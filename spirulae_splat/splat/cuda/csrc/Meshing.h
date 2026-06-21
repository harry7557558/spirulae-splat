/*
 * Meshing.h
 *
 * Surface mesh extraction from a trained 3D Gaussian Splatting model.
 *
 * Pipeline (see MeshingHost.cpp for the orchestration):
 *   1. Sample a point cloud: 7 points per (kept) Gaussian -- the center plus
 *      6 points at +/- k*sigma along each principal axis, where
 *      k = sqrt(2 ln(opac / ALPHA_THRESHOLD)).
 *   2. Delaunay tetrahedralization of the point cloud (Delaunay3D.h).
 *   3. Marching tetrahedra over an occupancy field (iso = 0.5), with the
 *      crossing point on each cut edge found by bisection (not lerp).
 *   4. Manifold-preserving merge of short edges (local adaptive threshold).
 *   5. Write the mesh to a binary PLY.
 *
 * The occupancy field and all per-point GPU work live in Meshing.cu; the
 * OccupancyEvaluator owns the device-side Gaussian buffers + spatial grid and
 * is reused for both vertex evaluation and edge bisection.
 */
#ifndef SPIRULAE_MESHING_H
#define SPIRULAE_MESHING_H

#include <cstdint>
#include <string>
#include <vector>

namespace meshing {

struct MeshingConfig {
    float iso = 0.5f;            // isosurface level on occupancy in [0,1]
    float merge_factor = 1.0f;   // local merge threshold multiplier; <=0 disables
    int   bisection_iters = 3;  // bisection steps per cut edge; a linear
                                // interpolation finishes the crossing estimate
    int   max_cameras = -1;      // cap on cameras used (k-means subset)
    int   max_grid_res = 512;    // per-axis cap on the acceleration grid
    float grid_cell_factor = 2.0f; // cell size = factor * mean Gaussian radius
    int   num_threads = 0;       // 0 = all hardware threads (Delaunay + host)
    bool  verbose = true;

    // --- rasterize-and-sample occupancy/color (dataset path) ---
    // When camera intrinsics are supplied, occupancy/color come from rendering
    // each camera once (3DGUT) and sampling, instead of LBVH ray traversal.
    int   render_width = 0;      // 0 = use dataset-native width
    int   render_height = 0;     // 0 = use dataset-native height
    int   fill_hole_max_edges = 2; // boundary loops with <= this many edges are
                                   // triangulated; larger holes (sky) stay open
    bool  cull_unseen = true;    // drop final mesh verts seen by no camera

    // --- occupancy estimate robustness (render path) ---                                                         
    // The surface is lost when min-over-cameras meets per-view underestimation:                                   
    // one view that underestimates occupancy at a true surface point carves it,                                   
    // and adding cameras adds such outliers. These knobs make it robust.                                          
    //   occ_normalize : occ = Phi((z-mean)/std) decoupled from the pixel's total                                  
    //       accumulated alpha m0 (so the iso=0.5 surface sits at the mean depth,                                  
    //       not behind it). When false, occ = m0 * Phi(...) (old behavior).                                       
    //   occ_m0_thresh : a view whose pixel accumulated less than this alpha is                                    
    //       uncertain (sky / grazing / thin) and ABSTAINS instead of carving.                                     
    //   carve_k : aggregate occupancy as the k-th SMALLEST over the cameras that                                  
    //       see the point (k=1 = strict min / space carving; k>1 ignores the                                      
    //       k-1 most-underestimating views, trading carving power for robustness).                                
    bool  occ_normalize = true;                                                                                    
    float occ_m0_thresh = 0.5f;
    int   carve_k = 1;
};

// Full camera intrinsics/extrinsics for the rasterize-and-sample (dataset)
// path. All host pointers; when !valid() the evaluator falls back to the static
// LBVH path. Distortion uses the engine [C,10] layout expected by the 3DGUT
// projection: k1,k2,p1,p2,k3,k4,k5,k6,sx1,sy1 (tail-padded with zeros).
struct CameraParams {
    const float* viewmats = nullptr;    // [C*16] row-major world->cam 4x4
    const float* intrins = nullptr;     // [C*4]  fx, fy, cx, cy
    const float* dist_coeffs = nullptr; // [C*10] (may be null => treated as 0)
    const int* widths  = nullptr;       // [C] per-camera image width
    const int* heights = nullptr;       // [C] per-camera image height
    std::string camera_model;           // engine name, e.g. "PINHOLE"/"FISHEYE"
    bool valid() const {
        return viewmats && intrins && widths && heights;
    }
};

// Owns the device-side Gaussian buffers and an LBVH built over their support,
// plus the (optional) camera positions. Constructed once from host-side raw
// activated/un-activated splat parameters, then reused.
//
// All inputs are host pointers; the constructor uploads + activates internally.
//   means      : [N*3]  (x,y,z)
//   quats      : [N*4]  (w,x,y,z), un-normalized ok
//   log_scales : [N*3]  log std-devs (exp activation)
//   logit_opac : [N]    logit opacities (sigmoid activation)
//   features_dc: [N*3]  SH band-0 (DC) color coefficients
//   cam_pos    : [C*3]  camera centers in the SAME frame as means (may be null)
class OccupancyEvaluator {
public:
    OccupancyEvaluator(
        const float* means, const float* quats,
        const float* log_scales, const float* logit_opac, const float* features_dc,
        int num_splats,
        const float* cam_pos, int num_cameras,
        const CameraParams& cams,
        const MeshingConfig& cfg);
    ~OccupancyEvaluator();

    OccupancyEvaluator(const OccupancyEvaluator&) = delete;
    OccupancyEvaluator& operator=(const OccupancyEvaluator&) = delete;

    // Generate the 7-points-per-kept-Gaussian point cloud, interleaved (x,y,z).
    // Size = 3 * num_points(). float here; the caller widens to double only at
    // the Delaunay boundary (the one place that needs the extra precision).
    void generate_point_cloud(std::vector<float>& xyz_out);

    // Evaluate the occupancy field at host points [n*3] floats -> occ[n].
    void evaluate(const float* xyz, int n, float* occ_out);

    // For each cut edge (a,b) with endpoint occupancies (occ_a,occ_b) straddling
    // iso, bisection-search the crossing point (finished with a linear
    // interpolation). Endpoints are point-cloud indices into the cloud used to
    // build this evaluator's point set. Writes crossing positions to xyz_out
    // [n_edges*3] floats.
    void bisect_edges(
        const float* cloud_xyz,          // [num_points*3], the MT point cloud
        const int32_t* edge_a, const int32_t* edge_b,
        const float* occ_a, const float* occ_b,
        int n_edges, float* xyz_out);

    // Per-vertex RGB (in [0,1], interleaved) for host points [n*3]. Uses the
    // camera-aware compositing when cameras are present, else the static
    // density-weighted average, with a nearest-splat fallback.
    void colorize(const float* verts, int n, float* rgb_out);

    int num_points() const { return num_points_; }
    int num_kept()   const { return num_kept_; }

    // Debug: render camera `cam_idx`'s occupancy moments (m0, mean_depth,
    // std_depth) into `out` [width*height*3], row-major. Returns false if the
    // rasterize-and-sample render path is unavailable (no camera intrinsics).
    bool debug_render_moments(int cam_idx, std::vector<float>& out,
                              int& width, int& height);

private:
    struct Impl;
    Impl* impl_;
    int num_points_ = 0;
    int num_kept_ = 0;
};

// Full pipeline entry point. Writes a binary PLY to output_path.
// Returns true on success.
bool generate_mesh(
    const float* means, const float* quats,
    const float* log_scales, const float* logit_opac, const float* features_dc,
    int num_splats,
    const float* cam_pos, int num_cameras,
    const CameraParams& cams,
    const MeshingConfig& cfg,
    const std::string& output_path);

} // namespace meshing

#endif
