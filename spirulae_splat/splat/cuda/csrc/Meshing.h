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
    int   bisection_iters = 20;  // bisection steps per cut edge
    int   max_cameras = 64;      // cap on cameras used (evenly subsampled)
    int   max_grid_res = 512;    // per-axis cap on the acceleration grid
    float grid_cell_factor = 2.0f; // cell size = factor * mean Gaussian radius
    int   num_threads = 0;       // 0 = all hardware threads (Delaunay + host)
    bool  verbose = true;
};

// Owns the device-side Gaussian buffers and a uniform spatial grid built over
// their support, plus the (optional) camera positions. Constructed once from
// host-side raw activated/un-activated splat parameters, then reused.
//
// All inputs are host pointers; the constructor uploads + activates internally.
//   means      : [N*3]  (x,y,z)
//   quats      : [N*4]  (w,x,y,z), un-normalized ok
//   log_scales : [N*3]  log std-devs (exp activation)
//   logit_opac : [N]    logit opacities (sigmoid activation)
//   cam_pos    : [C*3]  camera centers in the SAME frame as means (may be null)
class OccupancyEvaluator {
public:
    OccupancyEvaluator(
        const float* means, const float* quats,
        const float* log_scales, const float* logit_opac, int num_splats,
        const float* cam_pos, int num_cameras,
        const MeshingConfig& cfg);
    ~OccupancyEvaluator();

    OccupancyEvaluator(const OccupancyEvaluator&) = delete;
    OccupancyEvaluator& operator=(const OccupancyEvaluator&) = delete;

    // Generate the 7-points-per-kept-Gaussian point cloud. Returns interleaved
    // (x,y,z) doubles, ready for Delaunay3D. Size = 3 * num_points().
    void generate_point_cloud(std::vector<double>& xyz_out);

    // Evaluate the occupancy field at host points [n*3] doubles -> occ[n].
    void evaluate(const double* xyz, int n, float* occ_out);

    // For each cut edge (a,b) with endpoint occupancies (occ_a,occ_b) straddling
    // iso, bisection-search the crossing point. Endpoints are point-cloud
    // indices into the cloud used to build this evaluator's point set.
    // Writes crossing positions to xyz_out [n_edges*3] doubles.
    void bisect_edges(
        const double* cloud_xyz,         // [num_points*3], the MT point cloud
        const int32_t* edge_a, const int32_t* edge_b,
        const float* occ_a, const float* occ_b,
        int n_edges, double* xyz_out);

    int num_points() const { return num_points_; }
    int num_kept()   const { return num_kept_; }

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
    const float* log_scales, const float* logit_opac, int num_splats,
    const float* cam_pos, int num_cameras,
    const MeshingConfig& cfg,
    const std::string& output_path);

} // namespace meshing

#endif
