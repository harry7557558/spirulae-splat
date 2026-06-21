#pragma once

/*
 * MeshingRaster.cuh
 *
 * Rasterize-and-sample render driver for the meshing pipeline. Holds device
 * copies of the (raw) splat parameters and the camera intrinsics, and renders
 * one camera at a time via the existing 3DGUT projection + tile intersection +
 * (moment / color) rasterization. The per-camera render replaces the
 * per-point/per-camera LBVH traversal on the dataset path.
 *
 * Opaque RenderContext so Meshing.cu (which owns the OccupancyEvaluator) can
 * drive rendering without pulling in the heavy Tensor / projection headers.
 */

#include <cstdint>
#include <string>

namespace meshing {

struct RenderContext;

// Build a render context. All pointers are HOST arrays (uploaded internally).
//   means/quats/log_scales/logit_opac/features_dc : raw (un-activated) splat
//       params, [num_splats * {3,4,3,1,3}]; the projection activates them.
//   viewmats : [num_cameras*16] row-major world->cam 4x4
//   intrins  : [num_cameras*4]  fx, fy, cx, cy
//   dist     : [num_cameras*10] engine distortion layout (may be null => zeros)
//   widths   : [num_cameras] per-camera image width   (HOST array)
//   heights  : [num_cameras] per-camera image height  (HOST array)
RenderContext* render_context_create(
    const float* means, const float* quats, const float* log_scales,
    const float* logit_opac, const float* features_dc, int num_splats,
    const float* viewmats, const float* intrins, const float* dist,
    int num_cameras, const int* widths, const int* heights,
    const std::string& camera_model, int carve_k);

void render_context_destroy(RenderContext*);

int render_context_width(const RenderContext*, int cam_idx);
int render_context_height(const RenderContext*, int cam_idx);
int render_context_num_cameras(const RenderContext*);

// Render camera `cam_idx`'s occupancy moments into `d_moments` (device,
// [width*height] float3 = (m0, mean_depth, std_depth)). Caller owns d_moments.
void render_camera_moments(RenderContext*, int cam_idx, void* d_moments);

// Evaluate occupancy at device query points d_xyz [n*3], writing d_occ [n]
// (device). Renders each camera in cam_indices once, samples (project + bilinear
// + Gaussian-CDF), aggregates by min over the cameras that SEE each point;
// points seen by no camera get occupancy 0 (free / carved). Manages its own
// scratch internally.
void render_evaluate_occupancy(
    RenderContext*, const int* cam_indices, int num_cams,
    const float* d_xyz, int n, float* d_occ);

// Evaluate vertex color at device query points d_xyz [n*3], writing d_rgb [n*3]
// (device, RGB in [0,1]). Per camera, renders the full-ray DC color + moments,
// samples color (bilinear) weighted by transmittance T(z) until the point, and
// averages over cameras. Points with no usable view get rgb = (-1,-1,-1) so the
// caller can fall back to the static density-weighted color.
void render_evaluate_color(
    RenderContext*, const int* cam_indices, int num_cams,
    const float* d_xyz, int n, float* d_rgb);

} // namespace meshing
