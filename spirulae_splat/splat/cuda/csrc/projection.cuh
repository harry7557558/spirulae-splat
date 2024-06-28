#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"


/* == AUTO HEADER GENERATOR - DO NOT CHANGE THIS LINE == */




// kernel function for projecting each gaussian on device
// each thread processes one gaussian
__global__ void project_gaussians_forward_kernel(
    const int num_points,
    const float3* __restrict__ means3d,
    const float2* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ viewmat,
    const float4 intrins,
    const dim3 tile_bounds,
    const unsigned block_width,
    const float clip_thresh,
    int4* __restrict__ bounds,
    int32_t* __restrict__ num_tiles_hit,
    float3* __restrict__ positions,
    float3* __restrict__ axes_u,
    float3* __restrict__ axes_v,
    float2* __restrict__ depth_grads
);


__global__ void project_gaussians_backward_kernel(
    const int num_points,
    const float3* __restrict__ means3d,
    const float2* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ viewmat,  // 3x4 row major
    const float4 intrins,
    const int* __restrict__ num_tiles_hit,
    const float3* __restrict__ v_positions,
    const float3* __restrict__ v_axes_u,
    const float3* __restrict__ v_axes_v,
    const float2* __restrict__ v_depth_grads,
    float3* __restrict__ v_means3d,
    float2* __restrict__ v_scales,
    float4* __restrict__ v_quats,
    float* __restrict__ v_viewmat
);


// kernel to map each intersection from tile ID and depth to a gaussian
// writes output to isect_ids and gaussian_ids
__global__ void map_gaussian_to_intersects(
    const int num_points,
    const float3* __restrict__ positions,
    int4* __restrict__ bounds,
    const int32_t* __restrict__ cum_tiles_hit,
    const dim3 tile_bounds,
    const unsigned block_width,
    int64_t* __restrict__ isect_ids,
    int32_t* __restrict__ gaussian_ids
);


// kernel to map sorted intersection IDs to tile bins
// expect that intersection IDs are sorted by increasing tile ID
// i.e. intersections of a tile are in contiguous chunks
__global__ void get_tile_bin_edges(
    const int num_intersects, const int64_t* __restrict__ isect_ids_sorted, int2* __restrict__ tile_bins
);
