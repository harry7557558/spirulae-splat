#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"

// compute the 2d gaussian parameters from 3d gaussian parameters
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

__global__ void rasterize_simple_forward(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3& __restrict__ background,
    int* __restrict__ final_index,
    float3* __restrict__ out_img,
    float* __restrict__ out_alpha
);

// compute output color image from binned and sorted gaussians
__global__ void rasterize_forward(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const unsigned ch_degree_r,
    const unsigned ch_degree_phi,
    const float* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    const float3& __restrict__ background,
    const float2* __restrict__ depth_grads,
    const float2* __restrict__ depth_normal_ref_im,
    int* __restrict__ final_index,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float3* __restrict__ out_depth_grad,
    float* __restrict__ out_reg_depth,
    float* __restrict__ out_reg_normal
);

// device helper to approximate projected 2d cov from 3d mean and cov
__device__ void project_cov3d_ewa(
    const float3 &mean3d,
    const float *cov3d,
    const float *viewmat,
    const float fx,
    const float fy,
    const float tan_fovx,
    const float tan_fovy,
    float3 &cov2d,
    float &comp
);

// device helper to get 3D covariance from scale and quat parameters
__device__ void scale_rot_to_cov3d(
    const float2 scale, const float4 quat, float *cov3d
);

__device__ float2 projected_depth_grad(
    const float3 p, const glm::mat3 R,
    const float fx, const float fy
);

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

__global__ void get_tile_bin_edges(
    const int num_intersects, const int64_t* __restrict__ isect_ids_sorted, int2* __restrict__ tile_bins
);
