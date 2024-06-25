#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"

// for f : R(n) -> R(m), J in R(m, n),
// v is cotangent in R(m), e.g. dL/df in R(m),
// compute vjp i.e. vT J -> R(n)
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

__global__ void rasterize_simple_backward_kernel(
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
    const float2* __restrict__ anisotropies,
    const float3& __restrict__ background,
    const int* __restrict__ final_index,
    const float* __restrict__ output_alpha,
    const float3* __restrict__ v_output,
    const float* __restrict__ v_output_alpha,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float* __restrict__ v_opacities,
    float2* __restrict__ v_anisotropies
);

__global__ void rasterize_backward_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_phi,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    const float3& __restrict__ background,
    const float2* __restrict__ depth_grads,
    const float3* __restrict__ depth_ref_im,
    const int* __restrict__ final_index,
    const float* __restrict__ output_alpha,
    const float4* __restrict__ output_depth_grad,
    const float* __restrict__ v_output_alpha,
    const float3* __restrict__ v_output,
    const float4* __restrict__ v_output_depth_grad,
    const float* __restrict__ v_output_reg_depth,
    const float* __restrict__ v_output_reg_normal,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float3* __restrict__ v_ch_coeffs,
    float* __restrict__ v_ch_coeffs_abs,
    float* __restrict__ v_opacities,
    float2* __restrict__ v_anisotropies,
    float2* __restrict__ v_depth_grad,
    float3* __restrict__ v_depth_ref_im
);

__global__ void rasterize_depth_backward_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    const int* __restrict__ final_index,
    const float* __restrict__ out_depth,
    const float2* __restrict__ out_visibility,
    const float* __restrict__ v_out_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float* __restrict__ v_opacities,
    float2* __restrict__ v_anisotropies
);

__device__ void project_cov3d_ewa_vjp(
    const float3 &mean3d,
    const float *cov3d,
    const float *viewmat,
    const float fx,
    const float fy,
    const float3 &v_cov2d,
    float3 &v_mean3d,
    float *v_cov3d
);

__device__ void scale_rot_to_cov3d_vjp(
    const float2 scale,
    const float4 quat,
    const float *v_cov3d,
    float2 &v_scale,
    float4 &v_quat
);

__device__ void projected_depth_grad_vjp(
    const glm::vec3 p, const glm::mat3 R,
    const float fx, const float fy,
    const float2 v_depth_grad,
    glm::vec3 &v_p_view, glm::mat3 &v_R
);
