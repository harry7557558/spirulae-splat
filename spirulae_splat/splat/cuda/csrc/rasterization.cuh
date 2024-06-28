#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"


/* == AUTO HEADER GENERATOR - DO NOT CHANGE THIS LINE == */



__global__ void rasterize_simple_forward_kernel(
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
    int* __restrict__ final_index,
    float3* __restrict__ out_img,
    float* __restrict__ out_alpha
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


__global__ void rasterize_depth_forward_kernel(
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
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
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


__global__ void rasterize_forward_kernel(
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
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    const float3& __restrict__ background,
    const float2* __restrict__ depth_grads,
    const float3* __restrict__ depth_ref_im,
    int* __restrict__ final_index,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float4* __restrict__ out_depth_grad,
    float* __restrict__ out_reg_depth,
    float* __restrict__ out_reg_normal
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
