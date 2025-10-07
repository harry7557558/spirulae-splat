#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"


// refactored function arguments

#define _ARGS_rasterize_simple_forward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const int32_t* __restrict__ gaussian_ids_sorted, \
    const int2* __restrict__ tile_bins, \
    const float3* __restrict__ positions, \
    const float3* __restrict__ axes_u, \
    const float3* __restrict__ axes_v, \
    const float3* __restrict__ colors, \
    const float* __restrict__ opacities, \
    const float3 __restrict__ background, \
    int* __restrict__ final_index, \
    float3* __restrict__ out_img, \
    float* __restrict__ out_alpha

#define _ARGS_rasterize_depth_forward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const int32_t* __restrict__ gaussian_ids_sorted, \
    const int2* __restrict__ tile_bins, \
    const float3* __restrict__ positions, \
    const float3* __restrict__ axes_u, \
    const float3* __restrict__ axes_v, \
    const float* __restrict__ opacities, \
    int* __restrict__ final_index, \
    float* __restrict__ out_depth, \
    float2* __restrict__ out_visibility

#define _ARGS_rasterize_depth_backward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const int32_t* __restrict__ gaussian_ids_sorted, \
    const int2* __restrict__ tile_bins, \
    const float3* __restrict__ positions, \
    const float3* __restrict__ axes_u, \
    const float3* __restrict__ axes_v, \
    const float* __restrict__ opacities, \
    const int* __restrict__ final_index, \
    const float* __restrict__ out_depth, \
    const float2* __restrict__ out_visibility, \
    const float* __restrict__ v_out_depth, \
    float3* __restrict__ v_positions, \
    float2* __restrict__ v_positions_xy_abs, \
    float3* __restrict__ v_axes_u, \
    float3* __restrict__ v_axes_v, \
    float* __restrict__ v_opacities

#define _ARGS_rasterize_simplified_forward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const int32_t* __restrict__ gaussian_ids_sorted, \
    const int2* __restrict__ tile_bins, \
    const float3* __restrict__ positions, \
    const float3* __restrict__ axes_u, \
    const float3* __restrict__ axes_v, \
    const float3* __restrict__ colors, \
    const float* __restrict__ opacities, \
    int* __restrict__ final_index, \
    float* __restrict__ out_alpha, \
    float3* __restrict__ out_img, \
    float2* __restrict__ out_depth,  /* { depth, depth^2 } */ \
    float3* __restrict__ out_normal, \
    float* __restrict__ out_depth_reg

#define _ARGS_rasterize_simplified_backward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const int32_t* __restrict__ gaussian_ids_sorted, \
    const int2* __restrict__ tile_bins, \
    const float3* __restrict__ positions, \
    const float3* __restrict__ axes_u, \
    const float3* __restrict__ axes_v, \
    const float3* __restrict__ colors, \
    const float* __restrict__ opacities, \
    const int* __restrict__ final_index, \
    const float* __restrict__ output_alpha, \
    const float2* __restrict__ output_depth, \
    const float* __restrict__ v_output_alpha, \
    const float3* __restrict__ v_output_img, \
    const float2* __restrict__ v_output_depth, \
    const float3* __restrict__ v_output_normal, \
    const float* __restrict__ v_output_depth_reg, \
    float3* __restrict__ v_positions, \
    float2* __restrict__ v_positions_xy_abs, \
    float3* __restrict__ v_axes_u, \
    float3* __restrict__ v_axes_v, \
    float3* __restrict__ v_colors, \
    float* __restrict__ v_opacities

#define _ARGS_render_background_sh_forward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const float* rotation,  /* row major 3x3 */ \
    const unsigned sh_degree, \
    const float3* __restrict__ sh_coeffs_float3, \
    float3* __restrict__ out_img

#define _ARGS_render_background_sh_backward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const float* rotation,  /* row major 3x3 */ \
    const unsigned sh_degree, \
    const float3* __restrict__ sh_coeffs_float3, \
    const float3* __restrict__ out_color, \
    const float3* __restrict__ v_out_color, \
    float3* __restrict__ v_rotation, \
    float3* __restrict__ v_sh_coeffs


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_forward_kernel(
    _ARGS_render_background_sh_forward_kernel
);


template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_backward_kernel(
    _ARGS_render_background_sh_backward_kernel
);
