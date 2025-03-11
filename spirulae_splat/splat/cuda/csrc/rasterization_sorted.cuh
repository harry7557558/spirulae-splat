#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include "config.h"
#include "glm/glm/glm.hpp"


enum class PerPixelSortType {
	InsertionSort,
	QuickSort,
	HeapSort,
	RandomizedQuickSort,
};


/* == AUTO HEADER GENERATOR - DO NOT CHANGE THIS LINE == */



template<CameraType CAMERA_TYPE>
__global__ void rasterize_indices_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const float2* __restrict__ undistortion_map,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ num_intersects,
    int32_t* __restrict__ sorted_indices_,
    float* __restrict__ sorted_depths_
);


template <PerPixelSortType SORT_TYPE>
__global__ void sort_per_pixel_kernel(
    const unsigned num_pixels,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);


template<CameraType CAMERA_TYPE>
__global__ void rasterize_simple_sorted_forward_kernel(
    const int img_width,
    const int img_height,
    const float4 intrins,
    const float2* __restrict__ undistortion_map,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3 __restrict__ background,
    float3* __restrict__ out_img,
    float* __restrict__ out_alpha
);


__global__ void rasterize_simple_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3 __restrict__ background,
    const float* __restrict__ output_alpha,
    const float3* __restrict__ v_output,
    const float* __restrict__ v_output_alpha,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float* __restrict__ v_opacities
);


template <DepthMode DEPTH_MODE, CameraType CAMERA_TYPE>
__global__ void rasterize_depth_sorted_forward_kernel(
    const int img_width,
    const int img_height,
    const float4 intrins,
    const float2* __restrict__ undistortion_map,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
);


template <DepthMode DEPTH_MODE>
__global__ void rasterize_depth_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ final_index,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float* __restrict__ out_depth,
    const float2* __restrict__ out_visibility,
    const float* __restrict__ v_out_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float* __restrict__ v_opacities
);


__global__ void rasterize_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const float depth_reg_pairwise_factor,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    // const float3 __restrict__ background,
    const float* __restrict__ depth_ref_im,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float2* __restrict__ out_depth,
    float3* __restrict__ out_normal,
    float* __restrict__ out_reg_depth
);


__global__ void rasterize_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float depth_reg_pairwise_factor,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    // const float3 __restrict__ background,
    const float* __restrict__ depth_ref_im,
    const float* __restrict__ output_alpha,
    const float2* __restrict__ output_depth,
    const float* __restrict__ v_output_alpha,
    const float3* __restrict__ v_output,
    const float2* __restrict__ v_output_depth,
    const float3* __restrict__ v_output_normal,
    const float* __restrict__ v_output_reg_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float3* __restrict__ v_ch_coeffs,
    // float* __restrict__ v_ch_coeffs_abs,
    float* __restrict__ v_opacities,
    // float3* __restrict__ v_background,
    float* __restrict__ v_depth_ref_im
);


__global__ void rasterize_simplified_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float2* __restrict__ out_depth,  // { depth, depth^2 }
    float3* __restrict__ out_normal,
    float* __restrict__ out_depth_reg
);


__global__ void rasterize_simplified_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float* __restrict__ output_alpha,
    const float2* __restrict__ output_depth,
    const float* __restrict__ v_output_alpha,
    const float3* __restrict__ v_output_img,
    const float2* __restrict__ v_output_depth,
    const float3* __restrict__ v_output_normal,
    const float* __restrict__ v_output_depth_reg,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float* __restrict__ v_opacities
);
