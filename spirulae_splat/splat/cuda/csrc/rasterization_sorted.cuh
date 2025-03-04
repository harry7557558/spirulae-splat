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



__global__ void rasterize_indices_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
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
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);


__global__ void rasterize_simple_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3& __restrict__ background,
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
    const float3& __restrict__ background,
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


template <DepthMode DEPTH_MODE>
__global__ void rasterize_depth_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
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
