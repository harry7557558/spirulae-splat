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
