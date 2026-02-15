#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "PrimitiveOpaqueTriangle.cuh"

#include "helpers.cuh"
#include <algorithm>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>


template<typename T>
__device__ __forceinline__ unsigned upper_bound(
    const T *arr, unsigned n, T value
) {
    unsigned left = 0, right = n;
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] <= value)
            left = mid + 1;
        else
            right = mid;
    }
    return left;
}

template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_hetero_forward_kernel(
    const uint32_t C,
    const uint32_t nnz,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const int32_t* __restrict__ intersection_count_map,  // [C+1]
    const int32_t* __restrict__ intersection_splat_id,  // [nnz]
    // outputs
    int64_t *__restrict__ camera_ids,    // [nnz]
    int64_t *__restrict__ gaussian_ids,  // [nnz]
    int4 *__restrict__ aabbs,    // [nnz, 4]
    typename SplatPrimitive::Screen::Buffer splats_proj
) {
    int32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= nnz)
        return;
    int32_t camera_idx = upper_bound(intersection_count_map, C+1, thread_idx) - 1;
    int32_t splat_idx = intersection_splat_id[thread_idx];

    // Load camera
    viewmats += camera_idx * 16;
    float4 intrin = intrins[camera_idx];
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    typename SplatPrimitive::FwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
        near_plane, far_plane,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(camera_idx);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, splat_idx);

    // Projection
    int4 aabb;
    typename SplatPrimitive::Screen splat_proj;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_proj, aabb);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye(splat_world, cam, splat_proj, aabb);
        break;
    }

    // Save results
    // aabb.x = min(max(aabb.x, 0), tile_width-1);
    // aabb.y = min(max(aabb.y, 0), tile_height-1);
    // aabb.z = min(max(aabb.z, 0), tile_width-1);
    // aabb.w = min(max(aabb.w, 0), tile_height-1);
    // TODO: use actual original cx, cy
    int offset_x = (int)roundf(0.5f * image_width - cx);
    int offset_y = (int)roundf(0.5f * image_height - cy);
    aabb.x = min(max(aabb.x + offset_x, 0), image_width-1) - offset_x;
    aabb.y = min(max(aabb.y + offset_y, 0), image_height-1) - offset_y;
    aabb.z = min(max(aabb.z + offset_x, 0), image_width-1) - offset_x;
    aabb.w = min(max(aabb.w + offset_y, 0), image_height-1) - offset_y;
    camera_ids[thread_idx] = camera_idx;
    gaussian_ids[thread_idx] = splat_idx;
    aabbs[thread_idx] = aabb;
    splat_proj.saveParamsToBuffer(splats_proj, thread_idx);
}


// TODO: refactor stuff in Projection.*.cu


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */

