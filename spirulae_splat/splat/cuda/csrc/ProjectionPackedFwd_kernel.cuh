#include <cuda_runtime.h>
#include <cstdint>

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "types.cuh"

#include <cooperative_groups.h>
namespace cg = cooperative_groups;



template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_packed_mask_kernel(
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    // outputs
    bool *__restrict__ intersection_mask  // [B, C, N]
) {
    // parallelize over B * C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= B * C * N) {
        return;
    }
    const uint32_t bid = idx / (C * N); // batch id
    const uint32_t cid = (idx / N) % C; // camera id
    const uint32_t gid = idx % N; // gaussian id

    // Load camera
    viewmats += bid * C * 16 + cid * 16;
    float4 intrin = intrins[bid * C + cid];
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
    cam.dist_coeffs = dist_coeffs_buffer.load(bid * C + cid);

    // Load splat
    // TODO: verify that SH is not loaded after compiler optimization
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);

    // Projection
    float4 aabb;
    typename SplatPrimitive::Screen splat_screen;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_screen, aabb);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye(splat_world, cam, splat_screen, aabb);
        break;
    }

    // Save results
    aabb.x = fminf(fmaxf(aabb.x, 0.0f), image_width-1.0f);
    aabb.y = fminf(fmaxf(aabb.y, 0.0f), image_height-1.0f);
    aabb.z = fminf(fmaxf(aabb.z, 0.0f), image_width-1.0f);
    aabb.w = fminf(fmaxf(aabb.w, 0.0f), image_height-1.0f);
    intersection_mask[idx] = (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f);
}


template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_packed_fwd_kernel(
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const int64_t* __restrict__ intersection_mask_scan,  // [B, C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    typename SplatPrimitive::Screen::Buffer splats_screen  // [nnz, ...]
) {
    // parallelize over B * C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= B * C * N) {
        return;
    }

    int64_t out_idx = (idx == 0 ? 0 : intersection_mask_scan[idx-1]);
    int64_t out_idx_1 = intersection_mask_scan[idx];
    if (out_idx_1 == out_idx)
        return;

    const uint32_t bid = idx / (C * N); // batch id
    const uint32_t cid = (idx / N) % C; // camera id
    const uint32_t gid = idx % N; // gaussian id

    // Load camera
    viewmats += bid * C * 16 + cid * 16;
    float4 intrin = intrins[bid * C + cid];
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
    cam.dist_coeffs = dist_coeffs_buffer.load(bid * C + cid);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);

    // Projection
    float4 aabb;
    typename SplatPrimitive::Screen splat_screen;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_screen, aabb);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye(splat_world, cam, splat_screen, aabb);
        break;
    }

    // Save results
    camera_ids[out_idx] = (int32_t)cid;
    gaussian_ids[out_idx] = (int32_t)gid;
    aabb.x = fminf(fmaxf(aabb.x, 0.0f), image_width-1.0f);
    aabb.y = fminf(fmaxf(aabb.y, 0.0f), image_height-1.0f);
    aabb.z = fminf(fmaxf(aabb.z, 0.0f), image_width-1.0f);
    aabb.w = fminf(fmaxf(aabb.w, 0.0f), image_height-1.0f);
    aabbs[out_idx] = aabb;
    splat_screen.saveParamsToBuffer(splats_screen, out_idx);
}


template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
void projection_packed_mask_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    // outputs
    bool *__restrict__ intersection_mask  // [B, C, N]
) {
    constexpr uint block = 128;
    projection_packed_mask_kernel<SplatPrimitive, camera_model>
    <<<_CEIL_DIV(B*C*N, block), block, 0, stream>>>(
        B, C, N,
        splats_world, viewmats, intrins, dist_coeffs_buffer,
        image_width, image_height, near_plane, far_plane,
        intersection_mask
    );
}

template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
void projection_packed_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const int64_t* __restrict__ intersection_mask_scan,  // [B, C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    typename SplatPrimitive::Screen::Buffer splats_screen  // [nnz, ...]
) {
    constexpr uint block = 128;
    projection_packed_fwd_kernel<SplatPrimitive, camera_model>
    <<<_CEIL_DIV(B*C*N, block), block, 0, stream>>>(
        B, C, N,
        splats_world, viewmats, intrins, dist_coeffs_buffer,
        image_width, image_height, near_plane, far_plane,
        intersection_mask_scan,
        camera_ids, gaussian_ids, aabbs, splats_screen
    );
}
