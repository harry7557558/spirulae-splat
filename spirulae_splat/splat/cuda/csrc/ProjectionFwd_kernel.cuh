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



template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
__global__ void projection_fused_fwd_kernel(
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    float4 *__restrict__ aabbs,         // [B, C, N, 4]
    float *__restrict__ sorting_depths,  // [B, C, N, 1]
    float *__restrict__ radii,  // [N, 1]
    typename SplatPrimitive::ScreenBuffer splats_screen
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
    ProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(bid * C + cid);

    // Load splat
    typename SplatPrimitive::World splat_world;
    splat_world.load(splats_world, bid * N + gid);

    // Projection
    float sorting_depth;
    float4 aabb;
    float radius = 0.0f;
    typename SplatPrimitive::Screen splat_screen;
    splat_world.project<camera_model>(cam, splat_screen, aabb, sorting_depth, radius);

    // Save results
    aabb.x = fminf(fmaxf(aabb.x, 0.0f), image_width-1.0f);
    aabb.y = fminf(fmaxf(aabb.y, 0.0f), image_height-1.0f);
    aabb.z = fminf(fmaxf(aabb.z, 0.0f), image_width-1.0f);
    aabb.w = fminf(fmaxf(aabb.w, 0.0f), image_height-1.0f);
    if (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f) {
        splat_screen.store(splats_screen, idx);
        aabbs[idx] = aabb;
        sorting_depths[idx] = sorting_depth;
        atomicMax(&radii[idx%N], radius);
    } else {
        aabbs[idx] = {0.0f, 0.0f, 0.0f, 0.0f};
        sorting_depths[idx] = 0.0f;
    }
}


template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_fused_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    float4 *__restrict__ aabbs,         // [B, C, N, 4]
    float *__restrict__ sorting_depths,  // [B, C, N, 1]
    float *__restrict__ radii,  // [N, 1]
    typename SplatPrimitive::ScreenBuffer splats_screen
) {
    constexpr uint block = 128;
    projection_fused_fwd_kernel<SplatPrimitive, camera_model>
    <<<_CEIL_DIV(B*C*N, block), block, 0, stream>>>(
        B, C, N,
        splats_world, viewmats, intrins, dist_coeffs_buffer,
        image_width, image_height,
        aabbs, sorting_depths, radii, splats_screen
    );
}
