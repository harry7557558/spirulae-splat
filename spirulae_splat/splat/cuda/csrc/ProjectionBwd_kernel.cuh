#include <cuda_runtime.h>
#include <cstdint>

#include <Common.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif

#include "Common.cuh"

#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<
    typename SplatPrimitive,
    CameraModelType camera_model
>
__global__ void projection_fused_bwd_kernel(
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t *__restrict__ camera_ids,          // [nnz, 4]
    const int32_t *__restrict__ gaussian_ids,          // [nnz, 4]
    const float4 *__restrict__ aabb,          // [C, N, 4]
    // grad outputs
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splats_world,
    float *__restrict__ v_viewmats // [C, 4, 4] optional
) {
    uint32_t idx = cg::this_grid().thread_rank();
    uint32_t cid, gid;
    bool packed = (camera_ids != nullptr && gaussian_ids != nullptr);
    if (packed) {
        if (idx >= N) return;
        cid = camera_ids[idx];
        gid = gaussian_ids[idx];
    } else {
        // parallelize over C * N.
        if (idx >= C * N || (aabb[idx].z <= aabb[idx].x || aabb[idx].w <= aabb[idx].y)) {
            return;
        }
        cid = (idx / N) % C; // camera id
        gid = idx % N; // gaussian id
    }

    // Load camera
    viewmats += cid * 16;
    float4 intrin = intrins[cid];
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
    cam.dist_coeffs = dist_coeffs_buffer.load(cid);

    // Load splat
    typename SplatPrimitive::World splat_world;
    typename SplatPrimitive::Screen v_splat_screen;
    splat_world.load(splats_world, gid);
    v_splat_screen.load(v_splats_screen, idx);

    // Projection
    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero(v_splats_world, gid);
    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    splat_world.template project_vjp<camera_model>(cam, v_splat_screen, v_splat_world, v_R, v_t);

    // Save results
    v_splat_world.atomicStore(v_splats_world, gid);

    if (v_viewmats != nullptr) {
        auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
        auto warp_group_c = cg::labeled_partition(warp, cid);
        warpSum(v_R[0], warp_group_c);
        warpSum(v_R[1], warp_group_c);
        warpSum(v_R[2], warp_group_c);
        warpSum(v_t, warp_group_c);
        if (warp_group_c.thread_rank() == 0) {
            v_viewmats += cid * 16;
            #pragma unroll
            for (uint32_t i = 0; i < 3; i++) { // rows
                atomicAdd(v_viewmats + i * 4 + 0, v_R[i].x);
                atomicAdd(v_viewmats + i * 4 + 1, v_R[i].y);
                atomicAdd(v_viewmats + i * 4 + 2, v_R[i].z);
            }
            atomicAdd(v_viewmats + 0 * 4 + 3, v_t.x);
            atomicAdd(v_viewmats + 1 * 4 + 3, v_t.y);
            atomicAdd(v_viewmats + 2 * 4 + 3, v_t.z);
        }
    }

}

template<
    typename SplatPrimitive,
    CameraModelType camera_model
>
void projection_fused_bwd_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float * viewmats, // [C, 4, 4]
    const float4 * intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int32_t * camera_ids,          // [nnz, 4]
    const int32_t * gaussian_ids,          // [nnz, 4]
    const float4 * aabb,          // [C, N, 4]
    // grad outputs
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splats_world,
    float * v_viewmats // [C, 4, 4] optional
) {
    constexpr uint block = 128;
    projection_fused_bwd_kernel<SplatPrimitive, camera_model>
    <<<_CEIL_DIV(C*N, block), block, 0, stream>>>(
        C, N,
        splats_world, viewmats, intrins, dist_coeffs_buffer, image_width, image_height,
        camera_ids, gaussian_ids, aabb, v_splats_screen,
        v_splats_world, v_viewmats
    );
}
