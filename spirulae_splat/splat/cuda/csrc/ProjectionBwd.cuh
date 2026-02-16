#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include "types.cuh"

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include "common.cuh"


enum class HessianDiagonalOutputMode {
    None,
    Position,
    AllReasonable
};


template<
    typename SplatPrimitive,
    gsplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode
>
__global__ void projection_fused_bwd_kernel(
    // fwd inputs
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int4 *__restrict__ aabb,          // [B, C, N, 4]
    // grad outputs
    typename SplatPrimitive::Screen::Buffer v_splats_screen,
    typename SplatPrimitive::Screen::Buffer vr_splats_screen,
    typename SplatPrimitive::Screen::Buffer h_splats_screen,
    // grad inputs
    typename SplatPrimitive::World::Buffer v_splats_world,
    float3* vr_world_pos_buffer,
    float3* h_world_pos_buffer,
    typename SplatPrimitive::World::Buffer vr_splats_world,
    typename SplatPrimitive::World::Buffer h_splats_world,
    float *__restrict__ v_viewmats // [B, C, 4, 4] optional
) {
    // parallelize over B * C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= B * C * N || (aabb[idx].z-aabb[idx].x)*(aabb[idx].w-aabb[idx].y) <= 0) {
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
    typename SplatPrimitive::BwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(bid * C + cid);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);
    typename SplatPrimitive::Screen v_splat_screen =
        SplatPrimitive::Screen::load(v_splats_screen, idx);
    typename SplatPrimitive::Screen vr_splat_screen;
    typename SplatPrimitive::Screen h_splat_screen;
    if (hessian_diagonal_output_mode != HessianDiagonalOutputMode::None) {
        vr_splat_screen = SplatPrimitive::Screen::load(vr_splats_screen, idx);
        h_splat_screen = SplatPrimitive::Screen::load(h_splats_screen, idx);
    }

    // Projection
    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    float3 vr_world_pos = {0.f, 0.f, 0.f};
    float3 h_world_pos = {0.f, 0.f, 0.f};
    typename SplatPrimitive::World vr_splat_world = SplatPrimitive::World::zero();
    typename SplatPrimitive::World h_splat_world = SplatPrimitive::World::zero();
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::None) {
        switch (camera_model) {
        case gsplat::CameraModelType::PINHOLE: // perspective projection
            SplatPrimitive::project_persp_vjp(splat_world, cam, v_splat_screen, v_splat_world, v_R, v_t);
            break;
        // case gsplat::CameraModelType::ORTHO: // orthographic projection
        //     SplatPrimitive::project_ortho_vjp(splat_world, cam, v_splat_screen, v_splat_world, v_R, v_t);
        //     break;
        case gsplat::CameraModelType::FISHEYE: // fisheye projection
            SplatPrimitive::project_fisheye_vjp(splat_world, cam, v_splat_screen, v_splat_world, v_R, v_t);
            break;
        }
    } else if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        switch (camera_model) {
        case gsplat::CameraModelType::PINHOLE: // perspective projection
            SplatPrimitive::project_persp_vjp(splat_world, cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_world_pos, h_world_pos);
            break;
        case gsplat::CameraModelType::FISHEYE: // fisheye projection
            SplatPrimitive::project_fisheye_vjp(splat_world, cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_world_pos, h_world_pos);
            break;
        }
    } else if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        switch (camera_model) {
        case gsplat::CameraModelType::PINHOLE: // perspective projection
            SplatPrimitive::project_persp_vjp(splat_world, cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_splat_world, h_splat_world);
            break;
        case gsplat::CameraModelType::FISHEYE: // fisheye projection
            SplatPrimitive::project_fisheye_vjp(splat_world, cam, v_splat_screen,
                vr_splat_screen, h_splat_screen, v_splat_world, v_R, v_t, vr_splat_world, h_splat_world);
            break;
        }
    }

    // Save results
    v_splat_world.atomicAddGradientToBuffer(v_splats_world, bid * N + gid);
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        atomicAddFVec(&vr_world_pos_buffer[bid * N + gid], vr_world_pos);
        atomicAddFVec(&h_world_pos_buffer[bid * N + gid], h_world_pos);
    } else if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splat_world.atomicAddGradientToBuffer(vr_splats_world, bid * N + gid);
        h_splat_world.atomicAddGradientToBuffer(h_splats_world, bid * N + gid);
    }

    if (v_viewmats != nullptr) {
        auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
        auto warp_group_c = cg::labeled_partition(warp, cid);
        warpSum(v_R[0], warp_group_c);
        warpSum(v_R[1], warp_group_c);
        warpSum(v_R[2], warp_group_c);
        warpSum(v_t, warp_group_c);
        if (warp_group_c.thread_rank() == 0) {
            v_viewmats += bid * C * 16 + cid * 16;
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



template<typename SplatPrimitive, HessianDiagonalOutputMode hessian_diagonal_output_mode>
inline std::tuple<
    typename SplatPrimitive::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, typename SplatPrimitive::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, typename SplatPrimitive::World::TensorTuple>  // h_world_pos or h_splats
> _launch_projection_projection_fused_bwd_kernel(
    // fwd inputs
    const typename SplatPrimitive::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const typename SplatPrimitive::Screen::TensorTupleProj &v_splats_screen_tuple,
    const typename SplatPrimitive::Screen::TensorTupleProj *vr_splats_screen_tuple,
    const typename SplatPrimitive::Screen::TensorTupleProj *h_splats_screen_tuple,
    const bool viewmats_requires_grad
) {
    typename SplatPrimitive::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    typename SplatPrimitive::Screen::Tensor v_splats_screen(v_splats_screen_tuple);

    typename SplatPrimitive::World::Tensor v_splats_world = splats_world.allocProjBwd(false);

    typename SplatPrimitive::Screen::Tensor vr_splats_screen;
    typename SplatPrimitive::Screen::Tensor h_splats_screen;
    if (hessian_diagonal_output_mode != HessianDiagonalOutputMode::None) {
        vr_splats_screen = *vr_splats_screen_tuple;
        h_splats_screen = *h_splats_screen_tuple;
    }

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    at::Tensor vr_world_pos;
    at::Tensor h_world_pos;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        vr_world_pos = at::zeros({B, N, 3}, opt);
        h_world_pos = at::zeros({B, N, 3}, opt);
    }
    typename SplatPrimitive::World::Tensor vr_splats_world;
    typename SplatPrimitive::World::Tensor h_splats_world;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splats_world = vr_splats_world.allocProjBwd(true);
        h_splats_world = h_splats_world.allocProjBwd(true);
    }

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(B*C*N, block)>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_screen.buffer(), \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_screen.buffer() : typename SplatPrimitive::Screen::Buffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_screen.buffer() : typename SplatPrimitive::Screen::Buffer{}, \
            v_splats_world.buffer(), \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? (float3*)vr_world_pos.data_ptr<float>() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? (float3*)h_world_pos.data_ptr<float>() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? vr_splats_world.buffer() : typename SplatPrimitive::World::Buffer{}, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? h_splats_world.buffer() : typename SplatPrimitive::World::Buffer{}, \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_bwd_kernel<SplatPrimitive, gsplat::CameraModelType::PINHOLE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    // else if (camera_model == gsplat::CameraModelType::ORTHO)
    //     projection_fused_bwd_kernel<SplatPrimitive, gsplat::CameraModelType::ORTHO> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_fused_bwd_kernel<SplatPrimitive, gsplat::CameraModelType::FISHEYE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats, vr_world_pos, h_world_pos);
}

template<typename SplatPrimitive>
inline std::tuple<
    typename SplatPrimitive::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> launch_projection_projection_fused_bwd_kernel(
    // fwd inputs
    const typename SplatPrimitive::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const typename SplatPrimitive::Screen::TensorTupleProj &v_splats_screen_tuple,
    const bool viewmats_requires_grad
) {
    auto [v_splats, v_viewmats, vr_splats, h_splats] =
        _launch_projection_projection_fused_bwd_kernel
        <SplatPrimitive, HessianDiagonalOutputMode::None>
    (
        splats_world_tuple,
        viewmats,
        intrins,
        image_width,
        image_height,
        camera_model,
        dist_coeffs,
        aabb,
        v_splats_screen_tuple,
        nullptr,
        nullptr,
        viewmats_requires_grad
    );
    return std::make_tuple(v_splats, v_viewmats);
}
