#include "ProjectionEWA3DGS.cuh"

#include <ATen/Dispatch.h>
#include <ATen/core/Tensor.h>
#include <ATen/cuda/Atomic.cuh>
#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>


template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_fused_fwd_kernel(
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const bool antialiased,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    // outputs
    int2 *__restrict__ radii,         // [B, C, N, 2]
    float *__restrict__ depths,       // [B, C, N]
    typename SplatPrimitive::Screen::Buffer splats_screen
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
    Ks += bid * C * 9 + cid * 9;
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    typename SplatPrimitive::FwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height, antialiased,
        near_plane, far_plane, radius_clip
    };

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);

    // Projection
    int2 radius;
    float depth;
    typename SplatPrimitive::Screen splat_screen;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_screen, radius, depth);
        break;
    // case gsplat::CameraModelType::ORTHO: // orthographic projection
    //     SplatPrimitive::project_ortho(splat_world, cam, splat_screen, radius, depth);
    //     break;
    // case gsplat::CameraModelType::FISHEYE: // fisheye projection
    //     SplatPrimitive::project_fisheye(splat_world, cam, splat_screen, radius, depth);
    //     break;
    }

    // Save results
    radii[idx] = radius;
    if (radius.x * radius.y > 0) {
        depths[idx] = depth;
        splat_screen.saveBuffer(splats_screen, idx);
    }
}

std::tuple<
    at::Tensor,  // radii
    at::Tensor,  // depths
    Vanilla3DGS::Screen::TensorTuple  // out splats
> projection_ewa_3dgs_forward_tensor(
    const bool antialiased,
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    const gsplat::CameraModelType camera_model
) {
    Vanilla3DGS::World::Tensor splats_world(in_splats);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    auto opt = splats_world.options();
    at::Tensor radii = at::empty({C, N, 2}, opt.dtype(at::kInt));
    at::Tensor means2d = at::empty({C, N, 2}, opt);
    at::Tensor depths = at::empty({C, N}, opt);
    at::Tensor conics = at::empty({C, N, 3}, opt);
    at::Tensor out_opacities = at::empty({C, N}, opt);

    Vanilla3DGS::Screen::Tensor splats_screen = std::make_tuple(
        means2d, conics, out_opacities
    );

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            antialiased, \
            splats_world.buffer(), \
            viewmats.data_ptr<float>(), \
            Ks.data_ptr<float>(), \
            image_width, \
            image_height, \
            eps2d, \
            near_plane, \
            far_plane, \
            radius_clip, \
            (int2*)radii.data_ptr<int32_t>(), \
            depths.data_ptr<float>(), \
            splats_screen.buffer() \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_fwd_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    // else if (camera_model == gsplat::CameraModelType::ORTHO)
    //     projection_fused_fwd_kernel<Vanilla3DGS, gsplat::CameraModelType::ORTHO> _LAUNCH_ARGS;
    // else if (camera_model == gsplat::CameraModelType::FISHEYE)
    //     projection_fused_fwd_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(radii, depths, std::make_tuple(means2d, conics, out_opacities));
}

template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_fused_bwd_kernel(
    // fwd inputs
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const bool antialiased,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    // fwd outputs
    const int2 *__restrict__ radii,          // [B, C, N, 2]
    // grad outputs
    const float *__restrict__ v_depths,        // [B, C, N]
    typename SplatPrimitive::Screen::Buffer v_splats_screen,
    // grad inputs
    typename SplatPrimitive::World::Buffer v_splats_world,
    float *__restrict__ v_viewmats // [B, C, 4, 4] optional
) {
    // parallelize over B * C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= B * C * N || radii[idx].x * radii[idx].y <= 0) {
        return;
    }
    const uint32_t bid = idx / (C * N); // batch id
    const uint32_t cid = (idx / N) % C; // camera id
    const uint32_t gid = idx % N; // gaussian id

    // Load camera
    viewmats += bid * C * 16 + cid * 16;
    Ks += bid * C * 9 + cid * 9;
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    typename SplatPrimitive::BwdProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height, antialiased,
    };

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);
    typename SplatPrimitive::Screen v_splat_screen =
        SplatPrimitive::Screen::load(v_splats_screen, idx);

    // Projection
    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp_vjp(splat_world, cam, v_splat_screen, 0.f, v_splat_world, v_R, v_t);
        break;
    // case gsplat::CameraModelType::ORTHO: // orthographic projection
    //     SplatPrimitive::project_ortho_vjp(splat_world, cam, v_splat_screen, 0.f, v_splat_world, v_R, v_t);
    //     break;
    // case gsplat::CameraModelType::FISHEYE: // fisheye projection
    //     SplatPrimitive::project_fisheye_vjp(splat_world, cam, v_splat_screen, 0.f, v_splat_world, v_R, v_t);
    //     break;
    }

    // Save results
    auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
    auto warp_group_g = cg::labeled_partition(warp, gid);
    v_splat_world.atomicAddBuffer(v_splats_world, bid * N + gid);

    if (v_viewmats != nullptr) {
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

std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_ewa_3dgs_backward_tensor(
    const bool antialiased,
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const gsplat::CameraModelType camera_model,
    // fwd outputs
    const at::Tensor radii,                       // [..., C, N, 2]
    // grad outputs
    const at::Tensor v_depths,                      // [..., C, N]
    const Vanilla3DGS::Screen::TensorTuple &v_splats_screen_tuple,
    const bool viewmats_requires_grad
) {
    Vanilla3DGS::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    Vanilla3DGS::Screen::Tensor v_splats_screen(v_splats_screen_tuple);

    Vanilla3DGS::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            antialiased, \
            splats_world.buffer(), \
            viewmats.data_ptr<float>(), \
            Ks.data_ptr<float>(), \
            image_width, \
            image_height, \
            eps2d, \
            (int2*)radii.data_ptr<int32_t>(), \
            v_depths.data_ptr<float>(), \
            v_splats_screen.buffer(), \
            v_splats_world.buffer(), \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_bwd_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    // else if (camera_model == gsplat::CameraModelType::ORTHO)
    //     projection_fused_bwd_kernel<Vanilla3DGS, gsplat::CameraModelType::ORTHO> _LAUNCH_ARGS;
    // else if (camera_model == gsplat::CameraModelType::FISHEYE)
    //     projection_fused_bwd_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

