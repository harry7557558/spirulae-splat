#include "ProjectionEval3D.cuh"

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
__global__ void projection_eval3d_fwd_kernel(
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    // outputs
    int4 *__restrict__ aabbs,         // [B, C, N, 4]
    typename SplatPrimitive::WorldEval3D::Buffer splats_proj
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
        image_width, image_height, false,
        near_plane, far_plane
    };
    if (camera_model == gsplat::CameraModelType::FISHEYE) {
        if (dist_coeffs.radial_coeffs != nullptr)
            cam.radial_coeffs = dist_coeffs.radial_coeffs[bid * C + cid];
        if (dist_coeffs.tangential_coeffs != nullptr)
            cam.tangential_coeffs = dist_coeffs.tangential_coeffs[bid * C + cid];
        if (dist_coeffs.thin_prism_coeffs != nullptr)
            cam.thin_prism_coeffs = dist_coeffs.thin_prism_coeffs[bid * C + cid];
    }

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);

    // Projection
    int4 aabb;
    typename SplatPrimitive::WorldEval3D splat_proj;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp_eval3d(splat_world, cam, splat_proj, aabb);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye_eval3d(splat_world, cam, splat_proj, aabb);
        break;
    }

    // Save results
    if ((aabb.z-aabb.x)*(aabb.w-aabb.y) > 0) {
        splat_proj.saveBuffer(splats_proj, idx);
        aabbs[idx] = aabb;
    } else {
        aabbs[idx] = {0, 0, 0, 0};
    }
}

template<typename SplatPrimitive, gsplat::CameraModelType camera_model>
__global__ void projection_fused_eval3d_bwd_kernel(
    // fwd inputs
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int4 *__restrict__ aabb,          // [B, C, N, 4]
    // grad outputs
    typename SplatPrimitive::WorldEval3D::Buffer v_splats_proj,
    // grad inputs
    typename SplatPrimitive::World::Buffer v_splats_world,
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
        image_width, image_height, false,
    };
    if (camera_model == gsplat::CameraModelType::FISHEYE) {
        if (dist_coeffs.radial_coeffs != nullptr)
            cam.radial_coeffs = dist_coeffs.radial_coeffs[bid * C + cid];
        if (dist_coeffs.tangential_coeffs != nullptr)
            cam.tangential_coeffs = dist_coeffs.tangential_coeffs[bid * C + cid];
        if (dist_coeffs.thin_prism_coeffs != nullptr)
            cam.thin_prism_coeffs = dist_coeffs.thin_prism_coeffs[bid * C + cid];
    }

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);
    typename SplatPrimitive::WorldEval3D v_splat_proj =
        SplatPrimitive::WorldEval3D::load(v_splats_proj, idx);

    // Projection
    typename SplatPrimitive::World v_splat_world = SplatPrimitive::World::zero();
    float3x3 v_R = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    float3 v_t = {0.f, 0.f, 0.f};
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp_eval3d_vjp(splat_world, cam, v_splat_proj, v_splat_world, v_R, v_t);
        break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye_eval3d_vjp(splat_world, cam, v_splat_proj, v_splat_world, v_R, v_t);
        break;
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
    at::Tensor,  // aabb
    Vanilla3DGS::WorldEval3D::TensorTupleProj  // out splats
> projection_ewa_3dgs_eval3d_forward_tensor(
    const Vanilla3DGS::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    Vanilla3DGS::World::Tensor splats_world(in_splats);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    auto opt = splats_world.options();
    at::Tensor aabb = at::empty({C, N, 4}, opt.dtype(at::kInt));

    Vanilla3DGS::WorldEval3D::Tensor splats_proj =
        Vanilla3DGS::WorldEval3D::Tensor::empty(C, N, splats_world.options());

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_eval3d_fwd_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_eval3d_fwd_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(aabb, splats_proj.tupleProj());
}

std::tuple<
    at::Tensor,  // aabb
    OpaqueTriangle::WorldEval3D::TensorTupleProj  // out splats
> projection_opaque_triangle_eval3d_forward_tensor(
    const OpaqueTriangle::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    OpaqueTriangle::World::Tensor splats_world(in_splats);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    auto opt = splats_world.options();
    at::Tensor aabb = at::empty({C, N, 4}, opt.dtype(at::kInt));

    OpaqueTriangle::WorldEval3D::Tensor splats_proj =
        OpaqueTriangle::WorldEval3D::Tensor::empty(C, N, splats_world.options());

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_eval3d_fwd_kernel<OpaqueTriangle, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_eval3d_fwd_kernel<OpaqueTriangle, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(aabb, splats_proj.tupleProj());
}

std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_ewa_3dgs_eval3d_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::WorldEval3D::TensorTupleProj &v_splats_proj_tuple,
    const bool viewmats_requires_grad
) {
    Vanilla3DGS::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    Vanilla3DGS::WorldEval3D::Tensor v_splats_proj(v_splats_proj_tuple);

    Vanilla3DGS::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), v_splats_world.buffer(), \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_eval3d_bwd_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_fused_eval3d_bwd_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

std::tuple<
    OpaqueTriangle::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_opaque_triangle_eval3d_backward_tensor(
    // fwd inputs
    const OpaqueTriangle::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const OpaqueTriangle::WorldEval3D::TensorTupleProj &v_splats_proj_tuple,
    const bool viewmats_requires_grad
) {
    OpaqueTriangle::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    OpaqueTriangle::WorldEval3D::Tensor v_splats_proj(v_splats_proj_tuple);

    OpaqueTriangle::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    #define _LAUNCH_ARGS \
        <<<_CEIL_DIV(B*C*N, block), block>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), v_splats_world.buffer(), \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_eval3d_bwd_kernel<OpaqueTriangle, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_fused_eval3d_bwd_kernel<OpaqueTriangle, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}

