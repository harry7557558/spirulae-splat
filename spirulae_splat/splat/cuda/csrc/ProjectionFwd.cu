#include "ProjectionFwd.cuh"

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
    const typename SplatPrimitive::World::Buffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float *__restrict__ Ks,       // [B, C, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    // outputs
    int4 *__restrict__ aabbs,         // [B, C, N, 4]
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
        image_width, image_height,
        near_plane, far_plane,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(bid * C + cid);

    // Load splat
    typename SplatPrimitive::World splat_world =
        SplatPrimitive::World::load(splats_world, bid * N + gid);

    // Projection
    int4 aabb;
    typename SplatPrimitive::Screen splat_screen;
    switch (camera_model) {
    case gsplat::CameraModelType::PINHOLE: // perspective projection
        SplatPrimitive::project_persp(splat_world, cam, splat_screen, aabb);
        break;
    // case gsplat::CameraModelType::ORTHO: // orthographic projection
    //     SplatPrimitive::project_ortho(splat_world, cam, splat_screen, aabb);
    //     break;
    case gsplat::CameraModelType::FISHEYE: // fisheye projection
        SplatPrimitive::project_fisheye(splat_world, cam, splat_screen, aabb);
        break;
    }

    // Save results
    aabb.x = min(max(aabb.x, 0), image_width-1);
    aabb.y = min(max(aabb.y, 0), image_height-1);
    aabb.z = min(max(aabb.z, 0), image_width-1);
    aabb.w = min(max(aabb.w, 0), image_height-1);
    if ((aabb.z-aabb.x)*(aabb.w-aabb.y) > 0) {
        splat_screen.saveParamsToBuffer(splats_screen, idx);
        aabbs[idx] = aabb;
    } else {
        aabbs[idx] = {0, 0, 0, 0};
    }
}



template<typename SplatPrimitive>
std::tuple<
    at::Tensor,  // aabb
    typename SplatPrimitive::Screen::TensorTupleProj  // out splats
> launch_projection_fused_fwd_kernel(
    const typename SplatPrimitive::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    typename SplatPrimitive::World::Tensor splats_world(in_splats);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    uint32_t B = splats_world.batchSize();    // number of batches

    auto opt = splats_world.options();
    at::Tensor aabb = at::empty({C, N, 4}, opt.dtype(at::kInt));

    typename SplatPrimitive::Screen::Tensor splats_screen =
        SplatPrimitive::Screen::Tensor::allocProjFwd(C, N, splats_world.options());

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(B*C*N, block)>>>( \
            B, C, N, \
            splats_world.buffer(), viewmats.data_ptr<float>(), Ks.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, near_plane, far_plane, \
            (int4*)aabb.data_ptr<int32_t>(), splats_screen.buffer() \
        )

    constexpr uint block = 256;
    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_fwd_kernel<SplatPrimitive, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_fused_fwd_kernel<SplatPrimitive, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(aabb, splats_screen.tupleProjFwd());
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    Vanilla3DGS::Screen::TensorTupleProj  // out splats
> projection_3dgs_forward_tensor(
    // inputs
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
    return launch_projection_fused_fwd_kernel<Vanilla3DGS>(
        in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    MipSplatting::Screen::TensorTupleProj  // out splats
> projection_mip_forward_tensor(
    // inputs
    const MipSplatting::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<MipSplatting>(
        in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    Vanilla3DGUT::Screen::TensorTupleProj  // out splats
> projection_3dgut_forward_tensor(
    // inputs
    const Vanilla3DGUT::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<Vanilla3DGUT>(
        in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj  // out splats
> projection_3dgut_sv_forward_tensor(
    // inputs
    const SphericalVoronoi3DGUT_Default::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    int num_sv = std::get<5>(in_splats).size(-2);
    #define _CASE(n) \
        if (num_sv == n) return launch_projection_fused_fwd_kernel<SphericalVoronoi3DGUT<n>>( \
            in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs); \
    _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
    #undef _CASE
    throw std::invalid_argument("Unsupported num_sv");
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    OpaqueTriangle::Screen::TensorTupleProj  // out splats
> projection_opaque_triangle_forward_tensor(
    // inputs
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
    return launch_projection_fused_fwd_kernel<OpaqueTriangle>(
        in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    VoxelPrimitive::Screen::TensorTupleProj  // out splats
> projection_voxel_forward_tensor(
    // inputs
    const VoxelPrimitive::World::TensorTuple &in_splats,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<VoxelPrimitive>(
        in_splats, viewmats, Ks, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}
