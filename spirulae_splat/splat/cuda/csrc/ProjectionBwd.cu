#include "ProjectionBwd.cuh"

#include <gsplat/Utils.cuh>

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<
    typename SplatPrimitive,
    gsplat::CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode
>
void projection_fused_bwd_kernel_wrapper(
    cudaStream_t stream,
    // fwd inputs
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::World::Buffer splats_world,
    const float * viewmats, // [B, C, 4, 4]
    const float4 * intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // fwd outputs
    const int4 * aabb,          // [B, C, N, 4]
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
    float * v_viewmats // [B, C, 4, 4] optional
);


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
        v_viewmats = zeros_like<float>(viewmats);

    at::Tensor vr_world_pos;
    at::Tensor h_world_pos;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        vr_world_pos = at::empty({B, N, 3}, opt);
        h_world_pos = at::empty({B, N, 3}, opt);
        set_zero<float>(vr_world_pos);
        set_zero<float>(h_world_pos);
    }
    typename SplatPrimitive::World::Tensor vr_splats_world;
    typename SplatPrimitive::World::Tensor h_splats_world;
    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splats_world = splats_world.allocProjBwd(true);
        h_splats_world = splats_world.allocProjBwd(true);
    }

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), B, C, N, \
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

    if (camera_model == gsplat::CameraModelType::PINHOLE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, gsplat::CameraModelType::PINHOLE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else if (camera_model == gsplat::CameraModelType::FISHEYE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, gsplat::CameraModelType::FISHEYE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    if (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable)
        return std::make_tuple(v_splats_world.tuple(), v_viewmats, vr_splats_world.tuple(), h_splats_world.tuple());
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


// ================
// Vanilla3DGS
// ================


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,  // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<Vanilla3DGS>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgs_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgs_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<Vanilla3DGS, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}



// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_backward_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<MipSplatting>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>  // h_world_pos or h_splats
> projection_mip_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &vr_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>  // h_world_pos or h_splats
> projection_mip_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &vr_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<MipSplatting, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_backward_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgut_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, HessianDiagonalOutputMode::AllReasonable>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgut_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
) {
    return _launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT, HessianDiagonalOutputMode::Position>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, &vr_splats_screen, &h_splats_screen, viewmats_requires_grad);
}



// ================
// SphericalVoronoi3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    SphericalVoronoi3DGUT_Default::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_sv_backward_tensor(
    // fwd inputs
    const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    int num_sv = std::get<5>(splats_world).size(-2);
    #define _CASE(n) \
        if (num_sv == n) return launch_projection_projection_fused_bwd_kernel<SphericalVoronoi3DGUT<n>>( \
            splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs, \
            aabb, v_splats_screen, viewmats_requires_grad);
    _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
    #undef _CASE
    throw std::invalid_argument("Unsupported num_sv");
}



// ================
// OpaqueTriangle
// ================


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    OpaqueTriangle::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_opaque_triangle_backward_tensor(
    // fwd inputs
    const OpaqueTriangle::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const OpaqueTriangle::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<OpaqueTriangle>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}


// ================
// VoxelPrimitive
// ================


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    VoxelPrimitive::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_voxel_backward_tensor(
    // fwd inputs
    const VoxelPrimitive::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
) {
    return launch_projection_projection_fused_bwd_kernel<VoxelPrimitive>(
        splats_world, viewmats, intrins, image_width, image_height, camera_model, dist_coeffs,
        aabb, v_splats_screen, viewmats_requires_grad);
}

