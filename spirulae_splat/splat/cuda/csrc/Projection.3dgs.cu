#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionHeteroFwd.cuh"
#include "ProjectionHeteroBwd.cuh"


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    Vanilla3DGS::Screen::TensorTupleProj  // out splats
> projection_3dgs_forward_tensor(
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<Vanilla3DGS>(
        in_splats, viewmats, intrins, image_width, image_height, near_plane, far_plane, camera_model, dist_coeffs);
}


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
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    Vanilla3DGS::Screen::TensorTuple  // out splats
> projection_3dgs_hetero_forward_tensor(
    // inputs
    const Vanilla3DGS::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const float near_plane,
    const float far_plane,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    Vanilla3DGS::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kLong));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kInt));
    Vanilla3DGS::Screen::Tensor splats_proj =
        Vanilla3DGS::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, near_plane, far_plane, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), \
            (int4*)aabb.data_ptr<int32_t>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb,
        splats_proj.tupleProjFwdPacked()
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad,
    const bool sparse_grad
) {
    Vanilla3DGS::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    Vanilla3DGS::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    Vanilla3DGS::World::Tensor v_splats_world = splats_world.zeros_like();

    auto opt = splats_world.options();
    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = at::zeros_like(viewmats, opt);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int64_t>(), gaussian_ids.data_ptr<int64_t>(), (int4*)aabb.data_ptr<int32_t>(), \
            v_splats_proj.buffer(), sparse_grad, v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        constexpr uint block = 256;
        if (camera_model == gsplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGS, gsplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (camera_model == gsplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGS, gsplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}
