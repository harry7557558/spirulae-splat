#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionHeteroFwd.cuh"
#include "ProjectionHeteroBwd.cuh"

#include "Primitive3DGS.cuh"




/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    MipSplatting::Screen::TensorTuple  // out splats
> projection_mip_hetero_forward_tensor(
    // inputs
    const MipSplatting::World::TensorTuple &in_splats_tensor,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    MipSplatting::World::Tensor in_splats = in_splats_tensor;
    uint32_t N = in_splats.batchSize() * in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    auto opt = in_splats.options();
    at::Tensor camera_ids = at::empty({nnz}, opt.dtype(at::kInt));
    at::Tensor gaussian_ids = at::empty({nnz}, opt.dtype(at::kInt));
    at::Tensor aabb = at::empty({nnz, 4}, opt.dtype(at::kFloat));
    MipSplatting::Screen::Tensor splats_proj =
        MipSplatting::Screen::Tensor::allocProjFwdPacked(nnz, opt);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, block)>>>( \
            C, nnz, \
            in_splats.buffer(), viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int32_t>(), gaussian_ids.data_ptr<int32_t>(), \
            (float4*)aabb.data_ptr<float>(), splats_proj.buffer() \
        )

    if (nnz != 0) {
        constexpr uint block = 128;
        if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<MipSplatting, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<MipSplatting, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
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
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_hetero_backward_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world_tuple,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const MipSplatting::Screen::TensorTuple &v_splats_proj_tuple,
    const bool viewmats_requires_grad
) {
    MipSplatting::World::Tensor splats_world(splats_world_tuple);
    uint32_t N = splats_world.batchSize() * splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    MipSplatting::Screen::Tensor v_splats_proj(v_splats_proj_tuple);

    MipSplatting::World::Tensor v_splats_world = splats_world.allocProjBwd(false);

    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = zeros_like_tensor(viewmats);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, 128)>>>( \
            C, N, nnz, \
            splats_world.buffer(), viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int32_t>(), gaussian_ids.data_ptr<int32_t>(), (float4*)aabb.data_ptr<float>(), \
            v_splats_proj.buffer(), v_splats_world.buffer(),  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<MipSplatting, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<MipSplatting, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world.tuple(), v_viewmats);
}
