#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionHeteroFwd.cuh"
#include "ProjectionHeteroBwd.cuh"

#include "Primitive3DGUT.cuh"



/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_3dgut_hetero_forward_tensor(
    // inputs
    const TensorList &in_splats_tensor,
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
    Vanilla3DGUT::WorldBuffer in_splats = in_splats_tensor;
    uint32_t N = in_splats.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    at::Tensor camera_ids = at::empty({nnz}, kTensorOptionI32());
    at::Tensor gaussian_ids = at::empty({nnz}, kTensorOptionI32());
    at::Tensor aabb = at::empty({nnz, 4}, kTensorOptionF32());
    at::Tensor sorting_depths = at::empty({nnz}, kTensorOptionF32());
    at::Tensor radii = at::empty({nnz}, kTensorOptionF32());
    TensorList splats_proj = Vanilla3DGUT::ScreenBuffer::empty(nnz);

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, 128)>>>( \
            C, nnz, \
            in_splats, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            intersection_count_map.data_ptr<int32_t>(), intersection_splat_id.data_ptr<int32_t>(), \
            camera_ids.data_ptr<int32_t>(), gaussian_ids.data_ptr<int32_t>(), \
            (float4*)aabb.data_ptr<float>(), sorting_depths.data_ptr<float>(), radii.data_ptr<float>(), splats_proj \
        )

    if (nnz != 0) {
        if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE)
            projection_hetero_forward_kernel<Vanilla3DGUT, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE)
            projection_hetero_forward_kernel<Vanilla3DGUT, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::EQUISOLID)
            projection_hetero_forward_kernel<Vanilla3DGUT, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(
        camera_ids, gaussian_ids, aabb, sorting_depths, radii,
        splats_proj
    );
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    TensorList,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_hetero_backward_tensor(
    // fwd inputs
    const TensorList &splats_world_tuple,
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
    const TensorList &v_splats_proj_tuple,
    const bool viewmats_requires_grad
) {
    Vanilla3DGUT::WorldBuffer splats_world(splats_world_tuple);
    uint32_t N = splats_world.size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    Vanilla3DGUT::ScreenBuffer v_splats_proj(v_splats_proj_tuple);

    // Vanilla3DGUT::WorldBuffer v_splats_world = splats_world.allocProjBwd(false);
    TensorList v_splats_world = Vanilla3DGUT::WorldBuffer::zeros_like(splats_world);

    at::Tensor v_viewmats;
    if (viewmats_requires_grad)
        v_viewmats = zeros_like_tensor(viewmats);

    auto stream = at::cuda::getCurrentCUDAStream();

    #define _LAUNCH_ARGS \
        <<<_LAUNCH_ARGS_1D(nnz, 128)>>>( \
            C, N, nnz, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, tile_width, tile_height, \
            camera_ids.data_ptr<int32_t>(), gaussian_ids.data_ptr<int32_t>(), (float4*)aabb.data_ptr<float>(), \
            v_splats_proj, v_splats_world,  \
            viewmats_requires_grad ? v_viewmats.data_ptr<float>() : nullptr \
        )

    if (nnz != 0) {
        if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGUT, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGUT, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
        else if (cmt(camera_model) == ssplat::CameraModelType::EQUISOLID)
            projection_3dgs_hetero_backward_kernel<Vanilla3DGUT, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS;
        else
            throw std::runtime_error("Unsupported camera model");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world, v_viewmats);
}
