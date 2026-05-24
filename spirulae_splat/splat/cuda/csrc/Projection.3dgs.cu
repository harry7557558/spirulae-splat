#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionHeteroFwd.cuh"
#include "ProjectionHeteroBwd.cuh"

#include "Primitive3DGS.cuh"


#if 0


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting depths
    at::Tensor,  // radii
    std::vector<DeviceTensorFloatND>  // out splats
> projection_3dgs_hetero_forward(
    // inputs
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
) {
    uint32_t N = Vanilla3DGS<0>::WorldBuffer(in_splats).size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = intersection_splat_id.size(-1);  // number of intersections

    at::Tensor camera_ids = at::empty({nnz}, kTensorOptionI32());
    at::Tensor gaussian_ids = at::empty({nnz}, kTensorOptionI32());
    at::Tensor aabb = at::empty({nnz, 4}, kTensorOptionF32());
    at::Tensor sorting_depths = at::empty({nnz}, kTensorOptionF32());
    at::Tensor radii = at::empty({nnz}, kTensorOptionF32());
    std::vector<DeviceTensorFloatND> splats_proj = Vanilla3DGS<0>::ScreenBuffer::empty(nnz);

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
        int sh_degree = Vanilla3DGS<0>::WorldBuffer(in_splats).sh_degree();
        sh_degree = min(sh_degree, max_sh_degree);
        #define LAUNCH(n) if (sh_degree == n) { \
            if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE) \
                projection_hetero_forward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS; \
            else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE) \
                projection_hetero_forward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS; \
            else if (cmt(camera_model) == ssplat::CameraModelType::EQUISOLID) \
                projection_hetero_forward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS; \
            else \
                throw std::runtime_error("Unsupported camera model"); \
        }
        LAUNCH(3)
        else LAUNCH(2)
        else LAUNCH(1)
        else LAUNCH(0)
        else LAUNCH(4)
        else throw std::runtime_error("Unsupported SH degree");
        #undef LAUNCH
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
    std::vector<DeviceTensorFloatND>,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_hetero_backward(
    // fwd inputs
    const std::vector<DeviceTensorFloatND> splats_world,
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor aabb,  // [nnz, 4]
    // grad outputs
    const std::vector<DeviceTensorFloatND> v_splats_proj,
    const bool viewmats_requires_grad
) {
    uint32_t N = Vanilla3DGS<0>::WorldBuffer(splats_world).size();  // number of splats
    uint32_t C = viewmats.size(-3);  // number of cameras
    uint32_t nnz = camera_ids.size(0);  // number of intersections

    std::vector<DeviceTensorFloatND> v_splats_world = Vanilla3DGS<0>::WorldBuffer::zeros_like(splats_world);

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
        int sh_degree = Vanilla3DGS<0>::WorldBuffer(splats_world).sh_degree();
        #define LAUNCH(n) if (sh_degree == n) { \
            if (cmt(camera_model) == ssplat::CameraModelType::PINHOLE) \
                projection_3dgs_hetero_backward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS; \
            else if (cmt(camera_model) == ssplat::CameraModelType::FISHEYE) \
                projection_3dgs_hetero_backward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS; \
            else if (cmt(camera_model) == ssplat::CameraModelType::EQUISOLID) \
                projection_3dgs_hetero_backward_kernel<Vanilla3DGS<n>, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS; \
            else \
                throw std::runtime_error("Unsupported camera model"); \
        }
        LAUNCH(3)
        else LAUNCH(2)
        else LAUNCH(1)
        else LAUNCH(0)
        else LAUNCH(4)
        else throw std::runtime_error("Unsupported SH degree");
        #undef LAUNCH
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(v_splats_world, v_viewmats);
}

#endif
