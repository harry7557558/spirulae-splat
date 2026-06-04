#include "ProjectionBwd.cuh"

#include <Common.cuh>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<
    typename SplatPrimitive,
    CameraModelType camera_model,
    HessianDiagonalOutputMode hessian_diagonal_output_mode
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
    typename SplatPrimitive::ScreenBuffer vr_splats_screen,
    typename SplatPrimitive::ScreenBuffer h_splats_screen,
    // grad inputs
    typename SplatPrimitive::WorldBuffer v_splats_world,
    float3* vr_world_pos_buffer,
    float3* h_world_pos_buffer,
    typename SplatPrimitive::WorldBuffer vr_splats_world,
    typename SplatPrimitive::WorldBuffer h_splats_world,
    float * v_viewmats // [C, 4, 4] optional
);


template<typename SplatPrimitive, HessianDiagonalOutputMode hessian_diagonal_output_mode>
inline void launch_projection_projection_fused_bwd_kernel(
    // fwd inputs
    int64_t N,
    const std::vector<DeviceTensorFloatND> &splats_world_tuple,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> v_splats_screen,
    const std::vector<DeviceTensorFloatND> vr_splats_screen,
    const std::vector<DeviceTensorFloatND> h_splats_screen,
    // returns
    std::vector<DeviceTensorFloatND> v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    DeviceVector<float3>* vr_world_pos_arg,
    const std::vector<DeviceTensorFloatND>* vr_splats_world_arg,
    DeviceVector<float3>* h_world_pos_arg,
    const std::vector<DeviceTensorFloatND>* h_splats_world_arg
) {
    typename SplatPrimitive::WorldBuffer splats_world(splats_world_tuple);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0]; // number of cameras
    const float* viewmats_ptr = (const float*)std::get<0>(viewmats);
    const float4* intrins_ptr = (const float4*)std::get<0>(intrins);

    DeviceVector<float3> vr_world_pos;
    DeviceVector<float3> h_world_pos;
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position) {
        vr_world_pos = *vr_world_pos_arg;
        h_world_pos = *h_world_pos_arg;
    }
    typename SplatPrimitive::WorldBuffer vr_splats_world;
    typename SplatPrimitive::WorldBuffer h_splats_world;
    if constexpr (hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable) {
        vr_splats_world = typename SplatPrimitive::WorldBuffer(*vr_splats_world_arg);
        h_splats_world = typename SplatPrimitive::WorldBuffer(*h_splats_world_arg);
    }

    if (camera_ids.data_ptr() && gaussian_ids.data_ptr()) {  // packed
        N = camera_ids.size();
        C = 1;
    }

    if (C*N == 0)
        return;

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, N, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            camera_ids.data_ptr(), \
            gaussian_ids.data_ptr(), \
            aabb.data_ptr(), \
            v_splats_screen, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? vr_splats_screen : typename SplatPrimitive::ScreenBuffer{}, \
            hessian_diagonal_output_mode != HessianDiagonalOutputMode::None ? h_splats_screen : typename SplatPrimitive::ScreenBuffer{}, \
            v_splats_world, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? vr_world_pos.data_ptr() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::Position ? h_world_pos.data_ptr() : nullptr, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? vr_splats_world : typename SplatPrimitive::WorldBuffer{}, \
            hessian_diagonal_output_mode == HessianDiagonalOutputMode::AllReasonable ? h_splats_world : typename SplatPrimitive::WorldBuffer{}, \
            v_viewmats != nullptr ? v_viewmats->data_ptr() : nullptr \
        )

    if (camera_model == CameraModelType::PINHOLE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, CameraModelType::PINHOLE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::FISHEYE)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, CameraModelType::FISHEYE, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::EQUISOLID)
        projection_fused_bwd_kernel_wrapper<SplatPrimitive, CameraModelType::EQUISOLID, hessian_diagonal_output_mode> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

}



// ================
// Vanilla3DGS
// ================


/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats
) {
    int sh_degree = Vanilla3DGS<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGS<n>, HessianDiagonalOutputMode::None>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {}, \
            v_splats_world, v_viewmats, nullptr, nullptr, nullptr, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_with_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    const std::vector<DeviceTensorFloatND> &vr_splats_world,
    const std::vector<DeviceTensorFloatND> &h_splats_world
) {
    int sh_degree = Vanilla3DGS<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGS<n>, HessianDiagonalOutputMode::AllReasonable>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, nullptr, &vr_splats_world, nullptr, &h_splats_world);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_with_position_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    DeviceVector<float3>* vr_world_pos,
    DeviceVector<float3>* h_world_pos
) {
    int sh_degree = Vanilla3DGS<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGS<n>, HessianDiagonalOutputMode::Position>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, vr_world_pos, nullptr, h_world_pos, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}



// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats
) {
    int sh_degree = MipSplatting<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<MipSplatting<n>, HessianDiagonalOutputMode::None>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {}, \
            v_splats_world, v_viewmats, nullptr, nullptr, nullptr, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_with_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    const std::vector<DeviceTensorFloatND> &vr_splats_world,
    const std::vector<DeviceTensorFloatND> &h_splats_world
) {
    int sh_degree = MipSplatting<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<MipSplatting<n>, HessianDiagonalOutputMode::AllReasonable>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, nullptr, &vr_splats_world, nullptr, &h_splats_world);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_with_position_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    DeviceVector<float3>* vr_world_pos,
    DeviceVector<float3>* h_world_pos
) {
    int sh_degree = MipSplatting<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<MipSplatting<n>, HessianDiagonalOutputMode::Position>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, vr_world_pos, nullptr, h_world_pos, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats
) {
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT<n>, HessianDiagonalOutputMode::None>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, {}, {}, \
            v_splats_world, v_viewmats, nullptr, nullptr, nullptr, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_with_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    const std::vector<DeviceTensorFloatND> &vr_splats_world,
    const std::vector<DeviceTensorFloatND> &h_splats_world
) {
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT<n>, HessianDiagonalOutputMode::AllReasonable>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, nullptr, &vr_splats_world, nullptr, &h_splats_world);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}

/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_with_position_hessian_diagonal(
    // fwd inputs
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &splats_world,
    TorchTensorView viewmats,  // [C, 4, 4]
    TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    // fwd outputs
    const DeviceVector<int32_t> camera_ids,  // [nnz] or null
    const DeviceVector<int32_t> gaussian_ids,  // [nnz] or null
    const DeviceTensor2D<float4> aabb,  // [C, N]
    // grad outputs
    const std::vector<DeviceTensorFloatND> &v_splats_screen,
    const std::vector<DeviceTensorFloatND> &vr_splats_screen,
    const std::vector<DeviceTensorFloatND> &h_splats_screen,
    // returns
    const std::vector<DeviceTensorFloatND> &v_splats_world,
    DeviceTensor2D<float>* v_viewmats,
    DeviceVector<float3>* vr_world_pos,
    DeviceVector<float3>* h_world_pos
) {
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_projection_fused_bwd_kernel<Vanilla3DGUT<n>, HessianDiagonalOutputMode::Position>( \
            num_splats, splats_world, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, vr_splats_screen, h_splats_screen, \
            v_splats_world, v_viewmats, vr_world_pos, nullptr, h_world_pos, nullptr);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}
