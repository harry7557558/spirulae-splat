#include "ProjectionFwd.cuh"

#include <gsplat/Utils.cuh>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_fused_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    float4 *__restrict__ aabbs,         // [C, N, 4]
    float *__restrict__ sorting_depths,  // [C, N, 1]
    float *__restrict__ radii,  // [N, 1]
    typename SplatPrimitive::ScreenBuffer splats_screen
);



template<typename SplatPrimitive>
inline std::tuple<
    DeviceTensor2D<float4>,  // aabb [C, N]
    DeviceTensor2D<float>,   // sorting_depths [C, N]
    std::vector<DeviceTensorFloatND>  // out splats
> launch_projection_fused_fwd_kernel(
    const int64_t N,
    const std::vector<DeviceTensorFloatND> &in_splats,
    const float* viewmats_ptr,  // [C, 4, 4]
    const float4* intrins_ptr,  // [C, 4]
    const uint32_t C,
    const uint32_t image_width,
    const uint32_t image_height,
    const ssplat::CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
) {
    typename SplatPrimitive::WorldBuffer splats_world(in_splats);

    DeviceTensor2D<float4> aabb;
    aabb.resize("proj.aabb", C, N);
    DeviceTensor2D<float> sorting_depths;
    sorting_depths.resize("proj.depths", C, N);

    std::vector<DeviceTensorFloatND> splats_screen = SplatPrimitive::ScreenBuffer::empty_pool(C*splats_world.size(), "proj.screen");

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, N, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            aabb.data_ptr(), sorting_depths.data_ptr(), radii.data_ptr(), \
            splats_screen \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        projection_fused_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        projection_fused_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::EQUISOLID)
        projection_fused_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(aabb, sorting_depths, splats_screen);
}


// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgs_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
) {
    int sh_degree = Vanilla3DGS<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_fused_fwd_kernel<Vanilla3DGS<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs, radii);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}


// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_mip_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
) {
    int sh_degree = MipSplatting<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_fused_fwd_kernel<MipSplatting<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs, radii);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceTensor2D<float4>, DeviceTensor2D<float>, std::vector<DeviceTensorFloatND>
> projection_3dgut_forward(
    const int64_t num_splats, const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats, TorchTensorView intrins,
    const uint32_t image_width, const uint32_t image_height,
    const std::string camera_model, const TorchTensorView dist_coeffs,
    DeviceVector<float> radii
) {
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_fused_fwd_kernel<Vanilla3DGUT<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs, radii);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}
