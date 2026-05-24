#include "ProjectionPackedFwd.cuh"

#include <gsplat/Utils.cuh>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <cub/cub.cuh>

#include <Tensor.h>


template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_packed_mask_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    bool *__restrict__ intersection_mask  // [C, N]
);

template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_packed_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const int64_t* __restrict__ intersection_mask_scan,  // [C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    float *__restrict__ sorting_depths,         // [nnz]
    float *__restrict__ radii,  // [N]
    typename SplatPrimitive::ScreenBuffer splats_screen  // [nnz, ...]
);


template<typename SplatPrimitive>
inline std::tuple<
    DeviceVector<int32_t>,    // camera_ids [nnz]
    DeviceVector<int32_t>,    // gaussian_ids [nnz]
    DeviceVector<float4>,     // aabb [nnz]
    DeviceVector<float>,      // sorting_depths [nnz]
    DeviceVector<float>,      // radii [N]
    std::vector<DeviceTensorFloatND>  // out splats
> launch_projection_packed_fwd_kernel(
    const int64_t N,
    const std::vector<DeviceTensorFloatND> &in_splats,
    const float* viewmats_ptr,
    const float4* intrins_ptr,
    const uint32_t C,
    const uint32_t image_width,
    const uint32_t image_height,
    const ssplat::CameraModelType camera_model,
    const TorchTensorView dist_coeffs
) {
    typename SplatPrimitive::WorldBuffer splats_world(in_splats);

    // mask

    DeviceVector<bool> intersection_mask;
    intersection_mask.resize("proj_packed.mask", (int64_t)(C*N));

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, N, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            intersection_mask.data_ptr() \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        projection_packed_mask_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        projection_packed_mask_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::EQUISOLID)
        projection_packed_mask_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    // prefix sum
    DeviceVector<int64_t> intersection_mask_scan;
    intersection_mask_scan.resize("proj_packed.scan", (int64_t)(C*N));
    CUB_WRAPPER(cub::DeviceScan::InclusiveSum,
        intersection_mask.data_ptr(), intersection_mask_scan.data_ptr(), (int)(C*N));
    int64_t nnz = 0;
    if (C*N > 0)
        cudaMemcpy(&nnz, intersection_mask_scan.data_ptr() + (C*N - 1), sizeof(int64_t), cudaMemcpyDeviceToHost);

    // projection

    DeviceVector<int32_t> camera_ids; camera_ids.resize("proj_packed.camera_ids", nnz);
    DeviceVector<int32_t> gaussian_ids; gaussian_ids.resize("proj_packed.gaussian_ids", nnz);
    DeviceVector<float4> aabb; aabb.resize("proj_packed.aabb", nnz);
    DeviceVector<float> sorting_depths; sorting_depths.resize("proj_packed.depths", nnz);
    DeviceVector<float> radii; radii.resize("proj_packed.radii", splats_world.size()); radii.zero();

    std::vector<DeviceTensorFloatND> splats_screen = SplatPrimitive::ScreenBuffer::empty_pool(nnz, "proj_packed.screen");

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, N, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            intersection_mask_scan.data_ptr(), \
            camera_ids.data_ptr(), gaussian_ids.data_ptr(), \
            aabb.data_ptr(), sorting_depths.data_ptr(), radii.data_ptr(), \
            splats_screen \
        )

    if (camera_model == ssplat::CameraModelType::PINHOLE)
        projection_packed_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::FISHEYE)
        projection_packed_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else if (camera_model == ssplat::CameraModelType::EQUISOLID)
        projection_packed_fwd_kernel_wrapper<SplatPrimitive, ssplat::CameraModelType::EQUISOLID> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _LAUNCH_ARGS

    return std::make_tuple(camera_ids, gaussian_ids, aabb, sorting_depths, radii, splats_screen);
}


// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_3dgs_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs
) {
    int sh_degree = Vanilla3DGS<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_packed_fwd_kernel<Vanilla3DGS<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}


// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_mip_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs
) {
    int sh_degree = MipSplatting<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_packed_fwd_kernel<MipSplatting<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    DeviceVector<int32_t>, DeviceVector<int32_t>, DeviceVector<float4>,
    DeviceVector<float>, DeviceVector<float>, std::vector<DeviceTensorFloatND>
> projection_3dgut_packed_forward(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND> &in_splats,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs
) {
    int sh_degree = Vanilla3DGUT<0>::WorldBuffer(in_splats).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float* vm = (const float*)std::get<0>(viewmats);
    const float4* intr = (const float4*)std::get<0>(intrins);
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_packed_fwd_kernel<Vanilla3DGUT<n>>( \
            num_splats, in_splats, vm, intr, C, image_width, image_height, cmt(camera_model), dist_coeffs);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
    return {};
}
