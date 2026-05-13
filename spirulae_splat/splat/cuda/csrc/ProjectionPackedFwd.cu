#include "ProjectionFwd.cuh"

#include <gsplat/Utils.cuh>

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <cub/cub.cuh>


template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_packed_mask_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    bool *__restrict__ intersection_mask  // [B, C, N]
);

template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_packed_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,  // [B, N, ...]
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const int64_t* __restrict__ intersection_mask_scan,  // [B, C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    float *__restrict__ sorting_depths,         // [nnz]
    typename SplatPrimitive::ScreenBuffer splats_screen  // [nnz, ...]
);


template<typename SplatPrimitive>
inline std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    TensorList  // out splats
> launch_projection_packed_fwd_kernel(
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const ssplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    typename SplatPrimitive::WorldBuffer splats_world(in_splats);
    uint32_t N = splats_world.size();    // number of gaussians
    uint32_t C = viewmats.size(-3); // number of cameras
    // uint32_t B = splats_world.batchSize();    // number of batches
    uint32_t B = 1;  // TODO

    // mask

    at::Tensor intersection_mask = at::empty({B*C*N}, kTensorOptionBool());

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), B, C, N, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, \
            intersection_mask.data_ptr<bool>() \
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
    #if 0
    at::Tensor intersection_mask_scan = at::cumsum(intersection_mask, -1);
    #else
    at::Tensor intersection_mask_scan = at::empty({B*C*N}, kTensorOptionI64());
    {
        size_t temp_storage_bytes = 0;
        cub::DeviceScan::InclusiveSum(
            nullptr, temp_storage_bytes,
            intersection_mask.data_ptr<bool>(),
            intersection_mask_scan.data_ptr<int64_t>(),
            B*C*N);

        at::Tensor temp_storage = at::empty({(int64_t)temp_storage_bytes}, kTensorOptionByte());

        cub::DeviceScan::InclusiveSum(
            temp_storage.data_ptr<uint8_t>(), temp_storage_bytes,
            intersection_mask.data_ptr<bool>(),
            intersection_mask_scan.data_ptr<int64_t>(),
            B*C*N);
    }
    #endif
    int64_t nnz = intersection_mask_scan[-1].item<int64_t>();

    // projection

    at::Tensor camera_ids = at::empty({nnz,}, kTensorOptionI32());
    at::Tensor gaussian_ids = at::empty({nnz,}, kTensorOptionI32());
    at::Tensor aabb = at::empty({nnz, 4}, kTensorOptionF32());
    at::Tensor sorting_depths = at::empty({nnz,}, kTensorOptionF32());

    TensorList splats_screen = SplatPrimitive::ScreenBuffer::empty(nnz);

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), B, C, N, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, \
            intersection_mask_scan.data_ptr<int64_t>(), \
            camera_ids.data_ptr<int32_t>(), gaussian_ids.data_ptr<int32_t>(), \
            (float4*)aabb.data_ptr<float>(), sorting_depths.data_ptr<float>(), splats_screen \
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

    return std::make_tuple(camera_ids, gaussian_ids, aabb, sorting_depths, splats_screen);
}


// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    TensorList  // out splats
> projection_3dgs_packed_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_packed_fwd_kernel<Vanilla3DGS>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}


// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    TensorList  // out splats
> projection_mip_packed_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_packed_fwd_kernel<MipSplatting>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    TensorList  // out splats
> projection_3dgut_packed_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_packed_fwd_kernel<Vanilla3DGUT>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}


// // ================
// // SphericalVoronoi3DGUT
// // ================

// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     at::Tensor,  // camera_ids
//     at::Tensor,  // gaussian_ids
//     at::Tensor,  // aabb
//     SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj  // out splats
// > projection_3dgut_sv_packed_forward_tensor(
//     // inputs
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     int num_sv = std::get<5>(in_splats).size(-2);
//     #define _CASE(n) \
//         if (num_sv == n) return launch_projection_packed_fwd_kernel<SphericalVoronoi3DGUT<n>>( \
//             in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs); \
//     _CASE(2) _CASE(3) _CASE(4) _CASE(5) _CASE(6) _CASE(7) _CASE(8)
//     #undef _CASE
//     throw std::invalid_argument("Unsupported num_sv");
// }



// // ================
// // OpaqueTriangle
// // ================


// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     at::Tensor,  // camera_ids
//     at::Tensor,  // gaussian_ids
//     at::Tensor,  // aabb
//     OpaqueTriangle::Screen::TensorTupleProj  // out splats
// > projection_opaque_triangle_packed_forward_tensor(
//     // inputs
//     const OpaqueTriangle::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     return launch_projection_packed_fwd_kernel<OpaqueTriangle>(
//         in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
// }



// // ================
// // VoxelPrimitive
// // ================


// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     at::Tensor,  // camera_ids
//     at::Tensor,  // gaussian_ids
//     at::Tensor,  // aabb
//     VoxelPrimitive::Screen::TensorTupleProj  // out splats
// > projection_voxel_packed_forward_tensor(
//     // inputs
//     const VoxelPrimitive::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     return launch_projection_packed_fwd_kernel<VoxelPrimitive>(
//         in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
// }

