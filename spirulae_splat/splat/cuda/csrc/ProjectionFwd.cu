#include "ProjectionFwd.cuh"

#include <gsplat/Utils.cuh>

#include <c10/cuda/CUDAStream.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;


template<typename SplatPrimitive, ssplat::CameraModelType camera_model>
void projection_fused_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t B,
    const uint32_t C,
    const uint32_t N,
    const typename SplatPrimitive::WorldBuffer splats_world,
    const float *__restrict__ viewmats, // [B, C, 4, 4]
    const float4 *__restrict__ intrins,  // [B, C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    float4 *__restrict__ aabbs,         // [B, C, N, 4]
    float *__restrict__ sorting_depths,  // [B, C, N, 1]
    float *__restrict__ radii,  // [N, 1]
    typename SplatPrimitive::ScreenBuffer splats_screen
);



template<typename SplatPrimitive>
inline std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> launch_projection_fused_fwd_kernel(
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

    at::Tensor aabb = at::empty({C, N, 4}, kTensorOptionF32());
    at::Tensor sorting_depths = at::empty({C, N}, kTensorOptionF32());
    at::Tensor radii = at::empty({N}, kTensorOptionF32());
    set_zero_tensor(radii);

    TensorList splats_screen = SplatPrimitive::ScreenBuffer::empty(C*N);

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)at::cuda::getCurrentCUDAStream(), B, C, N, \
            splats_world, viewmats.data_ptr<float>(), (float4*)intrins.data_ptr<float>(), dist_coeffs, \
            image_width, image_height, \
            (float4*)aabb.data_ptr<float>(), sorting_depths.data_ptr<float>(), radii.data_ptr<float>(), \
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

    return std::make_tuple(aabb, sorting_depths, radii, splats_screen);
}


// ================
// Vanilla3DGS
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_3dgs_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<Vanilla3DGS>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}


// ================
// MipSplatting
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_mip_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<MipSplatting>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}



// ================
// Vanilla3DGUT
// ================

/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // aabb
    at::Tensor,  // sorting_depths
    at::Tensor,  // radii
    TensorList  // out splats
> projection_3dgut_forward_tensor(
    // inputs
    const TensorList &in_splats,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs
) {
    return launch_projection_fused_fwd_kernel<Vanilla3DGUT>(
        in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
}


// // ================
// // SphericalVoronoi3DGUT
// // ================

// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     at::Tensor,  // aabb
//     TensorList  // out splats
// > projection_3dgut_sv_forward_tensor(
//     // inputs
//     const TensorList &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     int num_sv = std::get<5>(in_splats).size(-2);
//     #define _CASE(n) \
//         if (num_sv == n) return launch_projection_fused_fwd_kernel<SphericalVoronoi3DGUT<n>>( \
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
//     at::Tensor,  // aabb
//     OpaqueTriangle::Screen::TensorTupleProj  // out splats
// > projection_opaque_triangle_forward_tensor(
//     // inputs
//     const OpaqueTriangle::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     return launch_projection_fused_fwd_kernel<OpaqueTriangle>(
//         in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
// }



// // ================
// // VoxelPrimitive
// // ================


// /*[AutoHeaderGeneratorExport]*/
// std::tuple<
//     at::Tensor,  // aabb
//     VoxelPrimitive::Screen::TensorTupleProj  // out splats
// > projection_voxel_forward_tensor(
//     // inputs
//     const VoxelPrimitive::World::TensorTuple &in_splats,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,   // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs
// ) {
//     return launch_projection_fused_fwd_kernel<VoxelPrimitive>(
//         in_splats, viewmats, intrins, image_width, image_height, cmt(camera_model), dist_coeffs);
// }

