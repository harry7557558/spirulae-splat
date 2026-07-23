#include "kernels/projection/ProjectionBwdQuantGrad.cuh"

#include <core/Common.cuh>
#include <optional>

#include "primitives/Primitive3DGS.cuh"
#include "primitives/Primitive3DGUT.cuh"

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include <cub/cub.cuh>


// Thin per-(primitive, camera_model) launcher, defined by the codegen'd
// instantiations in ins/*.cu (extern here). Dispatches VALUE_BITS internally.
template<
    typename SplatPrimitive,
    CameraModelType camera_model
>
void projection_bwd_quantgrad_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    const uint32_t num_sh_buffer,
    typename SplatPrimitive::WorldBuffer splats_world,
    const float * viewmats,
    const float4 * intrins,
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const int32_t * camera_id_bounds,
    const int32_t * camera_ids,
    const int32_t * perm,
    const float4 * aabb,
    typename SplatPrimitive::WorldBuffer v_splats_world,
    typename SplatPrimitive::ScreenBuffer v_splats_screen,
    GradQuantBuffers gq,
    const uint8_t* sh_value_packed,
    const float2* sh_value_bounds,
    const int64_t sh_value_bounds_stride,
    const int sh_value_bits
);


// Identity permutation + per-gaussian intersection ranges (copies of the FPBO
// helpers; static so they don't clash with FusedProjectionBwdOptim.cu's).
__global__ static void qg_iota_kernel(int64_t n, int32_t* __restrict__ buf) {
    int64_t i = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = (int32_t)i;
}

__global__ static void qg_camera_id_bounds_kernel(
    int64_t nnz, int64_t N,
    const int32_t* __restrict__ gaussian_ids,   // [nnz], sorted by gid
    int32_t* __restrict__ camera_id_bounds       // [N+1]
) {
    int64_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid > nnz) return;
    int32_t cid   = (gid == nnz) ? N : gaussian_ids[gid];
    int32_t cid_0 = (gid == 0) ? 0 : gaussian_ids[gid-1] + 1;
    for (int32_t i = cid_0; i <= cid; ++i)
        camera_id_bounds[i] = gid;
}


template<typename SplatPrimitive>
static inline void launch_projection_bwd_quantgrad(
    int64_t N,
    const uint32_t num_sh_buffer,
    const std::vector<DeviceTensorFloatND>& splats_world_tuple,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const CameraModelType camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    const DeviceTensor2D<float4> aabb,
    const std::vector<DeviceTensorFloatND>& v_splats_screen_tuple,
    const std::vector<DeviceTensorFloatND>& v_splats_world_tuple,
    GradQuantBuffers gq,
    const uint8_t* sh_value_packed,
    const float2* sh_value_bounds,
    const int64_t sh_value_bounds_stride,
    const int sh_value_bits
) {
    typename SplatPrimitive::WorldBuffer  splats_world(splats_world_tuple);
    typename SplatPrimitive::WorldBuffer  v_splats_world(v_splats_world_tuple);
    typename SplatPrimitive::ScreenBuffer v_splats_screen(v_splats_screen_tuple);

    uint32_t C = (uint32_t)std::get<2>(viewmats)[0];
    const float*  viewmats_ptr = (const float*)std::get<0>(viewmats);
    const float4* intrins_ptr  = (const float4*)std::get<0>(intrins);

    const bool is_packed = (camera_ids.data_ptr() && gaussian_ids.data_ptr());
    // NOTE: this kernel is SPLAT-parallel (one thread per splat, FPBO-style),
    // NOT intersection-parallel -- so N stays = num_splats even when packed
    // (unlike the C*N atomic projection backward). `nnz` (intersection count)
    // is used only to size/sort the per-splat camera_id_bounds below.
    if (N == 0) return;
    if (!is_packed && aabb.data_ptr() == nullptr) return;

    // Build per-gaussian intersection ranges for the packed camera loop, exactly
    // as the FPBO launcher does (sort a permutation alongside gaussian_ids so the
    // kernel can recover original out_idx for aabb / screen-grad loads).
    DeviceVector<int32_t> camera_id_bounds;
    DeviceVector<int32_t> gaussian_ids_sorted_buf, perm_buf, perm_sorted_buf;
    const int32_t* sorted_gaussian_ids_ptr = nullptr;
    const int32_t* sorted_perm_ptr = nullptr;
    if (is_packed) {
        long nnz = camera_ids.size();
        gaussian_ids_sorted_buf.resize(PoolSlot::FusedProjBwdGaussSorted, nnz);
        perm_buf.resize(PoolSlot::FusedProjBwdPerm, nnz);
        perm_sorted_buf.resize(PoolSlot::FusedProjBwdPermSorted, nnz);
        qg_iota_kernel<<<_LAUNCH_ARGS_1D(nnz, 256)>>>(nnz, perm_buf.data_ptr());
        CHECK_DEVICE_ERROR(cudaGetLastError());
        cub::DoubleBuffer<int32_t> d_keys(gaussian_ids.data_ptr(), gaussian_ids_sorted_buf.data_ptr());
        cub::DoubleBuffer<int32_t> d_values(perm_buf.data_ptr(), perm_sorted_buf.data_ptr());
        int n_bits = 0;
        while ((1U << n_bits) <= N) ++n_bits;
        CUB_WRAPPER(cub::DeviceRadixSort::SortPairs, d_keys, d_values, nnz, 0, n_bits);
        CHECK_DEVICE_ERROR(cudaGetLastError());
        sorted_gaussian_ids_ptr = d_keys.selector ? gaussian_ids_sorted_buf.data_ptr() : gaussian_ids.data_ptr();
        sorted_perm_ptr         = d_values.selector ? perm_sorted_buf.data_ptr() : perm_buf.data_ptr();
        camera_id_bounds.resize(PoolSlot::FusedProjBwdCamBounds, (int64_t)(N+1));
        qg_camera_id_bounds_kernel<<<_LAUNCH_ARGS_1D(nnz+1, 256)>>>(
            nnz, N, sorted_gaussian_ids_ptr, camera_id_bounds.data_ptr());
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    #define _LAUNCH_ARGS ( \
            (cudaStream_t)0, C, (uint32_t)N, num_sh_buffer, \
            splats_world, viewmats_ptr, intrins_ptr, dist_coeffs, \
            image_width, image_height, \
            is_packed ? camera_id_bounds.data_ptr() : nullptr, \
            is_packed ? camera_ids.data_ptr() : nullptr, \
            is_packed ? sorted_perm_ptr : nullptr, \
            aabb.data_ptr(), \
            v_splats_world, v_splats_screen, gq, \
            sh_value_packed, sh_value_bounds, sh_value_bounds_stride, sh_value_bits )

    if (camera_model == CameraModelType::PINHOLE)
        projection_bwd_quantgrad_kernel_wrapper<SplatPrimitive, CameraModelType::PINHOLE> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::FISHEYE)
        projection_bwd_quantgrad_kernel_wrapper<SplatPrimitive, CameraModelType::FISHEYE> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::EQUISOLID)
        projection_bwd_quantgrad_kernel_wrapper<SplatPrimitive, CameraModelType::EQUISOLID> _LAUNCH_ARGS;
    else if (camera_model == CameraModelType::EQUIRECTANGULAR)
        projection_bwd_quantgrad_kernel_wrapper<SplatPrimitive, CameraModelType::EQUIRECTANGULAR> _LAUNCH_ARGS;
    else
        throw std::runtime_error("Unsupported camera model");
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #undef _LAUNCH_ARGS
}


// Shared per-primitive dispatch: cap sh_degree at the runtime warmup degree and
// pick the PrimT<sh_degree> specialization.
template<template<int> class PrimT>
static inline void _projection_bwd_quantgrad_dispatch(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND>& splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    const DeviceTensor2D<float4> aabb,
    const std::vector<DeviceTensorFloatND>& v_splats_screen,
    const std::vector<DeviceTensorFloatND>& v_splats_world,
    GradQuantBuffers gq,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_value_bounds_stride
) {
    int sh_degree = typename PrimT<0>::WorldBuffer(splats_world).sh_degree();
    sh_degree = min(sh_degree, max_sh_degree);
    const uint8_t* vp = sh_value_packed.has_value() ? (const uint8_t*)std::get<0>(sh_value_packed.value()) : nullptr;
    const float2*  vb = sh_value_bounds.has_value() ? (const float2*)std::get<0>(sh_value_bounds.value()) : nullptr;
    #define LAUNCH(n) if (sh_degree == (n)) \
        return launch_projection_bwd_quantgrad<PrimT<n>>( \
            num_splats, num_sh_buffer, splats_world, viewmats, intrins, \
            image_width, image_height, cmt(camera_model), dist_coeffs, \
            camera_ids, gaussian_ids, aabb, v_splats_screen, v_splats_world, \
            gq, vp, vb, sh_value_bounds_stride, sh_value_bits);
    LAUNCH(3) LAUNCH(2) LAUNCH(1) LAUNCH(0) LAUNCH(4)
    #undef LAUNCH
}


/*[AutoHeaderGeneratorExport]*/
void projection_3dgs_backward_quantgrad(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND>& splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    const DeviceTensor2D<float4> aabb,
    const std::vector<DeviceTensorFloatND>& v_splats_screen,
    const std::vector<DeviceTensorFloatND>& v_splats_world,
    GradQuantBuffers gq,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_value_bounds_stride
) {
    _projection_bwd_quantgrad_dispatch<Vanilla3DGS>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins,
        image_width, image_height, camera_model, dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, v_splats_world,
        gq, sh_value_packed, sh_value_bounds, num_sh_buffer, sh_value_bits, sh_value_bounds_stride);
}


/*[AutoHeaderGeneratorExport]*/
void projection_mip_backward_quantgrad(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND>& splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    const DeviceTensor2D<float4> aabb,
    const std::vector<DeviceTensorFloatND>& v_splats_screen,
    const std::vector<DeviceTensorFloatND>& v_splats_world,
    GradQuantBuffers gq,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_value_bounds_stride
) {
    _projection_bwd_quantgrad_dispatch<MipSplatting>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins,
        image_width, image_height, camera_model, dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, v_splats_world,
        gq, sh_value_packed, sh_value_bounds, num_sh_buffer, sh_value_bits, sh_value_bounds_stride);
}


/*[AutoHeaderGeneratorExport]*/
void projection_3dgut_backward_quantgrad(
    const int64_t num_splats,
    const int max_sh_degree,
    const std::vector<DeviceTensorFloatND>& splats_world,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const TorchTensorView dist_coeffs,
    const DeviceVector<int32_t> camera_ids,
    const DeviceVector<int32_t> gaussian_ids,
    const DeviceTensor2D<float4> aabb,
    const std::vector<DeviceTensorFloatND>& v_splats_screen,
    const std::vector<DeviceTensorFloatND>& v_splats_world,
    GradQuantBuffers gq,
    const std::optional<TorchTensorView> sh_value_packed,
    const std::optional<TorchTensorView> sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_value_bounds_stride
) {
    _projection_bwd_quantgrad_dispatch<Vanilla3DGUT>(
        num_splats, max_sh_degree, splats_world, viewmats, intrins,
        image_width, image_height, camera_model, dist_coeffs,
        camera_ids, gaussian_ids, aabb, v_splats_screen, v_splats_world,
        gq, sh_value_packed, sh_value_bounds, num_sh_buffer, sh_value_bits, sh_value_bounds_stride);
}
