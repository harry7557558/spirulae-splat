#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"
#include "common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



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
);


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
);


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
);


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
);


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
);


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
);


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
);


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
);


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
);


// void projection_3dgut_sv_backward(
//     // fwd inputs
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
//     TorchTensorView viewmats,  // [C, 4, 4]
//     TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const TorchTensorView dist_coeffs,
//     // fwd outputs
//     const std::optional<DeviceVector<int32_t>> camera_ids,  // [nnz]
//     const std::optional<DeviceVector<int32_t>> gaussian_ids,  // [nnz]
//     const DeviceTensor2D<float4> aabb,  // [C, N]
//     // grad outputs
//     const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_opaque_triangle_backward(
//     // fwd inputs
//     const OpaqueTriangle::World::TensorTuple &splats_world,
//     TorchTensorView viewmats,  // [C, 4, 4]
//     TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const TorchTensorView dist_coeffs,
//     // fwd outputs
//     const std::optional<DeviceVector<int32_t>> camera_ids,  // [nnz]
//     const std::optional<DeviceVector<int32_t>> gaussian_ids,  // [nnz]
//     const DeviceTensor2D<float4> aabb,  // [C, N]
//     // grad outputs
//     const OpaqueTriangle::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const OpaqueTriangle::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_voxel_backward(
//     // fwd inputs
//     const VoxelPrimitive::World::TensorTuple &splats_world,
//     TorchTensorView viewmats,  // [C, 4, 4]
//     TorchTensorView intrins,   // [C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const TorchTensorView dist_coeffs,
//     // fwd outputs
//     const std::optional<DeviceVector<int32_t>> camera_ids,  // [nnz]
//     const std::optional<DeviceVector<int32_t>> gaussian_ids,  // [nnz]
//     const DeviceTensor2D<float4> aabb,  // [C, N]
//     // grad outputs
//     const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const VoxelPrimitive::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );
