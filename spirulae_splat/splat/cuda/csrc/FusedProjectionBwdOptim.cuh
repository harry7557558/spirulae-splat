#pragma once

#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "Primitive3DGS.cuh"
#include "Primitive3DGUT.cuh"
#include "Primitive3DGUT_SV.cuh"
#include "PrimitiveOpaqueTriangle.cuh"
#include "PrimitiveVoxel.cuh"

#include "types.cuh"
#include "common.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



// void projection_3dgs_backward_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,  // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_3dgs_backward_with_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const TensorList &vr_splats_world,
//     const TensorList &h_splats_world
// );


// void projection_3dgs_backward_with_position_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const at::Tensor &vr_splats_world,
//     const at::Tensor &h_splats_world
// );


// void projection_mip_backward_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_mip_backward_with_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const TensorList &vr_splats_world,
//     const TensorList &h_splats_world
// );


// void projection_mip_backward_with_position_hessian_diagonal_tensor(
//     // fwd inputs
//     const TensorList &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const TensorList &v_splats_screen,
//     const TensorList &vr_splats_screen,
//     const TensorList &h_splats_screen,
//     // returns
//     const TensorList &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats,
//     const at::Tensor &vr_splats_world,
//     const at::Tensor &h_splats_world
// );


void fused_projection_bwd_optimizer_3dgut_tensor(
    // fwd inputs
    const int64_t num_splats,
    TensorList splats_world,
    const at::Tensor viewmats,
    const at::Tensor intrins,
    const uint32_t image_width,
    const uint32_t image_height,
    const std::string camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,
    const std::optional<at::Tensor> gaussian_ids,
    const at::Tensor aabb,
    // grad outputs
    const TensorList v_splats_world,
    const std::optional<TensorList> vr_splats_world,
    const std::optional<TensorList> h_splats_world,
    const TensorList v_splats_screen,
    const std::optional<TensorList> vr_splats_screen,
    const std::optional<TensorList> h_splats_screen,
    // optimizer states
    const TensorList g1_splats_world,
    const TensorList g2_splats_world,
    const std::optional<at::Tensor> g1_features_sh,
    const std::optional<at::Tensor> g2_features_sh,
    const std::optional<at::Tensor> sh_quant_bounds,
    // optimizer params
    const at::Tensor radii,
    const float lr_means,
    const float lr_quats,
    const float lr_scales,
    const float lr_opacs,
    const float lr_features_dc,
    const float lr_features_sh,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    std::variant<int32_t, at::Tensor> step
);


// void projection_3dgut_sv_backward_tensor(
//     // fwd inputs
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const SphericalVoronoi3DGUT_Default::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_opaque_triangle_backward_tensor(
//     // fwd inputs
//     const OpaqueTriangle::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const OpaqueTriangle::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const OpaqueTriangle::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );


// void projection_voxel_backward_tensor(
//     // fwd inputs
//     const VoxelPrimitive::World::TensorTuple &splats_world,
//     const at::Tensor viewmats,  // [..., C, 4, 4]
//     const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
//     const uint32_t image_width,
//     const uint32_t image_height,
//     const std::string camera_model,
//     const CameraDistortionCoeffsTensor dist_coeffs,
//     // fwd outputs
//     const std::optional<at::Tensor> camera_ids,  // [nnz]
//     const std::optional<at::Tensor> gaussian_ids,  // [nnz]
//     const at::Tensor aabb,                       // [..., C, N, 2]
//     // grad outputs
//     const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
//     // returns
//     const VoxelPrimitive::World::TensorTuple &v_splats_world,
//     const std::optional<at::Tensor> &v_viewmats
// );
