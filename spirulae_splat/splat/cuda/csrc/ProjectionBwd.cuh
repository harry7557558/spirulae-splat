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



std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgs_backward_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,  // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgs_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    Vanilla3DGS::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGS::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgs_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGS::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGS::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGS::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_mip_backward_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>  // h_world_pos or h_splats
> projection_mip_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &vr_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    MipSplatting::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, MipSplatting::World::TensorTuple>  // h_world_pos or h_splats
> projection_mip_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const MipSplatting::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const MipSplatting::Screen::TensorTupleProj &v_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &vr_splats_screen,
    const MipSplatting::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_backward_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgut_backward_with_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    Vanilla3DGUT::World::TensorTuple,  // v_splats
    at::Tensor,  // v_viewmats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>,  // vr_world_pos or vr_splats
    std::variant<at::Tensor, Vanilla3DGUT::World::TensorTuple>  // h_world_pos or h_splats
> projection_3dgut_backward_with_position_hessian_diagonal_tensor(
    // fwd inputs
    const Vanilla3DGUT::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const Vanilla3DGUT::Screen::TensorTupleProj &v_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &vr_splats_screen,
    const Vanilla3DGUT::Screen::TensorTupleProj &h_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    SphericalVoronoi3DGUT_Default::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_3dgut_sv_backward_tensor(
    // fwd inputs
    const SphericalVoronoi3DGUT_Default::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const SphericalVoronoi3DGUT_Default::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    OpaqueTriangle::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_opaque_triangle_backward_tensor(
    // fwd inputs
    const OpaqueTriangle::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const OpaqueTriangle::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);


std::tuple<
    VoxelPrimitive::World::TensorTuple,  // v_splats
    at::Tensor  // v_viewmats
> projection_voxel_backward_tensor(
    // fwd inputs
    const VoxelPrimitive::World::TensorTuple &splats_world,
    const at::Tensor viewmats,  // [..., C, 4, 4]
    const at::Tensor intrins,  // [..., C, 4], fx, fy, cx, cy
    const uint32_t image_width,
    const uint32_t image_height,
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    // fwd outputs
    const std::optional<at::Tensor> camera_ids,  // [nnz]
    const std::optional<at::Tensor> gaussian_ids,  // [nnz]
    const at::Tensor aabb,                       // [..., C, N, 2]
    // grad outputs
    const VoxelPrimitive::Screen::TensorTupleProj &v_splats_screen,
    const bool viewmats_requires_grad
);
