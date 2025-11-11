#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include "Primitive3DGS.cuh"
#include "PrimitiveOpaqueTriangle.cuh"

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    Vanilla3DGS::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<Vanilla3DGS::RenderOutput::TensorTuple>,
    std::optional<Vanilla3DGS::RenderOutput::TensorTuple>
> rasterize_to_pixels_3dgs_eval3d_fwd(
    // Gaussian parameters
    Vanilla3DGS::WorldEval3D::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
);


std::tuple<
    OpaqueTriangle::RenderOutput::TensorTuple,
    at::Tensor,
    at::Tensor,
    std::optional<OpaqueTriangle::RenderOutput::TensorTuple>,
    std::optional<OpaqueTriangle::RenderOutput::TensorTuple>
> rasterize_to_pixels_opaque_triangle_eval3d_fwd(
    // Gaussian parameters
    OpaqueTriangle::WorldEval3D::TensorTuple splats_tuple,
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const gsplat::CameraModelType camera_model,
    const CameraDistortionCoeffsTensor dist_coeffs,
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids   // [n_isects]
);
