#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>

#include <gsplat/Common.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



std::tuple<
    at::Tensor,  // camera_ids
    at::Tensor,  // gaussian_ids
    at::Tensor,  // radii
    at::Tensor,  // means2d
    at::Tensor,  // depths
    at::Tensor,  // conics
    at::Tensor  // compensations
> projection_ewa_3dgs_hetero_forward_tensor(
    // inputs
    const at::Tensor means,                // [..., N, 3]
    const at::Tensor quats,  // [..., N, 4]
    const at::Tensor scales, // [..., N, 3]
    const at::optional<at::Tensor> opacities, // [..., N]
    const at::Tensor viewmats,             // [..., C, 4, 4]
    const at::Tensor Ks,                   // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const float near_plane,
    const float far_plane,
    const float radius_clip,
    const bool calc_compensations,
    const gsplat::CameraModelType camera_model,
    const at::Tensor intersection_count_map,  // [C+1]
    const at::Tensor intersection_splat_id  // [nnz]
);


std::tuple<
    at::Tensor,  // v_means
    at::Tensor,  // v_quats
    at::Tensor,  // v_scales
    at::Tensor  // v_viewmats
> projection_ewa_3dgs_hetero_backward_tensor(
    // fwd inputs
    const at::Tensor means, // [..., N, 3]
    const at::Tensor quats, // [..., N, 4]
    const at::Tensor scales, // [..., N, 3]
    const at::Tensor viewmats, // [..., C, 4, 4]
    const at::Tensor Ks, // [..., C, 3, 3]
    const uint32_t image_width,
    const uint32_t image_height,
    const float eps2d,
    const gsplat::CameraModelType camera_model,
    // fwd outputs
    const at::Tensor camera_ids, // [nnz]
    const at::Tensor gaussian_ids, // [nnz]
    const at::Tensor conics, // [nnz, 3]
    const at::optional<at::Tensor> compensations, // [nnz] optional
    // grad outputs
    const at::Tensor v_means2d, // [nnz, 2]
    const at::Tensor v_depths, // [nnz]
    const at::Tensor v_conics, // [nnz, 3]
    const at::optional<at::Tensor> v_compensations, // [nnz] optional
    const bool viewmats_requires_grad,
    const bool sparse_grad
);
