#include "cuda_runtime.h"
#include "forward.cuh"
#include <cstdio>
#include <iostream>
#include <math.h>
#include <torch/extension.h>
#include <tuple>
#include <c10/cuda/CUDAGuard.h>

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x)                                                    \
    TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x)                                                         \
    CHECK_CUDA(x);                                                             \
    CHECK_CONTIGUOUS(x)
#define DEVICE_GUARD(_ten) \
    const at::cuda::OptionalCUDAGuard device_guard(device_of(_ten));

std::tuple<
    torch::Tensor, // output conics
    torch::Tensor> // output radii
compute_cov2d_bounds_tensor(const int num_pts, torch::Tensor &A);

torch::Tensor compute_sh_forward_tensor(
    const std::string &method,
    unsigned num_points,
    unsigned degree,
    unsigned degrees_to_use,
    torch::Tensor &viewdirs,
    torch::Tensor &coeffs
);

torch::Tensor compute_sh_backward_tensor(
    const std::string &method,
    unsigned num_points,
    unsigned degree,
    unsigned degrees_to_use,
    torch::Tensor &viewdirs,
    torch::Tensor &v_colors
);

std::tuple<
    torch::Tensor,  // bounds
    torch::Tensor,  // num_tiles_hit
    torch::Tensor,  // positions
    torch::Tensor,  // axes_u
    torch::Tensor,  // axes_v
    torch::Tensor  // depth_grads
> project_gaussians_forward_tensor(
    const int num_points,
    torch::Tensor &means3d,
    torch::Tensor &scales,
    torch::Tensor &quats,
    torch::Tensor &viewmat,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const unsigned img_height,
    const unsigned img_width,
    const unsigned block_width,
    const float clip_thresh
);

std::tuple<
    torch::Tensor,  // v_means3d
    torch::Tensor,  // v_scales
    torch::Tensor  // v_quats
> project_gaussians_backward_tensor(
    const int num_points,
    torch::Tensor &means3d,
    torch::Tensor &scales,
    torch::Tensor &quats,
    torch::Tensor &viewmat,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    torch::Tensor &num_tiles_hit,
    torch::Tensor &v_positions,
    torch::Tensor &v_axes_u,
    torch::Tensor &v_axes_v,
    torch::Tensor &v_depth_grads
);


std::tuple<
    torch::Tensor, torch::Tensor
> map_gaussian_to_intersects_tensor(
    const int num_points,
    const int num_intersects,
    const torch::Tensor &positions,
    const torch::Tensor &bounds,
    const torch::Tensor &cum_tiles_hit,
    const std::tuple<int, int, int> tile_bounds,
    const unsigned block_width
);

torch::Tensor get_tile_bin_edges_tensor(
    int num_intersects,
    const torch::Tensor &isect_ids_sorted,
    const std::tuple<int, int, int> tile_bounds
);

std::tuple<
    torch::Tensor,  // final_index
    torch::Tensor,  // out_img
    torch::Tensor  // out_alpha
> rasterize_simple_forward_tensor(
    const std::tuple<int, int, int> tile_bounds,
    const std::tuple<int, int, int> block,
    const std::tuple<int, int, int> img_size,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &gaussian_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background
);

std::tuple<
    torch::Tensor, // final_idx
    torch::Tensor,  // out_alpha
    torch::Tensor,  // out_img
    torch::Tensor,  // out_depth
    torch::Tensor,  // out_reg_depth
    torch::Tensor  // out_reg_normal
> rasterize_forward_tensor(
    const std::tuple<int, int, int> tile_bounds,
    const std::tuple<int, int, int> block,
    const std::tuple<int, int, int> img_size,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &gaussian_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background,
    const torch::Tensor &depth_grads,
    const torch::Tensor &depth_normal_ref_im
);

std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor  // v_opacities
> rasterize_simple_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const unsigned block_width,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_alpha
);

std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_opacities
    torch::Tensor, // v_depth_grad
    torch::Tensor  // v_depth_normal_ref
> rasterize_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const unsigned block_width,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background,
    const torch::Tensor &depth_grads,
    const torch::Tensor &depth_normal_ref_im,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth_grad,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_depth_grad,
    const torch::Tensor &v_output_reg_depth,
    const torch::Tensor &v_output_reg_normal
);
