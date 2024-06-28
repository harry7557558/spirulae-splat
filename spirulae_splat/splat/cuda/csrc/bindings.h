#include "cuda_runtime.h"
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


/* == AUTO HEADER GENERATOR - DO NOT CHANGE THIS LINE == */



torch::Tensor compute_sh_forward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
);


torch::Tensor compute_sh_backward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &v_colors  // [..., 3]
);


std::tuple<
    torch::Tensor,  // bounds, [N, 4], int
    torch::Tensor,  // num_tiles_hit, [N]
    torch::Tensor,  // positions, [N, 3]
    torch::Tensor,  // axes_u, [N, 3]
    torch::Tensor,  // axes_v, [N, 4]
    torch::Tensor  // depth_gradsm, [N, 2]
> project_gaussians_forward_tensor(
    const int num_points,
    torch::Tensor &means3d,  // [N, 3]
    torch::Tensor &scales,  // [N, 2]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &viewmat,  // [3|4, 4]
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
    torch::Tensor,  // v_quats
    torch::Tensor  // v_viewmat, [4, 4]
> project_gaussians_backward_tensor(
    const int num_points,
    torch::Tensor &means3d,  // [N, 3]
    torch::Tensor &scales,  // [N, 2]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &viewmat,  // [3|4, 4]
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    torch::Tensor &num_tiles_hit,  // [N], int
    torch::Tensor &v_positions,  // [N, 3]
    torch::Tensor &v_axes_u,  // [N, 3]
    torch::Tensor &v_axes_v,  // [N, 3]
    torch::Tensor &v_depth_grads  // [N, 2]
);


std::tuple<
    torch::Tensor,
    torch::Tensor
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
    int num_intersects, const torch::Tensor &isect_ids_sorted, 
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
    const torch::Tensor &anisotropies,
    const torch::Tensor &background
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_opacities
    torch::Tensor  // v_anisotropies
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
    const torch::Tensor &anisotropies,
    const torch::Tensor &background,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_alpha
);


std::tuple<
    torch::Tensor,  // final_index
    torch::Tensor,  // out_img
    torch::Tensor  // out_visibility
> rasterize_depth_forward_tensor(
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
    const torch::Tensor &opacities,
    const torch::Tensor &anisotropies
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_opacities
    torch::Tensor  // v_anisotropies
> rasterize_depth_backward_tensor(
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
    const torch::Tensor &opacities,
    const torch::Tensor &anisotropies,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_depth,
    const torch::Tensor &output_visibility,
    const torch::Tensor &v_output_depth
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
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const torch::Tensor &ch_coeffs,
    const torch::Tensor &opacities,
    const torch::Tensor &anisotropies,
    const torch::Tensor &background,
    const torch::Tensor &depth_grads,
    const torch::Tensor &depth_ref_im
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_ch_coeffs
    torch::Tensor, // v_ch_coeffs_abs
    torch::Tensor, // v_opacities
    torch::Tensor, // v_anisotropies
    torch::Tensor, // v_depth_grad
    torch::Tensor  // v_depth_ref_im
> rasterize_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const unsigned block_width,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &ch_coeffs,
    const torch::Tensor &opacities,
    const torch::Tensor &anisotropies,
    const torch::Tensor &background,
    const torch::Tensor &depth_grads,
    const torch::Tensor &depth_ref_im,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth_grad,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_depth_grad,
    const torch::Tensor &v_output_reg_depth,
    const torch::Tensor &v_output_reg_normal
);
