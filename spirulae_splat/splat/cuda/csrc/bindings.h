#include "cuda_runtime.h"
#include <cstdio>
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



float4 tuple2float4(std::tuple<float, float, float, float> v);


torch::Tensor compute_sh_forward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs0,   // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
);


std::tuple<torch::Tensor, torch::Tensor>
compute_sh_backward_tensor(
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
    torch::Tensor  // axes_v, [N, 4]
    // torch::Tensor  // depth_grads, [N, 2]
> project_gaussians_forward_tensor(
    const int num_points,
    torch::Tensor &means3d,  // [N, 3]
    torch::Tensor &scales,  // [N, 2]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &viewmat,  // [3|4, 4]
    std::string camera_model,
    std::tuple<float, float, float, float> intrins,
    std::tuple<float, float, float, float> dist_coeffs,
    const unsigned img_height,
    const unsigned img_width,
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
    std::tuple<float, float, float, float> intrins,
    torch::Tensor &num_tiles_hit,  // [N], int
    torch::Tensor &v_positions,  // [N, 3]
    torch::Tensor &v_axes_u,  // [N, 3]
    torch::Tensor &v_axes_v  // [N, 3]
    // torch::Tensor &v_depth_grads  // [N, 2]
);


torch::Tensor render_undistortion_map_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::tuple<float, float, float, float> dist_coeffs
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
    const unsigned img_height,
    const unsigned img_width
);


torch::Tensor get_tile_bin_edges_tensor(
    int num_intersects,
    const torch::Tensor &isect_ids_sorted, 
    const unsigned img_height,
    const unsigned img_width
);


std::tuple<
    torch::Tensor,  // new_opacities
    torch::Tensor  // new_scales
> compute_relocation_tensor(
    torch::Tensor &opacities,
    torch::Tensor &scales,
    torch::Tensor &ratios,
    torch::Tensor &binoms,
    const int n_max
);


std::tuple<
    torch::Tensor,  // new_position_offsets
    torch::Tensor,  // new_opacities
    torch::Tensor  // new_scales
> compute_relocation_split_tensor(
    torch::Tensor &positions,
    torch::Tensor &quats,
    torch::Tensor &opacities,
    torch::Tensor &scales
);


std::tuple<
    torch::Tensor,  // final_index
    torch::Tensor,  // out_img
    torch::Tensor  // out_alpha
> rasterize_simple_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
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
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor // v_opacities
> rasterize_simple_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
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
    torch::Tensor,  // final_index
    torch::Tensor,  // out_img
    torch::Tensor  // out_visibility
> rasterize_depth_forward_tensor(
    const std::string depth_mode,
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &gaussian_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &opacities
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor // v_opacities
> rasterize_depth_backward_tensor(
    const std::string depth_mode,
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &opacities,
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
    torch::Tensor,  // out_normal
    torch::Tensor  // out_reg_depth
> rasterize_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const float depth_reg_pairwise_factor,
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
    // const torch::Tensor &background,
    const torch::Tensor &depth_ref_im
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_ch_coeffs
    // torch::Tensor, // v_ch_coeffs_abs
    torch::Tensor, // v_opacities
    // torch::Tensor, // v_background
    torch::Tensor  // v_depth_ref_im
> rasterize_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float depth_reg_pairwise_factor,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &ch_coeffs,
    const torch::Tensor &opacities,
    // const torch::Tensor &background,
    const torch::Tensor &depth_ref_im,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_depth,
    const torch::Tensor &v_output_normal,
    const torch::Tensor &v_output_reg_depth
);


std::tuple<
    torch::Tensor, // final_idx
    torch::Tensor,  // out_alpha
    torch::Tensor,  // out_img
    torch::Tensor,  // out_depth
    torch::Tensor,  // out_normal
    torch::Tensor  // out_depth_reg
> rasterize_simplified_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &gaussian_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor // v_opacities
> rasterize_simplified_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &final_idx,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output_img,
    const torch::Tensor &v_output_depth,
    const torch::Tensor &v_output_normal,
    const torch::Tensor &v_output_depth_reg
);


std::tuple<
    torch::Tensor,  // num_intersects
    torch::Tensor,  // sorted_indices
    torch::Tensor  // sorted_depths
> rasterize_indices_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &gaussian_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &opacities
);


void sort_per_pixel_tensor(
    const std::string &method,
    const unsigned num_pixels,
    torch::Tensor &num_intersects,  // [h, w]
    torch::Tensor &indices,  // [h, w, MAX_SORTED_SPLATS]
    torch::Tensor &depths  // [h, w, MAX_SORTED_SPLATS]
);


std::tuple<
    torch::Tensor,  // out_img
    torch::Tensor  // out_alpha
> rasterize_simple_sorted_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor // v_opacities
> rasterize_simple_sorted_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const torch::Tensor &num_intersects,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &background,
    const torch::Tensor &output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_alpha
);


std::tuple<
    torch::Tensor,  // final_idx
    torch::Tensor,  // out_img
    torch::Tensor  // out_visibility
> rasterize_depth_sorted_forward_tensor(
    const std::string depth_mode,
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &opacities
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor // v_opacities
> rasterize_depth_sorted_backward_tensor(
    const std::string depth_mode,
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const torch::Tensor &final_idx,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &opacities,
    const torch::Tensor &output_depth,
    const torch::Tensor &output_visibility,
    const torch::Tensor &v_output_depth
);


std::tuple<
    torch::Tensor,  // out_alpha
    torch::Tensor,  // out_img
    torch::Tensor,  // out_depth
    torch::Tensor,  // out_normal
    torch::Tensor  // out_reg_depth
> rasterize_sorted_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const float depth_reg_pairwise_factor,
    const torch::Tensor &sorted_indices,
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
    // const torch::Tensor &background,
    const torch::Tensor &depth_ref_im
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_ch_coeffs
    // torch::Tensor, // v_ch_coeffs_abs
    torch::Tensor, // v_opacities
    // torch::Tensor, // v_background
    torch::Tensor  // v_depth_ref_im
> rasterize_sorted_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::tuple<float, float, float, float> intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float depth_reg_pairwise_factor,
    const torch::Tensor &num_intersects,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &ch_coeffs,
    const torch::Tensor &opacities,
    // const torch::Tensor &background,
    const torch::Tensor &depth_ref_im,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output,
    const torch::Tensor &v_output_depth,
    const torch::Tensor &v_output_normal,
    const torch::Tensor &v_output_reg_depth
);


std::tuple<
    torch::Tensor,  // out_alpha
    torch::Tensor,  // out_img
    torch::Tensor,  // out_depth
    torch::Tensor,  // out_normal
    torch::Tensor  // out_depth_reg
> rasterize_simplified_sorted_forward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities
);


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor // v_opacities
> rasterize_simplified_sorted_backward_tensor(
    const unsigned img_height,
    const unsigned img_width,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &num_intersects,
    const torch::Tensor &sorted_indices,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &opacities,
    const torch::Tensor &output_alpha,
    const torch::Tensor &output_depth,
    const torch::Tensor &v_output_alpha,
    const torch::Tensor &v_output_img,
    const torch::Tensor &v_output_depth,
    const torch::Tensor &v_output_normal,
    const torch::Tensor &v_output_depth_reg
);


torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs
);


std::tuple<
    torch::Tensor,  // v_rotation
    torch::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs,
    const torch::Tensor &out_color,
    const torch::Tensor &v_out_color
);
