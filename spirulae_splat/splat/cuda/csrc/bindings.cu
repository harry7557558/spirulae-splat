#include "bindings.h"
#include "helpers.cuh"
#include "projection.cuh"
#include "rasterization.cuh"
#include "sh.cuh"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <math.h>
#include <torch/extension.h>
#include <tuple>

namespace cg = cooperative_groups;



torch::Tensor compute_sh_forward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
) {
    DEVICE_GUARD(viewdirs);
    unsigned num_bases = num_sh_bases(degree);
    if (coeffs.ndimension() != 3 || coeffs.size(0) != num_points ||
        coeffs.size(1) != num_bases || coeffs.size(2) != 3) {
        AT_ERROR("coeffs must have dimensions (N, D, 3)");
    }
    torch::Tensor colors = torch::empty({num_points, 3}, coeffs.options());    
    if (method == "poly") {
        compute_sh_forward_kernel<SHType::Poly><<<
            (num_points + N_THREADS - 1) / N_THREADS,
            N_THREADS>>>(
            num_points,
            degree,
            degrees_to_use,
            (float3 *)viewdirs.contiguous().data_ptr<float>(),
            coeffs.contiguous().data_ptr<float>(),
            colors.contiguous().data_ptr<float>()
        );
    } else if (method == "fast") {
        compute_sh_forward_kernel<SHType::Fast><<<
            (num_points + N_THREADS - 1) / N_THREADS,
            N_THREADS>>>(
            num_points,
            degree,
            degrees_to_use,
            (float3 *)viewdirs.contiguous().data_ptr<float>(),
            coeffs.contiguous().data_ptr<float>(),
            colors.contiguous().data_ptr<float>()
        );
    } else {
        AT_ERROR("Invalid method: ", method);
    }
    return colors;
}


torch::Tensor compute_sh_backward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &v_colors  // [..., 3]
) {
    DEVICE_GUARD(viewdirs);
    if (viewdirs.ndimension() != 2 || viewdirs.size(0) != num_points ||
        viewdirs.size(1) != 3) {
        AT_ERROR("viewdirs must have dimensions (N, 3)");
    }
    if (v_colors.ndimension() != 2 || v_colors.size(0) != num_points ||
        v_colors.size(1) != 3) {
        AT_ERROR("v_colors must have dimensions (N, 3)");
    }
    unsigned num_bases = num_sh_bases(degree);
    torch::Tensor v_coeffs =
        torch::zeros({num_points, num_bases, 3}, v_colors.options());
    if (method == "poly") {
        compute_sh_backward_kernel<SHType::Poly><<<
            (num_points + N_THREADS - 1) / N_THREADS,
            N_THREADS>>>(
            num_points,
            degree,
            degrees_to_use,
            (float3 *)viewdirs.contiguous().data_ptr<float>(),
            v_colors.contiguous().data_ptr<float>(),
            v_coeffs.contiguous().data_ptr<float>()
        );
    } else if (method == "fast") {
        compute_sh_backward_kernel<SHType::Fast><<<
            (num_points + N_THREADS - 1) / N_THREADS,
            N_THREADS>>>(
            num_points,
            degree,
            degrees_to_use,
            (float3 *)viewdirs.contiguous().data_ptr<float>(),
            v_colors.contiguous().data_ptr<float>(),
            v_coeffs.contiguous().data_ptr<float>()
        );
    } else {
        AT_ERROR("Invalid method: ", method);
    }
    return v_coeffs;
}


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
) {
    DEVICE_GUARD(means3d);

    dim3 tile_bounds_dim3;
    tile_bounds_dim3.x = int((img_width + block_width - 1) / block_width);
    tile_bounds_dim3.y = int((img_height + block_width - 1) / block_width);
    tile_bounds_dim3.z = 1;

    float4 intrins = {fx, fy, cx, cy};

    auto int32 = means3d.options().dtype(torch::kInt32);
    auto float32 = means3d.options().dtype(torch::kFloat32);
    torch::Tensor bounds_d = torch::zeros({num_points, 4}, int32);
    torch::Tensor num_tiles_hit_d = torch::zeros({num_points}, int32);
    torch::Tensor positions_d = torch::zeros({num_points, 3}, float32);
    torch::Tensor axes_u_d = torch::zeros({num_points, 3}, float32);
    torch::Tensor axes_v_d = torch::zeros({num_points, 3}, float32);
    torch::Tensor depth_grads_d = torch::zeros({num_points, 2}, float32);

    project_gaussians_forward_kernel<<<
        (num_points + N_THREADS - 1) / N_THREADS,
        N_THREADS>>>(
        num_points,
        (float3 *)means3d.contiguous().data_ptr<float>(),
        (float2 *)scales.contiguous().data_ptr<float>(),
        (float4 *)quats.contiguous().data_ptr<float>(),
        viewmat.contiguous().data_ptr<float>(),
        intrins,
        tile_bounds_dim3,
        block_width,
        clip_thresh,
        // Outputs.
        (int4 *)bounds_d.contiguous().data_ptr<int32_t>(),
        num_tiles_hit_d.contiguous().data_ptr<int32_t>(),
        (float3 *)positions_d.contiguous().data_ptr<float>(),
        (float3 *)axes_u_d.contiguous().data_ptr<float>(),
        (float3 *)axes_v_d.contiguous().data_ptr<float>(),
        (float2 *)depth_grads_d.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        bounds_d, num_tiles_hit_d,
        positions_d, axes_u_d, axes_v_d,
        depth_grads_d
    );
}

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
) {
    DEVICE_GUARD(means3d);

    float4 intrins = {fx, fy, cx, cy};

    auto int32 = means3d.options().dtype(torch::kInt32);
    auto float32 = means3d.options().dtype(torch::kFloat32);
    torch::Tensor v_means3d = torch::zeros({num_points, 3}, float32);
    torch::Tensor v_scales = torch::zeros({num_points, 2}, float32);
    torch::Tensor v_quats = torch::zeros({num_points, 4}, float32);
    torch::Tensor v_viewmat = torch::zeros({4, 4}, float32);

    project_gaussians_backward_kernel<<<
        (num_points + N_THREADS - 1) / N_THREADS,
        N_THREADS>>>(
        num_points,
        (float3 *)means3d.contiguous().data_ptr<float>(),
        (float2 *)scales.contiguous().data_ptr<float>(),
        (float4 *)quats.contiguous().data_ptr<float>(),
        viewmat.contiguous().data_ptr<float>(),
        intrins,
        num_tiles_hit.contiguous().data_ptr<int32_t>(),
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float2 *)v_depth_grads.contiguous().data_ptr<float>(),
        // Outputs.
        (float3 *)v_means3d.contiguous().data_ptr<float>(),
        (float2 *)v_scales.contiguous().data_ptr<float>(),
        (float4 *)v_quats.contiguous().data_ptr<float>(),
        (float *)v_viewmat.contiguous().data_ptr<float>()
    );

    return std::make_tuple(v_means3d, v_scales, v_quats, v_viewmat);
}


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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(bounds);
    CHECK_INPUT(cum_tiles_hit);

    dim3 tile_bounds_dim3;
    tile_bounds_dim3.x = std::get<0>(tile_bounds);
    tile_bounds_dim3.y = std::get<1>(tile_bounds);
    tile_bounds_dim3.z = std::get<2>(tile_bounds);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto int64 = positions.options().dtype(torch::kInt64);
    torch::Tensor gaussian_ids_unsorted =
        torch::zeros({num_intersects}, int32);
    torch::Tensor isect_ids_unsorted =
        torch::zeros({num_intersects}, int64);

    map_gaussian_to_intersects<<<
        (num_points + N_THREADS - 1) / N_THREADS,
        N_THREADS>>>(
        num_points,
        (float3 *)positions.contiguous().data_ptr<float>(),
        (int4 *)bounds.contiguous().data_ptr<int32_t>(),
        cum_tiles_hit.contiguous().data_ptr<int32_t>(),
        tile_bounds_dim3,
        block_width,
        // Outputs.
        isect_ids_unsorted.contiguous().data_ptr<int64_t>(),
        gaussian_ids_unsorted.contiguous().data_ptr<int32_t>()
    );

    return std::make_tuple(isect_ids_unsorted, gaussian_ids_unsorted);
}


torch::Tensor get_tile_bin_edges_tensor(
    int num_intersects, const torch::Tensor &isect_ids_sorted, 
    const std::tuple<int, int, int> tile_bounds
) {
    DEVICE_GUARD(isect_ids_sorted);
    CHECK_INPUT(isect_ids_sorted);
    int num_tiles = std::get<0>(tile_bounds) * std::get<1>(tile_bounds);
    torch::Tensor tile_bins = torch::zeros(
        {num_tiles, 2}, isect_ids_sorted.options().dtype(torch::kInt32)
    );
    get_tile_bin_edges<<<
        (num_intersects + N_THREADS - 1) / N_THREADS,
        N_THREADS>>>(
        num_intersects,
        isect_ids_sorted.contiguous().data_ptr<int64_t>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>()
    );
    return tile_bins;
}



std::tuple<
    torch::Tensor,  // new_opacities
    torch::Tensor  // new_scales
> compute_relocation_tensor(
    torch::Tensor &opacities,
    torch::Tensor &scales,
    torch::Tensor &ratios,
    torch::Tensor &binoms,
    const int n_max
) {
    DEVICE_GUARD(opacities);
    CHECK_INPUT(opacities);
    CHECK_INPUT(scales);
    CHECK_INPUT(ratios);
    CHECK_INPUT(binoms);
    torch::Tensor new_opacities = torch::empty_like(opacities);
    torch::Tensor new_scales = torch::empty_like(scales);

    uint32_t N = opacities.size(0);
    uint32_t num_scales = scales.size(1);
    if (N) {
        at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();
        compute_relocation_kernel<<<(N + N_THREADS - 1) / N_THREADS, N_THREADS, 0, stream>>>(
            (int)N, (int)num_scales,
            opacities.data_ptr<float>(),
            scales.data_ptr<float>(),
            ratios.data_ptr<int>(),
            binoms.data_ptr<float>(),
            n_max,
            new_opacities.data_ptr<float>(),
            new_scales.data_ptr<float>()
        );
    }
    return std::make_tuple(new_opacities, new_scales);
}

std::tuple<
    torch::Tensor,  // new_position_offsets
    torch::Tensor,  // new_opacities
    torch::Tensor  // new_scales
> compute_relocation_split_tensor(
    torch::Tensor &positions,
    torch::Tensor &quats,
    torch::Tensor &opacities,
    torch::Tensor &scales
) {
    DEVICE_GUARD(opacities);
    CHECK_INPUT(opacities);
    CHECK_INPUT(scales);
    CHECK_INPUT(positions);
    CHECK_INPUT(quats);

    torch::Tensor new_position_offsets = torch::empty_like(positions);
    torch::Tensor new_opacities = torch::empty_like(opacities);
    torch::Tensor new_scales = torch::empty_like(scales);

    uint32_t N = opacities.size(0);
    uint32_t num_scales = scales.size(1);
    if (N) {
        at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();
        compute_relocation_split_kernel<<<(N + N_THREADS - 1) / N_THREADS, N_THREADS, 0, stream>>>(
            (int)N,
            (float3*)positions.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            opacities.data_ptr<float>(),
            (float2*)scales.data_ptr<float>(),
            (float3*)new_position_offsets.data_ptr<float>(),
            new_opacities.data_ptr<float>(),
            (float2*)new_scales.data_ptr<float>()
        );
    }
    return std::make_tuple(
        new_position_offsets,
        new_opacities, new_scales
    );
}





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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    CHECK_INPUT(anisotropies);
    CHECK_INPUT(background);

    dim3 tile_bounds_dim3;
    tile_bounds_dim3.x = std::get<0>(tile_bounds);
    tile_bounds_dim3.y = std::get<1>(tile_bounds);
    tile_bounds_dim3.z = std::get<2>(tile_bounds);

    dim3 block_dim3;
    block_dim3.x = std::get<0>(block);
    block_dim3.y = std::get<1>(block);
    block_dim3.z = std::get<2>(block);

    dim3 img_size_dim3;
    img_size_dim3.x = std::get<0>(img_size);
    img_size_dim3.y = std::get<1>(img_size);
    img_size_dim3.z = std::get<2>(img_size);

    const int channels = colors.size(1);
    const int img_width = img_size_dim3.x;
    const int img_height = img_size_dim3.y;

    float4 intrins = {fx, fy, cx, cy};

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width}, float32
    );

    rasterize_simple_forward_kernel<<<tile_bounds_dim3, block_dim3>>>(
        tile_bounds_dim3,
        img_size_dim3,
        intrins,
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float2 *)anisotropies.contiguous().data_ptr<float>(),
        *(float3 *)background.contiguous().data_ptr<float>(),
        // outputs
        final_idx.contiguous().data_ptr<int>(),
        (float3 *)out_img.contiguous().data_ptr<float>(),
        out_alpha.contiguous().data_ptr<float>()
    );

    return std::make_tuple(final_idx, out_img, out_alpha);
}

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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(v_output);
    CHECK_INPUT(v_output_alpha);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("positions must have dimensions (num_points, 3)");
    }

    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }

    const int num_points = positions.size(0);
    const dim3 tile_bounds = {
        (img_width + block_width - 1) / block_width,
        (img_height + block_width - 1) / block_width,
        1
    };
    const dim3 block(block_width, block_width, 1);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);
    float4 intrins = {fx, fy, cx, cy};

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);
    torch::Tensor v_anisotropies = torch::zeros({num_points, 2}, options);

    rasterize_simple_backward_kernel<<<tile_bounds, block>>>(
        tile_bounds,
        img_size,
        intrins,
        gaussians_ids_sorted.contiguous().data_ptr<int>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float2 *)anisotropies.contiguous().data_ptr<float>(),
        *(float3 *)background.contiguous().data_ptr<float>(),
        final_idx.contiguous().data_ptr<int>(),
        output_alpha.contiguous().data_ptr<float>(),
        (float3 *)v_output.contiguous().data_ptr<float>(),
        v_output_alpha.contiguous().data_ptr<float>(),
        // outputs
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float3 *)v_colors.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        (float2 *)v_anisotropies.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_opacities, v_anisotropies
    );
}



std::tuple<
    torch::Tensor,  // final_index
    torch::Tensor,  // out_img
    torch::Tensor  // out_visibility
> rasterize_depth_forward_tensor(
    const int depth_mode,
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);
    CHECK_INPUT(anisotropies);

    dim3 tile_bounds_dim3;
    tile_bounds_dim3.x = std::get<0>(tile_bounds);
    tile_bounds_dim3.y = std::get<1>(tile_bounds);
    tile_bounds_dim3.z = std::get<2>(tile_bounds);

    dim3 block_dim3;
    block_dim3.x = std::get<0>(block);
    block_dim3.y = std::get<1>(block);
    block_dim3.z = std::get<2>(block);

    dim3 img_size_dim3;
    img_size_dim3.x = std::get<0>(img_size);
    img_size_dim3.y = std::get<1>(img_size);
    img_size_dim3.z = std::get<2>(img_size);

    const int img_width = img_size_dim3.x;
    const int img_height = img_size_dim3.y;

    float4 intrins = {fx, fy, cx, cy};

    auto dtype = positions.dtype();
    if (dtype != torch::kFloat32 && dtype != torch::kFloat16)
        AT_ERROR("Must be FP32 or FP16");
    if (axes_u.dtype() != positions.dtype() ||
        axes_v.dtype() != positions.dtype() ||
        opacities.dtype() != positions.dtype() ||
        anisotropies.dtype() != positions.dtype()
    ) AT_ERROR("dtype must match");

    auto int32 = positions.options().dtype(torch::kInt32);
    auto floatt = positions.options();
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 1}, floatt
    );
    torch::Tensor out_visibility = torch::zeros(
        {img_height, img_width, 2}, floatt
    );

    if (positions.dtype() == torch::kFloat32)
        rasterize_depth_forward_kernel<float><<<tile_bounds_dim3, block_dim3>>>(
            depth_mode,
            tile_bounds_dim3,
            img_size_dim3,
            intrins,
            gaussian_ids_sorted.contiguous().data_ptr<int32_t>(),
            (int2 *)tile_bins.contiguous().data_ptr<int>(),
            (vec3<float> *)positions.contiguous().data_ptr<float>(),
            (vec3<float> *)axes_u.contiguous().data_ptr<float>(),
            (vec3<float> *)axes_v.contiguous().data_ptr<float>(),
            opacities.contiguous().data_ptr<float>(),
            (vec2<float> *)anisotropies.contiguous().data_ptr<float>(),
            // outputs
            final_idx.contiguous().data_ptr<int>(),
            out_depth.contiguous().data_ptr<float>(),
            (vec2<float> *)out_visibility.contiguous().data_ptr<float>()
        );

    if (positions.dtype() == torch::kFloat16)
        rasterize_depth_forward_kernel<halfc><<<tile_bounds_dim3, block_dim3>>>(
            depth_mode,
            tile_bounds_dim3,
            img_size_dim3,
            intrins,
            gaussian_ids_sorted.contiguous().data_ptr<int32_t>(),
            (int2 *)tile_bins.contiguous().data_ptr<int>(),
            (vec3<halfc> *)positions.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)axes_u.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)axes_v.contiguous().data_ptr<at::Half>(),
            (halfc *)opacities.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)anisotropies.contiguous().data_ptr<at::Half>(),
            // outputs
            final_idx.contiguous().data_ptr<int>(),
            (halfc *)out_depth.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)out_visibility.contiguous().data_ptr<at::Half>()
        );

    return std::make_tuple(final_idx, out_depth, out_visibility);
}

std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_opacities
    torch::Tensor  // v_anisotropies
> rasterize_depth_backward_tensor(
    const int depth_mode,
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(output_depth);
    CHECK_INPUT(output_visibility);
    CHECK_INPUT(v_output_depth);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("positions must have dimensions (num_points, 3)");
    }
    if (output_visibility.ndimension() != 3 || output_visibility.size(2) != 2) {
        AT_ERROR("output_visibility must have dimensions (h, w, 2)");
    }

    const int num_points = positions.size(0);
    const dim3 tile_bounds = {
        (img_width + block_width - 1) / block_width,
        (img_height + block_width - 1) / block_width,
        1
    };
    const dim3 block(block_width, block_width, 1);
    const dim3 img_size = {img_width, img_height, 1};
    float4 intrins = {fx, fy, cx, cy};

    // rewritten to test if FP16 is faster (answer: no)

    auto dtype = positions.dtype();
    if (dtype != torch::kFloat32 && dtype != torch::kFloat16)
        AT_ERROR("Must be FP32 or FP16");
    if (axes_u.dtype() != positions.dtype() ||
        axes_v.dtype() != positions.dtype() ||
        opacities.dtype() != positions.dtype() ||
        anisotropies.dtype() != positions.dtype() ||
        output_depth.dtype() != positions.dtype() ||
        output_visibility.dtype() != positions.dtype() ||
        v_output_depth.dtype() != positions.dtype()
    ) AT_ERROR("dtype must match");

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);
    torch::Tensor v_anisotropies = torch::zeros({num_points, 2}, options);

    if (positions.dtype() == torch::kFloat32)
        rasterize_depth_backward_kernel<float><<<tile_bounds, block>>>(
            depth_mode,
            tile_bounds,
            img_size,
            intrins,
            gaussians_ids_sorted.contiguous().data_ptr<int>(),
            (int2 *)tile_bins.contiguous().data_ptr<int>(),
            (vec3<float> *)positions.contiguous().data_ptr<float>(),
            (vec3<float> *)axes_u.contiguous().data_ptr<float>(),
            (vec3<float> *)axes_v.contiguous().data_ptr<float>(),
            opacities.contiguous().data_ptr<float>(),
            (vec2<float> *)anisotropies.contiguous().data_ptr<float>(),
            final_idx.contiguous().data_ptr<int>(),
            output_depth.contiguous().data_ptr<float>(),
            (vec2<float> *)output_visibility.contiguous().data_ptr<float>(),
            v_output_depth.contiguous().data_ptr<float>(),
            // outputs
            (vec3<float> *)v_positions.contiguous().data_ptr<float>(),
            (vec2<float> *)v_positions_xy_abs.contiguous().data_ptr<float>(),
            (vec3<float> *)v_axes_u.contiguous().data_ptr<float>(),
            (vec3<float> *)v_axes_v.contiguous().data_ptr<float>(),
            v_opacities.contiguous().data_ptr<float>(),
            (vec2<float> *)v_anisotropies.contiguous().data_ptr<float>()
        );

    if (positions.dtype() == torch::kFloat16)
        rasterize_depth_backward_kernel<halfc><<<tile_bounds, block>>>(
            depth_mode,
            tile_bounds,
            img_size,
            intrins,
            gaussians_ids_sorted.contiguous().data_ptr<int>(),
            (int2 *)tile_bins.contiguous().data_ptr<int>(),
            (vec3<halfc> *)positions.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)axes_u.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)axes_v.contiguous().data_ptr<at::Half>(),
            (halfc *)opacities.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)anisotropies.contiguous().data_ptr<at::Half>(),
            final_idx.contiguous().data_ptr<int>(),
            (halfc *)output_depth.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)output_visibility.contiguous().data_ptr<at::Half>(),
            (halfc *)v_output_depth.contiguous().data_ptr<at::Half>(),
            // outputs
            (vec3<halfc> *)v_positions.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)v_positions_xy_abs.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)v_axes_u.contiguous().data_ptr<at::Half>(),
            (vec3<halfc> *)v_axes_v.contiguous().data_ptr<at::Half>(),
            (halfc *)v_opacities.contiguous().data_ptr<at::Half>(),
            (vec2<halfc> *)v_anisotropies.contiguous().data_ptr<at::Half>()
        );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_opacities, v_anisotropies
    );
}



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
    const torch::Tensor &anisotropies,
    // const torch::Tensor &background,
    const torch::Tensor &depth_grads,
    const torch::Tensor &depth_ref_im
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(ch_coeffs);
    CHECK_INPUT(opacities);
    CHECK_INPUT(anisotropies);
    // CHECK_INPUT(background);
    CHECK_INPUT(depth_grads);
    CHECK_INPUT(depth_ref_im);

    dim3 tile_bounds_dim3;
    tile_bounds_dim3.x = std::get<0>(tile_bounds);
    tile_bounds_dim3.y = std::get<1>(tile_bounds);
    tile_bounds_dim3.z = std::get<2>(tile_bounds);

    dim3 block_dim3;
    block_dim3.x = std::get<0>(block);
    block_dim3.y = std::get<1>(block);
    block_dim3.z = std::get<2>(block);

    dim3 img_size_dim3;
    img_size_dim3.x = std::get<0>(img_size);
    img_size_dim3.y = std::get<1>(img_size);
    img_size_dim3.z = std::get<2>(img_size);

    const int channels = colors.size(1);
    const int img_width = img_size_dim3.x;
    const int img_height = img_size_dim3.y;

    float4 intrins = {fx, fy, cx, cy};

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width}, float32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_depth_grad = torch::zeros(
        {img_height, img_width, 4}, float32
    );
    torch::Tensor out_reg_depth = torch::zeros(
        {img_height, img_width}, float32
    );
    torch::Tensor out_reg_normal = torch::zeros(
        {img_height, img_width}, float32
    );

    rasterize_forward_kernel<<<tile_bounds_dim3, block_dim3>>>(
        tile_bounds_dim3,
        img_size_dim3,
        intrins,
        depth_reg_pairwise_factor,
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        (unsigned)ch_degree_r, (unsigned)ch_degree_r_to_use,
        (unsigned)ch_degree_phi, (unsigned)ch_degree_phi_to_use,
        (float3 *)ch_coeffs.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float2 *)anisotropies.contiguous().data_ptr<float>(),
        // *(float3 *)background.contiguous().data_ptr<float>(),
        (float2 *)depth_grads.contiguous().data_ptr<float>(),
        (float3 *)depth_ref_im.contiguous().data_ptr<float>(),
        // outputs
        final_idx.contiguous().data_ptr<int>(),
        out_alpha.contiguous().data_ptr<float>(),
        (float3 *)out_img.contiguous().data_ptr<float>(),
        (float4 *)out_depth_grad.contiguous().data_ptr<float>(),
        out_reg_depth.contiguous().data_ptr<float>(),
        out_reg_normal.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        final_idx, out_alpha,
        out_img, out_depth_grad,
        out_reg_depth, out_reg_normal
    );
}


std::tuple<
    torch::Tensor, // v_positions
    torch::Tensor, // v_positions_xy_abs
    torch::Tensor, // v_axes_u
    torch::Tensor, // v_axes_v
    torch::Tensor, // v_colors
    torch::Tensor, // v_ch_coeffs
    // torch::Tensor, // v_ch_coeffs_abs
    torch::Tensor, // v_opacities
    torch::Tensor, // v_anisotropies
    // torch::Tensor, // v_background
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
    const float depth_reg_pairwise_factor,
    const torch::Tensor &gaussians_ids_sorted,
    const torch::Tensor &tile_bins,
    const torch::Tensor &positions,
    const torch::Tensor &axes_u,
    const torch::Tensor &axes_v,
    const torch::Tensor &colors,
    const torch::Tensor &ch_coeffs,
    const torch::Tensor &opacities,
    const torch::Tensor &anisotropies,
    // const torch::Tensor &background,
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(ch_coeffs);
    CHECK_INPUT(opacities);
    CHECK_INPUT(anisotropies);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(output_depth_grad);
    CHECK_INPUT(v_output_alpha);
    CHECK_INPUT(v_output);
    CHECK_INPUT(v_output_depth_grad);
    CHECK_INPUT(v_output_reg_depth);
    CHECK_INPUT(v_output_reg_normal);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("xys must have dimensions (num_points, 2)");
    }
    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }
    if (output_depth_grad.ndimension() != 3 || output_depth_grad.size(2) != 4) {
        AT_ERROR("output_depth_grad must have 3 dimensions");
    }
    if (v_output_depth_grad.ndimension() != 3 || v_output_depth_grad.size(2) != 4) {
        AT_ERROR("v_output_depth_grad must have 3 dimensions");
    }

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    const int num_points = positions.size(0);
    const dim3 tile_bounds = {
        (img_width + block_width - 1) / block_width,
        (img_height + block_width - 1) / block_width,
        1
    };
    const dim3 block(block_width, block_width, 1);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);
    float4 intrins = {fx, fy, cx, cy};

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_ch_coeffs = torch::zeros({num_points, dim_ch, channels}, options);
    // torch::Tensor v_ch_coeffs_abs = torch::zeros({num_points, 1}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);
    torch::Tensor v_anisotropies = torch::zeros({num_points, 2}, options);
    // torch::Tensor v_background = torch::zeros({3}, options);
    torch::Tensor v_depth_grad = torch::zeros({num_points, 2}, options);
    torch::Tensor v_depth_ref_im = torch::zeros({img_height, img_width, 3}, options);

    rasterize_backward_kernel<<<tile_bounds, block>>>(
        tile_bounds, img_size, intrins,
        (unsigned)ch_degree_r, (unsigned)ch_degree_r_to_use,
        (unsigned)ch_degree_phi, (unsigned)ch_degree_phi_to_use,
        depth_reg_pairwise_factor,
        gaussians_ids_sorted.contiguous().data_ptr<int>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        (float3 *)ch_coeffs.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        (float2 *)anisotropies.contiguous().data_ptr<float>(),
        // *(float3 *)background.contiguous().data_ptr<float>(),
        (float2 *)depth_grads.contiguous().data_ptr<float>(),
        (float3 *)depth_ref_im.contiguous().data_ptr<float>(),
        final_idx.contiguous().data_ptr<int>(),
        output_alpha.contiguous().data_ptr<float>(),
        (float4 *)output_depth_grad.contiguous().data_ptr<float>(),
        v_output_alpha.contiguous().data_ptr<float>(),
        (float3 *)v_output.contiguous().data_ptr<float>(),
        (float4 *)v_output_depth_grad.contiguous().data_ptr<float>(),
        v_output_reg_depth.contiguous().data_ptr<float>(),
        v_output_reg_normal.contiguous().data_ptr<float>(),
        // outputs
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float3 *)v_colors.contiguous().data_ptr<float>(),
        (float3 *)v_ch_coeffs.contiguous().data_ptr<float>(),
        // v_ch_coeffs_abs.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        (float2 *)v_anisotropies.contiguous().data_ptr<float>(),
        // (float3 *)v_background.contiguous().data_ptr<float>(),
        (float2 *)v_depth_grad.contiguous().data_ptr<float>(),
        (float3 *)v_depth_ref_im.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_ch_coeffs, //v_ch_coeffs_abs,
        v_opacities, v_anisotropies,
        // v_background,
        v_depth_grad, v_depth_ref_im
    );
}



torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    const unsigned block_width,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(rotation);

    if (rotation.numel() != 9) {
        AT_ERROR("rotation must be 3x3");
    }
    if (sh_coeffs.ndimension() != 2 ||
        sh_coeffs.size(0) != sh_degree*sh_degree ||
        sh_coeffs.size(1) != 3) {
        AT_ERROR("sh_coeffs must be (sh_regree**2, 3)");
    }

    const dim3 tile_bounds = {
        (w + block_width - 1) / block_width,
        (h + block_width - 1) / block_width,
        1
    };
    const dim3 block(block_width, block_width, 1);
    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor out_color = torch::empty({h, w, 3}, options);

    render_background_sh_forward_kernel<<<tile_bounds, block>>>(
        tile_bounds, img_size,
        fx, fy, cx, cy,
        rotation.contiguous().data_ptr<float>(),
        sh_degree,
        (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
        (float3 *)out_color.contiguous().data_ptr<float>()
    );

    return out_color;
}


std::tuple<
    torch::Tensor,  // v_rotation
    torch::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const unsigned block_width,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs,
    const torch::Tensor &out_color,
    const torch::Tensor &v_out_color
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(rotation);
    CHECK_INPUT(v_out_color);

    if (rotation.numel() != 9) {
        AT_ERROR("rotation must be 3x3");
    }
    if (sh_coeffs.ndimension() != 2 ||
        sh_coeffs.size(0) != sh_degree*sh_degree ||
        sh_coeffs.size(1) != 3) {
        AT_ERROR("sh_coeffs shape must be (sh_regree**2, 3)");
    }
    if (out_color.ndimension() != 3 ||
        out_color.size(0) != h ||
        out_color.size(1) != w ||
        out_color.size(2) != 3) {
        AT_ERROR("out_color shape must be (h, w, 3)");
    }
    if (v_out_color.ndimension() != 3 ||
        v_out_color.size(0) != h ||
        v_out_color.size(1) != w ||
        v_out_color.size(2) != 3) {
        AT_ERROR("v_out_color shape must be (h, w, 3)");
    }

    const dim3 tile_bounds = {
        (w + block_width - 1) / block_width,
        (h + block_width - 1) / block_width,
        1
    };
    const dim3 block(block_width, block_width, 1);
    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor v_rotation = torch::zeros({3, 3}, options);
    torch::Tensor v_sh_coeffs = torch::zeros({sh_degree*sh_degree, 3}, options);

    render_background_sh_backward_kernel<<<tile_bounds, block>>>(
        tile_bounds, img_size,
        fx, fy, cx, cy,
        rotation.contiguous().data_ptr<float>(),
        sh_degree,
        (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
        (float3 *)out_color.contiguous().data_ptr<float>(),
        (float3 *)v_out_color.contiguous().data_ptr<float>(),
        v_rotation.contiguous().data_ptr<float>(),
        (float3 *)v_sh_coeffs.contiguous().data_ptr<float>()
    );

    return std::make_tuple(v_rotation, v_sh_coeffs);
}
