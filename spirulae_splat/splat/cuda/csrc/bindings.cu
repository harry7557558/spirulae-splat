#include "bindings.h"
#include "helpers.cuh"
#include "projection.cuh"
#include "rasterization.cuh"
#include "rasterization_sorted.cuh"
#include "sh.cuh"

#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include <torch/extension.h>
#include <torch/types.h>
#include <tuple>
#include <optional>


#define BLOCK_DIM3 dim3(BLOCK_WIDTH, BLOCK_WIDTH, 1)


inline __host__ float4 tuple2float4(std::tuple<float, float, float, float> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)};
}

inline __host__ dim3 tuple2dim3(std::tuple<unsigned, unsigned, unsigned> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v)};
}

inline __host__ dim3 whb2tb(unsigned width, unsigned height) {
    return {
        (width + BLOCK_WIDTH - 1) / BLOCK_WIDTH,
        (height + BLOCK_WIDTH - 1) / BLOCK_WIDTH,
        1
    };
}


torch::Tensor compute_sh_forward_tensor(
    const std::string &method,
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs0,   // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
) {
    DEVICE_GUARD(viewdirs);
    unsigned num_bases = num_sh_bases(degree);
    if (coeffs0.ndimension() != 2 || coeffs0.size(0) != num_points ||
        coeffs0.size(1) != 3) {
        AT_ERROR("coeffs0 must have dimensions (N, 3)");
    }
    if (coeffs.ndimension() != 3 || coeffs.size(0) != num_points ||
        coeffs.size(1) != num_bases-1 || coeffs.size(2) != 3) {
        AT_ERROR("coeffs must have dimensions (N, D, 3)");
    }
    torch::Tensor colors = torch::empty({num_points, 3}, coeffs.options());

    #define _TEMP_ARGS  \
        num_points, degree, degrees_to_use, \
        (float3 *)viewdirs.contiguous().data_ptr<float>(), \
        coeffs0.contiguous().data_ptr<float>(), \
        coeffs.contiguous().data_ptr<float>(), \
        colors.contiguous().data_ptr<float>()

    if (method == "poly") {
        compute_sh_forward_kernel<SHType::Poly>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    } else if (method == "fast") {
        compute_sh_forward_kernel<SHType::Fast>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    } else {
        AT_ERROR("Invalid method: ", method);
    }

    #undef _TEMP_ARGS

    return colors;
}


std::tuple<torch::Tensor, torch::Tensor>
compute_sh_backward_tensor(
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
    torch::Tensor v_coeffs0 = torch::zeros({num_points, 3}, v_colors.options());
    torch::Tensor v_coeffs = torch::zeros({num_points, num_bases-1, 3}, v_colors.options());

    #define _TEMP_ARGS  \
        num_points, degree, degrees_to_use, \
        (float3 *)viewdirs.contiguous().data_ptr<float>(), \
        v_colors.contiguous().data_ptr<float>(), \
        v_coeffs0.contiguous().data_ptr<float>(), \
        v_coeffs.contiguous().data_ptr<float>()

    if (method == "poly") {
        compute_sh_backward_kernel<SHType::Poly>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    } else if (method == "fast") {
        compute_sh_backward_kernel<SHType::Fast>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    } else {
        AT_ERROR("Invalid method: ", method);
    }

    #undef _TEMP_ARGS

    return std::make_tuple(v_coeffs0, v_coeffs);
}


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
) {
    DEVICE_GUARD(means3d);

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    auto int32 = means3d.options().dtype(torch::kInt32);
    auto float32 = means3d.options().dtype(torch::kFloat32);
    torch::Tensor bounds_d = torch::zeros({num_points, 4}, int32);
    torch::Tensor num_tiles_hit_d = torch::zeros({num_points}, int32);
    torch::Tensor positions_d = torch::zeros({num_points, 3}, float32);
    torch::Tensor axes_u_d = torch::zeros({num_points, 3}, float32);
    torch::Tensor axes_v_d = torch::zeros({num_points, 3}, float32);
    // torch::Tensor depth_grads_d = torch::zeros({num_points, 2}, float32);

    #define _TEMP_ARGS \
        num_points, \
        (float3 *)means3d.contiguous().data_ptr<float>(), \
        (float2 *)scales.contiguous().data_ptr<float>(), \
        (float4 *)quats.contiguous().data_ptr<float>(), \
        viewmat.contiguous().data_ptr<float>(), \
        tuple2float4(intrins), tuple2float4(dist_coeffs), \
        tile_bounds, \
        clip_thresh, \
        /* Outputs */ \
        (int4 *)bounds_d.contiguous().data_ptr<int32_t>(), \
        num_tiles_hit_d.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions_d.contiguous().data_ptr<float>(), \
        (float3 *)axes_u_d.contiguous().data_ptr<float>(), \
        (float3 *)axes_v_d.contiguous().data_ptr<float>()

    if (camera_model == "")
        project_gaussians_forward_kernel<CameraType::Undistorted>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    else if (camera_model == "OPENCV_FISHEYE")
        project_gaussians_forward_kernel<CameraType::OPENCV_FISHEYE>
        <<<_LAUNGH_ARGS_1D(num_points)>>>(_TEMP_ARGS);
    else AT_ERROR("Invalid camera model: ", camera_model);

    #undef _TEMP_ARGS

    return std::make_tuple(
        bounds_d, num_tiles_hit_d,
        positions_d, axes_u_d, axes_v_d
        // depth_grads_d
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
    std::tuple<float, float, float, float> intrins,
    torch::Tensor &num_tiles_hit,  // [N], int
    torch::Tensor &v_positions,  // [N, 3]
    torch::Tensor &v_axes_u,  // [N, 3]
    torch::Tensor &v_axes_v  // [N, 3]
    // torch::Tensor &v_depth_grads  // [N, 2]
) {
    DEVICE_GUARD(means3d);

    auto int32 = means3d.options().dtype(torch::kInt32);
    auto float32 = means3d.options().dtype(torch::kFloat32);
    torch::Tensor v_means3d = torch::zeros({num_points, 3}, float32);
    torch::Tensor v_scales = torch::zeros({num_points, 2}, float32);
    torch::Tensor v_quats = torch::zeros({num_points, 4}, float32);
    torch::Tensor v_viewmat = torch::zeros({4, 4}, float32);

    project_gaussians_backward_kernel<<<_LAUNGH_ARGS_1D(num_points)>>>(
        num_points,
        (float3 *)means3d.contiguous().data_ptr<float>(),
        (float2 *)scales.contiguous().data_ptr<float>(),
        (float4 *)quats.contiguous().data_ptr<float>(),
        viewmat.contiguous().data_ptr<float>(),
        tuple2float4(intrins),
        num_tiles_hit.contiguous().data_ptr<int32_t>(),
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        // (float2 *)v_depth_grads.contiguous().data_ptr<float>(),
        // Outputs.
        (float3 *)v_means3d.contiguous().data_ptr<float>(),
        (float2 *)v_scales.contiguous().data_ptr<float>(),
        (float4 *)v_quats.contiguous().data_ptr<float>(),
        (float *)v_viewmat.contiguous().data_ptr<float>()
    );

    return std::make_tuple(v_means3d, v_scales, v_quats, v_viewmat);
}


torch::Tensor render_undistortion_map_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::tuple<float, float, float, float> dist_coeffs
) {
    const dim3 tile_bounds = whb2tb(w, h);
    const dim3 img_size = {w, h, 1};

    auto options = torch::dtype(torch::kFloat32).device(torch::kCUDA, -1);
    torch::Tensor out_img = torch::empty({h, w, 2}, options);

    if (camera_model == "OPENCV_FISHEYE")
        render_undistortion_map_kernel<CameraType::OPENCV_FISHEYE>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), tuple2float4(dist_coeffs),
            (float2 *)out_img.contiguous().data_ptr<float>()
        );

    else AT_ERROR("Invalid camera model: ", camera_model);

    return out_img;
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
    const unsigned img_height,
    const unsigned img_width
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(bounds);
    CHECK_INPUT(cum_tiles_hit);

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto int64 = positions.options().dtype(torch::kInt64);
    torch::Tensor gaussian_ids_unsorted =
        torch::zeros({num_intersects}, int32);
    torch::Tensor isect_ids_unsorted =
        torch::zeros({num_intersects}, int64);

    map_gaussian_to_intersects<<<_LAUNGH_ARGS_1D(num_points)>>>(
        num_points,
        (float3 *)positions.contiguous().data_ptr<float>(),
        (int4 *)bounds.contiguous().data_ptr<int32_t>(),
        cum_tiles_hit.contiguous().data_ptr<int32_t>(),
        tile_bounds,
        // Outputs.
        isect_ids_unsorted.contiguous().data_ptr<int64_t>(),
        gaussian_ids_unsorted.contiguous().data_ptr<int32_t>()
    );

    return std::make_tuple(isect_ids_unsorted, gaussian_ids_unsorted);
}


torch::Tensor get_tile_bin_edges_tensor(
    int num_intersects,
    const torch::Tensor &isect_ids_sorted, 
    const unsigned img_height,
    const unsigned img_width
) {
    DEVICE_GUARD(isect_ids_sorted);
    CHECK_INPUT(isect_ids_sorted);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    int num_tiles = tile_bounds.x * tile_bounds.y;
    torch::Tensor tile_bins = torch::zeros(
        {num_tiles, 2}, isect_ids_sorted.options().dtype(torch::kInt32)
    );
    get_tile_bin_edges<<<_LAUNGH_ARGS_1D(num_intersects)>>>(
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    // CHECK_INPUT(background);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    const int channels = colors.size(1);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    #define _TEMP_ARGS \
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        *(float3 *)background.contiguous().data_ptr<float>(), \
        /* outputs */ \
        final_idx.contiguous().data_ptr<int>(), \
        (float3 *)out_img.contiguous().data_ptr<float>(), \
        out_alpha.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simple_forward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }

    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        rasterize_simple_forward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(final_idx, out_img, out_alpha);
}

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
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    rasterize_simple_backward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        tile_bounds,
        img_size,
        tuple2float4(intrins),
        gaussians_ids_sorted.contiguous().data_ptr<int>(),
        (int2 *)tile_bins.contiguous().data_ptr<int>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
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
        v_opacities.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_opacities
    );
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_visibility = torch::zeros(
        {img_height, img_width, 2}, float32
    );

    #define _TEMP_ARGS \
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        /* outputs */ \
        final_idx.contiguous().data_ptr<int>(), \
        out_depth.contiguous().data_ptr<float>(), \
        (float2 *)out_visibility.contiguous().data_ptr<float>()

    if (camera_model == "") {
        if (depth_mode == "mean")
            rasterize_depth_forward_kernel<DepthMode::Mean, CameraType::Undistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                tile_bounds, img_size,
                tuple2float4(intrins), nullptr,
                _TEMP_ARGS
            );
        else if (depth_mode == "median")
            rasterize_depth_forward_kernel<DepthMode::Median, CameraType::Undistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                tile_bounds, img_size,
                tuple2float4(intrins), nullptr,
                _TEMP_ARGS
            );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        if (depth_mode == "mean")
            rasterize_depth_forward_kernel<DepthMode::Mean, CameraType::GenericDistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                tile_bounds, img_size,
                tuple2float4(intrins),
                (float2 *)undistortion_map.contiguous().data_ptr<float>(),
                _TEMP_ARGS
            );
        else if (depth_mode == "median")
            rasterize_depth_forward_kernel<DepthMode::Median, CameraType::GenericDistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                tile_bounds, img_size,
                tuple2float4(intrins),
                (float2 *)undistortion_map.contiguous().data_ptr<float>(),
                _TEMP_ARGS
            );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(final_idx, out_depth, out_visibility);
}

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
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    #define _TEMP_ARGS \
        tile_bounds, img_size, \
        tuple2float4(intrins), \
        gaussians_ids_sorted.contiguous().data_ptr<int>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        final_idx.contiguous().data_ptr<int>(), \
        output_depth.contiguous().data_ptr<float>(), \
        (float2 *)output_visibility.contiguous().data_ptr<float>(), \
        v_output_depth.contiguous().data_ptr<float>(), \
        /* outputs */ \
        (float3 *)v_positions.contiguous().data_ptr<float>(), \
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_u.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_v.contiguous().data_ptr<float>(), \
        v_opacities.contiguous().data_ptr<float>()

    if (depth_mode == "mean")
        rasterize_depth_backward_kernel<DepthMode::Mean><<<tile_bounds, BLOCK_DIM3>>>(
            _TEMP_ARGS
        );
    if (depth_mode == "median")
        rasterize_depth_backward_kernel<DepthMode::Median><<<tile_bounds, BLOCK_DIM3>>>(
            _TEMP_ARGS
        );

    #undef _TEMP_ARGS

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_opacities
    );
}



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
    // CHECK_INPUT(background);
    CHECK_INPUT(depth_ref_im);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    const int channels = colors.size(1);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 2}, float32
    );
    torch::Tensor out_normal = torch::zeros(
        {img_height, img_width, 3}, float32
    );
    torch::Tensor out_reg_depth = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    rasterize_forward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        tile_bounds, img_size,
        tuple2float4(intrins),
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
        // *(float3 *)background.contiguous().data_ptr<float>(),
        depth_ref_im.contiguous().data_ptr<float>(),
        // outputs
        final_idx.contiguous().data_ptr<int>(),
        out_alpha.contiguous().data_ptr<float>(),
        (float3 *)out_img.contiguous().data_ptr<float>(),
        (float2 *)out_depth.contiguous().data_ptr<float>(),
        (float3 *)out_normal.contiguous().data_ptr<float>(),
        out_reg_depth.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        final_idx, out_alpha,
        out_img, out_depth,
        out_normal, out_reg_depth
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(ch_coeffs);
    CHECK_INPUT(opacities);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(output_depth);
    CHECK_INPUT(v_output_alpha);
    CHECK_INPUT(v_output);
    CHECK_INPUT(v_output_depth);
    CHECK_INPUT(v_output_normal);
    CHECK_INPUT(v_output_reg_depth);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("xys must have dimensions (num_points, 2)");
    }
    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }
    if (output_depth.ndimension() != 3 || output_depth.size(2) != 2) {
        AT_ERROR("output_depth_grad must have 3 dimensions");
    }
    if (v_output_depth.ndimension() != 3 || v_output_depth.size(2) != 2) {
        AT_ERROR("v_output_depth must have 3 dimensions");
    }
    if (depth_ref_im.ndimension() != 3 || depth_ref_im.size(2) != 1) {
        AT_ERROR("depth_ref_im must have 3 dimensions");
    }

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    const int num_points = positions.size(0);
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_ch_coeffs = torch::zeros({num_points, dim_ch, channels}, options);
    // torch::Tensor v_ch_coeffs_abs = torch::zeros({num_points, 1}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);
    // torch::Tensor v_background = torch::zeros({3}, options);
    torch::Tensor v_depth_ref_im = torch::zeros({img_height, img_width, 1}, options);

    rasterize_backward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        tile_bounds, img_size,
        tuple2float4(intrins),
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
        // *(float3 *)background.contiguous().data_ptr<float>(),
        depth_ref_im.contiguous().data_ptr<float>(),
        final_idx.contiguous().data_ptr<int>(),
        output_alpha.contiguous().data_ptr<float>(),
        (float2 *)output_depth.contiguous().data_ptr<float>(),
        v_output_alpha.contiguous().data_ptr<float>(),
        (float3 *)v_output.contiguous().data_ptr<float>(),
        (float2 *)v_output_depth.contiguous().data_ptr<float>(),
        (float3 *)v_output_normal.contiguous().data_ptr<float>(),
        v_output_reg_depth.contiguous().data_ptr<float>(),
        // outputs
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float3 *)v_colors.contiguous().data_ptr<float>(),
        (float3 *)v_ch_coeffs.contiguous().data_ptr<float>(),
        // v_ch_coeffs_abs.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        // (float3 *)v_background.contiguous().data_ptr<float>(),
        v_depth_ref_im.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_ch_coeffs, //v_ch_coeffs_abs,
        v_opacities,
        // v_background,
        v_depth_ref_im
    );
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    const int channels = colors.size(1);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width}, int32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 2}, float32
    );
    torch::Tensor out_normal = torch::zeros(
        {img_height, img_width, 3}, float32
    );
    torch::Tensor out_depth_reg = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    #define _TEMP_ARGS \
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        /* outputs */ \
        final_idx.contiguous().data_ptr<int>(), \
        out_alpha.contiguous().data_ptr<float>(), \
        (float3 *)out_img.contiguous().data_ptr<float>(), \
        (float2 *)out_depth.contiguous().data_ptr<float>(), \
        (float3 *)out_normal.contiguous().data_ptr<float>(), \
        out_depth_reg.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simplified_forward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        rasterize_simplified_forward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(
        final_idx, out_alpha,
        out_img,
        out_depth, out_normal, out_depth_reg
    );
}


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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(output_depth);
    CHECK_INPUT(v_output_alpha);
    CHECK_INPUT(v_output_img);
    CHECK_INPUT(v_output_depth);
    CHECK_INPUT(v_output_normal);
    CHECK_INPUT(v_output_depth_reg);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("xys must have dimensions (num_points, 2)");
    }
    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }
    if (v_output_depth.ndimension() != 3 || v_output_depth.size(2) != 2) {
        AT_ERROR("v_output_depth must have 3 dimensions");
    }
    if (v_output_normal.ndimension() != 3 || v_output_normal.size(2) != 3) {
        AT_ERROR("v_output_normal must have 3 dimensions");
    }
    if (v_output_depth_reg.ndimension() != 3 || v_output_depth_reg.size(2) != 1) {
        AT_ERROR("v_output_depth_reg must have 3 dimensions");
    }

    const int num_points = positions.size(0);
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    #define _TEMP_ARGS \
        gaussians_ids_sorted.contiguous().data_ptr<int>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        final_idx.contiguous().data_ptr<int>(), \
        output_alpha.contiguous().data_ptr<float>(), \
        (float2 *)output_depth.contiguous().data_ptr<float>(), \
        v_output_alpha.contiguous().data_ptr<float>(), \
        (float3 *)v_output_img.contiguous().data_ptr<float>(), \
        (float2 *)v_output_depth.contiguous().data_ptr<float>(), \
        (float3 *)v_output_normal.contiguous().data_ptr<float>(), \
        v_output_depth_reg.contiguous().data_ptr<float>(), \
        /* outputs */ \
        (float3 *)v_positions.contiguous().data_ptr<float>(), \
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_u.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_v.contiguous().data_ptr<float>(), \
        (float3 *)v_colors.contiguous().data_ptr<float>(), \
        v_opacities.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simplified_backward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        rasterize_simplified_backward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_opacities
    );
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(gaussian_ids_sorted);
    CHECK_INPUT(tile_bins);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor num_intersects = torch::empty(
        {img_height, img_width}, int32
    );
    torch::Tensor sorted_indices = torch::empty(
        {img_height, img_width, MAX_SORTED_SPLATS}, int32
    );
    torch::Tensor sorted_depths = torch::empty(
        {img_height, img_width, MAX_SORTED_SPLATS}, float32
    );

    #define _TEMP_ARGS \
        gaussian_ids_sorted.contiguous().data_ptr<int32_t>(), \
        (int2 *)tile_bins.contiguous().data_ptr<int>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        /* outputs */ \
        num_intersects.contiguous().data_ptr<int>(), \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        sorted_depths.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_indices_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        rasterize_indices_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(
        num_intersects, sorted_indices, sorted_depths
    );
}


void sort_per_pixel_tensor(
    const std::string &method,
    const unsigned img_height,
    const unsigned img_width,
    torch::Tensor &num_intersects,  // [h, w]
    torch::Tensor &indices,  // [h, w, MAX_SORTED_SPLATS]
    torch::Tensor &depths  // [h, w, MAX_SORTED_SPLATS]
) {
    DEVICE_GUARD(num_intersects);
    CHECK_INPUT(num_intersects);
    CHECK_INPUT(indices);
    CHECK_INPUT(depths);

    const dim3 tile_bounds = {
        (img_width + BLOCK_WIDTH_PPS - 1) / BLOCK_WIDTH_PPS,
        (img_height + BLOCK_WIDTH_PPS - 1) / BLOCK_WIDTH_PPS,
        1
    };
    const dim3 block = {BLOCK_WIDTH_PPS, BLOCK_WIDTH_PPS, 1};

    #define _TEMP_ARGS \
        <<<tile_bounds, block>>>( \
            img_height, img_width, \
            num_intersects.contiguous().data_ptr<int>(), \
            indices.contiguous().data_ptr<int32_t>(), \
            depths.contiguous().data_ptr<float>() \
        )

    if (method == "insertion") {
        sort_per_pixel_kernel<PerPixelSortType::InsertionSort> _TEMP_ARGS;
    } else if (method == "quick") {
        sort_per_pixel_kernel<PerPixelSortType::QuickSort> _TEMP_ARGS;
    } else if (method == "heap") {
        sort_per_pixel_kernel<PerPixelSortType::HeapSort> _TEMP_ARGS;
    } else if (method == "random_quick") {
        sort_per_pixel_kernel<PerPixelSortType::RandomizedQuickSort> _TEMP_ARGS;
    } else if (method == "network") {
        sort_per_pixel_kernel<PerPixelSortType::NetworkSort> _TEMP_ARGS;
    } else {
        AT_ERROR("Invalid sorting method: ", method);
    }

    #undef _TEMP_ARGS
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    // CHECK_INPUT(background);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, 3}, float32
    );
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    #define _TEMP_ARGS \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        *(float3 *)background.contiguous().data_ptr<float>(), \
        /* outputs */ \
        (float3 *)out_img.contiguous().data_ptr<float>(), \
        out_alpha.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simple_sorted_forward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        rasterize_simple_sorted_forward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(out_img, out_alpha);
}


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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(num_intersects);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
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
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    rasterize_simple_sorted_backward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        img_size,
        tuple2float4(intrins),
        num_intersects.contiguous().data_ptr<int>(),
        sorted_indices.contiguous().data_ptr<int32_t>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        *(float3 *)background.contiguous().data_ptr<float>(),
        output_alpha.contiguous().data_ptr<float>(),
        (float3 *)v_output.contiguous().data_ptr<float>(),
        v_output_alpha.contiguous().data_ptr<float>(),
        // outputs
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float3 *)v_colors.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_opacities
    );
}


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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor final_idx = torch::zeros(
        {img_height, img_width, 2}, int32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_visibility = torch::zeros(
        {img_height, img_width, 2}, float32
    );

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    #define _TEMP_ARGS \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        /* outputs */ \
        final_idx.contiguous().data_ptr<int>(), \
        out_depth.contiguous().data_ptr<float>(), \
        (float2 *)out_visibility.contiguous().data_ptr<float>()

    if (camera_model == "") {
        if (depth_mode == "mean")
            rasterize_depth_sorted_forward_kernel<DepthMode::Mean, CameraType::Undistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                img_height, img_width,
                tuple2float4(intrins), nullptr,
                _TEMP_ARGS
            );
        else if (depth_mode == "median")
            rasterize_depth_sorted_forward_kernel<DepthMode::Median, CameraType::Undistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                img_height, img_width,
                tuple2float4(intrins), nullptr,
                _TEMP_ARGS
            );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        if (depth_mode == "mean")
            rasterize_depth_sorted_forward_kernel<DepthMode::Mean, CameraType::GenericDistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                img_height, img_width,
                tuple2float4(intrins),
                (float2 *)undistortion_map.contiguous().data_ptr<float>(),
                _TEMP_ARGS
            );
        else if (depth_mode == "median")
            rasterize_depth_sorted_forward_kernel<DepthMode::Median, CameraType::GenericDistorted>
            <<<tile_bounds, BLOCK_DIM3>>>(
                img_height, img_width,
                tuple2float4(intrins),
                (float2 *)undistortion_map.contiguous().data_ptr<float>(),
                _TEMP_ARGS
            );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(final_idx, out_depth, out_visibility);
}

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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(final_idx);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(opacities);
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
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    #define _TEMP_ARGS \
        img_size, \
        tuple2float4(intrins), \
        final_idx.contiguous().data_ptr<int>(), \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        output_depth.contiguous().data_ptr<float>(), \
        (float2 *)output_visibility.contiguous().data_ptr<float>(), \
        v_output_depth.contiguous().data_ptr<float>(), \
        /* outputs */ \
        (float3 *)v_positions.contiguous().data_ptr<float>(), \
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_u.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_v.contiguous().data_ptr<float>(), \
        v_opacities.contiguous().data_ptr<float>()

    if (depth_mode == "mean")
        rasterize_depth_sorted_backward_kernel<DepthMode::Mean><<<tile_bounds, BLOCK_DIM3>>>(
            _TEMP_ARGS
        );
    if (depth_mode == "median")
        rasterize_depth_sorted_backward_kernel<DepthMode::Median><<<tile_bounds, BLOCK_DIM3>>>(
            _TEMP_ARGS
        );

    #undef _TEMP_ARGS

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_opacities
    );
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(ch_coeffs);
    CHECK_INPUT(opacities);
    // CHECK_INPUT(background);
    CHECK_INPUT(depth_ref_im);

    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};

    const int channels = colors.size(1);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 2}, float32
    );
    torch::Tensor out_normal = torch::zeros(
        {img_height, img_width, 3}, float32
    );
    torch::Tensor out_reg_depth = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    rasterize_sorted_forward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        img_size,
        tuple2float4(intrins),
        depth_reg_pairwise_factor,
        sorted_indices.contiguous().data_ptr<int32_t>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        (unsigned)ch_degree_r, (unsigned)ch_degree_r_to_use,
        (unsigned)ch_degree_phi, (unsigned)ch_degree_phi_to_use,
        (float3 *)ch_coeffs.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        // *(float3 *)background.contiguous().data_ptr<float>(),
        depth_ref_im.contiguous().data_ptr<float>(),
        // outputs
        out_alpha.contiguous().data_ptr<float>(),
        (float3 *)out_img.contiguous().data_ptr<float>(),
        (float2 *)out_depth.contiguous().data_ptr<float>(),
        (float3 *)out_normal.contiguous().data_ptr<float>(),
        out_reg_depth.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        out_alpha,
        out_img, out_depth,
        out_normal, out_reg_depth
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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(num_intersects);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(ch_coeffs);
    CHECK_INPUT(opacities);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(output_depth);
    CHECK_INPUT(v_output_alpha);
    CHECK_INPUT(v_output);
    CHECK_INPUT(v_output_depth);
    CHECK_INPUT(v_output_normal);
    CHECK_INPUT(v_output_reg_depth);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("xys must have dimensions (num_points, 2)");
    }
    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }
    if (output_depth.ndimension() != 3 || output_depth.size(2) != 2) {
        AT_ERROR("output_depth_grad must have 3 dimensions");
    }
    if (v_output_depth.ndimension() != 3 || v_output_depth.size(2) != 2) {
        AT_ERROR("v_output_depth must have 3 dimensions");
    }
    if (depth_ref_im.ndimension() != 3 || depth_ref_im.size(2) != 1) {
        AT_ERROR("depth_ref_im must have 3 dimensions");
    }

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    const int num_points = positions.size(0);
    const dim3 tile_bounds = whb2tb(img_width, img_height);
    const dim3 img_size = {img_width, img_height, 1};
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_ch_coeffs = torch::zeros({num_points, dim_ch, channels}, options);
    // torch::Tensor v_ch_coeffs_abs = torch::zeros({num_points, 1}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);
    // torch::Tensor v_background = torch::zeros({3}, options);
    torch::Tensor v_depth_ref_im = torch::zeros({img_height, img_width, 1}, options);

    rasterize_sorted_backward_kernel<<<tile_bounds, BLOCK_DIM3>>>(
        img_size,
        tuple2float4(intrins),
        (unsigned)ch_degree_r, (unsigned)ch_degree_r_to_use,
        (unsigned)ch_degree_phi, (unsigned)ch_degree_phi_to_use,
        depth_reg_pairwise_factor,
        num_intersects.contiguous().data_ptr<int>(),
        sorted_indices.contiguous().data_ptr<int32_t>(),
        (float3 *)positions.contiguous().data_ptr<float>(),
        (float3 *)axes_u.contiguous().data_ptr<float>(),
        (float3 *)axes_v.contiguous().data_ptr<float>(),
        (float3 *)colors.contiguous().data_ptr<float>(),
        (float3 *)ch_coeffs.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        // *(float3 *)background.contiguous().data_ptr<float>(),
        depth_ref_im.contiguous().data_ptr<float>(),
        output_alpha.contiguous().data_ptr<float>(),
        (float2 *)output_depth.contiguous().data_ptr<float>(),
        v_output_alpha.contiguous().data_ptr<float>(),
        (float3 *)v_output.contiguous().data_ptr<float>(),
        (float2 *)v_output_depth.contiguous().data_ptr<float>(),
        (float3 *)v_output_normal.contiguous().data_ptr<float>(),
        v_output_reg_depth.contiguous().data_ptr<float>(),
        // outputs
        (float3 *)v_positions.contiguous().data_ptr<float>(),
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(),
        (float3 *)v_axes_u.contiguous().data_ptr<float>(),
        (float3 *)v_axes_v.contiguous().data_ptr<float>(),
        (float3 *)v_colors.contiguous().data_ptr<float>(),
        (float3 *)v_ch_coeffs.contiguous().data_ptr<float>(),
        // v_ch_coeffs_abs.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        // (float3 *)v_background.contiguous().data_ptr<float>(),
        v_depth_ref_im.contiguous().data_ptr<float>()
    );

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_ch_coeffs, //v_ch_coeffs_abs,
        v_opacities,
        // v_background,
        v_depth_ref_im
    );
}



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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);

    const int channels = colors.size(1);

    auto int32 = positions.options().dtype(torch::kInt32);
    auto float32 = positions.options().dtype(torch::kFloat32);
    torch::Tensor out_alpha = torch::zeros(
        {img_height, img_width, 1}, float32
    );
    torch::Tensor out_img = torch::zeros(
        {img_height, img_width, channels}, float32
    );
    torch::Tensor out_depth = torch::zeros(
        {img_height, img_width, 2}, float32
    );
    torch::Tensor out_normal = torch::zeros(
        {img_height, img_width, 3}, float32
    );
    torch::Tensor out_depth_reg = torch::zeros(
        {img_height, img_width, 1}, float32
    );

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    #define _TEMP_ARGS \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        /* outputs */ \
        out_alpha.contiguous().data_ptr<float>(), \
        (float3 *)out_img.contiguous().data_ptr<float>(), \
        (float2 *)out_depth.contiguous().data_ptr<float>(), \
        (float3 *)out_normal.contiguous().data_ptr<float>(), \
        out_depth_reg.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simplified_sorted_forward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        rasterize_simplified_sorted_forward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(
        out_alpha, out_img,
        out_depth, out_normal, out_depth_reg
    );
}


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
) {
    DEVICE_GUARD(positions);
    CHECK_INPUT(num_intersects);
    CHECK_INPUT(sorted_indices);
    CHECK_INPUT(positions);
    CHECK_INPUT(axes_u);
    CHECK_INPUT(axes_v);
    CHECK_INPUT(colors);
    CHECK_INPUT(opacities);
    CHECK_INPUT(output_alpha);
    CHECK_INPUT(output_depth);
    CHECK_INPUT(v_output_alpha);
    CHECK_INPUT(v_output_img);
    CHECK_INPUT(v_output_depth);
    CHECK_INPUT(v_output_normal);
    CHECK_INPUT(v_output_depth_reg);

    if (positions.ndimension() != 2 || positions.size(1) != 3) {
        AT_ERROR("xys must have dimensions (num_points, 2)");
    }
    if (colors.ndimension() != 2 || colors.size(1) != 3) {
        AT_ERROR("colors must have 2 dimensions");
    }
    if (v_output_depth.ndimension() != 3 || v_output_depth.size(2) != 2) {
        AT_ERROR("v_output_depth must have 3 dimensions");
    }
    if (v_output_normal.ndimension() != 3 || v_output_normal.size(2) != 3) {
        AT_ERROR("v_output_normal must have 3 dimensions");
    }
    if (v_output_depth_reg.ndimension() != 3 || v_output_depth_reg.size(2) != 1) {
        AT_ERROR("v_output_depth_reg must have 3 dimensions");
    }

    const int num_points = positions.size(0);
    const int channels = colors.size(1);

    auto options = positions.options();
    torch::Tensor v_positions = torch::zeros({num_points, 3}, options);
    torch::Tensor v_positions_xy_abs = torch::zeros({num_points, 2}, options);
    torch::Tensor v_axes_u = torch::zeros({num_points, 3}, options);
    torch::Tensor v_axes_v = torch::zeros({num_points, 3}, options);
    torch::Tensor v_colors = torch::zeros({num_points, channels}, options);
    torch::Tensor v_opacities = torch::zeros({num_points, 1}, options);

    const dim3 tile_bounds = whb2tb(img_width, img_height);

    #define _TEMP_ARGS \
        num_intersects.contiguous().data_ptr<int>(), \
        sorted_indices.contiguous().data_ptr<int32_t>(), \
        (float3 *)positions.contiguous().data_ptr<float>(), \
        (float3 *)axes_u.contiguous().data_ptr<float>(), \
        (float3 *)axes_v.contiguous().data_ptr<float>(), \
        (float3 *)colors.contiguous().data_ptr<float>(), \
        opacities.contiguous().data_ptr<float>(), \
        output_alpha.contiguous().data_ptr<float>(), \
        (float2 *)output_depth.contiguous().data_ptr<float>(), \
        v_output_alpha.contiguous().data_ptr<float>(), \
        (float3 *)v_output_img.contiguous().data_ptr<float>(), \
        (float2 *)v_output_depth.contiguous().data_ptr<float>(), \
        (float3 *)v_output_normal.contiguous().data_ptr<float>(), \
        v_output_depth_reg.contiguous().data_ptr<float>(), \
        /* outputs */ \
        (float3 *)v_positions.contiguous().data_ptr<float>(), \
        (float2 *)v_positions_xy_abs.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_u.contiguous().data_ptr<float>(), \
        (float3 *)v_axes_v.contiguous().data_ptr<float>(), \
        (float3 *)v_colors.contiguous().data_ptr<float>(), \
        v_opacities.contiguous().data_ptr<float>()

    if (camera_model == "") {
        rasterize_simplified_sorted_backward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        rasterize_simplified_sorted_backward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            img_height, img_width,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(
        v_positions, v_positions_xy_abs,
        v_axes_u, v_axes_v,
        v_colors, v_opacities
    );
}





torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
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

    const dim3 tile_bounds = whb2tb(w, h);
    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor out_color = torch::empty({h, w, 3}, options);

    if (camera_model == "") {
        render_background_sh_forward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            rotation.contiguous().data_ptr<float>(),
            sh_degree,
            (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
            (float3 *)out_color.contiguous().data_ptr<float>()
        );
    }

    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        render_background_sh_forward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            rotation.contiguous().data_ptr<float>(),
            sh_degree,
            (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
            (float3 *)out_color.contiguous().data_ptr<float>()
        );
    }

    return out_color;
}


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

    const dim3 tile_bounds = whb2tb(w, h);
    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor v_rotation = torch::zeros({3, 3}, options);
    torch::Tensor v_sh_coeffs = torch::zeros({sh_degree*sh_degree, 3}, options);

    #define _TEMP_ARGS \
        rotation.contiguous().data_ptr<float>(), \
        sh_degree, \
        (float3 *)sh_coeffs.contiguous().data_ptr<float>(), \
        (float3 *)out_color.contiguous().data_ptr<float>(), \
        (float3 *)v_out_color.contiguous().data_ptr<float>(), \
        v_rotation.contiguous().data_ptr<float>(), \
        (float3 *)v_sh_coeffs.contiguous().data_ptr<float>()

    if (camera_model == "") {
        render_background_sh_backward_kernel<CameraType::Undistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        render_background_sh_backward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, BLOCK_DIM3>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(v_rotation, v_sh_coeffs);
}
