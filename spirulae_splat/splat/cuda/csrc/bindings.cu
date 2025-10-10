#include "bindings.h"
#include "helpers.cuh"
#include "projection.cuh"
#include "rasterization.cuh"
#include "sh.cuh"
#include "misc.cuh"

#include "splat_tile_intersector.cuh"

#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include <tuple>
#include <optional>

#define TORCH_INDUCTOR_CPP_WRAPPER
#include <torch/extension.h>
#include <torch/types.h>


#define BLOCK_DIM3 dim3(BLOCK_WIDTH, BLOCK_WIDTH, 1)


inline __host__ float4 tuple2float4(std::tuple<float, float, float, float> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)};
}

inline __host__ dim3 tuple2dim3(std::tuple<unsigned, unsigned, unsigned> v) {
    return {std::get<0>(v), std::get<1>(v), std::get<2>(v)};
}

inline __host__ dim3 whb2tb(unsigned width, unsigned height, unsigned block_width=BLOCK_WIDTH) {
    return {
        (width + block_width - 1) / block_width,
        (height + block_width - 1) / block_width,
        1
    };
}


template<typename T, int ndim>
TensorView<T, ndim> tensor2view(torch::Tensor& tensor) {
    TensorView<T, ndim> view;
    view.data = tensor.data_ptr<T>();
    for (int i = 0; i < ndim; i++) {
        view.shape[i] = tensor.size(i);
        view.strides[i] = *(tensor.strides().begin() + i);
    }
    return view;
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

    // unsigned block_width = BLOCK_WIDTH;
    unsigned block_width = 32;  // 1024 threads
    const dim3 tile_bounds = whb2tb(w, h, block_width);
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
        (float3 *)v_rotation.contiguous().data_ptr<float>(), \
        (float3 *)v_sh_coeffs.contiguous().data_ptr<float>()

    if (camera_model == "") {
        render_background_sh_backward_kernel<CameraType::Undistorted>
        <<<tile_bounds, dim3(block_width, block_width, 1)>>>(
            tile_bounds, img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        render_background_sh_backward_kernel<CameraType::GenericDistorted>
        <<<tile_bounds, dim3(block_width, block_width, 1)>>>(
            tile_bounds, img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(v_rotation, v_sh_coeffs);
}


torch::Tensor compute_per_splat_losses_forward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(quats);

    const size_t num_points = opacities.size(0);
    const bool is_3dgs = (scales.size(-1) == 3);

    if (scales.ndimension() != 2 || scales.size(0) != num_points ||
        (scales.size(1) != 3 && scales.size(1) != 2))
        AT_ERROR("scales shape must be (n, 2) or (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");

    torch::Tensor loss = torch::zeros({kNumPerSplatLosses}, opacities.options());

    per_splat_losses_forward_kernel<<<_LAUNGH_ARGS_1D(num_points)>>>(
        is_3dgs,
        num_points,
        scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        quats.contiguous().data_ptr<float>(),
        loss.data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );

    return loss;
}


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
compute_per_splat_losses_backward_tensor(
    torch::Tensor &scales,  // [N, 3] or [N, 2]
    torch::Tensor &opacities,  // [N, 1]
    torch::Tensor &quats,  // [N, 4]
    torch::Tensor &v_losses,  // [kNumPerSplatLosses]
    float mcmc_opacity_reg_weight,
    float mcmc_scale_reg_weight,
    float max_gauss_ratio,
    float scale_regularization_weight,
    float erank_reg_weight,
    float erank_reg_weight_s3,
    float quat_norm_reg_weight
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(quats);
    CHECK_INPUT(v_losses);

    const size_t num_points = opacities.size(0);
    const bool is_3dgs = (scales.size(-1) == 3);

    if (scales.ndimension() != 2 || scales.size(0) != num_points ||
        (scales.size(1) != 3 && scales.size(1) != 2))
        AT_ERROR("scales shape must be (n, 2) or (n, 3)");
    if (opacities.ndimension() != 2 || opacities.size(0) != num_points || opacities.size(1) != 1)
        AT_ERROR("opacities shape must be (n, 1)");
    if (quats.ndimension() != 2 || quats.size(0) != num_points || quats.size(1) != 4)
        AT_ERROR("quats shape must be (n, 4)");
    if (v_losses.ndimension() != 1 || v_losses.size(0) != kNumPerSplatLosses)
        AT_ERROR("v_losses shape must be (kNumPerSplatLosses,)");

    torch::Tensor v_scales = torch::empty_like(scales);
    torch::Tensor v_opacities = torch::empty_like(opacities);
    torch::Tensor v_quats = torch::empty_like(quats);

    per_splat_losses_backward_kernel<<<_LAUNGH_ARGS_1D(num_points)>>>(
        is_3dgs,
        num_points,
        scales.contiguous().data_ptr<float>(),
        opacities.contiguous().data_ptr<float>(),
        quats.contiguous().data_ptr<float>(),
        v_losses.data_ptr<float>(),
        v_scales.contiguous().data_ptr<float>(),
        v_opacities.contiguous().data_ptr<float>(),
        v_quats.contiguous().data_ptr<float>(),
        mcmc_opacity_reg_weight,
        mcmc_scale_reg_weight,
        max_gauss_ratio,
        scale_regularization_weight,
        erank_reg_weight,
        erank_reg_weight_s3,
        quat_norm_reg_weight
    );

    return std::make_tuple(v_scales, v_opacities, v_quats);
}



torch::Tensor blend_background_forward_tensor(
    torch::Tensor &rgb,  // [H, W, 3]
    torch::Tensor &alpha,  // [H, W, 1]
    torch::Tensor &background  // [H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(alpha);
    CHECK_CUDA(background);

    if (rgb.ndimension() != 3 || rgb.size(2) != 3)
        AT_ERROR("rgb shape must be (h, w, 3)");
    long h = rgb.size(0), w = rgb.size(1);
    if (alpha.ndimension() != 3 || alpha.size(0) != h || alpha.size(1) != w || alpha.size(2) != 1)
        AT_ERROR("alpha shape must be (h, w, 1)");
    if (background.ndimension() != 3 || background.size(0) != h || background.size(1) != w || background.size(2) != 3)
        AT_ERROR("background shape must be (h, w, 3)");

    torch::Tensor out_rgb = torch::empty({h, w, 3}, rgb.options());

    blend_background_forward_kernel<<<_LAUNGH_ARGS_1D(h*w)>>>(
        tensor2view<float, 3>(rgb), tensor2view<float, 3>(alpha), tensor2view<float, 3>(background),
        tensor2view<float, 3>(out_rgb)
    );

    return out_rgb;
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
blend_background_backward_tensor(
    torch::Tensor &rgb,  // [H, W, 3]
    torch::Tensor &alpha,  // [H, W, 1]
    torch::Tensor &background,  // [H, W, 3]
    torch::Tensor &v_out_rgb  // [H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(alpha);
    CHECK_CUDA(background);
    CHECK_CUDA(v_out_rgb);

    if (rgb.ndimension() != 3 || rgb.size(2) != 3)
        AT_ERROR("rgb shape must be (h, w, 3)");
    long h = rgb.size(0), w = rgb.size(1);
    if (alpha.ndimension() != 3 || alpha.size(0) != h || alpha.size(1) != w || alpha.size(2) != 1)
        AT_ERROR("alpha shape must be (h, w, 1)");
    if (background.ndimension() != 3 || background.size(0) != h || background.size(1) != w || background.size(2) != 3)
        AT_ERROR("background shape must be (h, w, 3)");
    if (v_out_rgb.ndimension() != 3 || v_out_rgb.size(0) != h || v_out_rgb.size(1) != w || v_out_rgb.size(2) != 3)
        AT_ERROR("v_out_rgb shape must be (h, w, 3)");

    torch::Tensor v_rgb = torch::empty({h, w, 3}, rgb.options());
    torch::Tensor v_alpha = torch::empty({h, w, 1}, alpha.options());
    torch::Tensor v_background = torch::empty({h, w, 3}, background.options());

    blend_background_backward_kernel<<<_LAUNGH_ARGS_1D(h*w)>>>(
        tensor2view<float, 3>(rgb), tensor2view<float, 3>(alpha), tensor2view<float, 3>(background),
        tensor2view<float, 3>(v_out_rgb),
        tensor2view<float, 3>(v_rgb), tensor2view<float, 3>(v_alpha), tensor2view<float, 3>(v_background)
    );

    return std::make_tuple(v_rgb, v_alpha, v_background);
}


std::tuple<torch::Tensor, torch::Tensor>
intersect_splat_tile(
    torch::Tensor& means,
    torch::Tensor& scales,
    torch::Tensor& opacs,
    torch::Tensor& quats,
    torch::Tensor& viewmats,
    torch::Tensor& Ks
) {
    SplatBuffers splat_buffers = {means, scales, opacs, quats};
    TileBuffers tile_buffers = {viewmats, Ks};

    return SplatTileIntersector(
        means.options(),
        splat_buffers,
        tile_buffers
    ).getIntersections_lbvh();
}
