#include "PixelWise.cuh"

#define TensorView _Slang_TensorView
#include "generated/slang_all.cu"
#undef TensorView

#include "common.cuh"

__global__ void blend_background_forward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    TensorView<float, 3> out_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= in_rgb.shape[0]*in_rgb.shape[1])
        return;
    unsigned y = gid / in_rgb.shape[1];
    unsigned x = gid % in_rgb.shape[1];

    float3 rgb = in_rgb.load3f(y, x);
    float alpha = in_alpha.load1f(y, x);
    float3 background = in_background.load3f(y, x);

    rgb = blend_background(rgb, alpha, background);

    out_rgb.store3f(y, x, rgb);
}


__global__ void blend_background_backward_kernel(
    const TensorView<float, 3> in_rgb,
    const TensorView<float, 3> in_alpha,
    const TensorView<float, 3> in_background,
    const TensorView<float, 3> v_out_rgb,
    TensorView<float, 3> v_in_rgb,
    TensorView<float, 3> v_in_alpha,
    TensorView<float, 3> v_in_background
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= in_rgb.shape[0]*in_rgb.shape[1])
        return;
    unsigned y = gid / in_rgb.shape[1];
    unsigned x = gid % in_rgb.shape[1];

    float3 rgb = in_rgb.load3f(y, x);
    float alpha = in_alpha.load1f(y, x);
    float3 background = in_background.load3f(y, x);

    float3 v_out = v_out_rgb.load3f(y, x);

    float3 v_rgb; float v_alpha; float3 v_background;
    blend_background_bwd(
        rgb, alpha, background,
        v_out,
        &v_rgb, &v_alpha, &v_background
    );

    v_in_rgb.store3f(y, x, v_rgb);
    v_in_alpha.store1f(y, x, v_alpha);
    v_in_background.store3f(y, x, v_background);

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

