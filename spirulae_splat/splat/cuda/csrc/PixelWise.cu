#include "PixelWise.cuh"

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangPixelWise {
#include "generated/set_namespace.cuh"
#include "generated/pixel_wise.cuh"
}
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
namespace SlangPPISP {
#include "generated/set_namespace.cuh"
#include "generated/ppisp.cuh"
}
#endif

#include "common.cuh"


CameraDistortionCoeffsBuffer::CameraDistortionCoeffsBuffer(
    const CameraDistortionCoeffsTensor &tensors
) {
    coeffs = nullptr;

    if (tensors.has_value()) {
        CHECK_INPUT(tensors.value());
        if (tensors.value().size(-1) != 10)
            AT_ERROR("dist_coeffs must be (..., 10)");
        coeffs = (float*)tensors.value().data_ptr<float>();
    }
}

// ================
// Type Conversion
// ================

__global__ void uint8_image_to_float_kernel(
    const TensorView<uint8_t, 4> img_in,
    TensorView<float, 4> img_out
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = img_in.shape[0], H = img_in.shape[1], W = img_in.shape[2], C = img_in.shape[3];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    for (int i = 0; i < C; ++i) {
        uint8_t c_in = img_in.at(bid, y, x, i);
        float c_out = (float)c_in / 255.0f;
        img_out.at(bid, y, x, i) = c_out;
    }
}

__global__ void uint16_image_to_float_kernel(
    const TensorView<uint16_t, 4> img_in,
    TensorView<float, 4> img_out
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = img_in.shape[0], H = img_in.shape[1], W = img_in.shape[2], C = img_in.shape[3];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    for (int i = 0; i < C; ++i) {
        uint16_t c_in = img_in.at(bid, y, x, i);
        float c_out = (float)c_in / 65535.0f;
        img_out.at(bid, y, x, i) = c_out;
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor uint8_image_to_float_tensor(
    at::Tensor &img_in  // [B, H, W, C]
) {
    DEVICE_GUARD(img_in);
    CHECK_CUDA(img_in);

    long b = img_in.size(0), h = img_in.size(1), w = img_in.size(2), c = img_in.size(3);

    at::Tensor img_out = at::empty({b, h, w, c}, img_in.options().dtype(at::kFloat));

    uint8_image_to_float_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<uint8_t, 4>(img_in),
        tensor2view<float, 4>(img_out)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return img_out;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor uint16_image_to_float_tensor(
    at::Tensor &img_in  // [B, H, W, C]
) {
    DEVICE_GUARD(img_in);
    CHECK_CUDA(img_in);

    long b = img_in.size(0), h = img_in.size(1), w = img_in.size(2), c = img_in.size(3);

    at::Tensor img_out = at::empty({b, h, w, c}, img_in.options().dtype(at::kFloat));

    uint16_image_to_float_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<uint16_t, 4>(img_in),
        tensor2view<float, 4>(img_out)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return img_out;
}


// ================
// Rendered Depth to Expected Depth
// ================

__global__ void rendered_depth_to_expected_depth_forward_kernel(
    const TensorView<float, 4> in_depth,
    const TensorView<float, 4> in_transmittance,
    TensorView<float, 4> out_depth,
    float* __restrict__ max_out_depth
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_depth.shape[0], H = in_depth.shape[1], W = in_depth.shape[2];
    if (bid >= B)
        return;
    bool inside = (gid < H*W);

    float depth = 0.0f;
    if (inside) {
        unsigned y = gid / W;
        unsigned x = gid % W;

        depth = in_depth.load1(bid, y, x);
        float transmittance = in_transmittance.load1(bid, y, x);

        depth = SlangPixelWise::rendered_depth_to_expected_depth(depth, transmittance);

        out_depth.store1(bid, y, x, depth);
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    warpMax(depth, warp);
    if (warp.thread_rank() == 0) {
        atomicMax(&max_out_depth[bid], depth);
    }
}

__global__ void rendered_depth_to_expected_depth_filter_kernel(
    const TensorView<float, 4> in_transmittance,
    TensorView<float, 4> depth,
    const float* __restrict__ max_out_depth
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = depth.shape[0], H = depth.shape[1], W = depth.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float transmittance = in_transmittance.load1(bid, y, x);
    if (transmittance == 1.0f) {
        depth.store1(bid, y, x, max_out_depth[bid]);
    }
    // note that in backward, transmittance=1 automatically leads to zero output gradient
}


__global__ void rendered_depth_to_expected_depth_backward_kernel(
    const TensorView<float, 4> in_depth,
    const TensorView<float, 4> in_transmittance,
    const TensorView<float, 4> v_out_depth,
    TensorView<float, 4> v_in_depth,
    TensorView<float, 4> v_in_transmittance
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_depth.shape[0], H = in_depth.shape[1], W = in_depth.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float depth = in_depth.load1(bid, y, x);
    float transmittance = in_transmittance.load1(bid, y, x);

    float v_out = v_out_depth.load1(bid, y, x);

    float v_depth, v_transmittance;
    SlangPixelWise::rendered_depth_to_expected_depth_bwd(
        depth, transmittance,
        v_out,
        &v_depth, &v_transmittance
    );

    v_in_depth.store1(bid, y, x, v_depth);
    v_in_transmittance.store1(bid, y, x, v_transmittance);
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor rendered_depth_to_expected_depth_forward_tensor(
    at::Tensor &depth,  // [B, H, W, 1]
    at::Tensor &transmittance  // [B, H, W, 1]
) {
    DEVICE_GUARD(depth);
    CHECK_CUDA(depth);
    CHECK_CUDA(transmittance);

    if (depth.ndimension() != 4 || depth.size(-1) != 1)
        AT_ERROR("depth shape must be (b, h, w, 1)");
    long b = depth.size(0), h = depth.size(1), w = depth.size(2);
    if (transmittance.ndimension() != 4 || transmittance.size(0) != b || transmittance.size(1) != h || transmittance.size(2) != w || transmittance.size(3) != 1)
        AT_ERROR("transmittance shape must be (b, h, w, 1)");

    at::Tensor out_depth = at::empty_like(depth);
    at::Tensor max_depth = at::zeros({b,}, depth.options());

    rendered_depth_to_expected_depth_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(depth), tensor2view<float, 4>(transmittance),
        tensor2view<float, 4>(out_depth), max_depth.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    rendered_depth_to_expected_depth_filter_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(transmittance),
        tensor2view<float, 4>(out_depth),
        max_depth.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_depth;
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
rendered_depth_to_expected_depth_backward_tensor(
    at::Tensor &depth,  // [B, H, W, 1]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &v_out_depth  // [B, H, W, 1]
) {
    DEVICE_GUARD(depth);
    CHECK_CUDA(depth);
    CHECK_CUDA(transmittance);
    CHECK_CUDA(v_out_depth);

    long b = depth.size(0), h = depth.size(1), w = depth.size(2);

    at::Tensor v_depth = at::empty_like(depth);
    at::Tensor v_transmittance = at::empty_like(transmittance);

    rendered_depth_to_expected_depth_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(depth), tensor2view<float, 4>(transmittance),
        tensor2view<float, 4>(v_out_depth),
        tensor2view<float, 4>(v_depth), tensor2view<float, 4>(v_transmittance)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_depth, v_transmittance);
}



// ================
// Blend Background
// ================

__global__ void blend_background_forward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_transmittance,
    const TensorView<float, 4> in_background,
    TensorView<float, 4> out_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);
    float transmittance = in_transmittance.load1(bid, y, x);
    float3 background = in_background.load3(bid, y, x);

    rgb = SlangPixelWise::blend_background(rgb, transmittance, background);

    out_rgb.store3(bid, y, x, rgb);
}

__global__ void blend_background_backward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_transmittance,
    const TensorView<float, 4> in_background,
    const TensorView<float, 4> v_out_rgb,
    TensorView<float, 4> v_in_rgb,
    TensorView<float, 4> v_in_transmittance,
    TensorView<float, 4> v_in_background
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);
    float transmittance = in_transmittance.load1(bid, y, x);
    float3 background = in_background.load3(bid, y, x);

    float3 v_out = v_out_rgb.load3(bid, y, x);

    float3 v_rgb; float v_transmittance; float3 v_background;
    SlangPixelWise::blend_background_bwd(
        rgb, transmittance, background,
        v_out,
        &v_rgb, &v_transmittance, &v_background
    );

    v_in_rgb.store3(bid, y, x, v_rgb);
    v_in_transmittance.store1(bid, y, x, v_transmittance);
    v_in_background.store3(bid, y, x, v_background);

}

/*[AutoHeaderGeneratorExport]*/
at::Tensor blend_background_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &background  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(transmittance);
    CHECK_CUDA(background);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);
    if (transmittance.ndimension() != 4 || transmittance.size(0) != b || transmittance.size(1) != h || transmittance.size(2) != w || transmittance.size(3) != 1)
        AT_ERROR("transmittance shape must be (b, h, w, 1)");
    if (background.ndimension() != 4 || background.size(0) != b || background.size(1) != h || background.size(2) != w || background.size(3) != 3)
        AT_ERROR("background shape must be (b, h, w, 3)");

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    blend_background_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(transmittance), tensor2view<float, 4>(background),
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor, at::Tensor>
blend_background_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    at::Tensor &background,  // [B, H, W, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(transmittance);
    CHECK_CUDA(background);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());
    at::Tensor v_transmittance = at::empty({b, h, w, 1}, transmittance.options());
    at::Tensor v_background = at::empty({b, h, w, 3}, background.options());

    blend_background_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(transmittance), tensor2view<float, 4>(background),
        tensor2view<float, 4>(v_out_rgb),
        tensor2view<float, 4>(v_rgb), tensor2view<float, 4>(v_transmittance), tensor2view<float, 4>(v_background)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_rgb, v_transmittance, v_background);
}


// ================
// Blend Background with Random Noise
// ================

__forceinline__ __device__ uint32_t hash_uint3(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t hash = a;

    hash *= 0x01000193;
    hash = (hash << 16) | (hash >> 16);

    hash ^= b;
    hash *= 0x01000193;
    hash = (hash << 16) | (hash >> 16);

    hash ^= c;
    hash *= 0x01000193;

    hash ^= hash >> 15;
    hash *= 0x85ebca6b;
    hash ^= hash >> 13;
    hash *= 0xc2b2ae35;
    hash ^= hash >> 16;

    return hash;
}

template<bool is_linear>
__global__ void blend_background_noise_forward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_transmittance,
    const float randomize_weight,
    const uint32_t seed,
    TensorView<float, 4> out_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);
    float transmittance = in_transmittance.load1(bid, y, x);

    float3 background;
    background.x = (float)hash_uint3(seed + 0, gid, bid) * exp2f(-32.0f);
    background.y = (float)hash_uint3(seed + 1, gid, bid) * exp2f(-32.0f);
    background.z = (float)hash_uint3(seed + 2, gid, bid) * exp2f(-32.0f);
    background = 0.5 + 0.5*randomize_weight * background;
    if (is_linear) {
        background.x = SlangPixelWise::srgb_to_linear_rgb(background.x);
        background.y = SlangPixelWise::srgb_to_linear_rgb(background.y);
        background.z = SlangPixelWise::srgb_to_linear_rgb(background.z);
    }

    rgb = SlangPixelWise::blend_background(rgb, transmittance, background);

    out_rgb.store3(bid, y, x, rgb);
}

template<bool is_linear>
__global__ void blend_background_noise_backward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_transmittance,
    const float randomize_weight,
    const uint32_t seed,
    const TensorView<float, 4> v_out_rgb,
    TensorView<float, 4> v_in_rgb,
    TensorView<float, 4> v_in_transmittance
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);
    float transmittance = in_transmittance.load1(bid, y, x);

    float3 background;
    background.x = (float)hash_uint3(seed + 0, gid, bid) * exp2f(-32.0f);
    background.y = (float)hash_uint3(seed + 1, gid, bid) * exp2f(-32.0f);
    background.z = (float)hash_uint3(seed + 2, gid, bid) * exp2f(-32.0f);
    background = 0.5 + 0.5*randomize_weight * background;
    if (is_linear) {
        background.x = SlangPixelWise::srgb_to_linear_rgb(background.x);
        background.y = SlangPixelWise::srgb_to_linear_rgb(background.y);
        background.z = SlangPixelWise::srgb_to_linear_rgb(background.z);
    }

    float3 v_out = v_out_rgb.load3(bid, y, x);

    float3 v_rgb; float v_transmittance; float3 v_background;
    SlangPixelWise::blend_background_bwd(
        rgb, transmittance, background,
        v_out,
        &v_rgb, &v_transmittance, &v_background
    );

    v_in_rgb.store3(bid, y, x, v_rgb);
    v_in_transmittance.store1(bid, y, x, v_transmittance);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor blend_background_noise_forward_tensor(
    bool is_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(transmittance);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);
    if (transmittance.ndimension() != 4 || transmittance.size(0) != b || transmittance.size(1) != h || transmittance.size(2) != w || transmittance.size(3) != 1)
        AT_ERROR("transmittance shape must be (b, h, w, 1)");

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    (is_linear ? blend_background_noise_forward_kernel<true> : blend_background_noise_forward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(transmittance),
        randomize_weight, seed,
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
blend_background_noise_backward_tensor(
    bool is_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &transmittance,  // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed,
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(transmittance);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());
    at::Tensor v_transmittance = at::empty({b, h, w, 1}, transmittance.options());

    (is_linear ? blend_background_noise_backward_kernel<true> : blend_background_noise_backward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(transmittance),
        randomize_weight, seed,
        tensor2view<float, 4>(v_out_rgb),
        tensor2view<float, 4>(v_rgb), tensor2view<float, 4>(v_transmittance)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_rgb, v_transmittance);
}


// ================
// Log Map Image
// ================

template<bool is_input_linear>
__global__ void rgb_to_srgb_forward_kernel(
    const TensorView<float, 4> in_rgb,
    const float* __restrict__ color_matrix_buffer,
    TensorView<float, 4> out_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);

    float3x3 color_matrix;
    color_matrix[0].x = color_matrix_buffer[0];
    color_matrix[0].y = color_matrix_buffer[1];
    color_matrix[0].z = color_matrix_buffer[2];
    color_matrix[1].x = color_matrix_buffer[3];
    color_matrix[1].y = color_matrix_buffer[4];
    color_matrix[1].z = color_matrix_buffer[5];
    color_matrix[2].x = color_matrix_buffer[6];
    color_matrix[2].y = color_matrix_buffer[7];
    color_matrix[2].z = color_matrix_buffer[8];

    if (is_input_linear)
        rgb = SlangPixelWise::linear_rgb_to_srgb(rgb, color_matrix);
    else
        rgb = SlangPixelWise::rgb_to_srgb(rgb, color_matrix);

    out_rgb.store3(bid, y, x, rgb);
}

template<bool is_input_linear>
__global__ void rgb_to_srgb_backward_kernel(
    const TensorView<float, 4> in_rgb,
    const float* __restrict__ color_matrix_buffer,
    const TensorView<float, 4> v_out_rgb,
    TensorView<float, 4> v_in_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_rgb.shape[0], H = in_rgb.shape[1], W = in_rgb.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 rgb = in_rgb.load3(bid, y, x);

    float3x3 color_matrix;
    color_matrix[0].x = color_matrix_buffer[0];
    color_matrix[0].y = color_matrix_buffer[1];
    color_matrix[0].z = color_matrix_buffer[2];
    color_matrix[1].x = color_matrix_buffer[3];
    color_matrix[1].y = color_matrix_buffer[4];
    color_matrix[1].z = color_matrix_buffer[5];
    color_matrix[2].x = color_matrix_buffer[6];
    color_matrix[2].y = color_matrix_buffer[7];
    color_matrix[2].z = color_matrix_buffer[8];

    float3 v_out = v_out_rgb.load3(bid, y, x);

    float3 v_rgb;
    if (is_input_linear)
        v_rgb = SlangPixelWise::linear_rgb_to_srgb_bwd(rgb, color_matrix, v_out);
    else
        v_rgb = SlangPixelWise::rgb_to_srgb_bwd(rgb, color_matrix, v_out);

    v_in_rgb.store3(bid, y, x, v_rgb);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor rgb_to_srgb_forward_tensor(
    bool is_input_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &color_matrix   // [3, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_INPUT(color_matrix);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    (is_input_linear ? rgb_to_srgb_forward_kernel<true> : rgb_to_srgb_forward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb),
        color_matrix.data_ptr<float>(),
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor rgb_to_srgb_backward_tensor(
    bool is_input_linear,
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &color_matrix,   // [3, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());

    (is_input_linear ? rgb_to_srgb_backward_kernel<true> : rgb_to_srgb_backward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb),
        color_matrix.data_ptr<float>(),
        tensor2view<float, 4>(v_out_rgb),
        tensor2view<float, 4>(v_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_rgb;
}


// ================
// Depth to Points
// ================


__global__ void depth_to_points_forward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> in_depths,  // [B, H, W, 1]
    TensorView<float, 4> out_points  // [B, H, W, 3]
) {
    const int B = in_depths.shape[0],
        H = in_depths.shape[1],
        W = in_depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float in_depth = in_depths.load1(bid, j, i);
    float3 out_point = SlangPixelWise::depth_to_point(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE,
        is_ray_depth,
        in_depth
    );
    out_points.store3(bid, j, i, out_point);
}

__global__ void depth_to_points_backward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> in_depths,  // [B, H, W, 1]
    const TensorView<float, 4> v_out_points,  // [B, H, W, 3]
    TensorView<float, 4> v_in_depths  // [B, H, W, 1]
) {
    const int B = v_out_points.shape[0],
        H = v_out_points.shape[1],
        W = v_out_points.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float in_depth = in_depths.load1(bid, j, i);
    float3 v_out_point = v_out_points.load3(bid, j, i);
    float v_in_depth = SlangPixelWise::depth_to_point_vjp(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE,
        is_ray_depth,
        in_depth, v_out_point
    );
    v_in_depths.store1(bid, j, i, v_in_depth);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_points_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != intrins.size(0))
        AT_ERROR("depths and intrins batch dimension mismatch");

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor out_points = at::empty({b, h, w, 3}, depths.options());

    depth_to_points_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs, is_ray_depth,
        tensor2view<float, 4>(depths), tensor2view<float, 4>(out_points)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_points;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_points_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor in_depths,  // [B, H, W, 1]
    at::Tensor v_out_points  // [B, H, W, 3]
) {
    DEVICE_GUARD(in_depths);
    CHECK_CUDA(in_depths);
    CHECK_CUDA(v_out_points);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");

    int b = in_depths.size(0), h = in_depths.size(1), w = in_depths.size(2);
    at::Tensor v_in_depths = at::empty_like(in_depths);

    depth_to_points_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs, is_ray_depth,
        tensor2view<float, 4>(in_depths), tensor2view<float, 4>(v_out_points),
        tensor2view<float, 4>(v_in_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_in_depths;
}



// ================
// Depth to Normal
// ================


__global__ void depth_to_normal_forward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    TensorView<float, 4> normals  // [B, H, W, 3]
) {
    const int B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    constexpr int TILE = 16;  // blockDim.x and blockDim.y; blockDim.z should be 1
    uint32_t bid = blockIdx.z;
    uint32_t i = blockIdx.x * TILE + threadIdx.x;
    uint32_t j = blockIdx.y * TILE + threadIdx.y;

    bool inside = (bid < B && i < W && j < H);

    // Zero for border pixels (consistent with PyTorch implementation)
    if (inside && (i == 0 || i == W-1 || j == 0 || j == H-1)) {
        normals.store3(bid, j, i, make_float3(0.0f));
        inside = false;
    }

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
#if 0
    if (!inside) return;
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 normal = depth_to_normal(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth,
        depth
    );
#else
    __shared__ float3 shared_points[TILE+2][TILE+2];
    #pragma unroll 2
    for (int k = threadIdx.y * blockDim.x + threadIdx.x;
            k < (TILE+2)*(TILE+2); k += TILE*TILE) {
        int it = k % (TILE+2), jt = k / (TILE+2);
        int ig = int(blockIdx.x * TILE) + it - 1;
        int jg = int(blockIdx.y * TILE) + jt - 1;
        float depth = (ig >= 0 && ig < W && jg >= 0 && jg < H) ?
            depths.load1(bid, jg, ig) : 0.0f;
        float3 ray = SlangPixelWise::generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, dist_coeffs,
            camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth
        );
        shared_points[jt][it] = ray * depth;
    }
    __syncthreads();
    if (!inside) return;

    FixedArray<float3, 4> points;
    int it = threadIdx.x+1, jt = threadIdx.y+1;
    points[0] = shared_points[jt][it-1];
    points[1] = shared_points[jt][it+1];
    points[2] = shared_points[jt-1][it];
    points[3] = shared_points[jt+1][it];
    float3 normal = SlangPixelWise::points_to_normal(points);
#endif
    normals.store3(bid, j, i, normal);

}


__global__ void depth_to_normal_backward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    const TensorView<float, 4> v_normals,  // [B, H, W, 3]
    TensorView<float, 4> v_depths  // [B, H, W, 1]
) {
    const int B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    constexpr int TILE = 16;  // blockDim.x and blockDim.y; blockDim.z should be 1
    uint32_t bid = blockIdx.z;
    uint32_t i = blockIdx.x * TILE + threadIdx.x;
    uint32_t j = blockIdx.y * TILE + threadIdx.y;

    bool inside = (bid < B && i < W && j < H);

    // Zero for border pixels (consistent with PyTorch implementation)
    if (i == 0 || i == W-1 || j == 0 || j == H-1) {
        inside = false;
    }

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
#if 0
    if (!inside) return;
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 v_normal = v_normals.load3(bid, j, i);
    float4 v_depth;
    depth_to_normal_vjp(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth,
        depth, v_normal, &v_depth
    );
    v_depths.atomicStore1(bid, j, i-1, v_depth.x);
    v_depths.atomicStore1(bid, j, i+1, v_depth.y);
    v_depths.atomicStore1(bid, j-1, i, v_depth.z);
    v_depths.atomicStore1(bid, j+1, i, v_depth.w);
#else
    __shared__ float4 shared_points[TILE+2][TILE+2];
    #pragma unroll 2
    for (int k = threadIdx.y * blockDim.x + threadIdx.x;
            k < (TILE+2)*(TILE+2); k += TILE*TILE) {
        int it = k % (TILE+2), jt = k / (TILE+2);
        int ig = int(blockIdx.x * TILE) + it - 1;
        int jg = int(blockIdx.y * TILE) + jt - 1;
        float depth = (ig >= 0 && ig < W && jg >= 0 && jg < H) ?
            depths.load1(bid, jg, ig) : 0.0f;
        float3 ray = SlangPixelWise::generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, dist_coeffs,
            camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth
        );
        shared_points[jt][it] = make_float4(ray.x, ray.y, ray.z, depth);
    }
    __syncthreads();
    if (!inside) return;

    float3 v_normal = v_normals.load3(bid, j, i);

    FixedArray<float3, 4> rays;
    FixedArray<float3, 4> points;
    int it = threadIdx.x+1, jt = threadIdx.y+1;
    float4 t;
    t = shared_points[jt][it-1]; rays[0] = {t.x, t.y, t.z}; points[0] = rays[0] * t.w;
    t = shared_points[jt][it+1]; rays[1] = {t.x, t.y, t.z}; points[1] = rays[1] * t.w;
    t = shared_points[jt-1][it]; rays[2] = {t.x, t.y, t.z}; points[2] = rays[2] * t.w;
    t = shared_points[jt+1][it]; rays[3] = {t.x, t.y, t.z}; points[3] = rays[3] * t.w;
    FixedArray<float3, 4> v_points;
    SlangPixelWise::points_to_normal_vjp(points, v_normal, &v_points);

    v_depths.atomicStore1(bid, j, i-1, dot(v_points[0], rays[0]));
    v_depths.atomicStore1(bid, j, i+1, dot(v_points[1], rays[1]));
    v_depths.atomicStore1(bid, j-1, i, dot(v_points[2], rays[2]));
    v_depths.atomicStore1(bid, j+1, i, dot(v_points[3], rays[3]));
#endif
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_normal_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != intrins.size(0))
        AT_ERROR("depths and intrins batch dimension mismatch");

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor normals = at::empty({b, h, w, 3}, depths.options());

    depth_to_normal_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths), tensor2view<float, 4>(normals)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return normals;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_normal_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor v_normals  // [B, H, W, 3]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(intrins);
    CHECK_CUDA(v_normals);

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor v_depths = zeros_like<float>(depths);

    depth_to_normal_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths),
        tensor2view<float, 4>(v_normals),
        tensor2view<float, 4>(v_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_depths;
}



// ================
// Loss between Depth and Normal
// ================


__global__ void depth_normal_loss_forward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    const TensorView<float, 4> gt_normals,  // [B, H, W, 3]
    TensorView<float, 4> losses  // [B, H, W, 1]
) {
    const int B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    constexpr int TILE = 16;  // blockDim.x and blockDim.y; blockDim.z should be 1
    uint32_t bid = blockIdx.z;
    uint32_t i = blockIdx.x * TILE + threadIdx.x;
    uint32_t j = blockIdx.y * TILE + threadIdx.y;

    bool inside = (bid < B && i < W && j < H);

    // Zero for border pixels (consistent with PyTorch implementation)
    if (inside && (i == 0 || i == W-1 || j == 0 || j == H-1)) {
        losses.store1(bid, j, i, 0.0f);
        inside = false;
    }

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
#if 1
    if (!inside) return;
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 gt_normal = gt_normals.load3(bid, j, i);
    float loss = SlangPixelWise::depth_normal_loss(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth,
        depth, gt_normal
    );
#else
    // TODO: pre-undistort and store in shared memory
#endif
    losses.store1(bid, j, i, loss);
}


__global__ void depth_normal_loss_backward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    const TensorView<float, 4> gt_normals,  // [B, H, W, 3]
    const TensorView<float, 4> v_losses,  // [B, H, W, 1]
    TensorView<float, 4> v_depths,  // [B, H, W, 1]
    TensorView<float, 4> v_gt_normals  // [B, H, W, 3]
) {
    const int B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    constexpr int TILE = 16;  // blockDim.x and blockDim.y; blockDim.z should be 1
    uint32_t bid = blockIdx.z;
    uint32_t i = blockIdx.x * TILE + threadIdx.x;
    uint32_t j = blockIdx.y * TILE + threadIdx.y;

    bool inside = (bid < B && i < W && j < H);

    // Zero for border pixels (consistent with PyTorch implementation)
    if (i == 0 || i == W-1 || j == 0 || j == H-1) {
        inside = false;
    }

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
#if 1
    if (!inside) {
        v_gt_normals.store3(bid, j, i, float3{0.0f, 0.0f, 0.0f});
        return;
    }
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 gt_normal = gt_normals.load3(bid, j, i);
    float v_loss = v_losses.load1(bid, j, i);
    float4 v_depth = {0.0f, 0.0f, 0.0f, 0.0f};
    float3 v_gt_normal = float3{0.0f, 0.0f, 0.0f};
    SlangPixelWise::depth_normal_loss_vjp(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE, is_ray_depth,
        depth, gt_normal, v_loss, &v_depth, &v_gt_normal
    );
    v_depths.atomicStore1(bid, j, i-1, v_depth.x);
    v_depths.atomicStore1(bid, j, i+1, v_depth.y);
    v_depths.atomicStore1(bid, j-1, i, v_depth.z);
    v_depths.atomicStore1(bid, j+1, i, v_depth.w);
    v_gt_normals.store3(bid, j, i, v_gt_normal);
#else
    // TODO: pre-undistort and store in shared memory
#endif
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_normal_loss_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor gt_normals  // [B, H, W, 3]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_CUDA(gt_normals);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != intrins.size(0))
        AT_ERROR("depths and intrins batch dimension mismatch");
    if (gt_normals.ndimension() != 4 || gt_normals.size(-1) != 3)
        AT_ERROR("gt_normals shape must be (B, H, W, 3)");
    if (gt_normals.size(0) != intrins.size(0))
        AT_ERROR("gt_normals and intrins batch dimension mismatch");

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor losses = at::empty_like(depths);

    depth_normal_loss_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths), tensor2view<float, 4>(gt_normals),
        tensor2view<float, 4>(losses)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return losses;
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor> depth_normal_loss_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    at::Tensor depths,  // [B, H, W, 1]
    at::Tensor gt_normals,  // [B, H, W, 3]
    at::Tensor v_losses  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_CUDA(gt_normals);
    CHECK_CUDA(v_losses);
    CHECK_INPUT(intrins);

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor v_depths = zeros_like<float>(depths);
    at::Tensor v_gt_normals = at::empty_like(gt_normals);

    depth_normal_loss_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths), tensor2view<float, 4>(gt_normals),
        tensor2view<float, 4>(v_losses),
        tensor2view<float, 4>(v_depths), tensor2view<float, 4>(v_gt_normals)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_depths, v_gt_normals);
}


// ================
// Ray Depth To Linear Depth
// ================


__global__ void ray_depth_to_linear_depth_forward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> in_depths,  // [B, H, W, 1]
    TensorView<float, 4> out_depths  // [B, H, W, 1]
) {
    const int B = in_depths.shape[0],
        H = in_depths.shape[1],
        W = in_depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float in_depth = in_depths.load1(bid, j, i);
    float out_depth = in_depth * SlangPixelWise::ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE
    );
    out_depths.store1(bid, j, i, out_depth);
}

__global__ void ray_depth_to_linear_depth_backward_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> v_out_depths,  // [B, H, W, 1]
    TensorView<float, 4> v_in_depths  // [B, H, W, 1]
) {
    const int B = v_out_depths.shape[0],
        H = v_out_depths.shape[1],
        W = v_out_depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float v_out_depth = v_out_depths.load1(bid, j, i);
    float factor = SlangPixelWise::ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        camera_model == ssplat::CameraModelType::FISHEYE
    );
    float v_in_depth = factor * v_out_depth;
    v_in_depths.store1(bid, j, i, v_in_depth);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ray_depth_to_linear_depth_forward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != intrins.size(0))
        AT_ERROR("depths and intrins batch dimension mismatch");

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor out_depths = at::empty_like(depths);

    ray_depth_to_linear_depth_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(depths), tensor2view<float, 4>(out_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_depths;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ray_depth_to_linear_depth_backward_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor v_out_depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(v_out_depths);
    CHECK_CUDA(v_out_depths);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (v_out_depths.ndimension() != 4 || v_out_depths.size(-1) != 1)
        AT_ERROR("v_out_depths shape must be (B, H, W, 1)");
    if (v_out_depths.size(0) != intrins.size(0))
        AT_ERROR("v_out_depths and intrins batch dimension mismatch");

    int b = v_out_depths.size(0), h = v_out_depths.size(1), w = v_out_depths.size(2);
    at::Tensor v_in_depths = at::empty_like(v_out_depths);

    ray_depth_to_linear_depth_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(v_out_depths), tensor2view<float, 4>(v_in_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_in_depths;
}


// ================
// Distort / Undistort
// ================

inline __device__ float get_pixel_bilinear(
    const TensorView<float, 4> image,  // [B, H, W, C]
    uint32_t bid,
    uint32_t cid,
    float x,
    float y,
    float padding = 0.0f
) {
    const long W = image.shape[2],
        H = image.shape[1];
    long x0 = (long)floorf(x);
    long x1 = x0 + 1;
    long y0 = (long)floorf(y);
    long y1 = y0 + 1;
    float wx1 = x - x0;
    float wx0 = 1.0f - wx1;
    float wy1 = y - y0;
    float wy0 = 1.0f - wy1;

    float c00 = (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, y0, x0, cid) : padding;
    float c10 = (x1 >= 0 && x1 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, y0, x1, cid) : padding;
    float c01 = (x0 >= 0 && x0 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, y1, x0, cid) : padding;
    float c11 = (x1 >= 0 && x1 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, y1, x1, cid) : padding;

    float c = 0.0f;
    c += c00 * (wx0 * wy0);
    c += c10 * (wx1 * wy0);
    c += c01 * (wx0 * wy1);
    c += c11 * (wx1 * wy1);
    return c;
}

inline __device__ float get_pixel_bilinear(
    const TensorView<float, 5> image,  // [B, K, H, W, C]
    uint32_t bid,
    uint32_t kid,
    uint32_t cid,
    float x,
    float y,
    float padding = 0.0f
) {
    const long W = image.shape[3],
        H = image.shape[2];
    long x0 = (long)floorf(x);
    long x1 = x0 + 1;
    long y0 = (long)floorf(y);
    long y1 = y0 + 1;
    float wx1 = x - x0;
    float wx0 = 1.0f - wx1;
    float wy1 = y - y0;
    float wy0 = 1.0f - wy1;

    float c00 = (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, kid, y0, x0, cid) : padding;
    float c10 = (x1 >= 0 && x1 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, kid, y0, x1, cid) : padding;
    float c01 = (x0 >= 0 && x0 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, kid, y1, x0, cid) : padding;
    float c11 = (x1 >= 0 && x1 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, kid, y1, x1, cid) : padding;

    float c = 0.0f;
    c += c00 * (wx0 * wy0);
    c += c10 * (wx1 * wy0);
    c += c01 * (wx0 * wy1);
    c += c11 * (wx1 * wy1);
    return c;
}

template<bool is_undistort>
__global__ void distort_image_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    TensorView<float, 4> out_image  // [B, H, W, C]
) {
    const int B = in_image.shape[0],
        H = in_image.shape[1],
        W = in_image.shape[2],
        C = in_image.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Load camera
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Undistort point
    float2 uv = { (i+0.5f-cx) / fx, (j+0.5f-cy) / fy };
    if (is_undistort) {
        if (dot(uv, uv) > 0.0f && !SlangProjectionUtils::is_valid_distortion(
            camera_model == ssplat::CameraModelType::FISHEYE ?
                normalize(uv) * atanf(length(uv)) : uv,
            dist_coeffs
        ))
            return;
        uv = SlangProjectionUtils::distort_point(uv, camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs);
    }
    else {
        if (!SlangProjectionUtils::undistort_point(uv, camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs, &uv))
            return;
    }

    // Process
    for (int c = 0; c < C; c++) {
        out_image.at(bid, j, i, c) = get_pixel_bilinear(in_image, bid, c, uv.x*fx+cx, uv.y*fy+cy);
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor distort_image_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (in_image.ndimension() != 4 || in_image.size(0) != intrins.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with intrins");

    int b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = zeros_like<float>(in_image);

    distort_image_kernel<false><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_image;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor undistort_image_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor in_image  // [B, H, W, C]
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(intrins);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (in_image.ndimension() != 4 || in_image.size(0) != intrins.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with intrins");

    int b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = zeros_like<float>(in_image);

    distort_image_kernel<true><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_image;
}


// ================
// Warp / Unwarp
// ================

__global__ void warp_image_wide_to_pinhole_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<float, 4> wide_image,  // [B, H, W, C]
    TensorView<float, 5> pinhole_images,  // [B*K, H, W, C]
    const float* __restrict__ axes  // [K, 3, 3]
) {
    const int B = wide_image.shape[0],
        K = pinhole_images.shape[1],
        Hp = pinhole_images.shape[2],
        Wp = pinhole_images.shape[3],
        C = wide_image.shape[3];

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp)
        return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    float4 intrin = intrins[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[0], axes[1], axes[2]};
        float3 axis_y = {axes[3], axes[4], axes[5]};
        float3 axis_z = {axes[6], axes[7], axes[8]};
        axes += 9;

        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv;
        bool valid = camera_model == ssplat::CameraModelType::FISHEYE ?
            SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        if (valid) {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) = get_pixel_bilinear(wide_image, bid, c, uv.x, uv.y, 0.5f);
        } else {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) = 0.5f;
        }
    }

}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_image_wide_to_pinhole_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor wide_image,  // [B, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(wide_image);
    CHECK_CUDA(wide_image);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (wide_image.ndimension() != 4 || wide_image.size(0) != intrins.size(0))
        AT_ERROR("wide_image shape must be (B, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3)
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = wide_image.size(0), c = wide_image.size(3);
    int k = axes.size(0);
    at::Tensor pinhole_images = at::empty({b, k, out_h, out_w, c}, wide_image.options());

    warp_image_wide_to_pinhole_kernel<<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(wide_image), tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return pinhole_images;
}

enum class WarpImageType {
    Default,
    LinearDepth,
    RayDepth,
    Points,
};

__forceinline__ __device__ float3 solve3(float3 col0, float3 col1, float3 col2, float3 vec) {
    float invdet = 1.0f / dot(cross(col0, col1), col2);
    float x = dot(cross(vec, col1), col2) * invdet;
    float y = dot(cross(col0, vec), col2) * invdet;
    float z = dot(cross(col0, col1), vec) * invdet;
    return {x, y, z};
}

template<WarpImageType type>
__global__ void warp_image_pinhole_to_wide_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<float, 4> wide_image,  // [B, H, W, C]
    TensorView<float, 5> pinhole_images,  // [B, K, H, W, C]
    const float* __restrict__ axes  // [K, 3, 3]
) {
    const int B = wide_image.shape[0],
        K = pinhole_images.shape[1],
        Hw = wide_image.shape[1],
        Ww = wide_image.shape[2],
        Hp = pinhole_images.shape[2],
        Wp = pinhole_images.shape[3],
        C = (type == WarpImageType::LinearDepth ? 1 :
             type == WarpImageType::Points ? 3 :
             wide_image.shape[3]);

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Ww || j >= Hw)
        return;

    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    float total[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float total_count = 0;

    float2 uv = { (i+0.5f-cx) / fx, (j+0.5f-cy) / fy };
    float3 raydir;
    if (!SlangProjectionUtils::unproject_point(uv, camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs, &raydir)) {
        for (int c = 0; c < C; c++)
            wide_image.at(bid, j, i, c) = 0.0f;
        return;
    }

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[0], axes[1], axes[2]};
        float3 axis_y = {axes[3], axes[4], axes[5]};
        float3 axis_z = {axes[6], axes[7], axes[8]};
        axes += 9;

        // axis_z + axis_x * tx + axis_y * ty = raydir * t
        // [axis_x axis_y raydir] [tx ty -t] = -axis_z

        float invdet = 1.0f / dot(cross(axis_x, axis_y), raydir);
        float tx = dot(cross(raydir, axis_y), axis_z) * invdet;
        float ty = dot(cross(axis_x, raydir), axis_z) * invdet;
        float t = dot(cross(axis_x, axis_y), axis_z) * invdet;
        if (fabsf(tx) >= 1.0f || fabsf(ty) >= 1.0f || t < 0.0f)
            continue;

        float2 xy = {(0.5f+0.5f*tx)*Wp, (0.5f+0.5f*ty)*Hp};

        float weight = 1.0f;
        float wx = length(axis_z) / length(axis_x);
        float wy = length(axis_z) / length(axis_y);
        if (wx < 1.0f && wy < 1.0f) {
            float ux = (1.0f - fabsf(tx)) / (1.0f - wx);
            float uy = (1.0f - fabsf(ty)) / (1.0f - wy);
            weight = fmaxf(fminf(fminf(ux, uy), 1.0f), 0.0f);
        }

        if (type == WarpImageType::LinearDepth || type == WarpImageType::RayDepth) {
            weight = weight*weight*(3.0f-2.0f*weight); // smoothstep
            // weight = weight*weight*weight*(weight*(6.0f*weight-15.0f)+10.0f);
            // constexpr float k = 3.0f;
            // weight = powf(weight, k) / (powf(weight, k) + powf(1.0f-weight, k));
            // TODO: this assums multi view consistent depth; Handle relative depth?
            float depth = get_pixel_bilinear(pinhole_images, bid, ki, 0, xy.x, xy.y);
            float3 point = normalize(raydir) * length(float3{uv.x, uv.y, 1.0f}) * depth;
            total[0] += weight * (type == WarpImageType::LinearDepth ? point.z : length(point));
        }
        else if (type == WarpImageType::Points) {
            // assums axes are orthogonal
            float3 ax = normalize(axis_x), ay = normalize(axis_y), az = normalize(axis_z);
            float3 point = float3{
                get_pixel_bilinear(pinhole_images, bid, ki, 0, xy.x, xy.y),
                get_pixel_bilinear(pinhole_images, bid, ki, 1, xy.x, xy.y),
                get_pixel_bilinear(pinhole_images, bid, ki, 2, xy.x, xy.y)
            };
            point = weight * float3{
                // dot(ax, point), dot(ay, point), dot(az, point)
                dot(float3{ax.x, ay.x, az.x}, point),
                dot(float3{ax.y, ay.y, az.y}, point),
                dot(float3{ax.z, ay.z, az.z}, point)
            };
            total[0] += point.x, total[1] += point.y, total[2] += point.z;
        }
        else {
            for (int c = 0; c < C; c++)
                total[c] += weight * get_pixel_bilinear(pinhole_images, bid, ki, c, uv.x, uv.y);
        }

        total_count += weight;
    }

    for (int c = 0; c < C; c++)
        wide_image.at(bid, j, i, c) = (total_count == 0.0f ? 0.0f : total[c] / total_count);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_image_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(pinhole_images);
    CHECK_CUDA(pinhole_images);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (pinhole_images.ndimension() != 5 || pinhole_images.size(0) != intrins.size(0))
        AT_ERROR("pinhole_images shape must be (B, K, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3 || axes.size(0) != pinhole_images.size(1))
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = pinhole_images.size(0), k = pinhole_images.size(1), c = pinhole_images.size(4);
    if (c > 4)
        AT_ERROR("warp_image_pinhole_to_wide supports max 4 channels");
    at::Tensor wide_image = at::empty({b, out_h, out_w, c}, pinhole_images.options());

    warp_image_pinhole_to_wide_kernel<WarpImageType::Default><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(wide_image), tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return wide_image;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_linear_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(pinhole_images);
    CHECK_CUDA(pinhole_images);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (pinhole_images.ndimension() != 5 || pinhole_images.size(0) != intrins.size(0))
        AT_ERROR("pinhole_images shape must be (B, K, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3 || axes.size(0) != pinhole_images.size(1))
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = pinhole_images.size(0), k = pinhole_images.size(1), c = pinhole_images.size(4);
    if (c != 1)
        AT_ERROR("depth map must have 1 channel");
    at::Tensor wide_image = at::empty({b, out_h, out_w, c}, pinhole_images.options());

    warp_image_pinhole_to_wide_kernel<WarpImageType::LinearDepth><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(wide_image), tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return wide_image;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_ray_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(pinhole_images);
    CHECK_CUDA(pinhole_images);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (pinhole_images.ndimension() != 5 || pinhole_images.size(0) != intrins.size(0))
        AT_ERROR("pinhole_images shape must be (B, K, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3 || axes.size(0) != pinhole_images.size(1))
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = pinhole_images.size(0), k = pinhole_images.size(1), c = pinhole_images.size(4);
    if (c != 1)
        AT_ERROR("depth map must have 1 channel");
    at::Tensor wide_image = at::empty({b, out_h, out_w, c}, pinhole_images.options());

    warp_image_pinhole_to_wide_kernel<WarpImageType::RayDepth><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(wide_image), tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return wide_image;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_points_pinhole_to_wide_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(pinhole_images);
    CHECK_CUDA(pinhole_images);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (pinhole_images.ndimension() != 5 || pinhole_images.size(0) != intrins.size(0))
        AT_ERROR("pinhole_images shape must be (B, K, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3 || axes.size(0) != pinhole_images.size(1))
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = pinhole_images.size(0), k = pinhole_images.size(1), c = pinhole_images.size(4);
    if (c != 3)
        AT_ERROR("point map must have 3 channels");
    at::Tensor wide_image = at::empty({b, out_h, out_w, c}, pinhole_images.options());

    warp_image_pinhole_to_wide_kernel<WarpImageType::Points><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(wide_image), tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return wide_image;
}


// Resolve scale in relative depth -
// The eigenvector corresponding to the smallest eigenvalue of this matrix tells how much depth maps need to be scaled
__global__ void warp_depth_pinhole_to_wide_scale_matrix_kernel(
    ssplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const int out_w, const int out_h,
    const TensorView<float, 5> pinhole_images,  // [B, K, H, W, C]
    const float* __restrict__ axes,  // [K, 3, 3]
    float* __restrict__ out_matrix  // [B, K, K]
) {
    const int B = pinhole_images.shape[0],
        K = pinhole_images.shape[1],
        Hw = out_h,
        Ww = out_w,
        Hp = pinhole_images.shape[2],
        Wp = pinhole_images.shape[3];

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Ww || j >= Hw)
        return;

    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    float2 uv = { (i+0.5f-cx) / fx, (j+0.5f-cy) / fy };
    float3 raydir;
    if (!SlangProjectionUtils::unproject_point(uv, camera_model == ssplat::CameraModelType::FISHEYE, dist_coeffs, &raydir))
        return;

    constexpr int MAX_K = 12;
    float depths[MAX_K];

    #pragma unroll MAX_K
    for (int ki = 0; ki < K; ++ki) {
        const float* axis_i = axes + 9 * ki;
        float3 axis_x = {axis_i[0], axis_i[1], axis_i[2]};
        float3 axis_y = {axis_i[3], axis_i[4], axis_i[5]};
        float3 axis_z = {axis_i[6], axis_i[7], axis_i[8]};

        // axis_z + axis_x * tx + axis_y * ty = raydir * t
        // [axis_x axis_y raydir] [tx ty -t] = -axis_z
        float invdet = 1.0f / dot(cross(axis_x, axis_y), raydir);
        float tx = dot(cross(raydir, axis_y), axis_z) * invdet;
        float ty = dot(cross(axis_x, raydir), axis_z) * invdet;
        float t = dot(cross(axis_x, axis_y), axis_z) * invdet;

        if (fabsf(tx) >= 1.0f || fabsf(ty) >= 1.0f || t < 0.0f) {
            depths[ki] = 0.0f;
        } else {
            float2 xy = {(0.5f+0.5f*tx)*Wp, (0.5f+0.5f*ty)*Hp};
            float depth = get_pixel_bilinear(pinhole_images, bid, ki, 0, xy.x, xy.y);
            // depths[ki] = depth;
            depths[ki] = __logf(fmaxf(depth, 1e-10f));
        }
    }

    out_matrix += bid * K*K;

    #pragma unroll MAX_K
    for (int i = 0; i < K; ++i)
        #pragma unroll MAX_K
        for (int j = 0; j < K; ++j) {
            float zi = depths[i], zj = depths[j];
            float w = (i == j ? 1.0f : -1.0f) * (zi * zj);
            atomicAddFVec<WARP_SIZE>(&out_matrix[K*i+j], w / (Ww*Hw));
        }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor warp_depth_pinhole_to_wide_scale_matrix_tensor(
    std::string camera_model,
    at::Tensor intrins,  // fx, fy, cx, cy
    CameraDistortionCoeffsTensor dist_coeffs,
    at::Tensor pinhole_images,  // [B, K, H, W, C]
    at::Tensor axes,  // [K, 3, 3]
    int out_w, int out_h
) {
    DEVICE_GUARD(pinhole_images);
    CHECK_CUDA(pinhole_images);
    CHECK_INPUT(intrins);
    CHECK_INPUT(axes);

    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    if (pinhole_images.ndimension() != 5 || pinhole_images.size(0) != intrins.size(0))
        AT_ERROR("pinhole_images shape must be (B, K, H, W, C) where B is consistent with intrins");
    if (axes.ndimension() != 3 || axes.size(1) != 3 || axes.size(2) != 3 || axes.size(0) != pinhole_images.size(1))
        AT_ERROR("axes shape must be (K, 3, 3)");

    int b = pinhole_images.size(0), k = pinhole_images.size(1), c = pinhole_images.size(4);
    if (c != 1)
        AT_ERROR("depth map must have 1 channel");
    constexpr int MAX_K = 12;
    if (k > MAX_K)
        AT_ERROR("k is too large, max " + std::to_string(MAX_K));
    at::Tensor matrix = at::zeros({b, k, k}, pinhole_images.options());

    warp_depth_pinhole_to_wide_scale_matrix_kernel<<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)intrins.data_ptr<float>(), dist_coeffs,
        out_w, out_h,
        tensor2view<float, 5>(pinhole_images),
        axes.data_ptr<float>(),
        matrix.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return matrix;
}


// ================
// PPISP
// ================

static constexpr int kNumPPISPParams = 36;
static constexpr int kNumPPISPParamsRQS = 39;

enum class PPISPParamType : int {
    Original,
    RQS,
};

template<PPISPParamType param_type>
__global__ void ppisp_forward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    TensorView<float, 4> out_image  // [B, H, W, C]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    static constexpr int kNumParams = (param_type == PPISPParamType::Original) ?
        kNumPPISPParams : kNumPPISPParamsRQS;
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[bid * kNumParams + i];
    }

    float3 pixel = in_image.load3(bid, y, x);

    float3 out_pixel;
    if (param_type == PPISPParamType::Original)
        out_pixel = SlangPPISP::apply_ppisp(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params)
        );
    else
        out_pixel = SlangPPISP::apply_ppisp_rqs(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params)
        );

    out_image.store3(bid, y, x, out_pixel);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ppisp_forward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    std::string param_type
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(intrins);
    if (in_image.ndimension() != 4 || in_image.size(0) != ppisp_params.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with ppisp_params");
    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    long b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = at::empty_like(in_image);
    if (param_type == "original" || param_type == "") {
        if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParams)
            AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS)");
        ppisp_forward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            tensor2view<float, 4>(in_image),
            ppisp_params.data_ptr<float>(),
            (float4*)intrins.data_ptr<float>(),
            actual_image_width,
            actual_image_height,
            tensor2view<float, 4>(out_image)
        );
    }
    else if (param_type == "rqs") {
        if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParamsRQS)
            AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS_RQS)");
        ppisp_forward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            tensor2view<float, 4>(in_image),
            ppisp_params.data_ptr<float>(),
            (float4*)intrins.data_ptr<float>(),
            actual_image_width,
            actual_image_height,
            tensor2view<float, 4>(out_image)
        );
    }
    else {
        AT_ERROR("invalid PPISP param_type, must be \"original\" or \"rqs\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
    return out_image;
}

template<PPISPParamType param_type>
__global__ void ppisp_backward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    const TensorView<float, 4> v_out_image,  // [B, H, W, C]
    TensorView<float, 4> v_in_image,  // [B, H, W, C]
    float* __restrict__ v_ppisp_params  // [B, PPISP_NUM_PARAMS]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 pixel = in_image.load3(bid, y, x);
    float3 v_out_pixel = v_out_image.load3(bid, y, x);

    static constexpr int kNumParams = (param_type == PPISPParamType::Original) ?
        kNumPPISPParams : kNumPPISPParamsRQS;
#if 0
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[bid * kNumParams + i];
    }
#else
    __shared__ float params_shared[kNumParams];
    if (threadIdx.x < kNumParams) {  // assume blockDim.x >= kNumParams
        float value = ppisp_params[bid * kNumParams + threadIdx.x];
        params_shared[threadIdx.x] = value;
    }
    __syncthreads();
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = params_shared[i];
        // int j = (i + threadIdx.x) % kNumParams;
        // params[j] = params_shared[j];
    }
#endif

    float3 v_pixel;
    FixedArray<float, kNumParams> v_params;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::apply_ppisp_vjp(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            v_out_pixel,
            &v_pixel,
            reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&v_params)
        );
    else
        SlangPPISP::apply_ppisp_rqs_vjp(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
            v_out_pixel,
            &v_pixel,
            reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&v_params)
        );

    v_in_image.store3(bid, y, x, v_pixel);

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        param = cg::reduce(warp, param, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && param != 0.0f)
            atomicAdd(&v_ppisp_params[bid * kNumParams + i], param);
    }
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor> ppisp_backward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    at::Tensor &v_out_image,  // [B, H, W, C]
    std::string param_type
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(intrins);
    CHECK_CUDA(v_out_image);
    if (in_image.ndimension() != 4 || in_image.size(0) != ppisp_params.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with ppisp_params");
    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    long b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor v_in_image = at::empty_like(in_image);
    at::Tensor v_ppisp_params = zeros_like<float>(ppisp_params);
    if (param_type == "original" || param_type == "") {
        if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParams)
            AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS)");
        ppisp_backward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            tensor2view<float, 4>(in_image),
            ppisp_params.data_ptr<float>(),
            (float4*)intrins.data_ptr<float>(),
            actual_image_width,
            actual_image_height,
            tensor2view<float, 4>(v_out_image),
            tensor2view<float, 4>(v_in_image),
            v_ppisp_params.data_ptr<float>()
        );
    }
    else if (param_type == "rqs") {
        if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParamsRQS)
            AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS_RQS)");
        ppisp_backward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            tensor2view<float, 4>(in_image),
            ppisp_params.data_ptr<float>(),
            (float4*)intrins.data_ptr<float>(),
            actual_image_width,
            actual_image_height,
            tensor2view<float, 4>(v_out_image),
            tensor2view<float, 4>(v_in_image),
            v_ppisp_params.data_ptr<float>()
        );
    }
    else {
        AT_ERROR("invalid PPISP param_type, must be \"original\" or \"rqs\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
    return std::make_tuple(v_in_image, v_ppisp_params);
}

template<PPISPParamType param_type>
__global__ void compute_raw_ppisp_regularization_forward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    float* __restrict__ raw_losses  // [B+1, RawPPISPRegLossIndex::length]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= B)
        return;

    static constexpr int kNumParams = (param_type == PPISPParamType::Original) ?
        kNumPPISPParams : kNumPPISPParamsRQS;
    static constexpr int kNumRawLosses = (param_type == PPISPParamType::Original) ?
        (int)RawPPISPRegLossIndex::length : (int)RawPPISPRegLossIndexRQS::length;

    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[bid * kNumParams + i];
    }

    FixedArray<float, kNumRawLosses> losses;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_raw_ppisp_regularization_loss(
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&losses)
        );
    else
        SlangPPISP::compute_raw_ppisp_rqs_regularization_loss(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&losses)
        );

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        float loss = isfinite(losses[i]) ? losses[i] : 0.0f;
        raw_losses[bid*kNumRawLosses + i] = loss;
        loss = cg::reduce(warp, loss, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && loss != 0.0f)
            atomicAdd(&raw_losses[B*kNumRawLosses + i], loss);
    }
}

template<PPISPParamType param_type>
__global__ void compute_ppisp_regularization_forward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    float* __restrict__ losses_buffer  // [PPISPRegLossIndex::length]
) {
    static constexpr int kNumRawLosses = (param_type == PPISPParamType::Original) ?
        (int)RawPPISPRegLossIndex::length : (int)RawPPISPRegLossIndexRQS::length;

    FixedArray<float, kNumRawLosses> raw_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (int)PPISPRegLossIndex::length> losses;

    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_ppisp_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );
    else
        SlangPPISP::compute_ppisp_rqs_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );

    #pragma unroll
    for (int i = 0; i < (int)PPISPRegLossIndex::length; i++) {
        losses_buffer[i] = losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
compute_ppsip_regularization_forward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    std::string param_type
) {
    DEVICE_GUARD(ppisp_params);
    CHECK_INPUT(ppisp_params);

    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = ppisp_params.size(0);

    at::Tensor raw_losses;
    at::Tensor losses = at::zeros({(int)PPISPRegLossIndex::length}, ppisp_params.options());

    if (param_type == "original" || param_type == "") {
        raw_losses = at::zeros({B+1, (int)RawPPISPRegLossIndex::length}, ppisp_params.options());
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            ppisp_params.data_ptr<float>(),
            raw_losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            raw_losses.data_ptr<float>() + B*(int)RawPPISPRegLossIndex::length,
            loss_weights,
            losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        raw_losses = at::zeros({B+1, (int)RawPPISPRegLossIndexRQS::length}, ppisp_params.options());
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            ppisp_params.data_ptr<float>(),
            raw_losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            raw_losses.data_ptr<float>() + B*(int)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        AT_ERROR("invalid PPISP param_type, must be \"original\" or \"rqs\"");
    }

    return std::make_tuple(losses, raw_losses);
}

template<PPISPParamType param_type>
__global__ void compute_raw_ppisp_regularization_backward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const float* __restrict__ v_raw_losses,  // [RawPPISPRegLossIndex::length]
    float* __restrict__ v_ppisp_params  // [B, PPISP_NUM_PARAMS]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= B)
        return;

    static constexpr int kNumParams = (param_type == PPISPParamType::Original) ?
        kNumPPISPParams : kNumPPISPParamsRQS;
    static constexpr int kNumRawLosses = (param_type == PPISPParamType::Original) ?
        (int)RawPPISPRegLossIndex::length : (int)RawPPISPRegLossIndexRQS::length;

    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[bid * kNumParams + i];
    }

    FixedArray<float, kNumRawLosses> v_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        v_losses[i] = v_raw_losses[i];
    }

    FixedArray<float, kNumParams> v_params;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_raw_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&v_params)
        );
    else
        SlangPPISP::compute_raw_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&v_params)
        );

    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        v_ppisp_params[bid * kNumParams + i] = param;
    }
}

template<PPISPParamType param_type>
__global__ void compute_ppisp_regularization_backward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    const float* __restrict__ v_losses_buffer,  // [PPISPRegLossIndex::length]
    float* __restrict__ v_raw_losses_buffer  // [RawPPISPRegLossIndex::length]
) {
    static constexpr int kNumRawLosses = (param_type == PPISPParamType::Original) ?
        (int)RawPPISPRegLossIndex::length : (int)RawPPISPRegLossIndexRQS::length;

    FixedArray<float, kNumRawLosses> raw_losses;
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (int)PPISPRegLossIndex::length> v_losses;
    #pragma unroll
    for (int i = 0; i < (int)PPISPRegLossIndex::length; i++) {
        v_losses[i] = v_losses_buffer[i];
    }

    FixedArray<float, kNumRawLosses> v_raw_losses;
    if (param_type == PPISPParamType::Original)
        SlangPPISP::compute_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_raw_losses)
        );
    else
        SlangPPISP::compute_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_raw_losses)
        );

    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        v_raw_losses_buffer[i] = v_raw_losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor compute_ppsip_regularization_backward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    at::Tensor &raw_losses,  // [B+1, RawPPISPRegLossIndex::length]
    at::Tensor &v_losses,  // [PPISPRegLossIndex::length]
    std::string param_type
) {
    DEVICE_GUARD(ppisp_params);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(raw_losses);
    CHECK_INPUT(v_losses);

    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = ppisp_params.size(0);

    at::Tensor v_raw_losses;
    at::Tensor v_ppisp_params = zeros_like<float>(ppisp_params);

    if (param_type == "original" || param_type == "") {
        v_raw_losses = at::zeros({(int)RawPPISPRegLossIndex::length}, ppisp_params.options());

        compute_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            raw_losses.data_ptr<float>() + B*(uint)RawPPISPRegLossIndex::length,
            loss_weights,
            v_losses.data_ptr<float>(),
            v_raw_losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            ppisp_params.data_ptr<float>(),
            v_raw_losses.data_ptr<float>(),
            v_ppisp_params.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        v_raw_losses = at::zeros({(int)RawPPISPRegLossIndexRQS::length}, ppisp_params.options());

        compute_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            raw_losses.data_ptr<float>() + B*(uint)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            v_losses.data_ptr<float>(),
            v_raw_losses.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            ppisp_params.data_ptr<float>(),
            v_raw_losses.data_ptr<float>(),
            v_ppisp_params.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        AT_ERROR("invalid PPISP param_type, must be \"original\" or \"rqs\"");
    }

    return v_ppisp_params;
}
