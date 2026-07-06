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

#include "Common.cuh"


static TensorView<float, 4> _bhw1_view(const TorchTensorView& tv) {
    const auto& s = std::get<2>(tv);
    int64_t B = s[0], H = s[1], W = s[2];
    TensorView<float, 4> v;
    v.data = (float*)std::get<0>(tv);
    v.shape[0] = B; v.shape[1] = H; v.shape[2] = W; v.shape[3] = 1;
    v.strides[0] = H*W; v.strides[1] = W; v.strides[2] = 1; v.strides[3] = 1;
    return v;
}

// Convert DeviceTensor3D<T> to contiguous TensorView<U, 4> where T packs sizeof(T)/sizeof(U) channels.
template<typename U, typename T>
static TensorView<U, 4> _dt3d_to_tv4(const DeviceTensor3D<T>& dt) {
    static_assert(sizeof(T) % sizeof(U) == 0, "T must be a multiple of U");
    constexpr int C = sizeof(T) / sizeof(U);
    TensorView<U, 4> v;
    v.data = (U*)dt.data_ptr();
    v.shape[0] = dt.template size<0>();
    v.shape[1] = dt.template size<1>();
    v.shape[2] = dt.template size<2>();
    v.shape[3] = C;
    long HW = v.shape[1] * v.shape[2];
    v.strides[0] = HW * C;
    v.strides[1] = v.shape[2] * C;
    v.strides[2] = C;
    v.strides[3] = 1;
    return v;
}

// Convert DeviceTensor2D<T> to contiguous TensorView<U, 3>.
template<typename U, typename T>
static TensorView<U, 3> _dt2d_to_tv3(const DeviceTensor2D<T>& dt) {
    static_assert(sizeof(T) % sizeof(U) == 0, "T must be a multiple of U");
    constexpr int C = sizeof(T) / sizeof(U);
    TensorView<U, 3> v;
    v.data = (U*)dt.data_ptr();
    v.shape[0] = dt.template size<0>();
    v.shape[1] = dt.template size<1>();
    v.shape[2] = C;
    v.strides[0] = v.shape[1] * C;
    v.strides[1] = C;
    v.strides[2] = 1;
    return v;
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
void uint8_image_to_float_tensor(
    DeviceTensor3D<uint8_t> img_in,  // [B, H, W, C] (C packed as uint8_t)
    DeviceTensor3D<float> img_out    // [B, H, W, C]
) {
    long b = img_in.size<0>(), h = img_in.size<1>(), w = img_in.size<2>();

    uint8_image_to_float_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<uint8_t>(img_in),
        _dt3d_to_tv4<float>(img_out)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void uint16_image_to_float_tensor(
    DeviceTensor3D<uint16_t> img_in,  // [B, H, W, C] (C packed as uint16_t)
    DeviceTensor3D<float> img_out     // [B, H, W, C]
) {
    long b = img_in.size<0>(), h = img_in.size<1>(), w = img_in.size<2>();

    uint16_image_to_float_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<uint16_t>(img_in),
        _dt3d_to_tv4<float>(img_out)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

// Raw-pointer entry points so callers can drive the kernel without packing
// inputs into DeviceTensor3D (whose element-size constraint forces a particular
// shape carrier). Used by Engine.cpp to avoid an extra CPU float pass on
// gt_rgb: H2D-copy uint8 [B, H, W, C] and convert on the GPU into a float
// buffer of the same layout.
void uint8_image_to_float_raw(
    const uint8_t* d_in, float* d_out,
    int B, int H, int W, int C
) {
    TensorView<uint8_t, 4> in_v;
    in_v.data = const_cast<uint8_t*>(d_in);
    in_v.shape[0] = B; in_v.shape[1] = H; in_v.shape[2] = W; in_v.shape[3] = C;
    in_v.strides[0] = (long)H*W*C; in_v.strides[1] = (long)W*C; in_v.strides[2] = C; in_v.strides[3] = 1;

    TensorView<float, 4> out_v;
    out_v.data = d_out;
    out_v.shape[0] = B; out_v.shape[1] = H; out_v.shape[2] = W; out_v.shape[3] = C;
    out_v.strides[0] = (long)H*W*C; out_v.strides[1] = (long)W*C; out_v.strides[2] = C; out_v.strides[3] = 1;

    uint8_image_to_float_kernel<<<_LAUNCH_ARGS_2D(H*W, B, 256, 1)>>>(in_v, out_v);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

void uint16_image_to_float_raw(
    const uint16_t* d_in, float* d_out,
    int B, int H, int W, int C
) {
    TensorView<uint16_t, 4> in_v;
    in_v.data = const_cast<uint16_t*>(d_in);
    in_v.shape[0] = B; in_v.shape[1] = H; in_v.shape[2] = W; in_v.shape[3] = C;
    in_v.strides[0] = (long)H*W*C; in_v.strides[1] = (long)W*C; in_v.strides[2] = C; in_v.strides[3] = 1;

    TensorView<float, 4> out_v;
    out_v.data = d_out;
    out_v.shape[0] = B; out_v.shape[1] = H; out_v.shape[2] = W; out_v.shape[3] = C;
    out_v.strides[0] = (long)H*W*C; out_v.strides[1] = (long)W*C; out_v.strides[2] = C; out_v.strides[3] = 1;

    uint16_image_to_float_kernel<<<_LAUNCH_ARGS_2D(H*W, B, 256, 1)>>>(in_v, out_v);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

// ---- gt_normal: uint8 [0,255] -> float in [-1, 1] ----
// Mirrors model.py's `gt_normal.float() / (255/2) - 1.0`. Per-channel:
//   out = (float)in / 127.5f - 1.0f
__global__ void uint8_normal_to_float_kernel(
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
        float c_out = (float)c_in * (1.0f / 127.5f) - 1.0f;
        img_out.at(bid, y, x, i) = c_out;
    }
}

void uint8_normal_to_float_raw(
    const uint8_t* d_in, float* d_out,
    int B, int H, int W, int C
) {
    TensorView<uint8_t, 4> in_v;
    in_v.data = const_cast<uint8_t*>(d_in);
    in_v.shape[0] = B; in_v.shape[1] = H; in_v.shape[2] = W; in_v.shape[3] = C;
    in_v.strides[0] = (long)H*W*C; in_v.strides[1] = (long)W*C; in_v.strides[2] = C; in_v.strides[3] = 1;

    TensorView<float, 4> out_v;
    out_v.data = d_out;
    out_v.shape[0] = B; out_v.shape[1] = H; out_v.shape[2] = W; out_v.shape[3] = C;
    out_v.strides[0] = (long)H*W*C; out_v.strides[1] = (long)W*C; out_v.strides[2] = C; out_v.strides[3] = 1;

    uint8_normal_to_float_kernel<<<_LAUNCH_ARGS_2D(H*W, B, 256, 1)>>>(in_v, out_v);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

// ---- gt_depth: uint16 -> float (cast only, no scaling) ----
// Mirrors model.py's `gt_depth.float()` (which preserves the raw integer value).
__global__ void uint16_depth_to_float_kernel(
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
        img_out.at(bid, y, x, i) = (float)img_in.at(bid, y, x, i);
    }
}

void uint16_depth_to_float_raw(
    const uint16_t* d_in, float* d_out,
    int B, int H, int W, int C
) {
    TensorView<uint16_t, 4> in_v;
    in_v.data = const_cast<uint16_t*>(d_in);
    in_v.shape[0] = B; in_v.shape[1] = H; in_v.shape[2] = W; in_v.shape[3] = C;
    in_v.strides[0] = (long)H*W*C; in_v.strides[1] = (long)W*C; in_v.strides[2] = C; in_v.strides[3] = 1;

    TensorView<float, 4> out_v;
    out_v.data = d_out;
    out_v.shape[0] = B; out_v.shape[1] = H; out_v.shape[2] = W; out_v.shape[3] = C;
    out_v.strides[0] = (long)H*W*C; out_v.strides[1] = (long)W*C; out_v.strides[2] = C; out_v.strides[3] = 1;

    uint16_depth_to_float_kernel<<<_LAUNCH_ARGS_2D(H*W, B, 256, 1)>>>(in_v, out_v);
    CHECK_DEVICE_ERROR(cudaGetLastError());
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
void rendered_depth_to_expected_depth_forward(
    TorchTensorView depth,  // [B, H, W, 1]
    TorchTensorView transmittance,  // [B, H, W, 1]
    TorchTensorView out_depth  // [B, H, W, 1]
) {
    const auto& s = std::get<2>(depth);
    int64_t b = s[0], h = s[1], w = s[2];

    float* max_depth = DevicePool::global().acquire<float>(PoolSlot::RdtedMaxDepth, b);
    cudaMemset(max_depth, 0, b * sizeof(float));

    rendered_depth_to_expected_depth_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _bhw1_view(depth), _bhw1_view(transmittance),
        _bhw1_view(out_depth), max_depth
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    rendered_depth_to_expected_depth_filter_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _bhw1_view(transmittance),
        _bhw1_view(out_depth),
        max_depth
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void rendered_depth_to_expected_depth_backward(
    TorchTensorView depth,  // [B, H, W, 1]
    TorchTensorView transmittance,  // [B, H, W, 1]
    TorchTensorView v_out_depth,  // [B, H, W, 1]
    TorchTensorView v_depth,  // [B, H, W, 1]
    TorchTensorView v_transmittance  // [B, H, W, 1]
) {
    const auto& s = std::get<2>(depth);
    int64_t b = s[0], h = s[1], w = s[2];

    rendered_depth_to_expected_depth_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _bhw1_view(depth), _bhw1_view(transmittance),
        _bhw1_view(v_out_depth),
        _bhw1_view(v_depth), _bhw1_view(v_transmittance)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
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
void blend_background_forward(
    DeviceTensor3D<float3> rgb,           // [B, H, W, 3]
    DeviceTensor3D<float>  transmittance, // [B, H, W, 1]
    DeviceTensor3D<float3> background,    // [B, H, W, 3]
    DeviceTensor3D<float3> out_rgb        // [B, H, W, 3]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    blend_background_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb), _dt3d_to_tv4<float>(transmittance), _dt3d_to_tv4<float>(background),
        _dt3d_to_tv4<float>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void blend_background_backward(
    DeviceTensor3D<float3> rgb,              // [B, H, W, 3]
    DeviceTensor3D<float>  transmittance,    // [B, H, W, 1]
    DeviceTensor3D<float3> background,       // [B, H, W, 3]
    DeviceTensor3D<float3> v_out_rgb,        // [B, H, W, 3]
    DeviceTensor3D<float3> v_rgb,            // [B, H, W, 3]
    DeviceTensor3D<float>  v_transmittance,  // [B, H, W, 1]
    DeviceTensor3D<float3> v_background      // [B, H, W, 3]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    blend_background_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb), _dt3d_to_tv4<float>(transmittance), _dt3d_to_tv4<float>(background),
        _dt3d_to_tv4<float>(v_out_rgb),
        _dt3d_to_tv4<float>(v_rgb), _dt3d_to_tv4<float>(v_transmittance), _dt3d_to_tv4<float>(v_background)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Blend Background with Random Noise
// ================

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
void blend_background_noise_forward(
    bool is_linear,
    DeviceTensor3D<float3> rgb,           // [B, H, W, 3]
    DeviceTensor3D<float>  transmittance, // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed,
    DeviceTensor3D<float3> out_rgb        // [B, H, W, 3]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    (is_linear ? blend_background_noise_forward_kernel<true> : blend_background_noise_forward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb), _dt3d_to_tv4<float>(transmittance),
        randomize_weight, seed,
        _dt3d_to_tv4<float>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void blend_background_noise_backward(
    bool is_linear,
    DeviceTensor3D<float3> rgb,              // [B, H, W, 3]
    DeviceTensor3D<float>  transmittance,    // [B, H, W, 1]
    float randomize_weight,
    uint32_t seed,
    DeviceTensor3D<float3> v_out_rgb,        // [B, H, W, 3]
    DeviceTensor3D<float3> v_rgb,            // [B, H, W, 3]
    DeviceTensor3D<float>  v_transmittance   // [B, H, W, 1]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    (is_linear ? blend_background_noise_backward_kernel<true> : blend_background_noise_backward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb), _dt3d_to_tv4<float>(transmittance),
        randomize_weight, seed,
        _dt3d_to_tv4<float>(v_out_rgb),
        _dt3d_to_tv4<float>(v_rgb), _dt3d_to_tv4<float>(v_transmittance)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
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
void rgb_to_srgb_forward(
    bool is_input_linear,
    DeviceTensor3D<float3> rgb,          // [B, H, W, 3]
    DeviceTensor2D<float3> color_matrix, // [3, 3] stored as 3 float3
    DeviceTensor3D<float3> out_rgb       // [B, H, W, 3]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    (is_input_linear ? rgb_to_srgb_forward_kernel<true> : rgb_to_srgb_forward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb),
        (float*)color_matrix.data_ptr(),
        _dt3d_to_tv4<float>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void rgb_to_srgb_backward(
    bool is_input_linear,
    DeviceTensor3D<float3> rgb,          // [B, H, W, 3]
    DeviceTensor2D<float3> color_matrix, // [3, 3] stored as 3 float3
    DeviceTensor3D<float3> v_out_rgb,    // [B, H, W, 3]
    DeviceTensor3D<float3> v_rgb         // [B, H, W, 3]
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();

    (is_input_linear ? rgb_to_srgb_backward_kernel<true> : rgb_to_srgb_backward_kernel<false>)
    <<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb),
        (float*)color_matrix.data_ptr(),
        _dt3d_to_tv4<float>(v_out_rgb),
        _dt3d_to_tv4<float>(v_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Overexposure Regularization
// ================

// Penalizes rgb values outside [0,1] with an L2 loss
//   L = weight * mean_{b,y,x,c}( max(-x, x-1, 0)^2 )
// The actual scalar loss is never materialized; this kernel only adds the
// per-pixel gradient dL/dx into v_rgb in-place. With N = B*H*W*3 the
// per-channel grad is
//   x < 0  : 2 * weight * x / N
//   x > 1  : 2 * weight * (x - 1) / N
//   else   : 0
__global__ void overexposure_grad_add_kernel(
    const TensorView<float, 4> rgb,
    const float scale,                // 2 * weight / (B * H * W * 3)
    TensorView<float, 4> v_rgb
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = rgb.shape[0], H = rgb.shape[1], W = rgb.shape[2];
    if (bid >= B || gid >= H * W) return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    float3 c   = rgb.load3(bid, y, x);
    float3 v   = v_rgb.load3(bid, y, x);

    auto over = [scale](float v_in, float v_pix) -> float {
        float g = (v_pix < 0.0f) ? v_pix
                : (v_pix > 1.0f) ? (v_pix - 1.0f)
                : 0.0f;
        return v_in + scale * g;
    };
    v.x = over(v.x, c.x);
    v.y = over(v.y, c.y);
    v.z = over(v.z, c.z);

    v_rgb.store3(bid, y, x, v);
}

/*[AutoHeaderGeneratorExport]*/
void overexposure_grad_add(
    DeviceTensor3D<float3> rgb,    // [B, H, W, 3]
    float weight,                  // L = weight * mean(max(-x, x-1, 0)^2)
    DeviceTensor3D<float3> v_rgb   // [B, H, W, 3], in/out
) {
    long b = rgb.size<0>(), h = rgb.size<1>(), w = rgb.size<2>();
    if (b <= 0 || h <= 0 || w <= 0 || weight == 0.0f) return;
    double N = (double)b * (double)h * (double)w * 3.0;
    float scale = (float)(2.0 * (double)weight / N);

    overexposure_grad_add_kernel<<<_LAUNCH_ARGS_2D(h * w, b, 256, 1)>>>(
        _dt3d_to_tv4<float>(rgb),
        scale,
        _dt3d_to_tv4<float>(v_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Depth to Points
// ================


__global__ void depth_to_points_forward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model,
        is_ray_depth,
        in_depth
    );
    out_points.store3(bid, j, i, out_point);
}

__global__ void depth_to_points_backward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model,
        is_ray_depth,
        in_depth, v_out_point
    );
    v_in_depths.store1(bid, j, i, v_in_depth);
}

/*[AutoHeaderGeneratorExport]*/
void depth_to_points_forward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  depths,      // [B, H, W, 1]
    DeviceTensor3D<float3> out_points   // [B, H, W, 3]
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();

    depth_to_points_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs, is_ray_depth,
        _dt3d_to_tv4<float>(depths), _dt3d_to_tv4<float>(out_points)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void depth_to_points_backward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  in_depths,   // [B, H, W, 1]
    DeviceTensor3D<float3> v_out_points,// [B, H, W, 3]
    DeviceTensor3D<float>  v_in_depths  // [B, H, W, 1]
) {
    int b = in_depths.size<0>(), h = in_depths.size<1>(), w = in_depths.size<2>();

    depth_to_points_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs, is_ray_depth,
        _dt3d_to_tv4<float>(in_depths), _dt3d_to_tv4<float>(v_out_points),
        _dt3d_to_tv4<float>(v_in_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Depth to Normal
// ================


__global__ void depth_to_normal_forward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model, is_ray_depth,
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
            (int)camera_model, is_ray_depth
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
    CameraModelType camera_model,
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
        (int)camera_model, is_ray_depth,
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
            (int)camera_model, is_ray_depth
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
void depth_to_normal_forward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  depths,      // [B, H, W, 1]
    DeviceTensor3D<float3> normals      // [B, H, W, 3]
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();

    depth_to_normal_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        is_ray_depth, _dt3d_to_tv4<float>(depths), _dt3d_to_tv4<float>(normals)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void depth_to_normal_backward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  depths,      // [B, H, W, 1]
    DeviceTensor3D<float3> v_normals,   // [B, H, W, 3]
    DeviceTensor3D<float>  v_depths     // [B, H, W, 1] (accumulated in-place)
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();

    depth_to_normal_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        is_ray_depth, _dt3d_to_tv4<float>(depths),
        _dt3d_to_tv4<float>(v_normals),
        _dt3d_to_tv4<float>(v_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// Python-callable wrappers using only TorchTensorView (pybind11 can convert tuples).
/*[AutoHeaderGeneratorExport]*/
void depth_to_normal_forward_tv(
    std::string camera_model,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    bool is_ray_depth,
    TorchTensorView depths,    // [B, H, W, 1] float32, CUDA
    TorchTensorView normals    // [B, H, W, 3] float32, CUDA (pre-allocated output)
) {
    depth_to_normal_forward(
        camera_model, intrins, dist_coeffs, is_ray_depth,
        DeviceTensor3D<float>(depths),
        DeviceTensor3D<float3>(normals)
    );
}

/*[AutoHeaderGeneratorExport]*/
void depth_to_normal_backward_tv(
    std::string camera_model,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    bool is_ray_depth,
    TorchTensorView depths,    // [B, H, W, 1]
    TorchTensorView v_normals, // [B, H, W, 3]
    TorchTensorView v_depths   // [B, H, W, 1] accumulated in-place
) {
    depth_to_normal_backward(
        camera_model, intrins, dist_coeffs, is_ray_depth,
        DeviceTensor3D<float>(depths),
        DeviceTensor3D<float3>(v_normals),
        DeviceTensor3D<float>(v_depths)
    );
}


// ================
// Loss between Depth and Normal
// ================


__global__ void depth_normal_loss_forward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model, is_ray_depth,
        depth, gt_normal
    );
#else
    // TODO: pre-undistort and store in shared memory
#endif
    losses.store1(bid, j, i, loss);
}


__global__ void depth_normal_loss_backward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model, is_ray_depth,
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
void depth_normal_loss_forward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  depths,      // [B, H, W, 1]
    DeviceTensor3D<float3> gt_normals,  // [B, H, W, 3]
    DeviceTensor3D<float>  losses       // [B, H, W, 1]
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();

    depth_normal_loss_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        is_ray_depth, _dt3d_to_tv4<float>(depths), _dt3d_to_tv4<float>(gt_normals),
        _dt3d_to_tv4<float>(losses)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void depth_normal_loss_backward(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    bool is_ray_depth,
    DeviceTensor3D<float>  depths,      // [B, H, W, 1]
    DeviceTensor3D<float3> gt_normals,  // [B, H, W, 3]
    DeviceTensor3D<float>  v_losses,    // [B, H, W, 1]
    DeviceTensor3D<float>  v_depths,    // [B, H, W, 1] (must be pre-zeroed)
    DeviceTensor3D<float3> v_gt_normals // [B, H, W, 3]
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();

    depth_normal_loss_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        is_ray_depth, _dt3d_to_tv4<float>(depths), _dt3d_to_tv4<float>(gt_normals),
        _dt3d_to_tv4<float>(v_losses),
        _dt3d_to_tv4<float>(v_depths), _dt3d_to_tv4<float>(v_gt_normals)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Ray Depth To Linear Depth
// ================


__global__ void ray_depth_to_linear_depth_forward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model
    );
    out_depths.store1(bid, j, i, out_depth);
}

__global__ void ray_depth_to_linear_depth_backward_kernel(
    CameraModelType camera_model,
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
        (int)camera_model
    );
    float v_in_depth = factor * v_out_depth;
    v_in_depths.store1(bid, j, i, v_in_depth);
}

/*[AutoHeaderGeneratorExport]*/
void ray_depth_to_linear_depth_forward(
    std::string camera_model,
    TorchTensorView intrins,  // [B, 4]
    TorchTensorView dist_coeffs,  // [B, 10]
    TorchTensorView depths,  // [B, H, W, 1]
    TorchTensorView out_depths  // [B, H, W, 1]
) {
    const auto& s = std::get<2>(depths);
    int64_t b = s[0], h = s[1], w = s[2];

    ray_depth_to_linear_depth_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        _bhw1_view(depths), _bhw1_view(out_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void ray_depth_to_linear_depth_backward(
    std::string camera_model,
    TorchTensorView intrins,  // [B, 4]
    TorchTensorView dist_coeffs,  // [B, 10]
    TorchTensorView v_out_depths,  // [B, H, W, 1]
    TorchTensorView v_in_depths  // [B, H, W, 1]
) {
    const auto& s = std::get<2>(v_out_depths);
    int64_t b = s[0], h = s[1], w = s[2];

    ray_depth_to_linear_depth_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        _bhw1_view(v_out_depths), _bhw1_view(v_in_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// In-place conversion of LINEAR (z) depth to RAY depth (Euclidean distance
// along the camera ray). Used for supervision depth maps that store linear
// depth while the rasterizer renders ray depth. We reuse
// ray_depth_to_linear_depth_factor (linear = ray * factor, factor =
// sign(z)/length(raydir)); the inverse is ray = linear / factor. The
// intrinsics are at the *image* resolution; (sx, sy) rescale them to the depth
// map resolution (depth maps may differ in size from the rendered image). The
// zero sentinel (no GT) and degenerate rays (undistort failure -> factor 0)
// map to 0 (i.e. "no supervision here").
__global__ void linear_depth_to_ray_depth_inplace_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // image-res fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    float sx, float sy,                  // depth_res / image_res scale
    TensorView<float, 4> depths          // [B, Hd, Wd, 1] in/out
) {
    const int B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Intrinsics rescaled from image to depth-map resolution. Distortion
    // coeffs operate on normalized image coordinates and are resolution
    // independent, so they need no scaling.
    float4 intrin = intrins[bid];
    float fx = intrin.x * sx, fy = intrin.y * sy;
    float cx = intrin.z * sx, cy = intrin.w * sy;
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    float in_depth = depths.load1(bid, j, i);
    float factor = SlangPixelWise::ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
        (int)camera_model
    );
    float out_depth = (factor > 0.0f) ? (in_depth / factor) : 0.0f;
    depths.store1(bid, j, i, out_depth);
}

/*[AutoHeaderGeneratorExport]*/
void linear_depth_to_ray_depth_inplace(
    std::string camera_model,
    TorchTensorView intrins,        // [B, 4] at image resolution
    TorchTensorView dist_coeffs,    // [B, 10]
    int image_width, int image_height,
    DeviceTensor3D<float> depths    // [B, Hd, Wd, 1] in/out
) {
    int b = depths.size<0>(), h = depths.size<1>(), w = depths.size<2>();
    if (b <= 0 || h <= 0 || w <= 0) return;
    float sx = (image_width  > 0) ? (float)w / (float)image_width  : 1.0f;
    float sy = (image_height > 0) ? (float)h / (float)image_height : 1.0f;

    linear_depth_to_ray_depth_inplace_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        sx, sy, _dt3d_to_tv4<float>(depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Distort / Undistort
// ================

template<typename T>
inline __device__ T get_pixel_bilinear(
    const TensorView<T, 4> image,  // [B, H, W, C]
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
        (float)image.at(bid, y0, x0, cid) : padding;
    float c10 = (x1 >= 0 && x1 < W && y0 >= 0 && y0 < H) ?
        (float)image.at(bid, y0, x1, cid) : padding;
    float c01 = (x0 >= 0 && x0 < W && y1 >= 0 && y1 < H) ?
        (float)image.at(bid, y1, x0, cid) : padding;
    float c11 = (x1 >= 0 && x1 < W && y1 >= 0 && y1 < H) ?
        (float)image.at(bid, y1, x1, cid) : padding;

    float c = 0.0f;
    c += c00 * (wx0 * wy0);
    c += c10 * (wx1 * wy0);
    c += c01 * (wx0 * wy1);
    c += c11 * (wx1 * wy1);
    return (T)c;
}

template<typename T>
inline __device__ T get_pixel_bilinear(
    const TensorView<T, 5> image,  // [B, K, H, W, C]
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
        (float)image.at(bid, kid, y0, x0, cid) : padding;
    float c10 = (x1 >= 0 && x1 < W && y0 >= 0 && y0 < H) ?
        (float)image.at(bid, kid, y0, x1, cid) : padding;
    float c01 = (x0 >= 0 && x0 < W && y1 >= 0 && y1 < H) ?
        (float)image.at(bid, kid, y1, x0, cid) : padding;
    float c11 = (x1 >= 0 && x1 < W && y1 >= 0 && y1 < H) ?
        (float)image.at(bid, kid, y1, x1, cid) : padding;

    float c = 0.0f;
    c += c00 * (wx0 * wy0);
    c += c10 * (wx1 * wy0);
    c += c01 * (wx0 * wy1);
    c += c11 * (wx1 * wy1);
    return (T)c;
}

template<bool is_undistort>
__global__ void distort_image_kernel(
    CameraModelType camera_model,
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
            camera_model == CameraModelType::FISHEYE ?
                normalize(uv) * atanf(length(uv)) :
            camera_model == CameraModelType::EQUISOLID ?
                normalize(uv) * (2.0f * sinf(0.5f * atanf(length(uv)))) :
            uv,
            dist_coeffs
        ))
            return;
        uv = SlangProjectionUtils::distort_point(uv, (int)camera_model, dist_coeffs);
    }
    else {
        if (!SlangProjectionUtils::undistort_point(uv, (int)camera_model, dist_coeffs, &uv))
            return;
    }

    // Process
    for (int c = 0; c < C; c++) {
        out_image.at(bid, j, i, c) = get_pixel_bilinear(in_image, bid, c, uv.x*fx+cx, uv.y*fy+cy);
    }
}

// TorchTensorView -> contiguous TensorView<float, 4> for [B, H, W, C] images.
// Used by distort/undistort so the channel count is taken from the runtime shape
// rather than baked into a DeviceTensor3D<float> carrier (which would reject C != 1).
static TensorView<float, 4> _bhwc_view(const TorchTensorView& tv) {
    const auto& s = std::get<2>(tv);
    int64_t B = s[0], H = s[1], W = s[2], C = (s.size() >= 4 ? s[3] : 1);
    TensorView<float, 4> v;
    v.data = (float*)std::get<0>(tv);
    v.shape[0] = B; v.shape[1] = H; v.shape[2] = W; v.shape[3] = C;
    v.strides[0] = H*W*C; v.strides[1] = W*C; v.strides[2] = C; v.strides[3] = 1;
    return v;
}

/*[AutoHeaderGeneratorExport]*/
void distort_image_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView in_image,           // [B, H, W, C] float
    TorchTensorView out_image           // [B, H, W, C] float (must be pre-zeroed)
) {
    const auto& s = std::get<2>(in_image);
    int b = s[0], h = s[1], w = s[2];

    distort_image_kernel<false><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        _bhwc_view(in_image), _bhwc_view(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void undistort_image_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView in_image,           // [B, H, W, C] float
    TorchTensorView out_image           // [B, H, W, C] float (must be pre-zeroed)
) {
    const auto& s = std::get<2>(in_image);
    int b = s[0], h = s[1], w = s[2];

    distort_image_kernel<true><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        _bhwc_view(in_image), _bhwc_view(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Warp / Unwarp
// ================

template<typename T>
__global__ void warp_image_wide_to_pinhole_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<T, 4> wide_image,  // [B, H, W, C]
    TensorView<T, 5> pinhole_images,  // [B*K, H, W, C]
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
        bool valid = camera_model == CameraModelType::FISHEYE ?
                SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            camera_model == CameraModelType::EQUISOLID ?
                SlangProjectionUtils::equisolid_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        if (valid) {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) = get_pixel_bilinear<T>(wide_image, bid, c, uv.x, uv.y, 0.5f);
        } else {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) = 0.5f;
        }
    }

}

/*[AutoHeaderGeneratorExport]*/
void warp_image_wide_to_pinhole_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView wide_image,         // [B, H, W, C] (float)
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView pinhole_images      // [B, K, H, W, C] (float)
) {
    const auto& ws = std::get<2>(wide_image);
    int b = ws[0];
    // Construct TensorViews from TorchTensorView
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_wide_to_pinhole_kernel<<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        make_tv4(wide_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

template<typename T>
__global__ void warp_image_equirectangular_to_pinhole_kernel(
    TensorView<T, 4> equirectangular_image,  // [B, H, W, C]
    TensorView<T, 5> pinhole_images,  // [B*K, H, W, C]
    const float* __restrict__ axes,  // [K, 3, 3]
    float padding
) {
    const int B = equirectangular_image.shape[0],
        h = equirectangular_image.shape[1],
        w = equirectangular_image.shape[2],
        K = pinhole_images.shape[1],
        Hp = pinhole_images.shape[2],
        Wp = pinhole_images.shape[3],
        C = equirectangular_image.shape[3];

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp)
        return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    const float f = (float)w * (0.5f / (float)M_PI);
    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[0], axes[1], axes[2]};
        float3 axis_y = {axes[3], axes[4], axes[5]};
        float3 axis_z = {axes[6], axes[7], axes[8]};
        axes += 9;

        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv = {
            0.5f * (float)w + f * atan2f(raydir.x, raydir.z),
            0.5f * (float)h + f * atan2f(raydir.y, hypotf(raydir.x, raydir.z))
        };
        for (int c = 0; c < C; c++)
            pinhole_images.at(bid, ki, j, i, c) = get_pixel_bilinear(
                equirectangular_image, bid, c, uv.x, uv.y, padding);
    }

}

/*[AutoHeaderGeneratorExport]*/
void warp_image_equirectangular_to_pinhole_tensor(
    TorchTensorView equirectangular_image,  // [B, H, W, C] (float)
    TorchTensorView axes,                   // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView pinhole_images          // [B, K, H, W, C] (float)
) {
    const auto& es = std::get<2>(equirectangular_image);
    int b = es[0];
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_equirectangular_to_pinhole_kernel<<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        make_tv4(equirectangular_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes), 0.5f
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

// ================
// Byte-to-float fused warp kernels
// ================
//
// These read a uint8 / uint16 H2D staging buffer for the GT image and write
// directly into the engine's pre-split pinhole float buffer, fusing the
// byte->float conversion with the warp. This is the path used by the C++
// DataManager when warp_to_pinhole is enabled (for fisheye/equisolid) or
// when an equirectangular camera is present (always split). Saves an
// intermediate full-resolution float allocation vs the two-pass version
// (uint*_image_to_float_raw -> warp_image_*_tensor).
//
// Output values are normalized to [0, 1] via the `norm_inv` arg
// (1/255 for uint8, 1/65535 for uint16, 1.0 for float passthrough).

template<typename T_in>
inline __device__ float bilinear_byte_norm(
    const TensorView<T_in, 4>& image, uint32_t bid, uint32_t cid,
    float x, float y, float norm_inv, float padding_norm
) {
    const long W = image.shape[2], H = image.shape[1];
    long x0 = (long)floorf(x);
    long x1 = x0 + 1;
    long y0 = (long)floorf(y);
    long y1 = y0 + 1;
    float wx1 = x - x0, wx0 = 1.0f - wx1;
    float wy1 = y - y0, wy0 = 1.0f - wy1;
    auto fetch = [&](long y_, long x_) -> float {
        if (x_ < 0 || x_ >= W || y_ < 0 || y_ >= H) return padding_norm;
        return (float)image.at(bid, y_, x_, cid) * norm_inv;
    };
    return fetch(y0, x0) * (wx0 * wy0)
         + fetch(y0, x1) * (wx1 * wy0)
         + fetch(y1, x0) * (wx0 * wy1)
         + fetch(y1, x1) * (wx1 * wy1);
}

// Variant that wraps the x axis (for equirectangular sampling) and clamps y.
// At the panorama wrap u=0 / u=W the four bilinear taps stay valid -- without
// this the back-face of the cubemap (face 5) gets a full-height gray seam
// where atan2 flips sign across u=0 / u=W.
template<typename T_in>
inline __device__ float bilinear_byte_norm_wrap_u(
    const TensorView<T_in, 4>& image, uint32_t bid, uint32_t cid,
    float x, float y, float norm_inv, float padding_norm
) {
    const long W = image.shape[2], H = image.shape[1];
    long x0 = (long)floorf(x);
    long x1 = x0 + 1;
    long y0 = (long)floorf(y);
    long y1 = y0 + 1;
    float wx1 = x - x0, wx0 = 1.0f - wx1;
    float wy1 = y - y0, wy0 = 1.0f - wy1;
    auto wrap_x = [W](long v) {
        long r = v % W;
        if (r < 0) r += W;
        return r;
    };
    long x0w = wrap_x(x0), x1w = wrap_x(x1);
    auto fetch = [&](long y_, long xw) -> float {
        if (y_ < 0 || y_ >= H) return padding_norm;
        return (float)image.at(bid, y_, xw, cid) * norm_inv;
    };
    return fetch(y0, x0w) * (wx0 * wy0)
         + fetch(y0, x1w) * (wx1 * wy0)
         + fetch(y1, x0w) * (wx0 * wy1)
         + fetch(y1, x1w) * (wx1 * wy1);
}

// Wide -> pinhole, fused byte->float. Output is [B*K, H_out, W_out, C]
// laid out as [B, K, H_out, W_out, C] -- flatten via stride math.
template<typename T_in>
__global__ void warp_image_wide_to_pinhole_byte_to_float_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,                  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<T_in, 4>  wide_image,                     // [B, H, W, C]
    TensorView<float, 5> pinhole_images,                 // [B, K, H_out, W_out, C]
    const float* __restrict__ axes,                      // [K, 3, 3]
    float norm_inv                                       // 1/255 or 1/65535
) {
    const int B  = wide_image.shape[0],
              K  = pinhole_images.shape[1],
              Hp = pinhole_images.shape[2],
              Wp = pinhole_images.shape[3],
              C  = wide_image.shape[3];

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    float4 intrin = intrins[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv;
        bool valid = camera_model == CameraModelType::FISHEYE ?
                SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            camera_model == CameraModelType::EQUISOLID ?
                SlangProjectionUtils::equisolid_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        if (valid) {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) =
                    bilinear_byte_norm<T_in>(wide_image, bid, c, uv.x, uv.y, norm_inv, 0.5f);
        } else {
            for (int c = 0; c < C; c++)
                pinhole_images.at(bid, ki, j, i, c) = 0.5f;
        }
    }
}

// Equirectangular -> pinhole, fused byte->float.
template<typename T_in>
__global__ void warp_image_equirectangular_to_pinhole_byte_to_float_kernel(
    TensorView<T_in, 4>  equirectangular_image,           // [B, H, W, C]
    TensorView<float, 5> pinhole_images,                  // [B, K, H_out, W_out, C]
    const float* __restrict__ axes,                       // [K, 3, 3]
    float norm_inv
) {
    const int B = equirectangular_image.shape[0],
              h = equirectangular_image.shape[1],
              w = equirectangular_image.shape[2],
              K = pinhole_images.shape[1],
              Hp = pinhole_images.shape[2],
              Wp = pinhole_images.shape[3],
              C  = equirectangular_image.shape[3];

    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    const float f = (float)w * (0.5f / (float)M_PI);
    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv = {
            0.5f * (float)w + f * atan2f(raydir.x, raydir.z),
            0.5f * (float)h + f * atan2f(raydir.y, hypotf(raydir.x, raydir.z))
        };
        for (int c = 0; c < C; c++)
            pinhole_images.at(bid, ki, j, i, c) =
                bilinear_byte_norm_wrap_u<T_in>(
                    equirectangular_image, bid, c, uv.x, uv.y, norm_inv, 0.5f);
    }
}

// Nearest-neighbor mask warps. Mask is bool; we read uint8 (0 / nonzero)
// and emit uint8 0/1 at the post-split resolution. Camera-model dispatch
// matches the RGB kernels.
__global__ void warp_mask_wide_to_pinhole_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<uint8_t, 4> wide_mask,                    // [B, H, W, 1]
    TensorView<uint8_t, 5> pinhole_masks,                // [B, K, H_out, W_out, 1]
    const float* __restrict__ axes
) {
    const int B  = wide_mask.shape[0],
              H  = wide_mask.shape[1],
              W  = wide_mask.shape[2],
              K  = pinhole_masks.shape[1],
              Hp = pinhole_masks.shape[2],
              Wp = pinhole_masks.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    float4 intrin = intrins[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv;
        bool valid = camera_model == CameraModelType::FISHEYE ?
                SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            camera_model == CameraModelType::EQUISOLID ?
                SlangProjectionUtils::equisolid_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        uint8_t out = 0;
        if (valid) {
            int xs = (int)floorf(uv.x + 0.5f);
            int ys = (int)floorf(uv.y + 0.5f);
            if (xs >= 0 && xs < W && ys >= 0 && ys < H) {
                out = (wide_mask.at(bid, ys, xs, 0) != 0) ? 1 : 0;
            }
        }
        pinhole_masks.at(bid, ki, j, i, 0) = out;
    }
}

__global__ void warp_mask_equirectangular_to_pinhole_kernel(
    TensorView<uint8_t, 4> equi_mask,                    // [B, H, W, 1]
    TensorView<uint8_t, 5> pinhole_masks,                // [B, K, H_out, W_out, 1]
    const float* __restrict__ axes
) {
    const int B  = equi_mask.shape[0],
              h  = equi_mask.shape[1],
              w  = equi_mask.shape[2],
              K  = pinhole_masks.shape[1],
              Hp = pinhole_masks.shape[2],
              Wp = pinhole_masks.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    const float f = (float)w * (0.5f / (float)M_PI);
    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float u = 0.5f * (float)w + f * atan2f(raydir.x, raydir.z);
        float v = 0.5f * (float)h + f * atan2f(raydir.y, hypotf(raydir.x, raydir.z));
        int xs = (int)floorf(u + 0.5f);
        int ys = (int)floorf(v + 0.5f);
        // wrap u for equirectangular (panoramic continuity); clamp v.
        if (xs < 0) xs += w; else if (xs >= w) xs -= w;
        if (ys < 0) ys = 0; if (ys >= h) ys = h - 1;
        pinhole_masks.at(bid, ki, j, i, 0) =
            (equi_mask.at(bid, ys, xs, 0) != 0) ? 1 : 0;
    }
}


// ---- Host launchers (raw-pointer, since the engine drives these without
// going through DeviceTensor wrappers). Picks the right kernel by type tag.

/*[AutoHeaderGeneratorExport]*/
void launch_warp_byte_to_float_wide(
    std::string camera_model,
    const float* d_intrins,                // [B, 4]
    const float* d_dist_coeffs,            // [B, 10] (nullable -> all-zeros)
    const void* d_byte, bool input_is_u16,
    int B, int Hin, int Win, int C,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    auto make_5 = [&](float* p) {
        TensorView<float, 5> v;
        v.data = p;
        v.shape[0] = B; v.shape[1] = K; v.shape[2] = Hout; v.shape[3] = Wout; v.shape[4] = C;
        long s4 = C, s3 = (long)Wout * s4, s2 = (long)Hout * s3, s1 = (long)K * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };
    auto make_4_u8 = [&](const uint8_t* p) {
        TensorView<uint8_t, 4> v;
        v.data = const_cast<uint8_t*>(p);
        v.shape[0] = B; v.shape[1] = Hin; v.shape[2] = Win; v.shape[3] = C;
        long s3 = C, s2 = (long)Win * s3, s1 = (long)Hin * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_4_u16 = [&](const uint16_t* p) {
        TensorView<uint16_t, 4> v;
        v.data = const_cast<uint16_t*>(p);
        v.shape[0] = B; v.shape[1] = Hin; v.shape[2] = Win; v.shape[3] = C;
        long s3 = C, s2 = (long)Win * s3, s1 = (long)Hin * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    CameraModelType cm = cmt(camera_model);
    CameraDistortionCoeffsBuffer dcb(const_cast<float*>(d_dist_coeffs));
    const float4* intrins_f4 = (const float4*)d_intrins;
    if (input_is_u16) {
        warp_image_wide_to_pinhole_byte_to_float_kernel<uint16_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                make_4_u16((const uint16_t*)d_byte), make_5(d_float_out),
                d_axes, 1.0f / 65535.0f);
    } else {
        warp_image_wide_to_pinhole_byte_to_float_kernel<uint8_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                make_4_u8((const uint8_t*)d_byte), make_5(d_float_out),
                d_axes, 1.0f / 255.0f);
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_byte_to_float_equi(
    const void* d_byte, bool input_is_u16,
    int B, int Hin, int Win, int C,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    auto make_5 = [&](float* p) {
        TensorView<float, 5> v;
        v.data = p;
        v.shape[0] = B; v.shape[1] = K; v.shape[2] = Hout; v.shape[3] = Wout; v.shape[4] = C;
        long s4 = C, s3 = (long)Wout * s4, s2 = (long)Hout * s3, s1 = (long)K * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };
    auto make_4_u8 = [&](const uint8_t* p) {
        TensorView<uint8_t, 4> v;
        v.data = const_cast<uint8_t*>(p);
        v.shape[0] = B; v.shape[1] = Hin; v.shape[2] = Win; v.shape[3] = C;
        long s3 = C, s2 = (long)Win * s3, s1 = (long)Hin * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_4_u16 = [&](const uint16_t* p) {
        TensorView<uint16_t, 4> v;
        v.data = const_cast<uint16_t*>(p);
        v.shape[0] = B; v.shape[1] = Hin; v.shape[2] = Win; v.shape[3] = C;
        long s3 = C, s2 = (long)Win * s3, s1 = (long)Hin * s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    if (input_is_u16) {
        warp_image_equirectangular_to_pinhole_byte_to_float_kernel<uint16_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                make_4_u16((const uint16_t*)d_byte), make_5(d_float_out),
                d_axes, 1.0f / 65535.0f);
    } else {
        warp_image_equirectangular_to_pinhole_byte_to_float_kernel<uint8_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                make_4_u8((const uint8_t*)d_byte), make_5(d_float_out),
                d_axes, 1.0f / 255.0f);
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_mask_wide(
    std::string camera_model,
    const float* d_intrins,                // [B, 4]
    const float* d_dist_coeffs,            // [B, 10] (nullable)
    const uint8_t* d_byte_mask,
    int B, int Hin, int Win,
    uint8_t* d_byte_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    TensorView<uint8_t, 4> in_v;
    in_v.data = const_cast<uint8_t*>(d_byte_mask);
    in_v.shape[0] = B; in_v.shape[1] = Hin; in_v.shape[2] = Win; in_v.shape[3] = 1;
    in_v.strides[0] = (long)Hin*Win; in_v.strides[1] = Win; in_v.strides[2] = 1; in_v.strides[3] = 1;

    TensorView<uint8_t, 5> out_v;
    out_v.data = d_byte_out;
    out_v.shape[0] = B; out_v.shape[1] = K; out_v.shape[2] = Hout; out_v.shape[3] = Wout; out_v.shape[4] = 1;
    long s4 = 1, s3 = Wout * s4, s2 = Hout * s3, s1 = K * s2;
    out_v.strides[0] = s1; out_v.strides[1] = s2; out_v.strides[2] = s3; out_v.strides[3] = s4; out_v.strides[4] = 1;

    CameraModelType cm = cmt(camera_model);
    CameraDistortionCoeffsBuffer dcb(const_cast<float*>(d_dist_coeffs));
    const float4* intrins_f4 = (const float4*)d_intrins;
    warp_mask_wide_to_pinhole_kernel<<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
        cm, intrins_f4, dcb, in_v, out_v, d_axes);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_mask_equi(
    const uint8_t* d_byte_mask,
    int B, int Hin, int Win,
    uint8_t* d_byte_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    TensorView<uint8_t, 4> in_v;
    in_v.data = const_cast<uint8_t*>(d_byte_mask);
    in_v.shape[0] = B; in_v.shape[1] = Hin; in_v.shape[2] = Win; in_v.shape[3] = 1;
    in_v.strides[0] = (long)Hin*Win; in_v.strides[1] = Win; in_v.strides[2] = 1; in_v.strides[3] = 1;

    TensorView<uint8_t, 5> out_v;
    out_v.data = d_byte_out;
    out_v.shape[0] = B; out_v.shape[1] = K; out_v.shape[2] = Hout; out_v.shape[3] = Wout; out_v.shape[4] = 1;
    long s4 = 1, s3 = Wout * s4, s2 = Hout * s3, s1 = K * s2;
    out_v.strides[0] = s1; out_v.strides[1] = s2; out_v.strides[2] = s3; out_v.strides[3] = s4; out_v.strides[4] = 1;

    warp_mask_equirectangular_to_pinhole_kernel<<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
        in_v, out_v, d_axes);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// GT depth / normal wide -> pinhole warps (for split-mode supervision)
// ================
//
// These mirror the RGB byte->float warp above, but produce the per-face GT
// geometry buffers consumed by the per-pixel depth / normal supervision loss.
// The output layout matches gt.depth / gt.normal: [B*K, out_H, out_W, 1|3]
// float at the POST-split (render) resolution, so the loss kernel's bilinear
// sampler collapses to a 1:1 read (same projection as the rendered face).

// Convert a sampled wide-frame depth value to the per-face RAY depth (= the
// distance |P| from the camera center, which is rotation-invariant, so the
// same value serves every face). The rasterizer renders ray depth, so this is
// the convention the loss expects. Ray-depth input is copied as-is; linear (z)
// input is divided by cos(theta) w.r.t. the INPUT optical axis. That cosine is
// non-positive for rays at / behind the input +z axis (the side / back faces
// of a wide capture), where linear depth is ill-defined -- those are dropped
// (0 = "no GT"). d <= 0 is the invalid sentinel and passes through unchanged.
__forceinline__ __device__ float _wide_depth_to_face_ray_depth(
    float d, float3 raydir, bool input_is_ray_depth
) {
    if (d <= 0.0f) return 0.0f;
    if (input_is_ray_depth) return d;
    float rl = length(raydir);
    float rz = raydir.z;
    if (rz <= 1e-6f * rl) return 0.0f;
    return d * rl / rz;
}

template<typename T_in>
__global__ void warp_depth_wide_to_pinhole_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,                  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<T_in, 4>  wide_depth,                     // [B, Hd, Wd, 1]
    TensorView<float, 5> pinhole_depth,                  // [B, K, H_out, W_out, 1]
    const float* __restrict__ axes,                      // [K, 3, 3]
    int in_H, int in_W,                                  // intrinsics ref res
    float norm_inv,
    bool input_is_ray_depth
) {
    const int B  = wide_depth.shape[0],
              Hd = wide_depth.shape[1],
              Wd = wide_depth.shape[2],
              K  = pinhole_depth.shape[1],
              Hp = pinhole_depth.shape[2],
              Wp = pinhole_depth.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    float4 intrin = intrins[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);
    // Projected pixel is in intrinsics-reference (RGB input) pixel space;
    // rescale to the GT depth map's own resolution.
    float sx = (float)Wd / (float)in_W, sy = (float)Hd / (float)in_H;

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv;
        bool valid = camera_model == CameraModelType::FISHEYE ?
                SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            camera_model == CameraModelType::EQUISOLID ?
                SlangProjectionUtils::equisolid_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        float out = 0.0f;
        if (valid) {
            float d = bilinear_byte_norm<T_in>(
                wide_depth, bid, 0, uv.x * sx, uv.y * sy, norm_inv, 0.0f);
            out = _wide_depth_to_face_ray_depth(d, raydir, input_is_ray_depth);
        }
        pinhole_depth.at(bid, ki, j, i, 0) = out;
    }
}

template<typename T_in>
__global__ void warp_depth_equirectangular_to_pinhole_kernel(
    TensorView<T_in, 4>  wide_depth,                      // [B, Hd, Wd, 1]
    TensorView<float, 5> pinhole_depth,                   // [B, K, H_out, W_out, 1]
    const float* __restrict__ axes,                       // [K, 3, 3]
    float norm_inv,
    bool input_is_ray_depth
) {
    const int B  = wide_depth.shape[0],
              h  = wide_depth.shape[1],
              w  = wide_depth.shape[2],
              K  = pinhole_depth.shape[1],
              Hp = pinhole_depth.shape[2],
              Wp = pinhole_depth.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    const float f = (float)w * (0.5f / (float)M_PI);
    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv = {
            0.5f * (float)w + f * atan2f(raydir.x, raydir.z),
            0.5f * (float)h + f * atan2f(raydir.y, hypotf(raydir.x, raydir.z))
        };
        float d = bilinear_byte_norm_wrap_u<T_in>(
            wide_depth, bid, 0, uv.x, uv.y, norm_inv, 0.0f);
        pinhole_depth.at(bid, ki, j, i, 0) =
            _wide_depth_to_face_ray_depth(d, raydir, input_is_ray_depth);
    }
}

// Sample a wide-frame normal (in the INPUT camera frame) and rotate it into a
// pinhole face's camera frame. `axis_*` are the face axes expressed in the
// input frame (unit + orthogonal per the cubemap tables); the face-frame
// components are the dot products with the normalized axes. Returns the
// (-1,-1,-1) invalid sentinel when the sample is "no data" or degenerate.
template<typename T_in>
__forceinline__ __device__ float3 _warp_one_normal(
    const TensorView<T_in, 4>& wide_normal, uint32_t bid,
    float x, float y, float norm_inv, float decode_off, bool wrap_u,
    float3 axis_x, float3 axis_y, float3 axis_z
) {
    float3 n;
    if (wrap_u) {
        n.x = bilinear_byte_norm_wrap_u<T_in>(wide_normal, bid, 0, x, y, norm_inv, 0.0f) + decode_off;
        n.y = bilinear_byte_norm_wrap_u<T_in>(wide_normal, bid, 1, x, y, norm_inv, 0.0f) + decode_off;
        n.z = bilinear_byte_norm_wrap_u<T_in>(wide_normal, bid, 2, x, y, norm_inv, 0.0f) + decode_off;
    } else {
        n.x = bilinear_byte_norm<T_in>(wide_normal, bid, 0, x, y, norm_inv, 0.0f) + decode_off;
        n.y = bilinear_byte_norm<T_in>(wide_normal, bid, 1, x, y, norm_inv, 0.0f) + decode_off;
        n.z = bilinear_byte_norm<T_in>(wide_normal, bid, 2, x, y, norm_inv, 0.0f) + decode_off;
    }
    // "no data" sentinel (byte 0 -> -1 per channel; sum <= -2.366 == invalid,
    // matching the loss kernel's ref_normal validity test).
    if (n.x + n.y + n.z <= -2.366f) return make_float3(-1.0f, -1.0f, -1.0f);
    float3 ax = normalize(axis_x), ay = normalize(axis_y), az = normalize(axis_z);
    float3 r = make_float3(dot(ax, n), dot(ay, n), dot(az, n));
    float rl = length(r);
    if (rl <= 1e-8f) return make_float3(-1.0f, -1.0f, -1.0f);
    return r / rl;
}

template<typename T_in>
__global__ void warp_normal_wide_to_pinhole_kernel(
    CameraModelType camera_model,
    const float4 *__restrict__ intrins,                  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    TensorView<T_in, 4>  wide_normal,                    // [B, Hn, Wn, 3]
    TensorView<float, 5> pinhole_normal,                 // [B, K, H_out, W_out, 3]
    const float* __restrict__ axes,                      // [K, 3, 3]
    int in_H, int in_W,
    float norm_inv, float decode_off
) {
    const int B  = wide_normal.shape[0],
              Hn = wide_normal.shape[1],
              Wn = wide_normal.shape[2],
              K  = pinhole_normal.shape[1],
              Hp = pinhole_normal.shape[2],
              Wp = pinhole_normal.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    float4 intrin = intrins[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);
    float sx = (float)Wn / (float)in_W, sy = (float)Hn / (float)in_H;

    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv;
        bool valid = camera_model == CameraModelType::FISHEYE ?
                SlangProjectionUtils::fisheye_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            camera_model == CameraModelType::EQUISOLID ?
                SlangProjectionUtils::equisolid_proj_nav(raydir, intrin, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(raydir, intrin, dist_coeffs, &uv);
        float3 nf = make_float3(-1.0f, -1.0f, -1.0f);
        if (valid) {
            nf = _warp_one_normal<T_in>(
                wide_normal, bid, uv.x * sx, uv.y * sy, norm_inv, decode_off,
                /*wrap_u=*/false, axis_x, axis_y, axis_z);
        }
        pinhole_normal.at(bid, ki, j, i, 0) = nf.x;
        pinhole_normal.at(bid, ki, j, i, 1) = nf.y;
        pinhole_normal.at(bid, ki, j, i, 2) = nf.z;
    }
}

template<typename T_in>
__global__ void warp_normal_equirectangular_to_pinhole_kernel(
    TensorView<T_in, 4>  wide_normal,                     // [B, Hn, Wn, 3]
    TensorView<float, 5> pinhole_normal,                  // [B, K, H_out, W_out, 3]
    const float* __restrict__ axes,                       // [K, 3, 3]
    float norm_inv, float decode_off
) {
    const int B  = wide_normal.shape[0],
              h  = wide_normal.shape[1],
              w  = wide_normal.shape[2],
              K  = pinhole_normal.shape[1],
              Hp = pinhole_normal.shape[2],
              Wp = pinhole_normal.shape[3];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i   = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j   = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= Wp || j >= Hp) return;
    float tx = -1.0f + 2.0f * ((float)i + 0.5f) / (float)Wp;
    float ty = -1.0f + 2.0f * ((float)j + 0.5f) / (float)Hp;

    const float f = (float)w * (0.5f / (float)M_PI);
    for (int ki = 0; ki < K; ++ki) {
        float3 axis_x = {axes[9*ki + 0], axes[9*ki + 1], axes[9*ki + 2]};
        float3 axis_y = {axes[9*ki + 3], axes[9*ki + 4], axes[9*ki + 5]};
        float3 axis_z = {axes[9*ki + 6], axes[9*ki + 7], axes[9*ki + 8]};
        float3 raydir = axis_z + tx * axis_x + ty * axis_y;
        float2 uv = {
            0.5f * (float)w + f * atan2f(raydir.x, raydir.z),
            0.5f * (float)h + f * atan2f(raydir.y, hypotf(raydir.x, raydir.z))
        };
        float3 nf = _warp_one_normal<T_in>(
            wide_normal, bid, uv.x, uv.y, norm_inv, decode_off,
            /*wrap_u=*/true, axis_x, axis_y, axis_z);
        pinhole_normal.at(bid, ki, j, i, 0) = nf.x;
        pinhole_normal.at(bid, ki, j, i, 1) = nf.y;
        pinhole_normal.at(bid, ki, j, i, 2) = nf.z;
    }
}


// ---- Host launchers for the GT depth / normal warps -------------------------

// Build a [B, Hin, Win, C] input TensorView over a raw device pointer.
template<typename T>
static inline TensorView<T, 4> _make_tv4_in(const T* p, int B, int Hin, int Win, int C) {
    TensorView<T, 4> v;
    v.data = const_cast<T*>(p);
    v.shape[0] = B; v.shape[1] = Hin; v.shape[2] = Win; v.shape[3] = C;
    long s3 = C, s2 = (long)Win * s3, s1 = (long)Hin * s2;
    v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
    return v;
}

static inline TensorView<float, 5> _make_tv5_out(
    float* p, int B, int K, int Hout, int Wout, int C) {
    TensorView<float, 5> v;
    v.data = p;
    v.shape[0] = B; v.shape[1] = K; v.shape[2] = Hout; v.shape[3] = Wout; v.shape[4] = C;
    long s4 = C, s3 = (long)Wout * s4, s2 = (long)Hout * s3, s1 = (long)K * s2;
    v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
    return v;
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_depth_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const void* d_depth, uint32_t elem_size,
    int B, int Hin, int Win,
    int in_H, int in_W,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes, bool input_is_ray_depth)
{
    CameraModelType cm = cmt(camera_model);
    CameraDistortionCoeffsBuffer dcb(const_cast<float*>(d_dist_coeffs));
    const float4* intrins_f4 = (const float4*)d_intrins;
    auto out_v = _make_tv5_out(d_float_out, B, K, Hout, Wout, 1);
    if (elem_size == 2) {
        warp_depth_wide_to_pinhole_kernel<uint16_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                _make_tv4_in<uint16_t>((const uint16_t*)d_depth, B, Hin, Win, 1),
                out_v, d_axes, in_H, in_W, 1.0f, input_is_ray_depth);
    } else if (elem_size == 4) {
        warp_depth_wide_to_pinhole_kernel<float>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                _make_tv4_in<float>((const float*)d_depth, B, Hin, Win, 1),
                out_v, d_axes, in_H, in_W, 1.0f, input_is_ray_depth);
    } else {
        throw std::runtime_error("launch_warp_depth_wide: depth must be uint16 or float32");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_depth_equi(
    const void* d_depth, uint32_t elem_size,
    int B, int Hin, int Win,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes, bool input_is_ray_depth)
{
    auto out_v = _make_tv5_out(d_float_out, B, K, Hout, Wout, 1);
    if (elem_size == 2) {
        warp_depth_equirectangular_to_pinhole_kernel<uint16_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                _make_tv4_in<uint16_t>((const uint16_t*)d_depth, B, Hin, Win, 1),
                out_v, d_axes, 1.0f, input_is_ray_depth);
    } else if (elem_size == 4) {
        warp_depth_equirectangular_to_pinhole_kernel<float>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                _make_tv4_in<float>((const float*)d_depth, B, Hin, Win, 1),
                out_v, d_axes, 1.0f, input_is_ray_depth);
    } else {
        throw std::runtime_error("launch_warp_depth_equi: depth must be uint16 or float32");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_normal_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const void* d_normal, uint32_t elem_size,
    int B, int Hin, int Win,
    int in_H, int in_W,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    CameraModelType cm = cmt(camera_model);
    CameraDistortionCoeffsBuffer dcb(const_cast<float*>(d_dist_coeffs));
    const float4* intrins_f4 = (const float4*)d_intrins;
    auto out_v = _make_tv5_out(d_float_out, B, K, Hout, Wout, 3);
    if (elem_size == 1) {
        warp_normal_wide_to_pinhole_kernel<uint8_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                _make_tv4_in<uint8_t>((const uint8_t*)d_normal, B, Hin, Win, 3),
                out_v, d_axes, in_H, in_W, 1.0f / 127.5f, -1.0f);
    } else if (elem_size == 4) {
        warp_normal_wide_to_pinhole_kernel<float>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                cm, intrins_f4, dcb,
                _make_tv4_in<float>((const float*)d_normal, B, Hin, Win, 3),
                out_v, d_axes, in_H, in_W, 1.0f, 0.0f);
    } else {
        throw std::runtime_error("launch_warp_normal_wide: normal must be uint8 or float32");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void launch_warp_normal_equi(
    const void* d_normal, uint32_t elem_size,
    int B, int Hin, int Win,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    auto out_v = _make_tv5_out(d_float_out, B, K, Hout, Wout, 3);
    if (elem_size == 1) {
        warp_normal_equirectangular_to_pinhole_kernel<uint8_t>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                _make_tv4_in<uint8_t>((const uint8_t*)d_normal, B, Hin, Win, 3),
                out_v, d_axes, 1.0f / 127.5f, -1.0f);
    } else if (elem_size == 4) {
        warp_normal_equirectangular_to_pinhole_kernel<float>
            <<<_LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1)>>>(
                _make_tv4_in<float>((const float*)d_normal, B, Hin, Win, 3),
                out_v, d_axes, 1.0f, 0.0f);
    } else {
        throw std::runtime_error("launch_warp_normal_equi: normal must be uint8 or float32");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
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
    CameraModelType camera_model,
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
    if (!SlangProjectionUtils::unproject_point(uv, (int)camera_model, dist_coeffs, &raydir)) {
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
void warp_image_pinhole_to_wide_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView pinhole_images,     // [B, K, H, W, C]
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView wide_image          // [B, H, W, C]
) {
    const auto& ps = std::get<2>(pinhole_images);
    int b = ps[0];
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_pinhole_to_wide_kernel<WarpImageType::Default><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        make_tv4(wide_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void warp_linear_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView pinhole_images,     // [B, K, H, W, 1]
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView wide_image          // [B, H, W, 1]
) {
    const auto& ps = std::get<2>(pinhole_images);
    int b = ps[0];
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_pinhole_to_wide_kernel<WarpImageType::LinearDepth><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        make_tv4(wide_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void warp_ray_depth_pinhole_to_wide_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView pinhole_images,     // [B, K, H, W, 1]
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView wide_image          // [B, H, W, 1]
) {
    const auto& ps = std::get<2>(pinhole_images);
    int b = ps[0];
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_pinhole_to_wide_kernel<WarpImageType::RayDepth><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        make_tv4(wide_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void warp_points_pinhole_to_wide_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView pinhole_images,     // [B, K, H, W, 3]
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView wide_image          // [B, H, W, 3]
) {
    const auto& ps = std::get<2>(pinhole_images);
    int b = ps[0];
    auto make_tv4 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 4> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3];
        long s3 = s[3], s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = 1;
        return v;
    };
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_image_pinhole_to_wide_kernel<WarpImageType::Points><<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        make_tv4(wide_image), make_tv5(pinhole_images),
        (float*)std::get<0>(axes)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// Resolve scale in relative depth -
// The eigenvector corresponding to the smallest eigenvalue of this matrix tells how much depth maps need to be scaled
__global__ void warp_depth_pinhole_to_wide_scale_matrix_kernel(
    CameraModelType camera_model,
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
    if (!SlangProjectionUtils::unproject_point(uv, (int)camera_model, dist_coeffs, &raydir))
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
void warp_depth_pinhole_to_wide_scale_matrix_tensor(
    std::string camera_model,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 10]
    TorchTensorView pinhole_images,     // [B, K, H, W, 1]
    TorchTensorView axes,               // [K, 3, 3]
    int out_w, int out_h,
    TorchTensorView matrix              // [B, K, K] (must be pre-zeroed)
) {
    const auto& ps = std::get<2>(pinhole_images);
    int b = ps[0];
    auto make_tv5 = [](const TorchTensorView& tv) {
        const auto& s = std::get<2>(tv);
        TensorView<float, 5> v;
        v.data = (float*)std::get<0>(tv);
        v.shape[0] = s[0]; v.shape[1] = s[1]; v.shape[2] = s[2]; v.shape[3] = s[3]; v.shape[4] = s[4];
        long s4 = s[4], s3 = s[3]*s4, s2 = s[2]*s3, s1 = s[1]*s2;
        v.strides[0] = s1; v.strides[1] = s2; v.strides[2] = s3; v.strides[3] = s4; v.strides[4] = 1;
        return v;
    };

    warp_depth_pinhole_to_wide_scale_matrix_kernel<<<_LAUNCH_ARGS_3D(out_w, out_h, b, 16, 16, 1)>>>(
        cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs,
        out_w, out_h,
        make_tv5(pinhole_images),
        (float*)std::get<0>(axes),
        (float*)std::get<0>(matrix)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// PPISP
// ================

static constexpr int kNumPPISPParams = 36;
static constexpr int kNumPPISPParamsRQS = 39;
static constexpr int kNumPPISPParamsNoCRF = 24;

enum class PPISPParamType : int {
    Original,
    RQS,
    NoCRF,
};

template<PPISPParamType param_type>
__global__ void ppisp_forward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [N_cam or B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    const int* __restrict__ cam_indices,  // [B], or nullptr -> identity
    TensorView<float, 4> out_image  // [B, H, W, C]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    if (bid >= B || gid >= H*W)
        return;
    unsigned y = gid / W;
    unsigned x = gid % W;

    int p_id = cam_indices ? cam_indices[bid] : (int)bid;

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[p_id * kNumParams + i];
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
    else if (param_type == PPISPParamType::RQS)
        out_pixel = SlangPPISP::apply_ppisp_rqs(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params)
        );
    else
        out_pixel = SlangPPISP::apply_ppisp_no_crf(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params)
        );

    out_image.store3(bid, y, x, out_pixel);
}

/*[AutoHeaderGeneratorExport]*/
void ppisp_forward(
    DeviceTensor3D<float3> in_image,    // [B, H, W, 3]
    TorchTensorView ppisp_params,       // [N_cam or B, PPISP_NUM_PARAMS]
    TorchTensorView intrins,            // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    std::string param_type,
    TorchTensorView cam_indices,        // [B] int32, or null -> identity (ppisp_params is [B,P])
    DeviceTensor3D<float3> out_image    // [B, H, W, 3]
) {
    long b = in_image.size<0>(), h = in_image.size<1>(), w = in_image.size<2>();
    const int* cam_idx_ptr = (std::get<0>(cam_indices) == 0) ?
        nullptr : (const int*)std::get<0>(cam_indices);
    if (param_type == "original" || param_type == "") {
        ppisp_forward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else if (param_type == "rqs") {
        ppisp_forward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else if (param_type == "no_crf") {
        ppisp_forward_kernel<PPISPParamType::NoCRF><<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(out_image)
        );
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

template<PPISPParamType param_type>
__global__ void ppisp_backward_kernel(
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    const float* __restrict__ ppisp_params,  // [N_cam or B, PPISP_NUM_PARAMS]
    const float4 *__restrict__ intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    const int* __restrict__ cam_indices,  // [B], or nullptr -> identity
    const TensorView<float, 4> v_out_image,  // [B, H, W, C]
    TensorView<float, 4> v_in_image,  // [B, H, W, C]
    float* __restrict__ v_ppisp_params  // [N_cam or B, PPISP_NUM_PARAMS]
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned bid = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned B = in_image.shape[0], H = in_image.shape[1], W = in_image.shape[2];
    bool inside = (bid < B) && (gid < H*W);
    unsigned y = inside ? gid / W : 0u;
    unsigned x = inside ? gid % W : 0u;

    int p_id = (bid < B) ? (cam_indices ? cam_indices[bid] : (int)bid) : 0;

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
#if 0
    FixedArray<float, kNumParams> params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        params[i] = ppisp_params[p_id * kNumParams + i];
    }
#else
    __shared__ float params_shared[kNumParams];
    if (threadIdx.x < kNumParams) {  // assume blockDim.x >= kNumParams
        float value = ppisp_params[p_id * kNumParams + threadIdx.x];
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

    float3 v_pixel = make_float3(0.0f, 0.0f, 0.0f);
    FixedArray<float, kNumParams> v_params;
    #pragma unroll
    for (int i = 0; i < kNumParams; i++)
        v_params[i] = 0.0f;

    if (inside) {
        float3 pixel = in_image.load3(bid, y, x);
        float3 v_out_pixel = v_out_image.load3(bid, y, x);
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
        else if (param_type == PPISPParamType::RQS)
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
        else
            SlangPPISP::apply_ppisp_no_crf_vjp(
                pixel,
                make_float2((float)x, (float)y),
                make_float2(intrins[bid].z, intrins[bid].w),
                make_float2(actual_image_width, actual_image_height),
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
                v_out_pixel,
                &v_pixel,
                reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&v_params)
            );

        v_in_image.store3(bid, y, x, v_pixel);
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        param = cg::reduce(warp, param, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && param != 0.0f)
            atomicAdd(&v_ppisp_params[p_id * kNumParams + i], param);
    }
}

/*[AutoHeaderGeneratorExport]*/
void ppisp_backward(
    DeviceTensor3D<float3> in_image,    // [B, H, W, 3]
    TorchTensorView ppisp_params,       // [N_cam or B, PPISP_NUM_PARAMS]
    TorchTensorView intrins,            // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    DeviceTensor3D<float3> v_out_image, // [B, H, W, 3]
    std::string param_type,
    TorchTensorView cam_indices,        // [B] int32, or null -> identity
    DeviceTensor3D<float3> v_in_image,  // [B, H, W, 3]
    TorchTensorView v_ppisp_params      // [N_cam or B, PPISP_NUM_PARAMS] (must be pre-zeroed)
) {
    long b = in_image.size<0>(), h = in_image.size<1>(), w = in_image.size<2>();
    const int* cam_idx_ptr = (std::get<0>(cam_indices) == 0) ?
        nullptr : (const int*)std::get<0>(cam_indices);
    if (param_type == "original" || param_type == "") {
        ppisp_backward_kernel<PPISPParamType::Original><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else if (param_type == "rqs") {
        ppisp_backward_kernel<PPISPParamType::RQS><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else if (param_type == "no_crf") {
        ppisp_backward_kernel<PPISPParamType::NoCRF><<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
            _dt3d_to_tv4<float>(in_image),
            (float*)std::get<0>(ppisp_params),
            (float4*)std::get<0>(intrins),
            actual_image_width,
            actual_image_height,
            cam_idx_ptr,
            _dt3d_to_tv4<float>(v_out_image),
            _dt3d_to_tv4<float>(v_in_image),
            (float*)std::get<0>(v_ppisp_params)
        );
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

template<PPISPParamType param_type>
__global__ void compute_raw_ppisp_regularization_forward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    float* __restrict__ raw_losses  // [B+1, RawPPISPRegLossIndex::length]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    bool inside = (bid < (unsigned)B);

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

    FixedArray<float, kNumRawLosses> losses;
    if (inside) {
        FixedArray<float, kNumParams> params;
        #pragma unroll
        for (int i = 0; i < kNumParams; i++) {
            params[i] = ppisp_params[bid * kNumParams + i];
        }

        if (param_type == PPISPParamType::Original)
            SlangPPISP::compute_raw_ppisp_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&losses)
            );
        else if (param_type == PPISPParamType::RQS)
            SlangPPISP::compute_raw_ppisp_rqs_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&losses)
            );
        else
            SlangPPISP::compute_raw_ppisp_no_crf_regularization_loss(
                *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
                reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&losses)
            );
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        float loss = (inside && isfinite(losses[i])) ? losses[i] : 0.0f;
        if (inside)
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
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

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
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_ppisp_rqs_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );
    else
        SlangPPISP::compute_ppisp_no_crf_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );

    #pragma unroll
    for (int i = 0; i < (int)PPISPRegLossIndex::length; i++) {
        losses_buffer[i] = losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
void compute_ppsip_regularization_forward(
    TorchTensorView ppisp_params,       // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    std::string param_type,
    TorchTensorView losses,             // [PPISPRegLossIndex::length] (must be pre-zeroed)
    TorchTensorView raw_losses          // [B+1, RawPPISPRegLossIndex::length] (must be pre-zeroed)
) {
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = std::get<2>(ppisp_params)[0];

    if (param_type == "original" || param_type == "") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndex::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "no_crf") {
        compute_raw_ppisp_regularization_forward_kernel<PPISPParamType::NoCRF>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            (float*)std::get<0>(raw_losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_ppisp_regularization_forward_kernel<PPISPParamType::NoCRF>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(int)RawPPISPRegLossIndexNoCRF::length,
            loss_weights,
            (float*)std::get<0>(losses)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
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

    static constexpr int kNumParams =
        (param_type == PPISPParamType::Original) ? kNumPPISPParams :
        (param_type == PPISPParamType::RQS)      ? kNumPPISPParamsRQS :
                                                   kNumPPISPParamsNoCRF;
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

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
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_raw_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParamsRQS>*>(&v_params)
        );
    else
        SlangPPISP::compute_raw_ppisp_no_crf_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParamsNoCRF>*>(&v_params)
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
    static constexpr int kNumRawLosses =
        (param_type == PPISPParamType::Original) ? (int)RawPPISPRegLossIndex::length :
        (param_type == PPISPParamType::RQS)      ? (int)RawPPISPRegLossIndexRQS::length :
                                                   (int)RawPPISPRegLossIndexNoCRF::length;

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
    else if (param_type == PPISPParamType::RQS)
        SlangPPISP::compute_ppisp_rqs_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexRQS::length>*>(&v_raw_losses)
        );
    else
        SlangPPISP::compute_ppisp_no_crf_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndexNoCRF::length>*>(&v_raw_losses)
        );

    #pragma unroll
    for (int i = 0; i < kNumRawLosses; i++) {
        v_raw_losses_buffer[i] = v_raw_losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
void compute_ppsip_regularization_backward(
    TorchTensorView ppisp_params,       // [B, PPISP_NUM_PARAMS]
    const std::array<float, (int)PPISPRegLossIndex::length> loss_weights_0,
    TorchTensorView raw_losses,         // [B+1, RawPPISPRegLossIndex::length]
    TorchTensorView v_losses,           // [PPISPRegLossIndex::length]
    std::string param_type,
    TorchTensorView v_ppisp_params      // [B, PPISP_NUM_PARAMS] (must be pre-zeroed)
) {
    FixedArray<float, (int)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (int)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = std::get<2>(ppisp_params)[0];

    if (param_type == "original" || param_type == "") {
        // v_raw_losses is a small scratch buffer
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndex::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndex::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndex::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::Original>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "rqs") {
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndexRQS::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndexRQS::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndexRQS::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::RQS>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else if (param_type == "no_crf") {
        float* v_raw_losses = DevicePool::global().acquire<float>(
            PoolSlot::PpispVRawLosses, (int)RawPPISPRegLossIndexNoCRF::length);
        cudaMemset(v_raw_losses, 0, (int)RawPPISPRegLossIndexNoCRF::length * sizeof(float));

        compute_ppisp_regularization_backward_kernel<PPISPParamType::NoCRF>
        <<<1, 1>>>(
            B,
            (float*)std::get<0>(raw_losses) + B*(uint)RawPPISPRegLossIndexNoCRF::length,
            loss_weights,
            (float*)std::get<0>(v_losses),
            v_raw_losses
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());

        compute_raw_ppisp_regularization_backward_kernel<PPISPParamType::NoCRF>
        <<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
            B,
            (float*)std::get<0>(ppisp_params),
            v_raw_losses,
            (float*)std::get<0>(v_ppisp_params)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    else {
        throw std::runtime_error("invalid PPISP param_type, must be \"original\", \"rqs\", or \"no_crf\"");
    }
}
