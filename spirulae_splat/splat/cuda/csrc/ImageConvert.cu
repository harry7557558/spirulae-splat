// ImageConvert.cu -- uint8/uint16 image, normal and depth buffers to float;
// rendered depth -> expected depth.
//
// Part of the PixelWise family -- see PixelWiseCommon.cuh.

#include "PixelWiseCommon.cuh"


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



