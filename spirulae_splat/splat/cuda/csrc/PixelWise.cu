#include "PixelWise.cuh"

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
// Blend Background
// ================

__global__ void blend_background_forward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_alpha,
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
    float alpha = in_alpha.load1(bid, y, x);
    float3 background = in_background.load3(bid, y, x);

    rgb = blend_background(rgb, alpha, background);

    out_rgb.store3(bid, y, x, rgb);
}


__global__ void blend_background_backward_kernel(
    const TensorView<float, 4> in_rgb,
    const TensorView<float, 4> in_alpha,
    const TensorView<float, 4> in_background,
    const TensorView<float, 4> v_out_rgb,
    TensorView<float, 4> v_in_rgb,
    TensorView<float, 4> v_in_alpha,
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
    float alpha = in_alpha.load1(bid, y, x);
    float3 background = in_background.load3(bid, y, x);

    float3 v_out = v_out_rgb.load3(bid, y, x);

    float3 v_rgb; float v_alpha; float3 v_background;
    blend_background_bwd(
        rgb, alpha, background,
        v_out,
        &v_rgb, &v_alpha, &v_background
    );

    v_in_rgb.store3(bid, y, x, v_rgb);
    v_in_alpha.store1(bid, y, x, v_alpha);
    v_in_background.store3(bid, y, x, v_background);

}



/*[AutoHeaderGeneratorExport]*/
at::Tensor blend_background_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &alpha,  // [B, H, W, 1]
    at::Tensor &background  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(alpha);
    CHECK_CUDA(background);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);
    if (alpha.ndimension() != 4 || alpha.size(0) != b || alpha.size(1) != h || alpha.size(2) != w || alpha.size(3) != 1)
        AT_ERROR("alpha shape must be (b, h, w, 1)");
    if (background.ndimension() != 4 || background.size(0) != b || background.size(1) != h || background.size(2) != w || background.size(3) != 3)
        AT_ERROR("background shape must be (b, h, w, 3)");

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    blend_background_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(alpha), tensor2view<float, 4>(background),
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor, at::Tensor>
blend_background_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &alpha,  // [B, H, W, 1]
    at::Tensor &background,  // [B, H, W, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(alpha);
    CHECK_CUDA(background);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());
    at::Tensor v_alpha = at::empty({b, h, w, 1}, alpha.options());
    at::Tensor v_background = at::empty({b, h, w, 3}, background.options());

    blend_background_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), tensor2view<float, 4>(alpha), tensor2view<float, 4>(background),
        tensor2view<float, 4>(v_out_rgb),
        tensor2view<float, 4>(v_rgb), tensor2view<float, 4>(v_alpha), tensor2view<float, 4>(v_background)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(v_rgb, v_alpha, v_background);
}


// ================
// Log Map Image
// ================

__global__ void log_map_image_forward_kernel(
    const TensorView<float, 4> in_rgb,
    float t,
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

    rgb = log_map_image(rgb, t);

    out_rgb.store3(bid, y, x, rgb);
}


__global__ void log_map_image_backward_kernel(
    const TensorView<float, 4> in_rgb,
    float t,
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

    float3 v_out = v_out_rgb.load3(bid, y, x);

    float3 v_rgb = log_map_image_bwd(rgb, t, v_out);

    v_in_rgb.store3(bid, y, x, v_rgb);
}



/*[AutoHeaderGeneratorExport]*/
at::Tensor log_map_image_forward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    float t
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    log_map_image_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), t,
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor log_map_image_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    float t,
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());

    log_map_image_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb), t,
        tensor2view<float, 4>(v_out_rgb),
        tensor2view<float, 4>(v_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_rgb;
}



// ================
// Depth to Normal
// ================


__global__ void depth_to_normal_forward_kernel(
    gsplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    TensorView<float, 4> normals  // [B, H, W, 3]
) {
    const uint B = depths.shape[0],
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
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth,
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
        float3 ray = generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, &dist_coeffs,
            camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth
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
    float3 normal = points_to_normal(&points);
#endif
    normals.store3(bid, j, i, normal);

}


__global__ void depth_to_normal_backward_kernel(
    gsplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    const TensorView<float, 4> v_normals,  // [B, H, W, 3]
    TensorView<float, 4> v_depths  // [B, H, W, 1]
) {
    const uint B = depths.shape[0],
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
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth,
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
        float3 ray = generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, &dist_coeffs,
            camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth
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
    points_to_normal_vjp(&points, v_normal, &v_points);

    v_depths.atomicStore1(bid, j, i-1, dot(v_points[0], rays[0]));
    v_depths.atomicStore1(bid, j, i+1, dot(v_points[1], rays[1]));
    v_depths.atomicStore1(bid, j-1, i, dot(v_points[2], rays[2]));
    v_depths.atomicStore1(bid, j+1, i, dot(v_points[3], rays[3]));
#endif
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_normal_forward_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor normals = at::empty({b, h, w, 3}, depths.options());

    depth_to_normal_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths), tensor2view<float, 4>(normals)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return normals;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor depth_to_normal_backward_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor v_depths = at::zeros_like(depths);

    depth_to_normal_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths),
        tensor2view<float, 4>(v_normals),
        tensor2view<float, 4>(v_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_depths;
}


// ================
// Ray Depth To Linear Depth
// ================


__global__ void ray_depth_to_linear_depth_forward_kernel(
    gsplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> in_depths,  // [B, H, W, 1]
    TensorView<float, 4> out_depths  // [B, H, W, 1]
) {
    const uint B = in_depths.shape[0],
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
    float out_depth = in_depth * ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE
    );
    out_depths.store1(bid, j, i, out_depth);
}

__global__ void ray_depth_to_linear_depth_backward_kernel(
    gsplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> v_out_depths,  // [B, H, W, 1]
    TensorView<float, 4> v_in_depths  // [B, H, W, 1]
) {
    const uint B = v_out_depths.shape[0],
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
    float factor = ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE
    );
    float v_in_depth = factor * v_out_depth;
    v_in_depths.store1(bid, j, i, v_in_depth);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ray_depth_to_linear_depth_forward_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    at::Tensor out_depths = at::empty_like(depths);

    ray_depth_to_linear_depth_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(depths), tensor2view<float, 4>(out_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_depths;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ray_depth_to_linear_depth_backward_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = v_out_depths.size(0), h = v_out_depths.size(1), w = v_out_depths.size(2);
    at::Tensor v_in_depths = at::empty_like(v_out_depths);

    ray_depth_to_linear_depth_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(v_out_depths), tensor2view<float, 4>(v_in_depths)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_in_depths;
}


// ================
// Distort / Undistort
// ================

template<bool is_undistort>
__global__ void distort_image_kernel(
    gsplat::CameraModelType camera_model,
    const float4 *__restrict__ intrins,  // [B, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const TensorView<float, 4> in_image,  // [B, H, W, C]
    TensorView<float, 4> out_image  // [B, H, W, C]
) {
    const uint B = in_image.shape[0],
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
    if (is_undistort)
        uv = distort_point(uv, camera_model == gsplat::CameraModelType::FISHEYE, &dist_coeffs);
    else {
        if (!undistort_point(uv, camera_model == gsplat::CameraModelType::FISHEYE, &dist_coeffs, &uv))
            return;
    }
    int i1 = (int)floor(uv.x*fx+cx);
    int j1 = (int)floor(uv.y*fy+cy);

    // Process
    if (i1 < 0 || i1 >= W || j1 < 0 || j1 >= H)
        return;
    for (int c = 0; c < C; c++) {
        // TODO: bilinear/bicubic interpolation
        out_image.at(bid, j, i, c) = in_image.at(bid, j1, i1, c);
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor distort_image_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = at::zeros_like(in_image);

    distort_image_kernel<false><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_image;
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor undistort_image_tensor(
    gsplat::CameraModelType camera_model,
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

    uint b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = at::zeros_like(in_image);

    distort_image_kernel<true><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, (float4*)intrins.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_image;
}


// ================
// PPISP
// ================

static constexpr int kNumPPISPParams = 36;

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

    FixedArray<float, kNumPPISPParams> params;
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        params[i] = ppisp_params[bid * kNumPPISPParams + i];
    }

    float3 pixel = in_image.load3(bid, y, x);

    float3 out_pixel = apply_ppisp(
        pixel,
        (float2){(float)x, (float)y},
        (float2){intrins[bid].z, intrins[bid].w},
        (float2){actual_image_width, actual_image_height},
        &params
    );

    out_image.store3(bid, y, x, out_pixel);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor ppisp_forward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(intrins);
    if (in_image.ndimension() != 4 || in_image.size(0) != ppisp_params.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with ppisp_params");
    if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParams)
        AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS)");
    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    long b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor out_image = at::empty_like(in_image);
    ppisp_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(in_image),
        ppisp_params.data_ptr<float>(),
        (float4*)intrins.data_ptr<float>(),
        actual_image_width,
        actual_image_height,
        tensor2view<float, 4>(out_image)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    return out_image;
}

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

#if 0
    FixedArray<float, kNumPPISPParams> params;
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        params[i] = ppisp_params[bid * kNumPPISPParams + i];
    }
#else
    __shared__ float params_shared[kNumPPISPParams];
    if (threadIdx.x < kNumPPISPParams) {  // assume blockDim.x >= kNumPPISPParams
        float value = ppisp_params[bid * kNumPPISPParams + threadIdx.x];
        params_shared[threadIdx.x] = value;
    }
    __syncthreads();
    FixedArray<float, kNumPPISPParams> params;
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        params[i] = params_shared[i];
        // int j = (i + threadIdx.x) % kNumPPISPParams;
        // params[j] = params_shared[j];
    }
#endif

    float3 v_pixel;
    FixedArray<float, kNumPPISPParams> v_params;
    apply_ppisp_vjp(
        pixel,
        (float2){(float)x, (float)y},
        (float2){intrins[bid].z, intrins[bid].w},
        (float2){actual_image_width, actual_image_height},
        &params,
        v_out_pixel,
        &v_pixel,
        &v_params
    );

    v_in_image.store3(bid, y, x, v_pixel);

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        param = cg::reduce(warp, param, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && param != 0.0f)
            atomicAdd(&v_ppisp_params[bid * kNumPPISPParams + i], param);
    }
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor> ppisp_backward_tensor(
    at::Tensor &in_image,  // [B, H, W, C]
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    at::Tensor &intrins,  // [B, 4]
    const float actual_image_width,
    const float actual_image_height,
    at::Tensor &v_out_image  // [B, H, W, C]
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(intrins);
    CHECK_CUDA(v_out_image);
    if (in_image.ndimension() != 4 || in_image.size(0) != ppisp_params.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with ppisp_params");
    if (ppisp_params.ndimension() != 2 || ppisp_params.size(1) != kNumPPISPParams)
        AT_ERROR("ppisp_params shape must be (B, PPISP_NUM_PARAMS)");
    if (intrins.ndimension() != 2 || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (B, 4)");
    long b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    at::Tensor v_in_image = at::empty_like(in_image);
    at::Tensor v_ppisp_params = at::zeros_like(ppisp_params);
    ppisp_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 64, 1)>>>(
        tensor2view<float, 4>(in_image),
        ppisp_params.data_ptr<float>(),
        (float4*)intrins.data_ptr<float>(),
        actual_image_width,
        actual_image_height,
        tensor2view<float, 4>(v_out_image),
        tensor2view<float, 4>(v_in_image),
        v_ppisp_params.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    return std::make_tuple(v_in_image, v_ppisp_params);
}

__global__ void compute_raw_ppisp_regularization_forward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    float* __restrict__ raw_losses  // [B+1, RawPPISPRegLossIndex::length]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= B)
        return;

    FixedArray<float, kNumPPISPParams> params;
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        params[i] = ppisp_params[bid * kNumPPISPParams + i];
    }

    FixedArray<float, (uint)RawPPISPRegLossIndex::length> losses;
    compute_raw_ppisp_regularization_loss(&params, &losses);

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    #pragma unroll
    for (int i = 0; i < (uint)RawPPISPRegLossIndex::length; i++) {
        float loss = isfinite(losses[i]) ? losses[i] : 0.0f;
        raw_losses[bid*(uint)RawPPISPRegLossIndex::length + i] = loss;
        loss = cg::reduce(warp, loss, cg::plus<float>());
        if (threadIdx.x % WARP_SIZE == 0 && loss != 0.0f)
            atomicAdd(&raw_losses[B*(uint)RawPPISPRegLossIndex::length + i], loss);
    }
}

__global__ void compute_ppisp_regularization_forward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (uint)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    float* __restrict__ losses_buffer  // [PPISPRegLossIndex::length]
) {
    FixedArray<float, (uint)RawPPISPRegLossIndex::length> raw_losses;
    #pragma unroll
    for (int i = 0; i < (uint)RawPPISPRegLossIndex::length; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (uint)PPISPRegLossIndex::length> losses;

    compute_ppisp_regularization_loss(
        &raw_losses, num_train_images, &loss_weights, &losses
    );

    #pragma unroll
    for (int i = 0; i < (uint)PPISPRegLossIndex::length; i++) {
        losses_buffer[i] = losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
compute_ppsip_regularization_forward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (uint)PPISPRegLossIndex::length> loss_weights_0
) {
    DEVICE_GUARD(ppisp_params);
    CHECK_INPUT(ppisp_params);

    FixedArray<float, (uint)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = ppisp_params.size(0);

    at::Tensor raw_losses = at::zeros({B+1, (uint)RawPPISPRegLossIndex::length}, ppisp_params.options());
    at::Tensor losses = at::zeros({(uint)PPISPRegLossIndex::length}, ppisp_params.options());

    compute_raw_ppisp_regularization_forward_kernel<<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
        B,
        ppisp_params.data_ptr<float>(),
        raw_losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    compute_ppisp_regularization_forward_kernel<<<1, 1>>>(
        B,
        raw_losses.data_ptr<float>() + B*(uint)RawPPISPRegLossIndex::length,
        loss_weights,
        losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(losses, raw_losses);
}

__global__ void compute_raw_ppisp_regularization_backward_kernel(
    int B,  // number of images
    const float* __restrict__ ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const float* __restrict__ v_raw_losses,  // [RawPPISPRegLossIndex::length]
    float* __restrict__ v_ppisp_params  // [B, PPISP_NUM_PARAMS]
) {
    unsigned bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= B)
        return;

    FixedArray<float, kNumPPISPParams> params;
    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        params[i] = ppisp_params[bid * kNumPPISPParams + i];
    }

    FixedArray<float, (uint)RawPPISPRegLossIndex::length> v_losses;
    #pragma unroll
    for (int i = 0; i < (uint)RawPPISPRegLossIndex::length; i++) {
        v_losses[i] = v_raw_losses[i];
    }

    FixedArray<float, kNumPPISPParams> v_params;
    compute_raw_ppisp_regularization_loss_vjp(&params, &v_losses, &v_params);

    #pragma unroll
    for (int i = 0; i < kNumPPISPParams; i++) {
        float param = isfinite(v_params[i]) ? v_params[i] : 0.0f;
        v_ppisp_params[bid * kNumPPISPParams + i] = param;
    }
}

__global__ void compute_ppisp_regularization_backward_kernel(
    int num_train_images,
    const float* __restrict__ raw_losses_buffer,  // [RawPPISPRegLossIndex::length]
    FixedArray<float, (uint)PPISPRegLossIndex::length> loss_weights,  // [PPISPRegLossIndex::length]
    const float* __restrict__ v_losses_buffer,  // [PPISPRegLossIndex::length]
    float* __restrict__ v_raw_losses_buffer  // [RawPPISPRegLossIndex::length]
) {
    FixedArray<float, (uint)RawPPISPRegLossIndex::length> raw_losses;
    #pragma unroll
    for (int i = 0; i < (uint)RawPPISPRegLossIndex::length; i++) {
        raw_losses[i] = raw_losses_buffer[i];
    }

    FixedArray<float, (uint)PPISPRegLossIndex::length> v_losses;
    #pragma unroll
    for (int i = 0; i < (uint)PPISPRegLossIndex::length; i++) {
        v_losses[i] = v_losses_buffer[i];
    }

    FixedArray<float, (uint)RawPPISPRegLossIndex::length> v_raw_losses;
    compute_ppisp_regularization_loss_vjp(
        &raw_losses, num_train_images, &loss_weights, &v_losses, &v_raw_losses
    );

    #pragma unroll
    for (int i = 0; i < (uint)RawPPISPRegLossIndex::length; i++) {
        v_raw_losses_buffer[i] = v_raw_losses[i];
    }
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor compute_ppsip_regularization_backward_tensor(
    at::Tensor &ppisp_params,  // [B, PPISP_NUM_PARAMS]
    const std::array<float, (uint)PPISPRegLossIndex::length> loss_weights_0,
    at::Tensor &raw_losses,  // [B+1, RawPPISPRegLossIndex::length]
    at::Tensor &v_losses  // [PPISPRegLossIndex::length]
) {
    DEVICE_GUARD(ppisp_params);
    CHECK_INPUT(ppisp_params);
    CHECK_INPUT(raw_losses);
    CHECK_INPUT(v_losses);

    FixedArray<float, (uint)PPISPRegLossIndex::length> loss_weights =
        *reinterpret_cast<const FixedArray<float, (uint)PPISPRegLossIndex::length>*>(loss_weights_0.data());

    long B = ppisp_params.size(0);

    at::Tensor v_raw_losses = at::zeros({(uint)RawPPISPRegLossIndex::length}, ppisp_params.options());
    at::Tensor v_ppisp_params = at::zeros_like(ppisp_params);

    compute_ppisp_regularization_backward_kernel<<<1, 1>>>(
        B,
        raw_losses.data_ptr<float>() + B*(uint)RawPPISPRegLossIndex::length,
        loss_weights,
        v_losses.data_ptr<float>(),
        v_raw_losses.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    compute_raw_ppisp_regularization_backward_kernel<<<_LAUNCH_ARGS_1D(B, WARP_SIZE)>>>(
        B,
        ppisp_params.data_ptr<float>(),
        v_raw_losses.data_ptr<float>(),
        v_ppisp_params.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return v_ppisp_params;
}
