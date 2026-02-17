#include "PixelWise.cuh"

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangAll {
#include "generated/set_namespace.cuh"
#include "generated/slang_all.cuh"
}
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
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

    rgb = SlangAll::blend_background(rgb, alpha, background);

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
    SlangAll::blend_background_bwd(
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

__global__ void linear_rgb_to_srgb_forward_kernel(
    const TensorView<float, 4> in_rgb,
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

    rgb = SlangAll::linear_rgb_to_srgb(rgb);

    out_rgb.store3(bid, y, x, rgb);
}


__global__ void linear_rgb_to_srgb_backward_kernel(
    const TensorView<float, 4> in_rgb,
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

    float3 v_rgb = SlangAll::linear_rgb_to_srgb_bwd(rgb, v_out);

    v_in_rgb.store3(bid, y, x, v_rgb);
}



/*[AutoHeaderGeneratorExport]*/
at::Tensor linear_rgb_to_srgb_forward_tensor(
    at::Tensor &rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);

    if (rgb.ndimension() != 4 || rgb.size(-1) != 3)
        AT_ERROR("rgb shape must be (b, h, w, 3)");
    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor out_rgb = at::empty({b, h, w, 3}, rgb.options());

    linear_rgb_to_srgb_forward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb),
        tensor2view<float, 4>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor linear_rgb_to_srgb_backward_tensor(
    at::Tensor &rgb,  // [B, H, W, 3]
    at::Tensor &v_out_rgb  // [B, H, W, 3]
) {
    DEVICE_GUARD(rgb);
    CHECK_CUDA(rgb);
    CHECK_CUDA(v_out_rgb);

    long b = rgb.size(0), h = rgb.size(1), w = rgb.size(2);

    at::Tensor v_rgb = at::empty({b, h, w, 3}, rgb.options());

    linear_rgb_to_srgb_backward_kernel<<<_LAUNCH_ARGS_2D(h*w, b, 256, 1)>>>(
        tensor2view<float, 4>(rgb),
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
        float3 ray = SlangAll::generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, dist_coeffs,
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
    float3 normal = SlangAll::points_to_normal(points);
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
        float3 ray = SlangAll::generate_ray_d2n(
            {(float)ig+0.5f, (float)jg+0.5f},
            {fx, fy, cx, cy}, dist_coeffs,
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
    SlangAll::points_to_normal_vjp(points, v_normal, &v_points);

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

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
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

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
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
    float out_depth = in_depth * SlangAll::ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
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
    float factor = SlangAll::ray_depth_to_linear_depth_factor(
        {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, dist_coeffs,
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

    int b = depths.size(0), h = depths.size(1), w = depths.size(2);
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

    int b = v_out_depths.size(0), h = v_out_depths.size(1), w = v_out_depths.size(2);
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

__device__ float get_pixel_bilinear(
    const TensorView<float, 4> image,  // [B, H, W, C]
    uint32_t bid,
    uint32_t cid,
    float x,
    float y
) {
    const int W = image.shape[2],
        H = image.shape[1];
    int x0 = (int)floorf(x);
    int x1 = x0 + 1;
    int y0 = (int)floorf(y);
    int y1 = y0 + 1;
    float wx1 = x - x0;
    float wx0 = 1.0f - wx1;
    float wy1 = y - y0;
    float wy0 = 1.0f - wy1;

    float c00 = (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, y0, x0, cid) : 0.0f;
    float c10 = (x1 >= 0 && x1 < W && y0 >= 0 && y0 < H) ?
        image.at(bid, y0, x1, cid) : 0.0f;
    float c01 = (x0 >= 0 && x0 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, y1, x0, cid) : 0.0f;
    float c11 = (x1 >= 0 && x1 < W && y1 >= 0 && y1 < H) ?
        image.at(bid, y1, x1, cid) : 0.0f;

    float c = 0.0f;
    c += c00 * (wx0 * wy0);
    c += c10 * (wx1 * wy0);
    c += c01 * (wx0 * wy1);
    c += c11 * (wx1 * wy1);
    return c;
}

template<bool is_undistort>
__global__ void distort_image_kernel(
    gsplat::CameraModelType camera_model,
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
        if (!SlangProjectionUtils::is_valid_distortion(uv, dist_coeffs))
            return;
        uv = SlangProjectionUtils::distort_point(uv, camera_model == gsplat::CameraModelType::FISHEYE, dist_coeffs);
    }
    else {
        if (!SlangProjectionUtils::undistort_point(uv, camera_model == gsplat::CameraModelType::FISHEYE, dist_coeffs, &uv))
            return;
    }

    // Process
    for (int c = 0; c < C; c++) {
        out_image.at(bid, j, i, c) = get_pixel_bilinear(in_image, bid, c, uv.x*fx+cx, uv.y*fy+cy);
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

    int b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
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

    int b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
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
        out_pixel = SlangAll::apply_ppisp(
            pixel,
            make_float2((float)x, (float)y),
            make_float2(intrins[bid].z, intrins[bid].w),
            make_float2(actual_image_width, actual_image_height),
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params)
        );
    else
        out_pixel = SlangAll::apply_ppisp_rqs(
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
        SlangAll::apply_ppisp_vjp(
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
        SlangAll::apply_ppisp_rqs_vjp(
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
    at::Tensor v_ppisp_params = at::zeros_like(ppisp_params);
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
        SlangAll::compute_raw_ppisp_regularization_loss(
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&losses)
        );
    else
        SlangAll::compute_raw_ppisp_rqs_regularization_loss(
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
        SlangAll::compute_ppisp_regularization_loss(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, &losses
        );
    else
        SlangAll::compute_ppisp_rqs_regularization_loss(
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
        SlangAll::compute_raw_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&params),
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_losses),
            reinterpret_cast<FixedArray<float, kNumPPISPParams>*>(&v_params)
        );
    else
        SlangAll::compute_raw_ppisp_rqs_regularization_loss_vjp(
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
        SlangAll::compute_ppisp_regularization_loss_vjp(
            *reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&raw_losses),
            num_train_images, loss_weights, v_losses,
            reinterpret_cast<FixedArray<float, (int)RawPPISPRegLossIndex::length>*>(&v_raw_losses)
        );
    else
        SlangAll::compute_ppisp_rqs_regularization_loss_vjp(
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
    at::Tensor v_ppisp_params = at::zeros_like(ppisp_params);

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
