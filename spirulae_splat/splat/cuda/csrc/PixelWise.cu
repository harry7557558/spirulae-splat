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

    float3 rgb = in_rgb.load3(y, x);
    float alpha = in_alpha.load1(y, x);
    float3 background = in_background.load3(y, x);

    rgb = blend_background(rgb, alpha, background);

    out_rgb.store3(y, x, rgb);
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

    float3 rgb = in_rgb.load3(y, x);
    float alpha = in_alpha.load1(y, x);
    float3 background = in_background.load3(y, x);

    float3 v_out = v_out_rgb.load3(y, x);

    float3 v_rgb; float v_alpha; float3 v_background;
    blend_background_bwd(
        rgb, alpha, background,
        v_out,
        &v_rgb, &v_alpha, &v_background
    );

    v_in_rgb.store3(y, x, v_rgb);
    v_in_alpha.store1(y, x, v_alpha);
    v_in_background.store3(y, x, v_background);

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

    blend_background_forward_kernel<<<_LAUNCH_ARGS_1D(h*w, 256)>>>(
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

    blend_background_backward_kernel<<<_LAUNCH_ARGS_1D(h*w, 256)>>>(
        tensor2view<float, 3>(rgb), tensor2view<float, 3>(alpha), tensor2view<float, 3>(background),
        tensor2view<float, 3>(v_out_rgb),
        tensor2view<float, 3>(v_rgb), tensor2view<float, 3>(v_alpha), tensor2view<float, 3>(v_background)
    );

    return std::make_tuple(v_rgb, v_alpha, v_background);
}



// ================
// Depth to Normal
// ================


__global__ void depth_to_normal_forward_kernel(
    gsplat::CameraModelType camera_model,
    const float *__restrict__ Ks,  // [B, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    TensorView<float, 4> normals  // [B, H, W, 3]
) {
    const uint B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Zero for border pixels (consistent with Python)
    // TODO: possibly do a finite difference?
    if (i == 0 || i == W-1 || j == 0 || j == H-1) {
        normals.store3(bid, j, i, make_float3(0.0f));
        return;
    }

    // Load camera
    Ks += bid * 9;
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 normal;
    depth_to_normal(
        W, H, {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth,
        depth, &normal
    );
    normals.store3(bid, j, i, normal);

}


__global__ void depth_to_normal_backward_kernel(
    gsplat::CameraModelType camera_model,
    const float *__restrict__ Ks,  // [B, 3, 3]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const bool is_ray_depth,
    const TensorView<float, 4> depths,  // [B, H, W, 1]
    const TensorView<float, 4> v_normals,  // [B, H, W, 3]
    TensorView<float, 4> v_depths  // [B, H, W, 1]
) {
    const uint B = depths.shape[0],
        H = depths.shape[1],
        W = depths.shape[2];
    uint32_t bid = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (bid >= B || i >= W || j >= H)
        return;

    // Zero for border pixels (consistent with Python)
    // TODO: possibly do a finite difference?
    if (i == 0 || i == W-1 || j == 0 || j == H-1) {
        return;
    }

    // Load camera
    Ks += bid * 9;
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float4 depth = {
        depths.load1(bid, j, i-1),
        depths.load1(bid, j, i+1),
        depths.load1(bid, j-1, i),
        depths.load1(bid, j+1, i),
    };
    float3 v_normal = v_normals.load3(bid, j, i);
    float4 v_depth;
    depth_to_normal_vjp(
        W, H, {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE, is_ray_depth,
        depth, v_normal, &v_depth
    );
    v_depths.atomicStore1(bid, j, i-1, v_depth.x);
    v_depths.atomicStore1(bid, j, i+1, v_depth.y);
    v_depths.atomicStore1(bid, j-1, i, v_depth.z);
    v_depths.atomicStore1(bid, j+1, i, v_depth.w);
}


torch::Tensor depth_to_normal_forward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    torch::Tensor depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(Ks);

    if (Ks.ndimension() != 3 || Ks.size(-1) != 3 || Ks.size(-2) != 3)
        AT_ERROR("Ks shape must be (B, 3, 3)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != Ks.size(0))
        AT_ERROR("depths and Ks batch dimension mismatch");

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    torch::Tensor normals = torch::empty({b, h, w, 3}, depths.options());

    depth_to_normal_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths), tensor2view<float, 4>(normals)
    );

    return normals;
}

torch::Tensor depth_to_normal_backward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    bool is_ray_depth,
    torch::Tensor depths,  // [B, H, W, 1]
    torch::Tensor v_normals  // [B, H, W, 3]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(Ks);
    CHECK_CUDA(v_normals);

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    torch::Tensor v_depths = torch::zeros_like(depths);

    depth_to_normal_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        is_ray_depth, tensor2view<float, 4>(depths),
        tensor2view<float, 4>(v_normals),
        tensor2view<float, 4>(v_depths)
    );

    return v_depths;
}


// ================
// Ray Depth To Linear Depth
// ================


__global__ void ray_depth_to_linear_depth_forward_kernel(
    gsplat::CameraModelType camera_model,
    const float *__restrict__ Ks,  // [B, 3, 3]
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
    Ks += bid * 9;
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float in_depth = in_depths.load1(bid, j, i);
    float out_depth = in_depth * ray_depth_to_linear_depth_factor(
        W, H, {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE
    );
    out_depths.store1(bid, j, i, out_depth);
}

__global__ void ray_depth_to_linear_depth_backward_kernel(
    gsplat::CameraModelType camera_model,
    const float *__restrict__ Ks,  // [B, 3, 3]
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
    Ks += bid * 9;
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);

    // Process
    float v_out_depth = v_out_depths.load1(bid, j, i);
    float factor = ray_depth_to_linear_depth_factor(
        W, H, {(float)i+0.5f, (float)j+0.5f},
        {fx, fy, cx, cy}, &dist_coeffs,
        camera_model == gsplat::CameraModelType::FISHEYE
    );
    float v_in_depth = factor * v_out_depth;
    v_in_depths.store1(bid, j, i, v_in_depth);
}

torch::Tensor ray_depth_to_linear_depth_forward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(depths);
    CHECK_CUDA(depths);
    CHECK_INPUT(Ks);

    if (Ks.ndimension() != 3 || Ks.size(-1) != 3 || Ks.size(-2) != 3)
        AT_ERROR("Ks shape must be (B, 3, 3)");
    if (depths.ndimension() != 4 || depths.size(-1) != 1)
        AT_ERROR("depths shape must be (B, H, W, 1)");
    if (depths.size(0) != Ks.size(0))
        AT_ERROR("depths and Ks batch dimension mismatch");

    uint b = depths.size(0), h = depths.size(1), w = depths.size(2);
    torch::Tensor out_depths = torch::empty_like(depths);

    ray_depth_to_linear_depth_forward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(depths), tensor2view<float, 4>(out_depths)
    );

    return out_depths;
}

torch::Tensor ray_depth_to_linear_depth_backward_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor v_out_depths  // [B, H, W, 1]
) {
    DEVICE_GUARD(v_out_depths);
    CHECK_CUDA(v_out_depths);
    CHECK_INPUT(Ks);

    if (Ks.ndimension() != 3 || Ks.size(-1) != 3 || Ks.size(-2) != 3)
        AT_ERROR("Ks shape must be (B, 3, 3)");
    if (v_out_depths.ndimension() != 4 || v_out_depths.size(-1) != 1)
        AT_ERROR("v_out_depths shape must be (B, H, W, 1)");
    if (v_out_depths.size(0) != Ks.size(0))
        AT_ERROR("v_out_depths and Ks batch dimension mismatch");

    uint b = v_out_depths.size(0), h = v_out_depths.size(1), w = v_out_depths.size(2);
    torch::Tensor v_in_depths = torch::empty_like(v_out_depths);

    ray_depth_to_linear_depth_backward_kernel<<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(v_out_depths), tensor2view<float, 4>(v_in_depths)
    );

    return v_in_depths;
}


// ================
// Distort / Undistort
// ================

template<bool is_undistort>
__global__ void distort_image_kernel(
    gsplat::CameraModelType camera_model,
    const float *__restrict__ Ks,  // [B, 3, 3]
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
    Ks += bid * 9;
    float fx = Ks[0], fy = Ks[4], cx = Ks[2], cy = Ks[5];
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

torch::Tensor distort_image_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor in_image  // [B, H, W, C]
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(Ks);

    if (Ks.ndimension() != 3 || Ks.size(-1) != 3 || Ks.size(-2) != 3)
        AT_ERROR("Ks shape must be (B, 3, 3)");
    if (in_image.ndimension() != 4 || in_image.size(0) != Ks.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with Ks");

    uint b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    torch::Tensor out_image = torch::zeros_like(in_image);

    distort_image_kernel<false><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );

    return out_image;
}

torch::Tensor undistort_image_tensor(
    gsplat::CameraModelType camera_model,
    torch::Tensor Ks,  // [B, 3, 3]
    CameraDistortionCoeffsTensor dist_coeffs,
    torch::Tensor in_image  // [B, H, W, C]
) {
    DEVICE_GUARD(in_image);
    CHECK_CUDA(in_image);
    CHECK_INPUT(Ks);

    if (Ks.ndimension() != 3 || Ks.size(-1) != 3 || Ks.size(-2) != 3)
        AT_ERROR("Ks shape must be (B, 3, 3)");
    if (in_image.ndimension() != 4 || in_image.size(0) != Ks.size(0))
        AT_ERROR("in_image shape must be (B, H, W, C) where B is consistent with Ks");

    uint b = in_image.size(0), h = in_image.size(1), w = in_image.size(2);
    torch::Tensor out_image = torch::zeros_like(in_image);

    distort_image_kernel<true><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>(
        camera_model, Ks.data_ptr<float>(), dist_coeffs,
        tensor2view<float, 4>(in_image), tensor2view<float, 4>(out_image)
    );

    return out_image;
}
