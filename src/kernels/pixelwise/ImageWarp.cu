// ImageWarp.cu -- warps between wide (fisheye/equirectangular) and pinhole
// cameras, including the byte->float fused variants used by the DataManager.
//
// Part of the PixelWise family -- see PixelWiseCommon.cuh.

#include "kernels/pixelwise/BilinearSample.cuh"

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


