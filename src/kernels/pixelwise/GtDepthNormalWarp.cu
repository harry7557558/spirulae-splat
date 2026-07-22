// GtDepthNormalWarp.cu -- ground-truth depth / normal wide -> pinhole warps
// for split-mode supervision.
//
// Part of the PixelWise family -- see PixelWiseCommon.cuh.

#include "kernels/pixelwise/BilinearSample.cuh"

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


