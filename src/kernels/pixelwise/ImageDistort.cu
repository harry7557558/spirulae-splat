// ImageDistort.cu -- image distort / undistort.
//
// Part of the PixelWise family -- see PixelWiseCommon.cuh.

#include "kernels/pixelwise/BilinearSample.cuh"

// ================
// Distort / Undistort
// ================


template<CameraDistortionType distortion, bool is_undistort>
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
    using Dist = SlangDistortion<distortion>;
    typename Dist::Coeffs dist_coeffs = dist_coeffs_buffer.load<distortion>(bid);

    // Undistort point
    float2 uv = { (i+0.5f-cx) / fx, (j+0.5f-cy) / fy };
    if (is_undistort) {
        if (dot(uv, uv) > 0.0f && !Dist::is_valid_distortion(
            camera_model == CameraModelType::FISHEYE ?
                normalize(uv) * atanf(length(uv)) :
            camera_model == CameraModelType::EQUISOLID ?
                normalize(uv) * (2.0f * sinf(0.5f * atanf(length(uv)))) :
            uv,
            dist_coeffs
        ))
            return;
        uv = Dist::distort_point(uv, (int)camera_model, dist_coeffs);
    }
    else {
        if (!Dist::undistort_point(uv, (int)camera_model, dist_coeffs, &uv))
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
    std::string distortion,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 8]
    TorchTensorView in_image,           // [B, H, W, C] float
    TorchTensorView out_image           // [B, H, W, C] float (must be pre-zeroed)
) {
    const auto& s = std::get<2>(in_image);
    int b = s[0], h = s[1], w = s[2];

    #define LAUNCH(D) \
        distort_image_kernel<D, false><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>( \
            cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs, \
            _bhwc_view(in_image), _bhwc_view(out_image))
    _SS_DISPATCH_DISTORTION(distortion, LAUNCH);
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void undistort_image_tensor(
    std::string camera_model,
    std::string distortion,
    TorchTensorView intrins,            // [B, 4]
    TorchTensorView dist_coeffs,        // [B, 8]
    TorchTensorView in_image,           // [B, H, W, C] float
    TorchTensorView out_image           // [B, H, W, C] float (must be pre-zeroed)
) {
    const auto& s = std::get<2>(in_image);
    int b = s[0], h = s[1], w = s[2];

    #define LAUNCH(D) \
        distort_image_kernel<D, true><<<_LAUNCH_ARGS_3D(w, h, b, 16, 16, 1)>>>( \
            cmt(camera_model), (float4*)std::get<0>(intrins), dist_coeffs, \
            _bhwc_view(in_image), _bhwc_view(out_image))
    _SS_DISPATCH_DISTORTION(distortion, LAUNCH);
    #undef LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


