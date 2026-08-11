#pragma once

// Shared preamble for the image-space per-pixel kernels.
//
// PixelWise was one file (~3.7k lines) covering unrelated image-space
// operations. It is now split by function; the parts all declare into the
// single PixelWise.cuh via the explicit HEADER_SOURCES['PixelWise'] list in
// generate_headers.py. See docs/codegen.md.
//
//   ImageConvert.cu        uint8/uint16 -> float; rendered -> expected depth
//   ImageColorOps.cu       background blending, log map, overexposure reg
//   DepthGeometry.cu       depth -> points / normal, depth-normal loss,
//                          ray <-> linear depth
//   ImageDistort.cu        distort / undistort
//   ImageWarp.cu           wide <-> pinhole warps, incl. byte-fused
//   GtDepthNormalWarp.cu   GT depth / normal wide -> pinhole warps
//   Ppisp.cu               per-pixel image signal processing
//
// Device samplers shared across the distort/warp parts live in
// BilinearSample.cuh.

#include "kernels/pixelwise/PixelWise.cuh"

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
namespace SlangCameraSource {
#include "generated/set_namespace.cuh"
#include "generated/camera_source.cuh"
}
#endif

#include "core/Common.cuh"
#include "core/CameraDistortion.cuh"

#ifdef __CUDACC__

// Per-tier pixel_wise exports, same shape as SlangDistortion in
// core/CameraDistortion.cuh (slangc cannot export a generic).
template<CameraDistortionType D>
struct SlangPixelWiseDist;

#define _SS_DEF_SLANG_PIXEL_WISE(TIER, SUFFIX)                                     \
template<>                                                                         \
struct SlangPixelWiseDist<CameraDistortionType::TIER> {                            \
    using Coeffs = CameraDistortionCoeffsT<CameraDistortionType::TIER>;            \
                                                                                   \
    static __device__ __forceinline__                                              \
    float3 generate_ray_d2n(float2 pix_pos, float4 intrins, const Coeffs& c,        \
                            int camera_model, bool is_ray_depth) {                  \
        return SlangPixelWise::generate_ray_d2n##SUFFIX(                            \
            pix_pos, intrins, c.v, camera_model, is_ray_depth);                     \
    }                                                                              \
    static __device__ __forceinline__                                              \
    float3 depth_to_point(float2 pix_pos, float4 intrins, const Coeffs& c,          \
                          int camera_model, bool is_ray_depth, float depth) {       \
        return SlangPixelWise::depth_to_point##SUFFIX(                              \
            pix_pos, intrins, c.v, camera_model, is_ray_depth, depth);              \
    }                                                                              \
    static __device__ __forceinline__                                              \
    float depth_to_point_vjp(float2 pix_pos, float4 intrins, const Coeffs& c,       \
                             int camera_model, bool is_ray_depth,                   \
                             float depth, float3 v_point) {                         \
        return SlangPixelWise::depth_to_point_vjp##SUFFIX(                          \
            pix_pos, intrins, c.v, camera_model, is_ray_depth, depth, v_point);     \
    }                                                                              \
    static __device__ __forceinline__                                              \
    float3 depth_to_normal(float2 pix_center, float4 intrins, const Coeffs& c,      \
                           int camera_model, bool is_ray_depth, float4 depths) {    \
        return SlangPixelWise::depth_to_normal##SUFFIX(                             \
            pix_center, intrins, c.v, camera_model, is_ray_depth, depths);          \
    }                                                                              \
    static __device__ __forceinline__                                              \
    void depth_to_normal_vjp(float2 pix_center, float4 intrins, const Coeffs& c,    \
                             int camera_model, bool is_ray_depth, float4 depths,    \
                             float3 v_normal, float4* v_depths) {                   \
        SlangPixelWise::depth_to_normal_vjp##SUFFIX(                                \
            pix_center, intrins, c.v, camera_model, is_ray_depth, depths,           \
            v_normal, v_depths);                                                    \
    }                                                                              \
    static __device__ __forceinline__                                              \
    float ray_depth_to_linear_depth_factor(float2 pix_center, float4 intrins,       \
                                           const Coeffs& c, int camera_model) {     \
        return SlangPixelWise::ray_depth_to_linear_depth_factor##SUFFIX(            \
            pix_center, intrins, c.v, camera_model);                                \
    }                                                                              \
    static __device__ __forceinline__                                              \
    float depth_normal_loss(float2 pix_center, float4 intrins, const Coeffs& c,     \
                            int camera_model, bool is_ray_depth, float4 depths,     \
                            float3 gt_normal) {                                     \
        return SlangPixelWise::depth_normal_loss##SUFFIX(                           \
            pix_center, intrins, c.v, camera_model, is_ray_depth, depths,           \
            gt_normal);                                                             \
    }                                                                              \
    static __device__ __forceinline__                                              \
    void depth_normal_loss_vjp(float2 pix_center, float4 intrins, const Coeffs& c,  \
                               int camera_model, bool is_ray_depth, float4 depths,  \
                               float3 gt_normal, float v_loss,                      \
                               float4* v_depths, float3* v_gt_normal) {             \
        SlangPixelWise::depth_normal_loss_vjp##SUFFIX(                              \
            pix_center, intrins, c.v, camera_model, is_ray_depth, depths,           \
            gt_normal, v_loss, v_depths, v_gt_normal);                              \
    }                                                                              \
};

_SS_DEF_SLANG_PIXEL_WISE(None,      _none)
_SS_DEF_SLANG_PIXEL_WISE(OpenCV,    _opencv)
_SS_DEF_SLANG_PIXEL_WISE(ThinPrism, _prism)
_SS_DEF_SLANG_PIXEL_WISE(Rational,  _rational)

#undef _SS_DEF_SLANG_PIXEL_WISE

#endif  // __CUDACC__


// Run BODY(tier) for the tier `name` spells; the camera model stays runtime,
// so all four tiers are reachable for every model.
#define _SS_DISPATCH_DISTORTION(name, BODY)                                        \
    do { switch (cdt(name)) {                                                      \
        case CameraDistortionType::None:      BODY(CameraDistortionType::None);      break; \
        case CameraDistortionType::OpenCV:    BODY(CameraDistortionType::OpenCV);    break; \
        case CameraDistortionType::ThinPrism: BODY(CameraDistortionType::ThinPrism); break; \
        case CameraDistortionType::Rational:  BODY(CameraDistortionType::Rational);  break; \
        default: throw std::runtime_error(                                         \
            "Unknown camera distortion: " + std::string(name));                    \
    } } while (0)


static inline TensorView<float, 4> _bhw1_view(const TorchTensorView& tv) {
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
