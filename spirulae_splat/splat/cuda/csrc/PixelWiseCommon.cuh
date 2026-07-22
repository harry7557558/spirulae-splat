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
