#pragma once

// Device-side bilinear samplers shared by the distort/warp kernels
// (ImageDistort.cu, ImageWarp.cu, GtDepthNormalWarp.cu).
//
//   get_pixel_bilinear         float sampling of [B,H,W,C] / [B,K,H,W,C]
//   bilinear_byte_norm         uint8/uint16 sampling fused with the
//                              byte -> normalized-float conversion
//   bilinear_byte_norm_wrap_u  ditto, wrapping x for equirectangular input

#include "PixelWiseCommon.cuh"

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
