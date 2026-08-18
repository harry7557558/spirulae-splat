#pragma once

// Device-side bilinear samplers shared by the distort/warp kernels
// (ImageDistort.cu, ImageWarp.cu, GtDepthNormalWarp.cu). The `_valid` pair
// skips GT sentinel taps; the byte variants fuse the normalizing conversion.
//
// (x, y) follow this codebase's pixel convention: the centre of pixel i is at
// i + 0.5 (`px = j + 0.5f` in the rasterizer). Every sampler subtracts that
// half pixel, so a pixel centre reads one tap and an identity remap is exact.

#include "kernels/pixelwise/PixelWiseCommon.cuh"

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
    x -= 0.5f; y -= 0.5f;
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
    x -= 0.5f; y -= 0.5f;
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
    x -= 0.5f; y -= 0.5f;
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
    float x, float y, float norm_inv
) {
    const long W = image.shape[2], H = image.shape[1];
    x -= 0.5f; y -= 0.5f;
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
    // y clamps: the half-pixel offset puts the pole half-row outside the grid,
    // and padding it would band the top and bottom of every panorama.
    auto fetch = [&](long y_, long xw) -> float {
        y_ = y_ < 0 ? 0 : (y_ >= H ? H - 1 : y_);
        return (float)image.at(bid, y_, xw, cid) * norm_inv;
    };
    return fetch(y0, x0w) * (wx0 * wy0)
         + fetch(y0, x1w) * (wx1 * wy0)
         + fetch(y1, x0w) * (wx0 * wy1)
         + fetch(y1, x1w) * (wx1 * wy1);
}


// The two samplers below drop a GT sentinel tap instead of blending it, and
// report "nothing here" when no tap has data. The invariant is stated once at
// revalidate_weights in core/Interpolation.cuh.

namespace _bilinear_valid {

// The four taps at (x, y) plus their weights, indexed q = 2*dy + dx.
// `wrap_u` wraps the x column (equirectangular); y is left to the bounds test.
struct Taps {
    long  x[2], y[2];
    float w[4];
    bool  x_ok[2];
};

inline __device__ Taps resolve(float x, float y, long W, bool wrap_u) {
    Taps t;
    x -= 0.5f; y -= 0.5f;
    long x0 = (long)floorf(x), y0 = (long)floorf(y);
    float fx = x - (float)x0, fy = y - (float)y0;
    t.y[0] = y0; t.y[1] = y0 + 1;
    t.w[0] = (1.0f - fx) * (1.0f - fy);
    t.w[1] = fx * (1.0f - fy);
    t.w[2] = (1.0f - fx) * fy;
    t.w[3] = fx * fy;
    for (int k = 0; k < 2; k++) {
        long xv = x0 + k;
        if (wrap_u) {
            long r = xv % W;
            t.x[k] = r < 0 ? r + W : r;
            t.x_ok[k] = true;
        } else {
            t.x[k] = xv;
            t.x_ok[k] = (xv >= 0 && xv < W);
        }
    }
    return t;
}

}  // namespace _bilinear_valid

// False (and *out = 0) when no tap carries depth.
template<typename T_in>
inline __device__ bool bilinear_depth_valid(
    const TensorView<T_in, 4>& image, uint32_t bid,
    float x, float y, float norm_inv, bool wrap_u, float* out
) {
    const long W = image.shape[2], H = image.shape[1];
    _bilinear_valid::Taps t = _bilinear_valid::resolve(x, y, W, wrap_u);
    float acc = 0.0f, wsum = 0.0f;
    for (int q = 0; q < 4; q++) {
        long yv = t.y[q >> 1];
        int  xi = q & 1;
        if (!t.x_ok[xi] || yv < 0 || yv >= H) continue;
        float d = (float)image.at(bid, yv, t.x[xi], 0) * norm_inv;
        if (d <= 0.0f) continue;
        acc  += t.w[q] * d;
        wsum += t.w[q];
    }
    *out = wsum > 0.0f ? acc / wsum : 0.0f;
    return wsum > 0.0f;
}

// The sampled normal is NOT renormalized here -- the callers rotate it into
// their own frame first and normalize after. False (and *out = the
// (-1,-1,-1) sentinel) when no tap carries a normal.
template<typename T_in>
inline __device__ bool bilinear_normal_valid(
    const TensorView<T_in, 4>& image, uint32_t bid,
    float x, float y, float norm_inv, float decode_off, bool wrap_u,
    float3* out
) {
    const long W = image.shape[2], H = image.shape[1];
    _bilinear_valid::Taps t = _bilinear_valid::resolve(x, y, W, wrap_u);
    float3 acc = make_float3(0.0f, 0.0f, 0.0f);
    float wsum = 0.0f;
    for (int q = 0; q < 4; q++) {
        long yv = t.y[q >> 1];
        int  xi = q & 1;
        if (!t.x_ok[xi] || yv < 0 || yv >= H) continue;
        float3 n = make_float3(
            (float)image.at(bid, yv, t.x[xi], 0) * norm_inv + decode_off,
            (float)image.at(bid, yv, t.x[xi], 1) * norm_inv + decode_off,
            (float)image.at(bid, yv, t.x[xi], 2) * norm_inv + decode_off);
        // GT normal validity: core/Interpolation.cuh's gt_normal_valid.
        if (n.x + n.y + n.z <= -2.366f || dot(n, n) <= 0.25f) continue;
        acc  += t.w[q] * n;
        wsum += t.w[q];
    }
    if (wsum <= 0.0f) {
        *out = make_float3(-1.0f, -1.0f, -1.0f);
        return false;
    }
    *out = acc / wsum;
    return true;
}
