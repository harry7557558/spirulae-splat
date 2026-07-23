#pragma once

// Bilinear sampling helpers used by the per-pixel fused loss to read GT
// modalities (depth / normal / mask) whose resolution differs from the
// rendered output. Header-only, __device__-callable.
//
// Coordinate convention: integer pixel (x_dst, y_dst) at destination
// resolution (W_dst, H_dst) maps to source resolution (W_src, H_src) using
// the "align-corners=false" / "half-pixel-center" convention that matches
// torch.nn.functional.grid_sample / OpenCV warps:
//
//     u = (x_dst + 0.5) * W_src / W_dst - 0.5
//     v = (y_dst + 0.5) * H_src / H_dst - 0.5
//
// Sources of the same shape decode to u=x, v=y exactly, so equal-shape GT
// produces bit-identical reads to the previous direct-index path.
//
// Edge clamping: out-of-bound integer taps clamp to [0, src-1]. This is the
// standard "border" mode and keeps backward gradient delivery in-bounds.

#ifdef __CUDACC__

#include <cuda_runtime.h>


namespace _bilinear_detail {

// Resolve the four bilinear taps at destination (x_dst, y_dst) sampling
// source (W_src, H_src) from destination (W_dst, H_dst). Returns the four
// clamped integer taps + their bilinear weights.
__device__ __forceinline__ void resolve_taps(
    int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src,
    int& x0, int& y0, int& x1, int& y1,
    float& w00, float& w01, float& w10, float& w11
) {
    float u = ((float)x_dst + 0.5f) * (float)W_src / (float)W_dst - 0.5f;
    float v = ((float)y_dst + 0.5f) * (float)H_src / (float)H_dst - 0.5f;

    int x0_raw = (int)floorf(u);
    int y0_raw = (int)floorf(v);
    float fx = u - (float)x0_raw;
    float fy = v - (float)y0_raw;

    x0 = max(0, min(W_src - 1, x0_raw));
    y0 = max(0, min(H_src - 1, y0_raw));
    x1 = max(0, min(W_src - 1, x0_raw + 1));
    y1 = max(0, min(H_src - 1, y0_raw + 1));

    w00 = (1.0f - fx) * (1.0f - fy);
    w10 = fx        * (1.0f - fy);
    w01 = (1.0f - fx) * fy;
    w11 = fx        * fy;
}

} // namespace _bilinear_detail


// Bilinear sample of a [B, H, W, 1]-laid-out float source at the destination
// pixel (x_dst, y_dst) of a (H_dst, W_dst) target. b indexes the batch slice.
__device__ __forceinline__ float bilinear_sample_f(
    const float* __restrict__ src,
    int b, int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src
) {
    int x0, y0, x1, y1; float w00, w01, w10, w11;
    _bilinear_detail::resolve_taps(
        x_dst, y_dst, W_dst, H_dst, W_src, H_src,
        x0, y0, x1, y1, w00, w01, w10, w11);
    size_t plane = (size_t)b * (size_t)H_src * (size_t)W_src;
    float v00 = src[plane + (size_t)y0 * W_src + x0];
    float v10 = src[plane + (size_t)y0 * W_src + x1];
    float v01 = src[plane + (size_t)y1 * W_src + x0];
    float v11 = src[plane + (size_t)y1 * W_src + x1];
    return w00 * v00 + w10 * v10 + w01 * v01 + w11 * v11;
}

// Same for float3 (per-channel bilinear).
__device__ __forceinline__ float3 bilinear_sample_f3(
    const float3* __restrict__ src,
    int b, int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src
) {
    int x0, y0, x1, y1; float w00, w01, w10, w11;
    _bilinear_detail::resolve_taps(
        x_dst, y_dst, W_dst, H_dst, W_src, H_src,
        x0, y0, x1, y1, w00, w01, w10, w11);
    size_t plane = (size_t)b * (size_t)H_src * (size_t)W_src;
    float3 v00 = src[plane + (size_t)y0 * W_src + x0];
    float3 v10 = src[plane + (size_t)y0 * W_src + x1];
    float3 v01 = src[plane + (size_t)y1 * W_src + x0];
    float3 v11 = src[plane + (size_t)y1 * W_src + x1];
    float3 out;
    out.x = w00*v00.x + w10*v10.x + w01*v01.x + w11*v11.x;
    out.y = w00*v00.y + w10*v10.y + w01*v01.y + w11*v11.y;
    out.z = w00*v00.z + w10*v10.z + w01*v01.z + w11*v11.z;
    return out;
}

// Nearest-neighbor bool sampler — masks are boolean, bilinear has no
// meaningful interpolation. Returns the value of the nearest source pixel.
__device__ __forceinline__ bool nearest_sample_b(
    const bool* __restrict__ src,
    int b, int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src
) {
    float u = ((float)x_dst + 0.5f) * (float)W_src / (float)W_dst - 0.5f;
    float v = ((float)y_dst + 0.5f) * (float)H_src / (float)H_dst - 0.5f;
    int x = max(0, min(W_src - 1, (int)floorf(u + 0.5f)));
    int y = max(0, min(H_src - 1, (int)floorf(v + 0.5f)));
    size_t plane = (size_t)b * (size_t)H_src * (size_t)W_src;
    return src[plane + (size_t)y * W_src + x];
}


// Backward scatter: distribute `dvalue` across the four source taps using the
// same bilinear weights. Atomic — many destination pixels may share source
// taps when the destination is larger than the source.
__device__ __forceinline__ void bilinear_scatter_add_f(
    float* __restrict__ dst,
    int b, int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src,
    float dvalue
) {
    int x0, y0, x1, y1; float w00, w01, w10, w11;
    _bilinear_detail::resolve_taps(
        x_dst, y_dst, W_dst, H_dst, W_src, H_src,
        x0, y0, x1, y1, w00, w01, w10, w11);
    size_t plane = (size_t)b * (size_t)H_src * (size_t)W_src;
    if (w00 != 0.0f) atomicAdd(dst + plane + (size_t)y0 * W_src + x0, w00 * dvalue);
    if (w10 != 0.0f) atomicAdd(dst + plane + (size_t)y0 * W_src + x1, w10 * dvalue);
    if (w01 != 0.0f) atomicAdd(dst + plane + (size_t)y1 * W_src + x0, w01 * dvalue);
    if (w11 != 0.0f) atomicAdd(dst + plane + (size_t)y1 * W_src + x1, w11 * dvalue);
}

// Same for float3 — scatter per channel.
__device__ __forceinline__ void bilinear_scatter_add_f3(
    float3* __restrict__ dst,
    int b, int x_dst, int y_dst,
    int W_dst, int H_dst,
    int W_src, int H_src,
    float3 dvalue
) {
    int x0, y0, x1, y1; float w00, w01, w10, w11;
    _bilinear_detail::resolve_taps(
        x_dst, y_dst, W_dst, H_dst, W_src, H_src,
        x0, y0, x1, y1, w00, w01, w10, w11);
    size_t plane = (size_t)b * (size_t)H_src * (size_t)W_src;
    auto scatter = [&](size_t off, float w) {
        if (w == 0.0f) return;
        float* p = (float*)(dst + off);
        atomicAdd(p + 0, w * dvalue.x);
        atomicAdd(p + 1, w * dvalue.y);
        atomicAdd(p + 2, w * dvalue.z);
    };
    scatter(plane + (size_t)y0 * W_src + x0, w00);
    scatter(plane + (size_t)y0 * W_src + x1, w10);
    scatter(plane + (size_t)y1 * W_src + x0, w01);
    scatter(plane + (size_t)y1 * W_src + x1, w11);
}

#endif // __CUDACC__
