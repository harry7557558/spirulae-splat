// Dataset GT warp + byte->float conversion kernels (Vulkan backend),
// mirroring kernels/pixelwise/ImageWarp.cu + GtDepthNormalWarp.cu (launch_warp_*) and
// the raw converters (uint8/16_*_to_float_raw). The wide warps project each
// pinhole-face ray through the same canonical projection_utils.slang
// *_proj_nav exports the CUDA kernels call, so per-pixel outputs differ only
// by fast-math rounding (plus isolated valid-flip / nearest-rounding pixels
// at projection boundaries).
//
// Layout conventions match the CUDA launchers exactly: input byte/float
// buffers are contiguous [B, Hin, Win, C], warp outputs are contiguous
// [B, K, Hout, Wout, C] float. Sub-word (u8/u16) elements are extracted from
// u32 words -- 8/16-bit storage access is not baseline SPIR-V -- and the
// mask warp writes whole u32 words (4 output bytes per thread; allocations
// are 16-byte rounded so the trailing partial word stays in bounds).
//
// The camera model is a plain runtime params field, exactly like the CUDA
// kernels' `CameraModelType` argument (all loads go through the same
// pointers on every path, so the branch is speculation-safe on llvmpipe).

#define CudaDeviceExport ForceInline
#include "int64_compat.slang"
#include "int8_compat.slang"
#include "shaders/projection_utils.slang"
#undef CudaDeviceExport

// Input element kinds (the CUDA T_in template axis / elem_size dispatch).
static const int kElemU8 = 0;
static const int kElemU16 = 1;
static const int kElemF32 = 2;

static const float kWarpPi = 3.14159265358979323846f;

// Load element `idx` of the raw input buffer as an un-normalized float.
[ForceInline]
float _wp_load(uint* buf, gindex_t idx, int kind) {
    if (kind == kElemU8)
        return float(u8_load(buf, idx));
    if (kind == kElemU16) {
        uint w = buf[idx >> 1];
        return float((w >> ((uint(idx) & 1u) * 16u)) & 0xffffu);
    }
    return asfloat(buf[idx]);
}

[ForceInline]
float _wp_tap(uint* buf, int kind, gindex_t plane, int H, int W, int C,
              int cid, int x_, int y_, float norm_inv, float padding) {
    if (x_ < 0 || x_ >= W || y_ < 0 || y_ >= H) return padding;
    gindex_t idx = (plane + gindex_t(y_) * W + x_) * C + cid;
    return _wp_load(buf, idx, kind) * norm_inv;
}

[ForceInline]
int _wp_wrap_x(int v, int W) {
    int r = v % W;
    if (r < 0) r += W;
    return r;
}

// bilinear_byte_norm / bilinear_byte_norm_wrap_u mirror: half-open taps with
// a constant padding value outside the image (wrap_u wraps x and clamps y,
// keeping the equirectangular seam at u = 0 / u = W continuous). `plane` is
// the [bid] pixel offset (bid * H * W). Non-wrap coordinates are pre-clamped
// to a window that preserves the all-taps-padding result -- float->int
// conversion of an unbounded projected uv is undefined, and the CUDA kernels
// only avoid it by accident of the `valid` guard.
[ForceInline]
float _wp_bilinear(uint* buf, int kind, gindex_t plane, int H, int W, int C,
                   int cid, float x, float y, float norm_inv, float padding,
                   bool wrap_u) {
    if (!wrap_u) {
        x = clamp(x, -8.0f, float(W) + 8.0f);
        y = clamp(y, -8.0f, float(H) + 8.0f);
    }
    int x0 = int(floor(x)), x1 = x0 + 1;
    int y0 = int(floor(y)), y1 = y0 + 1;
    float wx1 = x - float(x0), wx0 = 1.0f - wx1;
    float wy1 = y - float(y0), wy0 = 1.0f - wy1;
    if (wrap_u) {
        x0 = _wp_wrap_x(x0, W);
        x1 = _wp_wrap_x(x1, W);
    }
    return _wp_tap(buf, kind, plane, H, W, C, cid, x0, y0, norm_inv, padding) * (wx0 * wy0)
         + _wp_tap(buf, kind, plane, H, W, C, cid, x1, y0, norm_inv, padding) * (wx1 * wy0)
         + _wp_tap(buf, kind, plane, H, W, C, cid, x0, y1, norm_inv, padding) * (wx0 * wy1)
         + _wp_tap(buf, kind, plane, H, W, C, cid, x1, y1, norm_inv, padding) * (wx1 * wy1);
}

[ForceInline]
float4 _wp_intrin(float* intrins, int bid) {
    return float4(intrins[bid * 4 + 0], intrins[bid * 4 + 1],
                  intrins[bid * 4 + 2], intrins[bid * 4 + 3]);
}

[ForceInline]
float[10] _wp_dist(float* dist_coeffs, int bid) {
    float[10] dc;
    [ForceUnroll]
    for (int k = 0; k < 10; k++) dc[k] = dist_coeffs[bid * 10 + k];
    return dc;
}

[ForceInline]
void _wp_axes(float* axes, int ki, out float3 ax, out float3 ay,
              out float3 az) {
    int b = 9 * ki;
    ax = float3(axes[b + 0], axes[b + 1], axes[b + 2]);
    ay = float3(axes[b + 3], axes[b + 4], axes[b + 5]);
    az = float3(axes[b + 6], axes[b + 7], axes[b + 8]);
}

// Camera-model dispatch matching the CUDA kernels' ternary chain.
[ForceInline]
bool _wp_proj(int camera_model, float3 raydir, float4 intrin,
              float[10] dist_coeffs, out float2 uv) {
    if (camera_model == int(CameraModelType::FISHEYE))
        return fisheye_proj_nav(raydir, intrin, dist_coeffs, uv);
    if (camera_model == int(CameraModelType::EQUISOLID))
        return equisolid_proj_nav(raydir, intrin, dist_coeffs, uv);
    return persp_proj_nav(raydir, intrin, dist_coeffs, uv);
}

// Equirectangular sample position for a face ray (see ImageWarp.cu; f is
// w * (0.5 / pi)).
[ForceInline]
float2 _wp_equi_uv(float3 raydir, int h, int w) {
    float f = float(w) * (0.5f / kWarpPi);
    return float2(
        0.5f * float(w) + f * atan2(raydir.x, raydir.z),
        0.5f * float(h) + f * atan2(raydir.y, length(float2(raydir.x, raydir.z))));
}

// _wide_depth_to_face_ray_depth mirror: convert a sampled wide-frame depth
// to per-face RAY depth. Linear (z) input is divided by cos(theta) w.r.t.
// the INPUT optical axis; rays at/behind the input +z plane are dropped.
// d <= 0 is the invalid sentinel.
[ForceInline]
float _wp_face_ray_depth(float d, float3 raydir, bool input_is_ray_depth) {
    if (d <= 0.0f) return 0.0f;
    if (input_is_ray_depth) return d;
    float rl = length(raydir);
    float rz = raydir.z;
    if (rz <= 1e-6f * rl) return 0.0f;
    return d * rl / rz;
}

// _warp_one_normal mirror: sample a wide-frame normal (INPUT camera frame)
// and rotate it into the face frame; (-1,-1,-1) sentinel for "no data" /
// degenerate samples.
[ForceInline]
float3 _wp_one_normal(uint* src, int kind, gindex_t plane, int H, int W,
                      float x, float y, float norm_inv, float decode_off,
                      bool wrap_u, float3 axis_x, float3 axis_y,
                      float3 axis_z) {
    float3 n;
    n.x = _wp_bilinear(src, kind, plane, H, W, 3, 0, x, y, norm_inv, 0.0f, wrap_u) + decode_off;
    n.y = _wp_bilinear(src, kind, plane, H, W, 3, 1, x, y, norm_inv, 0.0f, wrap_u) + decode_off;
    n.z = _wp_bilinear(src, kind, plane, H, W, 3, 2, x, y, norm_inv, 0.0f, wrap_u) + decode_off;
    if (n.x + n.y + n.z <= -2.366f) return float3(-1.0f);
    float3 ax = normalize(axis_x), ay = normalize(axis_y), az = normalize(axis_z);
    float3 r = float3(dot(ax, n), dot(ay, n), dot(az, n));
    float rl = length(r);
    if (rl <= 1e-8f) return float3(-1.0f);
    return r / rl;
}

// Shared params for the image / depth / normal warp entries. 96 bytes:
// pushed directly (by value). Pointer fields are never null (zeros
// fallback); the *_equi entries statically ignore intrins / dist_coeffs.
struct WarpParams {
    float* intrins;      // [B, 4]
    float* dist_coeffs;  // [B, 10]
    uint* src;           // [B, Hin, Win, C] u8 / u16 / f32 elements
    float* dst;          // [B, K, Hout, Wout, C] float
    float* axes;         // [K, 3, 3] face axes
    int B;
    int Hin;
    int Win;
    int C;
    int K;
    int Hout;
    int Wout;
    int elem_kind;     // kElemU8 / kElemU16 / kElemF32
    int camera_model;  // CameraModelType (wide entries)
    int in_H;          // intrinsics reference resolution (the projected uv is
    int in_W;          //   in this pixel space; depth/normal rescale to Hin/Win)
    int ray_depth;     // depth entries: input already stores ray depth
    float norm_inv;
    float decode_off;  // normal entries: post-normalization offset
};

// 16x16 logical tile decode: flat 256-wide workgroups (the pinned 32-wide
// subgroup divides 256), grid = (ceil(Wout/16), ceil(Hout/16), B).
[ForceInline]
bool _wp_tile(WarpParams p, uint3 wg, uint tid, out int i, out int j,
              out int bid, out float tx, out float ty) {
    i = int(wg.x) * 16 + int(tid % 16u);
    j = int(wg.y) * 16 + int(tid / 16u);
    bid = int(wg.z);
    tx = -1.0f + 2.0f * (float(i) + 0.5f) / float(p.Wout);
    ty = -1.0f + 2.0f * (float(j) + 0.5f) / float(p.Hout);
    return i < p.Wout && j < p.Hout;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_img_wide(uniform WarpParams p, uint3 wg: SV_GroupID,
                   uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    float4 intrin = _wp_intrin(p.intrins, bid);
    float[10] dc = _wp_dist(p.dist_coeffs, bid);
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv;
        bool valid = _wp_proj(p.camera_model, raydir, intrin, dc, uv);
        gindex_t base = (((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i) * p.C;
        for (int c = 0; c < p.C; c++) {
            float v = 0.5f;
            if (valid)
                v = _wp_bilinear(p.src, p.elem_kind, plane, p.Hin, p.Win,
                                 p.C, c, uv.x, uv.y, p.norm_inv, 0.5f, false);
            p.dst[base + c] = v;
        }
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_img_equi(uniform WarpParams p, uint3 wg: SV_GroupID,
                   uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv = _wp_equi_uv(raydir, p.Hin, p.Win);
        gindex_t base = (((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i) * p.C;
        for (int c = 0; c < p.C; c++)
            p.dst[base + c] = _wp_bilinear(p.src, p.elem_kind, plane, p.Hin,
                                           p.Win, p.C, c, uv.x, uv.y,
                                           p.norm_inv, 0.5f, true);
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_depth_wide(uniform WarpParams p, uint3 wg: SV_GroupID,
                     uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    float4 intrin = _wp_intrin(p.intrins, bid);
    float[10] dc = _wp_dist(p.dist_coeffs, bid);
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    // Projected uv is in intrinsics-reference (RGB input) pixel space;
    // rescale to the GT depth map's own resolution.
    float sx = float(p.Win) / float(p.in_W), sy = float(p.Hin) / float(p.in_H);
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv;
        bool valid = _wp_proj(p.camera_model, raydir, intrin, dc, uv);
        float outv = 0.0f;
        if (valid) {
            float d = _wp_bilinear(p.src, p.elem_kind, plane, p.Hin, p.Win, 1,
                                   0, uv.x * sx, uv.y * sy, p.norm_inv, 0.0f,
                                   false);
            outv = _wp_face_ray_depth(d, raydir, p.ray_depth != 0);
        }
        p.dst[((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i] = outv;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_depth_equi(uniform WarpParams p, uint3 wg: SV_GroupID,
                     uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv = _wp_equi_uv(raydir, p.Hin, p.Win);
        float d = _wp_bilinear(p.src, p.elem_kind, plane, p.Hin, p.Win, 1, 0,
                               uv.x, uv.y, p.norm_inv, 0.0f, true);
        p.dst[((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i] =
            _wp_face_ray_depth(d, raydir, p.ray_depth != 0);
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_normal_wide(uniform WarpParams p, uint3 wg: SV_GroupID,
                      uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    float4 intrin = _wp_intrin(p.intrins, bid);
    float[10] dc = _wp_dist(p.dist_coeffs, bid);
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    float sx = float(p.Win) / float(p.in_W), sy = float(p.Hin) / float(p.in_H);
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv;
        bool valid = _wp_proj(p.camera_model, raydir, intrin, dc, uv);
        float3 nf = float3(-1.0f);
        if (valid)
            nf = _wp_one_normal(p.src, p.elem_kind, plane, p.Hin, p.Win,
                                uv.x * sx, uv.y * sy, p.norm_inv,
                                p.decode_off, false, ax, ay, az);
        gindex_t base = (((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i) * 3;
        p.dst[base + 0] = nf.x;
        p.dst[base + 1] = nf.y;
        p.dst[base + 2] = nf.z;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_normal_equi(uniform WarpParams p, uint3 wg: SV_GroupID,
                      uint tid: SV_GroupIndex) {
    int i, j, bid;
    float tx, ty;
    if (!_wp_tile(p, wg, tid, i, j, bid, tx, ty)) return;
    gindex_t plane = gindex_t(bid) * p.Hin * p.Win;
    for (int ki = 0; ki < p.K; ki++) {
        float3 ax, ay, az;
        _wp_axes(p.axes, ki, ax, ay, az);
        float3 raydir = az + tx * ax + ty * ay;
        float2 uv = _wp_equi_uv(raydir, p.Hin, p.Win);
        float3 nf = _wp_one_normal(p.src, p.elem_kind, plane, p.Hin, p.Win,
                                   uv.x, uv.y, p.norm_inv, p.decode_off,
                                   true, ax, ay, az);
        gindex_t base = (((gindex_t(bid) * p.K + ki) * p.Hout + j) * p.Wout + i) * 3;
        p.dst[base + 0] = nf.x;
        p.dst[base + 1] = nf.y;
        p.dst[base + 2] = nf.z;
    }
}

// Nearest-neighbor mask warps (uint8 0/nonzero in, uint8 0/1 out). Threads
// map to whole u32 OUTPUT words in the flat [B, K, Hout, Wout] byte space:
// each thread decodes its 4 byte positions independently and stores one
// word, so there is exactly one writer per word and no atomics.
struct WarpMaskParams {  // 72 bytes: pushed directly
    float* intrins;      // [B, 4]
    float* dist_coeffs;  // [B, 10]
    uint* src;           // [B, Hin, Win] uint8
    uint* dst;           // [B, K, Hout, Wout] uint8 (word writes)
    float* axes;         // [K, 3, 3]
    int B;
    int Hin;
    int Win;
    int K;
    int Hout;
    int Wout;
    int camera_model;
    uint wgs_per_row;
};

[ForceInline]
uint _wp_mask_byte(WarpMaskParams p, uint n, bool equi) {
    uint i = n % uint(p.Wout);
    uint t = n / uint(p.Wout);
    uint j = t % uint(p.Hout);
    t /= uint(p.Hout);
    uint ki = t % uint(p.K);
    uint bid = t / uint(p.K);
    float tx = -1.0f + 2.0f * (float(i) + 0.5f) / float(p.Wout);
    float ty = -1.0f + 2.0f * (float(j) + 0.5f) / float(p.Hout);
    float3 ax, ay, az;
    _wp_axes(p.axes, int(ki), ax, ay, az);
    float3 raydir = az + tx * ax + ty * ay;
    int xs, ys;
    if (equi) {
        float2 uv = _wp_equi_uv(raydir, p.Hin, p.Win);
        xs = int(floor(uv.x + 0.5f));
        ys = int(floor(uv.y + 0.5f));
        // wrap u for equirectangular (panoramic continuity); clamp v.
        if (xs < 0) xs += p.Win; else if (xs >= p.Win) xs -= p.Win;
        ys = clamp(ys, 0, p.Hin - 1);
    } else {
        float4 intrin = _wp_intrin(p.intrins, int(bid));
        float[10] dc = _wp_dist(p.dist_coeffs, int(bid));
        float2 uv;
        if (!_wp_proj(p.camera_model, raydir, intrin, dc, uv)) return 0u;
        xs = int(floor(clamp(uv.x + 0.5f, -8.0f, float(p.Win) + 8.0f)));
        ys = int(floor(clamp(uv.y + 0.5f, -8.0f, float(p.Hin) + 8.0f)));
        if (xs < 0 || xs >= p.Win || ys < 0 || ys >= p.Hin) return 0u;
    }
    gindex_t idx = gindex_t(bid) * p.Hin * p.Win + gindex_t(ys) * p.Win + xs;
    return _wp_load(p.src, idx, kElemU8) != 0.0f ? 1u : 0u;
}

[ForceInline]
void _wp_mask_main(WarpMaskParams p, uint3 wg, uint tid, bool equi) {
    uint word = (wg.y * p.wgs_per_row + wg.x) * 256u + tid;
    uint total = uint(p.B * p.K) * uint(p.Hout) * uint(p.Wout);
    if (word >= (total + 3u) / 4u) return;
    uint outw = 0u;
    for (uint b = 0; b < 4u; b++) {
        uint n = word * 4u + b;
        if (n < total)
            outw |= _wp_mask_byte(p, n, equi) << (b * 8u);
    }
    p.dst[word] = outw;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_mask_wide(uniform WarpMaskParams p, uint3 wg: SV_GroupID,
                    uint tid: SV_GroupIndex) {
    _wp_mask_main(p, wg, tid, false);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void warp_mask_equi(uniform WarpMaskParams p, uint3 wg: SV_GroupID,
                    uint tid: SV_GroupIndex) {
    _wp_mask_main(p, wg, tid, true);
}

// Flat byte/u16 -> float conversion: out[i] = elem(i) * scale + offset.
// Covers uint8/16_image_to_float (1/255, 1/65535), uint8_normal_to_float
// (1/127.5, -1) and uint16_depth_to_float (1, 0).
struct BytesToFloatParams {  // 40 bytes: pushed directly
    uint* src;
    float* dst;
    uint total;
    int elem_kind;
    float scale;
    float offset;
    uint wgs_per_row;
    int _pad0;
};

[shader("compute")]
[numthreads(256, 1, 1)]
void bytes_to_float(uniform BytesToFloatParams p, uint3 wg: SV_GroupID,
                    uint tid: SV_GroupIndex) {
    uint idx = (wg.y * p.wgs_per_row + wg.x) * 256u + tid;
    if (idx >= p.total) return;
    p.dst[idx] = _wp_load(p.src, gindex_t(idx), p.elem_kind) * p.scale + p.offset;
}
