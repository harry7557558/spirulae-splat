// Decoded-video post-processing: YCbCr to RGB, and the sharpness metric that
// picks the least motion-blurred frame of a window.
//
// The decoder hands back multi-planar images; their planes are copied into
// plain buffers (one vkCmdCopyImageToBuffer per aspect) and read here through
// the same `uint*` convention as every other kernel -- 8-bit samples come out
// of a word four at a time, 10/12/16-bit two at a time, so no 8- or 16-bit
// storage feature is needed.

#include "_common.slang"

static const int WG = 256;

// flags bits, shared by every entry point here.
static const uint kF_TwoBytes = 1u << 0;   // 16-bit plane samples
static const uint kF_SubX     = 1u << 1;   // chroma is half-width
static const uint kF_SubY     = 1u << 2;   // chroma is half-height
static const uint kF_Planar   = 1u << 3;   // Cb and Cr in separate planes
static const uint kF_RotShift = 4u;        // 2 bits: 0/90/180/270 clockwise

// One luma/chroma sample. `stride` and `idx` are in texels, not bytes.
float load_sample(uint* buf, uint idx, uint flags, uint shift) {
    if ((flags & kF_TwoBytes) != 0u) {
        uint w = buf[idx >> 1u];
        uint v = (w >> ((idx & 1u) * 16u)) & 0xFFFFu;
        return float(v >> shift);
    }
    uint w = buf[idx >> 2u];
    return float((w >> ((idx & 3u) * 8u)) & 0xFFu);
}

// ---------------------------------------------------------------------------
// YCbCr -> RGB, with box downscale and rotation folded in
// ---------------------------------------------------------------------------

struct YuvParams {
    // Rows of the 3x4 conversion matrix, in [Y, Cb, Cr, 1] order. The host
    // bakes bit depth, studio/full range and the BT.601/709/2020 coefficients
    // into it, so the shader has one code path for every combination.
    float4 m0;
    float4 m1;
    float4 m2;
    uint*  out;      // ceil(out_w*out_h*3 / 4) words of packed RGB bytes
    uint*  luma;
    uint*  cb;
    uint*  cr;       // == cb for interleaved (2-plane) formats
    uint   out_w;
    uint   out_h;
    uint   src_w;
    uint   src_h;
    uint   luma_stride;    // texels per row
    uint   chroma_stride;  // texels per row (2 components per texel if !planar)
    uint   flags;
    uint   shift;
    float  inv_sx;         // source texels per scaled texel
    float  inv_sy;
    uint   groups_per_row;
};

// Four consecutive output pixels per thread: 12 bytes is exactly three words,
// so every store is a whole word (byte stores would need 8-bit storage).
[shader("compute")]
[numthreads(WG, 1, 1)]
void yuv_to_rgb(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                uniform YuvParams p) {
    uint quad = fold_index(gid, tid, p.groups_per_row, uint(WG));
    uint total = p.out_w * p.out_h;
    uint n_quads = (total + 3u) / 4u;
    if (quad >= n_quads) return;

    uint rot = (p.flags >> kF_RotShift) & 3u;
    uint sw = (rot & 1u) != 0u ? p.out_h : p.out_w;   // scaled-source dims
    uint sh = (rot & 1u) != 0u ? p.out_w : p.out_h;
    uint sub_x = (p.flags & kF_SubX) != 0u ? 1u : 0u;
    uint sub_y = (p.flags & kF_SubY) != 0u ? 1u : 0u;
    bool planar = (p.flags & kF_Planar) != 0u;

    uint words[3] = {0u, 0u, 0u};
    for (uint b = 0; b < 4u; ++b) {
        uint idx = quad * 4u + b;
        if (idx >= total) break;
        uint ox = idx % p.out_w;
        uint oy = idx / p.out_w;

        // Undo the rotation to reach coordinates in the scaled source.
        uint px, py;
        if (rot == 0u)      { px = ox;              py = oy; }
        else if (rot == 1u) { px = oy;              py = sh - 1u - ox; }
        else if (rot == 2u) { px = sw - 1u - ox;    py = sh - 1u - oy; }
        else                { px = sw - 1u - oy;    py = ox; }

        // Box-average the source footprint. At scale 1 this is a single texel.
        uint x0 = uint(float(px) * p.inv_sx);
        uint x1 = max(x0 + 1u, uint(float(px + 1u) * p.inv_sx));
        uint y0 = uint(float(py) * p.inv_sy);
        uint y1 = max(y0 + 1u, uint(float(py + 1u) * p.inv_sy));
        x1 = min(x1, p.src_w);
        y1 = min(y1, p.src_h);
        x0 = min(x0, x1 - 1u);
        y0 = min(y0, y1 - 1u);

        float y_sum = 0.0f;
        for (uint yy = y0; yy < y1; ++yy)
            for (uint xx = x0; xx < x1; ++xx)
                y_sum += load_sample(p.luma, yy * p.luma_stride + xx, p.flags, p.shift);
        float y_val = y_sum / float((x1 - x0) * (y1 - y0));

        uint cx0 = x0 >> sub_x, cy0 = y0 >> sub_y;
        uint cx1 = max(cx0 + 1u, (x1 + ((1u << sub_x) - 1u)) >> sub_x);
        uint cy1 = max(cy0 + 1u, (y1 + ((1u << sub_y) - 1u)) >> sub_y);
        float cb_sum = 0.0f, cr_sum = 0.0f;
        for (uint yy = cy0; yy < cy1; ++yy) {
            for (uint xx = cx0; xx < cx1; ++xx) {
                if (planar) {
                    cb_sum += load_sample(p.cb, yy * p.chroma_stride + xx, p.flags, p.shift);
                    cr_sum += load_sample(p.cr, yy * p.chroma_stride + xx, p.flags, p.shift);
                } else {
                    uint o = (yy * p.chroma_stride + xx) * 2u;
                    cb_sum += load_sample(p.cb, o, p.flags, p.shift);
                    cr_sum += load_sample(p.cb, o + 1u, p.flags, p.shift);
                }
            }
        }
        float inv_c = 1.0f / float((cx1 - cx0) * (cy1 - cy0));
        float4 v = float4(y_val, cb_sum * inv_c, cr_sum * inv_c, 1.0f);

        uint r = uint(clamp(dot(p.m0, v), 0.0f, 255.0f) + 0.5f);
        uint g = uint(clamp(dot(p.m1, v), 0.0f, 255.0f) + 0.5f);
        uint bl = uint(clamp(dot(p.m2, v), 0.0f, 255.0f) + 0.5f);

        uint byte0 = b * 3u;
        words[byte0 >> 2u] |= r << ((byte0 & 3u) * 8u);
        uint byte1 = byte0 + 1u;
        words[byte1 >> 2u] |= g << ((byte1 & 3u) * 8u);
        uint byte2 = byte0 + 2u;
        words[byte2 >> 2u] |= bl << ((byte2 & 3u) * 8u);
    }
    p.out[quad * 3u + 0u] = words[0];
    p.out[quad * 3u + 1u] = words[1];
    p.out[quad * 3u + 2u] = words[2];
}

// ---------------------------------------------------------------------------
// Sharpness: variance of the Laplacian on a fixed-size grey thumbnail
// ---------------------------------------------------------------------------
//
// The reference (reference/scripts/extract_frames.py) resizes to 512x512,
// takes the grey channel and reports var(Laplacian). Reproducing the size
// exactly keeps the numbers comparable; luma alone stands in for grey,
// which saves copying the chroma planes for a frame that may never be
// written out.

struct ThumbParams {
    float* out;      // S*S
    uint*  luma;
    uint   src_w;
    uint   src_h;
    uint   luma_stride;
    uint   size;     // S
    uint   flags;
    uint   shift;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void luma_thumbnail(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                    uniform ThumbParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.size * p.size) return;
    uint ox = i % p.size;
    uint oy = i / p.size;

    uint x0 = (ox * p.src_w) / p.size;
    uint x1 = max(x0 + 1u, ((ox + 1u) * p.src_w) / p.size);
    uint y0 = (oy * p.src_h) / p.size;
    uint y1 = max(y0 + 1u, ((oy + 1u) * p.src_h) / p.size);
    x1 = min(x1, p.src_w);
    y1 = min(y1, p.src_h);

    float sum = 0.0f;
    for (uint yy = y0; yy < y1; ++yy)
        for (uint xx = x0; xx < x1; ++xx)
            sum += load_sample(p.luma, yy * p.luma_stride + xx, p.flags, p.shift);
    p.out[i] = sum / float((x1 - x0) * (y1 - y0));
}

struct VarianceParams {
    float* out;      // one float, at index `slot`
    float* thumb;
    uint   size;
    uint   slot;
};

groupshared float g_sum[WG];
groupshared float g_sqr[WG];

// One workgroup for the whole 512x512 thumbnail: 256 threads x 1024 pixels is
// microseconds, and a single dispatch keeps the per-frame command count at two.
[shader("compute")]
[numthreads(WG, 1, 1)]
void laplacian_variance(uint3 tid : SV_GroupThreadID, uniform VarianceParams p) {
    uint t = tid.x;
    uint S = p.size;
    float sum = 0.0f, sqr = 0.0f;
    uint n = 0u;
    for (uint i = t; i < S * S; i += uint(WG)) {
        uint x = i % S;
        uint y = i / S;
        if (x == 0u || y == 0u || x + 1u >= S || y + 1u >= S) continue;
        float c = p.thumb[i];
        float lap = p.thumb[i - 1u] + p.thumb[i + 1u] + p.thumb[i - S] + p.thumb[i + S] -
                    4.0f * c;
        sum += lap;
        sqr += lap * lap;
        ++n;
    }
    g_sum[t] = sum;
    g_sqr[t] = sqr;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = uint(WG) / 2u; s > 0u; s >>= 1u) {
        if (t < s) {
            g_sum[t] += g_sum[t + s];
            g_sqr[t] += g_sqr[t + s];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (t == 0u) {
        float count = float((S - 2u) * (S - 2u));
        float mean = g_sum[0] / count;
        p.out[p.slot] = g_sqr[0] / count - mean * mean;
    }
}
