// Spatial resampling on channel-last [H, W, C] maps.
//
// The bilinear kernel follows torch's `align_corners=False` mapping exactly --
// src = (dst + 0.5) * scale - 0.5, clamped, with the low index taken by floor
// and the high index clamped to H-1. Mask upsampling runs through it on every
// frame of a video, and an off-by-half here shifts every boundary by a pixel.

#include "_common.slang"

// Spec ids 0..1.
[SpecializationConstant] const int kInputF16 = 0;
[SpecializationConstant] const int kAct      = 0;

static const int WG = 256;

struct ResizeParams {
    float* out;
    uint*  x;
    uint   Ho;
    uint   Wo;
    uint   Hi;
    uint   Wi;
    uint   C;
    uint   groups_per_row;
};

[shader("compute")]
[numthreads(WG, 1, 1)]
void resize_bilinear(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                     uniform ResizeParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.Ho * p.Wo * p.C) return;
    uint c = i % p.C;
    uint rest = i / p.C;
    uint xo = rest % p.Wo;
    uint yo = rest / p.Wo;

    float sy = (float(yo) + 0.5f) * (float(p.Hi) / float(p.Ho)) - 0.5f;
    float sx = (float(xo) + 0.5f) * (float(p.Wi) / float(p.Wo)) - 0.5f;
    sy = max(sy, 0.0f);
    sx = max(sx, 0.0f);
    p.out[i] = sample_bilinear(p.x, kInputF16 != 0, int(p.Hi), int(p.Wi), int(p.C),
                               int(c), sy, sx);
}

// Nearest-neighbour 2x upsample -- what the seg head's PixelDecoder uses
// between FPN levels (F.interpolate(mode="nearest")).
[shader("compute")]
[numthreads(WG, 1, 1)]
void upsample_nearest2x(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                        uniform ResizeParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.Ho * p.Wo * p.C) return;
    uint c = i % p.C;
    uint rest = i / p.C;
    uint xo = rest % p.Wo;
    uint yo = rest / p.Wo;
    uint ys = min(yo >> 1, p.Hi - 1u);
    uint xs = min(xo >> 1, p.Wi - 1u);
    p.out[i] = load_elem(p.x, (ys * p.Wi + xs) * p.C + c, kInputF16 != 0);
}

// 2x2 max pool, stride 2 -- the neck's 0.5x scale.
[shader("compute")]
[numthreads(WG, 1, 1)]
void maxpool2x2(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                uniform ResizeParams p) {
    uint i = fold_index(gid, tid, p.groups_per_row, uint(WG));
    if (i >= p.Ho * p.Wo * p.C) return;
    uint c = i % p.C;
    uint rest = i / p.C;
    uint xo = rest % p.Wo;
    uint yo = rest / p.Wo;
    float best = -3.402823466e+38f;
    for (uint dy = 0; dy < 2u; ++dy) {
        uint y = yo * 2u + dy;
        if (y >= p.Hi) continue;
        for (uint dx = 0; dx < 2u; ++dx) {
            uint x = xo * 2u + dx;
            if (x >= p.Wi) continue;
            best = max(best, load_elem(p.x, (y * p.Wi + x) * p.C + c, kInputF16 != 0));
        }
    }
    p.out[i] = best;
}

// ---------------------------------------------------------------------------
// Mask export: resize + threshold + pack to bytes, in one pass
// ---------------------------------------------------------------------------

struct MaskExportParams {
    uint*  out;    // ceil(Ho*Wo/4) u32 words, 4 mask bytes each
    uint*  x;      // [Hi, Wi] single-channel logits
    uint   Ho;
    uint   Wo;
    uint   Hi;
    uint   Wi;
    float  threshold;
    uint   groups_per_row;
};

// A tracked video frame can carry 200 instance masks at 1920x1080; pulling
// those back as floats would be 1.6 GB. Thresholding on device and returning
// packed bytes (0 / 255) makes it 400 MB of u8, and the host never sees a
// float mask. Four consecutive x per thread so the store is a whole word --
// byte stores would need the 8-bit-storage feature, which is off-baseline.
[shader("compute")]
[numthreads(WG, 1, 1)]
void resize_binarize(uint3 gid : SV_GroupID, uint3 tid : SV_GroupThreadID,
                     uniform MaskExportParams p) {
    uint w = fold_index(gid, tid, p.groups_per_row, uint(WG));
    uint total = p.Ho * p.Wo;
    uint n_words = (total + 3u) / 4u;
    if (w >= n_words) return;

    uint packed = 0u;
    for (uint b = 0; b < 4u; ++b) {
        uint idx = w * 4u + b;
        if (idx >= total) break;
        uint xo = idx % p.Wo;
        uint yo = idx / p.Wo;
        float sy = max((float(yo) + 0.5f) * (float(p.Hi) / float(p.Ho)) - 0.5f, 0.0f);
        float sx = max((float(xo) + 0.5f) * (float(p.Wi) / float(p.Wo)) - 0.5f, 0.0f);
        float v = sample_bilinear(p.x, kInputF16 != 0, int(p.Hi), int(p.Wi), 1, 0, sy, sx);
        if (v > p.threshold) packed |= 255u << (b * 8u);
    }
    p.out[w] = packed;
}
