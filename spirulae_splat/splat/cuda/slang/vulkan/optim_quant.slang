// Block-wise quantized-state codecs (Vulkan backend) mirroring the CUDA
// device codecs in csrc/Tensor.h (QuantizedAdamState, QuantizedTensor) and
// csrc/GradQuant.cuh (gradq::Codec<8>). All packed storage is read and
// written as u32 WORDS (no 8/16-bit storage features; README rules). Cell
// indices are u32 — launchers validate numel < 2^32.
//
// Word-write cooperation: cells narrower than a word are assembled in
// groupshared (_oq_pack) and stored as whole words by the lowest thread of
// each word group. A 256-cell quant block always starts word-aligned, and
// the up-to-3 tail bytes past the last cell land in the allocation's 4-byte
// rounding (backend::device_malloc guarantees it). Unlike CUDA — which
// leaves tail bytes untouched — the tail bytes inside the final word are
// written as zero; nothing reads them.
//
// log1p/expm1 use the classic compensated forms, tracking CUDA's
// log1pf/expm1f closely but not bit-exactly; parity tests allow +-1
// quantization step on packed cells.

static const float kAdamEps = 1e-15f;

[ForceInline]
float _oq_log1p(float x) {
    float u = 1.0f + x;
    if (u == 1.0f) return x;
    return log(u) * (x / (u - 1.0f));
}

[ForceInline]
float _oq_expm1(float x) {
    float u = exp(x);
    if (u == 1.0f) return x;
    return (u - 1.0f) * (x / log(u));
}

// floor(x + 0.5) matches CUDA roundf everywhere the callers' subsequent
// clamp-to-[0, qmax] doesn't already coalesce the difference (the only
// mismatch is exactly -0.5, which both paths clamp to 0).
[ForceInline]
float _oq_round(float x) { return floor(x + 0.5f); }

// Scratch shared by the word-assembly stores and the block reductions.
// Every helper both begins and ends its groupshared use with a barrier, so
// they are safely callable in sequence within one kernel.
groupshared uint _oq_pack[256];
groupshared float4 _oq_red4[256];

// Block-wide (min, max, min, max) reduction over all 256 threads,
// broadcast to every thread. Deterministic serial combine by thread 0 (no
// subgroup-id assumptions; the quantized-optimizer kernels are
// memory-bound). Pass sentinels (+1e30, -1e30, ...) for threads that must
// not contribute.
[ForceInline]
float4 oq_block_reduce_mm4(float4 v, uint tid) {
    GroupMemoryBarrierWithGroupSync();
    _oq_red4[tid] = v;
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0) {
        float4 m = _oq_red4[0];
        for (uint i = 1; i < 256; i++) {
            float4 t = _oq_red4[i];
            m = float4(min(m.x, t.x), max(m.y, t.y),
                       min(m.z, t.z), max(m.w, t.w));
        }
        _oq_red4[0] = m;
    }
    GroupMemoryBarrierWithGroupSync();
    return _oq_red4[0];
}

// float2 (min, max) variant over the same scratch.
[ForceInline]
float2 oq_block_reduce_mm2(float2 v, uint tid) {
    float4 r = oq_block_reduce_mm4(float4(v.x, v.y, 1e30f, -1e30f), tid);
    return float2(r.x, r.y);
}

// ---------------------------------------------------------------------------
// QuantizedAdamState<BITS, 256> — joint (u, log_s) Adam state, AoS.
//   BITS=16: u32 word/cell   BITS=8: halfword/cell   BITS=4: byte/cell
// ---------------------------------------------------------------------------

[ForceInline]
float qadam_qmax<let BITS : int>() { return float((1 << BITS) - 1); }

[ForceInline]
uint2 qadam_read_q<let BITS : int>(uint* packed, uint idx) {
    if (BITS == 16) {
        uint w = packed[idx];
        return uint2(w & 0xffffu, w >> 16);
    } else if (BITS == 8) {
        uint w = packed[idx >> 1];
        uint h = (w >> (16u * (idx & 1u))) & 0xffffu;
        return uint2(h & 0xffu, h >> 8);
    } else {  // BITS == 4
        uint w = packed[idx >> 2];
        uint b = (w >> (8u * (idx & 3u))) & 0xffu;
        return uint2(b & 0xfu, b >> 4);
    }
}

[ForceInline]
float2 qadam_decode_us<let BITS : int>(uint* packed, uint idx, float4 mm) {
    uint2 q = qadam_read_q<BITS>(packed, idx);
    float inv = 1.0f / qadam_qmax<BITS>();
    return float2(mm.x + (mm.y - mm.x) * (float(q.x) * inv),
                  mm.z + (mm.w - mm.z) * (float(q.y) * inv));
}

[ForceInline]
float2 qadam_decode_g1g2<let BITS : int>(uint* packed, uint idx, float4 mm) {
    float2 prim = qadam_decode_us<BITS>(packed, idx, mm);
    float sqrt_g2 = kAdamEps * _oq_expm1(prim.y);
    return float2(prim.x * (sqrt_g2 + kAdamEps), sqrt_g2 * sqrt_g2);
}

// (g1, g2) -> linearly-quantized primitives (u, log_s); the block bounds
// reduction operates on these.
[ForceInline]
float2 qadam_g1g2_to_us(float g1, float g2) {
    float sqrt_g2 = sqrt(max(g2, 0.0f));
    float u = g1 / (sqrt_g2 + kAdamEps);
    return float2(u, _oq_log1p(max(sqrt_g2, 0.0f) * (1.0f / kAdamEps)));
}

[ForceInline]
uint2 qadam_quantize_us<let BITS : int>(float u_val, float log_s_val,
                                        float4 mm) {
    float qmax = qadam_qmax<BITS>();
    float u_range = max(mm.y - mm.x, kAdamEps);
    float s_range = max(mm.w - mm.z, kAdamEps);
    float u_qf = clamp(_oq_round(qmax * (u_val - mm.x) / u_range), 0.0f, qmax);
    float s_qf = clamp(_oq_round(qmax * (log_s_val - mm.z) / s_range), 0.0f,
                       qmax);
    return uint2(uint(u_qf), uint(s_qf));
}

// Cooperative packed store of one cell's (u_q, s_q). Call from ALL 256
// threads (barriers inside); `inside` gates the actual memory writes.
[ForceInline]
void qadam_store<let BITS : int>(uint* packed, uint idx, uint tid,
                                 bool inside, uint2 q) {
    if (BITS == 16) {
        if (inside) packed[idx] = q.x | (q.y << 16);
    } else if (BITS == 8) {
        GroupMemoryBarrierWithGroupSync();
        _oq_pack[tid] = inside ? (q.x | (q.y << 8)) : 0u;
        GroupMemoryBarrierWithGroupSync();
        if ((tid & 1u) == 0u && inside)
            packed[idx >> 1] = _oq_pack[tid] | (_oq_pack[tid + 1] << 16);
        GroupMemoryBarrierWithGroupSync();
    } else {  // BITS == 4
        GroupMemoryBarrierWithGroupSync();
        _oq_pack[tid] = inside ? ((q.x & 0xfu) | ((q.y & 0xfu) << 4)) : 0u;
        GroupMemoryBarrierWithGroupSync();
        if ((tid & 3u) == 0u && inside)
            packed[idx >> 2] = _oq_pack[tid] | (_oq_pack[tid + 1] << 8) |
                               (_oq_pack[tid + 2] << 16) |
                               (_oq_pack[tid + 3] << 24);
        GroupMemoryBarrierWithGroupSync();
    }
}

// ---------------------------------------------------------------------------
// QuantizedTensor<BITS, 256> — linear value codec.
//   BITS=16: halfword/cell   BITS=8: byte/cell
// (Same codec the SH value-quant render path decodes via sh_quant.slang.)
// ---------------------------------------------------------------------------

[ForceInline]
float qtensor_decode_v<let BITS : int>(uint* packed, uint idx, float2 mm) {
    float q;
    if (BITS == 16) {
        uint w = packed[idx >> 1];
        q = float((w >> (16u * (idx & 1u))) & 0xffffu);
    } else {  // BITS == 8
        uint w = packed[idx >> 2];
        q = float((w >> (8u * (idx & 3u))) & 0xffu);
    }
    return mm.x + (mm.y - mm.x) * (q / float((1 << BITS) - 1));
}

[ForceInline]
uint qtensor_quantize_v<let BITS : int>(float v, float2 mm) {
    float qmax = float((1 << BITS) - 1);
    float range = max(mm.y - mm.x, 1e-30f);
    return uint(clamp(_oq_round(qmax * (v - mm.x) / range), 0.0f, qmax));
}

// Cooperative packed store of one cell's value code. Call from ALL 256
// threads (barriers inside); `inside` gates the actual memory writes.
[ForceInline]
void qtensor_store<let BITS : int>(uint* packed, uint idx, uint tid,
                                   bool inside, uint q) {
    GroupMemoryBarrierWithGroupSync();
    _oq_pack[tid] = inside ? q : 0u;
    GroupMemoryBarrierWithGroupSync();
    if (BITS == 16) {
        if ((tid & 1u) == 0u && inside)
            packed[idx >> 1] = _oq_pack[tid] | (_oq_pack[tid + 1] << 16);
    } else {  // BITS == 8
        if ((tid & 3u) == 0u && inside)
            packed[idx >> 2] = _oq_pack[tid] | (_oq_pack[tid + 1] << 8) |
                               (_oq_pack[tid + 2] << 16) |
                               (_oq_pack[tid + 3] << 24);
    }
    GroupMemoryBarrierWithGroupSync();
}

// ---------------------------------------------------------------------------
// Masked sub-word atomic stores. For packed layouts where one thread owns a
// CONTIGUOUS cell run that need not be word-aligned (per-splat runs in the
// grad-quant / FPBO writebacks), neighbor threads can share a word. And+Or
// with disjoint masks is correct under any interleaving: each atomic only
// touches the owner's bits, so the final word is
// (old & ~maskA & ~maskB) | bitsA | bitsB regardless of order.
// ---------------------------------------------------------------------------

[ForceInline]
void vkq_store_u16(uint* packed, uint idx, uint h) {  // idx in halfwords
    uint shift = 16u * (idx & 1u);
    InterlockedAnd(packed[idx >> 1], ~(0xffffu << shift));
    InterlockedOr(packed[idx >> 1], (h & 0xffffu) << shift);
}

[ForceInline]
void vkq_store_u8(uint* packed, uint idx, uint b) {  // idx in bytes
    uint shift = 8u * (idx & 3u);
    InterlockedAnd(packed[idx >> 2], ~(0xffu << shift));
    InterlockedOr(packed[idx >> 2], (b & 0xffu) << shift);
}

// ---------------------------------------------------------------------------
// gradq::Codec<BITS> — signed symmetric gradient codec (int8/int16 cells,
// codes in [-qmax, qmax] scaled by the block amplitude; code 0 decodes to
// exactly 0). Mirrors csrc/GradQuant.cuh.
// ---------------------------------------------------------------------------

[ForceInline]
float gradq_amp(float2 mm) { return max(abs(mm.x), abs(mm.y)); }

[ForceInline]
float gradq_qmax<let BITS : int>() { return BITS == 16 ? 32767.0f : 127.0f; }

[ForceInline]
float gradq_decode1<let BITS : int>(uint* packed, uint cell, float2 mm) {
    int s;
    if (BITS == 16) {
        uint w = packed[cell >> 1];
        uint h = (w >> (16u * (cell & 1u))) & 0xffffu;
        s = (int(h) << 16) >> 16;  // sign-extend int16
    } else {
        uint w = packed[cell >> 2];
        uint b = (w >> (8u * (cell & 3u))) & 0xffu;
        s = (int(b) << 24) >> 24;  // sign-extend int8
    }
    return float(s) * (gradq_amp(mm) * (1.0f / gradq_qmax<BITS>()));
}

[ForceInline]
float gradq8_decode1(uint* packed, uint cell, float2 mm) {
    return gradq_decode1<8>(packed, cell, mm);
}

// Round-half-even, matching CUDA rintf in gradq::Codec::encode1.
[ForceInline]
float _gq_rne(float x) {
    float f = floor(x);
    float d = x - f;
    if (d < 0.5f) return f;
    if (d > 0.5f) return f + 1.0f;
    return (f - 2.0f * floor(0.5f * f)) == 0.0f ? f : f + 1.0f;
}

// Encode one value to its signed code bits (caller stores via vkq_store_*).
[ForceInline]
uint gradq_encode1<let BITS : int>(float v, float2 mm) {
    float a = gradq_amp(mm);
    float qmax = gradq_qmax<BITS>();
    float qf = (a > 0.0f) ? _gq_rne(v * (qmax / a)) : 0.0f;
    qf = clamp(qf, -qmax, qmax);
    int q = int(qf);
    return uint(q) & (BITS == 16 ? 0xffffu : 0xffu);
}
