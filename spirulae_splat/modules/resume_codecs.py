"""NumPy host-side mirrors of the block quantization codecs in
src/core/Tensor.h, for adapting a checkpoint's quantized
buffers when resuming with a different layout (fewer Gaussians / different SH
degree). Decode -> transform (gather / band-resample) -> re-encode with fresh
per-block bounds. All host / CPU, vectorized (no per-element Python loops).

Layouts (see Tensor.h):
  QuantizedTensor<BITS>     -- linear, uint packed, float2 bounds (min,max)/block
  QuantizedTensorLog<BITS>  -- log1p(v/eps), uint packed, float2 bounds/block
  QuantizedAdamState<BITS>  -- joint (u, log_s) AoS pair, float4 bounds/block

The per-block "block" spans `bbc` CONTIGUOUS cells. Cell-block layout uses
bbc=256; the per-splat-block (FPBO) layout uses bbc = 256 * cells_per_splat
(a splat's cells are contiguous, so 256 splats == 256*cpp contiguous cells).
"""

import numpy as np

BLOCK = 256
EPS = 1e-15


# --- packed uint <-> value arrays -------------------------------------------

def unpack_uint(packed, n, bits):
    """Single-value packed stream -> uint array [n]."""
    p = np.asarray(packed, dtype=np.uint8)
    if bits == 8:
        return p[:n].astype(np.uint32)
    if bits == 16:
        return p.view(np.uint16)[:n].astype(np.uint32)
    if bits == 4:
        b = p[: (n + 1) // 2].astype(np.uint32)
        q = np.zeros(n, dtype=np.uint32)
        q[0::2] = b & 0x0F
        q[1::2] = (b[: n // 2] >> 4) & 0x0F
        return q
    raise ValueError(f"bits={bits}")


def pack_uint(q, bits):
    q = np.asarray(q, dtype=np.uint32)
    n = q.shape[0]
    if bits == 8:
        return q.astype(np.uint8)
    if bits == 16:
        return q.astype(np.uint16).view(np.uint8)
    if bits == 4:
        out = np.zeros((n + 1) // 2, dtype=np.uint8)
        out |= (q[0::2] & 0x0F).astype(np.uint8)
        out[: n // 2] |= ((q[1::2] & 0x0F) << 4).astype(np.uint8)
        return out
    raise ValueError(f"bits={bits}")


def unpack_pairs(packed, n, bits):
    """AoS (u, log_s) pair stream -> (u_q[n], s_q[n])."""
    p = np.asarray(packed, dtype=np.uint8)
    if bits == 16:
        w = p.view(np.uint16)
        return w[0 : 2 * n : 2].astype(np.uint32), w[1 : 2 * n : 2].astype(np.uint32)
    if bits == 8:
        return p[0 : 2 * n : 2].astype(np.uint32), p[1 : 2 * n : 2].astype(np.uint32)
    if bits == 4:
        b = p[:n].astype(np.uint32)
        return b & 0x0F, (b >> 4) & 0x0F
    raise ValueError(f"bits={bits}")


def pack_pairs(u, s, bits):
    u = np.asarray(u, dtype=np.uint32); s = np.asarray(s, dtype=np.uint32)
    n = u.shape[0]
    if bits == 16:
        out = np.empty(2 * n, dtype=np.uint16)
        out[0::2] = u; out[1::2] = s
        return out.view(np.uint8)
    if bits == 8:
        out = np.empty(2 * n, dtype=np.uint8)
        out[0::2] = u.astype(np.uint8); out[1::2] = s.astype(np.uint8)
        return out
    if bits == 4:
        return ((u & 0x0F) | ((s & 0x0F) << 4)).astype(np.uint8)
    raise ValueError(f"bits={bits}")


# --- per-block min/max over contiguous runs of `bbc` cells ------------------

def _block_minmax(vals, bbc, n):
    idx = np.arange(0, n, bbc)
    return np.minimum.reduceat(vals, idx), np.maximum.reduceat(vals, idx)


def _bidx(n, bbc):
    return np.arange(n, dtype=np.int64) // bbc


# --- QuantizedTensor (linear) -----------------------------------------------

def decode_linear(packed, bounds, n, bits, bbc):
    q = unpack_uint(packed, n, bits).astype(np.float64)
    qmax = float((1 << bits) - 1)
    b = np.asarray(bounds, dtype=np.float64).reshape(-1, 2)
    bi = _bidx(n, bbc)
    lo = b[bi, 0]; hi = b[bi, 1]
    return lo + (hi - lo) * (q / qmax)


def encode_linear(vals, bits, bbc):
    vals = np.asarray(vals, dtype=np.float64)
    n = vals.shape[0]
    lo, hi = _block_minmax(vals, bbc, n)
    qmax = float((1 << bits) - 1)
    bi = _bidx(n, bbc)
    rng = np.maximum(hi[bi] - lo[bi], 1e-30)
    q = np.clip(np.round(qmax * (vals - lo[bi]) / rng), 0, qmax)
    bounds = np.stack([lo, hi], axis=1).astype(np.float32)
    return pack_uint(q, bits), bounds


# --- QuantizedTensorLog (log domain) ----------------------------------------

def decode_log(packed, bounds, n, bits, bbc):
    log_v = decode_linear(packed, bounds, n, bits, bbc)  # bounds live in log domain
    return EPS * np.expm1(log_v)


def encode_log(vals, bits, bbc):
    log_v = np.log1p(np.maximum(np.asarray(vals, np.float64), 0.0) / EPS)
    # same as encode_linear but range floor is kEps (matches Tensor.h encode_log)
    n = log_v.shape[0]
    lo, hi = _block_minmax(log_v, bbc, n)
    qmax = float((1 << bits) - 1)
    bi = _bidx(n, bbc)
    rng = np.maximum(hi[bi] - lo[bi], EPS)
    q = np.clip(np.round(qmax * (log_v - lo[bi]) / rng), 0, qmax)
    return pack_uint(q, bits), np.stack([lo, hi], axis=1).astype(np.float32)


# --- QuantizedAdamState (joint u, log_s) ------------------------------------

def decode_adam(packed, bounds, n, bits, bbc):
    """Returns (g1, g2) arrays of length n."""
    u_q, s_q = unpack_pairs(packed, n, bits)
    qmax = float((1 << bits) - 1)
    b = np.asarray(bounds, dtype=np.float64).reshape(-1, 4)
    bi = _bidx(n, bbc)
    u = b[bi, 0] + (b[bi, 1] - b[bi, 0]) * (u_q.astype(np.float64) / qmax)
    log_s = b[bi, 2] + (b[bi, 3] - b[bi, 2]) * (s_q.astype(np.float64) / qmax)
    sqrt_g2 = EPS * np.expm1(log_s)
    g1 = u * (sqrt_g2 + EPS)
    g2 = sqrt_g2 * sqrt_g2
    return g1, g2


def encode_adam(g1, g2, bits, bbc):
    g1 = np.asarray(g1, np.float64); g2 = np.asarray(g2, np.float64)
    n = g1.shape[0]
    sqrt_g2 = np.sqrt(np.maximum(g2, 0.0))
    u = g1 / (sqrt_g2 + EPS)
    log_s = np.log1p(sqrt_g2 / EPS)
    umin, umax = _block_minmax(u, bbc, n)
    smin, smax = _block_minmax(log_s, bbc, n)
    qmax = float((1 << bits) - 1)
    bi = _bidx(n, bbc)
    u_rng = np.maximum(umax[bi] - umin[bi], EPS)
    s_rng = np.maximum(smax[bi] - smin[bi], EPS)
    u_q = np.clip(np.round(qmax * (u - umin[bi]) / u_rng), 0, qmax)
    s_q = np.clip(np.round(qmax * (log_s - smin[bi]) / s_rng), 0, qmax)
    bounds = np.stack([umin, umax, smin, smax], axis=1).astype(np.float32)
    return pack_pairs(u_q, s_q, bits), bounds
