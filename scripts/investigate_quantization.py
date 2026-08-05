#!/usr/bin/env python3
"""Investigate optimizer-state quantization for a Spirula Studio checkpoint.

Loads g1_*/g2_*/bilagrid_*/ppisp_* tensors from <ckpt>/full/ (the layout
written by engine_save_checkpoint(full_dump=True)), mirrors the current
on-device 8-bit + per-256-cell-block quantization (sqrt transform for g2),
and evaluates ALTERNATIVE block-ordering strategies to test the hypothesis
that grouping similar-magnitude parameters into the same block reduces
quantization error.

Strategies explored, per tensor:
  - identity (memory order: matches the current kernel)
  - all axis permutations (or a curated set for high-rank tensors)
  - sort-by-|x| within super-blocks of various sizes (4K, 64K)
  - sort-by-|x| globally  (oracle lower bound; not practical to deploy)

Quantitative output per tensor x strategy:
  - RMSE, MAE, max-abs error, dequant SNR (dB), mean / p90 / max block range
  - "wasted-precision" indicator (range CV across blocks)

Visualizations:
  - linear + log histograms of values and of |x|
  - rank-sorted |x| log-log
  - heatmap of block-range under each strategy (cumulative density)
  - SH g_features_sh: per-coefficient mean magnitude bar (K x 3 axes)
  - bilagrid g_rgb/depth/normal: per-channel mean magnitude, spatial L-tile
    mosaic (mean over cameras), camera magnitude variance
  - g1 vs g2 joint scatter (log-log) — Adam moment correlation
  - per-strategy block-range histogram overlay
  - SNR bar chart across strategies

The full report is written as a single self-contained HTML file with
embedded base64 PNGs (no external assets), easy to scp / view offline.
"""

from __future__ import annotations

import argparse
import base64
import io
import itertools
import os
import re
import sys
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Callable

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


BLOCK_SIZE = 256                # cells per quantization block (matches kernel)
QUANT_LEVELS = 256              # 8 bits/cell
SAMPLE_CAP = 4_000_000          # cap for histogram inputs (huge tensors)
SUPER_BLOCK_SIZES = (4096, 65536)  # for sort-within-superblock strategies


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _b64_fig(fig, dpi: int = 110) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _img(fig, alt: str = "") -> str:
    return f'<img alt="{alt}" src="data:image/png;base64,{_b64_fig(fig)}" />'


def _fmt(x: float, n: int = 4) -> str:
    if x == 0:
        return "0"
    if not np.isfinite(x):
        return str(x)
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        return f"{x:.{n}g}"
    return f"{x:.{n}g}"


def _subsample(x: np.ndarray, k: int = SAMPLE_CAP, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng(0)
    if x.size <= k:
        return x
    idx = rng.choice(x.size, size=k, replace=False)
    return x.reshape(-1)[idx]


# ---------------------------------------------------------------------------
# Field discovery
# ---------------------------------------------------------------------------

@dataclass
class Field:
    name: str                      # filename minus extension, e.g. "g1_means"
    kind: str                      # "g1" | "g2" | "param" | "bilagrid_g1" | "bilagrid_g2" | ...
    path: Path
    is_quantized_on_disk: bool     # checkpoint stored as q + bounds (lossy already)
    bounds_path: Optional[Path] = None
    array: Optional[np.ndarray] = None
    shape: tuple = ()
    dtype: str = ""


def _maybe_reshape_buffer(name: str, arr: np.ndarray, meta: dict) -> np.ndarray:
    """When the on-disk array is a flat 1D float (the case after decoding a
    packed AoS bilagrid or SH momentum buffer), re-impose the buffer's
    natural multi-dim shape using metadata so structural plots and axis
    permutations apply."""
    if arr.ndim != 1:
        return arr
    # SH momentum: shape (N_splats, K, 3) — channel last.
    if name.startswith("g1_features_sh") or name.startswith("g2_features_sh"):
        try:
            K = int(meta["num_sh"])
            N = int(meta.get("cur_num_splats", meta.get("max_num_splats", 0)))
            if K > 0 and arr.size % (K * 3) == 0:
                N = arr.size // (K * 3)
                return arr.reshape(N, K, 3)
        except (KeyError, ValueError):
            pass
        return arr
    # Bilagrid grids: shape (N_cam, L, H, W, C).
    if name.startswith("bilagrid_rgb"):
        prefix = "bilagrid_rgb"
    elif name.startswith("bilagrid_depth"):
        prefix = "bilagrid_depth"
    elif name.startswith("bilagrid_normal"):
        prefix = "bilagrid_normal"
    else:
        return arr
    try:
        L = int(meta[f"{prefix}_L"]); H = int(meta[f"{prefix}_H"]); W = int(meta[f"{prefix}_W"])
        if prefix == "bilagrid_rgb":
            C = int(meta.get("bilagrid_rgb_C", 0)) or (9 if meta.get("bilagrid_rgb_type") in ("ppisp", "loglinear") else 12)
        elif prefix == "bilagrid_depth":
            C = 2
        else:
            C = 3
        per_cam = L * H * W * C
        if arr.size % per_cam != 0:
            return arr
        N = arr.size // per_cam
        return arr.reshape(N, L, H, W, C)
    except (KeyError, ValueError):
        return arr


# Backwards-compat alias kept for callers that referenced the older name.
_maybe_reshape_bilagrid = _maybe_reshape_buffer


# ---------------------------------------------------------------------------
# Joint AoS packed dequantizer.
#
# Two codecs share the same on-disk layout (one AoS uint8 stream + one float4
# bounds array per 256-cell block); the meta.txt `quant_codec` tag picks the
# inverse map for the second half:
#
#   joint_u_sqrt_g2_v0   — legacy: both halves linearly quantized in raw units.
#                          u       = u_lo + (u_hi - u_lo) * u_q / 255
#                          sqrt_g2 = s_lo + (s_hi - s_lo) * s_q / 255
#
#   joint_u_log_sqrt_g2_v1 — current: u linear, sqrt_g2 = eps · expm1(log_s).
#                          The packed `s_q` byte encodes `log_s` linearly within
#                          the (log-space) per-block bounds. This is what
#                          QuantizedAdamState's forward_sqrt_g2 / inverse_sqrt_g2
#                          do on-device. Buys ~20–30 dB of post-step SNR(u_next)
#                          over the linear variant at no storage cost — see the
#                          per-pair table in the report's joint-pair section.
#
# In both cases the reconstructed Adam state is
#   g2 = sqrt_g2 ** 2
#   g1 = u * (sqrt_g2 + eps)
# matching the on-device decoder.
# ---------------------------------------------------------------------------

JOINT_AOS_EPS = 1e-15

QUANT_CODEC_LINEAR = "joint_u_sqrt_g2_v0"          # legacy linear sqrt_g2
QUANT_CODEC_LOG    = "joint_u_log_sqrt_g2_v1"      # current log_pos sqrt_g2


def _dequant_packed_aos(packed: np.ndarray, bounds: np.ndarray,
                         eps: float = JOINT_AOS_EPS,
                         codec: str = QUANT_CODEC_LOG,
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Decode a joint AoS-packed uint8 buffer into (g1, g2). `codec` picks
    the inverse map for the sqrt_g2 half."""
    flat_u8 = packed.reshape(-1)
    n = flat_u8.size
    if n % 2 != 0:
        raise ValueError(f"packed size {n} not even — expected interleaved bytes")
    n_cells = n // 2
    nb = (n_cells + BLOCK_SIZE - 1) // BLOCK_SIZE
    if bounds.shape[0] != nb:
        raise ValueError(f"bounds shape {bounds.shape} mismatch with packed cells "
                          f"{n_cells} → expected n_blocks={nb}")
    u_q = flat_u8[0::2].astype(np.float32)
    s_q = flat_u8[1::2].astype(np.float32)
    pad = nb * BLOCK_SIZE - n_cells
    if pad:
        u_q = np.concatenate([u_q, np.zeros(pad, dtype=np.float32)])
        s_q = np.concatenate([s_q, np.zeros(pad, dtype=np.float32)])
    u_q = u_q.reshape(nb, BLOCK_SIZE)
    s_q = s_q.reshape(nb, BLOCK_SIZE)
    u_lo = bounds[:, 0:1]; u_hi = bounds[:, 1:2]
    s_lo = bounds[:, 2:3]; s_hi = bounds[:, 3:4]
    u_val = u_lo + (u_hi - u_lo) * (u_q * (1.0 / 255.0))
    s_val = s_lo + (s_hi - s_lo) * (s_q * (1.0 / 255.0))
    if codec == QUANT_CODEC_LOG:
        # log_s -> sqrt_g2 = eps * expm1(log_s). expm1 in float64 to avoid
        # overflow at large log_s (typical training values land at ~20–35).
        sqrt_g2 = (eps * np.expm1(s_val.astype(np.float64))).astype(np.float32)
    elif codec == QUANT_CODEC_LINEAR:
        sqrt_g2 = s_val
    else:
        raise ValueError(f"unknown quant codec: {codec!r}")
    g2 = sqrt_g2 * sqrt_g2
    g1 = u_val * (sqrt_g2 + eps)
    g1 = g1.reshape(-1)[:n_cells].astype(np.float32)
    g2 = g2.reshape(-1)[:n_cells].astype(np.float32)
    return g1, g2


def discover_fields(full_dir: Path, meta: Optional[dict] = None) -> list[Field]:
    """Pick optimizer-state and parameter buffers worth analyzing.

    Three storage formats are recognised:
      1. Float SoA  : `<name>.npy` (regular per-buffer float arrays).
      2. Old quant  : `<name>_q.npy` (uint8) + sibling `*_quant_bounds.npy`
                      (per-component encoded; legacy from the pre-joint era).
      3. AoS packed : `<base>_packed.npy` (uint8, interleaved u and sqrt_g2)
                      + sibling `*_quant_bounds.npy` or `*_bounds.npy`. One
                      packed file decodes into BOTH g1 and g2 logical fields.
    """
    meta = meta or {}
    candidates = sorted(full_dir.glob("*.npy"))
    fields: list[Field] = []
    by_name = {p.stem: p for p in candidates}

    def looks_like_momentum(stem: str) -> Optional[str]:
        if stem.startswith("g1_"):
            return "g1"
        if stem.startswith("g2_"):
            return "g2"
        if stem.endswith("_g1"):
            return "bilagrid_g1"
        if stem.endswith("_g2"):
            return "bilagrid_g2"
        if stem == "ppisp_g1":
            return "ppisp_g1"
        if stem == "ppisp_g2":
            return "ppisp_g2"
        return None

    def find_bounds(base: str) -> Optional[Path]:
        # Try both naming conventions emitted by engine_save_checkpoint:
        # bilagrid: `<base>_quant_bounds.npy`; SH: `<base>_bounds.npy`.
        for cand in (f"{base}_quant_bounds", f"{base}_bounds"):
            if cand in by_name:
                return by_name[cand]
        return None

    # --- AoS packed files: decode into a (g1, g2) pair ----------------------
    consumed: set[str] = set()
    for stem, path in by_name.items():
        if not stem.endswith("_packed"):
            continue
        base = stem[:-len("_packed")]     # "sh_quant" / "bilagrid_rgb" / ...
        bounds_path = find_bounds(base)
        if bounds_path is None:
            print(f"warn: {stem} has no sibling bounds file; skipping",
                   file=sys.stderr)
            continue
        consumed.add(stem)
        consumed.add(bounds_path.stem)

        if base == "sh_quant":
            g1_name = "g1_features_sh"
            g2_name = "g2_features_sh"
            g1_kind = "g1"
            g2_kind = "g2"
        else:
            g1_name = f"{base}_g1"
            g2_name = f"{base}_g2"
            g1_kind = "bilagrid_g1" if base.startswith("bilagrid_") else "param"
            g2_kind = "bilagrid_g2" if base.startswith("bilagrid_") else "param"

        # We synthesize TWO Field entries per packed file. Both share the
        # same source `.path` (packed bytes) and `.bounds_path`, but `kind`
        # picks which half to surface as the field's `.array`.
        fields.append(Field(name=g1_name, kind=g1_kind, path=path,
                             is_quantized_on_disk=True, bounds_path=bounds_path))
        fields.append(Field(name=g2_name, kind=g2_kind, path=path,
                             is_quantized_on_disk=True, bounds_path=bounds_path))

    # --- Legacy per-component `*_q.npy` files ------------------------------
    for stem, path in by_name.items():
        if stem in consumed or not stem.endswith("_q"):
            continue
        base = stem[:-2]
        kind = looks_like_momentum(base) or "param"
        if base in ("g1_features_sh", "g2_features_sh"):
            bounds_key = "quant_bounds_sh"
        elif base.endswith("_g1") or base.endswith("_g2"):
            bounds_key = f"{base[:-3]}_quant_bounds"
        else:
            bounds_key = f"{base}_quant_bounds"
        bounds_path = by_name.get(bounds_key)
        consumed.add(stem)
        if bounds_path is not None:
            consumed.add(bounds_path.stem)
        if bounds_path is None:
            print(f"warn: {stem} has no bounds file (expected "
                   f"'{bounds_key}.npy'); loading raw uint8", file=sys.stderr)
        fields.append(Field(name=base, kind=kind, path=path,
                             is_quantized_on_disk=True,
                             bounds_path=bounds_path))

    # --- Float SoA `<name>.npy` files --------------------------------------
    for stem, path in by_name.items():
        if stem in consumed:
            continue
        if looks_like_momentum(stem):
            fields.append(Field(name=stem, kind=looks_like_momentum(stem),
                                 path=path, is_quantized_on_disk=False))
        elif stem.startswith("world_") or stem in ("accum_buffer", "bias_correction_steps"):
            fields.append(Field(name=stem, kind="param", path=path,
                                 is_quantized_on_disk=False))

    # --- Load arrays ------------------------------------------------------
    # Cache decoded (g1, g2) per packed source path so we only decode once.
    decoded_cache: dict[Path, tuple[np.ndarray, np.ndarray]] = {}
    # Old checkpoints (saved before the codec tag was added) lacked the
    # `quant_codec` meta line — those used the legacy linear sqrt_g2 mapping,
    # so we fall back to that when the tag is absent.
    codec_tag = meta.get("quant_codec", QUANT_CODEC_LINEAR)

    def _decode_packed(packed_path: Path, bounds_path: Path
                        ) -> tuple[np.ndarray, np.ndarray]:
        if packed_path in decoded_cache:
            return decoded_cache[packed_path]
        packed = np.load(packed_path)
        bounds = np.load(bounds_path)
        if packed.dtype != np.uint8:
            raise ValueError(f"{packed_path.name}: expected uint8, got {packed.dtype}")
        if bounds.dtype != np.float32:
            bounds = bounds.astype(np.float32)
        if bounds.ndim == 1 and bounds.size % 4 == 0:
            bounds = bounds.reshape(-1, 4)
        if bounds.ndim != 2 or bounds.shape[1] != 4:
            raise ValueError(f"{bounds_path.name}: expected shape [n_blocks, 4], got {bounds.shape}")
        g1, g2 = _dequant_packed_aos(packed, bounds, codec=codec_tag)
        decoded_cache[packed_path] = (g1, g2)
        return g1, g2

    for f in fields:
        if f.is_quantized_on_disk and f.bounds_path is not None and f.path.name.endswith("_packed.npy"):
            g1, g2 = _decode_packed(f.path, f.bounds_path)
            # Choose the right half by the field's logical role.
            half = g1 if (f.kind.startswith(("g1", "bilagrid_g1", "ppisp_g1"))
                          or f.name.startswith("g1_")
                          or f.name.endswith("_g1")) else g2
            f.array = _maybe_reshape_buffer(f.name, half, meta)
        else:
            a = np.load(f.path)
            if f.is_quantized_on_disk and f.bounds_path is not None:
                # Legacy per-component path: SoA uint8 + bounds.
                bounds = np.load(f.bounds_path)
                a = _dequant_uint8_blocks(a, bounds, kind=f.kind)
                a = _maybe_reshape_buffer(f.name, a, meta)
            f.array = a
        f.shape = tuple(f.array.shape)
        f.dtype = str(f.array.dtype)

    fields.sort(key=lambda f: (0 if f.kind.startswith(("g1", "bilagrid_g1", "ppisp_g1"))
                                else 1 if f.kind.startswith(("g2", "bilagrid_g2", "ppisp_g2"))
                                else 2, f.name))
    return fields


def _dequant_uint8_blocks(q: np.ndarray, bounds: np.ndarray, kind: str) -> np.ndarray:
    """Dequant a [..., 3?] uint8 array stored in 256-cell blocks (memory order)
    using float4 bounds = (g1lo, g1hi, sqrtg2lo, sqrtg2hi). Used to recover the
    on-disk-quantized SH momentum so the rest of the script can treat it as
    float. This is intrinsically lossy — we add a note in the report."""
    flat_u8 = q.reshape(-1).astype(np.float32)
    n = flat_u8.size
    nb = (n + BLOCK_SIZE - 1) // BLOCK_SIZE
    assert bounds.shape[0] == nb, (q.shape, bounds.shape)
    # Pick the (lo, hi) channel based on kind.
    if kind in ("g1", "bilagrid_g1", "ppisp_g1") or "g1_" in kind:
        lo = bounds[:, 0]
        hi = bounds[:, 1]
        apply_sq = False
    else:
        lo = bounds[:, 2]
        hi = bounds[:, 3]
        apply_sq = True

    out = np.empty(nb * BLOCK_SIZE, dtype=np.float32)
    out[:n] = ((flat_u8 + 0.5) / 256.0)
    # broadcast lo/hi per-block:
    out_blocks = out.reshape(nb, BLOCK_SIZE)
    out_blocks = lo[:, None] + (hi[:, None] - lo[:, None]) * out_blocks
    if apply_sq:
        out_blocks = out_blocks * out_blocks
    return out_blocks.reshape(-1)[:n].reshape(q.shape).astype(np.float32)


# ---------------------------------------------------------------------------
# Quantization simulation under arbitrary block orderings
# ---------------------------------------------------------------------------

def _pad_to_blocks(flat: np.ndarray) -> np.ndarray:
    n = flat.size
    nb = (n + BLOCK_SIZE - 1) // BLOCK_SIZE
    pad = nb * BLOCK_SIZE - n
    if pad:
        flat = np.concatenate([flat, np.zeros(pad, dtype=flat.dtype)])
    return flat.reshape(nb, BLOCK_SIZE)


# ---------------------------------------------------------------------------
# Value mappings applied before the per-block min/max + uint8 quantize step.
# A mapping is a (forward, inverse) pair acting elementwise on a float array.
# We map -> quantize linearly in the mapped domain -> dequantize -> map^{-1}.
#
# Choices:
#   linear      : identity. Good when in-block dynamic range is moderate.
#   sqrt        : x >= 0 -> sqrt(x). Matches the on-device g2 path: spreads
#                 [0, max] more uniformly when distribution is gradient-squared.
#   signed_log  : preserves sign; 0 -> 0; |x| -> log1p(|x|/eps).
#                 Useful for g1 buffers whose magnitude spans many orders of
#                 magnitude — concentrates resolution near zero and gracefully
#                 handles values close to zero (where Adam still cares about
#                 sign/direction). eps controls the resolution near 0.
#   log_pos     : non-negative input only. x -> log(max(x,0) + eps). Useful
#                 for sqrt(g2)-style or g2 buffers, similar motivation to
#                 signed_log without the sign-preservation cost.
#   log_sqrt    : composite for g2: y = log(sqrt(g2) + eps). Halves the log
#                 dynamic range vs log_pos applied to g2 directly.
# ---------------------------------------------------------------------------

VALUE_MAPPING_EPS = 1e-15        # matches the Adam denominator epsilon


def _map_forward(x: np.ndarray, mapping: str) -> np.ndarray:
    x = x.astype(np.float64, copy=False)
    if mapping == "linear":
        return x
    if mapping == "sqrt":
        return np.sqrt(np.maximum(x, 0.0))
    if mapping == "signed_log":
        return np.sign(x) * np.log1p(np.abs(x) / VALUE_MAPPING_EPS)
    if mapping == "log_pos":
        return np.log(np.maximum(x, 0.0) + VALUE_MAPPING_EPS)
    if mapping == "log_sqrt":
        return np.log(np.sqrt(np.maximum(x, 0.0)) + VALUE_MAPPING_EPS)
    raise ValueError(f"unknown mapping: {mapping}")


def _map_inverse(y: np.ndarray, mapping: str) -> np.ndarray:
    if mapping == "linear":
        return y
    if mapping == "sqrt":
        return y * y
    if mapping == "signed_log":
        return np.sign(y) * VALUE_MAPPING_EPS * np.expm1(np.abs(y))
    if mapping == "log_pos":
        return np.maximum(np.exp(y) - VALUE_MAPPING_EPS, 0.0)
    if mapping == "log_sqrt":
        s = np.maximum(np.exp(y) - VALUE_MAPPING_EPS, 0.0)
        return s * s
    raise ValueError(f"unknown mapping: {mapping}")


def _block_quantize(blocks: np.ndarray, mapping: str) -> tuple[np.ndarray, np.ndarray]:
    """Per-block 8-bit quantize under the given value mapping.
    Returns (dequantized_blocks_in_raw_space, per_block_mapped_range)."""
    t = _map_forward(blocks, mapping)
    lo = t.min(axis=1, keepdims=True)
    hi = t.max(axis=1, keepdims=True)
    rng = np.maximum(hi - lo, 1e-30)
    q = np.clip(np.floor(QUANT_LEVELS * (t - lo) / rng), 0, QUANT_LEVELS - 1)
    deq_t = lo + rng * ((q + 0.5) / QUANT_LEVELS)
    deq = _map_inverse(deq_t, mapping)
    return deq.astype(np.float32), (hi - lo).squeeze(1).astype(np.float64)


def _default_mapping(kind: str) -> str:
    """The mapping that matches the current on-device kernel for each kind."""
    if kind.endswith("g2") or "_g2" in kind:
        return "sqrt"
    return "linear"


def _is_positive_kind(kind: str) -> bool:
    """Whether values in this buffer are non-negative (so log_pos / log_sqrt
    are admissible)."""
    return kind.endswith("g2") or "_g2" in kind


# ---------------------------------------------------------------------------
# Spatial sorting: precompute per-splat permutations once per checkpoint,
# then reuse for every per-splat field (means/scales/quats/opacities/dc/sh
# and their g1/g2 momentum). Hypothesis: spatially nearby splats optimize
# similar features, so their momentum magnitudes are correlated; grouping
# them into the same 256-cell block tightens the block range.
# ---------------------------------------------------------------------------

NBITS_SPATIAL = 10               # grid size per axis = 1024 (3D = 1B cells)


def _norm_bbox(coords: np.ndarray) -> np.ndarray:
    """Min-max per axis -> [0, 1]^3."""
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    rng = np.maximum(hi - lo, 1e-30)
    return ((coords - lo) / rng).clip(0.0, 1.0)


def _norm_quantile(coords: np.ndarray) -> np.ndarray:
    """Independent per-axis rank/N -> [0, 1]^3 (true CDF normalization)."""
    out = np.empty(coords.shape, dtype=np.float64)
    n = coords.shape[0]
    for d in range(coords.shape[1]):
        # argsort(argsort(x)) gives the rank of each element.
        ranks = np.empty(n, dtype=np.int64)
        ranks[np.argsort(coords[:, d], kind="stable")] = np.arange(n)
        out[:, d] = ranks / max(n - 1, 1)
    return out


def _norm_moment_std(coords: np.ndarray, k: float) -> np.ndarray:
    """Per-axis (x - mean) / (k * std), shifted to [0, 1] with clip. Higher k
    keeps more tail inside [0, 1] (less clipping, lower spatial resolution)."""
    mu = coords.mean(axis=0)
    sigma = coords.std(axis=0) + 1e-30
    return np.clip((coords - mu) / (k * sigma) + 0.5, 0.0, 1.0)


def _norm_moment_mad(coords: np.ndarray, k: float) -> np.ndarray:
    """Robust variant: median + MAD (scaled so MAD == sigma for Gaussian)."""
    med = np.median(coords, axis=0)
    mad = np.median(np.abs(coords - med), axis=0) + 1e-30
    return np.clip((coords - med) / (k * mad * 1.4826) + 0.5, 0.0, 1.0)


# --- Morton (z-curve) 3D ---------------------------------------------------

def _morton_expand_3d(x: np.ndarray) -> np.ndarray:
    """Bit-expand low 21 bits of x into positions 0, 3, 6, ... of a uint64."""
    x = x.astype(np.uint64) & np.uint64(0x1fffff)
    x = (x | (x << np.uint64(32))) & np.uint64(0x1f00000000ffff)
    x = (x | (x << np.uint64(16))) & np.uint64(0x1f0000ff0000ff)
    x = (x | (x << np.uint64(8)))  & np.uint64(0x100f00f00f00f00f)
    x = (x | (x << np.uint64(4)))  & np.uint64(0x10c30c30c30c30c3)
    x = (x | (x << np.uint64(2)))  & np.uint64(0x1249249249249249)
    return x


def _morton_3d(coords_int: np.ndarray) -> np.ndarray:
    return (_morton_expand_3d(coords_int[:, 0])
            | (_morton_expand_3d(coords_int[:, 1]) << np.uint64(1))
            | (_morton_expand_3d(coords_int[:, 2]) << np.uint64(2)))


# --- Hilbert curve 3D (vectorized Skilling AxestoTranspose) ---------------

def _hilbert_axes_to_transpose(X: np.ndarray, nbits: int) -> None:
    """In-place; X is (N, 3) uint64."""
    ndim = X.shape[1]
    M = np.uint64(1) << np.uint64(nbits - 1)
    Q = M
    while Q > np.uint64(1):
        P = Q - np.uint64(1)
        for i in range(ndim):
            mask = (X[:, i] & Q) != 0
            X[mask, 0] ^= P
            # swap low-P bits of X[0] and X[i] where mask is False
            not_mask = ~mask
            t = (X[:, 0] ^ X[:, i]) & P
            X[not_mask, 0] ^= t[not_mask]
            X[not_mask, i] ^= t[not_mask]
        Q >>= np.uint64(1)
    # Gray-encode
    for i in range(1, ndim):
        X[:, i] ^= X[:, i - 1]
    t = np.zeros(X.shape[0], dtype=np.uint64)
    Q = M
    while Q > np.uint64(1):
        mask = (X[:, ndim - 1] & Q) != 0
        t[mask] ^= (Q - np.uint64(1))
        Q >>= np.uint64(1)
    for i in range(ndim):
        X[:, i] ^= t


def _hilbert_3d(coords_int: np.ndarray, nbits: int) -> np.ndarray:
    X = coords_int.astype(np.uint64).copy()
    _hilbert_axes_to_transpose(X, nbits)
    # Bit-interleave the transposed coordinates -> linear Hilbert position.
    out = np.zeros(coords_int.shape[0], dtype=np.uint64)
    for bit in range(nbits):
        for d in range(coords_int.shape[1]):
            b = (X[:, d] >> np.uint64(bit)) & np.uint64(1)
            out |= b << np.uint64(bit * coords_int.shape[1] + d)
    return out


NORMALIZATIONS: dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "bbox": _norm_bbox,
    "quantile": _norm_quantile,
    "std_3sig": lambda c: _norm_moment_std(c, 3.0),
    "std_5sig": lambda c: _norm_moment_std(c, 5.0),
    "mad_3sig": lambda c: _norm_moment_mad(c, 3.0),
}


def build_spatial_permutations(world_means: np.ndarray, n_active: int,
                               nbits: int = NBITS_SPATIAL,
                               curves: tuple[str, ...] = ("morton", "hilbert"),
                               norms: tuple[str, ...] = tuple(NORMALIZATIONS),
                              ) -> dict[str, np.ndarray]:
    """Return {strategy_name -> (n_active,) int64 permutation along axis 0}.

    n_active is the *current* number of live splats (means table can be
    over-allocated to max_num_splats); we only sort the first n_active rows
    and apply the resulting permutation downstream."""
    out: dict[str, np.ndarray] = {}
    coords_all = world_means[:n_active]
    levels = (1 << nbits) - 1
    for norm_name in norms:
        norm_fn = NORMALIZATIONS[norm_name]
        norm_coords = norm_fn(coords_all)
        coords_int = np.minimum(np.floor(norm_coords * levels), levels).astype(np.uint32)
        for curve in curves:
            if curve == "morton":
                codes = _morton_3d(coords_int)
            elif curve == "hilbert":
                codes = _hilbert_3d(coords_int, nbits)
            else:
                continue
            order = np.argsort(codes, kind="stable")
            out[f"{curve}_{norm_name}"] = order.astype(np.int64)
    return out


# ---------------------------------------------------------------------------
# Strategy descriptors and execution
# ---------------------------------------------------------------------------

@dataclass
class Strategy:
    """A reordering + value-mapping recipe applied to a tensor before the
    8-bit / 256-cell block quantizer. Components compose in this order:
        1. axis0_perm   (per-splat reorder along axis 0)
        2. axis_perm    (transpose of all axes)
        3. sort_within  (sort by |x| within each super-block)
        4. global_sort  (sort all values by |x|; oracle)
    Then the value mapping (linear/sqrt/signed_log/log_pos/log_sqrt) is
    applied per-block before the uint8 quantize.
    """
    name: str
    extra: str = ""
    axis0_perm: Optional[np.ndarray] = None      # (N,) int permutation of axis 0
    axis_perm: Optional[tuple] = None            # axes permutation tuple
    sort_within: Optional[int] = None            # super-block size
    global_sort: bool = False
    mapping: Optional[str] = None                # None = default for kind


@dataclass
class StrategyResult:
    name: str
    rmse: float
    mae: float
    maxerr: float
    snr_db: float
    block_range_mean: float
    block_range_p90: float
    block_range_max: float
    block_range_cv: float       # std(range) / mean(range)
    block_ranges: np.ndarray = field(default_factory=lambda: np.zeros(0))
    extra: str = ""


def _error_metrics(orig: np.ndarray, recon: np.ndarray) -> tuple[float, float, float, float]:
    o = orig.astype(np.float64)
    r = recon.astype(np.float64)
    d = o - r
    rmse = float(np.sqrt(np.mean(d * d)))
    mae = float(np.mean(np.abs(d)))
    maxerr = float(np.max(np.abs(d)))
    orig_rms = float(np.sqrt(np.mean(o * o)))
    snr = 20.0 * np.log10(max(orig_rms, 1e-30) / max(rmse, 1e-30))
    return rmse, mae, maxerr, snr


def _build_orderings(shape: tuple) -> list[tuple[str, Optional[tuple]]]:
    """Return (name, axis_permutation_or_None) for the strategies that just
    rearrange axes. Sort-based strategies are handled outside."""
    nd = len(shape)
    out = [("identity", tuple(range(nd)))]
    perms = list(itertools.permutations(range(nd)))
    if nd <= 4 and len(perms) <= 24:
        for p in perms:
            if p == tuple(range(nd)):
                continue
            label = f"perm({','.join(str(i) for i in p)})"
            out.append((label, p))
    else:
        # 5D bilagrid (N_cam, L, H, W, C): pick a few semantically meaningful
        # orderings. Identity is already channel-last; these explore moving
        # the channel axis around.
        labels_perms = []
        if nd == 5:
            labels_perms = [
                ("perm(N,C,L,H,W)", (0, 4, 1, 2, 3)),  # legacy channel-first
                ("perm(C,N,L,H,W)", (4, 0, 1, 2, 3)),  # channel outermost
                ("perm(L,H,W,N,C)", (1, 2, 3, 0, 4)),  # spatial outermost
                ("perm(C,L,H,W,N)", (4, 1, 2, 3, 0)),  # cam innermost-but-channel
            ]
        for lab, p in labels_perms:
            out.append((lab, p))
    return out


def _run_strategy(arr: np.ndarray, kind: str, s: Strategy) -> StrategyResult:
    """Apply Strategy s as a reordering + mapping of arr, run the quantizer,
    invert, measure error against the original arr in the raw value space."""
    a = arr
    inv_axis0 = None
    inv_axis = None
    mapping = s.mapping or _default_mapping(kind)

    if s.axis0_perm is not None:
        p = s.axis0_perm
        if p.shape[0] != a.shape[0]:
            p_full = np.arange(a.shape[0], dtype=np.int64)
            p_full[:p.shape[0]] = p
            p = p_full
        a = a[p]
        inv_axis0 = np.argsort(p, kind="stable")

    if s.axis_perm is not None:
        a = np.transpose(a, s.axis_perm)
        inv_axis = np.argsort(s.axis_perm)

    a = np.ascontiguousarray(a)
    flat = a.reshape(-1).astype(np.float32)
    n = flat.size

    # Sort key uses the mapped magnitude — for sqrt/log_pos this matches the
    # on-device kernel (which sorts after sqrt); for signed_log it groups
    # values by signed log-magnitude.
    def _sort_key(x: np.ndarray) -> np.ndarray:
        if mapping == "sqrt":
            return np.sqrt(np.maximum(x, 0.0))
        if mapping in ("log_pos", "log_sqrt"):
            return np.log(np.maximum(x, 0.0) + VALUE_MAPPING_EPS)
        if mapping == "signed_log":
            return np.sign(x) * np.log1p(np.abs(x) / VALUE_MAPPING_EPS)
        return x

    if s.global_sort:
        order = np.argsort(_sort_key(flat), kind="stable")
        flat_re = flat[order]
        blocks = _pad_to_blocks(flat_re)
        deq, ranges = _block_quantize(blocks, mapping)
        deq_flat = deq.reshape(-1)[:n]
        deq_back = np.empty_like(flat)
        deq_back[order] = deq_flat
    elif s.sort_within is not None:
        sb = s.sort_within
        key = _sort_key(flat)
        order = np.empty(n, dtype=np.int64)
        nsb = (n + sb - 1) // sb
        for i in range(nsb):
            lo, hi = i * sb, min(i * sb + sb, n)
            order[lo:hi] = lo + np.argsort(key[lo:hi], kind="stable")
        flat_re = flat[order]
        blocks = _pad_to_blocks(flat_re)
        deq, ranges = _block_quantize(blocks, mapping)
        deq_flat = deq.reshape(-1)[:n]
        deq_back = np.empty_like(flat)
        deq_back[order] = deq_flat
    else:
        blocks = _pad_to_blocks(flat)
        deq, ranges = _block_quantize(blocks, mapping)
        deq_back = deq.reshape(-1)[:n]

    out = deq_back.reshape(a.shape)
    if inv_axis is not None:
        out = np.ascontiguousarray(out.transpose(inv_axis))
    if inv_axis0 is not None:
        out = out[inv_axis0]

    rmse, mae, maxerr, snr = _error_metrics(arr, out)
    rng_mean = float(ranges.mean()) if ranges.size else 0.0
    return StrategyResult(
        name=s.name, rmse=rmse, mae=mae, maxerr=maxerr, snr_db=snr,
        block_range_mean=rng_mean,
        block_range_p90=float(np.percentile(ranges, 90)) if ranges.size else 0.0,
        block_range_max=float(ranges.max()) if ranges.size else 0.0,
        block_range_cv=float(ranges.std() / max(rng_mean, 1e-30)),
        block_ranges=ranges,
        extra=s.extra,
    )


def build_strategies(arr_shape: tuple, n_splats: Optional[int], kind: str,
                     spatial_perms: dict[str, np.ndarray]) -> list[Strategy]:
    """Enumerate every reordering + value-mapping strategy to evaluate."""
    strats: list[Strategy] = []
    default_map = _default_mapping(kind)
    positive = _is_positive_kind(kind)

    # --- 1. identity + axis permutations (default mapping) ---
    for label, perm in _build_orderings(arr_shape):
        if perm is None:
            continue
        if perm == tuple(range(len(arr_shape))):
            strats.append(Strategy(name="identity"))
        else:
            strats.append(Strategy(name=label, axis_perm=perm))

    # --- 2. sort within super-block ---
    n = int(np.prod(arr_shape))
    for sb in SUPER_BLOCK_SIZES:
        if sb < n:
            strats.append(Strategy(name=f"sort_within_{sb}", sort_within=sb,
                                    extra="needs per-cell index within super-block"))

    # --- 3. global oracle sort ---
    strats.append(Strategy(name="oracle_sorted_abs", global_sort=True,
                            extra="lower bound; needs full per-cell permutation at decode"))

    # --- 4. spatial sort along axis 0 (per-splat fields only) ---
    is_per_splat = (n_splats is not None
                    and len(arr_shape) >= 1
                    and arr_shape[0] == n_splats
                    and arr_shape[0] > 1)
    if is_per_splat:
        for sp_name, sp_perm in spatial_perms.items():
            strats.append(Strategy(name=f"spatial_{sp_name}",
                                    axis0_perm=sp_perm,
                                    extra="per-splat reorder; "
                                          "decode cost = 1 perm table shared by all per-splat fields"))
        # Combine the best curve+norm picks with an inner-axis swap (SH: K<->3).
        if len(arr_shape) == 3 and "hilbert_quantile" in spatial_perms:
            strats.append(Strategy(
                name="spatial_hilbert_quantile + perm(0,2,1)",
                axis0_perm=spatial_perms["hilbert_quantile"],
                axis_perm=(0, 2, 1),
                extra="combine spatial sort with inner-axis swap"))
        if len(arr_shape) == 3 and "morton_quantile" in spatial_perms:
            strats.append(Strategy(
                name="spatial_morton_quantile + perm(0,2,1)",
                axis0_perm=spatial_perms["morton_quantile"],
                axis_perm=(0, 2, 1)))

    # --- 5. value-mapping variants (decoder cost: same; bytes/cell: same) ---
    # signed_log is admissible for signed buffers (g1, params); log_pos /
    # log_sqrt require non-negative inputs (g2 only). The mapping replaces
    # the default sqrt-or-linear in the per-block min/max + uint8 step.
    mapping_variants = []
    if positive:
        # g2 default is sqrt; compare with two log-style alternatives.
        mapping_variants.extend(["log_pos", "log_sqrt"])
    else:
        # g1 default is linear; signed_log compresses long tails.
        mapping_variants.append("signed_log")

    for m in mapping_variants:
        if m == default_map:
            continue
        strats.append(Strategy(
            name=f"identity_{m}",
            mapping=m,
            extra=f"value mapping: {m} (default = {default_map})"))
        if SUPER_BLOCK_SIZES and SUPER_BLOCK_SIZES[-1] < n:
            sb = SUPER_BLOCK_SIZES[-1]
            strats.append(Strategy(
                name=f"sort_within_{sb}_{m}",
                sort_within=sb,
                mapping=m,
                extra=f"sort + value mapping: {m}"))
    return strats


def simulate_all_strategies(arr: np.ndarray, kind: str,
                             n_splats: Optional[int] = None,
                             spatial_perms: Optional[dict[str, np.ndarray]] = None,
                             max_workers: int = 1) -> list[StrategyResult]:
    """Run all strategies. The reported RMSE is in raw value space (after
    mapping inverse). Parallel across strategies via ThreadPoolExecutor."""
    strategies = build_strategies(arr.shape, n_splats, kind, spatial_perms or {})
    if max_workers > 1 and len(strategies) > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            fut_map = {ex.submit(_run_strategy, arr, kind, s): i
                        for i, s in enumerate(strategies)}
            results: list[Optional[StrategyResult]] = [None] * len(strategies)
            for f in as_completed(fut_map):
                results[fut_map[f]] = f.result()
        return [r for r in results if r is not None]
    else:
        return [_run_strategy(arr, kind, s) for s in strategies]


# ---------------------------------------------------------------------------
# Plotting helpers (used by per-field section)
# ---------------------------------------------------------------------------

def plot_value_histogram(arr: np.ndarray, title: str) -> str:
    sample = _subsample(arr.reshape(-1))
    finite = sample[np.isfinite(sample)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.5))

    axes[0].hist(finite, bins=200, color="C0")
    axes[0].set_title("values (linear y)")
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)

    abs_finite = np.abs(finite)
    abs_finite = abs_finite[abs_finite > 0]
    if abs_finite.size:
        axes[1].hist(np.log10(abs_finite), bins=200, color="C1")
        axes[1].set_title("log10(|x|)  (non-zero)")
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)

    # Rank-sorted |x| plot (log-log).
    s = np.sort(abs_finite)[::-1]
    ranks = np.arange(1, s.size + 1)
    axes[2].loglog(ranks, s, color="C2")
    axes[2].set_title("rank-sorted |x|")
    axes[2].set_xlabel("rank")
    axes[2].grid(True, which="both", alpha=0.3)

    fig.suptitle(title)
    return _img(fig, alt="value distribution")


def plot_block_ranges(results: list[StrategyResult], title: str) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(14, 3.8))
    for r in results:
        rng = r.block_ranges
        rng = rng[rng > 0]
        if rng.size:
            axes[0].hist(np.log10(rng), bins=80, histtype="step",
                         label=r.name, linewidth=1.2)
    axes[0].set_xlabel("log10(block range)")
    axes[0].set_ylabel("# blocks")
    axes[0].set_yscale("log")
    axes[0].set_title("block-range distribution")
    axes[0].legend(fontsize=7, loc="best")
    axes[0].grid(True, alpha=0.3)

    # Block-range CDF (linear).
    for r in results:
        rng = np.sort(r.block_ranges)
        if rng.size:
            cdf = np.arange(rng.size) / max(rng.size - 1, 1)
            axes[1].plot(rng, cdf, label=r.name, linewidth=1.2)
    axes[1].set_xlabel("block range")
    axes[1].set_ylabel("CDF")
    axes[1].set_xscale("log")
    axes[1].set_title("block-range CDF (smaller is better)")
    axes[1].legend(fontsize=7, loc="best")
    axes[1].grid(True, which="both", alpha=0.3)

    fig.suptitle(title)
    return _img(fig, alt="block ranges")


def plot_snr_bars(results: list[StrategyResult], title: str) -> str:
    names = [r.name for r in results]
    snr = [r.snr_db for r in results]
    rmse = [r.rmse for r in results]
    fig, axes = plt.subplots(1, 2, figsize=(14, max(3.0, 0.35 * len(results) + 1)))

    idx = np.arange(len(names))
    axes[0].barh(idx, snr, color="C0")
    axes[0].set_yticks(idx)
    axes[0].set_yticklabels(names, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("SNR (dB) — higher is better")
    axes[0].grid(True, axis="x", alpha=0.3)
    axes[0].set_title("dequant SNR per strategy")

    axes[1].barh(idx, rmse, color="C3", log=True)
    axes[1].set_yticks(idx)
    axes[1].set_yticklabels(names, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("RMSE (log) — lower is better")
    axes[1].grid(True, axis="x", which="both", alpha=0.3)
    axes[1].set_title("round-trip RMSE per strategy")

    fig.suptitle(title)
    return _img(fig, alt="SNR/RMSE bars")


def plot_g1_vs_g2(g1: np.ndarray, g2: np.ndarray, title: str) -> str:
    """Hex-bin scatter of log10|g1| vs log10(g2) (sqrt-domain). For Adam,
    g2 ~ g1^2 in expectation when gradient is concentrated."""
    if g1.shape != g2.shape:
        return ""
    rng = np.random.default_rng(0)
    s_g1 = _subsample(g1.reshape(-1), 500_000, rng)
    s_g2 = _subsample(g2.reshape(-1), 500_000, rng)
    if s_g1.size != s_g2.size:
        m = min(s_g1.size, s_g2.size)
        s_g1, s_g2 = s_g1[:m], s_g2[:m]
    a1 = np.abs(s_g1)
    a2 = np.sqrt(np.maximum(s_g2, 0))
    keep = (a1 > 0) & (a2 > 0)
    a1, a2 = a1[keep], a2[keep]
    if a1.size == 0:
        return ""
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 5.0))
    hb = ax.hexbin(np.log10(a1), np.log10(a2), gridsize=80, cmap="viridis",
                   bins="log")
    ax.set_xlabel("log10 |g1|")
    ax.set_ylabel("log10 sqrt(g2)")
    fig.colorbar(hb, ax=ax, label="density (log)")
    # Reference: g2 == g1^2 -> sqrt(g2) == |g1|.
    lo, hi = ax.get_xlim()
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1, label="sqrt(g2) = |g1|")
    ax.legend(fontsize=8)
    ax.set_title(title)
    return _img(fig, alt="g1 vs g2 joint")


def plot_sh_structure(arr: np.ndarray, title: str) -> str:
    """For SH momentum (N, K, 3): per-coefficient and per-channel statistics."""
    if arr.ndim != 3 or arr.shape[-1] != 3:
        return ""
    N, K, _ = arr.shape
    abs_arr = np.abs(arr).astype(np.float64)
    per_coef = abs_arr.mean(axis=0)        # [K, 3]
    per_coef_max = abs_arr.max(axis=0)     # [K, 3]
    per_coef_std = abs_arr.std(axis=0)     # [K, 3]

    fig, axes = plt.subplots(2, 2, figsize=(13, 6.5))
    coefs = np.arange(K)
    for c, color, lab in zip(range(3), "rgb", "RGB"):
        axes[0, 0].bar(coefs + 0.25 * (c - 1), per_coef[:, c], width=0.25,
                       color=color, label=lab)
    axes[0, 0].set_title("mean |x| by SH coeff x channel")
    axes[0, 0].set_xlabel("SH coefficient index"); axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    for c, color, lab in zip(range(3), "rgb", "RGB"):
        axes[0, 1].bar(coefs + 0.25 * (c - 1), per_coef_max[:, c], width=0.25,
                       color=color, label=lab)
    axes[0, 1].set_title("max |x| by SH coeff x channel  (log)")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel("SH coefficient index"); axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Per-splat magnitude (max |x| across K and 3).
    per_splat = abs_arr.max(axis=(1, 2))
    axes[1, 0].hist(np.log10(per_splat[per_splat > 0]), bins=120, color="C4")
    axes[1, 0].set_title("log10 per-splat max |x|")
    axes[1, 0].set_yscale("log")
    axes[1, 0].grid(True, alpha=0.3)

    # CV (std/mean) per coeff index — high CV means a single block at that
    # index will span very different magnitudes across splats.
    cv = per_coef_std / np.maximum(per_coef, 1e-30)
    axes[1, 1].plot(coefs, cv.mean(axis=1), "k-o", label="mean over RGB")
    for c, color, lab in zip(range(3), "rgb", "RGB"):
        axes[1, 1].plot(coefs, cv[:, c], color=color, alpha=0.4, label=lab)
    axes[1, 1].set_title("per-coeff magnitude CV (std/mean over splats)")
    axes[1, 1].set_xlabel("SH coefficient index"); axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)

    fig.suptitle(title)
    return _img(fig, alt="SH structure")


def plot_bilagrid_structure(arr: np.ndarray, title: str) -> str:
    """For bilagrid momentum (N_cam, L, H, W, C): channel/spatial views.
    Channel-last layout: C is the innermost axis."""
    if arr.ndim != 5:
        return ""
    N, L, H, W, C = arr.shape
    abs_arr = np.abs(arr).astype(np.float64)

    # Per-channel mean and max across batch/spatial.
    per_C_mean = abs_arr.mean(axis=(0, 1, 2, 3))   # [C]
    per_C_max  = abs_arr.max(axis=(0, 1, 2, 3))    # [C]
    # Per-camera mean across spatial+C.
    per_cam_mean = abs_arr.mean(axis=(1, 2, 3, 4)) # [N]
    # Spatial: mean over batch+C, leave (L, H, W).
    spatial = abs_arr.mean(axis=(0, 4))            # [L, H, W]

    fig = plt.figure(figsize=(15, 6.5))
    gs = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.35)

    ax = fig.add_subplot(gs[0, 0])
    ax.bar(np.arange(C), per_C_mean, color="C0")
    ax.set_title(f"mean |x| per channel (C={C})")
    ax.set_yscale("log"); ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    ax.bar(np.arange(C), per_C_max, color="C3")
    ax.set_title("max |x| per channel"); ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    # Per-cam mean magnitude (sorted) -- shows whether camera mix dominates.
    ax = fig.add_subplot(gs[0, 2])
    cm = np.sort(per_cam_mean)
    ax.plot(np.arange(N), cm, color="C2")
    ax.set_yscale("log"); ax.grid(True, which="both", alpha=0.3)
    ax.set_title("per-camera mean |x| (sorted)")

    # Per-cam mean histogram (log).
    ax = fig.add_subplot(gs[0, 3])
    pos = per_cam_mean[per_cam_mean > 0]
    if pos.size:
        ax.hist(np.log10(pos), bins=60, color="C2")
    ax.set_title("log10 per-camera mean |x|")
    ax.set_yscale("log"); ax.grid(True, alpha=0.3)

    # Spatial mosaic: tile L slices in a row, normalized per L.
    ax = fig.add_subplot(gs[1, :])
    tile_h, tile_w = H, W
    pad = 2
    mosaic = np.full((tile_h, L * tile_w + (L - 1) * pad), np.nan)
    vmax = spatial.max() if spatial.max() > 0 else 1.0
    for l in range(L):
        mosaic[:, l * (tile_w + pad):l * (tile_w + pad) + tile_w] = spatial[l]
    im = ax.imshow(mosaic, cmap="viridis", aspect="auto",
                   norm=LogNorm(vmin=max(spatial[spatial > 0].min(), vmax * 1e-6),
                                vmax=vmax) if (spatial > 0).any() else None)
    ax.set_title(f"mean |x| over (batch, C), tiled across L={L}  "
                  f"(each tile is HxW={H}x{W}, log scale)")
    ax.set_yticks([]); ax.set_xticks([])
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    fig.suptitle(title)
    return _img(fig, alt="bilagrid structure")


# ---------------------------------------------------------------------------
# Section rendering
# ---------------------------------------------------------------------------

CSS = """
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif;
       margin: 24px auto; max-width: 1300px; padding: 0 16px; color: #222; }
h1 { border-bottom: 2px solid #333; padding-bottom: 6px; }
h2 { background: #eef; padding: 8px 10px; border-left: 4px solid #66a; margin-top: 32px; }
h3 { color: #246; margin-top: 24px; }
table { border-collapse: collapse; margin: 8px 0 16px 0; font-size: 12px; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
th { background: #f0f0f0; }
td.left, th.left { text-align: left; }
img { display: block; max-width: 100%; height: auto; margin: 6px 0 16px 0; }
.note { color: #555; background: #fafaee; border-left: 3px solid #cc9;
        padding: 6px 10px; margin: 8px 0; font-size: 13px; }
.bad  { color: #b22; font-weight: bold; }
.good { color: #1a6; font-weight: bold; }
.small { font-size: 12px; color: #666; }
pre { background: #f7f7f7; padding: 6px 10px; border: 1px solid #ddd; }
</style>
"""


def render_meta(meta: dict, ckpt_path: Path) -> str:
    rows = "".join(
        f'<tr><td class="left">{k}</td><td class="left">{v}</td></tr>'
        for k, v in meta.items()
    )
    return (
        f"<h2>Checkpoint</h2>"
        f"<p class='small'>Loaded from <code>{ckpt_path}</code>.</p>"
        f"<table>{rows}</table>"
    )


def parse_meta(path: Path) -> dict:
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def render_overview(fields: list[Field]) -> str:
    header = (
        "<tr><th class='left'>field</th><th class='left'>kind</th><th>shape</th>"
        "<th>numel</th><th>%zero</th><th>min</th><th>max</th>"
        "<th>mean |x|</th><th>p99 |x|</th><th>on-disk q?</th></tr>"
    )
    rows = [header]
    for f in fields:
        a = f.array
        finite = a[np.isfinite(a)]
        abs_a = np.abs(finite).astype(np.float64)
        p99 = float(np.percentile(abs_a, 99)) if abs_a.size else 0.0
        rows.append(
            f"<tr>"
            f"<td class='left'><code>{f.name}</code></td>"
            f"<td class='left'>{f.kind}</td>"
            f"<td>{a.shape}</td>"
            f"<td>{a.size:,}</td>"
            f"<td>{100 * float(np.mean(a == 0)):.2f}</td>"
            f"<td>{_fmt(float(finite.min()) if finite.size else 0)}</td>"
            f"<td>{_fmt(float(finite.max()) if finite.size else 0)}</td>"
            f"<td>{_fmt(float(abs_a.mean()) if abs_a.size else 0)}</td>"
            f"<td>{_fmt(p99)}</td>"
            f"<td>{'YES' if f.is_quantized_on_disk else '—'}</td>"
            f"</tr>"
        )
    return ("<h2>Buffer overview</h2>"
            "<table>" + "".join(rows) + "</table>"
            "<div class='note'>"
            "Fields marked <b>on-disk q? YES</b> were stored as uint8 + per-block "
            "float4 bounds; we re-inflated them via the bounds and labelled the "
            "section accordingly. Their absolute error simulations are upper "
            "bounds — they measure re-quantizing already-quantized values."
            "</div>")


def render_field(f: Field, results: list[StrategyResult],
                  plot_arr: Optional[np.ndarray] = None) -> str:
    """plot_arr defaults to f.array; pass the (possibly subsampled) array
    actually used in the strategy simulation for consistency with the table."""
    a = plot_arr if plot_arr is not None else f.array
    # Quick stats.
    finite = a[np.isfinite(a)]
    abs_a = np.abs(finite).astype(np.float64)
    p = [50, 90, 99, 99.9, 99.99]
    percs = np.percentile(abs_a, p) if abs_a.size else np.zeros(len(p))
    perc_html = "".join(
        f"<tr><td>{pp}%</td><td>{_fmt(float(v))}</td></tr>" for pp, v in zip(p, percs)
    )

    # Strategy table.
    strat_html = (
        "<tr><th class='left'>strategy</th><th>RMSE</th><th>MAE</th>"
        "<th>max err</th><th>SNR (dB)</th><th>block range mean</th>"
        "<th>p90 range</th><th>max range</th><th>range CV</th><th>note</th></tr>"
    )
    rows = [strat_html]
    best_snr = max(r.snr_db for r in results)
    for r in results:
        cls = "good" if abs(r.snr_db - best_snr) < 1e-6 else ""
        rows.append(
            f"<tr>"
            f"<td class='left'>{r.name}</td>"
            f"<td>{_fmt(r.rmse)}</td>"
            f"<td>{_fmt(r.mae)}</td>"
            f"<td>{_fmt(r.maxerr)}</td>"
            f"<td class='{cls}'>{r.snr_db:.2f}</td>"
            f"<td>{_fmt(r.block_range_mean)}</td>"
            f"<td>{_fmt(r.block_range_p90)}</td>"
            f"<td>{_fmt(r.block_range_max)}</td>"
            f"<td>{_fmt(r.block_range_cv)}</td>"
            f"<td class='left small'>{r.extra}</td>"
            f"</tr>"
        )

    head = (
        f"<h3>{f.name}  <span class='small'>(kind: {f.kind}, shape: {f.shape}, "
        f"dtype: {f.dtype}, "
        f"{'<span class=bad>on-disk-quantized</span>' if f.is_quantized_on_disk else 'float'}"
        f")</span></h3>"
    )

    parts = [head]

    # Percentile table + value distribution.
    parts.append(
        "<table style='display:inline-block; vertical-align:top; margin-right:20px'>"
        "<tr><th>|x| pct</th><th>value</th></tr>" + perc_html + "</table>"
    )
    parts.append(plot_value_histogram(a, f"value distribution — {f.name}"))

    # g1 vs g2 joint (filled in cross_field pass).

    # Structural plots.
    if f.kind in ("g1", "g2") and "features_sh" in f.name:
        parts.append(plot_sh_structure(a, f"SH structure — {f.name}"))
    if a.ndim == 5:
        parts.append(plot_bilagrid_structure(a, f"bilagrid structure — {f.name}"))

    # Strategy comparison.
    parts.append("<table>" + "".join(rows) + "</table>")
    parts.append(plot_snr_bars(results, f"strategy comparison — {f.name}"))
    parts.append(plot_block_ranges(results, f"block-range distributions — {f.name}"))

    if f.is_quantized_on_disk:
        parts.append(
            "<div class='note'>This buffer was stored quantized in the "
            "checkpoint. The RMSE / SNR numbers above measure additional error "
            "from <i>re-quantizing</i> the dequantized values, not the true "
            "vs-ground-truth error. Train this run with the quantize-optim "
            "flag <i>off</i> to get an apples-to-apples comparison.</div>"
        )

    return "".join(parts)


def render_cross_field(fields: list[Field]) -> str:
    """g1 vs g2 joint distribution where a matched pair exists."""
    by_base: dict[str, dict[str, Field]] = {}
    for f in fields:
        if f.kind in ("g1", "bilagrid_g1", "ppisp_g1"):
            base = f.name[3:] if f.name.startswith("g1_") else f.name.replace("_g1", "")
            by_base.setdefault(base, {})["g1"] = f
        elif f.kind in ("g2", "bilagrid_g2", "ppisp_g2"):
            base = f.name[3:] if f.name.startswith("g2_") else f.name.replace("_g2", "")
            by_base.setdefault(base, {})["g2"] = f

    parts = ["<h2>Adam moment correlation (g1 vs g2)</h2>",
             "<p class='small'>For a stable Adam state with bursty gradients, "
             "<code>sqrt(g2)</code> traces |g1|; broad scatter indicates "
             "gradient variance dominates the second-moment estimate, while "
             "tight collapse to the diagonal suggests rare large updates have "
             "settled. Both regimes affect quantization differently: tight "
             "scatter means joint-quantization (one bounds vector for both) "
             "would be efficient, broad scatter argues for separate ranges."
             "</p>"]
    for base, pair in by_base.items():
        if "g1" in pair and "g2" in pair:
            parts.append(f"<h3>{base}</h3>")
            parts.append(plot_g1_vs_g2(pair["g1"].array, pair["g2"].array, base))
    return "".join(parts)


def render_spatial_intro(spatial_perms: dict[str, np.ndarray],
                          world_means: np.ndarray) -> str:
    """Show the spatial sort: a few normalization views of world_means and
    the order in which Morton/Hilbert visit splats."""
    sample = world_means
    if sample.shape[0] > 200_000:
        idx = np.random.default_rng(0).choice(sample.shape[0], 200_000, replace=False)
        sample = sample[idx]
        sample_idx = idx
    else:
        sample_idx = np.arange(sample.shape[0])

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.3)

    norm_demo = [
        ("bbox", _norm_bbox(world_means)),
        ("quantile", _norm_quantile(world_means)),
        ("std_3sig", _norm_moment_std(world_means, 3.0)),
    ]
    for col, (name, nc) in enumerate(norm_demo):
        ax = fig.add_subplot(gs[0, col])
        ax.scatter(nc[sample_idx, 0], nc[sample_idx, 1], s=0.5, alpha=0.3, c="C0")
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")
        ax.set_title(f"normalization: {name}  (XY projection)")
        ax.grid(True, alpha=0.3)

    # Color splats by visit order under the most "principled" combos.
    interesting = [k for k in ("morton_quantile", "hilbert_quantile",
                                "hilbert_std_3sig") if k in spatial_perms]
    for col, sname in enumerate(interesting[:3]):
        ax = fig.add_subplot(gs[1, col])
        perm = spatial_perms[sname]
        # visit_order[i] = position at which splat i is visited.
        visit_order = np.empty_like(perm)
        visit_order[perm] = np.arange(perm.shape[0])
        nc = _norm_bbox(world_means)
        ax.scatter(nc[sample_idx, 0], nc[sample_idx, 1],
                   c=visit_order[sample_idx], s=0.5, alpha=0.5, cmap="viridis")
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")
        ax.set_title(f"visit order: {sname}")
        ax.grid(True, alpha=0.3)

    intro = ("<h2>Spatial sorting</h2>"
             "<p class='small'>For per-splat tensors (means, scales, quats, "
             "opacities, dc, sh, plus their g1/g2 momentum), we precompute a "
             "single permutation along axis 0 derived from "
             "<code>world_means</code>, then apply the regular identity quantizer. "
             "Two space-filling curves (Morton/z-order and Hilbert) are paired "
             "with five normalizations of the coordinates: <code>bbox</code> "
             "(per-axis min/max), <code>quantile</code> (per-axis CDF), "
             "<code>std_Nsig</code> (mean ± N·sigma, clip), <code>mad_3sig</code> "
             "(robust to outliers via median + MAD). The same permutation is "
             "reused across all per-splat fields, so the decode cost amortizes — "
             "one extra shared index table for the whole optimizer state.</p>")
    return intro + _img(fig, alt="spatial-sort visualization")


def render_spatial_ranking(per_field: list[tuple[Field, list[StrategyResult]]]) -> str:
    """Heatmap-style table: rows=per-splat fields, columns=spatial strategies,
    cells coloured by Δ SNR vs identity. Identifies the winning (curve, norm)
    combination per field at a glance."""
    rows: list[tuple[str, dict[str, float], float, str, float]] = []
    for f, res in per_field:
        if not any(r.name.startswith("spatial_") for r in res):
            continue
        baseline = next((r.snr_db for r in res if r.name == "identity"),
                         res[0].snr_db)
        d = {r.name: r.snr_db - baseline for r in res
             if r.name.startswith("spatial_")}
        best = max(res, key=lambda r: r.snr_db)
        rows.append((f.name, d, baseline, best.name, best.snr_db - baseline))

    if not rows:
        return ""

    # Cell colouring: deepest green = biggest improvement.
    all_strats = sorted({k for _, d, *_ in rows for k in d})
    max_gain = max((max(d.values()) if d else 0.0) for _, d, *_ in rows)
    max_gain = max(max_gain, 1.0)

    def cell_style(v: float) -> str:
        if v <= 0:
            t = max(min(-v / max_gain, 1.0), 0.0)
            r, g, b = 255, int(255 * (1 - 0.5 * t)), int(255 * (1 - 0.5 * t))
        else:
            t = min(v / max_gain, 1.0)
            r, g, b = int(255 * (1 - 0.6 * t)), 255, int(255 * (1 - 0.6 * t))
        return f"background:rgb({r},{g},{b})"

    header = ("<tr><th class='left'>field</th><th>identity SNR</th>"
              + "".join(f"<th>{s.replace('spatial_', '')}</th>" for s in all_strats)
              + "<th>best overall</th><th>Δ SNR</th></tr>")
    body = []
    for name, d, baseline, best_name, best_gain in rows:
        cells = "".join(
            f"<td style='{cell_style(d.get(s, 0.0))}'>{d.get(s, 0.0):+.2f}</td>"
            for s in all_strats
        )
        body.append(
            f"<tr><td class='left'><code>{name}</code></td>"
            f"<td>{baseline:.2f}</td>{cells}"
            f"<td class='left small'>{best_name}</td>"
            f"<td><b>{best_gain:+.2f}</b></td></tr>"
        )
    return ("<h2>Spatial-sort ranking</h2>"
            "<p class='small'>Cell values are Δ SNR (dB) vs the identity "
            "baseline; deeper green = larger improvement, red = regression. "
            "<code>best overall</code> reports the winning strategy across "
            "all categories (axis perms, sort-within, oracle, spatial).</p>"
            "<table>" + header + "".join(body) + "</table>")


# ---------------------------------------------------------------------------
# Joint (g1, g2) pair analysis. Per-component SNR can mislead about training
# dynamics: Adam's actual update direction is u = g1 / (sqrt(g2) + eps), and a
# given quantization scheme can have great per-component SNR but waste bits on
# components of (g1, g2) that cancel in the ratio. Here we compare schemes by
# how well they preserve u.
#
# A scheme takes (g1, g2) and outputs two uint8 streams + per-block float4
# bounds (= 2 bytes/cell + 16 bytes/256-cell block, matching the on-device
# storage). The two streams need not be independent: in some schemes we
# encode (u, sqrt(g2)) directly, then reconstruct g1 from the pair at decode.
# ---------------------------------------------------------------------------

ADAM_EPS = 1e-15          # Adam denominator epsilon
ADAM_BETA1 = 0.9
ADAM_BETA2 = 0.999


def _quantize_field_simple(arr: np.ndarray, mapping: str) -> np.ndarray:
    """One-shot identity-ordering 256-block uint8 round-trip on a single
    tensor under `mapping`; returns dequantized values reshaped to `arr`."""
    flat = arr.reshape(-1).astype(np.float32)
    n = flat.size
    blocks = _pad_to_blocks(flat)
    deq, _ = _block_quantize(blocks, mapping)
    return deq.reshape(-1)[:n].reshape(arr.shape)


def _adam_update(g1: np.ndarray, g2: np.ndarray, eps: float = ADAM_EPS,
                  bc_step: Optional[int] = None) -> np.ndarray:
    """Reproduce Adam's per-cell update direction. When bc_step is given,
    apply the standard bias correction (g1_hat = g1 / (1 - beta1^t), etc.)."""
    g1_eff = g1
    g2_eff = g2
    if bc_step is not None and bc_step > 0:
        bc1 = 1.0 - ADAM_BETA1 ** bc_step
        bc2 = 1.0 - ADAM_BETA2 ** bc_step
        g1_eff = g1 / bc1
        g2_eff = g2 / bc2
    return g1_eff / (np.sqrt(np.maximum(g2_eff, 0.0)) + eps)


def _adam_one_step_u(g1: np.ndarray, g2: np.ndarray, v: np.ndarray,
                      eps: float = ADAM_EPS,
                      bc_step: Optional[int] = None) -> np.ndarray:
    """Apply ONE Adam EMA step (g1 → β1·g1 + (1-β1)·v, g2 → β2·g2 + (1-β2)·v²)
    and return the resulting update direction `u_next = g1_new / (sqrt(g2_new) + eps)`,
    optionally bias-corrected at step `bc_step`.

    Unlike the round-trip update direction `g1/(sqrt(g2)+eps)`, this metric
    depends on BOTH the g1 and g2 quantization errors because the gradient
    `v` mixes them through the EMA — it breaks the analytic identity that
    causes single-step SNR(u) to be invariant to the sqrt_g2 mapping under
    the joint (u, sqrt_g2) scheme."""
    g1_new = ADAM_BETA1 * g1 + (1.0 - ADAM_BETA1) * v
    g2_new = ADAM_BETA2 * g2 + (1.0 - ADAM_BETA2) * v * v
    if bc_step is not None and bc_step > 0:
        bc1 = 1.0 - ADAM_BETA1 ** bc_step
        bc2 = 1.0 - ADAM_BETA2 ** bc_step
        g1_new = g1_new / bc1
        g2_new = g2_new / bc2
    return g1_new / (np.sqrt(np.maximum(g2_new, 0.0)) + eps)


@dataclass
class PairSchemeResult:
    name: str
    extra: str
    snr_u_unbiased: float        # round-trip SNR(u) without bias correction
    snr_u_biased: float          # round-trip SNR(u) at a representative step
    snr_u_next_unbiased: float   # post-one-Adam-step SNR(u) unbiased
    snr_u_next_biased: float     # post-one-Adam-step SNR(u) bias-corrected
    snr_g1: float                # per-component (cross-reference)
    snr_g2: float
    rmse_u: float                # raw RMSE on u
    bytes_per_cell: int = 2      # 2 uint8 + amortized float4 bounds


def _snr_pair(orig: np.ndarray, recon: np.ndarray) -> tuple[float, float]:
    """Return (rmse, snr_db) computed in float64 on flattened inputs."""
    o = orig.reshape(-1).astype(np.float64)
    r = recon.reshape(-1).astype(np.float64)
    finite = np.isfinite(o) & np.isfinite(r)
    if not finite.any():
        return 0.0, 0.0
    o = o[finite]; r = r[finite]
    diff = o - r
    rmse = float(np.sqrt(np.mean(diff * diff)))
    orig_rms = float(np.sqrt(np.mean(o * o)))
    snr = 20.0 * np.log10(max(orig_rms, 1e-30) / max(rmse, 1e-30))
    return rmse, snr


def analyze_pair_schemes(g1: np.ndarray, g2: np.ndarray,
                          eps: float = ADAM_EPS,
                          bc_step: int = 10_000,
                          max_workers: int = 8,
                          ) -> list[PairSchemeResult]:
    """Try several encodings of an Adam (g1, g2) pair, score reconstruction
    of the Adam update direction u. Each scheme uses 16 bits per cell
    (= 2 uint8 streams) plus shared per-block bounds.

    Implementation note: a scheme's two halves are quantizations of one of the
    primitive tensors {g1, sqrt_g2, u} under one of {linear, signed_log,
    sqrt, log_pos, log_sqrt}. We compute the (primitive, mapping) → quantized
    tensor table once, in parallel, then assemble each scheme's (g1_q, g2_q)
    from cached entries. This roughly halves the work vs the obvious
    'one quantize per scheme half' loop."""

    sqrt_g2_true = np.sqrt(np.maximum(g2, 0.0))
    u_unb = g1 / (sqrt_g2_true + eps)
    u_true_unbiased = u_unb
    u_true_biased = _adam_update(g1, g2, eps, bc_step=bc_step)

    # Synthetic gradient for the post-Adam-step metric: per-cell Gaussian
    # with variance ≈ g2_true (matches steady-state Adam's E[v²] = g2 expectation).
    # Same v is reused across all schemes so differences reflect only the
    # quantization choice, not the gradient sample.
    rng = np.random.default_rng(0)
    v_synth = (rng.standard_normal(g1.shape).astype(np.float32)
               * sqrt_g2_true.astype(np.float32))
    u_next_true_unbiased = _adam_one_step_u(g1, g2, v_synth, eps, bc_step=None)
    u_next_true_biased   = _adam_one_step_u(g1, g2, v_synth, eps, bc_step=bc_step + 1)

    # Declare every (primitive, mapping) we'll need across all schemes.
    needed: dict[str, tuple[np.ndarray, str]] = {
        "g1_linear":         (g1, "linear"),
        "g1_signed_log":     (g1, "signed_log"),
        "g2_sqrt":           (g2, "sqrt"),
        "g2_log_sqrt":       (g2, "log_sqrt"),
        "u_linear":          (u_unb, "linear"),
        "u_signed_log":      (u_unb, "signed_log"),
        "sqrt_g2_linear":    (sqrt_g2_true, "linear"),
        "sqrt_g2_log_pos":   (sqrt_g2_true, "log_pos"),
    }
    cache: dict[str, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(_quantize_field_simple, arr, mapping): key
                 for key, (arr, mapping) in needed.items()}
        for f in as_completed(futs):
            cache[futs[f]] = f.result()

    # Scheme list as (name, extra, g1_q_factory, g2_q_factory). Factories take
    # the cache and return the reconstructed (g1_q, g2_q) tensors.
    indep_specs = [
        ("indep_linear_sqrt",         "current on-device scheme",
            "g1_linear",      "g2_sqrt"),
        ("indep_signed_log_sqrt",     "g1 via signed log; g2 via sqrt (current)",
            "g1_signed_log",  "g2_sqrt"),
        ("indep_signed_log_log_sqrt", "g1 via signed log; sqrt(g2) via log",
            "g1_signed_log",  "g2_log_sqrt"),
        ("indep_linear_log_sqrt",     "g1 linear; sqrt(g2) via log",
            "g1_linear",      "g2_log_sqrt"),
    ]
    joint_specs = [
        ("joint_u_sqrtg2_linear",     "store (u, sqrt(g2)) linearly",
            "u_linear",       "sqrt_g2_linear"),
        ("joint_u_sqrtg2_signed_log", "store (signed_log(u), sqrt(g2))",
            "u_signed_log",   "sqrt_g2_linear"),
        ("joint_u_sqrtg2_log_log",    "store (signed_log(u), log(sqrt(g2)))",
            "u_signed_log",   "sqrt_g2_log_pos"),
        ("joint_u_sqrtg2_linear_log", "store (u, log(sqrt(g2)))",
            "u_linear",       "sqrt_g2_log_pos"),
    ]

    def assemble_indep(g1_key: str, g2_key: str) -> tuple[np.ndarray, np.ndarray]:
        return cache[g1_key], cache[g2_key]

    def assemble_joint(u_key: str, s_key: str) -> tuple[np.ndarray, np.ndarray]:
        u_q = cache[u_key]
        sqrt_g2_q = cache[s_key]
        g2_q = sqrt_g2_q * sqrt_g2_q
        g1_q = u_q * (sqrt_g2_q + eps)
        return g1_q, g2_q

    pending: list[tuple[str, str, np.ndarray, np.ndarray]] = []
    for name, extra, k1, k2 in indep_specs:
        g1_q, g2_q = assemble_indep(k1, k2)
        pending.append((name, extra, g1_q, g2_q))
    for name, extra, k1, k2 in joint_specs:
        g1_q, g2_q = assemble_joint(k1, k2)
        pending.append((name, extra, g1_q, g2_q))

    def score_one(name: str, extra: str,
                   g1_q: np.ndarray, g2_q: np.ndarray) -> PairSchemeResult:
        u_q_unb = _adam_update(g1_q, g2_q, eps, bc_step=None)
        u_q_bia = _adam_update(g1_q, g2_q, eps, bc_step=bc_step)
        u_next_q_unb = _adam_one_step_u(g1_q, g2_q, v_synth, eps, bc_step=None)
        u_next_q_bia = _adam_one_step_u(g1_q, g2_q, v_synth, eps, bc_step=bc_step + 1)
        rmse_u, snr_u = _snr_pair(u_true_unbiased, u_q_unb)
        _, snr_u_bc = _snr_pair(u_true_biased, u_q_bia)
        _, snr_u_next_unb = _snr_pair(u_next_true_unbiased, u_next_q_unb)
        _, snr_u_next_bia = _snr_pair(u_next_true_biased, u_next_q_bia)
        _, snr_g1 = _snr_pair(g1, g1_q)
        _, snr_g2 = _snr_pair(g2, g2_q)
        return PairSchemeResult(
            name=name, extra=extra,
            snr_u_unbiased=snr_u, snr_u_biased=snr_u_bc,
            snr_u_next_unbiased=snr_u_next_unb, snr_u_next_biased=snr_u_next_bia,
            snr_g1=snr_g1, snr_g2=snr_g2,
            rmse_u=rmse_u,
        )

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(score_one, n, e, g1q, g2q)
                 for n, e, g1q, g2q in pending]
        results = [f.result() for f in futs]

    results.sort(key=lambda r: -r.snr_u_unbiased)
    return results


def _plot_pair_joint(g1: np.ndarray, g2: np.ndarray,
                      title: str, eps: float = ADAM_EPS) -> str:
    """Hex-bin (g1, sqrt(g2)) joint + (u, sqrt(g2)) joint side by side.
    Shows whether the Adam update direction has structure that joint
    encoding can exploit."""
    rng = np.random.default_rng(0)
    n = min(500_000, g1.size)
    idx = rng.choice(g1.size, size=n, replace=False)
    s_g1 = g1.reshape(-1)[idx].astype(np.float64)
    s_g2 = g2.reshape(-1)[idx].astype(np.float64)
    s_sqrt_g2 = np.sqrt(np.maximum(s_g2, 0.0))
    s_u = s_g1 / (s_sqrt_g2 + eps)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # Left: log|g1| vs log(sqrt(g2)). For decoupled gaussian noise these
    # would be ~independent; correlation is the structure to exploit.
    a1 = np.abs(s_g1); s2 = s_sqrt_g2
    keep = (a1 > 0) & (s2 > 0)
    a1, s2_l = a1[keep], s2[keep]
    if a1.size:
        hb = axes[0].hexbin(np.log10(a1), np.log10(s2_l),
                             gridsize=80, cmap="viridis", bins="log")
        fig.colorbar(hb, ax=axes[0], label="density (log)")
        lo, hi = axes[0].get_xlim()
        axes[0].plot([lo, hi], [lo, hi], "r--", linewidth=1,
                     label="sqrt(g2) = |g1|")
        axes[0].set_xlabel("log10 |g1|")
        axes[0].set_ylabel("log10 sqrt(g2)")
        axes[0].set_title("|g1| vs sqrt(g2) joint")
        axes[0].legend(fontsize=8)

    # Right: log(sqrt(g2)) vs Adam update direction u = g1/(sqrt(g2)+eps).
    # If u has small dynamic range, joint encoding wins.
    s2 = s_sqrt_g2
    a_u = np.abs(s_u)
    keep = (s2 > 0) & (a_u > 0)
    if keep.any():
        hb = axes[1].hexbin(np.log10(s2[keep]), np.log10(a_u[keep]),
                             gridsize=80, cmap="viridis", bins="log")
        fig.colorbar(hb, ax=axes[1], label="density (log)")
        axes[1].set_xlabel("log10 sqrt(g2)")
        axes[1].set_ylabel("log10 |u| = |g1| / (sqrt(g2)+eps)")
        axes[1].set_title("u vs sqrt(g2)  (what joint schemes exploit)")
    fig.suptitle(title)
    return _img(fig, alt="pair joint")


def _subsample_pair(g1: np.ndarray, g2: np.ndarray,
                     max_numel: int = 20_000_000) -> tuple[np.ndarray, np.ndarray]:
    """Cap pair analysis size by subsampling along the outer axis. The same
    rows are kept for both arrays so per-cell correspondence is preserved."""
    if g1.size <= max_numel or g1.ndim < 1 or g1.shape[0] <= 1:
        return g1, g2
    keep = max(int(max_numel / max(np.prod(g1.shape[1:]), 1)), 1)
    rng = np.random.default_rng(0)
    sel = rng.choice(g1.shape[0], size=min(keep, g1.shape[0]), replace=False)
    sel.sort()
    return g1[sel].copy(), g2[sel].copy()


def render_pair_section(fields: list[Field],
                         bc_step: int = 10_000,
                         compute_workers: int = 8,
                         plot_workers: int = 16,
                         max_numel: int = 20_000_000) -> str:
    """Pair every g1_* / g2_* across the checkpoint and report the joint
    schemes' update-direction SNR. Each pair is analyzed in its own thread;
    within a pair, schemes run concurrently. Plot rendering runs on a
    separate thread budget."""
    by_base: dict[str, dict[str, np.ndarray]] = {}
    for f in fields:
        if f.kind in ("g1", "bilagrid_g1", "ppisp_g1"):
            base = (f.name[3:] if f.name.startswith("g1_")
                     else f.name.replace("_g1", ""))
            by_base.setdefault(base, {})["g1"] = f.array
        elif f.kind in ("g2", "bilagrid_g2", "ppisp_g2"):
            base = (f.name[3:] if f.name.startswith("g2_")
                     else f.name.replace("_g2", ""))
            by_base.setdefault(base, {})["g2"] = f.array

    pairs = []
    for k, p in by_base.items():
        if "g1" not in p or "g2" not in p or p["g1"].shape != p["g2"].shape:
            continue
        g1_rms = float(np.sqrt(np.mean(p["g1"].astype(np.float64) ** 2)))
        g2_rms = float(np.sqrt(np.mean(p["g2"].astype(np.float64) ** 2)))
        if g1_rms < 1e-20 and g2_rms < 1e-20:
            continue
        g1_s, g2_s = _subsample_pair(p["g1"], p["g2"], max_numel=max_numel)
        pairs.append((k, g1_s, g2_s, p["g1"].shape != g1_s.shape, p["g1"].shape))

    if not pairs:
        return ""

    print(f"  pair analysis: {len(pairs)} pairs "
           f"(compute_j={compute_workers}, plot_j={plot_workers})")

    # Schemes per pair run inside analyze_pair_schemes via a thread pool
    # capped at compute_workers. To run pairs themselves in parallel without
    # multiplying that cost N times, we use a small top-level pool whose
    # workers each sequence their own analyze (which then internally fans
    # out to compute_workers). max_pair_concurrency is conservative: it lets
    # 2 pairs progress at once, giving overlap between the slow SH pairs
    # and the fast scalar ones without blowing the RAM budget.
    max_pair_concurrency = max(1, min(len(pairs), 2))

    def analyze_one(base, g1, g2):
        t0 = time.time()
        r = analyze_pair_schemes(g1, g2, bc_step=bc_step,
                                  max_workers=compute_workers)
        return base, g1, g2, r, time.time() - t0

    pair_results: list[tuple] = [None] * len(pairs)
    with ThreadPoolExecutor(max_workers=max_pair_concurrency) as ex:
        futs = {ex.submit(analyze_one, b, g1, g2): i
                 for i, (b, g1, g2, *_) in enumerate(pairs)}
        for f in as_completed(futs):
            i = futs[f]
            base, g1, g2, results, dt = f.result()
            pair_results[i] = (base, g1, g2, results, dt)
            print(f"    {base:24s} {len(results)} schemes in {dt:.1f}s")

    # --- Render plots in parallel ---
    plot_jobs = []   # (pair_idx, kind, args)
    for i, (base, g1, g2, results, _dt) in enumerate(pair_results):
        plot_jobs.append((i, "bars", (results, base)))
        if base in ("means", "scales", "features_dc", "features_sh"):
            plot_jobs.append((i, "joint", (g1, g2,
                              f"(g1, g2) joint distribution — {base}")))
    plots_html: dict[tuple[int, str], str] = {}
    def render_plot(kind, args):
        if kind == "bars":
            return _plot_pair_snr_bars(*args)
        elif kind == "joint":
            return _plot_pair_joint(*args)
        return ""
    with ThreadPoolExecutor(max_workers=plot_workers) as ex:
        futs = {ex.submit(render_plot, k, a): (i, k) for i, k, a in plot_jobs}
        for f in as_completed(futs):
            plots_html[futs[f]] = f.result()

    parts = ["<h2>Joint (g1, g2) pair analysis: Adam update direction</h2>",
             "<p class='small'>Each cell of the on-device optimizer state "
             "holds 2 bytes (one for g1, one for g2). Independent quantization "
             "(the current scheme) can have great per-component SNR but waste "
             "bits on (g1, g2) variations that cancel in Adam's update "
             "<code>u = g1 / (sqrt(g2) + eps)</code>. We compare independent "
             "vs joint encodings on how well they preserve u — the quantity "
             "training actually consumes. "
             f"<b>Bias correction</b> applied at step "
             f"<code>t={bc_step}</code>: "
             f"<code>(g1/(1-β1^t)) / (sqrt(g2/(1-β2^t)) + eps)</code> with "
             f"<code>eps={ADAM_EPS}</code>, "
             f"<code>β1=0.9, β2=0.999</code>. SNR is scale-invariant, so the "
             f"unbiased and bias-corrected columns differ only when "
             f"<code>eps</code> dominates one denominator but not the other "
             f"(typically only for extremely small <code>sqrt(g2)</code>).</p>"]

    for i, (base, g1, g2, results, _dt) in enumerate(pair_results):
        sub_note = ""
        # pairs[i] tuple was (base, g1, g2, was_subsampled, orig_shape)
        was_subsampled = pairs[i][3]
        orig_shape = pairs[i][4]
        if was_subsampled:
            sub_note = (f" <span class='small'>[subsampled to {g1.shape} "
                         f"from {orig_shape}]</span>")
        parts.append(f"<h3>{base}  <span class='small'>"
                      f"(shape {g1.shape})</span>{sub_note}</h3>")
        head = ("<tr><th class='left'>scheme</th>"
                "<th>SNR(u)  unbiased</th>"
                f"<th>SNR(u)  bc step {bc_step}</th>"
                "<th>SNR(u_next)  unbiased</th>"
                f"<th>SNR(u_next)  bc step {bc_step + 1}</th>"
                "<th>SNR(g1)</th><th>SNR(g2)</th>"
                "<th>RMSE(u)</th><th class='left'>note</th></tr>")
        rows = [head]
        # Sort the bar order by post-step SNR — the round-trip SNR(u) is
        # degenerate across joint schemes that share the u mapping (see note
        # below); SNR(u_next) propagates both halves through one Adam step.
        best_next = max(r.snr_u_next_unbiased for r in results)
        baseline_idx = next((j for j, r in enumerate(results)
                              if r.name == "indep_linear_sqrt"), None)
        baseline = (results[baseline_idx].snr_u_next_unbiased
                    if baseline_idx is not None else 0.0)
        for r in results:
            cls = "good" if abs(r.snr_u_next_unbiased - best_next) < 1e-6 else ""
            delta = r.snr_u_next_unbiased - baseline
            rows.append(
                f"<tr><td class='left'><code>{r.name}</code></td>"
                f"<td>{r.snr_u_unbiased:.2f}</td>"
                f"<td>{r.snr_u_biased:.2f}</td>"
                f"<td class='{cls}'>{r.snr_u_next_unbiased:.2f} "
                 f"<span class='small'>({delta:+.2f})</span></td>"
                f"<td>{r.snr_u_next_biased:.2f}</td>"
                f"<td>{r.snr_g1:.2f}</td>"
                f"<td>{r.snr_g2:.2f}</td>"
                f"<td>{_fmt(r.rmse_u)}</td>"
                f"<td class='left small'>{r.extra}</td></tr>"
            )
        parts.append("<table>" + "".join(rows) + "</table>")
        parts.append(
            "<div class='note'>"
            "<b>Why does SNR(u) collapse across joint variants?</b> When the "
            "joint scheme stores <code>u</code> directly and decodes "
            "<code>g1 = u_q · (sqrt_g2_q + ε)</code> / <code>g2 = sqrt_g2_q²</code>, "
            "the round-trip <code>u_recon = g1 / (sqrt(g2) + ε) = u_q · "
            "(sqrt_g2_q + ε) / (sqrt_g2_q + ε) = u_q</code> — exactly. So "
            "<b>SNR(u) is a function of the u mapping only</b> and is "
            "identical across (e.g.) <code>joint_u_sqrtg2_linear</code> and "
            "<code>joint_u_sqrtg2_linear_log</code>. To see the sqrt_g2 "
            "mapping's actual effect on training dynamics, look at "
            "<b>SNR(u_next)</b>: it computes one Adam EMA step "
            "<code>g1' = β1·g1 + (1-β1)·v</code>, "
            "<code>g2' = β2·g2 + (1-β2)·v²</code> with a synthetic gradient "
            "<code>v ~ N(0, g2_true)</code> shared across schemes, then "
            "evaluates SNR on <code>u' = g1' / (sqrt(g2') + ε)</code>. Because "
            "the new gradient term mixes <code>g1</code> and <code>g2</code> "
            "through the EMA, errors in either propagate to <code>u'</code> "
            "and the metric becomes sensitive to both quantization choices. "
            "Per-component <code>SNR(g1)</code> and <code>SNR(g2)</code> also "
            "distinguish the schemes."
            "</div>"
        )
        parts.append(plots_html.get((i, "bars"), ""))
        joint_html = plots_html.get((i, "joint"))
        if joint_html:
            parts.append(joint_html)
    return "".join(parts)


def _plot_pair_snr_bars(results: list[PairSchemeResult], base: str) -> str:
    """Three side-by-side bar groups:
      (left)   round-trip SNR(u)       — flat across joint variants that share u-mapping
      (middle) post-one-Adam-step SNR  — separates the variants
      (right)  per-component SNR(g1) / SNR(g2)."""
    names = [r.name for r in results]
    fig, axes = plt.subplots(1, 3, figsize=(17, max(2.8, 0.34 * len(results) + 1)))

    idx = np.arange(len(names))
    w = 0.35

    axes[0].barh(idx - w / 2, [r.snr_u_unbiased for r in results], w,
                  color="C0", label="unbiased")
    axes[0].barh(idx + w / 2, [r.snr_u_biased for r in results], w,
                  color="C2", label="bias-corrected")
    axes[0].set_yticks(idx); axes[0].set_yticklabels(names, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("SNR(u) round-trip (dB)")
    axes[0].grid(True, axis="x", alpha=0.3)
    axes[0].set_title(f"SNR(u) — {base}\n(degenerate across joint variants)")
    axes[0].legend(fontsize=8)

    axes[1].barh(idx - w / 2, [r.snr_u_next_unbiased for r in results], w,
                  color="C4", label="unbiased")
    axes[1].barh(idx + w / 2, [r.snr_u_next_biased for r in results], w,
                  color="C5", label="bias-corrected")
    axes[1].set_yticks(idx); axes[1].set_yticklabels(names, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("SNR(u_next) after one Adam step (dB)")
    axes[1].grid(True, axis="x", alpha=0.3)
    axes[1].set_title("SNR(u_next)\n(propagates BOTH g1 and g2 error)")
    axes[1].legend(fontsize=8)

    axes[2].barh(idx - w / 2, [r.snr_g1 for r in results], w,
                  color="C1", label="g1")
    axes[2].barh(idx + w / 2, [r.snr_g2 for r in results], w,
                  color="C3", label="g2")
    axes[2].set_yticks(idx); axes[2].set_yticklabels(names, fontsize=8)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("per-component SNR (dB)")
    axes[2].grid(True, axis="x", alpha=0.3)
    axes[2].set_title("per-component SNR  (cross-reference)")
    axes[2].legend(fontsize=8)

    return _img(fig, alt="pair snr bars")


def render_summary(per_field: list[tuple[Field, list[StrategyResult]]]) -> str:
    """Aggregated summary: SNR improvement vs identity for each field, side-by-side."""
    fig, ax = plt.subplots(1, 1, figsize=(12, max(3.0, 0.42 * len(per_field) + 1)))
    field_names = []
    improvements = []          # list of dicts {strat -> delta-snr}
    strats_seen: list[str] = []
    for f, res in per_field:
        if not res:
            continue
        baseline = next((r.snr_db for r in res if r.name == "identity"), res[0].snr_db)
        d = {}
        for r in res:
            d[r.name] = r.snr_db - baseline
            if r.name not in strats_seen:
                strats_seen.append(r.name)
        field_names.append(f.name)
        improvements.append(d)

    # Plot grouped bars (one bar per strategy per field).
    y = np.arange(len(field_names))
    h = 0.8 / max(len(strats_seen), 1)
    for i, s in enumerate(strats_seen):
        vals = [imp.get(s, 0.0) for imp in improvements]
        ax.barh(y + i * h - 0.4 + h / 2, vals, h, label=s)
    ax.set_yticks(y)
    ax.set_yticklabels(field_names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Δ SNR vs identity (dB) — positive = better than current scheme")
    ax.axvline(0, color="k", linewidth=0.6)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    ax.set_title("strategy improvement over identity (memory order), per field")
    return ("<h2>Summary: strategy improvement across fields</h2>"
            + _img(fig, alt="summary"))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ckpt", help="path to a step-XXXXXXXXX.ckpt directory "
                                       "(must contain a 'full/' subfolder)")
    parser.add_argument("--output", "-o", default=None,
                        help="output HTML path (default: <ckpt>/quantization_report.html)")
    parser.add_argument("--skip", nargs="*", default=[],
                        help="field names to skip (e.g. world_features_sh)")
    parser.add_argument("--only", nargs="*", default=None,
                        help="restrict to listed field names")
    parser.add_argument("--max-numel", type=int, default=20_000_000,
                        help="cap per-field numel for strategy simulation "
                             "(subsamples deterministically when exceeded)")
    parser.add_argument("--threads", "-j", type=int,
                        default=min(8, os.cpu_count() or 4),
                        help="strategies to evaluate in parallel per field "
                             "(numpy releases the GIL for the hot inner loops)")
    parser.add_argument("--no-spatial", action="store_true",
                        help="skip Morton/Hilbert spatial-sort strategies "
                             "(faster, smaller report)")
    parser.add_argument("--no-pair", action="store_true",
                        help="skip joint (g1, g2) pair analysis section")
    parser.add_argument("--bc-step", type=int, default=10_000,
                        help="step number at which to apply bias correction "
                             "in the pair-analysis update-SNR metric")
    parser.add_argument("--plot-threads", type=int,
                        default=min(16, (os.cpu_count() or 4) * 2),
                        help="thread budget for matplotlib plot rendering "
                             "(cheap on memory, can run wider than --threads)")
    parser.add_argument("--pair-max-numel", type=int, default=20_000_000,
                        help="cap per-pair numel for joint pair analysis "
                             "(outer-axis subsample when exceeded)")
    args = parser.parse_args(argv)

    ckpt_path = Path(args.ckpt).resolve()
    full = ckpt_path / "full"
    if not full.is_dir():
        sys.exit(f"error: {full} does not exist — pass a checkpoint dir written "
                  f"with --save-full-checkpoint")

    meta = parse_meta(ckpt_path / "meta.txt")
    fields = discover_fields(full, meta=meta)

    if args.only:
        fields = [f for f in fields if f.name in args.only]
    if args.skip:
        fields = [f for f in fields if f.name not in args.skip]

    # Skip params for strategy sim (just include in overview); momentum is the
    # interesting target. But we still show their stats for context.
    momentum_kinds = {"g1", "g2", "bilagrid_g1", "bilagrid_g2", "ppisp_g1", "ppisp_g2"}

    # --- Build spatial permutations once from world_means ----------------
    world_means_field = next((f for f in fields if f.name == "world_means"), None)
    n_splats = None
    spatial_perms_global: dict[str, np.ndarray] = {}
    if world_means_field is not None and world_means_field.array.ndim == 2:
        n_splats = int(meta.get("cur_num_splats", world_means_field.array.shape[0]))
        if not args.no_spatial and n_splats > 1:
            t0 = time.time()
            print(f"Building spatial sort permutations on {n_splats:,} splats "
                   f"(nbits={NBITS_SPATIAL}) ...")
            spatial_perms_global = build_spatial_permutations(
                world_means_field.array, n_splats)
            print(f"  built {len(spatial_perms_global)} perms in "
                   f"{time.time() - t0:.1f}s "
                   f"({', '.join(spatial_perms_global)})")

    print(f"Loaded {len(fields)} fields from {full}")
    out_html = [
        "<!DOCTYPE html><html><head>",
        "<meta charset='utf-8'><title>Quantization investigation</title>",
        CSS,
        "</head><body>",
        f"<h1>Optimizer-state quantization investigation</h1>",
        "<p class='small'>Simulating the current on-device scheme "
        "(<code>BLOCK_SIZE=256, 8-bit per cell, float4 per-block bounds, "
        "sqrt(g2) before quantize</code>) under alternative block orderings, "
        "to test whether grouping similar-magnitude parameters into the "
        "same block reduces quantization error. The Morton / Hilbert strategies "
        "reorder splats along a 3D space-filling curve over <code>world_means</code> "
        "before applying the regular identity quantizer — the hypothesis being that "
        "spatially-close splats optimize similar features, so their momentum "
        "magnitudes are correlated.</p>",
        render_meta(meta, ckpt_path),
        render_overview(fields),
    ]

    if spatial_perms_global:
        out_html.append(render_spatial_intro(spatial_perms_global,
                                              world_means_field.array[:n_splats]))

    per_field_results: list[tuple[Field, list[StrategyResult]]] = []
    field_subsampled_arrays: list[np.ndarray] = []
    out_html.append("<h2>Per-buffer analysis</h2>")
    for f in fields:
        if f.kind not in momentum_kinds:
            continue
        a = f.array
        is_per_splat = (n_splats is not None and a.ndim >= 1
                        and a.shape[0] == int(meta.get("max_num_splats", a.shape[0])))
        if is_per_splat and n_splats < a.shape[0]:
            a = a[:n_splats]
        local_spatial = spatial_perms_global
        used_sel: Optional[np.ndarray] = None
        if a.size > args.max_numel and a.ndim >= 1 and a.shape[0] > 1:
            keep = max(int(args.max_numel / max(np.prod(a.shape[1:]), 1)), 1)
            rng = np.random.default_rng(0)
            sel = rng.choice(a.shape[0], size=min(keep, a.shape[0]), replace=False)
            sel.sort()
            a = a[sel].copy()
            used_sel = sel
            if spatial_perms_global and world_means_field is not None:
                local_spatial = build_spatial_permutations(
                    world_means_field.array[sel], a.shape[0])
        field_is_per_splat = (n_splats is not None
                              and len(f.shape) >= 1
                              and f.shape[0] == n_splats)
        local_n_splats = a.shape[0] if field_is_per_splat else None
        msg = f"  analysing {f.name} ({tuple(a.shape)}, {f.kind})"
        if used_sel is not None:
            msg += f"  [subsampled to {a.shape[0]:,} of {f.shape[0]:,}]"
        print(msg)
        t0 = time.time()
        try:
            results = simulate_all_strategies(
                a, kind=f.kind,
                n_splats=local_n_splats,
                spatial_perms=local_spatial,
                max_workers=args.threads)
        except Exception as e:
            print(f"    SKIP {f.name}: {e}")
            continue
        print(f"    {len(results)} strategies in {time.time() - t0:.1f}s "
               f"(j={args.threads})")
        per_field_results.append((f, results))
        field_subsampled_arrays.append(a)

    # --- Render per-field HTML in parallel. Each render generates several
    # matplotlib figures (histograms + heatmaps + strategy bar charts). The
    # Agg backend is safe across threads as long as each figure object is
    # owned by one thread (which is the case here). ---
    if per_field_results:
        print(f"  rendering {len(per_field_results)} per-field sections "
               f"(plot_j={args.plot_threads})")
        t0 = time.time()
        field_html: list[Optional[str]] = [None] * len(per_field_results)
        def _render_one(i, f, res, a):
            return i, render_field(f, res, plot_arr=a)
        with ThreadPoolExecutor(max_workers=args.plot_threads) as ex:
            futs = [ex.submit(_render_one, i, f, res, a)
                     for i, ((f, res), a) in enumerate(
                         zip(per_field_results, field_subsampled_arrays))]
            for fut in as_completed(futs):
                i, html = fut.result()
                field_html[i] = html
        for h in field_html:
            if h is not None:
                out_html.append(h)
        print(f"    rendered in {time.time() - t0:.1f}s")

    out_html.append(render_cross_field(fields))
    if spatial_perms_global:
        out_html.append(render_spatial_ranking(per_field_results))
    if not args.no_pair:
        out_html.append(render_pair_section(
            fields, bc_step=args.bc_step,
            compute_workers=args.threads,
            plot_workers=args.plot_threads,
            max_numel=args.pair_max_numel))
    out_html.append(render_summary(per_field_results))

    out_html.append("</body></html>")

    out_path = Path(args.output) if args.output else (ckpt_path / "quantization_report.html")
    out_path.write_text("\n".join(out_html))
    print(f"\nReport written to: {out_path}")
    print(f"Size: {out_path.stat().st_size / 1024:.1f} KiB")


if __name__ == "__main__":
    main()
