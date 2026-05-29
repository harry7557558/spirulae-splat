#!/usr/bin/env python3
"""Investigate optimizer-state quantization for a spirulae-splat checkpoint.

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


def _maybe_reshape_bilagrid(name: str, arr: np.ndarray, meta: dict) -> np.ndarray:
    """When the on-disk array is a flat uint8 (the case for the *_g{1,2}_q.npy
    quantized bilagrid buffers), re-impose the (N_cam, C, L, H, W) shape using
    metadata so structural plots and axis permutations apply."""
    if arr.ndim != 1:
        return arr
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
        per_cam = C * L * H * W
        if arr.size % per_cam != 0:
            return arr
        N = arr.size // per_cam
        return arr.reshape(N, C, L, H, W)
    except (KeyError, ValueError):
        return arr


def discover_fields(full_dir: Path, meta: Optional[dict] = None) -> list[Field]:
    """Pick optimizer-state and parameter buffers worth analyzing."""
    meta = meta or {}
    candidates = sorted(full_dir.glob("*.npy"))
    # We want: g1_*, g2_*, bilagrid_*_g1, bilagrid_*_g2, ppisp_g1, ppisp_g2.
    # The quantized variants are *_q.npy (uint8) + sibling _quant_bounds.npy.
    # World params (world_*) are kept around for context (param-vs-momentum
    # comparison) but flagged as "param".
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

    # Resolve quantized pairs once: foo_q.npy plus matching bounds.
    # Bounds-name resolution by buffer family (g1 and g2 share one bounds vector
    # in the SH and bilagrid kernels — both halves of the float4 sit in the
    # same file).
    for stem, path in by_name.items():
        if stem.endswith("_q"):
            base = stem[:-2]               # "g1_features_sh_q" -> "g1_features_sh"
            kind = looks_like_momentum(base) or "param"
            if base in ("g1_features_sh", "g2_features_sh"):
                bounds_key = "quant_bounds_sh"
            elif base.endswith("_g1") or base.endswith("_g2"):
                # bilagrid_<type>_g{1,2} -> bilagrid_<type>_quant_bounds
                bounds_key = f"{base[:-3]}_quant_bounds"
            else:
                bounds_key = f"{base}_quant_bounds"
            bounds_path = by_name.get(bounds_key)
            if bounds_path is None:
                print(f"warn: {stem} has no bounds file (expected "
                       f"'{bounds_key}.npy'); loading raw uint8", file=sys.stderr)
            fields.append(Field(name=base, kind=kind, path=path,
                                 is_quantized_on_disk=True,
                                 bounds_path=bounds_path))
        elif looks_like_momentum(stem):
            # Skip if we already added the quantized pair for it.
            if f"{stem}_q" in by_name:
                continue
            fields.append(Field(name=stem, kind=looks_like_momentum(stem),
                                 path=path, is_quantized_on_disk=False))
        elif stem.startswith("world_") or stem in ("accum_buffer", "bias_correction_steps"):
            fields.append(Field(name=stem, kind="param", path=path,
                                 is_quantized_on_disk=False))

    # Load arrays. For quantized-on-disk fields, dequant using the bounds
    # so we have a usable float array for analysis.
    for f in fields:
        a = np.load(f.path)
        if f.is_quantized_on_disk and f.bounds_path is not None:
            bounds = np.load(f.bounds_path)   # [n_blocks, 4] -> (g1lo, g1hi, sqrtg2lo, sqrtg2hi)
            a = _dequant_uint8_blocks(a, bounds, kind=f.kind)
            a = _maybe_reshape_bilagrid(f.name, a, meta)
        f.array = a
        f.shape = tuple(a.shape)
        f.dtype = str(a.dtype)

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


def _block_quantize(blocks: np.ndarray, sqrt_xform: bool) -> tuple[np.ndarray, np.ndarray]:
    """Per-block uint8 quantize. Returns (dequant_blocks, block_range)."""
    if sqrt_xform:
        t = np.sqrt(np.maximum(blocks, 0.0)).astype(np.float64)
    else:
        t = blocks.astype(np.float64)
    lo = t.min(axis=1, keepdims=True)
    hi = t.max(axis=1, keepdims=True)
    rng = np.maximum(hi - lo, 1e-30)
    q = np.clip(np.floor(QUANT_LEVELS * (t - lo) / rng), 0, QUANT_LEVELS - 1)
    deq = lo + rng * ((q + 0.5) / QUANT_LEVELS)
    if sqrt_xform:
        deq = deq * deq
    return deq.astype(np.float32), (hi - lo).squeeze(1).astype(np.float64)


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
    """A reordering recipe applied to a tensor before the 8-bit / 256-cell
    block quantizer. Components compose in this order:
        1. axis0_perm   (per-splat reorder along axis 0)
        2. axis_perm    (transpose of all axes)
        3. sort_within  (sort by |x| within each super-block)
        4. global_sort  (sort all values by |x|; oracle)
    """
    name: str
    extra: str = ""
    axis0_perm: Optional[np.ndarray] = None      # (N,) int permutation of axis 0
    axis_perm: Optional[tuple] = None            # axes permutation tuple
    sort_within: Optional[int] = None            # super-block size
    global_sort: bool = False


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
        # 5D bilagrid (N_cam, C, L, H, W): pick a few semantically meaningful ones.
        labels_perms = []
        if nd == 5:
            labels_perms = [
                ("perm(C,N,L,H,W)", (1, 0, 2, 3, 4)),
                ("perm(L,H,W,C,N)", (2, 3, 4, 1, 0)),
                ("perm(N,L,H,W,C)", (0, 2, 3, 4, 1)),
                ("perm(C,L,H,W,N)", (1, 2, 3, 4, 0)),
            ]
        for lab, p in labels_perms:
            out.append((lab, p))
    return out


def _run_strategy(arr: np.ndarray, sqrt_xform: bool, s: Strategy) -> StrategyResult:
    """Apply Strategy s as a reordering of arr, run the quantizer, invert,
    measure error against the original arr."""
    a = arr
    inv_axis0 = None
    inv_axis = None

    if s.axis0_perm is not None:
        # Restrict perm to the leading n_active rows if shorter than axis 0.
        p = s.axis0_perm
        if p.shape[0] != a.shape[0]:
            # axis0_perm covers only the active prefix — keep tail in place.
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

    if s.global_sort:
        key = np.sqrt(np.maximum(flat, 0)) if sqrt_xform else flat
        order = np.argsort(key, kind="stable")
        flat_re = flat[order]
        blocks = _pad_to_blocks(flat_re)
        deq, ranges = _block_quantize(blocks, sqrt_xform)
        deq_flat = deq.reshape(-1)[:n]
        deq_back = np.empty_like(flat)
        deq_back[order] = deq_flat
    elif s.sort_within is not None:
        sb = s.sort_within
        key = np.sqrt(np.maximum(flat, 0)) if sqrt_xform else flat
        order = np.empty(n, dtype=np.int64)
        nsb = (n + sb - 1) // sb
        for i in range(nsb):
            lo, hi = i * sb, min(i * sb + sb, n)
            order[lo:hi] = lo + np.argsort(key[lo:hi], kind="stable")
        flat_re = flat[order]
        blocks = _pad_to_blocks(flat_re)
        deq, ranges = _block_quantize(blocks, sqrt_xform)
        deq_flat = deq.reshape(-1)[:n]
        deq_back = np.empty_like(flat)
        deq_back[order] = deq_flat
    else:
        blocks = _pad_to_blocks(flat)
        deq, ranges = _block_quantize(blocks, sqrt_xform)
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


def build_strategies(arr_shape: tuple, n_splats: Optional[int],
                     spatial_perms: dict[str, np.ndarray]) -> list[Strategy]:
    """Enumerate every reordering strategy to evaluate on a tensor."""
    strats: list[Strategy] = []

    # --- 1. identity + axis permutations ---
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
        # Only emit combos when there is something to permute.
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
    return strats


def simulate_all_strategies(arr: np.ndarray, kind: str,
                             n_splats: Optional[int] = None,
                             spatial_perms: Optional[dict[str, np.ndarray]] = None,
                             max_workers: int = 1) -> list[StrategyResult]:
    """Run all strategies. For g2-style fields the per-block sqrt transform is
    applied (matching the on-device kernel); reported RMSE is in raw units.
    Parallel across strategies via ThreadPoolExecutor (numpy releases the GIL
    inside the hot inner loops)."""
    sqrt_xform = kind.endswith("g2") or "_g2" in kind
    strategies = build_strategies(arr.shape, n_splats, spatial_perms or {})
    if max_workers > 1 and len(strategies) > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            fut_map = {ex.submit(_run_strategy, arr, sqrt_xform, s): i
                        for i, s in enumerate(strategies)}
            results: list[Optional[StrategyResult]] = [None] * len(strategies)
            for f in as_completed(fut_map):
                results[fut_map[f]] = f.result()
        return [r for r in results if r is not None]
    else:
        return [_run_strategy(arr, sqrt_xform, s) for s in strategies]


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
    """For bilagrid momentum (N_cam, C, L, H, W): channel/spatial views."""
    if arr.ndim != 5:
        return ""
    N, C, L, H, W = arr.shape
    abs_arr = np.abs(arr).astype(np.float64)

    # Per-channel mean and max across batch/spatial.
    per_C_mean = abs_arr.mean(axis=(0, 2, 3, 4))   # [C]
    per_C_max  = abs_arr.max(axis=(0, 2, 3, 4))    # [C]
    # Per-camera mean across C and spatial.
    per_cam_mean = abs_arr.mean(axis=(1, 2, 3, 4)) # [N]
    # Spatial: mean over batch+C, leave (L, H, W).
    spatial = abs_arr.mean(axis=(0, 1))            # [L, H, W]

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


def render_field(f: Field, results: list[StrategyResult]) -> str:
    a = f.array
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
    out_html.append("<h2>Per-buffer analysis</h2>")
    for f in fields:
        if f.kind not in momentum_kinds:
            continue
        a = f.array
        # Restrict to the active splats range (the table may be over-allocated
        # to max_num_splats with garbage past cur_num_splats).
        is_per_splat = (n_splats is not None and a.ndim >= 1
                        and a.shape[0] == int(meta.get("max_num_splats", a.shape[0])))
        if is_per_splat and n_splats < a.shape[0]:
            a = a[:n_splats]
        # Subsample huge tensors along the OUTER axis so structure is preserved.
        local_spatial = spatial_perms_global
        used_sel: Optional[np.ndarray] = None
        if a.size > args.max_numel and a.ndim >= 1 and a.shape[0] > 1:
            keep = max(int(args.max_numel / max(np.prod(a.shape[1:]), 1)), 1)
            rng = np.random.default_rng(0)
            sel = rng.choice(a.shape[0], size=min(keep, a.shape[0]), replace=False)
            sel.sort()
            a = a[sel].copy()
            used_sel = sel
            # Rebuild spatial perms for the subsample (cheap relative to sim).
            if spatial_perms_global and world_means_field is not None:
                local_spatial = build_spatial_permutations(
                    world_means_field.array[sel], a.shape[0])
        local_n_splats = a.shape[0] if (n_splats is not None
                                          and len(f.shape) >= 1
                                          and (f.shape[0] == n_splats
                                                or (used_sel is not None))) else None
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
        out_html.append(render_field(f, results))

    out_html.append(render_cross_field(fields))
    if spatial_perms_global:
        out_html.append(render_spatial_ranking(per_field_results))
    out_html.append(render_summary(per_field_results))

    out_html.append("</body></html>")

    out_path = Path(args.output) if args.output else (ckpt_path / "quantization_report.html")
    out_path.write_text("\n".join(out_html))
    print(f"\nReport written to: {out_path}")
    print(f"Size: {out_path.stat().st_size / 1024:.1f} KiB")


if __name__ == "__main__":
    main()
