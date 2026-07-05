"""Adapt a checkpoint's saved buffers to a different target layout so training
can resume with fewer Gaussians, a different SH degree, or with/without
bilagrid/PPISP. Produces a new state.tar at the TARGET layout which the
(unchanged, exact-match) C++ engine_load_checkpoint then loads.

All host / CPU, buffer-at-a-time (bounded host RAM, no extra VRAM), and
quantization-aware (decode -> transform -> re-encode via resume_codecs). Heavy
work is numpy-vectorized. A big motivation is reducing parameters after an OOM,
so the common path is: same dataset, smaller cap_max / SH degree, maybe dropping
appearance modules.

Splat reduction policy (target = min(cur_ckpt, cap_new)):
  1. truncate the tail (newest splats) by up to the checkpoint's "unsaturated
     slack" (max_ckpt - cur_ckpt), then
  2. drop the remaining count by lowest opacity.
"""

import io
import json
import os
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from spirulae_splat.modules import resume_codecs as codecs
from spirulae_splat.modules.resume import read_state_json

BLOCK = 256

# Per-splat, non-SH attributes and their channel counts.
_ATTR_CH = {"means": 3, "quats": 4, "scales": 3, "opacities": 1, "features_dc": 3}


# --- reduction ---------------------------------------------------------------

def keep_indices(opac, cur, max_ckpt, target):
    """Indices (ascending) of the `target` splats to keep out of the first `cur`.
    Tail-first (up to the unsaturated slack), then lowest-opacity."""
    if target >= cur:
        return np.arange(cur, dtype=np.int64)
    d = int(cur - target)
    keep = np.ones(cur, dtype=bool)
    n_tail = min(d, max(0, int(max_ckpt) - int(cur)))   # unsaturated slack
    if n_tail > 0:
        keep[cur - n_tail:] = False
    d -= n_tail
    if d > 0:
        idx = np.where(keep)[0]
        order = np.argsort(opac[idx], kind="stable")    # ascending opacity
        keep[idx[order[:d]]] = False
    out = np.where(keep)[0]
    assert out.shape[0] == target, (out.shape[0], target)
    return out


# --- per-splat gather (values already dequantized or plain float) -----------

def _gather_rows(a2d, keep, max_new):
    """a2d: [max_ckpt, ch] -> [max_new, ch] with rows[keep] placed at the front."""
    out = np.zeros((max_new,) + a2d.shape[1:], dtype=a2d.dtype)
    out[: keep.shape[0]] = a2d[keep]
    return out


def _resample_sh(a3d, keep, K_new, max_new):
    """a3d: [max_ckpt, K_ckpt, 3] -> [max_new, K_new, 3], band-clamped/zero-padded."""
    Kc = min(a3d.shape[1], K_new)
    out = np.zeros((max_new, K_new, 3), dtype=a3d.dtype)
    out[: keep.shape[0], :Kc] = a3d[keep][:, :Kc]
    return out


# --- quantized per-splat / SH reindex ---------------------------------------

def _adapt_quant(kind, bits, packed, bounds, keep, max_ckpt, max_new,
                 cpp_ckpt, cpp_new, is_sh, K_ckpt, K_new, fpbo):
    """Decode -> reshape per splat -> gather (+ SH band resample) -> re-encode.
    kind: 'adam' | 'linear' | 'log'. Returns (packed_new, bounds_new_flat)."""
    n_ckpt = int(max_ckpt) * cpp_ckpt
    bbc_ckpt = (BLOCK * cpp_ckpt) if fpbo else BLOCK

    def _reindex(v):
        if is_sh:
            a = v.reshape(max_ckpt, K_ckpt, 3)
            return _resample_sh(a, keep, K_new, max_new).reshape(-1)
        a = v.reshape(max_ckpt, cpp_ckpt)
        return _gather_rows(a, keep, max_new).reshape(-1)

    n_new = int(max_new) * cpp_new
    bbc_new = (BLOCK * cpp_new) if fpbo else BLOCK

    if kind == "adam":
        g1, g2 = codecs.decode_adam(packed, bounds, n_ckpt, bits, bbc_ckpt)
        g1n, g2n = _reindex(g1), _reindex(g2)
        p, b = codecs.encode_adam(g1n, g2n, bits, bbc_new)
    elif kind == "linear":
        v = codecs.decode_linear(packed, bounds, n_ckpt, bits, bbc_ckpt)
        p, b = codecs.encode_linear(_reindex(v), bits, bbc_new)
    elif kind == "log":
        v = codecs.decode_log(packed, bounds, n_ckpt, bits, bbc_ckpt)
        p, b = codecs.encode_log(_reindex(v), bits, bbc_new)
    else:
        raise ValueError(kind)
    return p, b.reshape(-1)


# --- per-splat slot classification ------------------------------------------

def classify(name):
    """Describe a per-splat slot's transform, or None if not per-splat.
    Returns dict: {sh: bool, ch|cpp, quant: None|('adam'|'linear',bits,fpbo)}."""
    # plain float per-splat (world + Adam fp32 moments)
    for pfx in ("world.", "eng.g1_", "eng.g2_"):
        if name.startswith(pfx):
            attr = name[len(pfx):]
            if attr == "features_sh":
                return {"sh": True, "dtype": "<f4"}
            if attr in _ATTR_CH:
                return {"sh": False, "ch": _ATTR_CH[attr], "dtype": "<f4"}
    if name == "eng.radii":
        return {"sh": False, "ch": 1, "dtype": "<f4"}
    if name == "eng.accum_buffer":
        return {"sh": False, "ch": 2, "dtype": "<f4"}
    if name == "eng.bias_correction_steps":
        return {"sh": False, "ch": 1, "dtype": "<i4"}
    # quantized Adam moments, FPBO (means/quats/.../features_dc)
    if name.endswith("_qfpbo"):
        attr = name[len("eng."):-len("_qfpbo")]
        if attr in _ATTR_CH:
            return {"sh": False, "cpp": _ATTR_CH[attr], "quant": ("adam", 16, True)}
    # quantized SH Adam state (cell-block / FPBO)
    if name in ("eng.sh_quant", "eng.sh_quant_fpbo"):
        return {"sh": True, "quant": ("adam", 8, name.endswith("_fpbo"))}
    # quantized SH VALUE store
    if name.startswith("eng.world.sh_vq"):
        bits = 16 if "16" in name else 8
        return {"sh": True, "quant": ("linear", bits, name.endswith("_fpbo"))}
    return None


# --- bilagrid spatial resample ----------------------------------------------

def _bilagrid_C(t, state):
    if t == "rgb":
        return int(state.get("bilagrid_rgb", {}).get("C", 12))
    return 2 if t == "depth" else 3            # depth C=2, normal C=3


def _resample_grid_float(flat, n_grids, lhw_ck, C, lhw_new):
    from scipy.ndimage import zoom
    L, H, W = lhw_ck; L2, H2, W2 = lhw_new
    g = flat.reshape(n_grids, L, H, W, C).astype(np.float64)
    z = zoom(g, (1, L2 / L, H2 / H, W2 / W, 1), order=1, mode="nearest")
    return z.astype(np.float32).reshape(-1)


# --- main --------------------------------------------------------------------

def _read_tar(ckpt_dir):
    arrs = {}
    with tarfile.open(Path(ckpt_dir) / "state.tar") as tf:
        for m in tf.getmembers():
            if m.name.endswith(".npy"):
                arrs[m.name[:-4]] = np.load(io.BytesIO(tf.extractfile(m).read()))
    return arrs


def _group_logical(arrs):
    logical = {}
    for name, a in arrs.items():
        if name.endswith(".qb"):
            logical.setdefault(name[:-3], {})["qb"] = a
        elif name.endswith(".q"):
            logical.setdefault(name[:-2], {})["q"] = a
        else:
            logical.setdefault(name, {})["plain"] = a
    return logical


def _add(tf, name, data):
    ti = tarfile.TarInfo(name); ti.size = len(data)
    tf.addfile(ti, io.BytesIO(data))


def _npy_bytes(arr):
    buf = io.BytesIO(); np.save(buf, np.ascontiguousarray(arr)); return buf.getvalue()


def needs_adapt(state, target):
    if int(state["max_num_splats"]) != int(target["max_num_splats"]):
        return True
    if int(state["num_sh"]) != int(target["num_sh"]):
        return True
    if int(state["cur_num_splats"]) > int(target["max_num_splats"]):
        return True
    for t in ("rgb", "depth", "normal"):
        bg = state.get(f"bilagrid_{t}", {})
        ck_on = bool(bg.get("enabled"))
        tg = target["bilagrid"].get(t)
        if ck_on != (tg is not None):
            return True
        if ck_on and tg is not None and (int(bg["L"]), int(bg["H"]), int(bg["W"])) != tuple(tg):
            return True
    if bool(state.get("ppisp", {}).get("enabled")) != bool(target["ppisp"]):
        return True
    return False


def adapt_checkpoint(ckpt_dir, target, out_dir):
    """Transform the checkpoint at `ckpt_dir` to `target` layout and write an
    adapted state.tar under `out_dir`. Returns out_dir, or None if no adaptation
    is needed (caller loads the original directly)."""
    state = read_state_json(ckpt_dir)
    if not needs_adapt(state, target):
        return None

    arrs = _read_tar(ckpt_dir)
    logical = _group_logical(arrs)

    max_ck = int(state["max_num_splats"]); cur_ck = int(state["cur_num_splats"])
    K_ck = int(state["num_sh"])
    max_new = int(target["max_num_splats"]); K_new = int(target["num_sh"])
    n_img = int(target["num_images"])

    # keep-indices from opacities (per-splat reduction)
    target_cur = min(cur_ck, max_new)
    if "world.opacities" not in logical:
        raise RuntimeError("adapt_checkpoint: world.opacities missing; a full "
                           "(save_full_checkpoint) checkpoint is required to resume "
                           "with a different number of Gaussians.")
    opac = logical["world.opacities"]["plain"].reshape(max_ck)[:cur_ck]
    keep = keep_indices(opac, cur_ck, max_ck, target_cur)

    emit = {}
    for base, parts in logical.items():
        cls = classify(base)
        if cls is not None and "quant" not in cls:                 # per-splat float
            a = parts["plain"]
            if cls["sh"]:
                out = _resample_sh(a.reshape(max_ck, K_ck, 3), keep, K_new, max_new)
            else:
                out = _gather_rows(a.reshape(max_ck, cls["ch"]), keep, max_new)
            emit[base] = out.reshape(-1)
        elif cls is not None:                                       # per-splat quant
            kindq, bits, fpbo = cls["quant"]
            if cls["sh"]:
                pk, bd = _adapt_quant(kindq, bits, parts["q"], parts["qb"], keep,
                                      max_ck, max_new, 3 * K_ck, 3 * K_new, True,
                                      K_ck, K_new, fpbo)
            else:
                cpp = cls["cpp"]
                pk, bd = _adapt_quant(kindq, bits, parts["q"], parts["qb"], keep,
                                      max_ck, max_new, cpp, cpp, False, K_ck, K_new, fpbo)
            emit[base + ".q"] = pk; emit[base + ".qb"] = bd
        elif base.startswith("eng.bg.") and not base.startswith("eng.bg_sky."):
            _adapt_bilagrid(base, parts, state, target, n_img, emit)
        elif base.startswith("eng.ppisp."):
            _adapt_ppisp(base, parts, state, target, n_img, emit)
        else:                                                       # bg_sky / color* fixed
            _copy_logical(base, parts, emit)

    st = dict(state)
    st["cur_num_splats"] = int(target_cur)
    st["max_num_splats"] = int(max_new)
    st["num_sh"] = int(K_new)

    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(out / "state.tar", "w") as tf:
        _add(tf, "state.json", json.dumps(st, indent=2).encode())
        for name, arr in emit.items():
            _add(tf, name + ".npy", _npy_bytes(arr))
    return str(out)


def _copy_logical(base, parts, emit):
    if "plain" in parts:
        emit[base] = parts["plain"]
    else:
        emit[base + ".q"] = parts["q"]; emit[base + ".qb"] = parts["qb"]


def _adapt_bilagrid(base, parts, state, target, n_img, emit):
    comps = base.split(".")                 # eng bg t field...
    t = comps[2]
    field = ".".join(comps[3:])
    tgt = target["bilagrid"].get(t)
    if tgt is None:
        return                              # config disabled this channel -> drop
    bg = state[f"bilagrid_{t}"]
    lhw_ck = (int(bg["L"]), int(bg["H"]), int(bg["W"]))
    C = _bilagrid_C(t, state)
    same = tuple(tgt) == lhw_ck

    if field == "scalars":                  # depth per-camera; per-image
        emit[base] = parts["plain"]; return
    if field == "grids":
        flat = parts["plain"]
        n_grids = flat.size // (lhw_ck[0] * lhw_ck[1] * lhw_ck[2] * C)
        if n_grids != n_img:
            raise RuntimeError(_img_err("bilagrid", n_grids, n_img))
        emit[base] = flat if same else _resample_grid_float(flat, n_grids, lhw_ck, C, tuple(tgt))
        return
    if field == "grids_q":                  # value-quant grids (linear 16-bit, cell-block)
        bits = 16 if int(bg.get("value_bits", 16)) == 16 else 8
        n_cells = lhw_ck[0] * lhw_ck[1] * lhw_ck[2] * C
        # derive n_grids from packed length
        n_grids = parts["q"].size // (n_cells * (2 if bits == 16 else 1))
        if n_grids != n_img:
            raise RuntimeError(_img_err("bilagrid", n_grids, n_img))
        if same:
            emit[base + ".q"] = parts["q"]; emit[base + ".qb"] = parts["qb"]; return
        v = codecs.decode_linear(parts["q"], parts["qb"], n_grids * n_cells, bits, BLOCK)
        vr = _resample_grid_float(v.astype(np.float32), n_grids, lhw_ck, C, tuple(tgt))
        pk, bd = codecs.encode_linear(vr, bits, BLOCK)
        emit[base + ".q"] = pk; emit[base + ".qb"] = bd.reshape(-1); return
    # optimizer state (g1/g2/accum/bg_ag/bg_quant): keep only if resolution
    # unchanged; otherwise reset (omit -> engine re-inits to zero).
    if same:
        _copy_logical(base, parts, emit)


def _adapt_ppisp(base, parts, state, target, n_img, emit):
    if not target["ppisp"]:
        return                              # config disabled -> drop
    P = int(state["ppisp"]["num_params"])
    flat = parts["plain"]
    n_cam = flat.size // P
    if n_cam != n_img:
        raise RuntimeError(_img_err("PPISP", n_cam, n_img))
    emit[base] = flat                       # per-image table, direct copy


def _img_err(what, ck, tgt):
    return (f"resume: {what} image count changed ({ck} -> {tgt}). The dataset, "
            f"train/eval split, or warp_to_pinhole appears to have changed in a way "
            f"that would require rewriting per-image {what} tables. This is not "
            f"supported (TODO); resume on the original dataset/split, or disable "
            f"{what} for this run.")
