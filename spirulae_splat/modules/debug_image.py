"""Visual-debug helpers for the engine's per-step image buffers.

These pull whatever is currently in ``engine().gt.rgb``, ``engine().gt.alpha``,
and ``engine().fwd.renders.rgb`` back to host and dump them as PNG tile grids
so you can actually *see* what the engine is training against -- especially
helpful on the warp-to-pinhole path where the GT is generated on GPU by a
fused warp kernel and never materialized in Python.

Typical usage from anywhere inside the training loop::

    from spirulae_splat.modules.debug_image import dump_engine_step
    dump_engine_step(f"debug/step{step:06d}")

That writes ``debug/step000000_gt_rgb.png``, ``debug/step000000_render_rgb.png``,
and ``debug/step000000_gt_alpha.png`` (the last only if a mask is present),
each tiled across the batch dimension.

Everything here is import-safe (no torch/cuda init side effects) so you can
drop it into a training script without slowing the no-debug path down.
"""

from __future__ import annotations

import math
import os
from typing import Optional, Tuple

import numpy as np

import spirulae_splat.csrc as _C


# ---------------------------------------------------------------------------
# Tile helpers
# ---------------------------------------------------------------------------

def _grid_layout(n: int) -> Tuple[int, int]:
    """Pick a near-square (rows, cols) layout for n tiles."""
    if n <= 0:
        return 0, 0
    cols = int(math.ceil(math.sqrt(n)))
    rows = int(math.ceil(n / cols))
    return rows, cols


def tile_images(
    images: np.ndarray,   # [B, H, W, C] float in [0, 1] or uint8
    pad: int = 2,
    pad_value: float = 0.25,
    annotate: bool = True,
) -> np.ndarray:
    """Arrange a batch of images into a single (R*H, C*W, C) grid.

    Adds a small inter-image padding gutter and optional indices in the top
    left of each tile (rendered with a simple bitmap font so we don't need
    a PIL / matplotlib dep here). Returns a uint8 image ready for saving
    via stbi / PIL / matplotlib.
    """
    if images.ndim == 3:
        images = images[None, ...]
    if images.ndim != 4:
        raise ValueError(f"expected 4-D [B, H, W, C], got {images.shape}")
    B, H, W, C = images.shape
    if B == 0:
        return np.zeros((1, 1, max(C, 3)), dtype=np.uint8)

    # Normalize to uint8 in [0, 255].
    if images.dtype != np.uint8:
        im = images
        if im.dtype != np.float32:
            im = im.astype(np.float32)
        im = np.clip(im, 0.0, 1.0) * 255.0 + 0.5
        im = im.astype(np.uint8)
    else:
        im = images

    rows, cols = _grid_layout(B)
    g_h = rows * H + (rows + 1) * pad
    g_w = cols * W + (cols + 1) * pad
    pad_u8 = int(np.clip(pad_value, 0.0, 1.0) * 255)
    out = np.full((g_h, g_w, max(C, 3)), pad_u8, dtype=np.uint8)

    for b in range(B):
        r = b // cols
        c = b % cols
        y = pad + r * (H + pad)
        x = pad + c * (W + pad)
        tile = im[b]
        if C == 1:
            tile = np.repeat(tile, 3, axis=-1)
        elif C == 4:
            tile = tile[..., :3]
        out[y:y+H, x:x+W, : tile.shape[-1]] = tile

        if annotate:
            _draw_label(out, x + 2, y + 2, str(b))

    return out


# ---------------------------------------------------------------------------
# 3x5 bitmap digits + a few letters, just enough to label tiles. Avoids a
# PIL / matplotlib dep so this module stays cheap to import.
# ---------------------------------------------------------------------------
_FONT = {
    '0': ["111","101","101","101","111"],
    '1': ["010","110","010","010","111"],
    '2': ["111","001","111","100","111"],
    '3': ["111","001","011","001","111"],
    '4': ["101","101","111","001","001"],
    '5': ["111","100","111","001","111"],
    '6': ["111","100","111","101","111"],
    '7': ["111","001","010","100","100"],
    '8': ["111","101","111","101","111"],
    '9': ["111","101","111","001","111"],
    ' ': ["000","000","000","000","000"],
    ':': ["000","010","000","010","000"],
    '-': ["000","000","111","000","000"],
}


def _draw_label(buf: np.ndarray, x: int, y: int, text: str,
                color=(255, 255, 0), bg=(0, 0, 0)) -> None:
    """In-place: stamp a tiny yellow-on-black label at (x, y)."""
    cw, ch = 4, 6
    label_h = ch
    label_w = cw * len(text)
    if y + label_h >= buf.shape[0] or x + label_w >= buf.shape[1]:
        return
    buf[y:y+label_h, x:x+label_w, :3] = bg
    for ci, ch_ in enumerate(text):
        glyph = _FONT.get(ch_, _FONT[' '])
        for gy, row in enumerate(glyph):
            for gx, p in enumerate(row):
                if p == '1':
                    buf[y + gy, x + ci * cw + gx, :3] = color


# ---------------------------------------------------------------------------
# Engine fetch helpers
# ---------------------------------------------------------------------------

def _tv(t):
    """Wrap a numpy ndarray as the (ptr, elem_size, shape) tuple the engine
    bindings expect."""
    if not t.flags.c_contiguous:
        t = np.ascontiguousarray(t)
    return (int(t.ctypes.data), int(t.dtype.itemsize), list(t.shape))


def get_engine_gt_rgb() -> Optional[np.ndarray]:
    """Pull engine().gt.rgb to host as [B, H, W, 3] float32. None if empty."""
    B, H, W, C = _C.engine_get_gt_rgb_shape()
    if B == 0:
        return None
    buf = np.zeros((B, H, W, C), dtype=np.float32)
    _C.engine_copy_gt_rgb_to_host(_tv(buf))
    return buf


def get_engine_gt_alpha() -> Optional[np.ndarray]:
    """Pull engine().gt.alpha to host as [B, H, W, 1] bool. None if empty."""
    B, H, W, C = _C.engine_get_gt_alpha_shape()
    if B == 0:
        return None
    buf = np.zeros((B, H, W, C), dtype=bool)
    _C.engine_copy_gt_alpha_to_host(_tv(buf))
    return buf


def get_engine_render_rgb() -> Optional[np.ndarray]:
    """Pull engine().fwd.renders.rgb to host as [B, H, W, 3] float32."""
    B, H, W, C = _C.engine_get_render_rgb_shape()
    if B == 0:
        return None
    rgb = np.zeros((B, H, W, C), dtype=np.float32)
    depth = np.zeros((B, H, W, 1), dtype=np.float32)  # unused, but the
    Ts    = np.zeros((B, H, W, 1), dtype=np.float32)  # binding requires args
    _C.engine_copy_render_to_host(_tv(rgb), _tv(depth), _tv(Ts), (0, 0, []))
    return rgb


# ---------------------------------------------------------------------------
# Public: dump one or more buffers as PNG tiles
# ---------------------------------------------------------------------------

def _save_png(path: str, arr_u8: np.ndarray) -> None:
    """Save a uint8 (H, W, C) array as PNG. Tries pillow first then cv2."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    try:
        from PIL import Image
        if arr_u8.shape[-1] == 1:
            Image.fromarray(arr_u8[..., 0], mode="L").save(path)
        else:
            Image.fromarray(arr_u8[..., :3], mode="RGB").save(path)
        return
    except ImportError:
        pass
    try:
        import cv2
        if arr_u8.shape[-1] >= 3:
            cv2.imwrite(path, cv2.cvtColor(arr_u8[..., :3], cv2.COLOR_RGB2BGR))
        else:
            cv2.imwrite(path, arr_u8[..., 0])
        return
    except ImportError:
        pass
    raise RuntimeError("debug_image: need PIL or cv2 installed to save PNGs")


def dump_engine_gt(path: str, **tile_kwargs) -> Optional[np.ndarray]:
    """Save the current engine GT RGB batch as a tile grid PNG."""
    rgb = get_engine_gt_rgb()
    if rgb is None:
        return None
    tile = tile_images(rgb, **tile_kwargs)
    _save_png(path, tile)
    return rgb


def dump_engine_render(path: str, **tile_kwargs) -> Optional[np.ndarray]:
    """Save the current engine RENDER RGB batch as a tile grid PNG."""
    rgb = get_engine_render_rgb()
    if rgb is None:
        return None
    tile = tile_images(rgb, **tile_kwargs)
    _save_png(path, tile)
    return rgb


def dump_engine_gt_alpha(path: str, **tile_kwargs) -> Optional[np.ndarray]:
    """Save the current engine GT mask batch as a tile grid PNG."""
    a = get_engine_gt_alpha()
    if a is None:
        return None
    af = a.astype(np.float32)
    tile = tile_images(af, **tile_kwargs)
    _save_png(path, tile)
    return a


def dump_engine_step(prefix: str,
                     gt: bool = True,
                     render: bool = True,
                     mask: bool = True,
                     **tile_kwargs) -> dict:
    """Convenience: dump all live buffers at once.

    Writes ``{prefix}_gt_rgb.png``, ``{prefix}_render_rgb.png``, and
    ``{prefix}_gt_alpha.png``. Returns a dict mapping suffix -> array
    actually pulled (for further programmatic inspection).
    """
    out = {}
    if gt:
        a = dump_engine_gt(f"{prefix}_gt_rgb.png", **tile_kwargs)
        if a is not None: out["gt_rgb"] = a
    if render:
        a = dump_engine_render(f"{prefix}_render_rgb.png", **tile_kwargs)
        if a is not None: out["render_rgb"] = a
    if mask:
        a = dump_engine_gt_alpha(f"{prefix}_gt_alpha.png", **tile_kwargs)
        if a is not None: out["gt_alpha"] = a
    return out


# ---------------------------------------------------------------------------
# Stats helper: catches the common "buffer is uniformly gray" failure mode
# ---------------------------------------------------------------------------

def buffer_stats(arr: np.ndarray) -> str:
    """One-line summary of a 4-D image buffer's value distribution. Useful
    to confirm a 'mess' is actually uniform gray vs structured noise vs ok."""
    if arr is None or arr.size == 0:
        return "(empty)"
    a = arr.astype(np.float32)
    return (f"shape={tuple(arr.shape)} dtype={arr.dtype} "
            f"min={a.min():.4f} max={a.max():.4f} "
            f"mean={a.mean():.4f} std={a.std():.4f}")
