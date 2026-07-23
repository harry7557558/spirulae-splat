#!/usr/bin/env python3
"""Regenerate tests/python/dataparser_golden.json.

The golden values were captured while `modules/dataparser.py` still had its
parsing implementation and `test_dataparser_parity.py` proved the native
parsers agreed with it on every frame, pose, intrinsic, distortion
coefficient, seed point and split index -- for 4 formats x 4 config variants x
both splits. The Python implementation is gone; the golden is what is left of
that proof.

Run this ONLY when a parser change is intentional, and say in the commit
message what moved.

    python tests/python/make_dataparser_golden.py
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

import dataset_fixtures as fixtures
from spirulae_splat.modules.native_dataparser import NativeParserConfig, parse_dataset

# Kept in sync with test_dataparser_parity.py.
FORMAT_WRITERS = {
    "colmap_text": fixtures.write_colmap_text,
    "colmap_binary": fixtures.write_colmap_binary,
    "nerfstudio": fixtures.write_nerfstudio,
    "metashape": fixtures.write_metashape,
}

CONFIG_VARIANTS = {
    "default": {},
    "eval_interval": {"eval_mode": "interval", "eval_interval": 3},
    "eval_fraction": {"eval_mode": "fraction", "train_split_fraction": 0.7},
    "validation": {"validation_fraction": 0.25},
}


def _arr(a):
    """Round-trippable, diffable representation of a float/int array."""
    a = np.asarray(a)
    if a.dtype.kind == "f":
        return [float(f"{v:.7g}") for v in a.reshape(-1).tolist()]
    return [int(v) for v in a.reshape(-1).tolist()]


def describe(ds) -> dict:
    """Everything the parity test used to compare, per parse."""
    return {
        "num_cameras": int(ds.num_cameras),
        "image_names": [Path(f).name for f in ds.image_filenames],
        "widths": _arr(ds.widths),
        "heights": _arr(ds.heights),
        "camera_models": _arr(ds.camera_models),
        "c2w": _arr(ds.c2w),
        "intrins": _arr(ds.intrins),
        "dist_coeffs": _arr(ds.dist_coeffs),
        "train_frame_scale": float(f"{ds.train_frame_scale:.7g}"),
        "train_to_normalized": _arr(ds.train_to_normalized),
        "val_indices": _arr(ds.val_indices),
        # The seed cloud is the same for every variant of a format, so it is
        # stored as a digest rather than 3x N floats.
        "points_xyz_sha256": hashlib.sha256(
            np.ascontiguousarray(ds.points_xyz, dtype=np.float32).tobytes()
        ).hexdigest(),
        "points_rgb_sha256": hashlib.sha256(
            np.ascontiguousarray(ds.points_rgb, dtype=np.uint8).tobytes()
        ).hexdigest(),
        "num_points": int(np.asarray(ds.points_xyz).reshape(-1, 3).shape[0]),
    }


def build(root: Path) -> dict:
    scene = fixtures.make_scene()
    datasets = {name: writer(root / name, scene)
                for name, writer in FORMAT_WRITERS.items()}
    golden = {}
    for fmt in sorted(FORMAT_WRITERS):
        golden[fmt] = {}
        for variant, overrides in sorted(CONFIG_VARIANTS.items()):
            golden[fmt][variant] = {
                split: describe(parse_dataset(
                    datasets[fmt], NativeParserConfig(split=split, **overrides)))
                for split in ("train", "eval")
            }
    return golden


def main():
    with tempfile.TemporaryDirectory() as tmp:
        golden = build(Path(tmp))
    path = Path(__file__).with_name("dataparser_golden.json")
    path.write_text(json.dumps(golden, indent=1, sort_keys=True) + "\n")
    n = len(FORMAT_WRITERS) * len(CONFIG_VARIANTS) * 2
    print(f"wrote {path} ({n} parses, {path.stat().st_size/1024:.0f} KiB)")


if __name__ == "__main__":
    main()
