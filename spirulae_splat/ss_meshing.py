"""Generate a surface mesh from an already-trained 3DGS model.

The heavy lifting (occupancy field, Delaunay, marching tetrahedra, bisection,
vertex merge, PLY writing) is implemented in C++/CUDA (Meshing.{h,cu},
MeshingHost.cpp). This script only handles interfacing: loading the checkpoint
PLY and, optionally, the dataset camera poses.

Usage:
    python -m spirulae_splat.ss_meshing <checkpoint_dir> [--data <dataset_dir>]

`checkpoint_dir` is either a run directory (containing config.json and a
`step-*.ckpt/` folder) or a `*.ckpt` directory directly. The mesh is written to
`mesh.ply` next to the checkpoint's `splat.ply`.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch

import spirulae_splat.csrc as _C


def _find_ckpt(checkpoint_dir: Path):
    """Return (splat_ply_path, run_dir) given a run dir or a *.ckpt dir."""
    if (checkpoint_dir / "splat.ply").is_file():
        return checkpoint_dir / "splat.ply", checkpoint_dir.parent
    ckpts = sorted(checkpoint_dir.glob("step-*.ckpt"))
    if not ckpts:
        ckpts = sorted(checkpoint_dir.glob("*.ckpt"))
    for ck in reversed(ckpts):
        if (ck / "splat.ply").is_file():
            return ck / "splat.ply", checkpoint_dir
    raise FileNotFoundError(f"No splat.ply found under {checkpoint_dir}")


def load_splats(ply_path: Path):
    """Load raw (un-activated) 3DGS parameters from an inria-format PLY."""
    from plyfile import PlyData

    ply = PlyData.read(str(ply_path))
    v = ply["vertex"]
    means = np.stack([v["x"], v["y"], v["z"]], axis=-1).astype(np.float32)
    quats = np.stack([v["rot_0"], v["rot_1"], v["rot_2"], v["rot_3"]], axis=-1).astype(np.float32)
    scales = np.stack([v["scale_0"], v["scale_1"], v["scale_2"]], axis=-1).astype(np.float32)
    opacities = np.asarray(v["opacity"]).astype(np.float32)
    return (
        torch.from_numpy(means),
        torch.from_numpy(quats),
        torch.from_numpy(scales),
        torch.from_numpy(opacities),
    )


def load_camera_positions(run_dir: Path, data_dir: Path):
    """Parse the dataset and return camera centers in the splat coordinate frame."""
    from spirulae_splat.modules.dataparser import (
        SpirulaeSplatDataparser,
        SpirulaeSplatDataParserConfig,
    )

    cfg_path = run_dir / "config.json"
    cfg_json = json.load(open(cfg_path)) if cfg_path.is_file() else {}
    dp_json = cfg_json.get("dataparser", {})

    cfg = SpirulaeSplatDataParserConfig()
    path_fields = {
        "colmap_recon_dir", "image_dir", "mask_dir", "depth_dir", "normal_dir",
        "metashape_xml", "metashape_ply", "metashape_psx",
    }
    for k, val in dp_json.items():
        if not hasattr(cfg, k):
            continue
        if k in path_fields and isinstance(val, str):
            val = Path(val)
        setattr(cfg, k, val)

    parser = SpirulaeSplatDataparser(cfg, data_dir)
    train_out, _eval_out = parser.parse()
    c2w = train_out["cameras"].camera_to_worlds  # (N, 3, 4)
    if c2w.dim() == 2:
        c2w = c2w.unsqueeze(0)
    positions = c2w[:, :3, 3].contiguous().float()  # (N, 3)

    # The splats live in a frame scaled by relative_scale relative to c2w.
    rel = cfg_json.get("model", {}).get("relative_scale", None)
    if rel is not None:
        positions = positions * float(rel)
    return positions


def main():
    ap = argparse.ArgumentParser(description="Generate a mesh from a trained 3DGS model.")
    ap.add_argument("checkpoint", type=Path, help="run dir or *.ckpt dir")
    ap.add_argument("--data", type=Path, default=None,
                    help="dataset dir for camera-based occupancy. Defaults to "
                         "config.json's `data` if present; pass --no-data to force "
                         "the static (density-only) occupancy field.")
    ap.add_argument("--no-data", action="store_true",
                    help="ignore the dataset; use the static density occupancy field")
    ap.add_argument("--output", type=Path, default=None, help="output PLY path")
    ap.add_argument("--iso", type=float, default=0.5)
    ap.add_argument("--merge-factor", type=float, default=1.0,
                    help="local short-edge merge threshold multiplier (0 disables)")
    ap.add_argument("--bisection-iters", type=int, default=20)
    ap.add_argument("--max-cameras", type=int, default=64)
    ap.add_argument("--max-grid-res", type=int, default=512)
    ap.add_argument("--grid-cell-factor", type=float, default=2.0)
    ap.add_argument("--num-threads", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    splat_ply, run_dir = _find_ckpt(args.checkpoint)
    print(f"[meshing] loading {splat_ply}")
    means, quats, scales, opacities = load_splats(splat_ply)
    print(f"[meshing] {means.shape[0]} Gaussians")

    cam_positions = torch.empty((0, 3), dtype=torch.float32)
    if not args.no_data:
        data_dir = args.data
        if data_dir is None:
            cfg_path = run_dir / "config.json"
            if cfg_path.is_file():
                data_str = json.load(open(cfg_path)).get("data", None)
                if data_str is not None:
                    cand = Path(data_str)
                    if not cand.is_absolute():
                        cand = (run_dir / data_str).resolve()
                    data_dir = cand if cand.exists() else None
        if data_dir is not None and Path(data_dir).exists():
            print(f"[meshing] loading cameras from {data_dir}")
            cam_positions = load_camera_positions(run_dir, Path(data_dir))
            print(f"[meshing] {cam_positions.shape[0]} cameras")
        else:
            print("[meshing] no dataset found; using static density occupancy field")

    out_path = args.output if args.output is not None else (splat_ply.parent / "mesh.ply")

    _C.generate_mesh(
        means, quats, scales, opacities, cam_positions, str(out_path),
        iso=args.iso,
        merge_factor=args.merge_factor,
        bisection_iters=args.bisection_iters,
        max_cameras=args.max_cameras,
        max_grid_res=args.max_grid_res,
        grid_cell_factor=args.grid_cell_factor,
        num_threads=args.num_threads,
        verbose=not args.quiet,
    )
    print(f"[meshing] done -> {out_path}")


if __name__ == "__main__":
    main()
