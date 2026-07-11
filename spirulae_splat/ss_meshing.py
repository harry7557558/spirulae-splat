"""Generate a surface mesh from an already-trained 3DGS model.

The heavy lifting (occupancy field, Delaunay, marching tetrahedra, bisection,
vertex merge, PLY writing) is implemented in C++/CUDA (Meshing.{h,cu},
MeshingHost.cpp). This script only handles interfacing: loading the checkpoint
PLY and, optionally, the dataset camera poses.

Usage:
    python -m spirulae_splat.ss_meshing <checkpoint_dir> [--data <dataset_dir>]

`checkpoint_dir` is either a run directory (containing config.json and a
`step-*.ckpt/` folder) or a `*.ckpt` directory directly. The mesh is written
next to the checkpoint's `splat.ply` as `mesh.<ext>` for each requested
--format (ply, obj, gltf, glb), colored per --color (none / vertex / texture).
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
    features_dc = np.stack([v["f_dc_0"], v["f_dc_1"], v["f_dc_2"]], axis=-1).astype(np.float32)
    return (
        torch.from_numpy(means),
        torch.from_numpy(quats),
        torch.from_numpy(scales),
        torch.from_numpy(opacities),
        torch.from_numpy(features_dc),
    )


# Dataparser distortion columns (see dataparser.DISTORTION_KEYS) -> engine [10]
# layout expected by the 3DGUT projection (Camera.h kCameraDistortionParams:
# k1,k2,p1,p2,k3,k4,k5,k6,sx1,sy1). Source order is k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2.
# k5,k6 have no source column (filled with zero); b1,b2 (metashape) are dropped.
# NOTE: only validated on near-zero-distortion scenes (mip360); revisit for a
# heavily distorted dataset.
_DIST_SRC = "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split()
_DIST_DST = ["k1", "k2", "p1", "p2", "k3", "k4", "k5", "k6", "sx1", "sy1"]

_IMG_EXTS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff", ".JPG", ".PNG"}


def _resolve_render_resolution(dp_json, data_dir: Path, native_w, native_h):
    """Pick the render resolution from config's image_dir/rescale_camera_to_fit.

    The nerfstudio parse keys image paths off transforms.json's `file_path`
    (typically full-res `images/`) and ignores `image_dir`, so a bool
    `rescale_camera_to_fit=True` measures the wrong (full-res) image. We instead
    honor `image_dir` directly: read the size of one actual image in
    data_dir/image_dir and scale intrinsics to it. Falls back to the native
    (parsed) resolution when image_dir is unset/missing. A numeric
    `rescale_camera_to_fit` (divide-by-N) is applied when present.

    `native_w`/`native_h` are per-camera int tensors (C,); datasets may mix
    resolutions, so the same downscale ratio (sx, sy) is applied to each camera's
    native size. Returns (widths(C,), heights(C,), sx, sy) where intrinsics scale
    by the scalar (sx, sy).
    """
    native_w = torch.as_tensor(native_w).reshape(-1)
    native_h = torch.as_tensor(native_h).reshape(-1)
    w0, h0 = int(native_w[0]), int(native_h[0])

    return native_w.to(torch.int32), native_h.to(torch.int32), 1.0, 1.0  # TODO

    rescale = dp_json.get("rescale_camera_to_fit", False)
    if not isinstance(rescale, bool) and isinstance(rescale, (int, float)) and rescale > 0:
        rfn = {"floor": math.floor, "ceil": math.ceil, "round": round}.get(
            dp_json.get("downscale_rounding_mode", "floor"), math.floor)
        widths = torch.tensor([int(rfn(int(w) / rescale)) for w in native_w], dtype=torch.int32)
        heights = torch.tensor([int(rfn(int(h) / rescale)) for h in native_h], dtype=torch.int32)
        return widths, heights, 1.0 / rescale, 1.0 / rescale

    image_dir = dp_json.get("image_dir", "images")
    if image_dir:
        d = Path(data_dir) / image_dir
        if d.is_dir():
            sample = next((p for p in sorted(d.iterdir())
                           if p.is_file() and p.suffix in _IMG_EXTS), None)
            if sample is not None:
                from PIL import Image
                with Image.open(sample) as im:
                    tw, th = im.size
                sx, sy = tw / w0, th / h0
                widths = torch.round(native_w.float() * sx).to(torch.int32)
                heights = torch.round(native_h.float() * sy).to(torch.int32)
                return widths, heights, sx, sy
    return native_w.to(torch.int32), native_h.to(torch.int32), 1.0, 1.0


def _camera_type_to_model(camera_type) -> str:
    """Map a nerfstudio/COLMAP CameraType value to an engine camera-model name."""
    from spirulae_splat.modules.camera import CameraType
    val = camera_type[0] if isinstance(camera_type, (list, tuple)) else camera_type
    if hasattr(val, "item"):
        val = val.item()
    if val == CameraType.EQUIDISTANT.value:
        return "FISHEYE"
    if val == CameraType.EQUIRECTANGULAR.value:
        return "EQUIRECTANGULAR"
    return "PINHOLE"


def load_cameras(run_dir: Path, data_dir: Path):
    """Parse the dataset; return camera intrinsics/extrinsics in the splat frame.

    Returns a dict with:
      positions   (C,3)      camera centers (splat frame)
      viewmats    (C,4,4)    world->camera (splat frame)
      intrins     (C,4)      fx, fy, cx, cy
      dist_coeffs (C,10)     engine distortion layout
      width,height (C,) int  per-camera render size
      camera_model str       engine model name
    """
    from spirulae_splat.modules.dataparser import (
        SpirulaeSplatDataparser,
        SpirulaeSplatDataParserConfig,
        DISTORTION_KEYS,
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

    # The bool auto-rescale opens a full-res image per frame (slow) and, on the
    # nerfstudio path, ignores image_dir -- we redo resolution from image_dir
    # below, so disable it here. A positive numeric rescale is left for the
    # parser; the C++ CLI stores "off" as the number 0 (and "auto" as -1),
    # which would otherwise divide intrinsics by zero.
    rcf = getattr(cfg, "rescale_camera_to_fit", False)
    if isinstance(rcf, bool) or (isinstance(rcf, (int, float)) and rcf <= 0):
        cfg.rescale_camera_to_fit = False

    parser = SpirulaeSplatDataparser(cfg, data_dir)
    train_out, _eval_out = parser.parse()
    cameras = train_out["cameras"]

    c2w = cameras.camera_to_worlds  # (C, 3, 4)
    if c2w.dim() == 2:
        c2w = c2w.unsqueeze(0)
    c2w = c2w.float()
    C = c2w.shape[0]

    rel = cfg_json.get("model", {}).get("relative_scale", None)
    rel = float(rel) if rel is not None else 1.0

    # The splats live in a frame scaled by relative_scale relative to c2w; scale
    # the translation, leave rotation alone.
    positions = (c2w[:, :3, 3].contiguous() * rel).float()  # (C, 3)

    # Build (C,4,4) world->camera. c2w is camera->world (OpenGL: y up, z back);
    # the projection expects OpenCV (y down, z forward), so flip the y,z columns
    # of R before inverting -- matching ss_viewer._camera_to_viewmat. The camera
    # translation is scaled by relative_scale (splat frame).
    R = c2w[:, :3, :3] * torch.tensor([1.0, -1.0, -1.0])[None, None, :]
    T = c2w[:, :3, 3:4] * rel
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.matmul(R_inv, T)
    viewmats = torch.eye(4).repeat(C, 1, 1)
    viewmats[:, :3, :3] = R_inv
    viewmats[:, :3, 3:4] = T_inv
    viewmats = viewmats.contiguous().float()  # (C,4,4)

    # intrins (fx, fy, cx, cy)
    intr = cameras.intrins
    if isinstance(intr, (tuple, list)):
        intr = torch.stack([torch.as_tensor(x).reshape(C).float() for x in intr], dim=-1)
    intrins = torch.as_tensor(intr).reshape(C, 4).contiguous().float()

    # distortion -> engine [C,10] layout
    dist_coeffs = torch.zeros(C, 10, dtype=torch.float32)
    dp = getattr(cameras, "distortion_params", None)
    if dp is not None and torch.as_tensor(dp).numel() > 0:
        dp = torch.as_tensor(dp).reshape(C, -1).float()
        src_idx = {k: i for i, k in enumerate(DISTORTION_KEYS)}
        for j, key in enumerate(_DIST_DST):
            if key in src_idx and src_idx[key] < dp.shape[1]:
                dist_coeffs[:, j] = dp[:, src_idx[key]]

    # Per-camera native size (datasets may mix resolutions); broadcast a scalar.
    native_w = torch.as_tensor(cameras.width).reshape(-1).expand(C).contiguous()
    native_h = torch.as_tensor(cameras.height).reshape(-1).expand(C).contiguous()

    # Resolution honoring config's image_dir / rescale_camera_to_fit.
    width, height, sx, sy = _resolve_render_resolution(dp_json, data_dir, native_w, native_h)
    if (sx, sy) != (1.0, 1.0):
        intrins[:, 0] *= sx  # fx
        intrins[:, 1] *= sy  # fy
        intrins[:, 2] *= sx  # cx
        intrins[:, 3] *= sy  # cy
        intrins = intrins.contiguous()

    camera_model = _camera_type_to_model(cameras.camera_type)

    return {
        "positions": positions,
        "viewmats": viewmats,
        "intrins": intrins,
        "dist_coeffs": dist_coeffs,
        "width": width,
        "height": height,
        "camera_model": camera_model,
    }


def entrypoint():
    ap = argparse.ArgumentParser(description="Generate a mesh from a trained 3DGS model.")
    ap.add_argument("checkpoint", type=Path, help="run dir or *.ckpt dir")
    ap.add_argument("--data", type=Path, default=None,
                    help="dataset dir for camera-based occupancy. Defaults to "
                         "config.json's `data` if present; pass --no-data to force "
                         "the static (density-only) occupancy field.")
    ap.add_argument("--no-data", action="store_true",
                    help="ignore the dataset; use the static density occupancy field")
    ap.add_argument("--output", type=Path, default=None,
                    help="output base path (a known mesh extension is stripped; "
                         "one file per requested format is written next to it)")
    ap.add_argument("--format", type=str, default="ply",
                    help="comma-separated output formats: ply, obj, gltf, glb. "
                         "With --color texture a format may carry a texture "
                         "encoding suffix: glb+png (default), glb+jpg "
                         "(JPEG quality 95), glb+jpeg75 (JPEG quality 75)")
    ap.add_argument("--color", type=str, default="vertex",
                    choices=["none", "vertex", "texture"],
                    help="mesh color: none, per-vertex color, or a baked texture "
                         "atlas. PLY+texture and OBJ+vertex are rejected.")
    ap.add_argument("--texture-size", type=int, default=2048,
                    help="texture atlas resolution (square)")
    ap.add_argument("--tex-gutter-px", type=int, default=4,
                    help="atlas spacing between UV charts, in texels")
    ap.add_argument("--chart-angle-deg", type=float, default=60.0,
                    help="max face-normal deviation within a UV chart")
    ap.add_argument("--iso", type=float, default=None,
                    help="isosurface level. Default: 0.5 with cameras (--data), "
                         "0.2 without.")
    ap.add_argument("--merge-factor", type=float, default=1.0,
                    help="local short-edge merge threshold multiplier (0 disables)")
    ap.add_argument("--bisection-iters", type=int, default=3)
    ap.add_argument("--max-cameras", type=int, default=-1)
    ap.add_argument("--max-grid-res", type=int, default=512)
    ap.add_argument("--grid-cell-factor", type=float, default=2.0)
    ap.add_argument("--num-threads", type=int, default=0)
    ap.add_argument("--occ-band", type=float, default=0.2,
                    help="half-width of the accumulated-alpha band around iso used "
                         "for the per-view occupancy crossing depths (render path)")
    ap.add_argument("--carve-k", type=int, default=1,
                    help="aggregate occupancy as the k-th smallest over cameras "
                         "(k=1 = strict min/space-carving; k>1 is robust to a few "
                         "underestimating views)")
    ap.add_argument("--no-cull-unseen", action="store_true",
                    help="keep mesh vertices that no training camera sees "
                         "(default: drop them; only applies with --data)")
    ap.add_argument("--merge-max-flip-deg", type=float, default=60.0,
                    help="reject a merge collapse that would rotate an incident "
                         "face normal past this many degrees (avoids folds / "
                         "slivers); <=0 or >=180 disables the guard")
    ap.add_argument("--floater-min-faces", type=int, default=10,
                    help="drop connected components with fewer than this many "
                         "faces (<=1 disables)")
    ap.add_argument("--fill-hole-ratio", type=float, default=0.05,
                    help="fill a boundary loop when its bounding-box diagonal is "
                         "less than this fraction of its connected component's "
                         "bounding-box diagonal (0 disables this criterion)")
    ap.add_argument("--fill-hole-max-edges", type=int, default=20,
                    help="also fill any boundary loop with at most this many edges "
                         "(so tiny holes always fill; 0 disables this criterion). "
                         "A loop is filled if either criterion is met.")
    ap.add_argument("--quality-iters", type=int, default=3,
                    help="valence-flip + tangential-relaxation iterations "
                         "(improves triangle shape/valence; 0 disables)")
    ap.add_argument("--degenerate-angle-deg", type=float, default=15.0,
                    help="remove triangles whose smallest angle is below this many "
                         "degrees by a fidelity-preserving local edit (collapse a "
                         "needle's short edge, flip a cap's long edge); <=0 disables")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    # Validate the export request before any heavy work. A token is
    # "fmt[+png|+jpg|+jpeg<quality>]"; the C++ side re-validates in depth.
    import re
    tokens = [f.strip().lower().lstrip(".") for f in args.format.split(",") if f.strip()]
    known = {"ply", "obj", "gltf", "glb"}
    formats, bases = [], []
    for tok in tokens:
        m = re.fullmatch(r"(ply|obj|gltf|glb)(?:\+(png|jpe?g(\d+)?))?", tok)
        if not m:
            ap.error(f"--format: invalid token {tok!r} (expected one of "
                     f"{', '.join(sorted(known))}, optionally with +png, "
                     f"+jpg, or +jpeg<quality>)")
        base, enc, quality = m.group(1), m.group(2), m.group(3)
        if base in bases:
            ap.error(f"--format: {base!r} requested more than once")
        if enc and args.color != "texture":
            ap.error(f"--format: {tok!r} has a texture encoding but --color "
                     f"is not texture")
        if quality is not None and not (1 <= int(quality) <= 100):
            ap.error(f"--format: {tok!r}: JPEG quality must be in 1..100")
        bases.append(base)
        formats.append(tok)
    if not formats:
        ap.error("--format: no valid format given")
    if args.color == "texture" and "ply" in bases:
        ap.error("PLY does not support textured meshes; use --color vertex "
                 "or drop ply from --format")
    if args.color == "vertex" and "obj" in bases:
        ap.error("OBJ has no standard per-vertex color; use --color texture "
                 "or drop obj from --format")

    splat_ply, run_dir = _find_ckpt(args.checkpoint)
    print(f"[meshing] loading {splat_ply}")
    means, quats, scales, opacities, features_dc = load_splats(splat_ply)
    print(f"[meshing] {means.shape[0]} Gaussians")

    cam_positions = torch.empty((0, 3), dtype=torch.float32)
    cams = None
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
            cams = load_cameras(run_dir, Path(data_dir))
            cam_positions = cams["positions"]
            cw, chh = cams["width"], cams["height"]
            if int(cw.min()) == int(cw.max()) and int(chh.min()) == int(chh.max()):
                res_str = f"{int(cw[0])}x{int(chh[0])}"
            else:
                res_str = (f"{int(cw.min())}x{int(chh.min())}..{int(cw.max())}x{int(chh.max())} "
                           "per-camera")
            print(f"[meshing] {cam_positions.shape[0]} cameras "
                  f"({res_str}, {cams['camera_model']})")

    using_cameras = cam_positions.shape[0] > 0
    if not using_cameras:
        # Bold yellow recommendation -- the static field is markedly lower quality.
        print("\033[1;33m"
              "[meshing] No camera dataset in use: meshing from Gaussian densities "
              "only.\n"
              "          Pass --data <dataset_dir> to use camera-based occupancy, "
              "which gives\n"
              "          significantly better surfaces."
              "\033[0m")

    # iso default depends on the occupancy field: the camera (visual-hull-like)
    # field is sharper, so a higher level sits on the surface.
    iso = args.iso if args.iso is not None else (0.5 if using_cameras else 0.2)

    out_path = args.output if args.output is not None else (splat_ply.parent / "mesh.ply")

    # Always pass the camera tensors explicitly: pybind cannot fall back on the
    # binding's `at::Tensor()` defaults when these are omitted, so the no-data
    # path must send empty tensors rather than leaving them unset.
    if cams is not None:
        cam_kwargs = dict(
            viewmats=cams["viewmats"],
            intrins=cams["intrins"],
            dist_coeffs=cams["dist_coeffs"],
            cam_widths=cams["width"],
            cam_heights=cams["height"],
            camera_model=cams["camera_model"],
        )
    else:
        empty = torch.empty(0, dtype=torch.float32)
        cam_kwargs = dict(
            viewmats=empty, intrins=empty, dist_coeffs=empty,
            cam_widths=empty, cam_heights=empty,
        )

    _C.generate_mesh(
        means, quats, scales, opacities, features_dc, cam_positions, str(out_path),
        iso=iso,
        merge_factor=args.merge_factor,
        bisection_iters=args.bisection_iters,
        max_cameras=args.max_cameras,
        max_grid_res=args.max_grid_res,
        grid_cell_factor=args.grid_cell_factor,
        num_threads=args.num_threads,
        verbose=not args.quiet,
        carve_k=args.carve_k,
        cull_unseen=not args.no_cull_unseen,
        merge_max_flip_deg=args.merge_max_flip_deg,
        floater_min_faces=args.floater_min_faces,
        fill_hole_ratio=args.fill_hole_ratio,
        fill_hole_max_edges=args.fill_hole_max_edges,
        degenerate_angle_deg=args.degenerate_angle_deg,
        quality_iters=args.quality_iters,
        color_mode=args.color,
        formats=",".join(formats),
        texture_size=args.texture_size,
        tex_gutter_px=args.tex_gutter_px,
        chart_angle_deg=args.chart_angle_deg,
        **cam_kwargs,
    )
    base = out_path.parent / out_path.stem if out_path.suffix.lower() in \
        {".ply", ".obj", ".gltf", ".glb"} else out_path
    print(f"[meshing] done -> {base}.{{{','.join(bases)}}}")


if __name__ == "__main__":
    entrypoint()
