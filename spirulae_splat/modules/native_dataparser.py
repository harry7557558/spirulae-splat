"""Adapter over the native dataset parsers (src/data/DatasetParser.h).

The C++ parsers are the ones the CLI trainer, the GUI and the WASM viewer
already use. This module is the Python client for them, and is the intended
replacement for the parsing implementation that used to live in
`dataparser.py` + `colmap_utils.py` + `metashape_utils.py` + `camera_utils.py`.
What became of each: the config dataclass stayed in `dataparser.py`; the COLMAP
and Metashape readers moved to `scripts/`, which still parses in Python for
preprocessing; and `camera_utils.py` stayed put but is now on no code path --
it is the reference implementation for the `orientation_method` /
`center_method` options the native parser does not implement yet (see
docs/notes/pose-normalization.md).

`parse_dataset()` imports neither nerfstudio nor torch. `to_dataparser_outputs()`
does import torch, because the dict it builds is made of torch tensors -- it is
the compatibility shim that keeps `dataset.py` / `datamanager.py` (the Python
eval pass's image loading) and `model.py` working off a native parse.

See docs/restructure-proposal.md §4.1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np


# Fields the native config accepts, in the order DatasetParserConfig declares
# them. Names match SpirulaeSplatDataParserConfig where the two overlap.
_STR_PATH_FIELDS = (
    "recon_dir", "image_dir", "mask_dir", "depth_dir", "normal_dir",
    "metashape_xml", "metashape_ply", "metashape_psx",
)


def _native():
    from spirulae_splat.splat.cuda import _C
    return _C


@dataclass
class NativeParserConfig:
    """Mirror of the C++ DatasetParserConfig, with Python-friendly types.

    Every field has the same name and default as its C++ counterpart, so
    `to_native()` is a copy rather than a translation.
    """

    data_format: Optional[Literal["colmap", "nerfstudio", "metashape"]] = None
    """None = auto-detect (nerfstudio, then COLMAP, then Metashape)."""

    split: Literal["train", "eval"] = "train"
    """Which side of the eval_mode split to return. Everything derived from
        the camera set (train_frame_scale, train_to_normalized, the outlier
        filter) is computed over ALL frames before the split, so the two
        parses agree frame-for-frame. eval_mode="all" gives the full set on
        both sides."""

    recon_dir: str = ""
    """COLMAP reconstruction dir relative to the dataset dir; "" = auto."""

    image_dir: str = "images"
    mask_dir: str = "masks"
    depth_dir: str = "depths"
    normal_dir: str = "normals"

    require_image_files: bool = True
    """False keeps frames whose image file is missing (cameras-only parse)."""

    validation_fraction: float = 0.0
    eval_mode: Literal["all", "fraction", "filename", "interval"] = "all"
    train_split_fraction: float = 0.9
    eval_interval: int = 8
    outlier_threshold: float = float("inf")

    rescale_camera_to_fit: float = 0.0
    """Divide stored intrinsics by this factor. 0 = off."""

    downscale_rounding_mode: Literal["floor", "ceil", "round"] = "floor"

    metashape_xml: str = ""
    metashape_ply: str = ""
    metashape_psx: str = ""
    metashape_component: int = -1

    def to_native(self):
        cfg = _native().DatasetParserConfig()
        for f in _STR_PATH_FIELDS:
            setattr(cfg, f, str(getattr(self, f) or ""))
        cfg.require_image_files = bool(self.require_image_files)
        cfg.validation_fraction = float(self.validation_fraction)
        cfg.eval_mode = str(self.eval_mode)
        cfg.train_split_fraction = float(self.train_split_fraction)
        cfg.eval_interval = int(self.eval_interval)
        cfg.outlier_threshold = float(self.outlier_threshold)
        cfg.rescale_camera_to_fit = float(self.rescale_camera_to_fit)
        cfg.downscale_rounding_mode = str(self.downscale_rounding_mode)
        cfg.metashape_component = int(self.metashape_component)
        cfg.split = str(self.split)
        return cfg


def to_native_parser_config(cfg, split: str = "train") -> NativeParserConfig:
    """`SpirulaeSplatDataParserConfig` -> `NativeParserConfig`.

    Field names match except `colmap_recon_dir` -> `recon_dir`.
    Deliberately dropped, because the native parsers do not implement them and
    `TrainerSession.check_config()` already rejects or warns about each:
    scene_scale, orientation_method, center_method, auto_scale_poses,
    train_frame (only "points"). depth_unit_scale_factor is a load-time
    scaling applied by dataset.py, not a parse-time one.
    """
    return NativeParserConfig(
        data_format=cfg.data_format,
        recon_dir=cfg.colmap_recon_dir or "",
        image_dir=cfg.image_dir,
        mask_dir=cfg.mask_dir,
        depth_dir=cfg.depth_dir,
        normal_dir=cfg.normal_dir,
        validation_fraction=cfg.validation_fraction,
        eval_mode=cfg.eval_mode,
        train_split_fraction=cfg.train_split_fraction,
        eval_interval=cfg.eval_interval,
        outlier_threshold=cfg.outlier_threshold,
        rescale_camera_to_fit=_rescale_to_native(cfg.rescale_camera_to_fit),
        downscale_rounding_mode=cfg.downscale_rounding_mode,
        metashape_xml=cfg.metashape_xml or "",
        metashape_ply=cfg.metashape_ply or "",
        metashape_psx=cfg.metashape_psx or "",
        split=split,
    )


def _rescale_to_native(v) -> float:
    """Union[bool, int] -> the native float (0 = off, -1 = auto-detect)."""
    if v is None or v is False:
        return 0.0
    if v is True:
        return -1.0
    return float(v)


@dataclass
class ParsedDataset:
    """Plain-numpy view of the native ParsedDataset.

    Cameras are PER-INPUT, sorted by image filename, in the raw
    `train_frame="points"` frame -- the same convention
    `modules/dataparser.py` produces.
    """

    num_cameras: int
    camera_models: np.ndarray          # [N] int32, CameraModelType
    image_filenames: List[str]
    mask_filenames: List[str]
    depth_filenames: List[str]
    normal_filenames: List[str]
    widths: np.ndarray                 # [N] int32
    heights: np.ndarray                # [N] int32
    c2w: np.ndarray                    # [N, 3, 4] float32, OpenGL convention
    intrins: np.ndarray                # [N, 4] fx, fy, cx, cy
    dist_coeffs: np.ndarray            # [N, 10] k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2
    train_indices: np.ndarray          # [.] int32
    val_indices: np.ndarray            # [.] int32
    points_xyz: np.ndarray             # [P, 3] float32
    points_rgb: np.ndarray             # [P, 3] uint8
    train_frame_scale: float
    train_to_normalized: np.ndarray    # [4, 4] float32, row-major

    @property
    def fx(self) -> np.ndarray: return self.intrins[:, 0]

    @property
    def fy(self) -> np.ndarray: return self.intrins[:, 1]

    @property
    def cx(self) -> np.ndarray: return self.intrins[:, 2]

    @property
    def cy(self) -> np.ndarray: return self.intrins[:, 3]


def _from_native(ds) -> ParsedDataset:
    return ParsedDataset(
        num_cameras=int(ds.num_cameras),
        camera_models=np.asarray(ds.camera_models),
        image_filenames=list(ds.image_filenames),
        mask_filenames=list(ds.mask_filenames),
        depth_filenames=list(ds.depth_filenames),
        normal_filenames=list(ds.normal_filenames),
        widths=np.asarray(ds.widths),
        heights=np.asarray(ds.heights),
        c2w=np.asarray(ds.c2w),
        intrins=np.asarray(ds.intrins),
        dist_coeffs=np.asarray(ds.dist_coeffs),
        train_indices=np.asarray(ds.train_indices),
        val_indices=np.asarray(ds.val_indices),
        points_xyz=np.asarray(ds.points_xyz),
        points_rgb=np.asarray(ds.points_rgb),
        train_frame_scale=float(ds.train_frame_scale),
        train_to_normalized=np.asarray(ds.train_to_normalized),
    )


def parse_dataset(
    dataset_dir: Union[str, Path],
    config: Optional[NativeParserConfig] = None,
) -> ParsedDataset:
    """Parse a COLMAP / Nerfstudio / Metashape dataset with the C++ parsers."""
    config = config or NativeParserConfig()
    fmt = config.data_format or ""
    ds = _native().parse_dataset(str(dataset_dir), config.to_native(), fmt)
    return _from_native(ds)


# ---------------------------------------------------------------------------
# dataparser_outputs adapter
# ---------------------------------------------------------------------------
# The shape `modules/dataparser.py` used to return. Still needed because
# `dataset.py` / `datamanager.py` load images for the Python eval pass, and
# `model.py` takes cameras + a seed cloud. Training does NOT go through this:
# it drives the native TrainerSession, which parses on its own.

_MODEL_NAMES = ("PINHOLE", "FISHEYE", "EQUISOLID", "EQUIRECTANGULAR")


def _model_id_to_name():
    _C = _native()
    return {int(_C.camera_model_from_name(n)): n for n in _MODEL_NAMES}


def to_dataparser_outputs(ds: ParsedDataset, depth_unit_scale_factor: float = 0.001,
                          train_frame: str = "points") -> dict:
    """`ParsedDataset` -> the dict `SpirulaeSplatDataset` / the model consume."""
    import torch
    from spirulae_splat.modules.camera import Cameras

    id_to_name = _model_id_to_name()
    n = ds.num_cameras

    intrins = torch.from_numpy(np.ascontiguousarray(ds.intrins)).float()
    cameras = Cameras(
        intrins=intrins,
        distortion_params=torch.from_numpy(
            np.ascontiguousarray(ds.dist_coeffs)).float(),
        height=torch.from_numpy(np.ascontiguousarray(ds.heights)).to(torch.int32),
        width=torch.from_numpy(np.ascontiguousarray(ds.widths)).to(torch.int32),
        camera_to_worlds=torch.from_numpy(np.ascontiguousarray(ds.c2w)).float(),
        camera_type=[id_to_name[int(m)] for m in ds.camera_models],
        metadata={},
    )

    # dataparser_transforms.json compatibility. The native parser reports the
    # train->normalized similarity as one 4x4; the legacy pair splits it into a
    # unit-rotation transform plus a scalar. Only ever consumed by external
    # tooling -- nothing in this repo reads the file back.
    T = np.asarray(ds.train_to_normalized, dtype=np.float64).reshape(4, 4)
    scale = float(np.linalg.norm(T[:3, 0])) or 1.0
    transform = torch.from_numpy(T[:3, :] / scale).float()

    def _or_none(seq):
        seq = [str(s) for s in seq]
        return seq if any(seq) else None

    return dict(
        cameras=cameras,
        image_filenames=[str(s) for s in ds.image_filenames],
        mask_filenames=_or_none(ds.mask_filenames),
        dataparser_scale=scale,
        dataparser_transform=transform,
        train_frame=train_frame,
        train_frame_scale=float(ds.train_frame_scale),
        train_to_normalized_transform=torch.from_numpy(T).float(),
        metadata={
            "depth_filenames": _or_none(ds.depth_filenames),
            "depth_unit_scale_factor": float(depth_unit_scale_factor),
            "normal_filenames": _or_none(ds.normal_filenames),
            "points3D_xyz": torch.from_numpy(
                np.ascontiguousarray(ds.points_xyz)).float(),
            "points3D_rgb": torch.from_numpy(
                np.ascontiguousarray(ds.points_rgb)).to(torch.uint8),
            "val_indices": [int(i) for i in ds.val_indices],
        },
    )
