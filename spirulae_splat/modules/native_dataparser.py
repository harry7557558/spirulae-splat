"""Adapter over the native dataset parsers (src/data/DatasetParser.h).

The C++ parsers are the ones the CLI trainer, the GUI and the WASM viewer
already use. This module is the Python client for them, and is the intended
replacement for `dataparser.py` + `colmap_utils.py` + `metashape_utils.py` and
the parsing half of `camera_utils.py` / `dataset.py`.

Nothing here imports nerfstudio or torch.

The deletion of the Python implementation is a separate, later commit; until
then `tests/python/test_dataparser_parity.py` asserts the two agree. See
docs/restructure-proposal.md §4.1.
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
        return cfg


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
