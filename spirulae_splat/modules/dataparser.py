"""The dataset-parsing config.

The parsing *implementation* that used to live here (plus `colmap_utils.py`
and `metashape_utils.py`) is gone: COLMAP / Nerfstudio / Metashape are parsed
by `src/data/parsers/` and reached from Python through
`modules/native_dataparser.py`. See docs/restructure-proposal.md §4.1.

The config dataclass stays, and stays HERE, because it is the source of truth
for two generated artifacts and one CLI:
  * `tools/codegen/generate_cli_config.py` AST-parses this file to emit
    `src/app/generated/cli_config.h` (and so the native CLI's flags);
  * `ss_trainer.py` builds its tyro CLI from it, and `--resume` reads it back
    out of a run's config.json;
  * `native_dataparser.to_native_parser_config()` maps it onto the native
    `DatasetParserConfig`.
Editing a field here therefore changes the native trainer too -- re-run the
generator.

Fields with no native counterpart (scene_scale, orientation_method,
center_method, auto_scale_poses, train_frame != "points") are kept so old
config.json files still load; `TrainerSession.check_config()` rejects or warns
about each.

The two preprocessing helpers moved to `scripts/`, which still parses in
Python: `scripts/colmap_utils.py`, `scripts/metashape_utils.py`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, Union


# Distortion columns, in the order both parsers emit them
# (native: ParsedDataset.dist_coeffs). Consumed by ss_meshing.py, which remaps
# to the engine's own ordering.
DISTORTION_KEYS = "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split()


@dataclass
class SpirulaeSplatDataParserConfig:
    """Spirulae-Splat dataset config"""

    data_format: Literal["colmap", "nerfstudio", "metashape", None] = None
    """Data format, leave None to auto detect"""
    colmap_recon_dir: Optional[Path] = None
    """Path to COLMAP reconstruction relative to dataset directory (e.g. sparse/0). Will auto detect if not specified (picking the model with the most registered images when several exist)."""
    image_dir: Path = Path("images")
    """Path to images relative to dataset directory, used for COLMAP and Metashape datasets"""
    mask_dir: Path = Path("masks")
    """Path to masks relative to dataset directory, used for COLMAP and Metashape datasets"""
    depth_dir: Path = Path("depths")
    """Path to depth maps relative to dataset directory, used for COLMAP and Metashape datasets"""
    normal_dir: Path = Path("normals")
    """Path to normal maps relative to dataset directory, used for COLMAP and Metashape datasets"""
    metashape_xml: Optional[Path] = None
    """Path to the Metashape xml file. Will automatically detect if not specified."""
    metashape_ply: Optional[Path] = None
    """Path to the Metashape point export ply file. Will automatically detect if not specified."""
    metashape_psx: Optional[Path] = None
    """Path to Metashape PSX file, used to resolve file name ambiguity when there are multiple images with the same file name"""
    rescale_camera_to_fit: Union[bool, int] = False
    """Whether to check if image resolution match camera resolution and scale camera intrinsics accordingly if not.
        Set this to a number to divide intrinsics by that number, e.g. Mip-NeRF 360 and Zip-NeRF with images_(2|4)
        Set this to True to detect resolution, e.g. tankt_db"""
    downscale_rounding_mode: Literal["floor", "ceil", "round"] = "floor"
    """Rounding mode applied to camera width/height when dividing by `rescale_camera_to_fit`.
        Use `round` to match the convention used by most image downscalers (e.g. Mip-NeRF 360 images_(2|4|8))."""

    scene_scale: float = 1.0
    """How much to scale the region of interest by."""
    orientation_method: Literal["pca", "up", "vertical", "none", "gsplat"] = "up"
    """The method to use for orientation."""
    center_method: Literal["poses", "focus", "none", "gsplat"] = "poses"
    """The method to use to center the poses."""
    auto_scale_poses: bool = True
    """Whether to automatically scale the poses to fit in +/- 1 bounding box."""
    outlier_threshold: float = float('inf')
    """Threshold to reject outlier camera poses."""
    train_frame: Literal["normalized", "camera", "points"] = "points"
    """Coordinate frame in which splats are trained."""
    eval_mode: Literal["fraction", "filename", "interval", "all"] = "all"
    """
    The method to use for splitting the dataset into train and eval.
    Fraction splits based on a percentage for train and the remaining for eval.
    Filename splits based on filenames containing train/eval.
    Interval uses every nth frame for eval.
    All uses all the images for any split.
    """
    train_split_fraction: float = 0.9
    """The percentage of the dataset to use for training. Only used when eval_mode is train-split-fraction."""
    eval_interval: int = 8
    """The interval between frames to use for eval. Only used when eval_mode is eval-interval."""
    depth_unit_scale_factor: float = 1e-3
    """Scales the depth values to meters. Default value is 0.001 for a millimeter to meter conversion."""

    validation_fraction: float = 0.0
    """Use this fraction of training images for validation. Stop training when performance on validation images start to drop."""


