"""Golden parity gate: the native dataset parsers vs the Python ones.

This is the verification gate for docs/restructure-proposal.md §4.1. It must
stay green until `modules/dataparser.py` (+ `colmap_utils.py`,
`metashape_utils.py`) are deleted; at that point it can be reduced to a
regression test of the native path alone.

Fixtures are generated from a fixed seed (dataset_fixtures.py) so the test
depends on no dataset that happens to be present on the machine. Set
`SSPLAT_TEST_DATASET=/path/to/a/real/dataset` to additionally run the same
comparison against a real one.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import dataset_fixtures as fixtures
from spirulae_splat.modules.native_dataparser import (
    NativeParserConfig, parse_dataset as parse_native,
)

nerfstudio = pytest.importorskip(
    "nerfstudio", reason="the Python dataparser half of the comparison needs it")

from spirulae_splat.modules.dataparser import (  # noqa: E402
    SpirulaeSplatDataParserConfig, SpirulaeSplatDataparser,
)

# CameraModelType in src/core/CameraModel.h
CAMERA_MODEL_IDS = {
    "PINHOLE": 0, "FISHEYE": 1, "EQUISOLID": 2, "EQUIRECTANGULAR": 3,
}

FORMAT_WRITERS = {
    "colmap_text": fixtures.write_colmap_text,
    "colmap_binary": fixtures.write_colmap_binary,
    "nerfstudio": fixtures.write_nerfstudio,
    "metashape": fixtures.write_metashape,
}

# Config variants exercised on every format: the split logic and the outlier
# filter are the parts most likely to drift between the two implementations.
CONFIG_VARIANTS = {
    "default": {},
    "eval_interval": {"eval_mode": "interval", "eval_interval": 3},
    "eval_fraction": {"eval_mode": "fraction", "train_split_fraction": 0.7},
    "validation": {"validation_fraction": 0.25},
}


def _np(x):
    """torch tensor / list / ndarray -> ndarray."""
    if hasattr(x, "detach"):
        x = x.detach().cpu()
    return np.asarray(x)


def _run_python(dataset_dir: Path, overrides: dict):
    cfg = SpirulaeSplatDataParserConfig(**overrides)
    # The Python parser returns (train_split, eval_split); the native parser
    # returns the train subset only, which is what we compare.
    return SpirulaeSplatDataparser(cfg, Path(dataset_dir)).parse()[0]


def _run_native(dataset_dir: Path, overrides: dict):
    return parse_native(dataset_dir, NativeParserConfig(**overrides))


@pytest.fixture(scope="module")
def scene():
    return fixtures.make_scene()


@pytest.fixture(scope="module")
def datasets(scene, tmp_path_factory):
    root = tmp_path_factory.mktemp("ssplat_fixtures")
    return {name: writer(root / name, scene)
            for name, writer in FORMAT_WRITERS.items()}


@pytest.mark.parametrize("fmt", sorted(FORMAT_WRITERS))
@pytest.mark.parametrize("variant", sorted(CONFIG_VARIANTS))
def test_parser_parity(datasets, fmt, variant):
    dataset_dir = datasets[fmt]
    overrides = CONFIG_VARIANTS[variant]

    py = _run_python(dataset_dir, overrides)
    cc = _run_native(dataset_dir, overrides)

    _assert_same(py, cc, f"{fmt}/{variant}")


def _assert_same(py, cc, tag: str):
    cams = py["cameras"]

    # ---- frame set and order -------------------------------------------
    py_names = [Path(str(f)).name for f in py["image_filenames"]]
    cc_names = [Path(f).name for f in cc.image_filenames]
    assert py_names == cc_names, f"{tag}: frame set/order differs"
    assert cc.num_cameras == len(py_names)

    # ---- resolution and camera model -----------------------------------
    np.testing.assert_array_equal(
        _np(cams.width).reshape(-1), cc.widths, err_msg=f"{tag}: width")
    np.testing.assert_array_equal(
        _np(cams.height).reshape(-1), cc.heights, err_msg=f"{tag}: height")
    py_models = np.array([CAMERA_MODEL_IDS[str(t)] for t in cams.camera_type])
    np.testing.assert_array_equal(
        py_models, cc.camera_models, err_msg=f"{tag}: camera model")

    # ---- geometry -------------------------------------------------------
    # float32 on both sides, but the two implementations reach the pose
    # through different orderings of the same algebra.
    np.testing.assert_allclose(
        _np(cams.camera_to_worlds), cc.c2w, atol=2e-5, rtol=1e-5,
        err_msg=f"{tag}: camera_to_worlds")
    np.testing.assert_allclose(
        _np(cams.intrins), cc.intrins, atol=1e-4, rtol=1e-6,
        err_msg=f"{tag}: intrinsics")
    np.testing.assert_allclose(
        _np(cams.distortion_params), cc.dist_coeffs, atol=1e-7, rtol=1e-5,
        err_msg=f"{tag}: distortion")

    # ---- seed point cloud ----------------------------------------------
    np.testing.assert_allclose(
        _np(py["metadata"]["points3D_xyz"]), cc.points_xyz,
        atol=1e-5, rtol=1e-5, err_msg=f"{tag}: point xyz")
    np.testing.assert_array_equal(
        _np(py["metadata"]["points3D_rgb"]), cc.points_rgb,
        err_msg=f"{tag}: point rgb")

    # ---- training-frame scalars ----------------------------------------
    np.testing.assert_allclose(
        py["train_frame_scale"], cc.train_frame_scale, rtol=1e-5,
        err_msg=f"{tag}: train_frame_scale")
    np.testing.assert_allclose(
        _np(py["train_to_normalized_transform"]).reshape(4, 4),
        cc.train_to_normalized, atol=2e-5, rtol=1e-5,
        err_msg=f"{tag}: train_to_normalized")

    # ---- validation split ----------------------------------------------
    py_val = np.asarray(py["metadata"].get("val_indices", []), dtype=np.int64)
    np.testing.assert_array_equal(
        np.sort(py_val), np.sort(np.asarray(cc.val_indices, dtype=np.int64)),
        err_msg=f"{tag}: val_indices")


def test_formats_describe_the_same_scene(datasets):
    """Cross-check the fixtures themselves.

    All four are written from one synthetic scene, so the native parser must
    recover the same cameras from each. A fixture bug would otherwise make the
    parity comparisons agree on garbage.
    """
    ref = _run_native(datasets["colmap_text"], {})
    for fmt in sorted(FORMAT_WRITERS):
        got = _run_native(datasets[fmt], {})
        assert got.num_cameras == ref.num_cameras, fmt
        np.testing.assert_allclose(got.c2w, ref.c2w, atol=1e-4,
                                   err_msg=f"{fmt} c2w vs colmap_text")
        np.testing.assert_allclose(got.intrins, ref.intrins, atol=1e-3,
                                   err_msg=f"{fmt} intrins vs colmap_text")
        np.testing.assert_allclose(got.train_frame_scale, ref.train_frame_scale,
                                   rtol=1e-4, err_msg=f"{fmt} scale")


def test_anisotropic_focal_colmap(tmp_path):
    """fx != fy, on the formats that can express it."""
    scene = fixtures.make_scene(fy=47.5)
    for fmt in ("colmap_text", "nerfstudio"):
        d = FORMAT_WRITERS[fmt](tmp_path / fmt, scene)
        _assert_same(_run_python(d, {}), _run_native(d, {}), f"{fmt}/aniso")


@pytest.mark.skipif(not os.environ.get("SSPLAT_TEST_DATASET"),
                    reason="set SSPLAT_TEST_DATASET to a dataset directory")
@pytest.mark.parametrize("variant", sorted(CONFIG_VARIANTS))
def test_parser_parity_real_dataset(variant):
    """Same comparison against a real dataset, opt-in via env var."""
    dataset_dir = Path(os.environ["SSPLAT_TEST_DATASET"])
    overrides = CONFIG_VARIANTS[variant]
    py = _run_python(dataset_dir, overrides)
    cc = _run_native(dataset_dir, overrides)
    _assert_same(py, cc, f"SSPLAT_TEST_DATASET/{variant}")
