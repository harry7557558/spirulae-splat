"""Regression gate on the native dataset parsers.

This started as the §4.1 parity gate: it ran COLMAP (text + binary),
Nerfstudio and Metashape fixtures through **both** the native parsers and
`modules/dataparser.py`, and asserted the frame set, poses, intrinsics,
distortion, seed cloud, train-frame scalars and both split sides all agreed.
They did, for 4 formats x 4 config variants x 2 splits -- and the Python
implementation has since been deleted (`dataparser.py` is now just the config
dataclass; `colmap_utils.py` / `metashape_utils.py` moved to `scripts/`, which
still parses in Python for preprocessing).

What survives is that proof, frozen: `dataparser_golden.json` holds the values
the comparison validated, and this file asserts the native parsers still
produce them. Regenerate deliberately with `make_dataparser_golden.py`.

Fixtures are generated from a fixed seed (`dataset_fixtures.py`), so the test
depends on no dataset present on the machine.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import dataset_fixtures as fixtures
from spirulae_splat.modules.native_dataparser import (
    NativeParserConfig, parse_dataset, to_dataparser_outputs,
)

import make_dataparser_golden as mk

GOLDEN = json.loads(
    (Path(__file__).with_name("dataparser_golden.json")).read_text())

# Per-field tolerances, carried over from the parity comparison: the two
# implementations reached the same pose through different orderings of the
# same float32 algebra, and the golden stores 7 significant digits.
TOLERANCES = {
    "c2w": 2e-5,
    "intrins": 1e-4,
    "dist_coeffs": 1e-7,
    "train_to_normalized": 2e-5,
    "train_frame_scale": 1e-5,
}
EXACT_FIELDS = ("num_cameras", "image_names", "widths", "heights",
                "camera_models", "val_indices", "num_points",
                "points_xyz_sha256", "points_rgb_sha256")


@pytest.fixture(scope="module")
def datasets(tmp_path_factory):
    root = tmp_path_factory.mktemp("ssplat_fixtures")
    scene = fixtures.make_scene()
    return {name: writer(root / name, scene)
            for name, writer in mk.FORMAT_WRITERS.items()}


@pytest.mark.parametrize("fmt", sorted(mk.FORMAT_WRITERS))
@pytest.mark.parametrize("variant", sorted(mk.CONFIG_VARIANTS))
@pytest.mark.parametrize("split", ["train", "eval"])
def test_parser_golden(datasets, fmt, variant, split):
    overrides = mk.CONFIG_VARIANTS[variant]
    ds = parse_dataset(datasets[fmt],
                       NativeParserConfig(split=split, **overrides))
    got = mk.describe(ds)
    want = GOLDEN[fmt][variant][split]
    tag = f"{fmt}/{variant}/{split}"

    for key in EXACT_FIELDS:
        assert got[key] == want[key], f"{tag}: {key}"
    for key, atol in TOLERANCES.items():
        np.testing.assert_allclose(
            np.asarray(got[key], dtype=np.float64),
            np.asarray(want[key], dtype=np.float64),
            atol=atol, rtol=1e-5, err_msg=f"{tag}: {key}")


def test_golden_covers_every_case():
    """Every fixture/variant/split combination must be pinned."""
    for fmt in mk.FORMAT_WRITERS:
        assert fmt in GOLDEN, fmt
        for variant in mk.CONFIG_VARIANTS:
            assert variant in GOLDEN[fmt], (fmt, variant)
            assert set(GOLDEN[fmt][variant]) == {"train", "eval"}, (fmt, variant)


def test_splits_partition_the_frames(datasets):
    """train and eval must be complementary -- and cover everything.

    A parser bug that dropped frames from both sides would keep each side
    self-consistent (so the golden would still pass) while silently shrinking
    the dataset.
    """
    for fmt in sorted(mk.FORMAT_WRITERS):
        for variant, overrides in sorted(mk.CONFIG_VARIANTS.items()):
            tr = parse_dataset(datasets[fmt], NativeParserConfig(**overrides))
            ev = parse_dataset(datasets[fmt],
                               NativeParserConfig(split="eval", **overrides))
            tr_names = {Path(f).name for f in tr.image_filenames}
            ev_names = {Path(f).name for f in ev.image_filenames}
            tag = f"{fmt}/{variant}"
            if overrides.get("eval_mode", "all") == "all":
                # "All uses all the images for any split."
                assert tr_names == ev_names, tag
            else:
                assert not (tr_names & ev_names), f"{tag}: splits overlap"
                assert len(tr_names) + len(ev_names) == fixtures.N_IMAGES, \
                    f"{tag}: splits do not cover every frame"
            # Anything derived from the full camera set is split-independent.
            assert tr.train_frame_scale == pytest.approx(
                ev.train_frame_scale, rel=1e-6), tag


def test_formats_describe_the_same_scene(datasets):
    """Cross-check the fixtures themselves.

    All four are written from one synthetic scene, so the native parser must
    recover the same cameras from each. A fixture bug would otherwise make the
    golden agree on garbage.
    """
    ref = parse_dataset(datasets["colmap_text"], NativeParserConfig())
    for fmt in sorted(mk.FORMAT_WRITERS):
        got = parse_dataset(datasets[fmt], NativeParserConfig())
        assert got.num_cameras == ref.num_cameras, fmt
        np.testing.assert_allclose(got.c2w, ref.c2w, atol=1e-4,
                                   err_msg=f"{fmt} c2w vs colmap_text")
        np.testing.assert_allclose(got.intrins, ref.intrins, atol=1e-3,
                                   err_msg=f"{fmt} intrins vs colmap_text")
        np.testing.assert_allclose(got.train_frame_scale, ref.train_frame_scale,
                                   rtol=1e-4, err_msg=f"{fmt} scale")


def test_anisotropic_focal(tmp_path):
    """fx != fy, on the formats that can express it.

    (Metashape stores a single <f>, so it cannot.)
    """
    scene = fixtures.make_scene(fy=47.5)
    for fmt in ("colmap_text", "nerfstudio"):
        d = mk.FORMAT_WRITERS[fmt](tmp_path / fmt, scene)
        ds = parse_dataset(d, NativeParserConfig())
        fx_, fy_ = ds.intrins[:, 0], ds.intrins[:, 1]
        assert not np.allclose(fx_, fy_), f"{fmt}: fy override was lost"
        np.testing.assert_allclose(fy_, 47.5, rtol=1e-4)


def test_dataparser_outputs_adapter(datasets):
    """The dict the eval loader and the model still consume."""
    ds = parse_dataset(datasets["colmap_text"], NativeParserConfig())
    dpo = to_dataparser_outputs(ds, depth_unit_scale_factor=1e-3)

    cams = dpo["cameras"]
    assert len(cams) == ds.num_cameras
    np.testing.assert_allclose(
        cams.camera_to_worlds.numpy().reshape(-1), np.asarray(ds.c2w).reshape(-1),
        atol=1e-6)
    assert len(dpo["image_filenames"]) == ds.num_cameras
    assert dpo["metadata"]["points3D_xyz"].shape[1] == 3
    assert dpo["metadata"]["depth_unit_scale_factor"] == 1e-3
    # dataparser_transform / _scale must recompose into train_to_normalized.
    T = np.eye(4)
    T[:3, :] = dpo["dataparser_transform"].numpy() * dpo["dataparser_scale"]
    np.testing.assert_allclose(
        T.reshape(-1), np.asarray(ds.train_to_normalized).reshape(-1), atol=1e-5)
