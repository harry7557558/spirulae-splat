"""Regression gate on the native training driver.

This started as the §4.3 parity gate: it compared the `EngineStepConfig` the
Python path built (`model.py::engine_train_step_managed` +
`core.py::_build_{optim,densify}_config`) against the one
`TrainerCore::build_step_config` builds, for 8 config variants x 4 run states
x 20 steps chosen to straddle every warmup and decay boundary. They agreed on
every field of all 640 -- and the Python implementation has since been deleted,
so there is no second implementation left to compare with.

What survives is that proof, frozen: `step_config_golden.json` holds the
values the comparison validated, and `test_step_config_golden` asserts the C++
builder still produces them. A schedule that changes by accident fails here
instead of as a quality regression 20k steps into a run. Regenerate
deliberately with `make_step_config_golden.py`, and say what moved.

The rest:
  * `test_native_config_matches_dataclass_defaults` -- the generated
    `SsplatConfig` (defaults + preset branches) must equal what
    `native_trainer.to_native_config()` reads off the Python dataclasses.
    These are still two live representations, so this stays a true parity
    test; it fails the moment the generated header goes stale.
  * `test_trainer_drives_the_session` -- the rewired `Trainer` really runs on
    a `TrainerSession` and reports its splats/losses.
  * `test_native_session_trains` (opt-in) -- a real short run end to end.

Fixtures are generated from a fixed seed (dataset_fixtures.py), so the test
depends on no dataset present on the machine. Set `SSPLAT_TEST_DATASET` to
additionally run against a real one.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

import dataset_fixtures as fixtures

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("needs a CUDA device", allow_module_level=True)

from spirulae_splat.splat.cuda import _C  # noqa: E402
from spirulae_splat.modules.native_trainer import (  # noqa: E402
    PRESET_BY_CLASS_NAME, make_session, preset_name, to_native_config,
)
from spirulae_splat.modules import trainer as trainer_mod  # noqa: E402


# Steps chosen to straddle every schedule boundary in build_step_config:
# sh_degree warmup, distortion / reg / supervision / median / alpha warmups,
# the LR decay curve, and refine_start_iter.
PARITY_STEPS = [0, 1, 7, 99, 100, 499, 500, 501, 999, 1000, 1001,
                1999, 2000, 2999, 3000, 4499, 4500, 9999, 24999, 29999]

# Fixed inputs for the golden, so it does not depend on a dataset.
GOLDEN_NUM_ITERATIONS = 30000
GOLDEN_TRAIN_FRAME_SCALE = 3.0041950345141890

# The four (bilagrid rgb/depth/normal, ppisp) enablement combinations that
# select different branches of build_step_config.
RUN_STATE_VARIANTS = {
    "none":      dict(),
    "rgb":       dict(bilagrid_rgb_init=True),
    "rgb+ppisp": dict(bilagrid_rgb_init=True, ppisp_init=True),
    "all":       dict(bilagrid_rgb_init=True, bilagrid_depth_init=True,
                      bilagrid_normal_init=True, ppisp_init=True),
}

# Config variants: the branches that depend on the config rather than the run
# state. Keys are SsplatConfig (== CLI) field names.
CONFIG_VARIANTS = {
    "default": {},
    "adagrad": {"use_adagrad_bilagrid_optim": True, "use_adagrad_ppisp_optim": True},
    "background_sh": {"background_mode": "sh"},
    "background_noise": {"background_mode": "noise"},
    "no_quantization": {"quantization_level": 0},
    "world_grad_score": {"densify_score_blend_world_grad": 1.0,
                         "use_revised_densification": True},
    "legacy_densify": {"use_revised_densification": False},
    "max_steps": {"max_steps": 7000},
}

STEP_CONFIG_SECTIONS = {
    "loss": ["weights", "w_ssim", "num_loss_scales", "loss_scale_min_pixels",
             "compute_loss_map", "loss_map_mode", "robust_edge_aware_quantile",
             "overexposure_reg_weight", "color_shift_reg_weight",
             "color_shift_reg_beta", "input_depth_is_ray_depth"],
    "optim": ["lr_means", "lr_quats", "lr_scales", "lr_opacities",
              "lr_features_dc", "lr_features_sh", "max_gauss_ratio",
              "scale_regularization_weight", "mcmc_opacity_reg_weight",
              "mcmc_scale_reg_weight", "erank_reg_weight", "erank_reg_weight_s3",
              "quat_norm_reg_weight", "sh_reg_weight", "use_scale_agnostic_mean",
              "sh_optim_bits", "sh_value_bits", "non_sh_optim_bits",
              "quantization_level", "use_per_splat_bias_correction",
              "use_fused_proj_bwd_optim", "write_densify_world_grad_score",
              "split_batch", "use_color_trust_region", "color_is_linear",
              "eps_tr"],
    "densify": ["refine_start_iter", "refine_stop_num_iter", "refine_stop_iter",
                "refine_every", "growth_factor", "min_opacity",
                "max_screen_size", "max_screen_size_clip_hardness",
                "max_world_size", "noise_lr", "noise_lr_final",
                "use_revised_densification", "score_mode",
                "score_blend_world_grad", "las_split_opacity_k_init",
                "las_split_opacity_k_final", "las_split_opacity_k_warmup"],
    "bilagrid": ["lr_rgb", "lr_depth", "lr_normal", "tv_weight_rgb",
                 "tv_weight_depth", "tv_weight_normal"],
    "ppisp": ["lr", "reg_weights", "run_before_bilagrid"],
    "background": ["lr_dc", "lr_sh", "randomize_weight", "seed"],
}

GOLDEN_PATH = Path(__file__).with_name("step_config_golden.json")


def _eq(a, b) -> bool:
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    if isinstance(a, float) or isinstance(b, float):
        a, b = float(a), float(b)
        if math.isnan(a) and math.isnan(b):
            return True
        if math.isinf(a) or math.isinf(b):
            return a == b
        return math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-12)
    return a == b


# ---------------------------------------------------------------------------
# 1. Config conversion -- still two live representations
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("class_name,preset", sorted(PRESET_BY_CLASS_NAME.items()))
def test_native_config_matches_dataclass_defaults(class_name, preset):
    """to_native_config(PresetClass()) == SsplatConfig() + apply_preset(name).

    Both sides claim to encode the same defaults: the generator bakes them
    into cli_config.h at codegen time, the adapter reads them off the live
    dataclasses. If they disagree, either the generator is stale or a preset
    branch is wrong -- and every native run would silently use a different
    config than the CLI flags describe.
    """
    cls = getattr(trainer_mod, class_name)
    assert preset_name(cls(data=Path("/nonexistent"))) == preset

    from_python = to_native_config(cls(data=Path("/nonexistent")))
    from_generator = _C.SsplatConfig()
    _C.ssplat_apply_preset(from_generator, preset)
    from_generator.data = "/nonexistent"

    mismatches = []
    for cli_key, _pyname, _group, _choices, _help in _C.ssplat_config_fields():
        va, vb = getattr(from_python, cli_key), getattr(from_generator, cli_key)
        if not _eq(va, vb):
            mismatches.append(f"  {cli_key}: python={va!r} generated={vb!r}")
    assert not mismatches, \
        f"{class_name}/{preset}: SsplatConfig differs:\n" + "\n".join(mismatches)


# ---------------------------------------------------------------------------
# 2. Per-step EngineStepConfig, against the frozen parity result
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("variant", sorted(CONFIG_VARIANTS))
def test_step_config_golden(variant):
    golden = json.loads(GOLDEN_PATH.read_text())
    assert variant in golden, \
        f"{GOLDEN_PATH.name} has no entry for '{variant}' -- regenerate it"

    native_cfg = _C.SsplatConfig()
    native_cfg.data = "/nonexistent"
    native_cfg.num_iterations = GOLDEN_NUM_ITERATIONS
    for key, value in CONFIG_VARIANTS[variant].items():
        setattr(native_cfg, key, value)

    mismatches = []
    for run_state, flags in sorted(RUN_STATE_VARIANTS.items()):
        st = _C.RunState()
        st.train_frame_scale = GOLDEN_TRAIN_FRAME_SCALE
        st.splat_linear = False
        for name, value in flags.items():
            setattr(st, name, value)

        entry = golden[variant][run_state]
        for step in PARITY_STEPS:
            # Fields that do not vary with the step are stored once.
            want = {**entry["constant"], **entry["per_step"][str(step)]}
            got = _C.build_step_config(native_cfg, st, step)
            for key, expected in want.items():
                section, field = key.split(".")
                actual = getattr(getattr(got, section), field)
                if isinstance(actual, (list, tuple)):
                    actual = list(actual)
                if not _eq(actual, expected):
                    mismatches.append(
                        f"  {run_state}/step={step} {key}: "
                        f"got={actual!r} golden={expected!r}")
    assert not mismatches, (
        f"{variant}: build_step_config drifted from the golden "
        f"({len(mismatches)} field(s)):\n" + "\n".join(mismatches[:40]))


def test_golden_covers_every_step_config_field():
    """The golden must pin every field, or drift can hide in an unpinned one."""
    golden = json.loads(GOLDEN_PATH.read_text())
    expected = {f"{s}.{f}" for s, fs in STEP_CONFIG_SECTIONS.items() for f in fs}
    for variant, per_variant in golden.items():
        for run_state, entry in per_variant.items():
            got = set(entry["constant"]) | set(
                entry["per_step"][str(PARITY_STEPS[0])])
            assert got == expected, (
                f"{variant}/{run_state}: golden field set differs from "
                f"STEP_CONFIG_SECTIONS; missing={sorted(expected - got)} "
                f"extra={sorted(got - expected)}")


# ---------------------------------------------------------------------------
# 3. The rewired Trainer really runs on the session
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    root = tmp_path_factory.mktemp("trainer_fixture")
    return fixtures.write_colmap_text(root / "colmap", fixtures.make_scene())


def test_trainer_drives_the_session(dataset, tmp_path):
    cfg = trainer_mod.TrainerConfig(data=Path(dataset))
    cfg.output_dir_prefix = Path(tmp_path)
    cfg.num_iterations = 4
    cfg.steps_per_save = 0

    t = trainer_mod.Trainer(cfg)
    try:
        # The engine world came from the session's seeding, and the model's
        # host parameters were synced from it rather than uploaded.
        assert t.model.core.cur_num_splats == int(_C.engine_get_cur_num_splats())
        assert float(t.model.gauss_params["means"].abs().sum()) > 0, \
            "attach mode did not pull the device world into host params"
        # The run state the session decided is what the model reports.
        assert t.model._bilagrid_rgb_init == t.session.run_state.bilagrid_rgb_init
        assert t.model._ppisp_init == t.session.run_state.ppisp_init

        t.train()
        assert t.session.cur_step == cfg.num_iterations
        assert (t.output_dir / "config.json").is_file()
        assert (t.output_dir / "dataparser_transforms.json").is_file()
    finally:
        _C.engine_reset()


def test_trainer_serves_the_native_viewer(dataset, tmp_path):
    import socket
    import urllib.request

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]

    cfg = trainer_mod.TrainerConfig(data=Path(dataset))
    cfg.output_dir_prefix = Path(tmp_path)
    cfg.num_iterations = 2
    cfg.steps_per_save = 0

    t = trainer_mod.Trainer(cfg)
    try:
        t.start_viewer(host="127.0.0.1", port=port)
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/progress",
                                    timeout=20) as r:
            progress = json.loads(r.read())
        assert progress["total_steps"] == cfg.num_iterations
        # toggle_pause drives the session's atomic, which the render worker
        # and the training loop both read.
        assert t.toggle_pause() is True and t.session.paused is True
        assert t.toggle_pause() is False
    finally:
        t.stop_viewer()
        _C.engine_reset()


# ---------------------------------------------------------------------------
# 4. End-to-end native run (opt-in: needs a real dataset)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.environ.get("SSPLAT_TEST_DATASET"),
                    reason="set SSPLAT_TEST_DATASET to a dataset directory")
def test_native_session_trains(tmp_path):
    """A short real run through TrainerSession, with densification on.

    refine_stop_num_iter counts back from the end, so the refine window is
    pulled in explicitly -- a short run with the defaults never densifies and
    would leave the densify half of the step config untested.
    """
    cfg = trainer_mod.TrainerConfig(data=Path(os.environ["SSPLAT_TEST_DATASET"]))
    cfg.output_dir_prefix = Path(tmp_path)
    cfg.num_iterations = 300
    cfg.steps_per_save = 0
    cfg.model.refine_start_iter = 50
    cfg.model.refine_stop_num_iter = 50
    cfg.model.refine_every = 50
    cfg.model.cap_max = 200000
    cfg.dataparser.image_dir = os.environ.get("SSPLAT_TEST_IMAGE_DIR", "images")
    cfg.dataparser.rescale_camera_to_fit = int(
        os.environ.get("SSPLAT_TEST_DOWNSCALE", "0"))
    cfg.dataparser.downscale_rounding_mode = "round"

    sess = make_session(cfg)
    sess.load_dataset()
    sess.setup_engine()
    assert sess.dataset.num_cameras > 0

    seen = []
    sess.train(on_step=lambda p: seen.append(
        (p.step, dict(p.losses), p.num_splats)))
    _C.engine_reset()

    assert len(seen) == cfg.num_iterations
    assert seen[0][0] == 0 and seen[-1][0] == cfg.num_iterations - 1

    def window(key, sl):
        vals = [v[key] for _, v, _ in seen[sl]]
        return sum(vals) / len(vals)

    # Step 0's metrics are not back from the device yet (the readback is one
    # step behind), so the "early" window starts at 1.
    for key, improves in (("rgb_loss", -1), ("psnr", +1)):
        first, last = window(key, slice(1, 21)), window(key, slice(-20, None))
        assert math.isfinite(first) and math.isfinite(last)
        assert improves * (last - first) > 0, \
            f"{key} did not improve: {first} -> {last}"

    # Densification ran: the refine window above is [50, 250].
    assert seen[-1][2] > seen[0][2], "no densification happened"
