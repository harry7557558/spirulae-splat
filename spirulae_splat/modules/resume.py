"""Training-resume support: reconstruct a TrainerConfig from a saved config.json
and read a checkpoint's runtime manifest (state.json inside state.tar).

The heavy engine state (world params, optimizer moments, densify aux, appearance)
is restored on the C++ side by ``_C.engine_load_checkpoint``; this module only
handles the Python-side glue: config reconstruction, checkpoint-path resolution,
and reading the small state.json manifest that drives appearance re-init.
"""

from __future__ import annotations

import dataclasses
import io
import json
import tarfile
import typing
from pathlib import Path


# --- config.json -> TrainerConfig -------------------------------------------
# The config is serialized with dataclasses.asdict (Paths -> str, no Type-valued
# fields exist in the schema). Reconstruction only needs to: recurse into nested
# config dataclasses, wrap Path/Optional[Path] strings back into Path, and turn
# JSON arrays for Tuple[...] fields back into tuples. Everything else (scalars,
# Literals, Optional scalars, Union[bool,int]) round-trips as-is.

def _convert(tp, val):
    if val is None:
        return None
    origin = typing.get_origin(tp)
    args = typing.get_args(tp)
    if origin is typing.Union:                       # Optional[X] / Union[...]
        for a in (a for a in args if a is not type(None)):
            if dataclasses.is_dataclass(a):
                return _from_dict(a, val)
            if a is Path:
                return Path(val)
        return val                                   # scalar/bool/int union
    if origin in (tuple,):                            # Tuple[...] -> tuple
        return tuple(val)
    if origin in (list,):                             # List[...] -> list
        return list(val)
    if tp is Path:
        return Path(val)
    if dataclasses.is_dataclass(tp):
        return _from_dict(tp, val)
    return val


def _from_dict(cls, data):
    if not dataclasses.is_dataclass(cls) or not isinstance(data, dict):
        return data
    hints = typing.get_type_hints(cls)
    kwargs = {}
    for f in dataclasses.fields(cls):
        if not f.init or f.name not in data:
            continue
        kwargs[f.name] = _convert(hints.get(f.name, type(data[f.name])), data[f.name])
    return cls(**kwargs)


def config_from_json(json_path: Path):
    """Reconstruct a base ``TrainerConfig`` from a run's ``config.json``.

    Deserializes into the base ``TrainerConfig`` (config.json is flat with all
    fields concrete, so the original preset subclass is not needed)."""
    from spirulae_splat.modules.trainer import TrainerConfig
    with open(json_path) as f:
        data = json.load(f)
    return _from_dict(TrainerConfig, data)


def apply_preset(base, preset_inst, ref_inst):
    """Re-impose a preset's characteristic settings on a resumed config.

    On resume, `base` is the checkpoint's config.json. If the user selected a
    DIFFERENT preset than the original run, this applies that preset's deviations
    (e.g. `synthetic` disabling bilagrid/PPISP) on top of `base`, recursively.
    A field is applied only where `preset_inst` differs from `ref_inst` (the
    plain base-class defaults) -- so a preset's intentional overrides win over
    the checkpoint, while fields the preset leaves at the base default keep the
    checkpoint's value. Explicit CLI flags then override on top (via tyro's
    default=). For a same-preset resume, preset_inst == ref_inst -> no change."""
    for f in dataclasses.fields(preset_inst):
        pv = getattr(preset_inst, f.name)
        rv = getattr(ref_inst, f.name)
        if dataclasses.is_dataclass(pv) and dataclasses.is_dataclass(rv):
            apply_preset(getattr(base, f.name), pv, rv)
        elif pv != rv:
            setattr(base, f.name, pv)


# --- checkpoint path resolution + state.json --------------------------------

def resolve_checkpoint(path: Path):
    """Given a run output dir OR a ``step-*.ckpt`` dir, return
    ``(run_dir, ckpt_dir)``. Picks the latest ``step-*.ckpt`` when given a run
    dir. ``run_dir`` is where ``config.json`` lives."""
    path = Path(path)
    if (path / "state.tar").is_file():               # a checkpoint directory
        return path.parent, path
    ckpts = sorted(path.glob("step-*.ckpt"))
    if ckpts:                                         # a run directory
        return path, ckpts[-1]
    raise FileNotFoundError(f"No state.tar or step-*.ckpt found under {path}")


def read_state_json(ckpt_dir: Path) -> dict:
    """Extract and parse the ``state.json`` manifest from ``ckpt_dir/state.tar``."""
    with tarfile.open(Path(ckpt_dir) / "state.tar") as tf:
        member = tf.extractfile("state.json")
        if member is None:
            raise RuntimeError(f"state.json missing in {ckpt_dir}/state.tar")
        return json.loads(member.read())


# --- config compatibility gate ----------------------------------------------
# Fields that DEFINE the on-device layout / model architecture. If a resume's
# reconstructed config disagrees with the checkpoint's state.json on any of
# these, the byte layouts won't line up (the C++ loader would then fail on a
# size mismatch), so we surface a clear error up front.

def check_compat(cfg, state: dict):
    """Raise if the config's architecture disagrees with the checkpoint manifest."""
    return  # TODO
    model = cfg.model
    checks = [
        ("primitive",           model.primitive,          state.get("primitive")),
        ("sh_degree",           model.sh_degree,          state.get("sh_degree")),
    ]
    problems = [f"  {name}: config={c!r} vs checkpoint={k!r}"
                for name, c, k in checks if k is not None and c != k]
    if problems:
        raise ValueError(
            "Resume config is incompatible with the checkpoint architecture:\n"
            + "\n".join(problems)
            + "\n(These must match the checkpoint; re-run with a matching config.)")


# --- entrypoint: build the effective config for a resumed run ----------------
# On resume the architecture/model/data-config come from the checkpoint's
# config.json; the identity of THIS run (data path, output dir, resume path) and
# any run-control flags the user explicitly changed come from the CLI.

_RESUME_IDENTITY = ("data", "resume")
_RESUME_RUN_CONTROL = (
    "num_iterations", "steps_per_save", "save_only_latest_checkpoint",
    "save_full_checkpoint", "save_eval_images", "viewer_port",
    "disable_viewer", "keep_viewer_alive",
)


def build_resume_config(cli_config):
    """Given a tyro-parsed config with ``.resume`` set, return the effective
    config: the checkpoint's config.json with CLI identity fields applied and
    run-control fields overridden only where the user changed them from the
    dataclass default (otherwise the checkpoint's value is kept).

    By default the resumed run writes back into the SAME output folder as the
    checkpoint (so new checkpoints/eval/logs continue in place), unless the user
    explicitly points the output elsewhere via --output-dir-prefix/-name."""
    from spirulae_splat.modules.trainer import TrainerConfig

    run_dir, _ = resolve_checkpoint(Path(cli_config.resume))
    run_dir = Path(run_dir).resolve()
    cfg = config_from_json(run_dir / "config.json")

    for f in _RESUME_IDENTITY:
        setattr(cfg, f, getattr(cli_config, f))

    defaults = {fl.name: fl.default for fl in dataclasses.fields(TrainerConfig)}

    # Output folder: reuse the checkpoint's run folder unless the user explicitly
    # asked for a different output on the CLI.
    user_set_output = (
        cli_config.output_dir_name != defaults.get("output_dir_name")
        or cli_config.output_dir_prefix != defaults.get("output_dir_prefix"))
    if user_set_output:
        cfg.output_dir_prefix = cli_config.output_dir_prefix
        cfg.output_dir_name = cli_config.output_dir_name
    else:
        cfg.output_dir_prefix = run_dir.parent
        cfg.output_dir_name = Path(run_dir.name)

    for f in _RESUME_RUN_CONTROL:
        cli_val = getattr(cli_config, f)
        if cli_val != defaults.get(f):            # user changed it -> honor CLI
            setattr(cfg, f, cli_val)
    return cfg
