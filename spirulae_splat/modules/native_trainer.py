"""Adapter over the native training session (src/app/TrainerCore.h).

`TrainerCore` is the training driver the standalone CLI (`ssplat-train`) and
the native GUI already run, and its header comment describes itself as a
line-by-line port of the Python managed path -- a comment doing the job a
shared implementation should do. This module is the Python client for it, and
is the intended replacement for `model.py::engine_train_step_managed` plus
`core.py::_build_optim_config` / `_build_densify_config`.

Nothing here imports nerfstudio. torch is imported only to load `_C`.

The deletion of the Python implementation is a separate, later commit; until
then `tests/python/test_trainer_parity.py` asserts the two agree, step by
step. See docs/restructure-proposal.md §4.3.

Config conversion
-----------------
`SsplatConfig` is generated from the Python dataclasses by
`tools/codegen/generate_cli_config.py`, and that generator also emits the
field table this module walks (`_C.ssplat_config_fields()` -> (cli_key,
pyname, group, choices, help)). So the flattening -- which Python sub-config a
field lives in, and what it is called after the rename pass -- exists exactly
once, in the generator. Adding a field to a dataclass and re-running the
generator is enough; nothing here needs editing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, Optional

# The preset subcommand names ss_trainer.py exposes, keyed by the config class
# that implements each. Used to record the preset in the run's config.json,
# exactly as the CLI does; the overrides themselves are already baked into the
# dataclass, so nothing is re-applied.
PRESET_BY_CLASS_NAME = {
    "TrainerConfig": "3dgs",
    "TrainerConfig360Camera": "360-camera",
    "TrainerConfigInTheWild": "in-the-wild",
    "TrainerConfigLinear": "linear-color",
    "TrainerConfigSynthetic": "synthetic",
    "TrainerConfigMeshing": "meshing",
    "TrainerConfigAcademicBaseline": "academic-baseline",
}

# Fields whose Python type the generator overrode (TYPE_OVERRIDES in
# generate_cli_config.py), so a straight copy is wrong.
#   rescale_camera_to_fit: Union[bool, int] on the Python side.
#     True  -> probe the image resolution        -> -1 natively
#     False -> off                               ->  0
#     int   -> divide intrinsics by this         ->  the number
def _rescale_camera_to_fit(v) -> float:
    if v is None or v is False:
        return 0.0
    if v is True:
        return -1.0
    return float(v)


_VALUE_OVERRIDES: Dict[str, Callable[[Any], Any]] = {
    "rescale_camera_to_fit": _rescale_camera_to_fit,
}


def _native():
    from spirulae_splat.splat.cuda import _C
    return _C


def preset_name(config) -> str:
    """The preset subcommand a TrainerConfig corresponds to ("3dgs" if none)."""
    return PRESET_BY_CLASS_NAME.get(type(config).__name__, "3dgs")


def to_native_config(config):
    """TrainerConfig (with its nested dataclasses) -> `_C.SsplatConfig`."""
    native = _native().SsplatConfig()

    for cli_key, pyname, group, choices, _help in _native().ssplat_config_fields():
        holder = config if group == "trainer" else getattr(config, group)
        try:
            value = getattr(holder, pyname)
        except AttributeError as e:
            raise AttributeError(
                f"{type(config).__name__}.{group}.{pyname} is missing; "
                f"cli_config.h is stale -- re-run "
                f"tools/codegen/generate_cli_config.py") from e

        override = _VALUE_OVERRIDES.get(cli_key)
        if override is not None:
            value = override(value)
        elif isinstance(value, Path):
            value = str(value)
        elif value is None:
            # A string field spells None as "" (the generator marks those with
            # a "none" choice); anything else is a genuine std::optional.
            if "none" in choices.split("|"):
                value = ""
        elif isinstance(value, tuple):
            value = list(value)

        try:
            setattr(native, cli_key, value)
        except TypeError as e:
            raise TypeError(
                f"cannot set SsplatConfig.{cli_key} from "
                f"{group}.{pyname}={value!r} ({type(value).__name__})") from e

    return native


def make_session(
    config,
    front_end_handles_resume: bool = True,
    front_end_handles_eval: bool = True,
    log_fn: Optional[Callable[[str], None]] = None,
):
    """A `_C.TrainerSession` configured from a Python `TrainerConfig`.

    Returns it after `check_config()`; the caller drives `load_dataset()` and
    `setup_engine()`, since a front-end usually wants to do its own work
    between the phases (the Python trainer restores a checkpoint after
    `setup_engine()`, which is the point at which the engine skeleton exists).

    The two `front_end_handles_*` flags default to True because the Python
    trainer does implement resume and eval; the standalone CLI does not, and
    its guards must keep firing.
    """
    sess = _native().TrainerSession()
    sess.config = to_native_config(config)
    sess.preset = preset_name(config)
    sess.front_end_handles_resume = front_end_handles_resume
    sess.front_end_handles_eval = front_end_handles_eval
    if log_fn is not None:
        sess.set_log_fn(log_fn)
    sess.check_config()
    return sess
