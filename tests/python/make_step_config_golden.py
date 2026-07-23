#!/usr/bin/env python3
"""Regenerate tests/python/step_config_golden.json.

The golden values were captured while `model.py::engine_train_step_managed`
still existed and `test_trainer_parity.py` proved the two implementations
emitted an identical EngineStepConfig for every entry in this sweep. The
Python implementation is gone; the golden is what is left of that proof, and
it turns the parity gate into a regression gate on the surviving C++ one.

Run this ONLY when a step-config change is intentional, and say in the commit
message which fields moved and why -- a regenerated golden with no explanation
is indistinguishable from a silently broken schedule.

    python tests/python/make_step_config_golden.py
"""

from __future__ import annotations

import json
from pathlib import Path

from spirulae_splat.splat.cuda import _C

import test_trainer_parity as tp


def dump_config(cfg) -> dict:
    out = {}
    for section, fields in tp.STEP_CONFIG_SECTIONS.items():
        sec = getattr(cfg, section)
        for f in fields:
            v = getattr(sec, f)
            if isinstance(v, (list, tuple)):
                v = [_num(x) for x in v]
            else:
                v = _num(v)
            out[f"{section}.{f}"] = v
    return out


def _num(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, float):
        # float32 values; 9 significant digits round-trips them exactly.
        return float(f"{v:.9g}")
    return v


def split_constant(per_step: dict) -> dict:
    """{step: {field: value}} -> {"constant": ..., "per_step": ...}.

    Most of the ~60 fields are config passthroughs that do not depend on the
    step; hoisting them out cuts the file ~8x and makes the remainder --
    exactly the scheduled quantities -- the thing you read in a diff.
    """
    steps = list(per_step)
    fields = list(per_step[steps[0]])
    constant, varying = {}, []
    for f in fields:
        v0 = per_step[steps[0]][f]
        if all(per_step[s][f] == v0 for s in steps[1:]):
            constant[f] = v0
        else:
            varying.append(f)
    return {
        "constant": constant,
        "per_step": {s: {f: per_step[s][f] for f in varying} for s in steps},
    }


def main():
    golden = {}
    for variant, overrides in sorted(tp.CONFIG_VARIANTS.items()):
        native_cfg = _C.SsplatConfig()
        native_cfg.data = "/nonexistent"
        native_cfg.num_iterations = tp.GOLDEN_NUM_ITERATIONS
        for key, value in overrides.items():
            setattr(native_cfg, key, value)

        per_variant = {}
        for run_state, flags in sorted(tp.RUN_STATE_VARIANTS.items()):
            st = _C.RunState()
            st.train_frame_scale = tp.GOLDEN_TRAIN_FRAME_SCALE
            st.splat_linear = False
            for name, value in flags.items():
                setattr(st, name, value)
            per_variant[run_state] = split_constant({
                str(step): dump_config(_C.build_step_config(native_cfg, st, step))
                for step in tp.PARITY_STEPS
            })
        golden[variant] = per_variant

    path = Path(__file__).with_name("step_config_golden.json")
    path.write_text(json.dumps(golden, indent=1, sort_keys=True) + "\n")
    n = len(tp.CONFIG_VARIANTS) * len(tp.RUN_STATE_VARIANTS) * len(tp.PARITY_STEPS)
    print(f"wrote {path} ({n} step configs, {path.stat().st_size/1024:.0f} KiB)")


if __name__ == "__main__":
    main()
