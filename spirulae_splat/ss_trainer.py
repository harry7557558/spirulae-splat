import threading

# import os
# os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
# os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
# os.environ["PYTORCH_NO_CUDA_MEMORY_CACHING"] = "1"

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated


def _find_resume(argv):
    """Extract the --resume value from raw argv before tyro parsing."""
    for i, a in enumerate(argv):
        if a == "--resume" and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith("--resume="):
            return a.split("=", 1)[1]
    return None


def entrypoint():
    import sys
    import tyro
    from pathlib import Path

    presets = {
        "3dgs": TrainerConfig,
        "360-camera": TrainerConfig360Camera,
        "in-the-wild": TrainerConfigInTheWild,
        "linear-color": TrainerConfigLinear,
        "synthetic": TrainerConfigSynthetic,
        "meshing": TrainerConfigMeshing,
        "academic-baseline": TrainerConfigAcademicBaseline,
    }

    argv = sys.argv[1:]
    resume = _find_resume(argv)

    if resume is not None:
        # Resume: the checkpoint's config.json is the BASE. If the user selected a
        # DIFFERENT preset than the original run, that preset's characteristic
        # settings (e.g. `synthetic` disabling bilagrid/PPISP) are re-imposed on
        # top; then flags the user explicitly passes (incl. --model.* overrides)
        # win, via tyro's `default=`. Output defaults back to the checkpoint's run
        # folder unless --output-dir-* is given. The Trainer restores + layout-
        # adapts the engine state in __init__.
        from spirulae_splat.modules.resume import (
            config_from_json, resolve_checkpoint, apply_preset, check_resumable)
        run_dir, ckpt_dir = resolve_checkpoint(Path(resume))
        check_resumable(ckpt_dir)          # fail fast before building the model
        base = config_from_json(run_dir / "config.json")
        base.resume = Path(resume)
        base.output_dir_prefix = run_dir.parent
        base.output_dir_name = Path(run_dir.name)
        # Leading token is the preset subcommand (e.g. "synthetic"); re-impose its
        # deviations from the base defaults, then parse the remaining flags.
        sub = argv[0] if (argv and not argv[0].startswith("-")) else None
        rest = argv[1:] if sub is not None else argv
        if sub in presets:
            apply_preset(base, presets[sub](data=base.data),
                         TrainerConfig(data=base.data))
        config = tyro.cli(TrainerConfig, default=base, args=rest)
    else:
        Config = Union[
            Annotated[TrainerConfig, tyro.conf.subcommand(name="3dgs")],
            Annotated[TrainerConfig360Camera, tyro.conf.subcommand(name="360-camera")],
            Annotated[TrainerConfigInTheWild, tyro.conf.subcommand(name="in-the-wild")],
            Annotated[TrainerConfigLinear, tyro.conf.subcommand(name="linear-color")],
            Annotated[TrainerConfigSynthetic, tyro.conf.subcommand(name="synthetic")],
            Annotated[TrainerConfigMeshing, tyro.conf.subcommand(name="meshing")],
            Annotated[TrainerConfigAcademicBaseline, tyro.conf.subcommand(name="academic-baseline")],
        ]
        config = tyro.cli(Config)

    trainer = Trainer(config)

    # The viewer is the native one (src/app/webviewer/): its HTTP server and
    # render worker are C++ threads reading the session directly, so there is
    # no Python thread to start and nothing to push to it.
    if not config.disable_viewer:
        trainer.start_viewer(port=config.viewer_port)

    trainer.train()
    # trainer._train_with_profiling()
    trainer.eval()

    if not config.disable_viewer and config.keep_viewer_alive:
        print("Training and eval complete. Viewer still running -- press Ctrl-C to exit.")
        try:
            threading.Event().wait()
        except KeyboardInterrupt:
            print("\nShutting down...")
    trainer.stop_viewer()

if __name__ == "__main__":
    entrypoint()
