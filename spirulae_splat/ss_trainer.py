import threading
import asyncio

# import os
# os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
# os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
# os.environ["PYTORCH_NO_CUDA_MEMORY_CACHING"] = "1"

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated

from spirulae_splat.viewer.server import ViewerServer



async def start_viewer_server(trainer: Trainer):

    server = ViewerServer(
        render_fn=trainer.render,
        progress_fn=trainer.get_progress,
        pause_toggle_fn=trainer.toggle_pause,
        # Set the "render desired" flag from the HTTP handler's submit() call,
        # not later from the worker thread's render_fn invocation. The earlier
        # this flips, the more reliably the training loop's next-iteration
        # boundary will see it and yield the lock to the viewer. The worker
        # clears it once it has fully drained pending requests.
        on_render_submit=trainer._render_pending.set,
        on_render_idle=trainer._render_pending.clear,
        http_host="0.0.0.0",
        http_port=trainer.config.viewer_port,
        open_browser=False,
    )
    print()

    server.start()
    server.wait()

async def start_viewer(trainer: Trainer):
    await asyncio.create_task(start_viewer_server(trainer))


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
            config_from_json, resolve_checkpoint, apply_preset)
        run_dir, _ = resolve_checkpoint(Path(resume))
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

    if not config.disable_viewer:
        thread = threading.Thread(
            target=lambda: asyncio.run(start_viewer(trainer)),
            daemon=True
        )
        thread.start()

    trainer.train()
    # trainer._train_with_profiling()
    trainer.eval()

    if not config.disable_viewer and config.keep_viewer_alive:
        print("Training and eval complete. Viewer still running -- press Ctrl-C to exit.")
        try:
            threading.Event().wait()
        except KeyboardInterrupt:
            print("\nShutting down...")

if __name__ == "__main__":
    entrypoint()
