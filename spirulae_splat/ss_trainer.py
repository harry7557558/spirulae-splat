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


def entrypoint():
    import tyro

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
