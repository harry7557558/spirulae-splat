import threading
import asyncio

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated

from spirulae_splat.viewer.server import ViewerServer



async def start_viewer_server(trainer: Trainer):

    server = ViewerServer(
        render_fn=trainer.render,
        progress_fn=trainer.get_progress,
        http_host="0.0.0.0",
        http_port=trainer.config.viewer_port,
        open_browser=False,
    )

    server.start()
    server.wait()

async def start_viewer(trainer: Trainer):
    await asyncio.create_task(start_viewer_server(trainer))


def entrypoint():
    import tyro

    Config = Union[
        Annotated[TrainerConfig, tyro.conf.subcommand(name="3dgs")],
        Annotated[TrainerConfigSquaredPos, tyro.conf.subcommand(name="3dgs^2-pos")],
        Annotated[TrainerConfigSquared, tyro.conf.subcommand(name="3dgs^2")],
        Annotated[TrainerConfigPatched, tyro.conf.subcommand(name="3dgs-patched")],
        Annotated[TrainerConfigTriangle, tyro.conf.subcommand(name="triangle")],
        Annotated[TrainerConfigTrianglePatched, tyro.conf.subcommand(name="triangle-patched")],
        Annotated[TrainerConfigVoxel, tyro.conf.subcommand(name="voxel")],
        Annotated[TrainerConfigConfinedLowTexture, tyro.conf.subcommand(name="3dgs-confined-low-texture")],
        Annotated[TrainerConfigConfined, tyro.conf.subcommand(name="3dgs-confined")],
        Annotated[TrainerConfigConfinedSquared, tyro.conf.subcommand(name="3dgs^2-confined")],
        Annotated[TrainerConfigOpenLowTexture, tyro.conf.subcommand(name="3dgs-open-low-texture")],
        Annotated[TrainerConfigOpen, tyro.conf.subcommand(name="3dgs-open")],
        Annotated[TrainerConfigOpenSquared, tyro.conf.subcommand(name="3dgs^2-open")],
        Annotated[TrainerConfigCenteredObject, tyro.conf.subcommand(name="3dgs-centered-object")],
        Annotated[TrainerConfigAcademicBaseline, tyro.conf.subcommand(name="academic-baseline")],
    ]

    config = tyro.cli(Config)
    trainer = Trainer(config)

    thread = threading.Thread(
        target=lambda: asyncio.run(start_viewer(trainer)),
        daemon=True
    )
    thread.start()

    trainer.train()
    trainer.eval()

if __name__ == "__main__":
    entrypoint()
