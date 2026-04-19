import threading
import asyncio

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated

from spirulae_splat.viewer.server import ViewerServer, SliderDef, DropdownDef



async def start_viewer(trainer: Trainer):

    server = ViewerServer(
        render_fn=trainer.render,
        http_host="localhost",
        http_port=7007,
        # Example extra controls:
        # extra_sliders=[
        #     SliderDef(id="brightness", label="Brightness", min=0.5, max=2.0, step=0.05, value=1.0, unit="×"),
        #     SliderDef(id="anim_speed", label="Anim Speed", min=0.0, max=5.0, step=0.1, value=1.0, unit="×"),
        # ],
        # extra_dropdowns=[
        #     DropdownDef(
        #         id="shading_mode",
        #         label="Shading",
        #         options=[
        #             {"value": "full", "label": "Full"},
        #             {"value": "diffuse", "label": "Diffuse"},
        #             {"value": "specular", "label": "Specular"},
        #         ],
        #         default="full",
        #     ),
        # ],
        open_browser=False,
    )

    server.start()
    server.wait()

async def main(trainer: Trainer):
    asyncio.create_task(start_viewer(trainer))


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

    thread = threading.Thread(target=trainer.train, daemon=True)
    thread.start()

    asyncio.run(main(trainer))

if __name__ == "__main__":
    entrypoint()
