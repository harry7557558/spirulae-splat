import threading
import asyncio

from spirulae_splat.modules.trainer import SpirulaeSplatTrainerConfig, SpirulaeSplatTrainer

from spirulae_splat.viewer.server import ViewerServer, SliderDef, DropdownDef



async def start_viewer(trainer: SpirulaeSplatTrainer):

    server = ViewerServer(
        render_fn=trainer.render,
        http_host="localhost",
        http_port=7007,
        jpeg_quality=75,
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

async def main(trainer: SpirulaeSplatTrainer):
    asyncio.create_task(start_viewer(trainer))


def entrypoint():
    import tyro
    config = tyro.cli(SpirulaeSplatTrainerConfig)
    trainer = SpirulaeSplatTrainer(config)

    thread = threading.Thread(target=trainer.train, daemon=True)
    thread.start()

    asyncio.run(main(trainer))

if __name__ == "__main__":
    entrypoint()
