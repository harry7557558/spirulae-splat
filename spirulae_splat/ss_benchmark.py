import threading
import asyncio

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated, Literal
import json



def bench_360_v2(config_class, path_to_360_v2: Path):
    if not (path_to_360_v2.exists() and path_to_360_v2.is_dir()):
        print("Dataset not found. Please download from http://storage.googleapis.com/gresearch/refraw360/360_v2.zip and unzip.")
        exit(0)

    all_metrics = []

    for scene, downscale in [
        ("bicycle", 4),
        ("garden", 4),
        ("stump", 4),
        ("bonsai", 2),
        ("counter", 2),
        ("kitchen", 2),
        ("room", 2),
    ]:
        print()
        print("Running:", scene)

        config = config_class(data=path_to_360_v2 / scene)  # type: TrainerConfig
        config.dataparser.data_format = "colmap"
        config.dataparser.colmap_recon_dir = Path("sparse/0")
        config.dataparser.rescale_camera_to_fit = downscale
        config.dataparser.image_dir = f"images_{downscale}"
        config.dataparser.eval_mode = "interval"

        # TODO: color correct
        config.model.use_bilateral_grid = False
        config.model.use_ppisp = False
        config.datamanager.load_depths = False
        config.datamanager.load_normals = False

        config.num_iterations = 1000

        from time import perf_counter
        try:
            trainer = Trainer(config)
            time0 = perf_counter()
            trainer.train()
            time1 = perf_counter()
            trainer.eval()
        except:
            import traceback
            traceback.print_exc()
            all_metrics.append((scene, None))
            continue

        with open(trainer.output_dir / "metrics.json", "r") as fp:
            metrics = json.load(fp)
        metrics['training_time'] = time1-time0
        # TODO: VRAM
        all_metrics.append((scene, metrics))

        del trainer
        torch.cuda.empty_cache()

    return all_metrics


@dataclass
class BenchmarkConfig:
    benchmark: Literal["360_v2", "tankt_db", "zipnerf", "blender", "llff", None]
    """Which benchmark to run.
        360_v2: 7 permissively released scenes from Mip-NeRF 360 dataset
        tankt_db: 2 scenes from Tank & Temples dataset, and 2 scenes from Deep Blending dataset
        zipnerf: 4 fisheye scenes from Zip-NeRF dataset
        None: run nothing, can be used to test if notification is working
    """

    preset: Literal[
        "3dgs", "3dgs^2-pos", "3dgs^2",
        "3dgs-confined-low-texture", "3dgs-confined", "3dgs^2-confined",
        "3dgs-open-low-texture", "3dgs-open", "3dgs^2-open",
        "3dgs-centered-object", "academic-baseline"
    ]
    """Preset to use to run the benchmark"""

    data: Path
    """Path to folder containing benchmark dataset."""

    notify_on_complete_at: Optional[str] = None
    """Set this if you wish to receive a notification when training completes.
        Currently, this supports Discord webhook URL."""


def entrypoint():

    import tyro
    benchmark_config = tyro.cli(BenchmarkConfig)

    config_class = {
        "3dgs": TrainerConfig,
        "3dgs^2-pos": TrainerConfigSquaredPos,
        "3dgs^2": TrainerConfigSquared,
        "3dgs-patched": TrainerConfigPatched,
        "triangle": TrainerConfigTriangle,
        "triangle-patched": TrainerConfigTrianglePatched,
        "voxel": TrainerConfigVoxel,
        "3dgs-confined-low-texture": TrainerConfigConfinedLowTexture,
        "3dgs-confined": TrainerConfigConfined,
        "3dgs^2-confined": TrainerConfigConfinedSquared,
        "3dgs-open-low-texture": TrainerConfigOpenLowTexture,
        "3dgs-open": TrainerConfigOpen,
        "3dgs^2-open": TrainerConfigOpenSquared,
        "3dgs-centered-object": TrainerConfigCenteredObject,
        "academic-baseline": TrainerConfigAcademicBaseline,
    }[benchmark_config.preset]

    if benchmark_config.benchmark == None:
        all_metrics = [
            (None, {"test": "If you see this message, notification is working."})
        ]
    elif benchmark_config.benchmark == "360_v2":
        all_metrics = bench_360_v2(config_class, benchmark_config.data)
    else:
        raise NotImplementedError()

    headers = []
    all_data = []
    failed = []

    for column, (scene_name, metrics) in enumerate(all_metrics):
        if scene_name is not None:
            print()
            print(scene_name)
        else:
            scene_name = "Notification Test"
        if metrics is None:
            print("Failed.")
            failed.append(scene_name)
            continue
        headers.append(scene_name)
        row = 0
        for key, value in metrics.items():
            if key.startswith("avg_"):
                key = key[len("avg_"):]
            if isinstance(value, list):
                continue
            if column == 0:
                all_data.append([key])
            if isinstance(value, float):
                value = f"{value:.7g}"
            all_data[row].append(value)
            if key != "test":
                print(key, value)
            row += 1
    print()

    if benchmark_config.notify_on_complete_at is None:
        return

    from tabulate import tabulate
    from urllib.request import Request, urlopen

    message = f"## Benchmark Complete (`{benchmark_config.benchmark}`, `{benchmark_config.preset}`)\n"
    message += f"Your benchmark has completed. Here are the results:\n"
    if len(failed) > 0:
        message += "\n⚠️⚠️⚠️ Failed scene" + 's'*(len(failed)>1) + ": " + ', '.join([f"`{scene}`" for scene in failed]) + "\n"
    message += f"```\n{tabulate(all_data, headers, numalign="left")}\n```"

    payload = {
        "content": message,
        "username": "spirulae-splat",
        "avatar_url": "https://avatars.githubusercontent.com/u/44143868",
        "attachments": []
    }
    data = json.dumps(payload).encode('utf-8')

    req = Request(benchmark_config.notify_on_complete_at, data=data, method='POST')
    req.add_header('Content-Type', 'application/json')
    # req.add_header('User-Agent', 'Python-urllib/3.x')
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:149.0) Gecko/20100101 Firefox/149.0')

    try:
        with urlopen(req) as response:
            print(f"Notification sent.")
    except Exception as e:
        print(f"Failed to send notification: {e}")

if __name__ == "__main__":
    entrypoint()
