import threading
import asyncio

from spirulae_splat.modules.trainer import *
from typing import Union, Annotated, Literal
import json


# Mip-NeRF 360 download command
"""
mkdir 360_v2 && cd 360_v2
wget http://storage.googleapis.com/gresearch/refraw360/360_v2.zip
unzip 360_v2.zip
rm 360_v2.zip
cd ..
"""

# Zip-NeRF download command
"""
mkdir zipnerf && cd zipnerf
wget https://storage.googleapis.com/gresearch/refraw360/zipnerf/alameda.zip
unzip alameda.zip
rm alameda.zip
wget https://storage.googleapis.com/gresearch/refraw360/zipnerf/berlin.zip
unzip berlin.zip
rm berlin.zip
wget https://storage.googleapis.com/gresearch/refraw360/zipnerf/london.zip
unzip london.zip
rm london.zip
wget https://storage.googleapis.com/gresearch/refraw360/zipnerf/nyc.zip
unzip nyc.zip
rm nyc.zip
cd ..
"""


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
        config.dataparser.rescale_camera_to_fit = downscale  # consistent with gsplat/Inria
        # config.dataparser.rescale_camera_to_fit = True
        config.dataparser.image_dir = f"images_{downscale}"
        config.dataparser.eval_mode = "interval"

        # TODO: color correct
        config.model.use_bilateral_grid = False
        config.model.use_ppisp = False
        config.datamanager.load_depths = False
        config.datamanager.load_normals = False
        config.model.use_bilateral_grid_for_geometry = False

        config.steps_per_save = 0

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


def bench_zipnerf(config_class, path_to_zipnerf: Path):
    if not (path_to_zipnerf.exists() and path_to_zipnerf.is_dir()):
        print("Dataset not found. Please download from https://smerf-3d.github.io/#data.")
        exit(0)

    all_metrics = []

    for scene, downscale in [
        ("alameda", 4),
        ("berlin", 4),
        ("london", 4),
        ("nyc", 4),
    ]:
        print()
        print("Running:", scene)

        config = config_class(data=path_to_zipnerf / scene)  # type: TrainerConfig
        config.dataparser.data_format = "colmap"
        config.dataparser.colmap_recon_dir = Path("sparse/0")
        config.dataparser.rescale_camera_to_fit = downscale
        # config.dataparser.rescale_camera_to_fit = True
        config.dataparser.image_dir = f"images_{downscale}"
        config.dataparser.eval_mode = "interval"

        # TODO: color correct
        config.model.use_bilateral_grid = False
        config.model.use_ppisp = False
        config.datamanager.load_depths = False
        config.datamanager.load_normals = False
        config.model.use_bilateral_grid_for_geometry = False

        config.steps_per_save = 0
        config.model.cap_max *= 3  # 3M Gaussians

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
    elif benchmark_config.benchmark == "zipnerf":
        all_metrics = bench_zipnerf(config_class, benchmark_config.data)
    else:
        raise NotImplementedError()

    headers = []
    all_data = []
    failed = []

    column = 0
    for scene_name, metrics in all_metrics:
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
        column += 1
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


# if __name__ == "__main__":
#     entrypoint()


def plot_results(path_to_report: str):
    if path_to_report.startswith("http"):
        import urllib.request
        with urllib.request.urlopen(path_to_report) as response:
            content = response.read().decode('utf-8')
    else:
        content = open(path_to_report).read()
    content = content.split('\n----\n')
    print(content)

plot_results("https://gist.githubusercontent.com/harry7557558/1709047429643a6bcd6a10fd0d76fb6a/raw/0e779a6904c1fb9e0d7a8ededc4b3cfa8d7e410a/260425.md")
