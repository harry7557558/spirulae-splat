"""Multi-scene benchmark driver.

Each scene runs in its own `spirulae-train` subprocess so the C++ engine state
(world splats, bilagrid/PPISP/background, optimizer moments, color-space
matrices, device pool) starts fresh. Running multiple scenes inside one
process leaks engine state between scenes, which silently degrades metrics on
later runs even after `engine_reset()` (any Python-side state we forgot to
reset stays). Subprocesses sidestep that entirely.

Each subprocess writes its own metrics.json (now includes `training_time` and
`engine_vram`); we read it back here and aggregate.
"""

import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple

# Importing trainer configs only to read preset defaults (e.g. cap_max);
# the actual training runs in a subprocess, not in this process.
from spirulae_splat.modules.trainer import (
    TrainerConfig,
    TrainerConfigSquaredPos,
    TrainerConfigSquared,
    TrainerConfigPatched,
    TrainerConfigConfinedLowTexture,
    TrainerConfigConfined,
    TrainerConfigConfinedSquared,
    TrainerConfigOpenLowTexture,
    TrainerConfigOpen,
    TrainerConfigOpenSquared,
    TrainerConfigCenteredObject,
    TrainerConfigAcademicBaseline,
)


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


PRESET_TO_CLASS = {
    "3dgs": TrainerConfig,
    "3dgs^2-pos": TrainerConfigSquaredPos,
    "3dgs^2": TrainerConfigSquared,
    "3dgs-patched": TrainerConfigPatched,
    "3dgs-confined-low-texture": TrainerConfigConfinedLowTexture,
    "3dgs-confined": TrainerConfigConfined,
    "3dgs^2-confined": TrainerConfigConfinedSquared,
    "3dgs-open-low-texture": TrainerConfigOpenLowTexture,
    "3dgs-open": TrainerConfigOpen,
    "3dgs^2-open": TrainerConfigOpenSquared,
    "3dgs-centered-object": TrainerConfigCenteredObject,
    "academic-baseline": TrainerConfigAcademicBaseline,
}


def _preset_cap_max(preset: str) -> int:
    """Read the preset class's default model.cap_max so we can apply the
    zipnerf-specific ×3 multiplier without hardcoding the base value."""
    cls = PRESET_TO_CLASS[preset]
    # Instantiate with a placeholder data path; we only read defaults from it.
    inst = cls(data=Path("/dev/null"))
    return int(inst.model.cap_max)


def _run_scene(
    preset: str,
    data_dir: Path,
    output_prefix: Path,
    output_name: str,
    downscale: int,
    extra_args: List[str],
) -> Optional[Dict]:
    """Launch one training run in a fresh subprocess and return its metrics.

    Inherits stdout/stderr so the user sees per-scene training progress.
    Returns the parsed metrics.json dict, or None if the subprocess failed or
    metrics.json was not written.
    """
    metrics_dir = output_prefix / output_name
    metrics_path = metrics_dir / "metrics.json"
    if metrics_path.exists():
        # Stale leftover from a prior run with the same name would silently
        # masquerade as this run's results. Remove it up front.
        metrics_path.unlink()

    cmd = [
        sys.executable, "-m", "spirulae_splat.ss_trainer", preset,
        "--data", str(data_dir),
        "--dataparser.data-format", "colmap",
        "--dataparser.colmap-recon-dir", "sparse/0",
        "--dataparser.rescale-camera-to-fit", str(downscale),
        "--dataparser.image-dir", f"images_{downscale}",
        "--dataparser.eval-mode", "interval",
        "--save-eval-images",
        "--datamanager.no-load-depths",
        "--datamanager.no-load-normals",
        "--steps-per-save", "0",
        "--disable-viewer",
        "--output-dir-prefix", str(output_prefix),
        "--output-dir-name", output_name,
        # "--num-iterations", str(2500),
        # "--model.refine-stop-num-iter", str(0),
    ] + extra_args

    print(">>>", " ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        print(f"!!! subprocess for {output_name} exited with code {rc}", flush=True)
        return None
    if not metrics_path.exists():
        print(f"!!! subprocess for {output_name} did not write metrics.json", flush=True)
        return None
    with open(metrics_path, "r") as fp:
        return json.load(fp)


def bench_360_v2(preset: str, path_to_360_v2: Path, output_prefix: Path,
                 run_tag: str, extra_cli_args: List[str]) -> List[Tuple[str, Optional[Dict]]]:
    if not (path_to_360_v2.exists() and path_to_360_v2.is_dir()):
        print("Dataset not found. Please download from http://storage.googleapis.com/gresearch/refraw360/360_v2.zip and unzip.")
        exit(0)

    scenes = [
        # ("bicycle", 4),
        # ("garden", 4),
        # ("stump", 4),
        # ("bonsai", 2),
        # ("counter", 2),
        # ("kitchen", 2),
        ("room", 2),
    ]

    results = []
    for scene, downscale in scenes:
        print()
        print("Running:", scene)
        output_name = f"benchmark-{run_tag}-{scene}"
        extra_args = list(extra_cli_args)[:]
        extra_args.extend([
            "--dataparser.downscale-rounding-mode",
            "round" if scene == "garden" else "ceil"
        ])
        metrics = _run_scene(
            preset=preset,
            data_dir=path_to_360_v2 / scene,
            output_prefix=output_prefix,
            output_name=output_name,
            downscale=downscale,
            extra_args=extra_args,
        )
        results.append((scene, metrics))
    return results


def bench_zipnerf(preset: str, path_to_zipnerf: Path, output_prefix: Path,
                  run_tag: str, extra_cli_args: List[str]) -> List[Tuple[str, Optional[Dict]]]:
    if not (path_to_zipnerf.exists() and path_to_zipnerf.is_dir()):
        print("Dataset not found. Please download from https://smerf-3d.github.io/#data.")
        exit(0)

    scenes = [
        ("alameda", 4),
        ("berlin", 4),
        ("london", 4),
        ("nyc", 4),
    ]
    cap_max_3x = _preset_cap_max(preset) * 3  # 3M Gaussians vs preset default

    results = []
    for scene, downscale in scenes:
        print()
        print("Running:", scene)
        output_name = f"benchmark-{run_tag}-{scene}"
        metrics = _run_scene(
            preset=preset,
            data_dir=path_to_zipnerf / scene,
            output_prefix=output_prefix,
            output_name=output_name,
            downscale=downscale,
            extra_args=[
                # "--model.no-use-bilateral-grid-for-geometry",
                "--model.cap-max", str(cap_max_3x),
                "--dataparser.downscale-rounding-mode", "ceil",
            ] + list(extra_cli_args),
        )
        results.append((scene, metrics))
    return results


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

    output_dir_prefix: Path = Path("outputs")
    """Where each per-scene subprocess writes its checkpoint + metrics.json.
        The benchmark reads metrics.json back from `<prefix>/benchmark-<tag>-<scene>/`."""

    notify_on_complete_at: Optional[str] = None
    """Set this if you wish to receive a notification when training completes.
        Currently, this supports Discord webhook URL."""


def entrypoint():
    import datetime

    import tyro
    benchmark_config, extra_cli_args = tyro.cli(
        BenchmarkConfig, return_unknown_args=True
    )

    # Tag every scene from this benchmark run so its output dir is unique and
    # the metrics.json files don't collide with prior runs.
    run_tag = (f"{benchmark_config.benchmark}-{benchmark_config.preset}-"
               f"{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}")
    output_prefix = benchmark_config.output_dir_prefix

    if benchmark_config.benchmark == None:
        all_metrics = [
            (None, {"test": "If you see this message, notification is working."})
        ]
    elif benchmark_config.benchmark == "360_v2":
        all_metrics = bench_360_v2(benchmark_config.preset, benchmark_config.data,
                                    output_prefix, run_tag, extra_cli_args)
    elif benchmark_config.benchmark == "zipnerf":
        all_metrics = bench_zipnerf(benchmark_config.preset, benchmark_config.data,
                                     output_prefix, run_tag, extra_cli_args)
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
                value = f"{value:.5g}"
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
    if extra_cli_args:
        message += "Extra training args: `" + " ".join(extra_cli_args) + "`\n"
    if len(failed) > 0:
        message += "\n⚠️⚠️⚠️ Failed scene" + 's'*(len(failed)>1) + ": " + ', '.join([f"`{scene}`" for scene in failed]) + "\n"
    message += f"```\n{tabulate(all_data, headers, numalign='left')}\n```"

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
