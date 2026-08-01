import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import torch
import tyro
from plyfile import PlyData

from spirulae_splat.modules.core import Renderer
from spirulae_splat.splat.cuda import _C


@dataclass
class ViewerConfig:
    ply: Path
    viewer_port: int = 7007
    primitive: Literal["3dgs", "mip", "3dgut"] = "3dgut"
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    open_browser: bool = False


class PlyViewer:
    def __init__(self, config: ViewerConfig):
        self.config = config
        self.device = torch.device(config.device)
        if self.device.type != "cuda":
            raise RuntimeError("SS viewer currently requires CUDA; please run with a CUDA-enabled PyTorch build.")

        self.splats_world, self.cur_num_splats = self._load_ply(config.ply, config.primitive)
        self.renderer = Renderer(config.primitive, self.splats_world, self.cur_num_splats)

    def _load_ply(self, ply_file_path: Path, primitive: str):
        plydata = PlyData.read(str(ply_file_path))
        if "vertex" not in plydata:
            raise ValueError(f"PLY file {ply_file_path} does not contain a vertex element.")

        vertex = plydata["vertex"]
        names = vertex.data.dtype.names
        count = vertex.count
        if count == 0:
            raise ValueError(f"PLY file {ply_file_path} contains no vertices.")

        def _has(*fields):
            return all(field in names for field in fields)

        def _read(*fields):
            return [np.asarray(vertex[field], dtype=np.float32) for field in fields]

        means = torch.from_numpy(np.stack(_read("x", "y", "z"), axis=-1)).contiguous().to(self.device)

        if _has("rot_0", "rot_1", "rot_2", "rot_3"):
            quats = torch.from_numpy(np.stack(_read("rot_0", "rot_1", "rot_2", "rot_3"), axis=-1))
        else:
            quats = torch.tensor([1.0, 0.0, 0.0, 0.0], dtype=torch.float32).repeat(count, 1)
        quats = quats.contiguous().to(self.device)

        if _has("scale_0", "scale_1", "scale_2"):
            scales = torch.from_numpy(np.stack(_read("scale_0", "scale_1", "scale_2"), axis=-1))
        else:
            scales = torch.ones((count, 3), dtype=torch.float32)
        scales = scales.contiguous().to(self.device)

        if _has("opacity"):
            opacities = torch.from_numpy(np.asarray(vertex["opacity"], dtype=np.float32)).unsqueeze(-1)
        elif _has("alpha"):
            opacities = torch.from_numpy(np.asarray(vertex["alpha"], dtype=np.float32)).unsqueeze(-1)
        else:
            opacities = torch.ones((count, 1), dtype=torch.float32)
        opacities = opacities.contiguous().to(self.device)

        if _has("f_dc_0", "f_dc_1", "f_dc_2"):
            features_dc = torch.from_numpy(np.stack(_read("f_dc_0", "f_dc_1", "f_dc_2"), axis=-1))
        elif _has("red", "green", "blue"):
            rgb = np.stack(_read("red", "green", "blue"), axis=-1)
            features_dc = torch.from_numpy(rgb.astype(np.float32) / 255.0)
        else:
            features_dc = torch.zeros((count, 3), dtype=torch.float32)
        features_dc = features_dc.contiguous().to(self.device)

        sh_keys = sorted([name for name in names if name.startswith("f_rest_")], key=lambda x: int(x.split("_")[-1]))
        if sh_keys:
            sh_values = np.stack([np.asarray(vertex[key], dtype=np.float32) for key in sh_keys], axis=-1)
            if sh_values.shape[1] % 3 != 0:
                raise ValueError("PLY SH fields must encode 3 channels in contiguous f_rest_<i> fields.")
            n_channels = sh_values.shape[1] // 3
            sh_values = sh_values.reshape(count, 3, n_channels)
            features_sh = torch.from_numpy(sh_values).permute(0, 2, 1)
        else:
            features_sh = torch.zeros((count, 0, 3), dtype=torch.float32)
        features_sh = features_sh.contiguous().to(self.device)

        splats_world = (means, quats, scales, opacities, features_dc, features_sh)

        return splats_world, count


def entrypoint_body(config: ViewerConfig) -> None:
    """Serve a PLY through the native web viewer.

    The Renderer here exists only to upload the PLY into the engine world; the
    rendering, HTTP serving and the browser client are all the native path
    (src/app/webviewer/), the same one `ssplat train` and the GUI use. An
    empty PostSplitCameras means the viewer has no training cameras to draw --
    a bare PLY has none.
    """
    viewer_obj = PlyViewer(config)   # uploads the splats via Renderer

    cfg = _C.ViewerRenderConfig()
    cfg.primitive = config.primitive

    server = _C.WebViewer()
    server.start("0.0.0.0", config.viewer_port, cfg, _C.PostSplitCameras())
    url = f"http://0.0.0.0:{config.viewer_port}/"
    print(f"Viewer at {url}")
    if config.open_browser:
        import webbrowser
        threading.Timer(0.5, webbrowser.open, args=[url]).start()
    try:
        threading.Event().wait()
    except KeyboardInterrupt:
        print("\nShutting down...")
    finally:
        server.stop()


def entrypoint() -> None:
    entrypoint_body(tyro.cli(ViewerConfig))


if __name__ == "__main__":
    entrypoint()
