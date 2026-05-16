import asyncio
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import torch
import tyro
from plyfile import PlyData

from spirulae_splat.modules.core import Renderer
from spirulae_splat.viewer.server import ViewerServer


@dataclass
class ViewerConfig:
    ply: Path
    viewer_port: int = 7007
    primitive: Literal["3dgs", "mip", "3dgut", "3dgut_sv", "opaque_triangle", "voxel"] = "3dgut"
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

        means = torch.from_numpy(np.stack(_read("x", "y", "z"), axis=-1)).to(self.device).contiguous()

        if _has("rot_0", "rot_1", "rot_2", "rot_3"):
            quats = torch.from_numpy(np.stack(_read("rot_0", "rot_1", "rot_2", "rot_3"), axis=-1))
        else:
            quats = torch.tensor([1.0, 0.0, 0.0, 0.0], dtype=torch.float32).repeat(count, 1)
        quats = quats.to(self.device).contiguous()

        if _has("scale_0", "scale_1", "scale_2"):
            scales = torch.from_numpy(np.stack(_read("scale_0", "scale_1", "scale_2"), axis=-1))
        else:
            scales = torch.ones((count, 3), dtype=torch.float32)
        scales = scales.to(self.device).contiguous()

        if _has("opacity"):
            opacities = torch.from_numpy(np.asarray(vertex["opacity"], dtype=np.float32)).unsqueeze(-1)
        elif _has("alpha"):
            opacities = torch.from_numpy(np.asarray(vertex["alpha"], dtype=np.float32)).unsqueeze(-1)
        else:
            opacities = torch.ones((count, 1), dtype=torch.float32)
        opacities = opacities.to(self.device).contiguous()

        if _has("f_dc_0", "f_dc_1", "f_dc_2"):
            features_dc = torch.from_numpy(np.stack(_read("f_dc_0", "f_dc_1", "f_dc_2"), axis=-1))
        elif _has("red", "green", "blue"):
            rgb = np.stack(_read("red", "green", "blue"), axis=-1)
            features_dc = torch.from_numpy(rgb.astype(np.float32) / 255.0)
        else:
            features_dc = torch.zeros((count, 3), dtype=torch.float32)
        features_dc = features_dc.to(self.device).contiguous()

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
        features_sh = features_sh.to(self.device).contiguous()

        if primitive == "opaque_triangle":
            features_ch = torch.zeros((count, 2, 3), dtype=torch.float32, device=self.device)
            splats_world = (means, quats, scales, opacities, features_dc, features_sh, features_ch)
        else:
            splats_world = (means, quats, scales, opacities, features_dc, features_sh)

        return splats_world, count

    def _camera_to_viewmat(self, c2w: np.ndarray) -> torch.Tensor:
        c2w = torch.from_numpy(c2w.astype(np.float32)).to(self.device)
        if c2w.shape != (3, 4):
            raise ValueError("Camera-to-world matrix must be 3x4.")

        R = c2w[:3, :3]
        T = c2w[:3, 3:4]
        R = R * torch.tensor([1.0, -1.0, -1.0], dtype=torch.float32, device=self.device)[None, :]
        R_inv = R.transpose(-1, -2)
        T_inv = -torch.matmul(R_inv, T)

        viewmat = torch.eye(4, dtype=torch.float32, device=self.device)
        viewmat[:3, :3] = R_inv
        viewmat[:3, 3:4] = T_inv
        return viewmat.unsqueeze(0)

    def render(self, c2w, fx, fy, cx, cy, width, height, camera_model):
        camera_model = camera_model.upper()
        viewmats = self._camera_to_viewmat(c2w)
        intrins = torch.tensor([[fx, fy, cx, cy]], dtype=torch.float32, device=self.device)

        self.renderer.set_params(
            viewmats=viewmats,
            intrins=intrins,
            width=int(width),
            height=int(height),
            packed=True,
            use_bvh=False,
            camera_model=camera_model,
        )
        self.renderer.forward()

        rgb = self.renderer.render_colors[0]
        return {
            "rgb": rgb[0],
            "_post_processor": lambda tensor, **kwargs: (255*torch.clamp(tensor, 0.0, 1.0)).to(torch.uint8),
        }


async def start_viewer_server(viewer: PlyViewer) -> None:
    server = ViewerServer(
        render_fn=viewer.render,
        progress_fn=None,
        http_host="0.0.0.0",
        http_port=viewer.config.viewer_port,
        open_browser=viewer.config.open_browser,
    )
    server.start()
    server.wait()


async def start_viewer(viewer: PlyViewer) -> None:
    await asyncio.create_task(start_viewer_server(viewer))


def entrypoint() -> None:
    config = tyro.cli(ViewerConfig)
    viewer = PlyViewer(config)
    thread = __import__("threading").Thread(
        target=lambda: asyncio.run(start_viewer(viewer)),
        daemon=False,
    )
    thread.start()
    thread.join()


if __name__ == "__main__":
    entrypoint()
