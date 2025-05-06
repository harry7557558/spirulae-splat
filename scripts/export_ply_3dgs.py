import numpy as np
import torch
from collections import OrderedDict
from plyfile import PlyData, PlyElement
from typing import Dict

from spirulae_splat.viewer.model import SplatModel


def write_ply(
    filename: str,
    count: int,
    tensors: Dict[str, np.ndarray],
):

    vertex_dtype = []
    for key, tensor in tensors.items():
        dtype = 'f4' if tensor.dtype.kind == 'f' else 'u1'
        shape = () if len(tensor.shape) == 1 else tensor.shape[1:]
        vertex_dtype.append((key, dtype, shape))
    
    vertex_data = np.zeros(count, dtype=vertex_dtype)

    for key, tensor in tensors.items():
        if len(tensor.shape) == 1:
            vertex_data[key] = tensor
        else:
            vertex_data[key] = tensor.reshape(count, -1)
    
    vertex_element = PlyElement.describe(vertex_data, 'vertex')
    
    ply_data = PlyData([vertex_element], text=False)
    
    ply_data.write(filename)


def process_model(model: SplatModel, output_path: str) -> None:

    map_to_tensors = OrderedDict()

    # position
    positions = model.means.cpu().numpy()
    count = positions.shape[0]
    n = count
    map_to_tensors["x"] = positions[:, 0]
    map_to_tensors["y"] = positions[:, 1]
    map_to_tensors["z"] = positions[:, 2]
    map_to_tensors["nx"] = np.zeros(n, dtype=np.float32)
    map_to_tensors["ny"] = np.zeros(n, dtype=np.float32)
    map_to_tensors["nz"] = np.zeros(n, dtype=np.float32)

    # color
    if False:
        colors = torch.clamp(model.colors, 0.0, 1.0).cpu().numpy()
        colors = (colors * 255).astype(np.uint8)
        map_to_tensors["red"] = colors[:, 0]
        map_to_tensors["green"] = colors[:, 1]
        map_to_tensors["blue"] = colors[:, 2]

    features_dc = model.features_dc.cpu().numpy()
    for i in range(3):
        map_to_tensors[f"f_dc_{i}"] = features_dc[:, i]

    # SH
    if model.sh_degree > 0:
        shs_rest = model.features_sh.transpose(1, 2).contiguous().cpu().numpy()
        shs_rest = shs_rest.reshape((n, -1))
        for i in range(shs_rest.shape[-1]):
            map_to_tensors[f"f_rest_{i}"] = shs_rest[:, i]

    # opac, scale, quat
    map_to_tensors["opacity"] = model.opacities.cpu().numpy().squeeze(-1)

    scales = model.scales.data.cpu().numpy()
    for i in range(3):
        map_to_tensors[f"scale_{i}"] = scales[:, i]

    quats = model.quats.data.cpu().numpy()
    for i in range(4):
        map_to_tensors[f"rot_{i}"] = quats[:, i]

    # filter out inf/nan
    select = np.ones(n, dtype=bool)
    for k, t in map_to_tensors.items():
        n_before = np.sum(select)
        select = np.logical_and(select, np.isfinite(t))
        n_after = np.sum(select)
        if n_after < n_before:
            print(f"{n_before - n_after} NaN/Inf elements in {k}")
    nan_count = np.sum(select) - n

    # filter out low opacity
    low_opacity_gaussians = map_to_tensors["opacity"] < -5.5373  # logit(1/255)
    lowopa_count = np.sum(low_opacity_gaussians)
    select[low_opacity_gaussians] = 0

    print(f"{nan_count} Gaussians have NaN/Inf and {lowopa_count} have low opacity, only export {np.sum(select)}/{n}")
    if np.sum(select) < n:
        for k, t in map_to_tensors.items():
            map_to_tensors[k] = map_to_tensors[k][select]
        count = np.sum(select)

    write_ply(output_path, count, map_to_tensors)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Export PLY, for 3DGS only.")
    parser.add_argument("input_file", nargs=1, help="The input config file.")
    parser.add_argument("--output", "-o", default="ssplat.ply", help="The output PLY file.")
    args = parser.parse_args()

    model = SplatModel(args.input_file[0])
    process_model(model, args.output)
