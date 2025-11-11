#!/usr/bin/env python3

import numpy as np
import torch
from collections import OrderedDict
from plyfile import PlyData, PlyElement
from typing import Dict
import os

from spirulae_splat.viewer.model import SplatModel
from spirulae_splat.viewer.camera import Camera

from scipy.spatial.transform import Rotation

from spirulae_splat.viewer.utils import quat_scale_to_triangle_verts

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
    
    vertex_data = np.zeros(3*count, dtype=vertex_dtype)

    for key, tensor in tensors.items():
        if len(tensor.shape) == 1:
            vertex_data[key] = tensor
        else:
            vertex_data[key] = tensor.reshape(count, -1)
    
    vertex_element = PlyElement.describe(vertex_data, 'vertex')

    faces = np.array(
        [(x,) for x in np.arange(3*count).reshape(-1, 3).tolist()],
        dtype=[('vertex_indices', 'i4', (3,))]
    )
    face_element = PlyElement.describe(faces, 'face')

    ply_data = PlyData([vertex_element, face_element], text=False)
    
    ply_data.write(filename)


def export_ply(model: SplatModel, output_path: str) -> None:

    map_to_tensors = OrderedDict()

    positions, colors = quat_scale_to_triangle_verts(
        model.quats, model.scales, model.means,
        model.features_dc, model.features_ch
    )
    vert0, vert1, vert2 = torch.unbind(positions, dim=-2)
    normals = torch.cross(vert1-vert0, vert2-vert0, dim=-1)
    normals = normals / torch.norm(normals, dim=-1, keepdim=True).clip(min=1e-12)
    positions = positions.detach().cpu().numpy()
    normals = normals.detach().cpu().numpy()

    # position
    print("Encoding position...")
    n = len(model.means)
    map_to_tensors["x"] = positions[:, :, 0].flatten()
    map_to_tensors["y"] = positions[:, :, 1].flatten()
    map_to_tensors["z"] = positions[:, :, 2].flatten()
    # map_to_tensors["nx"] = normals[:, 0].repeat(3)
    # map_to_tensors["ny"] = normals[:, 1].repeat(3)
    # map_to_tensors["nz"] = normals[:, 2].repeat(3)

    # color
    print("Encoding color...")
    colors = torch.clamp(colors, 0.0, 1.0).cpu().numpy()
    colors = (colors * 255).astype(np.uint8)
    map_to_tensors["red"] = colors[:, :, 0].flatten()
    map_to_tensors["green"] = colors[:, :, 1].flatten()
    map_to_tensors["blue"] = colors[:, :, 2].flatten()

    # filter out inf/nan
    print("Filtering...")
    select = np.ones(n, dtype=bool)
    for k, t in map_to_tensors.items():
        n_before = np.sum(select)
        select = np.logical_and(select, np.isfinite(t.reshape(-1, 3).sum(-1)))
        n_after = np.sum(select)
        if n_after < n_before:
            print(f"{n_before - n_after} NaN/Inf elements in {k}")
    nan_count = n - np.sum(select)

    print(f"NaN/Inf: {nan_count}; Export: {np.sum(select)}/{n}")
    count = np.sum(select)
    if count < n:
        for k, t in map_to_tensors.items():
            map_to_tensors[k] = map_to_tensors[k].reshape(-1, 3)[select].flatten()
        count = np.sum(select)

    print("Writing PLY...")
    write_ply(output_path, count, map_to_tensors)


from export_ply_3dgs import export_equirectangular

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Export PLY (and equirectangular map), for triangle only.")
    parser.add_argument("work_dir", nargs=1, help="Path to the work folder (the one named YYYY-MM-DD_hhmmss).")
    parser.add_argument("--dataset_dir", help="Path to dataset folder.")
    # parser.add_argument("--output", "-o", default="splat.ply", help="The output PLY file.")
    args = parser.parse_args()
    work_dir = args.work_dir[0]
    args.output = "splat.ply"

    print("Work directory:", work_dir)

    print("Loading model...")
    model = SplatModel(work_dir)

    print("Orienting model...")
    model.convert_to_input_frame()

    if args.dataset_dir is not None and os.path.exists(os.path.join(args.dataset_dir, 'markers.yaml')):
        print()
        print("markers.yaml detected in dataset directory, attempt to align using markers")
        print()
        from spirulae_splat.viewer.align_apriltag import get_alignment
        rot, tr, sc = get_alignment(args.dataset_dir, verbose=True)
        model.change_frame(rot, tr, sc)
        print("Alignment complete.")

    print()

    print("Start PLY export")
    output_path = args.output
    if not ('/' in args.output or '\\' in args.output or os.path.sep in args.output):
        output_path = os.path.join(work_dir, output_path)
    export_ply(model, output_path)
    print("PLY saved to", output_path)

    print()

    sky_path = os.path.join(os.path.dirname(output_path), "background.png")

    print("Exporting sky...")
    export_equirectangular(model, sky_path)
    print("Sky map saved to", sky_path)
