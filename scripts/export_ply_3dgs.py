import numpy as np
import torch
from collections import OrderedDict
from plyfile import PlyData, PlyElement
from typing import Dict
import os

from spirulae_splat.viewer.model import SplatModel
from spirulae_splat.viewer.camera import Camera

from scipy.spatial.transform import Rotation


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


def generate_equirectangular(
        render_fun, out_w=1024, out_h=512, face_res=512
        ):
    """
    Renders a 1024×512 equirectangular map by:
      1) rendering 6 cube faces at face_res×face_res, each with FoV=90°,
      2) remapping into equirectangular projection via vectorized bilinear interp.
    """

    # Prepare intrinsics for a 90° FoV pinhole camera at face_res
    fx = fy = face_res / 2.0
    cx = cy = face_res / 2.0
    face_intr = (face_res, face_res, fx, fy, cx, cy)

    # Define cube‐map face rotations (local → world)
    face_rotations = {
        'front': Rotation.from_rotvec([0, 0, 0]).as_matrix(),
        'right': Rotation.from_rotvec([0, np.pi/2, 0]).as_matrix(),
        'back':  Rotation.from_rotvec([0, np.pi, 0]).as_matrix(),
        'left':  Rotation.from_rotvec([0, -np.pi/2, 0]).as_matrix(),
        'up':    Rotation.from_rotvec([-np.pi/2, 0, 0]).as_matrix(),
        'down':  Rotation.from_rotvec([np.pi/2, 0, 0]).as_matrix(),
    }

    # 1) Render all six faces
    faces = {}
    for name, c2w in face_rotations.items():
        faces[name] = render_fun(face_intr, c2w)  # each → (face_res, face_res, 3)

    # 2) Build equirectangular map
    #   θ ∈ [−π, π), φ ∈ [ π/2, -π/2]
    u = np.linspace(-np.pi, np.pi,     out_w, endpoint=False)
    v = np.linspace( np.pi/2, -np.pi/2, out_h)
    θ, φ = np.meshgrid(u, v)

    # Spherical → Cartesian (world‐space unit directions)
    x = np.cos(φ) * np.sin(θ)
    y = np.sin(φ)
    z = np.cos(φ) * np.cos(θ)

    # Prepare output
    eq = np.zeros((out_h, out_w, 3), dtype=np.float32)

    # Vectorized bilinear sampler
    def sample(img, uvx, uvy):
        """
        img: H×W×3, uvx,uvy in [0,1]
        returns H×W×3 sampled bilinearly
        """
        H, W, _ = img.shape
        # Map [0,1]→[0, W-1] and [0, H-1]
        xs = uvx * (W - 1) + 0.5
        ys = uvy * (H - 1) + 0.5

        x0 = np.clip(np.floor(xs).astype(int), 0, W-1)
        x1 = np.clip(x0 + 1, 0, W-1)
        y0 = np.clip(np.floor(ys).astype(int), 0, H-1)
        y1 = np.clip(y0 + 1, 0, H-1)

        wa = np.maximum(x1 - xs, 0) * np.maximum(y1 - ys, 0)
        wb = np.maximum(xs - x0, 0) * np.maximum(y1 - ys, 0)
        wc = np.maximum(x1 - xs, 0) * np.maximum(ys - y0, 0)
        wd = np.maximum(xs - x0, 0) * np.maximum(ys - y0, 0)

        Ia = img[y0, x0]
        Ib = img[y0, x1]
        Ic = img[y1, x0]
        Id = img[y1, x1]

        return (wa[...,None] * Ia +
                wb[...,None] * Ib +
                wc[...,None] * Ic +
                wd[...,None] * Id) / (wa+wb+wc+wd)[..., None]

    # Decide which face each direction belongs to
    abs_x, abs_y, abs_z = np.abs(x), np.abs(y), np.abs(z)
    eps = 1e-3

    # +X face
    mask = (abs_x > abs_y-eps) & (abs_x > abs_z-eps) & (x > -eps)
    #   u =  -z/x, v =   y/x
    uvx = (-z[mask] / x[mask] + 1) * 0.5
    uvy = (  y[mask] / x[mask] + 1) * 0.5
    eq[mask] = sample(faces['right'], uvx, uvy)

    # -X face
    mask = (abs_x > abs_y-eps) & (abs_x > abs_z-eps) & (x < eps)
    #   u =   z/(-x), v =   y/(-x)
    uvx = ( z[mask] / -x[mask] + 1) * 0.5
    uvy = ( y[mask] / -x[mask] + 1) * 0.5
    eq[mask] = sample(faces['left'], uvx, uvy)

    # +Y face (up)
    mask = (abs_y > abs_x-eps) & (abs_y > abs_z-eps) & (y > -eps)
    #   u =   x/y, v =  -z/y
    uvx = ( x[mask] / y[mask] + 1) * 0.5
    uvy = (-z[mask] / y[mask] + 1) * 0.5
    eq[mask] = sample(faces['up'], uvx, uvy)

    # -Y face (down)
    mask = (abs_y > abs_x-eps) & (abs_y > abs_z-eps) & (y < eps)
    #   u =   x/(-y), v =   z/(-y)
    uvx = ( x[mask] / -y[mask] + 1) * 0.5
    uvy = ( z[mask] / -y[mask] + 1) * 0.5
    eq[mask] = sample(faces['down'], uvx, uvy)

    # +Z face (front)
    mask = (abs_z > abs_x-eps) & (abs_z > abs_y-eps) & (z > -eps)
    #   u =   x/z, v =   y/z
    uvx = ( x[mask] / z[mask] + 1) * 0.5
    uvy = ( y[mask] / z[mask] + 1) * 0.5
    eq[mask] = sample(faces['front'], uvx, uvy)

    # -Z face (back)
    mask = (abs_z > abs_x-eps) & (abs_z > abs_y-eps) & (z < eps)
    #   u =  -x/(-z), v =   y/(-z)   i.e. u = x/z, v = y/(-z)? but flip
    uvx = ( -x[mask] / -z[mask] + 1) * 0.5
    uvy = (  y[mask] / -z[mask] + 1) * 0.5
    eq[mask] = sample(faces['back'], uvx, uvy)

    eq = np.clip(eq, 0.0, 1.0)
    return eq


def export_ply(model: SplatModel, output_path: str) -> None:

    map_to_tensors = OrderedDict()

    # position
    print("Encoding position...")
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
    print("Encoding color...")
    if True:
        colors = torch.clamp(model.colors, 0.0, 1.0).cpu().numpy()
        colors = (colors * 255).astype(np.uint8)
        map_to_tensors["red"] = colors[:, 0]
        map_to_tensors["green"] = colors[:, 1]
        map_to_tensors["blue"] = colors[:, 2]

    features_dc = model.features_dc.cpu().numpy()
    for i in range(3):
        map_to_tensors[f"f_dc_{i}"] = features_dc[:, i]

    # SH
    print("Encoding SH...")
    if model.sh_degree > 0:
        shs_rest = model.features_sh.transpose(1, 2).contiguous().cpu().numpy()
        shs_rest = shs_rest.reshape((n, -1))
        for i in range(shs_rest.shape[-1]):
            map_to_tensors[f"f_rest_{i}"] = shs_rest[:, i]

    # opac, scale, quat
    print("Encoding opacity, scale, and rotation...")
    map_to_tensors["opacity"] = model.opacities.cpu().numpy().squeeze(-1)

    scales = model.scales.data.cpu().numpy()
    for i in range(3):
        map_to_tensors[f"scale_{i}"] = scales[:, i]

    quats = model.quats.data.cpu().numpy()
    for i in range(4):
        map_to_tensors[f"rot_{i}"] = quats[:, i]

    # filter out inf/nan
    print("Filtering...")
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

    print(f"NaN/Inf: {nan_count}; Low opacity: {lowopa_count}; Export: {np.sum(select)}/{n}")
    if np.sum(select) < n:
        for k, t in map_to_tensors.items():
            map_to_tensors[k] = map_to_tensors[k][select]
        count = np.sum(select)

    print("Writing PLY...")
    write_ply(output_path, count, map_to_tensors)


def export_equirectangular(model: SplatModel, output_path: str) -> None:
    camera_path = os.path.join(os.path.dirname(__file__), "../spirulae_splat/viewer/cameras/s21.yaml")
    camera = Camera(camera_path)

    def render(intrins, c2w):
        camera.w, camera.h, camera.fx, camera.fy, camera.cx, camera.cy = intrins
        return model._get_background_image(camera, c2w).cpu().numpy()

    im = generate_equirectangular(render)

    import imageio
    imageio.imwrite(output_path, (im * 255).astype(np.uint8))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Export PLY, for 3DGS only.")
    parser.add_argument("input_file", nargs=1, help="The input config file.")
    parser.add_argument("--output", "-o", default="splat.ply", help="The output PLY file.")
    args = parser.parse_args()

    print("Loading model...")
    model = SplatModel(args.input_file[0])

    print("Orienting model...")
    model.convert_to_input_frame()

    print()

    print("Start PLY export")
    export_ply(model, args.output)
    print("PLY saved to", args.output)

    print()

    sky_path = os.path.join(os.path.dirname(args.output), "background.png")

    print("Exporting sky...")
    export_equirectangular(model, sky_path)
    print("Sky map saved to", sky_path)
