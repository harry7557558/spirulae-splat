#!/usr/bin/python3

import torch
import numpy as np
import argparse
from io import BytesIO
import json
import os
from pathlib import Path

from typing import Literal


def quat_mult(q1, q2):
    w1, x1, y1, z1 = torch.unbind(q1, dim=-1)
    w2, x2, y2, z2 = torch.unbind(q2, dim=-1)
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    return torch.stack([w, x, y, z]).T.contiguous()

def quat_to_rotmat(quat):
    w, x, y, z = torch.unbind(quat, dim=-1)
    mat = torch.stack(
        [
            1 - 2 * (y**2 + z**2),
            2 * (x * y - w * z),
            2 * (x * z + w * y),
            2 * (x * y + w * z),
            1 - 2 * (x**2 + z**2),
            2 * (y * z - w * x),
            2 * (x * z - w * y),
            2 * (y * z + w * x),
            1 - 2 * (x**2 + y**2),
        ],
        dim=-1,
    )
    return mat.reshape(quat.shape[:-1] + (3, 3))


def generate_right_angle_rotation_matrices():
    """List all "right angle" rotation matrices"""
    import itertools
    valid_matrices = []
    for matrix in itertools.product([-1, 0, 1], repeat=9):
        R = np.array(matrix).reshape(3, 3)
        if np.linalg.det(R) == 1 and (np.dot(R.T, R) == np.eye(3)).all():
            valid_matrices.append(R)
    return valid_matrices


def generate_sh_basis(dirs: torch.Tensor, mode: Literal['nerfstudio', 'gsplat']):
    x, y, z = torch.unbind(dirs, -1)
    one = torch.ones_like(x)
    xx, yy, zz = x*x, y*y, z*z

    if mode == "nerfstudio":
        return torch.stack([
            0.28209479177387814 * one,
            0.4886025119029199 * y,
            0.4886025119029199 * z,
            0.4886025119029199 * x,
            1.0925484305920792 * x * y,
            1.0925484305920792 * y * z,
            (0.9461746957575601 * zz - 0.31539156525251999),
            1.0925484305920792 * x * z,
            0.5462742152960396 * (xx - yy),
            0.5900435899266435 * y * (3.0 * xx - yy),
            2.890611442640554 * x * y * z,
            0.4570457994644658 * y * (5.0 * zz - 1.0),
            0.3731763325901154 * z * (5.0 * zz - 3.0),
            0.4570457994644658 * x * (5.0 * zz - 1.0),
            1.445305721320277 * z * (xx - yy),
            0.5900435899266435 * x * (xx - 3.0 * yy),
            2.5033429417967046 * x * y * (xx - yy),
            1.7701307697799304 * y * z * (3.0 * xx - yy),
            0.9461746957575601 * x * y * (7.0 * zz - 1.0),
            0.6690465435572892 * y * z * (7.0 * zz - 3.0),
            0.10578554691520431 * (35.0 * zz * zz - 30.0 * zz + 3.0),
            0.6690465435572892 * x * z * (7.0 * zz - 3.0),
            0.47308734787878004 * (xx - yy) * (7.0 * zz - 1.0),
            1.7701307697799304 * x * z * (xx - 3.0 * yy),
            0.6258357354491761 * (xx * (xx - 3.0 * yy) - yy * (3.0 * xx - yy)),
        ], dim=1)

    elif mode == "gsplat":
        from spirulae_splat.splat._torch_impl import eval_sh_bases
        return eval_sh_bases(25, dirs)


def rotate_sh_coeffs(coeffs: torch.Tensor, R: torch.Tensor,
                     mode: Literal['nerfstudio', 'gsplat']) -> torch.Tensor:
    if isinstance(R, np.ndarray):
        R = torch.from_numpy(R).to(coeffs)

    # fibonacci sample on a sphere
    M = coeffs.shape[1]
    N = int(1.2*M+1)  # >= M, more is slower but more numerically stable
    idx = torch.arange(N, dtype=torch.float32, device=coeffs.device) + 0.5
    phi = 2 * np.pi * idx / ((1 + np.sqrt(5)) / 2)
    cos_theta = 1.0 - 2.0 * idx / N
    sin_theta = torch.sqrt(torch.relu(1.0 - cos_theta * cos_theta))
    dirs = torch.stack([sin_theta * torch.cos(phi),
                        sin_theta * torch.sin(phi),
                        cos_theta], dim=1)  # (N, 3)

    # evaluate basis
    Y_src = generate_sh_basis(dirs, mode)[:, 1:M+1]  # (N, M)
    dirs_rot = dirs @ R.t()  # (N, 3)
    Y_rot = generate_sh_basis(dirs_rot, mode)[:, 1:M+1]  # (N, M)
    # coeffs: (batch, M, 3)

    # direct linear transform
    f_vals = torch.einsum('mn,bnc->bmc', Y_src, coeffs)  # (batch, M, 3)
    Y_rot_inv = torch.linalg.inv(Y_rot.T @ Y_rot) @ Y_rot.T  # (N, M)
    coeffs_rot = torch.einsum('ni,bic->bnc', Y_rot_inv, f_vals)  # (batch, N, 3)

    # numerical check
    if False:
        print(torch.linalg.eigvalsh(Y_rot.T @ Y_rot))
        f_vals_rot = torch.einsum('ni,bic->bnc', Y_rot, coeffs_rot)
        print(abs(f_vals - f_vals_rot).max().item())

    return coeffs_rot


def custom_quantile(sorted_numbers, quantiles):
    n = len(sorted_numbers)
    indices = quantiles * (n - 1)
    i0 = torch.floor(indices).long()
    i1 = torch.ceil(indices).long()
    i0 = torch.clamp(i0, 0, n - 1)
    i1 = torch.clamp(i1, 0, n - 1)
    v0 = sorted_numbers[i0]
    v1 = sorted_numbers[i1]
    t = indices - i0.float()
    return v0 * (1 - t) + v1 * t


@torch.no_grad()
def quantize_tensor(name, tensor, n_bits, maxiter=0):
    device = tensor.device
    if tensor.numel() == 0:
        print(name, "empty tensor")
        centroids = torch.arange(0).float().to(device)
        bins = torch.zeros_like(tensor).long()
        return centroids, bins

    sorted_numbers, sorted_indices = torch.sort(tensor.flatten())
    quantiles = ((torch.arange(2**n_bits)+0.5) / (2**n_bits)).float().to(device)
    try:
        raise RuntimeError
        centroids = torch.quantile(sorted_numbers, quantiles)
    except RuntimeError:
        centroids = custom_quantile(sorted_numbers, quantiles)
    # initial_centroids = centroids.clone()
    
    initial_error = None
    prev_error = None

    for _ in range(maxiter):
        # divide into bins
        centers = (centroids[1:]+centroids[:-1])/2
        bins = torch.bucketize(sorted_numbers, centers)

        # monitor L1, which may increase as MSE decreases
        error = torch.mean(torch.abs(centroids[bins]-sorted_numbers)).item()
        # print(_, error)

        # Lloyd-Max quantization
        new_centroids = centroids.clone()
        for i in range(2**n_bits):
            m = sorted_numbers[bins == i]
            if len(m) > 0:
                new_centroids[i] += 1.0 * (m.mean()-new_centroids[i])
        centroids, indices = torch.sort(new_centroids)

        # convergence check
        if initial_error is None:
            initial_error = error
        if prev_error is not None and _ > 5 and \
            prev_error-error < 1e-3 * initial_error:
            break
        prev_error = error

    centers = (centroids[1:]+centroids[:-1])/2
    bins = torch.bucketize(tensor.contiguous(), centers)
    print(name, 'quantization loss:', torch.mean(torch.abs(centroids[bins]-tensor)).item())
    return centroids, bins

    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.kdeplot(data=tensor.cpu().numpy().flatten(), bw_adjust=0.2)
    # plt.plot(initial_centroids, 0.0*initial_centroids, '.')
    plt.plot(centroids.cpu().numpy(), 0.0*centroids.cpu().numpy(), '.')
    plt.show()


def pack_components_(config, components):
    components = components[:]
    for i, comp in enumerate(components):
        if isinstance(comp, torch.Tensor):
            components[i] = comp.cpu().numpy()

    length = config["length"]
    component_length = config["componentLength"]
    component_views = config["componentViews"]
    
    packed_data = np.zeros((length * component_length,), dtype=np.uint8)
    
    for i in range(length):
        for j, view in enumerate(component_views):
            data = components[j][i]
            bit_length = view["bitLength"]
            bit_offset = view["bitOffset"]
            data_type = view["type"]
            
            if data_type.startswith('quat'):
                num_elements = int(data_type[4:]) if len(data_type) > 4 else 1
                if num_elements == 0:
                    continue
                assert bit_length % num_elements == 0
                bits_per_element = bit_length // num_elements
                for k in range(num_elements):
                    value = data[k]
                    assert (value >= 0) and value < (1 << bits_per_element)
                    bit_start = bit_offset + k * bits_per_element
                    byte_start = bit_start // 8
                    bit_pos = bit_start % 8
                    for b in range(bits_per_element):
                        byte_idx = byte_start + (bit_pos + b) // 8
                        bit_in_byte = (bit_pos + b) % 8
                        if value & (1 << b):
                            packed_data[i*component_length+byte_idx] |= (1 << bit_in_byte)
            
            elif data_type.startswith('float'):
                num_elements = int(data_type[5:]) if len(data_type) > 5 else 1
                data = data.astype(np.float32).tobytes()
                assert bit_offset % 8 == 0
                for k in range(num_elements):
                    byte_start = (bit_offset // 8) + k * 4
                    i0 = i*component_length+byte_start
                    packed_data[i0:i0+4] = np.frombuffer(data[k*4:(k+1)*4], dtype=np.byte)
            
            elif data_type.startswith('byte'):
                num_elements = int(data_type[4:]) if len(data_type) > 4 else 1
                data = data.clip(0, 255).astype(np.uint8).tobytes()
                assert bit_offset % 8 == 0
                for k in range(num_elements):
                    byte_start = (bit_offset // 8) + k
                    packed_data[i*component_length+byte_start] = data[k]
    
    return packed_data.tobytes()


def pack_components(config, components):
    import subprocess
    import pybind11
    import importlib.util

    binary_path = "./pack_components.so"
    if not os.path.exists(binary_path):
        python_include = subprocess.check_output(["pkg-config", "--cflags", "python3"]).decode("utf-8").strip()
        compile_command = [
            "c++", "-O3", "-Wall", "-shared", "-std=c++11", "-fPIC", "pack_components.cpp", "-o", binary_path,
            "-I", pybind11.get_include(), *(python_include.split())
        ]
        subprocess.run(compile_command, check=True)

    spec = importlib.util.spec_from_file_location("pack_components", binary_path)
    pack_components_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pack_components_module)

    components = components[:]
    for i, comp in enumerate(components):
        if isinstance(comp, torch.Tensor):
            components[i] = comp.cpu().numpy().astype(np.int32)

    config_json = json.loads(json.dumps(config))
    packed_data = pack_components_module.pack_components(config_json, components)
    # print(packed_data)
    
    return bytearray(packed_data)


def component_view_psa(component_views):
    total = 0
    for view in component_views:
        view['bitOffset'] = total
        total += view['bitLength']
    return component_views, int(np.ceil(total/8))

def bufferview_psa(buffer_views):
    total = 0
    for view in buffer_views:
        view['byteOffset'] = total
        total += view['byteLength']
    return buffer_views


def load_nerfstudio_data(config):
    from nerfstudio.utils.io import load_from_json
    from nerfstudio.cameras import camera_utils
    from nerfstudio.data.utils.dataparsers_utils import get_train_eval_split_fraction

    data_dir = config.data
    meta = load_from_json(data_dir / "transforms.json")

    image_filenames = []
    poses = []

    fx_fixed = "fl_x" in meta
    fy_fixed = "fl_y" in meta
    cx_fixed = "cx" in meta
    cy_fixed = "cy" in meta
    height_fixed = "h" in meta
    width_fixed = "w" in meta
    distort_fixed = False
    for distort_key in ["k1", "k2", "k3", "p1", "p2", "distortion_params"]:
        if distort_key in meta:
            distort_fixed = True
            break
    fx = []
    fy = []
    cx = []
    cy = []
    height = []
    width = []
    distort = []

    # sort the frames by fname
    fnames = []
    for frame in meta["frames"]:
        filepath = Path(frame["file_path"])
        fname = data_dir / filepath
        fnames.append(fname)
    inds = np.argsort(fnames)
    frames = [meta["frames"][ind] for ind in inds]

    for frame in frames:
        filepath = Path(frame["file_path"])
        fname = data_dir / filepath

        if not fx_fixed:
            assert "fl_x" in frame, "fx not specified in frame"
            fx.append(float(frame["fl_x"]))
        if not fy_fixed:
            assert "fl_y" in frame, "fy not specified in frame"
            fy.append(float(frame["fl_y"]))
        if not cx_fixed:
            assert "cx" in frame, "cx not specified in frame"
            cx.append(float(frame["cx"]))
        if not cy_fixed:
            assert "cy" in frame, "cy not specified in frame"
            cy.append(float(frame["cy"]))
        if not height_fixed:
            assert "h" in frame, "height not specified in frame"
            height.append(int(frame["h"]))
        if not width_fixed:
            assert "w" in frame, "width not specified in frame"
            width.append(int(frame["w"]))
        if not distort_fixed:
            distort.append(
                np.array(frame["distortion_params"], dtype=np.float32)
                if "distortion_params" in frame
                else camera_utils.get_distortion_params(
                    k1=float(frame["k1"]) if "k1" in frame else 0.0,
                    k2=float(frame["k2"]) if "k2" in frame else 0.0,
                    k3=float(frame["k3"]) if "k3" in frame else 0.0,
                    k4=float(frame["k4"]) if "k4" in frame else 0.0,
                    p1=float(frame["p1"]) if "p1" in frame else 0.0,
                    p2=float(frame["p2"]) if "p2" in frame else 0.0,
                )
            )

        image_filenames.append(fname)
        poses.append(np.array(frame["transform_matrix"], dtype=np.float32))

    indices, _ = get_train_eval_split_fraction(image_filenames, config.train_split_fraction)

    poses = torch.from_numpy(np.array(poses).astype(np.float32))
    poses, transform_matrix = camera_utils.auto_orient_and_center_poses(
        poses,
        method=config.orientation_method,
        center_method=config.center_method,
    )
    poses = poses.numpy()

    # Scale poses
    scale_factor = 1.0
    if config.auto_scale_poses:
        scale_factor /= float(np.max(np.abs(poses[:, :3, 3])))
    scale_factor *= config.scale_factor

    poses[:, :3, 3] *= scale_factor

    idx_tensor = np.array(indices, dtype=np.int32)
    poses = poses[idx_tensor]

    camera_type = meta["camera_model"] if "camera_model" in meta else 'PERSPECTIVE'

    fx = float(meta["fl_x"]) if fx_fixed else np.array(fx, dtype=np.float32)[idx_tensor]
    fy = float(meta["fl_y"]) if fy_fixed else np.array(fy, dtype=np.float32)[idx_tensor]
    cx = float(meta["cx"]) if cx_fixed else np.array(cx, dtype=np.float32)[idx_tensor]
    cy = float(meta["cy"]) if cy_fixed else np.array(cy, dtype=np.float32)[idx_tensor]
    height = int(meta["h"]) if height_fixed else np.array(height, dtype=np.int32)[idx_tensor]
    width = int(meta["w"]) if width_fixed else np.array(width, dtype=np.int32)[idx_tensor]
    if distort_fixed:
        distortion_params = (
            np.array(meta["distortion_params"], dtype=np.float32)
            if "distortion_params" in meta
            else camera_utils.get_distortion_params(
                k1=float(meta["k1"]) if "k1" in meta else 0.0,
                k2=float(meta["k2"]) if "k2" in meta else 0.0,
                k3=float(meta["k3"]) if "k3" in meta else 0.0,
                k4=float(meta["k4"]) if "k4" in meta else 0.0,
                p1=float(meta["p1"]) if "p1" in meta else 0.0,
                p2=float(meta["p2"]) if "p2" in meta else 0.0,
            )
        )
    else:
        distortion_params = np.stack(distort, axis=0)[idx_tensor]

    return {
        'fx': fx,
        'fy': fy,
        'cx': cx,
        'cy': cy,
        'dist': distortion_params,
        'h': height,
        'w': width,
        'c2w': poses[:, :3, :4],
        'type': camera_type
    }


class SplatModel:
    def __init__(self, file_path):
        if file_path.endswith('.ckpt'):
            self.load_ckpt(file_path)
        elif file_path.endswith('config.yml'):
            self.load_config(file_path)
        elif file_path.endswith('.ply'):
            self.load_ply(file_path)

    def load_ckpt(self, file_path):
        checkpoint = torch.load(file_path)
        pipeline = checkpoint['pipeline']
        # print(pipeline.keys())
        # print(pipeline['_model.camera_optimizer.pose_adjustment'])

        self.features_dc = pipeline['_model.gauss_params.features_dc']  # 8 bits
        self.features_sh = pipeline['_model.gauss_params.features_sh']  # 4 bits
        self.features_ch = pipeline['_model.gauss_params.features_ch']  # 6 bits
        self.means = pipeline['_model.gauss_params.means']  # 12-16 bits
        self.opacities = pipeline['_model.gauss_params.opacities']  # 8 bits after sigmoid
        self.quats = pipeline['_model.gauss_params.quats']  # 8 bits
        self.scales = pipeline['_model.gauss_params.scales']  # 8-12 bits
        self.background_color = pipeline['_model.background_color']
        self.background_sh = pipeline['_model.background_sh'] \
            if '_model.background_sh' in pipeline else torch.zeros((0, 3))

    def load_config(self, file_path: str):
        save_dir = file_path[:file_path.rfind(os.path.sep)]
        ckpt_dir = os.path.join(save_dir, 'nerfstudio_models')
        for f in os.listdir(ckpt_dir):
            if f.endswith('.ckpt'):
                f = os.path.join(ckpt_dir, f)
                self.load_ckpt(f)
                break
        return

        # attempt to find camera poses
        import yaml
        with open(file_path, 'r') as fp:
            config = yaml.load(fp, yaml.UnsafeLoader)
        # print(config)
        launch_dir = file_path[:file_path.rfind(str(config.output_dir))]
        data_dir = os.path.join(launch_dir, str(config.data))
        dataparser = config.pipeline.datamanager.dataparser
        dataparser.data = Path(data_dir)
        print(dataparser)
        cameras = load_nerfstudio_data(dataparser)
        print(cameras)

    def load_ply(self, file_path):
        """for official 2DGS codebase"""
        from plyfile import PlyData
        plydata = PlyData.read(file_path)

        means = []
        scales = []
        quats = []
        features_dc = []
        opacities = []
        features_sh = []

        v = plydata["vertex"][0]
        sh_keys = []
        for i in range(72):
            key = f"f_rest_{i}"
            try:
                v[key]
                sh_keys.append(key)
            except ValueError:
                break

        for v in plydata["vertex"]:
            means.append([v["x"], v["y"], v["z"]])
            scales.append([v["scale_0"], v["scale_1"]])
            quats.append([v["rot_0"], v["rot_1"], v["rot_2"], v["rot_3"]])
            features_dc.append([v["f_dc_0"], v["f_dc_1"], v["f_dc_2"]])
            opacities.append([v["opacity"]])
            features_sh.append([v[key] for key in sh_keys])

        n = len(means)
        self.means = torch.tensor(means, dtype=torch.float)
        self.opacities = torch.tensor(opacities, dtype=torch.float)
        self.quats = torch.tensor(quats, dtype=torch.float)
        self.scales = torch.tensor(scales, dtype=torch.float)
        self.features_dc = torch.tensor(features_dc, dtype=torch.float)
        self.features_sh = torch.tensor(features_sh, dtype=torch.float).view(n, -1, 3)
        self.features_ch = torch.zeros((n, 0, 3))
        self.background_color = torch.zeros(3)
        self.background_sh = torch.zeros((0, 3))


def process_ckpt_to_ssplat(file_path, meta={}, cull_th=float('inf'), cull_n=1, exposure=1.0, scale=1.0, rotate=False):

    m = SplatModel(file_path)
    (features_dc, features_sh, features_ch,
        means, opacities, quats, scales,
        background_color, background_sh) = \
    (m.features_dc, m.features_sh, m.features_ch,
        m.means, m.opacities, m.quats, m.scales,
        m.background_color, m.background_sh)

    num_sh = features_sh.shape[1]
    num_ch = features_ch.shape[1]
    sh_degree = {
        1: 0, 4: 1, 9: 2, 16: 3, 25: 4
    }[num_sh+1]
    ch_degree_r, ch_degree_phi = {
        0: (0, 0), 1: (1, 0), 2: (2, 0),
        3: (1, 1), 5: (1, 2), 6: (2, 1),
        7: (1, 3), 9: (3, 1), 10: (2, 2),
        14: (2, 3), 15: (3, 2), 21: (3, 3)
    }[num_ch]
    background_sh_degree = {
        1: 0, 4: 1, 9: 2, 16: 3, 25: 4
    }[len(background_sh)+1]
    print("SH degree:", sh_degree)
    print("# of SH:", num_sh)
    print("CH degree r:", ch_degree_r)
    print("CH degree phi:", ch_degree_phi)
    print("# of CH:", num_ch)
    print("background:", background_color.cpu().numpy())
    print("background SH degree:", background_sh_degree)
    print()

    n_splat = len(means)
    n_floats = sum([
        x.numel() for x in [
            features_dc, features_sh, features_ch,
            means, opacities, quats, scales]
    ])
    print(n_splat, "splats")
    print(n_floats//n_splat, "floats per splat")
    print("{:.2f} MB floats".format(n_floats*4/1024**2))
    print()

    if 0 < cull_th < float('inf'):
        for _ in range(cull_n):
            center = torch.mean(means, axis=0)
            cov_matrix = torch.cov(means.T)
            mahalanobis_dist = torch.sqrt(torch.sum((
                (means-center) @ torch.linalg.inv(cov_matrix)) * (means-center), axis=1))
            mask = mahalanobis_dist < cull_th
            print("culled", torch.sum(~mask).item(), '/', len(means), "splats")

            (features_dc, features_sh, features_ch,
            means, opacities, quats, scales) = [
                tensor[mask].contiguous() for tensor in (
                    features_dc, features_sh, features_ch,
                    means, opacities, quats, scales)
            ]
        print()

    # normalize position
    means -= torch.mean(means, axis=0)
    # normalize scale
    scale = scale / ((means*means).sum(-1).mean()**0.5).item()
    means = means * scale
    scales = scales + np.log(scale)

    quats = quats / torch.norm(quats, dim=1, keepdim=True)
    opacities = torch.sigmoid(opacities)

    weight = torch.exp(scales[:,0]+scales[:,1]) * opacities[:,0]

    if rotate:
        from scipy.spatial.transform import Rotation
        # find rotation matrix
        print("Find rotation matrix...")
        if False:
            s = torch.cov(means.T, aweights=weight).cpu().numpy().astype(np.float64)
            s, r0 = np.linalg.eigh(s)
        else:
            normals = quat_to_rotmat(quats)[:, 2]
            m = torch.einsum('ki,kj->kij', normals, normals) - torch.eye(3).unsqueeze(0).to(normals)
            m = (m * weight.reshape(-1, 1, 1)).sum(0).cpu().numpy().astype(np.float64)
            s, r0 = np.linalg.eigh(m)
        if np.linalg.det(r0) < 0.0:
            r0[:, 2] *= -1.0
        possible_mats = generate_right_angle_rotation_matrices()
        diffs = [np.abs(r0-r).sum() for r in possible_mats]
        r1 = possible_mats[np.argmin(diffs)]
        dr = r1 @ r0.T

        # rotate position
        print("Rotate position and quaternion...")
        means = means @ torch.from_numpy(dr.T).to(means)

        # rotate quaternion
        dq = Rotation.from_matrix(dr).as_quat()
        dq = torch.from_numpy(dq[[3,0,1,2]]).to(quats)
        quats = quat_mult(dq, quats)

        # rotate SH coefficients
        print("Rotate SH...")
        features_sh = rotate_sh_coeffs(features_sh, dr, "gsplat")
        if background_sh_degree > 0:
            background_sh = rotate_sh_coeffs(background_sh[None], dr, "nerfstudio")[0]
        print()

    if exposure != 1.0:
        exposure = max(exposure, 0.0) ** (1.0/2.2)
        # adjust albedo
        color = features_dc*0.2820947917738781+0.5 if sh_degree > 0 else features_dc
        # color0 = torch.clip(color*exposure, max=1.0)
        color = 1.0 - torch.relu(1.0-color) ** exposure if exposure > 1.0 else color * exposure
        features_dc = (color - 0.5) / 0.2820947917738781 if sh_degree > 0 else color
        # adjust SH
        features_sh *= exposure * (1.0-color.unsqueeze(1)) if exposure > 1.0 else exposure
        # adjust background
        if background_sh_degree > 0:
            print(background_color)
            color = background_color*0.2820947917738781+0.5
            color = 1.0 - torch.relu(1.0-color) ** exposure if exposure > 1.0 else color * exposure
            background_color = (color - 0.5) / 0.2820947917738781
            # background_sh *= exposure * (1.0-color) if exposure > 1.0 else exposure
            background_sh *= exposure
        else:
            background_color = 1.0 - torch.relu(1.0-background_color) ** exposure if exposure > 1.0 else background_color * exposure
        # TODO: adjust CH


    sorted_indices = torch.argsort(-weight)
    sorted_indices = torch.arange(len(weight))

    from time import perf_counter
    time0 = perf_counter()
    with torch.no_grad():
        means_bins, means_q = quantize_tensor('means', means, 12)
        scales_bins, scales_q = quantize_tensor('scale', scales, 10)
        quats_bins, quats_q = quantize_tensor('quats', quats, 8)
        opacities_bins, opacities_q = quantize_tensor('opacs', opacities, 6)
        features_dc_bins, features_dc_q = quantize_tensor('color', features_dc, 6)
        features_ch_bins, features_ch_q = quantize_tensor('ch', features_ch, bit_ch, 100)
        features_sh_bins, features_sh_q = quantize_tensor('sh', features_sh, bit_sh, 100)
    time1 = perf_counter()
    print()

    buffer_views = [
        { "byteLength": 4*len(means_bins) },
        { "byteLength": 4*len(scales_bins) },
        { "byteLength": 4*len(quats_bins) },
        { "byteLength": 4*len(opacities_bins) },
        { "byteLength": 4*len(features_dc_bins) },
        { "byteLength": 4*len(features_ch_bins) },
        { "byteLength": 4*len(features_sh_bins) },
    ]

    print("Packing base...")
    componentViews, componentLength = component_view_psa([
        { "key": "means", "type": "quat3", "bitLength": 36, "quatBufferView": 0 },
        { "key": "scales", "type": "quat2", "bitLength": 20, "quatBufferView": 1 },
        { "key": "quats", "type": "quat4", "bitLength": 32, "quatBufferView": 2 },
        { "key": "opacities", "type": "quat", "bitLength": 6, "quatBufferView": 3 },
        { "key": "features_dc", "type": "quat3", "bitLength": 18, "quatBufferView": 4 },
    ])
    base_config = {
        "bufferView": len(buffer_views),
        "length": len(weight),
        "componentLength": componentLength,  # in bytes
        "componentViews": componentViews
    }
    base_buffer = pack_components(
        base_config,
        [means_q, scales_q, quats_q, opacities_q, features_dc_q]
    )

    print("Packing harmonics...")
    componentViews, componentLength = component_view_psa([
            { "key": "features_ch", "type": f"quat{3*num_ch}", "bitLength": 3*num_ch*bit_ch, "quatBufferView": 5 },
    ] * (num_ch > 0) + [
            { "key": "features_sh", "type": f"quat{3*num_sh}", "bitLength": 3*num_sh*bit_sh, "quatBufferView": 6 },
    ] * (num_sh > 0))
    harmonics_config = {
        "bufferView": len(buffer_views)+1,
        "length": len(weight),
        "componentLength": componentLength,  # in bytes
        "componentViews": componentViews
    }
    harmonics_buffer = pack_components(harmonics_config, [
        features_ch_q.reshape((len(features_ch), -1))
    ] * (num_ch > 0) + [
        features_sh_q.reshape((len(features_sh), -1))
    ] * (num_sh > 0))

    buffer = BytesIO()
    buffer.write(means_bins.cpu().numpy().tobytes())
    buffer.write(scales_bins.cpu().numpy().tobytes())
    buffer.write(quats_bins.cpu().numpy().tobytes())
    buffer.write(opacities_bins.cpu().numpy().tobytes())
    buffer.write(features_dc_bins.cpu().numpy().tobytes())
    buffer.write(features_ch_bins.cpu().numpy().tobytes())
    buffer.write(features_sh_bins.cpu().numpy().tobytes())
    buffer.write(base_buffer)
    buffer.write(harmonics_buffer)
    buffer = buffer.getvalue()
    while len(buffer) % 4:
        buffer += b'\0'

    header = {
        "asset": {
            "version": "0.0",
            "meta": meta,
        },
        "config": {
            "sh_degree": sh_degree,
            "ch_degree_r": ch_degree_r,
            "ch_degree_phi": ch_degree_phi,
            "background_color": background_color.cpu().numpy().flatten().tolist(),
            "background_sh_degree": background_sh_degree,
            "background_sh": background_sh.cpu().numpy().flatten().tolist()
        },
        "buffer": {
            "byteLength": len(buffer),
        },
        "primitives": {
            "base": base_config,
            "harmonics": harmonics_config,
        },
        "bufferViews": bufferview_psa(
            buffer_views + \
            [ { "byteLength": len(base_buffer) } ] + \
            [ { "byteLength": len(harmonics_buffer) } ] * (len(harmonics_buffer) > 0)
        )
    }

    header = json.dumps(header)
    while len(header) % 4:
        header += " "

    ssplat = BytesIO()
    ssplat.write(b"splt")
    ssplat.write(len(header).to_bytes(4, 'little'))
    ssplat.write(bytearray(header, 'utf-8'))
    ssplat.write(buffer)
    ssplat = ssplat.getvalue()
    return ssplat


def save_ssplat_file(ssplat_data, output_path):
    with open(output_path, "wb") as f:
        f.write(ssplat_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Nerfstudio spirulae-splat checkpoint files to SSPLAT format.")
    parser.add_argument(
        "input_files", nargs="+", help="The input Nerfstudio spirulae-splat checkpoint files to process."
    )
    parser.add_argument(
        "--output", "-o", default="model.ssplat", help="The output SSPLAT file."
    )
    parser.add_argument("--cull", "-c", default=float("inf"), type=float, help="Cull threshold.")
    parser.add_argument("--ncull", "-cn", default=1, type=int, help="Number of cull iterations.")
    parser.add_argument("--exposure", "-e", default=1.0, type=float, help="Relative exposure.")
    parser.add_argument("--scale", "-s", default=1.0, type=float, help="Relative scale.")
    parser.add_argument("--rotate", "-r", default=False, type=bool, help="Whether to rotate the scene to align with axis.")
    parser.add_argument("--bit_sh", "-bsh", default=3, type=int, help="Bits for each SH coefficient.")
    parser.add_argument("--bit_ch", "-bch", default=4, type=int, help="Bits for each CH coefficient.")
    args = parser.parse_args()

    bit_sh = args.bit_sh
    bit_ch = args.bit_ch

    meta = {
        'argv': __import__('sys').argv,
        'timestamp': __import__('datetime').datetime.now().astimezone().isoformat(),
    }

    cull_th = args.cull
    cull_n = args.ncull
    exposure = args.exposure
    rotate = args.rotate
    scale = args.scale

    for input_file in args.input_files:
        print(f"Processing {input_file}...", end='\n\n')
        ssplat_data = process_ckpt_to_ssplat(input_file, meta, cull_th, cull_n, exposure, scale, rotate)
        output_file = (
            args.output if len(args.input_files) == 1 else input_file + ".ssplat"
        )
        save_ssplat_file(ssplat_data, output_file)
        file_size = len(ssplat_data) / 2**20
        print(f"Saved {output_file} - {file_size:.2f} MB")
