#!/usr/bin/python

import torch
import numpy as np
import argparse
from io import BytesIO
import json


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
def quantize_tensor(tensor, n_bits, maxiter=0):
    device = tensor.device
    if tensor.numel() == 0:
        print("empty tensor")
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
    bins = torch.bucketize(tensor, centers)
    print('quatization loss:', torch.mean(torch.abs(centroids[bins]-tensor)).item())
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
    import os
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


class SplatModel:
    def __init__(self, file_path):
        if file_path.endswith('.ckpt'):
            self.load_ckpt(file_path)
        elif file_path.endswith('.ply'):
            self.load_ply(file_path)

    def load_ckpt(self, file_path):
        checkpoint = torch.load(file_path)
        pipeline = checkpoint['pipeline']
        # print(pipeline.keys())

        self.features_dc = pipeline['_model.gauss_params.features_dc']  # 8 bits
        self.features_sh = pipeline['_model.gauss_params.features_sh']  # 4 bits
        self.features_ch = pipeline['_model.gauss_params.features_ch']  # 6 bits
        self.means = pipeline['_model.gauss_params.means']  # 12-16 bits
        self.opacities = pipeline['_model.gauss_params.opacities']  # 8 bits after sigmoid
        self.anisotropies = pipeline['_model.gauss_params.anisotropies']  # 8 bits
        self.quats = pipeline['_model.gauss_params.quats']  # 8 bits
        self.scales = pipeline['_model.gauss_params.scales']  # 8-12 bits
        self.background_color = pipeline['_model.background_color']
        self.background_sh = pipeline['_model.background_sh'] \
            if '_model.background_sh' in pipeline else torch.zeros((0, 3))

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
        self.anisotropies = torch.zeros((n, 2))
        self.background_color = torch.zeros(3)
        self.background_sh = torch.zeros((0, 3))


def process_ckpt_to_ssplat(file_path):

    m = SplatModel(file_path)
    (features_dc, features_sh, features_ch,
        means, opacities, anisotropies, quats, scales,
        background_color, background_sh) = \
    (m.features_dc, m.features_sh, m.features_ch,
        m.means, m.opacities, m.anisotropies, m.quats, m.scales,
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
    use_anisotropy = (anisotropies != 0.0).any().item()
    print("SH degree:", sh_degree)
    print("# of SH:", num_sh)
    print("CH degree r:", ch_degree_r)
    print("CH degree phi:", ch_degree_phi)
    print("# of CH:", num_ch)
    print("background:", background_color.cpu().numpy())
    print("background SH degree:", background_sh_degree)
    print("anisotropy:", use_anisotropy)
    print()

    n_splat = len(means)
    n_floats = sum([
        x.numel() for x in [
            features_dc, features_sh, features_ch,
            means, opacities, anisotropies, quats, scales]
    ])
    print(n_splat, "splats")
    print(n_floats//n_splat, "floats per splat")
    print("{:.2f} MB floats".format(n_floats*4/1024**2))
    print()

    if False:
        center = torch.mean(means, axis=0)
        cov_matrix = torch.cov(means.T)
        mahalanobis_dist = torch.sqrt(torch.sum((
            (means-center) @ torch.linalg.inv(cov_matrix)) * (means-center), axis=1))
        mask = mahalanobis_dist < 1.5
        print("masked", torch.sum(mask).item(), "splats")
        print()

        (features_dc, features_sh, features_ch,
        means, opacities, anisotropies, quats, scales) = [
            tensor[mask].contiguous() for tensor in (
                features_dc, features_sh, features_ch,
                means, opacities, anisotropies, quats, scales)
        ]

    means -= torch.mean(means, axis=0)

    quats = quats / torch.norm(quats, dim=1, keepdim=True)
    opacities = torch.sigmoid(opacities)
    anisotropies = torch.arcsinh(anisotropies)

    weight = torch.exp(scales[:,0]+scales[:,1]) * opacities[:,0]
    sorted_indices = torch.argsort(-weight)
    sorted_indices = torch.arange(len(weight))

    from time import perf_counter
    time0 = perf_counter()
    with torch.no_grad():
        means_bins, means_q = quantize_tensor(means, 12)
        scales_bins, scales_q = quantize_tensor(scales, 10)
        quats_bins, quats_q = quantize_tensor(quats, 8)
        anisotropies_bins, anisotropies_q = quantize_tensor(anisotropies, 8*use_anisotropy)
        opacities_bins, opacities_q = quantize_tensor(opacities, 6)
        features_dc_bins, features_dc_q = quantize_tensor(features_dc, 6)
        features_ch_bins, features_ch_q = quantize_tensor(features_ch, 4, 100)
        features_sh_bins, features_sh_q = quantize_tensor(features_sh, 3, 100)
    time1 = perf_counter()
    print()

    buffer_views = [
        { "byteLength": 4*len(means_bins) },
        { "byteLength": 4*len(scales_bins) },
        { "byteLength": 4*len(quats_bins) },
        { "byteLength": 4*len(anisotropies_bins) },
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
    ] + [
        { "key": "anisotropies", "type": "quat2", "bitLength": 16, "quatBufferView": 3 },
    ] * use_anisotropy + [
        { "key": "opacities", "type": "quat", "bitLength": 6, "quatBufferView": 4 },
        { "key": "features_dc", "type": "quat3", "bitLength": 18, "quatBufferView": 5 },
    ])
    base_config = {
        "bufferView": len(buffer_views),
        "length": len(weight),
        "componentLength": componentLength,  # in bytes
        "componentViews": componentViews
    }
    base_buffer = pack_components(
        base_config,
        [means_q, scales_q, quats_q] + \
            [anisotropies_q] * use_anisotropy + \
            [opacities_q, features_dc_q]
    )

    print("Packing harmonics...")
    componentViews, componentLength = component_view_psa([
            { "key": "features_ch", "type": f"quat{3*num_ch}", "bitLength": 3*num_ch*4, "quatBufferView": 6 },
    ] * (num_ch > 0) + [
            { "key": "features_sh", "type": f"quat{3*num_sh}", "bitLength": 3*num_sh*3, "quatBufferView": 7 },
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
    buffer.write(anisotropies_bins.cpu().numpy().tobytes())
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
            "version": "0.0"
        },
        "config": {
            "sh_degree": sh_degree,
            "ch_degree_r": ch_degree_r,
            "ch_degree_phi": ch_degree_phi,
            "background_color": background_color.cpu().numpy().flatten().tolist(),
            "background_sh_degree": background_sh_degree,
            "background_sh": background_sh.cpu().numpy().flatten().tolist(),
            "use_anisotropy": use_anisotropy
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


def main():
    parser = argparse.ArgumentParser(
        description="Convert Nerfstudio spirulae-splat checkpoint files to SSPLAT format.")
    parser.add_argument(
        "input_files", nargs="+", help="The input Nerfstudio spirulae-splat checkpoint files to process."
    )
    parser.add_argument(
        "--output", "-o", default="model.ssplat", help="The output SSPLAT file."
    )
    args = parser.parse_args()
    for input_file in args.input_files:
        print(f"Processing {input_file}...", end='\n\n')
        ssplat_data = process_ckpt_to_ssplat(input_file)
        output_file = (
            args.output if len(args.input_files) == 1 else input_file + ".ssplat"
        )
        save_ssplat_file(ssplat_data, output_file)
        print(f"Saved {output_file}")


if __name__ == "__main__":
    main()
