#!/usr/bin/python

import torch
import numpy as np
import argparse
from io import BytesIO
import json


def quantize_tensor(tensor, n_bits, maxiter=100):
    device = tensor.device

    sorted_numbers, sorted_indices = torch.sort(tensor.flatten())

    centroids = torch.quantile(
        sorted_numbers,
        ((torch.arange(2**n_bits)+0.5) / (2**n_bits)).float().to(device)
    )
    # initial_centroids = centroids.clone()
    
    initial_error = None
    prev_error = None

    for _ in range(maxiter):
        # divide into bins
        centers = (centroids[1:]+centroids[:-1])/2
        bins = torch.bucketize(sorted_numbers, centers)

        # monitor L1, which may increase as MSE decreases
        error = torch.mean(torch.abs(centroids[bins]-sorted_numbers)).item()
        print(_, error)

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
    print('error:', torch.mean(torch.abs(centroids[bins]-tensor)).item())
    return centroids, bins

    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.kdeplot(data=tensor.cpu().numpy().flatten(), bw_adjust=0.2)
    # plt.plot(initial_centroids, 0.0*initial_centroids, '.')
    plt.plot(centroids.cpu().numpy(), 0.0*centroids.cpu().numpy(), '.')
    plt.show()


def pack_components(config, components):
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


def process_ckpt_to_splat(file_path):
    checkpoint = torch.load(file_path)
    pipeline = checkpoint['pipeline']
    
    # print(pipeline.keys())
    features_dc = pipeline['_model.gauss_params.features_dc']  # 8 bits
    features_sh = pipeline['_model.gauss_params.features_sh']  # 4 bits
    features_ch = pipeline['_model.gauss_params.features_ch']  # 6 bits
    means = pipeline['_model.gauss_params.means']  # 12-16 bits
    opacities = pipeline['_model.gauss_params.opacities']  # 8 bits after sigmoid
    anisotropies = pipeline['_model.gauss_params.anisotropies']  # 8 bits
    quats = pipeline['_model.gauss_params.quats']  # 8 bits
    scales = pipeline['_model.gauss_params.scales']  # 8-12 bits
    background_color = pipeline['_model.background_color']

    num_sh = features_sh.shape[1]
    num_ch = features_ch.shape[1]
    sh_degree = {
        1: 0, 4: 1, 9: 2, 16: 3, 25: 4
    }[num_sh+1]
    print("SH degree:", sh_degree)
    print("# of CH:", num_ch)
    print("background:", background_color.cpu().numpy())
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

    weight = torch.exp(scales[:,0]+scales[:,1]) * opacities[:,0]
    sorted_indices = torch.argsort(-weight)
    sorted_indices = torch.arange(len(weight))

    from time import perf_counter
    time0 = perf_counter()
    with torch.no_grad():
        means_bins, means_q = quantize_tensor(means, 12, maxiter=0)
        scales_bins, scales_q = quantize_tensor(scales, 10, maxiter=0)
    time1 = perf_counter()
    print("time elapsed:", time1-time0)

    SH_C0 = 0.28209479177387814
    if sh_degree > 0:
        features_dc = 0.5 + features_dc * SH_C0
    if num_ch > 0:  # has CH and not use CH
        features_dc *= 0.5
    opacities = torch.sigmoid(opacities)

    means_q = means_q.cpu().numpy()
    means_bins = means_bins.cpu().numpy()
    scales_q = scales_q.cpu().numpy()
    scales_bins = scales_bins.cpu().numpy()
    quats = quats.cpu().numpy()
    features_dc = features_dc.cpu().numpy()
    opacities = opacities.cpu().numpy()

    base_config = {
        "bufferView": 2,
        "length": len(weight),
        "componentLength": int(np.ceil((36+20+24+8+32)/8)),  # in bytes
        "componentViews": [
            { "key": "mean", "type": "quat3", "bitLength": 36, "bitOffset": 0, "quatBufferView": 0 },
            { "key": "scale", "type": "quat2", "bitLength": 20, "bitOffset": 36, "quatBufferView": 1 },
            { "key": "feature_dc", "type": "byte3", "bitLength": 24, "bitOffset": 36+20 },
            { "key": "opacity", "type": "byte", "bitLength": 8, "bitOffset": 36+20+24 },
            { "key": "quat", "type": "byte4", "bitLength": 32, "bitOffset": 36+20+24+8 },
        ]
    }
    base_buffer = pack_components(base_config, [
        means_q, scales_q,
        (features_dc*255+0.5),
        (opacities*255+0.5),
        ((quats / np.linalg.norm(quats, axis=1).reshape(-1,1)) * 128 + 128.5)])

    header = {
        "asset": {
            "version": "0.0"
        },
        "config": {
            "sh_degree": sh_degree,
            "ch_degree_r": 3,
            "ch_degree_phi": 3,
            "background_color": background_color.cpu().numpy().tolist()
        },
        "buffer": {
            "byteLength": len(base_buffer),
        },
        "primitives": {
            "base": base_config
        },
        "bufferViews": [
            {
                "byteLength": 4*len(means_bins),
                "byteOffset": 0
            },
            {
                "byteLength": 4*len(scales_bins),
                "byteOffset": 4*len(means_bins)
            },
            {
                "byteLength": len(base_buffer),
                "byteOffset": 4*len(means_bins)+4*len(scales_bins)
            }
        ]
    }
    header = json.dumps(header)
    while len(header) % 4:
        header += " "

    splat = BytesIO()
    splat.write(b"splt")
    splat.write(len(header).to_bytes(4, 'little'))
    splat.write(bytearray(header, 'utf-8'))
    splat.write(means_bins.tobytes())
    splat.write(scales_bins.tobytes())
    splat.write(base_buffer)

    splat = splat.getvalue()
    while len(splat) % 4:
        splat += b'\0'
    return splat


def save_splat_file(splat_data, output_path):
    with open(output_path, "wb") as f:
        f.write(splat_data)


def main():
    parser = argparse.ArgumentParser(
        description="Convert Nerfstudio spirulae-splat checkpoint files to SPLAT format.")
    parser.add_argument(
        "input_files", nargs="+", help="The input Nerfstudio spirulae-splat checkpoint files to process."
    )
    parser.add_argument(
        "--output", "-o", default="model.splat", help="The output SPLAT file."
    )
    args = parser.parse_args()
    for input_file in args.input_files:
        print(f"Processing {input_file}...", end='\n\n')
        splat_data = process_ckpt_to_splat(input_file)
        output_file = (
            args.output if len(args.input_files) == 1 else input_file + ".splat"
        )
        save_splat_file(splat_data, output_file)
        print(f"Saved {output_file}")


if __name__ == "__main__":
    main()
