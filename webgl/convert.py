#!/usr/bin/python

import torch
import numpy as np
import argparse
from io import BytesIO
import json


def quantize_tensor(tensor, n_bits, maxiter=0):
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


def pack_components(config, components):
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


def bufferview_psa(buffer_views):
    total = 0
    for view in buffer_views:
        view['byteOffset'] = total
        total += view['byteLength']
    return buffer_views


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
        anisotropies_bins, anisotropies_q = quantize_tensor(anisotropies, 8)
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
    base_config = {
        "bufferView": len(buffer_views),
        "length": len(weight),
        "componentLength": int(np.ceil((36+20+32+16+6+18)/8)),  # in bytes
        "componentViews": [
            { "key": "means", "type": "quat3", "bitLength": 36, "bitOffset": 0, "quatBufferView": 0 },
            { "key": "scales", "type": "quat2", "bitLength": 20, "bitOffset": 36, "quatBufferView": 1 },
            { "key": "quats", "type": "quat4", "bitLength": 32, "bitOffset": 36+20, "quatBufferView": 2 },
            { "key": "anisotropies", "type": "quat2", "bitLength": 16, "bitOffset": 36+20+32, "quatBufferView": 3 },
            { "key": "opacities", "type": "quat", "bitLength": 6, "bitOffset": 36+20+32+16, "quatBufferView": 4 },
            { "key": "features_dc", "type": "quat3", "bitLength": 18, "bitOffset": 36+20+32+16+6, "quatBufferView": 5 },
        ]
    }
    base_buffer = pack_components(base_config, [
        means_q, scales_q, quats_q, anisotropies_q,
        opacities_q, features_dc_q,
    ])

    print("Packing harmonics...")
    harmonics_config = {
        "bufferView": len(buffer_views)+1,
        "length": len(weight),
        "componentLength": int(np.ceil((63*4+45*3)/8)),  # in bytes
        "componentViews": [
            { "key": "features_ch", "type": "quat63", "bitLength": 63*4, "bitOffset": 0, "quatBufferView": 6 },
            { "key": "features_sh", "type": "quat45", "bitLength": 45*3, "bitOffset": 63*4, "quatBufferView": 7 },
        ]
    }
    harmonics_buffer = pack_components(harmonics_config, [
        features_ch_q.reshape((len(features_ch), -1)),
        features_sh_q.reshape((len(features_sh), -1))
    ])

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
            "sh_degree": 3,
            "ch_degree_r": 3,
            "ch_degree_phi": 3,
            "background_color": background_color.cpu().numpy().tolist()
        },
        "buffer": {
            "byteLength": len(buffer),
        },
        "primitives": {
            "base": base_config,
            "harmonics": harmonics_config,
        },
        "bufferViews": bufferview_psa(buffer_views + [
            { "byteLength": len(base_buffer) },
            { "byteLength": len(harmonics_buffer) },
        ])
    }

    header = json.dumps(header)
    while len(header) % 4:
        header += " "

    splat = BytesIO()
    splat.write(b"splt")
    splat.write(len(header).to_bytes(4, 'little'))
    splat.write(bytearray(header, 'utf-8'))
    splat.write(buffer)
    splat = splat.getvalue()
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
