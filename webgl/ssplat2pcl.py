import struct
import json
import numpy as np
from scipy.stats import qmc
import argparse

def load_model(filename):
    with open(filename, 'rb') as fs:
        fs.seek(4)
        size_buffer = fs.read(4)
        if len(size_buffer) != 4:
            raise IOError("Failed to read header size.")
        header_size = struct.unpack('<I', size_buffer)[0]
        json_buffer = fs.read(header_size)
        if len(json_buffer) != header_size:
            raise IOError("Failed to read JSON data.")
        json_string = json_buffer.decode('utf-8')
        json_object = json.loads(json_string)
        raw_bytes = fs.read()
    return unpack_model(json_object, raw_bytes)


def unpack_model(header, buffer):
    model = {
        "header": header
    }

    buffers = []
    for view in header["bufferViews"]:
        start = view["byteOffset"]
        end = start + view["byteLength"]
        buffers.append(buffer[start:end])

    primitives = header.get("primitives", {})
    for key, primitive_obj in primitives.items():
        if isinstance(primitive_obj, dict):
            buffer_slice = buffers[primitive_obj["bufferView"]]
            components = unpack_components(primitive_obj, buffer_slice, buffers)
            model[key] = components

    model["header"] = header
    return model


def unpack_components(config, packed_data, buffers):
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
    c_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(c_module)

    config_json = json.loads(json.dumps(config))
    unpacked_data = c_module.unpack_components(config_json, packed_data, buffers)

    return unpacked_data


def weighted_sample_low_discrepancy(weights, m, sampler):
    weights = np.array(weights)
    cumulative_weights = np.cumsum(weights)
    total_weight = cumulative_weights[-1]
    samples = sampler.random(n=m)
    scaled_samples = samples * total_weight
    indices = np.searchsorted(cumulative_weights, scaled_samples)
    return indices.flatten()

def quat_to_rotmat(quats):
    assert quats.shape[-1] == 4, quats.shape
    w, x, y, z = quats[:,0], quats[:,1], quats[:,2], quats[:,3]
    mat = np.stack([
        1 - 2 * (y**2 + z**2),
        2 * (x * y - w * z),
        2 * (x * z + w * y),
        2 * (x * y + w * z),
        1 - 2 * (x**2 + z**2),
        2 * (y * z - w * x),
        2 * (x * z - w * y),
        2 * (y * z + w * x),
        1 - 2 * (x**2 + y**2),
    ]).T
    return mat.reshape(quats.shape[:-1] + (3, 3))

def process_ssplat_to_pcl(input_file, num_points=1000000):
    model = load_model(input_file)
    base = model["base"]
    means = base["means"].reshape(-1, 3)
    quats = base["quats"].reshape(-1, 4)
    scales = np.exp(base["scales"].reshape(len(means), -1))
    opacs = 1.0/(1.0+np.exp(-base["opacities"]))
    quats /= np.linalg.norm(quats, axis=-1, keepdims=True)
    # weights = scales[:,0] * scales[:,1] * opacs
    weights = opacs

    sampler = qmc.Sobol(d=1, scramble=True, seed=42)
    indices = weighted_sample_low_discrepancy(weights, num_points, sampler)
    scales = scales[indices]
    R = quat_to_rotmat(quats[indices])
    means = means[indices]

    u = 2.0*np.pi * np.random.random(num_points)
    v = np.sqrt(1.0-np.sqrt(1.0-np.random.random(num_points)))
    p = np.stack((v*np.cos(u)*scales[:,0], v*np.sin(u)*scales[:,1], 0.0*v)).T
    p = np.einsum('nij,nj->ni', R, p) + means
    return p

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    ax = plt.axes(projection='3d')
    ax.scatter(p[:,0], p[:,1], p[:,2], s=1)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    plt.show()


def save_ply_file(points, filename):
    fp = open(filename, 'w')
    fp.write("ply\n")
    fp.write("format ascii 1.0\n")
    fp.write(f"element vertex {len(points)}\n", )
    for a in 'xyz':
        fp.write(f"property float {a}\n")
    fp.write("end_header\n")
    for p in points:
        fp.write("{:.6f} {:.6f} {:.6f}\n".format(*p))
    fp.close()


def main():
    parser = argparse.ArgumentParser(
        description="Convert SSPLAT file to PLY point cloud.")
    parser.add_argument("input_files", nargs="+", help="The input SSPLAT files.")
    parser.add_argument("--output", "-o", default="output.ply", help="The output PLY file.")
    parser.add_argument("--num_points", "-n", default="100000", help="Number of points in output PLY file.")
    args = parser.parse_args()
    for input_file in args.input_files:
        print(f"Processing {input_file}...", end='\n\n')
        points = process_ssplat_to_pcl(input_file)
        output_file = (
            args.output if len(args.input_files) == 1 else input_file + ".ply"
        )
        save_ply_file(points, output_file)
        print(f"Saved {output_file}")


if __name__ == "__main__":
    main()
