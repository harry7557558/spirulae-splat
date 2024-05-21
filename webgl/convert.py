#!/usr/bin/python

import torch
import numpy as np
import argparse
from io import BytesIO


def process_ckpt_to_splat(file_path):
    checkpoint = torch.load(file_path)
    pipeline = checkpoint['pipeline']
    
    features_dc = pipeline['_model.gauss_params.features_dc']
    features_rest = pipeline['_model.gauss_params.features_rest']
    means = pipeline['_model.gauss_params.means']
    opacities = pipeline['_model.gauss_params.opacities']
    quats = pipeline['_model.gauss_params.quats']
    scales = pipeline['_model.gauss_params.scales']

    # means -= torch.mean(means, axis=0)
    SH_C0 = 0.28209479177387814
    features_dc = 0.5 + features_dc * SH_C0
    opacities = torch.sigmoid(opacities)+0.5
    scales = torch.cat((torch.exp(scales), torch.zeros_like(scales[:,0:1])), dim=1)

    means = means.cpu().numpy()
    scales = scales.cpu().numpy()
    quats = quats.cpu().numpy()
    features_dc = features_dc.cpu().numpy()
    opacities = opacities.cpu().numpy()

    center = np.mean(means, axis=0)
    cov_matrix = np.cov(means, rowvar=False)
    mahalanobis_dist = np.sqrt(np.sum(np.dot(
        (means-center), np.linalg.inv(cov_matrix)) * (means-center), axis=1))
    mask = mahalanobis_dist < 1.5

    center = np.mean(means, axis=0)
    means -= center

    weight = scales[:,0] * scales[:,1] * opacities[:,0]
    sorted_indices = np.argsort(-weight)
    buffer = BytesIO()
    for idx in sorted_indices:
        if not mask[idx]:
            continue
        position = means[idx]
        scale = scales[idx]
        rot =  quats[idx]
        color = np.concatenate((features_dc[idx], opacities[idx]))
        buffer.write(position.tobytes())
        buffer.write(scale.tobytes())
        buffer.write((color*255+0.5).clip(0, 255).astype(np.uint8).tobytes())
        buffer.write(
            ((rot / np.linalg.norm(rot)) * 128 + 128.5)
            .clip(0, 255)
            .astype(np.uint8)
            .tobytes()
        )

    return buffer.getvalue()


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
        "--output", "-o", default="output.splat", help="The output SPLAT file."
    )
    args = parser.parse_args()
    for input_file in args.input_files:
        print(f"Processing {input_file}...")
        splat_data = process_ckpt_to_splat(input_file)
        output_file = (
            args.output if len(args.input_files) == 1 else input_file + ".splat"
        )
        save_splat_file(splat_data, output_file)
        print(f"Saved {output_file}")


if __name__ == "__main__":
    main()
