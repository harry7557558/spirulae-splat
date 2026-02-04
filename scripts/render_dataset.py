#!/usr/bin/env python3
"""Render a 3DGS PLY file given camera pose transforms.json file."""

import os
import json
import numpy as np
import cv2

from spirulae_splat.viewer.model import SplatModel
from spirulae_splat.viewer.camera import Camera

import threading
from concurrent.futures import ThreadPoolExecutor


def load_ply(filename: str) -> SplatModel:
    if not os.path.exists(filename):
        raise FileNotFoundError(f"PLY file '{filename}' not found.")
    if not filename.lower().endswith('.ply'):
        raise ValueError(f"File '{filename}' is not a PLY file.")

    model = SplatModel(filename)
    model.flip_yz = True
    return model


def load_transforms(transform_path: str):

    if os.path.isdir(transform_path):
        transform_path = os.path.join(transform_path, "transforms.json")
    if not os.path.exists(transform_path):
        raise FileNotFoundError(f"Transform file '{transform_path}' not found.")

    with open(transform_path, 'r') as f:
        transforms = json.load(f)

    applied_transform = np.array(transforms.get('applied_transform', np.eye(4).tolist()))
    if len(applied_transform) == 3:
        applied_transform = np.vstack([applied_transform, [0, 0, 0, 1]])

    def update_intrins(in_obj, out_obj):
        out_obj = {**out_obj}

        if 'camera_model' in in_obj:
            out_obj['camera_model'] = {
                'SIMPLE_PINHOLE': "OPENCV",
                'PINHOLE': "OPENCV",
                'SIMPLE_RADIAL': "OPENCV",
                'SIMPLE_RADIAL_FISHEYE': "OPENCV_FISHEYE",
                'RADIAL': "OPENCV",
                'RADIAL_FISHEYE': "OPENCV_FISHEYE",
                'OPENCV': "OPENCV",
                'OPENCV_FISHEYE': "OPENCV_FISHEYE",
                'THIN_PRISM_FISHEYE': "OPENCV_FISHEYE",
            }[in_obj['camera_model']]
        elif 'camera_model' not in out_obj:
            out_obj['camera_model'] = 'OPENCV'

        for key in 'w h fl_x fl_y cx cy k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2'.split():
            if key in in_obj:
                out_obj[key] = in_obj[key]
            elif key not in out_obj:
                out_obj[key] = 0.0

        return out_obj

    global_intrins = update_intrins(transforms, {})

    results = {}  # type: dict[str, dict]
    for frame in transforms.get('frames', []):
        file_path = frame.get('file_path', None)
        camera = update_intrins(frame, global_intrins)
        transform_matrix = np.array(frame.get('transform_matrix', np.eye(4).tolist()))
        # transform_matrix = np.linalg.inv(transform_matrix)
        camera['transform_matrix'] = transform_matrix.tolist()

        if file_path not in results:
            results[file_path] = camera
        else:
            raise ValueError(f"Duplicate frame file_path '{file_path}' in transforms.")
    return results


def render_dataset(
    ply_path: str,
    transform_path: str,
    output_dir: str,
    overwrite: bool=False,
):
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading model from '{ply_path}'...")
    model = load_ply(ply_path)

    print(f"Loading transforms from '{transform_path}'...")
    transforms = load_transforms(transform_path)

    model.return_torch = True
    lock = threading.Lock()

    def render_frame(args):
        file_path, frame = args
        if 'images' in os.path.normpath(file_path).split(os.sep):
            file_path = os.path.relpath(file_path, 'images')
        out_path = os.path.join(output_dir, file_path)
        if not overwrite and os.path.exists(out_path):
            print(f"Output '{out_path}' exists, skipping (use --overwrite to force).")
            return
        os.makedirs(os.path.dirname(out_path), exist_ok=True)

        camera = Camera(config_path=frame)
        c2w = frame['transform_matrix']
        with lock:
            image = model.render(camera, c2w)
        image = image.cpu().numpy()
        image = (np.clip(image, 0.0, 1.0) * 255).astype(np.uint8)

        cv2.imwrite(out_path, cv2.cvtColor(image, cv2.COLOR_RGB2BGR))

        print(f"Rendered frame from '{file_path}' to '{out_path}'")

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(render_frame, sorted(transforms.items()))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Render a 3DGS PLY file given camera pose transforms.json file.")
    parser.add_argument("--transforms", required=True, help="Path to the transforms.json file.")
    parser.add_argument("--output", "-o", required=True, help="The output image directory.")
    parser.add_argument("--ply", required=True, help="Path to the 3DGS PLY file.")
    parser.add_argument("--overwrite", action="store_true", help="Whether to overwrite existing output directory.")
    args = parser.parse_args()

    render_dataset(
        ply_path=args.ply,
        transform_path=args.transforms,
        output_dir=args.output,
        overwrite=args.overwrite,
    )
