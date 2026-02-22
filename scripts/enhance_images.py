#!/usr/bin/env python3

import os
from pathlib import Path
import numpy as np
import spirulae_splat.modules.enhancer as enhancer

from typing import Optional

import math

import cv2


def process_single_image(x: np.ndarray, max_tile_size: Optional[int] = None) -> np.ndarray:
    # x should be [H, W, C]
    num_splits_x, num_splits_y = 1, 1
    if max_tile_size is not None:
        if x.shape[0] * x.shape[1] > max_tile_size * max_tile_size:
            num_splits_x = math.ceil(x.shape[0] / max_tile_size)
            num_splits_y = math.ceil(x.shape[1] / max_tile_size)
    if num_splits_x > 1 or num_splits_y > 1:
        tile_size_x = math.ceil(x.shape[0] / num_splits_x)
        tile_size_y = math.ceil(x.shape[1] / num_splits_y)
        tiles = []
        for ix in range(num_splits_x):
            for iy in range(num_splits_y):
                x0 = ix * tile_size_x
                y0 = iy * tile_size_y
                x1 = min((ix+1) * tile_size_x, x.shape[0])
                y1 = min((iy+1) * tile_size_y, x.shape[1])
                tile = x[x0:x1, y0:y1, :]
                tiles.append(tile)
        processed_tiles = []
        for tile in tiles:
            processed_tile = enhancer.infer(tile)
            processed_tiles.append(processed_tile)
        # reconstruct
        y = np.zeros_like(x)
        index = 0
        for ix in range(num_splits_x):
            for iy in range(num_splits_y):
                x0 = ix * tile_size_x
                y0 = iy * tile_size_y
                x1 = min((ix+1) * tile_size_x, y.shape[0])
                y1 = min((iy+1) * tile_size_y, y.shape[1])
                y[x0:x1, y0:y1, :] = processed_tiles[index]
                index += 1
    else:
        y = enhancer.infer(x)
    return y

def process_folder(in_path, out_path, max_tile_size: Optional[int] = None):
    import os
    from pathlib import Path
    os.makedirs(out_path, exist_ok=True)
    from tqdm import tqdm
    
    # Recursively find all image files
    image_extensions = {'.jpg', '.jpeg', '.png', '.bmp', '.tiff', '.tif'}
    image_files = []
    
    for root, dirs, files in os.walk(in_path):
        for fn in files:
            if Path(fn).suffix.lower() in image_extensions:
                rel_path = os.path.relpath(os.path.join(root, fn), in_path)
                image_files.append(rel_path)
    
    # Process each image, maintaining directory structure
    for rel_path in tqdm(sorted(image_files)):
        src_file = os.path.join(in_path, rel_path)
        dst_file = os.path.join(out_path, rel_path)
        
        # Create output subdirectories if needed
        os.makedirs(os.path.dirname(dst_file), exist_ok=True)
        
        x = cv2.imread(src_file, cv2.IMREAD_UNCHANGED)
        if x is None:
            continue
        x = cv2.cvtColor(x, cv2.COLOR_BGR2RGB)
        y = process_single_image(x, max_tile_size=max_tile_size)
        y = cv2.cvtColor(y, cv2.COLOR_RGB2BGR)
        cv2.imwrite(dst_file, y)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Enhance a dataset of images using a pre-trained enhancer model.")
    parser.add_argument("input_dir", nargs=1, help="The dataset folder. Should contain a subfolder named `images`.")
    parser.add_argument("--backup_dir", type=str, default="images_original", help="The name of the sub-folder to back up original images to.")
    parser.add_argument("--max_tile_size", type=int, default=1024, help="The maximum tile size to use when enhancing large images. Larger values use more GPU memory.")
    args = parser.parse_args()

    already_processed_mark_filename = ".enhance_images_done"

    dataset_dir = args.input_dir[0]
    image_dir = Path(dataset_dir) / "images"
    if (Path(dataset_dir) / already_processed_mark_filename).exists():
        print(f"Dataset at {dataset_dir} already processed, skipping.")
        exit(0)

    backup_dir = None
    if args.backup_dir != "":
        backup_dir = Path(dataset_dir) / args.backup_dir
        os.makedirs(backup_dir, exist_ok=True)
        
        # Recursively move all images to backup directory, maintaining structure
        for root, dirs, files in os.walk(image_dir):
            for fn in sorted(files):
                src_file = os.path.join(root, fn)
                rel_path = os.path.relpath(src_file, image_dir)
                dst_file = os.path.join(backup_dir, rel_path)
                os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                os.rename(src_file, dst_file)
        
        # Remove empty directories
        for root, dirs, files in os.walk(image_dir, topdown=False):
            for d in dirs:
                dir_path = os.path.join(root, d)
                try:
                    os.rmdir(dir_path)
                except OSError:
                    pass

    in_dir = backup_dir if backup_dir != "" else image_dir
    out_dir = image_dir
    process_folder(in_dir, out_dir, max_tile_size=args.max_tile_size)

    Path(dataset_dir, already_processed_mark_filename).touch()
