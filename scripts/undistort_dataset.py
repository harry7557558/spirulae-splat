#!/usr/bin/env python3
"""
Nerfstudio Dataset Undistort Script

This script undistorts a Nerfstudio dataset.
For pinhole camera, it removes distortion from images and updates camera intrinsics accordingly.
For fisheye camera, it undistorts to an ideal equidistant fisheye model without distortion parameters, and updates intrinsics to match the new model.
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Literal
import random

from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

import cv2
import numpy as np

try:
    from plyfile import PlyData, PlyElement
except ImportError:
    print("Error: plyfile library is required. Install with: pip install plyfile")
    sys.exit(1)

import torch
from spirulae_splat.splat.cuda import undistort_image as ssplat_undistort_image


def validate_args(args):
    """Validate command line arguments."""
    if args.jpeg_quality is not None and not (1 <= args.jpeg_quality <= 100):
        raise ValueError("JPEG quality must be in range [1, 100]")
    
    if not os.path.exists(args.input_dir):
        raise ValueError(f"Input directory does not exist: {args.input_dir}")
    
    transforms_path = os.path.join(args.input_dir, "transforms.json")
    if not os.path.exists(transforms_path):
        raise ValueError(f"transforms.json not found in input directory: {args.input_dir}")
    
    # Check output directory
    if os.path.exists(args.output_dir):
        if os.listdir(args.output_dir):
            raise ValueError(f"Output directory is non-empty: {args.output_dir}")
    else:
        os.makedirs(args.output_dir, exist_ok=True)


def load_transforms(transforms_path: str) -> Dict[str, Any]:
    """Load transforms.json file."""
    with open(transforms_path, 'r') as f:
        return json.load(f)


def save_transforms(transforms: Dict[str, Any], output_path: str):
    """Save transforms.json file."""
    with open(output_path, 'w') as f:
        json.dump(transforms, f, indent=2)


def undistort_image(input_path: str, output_path: str, intrins: dict, 
                   jpeg_quality: Optional[int] = None):
    """Undistort an image and save it."""
    # Read image
    img = cv2.imread(input_path, -1)
    if img is None:
        raise ValueError(f"Could not read image: {input_path}")
    if img.ndim == 2:
        img = img[..., np.newaxis]  # [H, W, 1]
    
    # Undistort image
    sw = img.shape[1] / intrins['w']
    sh = img.shape[0] / intrins['h']
    intrins = (
        "pinhole",
        (intrins['fl_x']*sw, intrins['fl_y']*sh, intrins['cx']*sw, intrins['cy']*sh),
        tuple(intrins.get(key, 0.0) for key in "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split())
    )
    dtype = img.dtype
    img = torch.from_numpy(img).unsqueeze(0).float().cuda()  # [1, H, W, C]
    img = ssplat_undistort_image(img, *intrins)
    img = img.squeeze(0).cpu().numpy().astype(dtype)  # [H, W, C]

    # Determine output format and quality
    if jpeg_quality is not None:
        # Force JPEG output
        output_path = str(Path(output_path).with_suffix('.jpg'))
        cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, jpeg_quality])
    else:
        # Keep original format
        if input_path.lower().endswith(('.jpg', '.jpeg')):
            # Try to preserve original JPEG quality
            cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, 95])
        else:
            # Non-JPEG format
            if output_path.lower().endswith(('.jpg', '.jpeg')):
                cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, 90])
            else:
                cv2.imwrite(output_path, img)

    # TODO: properly generate from undistort function
    black_pixel_mask = (img == 0).all(axis=-1)
    mask = (~black_pixel_mask).astype(np.uint8)  # 1 for valid pixels, 0 for black borders
    if mask.ndim == 2:
        mask = mask[..., np.newaxis]  # [H, W, 1]
    if np.count_nonzero(mask) < mask.size:
        return mask
    return None


def adjust_intrinsics(intrinsics: Dict[str, Any]) -> None:
    """Adjust camera intrinsics in place."""
    
    for dist_coeffs in 'k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2'.split():
        if dist_coeffs in intrinsics:
            intrinsics.pop(dist_coeffs)

    if 'camera_model' in intrinsics:
        intrinsics['camera_model'] = {
            'SIMPLE_PINHOLE': "OPENCV",
            'PINHOLE': "OPENCV",
            'SIMPLE_RADIAL': "OPENCV",
            'SIMPLE_RADIAL_FISHEYE': "OPENCV_FISHEYE",
            'RADIAL': "OPENCV",
            'RADIAL_FISHEYE': "OPENCV_FISHEYE",
            'OPENCV': "OPENCV",
            'OPENCV_FISHEYE': "OPENCV_FISHEYE",
            'THIN_PRISM_FISHEYE': "OPENCV_FISHEYE",
        }[intrinsics['camera_model']]


def downscale_ply(input_path: str, output_path: str, scale_factor: float):
    """Downscale PLY point cloud by random sampling."""
    if scale_factor == 1.0:
        shutil.copy2(input_path, output_path)
        return
    
    # Read PLY file
    plydata = PlyData.read(input_path)
    vertex_data = plydata['vertex'].data
    
    # Calculate number of points to keep
    original_count = len(vertex_data)
    target_count = int(np.floor(original_count * scale_factor))
    
    if target_count >= original_count:
        shutil.copy2(input_path, output_path)
        return
    
    # Random sampling
    random.seed(42)
    indices = random.sample(range(original_count), target_count)
    sampled_data = vertex_data[indices]
    
    # Create new PLY element
    sampled_vertex = PlyElement.describe(sampled_data, 'vertex')
    
    # Write downscaled PLY
    PlyData([sampled_vertex]).write(output_path)
    
    print(f"PLY downscaled: {original_count} -> {target_count} points")


def main():
    parser = argparse.ArgumentParser(description="Undistort Nerfstudio dataset")
    parser.add_argument("input_dir", help="Input dataset directory")
    parser.add_argument("output_dir", help="Output dataset directory (must be empty or non-existing)")
    parser.add_argument("--jpeg_quality", type=int, help="JPEG quality (1-100)")
    # parser.add_argument("--max_workers", type=Optional[int], default=None, help="Maximum number of worker threads for parallel processing")
    
    args = parser.parse_args()
    
    try:
        validate_args(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Load transforms
    transforms_path = os.path.join(args.input_dir, "transforms.json")
    transforms = load_transforms(transforms_path)
    
    print(f"Processing dataset: {args.input_dir} -> {args.output_dir}")
    if args.jpeg_quality:
        print(f"JPEG quality: {args.jpeg_quality}")
    
    # Create output directory structure
    os.makedirs(args.output_dir, exist_ok=True)
    
    intrinsics_keys = ['w', 'h', 'fl_x', 'fl_y', 'cx', 'cy', 'k1', 'k2', 'k3', 'k4', 'p1', 'p2', 'camera_model']

    # Process images and update frames
    processed_frames = []

    def process_one_frame(frame, key: Literal['file_path', 'mask_path', 'depth_file_path', 'normal_file_path', None] = None):
        if key is None:
            return_val = None
            for key in ['file_path', 'mask_path', 'depth_file_path', 'normal_file_path'][::-1]:
                return_val = process_one_frame(frame, key)
    
            # handle black borders in undistorted images by generating masks and combining with existing masks if they exist
            mask = None
            if 'mask' in return_val:
                mask = return_val.pop('mask')
            if mask is not None:
                if 'mask_path' in return_val:
                    mask_input_path = os.path.join(args.input_dir, return_val['mask_path'])
                    loaded_mask = cv2.imread(mask_input_path, cv2.IMREAD_UNCHANGED)
                    if loaded_mask is not None and loaded_mask.size >= 4:
                        loaded_mask = cv2.resize(loaded_mask, (mask.shape[1], mask.shape[0]), interpolation=cv2.INTER_NEAREST)
                        if loaded_mask.ndim == 2:
                            loaded_mask = loaded_mask[..., np.newaxis]  # [H, W, 1]
                        if loaded_mask.shape[2] > 1:
                            loaded_mask = (loaded_mask != 0).all(axis=-1, keepdims=True).astype(np.uint8)
                        mask = np.logical_and(mask, loaded_mask).astype(np.uint8)
                if 'mask_path' not in return_val:
                    # images/cam0/00022.jpg -> masks/cam0/00022.jpg.png
                    mask_file_path = os.path.join('masks', return_val['file_path'].lstrip('.').lstrip(os.path.sep).lstrip('images').lstrip(os.path.sep) + ".png")
                    return_val['mask_path'] = mask_file_path
                mask_output_path = os.path.join(args.output_dir, return_val['mask_path'])
                os.makedirs(os.path.dirname(mask_output_path), exist_ok=True)
                cv2.imwrite(mask_output_path, mask.astype(np.uint8) * 255)

            return return_val

        if key not in frame:
            return frame
        file_path = frame[key]
        
        # Ensure file path is within dataset directory
        full_input_path = os.path.join(args.input_dir, file_path)
        if not os.path.exists(full_input_path):
            print(f"Warning: Image not found: {full_input_path}")
            return frame
        
        # Create output directory structure
        output_file_path = file_path
        if args.jpeg_quality is not None and not file_path.lower().endswith(('.jpg', '.jpeg')):
            # Convert to JPEG
            output_file_path = str(Path(file_path).with_suffix('.jpg'))
        
        full_output_path = os.path.join(args.output_dir, output_file_path)
        os.makedirs(os.path.dirname(full_output_path), exist_ok=True)
        
        # Update frame data
        new_frame = frame.copy()
        new_frame[key] = output_file_path
        
        # Adjust intrinsics (frame-specific or use global)
        intrinsics_for_undistortion = {}
        for key in intrinsics_keys:
            if key in frame:
                intrinsics_for_undistortion[key] = frame[key]
            elif key in transforms:
                intrinsics_for_undistortion[key] = transforms[key]
        
        # Undistort image
        mask = undistort_image(
            full_input_path, full_output_path, intrinsics_for_undistortion, args.jpeg_quality
        )
        if mask is not None:
            if 'mask' not in new_frame:
                new_frame['mask'] = mask
            else:
                new_frame['mask'] = np.logical_and(new_frame['mask'], mask).astype(np.uint8)
        
        # Apply adjustments
        adjust_intrinsics(new_frame)
        
        return new_frame
    
    # with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
    with ThreadPoolExecutor() as executor:
        for new_frame in tqdm(
            executor.map(process_one_frame, transforms['frames']),
            total=len(transforms["frames"]),
            desc="Undistorting images..."
        ):
            if new_frame is not None:
                processed_frames.append(new_frame)

    # Update global intrinsics if they exist
    adjust_intrinsics(transforms)
    
    # Update transforms with processed frames
    transforms['frames'] = processed_frames
    
    # Process PLY file if it exists
    if 'ply_file_path' in transforms:
        # Copy PLY file
        ply_input_path = os.path.join(args.input_dir, transforms['ply_file_path'])
        ply_output_path = os.path.join(args.output_dir, transforms['ply_file_path'])
        
        if os.path.exists(ply_input_path):
            os.makedirs(os.path.dirname(ply_output_path), exist_ok=True)
            shutil.copy2(ply_input_path, ply_output_path)
    
    # Save updated transforms
    output_transforms_path = os.path.join(args.output_dir, "transforms.json")
    save_transforms(transforms, output_transforms_path)
    
    print(f"Dataset undistorted successfully!")
    print(f"Processed {len(processed_frames)} frames")


if __name__ == "__main__":
    main()
