#!/usr/bin/env python3
"""
NeRF Studio Dataset Downscaler

This script downscales a NeRF Studio dataset by reducing image resolution,
adjusting camera intrinsics, and optionally downscaling point cloud data.
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, Any, Optional
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


def validate_args(args):
    """Validate command line arguments."""
    if not (0 < args.image_scale <= 1):
        raise ValueError("Image scale factor must be in range (0, 1]")
    
    if args.ply_scale is not None and not (0 < args.ply_scale <= 1):
        raise ValueError("PLY scale factor must be in range (0, 1]")
    
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


def get_jpeg_quality(image_path: str) -> Optional[int]:
    """Get JPEG quality from image file if it's a JPEG."""
    if not image_path.lower().endswith(('.jpg', '.jpeg')):
        return None
    
    try:
        # Read JPEG quality using OpenCV (approximation)
        # This is a rough estimation as OpenCV doesn't directly expose quality
        img = cv2.imread(image_path)
        if img is None:
            return None
        
        # Encode with different qualities and find closest match
        _, original_encoded = cv2.imencode('.jpg', img)
        original_size = len(original_encoded)
        
        for quality in range(95, 50, -5):
            _, test_encoded = cv2.imencode('.jpg', img, [cv2.IMWRITE_JPEG_QUALITY, quality])
            if len(test_encoded) <= original_size * 1.1:  # Allow 10% tolerance
                return quality
        
        return 90  # Default fallback
    except:
        return 90  # Default fallback


def downscale_image(input_path: str, output_path: str, scale_factor: float, 
                   jpeg_quality: Optional[int] = None):
    """Downscale an image and save it."""
    # Read image
    img = cv2.imread(input_path)
    if img is None:
        raise ValueError(f"Could not read image: {input_path}")
    
    if scale_factor == 1.0 and jpeg_quality is None:
        # Just copy the file
        shutil.copy2(input_path, output_path)
        return img.shape[1], img.shape[0]  # width, height
    
    # Calculate new dimensions
    original_height, original_width = img.shape[:2]
    new_width = int(original_width * scale_factor)
    new_height = int(original_height * scale_factor)
    
    # Downscale image
    if scale_factor != 1.0:
        img = cv2.resize(img, (new_width, new_height), interpolation=cv2.INTER_AREA)
    
    # Determine output format and quality
    if jpeg_quality is not None:
        # Force JPEG output
        output_path = str(Path(output_path).with_suffix('.jpg'))
        cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, jpeg_quality])
    else:
        # Keep original format
        if input_path.lower().endswith(('.jpg', '.jpeg')):
            # Try to preserve original JPEG quality
            original_quality = get_jpeg_quality(input_path)
            if original_quality:
                cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, original_quality])
            else:
                cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, 90])
        else:
            # Non-JPEG format
            if output_path.lower().endswith(('.jpg', '.jpeg')):
                cv2.imwrite(output_path, img, [cv2.IMWRITE_JPEG_QUALITY, 90])
            else:
                cv2.imwrite(output_path, img)
    
    return new_width, new_height


def adjust_intrinsics(intrinsics: Dict[str, Any], scale_factor: float, 
                     new_width: int, new_height: int) -> Dict[str, Any]:
    """Adjust camera intrinsics for downscaled image."""
    adjusted = intrinsics.copy()
    
    # Scale focal lengths and principal point
    if 'fl_x' in adjusted:
        adjusted['fl_x'] *= scale_factor
    if 'fl_y' in adjusted:
        adjusted['fl_y'] *= scale_factor
    if 'cx' in adjusted:
        adjusted['cx'] *= scale_factor
    if 'cy' in adjusted:
        adjusted['cy'] *= scale_factor
    
    # Update image dimensions
    adjusted['w'] = new_width
    adjusted['h'] = new_height

    return adjusted


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
    parser = argparse.ArgumentParser(description="Downscale NeRF Studio dataset")
    parser.add_argument("input_dir", help="Input dataset directory")
    parser.add_argument("output_dir", help="Output dataset directory (must be empty or non-existing)")
    parser.add_argument("image_scale", type=float, help="Image downscale factor (0, 1]")
    parser.add_argument("--ply_scale", type=float, default=1.0, help="PLY downscale factor (0, 1]")
    parser.add_argument("--jpeg_quality", type=int, help="JPEG quality (1-100)")
    
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
    print(f"Image scale factor: {args.image_scale}")
    if args.ply_scale:
        print(f"PLY scale factor: {args.ply_scale}")
    if args.jpeg_quality:
        print(f"JPEG quality: {args.jpeg_quality}")
    
    # Create output directory structure
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check for global intrinsics
    global_intrinsics = {}
    intrinsics_keys = ['w', 'h', 'fl_x', 'fl_y', 'cx', 'cy', 'k1', 'k2', 'k3', 'k4', 'p1', 'p2', 'camera_model']
    for key in intrinsics_keys:
        if key in transforms:
            global_intrinsics[key] = transforms[key]
    
    # Process images and update frames
    processed_frames = []

    def process_one_frame(frame):
        file_path = frame['file_path']
        
        # Ensure file path is within dataset directory
        full_input_path = os.path.join(args.input_dir, file_path)
        if not os.path.exists(full_input_path):
            print(f"Warning: Image not found: {full_input_path}")
            return
        
        # Create output directory structure
        output_file_path = file_path
        if args.jpeg_quality is not None and not file_path.lower().endswith(('.jpg', '.jpeg')):
            # Convert to JPEG
            output_file_path = str(Path(file_path).with_suffix('.jpg'))
        
        full_output_path = os.path.join(args.output_dir, output_file_path)
        os.makedirs(os.path.dirname(full_output_path), exist_ok=True)
        
        # Downscale image
        new_width, new_height = downscale_image(
            full_input_path, full_output_path, args.image_scale, args.jpeg_quality
        )
        
        # Update frame data
        new_frame = frame.copy()
        new_frame['file_path'] = output_file_path
        
        # Adjust intrinsics (frame-specific or use global)
        frame_intrinsics = {}
        for key in intrinsics_keys:
            if key in frame:
                frame_intrinsics[key] = frame[key]
            elif key in global_intrinsics:
                frame_intrinsics[key] = global_intrinsics[key]
        
        # Apply adjustments
        adjusted_intrinsics = adjust_intrinsics(frame_intrinsics, args.image_scale, new_width, new_height)
        
        # Update frame with adjusted intrinsics (only if they were in the original frame)
        is_non_global = any([k not in global_intrinsics for k in intrinsics_keys[:6]])
        for key, value in adjusted_intrinsics.items():
            if key in frame or is_non_global:
                new_frame[key] = value
        
        return new_frame
    
    with ThreadPoolExecutor() as executor:
        for new_frame in tqdm(
            executor.map(process_one_frame, transforms['frames']),
            total=len(transforms["frames"]),
            desc="Downscaling images..."
        ):
            if new_frame is not None:
                processed_frames.append(new_frame)

    # Update global intrinsics if they exist
    if global_intrinsics:
        # Use dimensions from first processed frame as reference
        if processed_frames:
            ref_frame = processed_frames[0]
            ref_width = ref_frame.get('w', global_intrinsics.get('w', 0))
            ref_height = ref_frame.get('h', global_intrinsics.get('h', 0))
            new_width = int(ref_width * args.image_scale)
            new_height = int(ref_height * args.image_scale)
            
            adjusted_global = adjust_intrinsics(global_intrinsics, args.image_scale, new_width, new_height)
            for key, value in adjusted_global.items():
                transforms[key] = value
    
    # Update transforms with processed frames
    transforms['frames'] = processed_frames
    
    # Process PLY file if it exists
    if 'ply_file_path' in transforms and args.ply_scale is not None:
        ply_input_path = os.path.join(args.input_dir, transforms['ply_file_path'])
        ply_output_path = os.path.join(args.output_dir, transforms['ply_file_path'])
        
        if os.path.exists(ply_input_path):
            os.makedirs(os.path.dirname(ply_output_path), exist_ok=True)
            downscale_ply(ply_input_path, ply_output_path, args.ply_scale)
        else:
            print(f"Warning: PLY file not found: {ply_input_path}")
    elif 'ply_file_path' in transforms:
        # Copy PLY file without downscaling
        ply_input_path = os.path.join(args.input_dir, transforms['ply_file_path'])
        ply_output_path = os.path.join(args.output_dir, transforms['ply_file_path'])
        
        if os.path.exists(ply_input_path):
            os.makedirs(os.path.dirname(ply_output_path), exist_ok=True)
            shutil.copy2(ply_input_path, ply_output_path)
    
    # Save updated transforms
    output_transforms_path = os.path.join(args.output_dir, "transforms.json")
    save_transforms(transforms, output_transforms_path)
    
    print(f"Dataset downscaled successfully!")
    print(f"Processed {len(processed_frames)} frames")


if __name__ == "__main__":
    main()
