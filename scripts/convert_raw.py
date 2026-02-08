#!/usr/bin/env python3
"""
DNG to JPG Batch Processor

Processes .dng files from an input folder and saves them as .jpg files
in organized subfolders with different scaling factors.
"""

import os
import sys
import shutil
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Tuple, Optional
import threading

try:
    import rawpy
    import numpy as np
    import cv2
    import piexif
except ImportError as e:
    print(f"Required package not installed: {e}")
    print("Install with: pip install rawpy pillow piexif")
    sys.exit(1)


def aces_tonemap(x: np.ndarray):
    a = 2.51
    b = 0.03
    c = 2.43
    d = 0.59
    e = 0.14
    y = (x * (a * x + b)) / (x * (c * x + d) + e)
    return np.clip(y, 0.0, 1.0)


class DNGProcessor:
    def __init__(self, input_folder: str, output_folder: str, 
                 color_space: str, save_format: str, max_workers: Optional[int] = None):
        self.input_folder = Path(input_folder)
        self.output_folder = Path(output_folder)
        self.color_space = color_space.lower()
        self.save_format = save_format.lower()
        self.max_workers = max_workers or os.cpu_count()
        self.lock = threading.Lock()
        self.processed_count = 0
        self.total_files = 0
        
        # Validate color space
        if self.color_space not in ['srgb', 'aces', 'linear-rgb', 'linear-aces']:
            raise ValueError("Color space must be 'sRGB', 'ACES', 'Linear-RGB', or 'Linear-ACES'")
    
    def setup_output_folders(self) -> None:
        """Setup output folder structure, ensuring it's empty or doesn't exist."""
        if self.output_folder.exists():
            if any(self.output_folder.iterdir()):
                response = input(f"Output folder is not empty. "
                               "Delete contents? (y/N): ")
                if response.lower() != 'y':
                    print("Operation cancelled.")
                    sys.exit(1)
                shutil.rmtree(self.output_folder)
        
        # Create main output folder and subfolders
        self.output_folder.mkdir(parents=True, exist_ok=True)
        for subfolder in ['images', 'images_2', 'images_4']:
            (self.output_folder / subfolder).mkdir(exist_ok=True)
    
    def get_dng_files(self) -> list:
        """Get list of .dng files from input folder."""
        dng_files = list(self.input_folder.glob('*.dng')) + \
                   list(self.input_folder.glob('*.DNG'))
        return sorted(dng_files)
    
    def process_raw_image(self, raw_path: Path) -> Tuple[np.ndarray, dict]:
        """Process a single RAW image with specified parameters."""
        with rawpy.imread(str(raw_path)) as raw:
            
            if self.color_space == 'srgb':
                output_color = rawpy.ColorSpace.sRGB
            elif self.color_space == 'linear-rgb':
                output_color = rawpy.ColorSpace.sRGB
            elif self.color_space == 'aces':
                output_color = rawpy.ColorSpace.ACES  # ACES2065-1
            elif self.color_space == 'linear-aces':
                output_color = rawpy.ColorSpace.ACES
            else:
                raise ValueError(f"Unknown color space: {self.color_space}")

            kwargs = {}
            if self.color_space.startswith('linear-'):
                kwargs['gamma'] = (1.0, 1.0)
            if self.save_format == 'png16':
                kwargs['output_bps'] = 16
            else:
                kwargs['output_bps'] = 8

            rgb_array = raw.postprocess(
                output_color=output_color,
                no_auto_bright=True,
                use_camera_wb=True,
                use_auto_wb=False,
                fbdd_noise_reduction=rawpy.FBDDNoiseReductionMode.Full,
                noise_thr=100,
                **kwargs
            )

            # brightness adjustment
            if False:
                target_l = 0.25
                if False:
                    # linear brightness, can lead to over-exposure
                    lut = target_l * np.arange(2**16) / np.mean(rgb_array)
                    lut = (255*aces_tonemap(lut)).astype(np.uint8)
                else:
                    # gamma, works well for dark area, may reduce contrast
                    k = np.log(target_l) / (np.mean(np.log(np.clip(rgb_array, min=1))) - np.log(2**16))
                    lut = (255*np.linspace(0, 1, 2**16)**k).astype(np.uint8)
                rgb_array = lut[rgb_array]
            
            metadata = {
                'camera_wb': raw.camera_whitebalance,
                'daylight_wb': raw.daylight_whitebalance,
                'color_matrix': raw.color_matrix,
                'raw_image_sizes': raw.raw_image.shape
            }
            
            return rgb_array, metadata
    
    def copy_exif_data(self, source_path: Path, target_path: Path) -> None:
        """Copy EXIF data from source DNG to target JPG."""
        try:
            # Read EXIF from source
            with rawpy.imread(str(source_path)) as raw:
                # Extract basic metadata that can be safely transferred
                exif_dict = {
                    "0th": {},
                    "Exif": {},
                    "GPS": {},
                    "1st": {},
                    "thumbnail": None
                }
                
                # Add basic camera information if available
                try:
                    # These are safe to copy and commonly supported
                    if hasattr(raw, 'color_desc'):
                        exif_dict["0th"][piexif.ImageIFD.Software] = "DNG Processor"
                    
                    # Add processing information
                    exif_dict["0th"][piexif.ImageIFD.ProcessingSoftware] = "rawpy + PIL"
                    
                except AttributeError:
                    pass
                
                # Generate EXIF bytes
                exif_bytes = piexif.dump(exif_dict)
                piexif.insert(exif_bytes, str(target_path))
                
        except Exception as e:
            print(f"Warning: Could not copy EXIF data for {source_path.name}: {e}")
            # Continue without EXIF data
    
    def resize_image(self, image: np.ndarray, scale_factor: int) -> np.ndarray:
        """Resize image by scale factor using high-quality resampling."""
        if scale_factor == 1:
            return image
        
        new_size = (image.shape[1] // scale_factor, image.shape[0] // scale_factor)
        
        resized = cv2.resize(image, new_size, interpolation=cv2.INTER_CUBIC)
        return resized
    
    def save_image_variants(self, image_array: np.ndarray, output_name: str) -> None:
        """Save image in different scales to appropriate subfolders."""
        scale_factors = [1, 2, 4]
        subfolder_names = ['images', 'images_2', 'images_4']
        
        for scale_factor, subfolder in zip(scale_factors, subfolder_names):
            scaled_array = self.resize_image(image_array, scale_factor)
            
            extension = 'jpg' if self.save_format == 'jpg' else 'png'
            output_path = self.output_folder / subfolder / f"{output_name}.{extension}"
            
            scaled_array = cv2.cvtColor(scaled_array, cv2.COLOR_RGB2BGR)
            if self.save_format == 'jpg':
                cv2.imwrite(str(output_path), scaled_array, [cv2.IMWRITE_JPEG_QUALITY, 95])
            else:
                cv2.imwrite(str(output_path), scaled_array)
    
    def process_single_file(self, dng_path: Path) -> bool:
        """Process a single DNG file."""
        try:
            rgb_array, metadata = self.process_raw_image(dng_path)

            output_name = dng_path.stem
            self.save_image_variants(rgb_array, output_name)
    
            if self.save_format == 'jpg':
                for subfolder in ['images', 'images_2', 'images_4']:
                    full_size_path = self.output_folder / subfolder / f"{output_name}.jpg"
                    self.copy_exif_data(dng_path, full_size_path)
            
            with self.lock:
                self.processed_count += 1
                print(f"Processed {self.processed_count}/{self.total_files}: {dng_path.name}")
            
            return True
            
        except Exception as e:
            print(f"Error processing {dng_path.name}: {e}")
            return False
    
    def process_all_files(self) -> None:
        """Process all DNG files using multithreading."""
        dng_files = self.get_dng_files()
        
        if not dng_files:
            print("No .dng files found in input folder.")
            return
        
        self.total_files = len(dng_files)
        print(f"Found {self.total_files} DNG files to process.")
        print(f"Using {self.max_workers} threads.")
        print(f"Color space: {self.color_space.upper()}")
        print(f"Output folder: {self.output_folder}")
        
        self.setup_output_folders()
        
        successful = 0
        failed = 0

        # for dng_file in dng_files:
        #     self.process_single_file(dng_file)
        # exit(0)
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_file = {
                executor.submit(self.process_single_file, dng_file): dng_file 
                for dng_file in dng_files
            }
            for future in as_completed(future_to_file):
                if future.result():
                    successful += 1
                else:
                    failed += 1
        
        print(f"\nProcessing complete!")
        print(f"Successfully processed: {successful}")
        print(f"Failed: {failed}")


def main():
    parser = argparse.ArgumentParser(description="Convert DNG files to JPG with multiple scales")
    parser.add_argument("input_folder", help="Input folder containing .dng files")
    parser.add_argument("output_folder", help="Output folder for processed images")
    parser.add_argument("--color-space", choices=['sRGB', 'ACES', 'linear-RGB', 'linear-ACES'], default='sRGB',
                       help="Color space for processing (default: sRGB)")
    parser.add_argument("--save-format", choices=['jpg', 'png', 'png16'], default='jpg',
                       help="Output image format (default: jpg)")
    parser.add_argument("--threads", type=int, default=None,
                       help="Number of processing threads (default: CPU count)")
    
    args = parser.parse_args()
    
    # Validate input folder
    if not Path(args.input_folder).exists():
        print(f"Error: Input folder '{args.input_folder}' does not exist.")
        sys.exit(1)
    
    # Create processor and run
    processor = DNGProcessor(
        input_folder=args.input_folder,
        output_folder=args.output_folder,
        color_space=args.color_space,
        save_format=args.save_format,
        max_workers=args.threads
    )
    
    try:
        processor.process_all_files()
    except KeyboardInterrupt:
        print("\nProcessing interrupted by user.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
