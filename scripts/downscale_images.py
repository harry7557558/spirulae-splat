#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

import cv2


IMAGE_EXTENSIONS = {
    ".jpg",
    ".jpeg",
    ".png",
    ".bmp",
    ".tif",
    ".tiff",
    ".webp",
    ".jp2",
}


def is_image_file(path: Path) -> bool:
    return path.suffix.lower() in IMAGE_EXTENSIONS


def find_images(root: Path):
    return [
        p
        for p in root.rglob("*")
        if p.is_file() and is_image_file(p)
    ]


def save_image(dst_path: Path, image, jpeg_quality: int):
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    ext = dst_path.suffix.lower()

    if ext in (".jpg", ".jpeg"):
        cv2.imwrite(
            str(dst_path),
            image,
            [cv2.IMWRITE_JPEG_QUALITY, jpeg_quality],
        )
    else:
        cv2.imwrite(str(dst_path), image)


def process_image(
    src_path: Path,
    src_root: Path,
    dst_roots: dict[int, Path],
    scales: list[int],
    jpeg_quality: int,
):
    try:
        img = cv2.imread(
            str(src_path),
            cv2.IMREAD_UNCHANGED,
        )

        if img is None:
            print(f"Failed to read: {src_path}")
            return

        h, w = img.shape[:2]
        rel_path = src_path.relative_to(src_root)

        for scale in scales:
            new_w = max(1, w // scale)
            new_h = max(1, h // scale)

            resized = cv2.resize(
                img,
                (new_w, new_h),
                interpolation=cv2.INTER_AREA,
            )

            dst_path = dst_roots[scale] / rel_path

            save_image(
                dst_path,
                resized,
                jpeg_quality,
            )

    except Exception as e:
        print(f"Error processing {src_path}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate downscaled image pyramids."
    )

    parser.add_argument(
        "work_folder",
        type=Path,
    )

    parser.add_argument(
        "--images-dir",
        default="images",
    )

    parser.add_argument(
        "--scales",
        type=int,
        nargs="+",
        default=[2, 4, 8],
    )

    parser.add_argument(
        "--jpeg-quality",
        type=int,
        default=95,
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count(),
    )

    args = parser.parse_args()

    work_folder = args.work_folder.resolve()
    src_root = work_folder / args.images_dir

    if not src_root.is_dir():
        raise RuntimeError(
            f"Source folder not found: {src_root}"
        )

    images = find_images(src_root)

    print(f"Found {len(images)} images.")

    dst_roots = {
        scale: work_folder / f"{args.images_dir}_{scale}"
        for scale in args.scales
    }

    # Prevent nested threading
    cv2.setNumThreads(1)

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        for _ in tqdm(
            executor.map(
                lambda p: process_image(
                    p,
                    src_root,
                    dst_roots,
                    args.scales,
                    args.jpeg_quality,
                ),
                images,
            ),
            total=len(images),
            desc="Processing images",
        ):
            pass

    print("Done.")


if __name__ == "__main__":
    main()
