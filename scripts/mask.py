import numpy as np
from PIL import Image

from dataclasses import dataclass
from typing import Optional

import os
from tqdm import tqdm
import json

from PIL import Image

from pathlib import Path

from io import StringIO
from contextlib import redirect_stdout


def process(predictor, dataset_dir: str, image_dir: str, mask_dir: str):

    image_dir = Path(dataset_dir) / image_dir
    mask_dir = Path(dataset_dir) / mask_dir

    image_files = []
    for file_path in image_dir.rglob("*"):
        if file_path.is_file():
            file_path = file_path.relative_to(image_dir)
            image_files.append(file_path)

    for image_path in tqdm(sorted(image_files)):
        try:
            image = Image.open(image_dir / image_path).convert("RGB")
        except:
            continue
        with redirect_stdout(StringIO()):
            outputs = predictor([image])
        result_mask = None
        for output in outputs:
            labels = output['labels']
            # print(labels)
            masks = output['masks']
            if not np.any(masks):
                continue
            masks = ~np.any(masks, axis=0)
            if result_mask is None:
                result_mask = masks
            else:
                result_mask &= masks
        if result_mask is None:
            result_mask = np.ones((1, 1), dtype=np.bool_)

        mask_path = mask_dir / (str(image_path) + ".png")
        os.makedirs(mask_path.parent, exist_ok=True)
        Image.fromarray(result_mask).save(mask_path)

        if False:
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(image)
            ax2.imshow(result_mask)
            plt.show()
            # break


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate masks with natural language prompts.")
    parser.add_argument("input_dir", nargs=1, help="The image directory.")
    parser.add_argument("--prompt", required=True,
                        help="Text prompt for mask objects, semicolon separated. "
                            "Example: \"people; cars; shadow; fisheye border\"")
    parser.add_argument("--images", default="images", help="Name of the subdirectory containing images.")
    parser.add_argument("--masks", default="masks", help="Name of the subdirectory to save masks.")
    parser.add_argument("--max_image_size", type=int, default=1600, help="Maximum image size.")
    parser.add_argument("--sam_type", default="sam2.1_hiera_small", help="SAM model to use.")
    parser.add_argument("--box_threshold", type=float, default=0.4, help="Box threshold for lang-sam model.")
    parser.add_argument("--text_threshold", type=float, default=0.25, help="Text threshold for lang-sam model.")
    args = parser.parse_args()

    try:
        from lang_sam.lang_sam import LangSAM
    except ImportError:
        print("lang-sam not found. Please install https://github.com/luca-medeiros/lang-segment-anything")
        exit(0)
    model = LangSAM(args.sam_type, device="cuda")

    def map_image(image: Image.Image):
        sc = args.max_image_size / max(image.size[0], image.size[1])
        if sc < 1.0:
            image = image.resize((int(image.size[0]*sc), int(image.size[1]*sc)))
        return image

    prompts = [s.strip() for s in args.prompt.split(';') if s.strip() != ""]
    # prompts = ["fisheye circle"]
    def predict(images_pil: list[Image.Image]):
        images_pil = [map_image(im) for im in images_pil]

        results = model.predict(images_pil*len(prompts), prompts, args.box_threshold, args.text_threshold)
        return results

        results = []
        for prompt in prompts:
            pred = model.predict(images_pil, [prompt], args.box_threshold, args.text_threshold)
            results.extend(pred)
        print(results)
        return results

    process(
        predict,
        args.input_dir[0],
        args.images,
        args.masks,
    )
