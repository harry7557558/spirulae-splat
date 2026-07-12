#!/usr/bin/env python3

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


def is_valid_mask(mask_path: Path) -> bool:
    if not mask_path.is_file():
        return False
    try:
        with Image.open(mask_path) as img:
            img.verify()
        # verify() doesn't decode pixel data; a truncated file can pass it
        with Image.open(mask_path) as img:
            img.load()
        return True
    except Exception:
        return False


def union_masks(outputs) -> Optional[np.ndarray]:
    """Union of all instance masks across outputs, as a 2D bool array.
    Returns None if nothing was detected."""
    result = None
    for output in outputs:
        masks = output['masks']
        if len(masks) == 0:
            continue
        if not isinstance(masks, np.ndarray):
            masks = masks.cpu().numpy()
        masks = np.any(masks, axis=tuple(range(masks.ndim - 2)))
        if not np.any(masks):
            continue
        result = masks if result is None else (result | masks)
    return result


def process(predictor, dataset_dir: str, image_dir: str, mask_dir: str):

    image_dir = Path(dataset_dir) / image_dir
    mask_dir = Path(dataset_dir) / mask_dir

    image_files = []
    for file_path in image_dir.rglob("*"):
        if file_path.is_file():
            file_path = file_path.relative_to(image_dir)
            image_files.append(file_path)

    for image_path in tqdm(sorted(image_files)):
        mask_path = mask_dir / (str(image_path) + ".png")
        if is_valid_mask(mask_path):
            continue
        try:
            image = Image.open(image_dir / image_path).convert("RGB")
        except:
            continue
        with redirect_stdout(StringIO()):
            pos_outputs, neg_outputs = predictor(image)
        pos_union = union_masks(pos_outputs)
        if pos_union is None:
            result_mask = np.ones((1, 1), dtype=np.bool_)
        else:
            neg_union = union_masks(neg_outputs)
            if neg_union is not None:
                # regions matching a negative prompt are kept even if they
                # also match a positive prompt
                pos_union &= ~neg_union
            result_mask = ~pos_union

        # resize mask to match original image resolution (model may run at reduced size)
        if result_mask.shape != (image.size[1], image.size[0]):
            mask_img = Image.fromarray(result_mask.astype(np.uint8) * 255)
            mask_img = mask_img.resize(image.size, Image.NEAREST)
            result_mask = np.asarray(mask_img) > 127

        os.makedirs(mask_path.parent, exist_ok=True)
        # write to a temp file then rename, so a crash mid-write can't leave a
        # corrupted file at the final path
        tmp_path = mask_path.with_name(mask_path.name + ".tmp")
        Image.fromarray(result_mask).save(tmp_path, format="PNG")
        os.replace(tmp_path, mask_path)

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
    parser.add_argument("input_dir", nargs=1, help="The dataset folder. Should contain a subfolder containing images. Masks will be saved to a different subfolder.")
    parser.add_argument("--prompt", required=True,
                        help="Text prompt for mask objects, semicolon separated. "
                            "Example: \"people; cars; shadow; fisheye border\"")
    parser.add_argument("--negative_prompt", "--negative-prompt", default="",
                        help="Text prompt for objects to exclude from masking, semicolon separated. "
                            "Regions matching these are kept even if they also match --prompt. "
                            "Example: --prompt \"person\" --negative_prompt \"person in painting\"")
    parser.add_argument("--images", default="images", help="Subfolder containing images. Default: images")
    parser.add_argument("--masks", default="masks", help="Subfolder to save masks. Default: masks")
    parser.add_argument("--max_image_size", type=int, default=1600, help="Maximum image size. Default: 1600")
    parser.add_argument("--model", default="sam2.1_hiera_large", help="SAM model to use.")
    parser.add_argument("--box_threshold", type=float, default=0.4, help="Box threshold for lang-sam model.")
    parser.add_argument("--text_threshold", type=float, default=0.25, help="Text threshold for lang-sam model.")
    args = parser.parse_args()

    prompts = [s.strip() for s in args.prompt.split(';') if s.strip() != ""]
    neg_prompts = [s.strip() for s in args.negative_prompt.split(';') if s.strip() != ""]

    def map_image(image: Image.Image):
        sc = args.max_image_size / max(image.size[0], image.size[1])
        if sc < 1.0:
            image = image.resize((int(image.size[0]*sc), int(image.size[1]*sc)))
        return image

    # lang-sam
    if args.model != "sam3":
        try:
            from lang_sam.lang_sam import LangSAM
        except ImportError:
            print("lang-sam not found or not installed properly. Please install https://github.com/luca-medeiros/lang-segment-anything")
            exit(0)
        model = LangSAM(args.model, device="cuda")

        all_prompts = prompts + neg_prompts

        def predict(image_pil: Image.Image):
            images_pil = [map_image(image_pil)]
            results = model.predict(images_pil*len(all_prompts), all_prompts, args.box_threshold, args.text_threshold)
            return results[:len(prompts)], results[len(prompts):]

    # SAM-3 (better in quality, need to request access)
    else:
        try:
            from sam3.model_builder import build_sam3_image_model
            from sam3.model.sam3_image_processor import Sam3Processor
        except ImportError:
            print("SAM-3 not found or not installed properly. Please install https://github.com/facebookresearch/sam3.git")
            exit(0)
        import torch
        model = build_sam3_image_model()
        processor = Sam3Processor(model)

        # text features are independent of the image; compute them once for
        # the whole dataset instead of per image
        with torch.inference_mode():
            text_features = {
                prompt: model.backbone.forward_text([prompt], device=processor.device)
                for prompt in prompts + neg_prompts
            }

        def predict_prompt(state, prompt):
            state["backbone_out"].update(text_features[prompt])
            if "geometric_prompt" not in state:
                state["geometric_prompt"] = model._get_dummy_prompt()
            processor._forward_grounding(state)
            # _forward_grounding overwrites state["masks"] on each call
            return {"masks": state["masks"]}

        def predict(image: Image.Image):
            image = map_image(image)
            # image embedding is independent of the prompt; run the backbone
            # once and reuse it for all prompts
            state = processor.set_image(image)
            pos = [predict_prompt(state, prompt) for prompt in prompts]
            neg = [predict_prompt(state, prompt) for prompt in neg_prompts]
            return pos, neg

    process(
        predict,
        args.input_dir[0],
        args.images,
        args.masks,
    )
