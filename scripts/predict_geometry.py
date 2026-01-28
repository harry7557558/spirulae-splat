#!/usr/bin/env python3

import torch
import torch.nn.functional as F
import numpy as np
from PIL import Image

from dataclasses import dataclass
from typing import Optional
from pathlib import Path

import os
from tqdm import tqdm
import json

from spirulae_splat.splat.cuda import undistort_image, distort_image

from typing import Tuple, Literal

from io import StringIO
from contextlib import redirect_stdout


@dataclass
class Config:
    dataset_dir: str
    max_size: int = 1600  # result not always better at high res


def expand_white_area(boolean_image, offset):
    """
    Expands the white area in a boolean image by a given offset using dilation via max pooling.

    Args:
        boolean_image (torch.Tensor): Input boolean image tensor (C, H, W or B, C, H, W).
                                      White pixels should be True/1 and black pixels False/0.
        offset (int): The number of pixels to expand the white areas by.
                      This will result in a kernel size of (2 * offset + 1).

    Returns:
        torch.Tensor: The image with expanded white areas.
    """
    if offset <= 0:
        return boolean_image
    
    # Ensure the input is a float tensor with a batch dimension (B, C, H, W) for max_pool2d
    if boolean_image.dim() == 3:
        image_with_batch = boolean_image.unsqueeze(0)
    else:
        image_with_batch = boolean_image.clone()

    # Convert to float for pooling operation
    image_float = image_with_batch.float()

    # Kernel size for dilation is (2 * offset + 1)
    kernel_size = 2 * offset + 1
    
    # Use max pooling with appropriate padding to implement dilation
    # A pixel becomes "white" if any pixel in its neighborhood (defined by the kernel) was white
    expanded_image_float = F.max_pool2d(
        image_float,
        kernel_size=kernel_size,
        stride=1,
        padding=offset
    )

    # Convert the result back to the original boolean type
    # Values > 0 will become True/1, and 0 will become False/0
    expanded_image_bool = expanded_image_float > 0.5
    
    if boolean_image.dim() == 3:
        return expanded_image_bool.squeeze(0)
    else:
        return expanded_image_bool


def process_image(
    model_type: Literal['metric3d', 'da3', 'da3+lang-sam', 'metric3d+da3+lang-sam'],
    models: Tuple,
    intrins: dict,
    image_path: str,
    depth_save_path: Optional[str]=None,
    normal_save_path: Optional[str]=None
):
    if (depth_save_path is None or os.path.exists(depth_save_path)) and \
        (normal_save_path is None or os.path.exists(normal_save_path)):
        return

    image = Image.open(image_path).convert("RGB")
    w, h = image.size
    sc = Config.max_size / max(w, h)
    if sc < 1.0:
        w, h = int(sc*w+0.5), int(sc*h+0.5)
        image = image.resize((w, h))

    image_original = image
    image = torch.from_numpy(np.array(image))[None].float().cuda() / 255.0

    sw = image.shape[2] / intrins['w']
    sh = image.shape[1] / intrins['h']
    intrins = (intrins['fl_x']*sw, intrins['fl_y']*sh, intrins['cx']*sw, intrins['cy']*sh)
    dist_coeffs = tuple(intrins.get(key, 0.0) for key in "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split())
    intrins = (intrins['camera_model'], intrins, dist_coeffs)
    image = undistort_image(image, *intrins)

    # pred_depth and pred_normal are distorted, pred_sky is original
    pred_depth, pred_normal, pred_sky = None, None, None

    if "metric3d" in model_type:
        model = models[model_type.split('+').index("metric3d")]
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            pred_depth, _, output_dict = model.inference({'input': image.permute((0, 3, 1, 2))})
        pred_normal = output_dict['prediction_normal']

        pred_depth = torch.nn.functional.interpolate(
            pred_depth, size=(h, w), mode='bilinear', align_corners=False
        ).permute(0, 2, 3, 1)  # [h, w, 1]
        pred_normal = torch.nn.functional.interpolate(
            pred_normal, size=(h, w), mode='bilinear', align_corners=False
        )[:, :3].permute(0, 2, 3, 1)  # [h, w, 3]

        pad = Config.max_size // 200 + 1  # fix bad values in case undistortion gives black border that messes up model
        sky_depth = torch.amax(pred_depth[0, pad:-pad, pad:-pad]).item()
        pred_depth = pred_depth.clip(max=sky_depth)

    # for fine sky only (depth is worse than metric3d)
    if "da3" in model_type:
        model = models[model_type.split('+').index("da3")]
        with redirect_stdout(StringIO()):
            outputs = model.inference([
                # image_original,
                Image.fromarray((255*image[0].clip(0,1)).to('cpu',torch.uint8).numpy())
            ])
        depth, sky = outputs.depth[-1], outputs.sky[-1]
        # sky_original = outputs.sky[0]  # doesn't work well with fisheye circle
        sky = torch.from_numpy(sky)[None][None].cuda()
        sky = torch.nn.functional.interpolate(
            sky.float(), size=(h, w), mode='bilinear', align_corners=False
        ).permute(0, 2, 3, 1) > 0.5  # [h, w, 1]

        assert pred_depth is not None
        sky_diffused = expand_white_area(sky, pad)  # TODO: might hurt distant details
        sky_depth = torch.quantile(pred_depth[~sky_diffused], 0.9999).item()
        pred_depth = torch.where(sky, sky_depth*torch.ones_like(pred_depth), pred_depth.clip(max=sky_depth))

    # for rough sky (can be inaccurate / miss high frequency details)
    # TODO: SAM-3
    if "lang-sam" in model_type:
        model = models[model_type.split('+').index("lang-sam")]
        box_threshold = 0.3
        text_threshold = 0.25
        with redirect_stdout(StringIO()):
            results = model.predict([image_original], ["sky"], box_threshold, text_threshold)
        for output in results:
            masks = output['masks']
            if not np.any(masks):
                continue
            masks = np.any(masks, axis=0)
            if pred_sky is None:
                pred_sky = masks
            else:
                pred_sky |= masks
        if pred_sky is not None:
            pred_sky = torch.from_numpy(pred_sky)[None][None].cuda()
            pred_sky = torch.nn.functional.interpolate(
                pred_sky.float(), size=(h, w), mode='bilinear', align_corners=False
            ).permute(0, 2, 3, 1)[0] > 0.5  # [h, w, 1]

    if depth_save_path is not None:
        pred_depth /= sky_depth
        pred_depth = distort_image(pred_depth, *intrins)[0]
        if pred_sky is not None:
            pred_depth = torch.where(
                pred_sky & (pred_depth == 0.0),
                torch.ones_like(pred_depth),
                pred_depth
            )
        pred_depth = torch.clip(65535*pred_depth, 0, 65535).cpu().numpy().astype(np.uint16)
        os.makedirs(Path(depth_save_path).parent, exist_ok=True)
        Image.fromarray(pred_depth.squeeze(-1), mode='I;16').save(depth_save_path)

    if normal_save_path is not None:
        # TODO: probably also need to distort values of normals
        pred_normal = 0.5 + 0.5 * pred_normal / torch.norm(pred_normal, dim=-1, keepdim=True)
        pred_normal = distort_image(pred_normal, *intrins)[0]
        pred_normal = torch.clip(255*pred_normal, 0, 255).to(torch.uint8).cpu().numpy()
        os.makedirs(Path(normal_save_path).parent, exist_ok=True)
        Image.fromarray(pred_normal).save(normal_save_path)


def update_intrins(in_obj, out_obj):
    out_obj = {**out_obj}

    if 'camera_model' in in_obj:
        out_obj['camera_model'] = {
            'SIMPLE_PINHOLE': "pinhole",
            'PINHOLE': "pinhole",
            'SIMPLE_RADIAL': "pinhole",
            'SIMPLE_RADIAL_FISHEYE': "fisheye",
            'RADIAL': "pinhole",
            'RADIAL_FISHEYE': "fisheye",
            'OPENCV': "pinhole",
            'OPENCV_FISHEYE': "fisheye",
            'THIN_PRISM_FISHEYE': "fisheye",
        }[in_obj['camera_model']]
    elif 'camera_model' not in out_obj:
        out_obj['camera_model'] = 'pinhole'

    for key in 'w h fl_x fl_y cx cy k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2'.split():
        if key in in_obj:
            out_obj[key] = in_obj[key]
        elif key not in out_obj:
            out_obj[key] = 0.0

    return out_obj


def process_dir(dataset_dir: str, include_normal: bool, include_sky: bool):
    image_dir = os.path.join(dataset_dir, "images")

    def with_ext(file_path, new_ext):
        base_name, _ = os.path.splitext(file_path)
        return base_name + '.' + new_ext

    image_filenames = []
    with open(os.path.join(dataset_dir, "transforms.json")) as fp:
        content = json.load(fp)
        global_intrins = update_intrins(content, {})
    for frame in content['frames']:
        file_path = frame['file_path'].lstrip('./')
        assert file_path.startswith("images")
        image_filename = file_path[len("images")+1:]
        depth_filename = with_ext(image_filename, 'png')
        normal_filename = with_ext(image_filename, 'png')
        image_filenames.append((
            image_filename, depth_filename, normal_filename,
            update_intrins(frame, global_intrins)
        ))
        frame["depth_file_path"] = os.path.join("depths", depth_filename)
        if include_normal:
            frame["normal_file_path"] = os.path.join("normals", normal_filename)
        elif 'normal_file_path' in frame and False:
            del frame['normal_file_path']
    if len(image_filenames) == 0:
        print("No image found")
        return
    image_filenames.sort()
    with open(os.path.join(dataset_dir, "transforms.json"), 'w') as fp:
        json.dump(content, fp, indent=4)

    depth_dir = os.path.join(dataset_dir, "depths")
    normal_dir = os.path.join(dataset_dir, "normals")
    os.makedirs(depth_dir, exist_ok=True)
    if include_normal:
        os.makedirs(normal_dir, exist_ok=True)

    # Metric3D v2 for good depth/normal
    print("Loading Metric3D v2 model...")
    try:
        model_metric3d = torch.hub.load('yvanyin/metric3d', 'metric3d_vit_large', pretrain=True)
    except ImportError:
        print("mmengine not found. Please run `pip install mmcv`.")
        exit(0)
    model_metric3d = model_metric3d.eval().cuda().bfloat16()
    print("Metric3D v2 model loaded")

    if include_sky:
        # DA2 for good sky
        print("Loading Depth Anything 3 model...")
        try:
            from depth_anything_3.api import DepthAnything3
        except ImportError:
            print("Depth Anything v3 not found. Please install following https://huggingface.co/depth-anything/DA3MONO-LARGE.")
            exit(0)
        model_da3 = DepthAnything3.from_pretrained("depth-anything/da3mono-large").cuda()
        print("Depth Anything 3 model loaded")

        # lang-sam for rough sky; TODO: SAM-3
        print("Loading lang-sam model...")
        try:
            from lang_sam.lang_sam import LangSAM
        except ImportError:
            print("lang-sam not found. Please install https://github.com/luca-medeiros/lang-segment-anything.")
            exit(0)
        model_langsam = LangSAM("sam2.1_hiera_large", device="cuda")
        models = (model_metric3d, model_da3, model_langsam)
        print("lang-sam model loaded")

    else:
        models = (model_metric3d,)

    for (image_filename, depth_filename, normal_filename, intrins) \
          in tqdm(image_filenames, "Predicting depths"):
        process_image(
            "metric3d+da3+lang-sam" if include_sky else "metric3d",
            models, intrins,
            os.path.join(image_dir, image_filename),
            os.path.join(depth_dir, depth_filename),
            os.path.join(normal_dir, normal_filename),
        )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Export PLY (and equirectangular map), for 3DGS only.")
    parser.add_argument("dataset_dir", nargs=1, help="Path to the dataset folder.")
    # parser.add_argument("--normal", action="store_true", help="Whether to predict normal in addition to depth.")
    parser.add_argument("--sky", action="store_true", help="Whether to predict sky for full images. Useful for highly distorted fisheye images.")
    args = parser.parse_args()
    # process_dir(args.dataset_dir[0], args.normal, args.sky)
    process_dir(args.dataset_dir[0], True, args.sky)
