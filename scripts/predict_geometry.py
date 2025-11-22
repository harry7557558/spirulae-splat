#!/usr/bin/env python3

import torch
import numpy as np
from PIL import Image

from dataclasses import dataclass
from typing import Optional
from pathlib import Path

import os
from tqdm import tqdm
import json

from spirulae_splat.splat.cuda import undistort_image, distort_image

from typing import Literal

from io import StringIO
from contextlib import redirect_stdout


@dataclass
class Config:
    dataset_dir: str
    max_size: int = 1600  # result not always better at high res


def process_image(
    model_type: Literal['metric3d', 'da3'],
    model,
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

    image = torch.from_numpy(np.array(image))[None].float().cuda() / 255.0

    sw = image.shape[2] / intrins['w']
    sh = image.shape[1] / intrins['h']
    Ks = (intrins['fl_x']*sw, intrins['fl_y']*sh, intrins['cx']*sw, intrins['cy']*sh)
    dist_coeffs = tuple(intrins.get(key, 0.0) for key in "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split())
    intrins = (intrins['camera_model'], Ks, dist_coeffs)
    image = undistort_image(image, *intrins)

    if model_type == "metric3d":
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            pred_depth, _, output_dict = model.inference({'input': image.permute((0, 3, 1, 2))})
        pred_normal = output_dict['prediction_normal']
    elif model_type == "da3":
        # TODO: provide intrinsics, or verify the model doesn't use intrinsics
        with redirect_stdout(StringIO()):
            outputs = model.inference([Image.fromarray((255*image[0].clip(0,1)).to('cpu',torch.uint8).numpy())])
        depth, sky = outputs.depth[0], outputs.sky[0]
        depth = np.where(sky, np.amax(depth)*np.ones_like(depth), depth)
        pred_depth = torch.from_numpy(depth)[None][None].cuda()

    if depth_save_path is not None:
        pred_depth = torch.nn.functional.interpolate(
            pred_depth, size=(h, w), mode='bilinear', align_corners=False
        ).permute(0, 2, 3, 1)  # [h, w, 1]
        pad = Config.max_size // 200 + 1  # fix bad values in case undistortion gives black border that messes up model
        pred_depth /= torch.amax(pred_depth[0, pad:-pad, pad:-pad])
        pred_depth = distort_image(pred_depth, *intrins)[0]
        pred_depth = torch.clip(65535*pred_depth, 0, 65535).cpu().numpy().astype(np.uint16)
        os.makedirs(Path(depth_save_path).parent, exist_ok=True)
        Image.fromarray(pred_depth.squeeze(-1), mode='I;16').save(depth_save_path)

    if normal_save_path is not None:
        pred_normal = torch.nn.functional.interpolate(
            pred_normal, size=(h, w), mode='bilinear', align_corners=False
        )[:, :3].permute(0, 2, 3, 1)  # [h, w, 3]
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


def process_dir(dataset_dir: str, include_normal: bool):
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
        elif 'normal_file_path' in frame:
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

    print("Loading Depth Anything 3 model...")
    try:
        from depth_anything_3.api import DepthAnything3
    except ImportError:
        print("Depth Anything v3 not found. Please install following https://huggingface.co/depth-anything/DA3MONO-LARGE.")
        exit(0)
    model = DepthAnything3.from_pretrained("depth-anything/da3mono-large")
    model = model.cuda()
    print("Depth Anything 3 model loaded")

    for (image_filename, depth_filename, normal_filename, intrins) \
          in tqdm(image_filenames, "Predicting depths"):
        process_image(
            "da3", model, intrins,
            os.path.join(image_dir, image_filename),
            os.path.join(depth_dir, depth_filename),
            # os.path.join(normal_dir, normal_filename)
            None
        )

    if not include_normal:
        return

    print("Loading Metric3D v2 model...")
    model = torch.hub.load('yvanyin/metric3d', 'metric3d_vit_large', pretrain=True)
    model = model.eval().cuda().bfloat16()
    print("Metric3D v2 model loaded")

    for (image_filename, depth_filename, normal_filename, intrins) \
          in tqdm(image_filenames, "Predicting normals"):
        process_image(
            "metric3d", model, intrins,
            os.path.join(image_dir, image_filename),
            # os.path.join(depth_dir, depth_filename),
            None,
            os.path.join(normal_dir, normal_filename)
        )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Export PLY (and equirectangular map), for 3DGS only.")
    parser.add_argument("dataset_dir", nargs=1, help="Path to the dataset folder.")
    parser.add_argument("--normal", action="store_true", help="Whether to predict normal in addition to depth.")
    args = parser.parse_args()
    process_dir(args.dataset_dir[0], args.normal)
