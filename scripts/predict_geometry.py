import torch
import numpy as np
from PIL import Image

from dataclasses import dataclass
from typing import Optional

import os
from tqdm import tqdm
import json


@dataclass
class Config:
    dataset_dir: str
    max_size: int = 99999999


def process_image(
    model,
    image_path: str,
    depth_save_path: Optional[str]=None,
    normal_save_path: Optional[str]=None
):
    image = Image.open(image_path).convert("RGB")
    w, h = image.size
    sc = Config.max_size / max(w, h)
    if sc < 1.0:
        w, h = int(sc*w+0.5), int(sc*h+0.5)
        image = image.resize((w, h))

    image = torch.from_numpy(np.array(image)).permute((2, 0, 1))[None].float().cuda() / 255.0
    with torch.autocast(device_type="cuda", dtype=torch.float32):
        pred_depth, _, output_dict = model.inference({'input': image})
    pred_normal = output_dict['prediction_normal']

    if depth_save_path is not None:
        pred_depth = torch.nn.functional.interpolate(
            pred_depth, size=(h, w), mode='bilinear', align_corners=False
        )[0].permute(1, 2, 0)  # [h, w, 1]
        # pred_depth = torch.log(torch.relu(pred_depth) + 1.0)
        pred_depth /= torch.amax(pred_depth)
        pred_depth = torch.clip(65535*pred_depth, 0, 65535).cpu().numpy().astype(np.uint16)
        Image.fromarray(pred_depth.squeeze(-1), mode='I;16').save(depth_save_path)

    if normal_save_path is not None:
        pred_normal = torch.nn.functional.interpolate(
            pred_normal, size=(h, w), mode='bilinear', align_corners=False
        )[0, :3].permute(1, 2, 0)  # [h, w, 3]
        pred_normal = 0.5 + 0.5 * pred_normal / torch.norm(pred_normal, dim=-1, keepdim=True)
        pred_normal = torch.clip(255*pred_normal, 0, 255).to(torch.uint8).cpu().numpy()
        Image.fromarray(pred_normal).save(normal_save_path)

def process_dir(dataset_dir: str):
    image_dir = os.path.join(dataset_dir, "images")

    def with_ext(file_path, new_ext):
        base_name, _ = os.path.splitext(file_path)
        return base_name + '.' + new_ext

    image_filenames = []
    with open(os.path.join(dataset_dir, "transforms.json")) as fp:
        content = json.load(fp)
    for frame in content['frames']:
        file_path = frame['file_path'].lstrip('./')
        assert file_path.startswith("images")
        image_filename = file_path[len("images")+1:]
        depth_filename = with_ext(image_filename, 'png')
        normal_filename = with_ext(image_filename, 'png')
        image_filenames.append((image_filename, depth_filename, normal_filename))
        frame["depth_file_path"] = os.path.join("depths", depth_filename)
        frame["normal_file_path"] = os.path.join("normals", normal_filename)
    if len(image_filenames) == 0:
        print("No image found")
        return
    image_filenames.sort()
    with open(os.path.join(dataset_dir, "transforms.json"), 'w') as fp:
        json.dump(content, fp, indent=4)

    depth_dir = os.path.join(dataset_dir, "depths")
    normal_dir = os.path.join(dataset_dir, "normals")
    os.makedirs(depth_dir, exist_ok=True)
    os.makedirs(normal_dir, exist_ok=True)

    print("Loading model...")
    model = torch.hub.load('yvanyin/metric3d', 'metric3d_vit_large', pretrain=True)
    model = model.eval().cuda()
    print("Model loaded")

    for (image_filename, depth_filename, normal_filename) \
          in tqdm(image_filenames, "Predicting depth/normal"):
        process_image(
            model,
            os.path.join(image_dir, image_filename),
            os.path.join(depth_dir, depth_filename),
            os.path.join(normal_dir, normal_filename)
        )


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python3 predict_geometry.py path/to/dataset/dir")
    process_dir(sys.argv[1])
