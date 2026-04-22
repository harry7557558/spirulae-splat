#!/usr/bin/env python3

"""Script to generate report for COLMAP dataset. (Distortion plot currently)"""

from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional, Union

import numpy as np
import matplotlib.pyplot as plt

from spirulae_splat.modules.colmap_utils import *

def generate_report(
    recon_dir: Path,
    save_path: Optional[Path] = None
):
    print("Recon dir:", recon_dir)
    if not recon_dir.exists() and recon_dir.is_dir():
        print(f"Error: Directory {recon_dir} does not exist.")
        exit(0)

    print("Loading points3D")
    colmap_points = load_colmap_points3D(recon_dir)
    print("Loading cameras")
    cam_id_to_camera = load_colmap_cameras(recon_dir)
    print("Loading images")
    im_id_to_image = load_colmap_images(recon_dir)

    camera_points = {}

    from tqdm import tqdm
    for im_id, im in tqdm(im_id_to_image.items(), "Processing"):

        camera = cam_id_to_camera[im.camera_id]
        camera_params = parse_colmap_camera_params(camera)

        (w, h, fx, fy, cx, cy, k1, k2, k3, k4, k5, k6, p1, p2, sx1, sy1) = \
            [camera_params.get(x, 0.0) for x in 'w h fl_x fl_y cx cy k1 k2 k3 k4 k5 k6 p1 p2 sx1 sy1'.split()]

        rotation = qvec2rotmat(im.qvec)
        translation = im.tvec.reshape(3, 1)

        xy_batch = im.xys[im.point3D_ids != -1]
        im_pids = im.point3D_ids[im.point3D_ids != -1]
        if len(im_pids) == 0:
            continue
        xyz_batch = np.stack([colmap_points[pid].xyz for pid in im_pids])    

        proj = xyz_batch @ rotation.T + translation.T

        x, y, z = proj.T
        if 'FISHEYE' in camera_params['camera_model']:
            r = np.hypot(x, y)
            theta = np.arctan2(r, z)
            x *= theta / (r + 1e-8)
            y *= theta / (r + 1e-8)
        else:
            x, y, z = x[z > 0], y[z > 0], z[z > 0]
            x, y = x / z, y / z
            xy_batch = xy_batch[np.where(z > 0)]

        r2 = x * x + y * y
        if camera_params['camera_model'] == "FULL_OPENCV":
            radial = (1+r2*(k1+r2*(k2+r2*k3))) / (1+r2*(k4+r2*(k5+r2*k6)))
        else:
            radial = 1.0 + r2*(k1 + r2*(k2 + r2*(k3 + r2*k4)))
        dx = 2.0*p1*x*y + p2*(r2+2.0*x*x) + sx1*r2
        dy = 2.0*p2*x*y + p1*(r2+2.0*y*y) + sy1*r2
        x = x * radial + dx
        y = y * radial + dy

        x = fx * x + cx
        y = fy * y + cy

        residual = np.stack((x, y), axis=-1) - xy_batch

        if im.camera_id not in camera_points:
            camera_points[im.camera_id] = {
                'gt_points': [],
                'residuals': [],
                'w': w,
                'h': h,
            }
        camera_points[im.camera_id]['gt_points'].append(xy_batch)
        camera_points[im.camera_id]['residuals'].append(residual)

    num_cameras = len(camera_points)
    print(num_cameras, "cameras")

    ncols = int(np.ceil(np.sqrt(num_cameras)))
    nrows = int(np.ceil(num_cameras / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(16, 12))
    if num_cameras == 1:
        axs = [axs]
    if nrows != 1 and nrows != 1:
        axs = sum([[*ax] for ax in axs], [])

    for axs_id, (camera_id, camera_point) in enumerate(sorted(camera_points.items())):
        print(f"Plotting camera {camera_id}")

        gt_points = camera_point['gt_points']
        residuals = camera_point['residuals']
        w, h = camera_point['w'], camera_point['h']

        gt_points = np.concatenate(gt_points)
        residuals = np.concatenate(residuals)

        sc = np.sqrt(w * h) / 64
        gw, gh = int(np.ceil(w / sc)), int(np.ceil(h / sc))
        # grid = np.zeros([gh, gw, 2])
        grid = np.zeros([gh, gw, 3])

        for (x, y), (dx, dy) in zip(gt_points, residuals):
            x, y = int(round(x / w * gw)), int(round(y / h * gh))
            if x < 0 or x >= gw or y < 0 or y >= gh:
                continue
            # if np.hypot(grid[y][x][0], grid[y][x][1]) < np.hypot(dx, dy):
            #     grid[y][x][0] = dx
            #     grid[y][x][1] = dy
            grid[y][x][0] += dx
            grid[y][x][1] += dy
            grid[y][x][2] += 1

        cnt = grid[:, :, 2].flatten()
        x, y = np.meshgrid(np.arange(gw), np.arange(gh))
        y = y.flatten()[cnt > 0] * (h / gh)
        x = x.flatten()[cnt > 0] * (w / gw)
        dy = grid[:, :, 1].flatten()[cnt > 0] / cnt[cnt > 0]
        dx = grid[:, :, 0].flatten()[cnt > 0] / cnt[cnt > 0]

        rms = (residuals[:, 0]**2 + residuals[:, 1]**2).mean()**0.5

        ax = axs[axs_id]
        ax.quiver(x, y, dx, dy, np.log10(np.hypot(dx, dy)+0.01),
                  headaxislength=0, headlength=0, headwidth=0)

        ax.set_aspect('equal')
        ax.set_xlim([0, w])
        ax.set_ylim([0, h])
        ax.invert_yaxis()
        ax.set_title(f"Camera {camera_id}, rms = {rms:.2f} px")

    plt.suptitle(f"Distortion plots for {recon_dir}")
    plt.tight_layout()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        print(f"Report saved to {save_path}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate report for COLMAP matching.")
    parser.add_argument("--dataset_dir", type=Path, required=True, help="Path to the dataset folder.")
    parser.add_argument("--colmap_model_path", type=Path, default="sparse/0", help="Path to the COLMAP reconstruction.")
    parser.add_argument("--save_report", help="If a path is specified, save report instead of display a window.")
    args = parser.parse_args()

    generate_report(
        (args.dataset_dir / args.colmap_model_path).absolute(),
        None if args.save_report is None else (args.dataset_dir / args.save_report).absolute()
    )
