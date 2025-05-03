#!/usr/bin/env python3

# Script for coloring point clouds for nerfstudio, especially those generated with GLOMAP

# To use: rename sparse_pc.ply to sparse_pc_raw.ply,
# then run this script in the same directory as nerfstudio dataset

import json
import numpy as np
import torch
import cv2
from plyfile import PlyData, PlyElement
from tqdm import tqdm
import os
import sys

def load_transforms(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data

def load_ply(ply_path):
    plydata = PlyData.read(ply_path)
    vertices = plydata['vertex']
    return torch.tensor(np.vstack([vertices['x'], vertices['y'], vertices['z']]).T, dtype=torch.float32).cuda()

def get_intrinsics(transforms):
    w, h = transforms['w'], transforms['h']
    K = np.array([
        [transforms['fl_x'], 0, transforms['cx']],
        [0, transforms['fl_y'], transforms['cy']],
        [0, 0, 1]
    ], dtype=np.float32)
    model = transforms['camera_model']
    dist_coeffs = [0.0, 0.0, 0.0, 0.0]
    if model == "OPENCV":
        dist_coeffs = [transforms.get(k, 0.0) for k in ['k1', 'k2', 'p1', 'p2']]
    elif model == "OPENCV_FISHEYE":
        dist_coeffs = [transforms.get(k, 0.0) for k in ['k1', 'k2', 'k3', 'k4']]
    dist_coeffs = np.array(dist_coeffs, dtype=np.float32)
    return w, h, K, model, dist_coeffs

def main(workdir):
    transforms_path = os.path.join(workdir, "transforms.json")
    input_ply_path = os.path.join(workdir, "sparse_pc_no_color.ply")
    output_ply_path = os.path.join(workdir, "sparse_pc.ply")
    # copy ply file
    open(input_ply_path, "wb").write(open(output_ply_path, "rb").read())

    transforms = load_transforms(transforms_path)
    points = load_ply(input_ply_path)
    
    colors_val = torch.zeros((points.shape[0], 3), dtype=torch.float32).cuda()
    colors_depth = 1e10 * torch.ones((points.shape[0], 1), dtype=torch.float32).cuda()

    wg, hg, Kg, modelG, distCoeffsG = [None]*5
    if 'fl_x' in transforms:
        wg, hg, Kg, modelG, distCoeffsG = get_intrinsics(transforms)

    b2w = transforms['applied_transform']
    if len(b2w) < 4:
        b2w.append([0.0, 0.0, 0.0, 1.0])
    b2w = np.array(b2w)

    transforms['frames'].sort(key=lambda _: _['file_path'])

    num_points = points.shape[0]
    num_images = len(transforms['frames'])

    for frame in tqdm(transforms['frames']):
        image = cv2.imread(os.path.join(workdir, frame['file_path']))
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        image_tensor = torch.from_numpy(image).cuda().float()

        w, h, K, model, distCoeffs = wg, hg, Kg, modelG, distCoeffsG
        if 'fl_x' in frame:
            w, h, K, model, distCoeffs = get_intrinsics(frame)
        c2w = frame['transform_matrix']
        w2c = np.linalg.inv(c2w)
        b2c = w2c #@ b2w
        b2c = np.diag([1, -1, -1, 1]) @ b2c
        b2c = torch.from_numpy(b2c).float().cuda()
        R = b2c[:3, :3]
        t = b2c[:3, 3:4]

        # projection
        pts_cam = points @ R.T + t.T
        depths = pts_cam[:, 2]
        x = pts_cam[:, 0] / depths
        y = pts_cam[:, 1] / depths

        if model == "OPENCV":
            r2 = x*x+y*y
            k1, k2, p1, p2 = distCoeffs
            radial = 1 + r2*(k1 + k2*r2)
            dx = 2*p1*x*y + p2*(r2 + 2*x*x)
            dy = p1*(r2 + 2*y*y) + 2*p2*x*y
            x = x * radial + dx
            y = y * radial + dy

        elif model == "OPENCV_FISHEYE":
            r = torch.hypot(x, y)
            theta = torch.atan(r)
            theta2 = theta**2
            k1, k2, k3, k4 = distCoeffs
            theta_d = theta * (1+theta2*(k1+theta2*(k2+theta2*(k3+theta2*k4))))
            scale = theta_d / torch.clip(r, min=1e-8)
            x = x * scale
            y = y * scale

        fx, fy, cx, cy = K[0,0], K[1,1], K[0,2], K[1,2]
        pts_2d_f = torch.stack([fx*x+cx, fy*y+cy], dim=1)
        pts_2d = (pts_2d_f+0.5).int()

        # filter out of bounds
        valid_mask = (pts_2d[:, 0] >= 0) & (pts_2d[:, 0] < w) & \
                     (pts_2d[:, 1] >= 0) & (pts_2d[:, 1] < h) & \
                     (depths > 0)

        # filter occluded - not perfect but good enough for initializing splats
        if True:
            # downscale
            pts_2d_filtered = pts_2d[valid_mask]
            # sc = min(0.1 * np.sqrt(len(pts_2d_filtered)/(w*h)), 1)
            sc = min(1.0 * np.sqrt((num_points/num_images)/(w*h)), 1)
            hs, ws = int(sc*h+1), int(sc*w+1)
            pts_2d_s = (pts_2d_f * sc + 0.5).long()
            valid_mask_s = (pts_2d_s[:, 0] >= 0) & (pts_2d_s[:, 0] < ws) & \
                        (pts_2d_s[:, 1] >= 0) & (pts_2d_s[:, 1] < hs) & \
                        (depths > 0)
            pts_2d_s = pts_2d_s[valid_mask_s]

            # render a flat depth image
            image = torch.zeros(hs*ws, device="cuda") + torch.amax(depths)
            raster_indices = pts_2d_s[:,1] * ws + pts_2d_s[:,0]
            raster_depths = depths[valid_mask_s]
            image.scatter_reduce_(0, raster_indices, raster_depths, reduce='amin', include_self=True)
            if False:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.imshow(image.reshape((hs, ws)).cpu().numpy())
                plt.show()

            # depth testing
            is_visible = 1.1 * image[raster_indices] > raster_depths
            valid_mask[valid_mask_s] &= is_visible

        # apply mask
        if 'mask_path' in frame:
            mask = cv2.imread(os.path.join(workdir, frame['mask_path']))[:,:,0] > 0
            mask_tensor = torch.from_numpy(mask).cuda()
            pts_2d_filtered = pts_2d[valid_mask]
            valid_mask[valid_mask.clone()] &= mask_tensor[pts_2d_filtered[:, 1], pts_2d_filtered[:, 0]]

        # accumulate color
        valid_pts_2d = pts_2d[valid_mask]
        colors = image_tensor[valid_pts_2d[:, 1], valid_pts_2d[:, 0]]

        colors_val[valid_mask] = torch.where(
            (pts_cam[:, 2:] < colors_depth)[valid_mask],
            colors, colors_val[valid_mask])
        colors_depth[valid_mask] = torch.fmin(colors_depth, pts_cam[:, 2:])[valid_mask]

    mean_colors = colors_val.clamp(0, 255).byte().cpu().numpy()
    
    vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                    ('red', 'u1'), ('green', 'u1'), ('blue', 'u1')]
    vertices = np.empty(points.shape[0], dtype=vertex_dtype)
    vertices['x'] = points.cpu().numpy()[:, 0]
    vertices['y'] = points.cpu().numpy()[:, 1]
    vertices['z'] = points.cpu().numpy()[:, 2]
    vertices['red'] = mean_colors[:, 0]
    vertices['green'] = mean_colors[:, 1]
    vertices['blue'] = mean_colors[:, 2]
    
    el = PlyElement.describe(vertices, 'vertex')
    PlyData([el]).write(output_ply_path)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 path/to/glomap_color.py path/to/colmap/directory")
        exit(-1)

    workdir = sys.argv[1]
    main(workdir)

