# Copyright 2022 the Regents of the University of California, Nerfstudio Team and contributors. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Modified from https://docs.nerf.studio/_modules/nerfstudio/data/dataparsers/nerfstudio_dataparser.html."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional, Tuple, Type

import numpy as np
import torch

import json

from spirulae_splat.modules.camera import (
    colmap_camera_model_to_type,
    Cameras,
)

from nerfstudio.cameras import camera_utils
from nerfstudio.data.utils.dataparsers_utils import (
    get_train_eval_split_all,
    get_train_eval_split_filename,
    get_train_eval_split_fraction,
    get_train_eval_split_interval,
)

MAX_AUTO_RESOLUTION = 1600


DISTORTION_KEYS = "k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2".split()


def geometric_median(X, eps=0.0, maxiter=10):
    """
    Calculates the geometric median of a set of 3D points.

    Args:
        X (np.ndarray): A NumPy array of shape (N, 3) where N is the number of points.
        eps (float): Tolerance for convergence.

    Returns:
        np.ndarray: The geometric median point (shape 3,).
    """
    from scipy.spatial.distance import cdist, euclidean
    y = np.median(X, 0)
    for iter in range(maxiter):
        D = cdist(X, [y])
        nonzeros = (D != 0)[:, 0]
        Dinv = 1 / D[nonzeros]
        Dinvs = np.sum(Dinv)
        W = Dinv / Dinvs
        T = np.sum(W * X[nonzeros], 0)
        num_zeros = len(X) - np.sum(nonzeros)

        if num_zeros == 0:
            y1 = T
        elif num_zeros == len(X):
            return y
        else:
            R = (T - y) * Dinvs
            r = np.linalg.norm(R)
            rinv = 0 if r == 0 else num_zeros / r
            y1 = max(0, 1 - rinv) * T + min(1, rinv) * y
        
        if euclidean(y, y1) <= eps:
            return y1
        y = y1
    return y


@dataclass
class SpirualeSplatDataParserConfig:
    """Spirulae-Splat dataset config"""

    scale_factor: float = 1.0
    """How much to scale the camera origins by."""
    downscale_factor: Optional[int] = None
    """How much to downscale images. If not set, images are chosen such that the max dimension is <1600px."""
    scene_scale: float = 1.0
    """How much to scale the region of interest by."""
    orientation_method: Literal["pca", "up", "vertical", "none"] = "up"
    """The method to use for orientation."""
    center_method: Literal["poses", "focus", "none"] = "poses"
    """The method to use to center the poses."""
    auto_scale_poses: bool = True
    """Whether to automatically scale the poses to fit in +/- 1 bounding box."""
    outlier_threshold: float = float('inf')
    """Threshold to reject outlier camera poses."""
    eval_mode: Literal["fraction", "filename", "interval", "all"] = "interval"
    """
    The method to use for splitting the dataset into train and eval.
    Fraction splits based on a percentage for train and the remaining for eval.
    Filename splits based on filenames containing train/eval.
    Interval uses every nth frame for eval.
    All uses all the images for any split.
    """
    train_split_fraction: float = 0.9
    """The percentage of the dataset to use for training. Only used when eval_mode is train-split-fraction."""
    eval_interval: int = 8
    """The interval between frames to use for eval. Only used when eval_mode is eval-interval."""
    depth_unit_scale_factor: float = 1e-3
    """Scales the depth values to meters. Default value is 0.001 for a millimeter to meter conversion."""
    mask_color: Optional[Tuple[float, float, float]] = None
    """Replace the unknown pixels with this color. Relevant if you have a mask but still sample everywhere."""

    load_thumbnails: bool = True
    """Whether to load thumbnails for viewer. If False, this can significantly speed up data loading for large datasets."""

    validation_fraction: float = 0.0
    """Use this fraction of training images for validation. Stop training when performance on validation images start to drop."""


@dataclass
class SpirulaeSplatDataparser:

    config: SpirualeSplatDataParserConfig

    def __init__(self, config: SpirualeSplatDataParserConfig, dataset_dir: Path):
        self.config = config
        self.dataset_dir = dataset_dir

    def parse(self):
        return self._parse_nerfstudio_data()

    def _parse_nerfstudio_data(self, split="train"):

        assert self.dataset_dir.exists() and self.dataset_dir.is_dir(), \
            f"Data directory {self.dataset_dir} does not exist."

        transforms_path = self.dataset_dir / "transforms.json"
        assert transforms_path.exists() and transforms_path.is_file(), \
            f"File {transforms_path} does not exist."
        with open(transforms_path, 'r') as fp:
            meta = json.load(fp)

        image_filenames = []
        mask_filenames = []
        depth_filenames = []
        normal_filenames = []
        poses = []

        fx_fixed = "fl_x" in meta
        fy_fixed = "fl_y" in meta
        cx_fixed = "cx" in meta
        cy_fixed = "cy" in meta
        height_fixed = "h" in meta
        width_fixed = "w" in meta
        default_camera_model = colmap_camera_model_to_type(meta.get("camera_model", "OPENCV"))
        distort_fixed = False
        for distort_key in DISTORTION_KEYS:
            if distort_key in meta:
                distort_fixed = True
                break
        fx = []
        fy = []
        cx = []
        cy = []
        height = []
        width = []
        camera_type = []
        distort = []

        # sort the frames by fname
        fnames = []
        for frame in meta["frames"]:
            filepath = Path(frame["file_path"])
            fname = self._get_fname(filepath)
            fnames.append(fname)
        inds = np.argsort(fnames)
        frames = [meta["frames"][ind] for ind in inds]

        if np.isfinite(self.config.outlier_threshold):
            camera_positions = np.array([np.array(frame["transform_matrix"])[:3, 3] for frame in frames])
            med_camera_positions = geometric_median(camera_positions)
            distances = np.linalg.norm(camera_positions - med_camera_positions[None], axis=-1)
            mad_camera_positions = np.median(distances)
            frames = [frame for dist, frame in zip(distances, frames)
                      if dist <= self.config.outlier_threshold * mad_camera_positions]

        for frame in frames:
            filepath = Path(frame["file_path"])
            fname = self._get_fname(filepath)

            if not fx_fixed:
                assert "fl_x" in frame, "fx not specified in frame"
                fx.append(float(frame["fl_x"]))
            if not fy_fixed:
                assert "fl_y" in frame, "fy not specified in frame"
                fy.append(float(frame["fl_y"]))
            if not cx_fixed:
                assert "cx" in frame, "cx not specified in frame"
                cx.append(float(frame["cx"]))
            if not cy_fixed:
                assert "cy" in frame, "cy not specified in frame"
                cy.append(float(frame["cy"]))
            if not height_fixed:
                assert "h" in frame, "height not specified in frame"
                height.append(int(frame["h"]))
            if not width_fixed:
                assert "w" in frame, "width not specified in frame"
                width.append(int(frame["w"]))
            camera_type.append(
                colmap_camera_model_to_type(frame["camera_model"])
                if 'camera_model' in frame else default_camera_model
            )
            if not distort_fixed:
                distort.append(
                    torch.tensor([
                        float(frame.get(key, 0.0)) for key in DISTORTION_KEYS
                    ])
                )

            image_filenames.append(fname)
            poses.append(np.array(frame["transform_matrix"]))
            if "mask_path" in frame:
                mask_filepath = Path(frame["mask_path"])
                mask_fname = self._get_fname(mask_filepath)
                mask_filenames.append(mask_fname)

            if "depth_file_path" in frame:
                depth_filepath = Path(frame["depth_file_path"])
                depth_fname = self._get_fname(depth_filepath)
                depth_filenames.append(depth_fname)

            if "normal_file_path" in frame:
                normal_filepath = Path(frame["normal_file_path"])
                normal_fname = self._get_fname(normal_filepath)
                normal_filenames.append(normal_fname)

        assert len(mask_filenames) == 0 or (len(mask_filenames) == len(image_filenames)), """
        Different number of image and mask filenames.
        You should check that mask_path is specified for every frame (or zero frames) in transforms.json.
        """
        assert len(depth_filenames) == 0 or (len(depth_filenames) == len(image_filenames)), """
        Different number of image and depth filenames.
        You should check that depth_file_path is specified for every frame (or zero frames) in transforms.json.
        """
        assert len(normal_filenames) == 0 or (len(normal_filenames) == len(image_filenames)), """
        Different number of image and normal filenames.
        You should check that normal_file_path is specified for every frame (or zero frames) in transforms.json.
        """

        has_split_files_spec = any(f"{split}_filenames" in meta for split in ("train", "val", "test"))
        if f"{split}_filenames" in meta:
            # Validate split first
            split_filenames = set(self._get_fname(Path(x)) for x in meta[f"{split}_filenames"])
            unmatched_filenames = split_filenames.difference(image_filenames)
            if unmatched_filenames:
                raise RuntimeError(f"Some filenames for split {split} were not found: {unmatched_filenames}.")

            indices = [i for i, path in enumerate(image_filenames) if path in split_filenames]
            print(f"WARNING: Dataset is overriding {split}_indices to {indices}")
            indices = np.array(indices, dtype=np.int32)
        elif has_split_files_spec:
            raise RuntimeError(f"The dataset's list of filenames for split {split} is missing.")
        else:
            # find train and eval indices based on the eval_mode specified
            if self.config.eval_mode == "fraction":
                i_train, i_eval = get_train_eval_split_fraction(image_filenames, self.config.train_split_fraction)
            elif self.config.eval_mode == "filename":
                i_train, i_eval = get_train_eval_split_filename(image_filenames)
            elif self.config.eval_mode == "interval":
                i_train, i_eval = get_train_eval_split_interval(image_filenames, self.config.eval_interval)
            elif self.config.eval_mode == "all":
                i_train, i_eval = get_train_eval_split_all(image_filenames)
            else:
                raise ValueError(f"Unknown eval mode {self.config.eval_mode}")

            if split == "train":
                indices = i_train
            elif split in ["val", "test"]:
                indices = i_eval
            else:
                raise ValueError(f"Unknown dataparser split {split}")

        orientation_method = self.config.orientation_method

        poses = torch.from_numpy(np.array(poses).astype(np.float32))
        poses, transform_matrix = camera_utils.auto_orient_and_center_poses(
            poses,
            method=orientation_method,
            center_method=self.config.center_method,
        )

        # Scale poses
        scale_factor = 1.0
        if self.config.auto_scale_poses:
            scale_factor /= float(torch.max(torch.abs(poses[:, :3, 3])))
        scale_factor *= self.config.scale_factor

        poses[:, :3, 3] *= scale_factor

        # Choose image_filenames and poses based on split, but after auto orient and scaling the poses.
        image_filenames = [image_filenames[i] for i in indices]
        mask_filenames = [mask_filenames[i] for i in indices] if len(mask_filenames) > 0 else []
        depth_filenames = [depth_filenames[i] for i in indices] if len(depth_filenames) > 0 else []
        normal_filenames = [normal_filenames[i] for i in indices] if len(normal_filenames) > 0 else []

        idx_tensor = torch.tensor(indices, dtype=torch.long)
        poses = poses[idx_tensor]

        fx = float(meta["fl_x"]) if fx_fixed else torch.tensor(fx, dtype=torch.float32)[idx_tensor]
        fy = float(meta["fl_y"]) if fy_fixed else torch.tensor(fy, dtype=torch.float32)[idx_tensor]
        cx = float(meta["cx"]) if cx_fixed else torch.tensor(cx, dtype=torch.float32)[idx_tensor]
        cy = float(meta["cy"]) if cy_fixed else torch.tensor(cy, dtype=torch.float32)[idx_tensor]
        height = int(meta["h"]) if height_fixed else torch.tensor(height, dtype=torch.int32)[idx_tensor]
        width = int(meta["w"]) if width_fixed else torch.tensor(width, dtype=torch.int32)[idx_tensor]
        camera_type = [camera_type[i] for i in idx_tensor]
        if distort_fixed:
            distortion_params = torch.tensor([
                float(meta.get(key, 0.0)) for key in DISTORTION_KEYS
            ])
        else:
            distortion_params = torch.stack(distort, dim=0)[idx_tensor]

        metadata = {}

        cameras = Cameras(
            intrins=(fx, fy, cx, cy),
            distortion_params=distortion_params,
            height=height,
            width=width,
            camera_to_worlds=poses[:, :3, :4],
            camera_type=camera_type,
            metadata=metadata,
        )

        # The naming is somewhat confusing, but:
        # - transform_matrix contains the transformation to dataparser output coordinates from saved coordinates.
        # - dataparser_transform_matrix contains the transformation to dataparser output coordinates from original data coordinates.
        # - applied_transform contains the transformation to saved coordinates from original data coordinates.
        applied_transform = None
        colmap_path = self.dataset_dir / "colmap/sparse/0"
        if "applied_transform" in meta:
            applied_transform = torch.tensor(meta["applied_transform"], dtype=transform_matrix.dtype)
        elif colmap_path.exists():
            # For converting from colmap, this was the effective value of applied_transform that was being
            # used before we added the applied_transform field to the output dataformat.
            meta["applied_transform"] = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0]]
            applied_transform = torch.tensor(meta["applied_transform"], dtype=transform_matrix.dtype)

        if applied_transform is not None:
            dataparser_transform_matrix = transform_matrix @ torch.cat(
                [applied_transform, torch.tensor([[0, 0, 0, 1]], dtype=transform_matrix.dtype)], 0
            )
        else:
            dataparser_transform_matrix = transform_matrix

        if "applied_scale" in meta:
            applied_scale = float(meta["applied_scale"])
            scale_factor *= applied_scale

        # reinitialize metadata for dataparser_outputs
        metadata = {}

        # _generate_dataparser_outputs might be called more than once so we check if we already loaded the point cloud
        try:
            self.prompted_user
        except AttributeError:
            self.prompted_user = False

        # Load 3D points
        assert "ply_file_path" in meta, "No initial point cloud found in transforms.json"
        ply_file_path = self.dataset_dir / meta["ply_file_path"]
        sparse_points = self._load_3D_points(ply_file_path, transform_matrix, scale_factor)
        if sparse_points is not None:
            metadata.update(sparse_points)

        dataparser_outputs = dict(
            cameras=cameras,
            image_filenames=image_filenames,
            mask_filenames=mask_filenames if len(mask_filenames) > 0 else None,
            dataparser_scale=scale_factor,
            dataparser_transform=dataparser_transform_matrix,
            metadata={
                "depth_filenames": depth_filenames if len(depth_filenames) > 0 else None,
                "depth_unit_scale_factor": self.config.depth_unit_scale_factor,
                "normal_filenames": normal_filenames if len(normal_filenames) > 0 else None,
                "mask_color": self.config.mask_color,
                "load_thumbnails": self.config.load_thumbnails,
                "val_indices": get_train_eval_split_fraction(image_filenames, 1-self.config.validation_fraction)[1].tolist(),
                **metadata,
            },
        )
        return dataparser_outputs

    def _load_3D_points(self, ply_file_path: Path, transform_matrix: torch.Tensor, scale_factor: float):
        """Loads point clouds positions and colors from .ply

        Args:
            ply_file_path: Path to .ply file
            transform_matrix: Matrix to transform world coordinates
            scale_factor: How much to scale the camera origins by.

        Returns:
            A dictionary of points: points3D_xyz and colors: points3D_rgb
        """
        from time import perf_counter
        import numpy as np
        import torch

        time0 = perf_counter()
        from plyfile import PlyData
        time1 = perf_counter()
        # print("Import plyfile:", time1 - time0)

        time0 = perf_counter()
        plydata = PlyData.read(str(ply_file_path))

        # Extract vertex data
        v = plydata['vertex']
        
        # Check if points exist
        if v.count == 0:
            return None

        # Stack x, y, z into a (N, 3) float32 array
        points_np = np.stack([v['x'], v['y'], v['z']], axis=-1).astype(np.float32)
        points3D = torch.from_numpy(points_np)

        # Apply transformation and scale
        points3D = (
            torch.cat(
                (
                    points3D,
                    torch.ones_like(points3D[..., :1]),
                ),
                -1,
            )
            @ transform_matrix.T
        )
        points3D *= scale_factor

        # Stack r, g, b into a (N, 3) uint8 array
        colors_np = np.stack([v['red'], v['green'], v['blue']], axis=-1).astype(np.uint8)
        points3D_rgb = torch.from_numpy(colors_np)

        out = {
            "points3D_xyz": points3D,
            "points3D_rgb": points3D_rgb,
        }
        time1 = perf_counter()
        # print("Load PLY:", time1-time0)
        return out

    def _get_fname(self, filepath: Path) -> Path:
        return self.dataset_dir / filepath  # TODO: downscale
