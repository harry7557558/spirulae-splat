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
from typing import Literal, Optional, Tuple, Type, List, Union

import numpy as np
import torch
import math
import os

import json

from spirulae_splat.modules.camera import (
    colmap_camera_model_to_type,
    Cameras,
)
from spirulae_splat.modules.colmap_utils import (
    qvec2rotmat,
    load_colmap_cameras,
    load_colmap_points3D,
    load_colmap_images,
    parse_colmap_camera_params,
)
from spirulae_splat.modules.metashape_utils import (
    metashape_to_json,
    find_metashape_cameras_dict,
    ET,
)

from spirulae_splat.modules import camera_utils


def get_train_eval_split_fraction(image_filenames: List, train_split_fraction: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the train/eval split fraction based on the number of images and the train split fraction.

    Args:
        image_filenames: list of image filenames
        train_split_fraction: fraction of images to use for training
    """

    # filter image_filenames and poses based on train/eval split percentage
    num_images = len(image_filenames)
    num_train_images = math.ceil(num_images * train_split_fraction)
    num_eval_images = num_images - num_train_images
    i_all = np.arange(num_images)
    i_train = np.linspace(
        0, num_images - 1, num_train_images, dtype=int
    )  # equally spaced training images starting and ending at 0 and num_images-1
    i_eval = np.setdiff1d(i_all, i_train)  # eval images are the remaining images
    assert len(i_eval) == num_eval_images

    return i_train, i_eval


def get_train_eval_split_filename(image_filenames: List) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the train/eval split based on the filename of the images.

    Args:
        image_filenames: list of image filenames
    """

    num_images = len(image_filenames)
    basenames = [os.path.basename(image_filename) for image_filename in image_filenames]
    i_all = np.arange(num_images)
    i_train = []
    i_eval = []
    for idx, basename in zip(i_all, basenames):
        # check the frame index
        if "train" in basename:
            i_train.append(idx)
        elif "eval" in basename:
            i_eval.append(idx)
        else:
            raise ValueError("frame should contain train/eval in its name to use this eval-frame-index eval mode")

    return np.array(i_train), np.array(i_eval)


def get_train_eval_split_interval(image_filenames: List, eval_interval: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the train/eval split based on the interval of the images.

    Args:
        image_filenames: list of image filenames
        eval_interval: interval of images to use for eval
    """

    num_images = len(image_filenames)
    all_indices = np.arange(num_images)
    train_indices = all_indices[all_indices % eval_interval != 0]
    eval_indices = all_indices[all_indices % eval_interval == 0]
    i_train = train_indices
    i_eval = eval_indices

    return i_train, i_eval


def get_train_eval_split_all(image_filenames: List) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the train/eval split where all indices are used for both train and eval.

    Args:
        image_filenames: list of image filenames
    """
    num_images = len(image_filenames)
    i_all = np.arange(num_images)
    i_train = i_all
    i_eval = i_all
    return i_train, i_eval


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
class SpirulaeSplatDataParserConfig:
    """Spirulae-Splat dataset config"""

    data_format: Literal["colmap", "nerfstudio", "metashape", None] = None
    """Data format, leave None to auto detect"""
    colmap_recon_dir: Optional[Path] = None
    """Path to COLMAP reconstruction relative to dataset directory (e.g. sparse/0). Will auto detect if not specified."""
    image_dir: Path = Path("images")
    """Path to images relative to dataset directory, used for COLMAP and Metashape datasets"""
    mask_dir: Path = Path("masks")
    """Path to masks relative to dataset directory, used for COLMAP and Metashape datasets"""
    depth_dir: Path = Path("depths")
    """Path to depth maps relative to dataset directory, used for COLMAP and Metashape datasets"""
    normal_dir: Path = Path("normals")
    """Path to normal maps relative to dataset directory, used for COLMAP and Metashape datasets"""
    metashape_xml: Optional[Path] = None
    """Path to the Metashape xml file. Will automatically detect if not specified."""
    metashape_ply: Optional[Path] = None
    """Path to the Metashape point export ply file. Will automatically detect if not specified."""
    metashape_psx: Optional[Path] = None
    """Path to Metashape PSX file, used to resolve file name ambiguity when there are multiple images with the same file name"""
    rescale_camera_to_fit: Union[bool, int] = False
    """Whether to check if image resolution match camera resolution and scale camera intrinsics accordingly if not.
        Set this to a number to divide intrinsics by that number, e.g. Mip-NeRF 360 and Zip-NeRF with images_(2|4)
        Set this to True to detect resolution, e.g. tankt_db"""

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
    eval_mode: Literal["fraction", "filename", "interval", "all"] = "all"
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

    validation_fraction: float = 0.0
    """Use this fraction of training images for validation. Stop training when performance on validation images start to drop."""


@dataclass
class SpirulaeSplatDataparser:

    config: SpirulaeSplatDataParserConfig

    def __init__(self, config: SpirulaeSplatDataParserConfig, dataset_dir: Path):
        self.config = config
        self.dataset_dir = dataset_dir

    def parse(self):
        if self.config.data_format == "colmap":
            return self._parse_colmap_data()
        if self.config.data_format == "nerfstudio":
            return self._parse_nerfstudio_data()
        if self.config.data_format == "metashape":
            return self._parser_metashape_data()

        try:
            print("Attempting to parse Nerfstudio data...")
            return self._parse_nerfstudio_data()
        except BaseException as e:
            print("Failed to parse Nerfstudio data:", ' '.join(map(str, e.args)))
        try:
            print("Attempting to parse COLMAP data...")
            return self._parse_colmap_data()
        except BaseException as e:
            print("Failed to parse COLMAP data:", ' '.join(map(str, e.args)))
        try:
            print("Attempting to parse Metashape data...")
            return self._parser_metashape_data()
        except BaseException as e:
            print("Failed to parse Metashape data:", ' '.join(map(str, e.args)))
        raise ValueError("No supported dataset format detected. Make sure you have a supported Nerfstudio, COLMAP, or Metashape dataset.")

    def _parse_nerfstudio_data(self, _meta=None, _points3D=None):

        if not (self.dataset_dir.exists() and self.dataset_dir.is_dir()):
            raise ValueError(f"Data directory {self.dataset_dir} does not exist.")

        if _meta is None:
            transforms_path = self.dataset_dir / "transforms.json"
            if not (transforms_path.exists() and transforms_path.is_file()):
                raise ValueError(f"File {transforms_path} does not exist.")
            with open(transforms_path, 'r') as fp:
                meta = json.load(fp)
        else:
            meta = _meta

        image_filenames = []
        mask_filenames = []
        depth_filenames = []
        normal_filenames = []
        poses = []

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

        if isinstance(self.config.rescale_camera_to_fit, bool) and self.config.rescale_camera_to_fit:
            from PIL import Image

        for frame in frames:
            filepath = Path(frame["file_path"])
            fname = self._get_fname(filepath)

            fx.append(float(frame["fl_x"] if 'fl_x' in frame else meta['fl_x']))
            fy.append(float(frame["fl_y"] if 'fl_y' in frame else meta['fl_y']))
            cx.append(float(frame["cx"] if 'cx' in frame else meta['cx']))
            cy.append(float(frame["cy"] if 'cy' in frame else meta['cy']))
            height.append(int(frame["h"] if 'h' in frame else meta['h']))
            width.append(int(frame["w"] if 'w' in frame else meta['w']))
            camera_type.append(
                colmap_camera_model_to_type(
                    frame["camera_model"]
                    if "camera_model" in frame else meta.get("camera_model", "OPENCV")
                )
            )
            distort.append(
                torch.tensor([
                    float(frame.get(key, meta.get(key, 0.0))) for key in DISTORTION_KEYS
                ])
            )

            if not isinstance(self.config.rescale_camera_to_fit, bool):
                fx[-1] /= self.config.rescale_camera_to_fit
                fy[-1] /= self.config.rescale_camera_to_fit
                cx[-1] /= self.config.rescale_camera_to_fit
                cy[-1] /= self.config.rescale_camera_to_fit
                height[-1] //= self.config.rescale_camera_to_fit
                width[-1] //= self.config.rescale_camera_to_fit
            elif isinstance(self.config.rescale_camera_to_fit, bool) and self.config.rescale_camera_to_fit:
                with Image.open(fname) as img:
                    w, h = img.size
                sx, sy = w / width[-1], h / height[-1]
                fx[-1] *= sx
                fy[-1] *= sy
                cx[-1] *= sx
                cy[-1] *= sy
                width[-1] = w
                height[-1] = h

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

        if not (len(mask_filenames) == 0 or (len(mask_filenames) == len(image_filenames))):
            raise ValueError("""
        Different number of image and mask filenames.
        You should check that mask_path is specified for every frame (or zero frames) in transforms.json.
        """)
        if not (len(depth_filenames) == 0 or (len(depth_filenames) == len(image_filenames))):
            raise ValueError("""
        Different number of image and depth filenames.
        You should check that depth_file_path is specified for every frame (or zero frames) in transforms.json.
        """)
        if not (len(normal_filenames) == 0 or (len(normal_filenames) == len(image_filenames))):
            raise ValueError("""
        Different number of image and normal filenames.
        You should check that normal_file_path is specified for every frame (or zero frames) in transforms.json.
        """)

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

        # Auto orient poses
        poses = torch.from_numpy(np.array(poses).astype(np.float32))
        poses, transform_matrix = camera_utils.auto_orient_and_center_poses(
            poses,
            method=self.config.orientation_method,
            center_method=self.config.center_method,
        )

        # Auto scale poses
        scale_factor = 1.0
        if self.config.auto_scale_poses:
            scale_factor /= float(torch.max(torch.abs(poses[:, :3, 3])))
        scale_factor *= self.config.scale_factor
        poses[:, :3, 3] *= scale_factor

        # Load 3D points
        metadata = {}
        if _points3D is None:
            if "ply_file_path" not in meta:
                raise ValueError("No initial point cloud found in transforms.json")
            ply_file_path = self.dataset_dir / meta["ply_file_path"]
            sparse_points = self._load_3D_points(ply_file_path, transform_matrix, scale_factor)
            metadata.update(sparse_points)
        else:
            metadata.update(_points3D)

        # Apply transformation to 3D points
        metadata['points3D_xyz'] = (
            torch.cat(
                (
                    metadata['points3D_xyz'],
                    torch.ones_like(metadata['points3D_xyz'][..., :1]),
                ),
                -1,
            )
            @ transform_matrix.T
        )
        metadata['points3D_xyz'] *= scale_factor

        all_dataparser_outputs = []
        for split_id, indices in enumerate([i_train, i_eval]):

            # Choose image_filenames and poses based on split, but after auto orient and scaling the poses.
            image_filenames_split = [image_filenames[i] for i in indices]
            mask_filenames_split = [mask_filenames[i] for i in indices] if len(mask_filenames) > 0 else []
            depth_filenames_split = [depth_filenames[i] for i in indices] if len(depth_filenames) > 0 else []
            normal_filenames_split = [normal_filenames[i] for i in indices] if len(normal_filenames) > 0 else []

            idx_tensor = torch.tensor(indices, dtype=torch.long)
            poses_split = poses[idx_tensor]

            fx_split = torch.tensor(fx, dtype=torch.float32)[idx_tensor]
            fy_split = torch.tensor(fy, dtype=torch.float32)[idx_tensor]
            cx_split = torch.tensor(cx, dtype=torch.float32)[idx_tensor]
            cy_split = torch.tensor(cy, dtype=torch.float32)[idx_tensor]
            height_split = torch.tensor(height, dtype=torch.int32)[idx_tensor]
            width_split = torch.tensor(width, dtype=torch.int32)[idx_tensor]
            camera_type_split = [camera_type[i] for i in idx_tensor]
            distortion_params_split = torch.stack(distort, dim=0)[idx_tensor]

            cameras = Cameras(
                intrins=(fx_split, fy_split, cx_split, cy_split),
                distortion_params=distortion_params_split,
                height=height_split,
                width=width_split,
                camera_to_worlds=poses_split[:, :3, :4],
                camera_type=camera_type_split,
                metadata={},
            )

            # The naming is somewhat confusing, but:
            # - transform_matrix contains the transformation to dataparser output coordinates from saved coordinates.
            # - dataparser_transform_matrix contains the transformation to dataparser output coordinates from original data coordinates.
            # - applied_transform contains the transformation to saved coordinates from original data coordinates.
            applied_transform = None
            if "applied_transform" in meta:
                applied_transform = torch.tensor(meta["applied_transform"], dtype=transform_matrix.dtype)
                dataparser_transform_matrix = transform_matrix @ torch.cat(
                    [applied_transform, torch.tensor([[0, 0, 0, 1]], dtype=transform_matrix.dtype)], 0
                )
            else:
                dataparser_transform_matrix = transform_matrix

            if "applied_scale" in meta:
                applied_scale = float(meta["applied_scale"])
                scale_factor *= applied_scale

            dataparser_outputs = dict(
                cameras=cameras,
                image_filenames=image_filenames_split,
                mask_filenames=mask_filenames_split if len(mask_filenames_split) > 0 else None,
                dataparser_scale=scale_factor,
                dataparser_transform=dataparser_transform_matrix,
                metadata={
                    "depth_filenames": depth_filenames_split if len(depth_filenames_split) > 0 else None,
                    "depth_unit_scale_factor": self.config.depth_unit_scale_factor,
                    "normal_filenames": normal_filenames_split if len(normal_filenames_split) > 0 else None,
                    "mask_color": self.config.mask_color,
                    **metadata,
                },
            )
            if split_id == 0:
                dataparser_outputs["metadata"]["val_indices"] = get_train_eval_split_fraction(
                    image_filenames_split, 1-self.config.validation_fraction)[1].tolist()

            all_dataparser_outputs.append(dataparser_outputs)

        return (*all_dataparser_outputs,)

    def _parse_colmap_data(self):
        if self.config.colmap_recon_dir is not None:
            recon_dir = self.dataset_dir / self.config.colmap_recon_dir
            colmap_points = load_colmap_points3D(recon_dir)
            cam_id_to_camera = load_colmap_cameras(recon_dir)
            im_id_to_image = load_colmap_images(recon_dir)
        else:
            okay = False
            for colmap_recon_dir in [
                Path("sparse/0"),
                Path("colmap/sparse/0"),
                Path("sparse"),
                Path("colmap"),
                Path("")
            ]:
                try:
                    recon_dir = self.dataset_dir / colmap_recon_dir
                    colmap_points = load_colmap_points3D(recon_dir)
                    cam_id_to_camera = load_colmap_cameras(recon_dir)
                    im_id_to_image = load_colmap_images(recon_dir)
                except:
                    continue
                okay = True
                print(f"Loaded COLMAP reconstruction from {recon_dir}")
            if not okay:
                raise ValueError("Could not find COLMAP reconstruction dir containing points, cameras, and images files. Specify --dataparser.colmap_recon_dir if needed.")

        frames = []
        for im in im_id_to_image.values():
            # intrinsics
            camera = cam_id_to_camera[im.camera_id]
            frame = parse_colmap_camera_params(camera)

            # extrinsics
            rotation = qvec2rotmat(im.qvec)
            translation = im.tvec.reshape(3, 1)
            c2w = np.eye(4)
            c2w[:3, :3] = rotation.T * [[1, -1, -1]]
            c2w[:3, 3:] = -rotation.T @ translation
            frame["transform_matrix"] = c2w

            # images
            image_filename = self.dataset_dir / self.config.image_dir / Path(im.name)
            if not image_filename.exists():
                raise ValueError(f"{image_filename} does not exist. Specify `--dataparser.image_dir` if needed.")
            frame['file_path'] = str(self.config.image_dir / Path(im.name))

            self._add_auxiliary_buffers(frame)
            frames.append(frame)

        points = [*colmap_points.values()]
        xyz = torch.from_numpy(np.stack([p.xyz for p in points])).float()
        rgb = torch.from_numpy(np.stack([p.rgb for p in points])).byte()
        return self._parse_nerfstudio_data(
            _meta={'frames': frames},
            _points3D={'points3D_xyz': xyz, 'points3D_rgb': rgb}
        )

    def _parser_metashape_data(self):

        if not os.path.exists(self.dataset_dir):
            raise ValueError("Work directory " + self.dataset_dir + " not found")
        dataset_dir_files = [Path(fn) for fn in os.listdir(self.dataset_dir)]

        # Load .xml file
        if self.config.metashape_xml is None:
            xml_files = [self.dataset_dir / fn for fn in dataset_dir_files if fn.suffix.lower() == ".xml"]
            if len(xml_files) == 0:
                raise ValueError("No XML file found in dataset_dir. Please specify using --dataparser.metashape_xml")
            if len(xml_files) > 1:
                raise ValueError("Multiple XML file found in dataset_dir. Please specify using --dataparser.metashape_xml")
            self.config.metashape_xml = xml_files[0]
            print("Using XML file found:", self.config.metashape_xml)
        elif not self.config.metashape_xml.exists():
            self.config.metashape_xml = self.dataset_dir / self.config.metashape_xml
        if not self.config.metashape_xml.exists():
            raise ValueError(f"XML file {self.config.metashape_xml} doesn't exist")
        if self.config.metashape_xml.suffix.lower() != ".xml":
            raise ValueError(f"XML file {self.config.metashape_xml} must have a .xml extension")

        # Load .ply file
        if self.config.metashape_ply is None:
            ply_files = [fn for fn in dataset_dir_files if fn.suffix.lower() == ".ply"]
            if len(ply_files) == 0:
                raise ValueError("No ply file found in dataset_dir. Please specify using --dataparser.metashape_ply")
            if len(ply_files) > 1:
                raise ValueError("Multiple ply file found in dataset_dir. Please specify using --dataparser.metashape_ply")
            self.config.metashape_ply = ply_files[0]
            print("Using PLY file found:", self.dataset_dir / self.config.metashape_ply)
        if not (self.dataset_dir / self.config.metashape_ply).exists():
            raise ValueError(f"ply file {self.dataset_dir / self.config.metashape_ply} doesn't exist")
        if self.config.metashape_ply.suffix.lower() != ".ply":
            raise ValueError(f"ply file {self.dataset_dir / self.config.metashape_ply} must have a .ply extension")

        # Load .psx file
        if self.config.metashape_psx is None:
            psx_files = [self.dataset_dir / fn for fn in dataset_dir_files if fn.suffix.lower() == ".psx"]
            if len(psx_files) > 1:
                raise ValueError("Multiple psx file found in dataset_dir. Please specify using --dataparser.metashape_psx")
            if len(psx_files) != 0:
                self.config.metashape_psx = psx_files[0]
                print("Using PSX file found:", self.config.metashape_psx)
        if self.config.metashape_psx is not None and not self.config.metashape_psx.exists():
            self.config.metashape_psx = self.dataset_dir / self.config.metashape_psx
            if not self.config.metashape_psx.exists():
                raise ValueError(f"psx file {self.config.metashape_psx} doesn't exist")
        camera_dict = None
        if self.config.metashape_psx is not None:
            if self.config.metashape_psx.suffix.lower() != ".psx":
                raise ValueError(f"psx file {self.config.metashape_psx} must have a .psx extension")
            root = ET.parse(self.config.metashape_psx).getroot()
            metashape_dir = self.config.metashape_psx.parent / root.get("path").replace("{projectname}", self.config.metashape_psx.name.rstrip(".psx"))
            camera_dict = find_metashape_cameras_dict(metashape_dir.parent)

        image_dir = self.dataset_dir / self.config.image_dir
        if not image_dir.exists():
            raise ValueError(f"Image directory `{image_dir}` does not exist. Please specify using --dataparser.image_dir")

        summary_log = []

        image_filenames = []
        for file_path in image_dir.rglob("*"):
            if file_path.is_file():
                file_path = file_path.relative_to(self.dataset_dir)
                image_filenames.append(str(file_path))

        transforms, summary = metashape_to_json(
            xml_filename=self.config.metashape_xml,
            ply_filename=self.config.metashape_ply,
            image_filenames=image_filenames,
            camera_dict=camera_dict,
        )
        summary_log.extend(summary)

        if len(transforms) > 1:
            print("WARNING: Multiple components found in Metashape export, which can have unintended effects. Clean up if needed.")

        for summary in summary_log:
            print(summary)

        for frame in transforms[0]['frames']:
            self._add_auxiliary_buffers(frame)
        return self._parse_nerfstudio_data(transforms[0])

    def _add_auxiliary_buffers(self, frame):
        name = str(Path(frame['file_path']).relative_to(self.config.image_dir))

        for mask_filename in [
            Path(name+".png"),
            Path(name).with_suffix(".png"),
            Path(name+".PNG"),
            Path(name).with_suffix(".PNG"),
            Path(name+".jpg"),
            Path(name).with_suffix(".jpg"),
            Path(name+".JPG"),
            Path(name).with_suffix(".JPG"),
            Path(name+".jpeg"),
            Path(name).with_suffix(".jpeg"),
            Path(name+".JPEG"),
            Path(name).with_suffix(".JPEG"),
            Path(name[:name.rfind('.')]+"_mask.png"),
        ]:
            if (self.dataset_dir / self.config.mask_dir / mask_filename).exists():
                frame['mask_path'] = self.config.mask_dir / mask_filename
                break

        # depths
        for depth_filename in [
            Path(name).with_suffix(".png"),
            Path(name).with_suffix(".jpg"),
            Path(name).with_suffix(".jpeg"),
            Path(name).with_suffix(".npy"),
            Path(name).with_suffix(".npz"),
            Path(name+".png"),
            Path(name+".jpg"),
            Path(name+".jpeg"),
            Path(name+".npy"),
            Path(name+".npz"),
            Path(name).with_suffix(".PNG"),
            Path(name).with_suffix(".JPG"),
            Path(name).with_suffix(".JPEG"),
            Path(name[:name.rfind('.')]+"_depth.png"),
            Path(name[:name.rfind('.')]+"_depth.jpg"),
            Path(name[:name.rfind('.')]+"_depth.jpeg"),
        ]:
            if (self.dataset_dir / self.config.depth_dir / depth_filename).exists():
                frame['depth_file_path'] = self.config.depth_dir / depth_filename
                break

        # normals
        for normal_filename in [
            Path(name).with_suffix(".png"),
            Path(name).with_suffix(".jpg"),
            Path(name).with_suffix(".jpeg"),
            Path(name).with_suffix(".npy"),
            Path(name).with_suffix(".npz"),
            Path(name+".png"),
            Path(name+".jpg"),
            Path(name+".jpeg"),
            Path(name+".npy"),
            Path(name+".npz"),
            Path(name).with_suffix(".PNG"),
            Path(name).with_suffix(".JPG"),
            Path(name).with_suffix(".JPEG"),
            Path(name[:name.rfind('.')]+"_normal.png"),
            Path(name[:name.rfind('.')]+"_normal.jpg"),
            Path(name[:name.rfind('.')]+"_normal.jpeg"),
        ]:
            if (self.dataset_dir / self.config.normal_dir / normal_filename).exists():
                frame['normal_file_path'] = self.config.normal_dir / normal_filename
                break


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
            raise ValueError(f"Failed to load points from {ply_file_path}.")

        # Stack x, y, z into a (N, 3) float32 array
        points_np = np.stack([v['x'], v['y'], v['z']], axis=-1).astype(np.float32)
        points3D = torch.from_numpy(points_np)

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
