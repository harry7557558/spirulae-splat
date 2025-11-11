#!/usr/bin/env python3

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

"""Helper utils for processing metashape data into the nerfstudio format."""

import json
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass
import tyro

import sys
import os

import numpy as np
import open3d as o3d

from nerfstudio.process_data.process_data_utils import CAMERA_MODELS
from nerfstudio.utils.rich_utils import CONSOLE

from difflib import SequenceMatcher


class StringMatcher:
    """Brute force solution, fast enough"""

    def __init__(self, strings):
        self.strings = strings
        # Precompute reversed strings for faster suffix matching
        self.reversed = [s[::-1] for s in strings]

    def query_suffix(self, s):
        """Find all strings with the longest common suffix with s."""
        rs = s[::-1]
        best_len, best_matches = 0, []
        for orig, rev in zip(self.strings, self.reversed):
            # find length of common prefix of reversed strings (i.e. suffix of original)
            i = 0
            while i < len(rs) and i < len(rev) and rs[i] == rev[i]:
                i += 1
            if i > best_len:
                best_len, best_matches = i, [orig]
            elif i == best_len and i > 0:
                best_matches.append(orig)
        return best_matches, best_len

    def query_substring(self, s):
        """Find all strings with the longest common substring with s."""
        best_len, best_matches = 0, []
        for orig in self.strings:
            match_len = SequenceMatcher(None, s, orig).find_longest_match(0, len(s), 0, len(orig)).size
            if match_len > best_len:
                best_len, best_matches = match_len, [orig]
            elif match_len == best_len and match_len > 0:
                best_matches.append(orig)
        return best_matches, best_len


def find_cameras_dict(root_dir):
    """
    Recursively search a directory (and one level inside .zip files) for an XML file
    containing a <cameras> element. If found, return a dict mapping:
        camera_id (str) -> photo path (str)
    Returns None if no valid <cameras> element is found.
    """
    def parse_cameras_element(root):
        cameras = root.find("cameras")
        if cameras is None:
            return None
        cam_dict = {}
        for cam in cameras.findall("camera"):
            cam_id = cam.get("camera_id")
            photo = cam.find("photo")
            if cam_id is not None and photo is not None:
                path = photo.get("path")
                if path is not None:
                    cam_dict[cam_id] = path
        return cam_dict if cam_dict else None

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)

            # Case 1: plain XML file
            if filename.lower().endswith(".xml"):
                try:
                    tree = ET.parse(filepath)
                    result = parse_cameras_element(tree.getroot())
                    if result is not None:
                        return result
                except ET.ParseError:
                    continue

            # Case 2: ZIP file containing XMLs
            elif filename.lower().endswith(".zip"):
                try:
                    with zipfile.ZipFile(filepath, "r") as zf:
                        for name in zf.namelist():
                            if name.lower().endswith(".xml"):
                                try:
                                    with zf.open(name) as xml_file:
                                        tree = ET.parse(xml_file)
                                        result = parse_cameras_element(tree.getroot())
                                        if result is not None:
                                            return result
                                except ET.ParseError:
                                    continue
                except zipfile.BadZipFile:
                    continue

    return None


def _find_param(calib_xml: ET.Element, param_name: str):
    param = calib_xml.find(param_name)
    if param is not None:
        return float(param.text)  # type: ignore
    return 0.0


def metashape_to_json(
    xml_filename: Path,
    output_dir: Path,
    ply_filename: Path,
    image_filenames: List[str],
    camera_dict: Optional[ET.ElementTree] = None
) -> List[str]:
    """Convert Metashape data into a nerfstudio dataset.

    Args:
        xml_filename: Path to the metashape cameras xml file.
        output_dir: Path to the output directory.
        ply_filename: Path to the exported ply file.
        verbose: Whether to print verbose output.

    Returns:
        Summary of the conversion.
    """

    xml_tree = ET.parse(xml_filename)
    root = xml_tree.getroot()
    chunk = root[0]
    sensors = chunk.find("sensors")

    if sensors is None:
        raise ValueError("No sensors found")

    calibrated_sensors = [
        sensor for sensor in sensors.iter("sensor") if sensor.get("type") == "spherical" or sensor.find("calibration") is not None
    ]
    if not calibrated_sensors:
        raise ValueError("No calibrated sensor found in Metashape XML")
    sensor_type = [s.get("type") for s in calibrated_sensors]
    if sensor_type.count(sensor_type[0]) != len(sensor_type):
        raise ValueError(
            "All Metashape sensors do not have the same sensor type. "
            "nerfstudio does not support per-frame camera_model types."
            "Only one camera type can be used: frame, fisheye or spherical (perspective, fisheye or equirectangular)"
        )
    data = {}
    if sensor_type[0] == "frame":
        data["camera_model"] = CAMERA_MODELS["perspective"].value
    elif sensor_type[0] == "fisheye":
        data["camera_model"] = CAMERA_MODELS["fisheye"].value
    elif sensor_type[0] == "spherical":
        data["camera_model"] = CAMERA_MODELS["equirectangular"].value
    else:
        # Cylindrical and RPC sensor types are not supported
        raise ValueError(f"Unsupported Metashape sensor type '{sensor_type[0]}'")

    sensor_dict = {}
    for sensor in calibrated_sensors:
        s = {}
        resolution = sensor.find("resolution")
        assert resolution is not None, "Resolution not found in Metashape xml"
        s["w"] = int(resolution.get("width"))  # type: ignore
        s["h"] = int(resolution.get("height"))  # type: ignore

        # https://github.com/facebookresearch/EyefulTower/issues/7
        # https://www.agisoft.com/pdf/metashape-pro_2_2_en.pdf Appendix D. Camera Models
        calib = sensor.find("calibration")
        if calib is None:
            assert sensor_type[0] == "spherical", "Only spherical sensors should have no intrinsics"
            s["fl_x"] = s["w"] / 2.0
            s["fl_y"] = s["h"]
            s["cx"] = s["w"] / 2.0
            s["cy"] = s["h"] / 2.0
        else:
            f = calib.find("f")
            assert f is not None, "Focal length not found in Metashape xml"
            s["fl_x"] = s["fl_y"] = float(f.text)  # type: ignore
            s["cx"] = _find_param(calib, "cx") + s["w"] / 2.0  # type: ignore
            s["cy"] = _find_param(calib, "cy") + s["h"] / 2.0  # type: ignore

            s["k1"] = _find_param(calib, "k1")
            s["k2"] = _find_param(calib, "k2")
            s["k3"] = _find_param(calib, "k3")
            s["k4"] = _find_param(calib, "k4")
            s["p1"] = _find_param(calib, "p2")
            s["p2"] = _find_param(calib, "p1")
            s["b1"] = _find_param(calib, "b1") / s["fl_x"]
            s["b2"] = _find_param(calib, "b2") / s["fl_x"]

        sensor_dict[sensor.get("id")] = s

    components = chunk.find("components")
    component_dict = {}
    if components is not None:
        for component in components.iter("component"):
            transform = component.find("transform")
            if transform is not None:
                rotation = transform.find("rotation")
                if rotation is None:
                    r = np.eye(3)
                else:
                    assert isinstance(rotation.text, str)
                    r = np.array([float(x) for x in rotation.text.split()]).reshape((3, 3))
                translation = transform.find("translation")
                if translation is None:
                    t = np.zeros(3)
                else:
                    assert isinstance(translation.text, str)
                    t = np.array([float(x) for x in translation.text.split()])
                scale = transform.find("scale")
                if scale is None:
                    s = 1.0
                else:
                    assert isinstance(scale.text, str)
                    s = float(scale.text)

                m = np.eye(4)
                m[:3, :3] = r
                m[:3, 3] = t / s
                component_dict[component.get("id")] = m

    image_filename_matcher = StringMatcher(image_filenames)

    frames = []
    cameras = chunk.find("cameras")
    assert cameras is not None, "Cameras not found in Metashape xml"
    num_skipped = 0
    for camera in cameras.iter("camera"):
        frame = {}
        if camera_dict is not None:
            image_filename = camera_dict[camera.get("id")]
            matches = image_filename_matcher.query_suffix(image_filename)[0]
            if len(matches) != 1:
                CONSOLE.log("WARNING: ambiguous filenames", matches, ", Skipping")
                num_skipped += 1
                continue
            frame["file_path"] = matches[0]
        else:
            camera_label = camera.get("label")
            matches = image_filename_matcher.query_substring(camera_label)[0]
            assert isinstance(camera_label, str)
            if len(matches) != 1:
                CONSOLE.log("WARNING: ambiguous filenames", matches, ", Skipping.\n"
                            "Specify Metashape .psx file using --psx to attempt resolving ambiguity.")
                num_skipped += 1
                continue
            frame["file_path"] = matches[0]

        sensor_id = camera.get("sensor_id")
        if sensor_id not in sensor_dict:
            # this should only happen when we have a sensor that doesn't have calibration
            CONSOLE.print(f"Missing sensor calibration for {camera.get('label')}, Skipping")
            num_skipped += 1
            continue
        # Add all sensor parameters to this frame.
        frame.update(sensor_dict[sensor_id])

        if camera.find("transform") is None:
            CONSOLE.print(f"Missing transforms data for {camera.get('label')}, Skipping")
            num_skipped += 1
            continue
        transform = np.array([float(x) for x in camera.find("transform").text.split()]).reshape((4, 4))  # type: ignore

        component_id = camera.get("component_id")
        if component_id in component_dict:
            transform = component_dict[component_id] @ transform

        # Metashape camera is looking towards -Z, +X is to the right and +Y is to the top/up of the first cam
        # Rotate the scene according to nerfstudio convention
        transform = transform[[2, 0, 1, 3], :]
        # Convert from Metashape's camera coordinate system (OpenCV) to ours (OpenGL)
        transform[:, 1:3] *= -1
        frame["transform_matrix"] = transform.tolist()
        frames.append(frame)

    data["frames"] = frames
    applied_transform = np.eye(4)[:3, :]
    applied_transform = applied_transform[np.array([2, 0, 1]), :]
    data["applied_transform"] = applied_transform.tolist()

    summary = []

    if ply_filename is not None:
        assert ply_filename.name != "sparse_pc.ply", "PLY file name must not be sparse_pc.ply"
        assert ply_filename.exists()
        pc = o3d.io.read_point_cloud(str(ply_filename))
        points3D = np.asarray(pc.points)
        points3D = np.einsum("ij,bj->bi", applied_transform[:3, :3], points3D) + applied_transform[:3, 3]
        pc.points = o3d.utility.Vector3dVector(points3D)
        o3d.io.write_point_cloud(str(output_dir / "sparse_pc.ply"), pc)
        data["ply_file_path"] = "sparse_pc.ply"
        summary.append(f"Imported {ply_filename} as starting points")

    with open(output_dir / "transforms.json", "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)

    if num_skipped == 1:
        summary.append(f"{num_skipped} image skipped.")
    if num_skipped > 1:
        summary.append(f"{num_skipped} images were skipped.")

    summary.append(f"Final dataset is {len(data['frames'])} frames.")

    return summary



@dataclass
class ProcessMetashape():
    """Process Metashape data into a nerfstudio dataset.

    This script assumes that cameras have been aligned using Metashape. After alignment, it is necessary to export the
    camera poses as a `.xml` and point clouds as a `.ply` file under the dataset folder. For point cloud, Make sure to
    export the data in non-binary format and exclude the normals. Optionally, provide the Metashape `.psx` file to 
    resolve file name ambiguity (Tested with Metashape Standard Version 2.2.2).

    """

    work_dir: Path
    """Path to the work directory."""

    xml: Optional[Path] = None
    """Path to the Metashape xml file. Will automatically find in work_dir if not specified."""

    images: Path = "images"
    """Path the image dir, relative to work_dir."""

    ply: Optional[Path] = None
    """Path to the Metashape point export ply file. Will automatically find in work_dir if not specified."""

    psx: Optional[Path] = None
    """Optional path to the Metashape psx file. Will automatically find in work_dir if not specified.
        Used to resolve ambiguity in file name that can't be resolved by xml file alone
        (This typically happens when you import multiple images with the same file name.)"""

    def main(self) -> None:
        """Process images into a nerfstudio dataset."""

        if not os.path.exists(self.work_dir):
            raise ValueError("Work directory " + self.work_dir + " not found")
        work_dir_files = [self.work_dir / fn for fn in os.listdir(self.work_dir)]

        # Load .xml file
        if self.xml is None:
            xml_files = [fn for fn in work_dir_files if fn.suffix == ".xml"]
            if len(xml_files) == 0:
                raise ValueError("No XML file found in work_dir. Please specify using --xml")
            if len(xml_files) > 1:
                raise ValueError("Multiple XML file found in work_dir. Please specify using --xml")
            self.xml = xml_files[0]
            CONSOLE.log("Using XML file found:", self.xml)
        elif not self.xml.exists():
            self.xml = self.work_dir / self.xml
        if not self.xml.exists():
            raise ValueError(f"XML file {self.xml} doesn't exist")
        if self.xml.suffix != ".xml":
            raise ValueError(f"XML file {self.xml} must have a .xml extension")

        # Load .ply file
        if self.ply is None:
            ply_files = [fn for fn in work_dir_files if fn.suffix == ".ply" and fn.name != "sparse_pc.ply"]
            if len(ply_files) == 0:
                raise ValueError("No ply file found in work_dir (other than sparse_pc.ply). Please specify using --ply")
            if len(ply_files) > 1:
                raise ValueError("Multiple ply file found in work_dir. Please specify using --ply")
            self.ply = ply_files[0]
            CONSOLE.log("Using PLY file found:", self.ply)
        elif not self.ply.exists():
            self.ply = self.work_dir / self.ply
        if not self.ply.exists():
            raise ValueError(f"ply file {self.ply} doesn't exist")
        if self.ply.suffix != ".ply":
            raise ValueError(f"ply file {self.ply} must have a .ply extension")

        # Load .psx file
        if self.psx is None:
            psx_files = [fn for fn in work_dir_files if fn.suffix == ".psx"]
            if len(psx_files) == 0:
                raise ValueError("No psx file found in work_dir. Please specify using --psx")
            if len(psx_files) > 1:
                raise ValueError("Multiple psx file found in work_dir. Please specify using --psx")
            self.psx = psx_files[0]
            CONSOLE.log("Using PSX file found:", self.psx)
        elif not self.psx.exists():
            self.psx = self.work_dir / self.psx
            if not self.psx.exists():
                raise ValueError(f"psx file {self.psx} doesn't exist")
        camera_dict = None
        if self.psx is not None:
            if self.psx.suffix != ".psx":
                raise ValueError(f"psx file {self.psx} must have a .psx extension")
            root = ET.parse(self.psx).getroot()
            metashape_dir = self.psx.parent / root.get("path").replace("{projectname}", self.psx.name.rstrip(".psx"))
            camera_dict = find_cameras_dict(metashape_dir.parent)

        image_dir = self.work_dir / self.images
        if not image_dir.exists():
            raise ValueError(f"Image directory {image_dir} does not exist")

        summary_log = []

        image_filenames = []
        for file_path in image_dir.rglob("*"):
            if file_path.is_file():
                file_path = file_path.relative_to(self.work_dir)
                image_filenames.append(str(file_path))

        summary_log.extend(
            metashape_to_json(
                xml_filename=self.xml,
                output_dir=self.work_dir,
                ply_filename=self.ply,
                image_filenames=image_filenames,
                camera_dict=camera_dict,
            )
        )

        CONSOLE.rule("[bold green]:tada: :tada: :tada: All DONE :tada: :tada: :tada:")

        for summary in summary_log:
            CONSOLE.print(summary, justify="center")
        CONSOLE.rule()


def entrypoint():
    """Entrypoint for use with pyproject scripts."""
    tyro.extras.set_accent_color("bright_yellow")
    try:
        tyro.cli(ProcessMetashape).main()
    except (RuntimeError, ValueError) as e:
        CONSOLE.log("[bold red]" + e.args[0])


if __name__ == "__main__":
    entrypoint()
