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


import sys
import os

import zipfile
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

import bisect

import numpy as np

from pathlib import Path

from spirulae_splat.modules.camera import CameraType


class StringMatcher:
    def __init__(self, strings):
        self.orig_strings = strings
        self.sep = '\0'
        self.corpus = self.sep.join(strings) + self.sep
        
        # Map corpus indices to original string index
        self.idx_map = []
        for i, s in enumerate(strings):
            self.idx_map.extend([i] * (len(s) + 1))
            
        self.sa = self._build_suffix_array(self.corpus)
        
        # Pre-sort reversed strings for suffix queries
        self.rev_sorted = sorted((s[::-1], i) for i, s in enumerate(strings))
        self.rev_keys = [x[0] for x in self.rev_sorted]

    def _build_suffix_array(self, s):
        n = len(s)
        # Initial rank based on single characters
        sa = list(range(n))
        rank = [ord(s[i]) for i in range(n)]
        k = 1
        while k < n:
            # Key for sorting: (current rank, rank of suffix k positions ahead)
            key = lambda x: (rank[x], rank[x + k] if x + k < n else -1)
            sa.sort(key=key)
            
            # Generate new ranks based on the sorted order
            new_rank = [0] * n
            for i in range(1, n):
                new_rank[sa[i]] = new_rank[sa[i-1]] + (1 if key(sa[i]) > key(sa[i-1]) else 0)
            rank = new_rank
            if rank[sa[n-1]] == n - 1: break # Optimization: all ranks unique
            k *= 2
        return sa

    def _get_lcp_len(self, s1, s2_start):
        # Efficiently compare query string to corpus slice without creating new strings
        match_len = 0
        n1, n2 = len(s1), len(self.corpus)
        while match_len < n1 and s2_start + match_len < n2:
            if s1[match_len] != self.corpus[s2_start + match_len]:
                break
            match_len += 1
        return match_len

    def query_suffix(self, s):
        rs = s[::-1]
        idx = bisect.bisect_left(self.rev_keys, rs)
        best_len, matches = 0, []
        for i in [idx - 1, idx]:
            if 0 <= i < len(self.rev_keys):
                # Simple prefix match on reversed strings
                lcp = 0
                r_match = self.rev_keys[i]
                for c1, c2 in zip(rs, r_match):
                    if c1 != c2: break
                    lcp += 1
                if lcp > best_len:
                    best_len, matches = lcp, [self.orig_strings[self.rev_sorted[i][1]]]
                elif lcp == len(rs) and lcp == best_len: # Handle multiple exact matches
                    matches.append(self.orig_strings[self.rev_sorted[i][1]])
        return list(set(matches)), best_len

    def query_substring(self, s):
        best_len, matches = 0, set()
        # Binary search for the query string in the Suffix Array
        low, high = 0, len(self.sa) - 1
        while low <= high:
            mid = (low + high) // 2
            # Compare query to corpus starting at sa[mid]
            suffix_start = self.sa[mid]
            current_lcp = self._get_lcp_len(s, suffix_start)
            
            if current_lcp > best_len:
                best_len = current_lcp
                matches = {self.orig_strings[self.idx_map[suffix_start]]}
            elif current_lcp == best_len and current_lcp > 0:
                matches.add(self.orig_strings[self.idx_map[suffix_start]])

            # Decide which way to move in binary search
            if s > self.corpus[suffix_start : suffix_start + len(s)]:
                low = mid + 1
            else:
                high = mid - 1
        return list(matches), best_len


def find_metashape_cameras_dict(root_dir):
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
        sensor for sensor in sensors.iter("sensor")
        if sensor.get("type") in ["spherical", "cylindrical"] or sensor.find("calibration") is not None
    ]
    sensor_dict = {}
    for sensor in calibrated_sensors:
        s = {}

        sensor_type = sensor.get("type")
        if sensor_type == "frame":
            s["camera_model"] = CameraType.PERSPECTIVE.value
        elif sensor_type in ["fisheye", "equidistant_fisheye"]:
            s["camera_model"] = CameraType.EQUIDISTANT.value
        elif sensor_type in ["equisolid_fisheye"]:
            s["camera_model"] = CameraType.EQUISOLID.value
        elif sensor_type == "spherical":
            s["camera_model"] = CameraType.EQUIRECTANGULAR.value
        elif sensor_type == "cylindrical":
            s["camera_model"] = CameraType.CYLINDRICAL.value
        else:
            # Unsupported: cylindrical, RPC, etc.
            raise ValueError(f"Unsupported Metashape sensor type '{sensor_type}'")

        resolution = sensor.find("resolution")
        assert resolution is not None, "Resolution not found in Metashape xml"
        s["w"] = int(resolution.get("width"))  # type: ignore
        s["h"] = int(resolution.get("height"))  # type: ignore

        # https://github.com/facebookresearch/EyefulTower/issues/7
        # https://www.agisoft.com/pdf/metashape-pro_2_2_en.pdf Appendix D. Camera Models
        calib = sensor.find("calibration[@class!='initial']")
        if calib is None:
            assert sensor_type in ["spherical", "cylindrical"], "Only spherical sensors should have no intrinsics"
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
    component_camera_ids_dict = {}
    component_transform_dict = {}
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
                m[:3, :3] = r * s
                m[:3, 3] = t
                component_transform_dict[component.get("id")] = m
            camera_ids = []
            for comp in component.iter("camera_ids"):
                camera_ids.extend(comp.text.strip().split())
            if len(camera_ids) > 0:
                component_camera_ids_dict[component.get("id")] = camera_ids

    image_filename_matcher = StringMatcher(image_filenames)

    cameras = chunk.find("cameras")
    assert cameras is not None, "Cameras not found in Metashape xml"
    num_skipped = 0
    valid_cameras = {}
    for camera in cameras.iter("camera"):
        frame = {}
        camera_id = camera.get("id")
        if camera_dict is not None:
            image_filename = camera_dict[camera_id]
            matches = image_filename_matcher.query_suffix(image_filename)[0]
            if len(matches) != 1:
                print("WARNING: ambiguous filenames", matches, ", Skipping")
                num_skipped += 1
                continue
            frame["file_path"] = matches[0]
        else:
            camera_label = camera.get("label")
            matches = image_filename_matcher.query_substring(camera_label)[0]
            assert isinstance(camera_label, str)
            if len(matches) != 1:
                print("WARNING: ambiguous filenames", matches, ", Skipping.\n"
                            "Specify Metashape .psx file using --psx or --dataparser.psx to attempt resolving ambiguity.")
                num_skipped += 1
                continue
            frame["file_path"] = matches[0]

        sensor_id = camera.get("sensor_id")
        if sensor_id not in sensor_dict:
            # this should only happen when we have a sensor that doesn't have calibration
            print(f"Missing sensor calibration for {camera.get('label')}, Skipping")
            num_skipped += 1
            continue
        # Add all sensor parameters to this frame.
        frame.update(sensor_dict[sensor_id])

        if camera.find("transform") is None:
            print(f"Missing transforms data for {camera.get('label')}, Skipping")
            num_skipped += 1
            continue

        valid_cameras[camera_id] = (camera, frame)

    for key, value in component_camera_ids_dict.items():
        value = [id for id in value if id in valid_cameras]
        component_camera_ids_dict[key] = value

    components = [*component_camera_ids_dict.items()]
    components.sort(key=lambda x: -len(x[1]))
    components = [set(c[1]) for c in components]
    all_transforms = []
    for component_index, component_camera_ids in [*enumerate(components)][::-1]:

        frames = []

        for camera_id, (camera, frame) in valid_cameras.items():
            if camera_id not in component_camera_ids:
                continue

            transform = np.array([float(x) for x in camera.find("transform").text.split()]).reshape((4, 4))  # type: ignore

            component_id = camera.get("component_id")
            if component_id in component_transform_dict:
                transform = component_transform_dict[component_id] @ transform
                transform[:3, :3] /= np.cbrt(np.linalg.det(transform[:3, :3]))

            # Metashape camera is looking towards -Z, +X is to the right and +Y is to the top/up of the first cam
            # Rotate the scene according to nerfstudio convention
            if False: transform = transform[[2, 0, 1, 3], :]
            # Convert from Metashape's camera coordinate system (OpenCV) to ours (OpenGL)
            transform[:, 1:3] *= -1
            frame["transform_matrix"] = transform.tolist()
            frames.append(frame)

        data = {}
        data["frames"] = frames
        applied_transform = np.eye(4)[:3, :]
        if False: applied_transform = applied_transform[np.array([2, 0, 1]), :]
        data["applied_transform"] = applied_transform.tolist()

        if ply_filename is not None:
            data["ply_file_path"] = ply_filename

        all_transforms.append(data)

    summary = []

    if num_skipped == 1:
        summary.append(f"{num_skipped} image skipped.")
    if num_skipped > 1:
        summary.append(f"{num_skipped} images were skipped.")

    summary.append(f"Final dataset is {len(data['frames'])} frames.")

    return all_transforms, summary
