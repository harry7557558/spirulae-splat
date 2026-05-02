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

"""Script to convert Metashape data into Nerfstudio format."""

import json
from pathlib import Path
from dataclasses import dataclass
import tyro

import os

import open3d as o3d

from typing import Optional

from spirulae_splat.modules.metashape_utils import *


@dataclass
class ProcessMetashape():
    """Process Metashape data into a nerfstudio dataset.

    This script assumes that cameras have been aligned using Metashape. After alignment, it is necessary to export the
    camera poses as a `.xml` and point clouds as a `.ply` file under the dataset folder. For point cloud, Make sure to
    export the data in non-binary format and exclude the normals. Optionally, provide the Metashape `.psx` file to 
    resolve file name ambiguity (Tested with Metashape Standard Version 2.2.2).

    """

    work_dir: Path
    """Path to the dataset directory."""

    xml: Optional[Path] = None
    """Path to the Metashape xml file. Will automatically find in work_dir if not specified."""

    images: Path = Path("images")
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
            xml_files = [fn for fn in work_dir_files if fn.suffix.lower() == ".xml"]
            if len(xml_files) == 0:
                raise ValueError("No XML file found in work_dir. Please specify using --xml")
            if len(xml_files) > 1:
                raise ValueError("Multiple XML file found in work_dir. Please specify using --xml")
            self.xml = xml_files[0]
            print("Using XML file found:", self.xml)
        elif not self.xml.exists():
            self.xml = self.work_dir / self.xml
        if not self.xml.exists():
            raise ValueError(f"XML file {self.xml} doesn't exist")
        if self.xml.suffix.lower() != ".xml":
            raise ValueError(f"XML file {self.xml} must have a .xml extension")

        # Load .ply file
        if self.ply is None:
            ply_files = [fn for fn in work_dir_files if fn.suffix.lower() == ".ply" and fn.name != "sparse_pc.ply"]
            if len(ply_files) == 0:
                raise ValueError("No ply file found in work_dir (other than sparse_pc.ply). Please specify using --ply")
            if len(ply_files) > 1:
                raise ValueError("Multiple ply file found in work_dir. Please specify using --ply")
            self.ply = ply_files[0]
            print("Using PLY file found:", self.ply)
        elif not self.ply.exists():
            self.ply = self.work_dir / self.ply
        if not self.ply.exists():
            raise ValueError(f"ply file {self.ply} doesn't exist")
        if self.ply.suffix.lower() != ".ply":
            raise ValueError(f"ply file {self.ply} must have a .ply extension")

        # Load .psx file
        if self.psx is None and False:
            psx_files = [fn for fn in work_dir_files if fn.suffix.lower() == ".psx"]
            if len(psx_files) > 1:
                raise ValueError("Multiple psx file found in work_dir. Please specify using --psx")
            if len(psx_files) != 0:
                self.psx = psx_files[0]
                print("Using PSX file found:", self.psx)
        if self.psx is not None and not self.psx.exists():
            self.psx = self.work_dir / self.psx
            if not self.psx.exists():
                raise ValueError(f"psx file {self.psx} doesn't exist")
        camera_dict = None
        if self.psx is not None:
            if self.psx.suffix.lower() != ".psx":
                raise ValueError(f"psx file {self.psx} must have a .psx extension")
            root = ET.parse(self.psx).getroot()
            metashape_dir = self.psx.parent / root.get("path").replace("{projectname}", self.psx.name.rstrip(".psx"))
            camera_dict = find_metashape_cameras_dict(metashape_dir.parent)

        image_dir = self.work_dir / self.images
        if not image_dir.exists():
            raise ValueError(f"Image directory `{image_dir}` does not exist")

        summary_log = []

        image_filenames = []
        for file_path in image_dir.rglob("*"):
            if file_path.is_file():
                file_path = file_path.relative_to(self.work_dir)
                image_filenames.append(str(file_path))

        transforms, summary = metashape_to_json(
            xml_filename=self.xml,
            ply_filename=self.ply,
            image_filenames=image_filenames,
            camera_dict=camera_dict,
        )
        summary_log.extend(summary)

        for component_index, data in enumerate(transforms):
            ply_filename = data.get("ply_file_path", None)
            applied_transform = np.array(data['applied_transform'])
            if ply_filename is not None:
                assert ply_filename.name != "sparse_pc.ply", "PLY file name must not be sparse_pc.ply"
                assert ply_filename.exists()
                pc = o3d.io.read_point_cloud(str(ply_filename))
                points3D = np.asarray(pc.points)
                points3D = np.einsum("ij,bj->bi", applied_transform[:3, :3], points3D) + applied_transform[:3, 3]
                pc.points = o3d.utility.Vector3dVector(points3D)
                o3d.io.write_point_cloud(str(self.work_dir / "sparse_pc.ply"), pc)
                data["ply_file_path"] = "sparse_pc.ply"
                summary.append(f"Imported {ply_filename} as starting points")

            transforms_filename = "transforms.json" if component_index == 0 else f"transforms_{component_index}.json"
            with open(self.work_dir / transforms_filename, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=4)

        print("All DONE.")

        for summary in summary_log:
            print(summary)


def entrypoint():
    """Entrypoint for use with pyproject scripts."""
    tyro.extras.set_accent_color("bright_yellow")
    try:
        tyro.cli(ProcessMetashape).main()
    except (RuntimeError, ValueError) as e:
        print(e.args[0])


if __name__ == "__main__":
    entrypoint()
