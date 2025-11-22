#!/usr/bin/env python3

"""Based on https://github.com/maturk/dn-splatter/blob/main/dn_splatter/export_mesh.py"""

import random
from dataclasses import dataclass
from pathlib import Path
import json
from typing import Literal, Optional, Tuple, Union, List
from torch import Tensor

import numpy as np
import open3d as o3d
import torch
import torch.nn.functional as F
import tyro
from tqdm import tqdm
from typing_extensions import Annotated

# from nerfstudio.cameras.cameras import Cameras
# from nerfstudio.models.splatfacto import SplatfactoModel
# from nerfstudio.utils.eval_utils import eval_setup
from nerfstudio.utils.rich_utils import CONSOLE

from spirulae_splat.viewer.model import SplatModel
from spirulae_splat.viewer.camera import Camera

from spirulae_splat.splat._torch_impl import depth_map, depth_inv_map
from spirulae_splat.splat import depth_to_normal


"""
Methods for extracting meshes from GS:

1) GaussiansToPoisson:
    - takes Gaussian means and predicted normals -> poisson
2) DepthAndNormalMapsPoisson
    - backproject rendered depth and normal maps -> poisson
3) LevelSetExtractor (SuGaR)
    - cast rays into scene from cameras
    - extract level sets based on gaussian density function
    - estimate normals (analytically or nearest gaussians)
    - poisson
4) Marching Cubes
    - voxelize scene and
    - compute isosurface level sets based on gaussian densities
    - run marching cubes algorithm
5) TSDF
    - voxelize scene and
    - backproject depths and integrate points into voxels for tsdf fusion
    - run marching cubes algorithm
"""



def get_camera_coords(img_size: tuple, pixel_offset: float = 0.5) -> Tensor:
    """Generates camera pixel coordinates [W,H]

    Returns:
        stacked coords [H*W,2] where [:,0] corresponds to W and [:,1] corresponds to H
    """

    # img size is (w,h)
    image_coords = torch.meshgrid(
        torch.arange(img_size[0]),
        torch.arange(img_size[1]),
        indexing="xy",  # W = u by H = v
    )
    image_coords = (
        torch.stack(image_coords, dim=-1) + pixel_offset
    )  # stored as (x, y) coordinates
    image_coords = image_coords.view(-1, 2)
    image_coords = image_coords.float()

    return image_coords


def get_means3d_backproj(
    depths: Tensor,
    fx: float,
    fy: float,
    cx: int,
    cy: int,
    img_size: tuple,
    c2w: Tensor,
    device: torch.device,
    mask: Optional[Tensor] = None,
) -> Tuple[Tensor, List]:
    """Backprojection using camera intrinsics and extrinsics

    image_coords -> (x,y,depth) -> (X, Y, depth)

    Returns:
        Tuple of (means: Tensor, image_coords: Tensor)
    """

    if depths.dim() == 3:
        depths = depths.view(-1, 1)
    elif depths.shape[-1] != 1:
        depths = depths.unsqueeze(-1).contiguous()
        depths = depths.view(-1, 1)
    if depths.dtype != torch.float:
        depths = depths.float()
        c2w = c2w.float()
    if c2w.device != device:
        c2w = c2w.to(device)

    image_coords = get_camera_coords(img_size)
    image_coords = image_coords.to(device)  # note image_coords is (H,W)

    # TODO: account for skew / radial distortion
    means3d = torch.empty(
        size=(img_size[0], img_size[1], 3), dtype=torch.float32, device=device
    ).view(-1, 3)
    means3d[:, 0] = (image_coords[:, 0] - cx) * depths[:, 0] / fx  # x
    means3d[:, 1] = (image_coords[:, 1] - cy) * depths[:, 0] / fy  # y
    means3d[:, 2] = depths[:, 0]  # z

    if mask is not None:
        if not torch.is_tensor(mask):
            mask = torch.tensor(mask, device=depths.device)
        means3d = means3d[mask]
        image_coords = image_coords[mask]

    if c2w is None:
        c2w = torch.eye((means3d.shape[0], 4, 4), device=device)

    # to world coords
    means3d = means3d @ torch.linalg.inv(c2w[..., :3, :3]) + c2w[..., :3, 3]
    return means3d, image_coords


def get_colored_points_from_depth(
    depths: Tensor,
    rgbs: Tensor,
    c2w: Tensor,
    fx: float,
    fy: float,
    cx: int,
    cy: int,
    img_size: tuple,
    mask: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor]:
    """Return colored pointclouds from depth and rgb frame and c2w. Optional masking.

    Returns:
        Tuple of (points, colors)
    """
    points, _ = get_means3d_backproj(
        depths=depths.float(),
        fx=fx,
        fy=fy,
        cx=cx,
        cy=cy,
        img_size=img_size,
        c2w=c2w.float(),
        device=depths.device,
    )
    points = points.squeeze(0)
    if mask is not None:
        if not torch.is_tensor(mask):
            mask = torch.tensor(mask, device=depths.device)
        colors = rgbs.view(-1, 3)[mask]
        points = points[mask]
    else:
        colors = rgbs.view(-1, 3)
        points = points
    return (points, colors)


def pick_indices_at_random(valid_mask, samples_per_frame):
    indices = torch.nonzero(torch.ravel(valid_mask))
    if samples_per_frame < len(indices):
        which = torch.randperm(len(indices))[:samples_per_frame]
        indices = indices[which]
    return torch.ravel(indices)


def find_depth_edges(depth_im, threshold=0.01, dilation_itr=3):
    laplacian_kernel = torch.tensor(
        [[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=depth_im.dtype, device=depth_im.device
    )
    laplacian_kernel = laplacian_kernel.unsqueeze(0).unsqueeze(0)
    depth_laplacian = (
        F.conv2d(
            (1.0 / (depth_im + 1e-6)).unsqueeze(0).unsqueeze(0).squeeze(-1),
            laplacian_kernel,
            padding=1,
        )
        .squeeze(0)
        .squeeze(0)
        .unsqueeze(-1)
    )

    edges = (depth_laplacian > threshold) * 1.0
    structure_el = laplacian_kernel * 0.0 + 1.0

    dilated_edges = edges
    for i in range(dilation_itr):
        dilated_edges = (
            F.conv2d(
                dilated_edges.unsqueeze(0).unsqueeze(0).squeeze(-1),
                structure_el,
                padding=1,
            )
            .squeeze(0)
            .squeeze(0)
            .unsqueeze(-1)
        )
    dilated_edges = (dilated_edges > 0.0) * 1.0
    return dilated_edges


@dataclass
class GSMeshExporter:
    """Base class for GS mesh exporters"""

    load_config: Path
    """Path to the work folder."""

    train_data: Path
    """Path to the training dataset."""

    output_dir: Optional[Path] = None
    """Path to the output directory."""

    max_image_size: int = 1440
    """Maximum size of rendered images"""

    depth_trunc_percentile: float = 0.95
    """Depth map truncation"""

    def get_depth_trunc(self, depth_map):
        depths = depth_map[~(depth_map > 0.999*depth_map.max())]
        if depths.numel() == 0:
            return 0.0
        return torch.quantile(depths, self.depth_trunc_percentile).item()

    def get_output_dir(self):
        if self.output_dir is not None:
            return self.output_dir
        return self.load_config.parent


@dataclass
class DepthAndNormalMapsPoisson(GSMeshExporter):
    """
    Idea: backproject depth and normal maps into 3D oriented point cloud -> Poisson
    """

    total_points: int = 4_000_000
    """Total target surface samples"""
    use_masks: bool = True
    """If dataset has masks, use these to auto crop gaussians within masked regions."""
    filter_edges_from_depth_maps: bool = False
    """Filter out edges when backprojecting from depth maps"""
    down_sample_voxel: Optional[float] = None
    """pcd down sample voxel size. Recommended value around 0.005"""
    outlier_removal: bool = True
    """Remove outliers"""
    std_ratio: float = 2.0
    """Threshold based on STD of the average distances across the point cloud to remove outliers."""
    edge_threshold: float = 0.004
    """Threshold for edge detection in depth maps (inverse depth Laplacian, resolution sensitive)"""
    edge_dilation_iterations: int = 10
    """Number of morphological dilation iterations for edge detection (swells edges)"""
    poisson_depth: int = 12
    """Poisson Octree max depth, higher values increase mesh detail"""

    @torch.no_grad()
    def main(self):
        if not self.get_output_dir().exists():
            self.get_output_dir().mkdir(parents=True)

        model = SplatModel(str(self.load_config))
        model.return_torch = True
        model.bgr = False

        applied_transform = np.array(transforms.get('applied_transform', np.eye(4)))
        if len(applied_transform) == 3:
            applied_transform = np.concatenate((applied_transform, [[0, 0, 0, 1]]))
        model.dataparser_transform = model.dataparser_transform @ np.linalg.inv(applied_transform)
        model.convert_to_input_frame()

        scales = torch.amin(model.scales, dim=-1, keepdim=True)
        model.gauss_params["opacities"] = torch.where(
            scales < torch.quantile(scales, 0.9),
            model.opacities,
            -10*torch.ones_like(model.opacities)
        )

        transforms_path = str(self.train_data / "transforms.json")
        with open(transforms_path, 'r') as fp:
            transforms = json.load(fp)

        num_frames = len(transforms["frames"])  # type: ignore
        samples_per_frame = (self.total_points + num_frames) // (num_frames)
        print("samples per frame: ", samples_per_frame)
        points = []
        normals = []
        colors = []
        for frame in tqdm(transforms["frames"]):

            for key in 'camera_model w h fl_x fl_y cx cy'.split():
                if key in transforms and key not in frame:
                    frame[key] = transforms[key]
            camera = Camera(frame)
            camera.model = "OPENCV"
            camera.distortion = (0, 0, 0, 0)

            sc = self.max_image_size / max(camera.w, camera.h)
            if sc < 1:
                camera.resize_rel(sc)

            c2w = np.array(frame['transform_matrix'])
            c2w = c2w @ np.diag([1, -1, -1, 1])

            rgb, depth_map = model.render(camera, c2w, return_depth=True)
            # depth_map = depth_inv_map(depth_map)

            H, W = camera.h, camera.w
            c2w = torch.from_numpy(c2w).float().cuda()

            if self.filter_edges_from_depth_maps:
                valid_depth = (
                    find_depth_edges(
                        depth_map,
                        threshold=self.edge_threshold,
                        dilation_itr=self.edge_dilation_iterations,
                    )
                    < 0.2
                )
            else:
                valid_depth = depth_map
            valid_mask = valid_depth
            valid_mask *= (depth_map < self.get_depth_trunc(depth_map)).float()
            valid_mask *= torch.isfinite(depth_map).float()

            indices = pick_indices_at_random(valid_mask, samples_per_frame)
            if len(indices) == 0:
                continue

            xyzs, rgbs = get_colored_points_from_depth(
                depths=depth_map,
                rgbs=rgb,
                fx=camera.fx,
                fy=camera.fy,
                cx=camera.cx,  # type: ignore
                cy=camera.cy,  # type: ignore
                img_size=(W, H),
                c2w=c2w,
                mask=indices,
            )

            # normals to OPENGL
            normal_map = depth_to_normal(depth_map, "pinhole", (camera.fx, camera.fy, camera.cx, camera.cy))

            # normals to World
            rot = c2w[:3, :3]
            normal_map = normal_map.permute(2, 0, 1).reshape(3, -1)
            normal_map = torch.nn.functional.normalize(normal_map, p=2, dim=0)
            normal_map = rot @ normal_map
            normal_map = normal_map.permute(1, 0).reshape(H, W, 3)

            normal_map = normal_map.view(-1, 3)[indices]

            points.append(xyzs)
            colors.append(rgbs)
            normals.append(normal_map)

        points = torch.cat(points, dim=0)
        colors = torch.cat(colors, dim=0)
        normals = torch.cat(normals, dim=0)

        points = points.cpu().numpy()
        normals = normals.cpu().numpy()
        colors = colors.cpu().numpy()

        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points)
        pcd.normals = o3d.utility.Vector3dVector(normals)
        pcd.colors = o3d.utility.Vector3dVector(colors)

        if self.outlier_removal:
            CONSOLE.print("Removing outliers...")
            cl, ind = pcd.remove_statistical_outlier(
                nb_neighbors=20, std_ratio=self.std_ratio
            )
            pcd = pcd.select_by_index(ind)

        CONSOLE.print("Writing point cloud...")
        o3d.io.write_point_cloud(
            str(self.get_output_dir() / "DepthAndNormalMapsPoisson_pcd.ply"), pcd
        )
        CONSOLE.print("Computing Mesh... this may take a while.")
        mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
            pcd, depth=self.poisson_depth
        )
        vertices_to_remove = densities < np.quantile(densities, 0.01)
        mesh.remove_vertices_by_mask(vertices_to_remove)
        CONSOLE.print("[bold green]:white_check_mark: Computing Mesh")

        CONSOLE.print(
            f"Saving Mesh to {str(self.get_output_dir() / 'DepthAndNormalMapsPoisson_poisson_mesh.ply')}"
        )
        o3d.io.write_triangle_mesh(
            str(self.get_output_dir() / "DepthAndNormalMapsPoisson_poisson_mesh.ply"),
            mesh,
        )


@dataclass
class TSDFFusion(GSMeshExporter):
    """
    Backproject depths and run TSDF fusion
    """

    voxel_size: float = 0.05
    """tsdf voxel size"""
    sdf_truc: float = 0.2
    """TSDF truncation"""
    total_points: int = 200000
    """Total target surface samples"""
    target_triangles: Optional[int] = None
    """Target number of triangles to simplify mesh to."""

    @torch.no_grad()
    def main(self):
        import vdbfusion

        if not self.get_output_dir().exists():
            self.get_output_dir().mkdir(parents=True)

        model = SplatModel(str(self.load_config))
        model.return_torch = True
        model.bgr = False

        transforms_path = str(self.train_data / "transforms.json")
        with open(transforms_path, 'r') as fp:
            transforms = json.load(fp)

        applied_transform = np.array(transforms.get('applied_transform', np.eye(4)))
        if len(applied_transform) == 3:
            applied_transform = np.concatenate((applied_transform, [[0, 0, 0, 1]]))
        model.dataparser_transform = model.dataparser_transform @ np.linalg.inv(applied_transform)
        model.convert_to_input_frame()

        TSDFvolume = vdbfusion.VDBVolume(
            voxel_size=self.voxel_size, sdf_trunc=self.sdf_truc, space_carving=True
        )

        num_frames = len(transforms)
        samples_per_frame = (self.total_points + num_frames) // (num_frames)
        print("samples per frame: ", samples_per_frame)
        points = []
        colors = []
        for frame in tqdm(transforms["frames"]):
            for key in 'camera_model w h fl_x fl_y cx cy'.split():
                if key in transforms and key not in frame:
                    frame[key] = transforms[key]
            camera = Camera(frame)

            c2w = np.array(frame['transform_matrix'])
            c2w = c2w @ np.diag([1, -1, -1, 1])

            rgb, depth_map = model.render(camera, c2w, return_depth=True)
            # depth_map = depth_inv_map(depth_map)

            H, W = camera.h, camera.w
            c2w = torch.from_numpy(c2w).float().cuda()

            depth_trunc = self.get_depth_trunc(depth_map)
            depth_map[~(depth_map < depth_trunc)] = 0

            xyzs, rgbs = get_colored_points_from_depth(
                depths=depth_map,
                rgbs=rgb,
                fx=camera.fx,
                fy=camera.fy,
                cx=camera.cx,
                cy=camera.cy,
                img_size=(W, H),
                c2w=c2w,
            )
            points.append(xyzs)
            colors.append(rgbs)
            TSDFvolume.integrate(
                xyzs.double().cpu().numpy(),
                extrinsic=c2w[:3, 3].double().cpu().numpy(),
            )

        vertices, faces = TSDFvolume.extract_triangle_mesh(min_weight=5)

        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(vertices)
        mesh.triangles = o3d.utility.Vector3iVector(faces)
        mesh.compute_vertex_normals()
        colors = torch.cat(colors, dim=0)
        colors = colors.cpu().numpy()
        mesh.vertex_colors = o3d.utility.Vector3dVector(colors)

        # simplify mesh
        if self.target_triangles is not None:
            mesh = mesh.simplify_quadric_decimation(self.target_triangles)

        o3d.io.write_triangle_mesh(
            str(self.get_output_dir() / "TSDFfusion_mesh.ply"),
            mesh,
        )
        CONSOLE.print(
            f"Finished computing mesh: {str(self.get_output_dir() / 'TSDFfusion.ply')}"
        )


@dataclass
class Open3DTSDFFusion(GSMeshExporter):
    """
    Backproject depths and run TSDF fusion
    """

    voxel_size: float = 0.0025
    """tsdf voxel size"""
    sdf_truc: float = 0.02
    """TSDF truncation"""

    @torch.no_grad()
    def main(self):
        import open3d as o3d

        model = SplatModel(str(self.load_config))
        model_scale = model.dataparser_scale
        model.return_torch = True
        model.bgr = False

        if False:
            scales = torch.amin(model.scales, dim=-1, keepdim=True)
            model.gauss_params["opacities"] = torch.where(
                scales < torch.quantile(scales, 0.9),
                model.opacities,
                -10*torch.ones_like(model.opacities)
            )

        transforms_path = str(self.train_data / "transforms.json")
        with open(transforms_path, 'r') as fp:
            transforms = json.load(fp)

        applied_transform = np.array(transforms.get('applied_transform', np.eye(4)))
        if len(applied_transform) == 3:
            applied_transform = np.concatenate((applied_transform, [[0, 0, 0, 1]]))
        model.dataparser_transform = model.dataparser_transform @ np.linalg.inv(applied_transform)
        model.convert_to_input_frame()

        volume = o3d.pipelines.integration.ScalableTSDFVolume(
            voxel_length=self.voxel_size/model_scale,
            sdf_trunc=self.sdf_truc/model_scale,
            color_type=o3d.pipelines.integration.TSDFVolumeColorType.RGB8,
        )

        for frame in tqdm(transforms["frames"]):

            for key in 'camera_model w h fl_x fl_y cx cy'.split():
                if key in transforms and key not in frame:
                    frame[key] = transforms[key]
            camera = Camera(frame)
            camera.model = "OPENCV"
            camera.distortion = (0, 0, 0, 0)

            sc = self.max_image_size / max(camera.w, camera.h)
            if sc < 1:
                camera.resize_rel(sc)

            c2w = np.array(frame['transform_matrix'])
            c2w = c2w @ np.diag([1, -1, -1, 1])

            rgb_map, depth_map = model.render(camera, c2w, return_depth=True)

            H, W = camera.h, camera.w
            c2w = torch.from_numpy(c2w).float().cuda()

            intrinsic = o3d.camera.PinholeCameraIntrinsic(
                width=W,
                height=H,
                fx=camera.fx,
                fy=camera.fy,
                cx=camera.cx,
                cy=camera.cy,
            )

            depth_trunc = self.get_depth_trunc(depth_map)

            rgbd = o3d.geometry.RGBDImage.create_from_color_and_depth(
                o3d.geometry.Image(
                    np.asarray(
                        rgb_map.cpu().numpy() * 255,
                        order="C",
                        dtype=np.uint8,
                    )
                ),
                o3d.geometry.Image(
                    np.asarray(depth_map.squeeze(-1).cpu().numpy(), order="C")
                ),
                depth_trunc=depth_trunc,
                convert_rgb_to_intensity=False,
                depth_scale=1.0,
            )

            volume.integrate(
                rgbd,
                intrinsic=intrinsic,
                extrinsic=np.linalg.inv(c2w.cpu().numpy()),
            )

        mesh = volume.extract_triangle_mesh()

        mesh_0 = mesh
        with o3d.utility.VerbosityContextManager(
            o3d.utility.VerbosityLevel.Debug
        ) as cm:
            (
                triangle_clusters,
                cluster_n_triangles,
                cluster_area,
            ) = mesh_0.cluster_connected_triangles()

        triangle_clusters = np.asarray(triangle_clusters)
        cluster_n_triangles = np.asarray(cluster_n_triangles)
        cluster_area = np.asarray(cluster_area)
        n_filter = min(50, len(cluster_n_triangles))
        n_cluster = np.sort(cluster_n_triangles.copy())[-n_filter]
        n_cluster = max(n_cluster, n_filter)  # filter meshes smaller than 50
        triangles_to_remove = cluster_n_triangles[triangle_clusters] < n_cluster
        mesh_0.remove_triangles_by_mask(triangles_to_remove)
        mesh_0.remove_unreferenced_vertices()
        mesh_0.remove_degenerate_triangles()

        if not self.get_output_dir().exists():
            self.get_output_dir().mkdir(parents=True)
        o3d.io.write_triangle_mesh(
            str(self.get_output_dir() / "Open3dTSDFfusion_mesh.ply"),
            mesh,
        )
        CONSOLE.print(
            f"Finished computing mesh: {str(self.get_output_dir() / 'Open3dTSDFfusion.ply')}"
        )


Commands = tyro.conf.FlagConversionOff[
    Union[
        Annotated[TSDFFusion, tyro.conf.subcommand(name="tsdf")],
        Annotated[Open3DTSDFFusion, tyro.conf.subcommand(name="o3dtsdf")],
        Annotated[DepthAndNormalMapsPoisson, tyro.conf.subcommand(name="dn")],
    ]
]


def entrypoint():
    """Entrypoint for use with pyproject scripts."""
    tyro.extras.set_accent_color("bright_yellow")
    tyro.cli(Commands).main()


if __name__ == "__main__":
    # tyro.cli(DepthAndNormalMapsPoisson).main()
    # tyro.cli(TSDFFusion).main()
    tyro.cli(Open3DTSDFFusion).main()
