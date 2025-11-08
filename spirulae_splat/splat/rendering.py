import torch
import torch.nn.functional as F
from torch import Tensor
from typing import Dict, Optional, Tuple, Literal
import math


from gsplat.cuda._wrapper import (
    RollingShutterType,
    # FThetaCameraDistortionParameters,
    # FThetaPolynomialType,
    # fully_fused_projection,
    fully_fused_projection_with_ut,
    isect_offset_encode,
    isect_tiles,
    # rasterize_to_pixels,
    # rasterize_to_pixels_eval3d,
    # spherical_harmonics,
)
from spirulae_splat.splat.cuda._wrapper import (
    intersect_splat_tile,
    fully_fused_projection,
    fully_fused_projection_hetero,
    rasterize_to_pixels,
    spherical_harmonics,
)
from gsplat.distributed import (
    all_gather_int32,
    all_gather_tensor_list,
    all_to_all_int32,
    all_to_all_tensor_list,
)

from .utils import depth_to_normal


# from gsplat

def rasterization(
    primitive: Literal["3dgs", "mip", "opaque_triangle"],
    splat_params: tuple[Tensor],  # means, quats, scales, opacities
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    width: int,
    height: int,
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    eps2d: float = 0.3,
    packed: bool = True,
    use_bvh: bool = False,
    tile_size: int = 16,
    backgrounds: Optional[Tensor] = None,
    render_mode: Literal["RGB", "D", "ED", "RGB+D", "RGB+ED"] = "RGB",
    sparse_grad: bool = False,
    absgrad: bool = False,
    # rasterize_mode: Literal["classic", "antialiased"] = "classic",
    channel_chunk: int = 32,
    distributed: bool = False,
    # camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    camera_model: Literal["pinhole", "ortho", "fisheye"] = "pinhole",
    segmented: bool = False,
    with_ut: bool = False,
    with_eval3d: bool = False,
    # distortion
    radial_coeffs: Optional[Tensor] = None,  # [..., C, 6] or [..., C, 4]
    tangential_coeffs: Optional[Tensor] = None,  # [..., C, 2]
    thin_prism_coeffs: Optional[Tensor] = None,  # [..., C, 4]
    # ftheta_coeffs: Optional[FThetaCameraDistortionParameters] = None,
    # rolling shutter
    rolling_shutter: RollingShutterType = RollingShutterType.GLOBAL,
    viewmats_rs: Optional[Tensor] = None,  # [..., C, 4, 4]
    world_rank: int = 0,
    world_size: int = 1
) -> Tuple[Tensor, Tensor, Dict]:
    """Rasterize a set of 3D Gaussians (N) to a batch of image planes (C).

    This function provides a handful features for 3D Gaussian rasterization, which
    we detail in the following notes. A complete profiling of the these features
    can be found in the :ref:`profiling` page.

    .. note::
        **Multi-GPU Distributed Rasterization**: This function can be used in a multi-GPU
        distributed scenario by setting `distributed` to True. When `distributed` is True,
        a subset of total Gaussians could be passed into this function in each rank, and
        the function will collaboratively render a set of images using Gaussians from all ranks. Note
        to achieve balanced computation, it is recommended (not enforced) to have similar number of
        Gaussians in each rank. But we do enforce that the number of cameras to be rendered
        in each rank is the same. The function will return the rendered images
        corresponds to the input cameras in each rank, and allows for gradients to flow back to the
        Gaussians living in other ranks. For the details, please refer to the paper
        `On Scaling Up 3D Gaussian Splatting Training <https://arxiv.org/abs/2406.18533>`_.

    .. note::
        **Batch Rasterization**: This function allows for rasterizing a set of 3D Gaussians
        to a batch of images in one go, by simplly providing the batched `viewmats` and `Ks`.

    .. note::
        **Depth Rendering**: This function supports colors or/and depths via `render_mode`.
        The supported modes are "RGB", "D", "ED", "RGB+D", and "RGB+ED". "RGB" renders the
        colored image that respects the `colors` argument. "D" renders the accumulated z-depth
        :math:`\\sum_i w_i z_i`. "ED" renders the expected z-depth
        :math:`\\frac{\\sum_i w_i z_i}{\\sum_i w_i}`. "RGB+D" and "RGB+ED" render both
        the colored image and the depth, in which the depth is the last channel of the output.

    .. note::
        **Memory-Speed Trade-off**: The `packed` argument provides a trade-off between
        memory footprint and runtime. If `packed` is True, the intermediate results are
        packed into sparse tensors, which is more memory efficient but might be slightly
        slower. This is especially helpful when the scene is large and each camera sees only
        a small portion of the scene. If `packed` is False, the intermediate results are
        with shape [..., C, N, ...], which is faster but might consume more memory.

    .. note::
        **Sparse Gradients**: If `sparse_grad` is True, the gradients for {means, quats, scales}
        will be stored in a `COO sparse layout <https://pytorch.org/docs/stable/generated/torch.sparse_coo_tensor.html>`_.
        This can be helpful for saving memory
        for training when the scene is large and each iteration only activates a small portion
        of the Gaussians. Usually a sparse optimizer is required to work with sparse gradients,
        such as `torch.optim.SparseAdam <https://pytorch.org/docs/stable/generated/torch.optim.SparseAdam.html#sparseadam>`_.
        This argument is only effective when `packed` is True.

    # .. note::
    #     **Antialiased Rendering**: If `rasterize_mode` is "antialiased", the function will
    #     apply a view-dependent compensation factor
    #     :math:`\\rho=\\sqrt{\\frac{Det(\\Sigma)}{Det(\\Sigma+ \\epsilon I)}}` to Gaussian
    #     opacities, where :math:`\\Sigma` is the projected 2D covariance matrix and :math:`\\epsilon`
    #     is the `eps2d`. This will make the rendered image more antialiased, as proposed in
    #     the paper `Mip-Splatting: Alias-free 3D Gaussian Splatting <https://arxiv.org/pdf/2311.16493>`_.

    .. note::
        **AbsGrad**: If `absgrad` is True, the absolute gradients of the projected
        2D means will be computed during the backward pass, which could be accessed by
        `meta["means2d"].absgrad`. This is an implementation of the paper
        `AbsGS: Recovering Fine Details for 3D Gaussian Splatting <https://arxiv.org/abs/2404.10484>`_,
        which is shown to be more effective for splitting Gaussians during training.

    .. note::
        **Camera Distortion and Rolling Shutter**: The function supports rendering with opencv
        distortion formula for pinhole and fisheye cameras (`radial_coeffs`, `tangential_coeffs`, `thin_prism_coeffs`).
        It also supports rolling shutter rendering with the `rolling_shutter` argument. We take
        reference from the paper `3DGUT: Enabling Distorted Cameras and Secondary Rays in Gaussian Splatting
        <https://arxiv.org/abs/2412.12507>`_.

    .. warning::
        This function is currently not differentiable w.r.t. the camera intrinsics `Ks`.

    Args:
        means: The 3D centers of the Gaussians. [..., N, 3]
        quats: The quaternions of the Gaussians (wxyz convension). It's not required to be normalized. [..., N, 4]
        scales: The scales of the Gaussians. [..., N, 3]
        opacities: The opacities of the Gaussians. [..., N]
        colors: The colors of the Gaussians. [..., (C,) N, D] or [..., (C,) N, K, 3] for SH coefficients.
        viewmats: The world-to-cam transformation of the cameras. [..., C, 4, 4]
        Ks: The camera intrinsics. [..., C, 3, 3]
        width: The width of the image.
        height: The height of the image.
        near_plane: The near plane for clipping. Default is 0.01.
        far_plane: The far plane for clipping. Default is 1e10.
        eps2d: An epsilon added to the egienvalues of projected 2D covariance matrices.
            This will prevents the projected GS to be too small. For example eps2d=0.3
            leads to minimal 3 pixel unit. Default is 0.3.
        packed: Whether to use packed mode which is more memory efficient but might or
            might not be as fast. Default is True.
        tile_size: The size of the tiles for rasterization. Default is 16.
            (Note: other values are not tested)
        backgrounds: The background colors. [..., C, D]. Default is None.
        render_mode: The rendering mode. Supported modes are "RGB", "D", "ED", "RGB+D",
            and "RGB+ED". "RGB" renders the colored image, "D" renders the accumulated depth, and
            "ED" renders the expected depth. Default is "RGB".
        sparse_grad: If true, the gradients for {means, quats, scales} will be stored in
            a COO sparse layout. This can be helpful for saving memory. Default is False.
        absgrad: If true, the absolute gradients of the projected 2D means
            will be computed during the backward pass, which could be accessed by
            `meta["means2d"].absgrad`. Default is False.
        # rasterize_mode: The rasterization mode. Supported modes are "classic" and
        #     "antialiased". Default is "classic".
        channel_chunk: The number of channels to render in one go. Default is 32.
            If the required rendering channels are larger than this value, the rendering
            will be done looply in chunks.
        distributed: Whether to use distributed rendering. Default is False. If True,
            The input Gaussians are expected to be a subset of scene in each rank, and
            the function will collaboratively render the images for all ranks.
        camera_model: The camera model to use. Supported models are "pinhole", "ortho",
            "fisheye", and "ftheta". Default is "pinhole".
        segmented: Whether to use segmented radix sort. Default is False.
            Segmented radix sort performs sorting in segments, which is more efficient for the sorting operation itself.
            However, since it requires offset indices as input, additional global memory access is needed, which results
            in slower overall performance in most use cases.
        with_ut: Whether to use Unscented Transform (UT) for projection. Default is False.
        with_eval3d: Whether to calculate Gaussian response in 3D world space, instead
            of 2D image space. Default is False.
        radial_coeffs: Opencv pinhole/fisheye radial distortion coefficients. Default is None.
            For pinhole camera, the shape should be [..., C, 6]. For fisheye camera, the shape
            should be [..., C, 4].
        tangential_coeffs: Opencv pinhole tangential distortion coefficients. Default is None.
            The shape should be [..., C, 2] if provided.
        thin_prism_coeffs: Opencv pinhole thin prism distortion coefficients. Default is None.
            The shape should be [..., C, 4] if provided.
        # ftheta_coeffs: F-Theta camera distortion coefficients shared for all cameras.
        #     Default is None. See `FThetaCameraDistortionParameters` for details.
        rolling_shutter: The rolling shutter type. Default `RollingShutterType.GLOBAL` means
            global shutter.
        viewmats_rs: The second viewmat when rolling shutter is used. Default is None.

    Returns:
        A tuple:

        **render_colors**: The rendered colors. [..., C, height, width, X].
        X depends on the `render_mode` and input `colors`. If `render_mode` is "RGB",
        X is D; if `render_mode` is "D" or "ED", X is 1; if `render_mode` is "RGB+D" or
        "RGB+ED", X is D+1.

        **render_alphas**: The rendered alphas. [..., C, height, width, 1].

        **meta**: A dictionary of intermediate results of the rasterization.

    Examples:

    .. code-block:: python

        >>> # define Gaussians
        >>> means = torch.randn((100, 3), device=device)
        >>> quats = torch.randn((100, 4), device=device)
        >>> scales = torch.rand((100, 3), device=device) * 0.1
        >>> colors = torch.rand((100, 3), device=device)
        >>> opacities = torch.rand((100,), device=device)
        >>> # define cameras
        >>> viewmats = torch.eye(4, device=device)[None, :, :]
        >>> Ks = torch.tensor([
        >>>    [300., 0., 150.], [0., 300., 100.], [0., 0., 1.]], device=device)[None, :, :]
        >>> width, height = 300, 200
        >>> # render
        >>> colors, alphas, meta = rasterization(
        >>>    means, quats, scales, opacities, colors, viewmats, Ks, width, height
        >>> )
        >>> print (colors.shape, alphas.shape)
        torch.Size([1, 200, 300, 3]) torch.Size([1, 200, 300, 1])
        >>> print (meta.keys())
        dict_keys(['camera_ids', 'gaussian_ids', 'radii', 'means2d', 'depths', 'conics',
        'opacities', 'tile_width', 'tile_height', 'tiles_per_gauss', 'isect_ids',
        'flatten_ids', 'isect_offsets', 'width', 'height', 'tile_size'])

    """
    meta = {}

    if primitive in ["3dgs", "mip", "opaque_triangle"]:
        if primitive in ["3dgs", "mip"]:
            assert len(splat_params) == 6, "3DGS requires 6 params (means, quats, scales, opacities, color, sh)"
            means, quats, scales, opacities, features_dc, features_sh = splat_params
        else:
            assert len(splat_params) == 7, "Opaque triangle requires 4 params (means, quats, scales, opacities, color, ch, sh)"
            means, quats, scales, opacities, features_dc, features_sh, features_ch = splat_params
        assert len(means.shape) >= 2, "means must have at least 2 dimensions"
        batch_dims = means.shape[:-2]
        N = means.shape[-2]
        device = means.device
        assert means.shape == batch_dims + (N, 3), means.shape
        assert quats.shape == batch_dims + (N, 4), quats.shape
        assert scales.shape == batch_dims + (N, 3), scales.shape
        if primitive in ["3dgs", "mip"]:
            assert opacities.shape == batch_dims + (N,), opacities.shape
            assert features_dc.shape == batch_dims + (N, 3), features_dc.shape
            assert features_sh.shape == batch_dims + (N, 3, 3) or \
                features_sh.shape == batch_dims + (N, 8, 3) or \
                features_sh.shape == batch_dims + (N, 15, 3), features_sh.shape
        else:
            assert opacities.shape == batch_dims + (N, 2), opacities.shape
            assert features_ch.shape == batch_dims + (N, 2, 3), features_ch.shape
    else:
        assert False, f"Invalid primitive ({primitive})"
    num_batch_dims = len(batch_dims)
    B = math.prod(batch_dims)
    C = viewmats.shape[-3]
    I = B * C
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    assert render_mode in ["RGB", "D", "ED", "RGB+D", "RGB+ED", "RGB+D+N", "RGB+ED+N"], render_mode

    def reshape_view(C: int, world_view: torch.Tensor, N_world: list) -> torch.Tensor:
        view_list = list(
            map(
                lambda x: x.split(int(x.shape[0] / C), dim=0),
                world_view.split([C * N_i for N_i in N_world], dim=0),
            )
        )
        return torch.stack([torch.cat(l, dim=0) for l in zip(*view_list)], dim=0)

    if absgrad:
        assert not distributed, "AbsGrad is not supported in distributed mode."

    if (
        # radial_coeffs is not None
        # or tangential_coeffs is not None
        thin_prism_coeffs is not None
        # or ftheta_coeffs is not None
        or rolling_shutter != RollingShutterType.GLOBAL
    ):
        assert (
            with_ut
        ), "Distortion and rolling shutter are only supported with `with_ut=True`."

    if rolling_shutter != RollingShutterType.GLOBAL:
        assert (
            viewmats_rs is not None
        ), "Rolling shutter requires to provide viewmats_rs."
    else:
        assert (
            viewmats_rs is None
        ), "viewmats_rs should be None for global rolling shutter."

    if with_ut or with_eval3d:
        assert (quats is not None) and (
            scales is not None
        ), "UT and eval3d requires to provide quats and scales."
        assert packed is False, "Packed mode is not supported with UT."
        assert sparse_grad is False, "Sparse grad is not supported with UT."

    if radial_coeffs is not None:
        radial_coeffs = radial_coeffs.contiguous()
    if tangential_coeffs is not None:
        tangential_coeffs = tangential_coeffs.contiguous()
    if thin_prism_coeffs is not None:
        thin_prism_coeffs = thin_prism_coeffs.contiguous()

    # Implement the multi-GPU strategy proposed in
    # `On Scaling Up 3D Gaussian Splatting Training <https://arxiv.org/abs/2406.18533>`.
    #
    # If in distributed mode, we distribute the projection computation over Gaussians
    # and the rasterize computation over cameras. So first we gather the cameras
    # from all ranks for projection.
    if distributed:
        raise NotImplementedError()
        assert batch_dims == (), "Distributed mode does not support batch dimensions"
        # world_rank = torch.distributed.get_rank()
        # world_size = torch.distributed.get_world_size()

        # Gather the number of Gaussians in each rank.
        N_world = all_gather_int32(world_size, N, device=device)

        # Enforce that the number of cameras is the same across all ranks.
        C_world = [C] * world_size
        viewmats, Ks = all_gather_tensor_list(world_size, [viewmats, Ks])
        if viewmats_rs is not None:
            (viewmats_rs,) = all_gather_tensor_list(world_size, [viewmats_rs])

        # Silently change C from local #Cameras to global #Cameras.
        C = len(viewmats)

    if use_bvh:
        raise NotImplementedError()
        assert not with_ut, "Not implemented"
        assert packed, "BVH must be packed"
        assert B == 1, "Not support batching"
        def dump_all():
            print(width, height)
            def dump(name, tensor):
                shape = '_'.join(map(str, tensor.shape))
                name = f"{name}_{shape}.bin"
                print("dump", name)
                tensor.contiguous().detach().cpu().numpy().tofile(f"/home/harry/temp/{name}")
            dump('means', means)
            dump('scales', scales)
            dump('opacities', opacities)
            dump('quats', quats)
            dump('viewmats', viewmats)
            dump('Ks', Ks)
            # exit(0)
        from time import perf_counter
        torch.cuda.synchronize()
        time0 = perf_counter()
        # dump_all()
        with torch.no_grad():
            intersection_count_map, intersection_splat_id = intersect_splat_tile(
                means.contiguous(),
                torch.exp(scales).contiguous(),
                torch.sigmoid(opacities).contiguous(),
                quats.contiguous(),
                width,
                height,
                viewmats.contiguous(),
                Ks.contiguous(),
            )
        torch.cuda.synchronize()
        time1 = perf_counter()
        # print(1e3*(time1-time0), 'ms')
        # if 1e3*(time1-time0) > 100:
        #     dump_all()
        #     exit(0)
        proj_results = fully_fused_projection_hetero(
            means,
            quats,
            torch.exp(scales),
            viewmats,
            Ks,
            width,
            height,
            intersection_count_map,
            intersection_splat_id,
            eps2d=eps2d,
            near_plane=near_plane,
            far_plane=far_plane,
            sparse_grad=sparse_grad,
            calc_compensations=(rasterize_mode == "antialiased"),
            camera_model=camera_model,
            opacities=torch.sigmoid(opacities),  # use opacities to compute a tigher bound for radii.
        )

    else:
        # Project splats to 2D
        proj_results = fully_fused_projection(
            primitive,
            with_eval3d,
            splat_params,
            viewmats,
            Ks,
            width,
            height,
            packed=packed,
            near_plane=near_plane,
            far_plane=far_plane,
            sparse_grad=sparse_grad,
            camera_model=camera_model,
            radial_coeffs=radial_coeffs,
            tangential_coeffs=tangential_coeffs,
            thin_prism_coeffs=thin_prism_coeffs,
        )
        # import gsplat.cuda._wrapper
        # proj_results_1 = gsplat.cuda._wrapper.fully_fused_projection(
        #     means,
        #     None,
        #     quats,
        #     torch.exp(scales),
        #     viewmats,
        #     Ks,
        #     width,
        #     height,
        #     eps2d=eps2d,
        #     packed=packed,
        #     near_plane=near_plane,
        #     far_plane=far_plane,
        #     sparse_grad=sparse_grad,
        #     calc_compensations=(primitive == "mip"),
        #     camera_model=camera_model,
        #     opacities=torch.sigmoid(opacities),  # use opacities to compute a tigher bound for radii.
        # )
        # for x0, x1 in zip(proj_results, proj_results_1):
        #     print(x0.shape, x0.dtype, x1.shape, x1.dtype)
        #     print(x0)
        #     print(x1)
        #     import matplotlib.pyplot as plt
        #     plt.scatter(x0.flatten().detach().cpu().numpy(), x1.flatten().detach().cpu().numpy())
        #     # plt.xscale('log')
        #     # plt.yscale('log')
        #     plt.show()

    if use_bvh:
        (
            camera_ids,
            gaussian_ids,
            radii,
            means2d,
            depths,
            conics,
            proj_opacities,
        ) = proj_results
        batch_ids, image_ids = 0, camera_ids
        proj_opacities = proj_opacities.view(B, N)[batch_ids, gaussian_ids]  # [nnz]
    elif packed:
        # The results are packed into shape [nnz, ...]. All elements are valid.
        (
            batch_ids,
            camera_ids,
            gaussian_ids,
            radii,
            means2d,
            depths,
            conics,
            proj_opacities,
        ) = proj_results
        image_ids = batch_ids * C + camera_ids
    else:
        # The results are with shape [..., C, N, ...]. Only the elements with radii > 0 are valid.
        aabb_xyxy, depths, proj_splats = proj_results
        batch_ids, camera_ids, gaussian_ids = None, None, None
        image_ids = None
    means2d = (aabb_xyxy[..., 2:] + aabb_xyxy[..., :2]).float() / 2
    radii = (aabb_xyxy[..., 2:] - aabb_xyxy[..., :2] + 1) // 2
    depths = torch.exp(depths)

    meta.update(
        {
            # global batch and camera ids
            "batch_ids": batch_ids,
            "camera_ids": camera_ids,
            # local gaussian_ids
            "gaussian_ids": gaussian_ids,
            "radii": radii,
            "depths": depths,
            # "normals": normals,
            "means2d": proj_splats[0],  # with grad
        }
    )
    # if heterogeneous:
    #     meta.update({
    #         "intersection_count": intersection_count_map[1:]-intersection_count_map[:-1]
    #     })

    # If in distributed mode, we need to scatter the GSs to the destination ranks, based
    # on which cameras they are visible to, which we already figured out in the projection
    # stage.
    if distributed:
        if packed:
            raise NotImplementedError()
            # count how many elements need to be sent to each rank
            cnts = torch.bincount(camera_ids, minlength=C)  # all cameras
            cnts = cnts.split(C_world, dim=0)
            cnts = [cuts.sum() for cuts in cnts]

            # all to all communication across all ranks. After this step, each rank
            # would have all the necessary GSs to render its own images.
            collected_splits = all_to_all_int32(world_size, cnts, device=device)
            (radii,) = all_to_all_tensor_list(
                world_size, [radii], cnts, output_splits=collected_splits
            )
            (means2d, depths, conics, proj_opacities, colors) = all_to_all_tensor_list(
                world_size,
                [means2d, depths, conics, proj_opacities, colors],
                cnts,
                output_splits=collected_splits,
            )

            # before sending the data, we should turn the camera_ids from global to local.
            # i.e. the camera_ids produced by the projection stage are over all cameras world-wide,
            # so we need to turn them into camera_ids that are local to each rank.
            offsets = torch.tensor(
                [0] + C_world[:-1], device=camera_ids.device, dtype=camera_ids.dtype
            )
            offsets = torch.cumsum(offsets, dim=0)
            offsets = offsets.repeat_interleave(torch.stack(cnts))
            camera_ids = camera_ids - offsets

            # and turn gaussian ids from local to global.
            offsets = torch.tensor(
                [0] + N_world[:-1],
                device=gaussian_ids.device,
                dtype=gaussian_ids.dtype,
            )
            offsets = torch.cumsum(offsets, dim=0)
            offsets = offsets.repeat_interleave(torch.stack(cnts))
            gaussian_ids = gaussian_ids + offsets

            # all to all communication across all ranks.
            (camera_ids, gaussian_ids) = all_to_all_tensor_list(
                world_size,
                [camera_ids, gaussian_ids],
                cnts,
                output_splits=collected_splits,
            )

            # Silently change C from global #Cameras to local #Cameras.
            C = C_world[world_rank]

        else:
            # Silently change C from global #Cameras to local #Cameras.
            C = C_world[world_rank]

            # all to all communication across all ranks. After this step, each rank
            # would have all the necessary GSs to render its own images.
            (radii,) = all_to_all_tensor_list(
                world_size,
                [radii.flatten(0, 1)],
                splits=[C_i * N for C_i in C_world],
                output_splits=[C * N_i for N_i in N_world],
            )
            radii = reshape_view(C, radii, N_world)

            is_batched = [(comp.shape[0] == C) for comp in proj_splats]  # TODO: splat count == C?
            proj_splats = all_to_all_tensor_list(
                world_size,
                [(comp.flatten(0, 1) if comp.shape[0] == C else comp)  # TODO
                 for comp in proj_splats],
                splits=[C_i * N for C_i in C_world],
                output_splits=[C * N_i for N_i in N_world],
            )
            proj_splats = [(reshape_view(C, comp, N_world) if ib else comp)
                           for (ib, comp) in zip(is_batched, proj_splats)]

    # Identify intersecting tiles
    tile_width = math.ceil(width / float(tile_size))
    tile_height = math.ceil(height / float(tile_size))
    tiles_per_gauss, isect_ids, flatten_ids = isect_tiles(
        means2d,
        radii,
        depths,
        tile_size,
        tile_width,
        tile_height,
        segmented=segmented,
        packed=packed,
        n_images=I,
        image_ids=image_ids,
        gaussian_ids=gaussian_ids,
    )
    # print("rank", world_rank, "Before isect_offset_encode")
    isect_offsets = isect_offset_encode(isect_ids, I, tile_width, tile_height)
    isect_offsets = isect_offsets.reshape(batch_dims + (C, tile_height, tile_width))

    meta.update(
        {
            "tile_width": tile_width,
            "tile_height": tile_height,
            "tiles_per_gauss": tiles_per_gauss,
            "isect_ids": isect_ids,
            "flatten_ids": flatten_ids,
            "isect_offsets": isect_offsets,
            "width": width,
            "height": height,
            "tile_size": tile_size,
            "n_batches": B,
            "n_cameras": C,
        }
    )

    # print("rank", world_rank, "Before rasterize_to_pixels")
    kwargs = {}
    if with_eval3d:
        kwargs = {
            "viewmats": viewmats,
            "Ks": Ks,
            "camera_model": camera_model,
            "radial_coeffs": radial_coeffs,
            "tangential_coeffs": tangential_coeffs,
            "thin_prism_coeffs": thin_prism_coeffs,
        }
    render_colors, render_alphas, render_meta = rasterize_to_pixels(
        primitive,
        with_eval3d,
        proj_splats,
        width,
        height,
        tile_size,
        isect_offsets,
        flatten_ids,
        backgrounds=backgrounds,
        packed=packed,
        absgrad=absgrad,
        **kwargs
    )
    meta.update(render_meta)

    # `rasterize_to_pixels` outputs log ray depth, map it to ray depth
    if len(render_colors) > 1:
        render_colors = (
            render_colors[0],
            torch.exp(render_colors[1] / render_alphas.clamp(min=1e-10)),
            *render_colors[2:]
        )

    # if "ED" in render_mode:
    #     # normalize the accumulated depth to get the expected depth
    #     depth_idx = 3 if "RGB" in render_mode else 0
    #     render_colors = torch.cat(
    #         [
    #             render_colors[..., :depth_idx],
    #             render_colors[..., depth_idx:depth_idx+1] / render_alphas.clamp(min=1e-10),
    #             render_colors[..., depth_idx+1:],
    #         ],
    #         dim=-1,
    #     )
    # if "N" in render_mode:
    #     render_colors = torch.cat(
    #         [
    #             render_colors[..., :-3],
    #             render_colors[..., -3:] / render_colors[..., -3:].norm(dim=-1, keepdim=True).clamp(min=1e-10),
    #         ],
    #         dim=-1,
    #     )

    return render_colors, render_alphas, meta

