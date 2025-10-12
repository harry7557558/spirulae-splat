import torch
import torch.nn.functional as F
from torch import Tensor
from typing import Dict, Optional, Tuple, Literal
import math


from gsplat.cuda._wrapper import (
    RollingShutterType,
    FThetaCameraDistortionParameters,
    FThetaPolynomialType,
    fully_fused_projection,
    fully_fused_projection_with_ut,
    isect_offset_encode,
    isect_tiles,
    rasterize_to_pixels,
    rasterize_to_pixels_eval3d,
    spherical_harmonics,
)
from spirulae_splat.splat.cuda._wrapper import (
    intersect_splat_tile,
    fully_fused_projection_hetero,
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
    means: Tensor,  # [..., N, 3]
    quats: Tensor,  # [..., N, 4]
    scales: Tensor,  # [..., N, 3]
    opacities: Tensor,  # [..., N]
    colors: Tensor,  # [..., (C,) N, D] or [..., (C,) N, K, 3]
    viewmats: Tensor,  # [..., C, 4, 4]
    Ks: Tensor,  # [..., C, 3, 3]
    width: int,
    height: int,
    near_plane: float = 0.01,
    far_plane: float = 1e10,
    radius_clip: float = 0.0,
    eps2d: float = 0.3,
    sh_degree: Optional[int] = None,
    packed: bool = True,
    use_bvh: bool = False,
    tile_size: int = 16,
    backgrounds: Optional[Tensor] = None,
    render_mode: Literal["RGB", "D", "ED", "RGB+D", "RGB+ED"] = "RGB",
    sparse_grad: bool = False,
    absgrad: bool = False,
    rasterize_mode: Literal["classic", "antialiased"] = "classic",
    channel_chunk: int = 32,
    distributed: bool = False,
    camera_model: Literal["pinhole", "ortho", "fisheye", "ftheta"] = "pinhole",
    segmented: bool = False,
    with_ut: bool = False,
    with_eval3d: bool = False,
    # distortion
    radial_coeffs: Optional[Tensor] = None,  # [..., C, 6] or [..., C, 4]
    tangential_coeffs: Optional[Tensor] = None,  # [..., C, 2]
    thin_prism_coeffs: Optional[Tensor] = None,  # [..., C, 4]
    ftheta_coeffs: Optional[FThetaCameraDistortionParameters] = None,
    # rolling shutter
    rolling_shutter: RollingShutterType = RollingShutterType.GLOBAL,
    viewmats_rs: Optional[Tensor] = None,  # [..., C, 4, 4]
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
        **Support N-D Features**: If `sh_degree` is None,
        the `colors` is expected to be with shape [..., N, D] or [..., C, N, D], in which D is the channel of
        the features to be rendered. The computation is slow when D > 32 at the moment.
        If `sh_degree` is set, the `colors` is expected to be the SH coefficients with
        shape [..., N, K, 3] or [..., C, N, K, 3], where K is the number of SH bases. In this case, it is expected
        that :math:`(\\textit{sh_degree} + 1) ^ 2 \\leq K`, where `sh_degree` controls the
        activated bases in the SH coefficients.

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

    .. note::
        **Speed-up for Large Scenes**: The `radius_clip` argument is extremely helpful for
        speeding up large scale scenes or scenes with large depth of fields. Gaussians with
        2D radius smaller or equal than this value (in pixel unit) will be skipped during rasterization.
        This will skip all the far-away Gaussians that are too small to be seen in the image.
        But be warned that if there are close-up Gaussians that are also below this threshold, they will
        also get skipped (which is rarely happened in practice). This is by default disabled by setting
        `radius_clip` to 0.0.

    .. note::
        **Antialiased Rendering**: If `rasterize_mode` is "antialiased", the function will
        apply a view-dependent compensation factor
        :math:`\\rho=\\sqrt{\\frac{Det(\\Sigma)}{Det(\\Sigma+ \\epsilon I)}}` to Gaussian
        opacities, where :math:`\\Sigma` is the projected 2D covariance matrix and :math:`\\epsilon`
        is the `eps2d`. This will make the rendered image more antialiased, as proposed in
        the paper `Mip-Splatting: Alias-free 3D Gaussian Splatting <https://arxiv.org/pdf/2311.16493>`_.

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
        radius_clip: Gaussians with 2D radius smaller or equal than this value will be
            skipped. This is extremely helpful for speeding up large scale scenes.
            Default is 0.0.
        eps2d: An epsilon added to the egienvalues of projected 2D covariance matrices.
            This will prevents the projected GS to be too small. For example eps2d=0.3
            leads to minimal 3 pixel unit. Default is 0.3.
        sh_degree: The SH degree to use, which can be smaller than the total
            number of bands. If set, the `colors` should be [..., (C,) N, K, 3] SH coefficients,
            else the `colors` should be [..., (C,) N, D] post-activation color values. Default is None.
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
        rasterize_mode: The rasterization mode. Supported modes are "classic" and
            "antialiased". Default is "classic".
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
        ftheta_coeffs: F-Theta camera distortion coefficients shared for all cameras.
            Default is None. See `FThetaCameraDistortionParameters` for details.
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

    batch_dims = means.shape[:-2]
    num_batch_dims = len(batch_dims)
    B = math.prod(batch_dims)
    N = means.shape[-2]
    C = viewmats.shape[-3]
    I = B * C
    device = means.device
    assert means.shape == batch_dims + (N, 3), means.shape
    assert quats.shape == batch_dims + (N, 4), quats.shape
    assert scales.shape == batch_dims + (N, 3), scales.shape
    assert opacities.shape == batch_dims + (N,), opacities.shape
    assert viewmats.shape == batch_dims + (C, 4, 4), viewmats.shape
    assert Ks.shape == batch_dims + (C, 3, 3), Ks.shape
    assert render_mode in ["RGB", "D", "ED", "RGB+D", "RGB+ED"], render_mode

    def reshape_view(C: int, world_view: torch.Tensor, N_world: list) -> torch.Tensor:
        view_list = list(
            map(
                lambda x: x.split(int(x.shape[0] / C), dim=0),
                world_view.split([C * N_i for N_i in N_world], dim=0),
            )
        )
        return torch.stack([torch.cat(l, dim=0) for l in zip(*view_list)], dim=0)

    if sh_degree is None:
        # treat colors as post-activation values, should be in shape [..., N, D] or [..., C, N, D]
        assert (
            colors.dim() == num_batch_dims + 2
            and colors.shape[:-1] == batch_dims + (N,)
        ) or (
            colors.dim() == num_batch_dims + 3
            and colors.shape[:-1] == batch_dims + (C, N)
        ), colors.shape
        if distributed:
            assert (
                colors.dim() == num_batch_dims + 2
            ), "Distributed mode only supports per-Gaussian colors."
    else:
        # treat colors as SH coefficients, should be in shape [..., N, K, 3] or [..., C, N, K, 3]
        # Allowing for activating partial SH bands
        assert (
            colors.dim() == num_batch_dims + 3
            and colors.shape[:-2] == batch_dims + (N,)
            and colors.shape[-1] == 3
        ) or (
            colors.dim() == num_batch_dims + 4
            and colors.shape[:-2] == batch_dims + (C, N)
            and colors.shape[-1] == 3
        ), colors.shape
        assert (sh_degree + 1) ** 2 <= colors.shape[-2], colors.shape
        if distributed:
            assert (
                colors.dim() == num_batch_dims + 3
            ), "Distributed mode only supports per-Gaussian colors."

    if absgrad:
        assert not distributed, "AbsGrad is not supported in distributed mode."

    if (
        radial_coeffs is not None
        or tangential_coeffs is not None
        or thin_prism_coeffs is not None
        or ftheta_coeffs is not None
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

    # Implement the multi-GPU strategy proposed in
    # `On Scaling Up 3D Gaussian Splatting Training <https://arxiv.org/abs/2406.18533>`.
    #
    # If in distributed mode, we distribute the projection computation over Gaussians
    # and the rasterize computation over cameras. So first we gather the cameras
    # from all ranks for projection.
    if distributed:
        assert batch_dims == (), "Distributed mode does not support batch dimensions"
        world_rank = torch.distributed.get_rank()
        world_size = torch.distributed.get_world_size()

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
                scales.contiguous(),
                opacities.contiguous(),
                quats.contiguous(),
                width,
                height,
                viewmats.contiguous(),
                Ks.contiguous(),
            )
        torch.cuda.synchronize()
        time1 = perf_counter()
        print(1e3*(time1-time0), 'ms')
        # if 1e3*(time1-time0) > 100:
        #     dump_all()
        #     exit(0)
        proj_results = fully_fused_projection_hetero(
            means,
            quats,
            scales,
            viewmats,
            Ks,
            width,
            height,
            intersection_count_map,
            intersection_splat_id,
            eps2d=eps2d,
            near_plane=near_plane,
            far_plane=far_plane,
            radius_clip=radius_clip,
            sparse_grad=sparse_grad,
            calc_compensations=(rasterize_mode == "antialiased"),
            camera_model=camera_model,
            opacities=opacities,  # use opacities to compute a tigher bound for radii.
        )

    elif with_ut:
        proj_results = fully_fused_projection_with_ut(
            means,
            quats,
            scales,
            opacities,  # use opacities to compute a tigher bound for radii.
            viewmats,
            Ks,
            width,
            height,
            eps2d=eps2d,
            near_plane=near_plane,
            far_plane=far_plane,
            radius_clip=radius_clip,
            calc_compensations=(rasterize_mode == "antialiased"),
            camera_model=camera_model,
            radial_coeffs=radial_coeffs,
            tangential_coeffs=tangential_coeffs,
            thin_prism_coeffs=thin_prism_coeffs,
            ftheta_coeffs=ftheta_coeffs,
            rolling_shutter=rolling_shutter,
            viewmats_rs=viewmats_rs,
        )

    else:
        # Project Gaussians to 2D. Directly pass in {quats, scales} is faster than precomputing covars.
        proj_results = fully_fused_projection(
            means,
            None,  # covars
            quats,
            scales,
            viewmats,
            Ks,
            width,
            height,
            eps2d=eps2d,
            packed=packed,
            near_plane=near_plane,
            far_plane=far_plane,
            radius_clip=radius_clip,
            sparse_grad=sparse_grad,
            calc_compensations=(rasterize_mode == "antialiased"),
            camera_model=camera_model,
            opacities=opacities,  # use opacities to compute a tigher bound for radii.
        )

    if use_bvh:
        (
            camera_ids,
            gaussian_ids,
            radii,
            means2d,
            depths,
            conics,
            compensations,
        ) = proj_results
        batch_ids, image_ids = 0, camera_ids
        opacities = opacities.view(B, N)[batch_ids, gaussian_ids]  # [nnz]
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
            compensations,
        ) = proj_results
        opacities = opacities.view(B, N)[batch_ids, gaussian_ids]  # [nnz]
        image_ids = batch_ids * C + camera_ids
    else:
        # The results are with shape [..., C, N, ...]. Only the elements with radii > 0 are valid.
        radii, means2d, depths, conics, compensations = proj_results
        opacities = torch.broadcast_to(
            opacities[..., None, :], batch_dims + (C, N)
        )  # [..., C, N]
        batch_ids, camera_ids, gaussian_ids = None, None, None
        image_ids = None

    if compensations is not None:
        opacities = opacities * compensations

    meta.update(
        {
            # global batch and camera ids
            "batch_ids": batch_ids,
            "camera_ids": camera_ids,
            # local gaussian_ids
            "gaussian_ids": gaussian_ids,
            "radii": radii,
            "means2d": means2d,
            "depths": depths,
            "conics": conics,
            "opacities": opacities,
        }
    )
    # if heterogeneous:
    #     meta.update({
    #         "intersection_count": intersection_count_map[1:]-intersection_count_map[:-1]
    #     })

    # Turn colors into [..., C, N, D] or [..., nnz, D] to pass into rasterize_to_pixels()
    if sh_degree is None:
        # Colors are post-activation values, with shape [..., N, D] or [..., C, N, D]
        if packed:
            if colors.dim() == num_batch_dims + 2:
                # Turn [..., N, D] into [nnz, D]
                colors = colors.view(B, N, -1)[batch_ids, gaussian_ids]
            else:
                # Turn [..., C, N, D] into [nnz, D]
                colors = colors.view(B, C, N, -1)[batch_ids, camera_ids, gaussian_ids]
        else:
            if colors.dim() == num_batch_dims + 2:
                # Turn [..., N, D] into [..., C, N, D]
                colors = torch.broadcast_to(
                    colors[..., None, :, :], batch_dims + (C, N, -1)
                )
            else:
                # colors is already [..., C, N, D]
                pass
    else:
        # Colors are SH coefficients, with shape [..., N, K, 3] or [..., C, N, K, 3]
        campos = torch.inverse(viewmats)[..., :3, 3]  # [..., C, 3]
        if viewmats_rs is not None:
            campos_rs = torch.inverse(viewmats_rs)[..., :3, 3]
            campos = 0.5 * (campos + campos_rs)  # [..., C, 3]
        if packed:
            dirs = (
                means.view(B, N, 3)[batch_ids, gaussian_ids]
                - campos.view(B, C, 3)[batch_ids, camera_ids]
            )  # [nnz, 3]
            masks = (radii > 0).all(dim=-1)  # [nnz]
            if colors.dim() == num_batch_dims + 3:
                # Turn [..., N, K, 3] into [nnz, 3]
                shs = colors.view(B, N, -1, 3)[batch_ids, gaussian_ids]  # [nnz, K, 3]
            else:
                # Turn [..., C, N, K, 3] into [nnz, 3]
                shs = colors.view(B, C, N, -1, 3)[
                    batch_ids, camera_ids, gaussian_ids
                ]  # [nnz, K, 3]
            colors = spherical_harmonics(sh_degree, dirs, shs, masks=masks)  # [nnz, 3]
        else:
            dirs = means[..., None, :, :] - campos[..., None, :]  # [..., C, N, 3]
            masks = (radii > 0).all(dim=-1)  # [..., C, N]
            if colors.dim() == num_batch_dims + 3:
                # Turn [..., N, K, 3] into [..., C, N, K, 3]
                shs = torch.broadcast_to(
                    colors[..., None, :, :, :], batch_dims + (C, N, -1, 3)
                )
            else:
                # colors is already [..., C, N, K, 3]
                shs = colors
            colors = spherical_harmonics(
                sh_degree, dirs, shs, masks=masks
            )  # [..., C, N, 3]
        # make it apple-to-apple with Inria's CUDA Backend.
        colors = torch.clamp_min(colors + 0.5, 0.0)

    # If in distributed mode, we need to scatter the GSs to the destination ranks, based
    # on which cameras they are visible to, which we already figured out in the projection
    # stage.
    if distributed:
        if packed:
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
            (means2d, depths, conics, opacities, colors) = all_to_all_tensor_list(
                world_size,
                [means2d, depths, conics, opacities, colors],
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

            (means2d, depths, conics, opacities, colors) = all_to_all_tensor_list(
                world_size,
                [
                    means2d.flatten(0, 1),
                    depths.flatten(0, 1),
                    conics.flatten(0, 1),
                    opacities.flatten(0, 1),
                    colors.flatten(0, 1),
                ],
                splits=[C_i * N for C_i in C_world],
                output_splits=[C * N_i for N_i in N_world],
            )
            means2d = reshape_view(C, means2d, N_world)
            depths = reshape_view(C, depths, N_world)
            conics = reshape_view(C, conics, N_world)
            opacities = reshape_view(C, opacities, N_world)
            colors = reshape_view(C, colors, N_world)

    # Rasterize to pixels
    if render_mode in ["RGB+D", "RGB+ED"]:
        colors = torch.cat((colors, depths[..., None]), dim=-1)
        if backgrounds is not None:
            backgrounds = torch.cat(
                [
                    backgrounds,
                    torch.zeros(batch_dims + (C, 1), device=backgrounds.device),
                ],
                dim=-1,
            )
    elif render_mode in ["D", "ED"]:
        colors = depths[..., None]
        if backgrounds is not None:
            backgrounds = torch.zeros(batch_dims + (C, 1), device=backgrounds.device)
    else:  # RGB
        pass

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
    if colors.shape[-1] > channel_chunk:
        # slice into chunks
        n_chunks = (colors.shape[-1] + channel_chunk - 1) // channel_chunk
        render_colors, render_alphas = [], []
        for i in range(n_chunks):
            colors_chunk = colors[..., i * channel_chunk : (i + 1) * channel_chunk]
            backgrounds_chunk = (
                backgrounds[..., i * channel_chunk : (i + 1) * channel_chunk]
                if backgrounds is not None
                else None
            )
            if with_eval3d:
                render_colors_, render_alphas_ = rasterize_to_pixels_eval3d(
                    means,
                    quats,
                    scales,
                    colors_chunk,
                    opacities,
                    viewmats,
                    Ks,
                    width,
                    height,
                    tile_size,
                    isect_offsets,
                    flatten_ids,
                    backgrounds=backgrounds_chunk,
                    camera_model=camera_model,
                    radial_coeffs=radial_coeffs,
                    tangential_coeffs=tangential_coeffs,
                    thin_prism_coeffs=thin_prism_coeffs,
                    ftheta_coeffs=ftheta_coeffs,
                    rolling_shutter=rolling_shutter,
                    viewmats_rs=viewmats_rs,
                )
            else:
                # TODO
                render_colors_, render_alphas_ = rasterize_to_pixels(
                    means2d,
                    conics,
                    colors_chunk,
                    opacities,
                    width,
                    height,
                    tile_size,
                    isect_offsets,
                    flatten_ids,
                    backgrounds=backgrounds_chunk,
                    packed=packed,
                    absgrad=absgrad,
                )
            render_colors.append(render_colors_)
            render_alphas.append(render_alphas_)
        render_colors = torch.cat(render_colors, dim=-1)
        render_alphas = render_alphas[0]  # discard the rest
    else:
        if with_eval3d:
            render_colors, render_alphas = rasterize_to_pixels_eval3d(
                means,
                quats,
                scales,
                colors,
                opacities,
                viewmats,
                Ks,
                width,
                height,
                tile_size,
                isect_offsets,
                flatten_ids,
                backgrounds=backgrounds,
                camera_model=camera_model,
                radial_coeffs=radial_coeffs,
                tangential_coeffs=tangential_coeffs,
                thin_prism_coeffs=thin_prism_coeffs,
                ftheta_coeffs=ftheta_coeffs,
                rolling_shutter=rolling_shutter,
                viewmats_rs=viewmats_rs,
            )
        else:
            render_colors, render_alphas = rasterize_to_pixels(
                means2d,
                conics,
                colors,
                opacities,
                width,
                height,
                tile_size,
                isect_offsets,
                flatten_ids,
                backgrounds=backgrounds,
                packed=packed,
                absgrad=absgrad,
            )
    if render_mode in ["ED", "RGB+ED"]:
        # normalize the accumulated depth to get the expected depth
        render_colors = torch.cat(
            [
                render_colors[..., :-1],
                render_colors[..., -1:] / render_alphas.clamp(min=1e-10),
            ],
            dim=-1,
        )

    return render_colors, render_alphas, meta

