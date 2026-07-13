"""Viewer-side annotation: thin wrappers around the engine-side BVH cache +
blit kernel.

Two entry points:

- ``ensure_viewer_initialized(cameras)`` -- idempotent; the first call uploads
  the dataset-wide camera arrays (intrins / dist_coeffs / c2w / model / W / H)
  to the engine's device-side viewer cache + allocates the per-camera
  thumbnail buffer + done-mask. Subsequent calls no-op.

- ``annotate_train_cameras(buffer_key, buffer, depth, alpha, view_camera, ...)``
  -- calls ``_C.engine_blit_view`` which (1) builds the BVH from the engine's
  cached camera arrays on the first ``show_training_cameras=True`` call and
  reuses it forever, (2) reads thumbnails out of the engine's cache (no
  thumbnails tensor argument), (3) writes a freshly-allocated CUDA uint8
  [H, W, 3] output.

The thumbnails themselves get populated automatically as training data
flows through ``set_training_data{,_warped}`` -- the engine's training
step hook downsamples each batch's POST-split RGB into the cache.
"""

from __future__ import annotations

import numpy as np
import torch

from spirulae_splat.modules.camera import Cameras, CameraType
from spirulae_splat.splat.cuda import _C


def _knn_dist(x: torch.Tensor, k: int = 4):
    if len(x) == 1:
        return 1.0

    # Deduplicate near-coincident positions before running knn. The
    # post-split camera arrays produced by warp_to_pinhole (K=5) and
    # equirectangular split (K=6) emit K sub-cameras at the *same*
    # translation per input image -- only the orientations differ. With
    # k=4 nearest neighbors and cluster size >= 5, every "nearest"
    # distance is 0, so the median collapses to 0 and the frustum scale
    # becomes 0 (invisible cameras, including pinhole entries in a mixed
    # dataset, since the whole visualization shares one camera_size).
    x_np = x.detach().cpu().numpy().astype(np.float32)
    scale = float(np.max(np.abs(x_np))) + 1e-9
    quantized = np.round(x_np / scale * 1e6).astype(np.int64)
    _, uniq_idx = np.unique(quantized, axis=0, return_index=True)
    x_unique = x_np[np.sort(uniq_idx)]
    n = len(x_unique)
    if n == 1:
        return 1.0
    if n <= k:
        from scipy.spatial.distance import pdist
        d = pdist(x_unique)
        d = d[d > 0.0]
        return float(np.median(d)) if d.size else 1.0
    from scipy.spatial import cKDTree
    tree = cKDTree(x_unique, balanced_tree=False)
    distances, _ = tree.query(x_unique, k=k + 1, workers=-1)
    return float(np.median(np.amax(distances[:, 1:], axis=-1)))


# Must match Common.cuh's CameraModelType.
_CAMERA_MODEL_INT = {
    "PINHOLE": 0,
    "FISHEYE": 1,
    "EQUISOLID": 2,
    "EQUIRECTANGULAR": 3,
}


_INIT_KEY = "_viewer_engine_initialized"


def _tv(t, keepalive):
    """TorchTensorView helper: (data_ptr, element_size, shape).

    `.cuda()` / `.contiguous()` may create a new temporary tensor; the caller
    MUST hold a reference to it through the kernel launch or the CUDA
    allocator will recycle the memory before the kernel reads it (the
    classic "rgb shows alpha" symptom). The `keepalive` list does that.
    """
    if t is None:
        return (0, 4, [0])
    if not t.is_cuda:
        t = t.cuda()
    t = t.contiguous()
    keepalive.append(t)
    return (t.data_ptr(), t.element_size(), list(t.shape))


def ensure_viewer_initialized(cameras: Cameras) -> None:
    """Idempotent: uploads POST-split camera arrays + allocates thumbnail cache
    on first call. cameras must already be the post-split layout (PINHOLE per
    cam for warp_to_pinhole; full set for EQUIRECTANGULAR split)."""
    if getattr(cameras, _INIT_KEY, False):
        return

    cuda_cameras = cameras.to("cuda")
    N = len(cuda_cameras)

    # Camera-to-world in the y/z-flipped form the blit kernel expects.
    R = cuda_cameras.camera_to_worlds[:, :3, :3]
    T = cuda_cameras.camera_to_worlds[:, :3, 3:4]
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]], device=R.device, dtype=R.dtype)
    c2w = torch.cat((R, T), dim=-1).contiguous()       # [N, 3, 4]

    intrins = cuda_cameras.intrins.contiguous()        # [N, 4]
    dist    = cuda_cameras.distortion_params.contiguous()  # [N, 10]
    widths  = cuda_cameras.width.to(torch.int32).contiguous().reshape(-1)
    heights = cuda_cameras.height.to(torch.int32).contiguous().reshape(-1)
    cmodels = torch.tensor(
        [_CAMERA_MODEL_INT[c.upper()] for c in cuda_cameras.camera_type],
        dtype=torch.int32, device="cuda",
    ).contiguous()

    camera_size = 0.2 * _knn_dist(T.squeeze(-1))

    # Pin to the cameras object so the tensors aren't GC'd while engine refs them.
    cameras._viewer_engine_keepalive = (intrins, dist, c2w, widths, heights, cmodels)

    keep = []
    _C.engine_viewer_init(
        _tv(cmodels, keep), _tv(intrins, keep), _tv(dist, keep), _tv(c2w, keep),
        _tv(widths, keep), _tv(heights, keep), float(camera_size),
    )
    # Axes/grid overlay (drawn when a render asks for show_grid). The grid is
    # axis-aligned in the engine's own frame -- the frame splats are trained
    # and saved in -- so its lines mark round coordinates of the exported
    # model; the cell size then adapts to the per-render grid_dist.
    radius = float(T.squeeze(-1).norm(dim=-1).max().clamp(min=1e-6))
    _C.engine_viewer_set_grid(radius, 0.0)
    setattr(cameras, _INIT_KEY, True)


@torch.no_grad()
def annotate_train_cameras(
    buffer_key: str,
    buffer: torch.Tensor,           # [H,W,1|3] float CUDA
    depth: torch.Tensor,            # [H,W,1] float CUDA (may be None)
    alpha: torch.Tensor,            # [H,W,1] float CUDA
    view_camera: Cameras,
    cameras: Cameras,
    **kwargs,
) -> torch.Tensor:
    """Annotated uint8 [H,W,3] CUDA tensor."""
    ensure_viewer_initialized(cameras)

    if buffer.ndim == 4:
        buffer = buffer.squeeze(0)
    if depth is not None and depth.ndim == 4:
        depth = depth.squeeze(0)
    if alpha.ndim == 4:
        alpha = alpha.squeeze(0)

    # Build view camera viewmat in the engine's frame.
    R = view_camera.camera_to_worlds[:, :3, :3]
    T = view_camera.camera_to_worlds[:, :3, 3:4]
    R = R * torch.tensor([[[1.0, -1.0, -1.0]]], device=R.device, dtype=R.dtype)
    R_inv = R.transpose(-1, -2)
    T_inv = -torch.bmm(R_inv, T)
    view_viewmat = torch.eye(4, dtype=view_camera.camera_to_worlds.dtype,
                             device=R.device)[None].repeat(len(view_camera), 1, 1)
    view_viewmat[:, :3, :3] = R_inv
    view_viewmat[:, :3, 3:4] = T_inv

    show_train = bool(kwargs.get("show_training_cameras", False))
    show_grid = bool(kwargs.get("show_grid", False))
    # Nav view distance + orbit target in the ENGINE frame (the trainer
    # already applied its similarity); dist 0 keeps the grid's current
    # cell size and center.
    grid_dist = float(kwargs.get("grid_dist", 0.0))
    grid_target = tuple(float(x) for x in
                        kwargs.get("grid_target", (0.0, 0.0, 0.0)))

    h, w = buffer.shape[0], buffer.shape[1]
    out_rgb = torch.empty(h, w, 3, dtype=torch.uint8,
                          device=buffer.device if buffer.is_cuda else torch.device("cuda"))

    view_intrins = view_camera.intrins.cuda().contiguous()
    view_dist    = view_camera.distortion_params.cuda().contiguous()
    view_viewmat = view_viewmat.cuda().contiguous()

    keep = []
    # depth may be missing for alpha-only render; pass an empty TV.
    depth_tv = _tv(depth, keep) if depth is not None else (0, 4, [0])

    _C.engine_blit_view(
        buffer_key,
        _tv(buffer, keep), depth_tv, _tv(alpha, keep),
        _CAMERA_MODEL_INT[view_camera.camera_type[0].upper()],
        _tv(view_intrins, keep), _tv(view_viewmat, keep), _tv(view_dist, keep),
        show_train,
        show_grid,
        grid_dist,
        grid_target[0], grid_target[1], grid_target[2],
        _tv(out_rgb, keep),
    )
    return out_rgb
