import pytest
import torch
from torch.func import vjp  # type: ignore

from spirulae_splat.splat._torch_impl import quat_to_rotmat

from spirulae_splat.modules.core import Renderer

# commit 48df61993fed8e4742d0758b6c8c2b4e599c4124
from gsplat.rendering import rasterization as gsplat_rasterization

from spirulae_splat.splat.cuda import (
    _C,
    ray_depth_to_linear_depth,
)

from utils import check_close, timeit

from typing import Literal

device = torch.device("cuda:0")


B, W, H = 4, 1440, 1080
N, SH_DEGREE, SH_DEGREE_TO_USE = 200000, 3, 3
PACKED = True
IS_FISHEYE = False
IS_ANTIALIASED = False
WITH_UT = False

def rasterize_ssplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    renderer = Renderer(
        primitive="3dgut" if WITH_UT else ["3dgs", "mip"][IS_ANTIALIASED],
        splats_world=(means, quats, scales, opacities.unsqueeze(-1), features_dc, features_sh),
        cur_num_splats=len(means),
        packed=PACKED,
    )
    camera_model = ["pinhole", "fisheye"][IS_FISHEYE]
    # quats = torch.nn.functional.normalize(quats, dim=-1)  # affects gradient
    intrins = torch.stack([Ks[..., 0, 0], Ks[..., 1, 1], Ks[..., 0, 2], Ks[..., 1, 2]], dim=-1)  # fx, fy, cx, cy
    renderer.set_params(
        viewmats=viewmats,  # [C, 4, 4]
        intrins=intrins,  # [C, 4]
        width=W,
        height=H,
        sh_degree_to_use=SH_DEGREE_TO_USE,
        camera_model=camera_model,
    )
    renderer.forward()
    rgbd = renderer.render_colors
    Ts = renderer.render_Ts
    rgbd = [*rgbd[:2]]
    # if WITH_UT:
    #     try:
    #         rgbd[1] = ray_depth_to_linear_depth(rgbd[1], camera_model, intrins)  # TODO: f(E[X]) != E[f(X)]
    #     except AttributeError:
    #         pass  # ray_depth_to_linear_depth_forward not compiled, depth comparison will differ
    return (*rgbd, 1.0 - Ts), renderer

def rasterize_gsplat(means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks):
    rgbd, alpha, meta = gsplat_rasterization(
        means=means,
        quats=quats,
        scales=torch.exp(scales),
        opacities=torch.sigmoid(opacities).squeeze(-1),
        colors=torch.concatenate([features_dc.unsqueeze(1), features_sh], dim=1),  # TODO: slow
        viewmats=viewmats,  # [C, 4, 4]
        Ks=Ks,  # [C, 3, 3]
        width=W,
        height=H,
        sh_degree=SH_DEGREE_TO_USE,
        packed=PACKED and not WITH_UT,
        sparse_grad=False,
        rasterize_mode=["classic", "antialiased"][IS_ANTIALIASED],
        eps2d=(0.1 if IS_ANTIALIASED else 0.3),
        distributed=False,
        camera_model=["pinhole", "fisheye"][IS_FISHEYE],
        with_ut=WITH_UT,
        with_eval3d=WITH_UT,
        render_mode="RGB+ED",
    )
    rgb, depth = rgbd[..., :3], rgbd[..., 3:]
    # depth = torch.where(alpha == 0.0, depth.amax(), depth)
    return rgb, depth, alpha


def get_inputs(seed=44):
    torch.manual_seed(seed)

    means = torch.randn((N, 3), device=device)
    quats = torch.randn((N, 4), device=device)
    scales = torch.randn((N, 3), device=device)*0.7 - 5.0
    # means[0] -= 1e4
    opacities = torch.randn((N,), device=device)
    features_dc = torch.rand((N, 3), device=device)
    features_sh = 0.2 * torch.randn((N, (SH_DEGREE+1)**2-1, 3), device=device)

    quats = torch.nn.functional.normalize(quats)

    viewmats = torch.eye(4)[None].repeat(B, 1, 1)
    viewmats[:, :3, 3] = 1.0 * torch.randn((B, 3))
    viewmats[:, 2, 3] = torch.abs(viewmats[:, 2, 3])
    rotation = torch.randn(B, 4)
    rotation[:, 1:] *= 0.1
    rotation /= torch.norm(rotation, dim=-1, keepdim=True)
    viewmats[:, :3, :3] = quat_to_rotmat(rotation)
    viewmats = viewmats.contiguous().to(device)

    cx, cy = 0.5*W, 0.5*H
    fx, fy = 0.4*W, 0.4*W  # 0.4
    Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]]).repeat(B, 1, 1).to(device)
    Ks *= torch.exp(0.2*torch.randn_like(Ks))

    inputs = (means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks)
    inputs = [torch.nn.Parameter(x.contiguous().clone()) for x in inputs]
    if WITH_UT:
        inputs = inputs[:-2] + [viewmats, Ks]
    return inputs


# ---------------------------------------------------------------------------
# Representative synthetic scenes for profiling across gaussian distributions.
# These are parametric stand-ins that reproduce the *projected* spatial density
# and depth skew of real captures (not photometric realism). Each returns the
# same 8-tuple as get_inputs().
# ---------------------------------------------------------------------------
import math


def _look_at(eye, target, up=(0.0, 1.0, 0.0)):
    eye = torch.tensor(eye, dtype=torch.float32)
    target = torch.tensor(target, dtype=torch.float32)
    up = torch.tensor(up, dtype=torch.float32)
    fwd = target - eye
    fwd = fwd / fwd.norm().clamp_min(1e-8)          # camera looks along +z
    right = torch.linalg.cross(up, fwd)
    right = right / right.norm().clamp_min(1e-8)
    tup = torch.linalg.cross(fwd, right)
    R = torch.stack([right, tup, fwd], dim=0)       # rows = camera axes in world
    vm = torch.eye(4)
    vm[:3, :3] = R
    vm[:3, 3] = -(R @ eye)
    return vm


def _plane(n, center, u, v, jitter, gen):
    # n points on a planar slab: center + a*u + b*v + c*normal, a,b in [-1,1]
    center = torch.tensor(center, dtype=torch.float32)
    u = torch.tensor(u, dtype=torch.float32)
    v = torch.tensor(v, dtype=torch.float32)
    nrm = torch.linalg.cross(u, v)
    nrm = nrm / nrm.norm().clamp_min(1e-8)
    a = torch.rand(n, 1, generator=gen) * 2 - 1
    b = torch.rand(n, 1, generator=gen) * 2 - 1
    c = (torch.rand(n, 1, generator=gen) * 2 - 1) * jitter
    return center + a * u + b * v + c * nrm


def _finish_scene(means, viewmats, scale_log_mean, opac_bias, seed):
    n = means.shape[0]
    g = torch.Generator().manual_seed(seed + 7)
    quats = torch.nn.functional.normalize(torch.randn(n, 4, generator=g))
    scales = torch.randn(n, 3, generator=g) * 0.5 + scale_log_mean
    opacities = torch.randn(n, generator=g) + opac_bias
    features_dc = torch.rand(n, 3, generator=g)
    features_sh = 0.2 * torch.randn(n, (SH_DEGREE + 1) ** 2 - 1, 3, generator=g)
    cx, cy = 0.5 * W, 0.5 * H
    fx = fy = 0.4 * W
    Ks = torch.tensor([[[fx, 0, cx], [0, fy, cy], [0, 0, 1]]]).repeat(viewmats.shape[0], 1, 1).float()
    inputs = (means, quats, scales, opacities, features_dc, features_sh, viewmats, Ks)
    inputs = [torch.nn.Parameter(x.contiguous().to(device).clone()) for x in inputs]
    if WITH_UT:
        inputs = inputs[:-2] + [viewmats.to(device), Ks.to(device)]
    return inputs


def gen_scene(kind, n, seed=44):
    g = torch.Generator().manual_seed(seed)

    if kind == "uniform":
        # baseline: gaussians spread through space, cameras around origin
        means = torch.randn(n, 3, generator=g)
        vms = [_look_at([3.0 * math.cos(2 * math.pi * k / B),
                         0.5 * (k - (B - 1) / 2),
                         3.0 * math.sin(2 * math.pi * k / B)], [0, 0, 0]) for k in range(B)]
        return _finish_scene(means, torch.stack(vms), -4.5, 0.0, seed)

    if kind == "object":
        # small object-centered capture: dense compact blob + sparse far shell
        n_obj = int(0.9 * n); n_bg = n - n_obj
        obj = 0.5 * torch.randn(n_obj, 3, generator=g)
        dirn = torch.nn.functional.normalize(torch.randn(n_bg, 3, generator=g))
        bg = dirn * (8 + 4 * torch.rand(n_bg, 1, generator=g))
        means = torch.cat([obj, bg], 0)
        vms = [_look_at([3.2 * math.cos(2 * math.pi * k / B),
                         0.7 * (k - (B - 1) / 2),
                         3.2 * math.sin(2 * math.pi * k / B)], [0, 0, 0]) for k in range(B)]
        return _finish_scene(means, torch.stack(vms), -5.3, 0.8, seed)

    if kind == "street":
        # narrow street: ground + two facades receding down +z, camera low.
        # grazing ground + converging facades -> depth grows toward the horizon.
        n3 = n // 3
        ground = _plane(n3, [0, 0, 20], [6, 0, 0], [0, 0, 20], 0.05, g); ground[:, 1] = 0.0
        facL = _plane(n3, [-3, 3, 20], [0, 3, 0], [0, 0, 20], 0.05, g)
        facR = _plane(n - 2 * n3, [3, 3, 20], [0, 3, 0], [0, 0, 20], 0.05, g)
        means = torch.cat([ground, facL, facR], 0)
        vms = [_look_at([0.4 * (k - (B - 1) / 2), 1.2, -1.0 - 0.3 * k],
                        [0.4 * (k - (B - 1) / 2), 1.0, 20]) for k in range(B)]
        return _finish_scene(means, torch.stack(vms), -4.8, 1.0, seed)

    if kind == "indoor":
        # multi-room: box of 6 planes, camera inside looking toward a corner so
        # two walls recede at a grazing angle (deep) and two are frontal (shallow).
        nf = n // 6
        planes = [
            _plane(nf, [0, -2, 0], [5, 0, 0], [0, 0, 5], 0.05, g),   # floor
            _plane(nf, [0, 3, 0], [5, 0, 0], [0, 0, 5], 0.05, g),    # ceiling
            _plane(nf, [-5, 0.5, 0], [0, 2.5, 0], [0, 0, 5], 0.05, g),  # wall x-
            _plane(nf, [5, 0.5, 0], [0, 2.5, 0], [0, 0, 5], 0.05, g),   # wall x+
            _plane(nf, [0, 0.5, -5], [5, 0, 0], [0, 2.5, 0], 0.05, g),  # wall z-
            _plane(n - 5 * nf, [0, 0.5, 5], [5, 0, 0], [0, 2.5, 0], 0.05, g),  # wall z+
        ]
        means = torch.cat(planes, 0)
        vms = [_look_at([-2 + 0.5 * k, 0.0, -2 + 0.3 * k], [4, 0, 4]) for k in range(B)]
        return _finish_scene(means, torch.stack(vms), -4.5, 0.5, seed)

    if kind == "garden":
        # unbounded outdoor: grazing ground + central foreground bushes + far shell
        n_g = n // 3; n_f = n // 3; n_b = n - n_g - n_f
        ground = _plane(n_g, [0, 0, 8], [12, 0, 0], [0, 0, 12], 0.03, g); ground[:, 1] = 0.0
        fg = torch.stack([2.0 * torch.randn(n_f, generator=g),
                          0.6 * torch.rand(n_f, generator=g),
                          8 + 2.0 * torch.randn(n_f, generator=g)], dim=1)
        dirn = torch.nn.functional.normalize(torch.randn(n_b, 3, generator=g))
        shell = dirn.abs() * torch.tensor([30.0, 15.0, 30.0]) + torch.tensor([0.0, 2.0, 8.0])
        means = torch.cat([ground, fg, shell], 0)
        vms = [_look_at([0.6 * (k - (B - 1) / 2), 1.4, -1.0], [0, 0.6, 8]) for k in range(B)]
        return _finish_scene(means, torch.stack(vms), -4.7, 0.7, seed)

    raise ValueError(f"unknown scene {kind}")


@torch.no_grad()
def scene_tile_stats(means, viewmats, Ks):
    # characterize skew: gaussian-center counts per tile (footprint ignored, so
    # absolute counts under-read true binning, but the *skew* is representative)
    m = means.detach().to(device).float()
    out = {}
    for name, (tx, ty) in [("macro", (64, 32)), ("micro", (8, 8))]:
        nx = (W + tx - 1) // tx; ny = (H + ty - 1) // ty
        mx = 0; nonempty = []
        for c in range(viewmats.shape[0]):
            R = viewmats[c, :3, :3].to(device).float(); t = viewmats[c, :3, 3].to(device).float()
            cam = m @ R.T + t
            z = cam[:, 2]
            u = Ks[c, 0, 0].item() * cam[:, 0] / z + Ks[c, 0, 2].item()
            v = Ks[c, 1, 1].item() * cam[:, 1] / z + Ks[c, 1, 2].item()
            ok = (z > 1e-4) & (u >= 0) & (u < W) & (v >= 0) & (v < H)
            idx = (v[ok] / ty).long() * nx + (u[ok] / tx).long()
            cnt = torch.bincount(idx, minlength=nx * ny)
            mx = max(mx, int(cnt.max().item()) if cnt.numel() else 0)
            nz = cnt[cnt > 0]
            if nz.numel():
                nonempty.append(nz.float())
        meanc = torch.cat(nonempty).mean().item() if nonempty else 0.0
        out[name + "_max"] = mx
        out[name + "_mean"] = meanc
    return out


def _time_ms(fn, warmup=3, repeat=30):
    from time import perf_counter
    for _ in range(warmup):
        fn()
    torch.cuda.synchronize()
    t0 = perf_counter()
    for _ in range(repeat):
        fn()
    torch.cuda.synchronize()
    return 1e3 * (perf_counter() - t0) / repeat


@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_scenes(kinds=None, n=200000):
    global N
    N = n
    if kinds is None:
        kinds = ["uniform", "object", "garden", "street", "indoor"]
    prim = "3dgut" if WITH_UT else ["3dgs", "mip"][IS_ANTIALIASED]

    def _isect_vram_mb():
        bd = _C.engine_get_pool_breakdown()
        return sum(e[1] for e in bd if e[0].lower().startswith("isect.")) / 1e6

    print(f"\n##### scene profiling (N={n}, primitive={prim}) #####")
    print(f"{'scene':8s} | {'macro max/mean':>16s} | {'micro max/mean':>15s} | {'fwd ms':>7s} | {'bwd ms':>7s} | {'isectMB':>7s}")
    for kind in kinds:
        try:
            inputs = gen_scene(kind, n)
            st = scene_tile_stats(inputs[0], inputs[6], inputs[7])
            cpu_inputs = [x.detach().cpu().contiguous() for x in inputs]
            _, renderer = rasterize_ssplat(*cpu_inputs)
            ft = _time_ms(lambda: _C.engine_forward_3dgs(
                renderer.primitive, renderer.sh_degree_to_use, renderer.packed))
            C, Hh, Ww = renderer.viewmats.shape[0], renderer.height, renderer.width
            v_rgb = torch.randn(C, Hh, Ww, 3, device=device)
            v_depth = torch.randn(C, Hh, Ww, 1, device=device)
            v_Ts = torch.randn(C, Hh, Ww, 1, device=device)
            bt = _time_ms(lambda: _C.engine_backward_from_render_grad(
                renderer._tv(v_rgb), renderer._tv(v_depth), renderer._tv(v_Ts)))
            vram = _isect_vram_mb()
            print(f"{kind:8s} | {st['macro_max']:7d}/{st['macro_mean']:7.0f} | "
                  f"{st['micro_max']:6d}/{st['micro_mean']:5.1f} | {ft:7.2f} | {bt:7.2f} | {vram:7.1f}")
            _C.engine_reset()
        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                _C.engine_reset()
            except Exception:
                pass



@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def test_rasterization():

    inputs = get_inputs()
    _inputs = get_inputs()

    inputs = [x.cpu() for x in inputs]
    outputs, renderer = rasterize_ssplat(*inputs)

    _outputs = rasterize_gsplat(*_inputs)
    _outputs = [x.cpu() for x in _outputs]

    print("test forward")
    tol = { 'atol': 1e-3, 'rtol': 1e-3 }
    check_close('rgb', outputs[0], _outputs[0], **tol)
    check_close('depth', outputs[1], _outputs[1], **tol)
    check_close('alpha', outputs[2], _outputs[2], **tol)
    print()

    if False:
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        rgb = outputs[0].detach().cpu().numpy()
        _rgb = _outputs[0].detach().cpu().numpy()
        depth = outputs[1].detach().cpu().numpy()
        _depth = _outputs[1].detach().cpu().numpy()
        alpha = outputs[2].detach().cpu().numpy()
        _alpha = _outputs[2].detach().cpu().numpy()
        ax1.imshow(rgb[0])
        ax2.imshow(alpha[0])
        # ax2.imshow(depth[0])
        # ax2.imshow(np.log10(np.abs(depth[0]-_depth[0])))
        # ax2.imshow(rgb[1])
        # ax3.imshow(rgb[2])
        # ax4.imshow(rgb[3])
        ax3.imshow(_rgb[0])
        ax4.imshow(_alpha[0])
        # ax4.imshow(_depth[0])
        # ax4.imshow(_rgb[1])
        plt.show()
        # plt.savefig("plot.png")
        exit(0)

    weights = [torch.randn_like(x.detach()) for x in _outputs]
    if WITH_UT:
        weights[1] *= 0.0  # depth
    def fun(outputs):
        return sum([(w*o).sum() for w, o in zip(weights, outputs)])
    fun(_outputs).backward()

    renderer.zero_grad()
    renderer.backward((*weights[:-1], None), -weights[-1])
    grads = renderer.v_splats_world

    # print(inputs[0].grad)
    # print(_inputs[0].grad)

    print("test backward")
    tol = { 'atol': 5e-2, 'rtol': 5e-2 }
    check_close('means', grads[0], _inputs[0].grad, **tol)
    check_close('quats', grads[1], _inputs[1].grad, **tol)
    check_close('scales', grads[2], _inputs[2].grad, **tol)
    check_close('opacs', grads[3].squeeze(-1), _inputs[3].grad, **tol)
    check_close('features_dc', grads[4], _inputs[4].grad, **tol)
    check_close('features_sh', grads[5].reshape(N, -1, 3), _inputs[5].grad, **tol)
    # viewmats gradient is not produced by the engine 3dgs / mip backward
    # (projection bwd is called with v_viewmats = null), so it is not compared.
    print()

    _C.engine_reset()

@pytest.mark.skipif(not torch.cuda.is_available(), reason="No CUDA device")
def profile_rasterization():

    # ssplat goes through the C++ engine, whose Python boundary is CPU tensors;
    # gsplat runs fully on CUDA.
    cpu_inputs = [x.cpu().contiguous() for x in get_inputs()]
    cuda_inputs = get_inputs()

    # --- Forward ---------------------------------------------------------
    # "full" includes the set_data H2D + copy_render_to_host D2H that the
    # engine boundary performs; "compute only" times just engine_forward_3dgs
    # after state is uploaded once (no host<->device transfer).
    print("profile forward")
    _, renderer = rasterize_ssplat(*cpu_inputs)  # primes engine state
    timeit(
        lambda: _C.engine_forward_3dgs(
            renderer.primitive, renderer.sh_degree_to_use, renderer.packed),
        "ssplat forward",
    )
    timeit(lambda: rasterize_gsplat(*cuda_inputs), "gsplat forward")
    print()

    # --- Backward --------------------------------------------------------
    # Output cotangents are built once on-device, so the timed region is the
    # engine raster + projection backward (the cotangent copy is a cheap d2d,
    # not a host<->device transfer).
    print("profile backward")
    C, H, W = renderer.viewmats.shape[0], renderer.height, renderer.width
    v_rgb   = torch.randn(C, H, W, 3, device=device)
    v_depth = torch.randn(C, H, W, 1, device=device)
    v_Ts    = torch.randn(C, H, W, 1, device=device)
    timeit(
        lambda: _C.engine_backward_from_render_grad(
            renderer._tv(v_rgb), renderer._tv(v_depth), renderer._tv(v_Ts)),
        "ssplat backward",
    )

    _outputs = rasterize_gsplat(*cuda_inputs)
    weights = [torch.randn_like(x.detach()) for x in _outputs]
    _loss = sum([(w*o).sum() for w, o in zip(weights, _outputs)])
    timeit(lambda: _loss.backward(retain_graph=True), "gsplat backward")
    print()

    _C.engine_reset()



if __name__ == "__main__":
    import itertools
    import sys

    # `python3 tests/test_rasterization.py scenes [ut]` -> benchmark the
    # representative scenes instead of the correctness/profile suite.
    if "scenes" in sys.argv:
        PACKED = True; IS_FISHEYE = False
        IS_ANTIALIASED = True; WITH_UT = False
        profile_scenes()
        if "ut" in sys.argv:
            IS_ANTIALIASED = False; WITH_UT = True
            profile_scenes()
        sys.exit(0)

    # packed x fisheye x (not_aa+no_ut, aa+no_ut, not_aa+with_ut)
    modes = [(False, False), (True, False), (False, True)]
    # modes = [(False, True)]
    combos = list(itertools.product([True, False], [False, True], modes))
    # modes = [(False, False), (True, False)]
    # combos = list(itertools.product([True, False], [False], modes))

    # combos = [(True, False, (False, True))]
    # combos = [(True, False, (True, False))]
    combos = [(True, False, (True, False)), (True, False, (False, True))]

    for packed_val, fisheye_val, (aa_val, ut_val) in combos:
        PACKED = packed_val
        IS_FISHEYE = fisheye_val
        IS_ANTIALIASED = aa_val
        WITH_UT = ut_val
        label = f"packed={PACKED} fisheye={IS_FISHEYE} antialiased={IS_ANTIALIASED} with_ut={WITH_UT}"
        print(f"=== {label} ===")

        N = 2000
        try:
            test_rasterization()
        except Exception as e:
            import traceback
            traceback.print_exc()

        N = 200000
        try:
            profile_rasterization()
        except Exception as e:
            import traceback
            traceback.print_exc()

