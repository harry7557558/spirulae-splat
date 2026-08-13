#!/usr/bin/env python3
"""Render assets/icon.png, the application icon.

The shape is the fractal spirula from one ofharry7557558/spirulae's
examples -- the u=1 level set of a logarithmic-spiral chart, iterated
the way a Mandelbulb iterates a sphere -- in the same turbo orbit-trap colour.

Run it when the icon changes. assets/icon.png is committed and
tools/package_macos.sh reads that, so neither the build nor packaging needs
what this imports.

    pip install taichi
    python3 tools/make_icon.py
"""

import os

import numpy as np
import taichi as ti
from PIL import Image

# ---- canvas -------------------------------------------------------------
SIZE = 1024
SS = 8               # supersampling: the whorls are finer than a pixel
INSET = 100          # macOS: 824 of artwork inside a 1024 canvas
RADIUS = 185         # corner radius of that 824 square

# ---- camera -------------------------------------------------------------
# Nearly down the spiral axis, from the side the whorls open away from: only
# from here do they stack into concentric blades around the pole, the reading
# worth having in an icon. The near side is a featureless dome. Well off the
# axis, so the disc is an ellipse and the shell has a thickness.
CAM_DIST = 5.0
CAM_TILT_DEG = 140.0
CAM_AZIM_DEG = 310.0
FOV_DEG = 60.0
KEY_DIR = (-0.5, 0.5, -0.7)

# ---- shape --------------------------------------------------------------
SP_B, SP_P, SP_S = 0.1, 1.5, 0.4
SP_NMAX, SP_K, SP_M, SP_N, SP_A = 2.0, 18.0, 6.0, 8.0, 0.5
SP_CSCALE, SP_NITER = 0.6, 20

BOUND_R = 5.0
MAX_STEPS = 12000
SHADOW_STEPS = 3000
MAX_STEP = 0.001
SURFACE_OFF = 0.002

ASSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "assets")
OUT = os.path.join(ASSETS, "icon.png")
ICO = os.path.join(ASSETS, "icon.ico")
SMALL = os.path.join(ASSETS, "icon_128.png")

# ---- README banner ------------------------------------------------------
# The same shape from the same camera, framed wide and close: the whorls run
# off all four edges, so the strip is subject rather than background. Artwork
# only -- every word on the page stays markdown, editable and translatable.
BANNER = os.path.join(ASSETS, "banner.png")
BAN_W, BAN_H = 1600, 480
BAN_FOV_DEG = 13.0        # vertical; tight enough that the strip is all shell
BAN_SHIFT = -0.30         # the pole off to the right of centre, not on it

try:
    ti.init(arch=ti.gpu, default_fp=ti.f32)
except Exception:
    ti.init(arch=ti.cpu, default_fp=ti.f32)

vec3 = ti.types.vector(3, ti.f32)
PI = 3.14159265358979

W = SIZE * SS
buf = ti.Vector.field(4, dtype=ti.f32, shape=(W, W))               # rgb + hit
ban = ti.Vector.field(4, dtype=ti.f32, shape=(BAN_H * SS, BAN_W * SS))


@ti.func
def xyz2uvw(p: vec3) -> vec3:
    r = ti.sqrt(p.x * p.x + p.y * p.y)
    a = ti.atan2(p.y, p.x)
    nf = ti.min((ti.log(ti.max(r, 1e-12)) / SP_B - a) / (2.0 * PI), SP_NMAX)
    v0 = 2.0 * PI * ti.floor(nf) + a
    v1 = 2.0 * PI * ti.ceil(nf) + a
    rc0, rc1 = ti.exp(SP_B * v0), ti.exp(SP_B * v1)
    v, rc = v0, rc0
    if ti.abs(r - rc1) < ti.abs(r - rc0):
        v, rc = v1, rc1
    d = r - rc
    l = ti.sqrt(d * d + p.z * p.z)
    return ti.Vector([l / (SP_S * ti.pow(rc, SP_P)), v, ti.atan2(d, p.z)])


@ti.func
def uvw2xyz(q: vec3) -> vec3:
    rc = ti.exp(SP_B * q.y)
    l = q.x * SP_S * ti.pow(rc, SP_P)
    r = rc + l * ti.sin(q.z)
    return ti.Vector([r * ti.cos(q.y), r * ti.sin(q.y), l * ti.cos(q.z)])


# Distance estimate, and the orbit trap that colours the surface.
@ti.func
def spirula_de(p_in: vec3):
    m = SP_M / ti.sqrt(1.0 + SP_A * SP_A)
    n = SP_N / ti.sqrt(1.0 + SP_A * SP_A)
    c = p_in * SP_CSCALE
    z = c
    dr = 1.0
    r = 0.0
    trap = 1e9
    for _ in range(SP_NITER):
        r = z.norm()
        trap = ti.min(trap, r)
        if r > 64.0:
            break
        q = xyz2uvw(z)
        u, v, w = q.x, q.y, q.z
        uk = ti.pow(u, SP_K)
        z1 = uvw2xyz(ti.Vector([uk, m * (v + SP_A * w), n * (-SP_A * v + w)])) + c

        rc, rc2 = ti.exp(SP_B * v), ti.exp(SP_B * m * v)
        uk1 = ti.pow(u, SP_K - 1.0) * ti.pow(rc2 / rc, SP_P)
        r_in = ti.max(ti.abs(rc + u * SP_S * ti.pow(rc, SP_P) * ti.sin(w)), 1e-8)
        r_out = ti.abs(rc2 + uk * SP_S * ti.pow(rc2, SP_P) * ti.sin(n * w))
        dr = ti.max(ti.max(SP_K * uk1, n * uk1), m * r_out / r_in) * dr + 1.0
        z = z1
    r = ti.max(z.norm(), 1e-8)
    return 0.5 * ti.log(r) * r / dr, trap


@ti.func
def normal_at(p: vec3) -> vec3:
    e = 1e-4
    dx = spirula_de(p + vec3([e, 0, 0]))[0] - spirula_de(p - vec3([e, 0, 0]))[0]
    dy = spirula_de(p + vec3([0, e, 0]))[0] - spirula_de(p - vec3([0, e, 0]))[0]
    dz = spirula_de(p + vec3([0, 0, e]))[0] - spirula_de(p - vec3([0, 0, e]))[0]
    return ti.Vector([dx, dy, dz]).normalized()


# Turbo, to color this fractal with.
@ti.func
def turbo(t: ti.f32) -> vec3:
    c = ti.Vector([-266.03287948, -63.82163439, -97.24397473])
    c = c * t + ti.Vector([874.85901369, 202.98596195, 447.99287080])
    c = c * t + ti.Vector([-1061.31232675, -245.26823636, -766.82678386])
    c = c * t + ti.Vector([550.44382083, 149.40813531, 604.07329733])
    c = c * t + ti.Vector([-90.46877071, -55.09180383, -203.76964931])
    c = c * t + ti.Vector([-10.15061040, 9.64581379, 9.26380342])
    c = c * t + ti.Vector([2.94711014, 2.06520532, 6.33792818])
    c = c * t + ti.Vector([0.14637796, 0.08594198, 0.23431523])
    return ti.max(ti.min(c, 1.0), 0.0)


# Three light levels and a rim, never a gradient: banding is what survives the
# downsample to 16 px.
@ti.func
def cel_shade(albedo: vec3, lum: ti.f32, rim: ti.f32, lit: ti.f32) -> vec3:
    band = 0.34
    if lum > 0.42:
        band = 0.68
    if lum > 0.74:
        band = 1.0
    band = ti.min(band, 0.34 + 0.66 * lit)
    col = albedo * band
    return col + (ti.Vector([1.0, 0.97, 0.90]) - col) * rim * lit


# Hard shadow, one ray, no softening: a cel shader wants the terminator sharp.
@ti.func
def shadow(p: vec3, ldir: vec3) -> ti.f32:
    t = 4.0 * SURFACE_OFF
    b = p.dot(ldir)
    tfar = -b + ti.sqrt(ti.max(b * b - (p.dot(p) - BOUND_R * BOUND_R), 0.0))
    lit = 1.0
    for _ in range(SHADOW_STEPS):
        if t > tfar:
            break
        de = spirula_de(p + ldir * t)[0]
        if de < SURFACE_OFF:
            lit = 0.0
            break
        t += ti.min(ti.max(de * 0.7, 5e-4), MAX_STEP)
    return lit


# tan_half is the VERTICAL half-FOV, so the two outputs frame the shell the
# same way down the short axis and the wide one only sees further sideways.
@ti.kernel
def render(out: ti.template(), eye: vec3, fwd: vec3, right: vec3, up: vec3,
           tan_half: ti.f32, shift: ti.f32, key_in: vec3):
    key = key_in.normalized()
    h, w = out.shape[0], out.shape[1]
    aspect = ti.cast(w, ti.f32) / ti.cast(h, ti.f32)
    for y, x in out:
        px = (2.0 * (x + 0.5) / w - 1.0 + shift) * tan_half * aspect
        py = (1.0 - 2.0 * (y + 0.5) / h) * tan_half
        rd = (fwd + right * px + up * py).normalized()

        # Bounding sphere first: outside it the estimator is meaningless.
        b = eye.dot(rd)
        disc = b * b - (eye.dot(eye) - BOUND_R * BOUND_R)
        col = ti.Vector([0.0, 0.0, 0.0])
        hit = 0.0
        if disc > 0.0:
            sq = ti.sqrt(disc)
            t = ti.max(-b - sq, 0.0)
            tfar = -b + sq
            for _ in range(MAX_STEPS):
                if t > tfar:
                    break
                p = eye + rd * t
                de, trap = spirula_de(p)
                if de < SURFACE_OFF:
                    nrm = normal_at(p)
                    if nrm.dot(rd) > 0.0:
                        nrm = -nrm
                    diff = 0.5 + 0.5 * nrm.dot(key)
                    rim = 0.45 * ti.pow(ti.max(1.0 + nrm.dot(rd), 0.0), 3.0)
                    col = cel_shade(turbo(0.5 - 0.5 * ti.cos(2.0 * trap)),
                                    diff, rim, shadow(p + nrm * 4.0 * SURFACE_OFF, key))
                    hit = 1.0
                    break
                t += ti.min(ti.max(de * 0.7, 5e-4), MAX_STEP)
        out[y, x] = ti.Vector([col.x, col.y, col.z, hit])


def camera():
    t = np.radians(CAM_TILT_DEG)
    a = np.radians(CAM_AZIM_DEG)
    eye = CAM_DIST * np.array([np.sin(t) * np.cos(a), np.sin(t) * np.sin(a), np.cos(t)])
    fwd = -eye / np.linalg.norm(eye)
    right = np.cross(fwd, [0.0, 1.0, 0.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    return eye, fwd, right, up


def ink(mask, rgb):
    """A silhouette line, drawn by eroding the hit mask."""
    m = mask > 0.5
    e = m.copy()
    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        e &= np.roll(m, (dy, dx), (0, 1))
    edge = (m & ~e)[..., None]
    return np.where(edge, rgb * 0.18, rgb)


def water(w, h, seed_x=0.72, seed_y=0.80):
    """Deep water, lit from one corner. Both outputs stand on it."""
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
    yy, xx = yy / max(h - 1, 1), xx / max(w - 1, 1)
    top = np.array([0.043, 0.055, 0.110], np.float32)
    bottom = np.array([0.075, 0.150, 0.210], np.float32)
    img = top + (bottom - top) * yy[..., None]
    glow = np.exp(-(((xx - seed_x) ** 2 + (yy - seed_y) ** 2)) / (2 * 0.35 ** 2))
    return img + glow[..., None] * np.array([0.015, 0.035, 0.045], np.float32)


def quantize(im, dither=Image.Dither.NONE):
    """256 adaptive colours: a third of the file size. Large flat areas need
    the dither -- the banner's water banded visibly without it."""
    return im.convert(mode="P", colors=256, palette=Image.Palette.ADAPTIVE,
                      dither=dither)


def rounded_mask():
    from PIL import ImageDraw
    s = 4
    m = Image.new("L", (SIZE * s, SIZE * s), 0)
    ImageDraw.Draw(m).rounded_rectangle(
        [INSET * s, INSET * s, (SIZE - INSET) * s - 1, (SIZE - INSET) * s - 1],
        radius=RADIUS * s, fill=255)
    return m.resize((SIZE, SIZE), Image.LANCZOS)


def shade(field, fov_deg, size, shift=0.0):
    """Trace into `field`, composite over the water, and downsample to `size`."""
    eye, fwd, right, up = camera()
    render(field, eye.tolist(), fwd.tolist(), right.tolist(), up.tolist(),
           float(np.tan(np.radians(fov_deg) * 0.5)), shift, KEY_DIR)
    out = field.to_numpy()
    rgb, mask = out[..., :3], out[..., 3]
    h, w = mask.shape
    rgb = ink(mask, np.where(mask[..., None] > 0.5, rgb, water(w, h)))
    img = Image.fromarray((np.clip(rgb, 0, 1) * 255 + 0.5).astype(np.uint8))
    return img.resize(size, Image.LANCZOS).convert("RGBA")


def main():
    img = shade(buf, FOV_DEG, (SIZE, SIZE))
    img.putalpha(rounded_mask())
    quantize(img).save(OUT, minimize=True)

    # Two derived forms, both committed because the build has no Python: the
    # .ico Windows compiles into the executable (src/app/app.rc), and the small
    # PNG GuiMain.cpp hands X11. Both are cropped past the macOS inset, which
    # would otherwise leave them a tenth smaller than every neighbour in a
    # taskbar.
    art = img.crop((INSET, INSET, SIZE - INSET, SIZE - INSET))
    art.save(ICO, sizes=[(16, 16), (24, 24), (32, 32), (48, 48),
                         (64, 64), (128, 128), (256, 256)])
    quantize(art.resize((128, 128), Image.LANCZOS)).save(SMALL, minimize=True)

    quantize(shade(ban, BAN_FOV_DEG, (BAN_W, BAN_H), BAN_SHIFT),
             Image.Dither.FLOYDSTEINBERG).save(BANNER, minimize=True)
    print(f"wrote {OUT}, {ICO}, {SMALL}, {BANNER} (traced at {SS}x)")


if __name__ == "__main__":
    main()
