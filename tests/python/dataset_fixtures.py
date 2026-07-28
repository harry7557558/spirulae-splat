"""Synthetic dataset fixtures for the dataparser parity test.

Everything is generated from a fixed seed into a temp directory: the tests must
not depend on any dataset that happens to live on the machine running them.
(`SSPLAT_TEST_DATASET` can additionally point the parity test at a real
dataset -- see test_dataparser_parity.py.)

Formats produced: COLMAP (text and binary), Nerfstudio, Metashape.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

WIDTH, HEIGHT = 64, 48
N_IMAGES = 7
N_POINTS = 40


# ---------------------------------------------------------------------------
# Shared synthetic scene
# ---------------------------------------------------------------------------

def _rng() -> np.random.Generator:
    return np.random.default_rng(20260722)


def _qvec_from_rotmat(R: np.ndarray) -> np.ndarray:
    """COLMAP-order (w, x, y, z) quaternion from a rotation matrix."""
    m = R
    t = np.trace(m)
    if t > 0:
        s = np.sqrt(t + 1.0) * 2
        w = 0.25 * s
        x = (m[2, 1] - m[1, 2]) / s
        y = (m[0, 2] - m[2, 0]) / s
        z = (m[1, 0] - m[0, 1]) / s
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        s = np.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
        w = (m[2, 1] - m[1, 2]) / s
        x = 0.25 * s
        y = (m[0, 1] + m[1, 0]) / s
        z = (m[0, 2] + m[2, 0]) / s
    elif m[1, 1] > m[2, 2]:
        s = np.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
        w = (m[0, 2] - m[2, 0]) / s
        x = (m[0, 1] + m[1, 0]) / s
        y = 0.25 * s
        z = (m[1, 2] + m[2, 1]) / s
    else:
        s = np.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2
        w = (m[1, 0] - m[0, 1]) / s
        x = (m[0, 2] + m[2, 0]) / s
        y = (m[1, 2] + m[2, 1]) / s
        z = 0.25 * s
    return np.array([w, x, y, z], dtype=np.float64)


def _look_at(eye: np.ndarray, target: np.ndarray) -> np.ndarray:
    """World-to-camera rotation, COLMAP/OpenCV convention (+Z forward)."""
    fwd = target - eye
    fwd /= np.linalg.norm(fwd)
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(up, fwd)) > 0.99:
        up = np.array([0.0, 1.0, 0.0])
    right = np.cross(fwd, up)
    right /= np.linalg.norm(right)
    down = np.cross(fwd, right)
    return np.stack([right, down, fwd], axis=0)   # rows = camera axes


def make_scene(fy: float = None) -> Dict:
    """Cameras on a ring looking at the origin, plus a point cloud.

    `fy` overrides the vertical focal length (default: equal to fx).
    """
    rng = _rng()
    angles = np.linspace(0.0, 2.0 * np.pi, N_IMAGES, endpoint=False)
    positions, rotations = [], []
    for i, a in enumerate(angles):
        # Slight per-camera jitter so the poses are not degenerate.
        eye = np.array([3.0 * np.cos(a), 3.0 * np.sin(a), 1.0 + 0.1 * i])
        R = _look_at(eye, np.zeros(3))
        positions.append(eye)
        rotations.append(R)
    positions = np.stack(positions)
    rotations = np.stack(rotations)
    # tvec = -R @ C
    tvecs = -np.einsum("nij,nj->ni", rotations, positions)

    xyz = rng.normal(scale=0.6, size=(N_POINTS, 3))
    rgb = rng.integers(0, 256, size=(N_POINTS, 3), dtype=np.uint8)

    return dict(
        positions=positions, rotations=rotations, tvecs=tvecs,
        xyz=xyz.astype(np.float64), rgb=rgb,
        names=[f"frame_{i:04d}.png" for i in range(N_IMAGES)],
        # OPENCV: fx fy cx cy k1 k2 p1 p2.
        # fx == fy on purpose: a Metashape frame calibration stores a single
        # <f>, so anisotropic focals cannot round-trip through that format and
        # the cross-format check below would compare apples to oranges.
        # test_anisotropic_focal_colmap covers fx != fy on the formats that
        # can express it.
        intrinsics=(52.0, 52.0 if fy is None else fy,
                    WIDTH / 2 - 0.5, HEIGHT / 2 + 0.5),
        distortion=(0.021, -0.0043, 0.00061, -0.00072),
    )


def make_equirect_scene(width: int = 64, height: int = 32) -> Dict:
    """The same scene shot with a COLMAP EQUIRECTANGULAR (model id 17) camera.

    Only the COLMAP writers honour "colmap_camera"; Nerfstudio and Metashape
    reach EQUIRECTANGULAR by their own routes (`camera_model` in
    transforms.json, sensor type "spherical" in the .xml).

    A spherical camera has no focal length and no distortion -- the params are
    exactly (w, h), which restate the resolution. Default 2:1, since anything
    else trips the parser's non-360x180 warning.
    """
    scene = make_scene()
    scene["size"] = (width, height)
    scene["colmap_camera"] = ("EQUIRECTANGULAR", 17,
                              (float(width), float(height)))
    return scene


def _write_images(image_dir: Path, names: List[str],
                  size: Tuple[int, int] = None) -> None:
    """Small real PNGs, so parsers that stat or open the image succeed."""
    from PIL import Image
    w, h = size if size is not None else (WIDTH, HEIGHT)
    image_dir.mkdir(parents=True, exist_ok=True)
    rng = _rng()
    for name in names:
        arr = rng.integers(0, 256, size=(h, w, 3), dtype=np.uint8)
        Image.fromarray(arr).save(image_dir / name)


# ---------------------------------------------------------------------------
# COLMAP
# ---------------------------------------------------------------------------

def _colmap_camera(scene: Dict) -> Tuple[str, int, int, int, Tuple[float, ...]]:
    """The cameras.{txt,bin} record: (model name, model id, w, h, params).

    Defaults to the OPENCV camera every golden-matrix fixture uses;
    `make_equirect_scene` overrides it via the "colmap_camera" key.
    """
    w, h = scene.get("size", (WIDTH, HEIGHT))
    model, model_id, params = scene.get(
        "colmap_camera",
        ("OPENCV", 4, (*scene["intrinsics"], *scene["distortion"])))
    return model, model_id, w, h, params


def write_colmap_text(root: Path, scene: Dict) -> Path:
    root = Path(root)
    recon = root / "sparse" / "0"
    recon.mkdir(parents=True, exist_ok=True)
    model, _, w, h, params = _colmap_camera(scene)

    with open(recon / "cameras.txt", "w") as f:
        f.write("# Camera list with one line of data per camera:\n")
        f.write("#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]\n")
        f.write(f"1 {model} {w} {h} "
                + " ".join(str(p) for p in params) + "\n")

    with open(recon / "images.txt", "w") as f:
        f.write("# Image list with two lines of data per image:\n")
        f.write("#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME\n")
        f.write("#   POINTS2D[] as (X, Y, POINT3D_ID)\n")
        for i, name in enumerate(scene["names"]):
            q = _qvec_from_rotmat(scene["rotations"][i])
            t = scene["tvecs"][i]
            f.write(f"{i + 1} {q[0]!r} {q[1]!r} {q[2]!r} {q[3]!r} "
                    f"{t[0]!r} {t[1]!r} {t[2]!r} 1 {name}\n")
            f.write("\n")   # empty POINTS2D line

    with open(recon / "points3D.txt", "w") as f:
        f.write("# 3D point list with one line of data per point:\n")
        f.write("#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[]\n")
        for i, (p, c) in enumerate(zip(scene["xyz"], scene["rgb"])):
            f.write(f"{i + 1} {p[0]!r} {p[1]!r} {p[2]!r} "
                    f"{int(c[0])} {int(c[1])} {int(c[2])} 0.5 1 0\n")

    _write_images(root / "images", scene["names"], scene.get("size"))
    return root


def write_colmap_binary(root: Path, scene: Dict) -> Path:
    root = Path(root)
    recon = root / "sparse" / "0"
    recon.mkdir(parents=True, exist_ok=True)
    _, model_id, w, h, params = _colmap_camera(scene)

    # cameras.bin -- model id 4 == OPENCV (8 params) unless overridden
    with open(recon / "cameras.bin", "wb") as f:
        f.write(struct.pack("<Q", 1))
        f.write(struct.pack("<iiQQ", 1, model_id, w, h))
        f.write(struct.pack("<%dd" % len(params), *params))

    with open(recon / "images.bin", "wb") as f:
        f.write(struct.pack("<Q", len(scene["names"])))
        for i, name in enumerate(scene["names"]):
            q = _qvec_from_rotmat(scene["rotations"][i])
            t = scene["tvecs"][i]
            f.write(struct.pack("<i", i + 1))
            f.write(struct.pack("<4d", *q))
            f.write(struct.pack("<3d", *t))
            f.write(struct.pack("<i", 1))
            f.write(name.encode() + b"\x00")
            f.write(struct.pack("<Q", 0))    # no 2D points

    with open(recon / "points3D.bin", "wb") as f:
        f.write(struct.pack("<Q", len(scene["xyz"])))
        for i, (p, c) in enumerate(zip(scene["xyz"], scene["rgb"])):
            f.write(struct.pack("<Q", i + 1))
            f.write(struct.pack("<3d", *p))
            f.write(struct.pack("<3B", int(c[0]), int(c[1]), int(c[2])))
            f.write(struct.pack("<d", 0.5))
            f.write(struct.pack("<Q", 1))         # track length
            f.write(struct.pack("<ii", 1, 0))     # (image_id, point2D_idx)

    _write_images(root / "images", scene["names"], scene.get("size"))
    return root


# ---------------------------------------------------------------------------
# Nerfstudio
# ---------------------------------------------------------------------------

def _colmap_to_nerfstudio_c2w(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """COLMAP world->cam (R, t) to nerfstudio camera-to-world (OpenGL axes)."""
    c2w = np.eye(4)
    c2w[:3, :3] = R.T * np.array([[1.0, -1.0, -1.0]])
    c2w[:3, 3] = -(R.T @ t)
    return c2w


def write_nerfstudio(root: Path, scene: Dict) -> Path:
    import json
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    fx, fy, cx, cy = scene["intrinsics"]
    k1, k2, p1, p2 = scene["distortion"]

    frames = []
    for i, name in enumerate(scene["names"]):
        c2w = _colmap_to_nerfstudio_c2w(scene["rotations"][i], scene["tvecs"][i])
        frames.append({
            "file_path": f"images/{name}",
            "transform_matrix": c2w.tolist(),
            "fl_x": fx, "fl_y": fy, "cx": cx, "cy": cy,
            "w": WIDTH, "h": HEIGHT,
            "camera_model": "OPENCV",
            "k1": k1, "k2": k2, "p1": p1, "p2": p2,
        })
    meta = {"frames": frames, "ply_file_path": "points3d.ply"}
    with open(root / "transforms.json", "w") as f:
        json.dump(meta, f, indent=1)

    write_ply(root / "points3d.ply", scene["xyz"], scene["rgb"])
    _write_images(root / "images", scene["names"])
    return root


def write_ply(path: Path, xyz: np.ndarray, rgb: np.ndarray) -> None:
    """Binary little-endian PLY with float x/y/z + uchar red/green/blue."""
    n = len(xyz)
    header = (
        "ply\n"
        "format binary_little_endian 1.0\n"
        f"element vertex {n}\n"
        "property float x\nproperty float y\nproperty float z\n"
        "property uchar red\nproperty uchar green\nproperty uchar blue\n"
        "end_header\n"
    )
    with open(path, "wb") as f:
        f.write(header.encode())
        for p, c in zip(xyz, rgb):
            f.write(struct.pack("<3f", *[float(v) for v in p]))
            f.write(struct.pack("<3B", int(c[0]), int(c[1]), int(c[2])))


# ---------------------------------------------------------------------------
# Metashape
# ---------------------------------------------------------------------------

def write_metashape(root: Path, scene: Dict) -> Path:
    """Metashape camera-export XML + PLY.

    Metashape stores camera-to-world transforms in its own chunk frame with
    +Y down / +Z forward, i.e. the COLMAP cam-to-world with no axis flip.
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    fx, fy, cx, cy = scene["intrinsics"]
    k1, k2, p1, p2 = scene["distortion"]

    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<document version="1.5.0">',
             '  <chunk label="Chunk 1" enabled="true">',
             '    <sensors next_id="1">',
             '      <sensor id="0" label="synthetic" type="frame">',
             f'        <resolution width="{WIDTH}" height="{HEIGHT}"/>',
             '        <calibration type="frame" class="adjusted">',
             f'          <resolution width="{WIDTH}" height="{HEIGHT}"/>',
             f'          <f>{fx}</f>',
             f'          <cx>{cx - WIDTH / 2}</cx>',
             f'          <cy>{cy - HEIGHT / 2}</cy>',
             f'          <k1>{k1}</k1>',
             f'          <k2>{k2}</k2>',
             f'          <p1>{p1}</p1>',
             f'          <p2>{p2}</p2>',
             '        </calibration>',
             '      </sensor>',
             '    </sensors>',
             f'    <cameras next_id="{N_IMAGES}" next_group_id="0">']
    for i, name in enumerate(scene["names"]):
        # Metashape transform = camera-to-world, row-major 4x4.
        R = scene["rotations"][i]
        t = scene["tvecs"][i]
        c2w = np.eye(4)
        c2w[:3, :3] = R.T
        c2w[:3, 3] = -(R.T @ t)
        flat = " ".join(repr(float(v)) for v in c2w.reshape(-1))
        lines += [f'      <camera id="{i}" sensor_id="0" component_id="0" '
                  f'label="{Path(name).stem}">',
                  f'        <transform>{flat}</transform>',
                  '      </camera>']
    camera_ids = " ".join(str(i) for i in range(N_IMAGES))
    lines += ['    </cameras>',
              '    <components next_id="1" active_id="0">',
              '      <component id="0" label="Component 1">',
              # Real exports always carry this; without it the Python
              # metashape_utils builds an empty component list and raises
              # UnboundLocalError on `data`.
              f'        <camera_ids>{camera_ids}</camera_ids>',
              '        <transform>',
              '          <rotation>1 0 0 0 1 0 0 0 1</rotation>',
              '          <translation>0 0 0</translation>',
              '          <scale>1</scale>',
              '        </transform>',
              '      </component>',
              '    </components>',
              '  </chunk>',
              '</document>']
    (root / "cameras.xml").write_text("\n".join(lines) + "\n")

    write_ply(root / "points.ply", scene["xyz"], scene["rgb"])
    _write_images(root / "images", scene["names"])
    return root
