"""Minimal COLMAP binary model reader (cameras.bin / images.bin / points3D.bin),
plus a Nerfstudio `transforms.json` reader that produces the same dicts.

Many of the benchmark datasets ship their reference poses as a Nerfstudio
`transforms.json` (from `ns-process-data`, Metashape or RealityCapture) and
never as a COLMAP model, so `read_model_any` accepts either and the scoring
tools do not care which.
"""
import json
import struct
import numpy as np
from pathlib import Path

CAMERA_MODELS = {
    0: ("SIMPLE_PINHOLE", 3), 1: ("PINHOLE", 4), 2: ("SIMPLE_RADIAL", 4),
    3: ("RADIAL", 5), 4: ("OPENCV", 8), 5: ("OPENCV_FISHEYE", 8),
    6: ("FULL_OPENCV", 12), 7: ("FOV", 5), 8: ("SIMPLE_RADIAL_FISHEYE", 4),
    9: ("RADIAL_FISHEYE", 5), 10: ("THIN_PRISM_FISHEYE", 12),
    11: ("RAD_TAN_THIN_PRISM_FISHEYE", 14), 12: ("SIMPLE_DIVISION", 4),
    13: ("DIVISION", 5), 14: ("SIMPLE_FISHEYE", 4), 15: ("FISHEYE", 5),
    16: ("EUCM", 6), 17: ("EQUIRECTANGULAR", 2),
}


class R:
    def __init__(self, path):
        self.f = open(path, "rb")

    def read(self, fmt):
        n = struct.calcsize(fmt)
        return struct.unpack(fmt, self.f.read(n))

    def string(self):
        out = b""
        while True:
            c = self.f.read(1)
            if c == b"\x00" or c == b"":
                return out.decode()
            out += c


def read_cameras(path):
    r = R(path)
    (n,) = r.read("<Q")
    cams = {}
    for _ in range(n):
        cid, model_id, w, h = r.read("<iiQQ")
        name, npar = CAMERA_MODELS[model_id]
        params = np.array(r.read("<" + "d" * npar))
        cams[cid] = dict(id=cid, model=name, width=w, height=h, params=params)
    return cams


def read_images(path, with_points=True):
    r = R(path)
    (n,) = r.read("<Q")
    imgs = {}
    for _ in range(n):
        iid = r.read("<I")[0]
        q = np.array(r.read("<4d"))          # qw qx qy qz  (world -> cam)
        t = np.array(r.read("<3d"))
        cid = r.read("<I")[0]
        name = r.string()
        (npts,) = r.read("<Q")
        if with_points:
            raw = np.frombuffer(r.f.read(npts * 24), dtype=np.uint8).reshape(npts, 24) if npts else None
            xy = np.frombuffer(raw.tobytes(), dtype="<f8").reshape(npts, 3)[:, :2].copy() if npts else np.zeros((0, 2))
            p3 = np.frombuffer(raw.tobytes(), dtype="<u8").reshape(npts, 3)[:, 2].copy() if npts else np.zeros((0,), dtype=np.uint64)
        else:
            r.f.seek(npts * 24, 1)
            xy, p3 = None, None
        imgs[iid] = dict(id=iid, qvec=q, tvec=t, camera_id=cid, name=name, xys=xy, p3D=p3)
    return imgs


def read_points3D(path):
    r = R(path)
    (n,) = r.read("<Q")
    ids = np.empty(n, dtype=np.uint64)
    xyz = np.empty((n, 3))
    err = np.empty(n)
    tracklen = np.empty(n, dtype=np.int64)
    for i in range(n):
        pid = r.read("<Q")[0]
        xyz[i] = r.read("<3d")
        r.read("<3B")
        err[i] = r.read("<d")[0]
        (tl,) = r.read("<Q")
        r.f.seek(tl * 8, 1)
        ids[i] = pid
        tracklen[i] = tl
    return dict(ids=ids, xyz=xyz, error=err, track_len=tracklen)


def qvec2rotmat(q):
    w, x, y, z = q
    return np.array([
        [1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * z * w, 2 * x * z + 2 * y * w],
        [2 * x * y + 2 * z * w, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * x * w],
        [2 * x * z - 2 * y * w, 2 * y * z + 2 * x * w, 1 - 2 * x * x - 2 * y * y]])


def _finish(m):
    for im in m["images"].values():
        im["R"] = qvec2rotmat(im["qvec"]) if "R" not in im else im["R"]
        im["C"] = -im["R"].T @ im["tvec"]     # camera centre in world coords
    return m


def read_model(d, with_points=True):
    d = Path(d)
    m = dict(cameras=read_cameras(d / "cameras.bin"),
             images=read_images(d / "images.bin", with_points))
    if with_points and (d / "points3D.bin").exists():
        m["points3D"] = read_points3D(d / "points3D.bin")
    return _finish(m)


# Nerfstudio camera_model -> (COLMAP name, ordered param keys). Only the
# parameters we compare on (focal, principal point) matter downstream, so the
# distortion tail is filled with whatever keys the file happens to carry.
_NS_MODELS = {
    "SIMPLE_PINHOLE": ("SIMPLE_PINHOLE", ["fl_x", "cx", "cy"]),
    "PINHOLE": ("PINHOLE", ["fl_x", "fl_y", "cx", "cy"]),
    "SIMPLE_RADIAL": ("SIMPLE_RADIAL", ["fl_x", "cx", "cy", "k1"]),
    "RADIAL": ("RADIAL", ["fl_x", "cx", "cy", "k1", "k2"]),
    "OPENCV": ("OPENCV", ["fl_x", "fl_y", "cx", "cy", "k1", "k2", "p1", "p2"]),
    "OPENCV_FISHEYE": ("OPENCV_FISHEYE", ["fl_x", "fl_y", "cx", "cy", "k1", "k2", "k3", "k4"]),
    "FULL_OPENCV": ("FULL_OPENCV", ["fl_x", "fl_y", "cx", "cy", "k1", "k2", "p1", "p2",
                                    "k3", "k4", "k5", "k6"]),
    "EQUIRECTANGULAR": ("EQUIRECTANGULAR", ["w", "h"]),
    # Metashape/RealityCapture exports spell Kannala-Brandt "FISHEYE".
    "FISHEYE": ("OPENCV_FISHEYE", ["fl_x", "fl_y", "cx", "cy", "k1", "k2", "k3", "k4"]),
    "EQUIDISTANT": ("OPENCV_FISHEYE", ["fl_x", "fl_y", "cx", "cy", "k1", "k2", "k3", "k4"]),
}

# Nerfstudio/Blender camera axes (x right, y up, z back) -> COLMAP (x right,
# y down, z forward). Right-multiplying the camera-to-world matrix by this
# flips the two camera axes; the world frame is untouched.
_NS_FLIP = np.diag([1.0, -1.0, -1.0])


def read_transforms(path, with_points=False):
    """Read a Nerfstudio `transforms.json` into the read_model() dicts.

    `transform_matrix` is camera-to-world in the OpenGL/Blender convention.
    `applied_transform`, when present, is the rigid world transform
    ns-process-data applied on import; it is deliberately *not* undone --
    every metric downstream is either alignment-free (relative pose, which a
    world rotation leaves untouched) or Sim(3)-aligned.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "transforms.json"
    d = json.loads(path.read_text())
    frames = d["frames"]

    cams, imgs = {}, {}
    cam_key = {}
    for i, fr in enumerate(frames):
        g = lambda k, dflt=None: fr.get(k, d.get(k, dflt))  # noqa: E731
        name = g("camera_model", "OPENCV")
        model, keys = _NS_MODELS.get(name, _NS_MODELS["OPENCV"])
        w, h = int(g("w", 0)), int(g("h", 0))
        params = tuple(float(g(k, w if k == "w" else h if k == "h" else 0.0)) for k in keys)
        ck = (model, w, h, params)
        if ck not in cam_key:
            cid = len(cam_key) + 1
            cam_key[ck] = cid
            cams[cid] = dict(id=cid, model=model, width=w, height=h,
                             params=np.array(params))
        c2w = np.array(fr["transform_matrix"], dtype=float)
        R_c2w = c2w[:3, :3] @ _NS_FLIP
        R = R_c2w.T                      # world -> camera
        t = -R @ c2w[:3, 3]
        imgs[i + 1] = dict(id=i + 1, qvec=None, tvec=t, camera_id=cam_key[ck],
                           name=fr["file_path"], xys=None, p3D=None, R=R)
    return _finish(dict(cameras=cams, images=imgs))


def read_metashape(path, largest_component_only=True):
    """Read an Agisoft Metashape camera export into the read_model() dicts.

    Several captures -- the 360 equirectangular sets in particular -- ship a
    Metashape `.xml` as their only reference. Metashape's per-component
    transforms place components in a common frame, but that placement is only
    as good as its component alignment, so by default we keep the largest
    component: an alignment-free relative-pose metric over frames Metashape
    itself never tied together would measure its guesswork, not ours.
    """
    import metashape_ref  # local import: only this path needs the XML parser
    ref = metashape_ref.load(path)
    keys = list(ref.poses)
    if largest_component_only and ref.comps:
        keys = next(iter(ref.comps.values()))

    cams, imgs = {}, {}
    for sid, s in ref.sensors.items():
        cams[len(cams) + 1] = dict(id=len(cams) + 1, model="OPENCV", width=s["width"],
                                   height=s["height"],
                                   params=np.array([s["f"], s["f"], s["cx"], s["cy"],
                                                    *s["k"][:2], *s["p"]]))
    sensor_cam = {sid: i + 1 for i, sid in enumerate(ref.sensors)}
    for i, k in enumerate(sorted(keys)):
        p = ref.poses[k]
        imgs[i + 1] = dict(id=i + 1, qvec=None, tvec=p["t"], name=k,
                           camera_id=sensor_cam.get(p["sensor"], 1),
                           xys=None, p3D=None, R=p["R"])
    return _finish(dict(cameras=cams, images=imgs))


def read_model_any(path, with_points=True):
    """A COLMAP model directory, a Nerfstudio transforms.json, a Metashape
    .xml, or a directory holding one of the latter two."""
    p = Path(path)
    if p.suffix == ".json":
        return read_transforms(p)
    if p.suffix == ".xml":
        return read_metashape(p)
    if (p / "cameras.bin").exists():
        return read_model(p, with_points)
    if (p / "transforms.json").exists():
        return read_transforms(p / "transforms.json")
    return read_model(p, with_points)  # let it raise with the COLMAP path
