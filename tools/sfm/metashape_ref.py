#!/usr/bin/env python3
"""Read an Agisoft Metashape camera-export .xml as a pose reference.

Same parse as `src/data/parsers/MetashapeParser.cpp`, reduced to what a pose
comparison needs: sensors, per-camera world poses, and the component grouping.
`colmap_io.read_metashape` is the only caller.

Camera labels are not file names -- a two-lens rig writes `00017` for both
`cam0/00017.jpg` and `cam1/00017.jpg`. The `.psx` project's `.files` tree holds
a camera_id -> photo-path table that resolves this; without it we fall back to
matching the label as a substring of the image paths, and drop the ambiguous
ones rather than scoring a pose against the wrong frame.
"""
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

IMG_EXT = {".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp"}


class Ref:
    """poses: name -> {R (world->cam), t, sensor}; sensors: id -> intrinsics;
    comps: component id -> [name, ...]."""

    def __init__(self, sensors, poses, comps):
        self.sensors, self.poses, self.comps = sensors, poses, comps


# ------------------------------------------------------------------ helpers
def _floats(text):
    return [float(x) for x in re.split(r"\s+", (text or "").strip()) if x]


def _child_float(node, tag):
    c = node.find(tag)
    return float(c.text) if c is not None and c.text else 0.0


def _camera_table(psx):
    """camera_id -> photo path, from the .psx project's .files tree."""
    root = ET.parse(psx).getroot()
    rel = (root.get("path") or "").replace("{projectname}", psx.stem)
    if not rel:
        return {}
    files_root = (psx.parent / rel).parent

    def from_xml_root(r):
        cams = r.find("cameras")
        if cams is None:
            return {}
        d = {}
        for cam in cams.findall("camera"):
            photo = cam.find("photo")
            if cam.get("camera_id") and photo is not None and photo.get("path"):
                d[cam.get("camera_id")] = photo.get("path")
        return d

    if not files_root.is_dir():
        return {}
    for p in sorted(files_root.rglob("*")):
        if not p.is_file():
            continue
        try:
            if p.suffix.lower() == ".xml":
                d = from_xml_root(ET.parse(p).getroot())
                if d:
                    return d
            elif p.suffix.lower() == ".zip":
                with zipfile.ZipFile(p) as z:
                    for n in z.namelist():
                        if not n.lower().endswith(".xml"):
                            continue
                        d = from_xml_root(ET.fromstring(z.read(n)))
                        if d:
                            return d
        except Exception:
            continue
    return {}


def _image_names(image_root):
    if image_root is None or not Path(image_root).is_dir():
        return []
    root = Path(image_root)
    return [p.relative_to(root).as_posix() for p in root.rglob("*")
            if p.is_file() and p.suffix.lower() in IMG_EXT]


def _match_suffix(names, path):
    """Names sharing the longest common suffix with `path`."""
    q = path.replace("\\", "/")
    best, out = 0, []
    for n in names:
        i = 0
        while i < min(len(q), len(n)) and q[-1 - i] == n[-1 - i]:
            i += 1
        if i > best:
            best, out = i, []
        if i == best and i > 0:
            out.append(n)
    return out


def _match_substring(names, label):
    q = label.replace("\\", "/")
    for ln in range(len(q), 0, -1):
        hits = [n for n in names if q[:ln] in n]
        if hits:
            return hits
    return []


# -------------------------------------------------------------------- entry
def load(path, image_root=None):
    path = Path(path)
    if image_root is None:
        for cand in ("images", "input", "omni"):
            if (path.parent / cand).is_dir():
                image_root = path.parent / cand
                break
    names = _image_names(image_root)

    table = {}
    psx = sorted(path.parent.glob("*.psx"))
    if psx:
        table = _camera_table(psx[0])

    chunk = ET.parse(path).getroot()[0]

    sensors = {}
    for s in chunk.find("sensors"):
        res = s.find("resolution")
        if res is None:
            continue
        w, h = int(res.get("width")), int(res.get("height"))
        calib = next((c for c in s.findall("calibration")
                      if c.get("class") and c.get("class") != "initial"), None)
        if calib is None:
            # spherical / cylindrical: the image is the calibration
            sensors[s.get("id")] = dict(width=w, height=h, f=w / 2.0,
                                        cx=w / 2.0, cy=h / 2.0, k=[0.0] * 4,
                                        p=[0.0, 0.0], type=s.get("type"))
            continue
        f = float(calib.find("f").text)
        sensors[s.get("id")] = dict(
            width=w, height=h, f=f,
            cx=_child_float(calib, "cx") + w / 2.0,
            cy=_child_float(calib, "cy") + h / 2.0,
            k=[_child_float(calib, k) for k in ("k1", "k2", "k3", "k4")],
            p=[_child_float(calib, "p2"), _child_float(calib, "p1")],
            type=s.get("type"))

    comp_xform = {}
    comps_el = chunk.find("components")
    if comps_el is not None:
        for comp in comps_el.iter("component"):
            tr = comp.find("transform")
            if tr is None:
                continue
            m = np.eye(4)
            rot = tr.find("rotation")
            if rot is not None and _floats(rot.text):
                m[:3, :3] = np.array(_floats(rot.text)[:9]).reshape(3, 3)
            sc = tr.find("scale")
            if sc is not None and sc.text:
                m[:3, :3] *= float(sc.text)
            t = tr.find("translation")
            if t is not None and _floats(t.text):
                m[:3, 3] = _floats(t.text)[:3]
            comp_xform[comp.get("id")] = m

    poses, comps = {}, {}
    for cam in chunk.find("cameras").iter("camera"):
        cid, tr = cam.get("id"), cam.find("transform")
        if cid is None or tr is None or len(_floats(tr.text)) < 16:
            continue
        if cam.get("sensor_id") not in sensors:
            continue
        if names:
            if table:
                hits = _match_suffix(names, table.get(cid, "")) if cid in table else []
            else:
                hits = _match_substring(names, cam.get("label") or "")
            if len(hits) != 1:
                continue
            name = hits[0]
        else:
            name = cam.get("label") or cid

        m = np.array(_floats(tr.text)[:16]).reshape(4, 4)
        comp_id = cam.get("component_id")
        if comp_id in comp_xform:
            m = comp_xform[comp_id] @ m
            det = np.linalg.det(m[:3, :3])
            if det != 0.0:
                m[:3, :3] /= np.cbrt(det)

        R = m[:3, :3].T                       # world -> camera
        poses[name] = dict(R=R, t=-R @ m[:3, 3], sensor=cam.get("sensor_id"))
        comps.setdefault(comp_id or "0", []).append(name)

    comps = dict(sorted(comps.items(), key=lambda kv: -len(kv[1])))
    return Ref(sensors, poses, comps)
