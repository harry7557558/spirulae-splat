#!/usr/bin/env python3
"""Compare our ALIKED extractor against COLMAP's ONNX one, on the same image.

This is the gate that says the port reproduces the reference, and it is the
reason the checkpoint is fetched from COLMAP's releases rather than converted:
both sides run the same weights, so a difference here is ours.

    # 1. COLMAP's side (needs a COLMAP built with ONNX support)
    colmap feature_extractor --database_path /tmp/db.db --image_path IMAGES \\
        --FeatureExtraction.type ALIKED_N16ROT

    # 2. ours
    ./build/aliked_test --image IMAGES/x.jpg --out /tmp/ours.bin

    # 3. compare
    python3 tools/aliked/compare_colmap.py /tmp/db.db /tmp/ours.bin

What it reports, and what each number means:

  matched      the fraction of our keypoints that have a COLMAP keypoint within
               `--tol` pixels. Detection is a threshold on a continuous score
               map, so a handful of borderline points differ between any two
               implementations; what would signal a real bug is a systematic
               miss, or an offset in the mean residual.
  offset       the mean signed (dx, dy) over matched pairs. This is the number
               that catches a half-pixel convention error, and it should be
               zero to several decimals, not merely small.
  cosine       descriptor agreement on matched pairs. Both are L2-normalized,
               so this is a dot product; anything below ~0.99 means the
               descriptor head diverges even where the detector agrees.
"""

import argparse
import sqlite3
import struct
import sys

import numpy as np


def read_ours(path):
    with open(path, "rb") as f:
        blob = f.read()
    if blob[:8] != b"ALIKEDFT":
        raise SystemExit(f"{path}: not an aliked_test dump")
    version, w, h, count, dim = struct.unpack_from("<IiiII", blob, 8)
    if version != 1:
        raise SystemExit(f"{path}: version {version}")
    off = 28
    kp = np.frombuffer(blob, np.float32, count * 3, off).reshape(count, 3)
    off += count * 3 * 4
    desc = np.frombuffer(blob, np.float32, count * dim, off).reshape(count, dim)
    return dict(width=w, height=h, xy=kp[:, :2], score=kp[:, 2], desc=desc)


def read_colmap(db_path, image_name=None):
    con = sqlite3.connect(db_path)
    rows = con.execute("SELECT image_id, name FROM images ORDER BY image_id").fetchall()
    if not rows:
        raise SystemExit(f"{db_path}: no images")
    if image_name is not None:
        rows = [r for r in rows if r[1] == image_name or r[1].endswith("/" + image_name)]
        if not rows:
            raise SystemExit(f"{db_path}: no image named {image_name}")
    image_id, name = rows[0]

    r, c, data = con.execute(
        "SELECT rows, cols, data FROM keypoints WHERE image_id=?", (image_id,)
    ).fetchone()
    kp = np.frombuffer(data, np.float32).reshape(r, c)

    r2, c2, data2 = con.execute(
        "SELECT rows, cols, data FROM descriptors WHERE image_id=?", (image_id,)
    ).fetchone()
    # ALIKED descriptors are float32 stored in COLMAP's uint8 descriptor blob
    # (feature/aliked.cc writes descriptor_dim * sizeof(float) bytes per row),
    # so the column count is 4x the descriptor width.
    desc = np.frombuffer(data2, np.uint8).reshape(r2, c2).view(np.float32)
    con.close()
    return dict(name=name, xy=kp[:, :2].copy(), desc=desc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("database")
    ap.add_argument("ours")
    ap.add_argument("--image", default=None, help="image name in the database")
    ap.add_argument("--tol", type=float, default=1.0, help="match radius, pixels")
    args = ap.parse_args()

    ours = read_ours(args.ours)
    ref = read_colmap(args.database, args.image)

    print(f"ours   : {len(ours['xy']):5d} keypoints, {ours['desc'].shape[1]}-D")
    print(f"colmap : {len(ref['xy']):5d} keypoints, {ref['desc'].shape[1]}-D  ({ref['name']})")
    if ours["desc"].shape[1] != ref["desc"].shape[1]:
        raise SystemExit("descriptor widths differ")

    # Nearest reference keypoint for each of ours. A few thousand points each
    # way, so the quadratic distance matrix is the simple thing that works.
    d = np.linalg.norm(ours["xy"][:, None, :] - ref["xy"][None, :, :], axis=2)
    nearest = d.argmin(axis=1)
    dist = d[np.arange(len(nearest)), nearest]
    ok = dist <= args.tol

    n_ok = int(ok.sum())
    print(f"\nmatched: {n_ok}/{len(ok)} ({100.0 * n_ok / max(len(ok), 1):.1f}%) "
          f"within {args.tol} px")
    if n_ok == 0:
        raise SystemExit("no keypoints matched -- check the coordinate convention")

    delta = ref["xy"][nearest[ok]] - ours["xy"][ok]
    print(f"offset : mean ({delta[:, 0].mean():+.4f}, {delta[:, 1].mean():+.4f}) px, "
          f"rms {np.sqrt((delta ** 2).sum(axis=1).mean()):.4f} px")
    print(f"         max |dx| {np.abs(delta[:, 0]).max():.4f}, "
          f"max |dy| {np.abs(delta[:, 1]).max():.4f}")

    cos = (ours["desc"][ok] * ref["desc"][nearest[ok]]).sum(axis=1)
    print(f"cosine : mean {cos.mean():.6f}, min {cos.min():.6f}, "
          f"5th pct {np.percentile(cos, 5):.6f}")
    below = int((cos < 0.99).sum())
    print(f"         {below} of {n_ok} below 0.99")

    # Not a pass/fail: what "close enough" means depends on what is being
    # changed, and a human reading the three numbers above is the gate.
    return 0


if __name__ == "__main__":
    sys.exit(main())
