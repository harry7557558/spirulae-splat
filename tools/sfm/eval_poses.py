#!/usr/bin/env python3
"""Compare a reconstructed COLMAP sparse model against ground-truth poses.

    tools/sfm/eval_poses.py EST_SPARSE_DIR GT_SPARSE_DIR [--json out.json]

Images are matched by file name (basename, extension-insensitive). Two families
of metrics are reported, because they answer different questions:

* Relative pose (alignment free). For every ground-truth image pair, the
  relative rotation and the direction of the relative translation are compared
  with the reconstruction's. The pair's error is the max of the two angles, and
  AUC@t is the normalised area under the cumulative error curve up to t
  degrees -- the usual SfM/pose-estimation metric. Pairs touching an
  unregistered image count as 180 deg failures, so registration rate is folded
  into the number and models cannot win by dropping hard images.

* Absolute pose. The reconstruction is aligned to the ground truth with a
  RANSAC Sim(3) fit on camera centres (COLMAP's model_aligner does the same),
  then per-camera rotation error and centre error are reported. Centre errors
  are also given as a fraction of the scene extent, since SfM scale is
  arbitrary.

No dependency beyond numpy; the model reader is tools/sfm/colmap_io.py.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from colmap_io import read_model, read_model_any  # noqa: E402


# ---------------------------------------------------------------- utilities
def pose_auc(errors, thresholds):
    """Normalised area under the cumulative error curve (SuperGlue convention)."""
    if len(errors) == 0:
        return [0.0 for _ in thresholds]
    errors = np.sort(np.asarray(errors, dtype=float))
    recall = np.arange(1, len(errors) + 1) / len(errors)
    errors = np.r_[0.0, errors]
    recall = np.r_[0.0, recall]
    aucs = []
    for t in thresholds:
        last = np.searchsorted(errors, t, side="right")
        r = np.r_[recall[:last], recall[last - 1]]
        e = np.r_[errors[:last], t]
        aucs.append(float(np.trapz(r, x=e) / t))
    return aucs


def rot_angle_deg(R):
    """Rotation angle of a batch of rotation matrices, in degrees."""
    tr = np.trace(R, axis1=-2, axis2=-1)
    return np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0)))


def umeyama(src, dst, with_scale=True):
    """Sim(3) taking src -> dst: returns (s, R, t) with dst ~ s * R @ src + t."""
    mu_s, mu_d = src.mean(0), dst.mean(0)
    S, D = src - mu_s, dst - mu_d
    cov = D.T @ S / len(src)
    U, sig, Vt = np.linalg.svd(cov)
    E = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        E[2, 2] = -1
    R = U @ E @ Vt
    s = float((sig * np.diag(E)).sum() / (S ** 2).sum() * len(src)) if with_scale else 1.0
    t = mu_d - s * R @ mu_s
    return s, R, t


def ransac_sim3(src, dst, thresh, iters=2000, seed=0):
    """RANSAC Sim(3) src -> dst on point correspondences; refit on inliers."""
    n = len(src)
    if n < 3:
        return umeyama(src, dst) + (np.ones(n, bool),)
    rng = np.random.default_rng(seed)
    best_inl = np.zeros(n, bool)
    for _ in range(iters):
        idx = rng.choice(n, 3, replace=False)
        try:
            s, R, t = umeyama(src[idx], dst[idx])
        except np.linalg.LinAlgError:
            continue
        if not np.isfinite(s) or s <= 0:
            continue
        err = np.linalg.norm((s * (R @ src.T).T + t) - dst, axis=1)
        inl = err < thresh
        if inl.sum() > best_inl.sum():
            best_inl = inl
    if best_inl.sum() < 3:
        best_inl = np.ones(n, bool)
    s, R, t = umeyama(src[best_inl], dst[best_inl])
    # one reweighting pass: re-classify with the refined transform, refit
    err = np.linalg.norm((s * (R @ src.T).T + t) - dst, axis=1)
    inl = err < thresh
    if inl.sum() >= 3:
        s, R, t = umeyama(src[inl], dst[inl])
        best_inl = inl
    return s, R, t, best_inl


def key(name):
    return Path(name).stem.lower()


def path_key(name):
    """Extension-less relative path, lowercased, with a leading image-root
    component (`images/`, `images_4/`, `input/`) dropped.

    Multi-folder captures (a two-camera rig, a set split over six folders)
    reuse frame numbers across folders, so the stem alone is ambiguous there.
    `evaluate` tries this key first and falls back to the stem when the two
    sides spell their paths differently."""
    p = Path(name.replace("\\", "/"))
    parts = list(p.parts)
    if parts and (parts[0].lower().startswith(("images", "input", "rgb"))):
        parts = parts[1:]
    if not parts:
        parts = [p.stem]
    else:
        parts[-1] = Path(parts[-1]).stem
    return "/".join(parts).lower()


def pick_keyer(gt_names, est_names):
    """The name->key function that matches the most images. Ambiguous stems
    (the same stem twice on either side) disqualify the stem keyer, so a rig
    capture cannot silently score frame_0001 of camera1 against camera2's."""
    def score(fn):
        g = [fn(n) for n in gt_names]
        e = [fn(n) for n in est_names]
        if len(set(g)) != len(g) or len(set(e)) != len(e):
            return -1
        return len(set(g) & set(e))
    return max((path_key, key), key=score)


# ------------------------------------------------------------------ metrics
def relative_pose_errors(Rs_gt, ts_gt, Rs_es, ts_es, ok, min_baseline_frac=1e-3):
    """Per-pair (rot_deg, trans_deg) over all pairs; failures are 180 deg.

    `ok[i]` says whether image i was registered by the reconstruction.
    Pairs whose ground-truth baseline is degenerate contribute rotation only.
    """
    n = len(Rs_gt)
    i, j = np.triu_indices(n, k=1)

    # ground truth relatives
    Rij_gt = np.einsum("nab,ncb->nac", Rs_gt[j], Rs_gt[i])
    tij_gt = ts_gt[j] - np.einsum("nab,nb->na", Rij_gt, ts_gt[i])
    nb_gt = np.linalg.norm(tij_gt, axis=1)

    rot = np.full(len(i), 180.0)
    tra = np.full(len(i), 180.0)
    both = ok[i] & ok[j]

    if both.any():
        bi, bj = i[both], j[both]
        Rij_es = np.einsum("nab,ncb->nac", Rs_es[bj], Rs_es[bi])
        tij_es = ts_es[bj] - np.einsum("nab,nb->na", Rij_es, ts_es[bi])
        rot[both] = rot_angle_deg(np.einsum("nab,nac->nbc", Rij_es, Rij_gt[both]))
        ng, ne = nb_gt[both], np.linalg.norm(tij_es, axis=1)
        good = (ng > 0) & (ne > 0)
        cos = np.zeros(both.sum())
        cos[good] = np.sum(tij_gt[both][good] / ng[good, None] * tij_es[good] / ne[good, None], 1)
        t_ang = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
        # A pair with (almost) no ground-truth baseline has no defined
        # translation direction -- score it on rotation alone.
        degenerate = ng < min_baseline_frac * np.median(nb_gt)
        t_ang[degenerate | ~good] = 0.0
        tra[both] = t_ang
    return rot, tra


def evaluate(est_dir, gt_dir, thresholds=(5.0, 10.0, 20.0), seed=0, only=None):
    """`only` restricts the ground truth to a set of image names (any of the
    key() spellings). Pass it whenever the run saw a subset of the ground
    truth's images -- otherwise the missing ones count as unregistered, and
    every pair involving them counts as a 180 deg error, so the AUC reports
    coverage rather than accuracy."""
    est = read_model_any(est_dir, with_points=False)
    gt = read_model_any(gt_dir, with_points=False)

    kf = pick_keyer([im["name"] for im in gt["images"].values()],
                    [im["name"] for im in est["images"].values()])
    gt_by_key = {kf(im["name"]): im for im in gt["images"].values()}
    if only:
        keep = {kf(n) for n in only}
        gt_by_key = {k: v for k, v in gt_by_key.items() if k in keep}
        if not gt_by_key:
            raise ValueError("none of the `only` names are in the ground truth")
    est_by_key = {kf(im["name"]): im for im in est["images"].values()}
    names = sorted(gt_by_key)
    ok = np.array([k in est_by_key for k in names])

    Rs_gt = np.stack([gt_by_key[k]["R"] for k in names])
    ts_gt = np.stack([gt_by_key[k]["tvec"] for k in names])
    Cs_gt = np.stack([gt_by_key[k]["C"] for k in names])
    ident = np.eye(3)
    Rs_es = np.stack([est_by_key[k]["R"] if o else ident for k, o in zip(names, ok)])
    ts_es = np.stack([est_by_key[k]["tvec"] if o else np.zeros(3) for k, o in zip(names, ok)])
    Cs_es = np.stack([est_by_key[k]["C"] if o else np.zeros(3) for k, o in zip(names, ok)])

    res = {
        "est_model": str(est_dir), "gt_model": str(gt_dir),
        "num_gt_images": len(names),
        "num_registered": int(ok.sum()),
        "registration_rate": float(ok.mean()),
        "num_extra_images": int(len(est_by_key) - sum(k in gt_by_key for k in est_by_key)),
        "thresholds_deg": list(thresholds),
    }

    # ---- relative pose, alignment free ----
    rot, tra = relative_pose_errors(Rs_gt, ts_gt, Rs_es, ts_es, ok)
    err = np.maximum(rot, tra)
    res["relative"] = {
        "num_pairs": int(len(err)),
        "auc_max": pose_auc(err, thresholds),
        "auc_rot": pose_auc(rot, thresholds),
        "auc_trans": pose_auc(tra, thresholds),
        "median_rot_deg": float(np.median(rot)),
        "median_trans_deg": float(np.median(tra)),
        "mean_rot_deg": float(np.mean(rot)),
        "mean_trans_deg": float(np.mean(tra)),
        "frac_under_deg": {str(t): float((err < t).mean()) for t in (1, 2, 5, 10, 20)},
    }

    # ---- absolute pose after robust Sim(3) alignment ----
    if ok.sum() >= 3:
        extent = float(np.linalg.norm(Cs_gt - Cs_gt.mean(0), axis=1).max()) * 2.0
        s, Rl, tl, inl = ransac_sim3(Cs_es[ok], Cs_gt[ok], thresh=0.02 * extent, seed=seed)
        C_al = s * (Rl @ Cs_es[ok].T).T + tl
        pos_err = np.linalg.norm(C_al - Cs_gt[ok], axis=1)
        # world rotation of the aligned model: R_est_aligned = R_est @ Rl^T
        R_al = np.einsum("nab,cb->nac", Rs_es[ok], Rl)
        rot_err = rot_angle_deg(np.einsum("nab,nac->nbc", R_al, Rs_gt[ok]))
        # unregistered images are failures for the absolute AUC too
        rot_all = np.full(len(names), 180.0)
        rot_all[ok] = rot_err
        res["absolute"] = {
            "scene_extent_gt_units": extent,
            "sim3_scale": float(s),
            "sim3_inliers": int(inl.sum()),
            "auc_rot": pose_auc(rot_all, thresholds),
            "rot_deg": {"median": float(np.median(rot_err)), "mean": float(np.mean(rot_err)),
                        "p95": float(np.percentile(rot_err, 95)), "max": float(rot_err.max())},
            "pos_gt_units": {"median": float(np.median(pos_err)), "mean": float(np.mean(pos_err)),
                             "rmse": float(np.sqrt((pos_err ** 2).mean())),
                             "p95": float(np.percentile(pos_err, 95)), "max": float(pos_err.max())},
            "pos_pct_of_extent": {"median": float(100 * np.median(pos_err) / extent),
                                  "mean": float(100 * np.mean(pos_err) / extent),
                                  "rmse": float(100 * np.sqrt((pos_err ** 2).mean()) / extent),
                                  "max": float(100 * pos_err.max() / extent)},
        }

    # ---- intrinsics ----
    # Focal lengths are compared as f/width: the extractor works on a
    # downscaled copy, so the model's pixel units are not the ground truth's.
    def foc(model):
        f = []
        for c in model["cameras"].values():
            if c["model"] == "EQUIRECTANGULAR" or not c["width"]:
                continue  # spherical: no focal length to compare
            p = c["params"]
            fpx = float(p[0] if c["model"] in ("SIMPLE_PINHOLE", "SIMPLE_RADIAL") else p[:2].mean())
            f.append((fpx, fpx / c["width"]))
        return f
    fg, fe = foc(gt), foc(est)
    res["intrinsics"] = {
        "gt_cameras": [(c["model"], int(c["width"]), int(c["height"])) for c in gt["cameras"].values()],
        "est_cameras": [(c["model"], int(c["width"]), int(c["height"])) for c in est["cameras"].values()],
        "gt_focal_mean_px": float(np.mean([x[0] for x in fg])) if fg else None,
        "est_focal_mean_px": float(np.mean([x[0] for x in fe])) if fe else None,
        "gt_focal_over_width": float(np.mean([x[1] for x in fg])) if fg else None,
        "est_focal_over_width": float(np.mean([x[1] for x in fe])) if fe else None,
    }
    if fg and fe:
        g, e = np.mean([x[1] for x in fg]), np.mean([x[1] for x in fe])
        res["intrinsics"]["focal_rel_err_pct"] = float(100 * (e - g) / g)
    return res


def format_report(r):
    t = r["thresholds_deg"]
    L = []
    L.append(f"images      : {r['num_registered']}/{r['num_gt_images']} registered "
             f"({100 * r['registration_rate']:.1f}%)")
    rel = r["relative"]
    L.append(f"relative pose AUC (max of rot & trans angle, {rel['num_pairs']} pairs):")
    L.append("    " + "  ".join(f"@{int(x)}deg {100 * a:6.2f}" for x, a in zip(t, rel["auc_max"])))
    L.append("    rot-only  " + "  ".join(f"@{int(x)}deg {100 * a:6.2f}" for x, a in zip(t, rel["auc_rot"])))
    L.append("    trans-only" + "  ".join(f"@{int(x)}deg {100 * a:6.2f}" for x, a in zip(t, rel["auc_trans"])))
    L.append(f"    median rot {rel['median_rot_deg']:.4f} deg, median trans {rel['median_trans_deg']:.4f} deg")
    if "absolute" in r:
        a = r["absolute"]
        L.append("absolute pose after RANSAC Sim(3) alignment:")
        L.append(f"    rot err  median {a['rot_deg']['median']:.4f} deg, mean {a['rot_deg']['mean']:.4f}, "
                 f"p95 {a['rot_deg']['p95']:.4f}, max {a['rot_deg']['max']:.4f}")
        L.append(f"    pos err  median {a['pos_pct_of_extent']['median']:.4f}%, "
                 f"rmse {a['pos_pct_of_extent']['rmse']:.4f}% of scene extent "
                 f"(rmse {a['pos_gt_units']['rmse']:.5f} gt units)")
        L.append("    AUC(rot) " + "  ".join(f"@{int(x)}deg {100 * v:6.2f}" for x, v in zip(t, a["auc_rot"])))
    i = r["intrinsics"]
    if i.get("focal_rel_err_pct") is not None:
        L.append(f"intrinsics  : f/width {i['est_focal_over_width']:.4f} vs gt {i['gt_focal_over_width']:.4f} "
                 f"({i['focal_rel_err_pct']:+.2f}%)")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("est")
    ap.add_argument("gt")
    ap.add_argument("--json", default=None)
    ap.add_argument("--thresholds", default="5,10,20")
    ap.add_argument("--only-dir", default=None,
                    help="score only ground-truth images that exist under this "
                         "directory (use when the run saw a subset)")
    args = ap.parse_args()
    th = tuple(float(x) for x in args.thresholds.split(","))
    only = None
    if args.only_dir:
        root = Path(args.only_dir)
        only = [str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()]
    r = evaluate(args.est, args.gt, th, only=only)
    print(format_report(r))
    if args.json:
        Path(args.json).write_text(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()
