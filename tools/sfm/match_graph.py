#!/usr/bin/env python3
"""Diagnose a matches.bin as a graph: components, degrees, inlier mass.

The mapper can only reconstruct what the view graph connects, so when a scene
comes back in several models (D41) the first question is whether the graph
really is disconnected -- a retrieval problem -- or whether it is connected but
too weakly for registration to cross the joins. This answers that from the
match database alone, before spending a mapping run:

    tools/sfm/match_graph.py matches.bin [--min-inliers N] [--top K]

`--min-inliers` restricts the graph to pairs with at least that many verified
inliers, which is roughly what the mapper needs from an edge to use it; sweeping
it shows how fast the graph falls apart as the edge requirement tightens.
"""
import argparse
import struct
import sys
from collections import Counter


def read_matches(path, verified_only=True, min_inliers=0):
    """(image names, [(i, j, num_inliers)]) from a VKMT v1/v2 file."""
    with open(path, "rb") as f:
        blob = f.read()
    if blob[:4] != b"VKMT":
        sys.exit(f"{path}: not a matches.bin (bad magic)")
    off = 4
    version, nimg = struct.unpack_from("<II", blob, off)
    off += 8
    names = []
    for _ in range(nimg):
        (ln,) = struct.unpack_from("<I", blob, off)
        off += 4
        names.append(blob[off:off + ln].decode())
        off += ln + 4  # + num_features
    (npairs,) = struct.unpack_from("<I", blob, off)
    off += 4
    edges = []
    for _ in range(npairs):
        i, j, config, nm = struct.unpack_from("<IIiI", blob, off)
        off += 16 + 8 * nm
        if verified_only and config == 0:
            continue
        if nm >= min_inliers:
            edges.append((i, j, nm))
    return names, edges


def components(n, edges):
    """Connected components as a list of member lists, largest first."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, _ in edges:
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=len, reverse=True)


def pct(v, total):
    return 100.0 * v / total if total else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matches")
    ap.add_argument("--min-inliers", type=int, default=0,
                    help="only count pairs with at least this many inliers")
    ap.add_argument("--top", type=int, default=8,
                    help="how many components to list")
    ap.add_argument("--sweep", action="store_true",
                    help="report the component split at several inlier thresholds")
    args = ap.parse_args()

    names, edges = read_matches(args.matches, min_inliers=args.min_inliers)
    n = len(names)
    deg = Counter()
    for i, j, _ in edges:
        deg[i] += 1
        deg[j] += 1
    degs = sorted(deg.get(i, 0) for i in range(n))
    isolated = sum(1 for d in degs if d == 0)

    print(f"{args.matches}: {n} images, {len(edges)} verified pairs "
          f"(>= {args.min_inliers} inliers)")
    if degs:
        q = lambda p: degs[min(len(degs) - 1, int(p * len(degs)))]
        print(f"  degree: min {degs[0]}, p10 {q(0.10)}, median {q(0.50)}, "
              f"p90 {q(0.90)}, max {degs[-1]}, mean {2*len(edges)/n:.1f}")
        print(f"  isolated images (degree 0): {isolated} ({pct(isolated, n):.1f}%)")

    comps = components(n, edges)
    big = [c for c in comps if len(c) > 1]
    print(f"  components: {len(comps)} total, {len(big)} with >1 image")
    for k, c in enumerate(comps[:args.top]):
        # Name prefix (the folder) is usually the rig camera; showing the mix
        # tells a two-lens rig apart from a genuinely separate part of a scene.
        pre = Counter(names[i].split("/")[0] for i in c)
        pre_s = ", ".join(f"{p}:{v}" for p, v in sorted(pre.items()))
        print(f"    [{k}] {len(c)} images ({pct(len(c), n):.1f}%)  {pre_s}")
    if len(comps) > args.top:
        rest = sum(len(c) for c in comps[args.top:])
        print(f"    ... {len(comps) - args.top} more holding {rest} images")

    if args.sweep:
        print("  component split vs edge strength:")
        print("    min inliers | pairs | largest comp | comps>1 | isolated")
        for t in (0, 15, 25, 50, 100, 200, 400):
            _, e = read_matches(args.matches, min_inliers=t)
            c = components(n, e)
            d = Counter()
            for i, j, _ in e:
                d[i] += 1
                d[j] += 1
            iso = sum(1 for i in range(n) if d.get(i, 0) == 0)
            nbig = sum(1 for x in c if len(x) > 1)
            print(f"    {t:>11} | {len(e):>5} | {len(c[0]):>5} "
                  f"({pct(len(c[0]), n):>4.1f}%) | {nbig:>7} | {iso:>8}")


if __name__ == "__main__":
    main()
