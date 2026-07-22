"""Quick test for the Delaunay3D python binding.

Requires the extension to be built (`bash build_develop.bash`), which exposes
`spirulae_splat.csrc.delaunay3d(points, num_threads=0)`.

Usage:
    python test_delaunay3d.py [path/to/splat.ply] [num_threads] [max_points]
"""
import sys
import time
import numpy as np
import torch

import spirulae_splat.csrc as _C


def read_ply_means(path, max_points=0):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"ply"
        fmt = f.readline().split()
        assert fmt[1] == b"binary_little_endian", fmt
        nverts, props = 0, []
        in_vertex = False
        while True:
            line = f.readline().split()
            if line[0] == b"element":
                in_vertex = line[1] == b"vertex"
                if in_vertex:
                    nverts = int(line[2])
            elif line[0] == b"property" and in_vertex:
                sizes = {b"float": 4, b"float32": 4, b"double": 8, b"int": 4,
                         b"uint": 4, b"uchar": 1, b"char": 1, b"short": 2,
                         b"ushort": 2}
                props.append((line[2], sizes[line[1]]))
            elif line[0] == b"end_header":
                break
        stride = sum(s for _, s in props)
        offs = {}
        o = 0
        for name, s in props:
            offs[name] = o
            o += s
        n = nverts if not max_points else min(nverts, max_points)
        raw = np.frombuffer(f.read(n * stride), dtype=np.uint8).reshape(n, stride)
        out = np.empty((n, 3), np.float32)
        for j, name in enumerate((b"x", b"y", b"z")):
            out[:, j] = raw[:, offs[name]:offs[name] + 4].copy().view(np.float32)[:, 0]
        return out


def main():
    if len(sys.argv) >= 2:
        path = sys.argv[1]
        nthreads = int(sys.argv[2]) if len(sys.argv) >= 3 else 0
        maxp = int(sys.argv[3]) if len(sys.argv) >= 4 else 0
        pts = torch.from_numpy(read_ply_means(path, maxp))
    else:
        nthreads = 0
        torch.manual_seed(0)
        # pts = torch.rand(50000, 3)
        # pts = torch.randn(1000000, 3)
        pts = (torch.randn(1000000, 3) + 100.0).half().float() - 100.0

    n = pts.shape[0]
    t0 = time.time()
    cells, adj = _C.delaunay3d(pts, nthreads)
    dt = time.time() - t0
    print(f"points={n}  tets={cells.shape[0]}  time={dt:.3f}s  "
          f"tets/points={cells.shape[0]/n:.2f}")
    assert cells.dtype == torch.int32 and cells.shape[1] == 4
    assert int(cells.min()) >= 0 and int(cells.max()) < n
    print("OK")


if __name__ == "__main__":
    main()
