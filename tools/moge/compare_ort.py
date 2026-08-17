#!/usr/bin/env python3
"""Compare src/moge/'s forward pass against onnxruntime, stage by stage.

The checkpoint is the same file both sides read, so a disagreement here is
ours. Run the C++ side with SS_MOGE_DUMP=<dir> to get one .npy per stage, then
point this at the same input.

    pip install onnx onnxruntime numpy
    SS_MOGE_F32_WEIGHTS=1 SS_MOGE_DUMP=/tmp/ours ./build/moge_test \
        --model moge2-vits --image street.jpg --max-size 640 --num-tokens 1200
    python3 tools/moge/compare_ort.py --onnx <the same .onnx> --ours /tmp/ours

`--list` prints the graph's tensor names when a new stage needs mapping.

SS_MOGE_F32_WEIGHTS=1 is not optional for a tight check: the default tolerance
is what fp32 accumulation order costs over a whole forward pass, and the f16
weights we ship are an order of magnitude above that on their own.
"""
import argparse
import os
import sys

import numpy as np

# Our stage name -> the ONNX tensor that should hold the same numbers. CHW is
# where onnxruntime returns channel-first and we are channel-last; the token
# sequences and the graph's own outputs agree and need no transpose.
CHW = {"encoder_feat", "neck0", "neck1", "neck2", "neck3", "neck4",
       "points_head_raw", "normal_head_raw", "mask_head_raw"}

STAGES = {
    "tokens_in":       "/encoder/Add_output_0",
    "tap0":            "/encoder/norm/Add_1_output_0",
    "tap1":            "/encoder/norm_1/Add_1_output_0",
    "tap2":            "/encoder/norm_2/Add_1_output_0",
    "tap3":            "/encoder/norm_3/Add_1_output_0",
    "encoder_feat":    "/encoder/ReduceSum_output_0",
    "neck0":           "/neck/input_blocks.0/Conv_output_0",
    # `{n}` is the LAST residual block at that level: vit-large's neck runs two
    # where the smaller variants run one, and the level output is the second's.
    "neck1":           "/neck/res_blocks.1/res_blocks.1.{n}/Add_output_0",
    "neck2":           "/neck/res_blocks.2/res_blocks.2.{n}/Add_output_0",
    "neck3":           "/neck/res_blocks.3/res_blocks.3.{n}/Add_output_0",
    "neck4":           "/neck/Add_3_output_0",
    "points_head_raw": "/points_head/output_blocks.4/Conv_output_0",
    "normal_head_raw": "/normal_head/output_blocks.4/Conv_output_0",
    "mask_head_raw":   "/mask_head/output_blocks.4/Conv_output_0",
    "points":          "points",
    "normal":          "normal",
    "mask":            "mask",
    # vit-small's export names this output `scale` and the other two
    # `metric_scale`; nothing else about it differs.
    "scale":           "scale|metric_scale",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--onnx", required=True)
    ap.add_argument("--ours", help="directory of .npy files written by SS_MOGE_DUMP")
    ap.add_argument("--list", action="store_true", help="print every tensor name")
    ap.add_argument("--tol", type=float, default=5e-3)
    args = ap.parse_args()

    import onnx
    import onnxruntime as ort

    model = onnx.load(args.onnx, load_external_data=False)
    if args.list:
        for n in model.graph.node:
            for o in n.output:
                print(o)
        return 0

    if not args.ours:
        ap.error("--ours is required unless --list is given")

    inp = os.path.join(args.ours, "input.npy")
    if not os.path.exists(inp):
        print(f"{inp} not found -- the C++ side writes the network input it used "
              f"so both sides are fed identical numbers", file=sys.stderr)
        return 2
    image = np.load(inp).astype(np.float32)
    num_tokens = np.array(int(np.load(os.path.join(args.ours, "num_tokens.npy"))[0]),
                          dtype=np.int64)

    ours = {}
    for f in sorted(os.listdir(args.ours)):
        if f.endswith(".npy") and f not in ("input.npy", "num_tokens.npy"):
            ours[f[:-4]] = np.load(os.path.join(args.ours, f))
    if not ours:
        print(f"no stage dumps in {args.ours}", file=sys.stderr)
        return 2

    # Re-export with the wanted intermediates promoted to graph outputs. ORT
    # will not fold anything a graph output depends on, so this changes the
    # numbers only where fusion would have changed them anyway.
    produced = {o for n in model.graph.node for o in n.output}
    graph_outputs = {o.name for o in model.graph.output}
    wanted = {}
    for name in ours:
        spec = STAGES.get(name)
        if spec is None:
            print(f"  {name}: no ONNX tensor mapped for this stage; skipping")
            continue
        if "|" in spec:
            spec = next((s for s in spec.split("|")
                         if s in produced or s in graph_outputs), spec)
        if "{n}" in spec:
            last = max((i for i in range(8) if spec.format(n=i) in produced),
                       default=None)
            if last is None:
                continue
            spec = spec.format(n=last)
        if spec not in produced and spec not in graph_outputs:
            continue          # a tap this variant does not have
        wanted[name] = spec
        if spec not in graph_outputs:
            model.graph.output.append(onnx.ValueInfoProto(name=spec))
            graph_outputs.add(spec)

    so = ort.SessionOptions()
    so.intra_op_num_threads = os.cpu_count() or 4
    sess = ort.InferenceSession(model.SerializeToString(), so,
                                providers=["CPUExecutionProvider"])
    names = [o.name for o in sess.get_outputs()]
    outs = dict(zip(names, sess.run(None, {"image": image, "num_tokens": num_tokens})))

    worst = 0.0
    for name in sorted(ours):
        tensor = wanted.get(name)
        if tensor is None or tensor not in outs:
            print(f"  {name:16s} not produced by the graph")
            continue
        a = np.asarray(outs[tensor], dtype=np.float64).squeeze()
        b = np.asarray(ours[name], dtype=np.float64).squeeze()
        if a.size != b.size:
            print(f"  {name:16s} SHAPE ort {a.shape} vs ours {b.shape}")
            worst = float("inf")
            continue
        # A single-channel map squeezes to [H, W] on both sides, so only a real
        # three-dimensional one needs the transpose.
        if name in CHW and a.ndim == 3:
            b = b.reshape(a.shape[1], a.shape[2], a.shape[0]).transpose(2, 0, 1)
        else:
            b = b.reshape(a.shape)
        scale = max(float(np.abs(a).max()), 1e-6)
        err = float(np.abs(a - b).max()) / scale
        # Relative L2 is what separates accumulated rounding from a bug: mixed
        # precision drifts the worst element without moving this much.
        l2 = float(np.linalg.norm(a - b) / max(np.linalg.norm(a), 1e-9))
        worst = max(worst, err)
        flag = "ok " if err <= args.tol else "BAD"
        print(f"  {flag} {name:16s} {str(a.shape):24s} rel {err:.3e}  L2 {l2:.3e}")

    print(f"\nworst relative error {worst:.3e} (tolerance {args.tol:g})")
    return 0 if worst <= args.tol else 1


if __name__ == "__main__":
    sys.exit(main())
