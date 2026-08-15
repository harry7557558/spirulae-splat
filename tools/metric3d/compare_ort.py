#!/usr/bin/env python3
"""Compare src/metric3d/'s forward pass against onnxruntime, stage by stage.

The checkpoint is the same file both sides read, so a disagreement here is
ours. Run the C++ side with SS_METRIC3D_DUMP=<dir> to get one .npy per stage,
then point this at the same input.

    pip install onnx onnxruntime numpy
    SS_METRIC3D_DUMP=/tmp/ours ./build/metric3d_test --model metric3d-vit-large \
        --image street.jpg --max-size 616
    python3 tools/metric3d/compare_ort.py --onnx ~/.cache/spirula-studio/models/\
metric3d-vit-large-fp16.onnx --ours /tmp/ours

`--list` prints the graph's tensor names when a new stage needs mapping.

Run the C++ side with SS_METRIC3D_F32_WEIGHTS=1 against an fp32 export for the
tight check: the default tolerance is what fp32 accumulation order costs over
a whole forward pass, and our f16 weights are ~20x that on their own.
"""
import argparse
import os
import re
import sys

import numpy as np

# Our stage name -> the ONNX tensor that should hold the same numbers.
# `{}` is filled with the encoder's last block index, which differs per variant.
# `chw` marks a tensor onnxruntime returns channel-first where ours is
# channel-last -- everything that has become a feature map by the time it is
# read. Token sequences are [N, C] on both sides and need no transpose.
CHW = {"state0", "net0", "net1", "net2", "read1", "read0", "upconv3", "upconv2", "ref_feat", "normal"}

STAGES = {
    "tokens_in":      "/depth_model/encoder/Concat_8_output_0",
    "tokens_block0":  "/depth_model/encoder/blocks.0/blocks.0.0/Add_1_output_0",
    "tokens":         "/depth_model/encoder/norm/Add_1_output_0",
    "tokens_tap0":    "/depth_model/encoder/blocks.0/blocks.0.{t0}/Add_1_output_0",
    "tokens_tap1":    "/depth_model/encoder/blocks.0/blocks.0.{t1}/Add_1_output_0",
    "tokens_tap2":    "/depth_model/encoder/blocks.0/blocks.0.{t2}/Add_1_output_0",
    "tokens_tap3":    "/depth_model/encoder/blocks.0/blocks.0.{t3}/Add_1_output_0",
    "state0":         "/depth_model/decoder/Sub_4_output_0",
    "net0":           "/depth_model/decoder/Tanh_output_0",
    "net1":           "/depth_model/decoder/Tanh_1_output_0",
    "net2":           "/depth_model/decoder/Tanh_2_output_0",
    "read3":          "/depth_model/decoder/token2feature/read_3/readoper/act/Mul_1_output_0",
    "read2":          "/depth_model/decoder/token2feature/read_2/readoper/act/Mul_1_output_0",
    "read1":          "/depth_model/decoder/token2feature/read_1/sample/ConvTranspose_output_0",
    "read0":          "/depth_model/decoder/token2feature/read_0/sample/sample.0/Conv_output_0",
    "upconv3":        "/depth_model/decoder/decoder_mono/upconv_3/out_conv/Conv_output_0",
    "upconv2":        "/depth_model/decoder/decoder_mono/upconv_2/out_conv/Conv_output_0",
    "ref_feat":       "/depth_model/decoder/decoder_mono/upconv_1/out_conv/Conv_output_0",
    "depth":          "predicted_depth",
    "normal":         "predicted_normal",
    "normal_conf":    "normal_confidence",
}


def load_graph(path):
    import onnx
    return onnx.load(path, load_external_data=False)


def encoder_taps(model):
    """(t0..t3) block indices, or None when the export ends in encoder.norm."""
    if any("encoder/norm/" in (n.name or "") for n in model.graph.node):
        return None
    idx = set()
    for n in model.graph.node:
        m = re.search(r"/blocks\.0\.(\d+)/", n.name or "")
        if m:
            idx.add(int(m.group(1)))
    depth = max(idx) + 1
    step = depth // 4
    return [step * (i + 1) - 1 for i in range(4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--onnx", required=True)
    ap.add_argument("--ours", help="directory of .npy files written by SS_METRIC3D_DUMP")
    ap.add_argument("--list", action="store_true", help="print every tensor name")
    ap.add_argument("--tol", type=float, default=5e-3)
    args = ap.parse_args()

    import onnx
    import onnxruntime as ort

    model = load_graph(args.onnx)
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
    pixel_values = np.load(inp)

    taps = encoder_taps(model)
    subs = {}
    if taps:
        subs = {f"t{i}": taps[i] for i in range(4)}

    ours = {}
    for f in sorted(os.listdir(args.ours)):
        if f.endswith(".npy") and f != "input.npy":
            ours[f[:-4]] = np.load(os.path.join(args.ours, f))
    if not ours:
        print(f"no stage dumps in {args.ours}", file=sys.stderr)
        return 2

    # Re-export with the wanted intermediates promoted to graph outputs. ORT
    # will not fold anything a graph output depends on, so this changes the
    # numbers only where fusion would have changed them anyway.
    wanted = {}
    for name in ours:
        spec = STAGES.get(name)
        if spec is None:
            print(f"  {name}: no ONNX tensor mapped for this stage; skipping")
            continue
        try:
            wanted[name] = spec.format(**subs)
        except KeyError:
            continue

    produced = {o for n in model.graph.node for o in n.output}
    graph_outputs = {o.name for o in model.graph.output}
    extra = [t for t in wanted.values() if t in produced and t not in graph_outputs]
    for t in extra:
        model.graph.output.append(onnx.ValueInfoProto(name=t))

    so = ort.SessionOptions()
    so.intra_op_num_threads = os.cpu_count() or 4
    if any(t.data_location == onnx.TensorProto.EXTERNAL
           for t in model.graph.initializer):
        # An external-data export keeps its weights in a sibling file named by
        # a RELATIVE path, so the augmented copy has to sit in the same
        # directory or the session cannot find them.
        tmp = os.path.join(os.path.dirname(os.path.abspath(args.onnx)),
                           "_compare_ort_tmp.onnx")
        onnx.save(model, tmp)
        try:
            sess = ort.InferenceSession(tmp, so, providers=["CPUExecutionProvider"])
        finally:
            os.remove(tmp)
    else:
        sess = ort.InferenceSession(model.SerializeToString(), so,
                                    providers=["CPUExecutionProvider"])
    names = [o.name for o in sess.get_outputs()]
    # An fp16 export takes an fp16 input. Feeding it the same numbers is what
    # makes the comparison tight: its weights are then bit-identical to the
    # ones we uploaded, so what is left is our arithmetic and not the rounding
    # of the checkpoint.
    if sess.get_inputs()[0].type == "tensor(float16)":
        pixel_values = pixel_values.astype(np.float16)
    outs = dict(zip(names, sess.run(None, {"pixel_values": pixel_values})))

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
        if name in CHW and b.ndim >= 2:
            # ours is [H, W, C]; the reference is [C, H, W]
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
