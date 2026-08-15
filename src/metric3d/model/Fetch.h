#pragma once
// Which Metric3D checkpoint to fetch, and from where.
//
// The artifacts are the onnx-community exports of Metric3D v2 -- the same
// weights `torch.hub.load('yvanyin/metric3d', 'metric3d_vit_large')` loads,
// exported once by someone else and fetched byte-for-byte. We host nothing and
// convert nothing, so there is no export to keep in step with upstream.
//
// Only the **fp16** exports are listed. The inference layer uploads every
// matrix weight as f16 regardless (that is what the tensor-core GEMM takes),
// so an fp32 download would be twice the bytes for numbers we round on the way
// in -- and giant2's fp32 export is 5.5 GB against 2.75. Pass a path to
// --model to run an fp32 file anyway, which is what the parity harness does.
//
// giant2 is two files: over 2 GB, protobuf's own message limit forces the
// weights out into a sibling `.onnx_data` that the reader resolves by name, so
// both have to land in the same directory.
//
// The licence is CC0-1.0 on the exports and BSD-2-Clause on Metric3D itself,
// so unlike the segmentation checkpoints (src/app/gui/ModelCache.cpp) this
// needs no consent gate. It still never downloads behind the user's back:
// `spirula geometry` prints what it is fetching and from where.

#include "nn/io/Fetch.h"

#include <string>

namespace metric3d {

struct ModelSource {
    const char*   id;    // "metric3d-vit-large" -- what --model spells
    nn::FetchFile onnx;
    // The external-data sibling, or a null `file` when the export is
    // self-contained. Fetched into the same directory and named exactly as the
    // graph's `location` string, which is the only name the reader will look
    // for.
    nn::FetchFile data;
};

// Null when `id` is not one of ours.
const ModelSource* find_model_source(const std::string& id);

// A path to a verified local copy, downloading through the system `curl` if
// needed. Throws nn::Error naming the URL when that fails.
std::string ensure_model(const ModelSource& src);

// Resolve what the user asked for: an explicit path is used as-is, otherwise
// the id is fetched. Throws with the list of known ids when it is neither.
std::string resolve_model(const std::string& id_or_path);

// Comma-separated known ids, for a usage message.
std::string model_id_list();

}  // namespace metric3d
