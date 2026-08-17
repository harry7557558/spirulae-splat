#pragma once
// Which MoGe-2 checkpoint to fetch, and from where.
//
// The artifacts are Ruicheng's own ONNX exports of MoGe-2, fetched
// byte-for-byte. There is no fp16 export to prefer as there is for Metric3D,
// so these are fp32 on disk and the loader rounds the matrices on the way to
// the device.
//
// The licence is MIT, so unlike the segmentation checkpoints
// (src/app/gui/ModelCache.cpp) this needs no consent gate. It still never
// downloads behind the user's back: `spirula geometry` prints what it is
// fetching and from where.

#include "nn/io/Fetch.h"

#include <string>

namespace moge {

struct ModelSource {
    const char*   id;    // "moge2-vitb" -- what --model spells
    nn::FetchFile onnx;
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

}  // namespace moge
