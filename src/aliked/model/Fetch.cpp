#include "aliked/model/Fetch.h"

#include "nn/core/Error.h"
#include "nn/io/Fetch.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace aliked {
namespace {

// The same URL / filename / SHA-256 triples COLMAP carries in
// src/colmap/feature/resources.h. Keep them identical: the point of fetching
// from COLMAP's release is that a parity run compares implementations, not
// checkpoints.
const ModelSource kSources[] = {
    {"aliked-n16rot",
     {"aliked-n16rot.onnx",
      "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n16rot.onnx",
      "39c423d0a6f03d39ec89d3d1d61853765c2fb6a8b8381376c703e5758778a547", 2997054ull}},
    {"aliked-n32",
     {"aliked-n32.onnx",
      "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n32.onnx",
      "a077728a02d2de1a775c66df6de8cfeb7c6b51ca57572c64c680131c988c8b3c", 4205634ull}},
    {"aliked-lightglue",
     {"aliked-lightglue.onnx",
      "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-lightglue.onnx",
      "b9a5de7204648b18a8cf5dcac819f9d30de1a5961ef03756803c8b86c2dceb8d", 0ull}},
};

}  // namespace

const ModelSource* find_model_source(const std::string& id) {
    for (const ModelSource& s : kSources)
        if (id == s.id) return &s;
    return nullptr;
}

std::string model_cache_path(const ModelSource& src) {
    return nn::cached_path(src.onnx);
}

std::string sha256_file(const std::string& path) {
    return nn::sha256_file(path);
}

std::string ensure_model(const ModelSource& src) {
    return nn::ensure_file(src.onnx, "aliked");
}

std::string resolve_model(const std::string& id_or_path) {
    if (const ModelSource* src = find_model_source(id_or_path)) return ensure_model(*src);

    std::error_code ec;
    NN_CHECK(fs::exists(id_or_path, ec),
             "'%s' is neither a known model id (aliked-n16rot, aliked-n32, "
             "aliked-lightglue) nor a file that exists",
             id_or_path.c_str());
    return id_or_path;
}

}  // namespace aliked
