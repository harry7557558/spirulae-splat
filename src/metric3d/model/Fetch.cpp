#include "metric3d/model/Fetch.h"

#include "nn/core/Error.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace metric3d {
namespace {

// SHA-256 is the LFS object id each file is published under, so a mismatch
// means the download broke rather than that upstream moved.
const ModelSource kSources[] = {
    {"metric3d-vit-small",
     {"metric3d-vit-small-fp16.onnx",
      "https://huggingface.co/onnx-community/metric3d-vit-small/resolve/main/onnx/model_fp16.onnx",
      "4afcc0893dbb3c0c63e270f8bb24bfa63ccf2dc68ab9c0c3601fdae4f0aafd9b", 75778144ull},
     {}},
    {"metric3d-vit-large",
     {"metric3d-vit-large-fp16.onnx",
      "https://huggingface.co/onnx-community/metric3d-vit-large/resolve/main/onnx/model_fp16.onnx",
      "bede661725a7e1808f2cd69d07f3d2a6e5b011a1bde4128e5df9d21232224d8e", 825204259ull},
     {}},
    // The graph's external-data `location` is the bare "model_fp16.onnx_data",
    // so the sibling has to be cached under exactly that name -- renaming it
    // the way the .onnx is renamed would make the reader look for a file that
    // is not there.
    {"metric3d-vit-giant2",
     {"metric3d-vit-giant2-fp16.onnx",
      "https://huggingface.co/onnx-community/metric3d-vit-giant2/resolve/main/onnx/model_fp16.onnx",
      "8ae9cb119397b42dedede351ea648d3cfcb26b6fa7d2eb537d43508c7fe4f72e", 1401125724ull},
     {"model_fp16.onnx_data",
      "https://huggingface.co/onnx-community/metric3d-vit-giant2/resolve/main/onnx/model_fp16.onnx_data",
      "7ebc3f5e95bd919d14b2aa5c1837cbb2c31f9276a5600adf4418b2b78353f4fd", 1355808768ull}},
};

}  // namespace

const ModelSource* find_model_source(const std::string& id) {
    for (const ModelSource& s : kSources)
        if (id == s.id) return &s;
    return nullptr;
}

std::string model_id_list() {
    std::string s;
    for (const ModelSource& m : kSources) {
        if (!s.empty()) s += ", ";
        s += m.id;
    }
    return s;
}

std::string ensure_model(const ModelSource& src) {
    const std::string path = nn::ensure_file(src.onnx, "metric3d");
    if (src.data.file) nn::ensure_file(src.data, "metric3d");
    return path;
}

std::string resolve_model(const std::string& id_or_path) {
    if (const ModelSource* src = find_model_source(id_or_path)) return ensure_model(*src);

    std::error_code ec;
    NN_CHECK(fs::exists(id_or_path, ec),
             "'%s' is neither a known model id (%s) nor a file that exists",
             id_or_path.c_str(), model_id_list().c_str());
    return id_or_path;
}

}  // namespace metric3d
