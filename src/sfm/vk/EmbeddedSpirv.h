// Lookup for the SPIR-V blobs compiled into the binary.
//
// The definitions live in the generated TU that cmake/SsplatSfm.cmake produces
// (sfm_shaders_embedded.cpp) and the `ssplat_sfm` library target compiles;
// blob names are the .spv file stems, e.g. "ba_double_huber". Which variants
// exist depends on SSPLAT_SFM_REALS / SSPLAT_SFM_LOSSES at configure time.
#pragma once

#include <cstddef>
#include <cstdint>

namespace sfm {

// Returns nullptr if the named variant was not built into this binary.
// *words (optional) receives the blob length in 32-bit words.
const uint32_t* findSpirv(const char* name, size_t* words);

// Enumerate blob names; returns nullptr past the end.
const char* spirvBlobName(size_t i);

}  // namespace sfm
