#pragma once
// Getting a checkpoint onto disk.
//
// The artifacts are COLMAP's, fetched from COLMAP's release URLs and verified
// against COLMAP's SHA-256. Nothing is
// bundled, nothing is committed, and nothing of ours is hosted anywhere -- so
// there is no converter to keep in step with an upstream re-export, and the
// bytes we run are the bytes the reference implementation runs.
//
// Both models are permissively licensed -- ALIKED is BSD-3-Clause, LightGlue
// Apache-2.0 -- so unlike the segmentation checkpoints
// (src/app/gui/ModelCache.cpp) this needs no consent gate. It still never
// downloads behind the user's back: `spirula sfm` prints what it is fetching
// and from where, and --aliked-model points at a file instead.

#include <cstdint>
#include <string>

namespace aliked {

// A checkpoint we know how to fetch. `sha256` is lowercase hex, and is checked
// after download and on every subsequent load of the cached file -- a
// truncated or tampered artifact must not reach the parser.
struct ModelSource {
    const char* id;      // "aliked-n16rot" -- what --features spells
    const char* file;    // basename in the cache directory
    const char* url;
    const char* sha256;
    uint64_t    bytes;   // approximate, for the "downloading N MB" line
};

// Null when `id` is not one of ours.
const ModelSource* find_model_source(const std::string& id);

// Where a cached checkpoint lives: <cache>/spirula-studio/models/<file>.
// Mirrors src/app/gui/AppPaths.cpp's cache_dir(); duplicated rather than
// shared because src/aliked/ sits below src/app/ in the layering and may not
// include it.
std::string model_cache_path(const ModelSource& src);

// Returns a path to a verified local copy, downloading through the system
// `curl` if needed. Throws nn::Error with an actionable message -- including
// the URL to fetch by hand -- when curl is missing, the download fails, or the
// hash does not match.
std::string ensure_model(const ModelSource& src);

// Resolve what the user asked for: an explicit path is used as-is (and, when
// it happens to be one of the known files, still hash-checked), otherwise the
// id is fetched.
std::string resolve_model(const std::string& id_or_path);

// Lowercase hex SHA-256 of a file's contents. Empty when it cannot be read.
std::string sha256_file(const std::string& path);

}  // namespace aliked
