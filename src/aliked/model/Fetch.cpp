#include "aliked/model/Fetch.h"

#include "core/Sha256.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

namespace fs = std::filesystem;

namespace aliked {
namespace {

// The same URL / filename / SHA-256 triples COLMAP carries in
// src/colmap/feature/resources.h. Keep them identical: the point of fetching
// from COLMAP's release is that a parity run compares implementations, not
// checkpoints.
const ModelSource kSources[] = {
    {"aliked-n16rot", "aliked-n16rot.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n16rot.onnx",
     "39c423d0a6f03d39ec89d3d1d61853765c2fb6a8b8381376c703e5758778a547", 2997054ull},
    {"aliked-n32", "aliked-n32.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n32.onnx",
     "a077728a02d2de1a775c66df6de8cfeb7c6b51ca57572c64c680131c988c8b3c", 4205634ull},
    {"aliked-lightglue", "aliked-lightglue.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-lightglue.onnx",
     "b9a5de7204648b18a8cf5dcac819f9d30de1a5961ef03756803c8b86c2dceb8d", 0ull},
};

std::string env_str(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::string(v) : std::string();
}

fs::path cache_root() {
#ifdef _WIN32
    fs::path dir = env_str("LOCALAPPDATA");
#else
    fs::path dir = env_str("XDG_CACHE_HOME");
    if (dir.empty()) {
        const std::string home = env_str("HOME");
        if (!home.empty()) dir = fs::path(home) / ".cache";
    }
#endif
    if (dir.empty()) dir = ".";
    // An existing spirulae-splat/ from before the rename is adopted where it
    // is -- the checkpoints in it are a big download to repeat.
    std::error_code ec;
    if (!fs::exists(dir / "spirula-studio", ec) &&
        fs::is_directory(dir / "spirulae-splat", ec))
        return dir / "spirulae-splat";
    return dir / "spirula-studio";
}

bool have_curl() {
#ifdef _WIN32
    return std::system("curl --version >NUL 2>&1") == 0;
#else
    return std::system("curl --version >/dev/null 2>&1") == 0;
#endif
}

}  // namespace

const ModelSource* find_model_source(const std::string& id) {
    for (const ModelSource& s : kSources)
        if (id == s.id) return &s;
    return nullptr;
}

std::string model_cache_path(const ModelSource& src) {
    return (cache_root() / "models" / src.file).string();
}

// The checkpoints are verified on every load, not only after download: a
// truncated or tampered artifact must not reach the parser.
std::string sha256_file(const std::string& path) {
    return spirula::sha256_file(path);
}

std::string ensure_model(const ModelSource& src) {
    const fs::path dst = model_cache_path(src);
    std::error_code ec;

    if (fs::exists(dst, ec)) {
        const std::string got = sha256_file(dst.string());
        if (got == src.sha256) return dst.string();
        // A cached file that does not hash is a failed or interrupted download
        // from a previous run, not a reason to stop: say so and refetch.
        NN_LOG_WARN("[aliked] cached %s does not match its checksum; re-downloading\n",
                    src.file);
        fs::remove(dst, ec);
    }

    fs::create_directories(dst.parent_path(), ec);
    NN_CHECK(!ec, "cannot create %s: %s", dst.parent_path().string().c_str(),
             ec.message().c_str());

    NN_CHECK(have_curl(),
             "curl was not found, and it is how checkpoints are fetched.\n"
             "  Install curl, or download\n    %s\n  to\n    %s\n  by hand.",
             src.url, dst.string().c_str());

    fs::path part = dst;
    part += ".part";

    NN_LOG_INFO("[aliked] fetching %s (%.1f MB) from %s\n", src.file,
                (double)src.bytes / 1e6, src.url);
    // -C - resumes a partial .part file; -f makes an HTTP error an exit code
    // rather than a saved error page. Downloading into .part and renaming is
    // what keeps an interrupted fetch from ever looking like a complete model.
    std::string cmd = "curl -L -f --progress-bar -C - -o \"" + part.string() + "\" \"" +
                      std::string(src.url) + "\"";
    const int rc = std::system(cmd.c_str());
    if (rc != 0) {
        fs::remove(part, ec);
        nn::fail("downloading %s failed (curl exit %d).\n"
                 "  Fetch it by hand from\n    %s\n  and save it as\n    %s",
                 src.file, rc, src.url, dst.string().c_str());
    }

    const std::string got = sha256_file(part.string());
    if (got != src.sha256) {
        fs::remove(part, ec);
        nn::fail("%s downloaded but its SHA-256 is\n    %s\n  expected\n    %s\n"
                 "  The file was discarded.",
                 src.file, got.c_str(), src.sha256);
    }

    fs::rename(part, dst, ec);
    NN_CHECK(!ec, "cannot move the download into place: %s", ec.message().c_str());
    NN_LOG_INFO("[aliked] saved %s\n", dst.string().c_str());
    return dst.string();
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
