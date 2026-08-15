#include "nn/io/Fetch.h"

#include "core/Sha256.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace nn {
namespace {

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

std::string model_cache_dir() {
    return (cache_root() / "models").string();
}

std::string cached_path(const FetchFile& f) {
    return (cache_root() / "models" / f.file).string();
}

std::string sha256_file(const std::string& path) {
    return spirula::sha256_file(path);
}

std::string ensure_file(const FetchFile& f, const char* tag) {
    const fs::path dst = cached_path(f);
    std::error_code ec;

    if (fs::exists(dst, ec)) {
        const std::string got = sha256_file(dst.string());
        if (got == f.sha256) return dst.string();
        // A cached file that does not hash is a failed or interrupted download
        // from a previous run, not a reason to stop: say so and refetch.
        NN_LOG_WARN("[%s] cached %s does not match its checksum; re-downloading\n", tag,
                    f.file);
        fs::remove(dst, ec);
    }

    fs::create_directories(dst.parent_path(), ec);
    NN_CHECK(!ec, "cannot create %s: %s", dst.parent_path().string().c_str(),
             ec.message().c_str());

    NN_CHECK(have_curl(),
             "curl was not found, and it is how checkpoints are fetched.\n"
             "  Install curl, or download\n    %s\n  to\n    %s\n  by hand.",
             f.url, dst.string().c_str());

    fs::path part = dst;
    part += ".part";

    NN_LOG_INFO("[%s] fetching %s (%.1f MB) from %s\n", tag, f.file,
                (double)f.bytes / 1e6, f.url);
    // -C - resumes a partial .part file; -f makes an HTTP error an exit code
    // rather than a saved error page. Downloading into .part and renaming is
    // what keeps an interrupted fetch from ever looking like a complete model.
    std::string cmd = "curl -L -f --progress-bar -C - -o \"" + part.string() + "\" \"" +
                      std::string(f.url) + "\"";
    const int rc = std::system(cmd.c_str());
    if (rc != 0) {
        fs::remove(part, ec);
        nn::fail("downloading %s failed (curl exit %d).\n"
                 "  Fetch it by hand from\n    %s\n  and save it as\n    %s",
                 f.file, rc, f.url, dst.string().c_str());
    }

    const std::string got = sha256_file(part.string());
    if (got != f.sha256) {
        fs::remove(part, ec);
        nn::fail("%s downloaded but its SHA-256 is\n    %s\n  expected\n    %s\n"
                 "  The file was discarded.",
                 f.file, got.c_str(), f.sha256);
    }

    fs::rename(part, dst, ec);
    NN_CHECK(!ec, "cannot move the download into place: %s", ec.message().c_str());
    NN_LOG_INFO("[%s] saved %s\n", tag, dst.string().c_str());
    return dst.string();
}

}  // namespace nn
