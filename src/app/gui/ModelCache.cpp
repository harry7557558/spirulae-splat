// ModelCache.cpp -- see ModelCache.h.

#include "app/gui/ModelCache.h"

#include "app/gui/AppPaths.h"
#include "app/gui/Subprocess.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

// The sam3.cpp conversions, which are byte-compatible with what this tree
// loads. Mirroring Meta's originals rather than converting them ourselves is
// deliberate: one published artifact, one checksum, one place to look when a
// checkpoint misbehaves.
const char* kBaseUrl = "https://huggingface.co/PABannier/sam3.cpp/resolve/main/";

std::string human_bytes(uint64_t b) {
    char buf[32];
    if (b >= (1ull << 30)) std::snprintf(buf, sizeof buf, "%.1f GB", b / 1073741824.0);
    else                   std::snprintf(buf, sizeof buf, "%.0f MB", b / 1048576.0);
    return buf;
}

}  // namespace

const std::vector<ModelEntry>& model_catalog() {
    // Quantized SAM 3 first: it is the only family that understands "mask out
    // the people", it is a third the size of the f16 weights, and the
    // quantization is of the *file* -- everything is dequantized to fp16 on
    // upload, so the accuracy difference is small and the download is not.
    //
    // The four SAM 2.1 sizes are all here because the choice is a real one and
    // nobody can make it from a name: they differ by ~2.5x in speed and the two
    // in the middle are where most captures want to be. The blurbs quote the
    // same measurement for each (one instance, 1080p frames, laptop GPU) so
    // they can actually be compared -- see src/sam/README.md for the table.
    static const std::vector<ModelEntry> kCatalog = {
        {"sam3-q4_0", "sam3-q4_0.ggml", "SAM 3 (recommended)",
         "Understands text prompts -- type what to mask out. 707 MB, ~2 GB "
         "VRAM, about 1 s per frame on a laptop GPU -- 3x slower than any "
         "SAM 2.1 below.",
         "sam3", 707ull << 20, true, "sam3"},
        {"sam3-f16", "sam3-f16.ggml", "SAM 3, full precision",
         "The same model without file quantization. Slightly better masks, "
         "much bigger download, same speed.",
         "sam3", 1884ull << 20, true, "sam3"},
        {"sam2.1-large", "sam2.1_hiera_large_f16.ggml", "SAM 2.1 Large",
         "Click or draw a box to select an object; no text prompts. The most "
         "accurate of the four and the one to pick for thin structure -- "
         "hair, railings, foliage. ~470 ms per frame. Apache-2.0.",
         "sam2", 430ull << 20, false, "sam2.1_hiera_large"},
        {"sam2.1-base-plus", "sam2.1_hiera_base_plus_f16.ggml", "SAM 2.1 Base+",
         "Clicks and boxes only. Close to Large on most subjects at two "
         "thirds the time, ~320 ms per frame. Apache-2.0.",
         "sam2", 156ull << 20, false, "sam2.1_hiera_base_plus"},
        {"sam2.1-small", "sam2.1_hiera_small_f16.ggml", "SAM 2.1 Small",
         "Clicks and boxes only, ~255 ms per frame. The best speed-for-quality "
         "of the four: below Large the frame is mostly tracking, which does "
         "not care how big the backbone is. Apache-2.0.",
         "sam2", 89ull << 20, false, "sam2.1_hiera_small"},
        {"sam2.1-tiny", "sam2.1_hiera_tiny_f16.ggml", "SAM 2.1 Tiny (fastest)",
         "Clicks and boxes only, 76 MB. Only 4% quicker than Small at ~245 ms "
         "per frame, and it loses thin structure first -- take it for the "
         "download size, not the speed. Apache-2.0.",
         "sam2", 76ull << 20, false, "sam2.1_hiera_tiny"},
    };
    return kCatalog;
}

const ModelEntry* find_model(const std::string& id) {
    for (const auto& e : model_catalog())
        if (id == e.id) return &e;
    return nullptr;
}

const LicenseInfo& license_for(const std::string& family) {
    // Written for someone who has not read a licence before. What they need to
    // know is (a) it is not ours, (b) whether they are agreeing to anything
    // beyond the ordinary, and (c) where the actual text is.
    static const LicenseInfo kSam3{
        "sam3", "SAM 3 License (Meta)",
        "SAM 3 is Meta's model, not part of Spirula Studio, and it comes with "
        "its own licence -- which is not a standard one. It is free to "
        "use, including commercially, but only on Meta's terms, so we cannot "
        "ship it with the app or accept them for you.\n\n"
        "Please read it before continuing -- it is short, and it is the actual "
        "agreement, not this summary of it.",
        "https://github.com/facebookresearch/sam3/blob/main/LICENSE", true};
    static const LicenseInfo kSam2{
        "sam2", "SAM 2.1 License (Meta, Apache-2.0)",
        "SAM 2.1 is Meta's model, released under the Apache 2.0 licence. "
        "Nothing unusual to agree to; it is downloaded rather than bundled "
        "only to keep the app small.",
        "https://github.com/facebookresearch/sam2/blob/main/LICENSE", false};
    return family == "sam2" ? kSam2 : kSam3;
}

std::string model_path(const ModelEntry& e) {
    return (fs::path(cache_dir()) / "models" / e.file).string();
}

bool model_is_cached(const ModelEntry& e) {
    std::error_code ec;
    const fs::path p = model_path(e);
    if (!fs::is_regular_file(p, ec)) return false;
    // A truncated file that somehow escaped the .part rename would fail deep
    // inside the weight loader with an unhelpful message; catch it here. The
    // catalog sizes are approximate, so this is a floor, not an equality.
    return fs::file_size(p, ec) > e.bytes / 2;
}

// ---------------------------------------------------------------------------
// Download
// ---------------------------------------------------------------------------

ModelDownload::~ModelDownload() {
    cancel();
    if (_worker.joinable()) _worker.join();
}

void ModelDownload::start(const ModelEntry& e) {
    if (_state.load() == State::Running) return;
    if (_worker.joinable()) _worker.join();
    _cancel = false;
    _progress = -1.0f;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _status.clear();
        _path.clear();
    }
    _state = State::Running;
    _worker = std::thread([this, e] { run(e); });
}

void ModelDownload::cancel() { _cancel = true; }

std::string ModelDownload::status() {
    std::lock_guard<std::mutex> lk(_mu);
    return _status;
}

std::string ModelDownload::path() {
    std::lock_guard<std::mutex> lk(_mu);
    return _path;
}

std::vector<std::string> ModelDownload::drain_log() {
    std::lock_guard<std::mutex> lk(_mu);
    std::vector<std::string> out;
    out.swap(_log);
    return out;
}

void ModelDownload::log(const std::string& line) {
    std::lock_guard<std::mutex> lk(_mu);
    _log.push_back(line);
    if (_log.size() > 500) _log.erase(_log.begin(), _log.begin() + 200);
}

void ModelDownload::run(ModelEntry e) {
    auto fail = [&](const std::string& why) {
        std::lock_guard<std::mutex> lk(_mu);
        _status = why;
        _state = _cancel.load() ? State::Cancelled : State::Failed;
    };

    const fs::path dst = model_path(e);
    std::error_code ec;
    fs::create_directories(dst.parent_path(), ec);
    if (ec) return fail("cannot create " + dst.parent_path().string());

    if (!command_exists("curl"))
        return fail("curl was not found. Install curl, or download\n" +
                    std::string(kBaseUrl) + e.file + "\nto " + dst.string() +
                    " by hand.");

    fs::path part = dst;
    part += ".part";

    log("Downloading " + std::string(e.file) + " (" + human_bytes(e.bytes) + ")");
    // -C - resumes a partial .part file, so a cancelled or dropped download
    // does not start over. --progress-bar writes "###...  42.0%" to stderr;
    // that percentage is the only progress the GUI needs.
    const int rc = run_process(
        {"curl", "-L", "-f", "--progress-bar", "-C", "-", "-o", part.string(),
         std::string(kBaseUrl) + e.file},
        "",
        [&](const std::string& line) {
            size_t pct = line.find('%');
            if (pct != std::string::npos) {
                size_t s = pct;
                while (s > 0 && (std::isdigit((unsigned char)line[s - 1]) ||
                                 line[s - 1] == '.'))
                    s--;
                if (s < pct) {
                    const float v = std::strtof(line.substr(s, pct - s).c_str(), nullptr);
                    _progress = v / 100.0f;
                    char buf[64];
                    std::snprintf(buf, sizeof buf, "%.0f%% of %s", v,
                                  human_bytes(e.bytes).c_str());
                    std::lock_guard<std::mutex> lk(_mu);
                    _status = buf;
                }
                return;
            }
            log(line);
        },
        _cancel);

    if (rc == kCancelled) {
        // The .part file stays: -C - picks it up if the user tries again.
        return fail("cancelled");
    }
    if (rc != 0) {
        fs::remove(part, ec);
        return fail("download failed (curl exit " + std::to_string(rc) +
                    "); see the log");
    }

    fs::rename(part, dst, ec);
    if (ec) return fail("cannot move the download into place: " + ec.message());

    _progress = 1.0f;
    {
        std::lock_guard<std::mutex> lk(_mu);
        _path = dst.string();
        _status = "ready";
    }
    log("Saved to " + dst.string());
    _state = State::Done;
}

}  // namespace gui
