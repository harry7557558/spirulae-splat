#pragma once

// ModelCache -- the segmentation checkpoints: which ones exist, whether one is
// already on disk, and fetching it with the user's consent.
//
// Why consent is a first-class thing here and not a download progress bar with
// a licence link under it: the SAM checkpoints are not ours and are not under
// this repository's licence. SAM 2 / 2.1 are Apache-2.0, which is easy. SAM 3
// is under Meta's own SAM 3 licence, which is not an open-source licence and
// is not compatible with GPLv3 -- so we cannot ship it, cannot relicense it,
// and must not fetch it on a user's behalf without them having seen the terms.
// The weights are downloaded at runtime into a cache directory, never bundled,
// and the acceptance is recorded once per family in the settings file.
//
// The goal is one clear sentence and one checkbox, not a wall of legal text: a
// beginner should understand *why* they are being asked, and an expert should
// be able to point --model at their own file and never see this at all.
//
// Downloads run on a worker thread through the system `curl` (Windows 10+ and
// every Linux/macOS ship it), resumable, into "<file>.part" so an interrupted
// fetch never looks like a complete model.

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

// A checkpoint the GUI can offer to fetch.
struct ModelEntry {
    const char* id;          // stable key used in the settings file
    const char* file;        // basename in the cache directory
    const char* label;       // what the combo box shows
    const char* blurb;       // one line under the combo box
    const char* family;      // "sam3" | "sam2" -- the licence unit
    uint64_t    bytes;       // expected download size, for the prompt
    bool        text_prompts;// false = clicks/boxes only (SAM 2 has no text tower)
    // What `scripts/mask.py` calls this checkpoint. Only the Python fallback
    // path uses it -- the built-in masker is handed a file -- but it has to name
    // the SAME model the user picked, or a run that falls back quietly changes
    // which network produced the masks.
    const char* legacy_name;
};

// Ordered best-default-first; index 0 is what a fresh install preselects.
const std::vector<ModelEntry>& model_catalog();
const ModelEntry* find_model(const std::string& id);

// Licence terms for a family, for the consent dialog.
struct LicenseInfo {
    const char* family;
    const char* title;       // "SAM 3 licence (Meta)"
    const char* summary;     // 2-3 short lines, plain language
    const char* url;
    bool        needs_tick;  // Apache-2.0 does not; the SAM 3 licence does
};
const LicenseInfo& license_for(const std::string& family);

// Where a model would live, whether or not it is there yet.
std::string model_path(const ModelEntry& e);
bool model_is_cached(const ModelEntry& e);

// A single background download. One at a time is enough for the GUI, so this
// is a plain object the screen owns rather than a queue.
class ModelDownload {
public:
    enum class State { Idle, Running, Done, Failed, Cancelled };

    ~ModelDownload();

    void start(const ModelEntry& e);
    void cancel();

    State state() const { return _state.load(); }
    // 0..1, or -1 when the server did not send a length.
    float progress() const { return _progress.load(); }
    std::string status();          // "412 MB of 707 MB", or the error
    std::string path();            // valid when Done
    std::vector<std::string> drain_log();

private:
    void run(ModelEntry e);
    void log(const std::string& line);

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<bool>  _cancel{false};
    std::atomic<float> _progress{-1.0f};
    std::mutex _mu;
    std::string _status, _path;
    std::vector<std::string> _log;
};

}  // namespace gui
