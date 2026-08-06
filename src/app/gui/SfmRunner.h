#pragma once

// SfmRunner -- turns raw images or a video into a trainable dataset using this
// repository's own structure-from-motion, with nothing else installed.
//
// Same public shape as ColmapRunner (start / cancel / state / stage / error /
// dataset_dir / image_dir / drain_log) plus a progress fraction, so GuiApp
// drives either one through the same code.
//
// Pipeline:
//   DatasetPrep  (frames, sharpest-frame selection, .insv track split, masks)
//   spirula sfm auto  ->  <workspace>/sparse/0/{cameras,images,points3D}.bin
//
// The reconstruction runs as a child process rather than in this one. That is
// a deliberate choice for now, not a shortcut left over from the COLMAP days:
//
//   * the SfM module is still a CLI at heart -- it prints to stdout and has no
//     cancellation token (docs/notes/sfm-port-plan.md phase 3), so in-process
//     it could neither be stopped nor reported on;
//   * global bundle adjustment on a large model and a live trainer must not
//     share a VRAM budget, and a child process gives that separation for free
//     -- every byte it held is gone when it exits;
//   * it keeps one Vulkan device live in the GUI process instead of two
//     (the port plan's own §10 risk).
//
// The child is this same executable run again (AppPaths::exe_path), so there
// is nothing to install and nothing to find: `spirula sfm` is one of the
// subcommands the binary already answers to. When phase 3 lands, only the
// run() body here changes.

#include "app/gui/DatasetPrep.h"
#include "i18n/catalog/Dataset.h"

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

// The camera models `spirula-sfm --camera-model` accepts that the dataset
// parser and renderer can also consume -- which now includes EQUIRECTANGULAR:
// the mapper could always write it, and ColmapParser reads model 17 (it also
// checks the 2:1 aspect a full sphere has to have).
//
// A capture can mix them. `--camera-model PREFIX=MODEL` sets one input's
// model, which is what lets a rig of a 360 camera and a phone reconstruct as
// one scene -- see append_camera_overrides().
inline const char* kSfmCameraModels[] = {
    "opencv", "pinhole", "simple-pinhole", "radial",
    "full-opencv", "opencv-fisheye", "thin-prism-fisheye", "equirectangular",
};
inline constexpr int kNumSfmCameraModels = 8;

// What the combo boxes show for each of the above, in the same order.
inline std::vector<const spirula::i18n::Msg*> sfm_camera_model_labels() {
    namespace m = spirula::i18n::msg::dataset;
    return {&m::lens_opencv, &m::lens_pinhole, &m::lens_simple_pinhole,
            &m::lens_radial, &m::lens_full_opencv, &m::lens_fisheye_kb,
            &m::lens_fisheye_thin_prism, &m::lens_equirectangular};
}

struct SfmJob {
    // ---- shared with ColmapRunner's path ----
    PrepJob prep;

    // ---- reconstruction ----
    int quality = 2;                  // 0 low, 1 medium, 2 high, 3 extreme
    int data_type = 0;                // 0 individual photos, 1 video, 2 internet
    std::string camera_model = "opencv";
    int camera_mode = 1;              // 0 single, 1 per folder, 2 per image
    int pairs = 0;                    // 0 auto, 1 exhaustive, 2 sequential, 3 prefilter
    int overlap = 10;                 // sequential neighbours
    // Sequential pairing only: also match the pairs GPU pair selection finds,
    // so a capture that comes back on itself links across the seam instead of
    // breaking into one model per unbroken run of frames.
    bool loop_closure = true;
    float init_focal_px = 0.0f;       // 0 = guess from EXIF / image size
    int max_features = 0;             // 0 = the quality preset's
    int max_image_size = 0;           // 0 = the quality preset's
    // 0 flat, 1 bottom-up. Flat for every capture, whatever its size: there is
    // no automatic switch, here or in `spirula sfm`.
    int mapper = 0;
    // 0 SIFT, 1 ALIKED-n16rot, 2 ALIKED-n32. The learned ones fetch a
    // checkpoint on first use and run on their own resolution ladder, so the
    // quality preset means something different for each -- which is why this
    // is a frontend choice and not a quality level.
    int features = 0;
    // 0 brute force, 1 LightGlue. Only meaningful with a learned frontend, and
    // an order of magnitude slower per pair -- the panel greys it out for SIFT
    // and the CLI refuses the combination outright.
    int matcher = 0;
    bool keep_intermediate = false;   // keep features/ and matches.bin

    // Extra flags typed by the user, appended verbatim. The escape hatch for
    // everything the panel does not surface -- `spirula-sfm auto --help` lists
    // the lot, and this is how an expert reaches it without us mirroring 130
    // flags into the GUI.
    std::string extra_args;
};

class SfmRunner {
public:
    enum class State { Idle, Running, Done, Failed, Cancelled };

    ~SfmRunner();

    // "" when spirula-sfm is available, otherwise why it is not.
    static std::string availability();

    void start(const SfmJob& job);
    void cancel();

    State state() const { return _state.load(); }
    std::string stage();
    std::string error();
    std::string dataset_dir();
    std::string image_dir();
    // What the trainer's mask_dir should be, or "" when the dataset has no
    // masks. Usually "masks" (the parser default); an absolute path when the
    // masks were only read, which is what photos used where they are do.
    std::string mask_dir();
    // 0..1 within the current stage, or -1 when it cannot be estimated.
    float progress() const { return _progress.load(); }
    // Done, but under half the images registered (or a high reprojection
    // error). The dataset is usable; the user should know it has gaps.
    bool partial() const { return _partial.load(); }
    std::vector<std::string> drain_log();

private:
    void run(SfmJob job);
    void log(const std::string& line);
    void set_stage(const std::string& s);
    // Stage changes driven by the child's output, which repeats a
    // stage's lines many times over.
    void set_stage_if_new(const char* s);
    // Reads one spirula-sfm output line for a progress fraction.
    void note_progress(const std::string& line);
    // Per-input `--camera-model DIR=MODEL` / `--focal DIR=PX`, resolved against
    // the frames that now exist.
    void append_camera_overrides(const SfmJob& job, const PrepResult& prep,
                                 std::vector<std::string>& argv);

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<bool> _cancel{false};
    std::atomic<float> _progress{-1.0f};
    std::atomic<bool> _partial{false};
    std::mutex _mu;
    std::string _stage, _error, _dataset_dir, _image_dir, _mask_dir;
    std::vector<std::string> _log;
};

}  // namespace gui
