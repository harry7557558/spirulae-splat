#pragma once

// DatasetPrep -- everything between "the user picked a video or a folder of
// photos" and "there is an image directory (and maybe a mask directory) ready
// for structure-from-motion".
//
// It was the half of ColmapRunner that had nothing to do with COLMAP: frame
// extraction, sharpest-frame selection, splitting a multi-track 360 file into
// one folder per lens, AI masking, and the resume rules that let an
// interrupted run reuse what it left behind. Both dataset paths -- the
// built-in SfM and the external COLMAP -- run it first, so a change to how
// frames are chosen cannot apply to only one of them.
//
// Each stage has a built-in implementation and an external fallback:
//
//   frames   in-process VK_KHR_video_decode_* (SSPLAT_ENABLE_PATENTED)
//            -> ffmpeg subprocess
//   masks    in-process SAM 2 / SAM 3          (SSPLAT_BUILD_SAM)
//            -> python + scripts/mask.py
//
// The built-in path is the default when it is compiled in and the device
// supports it; the fallback is picked automatically otherwise, and can be
// forced (a codec the driver cannot decode, an HDR transfer ffmpeg handles
// better). `Backends` reports what this build and this machine can actually
// do, so the GUI can say so instead of failing at run time.

#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace gui {

// A click made in the mask preview, kept because it is a prompt for the whole
// run and not just for the frame it was drawn on.
//
// Objects are separate on purpose: SAM segments ONE thing per prompt, and a
// single instance given a click on the dog and a click on the bicycle returns
// a mask that fits neither. Each object is tracked on its own and the masks
// are unioned at the end.
//
// The frame is recorded twice because the two masking paths number frames
// differently. `frame` is exact where it can be -- the decoded index for the
// path that reads the video itself, the position in the sorted list for a
// folder of photos -- while `position`, the same point as a fraction of the
// capture, is what the ffmpeg path has to fall back on, since it resamples the
// video to a frame rate the preview never saw.
struct MaskClick {
    float x = 0.0f, y = 0.0f;   // pixels of the source frame
    bool  positive = true;      // "this is it" vs "not this"
    int   object = 0;
    long long frame = 0;
    float position = 0.0f;      // 0..1 through the capture
};

struct PrepJob {
    std::string input_path;          // photo folder, or a video file
    bool is_video = false;
    std::string workspace;           // output dataset dir (created)
    bool resume = true;              // reuse what a previous run completed

    // ---- video extraction ----
    float video_fps = 2.0f;          // kept frames per second
    int   sharp_window = 3;          // keep the sharpest of N (1 = off)
    int   max_frames = 100000;
    bool  force_external_decode = false;
    std::string ffmpeg_exe = "ffmpeg";

    // ---- masking ----
    bool mask_enable = false;
    std::string mask_prompt;         // "people; cars; ..."
    std::string mask_negative_prompt;
    bool mask_keep_subject = false;  // prompt names what to KEEP, not remove
    int  mask_max_image_size = 1600;
    // Clicked objects. The only way to prompt a SAM 2 checkpoint, and usable
    // alongside a text prompt on a SAM 3 one.
    std::vector<MaskClick> mask_clicks;

    // Built-in: a checkpoint file (ModelCache resolves it). External: the
    // model name scripts/mask.py understands.
    std::string mask_model_path;
    std::string mask_model_name = "sam2.1_hiera_large";
    bool  force_external_masking = false;
    std::string python_exe = "python3";
};

struct PrepResult {
    std::string image_dir;           // absolute; what SfM should index
    std::string image_dir_cfg;       // what the trainer's image_dir should be
    std::string mask_dir;            // "" when there are no masks
    int  n_images = 0;
    bool multi_track = false;        // one camera per folder is required
};

// What this build, on this machine, can do without an external tool.
//
// Two strings per stage on purpose: `*_reason` names the option or the missing
// device feature and belongs in a log or a tooltip; `*_note` is the sentence
// shown on the screen, which should tell a user what will happen rather than
// which CMake flag was off when the binary was made.
struct Backends {
    bool builtin_video = false;
    std::string video_reason;
    std::string video_note;
    bool builtin_masking = false;
    std::string masking_reason;
    std::string masking_note;
};
// Probes the Vulkan device on first call for the video answer, so call it off
// the UI thread the first time if that matters.
const Backends& backends();

// Video container extensions the GUI offers, in the file dialog and for
// drag-and-drop. Sized here so a range-for over it works from another TU.
inline constexpr int kNumVideoExtensions = 12;
extern const char* const kVideoExtensions[kNumVideoExtensions];

class DatasetPrep {
public:
    using LogFn   = std::function<void(const std::string&)>;
    using StageFn = std::function<void(const std::string&)>;

    DatasetPrep(LogFn log, StageFn stage, const std::atomic<bool>& cancel)
        : _log(std::move(log)), _stage(std::move(stage)), _cancel(cancel) {}

    // False with `error` set on failure ("cancelled" when the token was set).
    bool run(const PrepJob& job, PrepResult& out, std::string& error);

    // Recursive, matching what COLMAP's feature_extractor indexes.
    static int count_images(const std::string& dir);
    // Dimensions of the first image found, for the focal-length prior.
    static bool first_image_dims(const std::string& dir, int& w, int& h);

private:
    bool extract_video(const PrepJob& job, PrepResult& out, std::string& error);
    bool extract_video_builtin(const PrepJob& job, PrepResult& out,
                               std::string& error);
    bool extract_video_ffmpeg(const PrepJob& job, PrepResult& out,
                              std::string& error);
    bool generate_masks(const PrepJob& job, const std::string& images,
                        const std::string& images_rel, std::string& error);
    bool generate_masks_builtin(const PrepJob& job, const std::string& images,
                                std::string& error);
    bool generate_masks_python(const PrepJob& job, const std::string& images_rel,
                               std::string& error);
    int exec(const std::vector<std::string>& argv);

    LogFn _log;
    StageFn _stage;
    const std::atomic<bool>& _cancel;
};

}  // namespace gui
