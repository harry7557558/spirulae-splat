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

// One thing the user picked: a video file, or a folder of photos. A job holds
// a list of them, because a capture is often shot as several clips, or on a rig
// whose lenses each write their own file -- and those only reconstruct together
// if they end up in one image tree.
struct PrepInput {
    std::string path;                // video file, or a folder of photos
    bool is_video = false;
    // Where this input's images land, relative to the dataset's images/. Empty
    // means images/ itself, which is what a lone input gets so a one-video
    // dataset keeps the layout it has always had. A multi-track video (an
    // Insta360 .insv carries two fisheye streams) adds cam0/, cam1/ ... below.
    std::string subdir;
    // The lens these images were shot with, and where its focal length starts.
    // Preparation reads neither -- they are here because this list is the one
    // place the inputs are enumerated, and the reconstruction has to be told
    // which folder each setting describes (SfmRunner turns them into
    // `--camera-model DIR=MODEL` / `--focal DIR=PX`). An empty model means the
    // job's dataset-wide one; a factor of 0 means no focal prior.
    std::string camera_model;
    float focal_factor = 0.0f;       // fx = fy = factor * image width
};

struct PrepJob {
    std::vector<PrepInput> inputs;   // in the order the user added them
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
    // Which input they were drawn on ("" = the only one). A click is a point on
    // a frame of one capture; carrying it into a second video would prompt the
    // model with whatever happens to be at those coordinates there, which is
    // the bug that looks like a working feature. Inputs it does not name are
    // masked from the text prompt alone, or skipped when there is none.
    std::string mask_clicks_source;

    // Built-in: a checkpoint file (ModelCache resolves it). External: the
    // model name scripts/mask.py understands.
    std::string mask_model_path;
    std::string mask_model_name = "sam2.1_hiera_large";
    bool  force_external_masking = false;
    std::string python_exe = "python3";
};

// The one case where images are read where they are instead of being gathered
// into the dataset's own images/ (see DatasetPrep::run): one folder of photos.
inline bool reads_photos_in_place(const std::vector<PrepInput>& inputs) {
    return inputs.size() == 1 && !inputs[0].is_video;
}

struct PrepResult {
    std::string image_dir;           // absolute; what SfM should index
    std::string image_dir_cfg;       // what the trainer's image_dir should be
    std::string mask_dir;            // "" when there are no masks
    int  n_images = 0;
    // images/ came out holding one sub-folder per camera -- several inputs, or
    // a multi-track video -- so intrinsics must not be shared across them.
    bool per_folder_cameras = false;
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

// Does this path name one of them? (Extension only; the file need not exist.)
bool is_video_path(const std::string& path);
// A dual-fisheye Insta360 file: two video tracks, one per lens, and a lens the
// default camera model does not fit.
bool is_dual_fisheye_path(const std::string& path);

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
    // One input's frames. `images` / `masks` are that input's own folders
    // (images/<subdir>, masks/<subdir>); `masked` comes back true when the
    // frames were masked on the way out of the decoder, which is the pass that
    // never has to read them back.
    bool extract_video(const PrepJob& job, const PrepInput& in,
                       const std::string& images, const std::string& masks,
                       PrepResult& out, bool& masked, std::string& error);
    bool extract_video_builtin(const PrepJob& job, const PrepInput& in,
                               const std::string& images, const std::string& masks,
                               PrepResult& out, bool& masked, std::string& error);
    bool extract_video_ffmpeg(const PrepJob& job, const PrepInput& in,
                              const std::string& images, PrepResult& out,
                              std::string& error);
    // Photos into the dataset's own images/<subdir>, when they cannot simply be
    // read where they are (see run()).
    bool gather_photos(const PrepJob& job, const PrepInput& in,
                       const std::string& images, std::string& error);
    // Masks for ONE input's images. Run per input rather than over the whole
    // tree so the tracker's memory bank never crosses from one capture into the
    // next, and so clicks reach only the input they were drawn on.
    bool generate_masks(const PrepJob& job, const PrepInput& in,
                        const std::string& images, const std::string& images_rel,
                        const std::string& masks, const std::string& masks_rel,
                        std::string& error);
    bool generate_masks_builtin(const PrepJob& job, const PrepInput& in,
                                const std::string& images, const std::string& masks,
                                std::string& error);
    bool generate_masks_python(const PrepJob& job, const std::string& images_rel,
                               const std::string& masks_rel, std::string& error);
    int exec(const std::vector<std::string>& argv);

    LogFn _log;
    StageFn _stage;
    const std::atomic<bool>& _cancel;
};

}  // namespace gui
