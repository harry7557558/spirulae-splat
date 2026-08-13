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
//   frames   in-process VK_KHR_video_decode_* (SS_ENABLE_PATENTED)
//            -> ffmpeg subprocess
//   masks    in-process SAM 2 / SAM 3          (SS_BUILD_SAM)
//            -> python + reference/scripts/mask.py
//
// The built-in path is the default when it is compiled in and the device
// supports it; the fallback is picked automatically otherwise, and can be
// forced (a codec the driver cannot decode, an HDR transfer ffmpeg handles
// better). `Backends` reports what this build and this machine can actually
// do, so the GUI can say so instead of failing at run time.

#include "app/FrameMask.h"

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
    // Masks that came WITH this input -- the `masks/` beside its photos, whose
    // tree mirrors the image tree (see resolve_photo_folder). Empty when it
    // brought none. An input that has its own masks is never AI-masked: the
    // files on disk are the answer, and generating over them would mean
    // writing into a folder the user only asked us to read.
    std::string mask_dir;
    // The lens these images were shot with, and where its focal length starts.
    // Preparation reads neither -- they are here because this list is the one
    // place the inputs are enumerated, and the reconstruction has to be told
    // which folder each setting describes (SfmRunner turns them into
    // `--camera-model DIR=MODEL` / `--focal DIR=PX`). An empty model means the
    // job's dataset-wide one; a factor of 0 means no focal prior.
    std::string camera_model;
    float focal_factor = 0.0f;       // fx = fy = factor * image width
    // Areas of the frame that are never scene -- the fisheye border, a
    // watermark, the rig in shot. Per input because it describes a lens, and
    // resolved per camera folder when it asks for the border to be fitted
    // (app::FrameStencil), which is what gives a dual-fisheye file two circles.
    app::FrameStencil stencil;
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
    // model name reference/scripts/mask.py understands.
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
    // ... and what the trainer's mask_dir should be: "masks" for masks the run
    // put in the dataset, an absolute path for masks it only read (photos used
    // where they are bring theirs with them).
    std::string mask_dir_cfg;
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

// ---- the ffmpeg fallback, for callers that are not a preparation run -------
//
// Everything the GUI does to a video without the built-in decoder goes through
// an external ffmpeg, and the mask preview needs the same two answers a
// preparation run does -- how long the capture is, and what one frame of it
// looks like -- without running one. Only `ffmpeg_exe` is needed (ffprobe is
// not assumed to be installed beside it): `ffmpeg -i` prints the stream table
// on its way to complaining that no output file was named.

// What an external ffmpeg says about a video. A zero means it did not say.
struct VideoFacts {
    double duration = 0.0;    // seconds
    double fps = 0.0;
    long long frames = 0;     // duration * fps; the container's own count is
                              // not printed by `ffmpeg -i`
    int width = 0, height = 0;   // one frame, before any scaling
};
bool ffmpeg_probe_video(const std::string& ffmpeg_exe, const std::string& path,
                        VideoFacts& out, const std::atomic<bool>& cancel);

// One frame, `seconds` into the file, written to `out_path` as a JPEG.
// False when ffmpeg is missing, was cancelled, or wrote nothing.
bool ffmpeg_extract_frame(const std::string& ffmpeg_exe, const std::string& video,
                          double seconds, const std::string& out_path,
                          const std::atomic<bool>& cancel);

// What a picked folder of photos actually means, by the layout conventions the
// rest of the project already uses -- `spirula sfm auto`'s own probing and the
// dataparsers' `mask_dir = "masks"`:
//
//   <picked>/images + <picked>/masks   a dataset folder: index images/, and the
//                                      masks beside it are already made
//   <picked> + <picked>/masks          photos at the top with their masks under
//                                      them (the masks folder is NOT indexed)
//   <picked> + <picked>/../masks       ... or the images folder was picked
//                                      directly and its sibling holds the masks
//
// `images` comes back as the folder to index and `masks` as the mask tree, or
// empty when there is none. Either half may be a symlink into the raw capture,
// which is why every walk in here follows directory symlinks.
void resolve_photo_folder(const std::string& picked, std::string& images,
                          std::string& masks);

// Does this folder hold any image at all, at any depth? Follows directory
// symlinks (a prepared capture's images/ is often a link into the raw one) and
// stops at the first hit, so it is cheap enough for the UI thread.
bool folder_has_images(const std::string& dir);

// Is this folder an already-reconstructed dataset -- something the trainer's
// dataparsers can read -- rather than raw input? True for a Nerfstudio
// transforms.json, a COLMAP sparse/ or colmap/, or a Metashape camera .xml
// next to a point-cloud .ply, which is what the parsers themselves probe for.
//
// Cheap and deliberately shallow: it decides where a dropped folder is routed
// and whether a batch row is worth starting, not whether the parse will
// succeed. It says nothing about raw photos -- a folder of JPEGs is not a
// dataset until something has reconstructed it.
bool folder_looks_like_dataset(const std::string& dir);

// What the output folder already holds, so the screen can say what a run would
// reuse and what it would replace.
//
// `inputs` is part of the question: the natural output folder for a capture
// whose images/ is already on disk is that folder itself, and its own images
// and masks are the input -- not leftovers from a previous run.
struct WorkspaceState {
    bool frames = false;    // images/ this run would extract into
    bool features = false;  // features/, matches.bin, database.db -- reusable
    bool masks = false;     // masks/ this run would generate into
    // sparse/: a reconstruction, which a new run replaces rather than resumes.
    // Deliberately NOT evidence of a resumable run: a folder holding one opens
    // in the trainer when it is dropped, so anything that gets this far with a
    // sparse/ in it is someone else's dataset or a finished one, and the honest
    // thing to do is warn that it is about to be written over.
    bool model = false;
    // Something a resumed run can pick up instead of redoing.
    bool resumable() const { return frames || features || masks; }
};
WorkspaceState probe_workspace(const std::string& workspace,
                               const std::vector<PrepInput>& inputs);

// Is this the mask half of one of those layouts, rather than an input of its
// own? By name, which is what makes it a convention: `--mask-dir masks` is the
// SfM default and `mask_dir = "masks"` the dataparsers'.
bool is_mask_folder(const std::string& path);

class DatasetPrep {
public:
    using LogFn   = std::function<void(const std::string&)>;
    using StageFn = std::function<void(const std::string&)>;

    DatasetPrep(LogFn log, StageFn stage, const std::atomic<bool>& cancel)
        : _log(std::move(log)), _stage(std::move(stage)), _cancel(cancel) {}

    // False with `error` set on failure ("cancelled" when the token was set).
    bool run(const PrepJob& job, PrepResult& out, std::string& error);

    // Recursive, matching what COLMAP's feature_extractor indexes. `skip` is a
    // sub-folder not to descend into: a masks/ nested under the images is full
    // of PNGs that are not views, and counting them doubles the dataset with
    // garbage (the guard `spirula sfm auto` has for the same layout).
    static int count_images(const std::string& dir, const std::string& skip = "");
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
    // read where they are (see run()) -- and the masks they came with into the
    // matching masks/<subdir>, so the two trees still mirror each other.
    bool gather_photos(const PrepJob& job, const PrepInput& in,
                       const std::string& images, const std::string& masks,
                       bool& have_masks, std::string& error);
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
    // The static stencil, intersected with whatever masking already wrote. No
    // model and no GPU, so it runs whether or not segmentation did.
    bool apply_stencil(const PrepInput& in, const std::string& images,
                       const std::string& masks, std::string& error);
    int exec(const std::vector<std::string>& argv);

    LogFn _log;
    StageFn _stage;
    const std::atomic<bool>& _cancel;
};

}  // namespace gui
