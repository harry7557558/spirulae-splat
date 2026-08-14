#pragma once

// SegmentPanel -- try a mask prompt on one real frame before committing to a
// run over the whole capture.
//
// Masking is the one dataset setting a beginner cannot reason about in the
// abstract: "people; cars" either catches the thing walking through the shot
// or it does not, and finding out after a twenty-minute reconstruction is the
// wrong time. This shows the actual mask, from the actual model, on a frame of
// the actual input, in about a second per attempt once the checkpoint is
// loaded -- and it is the same sam::Masker the dataset run uses, so what is on
// screen is what will be written.
//
// Clicks are supported too, which is the only way to prompt a SAM 2 checkpoint
// (it has no text tower): left-click marks the subject, right-click marks
// something to exclude. They belong to an object and to the frame they were
// drawn on -- see MaskClick -- and they are kept in the settings rather than
// in the panel, because they are prompts for the run and not a preview toy.
//
// The model lives on the GPU for as long as the panel is open and is dropped
// when it closes -- a reconstruction that follows should not be sharing VRAM
// with a 2 GB backbone that nobody is looking at.
//
// The same panel edits the input's static stencil (app::FrameStencil): the
// fisheye border it finds by itself, plus circles and boxes dragged over the
// picture. That half needs no model at all, and the panel stays usable -- and
// says so -- when there is none.

#include "app/FrameMask.h"
#include "app/gui/DatasetPrep.h"   // MaskClick
#include "app/gui/GlLoader.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

// The mask settings the panel edits, owned by the dataset screen so what is
// tried here is what runs.
struct MaskSettings {
    std::string prompt;              // "people; cars"
    std::string negative_prompt;
    bool keep_subject = false;       // prompt names what to KEEP
    int  max_image_size = 1600;
    float threshold = 0.5f;
    // Clicked objects, across every frame and every input the user visited;
    // each carries the input it was drawn on (MaskClick::source).
    std::vector<MaskClick> clicks;
    int object_count = 1;            // how many the user has opened
    int current_object = 0;          // which one a new click joins
};

class SegmentPanel {
public:
    // Both out of line: the pimpl'd Job is incomplete here.
    SegmentPanel();
    ~SegmentPanel();

    // `input` is a folder of images or a video file. Cheap: the frame list is
    // gathered here, nothing is decoded or loaded until the panel draws (a
    // video is asked how long it is, which is one probe either way).
    //
    // `ffmpeg_exe` and `force_ffmpeg` are the same two settings the dataset run
    // takes: a machine whose driver cannot decode video reads the preview frame
    // with an external ffmpeg, exactly as preparation would.
    void open(const std::string& input, bool is_video,
              const std::string& model_path, const std::string& ffmpeg_exe,
              bool force_ffmpeg);
    bool is_open() const { return _open; }
    void close();

    // Draws the modal window. Call once per frame from the dataset screen,
    // inside its ImGui frame. Both arguments are edited in place.
    void draw(MaskSettings& settings, app::FrameStencil& stencil);

    // Frees the GL textures; call while the GL context is current.
    void destroy_gl();

private:
    struct Rgb;
    struct Job;

    void start_job(const MaskSettings& s, const app::FrameMask& stencil);
    void start_detect();
    void collect_frames(const std::string& input, bool is_video);
    void upload_preview();
    void upload_stencil(const app::FrameMask& stencil);
    void draw_image(MaskSettings& settings, app::FrameStencil& stencil,
                    bool& edited);
    void draw_objects(MaskSettings& settings, bool& edited);
    void draw_stencil(app::FrameStencil& stencil, bool& edited);

    // The stencil as it will be written: the shapes plus, when asked for, the
    // border this panel found, shrunk by the current amount.
    app::FrameMask resolved(const app::FrameStencil& stencil) const;

    // Is the built-in decoder both compiled in and usable on this device, and
    // not switched off by the job's "always use ffmpeg"? False means every
    // frame this panel shows comes from ffmpeg.
    bool builtin_decode() const;

    // The settings hold every input's clicks; these are the ones on this picture.
    bool mine(const MaskClick& c) const { return c.source == _input; }

    bool _open = false;
    std::string _model_path;
    std::string _input;
    bool _is_video = false;
    std::string _ffmpeg_exe = "ffmpeg";
    bool _force_ffmpeg = false;
    // How long the video is, in seconds, when ffmpeg had to be asked. Zero for
    // a folder of photos, and for a video the built-in decoder can seek by
    // frame -- the ffmpeg path is the only one that needs a timestamp.
    double _video_seconds = 0.0;

    // One offer on the frame slider. `index` is what the masking run will call
    // this frame and `position` is where it falls in the capture; a click
    // records both, because which of the two survives depends on which
    // extraction path runs (see MaskClick).
    struct Frame {
        std::string path;       // empty for a video: decoded on demand
        long long   index = 0;
        float       position = 0.0f;
    };
    std::vector<Frame> _frames;
    // Every image of a photo input, which is what the border fit reads. The
    // slider offers a dozen of them; a fit wants a spread of two dozen.
    std::vector<std::string> _all_files;
    int  _frame_idx = 0;
    bool _frame_dirty = true;           // the chosen frame changed
    bool _needs_run = false;            // prompt edited; rerun on release

    // ---- the stencil ----
    // The border found for the frames this panel is showing, with no shrink
    // applied -- the slider re-applies it without another fit. The dataset run
    // fits one per camera; this is the one the preview can show.
    app::BorderDetect _border;          // UI thread
    app::BorderDetect _border_pending;  // guarded by _mu
    bool _border_ready = false;         // guarded by _mu
    std::atomic<bool> _detecting{false};
    bool _detect_asked = false;         // a fit has been asked for at least once
    int  _shape_sel = -1;               // which shape carries handles
    // Which of its handles is being dragged, kDragBody for the whole shape,
    // -1 for none. The body drag is relative, so it also remembers where the
    // pointer was last frame.
    static constexpr int kDragBody = -2;
    int   _drag_handle = -1;
    float _drag_from_u = 0.0f, _drag_from_v = 0.0f;
    GLuint _stencil_tex = 0;
    std::string _stencil_key;           // what _stencil_tex was built from

    // The composited RGB preview handed to GL, guarded by _mu.
    std::mutex _mu;
    std::vector<uint8_t> _preview;
    int _preview_w = 0, _preview_h = 0;
    bool _preview_dirty = false;
    std::string _status, _error;
    float _kept_fraction = -1.0f;

    GLuint _tex = 0;
    int _tex_w = 0, _tex_h = 0;

    std::thread _worker;
    std::atomic<bool> _busy{false};
    std::atomic<bool> _cancel{false};
    std::unique_ptr<Job> _job;          // the session, kept warm between runs
};

}  // namespace gui
