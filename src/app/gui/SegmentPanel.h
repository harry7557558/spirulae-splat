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
// something to exclude.
//
// The model lives on the GPU for as long as the panel is open and is dropped
// when it closes -- a reconstruction that follows should not be sharing VRAM
// with a 2 GB backbone that nobody is looking at.

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
};

class SegmentPanel {
public:
    // Both out of line: the pimpl'd Job is incomplete here.
    SegmentPanel();
    ~SegmentPanel();

    // `input` is a folder of images or a video file. Cheap: the frame list is
    // gathered here, nothing is decoded or loaded until the panel draws.
    void open(const std::string& input, bool is_video,
              const std::string& model_path);
    bool is_open() const { return _open; }
    void close();

    // Draws the modal window. Call once per frame from the dataset screen,
    // inside its ImGui frame. `settings` is edited in place.
    void draw(MaskSettings& settings);

    // Frees the GL texture; call while the GL context is current.
    void destroy_gl();

private:
    struct Job;

    void start_job(const MaskSettings& s);
    void collect_frames(const std::string& input, bool is_video);
    void upload_preview();
    void draw_image(MaskSettings& settings);

    bool _open = false;
    std::string _model_path;
    std::string _input;
    bool _is_video = false;

    std::vector<std::string> _frames;   // image paths, or "" for video frames
    int  _frame_idx = 0;
    bool _frame_dirty = true;           // the chosen frame changed
    bool _needs_run = false;            // prompt edited; rerun on release

    // Clicks, in pixels of the source frame.
    struct Click { float x, y; bool positive; };
    std::vector<Click> _clicks;

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
