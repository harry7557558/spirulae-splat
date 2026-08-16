#pragma once

// A handful of frames of one input, and the pixels of any of them: what the
// panels that try a setting on a real picture (SegmentPanel, GeometryPanel)
// both need before a dataset run has extracted anything.
//
// A video is read the way preparation will read it -- in process where the
// driver can, ffmpeg where it cannot or where the job said to -- so what the
// panel shows is what the run will see.

#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

namespace gui {

// One offer on a frame slider. A click records BOTH `index` and `position`,
// because which of the two survives depends on which extraction path the run
// takes (see MaskClick).
struct PreviewFrame {
    std::string path;       // empty for a video: decoded on demand
    long long   index = 0;
    float       position = 0.0f;   // 0..1 through the capture
};

// How an input is read, and how long it is. Built by collect_preview_frames
// and handed back to load_preview_frame, which is why it is one struct: the
// two answers have to agree about which decoder is in use.
struct PreviewSource {
    std::string input;      // video file, or a folder of photos
    bool is_video = false;
    std::string ffmpeg_exe = "ffmpeg";
    bool builtin_decode = true;
    // Seconds, when ffmpeg had to be asked. Zero for photos, and for a video
    // the built-in decoder can seek by frame -- the ffmpeg path is the only
    // one that needs a timestamp.
    double video_seconds = 0.0;
};

// `offers` frames spread through the input. `all_files` comes back holding
// every image of a photo folder (what a border fit reads), empty for a video.
void collect_preview_frames(PreviewSource& src, int offers,
                            std::vector<PreviewFrame>& frames,
                            std::vector<std::string>& all_files,
                            const std::atomic<bool>& cancel);

// One frame's pixels, w*h*3 interleaved. False with `error` set (already a
// sentence for the screen) when it cannot be read.
bool load_preview_frame(const PreviewSource& src, const PreviewFrame& frame,
                        int& w, int& h, std::vector<uint8_t>& rgb,
                        std::string& error, const std::atomic<bool>& cancel);

}  // namespace gui
