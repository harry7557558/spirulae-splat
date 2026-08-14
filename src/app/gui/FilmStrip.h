#pragma once

// FilmStrip -- the last few frames a dataset run has written, on the screen
// while it writes the next one.
//
// The point is that masking is the one stage whose output a user cannot judge
// from a counter: "1400 images masked" says nothing about whether the prompt
// caught the thing walking through the shot. This shows the frame with its
// mask composited exactly as the preview panel does.
//
// Lossy on purpose: a decode running at hundreds of frames a second must never
// wait for the screen, so `offer` declines rather than blocks, and the ring
// drops the oldest when it is full.

#include "app/gui/GlLoader.h"

#include <chrono>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

namespace gui {

// Longest side of a stored thumbnail, and how many are kept. 256 x 24 is under
// 5 MB and still reads at the size the panel draws.
inline constexpr int kFilmMaxSide = 256;
inline constexpr int kFilmSlots = 24;
inline constexpr size_t kMaxOverlay = 250;

// What to draw on top of a frame: feature positions in the source image's own
// pixels. Empty for the steps that have nothing to overlay.
struct FramePoints {
    const float* xy = nullptr;   // 2 floats per point
    size_t count = 0;
};

class FilmStrip;

// The strips a run feeds as it works, one per step that produces pictures.
// Either may be null: a caller with no screen wants the pipeline without them.
struct RunFilms {
    FilmStrip* frames = nullptr;
    FilmStrip* masks = nullptr;
};

class FilmStrip {
public:
    ~FilmStrip();

    // ---- producer (worker thread) ----
    // Is another frame wanted yet? False most of the time, which is what keeps
    // the downscale off the hot path entirely.
    bool wants(double min_interval_s = 0.12) const;
    // `rgb` is w*h*3 tightly packed; `mask`, when given, is w*h with 255 where
    // the pixel is kept. Downscales and copies, so the caller may reuse all
    // three.
    void offer(const uint8_t* rgb, int w, int h, const std::string& name,
               const uint8_t* mask = nullptr, FramePoints points = {});
    void clear();

    // ---- consumer (UI thread, GL context current) ----
    bool has_frames() const;
    // Draws the strip; returns false when there is nothing to draw. `height`
    // is the row height in pixels.
    bool draw(float height);
    void destroy_gl();

private:
    struct Slot {
        std::vector<uint8_t> rgb;    // downscaled, masked already composited
        int w = 0, h = 0;
        std::string name;
        // Overlay points as fractions of the thumbnail, already thinned to
        // kMaxOverlay: eight thousand dots on a 72-pixel picture is a
        // rectangle, and drawing them all costs a draw command each.
        std::vector<float> pts;
        uint64_t seq = 0;
        GLuint tex = 0;
        uint64_t uploaded = 0;
    };

    mutable std::mutex _mu;
    Slot _slots[kFilmSlots];
    uint64_t _seq = 0;
    // wants() is called far more often than offer(); the clock lives here so a
    // caller cannot forget to rate-limit itself.
    mutable std::chrono::steady_clock::time_point _last{};
};

}  // namespace gui
