#pragma once
// Prompt-to-mask policy, shared by ssam-cli and ssam-extract.
//
// This is the semantic layer spirulae-splat/scripts/mask.py implements in
// Python, kept in one place so the two tools cannot drift:
//
//   * several positive phrases, semicolon separated, unioned;
//   * negative phrases carved back out -- a region matching one is KEPT even
//     when it also matches a positive phrase;
//   * the output mask says what to KEEP, so by default the prompted objects
//     come out black and everything else white. `keep_prompted` flips that,
//     for the case where the prompt names the subject rather than a distractor;
//   * a longest-side cap on what the model sees, with the mask returned at the
//     source resolution.
//
// It also covers the visual path, where there is no text at all and one
// instance is seeded from clicks on the first frame.

#include "sam/Sam.h"

#include <memory>
#include <string>
#include <vector>

namespace sam {

struct MaskOptions {
    std::string model, device;
    std::string text, neg_text;
    // Clicks/box seeding one instance on the first frame. Used when there is no
    // text prompt, which is the only option for a SAM 2 checkpoint.
    sam::VisualPrompt seed;
    bool  has_seed = false;

    bool  video = true;           // track across frames vs. segment each alone
    bool  keep_prompted = false;  // white = the prompted objects
    float threshold = 0.5f, nms = 0.1f;
    int   detect_every = 1;
    int   memory_frames = 0;
    // Longest side handed to the model; 0 = off. Masks still come back at the
    // source resolution. mask.py's --max_image_size, same default.
    int   max_size = 1600;
    int   img_size = 0;
    bool  validate = false;
};

// "person; car;  bicycle" -> {"person", "car", "bicycle"}. Blank entries drop.
std::vector<std::string> split_phrases(const std::string& s);

// Longest-side cap. Returns `src` unchanged when it already fits.
nn::Image downscale_to_fit(const nn::Image& src, int max_size);

class Masker {
public:
    Masker();
    ~Masker();
    Masker(const Masker&) = delete;
    Masker& operator=(const Masker&) = delete;

    bool init(const MaskOptions& o, std::string& error);

    // Writes a mask the size of `image`: 255 where the pixel should be KEPT.
    // `overlay_out`, when given, receives the raw detections for diagnostics.
    bool run(const nn::Image& image, sam::Mask& out, sam::Result* overlay_out);

    const std::string& lastError() const;
    sam::Session& session();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace sam
