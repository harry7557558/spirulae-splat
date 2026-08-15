#pragma once

// Features drawn over a picture, the way OpenCV's drawKeypoints draws them: a
// circle over the region the descriptor was measured across rather than a dot
// on the pixel it sits on. The size is the thing a dot cannot say -- a frame
// whose features are all tiny is a frame matched on texture noise.
//
// Shared by the run's film reel and the match map's pair view, which draw the
// same overlay over pictures of different sizes.

#include "app/gui/Layout.h"         // px
#include "app/gui/SfmProgress.h"    // KeyPoint2D

#include "imgui.h"

#include <algorithm>
#include <vector>

namespace gui {

// `pts` are fractions of the picture, whose rectangle on screen is `at` and
// `size`. A circle smaller than a few pixels is a dot with extra draw calls,
// and one larger than the frame is a hoop, so the radius is clamped.
inline void draw_keypoints(ImDrawList* dl, const std::vector<KeyPoint2D>& pts,
                           const ImVec2& at, const ImVec2& size, ImU32 color) {
    const float lo = px(2.5f), hi = px(16.0f);
    const float span = std::max(size.x, size.y);
    // A dark ring under the bright one: the same green is invisible over grass
    // and over a white wall in turn.
    const ImU32 shade = IM_COL32(0, 0, 0, 90);
    for (const KeyPoint2D& k : pts) {
        const ImVec2 c(at.x + k.x * size.x, at.y + k.y * size.y);
        const float r = std::clamp(k.r * span, lo, hi);
        dl->AddCircle(c, r + 1.0f, shade, 0, 1.0f);
        dl->AddCircle(c, r, color, 0, px(1.4f));
    }
}

}  // namespace gui
