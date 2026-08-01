#include "sam/pipeline/PostProcess.h"

#include <algorithm>
#include <vector>

namespace sam {

float box_iou(const Box& a, const Box& b) {
    const float x0 = std::max(a.x0, b.x0);
    const float y0 = std::max(a.y0, b.y0);
    const float x1 = std::min(a.x1, b.x1);
    const float y1 = std::min(a.y1, b.y1);
    const float inter = std::max(0.0f, x1 - x0) * std::max(0.0f, y1 - y0);
    const float area_a = (a.x1 - a.x0) * (a.y1 - a.y0);
    const float area_b = (b.x1 - b.x0) * (b.y1 - b.y0);
    const float uni = area_a + area_b - inter;
    return uni > 0.0f ? inter / uni : 0.0f;
}

float mask_iou(const uint8_t* a, const uint8_t* b, size_t n) {
    size_t inter = 0, uni = 0;
    for (size_t i = 0; i < n; ++i) {
        const bool va = a[i] > 127, vb = b[i] > 127;
        inter += (va && vb) ? 1 : 0;
        uni += (va || vb) ? 1 : 0;
    }
    return uni ? (float)inter / (float)uni : 0.0f;
}

std::vector<int> nms(const std::vector<Detection>& dets, float iou_threshold) {
    std::vector<int> order((size_t)dets.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return dets[(size_t)a].score > dets[(size_t)b].score; });

    std::vector<bool> suppressed(dets.size(), false);
    std::vector<int> keep;
    for (int idx : order) {
        if (suppressed[(size_t)idx]) continue;
        keep.push_back(idx);
        for (int j : order) {
            if (j == idx || suppressed[(size_t)j]) continue;
            if (box_iou(dets[(size_t)idx].box, dets[(size_t)j].box) > iou_threshold)
                suppressed[(size_t)j] = true;
        }
    }
    return keep;
}

namespace {

// 4-connected flood fill over pixels matching `foreground`, labelling as it
// goes. Returns component sizes; `touches_border` is set per component.
void label_components(const uint8_t* mask, int w, int h, bool foreground,
                      std::vector<int>& labels, std::vector<int>& sizes,
                      std::vector<bool>& touches_border) {
    const int n = w * h;
    labels.assign((size_t)n, -1);
    sizes.clear();
    touches_border.clear();
    std::vector<int> queue;
    for (int start = 0; start < n; ++start) {
        const bool fg = mask[start] > 127;
        if (fg != foreground || labels[(size_t)start] >= 0) continue;
        const int label = (int)sizes.size();
        sizes.push_back(0);
        touches_border.push_back(false);
        queue.clear();
        queue.push_back(start);
        labels[(size_t)start] = label;
        for (size_t head = 0; head < queue.size(); ++head) {
            const int p = queue[head];
            ++sizes[(size_t)label];
            const int px = p % w, py = p / w;
            if (px == 0 || px == w - 1 || py == 0 || py == h - 1)
                touches_border[(size_t)label] = true;
            const int dx[4] = {-1, 1, 0, 0};
            const int dy[4] = {0, 0, -1, 1};
            for (int d = 0; d < 4; ++d) {
                const int nx = px + dx[d], ny = py + dy[d];
                if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
                const int ni = ny * w + nx;
                if (labels[(size_t)ni] >= 0) continue;
                if ((mask[ni] > 127) != foreground) continue;
                labels[(size_t)ni] = label;
                queue.push_back(ni);
            }
        }
    }
}

}  // namespace

void fill_holes(uint8_t* mask, int w, int h, int area_threshold) {
    if (area_threshold <= 0 || w <= 0 || h <= 0) return;
    std::vector<int> labels, sizes;
    std::vector<bool> border;
    label_components(mask, w, h, /*foreground=*/false, labels, sizes, border);
    for (int i = 0; i < w * h; ++i) {
        const int l = labels[(size_t)i];
        if (l < 0) continue;
        // A background component reaching the border is the outside, not a
        // hole, whatever its area.
        if (border[(size_t)l]) continue;
        if (sizes[(size_t)l] <= area_threshold) mask[i] = 255;
    }
}

void remove_sprinkles(uint8_t* mask, int w, int h, int area_threshold) {
    if (area_threshold <= 0 || w <= 0 || h <= 0) return;
    std::vector<int> labels, sizes;
    std::vector<bool> border;
    label_components(mask, w, h, /*foreground=*/true, labels, sizes, border);
    for (int i = 0; i < w * h; ++i) {
        const int l = labels[(size_t)i];
        if (l >= 0 && sizes[(size_t)l] <= area_threshold) mask[i] = 0;
    }
}

void resolve_overlaps(std::vector<Detection>& dets) {
    if (dets.size() <= 1) return;
    std::stable_sort(dets.begin(), dets.end(),
                     [](const Detection& a, const Detection& b) { return a.score > b.score; });
    const int w = dets[0].mask.width, h = dets[0].mask.height;
    if (w <= 0 || h <= 0) return;
    const size_t n = (size_t)w * h;
    for (size_t i = 0; i < n; ++i) {
        bool claimed = false;
        for (Detection& d : dets) {
            if (d.mask.data.size() != n) continue;
            if (claimed) d.mask.data[i] = 0;
            else if (d.mask.data[i] > 127) claimed = true;
        }
    }
}

void mask_bounding_box(const Mask& mask, Box& box) {
    float x0 = 1e30f, y0 = 1e30f, x1 = -1e30f, y1 = -1e30f;
    for (size_t i = 0; i < mask.data.size(); ++i) {
        if (mask.data[i] <= 127) continue;
        const float x = (float)(int)(i % (size_t)mask.width);
        const float y = (float)(int)(i / (size_t)mask.width);
        x0 = std::min(x0, x);
        y0 = std::min(y0, y);
        x1 = std::max(x1, x);
        y1 = std::max(y1, y);
    }
    if (x0 <= x1) box = Box{x0, y0, x1, y1};
}

Box cxcywh_to_xyxy(float cx, float cy, float w, float h, int img_w, int img_h) {
    Box b;
    b.x0 = (cx - w * 0.5f) * (float)img_w;
    b.y0 = (cy - h * 0.5f) * (float)img_h;
    b.x1 = (cx + w * 0.5f) * (float)img_w;
    b.y1 = (cy + h * 0.5f) * (float)img_h;
    b.x0 = std::max(0.0f, std::min(b.x0, (float)img_w));
    b.y0 = std::max(0.0f, std::min(b.y0, (float)img_h));
    b.x1 = std::max(0.0f, std::min(b.x1, (float)img_w));
    b.y1 = std::max(0.0f, std::min(b.y1, (float)img_h));
    return b;
}

}  // namespace sam
