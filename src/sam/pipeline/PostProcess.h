#pragma once
// Host-side mask and box post-processing.
//
// All of it is small, branchy, whole-image bookkeeping -- connected components,
// greedy NMS, priority resolution between overlapping instances. None of it
// belongs on the GPU: the masks are already bytes by the time they get here
// (the resize+threshold runs on device), and a 1920x1080 flood fill costs well
// under a millisecond.

#include "sam/Sam.h"

#include <vector>

namespace sam {

// Greedy NMS by box IoU, highest score first. Returns the kept indices.
std::vector<int> nms(const std::vector<Detection>& dets, float iou_threshold);

float box_iou(const Box& a, const Box& b);
float mask_iou(const uint8_t* a, const uint8_t* b, size_t n);

// Fills background components no larger than `area_threshold` that do not touch
// the image border. Border-touching components are the background itself and
// are always kept, however small the visible part is.
void fill_holes(uint8_t* mask, int w, int h, int area_threshold);

// Drops foreground components no larger than `area_threshold`.
void remove_sprinkles(uint8_t* mask, int w, int h, int area_threshold);

// Makes a set of instance masks mutually exclusive: at each pixel the
// highest-scoring instance keeps it. Sorts `dets` by score as a side effect,
// which is also the order results are returned in.
void resolve_overlaps(std::vector<Detection>& dets);

// Tight bounding box of a binary mask; leaves `box` untouched when empty.
void mask_bounding_box(const Mask& mask, Box& box);

// cxcywh in [0,1] -> xyxy in pixels, clamped to the image.
Box cxcywh_to_xyxy(float cx, float cy, float w, float h, int img_w, int img_h);

}  // namespace sam
