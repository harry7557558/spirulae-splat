#include "sam/Masking.h"

#include "nn/core/Parallel.h"

#include <algorithm>
#include <cstring>

namespace sam {

std::vector<std::string> split_phrases(const std::string& s) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size()) {
        const size_t end = s.find(';', start);
        std::string p = s.substr(start, end == std::string::npos ? std::string::npos
                                                                 : end - start);
        const size_t b = p.find_first_not_of(" \t");
        const size_t e = p.find_last_not_of(" \t");
        if (b != std::string::npos) out.push_back(p.substr(b, e - b + 1));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return out;
}

nn::Image downscale_to_fit(const nn::Image& src, int max_size) {
    const int longest = std::max(src.width, src.height);
    if (max_size <= 0 || longest <= max_size || src.empty()) return src;

    // PIL's Image.resize((w, h)) rounds the target down, as mask.py's
    // int(size * scale) does; the sampling below is the plain box-free bilinear
    // the model's own preprocessing would have applied anyway.
    const double sc = (double)max_size / (double)longest;
    nn::Image out;
    out.width = std::max(1, (int)(src.width * sc));
    out.height = std::max(1, (int)(src.height * sc));
    out.channels = src.channels;
    out.data.resize((size_t)out.width * out.height * 3);

    const double sx = (double)src.width / out.width;
    const double sy = (double)src.height / out.height;
    // Row-parallel: at 4K in, this is tens of milliseconds of the per-frame
    // budget on one thread and a rounding error on many, and every caller runs
    // it once per frame with the GPU idle behind it.
    nn::parallel_for(out.height, [&](int64_t y_lo, int64_t y_hi) {
    for (int y = (int)y_lo; y < (int)y_hi; ++y) {
        double fy = (y + 0.5) * sy - 0.5;
        if (fy < 0.0) fy = 0.0;
        const int y0 = (int)fy;
        const int y1 = (y0 < src.height - 1) ? y0 + 1 : y0;
        const double wy = fy - y0;
        for (int x = 0; x < out.width; ++x) {
            double fx = (x + 0.5) * sx - 0.5;
            if (fx < 0.0) fx = 0.0;
            const int x0 = (int)fx;
            const int x1 = (x0 < src.width - 1) ? x0 + 1 : x0;
            const double wx = fx - x0;
            for (int c = 0; c < 3; ++c) {
                const double p00 = src.data[((size_t)y0 * src.width + x0) * 3 + c];
                const double p01 = src.data[((size_t)y0 * src.width + x1) * 3 + c];
                const double p10 = src.data[((size_t)y1 * src.width + x0) * 3 + c];
                const double p11 = src.data[((size_t)y1 * src.width + x1) * 3 + c];
                const double v = (1 - wy) * ((1 - wx) * p00 + wx * p01) +
                                 wy * ((1 - wx) * p10 + wx * p11);
                out.data[((size_t)y * out.width + x) * 3 + c] =
                    (uint8_t)std::min(255, std::max(0, (int)(v + 0.5)));
            }
        }
    }
    }, /*min_chunk=*/16);
    return out;
}

namespace {

// Nearest-neighbour, matching mask.py's Image.NEAREST upscale back to the
// source resolution. Anything smoother would put grey along every edge of a
// mask that is by definition binary.
void upscale_nearest(const std::vector<uint8_t>& src, int sw, int sh,
                     std::vector<uint8_t>& dst, int dw, int dh) {
    dst.resize((size_t)dw * dh);
    nn::parallel_for(dh, [&](int64_t y_lo, int64_t y_hi) {
        for (int y = (int)y_lo; y < (int)y_hi; ++y) {
            const int sy = std::min(sh - 1, (int)((int64_t)y * sh / dh));
            for (int x = 0; x < dw; ++x) {
                const int sx = std::min(sw - 1, (int)((int64_t)x * sw / dw));
                dst[(size_t)y * dw + x] = src[(size_t)sy * sw + sx];
            }
        }
    }, /*min_chunk=*/64);
}

}  // namespace

struct Masker::Impl {
    sam::Session session;
    std::vector<std::unique_ptr<sam::Tracker>> trackers;
    std::vector<std::string> pos, neg;
    MaskOptions opts;
    int frame = 0;

    void accumulate(const sam::Result& r, std::vector<uint8_t>& hit, size_t n) {
        for (const auto& d : r.detections) {
            if (d.mask.data.size() != n) continue;
            for (size_t i = 0; i < n; ++i)
                if (d.mask.data[i] > 127) hit[i] = 1;
        }
    }
    static void append(sam::Result& all, const sam::Result& r) {
        all.detections.insert(all.detections.end(), r.detections.begin(),
                              r.detections.end());
    }
};

Masker::Masker() : impl_(new Impl()) {}
Masker::~Masker() = default;

const std::string& Masker::lastError() const { return impl_->session.lastError(); }
sam::Session& Masker::session() { return impl_->session; }

bool Masker::init(const MaskOptions& o, std::string& error) {
    impl_->opts = o;
    impl_->pos = split_phrases(o.text);
    impl_->neg = split_phrases(o.neg_text);
    if (impl_->pos.empty() && !o.has_seed) {
        error = "nothing to segment: give --text, or --point/--prompt-box to seed an "
                "instance on the first frame";
        return false;
    }

    sam::ModelParams mp;
    mp.model_path = o.model;
    mp.device_match = o.device;
    mp.validation = o.validate;
    mp.img_size = o.img_size;
    if (!impl_->session.loadModel(mp)) {
        error = impl_->session.lastError();
        return false;
    }
    if (!impl_->pos.empty() && !impl_->session.supportsTextPrompts()) {
        error = "this checkpoint has no text encoder, so --text cannot be used; seed an "
                "instance with --point instead (SAM 2 checkpoints are visual-only)";
        return false;
    }

    if (o.video) {
        // One tracker per positive phrase, all reading the same backbone pass.
        // A seed with no text gets one tracker that only ever propagates.
        auto make = [&](const std::string& phrase) {
            sam::VideoParams vp;
            vp.text_prompt = phrase;
            vp.score_threshold = o.threshold;
            vp.nms_threshold = o.nms;
            vp.detect_every = o.detect_every;
            vp.max_memory_frames = o.memory_frames;
            impl_->trackers.push_back(
                std::make_unique<sam::Tracker>(impl_->session, vp));
        };
        if (impl_->pos.empty()) make("");
        else for (const std::string& p : impl_->pos) make(p);
    }
    return true;
}

bool Masker::run(const nn::Image& image, sam::Mask& out, sam::Result* overlay_out) {
    Impl& s = *impl_;
    const nn::Image scaled = downscale_to_fit(image, s.opts.max_size);
    const size_t n = (size_t)scaled.width * scaled.height;
    std::vector<uint8_t> hit(n, 0);
    sam::Result all;

    if (!s.session.encodeImage(scaled)) return false;

    if (s.opts.video) {
        for (auto& t : s.trackers) {
            sam::Result r = t->trackEncoded();
            s.accumulate(r, hit, n);
            Impl::append(all, r);
        }
        // The seed lands after the first frame's propagation, which is a no-op
        // on an empty tracker -- the instance is then carried from frame 1 on.
        if (s.frame == 0 && s.opts.has_seed && !s.trackers.empty()) {
            sam::VisualPrompt vp = s.opts.seed;
            const double fx = (double)scaled.width / std::max(1, image.width);
            const double fy = (double)scaled.height / std::max(1, image.height);
            for (auto& p : vp.pos_points) { p.x = (float)(p.x * fx); p.y = (float)(p.y * fy); }
            for (auto& p : vp.neg_points) { p.x = (float)(p.x * fx); p.y = (float)(p.y * fy); }
            vp.box.x0 = (float)(vp.box.x0 * fx);
            vp.box.x1 = (float)(vp.box.x1 * fx);
            vp.box.y0 = (float)(vp.box.y0 * fy);
            vp.box.y1 = (float)(vp.box.y1 * fy);
            if (s.trackers.front()->addInstance(vp) < 0) return false;
        }
    } else {
        for (const std::string& phrase : s.pos) {
            sam::ConceptPrompt cp;
            cp.text = phrase;
            cp.score_threshold = s.opts.threshold;
            cp.nms_threshold = s.opts.nms;
            sam::Result r = s.session.segmentConcept(cp);
            s.accumulate(r, hit, n);
            Impl::append(all, r);
        }
        if (s.opts.has_seed) {
            sam::Result r = s.session.segmentVisual(s.opts.seed);
            s.accumulate(r, hit, n);
            Impl::append(all, r);
        }
    }

    // Negative phrases carve back out of what the positives found. The features
    // are already on the device, so this is a decoder pass only.
    for (const std::string& phrase : s.neg) {
        sam::ConceptPrompt cp;
        cp.text = phrase;
        cp.score_threshold = s.opts.threshold;
        cp.nms_threshold = s.opts.nms;
        sam::Result r = s.session.segmentConcept(cp);
        for (const auto& d : r.detections) {
            if (d.mask.data.size() != n) continue;
            for (size_t i = 0; i < n; ++i)
                if (d.mask.data[i] > 127) hit[i] = 0;
        }
    }

    std::vector<uint8_t> mask(n);
    for (size_t i = 0; i < n; ++i)
        mask[i] = (uint8_t)(((hit[i] != 0) == s.opts.keep_prompted) ? 255 : 0);

    out.width = image.width;
    out.height = image.height;
    if (scaled.width == image.width && scaled.height == image.height)
        out.data = std::move(mask);
    else
        upscale_nearest(mask, scaled.width, scaled.height, out.data, image.width,
                        image.height);

    if (overlay_out) {
        // The overlay is drawn over the caller's image, so per-instance masks
        // have to come back to its resolution too.
        if (scaled.width != image.width || scaled.height != image.height) {
            for (auto& d : all.detections) {
                if (d.mask.data.size() != n) continue;
                std::vector<uint8_t> up;
                upscale_nearest(d.mask.data, scaled.width, scaled.height, up, image.width,
                                image.height);
                d.mask.data = std::move(up);
                d.mask.width = image.width;
                d.mask.height = image.height;
                const double fx = (double)image.width / scaled.width;
                const double fy = (double)image.height / scaled.height;
                d.box.x0 = (float)(d.box.x0 * fx);
                d.box.x1 = (float)(d.box.x1 * fx);
                d.box.y0 = (float)(d.box.y0 * fy);
                d.box.y1 = (float)(d.box.y1 * fy);
            }
        }
        *overlay_out = std::move(all);
    }
    ++s.frame;
    return true;
}

}  // namespace sam
