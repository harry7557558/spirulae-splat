// FrameExtract.cpp -- see app/FrameExtract.h.
//
// Lifted out of the `spirula-sam extract` CLI when the GUI grew the same need:
// the decode/select/mask/write loop is the part worth having exactly once, and
// the two front ends differ only in how they report progress.

#include "app/FrameExtract.h"

#include "app/WriterPool.h"
#include "nn/core/Log.h"
#include "nn/io/Image.h"
#include "sam/Sam.h"
#include "video/Demuxer.h"
#include "video/VideoPipeline.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

namespace app {

namespace {

void log_line(const FrameExtractSinks& sinks, const std::string& s) {
    if (sinks.log) sinks.log(s);
    else NN_LOG_ERROR("%s\n", s.c_str());
}

// The encoder thread pool this depends on is app/WriterPool.h: masking and
// plain extraction both need it, so it does not live here.
using app::WriteJob;
using app::WriterPool;

// ---------------------------------------------------------------------------
// Extraction
// ---------------------------------------------------------------------------

bool extract_track(const FrameExtractJob& o, const FrameExtractSinks& sinks,
                   int track, const fs::path& image_dir,
                   const fs::path& mask_dir, sam::Masker* masker,
                   WriterPool& pool, FrameExtractStats& t,
                   std::string& error) {
    video::VideoPipeline pipe;
    // The blur window holds `keep` decoded pictures at once.
    if (!pipe.open(o.input, track, std::max(o.keep, 1), error)) return false;

    fs::create_directories(image_dir);
    if (masker) fs::create_directories(mask_dir);

    video::ConvertOpts conv;
    conv.scale = o.scale;
    conv.rotate = o.rotate;

    struct Buffered {
        video::FrameHandle h;
        int64_t index;
    };
    std::deque<Buffered> window;
    const int keep = o.keep;
    const int skip = o.skip;
    int64_t frame_count = 0;
    int written = 0;
    bool measured_pending = false;

    auto flush_window = [&](bool write) {
        if (window.empty()) return;
        if (write) {
            if (measured_pending) {
                const double t0 = nn::now_ms();
                pipe.flushSharpness();
                t.sharpness += nn::now_ms() - t0;
                measured_pending = false;
            }
            size_t best = 0;
            if (keep != 0) {
                float best_score = -1.0f;
                for (size_t i = 0; i < window.size(); ++i) {
                    const float s = pipe.sharpness(window[i].h);
                    if (s > best_score) {
                        best_score = s;
                        best = i;
                    }
                }
            }
            const Buffered& chosen = window[best];

            nn::Image image;
            double t0 = nn::now_ms();
            std::string err;
            if (!pipe.toImage(chosen.h, conv, image, err)) {
                log_line(sinks, "frame " + std::to_string(chosen.index) + ": " + err);
            } else {
                t.convert += nn::now_ms() - t0;

                char stem[64];
                std::snprintf(stem, sizeof(stem), "%05lld", (long long)chosen.index);
                if (masker) {
                    t0 = nn::now_ms();
                    sam::Mask mask;
                    sam::Result overlay;
                    // The decoded index, not the written one: a click was drawn
                    // on a frame of the video, and only one frame per sharpness
                    // window survives to be written.
                    if (masker->run(image, mask, o.write_overlay ? &overlay : nullptr,
                                    chosen.index)) {
                        t.mask += nn::now_ms() - t0;
                        WriteJob mj;
                        mj.mask = std::move(mask);
                        mj.path = (mask_dir / (std::string(stem) + ".png")).string();
                        pool.submit(std::move(mj));
                        if (o.write_overlay) {
                            sam::save_overlay_png(
                                image, overlay,
                                (mask_dir / (std::string(stem) + "_overlay.png")).string());
                        }
                    } else {
                        log_line(sinks, "frame " + std::to_string(chosen.index) +
                                                ": masking failed: " + masker->lastError());
                    }
                }

                const bool jpeg = o.quality >= 0 && o.quality <= 100;
                WriteJob job;
                job.path = (image_dir / (std::string(stem) + (jpeg ? ".jpg" : ".png"))).string();
                job.quality = o.quality;
                job.image = std::move(image);
                t0 = nn::now_ms();
                pool.submit(std::move(job));
                t.submit += nn::now_ms() - t0;
                ++written;
                ++t.written;
                if (sinks.progress) sinks.progress(t.written, t.decoded);
            }
        }
        for (auto& b : window) pipe.release(b.h);
        window.clear();
    };

    while (o.max_frames <= 0 || written < o.max_frames) {
        if (sinks.cancel && sinks.cancel->load()) {
            error = "cancelled";
            for (auto& b : window) pipe.release(b.h);
            return false;
        }
        video::FrameHandle h;
        const double t0 = nn::now_ms();
        if (!pipe.next(h, error)) {
            t.decode += nn::now_ms() - t0;
            if (!error.empty()) return false;
            break;  // end of stream
        }
        t.decode += nn::now_ms() - t0;
        ++t.decoded;

        const int64_t i = h.index;
        // A frame matters only if some write window can select it. Everything
        // else is decoded (inter prediction needs it) but never touched again.
        const bool candidate =
            (keep == 0) ? (i % skip == 0) : (((i + keep) % skip) < keep);
        if (!candidate) {
            pipe.release(h);
            ++frame_count;
            continue;
        }

        if (keep != 0) {
            const double t1 = nn::now_ms();
            pipe.queueSharpness(h);
            t.sharpness += nn::now_ms() - t1;
            measured_pending = true;
            ++t.measured;
        }
        window.push_back({h, i});
        if ((int)window.size() > std::max(keep, 1)) {
            pipe.release(window.front().h);
            window.pop_front();
        }

        // The Python increments its counter before the modulo test when a blur
        // window is in use, and after it when it is not; reproduce both so the
        // chosen frame indices match exactly.
        if (keep != 0) ++frame_count;
        const bool write_now = (frame_count % skip) == 0;
        if (keep == 0) ++frame_count;
        if (write_now) flush_window(true);
    }
    flush_window(false);
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

std::string video_decode_availability() {
    return video::VideoPipeline::availability();
}

int video_track_count(const std::string& path, std::string& error) {
    // The demuxer alone: creating a video session per track just to count them
    // would be several hundred megabytes of DPB for nothing.
    auto demux = video::open_demuxer(path, error);
    if (!demux) return 0;
    return (int)demux->tracks().size();
}

bool extract_frames(const FrameExtractJob& job_in, const FrameExtractSinks& sinks,
                    FrameExtractStats& stats, std::string& error) {
    FrameExtractJob job = job_in;
    if (job.skip < 1) job.skip = 1;
    if (job.keep < 0) job.keep = (int)(0.5 * job.skip + 0.5);
    if (job.rotate % 90 != 0) {
        error = "rotation must be a multiple of 90 degrees";
        return false;
    }

    error = video_decode_availability();
    if (!error.empty()) return false;

    std::vector<int> tracks;
    if (job.track >= 0) {
        tracks.push_back(job.track);
    } else {
        const int n = video_track_count(job.input, error);
        if (n <= 0) {
            if (error.empty()) error = "no video track in " + job.input;
            return false;
        }
        for (int i = 0; i < n; ++i) tracks.push_back(i);
    }
    stats.tracks = (int)tracks.size();

    std::unique_ptr<sam::Masker> masker;
    if (!job.mask.model.empty()) {
        masker = std::make_unique<sam::Masker>();
        if (!masker->init(job.mask, error)) return false;
    }

    stats.encoder_threads =
        job.threads > 0 ? job.threads
                        : std::max(1, (int)std::thread::hardware_concurrency() - 1);
    WriterPool pool(stats.encoder_threads);

    const fs::path base(job.image_dir);
    const fs::path mask_base(job.mask_dir.empty() ? (base.parent_path() / "masks")
                                                  : fs::path(job.mask_dir));

    const double t_start = nn::now_ms();
    bool ok = true;
    for (size_t ti = 0; ti < tracks.size() && ok; ++ti) {
        // A multi-track file (an Insta360 .insv carries two fisheye streams)
        // becomes cam0/, cam1/, ... -- one camera per folder downstream.
        const bool multi = tracks.size() > 1;
        const fs::path image_dir = multi ? base / ("cam" + std::to_string(ti)) : base;
        const fs::path mask_dir =
            multi ? mask_base / ("cam" + std::to_string(ti)) : mask_base;
        if (multi)
            log_line(sinks, "track " + std::to_string(tracks[ti]) + " -> " +
                                image_dir.string());
        ok = extract_track(job, sinks, tracks[ti], image_dir, mask_dir, masker.get(),
                           pool, stats, error);
    }

    const double t_drain = nn::now_ms();
    pool.finish();
    stats.drain = nn::now_ms() - t_drain;
    stats.total = nn::now_ms() - t_start;
    stats.encode_cpu = pool.busyMs();
    stats.write_failures = pool.failures();
    return ok;
}

std::string format_extract_stats(const FrameExtractStats& s,
                                 const std::string& out_dir, bool masked) {
    char b[2048];
    std::snprintf(b, sizeof b,
        "\nframes decoded      %lld\n"
        "frames measured     %lld\n"
        "frames written      %lld -> %s\n\n"
        "decode              %8.0f ms  (%.1f fps)\n"
        "sharpness           %8.0f ms\n"
        "convert + download  %8.0f ms\n"
        "%s"
        "encode (%2d threads) %8.0f ms of CPU across workers\n"
        "  stalled on queue  %8.0f ms\n"
        "  drain at exit     %8.0f ms\n"
        "total               %8.0f ms  (%.1f frames/s written)\n",
        (long long)s.decoded, (long long)s.measured, (long long)s.written,
        out_dir.c_str(),
        s.decode, s.decoded > 0 ? 1000.0 * (double)s.decoded / std::max(s.decode, 1e-6) : 0.0,
        s.sharpness, s.convert,
        masked ? "" : "",
        s.encoder_threads, s.encode_cpu, s.submit, s.drain,
        s.total, s.written > 0 ? 1000.0 * (double)s.written / std::max(s.total, 1e-6) : 0.0);
    std::string out = b;
    if (masked) {
        char m[64];
        std::snprintf(m, sizeof m, "segmentation        %8.0f ms\n", s.mask);
        const std::string anchor = "encode (";
        out.insert(out.find(anchor), m);
    }
    return out;
}

}  // namespace app
