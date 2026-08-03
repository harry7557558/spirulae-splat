// ssplat-sam -- segmentation, tracking and frame extraction from a shell.
//
//   ssplat-sam devices
//   ssplat-sam segment --model m.ggml --image cat.jpg --text "cat" --out out/
//   ssplat-sam segment --model m.ggml --image cat.jpg --point 315,250 --out out/
//   ssplat-sam track   --model m.ggml --frames frames/ --text "person" --out out/
//   ssplat-sam video   --info clip.mp4
//   ssplat-sam extract clip.mp4 --skip 30 --model m.ggml --text "person"
//
// The GUI drives the same library in-process (src/app/gui/DatasetPrep.cpp);
// this is the scriptable face of it, and how a masking or extraction problem
// gets reproduced without the GUI in the way.
//
// `video` and `extract` need the in-process decoder, which is only built with
// -DSSPLAT_ENABLE_PATENTED=ON (src/app/cli/sam_extract.cpp); without it they
// say so and point at ffmpeg.
//
// Diagnostics go to stderr, the result table to stdout, so a run pipes cleanly:
//   ssplat-sam segment ... 2>/dev/null > detections.tsv

#include "app/Tools.h"
#include "app/WriterPool.h"
#include "nn/core/Log.h"
#include "nn/Device.h"
#include "nn/io/Image.h"
#include "sam/Masking.h"
#include "sam/Sam.h"
#ifdef SSPLAT_HAVE_VIDEO
#include "video/Video.h"
#endif

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <future>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

void usage() {
    // Written against the historical name; app::help_text swaps in how this
    // tool was actually invoked ("ssplat sam", normally).
    static const char* kUsage =
        "ssplat-sam -- SAM 2 / SAM 3 segmentation and tracking on Vulkan\n"
        "\n"
        "  ssplat-sam devices\n"
        "        List Vulkan devices and whether each meets the baseline.\n"
        "\n"
        "  ssplat-sam segment --model <file> --image <file> [options]\n"
        "        --text <phrase>        concept prompt (all matching instances)\n"
        "        --box x0,y0,x1,y1      exemplar box; repeatable\n"
        "        --neg-box x0,y0,x1,y1  negative exemplar box; repeatable\n"
        "        --point x,y            positive click; repeatable (visual prompt)\n"
        "        --neg-point x,y        negative click; repeatable\n"
        "        --prompt-box x0,y0,x1,y1   box prompt for the visual path\n"
        "        --multimask            return the three ambiguity masks\n"
        "        --threshold <f>        detection score threshold (default 0.5)\n"
        "        --nms <f>              NMS IoU threshold (default 0.1)\n"
        "        --out <dir>            write mask PNGs and an overlay\n"
        "\n"
        "  ssplat-sam track --model <file> --frames <dir> [options]\n"
        "        --text <phrase>        detect and track matching instances;\n"
        "                               semicolon-separated for several concepts\n"
        "        --neg-text <phrase>    concepts to KEEP even where --text matches\n"
        "        --point x,y            click on an object to track; repeatable\n"
        "        --neg-point x,y        click on something that is NOT it\n"
        "        --object               end this object, start the next one --\n"
        "                               two things need two objects, since one\n"
        "                               instance prompted with both fits neither\n"
        "        --at-frame <n>         put the clicks that follow on frame n\n"
        "                               instead of the first, and use them to\n"
        "                               correct the object there. Frames are\n"
        "                               numbered as this command reads them:\n"
        "                                 --point 640,360 --at-frame 90 --point 700,300\n"
        "                                 --object --point 120,500\n"
        "                               is one object clicked once and corrected\n"
        "                               at frame 90, and a second object.\n"
        "        --detect-every <n>     run the detector every n frames (default 1);\n"
        "                               the memory bank carries tracks in between\n"
        "        --memory-frames <n>    cap spatial memory frames per instance\n"
        "        --max-frames <n>       stop after n frames\n"
        "        --out <dir>            write a per-frame binary mask PNG\n"
        "        --keep-prompted        white = the prompted objects. By default\n"
        "                               they are BLACK and everything else is\n"
        "                               white, which is what a reconstruction\n"
        "                               pipeline wants from \"mask out the people\"\n"
        "        --overlay              write a colour overlay instead\n"
        "\n"
        "  ssplat-sam video --info <file>\n"
        "        Probe a video file and report codec, geometry and decode support.\n"
        "\n"
        "  ssplat-sam extract <video> [options]\n"
        "        Write the sharpest frames of a video, optionally masked.\n"
        "        `ssplat-sam extract --help` lists its own options.\n"
        "\n"
        "Common: --device <index|name>  --vram  --profile  --validate  --img-size <n>\n"
        "        --max-size <n>         downscale inputs to fit (default 1600, 0 = off)\n"
        "Environment: SSPLAT_NN_LOG=0..3  SSPLAT_VK_DEVICE  SSPLAT_PROFILE=1\n"
        "             SSPLAT_VK_VALIDATION=1  SSPLAT_NN_DEBUG_SYNC=1\n";
    std::fprintf(stderr, "%s", app::help_text(kUsage, "ssplat-sam").c_str());
}

bool parse_floats(const char* s, float* out, int n) {
    for (int i = 0; i < n; ++i) {
        char* end = nullptr;
        out[i] = std::strtof(s, &end);
        if (end == s) return false;
        s = end;
        if (i + 1 < n) {
            while (*s == ' ') ++s;
            if (*s != ',' && *s != 'x') return false;
            ++s;
        }
    }
    return true;
}

struct Options {
    std::string command;
    std::string model, image, frames, out_dir, text, neg_text, video_path;
    std::string device;
    std::vector<sam::Box> pos_boxes, neg_boxes;
    std::vector<sam::Point> pos_points, neg_points;
    sam::Box prompt_box{};
    bool use_prompt_box = false;
    // Clicks for `track`, grouped into objects and frames as they are parsed.
    // `segment` works on one still image and reads the flat vectors above.
    std::vector<sam::SeedPrompt> seeds;
    int     cur_object = 0;
    int64_t cur_frame = 0;
    bool multimask = false, show_vram = false, profile = false, validate = false;
    bool overlay = false, keep_prompted = false;
    float threshold = 0.5f, nms = 0.1f;
    int max_frames = 0;
    int img_size = 0;
    int detect_every = 1, memory_frames = 0, max_size = 1600;
};

// The seed the next click goes into: one per (object, frame), created on
// demand so `--point a --point b` is two clicks on one object and not two
// objects, which is the distinction that matters to the model.
sam::SeedPrompt& current_seed(Options& o) {
    for (sam::SeedPrompt& s : o.seeds)
        if (s.object == o.cur_object && s.frame == o.cur_frame) return s;
    sam::SeedPrompt s;
    s.object = o.cur_object;
    s.frame = o.cur_frame;
    o.seeds.push_back(s);
    return o.seeds.back();
}

bool parse_args(int argc, char** argv, Options& o) {
    if (argc < 2) return false;
    o.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* what) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s needs a value\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        float v[4];
        if (a == "--model") o.model = next("--model");
        else if (a == "--image") o.image = next("--image");
        else if (a == "--frames") o.frames = next("--frames");
        else if (a == "--out") o.out_dir = next("--out");
        else if (a == "--text") o.text = next("--text");
        else if (a == "--neg-text") o.neg_text = next("--neg-text");
        else if (a == "--detect-every") o.detect_every = std::atoi(next("--detect-every"));
        else if (a == "--memory-frames") o.memory_frames = std::atoi(next("--memory-frames"));
        else if (a == "--max-size") o.max_size = std::atoi(next("--max-size"));
        else if (a == "--keep-prompted") o.keep_prompted = true;
        else if (a == "--info") o.video_path = next("--info");
        else if (a == "--device") o.device = next("--device");
        else if (a == "--threshold") o.threshold = std::strtof(next("--threshold"), nullptr);
        else if (a == "--nms") o.nms = std::strtof(next("--nms"), nullptr);
        else if (a == "--max-frames") o.max_frames = std::atoi(next("--max-frames"));
        else if (a == "--img-size") o.img_size = std::atoi(next("--img-size"));
        else if (a == "--multimask") o.multimask = true;
        else if (a == "--overlay") o.overlay = true;
        else if (a == "--vram") o.show_vram = true;
        else if (a == "--profile") o.profile = true;
        else if (a == "--validate") o.validate = true;
        else if (a == "--box" && parse_floats(next("--box"), v, 4))
            o.pos_boxes.push_back({v[0], v[1], v[2], v[3]});
        else if (a == "--neg-box" && parse_floats(next("--neg-box"), v, 4))
            o.neg_boxes.push_back({v[0], v[1], v[2], v[3]});
        else if (a == "--object") {
            // Objects are numbered from 0, so the first --object opens the
            // SECOND one and a command line that never says it still works.
            ++o.cur_object;
            o.cur_frame = 0;
        } else if (a == "--at-frame") o.cur_frame = std::atoll(next("--at-frame"));
        else if (a == "--prompt-box" && parse_floats(next("--prompt-box"), v, 4)) {
            o.prompt_box = {v[0], v[1], v[2], v[3]};
            o.use_prompt_box = true;
            sam::SeedPrompt& s = current_seed(o);
            s.prompt.box = o.prompt_box;
            s.prompt.use_box = true;
        } else if (a == "--point" && parse_floats(next("--point"), v, 2)) {
            o.pos_points.push_back({v[0], v[1]});
            current_seed(o).prompt.pos_points.push_back({v[0], v[1]});
        } else if (a == "--neg-point" && parse_floats(next("--neg-point"), v, 2)) {
            o.neg_points.push_back({v[0], v[1]});
            current_seed(o).prompt.neg_points.push_back({v[0], v[1]});
        } else {
            std::fprintf(stderr, "unknown or malformed option: %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

int cmd_devices() {
    auto devices = nn::list_devices();
    if (devices.empty()) {
        std::fprintf(stderr, "no Vulkan devices found -- is a driver installed?\n");
        return 1;
    }
    std::printf("%-3s %-42s %-11s %8s  %s\n", "idx", "name", "type", "vram", "status");
    for (const auto& d : devices)
        std::printf("%-3d %-42s %-11s %6.1f G  %s\n", d.index, d.name.c_str(),
                    d.type.c_str(), d.vram_bytes / 1073741824.0,
                    d.usable ? "ok" : d.unusable_reason.c_str());
    return 0;
}

void write_results(const Options& o, const nn::Image& image, const sam::Result& r,
                   const std::string& tag) {
    std::printf("#\tid\tscore\tiou\tx0\ty0\tx1\ty1\tarea\n");
    for (size_t i = 0; i < r.detections.size(); ++i) {
        const auto& d = r.detections[i];
        size_t area = 0;
        for (uint8_t m : d.mask.data) area += (m > 127) ? 1 : 0;
        std::printf("%zu\t%d\t%.4f\t%.4f\t%.1f\t%.1f\t%.1f\t%.1f\t%zu\n", i,
                    d.instance_id, d.score, d.iou_score, d.box.x0, d.box.y0, d.box.x1,
                    d.box.y1, area);
    }
    if (o.out_dir.empty()) return;
    std::error_code ec;
    fs::create_directories(o.out_dir, ec);
    for (size_t i = 0; i < r.detections.size(); ++i) {
        char name[512];
        std::snprintf(name, sizeof name, "%s/%s_mask_%02zu.png", o.out_dir.c_str(),
                      tag.c_str(), i);
        sam::save_mask_png(r.detections[i].mask, name);
    }
    char overlay[512];
    std::snprintf(overlay, sizeof overlay, "%s/%s_overlay.png", o.out_dir.c_str(),
                  tag.c_str());
    sam::save_overlay_png(image, r, overlay);
    std::fprintf(stderr, "wrote %zu masks + overlay to %s\n", r.detections.size(),
                 o.out_dir.c_str());
}

bool load_session(const Options& o, sam::Session& session) {
    sam::ModelParams mp;
    mp.model_path = o.model;
    mp.validation = o.validate;
    mp.profile = o.profile;
    mp.img_size = o.img_size;
    if (!o.device.empty()) {
        char* end = nullptr;
        long idx = std::strtol(o.device.c_str(), &end, 10);
        if (end && *end == '\0') mp.device_index = (int)idx;
        else mp.device_match = o.device;
    }
    if (!session.loadModel(mp)) {
        std::fprintf(stderr, "error: %s\n", session.lastError().c_str());
        return false;
    }
    return true;
}

int cmd_segment(const Options& o) {
    if (o.model.empty() || o.image.empty()) {
        std::fprintf(stderr, "segment needs --model and --image\n");
        return 2;
    }
    sam::Session session;
    if (!load_session(o, session)) return 1;

    nn::Image image = nn::load_image(o.image);
    if (image.empty()) return 1;
    if (!session.encodeImage(image)) {
        std::fprintf(stderr, "error: %s\n", session.lastError().c_str());
        return 1;
    }

    sam::Result r;
    const bool visual = !o.pos_points.empty() || !o.neg_points.empty() || o.use_prompt_box;
    if (visual) {
        sam::VisualPrompt vp;
        vp.pos_points = o.pos_points;
        vp.neg_points = o.neg_points;
        vp.box = o.prompt_box;
        vp.use_box = o.use_prompt_box;
        vp.multimask = o.multimask;
        r = session.segmentVisual(vp);
    } else {
        sam::ConceptPrompt cp;
        cp.text = o.text;
        cp.pos_exemplars = o.pos_boxes;
        cp.neg_exemplars = o.neg_boxes;
        cp.score_threshold = o.threshold;
        cp.nms_threshold = o.nms;
        r = session.segmentConcept(cp);
    }
    if (r.detections.empty() && !session.lastError().empty()) {
        std::fprintf(stderr, "error: %s\n", session.lastError().c_str());
        return 1;
    }
    write_results(o, image, r, "seg");
    if (o.show_vram) std::fprintf(stderr, "VRAM:\n%s", session.vramReport().c_str());
    if (o.profile) session.printProfile();
    return 0;
}

int cmd_track(const Options& o) {
    if (o.model.empty() || o.frames.empty()) {
        std::fprintf(stderr, "track needs --model and --frames <dir>\n");
        return 2;
    }
    std::vector<std::string> files;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(o.frames, ec)) {
        if (!e.is_regular_file()) continue;
        std::string ext = e.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".png" || ext == ".jpg" || ext == ".jpeg" || ext == ".bmp")
            files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) {
        std::fprintf(stderr, "no images in %s\n", o.frames.c_str());
        return 1;
    }
    if (o.max_frames > 0 && (int)files.size() > o.max_frames)
        files.resize((size_t)o.max_frames);

    // The mask policy -- several concepts, negatives that carve back out, and
    // which side of the mask is white -- lives in sam/Masking.h, shared with
    // `extract` and with the GUI so none of the three can disagree.
    sam::MaskOptions mo;
    mo.model = o.model;
    mo.device = o.device;
    mo.text = o.text;
    mo.neg_text = o.neg_text;
    mo.keep_prompted = o.keep_prompted;
    mo.threshold = o.threshold;
    mo.nms = o.nms;
    mo.detect_every = o.detect_every;
    mo.memory_frames = o.memory_frames;
    mo.max_size = o.max_size;
    mo.img_size = o.img_size;
    mo.validate = o.validate;
    mo.profile = o.profile;
    mo.seeds = o.seeds;

    sam::Masker masker;
    std::string error;
    if (!masker.init(mo, error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }

    // Neither reading the next frame nor writing the last mask needs the GPU,
    // and together they are about a third of a frame -- see app/WriterPool.h.
    app::WriterPool writers;
    std::future<nn::Image> ahead;
    auto load_at = [&](size_t i) {
        return std::async(std::launch::async,
                          [p = files[i]] { return nn::load_image(p); });
    };
    if (!files.empty()) ahead = load_at(0);

    double t_load = 0, t_track = 0, t_write = 0;
    const double t_all = nn::now_ms();
    for (size_t f = 0; f < files.size(); ++f) {
        double t0 = nn::now_ms();
        nn::Image frame = ahead.get();
        if (f + 1 < files.size()) ahead = load_at(f + 1);
        t_load += nn::now_ms() - t0;
        if (frame.empty()) continue;

        t0 = nn::now_ms();
        sam::Mask mask;
        sam::Result r;
        if (!masker.run(frame, mask, o.overlay ? &r : nullptr, (int64_t)f)) {
            std::fprintf(stderr, "error: %s\n", masker.lastError().c_str());
            return 1;
        }
        t_track += nn::now_ms() - t0;
        t0 = nn::now_ms();

        size_t white = 0;
        for (uint8_t v : mask.data) white += (v > 127) ? 1 : 0;
        std::fprintf(stderr, "[frame %3zu] %-40s %5.1f%% kept\n", f,
                     fs::path(files[f]).filename().string().c_str(),
                     100.0 * (double)white / (double)std::max<size_t>(mask.data.size(), 1));
        if (!o.out_dir.empty()) {
            std::error_code dec;
            fs::create_directories(o.out_dir, dec);
            char path[512];
            std::snprintf(path, sizeof path, "%s/frame_%05zu.png", o.out_dir.c_str(), f);
            if (o.overlay) {
                sam::save_overlay_png(frame, r, path);
            } else {
                app::WriteJob wj;
                wj.mask = std::move(mask);
                wj.path = path;
                writers.submit(std::move(wj));
            }
        }
        t_write += nn::now_ms() - t0;
    }
    writers.finish();
    const double total = nn::now_ms() - t_all;
    const double n = (double)files.size();
    std::fprintf(stderr,
                 "%zu frames in %.1f s -- %.0f ms/frame "
                 "(decode %.0f, model %.0f, write %.0f)\n",
                 files.size(), total / 1000.0, total / n, t_load / n, t_track / n,
                 t_write / n);
    if (o.show_vram)
        std::fprintf(stderr, "VRAM:\n%s", masker.session().vramReport().c_str());
    if (o.profile) masker.session().printProfile();
    return 0;
}

int cmd_video(const Options& o) {
#ifndef SSPLAT_HAVE_VIDEO
    (void)o;
    std::fprintf(stderr,
                 "this build has no in-process video decoder: it is compiled "
                 "only with -DSSPLAT_ENABLE_PATENTED=ON (see "
                 "cmake/SsplatOptions.cmake). Use ffmpeg to inspect a file or "
                 "extract frames instead.\n");
    return 1;
#else
    const std::string why = video::VideoReader::availability();
    if (why.empty()) std::fprintf(stderr, "Vulkan video decode: available\n");
    else std::fprintf(stderr, "Vulkan video decode: %s\n", why.c_str());

    if (o.video_path.empty()) return why.empty() ? 0 : 1;

    video::VideoReader reader;
    if (!reader.open(o.video_path)) {
        std::fprintf(stderr, "error: %s\n", reader.lastError().c_str());
        return 1;
    }
    const auto& i = reader.info();
    std::printf("codec\t%s\nwidth\t%d\nheight\t%d\nframes\t%d\nfps\t%.3f\n",
                i.codec.c_str(), i.width, i.height, i.frame_count, i.fps);

    // A decode smoke test: --max-frames N decodes the first N frames, writing
    // them as PNGs if --out is given. This is how a codec path gets compared
    // against a reference decoder.
    if (o.max_frames > 0) {
        const double t0 = nn::now_ms();
        int n = 0;
        for (; n < o.max_frames; ++n) {
            nn::Image frame = reader.readFrame();
            if (frame.empty()) {
                if (!reader.lastError().empty())
                    std::fprintf(stderr, "error: %s\n", reader.lastError().c_str());
                break;
            }
            if (o.out_dir.empty()) continue;
            char path[512];
            std::snprintf(path, sizeof(path), "%s/frame_%05d.png", o.out_dir.c_str(), n);
            nn::save_image(frame, path, -1);
            std::printf("wrote %s\n", path);
        }
        const double ms = nn::now_ms() - t0;
        std::printf("decoded\t%d frames in %.0f ms (%.1f fps)\n", n, ms,
                    n > 0 ? 1000.0 * n / ms : 0.0);
    }
    return 0;
#endif
}

}  // namespace

// Defined in src/app/cli/sam_extract.cpp; it has its own option set, so it
// parses its own argv rather than sharing Options above.
#ifdef SSPLAT_HAVE_VIDEO
int sam_cli_extract(int argc, char** argv);
#endif

int ssplat_sam_main(int argc, char** argv) {
    app::set_program_name(argc > 0 ? argv[0] : nullptr, "ssplat sam");
    if (argc >= 2 && std::strcmp(argv[1], "extract") == 0) {
#ifdef SSPLAT_HAVE_VIDEO
        int rc = sam_cli_extract(argc - 1, argv + 1);
        nn::shutdown();
        return rc;
#else
        std::fprintf(stderr,
                     "`extract` needs the in-process video decoder, which is "
                     "compiled only with -DSSPLAT_ENABLE_PATENTED=ON (see "
                     "cmake/SsplatOptions.cmake). Extract frames with ffmpeg "
                     "and mask them with `%s track` instead.\n",
                     app::program_name().c_str());
        return 1;
#endif
    }

    Options o;
    if (!parse_args(argc, argv, o)) {
        usage();
        return 2;
    }
    int rc = 2;
    if (o.command == "devices")      rc = cmd_devices();
    else if (o.command == "segment") rc = cmd_segment(o);
    else if (o.command == "track")   rc = cmd_track(o);
    else if (o.command == "video")   rc = cmd_video(o);
    else usage();
    nn::shutdown();
    return rc;
}
