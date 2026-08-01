// `ssplat-sam extract` -- pull the sharpest frames out of a video, and
// optionally mask them.
//
// The native replacement for scripts/extract_frames.py, with the
// OpenCV/ffmpeg/PyTorch stack replaced by this tree: decode is Vulkan Video,
// the motion-blur metric is a compute shader, and masking is SAM on the same
// device. Only the JPEG/PNG encode stays on the CPU, on a pool of worker
// threads, because that is what it is good at.
//
// The frame-selection arithmetic is deliberately identical to the Python's
// (a window of `keep` frames, one written every `skip`), so a dataset
// extracted either way has the same file names for the same source frames --
// which is also what lets the GUI fall back to ffmpeg + FrameSelect without
// the output changing.
//
// This file is compiled only with -DSSPLAT_ENABLE_PATENTED=ON; the dispatcher
// in sam_main.cpp explains the fallback when it is not.
//
//   ssplat-sam extract clip.mp4 --skip 10
//   ssplat-sam extract clip.mp4 --skip 10 --model sam3-f16.ggml --text "person; car"
//   ssplat-sam extract dish.mov --model sam3-f16.ggml --text food \
//                      --neg-text "cooked food" --mask-keep subject

#include "app/FrameExtract.h"
#include "app/Tools.h"
#include "nn/core/Log.h"
#include "sam/Masking.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace {

void set_env(const char* key, const char* value) {
#ifdef _WIN32
    _putenv_s(key, value);
#else
    setenv(key, value, 1);
#endif
}

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

void usage() {
    static const char* kUsage =
        "ssplat-sam extract <video> [options]\n"
        "\n"
        "Frame selection\n"
        "  -o, --out <dir>        output directory (default: <video-without-ext>/images)\n"
        "  -s, --skip <n>         write one frame every n source frames (default 1)\n"
        "  -k, --keep <n>         choose the sharpest of the last n frames;\n"
        "                         -1 = round(skip/2) (default), 0 = no selection\n"
        "  -n, --max-frames <n>   stop after writing n frames\n"
        "  -q, --quality <0..100> JPEG quality; outside that range writes PNG (default 95)\n"
        "  -r, --rotate <deg>     0, 90, 180 or 270, clockwise\n"
        "      --scale <f>        resize factor, <= 1\n"
        "      --track <i>        video track to read; default is every track,\n"
        "                         written to <out>/cam0, <out>/cam1, ...\n"
        "      --threads <n>      image-encoder threads (default: cores - 1)\n"
        "\n"
        "Masking (needs --model)\n"
        "      --model <file>     SAM 3 checkpoint\n"
        "      --text <phrases>   semicolon-separated noun phrases, e.g. \"person; car\"\n"
        "      --neg-text <ph>    phrases to subtract from the mask\n"
        "      --mask-mode <m>    video (default) tracks instances across frames;\n"
        "                         image treats every written frame independently\n"
        "      --mask-keep <w>    background (default): the prompt names distractors and\n"
        "                         everything else is kept; subject: the prompt names what\n"
        "                         to keep and everything else is masked out\n"
        "      --mask-out <dir>   default: 'masks' beside the image directory\n"
        "      --detect-every <n> run the detector every n frames (default 1); the\n"
        "                         memory bank carries instances in between\n"
        "      --memory-frames <n> cap spatial memory frames per instance; memory\n"
        "                         attention is linear in this and dominates the cost\n"
        "      --max-size <n>     longest side handed to the model (default 1600,\n"
        "                         0 = off); masks come back at frame resolution\n"
        "      --threshold <f>    detection score threshold (default 0.5)\n"
        "      --nms <f>          NMS IoU threshold (default 0.1)\n"
        "      --overlay          also write a colour overlay next to each mask\n"
        "\n"
        "Common: --device <index|name>  --profile  --validate\n";
    std::fprintf(stderr, "%s", app::help_text(kUsage, "ssplat-sam").c_str());
}

struct Options {
    std::string input;
    std::string out_dir, mask_dir;
    int    skip = 1, keep = -1, max_frames = 0;
    int    quality = 95, rotate = 0;
    float  scale = 1.0f;
    int    track = -1;
    int    threads = 0;

    std::string model, text, neg_text, device;
    std::string mask_mode = "video";
    bool   keep_subject = false;
    int    detect_every = 1, memory_frames = 0, max_size = 1600;
    float  threshold = 0.5f, nms = 0.1f;
    bool   overlay = false, profile = false, validate = false;
};

bool parse_args(int argc, char** argv, Options& o) {
    if (argc < 2) return false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* what) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s needs a value\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "-o" || a == "--out") o.out_dir = next("--out");
        else if (a == "-s" || a == "--skip") o.skip = std::atoi(next("--skip"));
        else if (a == "-k" || a == "--keep") o.keep = std::atoi(next("--keep"));
        else if (a == "-n" || a == "--max-frames") o.max_frames = std::atoi(next("--max-frames"));
        else if (a == "-q" || a == "--quality") o.quality = std::atoi(next("--quality"));
        else if (a == "-r" || a == "--rotate") o.rotate = std::atoi(next("--rotate"));
        else if (a == "--scale") o.scale = std::strtof(next("--scale"), nullptr);
        else if (a == "--track") o.track = std::atoi(next("--track"));
        else if (a == "--threads") o.threads = std::atoi(next("--threads"));
        else if (a == "--model") o.model = next("--model");
        else if (a == "--text") o.text = next("--text");
        else if (a == "--neg-text") o.neg_text = next("--neg-text");
        else if (a == "--mask-mode") o.mask_mode = next("--mask-mode");
        else if (a == "--mask-keep") o.keep_subject = std::strcmp(next("--mask-keep"), "subject") == 0;
        else if (a == "--mask-out") o.mask_dir = next("--mask-out");
        else if (a == "--detect-every") o.detect_every = std::atoi(next("--detect-every"));
        else if (a == "--memory-frames") o.memory_frames = std::atoi(next("--memory-frames"));
        else if (a == "--max-size") o.max_size = std::atoi(next("--max-size"));
        else if (a == "--threshold") o.threshold = std::strtof(next("--threshold"), nullptr);
        else if (a == "--nms") o.nms = std::strtof(next("--nms"), nullptr);
        else if (a == "--overlay") o.overlay = true;
        else if (a == "--device") o.device = next("--device");
        else if (a == "--profile") o.profile = true;
        else if (a == "--validate") o.validate = true;
        else if (a == "-h" || a == "--help") return false;
        else if (a[0] == '-') {
            std::fprintf(stderr, "unknown option %s\n", a.c_str());
            return false;
        } else if (o.input.empty()) {
            o.input = a;
        } else {
            std::fprintf(stderr, "unexpected argument %s\n", a.c_str());
            return false;
        }
    }
    if (o.input.empty()) return false;
    if (o.skip < 1) o.skip = 1;
    if (o.keep < 0) o.keep = (int)(0.5 * o.skip + 0.5);
    if (o.rotate % 90 != 0) {
        std::fprintf(stderr, "--rotate must be a multiple of 90\n");
        return false;
    }
    return true;
}

}  // namespace

// Declared in sam_main.cpp's dispatcher.
int sam_cli_extract(int argc, char** argv);

int sam_cli_extract(int argc, char** argv) {
    Options o;
    if (!parse_args(argc, argv, o)) {
        usage();
        return 2;
    }
    // The decoder and the masker each reach vk::Context::get() on their own,
    // so the device and the diagnostic switches are handed over the way the
    // context reads them rather than plumbed through two option structs.
    if (o.validate) set_env("SSPLAT_VK_VALIDATION", "1");
    if (o.profile) set_env("SSPLAT_PROFILE", "1");
    if (!o.device.empty()) set_env("SSPLAT_VK_DEVICE", o.device.c_str());

    const fs::path input(o.input);
    const fs::path base = o.out_dir.empty()
                              ? (input.parent_path() / input.stem() / "images")
                              : fs::path(o.out_dir);

    app::FrameExtractJob job;
    job.input = o.input;
    job.image_dir = base.string();
    job.mask_dir = o.mask_dir;
    job.skip = o.skip;
    job.keep = o.keep;
    job.max_frames = o.max_frames;
    job.quality = o.quality;
    job.rotate = o.rotate;
    job.scale = o.scale;
    job.track = o.track;
    job.threads = o.threads;
    job.write_overlay = o.overlay;
    if (!o.model.empty()) {
        job.mask.model = o.model;
        job.mask.device = o.device;
        job.mask.text = o.text;
        job.mask.neg_text = o.neg_text;
        job.mask.video = o.mask_mode != "image";
        job.mask.keep_prompted = o.keep_subject;
        job.mask.threshold = o.threshold;
        job.mask.nms = o.nms;
        job.mask.detect_every = o.detect_every;
        job.mask.memory_frames = o.memory_frames;
        job.mask.max_size = o.max_size;
        job.mask.validate = o.validate;
    }

    app::FrameExtractSinks sinks;
    sinks.log = [](const std::string& l) { std::fprintf(stderr, "%s\n", l.c_str()); };

    app::FrameExtractStats stats;
    std::string error;
    int rc = 0;
    if (!app::extract_frames(job, sinks, stats, error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        rc = 1;
    }
    std::printf("%s", app::format_extract_stats(stats, base.string(),
                                                !o.model.empty()).c_str());
    if (stats.write_failures) {
        std::fprintf(stderr, "warning: %d image(s) failed to write\n",
                     stats.write_failures);
        rc = 1;
    }
    return rc;
}
