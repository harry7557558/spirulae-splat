// Metric3D checkpoint and forward-pass gate.
//
// What it can check on its own is that a real checkpoint parses, that every
// weight the forward pass asks for exists at the width the rest of the model
// assumes, and that a prediction comes back finite and in range. It cannot
// check the numbers: we do not own these weights and cannot embed a golden
// copy. That is tools/metric3d/compare_ort.py's job, and it is the gate that
// matters -- run it against the same .onnx before believing a change.
//
//   ./build/metric3d_test                       # cached checkpoints, or SKIP
//   ./build/metric3d_test --fetch               # download them first
//   ./build/metric3d_test --model M --image IMG.jpg [--max-size N] [--repeat N]
//   SS_METRIC3D_DUMP=/tmp/ours ./build/metric3d_test --model M --image IMG.jpg

#include "metric3d/Metric3D.h"
#include "metric3d/model/Fetch.h"

#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "nn/io/Image.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

int g_failures = 0;

void check(bool ok, const char* what) {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) ++g_failures;
}

// Longest-side cap, rounded down to the patch stride so nothing is cropped
// after the resize. Nearest-neighbour on purpose: this is a test harness, and
// a resampling difference against the reference would be indistinguishable
// from a model difference.
std::vector<float> load_scaled(const std::string& path, int max_size, int& w, int& h) {
    const nn::Image img = nn::load_image(path);
    NN_CHECK(!img.empty(), "cannot read '%s'", path.c_str());
    const int P = metric3d::Predictor::sizeGranularity();
    double s = 1.0;
    if (max_size > 0 && std::max(img.width, img.height) > max_size)
        s = (double)max_size / std::max(img.width, img.height);
    w = std::max(2, (int)(img.width * s)) / P * P;
    h = std::max(2, (int)(img.height * s)) / P * P;

    std::vector<float> out((size_t)w * h * 3);
    for (int y = 0; y < h; ++y) {
        const int sy = std::min(img.height - 1, (int)((y + 0.5) / h * img.height));
        for (int x = 0; x < w; ++x) {
            const int sx = std::min(img.width - 1, (int)((x + 0.5) / w * img.width));
            const uint8_t* p = &img.data[((size_t)sy * img.width + sx) * 3];
            for (int c = 0; c < 3; ++c)
                out[((size_t)y * w + x) * 3 + c] = p[c] * (1.0f / 255.0f);
        }
    }
    return out;
}

void run_one(const std::string& model, const std::string& image, int max_size,
             int repeat) {
    std::printf("%s\n", model.c_str());
    const double t0 = nn::now_ms();
    metric3d::Predictor pred;
    pred.load(model);
    check(pred.loaded(), "checkpoint loads");
    std::printf("  weights %.1f MiB, loaded in %.2f s\n",
                (double)pred.deviceBytes() / 1048576.0, (nn::now_ms() - t0) / 1000.0);
    if (image.empty()) return;

    int w = 0, h = 0;
    const std::vector<float> rgb = load_scaled(image, max_size, w, h);
    std::printf("  %s at %dx%d\n", image.c_str(), w, h);
    metric3d::Prediction out = pred.predict(rgb.data(), w, h);
    // The first call builds pipelines and measures the GEMM tiling; the rate
    // a dataset actually runs at is the second one onwards.
    for (int i = 1; i < repeat; ++i) {
        const double t = nn::now_ms();
        out = pred.predict(rgb.data(), w, h);
        std::printf("  predict %.0f ms\n", nn::now_ms() - t);
    }

    check(out.width > 0 && out.height > 0, "prediction has a size");
    check((size_t)out.depth.size() == (size_t)out.width * out.height, "depth is full size");
    check((size_t)out.normal.size() == (size_t)out.width * out.height * 3,
          "normal is full size");

    double dmin = 1e30, dmax = -1e30, nmin = 1e30, nmax = -1e30;
    bool finite = true;
    for (float v : out.depth) {
        finite &= std::isfinite(v);
        dmin = std::min(dmin, (double)v);
        dmax = std::max(dmax, (double)v);
    }
    for (size_t i = 0; i + 2 < out.normal.size(); i += 3) {
        const double n = std::sqrt((double)out.normal[i] * out.normal[i] +
                                   (double)out.normal[i + 1] * out.normal[i + 1] +
                                   (double)out.normal[i + 2] * out.normal[i + 2]);
        finite &= std::isfinite(n);
        nmin = std::min(nmin, n);
        nmax = std::max(nmax, n);
    }
    check(finite, "every output is finite");
    check(dmin >= 0.0999 && dmax <= 200.001, "depth is inside the network's own range");
    check(nmin > 0.99 && nmax < 1.01, "normals are unit length");
    std::printf("  depth %.3f..%.3f (canonical; metres = d * fx/1000)\n", dmin, dmax);
}

}  // namespace

int main(int argc, char** argv) {
    std::string model, image;
    int max_size = 616, repeat = 1;
    bool fetch = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--model" && i + 1 < argc) model = argv[++i];
        else if (a == "--image" && i + 1 < argc) image = argv[++i];
        else if (a == "--max-size" && i + 1 < argc) max_size = std::atoi(argv[++i]);
        else if (a == "--repeat" && i + 1 < argc) repeat = std::atoi(argv[++i]);
        else if (a == "--fetch") fetch = true;
        else {
            std::fprintf(stderr, "unknown argument '%s'\n", a.c_str());
            return 2;
        }
    }

    try {
        if (!model.empty()) {
            run_one(model, image, max_size, repeat);
        } else {
            int ran = 0;
            for (const char* id : {"metric3d-vit-small", "metric3d-vit-large",
                                   "metric3d-vit-giant2"}) {
                const metric3d::ModelSource* src = metric3d::find_model_source(id);
                if (!src) continue;
                std::error_code ec;
                if (!fetch && !fs::exists(nn::cached_path(src->onnx), ec)) continue;
                run_one(id, image, max_size, repeat);
                ++ran;
            }
            if (ran == 0) {
                std::printf("SKIP: no cached checkpoint (--fetch downloads them, or "
                            "--model <path.onnx>)\n");
                return 0;
            }
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "\nerror: %s\n", e.what());
        return 1;
    }

    std::printf("\n%s\n", g_failures ? "FAILURES" : "all checks passed");
    return g_failures ? 1 : 0;
}
