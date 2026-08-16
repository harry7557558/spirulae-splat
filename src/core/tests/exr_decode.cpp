// exr_decode -- the EXR reader against a corpus OpenEXR itself produced.
//
//   python tools/gen_exr_cases.py /tmp/exr_cases && ./build/exr_decode /tmp/exr_cases
//
// Every case must decode BIT-EXACTLY to its .f32 companion, on one thread and
// on many; a case named err_* must be refused with a non-empty message.

#include "core/ExrImage.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

bool read_f32(const std::string& path, std::vector<float>& out) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return false;
    const std::streamsize n = f.tellg();
    f.seekg(0);
    out.resize((size_t)n / sizeof(float));
    return (bool)f.read((char*)out.data(), n);
}

// NaN compares unequal to itself, and PIZ round-trips NaN payloads exactly, so
// the two are compared by bit pattern rather than by value.
bool same(const std::vector<float>& a, const std::vector<float>& b) {
    if (a.size() != b.size()) return false;
    return std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: exr_decode <case-dir>\n");
        return 2;
    }
    const std::string dir = argv[1];
    std::ifstream manifest(dir + "/cases.txt");
    if (!manifest) {
        std::fprintf(stderr, "cannot read %s/cases.txt -- run tools/gen_exr_cases.py\n",
                     dir.c_str());
        return 2;
    }

    int total = 0, failed = 0;
    std::string name;
    int want_w = 0, want_h = 0;
    while (manifest >> name >> want_w >> want_h) {
        total++;
        const std::string path = dir + "/" + name + ".exr";
        const bool want_error = name.rfind("err_", 0) == 0;

        std::vector<float> want;
        if (!want_error && !read_f32(dir + "/" + name + ".f32", want)) {
            std::printf("BAD  %s: no reference\n", name.c_str());
            failed++;
            continue;
        }
        for (int threads : {1, 8}) {
            exr::Info info;
            exr::Options opt;
            opt.threads = threads;
            std::vector<float> got;
            const std::string err = exr::decode(path, opt, info, got);
            if (want_error) {
                if (err.empty()) {
                    std::printf("BAD  %s t%d: decoded, but should not have\n",
                                name.c_str(), threads);
                    failed++;
                } else if (threads == 1) {
                    std::printf("ok   %s: refused (%s)\n", name.c_str(), err.c_str());
                }
                continue;
            }
            if (!err.empty()) {
                std::printf("BAD  %s t%d: %s\n", name.c_str(), threads, err.c_str());
                failed++;
            } else if (info.width != want_w || info.height != want_h) {
                std::printf("BAD  %s t%d: %dx%d, expected %dx%d\n", name.c_str(),
                            threads, info.width, info.height, want_w, want_h);
                failed++;
            } else if (!same(got, want)) {
                size_t at = got.size();
                float a = 0.0f, b = 0.0f;
                for (size_t i = 0; i < got.size() && i < want.size(); i++)
                    if (std::memcmp(&got[i], &want[i], 4) != 0) {
                        at = i; a = got[i]; b = want[i];
                        break;
                    }
                std::printf("BAD  %s t%d: differs at %zu (%g vs %g)\n", name.c_str(),
                            threads, at, a, b);
                failed++;
            } else if (threads == 1) {
                std::printf("ok   %s: %dx%d %s\n", name.c_str(), info.width, info.height,
                            info.compression.c_str());
            }
        }
    }
    std::printf("\n%d cases, %d failures\n", total, failed);
    return failed ? 1 : 0;
}
