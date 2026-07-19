// Backend parity tool for the dataset GT warp + byte->float conversion
// launch APIs (csrc/PixelWise.cuh launch_warp_* + EngineInternal.h raw
// converters). The SAME source builds under both backends:
//
//   CUDA build:   ./warp_parity dump ref.bin
//   Vulkan build: ./warp_parity compare ref.bin   (per device)
//
// Ref format: [nf tight floats].
//
// Everything here is deterministic per-pixel math (no atomics), so all
// outputs live in one tight channel. Both backends run the same canonical
// slang projection math (CUDA via the 2026.2.1 emission, Vulkan via
// SPIR-V), so residual differences are fast-math rounding plus isolated
// boundary pixels: a proj_nav valid-flip or a nearest/floor rounding flip
// moves a whole sample, so the pass criterion is a small violation
// fraction rather than a max-error bound. Mask bytes (0/1) are compared
// through the same channel; a flipped boundary texel counts as one
// violation.
//
// Coverage: all three wide camera models (PINHOLE with null dist_coeffs ->
// zeros-fallback path), u8/u16 image, u8 mask, u16-ray + f32-linear depth
// (GT at a different resolution than the intrinsics reference, exercising
// the sx/sy rescale), u8/f32 normal (with all-zero "no data" sentinel
// pixels), the equirectangular variants (x-wrap sampling), and the four
// raw byte->float converters with odd element counts (word-tail handling).

#include <PixelWise.cuh>
#include <EngineInternal.h>
#include <Tensor.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <vector>

using backend::MemcpyKind;

static std::vector<float> g_tight;

template <typename T>
T* upload(const std::vector<T>& host) {
    T* d = (T*)backend::device_malloc(
        ((host.size() * sizeof(T) + 3) / 4) * 4);
    backend::memcpy_sync(d, host.data(), host.size() * sizeof(T),
                         MemcpyKind::HostToDevice);
    return d;
}

template <typename T>
T* alloc_out(int64_t n) {
    int64_t bytes = ((n * (int64_t)sizeof(T) + 3) / 4) * 4;
    T* d = (T*)backend::device_malloc(bytes);
    backend::memset_sync(d, 0, bytes);
    return d;
}

void readback_f(const float* d, int64_t n) {
    size_t off = g_tight.size();
    g_tight.resize(off + n);
    backend::memcpy_sync(g_tight.data() + off, d, n * sizeof(float),
                         MemcpyKind::DeviceToHost);
}

void readback_b(const uint8_t* d, int64_t n) {
    std::vector<uint8_t> h(n);
    backend::memcpy_sync(h.data(), d, n, MemcpyKind::DeviceToHost);
    for (uint8_t v : h) g_tight.push_back((float)v);
}

struct Rng {
    std::mt19937 rng;
    explicit Rng(uint32_t seed) : rng(seed) {}
    float uf(float lo, float hi) {
        return lo + (hi - lo) * (float)(rng() & 0xffffff) / 16777215.0f;
    }
    std::vector<uint8_t> bytes(int64_t n) {
        std::vector<uint8_t> v(n);
        for (auto& x : v) x = (uint8_t)(rng() & 0xff);
        return v;
    }
    std::vector<uint16_t> words(int64_t n, uint16_t max) {
        std::vector<uint16_t> v(n);
        for (auto& x : v) x = (uint16_t)(rng() % (uint32_t)(max + 1));
        return v;
    }
    std::vector<float> vec(int64_t n, float lo, float hi) {
        std::vector<float> v(n);
        for (auto& x : v) x = uf(lo, hi);
        return v;
    }
};

struct V3 {
    float x, y, z;
};
static V3 norm3(V3 a) {
    float l = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return {a.x / l, a.y / l, a.z / l};
}
static V3 cross3(V3 a, V3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

// [K, 3, 3] face-axes table: orthonormal frames with axis_x/axis_y scaled
// by the half-fov tangent (raydir = axis_z + tx*axis_x + ty*axis_y).
static std::vector<float> make_axes(int K, float s) {
    const V3 zs[4] = {{0, 0, 1},
                      {1, 0, 0},
                      {0.577f, 0.577f, 0.577f},
                      {0, -1, 0}};
    std::vector<float> axes;
    for (int k = 0; k < K; k++) {
        V3 z = norm3(zs[k % 4]);
        V3 up = std::fabs(z.y) > 0.9f ? V3{0, 0, 1} : V3{0, 1, 0};
        V3 x = norm3(cross3(up, z));
        V3 y = cross3(z, x);
        float row[9] = {s * x.x, s * x.y, s * x.z, s * y.x, s * y.y,
                        s * y.z, z.x,     z.y,     z.z};
        axes.insert(axes.end(), row, row + 9);
    }
    return axes;
}

int main(int argc, char** argv) {
    if (argc != 3 ||
        (std::strcmp(argv[1], "dump") && std::strcmp(argv[1], "compare"))) {
        std::fprintf(stderr, "usage: %s dump|compare <ref.bin>\n", argv[0]);
        return 2;
    }
    const bool dumping = std::strcmp(argv[1], "dump") == 0;

    Rng r(260722u);

    // ---- raw byte -> float converters (odd totals exercise word tails) ----
    {
        auto u8img = r.bytes(2 * 13 * 17 * 3);
        float* out = alloc_out<float>(u8img.size());
        uint8_image_to_float_raw(upload(u8img), out, 2, 13, 17, 3);
        readback_f(out, u8img.size());

        auto u16img = r.words(1 * 9 * 21, 65535);
        float* out2 = alloc_out<float>(u16img.size());
        uint16_image_to_float_raw(upload(u16img), out2, 1, 9, 21, 1);
        readback_f(out2, u16img.size());

        auto u8nrm = r.bytes(2 * 11 * 15 * 3);
        for (size_t i = 0; i < u8nrm.size(); i += 97) u8nrm[i] = 0;
        float* out3 = alloc_out<float>(u8nrm.size());
        uint8_normal_to_float_raw(upload(u8nrm), out3, 2, 11, 15, 3);
        readback_f(out3, u8nrm.size());

        auto u16dep = r.words(1 * 14 * 10, 5000);
        float* out4 = alloc_out<float>(u16dep.size());
        uint16_depth_to_float_raw(upload(u16dep), out4, 1, 14, 10, 1);
        readback_f(out4, u16dep.size());
    }

    // ---- wide warps (FISHEYE / EQUISOLID / PINHOLE) -----------------------
    {
        const int B = 2, Hin = 48, Win = 64, K = 3, Hout = 24, Wout = 32;
        const int Hd = 36, Wd = 48;  // GT depth/normal at their own res
        float* d_axes = upload(make_axes(K, 0.7f));

        struct Cam {
            const char* model;
            float fx, fy;
            bool null_dist;
        };
        const Cam cams[3] = {
            {"FISHEYE", 20.0f, 20.5f, false},
            {"EQUISOLID", 24.0f, 24.0f, false},
            {"PINHOLE", 30.0f, 31.0f, true},
        };
        for (const Cam& cam : cams) {
            std::vector<float> intr;
            for (int b = 0; b < B; b++) {
                intr.push_back(cam.fx + 0.5f * b);
                intr.push_back(cam.fy + 0.3f * b);
                intr.push_back(0.5f * Win + 0.4f * b);
                intr.push_back(0.5f * Hin - 0.2f * b);
            }
            float* d_intr = upload(intr);
            float* d_dist = nullptr;
            if (!cam.null_dist) {
                std::vector<float> dist = r.vec(B * 10, -0.02f, 0.02f);
                d_dist = upload(dist);
            }
            const int64_t n_out = (int64_t)B * K * Hout * Wout;

            // u8 image (+ u16 for the first model)
            {
                auto img = r.bytes((int64_t)B * Hin * Win * 3);
                float* out = alloc_out<float>(n_out * 3);
                launch_warp_byte_to_float_wide(
                    cam.model, d_intr, d_dist, upload(img), false, B, Hin,
                    Win, 3, out, K, Hout, Wout, d_axes);
                readback_f(out, n_out * 3);
            }
            if (cam.model[0] == 'F') {
                auto img = r.words((int64_t)B * Hin * Win * 3, 65535);
                float* out = alloc_out<float>(n_out * 3);
                launch_warp_byte_to_float_wide(
                    cam.model, d_intr, d_dist, upload(img), true, B, Hin,
                    Win, 3, out, K, Hout, Wout, d_axes);
                readback_f(out, n_out * 3);
            }

            // mask (mostly-set with holes)
            {
                auto m = r.bytes((int64_t)B * Hin * Win);
                for (auto& v : m) v = v < 200 ? 1 : 0;
                uint8_t* out = alloc_out<uint8_t>(n_out);
                launch_warp_mask_wide(cam.model, d_intr, d_dist, upload(m),
                                      B, Hin, Win, out, K, Hout, Wout,
                                      d_axes);
                readback_b(out, n_out);
            }

            // u16 ray depth + f32 linear depth at (Hd, Wd)
            {
                auto dep = r.words((int64_t)B * Hd * Wd, 5000);
                for (size_t i = 0; i < dep.size(); i += 53) dep[i] = 0;
                float* out = alloc_out<float>(n_out);
                launch_warp_depth_wide(cam.model, d_intr, d_dist,
                                       upload(dep), 2, B, Hd, Wd, Hin, Win,
                                       out, K, Hout, Wout, d_axes, true);
                readback_f(out, n_out);
            }
            {
                auto dep = r.vec((int64_t)B * Hd * Wd, 0.5f, 8.0f);
                for (size_t i = 0; i < dep.size(); i += 41) dep[i] = -0.5f;
                float* out = alloc_out<float>(n_out);
                launch_warp_depth_wide(cam.model, d_intr, d_dist,
                                       upload(dep), 4, B, Hd, Wd, Hin, Win,
                                       out, K, Hout, Wout, d_axes, false);
                readback_f(out, n_out);
            }

            // u8 normal at (Hd, Wd) with "no data" holes (+ f32 for the
            // first model)
            {
                auto nrm = r.bytes((int64_t)B * Hd * Wd * 3);
                for (size_t i = 0; i + 2 < nrm.size(); i += 87)
                    nrm[i] = nrm[i + 1] = nrm[i + 2] = 0;
                float* out = alloc_out<float>(n_out * 3);
                launch_warp_normal_wide(cam.model, d_intr, d_dist,
                                        upload(nrm), 1, B, Hd, Wd, Hin, Win,
                                        out, K, Hout, Wout, d_axes);
                readback_f(out, n_out * 3);
            }
            if (cam.model[0] == 'F') {
                auto nrm = r.vec((int64_t)B * Hd * Wd * 3, -1.0f, 1.0f);
                float* out = alloc_out<float>(n_out * 3);
                launch_warp_normal_wide(cam.model, d_intr, d_dist,
                                        upload(nrm), 4, B, Hd, Wd, Hin, Win,
                                        out, K, Hout, Wout, d_axes);
                readback_f(out, n_out * 3);
            }
        }
    }

    // ---- equirectangular warps -------------------------------------------
    {
        const int B = 1, Hin = 32, Win = 64, K = 4, Hout = 16, Wout = 16;
        float* d_axes = upload(make_axes(K, 0.8f));
        const int64_t n_out = (int64_t)B * K * Hout * Wout;

        {
            auto img = r.bytes((int64_t)B * Hin * Win * 3);
            float* out = alloc_out<float>(n_out * 3);
            launch_warp_byte_to_float_equi(upload(img), false, B, Hin, Win,
                                           3, out, K, Hout, Wout, d_axes);
            readback_f(out, n_out * 3);
        }
        {
            auto img = r.words((int64_t)B * Hin * Win * 3, 65535);
            float* out = alloc_out<float>(n_out * 3);
            launch_warp_byte_to_float_equi(upload(img), true, B, Hin, Win,
                                           3, out, K, Hout, Wout, d_axes);
            readback_f(out, n_out * 3);
        }
        {
            auto m = r.bytes((int64_t)B * Hin * Win);
            for (auto& v : m) v = v < 180 ? 1 : 0;
            uint8_t* out = alloc_out<uint8_t>(n_out);
            launch_warp_mask_equi(upload(m), B, Hin, Win, out, K, Hout,
                                  Wout, d_axes);
            readback_b(out, n_out);
        }
        {
            auto dep = r.words((int64_t)B * Hin * Win, 5000);
            for (size_t i = 0; i < dep.size(); i += 53) dep[i] = 0;
            float* out = alloc_out<float>(n_out);
            launch_warp_depth_equi(upload(dep), 2, B, Hin, Win, out, K,
                                   Hout, Wout, d_axes, true);
            readback_f(out, n_out);
        }
        {
            auto dep = r.vec((int64_t)B * Hin * Win, 0.5f, 8.0f);
            for (size_t i = 0; i < dep.size(); i += 41) dep[i] = -0.5f;
            float* out = alloc_out<float>(n_out);
            launch_warp_depth_equi(upload(dep), 4, B, Hin, Win, out, K,
                                   Hout, Wout, d_axes, false);
            readback_f(out, n_out);
        }
        {
            auto nrm = r.bytes((int64_t)B * Hin * Win * 3);
            for (size_t i = 0; i + 2 < nrm.size(); i += 87)
                nrm[i] = nrm[i + 1] = nrm[i + 2] = 0;
            float* out = alloc_out<float>(n_out * 3);
            launch_warp_normal_equi(upload(nrm), 1, B, Hin, Win, out, K,
                                    Hout, Wout, d_axes);
            readback_f(out, n_out * 3);
        }
        {
            auto nrm = r.vec((int64_t)B * Hin * Win * 3, -1.0f, 1.0f);
            float* out = alloc_out<float>(n_out * 3);
            launch_warp_normal_equi(upload(nrm), 4, B, Hin, Win, out, K,
                                    Hout, Wout, d_axes);
            readback_f(out, n_out * 3);
        }
    }

    backend::device_synchronize();
    if (const char* err = backend::last_error()) {
        std::fprintf(stderr, "backend error: %s\n", err);
        return 1;
    }

    const int64_t nf = (int64_t)g_tight.size();
    if (dumping) {
        std::ofstream f(argv[2], std::ios::binary);
        f.write((const char*)&nf, 8);
        f.write((const char*)g_tight.data(), nf * 4);
        std::printf("warp_parity: dumped %lld floats to %s\n",
                    (long long)nf, argv[2]);
        return 0;
    }

    std::ifstream f(argv[2], std::ios::binary);
    if (!f) {
        std::fprintf(stderr, "cannot open %s\n", argv[2]);
        return 2;
    }
    int64_t rn = 0;
    f.read((char*)&rn, 8);
    if (rn != nf) {
        std::fprintf(stderr, "float count mismatch: ref %lld vs got %lld\n",
                     (long long)rn, (long long)nf);
        return 1;
    }
    std::vector<float> ref(nf);
    f.read((char*)ref.data(), nf * 4);

    int64_t viol = 0;
    double max_abs = 0;
    for (int64_t i = 0; i < nf; i++) {
        double d = std::fabs((double)g_tight[i] - (double)ref[i]);
        double tol = 5e-3 + 1e-3 * std::fabs((double)ref[i]);
        max_abs = std::max(max_abs, d);
        if (d > tol) viol++;
    }
    double frac = nf ? (double)viol / (double)nf : 0.0;
    std::printf("warp_parity: %lld floats (max_abs %.3g, violations %lld = "
                "%.5f%%)\n",
                (long long)nf, max_abs, (long long)viol, 100.0 * frac);
    bool pass = frac <= 2e-3;
    std::printf(pass ? "warp_parity: PASSED\n" : "warp_parity: FAILED\n");
    return pass ? 0 : 1;
}
