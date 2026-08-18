// Self-checking test for the GT-side bilateral grids' two invariants: one
// depth scalar per camera in the batch, and the "no ground truth here"
// sentinels (depth 0, normal black) surviving the transform unchanged. Both
// failures are silent -- the loss just stops seeing (or starts inventing)
// supervision. Runs under either backend:
//
//   ./gt_bilagrid_sentinels

#include <kernels/bilagrid/BilagridBindings.h>
#include <core/Tensor.h>

#include <cmath>
#include <cstdio>
#include <vector>

using backend::MemcpyKind;

template <bool invert_quantile>
int batch_quantile_masked_radix_select(
    const float* d_x, int B, int N, float q, float* d_out, uint32_t* temp,
    backend::Stream stream);

static int g_failures = 0;

static void check(bool ok, const char* what) {
    std::printf("%-58s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) g_failures++;
}

template <typename T>
static T* upload(const std::vector<T>& host) {
    T* d = (T*)backend::device_malloc(((host.size() * sizeof(T) + 3) / 4) * 4);
    backend::memcpy_sync(d, host.data(), host.size() * sizeof(T),
                         MemcpyKind::HostToDevice);
    return d;
}

static std::vector<float> download(const float* d, int64_t n) {
    std::vector<float> h((size_t)n);
    backend::memcpy_sync(h.data(), d, (size_t)n * sizeof(float),
                         MemcpyKind::DeviceToHost);
    return h;
}

// EngineBilagrid.cpp's depth branch wants a median per camera over [C, h, w];
// the uniform sampler then reads scalars[cam]. `patched` fills slot 0 alone,
// so a batch of cameras must not ask for it.
static void test_depth_scalars() {
    const int C = 6, h = 32, w = 32, n = h * w;
    std::vector<float> depth((size_t)C * n);
    for (int c = 0; c < C; c++)
        for (int i = 0; i < n; i++)
            depth[(size_t)c * n + i] = 1.0f + (float)c + (float)(i % 7);

    float* d_depth = upload(depth);
    float* d_out = upload(std::vector<float>((size_t)C, -1.0f));
    uint32_t* temp = (uint32_t*)backend::device_malloc(
        (size_t)(256 + 5) * C * sizeof(uint32_t));
    backend::memset_sync(temp, 0, (size_t)(256 + 5) * C * sizeof(uint32_t));

    TorchTensorView depth_tv((uint64_t)d_depth, 4,
                             {(int64_t)C, (int64_t)h, (int64_t)w, 1LL, 1LL});
    TorchTensorView out_tv((uint64_t)d_out, 4, {(int64_t)C});

    compute_depth_scalars_tensor(depth_tv, /*patched=*/true, out_tv);
    backend::device_synchronize();
    std::vector<float> one = download(d_out, C);
    bool only_first = one[0] > 0.0f;
    for (int c = 1; c < C; c++) only_first &= (one[(size_t)c] == -1.0f);
    check(only_first, "depth scalars: `patched` writes slot 0 and nothing else");

    compute_depth_scalars_tensor(depth_tv, /*patched=*/false, out_tv);
    backend::device_synchronize();
    std::vector<float> got = download(d_out, C);
    bool all_written = true, distinct = false;
    for (int c = 0; c < C; c++) {
        if (!(got[(size_t)c] > 0.0f)) all_written = false;
        if (c > 0 && got[(size_t)c] != got[0]) distinct = true;
    }
    check(all_written, "depth scalars: every camera in the batch is written");
    check(distinct, "depth scalars: per-camera, not one median for the batch");
}

// A zeroed grid is the identity (exp(exp(0)*log(d)) == d), so anything but a
// round trip here is the transform mangling a value it was handed.
static void test_depth_sentinel() {
    const int C = 2, L = 4, H = 4, W = 4, h = 8, w = 8, n = h * w;
    std::vector<float> grid((size_t)C * L * H * W * 2, 0.0f);
    std::vector<float> depth((size_t)C * n);
    for (int c = 0; c < C; c++)
        for (int i = 0; i < n; i++)
            depth[(size_t)c * n + i] = (i % 5 == 0) ? 0.0f : 0.5f + 0.01f * i;

    float* d_grid = upload(grid);
    float* d_depth = upload(depth);
    float* d_scalars = upload(std::vector<float>((size_t)C, 1.0f));
    float* d_out = upload(std::vector<float>((size_t)C * n, -1.0f));

    bilagrid_depth_uniform_sample_forward(BilagridReader(d_grid), d_depth,
                                          d_scalars, d_out, C, L, H, W, h, w,
                                          backend::kDefaultStream, nullptr);
    backend::device_synchronize();

    std::vector<float> got = download(d_out, (int64_t)C * n);
    bool zeros_kept = true, values_kept = true;
    for (size_t i = 0; i < got.size(); i++) {
        if (depth[i] == 0.0f) zeros_kept &= (got[i] == 0.0f);
        else values_kept &= std::fabs(got[i] - depth[i]) < 1e-3f * depth[i];
    }
    check(zeros_kept, "depth bilagrid: 0 stays 0");
    check(values_kept, "depth bilagrid: identity grid round-trips a depth");
}

static void test_normal_sentinel() {
    const int C = 1, L = 4, H = 4, W = 4, h = 4, w = 4, n = h * w;
    std::vector<float> grid((size_t)C * L * H * W * 3, 0.0f);
    std::vector<float> normal((size_t)C * n * 3);
    // Black, a float map's 0, and mid-grey through byte/127.5 - 1 are all "no
    // normal here"; only the fourth is one.
    for (int i = 0; i < n; i++) {
        float* v = &normal[(size_t)i * 3];
        if (i % 4 == 0) { v[0] = v[1] = v[2] = -1.0f; }
        else if (i % 4 == 1) { v[0] = v[1] = v[2] = 0.0f; }
        else if (i % 4 == 2) { v[0] = v[1] = v[2] = 128.0f / 127.5f - 1.0f; }
        else { v[0] = 0.0f; v[1] = 0.0f; v[2] = 1.0f; }
    }

    float* d_grid = upload(grid);
    float* d_normal = upload(normal);
    float* d_out = upload(std::vector<float>((size_t)C * n * 3, -7.0f));

    bilagrid_normal_uniform_sample_forward(BilagridReader(d_grid), d_normal,
                                           d_out, C, L, H, W, h, w,
                                           backend::kDefaultStream, nullptr);
    backend::device_synchronize();

    std::vector<float> got = download(d_out, (int64_t)C * n * 3);
    bool kept = true, rotated = true;
    for (int i = 0; i < n; i++) {
        const float* v = &normal[(size_t)i * 3];
        const float* o = &got[(size_t)i * 3];
        const bool valid = v[0] + v[1] + v[2] > -2.366f &&
                           v[0]*v[0] + v[1]*v[1] + v[2]*v[2] > 0.25f;
        if (!valid)
            kept &= (o[0] == v[0] && o[1] == v[1] && o[2] == v[2]);
        else
            rotated &= std::fabs(o[2] - 1.0f) < 1e-5f;
    }
    check(kept, "normal bilagrid: the no-data sentinel is passed through");
    check(rotated, "normal bilagrid: identity grid round-trips a normal");
}

int main() {
    test_depth_scalars();
    test_depth_sentinel();
    test_normal_sentinel();
    std::printf("%s\n", g_failures ? "FAILED" : "all ok");
    return g_failures ? 1 : 0;
}
