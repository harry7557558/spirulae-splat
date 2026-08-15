// Progress.cpp -- see Progress.h.

#include "sfm/core/Progress.h"

#include "sfm/core/Model.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <vector>

namespace fs = std::filesystem;

namespace sfm {
namespace progress {

namespace {

using Clock = std::chrono::steady_clock;

// Long enough that a fast mapper does not spend its time serializing, short
// enough that a screen looks live. Registration on a small capture is
// milliseconds apart, so without this the file would be rewritten thousands of
// times a second.
constexpr double kInterval = 1.5;

struct State {
    std::mutex mu;
    std::string dir;
    Clock::time_point model_at{};
    Clock::time_point pairs_at{};
    bool model_started = false, pairs_started = false;

    // The pair matrix, binned down to kMatrixBins per side.
    uint32_t n_images = 0, bins = 0;
    std::vector<uint32_t> counts, planned, verified;
    bool pairs_dirty = false;
};

State& state() {
    static State s;
    return s;
}

// Whole file, then rename: a reader polling the directory either sees the
// previous snapshot or this one, never a prefix of one.
void write_atomic(const std::string& name, const std::string& bytes) {
    State& s = state();
    const fs::path dst = fs::path(s.dir) / name;
    const fs::path tmp = fs::path(s.dir) / (name + ".tmp");
    {
        std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
        if (!f) return;
        f.write(bytes.data(), (std::streamsize)bytes.size());
        if (!f) return;
    }
    std::error_code ec;
    fs::rename(tmp, dst, ec);
    if (ec) fs::remove(tmp, ec);
}

void put(std::string& b, const void* p, size_t n) {
    b.append((const char*)p, n);
}
void put_u32(std::string& b, uint32_t v) { put(b, &v, 4); }
void put_u64(std::string& b, uint64_t v) { put(b, &v, 8); }
void put_f32(std::string& b, float v) { put(b, &v, 4); }

// True when `last` is far enough behind now, and stamps it if so.
bool due(Clock::time_point& last, bool& started) {
    const auto now = Clock::now();
    if (started &&
        std::chrono::duration<double>(now - last).count() < kInterval)
        return false;
    started = true;
    last = now;
    return true;
}

void write_pairs_locked() {
    State& s = state();
    if (!s.pairs_dirty || s.counts.empty()) return;
    std::string b;
    b.reserve(16 + s.counts.size() * 12);
    put(b, "VKPP", 4);
    put_u32(b, 2);
    put_u32(b, s.n_images);
    put_u32(b, s.bins);
    put(b, s.counts.data(), s.counts.size() * 4);
    put(b, s.planned.data(), s.planned.size() * 4);
    put(b, s.verified.data(), s.verified.size() * 4);
    write_atomic("pairs.bin", b);
    s.pairs_dirty = false;
}

// Which cell of the matrix a pair of images falls in.
size_t cell_of(const State& s, uint32_t image1, uint32_t image2, size_t& mirror) {
    const uint32_t a = (uint32_t)((uint64_t)image1 * s.bins / s.n_images);
    const uint32_t b = (uint32_t)((uint64_t)image2 * s.bins / s.n_images);
    mirror = (size_t)b * s.bins + a;
    return (size_t)a * s.bins + b;
}

}  // namespace

void set_dir(const std::string& dir) {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    s.dir = dir;
    if (dir.empty()) return;
    std::error_code ec;
    fs::create_directories(dir, ec);
}

bool enabled() {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    return !s.dir.empty();
}

void model(const Reconstruction& rec, bool force, const PointColor& color) {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    if (s.dir.empty()) return;
    if (!due(s.model_at, s.model_started) && !force) return;

    // Registered images only: an unregistered one has no pose to draw.
    std::vector<const Image*> imgs;
    for (const auto& kv : rec.images)
        if (kv.second.registered) imgs.push_back(&kv.second);

    const uint64_t n_pts = rec.points3D.size();
    const uint64_t stride = n_pts > kMaxPoints ? (n_pts / kMaxPoints + 1) : 1;

    std::string b;
    b.reserve(64 + imgs.size() * 128 + (size_t)(n_pts / stride + 1) * 15);
    put(b, "VKPM", 4);
    put_u32(b, 2);
    put_u32(b, (uint32_t)rec.images.size());
    put_u32(b, (uint32_t)imgs.size());
    put_u64(b, n_pts);
    for (const Image* im : imgs) {
        // COLMAP world->camera (R, t) to nerfstudio/OpenGL camera->world:
        // R^T with columns 1 and 2 negated, translation -R^T t. Same
        // conversion ColmapParser does, done here so the reader needs none.
        const Mat3& R = im->pose.R;      // row-major, world -> camera
        const Vec3& t = im->pose.t;
        const double C[3] = {
            -(R[0] * t.x + R[3] * t.y + R[6] * t.z),
            -(R[1] * t.x + R[4] * t.y + R[7] * t.z),
            -(R[2] * t.x + R[5] * t.y + R[8] * t.z)};
        for (int r = 0; r < 3; r++) {
            put_f32(b, (float)R[r]);            // R^T row r, column 0
            put_f32(b, (float)-R[3 + r]);       // ... column 1, negated
            put_f32(b, (float)-R[6 + r]);       // ... column 2, negated
            put_f32(b, (float)C[r]);
        }
        // The camera as cameras.bin would hold it: a COLMAP model id and its
        // parameters. The reader hands them to the dataset parser's own
        // mapping rather than keeping a second copy of it, which is what makes
        // a fisheye frustum draw as a fisheye.
        const auto cam = rec.cameras.find(im->camera_id);
        const Camera c = cam == rec.cameras.end() ? Camera{} : cam->second;
        put_u32(b, (uint32_t)c.width);
        put_u32(b, (uint32_t)c.height);
        put_u32(b, (uint32_t)camColmapId(c.model));
        double ps[12] = {};
        packColmap(c, ps);
        const uint32_t np = (uint32_t)camColmapParams(c.model);
        put_u32(b, np);
        for (uint32_t k = 0; k < np; k++) put(b, &ps[k], 8);
    }

    // Count first, so the reader can size its buffers before the loop.
    uint32_t written = 0;
    for (uint64_t i = 0; i < n_pts; i += stride) written++;
    put_u32(b, written);
    uint64_t k = 0, next = 0;
    for (const auto& kv : rec.points3D) {
        if (k++ != next) continue;
        next += stride;
        put_f32(b, (float)kv.second.xyz.x);
        put_f32(b, (float)kv.second.xyz.y);
        put_f32(b, (float)kv.second.xyz.z);
        uint8_t rgb[3] = {kv.second.rgb[0], kv.second.rgb[1], kv.second.rgb[2]};
        if (color) color(kv.second, rgb);
        put(b, rgb, 3);
    }
    write_atomic("model.bin", b);
}

void begin_matching(uint32_t n_images,
                    const std::vector<std::pair<uint32_t, uint32_t>>& pairs) {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    if (s.dir.empty() || n_images == 0) return;
    s.n_images = n_images;
    s.bins = n_images < kMatrixBins ? n_images : kMatrixBins;
    s.counts.assign((size_t)s.bins * s.bins, 0);
    s.planned.assign((size_t)s.bins * s.bins, 0);
    s.verified.assign((size_t)s.bins * s.bins, 0);
    for (const auto& p : pairs) {
        if (p.first >= n_images || p.second >= n_images) continue;
        size_t mirror = 0;
        const size_t c = cell_of(s, p.first, p.second, mirror);
        s.planned[c]++;
        if (mirror != c) s.planned[mirror]++;
    }
    s.pairs_started = false;
    s.pairs_dirty = true;
}

void pair(uint32_t image1, uint32_t image2, uint32_t inliers) {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    if (s.dir.empty() || s.counts.empty()) return;
    if (image1 >= s.n_images || image2 >= s.n_images) return;
    size_t mirror = 0;
    const size_t c = cell_of(s, image1, image2, mirror);
    s.counts[c] += inliers;
    s.verified[c]++;
    if (mirror != c) {
        s.counts[mirror] += inliers;
        s.verified[mirror]++;
    }
    s.pairs_dirty = true;
    if (due(s.pairs_at, s.pairs_started)) write_pairs_locked();
}

void flush() {
    State& s = state();
    std::lock_guard<std::mutex> lk(s.mu);
    if (s.dir.empty()) return;
    write_pairs_locked();
}

}  // namespace progress
}  // namespace sfm
