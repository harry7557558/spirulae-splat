#include "app/GeometryWarp.h"

#include "data/CameraMath.h"
#include "data/SourceCamera.h"
#include "nn/core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace app {
namespace {

constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

// Face half-width as a multiple of the 45 degrees it owns: 1.25 spans 51.3,
// leaving neighbours a 6.3-degree strip to cross-fade over.
constexpr double kOverlap = 1.25;

// Tilt of the fan's four outer faces. At 1.27 rad the five reach 124 degrees
// off axis, which is a 220-degree lens with the overlap still intact.
constexpr double kFanTilt = 1.27;

// [K][3][3]: two in-plane axes whose LENGTH is the face's half-extent, then
// the unit optical axis. az + tx*ax + ty*ay over tx, ty in -1..1 is its ray.
struct AxisTable {
    double a[6][3][3];
};

AxisTable scaled(const double (*src)[3][3], int k) {
    AxisTable t{};
    for (int i = 0; i < k; ++i)
        for (int r = 0; r < 3; ++r)
            for (int m = 0; m < 3; ++m)
                t.a[i][r][m] = src[i][r][m] * (r < 2 ? kOverlap : 1.0);
    return t;
}

const double* face_axes(int k) {
    static const AxisTable fan = [] {
        const double s = std::sin(kFanTilt), c = std::cos(kFanTilt);
        const double base[5][3][3] = {
            {{ 1, 0, 0}, { 0,  1, 0}, { 0,  0, 1}},
            {{ 0, 1, 0}, {-c,  0, s}, { s,  0, c}},
            {{-1, 0, 0}, { 0, -c, s}, { 0,  s, c}},
            {{ 0,-1, 0}, { c,  0, s}, {-s,  0, c}},
            {{ 1, 0, 0}, { 0,  c, s}, { 0, -s, c}},
        };
        return scaled(base, 5);
    }();
    static const AxisTable cube = [] {
        const double base[6][3][3] = {
            {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
            {{ 0, 0,-1}, { 0, 1, 0}, { 1, 0, 0}},
            {{-1, 0, 0}, { 0, 0, 1}, { 0, 1, 0}},
            {{ 0,-1, 0}, { 0, 0, 1}, {-1, 0, 0}},
            {{ 1, 0, 0}, { 0, 0, 1}, { 0,-1, 0}},
            {{ 1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}},
        };
        return scaled(base, 6);
    }();
    if (k == 5) return &fan.a[0][0][0];
    if (k == 6) return &cube.a[0][0][0];
    return nullptr;
}

double dot3(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double len3(const double* a) { return std::sqrt(dot3(a, a)); }

void cross3(const double* a, const double* b, double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

// Solve az + tx*ax + ty*ay = t*r by Cramer. t <= 0 is the direction behind
// the camera, which meets the same plane at the same point.
bool face_coords(const double* a, const double r[3], double& tx, double& ty) {
    double c[3];
    cross3(a, a + 3, c);                       // ax x ay
    const double det = dot3(c, r);
    if (std::fabs(det) < 1e-12) return false;
    const double inv = 1.0 / det;
    if (dot3(c, a + 6) * inv <= 0.0) return false;
    double c1[3], c2[3];
    cross3(r, a + 3, c1);                      // r x ay
    cross3(a, r, c2);                          // ax x r
    tx = dot3(c1, a + 6) * inv;
    ty = dot3(c2, a + 6) * inv;
    return std::fabs(tx) < 1.0 && std::fabs(ty) < 1.0;
}

// 1 over the face's own 45 degrees, ramping to 0 at the edge it reaches past.
// Unit axes therefore weight the whole face 1 and blend nothing.
float face_weight(const double* a, double tx, double ty) {
    const double az = len3(a + 6);
    const double wx = az / len3(a), wy = az / len3(a + 3);
    double w = 1.0;
    if (wx < 1.0 && wy < 1.0)
        w = std::fmin(std::fmin((1.0 - std::fabs(tx)) / (1.0 - wx),
                                (1.0 - std::fabs(ty)) / (1.0 - wy)), 1.0);
    if (!(w > 0.0)) return 0.0f;
    return (float)(w * w * (3.0 - 2.0 * w));   // smoothstep: C1 across the seam
}

// Bilinear sample of a [h, w, C] map at a pixel coordinate whose 0.5 is the
// first pixel's centre, edges clamped.
void sample(const float* img, int w, int h, int C, float px, float py, float* out) {
    const float x = px - 0.5f, y = py - 0.5f;
    int x0 = (int)std::floor(x), y0 = (int)std::floor(y);
    const float fx = x - x0, fy = y - y0;
    const int x1 = std::min(std::max(x0 + 1, 0), w - 1);
    const int y1 = std::min(std::max(y0 + 1, 0), h - 1);
    x0 = std::min(std::max(x0, 0), w - 1);
    y0 = std::min(std::max(y0, 0), h - 1);
    const float* p00 = img + ((size_t)y0 * w + x0) * C;
    const float* p01 = img + ((size_t)y0 * w + x1) * C;
    const float* p10 = img + ((size_t)y1 * w + x0) * C;
    const float* p11 = img + ((size_t)y1 * w + x1) * C;
    for (int c = 0; c < C; ++c)
        out[c] = (p00[c] * (1 - fx) + p01[c] * fx) * (1 - fy) +
                 (p10[c] * (1 - fx) + p11[c] * fx) * fy;
}

}  // namespace

std::vector<float> resize_area(const uint8_t* src, int sw, int sh, int channels,
                               int dw, int dh) {
    std::vector<float> out((size_t)dw * dh * 3, 0.0f);
    if (!src || sw <= 0 || sh <= 0 || dw <= 0 || dh <= 0 || channels < 3)
        return out;
    nn::parallel_for(dh, [&](int64_t y0, int64_t y1) {
        for (int64_t y = y0; y < y1; ++y) {
            const int sy0 = (int)((double)y * sh / dh);
            const int sy1 = std::max(sy0 + 1, (int)((double)(y + 1) * sh / dh));
            for (int x = 0; x < dw; ++x) {
                const int sx0 = (int)((double)x * sw / dw);
                const int sx1 = std::max(sx0 + 1, (int)((double)(x + 1) * sw / dw));
                double acc[3] = {0, 0, 0};
                for (int sy = sy0; sy < sy1; ++sy)
                    for (int sx = sx0; sx < sx1; ++sx) {
                        const uint8_t* p = &src[((size_t)sy * sw + sx) * channels];
                        for (int c = 0; c < 3; ++c) acc[c] += p[c];
                    }
                const double n = (double)(sy1 - sy0) * (sx1 - sx0) * 255.0;
                for (int c = 0; c < 3; ++c)
                    out[((size_t)y * dw + x) * 3 + c] = (float)(acc[c] / n);
            }
        }
    });
    return out;
}

void GeometryWarp::plan(const GeometryCamera& cam, int out_w, int out_h, bool split,
                        int patch, int max_face) {
    out_w_ = out_w;
    out_h_ = out_h;

    const bool equirect = cam.model == (int)CameraModelType::EQUIRECTANGULAR;
    faces_ = split ? (equirect ? 6 : 5) : 1;
    axes_ = split ? face_axes(faces_) : nullptr;

    // The intrinsics at the output resolution: what the inverse map indexes.
    const double sx = (double)out_w / cam.width, sy = (double)out_h / cam.height;
    const double fx = cam.fx * sx, fy = cam.fy * sy;
    const double cx = cam.cx * sx, cy = cam.cy * sy;

    // The forward map reads at the sampling resolution instead, so a split at
    // native face size resamples the frame once rather than through a smaller
    // copy of itself.
    const int natural =
        split ? (int)std::lround(std::sqrt((double)cam.width * cam.height / faces_)) : 0;
    const double fs = (split && max_face > 0 && natural > max_face)
                          ? (double)max_face / natural
                          : 1.0;
    sw_ = split ? std::max(1, (int)std::lround(cam.width * fs)) : out_w;
    sh_ = split ? std::max(1, (int)std::lround(cam.height * fs)) : out_h;
    const double gx = (double)sw_ / cam.width, gy = (double)sh_ / cam.height;

    // A camera the parser had to fit reads through its TRUE lens; only that
    // direction has a closed form, so the inverse below still goes through the
    // fit. The one approximate step here, and sub-pixel.
    const bool from_source = cam.source_model >= 0;
    auto ray_to_pixel = [&](const double r[3], double& px, double& py) -> bool {
        double uv[2];
        if (from_source) {
            double u, v;
            if (!srccam::project(cam.source_model, cam.source_params, r[0], r[1], r[2],
                                 &u, &v))
                return false;
            px = u * gx;
            py = v * gy;
            return true;
        }
        if (!camhost::project_ray(r, cam.model, cam.distortion, cam.dist, uv))
            return false;
        px = (uv[0] * cam.fx + cam.cx) * gx;
        py = (uv[1] * cam.fy + cam.cy) * gy;
        return true;
    };

    if (faces_ == 1) {
        fw_ = out_w;
        fh_ = out_h;
        face_focal_ = fx;
        face_focal_y_ = fy;
        face_cx_ = cx;
        face_cy_ = cy;
    } else {
        int side = std::min(natural, max_face > 0 ? max_face : natural);
        side = std::max(patch, side / patch * patch);
        fw_ = fh_ = side;
        // tx spans -1..1 over the side, and the in-plane axes are longer than 1.
        face_focal_ = face_focal_y_ = side * 0.5 / kOverlap;
        face_cx_ = face_cy_ = side * 0.5;
    }

    // ---- forward: where each face pixel reads from ------------------------
    to_src_.assign((size_t)faces_, {});
    face_valid_.assign((size_t)faces_, {});
    for (int k = 0; k < faces_; ++k) {
        std::vector<float>& map = to_src_[(size_t)k];
        std::vector<float>& ok = face_valid_[(size_t)k];
        map.assign((size_t)fw_ * fh_ * 2, kNaN);
        ok.assign((size_t)fw_ * fh_, 0.0f);
        nn::parallel_for(fh_, [&](int64_t y0, int64_t y1) {
            for (int64_t y = y0; y < y1; ++y)
                for (int x = 0; x < fw_; ++x) {
                    // One face keeps the camera's own axis and intrinsics; a
                    // split face is a rotation with tx, ty spanning -1..1.
                    double r[3];
                    if (faces_ == 1) {
                        r[0] = ((double)x + 0.5 - cx) / fx;
                        r[1] = ((double)y + 0.5 - cy) / fy;
                        r[2] = 1.0;
                    } else {
                        const double tx = -1.0 + 2.0 * ((double)x + 0.5) / fw_;
                        const double ty = -1.0 + 2.0 * ((double)y + 0.5) / fh_;
                        const double* a = axes_ + (size_t)k * 9;
                        for (int i = 0; i < 3; ++i)
                            r[i] = a[6 + i] + tx * a[i] + ty * a[3 + i];
                    }
                    double px = 0, py = 0;
                    if (!ray_to_pixel(r, px, py)) continue;
                    const size_t i = (size_t)y * fw_ + x;
                    map[i * 2 + 0] = (float)px;
                    map[i * 2 + 1] = (float)py;
                    if (px >= 0 && py >= 0 && px <= sw_ && py <= sh_) ok[i] = 1.0f;
                }
        });
    }

    // ---- inverse: which faces reach each output pixel -----------------------
    const size_t n = (size_t)out_w * out_h;
    // Staged a row at a time: a slot per face per pixel would be six times
    // what this keeps, nine tenths of it empty.
    std::vector<std::vector<Contrib>> rows((size_t)out_h);
    std::vector<int32_t> count(n, 0);
    nn::parallel_for(out_h, [&](int64_t y0, int64_t y1) {
        for (int64_t y = y0; y < y1; ++y) {
            std::vector<Contrib>& row = rows[(size_t)y];
            row.reserve((size_t)out_w);
            for (int x = 0; x < out_w; ++x) {
                const size_t i = (size_t)y * out_w + x;
                const double u = ((double)x + 0.5 - cx) / fx;
                const double v = ((double)y + 0.5 - cy) / fy;
                double r[3];
                if (!camhost::generate_ray(u, v, cam.model, cam.distortion, cam.dist, r))
                    continue;

                if (faces_ == 1) {
                    // The pinhole face shares the optical axis, so the normal
                    // needs no rotation and the depth is already this camera's
                    // z. Ray depth is the secant of the undistorted angle.
                    if (r[2] <= 1e-9) continue;
                    const double un = r[0] / r[2], vn = r[1] / r[2];
                    const double px = un * fx + cx, py = vn * fy + cy;
                    if (px < 0 || py < 0 || px > fw_ || py > fh_) continue;
                    Contrib c;
                    c.face = 0;
                    c.px = (float)px;
                    c.py = (float)py;
                    c.w = 1.0f;
                    c.sz = 1.0f;
                    c.sr = (float)std::sqrt(un * un + vn * vn + 1.0);
                    row.push_back(c);
                    count[i] = 1;
                    continue;
                }

                for (int k = 0; k < faces_; ++k) {
                    const double* a = axes_ + (size_t)k * 9;
                    double tx = 0, ty = 0;
                    if (!face_coords(a, r, tx, ty)) continue;
                    const float w = face_weight(a, tx, ty);
                    if (w <= 0.0f) continue;
                    Contrib c;
                    c.face = k;
                    // 0.5-centred, which is what `sample` takes.
                    c.px = (float)((tx + 1.0) * 0.5 * fw_);
                    c.py = (float)((ty + 1.0) * 0.5 * fh_);
                    c.w = w;
                    // A point at face z-depth d sits at d * this vector, so
                    // its z and its length are the two conversions.
                    double rr[3];
                    for (int m = 0; m < 3; ++m)
                        rr[m] = a[6 + m] + tx * a[m] + ty * a[3 + m];
                    c.sz = (float)rr[2];
                    c.sr = (float)len3(rr);
                    row.push_back(c);
                    ++count[i];
                }
            }
        }
    });

    contrib_off_.assign(n + 1, 0);
    for (size_t i = 0; i < n; ++i) contrib_off_[i + 1] = contrib_off_[i] + count[i];
    contrib_.resize((size_t)contrib_off_[n]);
    if (contrib_.empty())
        throw std::runtime_error("no part of this camera's frame reaches a pinhole "
                                 "face; check the camera model and its coefficients");
    // A row's pixels are consecutive, so its stage lands in one run.
    nn::parallel_for(out_h, [&](int64_t y0, int64_t y1) {
        for (int64_t y = y0; y < y1; ++y)
            std::copy(rows[(size_t)y].begin(), rows[(size_t)y].end(),
                      contrib_.begin() + contrib_off_[(size_t)y * out_w]);
    });
}

void GeometryWarp::sampleFace(int k, const float* src, std::vector<float>& dst) const {
    dst.assign((size_t)fw_ * fh_ * 3, 0.5f);
    const std::vector<float>& map = to_src_[(size_t)k];
    const std::vector<float>& ok = face_valid_[(size_t)k];
    nn::parallel_for(fh_, [&](int64_t y0, int64_t y1) {
        for (int64_t y = y0; y < y1; ++y)
            for (int x = 0; x < fw_; ++x) {
                const size_t i = (size_t)y * fw_ + x;
                if (ok[i] == 0.0f) continue;   // outside the frame stays grey
                sample(src, sw_, sh_, 3, map[i * 2 + 0], map[i * 2 + 1], &dst[i * 3]);
            }
    });
}

// Per-face log offsets reconciling the median log depth ratio over every
// pair's shared strip -- a graph Laplacian on faces() nodes. Median because
// the disagreements are worst where one face resolves an edge and one does not.
std::vector<double> GeometryWarp::alignFaces(
    const std::vector<std::vector<float>>& depth) const {
    std::vector<double> out((size_t)faces_, 0.0);
    if (faces_ < 2) return out;

    std::vector<std::vector<float>> pairs((size_t)faces_ * faces_);
    const size_t n = (size_t)out_w_ * out_h_;
    for (size_t i = 0; i < n; ++i) {
        if (contrib_off_[i + 1] - contrib_off_[i] < 2) continue;
        double lv[8];
        int    fk[8], m = 0;
        for (int64_t c = contrib_off_[i]; c < contrib_off_[i + 1] && m < 8; ++c) {
            const Contrib& t = contrib_[(size_t)c];
            if (t.w <= 0.0f || depth[(size_t)t.face].empty()) continue;
            float ok = 0, d = 0;
            // What a face predicted over its mid-grey fill measures nothing.
            sample(face_valid_[(size_t)t.face].data(), fw_, fh_, 1, t.px, t.py, &ok);
            if (ok < 0.999f) continue;
            sample(depth[(size_t)t.face].data(), fw_, fh_, 1, t.px, t.py, &d);
            if (!(d > 0.0f)) continue;
            // As a distance from the camera: the one quantity two faces
            // agree on the meaning of.
            lv[m] = std::log((double)d * t.sr);
            fk[m] = t.face;
            ++m;
        }
        for (int a = 0; a < m; ++a)
            for (int b = a + 1; b < m; ++b) {
                int p = fk[a], q = fk[b];
                double d = lv[b] - lv[a];
                if (p > q) { std::swap(p, q); d = -d; }
                pairs[(size_t)p * faces_ + q].push_back((float)d);
            }
    }

    // L s = rhs, where s[p] - s[q] should be the median log ratio.
    double L[6][6] = {}, rhs[6] = {};
    bool any = false;
    for (int p = 0; p < faces_; ++p)
        for (int q = p + 1; q < faces_; ++q) {
            std::vector<float>& v = pairs[(size_t)p * faces_ + q];
            // A handful of shared pixels is a corner, not a seam.
            if (v.size() < 256) continue;
            std::nth_element(v.begin(), v.begin() + (long)(v.size() / 2), v.end());
            const double d = v[v.size() / 2];
            L[p][p] += 1.0;  L[q][q] += 1.0;
            L[p][q] -= 1.0;  L[q][p] -= 1.0;
            rhs[p]  += d;    rhs[q]  -= d;
            any = true;
        }
    if (!any) return out;
    // Singular along the constant: the ridge picks the smallest solution and
    // the mean below leaves the frame's overall scale alone.
    for (int p = 0; p < faces_; ++p) L[p][p] += 1e-3;

    for (int p = 0; p < faces_; ++p) {
        int piv = p;
        for (int q = p + 1; q < faces_; ++q)
            if (std::fabs(L[q][p]) > std::fabs(L[piv][p])) piv = q;
        if (std::fabs(L[piv][p]) < 1e-12) return out;
        std::swap_ranges(L[p], L[p] + faces_, L[piv]);
        std::swap(rhs[p], rhs[piv]);
        for (int q = p + 1; q < faces_; ++q) {
            const double f = L[q][p] / L[p][p];
            for (int m = p; m < faces_; ++m) L[q][m] -= f * L[p][m];
            rhs[q] -= f * rhs[p];
        }
    }
    for (int p = faces_ - 1; p >= 0; --p) {
        double s = rhs[p];
        for (int m = p + 1; m < faces_; ++m) s -= L[p][m] * out[(size_t)m];
        out[(size_t)p] = s / L[p][p];
    }

    double mean = 0;
    for (double v : out) mean += v;
    mean /= faces_;
    // Past 4x is not a scale disagreement, it is a face that saw something
    // else; leave that one where the network put it.
    for (double& v : out) v = std::fmin(std::fmax(v - mean, -1.386), 1.386);
    return out;
}

void GeometryWarp::gather(const std::vector<std::vector<float>>& depth,
                          const std::vector<std::vector<float>>& normal, bool ray_depth,
                          std::vector<float>* out_depth,
                          std::vector<float>* out_normal) const {
    const size_t n = (size_t)out_w_ * out_h_;
    const bool do_depth = out_depth && !depth.empty();
    const bool do_normal = out_normal && !normal.empty();
    if (out_depth) out_depth->assign(n, 0.0f);
    if (out_normal) out_normal->assign(n * 3, 0.0f);
    if (!do_depth && !do_normal) return;

    const std::vector<double> scale =
        do_depth ? alignFaces(depth) : std::vector<double>((size_t)faces_, 0.0);

    // Unit frames: the axis lengths are the face's extent, not a rotation.
    std::vector<double> frame((size_t)faces_ * 9, 0.0);
    for (int k = 0; k < faces_; ++k)
        for (int r = 0; r < 3; ++r) {
            const double* a = axes_ ? axes_ + (size_t)k * 9 + r * 3 : nullptr;
            double* f = &frame[(size_t)k * 9 + r * 3];
            if (!a) { f[r] = 1.0; continue; }
            const double l = len3(a);
            for (int m = 0; m < 3; ++m) f[m] = a[m] / (l > 1e-12 ? l : 1.0);
        }

    nn::parallel_for(out_h_, [&](int64_t y0, int64_t y1) {
        for (size_t i = (size_t)y0 * out_w_; i < (size_t)y1 * out_w_; ++i) {
            double wd = 0, ld = 0, wn = 0, an[3] = {0, 0, 0};
            for (int64_t c = contrib_off_[i]; c < contrib_off_[i + 1]; ++c) {
                const Contrib& t = contrib_[(size_t)c];
                const size_t k = (size_t)t.face;
                float ok = 0;
                sample(face_valid_[k].data(), fw_, fh_, 1, t.px, t.py, &ok);
                // Where the face saw the mid-grey fill instead of the frame,
                // the network still answered; it just answered about nothing.
                const double w = (double)t.w * ok;
                if (w <= 0.0) continue;

                // Linear depth is undefined past 90 degrees; the trainer's
                // own warp drops those rather than divide by a zero cosine.
                const double g = ray_depth ? t.sr : t.sz;
                if (do_depth && !depth[k].empty() && g > 1e-6) {
                    float d = 0;
                    sample(depth[k].data(), fw_, fh_, 1, t.px, t.py, &d);
                    if (d > 0.0f) {
                        ld += w * (std::log((double)d * g) + scale[k]);
                        wd += w;
                    }
                }
                if (do_normal && !normal[k].empty()) {
                    float nf[3];
                    sample(normal[k].data(), fw_, fh_, 3, t.px, t.py, nf);
                    const double* f = &frame[k * 9];
                    for (int m = 0; m < 3; ++m)
                        an[m] += w * (nf[0]*f[m] + nf[1]*f[3+m] + nf[2]*f[6+m]);
                    wn += w;
                }
            }
            // Blended in the log, where a scale error lives and where the
            // trainer's depth loss reads it.
            if (wd > 0.0) (*out_depth)[i] = (float)std::exp(ld / wd);
            if (wn > 0.0) {
                const double l = std::sqrt(an[0]*an[0] + an[1]*an[1] + an[2]*an[2]);
                if (l > 1e-12)
                    for (int m = 0; m < 3; ++m)
                        (*out_normal)[i * 3 + m] = (float)(an[m] / l);
            }
        }
    });
}

}  // namespace app
