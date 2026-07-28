// BAL problem loading and preprocessing into the generalized layout:
// per-image 6-DOF poses + per-group intrinsics (camera groups), observations
// sorted by point, per-camera-model observation lists, and the column layout
// of the reduced camera system.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Camera model registry; must match the entry points in sfm/shaders/ba/ba.slang.
// `n_intr` is how many parameters the kernel *reads* for this model (its
// ICameraModel::kNumIntr). How many of them bundle adjustment may *change* is
// a property of the problem, not the model -- see BAProblem::Group::n_intr and
// sfm/core/Camera.h camNumFreeParams (D50).
struct ModelDesc {
    const char* name;
    uint32_t n_intr;
    const char* cost_entry;
    const char* jac_entry;
};
static const ModelDesc kModels[] = {
    {"snavely", 3, "cost_snavely", "jac_snavely"},                       // BAL, log-focal
    {"snavely_f", 3, "cost_snavely_f", "jac_snavely_f"},                 // BAL, direct focal
    {"pinhole_radial", 5, "cost_pinhole_radial", "jac_pinhole_radial"},  // COLMAP RADIAL
    {"opencv", 8, "cost_opencv", "jac_opencv"},                          // COLMAP OPENCV (D29)
    {"simple_pinhole", 3, "cost_simple_pinhole", "jac_simple_pinhole"},  // COLMAP SIMPLE_PINHOLE
    {"pinhole", 4, "cost_pinhole", "jac_pinhole"},                       // COLMAP PINHOLE
    {"opencv_fisheye", 8, "cost_opencv_fisheye", "jac_opencv_fisheye"},  // COLMAP OPENCV_FISHEYE (D29-C)
    {"full_opencv", 9, "cost_full_opencv", "jac_full_opencv"},           // COLMAP FULL_OPENCV, reduced (D34)
    {"thin_prism_fisheye", 12, "cost_thin_prism_fisheye", "jac_thin_prism_fisheye"},  // COLMAP THIN_PRISM_FISHEYE (D34)
    {"equirect", 2, "cost_equirect", "jac_equirect"},                    // COLMAP EQUIRECTANGULAR (D49)
};
static const int kNumModels = sizeof(kModels) / sizeof(kModels[0]);
static const uint32_t kMaxCamDof = 18;  // 6 pose + up to 12 intrinsics (thin-prism); must match ba.slang

struct BAProblem {
    uint32_t num_images = 0, num_points = 0, num_obs = 0;

    // observations, sorted by (point, image)
    std::vector<uint32_t> obs_image, obs_point;
    std::vector<double> obs_xy;          // 2 per obs
    std::vector<uint32_t> obs_ranges;    // num_points+1

    // camera groups
    struct Group {
        uint32_t intr_offset;  // into flat intrinsics array (kModels[model].n_intr entries)
        uint32_t intr_col;     // column base in the reduced system
        uint32_t n_intr;       // *free* params: columns owned, and the dof beyond the
                               // pose. <= kModels[model].n_intr; the parameters
                               // past it are read by the kernels and held
                               // constant (D50). The kernels still read all
                               // kModels[model].n_intr from intr_offset.
        uint32_t model;
    };
    std::vector<uint32_t> image_group;   // per image
    std::vector<Group> groups;

    // parameters (host copy, double)
    std::vector<double> poses;   // 6 per image
    std::vector<double> intr;    // flat
    std::vector<double> points;  // 3 per point

    // per-model observation lists (concatenated) for specialized dispatches
    struct ModelRange { uint32_t model, offset, count; };
    std::vector<uint32_t> model_obs;
    std::vector<ModelRange> model_ranges;

    // A_cp block pool offsets (element index; block is dof x 3)
    std::vector<uint32_t> acp_off;
    uint64_t acp_total = 0;

    // pair-aggregated Schur tables: for each unordered image pair {a,b} seen
    // together by at least one point, the list of (obs_a, obs_b) index pairs
    // (one per shared point; a==b entries are the per-obs diagonal terms).
    // Grouped contiguously per pair, split into chunks; chunks with bit 31 of
    // the count set share their pair with other chunks and must use atomics.
    // Built on demand by buildPairTables (the solver skips it on the pure-CG
    // path, where the entry lists would only waste host+device memory).
    std::vector<uint32_t> pair_entries;  // 2 per entry
    std::vector<uint32_t> pair_chunks;   // 2 per chunk: offset, count|flag
    uint32_t num_pair_chunks = 0;
    bool use_pair_schur = false;

    // observations grouped by image (CSR), for the CG path's per-camera
    // kernels; built on demand by buildCamTables. cam_chunks splits each
    // camera's list into fixed-size pieces (one warp each; camera blocks are
    // accumulated with atomics) for occupancy and load balance.
    std::vector<uint32_t> cam_obs_ranges;  // num_images+1
    std::vector<uint32_t> cam_obs;         // num_obs
    std::vector<uint32_t> cam_chunks;      // 3 per chunk: image, start, count
    uint32_t num_cam_chunks = 0;

    uint32_t pose_dim = 0;   // 6 * num_images
    uint32_t total_intr = 0;  // intr.size(): parameters stored (free ones first)
    uint32_t free_intr = 0;   // of those, the ones with columns; n_dim - pose_dim
    uint32_t n_dim = 0;      // camera-side system dimension
};

// true if every intrinsics group is referenced by at most one image (required
// by both the pair-Schur and CG paths for exclusive column ownership)
inline bool exclusiveGroups(const BAProblem& P) {
    std::vector<uint32_t> guse(P.groups.size(), 0);
    for (uint32_t g : P.image_group)
        if (++guse[g] > 1) return false;
    return true;
}

// pair-Schur entry count = sum over points of t(t+1)/2 (cheap; used for VRAM
// estimates before deciding whether to build the tables at all)
inline uint64_t pairEntryCount(const BAProblem& P) {
    uint64_t total = 0;
    for (uint32_t p = 0; p < P.num_points; p++) {
        uint64_t t = P.obs_ranges[p + 1] - P.obs_ranges[p];
        total += t * (t + 1) / 2;
    }
    return total;
}

static const uint64_t kMaxPairEntries = 400ull << 20;  // 3.2 GB of entry data

// Group observations by image (CSR) for the CG path's per-camera kernels.
inline void buildCamTables(BAProblem& P) {
    P.cam_obs_ranges.assign(P.num_images + 1, 0);
    for (uint32_t o = 0; o < P.num_obs; o++) P.cam_obs_ranges[P.obs_image[o] + 1]++;
    for (uint32_t i = 0; i < P.num_images; i++) P.cam_obs_ranges[i + 1] += P.cam_obs_ranges[i];
    P.cam_obs.resize(P.num_obs);
    std::vector<uint32_t> fill(P.cam_obs_ranges.begin(), P.cam_obs_ranges.end() - 1);
    for (uint32_t o = 0; o < P.num_obs; o++) P.cam_obs[fill[P.obs_image[o]]++] = o;

    const uint32_t kChunk = 1024;
    P.cam_chunks.clear();
    for (uint32_t img = 0; img < P.num_images; img++)
        for (uint32_t o = P.cam_obs_ranges[img]; o < P.cam_obs_ranges[img + 1]; o += kChunk) {
            P.cam_chunks.push_back(img);
            P.cam_chunks.push_back(o);
            P.cam_chunks.push_back(std::min(kChunk, P.cam_obs_ranges[img + 1] - o));
        }
    P.num_cam_chunks = (uint32_t)(P.cam_chunks.size() / 3);
}

// Build per-model obs lists and A_cp offsets from obs/image/group tables.
inline void finalizeTables(BAProblem& P) {
    P.model_obs.clear();
    P.model_ranges.clear();
    P.model_obs.reserve(P.num_obs);
    for (int m = 0; m < kNumModels; m++) {
        uint32_t off = (uint32_t)P.model_obs.size();
        for (uint32_t o = 0; o < P.num_obs; o++)
            if (P.groups[P.image_group[P.obs_image[o]]].model == (uint32_t)m)
                P.model_obs.push_back(o);
        uint32_t cnt = (uint32_t)P.model_obs.size() - off;
        if (cnt) P.model_ranges.push_back({(uint32_t)m, off, cnt});
    }
    P.acp_off.resize(P.num_obs);
    uint64_t acc = 0;
    for (uint32_t o = 0; o < P.num_obs; o++) {
        P.acp_off[o] = (uint32_t)acc;
        uint32_t dof = 6 + P.groups[P.image_group[P.obs_image[o]]].n_intr;
        if (dof > kMaxCamDof) throw std::runtime_error("camera dof exceeds kMaxCamDof");
        acc += 3 * (size_t)dof;
    }
    P.acp_total = acc;
    if (acc > 0xFFFFFFFFull) throw std::runtime_error("Acp pool exceeds 32-bit indexing");
    // note: the packed-triangle 32-bit limit (n_dim <~ 65k) applies only to
    // the dense solver and is checked when that path is selected
}

// Build the pair-aggregated Schur tables (see BAProblem). Requires exclusive
// per-image ownership of intrinsics columns (no shared groups); falls back to
// the per-observation atomic kernel otherwise, or when the entry list would
// be unreasonably large (very long tracks).
inline void buildPairTables(BAProblem& P) {
    P.pair_entries.clear();
    P.pair_chunks.clear();
    P.num_pair_chunks = 0;
    P.use_pair_schur = false;
    if (P.num_obs == 0) return;

    // exclusivity: each group referenced by at most one image
    if (!exclusiveGroups(P)) return;

    uint64_t total = pairEntryCount(P);
    if (total > kMaxPairEntries) {
        fprintf(stderr, "[bal] pair-schur disabled (%llu entries)\n",
                (unsigned long long)total);
        return;
    }

    // counting sort by pair key; within a point's track images are strictly
    // increasing, so obs i >= j implies image_i >= image_j
    auto key = [&](uint32_t oi, uint32_t oj) {
        uint64_t a = P.obs_image[oi], b = P.obs_image[oj];
        return a * (a + 1) / 2 + b;
    };
    uint64_t nkeys = (uint64_t)P.num_images * (P.num_images + 1) / 2;
    std::vector<uint32_t> cnt(nkeys + 1, 0);
    for (uint32_t p = 0; p < P.num_points; p++)
        for (uint32_t i = P.obs_ranges[p]; i < P.obs_ranges[p + 1]; i++)
            for (uint32_t j = P.obs_ranges[p]; j <= i; j++)
                cnt[key(i, j) + 1]++;
    for (uint64_t k = 0; k < nkeys; k++) cnt[k + 1] += cnt[k];
    std::vector<uint32_t> fill(cnt.begin(), cnt.end() - 1);
    P.pair_entries.resize(2 * total);
    for (uint32_t p = 0; p < P.num_points; p++)
        for (uint32_t i = P.obs_ranges[p]; i < P.obs_ranges[p + 1]; i++)
            for (uint32_t j = P.obs_ranges[p]; j <= i; j++) {
                uint32_t e = fill[key(i, j)]++;
                P.pair_entries[2 * e] = i;
                P.pair_entries[2 * e + 1] = j;
            }

    // chunk pairs; pairs longer than kChunk are split and flagged for atomics
    const uint32_t kChunk = 1024;
    for (uint64_t k = 0; k < nkeys; k++) {
        uint32_t off = cnt[k], n = cnt[k + 1] - cnt[k];
        if (!n) continue;
        uint32_t flag = n > kChunk ? 0x80000000u : 0;
        for (uint32_t o = 0; o < n; o += kChunk) {
            P.pair_chunks.push_back(off + o);
            P.pair_chunks.push_back(std::min(kChunk, n - o) | flag);
        }
    }
    P.num_pair_chunks = (uint32_t)(P.pair_chunks.size() / 2);
    P.use_pair_schur = true;
    fprintf(stderr, "[bal] pair-schur: %llu entries, %u chunks\n",
            (unsigned long long)total, P.num_pair_chunks);
}


inline BAProblem loadBAL(const std::string& path, int model_id, bool shared_intrinsics) {
    if (model_id != 0 && model_id != 1)
        throw std::runtime_error("BAL loader supports snavely / snavely_f models");
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) throw std::runtime_error("cannot open " + path);
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    std::string text(len, 0);
    if (fread(text.data(), 1, len, fp) != (size_t)len) throw std::runtime_error("read failed");
    fclose(fp);

    char* s = text.data();
    auto nextInt = [&]() { return (uint32_t)strtoul(s, &s, 10); };
    auto nextDouble = [&]() { return strtod(s, &s); };

    BAProblem P;
    uint32_t nc = nextInt(), np = nextInt(), no = nextInt();
    P.num_images = nc;
    P.num_points = np;
    P.num_obs = no;
    fprintf(stderr, "[bal] %u cameras, %u points, %u observations\n", nc, np, no);

    std::vector<uint32_t> cam_idx(no), pnt_idx(no);
    std::vector<double> xy(2 * (size_t)no);
    for (uint32_t i = 0; i < no; i++) {
        cam_idx[i] = nextInt();
        pnt_idx[i] = nextInt();
        xy[2 * (size_t)i] = nextDouble();
        xy[2 * (size_t)i + 1] = nextDouble();
    }
    std::vector<double> cam9(9 * (size_t)nc);
    for (auto& v : cam9) v = nextDouble();
    P.points.resize(3 * (size_t)np);
    for (auto& v : P.points) v = nextDouble();

    // sort observations by (point, image)
    std::vector<uint32_t> order(no);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return pnt_idx[a] == pnt_idx[b] ? cam_idx[a] < cam_idx[b] : pnt_idx[a] < pnt_idx[b];
    });
    P.obs_image.resize(no);
    P.obs_point.resize(no);
    P.obs_xy.resize(2 * (size_t)no);
    for (uint32_t i = 0; i < no; i++) {
        uint32_t o = order[i];
        P.obs_image[i] = cam_idx[o];
        P.obs_point[i] = pnt_idx[o];
        P.obs_xy[2 * (size_t)i] = xy[2 * (size_t)o];
        P.obs_xy[2 * (size_t)i + 1] = xy[2 * (size_t)o + 1];
    }
    P.obs_ranges.assign(np + 1, 0);
    for (uint32_t i = 0; i < no; i++) P.obs_ranges[P.obs_point[i] + 1]++;
    for (uint32_t i = 0; i < np; i++) P.obs_ranges[i + 1] += P.obs_ranges[i];

    // poses + intrinsics; BAL camera = [angle-axis(3), t(3), f, k1, k2]
    P.poses.resize(6 * (size_t)nc);
    for (uint32_t c = 0; c < nc; c++)
        for (int j = 0; j < 6; j++) P.poses[6 * (size_t)c + j] = cam9[9 * (size_t)c + j];

    const uint32_t ni = kModels[model_id].n_intr;  // 3
    auto camIntr = [&](uint32_t c, int j) {
        double v = cam9[9 * (size_t)c + 6 + j];
        if (model_id == 0 && j == 0) v = std::log(v);  // log-focal parameterization
        return v;
    };
    if (shared_intrinsics) {
        P.groups.resize(1);
        P.image_group.assign(nc, 0);
        P.intr.assign(ni, 0.0);
        for (uint32_t c = 0; c < nc; c++)
            for (uint32_t j = 0; j < ni; j++) P.intr[j] += camIntr(c, j) / nc;
        P.groups[0] = {0, 0, ni, (uint32_t)model_id};
    } else {
        P.groups.resize(nc);
        P.image_group.resize(nc);
        P.intr.resize((size_t)ni * nc);
        for (uint32_t c = 0; c < nc; c++) {
            P.image_group[c] = c;
            for (uint32_t j = 0; j < ni; j++) P.intr[(size_t)ni * c + j] = camIntr(c, j);
            P.groups[c] = {ni * c, 0 /*fixed below*/, ni, (uint32_t)model_id};
        }
    }
    P.pose_dim = 6 * nc;
    P.total_intr = (uint32_t)P.intr.size();
    P.free_intr = P.total_intr;  // the BAL models refine everything they read
    P.n_dim = P.pose_dim + P.free_intr;
    for (auto& g : P.groups) g.intr_col = P.pose_dim + g.intr_offset;

    finalizeTables(P);
    return P;
}
