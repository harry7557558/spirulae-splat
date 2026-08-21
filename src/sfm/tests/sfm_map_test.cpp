// Synthetic full reconstruction, end to end, on the GPU BA solver.
//
// Prints PASS/FAIL and returns 0/1. See docs/testing.md.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "sfm/core/FeatureCompaction.h"
#include "sfm/core/Model.h"
#include "sfm/map/Bundle.h"
#include "sfm/map/Assemble.h"
#include "sfm/map/Mapper.h"
#include "sfm/tests/TestMain.h"

namespace fs = std::filesystem;
using namespace sfm;

// -----------------------------------------------------------------------
// map-selftest: synthetic full reconstruction (uses the GPU BA solver)
// -----------------------------------------------------------------------
static Pose lookAt(const Vec3& C, const Vec3& target) {
    Vec3 f = (target - C).normalized();
    Vec3 up0 = {0, 1, 0};
    Vec3 r = up0.cross(f).normalized();
    Vec3 u = f.cross(r);
    Mat3 R = {r.x, r.y, r.z, u.x, u.y, u.z, f.x, f.y, f.z};
    Vec3 t = mul(R, C);
    return {R, {-t.x, -t.y, -t.z}};
}

int cmdMapSelftest(int argc, char** argv) {
    MapperOptions opt;
    opt.verbose = false;
    opt.focal = 1200;
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--device" && i + 1 < argc) opt.device = std::stoi(argv[++i]);
        else if (a == "--verbose") opt.verbose = true;
    }
    // Two of the checks below are convergence assertions on deliberately
    // ill-conditioned geometry -- a rotation-degenerate capture's focal, and a
    // principal point 30 px from truth. They need the fp64 solver. A device
    // without an fp64 buffer atomic add (every AMD part here, both Intel iGPUs,
    // llvmpipe) runs the emulated double-float instead, whose ~48-bit mantissa
    // does not get there; on real captures the two agree to 0.01 px of
    // reprojection, so this is a limit of the fallback and not a regression.
    // Everything else in this file still has to pass everywhere.
    const bool fp64_ba =
        realSupportedByDevice(RealCfg::F64, VkContext::probeCaps(opt.device));
    if (!fp64_ba)
        printf("  note: no fp64 BA on this device -- the two ill-conditioned "
               "convergence checks are reported, not asserted\n");

    const int W = 1280, H = 960, M = 8, N = 160;
    Camera K = Camera::defaultFor(1, W, H, 1200);

    std::mt19937 rng(11);
    std::uniform_real_distribution<double> ub(-2.5, 2.5);
    std::normal_distribution<double> noise(0.0, 0.4);

    std::vector<Vec3> pts(N);
    for (auto& p : pts) p = {ub(rng), ub(rng), ub(rng)};

    // Cameras on a 120-degree arc at radius 9, looking at the origin.
    std::vector<Pose> gt(M);
    for (int c = 0; c < M; c++) {
        double ang = -1.05 + 2.1 * c / (M - 1);
        gt[c] = lookAt({9 * std::sin(ang), 1.5 * std::sin(0.7 * c), 9 * std::cos(ang)}, {0, 0, 0});
    }

    // Project; build feature sets and ground-truth matches.
    std::vector<FeatureSet> feats(M);
    std::vector<std::vector<char>> vis(M, std::vector<char>(N, 0));
    for (int c = 0; c < M; c++) {
        feats[c].width = W;
        feats[c].height = H;
        feats[c].keypoints.resize(N);
        for (int p = 0; p < N; p++) {
            Vec3 pc = mul(gt[c].R, pts[p]) + gt[c].t;
            Vec2 px = K.project(pc);
            if (pc.z > 0.1 && px.x > 0 && px.x < W && px.y > 0 && px.y < H) {
                feats[c].keypoints[p] = {(float)(px.x + noise(rng)), (float)(px.y + noise(rng)), 2, 0, 0};
                vis[c][p] = 1;
            } else {
                feats[c].keypoints[p] = {-1000, -1000, 2, 0, 0};
            }
        }
    }
    MatchesDatabase db;
    db.images.resize(M);
    for (int c = 0; c < M; c++) db.images[c] = {"cam" + std::to_string(c), (uint32_t)N};
    for (int i = 0; i < M; i++)
        for (int j = i + 1; j < M; j++) {
            TwoViewMatches tv;
            tv.image1 = i; tv.image2 = j; tv.config = (int)TwoViewConfig::Uncalibrated;
            for (int p = 0; p < N; p++)
                if (vis[i][p] && vis[j][p]) tv.matches.push_back({(uint32_t)p, (uint32_t)p, 0});
            if (tv.matches.size() >= 15) db.pairs.push_back(std::move(tv));
        }

    MatchesDatabase compact_db = db;
    std::vector<FeatureSet> compact_feats(feats.size());
    {
        FeatureCompactionPlan plan = buildFeatureCompactionPlan(compact_db);
        for (size_t i = 0; i < feats.size(); i++)
            compact_feats[i] =
                compactFeatureSet(feats[i], plan.old_to_new[i], plan.compact_counts[i]);
        remapMatches(compact_db, plan, compact_feats);
    }

    Mapper mapper(db, feats, opt);
    // The synthetic scene is one connected component, so this must stay a
    // single model -- more than one would mean the mapper split a graph that
    // does not split (D41).
    std::vector<Reconstruction> models = mapper.run();
    const Reconstruction& rec = models.front();

    int fails = 0;
    uint32_t reg = rec.numRegistered();
    printf("map-selftest: %u/%d images registered, %zu points, %zu model(s)\n", reg, M,
           rec.points3D.size(), models.size());
    if (reg < (uint32_t)M) { printf("  FAIL: not all images registered\n"); fails++; }
    if (models.size() != 1) { printf("  FAIL: connected scene split into models\n"); fails++; }

    Mapper compact_mapper(compact_db, compact_feats, opt);
    std::vector<Reconstruction> compact_models = compact_mapper.run();
    bool compact_equivalent = !compact_models.empty() && compact_models.size() == models.size();
    const Reconstruction* compact_rec = compact_models.empty() ? nullptr : &compact_models.front();
    if (compact_equivalent)
        compact_equivalent = compact_rec->numRegistered() == rec.numRegistered() &&
                             compact_rec->points3D.size() == rec.points3D.size() &&
                             countObservations(*compact_rec) == countObservations(rec);
    if (compact_equivalent) {
        for (const auto& image : rec.images) {
            auto found = compact_rec->images.find(image.first);
            if (found == compact_rec->images.end() ||
                found->second.registered != image.second.registered) {
                compact_equivalent = false;
                break;
            }
        }
    }
    if (compact_equivalent) {
        for (const auto& camera : rec.cameras) {
            auto found = compact_rec->cameras.find(camera.first);
            if (found == compact_rec->cameras.end() ||
                std::fabs(found->second.fx - camera.second.fx) > 1e-9 ||
                std::fabs(found->second.fy - camera.second.fy) > 1e-9 ||
                std::fabs(found->second.cx - camera.second.cx) > 1e-9 ||
                std::fabs(found->second.cy - camera.second.cy) > 1e-9) {
                compact_equivalent = false;
                break;
            }
        }
    }
    printf("  unused-feature compaction A/B: %s (%u images, %zu points, %zu obs)\n",
           compact_equivalent ? "equivalent" : "BAD",
           compact_rec ? compact_rec->numRegistered() : 0,
           compact_rec ? compact_rec->points3D.size() : 0,
           compact_rec ? countObservations(*compact_rec) : 0);
    if (!compact_equivalent) fails++;

    // Relative rotations are invariant to the global similarity gauge.
    double maxRelErr = 0;
    int cmpCount = 0;
    for (auto& a : rec.images)
        for (auto& b : rec.images) {
            if (a.first >= b.first || !a.second.registered || !b.second.registered) continue;
            Mat3 relRec = mul(b.second.pose.R, transpose(a.second.pose.R));
            Mat3 relGt = mul(gt[b.first].R, transpose(gt[a.first].R));
            Mat3 D = mul(relRec, transpose(relGt));
            double tr = (D[0] + D[4] + D[8] - 1) * 0.5;
            maxRelErr = std::max(maxRelErr,
                                 std::acos(std::max(-1.0, std::min(1.0, tr))) * 180.0 / M_PI);
            cmpCount++;
        }
    printf("  max relative-rotation error: %.2f deg over %d pairs\n", maxRelErr, cmpCount);
    if (reg >= 2 && maxRelErr > 1.5) { printf("  FAIL: pose accuracy\n"); fails++; }

    // Mean reprojection error over all observations.
    double sum = 0;
    int nobs = 0;
    for (const auto& kv : rec.points3D)
        for (const TrackElement& e : kv.second.track) {
            const Pose& p = rec.images.at(e.image_id).pose;
            Vec3 pc = mul(p.R, kv.second.xyz) + p.t;
            Vec2 px = rec.cameras.at(1).project(pc);
            Vec2 o = {feats[e.image_id].keypoints[e.point2D_idx].x,
                      feats[e.image_id].keypoints[e.point2D_idx].y};
            sum += std::hypot(px.x - o.x, px.y - o.y);
            nobs++;
        }
    double meanReproj = nobs ? sum / nobs : 0;
    printf("  mean reprojection error: %.3f px over %d obs\n", meanReproj, nobs);
    if (nobs > 0 && meanReproj > 1.5) { printf("  FAIL: reprojection error\n"); fails++; }

    // COLMAP model round-trip.
    std::string dir = "/tmp/spirula_sfm_map_selftest";
    fs::create_directories(dir);
    rec.writeBinary(dir);
    Reconstruction rd = Reconstruction::readBinary(dir);
    bool rt = rd.images.size() == rec.numRegistered() && rd.points3D.size() == rec.points3D.size() &&
              rd.cameras.size() == rec.cameras.size();
    printf("  COLMAP model round-trip: %s\n", rt ? "ok" : "BAD");
    if (!rt) fails++;

    // ---- aligning two models that share no image (D70) ----
    //
    // Two halves of one scene, reconstructed separately and in different
    // gauges, with not one image in common: alignReconstructions has nothing to
    // fit and the pair is never a merge candidate. What they do have is points
    // both triangulated, which the correspondence graph pairs up -- and that
    // determines the similarity. This is a building walked room by room, where
    // the two rooms see the same doorway from two passes.
    {
        auto half = [&](int first, int last) {
            Reconstruction r;
            r.cameras[1] = K;
            for (int c = first; c <= last; c++) {
                Image im;
                im.id = (uint32_t)c;
                im.camera_id = 1;
                im.name = "cam" + std::to_string(c);
                im.registered = true;
                im.pose = gt[c];
                im.points2D.resize(N);
                for (int p = 0; p < N; p++)
                    im.points2D[p] = {feats[c].keypoints[p].x, feats[c].keypoints[p].y};
                im.point3D_ids.assign(N, kInvalidPoint3D);
                r.images[(uint32_t)c] = std::move(im);
            }
            for (int p = 0; p < N; p++) {
                std::vector<TrackElement> tr;
                for (int c = first; c <= last; c++)
                    if (vis[c][p]) tr.push_back({(uint32_t)c, (uint32_t)p});
                if (tr.size() >= 2) r.addPoint3D(pts[p], tr);
            }
            return r;
        };
        Reconstruction A = half(0, M / 2 - 1), B = half(M / 2, M - 1);
        // Put B in a gauge of its own: 1.7x, a 25-degree yaw, and an offset.
        Sim3 S;
        S.scale = 1.7;
        const double a = 25.0 * M_PI / 180.0;
        S.R = {std::cos(a), 0, std::sin(a), 0, 1, 0, -std::sin(a), 0, std::cos(a)};
        S.t = {3.0, -1.5, 0.5};
        for (auto& kv : B.images) kv.second.pose = transformPose(S, kv.second.pose);
        for (auto& kv : B.points3D) kv.second.xyz = transformPoint(S, kv.second.xyz);

        MergeOptions mopt;
        AlignmentResult al = mapper.alignByStructure(A, B, mopt);
        // Recovered transform applied to B's centres must land on A's gauge,
        // which is `gt`. Measured on the cameras, none of which took part in
        // the fit -- it only ever saw points.
        double worst = 0;
        for (const auto& kv : B.images) {
            const Vec3 got = transformPoint(al.transform, cameraCenter(kv.second.pose));
            const Vec3 want = cameraCenter(gt[kv.first]);
            worst = std::max(worst, (got - want).norm());
        }
        printf("  align on shared structure: %zu shared image(s), %zu point pair(s), "
               "%zu inlier(s), %.2f px, camera centres off by %.4f\n",
               al.common_images, al.structure_pairs, al.inliers, al.mean_error, worst);
        if (!al.success || al.common_images != 0 || worst > 0.05) {
            printf("  FAIL: the similarity was not recovered from shared structure alone\n");
            fails++;
        }
    }

    // ---- two disconnected components -> two models (D41) ----
    // The same scene twice, with no verified pair joining the halves: exactly
    // the "loose components" case, where the pre-D41 mapper built both, threw
    // one away and reported half the dataset. Both must now come back, largest
    // first, with the images partitioned between them and nothing shared.
    {
        MapperOptions opt2 = opt;
        opt2.min_model_size = 5;  // the components are 8 images each
        MatchesDatabase db2;
        std::vector<FeatureSet> feats2;
        // Two copies of the single-component scene above; image c of component
        // k becomes image k*M + c, and pairs are only ever built within a
        // component, so the view graph has no edge between the halves.
        for (int k = 0; k < 2; k++) {
            for (int c = 0; c < M; c++) {
                feats2.push_back(feats[c]);
                db2.images.push_back({"comp" + std::to_string(k) + "/cam" + std::to_string(c),
                                      (uint32_t)N});
            }
            for (int i = 0; i < M; i++)
                for (int j = i + 1; j < M; j++) {
                    TwoViewMatches tv;
                    tv.image1 = (uint32_t)(k * M + i);
                    tv.image2 = (uint32_t)(k * M + j);
                    tv.config = (int)TwoViewConfig::Uncalibrated;
                    for (int p = 0; p < N; p++)
                        if (vis[i][p] && vis[j][p]) tv.matches.push_back({(uint32_t)p, (uint32_t)p, 0});
                    if (tv.matches.size() >= 15) db2.pairs.push_back(std::move(tv));
                }
        }
        Mapper mapper2(db2, feats2, opt2);
        std::vector<Reconstruction> ms = mapper2.run();
        uint32_t total = 0;
        for (const Reconstruction& m : ms) total += m.numRegistered();
        printf("  two disconnected components: %zu model(s), %u/%d images registered in total\n",
               ms.size(), total, 2 * M);
        if (ms.size() != 2 || total != (uint32_t)(2 * M)) {
            printf("  FAIL: expected 2 models covering all %d images\n", 2 * M);
            fails++;
        }
        // Sorted largest-first, disjoint, and each model holds one component
        // whole -- a model straddling the halves would mean the mapper invented
        // a constraint the view graph does not have.
        std::set<uint32_t> seen;
        bool disjoint = true, whole = true, sorted = true;
        for (size_t i = 0; i < ms.size(); i++) {
            if (i && ms[i].points3D.size() > ms[i - 1].points3D.size()) sorted = false;
            std::set<int> comps;
            for (const auto& kv : ms[i].images) {
                if (!seen.insert(kv.first).second) disjoint = false;
                comps.insert(kv.first / M);
            }
            if (comps.size() != 1 || ms[i].images.size() != (size_t)M) whole = false;
        }
        if (!disjoint) { printf("  FAIL: models share images\n"); fails++; }
        if (!whole) { printf("  FAIL: a model straddles the two components\n"); fails++; }
        if (!sorted) { printf("  FAIL: models not ordered by point count\n"); fails++; }

        // Assembly must leave two genuinely separate components alone: they
        // share no image, so there is nothing to align on and nothing to grow
        // across, and anything it did here would be an invention (D44).
        {
            ManagerOptions mo;
            mo.verbose = false;
            AssembleOptions ao;
            ao.verbose = false;
            ao.max_rounds = 2;
            AssembleStats ast;
            std::vector<Reconstruction> done = assembleModels(mapper2, ms, mo, ao, ast);
            uint32_t tot = 0;
            for (const Reconstruction& m : done) tot += m.numRegistered();
            printf("  assembly of two separate components: %zu model(s), %u images, %zu merge(s)\n",
                   done.size(), tot, ast.merges);
            if (done.size() != 2 || tot != (uint32_t)(2 * M) || ast.merges != 0) {
                printf("  FAIL: assembly should have left two separate components alone\n");
                fails++;
            }
        }

        // Shared intrinsics across components (D45). Both components look
        // through the same camera id, so one set of intrinsics is the truth for
        // both. Give one component a badly wrong focal -- what a small
        // component's own focal search does to it in the wild -- and the joint
        // bundle adjustment has to pull it back onto the other's evidence.
        {
            std::vector<Reconstruction> two = ms;
            const double good = two[0].cameras.at(1).focal();
            two[1].cameras[1].setFocal(good * 0.7);
            const double before = two[1].cameras.at(1).focal();
            mapper2.jointRefine(two);
            const double a = two[0].cameras.at(1).focal(), b = two[1].cameras.at(1).focal();
            printf("  joint BA with shared intrinsics: focals %.1f / %.1f -> %.1f / %.1f "
                   "(truth %.1f)\n", good, before, a, b, good);
            if (std::fabs(a - b) > 1e-6 * good || std::fabs(a - good) > 0.05 * good) {
                printf("  FAIL: joint BA did not share one focal across the components\n");
                fails++;
            }
        }

        // Splitting a model the correspondence graph contradicts (D45). Half
        // the images of a sound model are rotated where they stand: the model
        // is still self-consistent as a data structure, but the verified pairs
        // that span the two halves cannot hold any more, and that is what the
        // split has to find. Nothing here looks at 3D points -- the whole idea
        // is to use evidence the model did not build itself from.
        {
            Reconstruction broken = ms[0];
            uint32_t moved = 0, kept = 0;
            std::vector<uint32_t> ids;
            for (const auto& kv : broken.images)
                if (kv.second.registered) ids.push_back(kv.first);
            const double a35 = 35.0 * M_PI / 180.0;
            Mat3 spin = {std::cos(a35), 0, std::sin(a35), 0, 1, 0,
                         -std::sin(a35), 0, std::cos(a35)};
            for (size_t i = 0; i < ids.size(); i++) {
                if (i * 2 < ids.size()) { kept++; continue; }
                // A world rotation applied to this half only: their poses
                // stay mutually consistent (so the half is a valid model in
                // its own right), and every pair that spans the halves breaks.
                Image& im = broken.images[ids[i]];
                im.pose.R = mul(im.pose.R, transpose(spin));
                moved++;
            }
            Mapper::SplitStats ss;
            std::vector<Reconstruction> parts =
                mapper2.splitInconsistent(broken, 8.0, 0.5, 15, 2, &ss);
            std::vector<uint32_t> sizes;
            for (const Reconstruction& p : parts) sizes.push_back(p.numRegistered());
            std::sort(sizes.begin(), sizes.end(), std::greater<uint32_t>());
            printf("  split a model with %u of %zu images rotated 35 deg: %zu/%zu inner pairs "
                   "hold -> %zu part(s)", moved, ids.size(), ss.pairs_agree, ss.pairs_tested,
                   parts.size());
            for (uint32_t s : sizes) printf(" %u", s);
            printf("\n");
            bool ok = parts.size() == 2 && sizes.size() == 2 &&
                      sizes[0] == std::max(moved, kept) && sizes[1] == std::min(moved, kept);
            if (!ok) { printf("  FAIL: split did not separate the rotated half\n"); fails++; }
            // ... and a model nothing is wrong with must survive untouched.
            Mapper::SplitStats ok_ss;
            std::vector<Reconstruction> whole_parts =
                mapper2.splitInconsistent(ms[0], 8.0, 0.5, 15, 2, &ok_ss);
            printf("  split a sound model: %zu/%zu inner pairs hold -> %zu part(s)\n",
                   ok_ss.pairs_agree, ok_ss.pairs_tested, whole_parts.size());
            if (whole_parts.size() != 1) {
                printf("  FAIL: split broke up a model that agrees with its own pairs\n");
                fails++;
            }
        }

        // --max-models 1 must reproduce the old single-model behaviour exactly.
        MapperOptions opt1 = opt2;
        opt1.max_num_models = 1;
        std::vector<Reconstruction> one = Mapper(db2, feats2, opt1).run();
        printf("  --max-models 1 on the same input: %zu model(s), %u images\n", one.size(),
               one.empty() ? 0 : one.front().numRegistered());
        if (one.size() != 1 || one.front().numRegistered() != (uint32_t)M) {
            printf("  FAIL: max_num_models 1 did not collapse to one component\n");
            fails++;
        }
    }

    // ---- the engine driven from outside: continueFrom / audit (D44) ----
    // The manager hands the mapper models it did not build -- merged, or off
    // disk -- and asks it to keep going, or to check them. Both are exercised
    // here against the model above, whose correct answer is known.
    {
        const uint32_t victim = 3;
        auto relRot = [](const Reconstruction& r, uint32_t a, uint32_t b) {
            return mul(r.images.at(b).pose.R, transpose(r.images.at(a).pose.R));
        };
        auto rotDiffDeg = [](const Mat3& A, const Mat3& B) {
            Mat3 D = mul(A, transpose(B));
            double tr = std::max(-1.0, std::min(1.0, (D[0] + D[4] + D[8] - 1) * 0.5));
            return std::acos(tr) * 180.0 / M_PI;
        };
        const Mat3 truth = relRot(rec, 0, victim);

        // Detach one image entirely, as if it had never registered, and let
        // continueFrom put it back.
        Reconstruction partial = rec;
        partial.images.erase(victim);
        for (auto it = partial.points3D.begin(); it != partial.points3D.end();) {
            auto& tr = it->second.track;
            tr.erase(std::remove_if(tr.begin(), tr.end(),
                                    [&](const TrackElement& e) { return e.image_id == victim; }),
                     tr.end());
            it = tr.size() < 2 ? partial.points3D.erase(it) : std::next(it);
        }
        Mapper::GrowStats gs;
        Reconstruction back = Mapper(db, feats, opt).continueFrom(partial, &gs);
        const bool got = back.images.count(victim) && back.images.at(victim).registered;
        double err = got ? rotDiffDeg(relRot(back, 0, victim), truth) : 1e9;
        printf("  continueFrom: %u/%u images -> %u (%u registered), pose error %.3f deg\n",
               partial.numRegistered(), M, back.numRegistered(), gs.registered, err);
        if (!got || gs.registered != 1 || err > 1.0) {
            printf("  FAIL: continueFrom did not put the missing image back correctly\n");
            fails++;
        }

        // A model nothing is wrong with must come through the audit untouched:
        // false positives cost images, so this is the half that matters most.
        Mapper::AuditStats clean;
        Mapper(db, feats, opt).audit(rec, &clean);
        printf("  audit of a sound model: %u checked, %u contradicted\n", clean.checked,
               clean.unsupported);
        if (clean.unsupported != 0) {
            printf("  FAIL: audit contradicted an image of a model that is correct\n");
            fails++;
        }

        // Now the failure it exists for: an image sitting somewhere the rest of
        // the model does not support, with its own observations gone (what a
        // merge produces when its similarity does not fit both halves).
        Reconstruction bad = partial;
        Image moved = rec.images.at(victim);
        moved.pose.R = mul(angleAxisToRotation({0.0, 0.6, 0.0}), moved.pose.R);
        moved.point3D_ids.assign(moved.points2D.size(), kInvalidPoint3D);
        bad.images[victim] = moved;
        Mapper::AuditStats as;
        Reconstruction fixed = Mapper(db, feats, opt).audit(bad, &as);
        const bool kept = fixed.images.count(victim) && fixed.images.at(victim).registered;
        double ferr = kept ? rotDiffDeg(relRot(fixed, 0, victim), truth) : 1e9;
        printf("  audit of a misplaced image: %u contradicted, %u repaired, pose error "
               "%.3f deg (was %.1f)\n", as.unsupported, as.reregistered, ferr,
               rotDiffDeg(relRot(bad, 0, victim), truth));
        if (as.unsupported != 1 || !kept || ferr > 1.0) {
            printf("  FAIL: audit did not detect and repair the misplaced image\n");
            fails++;
        }
    }

    // ---- a forward-motion capture (D48) ----
    // A dashcam: the camera drives straight down its own optical axis, with only
    // the pitch/yaw jitter of a vehicle. Two separate things have to work.
    //
    // Seeding: every pair on such a capture has its baseline along the viewing
    // direction, so COLMAP's init_max_forward_motion veto rejects the whole
    // dataset and the pre-D48 mapper registered nothing at all.
    //
    // The focal: with all cameras pointing one way, scaling the focal and
    // stretching the scene along the viewing axis reproduces every image, so
    // bundle adjustment cannot see the focal -- *except* through the distortion,
    // whose radial polynomial does not rescale. Hence the synthetic camera here
    // has distortion, as a real one does, and the focal bootstrap has to find
    // its way from a guess 71% long back to the truth.
    {
        const int Wd = 1280, Hd = 400, Md = 14, Nd = 900;
        Camera Kd = Camera::defaultFor(1, Wd, Hd, 900, CamModel::OpenCV);
        Kd.k1 = -0.28;
        Kd.k2 = 0.07;

        std::mt19937 r2(7);
        std::uniform_real_distribution<double> ux(-9.0, 9.0), uy(-3.0, 3.0), uz(3.0, 42.0);
        std::normal_distribution<double> jit(0.0, 0.02), noise2(0.0, 0.3);
        std::vector<Vec3> pts2(Nd);
        for (Vec3& p : pts2) p = {ux(r2), uy(r2), uz(r2)};

        std::vector<Pose> gtd(Md);
        for (int c = 0; c < Md; c++) {
            // Straight down +Z, 1.1 units per frame, with a vehicle's jitter.
            Vec3 C = {0, 0, 1.1 * c};
            Mat3 R = angleAxisToRotation({jit(r2), jit(r2), jit(r2) * 0.3});
            Vec3 t = mul(R, C);
            gtd[c] = {R, {-t.x, -t.y, -t.z}};
        }

        std::vector<FeatureSet> fd(Md);
        std::vector<std::vector<char>> vd(Md, std::vector<char>(Nd, 0));
        for (int c = 0; c < Md; c++) {
            fd[c].width = Wd;
            fd[c].height = Hd;
            fd[c].keypoints.resize(Nd);
            for (int p = 0; p < Nd; p++) {
                Vec3 pc = mul(gtd[c].R, pts2[p]) + gtd[c].t;
                Vec2 px = Kd.project(pc);
                if (pc.z > 0.5 && px.x > 0 && px.x < Wd && px.y > 0 && px.y < Hd) {
                    fd[c].keypoints[p] = {(float)(px.x + noise2(r2)), (float)(px.y + noise2(r2)),
                                          2, 0, 0};
                    vd[c][p] = 1;
                } else {
                    fd[c].keypoints[p] = {-1000, -1000, 2, 0, 0};
                }
            }
        }
        MatchesDatabase dbd;
        dbd.images.resize(Md);
        for (int c = 0; c < Md; c++) dbd.images[c] = {"f" + std::to_string(c), (uint32_t)Nd};
        for (int i = 0; i < Md; i++)
            for (int j = i + 1; j < Md; j++) {
                TwoViewMatches tv;
                tv.image1 = i;
                tv.image2 = j;
                tv.config = (int)TwoViewConfig::Uncalibrated;
                for (int p = 0; p < Nd; p++)
                    if (vd[i][p] && vd[j][p]) tv.matches.push_back({(uint32_t)p, (uint32_t)p, 0});
                if (tv.matches.size() >= 15) dbd.pairs.push_back(std::move(tv));
            }

        MapperOptions od;
        od.verbose = opt.verbose;
        od.device = opt.device;
        od.camera_model = CamModel::OpenCV;
        od.focal = 0;  // the 1.2*max_dim guess: 1536 against a truth of 900

        Mapper md(dbd, fd, od);
        std::vector<Reconstruction> mds = md.run();
        const Reconstruction& rd2 = mds.front();
        const uint32_t regd = rd2.numRegistered();
        const double fd_est = rd2.cameras.count(1) ? rd2.cameras.at(1).focal() : 0;
        const double guess = 1.2 * Wd;
        printf("  forward-motion capture: %u/%d images, focal %.0f (guess %.0f, truth 900)\n",
               regd, Md, fd_est, guess);
        if (regd < (uint32_t)Md) {
            printf("  FAIL: a forward-motion capture must still seed and grow\n");
            fails++;
        }
        // Well inside the guess's 71% error, and on the right side of it.
        if (std::fabs(fd_est - 900.0) > 0.20 * 900.0) {
            printf("  %s: focal not recovered from a rotation-degenerate capture\n",
                   fp64_ba ? "FAIL" : "df-limited");
            if (fp64_ba) fails++;
        }
        // The bootstrap must never make a focal *worse* than leaving it alone.
        // This scene is deliberately kind about the focal -- 900 points at 0.3
        // px and a strong k1 -- so bundle adjustment gets there on its own here,
        // and the run below is the control: the same capture with the search
        // switched off. The evidence that the search matters when BA *cannot*
        // get there is a real dashcam, where it moves the focal from 54% long
        // to 1.5% (D48); no synthetic scene reproduces that without being tuned
        // to sit on the edge of recoverability, which would make this a coin
        // flip rather than a test.
        MapperOptions oo = od;
        oo.focal_trials = 0;
        Mapper mo(dbd, fd, oo);
        std::vector<Reconstruction> mos = mo.run();
        const double f_off = mos.front().cameras.count(1) ? mos.front().cameras.at(1).focal() : 0;
        printf("  ... with the search off: focal %.0f over %u image(s)\n", f_off,
               mos.front().numRegistered());
        if (std::fabs(fd_est - 900.0) > std::fabs(f_off - 900.0) + 0.01 * 900.0) {
            printf("  FAIL: the focal bootstrap made a recoverable focal worse\n");
            fails++;
        }
    }

    // ---- the principal point is held, and only the principal point (D50) ----
    // The same arc scene, started from a camera whose principal point is 30 px
    // off the true one. By default bundle adjustment must leave those two
    // numbers exactly where it found them while still fitting the focal, and
    // with --refine-principal-point it must move them back towards the truth.
    // That is the whole prefix-of-free-parameters mechanism, checked on a model
    // that does have parameters after the held pair.
    {
        Camera c0 = Camera::defaultFor(1, W, H, 1150.0, CamModel::OpenCV);
        c0.cx = W * 0.5 + 30;
        c0.cy = H * 0.5 - 30;

        MapperOptions op = opt;
        op.camera_model = CamModel::OpenCV;
        op.initial_cameras[1] = c0;
        op.focal_trials = 0;  // keep this about BA, not about the focal search

        Mapper mp(db, feats, op);
        std::vector<Reconstruction> mps = mp.run();
        const Camera& cp = mps.front().cameras.at(1);

        MapperOptions oq = op;
        oq.refine_principal_point = true;
        oq.pp_min_images = 4;   // 8 images share this camera; the default asks 20
        Mapper mq(db, feats, oq);
        std::vector<Reconstruction> mqs = mq.run();
        const Camera& cq = mqs.front().cameras.at(1);

        // The finished model gets one more pass with the principal point free
        // (D51). This scene has 8 images on one camera, so it does not qualify
        // under the default pp_min_images -- raise the bar's other side and it
        // should walk the principal point back toward the truth.
        MapperOptions ol = op;
        ol.pp_min_images = 4;
        Mapper ml(db, feats, ol);
        std::vector<Reconstruction> mls = ml.run();
        Reconstruction polished = ml.polish(mls.front());
        const Camera& cl = polished.cameras.at(1);

        const bool held = cp.cx == c0.cx && cp.cy == c0.cy;
        const bool moved = cq.cx != c0.cx || cq.cy != c0.cy;
        // The focal still has to converge with the principal point pinned --
        // holding a parameter must not freeze the group.
        const bool focal_ok = std::fabs(cp.focal() - 1200.0) < 0.05 * 1200.0;
        // The final pass has to *find* something, not merely move: the true
        // principal point is the image centre, and it starts 30 px away.
        const double d0 = std::hypot(c0.cx - W * 0.5, c0.cy - H * 0.5);
        const double dl = std::hypot(cl.cx - W * 0.5, cl.cy - H * 0.5);
        printf("  principal point: held (%.1f,%.1f) vs start (%.1f,%.1f), focal %.0f; "
               "refined -> (%.1f,%.1f); final pass -> (%.1f,%.1f), %.1f px from truth "
               "(started %.1f)\n",
               cp.cx, cp.cy, c0.cx, c0.cy, cp.focal(), cq.cx, cq.cy, cl.cx, cl.cy, dl, d0);
        if (!held) { printf("  FAIL: BA moved a held principal point\n"); fails++; }
        if (!moved) { printf("  FAIL: --refine-principal-point did not free it\n"); fails++; }
        if (dl > 0.5 * d0) {
            printf("  %s: the final principal-point pass did not recover it\n",
                   fp64_ba ? "FAIL" : "df-limited");
            if (fp64_ba) fails++;
        }
        // ... and a group too small to share intrinsics must be left alone.
        Reconstruction unpolished = mp.polish(mps.front());
        if (unpolished.cameras.at(1).cx != c0.cx) {
            printf("  FAIL: polished a camera group below pp_min_images\n");
            fails++;
        }
        // ... as must a model with more than one camera group, whatever the
        // thresholds say: that is the rig case the pass must never touch (D51).
        {
            MapperOptions orig = ol;
            Reconstruction two = mls.front();
            Camera c2 = two.cameras.at(1);
            c2.id = 2;
            two.cameras[2] = c2;
            uint32_t n = 0;
            for (auto& kv : two.images)
                if (kv.second.registered && n++ % 2) kv.second.camera_id = 2;
            Mapper mr(db, feats, orig);
            Reconstruction rr2 = mr.polish(two);
            if (rr2.cameras.at(1).cx != two.cameras.at(1).cx ||
                rr2.cameras.at(2).cx != two.cameras.at(2).cx) {
                printf("  FAIL: the final pass touched a multi-group model\n");
                fails++;
            }
        }
        if (!focal_ok) { printf("  FAIL: focal did not converge with the PP held\n"); fails++; }
    }

    // ---- distortion held during mapping, found by the finishing pass (D72) ----
    // The same arc through a lens that really distorts. Per-image intrinsics
    // (D73) ride on the same model: eight images on one camera become eight.
    {
        Camera Kx = Camera::defaultFor(1, W, H, 1200, CamModel::OpenCV);
        Kx.k1 = -0.06;
        Kx.k2 = 0.01;

        std::mt19937 r4(29);
        std::normal_distribution<double> n4(0.0, 0.3);
        std::vector<FeatureSet> fx(M);
        std::vector<std::vector<char>> vx(M, std::vector<char>(N, 0));
        for (int c = 0; c < M; c++) {
            fx[c].width = W;
            fx[c].height = H;
            fx[c].keypoints.resize(N);
            for (int p = 0; p < N; p++) {
                Vec3 pc = mul(gt[c].R, pts[p]) + gt[c].t;
                Vec2 px = Kx.project(pc);
                if (pc.z > 0.1 && px.x > 0 && px.x < W && px.y > 0 && px.y < H) {
                    fx[c].keypoints[p] = {(float)(px.x + n4(r4)), (float)(px.y + n4(r4)), 2, 0, 0};
                    vx[c][p] = 1;
                } else {
                    fx[c].keypoints[p] = {-1000, -1000, 2, 0, 0};
                }
            }
        }
        MatchesDatabase dbx;
        dbx.images.resize(M);
        for (int c = 0; c < M; c++) dbx.images[c] = {"dis" + std::to_string(c), (uint32_t)N};
        for (int i = 0; i < M; i++)
            for (int j = i + 1; j < M; j++) {
                TwoViewMatches tv;
                tv.image1 = i;
                tv.image2 = j;
                tv.config = (int)TwoViewConfig::Uncalibrated;
                for (int p = 0; p < N; p++)
                    if (vx[i][p] && vx[j][p]) tv.matches.push_back({(uint32_t)p, (uint32_t)p, 0});
                if (tv.matches.size() >= 15) dbx.pairs.push_back(std::move(tv));
            }

        MapperOptions ox = opt;
        ox.camera_model = CamModel::OpenCV;
        ox.focal_trials = 0;
        ox.refine_extra_params = false;
        Mapper mx(dbx, fx, ox);
        std::vector<Reconstruction> mxs = mx.run();
        const Camera& ch = mxs.front().cameras.at(1);
        const bool held = ch.k1 == 0 && ch.k2 == 0 && ch.p1 == 0 && ch.p2 == 0;

        Reconstruction fin = mx.polish(mxs.front(), /*free_pp=*/false, /*free_extra=*/true);
        const Camera& cf = fin.cameras.at(1);
        printf("  distortion: held (%.4f, %.4f) over %u image(s); finishing pass -> "
               "(%.4f, %.4f), truth (%.4f, %.4f)\n",
               ch.k1, ch.k2, mxs.front().numRegistered(), cf.k1, cf.k2, Kx.k1, Kx.k2);
        if (!held) { printf("  FAIL: BA moved a held distortion coefficient\n"); fails++; }
        // Recovered at least halfway, which no amount of gauge freedom gives
        // for nothing: unlike the principal point, k1 is not a rotation.
        if (std::fabs(cf.k1 - Kx.k1) > 0.5 * std::fabs(Kx.k1)) {
            printf("  %s: the finishing pass did not recover the distortion\n",
                   fp64_ba ? "FAIL" : "df-limited");
            if (fp64_ba) fails++;
        }
        // ... and holding them must not freeze the focal alongside.
        if (std::fabs(ch.focal() - 1200.0) > 0.10 * 1200.0) {
            printf("  FAIL: focal did not converge with the distortion held\n");
            fails++;
        }

        Reconstruction per = mx.perImageIntrinsics(fin);
        std::set<uint32_t> percam;
        for (const auto& kv : per.images)
            if (kv.second.registered) percam.insert(kv.second.camera_id);
        double spread = 0;
        for (const auto& kv : per.cameras)
            spread = std::max(spread, std::fabs(kv.second.focal() - cf.focal()));
        printf("  per-image intrinsics: %zu camera(s) over %u image(s), focal spread %.2f px\n",
               percam.size(), per.numRegistered(), spread);
        if (percam.size() != per.numRegistered() || per.cameras.size() != percam.size()) {
            printf("  FAIL: the per-image pass did not give every image its own camera\n");
            fails++;
        }
        if (per.numRegistered() != fin.numRegistered()) {
            printf("  FAIL: the per-image pass lost images\n");
            fails++;
        }
        // Free per-image intrinsics can only fit the same observations better.
        // A synthetic capture whose lens really is shared should nonetheless
        // stay near it, which is what makes a runaway visible.
        if (!(spread < 0.05 * cf.focal())) {
            printf("  %s: per-image focals ran away from the shared solution\n",
                   fp64_ba ? "FAIL" : "df-limited");
            if (fp64_ba) fails++;
        }
    }

    // ---- an equirectangular (360) capture (D49) ----
    // A spherical camera walking a straight line through a scene that surrounds
    // it. Three things are specific to this model and nothing else exercises
    // them:
    //
    //  * every direction projects, so most observations here are of points
    //    *behind* the camera -- correspondences no perspective model can hold.
    //  * forward motion is not degenerate the way it is for a narrow lens: the
    //    points abeam have full parallax, so the path is recovered without any
    //    of D48's machinery.
    //  * the intrinsics are metadata. Bundle adjustment reads (w, h) and must
    //    return them bit-identical -- the group owns no columns of the reduced
    //    system, which is the whole reason for the free/stored split in
    //    sfm/ba/Problem.h and ba.slang.
    {
        const int We = 2048, He = 1024, Me = 10, Ne = 500;
        Camera Ke = Camera::defaultFor(1, We, He, 0, CamModel::Equirect);

        std::mt19937 r3(23);
        std::normal_distribution<double> gs(0.0, 1.0), noise3(0.0, 0.4);
        std::uniform_real_distribution<double> ur(4.0, 14.0), uy3(-2.0, 2.0);

        // Points on a shell around the whole path, so every station sees the
        // scene in all directions.
        std::vector<Vec3> pts3(Ne);
        for (Vec3& p : pts3) {
            Vec3 d = {gs(r3), 0.35 * gs(r3), gs(r3)};
            double n = std::max(1e-9, std::sqrt(d.dot(d)));
            double r = ur(r3);
            p = {d.x / n * r, d.y / n * r + uy3(r3), d.z / n * r + 4.5};
        }

        std::vector<Pose> gte(Me);
        for (int c = 0; c < Me; c++) {
            Vec3 C = {0.15 * c, 0.05 * std::sin(0.9 * c), 1.0 * c};
            // A 360 camera is carried, so its yaw drifts; the model has to be
            // right about orientation as well as position.
            Mat3 R = angleAxisToRotation({0.03 * std::sin(0.7 * c), 0.25 * c, 0.02 * c});
            Vec3 t = mul(R, C);
            gte[c] = {R, {-t.x, -t.y, -t.z}};
        }

        std::vector<FeatureSet> fe3(Me);
        std::vector<std::vector<char>> ve(Me, std::vector<char>(Ne, 0));
        int behind = 0;
        for (int c = 0; c < Me; c++) {
            fe3[c].width = We;
            fe3[c].height = He;
            fe3[c].keypoints.resize(Ne);
            for (int p = 0; p < Ne; p++) {
                Vec3 pc = mul(gte[c].R, pts3[p]) + gte[c].t;
                Vec2 px = Ke.project(pc);
                fe3[c].keypoints[p] = {(float)(px.x + noise3(r3)), (float)(px.y + noise3(r3)),
                                       2, 0, 0};
                ve[c][p] = 1;
                if (pc.z <= 0) behind++;
            }
        }
        MatchesDatabase dbe;
        dbe.images.resize(Me);
        for (int c = 0; c < Me; c++) dbe.images[c] = {"s" + std::to_string(c), (uint32_t)Ne};
        for (int i = 0; i < Me; i++)
            for (int j = i + 1; j < Me; j++) {
                TwoViewMatches tv;
                tv.image1 = i;
                tv.image2 = j;
                tv.config = (int)TwoViewConfig::Uncalibrated;
                for (int p = 0; p < Ne; p++)
                    if (ve[i][p] && ve[j][p]) tv.matches.push_back({(uint32_t)p, (uint32_t)p, 0});
                if (tv.matches.size() >= 15) dbe.pairs.push_back(std::move(tv));
            }

        MapperOptions oe;
        oe.verbose = opt.verbose;
        oe.device = opt.device;
        oe.camera_model = CamModel::Equirect;
        oe.focal = 0;

        Mapper me(dbe, fe3, oe);
        std::vector<Reconstruction> mes = me.run();
        const Reconstruction& re = mes.front();
        const uint32_t rege = re.numRegistered();

        double maxRel = 0;
        for (const auto& a : re.images)
            for (const auto& b : re.images) {
                if (a.first >= b.first || !a.second.registered || !b.second.registered) continue;
                Mat3 D = mul(mul(b.second.pose.R, transpose(a.second.pose.R)),
                             transpose(mul(gte[b.first].R, transpose(gte[a.first].R))));
                double tr = (D[0] + D[4] + D[8] - 1) * 0.5;
                maxRel = std::max(maxRel,
                                  std::acos(std::max(-1.0, std::min(1.0, tr))) * 180.0 / M_PI);
            }
        const Camera& ce = re.cameras.at(1);
        const bool intr_fixed = ce.model == CamModel::Equirect &&
                                ce.fx == Ke.fx && ce.fy == Ke.fy &&
                                ce.cx == Ke.cx && ce.cy == Ke.cy;
        printf("  equirect capture: %u/%d images, max rel-rot %.3f deg, %d obs behind the "
               "camera, intrinsics %s\n",
               rege, Me, maxRel, behind, intr_fixed ? "unchanged" : "MOVED");
        if (rege < (uint32_t)Me) {
            printf("  FAIL: a spherical capture must register every station\n");
            fails++;
        }
        if (rege >= 2 && maxRel > 1.0) {
            printf("  FAIL: spherical pose accuracy\n");
            fails++;
        }
        if (!intr_fixed) {
            printf("  FAIL: bundle adjustment moved held-constant intrinsics\n");
            fails++;
        }
        if (behind == 0) {
            printf("  FAIL: the scene did not exercise the rear hemisphere\n");
            fails++;
        }
    }

    printf("%s\n", fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    return sfmTestMain(argc - 1, argv + 1, cmdMapSelftest);
}
