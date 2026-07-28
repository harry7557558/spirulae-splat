// Two-view geometric verification (src/sfm/README.md).
//
// Given putative correspondences, robustly estimate a fundamental matrix and a
// homography, pick the model that best explains the data (COLMAP's H-inlier
// ratio test), and return the surviving inliers and a config label. Optionally
// recover the calibrated relative pose from F + intrinsics.
//
// This is the filter that turns the phase-2 matcher's putative matches into
// trustworthy two-view geometries for the phase-4 mapper.
#pragma once

#include <vector>

#include "sfm/geometry/Essential.h"
#include "sfm/geometry/Fundamental.h"
#include "sfm/geometry/Homography.h"
#include "sfm/geometry/LinAlg.h"
#include "sfm/optim/Ransac.h"

namespace sfm {

enum class TwoViewConfig {
    Undefined = 0,
    Degenerate,          // too few inliers
    Uncalibrated,        // F best explains it
    PlanarOrPanoramic,   // H best explains it (planar scene or pure rotation)
};

inline const char* twoViewConfigName(TwoViewConfig c) {
    switch (c) {
        case TwoViewConfig::Degenerate: return "degenerate";
        case TwoViewConfig::Uncalibrated: return "uncalibrated";
        case TwoViewConfig::PlanarOrPanoramic: return "planar/panoramic";
        default: return "undefined";
    }
}

struct TwoViewOptions {
    // Confidence, trial counts and the rest are COLMAP's. The inlier radius is
    // not: it is 3 px of the image SIFT ran on, not 4 px of the source file
    // (D47). The unit differs from COLMAP's, so the number cannot match it; 3
    // is what a sweep over six scenes at both extraction scales chose.
    RansacOptions ransac;
    int min_num_inliers = 15;
    double max_H_inlier_ratio = 0.8;
    // Optional calibrated pose recovery from the F inliers.
    bool recover_pose = false;
    Mat3 K1 = mat3Identity(), K2 = mat3Identity();
    TwoViewOptions() {
        ransac.max_error = 3.0;
        ransac.min_inlier_ratio = 0.0;
    }
};

struct TwoViewGeometry {
    TwoViewConfig config = TwoViewConfig::Undefined;
    std::vector<char> inlier_mask;  // over the input correspondences
    int num_inliers = 0;
    // In the pixel path these are the fundamental and the homography, in
    // pixels. In the bearing path (D45) they are the same two maps expressed
    // on unit rays -- an essential matrix up to the two extra degrees of
    // freedom the 7-point fit leaves free, and a ray-to-ray homography.
    Mat3 F = {}, H = {};
    bool has_pose = false;
    Pose pose;  // relative pose (camera1 -> camera2), unit translation
};

namespace detail {
struct HModel {
    Mat3 H = mat3Identity(), Hinv = mat3Identity();
};

// Shared tail of both paths: pick F or H by COLMAP's H-inlier-ratio rule and
// copy the winner's inliers out.
inline void selectTwoViewModel(TwoViewGeometry& g, const RansacReport<Mat3>& fRep,
                               const RansacReport<HModel>& hRep, const TwoViewOptions& opt) {
    g.F = fRep.model;
    g.H = hRep.model.H;
    if (std::max(fRep.num_inliers, hRep.num_inliers) < opt.min_num_inliers) {
        g.config = TwoViewConfig::Degenerate;
        return;
    }
    // A high H-inlier ratio means a planar scene or pure rotation, where the
    // epipolar geometry is degenerate (COLMAP max_H_inlier_ratio = 0.8).
    double hRatio = (double)hRep.num_inliers / std::max(1, fRep.num_inliers);
    if (hRatio > opt.max_H_inlier_ratio) {
        g.config = TwoViewConfig::PlanarOrPanoramic;
        g.inlier_mask = hRep.inlier_mask;
        g.num_inliers = hRep.num_inliers;
    } else {
        g.config = TwoViewConfig::Uncalibrated;
        g.inlier_mask = fRep.inlier_mask;
        g.num_inliers = fRep.num_inliers;
    }
}
}  // namespace detail

// Verify one pair from unit bearings rather than pixels (D45).
//
// Everything the pixel path does, one coordinate system down: the epipolar
// constraint and the plane-induced homography both hold exactly for viewing
// rays whatever the field of view, whereas on raw pixels they hold only for a
// pinhole. `opt.ransac.max_error` is in **radians** here (the caller converts
// its pixel budget with the focal length), and residuals are angular.
//
// The config labels are deliberately the same as the pixel path's, so a
// verified pair means the same thing to matches.bin and to the mapper however
// it was verified.
inline TwoViewGeometry estimateTwoViewBearing(const std::vector<Vec3>& b1,
                                              const std::vector<Vec3>& b2,
                                              const TwoViewOptions& opt) {
    TwoViewGeometry g;
    int n = (int)b1.size();
    if (n < opt.min_num_inliers) return g;

    FitFn<Mat3> fFit = [&](const std::vector<int>& s) { return estimateEpipolar7Bearing(b1, b2, s); };
    FitFn<Mat3> fRefit = [&](const std::vector<int>& s) { return estimateEpipolar8Bearing(b1, b2, s); };
    ResidualFn<Mat3> fRes = [&](const Mat3& E, int i) { return sampsonSqBearing(E, b1[i], b2[i]); };
    RansacReport<Mat3> fRep = loransac<Mat3>(n, 7, fFit, fRefit, fRes, opt.ransac);

    FitFn<detail::HModel> hFit = [&](const std::vector<int>& s) {
        std::vector<detail::HModel> out;
        for (const Mat3& H : estimateHomographyBearing(b1, b2, s)) out.push_back({H, inverse3(H)});
        return out;
    };
    ResidualFn<detail::HModel> hRes = [&](const detail::HModel& m, int i) {
        return angularTransferSq(m.H, m.Hinv, b1[i], b2[i]);
    };
    RansacReport<detail::HModel> hRep = loransac<detail::HModel>(n, 4, hFit, hFit, hRes, opt.ransac);

    detail::selectTwoViewModel(g, fRep, hRep, opt);
    if (g.config == TwoViewConfig::Degenerate || g.config == TwoViewConfig::Undefined) return g;

    if (opt.recover_pose && g.config == TwoViewConfig::Uncalibrated && g.num_inliers >= 5) {
        std::vector<int> idx;
        for (int i = 0; i < n; i++)
            if (g.inlier_mask[i]) idx.push_back(i);
        // Only now impose the essential constraint the 7-point fit left free.
        int cnt = 0;
        g.pose = recoverRelativePoseBearing(b1, b2, idx, projectToEssential(g.F), cnt);
        g.has_pose = cnt > 0;
    }
    return g;
}

// Verify one image pair. p1[k] <-> p2[k] are matched keypoint pixel coords.
inline TwoViewGeometry estimateTwoView(const std::vector<Vec2>& p1, const std::vector<Vec2>& p2,
                                       const TwoViewOptions& opt) {
    TwoViewGeometry g;
    int n = (int)p1.size();
    if (n < opt.min_num_inliers) return g;

    // --- Fundamental ---
    FitFn<Mat3> fFit = [&](const std::vector<int>& s) { return estimateFundamental7(p1, p2, s); };
    FitFn<Mat3> fRefit = [&](const std::vector<int>& s) { return estimateFundamental8(p1, p2, s); };
    ResidualFn<Mat3> fRes = [&](const Mat3& F, int i) { return sampsonSq(F, p1[i], p2[i]); };
    RansacReport<Mat3> fRep = loransac<Mat3>(n, 7, fFit, fRefit, fRes, opt.ransac);

    // --- Homography ---
    FitFn<detail::HModel> hFit = [&](const std::vector<int>& s) {
        std::vector<detail::HModel> out;
        for (const Mat3& H : estimateHomography(p1, p2, s)) out.push_back({H, inverse3(H)});
        return out;
    };
    ResidualFn<detail::HModel> hRes = [&](const detail::HModel& m, int i) {
        return symmetricTransferSq(m.H, m.Hinv, p1[i], p2[i]);
    };
    RansacReport<detail::HModel> hRep = loransac<detail::HModel>(n, 4, hFit, hFit, hRes, opt.ransac);

    detail::selectTwoViewModel(g, fRep, hRep, opt);
    if (g.config == TwoViewConfig::Degenerate || g.config == TwoViewConfig::Undefined) return g;

    // Optional calibrated relative pose from the F inliers.
    if (opt.recover_pose && g.config == TwoViewConfig::Uncalibrated && g.num_inliers >= 5) {
        Mat3 K1inv = inverse3(opt.K1), K2inv = inverse3(opt.K2);
        std::vector<Vec2> n1(n), n2(n);
        std::vector<int> idx;
        for (int i = 0; i < n; i++) {
            Vec3 a = mul(K1inv, Vec3{p1[i].x, p1[i].y, 1});
            Vec3 b = mul(K2inv, Vec3{p2[i].x, p2[i].y, 1});
            n1[i] = {a.x / a.z, a.y / a.z};
            n2[i] = {b.x / b.z, b.y / b.z};
            if (g.inlier_mask[i]) idx.push_back(i);
        }
        Mat3 E = essentialFromFundamental(g.F, opt.K1, opt.K2);
        int cnt = 0;
        g.pose = recoverRelativePose(n1, n2, idx, E, cnt);
        g.has_pose = cnt > 0;
    }
    return g;
}

}  // namespace sfm
