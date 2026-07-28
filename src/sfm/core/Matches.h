// Feature matches and their flat-file interchange format.
//
// A match is a pair of feature indices into two images' FeatureSets. A
// MatchesDatabase collects the two-view match lists for a whole dataset (the
// input to phase-3 verification and the phase-4 correspondence graph). Like
// features.bin this is our own format (D4), self-describing and trivial to
// parse; it is not COLMAP's SQLite database.
#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sfm/core/Camera.h"

namespace sfm {

// One putative correspondence: feature idx1 in image1 <-> idx2 in image2.
struct FeatureMatch {
    uint32_t idx1 = 0, idx2 = 0;
    float distance = 0;   // L2 descriptor distance (not persisted)
};

// All matches between one ordered image pair. When geometric verification has
// run, `matches` holds only the inliers and `config` is the two-view config
// (0 = not verified / raw; otherwise the sfm/geometry/TwoView.h TwoViewConfig).
struct TwoViewMatches {
    uint32_t image1 = 0, image2 = 0;      // indices into MatchesDatabase::images
    int32_t config = 0;                   // 0 = unverified
    std::vector<FeatureMatch> matches;
};

// One image's identity in the match database.
struct ImageEntry {
    std::string name;         // feature-file stem (== image name)
    uint32_t num_features = 0;
};

struct MatchesDatabase {
    std::vector<ImageEntry> images;
    std::vector<TwoViewMatches> pairs;
    // The camera setup verification actually used (D47): which images share
    // intrinsics, what those intrinsics were, and whether the focal was a prior
    // or a guess. Recorded so the mapper does not have to reconstruct it -- its
    // only other option is to re-search the focal on the *verified inliers*,
    // which are biased towards whatever focal produced them. Empty in a file
    // written before this, or by `--no-verify`.
    std::vector<Camera> cameras;         // one per distinct camera id
    std::vector<uint32_t> camera_ids;    // per image, parallel to `images`
    std::vector<uint8_t> focal_prior;    // parallel to `cameras`; 1 = not a guess
    bool hasCameras() const {
        return !cameras.empty() && camera_ids.size() == images.size();
    }
};

// ---- matches.bin --------------------------------------------------------
// Header: char[4] "VKMT", u32 version=3, u32 num_images.
//   per image: u32 name_len, char[name_len], u32 num_features
// Then u32 num_pairs, per pair:
//   u32 image1, u32 image2, i32 config, u32 num_matches,
//   then num_matches*(u32 idx1, u32 idx2)
// Then (v3) the camera section: u32 num_cameras, per camera
//   u32 id, i32 width, i32 height, i32 colmap_model_id, u8 focal_prior,
//   f64 pixel_scale, u32 num_params, f64 params[num_params]
// followed by u32 num_camera_ids and that many u32 (one per image). A v2 file
// stops after the pairs and reads back with no cameras.

inline void writeMatches(const std::string& path, const MatchesDatabase& db) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    uint32_t version = 3, nimg = (uint32_t)db.images.size();
    f.write("VKMT", 4);
    f.write((const char*)&version, 4);
    f.write((const char*)&nimg, 4);
    for (const ImageEntry& im : db.images) {
        uint32_t len = (uint32_t)im.name.size();
        f.write((const char*)&len, 4);
        f.write(im.name.data(), len);
        f.write((const char*)&im.num_features, 4);
    }
    uint32_t npairs = (uint32_t)db.pairs.size();
    f.write((const char*)&npairs, 4);
    for (const TwoViewMatches& p : db.pairs) {
        uint32_t nm = (uint32_t)p.matches.size();
        f.write((const char*)&p.image1, 4);
        f.write((const char*)&p.image2, 4);
        f.write((const char*)&p.config, 4);
        f.write((const char*)&nm, 4);
        for (const FeatureMatch& m : p.matches) {
            f.write((const char*)&m.idx1, 4);
            f.write((const char*)&m.idx2, 4);
        }
    }
    uint32_t ncam = db.hasCameras() ? (uint32_t)db.cameras.size() : 0;
    f.write((const char*)&ncam, 4);
    for (uint32_t c = 0; c < ncam; c++) {
        const Camera& cam = db.cameras[c];
        int32_t model_id = camColmapId(cam.model);
        uint32_t np = (uint32_t)camColmapParams(cam.model);
        uint8_t prior = c < db.focal_prior.size() ? db.focal_prior[c] : 0;
        double params[16] = {0};
        packColmap(cam, params);
        f.write((const char*)&cam.id, 4);
        f.write((const char*)&cam.width, 4);
        f.write((const char*)&cam.height, 4);
        f.write((const char*)&model_id, 4);
        f.write((const char*)&prior, 1);
        f.write((const char*)&cam.pixel_scale, 8);
        f.write((const char*)&np, 4);
        f.write((const char*)params, np * 8);
    }
    if (ncam) {
        uint32_t nid = (uint32_t)db.camera_ids.size();
        f.write((const char*)&nid, 4);
        f.write((const char*)db.camera_ids.data(), (std::streamsize)nid * 4);
    }
}

inline MatchesDatabase readMatches(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    char magic[4];
    f.read(magic, 4);
    if (std::memcmp(magic, "VKMT", 4) != 0) throw std::runtime_error("bad magic in " + path);
    uint32_t version, nimg;
    f.read((char*)&version, 4);
    f.read((char*)&nimg, 4);
    MatchesDatabase db;
    db.images.resize(nimg);
    for (uint32_t i = 0; i < nimg; i++) {
        uint32_t len;
        f.read((char*)&len, 4);
        db.images[i].name.resize(len);
        f.read(&db.images[i].name[0], len);
        f.read((char*)&db.images[i].num_features, 4);
    }
    uint32_t npairs;
    f.read((char*)&npairs, 4);
    db.pairs.resize(npairs);
    for (uint32_t i = 0; i < npairs; i++) {
        TwoViewMatches& p = db.pairs[i];
        uint32_t nm;
        f.read((char*)&p.image1, 4);
        f.read((char*)&p.image2, 4);
        f.read((char*)&p.config, 4);
        f.read((char*)&nm, 4);
        p.matches.resize(nm);
        for (uint32_t j = 0; j < nm; j++) {
            f.read((char*)&p.matches[j].idx1, 4);
            f.read((char*)&p.matches[j].idx2, 4);
        }
    }
    if (version >= 3) {
        uint32_t ncam = 0;
        f.read((char*)&ncam, 4);
        if (f.gcount() != 4) ncam = 0;
        for (uint32_t c = 0; c < ncam; c++) {
            Camera cam;
            int32_t model_id = 0;
            uint32_t np = 0;
            uint8_t prior = 0;
            double params[16] = {0};
            f.read((char*)&cam.id, 4);
            f.read((char*)&cam.width, 4);
            f.read((char*)&cam.height, 4);
            f.read((char*)&model_id, 4);
            f.read((char*)&prior, 1);
            f.read((char*)&cam.pixel_scale, 8);
            f.read((char*)&np, 4);
            if (np > 16) throw std::runtime_error("bad camera in " + path);
            f.read((char*)params, (std::streamsize)np * 8);
            cam.model = camFromColmapId(model_id);
            unpackColmap(cam, params);
            db.cameras.push_back(cam);
            db.focal_prior.push_back(prior);
        }
        if (ncam) {
            uint32_t nid = 0;
            f.read((char*)&nid, 4);
            if (nid == db.images.size()) {
                db.camera_ids.resize(nid);
                f.read((char*)db.camera_ids.data(), (std::streamsize)nid * 4);
            }
        }
        if (!db.hasCameras()) {   // truncated section: no cameras, not bad ones
            db.cameras.clear();
            db.camera_ids.clear();
            db.focal_prior.clear();
        }
    }
    return db;
}

}  // namespace sfm
