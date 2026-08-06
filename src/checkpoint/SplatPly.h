#pragma once

// Reading a splat.ply back in.
//
// EngineCheckpoint.cpp writes one per checkpoint, in the property layout the
// original 3D Gaussian Splatting release established and every tool since has
// followed: x y z / nx ny nz / f_dc_0..2 / f_rest_0..3K-1 / opacity /
// scale_0..2 / rot_0..3, raw (un-activated) -- log scales, logit opacity,
// (w,x,y,z) quaternion, DC as the SH coefficient (colour - 0.5) / C0. Ours and
// theirs are the same file, so one reader serves both and a model trained
// anywhere opens here.
//
// The values land in exactly the layout set_data_3dgs() takes, with one
// transposition: a PLY stores the SH rest channel-major (all of red, then all
// of green, then blue) and the engine wants it coefficient-major.
//
// Used by `spirula mesh` (which meshes from a checkpoint's PLY) and by the
// GUI's viewer screen.

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace spirula {

// Raw Gaussians, as stored. Nothing here is activated: scales are logs,
// opacities logits, and the caller hands them to the engine unchanged.
struct SplatCloud {
    int64_t num = 0;
    int sh_degree = 0;               // 0 when the file carries no f_rest
    std::vector<float> means;        // [num, 3]
    std::vector<float> quats;        // [num, 4]  (w, x, y, z)
    std::vector<float> scales;       // [num, 3]  log
    std::vector<float> opacities;    // [num, 1]  logit
    std::vector<float> features_dc;  // [num, 3]
    // [num, dim_sh - 1, 3], coefficient-major -- what set_data_3dgs wants.
    // Empty when sh_degree == 0.
    std::vector<float> features_sh;

    int64_t dim_sh() const { return (int64_t)(sh_degree + 1) * (sh_degree + 1); }
};

// True when `path` is a PLY whose vertex element carries the 3DGS properties.
// False for a plain point cloud, a mesh, or a file that is not a PLY at all --
// it answers a question about the file, so it never throws.
bool is_splat_ply(const std::string& path);

// Throws std::runtime_error naming the file and what was wrong with it.
// `want_sh` false skips the f_rest columns, which is most of the file: the
// mesher only needs geometry and DC colour.
SplatCloud read_splat_ply(const std::string& path, bool want_sh = true);

// Resolve what a user pointed at into (splat.ply, run directory): a .ply
// directly, a step-*.ckpt / *.ckpt directory holding one, or a run directory
// whose newest checkpoint has one. The run directory is where config.json is
// looked for; it may not exist, and for a loose .ply it is a guess.
// Throws when nothing under `path` is a splat.ply.
std::pair<std::string, std::string> find_splat_ply(const std::string& path);

}  // namespace spirula
