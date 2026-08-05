/*
 * MeshUV.h
 *
 * UV-atlas generation + texture baking for the meshing pipeline. No external
 * dependencies.
 *
 * Pipeline (MeshUV.cpp):
 *   1. Chart segmentation: greedy region growing over the face-adjacency
 *      graph, bounded by the angle to the chart's (area-weighted) average
 *      normal. Charts are then reduced to topological disks (connected,
 *      chi = 1, one simple boundary loop) by recursive splitting -- this is
 *      what makes the atlas work for arbitrary topology.
 *   2. Per-chart parameterization: LSCM (least-squares conformal map, Levy
 *      2002) with two pinned boundary vertices, solved by Jacobi-
 *      preconditioned CGNR. Falls back to a best-fit-plane projection when
 *      the conformal solve degenerates. Charts are rescaled so UV area
 *      matches 3D area (uniform texel density across the atlas).
 *   3. Packing: minimal-area rotation per chart (rotating calipers over the
 *      UV convex hull), then skyline bottom-left packing with a gutter
 *      margin, normalized to [0,1]^2.
 *   4. Baking: per-chart rasterization of texel centers (with ~0.75px edge
 *      dilation), one batched color evaluation for all covered texels, then
 *      gutter dilation + background fill so bilinear/mip sampling never
 *      reads garbage.
 */
#ifndef SPIRULA_MESH_UV_H
#define SPIRULA_MESH_UV_H

#include <functional>
#include <vector>

#include "mesh/MeshExport.h"

namespace meshing {

struct UVAtlasConfig {
    int   texture_size = 0;        // square atlas resolution; 0 = auto (from
                                   // the face_weight texel budget, power of
                                   // two clamped to [1024, 8192])
    int   gutter_px = 4;           // min spacing between charts, in texels
    float chart_angle_deg = 60.0f; // max face-normal deviation within a chart
    int   max_chart_faces = 0;     // 0 = auto (clamp(nf/40, 1024, 65536))
    int   num_threads = 0;         // 0 = all hardware threads
    bool  verbose = true;
    // Optional per-face desired texel count ([nf], positive). Texel density
    // then follows the training views' observed resolution (weight = 3D area
    // x (best-view focal/z)^2 = the face's projected pixel area) instead of
    // giving every face an equal share. Empty = 1 per face.
    std::vector<float> face_weight;
};

// Compute the UV atlas for mesh.V/F: rewrites V and F (vertices on chart
// seams are duplicated, so UVs are per-vertex), duplicates N/C alongside, and
// fills mesh.UV with coordinates in [0,1]^2 (origin top-left). Returns the
// chart id per face (used by bake_texture for race-free parallelism).
// cfg is non-const: texture_size == 0 is resolved to the auto-chosen
// resolution (bake_texture must see the same cfg afterwards).
std::vector<int> build_uv_atlas(MeshData& mesh, UVAtlasConfig& cfg);

// Rasterize the atlas and bake mesh.texture (RGB8, cfg.texture_size square).
// `colorize(xyz, n, rgb_out)` evaluates the model color at n 3D points
// (OccupancyEvaluator::colorize). face_chart is build_uv_atlas's return.
void bake_texture(MeshData& mesh, const std::vector<int>& face_chart,
                  const UVAtlasConfig& cfg,
                  const std::function<void(const float*, int, float*)>& colorize);

} // namespace meshing

#endif
