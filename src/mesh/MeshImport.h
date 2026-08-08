#pragma once

/*
 * MeshImport.h
 *
 * The read side of MeshExport.h: a triangle mesh back out of the four formats
 * the pipeline writes (PLY, OBJ, glTF, GLB), into the same MeshData.
 *
 * It exists because the viewers have to open what the mesher produced -- the
 * GUI's mesh viewport and the "create mesh" preview, and (compiled to WASM)
 * the standalone web viewer. It is deliberately a READER FOR REAL FILES, not
 * only for ours: an OBJ from Blender, a vertex-colored PLY from MeshLab and a
 * textured GLB from anywhere else all open, because that is what a user who
 * drops a file on the window means.
 *
 * Scope, so the limits are known before they surprise anyone:
 *   PLY   ascii + binary little/big endian; any property order; float/double
 *         positions and normals; uchar or float vertex colors (red/green/blue
 *         or r/g/b, with or without alpha); s,t / u,v / texture_u,texture_v
 *         UVs; face lists of any size (polygons are fan-triangulated).
 *   OBJ   v / vn / vt / f with any of the v, v/vt, v//vn, v/vt/vn forms,
 *         negative (relative) indices, polygons (fan-triangulated), and the
 *         first map_Kd of the referenced .mtl as the texture.
 *   glTF  the .gltf + external .bin and the single-file .glb, POSITION /
 *   /GLB  NORMAL / COLOR_0 / TEXCOORD_0 / indices over every component type
 *         the spec allows for them, all primitives of all meshes merged into
 *         one, node transforms applied. The first baseColorTexture is the
 *         texture. Draco / sparse accessors / animation are NOT read.
 *
 * Everything is host-side and dependency-free apart from the vendored
 * stb_image (textures) and data/Json.h (glTF).
 *
 * WHY THIS IS NOT THE VIEWER'S READER. `viewer/` has its own mesh reader
 * (a streaming PLY parser in viewer/src/viewer.cpp sharing the splat path's
 * header parser, plus OBJ/glTF in JS), and it stays that way on purpose: it
 * parses multi-gigabyte files through a chunk buffer straight into WebGL
 * upload buffers, where this one loads a whole mesh into a MeshData. Two
 * readers with genuinely different constraints, not an accident -- but they
 * read the SAME files, so a format quirk fixed in one is worth checking
 * against the other.
 */

#include <string>

#include "mesh/MeshExport.h"   // MeshData

namespace meshing {

// True when `path`'s extension is one this reader handles. Says nothing about
// the contents -- a .ply is a container (see mesh_ply_has_faces).
bool is_mesh_path(const std::string& path);

// True when `path` is a PLY whose header declares a non-empty face element,
// i.e. a mesh rather than a Gaussian cloud or an SfM point cloud. False for
// anything unreadable, so a caller can use it as "should I treat this .ply as
// a mesh?" without a second check.
bool ply_is_mesh(const std::string& path);

// Read `path` into `out`. Returns false and fills `error` on failure; `out`
// is left in an unspecified state in that case.
//
// Positions are always filled. Normals are filled when the file has them (the
// caller can generate them otherwise -- see mesh_compute_normals). Colors,
// UVs and the texture are filled only when present.
bool read_mesh(const std::string& path, MeshData& out, std::string& error);

// Area-weighted vertex normals, for a file that carried none. No-op when
// `mesh.N` is already the right size unless `force`.
void mesh_compute_normals(MeshData& mesh, bool force = false);

}  // namespace meshing
