/*
 * MeshExport.h
 *
 * Dependency-free triangle-mesh writers for the meshing pipeline:
 *   PLY  (binary little-endian; positions + optional normals/vertex colors)
 *   OBJ  (+ .mtl and .png when textured)
 *   GLTF (.gltf JSON + external .bin, external .png when textured)
 *   GLB  (single binary; texture PNG embedded in the buffer)
 *
 * Color modes and format support:
 *   none    : all formats
 *   vertex  : PLY (uchar RGB), GLTF/GLB (COLOR_0, normalized ubyte4).
 *             OBJ is rejected -- it has no standard vertex-color storage.
 *   texture : OBJ, GLTF, GLB. PLY is rejected -- no standard texture support.
 *
 * PNG encoding uses the vendored stb_image_write (already an image dependency
 * of the engine); everything else is emitted by hand.
 */
#ifndef SPIRULA_MESH_EXPORT_H
#define SPIRULA_MESH_EXPORT_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace meshing {

enum class MeshColorMode { None = 0, Vertex = 1, Texture = 2 };

// Central mesh container shared by the pipeline, the UV atlas and the writers.
struct MeshData {
    std::vector<std::array<float, 3>> V;   // positions
    std::vector<std::array<int, 3>>   F;   // triangle vertex indices
    std::vector<std::array<float, 3>> N;   // per-vertex normals (optional)
    std::vector<std::array<unsigned char, 3>> C;  // per-vertex RGB (Vertex mode)

    // Texture mode: per-vertex UVs (seam vertices already duplicated by the
    // atlas), origin at the TOP-LEFT of the image (texel space; the OBJ writer
    // flips v). Plus the baked RGB8 texture, row 0 = top.
    std::vector<std::array<float, 2>> UV;
    int tex_width = 0, tex_height = 0;
    std::vector<unsigned char> texture;     // [tex_height * tex_width * 3]
};

// One requested output: a mesh format plus (texture mode only) the texture
// encoding. Parsed from tokens like "glb", "glb+jpg" (JPEG q95),
// "gltf+jpeg75" (JPEG q75), "obj+png".
struct MeshFormatSpec {
    std::string fmt;        // ply | obj | gltf | glb
    bool jpeg = false;      // texture encoding: false = PNG, true = JPEG
    int  quality = 95;      // JPEG quality (1..100)
    std::string token() const {          // canonical round-trip form
        if (!jpeg) return fmt;
        return fmt + "+jpeg" + std::to_string(quality);
    }
    const char* tex_ext() const { return jpeg ? "jpg" : "png"; }
};

// Parse one "fmt[+png|+jpg|+jpeg[quality]]" token. Throws std::runtime_error
// on unknown format names, bad encodings, or out-of-range quality.
MeshFormatSpec parse_one_mesh_format(const std::string& token);

// Parse a comma-separated format list ("ply,glb+jpg75"); throws on invalid
// tokens or a repeated base format (the outputs would overwrite each other).
// An empty string yields an empty vector.
std::vector<MeshFormatSpec> parse_mesh_formats(const std::string& csv);

// "" when `spec` can represent `mode`, else a human-readable error (also
// rejects a texture encoding suffix when mode is not Texture).
std::string check_export_support(const MeshFormatSpec& spec, MeshColorMode mode);

// Every file write_mesh(spec, mode, base_path) would create -- the mesh file
// first, then its sidecars. Callers use it to check, BEFORE doing any work,
// that the run is not about to write over one of its own inputs: meshing a
// `foo.ply` into a base of `foo` would silently replace the model it was
// asked to mesh, and by the time the write happens the model is already the
// only copy of itself. `check_mesh_outputs_safe` is that check.
std::vector<std::string> mesh_output_paths(const MeshFormatSpec& spec,
                                           MeshColorMode mode,
                                           const std::string& base_path);

// "" when none of the files `specs` would write at `base_path` is one of
// `inputs`, else a human-readable error naming the collision. Paths are
// compared by identity where both exist (so `./a.ply` and `a.ply` collide),
// and lexically otherwise.
std::string check_mesh_outputs_safe(const std::vector<MeshFormatSpec>& specs,
                                    MeshColorMode mode,
                                    const std::string& base_path,
                                    const std::vector<std::string>& inputs);

// Write `mesh` per `spec` at base_path + extension(s). base_path must have no
// extension; sidecar files (.mtl/.bin/.png/.jpg) are placed next to it and
// referenced by basename. Throws std::runtime_error on I/O failure or an
// unsupported format/mode combination.
void write_mesh(const MeshData& mesh, MeshColorMode mode,
                const MeshFormatSpec& spec, const std::string& base_path,
                bool verbose);

// Area-weighted smooth vertex normals from the (consistently oriented) faces.
void compute_vertex_normals(MeshData& mesh);

// Print a mesh-quality report (manifoldness, boundary loops, orientation
// consistency, connected components, min-angle stats, duplicate faces).
void print_mesh_stats(const MeshData& mesh);

} // namespace meshing

#endif
