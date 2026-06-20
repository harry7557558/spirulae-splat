/*
 * MeshingHost.cpp
 *
 * Host-side orchestration of the 3DGS meshing pipeline (see Meshing.h):
 *   point cloud -> Delaunay -> marching tetrahedra (bisection on cut edges) ->
 *   manifold-preserving short-edge merge -> binary PLY.
 *
 * The GPU work (occupancy field, bisection) lives in Meshing.cu behind
 * OccupancyEvaluator; everything here is plain multithreaded host C++.
 */

#include "Meshing.h"
#include "Delaunay3D.h"

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>

namespace meshing {

namespace {

using Clock = std::chrono::steady_clock;
static double secs_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// The 6 tet edges as corner-index pairs.
static const int kTetEdge[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

// Marching-tetrahedra triangle table. Case bit i is set when occ[i] < iso
// (the "empty" side). Each triangle is 3 edges; an edge is a corner pair.
// Winding follows Paul Bourke's PolygoniseTri so normals are consistent.
struct MTCase { int ntri; int e[2][3][2]; };
static const MTCase kMT[16] = {
    /*0000*/ {0, {}},
    /*0001*/ {1, {{{0,1},{0,2},{0,3}}}},
    /*0010*/ {1, {{{1,0},{1,3},{1,2}}}},
    /*0011*/ {2, {{{0,3},{0,2},{1,3}}, {{1,3},{1,2},{0,2}}}},
    /*0100*/ {1, {{{2,0},{2,1},{2,3}}}},
    /*0101*/ {2, {{{0,1},{2,3},{0,3}}, {{0,1},{1,2},{2,3}}}},
    /*0110*/ {2, {{{0,1},{1,3},{2,3}}, {{0,1},{0,2},{2,3}}}},
    /*0111*/ {1, {{{3,0},{3,2},{3,1}}}},
    /*1000*/ {1, {{{3,0},{3,2},{3,1}}}},
    /*1001*/ {2, {{{0,1},{1,3},{2,3}}, {{0,1},{0,2},{2,3}}}},
    /*1010*/ {2, {{{0,1},{2,3},{0,3}}, {{0,1},{1,2},{2,3}}}},
    /*1011*/ {1, {{{2,0},{2,1},{2,3}}}},
    /*1100*/ {2, {{{0,3},{0,2},{1,3}}, {{1,3},{1,2},{0,2}}}},
    /*1101*/ {1, {{{1,0},{1,3},{1,2}}}},
    /*1110*/ {1, {{{0,1},{0,2},{0,3}}}},
    /*1111*/ {0, {}},
};

static inline int64_t edge_key(int a, int b, int64_t P) {
    int64_t lo = a < b ? a : b, hi = a < b ? b : a;
    return lo * P + hi;
}

// ---------------------------------------------------------------------------
// Manifold-preserving short-edge merge.
//
// Collapses mesh edges shorter than a local threshold (merge_factor times the
// average incident edge length at the shorter end). A collapse is only applied
// when it satisfies the edge-collapse link condition, so anything locally
// manifold stays manifold.
// ---------------------------------------------------------------------------
struct Mesh {
    std::vector<std::array<double,3>> V;
    std::vector<std::array<int,3>> F;
};

static void merge_vertices(Mesh& mesh, float merge_factor, bool verbose) {
    const int nv = (int)mesh.V.size();
    if (nv == 0 || merge_factor <= 0.0f) return;

    std::vector<std::array<double,3>>& V = mesh.V;
    std::vector<std::unordered_set<int>> adj(nv);
    std::vector<std::unordered_set<int>> vt(nv);   // incident triangle ids
    std::vector<char> valive(nv, 1);
    std::vector<char> talive(mesh.F.size(), 1);

    auto add_edge = [&](int a, int b) { adj[a].insert(b); adj[b].insert(a); };
    for (int t = 0; t < (int)mesh.F.size(); ++t) {
        auto& f = mesh.F[t];
        if (f[0] == f[1] || f[1] == f[2] || f[0] == f[2]) { talive[t] = 0; continue; }
        add_edge(f[0], f[1]); add_edge(f[1], f[2]); add_edge(f[2], f[0]);
        vt[f[0]].insert(t); vt[f[1]].insert(t); vt[f[2]].insert(t);
    }

    auto dist = [&](int a, int b) {
        double dx = V[a][0]-V[b][0], dy = V[a][1]-V[b][1], dz = V[a][2]-V[b][2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    // per-vertex average incident edge length -> local feature size
    std::vector<double> Lavg(nv, 0.0);
    for (int v = 0; v < nv; ++v) {
        if (adj[v].empty()) continue;
        double s = 0; for (int w : adj[v]) s += dist(v, w);
        Lavg[v] = s / adj[v].size();
    }

    struct Cand { double len; int u, v; };
    std::vector<Cand> cands;
    for (int u = 0; u < nv; ++u)
        for (int w : adj[u]) {
            if (w <= u) continue;
            double l = dist(u, w);
            double thr = merge_factor * std::min(Lavg[u], Lavg[w]);
            if (l < thr) cands.push_back({l, u, w});
        }
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b){ return a.len < b.len; });

    // third vertices of triangles sharing edge (u,v)
    auto edge_opposites = [&](int u, int v, std::vector<int>& out) {
        out.clear();
        for (int t : vt[u]) {
            if (!talive[t]) continue;
            auto& f = mesh.F[t];
            bool hu=false, hv=false; int other=-1;
            for (int k=0;k<3;k++){ if(f[k]==u)hu=true; else if(f[k]==v)hv=true; else other=f[k]; }
            if (hu && hv) out.push_back(other);
        }
    };

    int collapses = 0;
    for (const Cand& c : cands) {
        int u = c.u, v = c.v;
        if (!valive[u] || !valive[v]) continue;
        if (adj[u].find(v) == adj[u].end()) continue;     // edge gone
        double thr = merge_factor * std::min(Lavg[u], Lavg[v]);
        if (dist(u, v) >= thr) continue;                  // no longer short

        // link condition: common neighbors == edge opposites
        std::vector<int> opp; edge_opposites(u, v, opp);
        std::unordered_set<int> oppset(opp.begin(), opp.end());
        bool ok = true;
        int common = 0;
        for (int w : adj[u]) {
            if (w == v) continue;
            if (adj[v].find(w) != adj[v].end()) {
                ++common;
                if (oppset.find(w) == oppset.end()) { ok = false; break; }
            }
        }
        if (!ok || common != (int)oppset.size()) continue;

        // ---- collapse v -> u (move u to midpoint) ----
        V[u][0] = 0.5*(V[u][0]+V[v][0]);
        V[u][1] = 0.5*(V[u][1]+V[v][1]);
        V[u][2] = 0.5*(V[u][2]+V[v][2]);

        // retarget triangles incident to v
        std::vector<int> vtris(vt[v].begin(), vt[v].end());
        for (int t : vtris) {
            if (!talive[t]) continue;
            auto& f = mesh.F[t];
            for (int k=0;k<3;k++) if (f[k]==v) f[k]=u;
            if (f[0]==f[1] || f[1]==f[2] || f[0]==f[2]) {
                talive[t] = 0;                  // degenerate -> drop
            } else {
                vt[u].insert(t);
            }
        }
        // update adjacency
        for (int w : adj[v]) {
            adj[w].erase(v);
            if (w != u) { adj[w].insert(u); adj[u].insert(w); }
        }
        adj[u].erase(v);
        adj[v].clear(); vt[v].clear(); valive[v] = 0;
        ++collapses;
    }

    // ---- compact ----
    std::vector<int> remap(nv, -1);
    std::vector<std::array<double,3>> newV;
    newV.reserve(nv);
    for (int v = 0; v < nv; ++v)
        if (valive[v]) { remap[v] = (int)newV.size(); newV.push_back(V[v]); }
    std::vector<std::array<int,3>> newF;
    newF.reserve(mesh.F.size());
    for (int t = 0; t < (int)mesh.F.size(); ++t) {
        if (!talive[t]) continue;
        auto& f = mesh.F[t];
        int a=remap[f[0]], b=remap[f[1]], c2=remap[f[2]];
        if (a<0||b<0||c2<0||a==b||b==c2||a==c2) continue;
        newF.push_back({a,b,c2});
    }
    mesh.V.swap(newV);
    mesh.F.swap(newF);
    if (verbose)
        printf("[meshing] merge: %d collapses -> %zu verts, %zu faces\n",
               collapses, mesh.V.size(), mesh.F.size());
}

// ---------------------------------------------------------------------------
// Binary little-endian PLY writer (float vertices, uchar/int faces).
// ---------------------------------------------------------------------------
static void write_ply(const Mesh& mesh, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("meshing: cannot open " + path);
    f << "ply\nformat binary_little_endian 1.0\n";
    f << "element vertex " << mesh.V.size() << "\n";
    f << "property float x\nproperty float y\nproperty float z\n";
    f << "element face " << mesh.F.size() << "\n";
    f << "property list uchar int vertex_indices\n";
    f << "end_header\n";
    for (const auto& v : mesh.V) {
        float xyz[3] = {(float)v[0], (float)v[1], (float)v[2]};
        f.write(reinterpret_cast<const char*>(xyz), sizeof(xyz));
    }
    for (const auto& tri : mesh.F) {
        unsigned char n = 3;
        f.write(reinterpret_cast<const char*>(&n), 1);
        int idx[3] = {tri[0], tri[1], tri[2]};
        f.write(reinterpret_cast<const char*>(idx), sizeof(idx));
    }
}

} // namespace

bool generate_mesh(
    const float* means, const float* quats,
    const float* log_scales, const float* logit_opac, int num_splats,
    const float* cam_pos, int num_cameras,
    const MeshingConfig& cfg,
    const std::string& output_path
) {
    const float iso = cfg.iso;
    auto t0 = Clock::now();

    OccupancyEvaluator ev(means, quats, log_scales, logit_opac, num_splats,
                          cam_pos, num_cameras, cfg);

    // ---- 1. point cloud ----
    std::vector<double> pts;
    ev.generate_point_cloud(pts);
    const int P = ev.num_points();
    if (cfg.verbose) printf("[meshing] point cloud: %d points (%.2fs)\n", P, secs_since(t0));

    // ---- 2. Delaunay tetrahedralization ----
    auto t1 = Clock::now();
    auto tri = delaunay3d::compute_delaunay_3d(pts.data(), P, cfg.num_threads, false);
    const int M = tri.nb_cells;
    if (cfg.verbose)
        printf("[meshing] Delaunay: %d tets, %d threads (%.2fs)\n",
               M, tri.num_threads, secs_since(t1));
    if (M == 0) { printf("[meshing] no tetrahedra produced\n"); return false; }

    // ---- 3. occupancy at all vertices ----
    auto t2 = Clock::now();
    std::vector<float> occ(P);
    ev.evaluate(pts.data(), P, occ.data());
    if (cfg.verbose) printf("[meshing] occupancy field (%.2fs)\n", secs_since(t2));

    // ---- 4. collect cut edges ----
    auto t3 = Clock::now();
    std::unordered_map<int64_t, int> edge_id;
    edge_id.reserve((size_t)M * 2);
    std::vector<int32_t> ea, eb;
    std::vector<float> oa, ob;
    auto get_edge = [&](int va, int vb) -> int {
        int64_t key = edge_key(va, vb, P);
        auto it = edge_id.find(key);
        if (it != edge_id.end()) return it->second;
        int id = (int)ea.size();
        edge_id.emplace(key, id);
        ea.push_back(va); eb.push_back(vb);
        oa.push_back(occ[va]); ob.push_back(occ[vb]);
        return id;
    };
    const int* cv = tri.cell_vertices.data();
    for (int t = 0; t < M; ++t) {
        const int* c = cv + 4 * t;
        for (int e = 0; e < 6; ++e) {
            int a = c[kTetEdge[e][0]], b = c[kTetEdge[e][1]];
            bool ea_empty = occ[a] < iso, eb_empty = occ[b] < iso;
            if (ea_empty != eb_empty) get_edge(a, b);
        }
    }
    const int E = (int)ea.size();
    if (cfg.verbose) printf("[meshing] %d cut edges (%.2fs)\n", E, secs_since(t3));
    if (E == 0) { printf("[meshing] isosurface empty (no edges cross %.3f)\n", iso); return false; }

    // ---- 5. bisection: one mesh vertex per cut edge ----
    auto t4 = Clock::now();
    Mesh mesh;
    mesh.V.resize(E);
    {
        std::vector<double> xyz(3 * E);
        ev.bisect_edges(pts.data(), ea.data(), eb.data(), oa.data(), ob.data(), E, xyz.data());
        for (int i = 0; i < E; ++i)
            mesh.V[i] = {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
    }
    if (cfg.verbose) printf("[meshing] bisection (%.2fs)\n", secs_since(t4));

    // ---- 6. emit triangles ----
    auto t5 = Clock::now();
    for (int t = 0; t < M; ++t) {
        const int* c = cv + 4 * t;
        int code = (occ[c[0]] < iso ? 1 : 0) | (occ[c[1]] < iso ? 2 : 0) |
                   (occ[c[2]] < iso ? 4 : 0) | (occ[c[3]] < iso ? 8 : 0);
        const MTCase& mc = kMT[code];
        for (int ti = 0; ti < mc.ntri; ++ti) {
            int vid[3];
            bool ok = true;
            for (int k = 0; k < 3; ++k) {
                int ca = c[mc.e[ti][k][0]], cb = c[mc.e[ti][k][1]];
                auto it = edge_id.find(edge_key(ca, cb, P));
                if (it == edge_id.end()) { ok = false; break; }
                vid[k] = it->second;
            }
            if (ok) mesh.F.push_back({vid[0], vid[1], vid[2]});
        }
    }
    if (cfg.verbose)
        printf("[meshing] marching tets: %zu raw faces (%.2fs)\n", mesh.F.size(), secs_since(t5));

    // ---- 7. manifold-preserving merge ----
    auto t6 = Clock::now();
    merge_vertices(mesh, cfg.merge_factor, cfg.verbose);
    if (cfg.verbose) printf("[meshing] merge (%.2fs)\n", secs_since(t6));

    // ---- 8. write ----
    write_ply(mesh, output_path);
    printf("[meshing] wrote %s: %zu vertices, %zu faces (total %.2fs)\n",
           output_path.c_str(), mesh.V.size(), mesh.F.size(), secs_since(t0));
    return true;
}

} // namespace meshing
