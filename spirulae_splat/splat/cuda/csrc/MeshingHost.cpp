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
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    std::vector<std::array<float,3>> V;
    std::vector<std::array<int,3>> F;
    std::vector<std::array<unsigned char,3>> C;  // per-vertex RGB (optional)
};

// Lock-free atomic min on a 64-bit slot (CAS loop).
static inline void atomic_min_u64(uint64_t* p, uint64_t val) {
    uint64_t old = *p;
    while (val < old) {
        uint64_t prev = __sync_val_compare_and_swap(p, old, val);
        if (prev == old) break;
        old = prev;
    }
}
// Order-preserving bit key for a non-negative float (priority by edge length).
static inline uint32_t float_key(float f) {
    uint32_t u; std::memcpy(&u, &f, 4); return u;
}

// Parallel manifold-preserving short-edge merge.
//
// Collapses mesh edges shorter than a local threshold (merge_factor times the
// average incident edge length at the shorter end), only when the edge-collapse
// link condition holds (so anything locally manifold stays manifold).
//
// Parallelism is round-based: each round selects a set of collapsible edges
// whose closed 1-rings are pairwise disjoint, then collapses them concurrently.
// Independence is decided by a "claim" pass -- every candidate writes a
// (length, id) priority into the claim slot of each vertex in its closed
// neighborhood via atomic-min; an edge wins iff it owns every slot it touched.
// The global-shortest candidate always wins, so each round makes progress, and
// the winners' disjoint neighborhoods make the concurrent collapses race-free.
// The outcome is deterministic regardless of thread count.
static void merge_vertices(Mesh& mesh, float merge_factor, float max_flip_deg,
                           bool verbose, int num_threads) {
    const int nv = (int)mesh.V.size();
    const long nf = (long)mesh.F.size();
    if (nv == 0 || merge_factor <= 0.0f) return;
#ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif

    std::vector<std::array<float,3>>& V = mesh.V;
    std::vector<std::array<int,3>>& F = mesh.F;
    std::vector<std::unordered_set<int>> adj(nv), vt(nv);
    std::vector<char> valive(nv, 1);
    std::vector<char> talive(nf, 1);

    #pragma omp parallel for schedule(static)
    for (long t = 0; t < nf; ++t) {
        const auto& f = F[t];
        if (f[0]==f[1] || f[1]==f[2] || f[0]==f[2]) talive[t] = 0;
    }

    // adjacency + incident triangles, built lock-free: each thread owns a
    // vertex range and only writes the sets it owns (scans all faces).
    #pragma omp parallel
    {
    #ifdef _OPENMP
        int nt = omp_get_num_threads(), tid = omp_get_thread_num();
    #else
        int nt = 1, tid = 0;
    #endif
        long lo = (long)nv * tid / nt, hi = (long)nv * (tid + 1) / nt;
        for (long t = 0; t < nf; ++t) {
            if (!talive[t]) continue;
            const auto& f = F[t];
            for (int k = 0; k < 3; ++k) {
                int w = f[k];
                if (w >= lo && w < hi) {
                    adj[w].insert(f[(k+1)%3]);
                    adj[w].insert(f[(k+2)%3]);
                    vt[w].insert((int)t);
                }
            }
        }
    }

    auto dist = [&](int a, int b) {
        float dx = V[a][0]-V[b][0], dy = V[a][1]-V[b][1], dz = V[a][2]-V[b][2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    std::vector<float> Lavg(nv, 0.0f);
    #pragma omp parallel for schedule(dynamic, 4096)
    for (int v = 0; v < nv; ++v) {
        if (adj[v].empty()) continue;
        float s = 0; for (int w : adj[v]) s += dist(v, w);
        Lavg[v] = s / (float)adj[v].size();
    }

    // link condition for edge (u,v): common neighbors of u and v must be
    // exactly the third vertices of triangles on that edge (1 boundary, 2
    // interior). Read-only on the current topology.
    auto link_ok = [&](int u, int v) -> bool {
        int opp0 = -1, opp1 = -1, opp_cnt = 0;
        for (int t : vt[u]) {
            if (!talive[t]) continue;
            const auto& f = F[t];
            bool hu=false, hv=false; int other=-1;
            for (int k=0;k<3;k++){ if(f[k]==u)hu=true; else if(f[k]==v)hv=true; else other=f[k]; }
            if (hu && hv) {
                if (opp_cnt == 0) opp0 = other;
                else if (opp_cnt == 1) opp1 = other;
                else return false;   // >2 faces on edge -> already non-manifold
                ++opp_cnt;
            }
        }
        const auto& au = adj[u]; const auto& av = adj[v];
        const auto& small = (au.size() < av.size()) ? au : av;
        const auto& big   = (au.size() < av.size()) ? av : au;
        int common = 0;
        for (int w : small) {
            if (w == u || w == v) continue;
            if (big.find(w) != big.end()) {
                ++common;
                if (w != opp0 && w != opp1) return false;
            }
        }
        return common == opp_cnt;
    };

    // Geometric guard: reject a collapse of edge (u,v) -> midpoint `np` if it
    // would fold the surface or create a degenerate face. For every alive face
    // incident to u or v (other than the two faces on the collapsing edge), the
    // reshaped face's normal must not rotate past max_flip_deg, and must not go
    // (near-)degenerate. Read-only on the current topology (winners' 1-rings are
    // disjoint within a round, so this is race-free in the parallel collapse).
    const float cos_flip = (max_flip_deg > 0.0f && max_flip_deg < 180.0f)
                               ? std::cos(max_flip_deg * 0.017453292519943295f)
                               : -2.0f;  // <= -1 sentinel: guard disabled
    auto tri_normal = [&](const float* a, const float* b, const float* c,
                          float& nx, float& ny, float& nz) {
        float e1x=b[0]-a[0], e1y=b[1]-a[1], e1z=b[2]-a[2];
        float e2x=c[0]-a[0], e2y=c[1]-a[1], e2z=c[2]-a[2];
        nx = e1y*e2z - e1z*e2y; ny = e1z*e2x - e1x*e2z; nz = e1x*e2y - e1y*e2x;
    };
    auto geom_ok = [&](int u, int v, const float* np) -> bool {
        if (cos_flip <= -1.0f) return true;                 // guard disabled
        for (int pass = 0; pass < 2; ++pass) {
            int center = pass == 0 ? u : v;                 // the vertex moving to np
            for (int t : vt[center]) {
                if (!talive[t]) continue;
                const auto& f = F[t];
                bool hu=false, hv=false;
                for (int k=0;k<3;k++){ if(f[k]==u)hu=true; else if(f[k]==v)hv=true; }
                if (hu && hv) continue;                     // face on the collapsed edge
                const float* q0 = (f[0]==center) ? np : V[f[0]].data();
                const float* q1 = (f[1]==center) ? np : V[f[1]].data();
                const float* q2 = (f[2]==center) ? np : V[f[2]].data();
                float ox,oy,oz, nx,ny,nz;
                tri_normal(V[f[0]].data(), V[f[1]].data(), V[f[2]].data(), ox,oy,oz);
                tri_normal(q0, q1, q2, nx, ny, nz);
                float lo = std::sqrt(ox*ox+oy*oy+oz*oz);
                float ln = std::sqrt(nx*nx+ny*ny+nz*nz);
                if (ln < 1e-12f) return false;              // collapse -> degenerate
                if (lo > 1e-12f &&
                    (ox*nx + oy*ny + oz*nz) / (lo*ln) < cos_flip) return false;  // fold
            }
        }
        return true;
    };

    std::vector<uint64_t> claim(nv);
    const uint64_t CLAIM_MAX = ~0ull;
    struct Cand { int u, v; uint32_t lk; };  // lk = length priority key

    // ---- gather the short-edge candidates ONCE (parallel) ----
    // The link condition is (re)checked when a candidate is actually selected,
    // not here -- like the original single pass, this avoids re-scanning the
    // whole mesh every round.
    std::vector<std::vector<Cand>> loc(nthreads);
    #pragma omp parallel
    {
    #ifdef _OPENMP
        int tid = omp_get_thread_num();
    #else
        int tid = 0;
    #endif
        auto& L = loc[tid];
        #pragma omp for schedule(dynamic, 4096)
        for (int u = 0; u < nv; ++u) {
            for (int w : adj[u]) {
                if (w <= u) continue;
                float len = dist(u, w);
                float thr = merge_factor * std::min(Lavg[u], Lavg[w]);
                if (len < thr) L.push_back({u, w, float_key(len)});
            }
        }
    }
    std::vector<Cand> cands;
    { size_t n=0; for (auto& l : loc) n += l.size(); cands.reserve(n); }
    for (auto& l : loc) cands.insert(cands.end(), l.begin(), l.end());

    // ---- process in rounds over a shrinking candidate list ----
    // Each round selects winners whose closed 1-rings are pairwise disjoint
    // (claim/atomic-min by length priority) and collapses them concurrently.
    // A winner is then removed whether it collapsed or was found invalid (dead
    // endpoint / edge gone / no longer short / link condition) -- the global
    // shortest candidate always wins, so every round removes >=1 candidate and
    // the loop terminates.
    long total_collapses = 0;
    std::vector<char> done;
    while (!cands.empty()) {
        const int NC = (int)cands.size();
        done.assign(NC, 0);

        #pragma omp parallel for schedule(static)
        for (int v = 0; v < nv; ++v) claim[v] = CLAIM_MAX;

        #pragma omp parallel for schedule(dynamic, 1024)
        for (int ci = 0; ci < NC; ++ci) {
            int u = cands[ci].u, v = cands[ci].v;
            if (!valive[u] || !valive[v]) continue;
            uint64_t key = ((uint64_t)cands[ci].lk << 32) | (uint32_t)ci;
            atomic_min_u64(&claim[u], key);
            atomic_min_u64(&claim[v], key);
            for (int w : adj[u]) atomic_min_u64(&claim[w], key);
            for (int w : adj[v]) atomic_min_u64(&claim[w], key);
        }

        long round_collapses = 0;
        #pragma omp parallel for schedule(dynamic, 1024) reduction(+:round_collapses)
        for (int ci = 0; ci < NC; ++ci) {
            int u = cands[ci].u, v = cands[ci].v;
            if (!valive[u] || !valive[v]) { done[ci] = 1; continue; }
            uint64_t key = ((uint64_t)cands[ci].lk << 32) | (uint32_t)ci;
            if (claim[u] != key || claim[v] != key) continue;  // not a winner; retry
            bool win = true;
            for (int w : adj[u]) if (claim[w] != key) { win = false; break; }
            if (win) for (int w : adj[v]) if (claim[w] != key) { win = false; break; }
            if (!win) continue;

            // winner: consumed this round regardless of the outcome below
            done[ci] = 1;
            if (adj[u].find(v) == adj[u].end()) continue;           // edge gone
            float len = dist(u, v);
            if (len >= merge_factor * std::min(Lavg[u], Lavg[v])) continue;  // not short now
            if (!link_ok(u, v)) continue;                           // would break manifold
            float np[3] = {0.5f*(V[u][0]+V[v][0]), 0.5f*(V[u][1]+V[v][1]),
                           0.5f*(V[u][2]+V[v][2])};
            if (!geom_ok(u, v, np)) continue;                       // would fold / degenerate

            // collapse v -> u (move u to midpoint); disjoint 1-rings => race-free
            V[u][0] = np[0]; V[u][1] = np[1]; V[u][2] = np[2];
            for (int t : vt[v]) {
                if (!talive[t]) continue;
                auto& f = F[t];
                for (int k=0;k<3;k++) if (f[k]==v) f[k]=u;
                if (f[0]==f[1] || f[1]==f[2] || f[0]==f[2]) talive[t] = 0;
                else vt[u].insert(t);
            }
            for (int w : adj[v]) {
                adj[w].erase(v);
                if (w != u) { adj[w].insert(u); adj[u].insert(w); }
            }
            adj[u].erase(v);
            adj[v].clear(); vt[v].clear(); valive[v] = 0;
            ++round_collapses;
        }
        total_collapses += round_collapses;

        // ---- compact: keep candidates that are neither consumed nor dead ----
        for (auto& l : loc) l.clear();
        #pragma omp parallel
        {
        #ifdef _OPENMP
            int tid = omp_get_thread_num();
        #else
            int tid = 0;
        #endif
            auto& L = loc[tid];
            #pragma omp for schedule(static)
            for (int ci = 0; ci < NC; ++ci) {
                if (done[ci]) continue;
                int u = cands[ci].u, v = cands[ci].v;
                if (valive[u] && valive[v]) L.push_back(cands[ci]);
            }
        }
        cands.clear();
        for (auto& l : loc) cands.insert(cands.end(), l.begin(), l.end());
    }

    // ---- compact ----
    std::vector<int> remap(nv, -1);
    std::vector<std::array<float,3>> newV; newV.reserve(nv);
    for (int v = 0; v < nv; ++v)
        if (valive[v]) { remap[v] = (int)newV.size(); newV.push_back(V[v]); }
    std::vector<std::array<int,3>> newF; newF.reserve(mesh.F.size());
    for (long t = 0; t < nf; ++t) {
        if (!talive[t]) continue;
        auto& f = F[t];
        int a=remap[f[0]], b=remap[f[1]], c2=remap[f[2]];
        if (a<0||b<0||c2<0||a==b||b==c2||a==c2) continue;
        newF.push_back({a,b,c2});
    }
    mesh.V.swap(newV);
    mesh.F.swap(newF);
    if (verbose)
        printf("[meshing] merge: %ld collapses -> %zu verts, %zu faces\n",
               total_collapses, mesh.V.size(), mesh.F.size());
}

// ---------------------------------------------------------------------------
// Binary little-endian PLY writer (float vertices, uchar/int faces).
// ---------------------------------------------------------------------------
static void write_ply(const Mesh& mesh, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("meshing: cannot open " + path);
    const bool has_color = mesh.C.size() == mesh.V.size();
    f << "ply\nformat binary_little_endian 1.0\n";
    f << "element vertex " << mesh.V.size() << "\n";
    f << "property float x\nproperty float y\nproperty float z\n";
    if (has_color)
        f << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
    f << "element face " << mesh.F.size() << "\n";
    f << "property list uchar int vertex_indices\n";
    f << "end_header\n";
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        f.write(reinterpret_cast<const char*>(mesh.V[i].data()), sizeof(float) * 3);
        if (has_color)
            f.write(reinterpret_cast<const char*>(mesh.C[i].data()), 3);
    }
    for (const auto& tri : mesh.F) {
        unsigned char n = 3;
        f.write(reinterpret_cast<const char*>(&n), 1);
        int idx[3] = {tri[0], tri[1], tri[2]};
        f.write(reinterpret_cast<const char*>(idx), sizeof(idx));
    }
}

// ---------------------------------------------------------------------------
// Global orientation: make the winding consistent across the whole (closed,
// manifold => orientable) surface. Flood-fill orientation across edge
// adjacency so neighbours traverse their shared edge in opposite directions,
// then flip each connected component as a whole if the consistent orientation
// disagrees with the per-triangle outward guess from emission (majority vote),
// so normals end up pointing outward.
// ---------------------------------------------------------------------------
static void orient_mesh(Mesh& mesh) {
    const int nf = (int)mesh.F.size();
    if (nf == 0) return;
    auto& F = mesh.F;
    const int64_t P = (int64_t)mesh.V.size();
    auto ekey = [P](int a, int b) -> int64_t {
        int64_t lo = a<b?a:b, hi = a<b?b:a; return lo*P + hi;
    };

    // undirected edge -> the (<=2) triangles on it
    std::unordered_map<int64_t, std::array<int,2>> etri;
    etri.reserve((size_t)nf * 3);
    for (int t = 0; t < nf; ++t) {
        const auto& f = F[t];
        for (int k = 0; k < 3; ++k) {
            auto& slot = etri.emplace(ekey(f[k], f[(k+1)%3]),
                                      std::array<int,2>{-1,-1}).first->second;
            if (slot[0] < 0) slot[0] = t; else if (slot[1] < 0) slot[1] = t;
        }
    }

    // Flood-fill consistent winding. Neighbours are looked up by undirected
    // edge each step (not a precomputed per-corner index) so they stay valid
    // even after a triangle's winding is flipped.
    //
    // BFS (a queue, via an index into `comp`), not DFS. For a clean orientable
    // surface either gives a fully consistent result, but the extracted mesh can
    // still contain small non-orientable defects (e.g. thin twisted sliver
    // strips). Around such a defect SOME edge stays inconsistent no matter what,
    // and the count of those leftover edges is the number of spanning-tree "back
    // edges" closing an odd cycle -- far smaller for BFS's shallow tree than for
    // DFS's long stringy one. BFS thus confines the residual front/back mixing
    // to the immediate neighbourhood of each defect instead of smearing it along
    // a deep traversal path across the surface.
    std::vector<char> vis(nf, 0), flipped(nf, 0);
    std::vector<int> comp;
    for (int seed = 0; seed < nf; ++seed) {
        if (vis[seed]) continue;
        vis[seed] = 1; comp.clear(); comp.push_back(seed);
        for (size_t h = 0; h < comp.size(); ++h) {
            int t = comp[h];
            const auto& f = F[t];
            for (int k = 0; k < 3; ++k) {
                int a = f[k], b = f[(k+1)%3];
                auto it = etri.find(ekey(a, b));
                int nb = (it->second[0] == t) ? it->second[1] : it->second[0];
                if (nb < 0 || vis[nb]) continue;
                auto& g = F[nb];
                bool same = false;  // neighbour traverses (a->b) the same way?
                for (int m = 0; m < 3; ++m)
                    if (g[m] == a && g[(m+1)%3] == b) { same = true; break; }
                if (same) { std::swap(g[1], g[2]); flipped[nb] = 1; }
                vis[nb] = 1; comp.push_back(nb);
            }
        }
        // align the component's global sign with the emission (outward) guess
        long fl = 0; for (int t : comp) fl += flipped[t];
        if (fl * 2 > (long)comp.size())
            for (int t : comp) std::swap(F[t][1], F[t][2]);
    }
}

// ---------------------------------------------------------------------------
// Duplicate-face removal.
//
// Marching tetrahedra over a dense point cloud occasionally emits two coincident
// triangles on the same vertex triple ("doublet" balloons of zero volume). They
// are degenerate junk and, being isolated 2-triangle shells, survive the merge.
// Drop every face whose (sorted) vertex triple was already seen, keeping one;
// the leftover singleton then falls to the floater pass.
// ---------------------------------------------------------------------------
static void dedup_faces(Mesh& mesh, bool verbose) {
    const long nf = (long)mesh.F.size();
    if (nf == 0) return;
    std::vector<std::array<int,4>> key(nf);
    for (long t = 0; t < nf; ++t) {
        int a = mesh.F[t][0], b = mesh.F[t][1], c = mesh.F[t][2];
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        key[t] = {a, b, c, (int)t};
    }
    std::sort(key.begin(), key.end());
    std::vector<char> drop(nf, 0);
    long ndup = 0;
    for (long i = 1; i < nf; ++i)
        if (key[i][0]==key[i-1][0] && key[i][1]==key[i-1][1] && key[i][2]==key[i-1][2]) {
            drop[key[i][3]] = 1; ++ndup;
        }
    if (ndup == 0) return;
    std::vector<std::array<int,3>> newF; newF.reserve(nf - ndup);
    for (long t = 0; t < nf; ++t) if (!drop[t]) newF.push_back(mesh.F[t]);
    mesh.F.swap(newF);
    if (verbose)
        printf("[meshing] dedup faces: removed %ld duplicate faces -> %zu faces\n",
               ndup, mesh.F.size());
}

// ---------------------------------------------------------------------------
// Floater removal: drop small disconnected components.
//
// Components are the connected components of the face graph (faces sharing a
// vertex are connected, via a union-find over vertices). A component is a
// "floater" when its face count is below min_faces; its faces are deleted and
// the mesh compacted (orphan vertices, i.e. vertices no longer used by any face,
// are dropped too). min_faces <= 1 disables the pass.
// ---------------------------------------------------------------------------
static void remove_floaters(Mesh& mesh, int min_faces, bool verbose) {
    const int nv = (int)mesh.V.size();
    const long nf = (long)mesh.F.size();
    if (min_faces <= 1 || nv == 0 || nf == 0) return;

    std::vector<int> parent(nv);
    for (int i = 0; i < nv; ++i) parent[i] = i;
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto unite = [&](int a, int b) { a=find(a); b=find(b); if (a!=b) parent[a]=b; };
    for (long t = 0; t < nf; ++t) {
        const auto& f = mesh.F[t];
        unite(f[0], f[1]); unite(f[1], f[2]);
    }

    std::unordered_map<int, long> comp_faces;        // root -> face count
    for (long t = 0; t < nf; ++t) comp_faces[find(mesh.F[t][0])] += 1;
    const long thresh = min_faces;

    std::vector<int> remap(nv, -1);
    std::vector<std::array<float,3>> newV; newV.reserve(nv);
    std::vector<std::array<int,3>> newF; newF.reserve(mesh.F.size());
    long dropped_f = 0;
    for (long t = 0; t < nf; ++t) {
        const auto& f = mesh.F[t];
        if (comp_faces[find(f[0])] < thresh) { ++dropped_f; continue; }
        int idx[3];
        for (int k = 0; k < 3; ++k) {
            int v = f[k];
            if (remap[v] < 0) { remap[v] = (int)newV.size(); newV.push_back(mesh.V[v]); }
            idx[k] = remap[v];
        }
        newF.push_back({idx[0], idx[1], idx[2]});
    }
    size_t removed_v = (size_t)nv - newV.size();
    mesh.V.swap(newV);
    mesh.F.swap(newF);
    if (verbose)
        printf("[meshing] floaters: removed %ld faces, %zu verts "
               "-> %zu verts, %zu faces\n",
               dropped_f, removed_v, mesh.V.size(), mesh.F.size());
}

// ---------------------------------------------------------------------------
// Hole filling: triangulate small boundary loops.
//
// A boundary edge is one used by a single triangle. Boundary edges are chained
// into loops (each boundary vertex on a clean manifold-with-boundary has exactly
// two boundary edges). A loop with <= max_edges edges is filled: a 3-loop adds
// one triangle; a longer loop gets a centroid fan (one new vertex at the loop
// centroid plus a triangle per boundary edge). Larger loops (open scene / sky /
// large unseen regions) are left open. Winding of the new faces is left to the
// subsequent orient pass. New centroid vertices are colored by the later pass.
// ---------------------------------------------------------------------------
static void fill_holes(Mesh& mesh, int max_edges, bool verbose) {
    const int nv = (int)mesh.V.size();
    const long nf = (long)mesh.F.size();
    if (max_edges <= 0 || nv == 0 || nf == 0) return;
    const int64_t P = (int64_t)nv;
    auto ekey = [P](int a, int b) -> int64_t {
        int64_t lo = a<b?a:b, hi = a<b?b:a; return lo*P + hi;
    };

    // boundary edges = undirected edges incident to exactly one triangle
    std::unordered_map<int64_t,int> ecount;
    ecount.reserve((size_t)nf * 3);
    for (long t = 0; t < nf; ++t) {
        const auto& f = mesh.F[t];
        for (int k = 0; k < 3; ++k) ecount[ekey(f[k], f[(k+1)%3])] += 1;
    }
    // boundary adjacency: vertex -> (neighbor, edge-id) for each boundary edge
    struct BE { int w; int64_t key; };
    std::unordered_map<int, std::vector<BE>> badj;
    std::unordered_map<int64_t,char> used;
    int n_bnd = 0;
    for (long t = 0; t < nf; ++t) {
        const auto& f = mesh.F[t];
        for (int k = 0; k < 3; ++k) {
            int a = f[k], b = f[(k+1)%3];
            int64_t key = ekey(a, b);
            if (ecount[key] != 1) continue;
            if (used.emplace(key, 0).second) {   // first time we see this bnd edge
                badj[a].push_back({b, key});
                badj[b].push_back({a, key});
                ++n_bnd;
            }
        }
    }
    if (n_bnd == 0) return;

    long n_filled = 0, n_tris = 0;
    std::vector<int> loop;
    for (auto& kv : badj) {
        for (BE& e0 : kv.second) {
            if (used[e0.key]) continue;
            // walk the loop starting along edge (start -> e0.w)
            int start = kv.first;
            loop.clear(); loop.push_back(start);
            used[e0.key] = 1;
            int prev = start, cur = e0.w;
            bool closed = false, ok = true;
            while ((int)loop.size() <= n_bnd) {
                if (cur == start) { closed = true; break; }
                loop.push_back(cur);
                // pick the next unused boundary edge out of cur (prefer != prev)
                int nxt = -1; int64_t nkey = -1;
                auto it = badj.find(cur);
                if (it != badj.end())
                    for (BE& e : it->second)
                        if (!used[e.key] && (nxt < 0 || e.w != prev)) { nxt = e.w; nkey = e.key; }
                if (nxt < 0) { ok = false; break; }
                used[nkey] = 1; prev = cur; cur = nxt;
            }
            if (!ok || !closed) continue;
            const int n = (int)loop.size();
            if (n < 3 || n > max_edges) continue;
            // Only fill simple loops. A repeated vertex (a pinch point with >2
            // boundary edges) would make a fan edge shared by >2 triangles, i.e.
            // non-manifold; leave those holes open.
            if (std::unordered_set<int>(loop.begin(), loop.end()).size() != (size_t)n)
                continue;

            // emit a triangle iff it is non-degenerate (skip collinear slivers)
            auto push_tri = [&](int a, int b, int c) {
                const auto& A = mesh.V[a]; const auto& B = mesh.V[b]; const auto& C = mesh.V[c];
                float e1[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]};
                float e2[3] = {C[0]-A[0], C[1]-A[1], C[2]-A[2]};
                float nx = e1[1]*e2[2]-e1[2]*e2[1], ny = e1[2]*e2[0]-e1[0]*e2[2],
                      nz = e1[0]*e2[1]-e1[1]*e2[0];
                if (nx*nx + ny*ny + nz*nz < 1e-24f) return;
                mesh.F.push_back({a, b, c});
                n_tris += 1;
            };
            if (n == 3) {
                push_tri(loop[0], loop[1], loop[2]);
            } else {
                std::array<float,3> c = {0,0,0};
                for (int v : loop) for (int a = 0; a < 3; ++a) c[a] += mesh.V[v][a];
                for (int a = 0; a < 3; ++a) c[a] /= (float)n;
                int ci = (int)mesh.V.size();
                mesh.V.push_back(c);
                for (int i = 0; i < n; ++i) push_tri(ci, loop[i], loop[(i+1)%n]);
            }
            ++n_filled;
        }
    }
    if (verbose)
        printf("[meshing] fill holes: filled %ld loops (<= %d edges), "
               "added %ld faces -> %zu verts, %zu faces\n",
               n_filled, max_edges, n_tris, mesh.V.size(), mesh.F.size());
}

// ---------------------------------------------------------------------------
// Split non-manifold vertices.
//
// A vertex is non-manifold ("bowtie") when its incident triangles form more
// than one edge-connected fan -- two surface sheets touching at a single point.
// Edges may all be manifold (<= 2 triangles) yet such a vertex still breaks the
// surface: the link is several cycles, not one, so orientation cannot be
// propagated consistently across it (the flood-fill in orient_mesh leaves mixed
// front/back faces). We duplicate the vertex into one copy per fan and reassign
// each fan's triangles to its own copy. Geometry is unchanged (copies share the
// position); the surface becomes a clean manifold-with-boundary, so the
// subsequent orient pass yields a single consistent winding per component.
// ---------------------------------------------------------------------------
static void split_nonmanifold_vertices(Mesh& mesh, bool verbose) {
    const int nv0 = (int)mesh.V.size();
    const long nf = (long)mesh.F.size();
    if (nv0 == 0 || nf == 0) return;

    std::vector<std::vector<int>> vt(nv0);
    for (long t = 0; t < nf; ++t)
        for (int k = 0; k < 3; ++k) vt[mesh.F[t][k]].push_back((int)t);

    std::unordered_map<int,int> parent;   // tri -> parent (local UF over incident tris)
    std::unordered_map<int,int> w2rep;    // other endpoint w -> a tri carrying edge (v,w)
    std::unordered_map<int,int> root2vid; // fan root -> assigned vertex id
    auto find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    long n_split = 0, n_new = 0;
    for (int v = 0; v < nv0; ++v) {
        auto& tris = vt[v];
        if (tris.size() < 2) continue;
        // union incident triangles that share an edge (v,w)
        parent.clear(); w2rep.clear();
        for (int t : tris) parent[t] = t;
        for (int t : tris) {
            const auto& f = mesh.F[t];
            for (int k = 0; k < 3; ++k) {
                int w = f[k];
                if (w == v) continue;
                auto it = w2rep.find(w);
                if (it == w2rep.end()) w2rep[w] = t;
                else { int a = find(it->second), b = find(t); if (a != b) parent[a] = b; }
            }
        }
        // first fan keeps v; each other fan gets a fresh copy
        root2vid.clear();
        int first_root = find(tris[0]);
        for (int t : tris) {
            int r = find(t);
            if (r == first_root) continue;
            int vid;
            auto it = root2vid.find(r);
            if (it == root2vid.end()) {
                vid = (int)mesh.V.size();
                mesh.V.push_back(mesh.V[v]);
                root2vid[r] = vid;
                ++n_new;
            } else vid = it->second;
            auto& f = mesh.F[t];
            for (int k = 0; k < 3; ++k) if (f[k] == v) f[k] = vid;
        }
        if (!root2vid.empty()) ++n_split;
    }
    if (verbose)
        printf("[meshing] split non-manifold verts: %ld split, %ld copies added "
               "-> %zu verts\n", n_split, n_new, mesh.V.size());
}

} // namespace

bool generate_mesh(
    const float* means, const float* quats,
    const float* log_scales, const float* logit_opac, const float* features_dc,
    int num_splats,
    const float* cam_pos, int num_cameras,
    const CameraParams& cams,
    const MeshingConfig& cfg,
    const std::string& output_path
) {
    const float iso = cfg.iso;
    auto t0 = Clock::now();

    OccupancyEvaluator ev(means, quats, log_scales, logit_opac, features_dc, num_splats,
                          cam_pos, num_cameras, cams, cfg);

    // ---- debug: dump one camera's occupancy moments and bail out ----
    // SS_MESH_DEBUG_RENDER=<cam_idx> -> /tmp/ss_moments.f32 (float32 [H,W,3]).
    if (const char* dbg = std::getenv("SS_MESH_DEBUG_RENDER")) {
        int cam = std::atoi(dbg);
        std::vector<float> mom; int W = 0, H = 0;
        if (!ev.debug_render_moments(cam, mom, W, H)) {
            printf("[meshing] debug render unavailable (no camera intrinsics)\n");
            return false;
        }
        double m0sum = 0; float m0max = 0, mn = 1e30f, mx = -1e30f;
        size_t seen = 0;
        for (size_t p = 0; p < (size_t)W * H; ++p) {
            float m0 = mom[3*p], me = mom[3*p+1];
            m0sum += m0; m0max = std::max(m0max, m0);
            if (m0 > 0.5f) { mn = std::min(mn, me); mx = std::max(mx, me); ++seen; }
        }
        printf("[meshing] debug moments cam %d: %dx%d, mean m0=%.4f max m0=%.4f; "
               "depth(m0>0.5) in [%.3f, %.3f] over %zu px\n",
               cam, W, H, m0sum / ((double)W * H), m0max, mn, mx, seen);
        FILE* f = std::fopen("/tmp/ss_moments.f32", "wb");
        if (f) {
            int hdr[3] = {H, W, 3};
            std::fwrite(hdr, sizeof(int), 3, f);
            std::fwrite(mom.data(), sizeof(float), mom.size(), f);
            std::fclose(f);
            printf("[meshing] wrote /tmp/ss_moments.f32 (header: H W 3, then floats)\n");
        }
        return true;
    }

    // ---- 1. point cloud (float everywhere except the Delaunay call) ----
    std::vector<float> pts;
    ev.generate_point_cloud(pts);
    const int P = ev.num_points();
    if (cfg.verbose) printf("[meshing] point cloud: %d points (%.2fs)\n", P, secs_since(t0));

    // ---- 2. Delaunay tetrahedralization (needs double precision) ----
    auto t1 = Clock::now();
    std::vector<double> pts_d(pts.begin(), pts.end());
    auto tri = delaunay3d::compute_delaunay_3d(pts_d.data(), P, cfg.num_threads, false);
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
        std::vector<float> xyz(3 * E);
        ev.bisect_edges(pts.data(), ea.data(), eb.data(), oa.data(), ob.data(), E, xyz.data());
        for (int i = 0; i < E; ++i)
            mesh.V[i] = {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
    }
    if (cfg.verbose) printf("[meshing] bisection (%.2fs)\n", secs_since(t4));

    // ---- 6. emit triangles ----
    // Delaunay tets have no consistent signed-volume orientation, so the MT
    // table's winding alone is not globally coherent. Orient every triangle by
    // the field instead: its normal must point toward the EMPTY side (occ<iso),
    // i.e. away from the occupied tet corners. This yields a consistently
    // outward-facing, watertight mesh (good for backface culling / normals).
    auto t5 = Clock::now();
    for (int t = 0; t < M; ++t) {
        const int* c = cv + 4 * t;
        int code = (occ[c[0]] < iso ? 1 : 0) | (occ[c[1]] < iso ? 2 : 0) |
                   (occ[c[2]] < iso ? 4 : 0) | (occ[c[3]] < iso ? 8 : 0);
        const MTCase& mc = kMT[code];
        if (mc.ntri == 0) continue;

        // outward direction = (centroid of empty corners) - (centroid of solid)
        float cin[3] = {0,0,0}, cout[3] = {0,0,0};
        int nin = 0, nout = 0;
        for (int k = 0; k < 4; ++k) {
            const float* pc = &pts[3 * c[k]];
            if (occ[c[k]] < iso) { cout[0]+=pc[0]; cout[1]+=pc[1]; cout[2]+=pc[2]; ++nout; }
            else                 { cin[0]+=pc[0];  cin[1]+=pc[1];  cin[2]+=pc[2];  ++nin; }
        }
        float dout[3];
        for (int a = 0; a < 3; ++a)
            dout[a] = (nout ? cout[a]/nout : 0.f) - (nin ? cin[a]/nin : 0.f);

        for (int ti = 0; ti < mc.ntri; ++ti) {
            int vid[3];
            bool ok = true;
            for (int k = 0; k < 3; ++k) {
                int ca = c[mc.e[ti][k][0]], cb = c[mc.e[ti][k][1]];
                auto it = edge_id.find(edge_key(ca, cb, P));
                if (it == edge_id.end()) { ok = false; break; }
                vid[k] = it->second;
            }
            if (!ok) continue;
            // flip winding if the triangle normal faces the solid side
            const auto& P0 = mesh.V[vid[0]];
            const auto& P1 = mesh.V[vid[1]];
            const auto& P2 = mesh.V[vid[2]];
            float e1[3] = {P1[0]-P0[0], P1[1]-P0[1], P1[2]-P0[2]};
            float e2[3] = {P2[0]-P0[0], P2[1]-P0[1], P2[2]-P0[2]};
            float nx = e1[1]*e2[2] - e1[2]*e2[1];
            float ny = e1[2]*e2[0] - e1[0]*e2[2];
            float nz = e1[0]*e2[1] - e1[1]*e2[0];
            if (nx*dout[0] + ny*dout[1] + nz*dout[2] < 0.0f)
                std::swap(vid[1], vid[2]);
            mesh.F.push_back({vid[0], vid[1], vid[2]});
        }
    }
    if (cfg.verbose)
        printf("[meshing] marching tets: %zu raw faces (%.2fs)\n", mesh.F.size(), secs_since(t5));

    // ---- 6b. cull vertices seen by no training camera (dataset path only) ----
    // A vertex is kept iff some camera both sees it in-frame and has an
    // unobstructed line of sight to it (occlusion tested against the mesh's own
    // triangles via a GPU LBVH; see OccupancyEvaluator::cull_unseen_vertices).
    // Removing a vertex drops every face incident to it. Since this only deletes
    // faces, each surviving edge keeps <= its original triangle count, so the
    // mesh stays manifold-with-boundary: the only new topology is boundary edges
    // (holes). The downstream merge (link condition handles 1-triangle boundary
    // edges) and orient (edge->triangle map allows the -1 "no neighbour" slot)
    // are already boundary-safe, so no further repair is needed.
    if (cfg.cull_unseen && ev.has_render_cameras() && !mesh.V.empty()) {
        auto tc = Clock::now();
        const int nv = (int)mesh.V.size();
        const int nf = (int)mesh.F.size();
        std::vector<unsigned char> vis(nv, 1);
        ev.cull_unseen_vertices(
            &mesh.V[0][0], nv, nf > 0 ? &mesh.F[0][0] : nullptr, nf, vis.data());

        std::vector<int> remap(nv, -1);
        std::vector<std::array<float,3>> newV; newV.reserve(nv);
        for (int v = 0; v < nv; ++v)
            if (vis[v]) { remap[v] = (int)newV.size(); newV.push_back(mesh.V[v]); }
        std::vector<std::array<int,3>> newF; newF.reserve(mesh.F.size());
        long dropped_f = 0;
        for (int t = 0; t < nf; ++t) {
            const auto& f = mesh.F[t];
            if (!vis[f[0]] || !vis[f[1]] || !vis[f[2]]) { ++dropped_f; continue; }
            newF.push_back({remap[f[0]], remap[f[1]], remap[f[2]]});
        }
        size_t removed_v = (size_t)nv - newV.size();
        mesh.V.swap(newV);
        mesh.F.swap(newF);
        if (cfg.verbose)
            printf("[meshing] cull unseen: removed %zu/%d verts, %ld faces "
                   "-> %zu verts, %zu faces (%.2fs)\n",
                   removed_v, nv, dropped_f, mesh.V.size(), mesh.F.size(), secs_since(tc));
    }

    // ---- 7. manifold-preserving merge (with fold/sliver guard) ----
    auto t6 = Clock::now();
    merge_vertices(mesh, cfg.merge_factor, cfg.merge_max_flip_deg, cfg.verbose, cfg.num_threads);
    if (cfg.verbose) printf("[meshing] merge (%.2fs)\n", secs_since(t6));

    // ---- 7a. post-merge cleanup: drop floaters, then fill small holes ----
    // Run after merge (final topology) and before orient so the orient pass
    // gives the new fill faces a consistent, outward winding. Floaters first so
    // we do not waste work filling holes in components about to be deleted.
    // Order matters: split_nonmanifold_vertices runs BEFORE fill_holes.
    //   - Culling deletes faces, which can split a vertex's fan and leave
    //     non-manifold (bowtie) vertices; at those the boundary has >2 edges, so
    //     a "butterfly" hole pinches through one point. Splitting first makes
    //     every boundary vertex have exactly two boundary edges, so the boundary
    //     becomes a set of disjoint simple loops -- a butterfly becomes two
    //     separate simple loops that fill_holes can then fill (instead of being
    //     skipped), and every filled loop is a true boundary circle. Filling a
    //     disk onto a true boundary circle preserves orientability; filling the
    //     bogus loops produced by ambiguous chaining at a bowtie did NOT (it was
    //     the sole source of non-orientable defects -- MT+merge+cull alone are
    //     orientable).
    //   - The split also isolates vertex-only-attached doublets into their own
    //     tiny components, which the final remove_floaters deletes.
    // A second split runs AFTER fill: fill can occasionally close a loop in a way
    // that leaves a single bowtie vertex, so re-splitting guarantees a fully
    // vertex-manifold result (split preserves orientability). The final floater
    // then removes whole isolated components (manifold-safe) left by either split.
    auto tcl = Clock::now();
    dedup_faces(mesh, cfg.verbose);
    remove_floaters(mesh, cfg.floater_min_faces, cfg.verbose);
    split_nonmanifold_vertices(mesh, cfg.verbose);
    fill_holes(mesh, cfg.fill_hole_max_edges, cfg.verbose);
    split_nonmanifold_vertices(mesh, cfg.verbose);
    remove_floaters(mesh, cfg.floater_min_faces, cfg.verbose);
    if (cfg.verbose) printf("[meshing] cleanup (%.2fs)\n", secs_since(tcl));

    // ---- 7b. globally consistent, outward-facing winding ----
    auto t7 = Clock::now();
    orient_mesh(mesh);
    if (cfg.verbose) printf("[meshing] orient (%.2fs)\n", secs_since(t7));

    // ---- 7c. per-vertex color from splats ----
    auto t8 = Clock::now();
    if (!mesh.V.empty()) {
        const int nv = (int)mesh.V.size();
        std::vector<float> rgb(3 * nv);
        ev.colorize(&mesh.V[0][0], nv, rgb.data());
        mesh.C.resize(nv);
        for (int i = 0; i < nv; ++i)
            for (int a = 0; a < 3; ++a)
                mesh.C[i][a] = (unsigned char)std::min(std::max(
                    (int)std::lround(rgb[3*i+a] * 255.0f), 0), 255);
    }
    if (cfg.verbose) printf("[meshing] color (%.2fs)\n", secs_since(t8));

    // ---- 8. write ----
    write_ply(mesh, output_path);
    printf("[meshing] wrote %s: %zu vertices, %zu faces (total %.2fs)\n",
           output_path.c_str(), mesh.V.size(), mesh.F.size(), secs_since(t0));
    return true;
}

} // namespace meshing
