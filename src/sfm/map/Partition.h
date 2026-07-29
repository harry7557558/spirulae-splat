// View-graph partitioning for the hierarchical mapper (src/sfm/README.md).
//
// The incremental mapper's cost per registration grows with the model it is
// registering into: every global bundle adjustment, every retriangulation pass
// and every filter sweep walks the whole reconstruction. Reconstructing a
// thousand images as one model therefore costs far more than reconstructing ten
// hundred-image models and gluing them together -- provided the cut is made
// where the capture is genuinely weakly connected, and provided the pieces
// overlap enough to be glued.
//
// That is what this file does, with no dependency beyond the standard library.
// COLMAP's SceneClustering asks METIS for the normalized min cut; the same cut
// comes out of the Fiedler vector of the normalized Laplacian, which a power
// iteration on the (sparse) view graph produces in a few milliseconds.
//
// The graph is the verified view graph: one node per image, edge weight = the
// number of inlier matches the pair kept. Nothing here looks at geometry, so
// it runs before any reconstruction exists.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <queue>
#include <vector>

#include "sfm/core/Matches.h"

namespace sfm {

// Undirected weighted graph in CSR form.
struct ViewGraph {
    std::vector<uint32_t> offs;   // n+1
    std::vector<uint32_t> adj;    // neighbours
    std::vector<double> w;        // parallel to adj
    size_t n() const { return offs.empty() ? 0 : offs.size() - 1; }
    double degree(uint32_t i) const {
        double d = 0;
        for (uint32_t k = offs[i]; k < offs[i + 1]; k++) d += w[k];
        return d;
    }
};

inline ViewGraph buildViewGraph(const MatchesDatabase& db) {
    const size_t n = db.images.size();
    ViewGraph g;
    g.offs.assign(n + 1, 0);
    for (const TwoViewMatches& p : db.pairs) {
        if (p.matches.empty() || p.image1 >= n || p.image2 >= n) continue;
        g.offs[p.image1 + 1]++;
        g.offs[p.image2 + 1]++;
    }
    for (size_t i = 1; i <= n; i++) g.offs[i] += g.offs[i - 1];
    g.adj.resize(g.offs[n]);
    g.w.resize(g.offs[n]);
    std::vector<uint32_t> fill(g.offs.begin(), g.offs.end() - 1);
    for (const TwoViewMatches& p : db.pairs) {
        if (p.matches.empty() || p.image1 >= n || p.image2 >= n) continue;
        const double weight = (double)p.matches.size();
        g.adj[fill[p.image1]] = p.image2;
        g.w[fill[p.image1]++] = weight;
        g.adj[fill[p.image2]] = p.image1;
        g.w[fill[p.image2]++] = weight;
    }
    return g;
}

// Connected components of the subgraph induced on `nodes`, largest first.
inline std::vector<std::vector<uint32_t>> connectedComponents(
    const ViewGraph& g, const std::vector<uint32_t>& nodes) {
    std::vector<char> inside(g.n(), 0), seen(g.n(), 0);
    for (uint32_t v : nodes) inside[v] = 1;
    std::vector<std::vector<uint32_t>> out;
    std::vector<uint32_t> stack;
    for (uint32_t s : nodes) {
        if (seen[s]) continue;
        std::vector<uint32_t> comp;
        stack.assign(1, s);
        seen[s] = 1;
        while (!stack.empty()) {
            const uint32_t v = stack.back();
            stack.pop_back();
            comp.push_back(v);
            for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++) {
                const uint32_t u = g.adj[k];
                if (inside[u] && !seen[u]) { seen[u] = 1; stack.push_back(u); }
            }
        }
        out.push_back(std::move(comp));
    }
    std::stable_sort(out.begin(), out.end(),
                     [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                         return a.size() > b.size();
                     });
    return out;
}

namespace detail {

// Fiedler vector of the induced subgraph's normalized Laplacian, by power
// iteration on M = D^-1/2 W D^-1/2 deflated against its top eigenvector
// D^1/2 * 1 (which belongs to eigenvalue 1 on a connected graph). Returns the
// per-node value in `local` index order; empty when it did not converge.
inline std::vector<double> fiedler(const ViewGraph& g, const std::vector<uint32_t>& nodes,
                                   const std::vector<uint32_t>& local_of, int iters = 300) {
    const size_t m = nodes.size();
    std::vector<double> deg(m, 0.0);
    for (size_t i = 0; i < m; i++) {
        const uint32_t v = nodes[i];
        for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++)
            if (local_of[g.adj[k]] != UINT32_MAX) deg[i] += g.w[k];
    }
    std::vector<double> isq(m);
    for (size_t i = 0; i < m; i++) isq[i] = deg[i] > 0 ? 1.0 / std::sqrt(deg[i]) : 0.0;
    // top eigenvector, normalized
    std::vector<double> v0(m);
    double n0 = 0;
    for (size_t i = 0; i < m; i++) { v0[i] = std::sqrt(deg[i]); n0 += v0[i] * v0[i]; }
    if (n0 <= 0) return {};
    n0 = std::sqrt(n0);
    for (double& x : v0) x /= n0;

    // A deterministic, spread-out start: alternating signs cannot be orthogonal
    // to the Fiedler vector by accident the way a constant start is.
    std::vector<double> x(m), y(m);
    for (size_t i = 0; i < m; i++) x[i] = (i % 2 ? -1.0 : 1.0) + 1e-3 * (double)(i % 7);
    auto orthonormalize = [&](std::vector<double>& z) {
        double dot = 0;
        for (size_t i = 0; i < m; i++) dot += z[i] * v0[i];
        double nz = 0;
        for (size_t i = 0; i < m; i++) { z[i] -= dot * v0[i]; nz += z[i] * z[i]; }
        return std::sqrt(nz);
    };
    if (orthonormalize(x) <= 0) return {};
    {
        double nz = 0;
        for (double q : x) nz += q * q;
        nz = std::sqrt(nz);
        for (double& q : x) q /= nz;
    }
    // M has eigenvalues in [-1, 1]; shifting to (M + I)/2 makes the one we want
    // the dominant mode of a positive operator, so the iteration converges to
    // it instead of to the most negative eigenvalue.
    for (int it = 0; it < iters; it++) {
        std::fill(y.begin(), y.end(), 0.0);
        for (size_t i = 0; i < m; i++) {
            const uint32_t v = nodes[i];
            double acc = 0;
            for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++) {
                const uint32_t j = local_of[g.adj[k]];
                if (j != UINT32_MAX) acc += g.w[k] * isq[j] * x[j];
            }
            y[i] = 0.5 * (isq[i] * acc + x[i]);
        }
        double nz = orthonormalize(y);
        if (nz <= 1e-300) return {};
        for (size_t i = 0; i < m; i++) y[i] /= nz;
        double diff = 0;
        for (size_t i = 0; i < m; i++) diff += std::fabs(y[i] - x[i]);
        x.swap(y);
        if (it > 10 && diff < 1e-7 * (double)m) break;
    }
    // Back to the Laplacian's coordinates: the cut is on D^-1/2 x.
    for (size_t i = 0; i < m; i++) x[i] *= isq[i];
    return x;
}

}  // namespace detail

struct PartitionOptions {
    size_t leaf_max_images = 160;  // split until every part is at most this big
    size_t overlap = 30;           // images each part borrows from its sibling
    size_t min_part = 20;          // a part smaller than this is not worth a model
};

// Split `nodes` in two along the normalized cut, then give each side the
// `overlap` nodes of the other side that are best connected to it. Returns
// {} when the split is not worth making.
inline std::vector<std::vector<uint32_t>> bisect(const ViewGraph& g,
                                                 const std::vector<uint32_t>& nodes,
                                                 const PartitionOptions& opt) {
    if (nodes.size() < 2 * opt.min_part) return {};
    std::vector<uint32_t> local_of(g.n(), UINT32_MAX);
    for (size_t i = 0; i < nodes.size(); i++) local_of[nodes[i]] = (uint32_t)i;

    std::vector<double> f = detail::fiedler(g, nodes, local_of);
    std::vector<uint32_t> order(nodes.size());
    std::iota(order.begin(), order.end(), 0u);
    if (f.empty()) {
        // No Fiedler vector (disconnected or degenerate): fall back to the
        // input order, which for a video capture is the capture order and is
        // exactly the cut one would draw by hand.
    } else {
        std::stable_sort(order.begin(), order.end(),
                         [&](uint32_t a, uint32_t b) { return f[a] < f[b]; });
    }

    // Sweep the sorted order for the best normalized cut, keeping both sides
    // above min_part. cut(S) / (vol(S) * vol(V\S)) is the usual objective; the
    // running sums make the whole sweep linear in the edge count.
    const size_t m = nodes.size();
    std::vector<uint32_t> rank(m);
    for (size_t i = 0; i < m; i++) rank[order[i]] = (uint32_t)i;
    double total_vol = 0;
    std::vector<double> deg(m, 0.0);
    for (size_t i = 0; i < m; i++) {
        const uint32_t v = nodes[i];
        for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++)
            if (local_of[g.adj[k]] != UINT32_MAX) deg[i] += g.w[k];
        total_vol += deg[i];
    }
    double vol = 0, cut = 0, best = 1e300;
    size_t best_split = 0;
    for (size_t i = 0; i + 1 < m; i++) {
        const uint32_t li = order[i];
        const uint32_t v = nodes[li];
        vol += deg[li];
        for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++) {
            const uint32_t j = local_of[g.adj[k]];
            if (j == UINT32_MAX) continue;
            cut += rank[j] <= i ? -g.w[k] : g.w[k];  // moved from cut to inside, or added
        }
        if (i + 1 < opt.min_part || m - i - 1 < opt.min_part) continue;
        const double other = total_vol - vol;
        if (vol <= 0 || other <= 0) continue;
        const double score = cut * (1.0 / vol + 1.0 / other);
        if (score < best) { best = score; best_split = i + 1; }
    }
    if (best_split == 0) return {};

    std::vector<std::vector<uint32_t>> parts(2);
    std::vector<char> side(m, 1);
    for (size_t i = 0; i < m; i++) {
        const bool left = rank[i] < best_split;
        side[i] = left ? 0 : 1;
        parts[left ? 0 : 1].push_back(nodes[i]);
    }

    // Overlap: each side takes the other side's nodes with the strongest
    // connection across the cut. That is what a later Sim(3) merge aligns on,
    // and it is also what lets a cluster triangulate the structure at its own
    // boundary instead of ending in a fringe of one-view features.
    for (int s = 0; s < 2; s++) {
        std::vector<std::pair<double, uint32_t>> cross;
        for (size_t i = 0; i < m; i++) {
            if (side[i] == s) continue;
            double link = 0;
            const uint32_t v = nodes[i];
            for (uint32_t k = g.offs[v]; k < g.offs[v + 1]; k++) {
                const uint32_t j = local_of[g.adj[k]];
                if (j != UINT32_MAX && side[j] == s) link += g.w[k];
            }
            if (link > 0) cross.emplace_back(link, v);
        }
        std::sort(cross.begin(), cross.end(), [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;
        });
        for (size_t k = 0; k < cross.size() && k < opt.overlap; k++)
            parts[s].push_back(cross[k].second);
        std::sort(parts[s].begin(), parts[s].end());
    }
    return parts;
}

// Recursive bisection down to leaves of at most `leaf_max_images` (before
// overlap). Disconnected inputs are separated first: a cut inside a component
// is meaningful, a cut between two is free and says nothing.
inline std::vector<std::vector<uint32_t>> partitionViewGraph(const ViewGraph& g,
                                                             const PartitionOptions& opt) {
    std::vector<uint32_t> all(g.n());
    std::iota(all.begin(), all.end(), 0u);
    std::vector<std::vector<uint32_t>> queue = connectedComponents(g, all), leaves;
    while (!queue.empty()) {
        std::vector<uint32_t> part = std::move(queue.back());
        queue.pop_back();
        // A part is over the leaf size only if it is over it *without* the
        // overlap it borrowed; otherwise the borrowing itself would force
        // another split, and every split borrows again.
        if (part.size() <= opt.leaf_max_images + opt.overlap) {
            if (part.size() >= 2) leaves.push_back(std::move(part));
            continue;
        }
        std::vector<std::vector<uint32_t>> halves = bisect(g, part, opt);
        if (halves.size() != 2) {
            leaves.push_back(std::move(part));
            continue;
        }
        queue.push_back(std::move(halves[0]));
        queue.push_back(std::move(halves[1]));
    }
    std::stable_sort(leaves.begin(), leaves.end(),
                     [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                         return a.size() > b.size();
                     });
    return leaves;
}

}  // namespace sfm
