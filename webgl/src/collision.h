#include "../../spirulae_splat/splat/cuda/csrc/glm/glm/glm.hpp"

#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <cassert>


struct AABB {
    glm::vec3 min;
    glm::vec3 max;

    AABB() {
        min = glm::vec3(std::numeric_limits<float>::infinity());
        max = glm::vec3(-std::numeric_limits<float>::infinity());
    }
    AABB(const glm::vec3 &mi, const glm::vec3 &ma) : min(mi), max(ma) {}

    bool intersects(const AABB &o) const {
        if (max.x < o.min.x || min.x > o.max.x) return false;
        if (max.y < o.min.y || min.y > o.max.y) return false;
        if (max.z < o.min.z || min.z > o.max.z) return false;
        return true;
    }
    bool contains(const AABB &o) const {
        // this contains o if this.min <= o.min and this.max >= o.max (component-wise)
        return (min.x <= o.min.x && min.y <= o.min.y && min.z <= o.min.z &&
                max.x >= o.max.x && max.y >= o.max.y && max.z >= o.max.z);
    }
    bool containsPoint(const glm::vec3 &p) const {
        return (p.x >= min.x && p.x <= max.x &&
                p.y >= min.y && p.y <= max.y &&
                p.z >= min.z && p.z <= max.z);
    }
    glm::vec3 center() const { return (min + max) * 0.5f; }
    glm::vec3 extent() const { return max - min; }
    void expandToInclude(const AABB &o) {
        min = glm::min(min, o.min);
        max = glm::max(max, o.max);
    }
};

// The user's primitive type must follow this interface:
struct Primitive {
    glm::vec3 mean;
    float opac;
    float inv_cov[6];

    AABB aabb;

    float eval(glm::vec3 p) const {
        float x = p.x-mean.x, y = p.y-mean.y, z = p.z-mean.z;
        float v = inv_cov[0] * x*x + inv_cov[1] * y*y + inv_cov[2] * z*z +
            inv_cov[3] * x*y + inv_cov[4] * x*z + inv_cov[5] * y*z;
        return opac * exp(v);
    }

    glm::vec3 integrateGradient(AABB query) const {
        query.min = glm::max(query.min, aabb.min);
        query.max = glm::min(query.max, aabb.max);
        if (query.min.x >= query.max.x || query.min.y >= query.max.y || query.min.z >= query.max.z)
            return glm::vec3(0.0f);
        glm::vec3 center = query.center();
        glm::vec3 extent = query.extent();
        float int_x = (
            eval(glm::vec3(query.max.x, center.y, center.z)) -
            eval(glm::vec3(query.min.x, center.y, center.z))
        ) * extent.y * extent.z;
        float int_y = (
            eval(glm::vec3(center.x, query.max.y, center.z)) -
            eval(glm::vec3(center.x, query.min.y, center.z))
        ) * extent.x * extent.z;
        float int_z = (
            eval(glm::vec3(center.x, center.y, query.max.z)) -
            eval(glm::vec3(center.x, center.y, query.min.z))
        ) * extent.x * extent.y;
        if (!isfinite(int_x + int_y + int_z))
            return glm::vec3(0);
        return glm::vec3(int_x, int_y, int_z);
    }
};

struct BVHNode {
    AABB bounds;
    std::unique_ptr<BVHNode> left;
    std::unique_ptr<BVHNode> right;

    // indices of primitives assigned to this node (i.e., primitives whose AABB
    // is contained in this node but not contained wholly in a single child)
    std::vector<int> assignedPrims;

    // Sum over assignedPrims of integrateGradient( primitive.aabb )
    glm::vec3 enclosedIntegralSum = glm::vec3(0.0f);
};

class BVHIntegral {
    const std::vector<Primitive> &primitives;
    std::unique_ptr<BVHNode> root;
    int maxLeafPrims;
    int maxDepth;

public:
    BVHIntegral(const std::vector<Primitive> &prims, int maxLeafPrims = 8, int maxDepth = 32)
        : primitives(prims), maxLeafPrims(maxLeafPrims), maxDepth(maxDepth)
    {
        std::vector<int> allIdx(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i) allIdx[i] = (int)i;

        // root bounds should cover all primitives
        AABB rootBounds;
        for (const auto &p : primitives) rootBounds.expandToInclude(p.aabb);

        root = std::make_unique<BVHNode>();
        root->bounds = rootBounds;
        buildNode(root.get(), allIdx, 0);
    }

    // Query: returns sum_i integrateGradient_i(queryBox)
    glm::vec3 query(const AABB &queryBox) const {
        return queryNode(root.get(), queryBox);
    }

private:

    void buildNode(BVHNode *node, const std::vector<int> &primIndices, int depth) {
        // Assign primitives: if node is leaf (by count or depth), store all here.
        if (depth >= maxDepth || primIndices.size() <= (size_t)maxLeafPrims) {
            node->assignedPrims = primIndices;
            // compute enclosed integral sum as sum over these primitives of integral over their own AABB
            node->enclosedIntegralSum = glm::vec3(0.0f);
            for (int idx : node->assignedPrims) {
                node->enclosedIntegralSum += primitives[idx].integrateGradient(primitives[idx].aabb);
            }
            return;
        }

        // Otherwise, split node along longest axis at centroid median
        glm::vec3 ext = node->bounds.extent();
        int axis = 0;
        if (ext.y >= ext.x && ext.y >= ext.z) axis = 1;
        else if (ext.z >= ext.x && ext.z >= ext.y) axis = 2;

        // compute centroids of primitives
        std::vector<std::pair<float,int>> centroids;
        centroids.reserve(primIndices.size());
        for (int idx : primIndices) {
            glm::vec3 c = primitives[idx].aabb.center();
            centroids.emplace_back(c[axis], idx);
        }
        // find median
        size_t mid = centroids.size() / 2;
        std::nth_element(centroids.begin(), centroids.begin() + mid, centroids.end(),
                         [](auto &a, auto &b){ return a.first < b.first; });
        float splitPos = centroids[mid].first;

        // create child bounds
        AABB leftB = node->bounds;
        AABB rightB = node->bounds;
        if (axis == 0) {
            leftB.max.x = splitPos;
            rightB.min.x = splitPos;
        } else if (axis == 1) {
            leftB.max.y = splitPos;
            rightB.min.y = splitPos;
        } else {
            leftB.max.z = splitPos;
            rightB.min.z = splitPos;
        }

        // Ensure children have non-empty extents. If splitPos collapses a child, perform equal split at center
        glm::vec3 leftExt = leftB.extent();
        glm::vec3 rightExt = rightB.extent();
        bool leftEmpty = (leftExt.x <= 0 || leftExt.y <= 0 || leftExt.z <= 0);
        bool rightEmpty = (rightExt.x <= 0 || rightExt.y <= 0 || rightExt.z <= 0);
        if (leftEmpty || rightEmpty) {
            // fallback: split at node center
            float centerPos = node->bounds.center()[axis];
            if (axis == 0) { leftB.max.x = centerPos; rightB.min.x = centerPos; }
            else if (axis == 1) { leftB.max.y = centerPos; rightB.min.y = centerPos; }
            else { leftB.max.z = centerPos; rightB.min.z = centerPos; }
        }

        node->left = std::make_unique<BVHNode>();
        node->right = std::make_unique<BVHNode>();
        node->left->bounds = leftB;
        node->right->bounds = rightB;

        // Partition primitives: if primitive.aabb is fully contained in a child, move it there.
        // Otherwise, keep it in current node.
        std::vector<int> leftIdx, rightIdx, stayIdx;
        leftIdx.reserve(primIndices.size()/2);
        rightIdx.reserve(primIndices.size()/2);
        stayIdx.reserve(primIndices.size()/4);

        for (int idx : primIndices) {
            const AABB &pa = primitives[idx].aabb;
            if (leftB.contains(pa)) leftIdx.push_back(idx);
            else if (rightB.contains(pa)) rightIdx.push_back(idx);
            else stayIdx.push_back(idx);
        }

        // store assigned (staying) primitives at this node
        node->assignedPrims = std::move(stayIdx);
        node->enclosedIntegralSum = glm::vec3(0.0f);
        for (int idx : node->assignedPrims) {
            node->enclosedIntegralSum += primitives[idx].integrateGradient(primitives[idx].aabb);
        }

        // If one child got zero primitives and the other got all, avoid infinite recursion:
        // if either child has no prims and node->assignedPrim count == 0, we should stop splitting.
        if ((leftIdx.empty() && rightIdx.empty()) || (leftIdx.size() + rightIdx.size() == primIndices.size() && node->assignedPrims.empty() && (leftIdx.empty() || rightIdx.empty()))) {
            // Make this a leaf (store all primitives here).
            // Combine all indices into node->assignedPrims and compute sum
            node->assignedPrims = primIndices;
            node->left.reset();
            node->right.reset();
            node->enclosedIntegralSum = glm::vec3(0.0f);
            for (int idx : node->assignedPrims) node->enclosedIntegralSum += primitives[idx].integrateGradient(primitives[idx].aabb);
            return;
        }

        // Recurse children with their lists
        if (!leftIdx.empty()) buildNode(node->left.get(), leftIdx, depth+1);
        else node->left.reset(); // no left child

        if (!rightIdx.empty()) buildNode(node->right.get(), rightIdx, depth+1);
        else node->right.reset(); // no right child
    }

    glm::vec3 queryNode(const BVHNode *node, const AABB &query) const {
        if (!node) return glm::vec3(0.0f);
        if (!node->bounds.intersects(query)) return glm::vec3(0.0f);

        // If node bounds are fully inside query, then every primitive whose AABB is contained in node->bounds
        // is fully inside query as well. We stored at this node only the primitives that were assigned to this node,
        // but children also have their own enclosed sums; however, children will be visited only if needed.
        if (query.contains(node->bounds)) {
            // Node bounds subset of query: we can take enclosed sum for this node,
            // BUT we must also add contributions from children (their enclosed sums) because children store primitives assigned to them.
            // Simpler approach: when node bounds subset of query, sum over:
            // - node->enclosedIntegralSum (primitives assigned here)
            // - plus recursively add children's enclosed sums (they will also be fully inside query)
            glm::vec3 sum = node->enclosedIntegralSum;
            if (node->left) sum += queryNode(node->left.get(), query);   // child top-level will see it's contained and return its enclosedSum + descendants
            if (node->right) sum += queryNode(node->right.get(), query);
            return sum;
        }

        // Partial overlap: we must:
        // 1) For every primitive assigned to this node (their AABBs are subset of node->bounds but node->bounds not subset of query),
        //    if primitive.aabb intersects query, add primitive.integrateGradient(query).
        // 2) Recurse into children that intersect query.
        glm::vec3 sum(0.0f);
        for (int idx : node->assignedPrims) {
            const AABB &pa = primitives[idx].aabb;
            if (pa.intersects(query)) {
                // primitive partially or fully intersects query; call the analytical integrator on the query box
                sum += primitives[idx].integrateGradient(query);
            }
        }

        if (node->left) sum += queryNode(node->left.get(), query);
        if (node->right) sum += queryNode(node->right.get(), query);

        return sum;
    }
};
