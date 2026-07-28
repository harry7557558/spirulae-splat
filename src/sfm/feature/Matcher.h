// Feature matching: a pluggable IFeatureMatcher interface and a GPU
// brute-force implementation (src/sfm/README.md).
//
// The interface is intentionally the whole contract a matcher must satisfy:
// given two FeatureSets, return the surviving putative correspondences
// (ratio-tested and, optionally, cross-checked). Nothing here assumes 128-D
// uint8 above the shader boundary, so a learned matcher (e.g. LightGlue,
// phase 6) can implement the same interface. Pair selection lives separately in
// sfm/feature/Pairing.h.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sfm/core/Features.h"
#include "sfm/core/Matches.h"
#include "sfm/vk/EmbeddedSpirv.h"
#include "sfm/vk/VkContext.h"

namespace sfm {

struct MatchOptions {
    float max_ratio = 0.8f;     // Lowe ratio (best_dist < ratio * second_dist)
    bool cross_check = true;    // keep only mutual nearest neighbours
    uint32_t max_num_matches = 32768;  // cap per pair (0 = unlimited)
    int device = -1;
    // Pairs per GPU submission. One command buffer holds every dispatch in a
    // batch, so the fence round trip is paid once per batch instead of per
    // pair; the batch's results are downloaded in one copy (D22).
    int batch_pairs = 64;
    // VRAM for resident descriptors. 8192 features/image is ~1 MB, so the
    // default holds ~1500 images; past that the block is recycled and images
    // are re-uploaded per batch, which is still far better than per pair.
    size_t descriptor_budget_bytes = 1536ull << 20;
};

struct IFeatureMatcher {
    virtual ~IFeatureMatcher() = default;
    // Putative matches between a and b (indices into a/b keypoints).
    virtual std::vector<FeatureMatch> match(const FeatureSet& a, const FeatureSet& b) = 0;
    virtual const char* name() const = 0;

    // Match pairs[begin..end) in one go, appending one entry per pair to `out`
    // in pair order. GPU matchers override this to amortize uploads and
    // submissions over the batch; the default just loops.
    virtual void matchBatch(const std::vector<FeatureSet>& feats,
                            const std::vector<std::pair<uint32_t, uint32_t>>& pairs, size_t begin,
                            size_t end, std::vector<std::vector<FeatureMatch>>& out) {
        out.resize(end - begin);
        for (size_t k = begin; k < end; k++)
            out[k - begin] = match(feats[pairs[k].first], feats[pairs[k].second]);
    }
};

// GPU brute-force matcher.
//
// Descriptors for many images stay resident in one device buffer and a whole
// batch of pairs is dispatched from a single command buffer, so an exhaustive
// pass uploads each image roughly once rather than once per pair it appears in
// (at 100 images that is 100 MB instead of 9.9 GB). Distances use the hardware
// packed uint8x4 dot product; see the shader for why that is exact.
class BruteForceMatcher : public IFeatureMatcher {
public:
    explicit BruteForceMatcher(const MatchOptions& opt = {}) : opt_(opt) {
        VkContextOptions vo;
        vo.deviceIndex = opt.device;
        vo.needIntDotProduct = true;
        ctx_.init(vo);
    }

    const char* name() const override { return "brute-force"; }

    // Convenience path (selftests, one-off pairs): a two-image session.
    std::vector<FeatureMatch> match(const FeatureSet& a, const FeatureSet& b) override {
        std::vector<FeatureSet> two;  // copies; this path is not the hot one
        two.push_back(a);
        two.push_back(b);
        std::vector<std::pair<uint32_t, uint32_t>> one{{0u, 1u}};
        std::vector<std::vector<FeatureMatch>> out;
        matchBatch(two, one, 0, 1, out);
        return out.empty() ? std::vector<FeatureMatch>() : std::move(out[0]);
    }

    void matchBatch(const std::vector<FeatureSet>& feats,
                    const std::vector<std::pair<uint32_t, uint32_t>>& pairs, size_t begin,
                    size_t end, std::vector<std::vector<FeatureMatch>>& out) override {
        out.assign(end - begin, {});
        if (begin >= end) return;
        for (size_t k = begin; k < end; k++) {
            checkDescriptors(feats[pairs[k].first]);
            checkDescriptors(feats[pairs[k].second]);
        }
        ensureSetup(feats);
        // The buffers are sized for batch_pairs; honour any range the caller
        // asks for by splitting it, rather than making the two settings agree.
        const size_t chunk = (size_t)std::max(1, opt_.batch_pairs);
        for (size_t b = begin; b < end; b += chunk)
            matchChunk(feats, pairs, b, std::min(b + chunk, end), out, begin);
    }

private:
    void matchChunk(const std::vector<FeatureSet>& feats,
                    const std::vector<std::pair<uint32_t, uint32_t>>& pairs, size_t begin,
                    size_t end, std::vector<std::vector<FeatureMatch>>& out, size_t outBase) {
        ensureResident(feats, pairs, begin, end);

        // Lay out this batch's results: each dispatch writes nQuery uint4s, and
        // the whole region comes back in one download. The column side (rB) is
        // only ever read by the cross-check, so without it neither the
        // reduce_cols dispatch nor its result slots are paid for.
        const bool cc = opt_.cross_check;
        struct Slot { uint32_t a, b, na, nb; uint32_t oa, ob; };
        std::vector<Slot> slots;
        slots.reserve(end - begin);
        uint64_t outCount = 0;
        for (size_t k = begin; k < end; k++) {
            uint32_t ia = pairs[k].first, ib = pairs[k].second;
            uint32_t na = feats[ia].count(), nb = feats[ib].count();
            Slot s{ia, ib, na, nb, (uint32_t)outCount, (uint32_t)(outCount + na)};
            if (na == 0 || nb == 0) s.na = s.nb = 0;  // nothing to dispatch
            else outCount += (uint64_t)na + (cc ? nb : 0);
            slots.push_back(s);
        }
        if (outCount == 0) return;
        if (outCount > resultCap_)  // chunking bounds this; a bug if it ever trips
            throw std::runtime_error("match result buffer too small for a chunk");

        VkCommandBuffer cb = ctx_.begin();
        if (normDirty_.first != normDirty_.second) {
            Push p;
            p.u0 = normDirty_.first;
            p.u1 = normDirty_.second - normDirty_.first;
            ctx_.dispatch(cb, "descriptor_norms", (p.u1 + 63) / 64, p);
            ctx_.barrier(cb);
            normDirty_ = {0, 0};
        }
        // One matrix per pair, then a fold of its column candidates. The two
        // share colPartial, so they are separated by a barrier and consecutive
        // pairs are too; each dispatch already fills the GPU, and the fold is
        // tiny next to the matrix.
        for (const Slot& s : slots) {
            if (s.na == 0 || s.nb == 0) continue;
            uint32_t numWG = (s.na + 63) / 64;
            Push p;
            p.u0 = resident_[s.a];
            p.u1 = s.na;
            p.u2 = resident_[s.b];
            p.u3 = s.nb;
            p.u4 = s.oa;
            p.u5 = s.ob;
            p.u6 = numWG;
            if (cc) {
                // colPartial is shared between pairs, so each pair's matrix and
                // column fold must retire before the next pair starts.
                ctx_.dispatch(cb, "match_pair", numWG, p);
                ctx_.barrier(cb);
                ctx_.dispatch(cb, "reduce_cols", (s.nb + 63) / 64, p);
                ctx_.barrier(cb);
            } else {
                // Row-only dispatches write disjoint result regions and read
                // nothing another pair writes: the whole batch runs with no
                // barriers, which is what keeps thousands of small scoring
                // dispatches from serializing on pipeline drains.
                ctx_.dispatch(cb, "match_rows", numWG, p);
            }
        }
        ctx_.submit(cb);

        std::vector<uint32_t> res((size_t)outCount * 4);
        ctx_.download(bResult_, res.data(), (VkDeviceSize)outCount * 16);

        for (size_t k = 0; k < slots.size(); k++) {
            const Slot& s = slots[k];
            if (s.na == 0 || s.nb == 0) continue;
            out[begin - outBase + k] =
                reduce(res.data() + (size_t)s.oa * 4, res.data() + (size_t)s.ob * 4, s.na, s.nb);
        }
    }

    static constexpr int DW = 32;  // uint words per descriptor

    static void checkDescriptors(const FeatureSet& f) {
        if (f.count() == 0) return;
        if (f.dtype != DType::U8 || f.dim != 128)
            throw std::runtime_error("brute-force matcher expects 128-D uint8 descriptors");
    }

    // Lowe ratio + cross-check, exactly as before -- the GPU only changed how
    // the distances are computed, not what counts as a match.
    std::vector<FeatureMatch> reduce(const uint32_t* rA, const uint32_t* rB, uint32_t na,
                                     uint32_t nb) const {
        std::vector<FeatureMatch> out;
        const float r2 = opt_.max_ratio * opt_.max_ratio;
        for (uint32_t i = 0; i < na; i++) {
            uint32_t j = rA[4 * i + 0];
            uint32_t bestD2 = rA[4 * i + 1], secondD2 = rA[4 * i + 2];
            if (j >= nb) continue;
            if (secondD2 != 0xffffffffu && (float)bestD2 >= r2 * (float)secondD2) continue;
            if (opt_.cross_check && rB[4 * j + 0] != i) continue;
            out.push_back({i, j, std::sqrt((float)bestD2)});
        }
        if (opt_.max_num_matches > 0 && out.size() > opt_.max_num_matches) {
            std::partial_sort(out.begin(), out.begin() + opt_.max_num_matches, out.end(),
                              [](const FeatureMatch& x, const FeatureMatch& y) {
                                  return x.distance < y.distance;
                              });
            out.resize(opt_.max_num_matches);
        }
        return out;
    }

    // All three buffers and the descriptor set must exist before the pipelines
    // load: createDescriptors is what builds the pipeline layout the compute
    // pipelines are created against. Everything is sized once, from the whole
    // feature set, so nothing is reallocated underneath a live pipeline.
    void ensureSetup(const std::vector<FeatureSet>& feats) {
        uint64_t total = 0, maxCount = 0;
        for (const FeatureSet& f : feats) {
            total += f.count();
            maxCount = std::max<uint64_t>(maxCount, f.count());
        }
        if (setup_) {
            // Buffers are allocated once and VkContext frees only at destruction
            // (no per-buffer regrow), so a second feature set that needs more
            // room cannot be served. Say so instead of overrunning: one
            // matcher per feature set.
            if (maxCount > setupMaxCount_)
                throw std::runtime_error(
                    "matcher was sized for a smaller feature set; use a fresh BruteForceMatcher");
            return;
        }
        setupMaxCount_ = maxCount;

        // A single batch can reference 2*batch_pairs distinct images, and all of
        // them must be resident at once, so that is the floor -- otherwise a
        // matcher first used on a two-image pair could never serve a real batch.
        const uint64_t batchFloor =
            (uint64_t)std::max(1, opt_.batch_pairs) * 2 * std::max<uint64_t>(maxCount, 1);
        const uint64_t budget = std::max<uint64_t>(opt_.descriptor_budget_bytes / 128, batchFloor);
        descCap_ = (uint32_t)std::min<uint64_t>(std::max<uint64_t>(total, batchFloor), budget);
        // Worst case a batch can ask for: every pair contributing both images'
        // results, all at the largest feature count.
        resultCap_ = (uint32_t)std::max<uint64_t>(
            1, (uint64_t)std::max(1, opt_.batch_pairs) * 2 * std::max<uint64_t>(maxCount, 1));
        // One pair's column candidates: ceil(nA/64) workgroups x nB columns.
        // Reused pair to pair, so it is sized for the largest single pair.
        const uint64_t mc = std::max<uint64_t>(maxCount, 1);
        uint64_t colCap = ((mc + 63) / 64) * mc;

        bDesc_ = ctx_.createBuffer((VkDeviceSize)descCap_ * DW * 4);
        bNorm_ = ctx_.createBuffer((VkDeviceSize)descCap_ * 4);
        bResult_ = ctx_.createBuffer((VkDeviceSize)resultCap_ * 16);
        bCol_ = ctx_.createBuffer((VkDeviceSize)colCap * 4);
        ctx_.createDescriptors({bDesc_.buf, bNorm_.buf, bResult_.buf, bCol_.buf});

        size_t words = 0;
        const uint32_t* code = findSpirv("match", &words);
        if (!code) throw std::runtime_error("match shader not built into this binary");
        ctx_.loadPipelines(code, words * 4,
                           {"match_pair", "match_rows", "reduce_cols", "descriptor_norms"});
        setup_ = true;
    }

    // Make every image this batch touches resident, uploading the newly added
    // ones as one contiguous block. The allocator bumps; when a batch does not
    // fit it recycles the whole block (exhaustive order means the working set
    // is the whole dataset, so partial eviction would buy nothing).
    void ensureResident(const std::vector<FeatureSet>& feats,
                        const std::vector<std::pair<uint32_t, uint32_t>>& pairs, size_t begin,
                        size_t end) {
        const uint32_t capDesc = descCap_;
        std::vector<uint32_t> need;
        for (size_t k = begin; k < end; k++)
            for (uint32_t img : {pairs[k].first, pairs[k].second})
                if (!resident_.count(img)) need.push_back(img);
        std::sort(need.begin(), need.end());
        need.erase(std::unique(need.begin(), need.end()), need.end());
        if (need.empty()) return;

        uint32_t want = 0;
        for (uint32_t img : need) want += feats[img].count();
        if (used_ + want > capDesc) {  // recycle
            resident_.clear();
            used_ = 0;
            normDirty_ = {0, 0};
            need.clear();
            for (size_t k = begin; k < end; k++)
                for (uint32_t img : {pairs[k].first, pairs[k].second}) need.push_back(img);
            std::sort(need.begin(), need.end());
            need.erase(std::unique(need.begin(), need.end()), need.end());
            want = 0;
            for (uint32_t img : need) want += feats[img].count();
            if (want > capDesc)
                throw std::runtime_error(  // batchFloor in ensureSetup rules this out
                    "descriptor budget too small for one batch of pairs");
        }

        // Gather into one staging-friendly blob: the new images occupy a
        // contiguous descriptor range, so this is a single upload and a single
        // norm dispatch instead of one of each per image.
        const uint32_t first = used_;
        std::vector<uint8_t> blob((size_t)want * 128);
        size_t off = 0;
        for (uint32_t img : need) {
            const FeatureSet& f = feats[img];
            resident_[img] = used_;
            used_ += f.count();
            if (f.count()) {
                memcpy(blob.data() + off, f.descriptors.data(), (size_t)f.count() * 128);
                off += (size_t)f.count() * 128;
            }
        }
        ctx_.upload(bDesc_, blob.data(), (VkDeviceSize)want * 128, (VkDeviceSize)first * 128);
        normDirty_ = {first, used_};
    }

    MatchOptions opt_;
    VkContext ctx_;
    GpuBuffer bDesc_, bNorm_, bResult_, bCol_;
    bool setup_ = false;
    std::map<uint32_t, uint32_t> resident_;   // image index -> first descriptor index
    std::pair<uint32_t, uint32_t> normDirty_{0, 0};  // descriptor range awaiting ||d||^2
    uint32_t used_ = 0, descCap_ = 0, resultCap_ = 0;
    uint64_t setupMaxCount_ = 0;  // largest feature count the buffers were sized for
};

}  // namespace sfm
