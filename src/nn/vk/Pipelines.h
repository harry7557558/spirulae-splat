#pragma once
// Compute-pipeline cache over the embedded SPIR-V registry.
//
// A kernel is identified by "<module>.<entry>" plus its specialization-constant
// tuple. Pipelines are created lazily and cached: the first dispatch of a
// variant pays the driver's compile, every later one is a hash lookup.
//
// There is exactly ONE pipeline layout in the process: a push-constant range of
// `maxPushConstantsSize` bytes and no descriptor sets. All data reaches a shader
// as 64-bit device addresses inside the pushed params struct (see AGENTS.md,
// "Kernel conventions"), which is what lets a launcher be a plain function of
// Tensors instead of a descriptor-binding ritual.

#include "nn/vk/Context.h"

#include <cassert>
#include <cstdint>
#include <initializer_list>

namespace nn {
namespace vk {

// Specialization constants, in declaration order within the entry's module
// (Slang assigns ids 0..n-1 in that order). A dispatch that only cares about a
// later axis must still supply the earlier ones -- omitting them silently
// leaves them at their defaults, which has cost this codebase's ancestors a
// day of debugging more than once.
struct SpecList {
    static constexpr uint32_t kMax = 12;
    uint32_t values[kMax] = {};
    uint32_t count = 0;

    SpecList() = default;
    SpecList(std::initializer_list<uint32_t> v) {
        assert(v.size() <= kMax && "SpecList: too many specialization constants");
        for (uint32_t x : v) {
            if (count >= kMax) break;
            values[count++] = x;
        }
    }
    bool operator==(const SpecList& o) const {
        if (count != o.count) return false;
        for (uint32_t i = 0; i < count; ++i)
            if (values[i] != o.values[i]) return false;
        return true;
    }
};

class Pipelines {
public:
    static Pipelines& get();

    // Cached pipeline for (entry, spec). Throws with a list of known entries
    // when the name does not resolve -- a typo in a kernel key is otherwise a
    // silent no-op.
    VkPipeline acquire(const char* entry_name, const SpecList& spec);

    VkPipelineLayout layout();

    // Destroys every pipeline, module and layout. Call at shutdown only.
    void shutdown();

    // Names of all entry points found in the embedded modules, for diagnostics.
    void dumpEntries();

private:
    Pipelines() = default;
    ~Pipelines();
    struct Impl;
    Impl* impl_ = nullptr;
    Impl& impl();
};

}  // namespace vk
}  // namespace nn
