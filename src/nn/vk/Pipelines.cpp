#include "nn/vk/Pipelines.h"

#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "nn/vk/EmbeddedSpirv.h"

NN_DECLARE_EMBEDDED_MODULES(nn)

#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace nn {
namespace vk {

namespace {

// Minimal SPIR-V scan for OpEntryPoint names. We only need the names, so the
// walk stops at the first instruction past the entry-point block.
//
//   word 0 : (wordCount << 16) | opcode
//   OpEntryPoint (15): [1] execution model, [2] id, [3..] literal name
std::vector<std::string> scanEntryPoints(const uint32_t* code, size_t words) {
    std::vector<std::string> out;
    if (words < 5 || code[0] != 0x07230203u) return out;
    size_t i = 5;
    while (i < words) {
        uint32_t inst = code[i];
        uint16_t len = (uint16_t)(inst >> 16);
        uint16_t op = (uint16_t)(inst & 0xFFFFu);
        if (len == 0) break;
        if (op == 15 /* OpEntryPoint */ && i + 3 < words) {
            const char* name = reinterpret_cast<const char*>(&code[i + 3]);
            size_t max_bytes = (words - (i + 3)) * 4;
            out.emplace_back(name, ::strnlen(name, max_bytes));
        }
        // Entry points always precede execution modes and debug info; stopping
        // early would be fragile across slangc versions, so walk the whole
        // module -- it is a few thousand words and happens once per process.
        i += len;
    }
    return out;
}

struct CacheKey {
    std::string entry;
    SpecList    spec;
    bool operator==(const CacheKey& o) const {
        return entry == o.entry && spec == o.spec;
    }
};

struct CacheKeyHash {
    size_t operator()(const CacheKey& k) const {
        size_t h = std::hash<std::string>{}(k.entry);
        for (uint32_t i = 0; i < k.spec.count; ++i)
            h = h * 1000003u ^ std::hash<uint32_t>{}(k.spec.values[i]);
        return h ^ (size_t)k.spec.count;
    }
};

}  // namespace

struct Pipelines::Impl {
    std::mutex mu;
    VkPipelineLayout layout = VK_NULL_HANDLE;
    // module name -> {VkShaderModule, entry names}
    struct Module {
        VkShaderModule          handle = VK_NULL_HANDLE;
        const uint32_t*         code = nullptr;
        size_t                  words = 0;
        std::vector<std::string> entries;
    };
    std::unordered_map<std::string, Module> modules;
    std::unordered_map<CacheKey, VkPipeline, CacheKeyHash> pipelines;
    size_t loaded_count = 0;   // registry entries already scanned

    void loadModules();
};

Pipelines& Pipelines::get() {
    static Pipelines inst;
    return inst;
}

Pipelines::~Pipelines() { shutdown(); }

Pipelines::Impl& Pipelines::impl() {
    if (!impl_) impl_ = new Impl();
    return *impl_;
}

// Additive: each library registers its blobs at its own entry point (see
// nn/vk/EmbeddedSpirv.h), so the registry can still grow after the first
// kernel of another library has been looked up. Only the new tail is scanned.
void Pipelines::Impl::loadModules() {
    NN_ENSURE_EMBEDDED_MODULES(nn);
    size_t n = 0;
    const EmbeddedModule* mods = embedded_modules(&n);
    if (n == loaded_count) return;
    for (size_t i = loaded_count; i < n; ++i) {
        Module m;
        m.code = mods[i].code;
        m.words = mods[i].words;
        m.entries = scanEntryPoints(m.code, m.words);
        modules.emplace(mods[i].name, std::move(m));
    }
    NN_LOG_DEBUG("[vk] %zu shader modules embedded (+%zu)\n", n,
                 n - loaded_count);
    loaded_count = n;
}

VkPipelineLayout Pipelines::layout() {
    Impl& s = impl();
    std::lock_guard<std::mutex> lock(s.mu);
    if (s.layout) return s.layout;
    Context& ctx = Context::get();
    VkPushConstantRange pcr{};
    pcr.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    pcr.offset = 0;
    // One layout for every kernel: a range as large as the device allows. A
    // shader may declare a smaller push block than the layout provides.
    pcr.size = std::min<uint32_t>(ctx.maxPushConstantsSize(), 256u);
    VkPipelineLayoutCreateInfo pli{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pli.pushConstantRangeCount = 1;
    pli.pPushConstantRanges = &pcr;
    NN_VK_CHECK(vkCreatePipelineLayout(ctx.device(), &pli, nullptr, &s.layout));
    return s.layout;
}

VkPipeline Pipelines::acquire(const char* entry_name, const SpecList& spec) {
    Impl& s = impl();
    VkPipelineLayout lay = layout();  // takes the lock itself

    std::lock_guard<std::mutex> lock(s.mu);
    s.loadModules();

    CacheKey key{entry_name, spec};
    auto it = s.pipelines.find(key);
    if (it != s.pipelines.end()) return it->second;

    // Split "<module>.<entry>".
    const char* dot = std::strchr(entry_name, '.');
    NN_CHECK(dot != nullptr,
               "kernel key '%s' must be \"<module>.<entry>\"", entry_name);
    std::string mod_name(entry_name, dot - entry_name);
    std::string ep(dot + 1);

    auto mit = s.modules.find(mod_name);
    if (mit == s.modules.end()) {
        std::string known;
        for (auto& kv : s.modules) { known += " "; known += kv.first; }
        fail("no shader module '%s' (kernel '%s'). Known modules:%s", mod_name.c_str(),
             entry_name, known.c_str());
    }
    Impl::Module& mod = mit->second;

    bool found = false;
    for (auto& e : mod.entries)
        if (e == ep) { found = true; break; }
    if (!found) {
        std::string known;
        for (auto& e : mod.entries) { known += " "; known += e; }
        fail("module '%s' has no entry point '%s'. Entries:%s", mod_name.c_str(),
             ep.c_str(), known.c_str());
    }

    Context& ctx = Context::get();
    if (mod.handle == VK_NULL_HANDLE) {
        VkShaderModuleCreateInfo smi{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        smi.codeSize = mod.words * 4;
        smi.pCode = mod.code;
        NN_VK_CHECK(vkCreateShaderModule(ctx.device(), &smi, nullptr, &mod.handle));
    }

    VkSpecializationMapEntry entries[SpecList::kMax];
    for (uint32_t i = 0; i < spec.count; ++i) {
        entries[i].constantID = i;
        entries[i].offset = i * 4;
        entries[i].size = 4;
    }
    VkSpecializationInfo si{};
    si.mapEntryCount = spec.count;
    si.pMapEntries = entries;
    si.dataSize = spec.count * 4;
    si.pData = spec.values;

    VkPipelineShaderStageCreateInfo stage{
        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
    stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    stage.module = mod.handle;
    stage.pName = ep.c_str();
    stage.pSpecializationInfo = spec.count ? &si : nullptr;

    // Pin the subgroup width where the device lets us. Intel/ANV's default is a
    // VARYING SIMD width, which silently breaks any `tid / WaveGetLaneCount()`
    // indexing. We do NOT set REQUIRE_FULL_SUBGROUPS: that would additionally
    // constrain local_size_x to a multiple of the pinned size (VUID-02757), and
    // pinning alone is what fixes the correctness problem.
    VkPipelineShaderStageRequiredSubgroupSizeCreateInfoEXT rss{
        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_REQUIRED_SUBGROUP_SIZE_CREATE_INFO_EXT};
    if (ctx.hasSubgroupSizeControl()) {
        rss.requiredSubgroupSize = ctx.preferredSubgroupSize();
        stage.pNext = &rss;
    }

    VkComputePipelineCreateInfo cpi{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    cpi.stage = stage;
    cpi.layout = lay;

    VkPipeline pipe = VK_NULL_HANDLE;
    double t0 = now_ms();
    NN_VK_CHECK(
        vkCreateComputePipelines(ctx.device(), VK_NULL_HANDLE, 1, &cpi, nullptr, &pipe));
    double dt = now_ms() - t0;
    if (dt > 100.0)
        NN_LOG_DEBUG("[vk] pipeline %s took %.0f ms to compile\n", entry_name, dt);

    s.pipelines.emplace(std::move(key), pipe);
    return pipe;
}

void Pipelines::dumpEntries() {
    Impl& s = impl();
    std::lock_guard<std::mutex> lock(s.mu);
    s.loadModules();
    for (auto& kv : s.modules)
        for (auto& e : kv.second.entries)
            NN_LOG_INFO("  %s.%s\n", kv.first.c_str(), e.c_str());
}

void Pipelines::shutdown() {
    if (!impl_) return;
    if (Context::initialized()) {
        Context& ctx = Context::get();
        vkDeviceWaitIdle(ctx.device());
        for (auto& kv : impl_->pipelines) vkDestroyPipeline(ctx.device(), kv.second, nullptr);
        for (auto& kv : impl_->modules)
            if (kv.second.handle) vkDestroyShaderModule(ctx.device(), kv.second.handle, nullptr);
        if (impl_->layout) vkDestroyPipelineLayout(ctx.device(), impl_->layout, nullptr);
    }
    delete impl_;
    impl_ = nullptr;
}

}  // namespace vk
}  // namespace nn
