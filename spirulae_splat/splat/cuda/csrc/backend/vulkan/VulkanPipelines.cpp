#include "VulkanPipelines.h"

#include <cstring>
#include <map>
#include <mutex>
#include <string>

namespace backend {
namespace vk {

// Defined by the generated embed TU (spirulae_splat/embed_spirv.py).
struct SpirvBlob {
    const char* name;
    const uint32_t* data;
    size_t size_bytes;
};
extern const SpirvBlob g_spirv_blobs[];
extern const size_t g_spirv_blob_count;

namespace {

std::mutex g_pipeline_mutex;
std::map<std::string, VkShaderModule> g_modules;         // by blob name
std::map<uint32_t, VkPipelineLayout> g_layouts;          // by push size
std::map<std::string, VkPipeline> g_pipelines;           // by full key

const SpirvBlob* find_blob(const char* name) {
    for (size_t i = 0; i < g_spirv_blob_count; i++)
        if (std::strcmp(g_spirv_blobs[i].name, name) == 0)
            return &g_spirv_blobs[i];
    return nullptr;
}

VkShaderModule get_module(const std::string& name) {
    auto it = g_modules.find(name);
    if (it != g_modules.end()) return it->second;
    const SpirvBlob* blob = find_blob(name.c_str());
    if (!blob) {
        set_error("no embedded SPIR-V blob with this entry name (rerun "
                  "slang/build_spirv.py?)", VK_SUCCESS);
        return VK_NULL_HANDLE;
    }
    VkShaderModuleCreateInfo sci{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    sci.codeSize = blob->size_bytes;
    sci.pCode = blob->data;
    VkShaderModule module = VK_NULL_HANDLE;
    VkResult r = vkCreateShaderModule(Context::get().device(), &sci, nullptr,
                                      &module);
    if (r != VK_SUCCESS) {
        set_error("vkCreateShaderModule failed", r);
        return VK_NULL_HANDLE;
    }
    g_modules[name] = module;
    return module;
}

VkPipelineLayout get_layout(uint32_t push_size) {
    auto it = g_layouts.find(push_size);
    if (it != g_layouts.end()) return it->second;
    VkPushConstantRange range{VK_SHADER_STAGE_COMPUTE_BIT, 0, push_size};
    VkPipelineLayoutCreateInfo lci{
        VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    if (push_size > 0) {
        lci.pushConstantRangeCount = 1;
        lci.pPushConstantRanges = &range;
    }
    VkPipelineLayout layout = VK_NULL_HANDLE;
    VkResult r = vkCreatePipelineLayout(Context::get().device(), &lci,
                                        nullptr, &layout);
    if (r != VK_SUCCESS) {
        set_error("vkCreatePipelineLayout failed", r);
        return VK_NULL_HANDLE;
    }
    g_layouts[push_size] = layout;
    return layout;
}

VkPipeline get_pipeline(const std::string& blob_name, const SpecList& spec,
                        uint32_t push_size, VkPipelineLayout* out_layout) {
    std::string key = blob_name;
    key += '|';
    key += std::to_string(push_size);
    for (uint32_t i = 0; i < spec.count; i++) {
        key += '|';
        key += std::to_string(spec.values[i]);
    }

    std::lock_guard<std::mutex> lock(g_pipeline_mutex);
    VkPipelineLayout layout = get_layout(push_size);
    if (layout == VK_NULL_HANDLE) return VK_NULL_HANDLE;
    *out_layout = layout;

    auto it = g_pipelines.find(key);
    if (it != g_pipelines.end()) return it->second;

    VkShaderModule module = get_module(blob_name);
    if (module == VK_NULL_HANDLE) return VK_NULL_HANDLE;

    VkSpecializationMapEntry entries[SpecList::kMax];
    for (uint32_t i = 0; i < spec.count; i++)
        entries[i] = {i, i * 4u, 4u};
    VkSpecializationInfo spec_info{};
    spec_info.mapEntryCount = spec.count;
    spec_info.pMapEntries = entries;
    spec_info.dataSize = spec.count * 4u;
    spec_info.pData = spec.values;

    VkComputePipelineCreateInfo pci{
        VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    pci.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    pci.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    pci.stage.module = module;
    // slangc emits each single-entry blob with the entry renamed to "main".
    pci.stage.pName = "main";
    // Pin the subgroup size wherever the device allows it: with the default
    // VARYING policy (notably Intel/ANV) WaveGetLaneCount() is not a
    // reliable launch-time constant and tid/WaveGetLaneCount() subgroup
    // indexing (radix sort ranking) miscounts. All kernel workgroup sizes
    // are multiples of any required size, as REQUIRE_FULL_SUBGROUPS needs.
    VkPipelineShaderStageRequiredSubgroupSizeCreateInfoEXT sgsz{
        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_REQUIRED_SUBGROUP_SIZE_CREATE_INFO_EXT};
    if (Context::get().caps().required_subgroup_size) {
        sgsz.requiredSubgroupSize =
            Context::get().caps().required_subgroup_size;
        pci.stage.pNext = &sgsz;
        pci.stage.flags |=
            VK_PIPELINE_SHADER_STAGE_CREATE_REQUIRE_FULL_SUBGROUPS_BIT_EXT;
    }
    if (spec.count) pci.stage.pSpecializationInfo = &spec_info;
    pci.layout = layout;

    VkPipeline pipeline = VK_NULL_HANDLE;
    VkResult r = vkCreateComputePipelines(Context::get().device(),
                                          VK_NULL_HANDLE, 1, &pci, nullptr,
                                          &pipeline);
    if (r != VK_SUCCESS) {
        set_error("vkCreateComputePipelines failed", r);
        return VK_NULL_HANDLE;
    }
    g_pipelines[key] = pipeline;
    return pipeline;
}

}  // namespace

bool dispatch(Stream stream, const char* entry_name, const SpecList& spec,
              uint32_t groups_x, uint32_t groups_y, uint32_t groups_z,
              const void* params, uint32_t params_size) {
    Context& ctx = Context::get();
    if (!ctx.ok()) return false;
    if (groups_x == 0 || groups_y == 0 || groups_z == 0) return true;
    if (params_size > ctx.caps().max_push_constants) {
        set_error("dispatch: params exceed maxPushConstantsSize (use a "
                  "params-ring indirection for this kernel)", VK_SUCCESS);
        return false;
    }
    // Vulkan guarantees only 65535 workgroups per dimension; kernels that
    // can exceed it must fold their grid (none do yet — fail loudly).
    if (groups_x > 65535 || groups_y > 65535 || groups_z > 65535) {
        set_error("dispatch: workgroup count exceeds 65535 (kernel needs a "
                  "folded grid)", VK_SUCCESS);
        return false;
    }

    uint32_t push_size = (params_size + 3u) / 4u * 4u;
    VkPipelineLayout layout = VK_NULL_HANDLE;
    VkPipeline pipeline = get_pipeline(entry_name, spec, push_size, &layout);
    if (pipeline == VK_NULL_HANDLE) return false;

    StreamImpl* s = stream_impl(stream);
    VkCommandBuffer cb = stream_begin(s);
    if (cb == VK_NULL_HANDLE) return false;

    vkCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
    if (params_size > 0) {
        if (push_size == params_size) {
            vkCmdPushConstants(cb, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0,
                               params_size, params);
        } else {  // pad the trailing bytes to the 4-byte range size
            unsigned char padded[256] = {};
            std::memcpy(padded, params, params_size);
            vkCmdPushConstants(cb, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0,
                               push_size, padded);
        }
    }
    vkCmdDispatch(cb, groups_x, groups_y, groups_z);
    stream_barrier(cb);
    return true;
}

void pipelines_shutdown() {
    VkDevice dev = Context::get().device();
    std::lock_guard<std::mutex> lock(g_pipeline_mutex);
    for (auto& kv : g_pipelines) vkDestroyPipeline(dev, kv.second, nullptr);
    for (auto& kv : g_layouts)
        vkDestroyPipelineLayout(dev, kv.second, nullptr);
    for (auto& kv : g_modules)
        vkDestroyShaderModule(dev, kv.second, nullptr);
    g_pipelines.clear();
    g_layouts.clear();
    g_modules.clear();
}

}  // namespace vk
}  // namespace backend
