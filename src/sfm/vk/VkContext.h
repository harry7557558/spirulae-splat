// Minimal Vulkan compute context: device setup, storage buffers with chunked
// staging upload/download, a single shared descriptor set, named compute
// pipelines from one SPIR-V module, and command recording helpers.
#pragma once

#include <vulkan/vulkan.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// The result codes worth naming: the ones a user can act on. Everything else
// prints its number. Out-of-memory in particular used to surface as a bare
// "Vulkan error -2" from a buffer allocation, which says nothing about the
// 3681-image model or the other job sharing the GPU that caused it.
inline const char* vkResultName(int r) {
    switch (r) {
        case -1: return "VK_ERROR_OUT_OF_HOST_MEMORY";
        case -2: return "VK_ERROR_OUT_OF_DEVICE_MEMORY (the GPU ran out of memory -- "
                        "a smaller --quality, a --vram-budget, or not sharing the device "
                        "with another job)";
        case -3: return "VK_ERROR_INITIALIZATION_FAILED";
        case -4: return "VK_ERROR_DEVICE_LOST";
        case -5: return "VK_ERROR_MEMORY_MAP_FAILED";
        case -9: return "VK_ERROR_FEATURE_NOT_PRESENT";
        case -11: return "VK_ERROR_FORMAT_NOT_SUPPORTED";
        default: return "";
    }
}

#define VK_CHECK(x) do { VkResult _r = (x); if (_r != VK_SUCCESS) { \
    fprintf(stderr, "Vulkan error %d %s at %s:%d\n", (int)_r, vkResultName((int)_r), \
            __FILE__, __LINE__); exit(1); } } while (0)

// Must match the push constant block in sfm/shaders/ba/ba.slang.
struct Push {
    uint32_t u0 = 0, u1 = 0, u2 = 0, u3 = 0;
    float f0 = 0, f1 = 0, f2 = 0, f3 = 0;
    // Appended, so the first 32 bytes keep the layout every existing shader
    // declares; the matcher needs more than four integer slots. Vulkan
    // guarantees at least 128 bytes of push constants, and a shader may declare
    // a smaller range than the pipeline layout provides.
    uint32_t u4 = 0, u5 = 0, u6 = 0, u7 = 0;
};

struct GpuBuffer {
    VkBuffer buf = VK_NULL_HANDLE;
    VkDeviceMemory mem = VK_NULL_HANDLE;
    VkDeviceSize size = 0;
};

struct VkContextOptions {
    bool needFloat64 = false;
    bool needFloatAtomics = false;   // VK_EXT_shader_atomic_float (f32 add, and f64 add if needFloat64)
    bool needInt64Atomics = false;   // buffer int64 atomics (emulated-double CAS)
    bool needIntDotProduct = false;  // VK_KHR_shader_integer_dot_product (packed uint8x4, matcher)
    int deviceIndex = -1;            // -1 = prefer discrete GPU
    bool validate = false;
    bool profile = false;            // per-dispatch GPU timestamps, aggregated by pipeline name
};

class VkContext {
public:
    VkContext() = default;
    // Owners hold VkContext by value and hand out GpuBuffers freely; a copy would
    // double-destroy every handle.
    VkContext(const VkContext&) = delete;
    VkContext& operator=(const VkContext&) = delete;

    // Full teardown. Every buffer ever created through this context is freed
    // here (createBufferRaw tracks them), so a scoped VkContext -- the mapper's
    // per-global-BA solver, the prefilter's matcher -- returns its VRAM when
    // it dies instead of leaking it for the life of the process. Before this,
    // a 1363-image run OOMed: each of the mapper's periodic global BAs leaked
    // a device plus problem-sized buffers.
    ~VkContext() {
        if (device_ == VK_NULL_HANDLE) {
            if (instance_ != VK_NULL_HANDLE) vkDestroyInstance(instance_, nullptr);
            return;
        }
        vkDeviceWaitIdle(device_);
        for (auto& p : pipelines_) vkDestroyPipeline(device_, p.second, nullptr);
        if (shaderModule_ != VK_NULL_HANDLE) vkDestroyShaderModule(device_, shaderModule_, nullptr);
        if (descPool_ != VK_NULL_HANDLE) vkDestroyDescriptorPool(device_, descPool_, nullptr);
        if (pipeLayout_ != VK_NULL_HANDLE) vkDestroyPipelineLayout(device_, pipeLayout_, nullptr);
        if (setLayout_ != VK_NULL_HANDLE) vkDestroyDescriptorSetLayout(device_, setLayout_, nullptr);
        if (queryPool_ != VK_NULL_HANDLE) vkDestroyQueryPool(device_, queryPool_, nullptr);
        if (stagingPtr_) vkUnmapMemory(device_, staging_.mem);
        if (stagingDlPtr_) vkUnmapMemory(device_, stagingDl_.mem);
        for (auto& a : allocations_) {
            vkDestroyBuffer(device_, a.first, nullptr);
            vkFreeMemory(device_, a.second, nullptr);
        }
        if (fence_ != VK_NULL_HANDLE) vkDestroyFence(device_, fence_, nullptr);
        if (cmdPool_ != VK_NULL_HANDLE) vkDestroyCommandPool(device_, cmdPool_, nullptr);
        vkDestroyDevice(device_, nullptr);
        vkDestroyInstance(instance_, nullptr);
    }

    void init(const VkContextOptions& opt) {
        // ---- instance ----
        VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
        app.pApplicationName = "vk_ba";
        app.apiVersion = VK_API_VERSION_1_2;
        VkInstanceCreateInfo ici{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
        ici.pApplicationInfo = &app;
        const char* layers[] = {"VK_LAYER_KHRONOS_validation"};
        if (opt.validate) {
            ici.enabledLayerCount = 1;
            ici.ppEnabledLayerNames = layers;
        }
        VK_CHECK(vkCreateInstance(&ici, nullptr, &instance_));

        // ---- physical device ----
        uint32_t n = 0;
        vkEnumeratePhysicalDevices(instance_, &n, nullptr);
        if (n == 0) throw std::runtime_error("no Vulkan devices");
        std::vector<VkPhysicalDevice> devs(n);
        vkEnumeratePhysicalDevices(instance_, &n, devs.data());
        int pick = opt.deviceIndex;
        if (pick < 0) {
            pick = 0;
            for (uint32_t i = 0; i < n; i++) {
                VkPhysicalDeviceProperties p;
                vkGetPhysicalDeviceProperties(devs[i], &p);
                if (p.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) { pick = (int)i; break; }
            }
        }
        phys_ = devs[pick];
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(phys_, &props);
        fprintf(stderr, "[vk] using device: %s\n", props.deviceName);

        // ---- queue family ----
        uint32_t qn = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(phys_, &qn, nullptr);
        std::vector<VkQueueFamilyProperties> qf(qn);
        vkGetPhysicalDeviceQueueFamilyProperties(phys_, &qn, qf.data());
        queueFamily_ = ~0u;
        for (uint32_t i = 0; i < qn; i++)
            if (qf[i].queueFlags & VK_QUEUE_COMPUTE_BIT) { queueFamily_ = i; break; }
        if (queueFamily_ == ~0u) throw std::runtime_error("no compute queue");

        // ---- device ----
        float prio = 1.0f;
        VkDeviceQueueCreateInfo qci{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        qci.queueFamilyIndex = queueFamily_;
        qci.queueCount = 1;
        qci.pQueuePriorities = &prio;

        VkPhysicalDeviceFeatures2 feat2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
        feat2.features.shaderFloat64 = opt.needFloat64 ? VK_TRUE : VK_FALSE;
        feat2.features.shaderInt64 = opt.needInt64Atomics ? VK_TRUE : VK_FALSE;

        VkPhysicalDeviceVulkan12Features f12{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
        f12.shaderBufferInt64Atomics = opt.needInt64Atomics ? VK_TRUE : VK_FALSE;
        feat2.pNext = &f12;

        VkPhysicalDeviceShaderAtomicFloatFeaturesEXT fAtom{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
        std::vector<const char*> exts;
        if (opt.needFloatAtomics) {
            fAtom.shaderBufferFloat32Atomics = VK_TRUE;
            fAtom.shaderBufferFloat32AtomicAdd = VK_TRUE;
            if (opt.needFloat64) {
                fAtom.shaderBufferFloat64Atomics = VK_TRUE;
                fAtom.shaderBufferFloat64AtomicAdd = VK_TRUE;
            }
            fAtom.pNext = feat2.pNext;
            feat2.pNext = &fAtom;
            exts.push_back(VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME);
        }

        // Packed uint8x4 dot product for the descriptor matcher (D21). Core in
        // Vulkan 1.3; we ask for 1.2, so take it as an extension.
        VkPhysicalDeviceShaderIntegerDotProductFeatures fDot{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_INTEGER_DOT_PRODUCT_FEATURES};
        if (opt.needIntDotProduct) {
            fDot.shaderIntegerDotProduct = VK_TRUE;
            fDot.pNext = feat2.pNext;
            feat2.pNext = &fDot;
            exts.push_back(VK_KHR_SHADER_INTEGER_DOT_PRODUCT_EXTENSION_NAME);
        }

        VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
        dci.pNext = &feat2;
        dci.queueCreateInfoCount = 1;
        dci.pQueueCreateInfos = &qci;
        dci.enabledExtensionCount = (uint32_t)exts.size();
        dci.ppEnabledExtensionNames = exts.data();
        VK_CHECK(vkCreateDevice(phys_, &dci, nullptr, &device_));
        vkGetDeviceQueue(device_, queueFamily_, 0, &queue_);

        vkGetPhysicalDeviceMemoryProperties(phys_, &memProps_);

        VkCommandPoolCreateInfo cpi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
        cpi.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
        cpi.queueFamilyIndex = queueFamily_;
        VK_CHECK(vkCreateCommandPool(device_, &cpi, nullptr, &cmdPool_));

        VkFenceCreateInfo fci{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
        VK_CHECK(vkCreateFence(device_, &fci, nullptr, &fence_));

        // Two staging buffers, because the two directions want opposite memory.
        // The first HOST_VISIBLE|HOST_COHERENT type on a discrete GPU is
        // write-combined: excellent to write through, and *uncached to read* --
        // a memcpy out of it runs at a small fraction of RAM speed. Every
        // readback in the pipeline goes through this path (match results, SIFT
        // keypoints and descriptors, BA parameters), so downloads get their own
        // buffer that asks for HOST_CACHED, and fall back to the shared type on
        // a device that has no such heap.
        staging_ = createBufferRaw(kStagingSize,
            VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
        VK_CHECK(vkMapMemory(device_, staging_.mem, 0, kStagingSize, 0, &stagingPtr_));
        stagingDl_ = createBufferRaw(kStagingSize,
            VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
            VK_MEMORY_PROPERTY_HOST_CACHED_BIT);
        VK_CHECK(vkMapMemory(device_, stagingDl_.mem, 0, kStagingSize, 0, &stagingDlPtr_));

        if (opt.profile) {
            profiling_ = true;
            timestampPeriod_ = props.limits.timestampPeriod;
            VkQueryPoolCreateInfo qpi{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
            qpi.queryType = VK_QUERY_TYPE_TIMESTAMP;
            qpi.queryCount = kMaxQueries;
            VK_CHECK(vkCreateQueryPool(device_, &qpi, nullptr, &queryPool_));
        }
    }

    GpuBuffer createBuffer(VkDeviceSize size) {
        GpuBuffer b = createBufferRaw(size,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        totalAllocated_ += size;
        return b;
    }

    // Free one buffer early (a long-lived context whose users come and go --
    // the persistent BA solver -- must not defer everything to ~VkContext). Safe
    // immediately after any submit(): submits are fenced synchronous, so no
    // GPU work can still reference the buffer.
    void destroyBuffer(GpuBuffer& b) {
        if (b.buf == VK_NULL_HANDLE) return;
        vkDestroyBuffer(device_, b.buf, nullptr);
        vkFreeMemory(device_, b.mem, nullptr);
        for (size_t i = 0; i < allocations_.size(); i++)
            if (allocations_[i].first == b.buf) {
                allocations_.erase(allocations_.begin() + i);
                break;
            }
        totalAllocated_ -= b.size;
        b = GpuBuffer{};
    }

    bool initialized() const { return device_ != VK_NULL_HANDLE; }
    bool hasPipeline(const std::string& name) const { return pipelines_.count(name) != 0; }

    void upload(const GpuBuffer& dst, const void* src, VkDeviceSize size, VkDeviceSize dstOff = 0) {
        const uint8_t* p = (const uint8_t*)src;
        for (VkDeviceSize off = 0; off < size; off += kStagingSize) {
            VkDeviceSize chunk = std::min<VkDeviceSize>(kStagingSize, size - off);
            memcpy(stagingPtr_, p + off, chunk);
            VkCommandBuffer cb = begin();
            VkBufferCopy c{0, dstOff + off, chunk};
            vkCmdCopyBuffer(cb, staging_.buf, dst.buf, 1, &c);
            submit(cb);
        }
    }

    void download(const GpuBuffer& src, void* dst, VkDeviceSize size, VkDeviceSize srcOff = 0) {
        uint8_t* p = (uint8_t*)dst;
        for (VkDeviceSize off = 0; off < size; off += kStagingSize) {
            VkDeviceSize chunk = std::min<VkDeviceSize>(kStagingSize, size - off);
            VkCommandBuffer cb = begin();
            VkBufferCopy c{srcOff + off, 0, chunk};
            vkCmdCopyBuffer(cb, src.buf, stagingDl_.buf, 1, &c);
            submit(cb);
            memcpy(p + off, stagingDlPtr_, chunk);
        }
    }

    // Fold a readback into a command buffer the caller is already recording:
    // record the copy with recordDownload(), submit once, then take the bytes
    // with stagingDownloadPtr(). Saves the extra submit-and-fence that a
    // separate download() costs after every dispatch batch. `size` must not
    // exceed stagingCapacity().
    void recordDownload(VkCommandBuffer cb, const GpuBuffer& src, VkDeviceSize size,
                        VkDeviceSize srcOff = 0) {
        VkBufferCopy c{srcOff, 0, size};
        vkCmdCopyBuffer(cb, src.buf, stagingDl_.buf, 1, &c);
    }
    const void* stagingDownloadPtr() const { return stagingDlPtr_; }
    static constexpr VkDeviceSize stagingCapacity() { return kStagingSize; }

    // ---- descriptor set (one set of N storage buffers shared by all pipelines) ----
    // Idempotent: the first call builds layout/pool/set, later calls (same
    // binding count) just rewrite the buffer bindings -- how a persistent
    // context re-targets the pipelines at a new problem's buffers. Safe
    // between submits (fenced synchronous, never while recording).
    void createDescriptors(const std::vector<VkBuffer>& buffers) {
        uint32_t nb = (uint32_t)buffers.size();
        if (setLayout_ != VK_NULL_HANDLE) {
            if (nb != descBindingCount_)
                throw std::runtime_error("createDescriptors: binding count changed");
            writeDescriptors(buffers);
            return;
        }
        descBindingCount_ = nb;
        std::vector<VkDescriptorSetLayoutBinding> binds(nb);
        for (uint32_t i = 0; i < nb; i++)
            binds[i] = {i, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr};
        VkDescriptorSetLayoutCreateInfo li{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        li.bindingCount = nb;
        li.pBindings = binds.data();
        VK_CHECK(vkCreateDescriptorSetLayout(device_, &li, nullptr, &setLayout_));

        VkPushConstantRange pcr{VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(Push)};
        VkPipelineLayoutCreateInfo pli{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pli.setLayoutCount = 1;
        pli.pSetLayouts = &setLayout_;
        pli.pushConstantRangeCount = 1;
        pli.pPushConstantRanges = &pcr;
        VK_CHECK(vkCreatePipelineLayout(device_, &pli, nullptr, &pipeLayout_));

        VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nb};
        VkDescriptorPoolCreateInfo dpi{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        dpi.maxSets = 1;
        dpi.poolSizeCount = 1;
        dpi.pPoolSizes = &ps;
        VK_CHECK(vkCreateDescriptorPool(device_, &dpi, nullptr, &descPool_));

        VkDescriptorSetAllocateInfo dai{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        dai.descriptorPool = descPool_;
        dai.descriptorSetCount = 1;
        dai.pSetLayouts = &setLayout_;
        VK_CHECK(vkAllocateDescriptorSets(device_, &dai, &descSet_));
        writeDescriptors(buffers);
    }

    void writeDescriptors(const std::vector<VkBuffer>& buffers) {
        uint32_t nb = (uint32_t)buffers.size();
        std::vector<VkDescriptorBufferInfo> infos(nb);
        std::vector<VkWriteDescriptorSet> writes(nb);
        for (uint32_t i = 0; i < nb; i++) {
            infos[i] = {buffers[i], 0, VK_WHOLE_SIZE};
            writes[i] = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
            writes[i].dstSet = descSet_;
            writes[i].dstBinding = i;
            writes[i].descriptorCount = 1;
            writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            writes[i].pBufferInfo = &infos[i];
        }
        vkUpdateDescriptorSets(device_, nb, writes.data(), 0, nullptr);
    }

    void loadPipelines(const std::string& spvPath, const std::vector<std::string>& entries) {
        std::ifstream f(spvPath, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("cannot open " + spvPath);
        size_t sz = (size_t)f.tellg();
        f.seekg(0);
        std::vector<char> code(sz);
        f.read(code.data(), sz);
        loadPipelines((const uint32_t*)code.data(), sz, entries);
    }

    // One module holding every entry point (slangc -fvk-use-entrypoint-name),
    // either read from disk above or taken from the embedded blobs.
    void loadPipelines(const uint32_t* code, size_t codeBytes,
                       const std::vector<std::string>& entries) {
        VkShaderModuleCreateInfo smi{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        smi.codeSize = codeBytes;
        smi.pCode = code;
        VK_CHECK(vkCreateShaderModule(device_, &smi, nullptr, &shaderModule_));

        for (const auto& e : entries) {
            VkComputePipelineCreateInfo cpi{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
            cpi.stage = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
            cpi.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
            cpi.stage.module = shaderModule_;
            cpi.stage.pName = e.c_str();
            cpi.layout = pipeLayout_;
            VkPipeline p;
            VK_CHECK(vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &cpi, nullptr, &p));
            pipelines_[e] = p;
        }
    }

    // ---- recording ----
    VkCommandBuffer begin() {
        VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        ai.commandPool = cmdPool_;
        ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        ai.commandBufferCount = 1;
        VkCommandBuffer cb;
        VK_CHECK(vkAllocateCommandBuffers(device_, &ai, &cb));
        VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
        bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        VK_CHECK(vkBeginCommandBuffer(cb, &bi));
        if (descSet_ != VK_NULL_HANDLE)
            vkCmdBindDescriptorSets(cb, VK_PIPELINE_BIND_POINT_COMPUTE, pipeLayout_, 0, 1,
                                    &descSet_, 0, nullptr);
        if (profiling_) {
            vkCmdResetQueryPool(cb, queryPool_, 0, kMaxQueries);
            queryNames_.clear();
        }
        return cb;
    }

    void dispatch(VkCommandBuffer cb, const std::string& name, uint32_t gx, const Push& push,
                  uint32_t gy = 1) {
        auto it = pipelines_.find(name);
        if (it == pipelines_.end()) throw std::runtime_error("unknown pipeline " + name);
        vkCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_COMPUTE, it->second);
        vkCmdPushConstants(cb, pipeLayout_, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(Push), &push);
        bool prof = profiling_ && 2 * queryNames_.size() + 2 <= kMaxQueries;
        if (prof)
            vkCmdWriteTimestamp(cb, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, queryPool_,
                                (uint32_t)(2 * queryNames_.size()));
        vkCmdDispatch(cb, gx, gy, 1);
        if (prof) {
            vkCmdWriteTimestamp(cb, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, queryPool_,
                                (uint32_t)(2 * queryNames_.size() + 1));
            queryNames_.push_back(name);
        }
    }

    void barrier(VkCommandBuffer cb,
                 VkPipelineStageFlags src = VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                 VkPipelineStageFlags dst = VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT) {
        VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        mb.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_TRANSFER_WRITE_BIT;
        mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT |
                           VK_ACCESS_TRANSFER_READ_BIT | VK_ACCESS_TRANSFER_WRITE_BIT;
        vkCmdPipelineBarrier(cb, src | VK_PIPELINE_STAGE_TRANSFER_BIT,
                             dst | VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 1, &mb, 0, nullptr, 0, nullptr);
    }

    void fillZero(VkCommandBuffer cb, const GpuBuffer& b) {
        vkCmdFillBuffer(cb, b.buf, 0, VK_WHOLE_SIZE, 0);
    }

    void copy(VkCommandBuffer cb, const GpuBuffer& src, const GpuBuffer& dst, VkDeviceSize size) {
        VkBufferCopy c{0, 0, size};
        vkCmdCopyBuffer(cb, src.buf, dst.buf, 1, &c);
    }

    void submit(VkCommandBuffer cb) {
        VK_CHECK(vkEndCommandBuffer(cb));
        VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
        si.commandBufferCount = 1;
        si.pCommandBuffers = &cb;
        VK_CHECK(vkQueueSubmit(queue_, 1, &si, fence_));
        VK_CHECK(vkWaitForFences(device_, 1, &fence_, VK_TRUE, ~0ull));
        VK_CHECK(vkResetFences(device_, 1, &fence_));
        vkFreeCommandBuffers(device_, cmdPool_, 1, &cb);
        if (profiling_ && !queryNames_.empty()) {
            std::vector<uint64_t> ts(2 * queryNames_.size());
            VK_CHECK(vkGetQueryPoolResults(device_, queryPool_, 0, (uint32_t)ts.size(),
                                           ts.size() * 8, ts.data(), 8,
                                           VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT));
            for (size_t i = 0; i < queryNames_.size(); i++) {
                auto& s = profStats_[queryNames_[i]];
                s.first += 1;
                s.second += (ts[2 * i + 1] - ts[2 * i]) * (double)timestampPeriod_ * 1e-6;
            }
            queryNames_.clear();
        }
    }

    void printProfile(FILE* f = stderr) const {
        if (!profiling_) return;
        std::vector<std::pair<std::string, std::pair<uint64_t, double>>> v(profStats_.begin(),
                                                                           profStats_.end());
        std::sort(v.begin(), v.end(),
                  [](auto& a, auto& b) { return a.second.second > b.second.second; });
        double total = 0;
        for (auto& e : v) total += e.second.second;
        fprintf(f, "[profile] %-18s %10s %8s %10s\n", "kernel", "total_ms", "count", "avg_us");
        for (auto& e : v)
            fprintf(f, "[profile] %-18s %10.2f %8llu %10.1f\n", e.first.c_str(), e.second.second,
                    (unsigned long long)e.second.first, 1e3 * e.second.second / e.second.first);
        fprintf(f, "[profile] %-18s %10.2f\n", "TOTAL", total);
    }

    VkDevice device() const { return device_; }
    double totalAllocatedMB() const { return totalAllocated_ / (1024.0 * 1024.0); }
    double deviceLocalHeapMB() const {
        VkDeviceSize best = 0;
        for (uint32_t i = 0; i < memProps_.memoryHeapCount; i++)
            if (memProps_.memoryHeaps[i].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT)
                best = std::max(best, memProps_.memoryHeaps[i].size);
        return best / (1024.0 * 1024.0);
    }

private:
    static constexpr VkDeviceSize kStagingSize = 64ull << 20;

    // `preferFlags` are tried on top of `memFlags` and dropped if no memory
    // type offers both.
    GpuBuffer createBufferRaw(VkDeviceSize size, VkBufferUsageFlags usage,
                              VkMemoryPropertyFlags memFlags,
                              VkMemoryPropertyFlags preferFlags = 0) {
        GpuBuffer b;
        b.size = size;
        VkBufferCreateInfo bci{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bci.size = size;
        bci.usage = usage;
        bci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        VK_CHECK(vkCreateBuffer(device_, &bci, nullptr, &b.buf));
        VkMemoryRequirements req;
        vkGetBufferMemoryRequirements(device_, b.buf, &req);
        uint32_t typeIdx = ~0u;
        for (int pass = preferFlags ? 0 : 1; pass < 2 && typeIdx == ~0u; pass++) {
            const VkMemoryPropertyFlags want = pass == 0 ? (memFlags | preferFlags) : memFlags;
            for (uint32_t i = 0; i < memProps_.memoryTypeCount; i++)
                if ((req.memoryTypeBits & (1u << i)) &&
                    (memProps_.memoryTypes[i].propertyFlags & want) == want) {
                    typeIdx = i;
                    break;
                }
        }
        if (typeIdx == ~0u) throw std::runtime_error("no suitable memory type");
        VkMemoryAllocateInfo mai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        mai.allocationSize = req.size;
        mai.memoryTypeIndex = typeIdx;
        VK_CHECK(vkAllocateMemory(device_, &mai, nullptr, &b.mem));
        VK_CHECK(vkBindBufferMemory(device_, b.buf, b.mem, 0));
        allocations_.emplace_back(b.buf, b.mem);
        return b;
    }

    VkInstance instance_ = VK_NULL_HANDLE;
    VkPhysicalDevice phys_ = VK_NULL_HANDLE;
    VkDevice device_ = VK_NULL_HANDLE;
    VkQueue queue_ = VK_NULL_HANDLE;
    uint32_t queueFamily_ = 0;
    VkPhysicalDeviceMemoryProperties memProps_{};
    VkCommandPool cmdPool_ = VK_NULL_HANDLE;
    VkFence fence_ = VK_NULL_HANDLE;
    VkDescriptorSetLayout setLayout_ = VK_NULL_HANDLE;
    VkPipelineLayout pipeLayout_ = VK_NULL_HANDLE;
    VkDescriptorPool descPool_ = VK_NULL_HANDLE;
    VkDescriptorSet descSet_ = VK_NULL_HANDLE;
    uint32_t descBindingCount_ = 0;
    VkShaderModule shaderModule_ = VK_NULL_HANDLE;
    std::map<std::string, VkPipeline> pipelines_;
    GpuBuffer staging_, stagingDl_;
    void* stagingPtr_ = nullptr;
    void* stagingDlPtr_ = nullptr;
    VkDeviceSize totalAllocated_ = 0;
    std::vector<std::pair<VkBuffer, VkDeviceMemory>> allocations_;  // freed in ~VkContext

    static constexpr uint32_t kMaxQueries = 16384;
    bool profiling_ = false;
    float timestampPeriod_ = 1.0f;
    VkQueryPool queryPool_ = VK_NULL_HANDLE;
    std::vector<std::string> queryNames_;
    std::map<std::string, std::pair<uint64_t, double>> profStats_;  // name -> (count, total ms)
};
