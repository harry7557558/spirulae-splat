#include "nn/vk/Context.h"

#include "nn/core/Error.h"
#include "nn/core/Log.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <mutex>

namespace nn {
namespace vk {

// ================
// Result names
// ================

const char* Context::resultName(VkResult r) {
    switch (r) {
        case VK_SUCCESS:                        return "VK_SUCCESS";
        case VK_NOT_READY:                      return "VK_NOT_READY";
        case VK_TIMEOUT:                        return "VK_TIMEOUT";
        case VK_INCOMPLETE:                     return "VK_INCOMPLETE";
        case VK_ERROR_OUT_OF_HOST_MEMORY:       return "VK_ERROR_OUT_OF_HOST_MEMORY";
        case VK_ERROR_OUT_OF_DEVICE_MEMORY:
            return "VK_ERROR_OUT_OF_DEVICE_MEMORY (the GPU ran out of memory -- "
                   "try a smaller --img-size, a smaller model, or a GPU not "
                   "shared with another job)";
        case VK_ERROR_INITIALIZATION_FAILED:    return "VK_ERROR_INITIALIZATION_FAILED";
        case VK_ERROR_DEVICE_LOST:
            return "VK_ERROR_DEVICE_LOST (a kernel faulted or hung; rerun with "
                   "SSPLAT_NN_DEBUG_SYNC=1 to bisect it to one dispatch)";
        case VK_ERROR_MEMORY_MAP_FAILED:        return "VK_ERROR_MEMORY_MAP_FAILED";
        case VK_ERROR_LAYER_NOT_PRESENT:        return "VK_ERROR_LAYER_NOT_PRESENT";
        case VK_ERROR_EXTENSION_NOT_PRESENT:    return "VK_ERROR_EXTENSION_NOT_PRESENT";
        case VK_ERROR_FEATURE_NOT_PRESENT:      return "VK_ERROR_FEATURE_NOT_PRESENT";
        case VK_ERROR_INCOMPATIBLE_DRIVER:      return "VK_ERROR_INCOMPATIBLE_DRIVER";
        case VK_ERROR_TOO_MANY_OBJECTS:         return "VK_ERROR_TOO_MANY_OBJECTS";
        case VK_ERROR_FORMAT_NOT_SUPPORTED:     return "VK_ERROR_FORMAT_NOT_SUPPORTED";
        case VK_ERROR_FRAGMENTED_POOL:          return "VK_ERROR_FRAGMENTED_POOL";
        case VK_ERROR_UNKNOWN:                  return "VK_ERROR_UNKNOWN";
        // Video decoding rejects an unsupported profile through these, and the
        // name is the whole diagnosis -- "this GPU cannot decode 4:2:2 HEVC".
        case VK_ERROR_VIDEO_PROFILE_OPERATION_NOT_SUPPORTED_KHR:
            return "VK_ERROR_VIDEO_PROFILE_OPERATION_NOT_SUPPORTED_KHR";
        case VK_ERROR_VIDEO_PROFILE_FORMAT_NOT_SUPPORTED_KHR:
            return "VK_ERROR_VIDEO_PROFILE_FORMAT_NOT_SUPPORTED_KHR";
        case VK_ERROR_VIDEO_PROFILE_CODEC_NOT_SUPPORTED_KHR:
            return "VK_ERROR_VIDEO_PROFILE_CODEC_NOT_SUPPORTED_KHR";
        case VK_ERROR_VIDEO_PICTURE_LAYOUT_NOT_SUPPORTED_KHR:
            return "VK_ERROR_VIDEO_PICTURE_LAYOUT_NOT_SUPPORTED_KHR";
        case VK_ERROR_VIDEO_STD_VERSION_NOT_SUPPORTED_KHR:
            return "VK_ERROR_VIDEO_STD_VERSION_NOT_SUPPORTED_KHR";
        default:                                return "VkResult";
    }
}

// ================
// Enumeration
// ================

namespace {

const char* deviceTypeName(VkPhysicalDeviceType t) {
    switch (t) {
        case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:   return "discrete";
        case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return "integrated";
        case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:    return "virtual";
        case VK_PHYSICAL_DEVICE_TYPE_CPU:            return "cpu";
        default:                                     return "other";
    }
}

int typeScore(const char* type) {
    if (std::strcmp(type, "discrete") == 0)   return 4;
    if (std::strcmp(type, "integrated") == 0) return 3;
    if (std::strcmp(type, "virtual") == 0)    return 2;
    if (std::strcmp(type, "cpu") == 0)        return 1;
    return 0;
}

uint64_t deviceLocalHeap(VkPhysicalDevice pd) {
    VkPhysicalDeviceMemoryProperties mp{};
    vkGetPhysicalDeviceMemoryProperties(pd, &mp);
    uint64_t best = 0;
    for (uint32_t i = 0; i < mp.memoryHeapCount; ++i)
        if (mp.memoryHeaps[i].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT)
            best = std::max<uint64_t>(best, mp.memoryHeaps[i].size);
    return best;
}

// Fills `usable` / `unusable_reason` against the Vulkan 1.2 baseline.
void probeBaseline(VkPhysicalDevice pd, DeviceInfo& out) {
    VkPhysicalDeviceProperties props{};
    vkGetPhysicalDeviceProperties(pd, &props);
    out.name = props.deviceName;
    out.type = deviceTypeName(props.deviceType);
    out.vram_bytes = deviceLocalHeap(pd);

    if (props.apiVersion < VK_API_VERSION_1_2) {
        out.unusable_reason = "driver reports Vulkan < 1.2";
        return;
    }

    VkPhysicalDeviceVulkan12Features f12{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    VkPhysicalDeviceFeatures2 f2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    f2.pNext = &f12;
    vkGetPhysicalDeviceFeatures2(pd, &f2);

    if (!f12.bufferDeviceAddress) {
        out.unusable_reason = "no bufferDeviceAddress";
        return;
    }
    if (!f12.timelineSemaphore) {
        out.unusable_reason = "no timelineSemaphore";
        return;
    }

    uint32_t nq = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(pd, &nq, nullptr);
    std::vector<VkQueueFamilyProperties> qf(nq);
    vkGetPhysicalDeviceQueueFamilyProperties(pd, &nq, qf.data());
    bool compute = false;
    for (auto& q : qf)
        if (q.queueFlags & VK_QUEUE_COMPUTE_BIT) compute = true;
    if (!compute) {
        out.unusable_reason = "no compute queue";
        return;
    }
    out.usable = true;
}

VkInstance createInstance(bool validation, VkDebugUtilsMessengerEXT* messenger);

VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
    VkDebugUtilsMessageSeverityFlagBitsEXT severity,
    VkDebugUtilsMessageTypeFlagsEXT, const VkDebugUtilsMessengerCallbackDataEXT* data,
    void*) {
    if (severity >= VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT)
        NN_LOG_WARN("[vk-validation] %s\n", data->pMessage);
    return VK_FALSE;
}

VkInstance createInstance(bool validation, VkDebugUtilsMessengerEXT* messenger) {
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "ssam";
    app.apiVersion = VK_API_VERSION_1_3;  // request 1.3, we only *require* 1.2

    std::vector<const char*> layers;
    std::vector<const char*> exts;
    if (validation) {
        uint32_t n = 0;
        vkEnumerateInstanceLayerProperties(&n, nullptr);
        std::vector<VkLayerProperties> props(n);
        vkEnumerateInstanceLayerProperties(&n, props.data());
        for (auto& p : props)
            if (std::strcmp(p.layerName, "VK_LAYER_KHRONOS_validation") == 0)
                layers.push_back("VK_LAYER_KHRONOS_validation");
        if (layers.empty())
            NN_LOG_WARN("[vk] validation requested but the layer is not installed\n");
        else
            exts.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    }

    VkInstanceCreateInfo ici{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ici.pApplicationInfo = &app;
    ici.enabledLayerCount = (uint32_t)layers.size();
    ici.ppEnabledLayerNames = layers.data();
    ici.enabledExtensionCount = (uint32_t)exts.size();
    ici.ppEnabledExtensionNames = exts.data();

    VkInstance inst = VK_NULL_HANDLE;
    VkResult r = vkCreateInstance(&ici, nullptr, &inst);
    if (r == VK_ERROR_INCOMPATIBLE_DRIVER) {
        // Retry at 1.2 for loaders that reject a 1.3 request outright.
        app.apiVersion = VK_API_VERSION_1_2;
        r = vkCreateInstance(&ici, nullptr, &inst);
    }
    NN_CHECK(r == VK_SUCCESS,
               "vkCreateInstance failed: %s (%d). Is a Vulkan driver installed?",
               Context::resultName(r), (int)r);

    if (messenger && !layers.empty()) {
        auto fn = (PFN_vkCreateDebugUtilsMessengerEXT)vkGetInstanceProcAddr(
            inst, "vkCreateDebugUtilsMessengerEXT");
        if (fn) {
            VkDebugUtilsMessengerCreateInfoEXT dci{
                VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT};
            dci.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT |
                                  VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
            dci.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT |
                              VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                              VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
            dci.pfnUserCallback = debugCallback;
            fn(inst, &dci, nullptr, messenger);
        }
    }
    return inst;
}

bool envFlag(const char* name) {
    const char* v = std::getenv(name);
    return v && v[0] && std::strcmp(v, "0") != 0;
}

std::string lower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

}  // namespace

std::vector<DeviceInfo> enumerate_devices() {
    std::vector<DeviceInfo> out;
    VkInstance inst = createInstance(false, nullptr);
    uint32_t n = 0;
    vkEnumeratePhysicalDevices(inst, &n, nullptr);
    std::vector<VkPhysicalDevice> devs(n);
    if (n) vkEnumeratePhysicalDevices(inst, &n, devs.data());
    for (uint32_t i = 0; i < n; ++i) {
        DeviceInfo di;
        di.index = (int)i;
        probeBaseline(devs[i], di);
        out.push_back(std::move(di));
    }
    vkDestroyInstance(inst, nullptr);
    return out;
}

// ================
// Lifetime
// ================

namespace {
Context*   g_ctx = nullptr;
std::mutex g_ctx_mu;
}  // namespace

bool Context::initialized() { return g_ctx != nullptr; }

Context& Context::get(const ContextOptions& opts) {
    std::lock_guard<std::mutex> lock(g_ctx_mu);
    if (!g_ctx) {
        g_ctx = new Context();
        g_ctx->init(opts);
    }
    return *g_ctx;
}

void Context::shutdown() {
    std::lock_guard<std::mutex> lock(g_ctx_mu);
    delete g_ctx;
    g_ctx = nullptr;
}

Context::~Context() {
    if (device_) {
        vkDeviceWaitIdle(device_);
        vkDestroyDevice(device_, nullptr);
    }
    if (debug_messenger_) {
        auto fn = (PFN_vkDestroyDebugUtilsMessengerEXT)vkGetInstanceProcAddr(
            instance_, "vkDestroyDebugUtilsMessengerEXT");
        if (fn) fn(instance_, debug_messenger_, nullptr);
    }
    if (instance_) vkDestroyInstance(instance_, nullptr);
}

void Context::init(const ContextOptions& opts) {
    const bool validation = opts.validation || envFlag("SSPLAT_VK_VALIDATION");
    profiling_ = opts.profile || envFlag("SSPLAT_PROFILE");

    instance_ = createInstance(validation, &debug_messenger_);
    pickPhysicalDevice(opts);
    createDevice(opts);

    NN_LOG_INFO("[vk] device: %s (%s, %.1f GiB)\n", info_.name.c_str(), info_.type,
                  info_.vram_bytes / (1024.0 * 1024.0 * 1024.0));
}

void Context::pickPhysicalDevice(const ContextOptions& opts) {
    uint32_t n = 0;
    vkEnumeratePhysicalDevices(instance_, &n, nullptr);
    NN_CHECK(n > 0, "no Vulkan physical devices found");
    std::vector<VkPhysicalDevice> devs(n);
    vkEnumeratePhysicalDevices(instance_, &n, devs.data());

    std::vector<DeviceInfo> infos(n);
    for (uint32_t i = 0; i < n; ++i) {
        infos[i].index = (int)i;
        probeBaseline(devs[i], infos[i]);
    }

    // Explicit selection: option, then $SSPLAT_VK_DEVICE (index or name substring).
    int chosen = -1;
    std::string match = opts.device_match;
    int want_index = opts.device_index;
    if (const char* env = std::getenv("SSPLAT_VK_DEVICE")) {
        char* end = nullptr;
        long v = std::strtol(env, &end, 10);
        if (end && *end == '\0') want_index = (int)v;
        else match = env;
    }
    if (want_index >= 0) {
        NN_CHECK(want_index < (int)n, "SSPLAT_VK_DEVICE index %d out of range (%u devices)",
                   want_index, n);
        chosen = want_index;
    } else if (!match.empty()) {
        std::string needle = lower(match);
        for (uint32_t i = 0; i < n; ++i)
            if (lower(infos[i].name).find(needle) != std::string::npos) {
                chosen = (int)i;
                break;
            }
        NN_CHECK(chosen >= 0, "no Vulkan device name contains '%s'", match.c_str());
    } else {
        int best_score = -1;
        uint64_t best_vram = 0;
        for (uint32_t i = 0; i < n; ++i) {
            if (!infos[i].usable) continue;
            int s = typeScore(infos[i].type);
            if (s > best_score || (s == best_score && infos[i].vram_bytes > best_vram)) {
                best_score = s;
                best_vram = infos[i].vram_bytes;
                chosen = (int)i;
            }
        }
        NN_CHECK(chosen >= 0,
                   "no Vulkan device meets the baseline (1.2 + bufferDeviceAddress "
                   "+ timelineSemaphore); first device said: %s",
                   infos[0].unusable_reason.c_str());
    }

    NN_CHECK(infos[chosen].usable, "device %d (%s) is unusable: %s", chosen,
               infos[chosen].name.c_str(), infos[chosen].unusable_reason.c_str());

    physical_ = devs[chosen];
    info_ = infos[chosen];
    VkPhysicalDeviceProperties props{};
    vkGetPhysicalDeviceProperties(physical_, &props);
    limits_ = props.limits;
    vkGetPhysicalDeviceMemoryProperties(physical_, &mem_props_);
}

void Context::createDevice(const ContextOptions& opts) {
    uint32_t nq = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(physical_, &nq, nullptr);
    std::vector<VkQueueFamilyProperties2> qf2(nq);
    std::vector<VkQueueFamilyVideoPropertiesKHR> qvid(nq);
    for (uint32_t i = 0; i < nq; ++i) {
        qf2[i] = {VK_STRUCTURE_TYPE_QUEUE_FAMILY_PROPERTIES_2};
        qvid[i] = {VK_STRUCTURE_TYPE_QUEUE_FAMILY_VIDEO_PROPERTIES_KHR};
        qf2[i].pNext = &qvid[i];
    }
    vkGetPhysicalDeviceQueueFamilyProperties2(physical_, &nq, qf2.data());

    // Prefer a compute family without graphics (async compute) when one exists.
    for (uint32_t i = 0; i < nq; ++i) {
        const auto& p = qf2[i].queueFamilyProperties;
        if ((p.queueFlags & VK_QUEUE_COMPUTE_BIT) &&
            !(p.queueFlags & VK_QUEUE_GRAPHICS_BIT)) {
            queue_family_ = i;
            break;
        }
    }
    if (queue_family_ == UINT32_MAX)
        for (uint32_t i = 0; i < nq; ++i)
            if (qf2[i].queueFamilyProperties.queueFlags & VK_QUEUE_COMPUTE_BIT) {
                queue_family_ = i;
                break;
            }
    NN_CHECK(queue_family_ != UINT32_MAX, "no compute queue family");

    // ---- optional extensions -------------------------------------------
    uint32_t ne = 0;
    vkEnumerateDeviceExtensionProperties(physical_, nullptr, &ne, nullptr);
    std::vector<VkExtensionProperties> avail(ne);
    vkEnumerateDeviceExtensionProperties(physical_, nullptr, &ne, avail.data());
    auto has = [&](const char* name) {
        for (auto& e : avail)
            if (std::strcmp(e.extensionName, name) == 0) return true;
        return false;
    };

    std::vector<const char*> exts;
    subgroup_size_control_ = has(VK_EXT_SUBGROUP_SIZE_CONTROL_EXTENSION_NAME);
    if (subgroup_size_control_)
        exts.push_back(VK_EXT_SUBGROUP_SIZE_CONTROL_EXTENSION_NAME);

    // ---- video decode ---------------------------------------------------
    bool want_video = opts.want_video;
#ifdef SSPLAT_HAVE_VIDEO
    want_video = true;  // always probe; the decoder is created lazily
#endif
    if (want_video) {
        const bool base = has(VK_KHR_VIDEO_QUEUE_EXTENSION_NAME) &&
                          has(VK_KHR_VIDEO_DECODE_QUEUE_EXTENSION_NAME);
        if (!base) {
            video_reason_ = "driver does not expose VK_KHR_video_queue / "
                            "VK_KHR_video_decode_queue";
        } else {
            for (uint32_t i = 0; i < nq; ++i) {
                if (!(qf2[i].queueFamilyProperties.queueFlags &
                      VK_QUEUE_VIDEO_DECODE_BIT_KHR))
                    continue;
                video_queue_family_ = i;
                video_codec_ops_ = qvid[i].videoCodecOperations;
                break;
            }
            if (video_queue_family_ == UINT32_MAX) {
                video_reason_ = "no queue family advertises VIDEO_DECODE";
            } else {
                exts.push_back(VK_KHR_VIDEO_QUEUE_EXTENSION_NAME);
                exts.push_back(VK_KHR_VIDEO_DECODE_QUEUE_EXTENSION_NAME);
                if (has(VK_KHR_VIDEO_DECODE_H264_EXTENSION_NAME))
                    exts.push_back(VK_KHR_VIDEO_DECODE_H264_EXTENSION_NAME);
                else
                    video_codec_ops_ &= ~(uint32_t)VK_VIDEO_CODEC_OPERATION_DECODE_H264_BIT_KHR;
                if (has(VK_KHR_VIDEO_DECODE_H265_EXTENSION_NAME))
                    exts.push_back(VK_KHR_VIDEO_DECODE_H265_EXTENSION_NAME);
                else
                    video_codec_ops_ &= ~(uint32_t)VK_VIDEO_CODEC_OPERATION_DECODE_H265_BIT_KHR;
#ifdef VK_KHR_video_decode_av1
                if (has(VK_KHR_VIDEO_DECODE_AV1_EXTENSION_NAME))
                    exts.push_back(VK_KHR_VIDEO_DECODE_AV1_EXTENSION_NAME);
                else
                    video_codec_ops_ &= ~(uint32_t)VK_VIDEO_CODEC_OPERATION_DECODE_AV1_BIT_KHR;
#endif
                if (video_codec_ops_ == 0) {
                    video_reason_ = "video-decode queue exists but no supported codec "
                                    "(H.264/H.265/AV1) is exposed";
                    video_queue_family_ = UINT32_MAX;
                    // Drop the extensions we just added: nothing will use them.
                    while (!exts.empty() &&
                           std::strncmp(exts.back(), "VK_KHR_video", 12) == 0)
                        exts.pop_back();
                } else {
                    video_reason_.clear();
                }
            }
        }
    }

    // ---- queues ---------------------------------------------------------
    const float prio = 1.0f;
    std::vector<VkDeviceQueueCreateInfo> qcis;
    {
        VkDeviceQueueCreateInfo q{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        q.queueFamilyIndex = queue_family_;
        q.queueCount = 1;
        q.pQueuePriorities = &prio;
        qcis.push_back(q);
        if (video_queue_family_ != UINT32_MAX && video_queue_family_ != queue_family_) {
            q.queueFamilyIndex = video_queue_family_;
            qcis.push_back(q);
        }
    }

    VkPhysicalDeviceVulkan12Features f12{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    f12.bufferDeviceAddress = VK_TRUE;
    f12.timelineSemaphore = VK_TRUE;

    VkPhysicalDeviceVulkan11Features f11{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES};
    f11.pNext = &f12;
    // Decoded pictures are multi-planar YCbCr images. We never sample them --
    // the planes are copied to buffers and unpacked in a shader -- but creating
    // an image in one of those formats is gated on this feature.
    if (video_queue_family_ != UINT32_MAX) f11.samplerYcbcrConversion = VK_TRUE;

    VkPhysicalDeviceSubgroupSizeControlFeaturesEXT fsg{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_SIZE_CONTROL_FEATURES_EXT};
    if (subgroup_size_control_) {
        fsg.subgroupSizeControl = VK_TRUE;
        fsg.computeFullSubgroups = VK_TRUE;
        fsg.pNext = &f11;
    }

    VkPhysicalDeviceFeatures2 f2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    f2.pNext = subgroup_size_control_ ? (void*)&fsg : (void*)&f11;

    VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.pNext = &f2;
    dci.queueCreateInfoCount = (uint32_t)qcis.size();
    dci.pQueueCreateInfos = qcis.data();
    dci.enabledExtensionCount = (uint32_t)exts.size();
    dci.ppEnabledExtensionNames = exts.data();
    NN_VK_CHECK(vkCreateDevice(physical_, &dci, nullptr, &device_));

    vkGetDeviceQueue(device_, queue_family_, 0, &queue_);
    if (video_queue_family_ != UINT32_MAX)
        vkGetDeviceQueue(device_, video_queue_family_,
                         video_queue_family_ == queue_family_ ? 0 : 0, &video_queue_);

    if (subgroup_size_control_) {
        VkPhysicalDeviceSubgroupSizeControlPropertiesEXT sp{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_SIZE_CONTROL_PROPERTIES_EXT};
        VkPhysicalDeviceProperties2 p2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
        p2.pNext = &sp;
        vkGetPhysicalDeviceProperties2(physical_, &p2);
        // Pin the widest size the device offers that our 64-wide-X kernels can
        // still fill; Intel/ANV's varying width is the failure this avoids.
        preferred_subgroup_ = std::min<uint32_t>(sp.maxSubgroupSize, 64u);
        preferred_subgroup_ = std::max<uint32_t>(preferred_subgroup_, sp.minSubgroupSize);
        NN_LOG_DEBUG("[vk] subgroup size pinned to %u (range %u..%u)\n",
                       preferred_subgroup_, sp.minSubgroupSize, sp.maxSubgroupSize);
    }
}

bool Context::hasVideoCodec(VkVideoCodecOperationFlagBitsKHR op) const {
    return (video_codec_ops_ & (uint32_t)op) != 0;
}

PFN_vkVoidFunction Context::deviceProc(const char* name) const {
    return vkGetDeviceProcAddr(device_, name);
}

}  // namespace vk
}  // namespace nn
