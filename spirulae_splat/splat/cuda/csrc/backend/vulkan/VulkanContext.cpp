#include "VulkanContext.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace backend {
namespace vk {

// --- sticky error (definition of backend::last_error lives in
// VulkanRuntime.cpp; it reads these) ---
std::mutex g_error_mutex;
std::string g_error;

void set_error(const char* what, VkResult result) {
    std::lock_guard<std::mutex> lock(g_error_mutex);
    if (!g_error.empty()) return;  // sticky: keep the FIRST error
    g_error = what;
    if (result != VK_SUCCESS) {
        g_error += " (VkResult ";
        g_error += std::to_string((int)result);
        g_error += ")";
    }
}

namespace {

const char* device_type_name(VkPhysicalDeviceType t) {
    switch (t) {
        case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU: return "discrete";
        case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return "integrated";
        case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU: return "virtual";
        case VK_PHYSICAL_DEVICE_TYPE_CPU: return "cpu";
        default: return "other";
    }
}

int device_type_score(VkPhysicalDeviceType t) {
    switch (t) {
        case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU: return 3;
        case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return 2;
        case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU: return 1;
        case VK_PHYSICAL_DEVICE_TYPE_CPU: return 0;
        default: return 0;
    }
}

VkDeviceSize device_local_vram(VkPhysicalDevice pd) {
    VkPhysicalDeviceMemoryProperties mp;
    vkGetPhysicalDeviceMemoryProperties(pd, &mp);
    VkDeviceSize total = 0;
    for (uint32_t i = 0; i < mp.memoryHeapCount; i++)
        if (mp.memoryHeaps[i].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT)
            total += mp.memoryHeaps[i].size;
    return total;
}

// Compute queue family with timestamp support preferred.
uint32_t pick_compute_queue_family(VkPhysicalDevice pd) {
    uint32_t n = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(pd, &n, nullptr);
    std::vector<VkQueueFamilyProperties> qf(n);
    vkGetPhysicalDeviceQueueFamilyProperties(pd, &n, qf.data());
    uint32_t best = UINT32_MAX;
    for (uint32_t i = 0; i < n; i++) {
        if (!(qf[i].queueFlags & VK_QUEUE_COMPUTE_BIT)) continue;
        if (best == UINT32_MAX) best = i;
        if (qf[i].timestampValidBits > 0) return i;
    }
    return best;
}

struct DeviceProbe {
    bool required_ok = false;
    bool atomic_float = false;
    uint32_t queue_family = UINT32_MAX;
    VkPhysicalDeviceProperties props{};
    VkPhysicalDeviceSubgroupProperties subgroup{};
};

DeviceProbe probe_device(VkPhysicalDevice pd) {
    DeviceProbe out;

    VkPhysicalDeviceVulkan12Features f12{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    VkPhysicalDeviceShaderAtomicFloatFeaturesEXT fatomic{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
    f12.pNext = &fatomic;
    VkPhysicalDeviceFeatures2 f2{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    f2.pNext = &f12;
    vkGetPhysicalDeviceFeatures2(pd, &f2);

    out.subgroup.sType =
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES;
    VkPhysicalDeviceProperties2 p2{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    p2.pNext = &out.subgroup;
    vkGetPhysicalDeviceProperties2(pd, &p2);
    out.props = p2.properties;

    out.queue_family = pick_compute_queue_family(pd);
    out.required_ok =
        out.props.apiVersion >= VK_API_VERSION_1_2 &&
        f2.features.shaderInt64 &&
        f12.bufferDeviceAddress &&
        f12.timelineSemaphore &&
        out.queue_family != UINT32_MAX;
    out.atomic_float = fatomic.shaderBufferFloat32AtomicAdd;
    return out;
}

bool has_extension(VkPhysicalDevice pd, const char* name) {
    uint32_t n = 0;
    vkEnumerateDeviceExtensionProperties(pd, nullptr, &n, nullptr);
    std::vector<VkExtensionProperties> exts(n);
    vkEnumerateDeviceExtensionProperties(pd, nullptr, &n, exts.data());
    for (const auto& e : exts)
        if (std::strcmp(e.extensionName, name) == 0) return true;
    return false;
}

}  // namespace

Context& Context::get() {
    static Context ctx;
    return ctx;
}

Context::Context() { init(); }

void Context::init() {
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "spirulae-splat";
    app.apiVersion = VK_API_VERSION_1_2;

    VkInstanceCreateInfo ici{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ici.pApplicationInfo = &app;
    const char* validation = "VK_LAYER_KHRONOS_validation";
    if (const char* env = std::getenv("SSPLAT_VK_VALIDATION");
        env && env[0] == '1') {
        ici.enabledLayerCount = 1;
        ici.ppEnabledLayerNames = &validation;
    }
    VkResult r = vkCreateInstance(&ici, nullptr, &_instance);
    if (r != VK_SUCCESS) {
        set_error("vkCreateInstance failed", r);
        return;
    }

    uint32_t n = 0;
    vkEnumeratePhysicalDevices(_instance, &n, nullptr);
    std::vector<VkPhysicalDevice> devices(n);
    vkEnumeratePhysicalDevices(_instance, &n, devices.data());
    if (n == 0) {
        set_error("no Vulkan devices found", VK_SUCCESS);
        return;
    }

    const char* want = std::getenv("SSPLAT_VK_DEVICE");
    int best = -1;
    long best_score = -1;
    std::vector<DeviceProbe> probes(n);
    for (uint32_t i = 0; i < n; i++) {
        probes[i] = probe_device(devices[i]);
        const DeviceProbe& p = probes[i];
        bool viable = p.required_ok;
        if (want && want[0]) {
            char* end = nullptr;
            long idx = std::strtol(want, &end, 10);
            bool matches = (end && *end == '\0')
                ? (idx == (long)i)
                : (std::strstr(p.props.deviceName, want) != nullptr);
            if (!matches) continue;
            if (!viable) {
                std::fprintf(stderr,
                    "[ssplat-vk] requested device '%s' lacks required "
                    "features (Vulkan 1.2 + shaderInt64 + "
                    "bufferDeviceAddress + timelineSemaphore)\n",
                    p.props.deviceName);
                continue;
            }
            best = (int)i;
            break;
        }
        if (!viable) continue;
        long score = device_type_score(p.props.deviceType) * 1000000 +
                     (long)(device_local_vram(devices[i]) >> 30);
        if (score > best_score) {
            best_score = score;
            best = (int)i;
        }
    }
    if (best < 0) {
        set_error("no viable Vulkan device (need Vulkan 1.2 + shaderInt64 + "
                  "bufferDeviceAddress + timelineSemaphore)", VK_SUCCESS);
        return;
    }

    _physical = devices[best];
    const DeviceProbe& probe = probes[best];
    _queue_family = probe.queue_family;
    _device_name = probe.props.deviceName;
    vkGetPhysicalDeviceMemoryProperties(_physical, &_mem_props);

    // Queue family timestamp support (for event timing).
    {
        uint32_t qn = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(_physical, &qn, nullptr);
        std::vector<VkQueueFamilyProperties> qf(qn);
        vkGetPhysicalDeviceQueueFamilyProperties(_physical, &qn, qf.data());
        _caps.timestamps = qf[_queue_family].timestampValidBits > 0;
    }
    _caps.timestamp_period_ns = probe.props.limits.timestampPeriod;
    _caps.subgroup_size = probe.subgroup.subgroupSize;
    _caps.max_push_constants = probe.props.limits.maxPushConstantsSize;
    _caps.max_workgroup_invocations =
        probe.props.limits.maxComputeWorkGroupInvocations;
    _caps.max_shared_memory = probe.props.limits.maxComputeSharedMemorySize;
    _caps.non_coherent_atom_size = probe.props.limits.nonCoherentAtomSize;
    _caps.float32_atomic_add = probe.atomic_float;

    float prio = 1.0f;
    VkDeviceQueueCreateInfo qci{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qci.queueFamilyIndex = _queue_family;
    qci.queueCount = 1;
    qci.pQueuePriorities = &prio;

    VkPhysicalDeviceVulkan12Features f12{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    f12.bufferDeviceAddress = VK_TRUE;
    f12.timelineSemaphore = VK_TRUE;
    // Optional but broadly available and useful for kernel ABIs.
    f12.shaderInt8 = VK_FALSE;

    VkPhysicalDeviceShaderAtomicFloatFeaturesEXT fatomic{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
    fatomic.shaderBufferFloat32AtomicAdd = VK_TRUE;

    std::vector<const char*> extensions;
    if (_caps.float32_atomic_add &&
        has_extension(_physical, VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME)) {
        extensions.push_back(VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME);
        f12.pNext = &fatomic;
    } else {
        _caps.float32_atomic_add = false;
    }

    // Pin the compute subgroup size when the device lets us (see
    // Capabilities::required_subgroup_size).
    VkPhysicalDeviceSubgroupSizeControlFeaturesEXT fsgc{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_SIZE_CONTROL_FEATURES_EXT};
    if (has_extension(_physical,
                      VK_EXT_SUBGROUP_SIZE_CONTROL_EXTENSION_NAME)) {
        VkPhysicalDeviceSubgroupSizeControlFeaturesEXT probe_sgc{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_SIZE_CONTROL_FEATURES_EXT};
        VkPhysicalDeviceFeatures2 probe_f2{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
        probe_f2.pNext = &probe_sgc;
        vkGetPhysicalDeviceFeatures2(_physical, &probe_f2);

        VkPhysicalDeviceSubgroupSizeControlPropertiesEXT sgc_props{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_SIZE_CONTROL_PROPERTIES_EXT};
        VkPhysicalDeviceProperties2 sgc_p2{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
        sgc_p2.pNext = &sgc_props;
        vkGetPhysicalDeviceProperties2(_physical, &sgc_p2);

        if (probe_sgc.subgroupSizeControl && probe_sgc.computeFullSubgroups &&
            (sgc_props.requiredSubgroupSizeStages &
             VK_SHADER_STAGE_COMPUTE_BIT)) {
            uint32_t want = _caps.subgroup_size;
            if (want < sgc_props.minSubgroupSize)
                want = sgc_props.minSubgroupSize;
            if (want > sgc_props.maxSubgroupSize)
                want = sgc_props.maxSubgroupSize;
            _caps.required_subgroup_size = want;
            extensions.push_back(VK_EXT_SUBGROUP_SIZE_CONTROL_EXTENSION_NAME);
            fsgc.subgroupSizeControl = VK_TRUE;
            fsgc.computeFullSubgroups = VK_TRUE;
            fsgc.pNext = f12.pNext;
            f12.pNext = &fsgc;
        }
    }

    VkPhysicalDeviceFeatures features{};
    features.shaderInt64 = VK_TRUE;

    VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.pNext = &f12;
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos = &qci;
    dci.pEnabledFeatures = &features;
    dci.enabledExtensionCount = (uint32_t)extensions.size();
    dci.ppEnabledExtensionNames = extensions.data();

    r = vkCreateDevice(_physical, &dci, nullptr, &_device);
    if (r != VK_SUCCESS) {
        set_error("vkCreateDevice failed", r);
        _device = VK_NULL_HANDLE;
        return;
    }
    vkGetDeviceQueue(_device, _queue_family, 0, &_queue);

    VkSemaphoreTypeCreateInfo tci{VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO};
    tci.semaphoreType = VK_SEMAPHORE_TYPE_TIMELINE;
    tci.initialValue = 0;
    VkSemaphoreCreateInfo sci{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    sci.pNext = &tci;
    r = vkCreateSemaphore(_device, &sci, nullptr, &_timeline);
    if (r != VK_SUCCESS) {
        set_error("vkCreateSemaphore (timeline) failed", r);
        vkDestroyDevice(_device, nullptr);
        _device = VK_NULL_HANDLE;
        return;
    }

    if (std::getenv("SSPLAT_VK_VERBOSE")) {
        std::fprintf(stderr,
            "[ssplat-vk] using %s (%s), subgroup %u, push %uB, "
            "float-atomic-add %s, timestamps %s\n",
            _device_name.c_str(), device_type_name(probe.props.deviceType),
            _caps.subgroup_size, _caps.max_push_constants,
            _caps.float32_atomic_add ? "native" : "EMULATED",
            _caps.timestamps ? "yes" : "no");
    }
}

Context::~Context() {
    if (_device != VK_NULL_HANDLE) {
        vkDeviceWaitIdle(_device);
        if (_shutdown_hook) _shutdown_hook();
        vkDestroySemaphore(_device, _timeline, nullptr);
        vkDestroyDevice(_device, nullptr);
    }
    if (_instance != VK_NULL_HANDLE) vkDestroyInstance(_instance, nullptr);
}

uint64_t Context::submit(VkCommandBuffer cb) {
    std::lock_guard<std::mutex> lock(_submit_mutex);
    uint64_t value = _last_value + 1;

    VkTimelineSemaphoreSubmitInfo tsi{
        VK_STRUCTURE_TYPE_TIMELINE_SEMAPHORE_SUBMIT_INFO};
    tsi.signalSemaphoreValueCount = 1;
    tsi.pSignalSemaphoreValues = &value;

    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.pNext = &tsi;
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cb;
    si.signalSemaphoreCount = 1;
    si.pSignalSemaphores = &_timeline;

    VkResult r = vkQueueSubmit(_queue, 1, &si, VK_NULL_HANDLE);
    if (r != VK_SUCCESS) {
        set_error("vkQueueSubmit failed", r);
        return 0;
    }
    _last_value = value;
    return value;
}

bool Context::wait(uint64_t value) {
    if (value == 0) return true;
    VkSemaphoreWaitInfo wi{VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
    wi.semaphoreCount = 1;
    wi.pSemaphores = &_timeline;
    wi.pValues = &value;
    VkResult r = vkWaitSemaphores(_device, &wi, UINT64_MAX);
    if (r != VK_SUCCESS) {
        set_error("vkWaitSemaphores failed", r);
        return false;
    }
    return true;
}

uint32_t Context::find_memory_type(uint32_t type_bits,
                                   VkMemoryPropertyFlags required) const {
    for (uint32_t i = 0; i < _mem_props.memoryTypeCount; i++) {
        if (!(type_bits & (1u << i))) continue;
        if ((_mem_props.memoryTypes[i].propertyFlags & required) == required)
            return i;
    }
    return UINT32_MAX;
}

}  // namespace vk
}  // namespace backend
