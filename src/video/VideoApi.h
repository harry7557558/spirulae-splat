#pragma once
// Entry points for VK_KHR_video_queue / VK_KHR_video_decode_*.
//
// The Vulkan loader does not export extension entry points, so they are
// resolved from the device and passed around as a table. One instance, filled
// by VideoPipeline.cpp on first use -- and re-filled whenever the process
// creates a new device, because the pointers belong to the device they came
// from (see "The entry-point table lives one device, not one process" in
// README.md).

#include <vulkan/vulkan.h>

namespace video {

struct VideoApi {
    PFN_vkGetPhysicalDeviceVideoCapabilitiesKHR       getCapabilities = nullptr;
    PFN_vkGetPhysicalDeviceVideoFormatPropertiesKHR   getFormatProperties = nullptr;
    PFN_vkCreateVideoSessionKHR                       createSession = nullptr;
    PFN_vkDestroyVideoSessionKHR                      destroySession = nullptr;
    PFN_vkGetVideoSessionMemoryRequirementsKHR        getSessionMemoryRequirements = nullptr;
    PFN_vkBindVideoSessionMemoryKHR                   bindSessionMemory = nullptr;
    PFN_vkCreateVideoSessionParametersKHR             createParameters = nullptr;
    PFN_vkUpdateVideoSessionParametersKHR             updateParameters = nullptr;
    PFN_vkDestroyVideoSessionParametersKHR            destroyParameters = nullptr;
    PFN_vkCmdBeginVideoCodingKHR                      cmdBeginCoding = nullptr;
    PFN_vkCmdEndVideoCodingKHR                        cmdEndCoding = nullptr;
    PFN_vkCmdControlVideoCodingKHR                    cmdControlCoding = nullptr;
    PFN_vkCmdDecodeVideoKHR                           cmdDecode = nullptr;

    bool complete() const {
        return getCapabilities && getFormatProperties && createSession && destroySession &&
               getSessionMemoryRequirements && bindSessionMemory && createParameters &&
               destroyParameters && cmdBeginCoding && cmdEndCoding && cmdControlCoding &&
               cmdDecode;
    }
};

// Resolved on first use per device; `complete()` is false when the extensions
// are absent (or when no device exists yet).
const VideoApi& video_api();

}  // namespace video
