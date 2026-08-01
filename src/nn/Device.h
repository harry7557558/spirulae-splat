#pragma once
// Device enumeration and process teardown for the inference layer.
//
// Separate from any model's header so an application can present a device
// picker before a checkpoint is chosen -- enumeration creates no device and
// allocates nothing.

#include <cstdint>
#include <string>
#include <vector>

namespace nn {

struct DeviceInfo {
    int         index = -1;
    std::string name;
    std::string type;              // discrete|integrated|virtual|cpu|other
    uint64_t    vram_bytes = 0;
    bool        usable = false;    // meets the Vulkan 1.2 baseline
    std::string unusable_reason;
};

std::vector<DeviceInfo> list_devices();

// Frees every device allocation and destroys the Vulkan device. Call at exit,
// once every Session, Tracker and VideoReader is gone.
void shutdown();

}  // namespace nn
