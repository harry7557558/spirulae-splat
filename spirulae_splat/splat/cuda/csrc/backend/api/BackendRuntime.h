#pragma once
// Backend-neutral device runtime, covering exactly the CUDA runtime surface
// the engine layer and Tensor.h use today (audited 2026-07: cudaMemcpy[Async],
// cudaMemset[Async], cudaMalloc/Free, cudaMallocHost/FreeHost, cudaEvent*,
// cudaDeviceSynchronize, streams). Engine .cpp files call these instead of
// cuda* directly; the ~54 existing call sites migrate mechanically.
//
// CUDA backend (the default; anything except -DSSPLAT_BACKEND_VULKAN):
// inline wrappers over cudart (zero-cost), included at the bottom of this
// header. Vulkan backend: suballocating pool + staging-buffer copies + fence
// waits (VkSplat model, see backend/README.md).
//
// Device memory is addressed by plain pointers. On Vulkan this presumes
// VK_KHR_buffer_device_address (core 1.2, supported by MoltenVK) so that a
// 64-bit device address behaves pointer-like end to end, including inside
// Slang kernels. The fallback — opaque buffer handle + offset — would ripple
// through Tensor.h and is documented as a rejected-unless-forced alternative.

#include <cstddef>
#include <cstdint>

namespace backend {

enum class MemcpyKind {
    HostToDevice,
    DeviceToHost,
    DeviceToDevice,
    Auto,  // cudaMemcpyDefault; Vulkan resolves by address range
};

// Opaque queue/stream token. Under the CUDA backend this IS cudaStream_t
// (forward-declared, so this header stays CUDA-include-free): declarations
// using backend::Stream and .cu definitions using cudaStream_t are the same
// type — no casts anywhere. Vulkan: command-batch handle (single compute
// queue; batches delimit submission like VkSplat's DEVICE_GUARD scopes).
#ifndef SSPLAT_BACKEND_VULKAN
using Stream = struct CUstream_st*;  // == cudaStream_t
#else
using Stream = void*;
#endif
inline constexpr Stream kDefaultStream = nullptr;

// --- device enumeration / selection ---
// A process runs on ONE device, picked lazily on the first device operation.
// Enumeration is side-effect-free (it does not initialize the backend), so
// apps can list devices and call device_select before starting work.
struct DeviceInfo {
    char name[256];       // human-readable device name ("" if index invalid)
    const char* type;     // "discrete"|"integrated"|"virtual"|"cpu"|"other"
                          // (static storage)
    uint64_t vram_bytes;  // device-local memory
    bool usable;          // meets the backend's feature requirements
};
// Number of devices visible to the backend (0: none / no driver).
int device_count();
// Info for device `index` in [0, device_count()).
DeviceInfo device_info(int index);
// Selects the device used by all subsequent backend work. Call before the
// first device operation. Returns false if `index` is out of range or not
// usable, or if the backend already initialized on a different device.
bool device_select(int index);
// Index of the device in use — or, before the backend initializes, the one
// it would pick (explicit selection, then backend env override, then
// auto-score). -1 if no usable device.
int device_current();

// --- memory ---
void* device_malloc(size_t bytes);
// CONTRACT: device_free must ensure all in-flight device work that may
// reference the allocation completes before the memory is reclaimed
// (cudaFree's implicit sync). DeviceScratch::acquire relies on this when
// growing under a live queue.
void  device_free(void* ptr);
void* host_malloc_pinned(size_t bytes);  // cudaMallocHost / persistently-mapped staging
void  host_free_pinned(void* ptr);

void memcpy_sync(void* dst, const void* src, size_t bytes, MemcpyKind kind);
void memcpy_async(void* dst, const void* src, size_t bytes, MemcpyKind kind,
                  Stream stream = kDefaultStream);
void memset_sync(void* dst, int value, size_t bytes);
void memset_async(void* dst, int value, size_t bytes,
                  Stream stream = kDefaultStream);

// True if `ptr` is a device (or managed) address; false for null and
// pageable-host pointers. Used to auto-route TorchTensorView inputs that may
// be either host or device memory (zero-copy view vs. staging upload).
// Vulkan: address-range check against the pool's buffer-device-address spans.
bool is_device_pointer(const void* ptr);

// Returns a human-readable description of the backend's pending sticky error
// and clears it, or nullptr if none. Like cudaGetLastError, this also surfaces
// failures from earlier async work (kernel launches) — callers that want
// graceful degradation (e.g. the viewer worker) check it after critical
// copies. The returned string has static storage duration.
const char* last_error();

// --- synchronization ---
void device_synchronize();
void stream_synchronize(Stream stream = kDefaultStream);

// --- events (async-readback fencing and perf timing) ---
struct Event;  // opaque per-backend
Event* event_create(bool enable_timing);
void   event_record(Event* event, Stream stream = kDefaultStream);
void   event_synchronize(Event* event);
void   event_destroy(Event* event);
float  event_elapsed_ms(Event* start, Event* end);  // requires enable_timing

}  // namespace backend

#ifndef SSPLAT_BACKEND_VULKAN
#include "../cuda/BackendRuntimeCuda.h"
#endif
