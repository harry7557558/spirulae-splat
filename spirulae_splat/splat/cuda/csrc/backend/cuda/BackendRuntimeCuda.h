#pragma once
// Inline cudart implementation of the backend runtime (the default backend).
// Included by BackendRuntime.h — do not include directly. Calls are
// unchecked, matching the raw cuda* call sites this layer replaced.

#include <cuda_runtime.h>

#include <cstdio>

namespace backend {

// --- device enumeration / selection ---
inline int device_count() {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) {
        cudaGetLastError();
        return 0;
    }
    return n;
}
inline DeviceInfo device_info(int index) {
    DeviceInfo info{};
    info.type = "other";
    cudaDeviceProp prop{};
    if (index < 0 || cudaGetDeviceProperties(&prop, index) != cudaSuccess) {
        cudaGetLastError();
        return info;
    }
    std::snprintf(info.name, sizeof(info.name), "%s", prop.name);
    info.type = prop.integrated ? "integrated" : "discrete";
    info.vram_bytes = (uint64_t)prop.totalGlobalMem;
    info.usable = true;
    return info;
}
inline bool device_select(int index) {
    if (index < 0 || index >= device_count()) return false;
    if (cudaSetDevice(index) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    return true;
}
inline int device_current() {
    int d = -1;
    if (cudaGetDevice(&d) != cudaSuccess) {
        cudaGetLastError();
        return -1;
    }
    return d;
}

inline cudaMemcpyKind _to_cuda(MemcpyKind kind) {
    switch (kind) {
        case MemcpyKind::HostToDevice:   return cudaMemcpyHostToDevice;
        case MemcpyKind::DeviceToHost:   return cudaMemcpyDeviceToHost;
        case MemcpyKind::DeviceToDevice: return cudaMemcpyDeviceToDevice;
        default:                         return cudaMemcpyDefault;
    }
}

inline cudaStream_t _to_cuda(Stream stream) {
    return stream;  // backend::Stream IS cudaStream_t under this backend
}

inline const char* last_error() {
    cudaError_t err = cudaGetLastError();
    return err == cudaSuccess ? nullptr : cudaGetErrorString(err);
}

inline bool is_device_pointer(const void* ptr) {
    if (ptr == nullptr) return false;
    cudaPointerAttributes attr{};
    cudaError_t err = cudaPointerGetAttributes(&attr, ptr);
    if (err != cudaSuccess) {
        // Pageable host pointer is unregistered -> reset error, treat as host.
        cudaGetLastError();
        return false;
    }
    return attr.type == cudaMemoryTypeDevice || attr.type == cudaMemoryTypeManaged;
}

// --- memory ---
inline void* device_malloc(size_t bytes) {
    void* ptr = nullptr;
    cudaMalloc(&ptr, bytes);
    return ptr;
}
inline void device_free(void* ptr) { cudaFree(ptr); }
inline void* host_malloc_pinned(size_t bytes) {
    void* ptr = nullptr;
    cudaMallocHost(&ptr, bytes);
    return ptr;
}
inline void host_free_pinned(void* ptr) { cudaFreeHost(ptr); }

inline void memcpy_sync(void* dst, const void* src, size_t bytes, MemcpyKind kind) {
    cudaMemcpy(dst, src, bytes, _to_cuda(kind));
}
inline void memcpy_async(void* dst, const void* src, size_t bytes, MemcpyKind kind,
                         Stream stream) {
    cudaMemcpyAsync(dst, src, bytes, _to_cuda(kind), _to_cuda(stream));
}
inline void memset_sync(void* dst, int value, size_t bytes) {
    cudaMemset(dst, value, bytes);
}
inline void memset_async(void* dst, int value, size_t bytes, Stream stream) {
    cudaMemsetAsync(dst, value, bytes, _to_cuda(stream));
}

// --- synchronization ---
inline void device_synchronize() { cudaDeviceSynchronize(); }
inline void stream_synchronize(Stream stream) {
    cudaStreamSynchronize(_to_cuda(stream));
}

// --- events ---
inline Event* event_create(bool enable_timing) {
    cudaEvent_t event;
    cudaEventCreateWithFlags(
        &event, enable_timing ? cudaEventDefault : cudaEventDisableTiming);
    return reinterpret_cast<Event*>(event);
}
inline void event_record(Event* event, Stream stream) {
    cudaEventRecord(reinterpret_cast<cudaEvent_t>(event), _to_cuda(stream));
}
inline void event_synchronize(Event* event) {
    cudaEventSynchronize(reinterpret_cast<cudaEvent_t>(event));
}
inline void event_destroy(Event* event) {
    cudaEventDestroy(reinterpret_cast<cudaEvent_t>(event));
}
inline float event_elapsed_ms(Event* start, Event* end) {
    float ms = 0.0f;
    cudaEventElapsedTime(&ms, reinterpret_cast<cudaEvent_t>(start),
                         reinterpret_cast<cudaEvent_t>(end));
    return ms;
}

}  // namespace backend
