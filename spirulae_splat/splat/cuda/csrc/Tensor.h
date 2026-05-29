#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <mutex>
#include <tuple>
#include <variant>
#include <vector>
#include <stdexcept>
#include <string>
#include <unordered_map>


// Use this to communicate with Python, without need to compile/link ATen/LibTorch
typedef std::tuple<
    uint64_t,  // data_ptr
    uint32_t,  // element_size
    std::vector<int64_t>  // shape
> TorchTensorView;


// AsyncReadout<T> — pinned-host ring buffer for non-blocking D->H scalar readouts.
//
// Usage pattern (one iteration of training):
//     static AsyncReadout<float> R(N);
//     const float* prev = R.read_previous();   // sync on prior issue (cheap)
//     LossValues out = prev ? pack(prev) : LossValues{};
//     R.issue(d_src, stream);                  // async D->H of this iter
//     return out;
//
// `read_previous()` returns the host pointer to the most-recently-issued slot
// AFTER synchronizing on its completion event. Because there's a full
// iteration of GPU work between issue() and the next read_previous(), the
// event is essentially always already complete — `cudaEventSynchronize`
// returns immediately, replacing what used to be a blocking `cudaMemcpy`.
//
// First call: read_previous() returns nullptr.
template<typename T>
class AsyncReadout {
public:
    explicit AsyncReadout(int n_elems, int n_slots = 2)
        : _n_elems(n_elems), _n_slots(n_slots),
          _h(n_slots, nullptr), _evt(n_slots),
          _valid(n_slots, false), _next_write(0)
    {
        for (int i = 0; i < _n_slots; ++i) {
            cudaMallocHost(&_h[i], (size_t)_n_elems * sizeof(T));
            cudaEventCreateWithFlags(&_evt[i], cudaEventDisableTiming);
        }
    }

    ~AsyncReadout() {
        for (int i = 0; i < _n_slots; ++i) {
            if (_h[i]) cudaFreeHost(_h[i]);
            cudaEventDestroy(_evt[i]);
        }
    }

    AsyncReadout(const AsyncReadout&) = delete;
    AsyncReadout& operator=(const AsyncReadout&) = delete;

    void issue(const T* d_src, cudaStream_t stream = 0) {
        int i = _next_write;
        cudaMemcpyAsync(_h[i], d_src, (size_t)_n_elems * sizeof(T),
                        cudaMemcpyDeviceToHost, stream);
        cudaEventRecord(_evt[i], stream);
        _valid[i] = true;
        _next_write = (_next_write + 1) % _n_slots;
    }

    // Pointer to the slot most recently filled by issue(), after syncing on
    // its completion event. nullptr if no issue() yet.
    const T* read_previous() {
        int i = (_next_write + _n_slots - 1) % _n_slots;
        if (!_valid[i]) return nullptr;
        cudaEventSynchronize(_evt[i]);
        return _h[i];
    }

    int n_elems() const { return _n_elems; }

private:
    int _n_elems;
    int _n_slots;
    std::vector<T*> _h;
    std::vector<cudaEvent_t> _evt;
    std::vector<bool> _valid;
    int _next_write;
};


// Global CUDA memory pool — keyed named buffers with high-water-mark semantics.
// Memory never shrinks on resize; only reallocates when the new size exceeds capacity.
// freeAll() is the only way to release memory, call it at VRAM pressure or exit.
// All DeviceVector/Tensor2D/3D resize() calls route through here.
class DevicePool {
    struct Slot {
        void*  ptr        = nullptr;
        size_t cap_bytes  = 0;  // allocated
        size_t used_bytes = 0;  // logical (last requested)
    };

    std::unordered_map<std::string, Slot> _slots;
    // Guards _slots. Multiple threads (e.g., trainer + viewer post-processor)
    // call acquire concurrently — unprotected unordered_map operations
    // (operator[] insert, rehash) are UB across threads.
    mutable std::mutex _mu;

    DevicePool() = default;
    DevicePool(const DevicePool&) = delete;
    DevicePool& operator=(const DevicePool&) = delete;

public:
    static DevicePool& global() {
        static DevicePool pool;
        return pool;
    }

    // Return a pointer to at least `n` elements of type T for the given key.
    // Reallocates only when current capacity is exceeded.
    template<typename T>
    T* acquire(const std::string& key, size_t n) {
        std::lock_guard<std::mutex> lock(_mu);
        size_t bytes = n * sizeof(T);
        auto& slot = _slots[key];
        if (bytes > slot.cap_bytes) {
            if (slot.ptr) cudaFree(slot.ptr);
            cudaMalloc(&slot.ptr, bytes);
            slot.cap_bytes = bytes;
        }
        slot.used_bytes = bytes;
        return static_cast<T*>(slot.ptr);
    }

    // Free a single named buffer and remove it from the pool.
    void free(const std::string& key) {
        std::lock_guard<std::mutex> lock(_mu);
        auto it = _slots.find(key);
        if (it != _slots.end()) {
            if (it->second.ptr) cudaFree(it->second.ptr);
            _slots.erase(it);
        }
    }

    // Free all managed CUDA memory.
    void freeAll() {
        std::lock_guard<std::mutex> lock(_mu);
        for (auto& [key, slot] : _slots)
            if (slot.ptr) cudaFree(slot.ptr);
        _slots.clear();
    }

    // Total bytes allocated (capacity, not logical size).
    size_t totalAllocBytes() const {
        std::lock_guard<std::mutex> lock(_mu);
        size_t total = 0;
        for (auto& [k, s] : _slots) total += s.cap_bytes;
        return total;
    }

    // Per-slot breakdown: returns [(key, used_bytes, cap_bytes), ...]
    std::vector<std::tuple<std::string, size_t, size_t>> getBreakdown() const {
        std::lock_guard<std::mutex> lock(_mu);
        std::vector<std::tuple<std::string, size_t, size_t>> result;
        result.reserve(_slots.size());
        for (auto& [k, s] : _slots)
            result.emplace_back(k, s.used_bytes, s.cap_bytes);
        return result;
    }

    ~DevicePool() { freeAll(); }
};


// Single reusable scratch buffer for CUB-style temporary allocations.
// Grows monotonically; never fragmented. Intended for single-stream use.
class DeviceScratch {
    void*  _ptr = nullptr;
    size_t _cap = 0;
    mutable std::mutex _mu;  // protects _ptr/_cap from concurrent threads

    DeviceScratch() = default;
    DeviceScratch(const DeviceScratch&) = delete;
    DeviceScratch& operator=(const DeviceScratch&) = delete;

public:
    static DeviceScratch& global() {
        static DeviceScratch s;
        return s;
    }

    // Returns a pointer to at least `bytes` bytes of device scratch space.
    // Note: callers using the returned pointer asynchronously across threads must
    // synchronize externally; cudaFree's implicit device sync inside (when
    // growing) ensures any in-flight kernel using the previous buffer completes
    // before the new allocation.
    void* acquire(size_t bytes) {
        std::lock_guard<std::mutex> lock(_mu);
        if (bytes > _cap) {
            if (_ptr) cudaFree(_ptr);
            cudaMalloc(&_ptr, bytes);
            _cap = bytes;
        }
        return _ptr;
    }

    size_t capBytes() const {
        std::lock_guard<std::mutex> lock(_mu);
        return _cap;
    }

    void freeAll() {
        std::lock_guard<std::mutex> lock(_mu);
        if (_ptr) { cudaFree(_ptr); _ptr = nullptr; _cap = 0; }
    }

    ~DeviceScratch() { freeAll(); }
};


// Free all device memory managed by the pool and scratch buffer.
// Safe to call at VRAM pressure or before process exit.
inline void freeAllDeviceMemory() {
    DevicePool::global().freeAll();
    DeviceScratch::global().freeAll();
}


// 1D device vector — non-owning view; memory owned by DevicePool or PyTorch.
template<typename T>
class DeviceVector {
protected:

    int64_t _numel;
    T* __restrict__ _data_ptr;

public:

    int64_t size() const {
        return _numel;
    }

    T* data_ptr() const {
        return _data_ptr;
    }

    DeviceVector() : _data_ptr(nullptr), _numel(0) {}

    DeviceVector(const TorchTensorView& view) {
        _data_ptr = (T*)std::get<0>(view);
        uint32_t element_size = std::get<1>(view);
        std::vector<int64_t> shape = std::get<2>(view);
        if (shape.size() == 1)
            shape.push_back(1);
        if (shape.size() != 2)
            throw std::runtime_error("Expected 2D tensor view with channel-last layout");
        if (shape[1] * element_size != sizeof(T))
            throw std::runtime_error("Element size mismatch");
        _numel = shape[0];
    }

    // shallow copy
    DeviceVector& operator=(const DeviceVector&) = default;

    // Pool-backed resize: memory owned by DevicePool, this is a non-owning view.
    void resize(const std::string& key, int64_t numel) {
        _data_ptr = DevicePool::global().acquire<T>(key, (size_t)numel);
        _numel = numel;
    }

    // fill zero
    void zero() {
        if (_data_ptr)
            cudaMemset(_data_ptr, 0, _numel * sizeof(T));
    }

#ifdef __CUDACC__
    __device__ T load(int64_t idx) const {
        return _data_ptr[idx];
    }
    __device__ void store(int64_t idx, T value) {
        _data_ptr[idx] = value;
    }
    __device__ T atomicAdd(int64_t idx, T value) {
        return ::atomicAdd(&_data_ptr[idx], value);
    }
    __device__ T atomicMax(int64_t idx, T value) {
        return ::atomicMax(&_data_ptr[idx], value);
    }
    __device__ T atomicMin(int64_t idx, T value) {
        return ::atomicMin(&_data_ptr[idx], value);
    }
#endif

};


// 2D device tensor — non-owning view; memory owned by DevicePool or PyTorch.
template<typename T>
class DeviceTensor2D {
protected:

    int64_t _numel_0;
    int64_t _numel_1;
    T* __restrict__ _data_ptr;

public:

    template<int dim>
    int64_t size() const {
        static_assert(dim == 0 || dim == 1, "Invalid dimension");
        if constexpr (dim == 0) return _numel_0;
        else if constexpr (dim == 1) return _numel_1;
    }

    int64_t numel() const {
        return _numel_0 * _numel_1;
    }

    T* data_ptr() const {
        return _data_ptr;
    }

    DeviceTensor2D() : _data_ptr(nullptr), _numel_0(0), _numel_1(0) {}

    DeviceTensor2D(const TorchTensorView& view) {
        _data_ptr = (T*)std::get<0>(view);
        uint32_t element_size = std::get<1>(view);
        std::vector<int64_t> shape = std::get<2>(view);
        if (shape.size() == 2)
            shape.push_back(1);
        if (shape.size() != 3)
            throw std::runtime_error("Expected 3D tensor view with channel-last layout");
        if (shape[2] * element_size != sizeof(T))
            throw std::runtime_error("Element size mismatch");
        _numel_0 = shape[0];
        _numel_1 = shape[1];
    }

    // shallow copy
    DeviceTensor2D& operator=(const DeviceTensor2D&) = default;

    // Pool-backed resize: memory owned by DevicePool, this is a non-owning view.
    void resize(const std::string& key, int64_t numel_0, int64_t numel_1) {
        _data_ptr = DevicePool::global().acquire<T>(key, (size_t)(numel_0 * numel_1));
        _numel_0 = numel_0;
        _numel_1 = numel_1;
    }

    // fill zero
    void zero() {
        if (_data_ptr)
            cudaMemset(_data_ptr, 0, _numel_0 * _numel_1 * sizeof(T));
    }

#ifdef __CUDACC__
    __device__ T load(int64_t idx0, int64_t idx1) const {
        return _data_ptr[idx0 * _numel_1 + idx1];
    }
    __device__ void store(int64_t idx0, int64_t idx1, T value) {
        _data_ptr[idx0 * _numel_1 + idx1] = value;
    }
    __device__ T atomicAdd(int64_t idx0, int64_t idx1, T value) {
        return ::atomicAdd(&_data_ptr[idx0 * _numel_1 + idx1], value);
    }
    __device__ T atomicMax(int64_t idx0, int64_t idx1, T value) {
        return ::atomicMax(&_data_ptr[idx0 * _numel_1 + idx1], value);
    }
    __device__ T atomicMin(int64_t idx0, int64_t idx1, T value) {
        return ::atomicMin(&_data_ptr[idx0 * _numel_1 + idx1], value);
    }
#endif

};


// 3D device tensor — non-owning view; memory owned by DevicePool or PyTorch.
template<typename T>
class DeviceTensor3D {
protected:

    int64_t _numel_0;
    int64_t _numel_1;
    int64_t _numel_2;
    T* __restrict__ _data_ptr;

public:

    template<int dim>
    int64_t size() const {
        static_assert(dim == 0 || dim == 1 || dim == 2, "Invalid dimension");
        if constexpr (dim == 0) return _numel_0;
        else if constexpr (dim == 1) return _numel_1;
        else if constexpr (dim == 2) return _numel_2;
    }

    int64_t numel() const {
        return _numel_0 * _numel_1 * _numel_2;
    }

    T* data_ptr() const {
        return _data_ptr;
    }

    DeviceTensor3D() : _data_ptr(nullptr), _numel_0(0), _numel_1(0), _numel_2(0) {}

    DeviceTensor3D(const TorchTensorView& view) {
        _data_ptr = (T*)std::get<0>(view);
        uint32_t element_size = std::get<1>(view);
        std::vector<int64_t> shape = std::get<2>(view);
        if (shape.size() == 3)
            shape.push_back(1);
        if (shape.size() != 4)
            throw std::runtime_error("Expected 4D tensor view with channel-last layout");
        if (shape[3] * element_size != sizeof(T))
            throw std::runtime_error("Element size mismatch");
        _numel_0 = shape[0];
        _numel_1 = shape[1];
        _numel_2 = shape[2];
    }

    // shallow copy
    DeviceTensor3D& operator=(const DeviceTensor3D&) = default;

    // Pool-backed resize: memory owned by DevicePool, this is a non-owning view.
    void resize(const std::string& key, int64_t numel_0, int64_t numel_1, int64_t numel_2) {
        _data_ptr = DevicePool::global().acquire<T>(key, (size_t)(numel_0 * numel_1 * numel_2));
        _numel_0 = numel_0;
        _numel_1 = numel_1;
        _numel_2 = numel_2;
    }

    // fill zero
    void zero() {
        if (_data_ptr)
            cudaMemset(_data_ptr, 0, _numel_0 * _numel_1 * _numel_2 * sizeof(T));
    }

#ifdef __CUDACC__
    __device__ T load(int64_t idx0, int64_t idx1, int64_t idx2) const {
        return _data_ptr[idx0 * _numel_1 * _numel_2 + idx1 * _numel_2 + idx2];
    }
    __device__ void store(int64_t idx0, int64_t idx1, int64_t idx2, T value) {
        _data_ptr[idx0 * _numel_1 * _numel_2 + idx1 * _numel_2 + idx2] = value;
    }
    __device__ T atomicAdd(int64_t idx0, int64_t idx1, int64_t idx2, T value) {
        return ::atomicAdd(&_data_ptr[idx0 * _numel_1 * _numel_2 + idx1 * _numel_2 + idx2], value);
    }
    __device__ T atomicMax(int64_t idx0, int64_t idx1, int64_t idx2, T value) {
        return ::atomicMax(&_data_ptr[idx0 * _numel_1 * _numel_2 + idx1 * _numel_2 + idx2], value);
    }
    __device__ T atomicMin(int64_t idx0, int64_t idx1, int64_t idx2, T value) {
        return ::atomicMin(&_data_ptr[idx0 * _numel_1 * _numel_2 + idx1 * _numel_2 + idx2], value);
    }
#endif

};


// 5D device tensor — non-owning view; memory owned by DevicePool or PyTorch.
// Channel-last expected (last dim packs into sizeof(T)); fits bilagrid [N,C,L,H,W]
// when each cell is a single float.
template<typename T>
class DeviceTensor5D {
protected:

    int64_t _numel_0;
    int64_t _numel_1;
    int64_t _numel_2;
    int64_t _numel_3;
    int64_t _numel_4;
    T* __restrict__ _data_ptr;

public:

    template<int dim>
    int64_t size() const {
        static_assert(dim >= 0 && dim < 5, "Invalid dimension");
        if constexpr (dim == 0) return _numel_0;
        else if constexpr (dim == 1) return _numel_1;
        else if constexpr (dim == 2) return _numel_2;
        else if constexpr (dim == 3) return _numel_3;
        else if constexpr (dim == 4) return _numel_4;
    }

    int64_t numel() const {
        return _numel_0 * _numel_1 * _numel_2 * _numel_3 * _numel_4;
    }

    T* data_ptr() const {
        return _data_ptr;
    }

    DeviceTensor5D() : _data_ptr(nullptr),
        _numel_0(0), _numel_1(0), _numel_2(0), _numel_3(0), _numel_4(0) {}

    DeviceTensor5D(const TorchTensorView& view) {
        _data_ptr = (T*)std::get<0>(view);
        uint32_t element_size = std::get<1>(view);
        std::vector<int64_t> shape = std::get<2>(view);
        if (shape.size() == 5)
            shape.push_back(1);
        if (shape.size() != 6)
            throw std::runtime_error("Expected 6D tensor view with channel-last layout");
        if (shape[5] * element_size != sizeof(T))
            throw std::runtime_error("Element size mismatch");
        _numel_0 = shape[0];
        _numel_1 = shape[1];
        _numel_2 = shape[2];
        _numel_3 = shape[3];
        _numel_4 = shape[4];
    }

    // shallow copy
    DeviceTensor5D& operator=(const DeviceTensor5D&) = default;

    // Pool-backed resize: memory owned by DevicePool, this is a non-owning view.
    void resize(const std::string& key, int64_t n0, int64_t n1, int64_t n2, int64_t n3, int64_t n4) {
        _data_ptr = DevicePool::global().acquire<T>(key, (size_t)(n0 * n1 * n2 * n3 * n4));
        _numel_0 = n0; _numel_1 = n1; _numel_2 = n2; _numel_3 = n3; _numel_4 = n4;
    }

    // fill zero
    void zero() {
        if (_data_ptr)
            cudaMemset(_data_ptr, 0, _numel_0 * _numel_1 * _numel_2 * _numel_3 * _numel_4 * sizeof(T));
    }

#ifdef __CUDACC__
    __device__ T load(int64_t i0, int64_t i1, int64_t i2, int64_t i3, int64_t i4) const {
        int64_t stride_3 = _numel_4;
        int64_t stride_2 = _numel_3 * stride_3;
        int64_t stride_1 = _numel_2 * stride_2;
        int64_t stride_0 = _numel_1 * stride_1;
        return _data_ptr[i0 * stride_0 + i1 * stride_1 + i2 * stride_2 + i3 * stride_3 + i4];
    }
    __device__ void store(int64_t i0, int64_t i1, int64_t i2, int64_t i3, int64_t i4, T value) {
        int64_t stride_3 = _numel_4;
        int64_t stride_2 = _numel_3 * stride_3;
        int64_t stride_1 = _numel_2 * stride_2;
        int64_t stride_0 = _numel_1 * stride_1;
        _data_ptr[i0 * stride_0 + i1 * stride_1 + i2 * stride_2 + i3 * stride_3 + i4] = value;
    }
#endif

};


// ND device tensor — non-owning view.
// Intended for heterogeneous containers (e.g. std::vector of tensors with different ranks).
// For fixed-rank device-friendly code, prefer DeviceVector/Tensor2D/Tensor3D.
template<typename T>
class DeviceTensorND {
protected:

    int64_t ndim;
    std::vector<int64_t> shape;
    T* __restrict__ _data_ptr;

public:

    int64_t size(int dim) const {
        if (dim < 0 || dim >= ndim)
            throw std::runtime_error("Invalid dimension");
        return shape[dim];
    }

    int64_t numel() const {
        int64_t n = 1;
        for (int i = 0; i < ndim; i++) n *= shape[i];
        return n;
    }

    T* data_ptr() const {
        return _data_ptr;
    }

    DeviceTensorND() : _data_ptr(nullptr), ndim(0) {}

    DeviceTensorND(const TorchTensorView& view) {
        _data_ptr = (T*)std::get<0>(view);
        uint32_t element_size = std::get<1>(view);
        std::vector<int64_t> s = std::get<2>(view);
        if (s.empty())
            throw std::runtime_error("Expected non-empty tensor view");
        if (s.back() * element_size != sizeof(T))
            throw std::runtime_error("Element size mismatch");
        ndim = (int64_t)s.size() - 1;
        shape.assign(s.begin(), s.begin() + ndim);
    }

    DeviceTensorND(const DeviceVector<T>& vec) {
        _data_ptr = vec.data_ptr();
        ndim = 1;
        shape = {vec.size()};
    }

    DeviceTensorND(const DeviceTensor2D<T>& tensor) {
        _data_ptr = tensor.data_ptr();
        ndim = 2;
        shape = {tensor.template size<0>(), tensor.template size<1>()};
    }

    DeviceTensorND(const DeviceTensor3D<T>& tensor) {
        _data_ptr = tensor.data_ptr();
        ndim = 3;
        shape = {tensor.template size<0>(), tensor.template size<1>(), tensor.template size<2>()};
    }
};


class DeviceTensorFloatND : public DeviceTensorND<float> {
public:
    using DeviceTensorND::DeviceTensorND;

    DeviceTensorFloatND(const DeviceVector<float2>& vec) {
        _data_ptr = (float*)vec.data_ptr();
        ndim = 2;
        shape = {vec.size(), 2};
    }
    DeviceTensorFloatND(const DeviceTensor2D<float2>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 3;
        shape = {tensor.template size<0>(), tensor.template size<1>(), 2};
    }
    DeviceTensorFloatND(const DeviceTensor3D<float2>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 4;
        shape = {tensor.template size<0>(), tensor.template size<1>(), tensor.template size<2>(), 2};
    }

    DeviceTensorFloatND(const DeviceVector<float3>& vec) {
        _data_ptr = (float*)vec.data_ptr();
        ndim = 2;
        shape = {vec.size(), 3};
    }
    DeviceTensorFloatND(const DeviceTensor2D<float3>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 3;
        shape = {tensor.template size<0>(), tensor.template size<1>(), 3};
    }
    DeviceTensorFloatND(const DeviceTensor3D<float3>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 4;
        shape = {tensor.template size<0>(), tensor.template size<1>(), tensor.template size<2>(), 3};
    }

    DeviceTensorFloatND(const DeviceVector<float4>& vec) {
        _data_ptr = (float*)vec.data_ptr();
        ndim = 2;
        shape = {vec.size(), 4};
    }
    DeviceTensorFloatND(const DeviceTensor2D<float4>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 3;
        shape = {tensor.template size<0>(), tensor.template size<1>(), 4};
    }
    DeviceTensorFloatND(const DeviceTensor3D<float4>& tensor) {
        _data_ptr = (float*)tensor.data_ptr();
        ndim = 4;
        shape = {tensor.template size<0>(), tensor.template size<1>(), tensor.template size<2>(), 4};
    }
};
