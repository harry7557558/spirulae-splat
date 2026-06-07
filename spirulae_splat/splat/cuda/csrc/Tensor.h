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


// AsyncReadout<T> -- pinned-host ring buffer for non-blocking D->H scalar readouts.
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
// event is essentially always already complete -- `cudaEventSynchronize`
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


// Global CUDA memory pool -- keyed named buffers with high-water-mark semantics.
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
    // call acquire concurrently -- unprotected unordered_map operations
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


// 1D device vector -- non-owning view; memory owned by DevicePool or PyTorch.
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


// 2D device tensor -- non-owning view; memory owned by DevicePool or PyTorch.
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

    // Shape-only descriptor: keeps the [n0, n1] dims but holds no device
    // memory. Used when the canonical storage moves to a packed/quantized
    // buffer elsewhere but downstream code still expects shape info (e.g.,
    // TensorArray's per-axis stride inference).
    void set_shape_no_alloc(int64_t numel_0, int64_t numel_1) {
        _data_ptr = nullptr;
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


// 3D device tensor -- non-owning view; memory owned by DevicePool or PyTorch.
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


// 5D device tensor -- non-owning view; memory owned by DevicePool or PyTorch.
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


// ND device tensor -- non-owning view.
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

    int64_t get_ndim() const {
        return ndim;
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


// ============================================================================
// QuantizedAdamState -- joint (u, log_s) Adam optimizer-state quantization.
//
// Per cell the codec stores TWO linearly-quantized "primitive" values:
//     primitives.x = u     = g1 / (sqrt(g2) + eps)
//     primitives.y = log_s = log1p(sqrt(g2) / eps)
// log1p(x/eps) is non-negative, monotone, sends 0 -> 0 (so the all-zero
// initial state of `packed` + `bounds` correctly decodes to sqrt_g2 = 0,
// matching the float-path zero-init invariant), and is approximately
// log(sqrt_g2 / eps) for sqrt_g2 >> eps. Its inverse is eps · expm1(log_s).
//
// Decode path:
//     u, log_s = linear-dequantize the two bytes within the per-block bounds
//     sqrt_g2  = eps · expm1(log_s)
//     g1       = u · (sqrt_g2 + eps)
//     g2       = sqrt_g2 · sqrt_g2
// Encode path is the same transformation reversed; the block-level min/max
// reduction inside the kernels operates on the linearly-quantized primitives
// (`u`, `log_s`), so the per-block bounds float4 means:
//     bounds.x, bounds.y = (u_min, u_max)             -- linear u domain
//     bounds.z, bounds.w = (log_s_min, log_s_max)     -- log-mapped sqrt(g2) domain
//
// Why log_s instead of plain sqrt(g2)? Empirically (see
// scripts/investigate_quantization.py "SNR(u_next)" column) the next-step
// Adam update u' = β1·g1 / sqrt(β2·g2 + ...) is dominated by errors in
// sqrt(g2) at small magnitudes -- exactly the region where linear quantization
// has the worst RELATIVE precision and where Adam's normalizer is most
// sensitive. Storing log1p(sqrt_g2 / eps) gives uniform relative resolution
// across all magnitudes within a block and consistently buys 20--30 dB on the
// post-step SNR metric at no storage cost.
//
// Storage is array-of-structure (AoS) interleaved over the two primitives:
//   BITS=16: 4 bytes/cell -- (uint16_t u_q) || (uint16_t log_s_q)
//   BITS=8:  2 bytes/cell -- byte[idx*2 + 0] = u_q, byte[idx*2 + 1] = log_s_q
//   BITS=4:  1 byte/cell  -- low nibble = u_q,    high nibble = log_s_q
// Quantization is endpoint-exact: q=0 -> lo, q=(kLevels-1) -> hi exactly.
//
// Hyperparameters (bit depth, block size, eps) are compile-time template
// arguments; all the math constants are constexpr so the codec compiles to
// straight-line arithmetic plus log1pf / expm1f (single CUDA math-lib calls).
// ============================================================================

template<int BITS, int BLOCK_SIZE = 256>
class QuantizedAdamState {
public:
    static_assert(BITS == 4 || BITS == 8 || BITS == 16,
                  "QuantizedAdamState: only 4-bit, 8-bit, and 16-bit are supported");
    static_assert(BLOCK_SIZE > 0 && (BLOCK_SIZE & (BLOCK_SIZE - 1)) == 0,
                  "QuantizedAdamState: BLOCK_SIZE must be a positive power of 2");

    static constexpr int   kBits         = BITS;
    static constexpr int   kBlockSize    = BLOCK_SIZE;
    static constexpr int   kLevels       = 1 << BITS;
    static constexpr float kQMax         = (float)(kLevels - 1);
    static constexpr float kInvQMax      = 1.0f / kQMax;
    static constexpr int   kBytesPerCell = (BITS == 16) ? 4 : (BITS == 8) ? 2 : 1;
    static constexpr float kEps          = 1e-15f;   // matches Adam denominator

    // Storage (non-owning views backed by the global DevicePool).
    DeviceVector<uint8_t> packed;     // size = n_cells * kBytesPerCell
    DeviceVector<float4>  bounds;     // size = n_bounds (caller-chosen)
    int64_t n_cells  = 0;
    int64_t n_bounds = 0;

    // ---- Host: storage management ------------------------------------------

    void resize(const std::string& key_prefix, int64_t cells, int64_t bnds) {
        n_cells  = cells;
        n_bounds = bnds;
        packed.resize(key_prefix + ".q",  (size_t)(cells * kBytesPerCell));
        bounds.resize(key_prefix + ".qb", (size_t)bnds);
    }

    void zero() {
        if (packed.data_ptr()) packed.zero();
        if (bounds.data_ptr()) bounds.zero();
    }

    bool initialized() const { return packed.data_ptr() != nullptr; }

    uint8_t* packed_ptr() const { return packed.data_ptr(); }
    float4*  bounds_ptr() const { return bounds.data_ptr(); }

    int64_t packed_bytes() const { return n_cells * kBytesPerCell; }
    int64_t bounds_bytes() const { return n_bounds * (int64_t)sizeof(float4); }
    int64_t total_bytes()  const { return packed_bytes() + bounds_bytes(); }

    // ---- Device: codec primitives ------------------------------------------
#ifdef __CUDACC__
    // Forward / inverse map for the sqrt_g2 half:
    //   forward:  sqrt_g2 -> log1pf(sqrt_g2 / eps)
    //   inverse:  log_s   -> eps * expm1f(log_s)
    // Both single CUDA math-lib calls; preserve the 0 <-> 0 fixed point so
    // the all-zero initial state of `packed` + `bounds` correctly decodes to
    // sqrt_g2 = 0 (no special init needed).
    __device__ static inline float forward_sqrt_g2(float sqrt_g2) {
        return log1pf(fmaxf(sqrt_g2, 0.0f) * (1.0f / kEps));
    }
    __device__ static inline float inverse_sqrt_g2(float log_s) {
        return kEps * expm1f(log_s);
    }

    // Decode the linearly-quantized primitives (u, log_s) at cell `idx`.
    // Note: returns the INTERNAL representation, not (u, sqrt_g2). Use
    // decode_g1g2 for the reconstructed Adam state.
    __device__ static inline float2 decode_us(
        const uint8_t* __restrict__ packed_ptr, int64_t idx, float4 mm
    ) {
        float u_q, s_q;
        if constexpr (BITS == 16) {
            const uint16_t* p = reinterpret_cast<const uint16_t*>(packed_ptr);
            u_q = (float)p[idx * 2 + 0];
            s_q = (float)p[idx * 2 + 1];
        } else if constexpr (BITS == 8) {
            u_q = (float)packed_ptr[idx * 2 + 0];
            s_q = (float)packed_ptr[idx * 2 + 1];
        } else {  // BITS == 4
            uint8_t b = packed_ptr[idx];
            u_q = (float)(b & 0x0Fu);
            s_q = (float)((b >> 4) & 0x0Fu);
        }
        return float2{
            mm.x + (mm.y - mm.x) * (u_q * kInvQMax),
            mm.z + (mm.w - mm.z) * (s_q * kInvQMax)
        };
    }

    // Decode + reconstruct Adam state (g1, g2) from the cell at `idx`.
    __device__ static inline float2 decode_g1g2(
        const uint8_t* __restrict__ packed_ptr, int64_t idx, float4 mm
    ) {
        float2 prim = decode_us(packed_ptr, idx, mm);
        float sqrt_g2 = inverse_sqrt_g2(prim.y);
        return float2{
            prim.x * (sqrt_g2 + kEps),   // g1 = u * (sqrt_g2 + eps)
            sqrt_g2 * sqrt_g2             // g2 = sqrt_g2 * sqrt_g2
        };
    }

    // Encode the linearly-quantized primitives (u, log_s) into the packed
    // byte stream at cell `idx`. Endpoint-exact 8-bit (or 4-bit) within mm.
    // Kernels typically obtain (u_val, log_s_val) by calling g1g2_to_us and
    // pass mm = the block-reduced bounds.
    __device__ static inline void encode_us(
        uint8_t* __restrict__ packed_ptr, int64_t idx,
        float u_val, float log_s_val, float4 mm
    ) {
        float u_range = fmaxf(mm.y - mm.x, kEps);
        float s_range = fmaxf(mm.w - mm.z, kEps);
        float u_qf = fminf(fmaxf(roundf(kQMax * (u_val     - mm.x) / u_range), 0.0f), kQMax);
        float s_qf = fminf(fmaxf(roundf(kQMax * (log_s_val - mm.z) / s_range), 0.0f), kQMax);
        if constexpr (BITS == 16) {
            uint16_t* p = reinterpret_cast<uint16_t*>(packed_ptr);
            p[idx * 2 + 0] = (uint16_t)u_qf;
            p[idx * 2 + 1] = (uint16_t)s_qf;
        } else if constexpr (BITS == 8) {
            packed_ptr[idx * 2 + 0] = (uint8_t)u_qf;
            packed_ptr[idx * 2 + 1] = (uint8_t)s_qf;
        } else {  // BITS == 4
            uint8_t u_q = (uint8_t)u_qf;
            uint8_t s_q = (uint8_t)s_qf;
            packed_ptr[idx] = (uint8_t)((u_q & 0x0Fu) | ((s_q & 0x0Fu) << 4));
        }
    }

    // Encode (g1, g2) -- internally derives (u, log_s) via g1g2_to_us.
    __device__ static inline void encode_g1g2(
        uint8_t* __restrict__ packed_ptr, int64_t idx,
        float g1, float g2, float4 mm
    ) {
        float2 prim = g1g2_to_us(g1, g2);
        encode_us(packed_ptr, idx, prim.x, prim.y, mm);
    }

    // (g1, g2) -> linearly-quantized primitives (u, log_s). Block reduction
    // inside host kernels operates on these -- both halves are min/max-reduced
    // in their own (linear) domain to produce the per-block float4 bounds.
    __device__ static inline float2 g1g2_to_us(float g1, float g2) {
        float sqrt_g2 = sqrtf(fmaxf(g2, 0.0f));
        float u       = g1 / (sqrt_g2 + kEps);
        return float2{u, forward_sqrt_g2(sqrt_g2)};
    }
#endif // __CUDACC__
};


// ============================================================================
// QuantizedTensorLog -- generic block-wise log-encoded quantized tensor.
//
// Stores a single non-negative scalar per cell as q ~ log1pf(v / eps), with
// per-block (min, max) float2 bounds in the log domain. Endpoint-exact within
// each block: q=0 -> decoded=lo, q=(kLevels-1) -> decoded=hi. Preserves the
// 0 <-> 0 fixed point so an all-zero (`packed` + `bounds`) initial state
// decodes to v = 0 exactly -- matches the "initial_accumulator_value=0" we use
// for AdaGrad accumulators and saves a separate init pass.
//
// Not Adam-specific. Anything that wants relative-precision storage for one
// non-negative value per cell (AdaGrad accum, running variance, occupancy
// likelihood, ...) can reuse it. The bounds float2 always means
// (log-domain min, log-domain max); the codec lives entirely in this class
// and a caller need only:
//     decode_v / encode_v       -- per-cell read / write
//     forward_v / inverse_v     -- to do their own block-reduction in the
//                                  log domain before writing bounds back
//
// Storage (AoS -- same per-cell stride conventions as QuantizedAdamState):
//   BITS=8: 1 byte/cell -- packed[idx] = q
//   BITS=4: 1 byte/2-cells -- low nibble = idx2k, high nibble = idx2k+1
// ============================================================================

template<int BITS, int BLOCK_SIZE = 256>
class QuantizedTensorLog {
public:
    static_assert(BITS == 4 || BITS == 8,
                  "QuantizedTensorLog: only 4-bit and 8-bit are supported");
    static_assert(BLOCK_SIZE > 0 && (BLOCK_SIZE & (BLOCK_SIZE - 1)) == 0,
                  "QuantizedTensorLog: BLOCK_SIZE must be a positive power of 2");

    static constexpr int   kBits         = BITS;
    static constexpr int   kBlockSize    = BLOCK_SIZE;
    static constexpr int   kLevels       = 1 << BITS;
    static constexpr float kQMax         = (float)(kLevels - 1);
    static constexpr float kInvQMax      = 1.0f / kQMax;
    // 8-bit: 1 byte per cell. 4-bit: 0.5 byte per cell (2 cells share a byte).
    // packed_bytes() rounds up to the nearest byte for 4-bit.
    static constexpr float kEps          = 1e-15f;

    // Storage (non-owning views backed by the global DevicePool).
    DeviceVector<uint8_t> packed;     // size = packed_bytes()
    DeviceVector<float2>  bounds;     // size = n_bounds (caller-chosen)
    int64_t n_cells  = 0;
    int64_t n_bounds = 0;

    // ---- Host: storage management ------------------------------------------

    void resize(const std::string& key_prefix, int64_t cells, int64_t bnds) {
        n_cells  = cells;
        n_bounds = bnds;
        packed.resize(key_prefix + ".q",  (size_t)packed_bytes_for(cells));
        bounds.resize(key_prefix + ".qb", (size_t)bnds);
    }

    void zero() {
        if (packed.data_ptr()) packed.zero();
        if (bounds.data_ptr()) bounds.zero();
    }

    bool initialized() const { return packed.data_ptr() != nullptr; }

    uint8_t* packed_ptr() const { return packed.data_ptr(); }
    float2*  bounds_ptr() const { return bounds.data_ptr(); }

    static int64_t packed_bytes_for(int64_t cells) {
        if constexpr (BITS == 8) return cells;
        else                     return (cells + 1) / 2;   // 4-bit, AoS
    }
    int64_t packed_bytes() const { return packed_bytes_for(n_cells); }
    int64_t bounds_bytes() const { return n_bounds * (int64_t)sizeof(float2); }
    int64_t total_bytes()  const { return packed_bytes() + bounds_bytes(); }

    // ---- Device: codec primitives ------------------------------------------
#ifdef __CUDACC__
    // Forward / inverse log map. 0 <-> 0 fixed point.
    __device__ static inline float forward_v(float v) {
        return log1pf(fmaxf(v, 0.0f) * (1.0f / kEps));
    }
    __device__ static inline float inverse_v(float log_v) {
        return kEps * expm1f(log_v);
    }

    // Decode the linearly-quantized log-primitive at cell `idx`. Returns
    // the value in the LOG domain -- use inverse_v(...) to get back to v.
    __device__ static inline float decode_log(
        const uint8_t* __restrict__ packed_ptr, int64_t idx, float2 mm
    ) {
        float q;
        if constexpr (BITS == 8) {
            q = (float)packed_ptr[idx];
        } else {  // BITS == 4
            uint8_t b = packed_ptr[idx >> 1];
            q = (idx & 1) ? (float)((b >> 4) & 0x0Fu) : (float)(b & 0x0Fu);
        }
        return mm.x + (mm.y - mm.x) * (q * kInvQMax);
    }

    // Decode + reconstruct v from the cell at `idx`.
    __device__ static inline float decode_v(
        const uint8_t* __restrict__ packed_ptr, int64_t idx, float2 mm
    ) {
        return inverse_v(decode_log(packed_ptr, idx, mm));
    }

    // Encode a log-domain value into the packed byte stream at cell `idx`.
    // Endpoint-exact within mm. Callers typically obtain `log_v` via
    // forward_v(v) and pass mm = the block-reduced log-domain bounds.
    //
    // NOTE: BITS==4 stores two cells per byte. Concurrent threads writing
    // odd/even partners of the same byte will race. The bilagrid kernels
    // serialize this with a __syncthreads() before the bytewise store --
    // see fused_bilagrid_tv_adagrad_kernel for the pattern.
    __device__ static inline void encode_log(
        uint8_t* __restrict__ packed_ptr, int64_t idx,
        float log_v, float2 mm
    ) {
        float range = fmaxf(mm.y - mm.x, kEps);
        uint8_t q = (uint8_t)fminf(
            fmaxf(roundf(kQMax * (log_v - mm.x) / range), 0.0f), kQMax);
        if constexpr (BITS == 8) {
            packed_ptr[idx] = q;
        } else {  // BITS == 4
            // 4-bit AoS: two cells per byte (idx&1 picks low/high nibble).
            // Caller-synchronized -- see note above.
            int64_t byte_idx = idx >> 1;
            uint8_t old = packed_ptr[byte_idx];
            uint8_t merged = (idx & 1)
                ? (uint8_t)((old & 0x0Fu) | ((q & 0x0Fu) << 4))
                : (uint8_t)((old & 0xF0u) | (q & 0x0Fu));
            packed_ptr[byte_idx] = merged;
        }
    }

    // Encode v directly -- internally calls forward_v.
    __device__ static inline void encode_v(
        uint8_t* __restrict__ packed_ptr, int64_t idx,
        float v, float2 mm
    ) {
        encode_log(packed_ptr, idx, forward_v(v), mm);
    }
#endif // __CUDACC__
};


// ============================================================================
// QuantizedTensor -- generic block-wise UNIFORM (linear) quantized tensor for a
// single scalar value per cell. Sibling of QuantizedTensorLog, but with linear
// (not log) quantization in the original value domain.
//
// Block-wise (min, max) float2 bounds per BLOCK_SIZE cells. Endpoint-exact
// within each block: q=0 -> decoded=lo, q=(kLevels-1) -> decoded=hi.
//
// Targets value (parameter) quantization where the value range is naturally
// bounded and small, so a linear codec captures it well without log-distortion
// of the precision-versus-magnitude trade-off. The codec is the building block
// for future SH-feature / bilagrid-grid value quantization (see "dequantize on
// load" plan): every consumer kernel takes BITS as a template parameter and
// calls decode_v on each read; the encoder writes once after the optimizer
// step alongside the per-block bound reduction.
//
// Storage (AoS):
//   BITS=8:  1 byte/cell  -- packed[idx]      = q
//   BITS=16: 2 bytes/cell -- packed16[idx]    = q   (reinterpret to uint16_t*)
//
// 16-bit is the high-precision option; 8-bit is the compact option. We do not
// support 4-bit here -- linear quantization at 4-bit is too coarse for the
// dynamic range of typical parameter values (use the log codec for AdaGrad-
// style accumulators where 4-bit pays off with log-encoding).
// ============================================================================

template<int BITS, int BLOCK_SIZE = 256>
class QuantizedTensor {
public:
    static_assert(BITS == 8 || BITS == 16,
                  "QuantizedTensor: only 8-bit and 16-bit are supported "
                  "(use QuantizedTensorLog for 4-bit log-domain storage)");
    static_assert(BLOCK_SIZE > 0 && (BLOCK_SIZE & (BLOCK_SIZE - 1)) == 0,
                  "QuantizedTensor: BLOCK_SIZE must be a positive power of 2");

    static constexpr int   kBits         = BITS;
    static constexpr int   kBlockSize    = BLOCK_SIZE;
    static constexpr int   kLevels       = 1 << BITS;
    static constexpr float kQMax         = (float)(kLevels - 1);
    static constexpr float kInvQMax      = 1.0f / kQMax;
    static constexpr int   kBytesPerCell = (BITS == 16) ? 2 : 1;

    // Storage (non-owning views backed by the global DevicePool).
    DeviceVector<uint8_t> packed;     // size = n_cells * kBytesPerCell
    DeviceVector<float2>  bounds;     // size = n_bounds (caller-chosen)
    int64_t n_cells  = 0;
    int64_t n_bounds = 0;

    // ---- Host: storage management ------------------------------------------

    void resize(const std::string& key_prefix, int64_t cells, int64_t bnds) {
        n_cells  = cells;
        n_bounds = bnds;
        packed.resize(key_prefix + ".q",  (size_t)(cells * kBytesPerCell));
        bounds.resize(key_prefix + ".qb", (size_t)bnds);
    }

    void zero() {
        if (packed.data_ptr()) packed.zero();
        if (bounds.data_ptr()) bounds.zero();
    }

    bool initialized() const { return packed.data_ptr() != nullptr; }

    uint8_t* packed_ptr() const { return packed.data_ptr(); }
    float2*  bounds_ptr() const { return bounds.data_ptr(); }

    int64_t packed_bytes() const { return n_cells * (int64_t)kBytesPerCell; }
    int64_t bounds_bytes() const { return n_bounds * (int64_t)sizeof(float2); }
    int64_t total_bytes()  const { return packed_bytes() + bounds_bytes(); }

    // ---- Device: codec primitives ------------------------------------------
#ifdef __CUDACC__
    // Decode the linearly-quantized value at cell `idx` directly (no log).
    __device__ static inline float decode_v(
        const uint8_t* __restrict__ packed_ptr, int64_t idx, float2 mm
    ) {
        float q;
        if constexpr (BITS == 8) {
            q = (float)packed_ptr[idx];
        } else {  // BITS == 16
            const uint16_t* p16 = reinterpret_cast<const uint16_t*>(packed_ptr);
            q = (float)p16[idx];
        }
        return mm.x + (mm.y - mm.x) * (q * kInvQMax);
    }

    // Encode v into the packed byte stream at cell `idx`. Endpoint-exact
    // within (mm.x, mm.y). Callers pass mm = the block-reduced bounds.
    __device__ static inline void encode_v(
        uint8_t* __restrict__ packed_ptr, int64_t idx,
        float v, float2 mm
    ) {
        float range = fmaxf(mm.y - mm.x, 1e-30f);
        float qf = roundf(kQMax * (v - mm.x) / range);
        qf = fminf(fmaxf(qf, 0.0f), kQMax);
        if constexpr (BITS == 8) {
            packed_ptr[idx] = (uint8_t)qf;
        } else {  // BITS == 16
            uint16_t* p16 = reinterpret_cast<uint16_t*>(packed_ptr);
            p16[idx] = (uint16_t)qf;
        }
    }
#endif // __CUDACC__
};
