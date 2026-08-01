#include "nn/Tensor.h"

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/vk/Stream.h"

#include <cmath>
#include <vector>

namespace nn {

namespace {
int infer_ndim(int ndim, int64_t d1, int64_t d2, int64_t d3) {
    if (ndim > 0) return ndim;
    if (d3 != 1) return 4;
    if (d2 != 1) return 3;
    if (d1 != 1) return 2;
    return 1;
}
}  // namespace

Tensor arena_tensor(vk::Arena& a, DType t, int64_t d0, int64_t d1, int64_t d2, int64_t d3,
                    int ndim) {
    Tensor out;
    out.dtype = t;
    out.ndim = infer_ndim(ndim, d1, d2, d3);
    out.shape[0] = d0; out.shape[1] = d1; out.shape[2] = d2; out.shape[3] = d3;
    out.ptr = a.alloc((vk::DevicePtr)out.numel() * dtype_size(t));
    return out;
}

Tensor pool_tensor(vk::PoolSlot slot, uint32_t sub, DType t, int64_t d0, int64_t d1,
                   int64_t d2, int64_t d3, int ndim) {
    Tensor out;
    out.dtype = t;
    out.ndim = infer_ndim(ndim, d1, d2, d3);
    out.shape[0] = d0; out.shape[1] = d1; out.shape[2] = d2; out.shape[3] = d3;
    out.ptr = vk::VramPool::get().acquire(slot, sub, (VkDeviceSize)out.bytes());
    return out;
}

void tensor_to_host(const Tensor& t, float* dst, int64_t count) {
    NN_CHECK(t.valid(), "tensor_to_host on an empty tensor");
    NN_CHECK(count <= t.numel(), "tensor_to_host: %lld > %lld elements",
               (long long)count, (long long)t.numel());
    if (t.dtype == DType::F32) {
        vk::Stream::get().download(dst, t.ptr, (VkDeviceSize)count * 4);
    } else if (t.dtype == DType::F16) {
        std::vector<uint16_t> tmp((size_t)count);
        vk::Stream::get().download(tmp.data(), t.ptr, (VkDeviceSize)count * 2);
        for (int64_t i = 0; i < count; ++i) dst[i] = half_to_float(tmp[(size_t)i]);
    } else if (t.dtype == DType::I32) {
        std::vector<int32_t> tmp((size_t)count);
        vk::Stream::get().download(tmp.data(), t.ptr, (VkDeviceSize)count * 4);
        for (int64_t i = 0; i < count; ++i) dst[i] = (float)tmp[(size_t)i];
    } else {
        std::vector<uint8_t> tmp((size_t)count);
        vk::Stream::get().download(tmp.data(), t.ptr, (VkDeviceSize)count);
        for (int64_t i = 0; i < count; ++i) dst[i] = (float)tmp[(size_t)i];
    }
}

void tensor_from_host(const Tensor& t, const float* src, int64_t count) {
    NN_CHECK(t.valid(), "tensor_from_host on an empty tensor");
    NN_CHECK(count <= t.numel(), "tensor_from_host: %lld > %lld elements",
               (long long)count, (long long)t.numel());
    if (t.dtype == DType::F32) {
        vk::Stream::get().upload(t.ptr, src, (VkDeviceSize)count * 4);
    } else if (t.dtype == DType::F16) {
        std::vector<uint16_t> tmp((size_t)count);
        for (int64_t i = 0; i < count; ++i) tmp[(size_t)i] = float_to_half(src[i]);
        vk::Stream::get().upload(t.ptr, tmp.data(), (VkDeviceSize)count * 2);
    } else {
        fail("tensor_from_host: unsupported dtype %s", dtype_name(t.dtype));
    }
}

void tensor_debug_dump(const char* label, const Tensor& t, int max_values) {
    if (log_level() < 3 || !t.valid()) return;
    const int64_t n = t.numel();
    std::vector<float> v((size_t)n);
    tensor_to_host(t, v.data(), n);
    double sum = 0, sq = 0, mn = 1e30, mx = -1e30;
    int64_t nan_count = 0;
    for (int64_t i = 0; i < n; ++i) {
        float x = v[(size_t)i];
        if (std::isnan(x) || std::isinf(x)) { ++nan_count; continue; }
        sum += x; sq += (double)x * x;
        mn = std::fmin(mn, x); mx = std::fmax(mx, x);
    }
    double mean = sum / (double)n;
    NN_LOG_DEBUG("[dbg] %-28s [", label);
    for (int i = 0; i < t.ndim; ++i)
        NN_LOG_DEBUG("%lld%s", (long long)t.shape[i], i + 1 < t.ndim ? "," : "");
    NN_LOG_DEBUG("] %s mean=%.5f std=%.5f min=%.5f max=%.5f%s\n      ",
                   dtype_name(t.dtype), mean,
                   std::sqrt(std::fmax(sq / (double)n - mean * mean, 0.0)), mn, mx,
                   nan_count ? " HAS-NAN" : "");
    for (int i = 0; i < max_values && i < n; ++i)
        NN_LOG_DEBUG("%.5f ", v[(size_t)i]);
    NN_LOG_DEBUG("\n");
}

}  // namespace nn
