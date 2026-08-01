// Launchers for the elementwise, resample, gather and shuffle kernels.

#include "nn/core/Error.h"
#include "nn/Ops.h"
#include "nn/vk/Stream.h"

#include <cstring>

namespace nn {

namespace {

struct BinaryParams {
    uint64_t out, a, b;
    uint32_t n, cols;
    float    alpha, beta;
    uint32_t groups_per_row;
};

struct UnaryParams {
    uint64_t out, x;
    uint32_t n;
    float    pre_scale, pre_bias, post_scale, post_bias;
    uint32_t groups_per_row;
};

struct ConvertParams {
    uint64_t out, x;
    uint32_t n, groups_per_row;
};

struct GatherParams {
    uint64_t out, table, ids;
    uint32_t n, cols, vocab, groups_per_row, _pad0;
};

struct StridedCopyParams {
    uint64_t out, x;
    uint32_t rows, cols, in_stride, out_stride, groups_per_row, _pad0;
};

struct WindowParams {
    uint64_t out, x;
    uint32_t H, W, C, ws, nwh, nww, groups_per_row, _pad0;
};

struct TiledAddParams {
    uint64_t out, x, tile;
    uint32_t H, W, C, th, tw, groups_per_row;
};

struct RoiAlignParams {
    uint64_t out, feat, boxes;
    uint32_t nboxes, H, W, C, S, groups_per_row, _pad0;
};

struct ResizeParams {
    uint64_t out, x;
    uint32_t Ho, Wo, Hi, Wi, C, groups_per_row;
};

struct MaskExportParams {
    uint64_t out, x;
    uint32_t Ho, Wo, Hi, Wi;
    float    threshold;
    uint32_t groups_per_row;
};

}  // namespace

// ================
// Elementwise
// ================

void fill(const Tensor& t, float value) {
    NN_CHECK(t.dtype == DType::F32 || value == 0.0f,
               "fill: only zero-fill is defined for %s tensors", dtype_name(t.dtype));
    uint32_t word;
    std::memcpy(&word, &value, 4);
    vk::Stream::get().fill(t.ptr, word, (VkDeviceSize)t.bytes());
}

void copy(const Tensor& dst, const Tensor& src) {
    NN_CHECK(dst.numel() == src.numel(), "copy: %lld vs %lld elements",
               (long long)dst.numel(), (long long)src.numel());
    if (dst.dtype == src.dtype) {
        vk::Stream::get().copy(dst.ptr, src.ptr, (VkDeviceSize)dst.bytes());
        return;
    }
    ConvertParams p{};
    p.out = dst.ptr;
    p.x = src.ptr;
    p.n = (uint32_t)src.numel();
    // kOp selects the direction: 0 -> f32 out, 1 -> packed f16 out.
    const uint32_t to_f16 = (dst.dtype == DType::F16) ? 1u : 0u;
    NN_CHECK(to_f16 || dst.dtype == DType::F32, "copy: unsupported target dtype %s",
               dtype_name(dst.dtype));
    vk::SpecList spec{to_f16, 0u, 0u, (uint32_t)(src.dtype == DType::F16)};
    const int64_t threads = to_f16 ? (src.numel() + 1) / 2 : src.numel();
    vk::Stream::get().dispatchFlat("elementwise.convert", spec, threads, 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

static void binary_op(const Tensor& out, const Tensor& a, const Tensor& b, uint32_t op,
                      float alpha, float beta, Act act) {
    BinaryParams p{};
    p.out = out.ptr;
    p.a = a.ptr;
    p.b = vk::or_fallback(b.ptr);
    p.n = (uint32_t)out.numel();
    p.cols = (uint32_t)out.cols();
    p.alpha = alpha;
    p.beta = beta;

    // Broadcast mode is inferred from the shape, which is unambiguous here:
    // matching element count is elementwise, a single row is per-column, a
    // single element is scalar.
    uint32_t bcast = 0;
    if (b.valid()) {
        if (b.numel() == out.numel())      bcast = 0;
        else if (b.numel() == out.cols())  bcast = 1;
        else if (b.numel() == 1)           bcast = 2;
        else fail("binary op: b has %lld elements, which broadcasts against neither "
                  "%lld (elementwise) nor %lld (per-column)",
                  (long long)b.numel(), (long long)out.numel(), (long long)out.cols());
    } else {
        bcast = 2;  // zeroed fallback
    }
    vk::SpecList spec{op, bcast, (uint32_t)act,
                      (uint32_t)(a.dtype == DType::F16)};
    vk::Stream::get().dispatchFlat("elementwise.binary", spec, out.numel(), 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

void add(const Tensor& out, const Tensor& a, const Tensor& b, float alpha, float beta,
         Act act) {
    binary_op(out, a, b, 0, alpha, beta, act);
}

void mul(const Tensor& out, const Tensor& a, const Tensor& b, Act act) {
    binary_op(out, a, b, 1, 1.0f, 1.0f, act);
}

void unary(const Tensor& out, const Tensor& x, Act act, float pre_scale, float pre_bias,
           float post_scale, float post_bias) {
    UnaryParams p{};
    p.out = out.ptr;
    p.x = x.ptr;
    p.n = (uint32_t)out.numel();
    p.pre_scale = pre_scale;
    p.pre_bias = pre_bias;
    p.post_scale = post_scale;
    p.post_bias = post_bias;
    vk::SpecList spec{0u, 0u, (uint32_t)act, (uint32_t)(x.dtype == DType::F16)};
    vk::Stream::get().dispatchFlat("elementwise.unary", spec, out.numel(), 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

// ================
// Gather / shuffle
// ================

void gather_rows(const Tensor& out, const Tensor& table, const Tensor& ids) {
    GatherParams p{};
    p.out = out.ptr;
    p.table = table.ptr;
    p.ids = ids.ptr;
    p.n = (uint32_t)out.rows();
    p.cols = (uint32_t)out.cols();
    p.vocab = (uint32_t)table.rows();
    vk::SpecList spec{(uint32_t)(table.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat("misc.gather_rows", spec, out.numel(), 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

void strided_copy(const Tensor& out, const Tensor& in, int64_t rows, int64_t cols,
                  int64_t in_stride, int64_t out_stride) {
    StridedCopyParams p{};
    p.out = out.ptr;
    p.x = in.ptr;
    p.rows = (uint32_t)rows;
    p.cols = (uint32_t)cols;
    p.in_stride = (uint32_t)in_stride;
    p.out_stride = (uint32_t)out_stride;
    vk::SpecList spec{(uint32_t)(in.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat("misc.strided_copy", spec, rows * cols, 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

static void window_op(const char* entry, const Tensor& out, const Tensor& in, int H,
                      int W, int C, int ws, int64_t total) {
    WindowParams p{};
    p.out = out.ptr;
    p.x = in.ptr;
    p.H = (uint32_t)H;
    p.W = (uint32_t)W;
    p.C = (uint32_t)C;
    p.ws = (uint32_t)ws;
    p.nwh = (uint32_t)((H + ws - 1) / ws);
    p.nww = (uint32_t)((W + ws - 1) / ws);
    vk::SpecList spec{(uint32_t)(in.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat(entry, spec, total, 256, &p, sizeof(p),
                                   &p.groups_per_row);
}

void window_partition(const Tensor& out, const Tensor& in, int H, int W, int C, int ws) {
    const int64_t nwh = (H + ws - 1) / ws, nww = (W + ws - 1) / ws;
    window_op("misc.window_partition", out, in, H, W, C, ws,
              nwh * nww * ws * ws * (int64_t)C);
}

void window_unpartition(const Tensor& out, const Tensor& in, int H, int W, int C,
                        int ws) {
    window_op("misc.window_unpartition", out, in, H, W, C, ws, (int64_t)H * W * C);
}

void add_tiled(const Tensor& out, const Tensor& in, const Tensor& tile, int H, int W,
               int C, int th, int tw) {
    TiledAddParams p{};
    p.out = out.ptr;
    p.x = in.ptr;
    p.tile = tile.ptr;
    p.H = (uint32_t)H;
    p.W = (uint32_t)W;
    p.C = (uint32_t)C;
    p.th = (uint32_t)th;
    p.tw = (uint32_t)tw;
    vk::SpecList spec{(uint32_t)(in.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat("misc.add_tiled", spec, (int64_t)H * W * C, 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

void roi_align(const Tensor& out, const Tensor& feat, const Tensor& boxes, int H, int W,
               int C, int S) {
    RoiAlignParams p{};
    p.out = out.ptr;
    p.feat = feat.ptr;
    p.boxes = boxes.ptr;
    p.nboxes = (uint32_t)boxes.rows();
    p.H = (uint32_t)H;
    p.W = (uint32_t)W;
    p.C = (uint32_t)C;
    p.S = (uint32_t)S;
    vk::SpecList spec{(uint32_t)(feat.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat("misc.roi_align", spec, out.numel(), 256, &p,
                                   sizeof(p), &p.groups_per_row);
}

// ================
// Resample
// ================

static void resize_op(const char* entry, const Tensor& out, const Tensor& in) {
    NN_CHECK(out.ndim == 3 && in.ndim == 3,
               "%s expects [H, W, C] tensors (got %dD and %dD)", entry, out.ndim,
               in.ndim);
    NN_CHECK(out.shape[2] == in.shape[2], "%s: channel counts differ (%lld vs %lld)",
               entry, (long long)out.shape[2], (long long)in.shape[2]);
    ResizeParams p{};
    p.out = out.ptr;
    p.x = in.ptr;
    p.Ho = (uint32_t)out.shape[0];
    p.Wo = (uint32_t)out.shape[1];
    p.Hi = (uint32_t)in.shape[0];
    p.Wi = (uint32_t)in.shape[1];
    p.C = (uint32_t)out.shape[2];
    vk::SpecList spec{(uint32_t)(in.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat(entry, spec, out.numel(), 256, &p, sizeof(p),
                                   &p.groups_per_row);
}

void resize_bilinear(const Tensor& out, const Tensor& in) {
    resize_op("resample.resize_bilinear", out, in);
}
void upsample_nearest2x(const Tensor& out, const Tensor& in) {
    resize_op("resample.upsample_nearest2x", out, in);
}
void maxpool2x2(const Tensor& out, const Tensor& in) {
    resize_op("resample.maxpool2x2", out, in);
}

void resize_binarize(const Tensor& out_u8, const Tensor& logits, int64_t Ho, int64_t Wo,
                     float threshold) {
    NN_CHECK(out_u8.dtype == DType::U8, "resize_binarize writes a u8 tensor");
    MaskExportParams p{};
    p.out = out_u8.ptr;
    p.x = logits.ptr;
    p.Ho = (uint32_t)Ho;
    p.Wo = (uint32_t)Wo;
    p.Hi = (uint32_t)logits.dim(0);
    p.Wi = (uint32_t)logits.dim(1);
    p.threshold = threshold;
    vk::SpecList spec{(uint32_t)(logits.dtype == DType::F16), 0u};
    vk::Stream::get().dispatchFlat("resample.resize_binarize", spec, (Ho * Wo + 3) / 4,
                                   256, &p, sizeof(p), &p.groups_per_row);
}

}  // namespace nn
