#pragma once
// Checkpoint loading: file -> a few device buffers -> named Tensors.
//
// Weights are packed into 256 MiB chunks, sub-allocated at 256-byte alignment.
// A SAM 3 checkpoint has ~1465 tensors; one allocation each would sit
// uncomfortably close to maxMemoryAllocationCount (4096 on most drivers) and
// buy nothing, since weights are never resized. One allocation for all of them
// does not work either: NVIDIA's Windows driver refuses a single device-local
// allocation of 1 GiB on an 8 GiB laptop part, with 7 GiB free -- so a 1.6 GiB
// f16 SAM 3 has to arrive in pieces regardless of what the heap size says.
//
// Device dtype policy (see AGENTS.md "Precision"):
//   * matmul / conv kernels          -> F16, halving VRAM at no measurable cost
//   * everything under 65536 elements-> F32; these are the learned tokens,
//                                       embeddings tables and reference points,
//                                       they cost nothing and some feed
//                                       coordinate math where fp16 shows
//   * RoPE frequency tables          -> F32, the kernel reads them as raw f32
//
// Quantized checkpoints (q4_0 / q4_1 / q8_0) are dequantized on the host during
// load. They exist to shrink the *file*; keeping the block codecs on the GPU
// would cost more in the inner loop than the bandwidth it saves at these matrix
// sizes, and it would put a second numeric path under every kernel.

#include "sam/Common.h"
#include "sam/model/Hparams.h"
#include "sam/model/Tokenizer.h"
#include "nn/Tensor.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace sam {

class WeightStore {
public:
    // Parses the header, uploads every tensor, and reads the embedded BPE
    // vocabulary. Throws nn::Error with a specific message on any mismatch.
    void load(const std::string& path);

    const Hparams&      hparams() const { return hparams_; }
    const BpeTokenizer& tokenizer() const { return tokenizer_; }
    bool                hasTokenizer() const { return tokenizer_.loaded(); }
    bool                visualOnly() const { return hparams_.visual_only != 0; }

    // Throws when the name is absent, listing near misses -- a typo in a weight
    // name is otherwise a null tensor that faults much later.
    nn::Tensor get(const std::string& name) const;
    nn::Tensor tryGet(const std::string& name) const;
    bool       has(const std::string& name) const { return tensors_.count(name) != 0; }

    // printf-style lookup for the per-layer names ("vit.blocks.%d.norm1.weight").
    nn::Tensor getf(const char* fmt, ...) const;
    nn::Tensor tryGetf(const char* fmt, ...) const;

    uint64_t deviceBytes() const { return device_bytes_; }
    size_t   count() const { return tensors_.size(); }

private:
    struct Record {
        std::string        name;
        std::vector<int64_t> shape;  // PyTorch order
        FileDType          dtype = FileDType::F32;
        uint64_t           file_offset = 0;
        uint64_t           file_bytes = 0;
        nn::DType          device_dtype = nn::DType::F16;
        uint32_t           chunk = 0;
        uint64_t           device_offset = 0;   // within its chunk
    };

    Hparams      hparams_;
    BpeTokenizer tokenizer_;
    std::unordered_map<std::string, nn::Tensor> tensors_;
    uint64_t     device_bytes_ = 0;
    std::vector<nn::DevicePtr> chunks_;
};

// ---------------------------------------------------------------------------
// Weight post-processing helpers, used by the modules at build time.
// ---------------------------------------------------------------------------

// ConvTranspose2d(k=2, s=2) weights arrive as [Cin, Cout, 2, 2]. The GEMM path
// wants [Cout*4, Cin] so the four kernel taps become four output-channel
// groups (see nn::conv_transpose2x2). Returns a new pool-backed tensor.
nn::Tensor repack_conv_transpose2x2(const nn::Tensor& w, int64_t cin, int64_t cout,
                                    uint32_t pool_sub);

// Reorder a [Cout, Cin, k, k] kernel into [Cout, k*k*Cin] so it can multiply a
// channel-last patch. Used for the geometry encoder's box pooling projection,
// which is a Conv2d(D, D, 7) applied to a 7x7 ROI-Align output -- i.e. a single
// matmul once the patch is laid out (ky, kx, ci) instead of (ci, ky, kx).
nn::Tensor repack_conv_to_hwc(const nn::Tensor& w, int64_t cout, int64_t cin, int64_t k,
                              uint32_t pool_sub);

}  // namespace sam
