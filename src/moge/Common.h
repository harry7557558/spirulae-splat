#pragma once
// The names moge/ borrows from the inference layer (src/nn/), exactly as
// metric3d/Common.h does. Include this at the top of a moge/ translation unit
// instead of qualifying nn:: at every call site.
//
// Nothing model-specific belongs here -- that is what moge/Moge.h and
// model/Hparams.h are for.

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"
#include "nn/io/Onnx.h"

namespace nn { namespace vk {} }

namespace moge {

namespace vk = ::nn::vk;
using ::nn::Error;
using ::nn::fail;
using ::nn::log_level;
using ::nn::now_ms;
using ::nn::OnnxFile;
using ::nn::OnnxNode;
using ::nn::OnnxTensor;
using ::nn::parallel_for;
using ::nn::ScopedTimer;

}  // namespace moge
