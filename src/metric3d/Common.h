#pragma once
// The names metric3d/ borrows from the inference layer (src/nn/), exactly as
// sam/Common.h and aliked/Common.h do. Include this at the top of a metric3d/
// translation unit instead of qualifying nn:: at every call site.
//
// Nothing model-specific belongs here -- that is what metric3d/Metric3D.h and
// model/Hparams.h are for.

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"
#include "nn/io/Onnx.h"

namespace nn { namespace vk {} }

namespace metric3d {

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

}  // namespace metric3d
