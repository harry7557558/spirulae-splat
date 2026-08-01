#pragma once
// The names video/ borrows from the inference layer (src/nn/) -- the same
// arrangement, and for the same reason, as src/sam/Common.h.

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"

namespace nn { namespace vk {} }

namespace video {

namespace vk = ::nn::vk;
using ::nn::Error;
using ::nn::fail;
using ::nn::now_ms;
using ::nn::log_level;
using ::nn::ScopedTimer;
using ::nn::parallel_for;
using ::nn::half_to_float;
using ::nn::float_to_half;

}  // namespace video
