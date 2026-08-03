#pragma once
// The names aliked/ borrows from the inference layer (src/nn/), exactly as
// sam/Common.h does for src/sam/. Include this at the top of an aliked/
// translation unit instead of qualifying nn:: at every call site.
//
// Nothing model-specific belongs here -- that is what aliked/Aliked.h is for.

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"

namespace nn { namespace vk {} }

namespace aliked {

namespace vk = ::nn::vk;
using ::nn::Error;
using ::nn::fail;
using ::nn::log_level;
using ::nn::now_ms;
using ::nn::ScopedTimer;
using ::nn::parallel_for;

}  // namespace aliked
