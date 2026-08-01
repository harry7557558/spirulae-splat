#pragma once
// The names sam/ borrows from the inference layer (src/nn/).
//
// Every module here is written against nn's Vulkan runtime, its exception type
// and its logging; qualifying all three at ~400 call sites would say nothing a
// reader does not already know from the directory. Include this once at the
// top of a sam/ translation unit or internal header instead. src/video/ opens
// the same way.
//
// Nothing model-specific belongs in here -- that is what sam/Sam.h is for.

#include "nn/core/Error.h"
#include "nn/core/Half.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"

// Opened empty so the alias below is well-formed without pulling vulkan.h into
// every sam/ header; the members are found through the alias as usual.
namespace nn { namespace vk {} }

namespace sam {

namespace vk = ::nn::vk;
using ::nn::Error;
using ::nn::fail;
using ::nn::now_ms;
using ::nn::log_level;
using ::nn::ScopedTimer;
using ::nn::parallel_for;
using ::nn::half_to_float;
using ::nn::float_to_half;

}  // namespace sam
