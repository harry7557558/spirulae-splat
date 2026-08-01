// The additive side of the SPIR-V registry -- see nn/vk/EmbeddedSpirv.h.
//
// The blob arrays themselves are generated; this is the small piece of hand
// code that lets several generated translation units share one table.

#include "nn/vk/EmbeddedSpirv.h"

#include <vector>

namespace nn {
namespace vk {

namespace {

// Function-local so registration from a static initializer in another TU
// cannot run before the vector is constructed.
std::vector<EmbeddedModule>& registry() {
    static std::vector<EmbeddedModule> r;
    return r;
}

}  // namespace

void register_modules(const EmbeddedModule* mods, size_t count) {
    registry().insert(registry().end(), mods, mods + count);
}

const EmbeddedModule* embedded_modules(size_t* count) {
    *count = registry().size();
    return registry().data();
}

}  // namespace vk
}  // namespace nn
