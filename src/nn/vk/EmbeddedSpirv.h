#pragma once
// The SPIR-V blob registry, produced at build time by spirv_tool's `embed
// --nn` mode (cmake/SsNn.cmake).
//
// One entry per .slang module; `name` is the file stem, which is the prefix of
// every kernel key the pipeline cache resolves ("<stem>.<entry>").
//
// The registry is additive: each library contributes its own generated
// translation unit, which registers its table before main(). That is what lets
// src/nn/, src/sam/ and src/video/ be three link units over one pipeline
// cache -- and what lets a future subsystem add shaders without touching this
// one.

#include <cstddef>
#include <cstdint>

namespace nn {
namespace vk {

struct EmbeddedModule {
    const char*     name;
    const uint32_t* code;
    size_t          words;
};

void register_modules(const EmbeddedModule* mods, size_t count);

// Every module registered so far, flattened. Valid for the process lifetime;
// the pipeline cache rescans only the tail, so registering late is fine.
const EmbeddedModule* embedded_modules(size_t* count);

}  // namespace vk
}  // namespace nn

// The two halves of a library's registration.
//
// Registration is an explicit call rather than a static initializer for one
// unglamorous reason: these libraries are static archives, and a member object
// that nothing references is not linked at all -- a static initializer in the
// generated blob TU would silently never run, and the kernels would come back
// "no shader module". Calling it from code the library's own consumers already
// pull in is what makes it certain.
//
// Put NN_ENSURE_EMBEDDED_MODULES at the library's entry point (whatever runs
// before its first dispatch); it registers once and costs a load after that.
#define NN_DECLARE_EMBEDDED_MODULES(tag)                                       \
    namespace nn {                                                             \
    namespace vk {                                                             \
    extern const EmbeddedModule kEmbeddedModules_##tag[];                      \
    extern const size_t kEmbeddedModuleCount_##tag;                            \
    /* inline, so the guard is one object across every TU that calls it */     \
    inline void ensure_modules_##tag() {                                       \
        static const bool once =                                               \
            (register_modules(kEmbeddedModules_##tag,                          \
                              kEmbeddedModuleCount_##tag),                     \
             true);                                                            \
        (void)once;                                                            \
    }                                                                          \
    }                                                                          \
    }

#define NN_ENSURE_EMBEDDED_MODULES(tag) ::nn::vk::ensure_modules_##tag()
