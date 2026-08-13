// Single translation unit that materializes the stb_image_write
// implementation (PNG encoding for the mesh texture export + the viewer's
// JPEG blit). All other TUs `#include "external/stb_image_write.h"` for the
// declarations only.

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsign-compare"
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
// Apple clang deprecates sprintf; the HDR writer uses it and stb is vendored.
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
