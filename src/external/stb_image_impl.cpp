// Single translation unit that materializes the stb_image implementation.
// All other TUs `#include "external/stb_image.h"` to get the declarations only.

#define STB_IMAGE_IMPLEMENTATION

// Compile-time slim-down: drop formats nothing here reads. PNG / JPEG / BMP /
// TGA cover everything the dataparser sees; PNM is kept because the SfM module
// accepts .ppm/.pgm inputs and masks (src/sfm/core/Image.h, Mask.h).
#define STBI_NO_PSD
#define STBI_NO_HDR
#define STBI_NO_PIC
#define STBI_NO_GIF

// Silence two warnings stb_image is noisy about under -Wsign-compare /
// -Wunused-but-set-variable on recent GCCs.
#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsign-compare"
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "external/stb_image.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
