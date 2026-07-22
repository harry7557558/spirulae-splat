// Hand-written companions to TrainingStubs.gen.cpp: throwing stubs for the
// symbols the generator cannot auto-stub (the meshing::OccupancyEvaluator
// class methods). Same contract — they satisfy the linker for the portable
// engine layer under the Vulkan backend; calling one throws. Remove entries
// as the corresponding modules get ported.

#include <mesh/Meshing.h>
#include <core/Tensor.h>

#include <cstdint>
#include <stdexcept>
#include <string>

namespace {
[[noreturn]] void _vk_stub(const char* name) {
    throw std::runtime_error(std::string("Vulkan backend: ") + name +
                             " is not implemented yet (training phase)");
}
}  // namespace

// Meshing occupancy evaluator (Meshing.cu). The constructor throws, so no
// other method can ever run on a live object; the destructor must stay
// non-throwing and has nothing to release.
namespace meshing {

OccupancyEvaluator::OccupancyEvaluator(
    const float*, const float*, const float*, const float*, const float*,
    int, const float*, int, const CameraParams&, const MeshingConfig&
) : impl_(nullptr) {
    _vk_stub("meshing::OccupancyEvaluator");
}

OccupancyEvaluator::~OccupancyEvaluator() {}

void OccupancyEvaluator::generate_point_cloud(std::vector<float>&) {
    _vk_stub("meshing::OccupancyEvaluator::generate_point_cloud");
}
void OccupancyEvaluator::evaluate(const float*, int, float*) {
    _vk_stub("meshing::OccupancyEvaluator::evaluate");
}
void OccupancyEvaluator::bisect_edges(const float*, const int32_t*,
                                      const int32_t*, const float*,
                                      const float*, int, float*) {
    _vk_stub("meshing::OccupancyEvaluator::bisect_edges");
}
void OccupancyEvaluator::colorize(const float*, int, float*) {
    _vk_stub("meshing::OccupancyEvaluator::colorize");
}
void OccupancyEvaluator::view_texel_density(const float*, int, float*) {
    _vk_stub("meshing::OccupancyEvaluator::view_texel_density");
}
bool OccupancyEvaluator::has_render_cameras() const {
    _vk_stub("meshing::OccupancyEvaluator::has_render_cameras");
}
void OccupancyEvaluator::cull_unseen_vertices(const float*, int, const int*,
                                              int, unsigned char*) {
    _vk_stub("meshing::OccupancyEvaluator::cull_unseen_vertices");
}
bool OccupancyEvaluator::debug_render_moments(int, std::vector<float>&, int&,
                                              int&) {
    _vk_stub("meshing::OccupancyEvaluator::debug_render_moments");
}

}  // namespace meshing
