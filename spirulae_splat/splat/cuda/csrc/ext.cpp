#include "config.h"
#include "bindings.h"
#include <torch/extension.h>

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {

    // constants
    m.attr("BLOCK_WIDTH") = py::int_(BLOCK_WIDTH);
    m.attr("MAX_SORTED_SPLATS") = py::int_(MAX_SORTED_SPLATS);
    m.attr("SORTED_INDEX_INF") = py::int_(SORTED_INDEX_INF);

    // splat rasterization
    m.def("rasterize_simple_forward", &rasterize_simple_forward_tensor);
    m.def("rasterize_simple_backward", &rasterize_simple_backward_tensor);
    m.def("rasterize_depth_forward", &rasterize_depth_forward_tensor);
    m.def("rasterize_depth_backward", &rasterize_depth_backward_tensor);
    m.def("rasterize_forward", &rasterize_forward_tensor);
    m.def("rasterize_backward", &rasterize_backward_tensor);
    m.def("rasterize_simplified_forward", &rasterize_simplified_forward_tensor);
    m.def("rasterize_simplified_backward", &rasterize_simplified_backward_tensor);

    // sorted splat rasterization
    m.def("rasterize_indices", &rasterize_indices_tensor);
    m.def("sort_per_pixel", &sort_per_pixel_tensor);
    m.def("rasterize_simple_sorted_forward", &rasterize_simple_sorted_forward_tensor);
    m.def("rasterize_simple_sorted_backward", &rasterize_simple_sorted_backward_tensor);
    m.def("rasterize_depth_sorted_forward", &rasterize_depth_sorted_forward_tensor);
    m.def("rasterize_depth_sorted_backward", &rasterize_depth_sorted_backward_tensor);
    m.def("rasterize_sorted_forward", &rasterize_sorted_forward_tensor);
    m.def("rasterize_sorted_backward", &rasterize_sorted_backward_tensor);
    m.def("rasterize_simplified_sorted_forward", &rasterize_simplified_sorted_forward_tensor);
    m.def("rasterize_simplified_sorted_backward", &rasterize_simplified_sorted_backward_tensor);

    // background
    m.def("render_background_sh_forward", &render_background_sh_forward_tensor);
    m.def("render_background_sh_backward", &render_background_sh_backward_tensor);
    m.def("render_undistortion_map", &render_undistortion_map_tensor);

    // projection
    m.def("map_gaussian_to_intersects", &map_gaussian_to_intersects_tensor);
    m.def("get_tile_bin_edges", &get_tile_bin_edges_tensor);
    m.def("project_gaussians_forward", &project_gaussians_forward_tensor);
    m.def("project_gaussians_backward", &project_gaussians_backward_tensor);
    m.def("compute_sh_forward", &compute_sh_forward_tensor);
    m.def("compute_sh_backward", &compute_sh_backward_tensor);

    // misc
    m.def("compute_relocation", &compute_relocation_tensor);
    m.def("compute_relocation_split", &compute_relocation_split_tensor);
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward_tensor);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward_tensor);
    m.def("blend_background_forward", &blend_background_forward_tensor);
    m.def("blend_background_backward", &blend_background_backward_tensor);
    
}
