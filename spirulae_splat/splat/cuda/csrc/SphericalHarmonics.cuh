#include <cuda_runtime.h>

#include <ATen/Tensor.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor compute_sh_forward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    at::Tensor &viewdirs,  // [..., 3]
    at::Tensor &coeffs0,   // [..., 3]
    at::Tensor &coeffs   // [..., K, 3]
);


std::tuple<at::Tensor, at::Tensor, at::Tensor>
compute_sh_backward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    at::Tensor &viewdirs,  // [..., 3]
    at::Tensor &coeffs,  // [..., 3]
    at::Tensor &colors,  // [..., 3]
    at::Tensor &v_colors  // [..., 3]
);
