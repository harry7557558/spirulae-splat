#include <cuda_runtime.h>

#include <torch/types.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor compute_sh_forward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs0,   // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
);


std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
compute_sh_backward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs,  // [..., 3]
    torch::Tensor &colors,  // [..., 3]
    torch::Tensor &v_colors  // [..., 3]
);
