#pragma once

#include <ATen/Tensor.h>

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void inplace_index_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


void inplace_scatter_add_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


void inplace_scatter_max_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
);


void mcmc_add_noise_tensor(
    std::string primitive,
    float scaler, float min_opacity,
    at::Tensor &means,
    at::Tensor &scales,
    at::Tensor &quats,
    at::Tensor &opacs
);


std::tuple<at::Tensor, at::Tensor>
compute_relocation_tensor(
    // inputs
    at::Tensor opacities, // [N]
    at::Tensor scales,    // [N, 3]
    at::Tensor ratios,    // [N]
    at::Tensor binoms,    // [n_max, n_max]
    const int n_max
);


std::tuple<at::Tensor, at::Tensor>
long_axis_split_tensor(
    std::string primitive,
    at::Tensor &scales,
    at::Tensor &quats
);
