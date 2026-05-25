#pragma once

#include "Tensor.h"

#include "IntersectTile.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
#include "SVHash.cuh"
// #include "Projection.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"
#include "Optimizer.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "Densify.cuh"
#include "BilagridUtils.cuh"
#include "Visualizer.cuh"


void set_data_3dgs(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
);

void set_camera_params(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
);

void forward_3dgs(
    std::string primitive,
    int sh_degree,
    bool packed,
    TorchTensorView out_rgb,
    TorchTensorView out_depth,
    TorchTensorView out_Ts
);

void backward_3dgs(
    TorchTensorView v_rgb,
    TorchTensorView v_depth,
    TorchTensorView v_Ts,
    TorchTensorView v_means,
    TorchTensorView v_quats,
    TorchTensorView v_scales,
    TorchTensorView v_opacities,
    TorchTensorView v_features_dc,
    TorchTensorView v_features_sh
);

void forward(std::string primitive, int sh_degree, bool packed);
