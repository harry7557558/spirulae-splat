#pragma once


#include "types.cuh"

enum class RawPPISPRegLossIndex {
    SumExposure,
    SumVignettingCrSquared,
    SumVignettingAlpha0Relu,
    SumVignettingAlpha1Relu,
    SumVignettingAlpha2Relu,
    SumVignettingCxChannelVariance,
    SumVignettingCyChannelVariance,
    SumVignettingAlpha0ChannelVariance,
    SumVignettingAlpha1ChannelVariance,
    SumVignettingAlpha2ChannelVariance,
    SumColorBx,
    SumColorBy,
    SumColorRx,
    SumColorRy,
    SumColorGx,
    SumColorGy,
    SumColorNx,
    SumColorNy,
    SumCRFToeChannelVariance,
    SumCRFShoulderChannelVariance,
    SumCRFGammaChannelVariance,
    SumCRFCenterChannelVariance,
    length
};

enum class RawPPISPRegLossIndexRQS {
    SumExposure,
    SumVignettingCrSquared,
    SumVignettingAlpha0Relu,
    SumVignettingAlpha1Relu,
    SumVignettingAlpha2Relu,
    SumVignettingCxChannelVariance,
    SumVignettingCyChannelVariance,
    SumVignettingAlpha0ChannelVariance,
    SumVignettingAlpha1ChannelVariance,
    SumVignettingAlpha2ChannelVariance,
    SumColorBx,
    SumColorBy,
    SumColorRx,
    SumColorRy,
    SumColorGx,
    SumColorGy,
    SumColorNx,
    SumColorNy,
    SumCRFG0ChannelVariance,
    SumCRFG1ChannelVariance,
    SumCRFX0ChannelVariance,
    SumCRFY0ChannelVariance,
    SumCRFGcChannelVariance,
    length
};

enum class PPISPRegLossIndex {
    ExposureMean,
    VignettingCenter,
    VignettingNonPositivity,
    VignettingChannelVariance,
    ColorMean,
    CRFChannelVariance,
    length
};


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void rendered_depth_to_expected_depth_forward(
    TorchTensorView depth,  // [B, H, W, 1]
    TorchTensorView transmittance,  // [B, H, W, 1]
    TorchTensorView out_depth  // [B, H, W, 1]
);


void rendered_depth_to_expected_depth_backward(
    TorchTensorView depth,  // [B, H, W, 1]
    TorchTensorView transmittance,  // [B, H, W, 1]
    TorchTensorView v_out_depth,  // [B, H, W, 1]
    TorchTensorView v_depth,  // [B, H, W, 1]
    TorchTensorView v_transmittance  // [B, H, W, 1]
);


void ray_depth_to_linear_depth_forward(
    std::string camera_model,
    TorchTensorView intrins,  // [B, 4]
    TorchTensorView dist_coeffs,  // [B, 10]
    TorchTensorView depths,  // [B, H, W, 1]
    TorchTensorView out_depths  // [B, H, W, 1]
);


void ray_depth_to_linear_depth_backward(
    std::string camera_model,
    TorchTensorView intrins,  // [B, 4]
    TorchTensorView dist_coeffs,  // [B, 10]
    TorchTensorView v_out_depths,  // [B, H, W, 1]
    TorchTensorView v_in_depths  // [B, H, W, 1]
);
