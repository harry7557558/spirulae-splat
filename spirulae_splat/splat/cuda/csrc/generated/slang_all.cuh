#pragma once

#include "slang.cuh"

struct VignettingChannelParams_0
{
    float cx_0;
    float cy_0;
    float alpha0_0;
    float alpha1_0;
    float alpha2_0;
};

inline __device__ VignettingChannelParams_0 VignettingChannelParams_x24_syn_dzero_0()
{
    VignettingChannelParams_0 result_0;
    (&result_0)->cx_0 = 0.0f;
    (&result_0)->cy_0 = 0.0f;
    (&result_0)->alpha0_0 = 0.0f;
    (&result_0)->alpha1_0 = 0.0f;
    (&result_0)->alpha2_0 = 0.0f;
    return result_0;
}

struct ColorPPISPParams_0
{
    float2  b_0;
    float2  r_0;
    float2  g_0;
    float2  n_0;
};

inline __device__ ColorPPISPParams_0 ColorPPISPParams_x24_syn_dzero_0()
{
    ColorPPISPParams_0 result_1;
    float2  _S1 = make_float2 (0.0f);
    (&result_1)->b_0 = _S1;
    (&result_1)->r_0 = _S1;
    (&result_1)->g_0 = _S1;
    (&result_1)->n_0 = _S1;
    return result_1;
}

struct RQSCRFPPISPChannelParams_0
{
    float g0_0;
    float g1_0;
    float x0_0;
    float y0_0;
    float gc_0;
};

inline __device__ RQSCRFPPISPChannelParams_0 RQSCRFPPISPChannelParams_x24_syn_dzero_0()
{
    RQSCRFPPISPChannelParams_0 result_2;
    (&result_2)->g0_0 = 0.0f;
    (&result_2)->g1_0 = 0.0f;
    (&result_2)->x0_0 = 0.0f;
    (&result_2)->y0_0 = 0.0f;
    (&result_2)->gc_0 = 0.0f;
    return result_2;
}

struct PPISPParamsRQS_0
{
    float exposure_0;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_0;
    ColorPPISPParams_0 color_params_0;
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  crf_params_0;
};

inline __device__ PPISPParamsRQS_0 PPISPParamsRQS_x24_syn_dzero_0()
{
    PPISPParamsRQS_0 result_3;
    (&result_3)->exposure_0 = 0.0f;
    VignettingChannelParams_0 _S2 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_3)->vignette_params_0[int(0)] = _S2;
    (&result_3)->vignette_params_0[int(1)] = _S2;
    (&result_3)->vignette_params_0[int(2)] = _S2;
    (&result_3)->color_params_0 = ColorPPISPParams_x24_syn_dzero_0();
    RQSCRFPPISPChannelParams_0 _S3 = RQSCRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_3)->crf_params_0[int(0)] = _S3;
    (&result_3)->crf_params_0[int(1)] = _S3;
    (&result_3)->crf_params_0[int(2)] = _S3;
    return result_3;
}

struct CRFPPISPChannelParams_0
{
    float toe_0;
    float shoulder_0;
    float gamma_0;
    float center_0;
};

inline __device__ CRFPPISPChannelParams_0 CRFPPISPChannelParams_x24_syn_dzero_0()
{
    CRFPPISPChannelParams_0 result_4;
    (&result_4)->toe_0 = 0.0f;
    (&result_4)->shoulder_0 = 0.0f;
    (&result_4)->gamma_0 = 0.0f;
    (&result_4)->center_0 = 0.0f;
    return result_4;
}

struct PPISPParams_0
{
    float exposure_1;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_1;
    ColorPPISPParams_0 color_params_1;
    FixedArray<CRFPPISPChannelParams_0, 3>  crf_params_1;
};

inline __device__ PPISPParams_0 PPISPParams_x24_syn_dzero_0()
{
    PPISPParams_0 result_5;
    (&result_5)->exposure_1 = 0.0f;
    VignettingChannelParams_0 _S4 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_5)->vignette_params_1[int(0)] = _S4;
    (&result_5)->vignette_params_1[int(1)] = _S4;
    (&result_5)->vignette_params_1[int(2)] = _S4;
    (&result_5)->color_params_1 = ColorPPISPParams_x24_syn_dzero_0();
    CRFPPISPChannelParams_0 _S5 = CRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_5)->crf_params_1[int(0)] = _S5;
    (&result_5)->crf_params_1[int(1)] = _S5;
    (&result_5)->crf_params_1[int(2)] = _S5;
    return result_5;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S6 = (F32_exp(((*dpx_0).primal_0))) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S6;
    return;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S7 = *dpx_1;
    float _S8;
    if(((*dpx_1).primal_0) > ((*dpy_0).primal_0))
    {
        _S8 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_0).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S7.primal_0;
    dpx_1->differential_0 = _S8;
    DiffPair_float_0 _S9 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S7.primal_0))
    {
        _S8 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_1).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S9.primal_0;
    dpy_0->differential_0 = _S8;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S10 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_2).primal_0)))))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S10;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_3)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_3;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_3).primal_0.x * dOut_3;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_3;
    *&((&y_d_result_0)->y) = (*dpx_3).primal_0.y * dOut_3;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_3;
    *&((&y_d_result_0)->z) = (*dpx_3).primal_0.z * dOut_3;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_0, float3  y_0)
{
    int i_0 = int(0);
    float result_6 = 0.0f;
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_7 = result_6 + _slang_vector_get_element(x_0, i_0) * _slang_vector_get_element(y_0, i_0);
        i_0 = i_0 + int(1);
        result_6 = result_7;
    }
    return result_6;
}

inline __device__ float dot_1(float4  x_1, float4  y_1)
{
    int i_1 = int(0);
    float result_8 = 0.0f;
    for(;;)
    {
        if(i_1 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_9 = result_8 + _slang_vector_get_element(x_1, i_1) * _slang_vector_get_element(y_1, i_1);
        i_1 = i_1 + int(1);
        result_8 = result_9;
    }
    return result_8;
}

inline __device__ float dot_2(float2  x_2, float2  y_2)
{
    int i_2 = int(0);
    float result_10 = 0.0f;
    for(;;)
    {
        if(i_2 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_11 = result_10 + _slang_vector_get_element(x_2, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_10 = result_11;
    }
    return result_10;
}

inline __device__ float length_0(float4  x_3)
{
    return (F32_sqrt((dot_1(x_3, x_3))));
}

inline __device__ float length_1(float2  x_4)
{
    return (F32_sqrt((dot_2(x_4, x_4))));
}

inline __device__ float length_2(float3  x_5)
{
    return (F32_sqrt((dot_0(x_5, x_5))));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_4, float dOut_4)
{
    float _S11 = 1.0f / (*dpx_4).primal_0 * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S11;
    return;
}

inline __device__ float3  exp_0(float3  x_6)
{
    float3  result_12;
    int i_3 = int(0);
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_12, i_3) = (F32_exp((_slang_vector_get_element(x_6, i_3))));
        i_3 = i_3 + int(1);
    }
    return result_12;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float3  dOut_5)
{
    float3  _S12 = exp_0((*dpx_5).primal_0) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S12;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ float2  exp_1(float2  x_7)
{
    float2  result_13;
    int i_4 = int(0);
    for(;;)
    {
        if(i_4 < int(2))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_13, i_4) = (F32_exp((_slang_vector_get_element(x_7, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_13;
}

inline __device__ void _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_6, float2  dOut_6)
{
    float2  _S13 = exp_1((*dpx_6).primal_0) * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S13;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpy_2, float dOut_7)
{
    DiffPair_float_0 _S14 = *dpx_7;
    float _S15;
    if(((*dpx_7).primal_0) < ((*dpy_2).primal_0))
    {
        _S15 = dOut_7;
    }
    else
    {
        if(((*dpx_7).primal_0) > ((*dpy_2).primal_0))
        {
            _S15 = 0.0f;
        }
        else
        {
            _S15 = 0.5f * dOut_7;
        }
    }
    dpx_7->primal_0 = _S14.primal_0;
    dpx_7->differential_0 = _S15;
    DiffPair_float_0 _S16 = *dpy_2;
    if(((*dpy_2).primal_0) < (_S14.primal_0))
    {
        _S15 = dOut_7;
    }
    else
    {
        if(((*dpy_2).primal_0) > ((*dpx_7).primal_0))
        {
            _S15 = 0.0f;
        }
        else
        {
            _S15 = 0.5f * dOut_7;
        }
    }
    dpy_2->primal_0 = _S16.primal_0;
    dpy_2->differential_0 = _S15;
    return;
}

inline __device__ void per_splat_losses(bool is_3dgs_0, float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S17)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    if(is_3dgs_0)
    {
        float3  _S18 = exp_0(scales_0);
        float _S19 = _S18.x;
        float _S20 = _S18.y;
        float _S21 = _S18.z;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S19 + _S20 + _S21) / 3.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S19), (_S20)))), (_S21))) / (F32_min(((F32_min((_S19), (_S20)))), (_S21)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float3  _S22 = exp_0(make_float3 (2.0f) * scales_0);
        float x_8 = _S22.x;
        float y_3 = _S22.y;
        float z_0 = _S22.z;
        float s_0 = x_8 + y_3 + z_0;
        float s1_0 = (F32_max(((F32_max((x_8), (y_3)))), (z_0))) / s_0;
        float s3_0 = (F32_min(((F32_min((x_8), (y_3)))), (z_0))) / s_0;
        float s2_0 = 1.0f - s1_0 - s3_0;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    }
    else
    {
        float2  _S23 = float2 {scales_0.x, scales_0.y};
        float2  _S24 = exp_1(_S23);
        float _S25 = _S24.x;
        float _S26 = _S24.y;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S25 + _S26) / 2.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max((_S25), (_S26))) / (F32_min((_S25), (_S26)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float2  _S27 = exp_1(make_float2 (2.0f) * _S23);
        float x_9 = _S27.x;
        float y_4 = _S27.y;
        float s_1 = x_9 + y_4;
        float s1_1 = (F32_max((x_9), (y_4))) / s_1;
        float s2_1 = (F32_min((x_9), (y_4))) / s_1;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_1 * (F32_log((s1_1))) - s2_1 * (F32_log((s2_1)))))) - 0.99998998641967773f)))), (0.0f)));
    }
    *_S17 = losses_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ float s_primal_ctx_exp_0(float _S28)
{
    return (F32_exp((_S28)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S29)
{
    return exp_0(_S29);
}

inline __device__ float s_primal_ctx_log_0(float _S30)
{
    return (F32_log((_S30)));
}

inline __device__ float2  s_primal_ctx_exp_2(float2  _S31)
{
    return exp_1(_S31);
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S32, float _S33)
{
    _d_log_0(_S32, _S33);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S34, float _S35)
{
    _d_exp_0(_S34, _S35);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S36, float3  _S37)
{
    _d_exp_vector_0(_S36, _S37);
    return;
}

inline __device__ void s_bwd_prop_exp_2(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S38, float2  _S39)
{
    _d_exp_vector_1(_S38, _S39);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S40, float _S41)
{
    _d_sqrt_0(_S40, _S41);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_8, float _s_dOut_0)
{
    float _S42 = (*dpx_8).primal_0.x;
    float _S43 = (*dpx_8).primal_0.y;
    float _S44 = (*dpx_8).primal_0.z;
    float _S45 = (*dpx_8).primal_0.w;
    DiffPair_float_0 _S46;
    (&_S46)->primal_0 = _S42 * _S42 + _S43 * _S43 + _S44 * _S44 + _S45 * _S45;
    (&_S46)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S46, _s_dOut_0);
    float _S47 = (*dpx_8).primal_0.w * _S46.differential_0;
    float _S48 = _S47 + _S47;
    float _S49 = (*dpx_8).primal_0.z * _S46.differential_0;
    float _S50 = _S49 + _S49;
    float _S51 = (*dpx_8).primal_0.y * _S46.differential_0;
    float _S52 = _S51 + _S51;
    float _S53 = (*dpx_8).primal_0.x * _S46.differential_0;
    float _S54 = _S53 + _S53;
    float4  _S55 = make_float4 (0.0f);
    *&((&_S55)->w) = _S48;
    *&((&_S55)->z) = _S50;
    *&((&_S55)->y) = _S52;
    *&((&_S55)->x) = _S54;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S55;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S56, float _S57)
{
    s_bwd_prop_length_impl_0(_S56, _S57);
    return;
}

inline __device__ void s_bwd_prop_per_splat_losses_0(bool is_3dgs_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1, FixedArray<float, 5>  * _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S58 = *dpscales_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S59 = *dpquat_0;
    float3  _S60 = make_float3 (0.0f);
    float2  _S61 = make_float2 (0.0f);
    float _S62 = - (*dpopacity_0).primal_0;
    float _S63 = 1.0f + s_primal_ctx_exp_0(_S62);
    float _S64 = _S63 * _S63;
    float _S65 = length_0((*dpquat_0).primal_0);
    float _S66;
    float _S67;
    float _S68;
    float _S69;
    float _S70;
    float _S71;
    float _S72;
    float _S73;
    float _S74;
    float _S75;
    float _S76;
    float _S77;
    float _S78;
    float _S79;
    float _S80;
    float _S81;
    float _S82;
    float _S83;
    float _S84;
    float _S85;
    float _S86;
    float _S87;
    float _S88;
    float _S89;
    float _S90;
    float _S91;
    float _S92;
    float _S93;
    float _S94;
    float _S95;
    float _S96;
    float _S97;
    float _S98;
    float _S99;
    float _S100;
    float _S101;
    float _S102;
    float _S103;
    float _S104;
    float _S105;
    float _S106;
    float _S107;
    float _S108;
    float _S109;
    float _S110;
    float _S111;
    float _S112;
    float _S113;
    float3  _S114;
    float2  _S115;
    float2  _S116;
    if(is_3dgs_1)
    {
        float3  _S117 = s_primal_ctx_exp_1(_S58.primal_0);
        float _S118 = _S117.x;
        float _S119 = _S117.y;
        float _S120 = _S117.z;
        float _S121 = (F32_max((_S118), (_S119)));
        float _S122 = (F32_max((_S121), (_S120)));
        float _S123 = (F32_min((_S118), (_S119)));
        float _S124 = (F32_min((_S123), (_S120)));
        float _S125 = _S122 / _S124;
        float _S126 = _S124 * _S124;
        float3  _S127 = make_float3 (2.0f) * _S58.primal_0;
        float3  _S128 = s_primal_ctx_exp_1(_S127);
        float x_10 = _S128.x;
        float y_5 = _S128.y;
        float z_1 = _S128.z;
        float s_2 = x_10 + y_5 + z_1;
        float _S129 = (F32_max((x_10), (y_5)));
        float _S130 = (F32_max((_S129), (z_1)));
        float s1_2 = _S130 / s_2;
        float _S131 = s_2 * s_2;
        float _S132 = (F32_min((x_10), (y_5)));
        float _S133 = (F32_min((_S132), (z_1)));
        float s3_1 = _S133 / s_2;
        float s2_2 = 1.0f - s1_2 - s3_1;
        float _S134 = - s1_2;
        float _S135 = s_primal_ctx_log_0(s1_2);
        float _S136 = s_primal_ctx_log_0(s2_2);
        float _S137 = s_primal_ctx_log_0(s3_1);
        float _S138 = _S134 * _S135 - s2_2 * _S136 - s3_1 * _S137;
        float _S139 = s_primal_ctx_exp_0(_S138) - 0.99998998641967773f;
        _S66 = - s_primal_ctx_log_0(_S139);
        _S67 = _S139;
        _S68 = _S138;
        _S69 = s3_1;
        _S70 = _S137;
        _S71 = s2_2;
        _S72 = _S136;
        _S73 = _S134;
        _S74 = _S135;
        _S75 = s1_2;
        _S76 = _S131;
        _S77 = _S133;
        _S78 = s_2;
        _S79 = _S132;
        _S80 = z_1;
        _S81 = x_10;
        _S82 = y_5;
        _S83 = _S130;
        _S84 = _S129;
        _S114 = _S127;
        _S85 = _S125;
        _S86 = _S126;
        _S87 = _S122;
        _S88 = _S124;
        _S89 = _S123;
        _S90 = _S120;
        _S91 = _S118;
        _S92 = _S119;
        _S93 = _S121;
        _S94 = 0.0f;
        _S95 = 0.0f;
        _S96 = 0.0f;
        _S97 = 0.0f;
        _S98 = 0.0f;
        _S99 = 0.0f;
        _S100 = 0.0f;
        _S101 = 0.0f;
        _S102 = 0.0f;
        _S103 = 0.0f;
        _S104 = 0.0f;
        _S105 = 0.0f;
        _S106 = 0.0f;
        _S107 = 0.0f;
        _S115 = _S61;
        _S108 = 0.0f;
        _S109 = 0.0f;
        _S110 = 0.0f;
        _S111 = 0.0f;
        _S112 = 0.0f;
        _S113 = 0.0f;
        _S116 = _S61;
    }
    else
    {
        float2  _S140 = float2 {_S58.primal_0.x, _S58.primal_0.y};
        float2  _S141 = s_primal_ctx_exp_2(_S140);
        float _S142 = _S141.x;
        float _S143 = _S141.y;
        float _S144 = (F32_max((_S142), (_S143)));
        float _S145 = (F32_min((_S142), (_S143)));
        float _S146 = _S144 / _S145;
        float _S147 = _S145 * _S145;
        float2  _S148 = make_float2 (2.0f) * _S140;
        float2  _S149 = s_primal_ctx_exp_2(_S148);
        float x_11 = _S149.x;
        float y_6 = _S149.y;
        float s_3 = x_11 + y_6;
        float _S150 = (F32_max((x_11), (y_6)));
        float s1_3 = _S150 / s_3;
        float _S151 = s_3 * s_3;
        float _S152 = (F32_min((x_11), (y_6)));
        float s2_3 = _S152 / s_3;
        float _S153 = - s1_3;
        float _S154 = s_primal_ctx_log_0(s1_3);
        float _S155 = s_primal_ctx_log_0(s2_3);
        float _S156 = _S153 * _S154 - s2_3 * _S155;
        float _S157 = s_primal_ctx_exp_0(_S156) - 0.99998998641967773f;
        float _S158 = - s_primal_ctx_log_0(_S157);
        _S66 = 0.0f;
        _S67 = 0.0f;
        _S68 = 0.0f;
        _S69 = 0.0f;
        _S70 = 0.0f;
        _S71 = 0.0f;
        _S72 = 0.0f;
        _S73 = 0.0f;
        _S74 = 0.0f;
        _S75 = 0.0f;
        _S76 = 0.0f;
        _S77 = 0.0f;
        _S78 = 0.0f;
        _S79 = 0.0f;
        _S80 = 0.0f;
        _S81 = 0.0f;
        _S82 = 0.0f;
        _S83 = 0.0f;
        _S84 = 0.0f;
        _S114 = _S60;
        _S85 = 0.0f;
        _S86 = 0.0f;
        _S87 = 0.0f;
        _S88 = 0.0f;
        _S89 = 0.0f;
        _S90 = 0.0f;
        _S91 = 0.0f;
        _S92 = 0.0f;
        _S93 = 0.0f;
        _S94 = _S158;
        _S95 = _S157;
        _S96 = _S156;
        _S97 = s2_3;
        _S98 = _S155;
        _S99 = _S153;
        _S100 = _S154;
        _S101 = s1_3;
        _S102 = _S151;
        _S103 = _S152;
        _S104 = s_3;
        _S105 = x_11;
        _S106 = y_6;
        _S107 = _S150;
        _S115 = _S148;
        _S108 = _S146;
        _S109 = _S147;
        _S110 = _S144;
        _S111 = _S145;
        _S112 = _S142;
        _S113 = _S143;
        _S116 = _S140;
    }
    if(is_3dgs_1)
    {
        float _S159 = erank_reg_weight_s3_1 * (*_s_dOut_1)[int(3)];
        float _S160 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S161;
        (&_S161)->primal_0 = _S66;
        (&_S161)->differential_0 = 0.0f;
        DiffPair_float_0 _S162;
        (&_S162)->primal_0 = 0.0f;
        (&_S162)->differential_0 = 0.0f;
        _d_max_0(&_S161, &_S162, _S160);
        float _S163 = - _S161.differential_0;
        DiffPair_float_0 _S164;
        (&_S164)->primal_0 = _S67;
        (&_S164)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S164, _S163);
        DiffPair_float_0 _S165;
        (&_S165)->primal_0 = _S68;
        (&_S165)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S165, _S164.differential_0);
        float _S166 = - _S165.differential_0;
        float _S167 = _S69 * _S166;
        float _S168 = _S70 * _S166;
        DiffPair_float_0 _S169;
        (&_S169)->primal_0 = _S69;
        (&_S169)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S169, _S167);
        float _S170 = _S71 * _S166;
        float _S171 = _S72 * _S166;
        DiffPair_float_0 _S172;
        (&_S172)->primal_0 = _S71;
        (&_S172)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S172, _S170);
        float _S173 = _S73 * _S165.differential_0;
        float _S174 = _S74 * _S165.differential_0;
        DiffPair_float_0 _S175;
        (&_S175)->primal_0 = _S75;
        (&_S175)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S175, _S173);
        float _S176 = - _S174;
        float _S177 = - (_S171 + _S172.differential_0);
        float _S178 = (_S159 + _S168 + _S169.differential_0 + _S177) / _S76;
        float _S179 = _S77 * - _S178;
        float _S180 = _S78 * _S178;
        DiffPair_float_0 _S181;
        (&_S181)->primal_0 = _S79;
        (&_S181)->differential_0 = 0.0f;
        DiffPair_float_0 _S182;
        (&_S182)->primal_0 = _S80;
        (&_S182)->differential_0 = 0.0f;
        _d_min_0(&_S181, &_S182, _S180);
        DiffPair_float_0 _S183;
        (&_S183)->primal_0 = _S81;
        (&_S183)->differential_0 = 0.0f;
        DiffPair_float_0 _S184;
        (&_S184)->primal_0 = _S82;
        (&_S184)->differential_0 = 0.0f;
        _d_min_0(&_S183, &_S184, _S181.differential_0);
        float _S185 = (_S175.differential_0 + _S176 + _S177) / _S76;
        float _S186 = _S83 * - _S185;
        float _S187 = _S78 * _S185;
        DiffPair_float_0 _S188;
        (&_S188)->primal_0 = _S84;
        (&_S188)->differential_0 = 0.0f;
        DiffPair_float_0 _S189;
        (&_S189)->primal_0 = _S80;
        (&_S189)->differential_0 = 0.0f;
        _d_max_0(&_S188, &_S189, _S187);
        DiffPair_float_0 _S190;
        (&_S190)->primal_0 = _S81;
        (&_S190)->differential_0 = 0.0f;
        DiffPair_float_0 _S191;
        (&_S191)->primal_0 = _S82;
        (&_S191)->differential_0 = 0.0f;
        _d_max_0(&_S190, &_S191, _S188.differential_0);
        float _S192 = _S179 + _S186;
        float3  _S193 = make_float3 (_S183.differential_0 + _S190.differential_0 + _S192, _S184.differential_0 + _S191.differential_0 + _S192, _S182.differential_0 + _S189.differential_0 + _S192);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S194;
        (&_S194)->primal_0 = _S114;
        (&_S194)->differential_0 = _S60;
        s_bwd_prop_exp_1(&_S194, _S193);
        float3  _S195 = make_float3 (2.0f) * _S194.differential_0;
        float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S196;
        (&_S196)->primal_0 = _S85;
        (&_S196)->differential_0 = 0.0f;
        DiffPair_float_0 _S197;
        (&_S197)->primal_0 = max_gauss_ratio_1;
        (&_S197)->differential_0 = 0.0f;
        _d_max_0(&_S196, &_S197, s_diff_scale_reg_T_0);
        float _S198 = _S196.differential_0 / _S86;
        float _S199 = _S87 * - _S198;
        float _S200 = _S88 * _S198;
        DiffPair_float_0 _S201;
        (&_S201)->primal_0 = _S89;
        (&_S201)->differential_0 = 0.0f;
        DiffPair_float_0 _S202;
        (&_S202)->primal_0 = _S90;
        (&_S202)->differential_0 = 0.0f;
        _d_min_0(&_S201, &_S202, _S199);
        DiffPair_float_0 _S203;
        (&_S203)->primal_0 = _S91;
        (&_S203)->differential_0 = 0.0f;
        DiffPair_float_0 _S204;
        (&_S204)->primal_0 = _S92;
        (&_S204)->differential_0 = 0.0f;
        _d_min_0(&_S203, &_S204, _S201.differential_0);
        DiffPair_float_0 _S205;
        (&_S205)->primal_0 = _S93;
        (&_S205)->differential_0 = 0.0f;
        DiffPair_float_0 _S206;
        (&_S206)->primal_0 = _S90;
        (&_S206)->differential_0 = 0.0f;
        _d_max_0(&_S205, &_S206, _S200);
        DiffPair_float_0 _S207;
        (&_S207)->primal_0 = _S91;
        (&_S207)->differential_0 = 0.0f;
        DiffPair_float_0 _S208;
        (&_S208)->primal_0 = _S92;
        (&_S208)->differential_0 = 0.0f;
        _d_max_0(&_S207, &_S208, _S205.differential_0);
        float _S209 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * (*_s_dOut_1)[int(1)]);
        float3  _S210 = make_float3 (_S203.differential_0 + _S207.differential_0 + _S209, _S204.differential_0 + _S208.differential_0 + _S209, _S202.differential_0 + _S206.differential_0 + _S209);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S211;
        (&_S211)->primal_0 = _S58.primal_0;
        (&_S211)->differential_0 = _S60;
        s_bwd_prop_exp_1(&_S211, _S210);
        float3  _S212 = _S195 + _S211.differential_0;
        _S66 = (*_s_dOut_1)[int(4)];
        _S67 = (*_s_dOut_1)[int(0)];
        _S114 = _S212;
    }
    else
    {
        float _S213 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S214;
        (&_S214)->primal_0 = _S94;
        (&_S214)->differential_0 = 0.0f;
        DiffPair_float_0 _S215;
        (&_S215)->primal_0 = 0.0f;
        (&_S215)->differential_0 = 0.0f;
        _d_max_0(&_S214, &_S215, _S213);
        float _S216 = - _S214.differential_0;
        DiffPair_float_0 _S217;
        (&_S217)->primal_0 = _S95;
        (&_S217)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S217, _S216);
        DiffPair_float_0 _S218;
        (&_S218)->primal_0 = _S96;
        (&_S218)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S218, _S217.differential_0);
        float _S219 = - _S218.differential_0;
        float _S220 = _S97 * _S219;
        float _S221 = _S98 * _S219;
        DiffPair_float_0 _S222;
        (&_S222)->primal_0 = _S97;
        (&_S222)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S222, _S220);
        float _S223 = _S99 * _S218.differential_0;
        float _S224 = _S100 * _S218.differential_0;
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = _S101;
        (&_S225)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S225, _S223);
        float _S226 = - _S224;
        float _S227 = (_S221 + _S222.differential_0) / _S102;
        float _S228 = _S103 * - _S227;
        float _S229 = _S104 * _S227;
        DiffPair_float_0 _S230;
        (&_S230)->primal_0 = _S105;
        (&_S230)->differential_0 = 0.0f;
        DiffPair_float_0 _S231;
        (&_S231)->primal_0 = _S106;
        (&_S231)->differential_0 = 0.0f;
        _d_min_0(&_S230, &_S231, _S229);
        float _S232 = (_S225.differential_0 + _S226) / _S102;
        float _S233 = _S107 * - _S232;
        float _S234 = _S104 * _S232;
        DiffPair_float_0 _S235;
        (&_S235)->primal_0 = _S105;
        (&_S235)->differential_0 = 0.0f;
        DiffPair_float_0 _S236;
        (&_S236)->primal_0 = _S106;
        (&_S236)->differential_0 = 0.0f;
        _d_max_0(&_S235, &_S236, _S234);
        float _S237 = _S228 + _S233;
        float2  _S238 = make_float2 (_S230.differential_0 + _S235.differential_0 + _S237, _S231.differential_0 + _S236.differential_0 + _S237);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S239;
        (&_S239)->primal_0 = _S115;
        (&_S239)->differential_0 = _S61;
        s_bwd_prop_exp_2(&_S239, _S238);
        float2  _S240 = make_float2 (2.0f) * _S239.differential_0;
        float s_diff_scale_reg_T_1 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S241;
        (&_S241)->primal_0 = _S108;
        (&_S241)->differential_0 = 0.0f;
        DiffPair_float_0 _S242;
        (&_S242)->primal_0 = max_gauss_ratio_1;
        (&_S242)->differential_0 = 0.0f;
        _d_max_0(&_S241, &_S242, s_diff_scale_reg_T_1);
        float _S243 = _S241.differential_0 / _S109;
        float _S244 = _S110 * - _S243;
        float _S245 = _S111 * _S243;
        DiffPair_float_0 _S246;
        (&_S246)->primal_0 = _S112;
        (&_S246)->differential_0 = 0.0f;
        DiffPair_float_0 _S247;
        (&_S247)->primal_0 = _S113;
        (&_S247)->differential_0 = 0.0f;
        _d_min_0(&_S246, &_S247, _S244);
        DiffPair_float_0 _S248;
        (&_S248)->primal_0 = _S112;
        (&_S248)->differential_0 = 0.0f;
        DiffPair_float_0 _S249;
        (&_S249)->primal_0 = _S113;
        (&_S249)->differential_0 = 0.0f;
        _d_max_0(&_S248, &_S249, _S245);
        float _S250 = mcmc_scale_reg_weight_1 * (0.5f * (*_s_dOut_1)[int(1)]);
        float2  _S251 = make_float2 (_S246.differential_0 + _S248.differential_0 + _S250, _S247.differential_0 + _S249.differential_0 + _S250);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S252;
        (&_S252)->primal_0 = _S116;
        (&_S252)->differential_0 = _S61;
        s_bwd_prop_exp_2(&_S252, _S251);
        float2  _S253 = _S240 + _S252.differential_0;
        float3  _S254 = make_float3 (_S253.x, _S253.y, 0.0f);
        _S66 = (*_s_dOut_1)[int(4)];
        _S67 = (*_s_dOut_1)[int(0)];
        _S114 = _S254;
    }
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * _S66;
    float _S255 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S256;
    (&_S256)->primal_0 = _S65;
    (&_S256)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S256, _S255);
    float _S257 = _S256.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S258 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S259;
    (&_S259)->primal_0 = _S59.primal_0;
    (&_S259)->differential_0 = _S258;
    s_bwd_length_impl_0(&_S259, _S257);
    float _S260 = - (mcmc_opacity_reg_weight_1 * _S67 / _S64);
    DiffPair_float_0 _S261;
    (&_S261)->primal_0 = _S62;
    (&_S261)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S261, _S260);
    float _S262 = - _S261.differential_0;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S259.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S262;
    dpscales_0->primal_0 = (*dpscales_0).primal_0;
    dpscales_0->differential_0 = _S114;
    return;
}

inline __device__ void s_bwd_per_splat_losses_0(bool _S263, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S264, DiffPair_float_0 * _S265, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S266, float _S267, float _S268, float _S269, float _S270, float _S271, float _S272, float _S273, FixedArray<float, 5>  * _S274)
{
    s_bwd_prop_per_splat_losses_0(_S263, _S264, _S265, _S266, _S267, _S268, _S269, _S270, _S271, _S272, _S273, _S274);
    return;
}

inline __device__ void per_splat_losses_bwd(bool is_3dgs_2, float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float3  _S275 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_scales_0;
    (&p_scales_0)->primal_0 = scales_1;
    (&p_scales_0)->differential_0 = _S275;
    DiffPair_float_0 p_opacity_0;
    (&p_opacity_0)->primal_0 = opacity_1;
    (&p_opacity_0)->differential_0 = 0.0f;
    float4  _S276 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 p_quat_0;
    (&p_quat_0)->primal_0 = quat_1;
    (&p_quat_0)->differential_0 = _S276;
    FixedArray<float, 5>  _S277 = v_loss_0;
    s_bwd_per_splat_losses_0(is_3dgs_2, &p_scales_0, &p_opacity_0, &p_quat_0, mcmc_opacity_reg_weight_2, mcmc_scale_reg_weight_2, max_gauss_ratio_2, scale_regularization_weight_2, erank_reg_weight_2, erank_reg_weight_s3_2, quat_norm_reg_weight_2, &_S277);
    *v_scales_0 = p_scales_0.differential_0;
    *v_opacity_0 = p_opacity_0.differential_0;
    *v_quat_0 = p_quat_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_12)
{
    Matrix<float, 3, 3>  result_14;
    int r_1 = int(0);
    for(;;)
    {
        if(r_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_14)->rows + (r_1)), c_0) = _slang_vector_get_element(x_12.rows[c_0], r_1);
            c_0 = c_0 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_14;
}

inline __device__ Matrix<float, 2, 2>  transpose_1(Matrix<float, 2, 2>  x_13)
{
    Matrix<float, 2, 2>  result_15;
    int r_2 = int(0);
    for(;;)
    {
        if(r_2 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_1 = int(0);
        for(;;)
        {
            if(c_1 < int(2))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_15)->rows + (r_2)), c_1) = _slang_vector_get_element(x_13.rows[c_1], r_2);
            c_1 = c_1 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_15;
}

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_2)
{
    float x_14 = quat_2.y;
    float x2_0 = x_14 * x_14;
    float y2_0 = quat_2.z * quat_2.z;
    float z2_0 = quat_2.w * quat_2.w;
    float xy_0 = quat_2.y * quat_2.z;
    float xz_0 = quat_2.y * quat_2.w;
    float yz_0 = quat_2.z * quat_2.w;
    float wx_0 = quat_2.x * quat_2.y;
    float wy_0 = quat_2.x * quat_2.z;
    float wz_0 = quat_2.x * quat_2.w;
    return transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0)));
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_0, Matrix<float, 3, 3>  dOut_8)
{
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_0;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_8.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(0)].x * dOut_8.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_8.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(0)].y * dOut_8.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_8.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(0)].z * dOut_8.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_8.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(0)].x * dOut_8.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_8.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(0)].y * dOut_8.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_8.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(0)].z * dOut_8.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_8.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(0)].x * dOut_8.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_8.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(0)].y * dOut_8.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_8.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(0)].z * dOut_8.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_8.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(1)].x * dOut_8.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_8.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(1)].y * dOut_8.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_8.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(1)].z * dOut_8.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_8.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(1)].x * dOut_8.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_8.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(1)].y * dOut_8.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_8.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(1)].z * dOut_8.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_8.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(1)].x * dOut_8.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_8.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(1)].y * dOut_8.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_8.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(1)].z * dOut_8.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_8.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(2)].x * dOut_8.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_8.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(2)].y * dOut_8.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_8.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(2)].z * dOut_8.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_8.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(2)].x * dOut_8.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_8.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(2)].y * dOut_8.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_8.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(2)].z * dOut_8.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_8.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(2)].x * dOut_8.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_8.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(2)].y * dOut_8.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_8.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(2)].z * dOut_8.rows[int(2)].z;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
{
    Matrix<float, 3, 3>  result_16;
    int r_3 = int(0);
    for(;;)
    {
        if(r_3 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_2 = int(0);
        for(;;)
        {
            if(c_2 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_5 = int(0);
            float sum_0 = 0.0f;
            for(;;)
            {
                if(i_5 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_1 = sum_0 + _slang_vector_get_element(left_1.rows[r_3], i_5) * _slang_vector_get_element(right_1.rows[i_5], c_2);
                i_5 = i_5 + int(1);
                sum_0 = sum_1;
            }
            *_slang_vector_get_element_ptr(((&result_16)->rows + (r_3)), c_2) = sum_0;
            c_2 = c_2 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_16;
}

inline __device__ void quat_scale_to_covar(float4  quat_3, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_15 = quat_3.y;
    float x2_1 = x_15 * x_15;
    float y2_1 = quat_3.z * quat_3.z;
    float z2_1 = quat_3.w * quat_3.w;
    float xy_1 = quat_3.y * quat_3.z;
    float xz_1 = quat_3.y * quat_3.w;
    float yz_1 = quat_3.z * quat_3.w;
    float wx_1 = quat_3.x * quat_3.y;
    float wy_1 = quat_3.x * quat_3.z;
    float wz_1 = quat_3.x * quat_3.w;
    Matrix<float, 3, 3>  M_0 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_1(M_0, transpose_0(M_0));
    return;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_0, FixedArray<float, 10>  dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float _S278 = 0.0f * v_0;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float s_diff_r2_0 = u_0 + u_0 + (_S278 + _S278);
    float _S279 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
    float _S280 = dist_coeffs_0[int(1)] + r2_0 * _S279;
    float _S281 = dist_coeffs_0[int(0)] + r2_0 * _S280;
    float radial_0 = 1.0f + r2_0 * _S281;
    float _S282 = 2.0f * dist_coeffs_0[int(4)];
    float _S283 = _S282 * u_0;
    float _S284 = 2.0f * u_0;
    float _S285 = 2.0f * dist_coeffs_0[int(5)];
    float _S286 = _S285 * u_0;
    float _S287 = 2.0f * v_0;
    float2  _S288 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S281 + (s_diff_r2_0 * _S280 + (s_diff_r2_0 * _S279 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (_S282 * v_0 + 0.0f * _S283 + (s_diff_r2_0 + (_S284 + _S284)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S285 * v_0 + 0.0f * _S286 + (s_diff_r2_0 + (_S278 + 0.0f * _S287)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
    float _S289 = 0.0f * u_0;
    float s_diff_r2_1 = _S289 + _S289 + (v_0 + v_0);
    float2  _S290 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S281 + (s_diff_r2_1 * _S280 + (s_diff_r2_1 * _S279 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (0.0f * _S282 * v_0 + _S283 + (s_diff_r2_1 + (_S289 + 0.0f * _S284)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S285 * v_0 + _S286 + (s_diff_r2_1 + (_S287 + _S287)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
    Matrix<float, 2, 2>  _S291 = transpose_1(makeMatrix<float, 2, 2> (_S288 + make_float2 (_S288.x * dist_coeffs_0[int(8)] + _S288.y * dist_coeffs_0[int(9)], 0.0f), _S290 + make_float2 (_S290.x * dist_coeffs_0[int(8)] + _S290.y * dist_coeffs_0[int(9)], 0.0f)));
    return (F32_min((determinant_0(_S291)), ((F32_min((_S291.rows[int(0)].x), (_S291.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_1, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_1)
{
    float2  _S292;
    if(is_fisheye_0)
    {
        float r_4 = length_1(uv_1);
        float theta_0 = (F32_atan((r_4)));
        float _S293;
        if(r_4 < 0.00100000004749745f)
        {
            _S293 = 1.0f - r_4 * r_4 / 3.0f;
        }
        else
        {
            _S293 = theta_0 / r_4;
        }
        _S292 = uv_1 * make_float2 (_S293);
    }
    else
    {
        _S292 = uv_1;
    }
    float u_1 = _S292.x;
    float v_1 = _S292.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S294 = _S292 * make_float2 (1.0f + r2_1 * (dist_coeffs_1[int(0)] + r2_1 * (dist_coeffs_1[int(1)] + r2_1 * (dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_1 * v_1 + dist_coeffs_1[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_1[int(6)] * r2_1, 2.0f * dist_coeffs_1[int(5)] * u_1 * v_1 + dist_coeffs_1[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_1[int(7)] * r2_1);
    return _S294 + make_float2 (dist_coeffs_1[int(8)] * _S294.x + dist_coeffs_1[int(9)] * _S294.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_2, FixedArray<float, 10>  * dist_coeffs_2, int maxiter_0, float2  * uv_undist_0)
{
    int i_6 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_6 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S295 = (*dist_coeffs_2)[int(3)];
        float _S296 = (*dist_coeffs_2)[int(4)];
        float _S297 = (*dist_coeffs_2)[int(5)];
        float _S298 = (*dist_coeffs_2)[int(6)];
        float _S299 = (*dist_coeffs_2)[int(7)];
        float _S300 = (*dist_coeffs_2)[int(8)];
        float _S301 = (*dist_coeffs_2)[int(9)];
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float _S302 = (*dist_coeffs_2)[int(2)] + r2_2 * (*dist_coeffs_2)[int(3)];
        float _S303 = (*dist_coeffs_2)[int(1)] + r2_2 * _S302;
        float _S304 = (*dist_coeffs_2)[int(0)] + r2_2 * _S303;
        float radial_1 = 1.0f + r2_2 * _S304;
        float _S305 = 2.0f * (*dist_coeffs_2)[int(4)];
        float _S306 = _S305 * u_2;
        float _S307 = 2.0f * u_2;
        float _S308 = 2.0f * (*dist_coeffs_2)[int(5)];
        float _S309 = _S308 * u_2;
        float _S310 = 2.0f * v_2;
        float2  _S311 = q_0 * make_float2 (radial_1) + make_float2 (_S306 * v_2 + (*dist_coeffs_2)[int(5)] * (r2_2 + _S307 * u_2) + (*dist_coeffs_2)[int(6)] * r2_2, _S309 * v_2 + (*dist_coeffs_2)[int(4)] * (r2_2 + _S310 * v_2) + (*dist_coeffs_2)[int(7)] * r2_2);
        float2  r_5 = _S311 + make_float2 ((*dist_coeffs_2)[int(8)] * _S311.x + (*dist_coeffs_2)[int(9)] * _S311.y, 0.0f) - uv_2;
        float _S312 = 0.0f * v_2;
        float s_diff_r2_2 = u_2 + u_2 + (_S312 + _S312);
        float2  _S313 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S304 + (s_diff_r2_2 * _S303 + (s_diff_r2_2 * _S302 + s_diff_r2_2 * _S295 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (_S305 * v_2 + 0.0f * _S306 + (s_diff_r2_2 + (_S307 + _S307)) * _S297 + s_diff_r2_2 * _S298, _S308 * v_2 + 0.0f * _S309 + (s_diff_r2_2 + (_S312 + 0.0f * _S310)) * _S296 + s_diff_r2_2 * _S299);
        float _S314 = 0.0f * u_2;
        float s_diff_r2_3 = _S314 + _S314 + (v_2 + v_2);
        float2  _S315 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S304 + (s_diff_r2_3 * _S303 + (s_diff_r2_3 * _S302 + s_diff_r2_3 * _S295 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (0.0f * _S305 * v_2 + _S306 + (s_diff_r2_3 + (_S314 + 0.0f * _S307)) * _S297 + s_diff_r2_3 * _S298, 0.0f * _S308 * v_2 + _S309 + (s_diff_r2_3 + (_S310 + _S310)) * _S296 + s_diff_r2_3 * _S299);
        Matrix<float, 2, 2>  _S316 = transpose_1(makeMatrix<float, 2, 2> (_S313 + make_float2 (_S313.x * _S300 + _S313.y * _S301, 0.0f), _S315 + make_float2 (_S315.x * _S300 + _S315.y * _S301, 0.0f)));
        float inv_det_0 = 1.0f / (_S316.rows[int(0)].x * _S316.rows[int(1)].y - _S316.rows[int(0)].y * _S316.rows[int(1)].x);
        float _S317 = r_5.x;
        float _S318 = r_5.y;
        float2  q_1 = q_0 - make_float2 ((_S317 * _S316.rows[int(1)].y - _S318 * _S316.rows[int(0)].y) * inv_det_0, (- _S317 * _S316.rows[int(1)].x + _S318 * _S316.rows[int(0)].x) * inv_det_0);
        i_6 = i_6 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S319 = (*dist_coeffs_2)[int(0)];
    float _S320 = (*dist_coeffs_2)[int(1)];
    float _S321 = (*dist_coeffs_2)[int(2)];
    float _S322 = (*dist_coeffs_2)[int(3)];
    float _S323 = (*dist_coeffs_2)[int(4)];
    float _S324 = (*dist_coeffs_2)[int(5)];
    float _S325 = (*dist_coeffs_2)[int(6)];
    float _S326 = (*dist_coeffs_2)[int(7)];
    float _S327 = (*dist_coeffs_2)[int(8)];
    float _S328 = (*dist_coeffs_2)[int(9)];
    float u_3 = q_0.x;
    float v_3 = q_0.y;
    float _S329 = 0.0f * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = u_3 + u_3 + (_S329 + _S329);
    float _S330 = (*dist_coeffs_2)[int(2)] + r2_3 * (*dist_coeffs_2)[int(3)];
    float _S331 = (*dist_coeffs_2)[int(1)] + r2_3 * _S330;
    float _S332 = (*dist_coeffs_2)[int(0)] + r2_3 * _S331;
    float radial_2 = 1.0f + r2_3 * _S332;
    float _S333 = 2.0f * (*dist_coeffs_2)[int(4)];
    float _S334 = _S333 * u_3;
    float _S335 = 2.0f * u_3;
    float _S336 = 2.0f * (*dist_coeffs_2)[int(5)];
    float _S337 = _S336 * u_3;
    float _S338 = 2.0f * v_3;
    float2  _S339 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S332 + (s_diff_r2_4 * _S331 + (s_diff_r2_4 * _S330 + s_diff_r2_4 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (_S333 * v_3 + 0.0f * _S334 + (s_diff_r2_4 + (_S335 + _S335)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_4 * (*dist_coeffs_2)[int(6)], _S336 * v_3 + 0.0f * _S337 + (s_diff_r2_4 + (_S329 + 0.0f * _S338)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_4 * (*dist_coeffs_2)[int(7)]);
    float _S340 = 0.0f * u_3;
    float s_diff_r2_5 = _S340 + _S340 + (v_3 + v_3);
    float2  _S341 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S332 + (s_diff_r2_5 * _S331 + (s_diff_r2_5 * _S330 + s_diff_r2_5 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (0.0f * _S333 * v_3 + _S334 + (s_diff_r2_5 + (_S340 + 0.0f * _S335)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_5 * (*dist_coeffs_2)[int(6)], 0.0f * _S336 * v_3 + _S337 + (s_diff_r2_5 + (_S338 + _S338)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_5 * (*dist_coeffs_2)[int(7)]);
    Matrix<float, 2, 2>  _S342 = transpose_1(makeMatrix<float, 2, 2> (_S339 + make_float2 (_S339.x * (*dist_coeffs_2)[int(8)] + _S339.y * (*dist_coeffs_2)[int(9)], 0.0f), _S341 + make_float2 (_S341.x * (*dist_coeffs_2)[int(8)] + _S341.y * (*dist_coeffs_2)[int(9)], 0.0f)));
    bool _S343;
    if((F32_min((determinant_0(_S342)), ((F32_min((_S342.rows[int(0)].x), (_S342.rows[int(1)].y)))))) > 0.0f)
    {
        float u_4 = (*uv_undist_0).x;
        float v_4 = (*uv_undist_0).y;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float2  _S344 = *uv_undist_0 * make_float2 (1.0f + r2_4 * (_S319 + r2_4 * (_S320 + r2_4 * (_S321 + r2_4 * _S322)))) + make_float2 (_S333 * u_4 * v_4 + _S324 * (r2_4 + 2.0f * u_4 * u_4) + _S325 * r2_4, _S336 * u_4 * v_4 + _S323 * (r2_4 + 2.0f * v_4 * v_4) + _S326 * r2_4);
        _S343 = (length_1(_S344 + make_float2 (_S327 * _S344.x + _S328 * _S344.y, 0.0f) - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S343 = false;
    }
    return _S343;
}

inline __device__ bool undistort_point(float2  uv_3, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_undist_1)
{
    float2  _S345 = uv_3;
    FixedArray<float, 10>  _S346 = dist_coeffs_3;
    bool _S347 = undistort_point_0(uv_3, &_S346, int(8), &_S345);
    if(!_S347)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S348 = _S345;
        float theta_1 = length_1(_S345);
        float _S349;
        if(theta_1 < 0.00100000004749745f)
        {
            _S349 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S349 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S350 = make_float3 ((_S348 * make_float2 (_S349)).x, (_S348 * make_float2 (_S349)).y, (F32_cos((theta_1))));
        raydir_0 = _S350;
    }
    else
    {
        raydir_0 = make_float3 (_S345.x, _S345.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_4, float3  * raydir_1)
{
    float2  _S351 = uv_4;
    int3  _S352 = make_int3 (int(0));
    float3  _S353 = make_float3 ((float)_S352.x, (float)_S352.y, (float)_S352.z);
    *raydir_1 = _S353;
    FixedArray<float, 10>  _S354 = dist_coeffs_4;
    bool _S355 = undistort_point_0(uv_4, &_S354, int(8), &_S351);
    if(!_S355)
    {
        return false;
    }
    float3  _S356;
    if(is_fisheye_2)
    {
        float2  _S357 = _S351;
        float theta_2 = length_1(_S351);
        float _S358;
        if(theta_2 < 0.00100000004749745f)
        {
            _S358 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S358 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S359 = make_float3 ((_S357 * make_float2 (_S358)).x, (_S357 * make_float2 (_S358)).y, (F32_cos((theta_2))));
        _S356 = _S359;
    }
    else
    {
        _S356 = make_float3 (_S351.x, _S351.y, 1.0f);
    }
    *raydir_1 = _S356;
    return true;
}

inline __device__ float3  normalize_0(float3  x_16)
{
    return x_16 / make_float3 (length_2(x_16));
}

inline __device__ float4  normalize_1(float4  x_17)
{
    return x_17 / make_float4 (length_0(x_17));
}

inline __device__ bool generate_ray(float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_5, float3  * raydir_2)
{
    float2  _S360 = uv_5;
    FixedArray<float, 10>  _S361 = dist_coeffs_5;
    bool _S362 = undistort_point_0(uv_5, &_S361, int(8), &_S360);
    if(!_S362)
    {
        int3  _S363 = make_int3 (int(0));
        float3  _S364 = make_float3 ((float)_S363.x, (float)_S363.y, (float)_S363.z);
        *raydir_2 = _S364;
        return false;
    }
    float3  _S365;
    if(is_fisheye_3)
    {
        float2  _S366 = _S360;
        float theta_3 = length_1(_S360);
        float _S367;
        if(theta_3 < 0.00100000004749745f)
        {
            _S367 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S367 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S368 = make_float3 ((_S366 * make_float2 (_S367)).x, (_S366 * make_float2 (_S367)).y, (F32_cos((theta_3))));
        _S365 = _S368;
    }
    else
    {
        _S365 = make_float3 (_S360.x, _S360.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S365);
    return true;
}

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_9)
{
    float _S369 = (*right_2).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_9.x;
    float sum_2 = _S369 + (*right_2).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_9.y;
    float sum_3 = sum_2 + (*right_2).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_9.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_3;
    float _S370 = (*right_2).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_9.x;
    float sum_4 = _S370 + (*right_2).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_9.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_9.z;
    *&((&left_d_result_1)->y) = sum_5;
    float _S371 = (*right_2).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_9.x;
    float sum_6 = _S371 + (*right_2).primal_0.rows[int(2)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_9.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(2)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_9.z;
    *&((&left_d_result_1)->z) = sum_7;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_2(float3  left_3, Matrix<float, 3, 3>  right_3)
{
    float3  result_17;
    int j_0 = int(0);
    for(;;)
    {
        if(j_0 < int(3))
        {
        }
        else
        {
            break;
        }
        int i_7 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(i_7 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_3, i_7) * _slang_vector_get_element(right_3.rows[i_7], j_0);
            i_7 = i_7 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_17, j_0) = sum_8;
        j_0 = j_0 + int(1);
    }
    return result_17;
}

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_0, float3  t_0)
{
    return - mul_2(t_0, R_0);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_1, float3  raydir_3)
{
    return mul_2(raydir_3, R_1);
}

inline __device__ float3  undo_transform_ray_d(Matrix<float, 3, 3>  R_2, float3  raydir_4)
{
    return mul_2(raydir_4, transpose_0(R_2));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S372, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S373, float3  _S374)
{
    _d_mul_0(_S372, _S373, _S374);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S375 = - _s_dOut_2;
    float3  _S376 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S377;
    (&_S377)->primal_0 = (*dpt_0).primal_0;
    (&_S377)->differential_0 = _S376;
    Matrix<float, 3, 3>  _S378 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S379;
    (&_S379)->primal_0 = (*dpR_0).primal_0;
    (&_S379)->differential_0 = _S378;
    s_bwd_prop_mul_0(&_S377, &_S379, _S375);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S377.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S379.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S380, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S381, float3  _S382)
{
    s_bwd_prop_transform_ray_o_0(_S380, _S381, _S382);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_3, float3  t_1, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S383 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S383;
    float3  _S384 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_1;
    (&dp_t_0)->differential_0 = _S384;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S385 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S386;
    (&_S386)->primal_0 = (*dpraydir_0).primal_0;
    (&_S386)->differential_0 = _S385;
    Matrix<float, 3, 3>  _S387 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S388;
    (&_S388)->primal_0 = (*dpR_1).primal_0;
    (&_S388)->differential_0 = _S387;
    s_bwd_prop_mul_0(&_S386, &_S388, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S386.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S388.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S389, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S390, float3  _S391)
{
    s_bwd_prop_transform_ray_d_0(_S389, _S390, _S391);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_4, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S392 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_4;
    (&dp_R_1)->differential_0 = _S392;
    float3  _S393 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S393;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_4, float3  scale_1)
{
    float x_18 = quat_4.y;
    float x2_2 = x_18 * x_18;
    float y2_2 = quat_4.z * quat_4.z;
    float z2_2 = quat_4.w * quat_4.w;
    float xy_2 = quat_4.y * quat_4.z;
    float xz_2 = quat_4.y * quat_4.w;
    float yz_2 = quat_4.z * quat_4.w;
    float wx_2 = quat_4.x * quat_4.y;
    float wy_2 = quat_4.x * quat_4.z;
    float wz_2 = quat_4.x * quat_4.w;
    return mul_1(makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2)))));
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S394, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S395, Matrix<float, 3, 3>  _S396)
{
    mul_0(_S394, _S395, _S396);
    return;
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_4)
{
    float _S397 = (*dpquat_1).primal_0.y;
    float x2_3 = _S397 * _S397;
    float y2_3 = (*dpquat_1).primal_0.z * (*dpquat_1).primal_0.z;
    float z2_3 = (*dpquat_1).primal_0.w * (*dpquat_1).primal_0.w;
    float xy_3 = (*dpquat_1).primal_0.y * (*dpquat_1).primal_0.z;
    float xz_3 = (*dpquat_1).primal_0.y * (*dpquat_1).primal_0.w;
    float yz_3 = (*dpquat_1).primal_0.z * (*dpquat_1).primal_0.w;
    float wx_3 = (*dpquat_1).primal_0.x * (*dpquat_1).primal_0.y;
    float wy_3 = (*dpquat_1).primal_0.x * (*dpquat_1).primal_0.z;
    float wz_3 = (*dpquat_1).primal_0.x * (*dpquat_1).primal_0.w;
    Matrix<float, 3, 3>  _S398 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))));
    Matrix<float, 3, 3>  _S399 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S400;
    (&_S400)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S400)->differential_0 = _S399;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S401;
    (&_S401)->primal_0 = _S398;
    (&_S401)->differential_0 = _S399;
    s_bwd_prop_mul_1(&_S400, &_S401, _s_dOut_4);
    Matrix<float, 3, 3>  _S402 = transpose_0(transpose_0(_S401.differential_0));
    float _S403 = 2.0f * - _S402.rows[int(2)].z;
    float _S404 = 2.0f * _S402.rows[int(2)].y;
    float _S405 = 2.0f * _S402.rows[int(2)].x;
    float _S406 = 2.0f * _S402.rows[int(1)].z;
    float _S407 = 2.0f * - _S402.rows[int(1)].y;
    float _S408 = 2.0f * _S402.rows[int(1)].x;
    float _S409 = 2.0f * _S402.rows[int(0)].z;
    float _S410 = 2.0f * _S402.rows[int(0)].y;
    float _S411 = 2.0f * - _S402.rows[int(0)].x;
    float _S412 = - _S408 + _S410;
    float _S413 = _S405 + - _S409;
    float _S414 = - _S404 + _S406;
    float _S415 = _S404 + _S406;
    float _S416 = _S405 + _S409;
    float _S417 = _S408 + _S410;
    float _S418 = (*dpquat_1).primal_0.w * (_S407 + _S411);
    float _S419 = (*dpquat_1).primal_0.z * (_S403 + _S411);
    float _S420 = (*dpquat_1).primal_0.y * (_S403 + _S407);
    float _S421 = (*dpquat_1).primal_0.x * _S412 + (*dpquat_1).primal_0.z * _S415 + (*dpquat_1).primal_0.y * _S416 + _S418 + _S418;
    float _S422 = (*dpquat_1).primal_0.x * _S413 + (*dpquat_1).primal_0.w * _S415 + (*dpquat_1).primal_0.y * _S417 + _S419 + _S419;
    float _S423 = (*dpquat_1).primal_0.x * _S414 + (*dpquat_1).primal_0.w * _S416 + (*dpquat_1).primal_0.z * _S417 + _S420 + _S420;
    float _S424 = (*dpquat_1).primal_0.w * _S412 + (*dpquat_1).primal_0.z * _S413 + (*dpquat_1).primal_0.y * _S414;
    float3  _S425 = make_float3 (_S400.differential_0.rows[int(0)].x, _S400.differential_0.rows[int(1)].y, _S400.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S425;
    float4  _S426 = make_float4 (0.0f);
    *&((&_S426)->w) = _S421;
    *&((&_S426)->z) = _S422;
    *&((&_S426)->y) = _S423;
    *&((&_S426)->x) = _S424;
    dpquat_1->primal_0 = (*dpquat_1).primal_0;
    dpquat_1->differential_0 = _S426;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S427, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S428, Matrix<float, 3, 3>  _S429)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S427, _S428, _S429);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_5, float3  scale_2, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_1, float3  * v_scale_0)
{
    float4  _S430 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_5;
    (&dp_quat_0)->differential_0 = _S430;
    float3  _S431 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_2;
    (&dp_scale_0)->differential_0 = _S431;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_1 = dp_quat_0.differential_0;
    *v_scale_0 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_4, float3  dOut_10)
{
    float _S432 = (*left_4).primal_0.rows[int(0)].x * dOut_10.x;
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = (*right_4).primal_0.x * dOut_10.x;
    float sum_10 = _S432 + (*left_4).primal_0.rows[int(1)].x * dOut_10.y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = (*right_4).primal_0.x * dOut_10.y;
    float sum_11 = sum_10 + (*left_4).primal_0.rows[int(2)].x * dOut_10.z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = (*right_4).primal_0.x * dOut_10.z;
    float3  right_d_result_2;
    *&((&right_d_result_2)->x) = sum_11;
    float _S433 = (*left_4).primal_0.rows[int(0)].y * dOut_10.x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = (*right_4).primal_0.y * dOut_10.x;
    float sum_12 = _S433 + (*left_4).primal_0.rows[int(1)].y * dOut_10.y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = (*right_4).primal_0.y * dOut_10.y;
    float sum_13 = sum_12 + (*left_4).primal_0.rows[int(2)].y * dOut_10.z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = (*right_4).primal_0.y * dOut_10.z;
    *&((&right_d_result_2)->y) = sum_13;
    float _S434 = (*left_4).primal_0.rows[int(0)].z * dOut_10.x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = (*right_4).primal_0.z * dOut_10.x;
    float sum_14 = _S434 + (*left_4).primal_0.rows[int(1)].z * dOut_10.y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = (*right_4).primal_0.z * dOut_10.y;
    float sum_15 = sum_14 + (*left_4).primal_0.rows[int(2)].z * dOut_10.z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = (*right_4).primal_0.z * dOut_10.z;
    *&((&right_d_result_2)->z) = sum_15;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_2;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_2;
    return;
}

struct DiffPair_matrixx3Cfloatx2C2x2C2x3E_0
{
    Matrix<float, 2, 2>  primal_0;
    Matrix<float, 2, 2>  differential_0;
};

inline __device__ void _d_mul_2(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * left_5, DiffPair_vectorx3Cfloatx2C2x3E_0 * right_5, float2  dOut_11)
{
    float _S435 = (*left_5).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 2, 2>  left_d_result_3;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = (*right_5).primal_0.x * dOut_11.x;
    float sum_16 = _S435 + (*left_5).primal_0.rows[int(1)].x * dOut_11.y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = (*right_5).primal_0.x * dOut_11.y;
    float2  right_d_result_3;
    *&((&right_d_result_3)->x) = sum_16;
    float _S436 = (*left_5).primal_0.rows[int(0)].y * dOut_11.x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = (*right_5).primal_0.y * dOut_11.x;
    float sum_17 = _S436 + (*left_5).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = (*right_5).primal_0.y * dOut_11.y;
    *&((&right_d_result_3)->y) = sum_17;
    left_5->primal_0 = (*left_5).primal_0;
    left_5->differential_0 = left_d_result_3;
    right_5->primal_0 = (*right_5).primal_0;
    right_5->differential_0 = right_d_result_3;
    return;
}

inline __device__ float3  mul_3(Matrix<float, 3, 3>  left_6, float3  right_6)
{
    float3  result_18;
    int i_8 = int(0);
    for(;;)
    {
        if(i_8 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_1 = int(0);
        float sum_18 = 0.0f;
        for(;;)
        {
            if(j_1 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_19 = sum_18 + _slang_vector_get_element(left_6.rows[i_8], j_1) * _slang_vector_get_element(right_6, j_1);
            j_1 = j_1 + int(1);
            sum_18 = sum_19;
        }
        *_slang_vector_get_element_ptr(&result_18, i_8) = sum_18;
        i_8 = i_8 + int(1);
    }
    return result_18;
}

inline __device__ float2  mul_4(Matrix<float, 2, 2>  left_7, float2  right_7)
{
    float2  result_19;
    int i_9 = int(0);
    for(;;)
    {
        if(i_9 < int(2))
        {
        }
        else
        {
            break;
        }
        int j_2 = int(0);
        float sum_20 = 0.0f;
        for(;;)
        {
            if(j_2 < int(2))
            {
            }
            else
            {
                break;
            }
            float sum_21 = sum_20 + _slang_vector_get_element(left_7.rows[i_9], j_2) * _slang_vector_get_element(right_7, j_2);
            j_2 = j_2 + int(1);
            sum_20 = sum_21;
        }
        *_slang_vector_get_element_ptr(&result_19, i_9) = sum_20;
        i_9 = i_9 + int(1);
    }
    return result_19;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_12)
{
    float _S437 = dOut_12.y;
    float _S438 = dOut_12.z;
    float _S439 = dOut_12.x;
    float _S440 = (*a_0).primal_0.z * _S437 + - (*a_0).primal_0.y * _S438;
    float _S441 = - (*a_0).primal_0.z * _S439 + (*a_0).primal_0.x * _S438;
    float _S442 = (*a_0).primal_0.y * _S439 + - (*a_0).primal_0.x * _S437;
    float3  _S443 = make_float3 (- (*b_1).primal_0.z * _S437 + (*b_1).primal_0.y * _S438, (*b_1).primal_0.z * _S439 + - (*b_1).primal_0.x * _S438, - (*b_1).primal_0.y * _S439 + (*b_1).primal_0.x * _S437);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S443;
    float3  _S444 = make_float3 (_S440, _S441, _S442);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S444;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S445 = left_8.y;
    float _S446 = right_8.z;
    float _S447 = left_8.z;
    float _S448 = right_8.y;
    float _S449 = right_8.x;
    float _S450 = left_8.x;
    return make_float3 (_S445 * _S446 - _S447 * _S448, _S447 * _S449 - _S450 * _S446, _S450 * _S448 - _S445 * _S449);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_0, Matrix<float, 3, 3>  iscl_rot_0, float opacity_2, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_3(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_3(iscl_rot_0, ray_o_0 - mean_0));
    return opacity_2 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S451, float3  _S452)
{
    return mul_3(_S451, _S452);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S453, float3  _S454)
{
    return cross_0(_S453, _S454);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S455, float3  _S456)
{
    return dot_0(_S455, _S456);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S457, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S458, float _S459)
{
    _d_dot_0(_S457, _S458, _S459);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S460, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S461, float3  _S462)
{
    _d_cross_0(_S460, _S461, _S462);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S463, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S464, float3  _S465)
{
    _d_mul_1(_S463, _S464, _S465);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_5)
{
    float3  _S466 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S467 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S466);
    float3  _S468 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S469 = s_primal_ctx_cross_0(_S468, _S467);
    float _S470 = -0.5f * s_primal_ctx_dot_0(_S469, _S469);
    float _S471 = s_primal_ctx_dot_0(_S468, _S468);
    float _S472 = _S470 / _S471;
    float _S473 = _S471 * _S471;
    float _S474 = (*dpopacity_1).primal_0 * _s_dOut_5;
    float _S475 = s_primal_ctx_exp_0(_S472) * _s_dOut_5;
    DiffPair_float_0 _S476;
    (&_S476)->primal_0 = _S472;
    (&_S476)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S476, _S474);
    float _S477 = _S476.differential_0 / _S473;
    float _S478 = _S470 * - _S477;
    float _S479 = _S471 * _S477;
    float3  _S480 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481;
    (&_S481)->primal_0 = _S468;
    (&_S481)->differential_0 = _S480;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S482;
    (&_S482)->primal_0 = _S468;
    (&_S482)->differential_0 = _S480;
    s_bwd_prop_dot_0(&_S481, &_S482, _S478);
    float _S483 = -0.5f * _S479;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S469;
    (&_S484)->differential_0 = _S480;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S485;
    (&_S485)->primal_0 = _S469;
    (&_S485)->differential_0 = _S480;
    s_bwd_prop_dot_0(&_S484, &_S485, _S483);
    float3  _S486 = _S485.differential_0 + _S484.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S487;
    (&_S487)->primal_0 = _S468;
    (&_S487)->differential_0 = _S480;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S488;
    (&_S488)->primal_0 = _S467;
    (&_S488)->differential_0 = _S480;
    s_bwd_prop_cross_0(&_S487, &_S488, _S486);
    float3  _S489 = _S482.differential_0 + _S481.differential_0 + _S487.differential_0;
    Matrix<float, 3, 3>  _S490 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S491;
    (&_S491)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S491)->differential_0 = _S490;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S492;
    (&_S492)->primal_0 = (*dpray_d_0).primal_0;
    (&_S492)->differential_0 = _S480;
    s_bwd_prop_mul_2(&_S491, &_S492, _S489);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S493;
    (&_S493)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S493)->differential_0 = _S490;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S494;
    (&_S494)->primal_0 = _S466;
    (&_S494)->differential_0 = _S480;
    s_bwd_prop_mul_2(&_S493, &_S494, _S488.differential_0);
    float3  _S495 = - _S494.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S492.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S494.differential_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = _S475;
    Matrix<float, 3, 3>  _S496 = _S491.differential_0 + _S493.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S496;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S495;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S497, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S498, DiffPair_float_0 * _S499, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S500, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S501, float _S502)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S497, _S498, _S499, _S500, _S501, _S502);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_1, Matrix<float, 3, 3>  iscl_rot_1, float opacity_3, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_0, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_1, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S503 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_1;
    (&dp_mean_0)->differential_0 = _S503;
    Matrix<float, 3, 3>  _S504 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S504;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_3;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S503;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_iscl_rot_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_0 = dp_mean_0.differential_0;
    *v_iscl_rot_1 = dp_iscl_rot_0.differential_0;
    *v_opacity_1 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_2, Matrix<float, 3, 3>  iscl_rot_2, float opacity_4, float3  rgb_0, float3  ray_o_2, float3  ray_d_2, float3  * out_rgb_0, float * depth_0)
{
    *out_rgb_0 = rgb_0;
    float3  grd_1 = mul_3(iscl_rot_2, ray_d_2);
    *depth_0 = - dot_0(mul_3(iscl_rot_2, ray_o_2 - mean_2), grd_1) / dot_0(grd_1, grd_1);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S505 = (*dpray_o_1).primal_0 - (*dpmean_1).primal_0;
    float3  _S506 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, _S505);
    float3  _S507 = s_primal_ctx_mul_0((*dpiscl_rot_1).primal_0, (*dpray_d_1).primal_0);
    float _S508 = s_primal_ctx_dot_0(_S507, _S507);
    float _S509 = dpdepth_0 / (_S508 * _S508);
    float _S510 = - s_primal_ctx_dot_0(_S506, _S507) * - _S509;
    float _S511 = _S508 * _S509;
    float3  _S512 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S513;
    (&_S513)->primal_0 = _S507;
    (&_S513)->differential_0 = _S512;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S514;
    (&_S514)->primal_0 = _S507;
    (&_S514)->differential_0 = _S512;
    s_bwd_prop_dot_0(&_S513, &_S514, _S510);
    float _S515 = - _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S516;
    (&_S516)->primal_0 = _S506;
    (&_S516)->differential_0 = _S512;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S517;
    (&_S517)->primal_0 = _S507;
    (&_S517)->differential_0 = _S512;
    s_bwd_prop_dot_0(&_S516, &_S517, _S515);
    float3  _S518 = _S514.differential_0 + _S513.differential_0 + _S517.differential_0;
    Matrix<float, 3, 3>  _S519 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S520;
    (&_S520)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S520)->differential_0 = _S519;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S521;
    (&_S521)->primal_0 = (*dpray_d_1).primal_0;
    (&_S521)->differential_0 = _S512;
    s_bwd_prop_mul_2(&_S520, &_S521, _S518);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S522;
    (&_S522)->primal_0 = (*dpiscl_rot_1).primal_0;
    (&_S522)->differential_0 = _S519;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S523;
    (&_S523)->primal_0 = _S505;
    (&_S523)->differential_0 = _S512;
    s_bwd_prop_mul_2(&_S522, &_S523, _S516.differential_0);
    float3  _S524 = - _S523.differential_0;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S521.differential_0;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S523.differential_0;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_2->primal_0 = (*dpopacity_2).primal_0;
    dpopacity_2->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S525 = _S520.differential_0 + _S522.differential_0;
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S525;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S524;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S526, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S527, DiffPair_float_0 * _S528, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S529, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S530, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S531, float3  _S532, float _S533)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S526, _S527, _S528, _S529, _S530, _S531, _S532, _S533);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_3, Matrix<float, 3, 3>  iscl_rot_3, float opacity_5, float3  rgb_1, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_0, float3  * v_mean_1, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_2, float3  * v_rgb_0, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S534 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_3;
    (&dp_mean_1)->differential_0 = _S534;
    Matrix<float, 3, 3>  _S535 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S535;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_5;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_1;
    (&dp_rgb_0)->differential_0 = _S534;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S534;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S534;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_0);
    *v_mean_1 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_2 = dp_opacity_1.differential_0;
    *v_rgb_0 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_4, float4  quat_6, float3  scale_3, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S536 = scale_3.x;
    float sx_0 = (F32_exp((_S536)));
    float _S537 = scale_3.y;
    float sy_0 = (F32_exp((_S537)));
    float sz_0 = scale_3.z - 0.5f * (_S536 + _S537);
    float x_19 = quat_6.y;
    float x2_4 = x_19 * x_19;
    float y2_4 = quat_6.z * quat_6.z;
    float z2_4 = quat_6.w * quat_6.w;
    float xy_4 = quat_6.y * quat_6.z;
    float xz_4 = quat_6.y * quat_6.w;
    float yz_4 = quat_6.z * quat_6.w;
    float wx_4 = quat_6.x * quat_6.y;
    float wy_4 = quat_6.x * quat_6.z;
    float wz_4 = quat_6.x * quat_6.w;
    Matrix<float, 3, 3>  _S538 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4)));
    *vert0_0 = mul_3(_S538, make_float3 (sx_0, 0.0f, 0.0f)) + mean_4;
    *vert1_0 = mul_3(_S538, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_4;
    *vert2_0 = mul_3(_S538, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_4;
    return;
}

inline __device__ float4  floor_0(float4  x_20)
{
    float4  result_20;
    int i_10 = int(0);
    for(;;)
    {
        if(i_10 < int(4))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_20, i_10) = (F32_floor((_slang_vector_get_element(x_20, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_20;
}

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_5, float3  scale_4, float4  quat_7, float opac_0)
{
    float4  _S539 = normalize_1(quat_7);
    float3  _S540 = exp_0(scale_4);
    float x_21 = _S539.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S539.z * _S539.z;
    float z2_5 = _S539.w * _S539.w;
    float xy_5 = _S539.y * _S539.z;
    float xz_5 = _S539.y * _S539.w;
    float yz_5 = _S539.z * _S539.w;
    float wx_5 = _S539.x * _S539.y;
    float wy_5 = _S539.x * _S539.z;
    float wz_5 = _S539.x * _S539.w;
    Matrix<float, 3, 3>  M_1 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S540.x, 0.0f, 0.0f, 0.0f, _S540.y, 0.0f, 0.0f, 0.0f, _S540.z));
    float4  _S541 = make_float4 (dot_0(*mean_5, *mean_5), dot_0(*mean_5, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S542 = _S541 - floor_0(_S541);
    float4  _S543 = _S542 + make_float4 (dot_1(_S542, float4 {_S542.w, _S542.z, _S542.x, _S542.y} + make_float4 (33.3300018310546875f)));
    float4  _S544 = (float4 {_S543.x, _S543.x, _S543.y, _S543.z} + float4 {_S543.y, _S543.z, _S543.z, _S543.w}) * float4 {_S543.z, _S543.y, _S543.w, _S543.x};
    float4  _S545 = _S544 - floor_0(_S544);
    float2  _S546 = float2 {_S545.x, _S545.z};
    float _S547 = 6.28318548202514648f * _S546.y;
    float2  _S548 = float2 {_S545.y, _S545.w};
    float _S549 = 6.28318548202514648f * _S548.y;
    *mean_5 = *mean_5 + mul_3(mul_1(M_1, transpose_0(M_1)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S546.x))))))) * make_float2 ((F32_cos((_S547))), (F32_sin((_S547))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S546.x))))))) * make_float2 ((F32_cos((_S547))), (F32_sin((_S547))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S548.x))))))) * make_float2 ((F32_cos((_S549))), (F32_sin((_S549))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_6, float3  scale_5, float4  quat_8, float opac_1)
{
    float4  _S550 = normalize_1(quat_8);
    float _S551 = scale_5.x;
    float sx_1 = (F32_exp((_S551)));
    float _S552 = scale_5.y;
    float sy_1 = (F32_exp((_S552)));
    float sz_1 = scale_5.z - 0.5f * (_S551 + _S552);
    float x_22 = _S550.y;
    float x2_6 = x_22 * x_22;
    float y2_6 = _S550.z * _S550.z;
    float z2_6 = _S550.w * _S550.w;
    float xy_6 = _S550.y * _S550.z;
    float xz_6 = _S550.y * _S550.w;
    float yz_6 = _S550.z * _S550.w;
    float wx_6 = _S550.x * _S550.y;
    float wy_6 = _S550.x * _S550.z;
    float wz_6 = _S550.x * _S550.w;
    Matrix<float, 3, 3>  _S553 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6)));
    float3  vert0_1 = mul_3(_S553, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_6;
    float3  vert1_1 = mul_3(_S553, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_6;
    float3  vert2_1 = mul_3(_S553, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_6;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S554 = make_float4 (dot_0(*mean_6, *mean_6), dot_0(*mean_6, scale_5), dot_0(scale_5, scale_5), dot_1(quat_8, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S555 = _S554 - floor_0(_S554);
    float4  _S556 = _S555 + make_float4 (dot_1(_S555, float4 {_S555.w, _S555.z, _S555.x, _S555.y} + make_float4 (33.3300018310546875f)));
    float4  _S557 = (float4 {_S556.x, _S556.x, _S556.y, _S556.z} + float4 {_S556.y, _S556.z, _S556.z, _S556.w}) * float4 {_S556.z, _S556.y, _S556.w, _S556.x};
    float4  _S558 = _S557 - floor_0(_S557);
    float2  _S559 = float2 {_S558.x, _S558.z};
    float _S560 = 6.28318548202514648f * _S559.y;
    float2  _S561 = float2 {_S558.y, _S558.w};
    float _S562 = 6.28318548202514648f * _S561.y;
    *mean_6 = *mean_6 + mul_3(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S559.x))))))) * make_float2 ((F32_cos((_S560))), (F32_sin((_S560))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S559.x))))))) * make_float2 ((F32_cos((_S560))), (F32_sin((_S560))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S561.x))))))) * make_float2 ((F32_cos((_S562))), (F32_sin((_S562))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_9, float dOut_13)
{
    float _S563 = _slang_select(((*dpx_9).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_9).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S563;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S564 = _slang_select(((*dpx_10).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_10).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S564;
    return;
}

inline __device__ float3  abs_0(float3  x_23)
{
    float3  result_21;
    int i_11 = int(0);
    for(;;)
    {
        if(i_11 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_21, i_11) = (F32_abs((_slang_vector_get_element(x_23, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_21;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_15)
{
    DiffPair_float_0 _S565 = *dpx_11;
    bool _S566;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S566 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S566 = false;
    }
    float _S567;
    if(_S566)
    {
        _S567 = dOut_15;
    }
    else
    {
        _S567 = 0.0f;
    }
    dpx_11->primal_0 = _S565.primal_0;
    dpx_11->differential_0 = _S567;
    DiffPair_float_0 _S568 = *dpMin_0;
    if((_S565.primal_0) < ((*dpMin_0).primal_0))
    {
        _S567 = dOut_15;
    }
    else
    {
        _S567 = 0.0f;
    }
    dpMin_0->primal_0 = _S568.primal_0;
    dpMin_0->differential_0 = _S567;
    DiffPair_float_0 _S569 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S567 = dOut_15;
    }
    else
    {
        _S567 = 0.0f;
    }
    dpMax_0->primal_0 = _S569.primal_0;
    dpMax_0->differential_0 = _S567;
    return;
}

inline __device__ float clamp_0(float x_24, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_24), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_12, float dOut_16)
{
    float _S570 = -0.5f / ((*dpx_12).primal_0 * (F32_sqrt(((*dpx_12).primal_0)))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S570;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_17)
{
    float _S571 = (1.0f - (*dps_0).primal_0) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S571;
    DiffPair_float_0 _S572 = *dpy_3;
    float _S573 = (*dps_0).primal_0 * dOut_17;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S573;
    float _S574 = (_S572.primal_0 - (*dpx_13).primal_0) * dOut_17;
    dps_0->primal_0 = _S572.primal_0;
    dps_0->differential_0 = _S574;
    return;
}

inline __device__ float lerp_0(float x_25, float y_7, float s_4)
{
    return x_25 + (y_7 - x_25) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  weights_0, FixedArray<float, 23>  * _S575)
{
    float3  _S576;
    bool _S577;
    bool _S578;
    FixedArray<float, 23>  losses_1;
    float _S579 = float(mask_0);
    float3  _S580 = ref_rgb_0 - render_rgb_0;
    float3  _S581 = abs_0(_S580);
    losses_1[int(0)] = weights_0[int(0)] * _S579 * ((_S581.x + _S581.y + _S581.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S579 * clamp_0(dot_0(_S580, _S580) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S582 = float(depth_mask_0 & mask_0);
    float _S583 = _S582 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S584 = _S582 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S583;
    losses_1[int(3)] = _S584;
    losses_1[int(4)] = _S583 * _S583;
    losses_1[int(5)] = _S584 * _S584;
    losses_1[int(6)] = _S583 * _S584;
    bool _S585 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S586 = norm2_0 == 0.0f;
        _S577 = _S586;
        if(_S586)
        {
            _S576 = make_float3 (0.0f);
            break;
        }
        _S576 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S587;
    bool _S588 = !_S577;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S589 = norm2_1 == 0.0f;
        _S578 = _S589;
        if(_S589)
        {
            _S587 = make_float3 (0.0f);
            break;
        }
        _S587 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S590;
    float3  _S591;
    bool _S592 = !_S578;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S591 = make_float3 (0.0f);
            _S590 = false;
            break;
        }
        _S591 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S590 = _S585;
        break;
    }
    float _S593 = float(_S588 & _S590);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S576, _S591);
    losses_1[int(7)] = weights_0[int(2)] * _S593 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S594 = float(_S592 & _S590);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S587, _S591);
    losses_1[int(8)] = weights_0[int(2)] * _S594 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S595 = float(_S588 & _S592);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S576, _S587);
    losses_1[int(11)] = weights_0[int(5)] * _S595 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S596 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S597 = float(alpha_mask_0);
    float _S598 = float(ref_alpha_0);
    float _S599 = (F32_max((_S596), (_S598)));
    losses_1[int(9)] = weights_0[int(3)] * _S597 * - lerp_0((F32_log(((F32_max((1.0f - _S599), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S599), (9.99999997475242708e-07f)))))), _S598);
    float _S600 = 1.0f - _S596;
    float _S601 = 1.0f - _S598;
    float _S602 = (F32_max((_S600), (_S601)));
    losses_1[int(10)] = weights_0[int(4)] * _S597 * - lerp_0((F32_log(((F32_max((1.0f - _S602), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S602), (9.99999997475242708e-07f)))))), _S601);
    losses_1[int(12)] = weights_0[int(6)] * 4.0f * _S596 * _S600;
    float _S603 = (F32_max((_S596), (9.999999960041972e-13f)));
    losses_1[int(13)] = weights_0[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S603;
    losses_1[int(14)] = weights_0[int(8)] * depth_dist_0 / _S603;
    losses_1[int(15)] = weights_0[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S603;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S579;
    losses_1[int(18)] = _S582;
    losses_1[int(19)] = _S593;
    losses_1[int(20)] = _S594;
    losses_1[int(21)] = _S595;
    losses_1[int(22)] = _S597;
    *_S575 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S604)
{
    return (F32_rsqrt((_S604)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S605, float _S606, float _S607)
{
    return clamp_0(_S605, _S606, _S607);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S608, DiffPair_float_0 * _S609, DiffPair_float_0 * _S610, float _S611)
{
    _d_lerp_0(_S608, _S609, _S610, _S611);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S612, DiffPair_float_0 * _S613, DiffPair_float_0 * _S614, float _S615)
{
    _d_clamp_0(_S612, _S613, _S614, _S615);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S616, float _S617)
{
    _d_rsqrt_0(_S616, _S617);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S618, float3  _S619)
{
    _d_abs_vector_0(_S618, _S619);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_6)
{
    DiffPair_float_0 _S620 = *dprender_depth_0;
    DiffPair_float_0 _S621 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S622 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S623 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S624 = *dpref_normal_0;
    DiffPair_float_0 _S625 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S626 = *dprgb_dist_0;
    DiffPair_float_0 _S627 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S628 = *dpnormal_dist_0;
    float3  _S629 = make_float3 (0.0f);
    float _S630 = float(mask_1);
    float _S631 = (*weights_1)[int(0)] * _S630;
    float3  _S632 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S633 = s_primal_ctx_dot_0(_S632, _S632) * 0.3333333432674408f;
    float _S634 = float(depth_mask_1 & mask_1);
    float _S635 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S636 = _S634 * s_primal_ctx_log_0(_S635);
    float _S637 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S638 = _S634 * s_primal_ctx_log_0(_S637);
    bool _S639 = normal_mask_1 & mask_1;
    float _S640 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S641 = _S640 == 0.0f;
    float3  _S642;
    if(_S641)
    {
        _S642 = make_float3 (0.0f);
    }
    bool _S643 = !_S641;
    float3  _S644;
    if(_S643)
    {
        float _S645 = s_primal_ctx_rsqrt_0(_S640);
        float3  _S646 = make_float3 (_S645);
        _S642 = _S622.primal_0 * make_float3 (_S645);
        _S644 = _S646;
    }
    else
    {
        _S644 = _S629;
    }
    float _S647 = s_primal_ctx_dot_0(_S623.primal_0, _S623.primal_0);
    bool _S648 = _S647 == 0.0f;
    float3  _S649;
    if(_S648)
    {
        _S649 = make_float3 (0.0f);
    }
    bool _S650 = !_S648;
    float3  _S651;
    if(_S650)
    {
        float _S652 = s_primal_ctx_rsqrt_0(_S647);
        float3  _S653 = make_float3 (_S652);
        _S649 = _S623.primal_0 * make_float3 (_S652);
        _S651 = _S653;
    }
    else
    {
        _S651 = _S629;
    }
    float _S654 = s_primal_ctx_dot_0(_S624.primal_0, _S624.primal_0);
    bool _S655 = _S654 == 0.0f;
    float3  _S656;
    bool _S657;
    if(_S655)
    {
        float3  _S658 = make_float3 (0.0f);
        _S657 = false;
        _S656 = _S658;
    }
    else
    {
        _S657 = _S639;
    }
    bool _S659 = !_S655;
    float3  _S660;
    if(_S659)
    {
        float _S661 = s_primal_ctx_rsqrt_0(_S654);
        float3  _S662 = make_float3 (_S661);
        _S656 = _S624.primal_0 * make_float3 (_S661);
        _S660 = _S662;
    }
    else
    {
        _S660 = _S629;
    }
    float _S663 = (*weights_1)[int(2)] * float(_S643 & _S657);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S642, _S656);
    float _S664 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S665 = (*weights_1)[int(2)] * float(_S650 & _S657);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S649, _S656);
    float _S666 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S667 = (*weights_1)[int(5)] * float(_S643 & _S650);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S642, _S649);
    float _S668 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S669 = s_primal_ctx_clamp_0(_S625.primal_0, 0.0f, 1.0f);
    float _S670 = float(alpha_mask_1);
    float _S671 = (*weights_1)[int(3)] * _S670;
    float _S672 = float(ref_alpha_1);
    float _S673 = (F32_max((_S669), (_S672)));
    float _S674 = 1.0f - _S673;
    float _S675 = (F32_max((_S674), (9.99999997475242708e-07f)));
    float _S676 = s_primal_ctx_log_0(_S675);
    float _S677 = (F32_max((_S673), (9.99999997475242708e-07f)));
    float _S678 = s_primal_ctx_log_0(_S677);
    float _S679 = (*weights_1)[int(4)] * _S670;
    float _S680 = 1.0f - _S669;
    float _S681 = 1.0f - _S672;
    float _S682 = (F32_max((_S680), (_S681)));
    float _S683 = 1.0f - _S682;
    float _S684 = (F32_max((_S683), (9.99999997475242708e-07f)));
    float _S685 = s_primal_ctx_log_0(_S684);
    float _S686 = (F32_max((_S682), (9.99999997475242708e-07f)));
    float _S687 = s_primal_ctx_log_0(_S686);
    float _S688 = (*weights_1)[int(6)] * 4.0f;
    float _S689 = _S688 * _S669;
    float _S690 = (F32_max((_S669), (9.999999960041972e-13f)));
    float _S691 = _S690 * _S690;
    float _S692 = (*_s_dOut_6)[int(0)];
    float _S693 = (*_s_dOut_6)[int(1)];
    float _S694 = (*_s_dOut_6)[int(2)];
    float _S695 = (*_s_dOut_6)[int(3)];
    float _S696 = (*_s_dOut_6)[int(4)];
    float _S697 = (*_s_dOut_6)[int(5)];
    float _S698 = (*_s_dOut_6)[int(6)];
    float _S699 = (*_s_dOut_6)[int(15)] / _S691;
    float _S700 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S690 * _S699));
    float _S701 = (*_s_dOut_6)[int(14)] / _S691;
    float _S702 = (*weights_1)[int(8)] * (_S690 * _S701);
    float _S703 = (*_s_dOut_6)[int(13)] / _S691;
    float _S704 = _S690 * _S703;
    float _S705 = (*weights_1)[int(9)] * ((_S628.primal_0.x + _S628.primal_0.y + _S628.primal_0.z) * 0.3333333432674408f) * - _S699 + (*weights_1)[int(8)] * _S627.primal_0 * - _S701 + (*weights_1)[int(7)] * ((_S626.primal_0.x + _S626.primal_0.y + _S626.primal_0.z) * 0.3333333432674408f) * - _S703;
    DiffPair_float_0 _S706;
    (&_S706)->primal_0 = _S669;
    (&_S706)->differential_0 = 0.0f;
    DiffPair_float_0 _S707;
    (&_S707)->primal_0 = 9.999999960041972e-13f;
    (&_S707)->differential_0 = 0.0f;
    _d_max_0(&_S706, &_S707, _S705);
    float _S708 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S704);
    float _S709 = _S689 * (*_s_dOut_6)[int(12)];
    float _S710 = _S688 * (_S680 * (*_s_dOut_6)[int(12)]);
    float _S711 = - (_S679 * (*_s_dOut_6)[int(10)]);
    DiffPair_float_0 _S712;
    (&_S712)->primal_0 = _S685;
    (&_S712)->differential_0 = 0.0f;
    DiffPair_float_0 _S713;
    (&_S713)->primal_0 = _S687;
    (&_S713)->differential_0 = 0.0f;
    DiffPair_float_0 _S714;
    (&_S714)->primal_0 = _S681;
    (&_S714)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S712, &_S713, &_S714, _S711);
    DiffPair_float_0 _S715;
    (&_S715)->primal_0 = _S686;
    (&_S715)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S715, _S713.differential_0);
    DiffPair_float_0 _S716;
    (&_S716)->primal_0 = _S682;
    (&_S716)->differential_0 = 0.0f;
    DiffPair_float_0 _S717;
    (&_S717)->primal_0 = 9.99999997475242708e-07f;
    (&_S717)->differential_0 = 0.0f;
    _d_max_0(&_S716, &_S717, _S715.differential_0);
    DiffPair_float_0 _S718;
    (&_S718)->primal_0 = _S684;
    (&_S718)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S718, _S712.differential_0);
    DiffPair_float_0 _S719;
    (&_S719)->primal_0 = _S683;
    (&_S719)->differential_0 = 0.0f;
    DiffPair_float_0 _S720;
    (&_S720)->primal_0 = 9.99999997475242708e-07f;
    (&_S720)->differential_0 = 0.0f;
    _d_max_0(&_S719, &_S720, _S718.differential_0);
    float _S721 = _S716.differential_0 + - _S719.differential_0;
    DiffPair_float_0 _S722;
    (&_S722)->primal_0 = _S680;
    (&_S722)->differential_0 = 0.0f;
    DiffPair_float_0 _S723;
    (&_S723)->primal_0 = _S681;
    (&_S723)->differential_0 = 0.0f;
    _d_max_0(&_S722, &_S723, _S721);
    float _S724 = - (_S709 + _S722.differential_0);
    float _S725 = - (_S671 * (*_s_dOut_6)[int(9)]);
    DiffPair_float_0 _S726;
    (&_S726)->primal_0 = _S676;
    (&_S726)->differential_0 = 0.0f;
    DiffPair_float_0 _S727;
    (&_S727)->primal_0 = _S678;
    (&_S727)->differential_0 = 0.0f;
    DiffPair_float_0 _S728;
    (&_S728)->primal_0 = _S672;
    (&_S728)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S726, &_S727, &_S728, _S725);
    DiffPair_float_0 _S729;
    (&_S729)->primal_0 = _S677;
    (&_S729)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S729, _S727.differential_0);
    DiffPair_float_0 _S730;
    (&_S730)->primal_0 = _S673;
    (&_S730)->differential_0 = 0.0f;
    DiffPair_float_0 _S731;
    (&_S731)->primal_0 = 9.99999997475242708e-07f;
    (&_S731)->differential_0 = 0.0f;
    _d_max_0(&_S730, &_S731, _S729.differential_0);
    DiffPair_float_0 _S732;
    (&_S732)->primal_0 = _S675;
    (&_S732)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S732, _S726.differential_0);
    DiffPair_float_0 _S733;
    (&_S733)->primal_0 = _S674;
    (&_S733)->differential_0 = 0.0f;
    DiffPair_float_0 _S734;
    (&_S734)->primal_0 = 9.99999997475242708e-07f;
    (&_S734)->differential_0 = 0.0f;
    _d_max_0(&_S733, &_S734, _S732.differential_0);
    float _S735 = _S730.differential_0 + - _S733.differential_0;
    DiffPair_float_0 _S736;
    (&_S736)->primal_0 = _S669;
    (&_S736)->differential_0 = 0.0f;
    DiffPair_float_0 _S737;
    (&_S737)->primal_0 = _S672;
    (&_S737)->differential_0 = 0.0f;
    _d_max_0(&_S736, &_S737, _S735);
    float _S738 = _S706.differential_0 + _S710 + _S724 + _S736.differential_0;
    DiffPair_float_0 _S739;
    (&_S739)->primal_0 = _S625.primal_0;
    (&_S739)->differential_0 = 0.0f;
    DiffPair_float_0 _S740;
    (&_S740)->primal_0 = 0.0f;
    (&_S740)->differential_0 = 0.0f;
    DiffPair_float_0 _S741;
    (&_S741)->primal_0 = 1.0f;
    (&_S741)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S739, &_S740, &_S741, _S738);
    DiffPair_float_0 _S742 = _S739;
    float _S743 = _S667 * (*_s_dOut_6)[int(11)];
    DiffPair_float_0 _S744;
    (&_S744)->primal_0 = _S668;
    (&_S744)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S744, _S743);
    DiffPair_float_0 _S745;
    (&_S745)->primal_0 = cos_sim_loss_5;
    (&_S745)->differential_0 = 0.0f;
    DiffPair_float_0 _S746;
    (&_S746)->primal_0 = 9.999999960041972e-13f;
    (&_S746)->differential_0 = 0.0f;
    _d_max_0(&_S745, &_S746, _S744.differential_0);
    float _S747 = 0.5f * - (_S743 + _S745.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S748;
    (&_S748)->primal_0 = _S642;
    (&_S748)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S749;
    (&_S749)->primal_0 = _S649;
    (&_S749)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S748, &_S749, _S747);
    float _S750 = _S665 * (*_s_dOut_6)[int(8)];
    DiffPair_float_0 _S751;
    (&_S751)->primal_0 = _S666;
    (&_S751)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S751, _S750);
    DiffPair_float_0 _S752;
    (&_S752)->primal_0 = cos_sim_loss_4;
    (&_S752)->differential_0 = 0.0f;
    DiffPair_float_0 _S753;
    (&_S753)->primal_0 = 9.999999960041972e-13f;
    (&_S753)->differential_0 = 0.0f;
    _d_max_0(&_S752, &_S753, _S751.differential_0);
    float _S754 = 0.5f * - (_S750 + _S752.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S755;
    (&_S755)->primal_0 = _S649;
    (&_S755)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S756;
    (&_S756)->primal_0 = _S656;
    (&_S756)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S755, &_S756, _S754);
    float _S757 = _S663 * (*_s_dOut_6)[int(7)];
    DiffPair_float_0 _S758;
    (&_S758)->primal_0 = _S664;
    (&_S758)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S758, _S757);
    DiffPair_float_0 _S759;
    (&_S759)->primal_0 = cos_sim_loss_3;
    (&_S759)->differential_0 = 0.0f;
    DiffPair_float_0 _S760;
    (&_S760)->primal_0 = 9.999999960041972e-13f;
    (&_S760)->differential_0 = 0.0f;
    _d_max_0(&_S759, &_S760, _S758.differential_0);
    float _S761 = 0.5f * - (_S757 + _S759.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S762;
    (&_S762)->primal_0 = _S642;
    (&_S762)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S763;
    (&_S763)->primal_0 = _S656;
    (&_S763)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S762, &_S763, _S761);
    float3  _S764 = _S756.differential_0 + _S763.differential_0;
    float3  _S765 = _S748.differential_0 + _S762.differential_0;
    float3  _S766 = make_float3 (_S700, _S700, _S700);
    float3  _S767 = make_float3 (_S708, _S708, _S708);
    float3  _S768 = _S749.differential_0 + _S755.differential_0;
    float _S769;
    if(_S659)
    {
        float3  _S770 = _S624.primal_0 * _S764;
        float3  _S771 = _S660 * _S764;
        float _S772 = _S770.x + _S770.y + _S770.z;
        DiffPair_float_0 _S773;
        (&_S773)->primal_0 = _S654;
        (&_S773)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S773, _S772);
        _S769 = _S773.differential_0;
        _S642 = _S771;
    }
    else
    {
        _S769 = 0.0f;
        _S642 = _S629;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S774;
    (&_S774)->primal_0 = _S624.primal_0;
    (&_S774)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S775;
    (&_S775)->primal_0 = _S624.primal_0;
    (&_S775)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S774, &_S775, _S769);
    float3  _S776 = _S775.differential_0 + _S774.differential_0 + _S642;
    if(_S650)
    {
        float3  _S777 = _S623.primal_0 * _S768;
        float3  _S778 = _S651 * _S768;
        float _S779 = _S777.x + _S777.y + _S777.z;
        DiffPair_float_0 _S780;
        (&_S780)->primal_0 = _S647;
        (&_S780)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S780, _S779);
        _S769 = _S780.differential_0;
        _S642 = _S778;
    }
    else
    {
        _S769 = 0.0f;
        _S642 = _S629;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S781;
    (&_S781)->primal_0 = _S623.primal_0;
    (&_S781)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S782;
    (&_S782)->primal_0 = _S623.primal_0;
    (&_S782)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S781, &_S782, _S769);
    float3  _S783 = _S782.differential_0 + _S781.differential_0 + _S642;
    if(_S643)
    {
        float3  _S784 = _S622.primal_0 * _S765;
        float3  _S785 = _S644 * _S765;
        float _S786 = _S784.x + _S784.y + _S784.z;
        DiffPair_float_0 _S787;
        (&_S787)->primal_0 = _S640;
        (&_S787)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S787, _S786);
        _S769 = _S787.differential_0;
        _S642 = _S785;
    }
    else
    {
        _S769 = 0.0f;
        _S642 = _S629;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S788;
    (&_S788)->primal_0 = _S622.primal_0;
    (&_S788)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S789;
    (&_S789)->primal_0 = _S622.primal_0;
    (&_S789)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S788, &_S789, _S769);
    float _S790 = _S638 * _S698;
    float _S791 = _S638 * _S697;
    float _S792 = _S636 * _S696;
    float _S793 = _S634 * (_S636 * _S698 + _S791 + _S791 + _S695);
    DiffPair_float_0 _S794;
    (&_S794)->primal_0 = _S637;
    (&_S794)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S794, _S793);
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = _S621.primal_0;
    (&_S795)->differential_0 = 0.0f;
    DiffPair_float_0 _S796;
    (&_S796)->primal_0 = 0.00009999999747379f;
    (&_S796)->differential_0 = 0.0f;
    _d_max_0(&_S795, &_S796, _S794.differential_0);
    float _S797 = _S634 * (_S790 + _S792 + _S792 + _S694);
    DiffPair_float_0 _S798;
    (&_S798)->primal_0 = _S635;
    (&_S798)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S798, _S797);
    DiffPair_float_0 _S799;
    (&_S799)->primal_0 = _S620.primal_0;
    (&_S799)->differential_0 = 0.0f;
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = 0.00009999999747379f;
    (&_S800)->differential_0 = 0.0f;
    _d_max_0(&_S799, &_S800, _S798.differential_0);
    float _S801 = _S630 * _S693;
    DiffPair_float_0 _S802;
    (&_S802)->primal_0 = _S633;
    (&_S802)->differential_0 = 0.0f;
    DiffPair_float_0 _S803;
    (&_S803)->primal_0 = 0.0f;
    (&_S803)->differential_0 = 0.0f;
    DiffPair_float_0 _S804;
    (&_S804)->primal_0 = 1.0f;
    (&_S804)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S802, &_S803, &_S804, _S801);
    float _S805 = 0.3333333432674408f * _S802.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S806;
    (&_S806)->primal_0 = _S632;
    (&_S806)->differential_0 = _S629;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S807;
    (&_S807)->primal_0 = _S632;
    (&_S807)->differential_0 = _S629;
    s_bwd_prop_dot_0(&_S806, &_S807, _S805);
    float _S808 = 0.3333333432674408f * (_S631 * _S692);
    float3  _S809 = make_float3 (_S808, _S808, _S808);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S810;
    (&_S810)->primal_0 = _S632;
    (&_S810)->differential_0 = _S629;
    s_bwd_prop_abs_0(&_S810, _S809);
    float3  _S811 = _S807.differential_0 + _S806.differential_0 + _S810.differential_0;
    float3  _S812 = - _S811;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S766;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S702;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S767;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S742.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S776;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S783;
    float3  _S813 = _S789.differential_0 + _S788.differential_0 + _S642;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S813;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S795.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S799.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S811;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S812;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S814, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S815, DiffPair_float_0 * _S816, DiffPair_float_0 * _S817, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S818, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S819, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S820, DiffPair_float_0 * _S821, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S822, DiffPair_float_0 * _S823, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S824, bool _S825, bool _S826, bool _S827, bool _S828, bool _S829, FixedArray<float, 10>  * _S830, FixedArray<float, 23>  * _S831)
{
    s_bwd_prop_per_pixel_losses_0(_S814, _S815, _S816, _S817, _S818, _S819, _S820, _S821, _S822, _S823, _S824, _S825, _S826, _S827, _S828, _S829, _S830, _S831);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S832 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S832;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S832;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S832;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S832;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S832;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S832;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S832;
    FixedArray<float, 10>  _S833 = weights_2;
    FixedArray<float, 23>  _S834 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S833, &_S834);
    *v_render_rgb_0 = dp_render_rgb_0.differential_0;
    *v_ref_rgb_0 = dp_ref_rgb_0.differential_0;
    *v_render_depth_0 = dp_render_depth_0.differential_0;
    *v_ref_depth_0 = dp_ref_depth_0.differential_0;
    *v_render_normal_0 = dp_render_normal_0.differential_0;
    *v_depth_normal_0 = dp_depth_normal_0.differential_0;
    *v_ref_normal_0 = dp_ref_normal_0.differential_0;
    *v_render_alpha_0 = dp_render_alpha_0.differential_0;
    *v_rgb_dist_0 = dp_rgb_dist_0.differential_0;
    *v_depth_dist_0 = dp_depth_dist_0.differential_0;
    *v_normal_dist_0 = dp_normal_dist_0.differential_0;
    return;
}

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_14, float dOut_18)
{
    float _S835 = 1.0f / ((*dpx_14).primal_0 * 2.30258512496948242f) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S835;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 10>  weights_3, FixedArray<float, 10>  * _S836)
{
    FixedArray<float, 10>  losses_2;
    float _S837 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_2[int(0)] = raw_losses_0[int(0)] / _S837;
    losses_2[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S837)));
    bool _S838;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S838 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S838 = false;
    }
    float _S839;
    if(_S838)
    {
        _S839 = weights_3[int(1)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S839 = 0.0f;
    }
    losses_2[int(2)] = _S839;
    losses_2[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_2[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S840 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_2[int(6)] = raw_losses_0[int(12)] / _S840;
    losses_2[int(7)] = raw_losses_0[int(13)] / _S840;
    losses_2[int(8)] = raw_losses_0[int(14)] / _S840;
    losses_2[int(9)] = raw_losses_0[int(15)] / _S840;
    *_S836 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S841)
{
    return (F32_sqrt((_S841)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S842, float _S843)
{
    _d_log10_0(_S842, _S843);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_7)
{
    FixedArray<float, 23>  _S844 = dpraw_losses_0->primal_0;
    float _S845 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S846 = _S845 * _S845;
    float _S847 = dpraw_losses_0->primal_0[int(1)] / _S845;
    bool _S848 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S849;
    if(_S848)
    {
        _S849 = (_S844[int(3)]) != 0.0f;
    }
    else
    {
        _S849 = false;
    }
    float _S850;
    float _S851;
    float _S852;
    float _S853;
    float _S854;
    float _S855;
    float _S856;
    float _S857;
    float _S858;
    float _S859;
    float _S860;
    float _S861;
    float _S862;
    float _S863;
    float _S864;
    if(_S849)
    {
        float _S865 = _S844[int(2)] * _S844[int(3)];
        float _S866 = _S844[int(18)] * _S844[int(18)];
        float _S867 = _S844[int(6)] - _S865 / _S844[int(18)];
        float _S868 = _S844[int(2)] * _S844[int(2)];
        float _S869 = _S844[int(4)] - _S868 / _S844[int(18)];
        float _S870 = _S844[int(3)] * _S844[int(3)];
        float _S871 = _S844[int(5)] - _S870 / _S844[int(18)];
        float _S872 = _S869 * _S871 + 1.0f;
        float _S873 = (F32_max((9.999999960041972e-13f), (_S872)));
        float _S874 = s_primal_ctx_sqrt_0(_S873);
        float _S875 = _S874 * _S874;
        float _S876 = 1.0f - _S867 / _S874;
        _S850 = (*weights_4)[int(1)];
        _S851 = _S876;
        _S852 = _S875;
        _S853 = _S867;
        _S854 = _S874;
        _S855 = _S873;
        _S856 = _S872;
        _S857 = _S869;
        _S858 = _S871;
        _S859 = _S866;
        _S860 = _S870;
        _S861 = _S844[int(3)];
        _S862 = _S868;
        _S863 = _S844[int(2)];
        _S864 = _S865;
    }
    else
    {
        _S850 = 0.0f;
        _S851 = 0.0f;
        _S852 = 0.0f;
        _S853 = 0.0f;
        _S854 = 0.0f;
        _S855 = 0.0f;
        _S856 = 0.0f;
        _S857 = 0.0f;
        _S858 = 0.0f;
        _S859 = 0.0f;
        _S860 = 0.0f;
        _S861 = 0.0f;
        _S862 = 0.0f;
        _S863 = 0.0f;
        _S864 = 0.0f;
    }
    float _S877 = (F32_max((_S844[int(19)]), (1.0f)));
    float _S878 = _S877 * _S877;
    float _S879 = (F32_max((_S844[int(20)]), (1.0f)));
    float _S880 = _S879 * _S879;
    float _S881 = float((I32_max((int((_S844[int(19)]) > 0.5f) + int((_S844[int(20)]) > 0.5f)), (int(1)))));
    float _S882 = _S844[int(9)] + _S844[int(10)];
    float _S883 = (F32_max((_S844[int(22)]), (1.0f)));
    float _S884 = _S883 * _S883;
    float _S885 = (F32_max((_S844[int(21)]), (1.0f)));
    float _S886 = _S885 * _S885;
    float _S887 = (F32_max((_S844[int(16)]), (1.0f)));
    float _S888 = _S887 * _S887;
    float _S889 = (*_s_dOut_7)[int(0)];
    float _S890 = (*_s_dOut_7)[int(1)];
    float _S891 = (*_s_dOut_7)[int(2)];
    float _S892 = (*_s_dOut_7)[int(9)] / _S888;
    float _S893 = _S887 * _S892;
    float _S894 = (*_s_dOut_7)[int(8)] / _S888;
    float _S895 = _S887 * _S894;
    float _S896 = (*_s_dOut_7)[int(7)] / _S888;
    float _S897 = _S887 * _S896;
    float _S898 = (*_s_dOut_7)[int(6)] / _S888;
    float _S899 = _S887 * _S898;
    float _S900 = _S844[int(15)] * - _S892 + _S844[int(14)] * - _S894 + _S844[int(13)] * - _S896 + _S844[int(12)] * - _S898;
    DiffPair_float_0 _S901;
    (&_S901)->primal_0 = _S844[int(16)];
    (&_S901)->differential_0 = 0.0f;
    DiffPair_float_0 _S902;
    (&_S902)->primal_0 = 1.0f;
    (&_S902)->differential_0 = 0.0f;
    _d_max_0(&_S901, &_S902, _S900);
    float _S903 = (*_s_dOut_7)[int(5)] / _S886;
    float _S904 = _S844[int(11)] * - _S903;
    float _S905 = _S885 * _S903;
    DiffPair_float_0 _S906;
    (&_S906)->primal_0 = _S844[int(21)];
    (&_S906)->differential_0 = 0.0f;
    DiffPair_float_0 _S907;
    (&_S907)->primal_0 = 1.0f;
    (&_S907)->differential_0 = 0.0f;
    _d_max_0(&_S906, &_S907, _S904);
    float _S908 = (*_s_dOut_7)[int(4)] / _S884;
    float _S909 = _S882 * - _S908;
    float _S910 = _S883 * _S908;
    DiffPair_float_0 _S911;
    (&_S911)->primal_0 = _S844[int(22)];
    (&_S911)->differential_0 = 0.0f;
    DiffPair_float_0 _S912;
    (&_S912)->primal_0 = 1.0f;
    (&_S912)->differential_0 = 0.0f;
    _d_max_0(&_S911, &_S912, _S909);
    float _S913 = (*_s_dOut_7)[int(3)] / _S881;
    float _S914 = _S913 / _S880;
    float _S915 = _S844[int(8)] * - _S914;
    float _S916 = _S879 * _S914;
    DiffPair_float_0 _S917;
    (&_S917)->primal_0 = _S844[int(20)];
    (&_S917)->differential_0 = 0.0f;
    DiffPair_float_0 _S918;
    (&_S918)->primal_0 = 1.0f;
    (&_S918)->differential_0 = 0.0f;
    _d_max_0(&_S917, &_S918, _S915);
    float _S919 = _S913 / _S878;
    float _S920 = _S844[int(7)] * - _S919;
    float _S921 = _S877 * _S919;
    DiffPair_float_0 _S922;
    (&_S922)->primal_0 = _S844[int(19)];
    (&_S922)->differential_0 = 0.0f;
    DiffPair_float_0 _S923;
    (&_S923)->primal_0 = 1.0f;
    (&_S923)->differential_0 = 0.0f;
    _d_max_0(&_S922, &_S923, _S920);
    FixedArray<float, 23>  _S924;
    _S924[int(0)] = 0.0f;
    _S924[int(1)] = 0.0f;
    _S924[int(2)] = 0.0f;
    _S924[int(3)] = 0.0f;
    _S924[int(4)] = 0.0f;
    _S924[int(5)] = 0.0f;
    _S924[int(6)] = 0.0f;
    _S924[int(7)] = 0.0f;
    _S924[int(8)] = 0.0f;
    _S924[int(9)] = 0.0f;
    _S924[int(10)] = 0.0f;
    _S924[int(11)] = 0.0f;
    _S924[int(12)] = 0.0f;
    _S924[int(13)] = 0.0f;
    _S924[int(14)] = 0.0f;
    _S924[int(15)] = 0.0f;
    _S924[int(16)] = 0.0f;
    _S924[int(17)] = 0.0f;
    _S924[int(18)] = 0.0f;
    _S924[int(19)] = 0.0f;
    _S924[int(20)] = 0.0f;
    _S924[int(21)] = 0.0f;
    _S924[int(22)] = 0.0f;
    _S924[int(15)] = _S893;
    _S924[int(14)] = _S895;
    _S924[int(13)] = _S897;
    _S924[int(16)] = _S901.differential_0;
    _S924[int(12)] = _S899;
    _S924[int(21)] = _S906.differential_0;
    _S924[int(11)] = _S905;
    _S924[int(22)] = _S911.differential_0;
    _S924[int(10)] = _S910;
    _S924[int(9)] = _S910;
    _S924[int(20)] = _S917.differential_0;
    _S924[int(8)] = _S916;
    _S924[int(19)] = _S922.differential_0;
    _S924[int(7)] = _S921;
    float _S925 = _S924[int(0)];
    float _S926 = _S924[int(1)];
    float _S927 = _S924[int(2)];
    float _S928 = _S924[int(3)];
    float _S929 = _S924[int(4)];
    float _S930 = _S924[int(5)];
    float _S931 = _S924[int(6)];
    float _S932 = _S924[int(7)];
    float _S933 = _S924[int(8)];
    float _S934 = _S924[int(9)];
    float _S935 = _S924[int(10)];
    float _S936 = _S924[int(11)];
    float _S937 = _S924[int(12)];
    float _S938 = _S924[int(13)];
    float _S939 = _S924[int(14)];
    float _S940 = _S924[int(15)];
    float _S941 = _S924[int(16)];
    float _S942 = _S924[int(17)];
    float _S943 = _S924[int(18)];
    float _S944 = _S924[int(19)];
    float _S945 = _S924[int(20)];
    float _S946 = _S924[int(21)];
    float _S947 = _S924[int(22)];
    FixedArray<float, 23>  _S948;
    if(_S849)
    {
        float _S949 = _S850 * _S891;
        DiffPair_float_0 _S950;
        (&_S950)->primal_0 = _S851;
        (&_S950)->differential_0 = 0.0f;
        DiffPair_float_0 _S951;
        (&_S951)->primal_0 = 0.0f;
        (&_S951)->differential_0 = 0.0f;
        DiffPair_float_0 _S952;
        (&_S952)->primal_0 = 2.0f;
        (&_S952)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S950, &_S951, &_S952, _S949);
        float _S953 = - _S950.differential_0 / _S852;
        float _S954 = _S853 * - _S953;
        float _S955 = _S854 * _S953;
        DiffPair_float_0 _S956;
        (&_S956)->primal_0 = _S855;
        (&_S956)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S956, _S954);
        DiffPair_float_0 _S957;
        (&_S957)->primal_0 = 9.999999960041972e-13f;
        (&_S957)->differential_0 = 0.0f;
        DiffPair_float_0 _S958;
        (&_S958)->primal_0 = _S856;
        (&_S958)->differential_0 = 0.0f;
        _d_max_0(&_S957, &_S958, _S956.differential_0);
        float _S959 = _S857 * _S958.differential_0;
        float _S960 = _S858 * _S958.differential_0;
        float _S961 = - _S959 / _S859;
        float _S962 = _S861 * (_S844[int(18)] * _S961);
        float _S963 = - _S960 / _S859;
        float _S964 = _S863 * (_S844[int(18)] * _S963);
        float _S965 = - _S955 / _S859;
        float _S966 = _S844[int(18)] * _S965;
        float _S967 = _S962 + _S962 + _S863 * _S966;
        float _S968 = _S964 + _S964 + _S861 * _S966;
        float _S969 = _S860 * - _S961 + _S862 * - _S963 + _S864 * - _S965;
        FixedArray<float, 23>  _S970;
        _S970[int(0)] = 0.0f;
        _S970[int(1)] = 0.0f;
        _S970[int(2)] = 0.0f;
        _S970[int(3)] = 0.0f;
        _S970[int(4)] = 0.0f;
        _S970[int(5)] = 0.0f;
        _S970[int(6)] = 0.0f;
        _S970[int(7)] = 0.0f;
        _S970[int(8)] = 0.0f;
        _S970[int(9)] = 0.0f;
        _S970[int(10)] = 0.0f;
        _S970[int(11)] = 0.0f;
        _S970[int(12)] = 0.0f;
        _S970[int(13)] = 0.0f;
        _S970[int(14)] = 0.0f;
        _S970[int(15)] = 0.0f;
        _S970[int(16)] = 0.0f;
        _S970[int(17)] = 0.0f;
        _S970[int(18)] = 0.0f;
        _S970[int(19)] = 0.0f;
        _S970[int(20)] = 0.0f;
        _S970[int(21)] = 0.0f;
        _S970[int(22)] = 0.0f;
        _S970[int(5)] = _S959;
        _S970[int(4)] = _S960;
        _S970[int(3)] = _S967;
        _S970[int(2)] = _S968;
        _S970[int(6)] = _S955;
        float _S971 = _S926 + _S970[int(1)];
        float _S972 = _S927 + _S970[int(2)];
        float _S973 = _S928 + _S970[int(3)];
        float _S974 = _S929 + _S970[int(4)];
        float _S975 = _S930 + _S970[int(5)];
        float _S976 = _S931 + _S970[int(6)];
        float _S977 = _S932 + _S970[int(7)];
        float _S978 = _S933 + _S970[int(8)];
        float _S979 = _S934 + _S970[int(9)];
        float _S980 = _S935 + _S970[int(10)];
        float _S981 = _S936 + _S970[int(11)];
        float _S982 = _S937 + _S970[int(12)];
        float _S983 = _S938 + _S970[int(13)];
        float _S984 = _S939 + _S970[int(14)];
        float _S985 = _S940 + _S970[int(15)];
        float _S986 = _S941 + _S970[int(16)];
        float _S987 = _S942 + _S970[int(17)];
        float _S988 = _S943 + _S970[int(18)];
        float _S989 = _S944 + _S970[int(19)];
        float _S990 = _S945 + _S970[int(20)];
        float _S991 = _S946 + _S970[int(21)];
        float _S992 = _S947 + _S970[int(22)];
        _S948[int(0)] = _S925 + _S970[int(0)];
        _S948[int(1)] = _S971;
        _S948[int(2)] = _S972;
        _S948[int(3)] = _S973;
        _S948[int(4)] = _S974;
        _S948[int(5)] = _S975;
        _S948[int(6)] = _S976;
        _S948[int(7)] = _S977;
        _S948[int(8)] = _S978;
        _S948[int(9)] = _S979;
        _S948[int(10)] = _S980;
        _S948[int(11)] = _S981;
        _S948[int(12)] = _S982;
        _S948[int(13)] = _S983;
        _S948[int(14)] = _S984;
        _S948[int(15)] = _S985;
        _S948[int(16)] = _S986;
        _S948[int(17)] = _S987;
        _S948[int(18)] = _S988;
        _S948[int(19)] = _S989;
        _S948[int(20)] = _S990;
        _S948[int(21)] = _S991;
        _S948[int(22)] = _S992;
        _S850 = _S969;
    }
    else
    {
        _S948[int(0)] = _S925;
        _S948[int(1)] = _S926;
        _S948[int(2)] = _S927;
        _S948[int(3)] = _S928;
        _S948[int(4)] = _S929;
        _S948[int(5)] = _S930;
        _S948[int(6)] = _S931;
        _S948[int(7)] = _S932;
        _S948[int(8)] = _S933;
        _S948[int(9)] = _S934;
        _S948[int(10)] = _S935;
        _S948[int(11)] = _S936;
        _S948[int(12)] = _S937;
        _S948[int(13)] = _S938;
        _S948[int(14)] = _S939;
        _S948[int(15)] = _S940;
        _S948[int(16)] = _S941;
        _S948[int(17)] = _S942;
        _S948[int(18)] = _S943;
        _S948[int(19)] = _S944;
        _S948[int(20)] = _S945;
        _S948[int(21)] = _S946;
        _S948[int(22)] = _S947;
        _S850 = 0.0f;
    }
    if(_S848)
    {
        FixedArray<float, 23>  _S993;
        _S993[int(0)] = 0.0f;
        _S993[int(1)] = 0.0f;
        _S993[int(2)] = 0.0f;
        _S993[int(3)] = 0.0f;
        _S993[int(4)] = 0.0f;
        _S993[int(5)] = 0.0f;
        _S993[int(6)] = 0.0f;
        _S993[int(7)] = 0.0f;
        _S993[int(8)] = 0.0f;
        _S993[int(9)] = 0.0f;
        _S993[int(10)] = 0.0f;
        _S993[int(11)] = 0.0f;
        _S993[int(12)] = 0.0f;
        _S993[int(13)] = 0.0f;
        _S993[int(14)] = 0.0f;
        _S993[int(15)] = 0.0f;
        _S993[int(16)] = 0.0f;
        _S993[int(17)] = 0.0f;
        _S993[int(18)] = 0.0f;
        _S993[int(19)] = 0.0f;
        _S993[int(20)] = 0.0f;
        _S993[int(21)] = 0.0f;
        _S993[int(22)] = 0.0f;
        _S993[int(3)] = 0.0f;
        float _S994 = _S948[int(1)] + _S993[int(1)];
        float _S995 = _S948[int(2)] + _S993[int(2)];
        float _S996 = _S948[int(3)] + _S993[int(3)];
        float _S997 = _S948[int(4)] + _S993[int(4)];
        float _S998 = _S948[int(5)] + _S993[int(5)];
        float _S999 = _S948[int(6)] + _S993[int(6)];
        float _S1000 = _S948[int(7)] + _S993[int(7)];
        float _S1001 = _S948[int(8)] + _S993[int(8)];
        float _S1002 = _S948[int(9)] + _S993[int(9)];
        float _S1003 = _S948[int(10)] + _S993[int(10)];
        float _S1004 = _S948[int(11)] + _S993[int(11)];
        float _S1005 = _S948[int(12)] + _S993[int(12)];
        float _S1006 = _S948[int(13)] + _S993[int(13)];
        float _S1007 = _S948[int(14)] + _S993[int(14)];
        float _S1008 = _S948[int(15)] + _S993[int(15)];
        float _S1009 = _S948[int(16)] + _S993[int(16)];
        float _S1010 = _S948[int(17)] + _S993[int(17)];
        float _S1011 = _S948[int(18)] + _S993[int(18)];
        float _S1012 = _S948[int(19)] + _S993[int(19)];
        float _S1013 = _S948[int(20)] + _S993[int(20)];
        float _S1014 = _S948[int(21)] + _S993[int(21)];
        float _S1015 = _S948[int(22)] + _S993[int(22)];
        _S948[int(0)] = _S948[int(0)] + _S993[int(0)];
        _S948[int(1)] = _S994;
        _S948[int(2)] = _S995;
        _S948[int(3)] = _S996;
        _S948[int(4)] = _S997;
        _S948[int(5)] = _S998;
        _S948[int(6)] = _S999;
        _S948[int(7)] = _S1000;
        _S948[int(8)] = _S1001;
        _S948[int(9)] = _S1002;
        _S948[int(10)] = _S1003;
        _S948[int(11)] = _S1004;
        _S948[int(12)] = _S1005;
        _S948[int(13)] = _S1006;
        _S948[int(14)] = _S1007;
        _S948[int(15)] = _S1008;
        _S948[int(16)] = _S1009;
        _S948[int(17)] = _S1010;
        _S948[int(18)] = _S1011;
        _S948[int(19)] = _S1012;
        _S948[int(20)] = _S1013;
        _S948[int(21)] = _S1014;
        _S948[int(22)] = _S1015;
    }
    float _S1016 = -10.0f * _S890;
    DiffPair_float_0 _S1017;
    (&_S1017)->primal_0 = _S847;
    (&_S1017)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S1017, _S1016);
    float _S1018 = _S1017.differential_0 / _S846;
    float _S1019 = _S845 * _S1018;
    float _S1020 = _S889 / _S846;
    float _S1021 = _S845 * _S1020;
    float _S1022 = _S844[int(1)] * - _S1018 + _S844[int(0)] * - _S1020;
    DiffPair_float_0 _S1023;
    (&_S1023)->primal_0 = _S844[int(17)];
    (&_S1023)->differential_0 = 0.0f;
    DiffPair_float_0 _S1024;
    (&_S1024)->primal_0 = 1.0f;
    (&_S1024)->differential_0 = 0.0f;
    _d_max_0(&_S1023, &_S1024, _S1022);
    FixedArray<float, 23>  _S1025;
    _S1025[int(0)] = 0.0f;
    _S1025[int(1)] = 0.0f;
    _S1025[int(2)] = 0.0f;
    _S1025[int(3)] = 0.0f;
    _S1025[int(4)] = 0.0f;
    _S1025[int(5)] = 0.0f;
    _S1025[int(6)] = 0.0f;
    _S1025[int(7)] = 0.0f;
    _S1025[int(8)] = 0.0f;
    _S1025[int(9)] = 0.0f;
    _S1025[int(10)] = 0.0f;
    _S1025[int(11)] = 0.0f;
    _S1025[int(12)] = 0.0f;
    _S1025[int(13)] = 0.0f;
    _S1025[int(14)] = 0.0f;
    _S1025[int(15)] = 0.0f;
    _S1025[int(16)] = 0.0f;
    _S1025[int(17)] = 0.0f;
    _S1025[int(18)] = 0.0f;
    _S1025[int(19)] = 0.0f;
    _S1025[int(20)] = 0.0f;
    _S1025[int(21)] = 0.0f;
    _S1025[int(22)] = 0.0f;
    _S1025[int(18)] = _S850;
    _S1025[int(1)] = _S1019;
    _S1025[int(17)] = _S1023.differential_0;
    _S1025[int(0)] = _S1021;
    FixedArray<float, 23>  _S1026 = {
        _S948[int(0)] + _S1025[int(0)], _S948[int(1)] + _S1025[int(1)], _S948[int(2)] + _S1025[int(2)], _S948[int(3)] + _S1025[int(3)], _S948[int(4)] + _S1025[int(4)], _S948[int(5)] + _S1025[int(5)], _S948[int(6)] + _S1025[int(6)], _S948[int(7)] + _S1025[int(7)], _S948[int(8)] + _S1025[int(8)], _S948[int(9)] + _S1025[int(9)], _S948[int(10)] + _S1025[int(10)], _S948[int(11)] + _S1025[int(11)], _S948[int(12)] + _S1025[int(12)], _S948[int(13)] + _S1025[int(13)], _S948[int(14)] + _S1025[int(14)], _S948[int(15)] + _S1025[int(15)], _S948[int(16)] + _S1025[int(16)], _S948[int(17)] + _S1025[int(17)], _S948[int(18)] + _S1025[int(18)], _S948[int(19)] + _S1025[int(19)], _S948[int(20)] + _S1025[int(20)], _S948[int(21)] + _S1025[int(21)], _S948[int(22)] + _S1025[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S1026;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S1027, FixedArray<float, 10>  * _S1028, FixedArray<float, 10>  * _S1029)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S1027, _S1028, _S1029);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 10>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S1030)
{
    FixedArray<float, 23>  _S1031 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S1031;
    FixedArray<float, 10>  _S1032 = weights_5;
    FixedArray<float, 10>  _S1033 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S1032, &_S1033);
    *_S1030 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_26, float3  y_8)
{
    float3  result_22;
    int i_12 = int(0);
    for(;;)
    {
        if(i_12 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_22, i_12) = (F32_min((_slang_vector_get_element(x_26, i_12)), (_slang_vector_get_element(y_8, i_12))));
        i_12 = i_12 + int(1);
    }
    return result_22;
}

inline __device__ float3  max_0(float3  x_27, float3  y_9)
{
    float3  result_23;
    int i_13 = int(0);
    for(;;)
    {
        if(i_13 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_23, i_13) = (F32_max((_slang_vector_get_element(x_27, i_13)), (_slang_vector_get_element(y_9, i_13))));
        i_13 = i_13 + int(1);
    }
    return result_23;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_19)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_15).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_4).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_19.x);
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_4;
    *&((&right_d_result_4)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_15).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_4).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_19.y);
    *&((&left_d_result_4)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_4)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_15).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_4).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_19.z);
    *&((&left_d_result_4)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_4)->z) = right_dp_2.differential_0;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = left_d_result_4;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  clamp_1(float3  x_28, float3  minBound_1, float3  maxBound_1)
{
    return min_0(max_0(x_28, minBound_1), maxBound_1);
}

inline __device__ float3  blend_background(float3  rgb_2, float alpha_0, float3  background_0)
{
    return clamp_1(rgb_2 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1034, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1035, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1036, float3  _S1037)
{
    _d_clamp_vector_0(_S1034, _S1035, _S1036, _S1037);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_8)
{
    float _S1038 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S1039 = make_float3 (_S1038);
    float3  _S1040 = make_float3 (0.0f);
    float3  _S1041 = make_float3 (1.0f);
    float3  _S1042 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1043;
    (&_S1043)->primal_0 = (*dprgb_1).primal_0 + make_float3 (_S1038) * (*dpbackground_0).primal_0;
    (&_S1043)->differential_0 = _S1042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1044;
    (&_S1044)->primal_0 = _S1040;
    (&_S1044)->differential_0 = _S1042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1045;
    (&_S1045)->primal_0 = _S1041;
    (&_S1045)->differential_0 = _S1042;
    s_bwd_prop_clamp_1(&_S1043, &_S1044, &_S1045, _s_dOut_8);
    float3  _S1046 = _S1039 * _S1043.differential_0;
    float3  _S1047 = (*dpbackground_0).primal_0 * _S1043.differential_0;
    float _S1048 = - (_S1047.x + _S1047.y + _S1047.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S1046;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S1048;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S1043.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1049, DiffPair_float_0 * _S1050, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1051, float3  _S1052)
{
    s_bwd_prop_blend_background_0(_S1049, _S1050, _S1051, _S1052);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_3, float alpha_1, float3  background_1, float3  v_out_rgb_1, float3  * v_rgb_1, float * v_alpha_1, float3  * v_background_0)
{
    float3  _S1053 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_3;
    (&p_rgb_0)->differential_0 = _S1053;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S1053;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_1);
    *v_rgb_1 = p_rgb_0.differential_0;
    *v_alpha_1 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpy_5, float dOut_20)
{
    if(((*dpx_16).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = 0.0f;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_16).primal_0), ((*dpy_5).primal_0)));
        DiffPair_float_0 _S1054 = *dpx_16;
        float _S1055 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_20;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S1055;
        float _S1056 = val_0 * (F32_log((_S1054.primal_0))) * dOut_20;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S1056;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_4)
{
    float3  _S1057 = rgb_4;
    float _S1058;
    if((rgb_4.x) < 0.00313080009073019f)
    {
        _S1058 = _S1057.x * 12.92000007629394531f;
    }
    else
    {
        _S1058 = 1.0549999475479126f * (F32_pow((_S1057.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1057)->x) = _S1058;
    if((_S1057.y) < 0.00313080009073019f)
    {
        _S1058 = _S1057.y * 12.92000007629394531f;
    }
    else
    {
        _S1058 = 1.0549999475479126f * (F32_pow((_S1057.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1057)->y) = _S1058;
    if((_S1057.z) < 0.00313080009073019f)
    {
        _S1058 = _S1057.z * 12.92000007629394531f;
    }
    else
    {
        _S1058 = 1.0549999475479126f * (F32_pow((_S1057.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1057)->z) = _S1058;
    return _S1057;
}

inline __device__ float s_primal_ctx_pow_0(float _S1059, float _S1060)
{
    return (F32_pow((_S1059), (_S1060)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1061, DiffPair_float_0 * _S1062, float _S1063)
{
    _d_pow_0(_S1061, _S1062, _S1063);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_2, float3  _s_dOut_9)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1064 = *dprgb_2;
    float _S1065 = (*dprgb_2).primal_0.x;
    bool _S1066 = _S1065 < 0.00313080009073019f;
    float _S1067;
    if(_S1066)
    {
        _S1067 = _S1065 * 12.92000007629394531f;
    }
    else
    {
        _S1067 = 1.0549999475479126f * s_primal_ctx_pow_0(_S1065, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S1068 = _S1064.primal_0;
    *&((&_S1068)->x) = _S1067;
    float _S1069 = _S1068.y;
    bool _S1070 = _S1069 < 0.00313080009073019f;
    if(_S1070)
    {
        _S1067 = _S1069 * 12.92000007629394531f;
    }
    else
    {
        _S1067 = 1.0549999475479126f * s_primal_ctx_pow_0(_S1069, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S1068)->y) = _S1067;
    float _S1071 = _S1068.z;
    bool _S1072 = _S1071 < 0.00313080009073019f;
    _S1068 = _s_dOut_9;
    *&((&_S1068)->z) = 0.0f;
    if(_S1072)
    {
        _S1067 = 12.92000007629394531f * _s_dOut_9.z;
    }
    else
    {
        float _S1073 = 1.0549999475479126f * _s_dOut_9.z;
        DiffPair_float_0 _S1074;
        (&_S1074)->primal_0 = _S1071;
        (&_S1074)->differential_0 = 0.0f;
        DiffPair_float_0 _S1075;
        (&_S1075)->primal_0 = 0.4166666567325592f;
        (&_S1075)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1074, &_S1075, _S1073);
        _S1067 = _S1074.differential_0;
    }
    float3  _S1076 = _S1068 + make_float3 (0.0f, 0.0f, _S1067);
    _S1068 = _S1076;
    *&((&_S1068)->y) = 0.0f;
    if(_S1070)
    {
        _S1067 = 12.92000007629394531f * _S1076.y;
    }
    else
    {
        float _S1077 = 1.0549999475479126f * _S1076.y;
        DiffPair_float_0 _S1078;
        (&_S1078)->primal_0 = _S1069;
        (&_S1078)->differential_0 = 0.0f;
        DiffPair_float_0 _S1079;
        (&_S1079)->primal_0 = 0.4166666567325592f;
        (&_S1079)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1078, &_S1079, _S1077);
        _S1067 = _S1078.differential_0;
    }
    float3  _S1080 = _S1068 + make_float3 (0.0f, _S1067, 0.0f);
    _S1068 = _S1080;
    *&((&_S1068)->x) = 0.0f;
    if(_S1066)
    {
        _S1067 = 12.92000007629394531f * _S1080.x;
    }
    else
    {
        float _S1081 = 1.0549999475479126f * _S1080.x;
        DiffPair_float_0 _S1082;
        (&_S1082)->primal_0 = _S1065;
        (&_S1082)->differential_0 = 0.0f;
        DiffPair_float_0 _S1083;
        (&_S1083)->primal_0 = 0.4166666567325592f;
        (&_S1083)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1082, &_S1083, _S1081);
        _S1067 = _S1082.differential_0;
    }
    float3  _S1084 = _S1068 + make_float3 (_S1067, 0.0f, 0.0f);
    dprgb_2->primal_0 = (*dprgb_2).primal_0;
    dprgb_2->differential_0 = _S1084;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1085, float3  _S1086)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S1085, _S1086);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_5, float3  v_out_rgb_2)
{
    float3  _S1087 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_5;
    (&p_rgb_1)->differential_0 = _S1087;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_2);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S1088 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_6 = _S1088;
    FixedArray<float, 10>  _S1089 = dist_coeffs_6;
    bool _S1090 = undistort_point_0(_S1088, &_S1089, int(12), &uv_6);
    if(!_S1090)
    {
        int3  _S1091 = make_int3 (int(0));
        float3  _S1092 = make_float3 ((float)_S1091.x, (float)_S1091.y, (float)_S1091.z);
        return _S1092;
    }
    float3  raydir_6;
    if(is_fisheye_4)
    {
        float theta_4 = length_1(uv_6);
        float3  raydir_7 = make_float3 ((uv_6 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_6 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_0)
        {
            raydir_6 = raydir_7 / make_float3 (raydir_7.z);
        }
        else
        {
            raydir_6 = raydir_7;
        }
    }
    else
    {
        float3  raydir_8 = make_float3 (uv_6.x, uv_6.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_6 = normalize_0(raydir_8);
        }
        else
        {
            raydir_6 = raydir_8;
        }
    }
    return raydir_6;
}

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  points_0)
{
    float3  normal_0 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    float3  normal_1;
    if((dot_0(normal_0, normal_0)) != 0.0f)
    {
        normal_1 = normal_0 / make_float3 (length_2(normal_0));
    }
    else
    {
        normal_1 = normal_0;
    }
    return normal_1;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0
{
    FixedArray<float3 , 4>  primal_0;
    FixedArray<float3 , 4>  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_10)
{
    float _S1093 = (*dpx_17).primal_0.x;
    float _S1094 = (*dpx_17).primal_0.y;
    float _S1095 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S1096;
    (&_S1096)->primal_0 = _S1093 * _S1093 + _S1094 * _S1094 + _S1095 * _S1095;
    (&_S1096)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1096, _s_dOut_10);
    float _S1097 = (*dpx_17).primal_0.z * _S1096.differential_0;
    float _S1098 = _S1097 + _S1097;
    float _S1099 = (*dpx_17).primal_0.y * _S1096.differential_0;
    float _S1100 = _S1099 + _S1099;
    float _S1101 = (*dpx_17).primal_0.x * _S1096.differential_0;
    float _S1102 = _S1101 + _S1101;
    float3  _S1103 = make_float3 (0.0f);
    *&((&_S1103)->z) = _S1098;
    *&((&_S1103)->y) = _S1100;
    *&((&_S1103)->x) = _S1102;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S1103;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1104, float _S1105)
{
    s_bwd_prop_length_impl_1(_S1104, _S1105);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_11)
{
    float3  _S1106 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S1107 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S1108 = s_primal_ctx_cross_0(dx_0, _S1107);
    bool _S1109 = (s_primal_ctx_dot_0(_S1108, _S1108)) != 0.0f;
    float3  _S1110;
    float3  _S1111;
    if(_S1109)
    {
        float _S1112 = length_2(_S1108);
        float3  _S1113 = make_float3 (_S1112);
        _S1110 = make_float3 (_S1112 * _S1112);
        _S1111 = _S1113;
    }
    else
    {
        _S1110 = _S1106;
        _S1111 = _S1106;
    }
    if(_S1109)
    {
        float3  _S1114 = _s_dOut_11 / _S1110;
        float3  _S1115 = _S1108 * - _S1114;
        float3  _S1116 = _S1111 * _S1114;
        float _S1117 = _S1115.x + _S1115.y + _S1115.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1118;
        (&_S1118)->primal_0 = _S1108;
        (&_S1118)->differential_0 = _S1106;
        s_bwd_length_impl_1(&_S1118, _S1117);
        _S1110 = _S1116 + _S1118.differential_0;
    }
    else
    {
        _S1110 = _s_dOut_11;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1119;
    (&_S1119)->primal_0 = _S1108;
    (&_S1119)->differential_0 = _S1106;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1120;
    (&_S1120)->primal_0 = _S1108;
    (&_S1120)->differential_0 = _S1106;
    s_bwd_prop_dot_0(&_S1119, &_S1120, 0.0f);
    float3  _S1121 = _S1120.differential_0 + _S1119.differential_0 + _S1110;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1122;
    (&_S1122)->primal_0 = dx_0;
    (&_S1122)->differential_0 = _S1106;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1123;
    (&_S1123)->primal_0 = _S1107;
    (&_S1123)->differential_0 = _S1106;
    s_bwd_prop_cross_0(&_S1122, &_S1123, _S1121);
    float3  s_diff_dy_T_0 = - _S1123.differential_0;
    float3  _S1124 = - s_diff_dy_T_0;
    float3  _S1125 = - _S1122.differential_0;
    FixedArray<float3 , 4>  _S1126;
    _S1126[int(0)] = _S1106;
    _S1126[int(1)] = _S1106;
    _S1126[int(2)] = _S1106;
    _S1126[int(3)] = _S1106;
    _S1126[int(2)] = _S1124;
    _S1126[int(3)] = s_diff_dy_T_0;
    _S1126[int(0)] = _S1125;
    _S1126[int(1)] = _S1122.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S1126;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S1127, float3  _S1128)
{
    s_bwd_prop_points_to_normal_0(_S1127, _S1128);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S1129 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S1129;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_7, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S1130 = float2 {intrins_1.z, intrins_1.w};
    float2  _S1131 = float2 {intrins_1.x, intrins_1.y};
    float2  _S1132 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S1130) / _S1131;
    float2  uv_7 = _S1132;
    FixedArray<float, 10>  _S1133 = dist_coeffs_7;
    bool _S1134 = undistort_point_0(_S1132, &_S1133, int(12), &uv_7);
    if(!_S1134)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_9;
    if(is_fisheye_5)
    {
        float theta_5 = length_1(uv_7);
        float3  raydir_10 = make_float3 ((uv_7 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_7 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        if(!is_ray_depth_1)
        {
            raydir_9 = raydir_10 / make_float3 (raydir_10.z);
        }
        else
        {
            raydir_9 = raydir_10;
        }
    }
    else
    {
        float3  raydir_11 = make_float3 (uv_7.x, uv_7.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_9 = normalize_0(raydir_11);
        }
        else
        {
            raydir_9 = raydir_11;
        }
    }
    points_2[int(0)] = make_float3 (depths_0.x) * raydir_9;
    float2  _S1135 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S1130) / _S1131;
    float2  uv_8 = _S1135;
    FixedArray<float, 10>  _S1136 = dist_coeffs_7;
    bool _S1137 = undistort_point_0(_S1135, &_S1136, int(12), &uv_8);
    if(!_S1137)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_6 = length_1(uv_8);
        float3  raydir_12 = make_float3 ((uv_8 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (uv_8 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
        if(!is_ray_depth_1)
        {
            raydir_9 = raydir_12 / make_float3 (raydir_12.z);
        }
        else
        {
            raydir_9 = raydir_12;
        }
    }
    else
    {
        float3  raydir_13 = make_float3 (uv_8.x, uv_8.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_9 = normalize_0(raydir_13);
        }
        else
        {
            raydir_9 = raydir_13;
        }
    }
    points_2[int(1)] = make_float3 (depths_0.y) * raydir_9;
    float2  _S1138 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S1130) / _S1131;
    float2  uv_9 = _S1138;
    FixedArray<float, 10>  _S1139 = dist_coeffs_7;
    bool _S1140 = undistort_point_0(_S1138, &_S1139, int(12), &uv_9);
    if(!_S1140)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_7 = length_1(uv_9);
        float3  raydir_14 = make_float3 ((uv_9 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (uv_9 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        if(!is_ray_depth_1)
        {
            raydir_9 = raydir_14 / make_float3 (raydir_14.z);
        }
        else
        {
            raydir_9 = raydir_14;
        }
    }
    else
    {
        float3  raydir_15 = make_float3 (uv_9.x, uv_9.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_9 = normalize_0(raydir_15);
        }
        else
        {
            raydir_9 = raydir_15;
        }
    }
    points_2[int(2)] = make_float3 (depths_0.z) * raydir_9;
    float2  _S1141 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S1130) / _S1131;
    float2  uv_10 = _S1141;
    FixedArray<float, 10>  _S1142 = dist_coeffs_7;
    bool _S1143 = undistort_point_0(_S1141, &_S1142, int(12), &uv_10);
    if(!_S1143)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_8 = length_1(uv_10);
        float3  raydir_16 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).x, (uv_10 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).y, (F32_cos((theta_8))));
        if(!is_ray_depth_1)
        {
            raydir_9 = raydir_16 / make_float3 (raydir_16.z);
        }
        else
        {
            raydir_9 = raydir_16;
        }
    }
    else
    {
        float3  raydir_17 = make_float3 (uv_10.x, uv_10.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_9 = normalize_0(raydir_17);
        }
        else
        {
            raydir_9 = raydir_17;
        }
    }
    points_2[int(3)] = make_float3 (depths_0.w) * raydir_9;
    float3  normal_2 = cross_0(points_2[int(1)] - points_2[int(0)], - (points_2[int(3)] - points_2[int(2)]));
    float3  normal_3;
    if((dot_0(normal_2, normal_2)) != 0.0f)
    {
        normal_3 = normal_2 / make_float3 (length_2(normal_2));
    }
    else
    {
        normal_3 = normal_2;
    }
    return normal_3;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S1144;
    bool _S1145;
    float2  _S1146;
    bool _S1147;
    float2  _S1148;
    bool _S1149;
    float2  _S1150;
    bool _S1151;
};

inline __device__ float s_primal_ctx_sin_0(float _S1152)
{
    return (F32_sin((_S1152)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1153)
{
    return (F32_cos((_S1153)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1154 = make_float2 (0.0f);
    _s_diff_ctx_0->_S1144 = _S1154;
    _s_diff_ctx_0->_S1145 = false;
    _s_diff_ctx_0->_S1146 = _S1154;
    _s_diff_ctx_0->_S1147 = false;
    _s_diff_ctx_0->_S1148 = _S1154;
    _s_diff_ctx_0->_S1149 = false;
    _s_diff_ctx_0->_S1150 = _S1154;
    _s_diff_ctx_0->_S1151 = false;
    _s_diff_ctx_0->_S1146 = _S1154;
    _s_diff_ctx_0->_S1147 = false;
    _s_diff_ctx_0->_S1148 = _S1154;
    _s_diff_ctx_0->_S1149 = false;
    _s_diff_ctx_0->_S1150 = _S1154;
    _s_diff_ctx_0->_S1151 = false;
    float3  _S1155 = make_float3 (0.0f);
    float2  _S1156 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1157 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1158 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1156) / _S1157;
    float2  _S1159 = _S1158;
    bool _S1160 = undistort_point_0(_S1158, dist_coeffs_8, int(12), &_S1159);
    _s_diff_ctx_0->_S1144 = _S1159;
    _s_diff_ctx_0->_S1145 = _S1160;
    float2  uv_11 = _S1159;
    bool _S1161 = !_S1160;
    float3  normal_4;
    if(_S1161)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1162 = !_S1161;
    int _S1163;
    FixedArray<float3 , 4>  points_3;
    if(_S1162)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1164 = length_1(uv_11);
            float3  raydir_19 = make_float3 ((uv_11 / make_float2 ((F32_max((_S1164), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1164))).x, (uv_11 / make_float2 ((F32_max((_S1164), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1164))).y, s_primal_ctx_cos_0(_S1164));
            if(!is_ray_depth_2)
            {
                raydir_18 = raydir_19 / make_float3 (raydir_19.z);
            }
            else
            {
                raydir_18 = raydir_19;
            }
        }
        else
        {
            float3  raydir_20 = make_float3 (uv_11.x, uv_11.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_18 = normalize_0(raydir_20);
            }
            else
            {
                raydir_18 = raydir_20;
            }
        }
        float3  _S1165 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1166 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1156) / _S1157;
        float2  _S1167 = _S1166;
        bool _S1168 = undistort_point_0(_S1166, dist_coeffs_8, int(12), &_S1167);
        _s_diff_ctx_0->_S1146 = _S1167;
        _s_diff_ctx_0->_S1147 = _S1168;
        float2  uv_12 = _S1167;
        bool _S1169 = !_S1168;
        if(_S1169)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1170 = !_S1169;
        if(_S1170)
        {
            if(is_fisheye_6)
            {
                float _S1171 = length_1(uv_12);
                float3  raydir_21 = make_float3 ((uv_12 / make_float2 ((F32_max((_S1171), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1171))).x, (uv_12 / make_float2 ((F32_max((_S1171), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1171))).y, s_primal_ctx_cos_0(_S1171));
                if(!is_ray_depth_2)
                {
                    raydir_18 = raydir_21 / make_float3 (raydir_21.z);
                }
                else
                {
                    raydir_18 = raydir_21;
                }
            }
            else
            {
                float3  raydir_22 = make_float3 (uv_12.x, uv_12.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_18 = normalize_0(raydir_22);
                }
                else
                {
                    raydir_18 = raydir_22;
                }
            }
            float3  _S1172 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1163 = int(2);
            points_3[int(0)] = _S1165;
            points_3[int(1)] = _S1172;
            points_3[int(2)] = _S1155;
            points_3[int(3)] = _S1155;
        }
        else
        {
            _S1163 = int(0);
            points_3[int(0)] = _S1165;
            points_3[int(1)] = _S1155;
            points_3[int(2)] = _S1155;
            points_3[int(3)] = _S1155;
        }
        bool _runFlag_0;
        if(_S1163 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1162;
            _S1163 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1173 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1156) / _S1157;
            float2  _S1174 = _S1173;
            bool _S1175 = undistort_point_0(_S1173, dist_coeffs_8, int(12), &_S1174);
            _s_diff_ctx_0->_S1148 = _S1174;
            _s_diff_ctx_0->_S1149 = _S1175;
            float2  uv_13 = _S1174;
            if(!_S1175)
            {
                float3  _S1176 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1163 = int(0);
                normal_4 = _S1176;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1177 = length_1(uv_13);
                    float3  raydir_23 = make_float3 ((uv_13 / make_float2 ((F32_max((_S1177), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1177))).x, (uv_13 / make_float2 ((F32_max((_S1177), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1177))).y, s_primal_ctx_cos_0(_S1177));
                    if(!is_ray_depth_2)
                    {
                        raydir_18 = raydir_23 / make_float3 (raydir_23.z);
                    }
                    else
                    {
                        raydir_18 = raydir_23;
                    }
                }
                else
                {
                    float3  raydir_24 = make_float3 (uv_13.x, uv_13.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_18 = normalize_0(raydir_24);
                    }
                    else
                    {
                        raydir_18 = raydir_24;
                    }
                }
                points_3[int(2)] = make_float3 (dpdepths_0.z) * raydir_18;
                float2  _S1178 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1156) / _S1157;
                float2  _S1179 = _S1178;
                bool _S1180 = undistort_point_0(_S1178, dist_coeffs_8, int(12), &_S1179);
                _s_diff_ctx_0->_S1150 = _S1179;
                _s_diff_ctx_0->_S1151 = _S1180;
                float2  uv_14 = _S1179;
                bool _S1181 = !_S1180;
                if(_S1181)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1182 = !_S1181;
                int _S1183;
                if(_S1182)
                {
                    if(is_fisheye_6)
                    {
                        float _S1184 = length_1(uv_14);
                        float3  raydir_25 = make_float3 ((uv_14 / make_float2 ((F32_max((_S1184), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1184))).x, (uv_14 / make_float2 ((F32_max((_S1184), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1184))).y, s_primal_ctx_cos_0(_S1184));
                        if(!is_ray_depth_2)
                        {
                            raydir_18 = raydir_25 / make_float3 (raydir_25.z);
                        }
                        else
                        {
                            raydir_18 = raydir_25;
                        }
                    }
                    else
                    {
                        float3  raydir_26 = make_float3 (uv_14.x, uv_14.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_18 = normalize_0(raydir_26);
                        }
                        else
                        {
                            raydir_18 = raydir_26;
                        }
                    }
                    points_3[int(3)] = make_float3 (dpdepths_0.w) * raydir_18;
                    _S1183 = int(2);
                }
                else
                {
                    _S1183 = int(0);
                }
                if(_S1183 != int(2))
                {
                    _runFlag_0 = false;
                    _S1163 = _S1183;
                }
                if(_runFlag_0)
                {
                    _S1163 = int(1);
                }
            }
        }
    }
    else
    {
        _S1163 = int(0);
        points_3[int(0)] = _S1155;
        points_3[int(1)] = _S1155;
        points_3[int(2)] = _S1155;
        points_3[int(3)] = _S1155;
    }
    if(!(_S1163 != int(1)))
    {
        float3  _S1185 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1185, _S1185)) != 0.0f)
        {
            normal_4 = _S1185 / make_float3 (length_2(_S1185));
        }
        else
        {
            normal_4 = _S1185;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_12, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1186 = *dpdepths_1;
    float3  _S1187 = make_float3 (0.0f);
    float2  _S1188 = _s_diff_ctx_1->_S1144;
    bool _S1189 = !!_s_diff_ctx_1->_S1145;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1190;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1191;
    if(_S1189)
    {
        if(is_fisheye_7)
        {
            float _S1192 = length_1(_S1188);
            float3  raydir_31 = make_float3 ((_S1188 / make_float2 ((F32_max((_S1192), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1192))).x, (_S1188 / make_float2 ((F32_max((_S1192), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1192))).y, s_primal_ctx_cos_0(_S1192));
            if(!is_ray_depth_3)
            {
                raydir_27 = raydir_31 / make_float3 (raydir_31.z);
            }
            else
            {
                raydir_27 = raydir_31;
            }
        }
        else
        {
            float3  raydir_32 = make_float3 (_S1188.x, _S1188.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1193 = make_float3 (_S1186.primal_0.x) * raydir_27;
        float2  _S1194 = _s_diff_ctx_1->_S1146;
        bool _S1195 = !!_s_diff_ctx_1->_S1147;
        if(_S1195)
        {
            if(is_fisheye_7)
            {
                float _S1196 = length_1(_S1194);
                float3  raydir_33 = make_float3 ((_S1194 / make_float2 ((F32_max((_S1196), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1196))).x, (_S1194 / make_float2 ((F32_max((_S1196), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1196))).y, s_primal_ctx_cos_0(_S1196));
                if(!is_ray_depth_3)
                {
                    raydir_28 = raydir_33 / make_float3 (raydir_33.z);
                }
                else
                {
                    raydir_28 = raydir_33;
                }
            }
            else
            {
                float3  raydir_34 = make_float3 (_S1194.x, _S1194.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1197 = make_float3 (_S1186.primal_0.y) * raydir_28;
            _S1190 = int(2);
            points_4[int(0)] = _S1193;
            points_4[int(1)] = _S1197;
            points_4[int(2)] = _S1187;
            points_4[int(3)] = _S1187;
        }
        else
        {
            _S1190 = int(0);
            points_4[int(0)] = _S1193;
            points_4[int(1)] = _S1187;
            points_4[int(2)] = _S1187;
            points_4[int(3)] = _S1187;
            raydir_28 = _S1187;
        }
        if(_S1190 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1189;
            _S1190 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S1198 = _s_diff_ctx_1->_S1148;
            if(!_s_diff_ctx_1->_S1149)
            {
                _runFlag_2 = false;
                _S1190 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1199 = length_1(_S1198);
                    float3  raydir_35 = make_float3 ((_S1198 / make_float2 ((F32_max((_S1199), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1199))).x, (_S1198 / make_float2 ((F32_max((_S1199), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1199))).y, s_primal_ctx_cos_0(_S1199));
                    if(!is_ray_depth_3)
                    {
                        raydir_29 = raydir_35 / make_float3 (raydir_35.z);
                    }
                    else
                    {
                        raydir_29 = raydir_35;
                    }
                }
                else
                {
                    float3  raydir_36 = make_float3 (_S1198.x, _S1198.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1186.primal_0.z) * raydir_29;
                float2  _S1200 = _s_diff_ctx_1->_S1150;
                bool _S1201 = !!_s_diff_ctx_1->_S1151;
                int _S1202;
                if(_S1201)
                {
                    if(is_fisheye_7)
                    {
                        float _S1203 = length_1(_S1200);
                        float3  raydir_37 = make_float3 ((_S1200 / make_float2 ((F32_max((_S1203), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1203))).x, (_S1200 / make_float2 ((F32_max((_S1203), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1203))).y, s_primal_ctx_cos_0(_S1203));
                        if(!is_ray_depth_3)
                        {
                            raydir_30 = raydir_37 / make_float3 (raydir_37.z);
                        }
                        else
                        {
                            raydir_30 = raydir_37;
                        }
                    }
                    else
                    {
                        float3  raydir_38 = make_float3 (_S1200.x, _S1200.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1186.primal_0.w) * raydir_30;
                    _S1202 = int(2);
                }
                else
                {
                    _S1202 = int(0);
                    raydir_30 = _S1187;
                }
                if(_S1202 != int(2))
                {
                    _runFlag_3 = false;
                    _S1190 = _S1202;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1190 = int(1);
                }
                float3  _S1204 = raydir_29;
                _runFlag_3 = _S1201;
                raydir_29 = raydir_30;
                raydir_30 = _S1204;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1187;
                raydir_30 = _S1187;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1187;
            raydir_30 = _S1187;
        }
        float3  _S1205 = raydir_27;
        float3  _S1206 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1191 = _S1195;
        raydir_29 = _S1206;
        raydir_30 = _S1205;
    }
    else
    {
        _S1190 = int(0);
        points_4[int(0)] = _S1187;
        points_4[int(1)] = _S1187;
        points_4[int(2)] = _S1187;
        points_4[int(3)] = _S1187;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1187;
        raydir_28 = _S1187;
        _S1191 = false;
        raydir_29 = _S1187;
        raydir_30 = _S1187;
    }
    bool _S1207 = !(_S1190 != int(1));
    float3  _S1208;
    float3  _S1209;
    float3  _S1210;
    float3  _S1211;
    float3  _S1212;
    bool _S1213;
    if(_S1207)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1214 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1215 = s_primal_ctx_cross_0(dx_1, _S1214);
        bool _S1216 = (s_primal_ctx_dot_0(_S1215, _S1215)) != 0.0f;
        if(_S1216)
        {
            float _S1217 = length_2(_S1215);
            float3  _S1218 = make_float3 (_S1217);
            _S1208 = make_float3 (_S1217 * _S1217);
            _S1209 = _S1218;
        }
        else
        {
            _S1208 = _S1187;
            _S1209 = _S1187;
        }
        float3  _S1219 = _S1209;
        _S1213 = _S1216;
        _S1209 = _S1215;
        _S1210 = _S1219;
        _S1211 = dx_1;
        _S1212 = _S1214;
    }
    else
    {
        _S1213 = false;
        _S1208 = _S1187;
        _S1209 = _S1187;
        _S1210 = _S1187;
        _S1211 = _S1187;
        _S1212 = _S1187;
    }
    float4  _S1220 = make_float4 (0.0f);
    if(_S1207)
    {
        if(_S1213)
        {
            float3  _S1221 = _s_dOut_12 / _S1208;
            float3  _S1222 = _S1209 * - _S1221;
            float3  _S1223 = _S1210 * _S1221;
            float _S1224 = _S1222.x + _S1222.y + _S1222.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1225;
            (&_S1225)->primal_0 = _S1209;
            (&_S1225)->differential_0 = _S1187;
            s_bwd_length_impl_1(&_S1225, _S1224);
            _S1208 = _S1223 + _S1225.differential_0;
        }
        else
        {
            _S1208 = _s_dOut_12;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1226;
        (&_S1226)->primal_0 = _S1209;
        (&_S1226)->differential_0 = _S1187;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1227;
        (&_S1227)->primal_0 = _S1209;
        (&_S1227)->differential_0 = _S1187;
        s_bwd_prop_dot_0(&_S1226, &_S1227, 0.0f);
        float3  _S1228 = _S1227.differential_0 + _S1226.differential_0 + _S1208;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1229;
        (&_S1229)->primal_0 = _S1211;
        (&_S1229)->differential_0 = _S1187;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1230;
        (&_S1230)->primal_0 = _S1212;
        (&_S1230)->differential_0 = _S1187;
        s_bwd_prop_cross_0(&_S1229, &_S1230, _S1228);
        float3  s_diff_dy_T_1 = - _S1230.differential_0;
        float3  _S1231 = - s_diff_dy_T_1;
        float3  _S1232 = - _S1229.differential_0;
        FixedArray<float3 , 4>  _S1233;
        _S1233[int(0)] = _S1187;
        _S1233[int(1)] = _S1187;
        _S1233[int(2)] = _S1187;
        _S1233[int(3)] = _S1187;
        _S1233[int(2)] = _S1231;
        _S1233[int(3)] = s_diff_dy_T_1;
        _S1233[int(0)] = _S1232;
        _S1233[int(1)] = _S1229.differential_0;
        points_4[int(0)] = _S1233[int(0)];
        points_4[int(1)] = _S1233[int(1)];
        points_4[int(2)] = _S1233[int(2)];
        points_4[int(3)] = _S1233[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1187;
        points_4[int(1)] = _S1187;
        points_4[int(2)] = _S1187;
        points_4[int(3)] = _S1187;
    }
    float4  _S1234;
    if(_S1189)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1235 = points_4;
                FixedArray<float3 , 4>  _S1236 = points_4;
                FixedArray<float3 , 4>  _S1237 = points_4;
                FixedArray<float3 , 4>  _S1238 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1239 = raydir_27 * _S1238[int(3)];
                    float _S1240 = _S1239.x + _S1239.y + _S1239.z;
                    float4  _S1241 = _S1220;
                    *&((&_S1241)->w) = _S1240;
                    points_4[int(0)] = _S1235[int(0)];
                    points_4[int(1)] = _S1236[int(1)];
                    points_4[int(2)] = _S1237[int(2)];
                    points_4[int(3)] = _S1187;
                    _S1234 = _S1241;
                }
                else
                {
                    points_4[int(0)] = _S1235[int(0)];
                    points_4[int(1)] = _S1236[int(1)];
                    points_4[int(2)] = _S1237[int(2)];
                    points_4[int(3)] = _S1238[int(3)];
                    _S1234 = _S1220;
                }
                float3  _S1242 = raydir_28 * points_4[int(2)];
                float _S1243 = _S1242.x + _S1242.y + _S1242.z;
                FixedArray<float3 , 4>  _S1244 = points_4;
                FixedArray<float3 , 4>  _S1245 = points_4;
                float4  _S1246 = _S1220;
                *&((&_S1246)->z) = _S1243;
                float4  _S1247 = _S1234 + _S1246;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1244[int(1)];
                points_4[int(2)] = _S1187;
                points_4[int(3)] = _S1245[int(3)];
                _S1234 = _S1247;
            }
            else
            {
                FixedArray<float3 , 4>  _S1248 = points_4;
                FixedArray<float3 , 4>  _S1249 = points_4;
                FixedArray<float3 , 4>  _S1250 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1248[int(1)];
                points_4[int(2)] = _S1249[int(2)];
                points_4[int(3)] = _S1250[int(3)];
                _S1234 = _S1220;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1251 = points_4;
            FixedArray<float3 , 4>  _S1252 = points_4;
            FixedArray<float3 , 4>  _S1253 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1251[int(1)];
            points_4[int(2)] = _S1252[int(2)];
            points_4[int(3)] = _S1253[int(3)];
            _S1234 = _S1220;
        }
        if(_S1191)
        {
            FixedArray<float3 , 4>  _S1254 = points_4;
            float3  _S1255 = raydir_29 * points_4[int(1)];
            float _S1256 = _S1255.x + _S1255.y + _S1255.z;
            float4  _S1257 = _S1220;
            *&((&_S1257)->y) = _S1256;
            float4  _S1258 = _S1234 + _S1257;
            points_4[int(0)] = _S1187;
            points_4[int(1)] = _S1187;
            points_4[int(2)] = _S1187;
            points_4[int(3)] = _S1187;
            raydir_27 = _S1254[int(0)];
            _S1234 = _S1258;
        }
        else
        {
            FixedArray<float3 , 4>  _S1259 = points_4;
            FixedArray<float3 , 4>  _S1260 = points_4;
            FixedArray<float3 , 4>  _S1261 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1259[int(1)];
            points_4[int(2)] = _S1260[int(2)];
            points_4[int(3)] = _S1261[int(3)];
            raydir_27 = _S1187;
        }
        float3  _S1262 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1263 = _S1262.x + _S1262.y + _S1262.z;
        float4  _S1264 = _S1220;
        *&((&_S1264)->x) = _S1263;
        _S1234 = _S1234 + _S1264;
    }
    else
    {
        _S1234 = _S1220;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1234;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1265, float4  _S1266, FixedArray<float, 10>  * _S1267, bool _S1268, bool _S1269, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1270, float3  _S1271)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1272;
    float3  _S1273 = s_primal_ctx_depth_to_normal_0(_S1265, _S1266, _S1267, _S1268, _S1269, (*_S1270).primal_0, &_S1272);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1274 = _S1272;
    s_bwd_prop_depth_to_normal_0(_S1265, _S1266, _S1267, _S1268, _S1269, _S1270, _S1271, &_S1274);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_10, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1275 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1275;
    FixedArray<float, 10>  _S1276 = dist_coeffs_10;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S1276, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_11, bool is_fisheye_9)
{
    float2  _S1277 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_15 = _S1277;
    FixedArray<float, 10>  _S1278 = dist_coeffs_11;
    bool _S1279 = undistort_point_0(_S1277, &_S1278, int(12), &uv_15);
    if(!_S1279)
    {
        return 0.0f;
    }
    float3  raydir_39;
    if(is_fisheye_9)
    {
        float theta_9 = length_1(uv_15);
        float3  raydir_40 = make_float3 ((uv_15 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).x, (uv_15 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).y, (F32_cos((theta_9))));
        raydir_39 = raydir_40 / make_float3 (raydir_40.z);
    }
    else
    {
        raydir_39 = make_float3 (uv_15.x, uv_15.y, 1.0f);
    }
    return float((F32_sign((raydir_39.z)))) / length_2(raydir_39);
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_18, float dOut_21)
{
    float _S1280 = (F32_exp2(((*dpx_18).primal_0))) * 0.69314718246459961f * dOut_21;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1280;
    return;
}

inline __device__ float3  apply_ppisp(float3  rgb_in_0, float2  pix_coord_0, float2  image_center_0, float2  img_size_0, FixedArray<float, 36>  params_0)
{
    PPISPParams_0 p_0;
    (&p_0)->exposure_1 = params_0[int(0)];
    (&(&p_0)->vignette_params_1[int(0)])->cx_0 = params_0[int(1)];
    (&(&p_0)->vignette_params_1[int(0)])->cy_0 = params_0[int(2)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha0_0 = params_0[int(3)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha1_0 = params_0[int(4)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha2_0 = params_0[int(5)];
    (&(&p_0)->vignette_params_1[int(1)])->cx_0 = params_0[int(6)];
    (&(&p_0)->vignette_params_1[int(1)])->cy_0 = params_0[int(7)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha0_0 = params_0[int(8)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha1_0 = params_0[int(9)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha2_0 = params_0[int(10)];
    (&(&p_0)->vignette_params_1[int(2)])->cx_0 = params_0[int(11)];
    (&(&p_0)->vignette_params_1[int(2)])->cy_0 = params_0[int(12)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha0_0 = params_0[int(13)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha1_0 = params_0[int(14)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha2_0 = params_0[int(15)];
    *&((&(&(&p_0)->color_params_1)->b_0)->x) = params_0[int(16)];
    *&((&(&(&p_0)->color_params_1)->b_0)->y) = params_0[int(17)];
    *&((&(&(&p_0)->color_params_1)->r_0)->x) = params_0[int(18)];
    *&((&(&(&p_0)->color_params_1)->r_0)->y) = params_0[int(19)];
    *&((&(&(&p_0)->color_params_1)->g_0)->x) = params_0[int(20)];
    *&((&(&(&p_0)->color_params_1)->g_0)->y) = params_0[int(21)];
    *&((&(&(&p_0)->color_params_1)->n_0)->x) = params_0[int(22)];
    *&((&(&(&p_0)->color_params_1)->n_0)->y) = params_0[int(23)];
    (&(&p_0)->crf_params_1[int(0)])->toe_0 = params_0[int(24)];
    (&(&p_0)->crf_params_1[int(0)])->shoulder_0 = params_0[int(25)];
    (&(&p_0)->crf_params_1[int(0)])->gamma_0 = params_0[int(26)];
    (&(&p_0)->crf_params_1[int(0)])->center_0 = params_0[int(27)];
    (&(&p_0)->crf_params_1[int(1)])->toe_0 = params_0[int(28)];
    (&(&p_0)->crf_params_1[int(1)])->shoulder_0 = params_0[int(29)];
    (&(&p_0)->crf_params_1[int(1)])->gamma_0 = params_0[int(30)];
    (&(&p_0)->crf_params_1[int(1)])->center_0 = params_0[int(31)];
    (&(&p_0)->crf_params_1[int(2)])->toe_0 = params_0[int(32)];
    (&(&p_0)->crf_params_1[int(2)])->shoulder_0 = params_0[int(33)];
    (&(&p_0)->crf_params_1[int(2)])->gamma_0 = params_0[int(34)];
    (&(&p_0)->crf_params_1[int(2)])->center_0 = params_0[int(35)];
    PPISPParams_0 _S1281 = p_0;
    float _S1282 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1283 = (pix_coord_0.x - image_center_0.x) / _S1282;
    float _S1284 = (pix_coord_0.y - image_center_0.y) / _S1282;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_2 = _S1283 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S1284 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_5 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_5 * r2_5;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_5) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1283 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S1284 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_6 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_6 * r2_6;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_6) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1283 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S1284 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_7 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_7 * r2_7;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_7) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_7 + 1.0f, 0.0f, 1.0f);
    float3  _S1285 = rgb_out_0;
    float2  bd_0 = mul_4(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_4(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_4(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_4(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S1286 = 0.3333333432674408f + nd_0.x;
    float _S1287 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_2 = mul_1(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1287, 1.0f, 0.0f, - _S1286, - _S1287, _S1286, 0.0f), T_0);
    float3  r0_0 = make_float3 (M_2.rows[int(0)].x, M_2.rows[int(0)].y, M_2.rows[int(0)].z);
    float3  r1_0 = make_float3 (M_2.rows[int(1)].x, M_2.rows[int(1)].y, M_2.rows[int(1)].z);
    float3  r2_8 = make_float3 (M_2.rows[int(2)].x, M_2.rows[int(2)].y, M_2.rows[int(2)].z);
    float3  lambda_v_0 = cross_0(r0_0, r1_0);
    float3  lambda_v_1;
    if((dot_0(lambda_v_0, lambda_v_0)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_2 = cross_0(r0_0, r2_8);
        if((dot_0(lambda_v_2, lambda_v_2)) < 9.99999968265522539e-21f)
        {
            lambda_v_1 = cross_0(r1_0, r2_8);
        }
        else
        {
            lambda_v_1 = lambda_v_2;
        }
    }
    else
    {
        lambda_v_1 = lambda_v_0;
    }
    Matrix<float, 3, 3>  H_0 = mul_1(mul_1(T_0, makeMatrix<float, 3, 3> (lambda_v_1.x, 0.0f, 0.0f, 0.0f, lambda_v_1.y, 0.0f, 0.0f, 0.0f, lambda_v_1.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_1;
    if((F32_abs((H_0.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_1 = H_0 * makeMatrix<float, 3, 3> (1.0f / H_0.rows[int(2)].z);
    }
    else
    {
        H_1 = H_0;
    }
    float _S1288 = _S1285.x;
    float _S1289 = _S1285.y;
    float intensity_0 = _S1288 + _S1289 + _S1285.z;
    float3  rgi_out_0 = mul_3(H_1, make_float3 (_S1288, _S1289, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1290 = rgi_out_1.x;
    float _S1291 = rgi_out_1.y;
    float3  _S1292 = clamp_1(make_float3 (_S1290, _S1291, rgi_out_1.z - _S1290 - _S1291), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1293 = _S1292.x;
    float _S1294 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(0)].toe_0))))));
    float _S1295 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(0)].shoulder_0))))));
    float _S1296 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(0)].gamma_0))))));
    float _S1297 = 1.0f / (1.0f + (F32_exp((- _S1281.crf_params_1[int(0)].center_0))));
    float a_1 = _S1295 * _S1297 / lerp_0(_S1294, _S1295, _S1297);
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1293 <= _S1297)
    {
        y_10 = a_1 * (F32_pow((_S1293 / _S1297), (_S1294)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1293) / (1.0f - _S1297)), (_S1295)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1296)));
    float _S1298 = _S1292.y;
    float _S1299 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(1)].toe_0))))));
    float _S1300 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(1)].shoulder_0))))));
    float _S1301 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(1)].gamma_0))))));
    float _S1302 = 1.0f / (1.0f + (F32_exp((- _S1281.crf_params_1[int(1)].center_0))));
    float a_2 = _S1300 * _S1302 / lerp_0(_S1299, _S1300, _S1302);
    float b_3 = 1.0f - a_2;
    if(_S1298 <= _S1302)
    {
        y_10 = a_2 * (F32_pow((_S1298 / _S1302), (_S1299)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1298) / (1.0f - _S1302)), (_S1300)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1301)));
    float _S1303 = _S1292.z;
    float _S1304 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(2)].toe_0))))));
    float _S1305 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(2)].shoulder_0))))));
    float _S1306 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1281.crf_params_1[int(2)].gamma_0))))));
    float _S1307 = 1.0f / (1.0f + (F32_exp((- _S1281.crf_params_1[int(2)].center_0))));
    float a_3 = _S1305 * _S1307 / lerp_0(_S1304, _S1305, _S1307);
    float b_4 = 1.0f - a_3;
    if(_S1303 <= _S1307)
    {
        y_10 = a_3 * (F32_pow((_S1303 / _S1307), (_S1304)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1303) / (1.0f - _S1307)), (_S1305)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1306)));
    return rgb_out_1;
}

inline __device__ float3  apply_ppisp_rqs(float3  rgb_in_1, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, FixedArray<float, 39>  params_1)
{
    PPISPParamsRQS_0 p_1;
    (&p_1)->exposure_0 = params_1[int(0)];
    (&(&p_1)->vignette_params_0[int(0)])->cx_0 = params_1[int(1)];
    (&(&p_1)->vignette_params_0[int(0)])->cy_0 = params_1[int(2)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha0_0 = params_1[int(3)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha1_0 = params_1[int(4)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha2_0 = params_1[int(5)];
    (&(&p_1)->vignette_params_0[int(1)])->cx_0 = params_1[int(6)];
    (&(&p_1)->vignette_params_0[int(1)])->cy_0 = params_1[int(7)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha0_0 = params_1[int(8)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha1_0 = params_1[int(9)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha2_0 = params_1[int(10)];
    (&(&p_1)->vignette_params_0[int(2)])->cx_0 = params_1[int(11)];
    (&(&p_1)->vignette_params_0[int(2)])->cy_0 = params_1[int(12)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha0_0 = params_1[int(13)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha1_0 = params_1[int(14)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha2_0 = params_1[int(15)];
    *&((&(&(&p_1)->color_params_0)->b_0)->x) = params_1[int(16)];
    *&((&(&(&p_1)->color_params_0)->b_0)->y) = params_1[int(17)];
    *&((&(&(&p_1)->color_params_0)->r_0)->x) = params_1[int(18)];
    *&((&(&(&p_1)->color_params_0)->r_0)->y) = params_1[int(19)];
    *&((&(&(&p_1)->color_params_0)->g_0)->x) = params_1[int(20)];
    *&((&(&(&p_1)->color_params_0)->g_0)->y) = params_1[int(21)];
    *&((&(&(&p_1)->color_params_0)->n_0)->x) = params_1[int(22)];
    *&((&(&(&p_1)->color_params_0)->n_0)->y) = params_1[int(23)];
    (&(&p_1)->crf_params_0[int(0)])->g0_0 = params_1[int(24)];
    (&(&p_1)->crf_params_0[int(0)])->g1_0 = params_1[int(25)];
    (&(&p_1)->crf_params_0[int(0)])->x0_0 = params_1[int(26)];
    (&(&p_1)->crf_params_0[int(0)])->y0_0 = params_1[int(27)];
    (&(&p_1)->crf_params_0[int(0)])->gc_0 = params_1[int(28)];
    (&(&p_1)->crf_params_0[int(1)])->g0_0 = params_1[int(29)];
    (&(&p_1)->crf_params_0[int(1)])->g1_0 = params_1[int(30)];
    (&(&p_1)->crf_params_0[int(1)])->x0_0 = params_1[int(31)];
    (&(&p_1)->crf_params_0[int(1)])->y0_0 = params_1[int(32)];
    (&(&p_1)->crf_params_0[int(1)])->gc_0 = params_1[int(33)];
    (&(&p_1)->crf_params_0[int(2)])->g0_0 = params_1[int(34)];
    (&(&p_1)->crf_params_0[int(2)])->g1_0 = params_1[int(35)];
    (&(&p_1)->crf_params_0[int(2)])->x0_0 = params_1[int(36)];
    (&(&p_1)->crf_params_0[int(2)])->y0_0 = params_1[int(37)];
    (&(&p_1)->crf_params_0[int(2)])->gc_0 = params_1[int(38)];
    PPISPParamsRQS_0 _S1308 = p_1;
    float _S1309 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S1310 = (pix_coord_1.x - image_center_1.x) / _S1309;
    float _S1311 = (pix_coord_1.y - image_center_1.y) / _S1309;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_5 = _S1310 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S1311 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_9 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_9 * r2_9;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_9) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float dx_6 = _S1310 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S1311 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_10 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_10 * r2_10;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_10) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_10 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S1310 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S1311 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_11 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_11 * r2_11;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_11) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_11 + 1.0f, 0.0f, 1.0f);
    float3  _S1312 = rgb_out_2;
    float2  bd_1 = mul_4(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_4(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_4(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_4(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S1313 = 0.3333333432674408f + nd_1.x;
    float _S1314 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_1(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1314, 1.0f, 0.0f, - _S1313, - _S1314, _S1313, 0.0f), T_1);
    float3  r0_1 = make_float3 (M_3.rows[int(0)].x, M_3.rows[int(0)].y, M_3.rows[int(0)].z);
    float3  r1_1 = make_float3 (M_3.rows[int(1)].x, M_3.rows[int(1)].y, M_3.rows[int(1)].z);
    float3  r2_12 = make_float3 (M_3.rows[int(2)].x, M_3.rows[int(2)].y, M_3.rows[int(2)].z);
    float3  lambda_v_3 = cross_0(r0_1, r1_1);
    float3  lambda_v_4;
    if((dot_0(lambda_v_3, lambda_v_3)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_5 = cross_0(r0_1, r2_12);
        if((dot_0(lambda_v_5, lambda_v_5)) < 9.99999968265522539e-21f)
        {
            lambda_v_4 = cross_0(r1_1, r2_12);
        }
        else
        {
            lambda_v_4 = lambda_v_5;
        }
    }
    else
    {
        lambda_v_4 = lambda_v_3;
    }
    Matrix<float, 3, 3>  H_2 = mul_1(mul_1(T_1, makeMatrix<float, 3, 3> (lambda_v_4.x, 0.0f, 0.0f, 0.0f, lambda_v_4.y, 0.0f, 0.0f, 0.0f, lambda_v_4.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_3;
    if((F32_abs((H_2.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_3 = H_2 * makeMatrix<float, 3, 3> (1.0f / H_2.rows[int(2)].z);
    }
    else
    {
        H_3 = H_2;
    }
    float _S1315 = _S1312.x;
    float _S1316 = _S1312.y;
    float intensity_1 = _S1315 + _S1316 + _S1312.z;
    float3  rgi_out_2 = mul_3(H_3, make_float3 (_S1315, _S1316, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / (rgi_out_2.z + 0.00000999999974738f));
    float _S1317 = rgi_out_3.x;
    float _S1318 = rgi_out_3.y;
    float3  _S1319 = clamp_1(make_float3 (_S1317, _S1318, rgi_out_3.z - _S1317 - _S1318), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S1320 = _S1319.x;
    float g0_1 = (F32_exp((_S1308.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S1308.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S1308.crf_params_0[int(0)].gc_0)));
    float y_11;
    if(_S1320 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S1320 / x0_1;
        float _S1321 = 1.0f - t0_0;
        y_11 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S1321) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S1321);
    }
    else
    {
        float _S1322 = 1.0f - y0_1;
        float _S1323 = 1.0f - x0_1;
        float s1_4 = _S1322 / _S1323;
        float t1_0 = (_S1320 - x0_1) / _S1323;
        float _S1324 = 1.0f - t1_0;
        y_11 = y0_1 + _S1322 * (s1_4 * t1_0 * t1_0 + gc_1 * t1_0 * _S1324) / (s1_4 + (gc_1 + g1_1 - 2.0f * s1_4) * t1_0 * _S1324);
    }
    *&((&rgb_out_3)->x) = y_11;
    float _S1325 = _S1319.y;
    float g0_2 = (F32_exp((_S1308.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S1308.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S1308.crf_params_0[int(1)].gc_0)));
    if(_S1325 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S1325 / x0_2;
        float _S1326 = 1.0f - t0_1;
        y_11 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S1326) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S1326);
    }
    else
    {
        float _S1327 = 1.0f - y0_2;
        float _S1328 = 1.0f - x0_2;
        float s1_5 = _S1327 / _S1328;
        float t1_1 = (_S1325 - x0_2) / _S1328;
        float _S1329 = 1.0f - t1_1;
        y_11 = y0_2 + _S1327 * (s1_5 * t1_1 * t1_1 + gc_2 * t1_1 * _S1329) / (s1_5 + (gc_2 + g1_2 - 2.0f * s1_5) * t1_1 * _S1329);
    }
    *&((&rgb_out_3)->y) = y_11;
    float _S1330 = _S1319.z;
    float g0_3 = (F32_exp((_S1308.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S1308.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S1308.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S1308.crf_params_0[int(2)].gc_0)));
    if(_S1330 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S1330 / x0_3;
        float _S1331 = 1.0f - t0_2;
        y_11 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S1331) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S1331);
    }
    else
    {
        float _S1332 = 1.0f - y0_3;
        float _S1333 = 1.0f - x0_3;
        float s1_6 = _S1332 / _S1333;
        float t1_2 = (_S1330 - x0_3) / _S1333;
        float _S1334 = 1.0f - t1_2;
        y_11 = y0_3 + _S1332 * (s1_6 * t1_2 * t1_2 + gc_3 * t1_2 * _S1334) / (s1_6 + (gc_3 + g1_3 - 2.0f * s1_6) * t1_2 * _S1334);
    }
    *&((&rgb_out_3)->z) = y_11;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1335)
{
    return (F32_exp2((_S1335)));
}

inline __device__ float2  s_primal_ctx_mul_1(Matrix<float, 2, 2>  _S1336, float2  _S1337)
{
    return mul_4(_S1336, _S1337);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1338, Matrix<float, 3, 3>  _S1339)
{
    return mul_1(_S1338, _S1339);
}

inline __device__ float s_primal_ctx_abs_0(float _S1340)
{
    return (F32_abs((_S1340)));
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1341, float3  _S1342, float3  _S1343)
{
    return clamp_1(_S1341, _S1342, _S1343);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1344, float _S1345, float _S1346)
{
    return lerp_0(_S1344, _S1345, _S1346);
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1347, float _S1348)
{
    _d_abs_0(_S1347, _S1348);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1349, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1350, float2  _S1351)
{
    _d_mul_2(_S1349, _S1350, _S1351);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1352, float _S1353)
{
    _d_exp2_0(_S1352, _S1353);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_13)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1354 = *dprgb_in_0;
    float3  _S1355 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1356 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1357 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1358 = {
        _S1357, _S1357, _S1357
    };
    float2  _S1359 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1360 = { _S1359, _S1359, _S1359, _S1359 };
    CRFPPISPChannelParams_0 _S1361 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1362 = {
        _S1361, _S1361, _S1361
    };
    PPISPParams_0 _S1363;
    (&_S1363)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S1363)->vignette_params_1 = _S1358;
    (&_S1363)->color_params_1 = _S1360;
    (&_S1363)->crf_params_1 = _S1362;
    (&(&_S1363)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1363)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1364 = dpparams_0->primal_0[int(3)];
    (&(&_S1363)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1365 = dpparams_0->primal_0[int(4)];
    (&(&_S1363)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1366 = dpparams_0->primal_0[int(5)];
    (&(&_S1363)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1363)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1363)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1367 = dpparams_0->primal_0[int(8)];
    (&(&_S1363)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1368 = dpparams_0->primal_0[int(9)];
    (&(&_S1363)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1369 = dpparams_0->primal_0[int(10)];
    (&(&_S1363)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1363)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1363)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1370 = dpparams_0->primal_0[int(13)];
    (&(&_S1363)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1371 = dpparams_0->primal_0[int(14)];
    (&(&_S1363)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1372 = dpparams_0->primal_0[int(15)];
    (&(&_S1363)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1363)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1363)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1363)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1363)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1363)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1363)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1363)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1363)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1373 = dpparams_0->primal_0[int(24)];
    (&(&_S1363)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1374 = dpparams_0->primal_0[int(25)];
    (&(&_S1363)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1375 = dpparams_0->primal_0[int(26)];
    (&(&_S1363)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1376 = dpparams_0->primal_0[int(27)];
    (&(&_S1363)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1377 = dpparams_0->primal_0[int(28)];
    (&(&_S1363)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1378 = dpparams_0->primal_0[int(29)];
    (&(&_S1363)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1379 = dpparams_0->primal_0[int(30)];
    (&(&_S1363)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1380 = dpparams_0->primal_0[int(31)];
    (&(&_S1363)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1381 = dpparams_0->primal_0[int(32)];
    (&(&_S1363)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1382 = dpparams_0->primal_0[int(33)];
    (&(&_S1363)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1383 = dpparams_0->primal_0[int(34)];
    (&(&_S1363)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1384 = dpparams_0->primal_0[int(35)];
    (&(&_S1363)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1385 = _S1363;
    float _S1386 = s_primal_ctx_exp2_0(_S1363.exposure_1);
    float3  _S1387 = make_float3 (_S1386);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S1386);
    float _S1388 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S1389 = (pix_coord_2.x - image_center_2.x) / _S1388;
    float _S1390 = (pix_coord_2.y - image_center_2.y) / _S1388;
    float dx_8 = _S1389 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S1390 - dpparams_0->primal_0[int(2)];
    float r2_13 = dx_8 * dx_8 + dy_6 * dy_6;
    float r4_6 = r2_13 * r2_13;
    float r6_0 = r4_6 * r2_13;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_13 + 1.0f;
    float _S1391 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1392 = rgb_out_4.x * _S1391;
    float3  _S1393 = rgb_out_4;
    *&((&_S1393)->x) = _S1392;
    float dx_9 = _S1389 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S1390 - dpparams_0->primal_0[int(7)];
    float r2_14 = dx_9 * dx_9 + dy_7 * dy_7;
    float r4_7 = r2_14 * r2_14;
    float r6_1 = r4_7 * r2_14;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_14 + 1.0f;
    float _S1394 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1393)->y) = rgb_out_4.y * _S1394;
    float dx_10 = _S1389 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S1390 - dpparams_0->primal_0[int(12)];
    float r2_15 = dx_10 * dx_10 + dy_8 * dy_8;
    float r4_8 = r2_15 * r2_15;
    float r6_2 = r4_8 * r2_15;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_15 + 1.0f;
    float _S1395 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1393)->z) = rgb_out_4.z * _S1395;
    PPISPParams_0 _S1396 = _S1363;
    float2  _S1397 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1363.color_params_1.b_0);
    float2  _S1398 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1363.color_params_1.r_0);
    float2  _S1399 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1363.color_params_1.g_0);
    float2  _S1400 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1363.color_params_1.n_0);
    float _S1401 = 0.3333333432674408f + _S1400.x;
    float _S1402 = 0.3333333432674408f + _S1400.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S1397.x, 1.0f + _S1398.x, _S1399.x, _S1397.y, _S1398.y, 1.0f + _S1399.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1402, 1.0f, 0.0f, - _S1401, - _S1402, _S1401, 0.0f);
    Matrix<float, 3, 3>  _S1403 = s_primal_ctx_mul_2(skew_0, T_2);
    float3  r0_2 = make_float3 (_S1403.rows[int(0)].x, _S1403.rows[int(0)].y, _S1403.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S1403.rows[int(1)].x, _S1403.rows[int(1)].y, _S1403.rows[int(1)].z);
    float3  r2_16 = make_float3 (_S1403.rows[int(2)].x, _S1403.rows[int(2)].y, _S1403.rows[int(2)].z);
    float3  _S1404 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S1405 = (s_primal_ctx_dot_0(_S1404, _S1404)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S1406;
    bool _S1407;
    if(_S1405)
    {
        float3  _S1408 = s_primal_ctx_cross_0(r0_2, r2_16);
        bool _S1409 = (s_primal_ctx_dot_0(_S1408, _S1408)) < 9.99999968265522539e-21f;
        if(_S1409)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_16);
        }
        else
        {
            lambda_v_6 = _S1408;
        }
        _S1407 = _S1409;
        _S1406 = _S1408;
    }
    else
    {
        lambda_v_6 = _S1404;
        _S1407 = false;
        _S1406 = _S1355;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S1410 = s_primal_ctx_mul_2(T_2, D_0);
    Matrix<float, 3, 3>  _S1411 = s_primal_ctx_mul_2(_S1410, S_inv_0);
    bool _S1412 = (s_primal_ctx_abs_0(_S1411.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S1413;
    float _S1414;
    if(_S1412)
    {
        float inv_s_0 = 1.0f / _S1411.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1415 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1416 = _S1411.rows[int(2)].z * _S1411.rows[int(2)].z;
        H_4 = _S1411 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1413 = _S1415;
        _S1414 = _S1416;
    }
    else
    {
        H_4 = _S1411;
        _S1413 = _S1356;
        _S1414 = 0.0f;
    }
    float _S1417 = _S1393.x;
    float _S1418 = _S1393.y;
    float intensity_2 = _S1417 + _S1418 + _S1393.z;
    float3  rgi_in_0 = make_float3 (_S1417, _S1418, intensity_2);
    float3  _S1419 = s_primal_ctx_mul_0(H_4, rgi_in_0);
    float _S1420 = _S1419.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S1420;
    float3  _S1421 = make_float3 (norm_factor_0);
    float _S1422 = _S1420 * _S1420;
    float3  rgi_out_4 = _S1419 * make_float3 (norm_factor_0);
    float _S1423 = rgi_out_4.x;
    float _S1424 = rgi_out_4.y;
    float3  _S1425 = make_float3 (_S1423, _S1424, rgi_out_4.z - _S1423 - _S1424);
    float3  _S1426 = make_float3 (0.0f);
    float3  _S1427 = make_float3 (1.0f);
    float3  _S1428 = s_primal_ctx_clamp_1(_S1425, _S1426, _S1427);
    float _S1429 = _S1428.x;
    float _S1430 = 1.0f + s_primal_ctx_exp_0(_S1373);
    float _S1431 = 0.30000001192092896f + s_primal_ctx_log_0(_S1430);
    float _S1432 = 1.0f + s_primal_ctx_exp_0(_S1374);
    float _S1433 = 0.30000001192092896f + s_primal_ctx_log_0(_S1432);
    float _S1434 = 1.0f + s_primal_ctx_exp_0(_S1375);
    float _S1435 = 0.10000000149011612f + s_primal_ctx_log_0(_S1434);
    float _S1436 = - _S1376;
    float _S1437 = 1.0f + s_primal_ctx_exp_0(_S1436);
    float _S1438 = 1.0f / _S1437;
    float _S1439 = _S1437 * _S1437;
    float _S1440 = s_primal_ctx_lerp_0(_S1431, _S1433, _S1438);
    float _S1441 = _S1433 * _S1438;
    float a_4 = _S1441 / _S1440;
    float _S1442 = _S1440 * _S1440;
    float b_5 = 1.0f - a_4;
    bool _S1443 = _S1429 <= _S1438;
    float y_12;
    float _S1444;
    float _S1445;
    float _S1446;
    float _S1447;
    float _S1448;
    float _S1449;
    float _S1450;
    float _S1451;
    if(_S1443)
    {
        float _S1452 = _S1429 / _S1438;
        float _S1453 = _S1438 * _S1438;
        float _S1454 = s_primal_ctx_pow_0(_S1452, _S1431);
        y_12 = a_4 * _S1454;
        _S1444 = _S1454;
        _S1445 = _S1452;
        _S1446 = _S1453;
        _S1447 = 0.0f;
        _S1448 = 0.0f;
        _S1449 = 0.0f;
        _S1450 = 0.0f;
        _S1451 = 0.0f;
    }
    else
    {
        float _S1455 = 1.0f - _S1429;
        float _S1456 = 1.0f - _S1438;
        float _S1457 = _S1455 / _S1456;
        float _S1458 = _S1456 * _S1456;
        float _S1459 = s_primal_ctx_pow_0(_S1457, _S1433);
        y_12 = 1.0f - b_5 * _S1459;
        _S1444 = 0.0f;
        _S1445 = 0.0f;
        _S1446 = 0.0f;
        _S1447 = _S1459;
        _S1448 = _S1457;
        _S1449 = _S1458;
        _S1450 = _S1455;
        _S1451 = _S1456;
    }
    float _S1460 = (F32_max((0.0f), (y_12)));
    float _S1461 = _S1428.y;
    float _S1462 = 1.0f + s_primal_ctx_exp_0(_S1377);
    float _S1463 = 0.30000001192092896f + s_primal_ctx_log_0(_S1462);
    float _S1464 = 1.0f + s_primal_ctx_exp_0(_S1378);
    float _S1465 = 0.30000001192092896f + s_primal_ctx_log_0(_S1464);
    float _S1466 = 1.0f + s_primal_ctx_exp_0(_S1379);
    float _S1467 = 0.10000000149011612f + s_primal_ctx_log_0(_S1466);
    float _S1468 = - _S1380;
    float _S1469 = 1.0f + s_primal_ctx_exp_0(_S1468);
    float _S1470 = 1.0f / _S1469;
    float _S1471 = _S1469 * _S1469;
    float _S1472 = s_primal_ctx_lerp_0(_S1463, _S1465, _S1470);
    float _S1473 = _S1465 * _S1470;
    float a_5 = _S1473 / _S1472;
    float _S1474 = _S1472 * _S1472;
    float b_6 = 1.0f - a_5;
    bool _S1475 = _S1461 <= _S1470;
    float y_13;
    float _S1476;
    float _S1477;
    float _S1478;
    float _S1479;
    float _S1480;
    float _S1481;
    float _S1482;
    float _S1483;
    if(_S1475)
    {
        float _S1484 = _S1461 / _S1470;
        float _S1485 = _S1470 * _S1470;
        float _S1486 = s_primal_ctx_pow_0(_S1484, _S1463);
        y_13 = a_5 * _S1486;
        _S1476 = _S1486;
        _S1477 = _S1484;
        _S1478 = _S1485;
        _S1479 = 0.0f;
        _S1480 = 0.0f;
        _S1481 = 0.0f;
        _S1482 = 0.0f;
        _S1483 = 0.0f;
    }
    else
    {
        float _S1487 = 1.0f - _S1461;
        float _S1488 = 1.0f - _S1470;
        float _S1489 = _S1487 / _S1488;
        float _S1490 = _S1488 * _S1488;
        float _S1491 = s_primal_ctx_pow_0(_S1489, _S1465);
        y_13 = 1.0f - b_6 * _S1491;
        _S1476 = 0.0f;
        _S1477 = 0.0f;
        _S1478 = 0.0f;
        _S1479 = _S1491;
        _S1480 = _S1489;
        _S1481 = _S1490;
        _S1482 = _S1487;
        _S1483 = _S1488;
    }
    float _S1492 = (F32_max((0.0f), (y_13)));
    float _S1493 = _S1428.z;
    float _S1494 = 1.0f + s_primal_ctx_exp_0(_S1381);
    float _S1495 = 0.30000001192092896f + s_primal_ctx_log_0(_S1494);
    float _S1496 = 1.0f + s_primal_ctx_exp_0(_S1382);
    float _S1497 = 0.30000001192092896f + s_primal_ctx_log_0(_S1496);
    float _S1498 = 1.0f + s_primal_ctx_exp_0(_S1383);
    float _S1499 = 0.10000000149011612f + s_primal_ctx_log_0(_S1498);
    float _S1500 = - _S1384;
    float _S1501 = 1.0f + s_primal_ctx_exp_0(_S1500);
    float _S1502 = 1.0f / _S1501;
    float _S1503 = _S1501 * _S1501;
    float _S1504 = s_primal_ctx_lerp_0(_S1495, _S1497, _S1502);
    float _S1505 = _S1497 * _S1502;
    float a_6 = _S1505 / _S1504;
    float _S1506 = _S1504 * _S1504;
    float b_7 = 1.0f - a_6;
    bool _S1507 = _S1493 <= _S1502;
    float y_14;
    float _S1508;
    float _S1509;
    float _S1510;
    float _S1511;
    float _S1512;
    float _S1513;
    float _S1514;
    float _S1515;
    if(_S1507)
    {
        float _S1516 = _S1493 / _S1502;
        float _S1517 = _S1502 * _S1502;
        float _S1518 = s_primal_ctx_pow_0(_S1516, _S1495);
        y_14 = a_6 * _S1518;
        _S1508 = _S1518;
        _S1509 = _S1516;
        _S1510 = _S1517;
        _S1511 = 0.0f;
        _S1512 = 0.0f;
        _S1513 = 0.0f;
        _S1514 = 0.0f;
        _S1515 = 0.0f;
    }
    else
    {
        float _S1519 = 1.0f - _S1493;
        float _S1520 = 1.0f - _S1502;
        float _S1521 = _S1519 / _S1520;
        float _S1522 = _S1520 * _S1520;
        float _S1523 = s_primal_ctx_pow_0(_S1521, _S1497);
        y_14 = 1.0f - b_7 * _S1523;
        _S1508 = 0.0f;
        _S1509 = 0.0f;
        _S1510 = 0.0f;
        _S1511 = _S1523;
        _S1512 = _S1521;
        _S1513 = _S1522;
        _S1514 = _S1519;
        _S1515 = _S1520;
    }
    float _S1524 = (F32_max((0.0f), (y_14)));
    DiffPair_float_0 _S1525;
    (&_S1525)->primal_0 = _S1524;
    (&_S1525)->differential_0 = 0.0f;
    DiffPair_float_0 _S1526;
    (&_S1526)->primal_0 = _S1499;
    (&_S1526)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1525, &_S1526, _s_dOut_13.z);
    DiffPair_float_0 _S1527 = _S1526;
    DiffPair_float_0 _S1528;
    (&_S1528)->primal_0 = 0.0f;
    (&_S1528)->differential_0 = 0.0f;
    DiffPair_float_0 _S1529;
    (&_S1529)->primal_0 = y_14;
    (&_S1529)->differential_0 = 0.0f;
    _d_max_0(&_S1528, &_S1529, _S1525.differential_0);
    DiffPair_float_0 _S1530 = _S1529;
    if(_S1507)
    {
        float _S1531 = a_6 * _S1530.differential_0;
        float _S1532 = _S1508 * _S1530.differential_0;
        DiffPair_float_0 _S1533;
        (&_S1533)->primal_0 = _S1509;
        (&_S1533)->differential_0 = 0.0f;
        DiffPair_float_0 _S1534;
        (&_S1534)->primal_0 = _S1495;
        (&_S1534)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1533, &_S1534, _S1531);
        float _S1535 = _S1533.differential_0 / _S1510;
        float _S1536 = _S1493 * - _S1535;
        float _S1537 = _S1502 * _S1535;
        y_14 = 0.0f;
        _S1508 = _S1532;
        _S1509 = _S1536;
        _S1510 = 0.0f;
        _S1511 = _S1534.differential_0;
        _S1512 = _S1537;
    }
    else
    {
        float _S1538 = - _S1530.differential_0;
        float _S1539 = b_7 * _S1538;
        float _S1540 = _S1511 * _S1538;
        DiffPair_float_0 _S1541;
        (&_S1541)->primal_0 = _S1512;
        (&_S1541)->differential_0 = 0.0f;
        DiffPair_float_0 _S1542;
        (&_S1542)->primal_0 = _S1497;
        (&_S1542)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1541, &_S1542, _S1539);
        float _S1543 = _S1541.differential_0 / _S1513;
        float _S1544 = - (_S1514 * - _S1543);
        float _S1545 = - (_S1515 * _S1543);
        y_14 = _S1540;
        _S1508 = 0.0f;
        _S1509 = _S1544;
        _S1510 = _S1542.differential_0;
        _S1511 = 0.0f;
        _S1512 = _S1545;
    }
    float _S1546 = (- y_14 + _S1508) / _S1506;
    float _S1547 = _S1505 * - _S1546;
    float _S1548 = _S1504 * _S1546;
    float _S1549 = _S1497 * _S1548;
    float _S1550 = _S1502 * _S1548;
    DiffPair_float_0 _S1551;
    (&_S1551)->primal_0 = _S1495;
    (&_S1551)->differential_0 = 0.0f;
    DiffPair_float_0 _S1552;
    (&_S1552)->primal_0 = _S1497;
    (&_S1552)->differential_0 = 0.0f;
    DiffPair_float_0 _S1553;
    (&_S1553)->primal_0 = _S1502;
    (&_S1553)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1551, &_S1552, &_S1553, _S1547);
    float _S1554 = - ((_S1549 + _S1553.differential_0 + _S1509) / _S1503);
    DiffPair_float_0 _S1555;
    (&_S1555)->primal_0 = _S1500;
    (&_S1555)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1555, _S1554);
    float _S1556 = - _S1555.differential_0;
    DiffPair_float_0 _S1557;
    (&_S1557)->primal_0 = _S1498;
    (&_S1557)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1557, _S1527.differential_0);
    DiffPair_float_0 _S1558;
    (&_S1558)->primal_0 = _S1383;
    (&_S1558)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1558, _S1557.differential_0);
    DiffPair_float_0 _S1559 = _S1558;
    float _S1560 = _S1550 + _S1552.differential_0 + _S1510;
    DiffPair_float_0 _S1561;
    (&_S1561)->primal_0 = _S1496;
    (&_S1561)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1561, _S1560);
    DiffPair_float_0 _S1562;
    (&_S1562)->primal_0 = _S1382;
    (&_S1562)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1562, _S1561.differential_0);
    DiffPair_float_0 _S1563 = _S1562;
    float _S1564 = _S1551.differential_0 + _S1511;
    DiffPair_float_0 _S1565;
    (&_S1565)->primal_0 = _S1494;
    (&_S1565)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1565, _S1564);
    DiffPair_float_0 _S1566;
    (&_S1566)->primal_0 = _S1381;
    (&_S1566)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1566, _S1565.differential_0);
    DiffPair_float_0 _S1567 = _S1566;
    float3  _S1568 = make_float3 (0.0f, 0.0f, _S1512);
    DiffPair_float_0 _S1569;
    (&_S1569)->primal_0 = _S1492;
    (&_S1569)->differential_0 = 0.0f;
    DiffPair_float_0 _S1570;
    (&_S1570)->primal_0 = _S1467;
    (&_S1570)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1569, &_S1570, _s_dOut_13.y);
    DiffPair_float_0 _S1571 = _S1570;
    DiffPair_float_0 _S1572;
    (&_S1572)->primal_0 = 0.0f;
    (&_S1572)->differential_0 = 0.0f;
    DiffPair_float_0 _S1573;
    (&_S1573)->primal_0 = y_13;
    (&_S1573)->differential_0 = 0.0f;
    _d_max_0(&_S1572, &_S1573, _S1569.differential_0);
    DiffPair_float_0 _S1574 = _S1573;
    if(_S1475)
    {
        float _S1575 = a_5 * _S1574.differential_0;
        float _S1576 = _S1476 * _S1574.differential_0;
        DiffPair_float_0 _S1577;
        (&_S1577)->primal_0 = _S1477;
        (&_S1577)->differential_0 = 0.0f;
        DiffPair_float_0 _S1578;
        (&_S1578)->primal_0 = _S1463;
        (&_S1578)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1577, &_S1578, _S1575);
        float _S1579 = _S1577.differential_0 / _S1478;
        float _S1580 = _S1461 * - _S1579;
        float _S1581 = _S1470 * _S1579;
        y_13 = 0.0f;
        _S1476 = _S1576;
        _S1477 = _S1580;
        _S1478 = 0.0f;
        _S1479 = _S1578.differential_0;
        _S1480 = _S1581;
    }
    else
    {
        float _S1582 = - _S1574.differential_0;
        float _S1583 = b_6 * _S1582;
        float _S1584 = _S1479 * _S1582;
        DiffPair_float_0 _S1585;
        (&_S1585)->primal_0 = _S1480;
        (&_S1585)->differential_0 = 0.0f;
        DiffPair_float_0 _S1586;
        (&_S1586)->primal_0 = _S1465;
        (&_S1586)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1585, &_S1586, _S1583);
        float _S1587 = _S1585.differential_0 / _S1481;
        float _S1588 = - (_S1482 * - _S1587);
        float _S1589 = - (_S1483 * _S1587);
        y_13 = _S1584;
        _S1476 = 0.0f;
        _S1477 = _S1588;
        _S1478 = _S1586.differential_0;
        _S1479 = 0.0f;
        _S1480 = _S1589;
    }
    float _S1590 = (- y_13 + _S1476) / _S1474;
    float _S1591 = _S1473 * - _S1590;
    float _S1592 = _S1472 * _S1590;
    float _S1593 = _S1465 * _S1592;
    float _S1594 = _S1470 * _S1592;
    DiffPair_float_0 _S1595;
    (&_S1595)->primal_0 = _S1463;
    (&_S1595)->differential_0 = 0.0f;
    DiffPair_float_0 _S1596;
    (&_S1596)->primal_0 = _S1465;
    (&_S1596)->differential_0 = 0.0f;
    DiffPair_float_0 _S1597;
    (&_S1597)->primal_0 = _S1470;
    (&_S1597)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1595, &_S1596, &_S1597, _S1591);
    float _S1598 = - ((_S1593 + _S1597.differential_0 + _S1477) / _S1471);
    DiffPair_float_0 _S1599;
    (&_S1599)->primal_0 = _S1468;
    (&_S1599)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1599, _S1598);
    float _S1600 = - _S1599.differential_0;
    DiffPair_float_0 _S1601;
    (&_S1601)->primal_0 = _S1466;
    (&_S1601)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1601, _S1571.differential_0);
    DiffPair_float_0 _S1602;
    (&_S1602)->primal_0 = _S1379;
    (&_S1602)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1602, _S1601.differential_0);
    DiffPair_float_0 _S1603 = _S1602;
    float _S1604 = _S1594 + _S1596.differential_0 + _S1478;
    DiffPair_float_0 _S1605;
    (&_S1605)->primal_0 = _S1464;
    (&_S1605)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1605, _S1604);
    DiffPair_float_0 _S1606;
    (&_S1606)->primal_0 = _S1378;
    (&_S1606)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1606, _S1605.differential_0);
    DiffPair_float_0 _S1607 = _S1606;
    float _S1608 = _S1595.differential_0 + _S1479;
    DiffPair_float_0 _S1609;
    (&_S1609)->primal_0 = _S1462;
    (&_S1609)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1609, _S1608);
    DiffPair_float_0 _S1610;
    (&_S1610)->primal_0 = _S1377;
    (&_S1610)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1610, _S1609.differential_0);
    DiffPair_float_0 _S1611 = _S1610;
    float3  _S1612 = _S1568 + make_float3 (0.0f, _S1480, 0.0f);
    DiffPair_float_0 _S1613;
    (&_S1613)->primal_0 = _S1460;
    (&_S1613)->differential_0 = 0.0f;
    DiffPair_float_0 _S1614;
    (&_S1614)->primal_0 = _S1435;
    (&_S1614)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1613, &_S1614, _s_dOut_13.x);
    DiffPair_float_0 _S1615 = _S1614;
    DiffPair_float_0 _S1616;
    (&_S1616)->primal_0 = 0.0f;
    (&_S1616)->differential_0 = 0.0f;
    DiffPair_float_0 _S1617;
    (&_S1617)->primal_0 = y_12;
    (&_S1617)->differential_0 = 0.0f;
    _d_max_0(&_S1616, &_S1617, _S1613.differential_0);
    DiffPair_float_0 _S1618 = _S1617;
    if(_S1443)
    {
        float _S1619 = a_4 * _S1618.differential_0;
        float _S1620 = _S1444 * _S1618.differential_0;
        DiffPair_float_0 _S1621;
        (&_S1621)->primal_0 = _S1445;
        (&_S1621)->differential_0 = 0.0f;
        DiffPair_float_0 _S1622;
        (&_S1622)->primal_0 = _S1431;
        (&_S1622)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1621, &_S1622, _S1619);
        float _S1623 = _S1621.differential_0 / _S1446;
        float _S1624 = _S1429 * - _S1623;
        float _S1625 = _S1438 * _S1623;
        y_12 = 0.0f;
        _S1444 = _S1620;
        _S1445 = _S1624;
        _S1446 = 0.0f;
        _S1447 = _S1622.differential_0;
        _S1448 = _S1625;
    }
    else
    {
        float _S1626 = - _S1618.differential_0;
        float _S1627 = b_5 * _S1626;
        float _S1628 = _S1447 * _S1626;
        DiffPair_float_0 _S1629;
        (&_S1629)->primal_0 = _S1448;
        (&_S1629)->differential_0 = 0.0f;
        DiffPair_float_0 _S1630;
        (&_S1630)->primal_0 = _S1433;
        (&_S1630)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1629, &_S1630, _S1627);
        float _S1631 = _S1629.differential_0 / _S1449;
        float _S1632 = - (_S1450 * - _S1631);
        float _S1633 = - (_S1451 * _S1631);
        y_12 = _S1628;
        _S1444 = 0.0f;
        _S1445 = _S1632;
        _S1446 = _S1630.differential_0;
        _S1447 = 0.0f;
        _S1448 = _S1633;
    }
    float _S1634 = (- y_12 + _S1444) / _S1442;
    float _S1635 = _S1441 * - _S1634;
    float _S1636 = _S1440 * _S1634;
    float _S1637 = _S1433 * _S1636;
    float _S1638 = _S1438 * _S1636;
    DiffPair_float_0 _S1639;
    (&_S1639)->primal_0 = _S1431;
    (&_S1639)->differential_0 = 0.0f;
    DiffPair_float_0 _S1640;
    (&_S1640)->primal_0 = _S1433;
    (&_S1640)->differential_0 = 0.0f;
    DiffPair_float_0 _S1641;
    (&_S1641)->primal_0 = _S1438;
    (&_S1641)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1639, &_S1640, &_S1641, _S1635);
    float _S1642 = - ((_S1637 + _S1641.differential_0 + _S1445) / _S1439);
    DiffPair_float_0 _S1643;
    (&_S1643)->primal_0 = _S1436;
    (&_S1643)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1643, _S1642);
    float _S1644 = - _S1643.differential_0;
    DiffPair_float_0 _S1645;
    (&_S1645)->primal_0 = _S1434;
    (&_S1645)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1645, _S1615.differential_0);
    DiffPair_float_0 _S1646;
    (&_S1646)->primal_0 = _S1375;
    (&_S1646)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1646, _S1645.differential_0);
    DiffPair_float_0 _S1647 = _S1646;
    float _S1648 = _S1638 + _S1640.differential_0 + _S1446;
    DiffPair_float_0 _S1649;
    (&_S1649)->primal_0 = _S1432;
    (&_S1649)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1649, _S1648);
    DiffPair_float_0 _S1650;
    (&_S1650)->primal_0 = _S1374;
    (&_S1650)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1650, _S1649.differential_0);
    DiffPair_float_0 _S1651 = _S1650;
    float _S1652 = _S1639.differential_0 + _S1447;
    DiffPair_float_0 _S1653;
    (&_S1653)->primal_0 = _S1430;
    (&_S1653)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1653, _S1652);
    DiffPair_float_0 _S1654;
    (&_S1654)->primal_0 = _S1373;
    (&_S1654)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1654, _S1653.differential_0);
    DiffPair_float_0 _S1655 = _S1654;
    float3  _S1656 = _S1612 + make_float3 (_S1448, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1657;
    (&_S1657)->primal_0 = _S1425;
    (&_S1657)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1658;
    (&_S1658)->primal_0 = _S1426;
    (&_S1658)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1659;
    (&_S1659)->primal_0 = _S1427;
    (&_S1659)->differential_0 = _S1355;
    s_bwd_prop_clamp_1(&_S1657, &_S1658, &_S1659, _S1656);
    float _S1660 = - _S1657.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1657.differential_0.x + _S1660, _S1657.differential_0.y + _S1660, _S1657.differential_0.z);
    float3  _S1661 = _S1419 * s_diff_rgi_out_T_0;
    float _S1662 = (_S1661.x + _S1661.y + _S1661.z) / _S1422;
    float _S1663 = _S1420 * _S1662;
    float3  _S1664 = _S1421 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_2 * - _S1662);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1665;
    (&_S1665)->primal_0 = H_4;
    (&_S1665)->differential_0 = _S1356;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1666;
    (&_S1666)->primal_0 = rgi_in_0;
    (&_S1666)->differential_0 = _S1355;
    s_bwd_prop_mul_2(&_S1665, &_S1666, _S1664);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1667 = _S1665;
    float _S1668 = _S1663 + _S1666.differential_0.z;
    float _S1669 = _S1666.differential_0.y + _S1668;
    float _S1670 = _S1666.differential_0.x + _S1668;
    float3  _S1671 = make_float3 (_S1670, _S1669, _S1668);
    if(_S1412)
    {
        Matrix<float, 3, 3>  _S1672 = _S1411 * _S1667.differential_0;
        Matrix<float, 3, 3>  _S1673 = _S1413 * _S1667.differential_0;
        _S1414 = - ((_S1672.rows[int(0)].x + _S1672.rows[int(0)].y + _S1672.rows[int(0)].z + _S1672.rows[int(1)].x + _S1672.rows[int(1)].y + _S1672.rows[int(1)].z + _S1672.rows[int(2)].x + _S1672.rows[int(2)].y + _S1672.rows[int(2)].z) / _S1414);
        H_4 = _S1673;
    }
    else
    {
        _S1414 = 0.0f;
        H_4 = _S1667.differential_0;
    }
    DiffPair_float_0 _S1674;
    (&_S1674)->primal_0 = _S1411.rows[int(2)].z;
    (&_S1674)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1674, 0.0f);
    float _S1675 = _S1674.differential_0 + _S1414;
    float3  _S1676 = _S1355;
    *&((&_S1676)->z) = _S1675;
    Matrix<float, 3, 3>  _S1677 = _S1356;
    _S1677[int(2)] = _S1676;
    Matrix<float, 3, 3>  _S1678 = H_4 + _S1677;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1679;
    (&_S1679)->primal_0 = _S1410;
    (&_S1679)->differential_0 = _S1356;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1680;
    (&_S1680)->primal_0 = S_inv_0;
    (&_S1680)->differential_0 = _S1356;
    s_bwd_prop_mul_1(&_S1679, &_S1680, _S1678);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1681;
    (&_S1681)->primal_0 = T_2;
    (&_S1681)->differential_0 = _S1356;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1682;
    (&_S1682)->primal_0 = D_0;
    (&_S1682)->differential_0 = _S1356;
    s_bwd_prop_mul_1(&_S1681, &_S1682, _S1679.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1683 = _S1681;
    float3  _S1684 = make_float3 (_S1682.differential_0.rows[int(0)].x, _S1682.differential_0.rows[int(1)].y, _S1682.differential_0.rows[int(2)].z);
    float3  _S1685;
    if(_S1405)
    {
        if(_S1407)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1686;
            (&_S1686)->primal_0 = r1_2;
            (&_S1686)->differential_0 = _S1355;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1687;
            (&_S1687)->primal_0 = r2_16;
            (&_S1687)->differential_0 = _S1355;
            s_bwd_prop_cross_0(&_S1686, &_S1687, _S1684);
            _S1393 = _S1355;
            lambda_v_6 = _S1687.differential_0;
            _S1685 = _S1686.differential_0;
        }
        else
        {
            _S1393 = _S1684;
            lambda_v_6 = _S1355;
            _S1685 = _S1355;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1688;
        (&_S1688)->primal_0 = _S1406;
        (&_S1688)->differential_0 = _S1355;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1689;
        (&_S1689)->primal_0 = _S1406;
        (&_S1689)->differential_0 = _S1355;
        s_bwd_prop_dot_0(&_S1688, &_S1689, 0.0f);
        float3  _S1690 = _S1689.differential_0 + _S1688.differential_0 + _S1393;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1691;
        (&_S1691)->primal_0 = r0_2;
        (&_S1691)->differential_0 = _S1355;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1692;
        (&_S1692)->primal_0 = r2_16;
        (&_S1692)->differential_0 = _S1355;
        s_bwd_prop_cross_0(&_S1691, &_S1692, _S1690);
        float3  _S1693 = _S1692.differential_0 + lambda_v_6;
        _S1393 = _S1355;
        lambda_v_6 = _S1693;
        _S1406 = _S1685;
        _S1685 = _S1691.differential_0;
    }
    else
    {
        _S1393 = _S1684;
        lambda_v_6 = _S1355;
        _S1406 = _S1355;
        _S1685 = _S1355;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1694;
    (&_S1694)->primal_0 = _S1404;
    (&_S1694)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1695;
    (&_S1695)->primal_0 = _S1404;
    (&_S1695)->differential_0 = _S1355;
    s_bwd_prop_dot_0(&_S1694, &_S1695, 0.0f);
    float3  _S1696 = _S1695.differential_0 + _S1694.differential_0 + _S1393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1697;
    (&_S1697)->primal_0 = r0_2;
    (&_S1697)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1698;
    (&_S1698)->primal_0 = r1_2;
    (&_S1698)->differential_0 = _S1355;
    s_bwd_prop_cross_0(&_S1697, &_S1698, _S1696);
    float3  _S1699 = _S1355;
    *&((&_S1699)->z) = lambda_v_6.z;
    *&((&_S1699)->y) = lambda_v_6.y;
    *&((&_S1699)->x) = lambda_v_6.x;
    float3  _S1700 = _S1698.differential_0 + _S1406;
    float3  _S1701 = _S1355;
    *&((&_S1701)->z) = _S1700.z;
    *&((&_S1701)->y) = _S1700.y;
    *&((&_S1701)->x) = _S1700.x;
    float3  _S1702 = _S1697.differential_0 + _S1685;
    float3  _S1703 = _S1355;
    *&((&_S1703)->z) = _S1702.z;
    *&((&_S1703)->y) = _S1702.y;
    *&((&_S1703)->x) = _S1702.x;
    Matrix<float, 3, 3>  _S1704 = _S1356;
    _S1704[int(2)] = _S1699;
    _S1704[int(1)] = _S1701;
    _S1704[int(0)] = _S1703;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1705;
    (&_S1705)->primal_0 = skew_0;
    (&_S1705)->differential_0 = _S1356;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1706;
    (&_S1706)->primal_0 = T_2;
    (&_S1706)->differential_0 = _S1356;
    s_bwd_prop_mul_1(&_S1705, &_S1706, _S1704);
    Matrix<float, 3, 3>  _S1707 = _S1706.differential_0 + _S1683.differential_0;
    float2  _S1708 = make_float2 (_S1705.differential_0.rows[int(2)].y + - _S1705.differential_0.rows[int(1)].z, _S1705.differential_0.rows[int(0)].z + - _S1705.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1709 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1710;
    (&_S1710)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1710)->differential_0 = _S1709;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1711;
    (&_S1711)->primal_0 = _S1396.color_params_1.n_0;
    (&_S1711)->differential_0 = _S1359;
    s_bwd_prop_mul_3(&_S1710, &_S1711, _S1708);
    float2  _S1712 = make_float2 (_S1707.rows[int(0)].z, _S1707.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1713;
    (&_S1713)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1713)->differential_0 = _S1709;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1714;
    (&_S1714)->primal_0 = _S1396.color_params_1.g_0;
    (&_S1714)->differential_0 = _S1359;
    s_bwd_prop_mul_3(&_S1713, &_S1714, _S1712);
    float2  _S1715 = make_float2 (_S1707.rows[int(0)].y, _S1707.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1716;
    (&_S1716)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1716)->differential_0 = _S1709;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1717;
    (&_S1717)->primal_0 = _S1396.color_params_1.r_0;
    (&_S1717)->differential_0 = _S1359;
    s_bwd_prop_mul_3(&_S1716, &_S1717, _S1715);
    float2  _S1718 = make_float2 (_S1707.rows[int(0)].x, _S1707.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1719;
    (&_S1719)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1719)->differential_0 = _S1709;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1720;
    (&_S1720)->primal_0 = _S1396.color_params_1.b_0;
    (&_S1720)->differential_0 = _S1359;
    s_bwd_prop_mul_3(&_S1719, &_S1720, _S1718);
    ColorPPISPParams_0 _S1721 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1721)->n_0 = _S1711.differential_0;
    (&_S1721)->g_0 = _S1714.differential_0;
    (&_S1721)->r_0 = _S1717.differential_0;
    (&_S1721)->b_0 = _S1720.differential_0;
    _S1393 = _S1671;
    *&((&_S1393)->z) = 0.0f;
    float _S1722 = rgb_out_4.z * _S1668;
    float _S1723 = _S1395 * _S1668;
    DiffPair_float_0 _S1724;
    (&_S1724)->primal_0 = falloff_2;
    (&_S1724)->differential_0 = 0.0f;
    DiffPair_float_0 _S1725;
    (&_S1725)->primal_0 = 0.0f;
    (&_S1725)->differential_0 = 0.0f;
    DiffPair_float_0 _S1726;
    (&_S1726)->primal_0 = 1.0f;
    (&_S1726)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1724, &_S1725, &_S1726, _S1722);
    float _S1727 = r2_15 * _S1724.differential_0;
    float _S1728 = r4_8 * _S1724.differential_0;
    float s_diff_r6_T_0 = _S1372 * _S1724.differential_0;
    float _S1729 = r6_2 * _S1724.differential_0;
    float _S1730 = r2_15 * (_S1371 * _S1724.differential_0 + r2_15 * s_diff_r6_T_0);
    float _S1731 = _S1370 * _S1724.differential_0 + r4_8 * s_diff_r6_T_0 + _S1730 + _S1730;
    float _S1732 = dy_8 * _S1731;
    float _S1733 = dx_10 * _S1731;
    float _S1734 = - (_S1732 + _S1732);
    float _S1735 = - (_S1733 + _S1733);
    *&((&_S1393)->y) = 0.0f;
    float _S1736 = rgb_out_4.y * _S1669;
    float _S1737 = _S1394 * _S1669;
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = falloff_1;
    (&_S1738)->differential_0 = 0.0f;
    DiffPair_float_0 _S1739;
    (&_S1739)->primal_0 = 0.0f;
    (&_S1739)->differential_0 = 0.0f;
    DiffPair_float_0 _S1740;
    (&_S1740)->primal_0 = 1.0f;
    (&_S1740)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1738, &_S1739, &_S1740, _S1736);
    float _S1741 = r2_14 * _S1738.differential_0;
    float _S1742 = r4_7 * _S1738.differential_0;
    float s_diff_r6_T_1 = _S1369 * _S1738.differential_0;
    float _S1743 = r6_1 * _S1738.differential_0;
    float _S1744 = r2_14 * (_S1368 * _S1738.differential_0 + r2_14 * s_diff_r6_T_1);
    float _S1745 = _S1367 * _S1738.differential_0 + r4_7 * s_diff_r6_T_1 + _S1744 + _S1744;
    float _S1746 = dy_7 * _S1745;
    float _S1747 = dx_9 * _S1745;
    float _S1748 = - (_S1746 + _S1746);
    float _S1749 = - (_S1747 + _S1747);
    *&((&_S1393)->x) = 0.0f;
    float _S1750 = rgb_out_4.x * _S1670;
    float _S1751 = _S1391 * _S1670;
    DiffPair_float_0 _S1752;
    (&_S1752)->primal_0 = falloff_0;
    (&_S1752)->differential_0 = 0.0f;
    DiffPair_float_0 _S1753;
    (&_S1753)->primal_0 = 0.0f;
    (&_S1753)->differential_0 = 0.0f;
    DiffPair_float_0 _S1754;
    (&_S1754)->primal_0 = 1.0f;
    (&_S1754)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1752, &_S1753, &_S1754, _S1750);
    float _S1755 = r2_13 * _S1752.differential_0;
    float _S1756 = r4_6 * _S1752.differential_0;
    float s_diff_r6_T_2 = _S1366 * _S1752.differential_0;
    float _S1757 = r6_0 * _S1752.differential_0;
    float _S1758 = r2_13 * (_S1365 * _S1752.differential_0 + r2_13 * s_diff_r6_T_2);
    float _S1759 = _S1364 * _S1752.differential_0 + r4_6 * s_diff_r6_T_2 + _S1758 + _S1758;
    float _S1760 = dy_6 * _S1759;
    float _S1761 = dx_8 * _S1759;
    float _S1762 = - (_S1760 + _S1760);
    float _S1763 = - (_S1761 + _S1761);
    float3  _S1764 = _S1355;
    *&((&_S1764)->z) = _S1723;
    *&((&_S1764)->y) = _S1737;
    *&((&_S1764)->x) = _S1751;
    float3  _S1765 = _S1393 + _S1764;
    float3  _S1766 = _S1354.primal_0 * _S1765;
    float3  _S1767 = _S1387 * _S1765;
    float _S1768 = _S1766.x + _S1766.y + _S1766.z;
    DiffPair_float_0 _S1769;
    (&_S1769)->primal_0 = _S1385.exposure_1;
    (&_S1769)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1769, _S1768);
    PPISPParams_0 _S1770 = PPISPParams_x24_syn_dzero_0();
    (&_S1770)->color_params_1 = _S1721;
    (&_S1770)->exposure_1 = _S1769.differential_0;
    _S1363 = _S1770;
    (&(&_S1363)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1771 = _S1770.crf_params_1[int(2)].center_0 + _S1556;
    (&(&_S1363)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1772 = _S1770.crf_params_1[int(2)].gamma_0 + _S1559.differential_0;
    (&(&_S1363)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1773 = _S1770.crf_params_1[int(2)].shoulder_0 + _S1563.differential_0;
    (&(&_S1363)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1774 = _S1770.crf_params_1[int(2)].toe_0 + _S1567.differential_0;
    (&(&_S1363)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1775 = _S1770.crf_params_1[int(1)].center_0 + _S1600;
    (&(&_S1363)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1776 = _S1770.crf_params_1[int(1)].gamma_0 + _S1603.differential_0;
    (&(&_S1363)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1777 = _S1770.crf_params_1[int(1)].shoulder_0 + _S1607.differential_0;
    (&(&_S1363)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1778 = _S1770.crf_params_1[int(1)].toe_0 + _S1611.differential_0;
    (&(&_S1363)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1779 = _S1770.crf_params_1[int(0)].center_0 + _S1644;
    (&(&_S1363)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1780 = _S1770.crf_params_1[int(0)].gamma_0 + _S1647.differential_0;
    (&(&_S1363)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1781 = _S1770.crf_params_1[int(0)].shoulder_0 + _S1651.differential_0;
    (&(&_S1363)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1782 = _S1770.crf_params_1[int(0)].toe_0 + _S1655.differential_0;
    *&((&(&(&_S1363)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1363)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1363)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1783 = _S1729 + _S1770.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1363)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1784 = _S1728 + _S1770.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1363)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1785 = _S1727 + _S1770.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1363)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1786 = _S1734 + _S1770.vignette_params_1[int(2)].cy_0;
    (&(&_S1363)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1787 = _S1735 + _S1770.vignette_params_1[int(2)].cx_0;
    (&(&_S1363)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1788 = _S1743 + _S1770.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1363)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1789 = _S1742 + _S1770.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1363)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1790 = _S1741 + _S1770.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1363)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1791 = _S1748 + _S1770.vignette_params_1[int(1)].cy_0;
    (&(&_S1363)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1792 = _S1749 + _S1770.vignette_params_1[int(1)].cx_0;
    (&(&_S1363)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1793 = _S1757 + _S1770.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1363)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1794 = _S1756 + _S1770.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1363)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1795 = _S1755 + _S1770.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1363)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1796 = _S1762 + _S1770.vignette_params_1[int(0)].cy_0;
    (&(&_S1363)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1797 = _S1763 + _S1770.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S1798;
    _S1798[int(0)] = 0.0f;
    _S1798[int(1)] = 0.0f;
    _S1798[int(2)] = 0.0f;
    _S1798[int(3)] = 0.0f;
    _S1798[int(4)] = 0.0f;
    _S1798[int(5)] = 0.0f;
    _S1798[int(6)] = 0.0f;
    _S1798[int(7)] = 0.0f;
    _S1798[int(8)] = 0.0f;
    _S1798[int(9)] = 0.0f;
    _S1798[int(10)] = 0.0f;
    _S1798[int(11)] = 0.0f;
    _S1798[int(12)] = 0.0f;
    _S1798[int(13)] = 0.0f;
    _S1798[int(14)] = 0.0f;
    _S1798[int(15)] = 0.0f;
    _S1798[int(16)] = 0.0f;
    _S1798[int(17)] = 0.0f;
    _S1798[int(18)] = 0.0f;
    _S1798[int(19)] = 0.0f;
    _S1798[int(20)] = 0.0f;
    _S1798[int(21)] = 0.0f;
    _S1798[int(22)] = 0.0f;
    _S1798[int(23)] = 0.0f;
    _S1798[int(24)] = 0.0f;
    _S1798[int(25)] = 0.0f;
    _S1798[int(26)] = 0.0f;
    _S1798[int(27)] = 0.0f;
    _S1798[int(28)] = 0.0f;
    _S1798[int(29)] = 0.0f;
    _S1798[int(30)] = 0.0f;
    _S1798[int(31)] = 0.0f;
    _S1798[int(32)] = 0.0f;
    _S1798[int(33)] = 0.0f;
    _S1798[int(34)] = 0.0f;
    _S1798[int(35)] = 0.0f;
    _S1798[int(8)] = _S1790;
    _S1798[int(16)] = _S1770.color_params_1.b_0.x;
    _S1798[int(15)] = _S1783;
    _S1798[int(14)] = _S1784;
    _S1798[int(13)] = _S1785;
    _S1798[int(12)] = _S1786;
    _S1798[int(11)] = _S1787;
    _S1798[int(10)] = _S1788;
    _S1798[int(9)] = _S1789;
    _S1798[int(17)] = _S1770.color_params_1.b_0.y;
    _S1798[int(7)] = _S1791;
    _S1798[int(6)] = _S1792;
    _S1798[int(5)] = _S1793;
    _S1798[int(4)] = _S1794;
    _S1798[int(3)] = _S1795;
    _S1798[int(2)] = _S1796;
    _S1798[int(1)] = _S1797;
    _S1798[int(0)] = _S1363.exposure_1;
    _S1798[int(26)] = _S1780;
    _S1798[int(34)] = _S1772;
    _S1798[int(33)] = _S1773;
    _S1798[int(32)] = _S1774;
    _S1798[int(31)] = _S1775;
    _S1798[int(30)] = _S1776;
    _S1798[int(29)] = _S1777;
    _S1798[int(28)] = _S1778;
    _S1798[int(27)] = _S1779;
    _S1798[int(35)] = _S1771;
    _S1798[int(25)] = _S1781;
    _S1798[int(24)] = _S1782;
    _S1798[int(23)] = _S1770.color_params_1.n_0.y;
    _S1798[int(22)] = _S1770.color_params_1.n_0.x;
    _S1798[int(21)] = _S1770.color_params_1.g_0.y;
    _S1798[int(20)] = _S1770.color_params_1.g_0.x;
    _S1798[int(19)] = _S1770.color_params_1.r_0.y;
    _S1798[int(18)] = _S1770.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1798;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1767;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1799, float2  _S1800, float2  _S1801, float2  _S1802, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1803, float3  _S1804)
{
    s_bwd_prop_apply_ppisp_0(_S1799, _S1800, _S1801, _S1802, _S1803, _S1804);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1805 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S1805;
    FixedArray<float, 36>  _S1806 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_2;
    (&dp_params_0)->differential_0 = _S1806;
    s_bwd_apply_ppisp_0(&dp_rgb_in_0, pix_coord_3, image_center_3, img_size_3, &dp_params_0, grad_out_0);
    *grad_rgb_in_0 = dp_rgb_in_0.differential_0;
    *grad_params_0 = (&dp_params_0)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C39x3E_0
{
    FixedArray<float, 39>  primal_0;
    FixedArray<float, 39>  differential_0;
};

inline __device__ void s_bwd_prop_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_1, float2  pix_coord_4, float2  image_center_4, float2  img_size_4, DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_1, float3  _s_dOut_14)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1807 = *dprgb_in_1;
    float3  _S1808 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1809 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1810 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1811 = {
        _S1810, _S1810, _S1810
    };
    float2  _S1812 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1813 = { _S1812, _S1812, _S1812, _S1812 };
    RQSCRFPPISPChannelParams_0 _S1814 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1815 = {
        _S1814, _S1814, _S1814
    };
    PPISPParamsRQS_0 _S1816;
    (&_S1816)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1816)->vignette_params_0 = _S1811;
    (&_S1816)->color_params_0 = _S1813;
    (&_S1816)->crf_params_0 = _S1815;
    (&(&_S1816)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1816)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S1817 = dpparams_1->primal_0[int(3)];
    (&(&_S1816)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S1818 = dpparams_1->primal_0[int(4)];
    (&(&_S1816)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S1819 = dpparams_1->primal_0[int(5)];
    (&(&_S1816)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1816)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1816)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S1820 = dpparams_1->primal_0[int(8)];
    (&(&_S1816)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S1821 = dpparams_1->primal_0[int(9)];
    (&(&_S1816)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S1822 = dpparams_1->primal_0[int(10)];
    (&(&_S1816)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1816)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1816)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S1823 = dpparams_1->primal_0[int(13)];
    (&(&_S1816)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S1824 = dpparams_1->primal_0[int(14)];
    (&(&_S1816)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S1825 = dpparams_1->primal_0[int(15)];
    (&(&_S1816)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1816)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1816)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1816)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1816)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1816)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1816)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1816)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1816)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S1826 = dpparams_1->primal_0[int(24)];
    (&(&_S1816)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S1827 = dpparams_1->primal_0[int(25)];
    (&(&_S1816)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S1828 = dpparams_1->primal_0[int(26)];
    (&(&_S1816)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S1829 = dpparams_1->primal_0[int(27)];
    (&(&_S1816)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S1830 = dpparams_1->primal_0[int(28)];
    (&(&_S1816)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S1831 = dpparams_1->primal_0[int(29)];
    (&(&_S1816)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S1832 = dpparams_1->primal_0[int(30)];
    (&(&_S1816)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S1833 = dpparams_1->primal_0[int(31)];
    (&(&_S1816)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S1834 = dpparams_1->primal_0[int(32)];
    (&(&_S1816)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S1835 = dpparams_1->primal_0[int(33)];
    (&(&_S1816)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S1836 = dpparams_1->primal_0[int(34)];
    (&(&_S1816)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S1837 = dpparams_1->primal_0[int(35)];
    (&(&_S1816)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S1838 = dpparams_1->primal_0[int(36)];
    (&(&_S1816)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S1839 = dpparams_1->primal_0[int(37)];
    (&(&_S1816)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S1840 = dpparams_1->primal_0[int(38)];
    (&(&_S1816)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S1841 = _S1816;
    float _S1842 = s_primal_ctx_exp2_0(_S1816.exposure_0);
    float3  _S1843 = make_float3 (_S1842);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S1842);
    float _S1844 = (F32_max((img_size_4.x), (img_size_4.y)));
    float _S1845 = (pix_coord_4.x - image_center_4.x) / _S1844;
    float _S1846 = (pix_coord_4.y - image_center_4.y) / _S1844;
    float dx_11 = _S1845 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S1846 - dpparams_1->primal_0[int(2)];
    float r2_17 = dx_11 * dx_11 + dy_9 * dy_9;
    float r4_9 = r2_17 * r2_17;
    float r6_3 = r4_9 * r2_17;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_17 + 1.0f;
    float _S1847 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S1848 = rgb_out_5.x * _S1847;
    float3  _S1849 = rgb_out_5;
    *&((&_S1849)->x) = _S1848;
    float dx_12 = _S1845 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S1846 - dpparams_1->primal_0[int(7)];
    float r2_18 = dx_12 * dx_12 + dy_10 * dy_10;
    float r4_10 = r2_18 * r2_18;
    float r6_4 = r4_10 * r2_18;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_18 + 1.0f;
    float _S1850 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S1849)->y) = rgb_out_5.y * _S1850;
    float dx_13 = _S1845 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S1846 - dpparams_1->primal_0[int(12)];
    float r2_19 = dx_13 * dx_13 + dy_11 * dy_11;
    float r4_11 = r2_19 * r2_19;
    float r6_5 = r4_11 * r2_19;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_19 + 1.0f;
    float _S1851 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S1849)->z) = rgb_out_5.z * _S1851;
    PPISPParamsRQS_0 _S1852 = _S1816;
    float2  _S1853 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1816.color_params_0.b_0);
    float2  _S1854 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1816.color_params_0.r_0);
    float2  _S1855 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1816.color_params_0.g_0);
    float2  _S1856 = s_primal_ctx_mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1816.color_params_0.n_0);
    float _S1857 = 0.3333333432674408f + _S1856.x;
    float _S1858 = 0.3333333432674408f + _S1856.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S1853.x, 1.0f + _S1854.x, _S1855.x, _S1853.y, _S1854.y, 1.0f + _S1855.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1858, 1.0f, 0.0f, - _S1857, - _S1858, _S1857, 0.0f);
    Matrix<float, 3, 3>  _S1859 = s_primal_ctx_mul_2(skew_1, T_3);
    float3  r0_3 = make_float3 (_S1859.rows[int(0)].x, _S1859.rows[int(0)].y, _S1859.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S1859.rows[int(1)].x, _S1859.rows[int(1)].y, _S1859.rows[int(1)].z);
    float3  r2_20 = make_float3 (_S1859.rows[int(2)].x, _S1859.rows[int(2)].y, _S1859.rows[int(2)].z);
    float3  _S1860 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S1861 = (s_primal_ctx_dot_0(_S1860, _S1860)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S1862;
    bool _S1863;
    if(_S1861)
    {
        float3  _S1864 = s_primal_ctx_cross_0(r0_3, r2_20);
        bool _S1865 = (s_primal_ctx_dot_0(_S1864, _S1864)) < 9.99999968265522539e-21f;
        if(_S1865)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_20);
        }
        else
        {
            lambda_v_7 = _S1864;
        }
        _S1863 = _S1865;
        _S1862 = _S1864;
    }
    else
    {
        lambda_v_7 = _S1860;
        _S1863 = false;
        _S1862 = _S1808;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S1866 = s_primal_ctx_mul_2(T_3, D_1);
    Matrix<float, 3, 3>  _S1867 = s_primal_ctx_mul_2(_S1866, S_inv_1);
    bool _S1868 = (s_primal_ctx_abs_0(_S1867.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S1869;
    float _S1870;
    if(_S1868)
    {
        float inv_s_1 = 1.0f / _S1867.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1871 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S1872 = _S1867.rows[int(2)].z * _S1867.rows[int(2)].z;
        H_5 = _S1867 * makeMatrix<float, 3, 3> (inv_s_1);
        _S1869 = _S1871;
        _S1870 = _S1872;
    }
    else
    {
        H_5 = _S1867;
        _S1869 = _S1809;
        _S1870 = 0.0f;
    }
    float _S1873 = _S1849.x;
    float _S1874 = _S1849.y;
    float intensity_3 = _S1873 + _S1874 + _S1849.z;
    float3  rgi_in_1 = make_float3 (_S1873, _S1874, intensity_3);
    float3  _S1875 = s_primal_ctx_mul_0(H_5, rgi_in_1);
    float _S1876 = _S1875.z + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S1876;
    float3  _S1877 = make_float3 (norm_factor_1);
    float _S1878 = _S1876 * _S1876;
    float3  rgi_out_5 = _S1875 * make_float3 (norm_factor_1);
    float _S1879 = rgi_out_5.x;
    float _S1880 = rgi_out_5.y;
    float3  _S1881 = make_float3 (_S1879, _S1880, rgi_out_5.z - _S1879 - _S1880);
    float3  _S1882 = make_float3 (0.0f);
    float3  _S1883 = make_float3 (1.0f);
    float3  _S1884 = s_primal_ctx_clamp_1(_S1881, _S1882, _S1883);
    float _S1885 = _S1884.x;
    float _S1886 = s_primal_ctx_exp_0(_S1826);
    float _S1887 = s_primal_ctx_exp_0(_S1827);
    float _S1888 = - _S1828;
    float _S1889 = 1.0f + s_primal_ctx_exp_0(_S1888);
    float x0_4 = 1.0f / _S1889;
    float _S1890 = _S1889 * _S1889;
    float _S1891 = - _S1829;
    float _S1892 = 1.0f + s_primal_ctx_exp_0(_S1891);
    float y0_4 = 1.0f / _S1892;
    float _S1893 = _S1892 * _S1892;
    float _S1894 = s_primal_ctx_exp_0(_S1830);
    bool _S1895 = _S1885 < x0_4;
    float _S1896;
    float _S1897;
    float _S1898;
    float _S1899;
    float _S1900;
    float _S1901;
    float _S1902;
    float _S1903;
    float _S1904;
    float _S1905;
    float _S1906;
    float _S1907;
    float _S1908;
    float _S1909;
    float _S1910;
    float _S1911;
    float _S1912;
    float _S1913;
    float _S1914;
    float _S1915;
    float _S1916;
    float _S1917;
    float _S1918;
    float _S1919;
    float _S1920;
    float _S1921;
    float _S1922;
    if(_S1895)
    {
        float s0_3 = y0_4 / x0_4;
        float _S1923 = x0_4 * x0_4;
        float t0_3 = _S1885 / x0_4;
        float _S1924 = s0_3 * t0_3;
        float _S1925 = _S1886 * t0_3;
        float _S1926 = 1.0f - t0_3;
        float _S1927 = _S1924 * t0_3 + _S1925 * _S1926;
        float _S1928 = y0_4 * _S1927;
        float _S1929 = _S1886 + _S1894 - 2.0f * s0_3;
        float _S1930 = _S1929 * t0_3;
        float _S1931 = s0_3 + _S1930 * _S1926;
        _S1896 = _S1931 * _S1931;
        _S1897 = _S1928;
        _S1898 = _S1931;
        _S1899 = _S1930;
        _S1900 = _S1926;
        _S1901 = _S1929;
        _S1902 = t0_3;
        _S1903 = _S1927;
        _S1904 = _S1925;
        _S1905 = _S1924;
        _S1906 = s0_3;
        _S1907 = _S1923;
        _S1908 = 0.0f;
        _S1909 = 0.0f;
        _S1910 = 0.0f;
        _S1911 = 0.0f;
        _S1912 = 0.0f;
        _S1913 = 0.0f;
        _S1914 = 0.0f;
        _S1915 = 0.0f;
        _S1916 = 0.0f;
        _S1917 = 0.0f;
        _S1918 = 0.0f;
        _S1919 = 0.0f;
        _S1920 = 0.0f;
        _S1921 = 0.0f;
        _S1922 = 0.0f;
    }
    else
    {
        float _S1932 = 1.0f - y0_4;
        float _S1933 = 1.0f - x0_4;
        float s1_7 = _S1932 / _S1933;
        float _S1934 = _S1933 * _S1933;
        float _S1935 = _S1885 - x0_4;
        float t1_3 = _S1935 / _S1933;
        float _S1936 = s1_7 * t1_3;
        float _S1937 = _S1894 * t1_3;
        float _S1938 = 1.0f - t1_3;
        float _S1939 = _S1936 * t1_3 + _S1937 * _S1938;
        float _S1940 = _S1932 * _S1939;
        float _S1941 = _S1894 + _S1887 - 2.0f * s1_7;
        float _S1942 = _S1941 * t1_3;
        float _S1943 = s1_7 + _S1942 * _S1938;
        float _S1944 = _S1943 * _S1943;
        _S1896 = 0.0f;
        _S1897 = 0.0f;
        _S1898 = 0.0f;
        _S1899 = 0.0f;
        _S1900 = 0.0f;
        _S1901 = 0.0f;
        _S1902 = 0.0f;
        _S1903 = 0.0f;
        _S1904 = 0.0f;
        _S1905 = 0.0f;
        _S1906 = 0.0f;
        _S1907 = 0.0f;
        _S1908 = _S1944;
        _S1909 = _S1940;
        _S1910 = _S1943;
        _S1911 = _S1942;
        _S1912 = _S1938;
        _S1913 = _S1941;
        _S1914 = t1_3;
        _S1915 = _S1932;
        _S1916 = _S1939;
        _S1917 = _S1937;
        _S1918 = _S1936;
        _S1919 = s1_7;
        _S1920 = _S1934;
        _S1921 = _S1935;
        _S1922 = _S1933;
    }
    float _S1945 = _S1884.y;
    float _S1946 = s_primal_ctx_exp_0(_S1831);
    float _S1947 = s_primal_ctx_exp_0(_S1832);
    float _S1948 = - _S1833;
    float _S1949 = 1.0f + s_primal_ctx_exp_0(_S1948);
    float x0_5 = 1.0f / _S1949;
    float _S1950 = _S1949 * _S1949;
    float _S1951 = - _S1834;
    float _S1952 = 1.0f + s_primal_ctx_exp_0(_S1951);
    float y0_5 = 1.0f / _S1952;
    float _S1953 = _S1952 * _S1952;
    float _S1954 = s_primal_ctx_exp_0(_S1835);
    bool _S1955 = _S1945 < x0_5;
    float _S1956;
    float _S1957;
    float _S1958;
    float _S1959;
    float _S1960;
    float _S1961;
    float _S1962;
    float _S1963;
    float _S1964;
    float _S1965;
    float _S1966;
    float _S1967;
    float _S1968;
    float _S1969;
    float _S1970;
    float _S1971;
    float _S1972;
    float _S1973;
    float _S1974;
    float _S1975;
    float _S1976;
    float _S1977;
    float _S1978;
    float _S1979;
    float _S1980;
    float _S1981;
    float _S1982;
    if(_S1955)
    {
        float s0_4 = y0_5 / x0_5;
        float _S1983 = x0_5 * x0_5;
        float t0_4 = _S1945 / x0_5;
        float _S1984 = s0_4 * t0_4;
        float _S1985 = _S1946 * t0_4;
        float _S1986 = 1.0f - t0_4;
        float _S1987 = _S1984 * t0_4 + _S1985 * _S1986;
        float _S1988 = y0_5 * _S1987;
        float _S1989 = _S1946 + _S1954 - 2.0f * s0_4;
        float _S1990 = _S1989 * t0_4;
        float _S1991 = s0_4 + _S1990 * _S1986;
        _S1956 = _S1991 * _S1991;
        _S1957 = _S1988;
        _S1958 = _S1991;
        _S1959 = _S1990;
        _S1960 = _S1986;
        _S1961 = _S1989;
        _S1962 = t0_4;
        _S1963 = _S1987;
        _S1964 = _S1985;
        _S1965 = _S1984;
        _S1966 = s0_4;
        _S1967 = _S1983;
        _S1968 = 0.0f;
        _S1969 = 0.0f;
        _S1970 = 0.0f;
        _S1971 = 0.0f;
        _S1972 = 0.0f;
        _S1973 = 0.0f;
        _S1974 = 0.0f;
        _S1975 = 0.0f;
        _S1976 = 0.0f;
        _S1977 = 0.0f;
        _S1978 = 0.0f;
        _S1979 = 0.0f;
        _S1980 = 0.0f;
        _S1981 = 0.0f;
        _S1982 = 0.0f;
    }
    else
    {
        float _S1992 = 1.0f - y0_5;
        float _S1993 = 1.0f - x0_5;
        float s1_8 = _S1992 / _S1993;
        float _S1994 = _S1993 * _S1993;
        float _S1995 = _S1945 - x0_5;
        float t1_4 = _S1995 / _S1993;
        float _S1996 = s1_8 * t1_4;
        float _S1997 = _S1954 * t1_4;
        float _S1998 = 1.0f - t1_4;
        float _S1999 = _S1996 * t1_4 + _S1997 * _S1998;
        float _S2000 = _S1992 * _S1999;
        float _S2001 = _S1954 + _S1947 - 2.0f * s1_8;
        float _S2002 = _S2001 * t1_4;
        float _S2003 = s1_8 + _S2002 * _S1998;
        float _S2004 = _S2003 * _S2003;
        _S1956 = 0.0f;
        _S1957 = 0.0f;
        _S1958 = 0.0f;
        _S1959 = 0.0f;
        _S1960 = 0.0f;
        _S1961 = 0.0f;
        _S1962 = 0.0f;
        _S1963 = 0.0f;
        _S1964 = 0.0f;
        _S1965 = 0.0f;
        _S1966 = 0.0f;
        _S1967 = 0.0f;
        _S1968 = _S2004;
        _S1969 = _S2000;
        _S1970 = _S2003;
        _S1971 = _S2002;
        _S1972 = _S1998;
        _S1973 = _S2001;
        _S1974 = t1_4;
        _S1975 = _S1992;
        _S1976 = _S1999;
        _S1977 = _S1997;
        _S1978 = _S1996;
        _S1979 = s1_8;
        _S1980 = _S1994;
        _S1981 = _S1995;
        _S1982 = _S1993;
    }
    float _S2005 = _S1884.z;
    float _S2006 = s_primal_ctx_exp_0(_S1836);
    float _S2007 = s_primal_ctx_exp_0(_S1837);
    float _S2008 = - _S1838;
    float _S2009 = 1.0f + s_primal_ctx_exp_0(_S2008);
    float x0_6 = 1.0f / _S2009;
    float _S2010 = _S2009 * _S2009;
    float _S2011 = - _S1839;
    float _S2012 = 1.0f + s_primal_ctx_exp_0(_S2011);
    float y0_6 = 1.0f / _S2012;
    float _S2013 = _S2012 * _S2012;
    float _S2014 = s_primal_ctx_exp_0(_S1840);
    bool _S2015 = _S2005 < x0_6;
    float _S2016;
    float _S2017;
    float _S2018;
    float _S2019;
    float _S2020;
    float _S2021;
    float _S2022;
    float _S2023;
    float _S2024;
    float _S2025;
    float _S2026;
    float _S2027;
    float _S2028;
    float _S2029;
    float _S2030;
    float _S2031;
    float _S2032;
    float _S2033;
    float _S2034;
    float _S2035;
    float _S2036;
    float _S2037;
    float _S2038;
    float _S2039;
    float _S2040;
    float _S2041;
    float _S2042;
    if(_S2015)
    {
        float s0_5 = y0_6 / x0_6;
        float _S2043 = x0_6 * x0_6;
        float t0_5 = _S2005 / x0_6;
        float _S2044 = s0_5 * t0_5;
        float _S2045 = _S2006 * t0_5;
        float _S2046 = 1.0f - t0_5;
        float _S2047 = _S2044 * t0_5 + _S2045 * _S2046;
        float _S2048 = y0_6 * _S2047;
        float _S2049 = _S2006 + _S2014 - 2.0f * s0_5;
        float _S2050 = _S2049 * t0_5;
        float _S2051 = s0_5 + _S2050 * _S2046;
        _S2016 = _S2051 * _S2051;
        _S2017 = _S2048;
        _S2018 = _S2051;
        _S2019 = _S2050;
        _S2020 = _S2046;
        _S2021 = _S2049;
        _S2022 = t0_5;
        _S2023 = _S2047;
        _S2024 = _S2045;
        _S2025 = _S2044;
        _S2026 = s0_5;
        _S2027 = _S2043;
        _S2028 = 0.0f;
        _S2029 = 0.0f;
        _S2030 = 0.0f;
        _S2031 = 0.0f;
        _S2032 = 0.0f;
        _S2033 = 0.0f;
        _S2034 = 0.0f;
        _S2035 = 0.0f;
        _S2036 = 0.0f;
        _S2037 = 0.0f;
        _S2038 = 0.0f;
        _S2039 = 0.0f;
        _S2040 = 0.0f;
        _S2041 = 0.0f;
        _S2042 = 0.0f;
    }
    else
    {
        float _S2052 = 1.0f - y0_6;
        float _S2053 = 1.0f - x0_6;
        float s1_9 = _S2052 / _S2053;
        float _S2054 = _S2053 * _S2053;
        float _S2055 = _S2005 - x0_6;
        float t1_5 = _S2055 / _S2053;
        float _S2056 = s1_9 * t1_5;
        float _S2057 = _S2014 * t1_5;
        float _S2058 = 1.0f - t1_5;
        float _S2059 = _S2056 * t1_5 + _S2057 * _S2058;
        float _S2060 = _S2052 * _S2059;
        float _S2061 = _S2014 + _S2007 - 2.0f * s1_9;
        float _S2062 = _S2061 * t1_5;
        float _S2063 = s1_9 + _S2062 * _S2058;
        float _S2064 = _S2063 * _S2063;
        _S2016 = 0.0f;
        _S2017 = 0.0f;
        _S2018 = 0.0f;
        _S2019 = 0.0f;
        _S2020 = 0.0f;
        _S2021 = 0.0f;
        _S2022 = 0.0f;
        _S2023 = 0.0f;
        _S2024 = 0.0f;
        _S2025 = 0.0f;
        _S2026 = 0.0f;
        _S2027 = 0.0f;
        _S2028 = _S2064;
        _S2029 = _S2060;
        _S2030 = _S2063;
        _S2031 = _S2062;
        _S2032 = _S2058;
        _S2033 = _S2061;
        _S2034 = t1_5;
        _S2035 = _S2052;
        _S2036 = _S2059;
        _S2037 = _S2057;
        _S2038 = _S2056;
        _S2039 = s1_9;
        _S2040 = _S2054;
        _S2041 = _S2055;
        _S2042 = _S2053;
    }
    if(_S2015)
    {
        float _S2065 = _s_dOut_14.z / _S2016;
        float _S2066 = _S2017 * - _S2065;
        float _S2067 = _S2018 * _S2065;
        float _S2068 = _S2020 * _S2066;
        float _S2069 = _S2022 * _S2068;
        float _S2070 = y0_6 * _S2067;
        float _S2071 = _S2020 * _S2070;
        float _S2072 = _S2022 * _S2070;
        float _S2073 = (_S2021 * _S2068 + - (_S2019 * _S2066 + _S2024 * _S2070) + _S2006 * _S2071 + _S2025 * _S2070 + _S2026 * _S2072) / _S2027;
        float _S2074 = x0_6 * _S2073;
        float _S2075 = (_S2066 + 2.0f * - _S2069 + _S2022 * _S2072) / _S2027;
        float _S2076 = _S2023 * _S2067 + x0_6 * _S2075;
        float _S2077 = _S2069 + _S2022 * _S2071;
        float _S2078 = _S2005 * - _S2073 + y0_6 * - _S2075;
        _S2016 = _S2069;
        _S2017 = _S2076;
        _S2018 = _S2078;
        _S2019 = 0.0f;
        _S2020 = _S2077;
        _S2021 = _S2074;
    }
    else
    {
        float _S2079 = _s_dOut_14.z / _S2028;
        float _S2080 = _S2029 * - _S2079;
        float _S2081 = _S2030 * _S2079;
        float _S2082 = _S2032 * _S2080;
        float _S2083 = _S2034 * _S2082;
        float _S2084 = _S2035 * _S2081;
        float _S2085 = _S2032 * _S2084;
        float _S2086 = _S2034 * _S2084;
        float _S2087 = (_S2033 * _S2082 + - (_S2031 * _S2080 + _S2037 * _S2084) + _S2014 * _S2085 + _S2038 * _S2084 + _S2039 * _S2086) / _S2040;
        float _S2088 = _S2042 * _S2087;
        float _S2089 = (_S2080 + 2.0f * - _S2083 + _S2034 * _S2086) / _S2040;
        float _S2090 = _s_dOut_14.z + - (_S2036 * _S2081 + _S2042 * _S2089);
        float _S2091 = - _S2088 + - (_S2041 * - _S2087 + _S2035 * - _S2089);
        _S2016 = _S2083 + _S2034 * _S2085;
        _S2017 = _S2090;
        _S2018 = _S2091;
        _S2019 = _S2083;
        _S2020 = 0.0f;
        _S2021 = _S2088;
    }
    DiffPair_float_0 _S2092;
    (&_S2092)->primal_0 = _S1840;
    (&_S2092)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2092, _S2016);
    DiffPair_float_0 _S2093 = _S2092;
    float _S2094 = - (_S2017 / _S2013);
    DiffPair_float_0 _S2095;
    (&_S2095)->primal_0 = _S2011;
    (&_S2095)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2095, _S2094);
    float _S2096 = - _S2095.differential_0;
    float _S2097 = - (_S2018 / _S2010);
    DiffPair_float_0 _S2098;
    (&_S2098)->primal_0 = _S2008;
    (&_S2098)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2098, _S2097);
    float _S2099 = - _S2098.differential_0;
    DiffPair_float_0 _S2100;
    (&_S2100)->primal_0 = _S1837;
    (&_S2100)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2100, _S2019);
    DiffPair_float_0 _S2101 = _S2100;
    DiffPair_float_0 _S2102;
    (&_S2102)->primal_0 = _S1836;
    (&_S2102)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2102, _S2020);
    DiffPair_float_0 _S2103 = _S2102;
    float3  _S2104 = make_float3 (0.0f, 0.0f, _S2021);
    if(_S1955)
    {
        float _S2105 = _s_dOut_14.y / _S1956;
        float _S2106 = _S1957 * - _S2105;
        float _S2107 = _S1958 * _S2105;
        float _S2108 = _S1960 * _S2106;
        float _S2109 = _S1962 * _S2108;
        float _S2110 = y0_5 * _S2107;
        float _S2111 = _S1960 * _S2110;
        float _S2112 = _S1962 * _S2110;
        float _S2113 = (_S1961 * _S2108 + - (_S1959 * _S2106 + _S1964 * _S2110) + _S1946 * _S2111 + _S1965 * _S2110 + _S1966 * _S2112) / _S1967;
        float _S2114 = x0_5 * _S2113;
        float _S2115 = (_S2106 + 2.0f * - _S2109 + _S1962 * _S2112) / _S1967;
        float _S2116 = _S1963 * _S2107 + x0_5 * _S2115;
        float _S2117 = _S2109 + _S1962 * _S2111;
        float _S2118 = _S1945 * - _S2113 + y0_5 * - _S2115;
        _S1956 = _S2109;
        _S1957 = _S2116;
        _S1958 = _S2118;
        _S1959 = 0.0f;
        _S1960 = _S2117;
        _S1961 = _S2114;
    }
    else
    {
        float _S2119 = _s_dOut_14.y / _S1968;
        float _S2120 = _S1969 * - _S2119;
        float _S2121 = _S1970 * _S2119;
        float _S2122 = _S1972 * _S2120;
        float _S2123 = _S1974 * _S2122;
        float _S2124 = _S1975 * _S2121;
        float _S2125 = _S1972 * _S2124;
        float _S2126 = _S1974 * _S2124;
        float _S2127 = (_S1973 * _S2122 + - (_S1971 * _S2120 + _S1977 * _S2124) + _S1954 * _S2125 + _S1978 * _S2124 + _S1979 * _S2126) / _S1980;
        float _S2128 = _S1982 * _S2127;
        float _S2129 = (_S2120 + 2.0f * - _S2123 + _S1974 * _S2126) / _S1980;
        float _S2130 = _s_dOut_14.y + - (_S1976 * _S2121 + _S1982 * _S2129);
        float _S2131 = - _S2128 + - (_S1981 * - _S2127 + _S1975 * - _S2129);
        _S1956 = _S2123 + _S1974 * _S2125;
        _S1957 = _S2130;
        _S1958 = _S2131;
        _S1959 = _S2123;
        _S1960 = 0.0f;
        _S1961 = _S2128;
    }
    DiffPair_float_0 _S2132;
    (&_S2132)->primal_0 = _S1835;
    (&_S2132)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2132, _S1956);
    DiffPair_float_0 _S2133 = _S2132;
    float _S2134 = - (_S1957 / _S1953);
    DiffPair_float_0 _S2135;
    (&_S2135)->primal_0 = _S1951;
    (&_S2135)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2135, _S2134);
    float _S2136 = - _S2135.differential_0;
    float _S2137 = - (_S1958 / _S1950);
    DiffPair_float_0 _S2138;
    (&_S2138)->primal_0 = _S1948;
    (&_S2138)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2138, _S2137);
    float _S2139 = - _S2138.differential_0;
    DiffPair_float_0 _S2140;
    (&_S2140)->primal_0 = _S1832;
    (&_S2140)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2140, _S1959);
    DiffPair_float_0 _S2141 = _S2140;
    DiffPair_float_0 _S2142;
    (&_S2142)->primal_0 = _S1831;
    (&_S2142)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2142, _S1960);
    DiffPair_float_0 _S2143 = _S2142;
    float3  _S2144 = _S2104 + make_float3 (0.0f, _S1961, 0.0f);
    if(_S1895)
    {
        float _S2145 = _s_dOut_14.x / _S1896;
        float _S2146 = _S1897 * - _S2145;
        float _S2147 = _S1898 * _S2145;
        float _S2148 = _S1900 * _S2146;
        float _S2149 = _S1902 * _S2148;
        float _S2150 = y0_4 * _S2147;
        float _S2151 = _S1900 * _S2150;
        float _S2152 = _S1902 * _S2150;
        float _S2153 = (_S1901 * _S2148 + - (_S1899 * _S2146 + _S1904 * _S2150) + _S1886 * _S2151 + _S1905 * _S2150 + _S1906 * _S2152) / _S1907;
        float _S2154 = x0_4 * _S2153;
        float _S2155 = (_S2146 + 2.0f * - _S2149 + _S1902 * _S2152) / _S1907;
        float _S2156 = _S1903 * _S2147 + x0_4 * _S2155;
        float _S2157 = _S2149 + _S1902 * _S2151;
        float _S2158 = _S1885 * - _S2153 + y0_4 * - _S2155;
        _S1896 = _S2149;
        _S1897 = _S2156;
        _S1898 = _S2158;
        _S1899 = 0.0f;
        _S1900 = _S2157;
        _S1901 = _S2154;
    }
    else
    {
        float _S2159 = _s_dOut_14.x / _S1908;
        float _S2160 = _S1909 * - _S2159;
        float _S2161 = _S1910 * _S2159;
        float _S2162 = _S1912 * _S2160;
        float _S2163 = _S1914 * _S2162;
        float _S2164 = _S1915 * _S2161;
        float _S2165 = _S1912 * _S2164;
        float _S2166 = _S1914 * _S2164;
        float _S2167 = (_S1913 * _S2162 + - (_S1911 * _S2160 + _S1917 * _S2164) + _S1894 * _S2165 + _S1918 * _S2164 + _S1919 * _S2166) / _S1920;
        float _S2168 = _S1922 * _S2167;
        float _S2169 = (_S2160 + 2.0f * - _S2163 + _S1914 * _S2166) / _S1920;
        float _S2170 = _s_dOut_14.x + - (_S1916 * _S2161 + _S1922 * _S2169);
        float _S2171 = - _S2168 + - (_S1921 * - _S2167 + _S1915 * - _S2169);
        _S1896 = _S2163 + _S1914 * _S2165;
        _S1897 = _S2170;
        _S1898 = _S2171;
        _S1899 = _S2163;
        _S1900 = 0.0f;
        _S1901 = _S2168;
    }
    DiffPair_float_0 _S2172;
    (&_S2172)->primal_0 = _S1830;
    (&_S2172)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2172, _S1896);
    DiffPair_float_0 _S2173 = _S2172;
    float _S2174 = - (_S1897 / _S1893);
    DiffPair_float_0 _S2175;
    (&_S2175)->primal_0 = _S1891;
    (&_S2175)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2175, _S2174);
    float _S2176 = - _S2175.differential_0;
    float _S2177 = - (_S1898 / _S1890);
    DiffPair_float_0 _S2178;
    (&_S2178)->primal_0 = _S1888;
    (&_S2178)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2178, _S2177);
    float _S2179 = - _S2178.differential_0;
    DiffPair_float_0 _S2180;
    (&_S2180)->primal_0 = _S1827;
    (&_S2180)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2180, _S1899);
    DiffPair_float_0 _S2181 = _S2180;
    DiffPair_float_0 _S2182;
    (&_S2182)->primal_0 = _S1826;
    (&_S2182)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2182, _S1900);
    DiffPair_float_0 _S2183 = _S2182;
    float3  _S2184 = _S2144 + make_float3 (_S1901, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2185;
    (&_S2185)->primal_0 = _S1881;
    (&_S2185)->differential_0 = _S1808;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2186;
    (&_S2186)->primal_0 = _S1882;
    (&_S2186)->differential_0 = _S1808;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2187;
    (&_S2187)->primal_0 = _S1883;
    (&_S2187)->differential_0 = _S1808;
    s_bwd_prop_clamp_1(&_S2185, &_S2186, &_S2187, _S2184);
    float _S2188 = - _S2185.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S2185.differential_0.x + _S2188, _S2185.differential_0.y + _S2188, _S2185.differential_0.z);
    float3  _S2189 = _S1875 * s_diff_rgi_out_T_1;
    float _S2190 = (_S2189.x + _S2189.y + _S2189.z) / _S1878;
    float _S2191 = _S1876 * _S2190;
    float3  _S2192 = _S1877 * s_diff_rgi_out_T_1 + make_float3 (0.0f, 0.0f, intensity_3 * - _S2190);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2193;
    (&_S2193)->primal_0 = H_5;
    (&_S2193)->differential_0 = _S1809;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2194;
    (&_S2194)->primal_0 = rgi_in_1;
    (&_S2194)->differential_0 = _S1808;
    s_bwd_prop_mul_2(&_S2193, &_S2194, _S2192);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2195 = _S2193;
    float _S2196 = _S2191 + _S2194.differential_0.z;
    float _S2197 = _S2194.differential_0.y + _S2196;
    float _S2198 = _S2194.differential_0.x + _S2196;
    float3  _S2199 = make_float3 (_S2198, _S2197, _S2196);
    if(_S1868)
    {
        Matrix<float, 3, 3>  _S2200 = _S1867 * _S2195.differential_0;
        Matrix<float, 3, 3>  _S2201 = _S1869 * _S2195.differential_0;
        _S1870 = - ((_S2200.rows[int(0)].x + _S2200.rows[int(0)].y + _S2200.rows[int(0)].z + _S2200.rows[int(1)].x + _S2200.rows[int(1)].y + _S2200.rows[int(1)].z + _S2200.rows[int(2)].x + _S2200.rows[int(2)].y + _S2200.rows[int(2)].z) / _S1870);
        H_5 = _S2201;
    }
    else
    {
        _S1870 = 0.0f;
        H_5 = _S2195.differential_0;
    }
    DiffPair_float_0 _S2202;
    (&_S2202)->primal_0 = _S1867.rows[int(2)].z;
    (&_S2202)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2202, 0.0f);
    float _S2203 = _S2202.differential_0 + _S1870;
    float3  _S2204 = _S1808;
    *&((&_S2204)->z) = _S2203;
    Matrix<float, 3, 3>  _S2205 = _S1809;
    _S2205[int(2)] = _S2204;
    Matrix<float, 3, 3>  _S2206 = H_5 + _S2205;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2207;
    (&_S2207)->primal_0 = _S1866;
    (&_S2207)->differential_0 = _S1809;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2208;
    (&_S2208)->primal_0 = S_inv_1;
    (&_S2208)->differential_0 = _S1809;
    s_bwd_prop_mul_1(&_S2207, &_S2208, _S2206);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2209;
    (&_S2209)->primal_0 = T_3;
    (&_S2209)->differential_0 = _S1809;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2210;
    (&_S2210)->primal_0 = D_1;
    (&_S2210)->differential_0 = _S1809;
    s_bwd_prop_mul_1(&_S2209, &_S2210, _S2207.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2211 = _S2209;
    float3  _S2212 = make_float3 (_S2210.differential_0.rows[int(0)].x, _S2210.differential_0.rows[int(1)].y, _S2210.differential_0.rows[int(2)].z);
    float3  _S2213;
    if(_S1861)
    {
        if(_S1863)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2214;
            (&_S2214)->primal_0 = r1_3;
            (&_S2214)->differential_0 = _S1808;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2215;
            (&_S2215)->primal_0 = r2_20;
            (&_S2215)->differential_0 = _S1808;
            s_bwd_prop_cross_0(&_S2214, &_S2215, _S2212);
            _S1849 = _S1808;
            lambda_v_7 = _S2215.differential_0;
            _S2213 = _S2214.differential_0;
        }
        else
        {
            _S1849 = _S2212;
            lambda_v_7 = _S1808;
            _S2213 = _S1808;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2216;
        (&_S2216)->primal_0 = _S1862;
        (&_S2216)->differential_0 = _S1808;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2217;
        (&_S2217)->primal_0 = _S1862;
        (&_S2217)->differential_0 = _S1808;
        s_bwd_prop_dot_0(&_S2216, &_S2217, 0.0f);
        float3  _S2218 = _S2217.differential_0 + _S2216.differential_0 + _S1849;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2219;
        (&_S2219)->primal_0 = r0_3;
        (&_S2219)->differential_0 = _S1808;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2220;
        (&_S2220)->primal_0 = r2_20;
        (&_S2220)->differential_0 = _S1808;
        s_bwd_prop_cross_0(&_S2219, &_S2220, _S2218);
        float3  _S2221 = _S2220.differential_0 + lambda_v_7;
        _S1849 = _S1808;
        lambda_v_7 = _S2221;
        _S1862 = _S2213;
        _S2213 = _S2219.differential_0;
    }
    else
    {
        _S1849 = _S2212;
        lambda_v_7 = _S1808;
        _S1862 = _S1808;
        _S2213 = _S1808;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2222;
    (&_S2222)->primal_0 = _S1860;
    (&_S2222)->differential_0 = _S1808;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2223;
    (&_S2223)->primal_0 = _S1860;
    (&_S2223)->differential_0 = _S1808;
    s_bwd_prop_dot_0(&_S2222, &_S2223, 0.0f);
    float3  _S2224 = _S2223.differential_0 + _S2222.differential_0 + _S1849;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2225;
    (&_S2225)->primal_0 = r0_3;
    (&_S2225)->differential_0 = _S1808;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2226;
    (&_S2226)->primal_0 = r1_3;
    (&_S2226)->differential_0 = _S1808;
    s_bwd_prop_cross_0(&_S2225, &_S2226, _S2224);
    float3  _S2227 = _S1808;
    *&((&_S2227)->z) = lambda_v_7.z;
    *&((&_S2227)->y) = lambda_v_7.y;
    *&((&_S2227)->x) = lambda_v_7.x;
    float3  _S2228 = _S2226.differential_0 + _S1862;
    float3  _S2229 = _S1808;
    *&((&_S2229)->z) = _S2228.z;
    *&((&_S2229)->y) = _S2228.y;
    *&((&_S2229)->x) = _S2228.x;
    float3  _S2230 = _S2225.differential_0 + _S2213;
    float3  _S2231 = _S1808;
    *&((&_S2231)->z) = _S2230.z;
    *&((&_S2231)->y) = _S2230.y;
    *&((&_S2231)->x) = _S2230.x;
    Matrix<float, 3, 3>  _S2232 = _S1809;
    _S2232[int(2)] = _S2227;
    _S2232[int(1)] = _S2229;
    _S2232[int(0)] = _S2231;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2233;
    (&_S2233)->primal_0 = skew_1;
    (&_S2233)->differential_0 = _S1809;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2234;
    (&_S2234)->primal_0 = T_3;
    (&_S2234)->differential_0 = _S1809;
    s_bwd_prop_mul_1(&_S2233, &_S2234, _S2232);
    Matrix<float, 3, 3>  _S2235 = _S2234.differential_0 + _S2211.differential_0;
    float2  _S2236 = make_float2 (_S2233.differential_0.rows[int(2)].y + - _S2233.differential_0.rows[int(1)].z, _S2233.differential_0.rows[int(0)].z + - _S2233.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S2237 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2238;
    (&_S2238)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2238)->differential_0 = _S2237;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2239;
    (&_S2239)->primal_0 = _S1852.color_params_0.n_0;
    (&_S2239)->differential_0 = _S1812;
    s_bwd_prop_mul_3(&_S2238, &_S2239, _S2236);
    float2  _S2240 = make_float2 (_S2235.rows[int(0)].z, _S2235.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2241;
    (&_S2241)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2241)->differential_0 = _S2237;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2242;
    (&_S2242)->primal_0 = _S1852.color_params_0.g_0;
    (&_S2242)->differential_0 = _S1812;
    s_bwd_prop_mul_3(&_S2241, &_S2242, _S2240);
    float2  _S2243 = make_float2 (_S2235.rows[int(0)].y, _S2235.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2244;
    (&_S2244)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2244)->differential_0 = _S2237;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2245;
    (&_S2245)->primal_0 = _S1852.color_params_0.r_0;
    (&_S2245)->differential_0 = _S1812;
    s_bwd_prop_mul_3(&_S2244, &_S2245, _S2243);
    float2  _S2246 = make_float2 (_S2235.rows[int(0)].x, _S2235.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2247;
    (&_S2247)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2247)->differential_0 = _S2237;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2248;
    (&_S2248)->primal_0 = _S1852.color_params_0.b_0;
    (&_S2248)->differential_0 = _S1812;
    s_bwd_prop_mul_3(&_S2247, &_S2248, _S2246);
    ColorPPISPParams_0 _S2249 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2249)->n_0 = _S2239.differential_0;
    (&_S2249)->g_0 = _S2242.differential_0;
    (&_S2249)->r_0 = _S2245.differential_0;
    (&_S2249)->b_0 = _S2248.differential_0;
    _S1849 = _S2199;
    *&((&_S1849)->z) = 0.0f;
    float _S2250 = rgb_out_5.z * _S2196;
    float _S2251 = _S1851 * _S2196;
    DiffPair_float_0 _S2252;
    (&_S2252)->primal_0 = falloff_5;
    (&_S2252)->differential_0 = 0.0f;
    DiffPair_float_0 _S2253;
    (&_S2253)->primal_0 = 0.0f;
    (&_S2253)->differential_0 = 0.0f;
    DiffPair_float_0 _S2254;
    (&_S2254)->primal_0 = 1.0f;
    (&_S2254)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2252, &_S2253, &_S2254, _S2250);
    float _S2255 = r2_19 * _S2252.differential_0;
    float _S2256 = r4_11 * _S2252.differential_0;
    float s_diff_r6_T_3 = _S1825 * _S2252.differential_0;
    float _S2257 = r6_5 * _S2252.differential_0;
    float _S2258 = r2_19 * (_S1824 * _S2252.differential_0 + r2_19 * s_diff_r6_T_3);
    float _S2259 = _S1823 * _S2252.differential_0 + r4_11 * s_diff_r6_T_3 + _S2258 + _S2258;
    float _S2260 = dy_11 * _S2259;
    float _S2261 = dx_13 * _S2259;
    float _S2262 = - (_S2260 + _S2260);
    float _S2263 = - (_S2261 + _S2261);
    *&((&_S1849)->y) = 0.0f;
    float _S2264 = rgb_out_5.y * _S2197;
    float _S2265 = _S1850 * _S2197;
    DiffPair_float_0 _S2266;
    (&_S2266)->primal_0 = falloff_4;
    (&_S2266)->differential_0 = 0.0f;
    DiffPair_float_0 _S2267;
    (&_S2267)->primal_0 = 0.0f;
    (&_S2267)->differential_0 = 0.0f;
    DiffPair_float_0 _S2268;
    (&_S2268)->primal_0 = 1.0f;
    (&_S2268)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2266, &_S2267, &_S2268, _S2264);
    float _S2269 = r2_18 * _S2266.differential_0;
    float _S2270 = r4_10 * _S2266.differential_0;
    float s_diff_r6_T_4 = _S1822 * _S2266.differential_0;
    float _S2271 = r6_4 * _S2266.differential_0;
    float _S2272 = r2_18 * (_S1821 * _S2266.differential_0 + r2_18 * s_diff_r6_T_4);
    float _S2273 = _S1820 * _S2266.differential_0 + r4_10 * s_diff_r6_T_4 + _S2272 + _S2272;
    float _S2274 = dy_10 * _S2273;
    float _S2275 = dx_12 * _S2273;
    float _S2276 = - (_S2274 + _S2274);
    float _S2277 = - (_S2275 + _S2275);
    *&((&_S1849)->x) = 0.0f;
    float _S2278 = rgb_out_5.x * _S2198;
    float _S2279 = _S1847 * _S2198;
    DiffPair_float_0 _S2280;
    (&_S2280)->primal_0 = falloff_3;
    (&_S2280)->differential_0 = 0.0f;
    DiffPair_float_0 _S2281;
    (&_S2281)->primal_0 = 0.0f;
    (&_S2281)->differential_0 = 0.0f;
    DiffPair_float_0 _S2282;
    (&_S2282)->primal_0 = 1.0f;
    (&_S2282)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2280, &_S2281, &_S2282, _S2278);
    float _S2283 = r2_17 * _S2280.differential_0;
    float _S2284 = r4_9 * _S2280.differential_0;
    float s_diff_r6_T_5 = _S1819 * _S2280.differential_0;
    float _S2285 = r6_3 * _S2280.differential_0;
    float _S2286 = r2_17 * (_S1818 * _S2280.differential_0 + r2_17 * s_diff_r6_T_5);
    float _S2287 = _S1817 * _S2280.differential_0 + r4_9 * s_diff_r6_T_5 + _S2286 + _S2286;
    float _S2288 = dy_9 * _S2287;
    float _S2289 = dx_11 * _S2287;
    float _S2290 = - (_S2288 + _S2288);
    float _S2291 = - (_S2289 + _S2289);
    float3  _S2292 = _S1808;
    *&((&_S2292)->z) = _S2251;
    *&((&_S2292)->y) = _S2265;
    *&((&_S2292)->x) = _S2279;
    float3  _S2293 = _S1849 + _S2292;
    float3  _S2294 = _S1807.primal_0 * _S2293;
    float3  _S2295 = _S1843 * _S2293;
    float _S2296 = _S2294.x + _S2294.y + _S2294.z;
    DiffPair_float_0 _S2297;
    (&_S2297)->primal_0 = _S1841.exposure_0;
    (&_S2297)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2297, _S2296);
    PPISPParamsRQS_0 _S2298 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2298)->color_params_0 = _S2249;
    (&_S2298)->exposure_0 = _S2297.differential_0;
    _S1816 = _S2298;
    (&(&_S1816)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2299 = _S2298.crf_params_0[int(2)].gc_0 + _S2093.differential_0;
    (&(&_S1816)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2300 = _S2298.crf_params_0[int(2)].y0_0 + _S2096;
    (&(&_S1816)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2301 = _S2298.crf_params_0[int(2)].x0_0 + _S2099;
    (&(&_S1816)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2302 = _S2298.crf_params_0[int(2)].g1_0 + _S2101.differential_0;
    (&(&_S1816)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2303 = _S2298.crf_params_0[int(2)].g0_0 + _S2103.differential_0;
    (&(&_S1816)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2304 = _S2298.crf_params_0[int(1)].gc_0 + _S2133.differential_0;
    (&(&_S1816)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2305 = _S2298.crf_params_0[int(1)].y0_0 + _S2136;
    (&(&_S1816)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2306 = _S2298.crf_params_0[int(1)].x0_0 + _S2139;
    (&(&_S1816)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2307 = _S2298.crf_params_0[int(1)].g1_0 + _S2141.differential_0;
    (&(&_S1816)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2308 = _S2298.crf_params_0[int(1)].g0_0 + _S2143.differential_0;
    (&(&_S1816)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2309 = _S2298.crf_params_0[int(0)].gc_0 + _S2173.differential_0;
    (&(&_S1816)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2310 = _S2298.crf_params_0[int(0)].y0_0 + _S2176;
    (&(&_S1816)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2311 = _S2298.crf_params_0[int(0)].x0_0 + _S2179;
    (&(&_S1816)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2312 = _S2298.crf_params_0[int(0)].g1_0 + _S2181.differential_0;
    (&(&_S1816)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2313 = _S2298.crf_params_0[int(0)].g0_0 + _S2183.differential_0;
    *&((&(&(&_S1816)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1816)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1816)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2314 = _S2257 + _S2298.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1816)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2315 = _S2256 + _S2298.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1816)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2316 = _S2255 + _S2298.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1816)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2317 = _S2262 + _S2298.vignette_params_0[int(2)].cy_0;
    (&(&_S1816)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2318 = _S2263 + _S2298.vignette_params_0[int(2)].cx_0;
    (&(&_S1816)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2319 = _S2271 + _S2298.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1816)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2320 = _S2270 + _S2298.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1816)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2321 = _S2269 + _S2298.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1816)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2322 = _S2276 + _S2298.vignette_params_0[int(1)].cy_0;
    (&(&_S1816)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2323 = _S2277 + _S2298.vignette_params_0[int(1)].cx_0;
    (&(&_S1816)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2324 = _S2285 + _S2298.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1816)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2325 = _S2284 + _S2298.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1816)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2326 = _S2283 + _S2298.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1816)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2327 = _S2290 + _S2298.vignette_params_0[int(0)].cy_0;
    (&(&_S1816)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2328 = _S2291 + _S2298.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2329;
    _S2329[int(0)] = 0.0f;
    _S2329[int(1)] = 0.0f;
    _S2329[int(2)] = 0.0f;
    _S2329[int(3)] = 0.0f;
    _S2329[int(4)] = 0.0f;
    _S2329[int(5)] = 0.0f;
    _S2329[int(6)] = 0.0f;
    _S2329[int(7)] = 0.0f;
    _S2329[int(8)] = 0.0f;
    _S2329[int(9)] = 0.0f;
    _S2329[int(10)] = 0.0f;
    _S2329[int(11)] = 0.0f;
    _S2329[int(12)] = 0.0f;
    _S2329[int(13)] = 0.0f;
    _S2329[int(14)] = 0.0f;
    _S2329[int(15)] = 0.0f;
    _S2329[int(16)] = 0.0f;
    _S2329[int(17)] = 0.0f;
    _S2329[int(18)] = 0.0f;
    _S2329[int(19)] = 0.0f;
    _S2329[int(20)] = 0.0f;
    _S2329[int(21)] = 0.0f;
    _S2329[int(22)] = 0.0f;
    _S2329[int(23)] = 0.0f;
    _S2329[int(24)] = 0.0f;
    _S2329[int(25)] = 0.0f;
    _S2329[int(26)] = 0.0f;
    _S2329[int(27)] = 0.0f;
    _S2329[int(28)] = 0.0f;
    _S2329[int(29)] = 0.0f;
    _S2329[int(30)] = 0.0f;
    _S2329[int(31)] = 0.0f;
    _S2329[int(32)] = 0.0f;
    _S2329[int(33)] = 0.0f;
    _S2329[int(34)] = 0.0f;
    _S2329[int(35)] = 0.0f;
    _S2329[int(36)] = 0.0f;
    _S2329[int(37)] = 0.0f;
    _S2329[int(38)] = 0.0f;
    _S2329[int(9)] = _S2320;
    _S2329[int(18)] = _S2298.color_params_0.r_0.x;
    _S2329[int(17)] = _S2298.color_params_0.b_0.y;
    _S2329[int(16)] = _S2298.color_params_0.b_0.x;
    _S2329[int(15)] = _S2314;
    _S2329[int(14)] = _S2315;
    _S2329[int(13)] = _S2316;
    _S2329[int(12)] = _S2317;
    _S2329[int(11)] = _S2318;
    _S2329[int(10)] = _S2319;
    _S2329[int(19)] = _S2298.color_params_0.r_0.y;
    _S2329[int(8)] = _S2321;
    _S2329[int(7)] = _S2322;
    _S2329[int(6)] = _S2323;
    _S2329[int(5)] = _S2324;
    _S2329[int(4)] = _S2325;
    _S2329[int(3)] = _S2326;
    _S2329[int(2)] = _S2327;
    _S2329[int(1)] = _S2328;
    _S2329[int(0)] = _S1816.exposure_0;
    _S2329[int(28)] = _S2309;
    _S2329[int(37)] = _S2300;
    _S2329[int(36)] = _S2301;
    _S2329[int(35)] = _S2302;
    _S2329[int(34)] = _S2303;
    _S2329[int(33)] = _S2304;
    _S2329[int(32)] = _S2305;
    _S2329[int(31)] = _S2306;
    _S2329[int(30)] = _S2307;
    _S2329[int(29)] = _S2308;
    _S2329[int(38)] = _S2299;
    _S2329[int(27)] = _S2310;
    _S2329[int(26)] = _S2311;
    _S2329[int(25)] = _S2312;
    _S2329[int(24)] = _S2313;
    _S2329[int(23)] = _S2298.color_params_0.n_0.y;
    _S2329[int(22)] = _S2298.color_params_0.n_0.x;
    _S2329[int(21)] = _S2298.color_params_0.g_0.y;
    _S2329[int(20)] = _S2298.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S2329;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S2295;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2330, float2  _S2331, float2  _S2332, float2  _S2333, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2334, float3  _S2335)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S2330, _S2331, _S2332, _S2333, _S2334, _S2335);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S2336 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S2336;
    FixedArray<float, 39>  _S2337 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_3;
    (&dp_params_1)->differential_0 = _S2337;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_4, FixedArray<float, 22>  * _S2338)
{
    PPISPParams_0 p_2;
    (&p_2)->exposure_1 = params_4[int(0)];
    (&(&p_2)->vignette_params_1[int(0)])->cx_0 = params_4[int(1)];
    (&(&p_2)->vignette_params_1[int(0)])->cy_0 = params_4[int(2)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha0_0 = params_4[int(3)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha1_0 = params_4[int(4)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha2_0 = params_4[int(5)];
    (&(&p_2)->vignette_params_1[int(1)])->cx_0 = params_4[int(6)];
    (&(&p_2)->vignette_params_1[int(1)])->cy_0 = params_4[int(7)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha0_0 = params_4[int(8)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha1_0 = params_4[int(9)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha2_0 = params_4[int(10)];
    (&(&p_2)->vignette_params_1[int(2)])->cx_0 = params_4[int(11)];
    (&(&p_2)->vignette_params_1[int(2)])->cy_0 = params_4[int(12)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha0_0 = params_4[int(13)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha1_0 = params_4[int(14)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha2_0 = params_4[int(15)];
    *&((&(&(&p_2)->color_params_1)->b_0)->x) = params_4[int(16)];
    *&((&(&(&p_2)->color_params_1)->b_0)->y) = params_4[int(17)];
    *&((&(&(&p_2)->color_params_1)->r_0)->x) = params_4[int(18)];
    *&((&(&(&p_2)->color_params_1)->r_0)->y) = params_4[int(19)];
    *&((&(&(&p_2)->color_params_1)->g_0)->x) = params_4[int(20)];
    *&((&(&(&p_2)->color_params_1)->g_0)->y) = params_4[int(21)];
    *&((&(&(&p_2)->color_params_1)->n_0)->x) = params_4[int(22)];
    *&((&(&(&p_2)->color_params_1)->n_0)->y) = params_4[int(23)];
    (&(&p_2)->crf_params_1[int(0)])->toe_0 = params_4[int(24)];
    (&(&p_2)->crf_params_1[int(0)])->shoulder_0 = params_4[int(25)];
    (&(&p_2)->crf_params_1[int(0)])->gamma_0 = params_4[int(26)];
    (&(&p_2)->crf_params_1[int(0)])->center_0 = params_4[int(27)];
    (&(&p_2)->crf_params_1[int(1)])->toe_0 = params_4[int(28)];
    (&(&p_2)->crf_params_1[int(1)])->shoulder_0 = params_4[int(29)];
    (&(&p_2)->crf_params_1[int(1)])->gamma_0 = params_4[int(30)];
    (&(&p_2)->crf_params_1[int(1)])->center_0 = params_4[int(31)];
    (&(&p_2)->crf_params_1[int(2)])->toe_0 = params_4[int(32)];
    (&(&p_2)->crf_params_1[int(2)])->shoulder_0 = params_4[int(33)];
    (&(&p_2)->crf_params_1[int(2)])->gamma_0 = params_4[int(34)];
    (&(&p_2)->crf_params_1[int(2)])->center_0 = params_4[int(35)];
    FixedArray<float, 22>  losses_3;
    losses_3[int(0)] = 0.0f;
    losses_3[int(1)] = 0.0f;
    losses_3[int(2)] = 0.0f;
    losses_3[int(3)] = 0.0f;
    losses_3[int(4)] = 0.0f;
    losses_3[int(5)] = 0.0f;
    losses_3[int(6)] = 0.0f;
    losses_3[int(7)] = 0.0f;
    losses_3[int(8)] = 0.0f;
    losses_3[int(9)] = 0.0f;
    losses_3[int(10)] = 0.0f;
    losses_3[int(11)] = 0.0f;
    losses_3[int(12)] = 0.0f;
    losses_3[int(13)] = 0.0f;
    losses_3[int(14)] = 0.0f;
    losses_3[int(15)] = 0.0f;
    losses_3[int(16)] = 0.0f;
    losses_3[int(17)] = 0.0f;
    losses_3[int(18)] = 0.0f;
    losses_3[int(19)] = 0.0f;
    losses_3[int(20)] = 0.0f;
    losses_3[int(21)] = 0.0f;
    losses_3[int(0)] = p_2.exposure_1;
    float _S2339 = p_2.vignette_params_1[int(0)].cx_0;
    float _S2340 = p_2.vignette_params_1[int(0)].cy_0;
    float _S2341 = p_2.vignette_params_1[int(1)].cx_0;
    float _S2342 = p_2.vignette_params_1[int(1)].cy_0;
    float _S2343 = p_2.vignette_params_1[int(2)].cx_0;
    float _S2344 = p_2.vignette_params_1[int(2)].cy_0;
    losses_3[int(1)] = _S2339 * _S2339 + _S2340 * _S2340 + _S2341 * _S2341 + _S2342 * _S2342 + _S2343 * _S2343 + _S2344 * _S2344;
    losses_3[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_3[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_3[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_7 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S2345 = p_2.vignette_params_1[int(0)].cx_0 - mean_7;
    float _S2346 = p_2.vignette_params_1[int(1)].cx_0 - mean_7;
    float _S2347 = p_2.vignette_params_1[int(2)].cx_0 - mean_7;
    losses_3[int(5)] = (_S2345 * _S2345 + _S2346 * _S2346 + _S2347 * _S2347) / 3.0f;
    float mean_8 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S2348 = p_2.vignette_params_1[int(0)].cy_0 - mean_8;
    float _S2349 = p_2.vignette_params_1[int(1)].cy_0 - mean_8;
    float _S2350 = p_2.vignette_params_1[int(2)].cy_0 - mean_8;
    losses_3[int(6)] = (_S2348 * _S2348 + _S2349 * _S2349 + _S2350 * _S2350) / 3.0f;
    float mean_9 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S2351 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_9;
    float _S2352 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_9;
    float _S2353 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_9;
    losses_3[int(7)] = (_S2351 * _S2351 + _S2352 * _S2352 + _S2353 * _S2353) / 3.0f;
    float mean_10 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S2354 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_10;
    float _S2355 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_10;
    float _S2356 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_10;
    losses_3[int(8)] = (_S2354 * _S2354 + _S2355 * _S2355 + _S2356 * _S2356) / 3.0f;
    float mean_11 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S2357 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_11;
    float _S2358 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_11;
    float _S2359 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_11;
    losses_3[int(9)] = (_S2357 * _S2357 + _S2358 * _S2358 + _S2359 * _S2359) / 3.0f;
    float2  bd_2 = mul_4(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_2.color_params_1.b_0);
    float2  rd_2 = mul_4(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_2.color_params_1.r_0);
    float2  gd_2 = mul_4(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_2.color_params_1.g_0);
    float2  nd_2 = mul_4(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_2.color_params_1.n_0);
    losses_3[int(10)] = bd_2.x;
    losses_3[int(11)] = bd_2.y;
    losses_3[int(12)] = rd_2.x;
    losses_3[int(13)] = rd_2.y;
    losses_3[int(14)] = gd_2.x;
    losses_3[int(15)] = gd_2.y;
    losses_3[int(16)] = nd_2.x;
    losses_3[int(17)] = nd_2.y;
    float mean_12 = (p_2.crf_params_1[int(0)].toe_0 + p_2.crf_params_1[int(1)].toe_0 + p_2.crf_params_1[int(2)].toe_0) / 3.0f;
    float _S2360 = p_2.crf_params_1[int(0)].toe_0 - mean_12;
    float _S2361 = p_2.crf_params_1[int(1)].toe_0 - mean_12;
    float _S2362 = p_2.crf_params_1[int(2)].toe_0 - mean_12;
    losses_3[int(18)] = (_S2360 * _S2360 + _S2361 * _S2361 + _S2362 * _S2362) / 3.0f;
    float mean_13 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S2363 = p_2.crf_params_1[int(0)].shoulder_0 - mean_13;
    float _S2364 = p_2.crf_params_1[int(1)].shoulder_0 - mean_13;
    float _S2365 = p_2.crf_params_1[int(2)].shoulder_0 - mean_13;
    losses_3[int(19)] = (_S2363 * _S2363 + _S2364 * _S2364 + _S2365 * _S2365) / 3.0f;
    float mean_14 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S2366 = p_2.crf_params_1[int(0)].gamma_0 - mean_14;
    float _S2367 = p_2.crf_params_1[int(1)].gamma_0 - mean_14;
    float _S2368 = p_2.crf_params_1[int(2)].gamma_0 - mean_14;
    losses_3[int(20)] = (_S2366 * _S2366 + _S2367 * _S2367 + _S2368 * _S2368) / 3.0f;
    float mean_15 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S2369 = p_2.crf_params_1[int(0)].center_0 - mean_15;
    float _S2370 = p_2.crf_params_1[int(1)].center_0 - mean_15;
    float _S2371 = p_2.crf_params_1[int(2)].center_0 - mean_15;
    losses_3[int(21)] = (_S2369 * _S2369 + _S2370 * _S2370 + _S2371 * _S2371) / 3.0f;
    *_S2338 = losses_3;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_5, FixedArray<float, 23>  * _S2372)
{
    PPISPParamsRQS_0 p_3;
    (&p_3)->exposure_0 = params_5[int(0)];
    (&(&p_3)->vignette_params_0[int(0)])->cx_0 = params_5[int(1)];
    (&(&p_3)->vignette_params_0[int(0)])->cy_0 = params_5[int(2)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha0_0 = params_5[int(3)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha1_0 = params_5[int(4)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha2_0 = params_5[int(5)];
    (&(&p_3)->vignette_params_0[int(1)])->cx_0 = params_5[int(6)];
    (&(&p_3)->vignette_params_0[int(1)])->cy_0 = params_5[int(7)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha0_0 = params_5[int(8)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha1_0 = params_5[int(9)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha2_0 = params_5[int(10)];
    (&(&p_3)->vignette_params_0[int(2)])->cx_0 = params_5[int(11)];
    (&(&p_3)->vignette_params_0[int(2)])->cy_0 = params_5[int(12)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha0_0 = params_5[int(13)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha1_0 = params_5[int(14)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha2_0 = params_5[int(15)];
    *&((&(&(&p_3)->color_params_0)->b_0)->x) = params_5[int(16)];
    *&((&(&(&p_3)->color_params_0)->b_0)->y) = params_5[int(17)];
    *&((&(&(&p_3)->color_params_0)->r_0)->x) = params_5[int(18)];
    *&((&(&(&p_3)->color_params_0)->r_0)->y) = params_5[int(19)];
    *&((&(&(&p_3)->color_params_0)->g_0)->x) = params_5[int(20)];
    *&((&(&(&p_3)->color_params_0)->g_0)->y) = params_5[int(21)];
    *&((&(&(&p_3)->color_params_0)->n_0)->x) = params_5[int(22)];
    *&((&(&(&p_3)->color_params_0)->n_0)->y) = params_5[int(23)];
    (&(&p_3)->crf_params_0[int(0)])->g0_0 = params_5[int(24)];
    (&(&p_3)->crf_params_0[int(0)])->g1_0 = params_5[int(25)];
    (&(&p_3)->crf_params_0[int(0)])->x0_0 = params_5[int(26)];
    (&(&p_3)->crf_params_0[int(0)])->y0_0 = params_5[int(27)];
    (&(&p_3)->crf_params_0[int(0)])->gc_0 = params_5[int(28)];
    (&(&p_3)->crf_params_0[int(1)])->g0_0 = params_5[int(29)];
    (&(&p_3)->crf_params_0[int(1)])->g1_0 = params_5[int(30)];
    (&(&p_3)->crf_params_0[int(1)])->x0_0 = params_5[int(31)];
    (&(&p_3)->crf_params_0[int(1)])->y0_0 = params_5[int(32)];
    (&(&p_3)->crf_params_0[int(1)])->gc_0 = params_5[int(33)];
    (&(&p_3)->crf_params_0[int(2)])->g0_0 = params_5[int(34)];
    (&(&p_3)->crf_params_0[int(2)])->g1_0 = params_5[int(35)];
    (&(&p_3)->crf_params_0[int(2)])->x0_0 = params_5[int(36)];
    (&(&p_3)->crf_params_0[int(2)])->y0_0 = params_5[int(37)];
    (&(&p_3)->crf_params_0[int(2)])->gc_0 = params_5[int(38)];
    FixedArray<float, 23>  losses_4;
    losses_4[int(0)] = 0.0f;
    losses_4[int(1)] = 0.0f;
    losses_4[int(2)] = 0.0f;
    losses_4[int(3)] = 0.0f;
    losses_4[int(4)] = 0.0f;
    losses_4[int(5)] = 0.0f;
    losses_4[int(6)] = 0.0f;
    losses_4[int(7)] = 0.0f;
    losses_4[int(8)] = 0.0f;
    losses_4[int(9)] = 0.0f;
    losses_4[int(10)] = 0.0f;
    losses_4[int(11)] = 0.0f;
    losses_4[int(12)] = 0.0f;
    losses_4[int(13)] = 0.0f;
    losses_4[int(14)] = 0.0f;
    losses_4[int(15)] = 0.0f;
    losses_4[int(16)] = 0.0f;
    losses_4[int(17)] = 0.0f;
    losses_4[int(18)] = 0.0f;
    losses_4[int(19)] = 0.0f;
    losses_4[int(20)] = 0.0f;
    losses_4[int(21)] = 0.0f;
    losses_4[int(22)] = 0.0f;
    losses_4[int(0)] = p_3.exposure_0;
    float _S2373 = p_3.vignette_params_0[int(0)].cx_0;
    float _S2374 = p_3.vignette_params_0[int(0)].cy_0;
    float _S2375 = p_3.vignette_params_0[int(1)].cx_0;
    float _S2376 = p_3.vignette_params_0[int(1)].cy_0;
    float _S2377 = p_3.vignette_params_0[int(2)].cx_0;
    float _S2378 = p_3.vignette_params_0[int(2)].cy_0;
    losses_4[int(1)] = _S2373 * _S2373 + _S2374 * _S2374 + _S2375 * _S2375 + _S2376 * _S2376 + _S2377 * _S2377 + _S2378 * _S2378;
    losses_4[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_4[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_4[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_16 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S2379 = p_3.vignette_params_0[int(0)].cx_0 - mean_16;
    float _S2380 = p_3.vignette_params_0[int(1)].cx_0 - mean_16;
    float _S2381 = p_3.vignette_params_0[int(2)].cx_0 - mean_16;
    losses_4[int(5)] = (_S2379 * _S2379 + _S2380 * _S2380 + _S2381 * _S2381) / 3.0f;
    float mean_17 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S2382 = p_3.vignette_params_0[int(0)].cy_0 - mean_17;
    float _S2383 = p_3.vignette_params_0[int(1)].cy_0 - mean_17;
    float _S2384 = p_3.vignette_params_0[int(2)].cy_0 - mean_17;
    losses_4[int(6)] = (_S2382 * _S2382 + _S2383 * _S2383 + _S2384 * _S2384) / 3.0f;
    float mean_18 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S2385 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_18;
    float _S2386 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_18;
    float _S2387 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_18;
    losses_4[int(7)] = (_S2385 * _S2385 + _S2386 * _S2386 + _S2387 * _S2387) / 3.0f;
    float mean_19 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S2388 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_19;
    float _S2389 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_19;
    float _S2390 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_19;
    losses_4[int(8)] = (_S2388 * _S2388 + _S2389 * _S2389 + _S2390 * _S2390) / 3.0f;
    float mean_20 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S2391 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_20;
    float _S2392 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_20;
    float _S2393 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_20;
    losses_4[int(9)] = (_S2391 * _S2391 + _S2392 * _S2392 + _S2393 * _S2393) / 3.0f;
    float2  bd_3 = mul_4(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_3.color_params_0.b_0);
    float2  rd_3 = mul_4(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_3.color_params_0.r_0);
    float2  gd_3 = mul_4(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_3.color_params_0.g_0);
    float2  nd_3 = mul_4(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_3.color_params_0.n_0);
    losses_4[int(10)] = bd_3.x;
    losses_4[int(11)] = bd_3.y;
    losses_4[int(12)] = rd_3.x;
    losses_4[int(13)] = rd_3.y;
    losses_4[int(14)] = gd_3.x;
    losses_4[int(15)] = gd_3.y;
    losses_4[int(16)] = nd_3.x;
    losses_4[int(17)] = nd_3.y;
    float mean_21 = (p_3.crf_params_0[int(0)].g0_0 + p_3.crf_params_0[int(1)].g0_0 + p_3.crf_params_0[int(2)].g0_0) / 3.0f;
    float _S2394 = p_3.crf_params_0[int(0)].g0_0 - mean_21;
    float _S2395 = p_3.crf_params_0[int(1)].g0_0 - mean_21;
    float _S2396 = p_3.crf_params_0[int(2)].g0_0 - mean_21;
    losses_4[int(18)] = (_S2394 * _S2394 + _S2395 * _S2395 + _S2396 * _S2396) / 3.0f;
    float mean_22 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S2397 = p_3.crf_params_0[int(0)].g1_0 - mean_22;
    float _S2398 = p_3.crf_params_0[int(1)].g1_0 - mean_22;
    float _S2399 = p_3.crf_params_0[int(2)].g1_0 - mean_22;
    losses_4[int(19)] = (_S2397 * _S2397 + _S2398 * _S2398 + _S2399 * _S2399) / 3.0f;
    float mean_23 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S2400 = p_3.crf_params_0[int(0)].x0_0 - mean_23;
    float _S2401 = p_3.crf_params_0[int(1)].x0_0 - mean_23;
    float _S2402 = p_3.crf_params_0[int(2)].x0_0 - mean_23;
    losses_4[int(20)] = (_S2400 * _S2400 + _S2401 * _S2401 + _S2402 * _S2402) / 3.0f;
    float mean_24 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S2403 = p_3.crf_params_0[int(0)].y0_0 - mean_24;
    float _S2404 = p_3.crf_params_0[int(1)].y0_0 - mean_24;
    float _S2405 = p_3.crf_params_0[int(2)].y0_0 - mean_24;
    losses_4[int(21)] = (_S2403 * _S2403 + _S2404 * _S2404 + _S2405 * _S2405) / 3.0f;
    float mean_25 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S2406 = p_3.crf_params_0[int(0)].gc_0 - mean_25;
    float _S2407 = p_3.crf_params_0[int(1)].gc_0 - mean_25;
    float _S2408 = p_3.crf_params_0[int(2)].gc_0 - mean_25;
    losses_4[int(22)] = (_S2406 * _S2406 + _S2407 * _S2407 + _S2408 * _S2408) / 3.0f;
    *_S2372 = losses_4;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_15)
{
    VignettingChannelParams_0 _S2409 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2410 = {
        _S2409, _S2409, _S2409
    };
    float2  _S2411 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2412 = { _S2411, _S2411, _S2411, _S2411 };
    CRFPPISPChannelParams_0 _S2413 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S2414 = {
        _S2413, _S2413, _S2413
    };
    PPISPParams_0 _S2415;
    (&_S2415)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S2415)->vignette_params_1 = _S2410;
    (&_S2415)->color_params_1 = _S2412;
    (&_S2415)->crf_params_1 = _S2414;
    (&(&_S2415)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S2415)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S2415)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S2415)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S2415)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S2415)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S2415)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S2415)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S2415)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S2415)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S2415)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S2415)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S2415)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S2415)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S2415)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S2415)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S2415)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S2415)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S2415)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S2415)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S2415)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S2415)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S2415)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S2415)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S2415)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S2415)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S2415)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S2415)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S2415)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S2415)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S2415)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S2415)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S2415)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S2415)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S2415)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_26 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S2416 = dpparams_2->primal_0[int(1)] - mean_26;
    float _S2417 = dpparams_2->primal_0[int(6)] - mean_26;
    float _S2418 = dpparams_2->primal_0[int(11)] - mean_26;
    float mean_27 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S2419 = dpparams_2->primal_0[int(2)] - mean_27;
    float _S2420 = dpparams_2->primal_0[int(7)] - mean_27;
    float _S2421 = dpparams_2->primal_0[int(12)] - mean_27;
    float mean_28 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S2422 = dpparams_2->primal_0[int(3)] - mean_28;
    float _S2423 = dpparams_2->primal_0[int(8)] - mean_28;
    float _S2424 = dpparams_2->primal_0[int(13)] - mean_28;
    float mean_29 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S2425 = dpparams_2->primal_0[int(4)] - mean_29;
    float _S2426 = dpparams_2->primal_0[int(9)] - mean_29;
    float _S2427 = dpparams_2->primal_0[int(14)] - mean_29;
    float mean_30 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S2428 = dpparams_2->primal_0[int(5)] - mean_30;
    float _S2429 = dpparams_2->primal_0[int(10)] - mean_30;
    float _S2430 = dpparams_2->primal_0[int(15)] - mean_30;
    float mean_31 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_32 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_33 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_34 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S2431 = 0.3333333432674408f * (*_s_dOut_15)[int(21)];
    float _S2432 = (dpparams_2->primal_0[int(35)] - mean_34) * _S2431;
    float _S2433 = _S2432 + _S2432;
    float _S2434 = (dpparams_2->primal_0[int(31)] - mean_34) * _S2431;
    float _S2435 = _S2434 + _S2434;
    float _S2436 = (dpparams_2->primal_0[int(27)] - mean_34) * _S2431;
    float _S2437 = _S2436 + _S2436;
    float _S2438 = 0.3333333432674408f * (- _S2433 + - _S2435 + - _S2437);
    float _S2439 = 0.3333333432674408f * (*_s_dOut_15)[int(20)];
    float _S2440 = (dpparams_2->primal_0[int(34)] - mean_33) * _S2439;
    float _S2441 = _S2440 + _S2440;
    float _S2442 = (dpparams_2->primal_0[int(30)] - mean_33) * _S2439;
    float _S2443 = _S2442 + _S2442;
    float _S2444 = (dpparams_2->primal_0[int(26)] - mean_33) * _S2439;
    float _S2445 = _S2444 + _S2444;
    float _S2446 = 0.3333333432674408f * (- _S2441 + - _S2443 + - _S2445);
    float _S2447 = 0.3333333432674408f * (*_s_dOut_15)[int(19)];
    float _S2448 = (dpparams_2->primal_0[int(33)] - mean_32) * _S2447;
    float _S2449 = _S2448 + _S2448;
    float _S2450 = (dpparams_2->primal_0[int(29)] - mean_32) * _S2447;
    float _S2451 = _S2450 + _S2450;
    float _S2452 = (dpparams_2->primal_0[int(25)] - mean_32) * _S2447;
    float _S2453 = _S2452 + _S2452;
    float _S2454 = 0.3333333432674408f * (- _S2449 + - _S2451 + - _S2453);
    float _S2455 = 0.3333333432674408f * (*_s_dOut_15)[int(18)];
    float _S2456 = (dpparams_2->primal_0[int(32)] - mean_31) * _S2455;
    float _S2457 = _S2456 + _S2456;
    float _S2458 = (dpparams_2->primal_0[int(28)] - mean_31) * _S2455;
    float _S2459 = _S2458 + _S2458;
    float _S2460 = (dpparams_2->primal_0[int(24)] - mean_31) * _S2455;
    float _S2461 = _S2460 + _S2460;
    float _S2462 = 0.3333333432674408f * (- _S2457 + - _S2459 + - _S2461);
    float2  _S2463 = make_float2 ((*_s_dOut_15)[int(16)], (*_s_dOut_15)[int(17)]);
    Matrix<float, 2, 2>  _S2464 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2465;
    (&_S2465)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2465)->differential_0 = _S2464;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2466;
    (&_S2466)->primal_0 = _S2415.color_params_1.n_0;
    (&_S2466)->differential_0 = _S2411;
    s_bwd_prop_mul_3(&_S2465, &_S2466, _S2463);
    float2  _S2467 = make_float2 ((*_s_dOut_15)[int(14)], (*_s_dOut_15)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2468;
    (&_S2468)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2468)->differential_0 = _S2464;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2469;
    (&_S2469)->primal_0 = _S2415.color_params_1.g_0;
    (&_S2469)->differential_0 = _S2411;
    s_bwd_prop_mul_3(&_S2468, &_S2469, _S2467);
    float2  _S2470 = make_float2 ((*_s_dOut_15)[int(12)], (*_s_dOut_15)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2471;
    (&_S2471)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2471)->differential_0 = _S2464;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2472;
    (&_S2472)->primal_0 = _S2415.color_params_1.r_0;
    (&_S2472)->differential_0 = _S2411;
    s_bwd_prop_mul_3(&_S2471, &_S2472, _S2470);
    float2  _S2473 = make_float2 ((*_s_dOut_15)[int(10)], (*_s_dOut_15)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2474;
    (&_S2474)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2474)->differential_0 = _S2464;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2475;
    (&_S2475)->primal_0 = _S2415.color_params_1.b_0;
    (&_S2475)->differential_0 = _S2411;
    s_bwd_prop_mul_3(&_S2474, &_S2475, _S2473);
    ColorPPISPParams_0 _S2476 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2476)->n_0 = _S2466.differential_0;
    (&_S2476)->g_0 = _S2469.differential_0;
    (&_S2476)->r_0 = _S2472.differential_0;
    (&_S2476)->b_0 = _S2475.differential_0;
    float _S2477 = 0.3333333432674408f * (*_s_dOut_15)[int(9)];
    float _S2478 = _S2430 * _S2477;
    float _S2479 = _S2478 + _S2478;
    float _S2480 = _S2429 * _S2477;
    float _S2481 = _S2480 + _S2480;
    float _S2482 = _S2428 * _S2477;
    float _S2483 = _S2482 + _S2482;
    float _S2484 = 0.3333333432674408f * (- _S2479 + - _S2481 + - _S2483);
    float _S2485 = 0.3333333432674408f * (*_s_dOut_15)[int(8)];
    float _S2486 = _S2427 * _S2485;
    float _S2487 = _S2486 + _S2486;
    float _S2488 = _S2426 * _S2485;
    float _S2489 = _S2488 + _S2488;
    float _S2490 = _S2425 * _S2485;
    float _S2491 = _S2490 + _S2490;
    float _S2492 = 0.3333333432674408f * (- _S2487 + - _S2489 + - _S2491);
    float _S2493 = 0.3333333432674408f * (*_s_dOut_15)[int(7)];
    float _S2494 = _S2424 * _S2493;
    float _S2495 = _S2494 + _S2494;
    float _S2496 = _S2423 * _S2493;
    float _S2497 = _S2496 + _S2496;
    float _S2498 = _S2422 * _S2493;
    float _S2499 = _S2498 + _S2498;
    float _S2500 = 0.3333333432674408f * (- _S2495 + - _S2497 + - _S2499);
    float _S2501 = 0.3333333432674408f * (*_s_dOut_15)[int(6)];
    float _S2502 = _S2421 * _S2501;
    float _S2503 = _S2502 + _S2502;
    float _S2504 = _S2420 * _S2501;
    float _S2505 = _S2504 + _S2504;
    float _S2506 = _S2419 * _S2501;
    float _S2507 = _S2506 + _S2506;
    float _S2508 = 0.3333333432674408f * (- _S2503 + - _S2505 + - _S2507);
    float _S2509 = 0.3333333432674408f * (*_s_dOut_15)[int(5)];
    float _S2510 = _S2418 * _S2509;
    float _S2511 = _S2510 + _S2510;
    float _S2512 = _S2417 * _S2509;
    float _S2513 = _S2512 + _S2512;
    float _S2514 = _S2416 * _S2509;
    float _S2515 = _S2514 + _S2514;
    float _S2516 = 0.3333333432674408f * (- _S2511 + - _S2513 + - _S2515);
    DiffPair_float_0 _S2517;
    (&_S2517)->primal_0 = 0.0f;
    (&_S2517)->differential_0 = 0.0f;
    DiffPair_float_0 _S2518;
    (&_S2518)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S2518)->differential_0 = 0.0f;
    _d_max_0(&_S2517, &_S2518, (*_s_dOut_15)[int(4)]);
    DiffPair_float_0 _S2519;
    (&_S2519)->primal_0 = 0.0f;
    (&_S2519)->differential_0 = 0.0f;
    DiffPair_float_0 _S2520;
    (&_S2520)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S2520)->differential_0 = 0.0f;
    _d_max_0(&_S2519, &_S2520, (*_s_dOut_15)[int(4)]);
    DiffPair_float_0 _S2521;
    (&_S2521)->primal_0 = 0.0f;
    (&_S2521)->differential_0 = 0.0f;
    DiffPair_float_0 _S2522;
    (&_S2522)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S2522)->differential_0 = 0.0f;
    _d_max_0(&_S2521, &_S2522, (*_s_dOut_15)[int(4)]);
    DiffPair_float_0 _S2523;
    (&_S2523)->primal_0 = 0.0f;
    (&_S2523)->differential_0 = 0.0f;
    DiffPair_float_0 _S2524;
    (&_S2524)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S2524)->differential_0 = 0.0f;
    _d_max_0(&_S2523, &_S2524, (*_s_dOut_15)[int(3)]);
    DiffPair_float_0 _S2525;
    (&_S2525)->primal_0 = 0.0f;
    (&_S2525)->differential_0 = 0.0f;
    DiffPair_float_0 _S2526;
    (&_S2526)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S2526)->differential_0 = 0.0f;
    _d_max_0(&_S2525, &_S2526, (*_s_dOut_15)[int(3)]);
    DiffPair_float_0 _S2527;
    (&_S2527)->primal_0 = 0.0f;
    (&_S2527)->differential_0 = 0.0f;
    DiffPair_float_0 _S2528;
    (&_S2528)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S2528)->differential_0 = 0.0f;
    _d_max_0(&_S2527, &_S2528, (*_s_dOut_15)[int(3)]);
    DiffPair_float_0 _S2529;
    (&_S2529)->primal_0 = 0.0f;
    (&_S2529)->differential_0 = 0.0f;
    DiffPair_float_0 _S2530;
    (&_S2530)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S2530)->differential_0 = 0.0f;
    _d_max_0(&_S2529, &_S2530, (*_s_dOut_15)[int(2)]);
    DiffPair_float_0 _S2531;
    (&_S2531)->primal_0 = 0.0f;
    (&_S2531)->differential_0 = 0.0f;
    DiffPair_float_0 _S2532;
    (&_S2532)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S2532)->differential_0 = 0.0f;
    _d_max_0(&_S2531, &_S2532, (*_s_dOut_15)[int(2)]);
    DiffPair_float_0 _S2533;
    (&_S2533)->primal_0 = 0.0f;
    (&_S2533)->differential_0 = 0.0f;
    DiffPair_float_0 _S2534;
    (&_S2534)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S2534)->differential_0 = 0.0f;
    _d_max_0(&_S2533, &_S2534, (*_s_dOut_15)[int(2)]);
    float _S2535 = dpparams_2->primal_0[int(12)] * (*_s_dOut_15)[int(1)];
    float _S2536 = dpparams_2->primal_0[int(11)] * (*_s_dOut_15)[int(1)];
    float _S2537 = dpparams_2->primal_0[int(7)] * (*_s_dOut_15)[int(1)];
    float _S2538 = dpparams_2->primal_0[int(6)] * (*_s_dOut_15)[int(1)];
    float _S2539 = dpparams_2->primal_0[int(2)] * (*_s_dOut_15)[int(1)];
    float _S2540 = dpparams_2->primal_0[int(1)] * (*_s_dOut_15)[int(1)];
    PPISPParams_0 _S2541 = PPISPParams_x24_syn_dzero_0();
    (&_S2541)->color_params_1 = _S2476;
    (&_S2541)->exposure_1 = (*_s_dOut_15)[int(0)];
    _S2415 = _S2541;
    (&(&_S2415)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2542 = _S2433 + _S2438 + _S2541.crf_params_1[int(2)].center_0;
    (&(&_S2415)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2543 = _S2441 + _S2446 + _S2541.crf_params_1[int(2)].gamma_0;
    (&(&_S2415)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2544 = _S2449 + _S2454 + _S2541.crf_params_1[int(2)].shoulder_0;
    (&(&_S2415)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2545 = _S2457 + _S2462 + _S2541.crf_params_1[int(2)].toe_0;
    (&(&_S2415)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2546 = _S2435 + _S2438 + _S2541.crf_params_1[int(1)].center_0;
    (&(&_S2415)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2547 = _S2443 + _S2446 + _S2541.crf_params_1[int(1)].gamma_0;
    (&(&_S2415)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2548 = _S2451 + _S2454 + _S2541.crf_params_1[int(1)].shoulder_0;
    (&(&_S2415)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2549 = _S2459 + _S2462 + _S2541.crf_params_1[int(1)].toe_0;
    (&(&_S2415)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2550 = _S2437 + _S2438 + _S2541.crf_params_1[int(0)].center_0;
    (&(&_S2415)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2551 = _S2445 + _S2446 + _S2541.crf_params_1[int(0)].gamma_0;
    (&(&_S2415)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2552 = _S2453 + _S2454 + _S2541.crf_params_1[int(0)].shoulder_0;
    (&(&_S2415)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2553 = _S2461 + _S2462 + _S2541.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S2415)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S2415)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S2415)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2554 = _S2479 + _S2484 + _S2518.differential_0 + _S2541.vignette_params_1[int(2)].alpha2_0;
    (&(&_S2415)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2555 = _S2487 + _S2492 + _S2524.differential_0 + _S2541.vignette_params_1[int(2)].alpha1_0;
    (&(&_S2415)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2556 = _S2495 + _S2500 + _S2530.differential_0 + _S2541.vignette_params_1[int(2)].alpha0_0;
    (&(&_S2415)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2557 = _S2503 + _S2508 + _S2535 + _S2535 + _S2541.vignette_params_1[int(2)].cy_0;
    (&(&_S2415)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2558 = _S2511 + _S2516 + _S2536 + _S2536 + _S2541.vignette_params_1[int(2)].cx_0;
    (&(&_S2415)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2559 = _S2481 + _S2484 + _S2520.differential_0 + _S2541.vignette_params_1[int(1)].alpha2_0;
    (&(&_S2415)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2560 = _S2489 + _S2492 + _S2526.differential_0 + _S2541.vignette_params_1[int(1)].alpha1_0;
    (&(&_S2415)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2561 = _S2497 + _S2500 + _S2532.differential_0 + _S2541.vignette_params_1[int(1)].alpha0_0;
    (&(&_S2415)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2562 = _S2505 + _S2508 + _S2537 + _S2537 + _S2541.vignette_params_1[int(1)].cy_0;
    (&(&_S2415)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2563 = _S2513 + _S2516 + _S2538 + _S2538 + _S2541.vignette_params_1[int(1)].cx_0;
    (&(&_S2415)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2564 = _S2483 + _S2484 + _S2522.differential_0 + _S2541.vignette_params_1[int(0)].alpha2_0;
    (&(&_S2415)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2565 = _S2491 + _S2492 + _S2528.differential_0 + _S2541.vignette_params_1[int(0)].alpha1_0;
    (&(&_S2415)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2566 = _S2499 + _S2500 + _S2534.differential_0 + _S2541.vignette_params_1[int(0)].alpha0_0;
    (&(&_S2415)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2567 = _S2507 + _S2508 + _S2539 + _S2539 + _S2541.vignette_params_1[int(0)].cy_0;
    (&(&_S2415)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2568 = _S2515 + _S2516 + _S2540 + _S2540 + _S2541.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2569;
    _S2569[int(0)] = 0.0f;
    _S2569[int(1)] = 0.0f;
    _S2569[int(2)] = 0.0f;
    _S2569[int(3)] = 0.0f;
    _S2569[int(4)] = 0.0f;
    _S2569[int(5)] = 0.0f;
    _S2569[int(6)] = 0.0f;
    _S2569[int(7)] = 0.0f;
    _S2569[int(8)] = 0.0f;
    _S2569[int(9)] = 0.0f;
    _S2569[int(10)] = 0.0f;
    _S2569[int(11)] = 0.0f;
    _S2569[int(12)] = 0.0f;
    _S2569[int(13)] = 0.0f;
    _S2569[int(14)] = 0.0f;
    _S2569[int(15)] = 0.0f;
    _S2569[int(16)] = 0.0f;
    _S2569[int(17)] = 0.0f;
    _S2569[int(18)] = 0.0f;
    _S2569[int(19)] = 0.0f;
    _S2569[int(20)] = 0.0f;
    _S2569[int(21)] = 0.0f;
    _S2569[int(22)] = 0.0f;
    _S2569[int(23)] = 0.0f;
    _S2569[int(24)] = 0.0f;
    _S2569[int(25)] = 0.0f;
    _S2569[int(26)] = 0.0f;
    _S2569[int(27)] = 0.0f;
    _S2569[int(28)] = 0.0f;
    _S2569[int(29)] = 0.0f;
    _S2569[int(30)] = 0.0f;
    _S2569[int(31)] = 0.0f;
    _S2569[int(32)] = 0.0f;
    _S2569[int(33)] = 0.0f;
    _S2569[int(34)] = 0.0f;
    _S2569[int(35)] = 0.0f;
    _S2569[int(8)] = _S2561;
    _S2569[int(16)] = _S2541.color_params_1.b_0.x;
    _S2569[int(15)] = _S2554;
    _S2569[int(14)] = _S2555;
    _S2569[int(13)] = _S2556;
    _S2569[int(12)] = _S2557;
    _S2569[int(11)] = _S2558;
    _S2569[int(10)] = _S2559;
    _S2569[int(9)] = _S2560;
    _S2569[int(17)] = _S2541.color_params_1.b_0.y;
    _S2569[int(7)] = _S2562;
    _S2569[int(6)] = _S2563;
    _S2569[int(5)] = _S2564;
    _S2569[int(4)] = _S2565;
    _S2569[int(3)] = _S2566;
    _S2569[int(2)] = _S2567;
    _S2569[int(1)] = _S2568;
    _S2569[int(0)] = _S2415.exposure_1;
    _S2569[int(26)] = _S2551;
    _S2569[int(34)] = _S2543;
    _S2569[int(33)] = _S2544;
    _S2569[int(32)] = _S2545;
    _S2569[int(31)] = _S2546;
    _S2569[int(30)] = _S2547;
    _S2569[int(29)] = _S2548;
    _S2569[int(28)] = _S2549;
    _S2569[int(27)] = _S2550;
    _S2569[int(35)] = _S2542;
    _S2569[int(25)] = _S2552;
    _S2569[int(24)] = _S2553;
    _S2569[int(23)] = _S2541.color_params_1.n_0.y;
    _S2569[int(22)] = _S2541.color_params_1.n_0.x;
    _S2569[int(21)] = _S2541.color_params_1.g_0.y;
    _S2569[int(20)] = _S2541.color_params_1.g_0.x;
    _S2569[int(19)] = _S2541.color_params_1.r_0.y;
    _S2569[int(18)] = _S2541.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S2569;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2570, FixedArray<float, 22>  * _S2571)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S2570, _S2571);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_6, FixedArray<float, 22>  grad_out_2, FixedArray<float, 36>  * _S2572)
{
    FixedArray<float, 36>  _S2573 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_6;
    (&dp_params_2)->differential_0 = _S2573;
    FixedArray<float, 22>  _S2574 = grad_out_2;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, &_S2574);
    *_S2572 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_16)
{
    VignettingChannelParams_0 _S2575 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2576 = {
        _S2575, _S2575, _S2575
    };
    float2  _S2577 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2578 = { _S2577, _S2577, _S2577, _S2577 };
    RQSCRFPPISPChannelParams_0 _S2579 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2580 = {
        _S2579, _S2579, _S2579
    };
    PPISPParamsRQS_0 _S2581;
    (&_S2581)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S2581)->vignette_params_0 = _S2576;
    (&_S2581)->color_params_0 = _S2578;
    (&_S2581)->crf_params_0 = _S2580;
    (&(&_S2581)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S2581)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S2581)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S2581)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S2581)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S2581)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S2581)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S2581)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S2581)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S2581)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S2581)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S2581)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S2581)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S2581)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S2581)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S2581)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S2581)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S2581)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S2581)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S2581)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S2581)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S2581)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S2581)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S2581)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S2581)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S2581)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S2581)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S2581)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S2581)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S2581)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S2581)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S2581)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S2581)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S2581)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S2581)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S2581)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S2581)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S2581)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_35 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S2582 = dpparams_3->primal_0[int(1)] - mean_35;
    float _S2583 = dpparams_3->primal_0[int(6)] - mean_35;
    float _S2584 = dpparams_3->primal_0[int(11)] - mean_35;
    float mean_36 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S2585 = dpparams_3->primal_0[int(2)] - mean_36;
    float _S2586 = dpparams_3->primal_0[int(7)] - mean_36;
    float _S2587 = dpparams_3->primal_0[int(12)] - mean_36;
    float mean_37 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S2588 = dpparams_3->primal_0[int(3)] - mean_37;
    float _S2589 = dpparams_3->primal_0[int(8)] - mean_37;
    float _S2590 = dpparams_3->primal_0[int(13)] - mean_37;
    float mean_38 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S2591 = dpparams_3->primal_0[int(4)] - mean_38;
    float _S2592 = dpparams_3->primal_0[int(9)] - mean_38;
    float _S2593 = dpparams_3->primal_0[int(14)] - mean_38;
    float mean_39 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S2594 = dpparams_3->primal_0[int(5)] - mean_39;
    float _S2595 = dpparams_3->primal_0[int(10)] - mean_39;
    float _S2596 = dpparams_3->primal_0[int(15)] - mean_39;
    float mean_40 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_41 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_42 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_43 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_44 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S2597 = 0.3333333432674408f * (*_s_dOut_16)[int(22)];
    float _S2598 = (dpparams_3->primal_0[int(38)] - mean_44) * _S2597;
    float _S2599 = _S2598 + _S2598;
    float _S2600 = (dpparams_3->primal_0[int(33)] - mean_44) * _S2597;
    float _S2601 = _S2600 + _S2600;
    float _S2602 = (dpparams_3->primal_0[int(28)] - mean_44) * _S2597;
    float _S2603 = _S2602 + _S2602;
    float _S2604 = 0.3333333432674408f * (- _S2599 + - _S2601 + - _S2603);
    float _S2605 = 0.3333333432674408f * (*_s_dOut_16)[int(21)];
    float _S2606 = (dpparams_3->primal_0[int(37)] - mean_43) * _S2605;
    float _S2607 = _S2606 + _S2606;
    float _S2608 = (dpparams_3->primal_0[int(32)] - mean_43) * _S2605;
    float _S2609 = _S2608 + _S2608;
    float _S2610 = (dpparams_3->primal_0[int(27)] - mean_43) * _S2605;
    float _S2611 = _S2610 + _S2610;
    float _S2612 = 0.3333333432674408f * (- _S2607 + - _S2609 + - _S2611);
    float _S2613 = 0.3333333432674408f * (*_s_dOut_16)[int(20)];
    float _S2614 = (dpparams_3->primal_0[int(36)] - mean_42) * _S2613;
    float _S2615 = _S2614 + _S2614;
    float _S2616 = (dpparams_3->primal_0[int(31)] - mean_42) * _S2613;
    float _S2617 = _S2616 + _S2616;
    float _S2618 = (dpparams_3->primal_0[int(26)] - mean_42) * _S2613;
    float _S2619 = _S2618 + _S2618;
    float _S2620 = 0.3333333432674408f * (- _S2615 + - _S2617 + - _S2619);
    float _S2621 = 0.3333333432674408f * (*_s_dOut_16)[int(19)];
    float _S2622 = (dpparams_3->primal_0[int(35)] - mean_41) * _S2621;
    float _S2623 = _S2622 + _S2622;
    float _S2624 = (dpparams_3->primal_0[int(30)] - mean_41) * _S2621;
    float _S2625 = _S2624 + _S2624;
    float _S2626 = (dpparams_3->primal_0[int(25)] - mean_41) * _S2621;
    float _S2627 = _S2626 + _S2626;
    float _S2628 = 0.3333333432674408f * (- _S2623 + - _S2625 + - _S2627);
    float _S2629 = 0.3333333432674408f * (*_s_dOut_16)[int(18)];
    float _S2630 = (dpparams_3->primal_0[int(34)] - mean_40) * _S2629;
    float _S2631 = _S2630 + _S2630;
    float _S2632 = (dpparams_3->primal_0[int(29)] - mean_40) * _S2629;
    float _S2633 = _S2632 + _S2632;
    float _S2634 = (dpparams_3->primal_0[int(24)] - mean_40) * _S2629;
    float _S2635 = _S2634 + _S2634;
    float _S2636 = 0.3333333432674408f * (- _S2631 + - _S2633 + - _S2635);
    float2  _S2637 = make_float2 ((*_s_dOut_16)[int(16)], (*_s_dOut_16)[int(17)]);
    Matrix<float, 2, 2>  _S2638 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2639;
    (&_S2639)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2639)->differential_0 = _S2638;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2640;
    (&_S2640)->primal_0 = _S2581.color_params_0.n_0;
    (&_S2640)->differential_0 = _S2577;
    s_bwd_prop_mul_3(&_S2639, &_S2640, _S2637);
    float2  _S2641 = make_float2 ((*_s_dOut_16)[int(14)], (*_s_dOut_16)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2642;
    (&_S2642)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2642)->differential_0 = _S2638;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2643;
    (&_S2643)->primal_0 = _S2581.color_params_0.g_0;
    (&_S2643)->differential_0 = _S2577;
    s_bwd_prop_mul_3(&_S2642, &_S2643, _S2641);
    float2  _S2644 = make_float2 ((*_s_dOut_16)[int(12)], (*_s_dOut_16)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2645;
    (&_S2645)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2645)->differential_0 = _S2638;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2646;
    (&_S2646)->primal_0 = _S2581.color_params_0.r_0;
    (&_S2646)->differential_0 = _S2577;
    s_bwd_prop_mul_3(&_S2645, &_S2646, _S2644);
    float2  _S2647 = make_float2 ((*_s_dOut_16)[int(10)], (*_s_dOut_16)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2648;
    (&_S2648)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2648)->differential_0 = _S2638;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2649;
    (&_S2649)->primal_0 = _S2581.color_params_0.b_0;
    (&_S2649)->differential_0 = _S2577;
    s_bwd_prop_mul_3(&_S2648, &_S2649, _S2647);
    ColorPPISPParams_0 _S2650 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2650)->n_0 = _S2640.differential_0;
    (&_S2650)->g_0 = _S2643.differential_0;
    (&_S2650)->r_0 = _S2646.differential_0;
    (&_S2650)->b_0 = _S2649.differential_0;
    float _S2651 = 0.3333333432674408f * (*_s_dOut_16)[int(9)];
    float _S2652 = _S2596 * _S2651;
    float _S2653 = _S2652 + _S2652;
    float _S2654 = _S2595 * _S2651;
    float _S2655 = _S2654 + _S2654;
    float _S2656 = _S2594 * _S2651;
    float _S2657 = _S2656 + _S2656;
    float _S2658 = 0.3333333432674408f * (- _S2653 + - _S2655 + - _S2657);
    float _S2659 = 0.3333333432674408f * (*_s_dOut_16)[int(8)];
    float _S2660 = _S2593 * _S2659;
    float _S2661 = _S2660 + _S2660;
    float _S2662 = _S2592 * _S2659;
    float _S2663 = _S2662 + _S2662;
    float _S2664 = _S2591 * _S2659;
    float _S2665 = _S2664 + _S2664;
    float _S2666 = 0.3333333432674408f * (- _S2661 + - _S2663 + - _S2665);
    float _S2667 = 0.3333333432674408f * (*_s_dOut_16)[int(7)];
    float _S2668 = _S2590 * _S2667;
    float _S2669 = _S2668 + _S2668;
    float _S2670 = _S2589 * _S2667;
    float _S2671 = _S2670 + _S2670;
    float _S2672 = _S2588 * _S2667;
    float _S2673 = _S2672 + _S2672;
    float _S2674 = 0.3333333432674408f * (- _S2669 + - _S2671 + - _S2673);
    float _S2675 = 0.3333333432674408f * (*_s_dOut_16)[int(6)];
    float _S2676 = _S2587 * _S2675;
    float _S2677 = _S2676 + _S2676;
    float _S2678 = _S2586 * _S2675;
    float _S2679 = _S2678 + _S2678;
    float _S2680 = _S2585 * _S2675;
    float _S2681 = _S2680 + _S2680;
    float _S2682 = 0.3333333432674408f * (- _S2677 + - _S2679 + - _S2681);
    float _S2683 = 0.3333333432674408f * (*_s_dOut_16)[int(5)];
    float _S2684 = _S2584 * _S2683;
    float _S2685 = _S2684 + _S2684;
    float _S2686 = _S2583 * _S2683;
    float _S2687 = _S2686 + _S2686;
    float _S2688 = _S2582 * _S2683;
    float _S2689 = _S2688 + _S2688;
    float _S2690 = 0.3333333432674408f * (- _S2685 + - _S2687 + - _S2689);
    DiffPair_float_0 _S2691;
    (&_S2691)->primal_0 = 0.0f;
    (&_S2691)->differential_0 = 0.0f;
    DiffPair_float_0 _S2692;
    (&_S2692)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S2692)->differential_0 = 0.0f;
    _d_max_0(&_S2691, &_S2692, (*_s_dOut_16)[int(4)]);
    DiffPair_float_0 _S2693;
    (&_S2693)->primal_0 = 0.0f;
    (&_S2693)->differential_0 = 0.0f;
    DiffPair_float_0 _S2694;
    (&_S2694)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S2694)->differential_0 = 0.0f;
    _d_max_0(&_S2693, &_S2694, (*_s_dOut_16)[int(4)]);
    DiffPair_float_0 _S2695;
    (&_S2695)->primal_0 = 0.0f;
    (&_S2695)->differential_0 = 0.0f;
    DiffPair_float_0 _S2696;
    (&_S2696)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S2696)->differential_0 = 0.0f;
    _d_max_0(&_S2695, &_S2696, (*_s_dOut_16)[int(4)]);
    DiffPair_float_0 _S2697;
    (&_S2697)->primal_0 = 0.0f;
    (&_S2697)->differential_0 = 0.0f;
    DiffPair_float_0 _S2698;
    (&_S2698)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S2698)->differential_0 = 0.0f;
    _d_max_0(&_S2697, &_S2698, (*_s_dOut_16)[int(3)]);
    DiffPair_float_0 _S2699;
    (&_S2699)->primal_0 = 0.0f;
    (&_S2699)->differential_0 = 0.0f;
    DiffPair_float_0 _S2700;
    (&_S2700)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S2700)->differential_0 = 0.0f;
    _d_max_0(&_S2699, &_S2700, (*_s_dOut_16)[int(3)]);
    DiffPair_float_0 _S2701;
    (&_S2701)->primal_0 = 0.0f;
    (&_S2701)->differential_0 = 0.0f;
    DiffPair_float_0 _S2702;
    (&_S2702)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S2702)->differential_0 = 0.0f;
    _d_max_0(&_S2701, &_S2702, (*_s_dOut_16)[int(3)]);
    DiffPair_float_0 _S2703;
    (&_S2703)->primal_0 = 0.0f;
    (&_S2703)->differential_0 = 0.0f;
    DiffPair_float_0 _S2704;
    (&_S2704)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S2704)->differential_0 = 0.0f;
    _d_max_0(&_S2703, &_S2704, (*_s_dOut_16)[int(2)]);
    DiffPair_float_0 _S2705;
    (&_S2705)->primal_0 = 0.0f;
    (&_S2705)->differential_0 = 0.0f;
    DiffPair_float_0 _S2706;
    (&_S2706)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S2706)->differential_0 = 0.0f;
    _d_max_0(&_S2705, &_S2706, (*_s_dOut_16)[int(2)]);
    DiffPair_float_0 _S2707;
    (&_S2707)->primal_0 = 0.0f;
    (&_S2707)->differential_0 = 0.0f;
    DiffPair_float_0 _S2708;
    (&_S2708)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S2708)->differential_0 = 0.0f;
    _d_max_0(&_S2707, &_S2708, (*_s_dOut_16)[int(2)]);
    float _S2709 = dpparams_3->primal_0[int(12)] * (*_s_dOut_16)[int(1)];
    float _S2710 = dpparams_3->primal_0[int(11)] * (*_s_dOut_16)[int(1)];
    float _S2711 = dpparams_3->primal_0[int(7)] * (*_s_dOut_16)[int(1)];
    float _S2712 = dpparams_3->primal_0[int(6)] * (*_s_dOut_16)[int(1)];
    float _S2713 = dpparams_3->primal_0[int(2)] * (*_s_dOut_16)[int(1)];
    float _S2714 = dpparams_3->primal_0[int(1)] * (*_s_dOut_16)[int(1)];
    PPISPParamsRQS_0 _S2715 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2715)->color_params_0 = _S2650;
    (&_S2715)->exposure_0 = (*_s_dOut_16)[int(0)];
    _S2581 = _S2715;
    (&(&_S2581)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2716 = _S2599 + _S2604 + _S2715.crf_params_0[int(2)].gc_0;
    (&(&_S2581)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2717 = _S2607 + _S2612 + _S2715.crf_params_0[int(2)].y0_0;
    (&(&_S2581)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2718 = _S2615 + _S2620 + _S2715.crf_params_0[int(2)].x0_0;
    (&(&_S2581)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2719 = _S2623 + _S2628 + _S2715.crf_params_0[int(2)].g1_0;
    (&(&_S2581)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2720 = _S2631 + _S2636 + _S2715.crf_params_0[int(2)].g0_0;
    (&(&_S2581)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2721 = _S2601 + _S2604 + _S2715.crf_params_0[int(1)].gc_0;
    (&(&_S2581)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2722 = _S2609 + _S2612 + _S2715.crf_params_0[int(1)].y0_0;
    (&(&_S2581)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2723 = _S2617 + _S2620 + _S2715.crf_params_0[int(1)].x0_0;
    (&(&_S2581)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2724 = _S2625 + _S2628 + _S2715.crf_params_0[int(1)].g1_0;
    (&(&_S2581)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2725 = _S2633 + _S2636 + _S2715.crf_params_0[int(1)].g0_0;
    (&(&_S2581)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2726 = _S2603 + _S2604 + _S2715.crf_params_0[int(0)].gc_0;
    (&(&_S2581)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2727 = _S2611 + _S2612 + _S2715.crf_params_0[int(0)].y0_0;
    (&(&_S2581)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2728 = _S2619 + _S2620 + _S2715.crf_params_0[int(0)].x0_0;
    (&(&_S2581)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2729 = _S2627 + _S2628 + _S2715.crf_params_0[int(0)].g1_0;
    (&(&_S2581)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2730 = _S2635 + _S2636 + _S2715.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S2581)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2581)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2581)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2731 = _S2653 + _S2658 + _S2692.differential_0 + _S2715.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2581)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2732 = _S2661 + _S2666 + _S2698.differential_0 + _S2715.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2581)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2733 = _S2669 + _S2674 + _S2704.differential_0 + _S2715.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2581)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2734 = _S2677 + _S2682 + _S2709 + _S2709 + _S2715.vignette_params_0[int(2)].cy_0;
    (&(&_S2581)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2735 = _S2685 + _S2690 + _S2710 + _S2710 + _S2715.vignette_params_0[int(2)].cx_0;
    (&(&_S2581)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2736 = _S2655 + _S2658 + _S2694.differential_0 + _S2715.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2581)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2737 = _S2663 + _S2666 + _S2700.differential_0 + _S2715.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2581)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2738 = _S2671 + _S2674 + _S2706.differential_0 + _S2715.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2581)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2739 = _S2679 + _S2682 + _S2711 + _S2711 + _S2715.vignette_params_0[int(1)].cy_0;
    (&(&_S2581)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2740 = _S2687 + _S2690 + _S2712 + _S2712 + _S2715.vignette_params_0[int(1)].cx_0;
    (&(&_S2581)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2741 = _S2657 + _S2658 + _S2696.differential_0 + _S2715.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2581)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2742 = _S2665 + _S2666 + _S2702.differential_0 + _S2715.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2581)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2743 = _S2673 + _S2674 + _S2708.differential_0 + _S2715.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2581)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2744 = _S2681 + _S2682 + _S2713 + _S2713 + _S2715.vignette_params_0[int(0)].cy_0;
    (&(&_S2581)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2745 = _S2689 + _S2690 + _S2714 + _S2714 + _S2715.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2746;
    _S2746[int(0)] = 0.0f;
    _S2746[int(1)] = 0.0f;
    _S2746[int(2)] = 0.0f;
    _S2746[int(3)] = 0.0f;
    _S2746[int(4)] = 0.0f;
    _S2746[int(5)] = 0.0f;
    _S2746[int(6)] = 0.0f;
    _S2746[int(7)] = 0.0f;
    _S2746[int(8)] = 0.0f;
    _S2746[int(9)] = 0.0f;
    _S2746[int(10)] = 0.0f;
    _S2746[int(11)] = 0.0f;
    _S2746[int(12)] = 0.0f;
    _S2746[int(13)] = 0.0f;
    _S2746[int(14)] = 0.0f;
    _S2746[int(15)] = 0.0f;
    _S2746[int(16)] = 0.0f;
    _S2746[int(17)] = 0.0f;
    _S2746[int(18)] = 0.0f;
    _S2746[int(19)] = 0.0f;
    _S2746[int(20)] = 0.0f;
    _S2746[int(21)] = 0.0f;
    _S2746[int(22)] = 0.0f;
    _S2746[int(23)] = 0.0f;
    _S2746[int(24)] = 0.0f;
    _S2746[int(25)] = 0.0f;
    _S2746[int(26)] = 0.0f;
    _S2746[int(27)] = 0.0f;
    _S2746[int(28)] = 0.0f;
    _S2746[int(29)] = 0.0f;
    _S2746[int(30)] = 0.0f;
    _S2746[int(31)] = 0.0f;
    _S2746[int(32)] = 0.0f;
    _S2746[int(33)] = 0.0f;
    _S2746[int(34)] = 0.0f;
    _S2746[int(35)] = 0.0f;
    _S2746[int(36)] = 0.0f;
    _S2746[int(37)] = 0.0f;
    _S2746[int(38)] = 0.0f;
    _S2746[int(9)] = _S2737;
    _S2746[int(18)] = _S2715.color_params_0.r_0.x;
    _S2746[int(17)] = _S2715.color_params_0.b_0.y;
    _S2746[int(16)] = _S2715.color_params_0.b_0.x;
    _S2746[int(15)] = _S2731;
    _S2746[int(14)] = _S2732;
    _S2746[int(13)] = _S2733;
    _S2746[int(12)] = _S2734;
    _S2746[int(11)] = _S2735;
    _S2746[int(10)] = _S2736;
    _S2746[int(19)] = _S2715.color_params_0.r_0.y;
    _S2746[int(8)] = _S2738;
    _S2746[int(7)] = _S2739;
    _S2746[int(6)] = _S2740;
    _S2746[int(5)] = _S2741;
    _S2746[int(4)] = _S2742;
    _S2746[int(3)] = _S2743;
    _S2746[int(2)] = _S2744;
    _S2746[int(1)] = _S2745;
    _S2746[int(0)] = _S2581.exposure_0;
    _S2746[int(28)] = _S2726;
    _S2746[int(37)] = _S2717;
    _S2746[int(36)] = _S2718;
    _S2746[int(35)] = _S2719;
    _S2746[int(34)] = _S2720;
    _S2746[int(33)] = _S2721;
    _S2746[int(32)] = _S2722;
    _S2746[int(31)] = _S2723;
    _S2746[int(30)] = _S2724;
    _S2746[int(29)] = _S2725;
    _S2746[int(38)] = _S2716;
    _S2746[int(27)] = _S2727;
    _S2746[int(26)] = _S2728;
    _S2746[int(25)] = _S2729;
    _S2746[int(24)] = _S2730;
    _S2746[int(23)] = _S2715.color_params_0.n_0.y;
    _S2746[int(22)] = _S2715.color_params_0.n_0.x;
    _S2746[int(21)] = _S2715.color_params_0.g_0.y;
    _S2746[int(20)] = _S2715.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S2746;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2747, FixedArray<float, 23>  * _S2748)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S2747, _S2748);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_7, FixedArray<float, 23>  grad_out_3, FixedArray<float, 39>  * _S2749)
{
    FixedArray<float, 39>  _S2750 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_7;
    (&dp_params_3)->differential_0 = _S2750;
    FixedArray<float, 23>  _S2751 = grad_out_3;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, &_S2751);
    *_S2749 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_2, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S2752)
{
    float _S2753;
    FixedArray<float, 6>  losses_5;
    float _S2754 = float(num_cameras_0);
    float _S2755 = raw_losses_2[int(0)] / _S2754;
    for(;;)
    {
        float _S2756 = (F32_abs((_S2755)));
        if(_S2756 < 0.10000000149011612f)
        {
            _S2753 = 0.5f * _S2755 * _S2755 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2753 = _S2756 - 0.05000000074505806f;
            break;
        }
    }
    losses_5[int(0)] = _S2753;
    losses_5[int(1)] = raw_losses_2[int(1)] / (3.0f * _S2754);
    losses_5[int(2)] = (raw_losses_2[int(2)] + raw_losses_2[int(3)] + raw_losses_2[int(4)]) / (9.0f * _S2754);
    losses_5[int(3)] = (raw_losses_2[int(5)] + raw_losses_2[int(6)] + raw_losses_2[int(7)] + raw_losses_2[int(8)] + raw_losses_2[int(9)]) / (5.0f * _S2754);
    float _S2757 = raw_losses_2[int(10)] / _S2754;
    for(;;)
    {
        float _S2758 = (F32_abs((_S2757)));
        if(_S2758 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2757 * _S2757 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2758 - 0.00249999994412065f;
            break;
        }
    }
    float _S2759;
    float _S2760 = raw_losses_2[int(11)] / _S2754;
    for(;;)
    {
        float _S2761 = (F32_abs((_S2760)));
        if(_S2761 < 0.00499999988824129f)
        {
            _S2759 = 0.5f * _S2760 * _S2760 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2759 = _S2761 - 0.00249999994412065f;
            break;
        }
    }
    float _S2762 = _S2753 + _S2759;
    float _S2763 = raw_losses_2[int(12)] / _S2754;
    for(;;)
    {
        float _S2764 = (F32_abs((_S2763)));
        if(_S2764 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2763 * _S2763 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2764 - 0.00249999994412065f;
            break;
        }
    }
    float _S2765 = _S2762 + _S2753;
    float _S2766 = raw_losses_2[int(13)] / _S2754;
    for(;;)
    {
        float _S2767 = (F32_abs((_S2766)));
        if(_S2767 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2766 * _S2766 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2767 - 0.00249999994412065f;
            break;
        }
    }
    float _S2768 = _S2765 + _S2753;
    float _S2769 = raw_losses_2[int(14)] / _S2754;
    for(;;)
    {
        float _S2770 = (F32_abs((_S2769)));
        if(_S2770 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2769 * _S2769 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2770 - 0.00249999994412065f;
            break;
        }
    }
    float _S2771 = _S2768 + _S2753;
    float _S2772 = raw_losses_2[int(15)] / _S2754;
    for(;;)
    {
        float _S2773 = (F32_abs((_S2772)));
        if(_S2773 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2772 * _S2772 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2773 - 0.00249999994412065f;
            break;
        }
    }
    float _S2774 = _S2771 + _S2753;
    float _S2775 = raw_losses_2[int(16)] / _S2754;
    for(;;)
    {
        float _S2776 = (F32_abs((_S2775)));
        if(_S2776 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2775 * _S2775 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2776 - 0.00249999994412065f;
            break;
        }
    }
    float _S2777 = _S2774 + _S2753;
    float _S2778 = raw_losses_2[int(17)] / _S2754;
    for(;;)
    {
        float _S2779 = (F32_abs((_S2778)));
        if(_S2779 < 0.00499999988824129f)
        {
            _S2753 = 0.5f * _S2778 * _S2778 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2753 = _S2779 - 0.00249999994412065f;
            break;
        }
    }
    float _S2780 = (_S2777 + _S2753) / 8.0f;
    float _S2781 = (raw_losses_2[int(18)] + raw_losses_2[int(19)] + raw_losses_2[int(20)] + raw_losses_2[int(21)]) / (4.0f * _S2754);
    losses_5[int(0)] = losses_5[int(0)] * loss_weights_0[int(0)];
    losses_5[int(1)] = losses_5[int(1)] * loss_weights_0[int(1)];
    losses_5[int(2)] = losses_5[int(2)] * loss_weights_0[int(2)];
    losses_5[int(3)] = losses_5[int(3)] * loss_weights_0[int(3)];
    losses_5[int(4)] = _S2780 * loss_weights_0[int(4)];
    losses_5[int(5)] = _S2781 * loss_weights_0[int(5)];
    *_S2752 = losses_5;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_3, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S2782)
{
    float _S2783;
    FixedArray<float, 6>  losses_6;
    float _S2784 = float(num_cameras_1);
    float _S2785 = raw_losses_3[int(0)] / _S2784;
    for(;;)
    {
        float _S2786 = (F32_abs((_S2785)));
        if(_S2786 < 0.10000000149011612f)
        {
            _S2783 = 0.5f * _S2785 * _S2785 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2783 = _S2786 - 0.05000000074505806f;
            break;
        }
    }
    losses_6[int(0)] = _S2783;
    losses_6[int(1)] = raw_losses_3[int(1)] / (3.0f * _S2784);
    losses_6[int(2)] = (raw_losses_3[int(2)] + raw_losses_3[int(3)] + raw_losses_3[int(4)]) / (9.0f * _S2784);
    float _S2787 = 5.0f * _S2784;
    losses_6[int(3)] = (raw_losses_3[int(5)] + raw_losses_3[int(6)] + raw_losses_3[int(7)] + raw_losses_3[int(8)] + raw_losses_3[int(9)]) / _S2787;
    float _S2788 = raw_losses_3[int(10)] / _S2784;
    for(;;)
    {
        float _S2789 = (F32_abs((_S2788)));
        if(_S2789 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2788 * _S2788 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2789 - 0.00249999994412065f;
            break;
        }
    }
    float _S2790;
    float _S2791 = raw_losses_3[int(11)] / _S2784;
    for(;;)
    {
        float _S2792 = (F32_abs((_S2791)));
        if(_S2792 < 0.00499999988824129f)
        {
            _S2790 = 0.5f * _S2791 * _S2791 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2790 = _S2792 - 0.00249999994412065f;
            break;
        }
    }
    float _S2793 = _S2783 + _S2790;
    float _S2794 = raw_losses_3[int(12)] / _S2784;
    for(;;)
    {
        float _S2795 = (F32_abs((_S2794)));
        if(_S2795 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2794 * _S2794 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2795 - 0.00249999994412065f;
            break;
        }
    }
    float _S2796 = _S2793 + _S2783;
    float _S2797 = raw_losses_3[int(13)] / _S2784;
    for(;;)
    {
        float _S2798 = (F32_abs((_S2797)));
        if(_S2798 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2797 * _S2797 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2798 - 0.00249999994412065f;
            break;
        }
    }
    float _S2799 = _S2796 + _S2783;
    float _S2800 = raw_losses_3[int(14)] / _S2784;
    for(;;)
    {
        float _S2801 = (F32_abs((_S2800)));
        if(_S2801 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2800 * _S2800 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2801 - 0.00249999994412065f;
            break;
        }
    }
    float _S2802 = _S2799 + _S2783;
    float _S2803 = raw_losses_3[int(15)] / _S2784;
    for(;;)
    {
        float _S2804 = (F32_abs((_S2803)));
        if(_S2804 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2803 * _S2803 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2804 - 0.00249999994412065f;
            break;
        }
    }
    float _S2805 = _S2802 + _S2783;
    float _S2806 = raw_losses_3[int(16)] / _S2784;
    for(;;)
    {
        float _S2807 = (F32_abs((_S2806)));
        if(_S2807 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2806 * _S2806 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2807 - 0.00249999994412065f;
            break;
        }
    }
    float _S2808 = _S2805 + _S2783;
    float _S2809 = raw_losses_3[int(17)] / _S2784;
    for(;;)
    {
        float _S2810 = (F32_abs((_S2809)));
        if(_S2810 < 0.00499999988824129f)
        {
            _S2783 = 0.5f * _S2809 * _S2809 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2783 = _S2810 - 0.00249999994412065f;
            break;
        }
    }
    float _S2811 = (_S2808 + _S2783) / 8.0f;
    float _S2812 = (raw_losses_3[int(18)] + raw_losses_3[int(19)] + raw_losses_3[int(20)] + raw_losses_3[int(21)] + raw_losses_3[int(22)]) / _S2787;
    losses_6[int(0)] = losses_6[int(0)] * loss_weights_1[int(0)];
    losses_6[int(1)] = losses_6[int(1)] * loss_weights_1[int(1)];
    losses_6[int(2)] = losses_6[int(2)] * loss_weights_1[int(2)];
    losses_6[int(3)] = losses_6[int(3)] * loss_weights_1[int(3)];
    losses_6[int(4)] = _S2811 * loss_weights_1[int(4)];
    losses_6[int(5)] = _S2812 * loss_weights_1[int(5)];
    *_S2782 = losses_6;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_17)
{
    FixedArray<float, 22>  _S2813 = dpraw_losses_1->primal_0;
    float _S2814 = float(num_cameras_2);
    float _S2815 = dpraw_losses_1->primal_0[int(0)] / _S2814;
    bool _S2816 = (s_primal_ctx_abs_0(_S2815)) < 0.10000000149011612f;
    float _S2817;
    if(_S2816)
    {
        _S2817 = 0.5f * _S2815;
    }
    else
    {
        _S2817 = 0.0f;
    }
    float _S2818 = 3.0f * _S2814;
    float _S2819 = 9.0f * _S2814;
    float _S2820 = 5.0f * _S2814;
    float _S2821 = _S2813[int(10)] / _S2814;
    bool _S2822 = (s_primal_ctx_abs_0(_S2821)) < 0.00499999988824129f;
    float _S2823;
    if(_S2822)
    {
        _S2823 = 0.5f * _S2821;
    }
    else
    {
        _S2823 = 0.0f;
    }
    float _S2824 = _S2813[int(11)] / _S2814;
    bool _S2825 = (s_primal_ctx_abs_0(_S2824)) < 0.00499999988824129f;
    float _S2826;
    if(_S2825)
    {
        _S2826 = 0.5f * _S2824;
    }
    else
    {
        _S2826 = 0.0f;
    }
    float _S2827 = _S2813[int(12)] / _S2814;
    bool _S2828 = (s_primal_ctx_abs_0(_S2827)) < 0.00499999988824129f;
    float _S2829;
    if(_S2828)
    {
        _S2829 = 0.5f * _S2827;
    }
    else
    {
        _S2829 = 0.0f;
    }
    float _S2830 = _S2813[int(13)] / _S2814;
    bool _S2831 = (s_primal_ctx_abs_0(_S2830)) < 0.00499999988824129f;
    float _S2832;
    if(_S2831)
    {
        _S2832 = 0.5f * _S2830;
    }
    else
    {
        _S2832 = 0.0f;
    }
    float _S2833 = _S2813[int(14)] / _S2814;
    bool _S2834 = (s_primal_ctx_abs_0(_S2833)) < 0.00499999988824129f;
    float _S2835;
    if(_S2834)
    {
        _S2835 = 0.5f * _S2833;
    }
    else
    {
        _S2835 = 0.0f;
    }
    float _S2836 = _S2813[int(15)] / _S2814;
    bool _S2837 = (s_primal_ctx_abs_0(_S2836)) < 0.00499999988824129f;
    float _S2838;
    if(_S2837)
    {
        _S2838 = 0.5f * _S2836;
    }
    else
    {
        _S2838 = 0.0f;
    }
    float _S2839 = _S2813[int(16)] / _S2814;
    bool _S2840 = (s_primal_ctx_abs_0(_S2839)) < 0.00499999988824129f;
    float _S2841;
    if(_S2840)
    {
        _S2841 = 0.5f * _S2839;
    }
    else
    {
        _S2841 = 0.0f;
    }
    float _S2842 = _S2813[int(17)] / _S2814;
    bool _S2843 = (s_primal_ctx_abs_0(_S2842)) < 0.00499999988824129f;
    float _S2844;
    if(_S2843)
    {
        _S2844 = 0.5f * _S2842;
    }
    else
    {
        _S2844 = 0.0f;
    }
    float _S2845 = (*loss_weights_2)[int(3)] * (*_s_dOut_17)[int(3)];
    float _S2846 = (*loss_weights_2)[int(2)] * (*_s_dOut_17)[int(2)];
    float _S2847 = (*loss_weights_2)[int(1)] * (*_s_dOut_17)[int(1)];
    float _S2848 = (*loss_weights_2)[int(0)] * (*_s_dOut_17)[int(0)];
    float _S2849 = (*loss_weights_2)[int(5)] * (*_s_dOut_17)[int(5)] / (4.0f * _S2814);
    float _S2850 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_17)[int(4)]);
    FixedArray<float, 22>  _S2851;
    _S2851[int(0)] = 0.0f;
    _S2851[int(1)] = 0.0f;
    _S2851[int(2)] = 0.0f;
    _S2851[int(3)] = 0.0f;
    _S2851[int(4)] = 0.0f;
    _S2851[int(5)] = 0.0f;
    _S2851[int(6)] = 0.0f;
    _S2851[int(7)] = 0.0f;
    _S2851[int(8)] = 0.0f;
    _S2851[int(9)] = 0.0f;
    _S2851[int(10)] = 0.0f;
    _S2851[int(11)] = 0.0f;
    _S2851[int(12)] = 0.0f;
    _S2851[int(13)] = 0.0f;
    _S2851[int(14)] = 0.0f;
    _S2851[int(15)] = 0.0f;
    _S2851[int(16)] = 0.0f;
    _S2851[int(17)] = 0.0f;
    _S2851[int(18)] = 0.0f;
    _S2851[int(19)] = 0.0f;
    _S2851[int(20)] = 0.0f;
    _S2851[int(21)] = 0.0f;
    _S2851[int(21)] = _S2849;
    _S2851[int(20)] = _S2849;
    _S2851[int(19)] = _S2849;
    _S2851[int(18)] = _S2849;
    float _S2852 = _S2851[int(0)];
    float _S2853 = _S2851[int(1)];
    float _S2854 = _S2851[int(2)];
    float _S2855 = _S2851[int(3)];
    float _S2856 = _S2851[int(4)];
    float _S2857 = _S2851[int(5)];
    float _S2858 = _S2851[int(6)];
    float _S2859 = _S2851[int(7)];
    float _S2860 = _S2851[int(8)];
    float _S2861 = _S2851[int(9)];
    float _S2862 = _S2851[int(10)];
    float _S2863 = _S2851[int(11)];
    float _S2864 = _S2851[int(12)];
    float _S2865 = _S2851[int(13)];
    float _S2866 = _S2851[int(14)];
    float _S2867 = _S2851[int(15)];
    float _S2868 = _S2851[int(16)];
    float _S2869 = _S2851[int(17)];
    float _S2870 = _S2851[int(18)];
    float _S2871 = _S2851[int(19)];
    float _S2872 = _S2851[int(20)];
    float _S2873 = _S2851[int(21)];
    float _S2874;
    if(_S2843)
    {
        float _S2875 = 200.0f * _S2850;
        float _S2876 = _S2844 * _S2875 + 0.5f * (_S2842 * _S2875);
        _S2844 = 0.0f;
        _S2874 = _S2876;
    }
    else
    {
        _S2844 = _S2850;
        _S2874 = 0.0f;
    }
    DiffPair_float_0 _S2877;
    (&_S2877)->primal_0 = _S2842;
    (&_S2877)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2877, _S2844);
    float _S2878 = (_S2877.differential_0 + _S2874) / _S2814;
    FixedArray<float, 22>  _S2879;
    _S2879[int(0)] = 0.0f;
    _S2879[int(1)] = 0.0f;
    _S2879[int(2)] = 0.0f;
    _S2879[int(3)] = 0.0f;
    _S2879[int(4)] = 0.0f;
    _S2879[int(5)] = 0.0f;
    _S2879[int(6)] = 0.0f;
    _S2879[int(7)] = 0.0f;
    _S2879[int(8)] = 0.0f;
    _S2879[int(9)] = 0.0f;
    _S2879[int(10)] = 0.0f;
    _S2879[int(11)] = 0.0f;
    _S2879[int(12)] = 0.0f;
    _S2879[int(13)] = 0.0f;
    _S2879[int(14)] = 0.0f;
    _S2879[int(15)] = 0.0f;
    _S2879[int(16)] = 0.0f;
    _S2879[int(17)] = 0.0f;
    _S2879[int(18)] = 0.0f;
    _S2879[int(19)] = 0.0f;
    _S2879[int(20)] = 0.0f;
    _S2879[int(21)] = 0.0f;
    _S2879[int(17)] = _S2878;
    float _S2880 = _S2852 + _S2879[int(0)];
    float _S2881 = _S2853 + _S2879[int(1)];
    float _S2882 = _S2854 + _S2879[int(2)];
    float _S2883 = _S2855 + _S2879[int(3)];
    float _S2884 = _S2856 + _S2879[int(4)];
    float _S2885 = _S2857 + _S2879[int(5)];
    float _S2886 = _S2858 + _S2879[int(6)];
    float _S2887 = _S2859 + _S2879[int(7)];
    float _S2888 = _S2860 + _S2879[int(8)];
    float _S2889 = _S2861 + _S2879[int(9)];
    float _S2890 = _S2862 + _S2879[int(10)];
    float _S2891 = _S2863 + _S2879[int(11)];
    float _S2892 = _S2864 + _S2879[int(12)];
    float _S2893 = _S2865 + _S2879[int(13)];
    float _S2894 = _S2866 + _S2879[int(14)];
    float _S2895 = _S2867 + _S2879[int(15)];
    float _S2896 = _S2868 + _S2879[int(16)];
    float _S2897 = _S2869 + _S2879[int(17)];
    float _S2898 = _S2870 + _S2879[int(18)];
    float _S2899 = _S2871 + _S2879[int(19)];
    float _S2900 = _S2872 + _S2879[int(20)];
    float _S2901 = _S2873 + _S2879[int(21)];
    if(_S2840)
    {
        float _S2902 = 200.0f * _S2850;
        float _S2903 = _S2841 * _S2902 + 0.5f * (_S2839 * _S2902);
        _S2841 = 0.0f;
        _S2844 = _S2903;
    }
    else
    {
        _S2841 = _S2850;
        _S2844 = 0.0f;
    }
    DiffPair_float_0 _S2904;
    (&_S2904)->primal_0 = _S2839;
    (&_S2904)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2904, _S2841);
    float _S2905 = (_S2904.differential_0 + _S2844) / _S2814;
    FixedArray<float, 22>  _S2906;
    _S2906[int(0)] = 0.0f;
    _S2906[int(1)] = 0.0f;
    _S2906[int(2)] = 0.0f;
    _S2906[int(3)] = 0.0f;
    _S2906[int(4)] = 0.0f;
    _S2906[int(5)] = 0.0f;
    _S2906[int(6)] = 0.0f;
    _S2906[int(7)] = 0.0f;
    _S2906[int(8)] = 0.0f;
    _S2906[int(9)] = 0.0f;
    _S2906[int(10)] = 0.0f;
    _S2906[int(11)] = 0.0f;
    _S2906[int(12)] = 0.0f;
    _S2906[int(13)] = 0.0f;
    _S2906[int(14)] = 0.0f;
    _S2906[int(15)] = 0.0f;
    _S2906[int(16)] = 0.0f;
    _S2906[int(17)] = 0.0f;
    _S2906[int(18)] = 0.0f;
    _S2906[int(19)] = 0.0f;
    _S2906[int(20)] = 0.0f;
    _S2906[int(21)] = 0.0f;
    _S2906[int(16)] = _S2905;
    float _S2907 = _S2880 + _S2906[int(0)];
    float _S2908 = _S2881 + _S2906[int(1)];
    float _S2909 = _S2882 + _S2906[int(2)];
    float _S2910 = _S2883 + _S2906[int(3)];
    float _S2911 = _S2884 + _S2906[int(4)];
    float _S2912 = _S2885 + _S2906[int(5)];
    float _S2913 = _S2886 + _S2906[int(6)];
    float _S2914 = _S2887 + _S2906[int(7)];
    float _S2915 = _S2888 + _S2906[int(8)];
    float _S2916 = _S2889 + _S2906[int(9)];
    float _S2917 = _S2890 + _S2906[int(10)];
    float _S2918 = _S2891 + _S2906[int(11)];
    float _S2919 = _S2892 + _S2906[int(12)];
    float _S2920 = _S2893 + _S2906[int(13)];
    float _S2921 = _S2894 + _S2906[int(14)];
    float _S2922 = _S2895 + _S2906[int(15)];
    float _S2923 = _S2896 + _S2906[int(16)];
    float _S2924 = _S2897 + _S2906[int(17)];
    float _S2925 = _S2898 + _S2906[int(18)];
    float _S2926 = _S2899 + _S2906[int(19)];
    float _S2927 = _S2900 + _S2906[int(20)];
    float _S2928 = _S2901 + _S2906[int(21)];
    if(_S2837)
    {
        float _S2929 = 200.0f * _S2850;
        float _S2930 = _S2838 * _S2929 + 0.5f * (_S2836 * _S2929);
        _S2838 = 0.0f;
        _S2841 = _S2930;
    }
    else
    {
        _S2838 = _S2850;
        _S2841 = 0.0f;
    }
    DiffPair_float_0 _S2931;
    (&_S2931)->primal_0 = _S2836;
    (&_S2931)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2931, _S2838);
    float _S2932 = (_S2931.differential_0 + _S2841) / _S2814;
    FixedArray<float, 22>  _S2933;
    _S2933[int(0)] = 0.0f;
    _S2933[int(1)] = 0.0f;
    _S2933[int(2)] = 0.0f;
    _S2933[int(3)] = 0.0f;
    _S2933[int(4)] = 0.0f;
    _S2933[int(5)] = 0.0f;
    _S2933[int(6)] = 0.0f;
    _S2933[int(7)] = 0.0f;
    _S2933[int(8)] = 0.0f;
    _S2933[int(9)] = 0.0f;
    _S2933[int(10)] = 0.0f;
    _S2933[int(11)] = 0.0f;
    _S2933[int(12)] = 0.0f;
    _S2933[int(13)] = 0.0f;
    _S2933[int(14)] = 0.0f;
    _S2933[int(15)] = 0.0f;
    _S2933[int(16)] = 0.0f;
    _S2933[int(17)] = 0.0f;
    _S2933[int(18)] = 0.0f;
    _S2933[int(19)] = 0.0f;
    _S2933[int(20)] = 0.0f;
    _S2933[int(21)] = 0.0f;
    _S2933[int(15)] = _S2932;
    float _S2934 = _S2907 + _S2933[int(0)];
    float _S2935 = _S2908 + _S2933[int(1)];
    float _S2936 = _S2909 + _S2933[int(2)];
    float _S2937 = _S2910 + _S2933[int(3)];
    float _S2938 = _S2911 + _S2933[int(4)];
    float _S2939 = _S2912 + _S2933[int(5)];
    float _S2940 = _S2913 + _S2933[int(6)];
    float _S2941 = _S2914 + _S2933[int(7)];
    float _S2942 = _S2915 + _S2933[int(8)];
    float _S2943 = _S2916 + _S2933[int(9)];
    float _S2944 = _S2917 + _S2933[int(10)];
    float _S2945 = _S2918 + _S2933[int(11)];
    float _S2946 = _S2919 + _S2933[int(12)];
    float _S2947 = _S2920 + _S2933[int(13)];
    float _S2948 = _S2921 + _S2933[int(14)];
    float _S2949 = _S2922 + _S2933[int(15)];
    float _S2950 = _S2923 + _S2933[int(16)];
    float _S2951 = _S2924 + _S2933[int(17)];
    float _S2952 = _S2925 + _S2933[int(18)];
    float _S2953 = _S2926 + _S2933[int(19)];
    float _S2954 = _S2927 + _S2933[int(20)];
    float _S2955 = _S2928 + _S2933[int(21)];
    if(_S2834)
    {
        float _S2956 = 200.0f * _S2850;
        float _S2957 = _S2835 * _S2956 + 0.5f * (_S2833 * _S2956);
        _S2835 = 0.0f;
        _S2838 = _S2957;
    }
    else
    {
        _S2835 = _S2850;
        _S2838 = 0.0f;
    }
    DiffPair_float_0 _S2958;
    (&_S2958)->primal_0 = _S2833;
    (&_S2958)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2958, _S2835);
    float _S2959 = (_S2958.differential_0 + _S2838) / _S2814;
    FixedArray<float, 22>  _S2960;
    _S2960[int(0)] = 0.0f;
    _S2960[int(1)] = 0.0f;
    _S2960[int(2)] = 0.0f;
    _S2960[int(3)] = 0.0f;
    _S2960[int(4)] = 0.0f;
    _S2960[int(5)] = 0.0f;
    _S2960[int(6)] = 0.0f;
    _S2960[int(7)] = 0.0f;
    _S2960[int(8)] = 0.0f;
    _S2960[int(9)] = 0.0f;
    _S2960[int(10)] = 0.0f;
    _S2960[int(11)] = 0.0f;
    _S2960[int(12)] = 0.0f;
    _S2960[int(13)] = 0.0f;
    _S2960[int(14)] = 0.0f;
    _S2960[int(15)] = 0.0f;
    _S2960[int(16)] = 0.0f;
    _S2960[int(17)] = 0.0f;
    _S2960[int(18)] = 0.0f;
    _S2960[int(19)] = 0.0f;
    _S2960[int(20)] = 0.0f;
    _S2960[int(21)] = 0.0f;
    _S2960[int(14)] = _S2959;
    float _S2961 = _S2934 + _S2960[int(0)];
    float _S2962 = _S2935 + _S2960[int(1)];
    float _S2963 = _S2936 + _S2960[int(2)];
    float _S2964 = _S2937 + _S2960[int(3)];
    float _S2965 = _S2938 + _S2960[int(4)];
    float _S2966 = _S2939 + _S2960[int(5)];
    float _S2967 = _S2940 + _S2960[int(6)];
    float _S2968 = _S2941 + _S2960[int(7)];
    float _S2969 = _S2942 + _S2960[int(8)];
    float _S2970 = _S2943 + _S2960[int(9)];
    float _S2971 = _S2944 + _S2960[int(10)];
    float _S2972 = _S2945 + _S2960[int(11)];
    float _S2973 = _S2946 + _S2960[int(12)];
    float _S2974 = _S2947 + _S2960[int(13)];
    float _S2975 = _S2948 + _S2960[int(14)];
    float _S2976 = _S2949 + _S2960[int(15)];
    float _S2977 = _S2950 + _S2960[int(16)];
    float _S2978 = _S2951 + _S2960[int(17)];
    float _S2979 = _S2952 + _S2960[int(18)];
    float _S2980 = _S2953 + _S2960[int(19)];
    float _S2981 = _S2954 + _S2960[int(20)];
    float _S2982 = _S2955 + _S2960[int(21)];
    if(_S2831)
    {
        float _S2983 = 200.0f * _S2850;
        float _S2984 = _S2832 * _S2983 + 0.5f * (_S2830 * _S2983);
        _S2832 = 0.0f;
        _S2835 = _S2984;
    }
    else
    {
        _S2832 = _S2850;
        _S2835 = 0.0f;
    }
    DiffPair_float_0 _S2985;
    (&_S2985)->primal_0 = _S2830;
    (&_S2985)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2985, _S2832);
    float _S2986 = (_S2985.differential_0 + _S2835) / _S2814;
    FixedArray<float, 22>  _S2987;
    _S2987[int(0)] = 0.0f;
    _S2987[int(1)] = 0.0f;
    _S2987[int(2)] = 0.0f;
    _S2987[int(3)] = 0.0f;
    _S2987[int(4)] = 0.0f;
    _S2987[int(5)] = 0.0f;
    _S2987[int(6)] = 0.0f;
    _S2987[int(7)] = 0.0f;
    _S2987[int(8)] = 0.0f;
    _S2987[int(9)] = 0.0f;
    _S2987[int(10)] = 0.0f;
    _S2987[int(11)] = 0.0f;
    _S2987[int(12)] = 0.0f;
    _S2987[int(13)] = 0.0f;
    _S2987[int(14)] = 0.0f;
    _S2987[int(15)] = 0.0f;
    _S2987[int(16)] = 0.0f;
    _S2987[int(17)] = 0.0f;
    _S2987[int(18)] = 0.0f;
    _S2987[int(19)] = 0.0f;
    _S2987[int(20)] = 0.0f;
    _S2987[int(21)] = 0.0f;
    _S2987[int(13)] = _S2986;
    float _S2988 = _S2961 + _S2987[int(0)];
    float _S2989 = _S2962 + _S2987[int(1)];
    float _S2990 = _S2963 + _S2987[int(2)];
    float _S2991 = _S2964 + _S2987[int(3)];
    float _S2992 = _S2965 + _S2987[int(4)];
    float _S2993 = _S2966 + _S2987[int(5)];
    float _S2994 = _S2967 + _S2987[int(6)];
    float _S2995 = _S2968 + _S2987[int(7)];
    float _S2996 = _S2969 + _S2987[int(8)];
    float _S2997 = _S2970 + _S2987[int(9)];
    float _S2998 = _S2971 + _S2987[int(10)];
    float _S2999 = _S2972 + _S2987[int(11)];
    float _S3000 = _S2973 + _S2987[int(12)];
    float _S3001 = _S2974 + _S2987[int(13)];
    float _S3002 = _S2975 + _S2987[int(14)];
    float _S3003 = _S2976 + _S2987[int(15)];
    float _S3004 = _S2977 + _S2987[int(16)];
    float _S3005 = _S2978 + _S2987[int(17)];
    float _S3006 = _S2979 + _S2987[int(18)];
    float _S3007 = _S2980 + _S2987[int(19)];
    float _S3008 = _S2981 + _S2987[int(20)];
    float _S3009 = _S2982 + _S2987[int(21)];
    if(_S2828)
    {
        float _S3010 = 200.0f * _S2850;
        float _S3011 = _S2829 * _S3010 + 0.5f * (_S2827 * _S3010);
        _S2829 = 0.0f;
        _S2832 = _S3011;
    }
    else
    {
        _S2829 = _S2850;
        _S2832 = 0.0f;
    }
    DiffPair_float_0 _S3012;
    (&_S3012)->primal_0 = _S2827;
    (&_S3012)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3012, _S2829);
    float _S3013 = (_S3012.differential_0 + _S2832) / _S2814;
    FixedArray<float, 22>  _S3014;
    _S3014[int(0)] = 0.0f;
    _S3014[int(1)] = 0.0f;
    _S3014[int(2)] = 0.0f;
    _S3014[int(3)] = 0.0f;
    _S3014[int(4)] = 0.0f;
    _S3014[int(5)] = 0.0f;
    _S3014[int(6)] = 0.0f;
    _S3014[int(7)] = 0.0f;
    _S3014[int(8)] = 0.0f;
    _S3014[int(9)] = 0.0f;
    _S3014[int(10)] = 0.0f;
    _S3014[int(11)] = 0.0f;
    _S3014[int(12)] = 0.0f;
    _S3014[int(13)] = 0.0f;
    _S3014[int(14)] = 0.0f;
    _S3014[int(15)] = 0.0f;
    _S3014[int(16)] = 0.0f;
    _S3014[int(17)] = 0.0f;
    _S3014[int(18)] = 0.0f;
    _S3014[int(19)] = 0.0f;
    _S3014[int(20)] = 0.0f;
    _S3014[int(21)] = 0.0f;
    _S3014[int(12)] = _S3013;
    float _S3015 = _S2988 + _S3014[int(0)];
    float _S3016 = _S2989 + _S3014[int(1)];
    float _S3017 = _S2990 + _S3014[int(2)];
    float _S3018 = _S2991 + _S3014[int(3)];
    float _S3019 = _S2992 + _S3014[int(4)];
    float _S3020 = _S2993 + _S3014[int(5)];
    float _S3021 = _S2994 + _S3014[int(6)];
    float _S3022 = _S2995 + _S3014[int(7)];
    float _S3023 = _S2996 + _S3014[int(8)];
    float _S3024 = _S2997 + _S3014[int(9)];
    float _S3025 = _S2998 + _S3014[int(10)];
    float _S3026 = _S2999 + _S3014[int(11)];
    float _S3027 = _S3000 + _S3014[int(12)];
    float _S3028 = _S3001 + _S3014[int(13)];
    float _S3029 = _S3002 + _S3014[int(14)];
    float _S3030 = _S3003 + _S3014[int(15)];
    float _S3031 = _S3004 + _S3014[int(16)];
    float _S3032 = _S3005 + _S3014[int(17)];
    float _S3033 = _S3006 + _S3014[int(18)];
    float _S3034 = _S3007 + _S3014[int(19)];
    float _S3035 = _S3008 + _S3014[int(20)];
    float _S3036 = _S3009 + _S3014[int(21)];
    if(_S2825)
    {
        float _S3037 = 200.0f * _S2850;
        float _S3038 = _S2826 * _S3037 + 0.5f * (_S2824 * _S3037);
        _S2826 = 0.0f;
        _S2829 = _S3038;
    }
    else
    {
        _S2826 = _S2850;
        _S2829 = 0.0f;
    }
    DiffPair_float_0 _S3039;
    (&_S3039)->primal_0 = _S2824;
    (&_S3039)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3039, _S2826);
    float _S3040 = (_S3039.differential_0 + _S2829) / _S2814;
    FixedArray<float, 22>  _S3041;
    _S3041[int(0)] = 0.0f;
    _S3041[int(1)] = 0.0f;
    _S3041[int(2)] = 0.0f;
    _S3041[int(3)] = 0.0f;
    _S3041[int(4)] = 0.0f;
    _S3041[int(5)] = 0.0f;
    _S3041[int(6)] = 0.0f;
    _S3041[int(7)] = 0.0f;
    _S3041[int(8)] = 0.0f;
    _S3041[int(9)] = 0.0f;
    _S3041[int(10)] = 0.0f;
    _S3041[int(11)] = 0.0f;
    _S3041[int(12)] = 0.0f;
    _S3041[int(13)] = 0.0f;
    _S3041[int(14)] = 0.0f;
    _S3041[int(15)] = 0.0f;
    _S3041[int(16)] = 0.0f;
    _S3041[int(17)] = 0.0f;
    _S3041[int(18)] = 0.0f;
    _S3041[int(19)] = 0.0f;
    _S3041[int(20)] = 0.0f;
    _S3041[int(21)] = 0.0f;
    _S3041[int(11)] = _S3040;
    float _S3042 = _S3015 + _S3041[int(0)];
    float _S3043 = _S3016 + _S3041[int(1)];
    float _S3044 = _S3017 + _S3041[int(2)];
    float _S3045 = _S3018 + _S3041[int(3)];
    float _S3046 = _S3019 + _S3041[int(4)];
    float _S3047 = _S3020 + _S3041[int(5)];
    float _S3048 = _S3021 + _S3041[int(6)];
    float _S3049 = _S3022 + _S3041[int(7)];
    float _S3050 = _S3023 + _S3041[int(8)];
    float _S3051 = _S3024 + _S3041[int(9)];
    float _S3052 = _S3025 + _S3041[int(10)];
    float _S3053 = _S3026 + _S3041[int(11)];
    float _S3054 = _S3027 + _S3041[int(12)];
    float _S3055 = _S3028 + _S3041[int(13)];
    float _S3056 = _S3029 + _S3041[int(14)];
    float _S3057 = _S3030 + _S3041[int(15)];
    float _S3058 = _S3031 + _S3041[int(16)];
    float _S3059 = _S3032 + _S3041[int(17)];
    float _S3060 = _S3033 + _S3041[int(18)];
    float _S3061 = _S3034 + _S3041[int(19)];
    float _S3062 = _S3035 + _S3041[int(20)];
    float _S3063 = _S3036 + _S3041[int(21)];
    if(_S2822)
    {
        float _S3064 = 200.0f * _S2850;
        float _S3065 = _S2823 * _S3064 + 0.5f * (_S2821 * _S3064);
        _S2823 = 0.0f;
        _S2826 = _S3065;
    }
    else
    {
        _S2823 = _S2850;
        _S2826 = 0.0f;
    }
    DiffPair_float_0 _S3066;
    (&_S3066)->primal_0 = _S2821;
    (&_S3066)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3066, _S2823);
    float _S3067 = (_S3066.differential_0 + _S2826) / _S2814;
    float _S3068 = _S2845 / _S2820;
    float _S3069 = _S2846 / _S2819;
    float _S3070 = _S2847 / _S2818;
    FixedArray<float, 22>  _S3071;
    _S3071[int(0)] = 0.0f;
    _S3071[int(1)] = 0.0f;
    _S3071[int(2)] = 0.0f;
    _S3071[int(3)] = 0.0f;
    _S3071[int(4)] = 0.0f;
    _S3071[int(5)] = 0.0f;
    _S3071[int(6)] = 0.0f;
    _S3071[int(7)] = 0.0f;
    _S3071[int(8)] = 0.0f;
    _S3071[int(9)] = 0.0f;
    _S3071[int(10)] = 0.0f;
    _S3071[int(11)] = 0.0f;
    _S3071[int(12)] = 0.0f;
    _S3071[int(13)] = 0.0f;
    _S3071[int(14)] = 0.0f;
    _S3071[int(15)] = 0.0f;
    _S3071[int(16)] = 0.0f;
    _S3071[int(17)] = 0.0f;
    _S3071[int(18)] = 0.0f;
    _S3071[int(19)] = 0.0f;
    _S3071[int(20)] = 0.0f;
    _S3071[int(21)] = 0.0f;
    _S3071[int(10)] = _S3067;
    _S3071[int(9)] = _S3068;
    _S3071[int(8)] = _S3068;
    _S3071[int(7)] = _S3068;
    _S3071[int(6)] = _S3068;
    _S3071[int(5)] = _S3068;
    _S3071[int(4)] = _S3069;
    _S3071[int(3)] = _S3069;
    _S3071[int(2)] = _S3069;
    _S3071[int(1)] = _S3070;
    float _S3072 = _S3042 + _S3071[int(0)];
    float _S3073 = _S3043 + _S3071[int(1)];
    float _S3074 = _S3044 + _S3071[int(2)];
    float _S3075 = _S3045 + _S3071[int(3)];
    float _S3076 = _S3046 + _S3071[int(4)];
    float _S3077 = _S3047 + _S3071[int(5)];
    float _S3078 = _S3048 + _S3071[int(6)];
    float _S3079 = _S3049 + _S3071[int(7)];
    float _S3080 = _S3050 + _S3071[int(8)];
    float _S3081 = _S3051 + _S3071[int(9)];
    float _S3082 = _S3052 + _S3071[int(10)];
    float _S3083 = _S3053 + _S3071[int(11)];
    float _S3084 = _S3054 + _S3071[int(12)];
    float _S3085 = _S3055 + _S3071[int(13)];
    float _S3086 = _S3056 + _S3071[int(14)];
    float _S3087 = _S3057 + _S3071[int(15)];
    float _S3088 = _S3058 + _S3071[int(16)];
    float _S3089 = _S3059 + _S3071[int(17)];
    float _S3090 = _S3060 + _S3071[int(18)];
    float _S3091 = _S3061 + _S3071[int(19)];
    float _S3092 = _S3062 + _S3071[int(20)];
    float _S3093 = _S3063 + _S3071[int(21)];
    if(_S2816)
    {
        float _S3094 = 10.0f * _S2848;
        float _S3095 = _S2817 * _S3094 + 0.5f * (_S2815 * _S3094);
        _S2817 = 0.0f;
        _S2823 = _S3095;
    }
    else
    {
        _S2817 = _S2848;
        _S2823 = 0.0f;
    }
    DiffPair_float_0 _S3096;
    (&_S3096)->primal_0 = _S2815;
    (&_S3096)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3096, _S2817);
    float _S3097 = (_S3096.differential_0 + _S2823) / _S2814;
    FixedArray<float, 22>  _S3098;
    _S3098[int(0)] = 0.0f;
    _S3098[int(1)] = 0.0f;
    _S3098[int(2)] = 0.0f;
    _S3098[int(3)] = 0.0f;
    _S3098[int(4)] = 0.0f;
    _S3098[int(5)] = 0.0f;
    _S3098[int(6)] = 0.0f;
    _S3098[int(7)] = 0.0f;
    _S3098[int(8)] = 0.0f;
    _S3098[int(9)] = 0.0f;
    _S3098[int(10)] = 0.0f;
    _S3098[int(11)] = 0.0f;
    _S3098[int(12)] = 0.0f;
    _S3098[int(13)] = 0.0f;
    _S3098[int(14)] = 0.0f;
    _S3098[int(15)] = 0.0f;
    _S3098[int(16)] = 0.0f;
    _S3098[int(17)] = 0.0f;
    _S3098[int(18)] = 0.0f;
    _S3098[int(19)] = 0.0f;
    _S3098[int(20)] = 0.0f;
    _S3098[int(21)] = 0.0f;
    _S3098[int(0)] = _S3097;
    FixedArray<float, 22>  _S3099 = {
        _S3072 + _S3098[int(0)], _S3073 + _S3098[int(1)], _S3074 + _S3098[int(2)], _S3075 + _S3098[int(3)], _S3076 + _S3098[int(4)], _S3077 + _S3098[int(5)], _S3078 + _S3098[int(6)], _S3079 + _S3098[int(7)], _S3080 + _S3098[int(8)], _S3081 + _S3098[int(9)], _S3082 + _S3098[int(10)], _S3083 + _S3098[int(11)], _S3084 + _S3098[int(12)], _S3085 + _S3098[int(13)], _S3086 + _S3098[int(14)], _S3087 + _S3098[int(15)], _S3088 + _S3098[int(16)], _S3089 + _S3098[int(17)], _S3090 + _S3098[int(18)], _S3091 + _S3098[int(19)], _S3092 + _S3098[int(20)], _S3093 + _S3098[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S3099;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S3100, int _S3101, FixedArray<float, 6>  * _S3102, FixedArray<float, 6>  * _S3103)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S3100, _S3101, _S3102, _S3103);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_4, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_4, FixedArray<float, 22>  * _S3104)
{
    FixedArray<float, 22>  _S3105 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_4;
    (&dp_raw_losses_1)->differential_0 = _S3105;
    FixedArray<float, 6>  _S3106 = loss_weights_3;
    FixedArray<float, 6>  _S3107 = grad_out_4;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_3, &_S3106, &_S3107);
    *_S3104 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_2, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_18)
{
    FixedArray<float, 23>  _S3108 = dpraw_losses_2->primal_0;
    float _S3109 = float(num_cameras_4);
    float _S3110 = dpraw_losses_2->primal_0[int(0)] / _S3109;
    bool _S3111 = (s_primal_ctx_abs_0(_S3110)) < 0.10000000149011612f;
    float _S3112;
    if(_S3111)
    {
        _S3112 = 0.5f * _S3110;
    }
    else
    {
        _S3112 = 0.0f;
    }
    float _S3113 = 3.0f * _S3109;
    float _S3114 = 9.0f * _S3109;
    float _S3115 = 5.0f * _S3109;
    float _S3116 = _S3108[int(10)] / _S3109;
    bool _S3117 = (s_primal_ctx_abs_0(_S3116)) < 0.00499999988824129f;
    float _S3118;
    if(_S3117)
    {
        _S3118 = 0.5f * _S3116;
    }
    else
    {
        _S3118 = 0.0f;
    }
    float _S3119 = _S3108[int(11)] / _S3109;
    bool _S3120 = (s_primal_ctx_abs_0(_S3119)) < 0.00499999988824129f;
    float _S3121;
    if(_S3120)
    {
        _S3121 = 0.5f * _S3119;
    }
    else
    {
        _S3121 = 0.0f;
    }
    float _S3122 = _S3108[int(12)] / _S3109;
    bool _S3123 = (s_primal_ctx_abs_0(_S3122)) < 0.00499999988824129f;
    float _S3124;
    if(_S3123)
    {
        _S3124 = 0.5f * _S3122;
    }
    else
    {
        _S3124 = 0.0f;
    }
    float _S3125 = _S3108[int(13)] / _S3109;
    bool _S3126 = (s_primal_ctx_abs_0(_S3125)) < 0.00499999988824129f;
    float _S3127;
    if(_S3126)
    {
        _S3127 = 0.5f * _S3125;
    }
    else
    {
        _S3127 = 0.0f;
    }
    float _S3128 = _S3108[int(14)] / _S3109;
    bool _S3129 = (s_primal_ctx_abs_0(_S3128)) < 0.00499999988824129f;
    float _S3130;
    if(_S3129)
    {
        _S3130 = 0.5f * _S3128;
    }
    else
    {
        _S3130 = 0.0f;
    }
    float _S3131 = _S3108[int(15)] / _S3109;
    bool _S3132 = (s_primal_ctx_abs_0(_S3131)) < 0.00499999988824129f;
    float _S3133;
    if(_S3132)
    {
        _S3133 = 0.5f * _S3131;
    }
    else
    {
        _S3133 = 0.0f;
    }
    float _S3134 = _S3108[int(16)] / _S3109;
    bool _S3135 = (s_primal_ctx_abs_0(_S3134)) < 0.00499999988824129f;
    float _S3136;
    if(_S3135)
    {
        _S3136 = 0.5f * _S3134;
    }
    else
    {
        _S3136 = 0.0f;
    }
    float _S3137 = _S3108[int(17)] / _S3109;
    bool _S3138 = (s_primal_ctx_abs_0(_S3137)) < 0.00499999988824129f;
    float _S3139;
    if(_S3138)
    {
        _S3139 = 0.5f * _S3137;
    }
    else
    {
        _S3139 = 0.0f;
    }
    float _S3140 = (*loss_weights_4)[int(3)] * (*_s_dOut_18)[int(3)];
    float _S3141 = (*loss_weights_4)[int(2)] * (*_s_dOut_18)[int(2)];
    float _S3142 = (*loss_weights_4)[int(1)] * (*_s_dOut_18)[int(1)];
    float _S3143 = (*loss_weights_4)[int(0)] * (*_s_dOut_18)[int(0)];
    float _S3144 = (*loss_weights_4)[int(5)] * (*_s_dOut_18)[int(5)] / _S3115;
    float _S3145 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_18)[int(4)]);
    FixedArray<float, 23>  _S3146;
    _S3146[int(0)] = 0.0f;
    _S3146[int(1)] = 0.0f;
    _S3146[int(2)] = 0.0f;
    _S3146[int(3)] = 0.0f;
    _S3146[int(4)] = 0.0f;
    _S3146[int(5)] = 0.0f;
    _S3146[int(6)] = 0.0f;
    _S3146[int(7)] = 0.0f;
    _S3146[int(8)] = 0.0f;
    _S3146[int(9)] = 0.0f;
    _S3146[int(10)] = 0.0f;
    _S3146[int(11)] = 0.0f;
    _S3146[int(12)] = 0.0f;
    _S3146[int(13)] = 0.0f;
    _S3146[int(14)] = 0.0f;
    _S3146[int(15)] = 0.0f;
    _S3146[int(16)] = 0.0f;
    _S3146[int(17)] = 0.0f;
    _S3146[int(18)] = 0.0f;
    _S3146[int(19)] = 0.0f;
    _S3146[int(20)] = 0.0f;
    _S3146[int(21)] = 0.0f;
    _S3146[int(22)] = 0.0f;
    _S3146[int(22)] = _S3144;
    _S3146[int(21)] = _S3144;
    _S3146[int(20)] = _S3144;
    _S3146[int(19)] = _S3144;
    _S3146[int(18)] = _S3144;
    float _S3147 = _S3146[int(0)];
    float _S3148 = _S3146[int(1)];
    float _S3149 = _S3146[int(2)];
    float _S3150 = _S3146[int(3)];
    float _S3151 = _S3146[int(4)];
    float _S3152 = _S3146[int(5)];
    float _S3153 = _S3146[int(6)];
    float _S3154 = _S3146[int(7)];
    float _S3155 = _S3146[int(8)];
    float _S3156 = _S3146[int(9)];
    float _S3157 = _S3146[int(10)];
    float _S3158 = _S3146[int(11)];
    float _S3159 = _S3146[int(12)];
    float _S3160 = _S3146[int(13)];
    float _S3161 = _S3146[int(14)];
    float _S3162 = _S3146[int(15)];
    float _S3163 = _S3146[int(16)];
    float _S3164 = _S3146[int(17)];
    float _S3165 = _S3146[int(18)];
    float _S3166 = _S3146[int(19)];
    float _S3167 = _S3146[int(20)];
    float _S3168 = _S3146[int(21)];
    float _S3169 = _S3146[int(22)];
    float _S3170;
    if(_S3138)
    {
        float _S3171 = 200.0f * _S3145;
        float _S3172 = _S3139 * _S3171 + 0.5f * (_S3137 * _S3171);
        _S3139 = 0.0f;
        _S3170 = _S3172;
    }
    else
    {
        _S3139 = _S3145;
        _S3170 = 0.0f;
    }
    DiffPair_float_0 _S3173;
    (&_S3173)->primal_0 = _S3137;
    (&_S3173)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3173, _S3139);
    float _S3174 = (_S3173.differential_0 + _S3170) / _S3109;
    FixedArray<float, 23>  _S3175;
    _S3175[int(0)] = 0.0f;
    _S3175[int(1)] = 0.0f;
    _S3175[int(2)] = 0.0f;
    _S3175[int(3)] = 0.0f;
    _S3175[int(4)] = 0.0f;
    _S3175[int(5)] = 0.0f;
    _S3175[int(6)] = 0.0f;
    _S3175[int(7)] = 0.0f;
    _S3175[int(8)] = 0.0f;
    _S3175[int(9)] = 0.0f;
    _S3175[int(10)] = 0.0f;
    _S3175[int(11)] = 0.0f;
    _S3175[int(12)] = 0.0f;
    _S3175[int(13)] = 0.0f;
    _S3175[int(14)] = 0.0f;
    _S3175[int(15)] = 0.0f;
    _S3175[int(16)] = 0.0f;
    _S3175[int(17)] = 0.0f;
    _S3175[int(18)] = 0.0f;
    _S3175[int(19)] = 0.0f;
    _S3175[int(20)] = 0.0f;
    _S3175[int(21)] = 0.0f;
    _S3175[int(22)] = 0.0f;
    _S3175[int(17)] = _S3174;
    float _S3176 = _S3147 + _S3175[int(0)];
    float _S3177 = _S3148 + _S3175[int(1)];
    float _S3178 = _S3149 + _S3175[int(2)];
    float _S3179 = _S3150 + _S3175[int(3)];
    float _S3180 = _S3151 + _S3175[int(4)];
    float _S3181 = _S3152 + _S3175[int(5)];
    float _S3182 = _S3153 + _S3175[int(6)];
    float _S3183 = _S3154 + _S3175[int(7)];
    float _S3184 = _S3155 + _S3175[int(8)];
    float _S3185 = _S3156 + _S3175[int(9)];
    float _S3186 = _S3157 + _S3175[int(10)];
    float _S3187 = _S3158 + _S3175[int(11)];
    float _S3188 = _S3159 + _S3175[int(12)];
    float _S3189 = _S3160 + _S3175[int(13)];
    float _S3190 = _S3161 + _S3175[int(14)];
    float _S3191 = _S3162 + _S3175[int(15)];
    float _S3192 = _S3163 + _S3175[int(16)];
    float _S3193 = _S3164 + _S3175[int(17)];
    float _S3194 = _S3165 + _S3175[int(18)];
    float _S3195 = _S3166 + _S3175[int(19)];
    float _S3196 = _S3167 + _S3175[int(20)];
    float _S3197 = _S3168 + _S3175[int(21)];
    float _S3198 = _S3169 + _S3175[int(22)];
    if(_S3135)
    {
        float _S3199 = 200.0f * _S3145;
        float _S3200 = _S3136 * _S3199 + 0.5f * (_S3134 * _S3199);
        _S3136 = 0.0f;
        _S3139 = _S3200;
    }
    else
    {
        _S3136 = _S3145;
        _S3139 = 0.0f;
    }
    DiffPair_float_0 _S3201;
    (&_S3201)->primal_0 = _S3134;
    (&_S3201)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3201, _S3136);
    float _S3202 = (_S3201.differential_0 + _S3139) / _S3109;
    FixedArray<float, 23>  _S3203;
    _S3203[int(0)] = 0.0f;
    _S3203[int(1)] = 0.0f;
    _S3203[int(2)] = 0.0f;
    _S3203[int(3)] = 0.0f;
    _S3203[int(4)] = 0.0f;
    _S3203[int(5)] = 0.0f;
    _S3203[int(6)] = 0.0f;
    _S3203[int(7)] = 0.0f;
    _S3203[int(8)] = 0.0f;
    _S3203[int(9)] = 0.0f;
    _S3203[int(10)] = 0.0f;
    _S3203[int(11)] = 0.0f;
    _S3203[int(12)] = 0.0f;
    _S3203[int(13)] = 0.0f;
    _S3203[int(14)] = 0.0f;
    _S3203[int(15)] = 0.0f;
    _S3203[int(16)] = 0.0f;
    _S3203[int(17)] = 0.0f;
    _S3203[int(18)] = 0.0f;
    _S3203[int(19)] = 0.0f;
    _S3203[int(20)] = 0.0f;
    _S3203[int(21)] = 0.0f;
    _S3203[int(22)] = 0.0f;
    _S3203[int(16)] = _S3202;
    float _S3204 = _S3176 + _S3203[int(0)];
    float _S3205 = _S3177 + _S3203[int(1)];
    float _S3206 = _S3178 + _S3203[int(2)];
    float _S3207 = _S3179 + _S3203[int(3)];
    float _S3208 = _S3180 + _S3203[int(4)];
    float _S3209 = _S3181 + _S3203[int(5)];
    float _S3210 = _S3182 + _S3203[int(6)];
    float _S3211 = _S3183 + _S3203[int(7)];
    float _S3212 = _S3184 + _S3203[int(8)];
    float _S3213 = _S3185 + _S3203[int(9)];
    float _S3214 = _S3186 + _S3203[int(10)];
    float _S3215 = _S3187 + _S3203[int(11)];
    float _S3216 = _S3188 + _S3203[int(12)];
    float _S3217 = _S3189 + _S3203[int(13)];
    float _S3218 = _S3190 + _S3203[int(14)];
    float _S3219 = _S3191 + _S3203[int(15)];
    float _S3220 = _S3192 + _S3203[int(16)];
    float _S3221 = _S3193 + _S3203[int(17)];
    float _S3222 = _S3194 + _S3203[int(18)];
    float _S3223 = _S3195 + _S3203[int(19)];
    float _S3224 = _S3196 + _S3203[int(20)];
    float _S3225 = _S3197 + _S3203[int(21)];
    float _S3226 = _S3198 + _S3203[int(22)];
    if(_S3132)
    {
        float _S3227 = 200.0f * _S3145;
        float _S3228 = _S3133 * _S3227 + 0.5f * (_S3131 * _S3227);
        _S3133 = 0.0f;
        _S3136 = _S3228;
    }
    else
    {
        _S3133 = _S3145;
        _S3136 = 0.0f;
    }
    DiffPair_float_0 _S3229;
    (&_S3229)->primal_0 = _S3131;
    (&_S3229)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3229, _S3133);
    float _S3230 = (_S3229.differential_0 + _S3136) / _S3109;
    FixedArray<float, 23>  _S3231;
    _S3231[int(0)] = 0.0f;
    _S3231[int(1)] = 0.0f;
    _S3231[int(2)] = 0.0f;
    _S3231[int(3)] = 0.0f;
    _S3231[int(4)] = 0.0f;
    _S3231[int(5)] = 0.0f;
    _S3231[int(6)] = 0.0f;
    _S3231[int(7)] = 0.0f;
    _S3231[int(8)] = 0.0f;
    _S3231[int(9)] = 0.0f;
    _S3231[int(10)] = 0.0f;
    _S3231[int(11)] = 0.0f;
    _S3231[int(12)] = 0.0f;
    _S3231[int(13)] = 0.0f;
    _S3231[int(14)] = 0.0f;
    _S3231[int(15)] = 0.0f;
    _S3231[int(16)] = 0.0f;
    _S3231[int(17)] = 0.0f;
    _S3231[int(18)] = 0.0f;
    _S3231[int(19)] = 0.0f;
    _S3231[int(20)] = 0.0f;
    _S3231[int(21)] = 0.0f;
    _S3231[int(22)] = 0.0f;
    _S3231[int(15)] = _S3230;
    float _S3232 = _S3204 + _S3231[int(0)];
    float _S3233 = _S3205 + _S3231[int(1)];
    float _S3234 = _S3206 + _S3231[int(2)];
    float _S3235 = _S3207 + _S3231[int(3)];
    float _S3236 = _S3208 + _S3231[int(4)];
    float _S3237 = _S3209 + _S3231[int(5)];
    float _S3238 = _S3210 + _S3231[int(6)];
    float _S3239 = _S3211 + _S3231[int(7)];
    float _S3240 = _S3212 + _S3231[int(8)];
    float _S3241 = _S3213 + _S3231[int(9)];
    float _S3242 = _S3214 + _S3231[int(10)];
    float _S3243 = _S3215 + _S3231[int(11)];
    float _S3244 = _S3216 + _S3231[int(12)];
    float _S3245 = _S3217 + _S3231[int(13)];
    float _S3246 = _S3218 + _S3231[int(14)];
    float _S3247 = _S3219 + _S3231[int(15)];
    float _S3248 = _S3220 + _S3231[int(16)];
    float _S3249 = _S3221 + _S3231[int(17)];
    float _S3250 = _S3222 + _S3231[int(18)];
    float _S3251 = _S3223 + _S3231[int(19)];
    float _S3252 = _S3224 + _S3231[int(20)];
    float _S3253 = _S3225 + _S3231[int(21)];
    float _S3254 = _S3226 + _S3231[int(22)];
    if(_S3129)
    {
        float _S3255 = 200.0f * _S3145;
        float _S3256 = _S3130 * _S3255 + 0.5f * (_S3128 * _S3255);
        _S3130 = 0.0f;
        _S3133 = _S3256;
    }
    else
    {
        _S3130 = _S3145;
        _S3133 = 0.0f;
    }
    DiffPair_float_0 _S3257;
    (&_S3257)->primal_0 = _S3128;
    (&_S3257)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3257, _S3130);
    float _S3258 = (_S3257.differential_0 + _S3133) / _S3109;
    FixedArray<float, 23>  _S3259;
    _S3259[int(0)] = 0.0f;
    _S3259[int(1)] = 0.0f;
    _S3259[int(2)] = 0.0f;
    _S3259[int(3)] = 0.0f;
    _S3259[int(4)] = 0.0f;
    _S3259[int(5)] = 0.0f;
    _S3259[int(6)] = 0.0f;
    _S3259[int(7)] = 0.0f;
    _S3259[int(8)] = 0.0f;
    _S3259[int(9)] = 0.0f;
    _S3259[int(10)] = 0.0f;
    _S3259[int(11)] = 0.0f;
    _S3259[int(12)] = 0.0f;
    _S3259[int(13)] = 0.0f;
    _S3259[int(14)] = 0.0f;
    _S3259[int(15)] = 0.0f;
    _S3259[int(16)] = 0.0f;
    _S3259[int(17)] = 0.0f;
    _S3259[int(18)] = 0.0f;
    _S3259[int(19)] = 0.0f;
    _S3259[int(20)] = 0.0f;
    _S3259[int(21)] = 0.0f;
    _S3259[int(22)] = 0.0f;
    _S3259[int(14)] = _S3258;
    float _S3260 = _S3232 + _S3259[int(0)];
    float _S3261 = _S3233 + _S3259[int(1)];
    float _S3262 = _S3234 + _S3259[int(2)];
    float _S3263 = _S3235 + _S3259[int(3)];
    float _S3264 = _S3236 + _S3259[int(4)];
    float _S3265 = _S3237 + _S3259[int(5)];
    float _S3266 = _S3238 + _S3259[int(6)];
    float _S3267 = _S3239 + _S3259[int(7)];
    float _S3268 = _S3240 + _S3259[int(8)];
    float _S3269 = _S3241 + _S3259[int(9)];
    float _S3270 = _S3242 + _S3259[int(10)];
    float _S3271 = _S3243 + _S3259[int(11)];
    float _S3272 = _S3244 + _S3259[int(12)];
    float _S3273 = _S3245 + _S3259[int(13)];
    float _S3274 = _S3246 + _S3259[int(14)];
    float _S3275 = _S3247 + _S3259[int(15)];
    float _S3276 = _S3248 + _S3259[int(16)];
    float _S3277 = _S3249 + _S3259[int(17)];
    float _S3278 = _S3250 + _S3259[int(18)];
    float _S3279 = _S3251 + _S3259[int(19)];
    float _S3280 = _S3252 + _S3259[int(20)];
    float _S3281 = _S3253 + _S3259[int(21)];
    float _S3282 = _S3254 + _S3259[int(22)];
    if(_S3126)
    {
        float _S3283 = 200.0f * _S3145;
        float _S3284 = _S3127 * _S3283 + 0.5f * (_S3125 * _S3283);
        _S3127 = 0.0f;
        _S3130 = _S3284;
    }
    else
    {
        _S3127 = _S3145;
        _S3130 = 0.0f;
    }
    DiffPair_float_0 _S3285;
    (&_S3285)->primal_0 = _S3125;
    (&_S3285)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3285, _S3127);
    float _S3286 = (_S3285.differential_0 + _S3130) / _S3109;
    FixedArray<float, 23>  _S3287;
    _S3287[int(0)] = 0.0f;
    _S3287[int(1)] = 0.0f;
    _S3287[int(2)] = 0.0f;
    _S3287[int(3)] = 0.0f;
    _S3287[int(4)] = 0.0f;
    _S3287[int(5)] = 0.0f;
    _S3287[int(6)] = 0.0f;
    _S3287[int(7)] = 0.0f;
    _S3287[int(8)] = 0.0f;
    _S3287[int(9)] = 0.0f;
    _S3287[int(10)] = 0.0f;
    _S3287[int(11)] = 0.0f;
    _S3287[int(12)] = 0.0f;
    _S3287[int(13)] = 0.0f;
    _S3287[int(14)] = 0.0f;
    _S3287[int(15)] = 0.0f;
    _S3287[int(16)] = 0.0f;
    _S3287[int(17)] = 0.0f;
    _S3287[int(18)] = 0.0f;
    _S3287[int(19)] = 0.0f;
    _S3287[int(20)] = 0.0f;
    _S3287[int(21)] = 0.0f;
    _S3287[int(22)] = 0.0f;
    _S3287[int(13)] = _S3286;
    float _S3288 = _S3260 + _S3287[int(0)];
    float _S3289 = _S3261 + _S3287[int(1)];
    float _S3290 = _S3262 + _S3287[int(2)];
    float _S3291 = _S3263 + _S3287[int(3)];
    float _S3292 = _S3264 + _S3287[int(4)];
    float _S3293 = _S3265 + _S3287[int(5)];
    float _S3294 = _S3266 + _S3287[int(6)];
    float _S3295 = _S3267 + _S3287[int(7)];
    float _S3296 = _S3268 + _S3287[int(8)];
    float _S3297 = _S3269 + _S3287[int(9)];
    float _S3298 = _S3270 + _S3287[int(10)];
    float _S3299 = _S3271 + _S3287[int(11)];
    float _S3300 = _S3272 + _S3287[int(12)];
    float _S3301 = _S3273 + _S3287[int(13)];
    float _S3302 = _S3274 + _S3287[int(14)];
    float _S3303 = _S3275 + _S3287[int(15)];
    float _S3304 = _S3276 + _S3287[int(16)];
    float _S3305 = _S3277 + _S3287[int(17)];
    float _S3306 = _S3278 + _S3287[int(18)];
    float _S3307 = _S3279 + _S3287[int(19)];
    float _S3308 = _S3280 + _S3287[int(20)];
    float _S3309 = _S3281 + _S3287[int(21)];
    float _S3310 = _S3282 + _S3287[int(22)];
    if(_S3123)
    {
        float _S3311 = 200.0f * _S3145;
        float _S3312 = _S3124 * _S3311 + 0.5f * (_S3122 * _S3311);
        _S3124 = 0.0f;
        _S3127 = _S3312;
    }
    else
    {
        _S3124 = _S3145;
        _S3127 = 0.0f;
    }
    DiffPair_float_0 _S3313;
    (&_S3313)->primal_0 = _S3122;
    (&_S3313)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3313, _S3124);
    float _S3314 = (_S3313.differential_0 + _S3127) / _S3109;
    FixedArray<float, 23>  _S3315;
    _S3315[int(0)] = 0.0f;
    _S3315[int(1)] = 0.0f;
    _S3315[int(2)] = 0.0f;
    _S3315[int(3)] = 0.0f;
    _S3315[int(4)] = 0.0f;
    _S3315[int(5)] = 0.0f;
    _S3315[int(6)] = 0.0f;
    _S3315[int(7)] = 0.0f;
    _S3315[int(8)] = 0.0f;
    _S3315[int(9)] = 0.0f;
    _S3315[int(10)] = 0.0f;
    _S3315[int(11)] = 0.0f;
    _S3315[int(12)] = 0.0f;
    _S3315[int(13)] = 0.0f;
    _S3315[int(14)] = 0.0f;
    _S3315[int(15)] = 0.0f;
    _S3315[int(16)] = 0.0f;
    _S3315[int(17)] = 0.0f;
    _S3315[int(18)] = 0.0f;
    _S3315[int(19)] = 0.0f;
    _S3315[int(20)] = 0.0f;
    _S3315[int(21)] = 0.0f;
    _S3315[int(22)] = 0.0f;
    _S3315[int(12)] = _S3314;
    float _S3316 = _S3288 + _S3315[int(0)];
    float _S3317 = _S3289 + _S3315[int(1)];
    float _S3318 = _S3290 + _S3315[int(2)];
    float _S3319 = _S3291 + _S3315[int(3)];
    float _S3320 = _S3292 + _S3315[int(4)];
    float _S3321 = _S3293 + _S3315[int(5)];
    float _S3322 = _S3294 + _S3315[int(6)];
    float _S3323 = _S3295 + _S3315[int(7)];
    float _S3324 = _S3296 + _S3315[int(8)];
    float _S3325 = _S3297 + _S3315[int(9)];
    float _S3326 = _S3298 + _S3315[int(10)];
    float _S3327 = _S3299 + _S3315[int(11)];
    float _S3328 = _S3300 + _S3315[int(12)];
    float _S3329 = _S3301 + _S3315[int(13)];
    float _S3330 = _S3302 + _S3315[int(14)];
    float _S3331 = _S3303 + _S3315[int(15)];
    float _S3332 = _S3304 + _S3315[int(16)];
    float _S3333 = _S3305 + _S3315[int(17)];
    float _S3334 = _S3306 + _S3315[int(18)];
    float _S3335 = _S3307 + _S3315[int(19)];
    float _S3336 = _S3308 + _S3315[int(20)];
    float _S3337 = _S3309 + _S3315[int(21)];
    float _S3338 = _S3310 + _S3315[int(22)];
    if(_S3120)
    {
        float _S3339 = 200.0f * _S3145;
        float _S3340 = _S3121 * _S3339 + 0.5f * (_S3119 * _S3339);
        _S3121 = 0.0f;
        _S3124 = _S3340;
    }
    else
    {
        _S3121 = _S3145;
        _S3124 = 0.0f;
    }
    DiffPair_float_0 _S3341;
    (&_S3341)->primal_0 = _S3119;
    (&_S3341)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3341, _S3121);
    float _S3342 = (_S3341.differential_0 + _S3124) / _S3109;
    FixedArray<float, 23>  _S3343;
    _S3343[int(0)] = 0.0f;
    _S3343[int(1)] = 0.0f;
    _S3343[int(2)] = 0.0f;
    _S3343[int(3)] = 0.0f;
    _S3343[int(4)] = 0.0f;
    _S3343[int(5)] = 0.0f;
    _S3343[int(6)] = 0.0f;
    _S3343[int(7)] = 0.0f;
    _S3343[int(8)] = 0.0f;
    _S3343[int(9)] = 0.0f;
    _S3343[int(10)] = 0.0f;
    _S3343[int(11)] = 0.0f;
    _S3343[int(12)] = 0.0f;
    _S3343[int(13)] = 0.0f;
    _S3343[int(14)] = 0.0f;
    _S3343[int(15)] = 0.0f;
    _S3343[int(16)] = 0.0f;
    _S3343[int(17)] = 0.0f;
    _S3343[int(18)] = 0.0f;
    _S3343[int(19)] = 0.0f;
    _S3343[int(20)] = 0.0f;
    _S3343[int(21)] = 0.0f;
    _S3343[int(22)] = 0.0f;
    _S3343[int(11)] = _S3342;
    float _S3344 = _S3316 + _S3343[int(0)];
    float _S3345 = _S3317 + _S3343[int(1)];
    float _S3346 = _S3318 + _S3343[int(2)];
    float _S3347 = _S3319 + _S3343[int(3)];
    float _S3348 = _S3320 + _S3343[int(4)];
    float _S3349 = _S3321 + _S3343[int(5)];
    float _S3350 = _S3322 + _S3343[int(6)];
    float _S3351 = _S3323 + _S3343[int(7)];
    float _S3352 = _S3324 + _S3343[int(8)];
    float _S3353 = _S3325 + _S3343[int(9)];
    float _S3354 = _S3326 + _S3343[int(10)];
    float _S3355 = _S3327 + _S3343[int(11)];
    float _S3356 = _S3328 + _S3343[int(12)];
    float _S3357 = _S3329 + _S3343[int(13)];
    float _S3358 = _S3330 + _S3343[int(14)];
    float _S3359 = _S3331 + _S3343[int(15)];
    float _S3360 = _S3332 + _S3343[int(16)];
    float _S3361 = _S3333 + _S3343[int(17)];
    float _S3362 = _S3334 + _S3343[int(18)];
    float _S3363 = _S3335 + _S3343[int(19)];
    float _S3364 = _S3336 + _S3343[int(20)];
    float _S3365 = _S3337 + _S3343[int(21)];
    float _S3366 = _S3338 + _S3343[int(22)];
    if(_S3117)
    {
        float _S3367 = 200.0f * _S3145;
        float _S3368 = _S3118 * _S3367 + 0.5f * (_S3116 * _S3367);
        _S3118 = 0.0f;
        _S3121 = _S3368;
    }
    else
    {
        _S3118 = _S3145;
        _S3121 = 0.0f;
    }
    DiffPair_float_0 _S3369;
    (&_S3369)->primal_0 = _S3116;
    (&_S3369)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3369, _S3118);
    float _S3370 = (_S3369.differential_0 + _S3121) / _S3109;
    float _S3371 = _S3140 / _S3115;
    float _S3372 = _S3141 / _S3114;
    float _S3373 = _S3142 / _S3113;
    FixedArray<float, 23>  _S3374;
    _S3374[int(0)] = 0.0f;
    _S3374[int(1)] = 0.0f;
    _S3374[int(2)] = 0.0f;
    _S3374[int(3)] = 0.0f;
    _S3374[int(4)] = 0.0f;
    _S3374[int(5)] = 0.0f;
    _S3374[int(6)] = 0.0f;
    _S3374[int(7)] = 0.0f;
    _S3374[int(8)] = 0.0f;
    _S3374[int(9)] = 0.0f;
    _S3374[int(10)] = 0.0f;
    _S3374[int(11)] = 0.0f;
    _S3374[int(12)] = 0.0f;
    _S3374[int(13)] = 0.0f;
    _S3374[int(14)] = 0.0f;
    _S3374[int(15)] = 0.0f;
    _S3374[int(16)] = 0.0f;
    _S3374[int(17)] = 0.0f;
    _S3374[int(18)] = 0.0f;
    _S3374[int(19)] = 0.0f;
    _S3374[int(20)] = 0.0f;
    _S3374[int(21)] = 0.0f;
    _S3374[int(22)] = 0.0f;
    _S3374[int(10)] = _S3370;
    _S3374[int(9)] = _S3371;
    _S3374[int(8)] = _S3371;
    _S3374[int(7)] = _S3371;
    _S3374[int(6)] = _S3371;
    _S3374[int(5)] = _S3371;
    _S3374[int(4)] = _S3372;
    _S3374[int(3)] = _S3372;
    _S3374[int(2)] = _S3372;
    _S3374[int(1)] = _S3373;
    float _S3375 = _S3344 + _S3374[int(0)];
    float _S3376 = _S3345 + _S3374[int(1)];
    float _S3377 = _S3346 + _S3374[int(2)];
    float _S3378 = _S3347 + _S3374[int(3)];
    float _S3379 = _S3348 + _S3374[int(4)];
    float _S3380 = _S3349 + _S3374[int(5)];
    float _S3381 = _S3350 + _S3374[int(6)];
    float _S3382 = _S3351 + _S3374[int(7)];
    float _S3383 = _S3352 + _S3374[int(8)];
    float _S3384 = _S3353 + _S3374[int(9)];
    float _S3385 = _S3354 + _S3374[int(10)];
    float _S3386 = _S3355 + _S3374[int(11)];
    float _S3387 = _S3356 + _S3374[int(12)];
    float _S3388 = _S3357 + _S3374[int(13)];
    float _S3389 = _S3358 + _S3374[int(14)];
    float _S3390 = _S3359 + _S3374[int(15)];
    float _S3391 = _S3360 + _S3374[int(16)];
    float _S3392 = _S3361 + _S3374[int(17)];
    float _S3393 = _S3362 + _S3374[int(18)];
    float _S3394 = _S3363 + _S3374[int(19)];
    float _S3395 = _S3364 + _S3374[int(20)];
    float _S3396 = _S3365 + _S3374[int(21)];
    float _S3397 = _S3366 + _S3374[int(22)];
    if(_S3111)
    {
        float _S3398 = 10.0f * _S3143;
        float _S3399 = _S3112 * _S3398 + 0.5f * (_S3110 * _S3398);
        _S3112 = 0.0f;
        _S3118 = _S3399;
    }
    else
    {
        _S3112 = _S3143;
        _S3118 = 0.0f;
    }
    DiffPair_float_0 _S3400;
    (&_S3400)->primal_0 = _S3110;
    (&_S3400)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3400, _S3112);
    float _S3401 = (_S3400.differential_0 + _S3118) / _S3109;
    FixedArray<float, 23>  _S3402;
    _S3402[int(0)] = 0.0f;
    _S3402[int(1)] = 0.0f;
    _S3402[int(2)] = 0.0f;
    _S3402[int(3)] = 0.0f;
    _S3402[int(4)] = 0.0f;
    _S3402[int(5)] = 0.0f;
    _S3402[int(6)] = 0.0f;
    _S3402[int(7)] = 0.0f;
    _S3402[int(8)] = 0.0f;
    _S3402[int(9)] = 0.0f;
    _S3402[int(10)] = 0.0f;
    _S3402[int(11)] = 0.0f;
    _S3402[int(12)] = 0.0f;
    _S3402[int(13)] = 0.0f;
    _S3402[int(14)] = 0.0f;
    _S3402[int(15)] = 0.0f;
    _S3402[int(16)] = 0.0f;
    _S3402[int(17)] = 0.0f;
    _S3402[int(18)] = 0.0f;
    _S3402[int(19)] = 0.0f;
    _S3402[int(20)] = 0.0f;
    _S3402[int(21)] = 0.0f;
    _S3402[int(22)] = 0.0f;
    _S3402[int(0)] = _S3401;
    FixedArray<float, 23>  _S3403 = {
        _S3375 + _S3402[int(0)], _S3376 + _S3402[int(1)], _S3377 + _S3402[int(2)], _S3378 + _S3402[int(3)], _S3379 + _S3402[int(4)], _S3380 + _S3402[int(5)], _S3381 + _S3402[int(6)], _S3382 + _S3402[int(7)], _S3383 + _S3402[int(8)], _S3384 + _S3402[int(9)], _S3385 + _S3402[int(10)], _S3386 + _S3402[int(11)], _S3387 + _S3402[int(12)], _S3388 + _S3402[int(13)], _S3389 + _S3402[int(14)], _S3390 + _S3402[int(15)], _S3391 + _S3402[int(16)], _S3392 + _S3402[int(17)], _S3393 + _S3402[int(18)], _S3394 + _S3402[int(19)], _S3395 + _S3402[int(20)], _S3396 + _S3402[int(21)], _S3397 + _S3402[int(22)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S3403;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S3404, int _S3405, FixedArray<float, 6>  * _S3406, FixedArray<float, 6>  * _S3407)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S3404, _S3405, _S3406, _S3407);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_5, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_5, FixedArray<float, 23>  * _S3408)
{
    FixedArray<float, 23>  _S3409 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S3409;
    FixedArray<float, 6>  _S3410 = loss_weights_5;
    FixedArray<float, 6>  _S3411 = grad_out_5;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_2, num_cameras_5, &_S3410, &_S3411);
    *_S3408 = (&dp_raw_losses_2)->differential_0;
    return;
}

