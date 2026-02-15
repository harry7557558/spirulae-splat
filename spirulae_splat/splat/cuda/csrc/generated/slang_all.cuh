#pragma once

// #if defined(CUDART_VERSION) && CUDART_VERSION >= 12000
// typedef double4_32a my_double4;
// typedef longlong4_32a my_longlong4;
// typedef ulonglong4_32a my_ulonglong4;
// #else
// typedef double4 my_double4;
// typedef longlong4 my_longlong4;
// typedef ulonglong4 my_ulonglong4;
// #endif

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

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_2)
{
    float x_13 = quat_2.y;
    float x2_0 = x_13 * x_13;
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_8)
{
    float _S278 = (*left_0).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_8.x;
    float sum_0 = _S278 + (*left_0).primal_0.rows[int(1)].x * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_8.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_8.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S279 = (*left_0).primal_0.rows[int(0)].y * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_8.x;
    float sum_2 = _S279 + (*left_0).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_8.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_8.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S280 = (*left_0).primal_0.rows[int(0)].z * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_8.x;
    float sum_4 = _S280 + (*left_0).primal_0.rows[int(1)].z * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_8.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_8.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

struct DiffPair_matrixx3Cfloatx2C2x2C2x3E_0
{
    Matrix<float, 2, 2>  primal_0;
    Matrix<float, 2, 2>  differential_0;
};

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * left_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * right_1, float2  dOut_9)
{
    float _S281 = (*left_1).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 2, 2>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_1).primal_0.x * dOut_9.x;
    float sum_6 = _S281 + (*left_1).primal_0.rows[int(1)].x * dOut_9.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_1).primal_0.x * dOut_9.y;
    float2  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_6;
    float _S282 = (*left_1).primal_0.rows[int(0)].y * dOut_9.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_1).primal_0.y * dOut_9.x;
    float sum_7 = _S282 + (*left_1).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_1).primal_0.y * dOut_9.y;
    *&((&right_d_result_1)->y) = sum_7;
    left_1->primal_0 = (*left_1).primal_0;
    left_1->differential_0 = left_d_result_1;
    right_1->primal_0 = (*right_1).primal_0;
    right_1->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_2, float3  right_2)
{
    float3  result_15;
    int i_5 = int(0);
    for(;;)
    {
        if(i_5 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_2.rows[i_5], j_0) * _slang_vector_get_element(right_2, j_0);
            j_0 = j_0 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_15, i_5) = sum_8;
        i_5 = i_5 + int(1);
    }
    return result_15;
}

inline __device__ float2  mul_1(Matrix<float, 2, 2>  left_3, float2  right_3)
{
    float2  result_16;
    int i_6 = int(0);
    for(;;)
    {
        if(i_6 < int(2))
        {
        }
        else
        {
            break;
        }
        int j_1 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(j_1 < int(2))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_3.rows[i_6], j_1) * _slang_vector_get_element(right_3, j_1);
            j_1 = j_1 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_16, i_6) = sum_10;
        i_6 = i_6 + int(1);
    }
    return result_16;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_4, Matrix<float, 3, 3>  dOut_10)
{
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_2;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_10.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(0)].x * dOut_10.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_10.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(0)].y * dOut_10.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_10.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(0)].z * dOut_10.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_10.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(0)].x * dOut_10.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_10.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(0)].y * dOut_10.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_10.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(0)].z * dOut_10.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_10.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(0)].x * dOut_10.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_10.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(0)].y * dOut_10.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_10.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(0)].z * dOut_10.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_10.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(1)].x * dOut_10.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_10.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(1)].y * dOut_10.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_10.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(1)].z * dOut_10.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_10.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(1)].x * dOut_10.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_10.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(1)].y * dOut_10.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_10.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(1)].z * dOut_10.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_10.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(1)].x * dOut_10.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_10.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(1)].y * dOut_10.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_10.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(1)].z * dOut_10.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_10.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(2)].x * dOut_10.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_10.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(2)].y * dOut_10.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_10.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(2)].z * dOut_10.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_10.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(2)].x * dOut_10.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_10.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(2)].y * dOut_10.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_10.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(2)].z * dOut_10.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_10.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(2)].x * dOut_10.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_10.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(2)].y * dOut_10.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_10.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(2)].z * dOut_10.rows[int(2)].z;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_2;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_2;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_3(Matrix<float, 3, 3>  left_5, Matrix<float, 3, 3>  right_5)
{
    Matrix<float, 3, 3>  result_17;
    int r_2 = int(0);
    for(;;)
    {
        if(r_2 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_1 = int(0);
        for(;;)
        {
            if(c_1 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_7 = int(0);
            float sum_12 = 0.0f;
            for(;;)
            {
                if(i_7 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_13 = sum_12 + _slang_vector_get_element(left_5.rows[r_2], i_7) * _slang_vector_get_element(right_5.rows[i_7], c_1);
                i_7 = i_7 + int(1);
                sum_12 = sum_13;
            }
            *_slang_vector_get_element_ptr(((&result_17)->rows + (r_2)), c_1) = sum_12;
            c_1 = c_1 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_17;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_3(mul_3(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_3, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_14 = quat_3.y;
    float x2_1 = x_14 * x_14;
    float y2_1 = quat_3.z * quat_3.z;
    float z2_1 = quat_3.w * quat_3.w;
    float xy_1 = quat_3.y * quat_3.z;
    float xz_1 = quat_3.y * quat_3.w;
    float yz_1 = quat_3.z * quat_3.w;
    float wx_1 = quat_3.x * quat_3.y;
    float wy_1 = quat_3.x * quat_3.z;
    float wz_1 = quat_3.x * quat_3.w;
    Matrix<float, 3, 3>  M_0 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_3(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_4, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_15 = quat_4.y;
    float x2_2 = x_15 * x_15;
    float y2_2 = quat_4.z * quat_4.z;
    float z2_2 = quat_4.w * quat_4.w;
    float xy_2 = quat_4.y * quat_4.z;
    float xz_2 = quat_4.y * quat_4.w;
    float yz_2 = quat_4.z * quat_4.w;
    float wx_2 = quat_4.x * quat_4.y;
    float wy_2 = quat_4.x * quat_4.z;
    float wz_2 = quat_4.x * quat_4.w;
    *M_1 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

inline __device__ Matrix<float, 2, 2>  camera_distortion_jac_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float _S283 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S284 = (*dist_coeffs_0)[int(1)] + r2_0 * _S283;
    float _S285 = (*dist_coeffs_0)[int(0)] + r2_0 * _S284;
    float2  _S286 = make_float2 (1.0f + r2_0 * _S285);
    float _S287 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S288 = _S287 * u_0;
    float _S289 = 2.0f * u_0;
    float _S290 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S291 = _S290 * u_0;
    float _S292 = 2.0f * v_0;
    float2  _S293 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S294 = uv_0 * _S293;
    float _S295 = (*dist_coeffs_0)[int(4)] * _S293.y;
    float _S296 = (*dist_coeffs_0)[int(5)] * _S293.x;
    float _S297 = _S294.x + _S294.y;
    float _S298 = r2_0 * _S297;
    float _S299 = r2_0 * _S298;
    float _S300 = (*dist_coeffs_0)[int(7)] * _S293.y + _S295 + (*dist_coeffs_0)[int(6)] * _S293.x + _S296 + _S285 * _S297 + _S284 * _S298 + _S283 * _S299 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S299);
    float _S301 = v_0 * _S300;
    float _S302 = u_0 * _S300;
    float2  _S303 = make_float2 (0.0f, 1.0f) + make_float2 ((*dist_coeffs_0)[int(8)] * 0.0f, (*dist_coeffs_0)[int(9)] * 0.0f);
    float2  _S304 = uv_0 * _S303;
    float _S305 = (*dist_coeffs_0)[int(4)] * _S303.y;
    float _S306 = (*dist_coeffs_0)[int(5)] * _S303.x;
    float _S307 = _S304.x + _S304.y;
    float _S308 = r2_0 * _S307;
    float _S309 = r2_0 * _S308;
    float _S310 = (*dist_coeffs_0)[int(7)] * _S303.y + _S305 + (*dist_coeffs_0)[int(6)] * _S303.x + _S306 + _S285 * _S307 + _S284 * _S308 + _S283 * _S309 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S309);
    float _S311 = v_0 * _S310;
    float _S312 = u_0 * _S310;
    return makeMatrix<float, 2, 2> (_S286 * _S293 + make_float2 (_S290 * (v_0 * _S293.y) + _S289 * _S296 + 2.0f * (u_0 * _S296) + _S287 * (v_0 * _S293.x) + _S302 + _S302, _S292 * _S295 + 2.0f * (v_0 * _S295) + _S291 * _S293.y + _S288 * _S293.x + _S301 + _S301), _S286 * _S303 + make_float2 (_S290 * (v_0 * _S303.y) + _S289 * _S306 + 2.0f * (u_0 * _S306) + _S287 * (v_0 * _S303.x) + _S312 + _S312, _S292 * _S305 + 2.0f * (v_0 * _S305) + _S291 * _S303.y + _S288 * _S303.x + _S311 + _S311));
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  dist_coeffs_1)
{
    FixedArray<float, 10>  _S313 = dist_coeffs_1;
    Matrix<float, 2, 2>  _S314 = camera_distortion_jac_0(uv_1, &_S313);
    return (F32_min((determinant_0(_S314)), ((F32_min((_S314.rows[int(0)].x), (_S314.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_2)
{
    float2  _S315;
    if(is_fisheye_0)
    {
        float r_3 = length_1(uv_2);
        float theta_0 = (F32_atan((r_3)));
        float _S316;
        if(r_3 < 0.00100000004749745f)
        {
            _S316 = 1.0f - r_3 * r_3 / 3.0f;
        }
        else
        {
            _S316 = theta_0 / r_3;
        }
        _S315 = uv_2 * make_float2 (_S316);
    }
    else
    {
        _S315 = uv_2;
    }
    float u_1 = _S315.x;
    float v_1 = _S315.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S317 = _S315 * make_float2 (1.0f + r2_1 * (dist_coeffs_2[int(0)] + r2_1 * (dist_coeffs_2[int(1)] + r2_1 * (dist_coeffs_2[int(2)] + r2_1 * dist_coeffs_2[int(3)])))) + make_float2 (2.0f * dist_coeffs_2[int(4)] * u_1 * v_1 + dist_coeffs_2[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_2[int(6)] * r2_1, 2.0f * dist_coeffs_2[int(5)] * u_1 * v_1 + dist_coeffs_2[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_2[int(7)] * r2_1);
    return _S317 + make_float2 (dist_coeffs_2[int(8)] * _S317.x + dist_coeffs_2[int(9)] * _S317.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_3, FixedArray<float, 10>  * dist_coeffs_3, int maxiter_0, float2  * uv_undist_0)
{
    int i_8 = int(0);
    float2  q_0 = uv_3;
    for(;;)
    {
        if(i_8 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S318 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_4 = _S318 + make_float2 ((*dist_coeffs_3)[int(8)] * _S318.x + (*dist_coeffs_3)[int(9)] * _S318.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S319 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S319.rows[int(0)].x * _S319.rows[int(1)].y - _S319.rows[int(0)].y * _S319.rows[int(1)].x);
        float _S320 = r_4.x;
        float _S321 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S320 * _S319.rows[int(1)].y - _S321 * _S319.rows[int(0)].y) * inv_det_0, (- _S320 * _S319.rows[int(1)].x + _S321 * _S319.rows[int(0)].x) * inv_det_0);
        i_8 = i_8 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S322 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S323;
    if((F32_min((determinant_0(_S322)), ((F32_min((_S322.rows[int(0)].x), (_S322.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S324 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S323 = (length_1(_S324 + make_float2 ((*dist_coeffs_3)[int(8)] * _S324.x + (*dist_coeffs_3)[int(9)] * _S324.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S323 = false;
    }
    return _S323;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S325 = uv_4;
    FixedArray<float, 10>  _S326 = dist_coeffs_4;
    bool _S327 = undistort_point_0(uv_4, &_S326, int(8), &_S325);
    if(!_S327)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S328 = _S325;
        float theta_1 = length_1(_S325);
        float _S329;
        if(theta_1 < 0.00100000004749745f)
        {
            _S329 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S329 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S330 = make_float3 ((_S328 * make_float2 (_S329)).x, (_S328 * make_float2 (_S329)).y, (F32_cos((theta_1))));
        raydir_0 = _S330;
    }
    else
    {
        raydir_0 = make_float3 (_S325.x, _S325.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_5, float3  * raydir_1)
{
    float2  _S331 = uv_5;
    int3  _S332 = make_int3 (int(0));
    float3  _S333 = make_float3 ((float)_S332.x, (float)_S332.y, (float)_S332.z);
    *raydir_1 = _S333;
    FixedArray<float, 10>  _S334 = dist_coeffs_5;
    bool _S335 = undistort_point_0(uv_5, &_S334, int(8), &_S331);
    if(!_S335)
    {
        return false;
    }
    float3  _S336;
    if(is_fisheye_2)
    {
        float2  _S337 = _S331;
        float theta_2 = length_1(_S331);
        float _S338;
        if(theta_2 < 0.00100000004749745f)
        {
            _S338 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S338 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S339 = make_float3 ((_S337 * make_float2 (_S338)).x, (_S337 * make_float2 (_S338)).y, (F32_cos((theta_2))));
        _S336 = _S339;
    }
    else
    {
        _S336 = make_float3 (_S331.x, _S331.y, 1.0f);
    }
    *raydir_1 = _S336;
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

inline __device__ bool generate_ray(float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_6, float3  * raydir_2)
{
    float2  _S340 = uv_6;
    FixedArray<float, 10>  _S341 = dist_coeffs_6;
    bool _S342 = undistort_point_0(uv_6, &_S341, int(8), &_S340);
    if(!_S342)
    {
        int3  _S343 = make_int3 (int(0));
        float3  _S344 = make_float3 ((float)_S343.x, (float)_S343.y, (float)_S343.z);
        *raydir_2 = _S344;
        return false;
    }
    float3  _S345;
    if(is_fisheye_3)
    {
        float2  _S346 = _S340;
        float theta_3 = length_1(_S340);
        float _S347;
        if(theta_3 < 0.00100000004749745f)
        {
            _S347 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S347 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S348 = make_float3 ((_S346 * make_float2 (_S347)).x, (_S346 * make_float2 (_S347)).y, (F32_cos((theta_3))));
        _S345 = _S348;
    }
    else
    {
        _S345 = make_float3 (_S340.x, _S340.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S345);
    return true;
}

inline __device__ void _d_mul_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_6, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_6, float3  dOut_11)
{
    float _S349 = (*right_6).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 3, 3>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = (*left_6).primal_0.x * dOut_11.x;
    float sum_14 = _S349 + (*right_6).primal_0.rows[int(0)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = (*left_6).primal_0.x * dOut_11.y;
    float sum_15 = sum_14 + (*right_6).primal_0.rows[int(0)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(0)))->z) = (*left_6).primal_0.x * dOut_11.z;
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = sum_15;
    float _S350 = (*right_6).primal_0.rows[int(1)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = (*left_6).primal_0.y * dOut_11.x;
    float sum_16 = _S350 + (*right_6).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = (*left_6).primal_0.y * dOut_11.y;
    float sum_17 = sum_16 + (*right_6).primal_0.rows[int(1)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(1)))->z) = (*left_6).primal_0.y * dOut_11.z;
    *&((&left_d_result_3)->y) = sum_17;
    float _S351 = (*right_6).primal_0.rows[int(2)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = (*left_6).primal_0.z * dOut_11.x;
    float sum_18 = _S351 + (*right_6).primal_0.rows[int(2)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = (*left_6).primal_0.z * dOut_11.y;
    float sum_19 = sum_18 + (*right_6).primal_0.rows[int(2)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(2)))->z) = (*left_6).primal_0.z * dOut_11.z;
    *&((&left_d_result_3)->z) = sum_19;
    left_6->primal_0 = (*left_6).primal_0;
    left_6->differential_0 = left_d_result_3;
    right_6->primal_0 = (*right_6).primal_0;
    right_6->differential_0 = right_d_result_3;
    return;
}

inline __device__ float3  mul_4(float3  left_7, Matrix<float, 3, 3>  right_7)
{
    float3  result_18;
    int j_2 = int(0);
    for(;;)
    {
        if(j_2 < int(3))
        {
        }
        else
        {
            break;
        }
        int i_9 = int(0);
        float sum_20 = 0.0f;
        for(;;)
        {
            if(i_9 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_21 = sum_20 + _slang_vector_get_element(left_7, i_9) * _slang_vector_get_element(right_7.rows[i_9], j_2);
            i_9 = i_9 + int(1);
            sum_20 = sum_21;
        }
        *_slang_vector_get_element_ptr(&result_18, j_2) = sum_20;
        j_2 = j_2 + int(1);
    }
    return result_18;
}

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_2, float3  t_1)
{
    return - mul_4(t_1, R_2);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_3, float3  raydir_3)
{
    return mul_4(raydir_3, R_3);
}

inline __device__ float3  undo_transform_ray_d(Matrix<float, 3, 3>  R_4, float3  raydir_4)
{
    return mul_4(raydir_4, transpose_0(R_4));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S352, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S353, float3  _S354)
{
    _d_mul_2(_S352, _S353, _S354);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S355 = - _s_dOut_2;
    float3  _S356 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S357;
    (&_S357)->primal_0 = (*dpt_0).primal_0;
    (&_S357)->differential_0 = _S356;
    Matrix<float, 3, 3>  _S358 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S359;
    (&_S359)->primal_0 = (*dpR_0).primal_0;
    (&_S359)->differential_0 = _S358;
    s_bwd_prop_mul_0(&_S357, &_S359, _S355);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S357.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S359.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S360, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S361, float3  _S362)
{
    s_bwd_prop_transform_ray_o_0(_S360, _S361, _S362);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_5, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S363 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_5;
    (&dp_R_0)->differential_0 = _S363;
    float3  _S364 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S364;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S365 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S366;
    (&_S366)->primal_0 = (*dpraydir_0).primal_0;
    (&_S366)->differential_0 = _S365;
    Matrix<float, 3, 3>  _S367 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S368;
    (&_S368)->primal_0 = (*dpR_1).primal_0;
    (&_S368)->differential_0 = _S367;
    s_bwd_prop_mul_0(&_S366, &_S368, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S366.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S368.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S369, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S370, float3  _S371)
{
    s_bwd_prop_transform_ray_d_0(_S369, _S370, _S371);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_6, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S372 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_6;
    (&dp_R_1)->differential_0 = _S372;
    float3  _S373 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S373;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S374 = scale_2.x;
    float sx_0 = (F32_exp((_S374)));
    float _S375 = scale_2.y;
    float sy_0 = (F32_exp((_S375)));
    float sz_0 = scale_2.z - 0.5f * (_S374 + _S375);
    float x_18 = quat_5.y;
    float x2_3 = x_18 * x_18;
    float y2_3 = quat_5.z * quat_5.z;
    float z2_3 = quat_5.w * quat_5.w;
    float xy_3 = quat_5.y * quat_5.z;
    float xz_3 = quat_5.y * quat_5.w;
    float yz_3 = quat_5.z * quat_5.w;
    float wx_3 = quat_5.x * quat_5.y;
    float wy_3 = quat_5.x * quat_5.z;
    float wz_3 = quat_5.x * quat_5.w;
    Matrix<float, 3, 3>  _S376 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S376, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S376, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S376, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
    return;
}

inline __device__ float4  floor_0(float4  x_19)
{
    float4  result_19;
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
        *_slang_vector_get_element_ptr(&result_19, i_10) = (F32_floor((_slang_vector_get_element(x_19, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_19;
}

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_1, float3  scale_3, float4  quat_6, float opac_0)
{
    float4  _S377 = normalize_1(quat_6);
    float3  _S378 = exp_0(scale_3);
    float x_20 = _S377.y;
    float x2_4 = x_20 * x_20;
    float y2_4 = _S377.z * _S377.z;
    float z2_4 = _S377.w * _S377.w;
    float xy_4 = _S377.y * _S377.z;
    float xz_4 = _S377.y * _S377.w;
    float yz_4 = _S377.z * _S377.w;
    float wx_4 = _S377.x * _S377.y;
    float wy_4 = _S377.x * _S377.z;
    float wz_4 = _S377.x * _S377.w;
    Matrix<float, 3, 3>  M_2 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S378.x, 0.0f, 0.0f, 0.0f, _S378.y, 0.0f, 0.0f, 0.0f, _S378.z));
    float4  _S379 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S380 = _S379 - floor_0(_S379);
    float4  _S381 = _S380 + make_float4 (dot_1(_S380, float4 {_S380.w, _S380.z, _S380.x, _S380.y} + make_float4 (33.3300018310546875f)));
    float4  _S382 = (float4 {_S381.x, _S381.x, _S381.y, _S381.z} + float4 {_S381.y, _S381.z, _S381.z, _S381.w}) * float4 {_S381.z, _S381.y, _S381.w, _S381.x};
    float4  _S383 = _S382 - floor_0(_S382);
    float2  _S384 = float2 {_S383.x, _S383.z};
    float _S385 = 6.28318548202514648f * _S384.y;
    float2  _S386 = float2 {_S383.y, _S383.w};
    float _S387 = 6.28318548202514648f * _S386.y;
    *mean_1 = *mean_1 + mul_0(mul_3(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S384.x))))))) * make_float2 ((F32_cos((_S385))), (F32_sin((_S385))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S384.x))))))) * make_float2 ((F32_cos((_S385))), (F32_sin((_S385))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S386.x))))))) * make_float2 ((F32_cos((_S387))), (F32_sin((_S387))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_12)
{
    float _S388 = dOut_12.y;
    float _S389 = dOut_12.z;
    float _S390 = dOut_12.x;
    float _S391 = (*a_0).primal_0.z * _S388 + - (*a_0).primal_0.y * _S389;
    float _S392 = - (*a_0).primal_0.z * _S390 + (*a_0).primal_0.x * _S389;
    float _S393 = (*a_0).primal_0.y * _S390 + - (*a_0).primal_0.x * _S388;
    float3  _S394 = make_float3 (- (*b_1).primal_0.z * _S388 + (*b_1).primal_0.y * _S389, (*b_1).primal_0.z * _S390 + - (*b_1).primal_0.x * _S389, - (*b_1).primal_0.y * _S390 + (*b_1).primal_0.x * _S388);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S394;
    float3  _S395 = make_float3 (_S391, _S392, _S393);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S395;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S396 = left_8.y;
    float _S397 = right_8.z;
    float _S398 = left_8.z;
    float _S399 = right_8.y;
    float _S400 = right_8.x;
    float _S401 = left_8.x;
    return make_float3 (_S396 * _S397 - _S398 * _S399, _S398 * _S400 - _S401 * _S397, _S401 * _S399 - _S396 * _S400);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S402 = normalize_1(quat_7);
    float _S403 = scale_4.x;
    float sx_1 = (F32_exp((_S403)));
    float _S404 = scale_4.y;
    float sy_1 = (F32_exp((_S404)));
    float sz_1 = scale_4.z - 0.5f * (_S403 + _S404);
    float x_21 = _S402.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S402.z * _S402.z;
    float z2_5 = _S402.w * _S402.w;
    float xy_5 = _S402.y * _S402.z;
    float xz_5 = _S402.y * _S402.w;
    float yz_5 = _S402.z * _S402.w;
    float wx_5 = _S402.x * _S402.y;
    float wy_5 = _S402.x * _S402.z;
    float wz_5 = _S402.x * _S402.w;
    Matrix<float, 3, 3>  _S405 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S405, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S405, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S405, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S406 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S407 = _S406 - floor_0(_S406);
    float4  _S408 = _S407 + make_float4 (dot_1(_S407, float4 {_S407.w, _S407.z, _S407.x, _S407.y} + make_float4 (33.3300018310546875f)));
    float4  _S409 = (float4 {_S408.x, _S408.x, _S408.y, _S408.z} + float4 {_S408.y, _S408.z, _S408.z, _S408.w}) * float4 {_S408.z, _S408.y, _S408.w, _S408.x};
    float4  _S410 = _S409 - floor_0(_S409);
    float2  _S411 = float2 {_S410.x, _S410.z};
    float _S412 = 6.28318548202514648f * _S411.y;
    float2  _S413 = float2 {_S410.y, _S410.w};
    float _S414 = 6.28318548202514648f * _S413.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S411.x))))))) * make_float2 ((F32_cos((_S412))), (F32_sin((_S412))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S411.x))))))) * make_float2 ((F32_cos((_S412))), (F32_sin((_S412))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S413.x))))))) * make_float2 ((F32_cos((_S414))), (F32_sin((_S414))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_9, float dOut_13)
{
    float _S415 = _slang_select(((*dpx_9).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_9).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S415;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S416 = _slang_select(((*dpx_10).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_10).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S416;
    return;
}

inline __device__ float3  abs_0(float3  x_22)
{
    float3  result_20;
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
        *_slang_vector_get_element_ptr(&result_20, i_11) = (F32_abs((_slang_vector_get_element(x_22, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_20;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_15)
{
    DiffPair_float_0 _S417 = *dpx_11;
    bool _S418;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S418 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S418 = false;
    }
    float _S419;
    if(_S418)
    {
        _S419 = dOut_15;
    }
    else
    {
        _S419 = 0.0f;
    }
    dpx_11->primal_0 = _S417.primal_0;
    dpx_11->differential_0 = _S419;
    DiffPair_float_0 _S420 = *dpMin_0;
    if((_S417.primal_0) < ((*dpMin_0).primal_0))
    {
        _S419 = dOut_15;
    }
    else
    {
        _S419 = 0.0f;
    }
    dpMin_0->primal_0 = _S420.primal_0;
    dpMin_0->differential_0 = _S419;
    DiffPair_float_0 _S421 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S419 = dOut_15;
    }
    else
    {
        _S419 = 0.0f;
    }
    dpMax_0->primal_0 = _S421.primal_0;
    dpMax_0->differential_0 = _S419;
    return;
}

inline __device__ float clamp_0(float x_23, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_23), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_12, float dOut_16)
{
    float _S422 = -0.5f / ((*dpx_12).primal_0 * (F32_sqrt(((*dpx_12).primal_0)))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S422;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_17)
{
    float _S423 = (1.0f - (*dps_0).primal_0) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S423;
    DiffPair_float_0 _S424 = *dpy_3;
    float _S425 = (*dps_0).primal_0 * dOut_17;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S425;
    float _S426 = (_S424.primal_0 - (*dpx_13).primal_0) * dOut_17;
    dps_0->primal_0 = _S424.primal_0;
    dps_0->differential_0 = _S426;
    return;
}

inline __device__ float lerp_0(float x_24, float y_7, float s_4)
{
    return x_24 + (y_7 - x_24) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  weights_0, FixedArray<float, 23>  * _S427)
{
    float3  _S428;
    bool _S429;
    bool _S430;
    FixedArray<float, 23>  losses_1;
    float _S431 = float(mask_0);
    float3  _S432 = ref_rgb_0 - render_rgb_0;
    float3  _S433 = abs_0(_S432);
    losses_1[int(0)] = weights_0[int(0)] * _S431 * ((_S433.x + _S433.y + _S433.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S431 * clamp_0(dot_0(_S432, _S432) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S434 = float(depth_mask_0 & mask_0);
    float _S435 = _S434 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S436 = _S434 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S435;
    losses_1[int(3)] = _S436;
    losses_1[int(4)] = _S435 * _S435;
    losses_1[int(5)] = _S436 * _S436;
    losses_1[int(6)] = _S435 * _S436;
    bool _S437 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S438 = norm2_0 == 0.0f;
        _S429 = _S438;
        if(_S438)
        {
            _S428 = make_float3 (0.0f);
            break;
        }
        _S428 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S439;
    bool _S440 = !_S429;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S441 = norm2_1 == 0.0f;
        _S430 = _S441;
        if(_S441)
        {
            _S439 = make_float3 (0.0f);
            break;
        }
        _S439 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S442;
    float3  _S443;
    bool _S444 = !_S430;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S443 = make_float3 (0.0f);
            _S442 = false;
            break;
        }
        _S443 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S442 = _S437;
        break;
    }
    float _S445 = float(_S440 & _S442);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S428, _S443);
    losses_1[int(7)] = weights_0[int(2)] * _S445 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S446 = float(_S444 & _S442);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S439, _S443);
    losses_1[int(8)] = weights_0[int(2)] * _S446 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S447 = float(_S440 & _S444);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S428, _S439);
    losses_1[int(11)] = weights_0[int(5)] * _S447 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S448 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S449 = float(alpha_mask_0);
    float _S450 = float(ref_alpha_0);
    float _S451 = (F32_max((_S448), (_S450)));
    losses_1[int(9)] = weights_0[int(3)] * _S449 * - lerp_0((F32_log(((F32_max((1.0f - _S451), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S451), (9.99999997475242708e-07f)))))), _S450);
    float _S452 = 1.0f - _S448;
    float _S453 = 1.0f - _S450;
    float _S454 = (F32_max((_S452), (_S453)));
    losses_1[int(10)] = weights_0[int(4)] * _S449 * - lerp_0((F32_log(((F32_max((1.0f - _S454), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S454), (9.99999997475242708e-07f)))))), _S453);
    losses_1[int(12)] = weights_0[int(6)] * 4.0f * _S448 * _S452;
    float _S455 = (F32_max((_S448), (9.999999960041972e-13f)));
    losses_1[int(13)] = weights_0[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S455;
    losses_1[int(14)] = weights_0[int(8)] * depth_dist_0 / _S455;
    losses_1[int(15)] = weights_0[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S455;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S431;
    losses_1[int(18)] = _S434;
    losses_1[int(19)] = _S445;
    losses_1[int(20)] = _S446;
    losses_1[int(21)] = _S447;
    losses_1[int(22)] = _S449;
    *_S427 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S456, float3  _S457)
{
    return dot_0(_S456, _S457);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S458)
{
    return (F32_rsqrt((_S458)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S459, float _S460, float _S461)
{
    return clamp_0(_S459, _S460, _S461);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S462, DiffPair_float_0 * _S463, DiffPair_float_0 * _S464, float _S465)
{
    _d_lerp_0(_S462, _S463, _S464, _S465);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S466, DiffPair_float_0 * _S467, DiffPair_float_0 * _S468, float _S469)
{
    _d_clamp_0(_S466, _S467, _S468, _S469);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S470, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S471, float _S472)
{
    _d_dot_0(_S470, _S471, _S472);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S473, float _S474)
{
    _d_rsqrt_0(_S473, _S474);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S475, float3  _S476)
{
    _d_abs_vector_0(_S475, _S476);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_float_0 _S477 = *dprender_depth_0;
    DiffPair_float_0 _S478 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S479 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S480 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481 = *dpref_normal_0;
    DiffPair_float_0 _S482 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S483 = *dprgb_dist_0;
    DiffPair_float_0 _S484 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S485 = *dpnormal_dist_0;
    float3  _S486 = make_float3 (0.0f);
    float _S487 = float(mask_1);
    float _S488 = (*weights_1)[int(0)] * _S487;
    float3  _S489 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S490 = s_primal_ctx_dot_0(_S489, _S489) * 0.3333333432674408f;
    float _S491 = float(depth_mask_1 & mask_1);
    float _S492 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S493 = _S491 * s_primal_ctx_log_0(_S492);
    float _S494 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S495 = _S491 * s_primal_ctx_log_0(_S494);
    bool _S496 = normal_mask_1 & mask_1;
    float _S497 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S498 = _S497 == 0.0f;
    float3  _S499;
    if(_S498)
    {
        _S499 = make_float3 (0.0f);
    }
    bool _S500 = !_S498;
    float3  _S501;
    if(_S500)
    {
        float _S502 = s_primal_ctx_rsqrt_0(_S497);
        float3  _S503 = make_float3 (_S502);
        _S499 = _S479.primal_0 * make_float3 (_S502);
        _S501 = _S503;
    }
    else
    {
        _S501 = _S486;
    }
    float _S504 = s_primal_ctx_dot_0(_S480.primal_0, _S480.primal_0);
    bool _S505 = _S504 == 0.0f;
    float3  _S506;
    if(_S505)
    {
        _S506 = make_float3 (0.0f);
    }
    bool _S507 = !_S505;
    float3  _S508;
    if(_S507)
    {
        float _S509 = s_primal_ctx_rsqrt_0(_S504);
        float3  _S510 = make_float3 (_S509);
        _S506 = _S480.primal_0 * make_float3 (_S509);
        _S508 = _S510;
    }
    else
    {
        _S508 = _S486;
    }
    float _S511 = s_primal_ctx_dot_0(_S481.primal_0, _S481.primal_0);
    bool _S512 = _S511 == 0.0f;
    float3  _S513;
    bool _S514;
    if(_S512)
    {
        float3  _S515 = make_float3 (0.0f);
        _S514 = false;
        _S513 = _S515;
    }
    else
    {
        _S514 = _S496;
    }
    bool _S516 = !_S512;
    float3  _S517;
    if(_S516)
    {
        float _S518 = s_primal_ctx_rsqrt_0(_S511);
        float3  _S519 = make_float3 (_S518);
        _S513 = _S481.primal_0 * make_float3 (_S518);
        _S517 = _S519;
    }
    else
    {
        _S517 = _S486;
    }
    float _S520 = (*weights_1)[int(2)] * float(_S500 & _S514);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S499, _S513);
    float _S521 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S522 = (*weights_1)[int(2)] * float(_S507 & _S514);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S506, _S513);
    float _S523 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S524 = (*weights_1)[int(5)] * float(_S500 & _S507);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S499, _S506);
    float _S525 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S526 = s_primal_ctx_clamp_0(_S482.primal_0, 0.0f, 1.0f);
    float _S527 = float(alpha_mask_1);
    float _S528 = (*weights_1)[int(3)] * _S527;
    float _S529 = float(ref_alpha_1);
    float _S530 = (F32_max((_S526), (_S529)));
    float _S531 = 1.0f - _S530;
    float _S532 = (F32_max((_S531), (9.99999997475242708e-07f)));
    float _S533 = s_primal_ctx_log_0(_S532);
    float _S534 = (F32_max((_S530), (9.99999997475242708e-07f)));
    float _S535 = s_primal_ctx_log_0(_S534);
    float _S536 = (*weights_1)[int(4)] * _S527;
    float _S537 = 1.0f - _S526;
    float _S538 = 1.0f - _S529;
    float _S539 = (F32_max((_S537), (_S538)));
    float _S540 = 1.0f - _S539;
    float _S541 = (F32_max((_S540), (9.99999997475242708e-07f)));
    float _S542 = s_primal_ctx_log_0(_S541);
    float _S543 = (F32_max((_S539), (9.99999997475242708e-07f)));
    float _S544 = s_primal_ctx_log_0(_S543);
    float _S545 = (*weights_1)[int(6)] * 4.0f;
    float _S546 = _S545 * _S526;
    float _S547 = (F32_max((_S526), (9.999999960041972e-13f)));
    float _S548 = _S547 * _S547;
    float _S549 = (*_s_dOut_4)[int(0)];
    float _S550 = (*_s_dOut_4)[int(1)];
    float _S551 = (*_s_dOut_4)[int(2)];
    float _S552 = (*_s_dOut_4)[int(3)];
    float _S553 = (*_s_dOut_4)[int(4)];
    float _S554 = (*_s_dOut_4)[int(5)];
    float _S555 = (*_s_dOut_4)[int(6)];
    float _S556 = (*_s_dOut_4)[int(15)] / _S548;
    float _S557 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S547 * _S556));
    float _S558 = (*_s_dOut_4)[int(14)] / _S548;
    float _S559 = (*weights_1)[int(8)] * (_S547 * _S558);
    float _S560 = (*_s_dOut_4)[int(13)] / _S548;
    float _S561 = _S547 * _S560;
    float _S562 = (*weights_1)[int(9)] * ((_S485.primal_0.x + _S485.primal_0.y + _S485.primal_0.z) * 0.3333333432674408f) * - _S556 + (*weights_1)[int(8)] * _S484.primal_0 * - _S558 + (*weights_1)[int(7)] * ((_S483.primal_0.x + _S483.primal_0.y + _S483.primal_0.z) * 0.3333333432674408f) * - _S560;
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = _S526;
    (&_S563)->differential_0 = 0.0f;
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = 9.999999960041972e-13f;
    (&_S564)->differential_0 = 0.0f;
    _d_max_0(&_S563, &_S564, _S562);
    float _S565 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S561);
    float _S566 = _S546 * (*_s_dOut_4)[int(12)];
    float _S567 = _S545 * (_S537 * (*_s_dOut_4)[int(12)]);
    float _S568 = - (_S536 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S542;
    (&_S569)->differential_0 = 0.0f;
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S544;
    (&_S570)->differential_0 = 0.0f;
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S538;
    (&_S571)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S569, &_S570, &_S571, _S568);
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S543;
    (&_S572)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S572, _S570.differential_0);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S539;
    (&_S573)->differential_0 = 0.0f;
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = 9.99999997475242708e-07f;
    (&_S574)->differential_0 = 0.0f;
    _d_max_0(&_S573, &_S574, _S572.differential_0);
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = _S541;
    (&_S575)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S575, _S569.differential_0);
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S540;
    (&_S576)->differential_0 = 0.0f;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = 9.99999997475242708e-07f;
    (&_S577)->differential_0 = 0.0f;
    _d_max_0(&_S576, &_S577, _S575.differential_0);
    float _S578 = _S573.differential_0 + - _S576.differential_0;
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S537;
    (&_S579)->differential_0 = 0.0f;
    DiffPair_float_0 _S580;
    (&_S580)->primal_0 = _S538;
    (&_S580)->differential_0 = 0.0f;
    _d_max_0(&_S579, &_S580, _S578);
    float _S581 = - (_S566 + _S579.differential_0);
    float _S582 = - (_S528 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S533;
    (&_S583)->differential_0 = 0.0f;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S535;
    (&_S584)->differential_0 = 0.0f;
    DiffPair_float_0 _S585;
    (&_S585)->primal_0 = _S529;
    (&_S585)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S583, &_S584, &_S585, _S582);
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = _S534;
    (&_S586)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S586, _S584.differential_0);
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = _S530;
    (&_S587)->differential_0 = 0.0f;
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = 9.99999997475242708e-07f;
    (&_S588)->differential_0 = 0.0f;
    _d_max_0(&_S587, &_S588, _S586.differential_0);
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = _S532;
    (&_S589)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S589, _S583.differential_0);
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = _S531;
    (&_S590)->differential_0 = 0.0f;
    DiffPair_float_0 _S591;
    (&_S591)->primal_0 = 9.99999997475242708e-07f;
    (&_S591)->differential_0 = 0.0f;
    _d_max_0(&_S590, &_S591, _S589.differential_0);
    float _S592 = _S587.differential_0 + - _S590.differential_0;
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = _S526;
    (&_S593)->differential_0 = 0.0f;
    DiffPair_float_0 _S594;
    (&_S594)->primal_0 = _S529;
    (&_S594)->differential_0 = 0.0f;
    _d_max_0(&_S593, &_S594, _S592);
    float _S595 = _S563.differential_0 + _S567 + _S581 + _S593.differential_0;
    DiffPair_float_0 _S596;
    (&_S596)->primal_0 = _S482.primal_0;
    (&_S596)->differential_0 = 0.0f;
    DiffPair_float_0 _S597;
    (&_S597)->primal_0 = 0.0f;
    (&_S597)->differential_0 = 0.0f;
    DiffPair_float_0 _S598;
    (&_S598)->primal_0 = 1.0f;
    (&_S598)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S596, &_S597, &_S598, _S595);
    DiffPair_float_0 _S599 = _S596;
    float _S600 = _S524 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = _S525;
    (&_S601)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S601, _S600);
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = cos_sim_loss_5;
    (&_S602)->differential_0 = 0.0f;
    DiffPair_float_0 _S603;
    (&_S603)->primal_0 = 9.999999960041972e-13f;
    (&_S603)->differential_0 = 0.0f;
    _d_max_0(&_S602, &_S603, _S601.differential_0);
    float _S604 = 0.5f * - (_S600 + _S602.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S605;
    (&_S605)->primal_0 = _S499;
    (&_S605)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S606;
    (&_S606)->primal_0 = _S506;
    (&_S606)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S605, &_S606, _S604);
    float _S607 = _S522 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S608;
    (&_S608)->primal_0 = _S523;
    (&_S608)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S608, _S607);
    DiffPair_float_0 _S609;
    (&_S609)->primal_0 = cos_sim_loss_4;
    (&_S609)->differential_0 = 0.0f;
    DiffPair_float_0 _S610;
    (&_S610)->primal_0 = 9.999999960041972e-13f;
    (&_S610)->differential_0 = 0.0f;
    _d_max_0(&_S609, &_S610, _S608.differential_0);
    float _S611 = 0.5f * - (_S607 + _S609.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S612;
    (&_S612)->primal_0 = _S506;
    (&_S612)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S613;
    (&_S613)->primal_0 = _S513;
    (&_S613)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S612, &_S613, _S611);
    float _S614 = _S520 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = _S521;
    (&_S615)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S615, _S614);
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = cos_sim_loss_3;
    (&_S616)->differential_0 = 0.0f;
    DiffPair_float_0 _S617;
    (&_S617)->primal_0 = 9.999999960041972e-13f;
    (&_S617)->differential_0 = 0.0f;
    _d_max_0(&_S616, &_S617, _S615.differential_0);
    float _S618 = 0.5f * - (_S614 + _S616.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S619;
    (&_S619)->primal_0 = _S499;
    (&_S619)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S620;
    (&_S620)->primal_0 = _S513;
    (&_S620)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S619, &_S620, _S618);
    float3  _S621 = _S613.differential_0 + _S620.differential_0;
    float3  _S622 = _S605.differential_0 + _S619.differential_0;
    float3  _S623 = make_float3 (_S557, _S557, _S557);
    float3  _S624 = make_float3 (_S565, _S565, _S565);
    float3  _S625 = _S606.differential_0 + _S612.differential_0;
    float _S626;
    if(_S516)
    {
        float3  _S627 = _S481.primal_0 * _S621;
        float3  _S628 = _S517 * _S621;
        float _S629 = _S627.x + _S627.y + _S627.z;
        DiffPair_float_0 _S630;
        (&_S630)->primal_0 = _S511;
        (&_S630)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S630, _S629);
        _S626 = _S630.differential_0;
        _S499 = _S628;
    }
    else
    {
        _S626 = 0.0f;
        _S499 = _S486;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S631;
    (&_S631)->primal_0 = _S481.primal_0;
    (&_S631)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
    (&_S632)->primal_0 = _S481.primal_0;
    (&_S632)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S631, &_S632, _S626);
    float3  _S633 = _S632.differential_0 + _S631.differential_0 + _S499;
    if(_S507)
    {
        float3  _S634 = _S480.primal_0 * _S625;
        float3  _S635 = _S508 * _S625;
        float _S636 = _S634.x + _S634.y + _S634.z;
        DiffPair_float_0 _S637;
        (&_S637)->primal_0 = _S504;
        (&_S637)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S637, _S636);
        _S626 = _S637.differential_0;
        _S499 = _S635;
    }
    else
    {
        _S626 = 0.0f;
        _S499 = _S486;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S480.primal_0;
    (&_S638)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S639;
    (&_S639)->primal_0 = _S480.primal_0;
    (&_S639)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S638, &_S639, _S626);
    float3  _S640 = _S639.differential_0 + _S638.differential_0 + _S499;
    if(_S500)
    {
        float3  _S641 = _S479.primal_0 * _S622;
        float3  _S642 = _S501 * _S622;
        float _S643 = _S641.x + _S641.y + _S641.z;
        DiffPair_float_0 _S644;
        (&_S644)->primal_0 = _S497;
        (&_S644)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S644, _S643);
        _S626 = _S644.differential_0;
        _S499 = _S642;
    }
    else
    {
        _S626 = 0.0f;
        _S499 = _S486;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S645;
    (&_S645)->primal_0 = _S479.primal_0;
    (&_S645)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S646;
    (&_S646)->primal_0 = _S479.primal_0;
    (&_S646)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S645, &_S646, _S626);
    float _S647 = _S495 * _S555;
    float _S648 = _S495 * _S554;
    float _S649 = _S493 * _S553;
    float _S650 = _S491 * (_S493 * _S555 + _S648 + _S648 + _S552);
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = _S494;
    (&_S651)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S651, _S650);
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = _S478.primal_0;
    (&_S652)->differential_0 = 0.0f;
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = 0.00009999999747379f;
    (&_S653)->differential_0 = 0.0f;
    _d_max_0(&_S652, &_S653, _S651.differential_0);
    float _S654 = _S491 * (_S647 + _S649 + _S649 + _S551);
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = _S492;
    (&_S655)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S655, _S654);
    DiffPair_float_0 _S656;
    (&_S656)->primal_0 = _S477.primal_0;
    (&_S656)->differential_0 = 0.0f;
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = 0.00009999999747379f;
    (&_S657)->differential_0 = 0.0f;
    _d_max_0(&_S656, &_S657, _S655.differential_0);
    float _S658 = _S487 * _S550;
    DiffPair_float_0 _S659;
    (&_S659)->primal_0 = _S490;
    (&_S659)->differential_0 = 0.0f;
    DiffPair_float_0 _S660;
    (&_S660)->primal_0 = 0.0f;
    (&_S660)->differential_0 = 0.0f;
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = 1.0f;
    (&_S661)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S659, &_S660, &_S661, _S658);
    float _S662 = 0.3333333432674408f * _S659.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S663;
    (&_S663)->primal_0 = _S489;
    (&_S663)->differential_0 = _S486;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S664;
    (&_S664)->primal_0 = _S489;
    (&_S664)->differential_0 = _S486;
    s_bwd_prop_dot_0(&_S663, &_S664, _S662);
    float _S665 = 0.3333333432674408f * (_S488 * _S549);
    float3  _S666 = make_float3 (_S665, _S665, _S665);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S667;
    (&_S667)->primal_0 = _S489;
    (&_S667)->differential_0 = _S486;
    s_bwd_prop_abs_0(&_S667, _S666);
    float3  _S668 = _S664.differential_0 + _S663.differential_0 + _S667.differential_0;
    float3  _S669 = - _S668;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S623;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S559;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S624;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S599.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S633;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S640;
    float3  _S670 = _S646.differential_0 + _S645.differential_0 + _S499;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S670;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S652.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S656.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S668;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S669;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S671, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S672, DiffPair_float_0 * _S673, DiffPair_float_0 * _S674, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S675, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S676, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S677, DiffPair_float_0 * _S678, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S679, DiffPair_float_0 * _S680, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S681, bool _S682, bool _S683, bool _S684, bool _S685, bool _S686, FixedArray<float, 10>  * _S687, FixedArray<float, 23>  * _S688)
{
    s_bwd_prop_per_pixel_losses_0(_S671, _S672, _S673, _S674, _S675, _S676, _S677, _S678, _S679, _S680, _S681, _S682, _S683, _S684, _S685, _S686, _S687, _S688);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S689 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S689;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S689;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S689;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S689;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S689;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S689;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S689;
    FixedArray<float, 10>  _S690 = weights_2;
    FixedArray<float, 23>  _S691 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S690, &_S691);
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
    float _S692 = 1.0f / ((*dpx_14).primal_0 * 2.30258512496948242f) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S692;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 10>  weights_3, FixedArray<float, 10>  * _S693)
{
    FixedArray<float, 10>  losses_2;
    float _S694 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_2[int(0)] = raw_losses_0[int(0)] / _S694;
    losses_2[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S694)));
    bool _S695;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S695 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S695 = false;
    }
    float _S696;
    if(_S695)
    {
        _S696 = weights_3[int(1)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S696 = 0.0f;
    }
    losses_2[int(2)] = _S696;
    losses_2[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_2[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S697 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_2[int(6)] = raw_losses_0[int(12)] / _S697;
    losses_2[int(7)] = raw_losses_0[int(13)] / _S697;
    losses_2[int(8)] = raw_losses_0[int(14)] / _S697;
    losses_2[int(9)] = raw_losses_0[int(15)] / _S697;
    *_S693 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S698)
{
    return (F32_sqrt((_S698)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S699, float _S700)
{
    _d_log10_0(_S699, _S700);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S701 = dpraw_losses_0->primal_0;
    float _S702 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S703 = _S702 * _S702;
    float _S704 = dpraw_losses_0->primal_0[int(1)] / _S702;
    bool _S705 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S706;
    if(_S705)
    {
        _S706 = (_S701[int(3)]) != 0.0f;
    }
    else
    {
        _S706 = false;
    }
    float _S707;
    float _S708;
    float _S709;
    float _S710;
    float _S711;
    float _S712;
    float _S713;
    float _S714;
    float _S715;
    float _S716;
    float _S717;
    float _S718;
    float _S719;
    float _S720;
    float _S721;
    if(_S706)
    {
        float _S722 = _S701[int(2)] * _S701[int(3)];
        float _S723 = _S701[int(18)] * _S701[int(18)];
        float _S724 = _S701[int(6)] - _S722 / _S701[int(18)];
        float _S725 = _S701[int(2)] * _S701[int(2)];
        float _S726 = _S701[int(4)] - _S725 / _S701[int(18)];
        float _S727 = _S701[int(3)] * _S701[int(3)];
        float _S728 = _S701[int(5)] - _S727 / _S701[int(18)];
        float _S729 = _S726 * _S728 + 1.0f;
        float _S730 = (F32_max((9.999999960041972e-13f), (_S729)));
        float _S731 = s_primal_ctx_sqrt_0(_S730);
        float _S732 = _S731 * _S731;
        float _S733 = 1.0f - _S724 / _S731;
        _S707 = (*weights_4)[int(1)];
        _S708 = _S733;
        _S709 = _S732;
        _S710 = _S724;
        _S711 = _S731;
        _S712 = _S730;
        _S713 = _S729;
        _S714 = _S726;
        _S715 = _S728;
        _S716 = _S723;
        _S717 = _S727;
        _S718 = _S701[int(3)];
        _S719 = _S725;
        _S720 = _S701[int(2)];
        _S721 = _S722;
    }
    else
    {
        _S707 = 0.0f;
        _S708 = 0.0f;
        _S709 = 0.0f;
        _S710 = 0.0f;
        _S711 = 0.0f;
        _S712 = 0.0f;
        _S713 = 0.0f;
        _S714 = 0.0f;
        _S715 = 0.0f;
        _S716 = 0.0f;
        _S717 = 0.0f;
        _S718 = 0.0f;
        _S719 = 0.0f;
        _S720 = 0.0f;
        _S721 = 0.0f;
    }
    float _S734 = (F32_max((_S701[int(19)]), (1.0f)));
    float _S735 = _S734 * _S734;
    float _S736 = (F32_max((_S701[int(20)]), (1.0f)));
    float _S737 = _S736 * _S736;
    float _S738 = float((I32_max((int((_S701[int(19)]) > 0.5f) + int((_S701[int(20)]) > 0.5f)), (int(1)))));
    float _S739 = _S701[int(9)] + _S701[int(10)];
    float _S740 = (F32_max((_S701[int(22)]), (1.0f)));
    float _S741 = _S740 * _S740;
    float _S742 = (F32_max((_S701[int(21)]), (1.0f)));
    float _S743 = _S742 * _S742;
    float _S744 = (F32_max((_S701[int(16)]), (1.0f)));
    float _S745 = _S744 * _S744;
    float _S746 = (*_s_dOut_5)[int(0)];
    float _S747 = (*_s_dOut_5)[int(1)];
    float _S748 = (*_s_dOut_5)[int(2)];
    float _S749 = (*_s_dOut_5)[int(9)] / _S745;
    float _S750 = _S744 * _S749;
    float _S751 = (*_s_dOut_5)[int(8)] / _S745;
    float _S752 = _S744 * _S751;
    float _S753 = (*_s_dOut_5)[int(7)] / _S745;
    float _S754 = _S744 * _S753;
    float _S755 = (*_s_dOut_5)[int(6)] / _S745;
    float _S756 = _S744 * _S755;
    float _S757 = _S701[int(15)] * - _S749 + _S701[int(14)] * - _S751 + _S701[int(13)] * - _S753 + _S701[int(12)] * - _S755;
    DiffPair_float_0 _S758;
    (&_S758)->primal_0 = _S701[int(16)];
    (&_S758)->differential_0 = 0.0f;
    DiffPair_float_0 _S759;
    (&_S759)->primal_0 = 1.0f;
    (&_S759)->differential_0 = 0.0f;
    _d_max_0(&_S758, &_S759, _S757);
    float _S760 = (*_s_dOut_5)[int(5)] / _S743;
    float _S761 = _S701[int(11)] * - _S760;
    float _S762 = _S742 * _S760;
    DiffPair_float_0 _S763;
    (&_S763)->primal_0 = _S701[int(21)];
    (&_S763)->differential_0 = 0.0f;
    DiffPair_float_0 _S764;
    (&_S764)->primal_0 = 1.0f;
    (&_S764)->differential_0 = 0.0f;
    _d_max_0(&_S763, &_S764, _S761);
    float _S765 = (*_s_dOut_5)[int(4)] / _S741;
    float _S766 = _S739 * - _S765;
    float _S767 = _S740 * _S765;
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = _S701[int(22)];
    (&_S768)->differential_0 = 0.0f;
    DiffPair_float_0 _S769;
    (&_S769)->primal_0 = 1.0f;
    (&_S769)->differential_0 = 0.0f;
    _d_max_0(&_S768, &_S769, _S766);
    float _S770 = (*_s_dOut_5)[int(3)] / _S738;
    float _S771 = _S770 / _S737;
    float _S772 = _S701[int(8)] * - _S771;
    float _S773 = _S736 * _S771;
    DiffPair_float_0 _S774;
    (&_S774)->primal_0 = _S701[int(20)];
    (&_S774)->differential_0 = 0.0f;
    DiffPair_float_0 _S775;
    (&_S775)->primal_0 = 1.0f;
    (&_S775)->differential_0 = 0.0f;
    _d_max_0(&_S774, &_S775, _S772);
    float _S776 = _S770 / _S735;
    float _S777 = _S701[int(7)] * - _S776;
    float _S778 = _S734 * _S776;
    DiffPair_float_0 _S779;
    (&_S779)->primal_0 = _S701[int(19)];
    (&_S779)->differential_0 = 0.0f;
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = 1.0f;
    (&_S780)->differential_0 = 0.0f;
    _d_max_0(&_S779, &_S780, _S777);
    FixedArray<float, 23>  _S781;
    _S781[int(0)] = 0.0f;
    _S781[int(1)] = 0.0f;
    _S781[int(2)] = 0.0f;
    _S781[int(3)] = 0.0f;
    _S781[int(4)] = 0.0f;
    _S781[int(5)] = 0.0f;
    _S781[int(6)] = 0.0f;
    _S781[int(7)] = 0.0f;
    _S781[int(8)] = 0.0f;
    _S781[int(9)] = 0.0f;
    _S781[int(10)] = 0.0f;
    _S781[int(11)] = 0.0f;
    _S781[int(12)] = 0.0f;
    _S781[int(13)] = 0.0f;
    _S781[int(14)] = 0.0f;
    _S781[int(15)] = 0.0f;
    _S781[int(16)] = 0.0f;
    _S781[int(17)] = 0.0f;
    _S781[int(18)] = 0.0f;
    _S781[int(19)] = 0.0f;
    _S781[int(20)] = 0.0f;
    _S781[int(21)] = 0.0f;
    _S781[int(22)] = 0.0f;
    _S781[int(15)] = _S750;
    _S781[int(14)] = _S752;
    _S781[int(13)] = _S754;
    _S781[int(16)] = _S758.differential_0;
    _S781[int(12)] = _S756;
    _S781[int(21)] = _S763.differential_0;
    _S781[int(11)] = _S762;
    _S781[int(22)] = _S768.differential_0;
    _S781[int(10)] = _S767;
    _S781[int(9)] = _S767;
    _S781[int(20)] = _S774.differential_0;
    _S781[int(8)] = _S773;
    _S781[int(19)] = _S779.differential_0;
    _S781[int(7)] = _S778;
    float _S782 = _S781[int(0)];
    float _S783 = _S781[int(1)];
    float _S784 = _S781[int(2)];
    float _S785 = _S781[int(3)];
    float _S786 = _S781[int(4)];
    float _S787 = _S781[int(5)];
    float _S788 = _S781[int(6)];
    float _S789 = _S781[int(7)];
    float _S790 = _S781[int(8)];
    float _S791 = _S781[int(9)];
    float _S792 = _S781[int(10)];
    float _S793 = _S781[int(11)];
    float _S794 = _S781[int(12)];
    float _S795 = _S781[int(13)];
    float _S796 = _S781[int(14)];
    float _S797 = _S781[int(15)];
    float _S798 = _S781[int(16)];
    float _S799 = _S781[int(17)];
    float _S800 = _S781[int(18)];
    float _S801 = _S781[int(19)];
    float _S802 = _S781[int(20)];
    float _S803 = _S781[int(21)];
    float _S804 = _S781[int(22)];
    FixedArray<float, 23>  _S805;
    if(_S706)
    {
        float _S806 = _S707 * _S748;
        DiffPair_float_0 _S807;
        (&_S807)->primal_0 = _S708;
        (&_S807)->differential_0 = 0.0f;
        DiffPair_float_0 _S808;
        (&_S808)->primal_0 = 0.0f;
        (&_S808)->differential_0 = 0.0f;
        DiffPair_float_0 _S809;
        (&_S809)->primal_0 = 2.0f;
        (&_S809)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S807, &_S808, &_S809, _S806);
        float _S810 = - _S807.differential_0 / _S709;
        float _S811 = _S710 * - _S810;
        float _S812 = _S711 * _S810;
        DiffPair_float_0 _S813;
        (&_S813)->primal_0 = _S712;
        (&_S813)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S813, _S811);
        DiffPair_float_0 _S814;
        (&_S814)->primal_0 = 9.999999960041972e-13f;
        (&_S814)->differential_0 = 0.0f;
        DiffPair_float_0 _S815;
        (&_S815)->primal_0 = _S713;
        (&_S815)->differential_0 = 0.0f;
        _d_max_0(&_S814, &_S815, _S813.differential_0);
        float _S816 = _S714 * _S815.differential_0;
        float _S817 = _S715 * _S815.differential_0;
        float _S818 = - _S816 / _S716;
        float _S819 = _S718 * (_S701[int(18)] * _S818);
        float _S820 = - _S817 / _S716;
        float _S821 = _S720 * (_S701[int(18)] * _S820);
        float _S822 = - _S812 / _S716;
        float _S823 = _S701[int(18)] * _S822;
        float _S824 = _S819 + _S819 + _S720 * _S823;
        float _S825 = _S821 + _S821 + _S718 * _S823;
        float _S826 = _S717 * - _S818 + _S719 * - _S820 + _S721 * - _S822;
        FixedArray<float, 23>  _S827;
        _S827[int(0)] = 0.0f;
        _S827[int(1)] = 0.0f;
        _S827[int(2)] = 0.0f;
        _S827[int(3)] = 0.0f;
        _S827[int(4)] = 0.0f;
        _S827[int(5)] = 0.0f;
        _S827[int(6)] = 0.0f;
        _S827[int(7)] = 0.0f;
        _S827[int(8)] = 0.0f;
        _S827[int(9)] = 0.0f;
        _S827[int(10)] = 0.0f;
        _S827[int(11)] = 0.0f;
        _S827[int(12)] = 0.0f;
        _S827[int(13)] = 0.0f;
        _S827[int(14)] = 0.0f;
        _S827[int(15)] = 0.0f;
        _S827[int(16)] = 0.0f;
        _S827[int(17)] = 0.0f;
        _S827[int(18)] = 0.0f;
        _S827[int(19)] = 0.0f;
        _S827[int(20)] = 0.0f;
        _S827[int(21)] = 0.0f;
        _S827[int(22)] = 0.0f;
        _S827[int(5)] = _S816;
        _S827[int(4)] = _S817;
        _S827[int(3)] = _S824;
        _S827[int(2)] = _S825;
        _S827[int(6)] = _S812;
        float _S828 = _S783 + _S827[int(1)];
        float _S829 = _S784 + _S827[int(2)];
        float _S830 = _S785 + _S827[int(3)];
        float _S831 = _S786 + _S827[int(4)];
        float _S832 = _S787 + _S827[int(5)];
        float _S833 = _S788 + _S827[int(6)];
        float _S834 = _S789 + _S827[int(7)];
        float _S835 = _S790 + _S827[int(8)];
        float _S836 = _S791 + _S827[int(9)];
        float _S837 = _S792 + _S827[int(10)];
        float _S838 = _S793 + _S827[int(11)];
        float _S839 = _S794 + _S827[int(12)];
        float _S840 = _S795 + _S827[int(13)];
        float _S841 = _S796 + _S827[int(14)];
        float _S842 = _S797 + _S827[int(15)];
        float _S843 = _S798 + _S827[int(16)];
        float _S844 = _S799 + _S827[int(17)];
        float _S845 = _S800 + _S827[int(18)];
        float _S846 = _S801 + _S827[int(19)];
        float _S847 = _S802 + _S827[int(20)];
        float _S848 = _S803 + _S827[int(21)];
        float _S849 = _S804 + _S827[int(22)];
        _S805[int(0)] = _S782 + _S827[int(0)];
        _S805[int(1)] = _S828;
        _S805[int(2)] = _S829;
        _S805[int(3)] = _S830;
        _S805[int(4)] = _S831;
        _S805[int(5)] = _S832;
        _S805[int(6)] = _S833;
        _S805[int(7)] = _S834;
        _S805[int(8)] = _S835;
        _S805[int(9)] = _S836;
        _S805[int(10)] = _S837;
        _S805[int(11)] = _S838;
        _S805[int(12)] = _S839;
        _S805[int(13)] = _S840;
        _S805[int(14)] = _S841;
        _S805[int(15)] = _S842;
        _S805[int(16)] = _S843;
        _S805[int(17)] = _S844;
        _S805[int(18)] = _S845;
        _S805[int(19)] = _S846;
        _S805[int(20)] = _S847;
        _S805[int(21)] = _S848;
        _S805[int(22)] = _S849;
        _S707 = _S826;
    }
    else
    {
        _S805[int(0)] = _S782;
        _S805[int(1)] = _S783;
        _S805[int(2)] = _S784;
        _S805[int(3)] = _S785;
        _S805[int(4)] = _S786;
        _S805[int(5)] = _S787;
        _S805[int(6)] = _S788;
        _S805[int(7)] = _S789;
        _S805[int(8)] = _S790;
        _S805[int(9)] = _S791;
        _S805[int(10)] = _S792;
        _S805[int(11)] = _S793;
        _S805[int(12)] = _S794;
        _S805[int(13)] = _S795;
        _S805[int(14)] = _S796;
        _S805[int(15)] = _S797;
        _S805[int(16)] = _S798;
        _S805[int(17)] = _S799;
        _S805[int(18)] = _S800;
        _S805[int(19)] = _S801;
        _S805[int(20)] = _S802;
        _S805[int(21)] = _S803;
        _S805[int(22)] = _S804;
        _S707 = 0.0f;
    }
    if(_S705)
    {
        FixedArray<float, 23>  _S850;
        _S850[int(0)] = 0.0f;
        _S850[int(1)] = 0.0f;
        _S850[int(2)] = 0.0f;
        _S850[int(3)] = 0.0f;
        _S850[int(4)] = 0.0f;
        _S850[int(5)] = 0.0f;
        _S850[int(6)] = 0.0f;
        _S850[int(7)] = 0.0f;
        _S850[int(8)] = 0.0f;
        _S850[int(9)] = 0.0f;
        _S850[int(10)] = 0.0f;
        _S850[int(11)] = 0.0f;
        _S850[int(12)] = 0.0f;
        _S850[int(13)] = 0.0f;
        _S850[int(14)] = 0.0f;
        _S850[int(15)] = 0.0f;
        _S850[int(16)] = 0.0f;
        _S850[int(17)] = 0.0f;
        _S850[int(18)] = 0.0f;
        _S850[int(19)] = 0.0f;
        _S850[int(20)] = 0.0f;
        _S850[int(21)] = 0.0f;
        _S850[int(22)] = 0.0f;
        _S850[int(3)] = 0.0f;
        float _S851 = _S805[int(1)] + _S850[int(1)];
        float _S852 = _S805[int(2)] + _S850[int(2)];
        float _S853 = _S805[int(3)] + _S850[int(3)];
        float _S854 = _S805[int(4)] + _S850[int(4)];
        float _S855 = _S805[int(5)] + _S850[int(5)];
        float _S856 = _S805[int(6)] + _S850[int(6)];
        float _S857 = _S805[int(7)] + _S850[int(7)];
        float _S858 = _S805[int(8)] + _S850[int(8)];
        float _S859 = _S805[int(9)] + _S850[int(9)];
        float _S860 = _S805[int(10)] + _S850[int(10)];
        float _S861 = _S805[int(11)] + _S850[int(11)];
        float _S862 = _S805[int(12)] + _S850[int(12)];
        float _S863 = _S805[int(13)] + _S850[int(13)];
        float _S864 = _S805[int(14)] + _S850[int(14)];
        float _S865 = _S805[int(15)] + _S850[int(15)];
        float _S866 = _S805[int(16)] + _S850[int(16)];
        float _S867 = _S805[int(17)] + _S850[int(17)];
        float _S868 = _S805[int(18)] + _S850[int(18)];
        float _S869 = _S805[int(19)] + _S850[int(19)];
        float _S870 = _S805[int(20)] + _S850[int(20)];
        float _S871 = _S805[int(21)] + _S850[int(21)];
        float _S872 = _S805[int(22)] + _S850[int(22)];
        _S805[int(0)] = _S805[int(0)] + _S850[int(0)];
        _S805[int(1)] = _S851;
        _S805[int(2)] = _S852;
        _S805[int(3)] = _S853;
        _S805[int(4)] = _S854;
        _S805[int(5)] = _S855;
        _S805[int(6)] = _S856;
        _S805[int(7)] = _S857;
        _S805[int(8)] = _S858;
        _S805[int(9)] = _S859;
        _S805[int(10)] = _S860;
        _S805[int(11)] = _S861;
        _S805[int(12)] = _S862;
        _S805[int(13)] = _S863;
        _S805[int(14)] = _S864;
        _S805[int(15)] = _S865;
        _S805[int(16)] = _S866;
        _S805[int(17)] = _S867;
        _S805[int(18)] = _S868;
        _S805[int(19)] = _S869;
        _S805[int(20)] = _S870;
        _S805[int(21)] = _S871;
        _S805[int(22)] = _S872;
    }
    float _S873 = -10.0f * _S747;
    DiffPair_float_0 _S874;
    (&_S874)->primal_0 = _S704;
    (&_S874)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S874, _S873);
    float _S875 = _S874.differential_0 / _S703;
    float _S876 = _S702 * _S875;
    float _S877 = _S746 / _S703;
    float _S878 = _S702 * _S877;
    float _S879 = _S701[int(1)] * - _S875 + _S701[int(0)] * - _S877;
    DiffPair_float_0 _S880;
    (&_S880)->primal_0 = _S701[int(17)];
    (&_S880)->differential_0 = 0.0f;
    DiffPair_float_0 _S881;
    (&_S881)->primal_0 = 1.0f;
    (&_S881)->differential_0 = 0.0f;
    _d_max_0(&_S880, &_S881, _S879);
    FixedArray<float, 23>  _S882;
    _S882[int(0)] = 0.0f;
    _S882[int(1)] = 0.0f;
    _S882[int(2)] = 0.0f;
    _S882[int(3)] = 0.0f;
    _S882[int(4)] = 0.0f;
    _S882[int(5)] = 0.0f;
    _S882[int(6)] = 0.0f;
    _S882[int(7)] = 0.0f;
    _S882[int(8)] = 0.0f;
    _S882[int(9)] = 0.0f;
    _S882[int(10)] = 0.0f;
    _S882[int(11)] = 0.0f;
    _S882[int(12)] = 0.0f;
    _S882[int(13)] = 0.0f;
    _S882[int(14)] = 0.0f;
    _S882[int(15)] = 0.0f;
    _S882[int(16)] = 0.0f;
    _S882[int(17)] = 0.0f;
    _S882[int(18)] = 0.0f;
    _S882[int(19)] = 0.0f;
    _S882[int(20)] = 0.0f;
    _S882[int(21)] = 0.0f;
    _S882[int(22)] = 0.0f;
    _S882[int(18)] = _S707;
    _S882[int(1)] = _S876;
    _S882[int(17)] = _S880.differential_0;
    _S882[int(0)] = _S878;
    FixedArray<float, 23>  _S883 = {
        _S805[int(0)] + _S882[int(0)], _S805[int(1)] + _S882[int(1)], _S805[int(2)] + _S882[int(2)], _S805[int(3)] + _S882[int(3)], _S805[int(4)] + _S882[int(4)], _S805[int(5)] + _S882[int(5)], _S805[int(6)] + _S882[int(6)], _S805[int(7)] + _S882[int(7)], _S805[int(8)] + _S882[int(8)], _S805[int(9)] + _S882[int(9)], _S805[int(10)] + _S882[int(10)], _S805[int(11)] + _S882[int(11)], _S805[int(12)] + _S882[int(12)], _S805[int(13)] + _S882[int(13)], _S805[int(14)] + _S882[int(14)], _S805[int(15)] + _S882[int(15)], _S805[int(16)] + _S882[int(16)], _S805[int(17)] + _S882[int(17)], _S805[int(18)] + _S882[int(18)], _S805[int(19)] + _S882[int(19)], _S805[int(20)] + _S882[int(20)], _S805[int(21)] + _S882[int(21)], _S805[int(22)] + _S882[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S883;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S884, FixedArray<float, 10>  * _S885, FixedArray<float, 10>  * _S886)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S884, _S885, _S886);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 10>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S887)
{
    FixedArray<float, 23>  _S888 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S888;
    FixedArray<float, 10>  _S889 = weights_5;
    FixedArray<float, 10>  _S890 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S889, &_S890);
    *_S887 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_25, float3  y_8)
{
    float3  result_21;
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
        *_slang_vector_get_element_ptr(&result_21, i_12) = (F32_min((_slang_vector_get_element(x_25, i_12)), (_slang_vector_get_element(y_8, i_12))));
        i_12 = i_12 + int(1);
    }
    return result_21;
}

inline __device__ float3  max_0(float3  x_26, float3  y_9)
{
    float3  result_22;
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
        *_slang_vector_get_element_ptr(&result_22, i_13) = (F32_max((_slang_vector_get_element(x_26, i_13)), (_slang_vector_get_element(y_9, i_13))));
        i_13 = i_13 + int(1);
    }
    return result_22;
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

inline __device__ float3  clamp_1(float3  x_27, float3  minBound_1, float3  maxBound_1)
{
    return min_0(max_0(x_27, minBound_1), maxBound_1);
}

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_1(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S891, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S892, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S893, float3  _S894)
{
    _d_clamp_vector_0(_S891, _S892, _S893, _S894);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S895 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S896 = make_float3 (_S895);
    float3  _S897 = make_float3 (0.0f);
    float3  _S898 = make_float3 (1.0f);
    float3  _S899 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S900;
    (&_S900)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S895) * (*dpbackground_0).primal_0;
    (&_S900)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S901;
    (&_S901)->primal_0 = _S897;
    (&_S901)->differential_0 = _S899;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S902;
    (&_S902)->primal_0 = _S898;
    (&_S902)->differential_0 = _S899;
    s_bwd_prop_clamp_1(&_S900, &_S901, &_S902, _s_dOut_6);
    float3  _S903 = _S896 * _S900.differential_0;
    float3  _S904 = (*dpbackground_0).primal_0 * _S900.differential_0;
    float _S905 = - (_S904.x + _S904.y + _S904.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S903;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S905;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S900.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S906, DiffPair_float_0 * _S907, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S908, float3  _S909)
{
    s_bwd_prop_blend_background_0(_S906, _S907, _S908, _S909);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S910 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S910;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S910;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
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
        DiffPair_float_0 _S911 = *dpx_16;
        float _S912 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_20;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S912;
        float _S913 = val_0 * (F32_log((_S911.primal_0))) * dOut_20;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S913;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S914 = rgb_2;
    float _S915;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S915 = _S914.x * 12.92000007629394531f;
    }
    else
    {
        _S915 = 1.0549999475479126f * (F32_pow((_S914.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S914)->x) = _S915;
    if((_S914.y) < 0.00313080009073019f)
    {
        _S915 = _S914.y * 12.92000007629394531f;
    }
    else
    {
        _S915 = 1.0549999475479126f * (F32_pow((_S914.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S914)->y) = _S915;
    if((_S914.z) < 0.00313080009073019f)
    {
        _S915 = _S914.z * 12.92000007629394531f;
    }
    else
    {
        _S915 = 1.0549999475479126f * (F32_pow((_S914.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S914)->z) = _S915;
    return _S914;
}

inline __device__ float s_primal_ctx_pow_0(float _S916, float _S917)
{
    return (F32_pow((_S916), (_S917)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S918, DiffPair_float_0 * _S919, float _S920)
{
    _d_pow_0(_S918, _S919, _S920);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_7)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S921 = *dprgb_1;
    float _S922 = (*dprgb_1).primal_0.x;
    bool _S923 = _S922 < 0.00313080009073019f;
    float _S924;
    if(_S923)
    {
        _S924 = _S922 * 12.92000007629394531f;
    }
    else
    {
        _S924 = 1.0549999475479126f * s_primal_ctx_pow_0(_S922, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S925 = _S921.primal_0;
    *&((&_S925)->x) = _S924;
    float _S926 = _S925.y;
    bool _S927 = _S926 < 0.00313080009073019f;
    if(_S927)
    {
        _S924 = _S926 * 12.92000007629394531f;
    }
    else
    {
        _S924 = 1.0549999475479126f * s_primal_ctx_pow_0(_S926, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S925)->y) = _S924;
    float _S928 = _S925.z;
    bool _S929 = _S928 < 0.00313080009073019f;
    _S925 = _s_dOut_7;
    *&((&_S925)->z) = 0.0f;
    if(_S929)
    {
        _S924 = 12.92000007629394531f * _s_dOut_7.z;
    }
    else
    {
        float _S930 = 1.0549999475479126f * _s_dOut_7.z;
        DiffPair_float_0 _S931;
        (&_S931)->primal_0 = _S928;
        (&_S931)->differential_0 = 0.0f;
        DiffPair_float_0 _S932;
        (&_S932)->primal_0 = 0.4166666567325592f;
        (&_S932)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S931, &_S932, _S930);
        _S924 = _S931.differential_0;
    }
    float3  _S933 = _S925 + make_float3 (0.0f, 0.0f, _S924);
    _S925 = _S933;
    *&((&_S925)->y) = 0.0f;
    if(_S927)
    {
        _S924 = 12.92000007629394531f * _S933.y;
    }
    else
    {
        float _S934 = 1.0549999475479126f * _S933.y;
        DiffPair_float_0 _S935;
        (&_S935)->primal_0 = _S926;
        (&_S935)->differential_0 = 0.0f;
        DiffPair_float_0 _S936;
        (&_S936)->primal_0 = 0.4166666567325592f;
        (&_S936)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S935, &_S936, _S934);
        _S924 = _S935.differential_0;
    }
    float3  _S937 = _S925 + make_float3 (0.0f, _S924, 0.0f);
    _S925 = _S937;
    *&((&_S925)->x) = 0.0f;
    if(_S923)
    {
        _S924 = 12.92000007629394531f * _S937.x;
    }
    else
    {
        float _S938 = 1.0549999475479126f * _S937.x;
        DiffPair_float_0 _S939;
        (&_S939)->primal_0 = _S922;
        (&_S939)->differential_0 = 0.0f;
        DiffPair_float_0 _S940;
        (&_S940)->primal_0 = 0.4166666567325592f;
        (&_S940)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S939, &_S940, _S938);
        _S924 = _S939.differential_0;
    }
    float3  _S941 = _S925 + make_float3 (_S924, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S941;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S942, float3  _S943)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S942, _S943);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S944 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S944;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_7, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S945 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = _S945;
    FixedArray<float, 10>  _S946 = dist_coeffs_7;
    bool _S947 = undistort_point_0(_S945, &_S946, int(12), &uv_7);
    if(!_S947)
    {
        int3  _S948 = make_int3 (int(0));
        float3  _S949 = make_float3 ((float)_S948.x, (float)_S948.y, (float)_S948.z);
        return _S949;
    }
    float3  raydir_6;
    if(is_fisheye_4)
    {
        float theta_4 = length_1(uv_7);
        float3  raydir_7 = make_float3 ((uv_7 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_7 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
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
        float3  raydir_8 = make_float3 (uv_7.x, uv_7.y, 1.0f);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S950, float3  _S951)
{
    return cross_0(_S950, _S951);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_8)
{
    float _S952 = (*dpx_17).primal_0.x;
    float _S953 = (*dpx_17).primal_0.y;
    float _S954 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S955;
    (&_S955)->primal_0 = _S952 * _S952 + _S953 * _S953 + _S954 * _S954;
    (&_S955)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S955, _s_dOut_8);
    float _S956 = (*dpx_17).primal_0.z * _S955.differential_0;
    float _S957 = _S956 + _S956;
    float _S958 = (*dpx_17).primal_0.y * _S955.differential_0;
    float _S959 = _S958 + _S958;
    float _S960 = (*dpx_17).primal_0.x * _S955.differential_0;
    float _S961 = _S960 + _S960;
    float3  _S962 = make_float3 (0.0f);
    *&((&_S962)->z) = _S957;
    *&((&_S962)->y) = _S959;
    *&((&_S962)->x) = _S961;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S962;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S963, float _S964)
{
    s_bwd_prop_length_impl_1(_S963, _S964);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S965, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S966, float3  _S967)
{
    _d_cross_0(_S965, _S966, _S967);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S968 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S969 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S970 = s_primal_ctx_cross_0(dx_0, _S969);
    bool _S971 = (s_primal_ctx_dot_0(_S970, _S970)) != 0.0f;
    float3  _S972;
    float3  _S973;
    if(_S971)
    {
        float _S974 = length_2(_S970);
        float3  _S975 = make_float3 (_S974);
        _S972 = make_float3 (_S974 * _S974);
        _S973 = _S975;
    }
    else
    {
        _S972 = _S968;
        _S973 = _S968;
    }
    if(_S971)
    {
        float3  _S976 = _s_dOut_9 / _S972;
        float3  _S977 = _S970 * - _S976;
        float3  _S978 = _S973 * _S976;
        float _S979 = _S977.x + _S977.y + _S977.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S980;
        (&_S980)->primal_0 = _S970;
        (&_S980)->differential_0 = _S968;
        s_bwd_length_impl_1(&_S980, _S979);
        _S972 = _S978 + _S980.differential_0;
    }
    else
    {
        _S972 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S981;
    (&_S981)->primal_0 = _S970;
    (&_S981)->differential_0 = _S968;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S982;
    (&_S982)->primal_0 = _S970;
    (&_S982)->differential_0 = _S968;
    s_bwd_prop_dot_0(&_S981, &_S982, 0.0f);
    float3  _S983 = _S982.differential_0 + _S981.differential_0 + _S972;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S984;
    (&_S984)->primal_0 = dx_0;
    (&_S984)->differential_0 = _S968;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S985;
    (&_S985)->primal_0 = _S969;
    (&_S985)->differential_0 = _S968;
    s_bwd_prop_cross_0(&_S984, &_S985, _S983);
    float3  s_diff_dy_T_0 = - _S985.differential_0;
    float3  _S986 = - s_diff_dy_T_0;
    float3  _S987 = - _S984.differential_0;
    FixedArray<float3 , 4>  _S988;
    _S988[int(0)] = _S968;
    _S988[int(1)] = _S968;
    _S988[int(2)] = _S968;
    _S988[int(3)] = _S968;
    _S988[int(2)] = _S986;
    _S988[int(3)] = s_diff_dy_T_0;
    _S988[int(0)] = _S987;
    _S988[int(1)] = _S984.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S988;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S989, float3  _S990)
{
    s_bwd_prop_points_to_normal_0(_S989, _S990);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S991 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S991;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_8, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S992 = float2 {intrins_1.z, intrins_1.w};
    float2  _S993 = float2 {intrins_1.x, intrins_1.y};
    float2  _S994 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S992) / _S993;
    float2  uv_8 = _S994;
    FixedArray<float, 10>  _S995 = dist_coeffs_8;
    bool _S996 = undistort_point_0(_S994, &_S995, int(12), &uv_8);
    if(!_S996)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_9;
    if(is_fisheye_5)
    {
        float theta_5 = length_1(uv_8);
        float3  raydir_10 = make_float3 ((uv_8 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_8 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
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
        float3  raydir_11 = make_float3 (uv_8.x, uv_8.y, 1.0f);
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
    float2  _S997 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S992) / _S993;
    float2  uv_9 = _S997;
    FixedArray<float, 10>  _S998 = dist_coeffs_8;
    bool _S999 = undistort_point_0(_S997, &_S998, int(12), &uv_9);
    if(!_S999)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_6 = length_1(uv_9);
        float3  raydir_12 = make_float3 ((uv_9 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (uv_9 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
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
        float3  raydir_13 = make_float3 (uv_9.x, uv_9.y, 1.0f);
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
    float2  _S1000 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S992) / _S993;
    float2  uv_10 = _S1000;
    FixedArray<float, 10>  _S1001 = dist_coeffs_8;
    bool _S1002 = undistort_point_0(_S1000, &_S1001, int(12), &uv_10);
    if(!_S1002)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_7 = length_1(uv_10);
        float3  raydir_14 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (uv_10 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
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
        float3  raydir_15 = make_float3 (uv_10.x, uv_10.y, 1.0f);
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
    float2  _S1003 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S992) / _S993;
    float2  uv_11 = _S1003;
    FixedArray<float, 10>  _S1004 = dist_coeffs_8;
    bool _S1005 = undistort_point_0(_S1003, &_S1004, int(12), &uv_11);
    if(!_S1005)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_8 = length_1(uv_11);
        float3  raydir_16 = make_float3 ((uv_11 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).x, (uv_11 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).y, (F32_cos((theta_8))));
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
        float3  raydir_17 = make_float3 (uv_11.x, uv_11.y, 1.0f);
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
    float2  _S1006;
    bool _S1007;
    float2  _S1008;
    bool _S1009;
    float2  _S1010;
    bool _S1011;
    float2  _S1012;
    bool _S1013;
};

inline __device__ float s_primal_ctx_sin_0(float _S1014)
{
    return (F32_sin((_S1014)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1015)
{
    return (F32_cos((_S1015)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1016 = make_float2 (0.0f);
    _s_diff_ctx_0->_S1006 = _S1016;
    _s_diff_ctx_0->_S1007 = false;
    _s_diff_ctx_0->_S1008 = _S1016;
    _s_diff_ctx_0->_S1009 = false;
    _s_diff_ctx_0->_S1010 = _S1016;
    _s_diff_ctx_0->_S1011 = false;
    _s_diff_ctx_0->_S1012 = _S1016;
    _s_diff_ctx_0->_S1013 = false;
    _s_diff_ctx_0->_S1008 = _S1016;
    _s_diff_ctx_0->_S1009 = false;
    _s_diff_ctx_0->_S1010 = _S1016;
    _s_diff_ctx_0->_S1011 = false;
    _s_diff_ctx_0->_S1012 = _S1016;
    _s_diff_ctx_0->_S1013 = false;
    float3  _S1017 = make_float3 (0.0f);
    float2  _S1018 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1019 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1020 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1018) / _S1019;
    float2  _S1021 = _S1020;
    bool _S1022 = undistort_point_0(_S1020, dist_coeffs_9, int(12), &_S1021);
    _s_diff_ctx_0->_S1006 = _S1021;
    _s_diff_ctx_0->_S1007 = _S1022;
    float2  uv_12 = _S1021;
    bool _S1023 = !_S1022;
    float3  normal_4;
    if(_S1023)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1024 = !_S1023;
    int _S1025;
    FixedArray<float3 , 4>  points_3;
    if(_S1024)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1026 = length_1(uv_12);
            float3  raydir_19 = make_float3 ((uv_12 / make_float2 ((F32_max((_S1026), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1026))).x, (uv_12 / make_float2 ((F32_max((_S1026), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1026))).y, s_primal_ctx_cos_0(_S1026));
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
            float3  raydir_20 = make_float3 (uv_12.x, uv_12.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_18 = normalize_0(raydir_20);
            }
            else
            {
                raydir_18 = raydir_20;
            }
        }
        float3  _S1027 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1028 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1018) / _S1019;
        float2  _S1029 = _S1028;
        bool _S1030 = undistort_point_0(_S1028, dist_coeffs_9, int(12), &_S1029);
        _s_diff_ctx_0->_S1008 = _S1029;
        _s_diff_ctx_0->_S1009 = _S1030;
        float2  uv_13 = _S1029;
        bool _S1031 = !_S1030;
        if(_S1031)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1032 = !_S1031;
        if(_S1032)
        {
            if(is_fisheye_6)
            {
                float _S1033 = length_1(uv_13);
                float3  raydir_21 = make_float3 ((uv_13 / make_float2 ((F32_max((_S1033), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1033))).x, (uv_13 / make_float2 ((F32_max((_S1033), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1033))).y, s_primal_ctx_cos_0(_S1033));
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
                float3  raydir_22 = make_float3 (uv_13.x, uv_13.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_18 = normalize_0(raydir_22);
                }
                else
                {
                    raydir_18 = raydir_22;
                }
            }
            float3  _S1034 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1025 = int(2);
            points_3[int(0)] = _S1027;
            points_3[int(1)] = _S1034;
            points_3[int(2)] = _S1017;
            points_3[int(3)] = _S1017;
        }
        else
        {
            _S1025 = int(0);
            points_3[int(0)] = _S1027;
            points_3[int(1)] = _S1017;
            points_3[int(2)] = _S1017;
            points_3[int(3)] = _S1017;
        }
        bool _runFlag_0;
        if(_S1025 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1024;
            _S1025 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1035 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1018) / _S1019;
            float2  _S1036 = _S1035;
            bool _S1037 = undistort_point_0(_S1035, dist_coeffs_9, int(12), &_S1036);
            _s_diff_ctx_0->_S1010 = _S1036;
            _s_diff_ctx_0->_S1011 = _S1037;
            float2  uv_14 = _S1036;
            if(!_S1037)
            {
                float3  _S1038 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1025 = int(0);
                normal_4 = _S1038;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1039 = length_1(uv_14);
                    float3  raydir_23 = make_float3 ((uv_14 / make_float2 ((F32_max((_S1039), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1039))).x, (uv_14 / make_float2 ((F32_max((_S1039), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1039))).y, s_primal_ctx_cos_0(_S1039));
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
                    float3  raydir_24 = make_float3 (uv_14.x, uv_14.y, 1.0f);
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
                float2  _S1040 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1018) / _S1019;
                float2  _S1041 = _S1040;
                bool _S1042 = undistort_point_0(_S1040, dist_coeffs_9, int(12), &_S1041);
                _s_diff_ctx_0->_S1012 = _S1041;
                _s_diff_ctx_0->_S1013 = _S1042;
                float2  uv_15 = _S1041;
                bool _S1043 = !_S1042;
                if(_S1043)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1044 = !_S1043;
                int _S1045;
                if(_S1044)
                {
                    if(is_fisheye_6)
                    {
                        float _S1046 = length_1(uv_15);
                        float3  raydir_25 = make_float3 ((uv_15 / make_float2 ((F32_max((_S1046), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1046))).x, (uv_15 / make_float2 ((F32_max((_S1046), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1046))).y, s_primal_ctx_cos_0(_S1046));
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
                        float3  raydir_26 = make_float3 (uv_15.x, uv_15.y, 1.0f);
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
                    _S1045 = int(2);
                }
                else
                {
                    _S1045 = int(0);
                }
                if(_S1045 != int(2))
                {
                    _runFlag_0 = false;
                    _S1025 = _S1045;
                }
                if(_runFlag_0)
                {
                    _S1025 = int(1);
                }
            }
        }
    }
    else
    {
        _S1025 = int(0);
        points_3[int(0)] = _S1017;
        points_3[int(1)] = _S1017;
        points_3[int(2)] = _S1017;
        points_3[int(3)] = _S1017;
    }
    if(!(_S1025 != int(1)))
    {
        float3  _S1047 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1047, _S1047)) != 0.0f)
        {
            normal_4 = _S1047 / make_float3 (length_2(_S1047));
        }
        else
        {
            normal_4 = _S1047;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1048 = *dpdepths_1;
    float3  _S1049 = make_float3 (0.0f);
    float2  _S1050 = _s_diff_ctx_1->_S1006;
    bool _S1051 = !!_s_diff_ctx_1->_S1007;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1052;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1053;
    if(_S1051)
    {
        if(is_fisheye_7)
        {
            float _S1054 = length_1(_S1050);
            float3  raydir_31 = make_float3 ((_S1050 / make_float2 ((F32_max((_S1054), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1054))).x, (_S1050 / make_float2 ((F32_max((_S1054), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1054))).y, s_primal_ctx_cos_0(_S1054));
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
            float3  raydir_32 = make_float3 (_S1050.x, _S1050.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1055 = make_float3 (_S1048.primal_0.x) * raydir_27;
        float2  _S1056 = _s_diff_ctx_1->_S1008;
        bool _S1057 = !!_s_diff_ctx_1->_S1009;
        if(_S1057)
        {
            if(is_fisheye_7)
            {
                float _S1058 = length_1(_S1056);
                float3  raydir_33 = make_float3 ((_S1056 / make_float2 ((F32_max((_S1058), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1058))).x, (_S1056 / make_float2 ((F32_max((_S1058), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1058))).y, s_primal_ctx_cos_0(_S1058));
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
                float3  raydir_34 = make_float3 (_S1056.x, _S1056.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1059 = make_float3 (_S1048.primal_0.y) * raydir_28;
            _S1052 = int(2);
            points_4[int(0)] = _S1055;
            points_4[int(1)] = _S1059;
            points_4[int(2)] = _S1049;
            points_4[int(3)] = _S1049;
        }
        else
        {
            _S1052 = int(0);
            points_4[int(0)] = _S1055;
            points_4[int(1)] = _S1049;
            points_4[int(2)] = _S1049;
            points_4[int(3)] = _S1049;
            raydir_28 = _S1049;
        }
        if(_S1052 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1051;
            _S1052 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S1060 = _s_diff_ctx_1->_S1010;
            if(!_s_diff_ctx_1->_S1011)
            {
                _runFlag_2 = false;
                _S1052 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1061 = length_1(_S1060);
                    float3  raydir_35 = make_float3 ((_S1060 / make_float2 ((F32_max((_S1061), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1061))).x, (_S1060 / make_float2 ((F32_max((_S1061), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1061))).y, s_primal_ctx_cos_0(_S1061));
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
                    float3  raydir_36 = make_float3 (_S1060.x, _S1060.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1048.primal_0.z) * raydir_29;
                float2  _S1062 = _s_diff_ctx_1->_S1012;
                bool _S1063 = !!_s_diff_ctx_1->_S1013;
                int _S1064;
                if(_S1063)
                {
                    if(is_fisheye_7)
                    {
                        float _S1065 = length_1(_S1062);
                        float3  raydir_37 = make_float3 ((_S1062 / make_float2 ((F32_max((_S1065), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1065))).x, (_S1062 / make_float2 ((F32_max((_S1065), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1065))).y, s_primal_ctx_cos_0(_S1065));
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
                        float3  raydir_38 = make_float3 (_S1062.x, _S1062.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1048.primal_0.w) * raydir_30;
                    _S1064 = int(2);
                }
                else
                {
                    _S1064 = int(0);
                    raydir_30 = _S1049;
                }
                if(_S1064 != int(2))
                {
                    _runFlag_3 = false;
                    _S1052 = _S1064;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1052 = int(1);
                }
                float3  _S1066 = raydir_29;
                _runFlag_3 = _S1063;
                raydir_29 = raydir_30;
                raydir_30 = _S1066;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1049;
                raydir_30 = _S1049;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1049;
            raydir_30 = _S1049;
        }
        float3  _S1067 = raydir_27;
        float3  _S1068 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1053 = _S1057;
        raydir_29 = _S1068;
        raydir_30 = _S1067;
    }
    else
    {
        _S1052 = int(0);
        points_4[int(0)] = _S1049;
        points_4[int(1)] = _S1049;
        points_4[int(2)] = _S1049;
        points_4[int(3)] = _S1049;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1049;
        raydir_28 = _S1049;
        _S1053 = false;
        raydir_29 = _S1049;
        raydir_30 = _S1049;
    }
    bool _S1069 = !(_S1052 != int(1));
    float3  _S1070;
    float3  _S1071;
    float3  _S1072;
    float3  _S1073;
    float3  _S1074;
    bool _S1075;
    if(_S1069)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1076 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1077 = s_primal_ctx_cross_0(dx_1, _S1076);
        bool _S1078 = (s_primal_ctx_dot_0(_S1077, _S1077)) != 0.0f;
        if(_S1078)
        {
            float _S1079 = length_2(_S1077);
            float3  _S1080 = make_float3 (_S1079);
            _S1070 = make_float3 (_S1079 * _S1079);
            _S1071 = _S1080;
        }
        else
        {
            _S1070 = _S1049;
            _S1071 = _S1049;
        }
        float3  _S1081 = _S1071;
        _S1075 = _S1078;
        _S1071 = _S1077;
        _S1072 = _S1081;
        _S1073 = dx_1;
        _S1074 = _S1076;
    }
    else
    {
        _S1075 = false;
        _S1070 = _S1049;
        _S1071 = _S1049;
        _S1072 = _S1049;
        _S1073 = _S1049;
        _S1074 = _S1049;
    }
    float4  _S1082 = make_float4 (0.0f);
    if(_S1069)
    {
        if(_S1075)
        {
            float3  _S1083 = _s_dOut_10 / _S1070;
            float3  _S1084 = _S1071 * - _S1083;
            float3  _S1085 = _S1072 * _S1083;
            float _S1086 = _S1084.x + _S1084.y + _S1084.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1087;
            (&_S1087)->primal_0 = _S1071;
            (&_S1087)->differential_0 = _S1049;
            s_bwd_length_impl_1(&_S1087, _S1086);
            _S1070 = _S1085 + _S1087.differential_0;
        }
        else
        {
            _S1070 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1088;
        (&_S1088)->primal_0 = _S1071;
        (&_S1088)->differential_0 = _S1049;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1089;
        (&_S1089)->primal_0 = _S1071;
        (&_S1089)->differential_0 = _S1049;
        s_bwd_prop_dot_0(&_S1088, &_S1089, 0.0f);
        float3  _S1090 = _S1089.differential_0 + _S1088.differential_0 + _S1070;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1091;
        (&_S1091)->primal_0 = _S1073;
        (&_S1091)->differential_0 = _S1049;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1092;
        (&_S1092)->primal_0 = _S1074;
        (&_S1092)->differential_0 = _S1049;
        s_bwd_prop_cross_0(&_S1091, &_S1092, _S1090);
        float3  s_diff_dy_T_1 = - _S1092.differential_0;
        float3  _S1093 = - s_diff_dy_T_1;
        float3  _S1094 = - _S1091.differential_0;
        FixedArray<float3 , 4>  _S1095;
        _S1095[int(0)] = _S1049;
        _S1095[int(1)] = _S1049;
        _S1095[int(2)] = _S1049;
        _S1095[int(3)] = _S1049;
        _S1095[int(2)] = _S1093;
        _S1095[int(3)] = s_diff_dy_T_1;
        _S1095[int(0)] = _S1094;
        _S1095[int(1)] = _S1091.differential_0;
        points_4[int(0)] = _S1095[int(0)];
        points_4[int(1)] = _S1095[int(1)];
        points_4[int(2)] = _S1095[int(2)];
        points_4[int(3)] = _S1095[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1049;
        points_4[int(1)] = _S1049;
        points_4[int(2)] = _S1049;
        points_4[int(3)] = _S1049;
    }
    float4  _S1096;
    if(_S1051)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1097 = points_4;
                FixedArray<float3 , 4>  _S1098 = points_4;
                FixedArray<float3 , 4>  _S1099 = points_4;
                FixedArray<float3 , 4>  _S1100 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1101 = raydir_27 * _S1100[int(3)];
                    float _S1102 = _S1101.x + _S1101.y + _S1101.z;
                    float4  _S1103 = _S1082;
                    *&((&_S1103)->w) = _S1102;
                    points_4[int(0)] = _S1097[int(0)];
                    points_4[int(1)] = _S1098[int(1)];
                    points_4[int(2)] = _S1099[int(2)];
                    points_4[int(3)] = _S1049;
                    _S1096 = _S1103;
                }
                else
                {
                    points_4[int(0)] = _S1097[int(0)];
                    points_4[int(1)] = _S1098[int(1)];
                    points_4[int(2)] = _S1099[int(2)];
                    points_4[int(3)] = _S1100[int(3)];
                    _S1096 = _S1082;
                }
                float3  _S1104 = raydir_28 * points_4[int(2)];
                float _S1105 = _S1104.x + _S1104.y + _S1104.z;
                FixedArray<float3 , 4>  _S1106 = points_4;
                FixedArray<float3 , 4>  _S1107 = points_4;
                float4  _S1108 = _S1082;
                *&((&_S1108)->z) = _S1105;
                float4  _S1109 = _S1096 + _S1108;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1106[int(1)];
                points_4[int(2)] = _S1049;
                points_4[int(3)] = _S1107[int(3)];
                _S1096 = _S1109;
            }
            else
            {
                FixedArray<float3 , 4>  _S1110 = points_4;
                FixedArray<float3 , 4>  _S1111 = points_4;
                FixedArray<float3 , 4>  _S1112 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1110[int(1)];
                points_4[int(2)] = _S1111[int(2)];
                points_4[int(3)] = _S1112[int(3)];
                _S1096 = _S1082;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1113 = points_4;
            FixedArray<float3 , 4>  _S1114 = points_4;
            FixedArray<float3 , 4>  _S1115 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1113[int(1)];
            points_4[int(2)] = _S1114[int(2)];
            points_4[int(3)] = _S1115[int(3)];
            _S1096 = _S1082;
        }
        if(_S1053)
        {
            FixedArray<float3 , 4>  _S1116 = points_4;
            float3  _S1117 = raydir_29 * points_4[int(1)];
            float _S1118 = _S1117.x + _S1117.y + _S1117.z;
            float4  _S1119 = _S1082;
            *&((&_S1119)->y) = _S1118;
            float4  _S1120 = _S1096 + _S1119;
            points_4[int(0)] = _S1049;
            points_4[int(1)] = _S1049;
            points_4[int(2)] = _S1049;
            points_4[int(3)] = _S1049;
            raydir_27 = _S1116[int(0)];
            _S1096 = _S1120;
        }
        else
        {
            FixedArray<float3 , 4>  _S1121 = points_4;
            FixedArray<float3 , 4>  _S1122 = points_4;
            FixedArray<float3 , 4>  _S1123 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1121[int(1)];
            points_4[int(2)] = _S1122[int(2)];
            points_4[int(3)] = _S1123[int(3)];
            raydir_27 = _S1049;
        }
        float3  _S1124 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1125 = _S1124.x + _S1124.y + _S1124.z;
        float4  _S1126 = _S1082;
        *&((&_S1126)->x) = _S1125;
        _S1096 = _S1096 + _S1126;
    }
    else
    {
        _S1096 = _S1082;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1096;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1127, float4  _S1128, FixedArray<float, 10>  * _S1129, bool _S1130, bool _S1131, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1132, float3  _S1133)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1134;
    float3  _S1135 = s_primal_ctx_depth_to_normal_0(_S1127, _S1128, _S1129, _S1130, _S1131, (*_S1132).primal_0, &_S1134);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1136 = _S1134;
    s_bwd_prop_depth_to_normal_0(_S1127, _S1128, _S1129, _S1130, _S1131, _S1132, _S1133, &_S1136);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1137 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1137;
    FixedArray<float, 10>  _S1138 = dist_coeffs_11;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S1138, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_12, bool is_fisheye_9)
{
    float2  _S1139 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_16 = _S1139;
    FixedArray<float, 10>  _S1140 = dist_coeffs_12;
    bool _S1141 = undistort_point_0(_S1139, &_S1140, int(12), &uv_16);
    if(!_S1141)
    {
        return 0.0f;
    }
    float3  raydir_39;
    if(is_fisheye_9)
    {
        float theta_9 = length_1(uv_16);
        float3  raydir_40 = make_float3 ((uv_16 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).x, (uv_16 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).y, (F32_cos((theta_9))));
        raydir_39 = raydir_40 / make_float3 (raydir_40.z);
    }
    else
    {
        raydir_39 = make_float3 (uv_16.x, uv_16.y, 1.0f);
    }
    return float((F32_sign((raydir_39.z)))) / length_2(raydir_39);
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_18, float dOut_21)
{
    float _S1142 = (F32_exp2(((*dpx_18).primal_0))) * 0.69314718246459961f * dOut_21;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1142;
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
    PPISPParams_0 _S1143 = p_0;
    float _S1144 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1145 = (pix_coord_0.x - image_center_0.x) / _S1144;
    float _S1146 = (pix_coord_0.y - image_center_0.y) / _S1144;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_2 = _S1145 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S1146 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_4 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_4 * r2_4;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_4) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1145 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S1146 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_5 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_5 * r2_5;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_5) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1145 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S1146 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_6 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_6 * r2_6;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_6) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S1147 = rgb_out_0;
    float2  bd_0 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S1148 = 0.3333333432674408f + nd_0.x;
    float _S1149 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1149, 1.0f, 0.0f, - _S1148, - _S1149, _S1148, 0.0f), T_0);
    float3  r0_0 = make_float3 (M_3.rows[int(0)].x, M_3.rows[int(0)].y, M_3.rows[int(0)].z);
    float3  r1_0 = make_float3 (M_3.rows[int(1)].x, M_3.rows[int(1)].y, M_3.rows[int(1)].z);
    float3  r2_7 = make_float3 (M_3.rows[int(2)].x, M_3.rows[int(2)].y, M_3.rows[int(2)].z);
    float3  lambda_v_0 = cross_0(r0_0, r1_0);
    float3  lambda_v_1;
    if((dot_0(lambda_v_0, lambda_v_0)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_2 = cross_0(r0_0, r2_7);
        if((dot_0(lambda_v_2, lambda_v_2)) < 9.99999968265522539e-21f)
        {
            lambda_v_1 = cross_0(r1_0, r2_7);
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
    Matrix<float, 3, 3>  H_0 = mul_3(mul_3(T_0, makeMatrix<float, 3, 3> (lambda_v_1.x, 0.0f, 0.0f, 0.0f, lambda_v_1.y, 0.0f, 0.0f, 0.0f, lambda_v_1.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_1;
    if((F32_abs((H_0.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_1 = H_0 * makeMatrix<float, 3, 3> (1.0f / H_0.rows[int(2)].z);
    }
    else
    {
        H_1 = H_0;
    }
    float _S1150 = _S1147.x;
    float _S1151 = _S1147.y;
    float intensity_0 = _S1150 + _S1151 + _S1147.z;
    float3  rgi_out_0 = mul_0(H_1, make_float3 (_S1150, _S1151, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1152 = rgi_out_1.x;
    float _S1153 = rgi_out_1.y;
    float3  _S1154 = clamp_1(make_float3 (_S1152, _S1153, rgi_out_1.z - _S1152 - _S1153), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1155 = _S1154.x;
    float _S1156 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(0)].toe_0))))));
    float _S1157 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(0)].shoulder_0))))));
    float _S1158 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(0)].gamma_0))))));
    float _S1159 = 1.0f / (1.0f + (F32_exp((- _S1143.crf_params_1[int(0)].center_0))));
    float a_1 = _S1157 * _S1159 / lerp_0(_S1156, _S1157, _S1159);
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1155 <= _S1159)
    {
        y_10 = a_1 * (F32_pow((_S1155 / _S1159), (_S1156)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1155) / (1.0f - _S1159)), (_S1157)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1158)));
    float _S1160 = _S1154.y;
    float _S1161 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(1)].toe_0))))));
    float _S1162 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(1)].shoulder_0))))));
    float _S1163 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(1)].gamma_0))))));
    float _S1164 = 1.0f / (1.0f + (F32_exp((- _S1143.crf_params_1[int(1)].center_0))));
    float a_2 = _S1162 * _S1164 / lerp_0(_S1161, _S1162, _S1164);
    float b_3 = 1.0f - a_2;
    if(_S1160 <= _S1164)
    {
        y_10 = a_2 * (F32_pow((_S1160 / _S1164), (_S1161)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1160) / (1.0f - _S1164)), (_S1162)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1163)));
    float _S1165 = _S1154.z;
    float _S1166 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(2)].toe_0))))));
    float _S1167 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(2)].shoulder_0))))));
    float _S1168 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1143.crf_params_1[int(2)].gamma_0))))));
    float _S1169 = 1.0f / (1.0f + (F32_exp((- _S1143.crf_params_1[int(2)].center_0))));
    float a_3 = _S1167 * _S1169 / lerp_0(_S1166, _S1167, _S1169);
    float b_4 = 1.0f - a_3;
    if(_S1165 <= _S1169)
    {
        y_10 = a_3 * (F32_pow((_S1165 / _S1169), (_S1166)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1165) / (1.0f - _S1169)), (_S1167)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1168)));
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
    PPISPParamsRQS_0 _S1170 = p_1;
    float _S1171 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S1172 = (pix_coord_1.x - image_center_1.x) / _S1171;
    float _S1173 = (pix_coord_1.y - image_center_1.y) / _S1171;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_5 = _S1172 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S1173 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_8 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_8 * r2_8;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_8) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_8 + 1.0f, 0.0f, 1.0f);
    float dx_6 = _S1172 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S1173 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_9 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_9 * r2_9;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_9) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S1172 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S1173 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_10 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_10 * r2_10;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_10) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_10 + 1.0f, 0.0f, 1.0f);
    float3  _S1174 = rgb_out_2;
    float2  bd_1 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S1175 = 0.3333333432674408f + nd_1.x;
    float _S1176 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_4 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1176, 1.0f, 0.0f, - _S1175, - _S1176, _S1175, 0.0f), T_1);
    float3  r0_1 = make_float3 (M_4.rows[int(0)].x, M_4.rows[int(0)].y, M_4.rows[int(0)].z);
    float3  r1_1 = make_float3 (M_4.rows[int(1)].x, M_4.rows[int(1)].y, M_4.rows[int(1)].z);
    float3  r2_11 = make_float3 (M_4.rows[int(2)].x, M_4.rows[int(2)].y, M_4.rows[int(2)].z);
    float3  lambda_v_3 = cross_0(r0_1, r1_1);
    float3  lambda_v_4;
    if((dot_0(lambda_v_3, lambda_v_3)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_5 = cross_0(r0_1, r2_11);
        if((dot_0(lambda_v_5, lambda_v_5)) < 9.99999968265522539e-21f)
        {
            lambda_v_4 = cross_0(r1_1, r2_11);
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
    Matrix<float, 3, 3>  H_2 = mul_3(mul_3(T_1, makeMatrix<float, 3, 3> (lambda_v_4.x, 0.0f, 0.0f, 0.0f, lambda_v_4.y, 0.0f, 0.0f, 0.0f, lambda_v_4.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_3;
    if((F32_abs((H_2.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_3 = H_2 * makeMatrix<float, 3, 3> (1.0f / H_2.rows[int(2)].z);
    }
    else
    {
        H_3 = H_2;
    }
    float _S1177 = _S1174.x;
    float _S1178 = _S1174.y;
    float intensity_1 = _S1177 + _S1178 + _S1174.z;
    float3  rgi_out_2 = mul_0(H_3, make_float3 (_S1177, _S1178, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / (rgi_out_2.z + 0.00000999999974738f));
    float _S1179 = rgi_out_3.x;
    float _S1180 = rgi_out_3.y;
    float3  _S1181 = clamp_1(make_float3 (_S1179, _S1180, rgi_out_3.z - _S1179 - _S1180), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S1182 = _S1181.x;
    float g0_1 = (F32_exp((_S1170.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S1170.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S1170.crf_params_0[int(0)].gc_0)));
    float y_11;
    if(_S1182 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S1182 / x0_1;
        float _S1183 = 1.0f - t0_0;
        y_11 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S1183) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S1183);
    }
    else
    {
        float _S1184 = 1.0f - y0_1;
        float _S1185 = 1.0f - x0_1;
        float s1_4 = _S1184 / _S1185;
        float t1_0 = (_S1182 - x0_1) / _S1185;
        float _S1186 = 1.0f - t1_0;
        y_11 = y0_1 + _S1184 * (s1_4 * t1_0 * t1_0 + gc_1 * t1_0 * _S1186) / (s1_4 + (gc_1 + g1_1 - 2.0f * s1_4) * t1_0 * _S1186);
    }
    *&((&rgb_out_3)->x) = y_11;
    float _S1187 = _S1181.y;
    float g0_2 = (F32_exp((_S1170.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S1170.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S1170.crf_params_0[int(1)].gc_0)));
    if(_S1187 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S1187 / x0_2;
        float _S1188 = 1.0f - t0_1;
        y_11 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S1188) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S1188);
    }
    else
    {
        float _S1189 = 1.0f - y0_2;
        float _S1190 = 1.0f - x0_2;
        float s1_5 = _S1189 / _S1190;
        float t1_1 = (_S1187 - x0_2) / _S1190;
        float _S1191 = 1.0f - t1_1;
        y_11 = y0_2 + _S1189 * (s1_5 * t1_1 * t1_1 + gc_2 * t1_1 * _S1191) / (s1_5 + (gc_2 + g1_2 - 2.0f * s1_5) * t1_1 * _S1191);
    }
    *&((&rgb_out_3)->y) = y_11;
    float _S1192 = _S1181.z;
    float g0_3 = (F32_exp((_S1170.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S1170.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S1170.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S1170.crf_params_0[int(2)].gc_0)));
    if(_S1192 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S1192 / x0_3;
        float _S1193 = 1.0f - t0_2;
        y_11 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S1193) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S1193);
    }
    else
    {
        float _S1194 = 1.0f - y0_3;
        float _S1195 = 1.0f - x0_3;
        float s1_6 = _S1194 / _S1195;
        float t1_2 = (_S1192 - x0_3) / _S1195;
        float _S1196 = 1.0f - t1_2;
        y_11 = y0_3 + _S1194 * (s1_6 * t1_2 * t1_2 + gc_3 * t1_2 * _S1196) / (s1_6 + (gc_3 + g1_3 - 2.0f * s1_6) * t1_2 * _S1196);
    }
    *&((&rgb_out_3)->z) = y_11;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1197)
{
    return (F32_exp2((_S1197)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1198, float2  _S1199)
{
    return mul_1(_S1198, _S1199);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1200, Matrix<float, 3, 3>  _S1201)
{
    return mul_3(_S1200, _S1201);
}

inline __device__ float s_primal_ctx_abs_0(float _S1202)
{
    return (F32_abs((_S1202)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1203, float3  _S1204)
{
    return mul_0(_S1203, _S1204);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1205, float3  _S1206, float3  _S1207)
{
    return clamp_1(_S1205, _S1206, _S1207);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1208, float _S1209, float _S1210)
{
    return lerp_0(_S1208, _S1209, _S1210);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1211, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1212, float3  _S1213)
{
    _d_mul_0(_S1211, _S1212, _S1213);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1214, float _S1215)
{
    _d_abs_0(_S1214, _S1215);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1216, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1217, Matrix<float, 3, 3>  _S1218)
{
    mul_2(_S1216, _S1217, _S1218);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1219, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1220, float2  _S1221)
{
    _d_mul_1(_S1219, _S1220, _S1221);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1222, float _S1223)
{
    _d_exp2_0(_S1222, _S1223);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1224 = *dprgb_in_0;
    float3  _S1225 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1226 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1227 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1228 = {
        _S1227, _S1227, _S1227
    };
    float2  _S1229 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1230 = { _S1229, _S1229, _S1229, _S1229 };
    CRFPPISPChannelParams_0 _S1231 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1232 = {
        _S1231, _S1231, _S1231
    };
    PPISPParams_0 _S1233;
    (&_S1233)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S1233)->vignette_params_1 = _S1228;
    (&_S1233)->color_params_1 = _S1230;
    (&_S1233)->crf_params_1 = _S1232;
    (&(&_S1233)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1233)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1234 = dpparams_0->primal_0[int(3)];
    (&(&_S1233)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1235 = dpparams_0->primal_0[int(4)];
    (&(&_S1233)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1236 = dpparams_0->primal_0[int(5)];
    (&(&_S1233)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1233)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1233)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1237 = dpparams_0->primal_0[int(8)];
    (&(&_S1233)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1238 = dpparams_0->primal_0[int(9)];
    (&(&_S1233)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1239 = dpparams_0->primal_0[int(10)];
    (&(&_S1233)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1233)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1233)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1240 = dpparams_0->primal_0[int(13)];
    (&(&_S1233)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1241 = dpparams_0->primal_0[int(14)];
    (&(&_S1233)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1242 = dpparams_0->primal_0[int(15)];
    (&(&_S1233)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1233)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1233)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1233)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1233)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1233)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1233)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1233)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1233)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1243 = dpparams_0->primal_0[int(24)];
    (&(&_S1233)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1244 = dpparams_0->primal_0[int(25)];
    (&(&_S1233)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1245 = dpparams_0->primal_0[int(26)];
    (&(&_S1233)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1246 = dpparams_0->primal_0[int(27)];
    (&(&_S1233)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1247 = dpparams_0->primal_0[int(28)];
    (&(&_S1233)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1248 = dpparams_0->primal_0[int(29)];
    (&(&_S1233)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1249 = dpparams_0->primal_0[int(30)];
    (&(&_S1233)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1250 = dpparams_0->primal_0[int(31)];
    (&(&_S1233)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1251 = dpparams_0->primal_0[int(32)];
    (&(&_S1233)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1252 = dpparams_0->primal_0[int(33)];
    (&(&_S1233)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1253 = dpparams_0->primal_0[int(34)];
    (&(&_S1233)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1254 = dpparams_0->primal_0[int(35)];
    (&(&_S1233)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1255 = _S1233;
    float _S1256 = s_primal_ctx_exp2_0(_S1233.exposure_1);
    float3  _S1257 = make_float3 (_S1256);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S1256);
    float _S1258 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S1259 = (pix_coord_2.x - image_center_2.x) / _S1258;
    float _S1260 = (pix_coord_2.y - image_center_2.y) / _S1258;
    float dx_8 = _S1259 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S1260 - dpparams_0->primal_0[int(2)];
    float r2_12 = dx_8 * dx_8 + dy_6 * dy_6;
    float r4_6 = r2_12 * r2_12;
    float r6_0 = r4_6 * r2_12;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_12 + 1.0f;
    float _S1261 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1262 = rgb_out_4.x * _S1261;
    float3  _S1263 = rgb_out_4;
    *&((&_S1263)->x) = _S1262;
    float dx_9 = _S1259 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S1260 - dpparams_0->primal_0[int(7)];
    float r2_13 = dx_9 * dx_9 + dy_7 * dy_7;
    float r4_7 = r2_13 * r2_13;
    float r6_1 = r4_7 * r2_13;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_13 + 1.0f;
    float _S1264 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1263)->y) = rgb_out_4.y * _S1264;
    float dx_10 = _S1259 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S1260 - dpparams_0->primal_0[int(12)];
    float r2_14 = dx_10 * dx_10 + dy_8 * dy_8;
    float r4_8 = r2_14 * r2_14;
    float r6_2 = r4_8 * r2_14;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_14 + 1.0f;
    float _S1265 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1263)->z) = rgb_out_4.z * _S1265;
    PPISPParams_0 _S1266 = _S1233;
    float2  _S1267 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1233.color_params_1.b_0);
    float2  _S1268 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1233.color_params_1.r_0);
    float2  _S1269 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1233.color_params_1.g_0);
    float2  _S1270 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1233.color_params_1.n_0);
    float _S1271 = 0.3333333432674408f + _S1270.x;
    float _S1272 = 0.3333333432674408f + _S1270.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S1267.x, 1.0f + _S1268.x, _S1269.x, _S1267.y, _S1268.y, 1.0f + _S1269.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1272, 1.0f, 0.0f, - _S1271, - _S1272, _S1271, 0.0f);
    Matrix<float, 3, 3>  _S1273 = s_primal_ctx_mul_1(skew_0, T_2);
    float3  r0_2 = make_float3 (_S1273.rows[int(0)].x, _S1273.rows[int(0)].y, _S1273.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S1273.rows[int(1)].x, _S1273.rows[int(1)].y, _S1273.rows[int(1)].z);
    float3  r2_15 = make_float3 (_S1273.rows[int(2)].x, _S1273.rows[int(2)].y, _S1273.rows[int(2)].z);
    float3  _S1274 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S1275 = (s_primal_ctx_dot_0(_S1274, _S1274)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S1276;
    bool _S1277;
    if(_S1275)
    {
        float3  _S1278 = s_primal_ctx_cross_0(r0_2, r2_15);
        bool _S1279 = (s_primal_ctx_dot_0(_S1278, _S1278)) < 9.99999968265522539e-21f;
        if(_S1279)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_15);
        }
        else
        {
            lambda_v_6 = _S1278;
        }
        _S1277 = _S1279;
        _S1276 = _S1278;
    }
    else
    {
        lambda_v_6 = _S1274;
        _S1277 = false;
        _S1276 = _S1225;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S1280 = s_primal_ctx_mul_1(T_2, D_0);
    Matrix<float, 3, 3>  _S1281 = s_primal_ctx_mul_1(_S1280, S_inv_0);
    bool _S1282 = (s_primal_ctx_abs_0(_S1281.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S1283;
    float _S1284;
    if(_S1282)
    {
        float inv_s_0 = 1.0f / _S1281.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1285 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1286 = _S1281.rows[int(2)].z * _S1281.rows[int(2)].z;
        H_4 = _S1281 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1283 = _S1285;
        _S1284 = _S1286;
    }
    else
    {
        H_4 = _S1281;
        _S1283 = _S1226;
        _S1284 = 0.0f;
    }
    float _S1287 = _S1263.x;
    float _S1288 = _S1263.y;
    float intensity_2 = _S1287 + _S1288 + _S1263.z;
    float3  rgi_in_0 = make_float3 (_S1287, _S1288, intensity_2);
    float3  _S1289 = s_primal_ctx_mul_2(H_4, rgi_in_0);
    float _S1290 = _S1289.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S1290;
    float3  _S1291 = make_float3 (norm_factor_0);
    float _S1292 = _S1290 * _S1290;
    float3  rgi_out_4 = _S1289 * make_float3 (norm_factor_0);
    float _S1293 = rgi_out_4.x;
    float _S1294 = rgi_out_4.y;
    float3  _S1295 = make_float3 (_S1293, _S1294, rgi_out_4.z - _S1293 - _S1294);
    float3  _S1296 = make_float3 (0.0f);
    float3  _S1297 = make_float3 (1.0f);
    float3  _S1298 = s_primal_ctx_clamp_1(_S1295, _S1296, _S1297);
    float _S1299 = _S1298.x;
    float _S1300 = 1.0f + s_primal_ctx_exp_0(_S1243);
    float _S1301 = 0.30000001192092896f + s_primal_ctx_log_0(_S1300);
    float _S1302 = 1.0f + s_primal_ctx_exp_0(_S1244);
    float _S1303 = 0.30000001192092896f + s_primal_ctx_log_0(_S1302);
    float _S1304 = 1.0f + s_primal_ctx_exp_0(_S1245);
    float _S1305 = 0.10000000149011612f + s_primal_ctx_log_0(_S1304);
    float _S1306 = - _S1246;
    float _S1307 = 1.0f + s_primal_ctx_exp_0(_S1306);
    float _S1308 = 1.0f / _S1307;
    float _S1309 = _S1307 * _S1307;
    float _S1310 = s_primal_ctx_lerp_0(_S1301, _S1303, _S1308);
    float _S1311 = _S1303 * _S1308;
    float a_4 = _S1311 / _S1310;
    float _S1312 = _S1310 * _S1310;
    float b_5 = 1.0f - a_4;
    bool _S1313 = _S1299 <= _S1308;
    float y_12;
    float _S1314;
    float _S1315;
    float _S1316;
    float _S1317;
    float _S1318;
    float _S1319;
    float _S1320;
    float _S1321;
    if(_S1313)
    {
        float _S1322 = _S1299 / _S1308;
        float _S1323 = _S1308 * _S1308;
        float _S1324 = s_primal_ctx_pow_0(_S1322, _S1301);
        y_12 = a_4 * _S1324;
        _S1314 = _S1324;
        _S1315 = _S1322;
        _S1316 = _S1323;
        _S1317 = 0.0f;
        _S1318 = 0.0f;
        _S1319 = 0.0f;
        _S1320 = 0.0f;
        _S1321 = 0.0f;
    }
    else
    {
        float _S1325 = 1.0f - _S1299;
        float _S1326 = 1.0f - _S1308;
        float _S1327 = _S1325 / _S1326;
        float _S1328 = _S1326 * _S1326;
        float _S1329 = s_primal_ctx_pow_0(_S1327, _S1303);
        y_12 = 1.0f - b_5 * _S1329;
        _S1314 = 0.0f;
        _S1315 = 0.0f;
        _S1316 = 0.0f;
        _S1317 = _S1329;
        _S1318 = _S1327;
        _S1319 = _S1328;
        _S1320 = _S1325;
        _S1321 = _S1326;
    }
    float _S1330 = (F32_max((0.0f), (y_12)));
    float _S1331 = _S1298.y;
    float _S1332 = 1.0f + s_primal_ctx_exp_0(_S1247);
    float _S1333 = 0.30000001192092896f + s_primal_ctx_log_0(_S1332);
    float _S1334 = 1.0f + s_primal_ctx_exp_0(_S1248);
    float _S1335 = 0.30000001192092896f + s_primal_ctx_log_0(_S1334);
    float _S1336 = 1.0f + s_primal_ctx_exp_0(_S1249);
    float _S1337 = 0.10000000149011612f + s_primal_ctx_log_0(_S1336);
    float _S1338 = - _S1250;
    float _S1339 = 1.0f + s_primal_ctx_exp_0(_S1338);
    float _S1340 = 1.0f / _S1339;
    float _S1341 = _S1339 * _S1339;
    float _S1342 = s_primal_ctx_lerp_0(_S1333, _S1335, _S1340);
    float _S1343 = _S1335 * _S1340;
    float a_5 = _S1343 / _S1342;
    float _S1344 = _S1342 * _S1342;
    float b_6 = 1.0f - a_5;
    bool _S1345 = _S1331 <= _S1340;
    float y_13;
    float _S1346;
    float _S1347;
    float _S1348;
    float _S1349;
    float _S1350;
    float _S1351;
    float _S1352;
    float _S1353;
    if(_S1345)
    {
        float _S1354 = _S1331 / _S1340;
        float _S1355 = _S1340 * _S1340;
        float _S1356 = s_primal_ctx_pow_0(_S1354, _S1333);
        y_13 = a_5 * _S1356;
        _S1346 = _S1356;
        _S1347 = _S1354;
        _S1348 = _S1355;
        _S1349 = 0.0f;
        _S1350 = 0.0f;
        _S1351 = 0.0f;
        _S1352 = 0.0f;
        _S1353 = 0.0f;
    }
    else
    {
        float _S1357 = 1.0f - _S1331;
        float _S1358 = 1.0f - _S1340;
        float _S1359 = _S1357 / _S1358;
        float _S1360 = _S1358 * _S1358;
        float _S1361 = s_primal_ctx_pow_0(_S1359, _S1335);
        y_13 = 1.0f - b_6 * _S1361;
        _S1346 = 0.0f;
        _S1347 = 0.0f;
        _S1348 = 0.0f;
        _S1349 = _S1361;
        _S1350 = _S1359;
        _S1351 = _S1360;
        _S1352 = _S1357;
        _S1353 = _S1358;
    }
    float _S1362 = (F32_max((0.0f), (y_13)));
    float _S1363 = _S1298.z;
    float _S1364 = 1.0f + s_primal_ctx_exp_0(_S1251);
    float _S1365 = 0.30000001192092896f + s_primal_ctx_log_0(_S1364);
    float _S1366 = 1.0f + s_primal_ctx_exp_0(_S1252);
    float _S1367 = 0.30000001192092896f + s_primal_ctx_log_0(_S1366);
    float _S1368 = 1.0f + s_primal_ctx_exp_0(_S1253);
    float _S1369 = 0.10000000149011612f + s_primal_ctx_log_0(_S1368);
    float _S1370 = - _S1254;
    float _S1371 = 1.0f + s_primal_ctx_exp_0(_S1370);
    float _S1372 = 1.0f / _S1371;
    float _S1373 = _S1371 * _S1371;
    float _S1374 = s_primal_ctx_lerp_0(_S1365, _S1367, _S1372);
    float _S1375 = _S1367 * _S1372;
    float a_6 = _S1375 / _S1374;
    float _S1376 = _S1374 * _S1374;
    float b_7 = 1.0f - a_6;
    bool _S1377 = _S1363 <= _S1372;
    float y_14;
    float _S1378;
    float _S1379;
    float _S1380;
    float _S1381;
    float _S1382;
    float _S1383;
    float _S1384;
    float _S1385;
    if(_S1377)
    {
        float _S1386 = _S1363 / _S1372;
        float _S1387 = _S1372 * _S1372;
        float _S1388 = s_primal_ctx_pow_0(_S1386, _S1365);
        y_14 = a_6 * _S1388;
        _S1378 = _S1388;
        _S1379 = _S1386;
        _S1380 = _S1387;
        _S1381 = 0.0f;
        _S1382 = 0.0f;
        _S1383 = 0.0f;
        _S1384 = 0.0f;
        _S1385 = 0.0f;
    }
    else
    {
        float _S1389 = 1.0f - _S1363;
        float _S1390 = 1.0f - _S1372;
        float _S1391 = _S1389 / _S1390;
        float _S1392 = _S1390 * _S1390;
        float _S1393 = s_primal_ctx_pow_0(_S1391, _S1367);
        y_14 = 1.0f - b_7 * _S1393;
        _S1378 = 0.0f;
        _S1379 = 0.0f;
        _S1380 = 0.0f;
        _S1381 = _S1393;
        _S1382 = _S1391;
        _S1383 = _S1392;
        _S1384 = _S1389;
        _S1385 = _S1390;
    }
    float _S1394 = (F32_max((0.0f), (y_14)));
    DiffPair_float_0 _S1395;
    (&_S1395)->primal_0 = _S1394;
    (&_S1395)->differential_0 = 0.0f;
    DiffPair_float_0 _S1396;
    (&_S1396)->primal_0 = _S1369;
    (&_S1396)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1395, &_S1396, _s_dOut_11.z);
    DiffPair_float_0 _S1397 = _S1396;
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = 0.0f;
    (&_S1398)->differential_0 = 0.0f;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = y_14;
    (&_S1399)->differential_0 = 0.0f;
    _d_max_0(&_S1398, &_S1399, _S1395.differential_0);
    DiffPair_float_0 _S1400 = _S1399;
    if(_S1377)
    {
        float _S1401 = a_6 * _S1400.differential_0;
        float _S1402 = _S1378 * _S1400.differential_0;
        DiffPair_float_0 _S1403;
        (&_S1403)->primal_0 = _S1379;
        (&_S1403)->differential_0 = 0.0f;
        DiffPair_float_0 _S1404;
        (&_S1404)->primal_0 = _S1365;
        (&_S1404)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1403, &_S1404, _S1401);
        float _S1405 = _S1403.differential_0 / _S1380;
        float _S1406 = _S1363 * - _S1405;
        float _S1407 = _S1372 * _S1405;
        y_14 = 0.0f;
        _S1378 = _S1402;
        _S1379 = _S1406;
        _S1380 = 0.0f;
        _S1381 = _S1404.differential_0;
        _S1382 = _S1407;
    }
    else
    {
        float _S1408 = - _S1400.differential_0;
        float _S1409 = b_7 * _S1408;
        float _S1410 = _S1381 * _S1408;
        DiffPair_float_0 _S1411;
        (&_S1411)->primal_0 = _S1382;
        (&_S1411)->differential_0 = 0.0f;
        DiffPair_float_0 _S1412;
        (&_S1412)->primal_0 = _S1367;
        (&_S1412)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1411, &_S1412, _S1409);
        float _S1413 = _S1411.differential_0 / _S1383;
        float _S1414 = - (_S1384 * - _S1413);
        float _S1415 = - (_S1385 * _S1413);
        y_14 = _S1410;
        _S1378 = 0.0f;
        _S1379 = _S1414;
        _S1380 = _S1412.differential_0;
        _S1381 = 0.0f;
        _S1382 = _S1415;
    }
    float _S1416 = (- y_14 + _S1378) / _S1376;
    float _S1417 = _S1375 * - _S1416;
    float _S1418 = _S1374 * _S1416;
    float _S1419 = _S1367 * _S1418;
    float _S1420 = _S1372 * _S1418;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1365;
    (&_S1421)->differential_0 = 0.0f;
    DiffPair_float_0 _S1422;
    (&_S1422)->primal_0 = _S1367;
    (&_S1422)->differential_0 = 0.0f;
    DiffPair_float_0 _S1423;
    (&_S1423)->primal_0 = _S1372;
    (&_S1423)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1421, &_S1422, &_S1423, _S1417);
    float _S1424 = - ((_S1419 + _S1423.differential_0 + _S1379) / _S1373);
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = _S1370;
    (&_S1425)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1425, _S1424);
    float _S1426 = - _S1425.differential_0;
    DiffPair_float_0 _S1427;
    (&_S1427)->primal_0 = _S1368;
    (&_S1427)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1427, _S1397.differential_0);
    DiffPair_float_0 _S1428;
    (&_S1428)->primal_0 = _S1253;
    (&_S1428)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1428, _S1427.differential_0);
    DiffPair_float_0 _S1429 = _S1428;
    float _S1430 = _S1420 + _S1422.differential_0 + _S1380;
    DiffPair_float_0 _S1431;
    (&_S1431)->primal_0 = _S1366;
    (&_S1431)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1431, _S1430);
    DiffPair_float_0 _S1432;
    (&_S1432)->primal_0 = _S1252;
    (&_S1432)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1432, _S1431.differential_0);
    DiffPair_float_0 _S1433 = _S1432;
    float _S1434 = _S1421.differential_0 + _S1381;
    DiffPair_float_0 _S1435;
    (&_S1435)->primal_0 = _S1364;
    (&_S1435)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1435, _S1434);
    DiffPair_float_0 _S1436;
    (&_S1436)->primal_0 = _S1251;
    (&_S1436)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1436, _S1435.differential_0);
    DiffPair_float_0 _S1437 = _S1436;
    float3  _S1438 = make_float3 (0.0f, 0.0f, _S1382);
    DiffPair_float_0 _S1439;
    (&_S1439)->primal_0 = _S1362;
    (&_S1439)->differential_0 = 0.0f;
    DiffPair_float_0 _S1440;
    (&_S1440)->primal_0 = _S1337;
    (&_S1440)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1439, &_S1440, _s_dOut_11.y);
    DiffPair_float_0 _S1441 = _S1440;
    DiffPair_float_0 _S1442;
    (&_S1442)->primal_0 = 0.0f;
    (&_S1442)->differential_0 = 0.0f;
    DiffPair_float_0 _S1443;
    (&_S1443)->primal_0 = y_13;
    (&_S1443)->differential_0 = 0.0f;
    _d_max_0(&_S1442, &_S1443, _S1439.differential_0);
    DiffPair_float_0 _S1444 = _S1443;
    if(_S1345)
    {
        float _S1445 = a_5 * _S1444.differential_0;
        float _S1446 = _S1346 * _S1444.differential_0;
        DiffPair_float_0 _S1447;
        (&_S1447)->primal_0 = _S1347;
        (&_S1447)->differential_0 = 0.0f;
        DiffPair_float_0 _S1448;
        (&_S1448)->primal_0 = _S1333;
        (&_S1448)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1447, &_S1448, _S1445);
        float _S1449 = _S1447.differential_0 / _S1348;
        float _S1450 = _S1331 * - _S1449;
        float _S1451 = _S1340 * _S1449;
        y_13 = 0.0f;
        _S1346 = _S1446;
        _S1347 = _S1450;
        _S1348 = 0.0f;
        _S1349 = _S1448.differential_0;
        _S1350 = _S1451;
    }
    else
    {
        float _S1452 = - _S1444.differential_0;
        float _S1453 = b_6 * _S1452;
        float _S1454 = _S1349 * _S1452;
        DiffPair_float_0 _S1455;
        (&_S1455)->primal_0 = _S1350;
        (&_S1455)->differential_0 = 0.0f;
        DiffPair_float_0 _S1456;
        (&_S1456)->primal_0 = _S1335;
        (&_S1456)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1455, &_S1456, _S1453);
        float _S1457 = _S1455.differential_0 / _S1351;
        float _S1458 = - (_S1352 * - _S1457);
        float _S1459 = - (_S1353 * _S1457);
        y_13 = _S1454;
        _S1346 = 0.0f;
        _S1347 = _S1458;
        _S1348 = _S1456.differential_0;
        _S1349 = 0.0f;
        _S1350 = _S1459;
    }
    float _S1460 = (- y_13 + _S1346) / _S1344;
    float _S1461 = _S1343 * - _S1460;
    float _S1462 = _S1342 * _S1460;
    float _S1463 = _S1335 * _S1462;
    float _S1464 = _S1340 * _S1462;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = _S1333;
    (&_S1465)->differential_0 = 0.0f;
    DiffPair_float_0 _S1466;
    (&_S1466)->primal_0 = _S1335;
    (&_S1466)->differential_0 = 0.0f;
    DiffPair_float_0 _S1467;
    (&_S1467)->primal_0 = _S1340;
    (&_S1467)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1465, &_S1466, &_S1467, _S1461);
    float _S1468 = - ((_S1463 + _S1467.differential_0 + _S1347) / _S1341);
    DiffPair_float_0 _S1469;
    (&_S1469)->primal_0 = _S1338;
    (&_S1469)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1469, _S1468);
    float _S1470 = - _S1469.differential_0;
    DiffPair_float_0 _S1471;
    (&_S1471)->primal_0 = _S1336;
    (&_S1471)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1471, _S1441.differential_0);
    DiffPair_float_0 _S1472;
    (&_S1472)->primal_0 = _S1249;
    (&_S1472)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1472, _S1471.differential_0);
    DiffPair_float_0 _S1473 = _S1472;
    float _S1474 = _S1464 + _S1466.differential_0 + _S1348;
    DiffPair_float_0 _S1475;
    (&_S1475)->primal_0 = _S1334;
    (&_S1475)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1475, _S1474);
    DiffPair_float_0 _S1476;
    (&_S1476)->primal_0 = _S1248;
    (&_S1476)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1476, _S1475.differential_0);
    DiffPair_float_0 _S1477 = _S1476;
    float _S1478 = _S1465.differential_0 + _S1349;
    DiffPair_float_0 _S1479;
    (&_S1479)->primal_0 = _S1332;
    (&_S1479)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1479, _S1478);
    DiffPair_float_0 _S1480;
    (&_S1480)->primal_0 = _S1247;
    (&_S1480)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1480, _S1479.differential_0);
    DiffPair_float_0 _S1481 = _S1480;
    float3  _S1482 = _S1438 + make_float3 (0.0f, _S1350, 0.0f);
    DiffPair_float_0 _S1483;
    (&_S1483)->primal_0 = _S1330;
    (&_S1483)->differential_0 = 0.0f;
    DiffPair_float_0 _S1484;
    (&_S1484)->primal_0 = _S1305;
    (&_S1484)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1483, &_S1484, _s_dOut_11.x);
    DiffPair_float_0 _S1485 = _S1484;
    DiffPair_float_0 _S1486;
    (&_S1486)->primal_0 = 0.0f;
    (&_S1486)->differential_0 = 0.0f;
    DiffPair_float_0 _S1487;
    (&_S1487)->primal_0 = y_12;
    (&_S1487)->differential_0 = 0.0f;
    _d_max_0(&_S1486, &_S1487, _S1483.differential_0);
    DiffPair_float_0 _S1488 = _S1487;
    if(_S1313)
    {
        float _S1489 = a_4 * _S1488.differential_0;
        float _S1490 = _S1314 * _S1488.differential_0;
        DiffPair_float_0 _S1491;
        (&_S1491)->primal_0 = _S1315;
        (&_S1491)->differential_0 = 0.0f;
        DiffPair_float_0 _S1492;
        (&_S1492)->primal_0 = _S1301;
        (&_S1492)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1491, &_S1492, _S1489);
        float _S1493 = _S1491.differential_0 / _S1316;
        float _S1494 = _S1299 * - _S1493;
        float _S1495 = _S1308 * _S1493;
        y_12 = 0.0f;
        _S1314 = _S1490;
        _S1315 = _S1494;
        _S1316 = 0.0f;
        _S1317 = _S1492.differential_0;
        _S1318 = _S1495;
    }
    else
    {
        float _S1496 = - _S1488.differential_0;
        float _S1497 = b_5 * _S1496;
        float _S1498 = _S1317 * _S1496;
        DiffPair_float_0 _S1499;
        (&_S1499)->primal_0 = _S1318;
        (&_S1499)->differential_0 = 0.0f;
        DiffPair_float_0 _S1500;
        (&_S1500)->primal_0 = _S1303;
        (&_S1500)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1499, &_S1500, _S1497);
        float _S1501 = _S1499.differential_0 / _S1319;
        float _S1502 = - (_S1320 * - _S1501);
        float _S1503 = - (_S1321 * _S1501);
        y_12 = _S1498;
        _S1314 = 0.0f;
        _S1315 = _S1502;
        _S1316 = _S1500.differential_0;
        _S1317 = 0.0f;
        _S1318 = _S1503;
    }
    float _S1504 = (- y_12 + _S1314) / _S1312;
    float _S1505 = _S1311 * - _S1504;
    float _S1506 = _S1310 * _S1504;
    float _S1507 = _S1303 * _S1506;
    float _S1508 = _S1308 * _S1506;
    DiffPair_float_0 _S1509;
    (&_S1509)->primal_0 = _S1301;
    (&_S1509)->differential_0 = 0.0f;
    DiffPair_float_0 _S1510;
    (&_S1510)->primal_0 = _S1303;
    (&_S1510)->differential_0 = 0.0f;
    DiffPair_float_0 _S1511;
    (&_S1511)->primal_0 = _S1308;
    (&_S1511)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1509, &_S1510, &_S1511, _S1505);
    float _S1512 = - ((_S1507 + _S1511.differential_0 + _S1315) / _S1309);
    DiffPair_float_0 _S1513;
    (&_S1513)->primal_0 = _S1306;
    (&_S1513)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1513, _S1512);
    float _S1514 = - _S1513.differential_0;
    DiffPair_float_0 _S1515;
    (&_S1515)->primal_0 = _S1304;
    (&_S1515)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1515, _S1485.differential_0);
    DiffPair_float_0 _S1516;
    (&_S1516)->primal_0 = _S1245;
    (&_S1516)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1516, _S1515.differential_0);
    DiffPair_float_0 _S1517 = _S1516;
    float _S1518 = _S1508 + _S1510.differential_0 + _S1316;
    DiffPair_float_0 _S1519;
    (&_S1519)->primal_0 = _S1302;
    (&_S1519)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1519, _S1518);
    DiffPair_float_0 _S1520;
    (&_S1520)->primal_0 = _S1244;
    (&_S1520)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1520, _S1519.differential_0);
    DiffPair_float_0 _S1521 = _S1520;
    float _S1522 = _S1509.differential_0 + _S1317;
    DiffPair_float_0 _S1523;
    (&_S1523)->primal_0 = _S1300;
    (&_S1523)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1523, _S1522);
    DiffPair_float_0 _S1524;
    (&_S1524)->primal_0 = _S1243;
    (&_S1524)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1524, _S1523.differential_0);
    DiffPair_float_0 _S1525 = _S1524;
    float3  _S1526 = _S1482 + make_float3 (_S1318, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1527;
    (&_S1527)->primal_0 = _S1295;
    (&_S1527)->differential_0 = _S1225;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1528;
    (&_S1528)->primal_0 = _S1296;
    (&_S1528)->differential_0 = _S1225;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1529;
    (&_S1529)->primal_0 = _S1297;
    (&_S1529)->differential_0 = _S1225;
    s_bwd_prop_clamp_1(&_S1527, &_S1528, &_S1529, _S1526);
    float _S1530 = - _S1527.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1527.differential_0.x + _S1530, _S1527.differential_0.y + _S1530, _S1527.differential_0.z);
    float3  _S1531 = _S1289 * s_diff_rgi_out_T_0;
    float _S1532 = (_S1531.x + _S1531.y + _S1531.z) / _S1292;
    float _S1533 = _S1290 * _S1532;
    float3  _S1534 = _S1291 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_2 * - _S1532);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1535;
    (&_S1535)->primal_0 = H_4;
    (&_S1535)->differential_0 = _S1226;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1536;
    (&_S1536)->primal_0 = rgi_in_0;
    (&_S1536)->differential_0 = _S1225;
    s_bwd_prop_mul_1(&_S1535, &_S1536, _S1534);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1537 = _S1535;
    float _S1538 = _S1533 + _S1536.differential_0.z;
    float _S1539 = _S1536.differential_0.y + _S1538;
    float _S1540 = _S1536.differential_0.x + _S1538;
    float3  _S1541 = make_float3 (_S1540, _S1539, _S1538);
    if(_S1282)
    {
        Matrix<float, 3, 3>  _S1542 = _S1281 * _S1537.differential_0;
        Matrix<float, 3, 3>  _S1543 = _S1283 * _S1537.differential_0;
        _S1284 = - ((_S1542.rows[int(0)].x + _S1542.rows[int(0)].y + _S1542.rows[int(0)].z + _S1542.rows[int(1)].x + _S1542.rows[int(1)].y + _S1542.rows[int(1)].z + _S1542.rows[int(2)].x + _S1542.rows[int(2)].y + _S1542.rows[int(2)].z) / _S1284);
        H_4 = _S1543;
    }
    else
    {
        _S1284 = 0.0f;
        H_4 = _S1537.differential_0;
    }
    DiffPair_float_0 _S1544;
    (&_S1544)->primal_0 = _S1281.rows[int(2)].z;
    (&_S1544)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1544, 0.0f);
    float _S1545 = _S1544.differential_0 + _S1284;
    float3  _S1546 = _S1225;
    *&((&_S1546)->z) = _S1545;
    Matrix<float, 3, 3>  _S1547 = _S1226;
    _S1547[int(2)] = _S1546;
    Matrix<float, 3, 3>  _S1548 = H_4 + _S1547;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1549;
    (&_S1549)->primal_0 = _S1280;
    (&_S1549)->differential_0 = _S1226;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1550;
    (&_S1550)->primal_0 = S_inv_0;
    (&_S1550)->differential_0 = _S1226;
    s_bwd_prop_mul_2(&_S1549, &_S1550, _S1548);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1551;
    (&_S1551)->primal_0 = T_2;
    (&_S1551)->differential_0 = _S1226;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1552;
    (&_S1552)->primal_0 = D_0;
    (&_S1552)->differential_0 = _S1226;
    s_bwd_prop_mul_2(&_S1551, &_S1552, _S1549.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1553 = _S1551;
    float3  _S1554 = make_float3 (_S1552.differential_0.rows[int(0)].x, _S1552.differential_0.rows[int(1)].y, _S1552.differential_0.rows[int(2)].z);
    float3  _S1555;
    if(_S1275)
    {
        if(_S1277)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1556;
            (&_S1556)->primal_0 = r1_2;
            (&_S1556)->differential_0 = _S1225;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1557;
            (&_S1557)->primal_0 = r2_15;
            (&_S1557)->differential_0 = _S1225;
            s_bwd_prop_cross_0(&_S1556, &_S1557, _S1554);
            _S1263 = _S1225;
            lambda_v_6 = _S1557.differential_0;
            _S1555 = _S1556.differential_0;
        }
        else
        {
            _S1263 = _S1554;
            lambda_v_6 = _S1225;
            _S1555 = _S1225;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1558;
        (&_S1558)->primal_0 = _S1276;
        (&_S1558)->differential_0 = _S1225;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1559;
        (&_S1559)->primal_0 = _S1276;
        (&_S1559)->differential_0 = _S1225;
        s_bwd_prop_dot_0(&_S1558, &_S1559, 0.0f);
        float3  _S1560 = _S1559.differential_0 + _S1558.differential_0 + _S1263;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1561;
        (&_S1561)->primal_0 = r0_2;
        (&_S1561)->differential_0 = _S1225;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1562;
        (&_S1562)->primal_0 = r2_15;
        (&_S1562)->differential_0 = _S1225;
        s_bwd_prop_cross_0(&_S1561, &_S1562, _S1560);
        float3  _S1563 = _S1562.differential_0 + lambda_v_6;
        _S1263 = _S1225;
        lambda_v_6 = _S1563;
        _S1276 = _S1555;
        _S1555 = _S1561.differential_0;
    }
    else
    {
        _S1263 = _S1554;
        lambda_v_6 = _S1225;
        _S1276 = _S1225;
        _S1555 = _S1225;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1564;
    (&_S1564)->primal_0 = _S1274;
    (&_S1564)->differential_0 = _S1225;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1565;
    (&_S1565)->primal_0 = _S1274;
    (&_S1565)->differential_0 = _S1225;
    s_bwd_prop_dot_0(&_S1564, &_S1565, 0.0f);
    float3  _S1566 = _S1565.differential_0 + _S1564.differential_0 + _S1263;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1567;
    (&_S1567)->primal_0 = r0_2;
    (&_S1567)->differential_0 = _S1225;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1568;
    (&_S1568)->primal_0 = r1_2;
    (&_S1568)->differential_0 = _S1225;
    s_bwd_prop_cross_0(&_S1567, &_S1568, _S1566);
    float3  _S1569 = _S1225;
    *&((&_S1569)->z) = lambda_v_6.z;
    *&((&_S1569)->y) = lambda_v_6.y;
    *&((&_S1569)->x) = lambda_v_6.x;
    float3  _S1570 = _S1568.differential_0 + _S1276;
    float3  _S1571 = _S1225;
    *&((&_S1571)->z) = _S1570.z;
    *&((&_S1571)->y) = _S1570.y;
    *&((&_S1571)->x) = _S1570.x;
    float3  _S1572 = _S1567.differential_0 + _S1555;
    float3  _S1573 = _S1225;
    *&((&_S1573)->z) = _S1572.z;
    *&((&_S1573)->y) = _S1572.y;
    *&((&_S1573)->x) = _S1572.x;
    Matrix<float, 3, 3>  _S1574 = _S1226;
    _S1574[int(2)] = _S1569;
    _S1574[int(1)] = _S1571;
    _S1574[int(0)] = _S1573;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1575;
    (&_S1575)->primal_0 = skew_0;
    (&_S1575)->differential_0 = _S1226;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1576;
    (&_S1576)->primal_0 = T_2;
    (&_S1576)->differential_0 = _S1226;
    s_bwd_prop_mul_2(&_S1575, &_S1576, _S1574);
    Matrix<float, 3, 3>  _S1577 = _S1576.differential_0 + _S1553.differential_0;
    float2  _S1578 = make_float2 (_S1575.differential_0.rows[int(2)].y + - _S1575.differential_0.rows[int(1)].z, _S1575.differential_0.rows[int(0)].z + - _S1575.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1579 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1580;
    (&_S1580)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1580)->differential_0 = _S1579;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1581;
    (&_S1581)->primal_0 = _S1266.color_params_1.n_0;
    (&_S1581)->differential_0 = _S1229;
    s_bwd_prop_mul_3(&_S1580, &_S1581, _S1578);
    float2  _S1582 = make_float2 (_S1577.rows[int(0)].z, _S1577.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1583;
    (&_S1583)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1583)->differential_0 = _S1579;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1584;
    (&_S1584)->primal_0 = _S1266.color_params_1.g_0;
    (&_S1584)->differential_0 = _S1229;
    s_bwd_prop_mul_3(&_S1583, &_S1584, _S1582);
    float2  _S1585 = make_float2 (_S1577.rows[int(0)].y, _S1577.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1586;
    (&_S1586)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1586)->differential_0 = _S1579;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1587;
    (&_S1587)->primal_0 = _S1266.color_params_1.r_0;
    (&_S1587)->differential_0 = _S1229;
    s_bwd_prop_mul_3(&_S1586, &_S1587, _S1585);
    float2  _S1588 = make_float2 (_S1577.rows[int(0)].x, _S1577.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1589;
    (&_S1589)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1589)->differential_0 = _S1579;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1590;
    (&_S1590)->primal_0 = _S1266.color_params_1.b_0;
    (&_S1590)->differential_0 = _S1229;
    s_bwd_prop_mul_3(&_S1589, &_S1590, _S1588);
    ColorPPISPParams_0 _S1591 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1591)->n_0 = _S1581.differential_0;
    (&_S1591)->g_0 = _S1584.differential_0;
    (&_S1591)->r_0 = _S1587.differential_0;
    (&_S1591)->b_0 = _S1590.differential_0;
    _S1263 = _S1541;
    *&((&_S1263)->z) = 0.0f;
    float _S1592 = rgb_out_4.z * _S1538;
    float _S1593 = _S1265 * _S1538;
    DiffPair_float_0 _S1594;
    (&_S1594)->primal_0 = falloff_2;
    (&_S1594)->differential_0 = 0.0f;
    DiffPair_float_0 _S1595;
    (&_S1595)->primal_0 = 0.0f;
    (&_S1595)->differential_0 = 0.0f;
    DiffPair_float_0 _S1596;
    (&_S1596)->primal_0 = 1.0f;
    (&_S1596)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1594, &_S1595, &_S1596, _S1592);
    float _S1597 = r2_14 * _S1594.differential_0;
    float _S1598 = r4_8 * _S1594.differential_0;
    float s_diff_r6_T_0 = _S1242 * _S1594.differential_0;
    float _S1599 = r6_2 * _S1594.differential_0;
    float _S1600 = r2_14 * (_S1241 * _S1594.differential_0 + r2_14 * s_diff_r6_T_0);
    float _S1601 = _S1240 * _S1594.differential_0 + r4_8 * s_diff_r6_T_0 + _S1600 + _S1600;
    float _S1602 = dy_8 * _S1601;
    float _S1603 = dx_10 * _S1601;
    float _S1604 = - (_S1602 + _S1602);
    float _S1605 = - (_S1603 + _S1603);
    *&((&_S1263)->y) = 0.0f;
    float _S1606 = rgb_out_4.y * _S1539;
    float _S1607 = _S1264 * _S1539;
    DiffPair_float_0 _S1608;
    (&_S1608)->primal_0 = falloff_1;
    (&_S1608)->differential_0 = 0.0f;
    DiffPair_float_0 _S1609;
    (&_S1609)->primal_0 = 0.0f;
    (&_S1609)->differential_0 = 0.0f;
    DiffPair_float_0 _S1610;
    (&_S1610)->primal_0 = 1.0f;
    (&_S1610)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1608, &_S1609, &_S1610, _S1606);
    float _S1611 = r2_13 * _S1608.differential_0;
    float _S1612 = r4_7 * _S1608.differential_0;
    float s_diff_r6_T_1 = _S1239 * _S1608.differential_0;
    float _S1613 = r6_1 * _S1608.differential_0;
    float _S1614 = r2_13 * (_S1238 * _S1608.differential_0 + r2_13 * s_diff_r6_T_1);
    float _S1615 = _S1237 * _S1608.differential_0 + r4_7 * s_diff_r6_T_1 + _S1614 + _S1614;
    float _S1616 = dy_7 * _S1615;
    float _S1617 = dx_9 * _S1615;
    float _S1618 = - (_S1616 + _S1616);
    float _S1619 = - (_S1617 + _S1617);
    *&((&_S1263)->x) = 0.0f;
    float _S1620 = rgb_out_4.x * _S1540;
    float _S1621 = _S1261 * _S1540;
    DiffPair_float_0 _S1622;
    (&_S1622)->primal_0 = falloff_0;
    (&_S1622)->differential_0 = 0.0f;
    DiffPair_float_0 _S1623;
    (&_S1623)->primal_0 = 0.0f;
    (&_S1623)->differential_0 = 0.0f;
    DiffPair_float_0 _S1624;
    (&_S1624)->primal_0 = 1.0f;
    (&_S1624)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1622, &_S1623, &_S1624, _S1620);
    float _S1625 = r2_12 * _S1622.differential_0;
    float _S1626 = r4_6 * _S1622.differential_0;
    float s_diff_r6_T_2 = _S1236 * _S1622.differential_0;
    float _S1627 = r6_0 * _S1622.differential_0;
    float _S1628 = r2_12 * (_S1235 * _S1622.differential_0 + r2_12 * s_diff_r6_T_2);
    float _S1629 = _S1234 * _S1622.differential_0 + r4_6 * s_diff_r6_T_2 + _S1628 + _S1628;
    float _S1630 = dy_6 * _S1629;
    float _S1631 = dx_8 * _S1629;
    float _S1632 = - (_S1630 + _S1630);
    float _S1633 = - (_S1631 + _S1631);
    float3  _S1634 = _S1225;
    *&((&_S1634)->z) = _S1593;
    *&((&_S1634)->y) = _S1607;
    *&((&_S1634)->x) = _S1621;
    float3  _S1635 = _S1263 + _S1634;
    float3  _S1636 = _S1224.primal_0 * _S1635;
    float3  _S1637 = _S1257 * _S1635;
    float _S1638 = _S1636.x + _S1636.y + _S1636.z;
    DiffPair_float_0 _S1639;
    (&_S1639)->primal_0 = _S1255.exposure_1;
    (&_S1639)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1639, _S1638);
    PPISPParams_0 _S1640 = PPISPParams_x24_syn_dzero_0();
    (&_S1640)->color_params_1 = _S1591;
    (&_S1640)->exposure_1 = _S1639.differential_0;
    _S1233 = _S1640;
    (&(&_S1233)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1641 = _S1640.crf_params_1[int(2)].center_0 + _S1426;
    (&(&_S1233)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1642 = _S1640.crf_params_1[int(2)].gamma_0 + _S1429.differential_0;
    (&(&_S1233)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1643 = _S1640.crf_params_1[int(2)].shoulder_0 + _S1433.differential_0;
    (&(&_S1233)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1644 = _S1640.crf_params_1[int(2)].toe_0 + _S1437.differential_0;
    (&(&_S1233)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1645 = _S1640.crf_params_1[int(1)].center_0 + _S1470;
    (&(&_S1233)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1646 = _S1640.crf_params_1[int(1)].gamma_0 + _S1473.differential_0;
    (&(&_S1233)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1647 = _S1640.crf_params_1[int(1)].shoulder_0 + _S1477.differential_0;
    (&(&_S1233)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1648 = _S1640.crf_params_1[int(1)].toe_0 + _S1481.differential_0;
    (&(&_S1233)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1649 = _S1640.crf_params_1[int(0)].center_0 + _S1514;
    (&(&_S1233)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1650 = _S1640.crf_params_1[int(0)].gamma_0 + _S1517.differential_0;
    (&(&_S1233)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1651 = _S1640.crf_params_1[int(0)].shoulder_0 + _S1521.differential_0;
    (&(&_S1233)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1652 = _S1640.crf_params_1[int(0)].toe_0 + _S1525.differential_0;
    *&((&(&(&_S1233)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1233)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1233)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1653 = _S1599 + _S1640.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1233)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1654 = _S1598 + _S1640.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1233)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1655 = _S1597 + _S1640.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1233)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1656 = _S1604 + _S1640.vignette_params_1[int(2)].cy_0;
    (&(&_S1233)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1657 = _S1605 + _S1640.vignette_params_1[int(2)].cx_0;
    (&(&_S1233)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1658 = _S1613 + _S1640.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1233)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1659 = _S1612 + _S1640.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1233)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1660 = _S1611 + _S1640.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1233)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1661 = _S1618 + _S1640.vignette_params_1[int(1)].cy_0;
    (&(&_S1233)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1662 = _S1619 + _S1640.vignette_params_1[int(1)].cx_0;
    (&(&_S1233)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1663 = _S1627 + _S1640.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1233)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1664 = _S1626 + _S1640.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1233)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1665 = _S1625 + _S1640.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1233)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1666 = _S1632 + _S1640.vignette_params_1[int(0)].cy_0;
    (&(&_S1233)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1667 = _S1633 + _S1640.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S1668;
    _S1668[int(0)] = 0.0f;
    _S1668[int(1)] = 0.0f;
    _S1668[int(2)] = 0.0f;
    _S1668[int(3)] = 0.0f;
    _S1668[int(4)] = 0.0f;
    _S1668[int(5)] = 0.0f;
    _S1668[int(6)] = 0.0f;
    _S1668[int(7)] = 0.0f;
    _S1668[int(8)] = 0.0f;
    _S1668[int(9)] = 0.0f;
    _S1668[int(10)] = 0.0f;
    _S1668[int(11)] = 0.0f;
    _S1668[int(12)] = 0.0f;
    _S1668[int(13)] = 0.0f;
    _S1668[int(14)] = 0.0f;
    _S1668[int(15)] = 0.0f;
    _S1668[int(16)] = 0.0f;
    _S1668[int(17)] = 0.0f;
    _S1668[int(18)] = 0.0f;
    _S1668[int(19)] = 0.0f;
    _S1668[int(20)] = 0.0f;
    _S1668[int(21)] = 0.0f;
    _S1668[int(22)] = 0.0f;
    _S1668[int(23)] = 0.0f;
    _S1668[int(24)] = 0.0f;
    _S1668[int(25)] = 0.0f;
    _S1668[int(26)] = 0.0f;
    _S1668[int(27)] = 0.0f;
    _S1668[int(28)] = 0.0f;
    _S1668[int(29)] = 0.0f;
    _S1668[int(30)] = 0.0f;
    _S1668[int(31)] = 0.0f;
    _S1668[int(32)] = 0.0f;
    _S1668[int(33)] = 0.0f;
    _S1668[int(34)] = 0.0f;
    _S1668[int(35)] = 0.0f;
    _S1668[int(8)] = _S1660;
    _S1668[int(16)] = _S1640.color_params_1.b_0.x;
    _S1668[int(15)] = _S1653;
    _S1668[int(14)] = _S1654;
    _S1668[int(13)] = _S1655;
    _S1668[int(12)] = _S1656;
    _S1668[int(11)] = _S1657;
    _S1668[int(10)] = _S1658;
    _S1668[int(9)] = _S1659;
    _S1668[int(17)] = _S1640.color_params_1.b_0.y;
    _S1668[int(7)] = _S1661;
    _S1668[int(6)] = _S1662;
    _S1668[int(5)] = _S1663;
    _S1668[int(4)] = _S1664;
    _S1668[int(3)] = _S1665;
    _S1668[int(2)] = _S1666;
    _S1668[int(1)] = _S1667;
    _S1668[int(0)] = _S1233.exposure_1;
    _S1668[int(26)] = _S1650;
    _S1668[int(34)] = _S1642;
    _S1668[int(33)] = _S1643;
    _S1668[int(32)] = _S1644;
    _S1668[int(31)] = _S1645;
    _S1668[int(30)] = _S1646;
    _S1668[int(29)] = _S1647;
    _S1668[int(28)] = _S1648;
    _S1668[int(27)] = _S1649;
    _S1668[int(35)] = _S1641;
    _S1668[int(25)] = _S1651;
    _S1668[int(24)] = _S1652;
    _S1668[int(23)] = _S1640.color_params_1.n_0.y;
    _S1668[int(22)] = _S1640.color_params_1.n_0.x;
    _S1668[int(21)] = _S1640.color_params_1.g_0.y;
    _S1668[int(20)] = _S1640.color_params_1.g_0.x;
    _S1668[int(19)] = _S1640.color_params_1.r_0.y;
    _S1668[int(18)] = _S1640.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1668;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1637;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1669, float2  _S1670, float2  _S1671, float2  _S1672, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1673, float3  _S1674)
{
    s_bwd_prop_apply_ppisp_0(_S1669, _S1670, _S1671, _S1672, _S1673, _S1674);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1675 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S1675;
    FixedArray<float, 36>  _S1676 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_2;
    (&dp_params_0)->differential_0 = _S1676;
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

inline __device__ void s_bwd_prop_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_1, float2  pix_coord_4, float2  image_center_4, float2  img_size_4, DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_1, float3  _s_dOut_12)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1677 = *dprgb_in_1;
    float3  _S1678 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1679 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1680 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1681 = {
        _S1680, _S1680, _S1680
    };
    float2  _S1682 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1683 = { _S1682, _S1682, _S1682, _S1682 };
    RQSCRFPPISPChannelParams_0 _S1684 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1685 = {
        _S1684, _S1684, _S1684
    };
    PPISPParamsRQS_0 _S1686;
    (&_S1686)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1686)->vignette_params_0 = _S1681;
    (&_S1686)->color_params_0 = _S1683;
    (&_S1686)->crf_params_0 = _S1685;
    (&(&_S1686)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1686)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S1687 = dpparams_1->primal_0[int(3)];
    (&(&_S1686)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S1688 = dpparams_1->primal_0[int(4)];
    (&(&_S1686)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S1689 = dpparams_1->primal_0[int(5)];
    (&(&_S1686)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1686)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1686)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S1690 = dpparams_1->primal_0[int(8)];
    (&(&_S1686)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S1691 = dpparams_1->primal_0[int(9)];
    (&(&_S1686)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S1692 = dpparams_1->primal_0[int(10)];
    (&(&_S1686)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1686)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1686)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S1693 = dpparams_1->primal_0[int(13)];
    (&(&_S1686)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S1694 = dpparams_1->primal_0[int(14)];
    (&(&_S1686)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S1695 = dpparams_1->primal_0[int(15)];
    (&(&_S1686)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1686)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1686)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1686)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1686)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1686)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1686)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1686)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1686)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S1696 = dpparams_1->primal_0[int(24)];
    (&(&_S1686)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S1697 = dpparams_1->primal_0[int(25)];
    (&(&_S1686)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S1698 = dpparams_1->primal_0[int(26)];
    (&(&_S1686)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S1699 = dpparams_1->primal_0[int(27)];
    (&(&_S1686)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S1700 = dpparams_1->primal_0[int(28)];
    (&(&_S1686)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S1701 = dpparams_1->primal_0[int(29)];
    (&(&_S1686)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S1702 = dpparams_1->primal_0[int(30)];
    (&(&_S1686)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S1703 = dpparams_1->primal_0[int(31)];
    (&(&_S1686)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S1704 = dpparams_1->primal_0[int(32)];
    (&(&_S1686)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S1705 = dpparams_1->primal_0[int(33)];
    (&(&_S1686)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S1706 = dpparams_1->primal_0[int(34)];
    (&(&_S1686)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S1707 = dpparams_1->primal_0[int(35)];
    (&(&_S1686)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S1708 = dpparams_1->primal_0[int(36)];
    (&(&_S1686)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S1709 = dpparams_1->primal_0[int(37)];
    (&(&_S1686)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S1710 = dpparams_1->primal_0[int(38)];
    (&(&_S1686)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S1711 = _S1686;
    float _S1712 = s_primal_ctx_exp2_0(_S1686.exposure_0);
    float3  _S1713 = make_float3 (_S1712);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S1712);
    float _S1714 = (F32_max((img_size_4.x), (img_size_4.y)));
    float _S1715 = (pix_coord_4.x - image_center_4.x) / _S1714;
    float _S1716 = (pix_coord_4.y - image_center_4.y) / _S1714;
    float dx_11 = _S1715 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S1716 - dpparams_1->primal_0[int(2)];
    float r2_16 = dx_11 * dx_11 + dy_9 * dy_9;
    float r4_9 = r2_16 * r2_16;
    float r6_3 = r4_9 * r2_16;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_16 + 1.0f;
    float _S1717 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S1718 = rgb_out_5.x * _S1717;
    float3  _S1719 = rgb_out_5;
    *&((&_S1719)->x) = _S1718;
    float dx_12 = _S1715 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S1716 - dpparams_1->primal_0[int(7)];
    float r2_17 = dx_12 * dx_12 + dy_10 * dy_10;
    float r4_10 = r2_17 * r2_17;
    float r6_4 = r4_10 * r2_17;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_17 + 1.0f;
    float _S1720 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S1719)->y) = rgb_out_5.y * _S1720;
    float dx_13 = _S1715 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S1716 - dpparams_1->primal_0[int(12)];
    float r2_18 = dx_13 * dx_13 + dy_11 * dy_11;
    float r4_11 = r2_18 * r2_18;
    float r6_5 = r4_11 * r2_18;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_18 + 1.0f;
    float _S1721 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S1719)->z) = rgb_out_5.z * _S1721;
    PPISPParamsRQS_0 _S1722 = _S1686;
    float2  _S1723 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1686.color_params_0.b_0);
    float2  _S1724 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1686.color_params_0.r_0);
    float2  _S1725 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1686.color_params_0.g_0);
    float2  _S1726 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1686.color_params_0.n_0);
    float _S1727 = 0.3333333432674408f + _S1726.x;
    float _S1728 = 0.3333333432674408f + _S1726.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S1723.x, 1.0f + _S1724.x, _S1725.x, _S1723.y, _S1724.y, 1.0f + _S1725.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1728, 1.0f, 0.0f, - _S1727, - _S1728, _S1727, 0.0f);
    Matrix<float, 3, 3>  _S1729 = s_primal_ctx_mul_1(skew_1, T_3);
    float3  r0_3 = make_float3 (_S1729.rows[int(0)].x, _S1729.rows[int(0)].y, _S1729.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S1729.rows[int(1)].x, _S1729.rows[int(1)].y, _S1729.rows[int(1)].z);
    float3  r2_19 = make_float3 (_S1729.rows[int(2)].x, _S1729.rows[int(2)].y, _S1729.rows[int(2)].z);
    float3  _S1730 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S1731 = (s_primal_ctx_dot_0(_S1730, _S1730)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S1732;
    bool _S1733;
    if(_S1731)
    {
        float3  _S1734 = s_primal_ctx_cross_0(r0_3, r2_19);
        bool _S1735 = (s_primal_ctx_dot_0(_S1734, _S1734)) < 9.99999968265522539e-21f;
        if(_S1735)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_19);
        }
        else
        {
            lambda_v_7 = _S1734;
        }
        _S1733 = _S1735;
        _S1732 = _S1734;
    }
    else
    {
        lambda_v_7 = _S1730;
        _S1733 = false;
        _S1732 = _S1678;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S1736 = s_primal_ctx_mul_1(T_3, D_1);
    Matrix<float, 3, 3>  _S1737 = s_primal_ctx_mul_1(_S1736, S_inv_1);
    bool _S1738 = (s_primal_ctx_abs_0(_S1737.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S1739;
    float _S1740;
    if(_S1738)
    {
        float inv_s_1 = 1.0f / _S1737.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1741 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S1742 = _S1737.rows[int(2)].z * _S1737.rows[int(2)].z;
        H_5 = _S1737 * makeMatrix<float, 3, 3> (inv_s_1);
        _S1739 = _S1741;
        _S1740 = _S1742;
    }
    else
    {
        H_5 = _S1737;
        _S1739 = _S1679;
        _S1740 = 0.0f;
    }
    float _S1743 = _S1719.x;
    float _S1744 = _S1719.y;
    float intensity_3 = _S1743 + _S1744 + _S1719.z;
    float3  rgi_in_1 = make_float3 (_S1743, _S1744, intensity_3);
    float3  _S1745 = s_primal_ctx_mul_2(H_5, rgi_in_1);
    float _S1746 = _S1745.z + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S1746;
    float3  _S1747 = make_float3 (norm_factor_1);
    float _S1748 = _S1746 * _S1746;
    float3  rgi_out_5 = _S1745 * make_float3 (norm_factor_1);
    float _S1749 = rgi_out_5.x;
    float _S1750 = rgi_out_5.y;
    float3  _S1751 = make_float3 (_S1749, _S1750, rgi_out_5.z - _S1749 - _S1750);
    float3  _S1752 = make_float3 (0.0f);
    float3  _S1753 = make_float3 (1.0f);
    float3  _S1754 = s_primal_ctx_clamp_1(_S1751, _S1752, _S1753);
    float _S1755 = _S1754.x;
    float _S1756 = s_primal_ctx_exp_0(_S1696);
    float _S1757 = s_primal_ctx_exp_0(_S1697);
    float _S1758 = - _S1698;
    float _S1759 = 1.0f + s_primal_ctx_exp_0(_S1758);
    float x0_4 = 1.0f / _S1759;
    float _S1760 = _S1759 * _S1759;
    float _S1761 = - _S1699;
    float _S1762 = 1.0f + s_primal_ctx_exp_0(_S1761);
    float y0_4 = 1.0f / _S1762;
    float _S1763 = _S1762 * _S1762;
    float _S1764 = s_primal_ctx_exp_0(_S1700);
    bool _S1765 = _S1755 < x0_4;
    float _S1766;
    float _S1767;
    float _S1768;
    float _S1769;
    float _S1770;
    float _S1771;
    float _S1772;
    float _S1773;
    float _S1774;
    float _S1775;
    float _S1776;
    float _S1777;
    float _S1778;
    float _S1779;
    float _S1780;
    float _S1781;
    float _S1782;
    float _S1783;
    float _S1784;
    float _S1785;
    float _S1786;
    float _S1787;
    float _S1788;
    float _S1789;
    float _S1790;
    float _S1791;
    float _S1792;
    if(_S1765)
    {
        float s0_3 = y0_4 / x0_4;
        float _S1793 = x0_4 * x0_4;
        float t0_3 = _S1755 / x0_4;
        float _S1794 = s0_3 * t0_3;
        float _S1795 = _S1756 * t0_3;
        float _S1796 = 1.0f - t0_3;
        float _S1797 = _S1794 * t0_3 + _S1795 * _S1796;
        float _S1798 = y0_4 * _S1797;
        float _S1799 = _S1756 + _S1764 - 2.0f * s0_3;
        float _S1800 = _S1799 * t0_3;
        float _S1801 = s0_3 + _S1800 * _S1796;
        _S1766 = _S1801 * _S1801;
        _S1767 = _S1798;
        _S1768 = _S1801;
        _S1769 = _S1800;
        _S1770 = _S1796;
        _S1771 = _S1799;
        _S1772 = t0_3;
        _S1773 = _S1797;
        _S1774 = _S1795;
        _S1775 = _S1794;
        _S1776 = s0_3;
        _S1777 = _S1793;
        _S1778 = 0.0f;
        _S1779 = 0.0f;
        _S1780 = 0.0f;
        _S1781 = 0.0f;
        _S1782 = 0.0f;
        _S1783 = 0.0f;
        _S1784 = 0.0f;
        _S1785 = 0.0f;
        _S1786 = 0.0f;
        _S1787 = 0.0f;
        _S1788 = 0.0f;
        _S1789 = 0.0f;
        _S1790 = 0.0f;
        _S1791 = 0.0f;
        _S1792 = 0.0f;
    }
    else
    {
        float _S1802 = 1.0f - y0_4;
        float _S1803 = 1.0f - x0_4;
        float s1_7 = _S1802 / _S1803;
        float _S1804 = _S1803 * _S1803;
        float _S1805 = _S1755 - x0_4;
        float t1_3 = _S1805 / _S1803;
        float _S1806 = s1_7 * t1_3;
        float _S1807 = _S1764 * t1_3;
        float _S1808 = 1.0f - t1_3;
        float _S1809 = _S1806 * t1_3 + _S1807 * _S1808;
        float _S1810 = _S1802 * _S1809;
        float _S1811 = _S1764 + _S1757 - 2.0f * s1_7;
        float _S1812 = _S1811 * t1_3;
        float _S1813 = s1_7 + _S1812 * _S1808;
        float _S1814 = _S1813 * _S1813;
        _S1766 = 0.0f;
        _S1767 = 0.0f;
        _S1768 = 0.0f;
        _S1769 = 0.0f;
        _S1770 = 0.0f;
        _S1771 = 0.0f;
        _S1772 = 0.0f;
        _S1773 = 0.0f;
        _S1774 = 0.0f;
        _S1775 = 0.0f;
        _S1776 = 0.0f;
        _S1777 = 0.0f;
        _S1778 = _S1814;
        _S1779 = _S1810;
        _S1780 = _S1813;
        _S1781 = _S1812;
        _S1782 = _S1808;
        _S1783 = _S1811;
        _S1784 = t1_3;
        _S1785 = _S1802;
        _S1786 = _S1809;
        _S1787 = _S1807;
        _S1788 = _S1806;
        _S1789 = s1_7;
        _S1790 = _S1804;
        _S1791 = _S1805;
        _S1792 = _S1803;
    }
    float _S1815 = _S1754.y;
    float _S1816 = s_primal_ctx_exp_0(_S1701);
    float _S1817 = s_primal_ctx_exp_0(_S1702);
    float _S1818 = - _S1703;
    float _S1819 = 1.0f + s_primal_ctx_exp_0(_S1818);
    float x0_5 = 1.0f / _S1819;
    float _S1820 = _S1819 * _S1819;
    float _S1821 = - _S1704;
    float _S1822 = 1.0f + s_primal_ctx_exp_0(_S1821);
    float y0_5 = 1.0f / _S1822;
    float _S1823 = _S1822 * _S1822;
    float _S1824 = s_primal_ctx_exp_0(_S1705);
    bool _S1825 = _S1815 < x0_5;
    float _S1826;
    float _S1827;
    float _S1828;
    float _S1829;
    float _S1830;
    float _S1831;
    float _S1832;
    float _S1833;
    float _S1834;
    float _S1835;
    float _S1836;
    float _S1837;
    float _S1838;
    float _S1839;
    float _S1840;
    float _S1841;
    float _S1842;
    float _S1843;
    float _S1844;
    float _S1845;
    float _S1846;
    float _S1847;
    float _S1848;
    float _S1849;
    float _S1850;
    float _S1851;
    float _S1852;
    if(_S1825)
    {
        float s0_4 = y0_5 / x0_5;
        float _S1853 = x0_5 * x0_5;
        float t0_4 = _S1815 / x0_5;
        float _S1854 = s0_4 * t0_4;
        float _S1855 = _S1816 * t0_4;
        float _S1856 = 1.0f - t0_4;
        float _S1857 = _S1854 * t0_4 + _S1855 * _S1856;
        float _S1858 = y0_5 * _S1857;
        float _S1859 = _S1816 + _S1824 - 2.0f * s0_4;
        float _S1860 = _S1859 * t0_4;
        float _S1861 = s0_4 + _S1860 * _S1856;
        _S1826 = _S1861 * _S1861;
        _S1827 = _S1858;
        _S1828 = _S1861;
        _S1829 = _S1860;
        _S1830 = _S1856;
        _S1831 = _S1859;
        _S1832 = t0_4;
        _S1833 = _S1857;
        _S1834 = _S1855;
        _S1835 = _S1854;
        _S1836 = s0_4;
        _S1837 = _S1853;
        _S1838 = 0.0f;
        _S1839 = 0.0f;
        _S1840 = 0.0f;
        _S1841 = 0.0f;
        _S1842 = 0.0f;
        _S1843 = 0.0f;
        _S1844 = 0.0f;
        _S1845 = 0.0f;
        _S1846 = 0.0f;
        _S1847 = 0.0f;
        _S1848 = 0.0f;
        _S1849 = 0.0f;
        _S1850 = 0.0f;
        _S1851 = 0.0f;
        _S1852 = 0.0f;
    }
    else
    {
        float _S1862 = 1.0f - y0_5;
        float _S1863 = 1.0f - x0_5;
        float s1_8 = _S1862 / _S1863;
        float _S1864 = _S1863 * _S1863;
        float _S1865 = _S1815 - x0_5;
        float t1_4 = _S1865 / _S1863;
        float _S1866 = s1_8 * t1_4;
        float _S1867 = _S1824 * t1_4;
        float _S1868 = 1.0f - t1_4;
        float _S1869 = _S1866 * t1_4 + _S1867 * _S1868;
        float _S1870 = _S1862 * _S1869;
        float _S1871 = _S1824 + _S1817 - 2.0f * s1_8;
        float _S1872 = _S1871 * t1_4;
        float _S1873 = s1_8 + _S1872 * _S1868;
        float _S1874 = _S1873 * _S1873;
        _S1826 = 0.0f;
        _S1827 = 0.0f;
        _S1828 = 0.0f;
        _S1829 = 0.0f;
        _S1830 = 0.0f;
        _S1831 = 0.0f;
        _S1832 = 0.0f;
        _S1833 = 0.0f;
        _S1834 = 0.0f;
        _S1835 = 0.0f;
        _S1836 = 0.0f;
        _S1837 = 0.0f;
        _S1838 = _S1874;
        _S1839 = _S1870;
        _S1840 = _S1873;
        _S1841 = _S1872;
        _S1842 = _S1868;
        _S1843 = _S1871;
        _S1844 = t1_4;
        _S1845 = _S1862;
        _S1846 = _S1869;
        _S1847 = _S1867;
        _S1848 = _S1866;
        _S1849 = s1_8;
        _S1850 = _S1864;
        _S1851 = _S1865;
        _S1852 = _S1863;
    }
    float _S1875 = _S1754.z;
    float _S1876 = s_primal_ctx_exp_0(_S1706);
    float _S1877 = s_primal_ctx_exp_0(_S1707);
    float _S1878 = - _S1708;
    float _S1879 = 1.0f + s_primal_ctx_exp_0(_S1878);
    float x0_6 = 1.0f / _S1879;
    float _S1880 = _S1879 * _S1879;
    float _S1881 = - _S1709;
    float _S1882 = 1.0f + s_primal_ctx_exp_0(_S1881);
    float y0_6 = 1.0f / _S1882;
    float _S1883 = _S1882 * _S1882;
    float _S1884 = s_primal_ctx_exp_0(_S1710);
    bool _S1885 = _S1875 < x0_6;
    float _S1886;
    float _S1887;
    float _S1888;
    float _S1889;
    float _S1890;
    float _S1891;
    float _S1892;
    float _S1893;
    float _S1894;
    float _S1895;
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
    if(_S1885)
    {
        float s0_5 = y0_6 / x0_6;
        float _S1913 = x0_6 * x0_6;
        float t0_5 = _S1875 / x0_6;
        float _S1914 = s0_5 * t0_5;
        float _S1915 = _S1876 * t0_5;
        float _S1916 = 1.0f - t0_5;
        float _S1917 = _S1914 * t0_5 + _S1915 * _S1916;
        float _S1918 = y0_6 * _S1917;
        float _S1919 = _S1876 + _S1884 - 2.0f * s0_5;
        float _S1920 = _S1919 * t0_5;
        float _S1921 = s0_5 + _S1920 * _S1916;
        _S1886 = _S1921 * _S1921;
        _S1887 = _S1918;
        _S1888 = _S1921;
        _S1889 = _S1920;
        _S1890 = _S1916;
        _S1891 = _S1919;
        _S1892 = t0_5;
        _S1893 = _S1917;
        _S1894 = _S1915;
        _S1895 = _S1914;
        _S1896 = s0_5;
        _S1897 = _S1913;
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
        _S1908 = 0.0f;
        _S1909 = 0.0f;
        _S1910 = 0.0f;
        _S1911 = 0.0f;
        _S1912 = 0.0f;
    }
    else
    {
        float _S1922 = 1.0f - y0_6;
        float _S1923 = 1.0f - x0_6;
        float s1_9 = _S1922 / _S1923;
        float _S1924 = _S1923 * _S1923;
        float _S1925 = _S1875 - x0_6;
        float t1_5 = _S1925 / _S1923;
        float _S1926 = s1_9 * t1_5;
        float _S1927 = _S1884 * t1_5;
        float _S1928 = 1.0f - t1_5;
        float _S1929 = _S1926 * t1_5 + _S1927 * _S1928;
        float _S1930 = _S1922 * _S1929;
        float _S1931 = _S1884 + _S1877 - 2.0f * s1_9;
        float _S1932 = _S1931 * t1_5;
        float _S1933 = s1_9 + _S1932 * _S1928;
        float _S1934 = _S1933 * _S1933;
        _S1886 = 0.0f;
        _S1887 = 0.0f;
        _S1888 = 0.0f;
        _S1889 = 0.0f;
        _S1890 = 0.0f;
        _S1891 = 0.0f;
        _S1892 = 0.0f;
        _S1893 = 0.0f;
        _S1894 = 0.0f;
        _S1895 = 0.0f;
        _S1896 = 0.0f;
        _S1897 = 0.0f;
        _S1898 = _S1934;
        _S1899 = _S1930;
        _S1900 = _S1933;
        _S1901 = _S1932;
        _S1902 = _S1928;
        _S1903 = _S1931;
        _S1904 = t1_5;
        _S1905 = _S1922;
        _S1906 = _S1929;
        _S1907 = _S1927;
        _S1908 = _S1926;
        _S1909 = s1_9;
        _S1910 = _S1924;
        _S1911 = _S1925;
        _S1912 = _S1923;
    }
    if(_S1885)
    {
        float _S1935 = _s_dOut_12.z / _S1886;
        float _S1936 = _S1887 * - _S1935;
        float _S1937 = _S1888 * _S1935;
        float _S1938 = _S1890 * _S1936;
        float _S1939 = _S1892 * _S1938;
        float _S1940 = y0_6 * _S1937;
        float _S1941 = _S1890 * _S1940;
        float _S1942 = _S1892 * _S1940;
        float _S1943 = (_S1891 * _S1938 + - (_S1889 * _S1936 + _S1894 * _S1940) + _S1876 * _S1941 + _S1895 * _S1940 + _S1896 * _S1942) / _S1897;
        float _S1944 = x0_6 * _S1943;
        float _S1945 = (_S1936 + 2.0f * - _S1939 + _S1892 * _S1942) / _S1897;
        float _S1946 = _S1893 * _S1937 + x0_6 * _S1945;
        float _S1947 = _S1939 + _S1892 * _S1941;
        float _S1948 = _S1875 * - _S1943 + y0_6 * - _S1945;
        _S1886 = _S1939;
        _S1887 = _S1946;
        _S1888 = _S1948;
        _S1889 = 0.0f;
        _S1890 = _S1947;
        _S1891 = _S1944;
    }
    else
    {
        float _S1949 = _s_dOut_12.z / _S1898;
        float _S1950 = _S1899 * - _S1949;
        float _S1951 = _S1900 * _S1949;
        float _S1952 = _S1902 * _S1950;
        float _S1953 = _S1904 * _S1952;
        float _S1954 = _S1905 * _S1951;
        float _S1955 = _S1902 * _S1954;
        float _S1956 = _S1904 * _S1954;
        float _S1957 = (_S1903 * _S1952 + - (_S1901 * _S1950 + _S1907 * _S1954) + _S1884 * _S1955 + _S1908 * _S1954 + _S1909 * _S1956) / _S1910;
        float _S1958 = _S1912 * _S1957;
        float _S1959 = (_S1950 + 2.0f * - _S1953 + _S1904 * _S1956) / _S1910;
        float _S1960 = _s_dOut_12.z + - (_S1906 * _S1951 + _S1912 * _S1959);
        float _S1961 = - _S1958 + - (_S1911 * - _S1957 + _S1905 * - _S1959);
        _S1886 = _S1953 + _S1904 * _S1955;
        _S1887 = _S1960;
        _S1888 = _S1961;
        _S1889 = _S1953;
        _S1890 = 0.0f;
        _S1891 = _S1958;
    }
    DiffPair_float_0 _S1962;
    (&_S1962)->primal_0 = _S1710;
    (&_S1962)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1962, _S1886);
    DiffPair_float_0 _S1963 = _S1962;
    float _S1964 = - (_S1887 / _S1883);
    DiffPair_float_0 _S1965;
    (&_S1965)->primal_0 = _S1881;
    (&_S1965)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1965, _S1964);
    float _S1966 = - _S1965.differential_0;
    float _S1967 = - (_S1888 / _S1880);
    DiffPair_float_0 _S1968;
    (&_S1968)->primal_0 = _S1878;
    (&_S1968)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1968, _S1967);
    float _S1969 = - _S1968.differential_0;
    DiffPair_float_0 _S1970;
    (&_S1970)->primal_0 = _S1707;
    (&_S1970)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1970, _S1889);
    DiffPair_float_0 _S1971 = _S1970;
    DiffPair_float_0 _S1972;
    (&_S1972)->primal_0 = _S1706;
    (&_S1972)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1972, _S1890);
    DiffPair_float_0 _S1973 = _S1972;
    float3  _S1974 = make_float3 (0.0f, 0.0f, _S1891);
    if(_S1825)
    {
        float _S1975 = _s_dOut_12.y / _S1826;
        float _S1976 = _S1827 * - _S1975;
        float _S1977 = _S1828 * _S1975;
        float _S1978 = _S1830 * _S1976;
        float _S1979 = _S1832 * _S1978;
        float _S1980 = y0_5 * _S1977;
        float _S1981 = _S1830 * _S1980;
        float _S1982 = _S1832 * _S1980;
        float _S1983 = (_S1831 * _S1978 + - (_S1829 * _S1976 + _S1834 * _S1980) + _S1816 * _S1981 + _S1835 * _S1980 + _S1836 * _S1982) / _S1837;
        float _S1984 = x0_5 * _S1983;
        float _S1985 = (_S1976 + 2.0f * - _S1979 + _S1832 * _S1982) / _S1837;
        float _S1986 = _S1833 * _S1977 + x0_5 * _S1985;
        float _S1987 = _S1979 + _S1832 * _S1981;
        float _S1988 = _S1815 * - _S1983 + y0_5 * - _S1985;
        _S1826 = _S1979;
        _S1827 = _S1986;
        _S1828 = _S1988;
        _S1829 = 0.0f;
        _S1830 = _S1987;
        _S1831 = _S1984;
    }
    else
    {
        float _S1989 = _s_dOut_12.y / _S1838;
        float _S1990 = _S1839 * - _S1989;
        float _S1991 = _S1840 * _S1989;
        float _S1992 = _S1842 * _S1990;
        float _S1993 = _S1844 * _S1992;
        float _S1994 = _S1845 * _S1991;
        float _S1995 = _S1842 * _S1994;
        float _S1996 = _S1844 * _S1994;
        float _S1997 = (_S1843 * _S1992 + - (_S1841 * _S1990 + _S1847 * _S1994) + _S1824 * _S1995 + _S1848 * _S1994 + _S1849 * _S1996) / _S1850;
        float _S1998 = _S1852 * _S1997;
        float _S1999 = (_S1990 + 2.0f * - _S1993 + _S1844 * _S1996) / _S1850;
        float _S2000 = _s_dOut_12.y + - (_S1846 * _S1991 + _S1852 * _S1999);
        float _S2001 = - _S1998 + - (_S1851 * - _S1997 + _S1845 * - _S1999);
        _S1826 = _S1993 + _S1844 * _S1995;
        _S1827 = _S2000;
        _S1828 = _S2001;
        _S1829 = _S1993;
        _S1830 = 0.0f;
        _S1831 = _S1998;
    }
    DiffPair_float_0 _S2002;
    (&_S2002)->primal_0 = _S1705;
    (&_S2002)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2002, _S1826);
    DiffPair_float_0 _S2003 = _S2002;
    float _S2004 = - (_S1827 / _S1823);
    DiffPair_float_0 _S2005;
    (&_S2005)->primal_0 = _S1821;
    (&_S2005)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2005, _S2004);
    float _S2006 = - _S2005.differential_0;
    float _S2007 = - (_S1828 / _S1820);
    DiffPair_float_0 _S2008;
    (&_S2008)->primal_0 = _S1818;
    (&_S2008)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2008, _S2007);
    float _S2009 = - _S2008.differential_0;
    DiffPair_float_0 _S2010;
    (&_S2010)->primal_0 = _S1702;
    (&_S2010)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2010, _S1829);
    DiffPair_float_0 _S2011 = _S2010;
    DiffPair_float_0 _S2012;
    (&_S2012)->primal_0 = _S1701;
    (&_S2012)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2012, _S1830);
    DiffPair_float_0 _S2013 = _S2012;
    float3  _S2014 = _S1974 + make_float3 (0.0f, _S1831, 0.0f);
    if(_S1765)
    {
        float _S2015 = _s_dOut_12.x / _S1766;
        float _S2016 = _S1767 * - _S2015;
        float _S2017 = _S1768 * _S2015;
        float _S2018 = _S1770 * _S2016;
        float _S2019 = _S1772 * _S2018;
        float _S2020 = y0_4 * _S2017;
        float _S2021 = _S1770 * _S2020;
        float _S2022 = _S1772 * _S2020;
        float _S2023 = (_S1771 * _S2018 + - (_S1769 * _S2016 + _S1774 * _S2020) + _S1756 * _S2021 + _S1775 * _S2020 + _S1776 * _S2022) / _S1777;
        float _S2024 = x0_4 * _S2023;
        float _S2025 = (_S2016 + 2.0f * - _S2019 + _S1772 * _S2022) / _S1777;
        float _S2026 = _S1773 * _S2017 + x0_4 * _S2025;
        float _S2027 = _S2019 + _S1772 * _S2021;
        float _S2028 = _S1755 * - _S2023 + y0_4 * - _S2025;
        _S1766 = _S2019;
        _S1767 = _S2026;
        _S1768 = _S2028;
        _S1769 = 0.0f;
        _S1770 = _S2027;
        _S1771 = _S2024;
    }
    else
    {
        float _S2029 = _s_dOut_12.x / _S1778;
        float _S2030 = _S1779 * - _S2029;
        float _S2031 = _S1780 * _S2029;
        float _S2032 = _S1782 * _S2030;
        float _S2033 = _S1784 * _S2032;
        float _S2034 = _S1785 * _S2031;
        float _S2035 = _S1782 * _S2034;
        float _S2036 = _S1784 * _S2034;
        float _S2037 = (_S1783 * _S2032 + - (_S1781 * _S2030 + _S1787 * _S2034) + _S1764 * _S2035 + _S1788 * _S2034 + _S1789 * _S2036) / _S1790;
        float _S2038 = _S1792 * _S2037;
        float _S2039 = (_S2030 + 2.0f * - _S2033 + _S1784 * _S2036) / _S1790;
        float _S2040 = _s_dOut_12.x + - (_S1786 * _S2031 + _S1792 * _S2039);
        float _S2041 = - _S2038 + - (_S1791 * - _S2037 + _S1785 * - _S2039);
        _S1766 = _S2033 + _S1784 * _S2035;
        _S1767 = _S2040;
        _S1768 = _S2041;
        _S1769 = _S2033;
        _S1770 = 0.0f;
        _S1771 = _S2038;
    }
    DiffPair_float_0 _S2042;
    (&_S2042)->primal_0 = _S1700;
    (&_S2042)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2042, _S1766);
    DiffPair_float_0 _S2043 = _S2042;
    float _S2044 = - (_S1767 / _S1763);
    DiffPair_float_0 _S2045;
    (&_S2045)->primal_0 = _S1761;
    (&_S2045)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2045, _S2044);
    float _S2046 = - _S2045.differential_0;
    float _S2047 = - (_S1768 / _S1760);
    DiffPair_float_0 _S2048;
    (&_S2048)->primal_0 = _S1758;
    (&_S2048)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2048, _S2047);
    float _S2049 = - _S2048.differential_0;
    DiffPair_float_0 _S2050;
    (&_S2050)->primal_0 = _S1697;
    (&_S2050)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2050, _S1769);
    DiffPair_float_0 _S2051 = _S2050;
    DiffPair_float_0 _S2052;
    (&_S2052)->primal_0 = _S1696;
    (&_S2052)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2052, _S1770);
    DiffPair_float_0 _S2053 = _S2052;
    float3  _S2054 = _S2014 + make_float3 (_S1771, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2055;
    (&_S2055)->primal_0 = _S1751;
    (&_S2055)->differential_0 = _S1678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2056;
    (&_S2056)->primal_0 = _S1752;
    (&_S2056)->differential_0 = _S1678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2057;
    (&_S2057)->primal_0 = _S1753;
    (&_S2057)->differential_0 = _S1678;
    s_bwd_prop_clamp_1(&_S2055, &_S2056, &_S2057, _S2054);
    float _S2058 = - _S2055.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S2055.differential_0.x + _S2058, _S2055.differential_0.y + _S2058, _S2055.differential_0.z);
    float3  _S2059 = _S1745 * s_diff_rgi_out_T_1;
    float _S2060 = (_S2059.x + _S2059.y + _S2059.z) / _S1748;
    float _S2061 = _S1746 * _S2060;
    float3  _S2062 = _S1747 * s_diff_rgi_out_T_1 + make_float3 (0.0f, 0.0f, intensity_3 * - _S2060);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2063;
    (&_S2063)->primal_0 = H_5;
    (&_S2063)->differential_0 = _S1679;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2064;
    (&_S2064)->primal_0 = rgi_in_1;
    (&_S2064)->differential_0 = _S1678;
    s_bwd_prop_mul_1(&_S2063, &_S2064, _S2062);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2065 = _S2063;
    float _S2066 = _S2061 + _S2064.differential_0.z;
    float _S2067 = _S2064.differential_0.y + _S2066;
    float _S2068 = _S2064.differential_0.x + _S2066;
    float3  _S2069 = make_float3 (_S2068, _S2067, _S2066);
    if(_S1738)
    {
        Matrix<float, 3, 3>  _S2070 = _S1737 * _S2065.differential_0;
        Matrix<float, 3, 3>  _S2071 = _S1739 * _S2065.differential_0;
        _S1740 = - ((_S2070.rows[int(0)].x + _S2070.rows[int(0)].y + _S2070.rows[int(0)].z + _S2070.rows[int(1)].x + _S2070.rows[int(1)].y + _S2070.rows[int(1)].z + _S2070.rows[int(2)].x + _S2070.rows[int(2)].y + _S2070.rows[int(2)].z) / _S1740);
        H_5 = _S2071;
    }
    else
    {
        _S1740 = 0.0f;
        H_5 = _S2065.differential_0;
    }
    DiffPair_float_0 _S2072;
    (&_S2072)->primal_0 = _S1737.rows[int(2)].z;
    (&_S2072)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2072, 0.0f);
    float _S2073 = _S2072.differential_0 + _S1740;
    float3  _S2074 = _S1678;
    *&((&_S2074)->z) = _S2073;
    Matrix<float, 3, 3>  _S2075 = _S1679;
    _S2075[int(2)] = _S2074;
    Matrix<float, 3, 3>  _S2076 = H_5 + _S2075;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2077;
    (&_S2077)->primal_0 = _S1736;
    (&_S2077)->differential_0 = _S1679;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2078;
    (&_S2078)->primal_0 = S_inv_1;
    (&_S2078)->differential_0 = _S1679;
    s_bwd_prop_mul_2(&_S2077, &_S2078, _S2076);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2079;
    (&_S2079)->primal_0 = T_3;
    (&_S2079)->differential_0 = _S1679;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2080;
    (&_S2080)->primal_0 = D_1;
    (&_S2080)->differential_0 = _S1679;
    s_bwd_prop_mul_2(&_S2079, &_S2080, _S2077.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2081 = _S2079;
    float3  _S2082 = make_float3 (_S2080.differential_0.rows[int(0)].x, _S2080.differential_0.rows[int(1)].y, _S2080.differential_0.rows[int(2)].z);
    float3  _S2083;
    if(_S1731)
    {
        if(_S1733)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2084;
            (&_S2084)->primal_0 = r1_3;
            (&_S2084)->differential_0 = _S1678;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2085;
            (&_S2085)->primal_0 = r2_19;
            (&_S2085)->differential_0 = _S1678;
            s_bwd_prop_cross_0(&_S2084, &_S2085, _S2082);
            _S1719 = _S1678;
            lambda_v_7 = _S2085.differential_0;
            _S2083 = _S2084.differential_0;
        }
        else
        {
            _S1719 = _S2082;
            lambda_v_7 = _S1678;
            _S2083 = _S1678;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2086;
        (&_S2086)->primal_0 = _S1732;
        (&_S2086)->differential_0 = _S1678;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2087;
        (&_S2087)->primal_0 = _S1732;
        (&_S2087)->differential_0 = _S1678;
        s_bwd_prop_dot_0(&_S2086, &_S2087, 0.0f);
        float3  _S2088 = _S2087.differential_0 + _S2086.differential_0 + _S1719;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2089;
        (&_S2089)->primal_0 = r0_3;
        (&_S2089)->differential_0 = _S1678;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2090;
        (&_S2090)->primal_0 = r2_19;
        (&_S2090)->differential_0 = _S1678;
        s_bwd_prop_cross_0(&_S2089, &_S2090, _S2088);
        float3  _S2091 = _S2090.differential_0 + lambda_v_7;
        _S1719 = _S1678;
        lambda_v_7 = _S2091;
        _S1732 = _S2083;
        _S2083 = _S2089.differential_0;
    }
    else
    {
        _S1719 = _S2082;
        lambda_v_7 = _S1678;
        _S1732 = _S1678;
        _S2083 = _S1678;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2092;
    (&_S2092)->primal_0 = _S1730;
    (&_S2092)->differential_0 = _S1678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2093;
    (&_S2093)->primal_0 = _S1730;
    (&_S2093)->differential_0 = _S1678;
    s_bwd_prop_dot_0(&_S2092, &_S2093, 0.0f);
    float3  _S2094 = _S2093.differential_0 + _S2092.differential_0 + _S1719;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2095;
    (&_S2095)->primal_0 = r0_3;
    (&_S2095)->differential_0 = _S1678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2096;
    (&_S2096)->primal_0 = r1_3;
    (&_S2096)->differential_0 = _S1678;
    s_bwd_prop_cross_0(&_S2095, &_S2096, _S2094);
    float3  _S2097 = _S1678;
    *&((&_S2097)->z) = lambda_v_7.z;
    *&((&_S2097)->y) = lambda_v_7.y;
    *&((&_S2097)->x) = lambda_v_7.x;
    float3  _S2098 = _S2096.differential_0 + _S1732;
    float3  _S2099 = _S1678;
    *&((&_S2099)->z) = _S2098.z;
    *&((&_S2099)->y) = _S2098.y;
    *&((&_S2099)->x) = _S2098.x;
    float3  _S2100 = _S2095.differential_0 + _S2083;
    float3  _S2101 = _S1678;
    *&((&_S2101)->z) = _S2100.z;
    *&((&_S2101)->y) = _S2100.y;
    *&((&_S2101)->x) = _S2100.x;
    Matrix<float, 3, 3>  _S2102 = _S1679;
    _S2102[int(2)] = _S2097;
    _S2102[int(1)] = _S2099;
    _S2102[int(0)] = _S2101;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2103;
    (&_S2103)->primal_0 = skew_1;
    (&_S2103)->differential_0 = _S1679;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2104;
    (&_S2104)->primal_0 = T_3;
    (&_S2104)->differential_0 = _S1679;
    s_bwd_prop_mul_2(&_S2103, &_S2104, _S2102);
    Matrix<float, 3, 3>  _S2105 = _S2104.differential_0 + _S2081.differential_0;
    float2  _S2106 = make_float2 (_S2103.differential_0.rows[int(2)].y + - _S2103.differential_0.rows[int(1)].z, _S2103.differential_0.rows[int(0)].z + - _S2103.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S2107 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2108;
    (&_S2108)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2108)->differential_0 = _S2107;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2109;
    (&_S2109)->primal_0 = _S1722.color_params_0.n_0;
    (&_S2109)->differential_0 = _S1682;
    s_bwd_prop_mul_3(&_S2108, &_S2109, _S2106);
    float2  _S2110 = make_float2 (_S2105.rows[int(0)].z, _S2105.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2111;
    (&_S2111)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2111)->differential_0 = _S2107;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2112;
    (&_S2112)->primal_0 = _S1722.color_params_0.g_0;
    (&_S2112)->differential_0 = _S1682;
    s_bwd_prop_mul_3(&_S2111, &_S2112, _S2110);
    float2  _S2113 = make_float2 (_S2105.rows[int(0)].y, _S2105.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2114;
    (&_S2114)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2114)->differential_0 = _S2107;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2115;
    (&_S2115)->primal_0 = _S1722.color_params_0.r_0;
    (&_S2115)->differential_0 = _S1682;
    s_bwd_prop_mul_3(&_S2114, &_S2115, _S2113);
    float2  _S2116 = make_float2 (_S2105.rows[int(0)].x, _S2105.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2117;
    (&_S2117)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2117)->differential_0 = _S2107;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2118;
    (&_S2118)->primal_0 = _S1722.color_params_0.b_0;
    (&_S2118)->differential_0 = _S1682;
    s_bwd_prop_mul_3(&_S2117, &_S2118, _S2116);
    ColorPPISPParams_0 _S2119 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2119)->n_0 = _S2109.differential_0;
    (&_S2119)->g_0 = _S2112.differential_0;
    (&_S2119)->r_0 = _S2115.differential_0;
    (&_S2119)->b_0 = _S2118.differential_0;
    _S1719 = _S2069;
    *&((&_S1719)->z) = 0.0f;
    float _S2120 = rgb_out_5.z * _S2066;
    float _S2121 = _S1721 * _S2066;
    DiffPair_float_0 _S2122;
    (&_S2122)->primal_0 = falloff_5;
    (&_S2122)->differential_0 = 0.0f;
    DiffPair_float_0 _S2123;
    (&_S2123)->primal_0 = 0.0f;
    (&_S2123)->differential_0 = 0.0f;
    DiffPair_float_0 _S2124;
    (&_S2124)->primal_0 = 1.0f;
    (&_S2124)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2122, &_S2123, &_S2124, _S2120);
    float _S2125 = r2_18 * _S2122.differential_0;
    float _S2126 = r4_11 * _S2122.differential_0;
    float s_diff_r6_T_3 = _S1695 * _S2122.differential_0;
    float _S2127 = r6_5 * _S2122.differential_0;
    float _S2128 = r2_18 * (_S1694 * _S2122.differential_0 + r2_18 * s_diff_r6_T_3);
    float _S2129 = _S1693 * _S2122.differential_0 + r4_11 * s_diff_r6_T_3 + _S2128 + _S2128;
    float _S2130 = dy_11 * _S2129;
    float _S2131 = dx_13 * _S2129;
    float _S2132 = - (_S2130 + _S2130);
    float _S2133 = - (_S2131 + _S2131);
    *&((&_S1719)->y) = 0.0f;
    float _S2134 = rgb_out_5.y * _S2067;
    float _S2135 = _S1720 * _S2067;
    DiffPair_float_0 _S2136;
    (&_S2136)->primal_0 = falloff_4;
    (&_S2136)->differential_0 = 0.0f;
    DiffPair_float_0 _S2137;
    (&_S2137)->primal_0 = 0.0f;
    (&_S2137)->differential_0 = 0.0f;
    DiffPair_float_0 _S2138;
    (&_S2138)->primal_0 = 1.0f;
    (&_S2138)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2136, &_S2137, &_S2138, _S2134);
    float _S2139 = r2_17 * _S2136.differential_0;
    float _S2140 = r4_10 * _S2136.differential_0;
    float s_diff_r6_T_4 = _S1692 * _S2136.differential_0;
    float _S2141 = r6_4 * _S2136.differential_0;
    float _S2142 = r2_17 * (_S1691 * _S2136.differential_0 + r2_17 * s_diff_r6_T_4);
    float _S2143 = _S1690 * _S2136.differential_0 + r4_10 * s_diff_r6_T_4 + _S2142 + _S2142;
    float _S2144 = dy_10 * _S2143;
    float _S2145 = dx_12 * _S2143;
    float _S2146 = - (_S2144 + _S2144);
    float _S2147 = - (_S2145 + _S2145);
    *&((&_S1719)->x) = 0.0f;
    float _S2148 = rgb_out_5.x * _S2068;
    float _S2149 = _S1717 * _S2068;
    DiffPair_float_0 _S2150;
    (&_S2150)->primal_0 = falloff_3;
    (&_S2150)->differential_0 = 0.0f;
    DiffPair_float_0 _S2151;
    (&_S2151)->primal_0 = 0.0f;
    (&_S2151)->differential_0 = 0.0f;
    DiffPair_float_0 _S2152;
    (&_S2152)->primal_0 = 1.0f;
    (&_S2152)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2150, &_S2151, &_S2152, _S2148);
    float _S2153 = r2_16 * _S2150.differential_0;
    float _S2154 = r4_9 * _S2150.differential_0;
    float s_diff_r6_T_5 = _S1689 * _S2150.differential_0;
    float _S2155 = r6_3 * _S2150.differential_0;
    float _S2156 = r2_16 * (_S1688 * _S2150.differential_0 + r2_16 * s_diff_r6_T_5);
    float _S2157 = _S1687 * _S2150.differential_0 + r4_9 * s_diff_r6_T_5 + _S2156 + _S2156;
    float _S2158 = dy_9 * _S2157;
    float _S2159 = dx_11 * _S2157;
    float _S2160 = - (_S2158 + _S2158);
    float _S2161 = - (_S2159 + _S2159);
    float3  _S2162 = _S1678;
    *&((&_S2162)->z) = _S2121;
    *&((&_S2162)->y) = _S2135;
    *&((&_S2162)->x) = _S2149;
    float3  _S2163 = _S1719 + _S2162;
    float3  _S2164 = _S1677.primal_0 * _S2163;
    float3  _S2165 = _S1713 * _S2163;
    float _S2166 = _S2164.x + _S2164.y + _S2164.z;
    DiffPair_float_0 _S2167;
    (&_S2167)->primal_0 = _S1711.exposure_0;
    (&_S2167)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2167, _S2166);
    PPISPParamsRQS_0 _S2168 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2168)->color_params_0 = _S2119;
    (&_S2168)->exposure_0 = _S2167.differential_0;
    _S1686 = _S2168;
    (&(&_S1686)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2169 = _S2168.crf_params_0[int(2)].gc_0 + _S1963.differential_0;
    (&(&_S1686)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2170 = _S2168.crf_params_0[int(2)].y0_0 + _S1966;
    (&(&_S1686)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2171 = _S2168.crf_params_0[int(2)].x0_0 + _S1969;
    (&(&_S1686)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2172 = _S2168.crf_params_0[int(2)].g1_0 + _S1971.differential_0;
    (&(&_S1686)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2173 = _S2168.crf_params_0[int(2)].g0_0 + _S1973.differential_0;
    (&(&_S1686)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2174 = _S2168.crf_params_0[int(1)].gc_0 + _S2003.differential_0;
    (&(&_S1686)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2175 = _S2168.crf_params_0[int(1)].y0_0 + _S2006;
    (&(&_S1686)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2176 = _S2168.crf_params_0[int(1)].x0_0 + _S2009;
    (&(&_S1686)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2177 = _S2168.crf_params_0[int(1)].g1_0 + _S2011.differential_0;
    (&(&_S1686)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2178 = _S2168.crf_params_0[int(1)].g0_0 + _S2013.differential_0;
    (&(&_S1686)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2179 = _S2168.crf_params_0[int(0)].gc_0 + _S2043.differential_0;
    (&(&_S1686)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2180 = _S2168.crf_params_0[int(0)].y0_0 + _S2046;
    (&(&_S1686)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2181 = _S2168.crf_params_0[int(0)].x0_0 + _S2049;
    (&(&_S1686)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2182 = _S2168.crf_params_0[int(0)].g1_0 + _S2051.differential_0;
    (&(&_S1686)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2183 = _S2168.crf_params_0[int(0)].g0_0 + _S2053.differential_0;
    *&((&(&(&_S1686)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1686)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1686)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2184 = _S2127 + _S2168.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1686)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2185 = _S2126 + _S2168.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1686)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2186 = _S2125 + _S2168.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1686)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2187 = _S2132 + _S2168.vignette_params_0[int(2)].cy_0;
    (&(&_S1686)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2188 = _S2133 + _S2168.vignette_params_0[int(2)].cx_0;
    (&(&_S1686)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2189 = _S2141 + _S2168.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1686)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2190 = _S2140 + _S2168.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1686)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2191 = _S2139 + _S2168.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1686)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2192 = _S2146 + _S2168.vignette_params_0[int(1)].cy_0;
    (&(&_S1686)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2193 = _S2147 + _S2168.vignette_params_0[int(1)].cx_0;
    (&(&_S1686)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2194 = _S2155 + _S2168.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1686)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2195 = _S2154 + _S2168.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1686)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2196 = _S2153 + _S2168.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1686)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2197 = _S2160 + _S2168.vignette_params_0[int(0)].cy_0;
    (&(&_S1686)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2198 = _S2161 + _S2168.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2199;
    _S2199[int(0)] = 0.0f;
    _S2199[int(1)] = 0.0f;
    _S2199[int(2)] = 0.0f;
    _S2199[int(3)] = 0.0f;
    _S2199[int(4)] = 0.0f;
    _S2199[int(5)] = 0.0f;
    _S2199[int(6)] = 0.0f;
    _S2199[int(7)] = 0.0f;
    _S2199[int(8)] = 0.0f;
    _S2199[int(9)] = 0.0f;
    _S2199[int(10)] = 0.0f;
    _S2199[int(11)] = 0.0f;
    _S2199[int(12)] = 0.0f;
    _S2199[int(13)] = 0.0f;
    _S2199[int(14)] = 0.0f;
    _S2199[int(15)] = 0.0f;
    _S2199[int(16)] = 0.0f;
    _S2199[int(17)] = 0.0f;
    _S2199[int(18)] = 0.0f;
    _S2199[int(19)] = 0.0f;
    _S2199[int(20)] = 0.0f;
    _S2199[int(21)] = 0.0f;
    _S2199[int(22)] = 0.0f;
    _S2199[int(23)] = 0.0f;
    _S2199[int(24)] = 0.0f;
    _S2199[int(25)] = 0.0f;
    _S2199[int(26)] = 0.0f;
    _S2199[int(27)] = 0.0f;
    _S2199[int(28)] = 0.0f;
    _S2199[int(29)] = 0.0f;
    _S2199[int(30)] = 0.0f;
    _S2199[int(31)] = 0.0f;
    _S2199[int(32)] = 0.0f;
    _S2199[int(33)] = 0.0f;
    _S2199[int(34)] = 0.0f;
    _S2199[int(35)] = 0.0f;
    _S2199[int(36)] = 0.0f;
    _S2199[int(37)] = 0.0f;
    _S2199[int(38)] = 0.0f;
    _S2199[int(9)] = _S2190;
    _S2199[int(18)] = _S2168.color_params_0.r_0.x;
    _S2199[int(17)] = _S2168.color_params_0.b_0.y;
    _S2199[int(16)] = _S2168.color_params_0.b_0.x;
    _S2199[int(15)] = _S2184;
    _S2199[int(14)] = _S2185;
    _S2199[int(13)] = _S2186;
    _S2199[int(12)] = _S2187;
    _S2199[int(11)] = _S2188;
    _S2199[int(10)] = _S2189;
    _S2199[int(19)] = _S2168.color_params_0.r_0.y;
    _S2199[int(8)] = _S2191;
    _S2199[int(7)] = _S2192;
    _S2199[int(6)] = _S2193;
    _S2199[int(5)] = _S2194;
    _S2199[int(4)] = _S2195;
    _S2199[int(3)] = _S2196;
    _S2199[int(2)] = _S2197;
    _S2199[int(1)] = _S2198;
    _S2199[int(0)] = _S1686.exposure_0;
    _S2199[int(28)] = _S2179;
    _S2199[int(37)] = _S2170;
    _S2199[int(36)] = _S2171;
    _S2199[int(35)] = _S2172;
    _S2199[int(34)] = _S2173;
    _S2199[int(33)] = _S2174;
    _S2199[int(32)] = _S2175;
    _S2199[int(31)] = _S2176;
    _S2199[int(30)] = _S2177;
    _S2199[int(29)] = _S2178;
    _S2199[int(38)] = _S2169;
    _S2199[int(27)] = _S2180;
    _S2199[int(26)] = _S2181;
    _S2199[int(25)] = _S2182;
    _S2199[int(24)] = _S2183;
    _S2199[int(23)] = _S2168.color_params_0.n_0.y;
    _S2199[int(22)] = _S2168.color_params_0.n_0.x;
    _S2199[int(21)] = _S2168.color_params_0.g_0.y;
    _S2199[int(20)] = _S2168.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S2199;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S2165;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2200, float2  _S2201, float2  _S2202, float2  _S2203, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2204, float3  _S2205)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S2200, _S2201, _S2202, _S2203, _S2204, _S2205);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S2206 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S2206;
    FixedArray<float, 39>  _S2207 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_3;
    (&dp_params_1)->differential_0 = _S2207;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_4, FixedArray<float, 22>  * _S2208)
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
    float _S2209 = p_2.vignette_params_1[int(0)].cx_0;
    float _S2210 = p_2.vignette_params_1[int(0)].cy_0;
    float _S2211 = p_2.vignette_params_1[int(1)].cx_0;
    float _S2212 = p_2.vignette_params_1[int(1)].cy_0;
    float _S2213 = p_2.vignette_params_1[int(2)].cx_0;
    float _S2214 = p_2.vignette_params_1[int(2)].cy_0;
    losses_3[int(1)] = _S2209 * _S2209 + _S2210 * _S2210 + _S2211 * _S2211 + _S2212 * _S2212 + _S2213 * _S2213 + _S2214 * _S2214;
    losses_3[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_3[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_3[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_3 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S2215 = p_2.vignette_params_1[int(0)].cx_0 - mean_3;
    float _S2216 = p_2.vignette_params_1[int(1)].cx_0 - mean_3;
    float _S2217 = p_2.vignette_params_1[int(2)].cx_0 - mean_3;
    losses_3[int(5)] = (_S2215 * _S2215 + _S2216 * _S2216 + _S2217 * _S2217) / 3.0f;
    float mean_4 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S2218 = p_2.vignette_params_1[int(0)].cy_0 - mean_4;
    float _S2219 = p_2.vignette_params_1[int(1)].cy_0 - mean_4;
    float _S2220 = p_2.vignette_params_1[int(2)].cy_0 - mean_4;
    losses_3[int(6)] = (_S2218 * _S2218 + _S2219 * _S2219 + _S2220 * _S2220) / 3.0f;
    float mean_5 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S2221 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_5;
    float _S2222 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_5;
    float _S2223 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_5;
    losses_3[int(7)] = (_S2221 * _S2221 + _S2222 * _S2222 + _S2223 * _S2223) / 3.0f;
    float mean_6 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S2224 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_6;
    float _S2225 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_6;
    float _S2226 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_6;
    losses_3[int(8)] = (_S2224 * _S2224 + _S2225 * _S2225 + _S2226 * _S2226) / 3.0f;
    float mean_7 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S2227 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_7;
    float _S2228 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_7;
    float _S2229 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_7;
    losses_3[int(9)] = (_S2227 * _S2227 + _S2228 * _S2228 + _S2229 * _S2229) / 3.0f;
    float2  bd_2 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_2.color_params_1.b_0);
    float2  rd_2 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_2.color_params_1.r_0);
    float2  gd_2 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_2.color_params_1.g_0);
    float2  nd_2 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_2.color_params_1.n_0);
    losses_3[int(10)] = bd_2.x;
    losses_3[int(11)] = bd_2.y;
    losses_3[int(12)] = rd_2.x;
    losses_3[int(13)] = rd_2.y;
    losses_3[int(14)] = gd_2.x;
    losses_3[int(15)] = gd_2.y;
    losses_3[int(16)] = nd_2.x;
    losses_3[int(17)] = nd_2.y;
    float mean_8 = (p_2.crf_params_1[int(0)].toe_0 + p_2.crf_params_1[int(1)].toe_0 + p_2.crf_params_1[int(2)].toe_0) / 3.0f;
    float _S2230 = p_2.crf_params_1[int(0)].toe_0 - mean_8;
    float _S2231 = p_2.crf_params_1[int(1)].toe_0 - mean_8;
    float _S2232 = p_2.crf_params_1[int(2)].toe_0 - mean_8;
    losses_3[int(18)] = (_S2230 * _S2230 + _S2231 * _S2231 + _S2232 * _S2232) / 3.0f;
    float mean_9 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S2233 = p_2.crf_params_1[int(0)].shoulder_0 - mean_9;
    float _S2234 = p_2.crf_params_1[int(1)].shoulder_0 - mean_9;
    float _S2235 = p_2.crf_params_1[int(2)].shoulder_0 - mean_9;
    losses_3[int(19)] = (_S2233 * _S2233 + _S2234 * _S2234 + _S2235 * _S2235) / 3.0f;
    float mean_10 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S2236 = p_2.crf_params_1[int(0)].gamma_0 - mean_10;
    float _S2237 = p_2.crf_params_1[int(1)].gamma_0 - mean_10;
    float _S2238 = p_2.crf_params_1[int(2)].gamma_0 - mean_10;
    losses_3[int(20)] = (_S2236 * _S2236 + _S2237 * _S2237 + _S2238 * _S2238) / 3.0f;
    float mean_11 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S2239 = p_2.crf_params_1[int(0)].center_0 - mean_11;
    float _S2240 = p_2.crf_params_1[int(1)].center_0 - mean_11;
    float _S2241 = p_2.crf_params_1[int(2)].center_0 - mean_11;
    losses_3[int(21)] = (_S2239 * _S2239 + _S2240 * _S2240 + _S2241 * _S2241) / 3.0f;
    *_S2208 = losses_3;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_5, FixedArray<float, 23>  * _S2242)
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
    float _S2243 = p_3.vignette_params_0[int(0)].cx_0;
    float _S2244 = p_3.vignette_params_0[int(0)].cy_0;
    float _S2245 = p_3.vignette_params_0[int(1)].cx_0;
    float _S2246 = p_3.vignette_params_0[int(1)].cy_0;
    float _S2247 = p_3.vignette_params_0[int(2)].cx_0;
    float _S2248 = p_3.vignette_params_0[int(2)].cy_0;
    losses_4[int(1)] = _S2243 * _S2243 + _S2244 * _S2244 + _S2245 * _S2245 + _S2246 * _S2246 + _S2247 * _S2247 + _S2248 * _S2248;
    losses_4[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_4[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_4[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_12 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S2249 = p_3.vignette_params_0[int(0)].cx_0 - mean_12;
    float _S2250 = p_3.vignette_params_0[int(1)].cx_0 - mean_12;
    float _S2251 = p_3.vignette_params_0[int(2)].cx_0 - mean_12;
    losses_4[int(5)] = (_S2249 * _S2249 + _S2250 * _S2250 + _S2251 * _S2251) / 3.0f;
    float mean_13 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S2252 = p_3.vignette_params_0[int(0)].cy_0 - mean_13;
    float _S2253 = p_3.vignette_params_0[int(1)].cy_0 - mean_13;
    float _S2254 = p_3.vignette_params_0[int(2)].cy_0 - mean_13;
    losses_4[int(6)] = (_S2252 * _S2252 + _S2253 * _S2253 + _S2254 * _S2254) / 3.0f;
    float mean_14 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S2255 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_14;
    float _S2256 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_14;
    float _S2257 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_14;
    losses_4[int(7)] = (_S2255 * _S2255 + _S2256 * _S2256 + _S2257 * _S2257) / 3.0f;
    float mean_15 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S2258 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_15;
    float _S2259 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_15;
    float _S2260 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_15;
    losses_4[int(8)] = (_S2258 * _S2258 + _S2259 * _S2259 + _S2260 * _S2260) / 3.0f;
    float mean_16 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S2261 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_16;
    float _S2262 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_16;
    float _S2263 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_16;
    losses_4[int(9)] = (_S2261 * _S2261 + _S2262 * _S2262 + _S2263 * _S2263) / 3.0f;
    float2  bd_3 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_3.color_params_0.b_0);
    float2  rd_3 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_3.color_params_0.r_0);
    float2  gd_3 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_3.color_params_0.g_0);
    float2  nd_3 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_3.color_params_0.n_0);
    losses_4[int(10)] = bd_3.x;
    losses_4[int(11)] = bd_3.y;
    losses_4[int(12)] = rd_3.x;
    losses_4[int(13)] = rd_3.y;
    losses_4[int(14)] = gd_3.x;
    losses_4[int(15)] = gd_3.y;
    losses_4[int(16)] = nd_3.x;
    losses_4[int(17)] = nd_3.y;
    float mean_17 = (p_3.crf_params_0[int(0)].g0_0 + p_3.crf_params_0[int(1)].g0_0 + p_3.crf_params_0[int(2)].g0_0) / 3.0f;
    float _S2264 = p_3.crf_params_0[int(0)].g0_0 - mean_17;
    float _S2265 = p_3.crf_params_0[int(1)].g0_0 - mean_17;
    float _S2266 = p_3.crf_params_0[int(2)].g0_0 - mean_17;
    losses_4[int(18)] = (_S2264 * _S2264 + _S2265 * _S2265 + _S2266 * _S2266) / 3.0f;
    float mean_18 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S2267 = p_3.crf_params_0[int(0)].g1_0 - mean_18;
    float _S2268 = p_3.crf_params_0[int(1)].g1_0 - mean_18;
    float _S2269 = p_3.crf_params_0[int(2)].g1_0 - mean_18;
    losses_4[int(19)] = (_S2267 * _S2267 + _S2268 * _S2268 + _S2269 * _S2269) / 3.0f;
    float mean_19 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S2270 = p_3.crf_params_0[int(0)].x0_0 - mean_19;
    float _S2271 = p_3.crf_params_0[int(1)].x0_0 - mean_19;
    float _S2272 = p_3.crf_params_0[int(2)].x0_0 - mean_19;
    losses_4[int(20)] = (_S2270 * _S2270 + _S2271 * _S2271 + _S2272 * _S2272) / 3.0f;
    float mean_20 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S2273 = p_3.crf_params_0[int(0)].y0_0 - mean_20;
    float _S2274 = p_3.crf_params_0[int(1)].y0_0 - mean_20;
    float _S2275 = p_3.crf_params_0[int(2)].y0_0 - mean_20;
    losses_4[int(21)] = (_S2273 * _S2273 + _S2274 * _S2274 + _S2275 * _S2275) / 3.0f;
    float mean_21 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S2276 = p_3.crf_params_0[int(0)].gc_0 - mean_21;
    float _S2277 = p_3.crf_params_0[int(1)].gc_0 - mean_21;
    float _S2278 = p_3.crf_params_0[int(2)].gc_0 - mean_21;
    losses_4[int(22)] = (_S2276 * _S2276 + _S2277 * _S2277 + _S2278 * _S2278) / 3.0f;
    *_S2242 = losses_4;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_13)
{
    VignettingChannelParams_0 _S2279 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2280 = {
        _S2279, _S2279, _S2279
    };
    float2  _S2281 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2282 = { _S2281, _S2281, _S2281, _S2281 };
    CRFPPISPChannelParams_0 _S2283 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S2284 = {
        _S2283, _S2283, _S2283
    };
    PPISPParams_0 _S2285;
    (&_S2285)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S2285)->vignette_params_1 = _S2280;
    (&_S2285)->color_params_1 = _S2282;
    (&_S2285)->crf_params_1 = _S2284;
    (&(&_S2285)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S2285)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S2285)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S2285)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S2285)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S2285)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S2285)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S2285)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S2285)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S2285)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S2285)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S2285)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S2285)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S2285)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S2285)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S2285)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S2285)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S2285)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S2285)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S2285)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S2285)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S2285)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S2285)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S2285)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S2285)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S2285)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S2285)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S2285)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S2285)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S2285)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S2285)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S2285)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S2285)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S2285)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S2285)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_22 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S2286 = dpparams_2->primal_0[int(1)] - mean_22;
    float _S2287 = dpparams_2->primal_0[int(6)] - mean_22;
    float _S2288 = dpparams_2->primal_0[int(11)] - mean_22;
    float mean_23 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S2289 = dpparams_2->primal_0[int(2)] - mean_23;
    float _S2290 = dpparams_2->primal_0[int(7)] - mean_23;
    float _S2291 = dpparams_2->primal_0[int(12)] - mean_23;
    float mean_24 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S2292 = dpparams_2->primal_0[int(3)] - mean_24;
    float _S2293 = dpparams_2->primal_0[int(8)] - mean_24;
    float _S2294 = dpparams_2->primal_0[int(13)] - mean_24;
    float mean_25 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S2295 = dpparams_2->primal_0[int(4)] - mean_25;
    float _S2296 = dpparams_2->primal_0[int(9)] - mean_25;
    float _S2297 = dpparams_2->primal_0[int(14)] - mean_25;
    float mean_26 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S2298 = dpparams_2->primal_0[int(5)] - mean_26;
    float _S2299 = dpparams_2->primal_0[int(10)] - mean_26;
    float _S2300 = dpparams_2->primal_0[int(15)] - mean_26;
    float mean_27 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_28 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_29 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_30 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S2301 = 0.3333333432674408f * (*_s_dOut_13)[int(21)];
    float _S2302 = (dpparams_2->primal_0[int(35)] - mean_30) * _S2301;
    float _S2303 = _S2302 + _S2302;
    float _S2304 = (dpparams_2->primal_0[int(31)] - mean_30) * _S2301;
    float _S2305 = _S2304 + _S2304;
    float _S2306 = (dpparams_2->primal_0[int(27)] - mean_30) * _S2301;
    float _S2307 = _S2306 + _S2306;
    float _S2308 = 0.3333333432674408f * (- _S2303 + - _S2305 + - _S2307);
    float _S2309 = 0.3333333432674408f * (*_s_dOut_13)[int(20)];
    float _S2310 = (dpparams_2->primal_0[int(34)] - mean_29) * _S2309;
    float _S2311 = _S2310 + _S2310;
    float _S2312 = (dpparams_2->primal_0[int(30)] - mean_29) * _S2309;
    float _S2313 = _S2312 + _S2312;
    float _S2314 = (dpparams_2->primal_0[int(26)] - mean_29) * _S2309;
    float _S2315 = _S2314 + _S2314;
    float _S2316 = 0.3333333432674408f * (- _S2311 + - _S2313 + - _S2315);
    float _S2317 = 0.3333333432674408f * (*_s_dOut_13)[int(19)];
    float _S2318 = (dpparams_2->primal_0[int(33)] - mean_28) * _S2317;
    float _S2319 = _S2318 + _S2318;
    float _S2320 = (dpparams_2->primal_0[int(29)] - mean_28) * _S2317;
    float _S2321 = _S2320 + _S2320;
    float _S2322 = (dpparams_2->primal_0[int(25)] - mean_28) * _S2317;
    float _S2323 = _S2322 + _S2322;
    float _S2324 = 0.3333333432674408f * (- _S2319 + - _S2321 + - _S2323);
    float _S2325 = 0.3333333432674408f * (*_s_dOut_13)[int(18)];
    float _S2326 = (dpparams_2->primal_0[int(32)] - mean_27) * _S2325;
    float _S2327 = _S2326 + _S2326;
    float _S2328 = (dpparams_2->primal_0[int(28)] - mean_27) * _S2325;
    float _S2329 = _S2328 + _S2328;
    float _S2330 = (dpparams_2->primal_0[int(24)] - mean_27) * _S2325;
    float _S2331 = _S2330 + _S2330;
    float _S2332 = 0.3333333432674408f * (- _S2327 + - _S2329 + - _S2331);
    float2  _S2333 = make_float2 ((*_s_dOut_13)[int(16)], (*_s_dOut_13)[int(17)]);
    Matrix<float, 2, 2>  _S2334 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2335;
    (&_S2335)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2335)->differential_0 = _S2334;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2336;
    (&_S2336)->primal_0 = _S2285.color_params_1.n_0;
    (&_S2336)->differential_0 = _S2281;
    s_bwd_prop_mul_3(&_S2335, &_S2336, _S2333);
    float2  _S2337 = make_float2 ((*_s_dOut_13)[int(14)], (*_s_dOut_13)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2338;
    (&_S2338)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2338)->differential_0 = _S2334;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2339;
    (&_S2339)->primal_0 = _S2285.color_params_1.g_0;
    (&_S2339)->differential_0 = _S2281;
    s_bwd_prop_mul_3(&_S2338, &_S2339, _S2337);
    float2  _S2340 = make_float2 ((*_s_dOut_13)[int(12)], (*_s_dOut_13)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2341;
    (&_S2341)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2341)->differential_0 = _S2334;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2342;
    (&_S2342)->primal_0 = _S2285.color_params_1.r_0;
    (&_S2342)->differential_0 = _S2281;
    s_bwd_prop_mul_3(&_S2341, &_S2342, _S2340);
    float2  _S2343 = make_float2 ((*_s_dOut_13)[int(10)], (*_s_dOut_13)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2344;
    (&_S2344)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2344)->differential_0 = _S2334;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2345;
    (&_S2345)->primal_0 = _S2285.color_params_1.b_0;
    (&_S2345)->differential_0 = _S2281;
    s_bwd_prop_mul_3(&_S2344, &_S2345, _S2343);
    ColorPPISPParams_0 _S2346 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2346)->n_0 = _S2336.differential_0;
    (&_S2346)->g_0 = _S2339.differential_0;
    (&_S2346)->r_0 = _S2342.differential_0;
    (&_S2346)->b_0 = _S2345.differential_0;
    float _S2347 = 0.3333333432674408f * (*_s_dOut_13)[int(9)];
    float _S2348 = _S2300 * _S2347;
    float _S2349 = _S2348 + _S2348;
    float _S2350 = _S2299 * _S2347;
    float _S2351 = _S2350 + _S2350;
    float _S2352 = _S2298 * _S2347;
    float _S2353 = _S2352 + _S2352;
    float _S2354 = 0.3333333432674408f * (- _S2349 + - _S2351 + - _S2353);
    float _S2355 = 0.3333333432674408f * (*_s_dOut_13)[int(8)];
    float _S2356 = _S2297 * _S2355;
    float _S2357 = _S2356 + _S2356;
    float _S2358 = _S2296 * _S2355;
    float _S2359 = _S2358 + _S2358;
    float _S2360 = _S2295 * _S2355;
    float _S2361 = _S2360 + _S2360;
    float _S2362 = 0.3333333432674408f * (- _S2357 + - _S2359 + - _S2361);
    float _S2363 = 0.3333333432674408f * (*_s_dOut_13)[int(7)];
    float _S2364 = _S2294 * _S2363;
    float _S2365 = _S2364 + _S2364;
    float _S2366 = _S2293 * _S2363;
    float _S2367 = _S2366 + _S2366;
    float _S2368 = _S2292 * _S2363;
    float _S2369 = _S2368 + _S2368;
    float _S2370 = 0.3333333432674408f * (- _S2365 + - _S2367 + - _S2369);
    float _S2371 = 0.3333333432674408f * (*_s_dOut_13)[int(6)];
    float _S2372 = _S2291 * _S2371;
    float _S2373 = _S2372 + _S2372;
    float _S2374 = _S2290 * _S2371;
    float _S2375 = _S2374 + _S2374;
    float _S2376 = _S2289 * _S2371;
    float _S2377 = _S2376 + _S2376;
    float _S2378 = 0.3333333432674408f * (- _S2373 + - _S2375 + - _S2377);
    float _S2379 = 0.3333333432674408f * (*_s_dOut_13)[int(5)];
    float _S2380 = _S2288 * _S2379;
    float _S2381 = _S2380 + _S2380;
    float _S2382 = _S2287 * _S2379;
    float _S2383 = _S2382 + _S2382;
    float _S2384 = _S2286 * _S2379;
    float _S2385 = _S2384 + _S2384;
    float _S2386 = 0.3333333432674408f * (- _S2381 + - _S2383 + - _S2385);
    DiffPair_float_0 _S2387;
    (&_S2387)->primal_0 = 0.0f;
    (&_S2387)->differential_0 = 0.0f;
    DiffPair_float_0 _S2388;
    (&_S2388)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S2388)->differential_0 = 0.0f;
    _d_max_0(&_S2387, &_S2388, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2389;
    (&_S2389)->primal_0 = 0.0f;
    (&_S2389)->differential_0 = 0.0f;
    DiffPair_float_0 _S2390;
    (&_S2390)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S2390)->differential_0 = 0.0f;
    _d_max_0(&_S2389, &_S2390, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2391;
    (&_S2391)->primal_0 = 0.0f;
    (&_S2391)->differential_0 = 0.0f;
    DiffPair_float_0 _S2392;
    (&_S2392)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S2392)->differential_0 = 0.0f;
    _d_max_0(&_S2391, &_S2392, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2393;
    (&_S2393)->primal_0 = 0.0f;
    (&_S2393)->differential_0 = 0.0f;
    DiffPair_float_0 _S2394;
    (&_S2394)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S2394)->differential_0 = 0.0f;
    _d_max_0(&_S2393, &_S2394, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2395;
    (&_S2395)->primal_0 = 0.0f;
    (&_S2395)->differential_0 = 0.0f;
    DiffPair_float_0 _S2396;
    (&_S2396)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S2396)->differential_0 = 0.0f;
    _d_max_0(&_S2395, &_S2396, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2397;
    (&_S2397)->primal_0 = 0.0f;
    (&_S2397)->differential_0 = 0.0f;
    DiffPair_float_0 _S2398;
    (&_S2398)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S2398)->differential_0 = 0.0f;
    _d_max_0(&_S2397, &_S2398, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2399;
    (&_S2399)->primal_0 = 0.0f;
    (&_S2399)->differential_0 = 0.0f;
    DiffPair_float_0 _S2400;
    (&_S2400)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S2400)->differential_0 = 0.0f;
    _d_max_0(&_S2399, &_S2400, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2401;
    (&_S2401)->primal_0 = 0.0f;
    (&_S2401)->differential_0 = 0.0f;
    DiffPair_float_0 _S2402;
    (&_S2402)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S2402)->differential_0 = 0.0f;
    _d_max_0(&_S2401, &_S2402, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2403;
    (&_S2403)->primal_0 = 0.0f;
    (&_S2403)->differential_0 = 0.0f;
    DiffPair_float_0 _S2404;
    (&_S2404)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S2404)->differential_0 = 0.0f;
    _d_max_0(&_S2403, &_S2404, (*_s_dOut_13)[int(2)]);
    float _S2405 = dpparams_2->primal_0[int(12)] * (*_s_dOut_13)[int(1)];
    float _S2406 = dpparams_2->primal_0[int(11)] * (*_s_dOut_13)[int(1)];
    float _S2407 = dpparams_2->primal_0[int(7)] * (*_s_dOut_13)[int(1)];
    float _S2408 = dpparams_2->primal_0[int(6)] * (*_s_dOut_13)[int(1)];
    float _S2409 = dpparams_2->primal_0[int(2)] * (*_s_dOut_13)[int(1)];
    float _S2410 = dpparams_2->primal_0[int(1)] * (*_s_dOut_13)[int(1)];
    PPISPParams_0 _S2411 = PPISPParams_x24_syn_dzero_0();
    (&_S2411)->color_params_1 = _S2346;
    (&_S2411)->exposure_1 = (*_s_dOut_13)[int(0)];
    _S2285 = _S2411;
    (&(&_S2285)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2412 = _S2303 + _S2308 + _S2411.crf_params_1[int(2)].center_0;
    (&(&_S2285)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2413 = _S2311 + _S2316 + _S2411.crf_params_1[int(2)].gamma_0;
    (&(&_S2285)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2414 = _S2319 + _S2324 + _S2411.crf_params_1[int(2)].shoulder_0;
    (&(&_S2285)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2415 = _S2327 + _S2332 + _S2411.crf_params_1[int(2)].toe_0;
    (&(&_S2285)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2416 = _S2305 + _S2308 + _S2411.crf_params_1[int(1)].center_0;
    (&(&_S2285)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2417 = _S2313 + _S2316 + _S2411.crf_params_1[int(1)].gamma_0;
    (&(&_S2285)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2418 = _S2321 + _S2324 + _S2411.crf_params_1[int(1)].shoulder_0;
    (&(&_S2285)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2419 = _S2329 + _S2332 + _S2411.crf_params_1[int(1)].toe_0;
    (&(&_S2285)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2420 = _S2307 + _S2308 + _S2411.crf_params_1[int(0)].center_0;
    (&(&_S2285)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2421 = _S2315 + _S2316 + _S2411.crf_params_1[int(0)].gamma_0;
    (&(&_S2285)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2422 = _S2323 + _S2324 + _S2411.crf_params_1[int(0)].shoulder_0;
    (&(&_S2285)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2423 = _S2331 + _S2332 + _S2411.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S2285)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S2285)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S2285)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2424 = _S2349 + _S2354 + _S2388.differential_0 + _S2411.vignette_params_1[int(2)].alpha2_0;
    (&(&_S2285)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2425 = _S2357 + _S2362 + _S2394.differential_0 + _S2411.vignette_params_1[int(2)].alpha1_0;
    (&(&_S2285)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2426 = _S2365 + _S2370 + _S2400.differential_0 + _S2411.vignette_params_1[int(2)].alpha0_0;
    (&(&_S2285)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2427 = _S2373 + _S2378 + _S2405 + _S2405 + _S2411.vignette_params_1[int(2)].cy_0;
    (&(&_S2285)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2428 = _S2381 + _S2386 + _S2406 + _S2406 + _S2411.vignette_params_1[int(2)].cx_0;
    (&(&_S2285)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2429 = _S2351 + _S2354 + _S2390.differential_0 + _S2411.vignette_params_1[int(1)].alpha2_0;
    (&(&_S2285)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2430 = _S2359 + _S2362 + _S2396.differential_0 + _S2411.vignette_params_1[int(1)].alpha1_0;
    (&(&_S2285)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2431 = _S2367 + _S2370 + _S2402.differential_0 + _S2411.vignette_params_1[int(1)].alpha0_0;
    (&(&_S2285)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2432 = _S2375 + _S2378 + _S2407 + _S2407 + _S2411.vignette_params_1[int(1)].cy_0;
    (&(&_S2285)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2433 = _S2383 + _S2386 + _S2408 + _S2408 + _S2411.vignette_params_1[int(1)].cx_0;
    (&(&_S2285)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2434 = _S2353 + _S2354 + _S2392.differential_0 + _S2411.vignette_params_1[int(0)].alpha2_0;
    (&(&_S2285)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2435 = _S2361 + _S2362 + _S2398.differential_0 + _S2411.vignette_params_1[int(0)].alpha1_0;
    (&(&_S2285)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2436 = _S2369 + _S2370 + _S2404.differential_0 + _S2411.vignette_params_1[int(0)].alpha0_0;
    (&(&_S2285)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2437 = _S2377 + _S2378 + _S2409 + _S2409 + _S2411.vignette_params_1[int(0)].cy_0;
    (&(&_S2285)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2438 = _S2385 + _S2386 + _S2410 + _S2410 + _S2411.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2439;
    _S2439[int(0)] = 0.0f;
    _S2439[int(1)] = 0.0f;
    _S2439[int(2)] = 0.0f;
    _S2439[int(3)] = 0.0f;
    _S2439[int(4)] = 0.0f;
    _S2439[int(5)] = 0.0f;
    _S2439[int(6)] = 0.0f;
    _S2439[int(7)] = 0.0f;
    _S2439[int(8)] = 0.0f;
    _S2439[int(9)] = 0.0f;
    _S2439[int(10)] = 0.0f;
    _S2439[int(11)] = 0.0f;
    _S2439[int(12)] = 0.0f;
    _S2439[int(13)] = 0.0f;
    _S2439[int(14)] = 0.0f;
    _S2439[int(15)] = 0.0f;
    _S2439[int(16)] = 0.0f;
    _S2439[int(17)] = 0.0f;
    _S2439[int(18)] = 0.0f;
    _S2439[int(19)] = 0.0f;
    _S2439[int(20)] = 0.0f;
    _S2439[int(21)] = 0.0f;
    _S2439[int(22)] = 0.0f;
    _S2439[int(23)] = 0.0f;
    _S2439[int(24)] = 0.0f;
    _S2439[int(25)] = 0.0f;
    _S2439[int(26)] = 0.0f;
    _S2439[int(27)] = 0.0f;
    _S2439[int(28)] = 0.0f;
    _S2439[int(29)] = 0.0f;
    _S2439[int(30)] = 0.0f;
    _S2439[int(31)] = 0.0f;
    _S2439[int(32)] = 0.0f;
    _S2439[int(33)] = 0.0f;
    _S2439[int(34)] = 0.0f;
    _S2439[int(35)] = 0.0f;
    _S2439[int(8)] = _S2431;
    _S2439[int(16)] = _S2411.color_params_1.b_0.x;
    _S2439[int(15)] = _S2424;
    _S2439[int(14)] = _S2425;
    _S2439[int(13)] = _S2426;
    _S2439[int(12)] = _S2427;
    _S2439[int(11)] = _S2428;
    _S2439[int(10)] = _S2429;
    _S2439[int(9)] = _S2430;
    _S2439[int(17)] = _S2411.color_params_1.b_0.y;
    _S2439[int(7)] = _S2432;
    _S2439[int(6)] = _S2433;
    _S2439[int(5)] = _S2434;
    _S2439[int(4)] = _S2435;
    _S2439[int(3)] = _S2436;
    _S2439[int(2)] = _S2437;
    _S2439[int(1)] = _S2438;
    _S2439[int(0)] = _S2285.exposure_1;
    _S2439[int(26)] = _S2421;
    _S2439[int(34)] = _S2413;
    _S2439[int(33)] = _S2414;
    _S2439[int(32)] = _S2415;
    _S2439[int(31)] = _S2416;
    _S2439[int(30)] = _S2417;
    _S2439[int(29)] = _S2418;
    _S2439[int(28)] = _S2419;
    _S2439[int(27)] = _S2420;
    _S2439[int(35)] = _S2412;
    _S2439[int(25)] = _S2422;
    _S2439[int(24)] = _S2423;
    _S2439[int(23)] = _S2411.color_params_1.n_0.y;
    _S2439[int(22)] = _S2411.color_params_1.n_0.x;
    _S2439[int(21)] = _S2411.color_params_1.g_0.y;
    _S2439[int(20)] = _S2411.color_params_1.g_0.x;
    _S2439[int(19)] = _S2411.color_params_1.r_0.y;
    _S2439[int(18)] = _S2411.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S2439;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2440, FixedArray<float, 22>  * _S2441)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S2440, _S2441);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_6, FixedArray<float, 22>  grad_out_2, FixedArray<float, 36>  * _S2442)
{
    FixedArray<float, 36>  _S2443 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_6;
    (&dp_params_2)->differential_0 = _S2443;
    FixedArray<float, 22>  _S2444 = grad_out_2;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, &_S2444);
    *_S2442 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_14)
{
    VignettingChannelParams_0 _S2445 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2446 = {
        _S2445, _S2445, _S2445
    };
    float2  _S2447 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2448 = { _S2447, _S2447, _S2447, _S2447 };
    RQSCRFPPISPChannelParams_0 _S2449 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2450 = {
        _S2449, _S2449, _S2449
    };
    PPISPParamsRQS_0 _S2451;
    (&_S2451)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S2451)->vignette_params_0 = _S2446;
    (&_S2451)->color_params_0 = _S2448;
    (&_S2451)->crf_params_0 = _S2450;
    (&(&_S2451)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S2451)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S2451)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S2451)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S2451)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S2451)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S2451)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S2451)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S2451)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S2451)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S2451)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S2451)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S2451)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S2451)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S2451)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S2451)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S2451)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S2451)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S2451)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S2451)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S2451)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S2451)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S2451)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S2451)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S2451)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S2451)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S2451)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S2451)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S2451)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S2451)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S2451)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S2451)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S2451)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S2451)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S2451)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S2451)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S2451)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S2451)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_31 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S2452 = dpparams_3->primal_0[int(1)] - mean_31;
    float _S2453 = dpparams_3->primal_0[int(6)] - mean_31;
    float _S2454 = dpparams_3->primal_0[int(11)] - mean_31;
    float mean_32 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S2455 = dpparams_3->primal_0[int(2)] - mean_32;
    float _S2456 = dpparams_3->primal_0[int(7)] - mean_32;
    float _S2457 = dpparams_3->primal_0[int(12)] - mean_32;
    float mean_33 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S2458 = dpparams_3->primal_0[int(3)] - mean_33;
    float _S2459 = dpparams_3->primal_0[int(8)] - mean_33;
    float _S2460 = dpparams_3->primal_0[int(13)] - mean_33;
    float mean_34 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S2461 = dpparams_3->primal_0[int(4)] - mean_34;
    float _S2462 = dpparams_3->primal_0[int(9)] - mean_34;
    float _S2463 = dpparams_3->primal_0[int(14)] - mean_34;
    float mean_35 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S2464 = dpparams_3->primal_0[int(5)] - mean_35;
    float _S2465 = dpparams_3->primal_0[int(10)] - mean_35;
    float _S2466 = dpparams_3->primal_0[int(15)] - mean_35;
    float mean_36 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_37 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_38 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_39 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_40 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S2467 = 0.3333333432674408f * (*_s_dOut_14)[int(22)];
    float _S2468 = (dpparams_3->primal_0[int(38)] - mean_40) * _S2467;
    float _S2469 = _S2468 + _S2468;
    float _S2470 = (dpparams_3->primal_0[int(33)] - mean_40) * _S2467;
    float _S2471 = _S2470 + _S2470;
    float _S2472 = (dpparams_3->primal_0[int(28)] - mean_40) * _S2467;
    float _S2473 = _S2472 + _S2472;
    float _S2474 = 0.3333333432674408f * (- _S2469 + - _S2471 + - _S2473);
    float _S2475 = 0.3333333432674408f * (*_s_dOut_14)[int(21)];
    float _S2476 = (dpparams_3->primal_0[int(37)] - mean_39) * _S2475;
    float _S2477 = _S2476 + _S2476;
    float _S2478 = (dpparams_3->primal_0[int(32)] - mean_39) * _S2475;
    float _S2479 = _S2478 + _S2478;
    float _S2480 = (dpparams_3->primal_0[int(27)] - mean_39) * _S2475;
    float _S2481 = _S2480 + _S2480;
    float _S2482 = 0.3333333432674408f * (- _S2477 + - _S2479 + - _S2481);
    float _S2483 = 0.3333333432674408f * (*_s_dOut_14)[int(20)];
    float _S2484 = (dpparams_3->primal_0[int(36)] - mean_38) * _S2483;
    float _S2485 = _S2484 + _S2484;
    float _S2486 = (dpparams_3->primal_0[int(31)] - mean_38) * _S2483;
    float _S2487 = _S2486 + _S2486;
    float _S2488 = (dpparams_3->primal_0[int(26)] - mean_38) * _S2483;
    float _S2489 = _S2488 + _S2488;
    float _S2490 = 0.3333333432674408f * (- _S2485 + - _S2487 + - _S2489);
    float _S2491 = 0.3333333432674408f * (*_s_dOut_14)[int(19)];
    float _S2492 = (dpparams_3->primal_0[int(35)] - mean_37) * _S2491;
    float _S2493 = _S2492 + _S2492;
    float _S2494 = (dpparams_3->primal_0[int(30)] - mean_37) * _S2491;
    float _S2495 = _S2494 + _S2494;
    float _S2496 = (dpparams_3->primal_0[int(25)] - mean_37) * _S2491;
    float _S2497 = _S2496 + _S2496;
    float _S2498 = 0.3333333432674408f * (- _S2493 + - _S2495 + - _S2497);
    float _S2499 = 0.3333333432674408f * (*_s_dOut_14)[int(18)];
    float _S2500 = (dpparams_3->primal_0[int(34)] - mean_36) * _S2499;
    float _S2501 = _S2500 + _S2500;
    float _S2502 = (dpparams_3->primal_0[int(29)] - mean_36) * _S2499;
    float _S2503 = _S2502 + _S2502;
    float _S2504 = (dpparams_3->primal_0[int(24)] - mean_36) * _S2499;
    float _S2505 = _S2504 + _S2504;
    float _S2506 = 0.3333333432674408f * (- _S2501 + - _S2503 + - _S2505);
    float2  _S2507 = make_float2 ((*_s_dOut_14)[int(16)], (*_s_dOut_14)[int(17)]);
    Matrix<float, 2, 2>  _S2508 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2509;
    (&_S2509)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2509)->differential_0 = _S2508;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2510;
    (&_S2510)->primal_0 = _S2451.color_params_0.n_0;
    (&_S2510)->differential_0 = _S2447;
    s_bwd_prop_mul_3(&_S2509, &_S2510, _S2507);
    float2  _S2511 = make_float2 ((*_s_dOut_14)[int(14)], (*_s_dOut_14)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2512;
    (&_S2512)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2512)->differential_0 = _S2508;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2513;
    (&_S2513)->primal_0 = _S2451.color_params_0.g_0;
    (&_S2513)->differential_0 = _S2447;
    s_bwd_prop_mul_3(&_S2512, &_S2513, _S2511);
    float2  _S2514 = make_float2 ((*_s_dOut_14)[int(12)], (*_s_dOut_14)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2515;
    (&_S2515)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2515)->differential_0 = _S2508;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2516;
    (&_S2516)->primal_0 = _S2451.color_params_0.r_0;
    (&_S2516)->differential_0 = _S2447;
    s_bwd_prop_mul_3(&_S2515, &_S2516, _S2514);
    float2  _S2517 = make_float2 ((*_s_dOut_14)[int(10)], (*_s_dOut_14)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2518;
    (&_S2518)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2518)->differential_0 = _S2508;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2519;
    (&_S2519)->primal_0 = _S2451.color_params_0.b_0;
    (&_S2519)->differential_0 = _S2447;
    s_bwd_prop_mul_3(&_S2518, &_S2519, _S2517);
    ColorPPISPParams_0 _S2520 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2520)->n_0 = _S2510.differential_0;
    (&_S2520)->g_0 = _S2513.differential_0;
    (&_S2520)->r_0 = _S2516.differential_0;
    (&_S2520)->b_0 = _S2519.differential_0;
    float _S2521 = 0.3333333432674408f * (*_s_dOut_14)[int(9)];
    float _S2522 = _S2466 * _S2521;
    float _S2523 = _S2522 + _S2522;
    float _S2524 = _S2465 * _S2521;
    float _S2525 = _S2524 + _S2524;
    float _S2526 = _S2464 * _S2521;
    float _S2527 = _S2526 + _S2526;
    float _S2528 = 0.3333333432674408f * (- _S2523 + - _S2525 + - _S2527);
    float _S2529 = 0.3333333432674408f * (*_s_dOut_14)[int(8)];
    float _S2530 = _S2463 * _S2529;
    float _S2531 = _S2530 + _S2530;
    float _S2532 = _S2462 * _S2529;
    float _S2533 = _S2532 + _S2532;
    float _S2534 = _S2461 * _S2529;
    float _S2535 = _S2534 + _S2534;
    float _S2536 = 0.3333333432674408f * (- _S2531 + - _S2533 + - _S2535);
    float _S2537 = 0.3333333432674408f * (*_s_dOut_14)[int(7)];
    float _S2538 = _S2460 * _S2537;
    float _S2539 = _S2538 + _S2538;
    float _S2540 = _S2459 * _S2537;
    float _S2541 = _S2540 + _S2540;
    float _S2542 = _S2458 * _S2537;
    float _S2543 = _S2542 + _S2542;
    float _S2544 = 0.3333333432674408f * (- _S2539 + - _S2541 + - _S2543);
    float _S2545 = 0.3333333432674408f * (*_s_dOut_14)[int(6)];
    float _S2546 = _S2457 * _S2545;
    float _S2547 = _S2546 + _S2546;
    float _S2548 = _S2456 * _S2545;
    float _S2549 = _S2548 + _S2548;
    float _S2550 = _S2455 * _S2545;
    float _S2551 = _S2550 + _S2550;
    float _S2552 = 0.3333333432674408f * (- _S2547 + - _S2549 + - _S2551);
    float _S2553 = 0.3333333432674408f * (*_s_dOut_14)[int(5)];
    float _S2554 = _S2454 * _S2553;
    float _S2555 = _S2554 + _S2554;
    float _S2556 = _S2453 * _S2553;
    float _S2557 = _S2556 + _S2556;
    float _S2558 = _S2452 * _S2553;
    float _S2559 = _S2558 + _S2558;
    float _S2560 = 0.3333333432674408f * (- _S2555 + - _S2557 + - _S2559);
    DiffPair_float_0 _S2561;
    (&_S2561)->primal_0 = 0.0f;
    (&_S2561)->differential_0 = 0.0f;
    DiffPair_float_0 _S2562;
    (&_S2562)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S2562)->differential_0 = 0.0f;
    _d_max_0(&_S2561, &_S2562, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2563;
    (&_S2563)->primal_0 = 0.0f;
    (&_S2563)->differential_0 = 0.0f;
    DiffPair_float_0 _S2564;
    (&_S2564)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S2564)->differential_0 = 0.0f;
    _d_max_0(&_S2563, &_S2564, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2565;
    (&_S2565)->primal_0 = 0.0f;
    (&_S2565)->differential_0 = 0.0f;
    DiffPair_float_0 _S2566;
    (&_S2566)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S2566)->differential_0 = 0.0f;
    _d_max_0(&_S2565, &_S2566, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2567;
    (&_S2567)->primal_0 = 0.0f;
    (&_S2567)->differential_0 = 0.0f;
    DiffPair_float_0 _S2568;
    (&_S2568)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S2568)->differential_0 = 0.0f;
    _d_max_0(&_S2567, &_S2568, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2569;
    (&_S2569)->primal_0 = 0.0f;
    (&_S2569)->differential_0 = 0.0f;
    DiffPair_float_0 _S2570;
    (&_S2570)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S2570)->differential_0 = 0.0f;
    _d_max_0(&_S2569, &_S2570, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2571;
    (&_S2571)->primal_0 = 0.0f;
    (&_S2571)->differential_0 = 0.0f;
    DiffPair_float_0 _S2572;
    (&_S2572)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S2572)->differential_0 = 0.0f;
    _d_max_0(&_S2571, &_S2572, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2573;
    (&_S2573)->primal_0 = 0.0f;
    (&_S2573)->differential_0 = 0.0f;
    DiffPair_float_0 _S2574;
    (&_S2574)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S2574)->differential_0 = 0.0f;
    _d_max_0(&_S2573, &_S2574, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2575;
    (&_S2575)->primal_0 = 0.0f;
    (&_S2575)->differential_0 = 0.0f;
    DiffPair_float_0 _S2576;
    (&_S2576)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S2576)->differential_0 = 0.0f;
    _d_max_0(&_S2575, &_S2576, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2577;
    (&_S2577)->primal_0 = 0.0f;
    (&_S2577)->differential_0 = 0.0f;
    DiffPair_float_0 _S2578;
    (&_S2578)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S2578)->differential_0 = 0.0f;
    _d_max_0(&_S2577, &_S2578, (*_s_dOut_14)[int(2)]);
    float _S2579 = dpparams_3->primal_0[int(12)] * (*_s_dOut_14)[int(1)];
    float _S2580 = dpparams_3->primal_0[int(11)] * (*_s_dOut_14)[int(1)];
    float _S2581 = dpparams_3->primal_0[int(7)] * (*_s_dOut_14)[int(1)];
    float _S2582 = dpparams_3->primal_0[int(6)] * (*_s_dOut_14)[int(1)];
    float _S2583 = dpparams_3->primal_0[int(2)] * (*_s_dOut_14)[int(1)];
    float _S2584 = dpparams_3->primal_0[int(1)] * (*_s_dOut_14)[int(1)];
    PPISPParamsRQS_0 _S2585 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2585)->color_params_0 = _S2520;
    (&_S2585)->exposure_0 = (*_s_dOut_14)[int(0)];
    _S2451 = _S2585;
    (&(&_S2451)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2586 = _S2469 + _S2474 + _S2585.crf_params_0[int(2)].gc_0;
    (&(&_S2451)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2587 = _S2477 + _S2482 + _S2585.crf_params_0[int(2)].y0_0;
    (&(&_S2451)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2588 = _S2485 + _S2490 + _S2585.crf_params_0[int(2)].x0_0;
    (&(&_S2451)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2589 = _S2493 + _S2498 + _S2585.crf_params_0[int(2)].g1_0;
    (&(&_S2451)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2590 = _S2501 + _S2506 + _S2585.crf_params_0[int(2)].g0_0;
    (&(&_S2451)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2591 = _S2471 + _S2474 + _S2585.crf_params_0[int(1)].gc_0;
    (&(&_S2451)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2592 = _S2479 + _S2482 + _S2585.crf_params_0[int(1)].y0_0;
    (&(&_S2451)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2593 = _S2487 + _S2490 + _S2585.crf_params_0[int(1)].x0_0;
    (&(&_S2451)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2594 = _S2495 + _S2498 + _S2585.crf_params_0[int(1)].g1_0;
    (&(&_S2451)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2595 = _S2503 + _S2506 + _S2585.crf_params_0[int(1)].g0_0;
    (&(&_S2451)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2596 = _S2473 + _S2474 + _S2585.crf_params_0[int(0)].gc_0;
    (&(&_S2451)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2597 = _S2481 + _S2482 + _S2585.crf_params_0[int(0)].y0_0;
    (&(&_S2451)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2598 = _S2489 + _S2490 + _S2585.crf_params_0[int(0)].x0_0;
    (&(&_S2451)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2599 = _S2497 + _S2498 + _S2585.crf_params_0[int(0)].g1_0;
    (&(&_S2451)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2600 = _S2505 + _S2506 + _S2585.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S2451)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2451)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2451)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2601 = _S2523 + _S2528 + _S2562.differential_0 + _S2585.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2451)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2602 = _S2531 + _S2536 + _S2568.differential_0 + _S2585.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2451)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2603 = _S2539 + _S2544 + _S2574.differential_0 + _S2585.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2451)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2604 = _S2547 + _S2552 + _S2579 + _S2579 + _S2585.vignette_params_0[int(2)].cy_0;
    (&(&_S2451)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2605 = _S2555 + _S2560 + _S2580 + _S2580 + _S2585.vignette_params_0[int(2)].cx_0;
    (&(&_S2451)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2606 = _S2525 + _S2528 + _S2564.differential_0 + _S2585.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2451)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2607 = _S2533 + _S2536 + _S2570.differential_0 + _S2585.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2451)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2608 = _S2541 + _S2544 + _S2576.differential_0 + _S2585.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2451)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2609 = _S2549 + _S2552 + _S2581 + _S2581 + _S2585.vignette_params_0[int(1)].cy_0;
    (&(&_S2451)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2610 = _S2557 + _S2560 + _S2582 + _S2582 + _S2585.vignette_params_0[int(1)].cx_0;
    (&(&_S2451)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2611 = _S2527 + _S2528 + _S2566.differential_0 + _S2585.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2451)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2612 = _S2535 + _S2536 + _S2572.differential_0 + _S2585.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2451)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2613 = _S2543 + _S2544 + _S2578.differential_0 + _S2585.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2451)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2614 = _S2551 + _S2552 + _S2583 + _S2583 + _S2585.vignette_params_0[int(0)].cy_0;
    (&(&_S2451)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2615 = _S2559 + _S2560 + _S2584 + _S2584 + _S2585.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2616;
    _S2616[int(0)] = 0.0f;
    _S2616[int(1)] = 0.0f;
    _S2616[int(2)] = 0.0f;
    _S2616[int(3)] = 0.0f;
    _S2616[int(4)] = 0.0f;
    _S2616[int(5)] = 0.0f;
    _S2616[int(6)] = 0.0f;
    _S2616[int(7)] = 0.0f;
    _S2616[int(8)] = 0.0f;
    _S2616[int(9)] = 0.0f;
    _S2616[int(10)] = 0.0f;
    _S2616[int(11)] = 0.0f;
    _S2616[int(12)] = 0.0f;
    _S2616[int(13)] = 0.0f;
    _S2616[int(14)] = 0.0f;
    _S2616[int(15)] = 0.0f;
    _S2616[int(16)] = 0.0f;
    _S2616[int(17)] = 0.0f;
    _S2616[int(18)] = 0.0f;
    _S2616[int(19)] = 0.0f;
    _S2616[int(20)] = 0.0f;
    _S2616[int(21)] = 0.0f;
    _S2616[int(22)] = 0.0f;
    _S2616[int(23)] = 0.0f;
    _S2616[int(24)] = 0.0f;
    _S2616[int(25)] = 0.0f;
    _S2616[int(26)] = 0.0f;
    _S2616[int(27)] = 0.0f;
    _S2616[int(28)] = 0.0f;
    _S2616[int(29)] = 0.0f;
    _S2616[int(30)] = 0.0f;
    _S2616[int(31)] = 0.0f;
    _S2616[int(32)] = 0.0f;
    _S2616[int(33)] = 0.0f;
    _S2616[int(34)] = 0.0f;
    _S2616[int(35)] = 0.0f;
    _S2616[int(36)] = 0.0f;
    _S2616[int(37)] = 0.0f;
    _S2616[int(38)] = 0.0f;
    _S2616[int(9)] = _S2607;
    _S2616[int(18)] = _S2585.color_params_0.r_0.x;
    _S2616[int(17)] = _S2585.color_params_0.b_0.y;
    _S2616[int(16)] = _S2585.color_params_0.b_0.x;
    _S2616[int(15)] = _S2601;
    _S2616[int(14)] = _S2602;
    _S2616[int(13)] = _S2603;
    _S2616[int(12)] = _S2604;
    _S2616[int(11)] = _S2605;
    _S2616[int(10)] = _S2606;
    _S2616[int(19)] = _S2585.color_params_0.r_0.y;
    _S2616[int(8)] = _S2608;
    _S2616[int(7)] = _S2609;
    _S2616[int(6)] = _S2610;
    _S2616[int(5)] = _S2611;
    _S2616[int(4)] = _S2612;
    _S2616[int(3)] = _S2613;
    _S2616[int(2)] = _S2614;
    _S2616[int(1)] = _S2615;
    _S2616[int(0)] = _S2451.exposure_0;
    _S2616[int(28)] = _S2596;
    _S2616[int(37)] = _S2587;
    _S2616[int(36)] = _S2588;
    _S2616[int(35)] = _S2589;
    _S2616[int(34)] = _S2590;
    _S2616[int(33)] = _S2591;
    _S2616[int(32)] = _S2592;
    _S2616[int(31)] = _S2593;
    _S2616[int(30)] = _S2594;
    _S2616[int(29)] = _S2595;
    _S2616[int(38)] = _S2586;
    _S2616[int(27)] = _S2597;
    _S2616[int(26)] = _S2598;
    _S2616[int(25)] = _S2599;
    _S2616[int(24)] = _S2600;
    _S2616[int(23)] = _S2585.color_params_0.n_0.y;
    _S2616[int(22)] = _S2585.color_params_0.n_0.x;
    _S2616[int(21)] = _S2585.color_params_0.g_0.y;
    _S2616[int(20)] = _S2585.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S2616;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2617, FixedArray<float, 23>  * _S2618)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S2617, _S2618);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_7, FixedArray<float, 23>  grad_out_3, FixedArray<float, 39>  * _S2619)
{
    FixedArray<float, 39>  _S2620 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_7;
    (&dp_params_3)->differential_0 = _S2620;
    FixedArray<float, 23>  _S2621 = grad_out_3;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, &_S2621);
    *_S2619 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_2, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S2622)
{
    float _S2623;
    FixedArray<float, 6>  losses_5;
    float _S2624 = float(num_cameras_0);
    float _S2625 = raw_losses_2[int(0)] / _S2624;
    for(;;)
    {
        float _S2626 = (F32_abs((_S2625)));
        if(_S2626 < 0.10000000149011612f)
        {
            _S2623 = 0.5f * _S2625 * _S2625 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2623 = _S2626 - 0.05000000074505806f;
            break;
        }
    }
    losses_5[int(0)] = _S2623;
    losses_5[int(1)] = raw_losses_2[int(1)] / (3.0f * _S2624);
    losses_5[int(2)] = (raw_losses_2[int(2)] + raw_losses_2[int(3)] + raw_losses_2[int(4)]) / (9.0f * _S2624);
    losses_5[int(3)] = (raw_losses_2[int(5)] + raw_losses_2[int(6)] + raw_losses_2[int(7)] + raw_losses_2[int(8)] + raw_losses_2[int(9)]) / (5.0f * _S2624);
    float _S2627 = raw_losses_2[int(10)] / _S2624;
    for(;;)
    {
        float _S2628 = (F32_abs((_S2627)));
        if(_S2628 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2627 * _S2627 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2628 - 0.00249999994412065f;
            break;
        }
    }
    float _S2629;
    float _S2630 = raw_losses_2[int(11)] / _S2624;
    for(;;)
    {
        float _S2631 = (F32_abs((_S2630)));
        if(_S2631 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2630 * _S2630 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2631 - 0.00249999994412065f;
            break;
        }
    }
    float _S2632 = _S2623 + _S2629;
    float _S2633 = raw_losses_2[int(12)] / _S2624;
    for(;;)
    {
        float _S2634 = (F32_abs((_S2633)));
        if(_S2634 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2633 * _S2633 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2634 - 0.00249999994412065f;
            break;
        }
    }
    float _S2635 = _S2632 + _S2623;
    float _S2636 = raw_losses_2[int(13)] / _S2624;
    for(;;)
    {
        float _S2637 = (F32_abs((_S2636)));
        if(_S2637 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2636 * _S2636 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2637 - 0.00249999994412065f;
            break;
        }
    }
    float _S2638 = _S2635 + _S2623;
    float _S2639 = raw_losses_2[int(14)] / _S2624;
    for(;;)
    {
        float _S2640 = (F32_abs((_S2639)));
        if(_S2640 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2639 * _S2639 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2640 - 0.00249999994412065f;
            break;
        }
    }
    float _S2641 = _S2638 + _S2623;
    float _S2642 = raw_losses_2[int(15)] / _S2624;
    for(;;)
    {
        float _S2643 = (F32_abs((_S2642)));
        if(_S2643 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2642 * _S2642 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2643 - 0.00249999994412065f;
            break;
        }
    }
    float _S2644 = _S2641 + _S2623;
    float _S2645 = raw_losses_2[int(16)] / _S2624;
    for(;;)
    {
        float _S2646 = (F32_abs((_S2645)));
        if(_S2646 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2645 * _S2645 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2646 - 0.00249999994412065f;
            break;
        }
    }
    float _S2647 = _S2644 + _S2623;
    float _S2648 = raw_losses_2[int(17)] / _S2624;
    for(;;)
    {
        float _S2649 = (F32_abs((_S2648)));
        if(_S2649 < 0.00499999988824129f)
        {
            _S2623 = 0.5f * _S2648 * _S2648 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2623 = _S2649 - 0.00249999994412065f;
            break;
        }
    }
    float _S2650 = (_S2647 + _S2623) / 8.0f;
    float _S2651 = (raw_losses_2[int(18)] + raw_losses_2[int(19)] + raw_losses_2[int(20)] + raw_losses_2[int(21)]) / (4.0f * _S2624);
    losses_5[int(0)] = losses_5[int(0)] * loss_weights_0[int(0)];
    losses_5[int(1)] = losses_5[int(1)] * loss_weights_0[int(1)];
    losses_5[int(2)] = losses_5[int(2)] * loss_weights_0[int(2)];
    losses_5[int(3)] = losses_5[int(3)] * loss_weights_0[int(3)];
    losses_5[int(4)] = _S2650 * loss_weights_0[int(4)];
    losses_5[int(5)] = _S2651 * loss_weights_0[int(5)];
    *_S2622 = losses_5;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_3, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S2652)
{
    float _S2653;
    FixedArray<float, 6>  losses_6;
    float _S2654 = float(num_cameras_1);
    float _S2655 = raw_losses_3[int(0)] / _S2654;
    for(;;)
    {
        float _S2656 = (F32_abs((_S2655)));
        if(_S2656 < 0.10000000149011612f)
        {
            _S2653 = 0.5f * _S2655 * _S2655 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2653 = _S2656 - 0.05000000074505806f;
            break;
        }
    }
    losses_6[int(0)] = _S2653;
    losses_6[int(1)] = raw_losses_3[int(1)] / (3.0f * _S2654);
    losses_6[int(2)] = (raw_losses_3[int(2)] + raw_losses_3[int(3)] + raw_losses_3[int(4)]) / (9.0f * _S2654);
    float _S2657 = 5.0f * _S2654;
    losses_6[int(3)] = (raw_losses_3[int(5)] + raw_losses_3[int(6)] + raw_losses_3[int(7)] + raw_losses_3[int(8)] + raw_losses_3[int(9)]) / _S2657;
    float _S2658 = raw_losses_3[int(10)] / _S2654;
    for(;;)
    {
        float _S2659 = (F32_abs((_S2658)));
        if(_S2659 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2658 * _S2658 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2659 - 0.00249999994412065f;
            break;
        }
    }
    float _S2660;
    float _S2661 = raw_losses_3[int(11)] / _S2654;
    for(;;)
    {
        float _S2662 = (F32_abs((_S2661)));
        if(_S2662 < 0.00499999988824129f)
        {
            _S2660 = 0.5f * _S2661 * _S2661 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2660 = _S2662 - 0.00249999994412065f;
            break;
        }
    }
    float _S2663 = _S2653 + _S2660;
    float _S2664 = raw_losses_3[int(12)] / _S2654;
    for(;;)
    {
        float _S2665 = (F32_abs((_S2664)));
        if(_S2665 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2664 * _S2664 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2665 - 0.00249999994412065f;
            break;
        }
    }
    float _S2666 = _S2663 + _S2653;
    float _S2667 = raw_losses_3[int(13)] / _S2654;
    for(;;)
    {
        float _S2668 = (F32_abs((_S2667)));
        if(_S2668 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2667 * _S2667 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2668 - 0.00249999994412065f;
            break;
        }
    }
    float _S2669 = _S2666 + _S2653;
    float _S2670 = raw_losses_3[int(14)] / _S2654;
    for(;;)
    {
        float _S2671 = (F32_abs((_S2670)));
        if(_S2671 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2670 * _S2670 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2671 - 0.00249999994412065f;
            break;
        }
    }
    float _S2672 = _S2669 + _S2653;
    float _S2673 = raw_losses_3[int(15)] / _S2654;
    for(;;)
    {
        float _S2674 = (F32_abs((_S2673)));
        if(_S2674 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2673 * _S2673 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2674 - 0.00249999994412065f;
            break;
        }
    }
    float _S2675 = _S2672 + _S2653;
    float _S2676 = raw_losses_3[int(16)] / _S2654;
    for(;;)
    {
        float _S2677 = (F32_abs((_S2676)));
        if(_S2677 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2676 * _S2676 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2677 - 0.00249999994412065f;
            break;
        }
    }
    float _S2678 = _S2675 + _S2653;
    float _S2679 = raw_losses_3[int(17)] / _S2654;
    for(;;)
    {
        float _S2680 = (F32_abs((_S2679)));
        if(_S2680 < 0.00499999988824129f)
        {
            _S2653 = 0.5f * _S2679 * _S2679 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2653 = _S2680 - 0.00249999994412065f;
            break;
        }
    }
    float _S2681 = (_S2678 + _S2653) / 8.0f;
    float _S2682 = (raw_losses_3[int(18)] + raw_losses_3[int(19)] + raw_losses_3[int(20)] + raw_losses_3[int(21)] + raw_losses_3[int(22)]) / _S2657;
    losses_6[int(0)] = losses_6[int(0)] * loss_weights_1[int(0)];
    losses_6[int(1)] = losses_6[int(1)] * loss_weights_1[int(1)];
    losses_6[int(2)] = losses_6[int(2)] * loss_weights_1[int(2)];
    losses_6[int(3)] = losses_6[int(3)] * loss_weights_1[int(3)];
    losses_6[int(4)] = _S2681 * loss_weights_1[int(4)];
    losses_6[int(5)] = _S2682 * loss_weights_1[int(5)];
    *_S2652 = losses_6;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_15)
{
    FixedArray<float, 22>  _S2683 = dpraw_losses_1->primal_0;
    float _S2684 = float(num_cameras_2);
    float _S2685 = dpraw_losses_1->primal_0[int(0)] / _S2684;
    bool _S2686 = (s_primal_ctx_abs_0(_S2685)) < 0.10000000149011612f;
    float _S2687;
    if(_S2686)
    {
        _S2687 = 0.5f * _S2685;
    }
    else
    {
        _S2687 = 0.0f;
    }
    float _S2688 = 3.0f * _S2684;
    float _S2689 = 9.0f * _S2684;
    float _S2690 = 5.0f * _S2684;
    float _S2691 = _S2683[int(10)] / _S2684;
    bool _S2692 = (s_primal_ctx_abs_0(_S2691)) < 0.00499999988824129f;
    float _S2693;
    if(_S2692)
    {
        _S2693 = 0.5f * _S2691;
    }
    else
    {
        _S2693 = 0.0f;
    }
    float _S2694 = _S2683[int(11)] / _S2684;
    bool _S2695 = (s_primal_ctx_abs_0(_S2694)) < 0.00499999988824129f;
    float _S2696;
    if(_S2695)
    {
        _S2696 = 0.5f * _S2694;
    }
    else
    {
        _S2696 = 0.0f;
    }
    float _S2697 = _S2683[int(12)] / _S2684;
    bool _S2698 = (s_primal_ctx_abs_0(_S2697)) < 0.00499999988824129f;
    float _S2699;
    if(_S2698)
    {
        _S2699 = 0.5f * _S2697;
    }
    else
    {
        _S2699 = 0.0f;
    }
    float _S2700 = _S2683[int(13)] / _S2684;
    bool _S2701 = (s_primal_ctx_abs_0(_S2700)) < 0.00499999988824129f;
    float _S2702;
    if(_S2701)
    {
        _S2702 = 0.5f * _S2700;
    }
    else
    {
        _S2702 = 0.0f;
    }
    float _S2703 = _S2683[int(14)] / _S2684;
    bool _S2704 = (s_primal_ctx_abs_0(_S2703)) < 0.00499999988824129f;
    float _S2705;
    if(_S2704)
    {
        _S2705 = 0.5f * _S2703;
    }
    else
    {
        _S2705 = 0.0f;
    }
    float _S2706 = _S2683[int(15)] / _S2684;
    bool _S2707 = (s_primal_ctx_abs_0(_S2706)) < 0.00499999988824129f;
    float _S2708;
    if(_S2707)
    {
        _S2708 = 0.5f * _S2706;
    }
    else
    {
        _S2708 = 0.0f;
    }
    float _S2709 = _S2683[int(16)] / _S2684;
    bool _S2710 = (s_primal_ctx_abs_0(_S2709)) < 0.00499999988824129f;
    float _S2711;
    if(_S2710)
    {
        _S2711 = 0.5f * _S2709;
    }
    else
    {
        _S2711 = 0.0f;
    }
    float _S2712 = _S2683[int(17)] / _S2684;
    bool _S2713 = (s_primal_ctx_abs_0(_S2712)) < 0.00499999988824129f;
    float _S2714;
    if(_S2713)
    {
        _S2714 = 0.5f * _S2712;
    }
    else
    {
        _S2714 = 0.0f;
    }
    float _S2715 = (*loss_weights_2)[int(3)] * (*_s_dOut_15)[int(3)];
    float _S2716 = (*loss_weights_2)[int(2)] * (*_s_dOut_15)[int(2)];
    float _S2717 = (*loss_weights_2)[int(1)] * (*_s_dOut_15)[int(1)];
    float _S2718 = (*loss_weights_2)[int(0)] * (*_s_dOut_15)[int(0)];
    float _S2719 = (*loss_weights_2)[int(5)] * (*_s_dOut_15)[int(5)] / (4.0f * _S2684);
    float _S2720 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_15)[int(4)]);
    FixedArray<float, 22>  _S2721;
    _S2721[int(0)] = 0.0f;
    _S2721[int(1)] = 0.0f;
    _S2721[int(2)] = 0.0f;
    _S2721[int(3)] = 0.0f;
    _S2721[int(4)] = 0.0f;
    _S2721[int(5)] = 0.0f;
    _S2721[int(6)] = 0.0f;
    _S2721[int(7)] = 0.0f;
    _S2721[int(8)] = 0.0f;
    _S2721[int(9)] = 0.0f;
    _S2721[int(10)] = 0.0f;
    _S2721[int(11)] = 0.0f;
    _S2721[int(12)] = 0.0f;
    _S2721[int(13)] = 0.0f;
    _S2721[int(14)] = 0.0f;
    _S2721[int(15)] = 0.0f;
    _S2721[int(16)] = 0.0f;
    _S2721[int(17)] = 0.0f;
    _S2721[int(18)] = 0.0f;
    _S2721[int(19)] = 0.0f;
    _S2721[int(20)] = 0.0f;
    _S2721[int(21)] = 0.0f;
    _S2721[int(21)] = _S2719;
    _S2721[int(20)] = _S2719;
    _S2721[int(19)] = _S2719;
    _S2721[int(18)] = _S2719;
    float _S2722 = _S2721[int(0)];
    float _S2723 = _S2721[int(1)];
    float _S2724 = _S2721[int(2)];
    float _S2725 = _S2721[int(3)];
    float _S2726 = _S2721[int(4)];
    float _S2727 = _S2721[int(5)];
    float _S2728 = _S2721[int(6)];
    float _S2729 = _S2721[int(7)];
    float _S2730 = _S2721[int(8)];
    float _S2731 = _S2721[int(9)];
    float _S2732 = _S2721[int(10)];
    float _S2733 = _S2721[int(11)];
    float _S2734 = _S2721[int(12)];
    float _S2735 = _S2721[int(13)];
    float _S2736 = _S2721[int(14)];
    float _S2737 = _S2721[int(15)];
    float _S2738 = _S2721[int(16)];
    float _S2739 = _S2721[int(17)];
    float _S2740 = _S2721[int(18)];
    float _S2741 = _S2721[int(19)];
    float _S2742 = _S2721[int(20)];
    float _S2743 = _S2721[int(21)];
    float _S2744;
    if(_S2713)
    {
        float _S2745 = 200.0f * _S2720;
        float _S2746 = _S2714 * _S2745 + 0.5f * (_S2712 * _S2745);
        _S2714 = 0.0f;
        _S2744 = _S2746;
    }
    else
    {
        _S2714 = _S2720;
        _S2744 = 0.0f;
    }
    DiffPair_float_0 _S2747;
    (&_S2747)->primal_0 = _S2712;
    (&_S2747)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2747, _S2714);
    float _S2748 = (_S2747.differential_0 + _S2744) / _S2684;
    FixedArray<float, 22>  _S2749;
    _S2749[int(0)] = 0.0f;
    _S2749[int(1)] = 0.0f;
    _S2749[int(2)] = 0.0f;
    _S2749[int(3)] = 0.0f;
    _S2749[int(4)] = 0.0f;
    _S2749[int(5)] = 0.0f;
    _S2749[int(6)] = 0.0f;
    _S2749[int(7)] = 0.0f;
    _S2749[int(8)] = 0.0f;
    _S2749[int(9)] = 0.0f;
    _S2749[int(10)] = 0.0f;
    _S2749[int(11)] = 0.0f;
    _S2749[int(12)] = 0.0f;
    _S2749[int(13)] = 0.0f;
    _S2749[int(14)] = 0.0f;
    _S2749[int(15)] = 0.0f;
    _S2749[int(16)] = 0.0f;
    _S2749[int(17)] = 0.0f;
    _S2749[int(18)] = 0.0f;
    _S2749[int(19)] = 0.0f;
    _S2749[int(20)] = 0.0f;
    _S2749[int(21)] = 0.0f;
    _S2749[int(17)] = _S2748;
    float _S2750 = _S2722 + _S2749[int(0)];
    float _S2751 = _S2723 + _S2749[int(1)];
    float _S2752 = _S2724 + _S2749[int(2)];
    float _S2753 = _S2725 + _S2749[int(3)];
    float _S2754 = _S2726 + _S2749[int(4)];
    float _S2755 = _S2727 + _S2749[int(5)];
    float _S2756 = _S2728 + _S2749[int(6)];
    float _S2757 = _S2729 + _S2749[int(7)];
    float _S2758 = _S2730 + _S2749[int(8)];
    float _S2759 = _S2731 + _S2749[int(9)];
    float _S2760 = _S2732 + _S2749[int(10)];
    float _S2761 = _S2733 + _S2749[int(11)];
    float _S2762 = _S2734 + _S2749[int(12)];
    float _S2763 = _S2735 + _S2749[int(13)];
    float _S2764 = _S2736 + _S2749[int(14)];
    float _S2765 = _S2737 + _S2749[int(15)];
    float _S2766 = _S2738 + _S2749[int(16)];
    float _S2767 = _S2739 + _S2749[int(17)];
    float _S2768 = _S2740 + _S2749[int(18)];
    float _S2769 = _S2741 + _S2749[int(19)];
    float _S2770 = _S2742 + _S2749[int(20)];
    float _S2771 = _S2743 + _S2749[int(21)];
    if(_S2710)
    {
        float _S2772 = 200.0f * _S2720;
        float _S2773 = _S2711 * _S2772 + 0.5f * (_S2709 * _S2772);
        _S2711 = 0.0f;
        _S2714 = _S2773;
    }
    else
    {
        _S2711 = _S2720;
        _S2714 = 0.0f;
    }
    DiffPair_float_0 _S2774;
    (&_S2774)->primal_0 = _S2709;
    (&_S2774)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2774, _S2711);
    float _S2775 = (_S2774.differential_0 + _S2714) / _S2684;
    FixedArray<float, 22>  _S2776;
    _S2776[int(0)] = 0.0f;
    _S2776[int(1)] = 0.0f;
    _S2776[int(2)] = 0.0f;
    _S2776[int(3)] = 0.0f;
    _S2776[int(4)] = 0.0f;
    _S2776[int(5)] = 0.0f;
    _S2776[int(6)] = 0.0f;
    _S2776[int(7)] = 0.0f;
    _S2776[int(8)] = 0.0f;
    _S2776[int(9)] = 0.0f;
    _S2776[int(10)] = 0.0f;
    _S2776[int(11)] = 0.0f;
    _S2776[int(12)] = 0.0f;
    _S2776[int(13)] = 0.0f;
    _S2776[int(14)] = 0.0f;
    _S2776[int(15)] = 0.0f;
    _S2776[int(16)] = 0.0f;
    _S2776[int(17)] = 0.0f;
    _S2776[int(18)] = 0.0f;
    _S2776[int(19)] = 0.0f;
    _S2776[int(20)] = 0.0f;
    _S2776[int(21)] = 0.0f;
    _S2776[int(16)] = _S2775;
    float _S2777 = _S2750 + _S2776[int(0)];
    float _S2778 = _S2751 + _S2776[int(1)];
    float _S2779 = _S2752 + _S2776[int(2)];
    float _S2780 = _S2753 + _S2776[int(3)];
    float _S2781 = _S2754 + _S2776[int(4)];
    float _S2782 = _S2755 + _S2776[int(5)];
    float _S2783 = _S2756 + _S2776[int(6)];
    float _S2784 = _S2757 + _S2776[int(7)];
    float _S2785 = _S2758 + _S2776[int(8)];
    float _S2786 = _S2759 + _S2776[int(9)];
    float _S2787 = _S2760 + _S2776[int(10)];
    float _S2788 = _S2761 + _S2776[int(11)];
    float _S2789 = _S2762 + _S2776[int(12)];
    float _S2790 = _S2763 + _S2776[int(13)];
    float _S2791 = _S2764 + _S2776[int(14)];
    float _S2792 = _S2765 + _S2776[int(15)];
    float _S2793 = _S2766 + _S2776[int(16)];
    float _S2794 = _S2767 + _S2776[int(17)];
    float _S2795 = _S2768 + _S2776[int(18)];
    float _S2796 = _S2769 + _S2776[int(19)];
    float _S2797 = _S2770 + _S2776[int(20)];
    float _S2798 = _S2771 + _S2776[int(21)];
    if(_S2707)
    {
        float _S2799 = 200.0f * _S2720;
        float _S2800 = _S2708 * _S2799 + 0.5f * (_S2706 * _S2799);
        _S2708 = 0.0f;
        _S2711 = _S2800;
    }
    else
    {
        _S2708 = _S2720;
        _S2711 = 0.0f;
    }
    DiffPair_float_0 _S2801;
    (&_S2801)->primal_0 = _S2706;
    (&_S2801)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2801, _S2708);
    float _S2802 = (_S2801.differential_0 + _S2711) / _S2684;
    FixedArray<float, 22>  _S2803;
    _S2803[int(0)] = 0.0f;
    _S2803[int(1)] = 0.0f;
    _S2803[int(2)] = 0.0f;
    _S2803[int(3)] = 0.0f;
    _S2803[int(4)] = 0.0f;
    _S2803[int(5)] = 0.0f;
    _S2803[int(6)] = 0.0f;
    _S2803[int(7)] = 0.0f;
    _S2803[int(8)] = 0.0f;
    _S2803[int(9)] = 0.0f;
    _S2803[int(10)] = 0.0f;
    _S2803[int(11)] = 0.0f;
    _S2803[int(12)] = 0.0f;
    _S2803[int(13)] = 0.0f;
    _S2803[int(14)] = 0.0f;
    _S2803[int(15)] = 0.0f;
    _S2803[int(16)] = 0.0f;
    _S2803[int(17)] = 0.0f;
    _S2803[int(18)] = 0.0f;
    _S2803[int(19)] = 0.0f;
    _S2803[int(20)] = 0.0f;
    _S2803[int(21)] = 0.0f;
    _S2803[int(15)] = _S2802;
    float _S2804 = _S2777 + _S2803[int(0)];
    float _S2805 = _S2778 + _S2803[int(1)];
    float _S2806 = _S2779 + _S2803[int(2)];
    float _S2807 = _S2780 + _S2803[int(3)];
    float _S2808 = _S2781 + _S2803[int(4)];
    float _S2809 = _S2782 + _S2803[int(5)];
    float _S2810 = _S2783 + _S2803[int(6)];
    float _S2811 = _S2784 + _S2803[int(7)];
    float _S2812 = _S2785 + _S2803[int(8)];
    float _S2813 = _S2786 + _S2803[int(9)];
    float _S2814 = _S2787 + _S2803[int(10)];
    float _S2815 = _S2788 + _S2803[int(11)];
    float _S2816 = _S2789 + _S2803[int(12)];
    float _S2817 = _S2790 + _S2803[int(13)];
    float _S2818 = _S2791 + _S2803[int(14)];
    float _S2819 = _S2792 + _S2803[int(15)];
    float _S2820 = _S2793 + _S2803[int(16)];
    float _S2821 = _S2794 + _S2803[int(17)];
    float _S2822 = _S2795 + _S2803[int(18)];
    float _S2823 = _S2796 + _S2803[int(19)];
    float _S2824 = _S2797 + _S2803[int(20)];
    float _S2825 = _S2798 + _S2803[int(21)];
    if(_S2704)
    {
        float _S2826 = 200.0f * _S2720;
        float _S2827 = _S2705 * _S2826 + 0.5f * (_S2703 * _S2826);
        _S2705 = 0.0f;
        _S2708 = _S2827;
    }
    else
    {
        _S2705 = _S2720;
        _S2708 = 0.0f;
    }
    DiffPair_float_0 _S2828;
    (&_S2828)->primal_0 = _S2703;
    (&_S2828)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2828, _S2705);
    float _S2829 = (_S2828.differential_0 + _S2708) / _S2684;
    FixedArray<float, 22>  _S2830;
    _S2830[int(0)] = 0.0f;
    _S2830[int(1)] = 0.0f;
    _S2830[int(2)] = 0.0f;
    _S2830[int(3)] = 0.0f;
    _S2830[int(4)] = 0.0f;
    _S2830[int(5)] = 0.0f;
    _S2830[int(6)] = 0.0f;
    _S2830[int(7)] = 0.0f;
    _S2830[int(8)] = 0.0f;
    _S2830[int(9)] = 0.0f;
    _S2830[int(10)] = 0.0f;
    _S2830[int(11)] = 0.0f;
    _S2830[int(12)] = 0.0f;
    _S2830[int(13)] = 0.0f;
    _S2830[int(14)] = 0.0f;
    _S2830[int(15)] = 0.0f;
    _S2830[int(16)] = 0.0f;
    _S2830[int(17)] = 0.0f;
    _S2830[int(18)] = 0.0f;
    _S2830[int(19)] = 0.0f;
    _S2830[int(20)] = 0.0f;
    _S2830[int(21)] = 0.0f;
    _S2830[int(14)] = _S2829;
    float _S2831 = _S2804 + _S2830[int(0)];
    float _S2832 = _S2805 + _S2830[int(1)];
    float _S2833 = _S2806 + _S2830[int(2)];
    float _S2834 = _S2807 + _S2830[int(3)];
    float _S2835 = _S2808 + _S2830[int(4)];
    float _S2836 = _S2809 + _S2830[int(5)];
    float _S2837 = _S2810 + _S2830[int(6)];
    float _S2838 = _S2811 + _S2830[int(7)];
    float _S2839 = _S2812 + _S2830[int(8)];
    float _S2840 = _S2813 + _S2830[int(9)];
    float _S2841 = _S2814 + _S2830[int(10)];
    float _S2842 = _S2815 + _S2830[int(11)];
    float _S2843 = _S2816 + _S2830[int(12)];
    float _S2844 = _S2817 + _S2830[int(13)];
    float _S2845 = _S2818 + _S2830[int(14)];
    float _S2846 = _S2819 + _S2830[int(15)];
    float _S2847 = _S2820 + _S2830[int(16)];
    float _S2848 = _S2821 + _S2830[int(17)];
    float _S2849 = _S2822 + _S2830[int(18)];
    float _S2850 = _S2823 + _S2830[int(19)];
    float _S2851 = _S2824 + _S2830[int(20)];
    float _S2852 = _S2825 + _S2830[int(21)];
    if(_S2701)
    {
        float _S2853 = 200.0f * _S2720;
        float _S2854 = _S2702 * _S2853 + 0.5f * (_S2700 * _S2853);
        _S2702 = 0.0f;
        _S2705 = _S2854;
    }
    else
    {
        _S2702 = _S2720;
        _S2705 = 0.0f;
    }
    DiffPair_float_0 _S2855;
    (&_S2855)->primal_0 = _S2700;
    (&_S2855)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2855, _S2702);
    float _S2856 = (_S2855.differential_0 + _S2705) / _S2684;
    FixedArray<float, 22>  _S2857;
    _S2857[int(0)] = 0.0f;
    _S2857[int(1)] = 0.0f;
    _S2857[int(2)] = 0.0f;
    _S2857[int(3)] = 0.0f;
    _S2857[int(4)] = 0.0f;
    _S2857[int(5)] = 0.0f;
    _S2857[int(6)] = 0.0f;
    _S2857[int(7)] = 0.0f;
    _S2857[int(8)] = 0.0f;
    _S2857[int(9)] = 0.0f;
    _S2857[int(10)] = 0.0f;
    _S2857[int(11)] = 0.0f;
    _S2857[int(12)] = 0.0f;
    _S2857[int(13)] = 0.0f;
    _S2857[int(14)] = 0.0f;
    _S2857[int(15)] = 0.0f;
    _S2857[int(16)] = 0.0f;
    _S2857[int(17)] = 0.0f;
    _S2857[int(18)] = 0.0f;
    _S2857[int(19)] = 0.0f;
    _S2857[int(20)] = 0.0f;
    _S2857[int(21)] = 0.0f;
    _S2857[int(13)] = _S2856;
    float _S2858 = _S2831 + _S2857[int(0)];
    float _S2859 = _S2832 + _S2857[int(1)];
    float _S2860 = _S2833 + _S2857[int(2)];
    float _S2861 = _S2834 + _S2857[int(3)];
    float _S2862 = _S2835 + _S2857[int(4)];
    float _S2863 = _S2836 + _S2857[int(5)];
    float _S2864 = _S2837 + _S2857[int(6)];
    float _S2865 = _S2838 + _S2857[int(7)];
    float _S2866 = _S2839 + _S2857[int(8)];
    float _S2867 = _S2840 + _S2857[int(9)];
    float _S2868 = _S2841 + _S2857[int(10)];
    float _S2869 = _S2842 + _S2857[int(11)];
    float _S2870 = _S2843 + _S2857[int(12)];
    float _S2871 = _S2844 + _S2857[int(13)];
    float _S2872 = _S2845 + _S2857[int(14)];
    float _S2873 = _S2846 + _S2857[int(15)];
    float _S2874 = _S2847 + _S2857[int(16)];
    float _S2875 = _S2848 + _S2857[int(17)];
    float _S2876 = _S2849 + _S2857[int(18)];
    float _S2877 = _S2850 + _S2857[int(19)];
    float _S2878 = _S2851 + _S2857[int(20)];
    float _S2879 = _S2852 + _S2857[int(21)];
    if(_S2698)
    {
        float _S2880 = 200.0f * _S2720;
        float _S2881 = _S2699 * _S2880 + 0.5f * (_S2697 * _S2880);
        _S2699 = 0.0f;
        _S2702 = _S2881;
    }
    else
    {
        _S2699 = _S2720;
        _S2702 = 0.0f;
    }
    DiffPair_float_0 _S2882;
    (&_S2882)->primal_0 = _S2697;
    (&_S2882)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2882, _S2699);
    float _S2883 = (_S2882.differential_0 + _S2702) / _S2684;
    FixedArray<float, 22>  _S2884;
    _S2884[int(0)] = 0.0f;
    _S2884[int(1)] = 0.0f;
    _S2884[int(2)] = 0.0f;
    _S2884[int(3)] = 0.0f;
    _S2884[int(4)] = 0.0f;
    _S2884[int(5)] = 0.0f;
    _S2884[int(6)] = 0.0f;
    _S2884[int(7)] = 0.0f;
    _S2884[int(8)] = 0.0f;
    _S2884[int(9)] = 0.0f;
    _S2884[int(10)] = 0.0f;
    _S2884[int(11)] = 0.0f;
    _S2884[int(12)] = 0.0f;
    _S2884[int(13)] = 0.0f;
    _S2884[int(14)] = 0.0f;
    _S2884[int(15)] = 0.0f;
    _S2884[int(16)] = 0.0f;
    _S2884[int(17)] = 0.0f;
    _S2884[int(18)] = 0.0f;
    _S2884[int(19)] = 0.0f;
    _S2884[int(20)] = 0.0f;
    _S2884[int(21)] = 0.0f;
    _S2884[int(12)] = _S2883;
    float _S2885 = _S2858 + _S2884[int(0)];
    float _S2886 = _S2859 + _S2884[int(1)];
    float _S2887 = _S2860 + _S2884[int(2)];
    float _S2888 = _S2861 + _S2884[int(3)];
    float _S2889 = _S2862 + _S2884[int(4)];
    float _S2890 = _S2863 + _S2884[int(5)];
    float _S2891 = _S2864 + _S2884[int(6)];
    float _S2892 = _S2865 + _S2884[int(7)];
    float _S2893 = _S2866 + _S2884[int(8)];
    float _S2894 = _S2867 + _S2884[int(9)];
    float _S2895 = _S2868 + _S2884[int(10)];
    float _S2896 = _S2869 + _S2884[int(11)];
    float _S2897 = _S2870 + _S2884[int(12)];
    float _S2898 = _S2871 + _S2884[int(13)];
    float _S2899 = _S2872 + _S2884[int(14)];
    float _S2900 = _S2873 + _S2884[int(15)];
    float _S2901 = _S2874 + _S2884[int(16)];
    float _S2902 = _S2875 + _S2884[int(17)];
    float _S2903 = _S2876 + _S2884[int(18)];
    float _S2904 = _S2877 + _S2884[int(19)];
    float _S2905 = _S2878 + _S2884[int(20)];
    float _S2906 = _S2879 + _S2884[int(21)];
    if(_S2695)
    {
        float _S2907 = 200.0f * _S2720;
        float _S2908 = _S2696 * _S2907 + 0.5f * (_S2694 * _S2907);
        _S2696 = 0.0f;
        _S2699 = _S2908;
    }
    else
    {
        _S2696 = _S2720;
        _S2699 = 0.0f;
    }
    DiffPair_float_0 _S2909;
    (&_S2909)->primal_0 = _S2694;
    (&_S2909)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2909, _S2696);
    float _S2910 = (_S2909.differential_0 + _S2699) / _S2684;
    FixedArray<float, 22>  _S2911;
    _S2911[int(0)] = 0.0f;
    _S2911[int(1)] = 0.0f;
    _S2911[int(2)] = 0.0f;
    _S2911[int(3)] = 0.0f;
    _S2911[int(4)] = 0.0f;
    _S2911[int(5)] = 0.0f;
    _S2911[int(6)] = 0.0f;
    _S2911[int(7)] = 0.0f;
    _S2911[int(8)] = 0.0f;
    _S2911[int(9)] = 0.0f;
    _S2911[int(10)] = 0.0f;
    _S2911[int(11)] = 0.0f;
    _S2911[int(12)] = 0.0f;
    _S2911[int(13)] = 0.0f;
    _S2911[int(14)] = 0.0f;
    _S2911[int(15)] = 0.0f;
    _S2911[int(16)] = 0.0f;
    _S2911[int(17)] = 0.0f;
    _S2911[int(18)] = 0.0f;
    _S2911[int(19)] = 0.0f;
    _S2911[int(20)] = 0.0f;
    _S2911[int(21)] = 0.0f;
    _S2911[int(11)] = _S2910;
    float _S2912 = _S2885 + _S2911[int(0)];
    float _S2913 = _S2886 + _S2911[int(1)];
    float _S2914 = _S2887 + _S2911[int(2)];
    float _S2915 = _S2888 + _S2911[int(3)];
    float _S2916 = _S2889 + _S2911[int(4)];
    float _S2917 = _S2890 + _S2911[int(5)];
    float _S2918 = _S2891 + _S2911[int(6)];
    float _S2919 = _S2892 + _S2911[int(7)];
    float _S2920 = _S2893 + _S2911[int(8)];
    float _S2921 = _S2894 + _S2911[int(9)];
    float _S2922 = _S2895 + _S2911[int(10)];
    float _S2923 = _S2896 + _S2911[int(11)];
    float _S2924 = _S2897 + _S2911[int(12)];
    float _S2925 = _S2898 + _S2911[int(13)];
    float _S2926 = _S2899 + _S2911[int(14)];
    float _S2927 = _S2900 + _S2911[int(15)];
    float _S2928 = _S2901 + _S2911[int(16)];
    float _S2929 = _S2902 + _S2911[int(17)];
    float _S2930 = _S2903 + _S2911[int(18)];
    float _S2931 = _S2904 + _S2911[int(19)];
    float _S2932 = _S2905 + _S2911[int(20)];
    float _S2933 = _S2906 + _S2911[int(21)];
    if(_S2692)
    {
        float _S2934 = 200.0f * _S2720;
        float _S2935 = _S2693 * _S2934 + 0.5f * (_S2691 * _S2934);
        _S2693 = 0.0f;
        _S2696 = _S2935;
    }
    else
    {
        _S2693 = _S2720;
        _S2696 = 0.0f;
    }
    DiffPair_float_0 _S2936;
    (&_S2936)->primal_0 = _S2691;
    (&_S2936)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2936, _S2693);
    float _S2937 = (_S2936.differential_0 + _S2696) / _S2684;
    float _S2938 = _S2715 / _S2690;
    float _S2939 = _S2716 / _S2689;
    float _S2940 = _S2717 / _S2688;
    FixedArray<float, 22>  _S2941;
    _S2941[int(0)] = 0.0f;
    _S2941[int(1)] = 0.0f;
    _S2941[int(2)] = 0.0f;
    _S2941[int(3)] = 0.0f;
    _S2941[int(4)] = 0.0f;
    _S2941[int(5)] = 0.0f;
    _S2941[int(6)] = 0.0f;
    _S2941[int(7)] = 0.0f;
    _S2941[int(8)] = 0.0f;
    _S2941[int(9)] = 0.0f;
    _S2941[int(10)] = 0.0f;
    _S2941[int(11)] = 0.0f;
    _S2941[int(12)] = 0.0f;
    _S2941[int(13)] = 0.0f;
    _S2941[int(14)] = 0.0f;
    _S2941[int(15)] = 0.0f;
    _S2941[int(16)] = 0.0f;
    _S2941[int(17)] = 0.0f;
    _S2941[int(18)] = 0.0f;
    _S2941[int(19)] = 0.0f;
    _S2941[int(20)] = 0.0f;
    _S2941[int(21)] = 0.0f;
    _S2941[int(10)] = _S2937;
    _S2941[int(9)] = _S2938;
    _S2941[int(8)] = _S2938;
    _S2941[int(7)] = _S2938;
    _S2941[int(6)] = _S2938;
    _S2941[int(5)] = _S2938;
    _S2941[int(4)] = _S2939;
    _S2941[int(3)] = _S2939;
    _S2941[int(2)] = _S2939;
    _S2941[int(1)] = _S2940;
    float _S2942 = _S2912 + _S2941[int(0)];
    float _S2943 = _S2913 + _S2941[int(1)];
    float _S2944 = _S2914 + _S2941[int(2)];
    float _S2945 = _S2915 + _S2941[int(3)];
    float _S2946 = _S2916 + _S2941[int(4)];
    float _S2947 = _S2917 + _S2941[int(5)];
    float _S2948 = _S2918 + _S2941[int(6)];
    float _S2949 = _S2919 + _S2941[int(7)];
    float _S2950 = _S2920 + _S2941[int(8)];
    float _S2951 = _S2921 + _S2941[int(9)];
    float _S2952 = _S2922 + _S2941[int(10)];
    float _S2953 = _S2923 + _S2941[int(11)];
    float _S2954 = _S2924 + _S2941[int(12)];
    float _S2955 = _S2925 + _S2941[int(13)];
    float _S2956 = _S2926 + _S2941[int(14)];
    float _S2957 = _S2927 + _S2941[int(15)];
    float _S2958 = _S2928 + _S2941[int(16)];
    float _S2959 = _S2929 + _S2941[int(17)];
    float _S2960 = _S2930 + _S2941[int(18)];
    float _S2961 = _S2931 + _S2941[int(19)];
    float _S2962 = _S2932 + _S2941[int(20)];
    float _S2963 = _S2933 + _S2941[int(21)];
    if(_S2686)
    {
        float _S2964 = 10.0f * _S2718;
        float _S2965 = _S2687 * _S2964 + 0.5f * (_S2685 * _S2964);
        _S2687 = 0.0f;
        _S2693 = _S2965;
    }
    else
    {
        _S2687 = _S2718;
        _S2693 = 0.0f;
    }
    DiffPair_float_0 _S2966;
    (&_S2966)->primal_0 = _S2685;
    (&_S2966)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2966, _S2687);
    float _S2967 = (_S2966.differential_0 + _S2693) / _S2684;
    FixedArray<float, 22>  _S2968;
    _S2968[int(0)] = 0.0f;
    _S2968[int(1)] = 0.0f;
    _S2968[int(2)] = 0.0f;
    _S2968[int(3)] = 0.0f;
    _S2968[int(4)] = 0.0f;
    _S2968[int(5)] = 0.0f;
    _S2968[int(6)] = 0.0f;
    _S2968[int(7)] = 0.0f;
    _S2968[int(8)] = 0.0f;
    _S2968[int(9)] = 0.0f;
    _S2968[int(10)] = 0.0f;
    _S2968[int(11)] = 0.0f;
    _S2968[int(12)] = 0.0f;
    _S2968[int(13)] = 0.0f;
    _S2968[int(14)] = 0.0f;
    _S2968[int(15)] = 0.0f;
    _S2968[int(16)] = 0.0f;
    _S2968[int(17)] = 0.0f;
    _S2968[int(18)] = 0.0f;
    _S2968[int(19)] = 0.0f;
    _S2968[int(20)] = 0.0f;
    _S2968[int(21)] = 0.0f;
    _S2968[int(0)] = _S2967;
    FixedArray<float, 22>  _S2969 = {
        _S2942 + _S2968[int(0)], _S2943 + _S2968[int(1)], _S2944 + _S2968[int(2)], _S2945 + _S2968[int(3)], _S2946 + _S2968[int(4)], _S2947 + _S2968[int(5)], _S2948 + _S2968[int(6)], _S2949 + _S2968[int(7)], _S2950 + _S2968[int(8)], _S2951 + _S2968[int(9)], _S2952 + _S2968[int(10)], _S2953 + _S2968[int(11)], _S2954 + _S2968[int(12)], _S2955 + _S2968[int(13)], _S2956 + _S2968[int(14)], _S2957 + _S2968[int(15)], _S2958 + _S2968[int(16)], _S2959 + _S2968[int(17)], _S2960 + _S2968[int(18)], _S2961 + _S2968[int(19)], _S2962 + _S2968[int(20)], _S2963 + _S2968[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2969;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S2970, int _S2971, FixedArray<float, 6>  * _S2972, FixedArray<float, 6>  * _S2973)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2970, _S2971, _S2972, _S2973);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_4, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_4, FixedArray<float, 22>  * _S2974)
{
    FixedArray<float, 22>  _S2975 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_4;
    (&dp_raw_losses_1)->differential_0 = _S2975;
    FixedArray<float, 6>  _S2976 = loss_weights_3;
    FixedArray<float, 6>  _S2977 = grad_out_4;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_3, &_S2976, &_S2977);
    *_S2974 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_2, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_16)
{
    FixedArray<float, 23>  _S2978 = dpraw_losses_2->primal_0;
    float _S2979 = float(num_cameras_4);
    float _S2980 = dpraw_losses_2->primal_0[int(0)] / _S2979;
    bool _S2981 = (s_primal_ctx_abs_0(_S2980)) < 0.10000000149011612f;
    float _S2982;
    if(_S2981)
    {
        _S2982 = 0.5f * _S2980;
    }
    else
    {
        _S2982 = 0.0f;
    }
    float _S2983 = 3.0f * _S2979;
    float _S2984 = 9.0f * _S2979;
    float _S2985 = 5.0f * _S2979;
    float _S2986 = _S2978[int(10)] / _S2979;
    bool _S2987 = (s_primal_ctx_abs_0(_S2986)) < 0.00499999988824129f;
    float _S2988;
    if(_S2987)
    {
        _S2988 = 0.5f * _S2986;
    }
    else
    {
        _S2988 = 0.0f;
    }
    float _S2989 = _S2978[int(11)] / _S2979;
    bool _S2990 = (s_primal_ctx_abs_0(_S2989)) < 0.00499999988824129f;
    float _S2991;
    if(_S2990)
    {
        _S2991 = 0.5f * _S2989;
    }
    else
    {
        _S2991 = 0.0f;
    }
    float _S2992 = _S2978[int(12)] / _S2979;
    bool _S2993 = (s_primal_ctx_abs_0(_S2992)) < 0.00499999988824129f;
    float _S2994;
    if(_S2993)
    {
        _S2994 = 0.5f * _S2992;
    }
    else
    {
        _S2994 = 0.0f;
    }
    float _S2995 = _S2978[int(13)] / _S2979;
    bool _S2996 = (s_primal_ctx_abs_0(_S2995)) < 0.00499999988824129f;
    float _S2997;
    if(_S2996)
    {
        _S2997 = 0.5f * _S2995;
    }
    else
    {
        _S2997 = 0.0f;
    }
    float _S2998 = _S2978[int(14)] / _S2979;
    bool _S2999 = (s_primal_ctx_abs_0(_S2998)) < 0.00499999988824129f;
    float _S3000;
    if(_S2999)
    {
        _S3000 = 0.5f * _S2998;
    }
    else
    {
        _S3000 = 0.0f;
    }
    float _S3001 = _S2978[int(15)] / _S2979;
    bool _S3002 = (s_primal_ctx_abs_0(_S3001)) < 0.00499999988824129f;
    float _S3003;
    if(_S3002)
    {
        _S3003 = 0.5f * _S3001;
    }
    else
    {
        _S3003 = 0.0f;
    }
    float _S3004 = _S2978[int(16)] / _S2979;
    bool _S3005 = (s_primal_ctx_abs_0(_S3004)) < 0.00499999988824129f;
    float _S3006;
    if(_S3005)
    {
        _S3006 = 0.5f * _S3004;
    }
    else
    {
        _S3006 = 0.0f;
    }
    float _S3007 = _S2978[int(17)] / _S2979;
    bool _S3008 = (s_primal_ctx_abs_0(_S3007)) < 0.00499999988824129f;
    float _S3009;
    if(_S3008)
    {
        _S3009 = 0.5f * _S3007;
    }
    else
    {
        _S3009 = 0.0f;
    }
    float _S3010 = (*loss_weights_4)[int(3)] * (*_s_dOut_16)[int(3)];
    float _S3011 = (*loss_weights_4)[int(2)] * (*_s_dOut_16)[int(2)];
    float _S3012 = (*loss_weights_4)[int(1)] * (*_s_dOut_16)[int(1)];
    float _S3013 = (*loss_weights_4)[int(0)] * (*_s_dOut_16)[int(0)];
    float _S3014 = (*loss_weights_4)[int(5)] * (*_s_dOut_16)[int(5)] / _S2985;
    float _S3015 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_16)[int(4)]);
    FixedArray<float, 23>  _S3016;
    _S3016[int(0)] = 0.0f;
    _S3016[int(1)] = 0.0f;
    _S3016[int(2)] = 0.0f;
    _S3016[int(3)] = 0.0f;
    _S3016[int(4)] = 0.0f;
    _S3016[int(5)] = 0.0f;
    _S3016[int(6)] = 0.0f;
    _S3016[int(7)] = 0.0f;
    _S3016[int(8)] = 0.0f;
    _S3016[int(9)] = 0.0f;
    _S3016[int(10)] = 0.0f;
    _S3016[int(11)] = 0.0f;
    _S3016[int(12)] = 0.0f;
    _S3016[int(13)] = 0.0f;
    _S3016[int(14)] = 0.0f;
    _S3016[int(15)] = 0.0f;
    _S3016[int(16)] = 0.0f;
    _S3016[int(17)] = 0.0f;
    _S3016[int(18)] = 0.0f;
    _S3016[int(19)] = 0.0f;
    _S3016[int(20)] = 0.0f;
    _S3016[int(21)] = 0.0f;
    _S3016[int(22)] = 0.0f;
    _S3016[int(22)] = _S3014;
    _S3016[int(21)] = _S3014;
    _S3016[int(20)] = _S3014;
    _S3016[int(19)] = _S3014;
    _S3016[int(18)] = _S3014;
    float _S3017 = _S3016[int(0)];
    float _S3018 = _S3016[int(1)];
    float _S3019 = _S3016[int(2)];
    float _S3020 = _S3016[int(3)];
    float _S3021 = _S3016[int(4)];
    float _S3022 = _S3016[int(5)];
    float _S3023 = _S3016[int(6)];
    float _S3024 = _S3016[int(7)];
    float _S3025 = _S3016[int(8)];
    float _S3026 = _S3016[int(9)];
    float _S3027 = _S3016[int(10)];
    float _S3028 = _S3016[int(11)];
    float _S3029 = _S3016[int(12)];
    float _S3030 = _S3016[int(13)];
    float _S3031 = _S3016[int(14)];
    float _S3032 = _S3016[int(15)];
    float _S3033 = _S3016[int(16)];
    float _S3034 = _S3016[int(17)];
    float _S3035 = _S3016[int(18)];
    float _S3036 = _S3016[int(19)];
    float _S3037 = _S3016[int(20)];
    float _S3038 = _S3016[int(21)];
    float _S3039 = _S3016[int(22)];
    float _S3040;
    if(_S3008)
    {
        float _S3041 = 200.0f * _S3015;
        float _S3042 = _S3009 * _S3041 + 0.5f * (_S3007 * _S3041);
        _S3009 = 0.0f;
        _S3040 = _S3042;
    }
    else
    {
        _S3009 = _S3015;
        _S3040 = 0.0f;
    }
    DiffPair_float_0 _S3043;
    (&_S3043)->primal_0 = _S3007;
    (&_S3043)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3043, _S3009);
    float _S3044 = (_S3043.differential_0 + _S3040) / _S2979;
    FixedArray<float, 23>  _S3045;
    _S3045[int(0)] = 0.0f;
    _S3045[int(1)] = 0.0f;
    _S3045[int(2)] = 0.0f;
    _S3045[int(3)] = 0.0f;
    _S3045[int(4)] = 0.0f;
    _S3045[int(5)] = 0.0f;
    _S3045[int(6)] = 0.0f;
    _S3045[int(7)] = 0.0f;
    _S3045[int(8)] = 0.0f;
    _S3045[int(9)] = 0.0f;
    _S3045[int(10)] = 0.0f;
    _S3045[int(11)] = 0.0f;
    _S3045[int(12)] = 0.0f;
    _S3045[int(13)] = 0.0f;
    _S3045[int(14)] = 0.0f;
    _S3045[int(15)] = 0.0f;
    _S3045[int(16)] = 0.0f;
    _S3045[int(17)] = 0.0f;
    _S3045[int(18)] = 0.0f;
    _S3045[int(19)] = 0.0f;
    _S3045[int(20)] = 0.0f;
    _S3045[int(21)] = 0.0f;
    _S3045[int(22)] = 0.0f;
    _S3045[int(17)] = _S3044;
    float _S3046 = _S3017 + _S3045[int(0)];
    float _S3047 = _S3018 + _S3045[int(1)];
    float _S3048 = _S3019 + _S3045[int(2)];
    float _S3049 = _S3020 + _S3045[int(3)];
    float _S3050 = _S3021 + _S3045[int(4)];
    float _S3051 = _S3022 + _S3045[int(5)];
    float _S3052 = _S3023 + _S3045[int(6)];
    float _S3053 = _S3024 + _S3045[int(7)];
    float _S3054 = _S3025 + _S3045[int(8)];
    float _S3055 = _S3026 + _S3045[int(9)];
    float _S3056 = _S3027 + _S3045[int(10)];
    float _S3057 = _S3028 + _S3045[int(11)];
    float _S3058 = _S3029 + _S3045[int(12)];
    float _S3059 = _S3030 + _S3045[int(13)];
    float _S3060 = _S3031 + _S3045[int(14)];
    float _S3061 = _S3032 + _S3045[int(15)];
    float _S3062 = _S3033 + _S3045[int(16)];
    float _S3063 = _S3034 + _S3045[int(17)];
    float _S3064 = _S3035 + _S3045[int(18)];
    float _S3065 = _S3036 + _S3045[int(19)];
    float _S3066 = _S3037 + _S3045[int(20)];
    float _S3067 = _S3038 + _S3045[int(21)];
    float _S3068 = _S3039 + _S3045[int(22)];
    if(_S3005)
    {
        float _S3069 = 200.0f * _S3015;
        float _S3070 = _S3006 * _S3069 + 0.5f * (_S3004 * _S3069);
        _S3006 = 0.0f;
        _S3009 = _S3070;
    }
    else
    {
        _S3006 = _S3015;
        _S3009 = 0.0f;
    }
    DiffPair_float_0 _S3071;
    (&_S3071)->primal_0 = _S3004;
    (&_S3071)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3071, _S3006);
    float _S3072 = (_S3071.differential_0 + _S3009) / _S2979;
    FixedArray<float, 23>  _S3073;
    _S3073[int(0)] = 0.0f;
    _S3073[int(1)] = 0.0f;
    _S3073[int(2)] = 0.0f;
    _S3073[int(3)] = 0.0f;
    _S3073[int(4)] = 0.0f;
    _S3073[int(5)] = 0.0f;
    _S3073[int(6)] = 0.0f;
    _S3073[int(7)] = 0.0f;
    _S3073[int(8)] = 0.0f;
    _S3073[int(9)] = 0.0f;
    _S3073[int(10)] = 0.0f;
    _S3073[int(11)] = 0.0f;
    _S3073[int(12)] = 0.0f;
    _S3073[int(13)] = 0.0f;
    _S3073[int(14)] = 0.0f;
    _S3073[int(15)] = 0.0f;
    _S3073[int(16)] = 0.0f;
    _S3073[int(17)] = 0.0f;
    _S3073[int(18)] = 0.0f;
    _S3073[int(19)] = 0.0f;
    _S3073[int(20)] = 0.0f;
    _S3073[int(21)] = 0.0f;
    _S3073[int(22)] = 0.0f;
    _S3073[int(16)] = _S3072;
    float _S3074 = _S3046 + _S3073[int(0)];
    float _S3075 = _S3047 + _S3073[int(1)];
    float _S3076 = _S3048 + _S3073[int(2)];
    float _S3077 = _S3049 + _S3073[int(3)];
    float _S3078 = _S3050 + _S3073[int(4)];
    float _S3079 = _S3051 + _S3073[int(5)];
    float _S3080 = _S3052 + _S3073[int(6)];
    float _S3081 = _S3053 + _S3073[int(7)];
    float _S3082 = _S3054 + _S3073[int(8)];
    float _S3083 = _S3055 + _S3073[int(9)];
    float _S3084 = _S3056 + _S3073[int(10)];
    float _S3085 = _S3057 + _S3073[int(11)];
    float _S3086 = _S3058 + _S3073[int(12)];
    float _S3087 = _S3059 + _S3073[int(13)];
    float _S3088 = _S3060 + _S3073[int(14)];
    float _S3089 = _S3061 + _S3073[int(15)];
    float _S3090 = _S3062 + _S3073[int(16)];
    float _S3091 = _S3063 + _S3073[int(17)];
    float _S3092 = _S3064 + _S3073[int(18)];
    float _S3093 = _S3065 + _S3073[int(19)];
    float _S3094 = _S3066 + _S3073[int(20)];
    float _S3095 = _S3067 + _S3073[int(21)];
    float _S3096 = _S3068 + _S3073[int(22)];
    if(_S3002)
    {
        float _S3097 = 200.0f * _S3015;
        float _S3098 = _S3003 * _S3097 + 0.5f * (_S3001 * _S3097);
        _S3003 = 0.0f;
        _S3006 = _S3098;
    }
    else
    {
        _S3003 = _S3015;
        _S3006 = 0.0f;
    }
    DiffPair_float_0 _S3099;
    (&_S3099)->primal_0 = _S3001;
    (&_S3099)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3099, _S3003);
    float _S3100 = (_S3099.differential_0 + _S3006) / _S2979;
    FixedArray<float, 23>  _S3101;
    _S3101[int(0)] = 0.0f;
    _S3101[int(1)] = 0.0f;
    _S3101[int(2)] = 0.0f;
    _S3101[int(3)] = 0.0f;
    _S3101[int(4)] = 0.0f;
    _S3101[int(5)] = 0.0f;
    _S3101[int(6)] = 0.0f;
    _S3101[int(7)] = 0.0f;
    _S3101[int(8)] = 0.0f;
    _S3101[int(9)] = 0.0f;
    _S3101[int(10)] = 0.0f;
    _S3101[int(11)] = 0.0f;
    _S3101[int(12)] = 0.0f;
    _S3101[int(13)] = 0.0f;
    _S3101[int(14)] = 0.0f;
    _S3101[int(15)] = 0.0f;
    _S3101[int(16)] = 0.0f;
    _S3101[int(17)] = 0.0f;
    _S3101[int(18)] = 0.0f;
    _S3101[int(19)] = 0.0f;
    _S3101[int(20)] = 0.0f;
    _S3101[int(21)] = 0.0f;
    _S3101[int(22)] = 0.0f;
    _S3101[int(15)] = _S3100;
    float _S3102 = _S3074 + _S3101[int(0)];
    float _S3103 = _S3075 + _S3101[int(1)];
    float _S3104 = _S3076 + _S3101[int(2)];
    float _S3105 = _S3077 + _S3101[int(3)];
    float _S3106 = _S3078 + _S3101[int(4)];
    float _S3107 = _S3079 + _S3101[int(5)];
    float _S3108 = _S3080 + _S3101[int(6)];
    float _S3109 = _S3081 + _S3101[int(7)];
    float _S3110 = _S3082 + _S3101[int(8)];
    float _S3111 = _S3083 + _S3101[int(9)];
    float _S3112 = _S3084 + _S3101[int(10)];
    float _S3113 = _S3085 + _S3101[int(11)];
    float _S3114 = _S3086 + _S3101[int(12)];
    float _S3115 = _S3087 + _S3101[int(13)];
    float _S3116 = _S3088 + _S3101[int(14)];
    float _S3117 = _S3089 + _S3101[int(15)];
    float _S3118 = _S3090 + _S3101[int(16)];
    float _S3119 = _S3091 + _S3101[int(17)];
    float _S3120 = _S3092 + _S3101[int(18)];
    float _S3121 = _S3093 + _S3101[int(19)];
    float _S3122 = _S3094 + _S3101[int(20)];
    float _S3123 = _S3095 + _S3101[int(21)];
    float _S3124 = _S3096 + _S3101[int(22)];
    if(_S2999)
    {
        float _S3125 = 200.0f * _S3015;
        float _S3126 = _S3000 * _S3125 + 0.5f * (_S2998 * _S3125);
        _S3000 = 0.0f;
        _S3003 = _S3126;
    }
    else
    {
        _S3000 = _S3015;
        _S3003 = 0.0f;
    }
    DiffPair_float_0 _S3127;
    (&_S3127)->primal_0 = _S2998;
    (&_S3127)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3127, _S3000);
    float _S3128 = (_S3127.differential_0 + _S3003) / _S2979;
    FixedArray<float, 23>  _S3129;
    _S3129[int(0)] = 0.0f;
    _S3129[int(1)] = 0.0f;
    _S3129[int(2)] = 0.0f;
    _S3129[int(3)] = 0.0f;
    _S3129[int(4)] = 0.0f;
    _S3129[int(5)] = 0.0f;
    _S3129[int(6)] = 0.0f;
    _S3129[int(7)] = 0.0f;
    _S3129[int(8)] = 0.0f;
    _S3129[int(9)] = 0.0f;
    _S3129[int(10)] = 0.0f;
    _S3129[int(11)] = 0.0f;
    _S3129[int(12)] = 0.0f;
    _S3129[int(13)] = 0.0f;
    _S3129[int(14)] = 0.0f;
    _S3129[int(15)] = 0.0f;
    _S3129[int(16)] = 0.0f;
    _S3129[int(17)] = 0.0f;
    _S3129[int(18)] = 0.0f;
    _S3129[int(19)] = 0.0f;
    _S3129[int(20)] = 0.0f;
    _S3129[int(21)] = 0.0f;
    _S3129[int(22)] = 0.0f;
    _S3129[int(14)] = _S3128;
    float _S3130 = _S3102 + _S3129[int(0)];
    float _S3131 = _S3103 + _S3129[int(1)];
    float _S3132 = _S3104 + _S3129[int(2)];
    float _S3133 = _S3105 + _S3129[int(3)];
    float _S3134 = _S3106 + _S3129[int(4)];
    float _S3135 = _S3107 + _S3129[int(5)];
    float _S3136 = _S3108 + _S3129[int(6)];
    float _S3137 = _S3109 + _S3129[int(7)];
    float _S3138 = _S3110 + _S3129[int(8)];
    float _S3139 = _S3111 + _S3129[int(9)];
    float _S3140 = _S3112 + _S3129[int(10)];
    float _S3141 = _S3113 + _S3129[int(11)];
    float _S3142 = _S3114 + _S3129[int(12)];
    float _S3143 = _S3115 + _S3129[int(13)];
    float _S3144 = _S3116 + _S3129[int(14)];
    float _S3145 = _S3117 + _S3129[int(15)];
    float _S3146 = _S3118 + _S3129[int(16)];
    float _S3147 = _S3119 + _S3129[int(17)];
    float _S3148 = _S3120 + _S3129[int(18)];
    float _S3149 = _S3121 + _S3129[int(19)];
    float _S3150 = _S3122 + _S3129[int(20)];
    float _S3151 = _S3123 + _S3129[int(21)];
    float _S3152 = _S3124 + _S3129[int(22)];
    if(_S2996)
    {
        float _S3153 = 200.0f * _S3015;
        float _S3154 = _S2997 * _S3153 + 0.5f * (_S2995 * _S3153);
        _S2997 = 0.0f;
        _S3000 = _S3154;
    }
    else
    {
        _S2997 = _S3015;
        _S3000 = 0.0f;
    }
    DiffPair_float_0 _S3155;
    (&_S3155)->primal_0 = _S2995;
    (&_S3155)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3155, _S2997);
    float _S3156 = (_S3155.differential_0 + _S3000) / _S2979;
    FixedArray<float, 23>  _S3157;
    _S3157[int(0)] = 0.0f;
    _S3157[int(1)] = 0.0f;
    _S3157[int(2)] = 0.0f;
    _S3157[int(3)] = 0.0f;
    _S3157[int(4)] = 0.0f;
    _S3157[int(5)] = 0.0f;
    _S3157[int(6)] = 0.0f;
    _S3157[int(7)] = 0.0f;
    _S3157[int(8)] = 0.0f;
    _S3157[int(9)] = 0.0f;
    _S3157[int(10)] = 0.0f;
    _S3157[int(11)] = 0.0f;
    _S3157[int(12)] = 0.0f;
    _S3157[int(13)] = 0.0f;
    _S3157[int(14)] = 0.0f;
    _S3157[int(15)] = 0.0f;
    _S3157[int(16)] = 0.0f;
    _S3157[int(17)] = 0.0f;
    _S3157[int(18)] = 0.0f;
    _S3157[int(19)] = 0.0f;
    _S3157[int(20)] = 0.0f;
    _S3157[int(21)] = 0.0f;
    _S3157[int(22)] = 0.0f;
    _S3157[int(13)] = _S3156;
    float _S3158 = _S3130 + _S3157[int(0)];
    float _S3159 = _S3131 + _S3157[int(1)];
    float _S3160 = _S3132 + _S3157[int(2)];
    float _S3161 = _S3133 + _S3157[int(3)];
    float _S3162 = _S3134 + _S3157[int(4)];
    float _S3163 = _S3135 + _S3157[int(5)];
    float _S3164 = _S3136 + _S3157[int(6)];
    float _S3165 = _S3137 + _S3157[int(7)];
    float _S3166 = _S3138 + _S3157[int(8)];
    float _S3167 = _S3139 + _S3157[int(9)];
    float _S3168 = _S3140 + _S3157[int(10)];
    float _S3169 = _S3141 + _S3157[int(11)];
    float _S3170 = _S3142 + _S3157[int(12)];
    float _S3171 = _S3143 + _S3157[int(13)];
    float _S3172 = _S3144 + _S3157[int(14)];
    float _S3173 = _S3145 + _S3157[int(15)];
    float _S3174 = _S3146 + _S3157[int(16)];
    float _S3175 = _S3147 + _S3157[int(17)];
    float _S3176 = _S3148 + _S3157[int(18)];
    float _S3177 = _S3149 + _S3157[int(19)];
    float _S3178 = _S3150 + _S3157[int(20)];
    float _S3179 = _S3151 + _S3157[int(21)];
    float _S3180 = _S3152 + _S3157[int(22)];
    if(_S2993)
    {
        float _S3181 = 200.0f * _S3015;
        float _S3182 = _S2994 * _S3181 + 0.5f * (_S2992 * _S3181);
        _S2994 = 0.0f;
        _S2997 = _S3182;
    }
    else
    {
        _S2994 = _S3015;
        _S2997 = 0.0f;
    }
    DiffPair_float_0 _S3183;
    (&_S3183)->primal_0 = _S2992;
    (&_S3183)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3183, _S2994);
    float _S3184 = (_S3183.differential_0 + _S2997) / _S2979;
    FixedArray<float, 23>  _S3185;
    _S3185[int(0)] = 0.0f;
    _S3185[int(1)] = 0.0f;
    _S3185[int(2)] = 0.0f;
    _S3185[int(3)] = 0.0f;
    _S3185[int(4)] = 0.0f;
    _S3185[int(5)] = 0.0f;
    _S3185[int(6)] = 0.0f;
    _S3185[int(7)] = 0.0f;
    _S3185[int(8)] = 0.0f;
    _S3185[int(9)] = 0.0f;
    _S3185[int(10)] = 0.0f;
    _S3185[int(11)] = 0.0f;
    _S3185[int(12)] = 0.0f;
    _S3185[int(13)] = 0.0f;
    _S3185[int(14)] = 0.0f;
    _S3185[int(15)] = 0.0f;
    _S3185[int(16)] = 0.0f;
    _S3185[int(17)] = 0.0f;
    _S3185[int(18)] = 0.0f;
    _S3185[int(19)] = 0.0f;
    _S3185[int(20)] = 0.0f;
    _S3185[int(21)] = 0.0f;
    _S3185[int(22)] = 0.0f;
    _S3185[int(12)] = _S3184;
    float _S3186 = _S3158 + _S3185[int(0)];
    float _S3187 = _S3159 + _S3185[int(1)];
    float _S3188 = _S3160 + _S3185[int(2)];
    float _S3189 = _S3161 + _S3185[int(3)];
    float _S3190 = _S3162 + _S3185[int(4)];
    float _S3191 = _S3163 + _S3185[int(5)];
    float _S3192 = _S3164 + _S3185[int(6)];
    float _S3193 = _S3165 + _S3185[int(7)];
    float _S3194 = _S3166 + _S3185[int(8)];
    float _S3195 = _S3167 + _S3185[int(9)];
    float _S3196 = _S3168 + _S3185[int(10)];
    float _S3197 = _S3169 + _S3185[int(11)];
    float _S3198 = _S3170 + _S3185[int(12)];
    float _S3199 = _S3171 + _S3185[int(13)];
    float _S3200 = _S3172 + _S3185[int(14)];
    float _S3201 = _S3173 + _S3185[int(15)];
    float _S3202 = _S3174 + _S3185[int(16)];
    float _S3203 = _S3175 + _S3185[int(17)];
    float _S3204 = _S3176 + _S3185[int(18)];
    float _S3205 = _S3177 + _S3185[int(19)];
    float _S3206 = _S3178 + _S3185[int(20)];
    float _S3207 = _S3179 + _S3185[int(21)];
    float _S3208 = _S3180 + _S3185[int(22)];
    if(_S2990)
    {
        float _S3209 = 200.0f * _S3015;
        float _S3210 = _S2991 * _S3209 + 0.5f * (_S2989 * _S3209);
        _S2991 = 0.0f;
        _S2994 = _S3210;
    }
    else
    {
        _S2991 = _S3015;
        _S2994 = 0.0f;
    }
    DiffPair_float_0 _S3211;
    (&_S3211)->primal_0 = _S2989;
    (&_S3211)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3211, _S2991);
    float _S3212 = (_S3211.differential_0 + _S2994) / _S2979;
    FixedArray<float, 23>  _S3213;
    _S3213[int(0)] = 0.0f;
    _S3213[int(1)] = 0.0f;
    _S3213[int(2)] = 0.0f;
    _S3213[int(3)] = 0.0f;
    _S3213[int(4)] = 0.0f;
    _S3213[int(5)] = 0.0f;
    _S3213[int(6)] = 0.0f;
    _S3213[int(7)] = 0.0f;
    _S3213[int(8)] = 0.0f;
    _S3213[int(9)] = 0.0f;
    _S3213[int(10)] = 0.0f;
    _S3213[int(11)] = 0.0f;
    _S3213[int(12)] = 0.0f;
    _S3213[int(13)] = 0.0f;
    _S3213[int(14)] = 0.0f;
    _S3213[int(15)] = 0.0f;
    _S3213[int(16)] = 0.0f;
    _S3213[int(17)] = 0.0f;
    _S3213[int(18)] = 0.0f;
    _S3213[int(19)] = 0.0f;
    _S3213[int(20)] = 0.0f;
    _S3213[int(21)] = 0.0f;
    _S3213[int(22)] = 0.0f;
    _S3213[int(11)] = _S3212;
    float _S3214 = _S3186 + _S3213[int(0)];
    float _S3215 = _S3187 + _S3213[int(1)];
    float _S3216 = _S3188 + _S3213[int(2)];
    float _S3217 = _S3189 + _S3213[int(3)];
    float _S3218 = _S3190 + _S3213[int(4)];
    float _S3219 = _S3191 + _S3213[int(5)];
    float _S3220 = _S3192 + _S3213[int(6)];
    float _S3221 = _S3193 + _S3213[int(7)];
    float _S3222 = _S3194 + _S3213[int(8)];
    float _S3223 = _S3195 + _S3213[int(9)];
    float _S3224 = _S3196 + _S3213[int(10)];
    float _S3225 = _S3197 + _S3213[int(11)];
    float _S3226 = _S3198 + _S3213[int(12)];
    float _S3227 = _S3199 + _S3213[int(13)];
    float _S3228 = _S3200 + _S3213[int(14)];
    float _S3229 = _S3201 + _S3213[int(15)];
    float _S3230 = _S3202 + _S3213[int(16)];
    float _S3231 = _S3203 + _S3213[int(17)];
    float _S3232 = _S3204 + _S3213[int(18)];
    float _S3233 = _S3205 + _S3213[int(19)];
    float _S3234 = _S3206 + _S3213[int(20)];
    float _S3235 = _S3207 + _S3213[int(21)];
    float _S3236 = _S3208 + _S3213[int(22)];
    if(_S2987)
    {
        float _S3237 = 200.0f * _S3015;
        float _S3238 = _S2988 * _S3237 + 0.5f * (_S2986 * _S3237);
        _S2988 = 0.0f;
        _S2991 = _S3238;
    }
    else
    {
        _S2988 = _S3015;
        _S2991 = 0.0f;
    }
    DiffPair_float_0 _S3239;
    (&_S3239)->primal_0 = _S2986;
    (&_S3239)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3239, _S2988);
    float _S3240 = (_S3239.differential_0 + _S2991) / _S2979;
    float _S3241 = _S3010 / _S2985;
    float _S3242 = _S3011 / _S2984;
    float _S3243 = _S3012 / _S2983;
    FixedArray<float, 23>  _S3244;
    _S3244[int(0)] = 0.0f;
    _S3244[int(1)] = 0.0f;
    _S3244[int(2)] = 0.0f;
    _S3244[int(3)] = 0.0f;
    _S3244[int(4)] = 0.0f;
    _S3244[int(5)] = 0.0f;
    _S3244[int(6)] = 0.0f;
    _S3244[int(7)] = 0.0f;
    _S3244[int(8)] = 0.0f;
    _S3244[int(9)] = 0.0f;
    _S3244[int(10)] = 0.0f;
    _S3244[int(11)] = 0.0f;
    _S3244[int(12)] = 0.0f;
    _S3244[int(13)] = 0.0f;
    _S3244[int(14)] = 0.0f;
    _S3244[int(15)] = 0.0f;
    _S3244[int(16)] = 0.0f;
    _S3244[int(17)] = 0.0f;
    _S3244[int(18)] = 0.0f;
    _S3244[int(19)] = 0.0f;
    _S3244[int(20)] = 0.0f;
    _S3244[int(21)] = 0.0f;
    _S3244[int(22)] = 0.0f;
    _S3244[int(10)] = _S3240;
    _S3244[int(9)] = _S3241;
    _S3244[int(8)] = _S3241;
    _S3244[int(7)] = _S3241;
    _S3244[int(6)] = _S3241;
    _S3244[int(5)] = _S3241;
    _S3244[int(4)] = _S3242;
    _S3244[int(3)] = _S3242;
    _S3244[int(2)] = _S3242;
    _S3244[int(1)] = _S3243;
    float _S3245 = _S3214 + _S3244[int(0)];
    float _S3246 = _S3215 + _S3244[int(1)];
    float _S3247 = _S3216 + _S3244[int(2)];
    float _S3248 = _S3217 + _S3244[int(3)];
    float _S3249 = _S3218 + _S3244[int(4)];
    float _S3250 = _S3219 + _S3244[int(5)];
    float _S3251 = _S3220 + _S3244[int(6)];
    float _S3252 = _S3221 + _S3244[int(7)];
    float _S3253 = _S3222 + _S3244[int(8)];
    float _S3254 = _S3223 + _S3244[int(9)];
    float _S3255 = _S3224 + _S3244[int(10)];
    float _S3256 = _S3225 + _S3244[int(11)];
    float _S3257 = _S3226 + _S3244[int(12)];
    float _S3258 = _S3227 + _S3244[int(13)];
    float _S3259 = _S3228 + _S3244[int(14)];
    float _S3260 = _S3229 + _S3244[int(15)];
    float _S3261 = _S3230 + _S3244[int(16)];
    float _S3262 = _S3231 + _S3244[int(17)];
    float _S3263 = _S3232 + _S3244[int(18)];
    float _S3264 = _S3233 + _S3244[int(19)];
    float _S3265 = _S3234 + _S3244[int(20)];
    float _S3266 = _S3235 + _S3244[int(21)];
    float _S3267 = _S3236 + _S3244[int(22)];
    if(_S2981)
    {
        float _S3268 = 10.0f * _S3013;
        float _S3269 = _S2982 * _S3268 + 0.5f * (_S2980 * _S3268);
        _S2982 = 0.0f;
        _S2988 = _S3269;
    }
    else
    {
        _S2982 = _S3013;
        _S2988 = 0.0f;
    }
    DiffPair_float_0 _S3270;
    (&_S3270)->primal_0 = _S2980;
    (&_S3270)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3270, _S2982);
    float _S3271 = (_S3270.differential_0 + _S2988) / _S2979;
    FixedArray<float, 23>  _S3272;
    _S3272[int(0)] = 0.0f;
    _S3272[int(1)] = 0.0f;
    _S3272[int(2)] = 0.0f;
    _S3272[int(3)] = 0.0f;
    _S3272[int(4)] = 0.0f;
    _S3272[int(5)] = 0.0f;
    _S3272[int(6)] = 0.0f;
    _S3272[int(7)] = 0.0f;
    _S3272[int(8)] = 0.0f;
    _S3272[int(9)] = 0.0f;
    _S3272[int(10)] = 0.0f;
    _S3272[int(11)] = 0.0f;
    _S3272[int(12)] = 0.0f;
    _S3272[int(13)] = 0.0f;
    _S3272[int(14)] = 0.0f;
    _S3272[int(15)] = 0.0f;
    _S3272[int(16)] = 0.0f;
    _S3272[int(17)] = 0.0f;
    _S3272[int(18)] = 0.0f;
    _S3272[int(19)] = 0.0f;
    _S3272[int(20)] = 0.0f;
    _S3272[int(21)] = 0.0f;
    _S3272[int(22)] = 0.0f;
    _S3272[int(0)] = _S3271;
    FixedArray<float, 23>  _S3273 = {
        _S3245 + _S3272[int(0)], _S3246 + _S3272[int(1)], _S3247 + _S3272[int(2)], _S3248 + _S3272[int(3)], _S3249 + _S3272[int(4)], _S3250 + _S3272[int(5)], _S3251 + _S3272[int(6)], _S3252 + _S3272[int(7)], _S3253 + _S3272[int(8)], _S3254 + _S3272[int(9)], _S3255 + _S3272[int(10)], _S3256 + _S3272[int(11)], _S3257 + _S3272[int(12)], _S3258 + _S3272[int(13)], _S3259 + _S3272[int(14)], _S3260 + _S3272[int(15)], _S3261 + _S3272[int(16)], _S3262 + _S3272[int(17)], _S3263 + _S3272[int(18)], _S3264 + _S3272[int(19)], _S3265 + _S3272[int(20)], _S3266 + _S3272[int(21)], _S3267 + _S3272[int(22)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S3273;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S3274, int _S3275, FixedArray<float, 6>  * _S3276, FixedArray<float, 6>  * _S3277)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S3274, _S3275, _S3276, _S3277);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_5, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_5, FixedArray<float, 23>  * _S3278)
{
    FixedArray<float, 23>  _S3279 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S3279;
    FixedArray<float, 6>  _S3280 = loss_weights_5;
    FixedArray<float, 6>  _S3281 = grad_out_5;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_2, num_cameras_5, &_S3280, &_S3281);
    *_S3278 = (&dp_raw_losses_2)->differential_0;
    return;
}

