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
    float3  result_16;
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
        *_slang_vector_get_element_ptr(&result_16, i_5) = sum_8;
        i_5 = i_5 + int(1);
    }
    return result_16;
}

inline __device__ float2  mul_1(Matrix<float, 2, 2>  left_3, float2  right_3)
{
    float2  result_17;
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
        *_slang_vector_get_element_ptr(&result_17, i_6) = sum_10;
        i_6 = i_6 + int(1);
    }
    return result_17;
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
    Matrix<float, 3, 3>  result_18;
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
                float sum_13 = sum_12 + _slang_vector_get_element(left_5.rows[r_3], i_7) * _slang_vector_get_element(right_5.rows[i_7], c_2);
                i_7 = i_7 + int(1);
                sum_12 = sum_13;
            }
            *_slang_vector_get_element_ptr(((&result_18)->rows + (r_3)), c_2) = sum_12;
            c_2 = c_2 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_18;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_3(mul_3(R_1, covarW_0), transpose_0(R_1));
    return;
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
    Matrix<float, 3, 3>  M_0 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_3(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_4, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_16 = quat_4.y;
    float x2_2 = x_16 * x_16;
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

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_0, FixedArray<float, 10>  dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float _S283 = 0.0f * v_0;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float s_diff_r2_0 = u_0 + u_0 + (_S283 + _S283);
    float _S284 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
    float _S285 = dist_coeffs_0[int(1)] + r2_0 * _S284;
    float _S286 = dist_coeffs_0[int(0)] + r2_0 * _S285;
    float radial_0 = 1.0f + r2_0 * _S286;
    float _S287 = 2.0f * dist_coeffs_0[int(4)];
    float _S288 = _S287 * u_0;
    float _S289 = 2.0f * u_0;
    float _S290 = 2.0f * dist_coeffs_0[int(5)];
    float _S291 = _S290 * u_0;
    float _S292 = 2.0f * v_0;
    float2  _S293 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S286 + (s_diff_r2_0 * _S285 + (s_diff_r2_0 * _S284 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (_S287 * v_0 + 0.0f * _S288 + (s_diff_r2_0 + (_S289 + _S289)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S290 * v_0 + 0.0f * _S291 + (s_diff_r2_0 + (_S283 + 0.0f * _S292)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
    float _S294 = 0.0f * u_0;
    float s_diff_r2_1 = _S294 + _S294 + (v_0 + v_0);
    float2  _S295 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S286 + (s_diff_r2_1 * _S285 + (s_diff_r2_1 * _S284 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * uv_0 + make_float2 (0.0f * _S287 * v_0 + _S288 + (s_diff_r2_1 + (_S294 + 0.0f * _S289)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S290 * v_0 + _S291 + (s_diff_r2_1 + (_S292 + _S292)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
    Matrix<float, 2, 2>  _S296 = transpose_1(makeMatrix<float, 2, 2> (_S293 + make_float2 (_S293.x * dist_coeffs_0[int(8)] + _S293.y * dist_coeffs_0[int(9)], 0.0f), _S295 + make_float2 (_S295.x * dist_coeffs_0[int(8)] + _S295.y * dist_coeffs_0[int(9)], 0.0f)));
    return (F32_min((determinant_0(_S296)), ((F32_min((_S296.rows[int(0)].x), (_S296.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_1, bool is_fisheye_0, FixedArray<float, 10>  dist_coeffs_1)
{
    float2  _S297;
    if(is_fisheye_0)
    {
        float r_4 = length_1(uv_1);
        float theta_0 = (F32_atan((r_4)));
        float _S298;
        if(r_4 < 0.00100000004749745f)
        {
            _S298 = 1.0f - r_4 * r_4 / 3.0f;
        }
        else
        {
            _S298 = theta_0 / r_4;
        }
        _S297 = uv_1 * make_float2 (_S298);
    }
    else
    {
        _S297 = uv_1;
    }
    float u_1 = _S297.x;
    float v_1 = _S297.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S299 = _S297 * make_float2 (1.0f + r2_1 * (dist_coeffs_1[int(0)] + r2_1 * (dist_coeffs_1[int(1)] + r2_1 * (dist_coeffs_1[int(2)] + r2_1 * dist_coeffs_1[int(3)])))) + make_float2 (2.0f * dist_coeffs_1[int(4)] * u_1 * v_1 + dist_coeffs_1[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_1[int(6)] * r2_1, 2.0f * dist_coeffs_1[int(5)] * u_1 * v_1 + dist_coeffs_1[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_1[int(7)] * r2_1);
    return _S299 + make_float2 (dist_coeffs_1[int(8)] * _S299.x + dist_coeffs_1[int(9)] * _S299.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_2, FixedArray<float, 10>  * dist_coeffs_2, int maxiter_0, float2  * uv_undist_0)
{
    int i_8 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_8 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S300 = (*dist_coeffs_2)[int(3)];
        float _S301 = (*dist_coeffs_2)[int(4)];
        float _S302 = (*dist_coeffs_2)[int(5)];
        float _S303 = (*dist_coeffs_2)[int(6)];
        float _S304 = (*dist_coeffs_2)[int(7)];
        float _S305 = (*dist_coeffs_2)[int(8)];
        float _S306 = (*dist_coeffs_2)[int(9)];
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float _S307 = (*dist_coeffs_2)[int(2)] + r2_2 * (*dist_coeffs_2)[int(3)];
        float _S308 = (*dist_coeffs_2)[int(1)] + r2_2 * _S307;
        float _S309 = (*dist_coeffs_2)[int(0)] + r2_2 * _S308;
        float radial_1 = 1.0f + r2_2 * _S309;
        float _S310 = 2.0f * (*dist_coeffs_2)[int(4)];
        float _S311 = _S310 * u_2;
        float _S312 = 2.0f * u_2;
        float _S313 = 2.0f * (*dist_coeffs_2)[int(5)];
        float _S314 = _S313 * u_2;
        float _S315 = 2.0f * v_2;
        float2  _S316 = q_0 * make_float2 (radial_1) + make_float2 (_S311 * v_2 + (*dist_coeffs_2)[int(5)] * (r2_2 + _S312 * u_2) + (*dist_coeffs_2)[int(6)] * r2_2, _S314 * v_2 + (*dist_coeffs_2)[int(4)] * (r2_2 + _S315 * v_2) + (*dist_coeffs_2)[int(7)] * r2_2);
        float2  r_5 = _S316 + make_float2 ((*dist_coeffs_2)[int(8)] * _S316.x + (*dist_coeffs_2)[int(9)] * _S316.y, 0.0f) - uv_2;
        float _S317 = 0.0f * v_2;
        float s_diff_r2_2 = u_2 + u_2 + (_S317 + _S317);
        float2  _S318 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S309 + (s_diff_r2_2 * _S308 + (s_diff_r2_2 * _S307 + s_diff_r2_2 * _S300 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (_S310 * v_2 + 0.0f * _S311 + (s_diff_r2_2 + (_S312 + _S312)) * _S302 + s_diff_r2_2 * _S303, _S313 * v_2 + 0.0f * _S314 + (s_diff_r2_2 + (_S317 + 0.0f * _S315)) * _S301 + s_diff_r2_2 * _S304);
        float _S319 = 0.0f * u_2;
        float s_diff_r2_3 = _S319 + _S319 + (v_2 + v_2);
        float2  _S320 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S309 + (s_diff_r2_3 * _S308 + (s_diff_r2_3 * _S307 + s_diff_r2_3 * _S300 * r2_2) * r2_2) * r2_2) * q_0 + make_float2 (0.0f * _S310 * v_2 + _S311 + (s_diff_r2_3 + (_S319 + 0.0f * _S312)) * _S302 + s_diff_r2_3 * _S303, 0.0f * _S313 * v_2 + _S314 + (s_diff_r2_3 + (_S315 + _S315)) * _S301 + s_diff_r2_3 * _S304);
        Matrix<float, 2, 2>  _S321 = transpose_1(makeMatrix<float, 2, 2> (_S318 + make_float2 (_S318.x * _S305 + _S318.y * _S306, 0.0f), _S320 + make_float2 (_S320.x * _S305 + _S320.y * _S306, 0.0f)));
        float inv_det_0 = 1.0f / (_S321.rows[int(0)].x * _S321.rows[int(1)].y - _S321.rows[int(0)].y * _S321.rows[int(1)].x);
        float _S322 = r_5.x;
        float _S323 = r_5.y;
        float2  q_1 = q_0 - make_float2 ((_S322 * _S321.rows[int(1)].y - _S323 * _S321.rows[int(0)].y) * inv_det_0, (- _S322 * _S321.rows[int(1)].x + _S323 * _S321.rows[int(0)].x) * inv_det_0);
        i_8 = i_8 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S324 = (*dist_coeffs_2)[int(0)];
    float _S325 = (*dist_coeffs_2)[int(1)];
    float _S326 = (*dist_coeffs_2)[int(2)];
    float _S327 = (*dist_coeffs_2)[int(3)];
    float _S328 = (*dist_coeffs_2)[int(4)];
    float _S329 = (*dist_coeffs_2)[int(5)];
    float _S330 = (*dist_coeffs_2)[int(6)];
    float _S331 = (*dist_coeffs_2)[int(7)];
    float _S332 = (*dist_coeffs_2)[int(8)];
    float _S333 = (*dist_coeffs_2)[int(9)];
    float u_3 = q_0.x;
    float v_3 = q_0.y;
    float _S334 = 0.0f * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_4 = u_3 + u_3 + (_S334 + _S334);
    float _S335 = (*dist_coeffs_2)[int(2)] + r2_3 * (*dist_coeffs_2)[int(3)];
    float _S336 = (*dist_coeffs_2)[int(1)] + r2_3 * _S335;
    float _S337 = (*dist_coeffs_2)[int(0)] + r2_3 * _S336;
    float radial_2 = 1.0f + r2_3 * _S337;
    float _S338 = 2.0f * (*dist_coeffs_2)[int(4)];
    float _S339 = _S338 * u_3;
    float _S340 = 2.0f * u_3;
    float _S341 = 2.0f * (*dist_coeffs_2)[int(5)];
    float _S342 = _S341 * u_3;
    float _S343 = 2.0f * v_3;
    float2  _S344 = make_float2 (1.0f, 0.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_4 * _S337 + (s_diff_r2_4 * _S336 + (s_diff_r2_4 * _S335 + s_diff_r2_4 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (_S338 * v_3 + 0.0f * _S339 + (s_diff_r2_4 + (_S340 + _S340)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_4 * (*dist_coeffs_2)[int(6)], _S341 * v_3 + 0.0f * _S342 + (s_diff_r2_4 + (_S334 + 0.0f * _S343)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_4 * (*dist_coeffs_2)[int(7)]);
    float _S345 = 0.0f * u_3;
    float s_diff_r2_5 = _S345 + _S345 + (v_3 + v_3);
    float2  _S346 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_5 * _S337 + (s_diff_r2_5 * _S336 + (s_diff_r2_5 * _S335 + s_diff_r2_5 * (*dist_coeffs_2)[int(3)] * r2_3) * r2_3) * r2_3) * q_0 + make_float2 (0.0f * _S338 * v_3 + _S339 + (s_diff_r2_5 + (_S345 + 0.0f * _S340)) * (*dist_coeffs_2)[int(5)] + s_diff_r2_5 * (*dist_coeffs_2)[int(6)], 0.0f * _S341 * v_3 + _S342 + (s_diff_r2_5 + (_S343 + _S343)) * (*dist_coeffs_2)[int(4)] + s_diff_r2_5 * (*dist_coeffs_2)[int(7)]);
    Matrix<float, 2, 2>  _S347 = transpose_1(makeMatrix<float, 2, 2> (_S344 + make_float2 (_S344.x * (*dist_coeffs_2)[int(8)] + _S344.y * (*dist_coeffs_2)[int(9)], 0.0f), _S346 + make_float2 (_S346.x * (*dist_coeffs_2)[int(8)] + _S346.y * (*dist_coeffs_2)[int(9)], 0.0f)));
    bool _S348;
    if((F32_min((determinant_0(_S347)), ((F32_min((_S347.rows[int(0)].x), (_S347.rows[int(1)].y)))))) > 0.0f)
    {
        float u_4 = (*uv_undist_0).x;
        float v_4 = (*uv_undist_0).y;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float2  _S349 = *uv_undist_0 * make_float2 (1.0f + r2_4 * (_S324 + r2_4 * (_S325 + r2_4 * (_S326 + r2_4 * _S327)))) + make_float2 (_S338 * u_4 * v_4 + _S329 * (r2_4 + 2.0f * u_4 * u_4) + _S330 * r2_4, _S341 * u_4 * v_4 + _S328 * (r2_4 + 2.0f * v_4 * v_4) + _S331 * r2_4);
        _S348 = (length_1(_S349 + make_float2 (_S332 * _S349.x + _S333 * _S349.y, 0.0f) - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S348 = false;
    }
    return _S348;
}

inline __device__ bool undistort_point(float2  uv_3, bool is_fisheye_1, FixedArray<float, 10>  dist_coeffs_3, float2  * uv_undist_1)
{
    float2  _S350 = uv_3;
    FixedArray<float, 10>  _S351 = dist_coeffs_3;
    bool _S352 = undistort_point_0(uv_3, &_S351, int(8), &_S350);
    if(!_S352)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S353 = _S350;
        float theta_1 = length_1(_S350);
        float _S354;
        if(theta_1 < 0.00100000004749745f)
        {
            _S354 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S354 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S355 = make_float3 ((_S353 * make_float2 (_S354)).x, (_S353 * make_float2 (_S354)).y, (F32_cos((theta_1))));
        raydir_0 = _S355;
    }
    else
    {
        raydir_0 = make_float3 (_S350.x, _S350.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  dist_coeffs_4, float3  * raydir_1)
{
    float2  _S356 = uv_4;
    int3  _S357 = make_int3 (int(0));
    float3  _S358 = make_float3 ((float)_S357.x, (float)_S357.y, (float)_S357.z);
    *raydir_1 = _S358;
    FixedArray<float, 10>  _S359 = dist_coeffs_4;
    bool _S360 = undistort_point_0(uv_4, &_S359, int(8), &_S356);
    if(!_S360)
    {
        return false;
    }
    float3  _S361;
    if(is_fisheye_2)
    {
        float2  _S362 = _S356;
        float theta_2 = length_1(_S356);
        float _S363;
        if(theta_2 < 0.00100000004749745f)
        {
            _S363 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S363 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S364 = make_float3 ((_S362 * make_float2 (_S363)).x, (_S362 * make_float2 (_S363)).y, (F32_cos((theta_2))));
        _S361 = _S364;
    }
    else
    {
        _S361 = make_float3 (_S356.x, _S356.y, 1.0f);
    }
    *raydir_1 = _S361;
    return true;
}

inline __device__ float3  normalize_0(float3  x_17)
{
    return x_17 / make_float3 (length_2(x_17));
}

inline __device__ float4  normalize_1(float4  x_18)
{
    return x_18 / make_float4 (length_0(x_18));
}

inline __device__ bool generate_ray(float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  dist_coeffs_5, float3  * raydir_2)
{
    float2  _S365 = uv_5;
    FixedArray<float, 10>  _S366 = dist_coeffs_5;
    bool _S367 = undistort_point_0(uv_5, &_S366, int(8), &_S365);
    if(!_S367)
    {
        int3  _S368 = make_int3 (int(0));
        float3  _S369 = make_float3 ((float)_S368.x, (float)_S368.y, (float)_S368.z);
        *raydir_2 = _S369;
        return false;
    }
    float3  _S370;
    if(is_fisheye_3)
    {
        float2  _S371 = _S365;
        float theta_3 = length_1(_S365);
        float _S372;
        if(theta_3 < 0.00100000004749745f)
        {
            _S372 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S372 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S373 = make_float3 ((_S371 * make_float2 (_S372)).x, (_S371 * make_float2 (_S372)).y, (F32_cos((theta_3))));
        _S370 = _S373;
    }
    else
    {
        _S370 = make_float3 (_S365.x, _S365.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S370);
    return true;
}

inline __device__ void _d_mul_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_6, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_6, float3  dOut_11)
{
    float _S374 = (*right_6).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 3, 3>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = (*left_6).primal_0.x * dOut_11.x;
    float sum_14 = _S374 + (*right_6).primal_0.rows[int(0)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = (*left_6).primal_0.x * dOut_11.y;
    float sum_15 = sum_14 + (*right_6).primal_0.rows[int(0)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(0)))->z) = (*left_6).primal_0.x * dOut_11.z;
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = sum_15;
    float _S375 = (*right_6).primal_0.rows[int(1)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = (*left_6).primal_0.y * dOut_11.x;
    float sum_16 = _S375 + (*right_6).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = (*left_6).primal_0.y * dOut_11.y;
    float sum_17 = sum_16 + (*right_6).primal_0.rows[int(1)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(1)))->z) = (*left_6).primal_0.y * dOut_11.z;
    *&((&left_d_result_3)->y) = sum_17;
    float _S376 = (*right_6).primal_0.rows[int(2)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = (*left_6).primal_0.z * dOut_11.x;
    float sum_18 = _S376 + (*right_6).primal_0.rows[int(2)].y * dOut_11.y;
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
    float3  result_19;
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
        *_slang_vector_get_element_ptr(&result_19, j_2) = sum_20;
        j_2 = j_2 + int(1);
    }
    return result_19;
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

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S377, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S378, float3  _S379)
{
    _d_mul_2(_S377, _S378, _S379);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S380 = - _s_dOut_2;
    float3  _S381 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S382;
    (&_S382)->primal_0 = (*dpt_0).primal_0;
    (&_S382)->differential_0 = _S381;
    Matrix<float, 3, 3>  _S383 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S384;
    (&_S384)->primal_0 = (*dpR_0).primal_0;
    (&_S384)->differential_0 = _S383;
    s_bwd_prop_mul_0(&_S382, &_S384, _S380);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S382.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S384.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S385, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S386, float3  _S387)
{
    s_bwd_prop_transform_ray_o_0(_S385, _S386, _S387);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_5, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S388 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_5;
    (&dp_R_0)->differential_0 = _S388;
    float3  _S389 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S389;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S390 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S391;
    (&_S391)->primal_0 = (*dpraydir_0).primal_0;
    (&_S391)->differential_0 = _S390;
    Matrix<float, 3, 3>  _S392 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S393;
    (&_S393)->primal_0 = (*dpR_1).primal_0;
    (&_S393)->differential_0 = _S392;
    s_bwd_prop_mul_0(&_S391, &_S393, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S391.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S393.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S394, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S395, float3  _S396)
{
    s_bwd_prop_transform_ray_d_0(_S394, _S395, _S396);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_6, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S397 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_6;
    (&dp_R_1)->differential_0 = _S397;
    float3  _S398 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S398;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S399 = scale_2.x;
    float sx_0 = (F32_exp((_S399)));
    float _S400 = scale_2.y;
    float sy_0 = (F32_exp((_S400)));
    float sz_0 = scale_2.z - 0.5f * (_S399 + _S400);
    float x_19 = quat_5.y;
    float x2_3 = x_19 * x_19;
    float y2_3 = quat_5.z * quat_5.z;
    float z2_3 = quat_5.w * quat_5.w;
    float xy_3 = quat_5.y * quat_5.z;
    float xz_3 = quat_5.y * quat_5.w;
    float yz_3 = quat_5.z * quat_5.w;
    float wx_3 = quat_5.x * quat_5.y;
    float wy_3 = quat_5.x * quat_5.z;
    float wz_3 = quat_5.x * quat_5.w;
    Matrix<float, 3, 3>  _S401 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S401, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S401, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S401, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
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

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_1, float3  scale_3, float4  quat_6, float opac_0)
{
    float4  _S402 = normalize_1(quat_6);
    float3  _S403 = exp_0(scale_3);
    float x_21 = _S402.y;
    float x2_4 = x_21 * x_21;
    float y2_4 = _S402.z * _S402.z;
    float z2_4 = _S402.w * _S402.w;
    float xy_4 = _S402.y * _S402.z;
    float xz_4 = _S402.y * _S402.w;
    float yz_4 = _S402.z * _S402.w;
    float wx_4 = _S402.x * _S402.y;
    float wy_4 = _S402.x * _S402.z;
    float wz_4 = _S402.x * _S402.w;
    Matrix<float, 3, 3>  M_2 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S403.x, 0.0f, 0.0f, 0.0f, _S403.y, 0.0f, 0.0f, 0.0f, _S403.z));
    float4  _S404 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S405 = _S404 - floor_0(_S404);
    float4  _S406 = _S405 + make_float4 (dot_1(_S405, float4 {_S405.w, _S405.z, _S405.x, _S405.y} + make_float4 (33.3300018310546875f)));
    float4  _S407 = (float4 {_S406.x, _S406.x, _S406.y, _S406.z} + float4 {_S406.y, _S406.z, _S406.z, _S406.w}) * float4 {_S406.z, _S406.y, _S406.w, _S406.x};
    float4  _S408 = _S407 - floor_0(_S407);
    float2  _S409 = float2 {_S408.x, _S408.z};
    float _S410 = 6.28318548202514648f * _S409.y;
    float2  _S411 = float2 {_S408.y, _S408.w};
    float _S412 = 6.28318548202514648f * _S411.y;
    *mean_1 = *mean_1 + mul_0(mul_3(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S409.x))))))) * make_float2 ((F32_cos((_S410))), (F32_sin((_S410))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S409.x))))))) * make_float2 ((F32_cos((_S410))), (F32_sin((_S410))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S411.x))))))) * make_float2 ((F32_cos((_S412))), (F32_sin((_S412))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_12)
{
    float _S413 = dOut_12.y;
    float _S414 = dOut_12.z;
    float _S415 = dOut_12.x;
    float _S416 = (*a_0).primal_0.z * _S413 + - (*a_0).primal_0.y * _S414;
    float _S417 = - (*a_0).primal_0.z * _S415 + (*a_0).primal_0.x * _S414;
    float _S418 = (*a_0).primal_0.y * _S415 + - (*a_0).primal_0.x * _S413;
    float3  _S419 = make_float3 (- (*b_1).primal_0.z * _S413 + (*b_1).primal_0.y * _S414, (*b_1).primal_0.z * _S415 + - (*b_1).primal_0.x * _S414, - (*b_1).primal_0.y * _S415 + (*b_1).primal_0.x * _S413);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S419;
    float3  _S420 = make_float3 (_S416, _S417, _S418);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S420;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S421 = left_8.y;
    float _S422 = right_8.z;
    float _S423 = left_8.z;
    float _S424 = right_8.y;
    float _S425 = right_8.x;
    float _S426 = left_8.x;
    return make_float3 (_S421 * _S422 - _S423 * _S424, _S423 * _S425 - _S426 * _S422, _S426 * _S424 - _S421 * _S425);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S427 = normalize_1(quat_7);
    float _S428 = scale_4.x;
    float sx_1 = (F32_exp((_S428)));
    float _S429 = scale_4.y;
    float sy_1 = (F32_exp((_S429)));
    float sz_1 = scale_4.z - 0.5f * (_S428 + _S429);
    float x_22 = _S427.y;
    float x2_5 = x_22 * x_22;
    float y2_5 = _S427.z * _S427.z;
    float z2_5 = _S427.w * _S427.w;
    float xy_5 = _S427.y * _S427.z;
    float xz_5 = _S427.y * _S427.w;
    float yz_5 = _S427.z * _S427.w;
    float wx_5 = _S427.x * _S427.y;
    float wy_5 = _S427.x * _S427.z;
    float wz_5 = _S427.x * _S427.w;
    Matrix<float, 3, 3>  _S430 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S430, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S430, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S430, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S431 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S432 = _S431 - floor_0(_S431);
    float4  _S433 = _S432 + make_float4 (dot_1(_S432, float4 {_S432.w, _S432.z, _S432.x, _S432.y} + make_float4 (33.3300018310546875f)));
    float4  _S434 = (float4 {_S433.x, _S433.x, _S433.y, _S433.z} + float4 {_S433.y, _S433.z, _S433.z, _S433.w}) * float4 {_S433.z, _S433.y, _S433.w, _S433.x};
    float4  _S435 = _S434 - floor_0(_S434);
    float2  _S436 = float2 {_S435.x, _S435.z};
    float _S437 = 6.28318548202514648f * _S436.y;
    float2  _S438 = float2 {_S435.y, _S435.w};
    float _S439 = 6.28318548202514648f * _S438.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S436.x))))))) * make_float2 ((F32_cos((_S437))), (F32_sin((_S437))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S436.x))))))) * make_float2 ((F32_cos((_S437))), (F32_sin((_S437))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S438.x))))))) * make_float2 ((F32_cos((_S439))), (F32_sin((_S439))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_9, float dOut_13)
{
    float _S440 = _slang_select(((*dpx_9).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_9).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S440;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S441 = _slang_select(((*dpx_10).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_10).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S441;
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
    DiffPair_float_0 _S442 = *dpx_11;
    bool _S443;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S443 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S443 = false;
    }
    float _S444;
    if(_S443)
    {
        _S444 = dOut_15;
    }
    else
    {
        _S444 = 0.0f;
    }
    dpx_11->primal_0 = _S442.primal_0;
    dpx_11->differential_0 = _S444;
    DiffPair_float_0 _S445 = *dpMin_0;
    if((_S442.primal_0) < ((*dpMin_0).primal_0))
    {
        _S444 = dOut_15;
    }
    else
    {
        _S444 = 0.0f;
    }
    dpMin_0->primal_0 = _S445.primal_0;
    dpMin_0->differential_0 = _S444;
    DiffPair_float_0 _S446 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S444 = dOut_15;
    }
    else
    {
        _S444 = 0.0f;
    }
    dpMax_0->primal_0 = _S446.primal_0;
    dpMax_0->differential_0 = _S444;
    return;
}

inline __device__ float clamp_0(float x_24, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_24), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_12, float dOut_16)
{
    float _S447 = -0.5f / ((*dpx_12).primal_0 * (F32_sqrt(((*dpx_12).primal_0)))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S447;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_17)
{
    float _S448 = (1.0f - (*dps_0).primal_0) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S448;
    DiffPair_float_0 _S449 = *dpy_3;
    float _S450 = (*dps_0).primal_0 * dOut_17;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S450;
    float _S451 = (_S449.primal_0 - (*dpx_13).primal_0) * dOut_17;
    dps_0->primal_0 = _S449.primal_0;
    dps_0->differential_0 = _S451;
    return;
}

inline __device__ float lerp_0(float x_25, float y_7, float s_4)
{
    return x_25 + (y_7 - x_25) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  weights_0, FixedArray<float, 23>  * _S452)
{
    float3  _S453;
    bool _S454;
    bool _S455;
    FixedArray<float, 23>  losses_1;
    float _S456 = float(mask_0);
    float3  _S457 = ref_rgb_0 - render_rgb_0;
    float3  _S458 = abs_0(_S457);
    losses_1[int(0)] = weights_0[int(0)] * _S456 * ((_S458.x + _S458.y + _S458.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S456 * clamp_0(dot_0(_S457, _S457) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S459 = float(depth_mask_0 & mask_0);
    float _S460 = _S459 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S461 = _S459 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S460;
    losses_1[int(3)] = _S461;
    losses_1[int(4)] = _S460 * _S460;
    losses_1[int(5)] = _S461 * _S461;
    losses_1[int(6)] = _S460 * _S461;
    bool _S462 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S463 = norm2_0 == 0.0f;
        _S454 = _S463;
        if(_S463)
        {
            _S453 = make_float3 (0.0f);
            break;
        }
        _S453 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S464;
    bool _S465 = !_S454;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S466 = norm2_1 == 0.0f;
        _S455 = _S466;
        if(_S466)
        {
            _S464 = make_float3 (0.0f);
            break;
        }
        _S464 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S467;
    float3  _S468;
    bool _S469 = !_S455;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S468 = make_float3 (0.0f);
            _S467 = false;
            break;
        }
        _S468 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S467 = _S462;
        break;
    }
    float _S470 = float(_S465 & _S467);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S453, _S468);
    losses_1[int(7)] = weights_0[int(2)] * _S470 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S471 = float(_S469 & _S467);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S464, _S468);
    losses_1[int(8)] = weights_0[int(2)] * _S471 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S472 = float(_S465 & _S469);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S453, _S464);
    losses_1[int(11)] = weights_0[int(5)] * _S472 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S473 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S474 = float(alpha_mask_0);
    float _S475 = float(ref_alpha_0);
    float _S476 = (F32_max((_S473), (_S475)));
    losses_1[int(9)] = weights_0[int(3)] * _S474 * - lerp_0((F32_log(((F32_max((1.0f - _S476), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S476), (9.99999997475242708e-07f)))))), _S475);
    float _S477 = 1.0f - _S473;
    float _S478 = 1.0f - _S475;
    float _S479 = (F32_max((_S477), (_S478)));
    losses_1[int(10)] = weights_0[int(4)] * _S474 * - lerp_0((F32_log(((F32_max((1.0f - _S479), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S479), (9.99999997475242708e-07f)))))), _S478);
    losses_1[int(12)] = weights_0[int(6)] * 4.0f * _S473 * _S477;
    float _S480 = (F32_max((_S473), (9.999999960041972e-13f)));
    losses_1[int(13)] = weights_0[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S480;
    losses_1[int(14)] = weights_0[int(8)] * depth_dist_0 / _S480;
    losses_1[int(15)] = weights_0[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S480;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S456;
    losses_1[int(18)] = _S459;
    losses_1[int(19)] = _S470;
    losses_1[int(20)] = _S471;
    losses_1[int(21)] = _S472;
    losses_1[int(22)] = _S474;
    *_S452 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S481, float3  _S482)
{
    return dot_0(_S481, _S482);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S483)
{
    return (F32_rsqrt((_S483)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S484, float _S485, float _S486)
{
    return clamp_0(_S484, _S485, _S486);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S487, DiffPair_float_0 * _S488, DiffPair_float_0 * _S489, float _S490)
{
    _d_lerp_0(_S487, _S488, _S489, _S490);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S491, DiffPair_float_0 * _S492, DiffPair_float_0 * _S493, float _S494)
{
    _d_clamp_0(_S491, _S492, _S493, _S494);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S495, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S496, float _S497)
{
    _d_dot_0(_S495, _S496, _S497);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S498, float _S499)
{
    _d_rsqrt_0(_S498, _S499);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S500, float3  _S501)
{
    _d_abs_vector_0(_S500, _S501);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_float_0 _S502 = *dprender_depth_0;
    DiffPair_float_0 _S503 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S504 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S505 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S506 = *dpref_normal_0;
    DiffPair_float_0 _S507 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S508 = *dprgb_dist_0;
    DiffPair_float_0 _S509 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S510 = *dpnormal_dist_0;
    float3  _S511 = make_float3 (0.0f);
    float _S512 = float(mask_1);
    float _S513 = (*weights_1)[int(0)] * _S512;
    float3  _S514 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S515 = s_primal_ctx_dot_0(_S514, _S514) * 0.3333333432674408f;
    float _S516 = float(depth_mask_1 & mask_1);
    float _S517 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S518 = _S516 * s_primal_ctx_log_0(_S517);
    float _S519 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S520 = _S516 * s_primal_ctx_log_0(_S519);
    bool _S521 = normal_mask_1 & mask_1;
    float _S522 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S523 = _S522 == 0.0f;
    float3  _S524;
    if(_S523)
    {
        _S524 = make_float3 (0.0f);
    }
    bool _S525 = !_S523;
    float3  _S526;
    if(_S525)
    {
        float _S527 = s_primal_ctx_rsqrt_0(_S522);
        float3  _S528 = make_float3 (_S527);
        _S524 = _S504.primal_0 * make_float3 (_S527);
        _S526 = _S528;
    }
    else
    {
        _S526 = _S511;
    }
    float _S529 = s_primal_ctx_dot_0(_S505.primal_0, _S505.primal_0);
    bool _S530 = _S529 == 0.0f;
    float3  _S531;
    if(_S530)
    {
        _S531 = make_float3 (0.0f);
    }
    bool _S532 = !_S530;
    float3  _S533;
    if(_S532)
    {
        float _S534 = s_primal_ctx_rsqrt_0(_S529);
        float3  _S535 = make_float3 (_S534);
        _S531 = _S505.primal_0 * make_float3 (_S534);
        _S533 = _S535;
    }
    else
    {
        _S533 = _S511;
    }
    float _S536 = s_primal_ctx_dot_0(_S506.primal_0, _S506.primal_0);
    bool _S537 = _S536 == 0.0f;
    float3  _S538;
    bool _S539;
    if(_S537)
    {
        float3  _S540 = make_float3 (0.0f);
        _S539 = false;
        _S538 = _S540;
    }
    else
    {
        _S539 = _S521;
    }
    bool _S541 = !_S537;
    float3  _S542;
    if(_S541)
    {
        float _S543 = s_primal_ctx_rsqrt_0(_S536);
        float3  _S544 = make_float3 (_S543);
        _S538 = _S506.primal_0 * make_float3 (_S543);
        _S542 = _S544;
    }
    else
    {
        _S542 = _S511;
    }
    float _S545 = (*weights_1)[int(2)] * float(_S525 & _S539);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S524, _S538);
    float _S546 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S547 = (*weights_1)[int(2)] * float(_S532 & _S539);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S531, _S538);
    float _S548 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S549 = (*weights_1)[int(5)] * float(_S525 & _S532);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S524, _S531);
    float _S550 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S551 = s_primal_ctx_clamp_0(_S507.primal_0, 0.0f, 1.0f);
    float _S552 = float(alpha_mask_1);
    float _S553 = (*weights_1)[int(3)] * _S552;
    float _S554 = float(ref_alpha_1);
    float _S555 = (F32_max((_S551), (_S554)));
    float _S556 = 1.0f - _S555;
    float _S557 = (F32_max((_S556), (9.99999997475242708e-07f)));
    float _S558 = s_primal_ctx_log_0(_S557);
    float _S559 = (F32_max((_S555), (9.99999997475242708e-07f)));
    float _S560 = s_primal_ctx_log_0(_S559);
    float _S561 = (*weights_1)[int(4)] * _S552;
    float _S562 = 1.0f - _S551;
    float _S563 = 1.0f - _S554;
    float _S564 = (F32_max((_S562), (_S563)));
    float _S565 = 1.0f - _S564;
    float _S566 = (F32_max((_S565), (9.99999997475242708e-07f)));
    float _S567 = s_primal_ctx_log_0(_S566);
    float _S568 = (F32_max((_S564), (9.99999997475242708e-07f)));
    float _S569 = s_primal_ctx_log_0(_S568);
    float _S570 = (*weights_1)[int(6)] * 4.0f;
    float _S571 = _S570 * _S551;
    float _S572 = (F32_max((_S551), (9.999999960041972e-13f)));
    float _S573 = _S572 * _S572;
    float _S574 = (*_s_dOut_4)[int(0)];
    float _S575 = (*_s_dOut_4)[int(1)];
    float _S576 = (*_s_dOut_4)[int(2)];
    float _S577 = (*_s_dOut_4)[int(3)];
    float _S578 = (*_s_dOut_4)[int(4)];
    float _S579 = (*_s_dOut_4)[int(5)];
    float _S580 = (*_s_dOut_4)[int(6)];
    float _S581 = (*_s_dOut_4)[int(15)] / _S573;
    float _S582 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S572 * _S581));
    float _S583 = (*_s_dOut_4)[int(14)] / _S573;
    float _S584 = (*weights_1)[int(8)] * (_S572 * _S583);
    float _S585 = (*_s_dOut_4)[int(13)] / _S573;
    float _S586 = _S572 * _S585;
    float _S587 = (*weights_1)[int(9)] * ((_S510.primal_0.x + _S510.primal_0.y + _S510.primal_0.z) * 0.3333333432674408f) * - _S581 + (*weights_1)[int(8)] * _S509.primal_0 * - _S583 + (*weights_1)[int(7)] * ((_S508.primal_0.x + _S508.primal_0.y + _S508.primal_0.z) * 0.3333333432674408f) * - _S585;
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = _S551;
    (&_S588)->differential_0 = 0.0f;
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = 9.999999960041972e-13f;
    (&_S589)->differential_0 = 0.0f;
    _d_max_0(&_S588, &_S589, _S587);
    float _S590 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S586);
    float _S591 = _S571 * (*_s_dOut_4)[int(12)];
    float _S592 = _S570 * (_S562 * (*_s_dOut_4)[int(12)]);
    float _S593 = - (_S561 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S594;
    (&_S594)->primal_0 = _S567;
    (&_S594)->differential_0 = 0.0f;
    DiffPair_float_0 _S595;
    (&_S595)->primal_0 = _S569;
    (&_S595)->differential_0 = 0.0f;
    DiffPair_float_0 _S596;
    (&_S596)->primal_0 = _S563;
    (&_S596)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S594, &_S595, &_S596, _S593);
    DiffPair_float_0 _S597;
    (&_S597)->primal_0 = _S568;
    (&_S597)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S597, _S595.differential_0);
    DiffPair_float_0 _S598;
    (&_S598)->primal_0 = _S564;
    (&_S598)->differential_0 = 0.0f;
    DiffPair_float_0 _S599;
    (&_S599)->primal_0 = 9.99999997475242708e-07f;
    (&_S599)->differential_0 = 0.0f;
    _d_max_0(&_S598, &_S599, _S597.differential_0);
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = _S566;
    (&_S600)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S600, _S594.differential_0);
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = _S565;
    (&_S601)->differential_0 = 0.0f;
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = 9.99999997475242708e-07f;
    (&_S602)->differential_0 = 0.0f;
    _d_max_0(&_S601, &_S602, _S600.differential_0);
    float _S603 = _S598.differential_0 + - _S601.differential_0;
    DiffPair_float_0 _S604;
    (&_S604)->primal_0 = _S562;
    (&_S604)->differential_0 = 0.0f;
    DiffPair_float_0 _S605;
    (&_S605)->primal_0 = _S563;
    (&_S605)->differential_0 = 0.0f;
    _d_max_0(&_S604, &_S605, _S603);
    float _S606 = - (_S591 + _S604.differential_0);
    float _S607 = - (_S553 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S608;
    (&_S608)->primal_0 = _S558;
    (&_S608)->differential_0 = 0.0f;
    DiffPair_float_0 _S609;
    (&_S609)->primal_0 = _S560;
    (&_S609)->differential_0 = 0.0f;
    DiffPair_float_0 _S610;
    (&_S610)->primal_0 = _S554;
    (&_S610)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S608, &_S609, &_S610, _S607);
    DiffPair_float_0 _S611;
    (&_S611)->primal_0 = _S559;
    (&_S611)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S611, _S609.differential_0);
    DiffPair_float_0 _S612;
    (&_S612)->primal_0 = _S555;
    (&_S612)->differential_0 = 0.0f;
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = 9.99999997475242708e-07f;
    (&_S613)->differential_0 = 0.0f;
    _d_max_0(&_S612, &_S613, _S611.differential_0);
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S557;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S614, _S608.differential_0);
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = _S556;
    (&_S615)->differential_0 = 0.0f;
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = 9.99999997475242708e-07f;
    (&_S616)->differential_0 = 0.0f;
    _d_max_0(&_S615, &_S616, _S614.differential_0);
    float _S617 = _S612.differential_0 + - _S615.differential_0;
    DiffPair_float_0 _S618;
    (&_S618)->primal_0 = _S551;
    (&_S618)->differential_0 = 0.0f;
    DiffPair_float_0 _S619;
    (&_S619)->primal_0 = _S554;
    (&_S619)->differential_0 = 0.0f;
    _d_max_0(&_S618, &_S619, _S617);
    float _S620 = _S588.differential_0 + _S592 + _S606 + _S618.differential_0;
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S507.primal_0;
    (&_S621)->differential_0 = 0.0f;
    DiffPair_float_0 _S622;
    (&_S622)->primal_0 = 0.0f;
    (&_S622)->differential_0 = 0.0f;
    DiffPair_float_0 _S623;
    (&_S623)->primal_0 = 1.0f;
    (&_S623)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S621, &_S622, &_S623, _S620);
    DiffPair_float_0 _S624 = _S621;
    float _S625 = _S549 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S626;
    (&_S626)->primal_0 = _S550;
    (&_S626)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S626, _S625);
    DiffPair_float_0 _S627;
    (&_S627)->primal_0 = cos_sim_loss_5;
    (&_S627)->differential_0 = 0.0f;
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = 9.999999960041972e-13f;
    (&_S628)->differential_0 = 0.0f;
    _d_max_0(&_S627, &_S628, _S626.differential_0);
    float _S629 = 0.5f * - (_S625 + _S627.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S630;
    (&_S630)->primal_0 = _S524;
    (&_S630)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S631;
    (&_S631)->primal_0 = _S531;
    (&_S631)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S630, &_S631, _S629);
    float _S632 = _S547 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S633;
    (&_S633)->primal_0 = _S548;
    (&_S633)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S633, _S632);
    DiffPair_float_0 _S634;
    (&_S634)->primal_0 = cos_sim_loss_4;
    (&_S634)->differential_0 = 0.0f;
    DiffPair_float_0 _S635;
    (&_S635)->primal_0 = 9.999999960041972e-13f;
    (&_S635)->differential_0 = 0.0f;
    _d_max_0(&_S634, &_S635, _S633.differential_0);
    float _S636 = 0.5f * - (_S632 + _S634.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S637;
    (&_S637)->primal_0 = _S531;
    (&_S637)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S538;
    (&_S638)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S637, &_S638, _S636);
    float _S639 = _S545 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S640;
    (&_S640)->primal_0 = _S546;
    (&_S640)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S640, _S639);
    DiffPair_float_0 _S641;
    (&_S641)->primal_0 = cos_sim_loss_3;
    (&_S641)->differential_0 = 0.0f;
    DiffPair_float_0 _S642;
    (&_S642)->primal_0 = 9.999999960041972e-13f;
    (&_S642)->differential_0 = 0.0f;
    _d_max_0(&_S641, &_S642, _S640.differential_0);
    float _S643 = 0.5f * - (_S639 + _S641.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S644;
    (&_S644)->primal_0 = _S524;
    (&_S644)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S645;
    (&_S645)->primal_0 = _S538;
    (&_S645)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S644, &_S645, _S643);
    float3  _S646 = _S638.differential_0 + _S645.differential_0;
    float3  _S647 = _S630.differential_0 + _S644.differential_0;
    float3  _S648 = make_float3 (_S582, _S582, _S582);
    float3  _S649 = make_float3 (_S590, _S590, _S590);
    float3  _S650 = _S631.differential_0 + _S637.differential_0;
    float _S651;
    if(_S541)
    {
        float3  _S652 = _S506.primal_0 * _S646;
        float3  _S653 = _S542 * _S646;
        float _S654 = _S652.x + _S652.y + _S652.z;
        DiffPair_float_0 _S655;
        (&_S655)->primal_0 = _S536;
        (&_S655)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S655, _S654);
        _S651 = _S655.differential_0;
        _S524 = _S653;
    }
    else
    {
        _S651 = 0.0f;
        _S524 = _S511;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S656;
    (&_S656)->primal_0 = _S506.primal_0;
    (&_S656)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S657;
    (&_S657)->primal_0 = _S506.primal_0;
    (&_S657)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S656, &_S657, _S651);
    float3  _S658 = _S657.differential_0 + _S656.differential_0 + _S524;
    if(_S532)
    {
        float3  _S659 = _S505.primal_0 * _S650;
        float3  _S660 = _S533 * _S650;
        float _S661 = _S659.x + _S659.y + _S659.z;
        DiffPair_float_0 _S662;
        (&_S662)->primal_0 = _S529;
        (&_S662)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S662, _S661);
        _S651 = _S662.differential_0;
        _S524 = _S660;
    }
    else
    {
        _S651 = 0.0f;
        _S524 = _S511;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S663;
    (&_S663)->primal_0 = _S505.primal_0;
    (&_S663)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S664;
    (&_S664)->primal_0 = _S505.primal_0;
    (&_S664)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S663, &_S664, _S651);
    float3  _S665 = _S664.differential_0 + _S663.differential_0 + _S524;
    if(_S525)
    {
        float3  _S666 = _S504.primal_0 * _S647;
        float3  _S667 = _S526 * _S647;
        float _S668 = _S666.x + _S666.y + _S666.z;
        DiffPair_float_0 _S669;
        (&_S669)->primal_0 = _S522;
        (&_S669)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S669, _S668);
        _S651 = _S669.differential_0;
        _S524 = _S667;
    }
    else
    {
        _S651 = 0.0f;
        _S524 = _S511;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S670;
    (&_S670)->primal_0 = _S504.primal_0;
    (&_S670)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S671;
    (&_S671)->primal_0 = _S504.primal_0;
    (&_S671)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S670, &_S671, _S651);
    float _S672 = _S520 * _S580;
    float _S673 = _S520 * _S579;
    float _S674 = _S518 * _S578;
    float _S675 = _S516 * (_S518 * _S580 + _S673 + _S673 + _S577);
    DiffPair_float_0 _S676;
    (&_S676)->primal_0 = _S519;
    (&_S676)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S676, _S675);
    DiffPair_float_0 _S677;
    (&_S677)->primal_0 = _S503.primal_0;
    (&_S677)->differential_0 = 0.0f;
    DiffPair_float_0 _S678;
    (&_S678)->primal_0 = 0.00009999999747379f;
    (&_S678)->differential_0 = 0.0f;
    _d_max_0(&_S677, &_S678, _S676.differential_0);
    float _S679 = _S516 * (_S672 + _S674 + _S674 + _S576);
    DiffPair_float_0 _S680;
    (&_S680)->primal_0 = _S517;
    (&_S680)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S680, _S679);
    DiffPair_float_0 _S681;
    (&_S681)->primal_0 = _S502.primal_0;
    (&_S681)->differential_0 = 0.0f;
    DiffPair_float_0 _S682;
    (&_S682)->primal_0 = 0.00009999999747379f;
    (&_S682)->differential_0 = 0.0f;
    _d_max_0(&_S681, &_S682, _S680.differential_0);
    float _S683 = _S512 * _S575;
    DiffPair_float_0 _S684;
    (&_S684)->primal_0 = _S515;
    (&_S684)->differential_0 = 0.0f;
    DiffPair_float_0 _S685;
    (&_S685)->primal_0 = 0.0f;
    (&_S685)->differential_0 = 0.0f;
    DiffPair_float_0 _S686;
    (&_S686)->primal_0 = 1.0f;
    (&_S686)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S684, &_S685, &_S686, _S683);
    float _S687 = 0.3333333432674408f * _S684.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S688;
    (&_S688)->primal_0 = _S514;
    (&_S688)->differential_0 = _S511;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S689;
    (&_S689)->primal_0 = _S514;
    (&_S689)->differential_0 = _S511;
    s_bwd_prop_dot_0(&_S688, &_S689, _S687);
    float _S690 = 0.3333333432674408f * (_S513 * _S574);
    float3  _S691 = make_float3 (_S690, _S690, _S690);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S692;
    (&_S692)->primal_0 = _S514;
    (&_S692)->differential_0 = _S511;
    s_bwd_prop_abs_0(&_S692, _S691);
    float3  _S693 = _S689.differential_0 + _S688.differential_0 + _S692.differential_0;
    float3  _S694 = - _S693;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S648;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S584;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S649;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S624.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S658;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S665;
    float3  _S695 = _S671.differential_0 + _S670.differential_0 + _S524;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S695;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S677.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S681.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S693;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S694;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S696, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S697, DiffPair_float_0 * _S698, DiffPair_float_0 * _S699, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S700, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S701, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S702, DiffPair_float_0 * _S703, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S704, DiffPair_float_0 * _S705, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S706, bool _S707, bool _S708, bool _S709, bool _S710, bool _S711, FixedArray<float, 10>  * _S712, FixedArray<float, 23>  * _S713)
{
    s_bwd_prop_per_pixel_losses_0(_S696, _S697, _S698, _S699, _S700, _S701, _S702, _S703, _S704, _S705, _S706, _S707, _S708, _S709, _S710, _S711, _S712, _S713);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S714 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S714;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S714;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S714;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S714;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S714;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S714;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S714;
    FixedArray<float, 10>  _S715 = weights_2;
    FixedArray<float, 23>  _S716 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S715, &_S716);
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
    float _S717 = 1.0f / ((*dpx_14).primal_0 * 2.30258512496948242f) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S717;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 10>  weights_3, FixedArray<float, 10>  * _S718)
{
    FixedArray<float, 10>  losses_2;
    float _S719 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_2[int(0)] = raw_losses_0[int(0)] / _S719;
    losses_2[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S719)));
    bool _S720;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S720 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S720 = false;
    }
    float _S721;
    if(_S720)
    {
        _S721 = weights_3[int(1)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S721 = 0.0f;
    }
    losses_2[int(2)] = _S721;
    losses_2[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_2[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S722 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_2[int(6)] = raw_losses_0[int(12)] / _S722;
    losses_2[int(7)] = raw_losses_0[int(13)] / _S722;
    losses_2[int(8)] = raw_losses_0[int(14)] / _S722;
    losses_2[int(9)] = raw_losses_0[int(15)] / _S722;
    *_S718 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S723)
{
    return (F32_sqrt((_S723)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S724, float _S725)
{
    _d_log10_0(_S724, _S725);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S726 = dpraw_losses_0->primal_0;
    float _S727 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S728 = _S727 * _S727;
    float _S729 = dpraw_losses_0->primal_0[int(1)] / _S727;
    bool _S730 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S731;
    if(_S730)
    {
        _S731 = (_S726[int(3)]) != 0.0f;
    }
    else
    {
        _S731 = false;
    }
    float _S732;
    float _S733;
    float _S734;
    float _S735;
    float _S736;
    float _S737;
    float _S738;
    float _S739;
    float _S740;
    float _S741;
    float _S742;
    float _S743;
    float _S744;
    float _S745;
    float _S746;
    if(_S731)
    {
        float _S747 = _S726[int(2)] * _S726[int(3)];
        float _S748 = _S726[int(18)] * _S726[int(18)];
        float _S749 = _S726[int(6)] - _S747 / _S726[int(18)];
        float _S750 = _S726[int(2)] * _S726[int(2)];
        float _S751 = _S726[int(4)] - _S750 / _S726[int(18)];
        float _S752 = _S726[int(3)] * _S726[int(3)];
        float _S753 = _S726[int(5)] - _S752 / _S726[int(18)];
        float _S754 = _S751 * _S753 + 1.0f;
        float _S755 = (F32_max((9.999999960041972e-13f), (_S754)));
        float _S756 = s_primal_ctx_sqrt_0(_S755);
        float _S757 = _S756 * _S756;
        float _S758 = 1.0f - _S749 / _S756;
        _S732 = (*weights_4)[int(1)];
        _S733 = _S758;
        _S734 = _S757;
        _S735 = _S749;
        _S736 = _S756;
        _S737 = _S755;
        _S738 = _S754;
        _S739 = _S751;
        _S740 = _S753;
        _S741 = _S748;
        _S742 = _S752;
        _S743 = _S726[int(3)];
        _S744 = _S750;
        _S745 = _S726[int(2)];
        _S746 = _S747;
    }
    else
    {
        _S732 = 0.0f;
        _S733 = 0.0f;
        _S734 = 0.0f;
        _S735 = 0.0f;
        _S736 = 0.0f;
        _S737 = 0.0f;
        _S738 = 0.0f;
        _S739 = 0.0f;
        _S740 = 0.0f;
        _S741 = 0.0f;
        _S742 = 0.0f;
        _S743 = 0.0f;
        _S744 = 0.0f;
        _S745 = 0.0f;
        _S746 = 0.0f;
    }
    float _S759 = (F32_max((_S726[int(19)]), (1.0f)));
    float _S760 = _S759 * _S759;
    float _S761 = (F32_max((_S726[int(20)]), (1.0f)));
    float _S762 = _S761 * _S761;
    float _S763 = float((I32_max((int((_S726[int(19)]) > 0.5f) + int((_S726[int(20)]) > 0.5f)), (int(1)))));
    float _S764 = _S726[int(9)] + _S726[int(10)];
    float _S765 = (F32_max((_S726[int(22)]), (1.0f)));
    float _S766 = _S765 * _S765;
    float _S767 = (F32_max((_S726[int(21)]), (1.0f)));
    float _S768 = _S767 * _S767;
    float _S769 = (F32_max((_S726[int(16)]), (1.0f)));
    float _S770 = _S769 * _S769;
    float _S771 = (*_s_dOut_5)[int(0)];
    float _S772 = (*_s_dOut_5)[int(1)];
    float _S773 = (*_s_dOut_5)[int(2)];
    float _S774 = (*_s_dOut_5)[int(9)] / _S770;
    float _S775 = _S769 * _S774;
    float _S776 = (*_s_dOut_5)[int(8)] / _S770;
    float _S777 = _S769 * _S776;
    float _S778 = (*_s_dOut_5)[int(7)] / _S770;
    float _S779 = _S769 * _S778;
    float _S780 = (*_s_dOut_5)[int(6)] / _S770;
    float _S781 = _S769 * _S780;
    float _S782 = _S726[int(15)] * - _S774 + _S726[int(14)] * - _S776 + _S726[int(13)] * - _S778 + _S726[int(12)] * - _S780;
    DiffPair_float_0 _S783;
    (&_S783)->primal_0 = _S726[int(16)];
    (&_S783)->differential_0 = 0.0f;
    DiffPair_float_0 _S784;
    (&_S784)->primal_0 = 1.0f;
    (&_S784)->differential_0 = 0.0f;
    _d_max_0(&_S783, &_S784, _S782);
    float _S785 = (*_s_dOut_5)[int(5)] / _S768;
    float _S786 = _S726[int(11)] * - _S785;
    float _S787 = _S767 * _S785;
    DiffPair_float_0 _S788;
    (&_S788)->primal_0 = _S726[int(21)];
    (&_S788)->differential_0 = 0.0f;
    DiffPair_float_0 _S789;
    (&_S789)->primal_0 = 1.0f;
    (&_S789)->differential_0 = 0.0f;
    _d_max_0(&_S788, &_S789, _S786);
    float _S790 = (*_s_dOut_5)[int(4)] / _S766;
    float _S791 = _S764 * - _S790;
    float _S792 = _S765 * _S790;
    DiffPair_float_0 _S793;
    (&_S793)->primal_0 = _S726[int(22)];
    (&_S793)->differential_0 = 0.0f;
    DiffPair_float_0 _S794;
    (&_S794)->primal_0 = 1.0f;
    (&_S794)->differential_0 = 0.0f;
    _d_max_0(&_S793, &_S794, _S791);
    float _S795 = (*_s_dOut_5)[int(3)] / _S763;
    float _S796 = _S795 / _S762;
    float _S797 = _S726[int(8)] * - _S796;
    float _S798 = _S761 * _S796;
    DiffPair_float_0 _S799;
    (&_S799)->primal_0 = _S726[int(20)];
    (&_S799)->differential_0 = 0.0f;
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = 1.0f;
    (&_S800)->differential_0 = 0.0f;
    _d_max_0(&_S799, &_S800, _S797);
    float _S801 = _S795 / _S760;
    float _S802 = _S726[int(7)] * - _S801;
    float _S803 = _S759 * _S801;
    DiffPair_float_0 _S804;
    (&_S804)->primal_0 = _S726[int(19)];
    (&_S804)->differential_0 = 0.0f;
    DiffPair_float_0 _S805;
    (&_S805)->primal_0 = 1.0f;
    (&_S805)->differential_0 = 0.0f;
    _d_max_0(&_S804, &_S805, _S802);
    FixedArray<float, 23>  _S806;
    _S806[int(0)] = 0.0f;
    _S806[int(1)] = 0.0f;
    _S806[int(2)] = 0.0f;
    _S806[int(3)] = 0.0f;
    _S806[int(4)] = 0.0f;
    _S806[int(5)] = 0.0f;
    _S806[int(6)] = 0.0f;
    _S806[int(7)] = 0.0f;
    _S806[int(8)] = 0.0f;
    _S806[int(9)] = 0.0f;
    _S806[int(10)] = 0.0f;
    _S806[int(11)] = 0.0f;
    _S806[int(12)] = 0.0f;
    _S806[int(13)] = 0.0f;
    _S806[int(14)] = 0.0f;
    _S806[int(15)] = 0.0f;
    _S806[int(16)] = 0.0f;
    _S806[int(17)] = 0.0f;
    _S806[int(18)] = 0.0f;
    _S806[int(19)] = 0.0f;
    _S806[int(20)] = 0.0f;
    _S806[int(21)] = 0.0f;
    _S806[int(22)] = 0.0f;
    _S806[int(15)] = _S775;
    _S806[int(14)] = _S777;
    _S806[int(13)] = _S779;
    _S806[int(16)] = _S783.differential_0;
    _S806[int(12)] = _S781;
    _S806[int(21)] = _S788.differential_0;
    _S806[int(11)] = _S787;
    _S806[int(22)] = _S793.differential_0;
    _S806[int(10)] = _S792;
    _S806[int(9)] = _S792;
    _S806[int(20)] = _S799.differential_0;
    _S806[int(8)] = _S798;
    _S806[int(19)] = _S804.differential_0;
    _S806[int(7)] = _S803;
    float _S807 = _S806[int(0)];
    float _S808 = _S806[int(1)];
    float _S809 = _S806[int(2)];
    float _S810 = _S806[int(3)];
    float _S811 = _S806[int(4)];
    float _S812 = _S806[int(5)];
    float _S813 = _S806[int(6)];
    float _S814 = _S806[int(7)];
    float _S815 = _S806[int(8)];
    float _S816 = _S806[int(9)];
    float _S817 = _S806[int(10)];
    float _S818 = _S806[int(11)];
    float _S819 = _S806[int(12)];
    float _S820 = _S806[int(13)];
    float _S821 = _S806[int(14)];
    float _S822 = _S806[int(15)];
    float _S823 = _S806[int(16)];
    float _S824 = _S806[int(17)];
    float _S825 = _S806[int(18)];
    float _S826 = _S806[int(19)];
    float _S827 = _S806[int(20)];
    float _S828 = _S806[int(21)];
    float _S829 = _S806[int(22)];
    FixedArray<float, 23>  _S830;
    if(_S731)
    {
        float _S831 = _S732 * _S773;
        DiffPair_float_0 _S832;
        (&_S832)->primal_0 = _S733;
        (&_S832)->differential_0 = 0.0f;
        DiffPair_float_0 _S833;
        (&_S833)->primal_0 = 0.0f;
        (&_S833)->differential_0 = 0.0f;
        DiffPair_float_0 _S834;
        (&_S834)->primal_0 = 2.0f;
        (&_S834)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S832, &_S833, &_S834, _S831);
        float _S835 = - _S832.differential_0 / _S734;
        float _S836 = _S735 * - _S835;
        float _S837 = _S736 * _S835;
        DiffPair_float_0 _S838;
        (&_S838)->primal_0 = _S737;
        (&_S838)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S838, _S836);
        DiffPair_float_0 _S839;
        (&_S839)->primal_0 = 9.999999960041972e-13f;
        (&_S839)->differential_0 = 0.0f;
        DiffPair_float_0 _S840;
        (&_S840)->primal_0 = _S738;
        (&_S840)->differential_0 = 0.0f;
        _d_max_0(&_S839, &_S840, _S838.differential_0);
        float _S841 = _S739 * _S840.differential_0;
        float _S842 = _S740 * _S840.differential_0;
        float _S843 = - _S841 / _S741;
        float _S844 = _S743 * (_S726[int(18)] * _S843);
        float _S845 = - _S842 / _S741;
        float _S846 = _S745 * (_S726[int(18)] * _S845);
        float _S847 = - _S837 / _S741;
        float _S848 = _S726[int(18)] * _S847;
        float _S849 = _S844 + _S844 + _S745 * _S848;
        float _S850 = _S846 + _S846 + _S743 * _S848;
        float _S851 = _S742 * - _S843 + _S744 * - _S845 + _S746 * - _S847;
        FixedArray<float, 23>  _S852;
        _S852[int(0)] = 0.0f;
        _S852[int(1)] = 0.0f;
        _S852[int(2)] = 0.0f;
        _S852[int(3)] = 0.0f;
        _S852[int(4)] = 0.0f;
        _S852[int(5)] = 0.0f;
        _S852[int(6)] = 0.0f;
        _S852[int(7)] = 0.0f;
        _S852[int(8)] = 0.0f;
        _S852[int(9)] = 0.0f;
        _S852[int(10)] = 0.0f;
        _S852[int(11)] = 0.0f;
        _S852[int(12)] = 0.0f;
        _S852[int(13)] = 0.0f;
        _S852[int(14)] = 0.0f;
        _S852[int(15)] = 0.0f;
        _S852[int(16)] = 0.0f;
        _S852[int(17)] = 0.0f;
        _S852[int(18)] = 0.0f;
        _S852[int(19)] = 0.0f;
        _S852[int(20)] = 0.0f;
        _S852[int(21)] = 0.0f;
        _S852[int(22)] = 0.0f;
        _S852[int(5)] = _S841;
        _S852[int(4)] = _S842;
        _S852[int(3)] = _S849;
        _S852[int(2)] = _S850;
        _S852[int(6)] = _S837;
        float _S853 = _S808 + _S852[int(1)];
        float _S854 = _S809 + _S852[int(2)];
        float _S855 = _S810 + _S852[int(3)];
        float _S856 = _S811 + _S852[int(4)];
        float _S857 = _S812 + _S852[int(5)];
        float _S858 = _S813 + _S852[int(6)];
        float _S859 = _S814 + _S852[int(7)];
        float _S860 = _S815 + _S852[int(8)];
        float _S861 = _S816 + _S852[int(9)];
        float _S862 = _S817 + _S852[int(10)];
        float _S863 = _S818 + _S852[int(11)];
        float _S864 = _S819 + _S852[int(12)];
        float _S865 = _S820 + _S852[int(13)];
        float _S866 = _S821 + _S852[int(14)];
        float _S867 = _S822 + _S852[int(15)];
        float _S868 = _S823 + _S852[int(16)];
        float _S869 = _S824 + _S852[int(17)];
        float _S870 = _S825 + _S852[int(18)];
        float _S871 = _S826 + _S852[int(19)];
        float _S872 = _S827 + _S852[int(20)];
        float _S873 = _S828 + _S852[int(21)];
        float _S874 = _S829 + _S852[int(22)];
        _S830[int(0)] = _S807 + _S852[int(0)];
        _S830[int(1)] = _S853;
        _S830[int(2)] = _S854;
        _S830[int(3)] = _S855;
        _S830[int(4)] = _S856;
        _S830[int(5)] = _S857;
        _S830[int(6)] = _S858;
        _S830[int(7)] = _S859;
        _S830[int(8)] = _S860;
        _S830[int(9)] = _S861;
        _S830[int(10)] = _S862;
        _S830[int(11)] = _S863;
        _S830[int(12)] = _S864;
        _S830[int(13)] = _S865;
        _S830[int(14)] = _S866;
        _S830[int(15)] = _S867;
        _S830[int(16)] = _S868;
        _S830[int(17)] = _S869;
        _S830[int(18)] = _S870;
        _S830[int(19)] = _S871;
        _S830[int(20)] = _S872;
        _S830[int(21)] = _S873;
        _S830[int(22)] = _S874;
        _S732 = _S851;
    }
    else
    {
        _S830[int(0)] = _S807;
        _S830[int(1)] = _S808;
        _S830[int(2)] = _S809;
        _S830[int(3)] = _S810;
        _S830[int(4)] = _S811;
        _S830[int(5)] = _S812;
        _S830[int(6)] = _S813;
        _S830[int(7)] = _S814;
        _S830[int(8)] = _S815;
        _S830[int(9)] = _S816;
        _S830[int(10)] = _S817;
        _S830[int(11)] = _S818;
        _S830[int(12)] = _S819;
        _S830[int(13)] = _S820;
        _S830[int(14)] = _S821;
        _S830[int(15)] = _S822;
        _S830[int(16)] = _S823;
        _S830[int(17)] = _S824;
        _S830[int(18)] = _S825;
        _S830[int(19)] = _S826;
        _S830[int(20)] = _S827;
        _S830[int(21)] = _S828;
        _S830[int(22)] = _S829;
        _S732 = 0.0f;
    }
    if(_S730)
    {
        FixedArray<float, 23>  _S875;
        _S875[int(0)] = 0.0f;
        _S875[int(1)] = 0.0f;
        _S875[int(2)] = 0.0f;
        _S875[int(3)] = 0.0f;
        _S875[int(4)] = 0.0f;
        _S875[int(5)] = 0.0f;
        _S875[int(6)] = 0.0f;
        _S875[int(7)] = 0.0f;
        _S875[int(8)] = 0.0f;
        _S875[int(9)] = 0.0f;
        _S875[int(10)] = 0.0f;
        _S875[int(11)] = 0.0f;
        _S875[int(12)] = 0.0f;
        _S875[int(13)] = 0.0f;
        _S875[int(14)] = 0.0f;
        _S875[int(15)] = 0.0f;
        _S875[int(16)] = 0.0f;
        _S875[int(17)] = 0.0f;
        _S875[int(18)] = 0.0f;
        _S875[int(19)] = 0.0f;
        _S875[int(20)] = 0.0f;
        _S875[int(21)] = 0.0f;
        _S875[int(22)] = 0.0f;
        _S875[int(3)] = 0.0f;
        float _S876 = _S830[int(1)] + _S875[int(1)];
        float _S877 = _S830[int(2)] + _S875[int(2)];
        float _S878 = _S830[int(3)] + _S875[int(3)];
        float _S879 = _S830[int(4)] + _S875[int(4)];
        float _S880 = _S830[int(5)] + _S875[int(5)];
        float _S881 = _S830[int(6)] + _S875[int(6)];
        float _S882 = _S830[int(7)] + _S875[int(7)];
        float _S883 = _S830[int(8)] + _S875[int(8)];
        float _S884 = _S830[int(9)] + _S875[int(9)];
        float _S885 = _S830[int(10)] + _S875[int(10)];
        float _S886 = _S830[int(11)] + _S875[int(11)];
        float _S887 = _S830[int(12)] + _S875[int(12)];
        float _S888 = _S830[int(13)] + _S875[int(13)];
        float _S889 = _S830[int(14)] + _S875[int(14)];
        float _S890 = _S830[int(15)] + _S875[int(15)];
        float _S891 = _S830[int(16)] + _S875[int(16)];
        float _S892 = _S830[int(17)] + _S875[int(17)];
        float _S893 = _S830[int(18)] + _S875[int(18)];
        float _S894 = _S830[int(19)] + _S875[int(19)];
        float _S895 = _S830[int(20)] + _S875[int(20)];
        float _S896 = _S830[int(21)] + _S875[int(21)];
        float _S897 = _S830[int(22)] + _S875[int(22)];
        _S830[int(0)] = _S830[int(0)] + _S875[int(0)];
        _S830[int(1)] = _S876;
        _S830[int(2)] = _S877;
        _S830[int(3)] = _S878;
        _S830[int(4)] = _S879;
        _S830[int(5)] = _S880;
        _S830[int(6)] = _S881;
        _S830[int(7)] = _S882;
        _S830[int(8)] = _S883;
        _S830[int(9)] = _S884;
        _S830[int(10)] = _S885;
        _S830[int(11)] = _S886;
        _S830[int(12)] = _S887;
        _S830[int(13)] = _S888;
        _S830[int(14)] = _S889;
        _S830[int(15)] = _S890;
        _S830[int(16)] = _S891;
        _S830[int(17)] = _S892;
        _S830[int(18)] = _S893;
        _S830[int(19)] = _S894;
        _S830[int(20)] = _S895;
        _S830[int(21)] = _S896;
        _S830[int(22)] = _S897;
    }
    float _S898 = -10.0f * _S772;
    DiffPair_float_0 _S899;
    (&_S899)->primal_0 = _S729;
    (&_S899)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S899, _S898);
    float _S900 = _S899.differential_0 / _S728;
    float _S901 = _S727 * _S900;
    float _S902 = _S771 / _S728;
    float _S903 = _S727 * _S902;
    float _S904 = _S726[int(1)] * - _S900 + _S726[int(0)] * - _S902;
    DiffPair_float_0 _S905;
    (&_S905)->primal_0 = _S726[int(17)];
    (&_S905)->differential_0 = 0.0f;
    DiffPair_float_0 _S906;
    (&_S906)->primal_0 = 1.0f;
    (&_S906)->differential_0 = 0.0f;
    _d_max_0(&_S905, &_S906, _S904);
    FixedArray<float, 23>  _S907;
    _S907[int(0)] = 0.0f;
    _S907[int(1)] = 0.0f;
    _S907[int(2)] = 0.0f;
    _S907[int(3)] = 0.0f;
    _S907[int(4)] = 0.0f;
    _S907[int(5)] = 0.0f;
    _S907[int(6)] = 0.0f;
    _S907[int(7)] = 0.0f;
    _S907[int(8)] = 0.0f;
    _S907[int(9)] = 0.0f;
    _S907[int(10)] = 0.0f;
    _S907[int(11)] = 0.0f;
    _S907[int(12)] = 0.0f;
    _S907[int(13)] = 0.0f;
    _S907[int(14)] = 0.0f;
    _S907[int(15)] = 0.0f;
    _S907[int(16)] = 0.0f;
    _S907[int(17)] = 0.0f;
    _S907[int(18)] = 0.0f;
    _S907[int(19)] = 0.0f;
    _S907[int(20)] = 0.0f;
    _S907[int(21)] = 0.0f;
    _S907[int(22)] = 0.0f;
    _S907[int(18)] = _S732;
    _S907[int(1)] = _S901;
    _S907[int(17)] = _S905.differential_0;
    _S907[int(0)] = _S903;
    FixedArray<float, 23>  _S908 = {
        _S830[int(0)] + _S907[int(0)], _S830[int(1)] + _S907[int(1)], _S830[int(2)] + _S907[int(2)], _S830[int(3)] + _S907[int(3)], _S830[int(4)] + _S907[int(4)], _S830[int(5)] + _S907[int(5)], _S830[int(6)] + _S907[int(6)], _S830[int(7)] + _S907[int(7)], _S830[int(8)] + _S907[int(8)], _S830[int(9)] + _S907[int(9)], _S830[int(10)] + _S907[int(10)], _S830[int(11)] + _S907[int(11)], _S830[int(12)] + _S907[int(12)], _S830[int(13)] + _S907[int(13)], _S830[int(14)] + _S907[int(14)], _S830[int(15)] + _S907[int(15)], _S830[int(16)] + _S907[int(16)], _S830[int(17)] + _S907[int(17)], _S830[int(18)] + _S907[int(18)], _S830[int(19)] + _S907[int(19)], _S830[int(20)] + _S907[int(20)], _S830[int(21)] + _S907[int(21)], _S830[int(22)] + _S907[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S908;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S909, FixedArray<float, 10>  * _S910, FixedArray<float, 10>  * _S911)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S909, _S910, _S911);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 10>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S912)
{
    FixedArray<float, 23>  _S913 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S913;
    FixedArray<float, 10>  _S914 = weights_5;
    FixedArray<float, 10>  _S915 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S914, &_S915);
    *_S912 = (&dp_raw_losses_0)->differential_0;
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

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_1(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S916, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S917, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S918, float3  _S919)
{
    _d_clamp_vector_0(_S916, _S917, _S918, _S919);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S920 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S921 = make_float3 (_S920);
    float3  _S922 = make_float3 (0.0f);
    float3  _S923 = make_float3 (1.0f);
    float3  _S924 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S925;
    (&_S925)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S920) * (*dpbackground_0).primal_0;
    (&_S925)->differential_0 = _S924;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S926;
    (&_S926)->primal_0 = _S922;
    (&_S926)->differential_0 = _S924;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S927;
    (&_S927)->primal_0 = _S923;
    (&_S927)->differential_0 = _S924;
    s_bwd_prop_clamp_1(&_S925, &_S926, &_S927, _s_dOut_6);
    float3  _S928 = _S921 * _S925.differential_0;
    float3  _S929 = (*dpbackground_0).primal_0 * _S925.differential_0;
    float _S930 = - (_S929.x + _S929.y + _S929.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S928;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S930;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S925.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S931, DiffPair_float_0 * _S932, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S933, float3  _S934)
{
    s_bwd_prop_blend_background_0(_S931, _S932, _S933, _S934);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S935 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S935;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S935;
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
        DiffPair_float_0 _S936 = *dpx_16;
        float _S937 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_20;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S937;
        float _S938 = val_0 * (F32_log((_S936.primal_0))) * dOut_20;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S938;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S939 = rgb_2;
    float _S940;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S940 = _S939.x * 12.92000007629394531f;
    }
    else
    {
        _S940 = 1.0549999475479126f * (F32_pow((_S939.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S939)->x) = _S940;
    if((_S939.y) < 0.00313080009073019f)
    {
        _S940 = _S939.y * 12.92000007629394531f;
    }
    else
    {
        _S940 = 1.0549999475479126f * (F32_pow((_S939.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S939)->y) = _S940;
    if((_S939.z) < 0.00313080009073019f)
    {
        _S940 = _S939.z * 12.92000007629394531f;
    }
    else
    {
        _S940 = 1.0549999475479126f * (F32_pow((_S939.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S939)->z) = _S940;
    return _S939;
}

inline __device__ float s_primal_ctx_pow_0(float _S941, float _S942)
{
    return (F32_pow((_S941), (_S942)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S943, DiffPair_float_0 * _S944, float _S945)
{
    _d_pow_0(_S943, _S944, _S945);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_7)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S946 = *dprgb_1;
    float _S947 = (*dprgb_1).primal_0.x;
    bool _S948 = _S947 < 0.00313080009073019f;
    float _S949;
    if(_S948)
    {
        _S949 = _S947 * 12.92000007629394531f;
    }
    else
    {
        _S949 = 1.0549999475479126f * s_primal_ctx_pow_0(_S947, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S950 = _S946.primal_0;
    *&((&_S950)->x) = _S949;
    float _S951 = _S950.y;
    bool _S952 = _S951 < 0.00313080009073019f;
    if(_S952)
    {
        _S949 = _S951 * 12.92000007629394531f;
    }
    else
    {
        _S949 = 1.0549999475479126f * s_primal_ctx_pow_0(_S951, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S950)->y) = _S949;
    float _S953 = _S950.z;
    bool _S954 = _S953 < 0.00313080009073019f;
    _S950 = _s_dOut_7;
    *&((&_S950)->z) = 0.0f;
    if(_S954)
    {
        _S949 = 12.92000007629394531f * _s_dOut_7.z;
    }
    else
    {
        float _S955 = 1.0549999475479126f * _s_dOut_7.z;
        DiffPair_float_0 _S956;
        (&_S956)->primal_0 = _S953;
        (&_S956)->differential_0 = 0.0f;
        DiffPair_float_0 _S957;
        (&_S957)->primal_0 = 0.4166666567325592f;
        (&_S957)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S956, &_S957, _S955);
        _S949 = _S956.differential_0;
    }
    float3  _S958 = _S950 + make_float3 (0.0f, 0.0f, _S949);
    _S950 = _S958;
    *&((&_S950)->y) = 0.0f;
    if(_S952)
    {
        _S949 = 12.92000007629394531f * _S958.y;
    }
    else
    {
        float _S959 = 1.0549999475479126f * _S958.y;
        DiffPair_float_0 _S960;
        (&_S960)->primal_0 = _S951;
        (&_S960)->differential_0 = 0.0f;
        DiffPair_float_0 _S961;
        (&_S961)->primal_0 = 0.4166666567325592f;
        (&_S961)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S960, &_S961, _S959);
        _S949 = _S960.differential_0;
    }
    float3  _S962 = _S950 + make_float3 (0.0f, _S949, 0.0f);
    _S950 = _S962;
    *&((&_S950)->x) = 0.0f;
    if(_S948)
    {
        _S949 = 12.92000007629394531f * _S962.x;
    }
    else
    {
        float _S963 = 1.0549999475479126f * _S962.x;
        DiffPair_float_0 _S964;
        (&_S964)->primal_0 = _S947;
        (&_S964)->differential_0 = 0.0f;
        DiffPair_float_0 _S965;
        (&_S965)->primal_0 = 0.4166666567325592f;
        (&_S965)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S964, &_S965, _S963);
        _S949 = _S964.differential_0;
    }
    float3  _S966 = _S950 + make_float3 (_S949, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S966;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S967, float3  _S968)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S967, _S968);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S969 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S969;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S970 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_6 = _S970;
    FixedArray<float, 10>  _S971 = dist_coeffs_6;
    bool _S972 = undistort_point_0(_S970, &_S971, int(12), &uv_6);
    if(!_S972)
    {
        int3  _S973 = make_int3 (int(0));
        float3  _S974 = make_float3 ((float)_S973.x, (float)_S973.y, (float)_S973.z);
        return _S974;
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S975, float3  _S976)
{
    return cross_0(_S975, _S976);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_8)
{
    float _S977 = (*dpx_17).primal_0.x;
    float _S978 = (*dpx_17).primal_0.y;
    float _S979 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S980;
    (&_S980)->primal_0 = _S977 * _S977 + _S978 * _S978 + _S979 * _S979;
    (&_S980)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S980, _s_dOut_8);
    float _S981 = (*dpx_17).primal_0.z * _S980.differential_0;
    float _S982 = _S981 + _S981;
    float _S983 = (*dpx_17).primal_0.y * _S980.differential_0;
    float _S984 = _S983 + _S983;
    float _S985 = (*dpx_17).primal_0.x * _S980.differential_0;
    float _S986 = _S985 + _S985;
    float3  _S987 = make_float3 (0.0f);
    *&((&_S987)->z) = _S982;
    *&((&_S987)->y) = _S984;
    *&((&_S987)->x) = _S986;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S987;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S988, float _S989)
{
    s_bwd_prop_length_impl_1(_S988, _S989);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S990, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S991, float3  _S992)
{
    _d_cross_0(_S990, _S991, _S992);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S993 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S994 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S995 = s_primal_ctx_cross_0(dx_0, _S994);
    bool _S996 = (s_primal_ctx_dot_0(_S995, _S995)) != 0.0f;
    float3  _S997;
    float3  _S998;
    if(_S996)
    {
        float _S999 = length_2(_S995);
        float3  _S1000 = make_float3 (_S999);
        _S997 = make_float3 (_S999 * _S999);
        _S998 = _S1000;
    }
    else
    {
        _S997 = _S993;
        _S998 = _S993;
    }
    if(_S996)
    {
        float3  _S1001 = _s_dOut_9 / _S997;
        float3  _S1002 = _S995 * - _S1001;
        float3  _S1003 = _S998 * _S1001;
        float _S1004 = _S1002.x + _S1002.y + _S1002.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1005;
        (&_S1005)->primal_0 = _S995;
        (&_S1005)->differential_0 = _S993;
        s_bwd_length_impl_1(&_S1005, _S1004);
        _S997 = _S1003 + _S1005.differential_0;
    }
    else
    {
        _S997 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1006;
    (&_S1006)->primal_0 = _S995;
    (&_S1006)->differential_0 = _S993;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1007;
    (&_S1007)->primal_0 = _S995;
    (&_S1007)->differential_0 = _S993;
    s_bwd_prop_dot_0(&_S1006, &_S1007, 0.0f);
    float3  _S1008 = _S1007.differential_0 + _S1006.differential_0 + _S997;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1009;
    (&_S1009)->primal_0 = dx_0;
    (&_S1009)->differential_0 = _S993;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1010;
    (&_S1010)->primal_0 = _S994;
    (&_S1010)->differential_0 = _S993;
    s_bwd_prop_cross_0(&_S1009, &_S1010, _S1008);
    float3  s_diff_dy_T_0 = - _S1010.differential_0;
    float3  _S1011 = - s_diff_dy_T_0;
    float3  _S1012 = - _S1009.differential_0;
    FixedArray<float3 , 4>  _S1013;
    _S1013[int(0)] = _S993;
    _S1013[int(1)] = _S993;
    _S1013[int(2)] = _S993;
    _S1013[int(3)] = _S993;
    _S1013[int(2)] = _S1011;
    _S1013[int(3)] = s_diff_dy_T_0;
    _S1013[int(0)] = _S1012;
    _S1013[int(1)] = _S1009.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S1013;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S1014, float3  _S1015)
{
    s_bwd_prop_points_to_normal_0(_S1014, _S1015);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S1016 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S1016;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_7, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S1017 = float2 {intrins_1.z, intrins_1.w};
    float2  _S1018 = float2 {intrins_1.x, intrins_1.y};
    float2  _S1019 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S1017) / _S1018;
    float2  uv_7 = _S1019;
    FixedArray<float, 10>  _S1020 = dist_coeffs_7;
    bool _S1021 = undistort_point_0(_S1019, &_S1020, int(12), &uv_7);
    if(!_S1021)
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
    float2  _S1022 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S1017) / _S1018;
    float2  uv_8 = _S1022;
    FixedArray<float, 10>  _S1023 = dist_coeffs_7;
    bool _S1024 = undistort_point_0(_S1022, &_S1023, int(12), &uv_8);
    if(!_S1024)
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
    float2  _S1025 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S1017) / _S1018;
    float2  uv_9 = _S1025;
    FixedArray<float, 10>  _S1026 = dist_coeffs_7;
    bool _S1027 = undistort_point_0(_S1025, &_S1026, int(12), &uv_9);
    if(!_S1027)
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
    float2  _S1028 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S1017) / _S1018;
    float2  uv_10 = _S1028;
    FixedArray<float, 10>  _S1029 = dist_coeffs_7;
    bool _S1030 = undistort_point_0(_S1028, &_S1029, int(12), &uv_10);
    if(!_S1030)
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
    float2  _S1031;
    bool _S1032;
    float2  _S1033;
    bool _S1034;
    float2  _S1035;
    bool _S1036;
    float2  _S1037;
    bool _S1038;
};

inline __device__ float s_primal_ctx_sin_0(float _S1039)
{
    return (F32_sin((_S1039)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1040)
{
    return (F32_cos((_S1040)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1041 = make_float2 (0.0f);
    _s_diff_ctx_0->_S1031 = _S1041;
    _s_diff_ctx_0->_S1032 = false;
    _s_diff_ctx_0->_S1033 = _S1041;
    _s_diff_ctx_0->_S1034 = false;
    _s_diff_ctx_0->_S1035 = _S1041;
    _s_diff_ctx_0->_S1036 = false;
    _s_diff_ctx_0->_S1037 = _S1041;
    _s_diff_ctx_0->_S1038 = false;
    _s_diff_ctx_0->_S1033 = _S1041;
    _s_diff_ctx_0->_S1034 = false;
    _s_diff_ctx_0->_S1035 = _S1041;
    _s_diff_ctx_0->_S1036 = false;
    _s_diff_ctx_0->_S1037 = _S1041;
    _s_diff_ctx_0->_S1038 = false;
    float3  _S1042 = make_float3 (0.0f);
    float2  _S1043 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1044 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1045 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1043) / _S1044;
    float2  _S1046 = _S1045;
    bool _S1047 = undistort_point_0(_S1045, dist_coeffs_8, int(12), &_S1046);
    _s_diff_ctx_0->_S1031 = _S1046;
    _s_diff_ctx_0->_S1032 = _S1047;
    float2  uv_11 = _S1046;
    bool _S1048 = !_S1047;
    float3  normal_4;
    if(_S1048)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1049 = !_S1048;
    int _S1050;
    FixedArray<float3 , 4>  points_3;
    if(_S1049)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1051 = length_1(uv_11);
            float3  raydir_19 = make_float3 ((uv_11 / make_float2 ((F32_max((_S1051), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1051))).x, (uv_11 / make_float2 ((F32_max((_S1051), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1051))).y, s_primal_ctx_cos_0(_S1051));
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
        float3  _S1052 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1053 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1043) / _S1044;
        float2  _S1054 = _S1053;
        bool _S1055 = undistort_point_0(_S1053, dist_coeffs_8, int(12), &_S1054);
        _s_diff_ctx_0->_S1033 = _S1054;
        _s_diff_ctx_0->_S1034 = _S1055;
        float2  uv_12 = _S1054;
        bool _S1056 = !_S1055;
        if(_S1056)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1057 = !_S1056;
        if(_S1057)
        {
            if(is_fisheye_6)
            {
                float _S1058 = length_1(uv_12);
                float3  raydir_21 = make_float3 ((uv_12 / make_float2 ((F32_max((_S1058), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1058))).x, (uv_12 / make_float2 ((F32_max((_S1058), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1058))).y, s_primal_ctx_cos_0(_S1058));
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
            float3  _S1059 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1050 = int(2);
            points_3[int(0)] = _S1052;
            points_3[int(1)] = _S1059;
            points_3[int(2)] = _S1042;
            points_3[int(3)] = _S1042;
        }
        else
        {
            _S1050 = int(0);
            points_3[int(0)] = _S1052;
            points_3[int(1)] = _S1042;
            points_3[int(2)] = _S1042;
            points_3[int(3)] = _S1042;
        }
        bool _runFlag_0;
        if(_S1050 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1049;
            _S1050 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1060 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1043) / _S1044;
            float2  _S1061 = _S1060;
            bool _S1062 = undistort_point_0(_S1060, dist_coeffs_8, int(12), &_S1061);
            _s_diff_ctx_0->_S1035 = _S1061;
            _s_diff_ctx_0->_S1036 = _S1062;
            float2  uv_13 = _S1061;
            if(!_S1062)
            {
                float3  _S1063 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1050 = int(0);
                normal_4 = _S1063;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1064 = length_1(uv_13);
                    float3  raydir_23 = make_float3 ((uv_13 / make_float2 ((F32_max((_S1064), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1064))).x, (uv_13 / make_float2 ((F32_max((_S1064), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1064))).y, s_primal_ctx_cos_0(_S1064));
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
                float2  _S1065 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1043) / _S1044;
                float2  _S1066 = _S1065;
                bool _S1067 = undistort_point_0(_S1065, dist_coeffs_8, int(12), &_S1066);
                _s_diff_ctx_0->_S1037 = _S1066;
                _s_diff_ctx_0->_S1038 = _S1067;
                float2  uv_14 = _S1066;
                bool _S1068 = !_S1067;
                if(_S1068)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1069 = !_S1068;
                int _S1070;
                if(_S1069)
                {
                    if(is_fisheye_6)
                    {
                        float _S1071 = length_1(uv_14);
                        float3  raydir_25 = make_float3 ((uv_14 / make_float2 ((F32_max((_S1071), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1071))).x, (uv_14 / make_float2 ((F32_max((_S1071), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1071))).y, s_primal_ctx_cos_0(_S1071));
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
                    _S1070 = int(2);
                }
                else
                {
                    _S1070 = int(0);
                }
                if(_S1070 != int(2))
                {
                    _runFlag_0 = false;
                    _S1050 = _S1070;
                }
                if(_runFlag_0)
                {
                    _S1050 = int(1);
                }
            }
        }
    }
    else
    {
        _S1050 = int(0);
        points_3[int(0)] = _S1042;
        points_3[int(1)] = _S1042;
        points_3[int(2)] = _S1042;
        points_3[int(3)] = _S1042;
    }
    if(!(_S1050 != int(1)))
    {
        float3  _S1072 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1072, _S1072)) != 0.0f)
        {
            normal_4 = _S1072 / make_float3 (length_2(_S1072));
        }
        else
        {
            normal_4 = _S1072;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1073 = *dpdepths_1;
    float3  _S1074 = make_float3 (0.0f);
    float2  _S1075 = _s_diff_ctx_1->_S1031;
    bool _S1076 = !!_s_diff_ctx_1->_S1032;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1077;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1078;
    if(_S1076)
    {
        if(is_fisheye_7)
        {
            float _S1079 = length_1(_S1075);
            float3  raydir_31 = make_float3 ((_S1075 / make_float2 ((F32_max((_S1079), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1079))).x, (_S1075 / make_float2 ((F32_max((_S1079), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1079))).y, s_primal_ctx_cos_0(_S1079));
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
            float3  raydir_32 = make_float3 (_S1075.x, _S1075.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1080 = make_float3 (_S1073.primal_0.x) * raydir_27;
        float2  _S1081 = _s_diff_ctx_1->_S1033;
        bool _S1082 = !!_s_diff_ctx_1->_S1034;
        if(_S1082)
        {
            if(is_fisheye_7)
            {
                float _S1083 = length_1(_S1081);
                float3  raydir_33 = make_float3 ((_S1081 / make_float2 ((F32_max((_S1083), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1083))).x, (_S1081 / make_float2 ((F32_max((_S1083), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1083))).y, s_primal_ctx_cos_0(_S1083));
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
                float3  raydir_34 = make_float3 (_S1081.x, _S1081.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1084 = make_float3 (_S1073.primal_0.y) * raydir_28;
            _S1077 = int(2);
            points_4[int(0)] = _S1080;
            points_4[int(1)] = _S1084;
            points_4[int(2)] = _S1074;
            points_4[int(3)] = _S1074;
        }
        else
        {
            _S1077 = int(0);
            points_4[int(0)] = _S1080;
            points_4[int(1)] = _S1074;
            points_4[int(2)] = _S1074;
            points_4[int(3)] = _S1074;
            raydir_28 = _S1074;
        }
        if(_S1077 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1076;
            _S1077 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S1085 = _s_diff_ctx_1->_S1035;
            if(!_s_diff_ctx_1->_S1036)
            {
                _runFlag_2 = false;
                _S1077 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1086 = length_1(_S1085);
                    float3  raydir_35 = make_float3 ((_S1085 / make_float2 ((F32_max((_S1086), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1086))).x, (_S1085 / make_float2 ((F32_max((_S1086), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1086))).y, s_primal_ctx_cos_0(_S1086));
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
                    float3  raydir_36 = make_float3 (_S1085.x, _S1085.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1073.primal_0.z) * raydir_29;
                float2  _S1087 = _s_diff_ctx_1->_S1037;
                bool _S1088 = !!_s_diff_ctx_1->_S1038;
                int _S1089;
                if(_S1088)
                {
                    if(is_fisheye_7)
                    {
                        float _S1090 = length_1(_S1087);
                        float3  raydir_37 = make_float3 ((_S1087 / make_float2 ((F32_max((_S1090), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1090))).x, (_S1087 / make_float2 ((F32_max((_S1090), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1090))).y, s_primal_ctx_cos_0(_S1090));
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
                        float3  raydir_38 = make_float3 (_S1087.x, _S1087.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1073.primal_0.w) * raydir_30;
                    _S1089 = int(2);
                }
                else
                {
                    _S1089 = int(0);
                    raydir_30 = _S1074;
                }
                if(_S1089 != int(2))
                {
                    _runFlag_3 = false;
                    _S1077 = _S1089;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1077 = int(1);
                }
                float3  _S1091 = raydir_29;
                _runFlag_3 = _S1088;
                raydir_29 = raydir_30;
                raydir_30 = _S1091;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1074;
                raydir_30 = _S1074;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1074;
            raydir_30 = _S1074;
        }
        float3  _S1092 = raydir_27;
        float3  _S1093 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1078 = _S1082;
        raydir_29 = _S1093;
        raydir_30 = _S1092;
    }
    else
    {
        _S1077 = int(0);
        points_4[int(0)] = _S1074;
        points_4[int(1)] = _S1074;
        points_4[int(2)] = _S1074;
        points_4[int(3)] = _S1074;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1074;
        raydir_28 = _S1074;
        _S1078 = false;
        raydir_29 = _S1074;
        raydir_30 = _S1074;
    }
    bool _S1094 = !(_S1077 != int(1));
    float3  _S1095;
    float3  _S1096;
    float3  _S1097;
    float3  _S1098;
    float3  _S1099;
    bool _S1100;
    if(_S1094)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1101 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1102 = s_primal_ctx_cross_0(dx_1, _S1101);
        bool _S1103 = (s_primal_ctx_dot_0(_S1102, _S1102)) != 0.0f;
        if(_S1103)
        {
            float _S1104 = length_2(_S1102);
            float3  _S1105 = make_float3 (_S1104);
            _S1095 = make_float3 (_S1104 * _S1104);
            _S1096 = _S1105;
        }
        else
        {
            _S1095 = _S1074;
            _S1096 = _S1074;
        }
        float3  _S1106 = _S1096;
        _S1100 = _S1103;
        _S1096 = _S1102;
        _S1097 = _S1106;
        _S1098 = dx_1;
        _S1099 = _S1101;
    }
    else
    {
        _S1100 = false;
        _S1095 = _S1074;
        _S1096 = _S1074;
        _S1097 = _S1074;
        _S1098 = _S1074;
        _S1099 = _S1074;
    }
    float4  _S1107 = make_float4 (0.0f);
    if(_S1094)
    {
        if(_S1100)
        {
            float3  _S1108 = _s_dOut_10 / _S1095;
            float3  _S1109 = _S1096 * - _S1108;
            float3  _S1110 = _S1097 * _S1108;
            float _S1111 = _S1109.x + _S1109.y + _S1109.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1112;
            (&_S1112)->primal_0 = _S1096;
            (&_S1112)->differential_0 = _S1074;
            s_bwd_length_impl_1(&_S1112, _S1111);
            _S1095 = _S1110 + _S1112.differential_0;
        }
        else
        {
            _S1095 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1113;
        (&_S1113)->primal_0 = _S1096;
        (&_S1113)->differential_0 = _S1074;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1114;
        (&_S1114)->primal_0 = _S1096;
        (&_S1114)->differential_0 = _S1074;
        s_bwd_prop_dot_0(&_S1113, &_S1114, 0.0f);
        float3  _S1115 = _S1114.differential_0 + _S1113.differential_0 + _S1095;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1116;
        (&_S1116)->primal_0 = _S1098;
        (&_S1116)->differential_0 = _S1074;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1117;
        (&_S1117)->primal_0 = _S1099;
        (&_S1117)->differential_0 = _S1074;
        s_bwd_prop_cross_0(&_S1116, &_S1117, _S1115);
        float3  s_diff_dy_T_1 = - _S1117.differential_0;
        float3  _S1118 = - s_diff_dy_T_1;
        float3  _S1119 = - _S1116.differential_0;
        FixedArray<float3 , 4>  _S1120;
        _S1120[int(0)] = _S1074;
        _S1120[int(1)] = _S1074;
        _S1120[int(2)] = _S1074;
        _S1120[int(3)] = _S1074;
        _S1120[int(2)] = _S1118;
        _S1120[int(3)] = s_diff_dy_T_1;
        _S1120[int(0)] = _S1119;
        _S1120[int(1)] = _S1116.differential_0;
        points_4[int(0)] = _S1120[int(0)];
        points_4[int(1)] = _S1120[int(1)];
        points_4[int(2)] = _S1120[int(2)];
        points_4[int(3)] = _S1120[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1074;
        points_4[int(1)] = _S1074;
        points_4[int(2)] = _S1074;
        points_4[int(3)] = _S1074;
    }
    float4  _S1121;
    if(_S1076)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1122 = points_4;
                FixedArray<float3 , 4>  _S1123 = points_4;
                FixedArray<float3 , 4>  _S1124 = points_4;
                FixedArray<float3 , 4>  _S1125 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1126 = raydir_27 * _S1125[int(3)];
                    float _S1127 = _S1126.x + _S1126.y + _S1126.z;
                    float4  _S1128 = _S1107;
                    *&((&_S1128)->w) = _S1127;
                    points_4[int(0)] = _S1122[int(0)];
                    points_4[int(1)] = _S1123[int(1)];
                    points_4[int(2)] = _S1124[int(2)];
                    points_4[int(3)] = _S1074;
                    _S1121 = _S1128;
                }
                else
                {
                    points_4[int(0)] = _S1122[int(0)];
                    points_4[int(1)] = _S1123[int(1)];
                    points_4[int(2)] = _S1124[int(2)];
                    points_4[int(3)] = _S1125[int(3)];
                    _S1121 = _S1107;
                }
                float3  _S1129 = raydir_28 * points_4[int(2)];
                float _S1130 = _S1129.x + _S1129.y + _S1129.z;
                FixedArray<float3 , 4>  _S1131 = points_4;
                FixedArray<float3 , 4>  _S1132 = points_4;
                float4  _S1133 = _S1107;
                *&((&_S1133)->z) = _S1130;
                float4  _S1134 = _S1121 + _S1133;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1131[int(1)];
                points_4[int(2)] = _S1074;
                points_4[int(3)] = _S1132[int(3)];
                _S1121 = _S1134;
            }
            else
            {
                FixedArray<float3 , 4>  _S1135 = points_4;
                FixedArray<float3 , 4>  _S1136 = points_4;
                FixedArray<float3 , 4>  _S1137 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1135[int(1)];
                points_4[int(2)] = _S1136[int(2)];
                points_4[int(3)] = _S1137[int(3)];
                _S1121 = _S1107;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1138 = points_4;
            FixedArray<float3 , 4>  _S1139 = points_4;
            FixedArray<float3 , 4>  _S1140 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1138[int(1)];
            points_4[int(2)] = _S1139[int(2)];
            points_4[int(3)] = _S1140[int(3)];
            _S1121 = _S1107;
        }
        if(_S1078)
        {
            FixedArray<float3 , 4>  _S1141 = points_4;
            float3  _S1142 = raydir_29 * points_4[int(1)];
            float _S1143 = _S1142.x + _S1142.y + _S1142.z;
            float4  _S1144 = _S1107;
            *&((&_S1144)->y) = _S1143;
            float4  _S1145 = _S1121 + _S1144;
            points_4[int(0)] = _S1074;
            points_4[int(1)] = _S1074;
            points_4[int(2)] = _S1074;
            points_4[int(3)] = _S1074;
            raydir_27 = _S1141[int(0)];
            _S1121 = _S1145;
        }
        else
        {
            FixedArray<float3 , 4>  _S1146 = points_4;
            FixedArray<float3 , 4>  _S1147 = points_4;
            FixedArray<float3 , 4>  _S1148 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1146[int(1)];
            points_4[int(2)] = _S1147[int(2)];
            points_4[int(3)] = _S1148[int(3)];
            raydir_27 = _S1074;
        }
        float3  _S1149 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1150 = _S1149.x + _S1149.y + _S1149.z;
        float4  _S1151 = _S1107;
        *&((&_S1151)->x) = _S1150;
        _S1121 = _S1121 + _S1151;
    }
    else
    {
        _S1121 = _S1107;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1121;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1152, float4  _S1153, FixedArray<float, 10>  * _S1154, bool _S1155, bool _S1156, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1157, float3  _S1158)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1159;
    float3  _S1160 = s_primal_ctx_depth_to_normal_0(_S1152, _S1153, _S1154, _S1155, _S1156, (*_S1157).primal_0, &_S1159);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1161 = _S1159;
    s_bwd_prop_depth_to_normal_0(_S1152, _S1153, _S1154, _S1155, _S1156, _S1157, _S1158, &_S1161);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_10, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1162 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1162;
    FixedArray<float, 10>  _S1163 = dist_coeffs_10;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S1163, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_11, bool is_fisheye_9)
{
    float2  _S1164 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_15 = _S1164;
    FixedArray<float, 10>  _S1165 = dist_coeffs_11;
    bool _S1166 = undistort_point_0(_S1164, &_S1165, int(12), &uv_15);
    if(!_S1166)
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
    float _S1167 = (F32_exp2(((*dpx_18).primal_0))) * 0.69314718246459961f * dOut_21;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1167;
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
    PPISPParams_0 _S1168 = p_0;
    float _S1169 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1170 = (pix_coord_0.x - image_center_0.x) / _S1169;
    float _S1171 = (pix_coord_0.y - image_center_0.y) / _S1169;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_2 = _S1170 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S1171 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_5 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_5 * r2_5;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_5) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1170 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S1171 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_6 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_6 * r2_6;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_6) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1170 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S1171 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_7 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_7 * r2_7;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_7) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_7 + 1.0f, 0.0f, 1.0f);
    float3  _S1172 = rgb_out_0;
    float2  bd_0 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S1173 = 0.3333333432674408f + nd_0.x;
    float _S1174 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1174, 1.0f, 0.0f, - _S1173, - _S1174, _S1173, 0.0f), T_0);
    float3  r0_0 = make_float3 (M_3.rows[int(0)].x, M_3.rows[int(0)].y, M_3.rows[int(0)].z);
    float3  r1_0 = make_float3 (M_3.rows[int(1)].x, M_3.rows[int(1)].y, M_3.rows[int(1)].z);
    float3  r2_8 = make_float3 (M_3.rows[int(2)].x, M_3.rows[int(2)].y, M_3.rows[int(2)].z);
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
    float _S1175 = _S1172.x;
    float _S1176 = _S1172.y;
    float intensity_0 = _S1175 + _S1176 + _S1172.z;
    float3  rgi_out_0 = mul_0(H_1, make_float3 (_S1175, _S1176, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1177 = rgi_out_1.x;
    float _S1178 = rgi_out_1.y;
    float3  _S1179 = clamp_1(make_float3 (_S1177, _S1178, rgi_out_1.z - _S1177 - _S1178), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1180 = _S1179.x;
    float _S1181 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(0)].toe_0))))));
    float _S1182 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(0)].shoulder_0))))));
    float _S1183 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(0)].gamma_0))))));
    float _S1184 = 1.0f / (1.0f + (F32_exp((- _S1168.crf_params_1[int(0)].center_0))));
    float a_1 = _S1182 * _S1184 / lerp_0(_S1181, _S1182, _S1184);
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1180 <= _S1184)
    {
        y_10 = a_1 * (F32_pow((_S1180 / _S1184), (_S1181)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1180) / (1.0f - _S1184)), (_S1182)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1183)));
    float _S1185 = _S1179.y;
    float _S1186 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(1)].toe_0))))));
    float _S1187 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(1)].shoulder_0))))));
    float _S1188 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(1)].gamma_0))))));
    float _S1189 = 1.0f / (1.0f + (F32_exp((- _S1168.crf_params_1[int(1)].center_0))));
    float a_2 = _S1187 * _S1189 / lerp_0(_S1186, _S1187, _S1189);
    float b_3 = 1.0f - a_2;
    if(_S1185 <= _S1189)
    {
        y_10 = a_2 * (F32_pow((_S1185 / _S1189), (_S1186)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1185) / (1.0f - _S1189)), (_S1187)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1188)));
    float _S1190 = _S1179.z;
    float _S1191 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(2)].toe_0))))));
    float _S1192 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(2)].shoulder_0))))));
    float _S1193 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1168.crf_params_1[int(2)].gamma_0))))));
    float _S1194 = 1.0f / (1.0f + (F32_exp((- _S1168.crf_params_1[int(2)].center_0))));
    float a_3 = _S1192 * _S1194 / lerp_0(_S1191, _S1192, _S1194);
    float b_4 = 1.0f - a_3;
    if(_S1190 <= _S1194)
    {
        y_10 = a_3 * (F32_pow((_S1190 / _S1194), (_S1191)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1190) / (1.0f - _S1194)), (_S1192)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1193)));
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
    PPISPParamsRQS_0 _S1195 = p_1;
    float _S1196 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S1197 = (pix_coord_1.x - image_center_1.x) / _S1196;
    float _S1198 = (pix_coord_1.y - image_center_1.y) / _S1196;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_5 = _S1197 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S1198 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_9 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_9 * r2_9;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_9) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float dx_6 = _S1197 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S1198 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_10 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_10 * r2_10;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_10) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_10 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S1197 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S1198 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_11 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_11 * r2_11;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_11) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_11 + 1.0f, 0.0f, 1.0f);
    float3  _S1199 = rgb_out_2;
    float2  bd_1 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S1200 = 0.3333333432674408f + nd_1.x;
    float _S1201 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_4 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1201, 1.0f, 0.0f, - _S1200, - _S1201, _S1200, 0.0f), T_1);
    float3  r0_1 = make_float3 (M_4.rows[int(0)].x, M_4.rows[int(0)].y, M_4.rows[int(0)].z);
    float3  r1_1 = make_float3 (M_4.rows[int(1)].x, M_4.rows[int(1)].y, M_4.rows[int(1)].z);
    float3  r2_12 = make_float3 (M_4.rows[int(2)].x, M_4.rows[int(2)].y, M_4.rows[int(2)].z);
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
    float _S1202 = _S1199.x;
    float _S1203 = _S1199.y;
    float intensity_1 = _S1202 + _S1203 + _S1199.z;
    float3  rgi_out_2 = mul_0(H_3, make_float3 (_S1202, _S1203, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / (rgi_out_2.z + 0.00000999999974738f));
    float _S1204 = rgi_out_3.x;
    float _S1205 = rgi_out_3.y;
    float3  _S1206 = clamp_1(make_float3 (_S1204, _S1205, rgi_out_3.z - _S1204 - _S1205), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S1207 = _S1206.x;
    float g0_1 = (F32_exp((_S1195.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S1195.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S1195.crf_params_0[int(0)].gc_0)));
    float y_11;
    if(_S1207 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S1207 / x0_1;
        float _S1208 = 1.0f - t0_0;
        y_11 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S1208) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S1208);
    }
    else
    {
        float _S1209 = 1.0f - y0_1;
        float _S1210 = 1.0f - x0_1;
        float s1_4 = _S1209 / _S1210;
        float t1_0 = (_S1207 - x0_1) / _S1210;
        float _S1211 = 1.0f - t1_0;
        y_11 = y0_1 + _S1209 * (s1_4 * t1_0 * t1_0 + gc_1 * t1_0 * _S1211) / (s1_4 + (gc_1 + g1_1 - 2.0f * s1_4) * t1_0 * _S1211);
    }
    *&((&rgb_out_3)->x) = y_11;
    float _S1212 = _S1206.y;
    float g0_2 = (F32_exp((_S1195.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S1195.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S1195.crf_params_0[int(1)].gc_0)));
    if(_S1212 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S1212 / x0_2;
        float _S1213 = 1.0f - t0_1;
        y_11 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S1213) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S1213);
    }
    else
    {
        float _S1214 = 1.0f - y0_2;
        float _S1215 = 1.0f - x0_2;
        float s1_5 = _S1214 / _S1215;
        float t1_1 = (_S1212 - x0_2) / _S1215;
        float _S1216 = 1.0f - t1_1;
        y_11 = y0_2 + _S1214 * (s1_5 * t1_1 * t1_1 + gc_2 * t1_1 * _S1216) / (s1_5 + (gc_2 + g1_2 - 2.0f * s1_5) * t1_1 * _S1216);
    }
    *&((&rgb_out_3)->y) = y_11;
    float _S1217 = _S1206.z;
    float g0_3 = (F32_exp((_S1195.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S1195.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S1195.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S1195.crf_params_0[int(2)].gc_0)));
    if(_S1217 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S1217 / x0_3;
        float _S1218 = 1.0f - t0_2;
        y_11 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S1218) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S1218);
    }
    else
    {
        float _S1219 = 1.0f - y0_3;
        float _S1220 = 1.0f - x0_3;
        float s1_6 = _S1219 / _S1220;
        float t1_2 = (_S1217 - x0_3) / _S1220;
        float _S1221 = 1.0f - t1_2;
        y_11 = y0_3 + _S1219 * (s1_6 * t1_2 * t1_2 + gc_3 * t1_2 * _S1221) / (s1_6 + (gc_3 + g1_3 - 2.0f * s1_6) * t1_2 * _S1221);
    }
    *&((&rgb_out_3)->z) = y_11;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1222)
{
    return (F32_exp2((_S1222)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1223, float2  _S1224)
{
    return mul_1(_S1223, _S1224);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1225, Matrix<float, 3, 3>  _S1226)
{
    return mul_3(_S1225, _S1226);
}

inline __device__ float s_primal_ctx_abs_0(float _S1227)
{
    return (F32_abs((_S1227)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1228, float3  _S1229)
{
    return mul_0(_S1228, _S1229);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1230, float3  _S1231, float3  _S1232)
{
    return clamp_1(_S1230, _S1231, _S1232);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1233, float _S1234, float _S1235)
{
    return lerp_0(_S1233, _S1234, _S1235);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1236, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1237, float3  _S1238)
{
    _d_mul_0(_S1236, _S1237, _S1238);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1239, float _S1240)
{
    _d_abs_0(_S1239, _S1240);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1241, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1242, Matrix<float, 3, 3>  _S1243)
{
    mul_2(_S1241, _S1242, _S1243);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1244, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1245, float2  _S1246)
{
    _d_mul_1(_S1244, _S1245, _S1246);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1247, float _S1248)
{
    _d_exp2_0(_S1247, _S1248);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1249 = *dprgb_in_0;
    float3  _S1250 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1251 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1252 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1253 = {
        _S1252, _S1252, _S1252
    };
    float2  _S1254 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1255 = { _S1254, _S1254, _S1254, _S1254 };
    CRFPPISPChannelParams_0 _S1256 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1257 = {
        _S1256, _S1256, _S1256
    };
    PPISPParams_0 _S1258;
    (&_S1258)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S1258)->vignette_params_1 = _S1253;
    (&_S1258)->color_params_1 = _S1255;
    (&_S1258)->crf_params_1 = _S1257;
    (&(&_S1258)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1258)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1259 = dpparams_0->primal_0[int(3)];
    (&(&_S1258)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1260 = dpparams_0->primal_0[int(4)];
    (&(&_S1258)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1261 = dpparams_0->primal_0[int(5)];
    (&(&_S1258)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1258)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1258)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1262 = dpparams_0->primal_0[int(8)];
    (&(&_S1258)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1263 = dpparams_0->primal_0[int(9)];
    (&(&_S1258)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1264 = dpparams_0->primal_0[int(10)];
    (&(&_S1258)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1258)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1258)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1265 = dpparams_0->primal_0[int(13)];
    (&(&_S1258)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1266 = dpparams_0->primal_0[int(14)];
    (&(&_S1258)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1267 = dpparams_0->primal_0[int(15)];
    (&(&_S1258)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1258)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1258)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1258)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1258)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1258)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1258)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1258)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1258)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1268 = dpparams_0->primal_0[int(24)];
    (&(&_S1258)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1269 = dpparams_0->primal_0[int(25)];
    (&(&_S1258)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1270 = dpparams_0->primal_0[int(26)];
    (&(&_S1258)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1271 = dpparams_0->primal_0[int(27)];
    (&(&_S1258)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1272 = dpparams_0->primal_0[int(28)];
    (&(&_S1258)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1273 = dpparams_0->primal_0[int(29)];
    (&(&_S1258)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1274 = dpparams_0->primal_0[int(30)];
    (&(&_S1258)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1275 = dpparams_0->primal_0[int(31)];
    (&(&_S1258)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1276 = dpparams_0->primal_0[int(32)];
    (&(&_S1258)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1277 = dpparams_0->primal_0[int(33)];
    (&(&_S1258)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1278 = dpparams_0->primal_0[int(34)];
    (&(&_S1258)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1279 = dpparams_0->primal_0[int(35)];
    (&(&_S1258)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1280 = _S1258;
    float _S1281 = s_primal_ctx_exp2_0(_S1258.exposure_1);
    float3  _S1282 = make_float3 (_S1281);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S1281);
    float _S1283 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S1284 = (pix_coord_2.x - image_center_2.x) / _S1283;
    float _S1285 = (pix_coord_2.y - image_center_2.y) / _S1283;
    float dx_8 = _S1284 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S1285 - dpparams_0->primal_0[int(2)];
    float r2_13 = dx_8 * dx_8 + dy_6 * dy_6;
    float r4_6 = r2_13 * r2_13;
    float r6_0 = r4_6 * r2_13;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_13 + 1.0f;
    float _S1286 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1287 = rgb_out_4.x * _S1286;
    float3  _S1288 = rgb_out_4;
    *&((&_S1288)->x) = _S1287;
    float dx_9 = _S1284 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S1285 - dpparams_0->primal_0[int(7)];
    float r2_14 = dx_9 * dx_9 + dy_7 * dy_7;
    float r4_7 = r2_14 * r2_14;
    float r6_1 = r4_7 * r2_14;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_14 + 1.0f;
    float _S1289 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1288)->y) = rgb_out_4.y * _S1289;
    float dx_10 = _S1284 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S1285 - dpparams_0->primal_0[int(12)];
    float r2_15 = dx_10 * dx_10 + dy_8 * dy_8;
    float r4_8 = r2_15 * r2_15;
    float r6_2 = r4_8 * r2_15;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_15 + 1.0f;
    float _S1290 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1288)->z) = rgb_out_4.z * _S1290;
    PPISPParams_0 _S1291 = _S1258;
    float2  _S1292 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1258.color_params_1.b_0);
    float2  _S1293 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1258.color_params_1.r_0);
    float2  _S1294 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1258.color_params_1.g_0);
    float2  _S1295 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1258.color_params_1.n_0);
    float _S1296 = 0.3333333432674408f + _S1295.x;
    float _S1297 = 0.3333333432674408f + _S1295.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S1292.x, 1.0f + _S1293.x, _S1294.x, _S1292.y, _S1293.y, 1.0f + _S1294.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1297, 1.0f, 0.0f, - _S1296, - _S1297, _S1296, 0.0f);
    Matrix<float, 3, 3>  _S1298 = s_primal_ctx_mul_1(skew_0, T_2);
    float3  r0_2 = make_float3 (_S1298.rows[int(0)].x, _S1298.rows[int(0)].y, _S1298.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S1298.rows[int(1)].x, _S1298.rows[int(1)].y, _S1298.rows[int(1)].z);
    float3  r2_16 = make_float3 (_S1298.rows[int(2)].x, _S1298.rows[int(2)].y, _S1298.rows[int(2)].z);
    float3  _S1299 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S1300 = (s_primal_ctx_dot_0(_S1299, _S1299)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S1301;
    bool _S1302;
    if(_S1300)
    {
        float3  _S1303 = s_primal_ctx_cross_0(r0_2, r2_16);
        bool _S1304 = (s_primal_ctx_dot_0(_S1303, _S1303)) < 9.99999968265522539e-21f;
        if(_S1304)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_16);
        }
        else
        {
            lambda_v_6 = _S1303;
        }
        _S1302 = _S1304;
        _S1301 = _S1303;
    }
    else
    {
        lambda_v_6 = _S1299;
        _S1302 = false;
        _S1301 = _S1250;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S1305 = s_primal_ctx_mul_1(T_2, D_0);
    Matrix<float, 3, 3>  _S1306 = s_primal_ctx_mul_1(_S1305, S_inv_0);
    bool _S1307 = (s_primal_ctx_abs_0(_S1306.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S1308;
    float _S1309;
    if(_S1307)
    {
        float inv_s_0 = 1.0f / _S1306.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1310 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1311 = _S1306.rows[int(2)].z * _S1306.rows[int(2)].z;
        H_4 = _S1306 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1308 = _S1310;
        _S1309 = _S1311;
    }
    else
    {
        H_4 = _S1306;
        _S1308 = _S1251;
        _S1309 = 0.0f;
    }
    float _S1312 = _S1288.x;
    float _S1313 = _S1288.y;
    float intensity_2 = _S1312 + _S1313 + _S1288.z;
    float3  rgi_in_0 = make_float3 (_S1312, _S1313, intensity_2);
    float3  _S1314 = s_primal_ctx_mul_2(H_4, rgi_in_0);
    float _S1315 = _S1314.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S1315;
    float3  _S1316 = make_float3 (norm_factor_0);
    float _S1317 = _S1315 * _S1315;
    float3  rgi_out_4 = _S1314 * make_float3 (norm_factor_0);
    float _S1318 = rgi_out_4.x;
    float _S1319 = rgi_out_4.y;
    float3  _S1320 = make_float3 (_S1318, _S1319, rgi_out_4.z - _S1318 - _S1319);
    float3  _S1321 = make_float3 (0.0f);
    float3  _S1322 = make_float3 (1.0f);
    float3  _S1323 = s_primal_ctx_clamp_1(_S1320, _S1321, _S1322);
    float _S1324 = _S1323.x;
    float _S1325 = 1.0f + s_primal_ctx_exp_0(_S1268);
    float _S1326 = 0.30000001192092896f + s_primal_ctx_log_0(_S1325);
    float _S1327 = 1.0f + s_primal_ctx_exp_0(_S1269);
    float _S1328 = 0.30000001192092896f + s_primal_ctx_log_0(_S1327);
    float _S1329 = 1.0f + s_primal_ctx_exp_0(_S1270);
    float _S1330 = 0.10000000149011612f + s_primal_ctx_log_0(_S1329);
    float _S1331 = - _S1271;
    float _S1332 = 1.0f + s_primal_ctx_exp_0(_S1331);
    float _S1333 = 1.0f / _S1332;
    float _S1334 = _S1332 * _S1332;
    float _S1335 = s_primal_ctx_lerp_0(_S1326, _S1328, _S1333);
    float _S1336 = _S1328 * _S1333;
    float a_4 = _S1336 / _S1335;
    float _S1337 = _S1335 * _S1335;
    float b_5 = 1.0f - a_4;
    bool _S1338 = _S1324 <= _S1333;
    float y_12;
    float _S1339;
    float _S1340;
    float _S1341;
    float _S1342;
    float _S1343;
    float _S1344;
    float _S1345;
    float _S1346;
    if(_S1338)
    {
        float _S1347 = _S1324 / _S1333;
        float _S1348 = _S1333 * _S1333;
        float _S1349 = s_primal_ctx_pow_0(_S1347, _S1326);
        y_12 = a_4 * _S1349;
        _S1339 = _S1349;
        _S1340 = _S1347;
        _S1341 = _S1348;
        _S1342 = 0.0f;
        _S1343 = 0.0f;
        _S1344 = 0.0f;
        _S1345 = 0.0f;
        _S1346 = 0.0f;
    }
    else
    {
        float _S1350 = 1.0f - _S1324;
        float _S1351 = 1.0f - _S1333;
        float _S1352 = _S1350 / _S1351;
        float _S1353 = _S1351 * _S1351;
        float _S1354 = s_primal_ctx_pow_0(_S1352, _S1328);
        y_12 = 1.0f - b_5 * _S1354;
        _S1339 = 0.0f;
        _S1340 = 0.0f;
        _S1341 = 0.0f;
        _S1342 = _S1354;
        _S1343 = _S1352;
        _S1344 = _S1353;
        _S1345 = _S1350;
        _S1346 = _S1351;
    }
    float _S1355 = (F32_max((0.0f), (y_12)));
    float _S1356 = _S1323.y;
    float _S1357 = 1.0f + s_primal_ctx_exp_0(_S1272);
    float _S1358 = 0.30000001192092896f + s_primal_ctx_log_0(_S1357);
    float _S1359 = 1.0f + s_primal_ctx_exp_0(_S1273);
    float _S1360 = 0.30000001192092896f + s_primal_ctx_log_0(_S1359);
    float _S1361 = 1.0f + s_primal_ctx_exp_0(_S1274);
    float _S1362 = 0.10000000149011612f + s_primal_ctx_log_0(_S1361);
    float _S1363 = - _S1275;
    float _S1364 = 1.0f + s_primal_ctx_exp_0(_S1363);
    float _S1365 = 1.0f / _S1364;
    float _S1366 = _S1364 * _S1364;
    float _S1367 = s_primal_ctx_lerp_0(_S1358, _S1360, _S1365);
    float _S1368 = _S1360 * _S1365;
    float a_5 = _S1368 / _S1367;
    float _S1369 = _S1367 * _S1367;
    float b_6 = 1.0f - a_5;
    bool _S1370 = _S1356 <= _S1365;
    float y_13;
    float _S1371;
    float _S1372;
    float _S1373;
    float _S1374;
    float _S1375;
    float _S1376;
    float _S1377;
    float _S1378;
    if(_S1370)
    {
        float _S1379 = _S1356 / _S1365;
        float _S1380 = _S1365 * _S1365;
        float _S1381 = s_primal_ctx_pow_0(_S1379, _S1358);
        y_13 = a_5 * _S1381;
        _S1371 = _S1381;
        _S1372 = _S1379;
        _S1373 = _S1380;
        _S1374 = 0.0f;
        _S1375 = 0.0f;
        _S1376 = 0.0f;
        _S1377 = 0.0f;
        _S1378 = 0.0f;
    }
    else
    {
        float _S1382 = 1.0f - _S1356;
        float _S1383 = 1.0f - _S1365;
        float _S1384 = _S1382 / _S1383;
        float _S1385 = _S1383 * _S1383;
        float _S1386 = s_primal_ctx_pow_0(_S1384, _S1360);
        y_13 = 1.0f - b_6 * _S1386;
        _S1371 = 0.0f;
        _S1372 = 0.0f;
        _S1373 = 0.0f;
        _S1374 = _S1386;
        _S1375 = _S1384;
        _S1376 = _S1385;
        _S1377 = _S1382;
        _S1378 = _S1383;
    }
    float _S1387 = (F32_max((0.0f), (y_13)));
    float _S1388 = _S1323.z;
    float _S1389 = 1.0f + s_primal_ctx_exp_0(_S1276);
    float _S1390 = 0.30000001192092896f + s_primal_ctx_log_0(_S1389);
    float _S1391 = 1.0f + s_primal_ctx_exp_0(_S1277);
    float _S1392 = 0.30000001192092896f + s_primal_ctx_log_0(_S1391);
    float _S1393 = 1.0f + s_primal_ctx_exp_0(_S1278);
    float _S1394 = 0.10000000149011612f + s_primal_ctx_log_0(_S1393);
    float _S1395 = - _S1279;
    float _S1396 = 1.0f + s_primal_ctx_exp_0(_S1395);
    float _S1397 = 1.0f / _S1396;
    float _S1398 = _S1396 * _S1396;
    float _S1399 = s_primal_ctx_lerp_0(_S1390, _S1392, _S1397);
    float _S1400 = _S1392 * _S1397;
    float a_6 = _S1400 / _S1399;
    float _S1401 = _S1399 * _S1399;
    float b_7 = 1.0f - a_6;
    bool _S1402 = _S1388 <= _S1397;
    float y_14;
    float _S1403;
    float _S1404;
    float _S1405;
    float _S1406;
    float _S1407;
    float _S1408;
    float _S1409;
    float _S1410;
    if(_S1402)
    {
        float _S1411 = _S1388 / _S1397;
        float _S1412 = _S1397 * _S1397;
        float _S1413 = s_primal_ctx_pow_0(_S1411, _S1390);
        y_14 = a_6 * _S1413;
        _S1403 = _S1413;
        _S1404 = _S1411;
        _S1405 = _S1412;
        _S1406 = 0.0f;
        _S1407 = 0.0f;
        _S1408 = 0.0f;
        _S1409 = 0.0f;
        _S1410 = 0.0f;
    }
    else
    {
        float _S1414 = 1.0f - _S1388;
        float _S1415 = 1.0f - _S1397;
        float _S1416 = _S1414 / _S1415;
        float _S1417 = _S1415 * _S1415;
        float _S1418 = s_primal_ctx_pow_0(_S1416, _S1392);
        y_14 = 1.0f - b_7 * _S1418;
        _S1403 = 0.0f;
        _S1404 = 0.0f;
        _S1405 = 0.0f;
        _S1406 = _S1418;
        _S1407 = _S1416;
        _S1408 = _S1417;
        _S1409 = _S1414;
        _S1410 = _S1415;
    }
    float _S1419 = (F32_max((0.0f), (y_14)));
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = _S1419;
    (&_S1420)->differential_0 = 0.0f;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1394;
    (&_S1421)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1420, &_S1421, _s_dOut_11.z);
    DiffPair_float_0 _S1422 = _S1421;
    DiffPair_float_0 _S1423;
    (&_S1423)->primal_0 = 0.0f;
    (&_S1423)->differential_0 = 0.0f;
    DiffPair_float_0 _S1424;
    (&_S1424)->primal_0 = y_14;
    (&_S1424)->differential_0 = 0.0f;
    _d_max_0(&_S1423, &_S1424, _S1420.differential_0);
    DiffPair_float_0 _S1425 = _S1424;
    if(_S1402)
    {
        float _S1426 = a_6 * _S1425.differential_0;
        float _S1427 = _S1403 * _S1425.differential_0;
        DiffPair_float_0 _S1428;
        (&_S1428)->primal_0 = _S1404;
        (&_S1428)->differential_0 = 0.0f;
        DiffPair_float_0 _S1429;
        (&_S1429)->primal_0 = _S1390;
        (&_S1429)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1428, &_S1429, _S1426);
        float _S1430 = _S1428.differential_0 / _S1405;
        float _S1431 = _S1388 * - _S1430;
        float _S1432 = _S1397 * _S1430;
        y_14 = 0.0f;
        _S1403 = _S1427;
        _S1404 = _S1431;
        _S1405 = 0.0f;
        _S1406 = _S1429.differential_0;
        _S1407 = _S1432;
    }
    else
    {
        float _S1433 = - _S1425.differential_0;
        float _S1434 = b_7 * _S1433;
        float _S1435 = _S1406 * _S1433;
        DiffPair_float_0 _S1436;
        (&_S1436)->primal_0 = _S1407;
        (&_S1436)->differential_0 = 0.0f;
        DiffPair_float_0 _S1437;
        (&_S1437)->primal_0 = _S1392;
        (&_S1437)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1436, &_S1437, _S1434);
        float _S1438 = _S1436.differential_0 / _S1408;
        float _S1439 = - (_S1409 * - _S1438);
        float _S1440 = - (_S1410 * _S1438);
        y_14 = _S1435;
        _S1403 = 0.0f;
        _S1404 = _S1439;
        _S1405 = _S1437.differential_0;
        _S1406 = 0.0f;
        _S1407 = _S1440;
    }
    float _S1441 = (- y_14 + _S1403) / _S1401;
    float _S1442 = _S1400 * - _S1441;
    float _S1443 = _S1399 * _S1441;
    float _S1444 = _S1392 * _S1443;
    float _S1445 = _S1397 * _S1443;
    DiffPair_float_0 _S1446;
    (&_S1446)->primal_0 = _S1390;
    (&_S1446)->differential_0 = 0.0f;
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1392;
    (&_S1447)->differential_0 = 0.0f;
    DiffPair_float_0 _S1448;
    (&_S1448)->primal_0 = _S1397;
    (&_S1448)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1446, &_S1447, &_S1448, _S1442);
    float _S1449 = - ((_S1444 + _S1448.differential_0 + _S1404) / _S1398);
    DiffPair_float_0 _S1450;
    (&_S1450)->primal_0 = _S1395;
    (&_S1450)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1450, _S1449);
    float _S1451 = - _S1450.differential_0;
    DiffPair_float_0 _S1452;
    (&_S1452)->primal_0 = _S1393;
    (&_S1452)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1452, _S1422.differential_0);
    DiffPair_float_0 _S1453;
    (&_S1453)->primal_0 = _S1278;
    (&_S1453)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1453, _S1452.differential_0);
    DiffPair_float_0 _S1454 = _S1453;
    float _S1455 = _S1445 + _S1447.differential_0 + _S1405;
    DiffPair_float_0 _S1456;
    (&_S1456)->primal_0 = _S1391;
    (&_S1456)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1456, _S1455);
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1277;
    (&_S1457)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1457, _S1456.differential_0);
    DiffPair_float_0 _S1458 = _S1457;
    float _S1459 = _S1446.differential_0 + _S1406;
    DiffPair_float_0 _S1460;
    (&_S1460)->primal_0 = _S1389;
    (&_S1460)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1460, _S1459);
    DiffPair_float_0 _S1461;
    (&_S1461)->primal_0 = _S1276;
    (&_S1461)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1461, _S1460.differential_0);
    DiffPair_float_0 _S1462 = _S1461;
    float3  _S1463 = make_float3 (0.0f, 0.0f, _S1407);
    DiffPair_float_0 _S1464;
    (&_S1464)->primal_0 = _S1387;
    (&_S1464)->differential_0 = 0.0f;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = _S1362;
    (&_S1465)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1464, &_S1465, _s_dOut_11.y);
    DiffPair_float_0 _S1466 = _S1465;
    DiffPair_float_0 _S1467;
    (&_S1467)->primal_0 = 0.0f;
    (&_S1467)->differential_0 = 0.0f;
    DiffPair_float_0 _S1468;
    (&_S1468)->primal_0 = y_13;
    (&_S1468)->differential_0 = 0.0f;
    _d_max_0(&_S1467, &_S1468, _S1464.differential_0);
    DiffPair_float_0 _S1469 = _S1468;
    if(_S1370)
    {
        float _S1470 = a_5 * _S1469.differential_0;
        float _S1471 = _S1371 * _S1469.differential_0;
        DiffPair_float_0 _S1472;
        (&_S1472)->primal_0 = _S1372;
        (&_S1472)->differential_0 = 0.0f;
        DiffPair_float_0 _S1473;
        (&_S1473)->primal_0 = _S1358;
        (&_S1473)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1472, &_S1473, _S1470);
        float _S1474 = _S1472.differential_0 / _S1373;
        float _S1475 = _S1356 * - _S1474;
        float _S1476 = _S1365 * _S1474;
        y_13 = 0.0f;
        _S1371 = _S1471;
        _S1372 = _S1475;
        _S1373 = 0.0f;
        _S1374 = _S1473.differential_0;
        _S1375 = _S1476;
    }
    else
    {
        float _S1477 = - _S1469.differential_0;
        float _S1478 = b_6 * _S1477;
        float _S1479 = _S1374 * _S1477;
        DiffPair_float_0 _S1480;
        (&_S1480)->primal_0 = _S1375;
        (&_S1480)->differential_0 = 0.0f;
        DiffPair_float_0 _S1481;
        (&_S1481)->primal_0 = _S1360;
        (&_S1481)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1480, &_S1481, _S1478);
        float _S1482 = _S1480.differential_0 / _S1376;
        float _S1483 = - (_S1377 * - _S1482);
        float _S1484 = - (_S1378 * _S1482);
        y_13 = _S1479;
        _S1371 = 0.0f;
        _S1372 = _S1483;
        _S1373 = _S1481.differential_0;
        _S1374 = 0.0f;
        _S1375 = _S1484;
    }
    float _S1485 = (- y_13 + _S1371) / _S1369;
    float _S1486 = _S1368 * - _S1485;
    float _S1487 = _S1367 * _S1485;
    float _S1488 = _S1360 * _S1487;
    float _S1489 = _S1365 * _S1487;
    DiffPair_float_0 _S1490;
    (&_S1490)->primal_0 = _S1358;
    (&_S1490)->differential_0 = 0.0f;
    DiffPair_float_0 _S1491;
    (&_S1491)->primal_0 = _S1360;
    (&_S1491)->differential_0 = 0.0f;
    DiffPair_float_0 _S1492;
    (&_S1492)->primal_0 = _S1365;
    (&_S1492)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1490, &_S1491, &_S1492, _S1486);
    float _S1493 = - ((_S1488 + _S1492.differential_0 + _S1372) / _S1366);
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = _S1363;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1494, _S1493);
    float _S1495 = - _S1494.differential_0;
    DiffPair_float_0 _S1496;
    (&_S1496)->primal_0 = _S1361;
    (&_S1496)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1496, _S1466.differential_0);
    DiffPair_float_0 _S1497;
    (&_S1497)->primal_0 = _S1274;
    (&_S1497)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1497, _S1496.differential_0);
    DiffPair_float_0 _S1498 = _S1497;
    float _S1499 = _S1489 + _S1491.differential_0 + _S1373;
    DiffPair_float_0 _S1500;
    (&_S1500)->primal_0 = _S1359;
    (&_S1500)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1500, _S1499);
    DiffPair_float_0 _S1501;
    (&_S1501)->primal_0 = _S1273;
    (&_S1501)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1501, _S1500.differential_0);
    DiffPair_float_0 _S1502 = _S1501;
    float _S1503 = _S1490.differential_0 + _S1374;
    DiffPair_float_0 _S1504;
    (&_S1504)->primal_0 = _S1357;
    (&_S1504)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1504, _S1503);
    DiffPair_float_0 _S1505;
    (&_S1505)->primal_0 = _S1272;
    (&_S1505)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1505, _S1504.differential_0);
    DiffPair_float_0 _S1506 = _S1505;
    float3  _S1507 = _S1463 + make_float3 (0.0f, _S1375, 0.0f);
    DiffPair_float_0 _S1508;
    (&_S1508)->primal_0 = _S1355;
    (&_S1508)->differential_0 = 0.0f;
    DiffPair_float_0 _S1509;
    (&_S1509)->primal_0 = _S1330;
    (&_S1509)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1508, &_S1509, _s_dOut_11.x);
    DiffPair_float_0 _S1510 = _S1509;
    DiffPair_float_0 _S1511;
    (&_S1511)->primal_0 = 0.0f;
    (&_S1511)->differential_0 = 0.0f;
    DiffPair_float_0 _S1512;
    (&_S1512)->primal_0 = y_12;
    (&_S1512)->differential_0 = 0.0f;
    _d_max_0(&_S1511, &_S1512, _S1508.differential_0);
    DiffPair_float_0 _S1513 = _S1512;
    if(_S1338)
    {
        float _S1514 = a_4 * _S1513.differential_0;
        float _S1515 = _S1339 * _S1513.differential_0;
        DiffPair_float_0 _S1516;
        (&_S1516)->primal_0 = _S1340;
        (&_S1516)->differential_0 = 0.0f;
        DiffPair_float_0 _S1517;
        (&_S1517)->primal_0 = _S1326;
        (&_S1517)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1516, &_S1517, _S1514);
        float _S1518 = _S1516.differential_0 / _S1341;
        float _S1519 = _S1324 * - _S1518;
        float _S1520 = _S1333 * _S1518;
        y_12 = 0.0f;
        _S1339 = _S1515;
        _S1340 = _S1519;
        _S1341 = 0.0f;
        _S1342 = _S1517.differential_0;
        _S1343 = _S1520;
    }
    else
    {
        float _S1521 = - _S1513.differential_0;
        float _S1522 = b_5 * _S1521;
        float _S1523 = _S1342 * _S1521;
        DiffPair_float_0 _S1524;
        (&_S1524)->primal_0 = _S1343;
        (&_S1524)->differential_0 = 0.0f;
        DiffPair_float_0 _S1525;
        (&_S1525)->primal_0 = _S1328;
        (&_S1525)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1524, &_S1525, _S1522);
        float _S1526 = _S1524.differential_0 / _S1344;
        float _S1527 = - (_S1345 * - _S1526);
        float _S1528 = - (_S1346 * _S1526);
        y_12 = _S1523;
        _S1339 = 0.0f;
        _S1340 = _S1527;
        _S1341 = _S1525.differential_0;
        _S1342 = 0.0f;
        _S1343 = _S1528;
    }
    float _S1529 = (- y_12 + _S1339) / _S1337;
    float _S1530 = _S1336 * - _S1529;
    float _S1531 = _S1335 * _S1529;
    float _S1532 = _S1328 * _S1531;
    float _S1533 = _S1333 * _S1531;
    DiffPair_float_0 _S1534;
    (&_S1534)->primal_0 = _S1326;
    (&_S1534)->differential_0 = 0.0f;
    DiffPair_float_0 _S1535;
    (&_S1535)->primal_0 = _S1328;
    (&_S1535)->differential_0 = 0.0f;
    DiffPair_float_0 _S1536;
    (&_S1536)->primal_0 = _S1333;
    (&_S1536)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1534, &_S1535, &_S1536, _S1530);
    float _S1537 = - ((_S1532 + _S1536.differential_0 + _S1340) / _S1334);
    DiffPair_float_0 _S1538;
    (&_S1538)->primal_0 = _S1331;
    (&_S1538)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1538, _S1537);
    float _S1539 = - _S1538.differential_0;
    DiffPair_float_0 _S1540;
    (&_S1540)->primal_0 = _S1329;
    (&_S1540)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1540, _S1510.differential_0);
    DiffPair_float_0 _S1541;
    (&_S1541)->primal_0 = _S1270;
    (&_S1541)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1541, _S1540.differential_0);
    DiffPair_float_0 _S1542 = _S1541;
    float _S1543 = _S1533 + _S1535.differential_0 + _S1341;
    DiffPair_float_0 _S1544;
    (&_S1544)->primal_0 = _S1327;
    (&_S1544)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1544, _S1543);
    DiffPair_float_0 _S1545;
    (&_S1545)->primal_0 = _S1269;
    (&_S1545)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1545, _S1544.differential_0);
    DiffPair_float_0 _S1546 = _S1545;
    float _S1547 = _S1534.differential_0 + _S1342;
    DiffPair_float_0 _S1548;
    (&_S1548)->primal_0 = _S1325;
    (&_S1548)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1548, _S1547);
    DiffPair_float_0 _S1549;
    (&_S1549)->primal_0 = _S1268;
    (&_S1549)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1549, _S1548.differential_0);
    DiffPair_float_0 _S1550 = _S1549;
    float3  _S1551 = _S1507 + make_float3 (_S1343, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1552;
    (&_S1552)->primal_0 = _S1320;
    (&_S1552)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1553;
    (&_S1553)->primal_0 = _S1321;
    (&_S1553)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1554;
    (&_S1554)->primal_0 = _S1322;
    (&_S1554)->differential_0 = _S1250;
    s_bwd_prop_clamp_1(&_S1552, &_S1553, &_S1554, _S1551);
    float _S1555 = - _S1552.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1552.differential_0.x + _S1555, _S1552.differential_0.y + _S1555, _S1552.differential_0.z);
    float3  _S1556 = _S1314 * s_diff_rgi_out_T_0;
    float _S1557 = (_S1556.x + _S1556.y + _S1556.z) / _S1317;
    float _S1558 = _S1315 * _S1557;
    float3  _S1559 = _S1316 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_2 * - _S1557);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1560;
    (&_S1560)->primal_0 = H_4;
    (&_S1560)->differential_0 = _S1251;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1561;
    (&_S1561)->primal_0 = rgi_in_0;
    (&_S1561)->differential_0 = _S1250;
    s_bwd_prop_mul_1(&_S1560, &_S1561, _S1559);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1562 = _S1560;
    float _S1563 = _S1558 + _S1561.differential_0.z;
    float _S1564 = _S1561.differential_0.y + _S1563;
    float _S1565 = _S1561.differential_0.x + _S1563;
    float3  _S1566 = make_float3 (_S1565, _S1564, _S1563);
    if(_S1307)
    {
        Matrix<float, 3, 3>  _S1567 = _S1306 * _S1562.differential_0;
        Matrix<float, 3, 3>  _S1568 = _S1308 * _S1562.differential_0;
        _S1309 = - ((_S1567.rows[int(0)].x + _S1567.rows[int(0)].y + _S1567.rows[int(0)].z + _S1567.rows[int(1)].x + _S1567.rows[int(1)].y + _S1567.rows[int(1)].z + _S1567.rows[int(2)].x + _S1567.rows[int(2)].y + _S1567.rows[int(2)].z) / _S1309);
        H_4 = _S1568;
    }
    else
    {
        _S1309 = 0.0f;
        H_4 = _S1562.differential_0;
    }
    DiffPair_float_0 _S1569;
    (&_S1569)->primal_0 = _S1306.rows[int(2)].z;
    (&_S1569)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1569, 0.0f);
    float _S1570 = _S1569.differential_0 + _S1309;
    float3  _S1571 = _S1250;
    *&((&_S1571)->z) = _S1570;
    Matrix<float, 3, 3>  _S1572 = _S1251;
    _S1572[int(2)] = _S1571;
    Matrix<float, 3, 3>  _S1573 = H_4 + _S1572;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1574;
    (&_S1574)->primal_0 = _S1305;
    (&_S1574)->differential_0 = _S1251;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1575;
    (&_S1575)->primal_0 = S_inv_0;
    (&_S1575)->differential_0 = _S1251;
    s_bwd_prop_mul_2(&_S1574, &_S1575, _S1573);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1576;
    (&_S1576)->primal_0 = T_2;
    (&_S1576)->differential_0 = _S1251;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1577;
    (&_S1577)->primal_0 = D_0;
    (&_S1577)->differential_0 = _S1251;
    s_bwd_prop_mul_2(&_S1576, &_S1577, _S1574.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1578 = _S1576;
    float3  _S1579 = make_float3 (_S1577.differential_0.rows[int(0)].x, _S1577.differential_0.rows[int(1)].y, _S1577.differential_0.rows[int(2)].z);
    float3  _S1580;
    if(_S1300)
    {
        if(_S1302)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1581;
            (&_S1581)->primal_0 = r1_2;
            (&_S1581)->differential_0 = _S1250;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1582;
            (&_S1582)->primal_0 = r2_16;
            (&_S1582)->differential_0 = _S1250;
            s_bwd_prop_cross_0(&_S1581, &_S1582, _S1579);
            _S1288 = _S1250;
            lambda_v_6 = _S1582.differential_0;
            _S1580 = _S1581.differential_0;
        }
        else
        {
            _S1288 = _S1579;
            lambda_v_6 = _S1250;
            _S1580 = _S1250;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1583;
        (&_S1583)->primal_0 = _S1301;
        (&_S1583)->differential_0 = _S1250;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1584;
        (&_S1584)->primal_0 = _S1301;
        (&_S1584)->differential_0 = _S1250;
        s_bwd_prop_dot_0(&_S1583, &_S1584, 0.0f);
        float3  _S1585 = _S1584.differential_0 + _S1583.differential_0 + _S1288;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1586;
        (&_S1586)->primal_0 = r0_2;
        (&_S1586)->differential_0 = _S1250;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1587;
        (&_S1587)->primal_0 = r2_16;
        (&_S1587)->differential_0 = _S1250;
        s_bwd_prop_cross_0(&_S1586, &_S1587, _S1585);
        float3  _S1588 = _S1587.differential_0 + lambda_v_6;
        _S1288 = _S1250;
        lambda_v_6 = _S1588;
        _S1301 = _S1580;
        _S1580 = _S1586.differential_0;
    }
    else
    {
        _S1288 = _S1579;
        lambda_v_6 = _S1250;
        _S1301 = _S1250;
        _S1580 = _S1250;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1589;
    (&_S1589)->primal_0 = _S1299;
    (&_S1589)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1590;
    (&_S1590)->primal_0 = _S1299;
    (&_S1590)->differential_0 = _S1250;
    s_bwd_prop_dot_0(&_S1589, &_S1590, 0.0f);
    float3  _S1591 = _S1590.differential_0 + _S1589.differential_0 + _S1288;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1592;
    (&_S1592)->primal_0 = r0_2;
    (&_S1592)->differential_0 = _S1250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1593;
    (&_S1593)->primal_0 = r1_2;
    (&_S1593)->differential_0 = _S1250;
    s_bwd_prop_cross_0(&_S1592, &_S1593, _S1591);
    float3  _S1594 = _S1250;
    *&((&_S1594)->z) = lambda_v_6.z;
    *&((&_S1594)->y) = lambda_v_6.y;
    *&((&_S1594)->x) = lambda_v_6.x;
    float3  _S1595 = _S1593.differential_0 + _S1301;
    float3  _S1596 = _S1250;
    *&((&_S1596)->z) = _S1595.z;
    *&((&_S1596)->y) = _S1595.y;
    *&((&_S1596)->x) = _S1595.x;
    float3  _S1597 = _S1592.differential_0 + _S1580;
    float3  _S1598 = _S1250;
    *&((&_S1598)->z) = _S1597.z;
    *&((&_S1598)->y) = _S1597.y;
    *&((&_S1598)->x) = _S1597.x;
    Matrix<float, 3, 3>  _S1599 = _S1251;
    _S1599[int(2)] = _S1594;
    _S1599[int(1)] = _S1596;
    _S1599[int(0)] = _S1598;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1600;
    (&_S1600)->primal_0 = skew_0;
    (&_S1600)->differential_0 = _S1251;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1601;
    (&_S1601)->primal_0 = T_2;
    (&_S1601)->differential_0 = _S1251;
    s_bwd_prop_mul_2(&_S1600, &_S1601, _S1599);
    Matrix<float, 3, 3>  _S1602 = _S1601.differential_0 + _S1578.differential_0;
    float2  _S1603 = make_float2 (_S1600.differential_0.rows[int(2)].y + - _S1600.differential_0.rows[int(1)].z, _S1600.differential_0.rows[int(0)].z + - _S1600.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1604 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1605;
    (&_S1605)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1605)->differential_0 = _S1604;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1606;
    (&_S1606)->primal_0 = _S1291.color_params_1.n_0;
    (&_S1606)->differential_0 = _S1254;
    s_bwd_prop_mul_3(&_S1605, &_S1606, _S1603);
    float2  _S1607 = make_float2 (_S1602.rows[int(0)].z, _S1602.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1608;
    (&_S1608)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1608)->differential_0 = _S1604;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1609;
    (&_S1609)->primal_0 = _S1291.color_params_1.g_0;
    (&_S1609)->differential_0 = _S1254;
    s_bwd_prop_mul_3(&_S1608, &_S1609, _S1607);
    float2  _S1610 = make_float2 (_S1602.rows[int(0)].y, _S1602.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1611;
    (&_S1611)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1611)->differential_0 = _S1604;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1612;
    (&_S1612)->primal_0 = _S1291.color_params_1.r_0;
    (&_S1612)->differential_0 = _S1254;
    s_bwd_prop_mul_3(&_S1611, &_S1612, _S1610);
    float2  _S1613 = make_float2 (_S1602.rows[int(0)].x, _S1602.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1614;
    (&_S1614)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1614)->differential_0 = _S1604;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1615;
    (&_S1615)->primal_0 = _S1291.color_params_1.b_0;
    (&_S1615)->differential_0 = _S1254;
    s_bwd_prop_mul_3(&_S1614, &_S1615, _S1613);
    ColorPPISPParams_0 _S1616 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1616)->n_0 = _S1606.differential_0;
    (&_S1616)->g_0 = _S1609.differential_0;
    (&_S1616)->r_0 = _S1612.differential_0;
    (&_S1616)->b_0 = _S1615.differential_0;
    _S1288 = _S1566;
    *&((&_S1288)->z) = 0.0f;
    float _S1617 = rgb_out_4.z * _S1563;
    float _S1618 = _S1290 * _S1563;
    DiffPair_float_0 _S1619;
    (&_S1619)->primal_0 = falloff_2;
    (&_S1619)->differential_0 = 0.0f;
    DiffPair_float_0 _S1620;
    (&_S1620)->primal_0 = 0.0f;
    (&_S1620)->differential_0 = 0.0f;
    DiffPair_float_0 _S1621;
    (&_S1621)->primal_0 = 1.0f;
    (&_S1621)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1619, &_S1620, &_S1621, _S1617);
    float _S1622 = r2_15 * _S1619.differential_0;
    float _S1623 = r4_8 * _S1619.differential_0;
    float s_diff_r6_T_0 = _S1267 * _S1619.differential_0;
    float _S1624 = r6_2 * _S1619.differential_0;
    float _S1625 = r2_15 * (_S1266 * _S1619.differential_0 + r2_15 * s_diff_r6_T_0);
    float _S1626 = _S1265 * _S1619.differential_0 + r4_8 * s_diff_r6_T_0 + _S1625 + _S1625;
    float _S1627 = dy_8 * _S1626;
    float _S1628 = dx_10 * _S1626;
    float _S1629 = - (_S1627 + _S1627);
    float _S1630 = - (_S1628 + _S1628);
    *&((&_S1288)->y) = 0.0f;
    float _S1631 = rgb_out_4.y * _S1564;
    float _S1632 = _S1289 * _S1564;
    DiffPair_float_0 _S1633;
    (&_S1633)->primal_0 = falloff_1;
    (&_S1633)->differential_0 = 0.0f;
    DiffPair_float_0 _S1634;
    (&_S1634)->primal_0 = 0.0f;
    (&_S1634)->differential_0 = 0.0f;
    DiffPair_float_0 _S1635;
    (&_S1635)->primal_0 = 1.0f;
    (&_S1635)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1633, &_S1634, &_S1635, _S1631);
    float _S1636 = r2_14 * _S1633.differential_0;
    float _S1637 = r4_7 * _S1633.differential_0;
    float s_diff_r6_T_1 = _S1264 * _S1633.differential_0;
    float _S1638 = r6_1 * _S1633.differential_0;
    float _S1639 = r2_14 * (_S1263 * _S1633.differential_0 + r2_14 * s_diff_r6_T_1);
    float _S1640 = _S1262 * _S1633.differential_0 + r4_7 * s_diff_r6_T_1 + _S1639 + _S1639;
    float _S1641 = dy_7 * _S1640;
    float _S1642 = dx_9 * _S1640;
    float _S1643 = - (_S1641 + _S1641);
    float _S1644 = - (_S1642 + _S1642);
    *&((&_S1288)->x) = 0.0f;
    float _S1645 = rgb_out_4.x * _S1565;
    float _S1646 = _S1286 * _S1565;
    DiffPair_float_0 _S1647;
    (&_S1647)->primal_0 = falloff_0;
    (&_S1647)->differential_0 = 0.0f;
    DiffPair_float_0 _S1648;
    (&_S1648)->primal_0 = 0.0f;
    (&_S1648)->differential_0 = 0.0f;
    DiffPair_float_0 _S1649;
    (&_S1649)->primal_0 = 1.0f;
    (&_S1649)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1647, &_S1648, &_S1649, _S1645);
    float _S1650 = r2_13 * _S1647.differential_0;
    float _S1651 = r4_6 * _S1647.differential_0;
    float s_diff_r6_T_2 = _S1261 * _S1647.differential_0;
    float _S1652 = r6_0 * _S1647.differential_0;
    float _S1653 = r2_13 * (_S1260 * _S1647.differential_0 + r2_13 * s_diff_r6_T_2);
    float _S1654 = _S1259 * _S1647.differential_0 + r4_6 * s_diff_r6_T_2 + _S1653 + _S1653;
    float _S1655 = dy_6 * _S1654;
    float _S1656 = dx_8 * _S1654;
    float _S1657 = - (_S1655 + _S1655);
    float _S1658 = - (_S1656 + _S1656);
    float3  _S1659 = _S1250;
    *&((&_S1659)->z) = _S1618;
    *&((&_S1659)->y) = _S1632;
    *&((&_S1659)->x) = _S1646;
    float3  _S1660 = _S1288 + _S1659;
    float3  _S1661 = _S1249.primal_0 * _S1660;
    float3  _S1662 = _S1282 * _S1660;
    float _S1663 = _S1661.x + _S1661.y + _S1661.z;
    DiffPair_float_0 _S1664;
    (&_S1664)->primal_0 = _S1280.exposure_1;
    (&_S1664)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1664, _S1663);
    PPISPParams_0 _S1665 = PPISPParams_x24_syn_dzero_0();
    (&_S1665)->color_params_1 = _S1616;
    (&_S1665)->exposure_1 = _S1664.differential_0;
    _S1258 = _S1665;
    (&(&_S1258)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1666 = _S1665.crf_params_1[int(2)].center_0 + _S1451;
    (&(&_S1258)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1667 = _S1665.crf_params_1[int(2)].gamma_0 + _S1454.differential_0;
    (&(&_S1258)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1668 = _S1665.crf_params_1[int(2)].shoulder_0 + _S1458.differential_0;
    (&(&_S1258)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1669 = _S1665.crf_params_1[int(2)].toe_0 + _S1462.differential_0;
    (&(&_S1258)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1670 = _S1665.crf_params_1[int(1)].center_0 + _S1495;
    (&(&_S1258)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1671 = _S1665.crf_params_1[int(1)].gamma_0 + _S1498.differential_0;
    (&(&_S1258)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1672 = _S1665.crf_params_1[int(1)].shoulder_0 + _S1502.differential_0;
    (&(&_S1258)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1673 = _S1665.crf_params_1[int(1)].toe_0 + _S1506.differential_0;
    (&(&_S1258)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1674 = _S1665.crf_params_1[int(0)].center_0 + _S1539;
    (&(&_S1258)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1675 = _S1665.crf_params_1[int(0)].gamma_0 + _S1542.differential_0;
    (&(&_S1258)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1676 = _S1665.crf_params_1[int(0)].shoulder_0 + _S1546.differential_0;
    (&(&_S1258)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1677 = _S1665.crf_params_1[int(0)].toe_0 + _S1550.differential_0;
    *&((&(&(&_S1258)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1258)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1258)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1678 = _S1624 + _S1665.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1258)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1679 = _S1623 + _S1665.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1258)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1680 = _S1622 + _S1665.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1258)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1681 = _S1629 + _S1665.vignette_params_1[int(2)].cy_0;
    (&(&_S1258)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1682 = _S1630 + _S1665.vignette_params_1[int(2)].cx_0;
    (&(&_S1258)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1683 = _S1638 + _S1665.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1258)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1684 = _S1637 + _S1665.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1258)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1685 = _S1636 + _S1665.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1258)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1686 = _S1643 + _S1665.vignette_params_1[int(1)].cy_0;
    (&(&_S1258)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1687 = _S1644 + _S1665.vignette_params_1[int(1)].cx_0;
    (&(&_S1258)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1688 = _S1652 + _S1665.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1258)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1689 = _S1651 + _S1665.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1258)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1690 = _S1650 + _S1665.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1258)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1691 = _S1657 + _S1665.vignette_params_1[int(0)].cy_0;
    (&(&_S1258)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1692 = _S1658 + _S1665.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S1693;
    _S1693[int(0)] = 0.0f;
    _S1693[int(1)] = 0.0f;
    _S1693[int(2)] = 0.0f;
    _S1693[int(3)] = 0.0f;
    _S1693[int(4)] = 0.0f;
    _S1693[int(5)] = 0.0f;
    _S1693[int(6)] = 0.0f;
    _S1693[int(7)] = 0.0f;
    _S1693[int(8)] = 0.0f;
    _S1693[int(9)] = 0.0f;
    _S1693[int(10)] = 0.0f;
    _S1693[int(11)] = 0.0f;
    _S1693[int(12)] = 0.0f;
    _S1693[int(13)] = 0.0f;
    _S1693[int(14)] = 0.0f;
    _S1693[int(15)] = 0.0f;
    _S1693[int(16)] = 0.0f;
    _S1693[int(17)] = 0.0f;
    _S1693[int(18)] = 0.0f;
    _S1693[int(19)] = 0.0f;
    _S1693[int(20)] = 0.0f;
    _S1693[int(21)] = 0.0f;
    _S1693[int(22)] = 0.0f;
    _S1693[int(23)] = 0.0f;
    _S1693[int(24)] = 0.0f;
    _S1693[int(25)] = 0.0f;
    _S1693[int(26)] = 0.0f;
    _S1693[int(27)] = 0.0f;
    _S1693[int(28)] = 0.0f;
    _S1693[int(29)] = 0.0f;
    _S1693[int(30)] = 0.0f;
    _S1693[int(31)] = 0.0f;
    _S1693[int(32)] = 0.0f;
    _S1693[int(33)] = 0.0f;
    _S1693[int(34)] = 0.0f;
    _S1693[int(35)] = 0.0f;
    _S1693[int(8)] = _S1685;
    _S1693[int(16)] = _S1665.color_params_1.b_0.x;
    _S1693[int(15)] = _S1678;
    _S1693[int(14)] = _S1679;
    _S1693[int(13)] = _S1680;
    _S1693[int(12)] = _S1681;
    _S1693[int(11)] = _S1682;
    _S1693[int(10)] = _S1683;
    _S1693[int(9)] = _S1684;
    _S1693[int(17)] = _S1665.color_params_1.b_0.y;
    _S1693[int(7)] = _S1686;
    _S1693[int(6)] = _S1687;
    _S1693[int(5)] = _S1688;
    _S1693[int(4)] = _S1689;
    _S1693[int(3)] = _S1690;
    _S1693[int(2)] = _S1691;
    _S1693[int(1)] = _S1692;
    _S1693[int(0)] = _S1258.exposure_1;
    _S1693[int(26)] = _S1675;
    _S1693[int(34)] = _S1667;
    _S1693[int(33)] = _S1668;
    _S1693[int(32)] = _S1669;
    _S1693[int(31)] = _S1670;
    _S1693[int(30)] = _S1671;
    _S1693[int(29)] = _S1672;
    _S1693[int(28)] = _S1673;
    _S1693[int(27)] = _S1674;
    _S1693[int(35)] = _S1666;
    _S1693[int(25)] = _S1676;
    _S1693[int(24)] = _S1677;
    _S1693[int(23)] = _S1665.color_params_1.n_0.y;
    _S1693[int(22)] = _S1665.color_params_1.n_0.x;
    _S1693[int(21)] = _S1665.color_params_1.g_0.y;
    _S1693[int(20)] = _S1665.color_params_1.g_0.x;
    _S1693[int(19)] = _S1665.color_params_1.r_0.y;
    _S1693[int(18)] = _S1665.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1693;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1662;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1694, float2  _S1695, float2  _S1696, float2  _S1697, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1698, float3  _S1699)
{
    s_bwd_prop_apply_ppisp_0(_S1694, _S1695, _S1696, _S1697, _S1698, _S1699);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1700 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S1700;
    FixedArray<float, 36>  _S1701 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_2;
    (&dp_params_0)->differential_0 = _S1701;
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
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1702 = *dprgb_in_1;
    float3  _S1703 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1704 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1705 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1706 = {
        _S1705, _S1705, _S1705
    };
    float2  _S1707 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1708 = { _S1707, _S1707, _S1707, _S1707 };
    RQSCRFPPISPChannelParams_0 _S1709 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1710 = {
        _S1709, _S1709, _S1709
    };
    PPISPParamsRQS_0 _S1711;
    (&_S1711)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1711)->vignette_params_0 = _S1706;
    (&_S1711)->color_params_0 = _S1708;
    (&_S1711)->crf_params_0 = _S1710;
    (&(&_S1711)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1711)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S1712 = dpparams_1->primal_0[int(3)];
    (&(&_S1711)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S1713 = dpparams_1->primal_0[int(4)];
    (&(&_S1711)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S1714 = dpparams_1->primal_0[int(5)];
    (&(&_S1711)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1711)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1711)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S1715 = dpparams_1->primal_0[int(8)];
    (&(&_S1711)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S1716 = dpparams_1->primal_0[int(9)];
    (&(&_S1711)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S1717 = dpparams_1->primal_0[int(10)];
    (&(&_S1711)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1711)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1711)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S1718 = dpparams_1->primal_0[int(13)];
    (&(&_S1711)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S1719 = dpparams_1->primal_0[int(14)];
    (&(&_S1711)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S1720 = dpparams_1->primal_0[int(15)];
    (&(&_S1711)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1711)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1711)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1711)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1711)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1711)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1711)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1711)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1711)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S1721 = dpparams_1->primal_0[int(24)];
    (&(&_S1711)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S1722 = dpparams_1->primal_0[int(25)];
    (&(&_S1711)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S1723 = dpparams_1->primal_0[int(26)];
    (&(&_S1711)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S1724 = dpparams_1->primal_0[int(27)];
    (&(&_S1711)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S1725 = dpparams_1->primal_0[int(28)];
    (&(&_S1711)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S1726 = dpparams_1->primal_0[int(29)];
    (&(&_S1711)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S1727 = dpparams_1->primal_0[int(30)];
    (&(&_S1711)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S1728 = dpparams_1->primal_0[int(31)];
    (&(&_S1711)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S1729 = dpparams_1->primal_0[int(32)];
    (&(&_S1711)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S1730 = dpparams_1->primal_0[int(33)];
    (&(&_S1711)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S1731 = dpparams_1->primal_0[int(34)];
    (&(&_S1711)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S1732 = dpparams_1->primal_0[int(35)];
    (&(&_S1711)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S1733 = dpparams_1->primal_0[int(36)];
    (&(&_S1711)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S1734 = dpparams_1->primal_0[int(37)];
    (&(&_S1711)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S1735 = dpparams_1->primal_0[int(38)];
    (&(&_S1711)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S1736 = _S1711;
    float _S1737 = s_primal_ctx_exp2_0(_S1711.exposure_0);
    float3  _S1738 = make_float3 (_S1737);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S1737);
    float _S1739 = (F32_max((img_size_4.x), (img_size_4.y)));
    float _S1740 = (pix_coord_4.x - image_center_4.x) / _S1739;
    float _S1741 = (pix_coord_4.y - image_center_4.y) / _S1739;
    float dx_11 = _S1740 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S1741 - dpparams_1->primal_0[int(2)];
    float r2_17 = dx_11 * dx_11 + dy_9 * dy_9;
    float r4_9 = r2_17 * r2_17;
    float r6_3 = r4_9 * r2_17;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_17 + 1.0f;
    float _S1742 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S1743 = rgb_out_5.x * _S1742;
    float3  _S1744 = rgb_out_5;
    *&((&_S1744)->x) = _S1743;
    float dx_12 = _S1740 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S1741 - dpparams_1->primal_0[int(7)];
    float r2_18 = dx_12 * dx_12 + dy_10 * dy_10;
    float r4_10 = r2_18 * r2_18;
    float r6_4 = r4_10 * r2_18;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_18 + 1.0f;
    float _S1745 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S1744)->y) = rgb_out_5.y * _S1745;
    float dx_13 = _S1740 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S1741 - dpparams_1->primal_0[int(12)];
    float r2_19 = dx_13 * dx_13 + dy_11 * dy_11;
    float r4_11 = r2_19 * r2_19;
    float r6_5 = r4_11 * r2_19;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_19 + 1.0f;
    float _S1746 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S1744)->z) = rgb_out_5.z * _S1746;
    PPISPParamsRQS_0 _S1747 = _S1711;
    float2  _S1748 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1711.color_params_0.b_0);
    float2  _S1749 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1711.color_params_0.r_0);
    float2  _S1750 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1711.color_params_0.g_0);
    float2  _S1751 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1711.color_params_0.n_0);
    float _S1752 = 0.3333333432674408f + _S1751.x;
    float _S1753 = 0.3333333432674408f + _S1751.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S1748.x, 1.0f + _S1749.x, _S1750.x, _S1748.y, _S1749.y, 1.0f + _S1750.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1753, 1.0f, 0.0f, - _S1752, - _S1753, _S1752, 0.0f);
    Matrix<float, 3, 3>  _S1754 = s_primal_ctx_mul_1(skew_1, T_3);
    float3  r0_3 = make_float3 (_S1754.rows[int(0)].x, _S1754.rows[int(0)].y, _S1754.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S1754.rows[int(1)].x, _S1754.rows[int(1)].y, _S1754.rows[int(1)].z);
    float3  r2_20 = make_float3 (_S1754.rows[int(2)].x, _S1754.rows[int(2)].y, _S1754.rows[int(2)].z);
    float3  _S1755 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S1756 = (s_primal_ctx_dot_0(_S1755, _S1755)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S1757;
    bool _S1758;
    if(_S1756)
    {
        float3  _S1759 = s_primal_ctx_cross_0(r0_3, r2_20);
        bool _S1760 = (s_primal_ctx_dot_0(_S1759, _S1759)) < 9.99999968265522539e-21f;
        if(_S1760)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_20);
        }
        else
        {
            lambda_v_7 = _S1759;
        }
        _S1758 = _S1760;
        _S1757 = _S1759;
    }
    else
    {
        lambda_v_7 = _S1755;
        _S1758 = false;
        _S1757 = _S1703;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S1761 = s_primal_ctx_mul_1(T_3, D_1);
    Matrix<float, 3, 3>  _S1762 = s_primal_ctx_mul_1(_S1761, S_inv_1);
    bool _S1763 = (s_primal_ctx_abs_0(_S1762.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S1764;
    float _S1765;
    if(_S1763)
    {
        float inv_s_1 = 1.0f / _S1762.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1766 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S1767 = _S1762.rows[int(2)].z * _S1762.rows[int(2)].z;
        H_5 = _S1762 * makeMatrix<float, 3, 3> (inv_s_1);
        _S1764 = _S1766;
        _S1765 = _S1767;
    }
    else
    {
        H_5 = _S1762;
        _S1764 = _S1704;
        _S1765 = 0.0f;
    }
    float _S1768 = _S1744.x;
    float _S1769 = _S1744.y;
    float intensity_3 = _S1768 + _S1769 + _S1744.z;
    float3  rgi_in_1 = make_float3 (_S1768, _S1769, intensity_3);
    float3  _S1770 = s_primal_ctx_mul_2(H_5, rgi_in_1);
    float _S1771 = _S1770.z + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S1771;
    float3  _S1772 = make_float3 (norm_factor_1);
    float _S1773 = _S1771 * _S1771;
    float3  rgi_out_5 = _S1770 * make_float3 (norm_factor_1);
    float _S1774 = rgi_out_5.x;
    float _S1775 = rgi_out_5.y;
    float3  _S1776 = make_float3 (_S1774, _S1775, rgi_out_5.z - _S1774 - _S1775);
    float3  _S1777 = make_float3 (0.0f);
    float3  _S1778 = make_float3 (1.0f);
    float3  _S1779 = s_primal_ctx_clamp_1(_S1776, _S1777, _S1778);
    float _S1780 = _S1779.x;
    float _S1781 = s_primal_ctx_exp_0(_S1721);
    float _S1782 = s_primal_ctx_exp_0(_S1722);
    float _S1783 = - _S1723;
    float _S1784 = 1.0f + s_primal_ctx_exp_0(_S1783);
    float x0_4 = 1.0f / _S1784;
    float _S1785 = _S1784 * _S1784;
    float _S1786 = - _S1724;
    float _S1787 = 1.0f + s_primal_ctx_exp_0(_S1786);
    float y0_4 = 1.0f / _S1787;
    float _S1788 = _S1787 * _S1787;
    float _S1789 = s_primal_ctx_exp_0(_S1725);
    bool _S1790 = _S1780 < x0_4;
    float _S1791;
    float _S1792;
    float _S1793;
    float _S1794;
    float _S1795;
    float _S1796;
    float _S1797;
    float _S1798;
    float _S1799;
    float _S1800;
    float _S1801;
    float _S1802;
    float _S1803;
    float _S1804;
    float _S1805;
    float _S1806;
    float _S1807;
    float _S1808;
    float _S1809;
    float _S1810;
    float _S1811;
    float _S1812;
    float _S1813;
    float _S1814;
    float _S1815;
    float _S1816;
    float _S1817;
    if(_S1790)
    {
        float s0_3 = y0_4 / x0_4;
        float _S1818 = x0_4 * x0_4;
        float t0_3 = _S1780 / x0_4;
        float _S1819 = s0_3 * t0_3;
        float _S1820 = _S1781 * t0_3;
        float _S1821 = 1.0f - t0_3;
        float _S1822 = _S1819 * t0_3 + _S1820 * _S1821;
        float _S1823 = y0_4 * _S1822;
        float _S1824 = _S1781 + _S1789 - 2.0f * s0_3;
        float _S1825 = _S1824 * t0_3;
        float _S1826 = s0_3 + _S1825 * _S1821;
        _S1791 = _S1826 * _S1826;
        _S1792 = _S1823;
        _S1793 = _S1826;
        _S1794 = _S1825;
        _S1795 = _S1821;
        _S1796 = _S1824;
        _S1797 = t0_3;
        _S1798 = _S1822;
        _S1799 = _S1820;
        _S1800 = _S1819;
        _S1801 = s0_3;
        _S1802 = _S1818;
        _S1803 = 0.0f;
        _S1804 = 0.0f;
        _S1805 = 0.0f;
        _S1806 = 0.0f;
        _S1807 = 0.0f;
        _S1808 = 0.0f;
        _S1809 = 0.0f;
        _S1810 = 0.0f;
        _S1811 = 0.0f;
        _S1812 = 0.0f;
        _S1813 = 0.0f;
        _S1814 = 0.0f;
        _S1815 = 0.0f;
        _S1816 = 0.0f;
        _S1817 = 0.0f;
    }
    else
    {
        float _S1827 = 1.0f - y0_4;
        float _S1828 = 1.0f - x0_4;
        float s1_7 = _S1827 / _S1828;
        float _S1829 = _S1828 * _S1828;
        float _S1830 = _S1780 - x0_4;
        float t1_3 = _S1830 / _S1828;
        float _S1831 = s1_7 * t1_3;
        float _S1832 = _S1789 * t1_3;
        float _S1833 = 1.0f - t1_3;
        float _S1834 = _S1831 * t1_3 + _S1832 * _S1833;
        float _S1835 = _S1827 * _S1834;
        float _S1836 = _S1789 + _S1782 - 2.0f * s1_7;
        float _S1837 = _S1836 * t1_3;
        float _S1838 = s1_7 + _S1837 * _S1833;
        float _S1839 = _S1838 * _S1838;
        _S1791 = 0.0f;
        _S1792 = 0.0f;
        _S1793 = 0.0f;
        _S1794 = 0.0f;
        _S1795 = 0.0f;
        _S1796 = 0.0f;
        _S1797 = 0.0f;
        _S1798 = 0.0f;
        _S1799 = 0.0f;
        _S1800 = 0.0f;
        _S1801 = 0.0f;
        _S1802 = 0.0f;
        _S1803 = _S1839;
        _S1804 = _S1835;
        _S1805 = _S1838;
        _S1806 = _S1837;
        _S1807 = _S1833;
        _S1808 = _S1836;
        _S1809 = t1_3;
        _S1810 = _S1827;
        _S1811 = _S1834;
        _S1812 = _S1832;
        _S1813 = _S1831;
        _S1814 = s1_7;
        _S1815 = _S1829;
        _S1816 = _S1830;
        _S1817 = _S1828;
    }
    float _S1840 = _S1779.y;
    float _S1841 = s_primal_ctx_exp_0(_S1726);
    float _S1842 = s_primal_ctx_exp_0(_S1727);
    float _S1843 = - _S1728;
    float _S1844 = 1.0f + s_primal_ctx_exp_0(_S1843);
    float x0_5 = 1.0f / _S1844;
    float _S1845 = _S1844 * _S1844;
    float _S1846 = - _S1729;
    float _S1847 = 1.0f + s_primal_ctx_exp_0(_S1846);
    float y0_5 = 1.0f / _S1847;
    float _S1848 = _S1847 * _S1847;
    float _S1849 = s_primal_ctx_exp_0(_S1730);
    bool _S1850 = _S1840 < x0_5;
    float _S1851;
    float _S1852;
    float _S1853;
    float _S1854;
    float _S1855;
    float _S1856;
    float _S1857;
    float _S1858;
    float _S1859;
    float _S1860;
    float _S1861;
    float _S1862;
    float _S1863;
    float _S1864;
    float _S1865;
    float _S1866;
    float _S1867;
    float _S1868;
    float _S1869;
    float _S1870;
    float _S1871;
    float _S1872;
    float _S1873;
    float _S1874;
    float _S1875;
    float _S1876;
    float _S1877;
    if(_S1850)
    {
        float s0_4 = y0_5 / x0_5;
        float _S1878 = x0_5 * x0_5;
        float t0_4 = _S1840 / x0_5;
        float _S1879 = s0_4 * t0_4;
        float _S1880 = _S1841 * t0_4;
        float _S1881 = 1.0f - t0_4;
        float _S1882 = _S1879 * t0_4 + _S1880 * _S1881;
        float _S1883 = y0_5 * _S1882;
        float _S1884 = _S1841 + _S1849 - 2.0f * s0_4;
        float _S1885 = _S1884 * t0_4;
        float _S1886 = s0_4 + _S1885 * _S1881;
        _S1851 = _S1886 * _S1886;
        _S1852 = _S1883;
        _S1853 = _S1886;
        _S1854 = _S1885;
        _S1855 = _S1881;
        _S1856 = _S1884;
        _S1857 = t0_4;
        _S1858 = _S1882;
        _S1859 = _S1880;
        _S1860 = _S1879;
        _S1861 = s0_4;
        _S1862 = _S1878;
        _S1863 = 0.0f;
        _S1864 = 0.0f;
        _S1865 = 0.0f;
        _S1866 = 0.0f;
        _S1867 = 0.0f;
        _S1868 = 0.0f;
        _S1869 = 0.0f;
        _S1870 = 0.0f;
        _S1871 = 0.0f;
        _S1872 = 0.0f;
        _S1873 = 0.0f;
        _S1874 = 0.0f;
        _S1875 = 0.0f;
        _S1876 = 0.0f;
        _S1877 = 0.0f;
    }
    else
    {
        float _S1887 = 1.0f - y0_5;
        float _S1888 = 1.0f - x0_5;
        float s1_8 = _S1887 / _S1888;
        float _S1889 = _S1888 * _S1888;
        float _S1890 = _S1840 - x0_5;
        float t1_4 = _S1890 / _S1888;
        float _S1891 = s1_8 * t1_4;
        float _S1892 = _S1849 * t1_4;
        float _S1893 = 1.0f - t1_4;
        float _S1894 = _S1891 * t1_4 + _S1892 * _S1893;
        float _S1895 = _S1887 * _S1894;
        float _S1896 = _S1849 + _S1842 - 2.0f * s1_8;
        float _S1897 = _S1896 * t1_4;
        float _S1898 = s1_8 + _S1897 * _S1893;
        float _S1899 = _S1898 * _S1898;
        _S1851 = 0.0f;
        _S1852 = 0.0f;
        _S1853 = 0.0f;
        _S1854 = 0.0f;
        _S1855 = 0.0f;
        _S1856 = 0.0f;
        _S1857 = 0.0f;
        _S1858 = 0.0f;
        _S1859 = 0.0f;
        _S1860 = 0.0f;
        _S1861 = 0.0f;
        _S1862 = 0.0f;
        _S1863 = _S1899;
        _S1864 = _S1895;
        _S1865 = _S1898;
        _S1866 = _S1897;
        _S1867 = _S1893;
        _S1868 = _S1896;
        _S1869 = t1_4;
        _S1870 = _S1887;
        _S1871 = _S1894;
        _S1872 = _S1892;
        _S1873 = _S1891;
        _S1874 = s1_8;
        _S1875 = _S1889;
        _S1876 = _S1890;
        _S1877 = _S1888;
    }
    float _S1900 = _S1779.z;
    float _S1901 = s_primal_ctx_exp_0(_S1731);
    float _S1902 = s_primal_ctx_exp_0(_S1732);
    float _S1903 = - _S1733;
    float _S1904 = 1.0f + s_primal_ctx_exp_0(_S1903);
    float x0_6 = 1.0f / _S1904;
    float _S1905 = _S1904 * _S1904;
    float _S1906 = - _S1734;
    float _S1907 = 1.0f + s_primal_ctx_exp_0(_S1906);
    float y0_6 = 1.0f / _S1907;
    float _S1908 = _S1907 * _S1907;
    float _S1909 = s_primal_ctx_exp_0(_S1735);
    bool _S1910 = _S1900 < x0_6;
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
    float _S1923;
    float _S1924;
    float _S1925;
    float _S1926;
    float _S1927;
    float _S1928;
    float _S1929;
    float _S1930;
    float _S1931;
    float _S1932;
    float _S1933;
    float _S1934;
    float _S1935;
    float _S1936;
    float _S1937;
    if(_S1910)
    {
        float s0_5 = y0_6 / x0_6;
        float _S1938 = x0_6 * x0_6;
        float t0_5 = _S1900 / x0_6;
        float _S1939 = s0_5 * t0_5;
        float _S1940 = _S1901 * t0_5;
        float _S1941 = 1.0f - t0_5;
        float _S1942 = _S1939 * t0_5 + _S1940 * _S1941;
        float _S1943 = y0_6 * _S1942;
        float _S1944 = _S1901 + _S1909 - 2.0f * s0_5;
        float _S1945 = _S1944 * t0_5;
        float _S1946 = s0_5 + _S1945 * _S1941;
        _S1911 = _S1946 * _S1946;
        _S1912 = _S1943;
        _S1913 = _S1946;
        _S1914 = _S1945;
        _S1915 = _S1941;
        _S1916 = _S1944;
        _S1917 = t0_5;
        _S1918 = _S1942;
        _S1919 = _S1940;
        _S1920 = _S1939;
        _S1921 = s0_5;
        _S1922 = _S1938;
        _S1923 = 0.0f;
        _S1924 = 0.0f;
        _S1925 = 0.0f;
        _S1926 = 0.0f;
        _S1927 = 0.0f;
        _S1928 = 0.0f;
        _S1929 = 0.0f;
        _S1930 = 0.0f;
        _S1931 = 0.0f;
        _S1932 = 0.0f;
        _S1933 = 0.0f;
        _S1934 = 0.0f;
        _S1935 = 0.0f;
        _S1936 = 0.0f;
        _S1937 = 0.0f;
    }
    else
    {
        float _S1947 = 1.0f - y0_6;
        float _S1948 = 1.0f - x0_6;
        float s1_9 = _S1947 / _S1948;
        float _S1949 = _S1948 * _S1948;
        float _S1950 = _S1900 - x0_6;
        float t1_5 = _S1950 / _S1948;
        float _S1951 = s1_9 * t1_5;
        float _S1952 = _S1909 * t1_5;
        float _S1953 = 1.0f - t1_5;
        float _S1954 = _S1951 * t1_5 + _S1952 * _S1953;
        float _S1955 = _S1947 * _S1954;
        float _S1956 = _S1909 + _S1902 - 2.0f * s1_9;
        float _S1957 = _S1956 * t1_5;
        float _S1958 = s1_9 + _S1957 * _S1953;
        float _S1959 = _S1958 * _S1958;
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
        _S1923 = _S1959;
        _S1924 = _S1955;
        _S1925 = _S1958;
        _S1926 = _S1957;
        _S1927 = _S1953;
        _S1928 = _S1956;
        _S1929 = t1_5;
        _S1930 = _S1947;
        _S1931 = _S1954;
        _S1932 = _S1952;
        _S1933 = _S1951;
        _S1934 = s1_9;
        _S1935 = _S1949;
        _S1936 = _S1950;
        _S1937 = _S1948;
    }
    if(_S1910)
    {
        float _S1960 = _s_dOut_12.z / _S1911;
        float _S1961 = _S1912 * - _S1960;
        float _S1962 = _S1913 * _S1960;
        float _S1963 = _S1915 * _S1961;
        float _S1964 = _S1917 * _S1963;
        float _S1965 = y0_6 * _S1962;
        float _S1966 = _S1915 * _S1965;
        float _S1967 = _S1917 * _S1965;
        float _S1968 = (_S1916 * _S1963 + - (_S1914 * _S1961 + _S1919 * _S1965) + _S1901 * _S1966 + _S1920 * _S1965 + _S1921 * _S1967) / _S1922;
        float _S1969 = x0_6 * _S1968;
        float _S1970 = (_S1961 + 2.0f * - _S1964 + _S1917 * _S1967) / _S1922;
        float _S1971 = _S1918 * _S1962 + x0_6 * _S1970;
        float _S1972 = _S1964 + _S1917 * _S1966;
        float _S1973 = _S1900 * - _S1968 + y0_6 * - _S1970;
        _S1911 = _S1964;
        _S1912 = _S1971;
        _S1913 = _S1973;
        _S1914 = 0.0f;
        _S1915 = _S1972;
        _S1916 = _S1969;
    }
    else
    {
        float _S1974 = _s_dOut_12.z / _S1923;
        float _S1975 = _S1924 * - _S1974;
        float _S1976 = _S1925 * _S1974;
        float _S1977 = _S1927 * _S1975;
        float _S1978 = _S1929 * _S1977;
        float _S1979 = _S1930 * _S1976;
        float _S1980 = _S1927 * _S1979;
        float _S1981 = _S1929 * _S1979;
        float _S1982 = (_S1928 * _S1977 + - (_S1926 * _S1975 + _S1932 * _S1979) + _S1909 * _S1980 + _S1933 * _S1979 + _S1934 * _S1981) / _S1935;
        float _S1983 = _S1937 * _S1982;
        float _S1984 = (_S1975 + 2.0f * - _S1978 + _S1929 * _S1981) / _S1935;
        float _S1985 = _s_dOut_12.z + - (_S1931 * _S1976 + _S1937 * _S1984);
        float _S1986 = - _S1983 + - (_S1936 * - _S1982 + _S1930 * - _S1984);
        _S1911 = _S1978 + _S1929 * _S1980;
        _S1912 = _S1985;
        _S1913 = _S1986;
        _S1914 = _S1978;
        _S1915 = 0.0f;
        _S1916 = _S1983;
    }
    DiffPair_float_0 _S1987;
    (&_S1987)->primal_0 = _S1735;
    (&_S1987)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1987, _S1911);
    DiffPair_float_0 _S1988 = _S1987;
    float _S1989 = - (_S1912 / _S1908);
    DiffPair_float_0 _S1990;
    (&_S1990)->primal_0 = _S1906;
    (&_S1990)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1990, _S1989);
    float _S1991 = - _S1990.differential_0;
    float _S1992 = - (_S1913 / _S1905);
    DiffPair_float_0 _S1993;
    (&_S1993)->primal_0 = _S1903;
    (&_S1993)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1993, _S1992);
    float _S1994 = - _S1993.differential_0;
    DiffPair_float_0 _S1995;
    (&_S1995)->primal_0 = _S1732;
    (&_S1995)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1995, _S1914);
    DiffPair_float_0 _S1996 = _S1995;
    DiffPair_float_0 _S1997;
    (&_S1997)->primal_0 = _S1731;
    (&_S1997)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1997, _S1915);
    DiffPair_float_0 _S1998 = _S1997;
    float3  _S1999 = make_float3 (0.0f, 0.0f, _S1916);
    if(_S1850)
    {
        float _S2000 = _s_dOut_12.y / _S1851;
        float _S2001 = _S1852 * - _S2000;
        float _S2002 = _S1853 * _S2000;
        float _S2003 = _S1855 * _S2001;
        float _S2004 = _S1857 * _S2003;
        float _S2005 = y0_5 * _S2002;
        float _S2006 = _S1855 * _S2005;
        float _S2007 = _S1857 * _S2005;
        float _S2008 = (_S1856 * _S2003 + - (_S1854 * _S2001 + _S1859 * _S2005) + _S1841 * _S2006 + _S1860 * _S2005 + _S1861 * _S2007) / _S1862;
        float _S2009 = x0_5 * _S2008;
        float _S2010 = (_S2001 + 2.0f * - _S2004 + _S1857 * _S2007) / _S1862;
        float _S2011 = _S1858 * _S2002 + x0_5 * _S2010;
        float _S2012 = _S2004 + _S1857 * _S2006;
        float _S2013 = _S1840 * - _S2008 + y0_5 * - _S2010;
        _S1851 = _S2004;
        _S1852 = _S2011;
        _S1853 = _S2013;
        _S1854 = 0.0f;
        _S1855 = _S2012;
        _S1856 = _S2009;
    }
    else
    {
        float _S2014 = _s_dOut_12.y / _S1863;
        float _S2015 = _S1864 * - _S2014;
        float _S2016 = _S1865 * _S2014;
        float _S2017 = _S1867 * _S2015;
        float _S2018 = _S1869 * _S2017;
        float _S2019 = _S1870 * _S2016;
        float _S2020 = _S1867 * _S2019;
        float _S2021 = _S1869 * _S2019;
        float _S2022 = (_S1868 * _S2017 + - (_S1866 * _S2015 + _S1872 * _S2019) + _S1849 * _S2020 + _S1873 * _S2019 + _S1874 * _S2021) / _S1875;
        float _S2023 = _S1877 * _S2022;
        float _S2024 = (_S2015 + 2.0f * - _S2018 + _S1869 * _S2021) / _S1875;
        float _S2025 = _s_dOut_12.y + - (_S1871 * _S2016 + _S1877 * _S2024);
        float _S2026 = - _S2023 + - (_S1876 * - _S2022 + _S1870 * - _S2024);
        _S1851 = _S2018 + _S1869 * _S2020;
        _S1852 = _S2025;
        _S1853 = _S2026;
        _S1854 = _S2018;
        _S1855 = 0.0f;
        _S1856 = _S2023;
    }
    DiffPair_float_0 _S2027;
    (&_S2027)->primal_0 = _S1730;
    (&_S2027)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2027, _S1851);
    DiffPair_float_0 _S2028 = _S2027;
    float _S2029 = - (_S1852 / _S1848);
    DiffPair_float_0 _S2030;
    (&_S2030)->primal_0 = _S1846;
    (&_S2030)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2030, _S2029);
    float _S2031 = - _S2030.differential_0;
    float _S2032 = - (_S1853 / _S1845);
    DiffPair_float_0 _S2033;
    (&_S2033)->primal_0 = _S1843;
    (&_S2033)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2033, _S2032);
    float _S2034 = - _S2033.differential_0;
    DiffPair_float_0 _S2035;
    (&_S2035)->primal_0 = _S1727;
    (&_S2035)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2035, _S1854);
    DiffPair_float_0 _S2036 = _S2035;
    DiffPair_float_0 _S2037;
    (&_S2037)->primal_0 = _S1726;
    (&_S2037)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2037, _S1855);
    DiffPair_float_0 _S2038 = _S2037;
    float3  _S2039 = _S1999 + make_float3 (0.0f, _S1856, 0.0f);
    if(_S1790)
    {
        float _S2040 = _s_dOut_12.x / _S1791;
        float _S2041 = _S1792 * - _S2040;
        float _S2042 = _S1793 * _S2040;
        float _S2043 = _S1795 * _S2041;
        float _S2044 = _S1797 * _S2043;
        float _S2045 = y0_4 * _S2042;
        float _S2046 = _S1795 * _S2045;
        float _S2047 = _S1797 * _S2045;
        float _S2048 = (_S1796 * _S2043 + - (_S1794 * _S2041 + _S1799 * _S2045) + _S1781 * _S2046 + _S1800 * _S2045 + _S1801 * _S2047) / _S1802;
        float _S2049 = x0_4 * _S2048;
        float _S2050 = (_S2041 + 2.0f * - _S2044 + _S1797 * _S2047) / _S1802;
        float _S2051 = _S1798 * _S2042 + x0_4 * _S2050;
        float _S2052 = _S2044 + _S1797 * _S2046;
        float _S2053 = _S1780 * - _S2048 + y0_4 * - _S2050;
        _S1791 = _S2044;
        _S1792 = _S2051;
        _S1793 = _S2053;
        _S1794 = 0.0f;
        _S1795 = _S2052;
        _S1796 = _S2049;
    }
    else
    {
        float _S2054 = _s_dOut_12.x / _S1803;
        float _S2055 = _S1804 * - _S2054;
        float _S2056 = _S1805 * _S2054;
        float _S2057 = _S1807 * _S2055;
        float _S2058 = _S1809 * _S2057;
        float _S2059 = _S1810 * _S2056;
        float _S2060 = _S1807 * _S2059;
        float _S2061 = _S1809 * _S2059;
        float _S2062 = (_S1808 * _S2057 + - (_S1806 * _S2055 + _S1812 * _S2059) + _S1789 * _S2060 + _S1813 * _S2059 + _S1814 * _S2061) / _S1815;
        float _S2063 = _S1817 * _S2062;
        float _S2064 = (_S2055 + 2.0f * - _S2058 + _S1809 * _S2061) / _S1815;
        float _S2065 = _s_dOut_12.x + - (_S1811 * _S2056 + _S1817 * _S2064);
        float _S2066 = - _S2063 + - (_S1816 * - _S2062 + _S1810 * - _S2064);
        _S1791 = _S2058 + _S1809 * _S2060;
        _S1792 = _S2065;
        _S1793 = _S2066;
        _S1794 = _S2058;
        _S1795 = 0.0f;
        _S1796 = _S2063;
    }
    DiffPair_float_0 _S2067;
    (&_S2067)->primal_0 = _S1725;
    (&_S2067)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2067, _S1791);
    DiffPair_float_0 _S2068 = _S2067;
    float _S2069 = - (_S1792 / _S1788);
    DiffPair_float_0 _S2070;
    (&_S2070)->primal_0 = _S1786;
    (&_S2070)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2070, _S2069);
    float _S2071 = - _S2070.differential_0;
    float _S2072 = - (_S1793 / _S1785);
    DiffPair_float_0 _S2073;
    (&_S2073)->primal_0 = _S1783;
    (&_S2073)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2073, _S2072);
    float _S2074 = - _S2073.differential_0;
    DiffPair_float_0 _S2075;
    (&_S2075)->primal_0 = _S1722;
    (&_S2075)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2075, _S1794);
    DiffPair_float_0 _S2076 = _S2075;
    DiffPair_float_0 _S2077;
    (&_S2077)->primal_0 = _S1721;
    (&_S2077)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2077, _S1795);
    DiffPair_float_0 _S2078 = _S2077;
    float3  _S2079 = _S2039 + make_float3 (_S1796, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2080;
    (&_S2080)->primal_0 = _S1776;
    (&_S2080)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2081;
    (&_S2081)->primal_0 = _S1777;
    (&_S2081)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2082;
    (&_S2082)->primal_0 = _S1778;
    (&_S2082)->differential_0 = _S1703;
    s_bwd_prop_clamp_1(&_S2080, &_S2081, &_S2082, _S2079);
    float _S2083 = - _S2080.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S2080.differential_0.x + _S2083, _S2080.differential_0.y + _S2083, _S2080.differential_0.z);
    float3  _S2084 = _S1770 * s_diff_rgi_out_T_1;
    float _S2085 = (_S2084.x + _S2084.y + _S2084.z) / _S1773;
    float _S2086 = _S1771 * _S2085;
    float3  _S2087 = _S1772 * s_diff_rgi_out_T_1 + make_float3 (0.0f, 0.0f, intensity_3 * - _S2085);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2088;
    (&_S2088)->primal_0 = H_5;
    (&_S2088)->differential_0 = _S1704;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2089;
    (&_S2089)->primal_0 = rgi_in_1;
    (&_S2089)->differential_0 = _S1703;
    s_bwd_prop_mul_1(&_S2088, &_S2089, _S2087);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2090 = _S2088;
    float _S2091 = _S2086 + _S2089.differential_0.z;
    float _S2092 = _S2089.differential_0.y + _S2091;
    float _S2093 = _S2089.differential_0.x + _S2091;
    float3  _S2094 = make_float3 (_S2093, _S2092, _S2091);
    if(_S1763)
    {
        Matrix<float, 3, 3>  _S2095 = _S1762 * _S2090.differential_0;
        Matrix<float, 3, 3>  _S2096 = _S1764 * _S2090.differential_0;
        _S1765 = - ((_S2095.rows[int(0)].x + _S2095.rows[int(0)].y + _S2095.rows[int(0)].z + _S2095.rows[int(1)].x + _S2095.rows[int(1)].y + _S2095.rows[int(1)].z + _S2095.rows[int(2)].x + _S2095.rows[int(2)].y + _S2095.rows[int(2)].z) / _S1765);
        H_5 = _S2096;
    }
    else
    {
        _S1765 = 0.0f;
        H_5 = _S2090.differential_0;
    }
    DiffPair_float_0 _S2097;
    (&_S2097)->primal_0 = _S1762.rows[int(2)].z;
    (&_S2097)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2097, 0.0f);
    float _S2098 = _S2097.differential_0 + _S1765;
    float3  _S2099 = _S1703;
    *&((&_S2099)->z) = _S2098;
    Matrix<float, 3, 3>  _S2100 = _S1704;
    _S2100[int(2)] = _S2099;
    Matrix<float, 3, 3>  _S2101 = H_5 + _S2100;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2102;
    (&_S2102)->primal_0 = _S1761;
    (&_S2102)->differential_0 = _S1704;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2103;
    (&_S2103)->primal_0 = S_inv_1;
    (&_S2103)->differential_0 = _S1704;
    s_bwd_prop_mul_2(&_S2102, &_S2103, _S2101);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2104;
    (&_S2104)->primal_0 = T_3;
    (&_S2104)->differential_0 = _S1704;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2105;
    (&_S2105)->primal_0 = D_1;
    (&_S2105)->differential_0 = _S1704;
    s_bwd_prop_mul_2(&_S2104, &_S2105, _S2102.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2106 = _S2104;
    float3  _S2107 = make_float3 (_S2105.differential_0.rows[int(0)].x, _S2105.differential_0.rows[int(1)].y, _S2105.differential_0.rows[int(2)].z);
    float3  _S2108;
    if(_S1756)
    {
        if(_S1758)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2109;
            (&_S2109)->primal_0 = r1_3;
            (&_S2109)->differential_0 = _S1703;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2110;
            (&_S2110)->primal_0 = r2_20;
            (&_S2110)->differential_0 = _S1703;
            s_bwd_prop_cross_0(&_S2109, &_S2110, _S2107);
            _S1744 = _S1703;
            lambda_v_7 = _S2110.differential_0;
            _S2108 = _S2109.differential_0;
        }
        else
        {
            _S1744 = _S2107;
            lambda_v_7 = _S1703;
            _S2108 = _S1703;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2111;
        (&_S2111)->primal_0 = _S1757;
        (&_S2111)->differential_0 = _S1703;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2112;
        (&_S2112)->primal_0 = _S1757;
        (&_S2112)->differential_0 = _S1703;
        s_bwd_prop_dot_0(&_S2111, &_S2112, 0.0f);
        float3  _S2113 = _S2112.differential_0 + _S2111.differential_0 + _S1744;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2114;
        (&_S2114)->primal_0 = r0_3;
        (&_S2114)->differential_0 = _S1703;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2115;
        (&_S2115)->primal_0 = r2_20;
        (&_S2115)->differential_0 = _S1703;
        s_bwd_prop_cross_0(&_S2114, &_S2115, _S2113);
        float3  _S2116 = _S2115.differential_0 + lambda_v_7;
        _S1744 = _S1703;
        lambda_v_7 = _S2116;
        _S1757 = _S2108;
        _S2108 = _S2114.differential_0;
    }
    else
    {
        _S1744 = _S2107;
        lambda_v_7 = _S1703;
        _S1757 = _S1703;
        _S2108 = _S1703;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2117;
    (&_S2117)->primal_0 = _S1755;
    (&_S2117)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2118;
    (&_S2118)->primal_0 = _S1755;
    (&_S2118)->differential_0 = _S1703;
    s_bwd_prop_dot_0(&_S2117, &_S2118, 0.0f);
    float3  _S2119 = _S2118.differential_0 + _S2117.differential_0 + _S1744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2120;
    (&_S2120)->primal_0 = r0_3;
    (&_S2120)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2121;
    (&_S2121)->primal_0 = r1_3;
    (&_S2121)->differential_0 = _S1703;
    s_bwd_prop_cross_0(&_S2120, &_S2121, _S2119);
    float3  _S2122 = _S1703;
    *&((&_S2122)->z) = lambda_v_7.z;
    *&((&_S2122)->y) = lambda_v_7.y;
    *&((&_S2122)->x) = lambda_v_7.x;
    float3  _S2123 = _S2121.differential_0 + _S1757;
    float3  _S2124 = _S1703;
    *&((&_S2124)->z) = _S2123.z;
    *&((&_S2124)->y) = _S2123.y;
    *&((&_S2124)->x) = _S2123.x;
    float3  _S2125 = _S2120.differential_0 + _S2108;
    float3  _S2126 = _S1703;
    *&((&_S2126)->z) = _S2125.z;
    *&((&_S2126)->y) = _S2125.y;
    *&((&_S2126)->x) = _S2125.x;
    Matrix<float, 3, 3>  _S2127 = _S1704;
    _S2127[int(2)] = _S2122;
    _S2127[int(1)] = _S2124;
    _S2127[int(0)] = _S2126;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2128;
    (&_S2128)->primal_0 = skew_1;
    (&_S2128)->differential_0 = _S1704;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2129;
    (&_S2129)->primal_0 = T_3;
    (&_S2129)->differential_0 = _S1704;
    s_bwd_prop_mul_2(&_S2128, &_S2129, _S2127);
    Matrix<float, 3, 3>  _S2130 = _S2129.differential_0 + _S2106.differential_0;
    float2  _S2131 = make_float2 (_S2128.differential_0.rows[int(2)].y + - _S2128.differential_0.rows[int(1)].z, _S2128.differential_0.rows[int(0)].z + - _S2128.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S2132 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2133;
    (&_S2133)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2133)->differential_0 = _S2132;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2134;
    (&_S2134)->primal_0 = _S1747.color_params_0.n_0;
    (&_S2134)->differential_0 = _S1707;
    s_bwd_prop_mul_3(&_S2133, &_S2134, _S2131);
    float2  _S2135 = make_float2 (_S2130.rows[int(0)].z, _S2130.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2136;
    (&_S2136)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2136)->differential_0 = _S2132;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2137;
    (&_S2137)->primal_0 = _S1747.color_params_0.g_0;
    (&_S2137)->differential_0 = _S1707;
    s_bwd_prop_mul_3(&_S2136, &_S2137, _S2135);
    float2  _S2138 = make_float2 (_S2130.rows[int(0)].y, _S2130.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2139;
    (&_S2139)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2139)->differential_0 = _S2132;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2140;
    (&_S2140)->primal_0 = _S1747.color_params_0.r_0;
    (&_S2140)->differential_0 = _S1707;
    s_bwd_prop_mul_3(&_S2139, &_S2140, _S2138);
    float2  _S2141 = make_float2 (_S2130.rows[int(0)].x, _S2130.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2142;
    (&_S2142)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2142)->differential_0 = _S2132;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2143;
    (&_S2143)->primal_0 = _S1747.color_params_0.b_0;
    (&_S2143)->differential_0 = _S1707;
    s_bwd_prop_mul_3(&_S2142, &_S2143, _S2141);
    ColorPPISPParams_0 _S2144 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2144)->n_0 = _S2134.differential_0;
    (&_S2144)->g_0 = _S2137.differential_0;
    (&_S2144)->r_0 = _S2140.differential_0;
    (&_S2144)->b_0 = _S2143.differential_0;
    _S1744 = _S2094;
    *&((&_S1744)->z) = 0.0f;
    float _S2145 = rgb_out_5.z * _S2091;
    float _S2146 = _S1746 * _S2091;
    DiffPair_float_0 _S2147;
    (&_S2147)->primal_0 = falloff_5;
    (&_S2147)->differential_0 = 0.0f;
    DiffPair_float_0 _S2148;
    (&_S2148)->primal_0 = 0.0f;
    (&_S2148)->differential_0 = 0.0f;
    DiffPair_float_0 _S2149;
    (&_S2149)->primal_0 = 1.0f;
    (&_S2149)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2147, &_S2148, &_S2149, _S2145);
    float _S2150 = r2_19 * _S2147.differential_0;
    float _S2151 = r4_11 * _S2147.differential_0;
    float s_diff_r6_T_3 = _S1720 * _S2147.differential_0;
    float _S2152 = r6_5 * _S2147.differential_0;
    float _S2153 = r2_19 * (_S1719 * _S2147.differential_0 + r2_19 * s_diff_r6_T_3);
    float _S2154 = _S1718 * _S2147.differential_0 + r4_11 * s_diff_r6_T_3 + _S2153 + _S2153;
    float _S2155 = dy_11 * _S2154;
    float _S2156 = dx_13 * _S2154;
    float _S2157 = - (_S2155 + _S2155);
    float _S2158 = - (_S2156 + _S2156);
    *&((&_S1744)->y) = 0.0f;
    float _S2159 = rgb_out_5.y * _S2092;
    float _S2160 = _S1745 * _S2092;
    DiffPair_float_0 _S2161;
    (&_S2161)->primal_0 = falloff_4;
    (&_S2161)->differential_0 = 0.0f;
    DiffPair_float_0 _S2162;
    (&_S2162)->primal_0 = 0.0f;
    (&_S2162)->differential_0 = 0.0f;
    DiffPair_float_0 _S2163;
    (&_S2163)->primal_0 = 1.0f;
    (&_S2163)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2161, &_S2162, &_S2163, _S2159);
    float _S2164 = r2_18 * _S2161.differential_0;
    float _S2165 = r4_10 * _S2161.differential_0;
    float s_diff_r6_T_4 = _S1717 * _S2161.differential_0;
    float _S2166 = r6_4 * _S2161.differential_0;
    float _S2167 = r2_18 * (_S1716 * _S2161.differential_0 + r2_18 * s_diff_r6_T_4);
    float _S2168 = _S1715 * _S2161.differential_0 + r4_10 * s_diff_r6_T_4 + _S2167 + _S2167;
    float _S2169 = dy_10 * _S2168;
    float _S2170 = dx_12 * _S2168;
    float _S2171 = - (_S2169 + _S2169);
    float _S2172 = - (_S2170 + _S2170);
    *&((&_S1744)->x) = 0.0f;
    float _S2173 = rgb_out_5.x * _S2093;
    float _S2174 = _S1742 * _S2093;
    DiffPair_float_0 _S2175;
    (&_S2175)->primal_0 = falloff_3;
    (&_S2175)->differential_0 = 0.0f;
    DiffPair_float_0 _S2176;
    (&_S2176)->primal_0 = 0.0f;
    (&_S2176)->differential_0 = 0.0f;
    DiffPair_float_0 _S2177;
    (&_S2177)->primal_0 = 1.0f;
    (&_S2177)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2175, &_S2176, &_S2177, _S2173);
    float _S2178 = r2_17 * _S2175.differential_0;
    float _S2179 = r4_9 * _S2175.differential_0;
    float s_diff_r6_T_5 = _S1714 * _S2175.differential_0;
    float _S2180 = r6_3 * _S2175.differential_0;
    float _S2181 = r2_17 * (_S1713 * _S2175.differential_0 + r2_17 * s_diff_r6_T_5);
    float _S2182 = _S1712 * _S2175.differential_0 + r4_9 * s_diff_r6_T_5 + _S2181 + _S2181;
    float _S2183 = dy_9 * _S2182;
    float _S2184 = dx_11 * _S2182;
    float _S2185 = - (_S2183 + _S2183);
    float _S2186 = - (_S2184 + _S2184);
    float3  _S2187 = _S1703;
    *&((&_S2187)->z) = _S2146;
    *&((&_S2187)->y) = _S2160;
    *&((&_S2187)->x) = _S2174;
    float3  _S2188 = _S1744 + _S2187;
    float3  _S2189 = _S1702.primal_0 * _S2188;
    float3  _S2190 = _S1738 * _S2188;
    float _S2191 = _S2189.x + _S2189.y + _S2189.z;
    DiffPair_float_0 _S2192;
    (&_S2192)->primal_0 = _S1736.exposure_0;
    (&_S2192)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2192, _S2191);
    PPISPParamsRQS_0 _S2193 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2193)->color_params_0 = _S2144;
    (&_S2193)->exposure_0 = _S2192.differential_0;
    _S1711 = _S2193;
    (&(&_S1711)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2194 = _S2193.crf_params_0[int(2)].gc_0 + _S1988.differential_0;
    (&(&_S1711)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2195 = _S2193.crf_params_0[int(2)].y0_0 + _S1991;
    (&(&_S1711)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2196 = _S2193.crf_params_0[int(2)].x0_0 + _S1994;
    (&(&_S1711)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2197 = _S2193.crf_params_0[int(2)].g1_0 + _S1996.differential_0;
    (&(&_S1711)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2198 = _S2193.crf_params_0[int(2)].g0_0 + _S1998.differential_0;
    (&(&_S1711)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2199 = _S2193.crf_params_0[int(1)].gc_0 + _S2028.differential_0;
    (&(&_S1711)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2200 = _S2193.crf_params_0[int(1)].y0_0 + _S2031;
    (&(&_S1711)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2201 = _S2193.crf_params_0[int(1)].x0_0 + _S2034;
    (&(&_S1711)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2202 = _S2193.crf_params_0[int(1)].g1_0 + _S2036.differential_0;
    (&(&_S1711)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2203 = _S2193.crf_params_0[int(1)].g0_0 + _S2038.differential_0;
    (&(&_S1711)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2204 = _S2193.crf_params_0[int(0)].gc_0 + _S2068.differential_0;
    (&(&_S1711)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2205 = _S2193.crf_params_0[int(0)].y0_0 + _S2071;
    (&(&_S1711)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2206 = _S2193.crf_params_0[int(0)].x0_0 + _S2074;
    (&(&_S1711)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2207 = _S2193.crf_params_0[int(0)].g1_0 + _S2076.differential_0;
    (&(&_S1711)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2208 = _S2193.crf_params_0[int(0)].g0_0 + _S2078.differential_0;
    *&((&(&(&_S1711)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1711)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1711)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2209 = _S2152 + _S2193.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1711)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2210 = _S2151 + _S2193.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1711)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2211 = _S2150 + _S2193.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1711)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2212 = _S2157 + _S2193.vignette_params_0[int(2)].cy_0;
    (&(&_S1711)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2213 = _S2158 + _S2193.vignette_params_0[int(2)].cx_0;
    (&(&_S1711)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2214 = _S2166 + _S2193.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1711)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2215 = _S2165 + _S2193.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1711)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2216 = _S2164 + _S2193.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1711)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2217 = _S2171 + _S2193.vignette_params_0[int(1)].cy_0;
    (&(&_S1711)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2218 = _S2172 + _S2193.vignette_params_0[int(1)].cx_0;
    (&(&_S1711)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2219 = _S2180 + _S2193.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1711)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2220 = _S2179 + _S2193.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1711)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2221 = _S2178 + _S2193.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1711)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2222 = _S2185 + _S2193.vignette_params_0[int(0)].cy_0;
    (&(&_S1711)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2223 = _S2186 + _S2193.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2224;
    _S2224[int(0)] = 0.0f;
    _S2224[int(1)] = 0.0f;
    _S2224[int(2)] = 0.0f;
    _S2224[int(3)] = 0.0f;
    _S2224[int(4)] = 0.0f;
    _S2224[int(5)] = 0.0f;
    _S2224[int(6)] = 0.0f;
    _S2224[int(7)] = 0.0f;
    _S2224[int(8)] = 0.0f;
    _S2224[int(9)] = 0.0f;
    _S2224[int(10)] = 0.0f;
    _S2224[int(11)] = 0.0f;
    _S2224[int(12)] = 0.0f;
    _S2224[int(13)] = 0.0f;
    _S2224[int(14)] = 0.0f;
    _S2224[int(15)] = 0.0f;
    _S2224[int(16)] = 0.0f;
    _S2224[int(17)] = 0.0f;
    _S2224[int(18)] = 0.0f;
    _S2224[int(19)] = 0.0f;
    _S2224[int(20)] = 0.0f;
    _S2224[int(21)] = 0.0f;
    _S2224[int(22)] = 0.0f;
    _S2224[int(23)] = 0.0f;
    _S2224[int(24)] = 0.0f;
    _S2224[int(25)] = 0.0f;
    _S2224[int(26)] = 0.0f;
    _S2224[int(27)] = 0.0f;
    _S2224[int(28)] = 0.0f;
    _S2224[int(29)] = 0.0f;
    _S2224[int(30)] = 0.0f;
    _S2224[int(31)] = 0.0f;
    _S2224[int(32)] = 0.0f;
    _S2224[int(33)] = 0.0f;
    _S2224[int(34)] = 0.0f;
    _S2224[int(35)] = 0.0f;
    _S2224[int(36)] = 0.0f;
    _S2224[int(37)] = 0.0f;
    _S2224[int(38)] = 0.0f;
    _S2224[int(9)] = _S2215;
    _S2224[int(18)] = _S2193.color_params_0.r_0.x;
    _S2224[int(17)] = _S2193.color_params_0.b_0.y;
    _S2224[int(16)] = _S2193.color_params_0.b_0.x;
    _S2224[int(15)] = _S2209;
    _S2224[int(14)] = _S2210;
    _S2224[int(13)] = _S2211;
    _S2224[int(12)] = _S2212;
    _S2224[int(11)] = _S2213;
    _S2224[int(10)] = _S2214;
    _S2224[int(19)] = _S2193.color_params_0.r_0.y;
    _S2224[int(8)] = _S2216;
    _S2224[int(7)] = _S2217;
    _S2224[int(6)] = _S2218;
    _S2224[int(5)] = _S2219;
    _S2224[int(4)] = _S2220;
    _S2224[int(3)] = _S2221;
    _S2224[int(2)] = _S2222;
    _S2224[int(1)] = _S2223;
    _S2224[int(0)] = _S1711.exposure_0;
    _S2224[int(28)] = _S2204;
    _S2224[int(37)] = _S2195;
    _S2224[int(36)] = _S2196;
    _S2224[int(35)] = _S2197;
    _S2224[int(34)] = _S2198;
    _S2224[int(33)] = _S2199;
    _S2224[int(32)] = _S2200;
    _S2224[int(31)] = _S2201;
    _S2224[int(30)] = _S2202;
    _S2224[int(29)] = _S2203;
    _S2224[int(38)] = _S2194;
    _S2224[int(27)] = _S2205;
    _S2224[int(26)] = _S2206;
    _S2224[int(25)] = _S2207;
    _S2224[int(24)] = _S2208;
    _S2224[int(23)] = _S2193.color_params_0.n_0.y;
    _S2224[int(22)] = _S2193.color_params_0.n_0.x;
    _S2224[int(21)] = _S2193.color_params_0.g_0.y;
    _S2224[int(20)] = _S2193.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S2224;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S2190;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2225, float2  _S2226, float2  _S2227, float2  _S2228, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2229, float3  _S2230)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S2225, _S2226, _S2227, _S2228, _S2229, _S2230);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S2231 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S2231;
    FixedArray<float, 39>  _S2232 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_3;
    (&dp_params_1)->differential_0 = _S2232;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_4, FixedArray<float, 22>  * _S2233)
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
    float _S2234 = p_2.vignette_params_1[int(0)].cx_0;
    float _S2235 = p_2.vignette_params_1[int(0)].cy_0;
    float _S2236 = p_2.vignette_params_1[int(1)].cx_0;
    float _S2237 = p_2.vignette_params_1[int(1)].cy_0;
    float _S2238 = p_2.vignette_params_1[int(2)].cx_0;
    float _S2239 = p_2.vignette_params_1[int(2)].cy_0;
    losses_3[int(1)] = _S2234 * _S2234 + _S2235 * _S2235 + _S2236 * _S2236 + _S2237 * _S2237 + _S2238 * _S2238 + _S2239 * _S2239;
    losses_3[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_3[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_3[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_3 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S2240 = p_2.vignette_params_1[int(0)].cx_0 - mean_3;
    float _S2241 = p_2.vignette_params_1[int(1)].cx_0 - mean_3;
    float _S2242 = p_2.vignette_params_1[int(2)].cx_0 - mean_3;
    losses_3[int(5)] = (_S2240 * _S2240 + _S2241 * _S2241 + _S2242 * _S2242) / 3.0f;
    float mean_4 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S2243 = p_2.vignette_params_1[int(0)].cy_0 - mean_4;
    float _S2244 = p_2.vignette_params_1[int(1)].cy_0 - mean_4;
    float _S2245 = p_2.vignette_params_1[int(2)].cy_0 - mean_4;
    losses_3[int(6)] = (_S2243 * _S2243 + _S2244 * _S2244 + _S2245 * _S2245) / 3.0f;
    float mean_5 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S2246 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_5;
    float _S2247 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_5;
    float _S2248 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_5;
    losses_3[int(7)] = (_S2246 * _S2246 + _S2247 * _S2247 + _S2248 * _S2248) / 3.0f;
    float mean_6 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S2249 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_6;
    float _S2250 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_6;
    float _S2251 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_6;
    losses_3[int(8)] = (_S2249 * _S2249 + _S2250 * _S2250 + _S2251 * _S2251) / 3.0f;
    float mean_7 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S2252 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_7;
    float _S2253 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_7;
    float _S2254 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_7;
    losses_3[int(9)] = (_S2252 * _S2252 + _S2253 * _S2253 + _S2254 * _S2254) / 3.0f;
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
    float _S2255 = p_2.crf_params_1[int(0)].toe_0 - mean_8;
    float _S2256 = p_2.crf_params_1[int(1)].toe_0 - mean_8;
    float _S2257 = p_2.crf_params_1[int(2)].toe_0 - mean_8;
    losses_3[int(18)] = (_S2255 * _S2255 + _S2256 * _S2256 + _S2257 * _S2257) / 3.0f;
    float mean_9 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S2258 = p_2.crf_params_1[int(0)].shoulder_0 - mean_9;
    float _S2259 = p_2.crf_params_1[int(1)].shoulder_0 - mean_9;
    float _S2260 = p_2.crf_params_1[int(2)].shoulder_0 - mean_9;
    losses_3[int(19)] = (_S2258 * _S2258 + _S2259 * _S2259 + _S2260 * _S2260) / 3.0f;
    float mean_10 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S2261 = p_2.crf_params_1[int(0)].gamma_0 - mean_10;
    float _S2262 = p_2.crf_params_1[int(1)].gamma_0 - mean_10;
    float _S2263 = p_2.crf_params_1[int(2)].gamma_0 - mean_10;
    losses_3[int(20)] = (_S2261 * _S2261 + _S2262 * _S2262 + _S2263 * _S2263) / 3.0f;
    float mean_11 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S2264 = p_2.crf_params_1[int(0)].center_0 - mean_11;
    float _S2265 = p_2.crf_params_1[int(1)].center_0 - mean_11;
    float _S2266 = p_2.crf_params_1[int(2)].center_0 - mean_11;
    losses_3[int(21)] = (_S2264 * _S2264 + _S2265 * _S2265 + _S2266 * _S2266) / 3.0f;
    *_S2233 = losses_3;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_5, FixedArray<float, 23>  * _S2267)
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
    float _S2268 = p_3.vignette_params_0[int(0)].cx_0;
    float _S2269 = p_3.vignette_params_0[int(0)].cy_0;
    float _S2270 = p_3.vignette_params_0[int(1)].cx_0;
    float _S2271 = p_3.vignette_params_0[int(1)].cy_0;
    float _S2272 = p_3.vignette_params_0[int(2)].cx_0;
    float _S2273 = p_3.vignette_params_0[int(2)].cy_0;
    losses_4[int(1)] = _S2268 * _S2268 + _S2269 * _S2269 + _S2270 * _S2270 + _S2271 * _S2271 + _S2272 * _S2272 + _S2273 * _S2273;
    losses_4[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_4[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_4[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_12 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S2274 = p_3.vignette_params_0[int(0)].cx_0 - mean_12;
    float _S2275 = p_3.vignette_params_0[int(1)].cx_0 - mean_12;
    float _S2276 = p_3.vignette_params_0[int(2)].cx_0 - mean_12;
    losses_4[int(5)] = (_S2274 * _S2274 + _S2275 * _S2275 + _S2276 * _S2276) / 3.0f;
    float mean_13 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S2277 = p_3.vignette_params_0[int(0)].cy_0 - mean_13;
    float _S2278 = p_3.vignette_params_0[int(1)].cy_0 - mean_13;
    float _S2279 = p_3.vignette_params_0[int(2)].cy_0 - mean_13;
    losses_4[int(6)] = (_S2277 * _S2277 + _S2278 * _S2278 + _S2279 * _S2279) / 3.0f;
    float mean_14 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S2280 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_14;
    float _S2281 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_14;
    float _S2282 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_14;
    losses_4[int(7)] = (_S2280 * _S2280 + _S2281 * _S2281 + _S2282 * _S2282) / 3.0f;
    float mean_15 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S2283 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_15;
    float _S2284 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_15;
    float _S2285 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_15;
    losses_4[int(8)] = (_S2283 * _S2283 + _S2284 * _S2284 + _S2285 * _S2285) / 3.0f;
    float mean_16 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S2286 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_16;
    float _S2287 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_16;
    float _S2288 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_16;
    losses_4[int(9)] = (_S2286 * _S2286 + _S2287 * _S2287 + _S2288 * _S2288) / 3.0f;
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
    float _S2289 = p_3.crf_params_0[int(0)].g0_0 - mean_17;
    float _S2290 = p_3.crf_params_0[int(1)].g0_0 - mean_17;
    float _S2291 = p_3.crf_params_0[int(2)].g0_0 - mean_17;
    losses_4[int(18)] = (_S2289 * _S2289 + _S2290 * _S2290 + _S2291 * _S2291) / 3.0f;
    float mean_18 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S2292 = p_3.crf_params_0[int(0)].g1_0 - mean_18;
    float _S2293 = p_3.crf_params_0[int(1)].g1_0 - mean_18;
    float _S2294 = p_3.crf_params_0[int(2)].g1_0 - mean_18;
    losses_4[int(19)] = (_S2292 * _S2292 + _S2293 * _S2293 + _S2294 * _S2294) / 3.0f;
    float mean_19 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S2295 = p_3.crf_params_0[int(0)].x0_0 - mean_19;
    float _S2296 = p_3.crf_params_0[int(1)].x0_0 - mean_19;
    float _S2297 = p_3.crf_params_0[int(2)].x0_0 - mean_19;
    losses_4[int(20)] = (_S2295 * _S2295 + _S2296 * _S2296 + _S2297 * _S2297) / 3.0f;
    float mean_20 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S2298 = p_3.crf_params_0[int(0)].y0_0 - mean_20;
    float _S2299 = p_3.crf_params_0[int(1)].y0_0 - mean_20;
    float _S2300 = p_3.crf_params_0[int(2)].y0_0 - mean_20;
    losses_4[int(21)] = (_S2298 * _S2298 + _S2299 * _S2299 + _S2300 * _S2300) / 3.0f;
    float mean_21 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S2301 = p_3.crf_params_0[int(0)].gc_0 - mean_21;
    float _S2302 = p_3.crf_params_0[int(1)].gc_0 - mean_21;
    float _S2303 = p_3.crf_params_0[int(2)].gc_0 - mean_21;
    losses_4[int(22)] = (_S2301 * _S2301 + _S2302 * _S2302 + _S2303 * _S2303) / 3.0f;
    *_S2267 = losses_4;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_13)
{
    VignettingChannelParams_0 _S2304 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2305 = {
        _S2304, _S2304, _S2304
    };
    float2  _S2306 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2307 = { _S2306, _S2306, _S2306, _S2306 };
    CRFPPISPChannelParams_0 _S2308 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S2309 = {
        _S2308, _S2308, _S2308
    };
    PPISPParams_0 _S2310;
    (&_S2310)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S2310)->vignette_params_1 = _S2305;
    (&_S2310)->color_params_1 = _S2307;
    (&_S2310)->crf_params_1 = _S2309;
    (&(&_S2310)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S2310)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S2310)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S2310)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S2310)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S2310)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S2310)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S2310)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S2310)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S2310)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S2310)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S2310)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S2310)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S2310)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S2310)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S2310)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S2310)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S2310)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S2310)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S2310)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S2310)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S2310)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S2310)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S2310)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S2310)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S2310)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S2310)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S2310)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S2310)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S2310)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S2310)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S2310)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S2310)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S2310)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S2310)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_22 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S2311 = dpparams_2->primal_0[int(1)] - mean_22;
    float _S2312 = dpparams_2->primal_0[int(6)] - mean_22;
    float _S2313 = dpparams_2->primal_0[int(11)] - mean_22;
    float mean_23 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S2314 = dpparams_2->primal_0[int(2)] - mean_23;
    float _S2315 = dpparams_2->primal_0[int(7)] - mean_23;
    float _S2316 = dpparams_2->primal_0[int(12)] - mean_23;
    float mean_24 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S2317 = dpparams_2->primal_0[int(3)] - mean_24;
    float _S2318 = dpparams_2->primal_0[int(8)] - mean_24;
    float _S2319 = dpparams_2->primal_0[int(13)] - mean_24;
    float mean_25 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S2320 = dpparams_2->primal_0[int(4)] - mean_25;
    float _S2321 = dpparams_2->primal_0[int(9)] - mean_25;
    float _S2322 = dpparams_2->primal_0[int(14)] - mean_25;
    float mean_26 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S2323 = dpparams_2->primal_0[int(5)] - mean_26;
    float _S2324 = dpparams_2->primal_0[int(10)] - mean_26;
    float _S2325 = dpparams_2->primal_0[int(15)] - mean_26;
    float mean_27 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_28 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_29 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_30 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S2326 = 0.3333333432674408f * (*_s_dOut_13)[int(21)];
    float _S2327 = (dpparams_2->primal_0[int(35)] - mean_30) * _S2326;
    float _S2328 = _S2327 + _S2327;
    float _S2329 = (dpparams_2->primal_0[int(31)] - mean_30) * _S2326;
    float _S2330 = _S2329 + _S2329;
    float _S2331 = (dpparams_2->primal_0[int(27)] - mean_30) * _S2326;
    float _S2332 = _S2331 + _S2331;
    float _S2333 = 0.3333333432674408f * (- _S2328 + - _S2330 + - _S2332);
    float _S2334 = 0.3333333432674408f * (*_s_dOut_13)[int(20)];
    float _S2335 = (dpparams_2->primal_0[int(34)] - mean_29) * _S2334;
    float _S2336 = _S2335 + _S2335;
    float _S2337 = (dpparams_2->primal_0[int(30)] - mean_29) * _S2334;
    float _S2338 = _S2337 + _S2337;
    float _S2339 = (dpparams_2->primal_0[int(26)] - mean_29) * _S2334;
    float _S2340 = _S2339 + _S2339;
    float _S2341 = 0.3333333432674408f * (- _S2336 + - _S2338 + - _S2340);
    float _S2342 = 0.3333333432674408f * (*_s_dOut_13)[int(19)];
    float _S2343 = (dpparams_2->primal_0[int(33)] - mean_28) * _S2342;
    float _S2344 = _S2343 + _S2343;
    float _S2345 = (dpparams_2->primal_0[int(29)] - mean_28) * _S2342;
    float _S2346 = _S2345 + _S2345;
    float _S2347 = (dpparams_2->primal_0[int(25)] - mean_28) * _S2342;
    float _S2348 = _S2347 + _S2347;
    float _S2349 = 0.3333333432674408f * (- _S2344 + - _S2346 + - _S2348);
    float _S2350 = 0.3333333432674408f * (*_s_dOut_13)[int(18)];
    float _S2351 = (dpparams_2->primal_0[int(32)] - mean_27) * _S2350;
    float _S2352 = _S2351 + _S2351;
    float _S2353 = (dpparams_2->primal_0[int(28)] - mean_27) * _S2350;
    float _S2354 = _S2353 + _S2353;
    float _S2355 = (dpparams_2->primal_0[int(24)] - mean_27) * _S2350;
    float _S2356 = _S2355 + _S2355;
    float _S2357 = 0.3333333432674408f * (- _S2352 + - _S2354 + - _S2356);
    float2  _S2358 = make_float2 ((*_s_dOut_13)[int(16)], (*_s_dOut_13)[int(17)]);
    Matrix<float, 2, 2>  _S2359 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2360;
    (&_S2360)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2360)->differential_0 = _S2359;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2361;
    (&_S2361)->primal_0 = _S2310.color_params_1.n_0;
    (&_S2361)->differential_0 = _S2306;
    s_bwd_prop_mul_3(&_S2360, &_S2361, _S2358);
    float2  _S2362 = make_float2 ((*_s_dOut_13)[int(14)], (*_s_dOut_13)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2363;
    (&_S2363)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2363)->differential_0 = _S2359;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2364;
    (&_S2364)->primal_0 = _S2310.color_params_1.g_0;
    (&_S2364)->differential_0 = _S2306;
    s_bwd_prop_mul_3(&_S2363, &_S2364, _S2362);
    float2  _S2365 = make_float2 ((*_s_dOut_13)[int(12)], (*_s_dOut_13)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2366;
    (&_S2366)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2366)->differential_0 = _S2359;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2367;
    (&_S2367)->primal_0 = _S2310.color_params_1.r_0;
    (&_S2367)->differential_0 = _S2306;
    s_bwd_prop_mul_3(&_S2366, &_S2367, _S2365);
    float2  _S2368 = make_float2 ((*_s_dOut_13)[int(10)], (*_s_dOut_13)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2369;
    (&_S2369)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2369)->differential_0 = _S2359;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2370;
    (&_S2370)->primal_0 = _S2310.color_params_1.b_0;
    (&_S2370)->differential_0 = _S2306;
    s_bwd_prop_mul_3(&_S2369, &_S2370, _S2368);
    ColorPPISPParams_0 _S2371 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2371)->n_0 = _S2361.differential_0;
    (&_S2371)->g_0 = _S2364.differential_0;
    (&_S2371)->r_0 = _S2367.differential_0;
    (&_S2371)->b_0 = _S2370.differential_0;
    float _S2372 = 0.3333333432674408f * (*_s_dOut_13)[int(9)];
    float _S2373 = _S2325 * _S2372;
    float _S2374 = _S2373 + _S2373;
    float _S2375 = _S2324 * _S2372;
    float _S2376 = _S2375 + _S2375;
    float _S2377 = _S2323 * _S2372;
    float _S2378 = _S2377 + _S2377;
    float _S2379 = 0.3333333432674408f * (- _S2374 + - _S2376 + - _S2378);
    float _S2380 = 0.3333333432674408f * (*_s_dOut_13)[int(8)];
    float _S2381 = _S2322 * _S2380;
    float _S2382 = _S2381 + _S2381;
    float _S2383 = _S2321 * _S2380;
    float _S2384 = _S2383 + _S2383;
    float _S2385 = _S2320 * _S2380;
    float _S2386 = _S2385 + _S2385;
    float _S2387 = 0.3333333432674408f * (- _S2382 + - _S2384 + - _S2386);
    float _S2388 = 0.3333333432674408f * (*_s_dOut_13)[int(7)];
    float _S2389 = _S2319 * _S2388;
    float _S2390 = _S2389 + _S2389;
    float _S2391 = _S2318 * _S2388;
    float _S2392 = _S2391 + _S2391;
    float _S2393 = _S2317 * _S2388;
    float _S2394 = _S2393 + _S2393;
    float _S2395 = 0.3333333432674408f * (- _S2390 + - _S2392 + - _S2394);
    float _S2396 = 0.3333333432674408f * (*_s_dOut_13)[int(6)];
    float _S2397 = _S2316 * _S2396;
    float _S2398 = _S2397 + _S2397;
    float _S2399 = _S2315 * _S2396;
    float _S2400 = _S2399 + _S2399;
    float _S2401 = _S2314 * _S2396;
    float _S2402 = _S2401 + _S2401;
    float _S2403 = 0.3333333432674408f * (- _S2398 + - _S2400 + - _S2402);
    float _S2404 = 0.3333333432674408f * (*_s_dOut_13)[int(5)];
    float _S2405 = _S2313 * _S2404;
    float _S2406 = _S2405 + _S2405;
    float _S2407 = _S2312 * _S2404;
    float _S2408 = _S2407 + _S2407;
    float _S2409 = _S2311 * _S2404;
    float _S2410 = _S2409 + _S2409;
    float _S2411 = 0.3333333432674408f * (- _S2406 + - _S2408 + - _S2410);
    DiffPair_float_0 _S2412;
    (&_S2412)->primal_0 = 0.0f;
    (&_S2412)->differential_0 = 0.0f;
    DiffPair_float_0 _S2413;
    (&_S2413)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S2413)->differential_0 = 0.0f;
    _d_max_0(&_S2412, &_S2413, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2414;
    (&_S2414)->primal_0 = 0.0f;
    (&_S2414)->differential_0 = 0.0f;
    DiffPair_float_0 _S2415;
    (&_S2415)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S2415)->differential_0 = 0.0f;
    _d_max_0(&_S2414, &_S2415, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2416;
    (&_S2416)->primal_0 = 0.0f;
    (&_S2416)->differential_0 = 0.0f;
    DiffPair_float_0 _S2417;
    (&_S2417)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S2417)->differential_0 = 0.0f;
    _d_max_0(&_S2416, &_S2417, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2418;
    (&_S2418)->primal_0 = 0.0f;
    (&_S2418)->differential_0 = 0.0f;
    DiffPair_float_0 _S2419;
    (&_S2419)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S2419)->differential_0 = 0.0f;
    _d_max_0(&_S2418, &_S2419, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2420;
    (&_S2420)->primal_0 = 0.0f;
    (&_S2420)->differential_0 = 0.0f;
    DiffPair_float_0 _S2421;
    (&_S2421)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S2421)->differential_0 = 0.0f;
    _d_max_0(&_S2420, &_S2421, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2422;
    (&_S2422)->primal_0 = 0.0f;
    (&_S2422)->differential_0 = 0.0f;
    DiffPair_float_0 _S2423;
    (&_S2423)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S2423)->differential_0 = 0.0f;
    _d_max_0(&_S2422, &_S2423, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2424;
    (&_S2424)->primal_0 = 0.0f;
    (&_S2424)->differential_0 = 0.0f;
    DiffPair_float_0 _S2425;
    (&_S2425)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S2425)->differential_0 = 0.0f;
    _d_max_0(&_S2424, &_S2425, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2426;
    (&_S2426)->primal_0 = 0.0f;
    (&_S2426)->differential_0 = 0.0f;
    DiffPair_float_0 _S2427;
    (&_S2427)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S2427)->differential_0 = 0.0f;
    _d_max_0(&_S2426, &_S2427, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2428;
    (&_S2428)->primal_0 = 0.0f;
    (&_S2428)->differential_0 = 0.0f;
    DiffPair_float_0 _S2429;
    (&_S2429)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S2429)->differential_0 = 0.0f;
    _d_max_0(&_S2428, &_S2429, (*_s_dOut_13)[int(2)]);
    float _S2430 = dpparams_2->primal_0[int(12)] * (*_s_dOut_13)[int(1)];
    float _S2431 = dpparams_2->primal_0[int(11)] * (*_s_dOut_13)[int(1)];
    float _S2432 = dpparams_2->primal_0[int(7)] * (*_s_dOut_13)[int(1)];
    float _S2433 = dpparams_2->primal_0[int(6)] * (*_s_dOut_13)[int(1)];
    float _S2434 = dpparams_2->primal_0[int(2)] * (*_s_dOut_13)[int(1)];
    float _S2435 = dpparams_2->primal_0[int(1)] * (*_s_dOut_13)[int(1)];
    PPISPParams_0 _S2436 = PPISPParams_x24_syn_dzero_0();
    (&_S2436)->color_params_1 = _S2371;
    (&_S2436)->exposure_1 = (*_s_dOut_13)[int(0)];
    _S2310 = _S2436;
    (&(&_S2310)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2437 = _S2328 + _S2333 + _S2436.crf_params_1[int(2)].center_0;
    (&(&_S2310)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2438 = _S2336 + _S2341 + _S2436.crf_params_1[int(2)].gamma_0;
    (&(&_S2310)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2439 = _S2344 + _S2349 + _S2436.crf_params_1[int(2)].shoulder_0;
    (&(&_S2310)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2440 = _S2352 + _S2357 + _S2436.crf_params_1[int(2)].toe_0;
    (&(&_S2310)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2441 = _S2330 + _S2333 + _S2436.crf_params_1[int(1)].center_0;
    (&(&_S2310)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2442 = _S2338 + _S2341 + _S2436.crf_params_1[int(1)].gamma_0;
    (&(&_S2310)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2443 = _S2346 + _S2349 + _S2436.crf_params_1[int(1)].shoulder_0;
    (&(&_S2310)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2444 = _S2354 + _S2357 + _S2436.crf_params_1[int(1)].toe_0;
    (&(&_S2310)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2445 = _S2332 + _S2333 + _S2436.crf_params_1[int(0)].center_0;
    (&(&_S2310)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2446 = _S2340 + _S2341 + _S2436.crf_params_1[int(0)].gamma_0;
    (&(&_S2310)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2447 = _S2348 + _S2349 + _S2436.crf_params_1[int(0)].shoulder_0;
    (&(&_S2310)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2448 = _S2356 + _S2357 + _S2436.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S2310)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S2310)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S2310)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2449 = _S2374 + _S2379 + _S2413.differential_0 + _S2436.vignette_params_1[int(2)].alpha2_0;
    (&(&_S2310)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2450 = _S2382 + _S2387 + _S2419.differential_0 + _S2436.vignette_params_1[int(2)].alpha1_0;
    (&(&_S2310)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2451 = _S2390 + _S2395 + _S2425.differential_0 + _S2436.vignette_params_1[int(2)].alpha0_0;
    (&(&_S2310)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2452 = _S2398 + _S2403 + _S2430 + _S2430 + _S2436.vignette_params_1[int(2)].cy_0;
    (&(&_S2310)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2453 = _S2406 + _S2411 + _S2431 + _S2431 + _S2436.vignette_params_1[int(2)].cx_0;
    (&(&_S2310)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2454 = _S2376 + _S2379 + _S2415.differential_0 + _S2436.vignette_params_1[int(1)].alpha2_0;
    (&(&_S2310)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2455 = _S2384 + _S2387 + _S2421.differential_0 + _S2436.vignette_params_1[int(1)].alpha1_0;
    (&(&_S2310)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2456 = _S2392 + _S2395 + _S2427.differential_0 + _S2436.vignette_params_1[int(1)].alpha0_0;
    (&(&_S2310)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2457 = _S2400 + _S2403 + _S2432 + _S2432 + _S2436.vignette_params_1[int(1)].cy_0;
    (&(&_S2310)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2458 = _S2408 + _S2411 + _S2433 + _S2433 + _S2436.vignette_params_1[int(1)].cx_0;
    (&(&_S2310)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2459 = _S2378 + _S2379 + _S2417.differential_0 + _S2436.vignette_params_1[int(0)].alpha2_0;
    (&(&_S2310)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2460 = _S2386 + _S2387 + _S2423.differential_0 + _S2436.vignette_params_1[int(0)].alpha1_0;
    (&(&_S2310)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2461 = _S2394 + _S2395 + _S2429.differential_0 + _S2436.vignette_params_1[int(0)].alpha0_0;
    (&(&_S2310)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2462 = _S2402 + _S2403 + _S2434 + _S2434 + _S2436.vignette_params_1[int(0)].cy_0;
    (&(&_S2310)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2463 = _S2410 + _S2411 + _S2435 + _S2435 + _S2436.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2464;
    _S2464[int(0)] = 0.0f;
    _S2464[int(1)] = 0.0f;
    _S2464[int(2)] = 0.0f;
    _S2464[int(3)] = 0.0f;
    _S2464[int(4)] = 0.0f;
    _S2464[int(5)] = 0.0f;
    _S2464[int(6)] = 0.0f;
    _S2464[int(7)] = 0.0f;
    _S2464[int(8)] = 0.0f;
    _S2464[int(9)] = 0.0f;
    _S2464[int(10)] = 0.0f;
    _S2464[int(11)] = 0.0f;
    _S2464[int(12)] = 0.0f;
    _S2464[int(13)] = 0.0f;
    _S2464[int(14)] = 0.0f;
    _S2464[int(15)] = 0.0f;
    _S2464[int(16)] = 0.0f;
    _S2464[int(17)] = 0.0f;
    _S2464[int(18)] = 0.0f;
    _S2464[int(19)] = 0.0f;
    _S2464[int(20)] = 0.0f;
    _S2464[int(21)] = 0.0f;
    _S2464[int(22)] = 0.0f;
    _S2464[int(23)] = 0.0f;
    _S2464[int(24)] = 0.0f;
    _S2464[int(25)] = 0.0f;
    _S2464[int(26)] = 0.0f;
    _S2464[int(27)] = 0.0f;
    _S2464[int(28)] = 0.0f;
    _S2464[int(29)] = 0.0f;
    _S2464[int(30)] = 0.0f;
    _S2464[int(31)] = 0.0f;
    _S2464[int(32)] = 0.0f;
    _S2464[int(33)] = 0.0f;
    _S2464[int(34)] = 0.0f;
    _S2464[int(35)] = 0.0f;
    _S2464[int(8)] = _S2456;
    _S2464[int(16)] = _S2436.color_params_1.b_0.x;
    _S2464[int(15)] = _S2449;
    _S2464[int(14)] = _S2450;
    _S2464[int(13)] = _S2451;
    _S2464[int(12)] = _S2452;
    _S2464[int(11)] = _S2453;
    _S2464[int(10)] = _S2454;
    _S2464[int(9)] = _S2455;
    _S2464[int(17)] = _S2436.color_params_1.b_0.y;
    _S2464[int(7)] = _S2457;
    _S2464[int(6)] = _S2458;
    _S2464[int(5)] = _S2459;
    _S2464[int(4)] = _S2460;
    _S2464[int(3)] = _S2461;
    _S2464[int(2)] = _S2462;
    _S2464[int(1)] = _S2463;
    _S2464[int(0)] = _S2310.exposure_1;
    _S2464[int(26)] = _S2446;
    _S2464[int(34)] = _S2438;
    _S2464[int(33)] = _S2439;
    _S2464[int(32)] = _S2440;
    _S2464[int(31)] = _S2441;
    _S2464[int(30)] = _S2442;
    _S2464[int(29)] = _S2443;
    _S2464[int(28)] = _S2444;
    _S2464[int(27)] = _S2445;
    _S2464[int(35)] = _S2437;
    _S2464[int(25)] = _S2447;
    _S2464[int(24)] = _S2448;
    _S2464[int(23)] = _S2436.color_params_1.n_0.y;
    _S2464[int(22)] = _S2436.color_params_1.n_0.x;
    _S2464[int(21)] = _S2436.color_params_1.g_0.y;
    _S2464[int(20)] = _S2436.color_params_1.g_0.x;
    _S2464[int(19)] = _S2436.color_params_1.r_0.y;
    _S2464[int(18)] = _S2436.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S2464;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2465, FixedArray<float, 22>  * _S2466)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S2465, _S2466);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_6, FixedArray<float, 22>  grad_out_2, FixedArray<float, 36>  * _S2467)
{
    FixedArray<float, 36>  _S2468 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_6;
    (&dp_params_2)->differential_0 = _S2468;
    FixedArray<float, 22>  _S2469 = grad_out_2;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, &_S2469);
    *_S2467 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_14)
{
    VignettingChannelParams_0 _S2470 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2471 = {
        _S2470, _S2470, _S2470
    };
    float2  _S2472 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2473 = { _S2472, _S2472, _S2472, _S2472 };
    RQSCRFPPISPChannelParams_0 _S2474 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2475 = {
        _S2474, _S2474, _S2474
    };
    PPISPParamsRQS_0 _S2476;
    (&_S2476)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S2476)->vignette_params_0 = _S2471;
    (&_S2476)->color_params_0 = _S2473;
    (&_S2476)->crf_params_0 = _S2475;
    (&(&_S2476)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S2476)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S2476)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S2476)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S2476)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S2476)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S2476)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S2476)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S2476)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S2476)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S2476)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S2476)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S2476)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S2476)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S2476)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S2476)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S2476)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S2476)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S2476)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S2476)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S2476)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S2476)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S2476)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S2476)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S2476)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S2476)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S2476)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S2476)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S2476)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S2476)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S2476)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S2476)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S2476)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S2476)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S2476)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S2476)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S2476)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S2476)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_31 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S2477 = dpparams_3->primal_0[int(1)] - mean_31;
    float _S2478 = dpparams_3->primal_0[int(6)] - mean_31;
    float _S2479 = dpparams_3->primal_0[int(11)] - mean_31;
    float mean_32 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S2480 = dpparams_3->primal_0[int(2)] - mean_32;
    float _S2481 = dpparams_3->primal_0[int(7)] - mean_32;
    float _S2482 = dpparams_3->primal_0[int(12)] - mean_32;
    float mean_33 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S2483 = dpparams_3->primal_0[int(3)] - mean_33;
    float _S2484 = dpparams_3->primal_0[int(8)] - mean_33;
    float _S2485 = dpparams_3->primal_0[int(13)] - mean_33;
    float mean_34 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S2486 = dpparams_3->primal_0[int(4)] - mean_34;
    float _S2487 = dpparams_3->primal_0[int(9)] - mean_34;
    float _S2488 = dpparams_3->primal_0[int(14)] - mean_34;
    float mean_35 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S2489 = dpparams_3->primal_0[int(5)] - mean_35;
    float _S2490 = dpparams_3->primal_0[int(10)] - mean_35;
    float _S2491 = dpparams_3->primal_0[int(15)] - mean_35;
    float mean_36 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_37 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_38 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_39 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_40 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S2492 = 0.3333333432674408f * (*_s_dOut_14)[int(22)];
    float _S2493 = (dpparams_3->primal_0[int(38)] - mean_40) * _S2492;
    float _S2494 = _S2493 + _S2493;
    float _S2495 = (dpparams_3->primal_0[int(33)] - mean_40) * _S2492;
    float _S2496 = _S2495 + _S2495;
    float _S2497 = (dpparams_3->primal_0[int(28)] - mean_40) * _S2492;
    float _S2498 = _S2497 + _S2497;
    float _S2499 = 0.3333333432674408f * (- _S2494 + - _S2496 + - _S2498);
    float _S2500 = 0.3333333432674408f * (*_s_dOut_14)[int(21)];
    float _S2501 = (dpparams_3->primal_0[int(37)] - mean_39) * _S2500;
    float _S2502 = _S2501 + _S2501;
    float _S2503 = (dpparams_3->primal_0[int(32)] - mean_39) * _S2500;
    float _S2504 = _S2503 + _S2503;
    float _S2505 = (dpparams_3->primal_0[int(27)] - mean_39) * _S2500;
    float _S2506 = _S2505 + _S2505;
    float _S2507 = 0.3333333432674408f * (- _S2502 + - _S2504 + - _S2506);
    float _S2508 = 0.3333333432674408f * (*_s_dOut_14)[int(20)];
    float _S2509 = (dpparams_3->primal_0[int(36)] - mean_38) * _S2508;
    float _S2510 = _S2509 + _S2509;
    float _S2511 = (dpparams_3->primal_0[int(31)] - mean_38) * _S2508;
    float _S2512 = _S2511 + _S2511;
    float _S2513 = (dpparams_3->primal_0[int(26)] - mean_38) * _S2508;
    float _S2514 = _S2513 + _S2513;
    float _S2515 = 0.3333333432674408f * (- _S2510 + - _S2512 + - _S2514);
    float _S2516 = 0.3333333432674408f * (*_s_dOut_14)[int(19)];
    float _S2517 = (dpparams_3->primal_0[int(35)] - mean_37) * _S2516;
    float _S2518 = _S2517 + _S2517;
    float _S2519 = (dpparams_3->primal_0[int(30)] - mean_37) * _S2516;
    float _S2520 = _S2519 + _S2519;
    float _S2521 = (dpparams_3->primal_0[int(25)] - mean_37) * _S2516;
    float _S2522 = _S2521 + _S2521;
    float _S2523 = 0.3333333432674408f * (- _S2518 + - _S2520 + - _S2522);
    float _S2524 = 0.3333333432674408f * (*_s_dOut_14)[int(18)];
    float _S2525 = (dpparams_3->primal_0[int(34)] - mean_36) * _S2524;
    float _S2526 = _S2525 + _S2525;
    float _S2527 = (dpparams_3->primal_0[int(29)] - mean_36) * _S2524;
    float _S2528 = _S2527 + _S2527;
    float _S2529 = (dpparams_3->primal_0[int(24)] - mean_36) * _S2524;
    float _S2530 = _S2529 + _S2529;
    float _S2531 = 0.3333333432674408f * (- _S2526 + - _S2528 + - _S2530);
    float2  _S2532 = make_float2 ((*_s_dOut_14)[int(16)], (*_s_dOut_14)[int(17)]);
    Matrix<float, 2, 2>  _S2533 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2534;
    (&_S2534)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2534)->differential_0 = _S2533;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2535;
    (&_S2535)->primal_0 = _S2476.color_params_0.n_0;
    (&_S2535)->differential_0 = _S2472;
    s_bwd_prop_mul_3(&_S2534, &_S2535, _S2532);
    float2  _S2536 = make_float2 ((*_s_dOut_14)[int(14)], (*_s_dOut_14)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2537;
    (&_S2537)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2537)->differential_0 = _S2533;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2538;
    (&_S2538)->primal_0 = _S2476.color_params_0.g_0;
    (&_S2538)->differential_0 = _S2472;
    s_bwd_prop_mul_3(&_S2537, &_S2538, _S2536);
    float2  _S2539 = make_float2 ((*_s_dOut_14)[int(12)], (*_s_dOut_14)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2540;
    (&_S2540)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2540)->differential_0 = _S2533;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2541;
    (&_S2541)->primal_0 = _S2476.color_params_0.r_0;
    (&_S2541)->differential_0 = _S2472;
    s_bwd_prop_mul_3(&_S2540, &_S2541, _S2539);
    float2  _S2542 = make_float2 ((*_s_dOut_14)[int(10)], (*_s_dOut_14)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2543;
    (&_S2543)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2543)->differential_0 = _S2533;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2544;
    (&_S2544)->primal_0 = _S2476.color_params_0.b_0;
    (&_S2544)->differential_0 = _S2472;
    s_bwd_prop_mul_3(&_S2543, &_S2544, _S2542);
    ColorPPISPParams_0 _S2545 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2545)->n_0 = _S2535.differential_0;
    (&_S2545)->g_0 = _S2538.differential_0;
    (&_S2545)->r_0 = _S2541.differential_0;
    (&_S2545)->b_0 = _S2544.differential_0;
    float _S2546 = 0.3333333432674408f * (*_s_dOut_14)[int(9)];
    float _S2547 = _S2491 * _S2546;
    float _S2548 = _S2547 + _S2547;
    float _S2549 = _S2490 * _S2546;
    float _S2550 = _S2549 + _S2549;
    float _S2551 = _S2489 * _S2546;
    float _S2552 = _S2551 + _S2551;
    float _S2553 = 0.3333333432674408f * (- _S2548 + - _S2550 + - _S2552);
    float _S2554 = 0.3333333432674408f * (*_s_dOut_14)[int(8)];
    float _S2555 = _S2488 * _S2554;
    float _S2556 = _S2555 + _S2555;
    float _S2557 = _S2487 * _S2554;
    float _S2558 = _S2557 + _S2557;
    float _S2559 = _S2486 * _S2554;
    float _S2560 = _S2559 + _S2559;
    float _S2561 = 0.3333333432674408f * (- _S2556 + - _S2558 + - _S2560);
    float _S2562 = 0.3333333432674408f * (*_s_dOut_14)[int(7)];
    float _S2563 = _S2485 * _S2562;
    float _S2564 = _S2563 + _S2563;
    float _S2565 = _S2484 * _S2562;
    float _S2566 = _S2565 + _S2565;
    float _S2567 = _S2483 * _S2562;
    float _S2568 = _S2567 + _S2567;
    float _S2569 = 0.3333333432674408f * (- _S2564 + - _S2566 + - _S2568);
    float _S2570 = 0.3333333432674408f * (*_s_dOut_14)[int(6)];
    float _S2571 = _S2482 * _S2570;
    float _S2572 = _S2571 + _S2571;
    float _S2573 = _S2481 * _S2570;
    float _S2574 = _S2573 + _S2573;
    float _S2575 = _S2480 * _S2570;
    float _S2576 = _S2575 + _S2575;
    float _S2577 = 0.3333333432674408f * (- _S2572 + - _S2574 + - _S2576);
    float _S2578 = 0.3333333432674408f * (*_s_dOut_14)[int(5)];
    float _S2579 = _S2479 * _S2578;
    float _S2580 = _S2579 + _S2579;
    float _S2581 = _S2478 * _S2578;
    float _S2582 = _S2581 + _S2581;
    float _S2583 = _S2477 * _S2578;
    float _S2584 = _S2583 + _S2583;
    float _S2585 = 0.3333333432674408f * (- _S2580 + - _S2582 + - _S2584);
    DiffPair_float_0 _S2586;
    (&_S2586)->primal_0 = 0.0f;
    (&_S2586)->differential_0 = 0.0f;
    DiffPair_float_0 _S2587;
    (&_S2587)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S2587)->differential_0 = 0.0f;
    _d_max_0(&_S2586, &_S2587, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2588;
    (&_S2588)->primal_0 = 0.0f;
    (&_S2588)->differential_0 = 0.0f;
    DiffPair_float_0 _S2589;
    (&_S2589)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S2589)->differential_0 = 0.0f;
    _d_max_0(&_S2588, &_S2589, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2590;
    (&_S2590)->primal_0 = 0.0f;
    (&_S2590)->differential_0 = 0.0f;
    DiffPair_float_0 _S2591;
    (&_S2591)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S2591)->differential_0 = 0.0f;
    _d_max_0(&_S2590, &_S2591, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2592;
    (&_S2592)->primal_0 = 0.0f;
    (&_S2592)->differential_0 = 0.0f;
    DiffPair_float_0 _S2593;
    (&_S2593)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S2593)->differential_0 = 0.0f;
    _d_max_0(&_S2592, &_S2593, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2594;
    (&_S2594)->primal_0 = 0.0f;
    (&_S2594)->differential_0 = 0.0f;
    DiffPair_float_0 _S2595;
    (&_S2595)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S2595)->differential_0 = 0.0f;
    _d_max_0(&_S2594, &_S2595, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2596;
    (&_S2596)->primal_0 = 0.0f;
    (&_S2596)->differential_0 = 0.0f;
    DiffPair_float_0 _S2597;
    (&_S2597)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S2597)->differential_0 = 0.0f;
    _d_max_0(&_S2596, &_S2597, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2598;
    (&_S2598)->primal_0 = 0.0f;
    (&_S2598)->differential_0 = 0.0f;
    DiffPair_float_0 _S2599;
    (&_S2599)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S2599)->differential_0 = 0.0f;
    _d_max_0(&_S2598, &_S2599, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2600;
    (&_S2600)->primal_0 = 0.0f;
    (&_S2600)->differential_0 = 0.0f;
    DiffPair_float_0 _S2601;
    (&_S2601)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S2601)->differential_0 = 0.0f;
    _d_max_0(&_S2600, &_S2601, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2602;
    (&_S2602)->primal_0 = 0.0f;
    (&_S2602)->differential_0 = 0.0f;
    DiffPair_float_0 _S2603;
    (&_S2603)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S2603)->differential_0 = 0.0f;
    _d_max_0(&_S2602, &_S2603, (*_s_dOut_14)[int(2)]);
    float _S2604 = dpparams_3->primal_0[int(12)] * (*_s_dOut_14)[int(1)];
    float _S2605 = dpparams_3->primal_0[int(11)] * (*_s_dOut_14)[int(1)];
    float _S2606 = dpparams_3->primal_0[int(7)] * (*_s_dOut_14)[int(1)];
    float _S2607 = dpparams_3->primal_0[int(6)] * (*_s_dOut_14)[int(1)];
    float _S2608 = dpparams_3->primal_0[int(2)] * (*_s_dOut_14)[int(1)];
    float _S2609 = dpparams_3->primal_0[int(1)] * (*_s_dOut_14)[int(1)];
    PPISPParamsRQS_0 _S2610 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2610)->color_params_0 = _S2545;
    (&_S2610)->exposure_0 = (*_s_dOut_14)[int(0)];
    _S2476 = _S2610;
    (&(&_S2476)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2611 = _S2494 + _S2499 + _S2610.crf_params_0[int(2)].gc_0;
    (&(&_S2476)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2612 = _S2502 + _S2507 + _S2610.crf_params_0[int(2)].y0_0;
    (&(&_S2476)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2613 = _S2510 + _S2515 + _S2610.crf_params_0[int(2)].x0_0;
    (&(&_S2476)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2614 = _S2518 + _S2523 + _S2610.crf_params_0[int(2)].g1_0;
    (&(&_S2476)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2615 = _S2526 + _S2531 + _S2610.crf_params_0[int(2)].g0_0;
    (&(&_S2476)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2616 = _S2496 + _S2499 + _S2610.crf_params_0[int(1)].gc_0;
    (&(&_S2476)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2617 = _S2504 + _S2507 + _S2610.crf_params_0[int(1)].y0_0;
    (&(&_S2476)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2618 = _S2512 + _S2515 + _S2610.crf_params_0[int(1)].x0_0;
    (&(&_S2476)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2619 = _S2520 + _S2523 + _S2610.crf_params_0[int(1)].g1_0;
    (&(&_S2476)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2620 = _S2528 + _S2531 + _S2610.crf_params_0[int(1)].g0_0;
    (&(&_S2476)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2621 = _S2498 + _S2499 + _S2610.crf_params_0[int(0)].gc_0;
    (&(&_S2476)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2622 = _S2506 + _S2507 + _S2610.crf_params_0[int(0)].y0_0;
    (&(&_S2476)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2623 = _S2514 + _S2515 + _S2610.crf_params_0[int(0)].x0_0;
    (&(&_S2476)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2624 = _S2522 + _S2523 + _S2610.crf_params_0[int(0)].g1_0;
    (&(&_S2476)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2625 = _S2530 + _S2531 + _S2610.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S2476)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2476)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2476)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2626 = _S2548 + _S2553 + _S2587.differential_0 + _S2610.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2476)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2627 = _S2556 + _S2561 + _S2593.differential_0 + _S2610.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2476)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2628 = _S2564 + _S2569 + _S2599.differential_0 + _S2610.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2476)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2629 = _S2572 + _S2577 + _S2604 + _S2604 + _S2610.vignette_params_0[int(2)].cy_0;
    (&(&_S2476)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2630 = _S2580 + _S2585 + _S2605 + _S2605 + _S2610.vignette_params_0[int(2)].cx_0;
    (&(&_S2476)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2631 = _S2550 + _S2553 + _S2589.differential_0 + _S2610.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2476)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2632 = _S2558 + _S2561 + _S2595.differential_0 + _S2610.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2476)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2633 = _S2566 + _S2569 + _S2601.differential_0 + _S2610.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2476)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2634 = _S2574 + _S2577 + _S2606 + _S2606 + _S2610.vignette_params_0[int(1)].cy_0;
    (&(&_S2476)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2635 = _S2582 + _S2585 + _S2607 + _S2607 + _S2610.vignette_params_0[int(1)].cx_0;
    (&(&_S2476)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2636 = _S2552 + _S2553 + _S2591.differential_0 + _S2610.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2476)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2637 = _S2560 + _S2561 + _S2597.differential_0 + _S2610.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2476)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2638 = _S2568 + _S2569 + _S2603.differential_0 + _S2610.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2476)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2639 = _S2576 + _S2577 + _S2608 + _S2608 + _S2610.vignette_params_0[int(0)].cy_0;
    (&(&_S2476)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2640 = _S2584 + _S2585 + _S2609 + _S2609 + _S2610.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2641;
    _S2641[int(0)] = 0.0f;
    _S2641[int(1)] = 0.0f;
    _S2641[int(2)] = 0.0f;
    _S2641[int(3)] = 0.0f;
    _S2641[int(4)] = 0.0f;
    _S2641[int(5)] = 0.0f;
    _S2641[int(6)] = 0.0f;
    _S2641[int(7)] = 0.0f;
    _S2641[int(8)] = 0.0f;
    _S2641[int(9)] = 0.0f;
    _S2641[int(10)] = 0.0f;
    _S2641[int(11)] = 0.0f;
    _S2641[int(12)] = 0.0f;
    _S2641[int(13)] = 0.0f;
    _S2641[int(14)] = 0.0f;
    _S2641[int(15)] = 0.0f;
    _S2641[int(16)] = 0.0f;
    _S2641[int(17)] = 0.0f;
    _S2641[int(18)] = 0.0f;
    _S2641[int(19)] = 0.0f;
    _S2641[int(20)] = 0.0f;
    _S2641[int(21)] = 0.0f;
    _S2641[int(22)] = 0.0f;
    _S2641[int(23)] = 0.0f;
    _S2641[int(24)] = 0.0f;
    _S2641[int(25)] = 0.0f;
    _S2641[int(26)] = 0.0f;
    _S2641[int(27)] = 0.0f;
    _S2641[int(28)] = 0.0f;
    _S2641[int(29)] = 0.0f;
    _S2641[int(30)] = 0.0f;
    _S2641[int(31)] = 0.0f;
    _S2641[int(32)] = 0.0f;
    _S2641[int(33)] = 0.0f;
    _S2641[int(34)] = 0.0f;
    _S2641[int(35)] = 0.0f;
    _S2641[int(36)] = 0.0f;
    _S2641[int(37)] = 0.0f;
    _S2641[int(38)] = 0.0f;
    _S2641[int(9)] = _S2632;
    _S2641[int(18)] = _S2610.color_params_0.r_0.x;
    _S2641[int(17)] = _S2610.color_params_0.b_0.y;
    _S2641[int(16)] = _S2610.color_params_0.b_0.x;
    _S2641[int(15)] = _S2626;
    _S2641[int(14)] = _S2627;
    _S2641[int(13)] = _S2628;
    _S2641[int(12)] = _S2629;
    _S2641[int(11)] = _S2630;
    _S2641[int(10)] = _S2631;
    _S2641[int(19)] = _S2610.color_params_0.r_0.y;
    _S2641[int(8)] = _S2633;
    _S2641[int(7)] = _S2634;
    _S2641[int(6)] = _S2635;
    _S2641[int(5)] = _S2636;
    _S2641[int(4)] = _S2637;
    _S2641[int(3)] = _S2638;
    _S2641[int(2)] = _S2639;
    _S2641[int(1)] = _S2640;
    _S2641[int(0)] = _S2476.exposure_0;
    _S2641[int(28)] = _S2621;
    _S2641[int(37)] = _S2612;
    _S2641[int(36)] = _S2613;
    _S2641[int(35)] = _S2614;
    _S2641[int(34)] = _S2615;
    _S2641[int(33)] = _S2616;
    _S2641[int(32)] = _S2617;
    _S2641[int(31)] = _S2618;
    _S2641[int(30)] = _S2619;
    _S2641[int(29)] = _S2620;
    _S2641[int(38)] = _S2611;
    _S2641[int(27)] = _S2622;
    _S2641[int(26)] = _S2623;
    _S2641[int(25)] = _S2624;
    _S2641[int(24)] = _S2625;
    _S2641[int(23)] = _S2610.color_params_0.n_0.y;
    _S2641[int(22)] = _S2610.color_params_0.n_0.x;
    _S2641[int(21)] = _S2610.color_params_0.g_0.y;
    _S2641[int(20)] = _S2610.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S2641;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2642, FixedArray<float, 23>  * _S2643)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S2642, _S2643);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_7, FixedArray<float, 23>  grad_out_3, FixedArray<float, 39>  * _S2644)
{
    FixedArray<float, 39>  _S2645 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_7;
    (&dp_params_3)->differential_0 = _S2645;
    FixedArray<float, 23>  _S2646 = grad_out_3;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, &_S2646);
    *_S2644 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_2, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S2647)
{
    float _S2648;
    FixedArray<float, 6>  losses_5;
    float _S2649 = float(num_cameras_0);
    float _S2650 = raw_losses_2[int(0)] / _S2649;
    for(;;)
    {
        float _S2651 = (F32_abs((_S2650)));
        if(_S2651 < 0.10000000149011612f)
        {
            _S2648 = 0.5f * _S2650 * _S2650 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2648 = _S2651 - 0.05000000074505806f;
            break;
        }
    }
    losses_5[int(0)] = _S2648;
    losses_5[int(1)] = raw_losses_2[int(1)] / (3.0f * _S2649);
    losses_5[int(2)] = (raw_losses_2[int(2)] + raw_losses_2[int(3)] + raw_losses_2[int(4)]) / (9.0f * _S2649);
    losses_5[int(3)] = (raw_losses_2[int(5)] + raw_losses_2[int(6)] + raw_losses_2[int(7)] + raw_losses_2[int(8)] + raw_losses_2[int(9)]) / (5.0f * _S2649);
    float _S2652 = raw_losses_2[int(10)] / _S2649;
    for(;;)
    {
        float _S2653 = (F32_abs((_S2652)));
        if(_S2653 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2652 * _S2652 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2653 - 0.00249999994412065f;
            break;
        }
    }
    float _S2654;
    float _S2655 = raw_losses_2[int(11)] / _S2649;
    for(;;)
    {
        float _S2656 = (F32_abs((_S2655)));
        if(_S2656 < 0.00499999988824129f)
        {
            _S2654 = 0.5f * _S2655 * _S2655 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2654 = _S2656 - 0.00249999994412065f;
            break;
        }
    }
    float _S2657 = _S2648 + _S2654;
    float _S2658 = raw_losses_2[int(12)] / _S2649;
    for(;;)
    {
        float _S2659 = (F32_abs((_S2658)));
        if(_S2659 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2658 * _S2658 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2659 - 0.00249999994412065f;
            break;
        }
    }
    float _S2660 = _S2657 + _S2648;
    float _S2661 = raw_losses_2[int(13)] / _S2649;
    for(;;)
    {
        float _S2662 = (F32_abs((_S2661)));
        if(_S2662 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2661 * _S2661 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2662 - 0.00249999994412065f;
            break;
        }
    }
    float _S2663 = _S2660 + _S2648;
    float _S2664 = raw_losses_2[int(14)] / _S2649;
    for(;;)
    {
        float _S2665 = (F32_abs((_S2664)));
        if(_S2665 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2664 * _S2664 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2665 - 0.00249999994412065f;
            break;
        }
    }
    float _S2666 = _S2663 + _S2648;
    float _S2667 = raw_losses_2[int(15)] / _S2649;
    for(;;)
    {
        float _S2668 = (F32_abs((_S2667)));
        if(_S2668 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2667 * _S2667 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2668 - 0.00249999994412065f;
            break;
        }
    }
    float _S2669 = _S2666 + _S2648;
    float _S2670 = raw_losses_2[int(16)] / _S2649;
    for(;;)
    {
        float _S2671 = (F32_abs((_S2670)));
        if(_S2671 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2670 * _S2670 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2671 - 0.00249999994412065f;
            break;
        }
    }
    float _S2672 = _S2669 + _S2648;
    float _S2673 = raw_losses_2[int(17)] / _S2649;
    for(;;)
    {
        float _S2674 = (F32_abs((_S2673)));
        if(_S2674 < 0.00499999988824129f)
        {
            _S2648 = 0.5f * _S2673 * _S2673 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2648 = _S2674 - 0.00249999994412065f;
            break;
        }
    }
    float _S2675 = (_S2672 + _S2648) / 8.0f;
    float _S2676 = (raw_losses_2[int(18)] + raw_losses_2[int(19)] + raw_losses_2[int(20)] + raw_losses_2[int(21)]) / (4.0f * _S2649);
    losses_5[int(0)] = losses_5[int(0)] * loss_weights_0[int(0)];
    losses_5[int(1)] = losses_5[int(1)] * loss_weights_0[int(1)];
    losses_5[int(2)] = losses_5[int(2)] * loss_weights_0[int(2)];
    losses_5[int(3)] = losses_5[int(3)] * loss_weights_0[int(3)];
    losses_5[int(4)] = _S2675 * loss_weights_0[int(4)];
    losses_5[int(5)] = _S2676 * loss_weights_0[int(5)];
    *_S2647 = losses_5;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_3, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S2677)
{
    float _S2678;
    FixedArray<float, 6>  losses_6;
    float _S2679 = float(num_cameras_1);
    float _S2680 = raw_losses_3[int(0)] / _S2679;
    for(;;)
    {
        float _S2681 = (F32_abs((_S2680)));
        if(_S2681 < 0.10000000149011612f)
        {
            _S2678 = 0.5f * _S2680 * _S2680 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2678 = _S2681 - 0.05000000074505806f;
            break;
        }
    }
    losses_6[int(0)] = _S2678;
    losses_6[int(1)] = raw_losses_3[int(1)] / (3.0f * _S2679);
    losses_6[int(2)] = (raw_losses_3[int(2)] + raw_losses_3[int(3)] + raw_losses_3[int(4)]) / (9.0f * _S2679);
    float _S2682 = 5.0f * _S2679;
    losses_6[int(3)] = (raw_losses_3[int(5)] + raw_losses_3[int(6)] + raw_losses_3[int(7)] + raw_losses_3[int(8)] + raw_losses_3[int(9)]) / _S2682;
    float _S2683 = raw_losses_3[int(10)] / _S2679;
    for(;;)
    {
        float _S2684 = (F32_abs((_S2683)));
        if(_S2684 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2683 * _S2683 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2684 - 0.00249999994412065f;
            break;
        }
    }
    float _S2685;
    float _S2686 = raw_losses_3[int(11)] / _S2679;
    for(;;)
    {
        float _S2687 = (F32_abs((_S2686)));
        if(_S2687 < 0.00499999988824129f)
        {
            _S2685 = 0.5f * _S2686 * _S2686 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2685 = _S2687 - 0.00249999994412065f;
            break;
        }
    }
    float _S2688 = _S2678 + _S2685;
    float _S2689 = raw_losses_3[int(12)] / _S2679;
    for(;;)
    {
        float _S2690 = (F32_abs((_S2689)));
        if(_S2690 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2689 * _S2689 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2690 - 0.00249999994412065f;
            break;
        }
    }
    float _S2691 = _S2688 + _S2678;
    float _S2692 = raw_losses_3[int(13)] / _S2679;
    for(;;)
    {
        float _S2693 = (F32_abs((_S2692)));
        if(_S2693 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2692 * _S2692 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2693 - 0.00249999994412065f;
            break;
        }
    }
    float _S2694 = _S2691 + _S2678;
    float _S2695 = raw_losses_3[int(14)] / _S2679;
    for(;;)
    {
        float _S2696 = (F32_abs((_S2695)));
        if(_S2696 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2695 * _S2695 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2696 - 0.00249999994412065f;
            break;
        }
    }
    float _S2697 = _S2694 + _S2678;
    float _S2698 = raw_losses_3[int(15)] / _S2679;
    for(;;)
    {
        float _S2699 = (F32_abs((_S2698)));
        if(_S2699 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2698 * _S2698 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2699 - 0.00249999994412065f;
            break;
        }
    }
    float _S2700 = _S2697 + _S2678;
    float _S2701 = raw_losses_3[int(16)] / _S2679;
    for(;;)
    {
        float _S2702 = (F32_abs((_S2701)));
        if(_S2702 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2701 * _S2701 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2702 - 0.00249999994412065f;
            break;
        }
    }
    float _S2703 = _S2700 + _S2678;
    float _S2704 = raw_losses_3[int(17)] / _S2679;
    for(;;)
    {
        float _S2705 = (F32_abs((_S2704)));
        if(_S2705 < 0.00499999988824129f)
        {
            _S2678 = 0.5f * _S2704 * _S2704 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2678 = _S2705 - 0.00249999994412065f;
            break;
        }
    }
    float _S2706 = (_S2703 + _S2678) / 8.0f;
    float _S2707 = (raw_losses_3[int(18)] + raw_losses_3[int(19)] + raw_losses_3[int(20)] + raw_losses_3[int(21)] + raw_losses_3[int(22)]) / _S2682;
    losses_6[int(0)] = losses_6[int(0)] * loss_weights_1[int(0)];
    losses_6[int(1)] = losses_6[int(1)] * loss_weights_1[int(1)];
    losses_6[int(2)] = losses_6[int(2)] * loss_weights_1[int(2)];
    losses_6[int(3)] = losses_6[int(3)] * loss_weights_1[int(3)];
    losses_6[int(4)] = _S2706 * loss_weights_1[int(4)];
    losses_6[int(5)] = _S2707 * loss_weights_1[int(5)];
    *_S2677 = losses_6;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_15)
{
    FixedArray<float, 22>  _S2708 = dpraw_losses_1->primal_0;
    float _S2709 = float(num_cameras_2);
    float _S2710 = dpraw_losses_1->primal_0[int(0)] / _S2709;
    bool _S2711 = (s_primal_ctx_abs_0(_S2710)) < 0.10000000149011612f;
    float _S2712;
    if(_S2711)
    {
        _S2712 = 0.5f * _S2710;
    }
    else
    {
        _S2712 = 0.0f;
    }
    float _S2713 = 3.0f * _S2709;
    float _S2714 = 9.0f * _S2709;
    float _S2715 = 5.0f * _S2709;
    float _S2716 = _S2708[int(10)] / _S2709;
    bool _S2717 = (s_primal_ctx_abs_0(_S2716)) < 0.00499999988824129f;
    float _S2718;
    if(_S2717)
    {
        _S2718 = 0.5f * _S2716;
    }
    else
    {
        _S2718 = 0.0f;
    }
    float _S2719 = _S2708[int(11)] / _S2709;
    bool _S2720 = (s_primal_ctx_abs_0(_S2719)) < 0.00499999988824129f;
    float _S2721;
    if(_S2720)
    {
        _S2721 = 0.5f * _S2719;
    }
    else
    {
        _S2721 = 0.0f;
    }
    float _S2722 = _S2708[int(12)] / _S2709;
    bool _S2723 = (s_primal_ctx_abs_0(_S2722)) < 0.00499999988824129f;
    float _S2724;
    if(_S2723)
    {
        _S2724 = 0.5f * _S2722;
    }
    else
    {
        _S2724 = 0.0f;
    }
    float _S2725 = _S2708[int(13)] / _S2709;
    bool _S2726 = (s_primal_ctx_abs_0(_S2725)) < 0.00499999988824129f;
    float _S2727;
    if(_S2726)
    {
        _S2727 = 0.5f * _S2725;
    }
    else
    {
        _S2727 = 0.0f;
    }
    float _S2728 = _S2708[int(14)] / _S2709;
    bool _S2729 = (s_primal_ctx_abs_0(_S2728)) < 0.00499999988824129f;
    float _S2730;
    if(_S2729)
    {
        _S2730 = 0.5f * _S2728;
    }
    else
    {
        _S2730 = 0.0f;
    }
    float _S2731 = _S2708[int(15)] / _S2709;
    bool _S2732 = (s_primal_ctx_abs_0(_S2731)) < 0.00499999988824129f;
    float _S2733;
    if(_S2732)
    {
        _S2733 = 0.5f * _S2731;
    }
    else
    {
        _S2733 = 0.0f;
    }
    float _S2734 = _S2708[int(16)] / _S2709;
    bool _S2735 = (s_primal_ctx_abs_0(_S2734)) < 0.00499999988824129f;
    float _S2736;
    if(_S2735)
    {
        _S2736 = 0.5f * _S2734;
    }
    else
    {
        _S2736 = 0.0f;
    }
    float _S2737 = _S2708[int(17)] / _S2709;
    bool _S2738 = (s_primal_ctx_abs_0(_S2737)) < 0.00499999988824129f;
    float _S2739;
    if(_S2738)
    {
        _S2739 = 0.5f * _S2737;
    }
    else
    {
        _S2739 = 0.0f;
    }
    float _S2740 = (*loss_weights_2)[int(3)] * (*_s_dOut_15)[int(3)];
    float _S2741 = (*loss_weights_2)[int(2)] * (*_s_dOut_15)[int(2)];
    float _S2742 = (*loss_weights_2)[int(1)] * (*_s_dOut_15)[int(1)];
    float _S2743 = (*loss_weights_2)[int(0)] * (*_s_dOut_15)[int(0)];
    float _S2744 = (*loss_weights_2)[int(5)] * (*_s_dOut_15)[int(5)] / (4.0f * _S2709);
    float _S2745 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_15)[int(4)]);
    FixedArray<float, 22>  _S2746;
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
    _S2746[int(21)] = _S2744;
    _S2746[int(20)] = _S2744;
    _S2746[int(19)] = _S2744;
    _S2746[int(18)] = _S2744;
    float _S2747 = _S2746[int(0)];
    float _S2748 = _S2746[int(1)];
    float _S2749 = _S2746[int(2)];
    float _S2750 = _S2746[int(3)];
    float _S2751 = _S2746[int(4)];
    float _S2752 = _S2746[int(5)];
    float _S2753 = _S2746[int(6)];
    float _S2754 = _S2746[int(7)];
    float _S2755 = _S2746[int(8)];
    float _S2756 = _S2746[int(9)];
    float _S2757 = _S2746[int(10)];
    float _S2758 = _S2746[int(11)];
    float _S2759 = _S2746[int(12)];
    float _S2760 = _S2746[int(13)];
    float _S2761 = _S2746[int(14)];
    float _S2762 = _S2746[int(15)];
    float _S2763 = _S2746[int(16)];
    float _S2764 = _S2746[int(17)];
    float _S2765 = _S2746[int(18)];
    float _S2766 = _S2746[int(19)];
    float _S2767 = _S2746[int(20)];
    float _S2768 = _S2746[int(21)];
    float _S2769;
    if(_S2738)
    {
        float _S2770 = 200.0f * _S2745;
        float _S2771 = _S2739 * _S2770 + 0.5f * (_S2737 * _S2770);
        _S2739 = 0.0f;
        _S2769 = _S2771;
    }
    else
    {
        _S2739 = _S2745;
        _S2769 = 0.0f;
    }
    DiffPair_float_0 _S2772;
    (&_S2772)->primal_0 = _S2737;
    (&_S2772)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2772, _S2739);
    float _S2773 = (_S2772.differential_0 + _S2769) / _S2709;
    FixedArray<float, 22>  _S2774;
    _S2774[int(0)] = 0.0f;
    _S2774[int(1)] = 0.0f;
    _S2774[int(2)] = 0.0f;
    _S2774[int(3)] = 0.0f;
    _S2774[int(4)] = 0.0f;
    _S2774[int(5)] = 0.0f;
    _S2774[int(6)] = 0.0f;
    _S2774[int(7)] = 0.0f;
    _S2774[int(8)] = 0.0f;
    _S2774[int(9)] = 0.0f;
    _S2774[int(10)] = 0.0f;
    _S2774[int(11)] = 0.0f;
    _S2774[int(12)] = 0.0f;
    _S2774[int(13)] = 0.0f;
    _S2774[int(14)] = 0.0f;
    _S2774[int(15)] = 0.0f;
    _S2774[int(16)] = 0.0f;
    _S2774[int(17)] = 0.0f;
    _S2774[int(18)] = 0.0f;
    _S2774[int(19)] = 0.0f;
    _S2774[int(20)] = 0.0f;
    _S2774[int(21)] = 0.0f;
    _S2774[int(17)] = _S2773;
    float _S2775 = _S2747 + _S2774[int(0)];
    float _S2776 = _S2748 + _S2774[int(1)];
    float _S2777 = _S2749 + _S2774[int(2)];
    float _S2778 = _S2750 + _S2774[int(3)];
    float _S2779 = _S2751 + _S2774[int(4)];
    float _S2780 = _S2752 + _S2774[int(5)];
    float _S2781 = _S2753 + _S2774[int(6)];
    float _S2782 = _S2754 + _S2774[int(7)];
    float _S2783 = _S2755 + _S2774[int(8)];
    float _S2784 = _S2756 + _S2774[int(9)];
    float _S2785 = _S2757 + _S2774[int(10)];
    float _S2786 = _S2758 + _S2774[int(11)];
    float _S2787 = _S2759 + _S2774[int(12)];
    float _S2788 = _S2760 + _S2774[int(13)];
    float _S2789 = _S2761 + _S2774[int(14)];
    float _S2790 = _S2762 + _S2774[int(15)];
    float _S2791 = _S2763 + _S2774[int(16)];
    float _S2792 = _S2764 + _S2774[int(17)];
    float _S2793 = _S2765 + _S2774[int(18)];
    float _S2794 = _S2766 + _S2774[int(19)];
    float _S2795 = _S2767 + _S2774[int(20)];
    float _S2796 = _S2768 + _S2774[int(21)];
    if(_S2735)
    {
        float _S2797 = 200.0f * _S2745;
        float _S2798 = _S2736 * _S2797 + 0.5f * (_S2734 * _S2797);
        _S2736 = 0.0f;
        _S2739 = _S2798;
    }
    else
    {
        _S2736 = _S2745;
        _S2739 = 0.0f;
    }
    DiffPair_float_0 _S2799;
    (&_S2799)->primal_0 = _S2734;
    (&_S2799)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2799, _S2736);
    float _S2800 = (_S2799.differential_0 + _S2739) / _S2709;
    FixedArray<float, 22>  _S2801;
    _S2801[int(0)] = 0.0f;
    _S2801[int(1)] = 0.0f;
    _S2801[int(2)] = 0.0f;
    _S2801[int(3)] = 0.0f;
    _S2801[int(4)] = 0.0f;
    _S2801[int(5)] = 0.0f;
    _S2801[int(6)] = 0.0f;
    _S2801[int(7)] = 0.0f;
    _S2801[int(8)] = 0.0f;
    _S2801[int(9)] = 0.0f;
    _S2801[int(10)] = 0.0f;
    _S2801[int(11)] = 0.0f;
    _S2801[int(12)] = 0.0f;
    _S2801[int(13)] = 0.0f;
    _S2801[int(14)] = 0.0f;
    _S2801[int(15)] = 0.0f;
    _S2801[int(16)] = 0.0f;
    _S2801[int(17)] = 0.0f;
    _S2801[int(18)] = 0.0f;
    _S2801[int(19)] = 0.0f;
    _S2801[int(20)] = 0.0f;
    _S2801[int(21)] = 0.0f;
    _S2801[int(16)] = _S2800;
    float _S2802 = _S2775 + _S2801[int(0)];
    float _S2803 = _S2776 + _S2801[int(1)];
    float _S2804 = _S2777 + _S2801[int(2)];
    float _S2805 = _S2778 + _S2801[int(3)];
    float _S2806 = _S2779 + _S2801[int(4)];
    float _S2807 = _S2780 + _S2801[int(5)];
    float _S2808 = _S2781 + _S2801[int(6)];
    float _S2809 = _S2782 + _S2801[int(7)];
    float _S2810 = _S2783 + _S2801[int(8)];
    float _S2811 = _S2784 + _S2801[int(9)];
    float _S2812 = _S2785 + _S2801[int(10)];
    float _S2813 = _S2786 + _S2801[int(11)];
    float _S2814 = _S2787 + _S2801[int(12)];
    float _S2815 = _S2788 + _S2801[int(13)];
    float _S2816 = _S2789 + _S2801[int(14)];
    float _S2817 = _S2790 + _S2801[int(15)];
    float _S2818 = _S2791 + _S2801[int(16)];
    float _S2819 = _S2792 + _S2801[int(17)];
    float _S2820 = _S2793 + _S2801[int(18)];
    float _S2821 = _S2794 + _S2801[int(19)];
    float _S2822 = _S2795 + _S2801[int(20)];
    float _S2823 = _S2796 + _S2801[int(21)];
    if(_S2732)
    {
        float _S2824 = 200.0f * _S2745;
        float _S2825 = _S2733 * _S2824 + 0.5f * (_S2731 * _S2824);
        _S2733 = 0.0f;
        _S2736 = _S2825;
    }
    else
    {
        _S2733 = _S2745;
        _S2736 = 0.0f;
    }
    DiffPair_float_0 _S2826;
    (&_S2826)->primal_0 = _S2731;
    (&_S2826)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2826, _S2733);
    float _S2827 = (_S2826.differential_0 + _S2736) / _S2709;
    FixedArray<float, 22>  _S2828;
    _S2828[int(0)] = 0.0f;
    _S2828[int(1)] = 0.0f;
    _S2828[int(2)] = 0.0f;
    _S2828[int(3)] = 0.0f;
    _S2828[int(4)] = 0.0f;
    _S2828[int(5)] = 0.0f;
    _S2828[int(6)] = 0.0f;
    _S2828[int(7)] = 0.0f;
    _S2828[int(8)] = 0.0f;
    _S2828[int(9)] = 0.0f;
    _S2828[int(10)] = 0.0f;
    _S2828[int(11)] = 0.0f;
    _S2828[int(12)] = 0.0f;
    _S2828[int(13)] = 0.0f;
    _S2828[int(14)] = 0.0f;
    _S2828[int(15)] = 0.0f;
    _S2828[int(16)] = 0.0f;
    _S2828[int(17)] = 0.0f;
    _S2828[int(18)] = 0.0f;
    _S2828[int(19)] = 0.0f;
    _S2828[int(20)] = 0.0f;
    _S2828[int(21)] = 0.0f;
    _S2828[int(15)] = _S2827;
    float _S2829 = _S2802 + _S2828[int(0)];
    float _S2830 = _S2803 + _S2828[int(1)];
    float _S2831 = _S2804 + _S2828[int(2)];
    float _S2832 = _S2805 + _S2828[int(3)];
    float _S2833 = _S2806 + _S2828[int(4)];
    float _S2834 = _S2807 + _S2828[int(5)];
    float _S2835 = _S2808 + _S2828[int(6)];
    float _S2836 = _S2809 + _S2828[int(7)];
    float _S2837 = _S2810 + _S2828[int(8)];
    float _S2838 = _S2811 + _S2828[int(9)];
    float _S2839 = _S2812 + _S2828[int(10)];
    float _S2840 = _S2813 + _S2828[int(11)];
    float _S2841 = _S2814 + _S2828[int(12)];
    float _S2842 = _S2815 + _S2828[int(13)];
    float _S2843 = _S2816 + _S2828[int(14)];
    float _S2844 = _S2817 + _S2828[int(15)];
    float _S2845 = _S2818 + _S2828[int(16)];
    float _S2846 = _S2819 + _S2828[int(17)];
    float _S2847 = _S2820 + _S2828[int(18)];
    float _S2848 = _S2821 + _S2828[int(19)];
    float _S2849 = _S2822 + _S2828[int(20)];
    float _S2850 = _S2823 + _S2828[int(21)];
    if(_S2729)
    {
        float _S2851 = 200.0f * _S2745;
        float _S2852 = _S2730 * _S2851 + 0.5f * (_S2728 * _S2851);
        _S2730 = 0.0f;
        _S2733 = _S2852;
    }
    else
    {
        _S2730 = _S2745;
        _S2733 = 0.0f;
    }
    DiffPair_float_0 _S2853;
    (&_S2853)->primal_0 = _S2728;
    (&_S2853)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2853, _S2730);
    float _S2854 = (_S2853.differential_0 + _S2733) / _S2709;
    FixedArray<float, 22>  _S2855;
    _S2855[int(0)] = 0.0f;
    _S2855[int(1)] = 0.0f;
    _S2855[int(2)] = 0.0f;
    _S2855[int(3)] = 0.0f;
    _S2855[int(4)] = 0.0f;
    _S2855[int(5)] = 0.0f;
    _S2855[int(6)] = 0.0f;
    _S2855[int(7)] = 0.0f;
    _S2855[int(8)] = 0.0f;
    _S2855[int(9)] = 0.0f;
    _S2855[int(10)] = 0.0f;
    _S2855[int(11)] = 0.0f;
    _S2855[int(12)] = 0.0f;
    _S2855[int(13)] = 0.0f;
    _S2855[int(14)] = 0.0f;
    _S2855[int(15)] = 0.0f;
    _S2855[int(16)] = 0.0f;
    _S2855[int(17)] = 0.0f;
    _S2855[int(18)] = 0.0f;
    _S2855[int(19)] = 0.0f;
    _S2855[int(20)] = 0.0f;
    _S2855[int(21)] = 0.0f;
    _S2855[int(14)] = _S2854;
    float _S2856 = _S2829 + _S2855[int(0)];
    float _S2857 = _S2830 + _S2855[int(1)];
    float _S2858 = _S2831 + _S2855[int(2)];
    float _S2859 = _S2832 + _S2855[int(3)];
    float _S2860 = _S2833 + _S2855[int(4)];
    float _S2861 = _S2834 + _S2855[int(5)];
    float _S2862 = _S2835 + _S2855[int(6)];
    float _S2863 = _S2836 + _S2855[int(7)];
    float _S2864 = _S2837 + _S2855[int(8)];
    float _S2865 = _S2838 + _S2855[int(9)];
    float _S2866 = _S2839 + _S2855[int(10)];
    float _S2867 = _S2840 + _S2855[int(11)];
    float _S2868 = _S2841 + _S2855[int(12)];
    float _S2869 = _S2842 + _S2855[int(13)];
    float _S2870 = _S2843 + _S2855[int(14)];
    float _S2871 = _S2844 + _S2855[int(15)];
    float _S2872 = _S2845 + _S2855[int(16)];
    float _S2873 = _S2846 + _S2855[int(17)];
    float _S2874 = _S2847 + _S2855[int(18)];
    float _S2875 = _S2848 + _S2855[int(19)];
    float _S2876 = _S2849 + _S2855[int(20)];
    float _S2877 = _S2850 + _S2855[int(21)];
    if(_S2726)
    {
        float _S2878 = 200.0f * _S2745;
        float _S2879 = _S2727 * _S2878 + 0.5f * (_S2725 * _S2878);
        _S2727 = 0.0f;
        _S2730 = _S2879;
    }
    else
    {
        _S2727 = _S2745;
        _S2730 = 0.0f;
    }
    DiffPair_float_0 _S2880;
    (&_S2880)->primal_0 = _S2725;
    (&_S2880)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2880, _S2727);
    float _S2881 = (_S2880.differential_0 + _S2730) / _S2709;
    FixedArray<float, 22>  _S2882;
    _S2882[int(0)] = 0.0f;
    _S2882[int(1)] = 0.0f;
    _S2882[int(2)] = 0.0f;
    _S2882[int(3)] = 0.0f;
    _S2882[int(4)] = 0.0f;
    _S2882[int(5)] = 0.0f;
    _S2882[int(6)] = 0.0f;
    _S2882[int(7)] = 0.0f;
    _S2882[int(8)] = 0.0f;
    _S2882[int(9)] = 0.0f;
    _S2882[int(10)] = 0.0f;
    _S2882[int(11)] = 0.0f;
    _S2882[int(12)] = 0.0f;
    _S2882[int(13)] = 0.0f;
    _S2882[int(14)] = 0.0f;
    _S2882[int(15)] = 0.0f;
    _S2882[int(16)] = 0.0f;
    _S2882[int(17)] = 0.0f;
    _S2882[int(18)] = 0.0f;
    _S2882[int(19)] = 0.0f;
    _S2882[int(20)] = 0.0f;
    _S2882[int(21)] = 0.0f;
    _S2882[int(13)] = _S2881;
    float _S2883 = _S2856 + _S2882[int(0)];
    float _S2884 = _S2857 + _S2882[int(1)];
    float _S2885 = _S2858 + _S2882[int(2)];
    float _S2886 = _S2859 + _S2882[int(3)];
    float _S2887 = _S2860 + _S2882[int(4)];
    float _S2888 = _S2861 + _S2882[int(5)];
    float _S2889 = _S2862 + _S2882[int(6)];
    float _S2890 = _S2863 + _S2882[int(7)];
    float _S2891 = _S2864 + _S2882[int(8)];
    float _S2892 = _S2865 + _S2882[int(9)];
    float _S2893 = _S2866 + _S2882[int(10)];
    float _S2894 = _S2867 + _S2882[int(11)];
    float _S2895 = _S2868 + _S2882[int(12)];
    float _S2896 = _S2869 + _S2882[int(13)];
    float _S2897 = _S2870 + _S2882[int(14)];
    float _S2898 = _S2871 + _S2882[int(15)];
    float _S2899 = _S2872 + _S2882[int(16)];
    float _S2900 = _S2873 + _S2882[int(17)];
    float _S2901 = _S2874 + _S2882[int(18)];
    float _S2902 = _S2875 + _S2882[int(19)];
    float _S2903 = _S2876 + _S2882[int(20)];
    float _S2904 = _S2877 + _S2882[int(21)];
    if(_S2723)
    {
        float _S2905 = 200.0f * _S2745;
        float _S2906 = _S2724 * _S2905 + 0.5f * (_S2722 * _S2905);
        _S2724 = 0.0f;
        _S2727 = _S2906;
    }
    else
    {
        _S2724 = _S2745;
        _S2727 = 0.0f;
    }
    DiffPair_float_0 _S2907;
    (&_S2907)->primal_0 = _S2722;
    (&_S2907)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2907, _S2724);
    float _S2908 = (_S2907.differential_0 + _S2727) / _S2709;
    FixedArray<float, 22>  _S2909;
    _S2909[int(0)] = 0.0f;
    _S2909[int(1)] = 0.0f;
    _S2909[int(2)] = 0.0f;
    _S2909[int(3)] = 0.0f;
    _S2909[int(4)] = 0.0f;
    _S2909[int(5)] = 0.0f;
    _S2909[int(6)] = 0.0f;
    _S2909[int(7)] = 0.0f;
    _S2909[int(8)] = 0.0f;
    _S2909[int(9)] = 0.0f;
    _S2909[int(10)] = 0.0f;
    _S2909[int(11)] = 0.0f;
    _S2909[int(12)] = 0.0f;
    _S2909[int(13)] = 0.0f;
    _S2909[int(14)] = 0.0f;
    _S2909[int(15)] = 0.0f;
    _S2909[int(16)] = 0.0f;
    _S2909[int(17)] = 0.0f;
    _S2909[int(18)] = 0.0f;
    _S2909[int(19)] = 0.0f;
    _S2909[int(20)] = 0.0f;
    _S2909[int(21)] = 0.0f;
    _S2909[int(12)] = _S2908;
    float _S2910 = _S2883 + _S2909[int(0)];
    float _S2911 = _S2884 + _S2909[int(1)];
    float _S2912 = _S2885 + _S2909[int(2)];
    float _S2913 = _S2886 + _S2909[int(3)];
    float _S2914 = _S2887 + _S2909[int(4)];
    float _S2915 = _S2888 + _S2909[int(5)];
    float _S2916 = _S2889 + _S2909[int(6)];
    float _S2917 = _S2890 + _S2909[int(7)];
    float _S2918 = _S2891 + _S2909[int(8)];
    float _S2919 = _S2892 + _S2909[int(9)];
    float _S2920 = _S2893 + _S2909[int(10)];
    float _S2921 = _S2894 + _S2909[int(11)];
    float _S2922 = _S2895 + _S2909[int(12)];
    float _S2923 = _S2896 + _S2909[int(13)];
    float _S2924 = _S2897 + _S2909[int(14)];
    float _S2925 = _S2898 + _S2909[int(15)];
    float _S2926 = _S2899 + _S2909[int(16)];
    float _S2927 = _S2900 + _S2909[int(17)];
    float _S2928 = _S2901 + _S2909[int(18)];
    float _S2929 = _S2902 + _S2909[int(19)];
    float _S2930 = _S2903 + _S2909[int(20)];
    float _S2931 = _S2904 + _S2909[int(21)];
    if(_S2720)
    {
        float _S2932 = 200.0f * _S2745;
        float _S2933 = _S2721 * _S2932 + 0.5f * (_S2719 * _S2932);
        _S2721 = 0.0f;
        _S2724 = _S2933;
    }
    else
    {
        _S2721 = _S2745;
        _S2724 = 0.0f;
    }
    DiffPair_float_0 _S2934;
    (&_S2934)->primal_0 = _S2719;
    (&_S2934)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2934, _S2721);
    float _S2935 = (_S2934.differential_0 + _S2724) / _S2709;
    FixedArray<float, 22>  _S2936;
    _S2936[int(0)] = 0.0f;
    _S2936[int(1)] = 0.0f;
    _S2936[int(2)] = 0.0f;
    _S2936[int(3)] = 0.0f;
    _S2936[int(4)] = 0.0f;
    _S2936[int(5)] = 0.0f;
    _S2936[int(6)] = 0.0f;
    _S2936[int(7)] = 0.0f;
    _S2936[int(8)] = 0.0f;
    _S2936[int(9)] = 0.0f;
    _S2936[int(10)] = 0.0f;
    _S2936[int(11)] = 0.0f;
    _S2936[int(12)] = 0.0f;
    _S2936[int(13)] = 0.0f;
    _S2936[int(14)] = 0.0f;
    _S2936[int(15)] = 0.0f;
    _S2936[int(16)] = 0.0f;
    _S2936[int(17)] = 0.0f;
    _S2936[int(18)] = 0.0f;
    _S2936[int(19)] = 0.0f;
    _S2936[int(20)] = 0.0f;
    _S2936[int(21)] = 0.0f;
    _S2936[int(11)] = _S2935;
    float _S2937 = _S2910 + _S2936[int(0)];
    float _S2938 = _S2911 + _S2936[int(1)];
    float _S2939 = _S2912 + _S2936[int(2)];
    float _S2940 = _S2913 + _S2936[int(3)];
    float _S2941 = _S2914 + _S2936[int(4)];
    float _S2942 = _S2915 + _S2936[int(5)];
    float _S2943 = _S2916 + _S2936[int(6)];
    float _S2944 = _S2917 + _S2936[int(7)];
    float _S2945 = _S2918 + _S2936[int(8)];
    float _S2946 = _S2919 + _S2936[int(9)];
    float _S2947 = _S2920 + _S2936[int(10)];
    float _S2948 = _S2921 + _S2936[int(11)];
    float _S2949 = _S2922 + _S2936[int(12)];
    float _S2950 = _S2923 + _S2936[int(13)];
    float _S2951 = _S2924 + _S2936[int(14)];
    float _S2952 = _S2925 + _S2936[int(15)];
    float _S2953 = _S2926 + _S2936[int(16)];
    float _S2954 = _S2927 + _S2936[int(17)];
    float _S2955 = _S2928 + _S2936[int(18)];
    float _S2956 = _S2929 + _S2936[int(19)];
    float _S2957 = _S2930 + _S2936[int(20)];
    float _S2958 = _S2931 + _S2936[int(21)];
    if(_S2717)
    {
        float _S2959 = 200.0f * _S2745;
        float _S2960 = _S2718 * _S2959 + 0.5f * (_S2716 * _S2959);
        _S2718 = 0.0f;
        _S2721 = _S2960;
    }
    else
    {
        _S2718 = _S2745;
        _S2721 = 0.0f;
    }
    DiffPair_float_0 _S2961;
    (&_S2961)->primal_0 = _S2716;
    (&_S2961)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2961, _S2718);
    float _S2962 = (_S2961.differential_0 + _S2721) / _S2709;
    float _S2963 = _S2740 / _S2715;
    float _S2964 = _S2741 / _S2714;
    float _S2965 = _S2742 / _S2713;
    FixedArray<float, 22>  _S2966;
    _S2966[int(0)] = 0.0f;
    _S2966[int(1)] = 0.0f;
    _S2966[int(2)] = 0.0f;
    _S2966[int(3)] = 0.0f;
    _S2966[int(4)] = 0.0f;
    _S2966[int(5)] = 0.0f;
    _S2966[int(6)] = 0.0f;
    _S2966[int(7)] = 0.0f;
    _S2966[int(8)] = 0.0f;
    _S2966[int(9)] = 0.0f;
    _S2966[int(10)] = 0.0f;
    _S2966[int(11)] = 0.0f;
    _S2966[int(12)] = 0.0f;
    _S2966[int(13)] = 0.0f;
    _S2966[int(14)] = 0.0f;
    _S2966[int(15)] = 0.0f;
    _S2966[int(16)] = 0.0f;
    _S2966[int(17)] = 0.0f;
    _S2966[int(18)] = 0.0f;
    _S2966[int(19)] = 0.0f;
    _S2966[int(20)] = 0.0f;
    _S2966[int(21)] = 0.0f;
    _S2966[int(10)] = _S2962;
    _S2966[int(9)] = _S2963;
    _S2966[int(8)] = _S2963;
    _S2966[int(7)] = _S2963;
    _S2966[int(6)] = _S2963;
    _S2966[int(5)] = _S2963;
    _S2966[int(4)] = _S2964;
    _S2966[int(3)] = _S2964;
    _S2966[int(2)] = _S2964;
    _S2966[int(1)] = _S2965;
    float _S2967 = _S2937 + _S2966[int(0)];
    float _S2968 = _S2938 + _S2966[int(1)];
    float _S2969 = _S2939 + _S2966[int(2)];
    float _S2970 = _S2940 + _S2966[int(3)];
    float _S2971 = _S2941 + _S2966[int(4)];
    float _S2972 = _S2942 + _S2966[int(5)];
    float _S2973 = _S2943 + _S2966[int(6)];
    float _S2974 = _S2944 + _S2966[int(7)];
    float _S2975 = _S2945 + _S2966[int(8)];
    float _S2976 = _S2946 + _S2966[int(9)];
    float _S2977 = _S2947 + _S2966[int(10)];
    float _S2978 = _S2948 + _S2966[int(11)];
    float _S2979 = _S2949 + _S2966[int(12)];
    float _S2980 = _S2950 + _S2966[int(13)];
    float _S2981 = _S2951 + _S2966[int(14)];
    float _S2982 = _S2952 + _S2966[int(15)];
    float _S2983 = _S2953 + _S2966[int(16)];
    float _S2984 = _S2954 + _S2966[int(17)];
    float _S2985 = _S2955 + _S2966[int(18)];
    float _S2986 = _S2956 + _S2966[int(19)];
    float _S2987 = _S2957 + _S2966[int(20)];
    float _S2988 = _S2958 + _S2966[int(21)];
    if(_S2711)
    {
        float _S2989 = 10.0f * _S2743;
        float _S2990 = _S2712 * _S2989 + 0.5f * (_S2710 * _S2989);
        _S2712 = 0.0f;
        _S2718 = _S2990;
    }
    else
    {
        _S2712 = _S2743;
        _S2718 = 0.0f;
    }
    DiffPair_float_0 _S2991;
    (&_S2991)->primal_0 = _S2710;
    (&_S2991)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2991, _S2712);
    float _S2992 = (_S2991.differential_0 + _S2718) / _S2709;
    FixedArray<float, 22>  _S2993;
    _S2993[int(0)] = 0.0f;
    _S2993[int(1)] = 0.0f;
    _S2993[int(2)] = 0.0f;
    _S2993[int(3)] = 0.0f;
    _S2993[int(4)] = 0.0f;
    _S2993[int(5)] = 0.0f;
    _S2993[int(6)] = 0.0f;
    _S2993[int(7)] = 0.0f;
    _S2993[int(8)] = 0.0f;
    _S2993[int(9)] = 0.0f;
    _S2993[int(10)] = 0.0f;
    _S2993[int(11)] = 0.0f;
    _S2993[int(12)] = 0.0f;
    _S2993[int(13)] = 0.0f;
    _S2993[int(14)] = 0.0f;
    _S2993[int(15)] = 0.0f;
    _S2993[int(16)] = 0.0f;
    _S2993[int(17)] = 0.0f;
    _S2993[int(18)] = 0.0f;
    _S2993[int(19)] = 0.0f;
    _S2993[int(20)] = 0.0f;
    _S2993[int(21)] = 0.0f;
    _S2993[int(0)] = _S2992;
    FixedArray<float, 22>  _S2994 = {
        _S2967 + _S2993[int(0)], _S2968 + _S2993[int(1)], _S2969 + _S2993[int(2)], _S2970 + _S2993[int(3)], _S2971 + _S2993[int(4)], _S2972 + _S2993[int(5)], _S2973 + _S2993[int(6)], _S2974 + _S2993[int(7)], _S2975 + _S2993[int(8)], _S2976 + _S2993[int(9)], _S2977 + _S2993[int(10)], _S2978 + _S2993[int(11)], _S2979 + _S2993[int(12)], _S2980 + _S2993[int(13)], _S2981 + _S2993[int(14)], _S2982 + _S2993[int(15)], _S2983 + _S2993[int(16)], _S2984 + _S2993[int(17)], _S2985 + _S2993[int(18)], _S2986 + _S2993[int(19)], _S2987 + _S2993[int(20)], _S2988 + _S2993[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2994;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S2995, int _S2996, FixedArray<float, 6>  * _S2997, FixedArray<float, 6>  * _S2998)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2995, _S2996, _S2997, _S2998);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_4, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_4, FixedArray<float, 22>  * _S2999)
{
    FixedArray<float, 22>  _S3000 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_4;
    (&dp_raw_losses_1)->differential_0 = _S3000;
    FixedArray<float, 6>  _S3001 = loss_weights_3;
    FixedArray<float, 6>  _S3002 = grad_out_4;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_3, &_S3001, &_S3002);
    *_S2999 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_2, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_16)
{
    FixedArray<float, 23>  _S3003 = dpraw_losses_2->primal_0;
    float _S3004 = float(num_cameras_4);
    float _S3005 = dpraw_losses_2->primal_0[int(0)] / _S3004;
    bool _S3006 = (s_primal_ctx_abs_0(_S3005)) < 0.10000000149011612f;
    float _S3007;
    if(_S3006)
    {
        _S3007 = 0.5f * _S3005;
    }
    else
    {
        _S3007 = 0.0f;
    }
    float _S3008 = 3.0f * _S3004;
    float _S3009 = 9.0f * _S3004;
    float _S3010 = 5.0f * _S3004;
    float _S3011 = _S3003[int(10)] / _S3004;
    bool _S3012 = (s_primal_ctx_abs_0(_S3011)) < 0.00499999988824129f;
    float _S3013;
    if(_S3012)
    {
        _S3013 = 0.5f * _S3011;
    }
    else
    {
        _S3013 = 0.0f;
    }
    float _S3014 = _S3003[int(11)] / _S3004;
    bool _S3015 = (s_primal_ctx_abs_0(_S3014)) < 0.00499999988824129f;
    float _S3016;
    if(_S3015)
    {
        _S3016 = 0.5f * _S3014;
    }
    else
    {
        _S3016 = 0.0f;
    }
    float _S3017 = _S3003[int(12)] / _S3004;
    bool _S3018 = (s_primal_ctx_abs_0(_S3017)) < 0.00499999988824129f;
    float _S3019;
    if(_S3018)
    {
        _S3019 = 0.5f * _S3017;
    }
    else
    {
        _S3019 = 0.0f;
    }
    float _S3020 = _S3003[int(13)] / _S3004;
    bool _S3021 = (s_primal_ctx_abs_0(_S3020)) < 0.00499999988824129f;
    float _S3022;
    if(_S3021)
    {
        _S3022 = 0.5f * _S3020;
    }
    else
    {
        _S3022 = 0.0f;
    }
    float _S3023 = _S3003[int(14)] / _S3004;
    bool _S3024 = (s_primal_ctx_abs_0(_S3023)) < 0.00499999988824129f;
    float _S3025;
    if(_S3024)
    {
        _S3025 = 0.5f * _S3023;
    }
    else
    {
        _S3025 = 0.0f;
    }
    float _S3026 = _S3003[int(15)] / _S3004;
    bool _S3027 = (s_primal_ctx_abs_0(_S3026)) < 0.00499999988824129f;
    float _S3028;
    if(_S3027)
    {
        _S3028 = 0.5f * _S3026;
    }
    else
    {
        _S3028 = 0.0f;
    }
    float _S3029 = _S3003[int(16)] / _S3004;
    bool _S3030 = (s_primal_ctx_abs_0(_S3029)) < 0.00499999988824129f;
    float _S3031;
    if(_S3030)
    {
        _S3031 = 0.5f * _S3029;
    }
    else
    {
        _S3031 = 0.0f;
    }
    float _S3032 = _S3003[int(17)] / _S3004;
    bool _S3033 = (s_primal_ctx_abs_0(_S3032)) < 0.00499999988824129f;
    float _S3034;
    if(_S3033)
    {
        _S3034 = 0.5f * _S3032;
    }
    else
    {
        _S3034 = 0.0f;
    }
    float _S3035 = (*loss_weights_4)[int(3)] * (*_s_dOut_16)[int(3)];
    float _S3036 = (*loss_weights_4)[int(2)] * (*_s_dOut_16)[int(2)];
    float _S3037 = (*loss_weights_4)[int(1)] * (*_s_dOut_16)[int(1)];
    float _S3038 = (*loss_weights_4)[int(0)] * (*_s_dOut_16)[int(0)];
    float _S3039 = (*loss_weights_4)[int(5)] * (*_s_dOut_16)[int(5)] / _S3010;
    float _S3040 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_16)[int(4)]);
    FixedArray<float, 23>  _S3041;
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
    _S3041[int(22)] = 0.0f;
    _S3041[int(22)] = _S3039;
    _S3041[int(21)] = _S3039;
    _S3041[int(20)] = _S3039;
    _S3041[int(19)] = _S3039;
    _S3041[int(18)] = _S3039;
    float _S3042 = _S3041[int(0)];
    float _S3043 = _S3041[int(1)];
    float _S3044 = _S3041[int(2)];
    float _S3045 = _S3041[int(3)];
    float _S3046 = _S3041[int(4)];
    float _S3047 = _S3041[int(5)];
    float _S3048 = _S3041[int(6)];
    float _S3049 = _S3041[int(7)];
    float _S3050 = _S3041[int(8)];
    float _S3051 = _S3041[int(9)];
    float _S3052 = _S3041[int(10)];
    float _S3053 = _S3041[int(11)];
    float _S3054 = _S3041[int(12)];
    float _S3055 = _S3041[int(13)];
    float _S3056 = _S3041[int(14)];
    float _S3057 = _S3041[int(15)];
    float _S3058 = _S3041[int(16)];
    float _S3059 = _S3041[int(17)];
    float _S3060 = _S3041[int(18)];
    float _S3061 = _S3041[int(19)];
    float _S3062 = _S3041[int(20)];
    float _S3063 = _S3041[int(21)];
    float _S3064 = _S3041[int(22)];
    float _S3065;
    if(_S3033)
    {
        float _S3066 = 200.0f * _S3040;
        float _S3067 = _S3034 * _S3066 + 0.5f * (_S3032 * _S3066);
        _S3034 = 0.0f;
        _S3065 = _S3067;
    }
    else
    {
        _S3034 = _S3040;
        _S3065 = 0.0f;
    }
    DiffPair_float_0 _S3068;
    (&_S3068)->primal_0 = _S3032;
    (&_S3068)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3068, _S3034);
    float _S3069 = (_S3068.differential_0 + _S3065) / _S3004;
    FixedArray<float, 23>  _S3070;
    _S3070[int(0)] = 0.0f;
    _S3070[int(1)] = 0.0f;
    _S3070[int(2)] = 0.0f;
    _S3070[int(3)] = 0.0f;
    _S3070[int(4)] = 0.0f;
    _S3070[int(5)] = 0.0f;
    _S3070[int(6)] = 0.0f;
    _S3070[int(7)] = 0.0f;
    _S3070[int(8)] = 0.0f;
    _S3070[int(9)] = 0.0f;
    _S3070[int(10)] = 0.0f;
    _S3070[int(11)] = 0.0f;
    _S3070[int(12)] = 0.0f;
    _S3070[int(13)] = 0.0f;
    _S3070[int(14)] = 0.0f;
    _S3070[int(15)] = 0.0f;
    _S3070[int(16)] = 0.0f;
    _S3070[int(17)] = 0.0f;
    _S3070[int(18)] = 0.0f;
    _S3070[int(19)] = 0.0f;
    _S3070[int(20)] = 0.0f;
    _S3070[int(21)] = 0.0f;
    _S3070[int(22)] = 0.0f;
    _S3070[int(17)] = _S3069;
    float _S3071 = _S3042 + _S3070[int(0)];
    float _S3072 = _S3043 + _S3070[int(1)];
    float _S3073 = _S3044 + _S3070[int(2)];
    float _S3074 = _S3045 + _S3070[int(3)];
    float _S3075 = _S3046 + _S3070[int(4)];
    float _S3076 = _S3047 + _S3070[int(5)];
    float _S3077 = _S3048 + _S3070[int(6)];
    float _S3078 = _S3049 + _S3070[int(7)];
    float _S3079 = _S3050 + _S3070[int(8)];
    float _S3080 = _S3051 + _S3070[int(9)];
    float _S3081 = _S3052 + _S3070[int(10)];
    float _S3082 = _S3053 + _S3070[int(11)];
    float _S3083 = _S3054 + _S3070[int(12)];
    float _S3084 = _S3055 + _S3070[int(13)];
    float _S3085 = _S3056 + _S3070[int(14)];
    float _S3086 = _S3057 + _S3070[int(15)];
    float _S3087 = _S3058 + _S3070[int(16)];
    float _S3088 = _S3059 + _S3070[int(17)];
    float _S3089 = _S3060 + _S3070[int(18)];
    float _S3090 = _S3061 + _S3070[int(19)];
    float _S3091 = _S3062 + _S3070[int(20)];
    float _S3092 = _S3063 + _S3070[int(21)];
    float _S3093 = _S3064 + _S3070[int(22)];
    if(_S3030)
    {
        float _S3094 = 200.0f * _S3040;
        float _S3095 = _S3031 * _S3094 + 0.5f * (_S3029 * _S3094);
        _S3031 = 0.0f;
        _S3034 = _S3095;
    }
    else
    {
        _S3031 = _S3040;
        _S3034 = 0.0f;
    }
    DiffPair_float_0 _S3096;
    (&_S3096)->primal_0 = _S3029;
    (&_S3096)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3096, _S3031);
    float _S3097 = (_S3096.differential_0 + _S3034) / _S3004;
    FixedArray<float, 23>  _S3098;
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
    _S3098[int(22)] = 0.0f;
    _S3098[int(16)] = _S3097;
    float _S3099 = _S3071 + _S3098[int(0)];
    float _S3100 = _S3072 + _S3098[int(1)];
    float _S3101 = _S3073 + _S3098[int(2)];
    float _S3102 = _S3074 + _S3098[int(3)];
    float _S3103 = _S3075 + _S3098[int(4)];
    float _S3104 = _S3076 + _S3098[int(5)];
    float _S3105 = _S3077 + _S3098[int(6)];
    float _S3106 = _S3078 + _S3098[int(7)];
    float _S3107 = _S3079 + _S3098[int(8)];
    float _S3108 = _S3080 + _S3098[int(9)];
    float _S3109 = _S3081 + _S3098[int(10)];
    float _S3110 = _S3082 + _S3098[int(11)];
    float _S3111 = _S3083 + _S3098[int(12)];
    float _S3112 = _S3084 + _S3098[int(13)];
    float _S3113 = _S3085 + _S3098[int(14)];
    float _S3114 = _S3086 + _S3098[int(15)];
    float _S3115 = _S3087 + _S3098[int(16)];
    float _S3116 = _S3088 + _S3098[int(17)];
    float _S3117 = _S3089 + _S3098[int(18)];
    float _S3118 = _S3090 + _S3098[int(19)];
    float _S3119 = _S3091 + _S3098[int(20)];
    float _S3120 = _S3092 + _S3098[int(21)];
    float _S3121 = _S3093 + _S3098[int(22)];
    if(_S3027)
    {
        float _S3122 = 200.0f * _S3040;
        float _S3123 = _S3028 * _S3122 + 0.5f * (_S3026 * _S3122);
        _S3028 = 0.0f;
        _S3031 = _S3123;
    }
    else
    {
        _S3028 = _S3040;
        _S3031 = 0.0f;
    }
    DiffPair_float_0 _S3124;
    (&_S3124)->primal_0 = _S3026;
    (&_S3124)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3124, _S3028);
    float _S3125 = (_S3124.differential_0 + _S3031) / _S3004;
    FixedArray<float, 23>  _S3126;
    _S3126[int(0)] = 0.0f;
    _S3126[int(1)] = 0.0f;
    _S3126[int(2)] = 0.0f;
    _S3126[int(3)] = 0.0f;
    _S3126[int(4)] = 0.0f;
    _S3126[int(5)] = 0.0f;
    _S3126[int(6)] = 0.0f;
    _S3126[int(7)] = 0.0f;
    _S3126[int(8)] = 0.0f;
    _S3126[int(9)] = 0.0f;
    _S3126[int(10)] = 0.0f;
    _S3126[int(11)] = 0.0f;
    _S3126[int(12)] = 0.0f;
    _S3126[int(13)] = 0.0f;
    _S3126[int(14)] = 0.0f;
    _S3126[int(15)] = 0.0f;
    _S3126[int(16)] = 0.0f;
    _S3126[int(17)] = 0.0f;
    _S3126[int(18)] = 0.0f;
    _S3126[int(19)] = 0.0f;
    _S3126[int(20)] = 0.0f;
    _S3126[int(21)] = 0.0f;
    _S3126[int(22)] = 0.0f;
    _S3126[int(15)] = _S3125;
    float _S3127 = _S3099 + _S3126[int(0)];
    float _S3128 = _S3100 + _S3126[int(1)];
    float _S3129 = _S3101 + _S3126[int(2)];
    float _S3130 = _S3102 + _S3126[int(3)];
    float _S3131 = _S3103 + _S3126[int(4)];
    float _S3132 = _S3104 + _S3126[int(5)];
    float _S3133 = _S3105 + _S3126[int(6)];
    float _S3134 = _S3106 + _S3126[int(7)];
    float _S3135 = _S3107 + _S3126[int(8)];
    float _S3136 = _S3108 + _S3126[int(9)];
    float _S3137 = _S3109 + _S3126[int(10)];
    float _S3138 = _S3110 + _S3126[int(11)];
    float _S3139 = _S3111 + _S3126[int(12)];
    float _S3140 = _S3112 + _S3126[int(13)];
    float _S3141 = _S3113 + _S3126[int(14)];
    float _S3142 = _S3114 + _S3126[int(15)];
    float _S3143 = _S3115 + _S3126[int(16)];
    float _S3144 = _S3116 + _S3126[int(17)];
    float _S3145 = _S3117 + _S3126[int(18)];
    float _S3146 = _S3118 + _S3126[int(19)];
    float _S3147 = _S3119 + _S3126[int(20)];
    float _S3148 = _S3120 + _S3126[int(21)];
    float _S3149 = _S3121 + _S3126[int(22)];
    if(_S3024)
    {
        float _S3150 = 200.0f * _S3040;
        float _S3151 = _S3025 * _S3150 + 0.5f * (_S3023 * _S3150);
        _S3025 = 0.0f;
        _S3028 = _S3151;
    }
    else
    {
        _S3025 = _S3040;
        _S3028 = 0.0f;
    }
    DiffPair_float_0 _S3152;
    (&_S3152)->primal_0 = _S3023;
    (&_S3152)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3152, _S3025);
    float _S3153 = (_S3152.differential_0 + _S3028) / _S3004;
    FixedArray<float, 23>  _S3154;
    _S3154[int(0)] = 0.0f;
    _S3154[int(1)] = 0.0f;
    _S3154[int(2)] = 0.0f;
    _S3154[int(3)] = 0.0f;
    _S3154[int(4)] = 0.0f;
    _S3154[int(5)] = 0.0f;
    _S3154[int(6)] = 0.0f;
    _S3154[int(7)] = 0.0f;
    _S3154[int(8)] = 0.0f;
    _S3154[int(9)] = 0.0f;
    _S3154[int(10)] = 0.0f;
    _S3154[int(11)] = 0.0f;
    _S3154[int(12)] = 0.0f;
    _S3154[int(13)] = 0.0f;
    _S3154[int(14)] = 0.0f;
    _S3154[int(15)] = 0.0f;
    _S3154[int(16)] = 0.0f;
    _S3154[int(17)] = 0.0f;
    _S3154[int(18)] = 0.0f;
    _S3154[int(19)] = 0.0f;
    _S3154[int(20)] = 0.0f;
    _S3154[int(21)] = 0.0f;
    _S3154[int(22)] = 0.0f;
    _S3154[int(14)] = _S3153;
    float _S3155 = _S3127 + _S3154[int(0)];
    float _S3156 = _S3128 + _S3154[int(1)];
    float _S3157 = _S3129 + _S3154[int(2)];
    float _S3158 = _S3130 + _S3154[int(3)];
    float _S3159 = _S3131 + _S3154[int(4)];
    float _S3160 = _S3132 + _S3154[int(5)];
    float _S3161 = _S3133 + _S3154[int(6)];
    float _S3162 = _S3134 + _S3154[int(7)];
    float _S3163 = _S3135 + _S3154[int(8)];
    float _S3164 = _S3136 + _S3154[int(9)];
    float _S3165 = _S3137 + _S3154[int(10)];
    float _S3166 = _S3138 + _S3154[int(11)];
    float _S3167 = _S3139 + _S3154[int(12)];
    float _S3168 = _S3140 + _S3154[int(13)];
    float _S3169 = _S3141 + _S3154[int(14)];
    float _S3170 = _S3142 + _S3154[int(15)];
    float _S3171 = _S3143 + _S3154[int(16)];
    float _S3172 = _S3144 + _S3154[int(17)];
    float _S3173 = _S3145 + _S3154[int(18)];
    float _S3174 = _S3146 + _S3154[int(19)];
    float _S3175 = _S3147 + _S3154[int(20)];
    float _S3176 = _S3148 + _S3154[int(21)];
    float _S3177 = _S3149 + _S3154[int(22)];
    if(_S3021)
    {
        float _S3178 = 200.0f * _S3040;
        float _S3179 = _S3022 * _S3178 + 0.5f * (_S3020 * _S3178);
        _S3022 = 0.0f;
        _S3025 = _S3179;
    }
    else
    {
        _S3022 = _S3040;
        _S3025 = 0.0f;
    }
    DiffPair_float_0 _S3180;
    (&_S3180)->primal_0 = _S3020;
    (&_S3180)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3180, _S3022);
    float _S3181 = (_S3180.differential_0 + _S3025) / _S3004;
    FixedArray<float, 23>  _S3182;
    _S3182[int(0)] = 0.0f;
    _S3182[int(1)] = 0.0f;
    _S3182[int(2)] = 0.0f;
    _S3182[int(3)] = 0.0f;
    _S3182[int(4)] = 0.0f;
    _S3182[int(5)] = 0.0f;
    _S3182[int(6)] = 0.0f;
    _S3182[int(7)] = 0.0f;
    _S3182[int(8)] = 0.0f;
    _S3182[int(9)] = 0.0f;
    _S3182[int(10)] = 0.0f;
    _S3182[int(11)] = 0.0f;
    _S3182[int(12)] = 0.0f;
    _S3182[int(13)] = 0.0f;
    _S3182[int(14)] = 0.0f;
    _S3182[int(15)] = 0.0f;
    _S3182[int(16)] = 0.0f;
    _S3182[int(17)] = 0.0f;
    _S3182[int(18)] = 0.0f;
    _S3182[int(19)] = 0.0f;
    _S3182[int(20)] = 0.0f;
    _S3182[int(21)] = 0.0f;
    _S3182[int(22)] = 0.0f;
    _S3182[int(13)] = _S3181;
    float _S3183 = _S3155 + _S3182[int(0)];
    float _S3184 = _S3156 + _S3182[int(1)];
    float _S3185 = _S3157 + _S3182[int(2)];
    float _S3186 = _S3158 + _S3182[int(3)];
    float _S3187 = _S3159 + _S3182[int(4)];
    float _S3188 = _S3160 + _S3182[int(5)];
    float _S3189 = _S3161 + _S3182[int(6)];
    float _S3190 = _S3162 + _S3182[int(7)];
    float _S3191 = _S3163 + _S3182[int(8)];
    float _S3192 = _S3164 + _S3182[int(9)];
    float _S3193 = _S3165 + _S3182[int(10)];
    float _S3194 = _S3166 + _S3182[int(11)];
    float _S3195 = _S3167 + _S3182[int(12)];
    float _S3196 = _S3168 + _S3182[int(13)];
    float _S3197 = _S3169 + _S3182[int(14)];
    float _S3198 = _S3170 + _S3182[int(15)];
    float _S3199 = _S3171 + _S3182[int(16)];
    float _S3200 = _S3172 + _S3182[int(17)];
    float _S3201 = _S3173 + _S3182[int(18)];
    float _S3202 = _S3174 + _S3182[int(19)];
    float _S3203 = _S3175 + _S3182[int(20)];
    float _S3204 = _S3176 + _S3182[int(21)];
    float _S3205 = _S3177 + _S3182[int(22)];
    if(_S3018)
    {
        float _S3206 = 200.0f * _S3040;
        float _S3207 = _S3019 * _S3206 + 0.5f * (_S3017 * _S3206);
        _S3019 = 0.0f;
        _S3022 = _S3207;
    }
    else
    {
        _S3019 = _S3040;
        _S3022 = 0.0f;
    }
    DiffPair_float_0 _S3208;
    (&_S3208)->primal_0 = _S3017;
    (&_S3208)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3208, _S3019);
    float _S3209 = (_S3208.differential_0 + _S3022) / _S3004;
    FixedArray<float, 23>  _S3210;
    _S3210[int(0)] = 0.0f;
    _S3210[int(1)] = 0.0f;
    _S3210[int(2)] = 0.0f;
    _S3210[int(3)] = 0.0f;
    _S3210[int(4)] = 0.0f;
    _S3210[int(5)] = 0.0f;
    _S3210[int(6)] = 0.0f;
    _S3210[int(7)] = 0.0f;
    _S3210[int(8)] = 0.0f;
    _S3210[int(9)] = 0.0f;
    _S3210[int(10)] = 0.0f;
    _S3210[int(11)] = 0.0f;
    _S3210[int(12)] = 0.0f;
    _S3210[int(13)] = 0.0f;
    _S3210[int(14)] = 0.0f;
    _S3210[int(15)] = 0.0f;
    _S3210[int(16)] = 0.0f;
    _S3210[int(17)] = 0.0f;
    _S3210[int(18)] = 0.0f;
    _S3210[int(19)] = 0.0f;
    _S3210[int(20)] = 0.0f;
    _S3210[int(21)] = 0.0f;
    _S3210[int(22)] = 0.0f;
    _S3210[int(12)] = _S3209;
    float _S3211 = _S3183 + _S3210[int(0)];
    float _S3212 = _S3184 + _S3210[int(1)];
    float _S3213 = _S3185 + _S3210[int(2)];
    float _S3214 = _S3186 + _S3210[int(3)];
    float _S3215 = _S3187 + _S3210[int(4)];
    float _S3216 = _S3188 + _S3210[int(5)];
    float _S3217 = _S3189 + _S3210[int(6)];
    float _S3218 = _S3190 + _S3210[int(7)];
    float _S3219 = _S3191 + _S3210[int(8)];
    float _S3220 = _S3192 + _S3210[int(9)];
    float _S3221 = _S3193 + _S3210[int(10)];
    float _S3222 = _S3194 + _S3210[int(11)];
    float _S3223 = _S3195 + _S3210[int(12)];
    float _S3224 = _S3196 + _S3210[int(13)];
    float _S3225 = _S3197 + _S3210[int(14)];
    float _S3226 = _S3198 + _S3210[int(15)];
    float _S3227 = _S3199 + _S3210[int(16)];
    float _S3228 = _S3200 + _S3210[int(17)];
    float _S3229 = _S3201 + _S3210[int(18)];
    float _S3230 = _S3202 + _S3210[int(19)];
    float _S3231 = _S3203 + _S3210[int(20)];
    float _S3232 = _S3204 + _S3210[int(21)];
    float _S3233 = _S3205 + _S3210[int(22)];
    if(_S3015)
    {
        float _S3234 = 200.0f * _S3040;
        float _S3235 = _S3016 * _S3234 + 0.5f * (_S3014 * _S3234);
        _S3016 = 0.0f;
        _S3019 = _S3235;
    }
    else
    {
        _S3016 = _S3040;
        _S3019 = 0.0f;
    }
    DiffPair_float_0 _S3236;
    (&_S3236)->primal_0 = _S3014;
    (&_S3236)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3236, _S3016);
    float _S3237 = (_S3236.differential_0 + _S3019) / _S3004;
    FixedArray<float, 23>  _S3238;
    _S3238[int(0)] = 0.0f;
    _S3238[int(1)] = 0.0f;
    _S3238[int(2)] = 0.0f;
    _S3238[int(3)] = 0.0f;
    _S3238[int(4)] = 0.0f;
    _S3238[int(5)] = 0.0f;
    _S3238[int(6)] = 0.0f;
    _S3238[int(7)] = 0.0f;
    _S3238[int(8)] = 0.0f;
    _S3238[int(9)] = 0.0f;
    _S3238[int(10)] = 0.0f;
    _S3238[int(11)] = 0.0f;
    _S3238[int(12)] = 0.0f;
    _S3238[int(13)] = 0.0f;
    _S3238[int(14)] = 0.0f;
    _S3238[int(15)] = 0.0f;
    _S3238[int(16)] = 0.0f;
    _S3238[int(17)] = 0.0f;
    _S3238[int(18)] = 0.0f;
    _S3238[int(19)] = 0.0f;
    _S3238[int(20)] = 0.0f;
    _S3238[int(21)] = 0.0f;
    _S3238[int(22)] = 0.0f;
    _S3238[int(11)] = _S3237;
    float _S3239 = _S3211 + _S3238[int(0)];
    float _S3240 = _S3212 + _S3238[int(1)];
    float _S3241 = _S3213 + _S3238[int(2)];
    float _S3242 = _S3214 + _S3238[int(3)];
    float _S3243 = _S3215 + _S3238[int(4)];
    float _S3244 = _S3216 + _S3238[int(5)];
    float _S3245 = _S3217 + _S3238[int(6)];
    float _S3246 = _S3218 + _S3238[int(7)];
    float _S3247 = _S3219 + _S3238[int(8)];
    float _S3248 = _S3220 + _S3238[int(9)];
    float _S3249 = _S3221 + _S3238[int(10)];
    float _S3250 = _S3222 + _S3238[int(11)];
    float _S3251 = _S3223 + _S3238[int(12)];
    float _S3252 = _S3224 + _S3238[int(13)];
    float _S3253 = _S3225 + _S3238[int(14)];
    float _S3254 = _S3226 + _S3238[int(15)];
    float _S3255 = _S3227 + _S3238[int(16)];
    float _S3256 = _S3228 + _S3238[int(17)];
    float _S3257 = _S3229 + _S3238[int(18)];
    float _S3258 = _S3230 + _S3238[int(19)];
    float _S3259 = _S3231 + _S3238[int(20)];
    float _S3260 = _S3232 + _S3238[int(21)];
    float _S3261 = _S3233 + _S3238[int(22)];
    if(_S3012)
    {
        float _S3262 = 200.0f * _S3040;
        float _S3263 = _S3013 * _S3262 + 0.5f * (_S3011 * _S3262);
        _S3013 = 0.0f;
        _S3016 = _S3263;
    }
    else
    {
        _S3013 = _S3040;
        _S3016 = 0.0f;
    }
    DiffPair_float_0 _S3264;
    (&_S3264)->primal_0 = _S3011;
    (&_S3264)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3264, _S3013);
    float _S3265 = (_S3264.differential_0 + _S3016) / _S3004;
    float _S3266 = _S3035 / _S3010;
    float _S3267 = _S3036 / _S3009;
    float _S3268 = _S3037 / _S3008;
    FixedArray<float, 23>  _S3269;
    _S3269[int(0)] = 0.0f;
    _S3269[int(1)] = 0.0f;
    _S3269[int(2)] = 0.0f;
    _S3269[int(3)] = 0.0f;
    _S3269[int(4)] = 0.0f;
    _S3269[int(5)] = 0.0f;
    _S3269[int(6)] = 0.0f;
    _S3269[int(7)] = 0.0f;
    _S3269[int(8)] = 0.0f;
    _S3269[int(9)] = 0.0f;
    _S3269[int(10)] = 0.0f;
    _S3269[int(11)] = 0.0f;
    _S3269[int(12)] = 0.0f;
    _S3269[int(13)] = 0.0f;
    _S3269[int(14)] = 0.0f;
    _S3269[int(15)] = 0.0f;
    _S3269[int(16)] = 0.0f;
    _S3269[int(17)] = 0.0f;
    _S3269[int(18)] = 0.0f;
    _S3269[int(19)] = 0.0f;
    _S3269[int(20)] = 0.0f;
    _S3269[int(21)] = 0.0f;
    _S3269[int(22)] = 0.0f;
    _S3269[int(10)] = _S3265;
    _S3269[int(9)] = _S3266;
    _S3269[int(8)] = _S3266;
    _S3269[int(7)] = _S3266;
    _S3269[int(6)] = _S3266;
    _S3269[int(5)] = _S3266;
    _S3269[int(4)] = _S3267;
    _S3269[int(3)] = _S3267;
    _S3269[int(2)] = _S3267;
    _S3269[int(1)] = _S3268;
    float _S3270 = _S3239 + _S3269[int(0)];
    float _S3271 = _S3240 + _S3269[int(1)];
    float _S3272 = _S3241 + _S3269[int(2)];
    float _S3273 = _S3242 + _S3269[int(3)];
    float _S3274 = _S3243 + _S3269[int(4)];
    float _S3275 = _S3244 + _S3269[int(5)];
    float _S3276 = _S3245 + _S3269[int(6)];
    float _S3277 = _S3246 + _S3269[int(7)];
    float _S3278 = _S3247 + _S3269[int(8)];
    float _S3279 = _S3248 + _S3269[int(9)];
    float _S3280 = _S3249 + _S3269[int(10)];
    float _S3281 = _S3250 + _S3269[int(11)];
    float _S3282 = _S3251 + _S3269[int(12)];
    float _S3283 = _S3252 + _S3269[int(13)];
    float _S3284 = _S3253 + _S3269[int(14)];
    float _S3285 = _S3254 + _S3269[int(15)];
    float _S3286 = _S3255 + _S3269[int(16)];
    float _S3287 = _S3256 + _S3269[int(17)];
    float _S3288 = _S3257 + _S3269[int(18)];
    float _S3289 = _S3258 + _S3269[int(19)];
    float _S3290 = _S3259 + _S3269[int(20)];
    float _S3291 = _S3260 + _S3269[int(21)];
    float _S3292 = _S3261 + _S3269[int(22)];
    if(_S3006)
    {
        float _S3293 = 10.0f * _S3038;
        float _S3294 = _S3007 * _S3293 + 0.5f * (_S3005 * _S3293);
        _S3007 = 0.0f;
        _S3013 = _S3294;
    }
    else
    {
        _S3007 = _S3038;
        _S3013 = 0.0f;
    }
    DiffPair_float_0 _S3295;
    (&_S3295)->primal_0 = _S3005;
    (&_S3295)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3295, _S3007);
    float _S3296 = (_S3295.differential_0 + _S3013) / _S3004;
    FixedArray<float, 23>  _S3297;
    _S3297[int(0)] = 0.0f;
    _S3297[int(1)] = 0.0f;
    _S3297[int(2)] = 0.0f;
    _S3297[int(3)] = 0.0f;
    _S3297[int(4)] = 0.0f;
    _S3297[int(5)] = 0.0f;
    _S3297[int(6)] = 0.0f;
    _S3297[int(7)] = 0.0f;
    _S3297[int(8)] = 0.0f;
    _S3297[int(9)] = 0.0f;
    _S3297[int(10)] = 0.0f;
    _S3297[int(11)] = 0.0f;
    _S3297[int(12)] = 0.0f;
    _S3297[int(13)] = 0.0f;
    _S3297[int(14)] = 0.0f;
    _S3297[int(15)] = 0.0f;
    _S3297[int(16)] = 0.0f;
    _S3297[int(17)] = 0.0f;
    _S3297[int(18)] = 0.0f;
    _S3297[int(19)] = 0.0f;
    _S3297[int(20)] = 0.0f;
    _S3297[int(21)] = 0.0f;
    _S3297[int(22)] = 0.0f;
    _S3297[int(0)] = _S3296;
    FixedArray<float, 23>  _S3298 = {
        _S3270 + _S3297[int(0)], _S3271 + _S3297[int(1)], _S3272 + _S3297[int(2)], _S3273 + _S3297[int(3)], _S3274 + _S3297[int(4)], _S3275 + _S3297[int(5)], _S3276 + _S3297[int(6)], _S3277 + _S3297[int(7)], _S3278 + _S3297[int(8)], _S3279 + _S3297[int(9)], _S3280 + _S3297[int(10)], _S3281 + _S3297[int(11)], _S3282 + _S3297[int(12)], _S3283 + _S3297[int(13)], _S3284 + _S3297[int(14)], _S3285 + _S3297[int(15)], _S3286 + _S3297[int(16)], _S3287 + _S3297[int(17)], _S3288 + _S3297[int(18)], _S3289 + _S3297[int(19)], _S3290 + _S3297[int(20)], _S3291 + _S3297[int(21)], _S3292 + _S3297[int(22)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S3298;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S3299, int _S3300, FixedArray<float, 6>  * _S3301, FixedArray<float, 6>  * _S3302)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S3299, _S3300, _S3301, _S3302);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_5, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_5, FixedArray<float, 23>  * _S3303)
{
    FixedArray<float, 23>  _S3304 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S3304;
    FixedArray<float, 6>  _S3305 = loss_weights_5;
    FixedArray<float, 6>  _S3306 = grad_out_5;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_2, num_cameras_5, &_S3305, &_S3306);
    *_S3303 = (&dp_raw_losses_2)->differential_0;
    return;
}

