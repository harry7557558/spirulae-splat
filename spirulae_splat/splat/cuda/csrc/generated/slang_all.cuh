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

inline __device__ float2  s_primal_ctx_exp_1(float2  _S29)
{
    return exp_1(_S29);
}

inline __device__ float s_primal_ctx_max_0(float _S30, float _S31)
{
    return (F32_max((_S30), (_S31)));
}

inline __device__ float s_primal_ctx_min_0(float _S32, float _S33)
{
    return (F32_min((_S32), (_S33)));
}

inline __device__ float s_primal_ctx_log_0(float _S34)
{
    return (F32_log((_S34)));
}

inline __device__ float3  s_primal_ctx_exp_2(float3  _S35)
{
    return exp_0(_S35);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S36, DiffPair_float_0 * _S37, float _S38)
{
    _d_max_0(_S36, _S37, _S38);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S39, float _S40)
{
    _d_log_0(_S39, _S40);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S41, float _S42)
{
    _d_exp_0(_S41, _S42);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S43, DiffPair_float_0 * _S44, float _S45)
{
    _d_min_0(_S43, _S44, _S45);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S46, float2  _S47)
{
    _d_exp_vector_1(_S46, _S47);
    return;
}

inline __device__ void s_bwd_prop_exp_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S48, float3  _S49)
{
    _d_exp_vector_0(_S48, _S49);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S50, float _S51)
{
    _d_sqrt_0(_S50, _S51);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_8, float _s_dOut_0)
{
    float _S52 = (*dpx_8).primal_0.x;
    float _S53 = (*dpx_8).primal_0.y;
    float _S54 = (*dpx_8).primal_0.z;
    float _S55 = (*dpx_8).primal_0.w;
    DiffPair_float_0 _S56;
    (&_S56)->primal_0 = _S52 * _S52 + _S53 * _S53 + _S54 * _S54 + _S55 * _S55;
    (&_S56)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S56, _s_dOut_0);
    float _S57 = (*dpx_8).primal_0.w * _S56.differential_0;
    float _S58 = _S57 + _S57;
    float _S59 = (*dpx_8).primal_0.z * _S56.differential_0;
    float _S60 = _S59 + _S59;
    float _S61 = (*dpx_8).primal_0.y * _S56.differential_0;
    float _S62 = _S61 + _S61;
    float _S63 = (*dpx_8).primal_0.x * _S56.differential_0;
    float _S64 = _S63 + _S63;
    float4  _S65 = make_float4 (0.0f);
    *&((&_S65)->w) = _S58;
    *&((&_S65)->z) = _S60;
    *&((&_S65)->y) = _S62;
    *&((&_S65)->x) = _S64;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S65;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S66, float _S67)
{
    s_bwd_prop_length_impl_0(_S66, _S67);
    return;
}

inline __device__ void s_bwd_prop_per_splat_losses_0(bool is_3dgs_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1, FixedArray<float, 5>  * _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S68 = *dpscales_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S69 = *dpquat_0;
    float2  _S70 = make_float2 (0.0f);
    float3  _S71 = make_float3 (0.0f);
    float _S72 = - (*dpopacity_0).primal_0;
    float _S73 = 1.0f + s_primal_ctx_exp_0(_S72);
    float _S74 = _S73 * _S73;
    float _S75 = length_0((*dpquat_0).primal_0);
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
    float _S114;
    float _S115;
    float _S116;
    float _S117;
    float _S118;
    float _S119;
    float _S120;
    float _S121;
    float _S122;
    float _S123;
    float2  _S124;
    float2  _S125;
    float3  _S126;
    if(is_3dgs_1)
    {
        float3  _S127 = s_primal_ctx_exp_2(_S68.primal_0);
        float _S128 = _S127.x;
        float _S129 = _S127.y;
        float _S130 = _S127.z;
        float _S131 = s_primal_ctx_max_0(_S128, _S129);
        float _S132 = s_primal_ctx_max_0(_S131, _S130);
        float _S133 = s_primal_ctx_min_0(_S128, _S129);
        float _S134 = s_primal_ctx_min_0(_S133, _S130);
        float _S135 = _S132 / _S134;
        float _S136 = _S134 * _S134;
        float3  _S137 = make_float3 (2.0f) * _S68.primal_0;
        float3  _S138 = s_primal_ctx_exp_2(_S137);
        float x_10 = _S138.x;
        float y_5 = _S138.y;
        float z_1 = _S138.z;
        float s_2 = x_10 + y_5 + z_1;
        float _S139 = s_primal_ctx_max_0(x_10, y_5);
        float _S140 = s_primal_ctx_max_0(_S139, z_1);
        float s1_2 = _S140 / s_2;
        float _S141 = s_2 * s_2;
        float _S142 = s_primal_ctx_min_0(x_10, y_5);
        float _S143 = s_primal_ctx_min_0(_S142, z_1);
        float s3_1 = _S143 / s_2;
        float s2_2 = 1.0f - s1_2 - s3_1;
        float _S144 = - s1_2;
        float _S145 = s_primal_ctx_log_0(s1_2);
        float _S146 = s_primal_ctx_log_0(s2_2);
        float _S147 = s_primal_ctx_log_0(s3_1);
        float _S148 = _S144 * _S145 - s2_2 * _S146 - s3_1 * _S147;
        float _S149 = s_primal_ctx_exp_0(_S148) - 0.99998998641967773f;
        float _S150 = - s_primal_ctx_log_0(_S149);
        _S76 = 0.0f;
        _S77 = 0.0f;
        _S78 = 0.0f;
        _S79 = 0.0f;
        _S80 = 0.0f;
        _S81 = 0.0f;
        _S82 = 0.0f;
        _S83 = 0.0f;
        _S84 = 0.0f;
        _S85 = 0.0f;
        _S86 = 0.0f;
        _S87 = 0.0f;
        _S88 = 0.0f;
        _S89 = 0.0f;
        _S124 = _S70;
        _S90 = 0.0f;
        _S91 = 0.0f;
        _S92 = 0.0f;
        _S93 = 0.0f;
        _S94 = 0.0f;
        _S95 = 0.0f;
        _S125 = _S70;
        _S96 = _S150;
        _S97 = _S149;
        _S98 = _S148;
        _S99 = s3_1;
        _S100 = _S147;
        _S101 = s2_2;
        _S102 = _S146;
        _S103 = _S144;
        _S104 = _S145;
        _S105 = s1_2;
        _S106 = _S141;
        _S107 = _S143;
        _S108 = s_2;
        _S109 = _S142;
        _S110 = z_1;
        _S111 = x_10;
        _S112 = y_5;
        _S113 = _S140;
        _S114 = _S139;
        _S126 = _S137;
        _S115 = _S135;
        _S116 = _S136;
        _S117 = _S132;
        _S118 = _S134;
        _S119 = _S133;
        _S120 = _S130;
        _S121 = _S128;
        _S122 = _S129;
        _S123 = _S131;
    }
    else
    {
        float2  _S151 = float2 {_S68.primal_0.x, _S68.primal_0.y};
        float2  _S152 = s_primal_ctx_exp_1(_S151);
        float _S153 = _S152.x;
        float _S154 = _S152.y;
        float _S155 = s_primal_ctx_max_0(_S153, _S154);
        float _S156 = s_primal_ctx_min_0(_S153, _S154);
        float _S157 = _S155 / _S156;
        float _S158 = _S156 * _S156;
        float2  _S159 = make_float2 (2.0f) * _S151;
        float2  _S160 = s_primal_ctx_exp_1(_S159);
        float x_11 = _S160.x;
        float y_6 = _S160.y;
        float s_3 = x_11 + y_6;
        float _S161 = s_primal_ctx_max_0(x_11, y_6);
        float s1_3 = _S161 / s_3;
        float _S162 = s_3 * s_3;
        float _S163 = s_primal_ctx_min_0(x_11, y_6);
        float s2_3 = _S163 / s_3;
        float _S164 = - s1_3;
        float _S165 = s_primal_ctx_log_0(s1_3);
        float _S166 = s_primal_ctx_log_0(s2_3);
        float _S167 = _S164 * _S165 - s2_3 * _S166;
        float _S168 = s_primal_ctx_exp_0(_S167) - 0.99998998641967773f;
        _S76 = - s_primal_ctx_log_0(_S168);
        _S77 = _S168;
        _S78 = _S167;
        _S79 = s2_3;
        _S80 = _S166;
        _S81 = _S164;
        _S82 = _S165;
        _S83 = s1_3;
        _S84 = _S162;
        _S85 = _S163;
        _S86 = s_3;
        _S87 = x_11;
        _S88 = y_6;
        _S89 = _S161;
        _S124 = _S159;
        _S90 = _S157;
        _S91 = _S158;
        _S92 = _S155;
        _S93 = _S156;
        _S94 = _S153;
        _S95 = _S154;
        _S125 = _S151;
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
        _S108 = 0.0f;
        _S109 = 0.0f;
        _S110 = 0.0f;
        _S111 = 0.0f;
        _S112 = 0.0f;
        _S113 = 0.0f;
        _S114 = 0.0f;
        _S126 = _S71;
        _S115 = 0.0f;
        _S116 = 0.0f;
        _S117 = 0.0f;
        _S118 = 0.0f;
        _S119 = 0.0f;
        _S120 = 0.0f;
        _S121 = 0.0f;
        _S122 = 0.0f;
        _S123 = 0.0f;
    }
    if(is_3dgs_1)
    {
        float _S169 = erank_reg_weight_s3_1 * (*_s_dOut_1)[int(3)];
        float _S170 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S171;
        (&_S171)->primal_0 = _S96;
        (&_S171)->differential_0 = 0.0f;
        DiffPair_float_0 _S172;
        (&_S172)->primal_0 = 0.0f;
        (&_S172)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S171, &_S172, _S170);
        float _S173 = - _S171.differential_0;
        DiffPair_float_0 _S174;
        (&_S174)->primal_0 = _S97;
        (&_S174)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S174, _S173);
        DiffPair_float_0 _S175;
        (&_S175)->primal_0 = _S98;
        (&_S175)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S175, _S174.differential_0);
        float _S176 = - _S175.differential_0;
        float _S177 = _S99 * _S176;
        float _S178 = _S100 * _S176;
        DiffPair_float_0 _S179;
        (&_S179)->primal_0 = _S99;
        (&_S179)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S179, _S177);
        float _S180 = _S101 * _S176;
        float _S181 = _S102 * _S176;
        DiffPair_float_0 _S182;
        (&_S182)->primal_0 = _S101;
        (&_S182)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S182, _S180);
        float _S183 = _S103 * _S175.differential_0;
        float _S184 = _S104 * _S175.differential_0;
        DiffPair_float_0 _S185;
        (&_S185)->primal_0 = _S105;
        (&_S185)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S185, _S183);
        float _S186 = - _S184;
        float _S187 = - (_S181 + _S182.differential_0);
        float _S188 = (_S169 + _S178 + _S179.differential_0 + _S187) / _S106;
        float _S189 = _S107 * - _S188;
        float _S190 = _S108 * _S188;
        DiffPair_float_0 _S191;
        (&_S191)->primal_0 = _S109;
        (&_S191)->differential_0 = 0.0f;
        DiffPair_float_0 _S192;
        (&_S192)->primal_0 = _S110;
        (&_S192)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S191, &_S192, _S190);
        DiffPair_float_0 _S193;
        (&_S193)->primal_0 = _S111;
        (&_S193)->differential_0 = 0.0f;
        DiffPair_float_0 _S194;
        (&_S194)->primal_0 = _S112;
        (&_S194)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S193, &_S194, _S191.differential_0);
        float _S195 = (_S185.differential_0 + _S186 + _S187) / _S106;
        float _S196 = _S113 * - _S195;
        float _S197 = _S108 * _S195;
        DiffPair_float_0 _S198;
        (&_S198)->primal_0 = _S114;
        (&_S198)->differential_0 = 0.0f;
        DiffPair_float_0 _S199;
        (&_S199)->primal_0 = _S110;
        (&_S199)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S198, &_S199, _S197);
        DiffPair_float_0 _S200;
        (&_S200)->primal_0 = _S111;
        (&_S200)->differential_0 = 0.0f;
        DiffPair_float_0 _S201;
        (&_S201)->primal_0 = _S112;
        (&_S201)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S200, &_S201, _S198.differential_0);
        float _S202 = _S189 + _S196;
        float3  _S203 = make_float3 (_S193.differential_0 + _S200.differential_0 + _S202, _S194.differential_0 + _S201.differential_0 + _S202, _S192.differential_0 + _S199.differential_0 + _S202);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S204;
        (&_S204)->primal_0 = _S126;
        (&_S204)->differential_0 = _S71;
        s_bwd_prop_exp_2(&_S204, _S203);
        float3  _S205 = make_float3 (2.0f) * _S204.differential_0;
        float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S206;
        (&_S206)->primal_0 = _S115;
        (&_S206)->differential_0 = 0.0f;
        DiffPair_float_0 _S207;
        (&_S207)->primal_0 = max_gauss_ratio_1;
        (&_S207)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S206, &_S207, s_diff_scale_reg_T_0);
        float _S208 = _S206.differential_0 / _S116;
        float _S209 = _S117 * - _S208;
        float _S210 = _S118 * _S208;
        DiffPair_float_0 _S211;
        (&_S211)->primal_0 = _S119;
        (&_S211)->differential_0 = 0.0f;
        DiffPair_float_0 _S212;
        (&_S212)->primal_0 = _S120;
        (&_S212)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S211, &_S212, _S209);
        DiffPair_float_0 _S213;
        (&_S213)->primal_0 = _S121;
        (&_S213)->differential_0 = 0.0f;
        DiffPair_float_0 _S214;
        (&_S214)->primal_0 = _S122;
        (&_S214)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S213, &_S214, _S211.differential_0);
        DiffPair_float_0 _S215;
        (&_S215)->primal_0 = _S123;
        (&_S215)->differential_0 = 0.0f;
        DiffPair_float_0 _S216;
        (&_S216)->primal_0 = _S120;
        (&_S216)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S215, &_S216, _S210);
        DiffPair_float_0 _S217;
        (&_S217)->primal_0 = _S121;
        (&_S217)->differential_0 = 0.0f;
        DiffPair_float_0 _S218;
        (&_S218)->primal_0 = _S122;
        (&_S218)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S217, &_S218, _S215.differential_0);
        float _S219 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * (*_s_dOut_1)[int(1)]);
        float3  _S220 = make_float3 (_S213.differential_0 + _S217.differential_0 + _S219, _S214.differential_0 + _S218.differential_0 + _S219, _S212.differential_0 + _S216.differential_0 + _S219);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S221;
        (&_S221)->primal_0 = _S68.primal_0;
        (&_S221)->differential_0 = _S71;
        s_bwd_prop_exp_2(&_S221, _S220);
        float3  _S222 = _S205 + _S221.differential_0;
        _S76 = (*_s_dOut_1)[int(4)];
        _S77 = (*_s_dOut_1)[int(0)];
        _S126 = _S222;
    }
    else
    {
        float _S223 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S224;
        (&_S224)->primal_0 = _S76;
        (&_S224)->differential_0 = 0.0f;
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = 0.0f;
        (&_S225)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S224, &_S225, _S223);
        float _S226 = - _S224.differential_0;
        DiffPair_float_0 _S227;
        (&_S227)->primal_0 = _S77;
        (&_S227)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S227, _S226);
        DiffPair_float_0 _S228;
        (&_S228)->primal_0 = _S78;
        (&_S228)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S228, _S227.differential_0);
        float _S229 = - _S228.differential_0;
        float _S230 = _S79 * _S229;
        float _S231 = _S80 * _S229;
        DiffPair_float_0 _S232;
        (&_S232)->primal_0 = _S79;
        (&_S232)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S232, _S230);
        float _S233 = _S81 * _S228.differential_0;
        float _S234 = _S82 * _S228.differential_0;
        DiffPair_float_0 _S235;
        (&_S235)->primal_0 = _S83;
        (&_S235)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S235, _S233);
        float _S236 = - _S234;
        float _S237 = (_S231 + _S232.differential_0) / _S84;
        float _S238 = _S85 * - _S237;
        float _S239 = _S86 * _S237;
        DiffPair_float_0 _S240;
        (&_S240)->primal_0 = _S87;
        (&_S240)->differential_0 = 0.0f;
        DiffPair_float_0 _S241;
        (&_S241)->primal_0 = _S88;
        (&_S241)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S240, &_S241, _S239);
        float _S242 = (_S235.differential_0 + _S236) / _S84;
        float _S243 = _S89 * - _S242;
        float _S244 = _S86 * _S242;
        DiffPair_float_0 _S245;
        (&_S245)->primal_0 = _S87;
        (&_S245)->differential_0 = 0.0f;
        DiffPair_float_0 _S246;
        (&_S246)->primal_0 = _S88;
        (&_S246)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S245, &_S246, _S244);
        float _S247 = _S238 + _S243;
        float2  _S248 = make_float2 (_S240.differential_0 + _S245.differential_0 + _S247, _S241.differential_0 + _S246.differential_0 + _S247);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S249;
        (&_S249)->primal_0 = _S124;
        (&_S249)->differential_0 = _S70;
        s_bwd_prop_exp_1(&_S249, _S248);
        float2  _S250 = make_float2 (2.0f) * _S249.differential_0;
        float s_diff_scale_reg_T_1 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S251;
        (&_S251)->primal_0 = _S90;
        (&_S251)->differential_0 = 0.0f;
        DiffPair_float_0 _S252;
        (&_S252)->primal_0 = max_gauss_ratio_1;
        (&_S252)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S251, &_S252, s_diff_scale_reg_T_1);
        float _S253 = _S251.differential_0 / _S91;
        float _S254 = _S92 * - _S253;
        float _S255 = _S93 * _S253;
        DiffPair_float_0 _S256;
        (&_S256)->primal_0 = _S94;
        (&_S256)->differential_0 = 0.0f;
        DiffPair_float_0 _S257;
        (&_S257)->primal_0 = _S95;
        (&_S257)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S256, &_S257, _S254);
        DiffPair_float_0 _S258;
        (&_S258)->primal_0 = _S94;
        (&_S258)->differential_0 = 0.0f;
        DiffPair_float_0 _S259;
        (&_S259)->primal_0 = _S95;
        (&_S259)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S258, &_S259, _S255);
        float _S260 = mcmc_scale_reg_weight_1 * (0.5f * (*_s_dOut_1)[int(1)]);
        float2  _S261 = make_float2 (_S256.differential_0 + _S258.differential_0 + _S260, _S257.differential_0 + _S259.differential_0 + _S260);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S262;
        (&_S262)->primal_0 = _S125;
        (&_S262)->differential_0 = _S70;
        s_bwd_prop_exp_1(&_S262, _S261);
        float2  _S263 = _S250 + _S262.differential_0;
        float3  _S264 = make_float3 (_S263.x, _S263.y, 0.0f);
        _S76 = (*_s_dOut_1)[int(4)];
        _S77 = (*_s_dOut_1)[int(0)];
        _S126 = _S264;
    }
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * _S76;
    float _S265 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S266;
    (&_S266)->primal_0 = _S75;
    (&_S266)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S266, _S265);
    float _S267 = _S266.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S268 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S269;
    (&_S269)->primal_0 = _S69.primal_0;
    (&_S269)->differential_0 = _S268;
    s_bwd_length_impl_0(&_S269, _S267);
    float _S270 = - (mcmc_opacity_reg_weight_1 * _S77 / _S74);
    DiffPair_float_0 _S271;
    (&_S271)->primal_0 = _S72;
    (&_S271)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S271, _S270);
    float _S272 = - _S271.differential_0;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S269.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S272;
    dpscales_0->primal_0 = (*dpscales_0).primal_0;
    dpscales_0->differential_0 = _S126;
    return;
}

inline __device__ void s_bwd_per_splat_losses_0(bool _S273, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S274, DiffPair_float_0 * _S275, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S276, float _S277, float _S278, float _S279, float _S280, float _S281, float _S282, float _S283, FixedArray<float, 5>  * _S284)
{
    s_bwd_prop_per_splat_losses_0(_S273, _S274, _S275, _S276, _S277, _S278, _S279, _S280, _S281, _S282, _S283, _S284);
    return;
}

inline __device__ void per_splat_losses_bwd(bool is_3dgs_2, float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  * v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float3  _S285 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_scales_0;
    (&p_scales_0)->primal_0 = scales_1;
    (&p_scales_0)->differential_0 = _S285;
    DiffPair_float_0 p_opacity_0;
    (&p_opacity_0)->primal_0 = opacity_1;
    (&p_opacity_0)->differential_0 = 0.0f;
    float4  _S286 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 p_quat_0;
    (&p_quat_0)->primal_0 = quat_1;
    (&p_quat_0)->differential_0 = _S286;
    s_bwd_per_splat_losses_0(is_3dgs_2, &p_scales_0, &p_opacity_0, &p_quat_0, mcmc_opacity_reg_weight_2, mcmc_scale_reg_weight_2, max_gauss_ratio_2, scale_regularization_weight_2, erank_reg_weight_2, erank_reg_weight_s3_2, quat_norm_reg_weight_2, v_loss_0);
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
    float _S287 = (*left_0).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_8.x;
    float sum_0 = _S287 + (*left_0).primal_0.rows[int(1)].x * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_8.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_8.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S288 = (*left_0).primal_0.rows[int(0)].y * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_8.x;
    float sum_2 = _S288 + (*left_0).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_8.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_8.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S289 = (*left_0).primal_0.rows[int(0)].z * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_8.x;
    float sum_4 = _S289 + (*left_0).primal_0.rows[int(1)].z * dOut_8.y;
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
    float _S290 = (*left_1).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 2, 2>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_1).primal_0.x * dOut_9.x;
    float sum_6 = _S290 + (*left_1).primal_0.rows[int(1)].x * dOut_9.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_1).primal_0.x * dOut_9.y;
    float2  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_6;
    float _S291 = (*left_1).primal_0.rows[int(0)].y * dOut_9.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_1).primal_0.y * dOut_9.x;
    float sum_7 = _S291 + (*left_1).primal_0.rows[int(1)].y * dOut_9.y;
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
    float _S292 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S293 = (*dist_coeffs_0)[int(1)] + r2_0 * _S292;
    float _S294 = (*dist_coeffs_0)[int(0)] + r2_0 * _S293;
    float2  _S295 = make_float2 (1.0f + r2_0 * _S294);
    float _S296 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S297 = _S296 * u_0;
    float _S298 = 2.0f * u_0;
    float _S299 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S300 = _S299 * u_0;
    float _S301 = 2.0f * v_0;
    float2  _S302 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S303 = uv_0 * _S302;
    float _S304 = (*dist_coeffs_0)[int(4)] * _S302.y;
    float _S305 = (*dist_coeffs_0)[int(5)] * _S302.x;
    float _S306 = _S303.x + _S303.y;
    float _S307 = r2_0 * _S306;
    float _S308 = r2_0 * _S307;
    float _S309 = (*dist_coeffs_0)[int(7)] * _S302.y + _S304 + (*dist_coeffs_0)[int(6)] * _S302.x + _S305 + _S294 * _S306 + _S293 * _S307 + _S292 * _S308 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S308);
    float _S310 = v_0 * _S309;
    float _S311 = u_0 * _S309;
    float2  _S312 = make_float2 (0.0f, 1.0f) + make_float2 ((*dist_coeffs_0)[int(8)] * 0.0f, (*dist_coeffs_0)[int(9)] * 0.0f);
    float2  _S313 = uv_0 * _S312;
    float _S314 = (*dist_coeffs_0)[int(4)] * _S312.y;
    float _S315 = (*dist_coeffs_0)[int(5)] * _S312.x;
    float _S316 = _S313.x + _S313.y;
    float _S317 = r2_0 * _S316;
    float _S318 = r2_0 * _S317;
    float _S319 = (*dist_coeffs_0)[int(7)] * _S312.y + _S314 + (*dist_coeffs_0)[int(6)] * _S312.x + _S315 + _S294 * _S316 + _S293 * _S317 + _S292 * _S318 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S318);
    float _S320 = v_0 * _S319;
    float _S321 = u_0 * _S319;
    return makeMatrix<float, 2, 2> (_S295 * _S302 + make_float2 (_S299 * (v_0 * _S302.y) + _S298 * _S305 + 2.0f * (u_0 * _S305) + _S296 * (v_0 * _S302.x) + _S311 + _S311, _S301 * _S304 + 2.0f * (v_0 * _S304) + _S300 * _S302.y + _S297 * _S302.x + _S310 + _S310), _S295 * _S312 + make_float2 (_S299 * (v_0 * _S312.y) + _S298 * _S315 + 2.0f * (u_0 * _S315) + _S296 * (v_0 * _S312.x) + _S321 + _S321, _S301 * _S314 + 2.0f * (v_0 * _S314) + _S300 * _S312.y + _S297 * _S312.x + _S320 + _S320));
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  * dist_coeffs_1)
{
    Matrix<float, 2, 2>  _S322 = camera_distortion_jac_0(uv_1, dist_coeffs_1);
    return (F32_min((determinant_0(_S322)), ((F32_min((_S322.rows[int(0)].x), (_S322.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_2)
{
    float2  _S323;
    if(is_fisheye_0)
    {
        float r_3 = length_1(uv_2);
        float theta_0 = (F32_atan((r_3)));
        float _S324;
        if(r_3 < 0.00100000004749745f)
        {
            _S324 = 1.0f - r_3 * r_3 / 3.0f;
        }
        else
        {
            _S324 = theta_0 / r_3;
        }
        _S323 = uv_2 * make_float2 (_S324);
    }
    else
    {
        _S323 = uv_2;
    }
    float u_1 = _S323.x;
    float v_1 = _S323.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S325 = _S323 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_2)[int(0)] + r2_1 * ((*dist_coeffs_2)[int(1)] + r2_1 * ((*dist_coeffs_2)[int(2)] + r2_1 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_1 * v_1 + (*dist_coeffs_2)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_2)[int(6)] * r2_1, 2.0f * (*dist_coeffs_2)[int(5)] * u_1 * v_1 + (*dist_coeffs_2)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_2)[int(7)] * r2_1);
    return _S325 + make_float2 ((*dist_coeffs_2)[int(8)] * _S325.x + (*dist_coeffs_2)[int(9)] * _S325.y, 0.0f);
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
        float2  _S326 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_4 = _S326 + make_float2 ((*dist_coeffs_3)[int(8)] * _S326.x + (*dist_coeffs_3)[int(9)] * _S326.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S327 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S327.rows[int(0)].x * _S327.rows[int(1)].y - _S327.rows[int(0)].y * _S327.rows[int(1)].x);
        float _S328 = r_4.x;
        float _S329 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S328 * _S327.rows[int(1)].y - _S329 * _S327.rows[int(0)].y) * inv_det_0, (- _S328 * _S327.rows[int(1)].x + _S329 * _S327.rows[int(0)].x) * inv_det_0);
        i_8 = i_8 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S330 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S331;
    if((F32_min((determinant_0(_S330)), ((F32_min((_S330.rows[int(0)].x), (_S330.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S332 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S331 = (length_1(_S332 + make_float2 ((*dist_coeffs_3)[int(8)] * _S332.x + (*dist_coeffs_3)[int(9)] * _S332.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S331 = false;
    }
    return _S331;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S333 = uv_4;
    bool _S334 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S333);
    if(!_S334)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S335 = _S333;
        float theta_1 = length_1(_S333);
        float _S336;
        if(theta_1 < 0.00100000004749745f)
        {
            _S336 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S336 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S337 = make_float3 ((_S335 * make_float2 (_S336)).x, (_S335 * make_float2 (_S336)).y, (F32_cos((theta_1))));
        raydir_0 = _S337;
    }
    else
    {
        raydir_0 = make_float3 (_S333.x, _S333.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_5, float3  * raydir_1)
{
    float2  _S338 = uv_5;
    int3  _S339 = make_int3 (int(0));
    float3  _S340 = make_float3 ((float)_S339.x, (float)_S339.y, (float)_S339.z);
    *raydir_1 = _S340;
    bool _S341 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S338);
    if(!_S341)
    {
        return false;
    }
    float3  _S342;
    if(is_fisheye_2)
    {
        float2  _S343 = _S338;
        float theta_2 = length_1(_S338);
        float _S344;
        if(theta_2 < 0.00100000004749745f)
        {
            _S344 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S344 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S345 = make_float3 ((_S343 * make_float2 (_S344)).x, (_S343 * make_float2 (_S344)).y, (F32_cos((theta_2))));
        _S342 = _S345;
    }
    else
    {
        _S342 = make_float3 (_S338.x, _S338.y, 1.0f);
    }
    *raydir_1 = _S342;
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

inline __device__ bool generate_ray(float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_6, float3  * raydir_2)
{
    float2  _S346 = uv_6;
    bool _S347 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S346);
    if(!_S347)
    {
        int3  _S348 = make_int3 (int(0));
        float3  _S349 = make_float3 ((float)_S348.x, (float)_S348.y, (float)_S348.z);
        *raydir_2 = _S349;
        return false;
    }
    float3  _S350;
    if(is_fisheye_3)
    {
        float2  _S351 = _S346;
        float theta_3 = length_1(_S346);
        float _S352;
        if(theta_3 < 0.00100000004749745f)
        {
            _S352 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S352 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S353 = make_float3 ((_S351 * make_float2 (_S352)).x, (_S351 * make_float2 (_S352)).y, (F32_cos((theta_3))));
        _S350 = _S353;
    }
    else
    {
        _S350 = make_float3 (_S346.x, _S346.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S350);
    return true;
}

inline __device__ void _d_mul_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_6, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_6, float3  dOut_11)
{
    float _S354 = (*right_6).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 3, 3>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = (*left_6).primal_0.x * dOut_11.x;
    float sum_14 = _S354 + (*right_6).primal_0.rows[int(0)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = (*left_6).primal_0.x * dOut_11.y;
    float sum_15 = sum_14 + (*right_6).primal_0.rows[int(0)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(0)))->z) = (*left_6).primal_0.x * dOut_11.z;
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = sum_15;
    float _S355 = (*right_6).primal_0.rows[int(1)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = (*left_6).primal_0.y * dOut_11.x;
    float sum_16 = _S355 + (*right_6).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = (*left_6).primal_0.y * dOut_11.y;
    float sum_17 = sum_16 + (*right_6).primal_0.rows[int(1)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(1)))->z) = (*left_6).primal_0.y * dOut_11.z;
    *&((&left_d_result_3)->y) = sum_17;
    float _S356 = (*right_6).primal_0.rows[int(2)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = (*left_6).primal_0.z * dOut_11.x;
    float sum_18 = _S356 + (*right_6).primal_0.rows[int(2)].y * dOut_11.y;
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

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S357, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S358, float3  _S359)
{
    _d_mul_2(_S357, _S358, _S359);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S360 = - _s_dOut_2;
    float3  _S361 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S362;
    (&_S362)->primal_0 = (*dpt_0).primal_0;
    (&_S362)->differential_0 = _S361;
    Matrix<float, 3, 3>  _S363 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S364;
    (&_S364)->primal_0 = (*dpR_0).primal_0;
    (&_S364)->differential_0 = _S363;
    s_bwd_prop_mul_0(&_S362, &_S364, _S360);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S362.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S364.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S365, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S366, float3  _S367)
{
    s_bwd_prop_transform_ray_o_0(_S365, _S366, _S367);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_5, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S368 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_5;
    (&dp_R_0)->differential_0 = _S368;
    float3  _S369 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S369;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S370 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S371;
    (&_S371)->primal_0 = (*dpraydir_0).primal_0;
    (&_S371)->differential_0 = _S370;
    Matrix<float, 3, 3>  _S372 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S373;
    (&_S373)->primal_0 = (*dpR_1).primal_0;
    (&_S373)->differential_0 = _S372;
    s_bwd_prop_mul_0(&_S371, &_S373, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S371.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S373.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S374, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S375, float3  _S376)
{
    s_bwd_prop_transform_ray_d_0(_S374, _S375, _S376);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_6, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S377 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_6;
    (&dp_R_1)->differential_0 = _S377;
    float3  _S378 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S378;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S379 = scale_2.x;
    float sx_0 = (F32_exp((_S379)));
    float _S380 = scale_2.y;
    float sy_0 = (F32_exp((_S380)));
    float sz_0 = scale_2.z - 0.5f * (_S379 + _S380);
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
    Matrix<float, 3, 3>  _S381 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S381, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S381, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S381, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
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
    float4  _S382 = normalize_1(quat_6);
    float3  _S383 = exp_0(scale_3);
    float x_20 = _S382.y;
    float x2_4 = x_20 * x_20;
    float y2_4 = _S382.z * _S382.z;
    float z2_4 = _S382.w * _S382.w;
    float xy_4 = _S382.y * _S382.z;
    float xz_4 = _S382.y * _S382.w;
    float yz_4 = _S382.z * _S382.w;
    float wx_4 = _S382.x * _S382.y;
    float wy_4 = _S382.x * _S382.z;
    float wz_4 = _S382.x * _S382.w;
    Matrix<float, 3, 3>  M_2 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S383.x, 0.0f, 0.0f, 0.0f, _S383.y, 0.0f, 0.0f, 0.0f, _S383.z));
    float4  _S384 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S385 = _S384 - floor_0(_S384);
    float4  _S386 = _S385 + make_float4 (dot_1(_S385, float4 {_S385.w, _S385.z, _S385.x, _S385.y} + make_float4 (33.3300018310546875f)));
    float4  _S387 = (float4 {_S386.x, _S386.x, _S386.y, _S386.z} + float4 {_S386.y, _S386.z, _S386.z, _S386.w}) * float4 {_S386.z, _S386.y, _S386.w, _S386.x};
    float4  _S388 = _S387 - floor_0(_S387);
    float2  _S389 = float2 {_S388.x, _S388.z};
    float _S390 = 6.28318548202514648f * _S389.y;
    float2  _S391 = float2 {_S388.y, _S388.w};
    float _S392 = 6.28318548202514648f * _S391.y;
    *mean_1 = *mean_1 + mul_0(mul_3(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S389.x))))))) * make_float2 ((F32_cos((_S390))), (F32_sin((_S390))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S389.x))))))) * make_float2 ((F32_cos((_S390))), (F32_sin((_S390))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S391.x))))))) * make_float2 ((F32_cos((_S392))), (F32_sin((_S392))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_12)
{
    float _S393 = dOut_12.y;
    float _S394 = dOut_12.z;
    float _S395 = dOut_12.x;
    float _S396 = (*a_0).primal_0.z * _S393 + - (*a_0).primal_0.y * _S394;
    float _S397 = - (*a_0).primal_0.z * _S395 + (*a_0).primal_0.x * _S394;
    float _S398 = (*a_0).primal_0.y * _S395 + - (*a_0).primal_0.x * _S393;
    float3  _S399 = make_float3 (- (*b_1).primal_0.z * _S393 + (*b_1).primal_0.y * _S394, (*b_1).primal_0.z * _S395 + - (*b_1).primal_0.x * _S394, - (*b_1).primal_0.y * _S395 + (*b_1).primal_0.x * _S393);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S399;
    float3  _S400 = make_float3 (_S396, _S397, _S398);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S400;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S401 = left_8.y;
    float _S402 = right_8.z;
    float _S403 = left_8.z;
    float _S404 = right_8.y;
    float _S405 = right_8.x;
    float _S406 = left_8.x;
    return make_float3 (_S401 * _S402 - _S403 * _S404, _S403 * _S405 - _S406 * _S402, _S406 * _S404 - _S401 * _S405);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S407 = normalize_1(quat_7);
    float _S408 = scale_4.x;
    float sx_1 = (F32_exp((_S408)));
    float _S409 = scale_4.y;
    float sy_1 = (F32_exp((_S409)));
    float sz_1 = scale_4.z - 0.5f * (_S408 + _S409);
    float x_21 = _S407.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S407.z * _S407.z;
    float z2_5 = _S407.w * _S407.w;
    float xy_5 = _S407.y * _S407.z;
    float xz_5 = _S407.y * _S407.w;
    float yz_5 = _S407.z * _S407.w;
    float wx_5 = _S407.x * _S407.y;
    float wy_5 = _S407.x * _S407.z;
    float wz_5 = _S407.x * _S407.w;
    Matrix<float, 3, 3>  _S410 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S410, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S410, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S410, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S411 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S412 = _S411 - floor_0(_S411);
    float4  _S413 = _S412 + make_float4 (dot_1(_S412, float4 {_S412.w, _S412.z, _S412.x, _S412.y} + make_float4 (33.3300018310546875f)));
    float4  _S414 = (float4 {_S413.x, _S413.x, _S413.y, _S413.z} + float4 {_S413.y, _S413.z, _S413.z, _S413.w}) * float4 {_S413.z, _S413.y, _S413.w, _S413.x};
    float4  _S415 = _S414 - floor_0(_S414);
    float2  _S416 = float2 {_S415.x, _S415.z};
    float _S417 = 6.28318548202514648f * _S416.y;
    float2  _S418 = float2 {_S415.y, _S415.w};
    float _S419 = 6.28318548202514648f * _S418.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S416.x))))))) * make_float2 ((F32_cos((_S417))), (F32_sin((_S417))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S416.x))))))) * make_float2 ((F32_cos((_S417))), (F32_sin((_S417))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S418.x))))))) * make_float2 ((F32_cos((_S419))), (F32_sin((_S419))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_9, float dOut_13)
{
    float _S420 = _slang_select(((*dpx_9).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_9).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S420;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S421 = _slang_select(((*dpx_10).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_10).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S421;
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
    DiffPair_float_0 _S422 = *dpx_11;
    bool _S423;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S423 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S423 = false;
    }
    float _S424;
    if(_S423)
    {
        _S424 = dOut_15;
    }
    else
    {
        _S424 = 0.0f;
    }
    dpx_11->primal_0 = _S422.primal_0;
    dpx_11->differential_0 = _S424;
    DiffPair_float_0 _S425 = *dpMin_0;
    if((_S422.primal_0) < ((*dpMin_0).primal_0))
    {
        _S424 = dOut_15;
    }
    else
    {
        _S424 = 0.0f;
    }
    dpMin_0->primal_0 = _S425.primal_0;
    dpMin_0->differential_0 = _S424;
    DiffPair_float_0 _S426 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S424 = dOut_15;
    }
    else
    {
        _S424 = 0.0f;
    }
    dpMax_0->primal_0 = _S426.primal_0;
    dpMax_0->differential_0 = _S424;
    return;
}

inline __device__ float clamp_0(float x_23, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_23), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_12, float dOut_16)
{
    float _S427 = -0.5f / ((*dpx_12).primal_0 * (F32_sqrt(((*dpx_12).primal_0)))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S427;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_17)
{
    float _S428 = (1.0f - (*dps_0).primal_0) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S428;
    DiffPair_float_0 _S429 = *dpy_3;
    float _S430 = (*dps_0).primal_0 * dOut_17;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S430;
    float _S431 = (_S429.primal_0 - (*dpx_13).primal_0) * dOut_17;
    dps_0->primal_0 = _S429.primal_0;
    dps_0->differential_0 = _S431;
    return;
}

inline __device__ float lerp_0(float x_24, float y_7, float s_4)
{
    return x_24 + (y_7 - x_24) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  * weights_0, FixedArray<float, 23>  * _S432)
{
    float3  _S433;
    bool _S434;
    bool _S435;
    FixedArray<float, 23>  losses_1;
    float _S436 = float(mask_0);
    float3  _S437 = ref_rgb_0 - render_rgb_0;
    float3  _S438 = abs_0(_S437);
    losses_1[int(0)] = (*weights_0)[int(0)] * _S436 * ((_S438.x + _S438.y + _S438.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S436 * clamp_0(dot_0(_S437, _S437) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S439 = float(depth_mask_0 & mask_0);
    float _S440 = _S439 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S441 = _S439 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S440;
    losses_1[int(3)] = _S441;
    losses_1[int(4)] = _S440 * _S440;
    losses_1[int(5)] = _S441 * _S441;
    losses_1[int(6)] = _S440 * _S441;
    bool _S442 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S443 = norm2_0 == 0.0f;
        _S434 = _S443;
        if(_S443)
        {
            _S433 = make_float3 (0.0f);
            break;
        }
        _S433 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S444;
    bool _S445 = !_S434;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S446 = norm2_1 == 0.0f;
        _S435 = _S446;
        if(_S446)
        {
            _S444 = make_float3 (0.0f);
            break;
        }
        _S444 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S447;
    float3  _S448;
    bool _S449 = !_S435;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S448 = make_float3 (0.0f);
            _S447 = false;
            break;
        }
        _S448 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S447 = _S442;
        break;
    }
    float _S450 = float(_S445 & _S447);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S433, _S448);
    losses_1[int(7)] = (*weights_0)[int(2)] * _S450 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S451 = float(_S449 & _S447);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S444, _S448);
    losses_1[int(8)] = (*weights_0)[int(2)] * _S451 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S452 = float(_S445 & _S449);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S433, _S444);
    losses_1[int(11)] = (*weights_0)[int(5)] * _S452 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S453 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S454 = float(alpha_mask_0);
    float _S455 = float(ref_alpha_0);
    float _S456 = (F32_max((_S453), (_S455)));
    losses_1[int(9)] = (*weights_0)[int(3)] * _S454 * - lerp_0((F32_log(((F32_max((1.0f - _S456), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S456), (9.99999997475242708e-07f)))))), _S455);
    float _S457 = 1.0f - _S453;
    float _S458 = 1.0f - _S455;
    float _S459 = (F32_max((_S457), (_S458)));
    losses_1[int(10)] = (*weights_0)[int(4)] * _S454 * - lerp_0((F32_log(((F32_max((1.0f - _S459), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S459), (9.99999997475242708e-07f)))))), _S458);
    losses_1[int(12)] = (*weights_0)[int(6)] * 4.0f * _S453 * _S457;
    float _S460 = (F32_max((_S453), (9.999999960041972e-13f)));
    losses_1[int(13)] = (*weights_0)[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S460;
    losses_1[int(14)] = (*weights_0)[int(8)] * depth_dist_0 / _S460;
    losses_1[int(15)] = (*weights_0)[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S460;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S436;
    losses_1[int(18)] = _S439;
    losses_1[int(19)] = _S450;
    losses_1[int(20)] = _S451;
    losses_1[int(21)] = _S452;
    losses_1[int(22)] = _S454;
    *_S432 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S461, float3  _S462)
{
    return dot_0(_S461, _S462);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S463)
{
    return (F32_rsqrt((_S463)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S464, float _S465, float _S466)
{
    return clamp_0(_S464, _S465, _S466);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S467, DiffPair_float_0 * _S468, DiffPair_float_0 * _S469, float _S470)
{
    _d_lerp_0(_S467, _S468, _S469, _S470);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S471, DiffPair_float_0 * _S472, DiffPair_float_0 * _S473, float _S474)
{
    _d_clamp_0(_S471, _S472, _S473, _S474);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S475, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S476, float _S477)
{
    _d_dot_0(_S475, _S476, _S477);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S478, float _S479)
{
    _d_rsqrt_0(_S478, _S479);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S480, float3  _S481)
{
    _d_abs_vector_0(_S480, _S481);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_float_0 _S482 = *dprender_depth_0;
    DiffPair_float_0 _S483 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S485 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S486 = *dpref_normal_0;
    DiffPair_float_0 _S487 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S488 = *dprgb_dist_0;
    DiffPair_float_0 _S489 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S490 = *dpnormal_dist_0;
    float3  _S491 = make_float3 (0.0f);
    float _S492 = float(mask_1);
    float _S493 = (*weights_1)[int(0)] * _S492;
    float3  _S494 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S495 = s_primal_ctx_dot_0(_S494, _S494) * 0.3333333432674408f;
    float _S496 = float(depth_mask_1 & mask_1);
    float _S497 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S498 = _S496 * s_primal_ctx_log_0(_S497);
    float _S499 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S500 = _S496 * s_primal_ctx_log_0(_S499);
    bool _S501 = normal_mask_1 & mask_1;
    float _S502 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S503 = _S502 == 0.0f;
    float3  _S504;
    if(_S503)
    {
        _S504 = make_float3 (0.0f);
    }
    bool _S505 = !_S503;
    float3  _S506;
    if(_S505)
    {
        float _S507 = s_primal_ctx_rsqrt_0(_S502);
        float3  _S508 = make_float3 (_S507);
        _S504 = _S484.primal_0 * make_float3 (_S507);
        _S506 = _S508;
    }
    else
    {
        _S506 = _S491;
    }
    float _S509 = s_primal_ctx_dot_0(_S485.primal_0, _S485.primal_0);
    bool _S510 = _S509 == 0.0f;
    float3  _S511;
    if(_S510)
    {
        _S511 = make_float3 (0.0f);
    }
    bool _S512 = !_S510;
    float3  _S513;
    if(_S512)
    {
        float _S514 = s_primal_ctx_rsqrt_0(_S509);
        float3  _S515 = make_float3 (_S514);
        _S511 = _S485.primal_0 * make_float3 (_S514);
        _S513 = _S515;
    }
    else
    {
        _S513 = _S491;
    }
    float _S516 = s_primal_ctx_dot_0(_S486.primal_0, _S486.primal_0);
    bool _S517 = _S516 == 0.0f;
    float3  _S518;
    bool _S519;
    if(_S517)
    {
        float3  _S520 = make_float3 (0.0f);
        _S519 = false;
        _S518 = _S520;
    }
    else
    {
        _S519 = _S501;
    }
    bool _S521 = !_S517;
    float3  _S522;
    if(_S521)
    {
        float _S523 = s_primal_ctx_rsqrt_0(_S516);
        float3  _S524 = make_float3 (_S523);
        _S518 = _S486.primal_0 * make_float3 (_S523);
        _S522 = _S524;
    }
    else
    {
        _S522 = _S491;
    }
    float _S525 = (*weights_1)[int(2)] * float(_S505 & _S519);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S504, _S518);
    float _S526 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S527 = (*weights_1)[int(2)] * float(_S512 & _S519);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S511, _S518);
    float _S528 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S529 = (*weights_1)[int(5)] * float(_S505 & _S512);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S504, _S511);
    float _S530 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S531 = s_primal_ctx_clamp_0(_S487.primal_0, 0.0f, 1.0f);
    float _S532 = float(alpha_mask_1);
    float _S533 = (*weights_1)[int(3)] * _S532;
    float _S534 = float(ref_alpha_1);
    float _S535 = s_primal_ctx_max_0(_S531, _S534);
    float _S536 = 1.0f - _S535;
    float _S537 = s_primal_ctx_max_0(_S536, 9.99999997475242708e-07f);
    float _S538 = s_primal_ctx_log_0(_S537);
    float _S539 = s_primal_ctx_max_0(_S535, 9.99999997475242708e-07f);
    float _S540 = s_primal_ctx_log_0(_S539);
    float _S541 = (*weights_1)[int(4)] * _S532;
    float _S542 = 1.0f - _S531;
    float _S543 = 1.0f - _S534;
    float _S544 = s_primal_ctx_max_0(_S542, _S543);
    float _S545 = 1.0f - _S544;
    float _S546 = s_primal_ctx_max_0(_S545, 9.99999997475242708e-07f);
    float _S547 = s_primal_ctx_log_0(_S546);
    float _S548 = s_primal_ctx_max_0(_S544, 9.99999997475242708e-07f);
    float _S549 = s_primal_ctx_log_0(_S548);
    float _S550 = (*weights_1)[int(6)] * 4.0f;
    float _S551 = _S550 * _S531;
    float _S552 = s_primal_ctx_max_0(_S531, 9.999999960041972e-13f);
    float _S553 = _S552 * _S552;
    float _S554 = (*_s_dOut_4)[int(15)] / _S553;
    float _S555 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S552 * _S554));
    float _S556 = (*_s_dOut_4)[int(14)] / _S553;
    float _S557 = (*weights_1)[int(8)] * (_S552 * _S556);
    float _S558 = (*_s_dOut_4)[int(13)] / _S553;
    float _S559 = _S552 * _S558;
    float _S560 = (*weights_1)[int(9)] * ((_S490.primal_0.x + _S490.primal_0.y + _S490.primal_0.z) * 0.3333333432674408f) * - _S554 + (*weights_1)[int(8)] * _S489.primal_0 * - _S556 + (*weights_1)[int(7)] * ((_S488.primal_0.x + _S488.primal_0.y + _S488.primal_0.z) * 0.3333333432674408f) * - _S558;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S531;
    (&_S561)->differential_0 = 0.0f;
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = 9.999999960041972e-13f;
    (&_S562)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S561, &_S562, _S560);
    float _S563 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S559);
    float _S564 = _S551 * (*_s_dOut_4)[int(12)];
    float _S565 = _S550 * (_S542 * (*_s_dOut_4)[int(12)]);
    float _S566 = - (_S541 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S547;
    (&_S567)->differential_0 = 0.0f;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S549;
    (&_S568)->differential_0 = 0.0f;
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S543;
    (&_S569)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S567, &_S568, &_S569, _S566);
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S548;
    (&_S570)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S570, _S568.differential_0);
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S544;
    (&_S571)->differential_0 = 0.0f;
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = 9.99999997475242708e-07f;
    (&_S572)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S571, &_S572, _S570.differential_0);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S546;
    (&_S573)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S573, _S567.differential_0);
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = _S545;
    (&_S574)->differential_0 = 0.0f;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = 9.99999997475242708e-07f;
    (&_S575)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S574, &_S575, _S573.differential_0);
    float _S576 = _S571.differential_0 + - _S574.differential_0;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = _S542;
    (&_S577)->differential_0 = 0.0f;
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = _S543;
    (&_S578)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S577, &_S578, _S576);
    float _S579 = - (_S564 + _S577.differential_0);
    float _S580 = - (_S533 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S581;
    (&_S581)->primal_0 = _S538;
    (&_S581)->differential_0 = 0.0f;
    DiffPair_float_0 _S582;
    (&_S582)->primal_0 = _S540;
    (&_S582)->differential_0 = 0.0f;
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S534;
    (&_S583)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S581, &_S582, &_S583, _S580);
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S539;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S584, _S582.differential_0);
    DiffPair_float_0 _S585;
    (&_S585)->primal_0 = _S535;
    (&_S585)->differential_0 = 0.0f;
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = 9.99999997475242708e-07f;
    (&_S586)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S585, &_S586, _S584.differential_0);
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = _S537;
    (&_S587)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S587, _S581.differential_0);
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = _S536;
    (&_S588)->differential_0 = 0.0f;
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = 9.99999997475242708e-07f;
    (&_S589)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S588, &_S589, _S587.differential_0);
    float _S590 = _S585.differential_0 + - _S588.differential_0;
    DiffPair_float_0 _S591;
    (&_S591)->primal_0 = _S531;
    (&_S591)->differential_0 = 0.0f;
    DiffPair_float_0 _S592;
    (&_S592)->primal_0 = _S534;
    (&_S592)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S591, &_S592, _S590);
    float _S593 = _S561.differential_0 + _S565 + _S579 + _S591.differential_0;
    DiffPair_float_0 _S594;
    (&_S594)->primal_0 = _S487.primal_0;
    (&_S594)->differential_0 = 0.0f;
    DiffPair_float_0 _S595;
    (&_S595)->primal_0 = 0.0f;
    (&_S595)->differential_0 = 0.0f;
    DiffPair_float_0 _S596;
    (&_S596)->primal_0 = 1.0f;
    (&_S596)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S594, &_S595, &_S596, _S593);
    DiffPair_float_0 _S597 = _S594;
    float _S598 = _S529 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S599;
    (&_S599)->primal_0 = _S530;
    (&_S599)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S599, _S598);
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = cos_sim_loss_5;
    (&_S600)->differential_0 = 0.0f;
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = 9.999999960041972e-13f;
    (&_S601)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S600, &_S601, _S599.differential_0);
    float _S602 = 0.5f * - (_S598 + _S600.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S603;
    (&_S603)->primal_0 = _S504;
    (&_S603)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S604;
    (&_S604)->primal_0 = _S511;
    (&_S604)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S603, &_S604, _S602);
    float _S605 = _S527 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S606;
    (&_S606)->primal_0 = _S528;
    (&_S606)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S606, _S605);
    DiffPair_float_0 _S607;
    (&_S607)->primal_0 = cos_sim_loss_4;
    (&_S607)->differential_0 = 0.0f;
    DiffPair_float_0 _S608;
    (&_S608)->primal_0 = 9.999999960041972e-13f;
    (&_S608)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S607, &_S608, _S606.differential_0);
    float _S609 = 0.5f * - (_S605 + _S607.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S610;
    (&_S610)->primal_0 = _S511;
    (&_S610)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S611;
    (&_S611)->primal_0 = _S518;
    (&_S611)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S610, &_S611, _S609);
    float _S612 = _S525 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = _S526;
    (&_S613)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S613, _S612);
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = cos_sim_loss_3;
    (&_S614)->differential_0 = 0.0f;
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = 9.999999960041972e-13f;
    (&_S615)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S614, &_S615, _S613.differential_0);
    float _S616 = 0.5f * - (_S612 + _S614.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S617;
    (&_S617)->primal_0 = _S504;
    (&_S617)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S618;
    (&_S618)->primal_0 = _S518;
    (&_S618)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S617, &_S618, _S616);
    float3  _S619 = _S611.differential_0 + _S618.differential_0;
    float3  _S620 = _S603.differential_0 + _S617.differential_0;
    float3  _S621 = make_float3 (_S555, _S555, _S555);
    float3  _S622 = make_float3 (_S563, _S563, _S563);
    float3  _S623 = _S604.differential_0 + _S610.differential_0;
    float _S624;
    if(_S521)
    {
        float3  _S625 = _S486.primal_0 * _S619;
        float3  _S626 = _S522 * _S619;
        float _S627 = _S625.x + _S625.y + _S625.z;
        DiffPair_float_0 _S628;
        (&_S628)->primal_0 = _S516;
        (&_S628)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S628, _S627);
        _S624 = _S628.differential_0;
        _S504 = _S626;
    }
    else
    {
        _S624 = 0.0f;
        _S504 = _S491;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S629;
    (&_S629)->primal_0 = _S486.primal_0;
    (&_S629)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S630;
    (&_S630)->primal_0 = _S486.primal_0;
    (&_S630)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S629, &_S630, _S624);
    float3  _S631 = _S630.differential_0 + _S629.differential_0 + _S504;
    if(_S512)
    {
        float3  _S632 = _S485.primal_0 * _S623;
        float3  _S633 = _S513 * _S623;
        float _S634 = _S632.x + _S632.y + _S632.z;
        DiffPair_float_0 _S635;
        (&_S635)->primal_0 = _S509;
        (&_S635)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S635, _S634);
        _S624 = _S635.differential_0;
        _S504 = _S633;
    }
    else
    {
        _S624 = 0.0f;
        _S504 = _S491;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S636;
    (&_S636)->primal_0 = _S485.primal_0;
    (&_S636)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S637;
    (&_S637)->primal_0 = _S485.primal_0;
    (&_S637)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S636, &_S637, _S624);
    float3  _S638 = _S637.differential_0 + _S636.differential_0 + _S504;
    if(_S505)
    {
        float3  _S639 = _S484.primal_0 * _S620;
        float3  _S640 = _S506 * _S620;
        float _S641 = _S639.x + _S639.y + _S639.z;
        DiffPair_float_0 _S642;
        (&_S642)->primal_0 = _S502;
        (&_S642)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S642, _S641);
        _S624 = _S642.differential_0;
        _S504 = _S640;
    }
    else
    {
        _S624 = 0.0f;
        _S504 = _S491;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S643;
    (&_S643)->primal_0 = _S484.primal_0;
    (&_S643)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S644;
    (&_S644)->primal_0 = _S484.primal_0;
    (&_S644)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S643, &_S644, _S624);
    float _S645 = _S500 * (*_s_dOut_4)[int(6)];
    float _S646 = _S500 * (*_s_dOut_4)[int(5)];
    float _S647 = _S498 * (*_s_dOut_4)[int(4)];
    float _S648 = _S496 * (_S498 * (*_s_dOut_4)[int(6)] + _S646 + _S646 + (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S649;
    (&_S649)->primal_0 = _S499;
    (&_S649)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S649, _S648);
    DiffPair_float_0 _S650;
    (&_S650)->primal_0 = _S483.primal_0;
    (&_S650)->differential_0 = 0.0f;
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = 0.00009999999747379f;
    (&_S651)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S650, &_S651, _S649.differential_0);
    float _S652 = _S496 * (_S645 + _S647 + _S647 + (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = _S497;
    (&_S653)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S653, _S652);
    DiffPair_float_0 _S654;
    (&_S654)->primal_0 = _S482.primal_0;
    (&_S654)->differential_0 = 0.0f;
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = 0.00009999999747379f;
    (&_S655)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S654, &_S655, _S653.differential_0);
    float _S656 = _S492 * (*_s_dOut_4)[int(1)];
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = _S495;
    (&_S657)->differential_0 = 0.0f;
    DiffPair_float_0 _S658;
    (&_S658)->primal_0 = 0.0f;
    (&_S658)->differential_0 = 0.0f;
    DiffPair_float_0 _S659;
    (&_S659)->primal_0 = 1.0f;
    (&_S659)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S657, &_S658, &_S659, _S656);
    float _S660 = 0.3333333432674408f * _S657.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S661;
    (&_S661)->primal_0 = _S494;
    (&_S661)->differential_0 = _S491;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S662;
    (&_S662)->primal_0 = _S494;
    (&_S662)->differential_0 = _S491;
    s_bwd_prop_dot_0(&_S661, &_S662, _S660);
    float _S663 = 0.3333333432674408f * (_S493 * (*_s_dOut_4)[int(0)]);
    float3  _S664 = make_float3 (_S663, _S663, _S663);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S665;
    (&_S665)->primal_0 = _S494;
    (&_S665)->differential_0 = _S491;
    s_bwd_prop_abs_0(&_S665, _S664);
    float3  _S666 = _S662.differential_0 + _S661.differential_0 + _S665.differential_0;
    float3  _S667 = - _S666;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S621;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S557;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S622;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S597.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S631;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S638;
    float3  _S668 = _S644.differential_0 + _S643.differential_0 + _S504;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S668;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S650.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S654.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S666;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S667;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S669, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S670, DiffPair_float_0 * _S671, DiffPair_float_0 * _S672, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S673, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S674, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S675, DiffPair_float_0 * _S676, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S677, DiffPair_float_0 * _S678, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S679, bool _S680, bool _S681, bool _S682, bool _S683, bool _S684, FixedArray<float, 10>  * _S685, FixedArray<float, 23>  * _S686)
{
    s_bwd_prop_per_pixel_losses_0(_S669, _S670, _S671, _S672, _S673, _S674, _S675, _S676, _S677, _S678, _S679, _S680, _S681, _S682, _S683, _S684, _S685, _S686);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  * weights_2, FixedArray<float, 23>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S687 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S687;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S687;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S687;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S687;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S687;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S687;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S687;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, weights_2, v_losses_0);
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
    float _S688 = 1.0f / ((*dpx_14).primal_0 * 52.30258560180664062f) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S688;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  * raw_losses_0, FixedArray<float, 10>  * weights_3, FixedArray<float, 10>  * _S689)
{
    FixedArray<float, 10>  losses_2;
    float _S690 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*raw_losses_0)[int(0)] / _S690;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S690)));
    bool _S691;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S691 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S691 = false;
    }
    float _S692;
    if(_S691)
    {
        _S692 = (*weights_3)[int(1)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S692 = 0.0f;
    }
    losses_2[int(2)] = _S692;
    losses_2[int(3)] = ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = ((*raw_losses_0)[int(9)] + (*raw_losses_0)[int(10)]) / (F32_max(((*raw_losses_0)[int(22)]), (1.0f)));
    losses_2[int(5)] = (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    float _S693 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(6)] = (*raw_losses_0)[int(12)] / _S693;
    losses_2[int(7)] = (*raw_losses_0)[int(13)] / _S693;
    losses_2[int(8)] = (*raw_losses_0)[int(14)] / _S693;
    losses_2[int(9)] = (*raw_losses_0)[int(15)] / _S693;
    *_S689 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S694)
{
    return (F32_sqrt((_S694)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S695, float _S696)
{
    _d_log10_0(_S695, _S696);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S697 = dpraw_losses_0->primal_0;
    float _S698 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S699 = _S698 * _S698;
    float _S700 = dpraw_losses_0->primal_0[int(1)] / _S698;
    bool _S701 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S702;
    if(_S701)
    {
        _S702 = (_S697[int(3)]) != 0.0f;
    }
    else
    {
        _S702 = false;
    }
    float _S703;
    float _S704;
    float _S705;
    float _S706;
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
    if(_S702)
    {
        float _S718 = _S697[int(2)] * _S697[int(3)];
        float _S719 = _S697[int(18)] * _S697[int(18)];
        float _S720 = _S697[int(6)] - _S718 / _S697[int(18)];
        float _S721 = _S697[int(2)] * _S697[int(2)];
        float _S722 = _S697[int(4)] - _S721 / _S697[int(18)];
        float _S723 = _S697[int(3)] * _S697[int(3)];
        float _S724 = _S697[int(5)] - _S723 / _S697[int(18)];
        float _S725 = _S722 * _S724 + 1.0f;
        float _S726 = s_primal_ctx_max_0(9.999999960041972e-13f, _S725);
        float _S727 = s_primal_ctx_sqrt_0(_S726);
        float _S728 = _S727 * _S727;
        float _S729 = 1.0f - _S720 / _S727;
        _S703 = (*weights_4)[int(1)];
        _S704 = _S729;
        _S705 = _S728;
        _S706 = _S720;
        _S707 = _S727;
        _S708 = _S726;
        _S709 = _S725;
        _S710 = _S722;
        _S711 = _S724;
        _S712 = _S719;
        _S713 = _S723;
        _S714 = _S697[int(3)];
        _S715 = _S721;
        _S716 = _S697[int(2)];
        _S717 = _S718;
    }
    else
    {
        _S703 = 0.0f;
        _S704 = 0.0f;
        _S705 = 0.0f;
        _S706 = 0.0f;
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
    }
    float _S730 = s_primal_ctx_max_0(_S697[int(19)], 1.0f);
    float _S731 = _S730 * _S730;
    float _S732 = s_primal_ctx_max_0(_S697[int(20)], 1.0f);
    float _S733 = _S732 * _S732;
    float _S734 = float((I32_max((int((_S697[int(19)]) > 0.5f) + int((_S697[int(20)]) > 0.5f)), (int(1)))));
    float _S735 = _S697[int(9)] + _S697[int(10)];
    float _S736 = s_primal_ctx_max_0(_S697[int(22)], 1.0f);
    float _S737 = _S736 * _S736;
    float _S738 = s_primal_ctx_max_0(_S697[int(21)], 1.0f);
    float _S739 = _S738 * _S738;
    float _S740 = s_primal_ctx_max_0(_S697[int(16)], 1.0f);
    float _S741 = _S740 * _S740;
    float _S742 = (*_s_dOut_5)[int(9)] / _S741;
    float _S743 = _S740 * _S742;
    float _S744 = (*_s_dOut_5)[int(8)] / _S741;
    float _S745 = _S740 * _S744;
    float _S746 = (*_s_dOut_5)[int(7)] / _S741;
    float _S747 = _S740 * _S746;
    float _S748 = (*_s_dOut_5)[int(6)] / _S741;
    float _S749 = _S740 * _S748;
    float _S750 = _S697[int(15)] * - _S742 + _S697[int(14)] * - _S744 + _S697[int(13)] * - _S746 + _S697[int(12)] * - _S748;
    DiffPair_float_0 _S751;
    (&_S751)->primal_0 = _S697[int(16)];
    (&_S751)->differential_0 = 0.0f;
    DiffPair_float_0 _S752;
    (&_S752)->primal_0 = 1.0f;
    (&_S752)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S751, &_S752, _S750);
    float _S753 = (*_s_dOut_5)[int(5)] / _S739;
    float _S754 = _S697[int(11)] * - _S753;
    float _S755 = _S738 * _S753;
    DiffPair_float_0 _S756;
    (&_S756)->primal_0 = _S697[int(21)];
    (&_S756)->differential_0 = 0.0f;
    DiffPair_float_0 _S757;
    (&_S757)->primal_0 = 1.0f;
    (&_S757)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S756, &_S757, _S754);
    float _S758 = (*_s_dOut_5)[int(4)] / _S737;
    float _S759 = _S735 * - _S758;
    float _S760 = _S736 * _S758;
    DiffPair_float_0 _S761;
    (&_S761)->primal_0 = _S697[int(22)];
    (&_S761)->differential_0 = 0.0f;
    DiffPair_float_0 _S762;
    (&_S762)->primal_0 = 1.0f;
    (&_S762)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S761, &_S762, _S759);
    float _S763 = (*_s_dOut_5)[int(3)] / _S734;
    float _S764 = _S763 / _S733;
    float _S765 = _S697[int(8)] * - _S764;
    float _S766 = _S732 * _S764;
    DiffPair_float_0 _S767;
    (&_S767)->primal_0 = _S697[int(20)];
    (&_S767)->differential_0 = 0.0f;
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = 1.0f;
    (&_S768)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S767, &_S768, _S765);
    float _S769 = _S763 / _S731;
    float _S770 = _S697[int(7)] * - _S769;
    float _S771 = _S730 * _S769;
    DiffPair_float_0 _S772;
    (&_S772)->primal_0 = _S697[int(19)];
    (&_S772)->differential_0 = 0.0f;
    DiffPair_float_0 _S773;
    (&_S773)->primal_0 = 1.0f;
    (&_S773)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S772, &_S773, _S770);
    FixedArray<float, 23>  _S774;
    _S774[int(0)] = 0.0f;
    _S774[int(1)] = 0.0f;
    _S774[int(2)] = 0.0f;
    _S774[int(3)] = 0.0f;
    _S774[int(4)] = 0.0f;
    _S774[int(5)] = 0.0f;
    _S774[int(6)] = 0.0f;
    _S774[int(7)] = 0.0f;
    _S774[int(8)] = 0.0f;
    _S774[int(9)] = 0.0f;
    _S774[int(10)] = 0.0f;
    _S774[int(11)] = 0.0f;
    _S774[int(12)] = 0.0f;
    _S774[int(13)] = 0.0f;
    _S774[int(14)] = 0.0f;
    _S774[int(15)] = 0.0f;
    _S774[int(16)] = 0.0f;
    _S774[int(17)] = 0.0f;
    _S774[int(18)] = 0.0f;
    _S774[int(19)] = 0.0f;
    _S774[int(20)] = 0.0f;
    _S774[int(21)] = 0.0f;
    _S774[int(22)] = 0.0f;
    _S774[int(15)] = _S743;
    _S774[int(14)] = _S745;
    _S774[int(13)] = _S747;
    _S774[int(16)] = _S751.differential_0;
    _S774[int(12)] = _S749;
    _S774[int(21)] = _S756.differential_0;
    _S774[int(11)] = _S755;
    _S774[int(22)] = _S761.differential_0;
    _S774[int(10)] = _S760;
    _S774[int(9)] = _S760;
    _S774[int(20)] = _S767.differential_0;
    _S774[int(8)] = _S766;
    _S774[int(19)] = _S772.differential_0;
    _S774[int(7)] = _S771;
    float _S775 = _S774[int(0)];
    float _S776 = _S774[int(1)];
    float _S777 = _S774[int(2)];
    float _S778 = _S774[int(3)];
    float _S779 = _S774[int(4)];
    float _S780 = _S774[int(5)];
    float _S781 = _S774[int(6)];
    float _S782 = _S774[int(7)];
    float _S783 = _S774[int(8)];
    float _S784 = _S774[int(9)];
    float _S785 = _S774[int(10)];
    float _S786 = _S774[int(11)];
    float _S787 = _S774[int(12)];
    float _S788 = _S774[int(13)];
    float _S789 = _S774[int(14)];
    float _S790 = _S774[int(15)];
    float _S791 = _S774[int(16)];
    float _S792 = _S774[int(17)];
    float _S793 = _S774[int(18)];
    float _S794 = _S774[int(19)];
    float _S795 = _S774[int(20)];
    float _S796 = _S774[int(21)];
    float _S797 = _S774[int(22)];
    FixedArray<float, 23>  _S798;
    if(_S702)
    {
        float _S799 = _S703 * (*_s_dOut_5)[int(2)];
        DiffPair_float_0 _S800;
        (&_S800)->primal_0 = _S704;
        (&_S800)->differential_0 = 0.0f;
        DiffPair_float_0 _S801;
        (&_S801)->primal_0 = 0.0f;
        (&_S801)->differential_0 = 0.0f;
        DiffPair_float_0 _S802;
        (&_S802)->primal_0 = 2.0f;
        (&_S802)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S800, &_S801, &_S802, _S799);
        float _S803 = - _S800.differential_0 / _S705;
        float _S804 = _S706 * - _S803;
        float _S805 = _S707 * _S803;
        DiffPair_float_0 _S806;
        (&_S806)->primal_0 = _S708;
        (&_S806)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S806, _S804);
        DiffPair_float_0 _S807;
        (&_S807)->primal_0 = 9.999999960041972e-13f;
        (&_S807)->differential_0 = 0.0f;
        DiffPair_float_0 _S808;
        (&_S808)->primal_0 = _S709;
        (&_S808)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S807, &_S808, _S806.differential_0);
        float _S809 = _S710 * _S808.differential_0;
        float _S810 = _S711 * _S808.differential_0;
        float _S811 = - _S809 / _S712;
        float _S812 = _S714 * (_S697[int(18)] * _S811);
        float _S813 = - _S810 / _S712;
        float _S814 = _S716 * (_S697[int(18)] * _S813);
        float _S815 = - _S805 / _S712;
        float _S816 = _S697[int(18)] * _S815;
        float _S817 = _S812 + _S812 + _S716 * _S816;
        float _S818 = _S814 + _S814 + _S714 * _S816;
        float _S819 = _S713 * - _S811 + _S715 * - _S813 + _S717 * - _S815;
        FixedArray<float, 23>  _S820;
        _S820[int(0)] = 0.0f;
        _S820[int(1)] = 0.0f;
        _S820[int(2)] = 0.0f;
        _S820[int(3)] = 0.0f;
        _S820[int(4)] = 0.0f;
        _S820[int(5)] = 0.0f;
        _S820[int(6)] = 0.0f;
        _S820[int(7)] = 0.0f;
        _S820[int(8)] = 0.0f;
        _S820[int(9)] = 0.0f;
        _S820[int(10)] = 0.0f;
        _S820[int(11)] = 0.0f;
        _S820[int(12)] = 0.0f;
        _S820[int(13)] = 0.0f;
        _S820[int(14)] = 0.0f;
        _S820[int(15)] = 0.0f;
        _S820[int(16)] = 0.0f;
        _S820[int(17)] = 0.0f;
        _S820[int(18)] = 0.0f;
        _S820[int(19)] = 0.0f;
        _S820[int(20)] = 0.0f;
        _S820[int(21)] = 0.0f;
        _S820[int(22)] = 0.0f;
        _S820[int(5)] = _S809;
        _S820[int(4)] = _S810;
        _S820[int(3)] = _S817;
        _S820[int(2)] = _S818;
        _S820[int(6)] = _S805;
        float _S821 = _S776 + _S820[int(1)];
        float _S822 = _S777 + _S820[int(2)];
        float _S823 = _S778 + _S820[int(3)];
        float _S824 = _S779 + _S820[int(4)];
        float _S825 = _S780 + _S820[int(5)];
        float _S826 = _S781 + _S820[int(6)];
        float _S827 = _S782 + _S820[int(7)];
        float _S828 = _S783 + _S820[int(8)];
        float _S829 = _S784 + _S820[int(9)];
        float _S830 = _S785 + _S820[int(10)];
        float _S831 = _S786 + _S820[int(11)];
        float _S832 = _S787 + _S820[int(12)];
        float _S833 = _S788 + _S820[int(13)];
        float _S834 = _S789 + _S820[int(14)];
        float _S835 = _S790 + _S820[int(15)];
        float _S836 = _S791 + _S820[int(16)];
        float _S837 = _S792 + _S820[int(17)];
        float _S838 = _S793 + _S820[int(18)];
        float _S839 = _S794 + _S820[int(19)];
        float _S840 = _S795 + _S820[int(20)];
        float _S841 = _S796 + _S820[int(21)];
        float _S842 = _S797 + _S820[int(22)];
        _S798[int(0)] = _S775 + _S820[int(0)];
        _S798[int(1)] = _S821;
        _S798[int(2)] = _S822;
        _S798[int(3)] = _S823;
        _S798[int(4)] = _S824;
        _S798[int(5)] = _S825;
        _S798[int(6)] = _S826;
        _S798[int(7)] = _S827;
        _S798[int(8)] = _S828;
        _S798[int(9)] = _S829;
        _S798[int(10)] = _S830;
        _S798[int(11)] = _S831;
        _S798[int(12)] = _S832;
        _S798[int(13)] = _S833;
        _S798[int(14)] = _S834;
        _S798[int(15)] = _S835;
        _S798[int(16)] = _S836;
        _S798[int(17)] = _S837;
        _S798[int(18)] = _S838;
        _S798[int(19)] = _S839;
        _S798[int(20)] = _S840;
        _S798[int(21)] = _S841;
        _S798[int(22)] = _S842;
        _S703 = _S819;
    }
    else
    {
        _S798[int(0)] = _S775;
        _S798[int(1)] = _S776;
        _S798[int(2)] = _S777;
        _S798[int(3)] = _S778;
        _S798[int(4)] = _S779;
        _S798[int(5)] = _S780;
        _S798[int(6)] = _S781;
        _S798[int(7)] = _S782;
        _S798[int(8)] = _S783;
        _S798[int(9)] = _S784;
        _S798[int(10)] = _S785;
        _S798[int(11)] = _S786;
        _S798[int(12)] = _S787;
        _S798[int(13)] = _S788;
        _S798[int(14)] = _S789;
        _S798[int(15)] = _S790;
        _S798[int(16)] = _S791;
        _S798[int(17)] = _S792;
        _S798[int(18)] = _S793;
        _S798[int(19)] = _S794;
        _S798[int(20)] = _S795;
        _S798[int(21)] = _S796;
        _S798[int(22)] = _S797;
        _S703 = 0.0f;
    }
    if(_S701)
    {
        FixedArray<float, 23>  _S843;
        _S843[int(0)] = 0.0f;
        _S843[int(1)] = 0.0f;
        _S843[int(2)] = 0.0f;
        _S843[int(3)] = 0.0f;
        _S843[int(4)] = 0.0f;
        _S843[int(5)] = 0.0f;
        _S843[int(6)] = 0.0f;
        _S843[int(7)] = 0.0f;
        _S843[int(8)] = 0.0f;
        _S843[int(9)] = 0.0f;
        _S843[int(10)] = 0.0f;
        _S843[int(11)] = 0.0f;
        _S843[int(12)] = 0.0f;
        _S843[int(13)] = 0.0f;
        _S843[int(14)] = 0.0f;
        _S843[int(15)] = 0.0f;
        _S843[int(16)] = 0.0f;
        _S843[int(17)] = 0.0f;
        _S843[int(18)] = 0.0f;
        _S843[int(19)] = 0.0f;
        _S843[int(20)] = 0.0f;
        _S843[int(21)] = 0.0f;
        _S843[int(22)] = 0.0f;
        _S843[int(3)] = 0.0f;
        float _S844 = _S798[int(1)] + _S843[int(1)];
        float _S845 = _S798[int(2)] + _S843[int(2)];
        float _S846 = _S798[int(3)] + _S843[int(3)];
        float _S847 = _S798[int(4)] + _S843[int(4)];
        float _S848 = _S798[int(5)] + _S843[int(5)];
        float _S849 = _S798[int(6)] + _S843[int(6)];
        float _S850 = _S798[int(7)] + _S843[int(7)];
        float _S851 = _S798[int(8)] + _S843[int(8)];
        float _S852 = _S798[int(9)] + _S843[int(9)];
        float _S853 = _S798[int(10)] + _S843[int(10)];
        float _S854 = _S798[int(11)] + _S843[int(11)];
        float _S855 = _S798[int(12)] + _S843[int(12)];
        float _S856 = _S798[int(13)] + _S843[int(13)];
        float _S857 = _S798[int(14)] + _S843[int(14)];
        float _S858 = _S798[int(15)] + _S843[int(15)];
        float _S859 = _S798[int(16)] + _S843[int(16)];
        float _S860 = _S798[int(17)] + _S843[int(17)];
        float _S861 = _S798[int(18)] + _S843[int(18)];
        float _S862 = _S798[int(19)] + _S843[int(19)];
        float _S863 = _S798[int(20)] + _S843[int(20)];
        float _S864 = _S798[int(21)] + _S843[int(21)];
        float _S865 = _S798[int(22)] + _S843[int(22)];
        _S798[int(0)] = _S798[int(0)] + _S843[int(0)];
        _S798[int(1)] = _S844;
        _S798[int(2)] = _S845;
        _S798[int(3)] = _S846;
        _S798[int(4)] = _S847;
        _S798[int(5)] = _S848;
        _S798[int(6)] = _S849;
        _S798[int(7)] = _S850;
        _S798[int(8)] = _S851;
        _S798[int(9)] = _S852;
        _S798[int(10)] = _S853;
        _S798[int(11)] = _S854;
        _S798[int(12)] = _S855;
        _S798[int(13)] = _S856;
        _S798[int(14)] = _S857;
        _S798[int(15)] = _S858;
        _S798[int(16)] = _S859;
        _S798[int(17)] = _S860;
        _S798[int(18)] = _S861;
        _S798[int(19)] = _S862;
        _S798[int(20)] = _S863;
        _S798[int(21)] = _S864;
        _S798[int(22)] = _S865;
    }
    float _S866 = -10.0f * (*_s_dOut_5)[int(1)];
    DiffPair_float_0 _S867;
    (&_S867)->primal_0 = _S700;
    (&_S867)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S867, _S866);
    float _S868 = _S867.differential_0 / _S699;
    float _S869 = _S698 * _S868;
    float _S870 = (*_s_dOut_5)[int(0)] / _S699;
    float _S871 = _S698 * _S870;
    float _S872 = _S697[int(1)] * - _S868 + _S697[int(0)] * - _S870;
    DiffPair_float_0 _S873;
    (&_S873)->primal_0 = _S697[int(17)];
    (&_S873)->differential_0 = 0.0f;
    DiffPair_float_0 _S874;
    (&_S874)->primal_0 = 1.0f;
    (&_S874)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S873, &_S874, _S872);
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
    _S875[int(18)] = _S703;
    _S875[int(1)] = _S869;
    _S875[int(17)] = _S873.differential_0;
    _S875[int(0)] = _S871;
    FixedArray<float, 23>  _S876 = {
        _S798[int(0)] + _S875[int(0)], _S798[int(1)] + _S875[int(1)], _S798[int(2)] + _S875[int(2)], _S798[int(3)] + _S875[int(3)], _S798[int(4)] + _S875[int(4)], _S798[int(5)] + _S875[int(5)], _S798[int(6)] + _S875[int(6)], _S798[int(7)] + _S875[int(7)], _S798[int(8)] + _S875[int(8)], _S798[int(9)] + _S875[int(9)], _S798[int(10)] + _S875[int(10)], _S798[int(11)] + _S875[int(11)], _S798[int(12)] + _S875[int(12)], _S798[int(13)] + _S875[int(13)], _S798[int(14)] + _S875[int(14)], _S798[int(15)] + _S875[int(15)], _S798[int(16)] + _S875[int(16)], _S798[int(17)] + _S875[int(17)], _S798[int(18)] + _S875[int(18)], _S798[int(19)] + _S875[int(19)], _S798[int(20)] + _S875[int(20)], _S798[int(21)] + _S875[int(21)], _S798[int(22)] + _S875[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S876;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S877, FixedArray<float, 10>  * _S878, FixedArray<float, 10>  * _S879)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S877, _S878, _S879);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  * raw_losses_1, FixedArray<float, 10>  * weights_5, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 23>  * _S880)
{
    FixedArray<float, 23>  _S881 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S881;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_5, v_losses_1);
    *_S880 = (&dp_raw_losses_0)->differential_0;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S882, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S883, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S884, float3  _S885)
{
    _d_clamp_vector_0(_S882, _S883, _S884, _S885);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S886 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S887 = make_float3 (_S886);
    float3  _S888 = make_float3 (0.0f);
    float3  _S889 = make_float3 (1.0f);
    float3  _S890 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S891;
    (&_S891)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S886) * (*dpbackground_0).primal_0;
    (&_S891)->differential_0 = _S890;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S892;
    (&_S892)->primal_0 = _S888;
    (&_S892)->differential_0 = _S890;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S893;
    (&_S893)->primal_0 = _S889;
    (&_S893)->differential_0 = _S890;
    s_bwd_prop_clamp_1(&_S891, &_S892, &_S893, _s_dOut_6);
    float3  _S894 = _S887 * _S891.differential_0;
    float3  _S895 = (*dpbackground_0).primal_0 * _S891.differential_0;
    float _S896 = - (_S895.x + _S895.y + _S895.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S894;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S896;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S891.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S897, DiffPair_float_0 * _S898, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S899, float3  _S900)
{
    s_bwd_prop_blend_background_0(_S897, _S898, _S899, _S900);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S901 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S901;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S901;
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
        DiffPair_float_0 _S902 = *dpx_16;
        float _S903 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_20;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S903;
        float _S904 = val_0 * (F32_log((_S902.primal_0))) * dOut_20;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S904;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S905 = rgb_2;
    float _S906;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S906 = _S905.x * 12.92000007629394531f;
    }
    else
    {
        _S906 = 1.0549999475479126f * (F32_pow((_S905.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S905)->x) = _S906;
    if((_S905.y) < 0.00313080009073019f)
    {
        _S906 = _S905.y * 12.92000007629394531f;
    }
    else
    {
        _S906 = 1.0549999475479126f * (F32_pow((_S905.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S905)->y) = _S906;
    if((_S905.z) < 0.00313080009073019f)
    {
        _S906 = _S905.z * 12.92000007629394531f;
    }
    else
    {
        _S906 = 1.0549999475479126f * (F32_pow((_S905.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S905)->z) = _S906;
    return _S905;
}

inline __device__ float s_primal_ctx_pow_0(float _S907, float _S908)
{
    return (F32_pow((_S907), (_S908)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S909, DiffPair_float_0 * _S910, float _S911)
{
    _d_pow_0(_S909, _S910, _S911);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_7)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S912 = *dprgb_1;
    float _S913 = (*dprgb_1).primal_0.x;
    bool _S914 = _S913 < 0.00313080009073019f;
    float _S915;
    if(_S914)
    {
        _S915 = _S913 * 12.92000007629394531f;
    }
    else
    {
        _S915 = 1.0549999475479126f * s_primal_ctx_pow_0(_S913, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S916 = _S912.primal_0;
    *&((&_S916)->x) = _S915;
    float _S917 = _S916.y;
    bool _S918 = _S917 < 0.00313080009073019f;
    if(_S918)
    {
        _S915 = _S917 * 12.92000007629394531f;
    }
    else
    {
        _S915 = 1.0549999475479126f * s_primal_ctx_pow_0(_S917, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S916)->y) = _S915;
    float _S919 = _S916.z;
    bool _S920 = _S919 < 0.00313080009073019f;
    _S916 = _s_dOut_7;
    *&((&_S916)->z) = 0.0f;
    if(_S920)
    {
        _S915 = 12.92000007629394531f * _s_dOut_7.z;
    }
    else
    {
        float _S921 = 1.0549999475479126f * _s_dOut_7.z;
        DiffPair_float_0 _S922;
        (&_S922)->primal_0 = _S919;
        (&_S922)->differential_0 = 0.0f;
        DiffPair_float_0 _S923;
        (&_S923)->primal_0 = 0.4166666567325592f;
        (&_S923)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S922, &_S923, _S921);
        _S915 = _S922.differential_0;
    }
    float3  _S924 = _S916 + make_float3 (0.0f, 0.0f, _S915);
    _S916 = _S924;
    *&((&_S916)->y) = 0.0f;
    if(_S918)
    {
        _S915 = 12.92000007629394531f * _S924.y;
    }
    else
    {
        float _S925 = 1.0549999475479126f * _S924.y;
        DiffPair_float_0 _S926;
        (&_S926)->primal_0 = _S917;
        (&_S926)->differential_0 = 0.0f;
        DiffPair_float_0 _S927;
        (&_S927)->primal_0 = 0.4166666567325592f;
        (&_S927)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S926, &_S927, _S925);
        _S915 = _S926.differential_0;
    }
    float3  _S928 = _S916 + make_float3 (0.0f, _S915, 0.0f);
    _S916 = _S928;
    *&((&_S916)->x) = 0.0f;
    if(_S914)
    {
        _S915 = 12.92000007629394531f * _S928.x;
    }
    else
    {
        float _S929 = 1.0549999475479126f * _S928.x;
        DiffPair_float_0 _S930;
        (&_S930)->primal_0 = _S913;
        (&_S930)->differential_0 = 0.0f;
        DiffPair_float_0 _S931;
        (&_S931)->primal_0 = 0.4166666567325592f;
        (&_S931)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S930, &_S931, _S929);
        _S915 = _S930.differential_0;
    }
    float3  _S932 = _S916 + make_float3 (_S915, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S932;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S933, float3  _S934)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S933, _S934);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S935 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S935;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_7, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S936 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = _S936;
    bool _S937 = undistort_point_0(_S936, dist_coeffs_7, int(12), &uv_7);
    if(!_S937)
    {
        int3  _S938 = make_int3 (int(0));
        float3  _S939 = make_float3 ((float)_S938.x, (float)_S938.y, (float)_S938.z);
        return _S939;
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

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  * points_0)
{
    float3  normal_0 = cross_0((*points_0)[int(1)] - (*points_0)[int(0)], - ((*points_0)[int(3)] - (*points_0)[int(2)]));
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S940, float3  _S941)
{
    return cross_0(_S940, _S941);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_8)
{
    float _S942 = (*dpx_17).primal_0.x;
    float _S943 = (*dpx_17).primal_0.y;
    float _S944 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S942 * _S942 + _S943 * _S943 + _S944 * _S944;
    (&_S945)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S945, _s_dOut_8);
    float _S946 = (*dpx_17).primal_0.z * _S945.differential_0;
    float _S947 = _S946 + _S946;
    float _S948 = (*dpx_17).primal_0.y * _S945.differential_0;
    float _S949 = _S948 + _S948;
    float _S950 = (*dpx_17).primal_0.x * _S945.differential_0;
    float _S951 = _S950 + _S950;
    float3  _S952 = make_float3 (0.0f);
    *&((&_S952)->z) = _S947;
    *&((&_S952)->y) = _S949;
    *&((&_S952)->x) = _S951;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S952;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S953, float _S954)
{
    s_bwd_prop_length_impl_1(_S953, _S954);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S955, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S956, float3  _S957)
{
    _d_cross_0(_S955, _S956, _S957);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S958 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S959 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S960 = s_primal_ctx_cross_0(dx_0, _S959);
    bool _S961 = (s_primal_ctx_dot_0(_S960, _S960)) != 0.0f;
    float3  _S962;
    float3  _S963;
    if(_S961)
    {
        float _S964 = length_2(_S960);
        float3  _S965 = make_float3 (_S964);
        _S962 = make_float3 (_S964 * _S964);
        _S963 = _S965;
    }
    else
    {
        _S962 = _S958;
        _S963 = _S958;
    }
    if(_S961)
    {
        float3  _S966 = _s_dOut_9 / _S962;
        float3  _S967 = _S960 * - _S966;
        float3  _S968 = _S963 * _S966;
        float _S969 = _S967.x + _S967.y + _S967.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S970;
        (&_S970)->primal_0 = _S960;
        (&_S970)->differential_0 = _S958;
        s_bwd_length_impl_1(&_S970, _S969);
        _S962 = _S968 + _S970.differential_0;
    }
    else
    {
        _S962 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S971;
    (&_S971)->primal_0 = _S960;
    (&_S971)->differential_0 = _S958;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S972;
    (&_S972)->primal_0 = _S960;
    (&_S972)->differential_0 = _S958;
    s_bwd_prop_dot_0(&_S971, &_S972, 0.0f);
    float3  _S973 = _S972.differential_0 + _S971.differential_0 + _S962;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S974;
    (&_S974)->primal_0 = dx_0;
    (&_S974)->differential_0 = _S958;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S975;
    (&_S975)->primal_0 = _S959;
    (&_S975)->differential_0 = _S958;
    s_bwd_prop_cross_0(&_S974, &_S975, _S973);
    float3  s_diff_dy_T_0 = - _S975.differential_0;
    float3  _S976 = - s_diff_dy_T_0;
    float3  _S977 = - _S974.differential_0;
    FixedArray<float3 , 4>  _S978;
    _S978[int(0)] = _S958;
    _S978[int(1)] = _S958;
    _S978[int(2)] = _S958;
    _S978[int(3)] = _S958;
    _S978[int(2)] = _S976;
    _S978[int(3)] = s_diff_dy_T_0;
    _S978[int(0)] = _S977;
    _S978[int(1)] = _S974.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S978;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S979, float3  _S980)
{
    s_bwd_prop_points_to_normal_0(_S979, _S980);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  * points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S981 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = *points_1;
    (&dp_points_0)->differential_0 = _S981;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S982 = float2 {intrins_1.z, intrins_1.w};
    float2  _S983 = float2 {intrins_1.x, intrins_1.y};
    float2  _S984 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S982) / _S983;
    float2  uv_8 = _S984;
    bool _S985 = undistort_point_0(_S984, dist_coeffs_8, int(12), &uv_8);
    if(!_S985)
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
    float2  _S986 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S982) / _S983;
    float2  uv_9 = _S986;
    bool _S987 = undistort_point_0(_S986, dist_coeffs_8, int(12), &uv_9);
    if(!_S987)
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
    float2  _S988 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S982) / _S983;
    float2  uv_10 = _S988;
    bool _S989 = undistort_point_0(_S988, dist_coeffs_8, int(12), &uv_10);
    if(!_S989)
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
    float2  _S990 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S982) / _S983;
    float2  uv_11 = _S990;
    bool _S991 = undistort_point_0(_S990, dist_coeffs_8, int(12), &uv_11);
    if(!_S991)
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
    float2  _S992;
    bool _S993;
    float2  _S994;
    bool _S995;
    float2  _S996;
    bool _S997;
    float2  _S998;
    bool _S999;
};

inline __device__ float s_primal_ctx_sin_0(float _S1000)
{
    return (F32_sin((_S1000)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1001)
{
    return (F32_cos((_S1001)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1002 = make_float2 (0.0f);
    _s_diff_ctx_0->_S992 = _S1002;
    _s_diff_ctx_0->_S993 = false;
    _s_diff_ctx_0->_S994 = _S1002;
    _s_diff_ctx_0->_S995 = false;
    _s_diff_ctx_0->_S996 = _S1002;
    _s_diff_ctx_0->_S997 = false;
    _s_diff_ctx_0->_S998 = _S1002;
    _s_diff_ctx_0->_S999 = false;
    _s_diff_ctx_0->_S994 = _S1002;
    _s_diff_ctx_0->_S995 = false;
    _s_diff_ctx_0->_S996 = _S1002;
    _s_diff_ctx_0->_S997 = false;
    _s_diff_ctx_0->_S998 = _S1002;
    _s_diff_ctx_0->_S999 = false;
    float3  _S1003 = make_float3 (0.0f);
    float2  _S1004 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1005 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1006 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1004) / _S1005;
    float2  _S1007 = _S1006;
    bool _S1008 = undistort_point_0(_S1006, dist_coeffs_9, int(12), &_S1007);
    _s_diff_ctx_0->_S992 = _S1007;
    _s_diff_ctx_0->_S993 = _S1008;
    float2  uv_12 = _S1007;
    bool _S1009 = !_S1008;
    float3  normal_4;
    if(_S1009)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1010 = !_S1009;
    int _S1011;
    FixedArray<float3 , 4>  points_3;
    if(_S1010)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1012 = length_1(uv_12);
            float3  raydir_19 = make_float3 ((uv_12 / make_float2 (s_primal_ctx_max_0(_S1012, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1012))).x, (uv_12 / make_float2 (s_primal_ctx_max_0(_S1012, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1012))).y, s_primal_ctx_cos_0(_S1012));
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
        float3  _S1013 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1014 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1004) / _S1005;
        float2  _S1015 = _S1014;
        bool _S1016 = undistort_point_0(_S1014, dist_coeffs_9, int(12), &_S1015);
        _s_diff_ctx_0->_S994 = _S1015;
        _s_diff_ctx_0->_S995 = _S1016;
        float2  uv_13 = _S1015;
        bool _S1017 = !_S1016;
        if(_S1017)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1018 = !_S1017;
        if(_S1018)
        {
            if(is_fisheye_6)
            {
                float _S1019 = length_1(uv_13);
                float3  raydir_21 = make_float3 ((uv_13 / make_float2 (s_primal_ctx_max_0(_S1019, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1019))).x, (uv_13 / make_float2 (s_primal_ctx_max_0(_S1019, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1019))).y, s_primal_ctx_cos_0(_S1019));
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
            float3  _S1020 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1011 = int(2);
            points_3[int(0)] = _S1013;
            points_3[int(1)] = _S1020;
            points_3[int(2)] = _S1003;
            points_3[int(3)] = _S1003;
        }
        else
        {
            _S1011 = int(0);
            points_3[int(0)] = _S1013;
            points_3[int(1)] = _S1003;
            points_3[int(2)] = _S1003;
            points_3[int(3)] = _S1003;
        }
        bool _runFlag_0;
        if(_S1011 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1010;
            _S1011 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1021 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1004) / _S1005;
            float2  _S1022 = _S1021;
            bool _S1023 = undistort_point_0(_S1021, dist_coeffs_9, int(12), &_S1022);
            _s_diff_ctx_0->_S996 = _S1022;
            _s_diff_ctx_0->_S997 = _S1023;
            float2  uv_14 = _S1022;
            if(!_S1023)
            {
                float3  _S1024 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1011 = int(0);
                normal_4 = _S1024;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1025 = length_1(uv_14);
                    float3  raydir_23 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S1025, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1025))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S1025, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1025))).y, s_primal_ctx_cos_0(_S1025));
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
                float2  _S1026 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1004) / _S1005;
                float2  _S1027 = _S1026;
                bool _S1028 = undistort_point_0(_S1026, dist_coeffs_9, int(12), &_S1027);
                _s_diff_ctx_0->_S998 = _S1027;
                _s_diff_ctx_0->_S999 = _S1028;
                float2  uv_15 = _S1027;
                bool _S1029 = !_S1028;
                if(_S1029)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1030 = !_S1029;
                int _S1031;
                if(_S1030)
                {
                    if(is_fisheye_6)
                    {
                        float _S1032 = length_1(uv_15);
                        float3  raydir_25 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S1032, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1032))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S1032, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1032))).y, s_primal_ctx_cos_0(_S1032));
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
                    _S1031 = int(2);
                }
                else
                {
                    _S1031 = int(0);
                }
                if(_S1031 != int(2))
                {
                    _runFlag_0 = false;
                    _S1011 = _S1031;
                }
                if(_runFlag_0)
                {
                    _S1011 = int(1);
                }
            }
        }
    }
    else
    {
        _S1011 = int(0);
        points_3[int(0)] = _S1003;
        points_3[int(1)] = _S1003;
        points_3[int(2)] = _S1003;
        points_3[int(3)] = _S1003;
    }
    if(!(_S1011 != int(1)))
    {
        float3  _S1033 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1033, _S1033)) != 0.0f)
        {
            normal_4 = _S1033 / make_float3 (length_2(_S1033));
        }
        else
        {
            normal_4 = _S1033;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1034 = *dpdepths_1;
    float3  _S1035 = make_float3 (0.0f);
    bool _S1036 = !!_s_diff_ctx_1->_S993;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1037;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1038;
    if(_S1036)
    {
        if(is_fisheye_7)
        {
            float _S1039 = length_1(_s_diff_ctx_1->_S992);
            float3  raydir_31 = make_float3 ((_s_diff_ctx_1->_S992 / make_float2 (s_primal_ctx_max_0(_S1039, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1039))).x, (_s_diff_ctx_1->_S992 / make_float2 (s_primal_ctx_max_0(_S1039, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1039))).y, s_primal_ctx_cos_0(_S1039));
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
            float3  raydir_32 = make_float3 (_s_diff_ctx_1->_S992.x, _s_diff_ctx_1->_S992.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1040 = make_float3 (_S1034.primal_0.x) * raydir_27;
        bool _S1041 = !!_s_diff_ctx_1->_S995;
        if(_S1041)
        {
            if(is_fisheye_7)
            {
                float _S1042 = length_1(_s_diff_ctx_1->_S994);
                float3  raydir_33 = make_float3 ((_s_diff_ctx_1->_S994 / make_float2 (s_primal_ctx_max_0(_S1042, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1042))).x, (_s_diff_ctx_1->_S994 / make_float2 (s_primal_ctx_max_0(_S1042, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1042))).y, s_primal_ctx_cos_0(_S1042));
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
                float3  raydir_34 = make_float3 (_s_diff_ctx_1->_S994.x, _s_diff_ctx_1->_S994.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1043 = make_float3 (_S1034.primal_0.y) * raydir_28;
            _S1037 = int(2);
            points_4[int(0)] = _S1040;
            points_4[int(1)] = _S1043;
            points_4[int(2)] = _S1035;
            points_4[int(3)] = _S1035;
        }
        else
        {
            _S1037 = int(0);
            points_4[int(0)] = _S1040;
            points_4[int(1)] = _S1035;
            points_4[int(2)] = _S1035;
            points_4[int(3)] = _S1035;
            raydir_28 = _S1035;
        }
        if(_S1037 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1036;
            _S1037 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_1->_S997)
            {
                _runFlag_2 = false;
                _S1037 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1044 = length_1(_s_diff_ctx_1->_S996);
                    float3  raydir_35 = make_float3 ((_s_diff_ctx_1->_S996 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).x, (_s_diff_ctx_1->_S996 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).y, s_primal_ctx_cos_0(_S1044));
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
                    float3  raydir_36 = make_float3 (_s_diff_ctx_1->_S996.x, _s_diff_ctx_1->_S996.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1034.primal_0.z) * raydir_29;
                bool _S1045 = !!_s_diff_ctx_1->_S999;
                int _S1046;
                if(_S1045)
                {
                    if(is_fisheye_7)
                    {
                        float _S1047 = length_1(_s_diff_ctx_1->_S998);
                        float3  raydir_37 = make_float3 ((_s_diff_ctx_1->_S998 / make_float2 (s_primal_ctx_max_0(_S1047, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1047))).x, (_s_diff_ctx_1->_S998 / make_float2 (s_primal_ctx_max_0(_S1047, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1047))).y, s_primal_ctx_cos_0(_S1047));
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
                        float3  raydir_38 = make_float3 (_s_diff_ctx_1->_S998.x, _s_diff_ctx_1->_S998.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1034.primal_0.w) * raydir_30;
                    _S1046 = int(2);
                }
                else
                {
                    _S1046 = int(0);
                    raydir_30 = _S1035;
                }
                if(_S1046 != int(2))
                {
                    _runFlag_3 = false;
                    _S1037 = _S1046;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1037 = int(1);
                }
                float3  _S1048 = raydir_29;
                _runFlag_3 = _S1045;
                raydir_29 = raydir_30;
                raydir_30 = _S1048;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1035;
                raydir_30 = _S1035;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1035;
            raydir_30 = _S1035;
        }
        float3  _S1049 = raydir_27;
        float3  _S1050 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1038 = _S1041;
        raydir_29 = _S1050;
        raydir_30 = _S1049;
    }
    else
    {
        _S1037 = int(0);
        points_4[int(0)] = _S1035;
        points_4[int(1)] = _S1035;
        points_4[int(2)] = _S1035;
        points_4[int(3)] = _S1035;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1035;
        raydir_28 = _S1035;
        _S1038 = false;
        raydir_29 = _S1035;
        raydir_30 = _S1035;
    }
    bool _S1051 = !(_S1037 != int(1));
    float3  _S1052;
    float3  _S1053;
    float3  _S1054;
    float3  _S1055;
    float3  _S1056;
    bool _S1057;
    if(_S1051)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1058 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1059 = s_primal_ctx_cross_0(dx_1, _S1058);
        bool _S1060 = (s_primal_ctx_dot_0(_S1059, _S1059)) != 0.0f;
        if(_S1060)
        {
            float _S1061 = length_2(_S1059);
            float3  _S1062 = make_float3 (_S1061);
            _S1052 = make_float3 (_S1061 * _S1061);
            _S1053 = _S1062;
        }
        else
        {
            _S1052 = _S1035;
            _S1053 = _S1035;
        }
        float3  _S1063 = _S1053;
        _S1057 = _S1060;
        _S1053 = _S1059;
        _S1054 = _S1063;
        _S1055 = dx_1;
        _S1056 = _S1058;
    }
    else
    {
        _S1057 = false;
        _S1052 = _S1035;
        _S1053 = _S1035;
        _S1054 = _S1035;
        _S1055 = _S1035;
        _S1056 = _S1035;
    }
    float4  _S1064 = make_float4 (0.0f);
    if(_S1051)
    {
        if(_S1057)
        {
            float3  _S1065 = _s_dOut_10 / _S1052;
            float3  _S1066 = _S1053 * - _S1065;
            float3  _S1067 = _S1054 * _S1065;
            float _S1068 = _S1066.x + _S1066.y + _S1066.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1069;
            (&_S1069)->primal_0 = _S1053;
            (&_S1069)->differential_0 = _S1035;
            s_bwd_length_impl_1(&_S1069, _S1068);
            _S1052 = _S1067 + _S1069.differential_0;
        }
        else
        {
            _S1052 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1070;
        (&_S1070)->primal_0 = _S1053;
        (&_S1070)->differential_0 = _S1035;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1071;
        (&_S1071)->primal_0 = _S1053;
        (&_S1071)->differential_0 = _S1035;
        s_bwd_prop_dot_0(&_S1070, &_S1071, 0.0f);
        float3  _S1072 = _S1071.differential_0 + _S1070.differential_0 + _S1052;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1073;
        (&_S1073)->primal_0 = _S1055;
        (&_S1073)->differential_0 = _S1035;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1074;
        (&_S1074)->primal_0 = _S1056;
        (&_S1074)->differential_0 = _S1035;
        s_bwd_prop_cross_0(&_S1073, &_S1074, _S1072);
        float3  s_diff_dy_T_1 = - _S1074.differential_0;
        float3  _S1075 = - s_diff_dy_T_1;
        float3  _S1076 = - _S1073.differential_0;
        FixedArray<float3 , 4>  _S1077;
        _S1077[int(0)] = _S1035;
        _S1077[int(1)] = _S1035;
        _S1077[int(2)] = _S1035;
        _S1077[int(3)] = _S1035;
        _S1077[int(2)] = _S1075;
        _S1077[int(3)] = s_diff_dy_T_1;
        _S1077[int(0)] = _S1076;
        _S1077[int(1)] = _S1073.differential_0;
        points_4[int(0)] = _S1077[int(0)];
        points_4[int(1)] = _S1077[int(1)];
        points_4[int(2)] = _S1077[int(2)];
        points_4[int(3)] = _S1077[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1035;
        points_4[int(1)] = _S1035;
        points_4[int(2)] = _S1035;
        points_4[int(3)] = _S1035;
    }
    float4  _S1078;
    if(_S1036)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1079 = points_4;
                FixedArray<float3 , 4>  _S1080 = points_4;
                FixedArray<float3 , 4>  _S1081 = points_4;
                FixedArray<float3 , 4>  _S1082 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1083 = raydir_27 * _S1082[int(3)];
                    float _S1084 = _S1083.x + _S1083.y + _S1083.z;
                    float4  _S1085 = _S1064;
                    *&((&_S1085)->w) = _S1084;
                    points_4[int(0)] = _S1079[int(0)];
                    points_4[int(1)] = _S1080[int(1)];
                    points_4[int(2)] = _S1081[int(2)];
                    points_4[int(3)] = _S1035;
                    _S1078 = _S1085;
                }
                else
                {
                    points_4[int(0)] = _S1079[int(0)];
                    points_4[int(1)] = _S1080[int(1)];
                    points_4[int(2)] = _S1081[int(2)];
                    points_4[int(3)] = _S1082[int(3)];
                    _S1078 = _S1064;
                }
                float3  _S1086 = raydir_28 * points_4[int(2)];
                float _S1087 = _S1086.x + _S1086.y + _S1086.z;
                FixedArray<float3 , 4>  _S1088 = points_4;
                FixedArray<float3 , 4>  _S1089 = points_4;
                float4  _S1090 = _S1064;
                *&((&_S1090)->z) = _S1087;
                float4  _S1091 = _S1078 + _S1090;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1088[int(1)];
                points_4[int(2)] = _S1035;
                points_4[int(3)] = _S1089[int(3)];
                _S1078 = _S1091;
            }
            else
            {
                FixedArray<float3 , 4>  _S1092 = points_4;
                FixedArray<float3 , 4>  _S1093 = points_4;
                FixedArray<float3 , 4>  _S1094 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1092[int(1)];
                points_4[int(2)] = _S1093[int(2)];
                points_4[int(3)] = _S1094[int(3)];
                _S1078 = _S1064;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1095 = points_4;
            FixedArray<float3 , 4>  _S1096 = points_4;
            FixedArray<float3 , 4>  _S1097 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1095[int(1)];
            points_4[int(2)] = _S1096[int(2)];
            points_4[int(3)] = _S1097[int(3)];
            _S1078 = _S1064;
        }
        if(_S1038)
        {
            FixedArray<float3 , 4>  _S1098 = points_4;
            float3  _S1099 = raydir_29 * points_4[int(1)];
            float _S1100 = _S1099.x + _S1099.y + _S1099.z;
            float4  _S1101 = _S1064;
            *&((&_S1101)->y) = _S1100;
            float4  _S1102 = _S1078 + _S1101;
            points_4[int(0)] = _S1035;
            points_4[int(1)] = _S1035;
            points_4[int(2)] = _S1035;
            points_4[int(3)] = _S1035;
            raydir_27 = _S1098[int(0)];
            _S1078 = _S1102;
        }
        else
        {
            FixedArray<float3 , 4>  _S1103 = points_4;
            FixedArray<float3 , 4>  _S1104 = points_4;
            FixedArray<float3 , 4>  _S1105 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1103[int(1)];
            points_4[int(2)] = _S1104[int(2)];
            points_4[int(3)] = _S1105[int(3)];
            raydir_27 = _S1035;
        }
        float3  _S1106 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1107 = _S1106.x + _S1106.y + _S1106.z;
        float4  _S1108 = _S1064;
        *&((&_S1108)->x) = _S1107;
        _S1078 = _S1078 + _S1108;
    }
    else
    {
        _S1078 = _S1064;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1078;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1109, float4  _S1110, FixedArray<float, 10>  * _S1111, bool _S1112, bool _S1113, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1114, float3  _S1115)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1116;
    float3  _S1117 = s_primal_ctx_depth_to_normal_0(_S1109, _S1110, _S1111, _S1112, _S1113, (*_S1114).primal_0, &_S1116);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1118 = _S1116;
    s_bwd_prop_depth_to_normal_0(_S1109, _S1110, _S1111, _S1112, _S1113, _S1114, _S1115, &_S1118);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1119 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1119;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, dist_coeffs_11, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9)
{
    float2  _S1120 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_16 = _S1120;
    bool _S1121 = undistort_point_0(_S1120, dist_coeffs_12, int(12), &uv_16);
    if(!_S1121)
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
    float _S1122 = (F32_exp2(((*dpx_18).primal_0))) * 0.693145751953125f * dOut_21;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1122;
    return;
}

inline __device__ float3  apply_ppisp(float3  rgb_in_0, float2  pix_coord_0, float2  image_center_0, float2  img_size_0, FixedArray<float, 36>  * params_0)
{
    PPISPParams_0 p_0;
    (&p_0)->exposure_1 = (*params_0)[int(0)];
    (&(&p_0)->vignette_params_1[int(0)])->cx_0 = (*params_0)[int(1)];
    (&(&p_0)->vignette_params_1[int(0)])->cy_0 = (*params_0)[int(2)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha0_0 = (*params_0)[int(3)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha1_0 = (*params_0)[int(4)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha2_0 = (*params_0)[int(5)];
    (&(&p_0)->vignette_params_1[int(1)])->cx_0 = (*params_0)[int(6)];
    (&(&p_0)->vignette_params_1[int(1)])->cy_0 = (*params_0)[int(7)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha0_0 = (*params_0)[int(8)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha1_0 = (*params_0)[int(9)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha2_0 = (*params_0)[int(10)];
    (&(&p_0)->vignette_params_1[int(2)])->cx_0 = (*params_0)[int(11)];
    (&(&p_0)->vignette_params_1[int(2)])->cy_0 = (*params_0)[int(12)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha0_0 = (*params_0)[int(13)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha1_0 = (*params_0)[int(14)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha2_0 = (*params_0)[int(15)];
    *&((&(&(&p_0)->color_params_1)->b_0)->x) = (*params_0)[int(16)];
    *&((&(&(&p_0)->color_params_1)->b_0)->y) = (*params_0)[int(17)];
    *&((&(&(&p_0)->color_params_1)->r_0)->x) = (*params_0)[int(18)];
    *&((&(&(&p_0)->color_params_1)->r_0)->y) = (*params_0)[int(19)];
    *&((&(&(&p_0)->color_params_1)->g_0)->x) = (*params_0)[int(20)];
    *&((&(&(&p_0)->color_params_1)->g_0)->y) = (*params_0)[int(21)];
    *&((&(&(&p_0)->color_params_1)->n_0)->x) = (*params_0)[int(22)];
    *&((&(&(&p_0)->color_params_1)->n_0)->y) = (*params_0)[int(23)];
    (&(&p_0)->crf_params_1[int(0)])->toe_0 = (*params_0)[int(24)];
    (&(&p_0)->crf_params_1[int(0)])->shoulder_0 = (*params_0)[int(25)];
    (&(&p_0)->crf_params_1[int(0)])->gamma_0 = (*params_0)[int(26)];
    (&(&p_0)->crf_params_1[int(0)])->center_0 = (*params_0)[int(27)];
    (&(&p_0)->crf_params_1[int(1)])->toe_0 = (*params_0)[int(28)];
    (&(&p_0)->crf_params_1[int(1)])->shoulder_0 = (*params_0)[int(29)];
    (&(&p_0)->crf_params_1[int(1)])->gamma_0 = (*params_0)[int(30)];
    (&(&p_0)->crf_params_1[int(1)])->center_0 = (*params_0)[int(31)];
    (&(&p_0)->crf_params_1[int(2)])->toe_0 = (*params_0)[int(32)];
    (&(&p_0)->crf_params_1[int(2)])->shoulder_0 = (*params_0)[int(33)];
    (&(&p_0)->crf_params_1[int(2)])->gamma_0 = (*params_0)[int(34)];
    (&(&p_0)->crf_params_1[int(2)])->center_0 = (*params_0)[int(35)];
    PPISPParams_0 _S1123 = p_0;
    float max_res_0 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1124 = (pix_coord_0.x - image_center_0.x) / max_res_0;
    float _S1125 = (pix_coord_0.y - image_center_0.y) / max_res_0;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_2 = _S1124 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S1125 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_4 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_4 * r2_4;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_4) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1124 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S1125 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_5 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_5 * r2_5;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_5) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1124 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S1125 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_6 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_6 * r2_6;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_6) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S1126 = rgb_out_0;
    float2  bd_0 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S1127 = 0.3333333432674408f + nd_0.x;
    float _S1128 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1128, 1.0f, 0.0f, - _S1127, - _S1128, _S1127, 0.0f), T_0);
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
    float _S1129 = _S1126.x;
    float _S1130 = _S1126.y;
    float intensity_0 = _S1129 + _S1130 + _S1126.z;
    float3  rgi_out_0 = mul_0(H_1, make_float3 (_S1129, _S1130, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1131 = rgi_out_1.x;
    float _S1132 = rgi_out_1.y;
    float3  _S1133 = clamp_1(make_float3 (_S1131, _S1132, rgi_out_1.z - _S1131 - _S1132), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1134 = _S1133.x;
    float _S1135 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(0)].toe_0))))));
    float _S1136 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(0)].shoulder_0))))));
    float _S1137 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(0)].gamma_0))))));
    float _S1138 = 1.0f / (1.0f + (F32_exp((- _S1123.crf_params_1[int(0)].center_0))));
    float a_1 = _S1136 * _S1138 / lerp_0(_S1135, _S1136, _S1138);
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1134 <= _S1138)
    {
        y_10 = a_1 * (F32_pow((_S1134 / _S1138), (_S1135)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1134) / (1.0f - _S1138)), (_S1136)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1137)));
    float _S1139 = _S1133.y;
    float _S1140 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(1)].toe_0))))));
    float _S1141 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(1)].shoulder_0))))));
    float _S1142 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(1)].gamma_0))))));
    float _S1143 = 1.0f / (1.0f + (F32_exp((- _S1123.crf_params_1[int(1)].center_0))));
    float a_2 = _S1141 * _S1143 / lerp_0(_S1140, _S1141, _S1143);
    float b_3 = 1.0f - a_2;
    if(_S1139 <= _S1143)
    {
        y_10 = a_2 * (F32_pow((_S1139 / _S1143), (_S1140)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1139) / (1.0f - _S1143)), (_S1141)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1142)));
    float _S1144 = _S1133.z;
    float _S1145 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(2)].toe_0))))));
    float _S1146 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(2)].shoulder_0))))));
    float _S1147 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1123.crf_params_1[int(2)].gamma_0))))));
    float _S1148 = 1.0f / (1.0f + (F32_exp((- _S1123.crf_params_1[int(2)].center_0))));
    float a_3 = _S1146 * _S1148 / lerp_0(_S1145, _S1146, _S1148);
    float b_4 = 1.0f - a_3;
    if(_S1144 <= _S1148)
    {
        y_10 = a_3 * (F32_pow((_S1144 / _S1148), (_S1145)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1144) / (1.0f - _S1148)), (_S1146)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1147)));
    return rgb_out_1;
}

inline __device__ float3  apply_ppisp_rqs(float3  rgb_in_1, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, FixedArray<float, 39>  * params_1)
{
    PPISPParamsRQS_0 p_1;
    (&p_1)->exposure_0 = (*params_1)[int(0)];
    (&(&p_1)->vignette_params_0[int(0)])->cx_0 = (*params_1)[int(1)];
    (&(&p_1)->vignette_params_0[int(0)])->cy_0 = (*params_1)[int(2)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha0_0 = (*params_1)[int(3)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha1_0 = (*params_1)[int(4)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha2_0 = (*params_1)[int(5)];
    (&(&p_1)->vignette_params_0[int(1)])->cx_0 = (*params_1)[int(6)];
    (&(&p_1)->vignette_params_0[int(1)])->cy_0 = (*params_1)[int(7)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha0_0 = (*params_1)[int(8)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha1_0 = (*params_1)[int(9)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha2_0 = (*params_1)[int(10)];
    (&(&p_1)->vignette_params_0[int(2)])->cx_0 = (*params_1)[int(11)];
    (&(&p_1)->vignette_params_0[int(2)])->cy_0 = (*params_1)[int(12)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha0_0 = (*params_1)[int(13)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha1_0 = (*params_1)[int(14)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha2_0 = (*params_1)[int(15)];
    *&((&(&(&p_1)->color_params_0)->b_0)->x) = (*params_1)[int(16)];
    *&((&(&(&p_1)->color_params_0)->b_0)->y) = (*params_1)[int(17)];
    *&((&(&(&p_1)->color_params_0)->r_0)->x) = (*params_1)[int(18)];
    *&((&(&(&p_1)->color_params_0)->r_0)->y) = (*params_1)[int(19)];
    *&((&(&(&p_1)->color_params_0)->g_0)->x) = (*params_1)[int(20)];
    *&((&(&(&p_1)->color_params_0)->g_0)->y) = (*params_1)[int(21)];
    *&((&(&(&p_1)->color_params_0)->n_0)->x) = (*params_1)[int(22)];
    *&((&(&(&p_1)->color_params_0)->n_0)->y) = (*params_1)[int(23)];
    (&(&p_1)->crf_params_0[int(0)])->g0_0 = (*params_1)[int(24)];
    (&(&p_1)->crf_params_0[int(0)])->g1_0 = (*params_1)[int(25)];
    (&(&p_1)->crf_params_0[int(0)])->x0_0 = (*params_1)[int(26)];
    (&(&p_1)->crf_params_0[int(0)])->y0_0 = (*params_1)[int(27)];
    (&(&p_1)->crf_params_0[int(0)])->gc_0 = (*params_1)[int(28)];
    (&(&p_1)->crf_params_0[int(1)])->g0_0 = (*params_1)[int(29)];
    (&(&p_1)->crf_params_0[int(1)])->g1_0 = (*params_1)[int(30)];
    (&(&p_1)->crf_params_0[int(1)])->x0_0 = (*params_1)[int(31)];
    (&(&p_1)->crf_params_0[int(1)])->y0_0 = (*params_1)[int(32)];
    (&(&p_1)->crf_params_0[int(1)])->gc_0 = (*params_1)[int(33)];
    (&(&p_1)->crf_params_0[int(2)])->g0_0 = (*params_1)[int(34)];
    (&(&p_1)->crf_params_0[int(2)])->g1_0 = (*params_1)[int(35)];
    (&(&p_1)->crf_params_0[int(2)])->x0_0 = (*params_1)[int(36)];
    (&(&p_1)->crf_params_0[int(2)])->y0_0 = (*params_1)[int(37)];
    (&(&p_1)->crf_params_0[int(2)])->gc_0 = (*params_1)[int(38)];
    PPISPParamsRQS_0 _S1149 = p_1;
    float max_res_1 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S1150 = (pix_coord_1.x - image_center_1.x) / max_res_1;
    float _S1151 = (pix_coord_1.y - image_center_1.y) / max_res_1;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_5 = _S1150 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S1151 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_8 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_8 * r2_8;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_8) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_8 + 1.0f, 0.0f, 1.0f);
    float dx_6 = _S1150 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S1151 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_9 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_9 * r2_9;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_9) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S1150 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S1151 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_10 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_10 * r2_10;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_10) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_10 + 1.0f, 0.0f, 1.0f);
    float3  _S1152 = rgb_out_2;
    float2  bd_1 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S1153 = 0.3333333432674408f + nd_1.x;
    float _S1154 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_4 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1154, 1.0f, 0.0f, - _S1153, - _S1154, _S1153, 0.0f), T_1);
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
    float _S1155 = _S1152.x;
    float _S1156 = _S1152.y;
    float intensity_1 = _S1155 + _S1156 + _S1152.z;
    float3  rgi_out_2 = mul_0(H_3, make_float3 (_S1155, _S1156, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / (rgi_out_2.z + 0.00000999999974738f));
    float _S1157 = rgi_out_3.x;
    float _S1158 = rgi_out_3.y;
    float3  _S1159 = clamp_1(make_float3 (_S1157, _S1158, rgi_out_3.z - _S1157 - _S1158), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S1160 = _S1159.x;
    float g0_1 = (F32_exp((_S1149.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S1149.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S1149.crf_params_0[int(0)].gc_0)));
    float y_11;
    if(_S1160 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S1160 / x0_1;
        float _S1161 = 1.0f - t0_0;
        y_11 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S1161) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S1161);
    }
    else
    {
        float _S1162 = 1.0f - y0_1;
        float _S1163 = 1.0f - x0_1;
        float s1_4 = _S1162 / _S1163;
        float t1_0 = (_S1160 - x0_1) / _S1163;
        float _S1164 = 1.0f - t1_0;
        y_11 = y0_1 + _S1162 * (s1_4 * t1_0 * t1_0 + gc_1 * t1_0 * _S1164) / (s1_4 + (gc_1 + g1_1 - 2.0f * s1_4) * t1_0 * _S1164);
    }
    *&((&rgb_out_3)->x) = y_11;
    float _S1165 = _S1159.y;
    float g0_2 = (F32_exp((_S1149.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S1149.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S1149.crf_params_0[int(1)].gc_0)));
    if(_S1165 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S1165 / x0_2;
        float _S1166 = 1.0f - t0_1;
        y_11 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S1166) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S1166);
    }
    else
    {
        float _S1167 = 1.0f - y0_2;
        float _S1168 = 1.0f - x0_2;
        float s1_5 = _S1167 / _S1168;
        float t1_1 = (_S1165 - x0_2) / _S1168;
        float _S1169 = 1.0f - t1_1;
        y_11 = y0_2 + _S1167 * (s1_5 * t1_1 * t1_1 + gc_2 * t1_1 * _S1169) / (s1_5 + (gc_2 + g1_2 - 2.0f * s1_5) * t1_1 * _S1169);
    }
    *&((&rgb_out_3)->y) = y_11;
    float _S1170 = _S1159.z;
    float g0_3 = (F32_exp((_S1149.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S1149.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S1149.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S1149.crf_params_0[int(2)].gc_0)));
    if(_S1170 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S1170 / x0_3;
        float _S1171 = 1.0f - t0_2;
        y_11 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S1171) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S1171);
    }
    else
    {
        float _S1172 = 1.0f - y0_3;
        float _S1173 = 1.0f - x0_3;
        float s1_6 = _S1172 / _S1173;
        float t1_2 = (_S1170 - x0_3) / _S1173;
        float _S1174 = 1.0f - t1_2;
        y_11 = y0_3 + _S1172 * (s1_6 * t1_2 * t1_2 + gc_3 * t1_2 * _S1174) / (s1_6 + (gc_3 + g1_3 - 2.0f * s1_6) * t1_2 * _S1174);
    }
    *&((&rgb_out_3)->z) = y_11;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1175)
{
    return (F32_exp2((_S1175)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1176, float2  _S1177)
{
    return mul_1(_S1176, _S1177);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1178, Matrix<float, 3, 3>  _S1179)
{
    return mul_3(_S1178, _S1179);
}

inline __device__ float s_primal_ctx_abs_0(float _S1180)
{
    return (F32_abs((_S1180)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1181, float3  _S1182)
{
    return mul_0(_S1181, _S1182);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1183, float3  _S1184, float3  _S1185)
{
    return clamp_1(_S1183, _S1184, _S1185);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1186, float _S1187, float _S1188)
{
    return lerp_0(_S1186, _S1187, _S1188);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1189, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1190, float3  _S1191)
{
    _d_mul_0(_S1189, _S1190, _S1191);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1192, float _S1193)
{
    _d_abs_0(_S1192, _S1193);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1194, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1195, Matrix<float, 3, 3>  _S1196)
{
    mul_2(_S1194, _S1195, _S1196);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1197, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1198, float2  _S1199)
{
    _d_mul_1(_S1197, _S1198, _S1199);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1200, float _S1201)
{
    _d_exp2_0(_S1200, _S1201);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1202 = *dprgb_in_0;
    float3  _S1203 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1204 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1205 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1206 = {
        _S1205, _S1205, _S1205
    };
    float2  _S1207 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1208 = { _S1207, _S1207, _S1207, _S1207 };
    CRFPPISPChannelParams_0 _S1209 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1210 = {
        _S1209, _S1209, _S1209
    };
    PPISPParams_0 _S1211;
    (&_S1211)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S1211)->vignette_params_1 = _S1206;
    (&_S1211)->color_params_1 = _S1208;
    (&_S1211)->crf_params_1 = _S1210;
    (&(&_S1211)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1211)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1212 = dpparams_0->primal_0[int(3)];
    (&(&_S1211)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1213 = dpparams_0->primal_0[int(4)];
    (&(&_S1211)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1214 = dpparams_0->primal_0[int(5)];
    (&(&_S1211)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1211)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1211)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1215 = dpparams_0->primal_0[int(8)];
    (&(&_S1211)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1216 = dpparams_0->primal_0[int(9)];
    (&(&_S1211)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1217 = dpparams_0->primal_0[int(10)];
    (&(&_S1211)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1211)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1211)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1218 = dpparams_0->primal_0[int(13)];
    (&(&_S1211)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1219 = dpparams_0->primal_0[int(14)];
    (&(&_S1211)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1220 = dpparams_0->primal_0[int(15)];
    (&(&_S1211)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1211)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1211)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1211)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1211)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1211)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1211)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1211)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1211)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1221 = dpparams_0->primal_0[int(24)];
    (&(&_S1211)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1222 = dpparams_0->primal_0[int(25)];
    (&(&_S1211)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1223 = dpparams_0->primal_0[int(26)];
    (&(&_S1211)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1224 = dpparams_0->primal_0[int(27)];
    (&(&_S1211)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1225 = dpparams_0->primal_0[int(28)];
    (&(&_S1211)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1226 = dpparams_0->primal_0[int(29)];
    (&(&_S1211)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1227 = dpparams_0->primal_0[int(30)];
    (&(&_S1211)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1228 = dpparams_0->primal_0[int(31)];
    (&(&_S1211)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1229 = dpparams_0->primal_0[int(32)];
    (&(&_S1211)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1230 = dpparams_0->primal_0[int(33)];
    (&(&_S1211)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1231 = dpparams_0->primal_0[int(34)];
    (&(&_S1211)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1232 = dpparams_0->primal_0[int(35)];
    (&(&_S1211)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1233 = _S1211;
    float _S1234 = s_primal_ctx_exp2_0(_S1211.exposure_1);
    float3  _S1235 = make_float3 (_S1234);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S1234);
    float _S1236 = s_primal_ctx_max_0(img_size_2.x, img_size_2.y);
    float _S1237 = (pix_coord_2.x - image_center_2.x) / _S1236;
    float _S1238 = (pix_coord_2.y - image_center_2.y) / _S1236;
    float dx_8 = _S1237 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S1238 - dpparams_0->primal_0[int(2)];
    float r2_12 = dx_8 * dx_8 + dy_6 * dy_6;
    float r4_6 = r2_12 * r2_12;
    float r6_0 = r4_6 * r2_12;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_12 + 1.0f;
    float _S1239 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1240 = rgb_out_4.x * _S1239;
    float3  _S1241 = rgb_out_4;
    *&((&_S1241)->x) = _S1240;
    float dx_9 = _S1237 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S1238 - dpparams_0->primal_0[int(7)];
    float r2_13 = dx_9 * dx_9 + dy_7 * dy_7;
    float r4_7 = r2_13 * r2_13;
    float r6_1 = r4_7 * r2_13;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_13 + 1.0f;
    float _S1242 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1241)->y) = rgb_out_4.y * _S1242;
    float dx_10 = _S1237 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S1238 - dpparams_0->primal_0[int(12)];
    float r2_14 = dx_10 * dx_10 + dy_8 * dy_8;
    float r4_8 = r2_14 * r2_14;
    float r6_2 = r4_8 * r2_14;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_14 + 1.0f;
    float _S1243 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1241)->z) = rgb_out_4.z * _S1243;
    PPISPParams_0 _S1244 = _S1211;
    float2  _S1245 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1211.color_params_1.b_0);
    float2  _S1246 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1211.color_params_1.r_0);
    float2  _S1247 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1211.color_params_1.g_0);
    float2  _S1248 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1211.color_params_1.n_0);
    float _S1249 = 0.3333333432674408f + _S1248.x;
    float _S1250 = 0.3333333432674408f + _S1248.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S1245.x, 1.0f + _S1246.x, _S1247.x, _S1245.y, _S1246.y, 1.0f + _S1247.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1250, 1.0f, 0.0f, - _S1249, - _S1250, _S1249, 0.0f);
    Matrix<float, 3, 3>  _S1251 = s_primal_ctx_mul_1(skew_0, T_2);
    float3  r0_2 = make_float3 (_S1251.rows[int(0)].x, _S1251.rows[int(0)].y, _S1251.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S1251.rows[int(1)].x, _S1251.rows[int(1)].y, _S1251.rows[int(1)].z);
    float3  r2_15 = make_float3 (_S1251.rows[int(2)].x, _S1251.rows[int(2)].y, _S1251.rows[int(2)].z);
    float3  _S1252 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S1253 = (s_primal_ctx_dot_0(_S1252, _S1252)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S1254;
    bool _S1255;
    if(_S1253)
    {
        float3  _S1256 = s_primal_ctx_cross_0(r0_2, r2_15);
        bool _S1257 = (s_primal_ctx_dot_0(_S1256, _S1256)) < 9.99999968265522539e-21f;
        if(_S1257)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_15);
        }
        else
        {
            lambda_v_6 = _S1256;
        }
        _S1255 = _S1257;
        _S1254 = _S1256;
    }
    else
    {
        lambda_v_6 = _S1252;
        _S1255 = false;
        _S1254 = _S1203;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S1258 = s_primal_ctx_mul_1(T_2, D_0);
    Matrix<float, 3, 3>  _S1259 = s_primal_ctx_mul_1(_S1258, S_inv_0);
    bool _S1260 = (s_primal_ctx_abs_0(_S1259.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S1261;
    float _S1262;
    if(_S1260)
    {
        float inv_s_0 = 1.0f / _S1259.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1263 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1264 = _S1259.rows[int(2)].z * _S1259.rows[int(2)].z;
        H_4 = _S1259 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1261 = _S1263;
        _S1262 = _S1264;
    }
    else
    {
        H_4 = _S1259;
        _S1261 = _S1204;
        _S1262 = 0.0f;
    }
    float _S1265 = _S1241.x;
    float _S1266 = _S1241.y;
    float intensity_2 = _S1265 + _S1266 + _S1241.z;
    float3  rgi_in_0 = make_float3 (_S1265, _S1266, intensity_2);
    float3  _S1267 = s_primal_ctx_mul_2(H_4, rgi_in_0);
    float _S1268 = _S1267.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S1268;
    float3  _S1269 = make_float3 (norm_factor_0);
    float _S1270 = _S1268 * _S1268;
    float3  rgi_out_4 = _S1267 * make_float3 (norm_factor_0);
    float _S1271 = rgi_out_4.x;
    float _S1272 = rgi_out_4.y;
    float3  _S1273 = make_float3 (_S1271, _S1272, rgi_out_4.z - _S1271 - _S1272);
    float3  _S1274 = make_float3 (0.0f);
    float3  _S1275 = make_float3 (1.0f);
    float3  _S1276 = s_primal_ctx_clamp_1(_S1273, _S1274, _S1275);
    float _S1277 = _S1276.x;
    float _S1278 = 1.0f + s_primal_ctx_exp_0(_S1221);
    float _S1279 = 0.30000001192092896f + s_primal_ctx_log_0(_S1278);
    float _S1280 = 1.0f + s_primal_ctx_exp_0(_S1222);
    float _S1281 = 0.30000001192092896f + s_primal_ctx_log_0(_S1280);
    float _S1282 = 1.0f + s_primal_ctx_exp_0(_S1223);
    float _S1283 = 0.10000000149011612f + s_primal_ctx_log_0(_S1282);
    float _S1284 = - _S1224;
    float _S1285 = 1.0f + s_primal_ctx_exp_0(_S1284);
    float _S1286 = 1.0f / _S1285;
    float _S1287 = _S1285 * _S1285;
    float _S1288 = s_primal_ctx_lerp_0(_S1279, _S1281, _S1286);
    float _S1289 = _S1281 * _S1286;
    float a_4 = _S1289 / _S1288;
    float _S1290 = _S1288 * _S1288;
    float b_5 = 1.0f - a_4;
    bool _S1291 = _S1277 <= _S1286;
    float y_12;
    float _S1292;
    float _S1293;
    float _S1294;
    float _S1295;
    float _S1296;
    float _S1297;
    float _S1298;
    float _S1299;
    if(_S1291)
    {
        float _S1300 = _S1277 / _S1286;
        float _S1301 = _S1286 * _S1286;
        float _S1302 = s_primal_ctx_pow_0(_S1300, _S1279);
        y_12 = a_4 * _S1302;
        _S1292 = 0.0f;
        _S1293 = 0.0f;
        _S1294 = 0.0f;
        _S1295 = 0.0f;
        _S1296 = 0.0f;
        _S1297 = _S1302;
        _S1298 = _S1300;
        _S1299 = _S1301;
    }
    else
    {
        float _S1303 = 1.0f - _S1277;
        float _S1304 = 1.0f - _S1286;
        float _S1305 = _S1303 / _S1304;
        float _S1306 = _S1304 * _S1304;
        float _S1307 = s_primal_ctx_pow_0(_S1305, _S1281);
        y_12 = 1.0f - b_5 * _S1307;
        _S1292 = _S1307;
        _S1293 = _S1305;
        _S1294 = _S1306;
        _S1295 = _S1303;
        _S1296 = _S1304;
        _S1297 = 0.0f;
        _S1298 = 0.0f;
        _S1299 = 0.0f;
    }
    float _S1308 = s_primal_ctx_max_0(0.0f, y_12);
    float _S1309 = _S1276.y;
    float _S1310 = 1.0f + s_primal_ctx_exp_0(_S1225);
    float _S1311 = 0.30000001192092896f + s_primal_ctx_log_0(_S1310);
    float _S1312 = 1.0f + s_primal_ctx_exp_0(_S1226);
    float _S1313 = 0.30000001192092896f + s_primal_ctx_log_0(_S1312);
    float _S1314 = 1.0f + s_primal_ctx_exp_0(_S1227);
    float _S1315 = 0.10000000149011612f + s_primal_ctx_log_0(_S1314);
    float _S1316 = - _S1228;
    float _S1317 = 1.0f + s_primal_ctx_exp_0(_S1316);
    float _S1318 = 1.0f / _S1317;
    float _S1319 = _S1317 * _S1317;
    float _S1320 = s_primal_ctx_lerp_0(_S1311, _S1313, _S1318);
    float _S1321 = _S1313 * _S1318;
    float a_5 = _S1321 / _S1320;
    float _S1322 = _S1320 * _S1320;
    float b_6 = 1.0f - a_5;
    bool _S1323 = _S1309 <= _S1318;
    float y_13;
    float _S1324;
    float _S1325;
    float _S1326;
    float _S1327;
    float _S1328;
    float _S1329;
    float _S1330;
    float _S1331;
    if(_S1323)
    {
        float _S1332 = _S1309 / _S1318;
        float _S1333 = _S1318 * _S1318;
        float _S1334 = s_primal_ctx_pow_0(_S1332, _S1311);
        y_13 = a_5 * _S1334;
        _S1324 = 0.0f;
        _S1325 = 0.0f;
        _S1326 = 0.0f;
        _S1327 = 0.0f;
        _S1328 = 0.0f;
        _S1329 = _S1334;
        _S1330 = _S1332;
        _S1331 = _S1333;
    }
    else
    {
        float _S1335 = 1.0f - _S1309;
        float _S1336 = 1.0f - _S1318;
        float _S1337 = _S1335 / _S1336;
        float _S1338 = _S1336 * _S1336;
        float _S1339 = s_primal_ctx_pow_0(_S1337, _S1313);
        y_13 = 1.0f - b_6 * _S1339;
        _S1324 = _S1339;
        _S1325 = _S1337;
        _S1326 = _S1338;
        _S1327 = _S1335;
        _S1328 = _S1336;
        _S1329 = 0.0f;
        _S1330 = 0.0f;
        _S1331 = 0.0f;
    }
    float _S1340 = s_primal_ctx_max_0(0.0f, y_13);
    float _S1341 = _S1276.z;
    float _S1342 = 1.0f + s_primal_ctx_exp_0(_S1229);
    float _S1343 = 0.30000001192092896f + s_primal_ctx_log_0(_S1342);
    float _S1344 = 1.0f + s_primal_ctx_exp_0(_S1230);
    float _S1345 = 0.30000001192092896f + s_primal_ctx_log_0(_S1344);
    float _S1346 = 1.0f + s_primal_ctx_exp_0(_S1231);
    float _S1347 = 0.10000000149011612f + s_primal_ctx_log_0(_S1346);
    float _S1348 = - _S1232;
    float _S1349 = 1.0f + s_primal_ctx_exp_0(_S1348);
    float _S1350 = 1.0f / _S1349;
    float _S1351 = _S1349 * _S1349;
    float _S1352 = s_primal_ctx_lerp_0(_S1343, _S1345, _S1350);
    float _S1353 = _S1345 * _S1350;
    float a_6 = _S1353 / _S1352;
    float _S1354 = _S1352 * _S1352;
    float b_7 = 1.0f - a_6;
    bool _S1355 = _S1341 <= _S1350;
    float y_14;
    float _S1356;
    float _S1357;
    float _S1358;
    float _S1359;
    float _S1360;
    float _S1361;
    float _S1362;
    float _S1363;
    if(_S1355)
    {
        float _S1364 = _S1341 / _S1350;
        float _S1365 = _S1350 * _S1350;
        float _S1366 = s_primal_ctx_pow_0(_S1364, _S1343);
        y_14 = a_6 * _S1366;
        _S1356 = 0.0f;
        _S1357 = 0.0f;
        _S1358 = 0.0f;
        _S1359 = 0.0f;
        _S1360 = 0.0f;
        _S1361 = _S1366;
        _S1362 = _S1364;
        _S1363 = _S1365;
    }
    else
    {
        float _S1367 = 1.0f - _S1341;
        float _S1368 = 1.0f - _S1350;
        float _S1369 = _S1367 / _S1368;
        float _S1370 = _S1368 * _S1368;
        float _S1371 = s_primal_ctx_pow_0(_S1369, _S1345);
        y_14 = 1.0f - b_7 * _S1371;
        _S1356 = _S1371;
        _S1357 = _S1369;
        _S1358 = _S1370;
        _S1359 = _S1367;
        _S1360 = _S1368;
        _S1361 = 0.0f;
        _S1362 = 0.0f;
        _S1363 = 0.0f;
    }
    float _S1372 = s_primal_ctx_max_0(0.0f, y_14);
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S1372;
    (&_S1373)->differential_0 = 0.0f;
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = _S1347;
    (&_S1374)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1373, &_S1374, _s_dOut_11.z);
    DiffPair_float_0 _S1375 = _S1374;
    DiffPair_float_0 _S1376;
    (&_S1376)->primal_0 = 0.0f;
    (&_S1376)->differential_0 = 0.0f;
    DiffPair_float_0 _S1377;
    (&_S1377)->primal_0 = y_14;
    (&_S1377)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1376, &_S1377, _S1373.differential_0);
    DiffPair_float_0 _S1378 = _S1377;
    if(_S1355)
    {
        float _S1379 = a_6 * _S1378.differential_0;
        float _S1380 = _S1361 * _S1378.differential_0;
        DiffPair_float_0 _S1381;
        (&_S1381)->primal_0 = _S1362;
        (&_S1381)->differential_0 = 0.0f;
        DiffPair_float_0 _S1382;
        (&_S1382)->primal_0 = _S1343;
        (&_S1382)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1381, &_S1382, _S1379);
        float _S1383 = _S1381.differential_0 / _S1363;
        float _S1384 = _S1341 * - _S1383;
        float _S1385 = _S1350 * _S1383;
        y_14 = 0.0f;
        _S1356 = _S1380;
        _S1357 = _S1384;
        _S1358 = 0.0f;
        _S1359 = _S1382.differential_0;
        _S1360 = _S1385;
    }
    else
    {
        float _S1386 = - _S1378.differential_0;
        float _S1387 = b_7 * _S1386;
        float _S1388 = _S1356 * _S1386;
        DiffPair_float_0 _S1389;
        (&_S1389)->primal_0 = _S1357;
        (&_S1389)->differential_0 = 0.0f;
        DiffPair_float_0 _S1390;
        (&_S1390)->primal_0 = _S1345;
        (&_S1390)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1389, &_S1390, _S1387);
        float _S1391 = _S1389.differential_0 / _S1358;
        float _S1392 = - (_S1359 * - _S1391);
        float _S1393 = - (_S1360 * _S1391);
        y_14 = _S1388;
        _S1356 = 0.0f;
        _S1357 = _S1392;
        _S1358 = _S1390.differential_0;
        _S1359 = 0.0f;
        _S1360 = _S1393;
    }
    float _S1394 = (- y_14 + _S1356) / _S1354;
    float _S1395 = _S1353 * - _S1394;
    float _S1396 = _S1352 * _S1394;
    float _S1397 = _S1345 * _S1396;
    float _S1398 = _S1350 * _S1396;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1343;
    (&_S1399)->differential_0 = 0.0f;
    DiffPair_float_0 _S1400;
    (&_S1400)->primal_0 = _S1345;
    (&_S1400)->differential_0 = 0.0f;
    DiffPair_float_0 _S1401;
    (&_S1401)->primal_0 = _S1350;
    (&_S1401)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1399, &_S1400, &_S1401, _S1395);
    float _S1402 = - ((_S1397 + _S1401.differential_0 + _S1357) / _S1351);
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = _S1348;
    (&_S1403)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1403, _S1402);
    float _S1404 = - _S1403.differential_0;
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = _S1346;
    (&_S1405)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1405, _S1375.differential_0);
    DiffPair_float_0 _S1406;
    (&_S1406)->primal_0 = _S1231;
    (&_S1406)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1406, _S1405.differential_0);
    DiffPair_float_0 _S1407 = _S1406;
    float _S1408 = _S1398 + _S1400.differential_0 + _S1358;
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1344;
    (&_S1409)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1409, _S1408);
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = _S1230;
    (&_S1410)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1410, _S1409.differential_0);
    DiffPair_float_0 _S1411 = _S1410;
    float _S1412 = _S1399.differential_0 + _S1359;
    DiffPair_float_0 _S1413;
    (&_S1413)->primal_0 = _S1342;
    (&_S1413)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1413, _S1412);
    DiffPair_float_0 _S1414;
    (&_S1414)->primal_0 = _S1229;
    (&_S1414)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1414, _S1413.differential_0);
    DiffPair_float_0 _S1415 = _S1414;
    float3  _S1416 = make_float3 (0.0f, 0.0f, _S1360);
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1340;
    (&_S1417)->differential_0 = 0.0f;
    DiffPair_float_0 _S1418;
    (&_S1418)->primal_0 = _S1315;
    (&_S1418)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1417, &_S1418, _s_dOut_11.y);
    DiffPair_float_0 _S1419 = _S1418;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = 0.0f;
    (&_S1420)->differential_0 = 0.0f;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = y_13;
    (&_S1421)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1420, &_S1421, _S1417.differential_0);
    DiffPair_float_0 _S1422 = _S1421;
    if(_S1323)
    {
        float _S1423 = a_5 * _S1422.differential_0;
        float _S1424 = _S1329 * _S1422.differential_0;
        DiffPair_float_0 _S1425;
        (&_S1425)->primal_0 = _S1330;
        (&_S1425)->differential_0 = 0.0f;
        DiffPair_float_0 _S1426;
        (&_S1426)->primal_0 = _S1311;
        (&_S1426)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1425, &_S1426, _S1423);
        float _S1427 = _S1425.differential_0 / _S1331;
        float _S1428 = _S1309 * - _S1427;
        float _S1429 = _S1318 * _S1427;
        y_13 = 0.0f;
        _S1324 = _S1424;
        _S1325 = _S1428;
        _S1326 = 0.0f;
        _S1327 = _S1426.differential_0;
        _S1328 = _S1429;
    }
    else
    {
        float _S1430 = - _S1422.differential_0;
        float _S1431 = b_6 * _S1430;
        float _S1432 = _S1324 * _S1430;
        DiffPair_float_0 _S1433;
        (&_S1433)->primal_0 = _S1325;
        (&_S1433)->differential_0 = 0.0f;
        DiffPair_float_0 _S1434;
        (&_S1434)->primal_0 = _S1313;
        (&_S1434)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1433, &_S1434, _S1431);
        float _S1435 = _S1433.differential_0 / _S1326;
        float _S1436 = - (_S1327 * - _S1435);
        float _S1437 = - (_S1328 * _S1435);
        y_13 = _S1432;
        _S1324 = 0.0f;
        _S1325 = _S1436;
        _S1326 = _S1434.differential_0;
        _S1327 = 0.0f;
        _S1328 = _S1437;
    }
    float _S1438 = (- y_13 + _S1324) / _S1322;
    float _S1439 = _S1321 * - _S1438;
    float _S1440 = _S1320 * _S1438;
    float _S1441 = _S1313 * _S1440;
    float _S1442 = _S1318 * _S1440;
    DiffPair_float_0 _S1443;
    (&_S1443)->primal_0 = _S1311;
    (&_S1443)->differential_0 = 0.0f;
    DiffPair_float_0 _S1444;
    (&_S1444)->primal_0 = _S1313;
    (&_S1444)->differential_0 = 0.0f;
    DiffPair_float_0 _S1445;
    (&_S1445)->primal_0 = _S1318;
    (&_S1445)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1443, &_S1444, &_S1445, _S1439);
    float _S1446 = - ((_S1441 + _S1445.differential_0 + _S1325) / _S1319);
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1316;
    (&_S1447)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1447, _S1446);
    float _S1448 = - _S1447.differential_0;
    DiffPair_float_0 _S1449;
    (&_S1449)->primal_0 = _S1314;
    (&_S1449)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1449, _S1419.differential_0);
    DiffPair_float_0 _S1450;
    (&_S1450)->primal_0 = _S1227;
    (&_S1450)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1450, _S1449.differential_0);
    DiffPair_float_0 _S1451 = _S1450;
    float _S1452 = _S1442 + _S1444.differential_0 + _S1326;
    DiffPair_float_0 _S1453;
    (&_S1453)->primal_0 = _S1312;
    (&_S1453)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1453, _S1452);
    DiffPair_float_0 _S1454;
    (&_S1454)->primal_0 = _S1226;
    (&_S1454)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1454, _S1453.differential_0);
    DiffPair_float_0 _S1455 = _S1454;
    float _S1456 = _S1443.differential_0 + _S1327;
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1310;
    (&_S1457)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1457, _S1456);
    DiffPair_float_0 _S1458;
    (&_S1458)->primal_0 = _S1225;
    (&_S1458)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1458, _S1457.differential_0);
    DiffPair_float_0 _S1459 = _S1458;
    float3  _S1460 = _S1416 + make_float3 (0.0f, _S1328, 0.0f);
    DiffPair_float_0 _S1461;
    (&_S1461)->primal_0 = _S1308;
    (&_S1461)->differential_0 = 0.0f;
    DiffPair_float_0 _S1462;
    (&_S1462)->primal_0 = _S1283;
    (&_S1462)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1461, &_S1462, _s_dOut_11.x);
    DiffPair_float_0 _S1463 = _S1462;
    DiffPair_float_0 _S1464;
    (&_S1464)->primal_0 = 0.0f;
    (&_S1464)->differential_0 = 0.0f;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = y_12;
    (&_S1465)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1464, &_S1465, _S1461.differential_0);
    DiffPair_float_0 _S1466 = _S1465;
    if(_S1291)
    {
        float _S1467 = a_4 * _S1466.differential_0;
        float _S1468 = _S1297 * _S1466.differential_0;
        DiffPair_float_0 _S1469;
        (&_S1469)->primal_0 = _S1298;
        (&_S1469)->differential_0 = 0.0f;
        DiffPair_float_0 _S1470;
        (&_S1470)->primal_0 = _S1279;
        (&_S1470)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1469, &_S1470, _S1467);
        float _S1471 = _S1469.differential_0 / _S1299;
        float _S1472 = _S1277 * - _S1471;
        float _S1473 = _S1286 * _S1471;
        y_12 = 0.0f;
        _S1292 = _S1468;
        _S1293 = _S1472;
        _S1294 = 0.0f;
        _S1295 = _S1470.differential_0;
        _S1296 = _S1473;
    }
    else
    {
        float _S1474 = - _S1466.differential_0;
        float _S1475 = b_5 * _S1474;
        float _S1476 = _S1292 * _S1474;
        DiffPair_float_0 _S1477;
        (&_S1477)->primal_0 = _S1293;
        (&_S1477)->differential_0 = 0.0f;
        DiffPair_float_0 _S1478;
        (&_S1478)->primal_0 = _S1281;
        (&_S1478)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1477, &_S1478, _S1475);
        float _S1479 = _S1477.differential_0 / _S1294;
        float _S1480 = - (_S1295 * - _S1479);
        float _S1481 = - (_S1296 * _S1479);
        y_12 = _S1476;
        _S1292 = 0.0f;
        _S1293 = _S1480;
        _S1294 = _S1478.differential_0;
        _S1295 = 0.0f;
        _S1296 = _S1481;
    }
    float _S1482 = (- y_12 + _S1292) / _S1290;
    float _S1483 = _S1289 * - _S1482;
    float _S1484 = _S1288 * _S1482;
    float _S1485 = _S1281 * _S1484;
    float _S1486 = _S1286 * _S1484;
    DiffPair_float_0 _S1487;
    (&_S1487)->primal_0 = _S1279;
    (&_S1487)->differential_0 = 0.0f;
    DiffPair_float_0 _S1488;
    (&_S1488)->primal_0 = _S1281;
    (&_S1488)->differential_0 = 0.0f;
    DiffPair_float_0 _S1489;
    (&_S1489)->primal_0 = _S1286;
    (&_S1489)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1487, &_S1488, &_S1489, _S1483);
    float _S1490 = - ((_S1485 + _S1489.differential_0 + _S1293) / _S1287);
    DiffPair_float_0 _S1491;
    (&_S1491)->primal_0 = _S1284;
    (&_S1491)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1491, _S1490);
    float _S1492 = - _S1491.differential_0;
    DiffPair_float_0 _S1493;
    (&_S1493)->primal_0 = _S1282;
    (&_S1493)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1493, _S1463.differential_0);
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = _S1223;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1494, _S1493.differential_0);
    DiffPair_float_0 _S1495 = _S1494;
    float _S1496 = _S1486 + _S1488.differential_0 + _S1294;
    DiffPair_float_0 _S1497;
    (&_S1497)->primal_0 = _S1280;
    (&_S1497)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1497, _S1496);
    DiffPair_float_0 _S1498;
    (&_S1498)->primal_0 = _S1222;
    (&_S1498)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1498, _S1497.differential_0);
    DiffPair_float_0 _S1499 = _S1498;
    float _S1500 = _S1487.differential_0 + _S1295;
    DiffPair_float_0 _S1501;
    (&_S1501)->primal_0 = _S1278;
    (&_S1501)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1501, _S1500);
    DiffPair_float_0 _S1502;
    (&_S1502)->primal_0 = _S1221;
    (&_S1502)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1502, _S1501.differential_0);
    DiffPair_float_0 _S1503 = _S1502;
    float3  _S1504 = _S1460 + make_float3 (_S1296, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1505;
    (&_S1505)->primal_0 = _S1273;
    (&_S1505)->differential_0 = _S1203;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1506;
    (&_S1506)->primal_0 = _S1274;
    (&_S1506)->differential_0 = _S1203;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1507;
    (&_S1507)->primal_0 = _S1275;
    (&_S1507)->differential_0 = _S1203;
    s_bwd_prop_clamp_1(&_S1505, &_S1506, &_S1507, _S1504);
    float _S1508 = - _S1505.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1505.differential_0.x + _S1508, _S1505.differential_0.y + _S1508, _S1505.differential_0.z);
    float3  _S1509 = _S1267 * s_diff_rgi_out_T_0;
    float _S1510 = (_S1509.x + _S1509.y + _S1509.z) / _S1270;
    float _S1511 = _S1268 * _S1510;
    float3  _S1512 = _S1269 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_2 * - _S1510);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1513;
    (&_S1513)->primal_0 = H_4;
    (&_S1513)->differential_0 = _S1204;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1514;
    (&_S1514)->primal_0 = rgi_in_0;
    (&_S1514)->differential_0 = _S1203;
    s_bwd_prop_mul_1(&_S1513, &_S1514, _S1512);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1515 = _S1513;
    float _S1516 = _S1511 + _S1514.differential_0.z;
    float _S1517 = _S1514.differential_0.y + _S1516;
    float _S1518 = _S1514.differential_0.x + _S1516;
    float3  _S1519 = make_float3 (_S1518, _S1517, _S1516);
    if(_S1260)
    {
        Matrix<float, 3, 3>  _S1520 = _S1259 * _S1515.differential_0;
        Matrix<float, 3, 3>  _S1521 = _S1261 * _S1515.differential_0;
        _S1262 = - ((_S1520.rows[int(0)].x + _S1520.rows[int(0)].y + _S1520.rows[int(0)].z + _S1520.rows[int(1)].x + _S1520.rows[int(1)].y + _S1520.rows[int(1)].z + _S1520.rows[int(2)].x + _S1520.rows[int(2)].y + _S1520.rows[int(2)].z) / _S1262);
        H_4 = _S1521;
    }
    else
    {
        _S1262 = 0.0f;
        H_4 = _S1515.differential_0;
    }
    DiffPair_float_0 _S1522;
    (&_S1522)->primal_0 = _S1259.rows[int(2)].z;
    (&_S1522)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1522, 0.0f);
    float _S1523 = _S1522.differential_0 + _S1262;
    float3  _S1524 = _S1203;
    *&((&_S1524)->z) = _S1523;
    Matrix<float, 3, 3>  _S1525 = _S1204;
    _S1525[int(2)] = _S1524;
    Matrix<float, 3, 3>  _S1526 = H_4 + _S1525;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1527;
    (&_S1527)->primal_0 = _S1258;
    (&_S1527)->differential_0 = _S1204;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1528;
    (&_S1528)->primal_0 = S_inv_0;
    (&_S1528)->differential_0 = _S1204;
    s_bwd_prop_mul_2(&_S1527, &_S1528, _S1526);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1529;
    (&_S1529)->primal_0 = T_2;
    (&_S1529)->differential_0 = _S1204;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1530;
    (&_S1530)->primal_0 = D_0;
    (&_S1530)->differential_0 = _S1204;
    s_bwd_prop_mul_2(&_S1529, &_S1530, _S1527.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1531 = _S1529;
    float3  _S1532 = make_float3 (_S1530.differential_0.rows[int(0)].x, _S1530.differential_0.rows[int(1)].y, _S1530.differential_0.rows[int(2)].z);
    float3  _S1533;
    if(_S1253)
    {
        if(_S1255)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1534;
            (&_S1534)->primal_0 = r1_2;
            (&_S1534)->differential_0 = _S1203;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1535;
            (&_S1535)->primal_0 = r2_15;
            (&_S1535)->differential_0 = _S1203;
            s_bwd_prop_cross_0(&_S1534, &_S1535, _S1532);
            _S1241 = _S1203;
            lambda_v_6 = _S1535.differential_0;
            _S1533 = _S1534.differential_0;
        }
        else
        {
            _S1241 = _S1532;
            lambda_v_6 = _S1203;
            _S1533 = _S1203;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1536;
        (&_S1536)->primal_0 = _S1254;
        (&_S1536)->differential_0 = _S1203;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1537;
        (&_S1537)->primal_0 = _S1254;
        (&_S1537)->differential_0 = _S1203;
        s_bwd_prop_dot_0(&_S1536, &_S1537, 0.0f);
        float3  _S1538 = _S1537.differential_0 + _S1536.differential_0 + _S1241;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1539;
        (&_S1539)->primal_0 = r0_2;
        (&_S1539)->differential_0 = _S1203;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1540;
        (&_S1540)->primal_0 = r2_15;
        (&_S1540)->differential_0 = _S1203;
        s_bwd_prop_cross_0(&_S1539, &_S1540, _S1538);
        float3  _S1541 = _S1540.differential_0 + lambda_v_6;
        _S1241 = _S1203;
        lambda_v_6 = _S1541;
        _S1254 = _S1533;
        _S1533 = _S1539.differential_0;
    }
    else
    {
        _S1241 = _S1532;
        lambda_v_6 = _S1203;
        _S1254 = _S1203;
        _S1533 = _S1203;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1542;
    (&_S1542)->primal_0 = _S1252;
    (&_S1542)->differential_0 = _S1203;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1543;
    (&_S1543)->primal_0 = _S1252;
    (&_S1543)->differential_0 = _S1203;
    s_bwd_prop_dot_0(&_S1542, &_S1543, 0.0f);
    float3  _S1544 = _S1543.differential_0 + _S1542.differential_0 + _S1241;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1545;
    (&_S1545)->primal_0 = r0_2;
    (&_S1545)->differential_0 = _S1203;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1546;
    (&_S1546)->primal_0 = r1_2;
    (&_S1546)->differential_0 = _S1203;
    s_bwd_prop_cross_0(&_S1545, &_S1546, _S1544);
    float3  _S1547 = _S1203;
    *&((&_S1547)->z) = lambda_v_6.z;
    *&((&_S1547)->y) = lambda_v_6.y;
    *&((&_S1547)->x) = lambda_v_6.x;
    float3  _S1548 = _S1546.differential_0 + _S1254;
    float3  _S1549 = _S1203;
    *&((&_S1549)->z) = _S1548.z;
    *&((&_S1549)->y) = _S1548.y;
    *&((&_S1549)->x) = _S1548.x;
    float3  _S1550 = _S1545.differential_0 + _S1533;
    float3  _S1551 = _S1203;
    *&((&_S1551)->z) = _S1550.z;
    *&((&_S1551)->y) = _S1550.y;
    *&((&_S1551)->x) = _S1550.x;
    Matrix<float, 3, 3>  _S1552 = _S1204;
    _S1552[int(2)] = _S1547;
    _S1552[int(1)] = _S1549;
    _S1552[int(0)] = _S1551;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1553;
    (&_S1553)->primal_0 = skew_0;
    (&_S1553)->differential_0 = _S1204;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1554;
    (&_S1554)->primal_0 = T_2;
    (&_S1554)->differential_0 = _S1204;
    s_bwd_prop_mul_2(&_S1553, &_S1554, _S1552);
    Matrix<float, 3, 3>  _S1555 = _S1554.differential_0 + _S1531.differential_0;
    float2  _S1556 = make_float2 (_S1553.differential_0.rows[int(2)].y + - _S1553.differential_0.rows[int(1)].z, _S1553.differential_0.rows[int(0)].z + - _S1553.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1557 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1558;
    (&_S1558)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1558)->differential_0 = _S1557;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1559;
    (&_S1559)->primal_0 = _S1244.color_params_1.n_0;
    (&_S1559)->differential_0 = _S1207;
    s_bwd_prop_mul_3(&_S1558, &_S1559, _S1556);
    float2  _S1560 = make_float2 (_S1555.rows[int(0)].z, _S1555.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1561;
    (&_S1561)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1561)->differential_0 = _S1557;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1562;
    (&_S1562)->primal_0 = _S1244.color_params_1.g_0;
    (&_S1562)->differential_0 = _S1207;
    s_bwd_prop_mul_3(&_S1561, &_S1562, _S1560);
    float2  _S1563 = make_float2 (_S1555.rows[int(0)].y, _S1555.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1564;
    (&_S1564)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1564)->differential_0 = _S1557;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1565;
    (&_S1565)->primal_0 = _S1244.color_params_1.r_0;
    (&_S1565)->differential_0 = _S1207;
    s_bwd_prop_mul_3(&_S1564, &_S1565, _S1563);
    float2  _S1566 = make_float2 (_S1555.rows[int(0)].x, _S1555.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1567;
    (&_S1567)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1567)->differential_0 = _S1557;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1568;
    (&_S1568)->primal_0 = _S1244.color_params_1.b_0;
    (&_S1568)->differential_0 = _S1207;
    s_bwd_prop_mul_3(&_S1567, &_S1568, _S1566);
    ColorPPISPParams_0 _S1569 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1569)->n_0 = _S1559.differential_0;
    (&_S1569)->g_0 = _S1562.differential_0;
    (&_S1569)->r_0 = _S1565.differential_0;
    (&_S1569)->b_0 = _S1568.differential_0;
    _S1241 = _S1519;
    *&((&_S1241)->z) = 0.0f;
    float _S1570 = rgb_out_4.z * _S1516;
    float _S1571 = _S1243 * _S1516;
    DiffPair_float_0 _S1572;
    (&_S1572)->primal_0 = falloff_2;
    (&_S1572)->differential_0 = 0.0f;
    DiffPair_float_0 _S1573;
    (&_S1573)->primal_0 = 0.0f;
    (&_S1573)->differential_0 = 0.0f;
    DiffPair_float_0 _S1574;
    (&_S1574)->primal_0 = 1.0f;
    (&_S1574)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1572, &_S1573, &_S1574, _S1570);
    float _S1575 = r2_14 * _S1572.differential_0;
    float _S1576 = r4_8 * _S1572.differential_0;
    float s_diff_r6_T_0 = _S1220 * _S1572.differential_0;
    float _S1577 = r6_2 * _S1572.differential_0;
    float _S1578 = r2_14 * (_S1219 * _S1572.differential_0 + r2_14 * s_diff_r6_T_0);
    float _S1579 = _S1218 * _S1572.differential_0 + r4_8 * s_diff_r6_T_0 + _S1578 + _S1578;
    float _S1580 = dy_8 * _S1579;
    float _S1581 = dx_10 * _S1579;
    float _S1582 = - (_S1580 + _S1580);
    float _S1583 = - (_S1581 + _S1581);
    *&((&_S1241)->y) = 0.0f;
    float _S1584 = rgb_out_4.y * _S1517;
    float _S1585 = _S1242 * _S1517;
    DiffPair_float_0 _S1586;
    (&_S1586)->primal_0 = falloff_1;
    (&_S1586)->differential_0 = 0.0f;
    DiffPair_float_0 _S1587;
    (&_S1587)->primal_0 = 0.0f;
    (&_S1587)->differential_0 = 0.0f;
    DiffPair_float_0 _S1588;
    (&_S1588)->primal_0 = 1.0f;
    (&_S1588)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1586, &_S1587, &_S1588, _S1584);
    float _S1589 = r2_13 * _S1586.differential_0;
    float _S1590 = r4_7 * _S1586.differential_0;
    float s_diff_r6_T_1 = _S1217 * _S1586.differential_0;
    float _S1591 = r6_1 * _S1586.differential_0;
    float _S1592 = r2_13 * (_S1216 * _S1586.differential_0 + r2_13 * s_diff_r6_T_1);
    float _S1593 = _S1215 * _S1586.differential_0 + r4_7 * s_diff_r6_T_1 + _S1592 + _S1592;
    float _S1594 = dy_7 * _S1593;
    float _S1595 = dx_9 * _S1593;
    float _S1596 = - (_S1594 + _S1594);
    float _S1597 = - (_S1595 + _S1595);
    *&((&_S1241)->x) = 0.0f;
    float _S1598 = rgb_out_4.x * _S1518;
    float _S1599 = _S1239 * _S1518;
    DiffPair_float_0 _S1600;
    (&_S1600)->primal_0 = falloff_0;
    (&_S1600)->differential_0 = 0.0f;
    DiffPair_float_0 _S1601;
    (&_S1601)->primal_0 = 0.0f;
    (&_S1601)->differential_0 = 0.0f;
    DiffPair_float_0 _S1602;
    (&_S1602)->primal_0 = 1.0f;
    (&_S1602)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1600, &_S1601, &_S1602, _S1598);
    float _S1603 = r2_12 * _S1600.differential_0;
    float _S1604 = r4_6 * _S1600.differential_0;
    float s_diff_r6_T_2 = _S1214 * _S1600.differential_0;
    float _S1605 = r6_0 * _S1600.differential_0;
    float _S1606 = r2_12 * (_S1213 * _S1600.differential_0 + r2_12 * s_diff_r6_T_2);
    float _S1607 = _S1212 * _S1600.differential_0 + r4_6 * s_diff_r6_T_2 + _S1606 + _S1606;
    float _S1608 = dy_6 * _S1607;
    float _S1609 = dx_8 * _S1607;
    float _S1610 = - (_S1608 + _S1608);
    float _S1611 = - (_S1609 + _S1609);
    float3  _S1612 = _S1203;
    *&((&_S1612)->z) = _S1571;
    *&((&_S1612)->y) = _S1585;
    *&((&_S1612)->x) = _S1599;
    float3  _S1613 = _S1241 + _S1612;
    float3  _S1614 = _S1202.primal_0 * _S1613;
    float3  _S1615 = _S1235 * _S1613;
    float _S1616 = _S1614.x + _S1614.y + _S1614.z;
    DiffPair_float_0 _S1617;
    (&_S1617)->primal_0 = _S1233.exposure_1;
    (&_S1617)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1617, _S1616);
    PPISPParams_0 _S1618 = PPISPParams_x24_syn_dzero_0();
    (&_S1618)->color_params_1 = _S1569;
    (&_S1618)->exposure_1 = _S1617.differential_0;
    _S1211 = _S1618;
    (&(&_S1211)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1619 = _S1618.crf_params_1[int(2)].center_0 + _S1404;
    (&(&_S1211)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1620 = _S1618.crf_params_1[int(2)].gamma_0 + _S1407.differential_0;
    (&(&_S1211)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1621 = _S1618.crf_params_1[int(2)].shoulder_0 + _S1411.differential_0;
    (&(&_S1211)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1622 = _S1618.crf_params_1[int(2)].toe_0 + _S1415.differential_0;
    (&(&_S1211)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1623 = _S1618.crf_params_1[int(1)].center_0 + _S1448;
    (&(&_S1211)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1624 = _S1618.crf_params_1[int(1)].gamma_0 + _S1451.differential_0;
    (&(&_S1211)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1625 = _S1618.crf_params_1[int(1)].shoulder_0 + _S1455.differential_0;
    (&(&_S1211)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1626 = _S1618.crf_params_1[int(1)].toe_0 + _S1459.differential_0;
    (&(&_S1211)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1627 = _S1618.crf_params_1[int(0)].center_0 + _S1492;
    (&(&_S1211)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1628 = _S1618.crf_params_1[int(0)].gamma_0 + _S1495.differential_0;
    (&(&_S1211)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1629 = _S1618.crf_params_1[int(0)].shoulder_0 + _S1499.differential_0;
    (&(&_S1211)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1630 = _S1618.crf_params_1[int(0)].toe_0 + _S1503.differential_0;
    *&((&(&(&_S1211)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1211)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1211)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1631 = _S1577 + _S1618.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1211)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1632 = _S1576 + _S1618.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1211)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1633 = _S1575 + _S1618.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1211)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1634 = _S1582 + _S1618.vignette_params_1[int(2)].cy_0;
    (&(&_S1211)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1635 = _S1583 + _S1618.vignette_params_1[int(2)].cx_0;
    (&(&_S1211)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1636 = _S1591 + _S1618.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1211)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1637 = _S1590 + _S1618.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1211)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1638 = _S1589 + _S1618.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1211)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1639 = _S1596 + _S1618.vignette_params_1[int(1)].cy_0;
    (&(&_S1211)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1640 = _S1597 + _S1618.vignette_params_1[int(1)].cx_0;
    (&(&_S1211)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1641 = _S1605 + _S1618.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1211)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1642 = _S1604 + _S1618.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1211)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1643 = _S1603 + _S1618.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1211)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1644 = _S1610 + _S1618.vignette_params_1[int(0)].cy_0;
    (&(&_S1211)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1645 = _S1611 + _S1618.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S1646;
    _S1646[int(0)] = 0.0f;
    _S1646[int(1)] = 0.0f;
    _S1646[int(2)] = 0.0f;
    _S1646[int(3)] = 0.0f;
    _S1646[int(4)] = 0.0f;
    _S1646[int(5)] = 0.0f;
    _S1646[int(6)] = 0.0f;
    _S1646[int(7)] = 0.0f;
    _S1646[int(8)] = 0.0f;
    _S1646[int(9)] = 0.0f;
    _S1646[int(10)] = 0.0f;
    _S1646[int(11)] = 0.0f;
    _S1646[int(12)] = 0.0f;
    _S1646[int(13)] = 0.0f;
    _S1646[int(14)] = 0.0f;
    _S1646[int(15)] = 0.0f;
    _S1646[int(16)] = 0.0f;
    _S1646[int(17)] = 0.0f;
    _S1646[int(18)] = 0.0f;
    _S1646[int(19)] = 0.0f;
    _S1646[int(20)] = 0.0f;
    _S1646[int(21)] = 0.0f;
    _S1646[int(22)] = 0.0f;
    _S1646[int(23)] = 0.0f;
    _S1646[int(24)] = 0.0f;
    _S1646[int(25)] = 0.0f;
    _S1646[int(26)] = 0.0f;
    _S1646[int(27)] = 0.0f;
    _S1646[int(28)] = 0.0f;
    _S1646[int(29)] = 0.0f;
    _S1646[int(30)] = 0.0f;
    _S1646[int(31)] = 0.0f;
    _S1646[int(32)] = 0.0f;
    _S1646[int(33)] = 0.0f;
    _S1646[int(34)] = 0.0f;
    _S1646[int(35)] = 0.0f;
    _S1646[int(8)] = _S1638;
    _S1646[int(16)] = _S1618.color_params_1.b_0.x;
    _S1646[int(15)] = _S1631;
    _S1646[int(14)] = _S1632;
    _S1646[int(13)] = _S1633;
    _S1646[int(12)] = _S1634;
    _S1646[int(11)] = _S1635;
    _S1646[int(10)] = _S1636;
    _S1646[int(9)] = _S1637;
    _S1646[int(17)] = _S1618.color_params_1.b_0.y;
    _S1646[int(7)] = _S1639;
    _S1646[int(6)] = _S1640;
    _S1646[int(5)] = _S1641;
    _S1646[int(4)] = _S1642;
    _S1646[int(3)] = _S1643;
    _S1646[int(2)] = _S1644;
    _S1646[int(1)] = _S1645;
    _S1646[int(0)] = _S1211.exposure_1;
    _S1646[int(26)] = _S1628;
    _S1646[int(34)] = _S1620;
    _S1646[int(33)] = _S1621;
    _S1646[int(32)] = _S1622;
    _S1646[int(31)] = _S1623;
    _S1646[int(30)] = _S1624;
    _S1646[int(29)] = _S1625;
    _S1646[int(28)] = _S1626;
    _S1646[int(27)] = _S1627;
    _S1646[int(35)] = _S1619;
    _S1646[int(25)] = _S1629;
    _S1646[int(24)] = _S1630;
    _S1646[int(23)] = _S1618.color_params_1.n_0.y;
    _S1646[int(22)] = _S1618.color_params_1.n_0.x;
    _S1646[int(21)] = _S1618.color_params_1.g_0.y;
    _S1646[int(20)] = _S1618.color_params_1.g_0.x;
    _S1646[int(19)] = _S1618.color_params_1.r_0.y;
    _S1646[int(18)] = _S1618.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1646;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1615;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1647, float2  _S1648, float2  _S1649, float2  _S1650, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1651, float3  _S1652)
{
    s_bwd_prop_apply_ppisp_0(_S1647, _S1648, _S1649, _S1650, _S1651, _S1652);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  * params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1653 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S1653;
    FixedArray<float, 36>  _S1654 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = *params_2;
    (&dp_params_0)->differential_0 = _S1654;
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
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1655 = *dprgb_in_1;
    float3  _S1656 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1657 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1658 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1659 = {
        _S1658, _S1658, _S1658
    };
    float2  _S1660 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1661 = { _S1660, _S1660, _S1660, _S1660 };
    RQSCRFPPISPChannelParams_0 _S1662 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1663 = {
        _S1662, _S1662, _S1662
    };
    PPISPParamsRQS_0 _S1664;
    (&_S1664)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1664)->vignette_params_0 = _S1659;
    (&_S1664)->color_params_0 = _S1661;
    (&_S1664)->crf_params_0 = _S1663;
    (&(&_S1664)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1664)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S1665 = dpparams_1->primal_0[int(3)];
    (&(&_S1664)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S1666 = dpparams_1->primal_0[int(4)];
    (&(&_S1664)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S1667 = dpparams_1->primal_0[int(5)];
    (&(&_S1664)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1664)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1664)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S1668 = dpparams_1->primal_0[int(8)];
    (&(&_S1664)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S1669 = dpparams_1->primal_0[int(9)];
    (&(&_S1664)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S1670 = dpparams_1->primal_0[int(10)];
    (&(&_S1664)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1664)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1664)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S1671 = dpparams_1->primal_0[int(13)];
    (&(&_S1664)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S1672 = dpparams_1->primal_0[int(14)];
    (&(&_S1664)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S1673 = dpparams_1->primal_0[int(15)];
    (&(&_S1664)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1664)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1664)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1664)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1664)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1664)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1664)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1664)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1664)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S1674 = dpparams_1->primal_0[int(24)];
    (&(&_S1664)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S1675 = dpparams_1->primal_0[int(25)];
    (&(&_S1664)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S1676 = dpparams_1->primal_0[int(26)];
    (&(&_S1664)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S1677 = dpparams_1->primal_0[int(27)];
    (&(&_S1664)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S1678 = dpparams_1->primal_0[int(28)];
    (&(&_S1664)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S1679 = dpparams_1->primal_0[int(29)];
    (&(&_S1664)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S1680 = dpparams_1->primal_0[int(30)];
    (&(&_S1664)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S1681 = dpparams_1->primal_0[int(31)];
    (&(&_S1664)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S1682 = dpparams_1->primal_0[int(32)];
    (&(&_S1664)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S1683 = dpparams_1->primal_0[int(33)];
    (&(&_S1664)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S1684 = dpparams_1->primal_0[int(34)];
    (&(&_S1664)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S1685 = dpparams_1->primal_0[int(35)];
    (&(&_S1664)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S1686 = dpparams_1->primal_0[int(36)];
    (&(&_S1664)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S1687 = dpparams_1->primal_0[int(37)];
    (&(&_S1664)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S1688 = dpparams_1->primal_0[int(38)];
    (&(&_S1664)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S1689 = _S1664;
    float _S1690 = s_primal_ctx_exp2_0(_S1664.exposure_0);
    float3  _S1691 = make_float3 (_S1690);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S1690);
    float _S1692 = s_primal_ctx_max_0(img_size_4.x, img_size_4.y);
    float _S1693 = (pix_coord_4.x - image_center_4.x) / _S1692;
    float _S1694 = (pix_coord_4.y - image_center_4.y) / _S1692;
    float dx_11 = _S1693 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S1694 - dpparams_1->primal_0[int(2)];
    float r2_16 = dx_11 * dx_11 + dy_9 * dy_9;
    float r4_9 = r2_16 * r2_16;
    float r6_3 = r4_9 * r2_16;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_16 + 1.0f;
    float _S1695 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S1696 = rgb_out_5.x * _S1695;
    float3  _S1697 = rgb_out_5;
    *&((&_S1697)->x) = _S1696;
    float dx_12 = _S1693 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S1694 - dpparams_1->primal_0[int(7)];
    float r2_17 = dx_12 * dx_12 + dy_10 * dy_10;
    float r4_10 = r2_17 * r2_17;
    float r6_4 = r4_10 * r2_17;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_17 + 1.0f;
    float _S1698 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S1697)->y) = rgb_out_5.y * _S1698;
    float dx_13 = _S1693 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S1694 - dpparams_1->primal_0[int(12)];
    float r2_18 = dx_13 * dx_13 + dy_11 * dy_11;
    float r4_11 = r2_18 * r2_18;
    float r6_5 = r4_11 * r2_18;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_18 + 1.0f;
    float _S1699 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S1697)->z) = rgb_out_5.z * _S1699;
    PPISPParamsRQS_0 _S1700 = _S1664;
    float2  _S1701 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1664.color_params_0.b_0);
    float2  _S1702 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1664.color_params_0.r_0);
    float2  _S1703 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1664.color_params_0.g_0);
    float2  _S1704 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1664.color_params_0.n_0);
    float _S1705 = 0.3333333432674408f + _S1704.x;
    float _S1706 = 0.3333333432674408f + _S1704.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S1701.x, 1.0f + _S1702.x, _S1703.x, _S1701.y, _S1702.y, 1.0f + _S1703.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1706, 1.0f, 0.0f, - _S1705, - _S1706, _S1705, 0.0f);
    Matrix<float, 3, 3>  _S1707 = s_primal_ctx_mul_1(skew_1, T_3);
    float3  r0_3 = make_float3 (_S1707.rows[int(0)].x, _S1707.rows[int(0)].y, _S1707.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S1707.rows[int(1)].x, _S1707.rows[int(1)].y, _S1707.rows[int(1)].z);
    float3  r2_19 = make_float3 (_S1707.rows[int(2)].x, _S1707.rows[int(2)].y, _S1707.rows[int(2)].z);
    float3  _S1708 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S1709 = (s_primal_ctx_dot_0(_S1708, _S1708)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S1710;
    bool _S1711;
    if(_S1709)
    {
        float3  _S1712 = s_primal_ctx_cross_0(r0_3, r2_19);
        bool _S1713 = (s_primal_ctx_dot_0(_S1712, _S1712)) < 9.99999968265522539e-21f;
        if(_S1713)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_19);
        }
        else
        {
            lambda_v_7 = _S1712;
        }
        _S1711 = _S1713;
        _S1710 = _S1712;
    }
    else
    {
        lambda_v_7 = _S1708;
        _S1711 = false;
        _S1710 = _S1656;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S1714 = s_primal_ctx_mul_1(T_3, D_1);
    Matrix<float, 3, 3>  _S1715 = s_primal_ctx_mul_1(_S1714, S_inv_1);
    bool _S1716 = (s_primal_ctx_abs_0(_S1715.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S1717;
    float _S1718;
    if(_S1716)
    {
        float inv_s_1 = 1.0f / _S1715.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1719 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S1720 = _S1715.rows[int(2)].z * _S1715.rows[int(2)].z;
        H_5 = _S1715 * makeMatrix<float, 3, 3> (inv_s_1);
        _S1717 = _S1719;
        _S1718 = _S1720;
    }
    else
    {
        H_5 = _S1715;
        _S1717 = _S1657;
        _S1718 = 0.0f;
    }
    float _S1721 = _S1697.x;
    float _S1722 = _S1697.y;
    float intensity_3 = _S1721 + _S1722 + _S1697.z;
    float3  rgi_in_1 = make_float3 (_S1721, _S1722, intensity_3);
    float3  _S1723 = s_primal_ctx_mul_2(H_5, rgi_in_1);
    float _S1724 = _S1723.z + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S1724;
    float3  _S1725 = make_float3 (norm_factor_1);
    float _S1726 = _S1724 * _S1724;
    float3  rgi_out_5 = _S1723 * make_float3 (norm_factor_1);
    float _S1727 = rgi_out_5.x;
    float _S1728 = rgi_out_5.y;
    float3  _S1729 = make_float3 (_S1727, _S1728, rgi_out_5.z - _S1727 - _S1728);
    float3  _S1730 = make_float3 (0.0f);
    float3  _S1731 = make_float3 (1.0f);
    float3  _S1732 = s_primal_ctx_clamp_1(_S1729, _S1730, _S1731);
    float _S1733 = _S1732.x;
    float _S1734 = s_primal_ctx_exp_0(_S1674);
    float _S1735 = s_primal_ctx_exp_0(_S1675);
    float _S1736 = - _S1676;
    float _S1737 = 1.0f + s_primal_ctx_exp_0(_S1736);
    float x0_4 = 1.0f / _S1737;
    float _S1738 = _S1737 * _S1737;
    float _S1739 = - _S1677;
    float _S1740 = 1.0f + s_primal_ctx_exp_0(_S1739);
    float y0_4 = 1.0f / _S1740;
    float _S1741 = _S1740 * _S1740;
    float _S1742 = s_primal_ctx_exp_0(_S1678);
    bool _S1743 = _S1733 < x0_4;
    float _S1744;
    float _S1745;
    float _S1746;
    float _S1747;
    float _S1748;
    float _S1749;
    float _S1750;
    float _S1751;
    float _S1752;
    float _S1753;
    float _S1754;
    float _S1755;
    float _S1756;
    float _S1757;
    float _S1758;
    float _S1759;
    float _S1760;
    float _S1761;
    float _S1762;
    float _S1763;
    float _S1764;
    float _S1765;
    float _S1766;
    float _S1767;
    float _S1768;
    float _S1769;
    float _S1770;
    if(_S1743)
    {
        float s0_3 = y0_4 / x0_4;
        float _S1771 = x0_4 * x0_4;
        float t0_3 = _S1733 / x0_4;
        float _S1772 = s0_3 * t0_3;
        float _S1773 = _S1734 * t0_3;
        float _S1774 = 1.0f - t0_3;
        float _S1775 = _S1772 * t0_3 + _S1773 * _S1774;
        float _S1776 = y0_4 * _S1775;
        float _S1777 = _S1734 + _S1742 - 2.0f * s0_3;
        float _S1778 = _S1777 * t0_3;
        float _S1779 = s0_3 + _S1778 * _S1774;
        float _S1780 = _S1779 * _S1779;
        _S1744 = 0.0f;
        _S1745 = 0.0f;
        _S1746 = 0.0f;
        _S1747 = 0.0f;
        _S1748 = 0.0f;
        _S1749 = 0.0f;
        _S1750 = 0.0f;
        _S1751 = 0.0f;
        _S1752 = 0.0f;
        _S1753 = 0.0f;
        _S1754 = 0.0f;
        _S1755 = 0.0f;
        _S1756 = 0.0f;
        _S1757 = 0.0f;
        _S1758 = 0.0f;
        _S1759 = _S1780;
        _S1760 = _S1776;
        _S1761 = _S1779;
        _S1762 = _S1778;
        _S1763 = _S1774;
        _S1764 = _S1777;
        _S1765 = t0_3;
        _S1766 = _S1775;
        _S1767 = _S1773;
        _S1768 = _S1772;
        _S1769 = s0_3;
        _S1770 = _S1771;
    }
    else
    {
        float _S1781 = 1.0f - y0_4;
        float _S1782 = 1.0f - x0_4;
        float s1_7 = _S1781 / _S1782;
        float _S1783 = _S1782 * _S1782;
        float _S1784 = _S1733 - x0_4;
        float t1_3 = _S1784 / _S1782;
        float _S1785 = s1_7 * t1_3;
        float _S1786 = _S1742 * t1_3;
        float _S1787 = 1.0f - t1_3;
        float _S1788 = _S1785 * t1_3 + _S1786 * _S1787;
        float _S1789 = _S1781 * _S1788;
        float _S1790 = _S1742 + _S1735 - 2.0f * s1_7;
        float _S1791 = _S1790 * t1_3;
        float _S1792 = s1_7 + _S1791 * _S1787;
        _S1744 = _S1792 * _S1792;
        _S1745 = _S1789;
        _S1746 = _S1792;
        _S1747 = _S1791;
        _S1748 = _S1787;
        _S1749 = _S1790;
        _S1750 = t1_3;
        _S1751 = _S1781;
        _S1752 = _S1788;
        _S1753 = _S1786;
        _S1754 = _S1785;
        _S1755 = s1_7;
        _S1756 = _S1783;
        _S1757 = _S1784;
        _S1758 = _S1782;
        _S1759 = 0.0f;
        _S1760 = 0.0f;
        _S1761 = 0.0f;
        _S1762 = 0.0f;
        _S1763 = 0.0f;
        _S1764 = 0.0f;
        _S1765 = 0.0f;
        _S1766 = 0.0f;
        _S1767 = 0.0f;
        _S1768 = 0.0f;
        _S1769 = 0.0f;
        _S1770 = 0.0f;
    }
    float _S1793 = _S1732.y;
    float _S1794 = s_primal_ctx_exp_0(_S1679);
    float _S1795 = s_primal_ctx_exp_0(_S1680);
    float _S1796 = - _S1681;
    float _S1797 = 1.0f + s_primal_ctx_exp_0(_S1796);
    float x0_5 = 1.0f / _S1797;
    float _S1798 = _S1797 * _S1797;
    float _S1799 = - _S1682;
    float _S1800 = 1.0f + s_primal_ctx_exp_0(_S1799);
    float y0_5 = 1.0f / _S1800;
    float _S1801 = _S1800 * _S1800;
    float _S1802 = s_primal_ctx_exp_0(_S1683);
    bool _S1803 = _S1793 < x0_5;
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
    float _S1818;
    float _S1819;
    float _S1820;
    float _S1821;
    float _S1822;
    float _S1823;
    float _S1824;
    float _S1825;
    float _S1826;
    float _S1827;
    float _S1828;
    float _S1829;
    float _S1830;
    if(_S1803)
    {
        float s0_4 = y0_5 / x0_5;
        float _S1831 = x0_5 * x0_5;
        float t0_4 = _S1793 / x0_5;
        float _S1832 = s0_4 * t0_4;
        float _S1833 = _S1794 * t0_4;
        float _S1834 = 1.0f - t0_4;
        float _S1835 = _S1832 * t0_4 + _S1833 * _S1834;
        float _S1836 = y0_5 * _S1835;
        float _S1837 = _S1794 + _S1802 - 2.0f * s0_4;
        float _S1838 = _S1837 * t0_4;
        float _S1839 = s0_4 + _S1838 * _S1834;
        float _S1840 = _S1839 * _S1839;
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
        _S1818 = 0.0f;
        _S1819 = _S1840;
        _S1820 = _S1836;
        _S1821 = _S1839;
        _S1822 = _S1838;
        _S1823 = _S1834;
        _S1824 = _S1837;
        _S1825 = t0_4;
        _S1826 = _S1835;
        _S1827 = _S1833;
        _S1828 = _S1832;
        _S1829 = s0_4;
        _S1830 = _S1831;
    }
    else
    {
        float _S1841 = 1.0f - y0_5;
        float _S1842 = 1.0f - x0_5;
        float s1_8 = _S1841 / _S1842;
        float _S1843 = _S1842 * _S1842;
        float _S1844 = _S1793 - x0_5;
        float t1_4 = _S1844 / _S1842;
        float _S1845 = s1_8 * t1_4;
        float _S1846 = _S1802 * t1_4;
        float _S1847 = 1.0f - t1_4;
        float _S1848 = _S1845 * t1_4 + _S1846 * _S1847;
        float _S1849 = _S1841 * _S1848;
        float _S1850 = _S1802 + _S1795 - 2.0f * s1_8;
        float _S1851 = _S1850 * t1_4;
        float _S1852 = s1_8 + _S1851 * _S1847;
        _S1804 = _S1852 * _S1852;
        _S1805 = _S1849;
        _S1806 = _S1852;
        _S1807 = _S1851;
        _S1808 = _S1847;
        _S1809 = _S1850;
        _S1810 = t1_4;
        _S1811 = _S1841;
        _S1812 = _S1848;
        _S1813 = _S1846;
        _S1814 = _S1845;
        _S1815 = s1_8;
        _S1816 = _S1843;
        _S1817 = _S1844;
        _S1818 = _S1842;
        _S1819 = 0.0f;
        _S1820 = 0.0f;
        _S1821 = 0.0f;
        _S1822 = 0.0f;
        _S1823 = 0.0f;
        _S1824 = 0.0f;
        _S1825 = 0.0f;
        _S1826 = 0.0f;
        _S1827 = 0.0f;
        _S1828 = 0.0f;
        _S1829 = 0.0f;
        _S1830 = 0.0f;
    }
    float _S1853 = _S1732.z;
    float _S1854 = s_primal_ctx_exp_0(_S1684);
    float _S1855 = s_primal_ctx_exp_0(_S1685);
    float _S1856 = - _S1686;
    float _S1857 = 1.0f + s_primal_ctx_exp_0(_S1856);
    float x0_6 = 1.0f / _S1857;
    float _S1858 = _S1857 * _S1857;
    float _S1859 = - _S1687;
    float _S1860 = 1.0f + s_primal_ctx_exp_0(_S1859);
    float y0_6 = 1.0f / _S1860;
    float _S1861 = _S1860 * _S1860;
    float _S1862 = s_primal_ctx_exp_0(_S1688);
    bool _S1863 = _S1853 < x0_6;
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
    float _S1878;
    float _S1879;
    float _S1880;
    float _S1881;
    float _S1882;
    float _S1883;
    float _S1884;
    float _S1885;
    float _S1886;
    float _S1887;
    float _S1888;
    float _S1889;
    float _S1890;
    if(_S1863)
    {
        float s0_5 = y0_6 / x0_6;
        float _S1891 = x0_6 * x0_6;
        float t0_5 = _S1853 / x0_6;
        float _S1892 = s0_5 * t0_5;
        float _S1893 = _S1854 * t0_5;
        float _S1894 = 1.0f - t0_5;
        float _S1895 = _S1892 * t0_5 + _S1893 * _S1894;
        float _S1896 = y0_6 * _S1895;
        float _S1897 = _S1854 + _S1862 - 2.0f * s0_5;
        float _S1898 = _S1897 * t0_5;
        float _S1899 = s0_5 + _S1898 * _S1894;
        float _S1900 = _S1899 * _S1899;
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
        _S1878 = 0.0f;
        _S1879 = _S1900;
        _S1880 = _S1896;
        _S1881 = _S1899;
        _S1882 = _S1898;
        _S1883 = _S1894;
        _S1884 = _S1897;
        _S1885 = t0_5;
        _S1886 = _S1895;
        _S1887 = _S1893;
        _S1888 = _S1892;
        _S1889 = s0_5;
        _S1890 = _S1891;
    }
    else
    {
        float _S1901 = 1.0f - y0_6;
        float _S1902 = 1.0f - x0_6;
        float s1_9 = _S1901 / _S1902;
        float _S1903 = _S1902 * _S1902;
        float _S1904 = _S1853 - x0_6;
        float t1_5 = _S1904 / _S1902;
        float _S1905 = s1_9 * t1_5;
        float _S1906 = _S1862 * t1_5;
        float _S1907 = 1.0f - t1_5;
        float _S1908 = _S1905 * t1_5 + _S1906 * _S1907;
        float _S1909 = _S1901 * _S1908;
        float _S1910 = _S1862 + _S1855 - 2.0f * s1_9;
        float _S1911 = _S1910 * t1_5;
        float _S1912 = s1_9 + _S1911 * _S1907;
        _S1864 = _S1912 * _S1912;
        _S1865 = _S1909;
        _S1866 = _S1912;
        _S1867 = _S1911;
        _S1868 = _S1907;
        _S1869 = _S1910;
        _S1870 = t1_5;
        _S1871 = _S1901;
        _S1872 = _S1908;
        _S1873 = _S1906;
        _S1874 = _S1905;
        _S1875 = s1_9;
        _S1876 = _S1903;
        _S1877 = _S1904;
        _S1878 = _S1902;
        _S1879 = 0.0f;
        _S1880 = 0.0f;
        _S1881 = 0.0f;
        _S1882 = 0.0f;
        _S1883 = 0.0f;
        _S1884 = 0.0f;
        _S1885 = 0.0f;
        _S1886 = 0.0f;
        _S1887 = 0.0f;
        _S1888 = 0.0f;
        _S1889 = 0.0f;
        _S1890 = 0.0f;
    }
    if(_S1863)
    {
        float _S1913 = _s_dOut_12.z / _S1879;
        float _S1914 = _S1880 * - _S1913;
        float _S1915 = _S1881 * _S1913;
        float _S1916 = _S1883 * _S1914;
        float _S1917 = _S1885 * _S1916;
        float _S1918 = y0_6 * _S1915;
        float _S1919 = _S1883 * _S1918;
        float _S1920 = _S1885 * _S1918;
        float _S1921 = (_S1884 * _S1916 + - (_S1882 * _S1914 + _S1887 * _S1918) + _S1854 * _S1919 + _S1888 * _S1918 + _S1889 * _S1920) / _S1890;
        float _S1922 = x0_6 * _S1921;
        float _S1923 = (_S1914 + 2.0f * - _S1917 + _S1885 * _S1920) / _S1890;
        float _S1924 = _S1886 * _S1915 + x0_6 * _S1923;
        float _S1925 = _S1917 + _S1885 * _S1919;
        float _S1926 = _S1853 * - _S1921 + y0_6 * - _S1923;
        _S1864 = _S1917;
        _S1865 = _S1924;
        _S1866 = _S1926;
        _S1867 = 0.0f;
        _S1868 = _S1925;
        _S1869 = _S1922;
    }
    else
    {
        float _S1927 = _s_dOut_12.z / _S1864;
        float _S1928 = _S1865 * - _S1927;
        float _S1929 = _S1866 * _S1927;
        float _S1930 = _S1868 * _S1928;
        float _S1931 = _S1870 * _S1930;
        float _S1932 = _S1871 * _S1929;
        float _S1933 = _S1868 * _S1932;
        float _S1934 = _S1870 * _S1932;
        float _S1935 = (_S1869 * _S1930 + - (_S1867 * _S1928 + _S1873 * _S1932) + _S1862 * _S1933 + _S1874 * _S1932 + _S1875 * _S1934) / _S1876;
        float _S1936 = _S1878 * _S1935;
        float _S1937 = (_S1928 + 2.0f * - _S1931 + _S1870 * _S1934) / _S1876;
        float _S1938 = _s_dOut_12.z + - (_S1872 * _S1929 + _S1878 * _S1937);
        float _S1939 = - _S1936 + - (_S1877 * - _S1935 + _S1871 * - _S1937);
        _S1864 = _S1931 + _S1870 * _S1933;
        _S1865 = _S1938;
        _S1866 = _S1939;
        _S1867 = _S1931;
        _S1868 = 0.0f;
        _S1869 = _S1936;
    }
    DiffPair_float_0 _S1940;
    (&_S1940)->primal_0 = _S1688;
    (&_S1940)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1940, _S1864);
    DiffPair_float_0 _S1941 = _S1940;
    float _S1942 = - (_S1865 / _S1861);
    DiffPair_float_0 _S1943;
    (&_S1943)->primal_0 = _S1859;
    (&_S1943)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1943, _S1942);
    float _S1944 = - _S1943.differential_0;
    float _S1945 = - (_S1866 / _S1858);
    DiffPair_float_0 _S1946;
    (&_S1946)->primal_0 = _S1856;
    (&_S1946)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1946, _S1945);
    float _S1947 = - _S1946.differential_0;
    DiffPair_float_0 _S1948;
    (&_S1948)->primal_0 = _S1685;
    (&_S1948)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1948, _S1867);
    DiffPair_float_0 _S1949 = _S1948;
    DiffPair_float_0 _S1950;
    (&_S1950)->primal_0 = _S1684;
    (&_S1950)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1950, _S1868);
    DiffPair_float_0 _S1951 = _S1950;
    float3  _S1952 = make_float3 (0.0f, 0.0f, _S1869);
    if(_S1803)
    {
        float _S1953 = _s_dOut_12.y / _S1819;
        float _S1954 = _S1820 * - _S1953;
        float _S1955 = _S1821 * _S1953;
        float _S1956 = _S1823 * _S1954;
        float _S1957 = _S1825 * _S1956;
        float _S1958 = y0_5 * _S1955;
        float _S1959 = _S1823 * _S1958;
        float _S1960 = _S1825 * _S1958;
        float _S1961 = (_S1824 * _S1956 + - (_S1822 * _S1954 + _S1827 * _S1958) + _S1794 * _S1959 + _S1828 * _S1958 + _S1829 * _S1960) / _S1830;
        float _S1962 = x0_5 * _S1961;
        float _S1963 = (_S1954 + 2.0f * - _S1957 + _S1825 * _S1960) / _S1830;
        float _S1964 = _S1826 * _S1955 + x0_5 * _S1963;
        float _S1965 = _S1957 + _S1825 * _S1959;
        float _S1966 = _S1793 * - _S1961 + y0_5 * - _S1963;
        _S1804 = _S1957;
        _S1805 = _S1964;
        _S1806 = _S1966;
        _S1807 = 0.0f;
        _S1808 = _S1965;
        _S1809 = _S1962;
    }
    else
    {
        float _S1967 = _s_dOut_12.y / _S1804;
        float _S1968 = _S1805 * - _S1967;
        float _S1969 = _S1806 * _S1967;
        float _S1970 = _S1808 * _S1968;
        float _S1971 = _S1810 * _S1970;
        float _S1972 = _S1811 * _S1969;
        float _S1973 = _S1808 * _S1972;
        float _S1974 = _S1810 * _S1972;
        float _S1975 = (_S1809 * _S1970 + - (_S1807 * _S1968 + _S1813 * _S1972) + _S1802 * _S1973 + _S1814 * _S1972 + _S1815 * _S1974) / _S1816;
        float _S1976 = _S1818 * _S1975;
        float _S1977 = (_S1968 + 2.0f * - _S1971 + _S1810 * _S1974) / _S1816;
        float _S1978 = _s_dOut_12.y + - (_S1812 * _S1969 + _S1818 * _S1977);
        float _S1979 = - _S1976 + - (_S1817 * - _S1975 + _S1811 * - _S1977);
        _S1804 = _S1971 + _S1810 * _S1973;
        _S1805 = _S1978;
        _S1806 = _S1979;
        _S1807 = _S1971;
        _S1808 = 0.0f;
        _S1809 = _S1976;
    }
    DiffPair_float_0 _S1980;
    (&_S1980)->primal_0 = _S1683;
    (&_S1980)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1980, _S1804);
    DiffPair_float_0 _S1981 = _S1980;
    float _S1982 = - (_S1805 / _S1801);
    DiffPair_float_0 _S1983;
    (&_S1983)->primal_0 = _S1799;
    (&_S1983)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1983, _S1982);
    float _S1984 = - _S1983.differential_0;
    float _S1985 = - (_S1806 / _S1798);
    DiffPair_float_0 _S1986;
    (&_S1986)->primal_0 = _S1796;
    (&_S1986)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1986, _S1985);
    float _S1987 = - _S1986.differential_0;
    DiffPair_float_0 _S1988;
    (&_S1988)->primal_0 = _S1680;
    (&_S1988)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1988, _S1807);
    DiffPair_float_0 _S1989 = _S1988;
    DiffPair_float_0 _S1990;
    (&_S1990)->primal_0 = _S1679;
    (&_S1990)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1990, _S1808);
    DiffPair_float_0 _S1991 = _S1990;
    float3  _S1992 = _S1952 + make_float3 (0.0f, _S1809, 0.0f);
    if(_S1743)
    {
        float _S1993 = _s_dOut_12.x / _S1759;
        float _S1994 = _S1760 * - _S1993;
        float _S1995 = _S1761 * _S1993;
        float _S1996 = _S1763 * _S1994;
        float _S1997 = _S1765 * _S1996;
        float _S1998 = y0_4 * _S1995;
        float _S1999 = _S1763 * _S1998;
        float _S2000 = _S1765 * _S1998;
        float _S2001 = (_S1764 * _S1996 + - (_S1762 * _S1994 + _S1767 * _S1998) + _S1734 * _S1999 + _S1768 * _S1998 + _S1769 * _S2000) / _S1770;
        float _S2002 = x0_4 * _S2001;
        float _S2003 = (_S1994 + 2.0f * - _S1997 + _S1765 * _S2000) / _S1770;
        float _S2004 = _S1766 * _S1995 + x0_4 * _S2003;
        float _S2005 = _S1997 + _S1765 * _S1999;
        float _S2006 = _S1733 * - _S2001 + y0_4 * - _S2003;
        _S1744 = _S1997;
        _S1745 = _S2004;
        _S1746 = _S2006;
        _S1747 = 0.0f;
        _S1748 = _S2005;
        _S1749 = _S2002;
    }
    else
    {
        float _S2007 = _s_dOut_12.x / _S1744;
        float _S2008 = _S1745 * - _S2007;
        float _S2009 = _S1746 * _S2007;
        float _S2010 = _S1748 * _S2008;
        float _S2011 = _S1750 * _S2010;
        float _S2012 = _S1751 * _S2009;
        float _S2013 = _S1748 * _S2012;
        float _S2014 = _S1750 * _S2012;
        float _S2015 = (_S1749 * _S2010 + - (_S1747 * _S2008 + _S1753 * _S2012) + _S1742 * _S2013 + _S1754 * _S2012 + _S1755 * _S2014) / _S1756;
        float _S2016 = _S1758 * _S2015;
        float _S2017 = (_S2008 + 2.0f * - _S2011 + _S1750 * _S2014) / _S1756;
        float _S2018 = _s_dOut_12.x + - (_S1752 * _S2009 + _S1758 * _S2017);
        float _S2019 = - _S2016 + - (_S1757 * - _S2015 + _S1751 * - _S2017);
        _S1744 = _S2011 + _S1750 * _S2013;
        _S1745 = _S2018;
        _S1746 = _S2019;
        _S1747 = _S2011;
        _S1748 = 0.0f;
        _S1749 = _S2016;
    }
    DiffPair_float_0 _S2020;
    (&_S2020)->primal_0 = _S1678;
    (&_S2020)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2020, _S1744);
    DiffPair_float_0 _S2021 = _S2020;
    float _S2022 = - (_S1745 / _S1741);
    DiffPair_float_0 _S2023;
    (&_S2023)->primal_0 = _S1739;
    (&_S2023)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2023, _S2022);
    float _S2024 = - _S2023.differential_0;
    float _S2025 = - (_S1746 / _S1738);
    DiffPair_float_0 _S2026;
    (&_S2026)->primal_0 = _S1736;
    (&_S2026)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2026, _S2025);
    float _S2027 = - _S2026.differential_0;
    DiffPair_float_0 _S2028;
    (&_S2028)->primal_0 = _S1675;
    (&_S2028)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2028, _S1747);
    DiffPair_float_0 _S2029 = _S2028;
    DiffPair_float_0 _S2030;
    (&_S2030)->primal_0 = _S1674;
    (&_S2030)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2030, _S1748);
    DiffPair_float_0 _S2031 = _S2030;
    float3  _S2032 = _S1992 + make_float3 (_S1749, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2033;
    (&_S2033)->primal_0 = _S1729;
    (&_S2033)->differential_0 = _S1656;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2034;
    (&_S2034)->primal_0 = _S1730;
    (&_S2034)->differential_0 = _S1656;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2035;
    (&_S2035)->primal_0 = _S1731;
    (&_S2035)->differential_0 = _S1656;
    s_bwd_prop_clamp_1(&_S2033, &_S2034, &_S2035, _S2032);
    float _S2036 = - _S2033.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S2033.differential_0.x + _S2036, _S2033.differential_0.y + _S2036, _S2033.differential_0.z);
    float3  _S2037 = _S1723 * s_diff_rgi_out_T_1;
    float _S2038 = (_S2037.x + _S2037.y + _S2037.z) / _S1726;
    float _S2039 = _S1724 * _S2038;
    float3  _S2040 = _S1725 * s_diff_rgi_out_T_1 + make_float3 (0.0f, 0.0f, intensity_3 * - _S2038);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2041;
    (&_S2041)->primal_0 = H_5;
    (&_S2041)->differential_0 = _S1657;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2042;
    (&_S2042)->primal_0 = rgi_in_1;
    (&_S2042)->differential_0 = _S1656;
    s_bwd_prop_mul_1(&_S2041, &_S2042, _S2040);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2043 = _S2041;
    float _S2044 = _S2039 + _S2042.differential_0.z;
    float _S2045 = _S2042.differential_0.y + _S2044;
    float _S2046 = _S2042.differential_0.x + _S2044;
    float3  _S2047 = make_float3 (_S2046, _S2045, _S2044);
    if(_S1716)
    {
        Matrix<float, 3, 3>  _S2048 = _S1715 * _S2043.differential_0;
        Matrix<float, 3, 3>  _S2049 = _S1717 * _S2043.differential_0;
        _S1718 = - ((_S2048.rows[int(0)].x + _S2048.rows[int(0)].y + _S2048.rows[int(0)].z + _S2048.rows[int(1)].x + _S2048.rows[int(1)].y + _S2048.rows[int(1)].z + _S2048.rows[int(2)].x + _S2048.rows[int(2)].y + _S2048.rows[int(2)].z) / _S1718);
        H_5 = _S2049;
    }
    else
    {
        _S1718 = 0.0f;
        H_5 = _S2043.differential_0;
    }
    DiffPair_float_0 _S2050;
    (&_S2050)->primal_0 = _S1715.rows[int(2)].z;
    (&_S2050)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2050, 0.0f);
    float _S2051 = _S2050.differential_0 + _S1718;
    float3  _S2052 = _S1656;
    *&((&_S2052)->z) = _S2051;
    Matrix<float, 3, 3>  _S2053 = _S1657;
    _S2053[int(2)] = _S2052;
    Matrix<float, 3, 3>  _S2054 = H_5 + _S2053;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2055;
    (&_S2055)->primal_0 = _S1714;
    (&_S2055)->differential_0 = _S1657;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2056;
    (&_S2056)->primal_0 = S_inv_1;
    (&_S2056)->differential_0 = _S1657;
    s_bwd_prop_mul_2(&_S2055, &_S2056, _S2054);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2057;
    (&_S2057)->primal_0 = T_3;
    (&_S2057)->differential_0 = _S1657;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2058;
    (&_S2058)->primal_0 = D_1;
    (&_S2058)->differential_0 = _S1657;
    s_bwd_prop_mul_2(&_S2057, &_S2058, _S2055.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2059 = _S2057;
    float3  _S2060 = make_float3 (_S2058.differential_0.rows[int(0)].x, _S2058.differential_0.rows[int(1)].y, _S2058.differential_0.rows[int(2)].z);
    float3  _S2061;
    if(_S1709)
    {
        if(_S1711)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2062;
            (&_S2062)->primal_0 = r1_3;
            (&_S2062)->differential_0 = _S1656;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2063;
            (&_S2063)->primal_0 = r2_19;
            (&_S2063)->differential_0 = _S1656;
            s_bwd_prop_cross_0(&_S2062, &_S2063, _S2060);
            _S1697 = _S1656;
            lambda_v_7 = _S2063.differential_0;
            _S2061 = _S2062.differential_0;
        }
        else
        {
            _S1697 = _S2060;
            lambda_v_7 = _S1656;
            _S2061 = _S1656;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2064;
        (&_S2064)->primal_0 = _S1710;
        (&_S2064)->differential_0 = _S1656;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2065;
        (&_S2065)->primal_0 = _S1710;
        (&_S2065)->differential_0 = _S1656;
        s_bwd_prop_dot_0(&_S2064, &_S2065, 0.0f);
        float3  _S2066 = _S2065.differential_0 + _S2064.differential_0 + _S1697;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2067;
        (&_S2067)->primal_0 = r0_3;
        (&_S2067)->differential_0 = _S1656;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2068;
        (&_S2068)->primal_0 = r2_19;
        (&_S2068)->differential_0 = _S1656;
        s_bwd_prop_cross_0(&_S2067, &_S2068, _S2066);
        float3  _S2069 = _S2068.differential_0 + lambda_v_7;
        _S1697 = _S1656;
        lambda_v_7 = _S2069;
        _S1710 = _S2061;
        _S2061 = _S2067.differential_0;
    }
    else
    {
        _S1697 = _S2060;
        lambda_v_7 = _S1656;
        _S1710 = _S1656;
        _S2061 = _S1656;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2070;
    (&_S2070)->primal_0 = _S1708;
    (&_S2070)->differential_0 = _S1656;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2071;
    (&_S2071)->primal_0 = _S1708;
    (&_S2071)->differential_0 = _S1656;
    s_bwd_prop_dot_0(&_S2070, &_S2071, 0.0f);
    float3  _S2072 = _S2071.differential_0 + _S2070.differential_0 + _S1697;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2073;
    (&_S2073)->primal_0 = r0_3;
    (&_S2073)->differential_0 = _S1656;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2074;
    (&_S2074)->primal_0 = r1_3;
    (&_S2074)->differential_0 = _S1656;
    s_bwd_prop_cross_0(&_S2073, &_S2074, _S2072);
    float3  _S2075 = _S1656;
    *&((&_S2075)->z) = lambda_v_7.z;
    *&((&_S2075)->y) = lambda_v_7.y;
    *&((&_S2075)->x) = lambda_v_7.x;
    float3  _S2076 = _S2074.differential_0 + _S1710;
    float3  _S2077 = _S1656;
    *&((&_S2077)->z) = _S2076.z;
    *&((&_S2077)->y) = _S2076.y;
    *&((&_S2077)->x) = _S2076.x;
    float3  _S2078 = _S2073.differential_0 + _S2061;
    float3  _S2079 = _S1656;
    *&((&_S2079)->z) = _S2078.z;
    *&((&_S2079)->y) = _S2078.y;
    *&((&_S2079)->x) = _S2078.x;
    Matrix<float, 3, 3>  _S2080 = _S1657;
    _S2080[int(2)] = _S2075;
    _S2080[int(1)] = _S2077;
    _S2080[int(0)] = _S2079;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2081;
    (&_S2081)->primal_0 = skew_1;
    (&_S2081)->differential_0 = _S1657;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2082;
    (&_S2082)->primal_0 = T_3;
    (&_S2082)->differential_0 = _S1657;
    s_bwd_prop_mul_2(&_S2081, &_S2082, _S2080);
    Matrix<float, 3, 3>  _S2083 = _S2082.differential_0 + _S2059.differential_0;
    float2  _S2084 = make_float2 (_S2081.differential_0.rows[int(2)].y + - _S2081.differential_0.rows[int(1)].z, _S2081.differential_0.rows[int(0)].z + - _S2081.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S2085 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2086;
    (&_S2086)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2086)->differential_0 = _S2085;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2087;
    (&_S2087)->primal_0 = _S1700.color_params_0.n_0;
    (&_S2087)->differential_0 = _S1660;
    s_bwd_prop_mul_3(&_S2086, &_S2087, _S2084);
    float2  _S2088 = make_float2 (_S2083.rows[int(0)].z, _S2083.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2089;
    (&_S2089)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2089)->differential_0 = _S2085;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2090;
    (&_S2090)->primal_0 = _S1700.color_params_0.g_0;
    (&_S2090)->differential_0 = _S1660;
    s_bwd_prop_mul_3(&_S2089, &_S2090, _S2088);
    float2  _S2091 = make_float2 (_S2083.rows[int(0)].y, _S2083.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2092;
    (&_S2092)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2092)->differential_0 = _S2085;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2093;
    (&_S2093)->primal_0 = _S1700.color_params_0.r_0;
    (&_S2093)->differential_0 = _S1660;
    s_bwd_prop_mul_3(&_S2092, &_S2093, _S2091);
    float2  _S2094 = make_float2 (_S2083.rows[int(0)].x, _S2083.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2095;
    (&_S2095)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2095)->differential_0 = _S2085;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2096;
    (&_S2096)->primal_0 = _S1700.color_params_0.b_0;
    (&_S2096)->differential_0 = _S1660;
    s_bwd_prop_mul_3(&_S2095, &_S2096, _S2094);
    ColorPPISPParams_0 _S2097 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2097)->n_0 = _S2087.differential_0;
    (&_S2097)->g_0 = _S2090.differential_0;
    (&_S2097)->r_0 = _S2093.differential_0;
    (&_S2097)->b_0 = _S2096.differential_0;
    _S1697 = _S2047;
    *&((&_S1697)->z) = 0.0f;
    float _S2098 = rgb_out_5.z * _S2044;
    float _S2099 = _S1699 * _S2044;
    DiffPair_float_0 _S2100;
    (&_S2100)->primal_0 = falloff_5;
    (&_S2100)->differential_0 = 0.0f;
    DiffPair_float_0 _S2101;
    (&_S2101)->primal_0 = 0.0f;
    (&_S2101)->differential_0 = 0.0f;
    DiffPair_float_0 _S2102;
    (&_S2102)->primal_0 = 1.0f;
    (&_S2102)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2100, &_S2101, &_S2102, _S2098);
    float _S2103 = r2_18 * _S2100.differential_0;
    float _S2104 = r4_11 * _S2100.differential_0;
    float s_diff_r6_T_3 = _S1673 * _S2100.differential_0;
    float _S2105 = r6_5 * _S2100.differential_0;
    float _S2106 = r2_18 * (_S1672 * _S2100.differential_0 + r2_18 * s_diff_r6_T_3);
    float _S2107 = _S1671 * _S2100.differential_0 + r4_11 * s_diff_r6_T_3 + _S2106 + _S2106;
    float _S2108 = dy_11 * _S2107;
    float _S2109 = dx_13 * _S2107;
    float _S2110 = - (_S2108 + _S2108);
    float _S2111 = - (_S2109 + _S2109);
    *&((&_S1697)->y) = 0.0f;
    float _S2112 = rgb_out_5.y * _S2045;
    float _S2113 = _S1698 * _S2045;
    DiffPair_float_0 _S2114;
    (&_S2114)->primal_0 = falloff_4;
    (&_S2114)->differential_0 = 0.0f;
    DiffPair_float_0 _S2115;
    (&_S2115)->primal_0 = 0.0f;
    (&_S2115)->differential_0 = 0.0f;
    DiffPair_float_0 _S2116;
    (&_S2116)->primal_0 = 1.0f;
    (&_S2116)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2114, &_S2115, &_S2116, _S2112);
    float _S2117 = r2_17 * _S2114.differential_0;
    float _S2118 = r4_10 * _S2114.differential_0;
    float s_diff_r6_T_4 = _S1670 * _S2114.differential_0;
    float _S2119 = r6_4 * _S2114.differential_0;
    float _S2120 = r2_17 * (_S1669 * _S2114.differential_0 + r2_17 * s_diff_r6_T_4);
    float _S2121 = _S1668 * _S2114.differential_0 + r4_10 * s_diff_r6_T_4 + _S2120 + _S2120;
    float _S2122 = dy_10 * _S2121;
    float _S2123 = dx_12 * _S2121;
    float _S2124 = - (_S2122 + _S2122);
    float _S2125 = - (_S2123 + _S2123);
    *&((&_S1697)->x) = 0.0f;
    float _S2126 = rgb_out_5.x * _S2046;
    float _S2127 = _S1695 * _S2046;
    DiffPair_float_0 _S2128;
    (&_S2128)->primal_0 = falloff_3;
    (&_S2128)->differential_0 = 0.0f;
    DiffPair_float_0 _S2129;
    (&_S2129)->primal_0 = 0.0f;
    (&_S2129)->differential_0 = 0.0f;
    DiffPair_float_0 _S2130;
    (&_S2130)->primal_0 = 1.0f;
    (&_S2130)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2128, &_S2129, &_S2130, _S2126);
    float _S2131 = r2_16 * _S2128.differential_0;
    float _S2132 = r4_9 * _S2128.differential_0;
    float s_diff_r6_T_5 = _S1667 * _S2128.differential_0;
    float _S2133 = r6_3 * _S2128.differential_0;
    float _S2134 = r2_16 * (_S1666 * _S2128.differential_0 + r2_16 * s_diff_r6_T_5);
    float _S2135 = _S1665 * _S2128.differential_0 + r4_9 * s_diff_r6_T_5 + _S2134 + _S2134;
    float _S2136 = dy_9 * _S2135;
    float _S2137 = dx_11 * _S2135;
    float _S2138 = - (_S2136 + _S2136);
    float _S2139 = - (_S2137 + _S2137);
    float3  _S2140 = _S1656;
    *&((&_S2140)->z) = _S2099;
    *&((&_S2140)->y) = _S2113;
    *&((&_S2140)->x) = _S2127;
    float3  _S2141 = _S1697 + _S2140;
    float3  _S2142 = _S1655.primal_0 * _S2141;
    float3  _S2143 = _S1691 * _S2141;
    float _S2144 = _S2142.x + _S2142.y + _S2142.z;
    DiffPair_float_0 _S2145;
    (&_S2145)->primal_0 = _S1689.exposure_0;
    (&_S2145)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2145, _S2144);
    PPISPParamsRQS_0 _S2146 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2146)->color_params_0 = _S2097;
    (&_S2146)->exposure_0 = _S2145.differential_0;
    _S1664 = _S2146;
    (&(&_S1664)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2147 = _S2146.crf_params_0[int(2)].gc_0 + _S1941.differential_0;
    (&(&_S1664)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2148 = _S2146.crf_params_0[int(2)].y0_0 + _S1944;
    (&(&_S1664)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2149 = _S2146.crf_params_0[int(2)].x0_0 + _S1947;
    (&(&_S1664)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2150 = _S2146.crf_params_0[int(2)].g1_0 + _S1949.differential_0;
    (&(&_S1664)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2151 = _S2146.crf_params_0[int(2)].g0_0 + _S1951.differential_0;
    (&(&_S1664)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2152 = _S2146.crf_params_0[int(1)].gc_0 + _S1981.differential_0;
    (&(&_S1664)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2153 = _S2146.crf_params_0[int(1)].y0_0 + _S1984;
    (&(&_S1664)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2154 = _S2146.crf_params_0[int(1)].x0_0 + _S1987;
    (&(&_S1664)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2155 = _S2146.crf_params_0[int(1)].g1_0 + _S1989.differential_0;
    (&(&_S1664)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2156 = _S2146.crf_params_0[int(1)].g0_0 + _S1991.differential_0;
    (&(&_S1664)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2157 = _S2146.crf_params_0[int(0)].gc_0 + _S2021.differential_0;
    (&(&_S1664)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2158 = _S2146.crf_params_0[int(0)].y0_0 + _S2024;
    (&(&_S1664)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2159 = _S2146.crf_params_0[int(0)].x0_0 + _S2027;
    (&(&_S1664)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2160 = _S2146.crf_params_0[int(0)].g1_0 + _S2029.differential_0;
    (&(&_S1664)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2161 = _S2146.crf_params_0[int(0)].g0_0 + _S2031.differential_0;
    *&((&(&(&_S1664)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1664)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1664)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2162 = _S2105 + _S2146.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1664)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2163 = _S2104 + _S2146.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1664)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2164 = _S2103 + _S2146.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1664)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2165 = _S2110 + _S2146.vignette_params_0[int(2)].cy_0;
    (&(&_S1664)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2166 = _S2111 + _S2146.vignette_params_0[int(2)].cx_0;
    (&(&_S1664)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2167 = _S2119 + _S2146.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1664)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2168 = _S2118 + _S2146.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1664)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2169 = _S2117 + _S2146.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1664)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2170 = _S2124 + _S2146.vignette_params_0[int(1)].cy_0;
    (&(&_S1664)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2171 = _S2125 + _S2146.vignette_params_0[int(1)].cx_0;
    (&(&_S1664)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2172 = _S2133 + _S2146.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1664)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2173 = _S2132 + _S2146.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1664)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2174 = _S2131 + _S2146.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1664)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2175 = _S2138 + _S2146.vignette_params_0[int(0)].cy_0;
    (&(&_S1664)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2176 = _S2139 + _S2146.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2177;
    _S2177[int(0)] = 0.0f;
    _S2177[int(1)] = 0.0f;
    _S2177[int(2)] = 0.0f;
    _S2177[int(3)] = 0.0f;
    _S2177[int(4)] = 0.0f;
    _S2177[int(5)] = 0.0f;
    _S2177[int(6)] = 0.0f;
    _S2177[int(7)] = 0.0f;
    _S2177[int(8)] = 0.0f;
    _S2177[int(9)] = 0.0f;
    _S2177[int(10)] = 0.0f;
    _S2177[int(11)] = 0.0f;
    _S2177[int(12)] = 0.0f;
    _S2177[int(13)] = 0.0f;
    _S2177[int(14)] = 0.0f;
    _S2177[int(15)] = 0.0f;
    _S2177[int(16)] = 0.0f;
    _S2177[int(17)] = 0.0f;
    _S2177[int(18)] = 0.0f;
    _S2177[int(19)] = 0.0f;
    _S2177[int(20)] = 0.0f;
    _S2177[int(21)] = 0.0f;
    _S2177[int(22)] = 0.0f;
    _S2177[int(23)] = 0.0f;
    _S2177[int(24)] = 0.0f;
    _S2177[int(25)] = 0.0f;
    _S2177[int(26)] = 0.0f;
    _S2177[int(27)] = 0.0f;
    _S2177[int(28)] = 0.0f;
    _S2177[int(29)] = 0.0f;
    _S2177[int(30)] = 0.0f;
    _S2177[int(31)] = 0.0f;
    _S2177[int(32)] = 0.0f;
    _S2177[int(33)] = 0.0f;
    _S2177[int(34)] = 0.0f;
    _S2177[int(35)] = 0.0f;
    _S2177[int(36)] = 0.0f;
    _S2177[int(37)] = 0.0f;
    _S2177[int(38)] = 0.0f;
    _S2177[int(9)] = _S2168;
    _S2177[int(18)] = _S2146.color_params_0.r_0.x;
    _S2177[int(17)] = _S2146.color_params_0.b_0.y;
    _S2177[int(16)] = _S2146.color_params_0.b_0.x;
    _S2177[int(15)] = _S2162;
    _S2177[int(14)] = _S2163;
    _S2177[int(13)] = _S2164;
    _S2177[int(12)] = _S2165;
    _S2177[int(11)] = _S2166;
    _S2177[int(10)] = _S2167;
    _S2177[int(19)] = _S2146.color_params_0.r_0.y;
    _S2177[int(8)] = _S2169;
    _S2177[int(7)] = _S2170;
    _S2177[int(6)] = _S2171;
    _S2177[int(5)] = _S2172;
    _S2177[int(4)] = _S2173;
    _S2177[int(3)] = _S2174;
    _S2177[int(2)] = _S2175;
    _S2177[int(1)] = _S2176;
    _S2177[int(0)] = _S1664.exposure_0;
    _S2177[int(28)] = _S2157;
    _S2177[int(37)] = _S2148;
    _S2177[int(36)] = _S2149;
    _S2177[int(35)] = _S2150;
    _S2177[int(34)] = _S2151;
    _S2177[int(33)] = _S2152;
    _S2177[int(32)] = _S2153;
    _S2177[int(31)] = _S2154;
    _S2177[int(30)] = _S2155;
    _S2177[int(29)] = _S2156;
    _S2177[int(38)] = _S2147;
    _S2177[int(27)] = _S2158;
    _S2177[int(26)] = _S2159;
    _S2177[int(25)] = _S2160;
    _S2177[int(24)] = _S2161;
    _S2177[int(23)] = _S2146.color_params_0.n_0.y;
    _S2177[int(22)] = _S2146.color_params_0.n_0.x;
    _S2177[int(21)] = _S2146.color_params_0.g_0.y;
    _S2177[int(20)] = _S2146.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S2177;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S2143;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2178, float2  _S2179, float2  _S2180, float2  _S2181, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2182, float3  _S2183)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S2178, _S2179, _S2180, _S2181, _S2182, _S2183);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  * params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S2184 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S2184;
    FixedArray<float, 39>  _S2185 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = *params_3;
    (&dp_params_1)->differential_0 = _S2185;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  * params_4, FixedArray<float, 22>  * _S2186)
{
    PPISPParams_0 p_2;
    (&p_2)->exposure_1 = (*params_4)[int(0)];
    (&(&p_2)->vignette_params_1[int(0)])->cx_0 = (*params_4)[int(1)];
    (&(&p_2)->vignette_params_1[int(0)])->cy_0 = (*params_4)[int(2)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha0_0 = (*params_4)[int(3)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha1_0 = (*params_4)[int(4)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha2_0 = (*params_4)[int(5)];
    (&(&p_2)->vignette_params_1[int(1)])->cx_0 = (*params_4)[int(6)];
    (&(&p_2)->vignette_params_1[int(1)])->cy_0 = (*params_4)[int(7)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha0_0 = (*params_4)[int(8)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha1_0 = (*params_4)[int(9)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha2_0 = (*params_4)[int(10)];
    (&(&p_2)->vignette_params_1[int(2)])->cx_0 = (*params_4)[int(11)];
    (&(&p_2)->vignette_params_1[int(2)])->cy_0 = (*params_4)[int(12)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha0_0 = (*params_4)[int(13)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha1_0 = (*params_4)[int(14)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha2_0 = (*params_4)[int(15)];
    *&((&(&(&p_2)->color_params_1)->b_0)->x) = (*params_4)[int(16)];
    *&((&(&(&p_2)->color_params_1)->b_0)->y) = (*params_4)[int(17)];
    *&((&(&(&p_2)->color_params_1)->r_0)->x) = (*params_4)[int(18)];
    *&((&(&(&p_2)->color_params_1)->r_0)->y) = (*params_4)[int(19)];
    *&((&(&(&p_2)->color_params_1)->g_0)->x) = (*params_4)[int(20)];
    *&((&(&(&p_2)->color_params_1)->g_0)->y) = (*params_4)[int(21)];
    *&((&(&(&p_2)->color_params_1)->n_0)->x) = (*params_4)[int(22)];
    *&((&(&(&p_2)->color_params_1)->n_0)->y) = (*params_4)[int(23)];
    (&(&p_2)->crf_params_1[int(0)])->toe_0 = (*params_4)[int(24)];
    (&(&p_2)->crf_params_1[int(0)])->shoulder_0 = (*params_4)[int(25)];
    (&(&p_2)->crf_params_1[int(0)])->gamma_0 = (*params_4)[int(26)];
    (&(&p_2)->crf_params_1[int(0)])->center_0 = (*params_4)[int(27)];
    (&(&p_2)->crf_params_1[int(1)])->toe_0 = (*params_4)[int(28)];
    (&(&p_2)->crf_params_1[int(1)])->shoulder_0 = (*params_4)[int(29)];
    (&(&p_2)->crf_params_1[int(1)])->gamma_0 = (*params_4)[int(30)];
    (&(&p_2)->crf_params_1[int(1)])->center_0 = (*params_4)[int(31)];
    (&(&p_2)->crf_params_1[int(2)])->toe_0 = (*params_4)[int(32)];
    (&(&p_2)->crf_params_1[int(2)])->shoulder_0 = (*params_4)[int(33)];
    (&(&p_2)->crf_params_1[int(2)])->gamma_0 = (*params_4)[int(34)];
    (&(&p_2)->crf_params_1[int(2)])->center_0 = (*params_4)[int(35)];
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
    float _S2187 = p_2.vignette_params_1[int(0)].cx_0;
    float _S2188 = p_2.vignette_params_1[int(0)].cy_0;
    float _S2189 = p_2.vignette_params_1[int(1)].cx_0;
    float _S2190 = p_2.vignette_params_1[int(1)].cy_0;
    float _S2191 = p_2.vignette_params_1[int(2)].cx_0;
    float _S2192 = p_2.vignette_params_1[int(2)].cy_0;
    losses_3[int(1)] = _S2187 * _S2187 + _S2188 * _S2188 + _S2189 * _S2189 + _S2190 * _S2190 + _S2191 * _S2191 + _S2192 * _S2192;
    losses_3[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_3[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_3[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_3 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S2193 = p_2.vignette_params_1[int(0)].cx_0 - mean_3;
    float _S2194 = p_2.vignette_params_1[int(1)].cx_0 - mean_3;
    float _S2195 = p_2.vignette_params_1[int(2)].cx_0 - mean_3;
    losses_3[int(5)] = (_S2193 * _S2193 + _S2194 * _S2194 + _S2195 * _S2195) / 3.0f;
    float mean_4 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S2196 = p_2.vignette_params_1[int(0)].cy_0 - mean_4;
    float _S2197 = p_2.vignette_params_1[int(1)].cy_0 - mean_4;
    float _S2198 = p_2.vignette_params_1[int(2)].cy_0 - mean_4;
    losses_3[int(6)] = (_S2196 * _S2196 + _S2197 * _S2197 + _S2198 * _S2198) / 3.0f;
    float mean_5 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S2199 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_5;
    float _S2200 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_5;
    float _S2201 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_5;
    losses_3[int(7)] = (_S2199 * _S2199 + _S2200 * _S2200 + _S2201 * _S2201) / 3.0f;
    float mean_6 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S2202 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_6;
    float _S2203 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_6;
    float _S2204 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_6;
    losses_3[int(8)] = (_S2202 * _S2202 + _S2203 * _S2203 + _S2204 * _S2204) / 3.0f;
    float mean_7 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S2205 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_7;
    float _S2206 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_7;
    float _S2207 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_7;
    losses_3[int(9)] = (_S2205 * _S2205 + _S2206 * _S2206 + _S2207 * _S2207) / 3.0f;
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
    float _S2208 = p_2.crf_params_1[int(0)].toe_0 - mean_8;
    float _S2209 = p_2.crf_params_1[int(1)].toe_0 - mean_8;
    float _S2210 = p_2.crf_params_1[int(2)].toe_0 - mean_8;
    losses_3[int(18)] = (_S2208 * _S2208 + _S2209 * _S2209 + _S2210 * _S2210) / 3.0f;
    float mean_9 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S2211 = p_2.crf_params_1[int(0)].shoulder_0 - mean_9;
    float _S2212 = p_2.crf_params_1[int(1)].shoulder_0 - mean_9;
    float _S2213 = p_2.crf_params_1[int(2)].shoulder_0 - mean_9;
    losses_3[int(19)] = (_S2211 * _S2211 + _S2212 * _S2212 + _S2213 * _S2213) / 3.0f;
    float mean_10 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S2214 = p_2.crf_params_1[int(0)].gamma_0 - mean_10;
    float _S2215 = p_2.crf_params_1[int(1)].gamma_0 - mean_10;
    float _S2216 = p_2.crf_params_1[int(2)].gamma_0 - mean_10;
    losses_3[int(20)] = (_S2214 * _S2214 + _S2215 * _S2215 + _S2216 * _S2216) / 3.0f;
    float mean_11 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S2217 = p_2.crf_params_1[int(0)].center_0 - mean_11;
    float _S2218 = p_2.crf_params_1[int(1)].center_0 - mean_11;
    float _S2219 = p_2.crf_params_1[int(2)].center_0 - mean_11;
    losses_3[int(21)] = (_S2217 * _S2217 + _S2218 * _S2218 + _S2219 * _S2219) / 3.0f;
    *_S2186 = losses_3;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  * params_5, FixedArray<float, 23>  * _S2220)
{
    PPISPParamsRQS_0 p_3;
    (&p_3)->exposure_0 = (*params_5)[int(0)];
    (&(&p_3)->vignette_params_0[int(0)])->cx_0 = (*params_5)[int(1)];
    (&(&p_3)->vignette_params_0[int(0)])->cy_0 = (*params_5)[int(2)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha0_0 = (*params_5)[int(3)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha1_0 = (*params_5)[int(4)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha2_0 = (*params_5)[int(5)];
    (&(&p_3)->vignette_params_0[int(1)])->cx_0 = (*params_5)[int(6)];
    (&(&p_3)->vignette_params_0[int(1)])->cy_0 = (*params_5)[int(7)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha0_0 = (*params_5)[int(8)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha1_0 = (*params_5)[int(9)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha2_0 = (*params_5)[int(10)];
    (&(&p_3)->vignette_params_0[int(2)])->cx_0 = (*params_5)[int(11)];
    (&(&p_3)->vignette_params_0[int(2)])->cy_0 = (*params_5)[int(12)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha0_0 = (*params_5)[int(13)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha1_0 = (*params_5)[int(14)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha2_0 = (*params_5)[int(15)];
    *&((&(&(&p_3)->color_params_0)->b_0)->x) = (*params_5)[int(16)];
    *&((&(&(&p_3)->color_params_0)->b_0)->y) = (*params_5)[int(17)];
    *&((&(&(&p_3)->color_params_0)->r_0)->x) = (*params_5)[int(18)];
    *&((&(&(&p_3)->color_params_0)->r_0)->y) = (*params_5)[int(19)];
    *&((&(&(&p_3)->color_params_0)->g_0)->x) = (*params_5)[int(20)];
    *&((&(&(&p_3)->color_params_0)->g_0)->y) = (*params_5)[int(21)];
    *&((&(&(&p_3)->color_params_0)->n_0)->x) = (*params_5)[int(22)];
    *&((&(&(&p_3)->color_params_0)->n_0)->y) = (*params_5)[int(23)];
    (&(&p_3)->crf_params_0[int(0)])->g0_0 = (*params_5)[int(24)];
    (&(&p_3)->crf_params_0[int(0)])->g1_0 = (*params_5)[int(25)];
    (&(&p_3)->crf_params_0[int(0)])->x0_0 = (*params_5)[int(26)];
    (&(&p_3)->crf_params_0[int(0)])->y0_0 = (*params_5)[int(27)];
    (&(&p_3)->crf_params_0[int(0)])->gc_0 = (*params_5)[int(28)];
    (&(&p_3)->crf_params_0[int(1)])->g0_0 = (*params_5)[int(29)];
    (&(&p_3)->crf_params_0[int(1)])->g1_0 = (*params_5)[int(30)];
    (&(&p_3)->crf_params_0[int(1)])->x0_0 = (*params_5)[int(31)];
    (&(&p_3)->crf_params_0[int(1)])->y0_0 = (*params_5)[int(32)];
    (&(&p_3)->crf_params_0[int(1)])->gc_0 = (*params_5)[int(33)];
    (&(&p_3)->crf_params_0[int(2)])->g0_0 = (*params_5)[int(34)];
    (&(&p_3)->crf_params_0[int(2)])->g1_0 = (*params_5)[int(35)];
    (&(&p_3)->crf_params_0[int(2)])->x0_0 = (*params_5)[int(36)];
    (&(&p_3)->crf_params_0[int(2)])->y0_0 = (*params_5)[int(37)];
    (&(&p_3)->crf_params_0[int(2)])->gc_0 = (*params_5)[int(38)];
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
    float _S2221 = p_3.vignette_params_0[int(0)].cx_0;
    float _S2222 = p_3.vignette_params_0[int(0)].cy_0;
    float _S2223 = p_3.vignette_params_0[int(1)].cx_0;
    float _S2224 = p_3.vignette_params_0[int(1)].cy_0;
    float _S2225 = p_3.vignette_params_0[int(2)].cx_0;
    float _S2226 = p_3.vignette_params_0[int(2)].cy_0;
    losses_4[int(1)] = _S2221 * _S2221 + _S2222 * _S2222 + _S2223 * _S2223 + _S2224 * _S2224 + _S2225 * _S2225 + _S2226 * _S2226;
    losses_4[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_4[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_4[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_12 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S2227 = p_3.vignette_params_0[int(0)].cx_0 - mean_12;
    float _S2228 = p_3.vignette_params_0[int(1)].cx_0 - mean_12;
    float _S2229 = p_3.vignette_params_0[int(2)].cx_0 - mean_12;
    losses_4[int(5)] = (_S2227 * _S2227 + _S2228 * _S2228 + _S2229 * _S2229) / 3.0f;
    float mean_13 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S2230 = p_3.vignette_params_0[int(0)].cy_0 - mean_13;
    float _S2231 = p_3.vignette_params_0[int(1)].cy_0 - mean_13;
    float _S2232 = p_3.vignette_params_0[int(2)].cy_0 - mean_13;
    losses_4[int(6)] = (_S2230 * _S2230 + _S2231 * _S2231 + _S2232 * _S2232) / 3.0f;
    float mean_14 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S2233 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_14;
    float _S2234 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_14;
    float _S2235 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_14;
    losses_4[int(7)] = (_S2233 * _S2233 + _S2234 * _S2234 + _S2235 * _S2235) / 3.0f;
    float mean_15 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S2236 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_15;
    float _S2237 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_15;
    float _S2238 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_15;
    losses_4[int(8)] = (_S2236 * _S2236 + _S2237 * _S2237 + _S2238 * _S2238) / 3.0f;
    float mean_16 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S2239 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_16;
    float _S2240 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_16;
    float _S2241 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_16;
    losses_4[int(9)] = (_S2239 * _S2239 + _S2240 * _S2240 + _S2241 * _S2241) / 3.0f;
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
    float _S2242 = p_3.crf_params_0[int(0)].g0_0 - mean_17;
    float _S2243 = p_3.crf_params_0[int(1)].g0_0 - mean_17;
    float _S2244 = p_3.crf_params_0[int(2)].g0_0 - mean_17;
    losses_4[int(18)] = (_S2242 * _S2242 + _S2243 * _S2243 + _S2244 * _S2244) / 3.0f;
    float mean_18 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S2245 = p_3.crf_params_0[int(0)].g1_0 - mean_18;
    float _S2246 = p_3.crf_params_0[int(1)].g1_0 - mean_18;
    float _S2247 = p_3.crf_params_0[int(2)].g1_0 - mean_18;
    losses_4[int(19)] = (_S2245 * _S2245 + _S2246 * _S2246 + _S2247 * _S2247) / 3.0f;
    float mean_19 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S2248 = p_3.crf_params_0[int(0)].x0_0 - mean_19;
    float _S2249 = p_3.crf_params_0[int(1)].x0_0 - mean_19;
    float _S2250 = p_3.crf_params_0[int(2)].x0_0 - mean_19;
    losses_4[int(20)] = (_S2248 * _S2248 + _S2249 * _S2249 + _S2250 * _S2250) / 3.0f;
    float mean_20 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S2251 = p_3.crf_params_0[int(0)].y0_0 - mean_20;
    float _S2252 = p_3.crf_params_0[int(1)].y0_0 - mean_20;
    float _S2253 = p_3.crf_params_0[int(2)].y0_0 - mean_20;
    losses_4[int(21)] = (_S2251 * _S2251 + _S2252 * _S2252 + _S2253 * _S2253) / 3.0f;
    float mean_21 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S2254 = p_3.crf_params_0[int(0)].gc_0 - mean_21;
    float _S2255 = p_3.crf_params_0[int(1)].gc_0 - mean_21;
    float _S2256 = p_3.crf_params_0[int(2)].gc_0 - mean_21;
    losses_4[int(22)] = (_S2254 * _S2254 + _S2255 * _S2255 + _S2256 * _S2256) / 3.0f;
    *_S2220 = losses_4;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_13)
{
    VignettingChannelParams_0 _S2257 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2258 = {
        _S2257, _S2257, _S2257
    };
    float2  _S2259 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2260 = { _S2259, _S2259, _S2259, _S2259 };
    CRFPPISPChannelParams_0 _S2261 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S2262 = {
        _S2261, _S2261, _S2261
    };
    PPISPParams_0 _S2263;
    (&_S2263)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S2263)->vignette_params_1 = _S2258;
    (&_S2263)->color_params_1 = _S2260;
    (&_S2263)->crf_params_1 = _S2262;
    (&(&_S2263)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S2263)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S2263)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S2263)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S2263)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S2263)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S2263)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S2263)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S2263)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S2263)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S2263)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S2263)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S2263)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S2263)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S2263)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S2263)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S2263)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S2263)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S2263)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S2263)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S2263)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S2263)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S2263)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S2263)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S2263)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S2263)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S2263)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S2263)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S2263)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S2263)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S2263)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S2263)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S2263)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S2263)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S2263)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_22 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S2264 = dpparams_2->primal_0[int(1)] - mean_22;
    float _S2265 = dpparams_2->primal_0[int(6)] - mean_22;
    float _S2266 = dpparams_2->primal_0[int(11)] - mean_22;
    float mean_23 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S2267 = dpparams_2->primal_0[int(2)] - mean_23;
    float _S2268 = dpparams_2->primal_0[int(7)] - mean_23;
    float _S2269 = dpparams_2->primal_0[int(12)] - mean_23;
    float mean_24 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S2270 = dpparams_2->primal_0[int(3)] - mean_24;
    float _S2271 = dpparams_2->primal_0[int(8)] - mean_24;
    float _S2272 = dpparams_2->primal_0[int(13)] - mean_24;
    float mean_25 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S2273 = dpparams_2->primal_0[int(4)] - mean_25;
    float _S2274 = dpparams_2->primal_0[int(9)] - mean_25;
    float _S2275 = dpparams_2->primal_0[int(14)] - mean_25;
    float mean_26 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S2276 = dpparams_2->primal_0[int(5)] - mean_26;
    float _S2277 = dpparams_2->primal_0[int(10)] - mean_26;
    float _S2278 = dpparams_2->primal_0[int(15)] - mean_26;
    float mean_27 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_28 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_29 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_30 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S2279 = 0.3333333432674408f * (*_s_dOut_13)[int(21)];
    float _S2280 = (dpparams_2->primal_0[int(35)] - mean_30) * _S2279;
    float _S2281 = _S2280 + _S2280;
    float _S2282 = (dpparams_2->primal_0[int(31)] - mean_30) * _S2279;
    float _S2283 = _S2282 + _S2282;
    float _S2284 = (dpparams_2->primal_0[int(27)] - mean_30) * _S2279;
    float _S2285 = _S2284 + _S2284;
    float _S2286 = 0.3333333432674408f * (- _S2281 + - _S2283 + - _S2285);
    float _S2287 = 0.3333333432674408f * (*_s_dOut_13)[int(20)];
    float _S2288 = (dpparams_2->primal_0[int(34)] - mean_29) * _S2287;
    float _S2289 = _S2288 + _S2288;
    float _S2290 = (dpparams_2->primal_0[int(30)] - mean_29) * _S2287;
    float _S2291 = _S2290 + _S2290;
    float _S2292 = (dpparams_2->primal_0[int(26)] - mean_29) * _S2287;
    float _S2293 = _S2292 + _S2292;
    float _S2294 = 0.3333333432674408f * (- _S2289 + - _S2291 + - _S2293);
    float _S2295 = 0.3333333432674408f * (*_s_dOut_13)[int(19)];
    float _S2296 = (dpparams_2->primal_0[int(33)] - mean_28) * _S2295;
    float _S2297 = _S2296 + _S2296;
    float _S2298 = (dpparams_2->primal_0[int(29)] - mean_28) * _S2295;
    float _S2299 = _S2298 + _S2298;
    float _S2300 = (dpparams_2->primal_0[int(25)] - mean_28) * _S2295;
    float _S2301 = _S2300 + _S2300;
    float _S2302 = 0.3333333432674408f * (- _S2297 + - _S2299 + - _S2301);
    float _S2303 = 0.3333333432674408f * (*_s_dOut_13)[int(18)];
    float _S2304 = (dpparams_2->primal_0[int(32)] - mean_27) * _S2303;
    float _S2305 = _S2304 + _S2304;
    float _S2306 = (dpparams_2->primal_0[int(28)] - mean_27) * _S2303;
    float _S2307 = _S2306 + _S2306;
    float _S2308 = (dpparams_2->primal_0[int(24)] - mean_27) * _S2303;
    float _S2309 = _S2308 + _S2308;
    float _S2310 = 0.3333333432674408f * (- _S2305 + - _S2307 + - _S2309);
    float2  _S2311 = make_float2 ((*_s_dOut_13)[int(16)], (*_s_dOut_13)[int(17)]);
    Matrix<float, 2, 2>  _S2312 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2313;
    (&_S2313)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2313)->differential_0 = _S2312;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2314;
    (&_S2314)->primal_0 = _S2263.color_params_1.n_0;
    (&_S2314)->differential_0 = _S2259;
    s_bwd_prop_mul_3(&_S2313, &_S2314, _S2311);
    float2  _S2315 = make_float2 ((*_s_dOut_13)[int(14)], (*_s_dOut_13)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2316;
    (&_S2316)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2316)->differential_0 = _S2312;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2317;
    (&_S2317)->primal_0 = _S2263.color_params_1.g_0;
    (&_S2317)->differential_0 = _S2259;
    s_bwd_prop_mul_3(&_S2316, &_S2317, _S2315);
    float2  _S2318 = make_float2 ((*_s_dOut_13)[int(12)], (*_s_dOut_13)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2319;
    (&_S2319)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2319)->differential_0 = _S2312;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2320;
    (&_S2320)->primal_0 = _S2263.color_params_1.r_0;
    (&_S2320)->differential_0 = _S2259;
    s_bwd_prop_mul_3(&_S2319, &_S2320, _S2318);
    float2  _S2321 = make_float2 ((*_s_dOut_13)[int(10)], (*_s_dOut_13)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2322;
    (&_S2322)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2322)->differential_0 = _S2312;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2323;
    (&_S2323)->primal_0 = _S2263.color_params_1.b_0;
    (&_S2323)->differential_0 = _S2259;
    s_bwd_prop_mul_3(&_S2322, &_S2323, _S2321);
    ColorPPISPParams_0 _S2324 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2324)->n_0 = _S2314.differential_0;
    (&_S2324)->g_0 = _S2317.differential_0;
    (&_S2324)->r_0 = _S2320.differential_0;
    (&_S2324)->b_0 = _S2323.differential_0;
    float _S2325 = 0.3333333432674408f * (*_s_dOut_13)[int(9)];
    float _S2326 = _S2278 * _S2325;
    float _S2327 = _S2326 + _S2326;
    float _S2328 = _S2277 * _S2325;
    float _S2329 = _S2328 + _S2328;
    float _S2330 = _S2276 * _S2325;
    float _S2331 = _S2330 + _S2330;
    float _S2332 = 0.3333333432674408f * (- _S2327 + - _S2329 + - _S2331);
    float _S2333 = 0.3333333432674408f * (*_s_dOut_13)[int(8)];
    float _S2334 = _S2275 * _S2333;
    float _S2335 = _S2334 + _S2334;
    float _S2336 = _S2274 * _S2333;
    float _S2337 = _S2336 + _S2336;
    float _S2338 = _S2273 * _S2333;
    float _S2339 = _S2338 + _S2338;
    float _S2340 = 0.3333333432674408f * (- _S2335 + - _S2337 + - _S2339);
    float _S2341 = 0.3333333432674408f * (*_s_dOut_13)[int(7)];
    float _S2342 = _S2272 * _S2341;
    float _S2343 = _S2342 + _S2342;
    float _S2344 = _S2271 * _S2341;
    float _S2345 = _S2344 + _S2344;
    float _S2346 = _S2270 * _S2341;
    float _S2347 = _S2346 + _S2346;
    float _S2348 = 0.3333333432674408f * (- _S2343 + - _S2345 + - _S2347);
    float _S2349 = 0.3333333432674408f * (*_s_dOut_13)[int(6)];
    float _S2350 = _S2269 * _S2349;
    float _S2351 = _S2350 + _S2350;
    float _S2352 = _S2268 * _S2349;
    float _S2353 = _S2352 + _S2352;
    float _S2354 = _S2267 * _S2349;
    float _S2355 = _S2354 + _S2354;
    float _S2356 = 0.3333333432674408f * (- _S2351 + - _S2353 + - _S2355);
    float _S2357 = 0.3333333432674408f * (*_s_dOut_13)[int(5)];
    float _S2358 = _S2266 * _S2357;
    float _S2359 = _S2358 + _S2358;
    float _S2360 = _S2265 * _S2357;
    float _S2361 = _S2360 + _S2360;
    float _S2362 = _S2264 * _S2357;
    float _S2363 = _S2362 + _S2362;
    float _S2364 = 0.3333333432674408f * (- _S2359 + - _S2361 + - _S2363);
    DiffPair_float_0 _S2365;
    (&_S2365)->primal_0 = 0.0f;
    (&_S2365)->differential_0 = 0.0f;
    DiffPair_float_0 _S2366;
    (&_S2366)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S2366)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2365, &_S2366, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2367;
    (&_S2367)->primal_0 = 0.0f;
    (&_S2367)->differential_0 = 0.0f;
    DiffPair_float_0 _S2368;
    (&_S2368)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S2368)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2367, &_S2368, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2369;
    (&_S2369)->primal_0 = 0.0f;
    (&_S2369)->differential_0 = 0.0f;
    DiffPair_float_0 _S2370;
    (&_S2370)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S2370)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2369, &_S2370, (*_s_dOut_13)[int(4)]);
    DiffPair_float_0 _S2371;
    (&_S2371)->primal_0 = 0.0f;
    (&_S2371)->differential_0 = 0.0f;
    DiffPair_float_0 _S2372;
    (&_S2372)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S2372)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2371, &_S2372, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2373;
    (&_S2373)->primal_0 = 0.0f;
    (&_S2373)->differential_0 = 0.0f;
    DiffPair_float_0 _S2374;
    (&_S2374)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S2374)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2373, &_S2374, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2375;
    (&_S2375)->primal_0 = 0.0f;
    (&_S2375)->differential_0 = 0.0f;
    DiffPair_float_0 _S2376;
    (&_S2376)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S2376)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2375, &_S2376, (*_s_dOut_13)[int(3)]);
    DiffPair_float_0 _S2377;
    (&_S2377)->primal_0 = 0.0f;
    (&_S2377)->differential_0 = 0.0f;
    DiffPair_float_0 _S2378;
    (&_S2378)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S2378)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2377, &_S2378, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2379;
    (&_S2379)->primal_0 = 0.0f;
    (&_S2379)->differential_0 = 0.0f;
    DiffPair_float_0 _S2380;
    (&_S2380)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S2380)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2379, &_S2380, (*_s_dOut_13)[int(2)]);
    DiffPair_float_0 _S2381;
    (&_S2381)->primal_0 = 0.0f;
    (&_S2381)->differential_0 = 0.0f;
    DiffPair_float_0 _S2382;
    (&_S2382)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S2382)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2381, &_S2382, (*_s_dOut_13)[int(2)]);
    float _S2383 = dpparams_2->primal_0[int(12)] * (*_s_dOut_13)[int(1)];
    float _S2384 = dpparams_2->primal_0[int(11)] * (*_s_dOut_13)[int(1)];
    float _S2385 = dpparams_2->primal_0[int(7)] * (*_s_dOut_13)[int(1)];
    float _S2386 = dpparams_2->primal_0[int(6)] * (*_s_dOut_13)[int(1)];
    float _S2387 = dpparams_2->primal_0[int(2)] * (*_s_dOut_13)[int(1)];
    float _S2388 = dpparams_2->primal_0[int(1)] * (*_s_dOut_13)[int(1)];
    PPISPParams_0 _S2389 = PPISPParams_x24_syn_dzero_0();
    (&_S2389)->color_params_1 = _S2324;
    (&_S2389)->exposure_1 = (*_s_dOut_13)[int(0)];
    _S2263 = _S2389;
    (&(&_S2263)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2390 = _S2281 + _S2286 + _S2389.crf_params_1[int(2)].center_0;
    (&(&_S2263)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2391 = _S2289 + _S2294 + _S2389.crf_params_1[int(2)].gamma_0;
    (&(&_S2263)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2392 = _S2297 + _S2302 + _S2389.crf_params_1[int(2)].shoulder_0;
    (&(&_S2263)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2393 = _S2305 + _S2310 + _S2389.crf_params_1[int(2)].toe_0;
    (&(&_S2263)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2394 = _S2283 + _S2286 + _S2389.crf_params_1[int(1)].center_0;
    (&(&_S2263)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2395 = _S2291 + _S2294 + _S2389.crf_params_1[int(1)].gamma_0;
    (&(&_S2263)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2396 = _S2299 + _S2302 + _S2389.crf_params_1[int(1)].shoulder_0;
    (&(&_S2263)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2397 = _S2307 + _S2310 + _S2389.crf_params_1[int(1)].toe_0;
    (&(&_S2263)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2398 = _S2285 + _S2286 + _S2389.crf_params_1[int(0)].center_0;
    (&(&_S2263)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2399 = _S2293 + _S2294 + _S2389.crf_params_1[int(0)].gamma_0;
    (&(&_S2263)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2400 = _S2301 + _S2302 + _S2389.crf_params_1[int(0)].shoulder_0;
    (&(&_S2263)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2401 = _S2309 + _S2310 + _S2389.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S2263)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S2263)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S2263)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2402 = _S2327 + _S2332 + _S2366.differential_0 + _S2389.vignette_params_1[int(2)].alpha2_0;
    (&(&_S2263)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2403 = _S2335 + _S2340 + _S2372.differential_0 + _S2389.vignette_params_1[int(2)].alpha1_0;
    (&(&_S2263)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2404 = _S2343 + _S2348 + _S2378.differential_0 + _S2389.vignette_params_1[int(2)].alpha0_0;
    (&(&_S2263)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2405 = _S2351 + _S2356 + _S2383 + _S2383 + _S2389.vignette_params_1[int(2)].cy_0;
    (&(&_S2263)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2406 = _S2359 + _S2364 + _S2384 + _S2384 + _S2389.vignette_params_1[int(2)].cx_0;
    (&(&_S2263)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2407 = _S2329 + _S2332 + _S2368.differential_0 + _S2389.vignette_params_1[int(1)].alpha2_0;
    (&(&_S2263)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2408 = _S2337 + _S2340 + _S2374.differential_0 + _S2389.vignette_params_1[int(1)].alpha1_0;
    (&(&_S2263)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2409 = _S2345 + _S2348 + _S2380.differential_0 + _S2389.vignette_params_1[int(1)].alpha0_0;
    (&(&_S2263)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2410 = _S2353 + _S2356 + _S2385 + _S2385 + _S2389.vignette_params_1[int(1)].cy_0;
    (&(&_S2263)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2411 = _S2361 + _S2364 + _S2386 + _S2386 + _S2389.vignette_params_1[int(1)].cx_0;
    (&(&_S2263)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2412 = _S2331 + _S2332 + _S2370.differential_0 + _S2389.vignette_params_1[int(0)].alpha2_0;
    (&(&_S2263)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2413 = _S2339 + _S2340 + _S2376.differential_0 + _S2389.vignette_params_1[int(0)].alpha1_0;
    (&(&_S2263)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2414 = _S2347 + _S2348 + _S2382.differential_0 + _S2389.vignette_params_1[int(0)].alpha0_0;
    (&(&_S2263)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2415 = _S2355 + _S2356 + _S2387 + _S2387 + _S2389.vignette_params_1[int(0)].cy_0;
    (&(&_S2263)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2416 = _S2363 + _S2364 + _S2388 + _S2388 + _S2389.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2417;
    _S2417[int(0)] = 0.0f;
    _S2417[int(1)] = 0.0f;
    _S2417[int(2)] = 0.0f;
    _S2417[int(3)] = 0.0f;
    _S2417[int(4)] = 0.0f;
    _S2417[int(5)] = 0.0f;
    _S2417[int(6)] = 0.0f;
    _S2417[int(7)] = 0.0f;
    _S2417[int(8)] = 0.0f;
    _S2417[int(9)] = 0.0f;
    _S2417[int(10)] = 0.0f;
    _S2417[int(11)] = 0.0f;
    _S2417[int(12)] = 0.0f;
    _S2417[int(13)] = 0.0f;
    _S2417[int(14)] = 0.0f;
    _S2417[int(15)] = 0.0f;
    _S2417[int(16)] = 0.0f;
    _S2417[int(17)] = 0.0f;
    _S2417[int(18)] = 0.0f;
    _S2417[int(19)] = 0.0f;
    _S2417[int(20)] = 0.0f;
    _S2417[int(21)] = 0.0f;
    _S2417[int(22)] = 0.0f;
    _S2417[int(23)] = 0.0f;
    _S2417[int(24)] = 0.0f;
    _S2417[int(25)] = 0.0f;
    _S2417[int(26)] = 0.0f;
    _S2417[int(27)] = 0.0f;
    _S2417[int(28)] = 0.0f;
    _S2417[int(29)] = 0.0f;
    _S2417[int(30)] = 0.0f;
    _S2417[int(31)] = 0.0f;
    _S2417[int(32)] = 0.0f;
    _S2417[int(33)] = 0.0f;
    _S2417[int(34)] = 0.0f;
    _S2417[int(35)] = 0.0f;
    _S2417[int(8)] = _S2409;
    _S2417[int(16)] = _S2389.color_params_1.b_0.x;
    _S2417[int(15)] = _S2402;
    _S2417[int(14)] = _S2403;
    _S2417[int(13)] = _S2404;
    _S2417[int(12)] = _S2405;
    _S2417[int(11)] = _S2406;
    _S2417[int(10)] = _S2407;
    _S2417[int(9)] = _S2408;
    _S2417[int(17)] = _S2389.color_params_1.b_0.y;
    _S2417[int(7)] = _S2410;
    _S2417[int(6)] = _S2411;
    _S2417[int(5)] = _S2412;
    _S2417[int(4)] = _S2413;
    _S2417[int(3)] = _S2414;
    _S2417[int(2)] = _S2415;
    _S2417[int(1)] = _S2416;
    _S2417[int(0)] = _S2263.exposure_1;
    _S2417[int(26)] = _S2399;
    _S2417[int(34)] = _S2391;
    _S2417[int(33)] = _S2392;
    _S2417[int(32)] = _S2393;
    _S2417[int(31)] = _S2394;
    _S2417[int(30)] = _S2395;
    _S2417[int(29)] = _S2396;
    _S2417[int(28)] = _S2397;
    _S2417[int(27)] = _S2398;
    _S2417[int(35)] = _S2390;
    _S2417[int(25)] = _S2400;
    _S2417[int(24)] = _S2401;
    _S2417[int(23)] = _S2389.color_params_1.n_0.y;
    _S2417[int(22)] = _S2389.color_params_1.n_0.x;
    _S2417[int(21)] = _S2389.color_params_1.g_0.y;
    _S2417[int(20)] = _S2389.color_params_1.g_0.x;
    _S2417[int(19)] = _S2389.color_params_1.r_0.y;
    _S2417[int(18)] = _S2389.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S2417;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2418, FixedArray<float, 22>  * _S2419)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S2418, _S2419);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  * params_6, FixedArray<float, 22>  * grad_out_2, FixedArray<float, 36>  * _S2420)
{
    FixedArray<float, 36>  _S2421 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = *params_6;
    (&dp_params_2)->differential_0 = _S2421;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, grad_out_2);
    *_S2420 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_14)
{
    VignettingChannelParams_0 _S2422 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2423 = {
        _S2422, _S2422, _S2422
    };
    float2  _S2424 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2425 = { _S2424, _S2424, _S2424, _S2424 };
    RQSCRFPPISPChannelParams_0 _S2426 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2427 = {
        _S2426, _S2426, _S2426
    };
    PPISPParamsRQS_0 _S2428;
    (&_S2428)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S2428)->vignette_params_0 = _S2423;
    (&_S2428)->color_params_0 = _S2425;
    (&_S2428)->crf_params_0 = _S2427;
    (&(&_S2428)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S2428)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S2428)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S2428)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S2428)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S2428)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S2428)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S2428)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S2428)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S2428)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S2428)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S2428)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S2428)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S2428)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S2428)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S2428)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S2428)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S2428)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S2428)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S2428)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S2428)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S2428)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S2428)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S2428)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S2428)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S2428)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S2428)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S2428)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S2428)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S2428)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S2428)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S2428)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S2428)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S2428)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S2428)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S2428)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S2428)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S2428)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_31 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S2429 = dpparams_3->primal_0[int(1)] - mean_31;
    float _S2430 = dpparams_3->primal_0[int(6)] - mean_31;
    float _S2431 = dpparams_3->primal_0[int(11)] - mean_31;
    float mean_32 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S2432 = dpparams_3->primal_0[int(2)] - mean_32;
    float _S2433 = dpparams_3->primal_0[int(7)] - mean_32;
    float _S2434 = dpparams_3->primal_0[int(12)] - mean_32;
    float mean_33 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S2435 = dpparams_3->primal_0[int(3)] - mean_33;
    float _S2436 = dpparams_3->primal_0[int(8)] - mean_33;
    float _S2437 = dpparams_3->primal_0[int(13)] - mean_33;
    float mean_34 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S2438 = dpparams_3->primal_0[int(4)] - mean_34;
    float _S2439 = dpparams_3->primal_0[int(9)] - mean_34;
    float _S2440 = dpparams_3->primal_0[int(14)] - mean_34;
    float mean_35 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S2441 = dpparams_3->primal_0[int(5)] - mean_35;
    float _S2442 = dpparams_3->primal_0[int(10)] - mean_35;
    float _S2443 = dpparams_3->primal_0[int(15)] - mean_35;
    float mean_36 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_37 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_38 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_39 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_40 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S2444 = 0.3333333432674408f * (*_s_dOut_14)[int(22)];
    float _S2445 = (dpparams_3->primal_0[int(38)] - mean_40) * _S2444;
    float _S2446 = _S2445 + _S2445;
    float _S2447 = (dpparams_3->primal_0[int(33)] - mean_40) * _S2444;
    float _S2448 = _S2447 + _S2447;
    float _S2449 = (dpparams_3->primal_0[int(28)] - mean_40) * _S2444;
    float _S2450 = _S2449 + _S2449;
    float _S2451 = 0.3333333432674408f * (- _S2446 + - _S2448 + - _S2450);
    float _S2452 = 0.3333333432674408f * (*_s_dOut_14)[int(21)];
    float _S2453 = (dpparams_3->primal_0[int(37)] - mean_39) * _S2452;
    float _S2454 = _S2453 + _S2453;
    float _S2455 = (dpparams_3->primal_0[int(32)] - mean_39) * _S2452;
    float _S2456 = _S2455 + _S2455;
    float _S2457 = (dpparams_3->primal_0[int(27)] - mean_39) * _S2452;
    float _S2458 = _S2457 + _S2457;
    float _S2459 = 0.3333333432674408f * (- _S2454 + - _S2456 + - _S2458);
    float _S2460 = 0.3333333432674408f * (*_s_dOut_14)[int(20)];
    float _S2461 = (dpparams_3->primal_0[int(36)] - mean_38) * _S2460;
    float _S2462 = _S2461 + _S2461;
    float _S2463 = (dpparams_3->primal_0[int(31)] - mean_38) * _S2460;
    float _S2464 = _S2463 + _S2463;
    float _S2465 = (dpparams_3->primal_0[int(26)] - mean_38) * _S2460;
    float _S2466 = _S2465 + _S2465;
    float _S2467 = 0.3333333432674408f * (- _S2462 + - _S2464 + - _S2466);
    float _S2468 = 0.3333333432674408f * (*_s_dOut_14)[int(19)];
    float _S2469 = (dpparams_3->primal_0[int(35)] - mean_37) * _S2468;
    float _S2470 = _S2469 + _S2469;
    float _S2471 = (dpparams_3->primal_0[int(30)] - mean_37) * _S2468;
    float _S2472 = _S2471 + _S2471;
    float _S2473 = (dpparams_3->primal_0[int(25)] - mean_37) * _S2468;
    float _S2474 = _S2473 + _S2473;
    float _S2475 = 0.3333333432674408f * (- _S2470 + - _S2472 + - _S2474);
    float _S2476 = 0.3333333432674408f * (*_s_dOut_14)[int(18)];
    float _S2477 = (dpparams_3->primal_0[int(34)] - mean_36) * _S2476;
    float _S2478 = _S2477 + _S2477;
    float _S2479 = (dpparams_3->primal_0[int(29)] - mean_36) * _S2476;
    float _S2480 = _S2479 + _S2479;
    float _S2481 = (dpparams_3->primal_0[int(24)] - mean_36) * _S2476;
    float _S2482 = _S2481 + _S2481;
    float _S2483 = 0.3333333432674408f * (- _S2478 + - _S2480 + - _S2482);
    float2  _S2484 = make_float2 ((*_s_dOut_14)[int(16)], (*_s_dOut_14)[int(17)]);
    Matrix<float, 2, 2>  _S2485 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2486;
    (&_S2486)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2486)->differential_0 = _S2485;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2487;
    (&_S2487)->primal_0 = _S2428.color_params_0.n_0;
    (&_S2487)->differential_0 = _S2424;
    s_bwd_prop_mul_3(&_S2486, &_S2487, _S2484);
    float2  _S2488 = make_float2 ((*_s_dOut_14)[int(14)], (*_s_dOut_14)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2489;
    (&_S2489)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2489)->differential_0 = _S2485;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2490;
    (&_S2490)->primal_0 = _S2428.color_params_0.g_0;
    (&_S2490)->differential_0 = _S2424;
    s_bwd_prop_mul_3(&_S2489, &_S2490, _S2488);
    float2  _S2491 = make_float2 ((*_s_dOut_14)[int(12)], (*_s_dOut_14)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2492;
    (&_S2492)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2492)->differential_0 = _S2485;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2493;
    (&_S2493)->primal_0 = _S2428.color_params_0.r_0;
    (&_S2493)->differential_0 = _S2424;
    s_bwd_prop_mul_3(&_S2492, &_S2493, _S2491);
    float2  _S2494 = make_float2 ((*_s_dOut_14)[int(10)], (*_s_dOut_14)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2495;
    (&_S2495)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2495)->differential_0 = _S2485;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2496;
    (&_S2496)->primal_0 = _S2428.color_params_0.b_0;
    (&_S2496)->differential_0 = _S2424;
    s_bwd_prop_mul_3(&_S2495, &_S2496, _S2494);
    ColorPPISPParams_0 _S2497 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2497)->n_0 = _S2487.differential_0;
    (&_S2497)->g_0 = _S2490.differential_0;
    (&_S2497)->r_0 = _S2493.differential_0;
    (&_S2497)->b_0 = _S2496.differential_0;
    float _S2498 = 0.3333333432674408f * (*_s_dOut_14)[int(9)];
    float _S2499 = _S2443 * _S2498;
    float _S2500 = _S2499 + _S2499;
    float _S2501 = _S2442 * _S2498;
    float _S2502 = _S2501 + _S2501;
    float _S2503 = _S2441 * _S2498;
    float _S2504 = _S2503 + _S2503;
    float _S2505 = 0.3333333432674408f * (- _S2500 + - _S2502 + - _S2504);
    float _S2506 = 0.3333333432674408f * (*_s_dOut_14)[int(8)];
    float _S2507 = _S2440 * _S2506;
    float _S2508 = _S2507 + _S2507;
    float _S2509 = _S2439 * _S2506;
    float _S2510 = _S2509 + _S2509;
    float _S2511 = _S2438 * _S2506;
    float _S2512 = _S2511 + _S2511;
    float _S2513 = 0.3333333432674408f * (- _S2508 + - _S2510 + - _S2512);
    float _S2514 = 0.3333333432674408f * (*_s_dOut_14)[int(7)];
    float _S2515 = _S2437 * _S2514;
    float _S2516 = _S2515 + _S2515;
    float _S2517 = _S2436 * _S2514;
    float _S2518 = _S2517 + _S2517;
    float _S2519 = _S2435 * _S2514;
    float _S2520 = _S2519 + _S2519;
    float _S2521 = 0.3333333432674408f * (- _S2516 + - _S2518 + - _S2520);
    float _S2522 = 0.3333333432674408f * (*_s_dOut_14)[int(6)];
    float _S2523 = _S2434 * _S2522;
    float _S2524 = _S2523 + _S2523;
    float _S2525 = _S2433 * _S2522;
    float _S2526 = _S2525 + _S2525;
    float _S2527 = _S2432 * _S2522;
    float _S2528 = _S2527 + _S2527;
    float _S2529 = 0.3333333432674408f * (- _S2524 + - _S2526 + - _S2528);
    float _S2530 = 0.3333333432674408f * (*_s_dOut_14)[int(5)];
    float _S2531 = _S2431 * _S2530;
    float _S2532 = _S2531 + _S2531;
    float _S2533 = _S2430 * _S2530;
    float _S2534 = _S2533 + _S2533;
    float _S2535 = _S2429 * _S2530;
    float _S2536 = _S2535 + _S2535;
    float _S2537 = 0.3333333432674408f * (- _S2532 + - _S2534 + - _S2536);
    DiffPair_float_0 _S2538;
    (&_S2538)->primal_0 = 0.0f;
    (&_S2538)->differential_0 = 0.0f;
    DiffPair_float_0 _S2539;
    (&_S2539)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S2539)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2538, &_S2539, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2540;
    (&_S2540)->primal_0 = 0.0f;
    (&_S2540)->differential_0 = 0.0f;
    DiffPair_float_0 _S2541;
    (&_S2541)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S2541)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2540, &_S2541, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2542;
    (&_S2542)->primal_0 = 0.0f;
    (&_S2542)->differential_0 = 0.0f;
    DiffPair_float_0 _S2543;
    (&_S2543)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S2543)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2542, &_S2543, (*_s_dOut_14)[int(4)]);
    DiffPair_float_0 _S2544;
    (&_S2544)->primal_0 = 0.0f;
    (&_S2544)->differential_0 = 0.0f;
    DiffPair_float_0 _S2545;
    (&_S2545)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S2545)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2544, &_S2545, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2546;
    (&_S2546)->primal_0 = 0.0f;
    (&_S2546)->differential_0 = 0.0f;
    DiffPair_float_0 _S2547;
    (&_S2547)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S2547)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2546, &_S2547, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2548;
    (&_S2548)->primal_0 = 0.0f;
    (&_S2548)->differential_0 = 0.0f;
    DiffPair_float_0 _S2549;
    (&_S2549)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S2549)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2548, &_S2549, (*_s_dOut_14)[int(3)]);
    DiffPair_float_0 _S2550;
    (&_S2550)->primal_0 = 0.0f;
    (&_S2550)->differential_0 = 0.0f;
    DiffPair_float_0 _S2551;
    (&_S2551)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S2551)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2550, &_S2551, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2552;
    (&_S2552)->primal_0 = 0.0f;
    (&_S2552)->differential_0 = 0.0f;
    DiffPair_float_0 _S2553;
    (&_S2553)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S2553)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2552, &_S2553, (*_s_dOut_14)[int(2)]);
    DiffPair_float_0 _S2554;
    (&_S2554)->primal_0 = 0.0f;
    (&_S2554)->differential_0 = 0.0f;
    DiffPair_float_0 _S2555;
    (&_S2555)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S2555)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2554, &_S2555, (*_s_dOut_14)[int(2)]);
    float _S2556 = dpparams_3->primal_0[int(12)] * (*_s_dOut_14)[int(1)];
    float _S2557 = dpparams_3->primal_0[int(11)] * (*_s_dOut_14)[int(1)];
    float _S2558 = dpparams_3->primal_0[int(7)] * (*_s_dOut_14)[int(1)];
    float _S2559 = dpparams_3->primal_0[int(6)] * (*_s_dOut_14)[int(1)];
    float _S2560 = dpparams_3->primal_0[int(2)] * (*_s_dOut_14)[int(1)];
    float _S2561 = dpparams_3->primal_0[int(1)] * (*_s_dOut_14)[int(1)];
    PPISPParamsRQS_0 _S2562 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2562)->color_params_0 = _S2497;
    (&_S2562)->exposure_0 = (*_s_dOut_14)[int(0)];
    _S2428 = _S2562;
    (&(&_S2428)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2563 = _S2446 + _S2451 + _S2562.crf_params_0[int(2)].gc_0;
    (&(&_S2428)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2564 = _S2454 + _S2459 + _S2562.crf_params_0[int(2)].y0_0;
    (&(&_S2428)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2565 = _S2462 + _S2467 + _S2562.crf_params_0[int(2)].x0_0;
    (&(&_S2428)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2566 = _S2470 + _S2475 + _S2562.crf_params_0[int(2)].g1_0;
    (&(&_S2428)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2567 = _S2478 + _S2483 + _S2562.crf_params_0[int(2)].g0_0;
    (&(&_S2428)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2568 = _S2448 + _S2451 + _S2562.crf_params_0[int(1)].gc_0;
    (&(&_S2428)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2569 = _S2456 + _S2459 + _S2562.crf_params_0[int(1)].y0_0;
    (&(&_S2428)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2570 = _S2464 + _S2467 + _S2562.crf_params_0[int(1)].x0_0;
    (&(&_S2428)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2571 = _S2472 + _S2475 + _S2562.crf_params_0[int(1)].g1_0;
    (&(&_S2428)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2572 = _S2480 + _S2483 + _S2562.crf_params_0[int(1)].g0_0;
    (&(&_S2428)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2573 = _S2450 + _S2451 + _S2562.crf_params_0[int(0)].gc_0;
    (&(&_S2428)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2574 = _S2458 + _S2459 + _S2562.crf_params_0[int(0)].y0_0;
    (&(&_S2428)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2575 = _S2466 + _S2467 + _S2562.crf_params_0[int(0)].x0_0;
    (&(&_S2428)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2576 = _S2474 + _S2475 + _S2562.crf_params_0[int(0)].g1_0;
    (&(&_S2428)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2577 = _S2482 + _S2483 + _S2562.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S2428)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2428)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2428)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2578 = _S2500 + _S2505 + _S2539.differential_0 + _S2562.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2428)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2579 = _S2508 + _S2513 + _S2545.differential_0 + _S2562.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2428)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2580 = _S2516 + _S2521 + _S2551.differential_0 + _S2562.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2428)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2581 = _S2524 + _S2529 + _S2556 + _S2556 + _S2562.vignette_params_0[int(2)].cy_0;
    (&(&_S2428)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2582 = _S2532 + _S2537 + _S2557 + _S2557 + _S2562.vignette_params_0[int(2)].cx_0;
    (&(&_S2428)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2583 = _S2502 + _S2505 + _S2541.differential_0 + _S2562.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2428)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2584 = _S2510 + _S2513 + _S2547.differential_0 + _S2562.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2428)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2585 = _S2518 + _S2521 + _S2553.differential_0 + _S2562.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2428)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2586 = _S2526 + _S2529 + _S2558 + _S2558 + _S2562.vignette_params_0[int(1)].cy_0;
    (&(&_S2428)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2587 = _S2534 + _S2537 + _S2559 + _S2559 + _S2562.vignette_params_0[int(1)].cx_0;
    (&(&_S2428)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2588 = _S2504 + _S2505 + _S2543.differential_0 + _S2562.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2428)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2589 = _S2512 + _S2513 + _S2549.differential_0 + _S2562.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2428)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2590 = _S2520 + _S2521 + _S2555.differential_0 + _S2562.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2428)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2591 = _S2528 + _S2529 + _S2560 + _S2560 + _S2562.vignette_params_0[int(0)].cy_0;
    (&(&_S2428)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2592 = _S2536 + _S2537 + _S2561 + _S2561 + _S2562.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2593;
    _S2593[int(0)] = 0.0f;
    _S2593[int(1)] = 0.0f;
    _S2593[int(2)] = 0.0f;
    _S2593[int(3)] = 0.0f;
    _S2593[int(4)] = 0.0f;
    _S2593[int(5)] = 0.0f;
    _S2593[int(6)] = 0.0f;
    _S2593[int(7)] = 0.0f;
    _S2593[int(8)] = 0.0f;
    _S2593[int(9)] = 0.0f;
    _S2593[int(10)] = 0.0f;
    _S2593[int(11)] = 0.0f;
    _S2593[int(12)] = 0.0f;
    _S2593[int(13)] = 0.0f;
    _S2593[int(14)] = 0.0f;
    _S2593[int(15)] = 0.0f;
    _S2593[int(16)] = 0.0f;
    _S2593[int(17)] = 0.0f;
    _S2593[int(18)] = 0.0f;
    _S2593[int(19)] = 0.0f;
    _S2593[int(20)] = 0.0f;
    _S2593[int(21)] = 0.0f;
    _S2593[int(22)] = 0.0f;
    _S2593[int(23)] = 0.0f;
    _S2593[int(24)] = 0.0f;
    _S2593[int(25)] = 0.0f;
    _S2593[int(26)] = 0.0f;
    _S2593[int(27)] = 0.0f;
    _S2593[int(28)] = 0.0f;
    _S2593[int(29)] = 0.0f;
    _S2593[int(30)] = 0.0f;
    _S2593[int(31)] = 0.0f;
    _S2593[int(32)] = 0.0f;
    _S2593[int(33)] = 0.0f;
    _S2593[int(34)] = 0.0f;
    _S2593[int(35)] = 0.0f;
    _S2593[int(36)] = 0.0f;
    _S2593[int(37)] = 0.0f;
    _S2593[int(38)] = 0.0f;
    _S2593[int(9)] = _S2584;
    _S2593[int(18)] = _S2562.color_params_0.r_0.x;
    _S2593[int(17)] = _S2562.color_params_0.b_0.y;
    _S2593[int(16)] = _S2562.color_params_0.b_0.x;
    _S2593[int(15)] = _S2578;
    _S2593[int(14)] = _S2579;
    _S2593[int(13)] = _S2580;
    _S2593[int(12)] = _S2581;
    _S2593[int(11)] = _S2582;
    _S2593[int(10)] = _S2583;
    _S2593[int(19)] = _S2562.color_params_0.r_0.y;
    _S2593[int(8)] = _S2585;
    _S2593[int(7)] = _S2586;
    _S2593[int(6)] = _S2587;
    _S2593[int(5)] = _S2588;
    _S2593[int(4)] = _S2589;
    _S2593[int(3)] = _S2590;
    _S2593[int(2)] = _S2591;
    _S2593[int(1)] = _S2592;
    _S2593[int(0)] = _S2428.exposure_0;
    _S2593[int(28)] = _S2573;
    _S2593[int(37)] = _S2564;
    _S2593[int(36)] = _S2565;
    _S2593[int(35)] = _S2566;
    _S2593[int(34)] = _S2567;
    _S2593[int(33)] = _S2568;
    _S2593[int(32)] = _S2569;
    _S2593[int(31)] = _S2570;
    _S2593[int(30)] = _S2571;
    _S2593[int(29)] = _S2572;
    _S2593[int(38)] = _S2563;
    _S2593[int(27)] = _S2574;
    _S2593[int(26)] = _S2575;
    _S2593[int(25)] = _S2576;
    _S2593[int(24)] = _S2577;
    _S2593[int(23)] = _S2562.color_params_0.n_0.y;
    _S2593[int(22)] = _S2562.color_params_0.n_0.x;
    _S2593[int(21)] = _S2562.color_params_0.g_0.y;
    _S2593[int(20)] = _S2562.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S2593;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2594, FixedArray<float, 23>  * _S2595)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S2594, _S2595);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  * params_7, FixedArray<float, 23>  * grad_out_3, FixedArray<float, 39>  * _S2596)
{
    FixedArray<float, 39>  _S2597 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = *params_7;
    (&dp_params_3)->differential_0 = _S2597;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, grad_out_3);
    *_S2596 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  * raw_losses_2, int num_cameras_0, FixedArray<float, 6>  * loss_weights_0, FixedArray<float, 6>  * _S2598)
{
    float _S2599;
    FixedArray<float, 6>  losses_5;
    float _S2600 = float(num_cameras_0);
    float _S2601 = (*raw_losses_2)[int(0)] / _S2600;
    for(;;)
    {
        float _S2602 = (F32_abs((_S2601)));
        if(_S2602 < 0.10000000149011612f)
        {
            _S2599 = 0.5f * _S2601 * _S2601 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2599 = _S2602 - 0.05000000074505806f;
            break;
        }
    }
    losses_5[int(0)] = _S2599;
    losses_5[int(1)] = (*raw_losses_2)[int(1)] / (3.0f * _S2600);
    losses_5[int(2)] = ((*raw_losses_2)[int(2)] + (*raw_losses_2)[int(3)] + (*raw_losses_2)[int(4)]) / (9.0f * _S2600);
    losses_5[int(3)] = ((*raw_losses_2)[int(5)] + (*raw_losses_2)[int(6)] + (*raw_losses_2)[int(7)] + (*raw_losses_2)[int(8)] + (*raw_losses_2)[int(9)]) / (5.0f * _S2600);
    float _S2603 = (*raw_losses_2)[int(10)] / _S2600;
    for(;;)
    {
        float _S2604 = (F32_abs((_S2603)));
        if(_S2604 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2603 * _S2603 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2604 - 0.00249999994412065f;
            break;
        }
    }
    float _S2605;
    float _S2606 = (*raw_losses_2)[int(11)] / _S2600;
    for(;;)
    {
        float _S2607 = (F32_abs((_S2606)));
        if(_S2607 < 0.00499999988824129f)
        {
            _S2605 = 0.5f * _S2606 * _S2606 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2605 = _S2607 - 0.00249999994412065f;
            break;
        }
    }
    float _S2608 = _S2599 + _S2605;
    float _S2609 = (*raw_losses_2)[int(12)] / _S2600;
    for(;;)
    {
        float _S2610 = (F32_abs((_S2609)));
        if(_S2610 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2609 * _S2609 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2610 - 0.00249999994412065f;
            break;
        }
    }
    float _S2611 = _S2608 + _S2599;
    float _S2612 = (*raw_losses_2)[int(13)] / _S2600;
    for(;;)
    {
        float _S2613 = (F32_abs((_S2612)));
        if(_S2613 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2612 * _S2612 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2613 - 0.00249999994412065f;
            break;
        }
    }
    float _S2614 = _S2611 + _S2599;
    float _S2615 = (*raw_losses_2)[int(14)] / _S2600;
    for(;;)
    {
        float _S2616 = (F32_abs((_S2615)));
        if(_S2616 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2615 * _S2615 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2616 - 0.00249999994412065f;
            break;
        }
    }
    float _S2617 = _S2614 + _S2599;
    float _S2618 = (*raw_losses_2)[int(15)] / _S2600;
    for(;;)
    {
        float _S2619 = (F32_abs((_S2618)));
        if(_S2619 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2618 * _S2618 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2619 - 0.00249999994412065f;
            break;
        }
    }
    float _S2620 = _S2617 + _S2599;
    float _S2621 = (*raw_losses_2)[int(16)] / _S2600;
    for(;;)
    {
        float _S2622 = (F32_abs((_S2621)));
        if(_S2622 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2621 * _S2621 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2622 - 0.00249999994412065f;
            break;
        }
    }
    float _S2623 = _S2620 + _S2599;
    float _S2624 = (*raw_losses_2)[int(17)] / _S2600;
    for(;;)
    {
        float _S2625 = (F32_abs((_S2624)));
        if(_S2625 < 0.00499999988824129f)
        {
            _S2599 = 0.5f * _S2624 * _S2624 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2599 = _S2625 - 0.00249999994412065f;
            break;
        }
    }
    float _S2626 = (_S2623 + _S2599) / 8.0f;
    float _S2627 = ((*raw_losses_2)[int(18)] + (*raw_losses_2)[int(19)] + (*raw_losses_2)[int(20)] + (*raw_losses_2)[int(21)]) / (4.0f * _S2600);
    losses_5[int(0)] = losses_5[int(0)] * (*loss_weights_0)[int(0)];
    losses_5[int(1)] = losses_5[int(1)] * (*loss_weights_0)[int(1)];
    losses_5[int(2)] = losses_5[int(2)] * (*loss_weights_0)[int(2)];
    losses_5[int(3)] = losses_5[int(3)] * (*loss_weights_0)[int(3)];
    losses_5[int(4)] = _S2626 * (*loss_weights_0)[int(4)];
    losses_5[int(5)] = _S2627 * (*loss_weights_0)[int(5)];
    *_S2598 = losses_5;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  * raw_losses_3, int num_cameras_1, FixedArray<float, 6>  * loss_weights_1, FixedArray<float, 6>  * _S2628)
{
    float _S2629;
    FixedArray<float, 6>  losses_6;
    float _S2630 = float(num_cameras_1);
    float _S2631 = (*raw_losses_3)[int(0)] / _S2630;
    for(;;)
    {
        float _S2632 = (F32_abs((_S2631)));
        if(_S2632 < 0.10000000149011612f)
        {
            _S2629 = 0.5f * _S2631 * _S2631 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2629 = _S2632 - 0.05000000074505806f;
            break;
        }
    }
    losses_6[int(0)] = _S2629;
    losses_6[int(1)] = (*raw_losses_3)[int(1)] / (3.0f * _S2630);
    losses_6[int(2)] = ((*raw_losses_3)[int(2)] + (*raw_losses_3)[int(3)] + (*raw_losses_3)[int(4)]) / (9.0f * _S2630);
    float _S2633 = 5.0f * _S2630;
    losses_6[int(3)] = ((*raw_losses_3)[int(5)] + (*raw_losses_3)[int(6)] + (*raw_losses_3)[int(7)] + (*raw_losses_3)[int(8)] + (*raw_losses_3)[int(9)]) / _S2633;
    float _S2634 = (*raw_losses_3)[int(10)] / _S2630;
    for(;;)
    {
        float _S2635 = (F32_abs((_S2634)));
        if(_S2635 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2634 * _S2634 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2635 - 0.00249999994412065f;
            break;
        }
    }
    float _S2636;
    float _S2637 = (*raw_losses_3)[int(11)] / _S2630;
    for(;;)
    {
        float _S2638 = (F32_abs((_S2637)));
        if(_S2638 < 0.00499999988824129f)
        {
            _S2636 = 0.5f * _S2637 * _S2637 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2636 = _S2638 - 0.00249999994412065f;
            break;
        }
    }
    float _S2639 = _S2629 + _S2636;
    float _S2640 = (*raw_losses_3)[int(12)] / _S2630;
    for(;;)
    {
        float _S2641 = (F32_abs((_S2640)));
        if(_S2641 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2640 * _S2640 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2641 - 0.00249999994412065f;
            break;
        }
    }
    float _S2642 = _S2639 + _S2629;
    float _S2643 = (*raw_losses_3)[int(13)] / _S2630;
    for(;;)
    {
        float _S2644 = (F32_abs((_S2643)));
        if(_S2644 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2643 * _S2643 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2644 - 0.00249999994412065f;
            break;
        }
    }
    float _S2645 = _S2642 + _S2629;
    float _S2646 = (*raw_losses_3)[int(14)] / _S2630;
    for(;;)
    {
        float _S2647 = (F32_abs((_S2646)));
        if(_S2647 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2646 * _S2646 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2647 - 0.00249999994412065f;
            break;
        }
    }
    float _S2648 = _S2645 + _S2629;
    float _S2649 = (*raw_losses_3)[int(15)] / _S2630;
    for(;;)
    {
        float _S2650 = (F32_abs((_S2649)));
        if(_S2650 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2649 * _S2649 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2650 - 0.00249999994412065f;
            break;
        }
    }
    float _S2651 = _S2648 + _S2629;
    float _S2652 = (*raw_losses_3)[int(16)] / _S2630;
    for(;;)
    {
        float _S2653 = (F32_abs((_S2652)));
        if(_S2653 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2652 * _S2652 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2653 - 0.00249999994412065f;
            break;
        }
    }
    float _S2654 = _S2651 + _S2629;
    float _S2655 = (*raw_losses_3)[int(17)] / _S2630;
    for(;;)
    {
        float _S2656 = (F32_abs((_S2655)));
        if(_S2656 < 0.00499999988824129f)
        {
            _S2629 = 0.5f * _S2655 * _S2655 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2629 = _S2656 - 0.00249999994412065f;
            break;
        }
    }
    float _S2657 = (_S2654 + _S2629) / 8.0f;
    float _S2658 = ((*raw_losses_3)[int(18)] + (*raw_losses_3)[int(19)] + (*raw_losses_3)[int(20)] + (*raw_losses_3)[int(21)] + (*raw_losses_3)[int(22)]) / _S2633;
    losses_6[int(0)] = losses_6[int(0)] * (*loss_weights_1)[int(0)];
    losses_6[int(1)] = losses_6[int(1)] * (*loss_weights_1)[int(1)];
    losses_6[int(2)] = losses_6[int(2)] * (*loss_weights_1)[int(2)];
    losses_6[int(3)] = losses_6[int(3)] * (*loss_weights_1)[int(3)];
    losses_6[int(4)] = _S2657 * (*loss_weights_1)[int(4)];
    losses_6[int(5)] = _S2658 * (*loss_weights_1)[int(5)];
    *_S2628 = losses_6;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_15)
{
    FixedArray<float, 22>  _S2659 = dpraw_losses_1->primal_0;
    float _S2660 = float(num_cameras_2);
    float _S2661 = dpraw_losses_1->primal_0[int(0)] / _S2660;
    bool _S2662 = (s_primal_ctx_abs_0(_S2661)) < 0.10000000149011612f;
    float _S2663;
    if(_S2662)
    {
        _S2663 = 0.5f * _S2661;
    }
    else
    {
        _S2663 = 0.0f;
    }
    float _S2664 = 3.0f * _S2660;
    float _S2665 = 9.0f * _S2660;
    float _S2666 = 5.0f * _S2660;
    float _S2667 = _S2659[int(10)] / _S2660;
    bool _S2668 = (s_primal_ctx_abs_0(_S2667)) < 0.00499999988824129f;
    float _S2669;
    if(_S2668)
    {
        _S2669 = 0.5f * _S2667;
    }
    else
    {
        _S2669 = 0.0f;
    }
    float _S2670 = _S2659[int(11)] / _S2660;
    bool _S2671 = (s_primal_ctx_abs_0(_S2670)) < 0.00499999988824129f;
    float _S2672;
    if(_S2671)
    {
        _S2672 = 0.5f * _S2670;
    }
    else
    {
        _S2672 = 0.0f;
    }
    float _S2673 = _S2659[int(12)] / _S2660;
    bool _S2674 = (s_primal_ctx_abs_0(_S2673)) < 0.00499999988824129f;
    float _S2675;
    if(_S2674)
    {
        _S2675 = 0.5f * _S2673;
    }
    else
    {
        _S2675 = 0.0f;
    }
    float _S2676 = _S2659[int(13)] / _S2660;
    bool _S2677 = (s_primal_ctx_abs_0(_S2676)) < 0.00499999988824129f;
    float _S2678;
    if(_S2677)
    {
        _S2678 = 0.5f * _S2676;
    }
    else
    {
        _S2678 = 0.0f;
    }
    float _S2679 = _S2659[int(14)] / _S2660;
    bool _S2680 = (s_primal_ctx_abs_0(_S2679)) < 0.00499999988824129f;
    float _S2681;
    if(_S2680)
    {
        _S2681 = 0.5f * _S2679;
    }
    else
    {
        _S2681 = 0.0f;
    }
    float _S2682 = _S2659[int(15)] / _S2660;
    bool _S2683 = (s_primal_ctx_abs_0(_S2682)) < 0.00499999988824129f;
    float _S2684;
    if(_S2683)
    {
        _S2684 = 0.5f * _S2682;
    }
    else
    {
        _S2684 = 0.0f;
    }
    float _S2685 = _S2659[int(16)] / _S2660;
    bool _S2686 = (s_primal_ctx_abs_0(_S2685)) < 0.00499999988824129f;
    float _S2687;
    if(_S2686)
    {
        _S2687 = 0.5f * _S2685;
    }
    else
    {
        _S2687 = 0.0f;
    }
    float _S2688 = _S2659[int(17)] / _S2660;
    bool _S2689 = (s_primal_ctx_abs_0(_S2688)) < 0.00499999988824129f;
    float _S2690;
    if(_S2689)
    {
        _S2690 = 0.5f * _S2688;
    }
    else
    {
        _S2690 = 0.0f;
    }
    float _S2691 = (*loss_weights_2)[int(3)] * (*_s_dOut_15)[int(3)];
    float _S2692 = (*loss_weights_2)[int(2)] * (*_s_dOut_15)[int(2)];
    float _S2693 = (*loss_weights_2)[int(1)] * (*_s_dOut_15)[int(1)];
    float _S2694 = (*loss_weights_2)[int(0)] * (*_s_dOut_15)[int(0)];
    float _S2695 = (*loss_weights_2)[int(5)] * (*_s_dOut_15)[int(5)] / (4.0f * _S2660);
    float _S2696 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_15)[int(4)]);
    FixedArray<float, 22>  _S2697;
    _S2697[int(0)] = 0.0f;
    _S2697[int(1)] = 0.0f;
    _S2697[int(2)] = 0.0f;
    _S2697[int(3)] = 0.0f;
    _S2697[int(4)] = 0.0f;
    _S2697[int(5)] = 0.0f;
    _S2697[int(6)] = 0.0f;
    _S2697[int(7)] = 0.0f;
    _S2697[int(8)] = 0.0f;
    _S2697[int(9)] = 0.0f;
    _S2697[int(10)] = 0.0f;
    _S2697[int(11)] = 0.0f;
    _S2697[int(12)] = 0.0f;
    _S2697[int(13)] = 0.0f;
    _S2697[int(14)] = 0.0f;
    _S2697[int(15)] = 0.0f;
    _S2697[int(16)] = 0.0f;
    _S2697[int(17)] = 0.0f;
    _S2697[int(18)] = 0.0f;
    _S2697[int(19)] = 0.0f;
    _S2697[int(20)] = 0.0f;
    _S2697[int(21)] = 0.0f;
    _S2697[int(21)] = _S2695;
    _S2697[int(20)] = _S2695;
    _S2697[int(19)] = _S2695;
    _S2697[int(18)] = _S2695;
    float _S2698 = _S2697[int(0)];
    float _S2699 = _S2697[int(1)];
    float _S2700 = _S2697[int(2)];
    float _S2701 = _S2697[int(3)];
    float _S2702 = _S2697[int(4)];
    float _S2703 = _S2697[int(5)];
    float _S2704 = _S2697[int(6)];
    float _S2705 = _S2697[int(7)];
    float _S2706 = _S2697[int(8)];
    float _S2707 = _S2697[int(9)];
    float _S2708 = _S2697[int(10)];
    float _S2709 = _S2697[int(11)];
    float _S2710 = _S2697[int(12)];
    float _S2711 = _S2697[int(13)];
    float _S2712 = _S2697[int(14)];
    float _S2713 = _S2697[int(15)];
    float _S2714 = _S2697[int(16)];
    float _S2715 = _S2697[int(17)];
    float _S2716 = _S2697[int(18)];
    float _S2717 = _S2697[int(19)];
    float _S2718 = _S2697[int(20)];
    float _S2719 = _S2697[int(21)];
    float _S2720;
    if(_S2689)
    {
        float _S2721 = 200.0f * _S2696;
        float _S2722 = _S2690 * _S2721 + 0.5f * (_S2688 * _S2721);
        _S2690 = 0.0f;
        _S2720 = _S2722;
    }
    else
    {
        _S2690 = _S2696;
        _S2720 = 0.0f;
    }
    DiffPair_float_0 _S2723;
    (&_S2723)->primal_0 = _S2688;
    (&_S2723)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2723, _S2690);
    float _S2724 = (_S2723.differential_0 + _S2720) / _S2660;
    FixedArray<float, 22>  _S2725;
    _S2725[int(0)] = 0.0f;
    _S2725[int(1)] = 0.0f;
    _S2725[int(2)] = 0.0f;
    _S2725[int(3)] = 0.0f;
    _S2725[int(4)] = 0.0f;
    _S2725[int(5)] = 0.0f;
    _S2725[int(6)] = 0.0f;
    _S2725[int(7)] = 0.0f;
    _S2725[int(8)] = 0.0f;
    _S2725[int(9)] = 0.0f;
    _S2725[int(10)] = 0.0f;
    _S2725[int(11)] = 0.0f;
    _S2725[int(12)] = 0.0f;
    _S2725[int(13)] = 0.0f;
    _S2725[int(14)] = 0.0f;
    _S2725[int(15)] = 0.0f;
    _S2725[int(16)] = 0.0f;
    _S2725[int(17)] = 0.0f;
    _S2725[int(18)] = 0.0f;
    _S2725[int(19)] = 0.0f;
    _S2725[int(20)] = 0.0f;
    _S2725[int(21)] = 0.0f;
    _S2725[int(17)] = _S2724;
    float _S2726 = _S2698 + _S2725[int(0)];
    float _S2727 = _S2699 + _S2725[int(1)];
    float _S2728 = _S2700 + _S2725[int(2)];
    float _S2729 = _S2701 + _S2725[int(3)];
    float _S2730 = _S2702 + _S2725[int(4)];
    float _S2731 = _S2703 + _S2725[int(5)];
    float _S2732 = _S2704 + _S2725[int(6)];
    float _S2733 = _S2705 + _S2725[int(7)];
    float _S2734 = _S2706 + _S2725[int(8)];
    float _S2735 = _S2707 + _S2725[int(9)];
    float _S2736 = _S2708 + _S2725[int(10)];
    float _S2737 = _S2709 + _S2725[int(11)];
    float _S2738 = _S2710 + _S2725[int(12)];
    float _S2739 = _S2711 + _S2725[int(13)];
    float _S2740 = _S2712 + _S2725[int(14)];
    float _S2741 = _S2713 + _S2725[int(15)];
    float _S2742 = _S2714 + _S2725[int(16)];
    float _S2743 = _S2715 + _S2725[int(17)];
    float _S2744 = _S2716 + _S2725[int(18)];
    float _S2745 = _S2717 + _S2725[int(19)];
    float _S2746 = _S2718 + _S2725[int(20)];
    float _S2747 = _S2719 + _S2725[int(21)];
    if(_S2686)
    {
        float _S2748 = 200.0f * _S2696;
        float _S2749 = _S2687 * _S2748 + 0.5f * (_S2685 * _S2748);
        _S2687 = 0.0f;
        _S2690 = _S2749;
    }
    else
    {
        _S2687 = _S2696;
        _S2690 = 0.0f;
    }
    DiffPair_float_0 _S2750;
    (&_S2750)->primal_0 = _S2685;
    (&_S2750)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2750, _S2687);
    float _S2751 = (_S2750.differential_0 + _S2690) / _S2660;
    FixedArray<float, 22>  _S2752;
    _S2752[int(0)] = 0.0f;
    _S2752[int(1)] = 0.0f;
    _S2752[int(2)] = 0.0f;
    _S2752[int(3)] = 0.0f;
    _S2752[int(4)] = 0.0f;
    _S2752[int(5)] = 0.0f;
    _S2752[int(6)] = 0.0f;
    _S2752[int(7)] = 0.0f;
    _S2752[int(8)] = 0.0f;
    _S2752[int(9)] = 0.0f;
    _S2752[int(10)] = 0.0f;
    _S2752[int(11)] = 0.0f;
    _S2752[int(12)] = 0.0f;
    _S2752[int(13)] = 0.0f;
    _S2752[int(14)] = 0.0f;
    _S2752[int(15)] = 0.0f;
    _S2752[int(16)] = 0.0f;
    _S2752[int(17)] = 0.0f;
    _S2752[int(18)] = 0.0f;
    _S2752[int(19)] = 0.0f;
    _S2752[int(20)] = 0.0f;
    _S2752[int(21)] = 0.0f;
    _S2752[int(16)] = _S2751;
    float _S2753 = _S2726 + _S2752[int(0)];
    float _S2754 = _S2727 + _S2752[int(1)];
    float _S2755 = _S2728 + _S2752[int(2)];
    float _S2756 = _S2729 + _S2752[int(3)];
    float _S2757 = _S2730 + _S2752[int(4)];
    float _S2758 = _S2731 + _S2752[int(5)];
    float _S2759 = _S2732 + _S2752[int(6)];
    float _S2760 = _S2733 + _S2752[int(7)];
    float _S2761 = _S2734 + _S2752[int(8)];
    float _S2762 = _S2735 + _S2752[int(9)];
    float _S2763 = _S2736 + _S2752[int(10)];
    float _S2764 = _S2737 + _S2752[int(11)];
    float _S2765 = _S2738 + _S2752[int(12)];
    float _S2766 = _S2739 + _S2752[int(13)];
    float _S2767 = _S2740 + _S2752[int(14)];
    float _S2768 = _S2741 + _S2752[int(15)];
    float _S2769 = _S2742 + _S2752[int(16)];
    float _S2770 = _S2743 + _S2752[int(17)];
    float _S2771 = _S2744 + _S2752[int(18)];
    float _S2772 = _S2745 + _S2752[int(19)];
    float _S2773 = _S2746 + _S2752[int(20)];
    float _S2774 = _S2747 + _S2752[int(21)];
    if(_S2683)
    {
        float _S2775 = 200.0f * _S2696;
        float _S2776 = _S2684 * _S2775 + 0.5f * (_S2682 * _S2775);
        _S2684 = 0.0f;
        _S2687 = _S2776;
    }
    else
    {
        _S2684 = _S2696;
        _S2687 = 0.0f;
    }
    DiffPair_float_0 _S2777;
    (&_S2777)->primal_0 = _S2682;
    (&_S2777)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2777, _S2684);
    float _S2778 = (_S2777.differential_0 + _S2687) / _S2660;
    FixedArray<float, 22>  _S2779;
    _S2779[int(0)] = 0.0f;
    _S2779[int(1)] = 0.0f;
    _S2779[int(2)] = 0.0f;
    _S2779[int(3)] = 0.0f;
    _S2779[int(4)] = 0.0f;
    _S2779[int(5)] = 0.0f;
    _S2779[int(6)] = 0.0f;
    _S2779[int(7)] = 0.0f;
    _S2779[int(8)] = 0.0f;
    _S2779[int(9)] = 0.0f;
    _S2779[int(10)] = 0.0f;
    _S2779[int(11)] = 0.0f;
    _S2779[int(12)] = 0.0f;
    _S2779[int(13)] = 0.0f;
    _S2779[int(14)] = 0.0f;
    _S2779[int(15)] = 0.0f;
    _S2779[int(16)] = 0.0f;
    _S2779[int(17)] = 0.0f;
    _S2779[int(18)] = 0.0f;
    _S2779[int(19)] = 0.0f;
    _S2779[int(20)] = 0.0f;
    _S2779[int(21)] = 0.0f;
    _S2779[int(15)] = _S2778;
    float _S2780 = _S2753 + _S2779[int(0)];
    float _S2781 = _S2754 + _S2779[int(1)];
    float _S2782 = _S2755 + _S2779[int(2)];
    float _S2783 = _S2756 + _S2779[int(3)];
    float _S2784 = _S2757 + _S2779[int(4)];
    float _S2785 = _S2758 + _S2779[int(5)];
    float _S2786 = _S2759 + _S2779[int(6)];
    float _S2787 = _S2760 + _S2779[int(7)];
    float _S2788 = _S2761 + _S2779[int(8)];
    float _S2789 = _S2762 + _S2779[int(9)];
    float _S2790 = _S2763 + _S2779[int(10)];
    float _S2791 = _S2764 + _S2779[int(11)];
    float _S2792 = _S2765 + _S2779[int(12)];
    float _S2793 = _S2766 + _S2779[int(13)];
    float _S2794 = _S2767 + _S2779[int(14)];
    float _S2795 = _S2768 + _S2779[int(15)];
    float _S2796 = _S2769 + _S2779[int(16)];
    float _S2797 = _S2770 + _S2779[int(17)];
    float _S2798 = _S2771 + _S2779[int(18)];
    float _S2799 = _S2772 + _S2779[int(19)];
    float _S2800 = _S2773 + _S2779[int(20)];
    float _S2801 = _S2774 + _S2779[int(21)];
    if(_S2680)
    {
        float _S2802 = 200.0f * _S2696;
        float _S2803 = _S2681 * _S2802 + 0.5f * (_S2679 * _S2802);
        _S2681 = 0.0f;
        _S2684 = _S2803;
    }
    else
    {
        _S2681 = _S2696;
        _S2684 = 0.0f;
    }
    DiffPair_float_0 _S2804;
    (&_S2804)->primal_0 = _S2679;
    (&_S2804)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2804, _S2681);
    float _S2805 = (_S2804.differential_0 + _S2684) / _S2660;
    FixedArray<float, 22>  _S2806;
    _S2806[int(0)] = 0.0f;
    _S2806[int(1)] = 0.0f;
    _S2806[int(2)] = 0.0f;
    _S2806[int(3)] = 0.0f;
    _S2806[int(4)] = 0.0f;
    _S2806[int(5)] = 0.0f;
    _S2806[int(6)] = 0.0f;
    _S2806[int(7)] = 0.0f;
    _S2806[int(8)] = 0.0f;
    _S2806[int(9)] = 0.0f;
    _S2806[int(10)] = 0.0f;
    _S2806[int(11)] = 0.0f;
    _S2806[int(12)] = 0.0f;
    _S2806[int(13)] = 0.0f;
    _S2806[int(14)] = 0.0f;
    _S2806[int(15)] = 0.0f;
    _S2806[int(16)] = 0.0f;
    _S2806[int(17)] = 0.0f;
    _S2806[int(18)] = 0.0f;
    _S2806[int(19)] = 0.0f;
    _S2806[int(20)] = 0.0f;
    _S2806[int(21)] = 0.0f;
    _S2806[int(14)] = _S2805;
    float _S2807 = _S2780 + _S2806[int(0)];
    float _S2808 = _S2781 + _S2806[int(1)];
    float _S2809 = _S2782 + _S2806[int(2)];
    float _S2810 = _S2783 + _S2806[int(3)];
    float _S2811 = _S2784 + _S2806[int(4)];
    float _S2812 = _S2785 + _S2806[int(5)];
    float _S2813 = _S2786 + _S2806[int(6)];
    float _S2814 = _S2787 + _S2806[int(7)];
    float _S2815 = _S2788 + _S2806[int(8)];
    float _S2816 = _S2789 + _S2806[int(9)];
    float _S2817 = _S2790 + _S2806[int(10)];
    float _S2818 = _S2791 + _S2806[int(11)];
    float _S2819 = _S2792 + _S2806[int(12)];
    float _S2820 = _S2793 + _S2806[int(13)];
    float _S2821 = _S2794 + _S2806[int(14)];
    float _S2822 = _S2795 + _S2806[int(15)];
    float _S2823 = _S2796 + _S2806[int(16)];
    float _S2824 = _S2797 + _S2806[int(17)];
    float _S2825 = _S2798 + _S2806[int(18)];
    float _S2826 = _S2799 + _S2806[int(19)];
    float _S2827 = _S2800 + _S2806[int(20)];
    float _S2828 = _S2801 + _S2806[int(21)];
    if(_S2677)
    {
        float _S2829 = 200.0f * _S2696;
        float _S2830 = _S2678 * _S2829 + 0.5f * (_S2676 * _S2829);
        _S2678 = 0.0f;
        _S2681 = _S2830;
    }
    else
    {
        _S2678 = _S2696;
        _S2681 = 0.0f;
    }
    DiffPair_float_0 _S2831;
    (&_S2831)->primal_0 = _S2676;
    (&_S2831)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2831, _S2678);
    float _S2832 = (_S2831.differential_0 + _S2681) / _S2660;
    FixedArray<float, 22>  _S2833;
    _S2833[int(0)] = 0.0f;
    _S2833[int(1)] = 0.0f;
    _S2833[int(2)] = 0.0f;
    _S2833[int(3)] = 0.0f;
    _S2833[int(4)] = 0.0f;
    _S2833[int(5)] = 0.0f;
    _S2833[int(6)] = 0.0f;
    _S2833[int(7)] = 0.0f;
    _S2833[int(8)] = 0.0f;
    _S2833[int(9)] = 0.0f;
    _S2833[int(10)] = 0.0f;
    _S2833[int(11)] = 0.0f;
    _S2833[int(12)] = 0.0f;
    _S2833[int(13)] = 0.0f;
    _S2833[int(14)] = 0.0f;
    _S2833[int(15)] = 0.0f;
    _S2833[int(16)] = 0.0f;
    _S2833[int(17)] = 0.0f;
    _S2833[int(18)] = 0.0f;
    _S2833[int(19)] = 0.0f;
    _S2833[int(20)] = 0.0f;
    _S2833[int(21)] = 0.0f;
    _S2833[int(13)] = _S2832;
    float _S2834 = _S2807 + _S2833[int(0)];
    float _S2835 = _S2808 + _S2833[int(1)];
    float _S2836 = _S2809 + _S2833[int(2)];
    float _S2837 = _S2810 + _S2833[int(3)];
    float _S2838 = _S2811 + _S2833[int(4)];
    float _S2839 = _S2812 + _S2833[int(5)];
    float _S2840 = _S2813 + _S2833[int(6)];
    float _S2841 = _S2814 + _S2833[int(7)];
    float _S2842 = _S2815 + _S2833[int(8)];
    float _S2843 = _S2816 + _S2833[int(9)];
    float _S2844 = _S2817 + _S2833[int(10)];
    float _S2845 = _S2818 + _S2833[int(11)];
    float _S2846 = _S2819 + _S2833[int(12)];
    float _S2847 = _S2820 + _S2833[int(13)];
    float _S2848 = _S2821 + _S2833[int(14)];
    float _S2849 = _S2822 + _S2833[int(15)];
    float _S2850 = _S2823 + _S2833[int(16)];
    float _S2851 = _S2824 + _S2833[int(17)];
    float _S2852 = _S2825 + _S2833[int(18)];
    float _S2853 = _S2826 + _S2833[int(19)];
    float _S2854 = _S2827 + _S2833[int(20)];
    float _S2855 = _S2828 + _S2833[int(21)];
    if(_S2674)
    {
        float _S2856 = 200.0f * _S2696;
        float _S2857 = _S2675 * _S2856 + 0.5f * (_S2673 * _S2856);
        _S2675 = 0.0f;
        _S2678 = _S2857;
    }
    else
    {
        _S2675 = _S2696;
        _S2678 = 0.0f;
    }
    DiffPair_float_0 _S2858;
    (&_S2858)->primal_0 = _S2673;
    (&_S2858)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2858, _S2675);
    float _S2859 = (_S2858.differential_0 + _S2678) / _S2660;
    FixedArray<float, 22>  _S2860;
    _S2860[int(0)] = 0.0f;
    _S2860[int(1)] = 0.0f;
    _S2860[int(2)] = 0.0f;
    _S2860[int(3)] = 0.0f;
    _S2860[int(4)] = 0.0f;
    _S2860[int(5)] = 0.0f;
    _S2860[int(6)] = 0.0f;
    _S2860[int(7)] = 0.0f;
    _S2860[int(8)] = 0.0f;
    _S2860[int(9)] = 0.0f;
    _S2860[int(10)] = 0.0f;
    _S2860[int(11)] = 0.0f;
    _S2860[int(12)] = 0.0f;
    _S2860[int(13)] = 0.0f;
    _S2860[int(14)] = 0.0f;
    _S2860[int(15)] = 0.0f;
    _S2860[int(16)] = 0.0f;
    _S2860[int(17)] = 0.0f;
    _S2860[int(18)] = 0.0f;
    _S2860[int(19)] = 0.0f;
    _S2860[int(20)] = 0.0f;
    _S2860[int(21)] = 0.0f;
    _S2860[int(12)] = _S2859;
    float _S2861 = _S2834 + _S2860[int(0)];
    float _S2862 = _S2835 + _S2860[int(1)];
    float _S2863 = _S2836 + _S2860[int(2)];
    float _S2864 = _S2837 + _S2860[int(3)];
    float _S2865 = _S2838 + _S2860[int(4)];
    float _S2866 = _S2839 + _S2860[int(5)];
    float _S2867 = _S2840 + _S2860[int(6)];
    float _S2868 = _S2841 + _S2860[int(7)];
    float _S2869 = _S2842 + _S2860[int(8)];
    float _S2870 = _S2843 + _S2860[int(9)];
    float _S2871 = _S2844 + _S2860[int(10)];
    float _S2872 = _S2845 + _S2860[int(11)];
    float _S2873 = _S2846 + _S2860[int(12)];
    float _S2874 = _S2847 + _S2860[int(13)];
    float _S2875 = _S2848 + _S2860[int(14)];
    float _S2876 = _S2849 + _S2860[int(15)];
    float _S2877 = _S2850 + _S2860[int(16)];
    float _S2878 = _S2851 + _S2860[int(17)];
    float _S2879 = _S2852 + _S2860[int(18)];
    float _S2880 = _S2853 + _S2860[int(19)];
    float _S2881 = _S2854 + _S2860[int(20)];
    float _S2882 = _S2855 + _S2860[int(21)];
    if(_S2671)
    {
        float _S2883 = 200.0f * _S2696;
        float _S2884 = _S2672 * _S2883 + 0.5f * (_S2670 * _S2883);
        _S2672 = 0.0f;
        _S2675 = _S2884;
    }
    else
    {
        _S2672 = _S2696;
        _S2675 = 0.0f;
    }
    DiffPair_float_0 _S2885;
    (&_S2885)->primal_0 = _S2670;
    (&_S2885)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2885, _S2672);
    float _S2886 = (_S2885.differential_0 + _S2675) / _S2660;
    FixedArray<float, 22>  _S2887;
    _S2887[int(0)] = 0.0f;
    _S2887[int(1)] = 0.0f;
    _S2887[int(2)] = 0.0f;
    _S2887[int(3)] = 0.0f;
    _S2887[int(4)] = 0.0f;
    _S2887[int(5)] = 0.0f;
    _S2887[int(6)] = 0.0f;
    _S2887[int(7)] = 0.0f;
    _S2887[int(8)] = 0.0f;
    _S2887[int(9)] = 0.0f;
    _S2887[int(10)] = 0.0f;
    _S2887[int(11)] = 0.0f;
    _S2887[int(12)] = 0.0f;
    _S2887[int(13)] = 0.0f;
    _S2887[int(14)] = 0.0f;
    _S2887[int(15)] = 0.0f;
    _S2887[int(16)] = 0.0f;
    _S2887[int(17)] = 0.0f;
    _S2887[int(18)] = 0.0f;
    _S2887[int(19)] = 0.0f;
    _S2887[int(20)] = 0.0f;
    _S2887[int(21)] = 0.0f;
    _S2887[int(11)] = _S2886;
    float _S2888 = _S2861 + _S2887[int(0)];
    float _S2889 = _S2862 + _S2887[int(1)];
    float _S2890 = _S2863 + _S2887[int(2)];
    float _S2891 = _S2864 + _S2887[int(3)];
    float _S2892 = _S2865 + _S2887[int(4)];
    float _S2893 = _S2866 + _S2887[int(5)];
    float _S2894 = _S2867 + _S2887[int(6)];
    float _S2895 = _S2868 + _S2887[int(7)];
    float _S2896 = _S2869 + _S2887[int(8)];
    float _S2897 = _S2870 + _S2887[int(9)];
    float _S2898 = _S2871 + _S2887[int(10)];
    float _S2899 = _S2872 + _S2887[int(11)];
    float _S2900 = _S2873 + _S2887[int(12)];
    float _S2901 = _S2874 + _S2887[int(13)];
    float _S2902 = _S2875 + _S2887[int(14)];
    float _S2903 = _S2876 + _S2887[int(15)];
    float _S2904 = _S2877 + _S2887[int(16)];
    float _S2905 = _S2878 + _S2887[int(17)];
    float _S2906 = _S2879 + _S2887[int(18)];
    float _S2907 = _S2880 + _S2887[int(19)];
    float _S2908 = _S2881 + _S2887[int(20)];
    float _S2909 = _S2882 + _S2887[int(21)];
    if(_S2668)
    {
        float _S2910 = 200.0f * _S2696;
        float _S2911 = _S2669 * _S2910 + 0.5f * (_S2667 * _S2910);
        _S2669 = 0.0f;
        _S2672 = _S2911;
    }
    else
    {
        _S2669 = _S2696;
        _S2672 = 0.0f;
    }
    DiffPair_float_0 _S2912;
    (&_S2912)->primal_0 = _S2667;
    (&_S2912)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2912, _S2669);
    float _S2913 = (_S2912.differential_0 + _S2672) / _S2660;
    float _S2914 = _S2691 / _S2666;
    float _S2915 = _S2692 / _S2665;
    float _S2916 = _S2693 / _S2664;
    FixedArray<float, 22>  _S2917;
    _S2917[int(0)] = 0.0f;
    _S2917[int(1)] = 0.0f;
    _S2917[int(2)] = 0.0f;
    _S2917[int(3)] = 0.0f;
    _S2917[int(4)] = 0.0f;
    _S2917[int(5)] = 0.0f;
    _S2917[int(6)] = 0.0f;
    _S2917[int(7)] = 0.0f;
    _S2917[int(8)] = 0.0f;
    _S2917[int(9)] = 0.0f;
    _S2917[int(10)] = 0.0f;
    _S2917[int(11)] = 0.0f;
    _S2917[int(12)] = 0.0f;
    _S2917[int(13)] = 0.0f;
    _S2917[int(14)] = 0.0f;
    _S2917[int(15)] = 0.0f;
    _S2917[int(16)] = 0.0f;
    _S2917[int(17)] = 0.0f;
    _S2917[int(18)] = 0.0f;
    _S2917[int(19)] = 0.0f;
    _S2917[int(20)] = 0.0f;
    _S2917[int(21)] = 0.0f;
    _S2917[int(10)] = _S2913;
    _S2917[int(9)] = _S2914;
    _S2917[int(8)] = _S2914;
    _S2917[int(7)] = _S2914;
    _S2917[int(6)] = _S2914;
    _S2917[int(5)] = _S2914;
    _S2917[int(4)] = _S2915;
    _S2917[int(3)] = _S2915;
    _S2917[int(2)] = _S2915;
    _S2917[int(1)] = _S2916;
    float _S2918 = _S2888 + _S2917[int(0)];
    float _S2919 = _S2889 + _S2917[int(1)];
    float _S2920 = _S2890 + _S2917[int(2)];
    float _S2921 = _S2891 + _S2917[int(3)];
    float _S2922 = _S2892 + _S2917[int(4)];
    float _S2923 = _S2893 + _S2917[int(5)];
    float _S2924 = _S2894 + _S2917[int(6)];
    float _S2925 = _S2895 + _S2917[int(7)];
    float _S2926 = _S2896 + _S2917[int(8)];
    float _S2927 = _S2897 + _S2917[int(9)];
    float _S2928 = _S2898 + _S2917[int(10)];
    float _S2929 = _S2899 + _S2917[int(11)];
    float _S2930 = _S2900 + _S2917[int(12)];
    float _S2931 = _S2901 + _S2917[int(13)];
    float _S2932 = _S2902 + _S2917[int(14)];
    float _S2933 = _S2903 + _S2917[int(15)];
    float _S2934 = _S2904 + _S2917[int(16)];
    float _S2935 = _S2905 + _S2917[int(17)];
    float _S2936 = _S2906 + _S2917[int(18)];
    float _S2937 = _S2907 + _S2917[int(19)];
    float _S2938 = _S2908 + _S2917[int(20)];
    float _S2939 = _S2909 + _S2917[int(21)];
    if(_S2662)
    {
        float _S2940 = 10.0f * _S2694;
        float _S2941 = _S2663 * _S2940 + 0.5f * (_S2661 * _S2940);
        _S2663 = 0.0f;
        _S2669 = _S2941;
    }
    else
    {
        _S2663 = _S2694;
        _S2669 = 0.0f;
    }
    DiffPair_float_0 _S2942;
    (&_S2942)->primal_0 = _S2661;
    (&_S2942)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2942, _S2663);
    float _S2943 = (_S2942.differential_0 + _S2669) / _S2660;
    FixedArray<float, 22>  _S2944;
    _S2944[int(0)] = 0.0f;
    _S2944[int(1)] = 0.0f;
    _S2944[int(2)] = 0.0f;
    _S2944[int(3)] = 0.0f;
    _S2944[int(4)] = 0.0f;
    _S2944[int(5)] = 0.0f;
    _S2944[int(6)] = 0.0f;
    _S2944[int(7)] = 0.0f;
    _S2944[int(8)] = 0.0f;
    _S2944[int(9)] = 0.0f;
    _S2944[int(10)] = 0.0f;
    _S2944[int(11)] = 0.0f;
    _S2944[int(12)] = 0.0f;
    _S2944[int(13)] = 0.0f;
    _S2944[int(14)] = 0.0f;
    _S2944[int(15)] = 0.0f;
    _S2944[int(16)] = 0.0f;
    _S2944[int(17)] = 0.0f;
    _S2944[int(18)] = 0.0f;
    _S2944[int(19)] = 0.0f;
    _S2944[int(20)] = 0.0f;
    _S2944[int(21)] = 0.0f;
    _S2944[int(0)] = _S2943;
    FixedArray<float, 22>  _S2945 = {
        _S2918 + _S2944[int(0)], _S2919 + _S2944[int(1)], _S2920 + _S2944[int(2)], _S2921 + _S2944[int(3)], _S2922 + _S2944[int(4)], _S2923 + _S2944[int(5)], _S2924 + _S2944[int(6)], _S2925 + _S2944[int(7)], _S2926 + _S2944[int(8)], _S2927 + _S2944[int(9)], _S2928 + _S2944[int(10)], _S2929 + _S2944[int(11)], _S2930 + _S2944[int(12)], _S2931 + _S2944[int(13)], _S2932 + _S2944[int(14)], _S2933 + _S2944[int(15)], _S2934 + _S2944[int(16)], _S2935 + _S2944[int(17)], _S2936 + _S2944[int(18)], _S2937 + _S2944[int(19)], _S2938 + _S2944[int(20)], _S2939 + _S2944[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2945;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S2946, int _S2947, FixedArray<float, 6>  * _S2948, FixedArray<float, 6>  * _S2949)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2946, _S2947, _S2948, _S2949);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  * raw_losses_4, int num_cameras_3, FixedArray<float, 6>  * loss_weights_3, FixedArray<float, 6>  * grad_out_4, FixedArray<float, 22>  * _S2950)
{
    FixedArray<float, 22>  _S2951 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = *raw_losses_4;
    (&dp_raw_losses_1)->differential_0 = _S2951;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_3, loss_weights_3, grad_out_4);
    *_S2950 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_2, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_16)
{
    FixedArray<float, 23>  _S2952 = dpraw_losses_2->primal_0;
    float _S2953 = float(num_cameras_4);
    float _S2954 = dpraw_losses_2->primal_0[int(0)] / _S2953;
    bool _S2955 = (s_primal_ctx_abs_0(_S2954)) < 0.10000000149011612f;
    float _S2956;
    if(_S2955)
    {
        _S2956 = 0.5f * _S2954;
    }
    else
    {
        _S2956 = 0.0f;
    }
    float _S2957 = 3.0f * _S2953;
    float _S2958 = 9.0f * _S2953;
    float _S2959 = 5.0f * _S2953;
    float _S2960 = _S2952[int(10)] / _S2953;
    bool _S2961 = (s_primal_ctx_abs_0(_S2960)) < 0.00499999988824129f;
    float _S2962;
    if(_S2961)
    {
        _S2962 = 0.5f * _S2960;
    }
    else
    {
        _S2962 = 0.0f;
    }
    float _S2963 = _S2952[int(11)] / _S2953;
    bool _S2964 = (s_primal_ctx_abs_0(_S2963)) < 0.00499999988824129f;
    float _S2965;
    if(_S2964)
    {
        _S2965 = 0.5f * _S2963;
    }
    else
    {
        _S2965 = 0.0f;
    }
    float _S2966 = _S2952[int(12)] / _S2953;
    bool _S2967 = (s_primal_ctx_abs_0(_S2966)) < 0.00499999988824129f;
    float _S2968;
    if(_S2967)
    {
        _S2968 = 0.5f * _S2966;
    }
    else
    {
        _S2968 = 0.0f;
    }
    float _S2969 = _S2952[int(13)] / _S2953;
    bool _S2970 = (s_primal_ctx_abs_0(_S2969)) < 0.00499999988824129f;
    float _S2971;
    if(_S2970)
    {
        _S2971 = 0.5f * _S2969;
    }
    else
    {
        _S2971 = 0.0f;
    }
    float _S2972 = _S2952[int(14)] / _S2953;
    bool _S2973 = (s_primal_ctx_abs_0(_S2972)) < 0.00499999988824129f;
    float _S2974;
    if(_S2973)
    {
        _S2974 = 0.5f * _S2972;
    }
    else
    {
        _S2974 = 0.0f;
    }
    float _S2975 = _S2952[int(15)] / _S2953;
    bool _S2976 = (s_primal_ctx_abs_0(_S2975)) < 0.00499999988824129f;
    float _S2977;
    if(_S2976)
    {
        _S2977 = 0.5f * _S2975;
    }
    else
    {
        _S2977 = 0.0f;
    }
    float _S2978 = _S2952[int(16)] / _S2953;
    bool _S2979 = (s_primal_ctx_abs_0(_S2978)) < 0.00499999988824129f;
    float _S2980;
    if(_S2979)
    {
        _S2980 = 0.5f * _S2978;
    }
    else
    {
        _S2980 = 0.0f;
    }
    float _S2981 = _S2952[int(17)] / _S2953;
    bool _S2982 = (s_primal_ctx_abs_0(_S2981)) < 0.00499999988824129f;
    float _S2983;
    if(_S2982)
    {
        _S2983 = 0.5f * _S2981;
    }
    else
    {
        _S2983 = 0.0f;
    }
    float _S2984 = (*loss_weights_4)[int(3)] * (*_s_dOut_16)[int(3)];
    float _S2985 = (*loss_weights_4)[int(2)] * (*_s_dOut_16)[int(2)];
    float _S2986 = (*loss_weights_4)[int(1)] * (*_s_dOut_16)[int(1)];
    float _S2987 = (*loss_weights_4)[int(0)] * (*_s_dOut_16)[int(0)];
    float _S2988 = (*loss_weights_4)[int(5)] * (*_s_dOut_16)[int(5)] / _S2959;
    float _S2989 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_16)[int(4)]);
    FixedArray<float, 23>  _S2990;
    _S2990[int(0)] = 0.0f;
    _S2990[int(1)] = 0.0f;
    _S2990[int(2)] = 0.0f;
    _S2990[int(3)] = 0.0f;
    _S2990[int(4)] = 0.0f;
    _S2990[int(5)] = 0.0f;
    _S2990[int(6)] = 0.0f;
    _S2990[int(7)] = 0.0f;
    _S2990[int(8)] = 0.0f;
    _S2990[int(9)] = 0.0f;
    _S2990[int(10)] = 0.0f;
    _S2990[int(11)] = 0.0f;
    _S2990[int(12)] = 0.0f;
    _S2990[int(13)] = 0.0f;
    _S2990[int(14)] = 0.0f;
    _S2990[int(15)] = 0.0f;
    _S2990[int(16)] = 0.0f;
    _S2990[int(17)] = 0.0f;
    _S2990[int(18)] = 0.0f;
    _S2990[int(19)] = 0.0f;
    _S2990[int(20)] = 0.0f;
    _S2990[int(21)] = 0.0f;
    _S2990[int(22)] = 0.0f;
    _S2990[int(22)] = _S2988;
    _S2990[int(21)] = _S2988;
    _S2990[int(20)] = _S2988;
    _S2990[int(19)] = _S2988;
    _S2990[int(18)] = _S2988;
    float _S2991 = _S2990[int(0)];
    float _S2992 = _S2990[int(1)];
    float _S2993 = _S2990[int(2)];
    float _S2994 = _S2990[int(3)];
    float _S2995 = _S2990[int(4)];
    float _S2996 = _S2990[int(5)];
    float _S2997 = _S2990[int(6)];
    float _S2998 = _S2990[int(7)];
    float _S2999 = _S2990[int(8)];
    float _S3000 = _S2990[int(9)];
    float _S3001 = _S2990[int(10)];
    float _S3002 = _S2990[int(11)];
    float _S3003 = _S2990[int(12)];
    float _S3004 = _S2990[int(13)];
    float _S3005 = _S2990[int(14)];
    float _S3006 = _S2990[int(15)];
    float _S3007 = _S2990[int(16)];
    float _S3008 = _S2990[int(17)];
    float _S3009 = _S2990[int(18)];
    float _S3010 = _S2990[int(19)];
    float _S3011 = _S2990[int(20)];
    float _S3012 = _S2990[int(21)];
    float _S3013 = _S2990[int(22)];
    float _S3014;
    if(_S2982)
    {
        float _S3015 = 200.0f * _S2989;
        float _S3016 = _S2983 * _S3015 + 0.5f * (_S2981 * _S3015);
        _S2983 = 0.0f;
        _S3014 = _S3016;
    }
    else
    {
        _S2983 = _S2989;
        _S3014 = 0.0f;
    }
    DiffPair_float_0 _S3017;
    (&_S3017)->primal_0 = _S2981;
    (&_S3017)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3017, _S2983);
    float _S3018 = (_S3017.differential_0 + _S3014) / _S2953;
    FixedArray<float, 23>  _S3019;
    _S3019[int(0)] = 0.0f;
    _S3019[int(1)] = 0.0f;
    _S3019[int(2)] = 0.0f;
    _S3019[int(3)] = 0.0f;
    _S3019[int(4)] = 0.0f;
    _S3019[int(5)] = 0.0f;
    _S3019[int(6)] = 0.0f;
    _S3019[int(7)] = 0.0f;
    _S3019[int(8)] = 0.0f;
    _S3019[int(9)] = 0.0f;
    _S3019[int(10)] = 0.0f;
    _S3019[int(11)] = 0.0f;
    _S3019[int(12)] = 0.0f;
    _S3019[int(13)] = 0.0f;
    _S3019[int(14)] = 0.0f;
    _S3019[int(15)] = 0.0f;
    _S3019[int(16)] = 0.0f;
    _S3019[int(17)] = 0.0f;
    _S3019[int(18)] = 0.0f;
    _S3019[int(19)] = 0.0f;
    _S3019[int(20)] = 0.0f;
    _S3019[int(21)] = 0.0f;
    _S3019[int(22)] = 0.0f;
    _S3019[int(17)] = _S3018;
    float _S3020 = _S2991 + _S3019[int(0)];
    float _S3021 = _S2992 + _S3019[int(1)];
    float _S3022 = _S2993 + _S3019[int(2)];
    float _S3023 = _S2994 + _S3019[int(3)];
    float _S3024 = _S2995 + _S3019[int(4)];
    float _S3025 = _S2996 + _S3019[int(5)];
    float _S3026 = _S2997 + _S3019[int(6)];
    float _S3027 = _S2998 + _S3019[int(7)];
    float _S3028 = _S2999 + _S3019[int(8)];
    float _S3029 = _S3000 + _S3019[int(9)];
    float _S3030 = _S3001 + _S3019[int(10)];
    float _S3031 = _S3002 + _S3019[int(11)];
    float _S3032 = _S3003 + _S3019[int(12)];
    float _S3033 = _S3004 + _S3019[int(13)];
    float _S3034 = _S3005 + _S3019[int(14)];
    float _S3035 = _S3006 + _S3019[int(15)];
    float _S3036 = _S3007 + _S3019[int(16)];
    float _S3037 = _S3008 + _S3019[int(17)];
    float _S3038 = _S3009 + _S3019[int(18)];
    float _S3039 = _S3010 + _S3019[int(19)];
    float _S3040 = _S3011 + _S3019[int(20)];
    float _S3041 = _S3012 + _S3019[int(21)];
    float _S3042 = _S3013 + _S3019[int(22)];
    if(_S2979)
    {
        float _S3043 = 200.0f * _S2989;
        float _S3044 = _S2980 * _S3043 + 0.5f * (_S2978 * _S3043);
        _S2980 = 0.0f;
        _S2983 = _S3044;
    }
    else
    {
        _S2980 = _S2989;
        _S2983 = 0.0f;
    }
    DiffPair_float_0 _S3045;
    (&_S3045)->primal_0 = _S2978;
    (&_S3045)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3045, _S2980);
    float _S3046 = (_S3045.differential_0 + _S2983) / _S2953;
    FixedArray<float, 23>  _S3047;
    _S3047[int(0)] = 0.0f;
    _S3047[int(1)] = 0.0f;
    _S3047[int(2)] = 0.0f;
    _S3047[int(3)] = 0.0f;
    _S3047[int(4)] = 0.0f;
    _S3047[int(5)] = 0.0f;
    _S3047[int(6)] = 0.0f;
    _S3047[int(7)] = 0.0f;
    _S3047[int(8)] = 0.0f;
    _S3047[int(9)] = 0.0f;
    _S3047[int(10)] = 0.0f;
    _S3047[int(11)] = 0.0f;
    _S3047[int(12)] = 0.0f;
    _S3047[int(13)] = 0.0f;
    _S3047[int(14)] = 0.0f;
    _S3047[int(15)] = 0.0f;
    _S3047[int(16)] = 0.0f;
    _S3047[int(17)] = 0.0f;
    _S3047[int(18)] = 0.0f;
    _S3047[int(19)] = 0.0f;
    _S3047[int(20)] = 0.0f;
    _S3047[int(21)] = 0.0f;
    _S3047[int(22)] = 0.0f;
    _S3047[int(16)] = _S3046;
    float _S3048 = _S3020 + _S3047[int(0)];
    float _S3049 = _S3021 + _S3047[int(1)];
    float _S3050 = _S3022 + _S3047[int(2)];
    float _S3051 = _S3023 + _S3047[int(3)];
    float _S3052 = _S3024 + _S3047[int(4)];
    float _S3053 = _S3025 + _S3047[int(5)];
    float _S3054 = _S3026 + _S3047[int(6)];
    float _S3055 = _S3027 + _S3047[int(7)];
    float _S3056 = _S3028 + _S3047[int(8)];
    float _S3057 = _S3029 + _S3047[int(9)];
    float _S3058 = _S3030 + _S3047[int(10)];
    float _S3059 = _S3031 + _S3047[int(11)];
    float _S3060 = _S3032 + _S3047[int(12)];
    float _S3061 = _S3033 + _S3047[int(13)];
    float _S3062 = _S3034 + _S3047[int(14)];
    float _S3063 = _S3035 + _S3047[int(15)];
    float _S3064 = _S3036 + _S3047[int(16)];
    float _S3065 = _S3037 + _S3047[int(17)];
    float _S3066 = _S3038 + _S3047[int(18)];
    float _S3067 = _S3039 + _S3047[int(19)];
    float _S3068 = _S3040 + _S3047[int(20)];
    float _S3069 = _S3041 + _S3047[int(21)];
    float _S3070 = _S3042 + _S3047[int(22)];
    if(_S2976)
    {
        float _S3071 = 200.0f * _S2989;
        float _S3072 = _S2977 * _S3071 + 0.5f * (_S2975 * _S3071);
        _S2977 = 0.0f;
        _S2980 = _S3072;
    }
    else
    {
        _S2977 = _S2989;
        _S2980 = 0.0f;
    }
    DiffPair_float_0 _S3073;
    (&_S3073)->primal_0 = _S2975;
    (&_S3073)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3073, _S2977);
    float _S3074 = (_S3073.differential_0 + _S2980) / _S2953;
    FixedArray<float, 23>  _S3075;
    _S3075[int(0)] = 0.0f;
    _S3075[int(1)] = 0.0f;
    _S3075[int(2)] = 0.0f;
    _S3075[int(3)] = 0.0f;
    _S3075[int(4)] = 0.0f;
    _S3075[int(5)] = 0.0f;
    _S3075[int(6)] = 0.0f;
    _S3075[int(7)] = 0.0f;
    _S3075[int(8)] = 0.0f;
    _S3075[int(9)] = 0.0f;
    _S3075[int(10)] = 0.0f;
    _S3075[int(11)] = 0.0f;
    _S3075[int(12)] = 0.0f;
    _S3075[int(13)] = 0.0f;
    _S3075[int(14)] = 0.0f;
    _S3075[int(15)] = 0.0f;
    _S3075[int(16)] = 0.0f;
    _S3075[int(17)] = 0.0f;
    _S3075[int(18)] = 0.0f;
    _S3075[int(19)] = 0.0f;
    _S3075[int(20)] = 0.0f;
    _S3075[int(21)] = 0.0f;
    _S3075[int(22)] = 0.0f;
    _S3075[int(15)] = _S3074;
    float _S3076 = _S3048 + _S3075[int(0)];
    float _S3077 = _S3049 + _S3075[int(1)];
    float _S3078 = _S3050 + _S3075[int(2)];
    float _S3079 = _S3051 + _S3075[int(3)];
    float _S3080 = _S3052 + _S3075[int(4)];
    float _S3081 = _S3053 + _S3075[int(5)];
    float _S3082 = _S3054 + _S3075[int(6)];
    float _S3083 = _S3055 + _S3075[int(7)];
    float _S3084 = _S3056 + _S3075[int(8)];
    float _S3085 = _S3057 + _S3075[int(9)];
    float _S3086 = _S3058 + _S3075[int(10)];
    float _S3087 = _S3059 + _S3075[int(11)];
    float _S3088 = _S3060 + _S3075[int(12)];
    float _S3089 = _S3061 + _S3075[int(13)];
    float _S3090 = _S3062 + _S3075[int(14)];
    float _S3091 = _S3063 + _S3075[int(15)];
    float _S3092 = _S3064 + _S3075[int(16)];
    float _S3093 = _S3065 + _S3075[int(17)];
    float _S3094 = _S3066 + _S3075[int(18)];
    float _S3095 = _S3067 + _S3075[int(19)];
    float _S3096 = _S3068 + _S3075[int(20)];
    float _S3097 = _S3069 + _S3075[int(21)];
    float _S3098 = _S3070 + _S3075[int(22)];
    if(_S2973)
    {
        float _S3099 = 200.0f * _S2989;
        float _S3100 = _S2974 * _S3099 + 0.5f * (_S2972 * _S3099);
        _S2974 = 0.0f;
        _S2977 = _S3100;
    }
    else
    {
        _S2974 = _S2989;
        _S2977 = 0.0f;
    }
    DiffPair_float_0 _S3101;
    (&_S3101)->primal_0 = _S2972;
    (&_S3101)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3101, _S2974);
    float _S3102 = (_S3101.differential_0 + _S2977) / _S2953;
    FixedArray<float, 23>  _S3103;
    _S3103[int(0)] = 0.0f;
    _S3103[int(1)] = 0.0f;
    _S3103[int(2)] = 0.0f;
    _S3103[int(3)] = 0.0f;
    _S3103[int(4)] = 0.0f;
    _S3103[int(5)] = 0.0f;
    _S3103[int(6)] = 0.0f;
    _S3103[int(7)] = 0.0f;
    _S3103[int(8)] = 0.0f;
    _S3103[int(9)] = 0.0f;
    _S3103[int(10)] = 0.0f;
    _S3103[int(11)] = 0.0f;
    _S3103[int(12)] = 0.0f;
    _S3103[int(13)] = 0.0f;
    _S3103[int(14)] = 0.0f;
    _S3103[int(15)] = 0.0f;
    _S3103[int(16)] = 0.0f;
    _S3103[int(17)] = 0.0f;
    _S3103[int(18)] = 0.0f;
    _S3103[int(19)] = 0.0f;
    _S3103[int(20)] = 0.0f;
    _S3103[int(21)] = 0.0f;
    _S3103[int(22)] = 0.0f;
    _S3103[int(14)] = _S3102;
    float _S3104 = _S3076 + _S3103[int(0)];
    float _S3105 = _S3077 + _S3103[int(1)];
    float _S3106 = _S3078 + _S3103[int(2)];
    float _S3107 = _S3079 + _S3103[int(3)];
    float _S3108 = _S3080 + _S3103[int(4)];
    float _S3109 = _S3081 + _S3103[int(5)];
    float _S3110 = _S3082 + _S3103[int(6)];
    float _S3111 = _S3083 + _S3103[int(7)];
    float _S3112 = _S3084 + _S3103[int(8)];
    float _S3113 = _S3085 + _S3103[int(9)];
    float _S3114 = _S3086 + _S3103[int(10)];
    float _S3115 = _S3087 + _S3103[int(11)];
    float _S3116 = _S3088 + _S3103[int(12)];
    float _S3117 = _S3089 + _S3103[int(13)];
    float _S3118 = _S3090 + _S3103[int(14)];
    float _S3119 = _S3091 + _S3103[int(15)];
    float _S3120 = _S3092 + _S3103[int(16)];
    float _S3121 = _S3093 + _S3103[int(17)];
    float _S3122 = _S3094 + _S3103[int(18)];
    float _S3123 = _S3095 + _S3103[int(19)];
    float _S3124 = _S3096 + _S3103[int(20)];
    float _S3125 = _S3097 + _S3103[int(21)];
    float _S3126 = _S3098 + _S3103[int(22)];
    if(_S2970)
    {
        float _S3127 = 200.0f * _S2989;
        float _S3128 = _S2971 * _S3127 + 0.5f * (_S2969 * _S3127);
        _S2971 = 0.0f;
        _S2974 = _S3128;
    }
    else
    {
        _S2971 = _S2989;
        _S2974 = 0.0f;
    }
    DiffPair_float_0 _S3129;
    (&_S3129)->primal_0 = _S2969;
    (&_S3129)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3129, _S2971);
    float _S3130 = (_S3129.differential_0 + _S2974) / _S2953;
    FixedArray<float, 23>  _S3131;
    _S3131[int(0)] = 0.0f;
    _S3131[int(1)] = 0.0f;
    _S3131[int(2)] = 0.0f;
    _S3131[int(3)] = 0.0f;
    _S3131[int(4)] = 0.0f;
    _S3131[int(5)] = 0.0f;
    _S3131[int(6)] = 0.0f;
    _S3131[int(7)] = 0.0f;
    _S3131[int(8)] = 0.0f;
    _S3131[int(9)] = 0.0f;
    _S3131[int(10)] = 0.0f;
    _S3131[int(11)] = 0.0f;
    _S3131[int(12)] = 0.0f;
    _S3131[int(13)] = 0.0f;
    _S3131[int(14)] = 0.0f;
    _S3131[int(15)] = 0.0f;
    _S3131[int(16)] = 0.0f;
    _S3131[int(17)] = 0.0f;
    _S3131[int(18)] = 0.0f;
    _S3131[int(19)] = 0.0f;
    _S3131[int(20)] = 0.0f;
    _S3131[int(21)] = 0.0f;
    _S3131[int(22)] = 0.0f;
    _S3131[int(13)] = _S3130;
    float _S3132 = _S3104 + _S3131[int(0)];
    float _S3133 = _S3105 + _S3131[int(1)];
    float _S3134 = _S3106 + _S3131[int(2)];
    float _S3135 = _S3107 + _S3131[int(3)];
    float _S3136 = _S3108 + _S3131[int(4)];
    float _S3137 = _S3109 + _S3131[int(5)];
    float _S3138 = _S3110 + _S3131[int(6)];
    float _S3139 = _S3111 + _S3131[int(7)];
    float _S3140 = _S3112 + _S3131[int(8)];
    float _S3141 = _S3113 + _S3131[int(9)];
    float _S3142 = _S3114 + _S3131[int(10)];
    float _S3143 = _S3115 + _S3131[int(11)];
    float _S3144 = _S3116 + _S3131[int(12)];
    float _S3145 = _S3117 + _S3131[int(13)];
    float _S3146 = _S3118 + _S3131[int(14)];
    float _S3147 = _S3119 + _S3131[int(15)];
    float _S3148 = _S3120 + _S3131[int(16)];
    float _S3149 = _S3121 + _S3131[int(17)];
    float _S3150 = _S3122 + _S3131[int(18)];
    float _S3151 = _S3123 + _S3131[int(19)];
    float _S3152 = _S3124 + _S3131[int(20)];
    float _S3153 = _S3125 + _S3131[int(21)];
    float _S3154 = _S3126 + _S3131[int(22)];
    if(_S2967)
    {
        float _S3155 = 200.0f * _S2989;
        float _S3156 = _S2968 * _S3155 + 0.5f * (_S2966 * _S3155);
        _S2968 = 0.0f;
        _S2971 = _S3156;
    }
    else
    {
        _S2968 = _S2989;
        _S2971 = 0.0f;
    }
    DiffPair_float_0 _S3157;
    (&_S3157)->primal_0 = _S2966;
    (&_S3157)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3157, _S2968);
    float _S3158 = (_S3157.differential_0 + _S2971) / _S2953;
    FixedArray<float, 23>  _S3159;
    _S3159[int(0)] = 0.0f;
    _S3159[int(1)] = 0.0f;
    _S3159[int(2)] = 0.0f;
    _S3159[int(3)] = 0.0f;
    _S3159[int(4)] = 0.0f;
    _S3159[int(5)] = 0.0f;
    _S3159[int(6)] = 0.0f;
    _S3159[int(7)] = 0.0f;
    _S3159[int(8)] = 0.0f;
    _S3159[int(9)] = 0.0f;
    _S3159[int(10)] = 0.0f;
    _S3159[int(11)] = 0.0f;
    _S3159[int(12)] = 0.0f;
    _S3159[int(13)] = 0.0f;
    _S3159[int(14)] = 0.0f;
    _S3159[int(15)] = 0.0f;
    _S3159[int(16)] = 0.0f;
    _S3159[int(17)] = 0.0f;
    _S3159[int(18)] = 0.0f;
    _S3159[int(19)] = 0.0f;
    _S3159[int(20)] = 0.0f;
    _S3159[int(21)] = 0.0f;
    _S3159[int(22)] = 0.0f;
    _S3159[int(12)] = _S3158;
    float _S3160 = _S3132 + _S3159[int(0)];
    float _S3161 = _S3133 + _S3159[int(1)];
    float _S3162 = _S3134 + _S3159[int(2)];
    float _S3163 = _S3135 + _S3159[int(3)];
    float _S3164 = _S3136 + _S3159[int(4)];
    float _S3165 = _S3137 + _S3159[int(5)];
    float _S3166 = _S3138 + _S3159[int(6)];
    float _S3167 = _S3139 + _S3159[int(7)];
    float _S3168 = _S3140 + _S3159[int(8)];
    float _S3169 = _S3141 + _S3159[int(9)];
    float _S3170 = _S3142 + _S3159[int(10)];
    float _S3171 = _S3143 + _S3159[int(11)];
    float _S3172 = _S3144 + _S3159[int(12)];
    float _S3173 = _S3145 + _S3159[int(13)];
    float _S3174 = _S3146 + _S3159[int(14)];
    float _S3175 = _S3147 + _S3159[int(15)];
    float _S3176 = _S3148 + _S3159[int(16)];
    float _S3177 = _S3149 + _S3159[int(17)];
    float _S3178 = _S3150 + _S3159[int(18)];
    float _S3179 = _S3151 + _S3159[int(19)];
    float _S3180 = _S3152 + _S3159[int(20)];
    float _S3181 = _S3153 + _S3159[int(21)];
    float _S3182 = _S3154 + _S3159[int(22)];
    if(_S2964)
    {
        float _S3183 = 200.0f * _S2989;
        float _S3184 = _S2965 * _S3183 + 0.5f * (_S2963 * _S3183);
        _S2965 = 0.0f;
        _S2968 = _S3184;
    }
    else
    {
        _S2965 = _S2989;
        _S2968 = 0.0f;
    }
    DiffPair_float_0 _S3185;
    (&_S3185)->primal_0 = _S2963;
    (&_S3185)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3185, _S2965);
    float _S3186 = (_S3185.differential_0 + _S2968) / _S2953;
    FixedArray<float, 23>  _S3187;
    _S3187[int(0)] = 0.0f;
    _S3187[int(1)] = 0.0f;
    _S3187[int(2)] = 0.0f;
    _S3187[int(3)] = 0.0f;
    _S3187[int(4)] = 0.0f;
    _S3187[int(5)] = 0.0f;
    _S3187[int(6)] = 0.0f;
    _S3187[int(7)] = 0.0f;
    _S3187[int(8)] = 0.0f;
    _S3187[int(9)] = 0.0f;
    _S3187[int(10)] = 0.0f;
    _S3187[int(11)] = 0.0f;
    _S3187[int(12)] = 0.0f;
    _S3187[int(13)] = 0.0f;
    _S3187[int(14)] = 0.0f;
    _S3187[int(15)] = 0.0f;
    _S3187[int(16)] = 0.0f;
    _S3187[int(17)] = 0.0f;
    _S3187[int(18)] = 0.0f;
    _S3187[int(19)] = 0.0f;
    _S3187[int(20)] = 0.0f;
    _S3187[int(21)] = 0.0f;
    _S3187[int(22)] = 0.0f;
    _S3187[int(11)] = _S3186;
    float _S3188 = _S3160 + _S3187[int(0)];
    float _S3189 = _S3161 + _S3187[int(1)];
    float _S3190 = _S3162 + _S3187[int(2)];
    float _S3191 = _S3163 + _S3187[int(3)];
    float _S3192 = _S3164 + _S3187[int(4)];
    float _S3193 = _S3165 + _S3187[int(5)];
    float _S3194 = _S3166 + _S3187[int(6)];
    float _S3195 = _S3167 + _S3187[int(7)];
    float _S3196 = _S3168 + _S3187[int(8)];
    float _S3197 = _S3169 + _S3187[int(9)];
    float _S3198 = _S3170 + _S3187[int(10)];
    float _S3199 = _S3171 + _S3187[int(11)];
    float _S3200 = _S3172 + _S3187[int(12)];
    float _S3201 = _S3173 + _S3187[int(13)];
    float _S3202 = _S3174 + _S3187[int(14)];
    float _S3203 = _S3175 + _S3187[int(15)];
    float _S3204 = _S3176 + _S3187[int(16)];
    float _S3205 = _S3177 + _S3187[int(17)];
    float _S3206 = _S3178 + _S3187[int(18)];
    float _S3207 = _S3179 + _S3187[int(19)];
    float _S3208 = _S3180 + _S3187[int(20)];
    float _S3209 = _S3181 + _S3187[int(21)];
    float _S3210 = _S3182 + _S3187[int(22)];
    if(_S2961)
    {
        float _S3211 = 200.0f * _S2989;
        float _S3212 = _S2962 * _S3211 + 0.5f * (_S2960 * _S3211);
        _S2962 = 0.0f;
        _S2965 = _S3212;
    }
    else
    {
        _S2962 = _S2989;
        _S2965 = 0.0f;
    }
    DiffPair_float_0 _S3213;
    (&_S3213)->primal_0 = _S2960;
    (&_S3213)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3213, _S2962);
    float _S3214 = (_S3213.differential_0 + _S2965) / _S2953;
    float _S3215 = _S2984 / _S2959;
    float _S3216 = _S2985 / _S2958;
    float _S3217 = _S2986 / _S2957;
    FixedArray<float, 23>  _S3218;
    _S3218[int(0)] = 0.0f;
    _S3218[int(1)] = 0.0f;
    _S3218[int(2)] = 0.0f;
    _S3218[int(3)] = 0.0f;
    _S3218[int(4)] = 0.0f;
    _S3218[int(5)] = 0.0f;
    _S3218[int(6)] = 0.0f;
    _S3218[int(7)] = 0.0f;
    _S3218[int(8)] = 0.0f;
    _S3218[int(9)] = 0.0f;
    _S3218[int(10)] = 0.0f;
    _S3218[int(11)] = 0.0f;
    _S3218[int(12)] = 0.0f;
    _S3218[int(13)] = 0.0f;
    _S3218[int(14)] = 0.0f;
    _S3218[int(15)] = 0.0f;
    _S3218[int(16)] = 0.0f;
    _S3218[int(17)] = 0.0f;
    _S3218[int(18)] = 0.0f;
    _S3218[int(19)] = 0.0f;
    _S3218[int(20)] = 0.0f;
    _S3218[int(21)] = 0.0f;
    _S3218[int(22)] = 0.0f;
    _S3218[int(10)] = _S3214;
    _S3218[int(9)] = _S3215;
    _S3218[int(8)] = _S3215;
    _S3218[int(7)] = _S3215;
    _S3218[int(6)] = _S3215;
    _S3218[int(5)] = _S3215;
    _S3218[int(4)] = _S3216;
    _S3218[int(3)] = _S3216;
    _S3218[int(2)] = _S3216;
    _S3218[int(1)] = _S3217;
    float _S3219 = _S3188 + _S3218[int(0)];
    float _S3220 = _S3189 + _S3218[int(1)];
    float _S3221 = _S3190 + _S3218[int(2)];
    float _S3222 = _S3191 + _S3218[int(3)];
    float _S3223 = _S3192 + _S3218[int(4)];
    float _S3224 = _S3193 + _S3218[int(5)];
    float _S3225 = _S3194 + _S3218[int(6)];
    float _S3226 = _S3195 + _S3218[int(7)];
    float _S3227 = _S3196 + _S3218[int(8)];
    float _S3228 = _S3197 + _S3218[int(9)];
    float _S3229 = _S3198 + _S3218[int(10)];
    float _S3230 = _S3199 + _S3218[int(11)];
    float _S3231 = _S3200 + _S3218[int(12)];
    float _S3232 = _S3201 + _S3218[int(13)];
    float _S3233 = _S3202 + _S3218[int(14)];
    float _S3234 = _S3203 + _S3218[int(15)];
    float _S3235 = _S3204 + _S3218[int(16)];
    float _S3236 = _S3205 + _S3218[int(17)];
    float _S3237 = _S3206 + _S3218[int(18)];
    float _S3238 = _S3207 + _S3218[int(19)];
    float _S3239 = _S3208 + _S3218[int(20)];
    float _S3240 = _S3209 + _S3218[int(21)];
    float _S3241 = _S3210 + _S3218[int(22)];
    if(_S2955)
    {
        float _S3242 = 10.0f * _S2987;
        float _S3243 = _S2956 * _S3242 + 0.5f * (_S2954 * _S3242);
        _S2956 = 0.0f;
        _S2962 = _S3243;
    }
    else
    {
        _S2956 = _S2987;
        _S2962 = 0.0f;
    }
    DiffPair_float_0 _S3244;
    (&_S3244)->primal_0 = _S2954;
    (&_S3244)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3244, _S2956);
    float _S3245 = (_S3244.differential_0 + _S2962) / _S2953;
    FixedArray<float, 23>  _S3246;
    _S3246[int(0)] = 0.0f;
    _S3246[int(1)] = 0.0f;
    _S3246[int(2)] = 0.0f;
    _S3246[int(3)] = 0.0f;
    _S3246[int(4)] = 0.0f;
    _S3246[int(5)] = 0.0f;
    _S3246[int(6)] = 0.0f;
    _S3246[int(7)] = 0.0f;
    _S3246[int(8)] = 0.0f;
    _S3246[int(9)] = 0.0f;
    _S3246[int(10)] = 0.0f;
    _S3246[int(11)] = 0.0f;
    _S3246[int(12)] = 0.0f;
    _S3246[int(13)] = 0.0f;
    _S3246[int(14)] = 0.0f;
    _S3246[int(15)] = 0.0f;
    _S3246[int(16)] = 0.0f;
    _S3246[int(17)] = 0.0f;
    _S3246[int(18)] = 0.0f;
    _S3246[int(19)] = 0.0f;
    _S3246[int(20)] = 0.0f;
    _S3246[int(21)] = 0.0f;
    _S3246[int(22)] = 0.0f;
    _S3246[int(0)] = _S3245;
    FixedArray<float, 23>  _S3247 = {
        _S3219 + _S3246[int(0)], _S3220 + _S3246[int(1)], _S3221 + _S3246[int(2)], _S3222 + _S3246[int(3)], _S3223 + _S3246[int(4)], _S3224 + _S3246[int(5)], _S3225 + _S3246[int(6)], _S3226 + _S3246[int(7)], _S3227 + _S3246[int(8)], _S3228 + _S3246[int(9)], _S3229 + _S3246[int(10)], _S3230 + _S3246[int(11)], _S3231 + _S3246[int(12)], _S3232 + _S3246[int(13)], _S3233 + _S3246[int(14)], _S3234 + _S3246[int(15)], _S3235 + _S3246[int(16)], _S3236 + _S3246[int(17)], _S3237 + _S3246[int(18)], _S3238 + _S3246[int(19)], _S3239 + _S3246[int(20)], _S3240 + _S3246[int(21)], _S3241 + _S3246[int(22)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S3247;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S3248, int _S3249, FixedArray<float, 6>  * _S3250, FixedArray<float, 6>  * _S3251)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S3248, _S3249, _S3250, _S3251);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  * raw_losses_5, int num_cameras_5, FixedArray<float, 6>  * loss_weights_5, FixedArray<float, 6>  * grad_out_5, FixedArray<float, 23>  * _S3252)
{
    FixedArray<float, 23>  _S3253 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = *raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S3253;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_2, num_cameras_5, loss_weights_5, grad_out_5);
    *_S3252 = (&dp_raw_losses_2)->differential_0;
    return;
}

