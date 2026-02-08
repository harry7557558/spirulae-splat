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

struct CRFPPISPChannelParams_0
{
    float toe_0;
    float shoulder_0;
    float gamma_0;
    float center_0;
};

inline __device__ CRFPPISPChannelParams_0 CRFPPISPChannelParams_x24_syn_dzero_0()
{
    CRFPPISPChannelParams_0 result_2;
    (&result_2)->toe_0 = 0.0f;
    (&result_2)->shoulder_0 = 0.0f;
    (&result_2)->gamma_0 = 0.0f;
    (&result_2)->center_0 = 0.0f;
    return result_2;
}

struct PPISPParams_0
{
    float exposure_0;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_0;
    ColorPPISPParams_0 color_params_0;
    FixedArray<CRFPPISPChannelParams_0, 3>  crf_params_0;
};

inline __device__ PPISPParams_0 PPISPParams_x24_syn_dzero_0()
{
    PPISPParams_0 result_3;
    (&result_3)->exposure_0 = 0.0f;
    VignettingChannelParams_0 _S2 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_3)->vignette_params_0[int(0)] = _S2;
    (&result_3)->vignette_params_0[int(1)] = _S2;
    (&result_3)->vignette_params_0[int(2)] = _S2;
    (&result_3)->color_params_0 = ColorPPISPParams_x24_syn_dzero_0();
    CRFPPISPChannelParams_0 _S3 = CRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_3)->crf_params_0[int(0)] = _S3;
    (&result_3)->crf_params_0[int(1)] = _S3;
    (&result_3)->crf_params_0[int(2)] = _S3;
    return result_3;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S4 = (F32_exp(((*dpx_0).primal_0))) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S4;
    return;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S5 = *dpx_1;
    float _S6;
    if(((*dpx_1).primal_0) > ((*dpy_0).primal_0))
    {
        _S6 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_0).primal_0))
        {
            _S6 = 0.0f;
        }
        else
        {
            _S6 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S5.primal_0;
    dpx_1->differential_0 = _S6;
    DiffPair_float_0 _S7 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S5.primal_0))
    {
        _S6 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_1).primal_0))
        {
            _S6 = 0.0f;
        }
        else
        {
            _S6 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S7.primal_0;
    dpy_0->differential_0 = _S6;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S8 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_2).primal_0)))))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S8;
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
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_5 = result_4 + _slang_vector_get_element(x_0, i_0) * _slang_vector_get_element(y_0, i_0);
        i_0 = i_0 + int(1);
        result_4 = result_5;
    }
    return result_4;
}

inline __device__ float dot_1(float4  x_1, float4  y_1)
{
    int i_1 = int(0);
    float result_6 = 0.0f;
    for(;;)
    {
        if(i_1 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_7 = result_6 + _slang_vector_get_element(x_1, i_1) * _slang_vector_get_element(y_1, i_1);
        i_1 = i_1 + int(1);
        result_6 = result_7;
    }
    return result_6;
}

inline __device__ float dot_2(float2  x_2, float2  y_2)
{
    int i_2 = int(0);
    float result_8 = 0.0f;
    for(;;)
    {
        if(i_2 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_9 = result_8 + _slang_vector_get_element(x_2, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_8 = result_9;
    }
    return result_8;
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
    float _S9 = 1.0f / (*dpx_4).primal_0 * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S9;
    return;
}

inline __device__ float3  exp_0(float3  x_6)
{
    float3  result_10;
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
        *_slang_vector_get_element_ptr(&result_10, i_3) = (F32_exp((_slang_vector_get_element(x_6, i_3))));
        i_3 = i_3 + int(1);
    }
    return result_10;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float3  dOut_5)
{
    float3  _S10 = exp_0((*dpx_5).primal_0) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S10;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ float2  exp_1(float2  x_7)
{
    float2  result_11;
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
        *_slang_vector_get_element_ptr(&result_11, i_4) = (F32_exp((_slang_vector_get_element(x_7, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_11;
}

inline __device__ void _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_6, float2  dOut_6)
{
    float2  _S11 = exp_1((*dpx_6).primal_0) * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S11;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpy_2, float dOut_7)
{
    DiffPair_float_0 _S12 = *dpx_7;
    float _S13;
    if(((*dpx_7).primal_0) < ((*dpy_2).primal_0))
    {
        _S13 = dOut_7;
    }
    else
    {
        if(((*dpx_7).primal_0) > ((*dpy_2).primal_0))
        {
            _S13 = 0.0f;
        }
        else
        {
            _S13 = 0.5f * dOut_7;
        }
    }
    dpx_7->primal_0 = _S12.primal_0;
    dpx_7->differential_0 = _S13;
    DiffPair_float_0 _S14 = *dpy_2;
    if(((*dpy_2).primal_0) < (_S12.primal_0))
    {
        _S13 = dOut_7;
    }
    else
    {
        if(((*dpy_2).primal_0) > ((*dpx_7).primal_0))
        {
            _S13 = 0.0f;
        }
        else
        {
            _S13 = 0.5f * dOut_7;
        }
    }
    dpy_2->primal_0 = _S14.primal_0;
    dpy_2->differential_0 = _S13;
    return;
}

inline __device__ void per_splat_losses(bool is_3dgs_0, float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S15)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    if(is_3dgs_0)
    {
        float3  _S16 = exp_0(scales_0);
        float _S17 = _S16.x;
        float _S18 = _S16.y;
        float _S19 = _S16.z;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S17 + _S18 + _S19) / 3.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S17), (_S18)))), (_S19))) / (F32_min(((F32_min((_S17), (_S18)))), (_S19)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float3  _S20 = exp_0(make_float3 (2.0f) * scales_0);
        float x_8 = _S20.x;
        float y_3 = _S20.y;
        float z_0 = _S20.z;
        float s_0 = x_8 + y_3 + z_0;
        float s1_0 = (F32_max(((F32_max((x_8), (y_3)))), (z_0))) / s_0;
        float s3_0 = (F32_min(((F32_min((x_8), (y_3)))), (z_0))) / s_0;
        float s2_0 = 1.0f - s1_0 - s3_0;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    }
    else
    {
        float2  _S21 = float2 {scales_0.x, scales_0.y};
        float2  _S22 = exp_1(_S21);
        float _S23 = _S22.x;
        float _S24 = _S22.y;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S23 + _S24) / 2.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max((_S23), (_S24))) / (F32_min((_S23), (_S24)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float2  _S25 = exp_1(make_float2 (2.0f) * _S21);
        float x_9 = _S25.x;
        float y_4 = _S25.y;
        float s_1 = x_9 + y_4;
        float s1_1 = (F32_max((x_9), (y_4))) / s_1;
        float s2_1 = (F32_min((x_9), (y_4))) / s_1;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_1 * (F32_log((s1_1))) - s2_1 * (F32_log((s2_1)))))) - 0.99998998641967773f)))), (0.0f)));
    }
    *_S15 = losses_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ float s_primal_ctx_exp_0(float _S26)
{
    return (F32_exp((_S26)));
}

inline __device__ float2  s_primal_ctx_exp_1(float2  _S27)
{
    return exp_1(_S27);
}

inline __device__ float s_primal_ctx_max_0(float _S28, float _S29)
{
    return (F32_max((_S28), (_S29)));
}

inline __device__ float s_primal_ctx_min_0(float _S30, float _S31)
{
    return (F32_min((_S30), (_S31)));
}

inline __device__ float s_primal_ctx_log_0(float _S32)
{
    return (F32_log((_S32)));
}

inline __device__ float3  s_primal_ctx_exp_2(float3  _S33)
{
    return exp_0(_S33);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S34, DiffPair_float_0 * _S35, float _S36)
{
    _d_max_0(_S34, _S35, _S36);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S37, float _S38)
{
    _d_log_0(_S37, _S38);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S39, float _S40)
{
    _d_exp_0(_S39, _S40);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S41, DiffPair_float_0 * _S42, float _S43)
{
    _d_min_0(_S41, _S42, _S43);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S44, float2  _S45)
{
    _d_exp_vector_1(_S44, _S45);
    return;
}

inline __device__ void s_bwd_prop_exp_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S46, float3  _S47)
{
    _d_exp_vector_0(_S46, _S47);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S48, float _S49)
{
    _d_sqrt_0(_S48, _S49);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_8, float _s_dOut_0)
{
    float _S50 = (*dpx_8).primal_0.x;
    float _S51 = (*dpx_8).primal_0.y;
    float _S52 = (*dpx_8).primal_0.z;
    float _S53 = (*dpx_8).primal_0.w;
    DiffPair_float_0 _S54;
    (&_S54)->primal_0 = _S50 * _S50 + _S51 * _S51 + _S52 * _S52 + _S53 * _S53;
    (&_S54)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S54, _s_dOut_0);
    float _S55 = (*dpx_8).primal_0.w * _S54.differential_0;
    float _S56 = _S55 + _S55;
    float _S57 = (*dpx_8).primal_0.z * _S54.differential_0;
    float _S58 = _S57 + _S57;
    float _S59 = (*dpx_8).primal_0.y * _S54.differential_0;
    float _S60 = _S59 + _S59;
    float _S61 = (*dpx_8).primal_0.x * _S54.differential_0;
    float _S62 = _S61 + _S61;
    float4  _S63 = make_float4 (0.0f);
    *&((&_S63)->w) = _S56;
    *&((&_S63)->z) = _S58;
    *&((&_S63)->y) = _S60;
    *&((&_S63)->x) = _S62;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S63;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S64, float _S65)
{
    s_bwd_prop_length_impl_0(_S64, _S65);
    return;
}

inline __device__ void s_bwd_prop_per_splat_losses_0(bool is_3dgs_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1, FixedArray<float, 5>  * _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S66 = *dpscales_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S67 = *dpquat_0;
    float2  _S68 = make_float2 (0.0f);
    float3  _S69 = make_float3 (0.0f);
    float _S70 = - (*dpopacity_0).primal_0;
    float _S71 = 1.0f + s_primal_ctx_exp_0(_S70);
    float _S72 = _S71 * _S71;
    float _S73 = length_0((*dpquat_0).primal_0);
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
    float _S114;
    float _S115;
    float _S116;
    float _S117;
    float _S118;
    float _S119;
    float _S120;
    float _S121;
    float2  _S122;
    float2  _S123;
    float3  _S124;
    if(is_3dgs_1)
    {
        float3  _S125 = s_primal_ctx_exp_2(_S66.primal_0);
        float _S126 = _S125.x;
        float _S127 = _S125.y;
        float _S128 = _S125.z;
        float _S129 = s_primal_ctx_max_0(_S126, _S127);
        float _S130 = s_primal_ctx_max_0(_S129, _S128);
        float _S131 = s_primal_ctx_min_0(_S126, _S127);
        float _S132 = s_primal_ctx_min_0(_S131, _S128);
        float _S133 = _S130 / _S132;
        float _S134 = _S132 * _S132;
        float3  _S135 = make_float3 (2.0f) * _S66.primal_0;
        float3  _S136 = s_primal_ctx_exp_2(_S135);
        float x_10 = _S136.x;
        float y_5 = _S136.y;
        float z_1 = _S136.z;
        float s_2 = x_10 + y_5 + z_1;
        float _S137 = s_primal_ctx_max_0(x_10, y_5);
        float _S138 = s_primal_ctx_max_0(_S137, z_1);
        float s1_2 = _S138 / s_2;
        float _S139 = s_2 * s_2;
        float _S140 = s_primal_ctx_min_0(x_10, y_5);
        float _S141 = s_primal_ctx_min_0(_S140, z_1);
        float s3_1 = _S141 / s_2;
        float s2_2 = 1.0f - s1_2 - s3_1;
        float _S142 = - s1_2;
        float _S143 = s_primal_ctx_log_0(s1_2);
        float _S144 = s_primal_ctx_log_0(s2_2);
        float _S145 = s_primal_ctx_log_0(s3_1);
        float _S146 = _S142 * _S143 - s2_2 * _S144 - s3_1 * _S145;
        float _S147 = s_primal_ctx_exp_0(_S146) - 0.99998998641967773f;
        float _S148 = - s_primal_ctx_log_0(_S147);
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
        _S85 = 0.0f;
        _S86 = 0.0f;
        _S87 = 0.0f;
        _S122 = _S68;
        _S88 = 0.0f;
        _S89 = 0.0f;
        _S90 = 0.0f;
        _S91 = 0.0f;
        _S92 = 0.0f;
        _S93 = 0.0f;
        _S123 = _S68;
        _S94 = _S148;
        _S95 = _S147;
        _S96 = _S146;
        _S97 = s3_1;
        _S98 = _S145;
        _S99 = s2_2;
        _S100 = _S144;
        _S101 = _S142;
        _S102 = _S143;
        _S103 = s1_2;
        _S104 = _S139;
        _S105 = _S141;
        _S106 = s_2;
        _S107 = _S140;
        _S108 = z_1;
        _S109 = x_10;
        _S110 = y_5;
        _S111 = _S138;
        _S112 = _S137;
        _S124 = _S135;
        _S113 = _S133;
        _S114 = _S134;
        _S115 = _S130;
        _S116 = _S132;
        _S117 = _S131;
        _S118 = _S128;
        _S119 = _S126;
        _S120 = _S127;
        _S121 = _S129;
    }
    else
    {
        float2  _S149 = float2 {_S66.primal_0.x, _S66.primal_0.y};
        float2  _S150 = s_primal_ctx_exp_1(_S149);
        float _S151 = _S150.x;
        float _S152 = _S150.y;
        float _S153 = s_primal_ctx_max_0(_S151, _S152);
        float _S154 = s_primal_ctx_min_0(_S151, _S152);
        float _S155 = _S153 / _S154;
        float _S156 = _S154 * _S154;
        float2  _S157 = make_float2 (2.0f) * _S149;
        float2  _S158 = s_primal_ctx_exp_1(_S157);
        float x_11 = _S158.x;
        float y_6 = _S158.y;
        float s_3 = x_11 + y_6;
        float _S159 = s_primal_ctx_max_0(x_11, y_6);
        float s1_3 = _S159 / s_3;
        float _S160 = s_3 * s_3;
        float _S161 = s_primal_ctx_min_0(x_11, y_6);
        float s2_3 = _S161 / s_3;
        float _S162 = - s1_3;
        float _S163 = s_primal_ctx_log_0(s1_3);
        float _S164 = s_primal_ctx_log_0(s2_3);
        float _S165 = _S162 * _S163 - s2_3 * _S164;
        float _S166 = s_primal_ctx_exp_0(_S165) - 0.99998998641967773f;
        _S74 = - s_primal_ctx_log_0(_S166);
        _S75 = _S166;
        _S76 = _S165;
        _S77 = s2_3;
        _S78 = _S164;
        _S79 = _S162;
        _S80 = _S163;
        _S81 = s1_3;
        _S82 = _S160;
        _S83 = _S161;
        _S84 = s_3;
        _S85 = x_11;
        _S86 = y_6;
        _S87 = _S159;
        _S122 = _S157;
        _S88 = _S155;
        _S89 = _S156;
        _S90 = _S153;
        _S91 = _S154;
        _S92 = _S151;
        _S93 = _S152;
        _S123 = _S149;
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
        _S108 = 0.0f;
        _S109 = 0.0f;
        _S110 = 0.0f;
        _S111 = 0.0f;
        _S112 = 0.0f;
        _S124 = _S69;
        _S113 = 0.0f;
        _S114 = 0.0f;
        _S115 = 0.0f;
        _S116 = 0.0f;
        _S117 = 0.0f;
        _S118 = 0.0f;
        _S119 = 0.0f;
        _S120 = 0.0f;
        _S121 = 0.0f;
    }
    if(is_3dgs_1)
    {
        float _S167 = erank_reg_weight_s3_1 * (*_s_dOut_1)[int(3)];
        float _S168 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S169;
        (&_S169)->primal_0 = _S94;
        (&_S169)->differential_0 = 0.0f;
        DiffPair_float_0 _S170;
        (&_S170)->primal_0 = 0.0f;
        (&_S170)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S169, &_S170, _S168);
        float _S171 = - _S169.differential_0;
        DiffPair_float_0 _S172;
        (&_S172)->primal_0 = _S95;
        (&_S172)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S172, _S171);
        DiffPair_float_0 _S173;
        (&_S173)->primal_0 = _S96;
        (&_S173)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S173, _S172.differential_0);
        float _S174 = - _S173.differential_0;
        float _S175 = _S97 * _S174;
        float _S176 = _S98 * _S174;
        DiffPair_float_0 _S177;
        (&_S177)->primal_0 = _S97;
        (&_S177)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S177, _S175);
        float _S178 = _S99 * _S174;
        float _S179 = _S100 * _S174;
        DiffPair_float_0 _S180;
        (&_S180)->primal_0 = _S99;
        (&_S180)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S180, _S178);
        float _S181 = _S101 * _S173.differential_0;
        float _S182 = _S102 * _S173.differential_0;
        DiffPair_float_0 _S183;
        (&_S183)->primal_0 = _S103;
        (&_S183)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S183, _S181);
        float _S184 = - _S182;
        float _S185 = - (_S179 + _S180.differential_0);
        float _S186 = (_S167 + _S176 + _S177.differential_0 + _S185) / _S104;
        float _S187 = _S105 * - _S186;
        float _S188 = _S106 * _S186;
        DiffPair_float_0 _S189;
        (&_S189)->primal_0 = _S107;
        (&_S189)->differential_0 = 0.0f;
        DiffPair_float_0 _S190;
        (&_S190)->primal_0 = _S108;
        (&_S190)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S189, &_S190, _S188);
        DiffPair_float_0 _S191;
        (&_S191)->primal_0 = _S109;
        (&_S191)->differential_0 = 0.0f;
        DiffPair_float_0 _S192;
        (&_S192)->primal_0 = _S110;
        (&_S192)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S191, &_S192, _S189.differential_0);
        float _S193 = (_S183.differential_0 + _S184 + _S185) / _S104;
        float _S194 = _S111 * - _S193;
        float _S195 = _S106 * _S193;
        DiffPair_float_0 _S196;
        (&_S196)->primal_0 = _S112;
        (&_S196)->differential_0 = 0.0f;
        DiffPair_float_0 _S197;
        (&_S197)->primal_0 = _S108;
        (&_S197)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S196, &_S197, _S195);
        DiffPair_float_0 _S198;
        (&_S198)->primal_0 = _S109;
        (&_S198)->differential_0 = 0.0f;
        DiffPair_float_0 _S199;
        (&_S199)->primal_0 = _S110;
        (&_S199)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S198, &_S199, _S196.differential_0);
        float _S200 = _S187 + _S194;
        float3  _S201 = make_float3 (_S191.differential_0 + _S198.differential_0 + _S200, _S192.differential_0 + _S199.differential_0 + _S200, _S190.differential_0 + _S197.differential_0 + _S200);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S202;
        (&_S202)->primal_0 = _S124;
        (&_S202)->differential_0 = _S69;
        s_bwd_prop_exp_2(&_S202, _S201);
        float3  _S203 = make_float3 (2.0f) * _S202.differential_0;
        float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S204;
        (&_S204)->primal_0 = _S113;
        (&_S204)->differential_0 = 0.0f;
        DiffPair_float_0 _S205;
        (&_S205)->primal_0 = max_gauss_ratio_1;
        (&_S205)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S204, &_S205, s_diff_scale_reg_T_0);
        float _S206 = _S204.differential_0 / _S114;
        float _S207 = _S115 * - _S206;
        float _S208 = _S116 * _S206;
        DiffPair_float_0 _S209;
        (&_S209)->primal_0 = _S117;
        (&_S209)->differential_0 = 0.0f;
        DiffPair_float_0 _S210;
        (&_S210)->primal_0 = _S118;
        (&_S210)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S209, &_S210, _S207);
        DiffPair_float_0 _S211;
        (&_S211)->primal_0 = _S119;
        (&_S211)->differential_0 = 0.0f;
        DiffPair_float_0 _S212;
        (&_S212)->primal_0 = _S120;
        (&_S212)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S211, &_S212, _S209.differential_0);
        DiffPair_float_0 _S213;
        (&_S213)->primal_0 = _S121;
        (&_S213)->differential_0 = 0.0f;
        DiffPair_float_0 _S214;
        (&_S214)->primal_0 = _S118;
        (&_S214)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S213, &_S214, _S208);
        DiffPair_float_0 _S215;
        (&_S215)->primal_0 = _S119;
        (&_S215)->differential_0 = 0.0f;
        DiffPair_float_0 _S216;
        (&_S216)->primal_0 = _S120;
        (&_S216)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S215, &_S216, _S213.differential_0);
        float _S217 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * (*_s_dOut_1)[int(1)]);
        float3  _S218 = make_float3 (_S211.differential_0 + _S215.differential_0 + _S217, _S212.differential_0 + _S216.differential_0 + _S217, _S210.differential_0 + _S214.differential_0 + _S217);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S219;
        (&_S219)->primal_0 = _S66.primal_0;
        (&_S219)->differential_0 = _S69;
        s_bwd_prop_exp_2(&_S219, _S218);
        float3  _S220 = _S203 + _S219.differential_0;
        _S74 = (*_s_dOut_1)[int(4)];
        _S75 = (*_s_dOut_1)[int(0)];
        _S124 = _S220;
    }
    else
    {
        float _S221 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S222;
        (&_S222)->primal_0 = _S74;
        (&_S222)->differential_0 = 0.0f;
        DiffPair_float_0 _S223;
        (&_S223)->primal_0 = 0.0f;
        (&_S223)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S222, &_S223, _S221);
        float _S224 = - _S222.differential_0;
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = _S75;
        (&_S225)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S225, _S224);
        DiffPair_float_0 _S226;
        (&_S226)->primal_0 = _S76;
        (&_S226)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S226, _S225.differential_0);
        float _S227 = - _S226.differential_0;
        float _S228 = _S77 * _S227;
        float _S229 = _S78 * _S227;
        DiffPair_float_0 _S230;
        (&_S230)->primal_0 = _S77;
        (&_S230)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S230, _S228);
        float _S231 = _S79 * _S226.differential_0;
        float _S232 = _S80 * _S226.differential_0;
        DiffPair_float_0 _S233;
        (&_S233)->primal_0 = _S81;
        (&_S233)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S233, _S231);
        float _S234 = - _S232;
        float _S235 = (_S229 + _S230.differential_0) / _S82;
        float _S236 = _S83 * - _S235;
        float _S237 = _S84 * _S235;
        DiffPair_float_0 _S238;
        (&_S238)->primal_0 = _S85;
        (&_S238)->differential_0 = 0.0f;
        DiffPair_float_0 _S239;
        (&_S239)->primal_0 = _S86;
        (&_S239)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S238, &_S239, _S237);
        float _S240 = (_S233.differential_0 + _S234) / _S82;
        float _S241 = _S87 * - _S240;
        float _S242 = _S84 * _S240;
        DiffPair_float_0 _S243;
        (&_S243)->primal_0 = _S85;
        (&_S243)->differential_0 = 0.0f;
        DiffPair_float_0 _S244;
        (&_S244)->primal_0 = _S86;
        (&_S244)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S243, &_S244, _S242);
        float _S245 = _S236 + _S241;
        float2  _S246 = make_float2 (_S238.differential_0 + _S243.differential_0 + _S245, _S239.differential_0 + _S244.differential_0 + _S245);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S247;
        (&_S247)->primal_0 = _S122;
        (&_S247)->differential_0 = _S68;
        s_bwd_prop_exp_1(&_S247, _S246);
        float2  _S248 = make_float2 (2.0f) * _S247.differential_0;
        float s_diff_scale_reg_T_1 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S249;
        (&_S249)->primal_0 = _S88;
        (&_S249)->differential_0 = 0.0f;
        DiffPair_float_0 _S250;
        (&_S250)->primal_0 = max_gauss_ratio_1;
        (&_S250)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S249, &_S250, s_diff_scale_reg_T_1);
        float _S251 = _S249.differential_0 / _S89;
        float _S252 = _S90 * - _S251;
        float _S253 = _S91 * _S251;
        DiffPair_float_0 _S254;
        (&_S254)->primal_0 = _S92;
        (&_S254)->differential_0 = 0.0f;
        DiffPair_float_0 _S255;
        (&_S255)->primal_0 = _S93;
        (&_S255)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S254, &_S255, _S252);
        DiffPair_float_0 _S256;
        (&_S256)->primal_0 = _S92;
        (&_S256)->differential_0 = 0.0f;
        DiffPair_float_0 _S257;
        (&_S257)->primal_0 = _S93;
        (&_S257)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S256, &_S257, _S253);
        float _S258 = mcmc_scale_reg_weight_1 * (0.5f * (*_s_dOut_1)[int(1)]);
        float2  _S259 = make_float2 (_S254.differential_0 + _S256.differential_0 + _S258, _S255.differential_0 + _S257.differential_0 + _S258);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S260;
        (&_S260)->primal_0 = _S123;
        (&_S260)->differential_0 = _S68;
        s_bwd_prop_exp_1(&_S260, _S259);
        float2  _S261 = _S248 + _S260.differential_0;
        float3  _S262 = make_float3 (_S261.x, _S261.y, 0.0f);
        _S74 = (*_s_dOut_1)[int(4)];
        _S75 = (*_s_dOut_1)[int(0)];
        _S124 = _S262;
    }
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * _S74;
    float _S263 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S264;
    (&_S264)->primal_0 = _S73;
    (&_S264)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S264, _S263);
    float _S265 = _S264.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S266 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S267;
    (&_S267)->primal_0 = _S67.primal_0;
    (&_S267)->differential_0 = _S266;
    s_bwd_length_impl_0(&_S267, _S265);
    float _S268 = - (mcmc_opacity_reg_weight_1 * _S75 / _S72);
    DiffPair_float_0 _S269;
    (&_S269)->primal_0 = _S70;
    (&_S269)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S269, _S268);
    float _S270 = - _S269.differential_0;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S267.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S270;
    dpscales_0->primal_0 = (*dpscales_0).primal_0;
    dpscales_0->differential_0 = _S124;
    return;
}

inline __device__ void s_bwd_per_splat_losses_0(bool _S271, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S272, DiffPair_float_0 * _S273, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S274, float _S275, float _S276, float _S277, float _S278, float _S279, float _S280, float _S281, FixedArray<float, 5>  * _S282)
{
    s_bwd_prop_per_splat_losses_0(_S271, _S272, _S273, _S274, _S275, _S276, _S277, _S278, _S279, _S280, _S281, _S282);
    return;
}

inline __device__ void per_splat_losses_bwd(bool is_3dgs_2, float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  * v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float3  _S283 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_scales_0;
    (&p_scales_0)->primal_0 = scales_1;
    (&p_scales_0)->differential_0 = _S283;
    DiffPair_float_0 p_opacity_0;
    (&p_opacity_0)->primal_0 = opacity_1;
    (&p_opacity_0)->differential_0 = 0.0f;
    float4  _S284 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 p_quat_0;
    (&p_quat_0)->primal_0 = quat_1;
    (&p_quat_0)->differential_0 = _S284;
    s_bwd_per_splat_losses_0(is_3dgs_2, &p_scales_0, &p_opacity_0, &p_quat_0, mcmc_opacity_reg_weight_2, mcmc_scale_reg_weight_2, max_gauss_ratio_2, scale_regularization_weight_2, erank_reg_weight_2, erank_reg_weight_s3_2, quat_norm_reg_weight_2, v_loss_0);
    *v_scales_0 = p_scales_0.differential_0;
    *v_opacity_0 = p_opacity_0.differential_0;
    *v_quat_0 = p_quat_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_12)
{
    Matrix<float, 3, 3>  result_12;
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
            *_slang_vector_get_element_ptr(((&result_12)->rows + (r_1)), c_0) = _slang_vector_get_element(x_12.rows[c_0], r_1);
            c_0 = c_0 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_12;
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
    float _S285 = (*left_0).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_8.x;
    float sum_0 = _S285 + (*left_0).primal_0.rows[int(1)].x * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_8.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_8.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S286 = (*left_0).primal_0.rows[int(0)].y * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_8.x;
    float sum_2 = _S286 + (*left_0).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_8.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_8.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_8.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S287 = (*left_0).primal_0.rows[int(0)].z * dOut_8.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_8.x;
    float sum_4 = _S287 + (*left_0).primal_0.rows[int(1)].z * dOut_8.y;
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
    float _S288 = (*left_1).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 2, 2>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_1).primal_0.x * dOut_9.x;
    float sum_6 = _S288 + (*left_1).primal_0.rows[int(1)].x * dOut_9.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_1).primal_0.x * dOut_9.y;
    float2  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_6;
    float _S289 = (*left_1).primal_0.rows[int(0)].y * dOut_9.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_1).primal_0.y * dOut_9.x;
    float sum_7 = _S289 + (*left_1).primal_0.rows[int(1)].y * dOut_9.y;
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
    float3  result_13;
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
        *_slang_vector_get_element_ptr(&result_13, i_5) = sum_8;
        i_5 = i_5 + int(1);
    }
    return result_13;
}

inline __device__ float2  mul_1(Matrix<float, 2, 2>  left_3, float2  right_3)
{
    float2  result_14;
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
        *_slang_vector_get_element_ptr(&result_14, i_6) = sum_10;
        i_6 = i_6 + int(1);
    }
    return result_14;
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
    Matrix<float, 3, 3>  result_15;
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
            *_slang_vector_get_element_ptr(((&result_15)->rows + (r_2)), c_1) = sum_12;
            c_1 = c_1 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_15;
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
    float _S290 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S291 = (*dist_coeffs_0)[int(1)] + r2_0 * _S290;
    float _S292 = (*dist_coeffs_0)[int(0)] + r2_0 * _S291;
    float2  _S293 = make_float2 (1.0f + r2_0 * _S292);
    float _S294 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S295 = _S294 * u_0;
    float _S296 = 2.0f * u_0;
    float _S297 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S298 = _S297 * u_0;
    float _S299 = 2.0f * v_0;
    float2  _S300 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S301 = uv_0 * _S300;
    float _S302 = (*dist_coeffs_0)[int(4)] * _S300.y;
    float _S303 = (*dist_coeffs_0)[int(5)] * _S300.x;
    float _S304 = _S301.x + _S301.y;
    float _S305 = r2_0 * _S304;
    float _S306 = r2_0 * _S305;
    float _S307 = (*dist_coeffs_0)[int(7)] * _S300.y + _S302 + (*dist_coeffs_0)[int(6)] * _S300.x + _S303 + _S292 * _S304 + _S291 * _S305 + _S290 * _S306 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S306);
    float _S308 = v_0 * _S307;
    float _S309 = u_0 * _S307;
    float2  _S310 = make_float2 (0.0f, 1.0f) + make_float2 ((*dist_coeffs_0)[int(8)] * 0.0f, (*dist_coeffs_0)[int(9)] * 0.0f);
    float2  _S311 = uv_0 * _S310;
    float _S312 = (*dist_coeffs_0)[int(4)] * _S310.y;
    float _S313 = (*dist_coeffs_0)[int(5)] * _S310.x;
    float _S314 = _S311.x + _S311.y;
    float _S315 = r2_0 * _S314;
    float _S316 = r2_0 * _S315;
    float _S317 = (*dist_coeffs_0)[int(7)] * _S310.y + _S312 + (*dist_coeffs_0)[int(6)] * _S310.x + _S313 + _S292 * _S314 + _S291 * _S315 + _S290 * _S316 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S316);
    float _S318 = v_0 * _S317;
    float _S319 = u_0 * _S317;
    return makeMatrix<float, 2, 2> (_S293 * _S300 + make_float2 (_S297 * (v_0 * _S300.y) + _S296 * _S303 + 2.0f * (u_0 * _S303) + _S294 * (v_0 * _S300.x) + _S309 + _S309, _S299 * _S302 + 2.0f * (v_0 * _S302) + _S298 * _S300.y + _S295 * _S300.x + _S308 + _S308), _S293 * _S310 + make_float2 (_S297 * (v_0 * _S310.y) + _S296 * _S313 + 2.0f * (u_0 * _S313) + _S294 * (v_0 * _S310.x) + _S319 + _S319, _S299 * _S312 + 2.0f * (v_0 * _S312) + _S298 * _S310.y + _S295 * _S310.x + _S318 + _S318));
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  * dist_coeffs_1)
{
    Matrix<float, 2, 2>  _S320 = camera_distortion_jac_0(uv_1, dist_coeffs_1);
    return (F32_min((determinant_0(_S320)), ((F32_min((_S320.rows[int(0)].x), (_S320.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_2)
{
    float2  _S321;
    if(is_fisheye_0)
    {
        float r_3 = length_1(uv_2);
        float theta_0 = (F32_atan((r_3)));
        float _S322;
        if(r_3 < 0.00100000004749745f)
        {
            _S322 = 1.0f - r_3 * r_3 / 3.0f;
        }
        else
        {
            _S322 = theta_0 / r_3;
        }
        _S321 = uv_2 * make_float2 (_S322);
    }
    else
    {
        _S321 = uv_2;
    }
    float u_1 = _S321.x;
    float v_1 = _S321.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S323 = _S321 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_2)[int(0)] + r2_1 * ((*dist_coeffs_2)[int(1)] + r2_1 * ((*dist_coeffs_2)[int(2)] + r2_1 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_1 * v_1 + (*dist_coeffs_2)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_2)[int(6)] * r2_1, 2.0f * (*dist_coeffs_2)[int(5)] * u_1 * v_1 + (*dist_coeffs_2)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_2)[int(7)] * r2_1);
    return _S323 + make_float2 ((*dist_coeffs_2)[int(8)] * _S323.x + (*dist_coeffs_2)[int(9)] * _S323.y, 0.0f);
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
        float2  _S324 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_4 = _S324 + make_float2 ((*dist_coeffs_3)[int(8)] * _S324.x + (*dist_coeffs_3)[int(9)] * _S324.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S325 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S325.rows[int(0)].x * _S325.rows[int(1)].y - _S325.rows[int(0)].y * _S325.rows[int(1)].x);
        float _S326 = r_4.x;
        float _S327 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S326 * _S325.rows[int(1)].y - _S327 * _S325.rows[int(0)].y) * inv_det_0, (- _S326 * _S325.rows[int(1)].x + _S327 * _S325.rows[int(0)].x) * inv_det_0);
        i_8 = i_8 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S328 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S329;
    if((F32_min((determinant_0(_S328)), ((F32_min((_S328.rows[int(0)].x), (_S328.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S330 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S329 = (length_1(_S330 + make_float2 ((*dist_coeffs_3)[int(8)] * _S330.x + (*dist_coeffs_3)[int(9)] * _S330.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S329 = false;
    }
    return _S329;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S331 = uv_4;
    bool _S332 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S331);
    if(!_S332)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S333 = _S331;
        float theta_1 = length_1(_S331);
        float _S334;
        if(theta_1 < 0.00100000004749745f)
        {
            _S334 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S334 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S335 = make_float3 ((_S333 * make_float2 (_S334)).x, (_S333 * make_float2 (_S334)).y, (F32_cos((theta_1))));
        raydir_0 = _S335;
    }
    else
    {
        raydir_0 = make_float3 (_S331.x, _S331.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_5, float3  * raydir_1)
{
    float2  _S336 = uv_5;
    int3  _S337 = make_int3 (int(0));
    float3  _S338 = make_float3 ((float)_S337.x, (float)_S337.y, (float)_S337.z);
    *raydir_1 = _S338;
    bool _S339 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S336);
    if(!_S339)
    {
        return false;
    }
    float3  _S340;
    if(is_fisheye_2)
    {
        float2  _S341 = _S336;
        float theta_2 = length_1(_S336);
        float _S342;
        if(theta_2 < 0.00100000004749745f)
        {
            _S342 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S342 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S343 = make_float3 ((_S341 * make_float2 (_S342)).x, (_S341 * make_float2 (_S342)).y, (F32_cos((theta_2))));
        _S340 = _S343;
    }
    else
    {
        _S340 = make_float3 (_S336.x, _S336.y, 1.0f);
    }
    *raydir_1 = _S340;
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
    float2  _S344 = uv_6;
    bool _S345 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S344);
    if(!_S345)
    {
        int3  _S346 = make_int3 (int(0));
        float3  _S347 = make_float3 ((float)_S346.x, (float)_S346.y, (float)_S346.z);
        *raydir_2 = _S347;
        return false;
    }
    float3  _S348;
    if(is_fisheye_3)
    {
        float2  _S349 = _S344;
        float theta_3 = length_1(_S344);
        float _S350;
        if(theta_3 < 0.00100000004749745f)
        {
            _S350 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S350 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S351 = make_float3 ((_S349 * make_float2 (_S350)).x, (_S349 * make_float2 (_S350)).y, (F32_cos((theta_3))));
        _S348 = _S351;
    }
    else
    {
        _S348 = make_float3 (_S344.x, _S344.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S348);
    return true;
}

inline __device__ void _d_mul_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_6, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_6, float3  dOut_11)
{
    float _S352 = (*right_6).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 3, 3>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = (*left_6).primal_0.x * dOut_11.x;
    float sum_14 = _S352 + (*right_6).primal_0.rows[int(0)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = (*left_6).primal_0.x * dOut_11.y;
    float sum_15 = sum_14 + (*right_6).primal_0.rows[int(0)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(0)))->z) = (*left_6).primal_0.x * dOut_11.z;
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = sum_15;
    float _S353 = (*right_6).primal_0.rows[int(1)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = (*left_6).primal_0.y * dOut_11.x;
    float sum_16 = _S353 + (*right_6).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = (*left_6).primal_0.y * dOut_11.y;
    float sum_17 = sum_16 + (*right_6).primal_0.rows[int(1)].z * dOut_11.z;
    *&(((&right_d_result_3)->rows + (int(1)))->z) = (*left_6).primal_0.y * dOut_11.z;
    *&((&left_d_result_3)->y) = sum_17;
    float _S354 = (*right_6).primal_0.rows[int(2)].x * dOut_11.x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = (*left_6).primal_0.z * dOut_11.x;
    float sum_18 = _S354 + (*right_6).primal_0.rows[int(2)].y * dOut_11.y;
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
    float3  result_16;
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
        *_slang_vector_get_element_ptr(&result_16, j_2) = sum_20;
        j_2 = j_2 + int(1);
    }
    return result_16;
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

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S355, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S356, float3  _S357)
{
    _d_mul_2(_S355, _S356, _S357);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S358 = - _s_dOut_2;
    float3  _S359 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S360;
    (&_S360)->primal_0 = (*dpt_0).primal_0;
    (&_S360)->differential_0 = _S359;
    Matrix<float, 3, 3>  _S361 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S362;
    (&_S362)->primal_0 = (*dpR_0).primal_0;
    (&_S362)->differential_0 = _S361;
    s_bwd_prop_mul_0(&_S360, &_S362, _S358);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S360.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S362.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S363, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S364, float3  _S365)
{
    s_bwd_prop_transform_ray_o_0(_S363, _S364, _S365);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_5, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S366 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_5;
    (&dp_R_0)->differential_0 = _S366;
    float3  _S367 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S367;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S368 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S369;
    (&_S369)->primal_0 = (*dpraydir_0).primal_0;
    (&_S369)->differential_0 = _S368;
    Matrix<float, 3, 3>  _S370 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S371;
    (&_S371)->primal_0 = (*dpR_1).primal_0;
    (&_S371)->differential_0 = _S370;
    s_bwd_prop_mul_0(&_S369, &_S371, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S369.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S371.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S372, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S373, float3  _S374)
{
    s_bwd_prop_transform_ray_d_0(_S372, _S373, _S374);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_6, float3  raydir_5, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S375 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_6;
    (&dp_R_1)->differential_0 = _S375;
    float3  _S376 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_5;
    (&dp_raydir_0)->differential_0 = _S376;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S377 = scale_2.x;
    float sx_0 = (F32_exp((_S377)));
    float _S378 = scale_2.y;
    float sy_0 = (F32_exp((_S378)));
    float sz_0 = scale_2.z - 0.5f * (_S377 + _S378);
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
    Matrix<float, 3, 3>  _S379 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S379, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S379, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S379, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
    return;
}

inline __device__ float4  floor_0(float4  x_19)
{
    float4  result_17;
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
        *_slang_vector_get_element_ptr(&result_17, i_10) = (F32_floor((_slang_vector_get_element(x_19, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_17;
}

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_1, float3  scale_3, float4  quat_6, float opac_0)
{
    float4  _S380 = normalize_1(quat_6);
    float3  _S381 = exp_0(scale_3);
    float x_20 = _S380.y;
    float x2_4 = x_20 * x_20;
    float y2_4 = _S380.z * _S380.z;
    float z2_4 = _S380.w * _S380.w;
    float xy_4 = _S380.y * _S380.z;
    float xz_4 = _S380.y * _S380.w;
    float yz_4 = _S380.z * _S380.w;
    float wx_4 = _S380.x * _S380.y;
    float wy_4 = _S380.x * _S380.z;
    float wz_4 = _S380.x * _S380.w;
    Matrix<float, 3, 3>  M_2 = mul_3(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S381.x, 0.0f, 0.0f, 0.0f, _S381.y, 0.0f, 0.0f, 0.0f, _S381.z));
    float4  _S382 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S383 = _S382 - floor_0(_S382);
    float4  _S384 = _S383 + make_float4 (dot_1(_S383, float4 {_S383.w, _S383.z, _S383.x, _S383.y} + make_float4 (33.3300018310546875f)));
    float4  _S385 = (float4 {_S384.x, _S384.x, _S384.y, _S384.z} + float4 {_S384.y, _S384.z, _S384.z, _S384.w}) * float4 {_S384.z, _S384.y, _S384.w, _S384.x};
    float4  _S386 = _S385 - floor_0(_S385);
    float2  _S387 = float2 {_S386.x, _S386.z};
    float _S388 = 6.28318548202514648f * _S387.y;
    float2  _S389 = float2 {_S386.y, _S386.w};
    float _S390 = 6.28318548202514648f * _S389.y;
    *mean_1 = *mean_1 + mul_0(mul_3(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S387.x))))))) * make_float2 ((F32_cos((_S388))), (F32_sin((_S388))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S387.x))))))) * make_float2 ((F32_cos((_S388))), (F32_sin((_S388))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S389.x))))))) * make_float2 ((F32_cos((_S390))), (F32_sin((_S390))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_12)
{
    float _S391 = dOut_12.y;
    float _S392 = dOut_12.z;
    float _S393 = dOut_12.x;
    float _S394 = (*a_0).primal_0.z * _S391 + - (*a_0).primal_0.y * _S392;
    float _S395 = - (*a_0).primal_0.z * _S393 + (*a_0).primal_0.x * _S392;
    float _S396 = (*a_0).primal_0.y * _S393 + - (*a_0).primal_0.x * _S391;
    float3  _S397 = make_float3 (- (*b_1).primal_0.z * _S391 + (*b_1).primal_0.y * _S392, (*b_1).primal_0.z * _S393 + - (*b_1).primal_0.x * _S392, - (*b_1).primal_0.y * _S393 + (*b_1).primal_0.x * _S391);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S397;
    float3  _S398 = make_float3 (_S394, _S395, _S396);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S398;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S399 = left_8.y;
    float _S400 = right_8.z;
    float _S401 = left_8.z;
    float _S402 = right_8.y;
    float _S403 = right_8.x;
    float _S404 = left_8.x;
    return make_float3 (_S399 * _S400 - _S401 * _S402, _S401 * _S403 - _S404 * _S400, _S404 * _S402 - _S399 * _S403);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S405 = normalize_1(quat_7);
    float _S406 = scale_4.x;
    float sx_1 = (F32_exp((_S406)));
    float _S407 = scale_4.y;
    float sy_1 = (F32_exp((_S407)));
    float sz_1 = scale_4.z - 0.5f * (_S406 + _S407);
    float x_21 = _S405.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S405.z * _S405.z;
    float z2_5 = _S405.w * _S405.w;
    float xy_5 = _S405.y * _S405.z;
    float xz_5 = _S405.y * _S405.w;
    float yz_5 = _S405.z * _S405.w;
    float wx_5 = _S405.x * _S405.y;
    float wy_5 = _S405.x * _S405.z;
    float wz_5 = _S405.x * _S405.w;
    Matrix<float, 3, 3>  _S408 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S408, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S408, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S408, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S409 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S410 = _S409 - floor_0(_S409);
    float4  _S411 = _S410 + make_float4 (dot_1(_S410, float4 {_S410.w, _S410.z, _S410.x, _S410.y} + make_float4 (33.3300018310546875f)));
    float4  _S412 = (float4 {_S411.x, _S411.x, _S411.y, _S411.z} + float4 {_S411.y, _S411.z, _S411.z, _S411.w}) * float4 {_S411.z, _S411.y, _S411.w, _S411.x};
    float4  _S413 = _S412 - floor_0(_S412);
    float2  _S414 = float2 {_S413.x, _S413.z};
    float _S415 = 6.28318548202514648f * _S414.y;
    float2  _S416 = float2 {_S413.y, _S413.w};
    float _S417 = 6.28318548202514648f * _S416.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S414.x))))))) * make_float2 ((F32_cos((_S415))), (F32_sin((_S415))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S414.x))))))) * make_float2 ((F32_cos((_S415))), (F32_sin((_S415))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S416.x))))))) * make_float2 ((F32_cos((_S417))), (F32_sin((_S417))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_9, float dOut_13)
{
    float _S418 = _slang_select(((*dpx_9).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_9).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S418;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S419 = _slang_select(((*dpx_10).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_10).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S419;
    return;
}

inline __device__ float3  abs_0(float3  x_22)
{
    float3  result_18;
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
        *_slang_vector_get_element_ptr(&result_18, i_11) = (F32_abs((_slang_vector_get_element(x_22, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_18;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_15)
{
    DiffPair_float_0 _S420 = *dpx_11;
    bool _S421;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S421 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S421 = false;
    }
    float _S422;
    if(_S421)
    {
        _S422 = dOut_15;
    }
    else
    {
        _S422 = 0.0f;
    }
    dpx_11->primal_0 = _S420.primal_0;
    dpx_11->differential_0 = _S422;
    DiffPair_float_0 _S423 = *dpMin_0;
    if((_S420.primal_0) < ((*dpMin_0).primal_0))
    {
        _S422 = dOut_15;
    }
    else
    {
        _S422 = 0.0f;
    }
    dpMin_0->primal_0 = _S423.primal_0;
    dpMin_0->differential_0 = _S422;
    DiffPair_float_0 _S424 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S422 = dOut_15;
    }
    else
    {
        _S422 = 0.0f;
    }
    dpMax_0->primal_0 = _S424.primal_0;
    dpMax_0->differential_0 = _S422;
    return;
}

inline __device__ float clamp_0(float x_23, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_23), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_12, float dOut_16)
{
    float _S425 = -0.5f / ((*dpx_12).primal_0 * (F32_sqrt(((*dpx_12).primal_0)))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S425;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_17)
{
    float _S426 = (1.0f - (*dps_0).primal_0) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S426;
    DiffPair_float_0 _S427 = *dpy_3;
    float _S428 = (*dps_0).primal_0 * dOut_17;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S428;
    float _S429 = (_S427.primal_0 - (*dpx_13).primal_0) * dOut_17;
    dps_0->primal_0 = _S427.primal_0;
    dps_0->differential_0 = _S429;
    return;
}

inline __device__ float lerp_0(float x_24, float y_7, float s_4)
{
    return x_24 + (y_7 - x_24) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  * weights_0, FixedArray<float, 23>  * _S430)
{
    float3  _S431;
    bool _S432;
    bool _S433;
    FixedArray<float, 23>  losses_1;
    float _S434 = float(mask_0);
    float3  _S435 = ref_rgb_0 - render_rgb_0;
    float3  _S436 = abs_0(_S435);
    losses_1[int(0)] = (*weights_0)[int(0)] * _S434 * ((_S436.x + _S436.y + _S436.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S434 * clamp_0(dot_0(_S435, _S435) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S437 = float(depth_mask_0 & mask_0);
    float _S438 = _S437 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S439 = _S437 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S438;
    losses_1[int(3)] = _S439;
    losses_1[int(4)] = _S438 * _S438;
    losses_1[int(5)] = _S439 * _S439;
    losses_1[int(6)] = _S438 * _S439;
    bool _S440 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S441 = norm2_0 == 0.0f;
        _S432 = _S441;
        if(_S441)
        {
            _S431 = make_float3 (0.0f);
            break;
        }
        _S431 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S442;
    bool _S443 = !_S432;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S444 = norm2_1 == 0.0f;
        _S433 = _S444;
        if(_S444)
        {
            _S442 = make_float3 (0.0f);
            break;
        }
        _S442 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S445;
    float3  _S446;
    bool _S447 = !_S433;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S446 = make_float3 (0.0f);
            _S445 = false;
            break;
        }
        _S446 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S445 = _S440;
        break;
    }
    float _S448 = float(_S443 & _S445);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S431, _S446);
    losses_1[int(7)] = (*weights_0)[int(2)] * _S448 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S449 = float(_S447 & _S445);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S442, _S446);
    losses_1[int(8)] = (*weights_0)[int(2)] * _S449 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S450 = float(_S443 & _S447);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S431, _S442);
    losses_1[int(11)] = (*weights_0)[int(5)] * _S450 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S451 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S452 = float(alpha_mask_0);
    float _S453 = float(ref_alpha_0);
    float _S454 = (F32_max((_S451), (_S453)));
    losses_1[int(9)] = (*weights_0)[int(3)] * _S452 * - lerp_0((F32_log(((F32_max((1.0f - _S454), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S454), (9.99999997475242708e-07f)))))), _S453);
    float _S455 = 1.0f - _S451;
    float _S456 = 1.0f - _S453;
    float _S457 = (F32_max((_S455), (_S456)));
    losses_1[int(10)] = (*weights_0)[int(4)] * _S452 * - lerp_0((F32_log(((F32_max((1.0f - _S457), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S457), (9.99999997475242708e-07f)))))), _S456);
    losses_1[int(12)] = (*weights_0)[int(6)] * 4.0f * _S451 * _S455;
    float _S458 = (F32_max((_S451), (9.999999960041972e-13f)));
    losses_1[int(13)] = (*weights_0)[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S458;
    losses_1[int(14)] = (*weights_0)[int(8)] * depth_dist_0 / _S458;
    losses_1[int(15)] = (*weights_0)[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S458;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S434;
    losses_1[int(18)] = _S437;
    losses_1[int(19)] = _S448;
    losses_1[int(20)] = _S449;
    losses_1[int(21)] = _S450;
    losses_1[int(22)] = _S452;
    *_S430 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S459, float3  _S460)
{
    return dot_0(_S459, _S460);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S461)
{
    return (F32_rsqrt((_S461)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S462, float _S463, float _S464)
{
    return clamp_0(_S462, _S463, _S464);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S465, DiffPair_float_0 * _S466, DiffPair_float_0 * _S467, float _S468)
{
    _d_lerp_0(_S465, _S466, _S467, _S468);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S469, DiffPair_float_0 * _S470, DiffPair_float_0 * _S471, float _S472)
{
    _d_clamp_0(_S469, _S470, _S471, _S472);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S473, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S474, float _S475)
{
    _d_dot_0(_S473, _S474, _S475);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S476, float _S477)
{
    _d_rsqrt_0(_S476, _S477);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S478, float3  _S479)
{
    _d_abs_vector_0(_S478, _S479);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_float_0 _S480 = *dprender_depth_0;
    DiffPair_float_0 _S481 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S482 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S483 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484 = *dpref_normal_0;
    DiffPair_float_0 _S485 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S486 = *dprgb_dist_0;
    DiffPair_float_0 _S487 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S488 = *dpnormal_dist_0;
    float3  _S489 = make_float3 (0.0f);
    float _S490 = float(mask_1);
    float _S491 = (*weights_1)[int(0)] * _S490;
    float3  _S492 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S493 = s_primal_ctx_dot_0(_S492, _S492) * 0.3333333432674408f;
    float _S494 = float(depth_mask_1 & mask_1);
    float _S495 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S496 = _S494 * s_primal_ctx_log_0(_S495);
    float _S497 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S498 = _S494 * s_primal_ctx_log_0(_S497);
    bool _S499 = normal_mask_1 & mask_1;
    float _S500 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S501 = _S500 == 0.0f;
    float3  _S502;
    if(_S501)
    {
        _S502 = make_float3 (0.0f);
    }
    bool _S503 = !_S501;
    float3  _S504;
    if(_S503)
    {
        float _S505 = s_primal_ctx_rsqrt_0(_S500);
        float3  _S506 = make_float3 (_S505);
        _S502 = _S482.primal_0 * make_float3 (_S505);
        _S504 = _S506;
    }
    else
    {
        _S504 = _S489;
    }
    float _S507 = s_primal_ctx_dot_0(_S483.primal_0, _S483.primal_0);
    bool _S508 = _S507 == 0.0f;
    float3  _S509;
    if(_S508)
    {
        _S509 = make_float3 (0.0f);
    }
    bool _S510 = !_S508;
    float3  _S511;
    if(_S510)
    {
        float _S512 = s_primal_ctx_rsqrt_0(_S507);
        float3  _S513 = make_float3 (_S512);
        _S509 = _S483.primal_0 * make_float3 (_S512);
        _S511 = _S513;
    }
    else
    {
        _S511 = _S489;
    }
    float _S514 = s_primal_ctx_dot_0(_S484.primal_0, _S484.primal_0);
    bool _S515 = _S514 == 0.0f;
    float3  _S516;
    bool _S517;
    if(_S515)
    {
        float3  _S518 = make_float3 (0.0f);
        _S517 = false;
        _S516 = _S518;
    }
    else
    {
        _S517 = _S499;
    }
    bool _S519 = !_S515;
    float3  _S520;
    if(_S519)
    {
        float _S521 = s_primal_ctx_rsqrt_0(_S514);
        float3  _S522 = make_float3 (_S521);
        _S516 = _S484.primal_0 * make_float3 (_S521);
        _S520 = _S522;
    }
    else
    {
        _S520 = _S489;
    }
    float _S523 = (*weights_1)[int(2)] * float(_S503 & _S517);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S502, _S516);
    float _S524 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S525 = (*weights_1)[int(2)] * float(_S510 & _S517);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S509, _S516);
    float _S526 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S527 = (*weights_1)[int(5)] * float(_S503 & _S510);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S502, _S509);
    float _S528 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S529 = s_primal_ctx_clamp_0(_S485.primal_0, 0.0f, 1.0f);
    float _S530 = float(alpha_mask_1);
    float _S531 = (*weights_1)[int(3)] * _S530;
    float _S532 = float(ref_alpha_1);
    float _S533 = s_primal_ctx_max_0(_S529, _S532);
    float _S534 = 1.0f - _S533;
    float _S535 = s_primal_ctx_max_0(_S534, 9.99999997475242708e-07f);
    float _S536 = s_primal_ctx_log_0(_S535);
    float _S537 = s_primal_ctx_max_0(_S533, 9.99999997475242708e-07f);
    float _S538 = s_primal_ctx_log_0(_S537);
    float _S539 = (*weights_1)[int(4)] * _S530;
    float _S540 = 1.0f - _S529;
    float _S541 = 1.0f - _S532;
    float _S542 = s_primal_ctx_max_0(_S540, _S541);
    float _S543 = 1.0f - _S542;
    float _S544 = s_primal_ctx_max_0(_S543, 9.99999997475242708e-07f);
    float _S545 = s_primal_ctx_log_0(_S544);
    float _S546 = s_primal_ctx_max_0(_S542, 9.99999997475242708e-07f);
    float _S547 = s_primal_ctx_log_0(_S546);
    float _S548 = (*weights_1)[int(6)] * 4.0f;
    float _S549 = _S548 * _S529;
    float _S550 = s_primal_ctx_max_0(_S529, 9.999999960041972e-13f);
    float _S551 = _S550 * _S550;
    float _S552 = (*_s_dOut_4)[int(15)] / _S551;
    float _S553 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S550 * _S552));
    float _S554 = (*_s_dOut_4)[int(14)] / _S551;
    float _S555 = (*weights_1)[int(8)] * (_S550 * _S554);
    float _S556 = (*_s_dOut_4)[int(13)] / _S551;
    float _S557 = _S550 * _S556;
    float _S558 = (*weights_1)[int(9)] * ((_S488.primal_0.x + _S488.primal_0.y + _S488.primal_0.z) * 0.3333333432674408f) * - _S552 + (*weights_1)[int(8)] * _S487.primal_0 * - _S554 + (*weights_1)[int(7)] * ((_S486.primal_0.x + _S486.primal_0.y + _S486.primal_0.z) * 0.3333333432674408f) * - _S556;
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = _S529;
    (&_S559)->differential_0 = 0.0f;
    DiffPair_float_0 _S560;
    (&_S560)->primal_0 = 9.999999960041972e-13f;
    (&_S560)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S559, &_S560, _S558);
    float _S561 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S557);
    float _S562 = _S549 * (*_s_dOut_4)[int(12)];
    float _S563 = _S548 * (_S540 * (*_s_dOut_4)[int(12)]);
    float _S564 = - (_S539 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S565;
    (&_S565)->primal_0 = _S545;
    (&_S565)->differential_0 = 0.0f;
    DiffPair_float_0 _S566;
    (&_S566)->primal_0 = _S547;
    (&_S566)->differential_0 = 0.0f;
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S541;
    (&_S567)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S565, &_S566, &_S567, _S564);
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S546;
    (&_S568)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S568, _S566.differential_0);
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S542;
    (&_S569)->differential_0 = 0.0f;
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = 9.99999997475242708e-07f;
    (&_S570)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S569, &_S570, _S568.differential_0);
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S544;
    (&_S571)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S571, _S565.differential_0);
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S543;
    (&_S572)->differential_0 = 0.0f;
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = 9.99999997475242708e-07f;
    (&_S573)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S572, &_S573, _S571.differential_0);
    float _S574 = _S569.differential_0 + - _S572.differential_0;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = _S540;
    (&_S575)->differential_0 = 0.0f;
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S541;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S575, &_S576, _S574);
    float _S577 = - (_S562 + _S575.differential_0);
    float _S578 = - (_S531 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S536;
    (&_S579)->differential_0 = 0.0f;
    DiffPair_float_0 _S580;
    (&_S580)->primal_0 = _S538;
    (&_S580)->differential_0 = 0.0f;
    DiffPair_float_0 _S581;
    (&_S581)->primal_0 = _S532;
    (&_S581)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S579, &_S580, &_S581, _S578);
    DiffPair_float_0 _S582;
    (&_S582)->primal_0 = _S537;
    (&_S582)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S582, _S580.differential_0);
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S533;
    (&_S583)->differential_0 = 0.0f;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = 9.99999997475242708e-07f;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S583, &_S584, _S582.differential_0);
    DiffPair_float_0 _S585;
    (&_S585)->primal_0 = _S535;
    (&_S585)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S585, _S579.differential_0);
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = _S534;
    (&_S586)->differential_0 = 0.0f;
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = 9.99999997475242708e-07f;
    (&_S587)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S586, &_S587, _S585.differential_0);
    float _S588 = _S583.differential_0 + - _S586.differential_0;
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = _S529;
    (&_S589)->differential_0 = 0.0f;
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = _S532;
    (&_S590)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S589, &_S590, _S588);
    float _S591 = _S559.differential_0 + _S563 + _S577 + _S589.differential_0;
    DiffPair_float_0 _S592;
    (&_S592)->primal_0 = _S485.primal_0;
    (&_S592)->differential_0 = 0.0f;
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = 0.0f;
    (&_S593)->differential_0 = 0.0f;
    DiffPair_float_0 _S594;
    (&_S594)->primal_0 = 1.0f;
    (&_S594)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S592, &_S593, &_S594, _S591);
    DiffPair_float_0 _S595 = _S592;
    float _S596 = _S527 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S597;
    (&_S597)->primal_0 = _S528;
    (&_S597)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S597, _S596);
    DiffPair_float_0 _S598;
    (&_S598)->primal_0 = cos_sim_loss_5;
    (&_S598)->differential_0 = 0.0f;
    DiffPair_float_0 _S599;
    (&_S599)->primal_0 = 9.999999960041972e-13f;
    (&_S599)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S598, &_S599, _S597.differential_0);
    float _S600 = 0.5f * - (_S596 + _S598.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S601;
    (&_S601)->primal_0 = _S502;
    (&_S601)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S602;
    (&_S602)->primal_0 = _S509;
    (&_S602)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S601, &_S602, _S600);
    float _S603 = _S525 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S604;
    (&_S604)->primal_0 = _S526;
    (&_S604)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S604, _S603);
    DiffPair_float_0 _S605;
    (&_S605)->primal_0 = cos_sim_loss_4;
    (&_S605)->differential_0 = 0.0f;
    DiffPair_float_0 _S606;
    (&_S606)->primal_0 = 9.999999960041972e-13f;
    (&_S606)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S605, &_S606, _S604.differential_0);
    float _S607 = 0.5f * - (_S603 + _S605.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S608;
    (&_S608)->primal_0 = _S509;
    (&_S608)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S609;
    (&_S609)->primal_0 = _S516;
    (&_S609)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S608, &_S609, _S607);
    float _S610 = _S523 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S611;
    (&_S611)->primal_0 = _S524;
    (&_S611)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S611, _S610);
    DiffPair_float_0 _S612;
    (&_S612)->primal_0 = cos_sim_loss_3;
    (&_S612)->differential_0 = 0.0f;
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = 9.999999960041972e-13f;
    (&_S613)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S612, &_S613, _S611.differential_0);
    float _S614 = 0.5f * - (_S610 + _S612.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S615;
    (&_S615)->primal_0 = _S502;
    (&_S615)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S616;
    (&_S616)->primal_0 = _S516;
    (&_S616)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S615, &_S616, _S614);
    float3  _S617 = _S609.differential_0 + _S616.differential_0;
    float3  _S618 = _S601.differential_0 + _S615.differential_0;
    float3  _S619 = make_float3 (_S553, _S553, _S553);
    float3  _S620 = make_float3 (_S561, _S561, _S561);
    float3  _S621 = _S602.differential_0 + _S608.differential_0;
    float _S622;
    if(_S519)
    {
        float3  _S623 = _S484.primal_0 * _S617;
        float3  _S624 = _S520 * _S617;
        float _S625 = _S623.x + _S623.y + _S623.z;
        DiffPair_float_0 _S626;
        (&_S626)->primal_0 = _S514;
        (&_S626)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S626, _S625);
        _S622 = _S626.differential_0;
        _S502 = _S624;
    }
    else
    {
        _S622 = 0.0f;
        _S502 = _S489;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S627;
    (&_S627)->primal_0 = _S484.primal_0;
    (&_S627)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S628;
    (&_S628)->primal_0 = _S484.primal_0;
    (&_S628)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S627, &_S628, _S622);
    float3  _S629 = _S628.differential_0 + _S627.differential_0 + _S502;
    if(_S510)
    {
        float3  _S630 = _S483.primal_0 * _S621;
        float3  _S631 = _S511 * _S621;
        float _S632 = _S630.x + _S630.y + _S630.z;
        DiffPair_float_0 _S633;
        (&_S633)->primal_0 = _S507;
        (&_S633)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S633, _S632);
        _S622 = _S633.differential_0;
        _S502 = _S631;
    }
    else
    {
        _S622 = 0.0f;
        _S502 = _S489;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S634;
    (&_S634)->primal_0 = _S483.primal_0;
    (&_S634)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S635;
    (&_S635)->primal_0 = _S483.primal_0;
    (&_S635)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S634, &_S635, _S622);
    float3  _S636 = _S635.differential_0 + _S634.differential_0 + _S502;
    if(_S503)
    {
        float3  _S637 = _S482.primal_0 * _S618;
        float3  _S638 = _S504 * _S618;
        float _S639 = _S637.x + _S637.y + _S637.z;
        DiffPair_float_0 _S640;
        (&_S640)->primal_0 = _S500;
        (&_S640)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S640, _S639);
        _S622 = _S640.differential_0;
        _S502 = _S638;
    }
    else
    {
        _S622 = 0.0f;
        _S502 = _S489;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S641;
    (&_S641)->primal_0 = _S482.primal_0;
    (&_S641)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S642;
    (&_S642)->primal_0 = _S482.primal_0;
    (&_S642)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S641, &_S642, _S622);
    float _S643 = _S498 * (*_s_dOut_4)[int(6)];
    float _S644 = _S498 * (*_s_dOut_4)[int(5)];
    float _S645 = _S496 * (*_s_dOut_4)[int(4)];
    float _S646 = _S494 * (_S496 * (*_s_dOut_4)[int(6)] + _S644 + _S644 + (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S647;
    (&_S647)->primal_0 = _S497;
    (&_S647)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S647, _S646);
    DiffPair_float_0 _S648;
    (&_S648)->primal_0 = _S481.primal_0;
    (&_S648)->differential_0 = 0.0f;
    DiffPair_float_0 _S649;
    (&_S649)->primal_0 = 0.00009999999747379f;
    (&_S649)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S648, &_S649, _S647.differential_0);
    float _S650 = _S494 * (_S643 + _S645 + _S645 + (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = _S495;
    (&_S651)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S651, _S650);
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = _S480.primal_0;
    (&_S652)->differential_0 = 0.0f;
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = 0.00009999999747379f;
    (&_S653)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S652, &_S653, _S651.differential_0);
    float _S654 = _S490 * (*_s_dOut_4)[int(1)];
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = _S493;
    (&_S655)->differential_0 = 0.0f;
    DiffPair_float_0 _S656;
    (&_S656)->primal_0 = 0.0f;
    (&_S656)->differential_0 = 0.0f;
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = 1.0f;
    (&_S657)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S655, &_S656, &_S657, _S654);
    float _S658 = 0.3333333432674408f * _S655.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S659;
    (&_S659)->primal_0 = _S492;
    (&_S659)->differential_0 = _S489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S660;
    (&_S660)->primal_0 = _S492;
    (&_S660)->differential_0 = _S489;
    s_bwd_prop_dot_0(&_S659, &_S660, _S658);
    float _S661 = 0.3333333432674408f * (_S491 * (*_s_dOut_4)[int(0)]);
    float3  _S662 = make_float3 (_S661, _S661, _S661);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S663;
    (&_S663)->primal_0 = _S492;
    (&_S663)->differential_0 = _S489;
    s_bwd_prop_abs_0(&_S663, _S662);
    float3  _S664 = _S660.differential_0 + _S659.differential_0 + _S663.differential_0;
    float3  _S665 = - _S664;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S619;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S555;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S620;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S595.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S629;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S636;
    float3  _S666 = _S642.differential_0 + _S641.differential_0 + _S502;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S666;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S648.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S652.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S664;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S665;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S667, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S668, DiffPair_float_0 * _S669, DiffPair_float_0 * _S670, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S671, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S672, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S673, DiffPair_float_0 * _S674, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S675, DiffPair_float_0 * _S676, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S677, bool _S678, bool _S679, bool _S680, bool _S681, bool _S682, FixedArray<float, 10>  * _S683, FixedArray<float, 23>  * _S684)
{
    s_bwd_prop_per_pixel_losses_0(_S667, _S668, _S669, _S670, _S671, _S672, _S673, _S674, _S675, _S676, _S677, _S678, _S679, _S680, _S681, _S682, _S683, _S684);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  * weights_2, FixedArray<float, 23>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S685 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S685;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S685;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S685;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S685;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S685;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S685;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S685;
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
    float _S686 = 1.0f / ((*dpx_14).primal_0 * 52.30258560180664062f) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S686;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  * raw_losses_0, FixedArray<float, 10>  * weights_3, FixedArray<float, 10>  * _S687)
{
    FixedArray<float, 10>  losses_2;
    float _S688 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*raw_losses_0)[int(0)] / _S688;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S688)));
    bool _S689;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S689 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S689 = false;
    }
    float _S690;
    if(_S689)
    {
        _S690 = (*weights_3)[int(1)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S690 = 0.0f;
    }
    losses_2[int(2)] = _S690;
    losses_2[int(3)] = ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = ((*raw_losses_0)[int(9)] + (*raw_losses_0)[int(10)]) / (F32_max(((*raw_losses_0)[int(22)]), (1.0f)));
    losses_2[int(5)] = (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    float _S691 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(6)] = (*raw_losses_0)[int(12)] / _S691;
    losses_2[int(7)] = (*raw_losses_0)[int(13)] / _S691;
    losses_2[int(8)] = (*raw_losses_0)[int(14)] / _S691;
    losses_2[int(9)] = (*raw_losses_0)[int(15)] / _S691;
    *_S687 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S692)
{
    return (F32_sqrt((_S692)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S693, float _S694)
{
    _d_log10_0(_S693, _S694);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S695 = dpraw_losses_0->primal_0;
    float _S696 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S697 = _S696 * _S696;
    float _S698 = dpraw_losses_0->primal_0[int(1)] / _S696;
    bool _S699 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S700;
    if(_S699)
    {
        _S700 = (_S695[int(3)]) != 0.0f;
    }
    else
    {
        _S700 = false;
    }
    float _S701;
    float _S702;
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
    if(_S700)
    {
        float _S716 = _S695[int(2)] * _S695[int(3)];
        float _S717 = _S695[int(18)] * _S695[int(18)];
        float _S718 = _S695[int(6)] - _S716 / _S695[int(18)];
        float _S719 = _S695[int(2)] * _S695[int(2)];
        float _S720 = _S695[int(4)] - _S719 / _S695[int(18)];
        float _S721 = _S695[int(3)] * _S695[int(3)];
        float _S722 = _S695[int(5)] - _S721 / _S695[int(18)];
        float _S723 = _S720 * _S722 + 1.0f;
        float _S724 = s_primal_ctx_max_0(9.999999960041972e-13f, _S723);
        float _S725 = s_primal_ctx_sqrt_0(_S724);
        float _S726 = _S725 * _S725;
        float _S727 = 1.0f - _S718 / _S725;
        _S701 = (*weights_4)[int(1)];
        _S702 = _S727;
        _S703 = _S726;
        _S704 = _S718;
        _S705 = _S725;
        _S706 = _S724;
        _S707 = _S723;
        _S708 = _S720;
        _S709 = _S722;
        _S710 = _S717;
        _S711 = _S721;
        _S712 = _S695[int(3)];
        _S713 = _S719;
        _S714 = _S695[int(2)];
        _S715 = _S716;
    }
    else
    {
        _S701 = 0.0f;
        _S702 = 0.0f;
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
    }
    float _S728 = s_primal_ctx_max_0(_S695[int(19)], 1.0f);
    float _S729 = _S728 * _S728;
    float _S730 = s_primal_ctx_max_0(_S695[int(20)], 1.0f);
    float _S731 = _S730 * _S730;
    float _S732 = float((I32_max((int((_S695[int(19)]) > 0.5f) + int((_S695[int(20)]) > 0.5f)), (int(1)))));
    float _S733 = _S695[int(9)] + _S695[int(10)];
    float _S734 = s_primal_ctx_max_0(_S695[int(22)], 1.0f);
    float _S735 = _S734 * _S734;
    float _S736 = s_primal_ctx_max_0(_S695[int(21)], 1.0f);
    float _S737 = _S736 * _S736;
    float _S738 = s_primal_ctx_max_0(_S695[int(16)], 1.0f);
    float _S739 = _S738 * _S738;
    float _S740 = (*_s_dOut_5)[int(9)] / _S739;
    float _S741 = _S738 * _S740;
    float _S742 = (*_s_dOut_5)[int(8)] / _S739;
    float _S743 = _S738 * _S742;
    float _S744 = (*_s_dOut_5)[int(7)] / _S739;
    float _S745 = _S738 * _S744;
    float _S746 = (*_s_dOut_5)[int(6)] / _S739;
    float _S747 = _S738 * _S746;
    float _S748 = _S695[int(15)] * - _S740 + _S695[int(14)] * - _S742 + _S695[int(13)] * - _S744 + _S695[int(12)] * - _S746;
    DiffPair_float_0 _S749;
    (&_S749)->primal_0 = _S695[int(16)];
    (&_S749)->differential_0 = 0.0f;
    DiffPair_float_0 _S750;
    (&_S750)->primal_0 = 1.0f;
    (&_S750)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S749, &_S750, _S748);
    float _S751 = (*_s_dOut_5)[int(5)] / _S737;
    float _S752 = _S695[int(11)] * - _S751;
    float _S753 = _S736 * _S751;
    DiffPair_float_0 _S754;
    (&_S754)->primal_0 = _S695[int(21)];
    (&_S754)->differential_0 = 0.0f;
    DiffPair_float_0 _S755;
    (&_S755)->primal_0 = 1.0f;
    (&_S755)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S754, &_S755, _S752);
    float _S756 = (*_s_dOut_5)[int(4)] / _S735;
    float _S757 = _S733 * - _S756;
    float _S758 = _S734 * _S756;
    DiffPair_float_0 _S759;
    (&_S759)->primal_0 = _S695[int(22)];
    (&_S759)->differential_0 = 0.0f;
    DiffPair_float_0 _S760;
    (&_S760)->primal_0 = 1.0f;
    (&_S760)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S759, &_S760, _S757);
    float _S761 = (*_s_dOut_5)[int(3)] / _S732;
    float _S762 = _S761 / _S731;
    float _S763 = _S695[int(8)] * - _S762;
    float _S764 = _S730 * _S762;
    DiffPair_float_0 _S765;
    (&_S765)->primal_0 = _S695[int(20)];
    (&_S765)->differential_0 = 0.0f;
    DiffPair_float_0 _S766;
    (&_S766)->primal_0 = 1.0f;
    (&_S766)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S765, &_S766, _S763);
    float _S767 = _S761 / _S729;
    float _S768 = _S695[int(7)] * - _S767;
    float _S769 = _S728 * _S767;
    DiffPair_float_0 _S770;
    (&_S770)->primal_0 = _S695[int(19)];
    (&_S770)->differential_0 = 0.0f;
    DiffPair_float_0 _S771;
    (&_S771)->primal_0 = 1.0f;
    (&_S771)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S770, &_S771, _S768);
    FixedArray<float, 23>  _S772;
    _S772[int(0)] = 0.0f;
    _S772[int(1)] = 0.0f;
    _S772[int(2)] = 0.0f;
    _S772[int(3)] = 0.0f;
    _S772[int(4)] = 0.0f;
    _S772[int(5)] = 0.0f;
    _S772[int(6)] = 0.0f;
    _S772[int(7)] = 0.0f;
    _S772[int(8)] = 0.0f;
    _S772[int(9)] = 0.0f;
    _S772[int(10)] = 0.0f;
    _S772[int(11)] = 0.0f;
    _S772[int(12)] = 0.0f;
    _S772[int(13)] = 0.0f;
    _S772[int(14)] = 0.0f;
    _S772[int(15)] = 0.0f;
    _S772[int(16)] = 0.0f;
    _S772[int(17)] = 0.0f;
    _S772[int(18)] = 0.0f;
    _S772[int(19)] = 0.0f;
    _S772[int(20)] = 0.0f;
    _S772[int(21)] = 0.0f;
    _S772[int(22)] = 0.0f;
    _S772[int(15)] = _S741;
    _S772[int(14)] = _S743;
    _S772[int(13)] = _S745;
    _S772[int(16)] = _S749.differential_0;
    _S772[int(12)] = _S747;
    _S772[int(21)] = _S754.differential_0;
    _S772[int(11)] = _S753;
    _S772[int(22)] = _S759.differential_0;
    _S772[int(10)] = _S758;
    _S772[int(9)] = _S758;
    _S772[int(20)] = _S765.differential_0;
    _S772[int(8)] = _S764;
    _S772[int(19)] = _S770.differential_0;
    _S772[int(7)] = _S769;
    float _S773 = _S772[int(0)];
    float _S774 = _S772[int(1)];
    float _S775 = _S772[int(2)];
    float _S776 = _S772[int(3)];
    float _S777 = _S772[int(4)];
    float _S778 = _S772[int(5)];
    float _S779 = _S772[int(6)];
    float _S780 = _S772[int(7)];
    float _S781 = _S772[int(8)];
    float _S782 = _S772[int(9)];
    float _S783 = _S772[int(10)];
    float _S784 = _S772[int(11)];
    float _S785 = _S772[int(12)];
    float _S786 = _S772[int(13)];
    float _S787 = _S772[int(14)];
    float _S788 = _S772[int(15)];
    float _S789 = _S772[int(16)];
    float _S790 = _S772[int(17)];
    float _S791 = _S772[int(18)];
    float _S792 = _S772[int(19)];
    float _S793 = _S772[int(20)];
    float _S794 = _S772[int(21)];
    float _S795 = _S772[int(22)];
    FixedArray<float, 23>  _S796;
    if(_S700)
    {
        float _S797 = _S701 * (*_s_dOut_5)[int(2)];
        DiffPair_float_0 _S798;
        (&_S798)->primal_0 = _S702;
        (&_S798)->differential_0 = 0.0f;
        DiffPair_float_0 _S799;
        (&_S799)->primal_0 = 0.0f;
        (&_S799)->differential_0 = 0.0f;
        DiffPair_float_0 _S800;
        (&_S800)->primal_0 = 2.0f;
        (&_S800)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S798, &_S799, &_S800, _S797);
        float _S801 = - _S798.differential_0 / _S703;
        float _S802 = _S704 * - _S801;
        float _S803 = _S705 * _S801;
        DiffPair_float_0 _S804;
        (&_S804)->primal_0 = _S706;
        (&_S804)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S804, _S802);
        DiffPair_float_0 _S805;
        (&_S805)->primal_0 = 9.999999960041972e-13f;
        (&_S805)->differential_0 = 0.0f;
        DiffPair_float_0 _S806;
        (&_S806)->primal_0 = _S707;
        (&_S806)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S805, &_S806, _S804.differential_0);
        float _S807 = _S708 * _S806.differential_0;
        float _S808 = _S709 * _S806.differential_0;
        float _S809 = - _S807 / _S710;
        float _S810 = _S712 * (_S695[int(18)] * _S809);
        float _S811 = - _S808 / _S710;
        float _S812 = _S714 * (_S695[int(18)] * _S811);
        float _S813 = - _S803 / _S710;
        float _S814 = _S695[int(18)] * _S813;
        float _S815 = _S810 + _S810 + _S714 * _S814;
        float _S816 = _S812 + _S812 + _S712 * _S814;
        float _S817 = _S711 * - _S809 + _S713 * - _S811 + _S715 * - _S813;
        FixedArray<float, 23>  _S818;
        _S818[int(0)] = 0.0f;
        _S818[int(1)] = 0.0f;
        _S818[int(2)] = 0.0f;
        _S818[int(3)] = 0.0f;
        _S818[int(4)] = 0.0f;
        _S818[int(5)] = 0.0f;
        _S818[int(6)] = 0.0f;
        _S818[int(7)] = 0.0f;
        _S818[int(8)] = 0.0f;
        _S818[int(9)] = 0.0f;
        _S818[int(10)] = 0.0f;
        _S818[int(11)] = 0.0f;
        _S818[int(12)] = 0.0f;
        _S818[int(13)] = 0.0f;
        _S818[int(14)] = 0.0f;
        _S818[int(15)] = 0.0f;
        _S818[int(16)] = 0.0f;
        _S818[int(17)] = 0.0f;
        _S818[int(18)] = 0.0f;
        _S818[int(19)] = 0.0f;
        _S818[int(20)] = 0.0f;
        _S818[int(21)] = 0.0f;
        _S818[int(22)] = 0.0f;
        _S818[int(5)] = _S807;
        _S818[int(4)] = _S808;
        _S818[int(3)] = _S815;
        _S818[int(2)] = _S816;
        _S818[int(6)] = _S803;
        float _S819 = _S774 + _S818[int(1)];
        float _S820 = _S775 + _S818[int(2)];
        float _S821 = _S776 + _S818[int(3)];
        float _S822 = _S777 + _S818[int(4)];
        float _S823 = _S778 + _S818[int(5)];
        float _S824 = _S779 + _S818[int(6)];
        float _S825 = _S780 + _S818[int(7)];
        float _S826 = _S781 + _S818[int(8)];
        float _S827 = _S782 + _S818[int(9)];
        float _S828 = _S783 + _S818[int(10)];
        float _S829 = _S784 + _S818[int(11)];
        float _S830 = _S785 + _S818[int(12)];
        float _S831 = _S786 + _S818[int(13)];
        float _S832 = _S787 + _S818[int(14)];
        float _S833 = _S788 + _S818[int(15)];
        float _S834 = _S789 + _S818[int(16)];
        float _S835 = _S790 + _S818[int(17)];
        float _S836 = _S791 + _S818[int(18)];
        float _S837 = _S792 + _S818[int(19)];
        float _S838 = _S793 + _S818[int(20)];
        float _S839 = _S794 + _S818[int(21)];
        float _S840 = _S795 + _S818[int(22)];
        _S796[int(0)] = _S773 + _S818[int(0)];
        _S796[int(1)] = _S819;
        _S796[int(2)] = _S820;
        _S796[int(3)] = _S821;
        _S796[int(4)] = _S822;
        _S796[int(5)] = _S823;
        _S796[int(6)] = _S824;
        _S796[int(7)] = _S825;
        _S796[int(8)] = _S826;
        _S796[int(9)] = _S827;
        _S796[int(10)] = _S828;
        _S796[int(11)] = _S829;
        _S796[int(12)] = _S830;
        _S796[int(13)] = _S831;
        _S796[int(14)] = _S832;
        _S796[int(15)] = _S833;
        _S796[int(16)] = _S834;
        _S796[int(17)] = _S835;
        _S796[int(18)] = _S836;
        _S796[int(19)] = _S837;
        _S796[int(20)] = _S838;
        _S796[int(21)] = _S839;
        _S796[int(22)] = _S840;
        _S701 = _S817;
    }
    else
    {
        _S796[int(0)] = _S773;
        _S796[int(1)] = _S774;
        _S796[int(2)] = _S775;
        _S796[int(3)] = _S776;
        _S796[int(4)] = _S777;
        _S796[int(5)] = _S778;
        _S796[int(6)] = _S779;
        _S796[int(7)] = _S780;
        _S796[int(8)] = _S781;
        _S796[int(9)] = _S782;
        _S796[int(10)] = _S783;
        _S796[int(11)] = _S784;
        _S796[int(12)] = _S785;
        _S796[int(13)] = _S786;
        _S796[int(14)] = _S787;
        _S796[int(15)] = _S788;
        _S796[int(16)] = _S789;
        _S796[int(17)] = _S790;
        _S796[int(18)] = _S791;
        _S796[int(19)] = _S792;
        _S796[int(20)] = _S793;
        _S796[int(21)] = _S794;
        _S796[int(22)] = _S795;
        _S701 = 0.0f;
    }
    if(_S699)
    {
        FixedArray<float, 23>  _S841;
        _S841[int(0)] = 0.0f;
        _S841[int(1)] = 0.0f;
        _S841[int(2)] = 0.0f;
        _S841[int(3)] = 0.0f;
        _S841[int(4)] = 0.0f;
        _S841[int(5)] = 0.0f;
        _S841[int(6)] = 0.0f;
        _S841[int(7)] = 0.0f;
        _S841[int(8)] = 0.0f;
        _S841[int(9)] = 0.0f;
        _S841[int(10)] = 0.0f;
        _S841[int(11)] = 0.0f;
        _S841[int(12)] = 0.0f;
        _S841[int(13)] = 0.0f;
        _S841[int(14)] = 0.0f;
        _S841[int(15)] = 0.0f;
        _S841[int(16)] = 0.0f;
        _S841[int(17)] = 0.0f;
        _S841[int(18)] = 0.0f;
        _S841[int(19)] = 0.0f;
        _S841[int(20)] = 0.0f;
        _S841[int(21)] = 0.0f;
        _S841[int(22)] = 0.0f;
        _S841[int(3)] = 0.0f;
        float _S842 = _S796[int(1)] + _S841[int(1)];
        float _S843 = _S796[int(2)] + _S841[int(2)];
        float _S844 = _S796[int(3)] + _S841[int(3)];
        float _S845 = _S796[int(4)] + _S841[int(4)];
        float _S846 = _S796[int(5)] + _S841[int(5)];
        float _S847 = _S796[int(6)] + _S841[int(6)];
        float _S848 = _S796[int(7)] + _S841[int(7)];
        float _S849 = _S796[int(8)] + _S841[int(8)];
        float _S850 = _S796[int(9)] + _S841[int(9)];
        float _S851 = _S796[int(10)] + _S841[int(10)];
        float _S852 = _S796[int(11)] + _S841[int(11)];
        float _S853 = _S796[int(12)] + _S841[int(12)];
        float _S854 = _S796[int(13)] + _S841[int(13)];
        float _S855 = _S796[int(14)] + _S841[int(14)];
        float _S856 = _S796[int(15)] + _S841[int(15)];
        float _S857 = _S796[int(16)] + _S841[int(16)];
        float _S858 = _S796[int(17)] + _S841[int(17)];
        float _S859 = _S796[int(18)] + _S841[int(18)];
        float _S860 = _S796[int(19)] + _S841[int(19)];
        float _S861 = _S796[int(20)] + _S841[int(20)];
        float _S862 = _S796[int(21)] + _S841[int(21)];
        float _S863 = _S796[int(22)] + _S841[int(22)];
        _S796[int(0)] = _S796[int(0)] + _S841[int(0)];
        _S796[int(1)] = _S842;
        _S796[int(2)] = _S843;
        _S796[int(3)] = _S844;
        _S796[int(4)] = _S845;
        _S796[int(5)] = _S846;
        _S796[int(6)] = _S847;
        _S796[int(7)] = _S848;
        _S796[int(8)] = _S849;
        _S796[int(9)] = _S850;
        _S796[int(10)] = _S851;
        _S796[int(11)] = _S852;
        _S796[int(12)] = _S853;
        _S796[int(13)] = _S854;
        _S796[int(14)] = _S855;
        _S796[int(15)] = _S856;
        _S796[int(16)] = _S857;
        _S796[int(17)] = _S858;
        _S796[int(18)] = _S859;
        _S796[int(19)] = _S860;
        _S796[int(20)] = _S861;
        _S796[int(21)] = _S862;
        _S796[int(22)] = _S863;
    }
    float _S864 = -10.0f * (*_s_dOut_5)[int(1)];
    DiffPair_float_0 _S865;
    (&_S865)->primal_0 = _S698;
    (&_S865)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S865, _S864);
    float _S866 = _S865.differential_0 / _S697;
    float _S867 = _S696 * _S866;
    float _S868 = (*_s_dOut_5)[int(0)] / _S697;
    float _S869 = _S696 * _S868;
    float _S870 = _S695[int(1)] * - _S866 + _S695[int(0)] * - _S868;
    DiffPair_float_0 _S871;
    (&_S871)->primal_0 = _S695[int(17)];
    (&_S871)->differential_0 = 0.0f;
    DiffPair_float_0 _S872;
    (&_S872)->primal_0 = 1.0f;
    (&_S872)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S871, &_S872, _S870);
    FixedArray<float, 23>  _S873;
    _S873[int(0)] = 0.0f;
    _S873[int(1)] = 0.0f;
    _S873[int(2)] = 0.0f;
    _S873[int(3)] = 0.0f;
    _S873[int(4)] = 0.0f;
    _S873[int(5)] = 0.0f;
    _S873[int(6)] = 0.0f;
    _S873[int(7)] = 0.0f;
    _S873[int(8)] = 0.0f;
    _S873[int(9)] = 0.0f;
    _S873[int(10)] = 0.0f;
    _S873[int(11)] = 0.0f;
    _S873[int(12)] = 0.0f;
    _S873[int(13)] = 0.0f;
    _S873[int(14)] = 0.0f;
    _S873[int(15)] = 0.0f;
    _S873[int(16)] = 0.0f;
    _S873[int(17)] = 0.0f;
    _S873[int(18)] = 0.0f;
    _S873[int(19)] = 0.0f;
    _S873[int(20)] = 0.0f;
    _S873[int(21)] = 0.0f;
    _S873[int(22)] = 0.0f;
    _S873[int(18)] = _S701;
    _S873[int(1)] = _S867;
    _S873[int(17)] = _S871.differential_0;
    _S873[int(0)] = _S869;
    FixedArray<float, 23>  _S874 = {
        _S796[int(0)] + _S873[int(0)], _S796[int(1)] + _S873[int(1)], _S796[int(2)] + _S873[int(2)], _S796[int(3)] + _S873[int(3)], _S796[int(4)] + _S873[int(4)], _S796[int(5)] + _S873[int(5)], _S796[int(6)] + _S873[int(6)], _S796[int(7)] + _S873[int(7)], _S796[int(8)] + _S873[int(8)], _S796[int(9)] + _S873[int(9)], _S796[int(10)] + _S873[int(10)], _S796[int(11)] + _S873[int(11)], _S796[int(12)] + _S873[int(12)], _S796[int(13)] + _S873[int(13)], _S796[int(14)] + _S873[int(14)], _S796[int(15)] + _S873[int(15)], _S796[int(16)] + _S873[int(16)], _S796[int(17)] + _S873[int(17)], _S796[int(18)] + _S873[int(18)], _S796[int(19)] + _S873[int(19)], _S796[int(20)] + _S873[int(20)], _S796[int(21)] + _S873[int(21)], _S796[int(22)] + _S873[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S874;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S875, FixedArray<float, 10>  * _S876, FixedArray<float, 10>  * _S877)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S875, _S876, _S877);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  * raw_losses_1, FixedArray<float, 10>  * weights_5, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 23>  * _S878)
{
    FixedArray<float, 23>  _S879 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S879;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_5, v_losses_1);
    *_S878 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_25, float3  y_8)
{
    float3  result_19;
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
        *_slang_vector_get_element_ptr(&result_19, i_12) = (F32_min((_slang_vector_get_element(x_25, i_12)), (_slang_vector_get_element(y_8, i_12))));
        i_12 = i_12 + int(1);
    }
    return result_19;
}

inline __device__ float3  max_0(float3  x_26, float3  y_9)
{
    float3  result_20;
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
        *_slang_vector_get_element_ptr(&result_20, i_13) = (F32_max((_slang_vector_get_element(x_26, i_13)), (_slang_vector_get_element(y_9, i_13))));
        i_13 = i_13 + int(1);
    }
    return result_20;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S880, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S881, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S882, float3  _S883)
{
    _d_clamp_vector_0(_S880, _S881, _S882, _S883);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S884 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S885 = make_float3 (_S884);
    float3  _S886 = make_float3 (0.0f);
    float3  _S887 = make_float3 (1.0f);
    float3  _S888 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S889;
    (&_S889)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S884) * (*dpbackground_0).primal_0;
    (&_S889)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S890;
    (&_S890)->primal_0 = _S886;
    (&_S890)->differential_0 = _S888;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S891;
    (&_S891)->primal_0 = _S887;
    (&_S891)->differential_0 = _S888;
    s_bwd_prop_clamp_1(&_S889, &_S890, &_S891, _s_dOut_6);
    float3  _S892 = _S885 * _S889.differential_0;
    float3  _S893 = (*dpbackground_0).primal_0 * _S889.differential_0;
    float _S894 = - (_S893.x + _S893.y + _S893.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S892;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S894;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S889.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S895, DiffPair_float_0 * _S896, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S897, float3  _S898)
{
    s_bwd_prop_blend_background_0(_S895, _S896, _S897, _S898);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S899 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S899;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S899;
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
        DiffPair_float_0 _S900 = *dpx_16;
        float _S901 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_20;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S901;
        float _S902 = val_0 * (F32_log((_S900.primal_0))) * dOut_20;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S902;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S903 = rgb_2;
    float _S904;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S904 = _S903.x * 12.92000007629394531f;
    }
    else
    {
        _S904 = 1.0549999475479126f * (F32_pow((_S903.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S903)->x) = _S904;
    if((_S903.y) < 0.00313080009073019f)
    {
        _S904 = _S903.y * 12.92000007629394531f;
    }
    else
    {
        _S904 = 1.0549999475479126f * (F32_pow((_S903.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S903)->y) = _S904;
    if((_S903.z) < 0.00313080009073019f)
    {
        _S904 = _S903.z * 12.92000007629394531f;
    }
    else
    {
        _S904 = 1.0549999475479126f * (F32_pow((_S903.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S903)->z) = _S904;
    return _S903;
}

inline __device__ float s_primal_ctx_pow_0(float _S905, float _S906)
{
    return (F32_pow((_S905), (_S906)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S907, DiffPair_float_0 * _S908, float _S909)
{
    _d_pow_0(_S907, _S908, _S909);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_7)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S910 = *dprgb_1;
    float _S911 = (*dprgb_1).primal_0.x;
    bool _S912 = _S911 < 0.00313080009073019f;
    float _S913;
    if(_S912)
    {
        _S913 = _S911 * 12.92000007629394531f;
    }
    else
    {
        _S913 = 1.0549999475479126f * s_primal_ctx_pow_0(_S911, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S914 = _S910.primal_0;
    *&((&_S914)->x) = _S913;
    float _S915 = _S914.y;
    bool _S916 = _S915 < 0.00313080009073019f;
    if(_S916)
    {
        _S913 = _S915 * 12.92000007629394531f;
    }
    else
    {
        _S913 = 1.0549999475479126f * s_primal_ctx_pow_0(_S915, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S914)->y) = _S913;
    float _S917 = _S914.z;
    bool _S918 = _S917 < 0.00313080009073019f;
    _S914 = _s_dOut_7;
    *&((&_S914)->z) = 0.0f;
    if(_S918)
    {
        _S913 = 12.92000007629394531f * _s_dOut_7.z;
    }
    else
    {
        float _S919 = 1.0549999475479126f * _s_dOut_7.z;
        DiffPair_float_0 _S920;
        (&_S920)->primal_0 = _S917;
        (&_S920)->differential_0 = 0.0f;
        DiffPair_float_0 _S921;
        (&_S921)->primal_0 = 0.4166666567325592f;
        (&_S921)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S920, &_S921, _S919);
        _S913 = _S920.differential_0;
    }
    float3  _S922 = _S914 + make_float3 (0.0f, 0.0f, _S913);
    _S914 = _S922;
    *&((&_S914)->y) = 0.0f;
    if(_S916)
    {
        _S913 = 12.92000007629394531f * _S922.y;
    }
    else
    {
        float _S923 = 1.0549999475479126f * _S922.y;
        DiffPair_float_0 _S924;
        (&_S924)->primal_0 = _S915;
        (&_S924)->differential_0 = 0.0f;
        DiffPair_float_0 _S925;
        (&_S925)->primal_0 = 0.4166666567325592f;
        (&_S925)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S924, &_S925, _S923);
        _S913 = _S924.differential_0;
    }
    float3  _S926 = _S914 + make_float3 (0.0f, _S913, 0.0f);
    _S914 = _S926;
    *&((&_S914)->x) = 0.0f;
    if(_S912)
    {
        _S913 = 12.92000007629394531f * _S926.x;
    }
    else
    {
        float _S927 = 1.0549999475479126f * _S926.x;
        DiffPair_float_0 _S928;
        (&_S928)->primal_0 = _S911;
        (&_S928)->differential_0 = 0.0f;
        DiffPair_float_0 _S929;
        (&_S929)->primal_0 = 0.4166666567325592f;
        (&_S929)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S928, &_S929, _S927);
        _S913 = _S928.differential_0;
    }
    float3  _S930 = _S914 + make_float3 (_S913, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S930;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S931, float3  _S932)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S931, _S932);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S933 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S933;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_7, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S934 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = _S934;
    bool _S935 = undistort_point_0(_S934, dist_coeffs_7, int(12), &uv_7);
    if(!_S935)
    {
        int3  _S936 = make_int3 (int(0));
        float3  _S937 = make_float3 ((float)_S936.x, (float)_S936.y, (float)_S936.z);
        return _S937;
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S938, float3  _S939)
{
    return cross_0(_S938, _S939);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_8)
{
    float _S940 = (*dpx_17).primal_0.x;
    float _S941 = (*dpx_17).primal_0.y;
    float _S942 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S943;
    (&_S943)->primal_0 = _S940 * _S940 + _S941 * _S941 + _S942 * _S942;
    (&_S943)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S943, _s_dOut_8);
    float _S944 = (*dpx_17).primal_0.z * _S943.differential_0;
    float _S945 = _S944 + _S944;
    float _S946 = (*dpx_17).primal_0.y * _S943.differential_0;
    float _S947 = _S946 + _S946;
    float _S948 = (*dpx_17).primal_0.x * _S943.differential_0;
    float _S949 = _S948 + _S948;
    float3  _S950 = make_float3 (0.0f);
    *&((&_S950)->z) = _S945;
    *&((&_S950)->y) = _S947;
    *&((&_S950)->x) = _S949;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S950;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S951, float _S952)
{
    s_bwd_prop_length_impl_1(_S951, _S952);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S953, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S954, float3  _S955)
{
    _d_cross_0(_S953, _S954, _S955);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S956 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S957 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S958 = s_primal_ctx_cross_0(dx_0, _S957);
    bool _S959 = (s_primal_ctx_dot_0(_S958, _S958)) != 0.0f;
    float3  _S960;
    float3  _S961;
    if(_S959)
    {
        float _S962 = length_2(_S958);
        float3  _S963 = make_float3 (_S962);
        _S960 = make_float3 (_S962 * _S962);
        _S961 = _S963;
    }
    else
    {
        _S960 = _S956;
        _S961 = _S956;
    }
    if(_S959)
    {
        float3  _S964 = _s_dOut_9 / _S960;
        float3  _S965 = _S958 * - _S964;
        float3  _S966 = _S961 * _S964;
        float _S967 = _S965.x + _S965.y + _S965.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S968;
        (&_S968)->primal_0 = _S958;
        (&_S968)->differential_0 = _S956;
        s_bwd_length_impl_1(&_S968, _S967);
        _S960 = _S966 + _S968.differential_0;
    }
    else
    {
        _S960 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S969;
    (&_S969)->primal_0 = _S958;
    (&_S969)->differential_0 = _S956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S970;
    (&_S970)->primal_0 = _S958;
    (&_S970)->differential_0 = _S956;
    s_bwd_prop_dot_0(&_S969, &_S970, 0.0f);
    float3  _S971 = _S970.differential_0 + _S969.differential_0 + _S960;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S972;
    (&_S972)->primal_0 = dx_0;
    (&_S972)->differential_0 = _S956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S973;
    (&_S973)->primal_0 = _S957;
    (&_S973)->differential_0 = _S956;
    s_bwd_prop_cross_0(&_S972, &_S973, _S971);
    float3  s_diff_dy_T_0 = - _S973.differential_0;
    float3  _S974 = - s_diff_dy_T_0;
    float3  _S975 = - _S972.differential_0;
    FixedArray<float3 , 4>  _S976;
    _S976[int(0)] = _S956;
    _S976[int(1)] = _S956;
    _S976[int(2)] = _S956;
    _S976[int(3)] = _S956;
    _S976[int(2)] = _S974;
    _S976[int(3)] = s_diff_dy_T_0;
    _S976[int(0)] = _S975;
    _S976[int(1)] = _S972.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S976;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S977, float3  _S978)
{
    s_bwd_prop_points_to_normal_0(_S977, _S978);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  * points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S979 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = *points_1;
    (&dp_points_0)->differential_0 = _S979;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S980 = float2 {intrins_1.z, intrins_1.w};
    float2  _S981 = float2 {intrins_1.x, intrins_1.y};
    float2  _S982 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S980) / _S981;
    float2  uv_8 = _S982;
    bool _S983 = undistort_point_0(_S982, dist_coeffs_8, int(12), &uv_8);
    if(!_S983)
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
    float2  _S984 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S980) / _S981;
    float2  uv_9 = _S984;
    bool _S985 = undistort_point_0(_S984, dist_coeffs_8, int(12), &uv_9);
    if(!_S985)
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
    float2  _S986 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S980) / _S981;
    float2  uv_10 = _S986;
    bool _S987 = undistort_point_0(_S986, dist_coeffs_8, int(12), &uv_10);
    if(!_S987)
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
    float2  _S988 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S980) / _S981;
    float2  uv_11 = _S988;
    bool _S989 = undistort_point_0(_S988, dist_coeffs_8, int(12), &uv_11);
    if(!_S989)
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
    float2  _S990;
    bool _S991;
    float2  _S992;
    bool _S993;
    float2  _S994;
    bool _S995;
    float2  _S996;
    bool _S997;
};

inline __device__ float s_primal_ctx_sin_0(float _S998)
{
    return (F32_sin((_S998)));
}

inline __device__ float s_primal_ctx_cos_0(float _S999)
{
    return (F32_cos((_S999)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1000 = make_float2 (0.0f);
    _s_diff_ctx_0->_S990 = _S1000;
    _s_diff_ctx_0->_S991 = false;
    _s_diff_ctx_0->_S992 = _S1000;
    _s_diff_ctx_0->_S993 = false;
    _s_diff_ctx_0->_S994 = _S1000;
    _s_diff_ctx_0->_S995 = false;
    _s_diff_ctx_0->_S996 = _S1000;
    _s_diff_ctx_0->_S997 = false;
    _s_diff_ctx_0->_S992 = _S1000;
    _s_diff_ctx_0->_S993 = false;
    _s_diff_ctx_0->_S994 = _S1000;
    _s_diff_ctx_0->_S995 = false;
    _s_diff_ctx_0->_S996 = _S1000;
    _s_diff_ctx_0->_S997 = false;
    float3  _S1001 = make_float3 (0.0f);
    float2  _S1002 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1003 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1004 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1002) / _S1003;
    float2  _S1005 = _S1004;
    bool _S1006 = undistort_point_0(_S1004, dist_coeffs_9, int(12), &_S1005);
    _s_diff_ctx_0->_S990 = _S1005;
    _s_diff_ctx_0->_S991 = _S1006;
    float2  uv_12 = _S1005;
    bool _S1007 = !_S1006;
    float3  normal_4;
    if(_S1007)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1008 = !_S1007;
    int _S1009;
    FixedArray<float3 , 4>  points_3;
    if(_S1008)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1010 = length_1(uv_12);
            float3  raydir_19 = make_float3 ((uv_12 / make_float2 (s_primal_ctx_max_0(_S1010, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1010))).x, (uv_12 / make_float2 (s_primal_ctx_max_0(_S1010, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1010))).y, s_primal_ctx_cos_0(_S1010));
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
        float3  _S1011 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1012 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1002) / _S1003;
        float2  _S1013 = _S1012;
        bool _S1014 = undistort_point_0(_S1012, dist_coeffs_9, int(12), &_S1013);
        _s_diff_ctx_0->_S992 = _S1013;
        _s_diff_ctx_0->_S993 = _S1014;
        float2  uv_13 = _S1013;
        bool _S1015 = !_S1014;
        if(_S1015)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1016 = !_S1015;
        if(_S1016)
        {
            if(is_fisheye_6)
            {
                float _S1017 = length_1(uv_13);
                float3  raydir_21 = make_float3 ((uv_13 / make_float2 (s_primal_ctx_max_0(_S1017, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1017))).x, (uv_13 / make_float2 (s_primal_ctx_max_0(_S1017, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1017))).y, s_primal_ctx_cos_0(_S1017));
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
            float3  _S1018 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1009 = int(2);
            points_3[int(0)] = _S1011;
            points_3[int(1)] = _S1018;
            points_3[int(2)] = _S1001;
            points_3[int(3)] = _S1001;
        }
        else
        {
            _S1009 = int(0);
            points_3[int(0)] = _S1011;
            points_3[int(1)] = _S1001;
            points_3[int(2)] = _S1001;
            points_3[int(3)] = _S1001;
        }
        bool _runFlag_0;
        if(_S1009 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1008;
            _S1009 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1019 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1002) / _S1003;
            float2  _S1020 = _S1019;
            bool _S1021 = undistort_point_0(_S1019, dist_coeffs_9, int(12), &_S1020);
            _s_diff_ctx_0->_S994 = _S1020;
            _s_diff_ctx_0->_S995 = _S1021;
            float2  uv_14 = _S1020;
            if(!_S1021)
            {
                float3  _S1022 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1009 = int(0);
                normal_4 = _S1022;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1023 = length_1(uv_14);
                    float3  raydir_23 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S1023, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1023))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S1023, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1023))).y, s_primal_ctx_cos_0(_S1023));
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
                float2  _S1024 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1002) / _S1003;
                float2  _S1025 = _S1024;
                bool _S1026 = undistort_point_0(_S1024, dist_coeffs_9, int(12), &_S1025);
                _s_diff_ctx_0->_S996 = _S1025;
                _s_diff_ctx_0->_S997 = _S1026;
                float2  uv_15 = _S1025;
                bool _S1027 = !_S1026;
                if(_S1027)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1028 = !_S1027;
                int _S1029;
                if(_S1028)
                {
                    if(is_fisheye_6)
                    {
                        float _S1030 = length_1(uv_15);
                        float3  raydir_25 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S1030, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1030))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S1030, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1030))).y, s_primal_ctx_cos_0(_S1030));
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
                    _S1029 = int(2);
                }
                else
                {
                    _S1029 = int(0);
                }
                if(_S1029 != int(2))
                {
                    _runFlag_0 = false;
                    _S1009 = _S1029;
                }
                if(_runFlag_0)
                {
                    _S1009 = int(1);
                }
            }
        }
    }
    else
    {
        _S1009 = int(0);
        points_3[int(0)] = _S1001;
        points_3[int(1)] = _S1001;
        points_3[int(2)] = _S1001;
        points_3[int(3)] = _S1001;
    }
    if(!(_S1009 != int(1)))
    {
        float3  _S1031 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1031, _S1031)) != 0.0f)
        {
            normal_4 = _S1031 / make_float3 (length_2(_S1031));
        }
        else
        {
            normal_4 = _S1031;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1032 = *dpdepths_1;
    float3  _S1033 = make_float3 (0.0f);
    bool _S1034 = !!_s_diff_ctx_1->_S991;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1035;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1036;
    if(_S1034)
    {
        if(is_fisheye_7)
        {
            float _S1037 = length_1(_s_diff_ctx_1->_S990);
            float3  raydir_31 = make_float3 ((_s_diff_ctx_1->_S990 / make_float2 (s_primal_ctx_max_0(_S1037, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1037))).x, (_s_diff_ctx_1->_S990 / make_float2 (s_primal_ctx_max_0(_S1037, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1037))).y, s_primal_ctx_cos_0(_S1037));
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
            float3  raydir_32 = make_float3 (_s_diff_ctx_1->_S990.x, _s_diff_ctx_1->_S990.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1038 = make_float3 (_S1032.primal_0.x) * raydir_27;
        bool _S1039 = !!_s_diff_ctx_1->_S993;
        if(_S1039)
        {
            if(is_fisheye_7)
            {
                float _S1040 = length_1(_s_diff_ctx_1->_S992);
                float3  raydir_33 = make_float3 ((_s_diff_ctx_1->_S992 / make_float2 (s_primal_ctx_max_0(_S1040, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1040))).x, (_s_diff_ctx_1->_S992 / make_float2 (s_primal_ctx_max_0(_S1040, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1040))).y, s_primal_ctx_cos_0(_S1040));
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
                float3  raydir_34 = make_float3 (_s_diff_ctx_1->_S992.x, _s_diff_ctx_1->_S992.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1041 = make_float3 (_S1032.primal_0.y) * raydir_28;
            _S1035 = int(2);
            points_4[int(0)] = _S1038;
            points_4[int(1)] = _S1041;
            points_4[int(2)] = _S1033;
            points_4[int(3)] = _S1033;
        }
        else
        {
            _S1035 = int(0);
            points_4[int(0)] = _S1038;
            points_4[int(1)] = _S1033;
            points_4[int(2)] = _S1033;
            points_4[int(3)] = _S1033;
            raydir_28 = _S1033;
        }
        if(_S1035 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1034;
            _S1035 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_1->_S995)
            {
                _runFlag_2 = false;
                _S1035 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1042 = length_1(_s_diff_ctx_1->_S994);
                    float3  raydir_35 = make_float3 ((_s_diff_ctx_1->_S994 / make_float2 (s_primal_ctx_max_0(_S1042, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1042))).x, (_s_diff_ctx_1->_S994 / make_float2 (s_primal_ctx_max_0(_S1042, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1042))).y, s_primal_ctx_cos_0(_S1042));
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
                    float3  raydir_36 = make_float3 (_s_diff_ctx_1->_S994.x, _s_diff_ctx_1->_S994.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1032.primal_0.z) * raydir_29;
                bool _S1043 = !!_s_diff_ctx_1->_S997;
                int _S1044;
                if(_S1043)
                {
                    if(is_fisheye_7)
                    {
                        float _S1045 = length_1(_s_diff_ctx_1->_S996);
                        float3  raydir_37 = make_float3 ((_s_diff_ctx_1->_S996 / make_float2 (s_primal_ctx_max_0(_S1045, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1045))).x, (_s_diff_ctx_1->_S996 / make_float2 (s_primal_ctx_max_0(_S1045, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1045))).y, s_primal_ctx_cos_0(_S1045));
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
                        float3  raydir_38 = make_float3 (_s_diff_ctx_1->_S996.x, _s_diff_ctx_1->_S996.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1032.primal_0.w) * raydir_30;
                    _S1044 = int(2);
                }
                else
                {
                    _S1044 = int(0);
                    raydir_30 = _S1033;
                }
                if(_S1044 != int(2))
                {
                    _runFlag_3 = false;
                    _S1035 = _S1044;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1035 = int(1);
                }
                float3  _S1046 = raydir_29;
                _runFlag_3 = _S1043;
                raydir_29 = raydir_30;
                raydir_30 = _S1046;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1033;
                raydir_30 = _S1033;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1033;
            raydir_30 = _S1033;
        }
        float3  _S1047 = raydir_27;
        float3  _S1048 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1036 = _S1039;
        raydir_29 = _S1048;
        raydir_30 = _S1047;
    }
    else
    {
        _S1035 = int(0);
        points_4[int(0)] = _S1033;
        points_4[int(1)] = _S1033;
        points_4[int(2)] = _S1033;
        points_4[int(3)] = _S1033;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1033;
        raydir_28 = _S1033;
        _S1036 = false;
        raydir_29 = _S1033;
        raydir_30 = _S1033;
    }
    bool _S1049 = !(_S1035 != int(1));
    float3  _S1050;
    float3  _S1051;
    float3  _S1052;
    float3  _S1053;
    float3  _S1054;
    bool _S1055;
    if(_S1049)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1056 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1057 = s_primal_ctx_cross_0(dx_1, _S1056);
        bool _S1058 = (s_primal_ctx_dot_0(_S1057, _S1057)) != 0.0f;
        if(_S1058)
        {
            float _S1059 = length_2(_S1057);
            float3  _S1060 = make_float3 (_S1059);
            _S1050 = make_float3 (_S1059 * _S1059);
            _S1051 = _S1060;
        }
        else
        {
            _S1050 = _S1033;
            _S1051 = _S1033;
        }
        float3  _S1061 = _S1051;
        _S1055 = _S1058;
        _S1051 = _S1057;
        _S1052 = _S1061;
        _S1053 = dx_1;
        _S1054 = _S1056;
    }
    else
    {
        _S1055 = false;
        _S1050 = _S1033;
        _S1051 = _S1033;
        _S1052 = _S1033;
        _S1053 = _S1033;
        _S1054 = _S1033;
    }
    float4  _S1062 = make_float4 (0.0f);
    if(_S1049)
    {
        if(_S1055)
        {
            float3  _S1063 = _s_dOut_10 / _S1050;
            float3  _S1064 = _S1051 * - _S1063;
            float3  _S1065 = _S1052 * _S1063;
            float _S1066 = _S1064.x + _S1064.y + _S1064.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1067;
            (&_S1067)->primal_0 = _S1051;
            (&_S1067)->differential_0 = _S1033;
            s_bwd_length_impl_1(&_S1067, _S1066);
            _S1050 = _S1065 + _S1067.differential_0;
        }
        else
        {
            _S1050 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1068;
        (&_S1068)->primal_0 = _S1051;
        (&_S1068)->differential_0 = _S1033;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1069;
        (&_S1069)->primal_0 = _S1051;
        (&_S1069)->differential_0 = _S1033;
        s_bwd_prop_dot_0(&_S1068, &_S1069, 0.0f);
        float3  _S1070 = _S1069.differential_0 + _S1068.differential_0 + _S1050;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1071;
        (&_S1071)->primal_0 = _S1053;
        (&_S1071)->differential_0 = _S1033;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1072;
        (&_S1072)->primal_0 = _S1054;
        (&_S1072)->differential_0 = _S1033;
        s_bwd_prop_cross_0(&_S1071, &_S1072, _S1070);
        float3  s_diff_dy_T_1 = - _S1072.differential_0;
        float3  _S1073 = - s_diff_dy_T_1;
        float3  _S1074 = - _S1071.differential_0;
        FixedArray<float3 , 4>  _S1075;
        _S1075[int(0)] = _S1033;
        _S1075[int(1)] = _S1033;
        _S1075[int(2)] = _S1033;
        _S1075[int(3)] = _S1033;
        _S1075[int(2)] = _S1073;
        _S1075[int(3)] = s_diff_dy_T_1;
        _S1075[int(0)] = _S1074;
        _S1075[int(1)] = _S1071.differential_0;
        points_4[int(0)] = _S1075[int(0)];
        points_4[int(1)] = _S1075[int(1)];
        points_4[int(2)] = _S1075[int(2)];
        points_4[int(3)] = _S1075[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1033;
        points_4[int(1)] = _S1033;
        points_4[int(2)] = _S1033;
        points_4[int(3)] = _S1033;
    }
    float4  _S1076;
    if(_S1034)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1077 = points_4;
                FixedArray<float3 , 4>  _S1078 = points_4;
                FixedArray<float3 , 4>  _S1079 = points_4;
                FixedArray<float3 , 4>  _S1080 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1081 = raydir_27 * _S1080[int(3)];
                    float _S1082 = _S1081.x + _S1081.y + _S1081.z;
                    float4  _S1083 = _S1062;
                    *&((&_S1083)->w) = _S1082;
                    points_4[int(0)] = _S1077[int(0)];
                    points_4[int(1)] = _S1078[int(1)];
                    points_4[int(2)] = _S1079[int(2)];
                    points_4[int(3)] = _S1033;
                    _S1076 = _S1083;
                }
                else
                {
                    points_4[int(0)] = _S1077[int(0)];
                    points_4[int(1)] = _S1078[int(1)];
                    points_4[int(2)] = _S1079[int(2)];
                    points_4[int(3)] = _S1080[int(3)];
                    _S1076 = _S1062;
                }
                float3  _S1084 = raydir_28 * points_4[int(2)];
                float _S1085 = _S1084.x + _S1084.y + _S1084.z;
                FixedArray<float3 , 4>  _S1086 = points_4;
                FixedArray<float3 , 4>  _S1087 = points_4;
                float4  _S1088 = _S1062;
                *&((&_S1088)->z) = _S1085;
                float4  _S1089 = _S1076 + _S1088;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1086[int(1)];
                points_4[int(2)] = _S1033;
                points_4[int(3)] = _S1087[int(3)];
                _S1076 = _S1089;
            }
            else
            {
                FixedArray<float3 , 4>  _S1090 = points_4;
                FixedArray<float3 , 4>  _S1091 = points_4;
                FixedArray<float3 , 4>  _S1092 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1090[int(1)];
                points_4[int(2)] = _S1091[int(2)];
                points_4[int(3)] = _S1092[int(3)];
                _S1076 = _S1062;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1093 = points_4;
            FixedArray<float3 , 4>  _S1094 = points_4;
            FixedArray<float3 , 4>  _S1095 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1093[int(1)];
            points_4[int(2)] = _S1094[int(2)];
            points_4[int(3)] = _S1095[int(3)];
            _S1076 = _S1062;
        }
        if(_S1036)
        {
            FixedArray<float3 , 4>  _S1096 = points_4;
            float3  _S1097 = raydir_29 * points_4[int(1)];
            float _S1098 = _S1097.x + _S1097.y + _S1097.z;
            float4  _S1099 = _S1062;
            *&((&_S1099)->y) = _S1098;
            float4  _S1100 = _S1076 + _S1099;
            points_4[int(0)] = _S1033;
            points_4[int(1)] = _S1033;
            points_4[int(2)] = _S1033;
            points_4[int(3)] = _S1033;
            raydir_27 = _S1096[int(0)];
            _S1076 = _S1100;
        }
        else
        {
            FixedArray<float3 , 4>  _S1101 = points_4;
            FixedArray<float3 , 4>  _S1102 = points_4;
            FixedArray<float3 , 4>  _S1103 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1101[int(1)];
            points_4[int(2)] = _S1102[int(2)];
            points_4[int(3)] = _S1103[int(3)];
            raydir_27 = _S1033;
        }
        float3  _S1104 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1105 = _S1104.x + _S1104.y + _S1104.z;
        float4  _S1106 = _S1062;
        *&((&_S1106)->x) = _S1105;
        _S1076 = _S1076 + _S1106;
    }
    else
    {
        _S1076 = _S1062;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1076;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1107, float4  _S1108, FixedArray<float, 10>  * _S1109, bool _S1110, bool _S1111, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1112, float3  _S1113)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1114;
    float3  _S1115 = s_primal_ctx_depth_to_normal_0(_S1107, _S1108, _S1109, _S1110, _S1111, (*_S1112).primal_0, &_S1114);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1116 = _S1114;
    s_bwd_prop_depth_to_normal_0(_S1107, _S1108, _S1109, _S1110, _S1111, _S1112, _S1113, &_S1116);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1117 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1117;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, dist_coeffs_11, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9)
{
    float2  _S1118 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_16 = _S1118;
    bool _S1119 = undistort_point_0(_S1118, dist_coeffs_12, int(12), &uv_16);
    if(!_S1119)
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
    float _S1120 = (F32_exp2(((*dpx_18).primal_0))) * 0.693145751953125f * dOut_21;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1120;
    return;
}

inline __device__ float3  apply_ppisp(float3  rgb_in_0, float2  pix_coord_0, float2  image_center_0, float2  img_size_0, FixedArray<float, 36>  * params_0)
{
    PPISPParams_0 p_0;
    (&p_0)->exposure_0 = (*params_0)[int(0)];
    (&(&p_0)->vignette_params_0[int(0)])->cx_0 = (*params_0)[int(1)];
    (&(&p_0)->vignette_params_0[int(0)])->cy_0 = (*params_0)[int(2)];
    (&(&p_0)->vignette_params_0[int(0)])->alpha0_0 = (*params_0)[int(3)];
    (&(&p_0)->vignette_params_0[int(0)])->alpha1_0 = (*params_0)[int(4)];
    (&(&p_0)->vignette_params_0[int(0)])->alpha2_0 = (*params_0)[int(5)];
    (&(&p_0)->vignette_params_0[int(1)])->cx_0 = (*params_0)[int(6)];
    (&(&p_0)->vignette_params_0[int(1)])->cy_0 = (*params_0)[int(7)];
    (&(&p_0)->vignette_params_0[int(1)])->alpha0_0 = (*params_0)[int(8)];
    (&(&p_0)->vignette_params_0[int(1)])->alpha1_0 = (*params_0)[int(9)];
    (&(&p_0)->vignette_params_0[int(1)])->alpha2_0 = (*params_0)[int(10)];
    (&(&p_0)->vignette_params_0[int(2)])->cx_0 = (*params_0)[int(11)];
    (&(&p_0)->vignette_params_0[int(2)])->cy_0 = (*params_0)[int(12)];
    (&(&p_0)->vignette_params_0[int(2)])->alpha0_0 = (*params_0)[int(13)];
    (&(&p_0)->vignette_params_0[int(2)])->alpha1_0 = (*params_0)[int(14)];
    (&(&p_0)->vignette_params_0[int(2)])->alpha2_0 = (*params_0)[int(15)];
    *&((&(&(&p_0)->color_params_0)->b_0)->x) = (*params_0)[int(16)];
    *&((&(&(&p_0)->color_params_0)->b_0)->y) = (*params_0)[int(17)];
    *&((&(&(&p_0)->color_params_0)->r_0)->x) = (*params_0)[int(18)];
    *&((&(&(&p_0)->color_params_0)->r_0)->y) = (*params_0)[int(19)];
    *&((&(&(&p_0)->color_params_0)->g_0)->x) = (*params_0)[int(20)];
    *&((&(&(&p_0)->color_params_0)->g_0)->y) = (*params_0)[int(21)];
    *&((&(&(&p_0)->color_params_0)->n_0)->x) = (*params_0)[int(22)];
    *&((&(&(&p_0)->color_params_0)->n_0)->y) = (*params_0)[int(23)];
    (&(&p_0)->crf_params_0[int(0)])->toe_0 = (*params_0)[int(24)];
    (&(&p_0)->crf_params_0[int(0)])->shoulder_0 = (*params_0)[int(25)];
    (&(&p_0)->crf_params_0[int(0)])->gamma_0 = (*params_0)[int(26)];
    (&(&p_0)->crf_params_0[int(0)])->center_0 = (*params_0)[int(27)];
    (&(&p_0)->crf_params_0[int(1)])->toe_0 = (*params_0)[int(28)];
    (&(&p_0)->crf_params_0[int(1)])->shoulder_0 = (*params_0)[int(29)];
    (&(&p_0)->crf_params_0[int(1)])->gamma_0 = (*params_0)[int(30)];
    (&(&p_0)->crf_params_0[int(1)])->center_0 = (*params_0)[int(31)];
    (&(&p_0)->crf_params_0[int(2)])->toe_0 = (*params_0)[int(32)];
    (&(&p_0)->crf_params_0[int(2)])->shoulder_0 = (*params_0)[int(33)];
    (&(&p_0)->crf_params_0[int(2)])->gamma_0 = (*params_0)[int(34)];
    (&(&p_0)->crf_params_0[int(2)])->center_0 = (*params_0)[int(35)];
    PPISPParams_0 _S1121 = p_0;
    float max_res_0 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1122 = (pix_coord_0.x - image_center_0.x) / max_res_0;
    float _S1123 = (pix_coord_0.y - image_center_0.y) / max_res_0;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_0))));
    float dx_2 = _S1122 - p_0.vignette_params_0[int(0)].cx_0;
    float dy_0 = _S1123 - p_0.vignette_params_0[int(0)].cy_0;
    float r2_4 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_4 * r2_4;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_0[int(0)].alpha2_0 * (r4_0 * r2_4) + p_0.vignette_params_0[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_0[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1122 - p_0.vignette_params_0[int(1)].cx_0;
    float dy_1 = _S1123 - p_0.vignette_params_0[int(1)].cy_0;
    float r2_5 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_5 * r2_5;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_0[int(1)].alpha2_0 * (r4_1 * r2_5) + p_0.vignette_params_0[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_0[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1122 - p_0.vignette_params_0[int(2)].cx_0;
    float dy_2 = _S1123 - p_0.vignette_params_0[int(2)].cy_0;
    float r2_6 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_6 * r2_6;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_0[int(2)].alpha2_0 * (r4_2 * r2_6) + p_0.vignette_params_0[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_0[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S1124 = rgb_out_0;
    float2  bd_0 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_0.b_0);
    float2  rd_0 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_0.r_0);
    float2  gd_0 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_0.g_0);
    float2  nd_0 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_0.n_0);
    float _S1125 = 0.3333333432674408f + nd_0.x;
    float _S1126 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1126, 1.0f, 0.0f, - _S1125, - _S1126, _S1125, 0.0f), T_0);
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
    float _S1127 = _S1124.x;
    float _S1128 = _S1124.y;
    float intensity_0 = _S1127 + _S1128 + _S1124.z;
    float3  rgi_out_0 = mul_0(H_1, make_float3 (_S1127, _S1128, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1129 = rgi_out_1.x;
    float _S1130 = rgi_out_1.y;
    float3  _S1131 = clamp_1(make_float3 (_S1129, _S1130, rgi_out_1.z - _S1129 - _S1130), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1132 = _S1131.x;
    float _S1133 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(0)].toe_0))))));
    float _S1134 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(0)].shoulder_0))))));
    float _S1135 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(0)].gamma_0))))));
    float _S1136 = 1.0f / (1.0f + (F32_exp((- _S1121.crf_params_0[int(0)].center_0))));
    float a_1 = _S1134 * _S1136 / lerp_0(_S1133, _S1134, _S1136);
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1132 <= _S1136)
    {
        y_10 = a_1 * (F32_pow((_S1132 / _S1136), (_S1133)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1132) / (1.0f - _S1136)), (_S1134)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1135)));
    float _S1137 = _S1131.y;
    float _S1138 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(1)].toe_0))))));
    float _S1139 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(1)].shoulder_0))))));
    float _S1140 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(1)].gamma_0))))));
    float _S1141 = 1.0f / (1.0f + (F32_exp((- _S1121.crf_params_0[int(1)].center_0))));
    float a_2 = _S1139 * _S1141 / lerp_0(_S1138, _S1139, _S1141);
    float b_3 = 1.0f - a_2;
    if(_S1137 <= _S1141)
    {
        y_10 = a_2 * (F32_pow((_S1137 / _S1141), (_S1138)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1137) / (1.0f - _S1141)), (_S1139)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1140)));
    float _S1142 = _S1131.z;
    float _S1143 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(2)].toe_0))))));
    float _S1144 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(2)].shoulder_0))))));
    float _S1145 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1121.crf_params_0[int(2)].gamma_0))))));
    float _S1146 = 1.0f / (1.0f + (F32_exp((- _S1121.crf_params_0[int(2)].center_0))));
    float a_3 = _S1144 * _S1146 / lerp_0(_S1143, _S1144, _S1146);
    float b_4 = 1.0f - a_3;
    if(_S1142 <= _S1146)
    {
        y_10 = a_3 * (F32_pow((_S1142 / _S1146), (_S1143)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1142) / (1.0f - _S1146)), (_S1144)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1145)));
    return rgb_out_1;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1147)
{
    return (F32_exp2((_S1147)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1148, float2  _S1149)
{
    return mul_1(_S1148, _S1149);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1150, Matrix<float, 3, 3>  _S1151)
{
    return mul_3(_S1150, _S1151);
}

inline __device__ float s_primal_ctx_abs_0(float _S1152)
{
    return (F32_abs((_S1152)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1153, float3  _S1154)
{
    return mul_0(_S1153, _S1154);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1155, float3  _S1156, float3  _S1157)
{
    return clamp_1(_S1155, _S1156, _S1157);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1158, float _S1159, float _S1160)
{
    return lerp_0(_S1158, _S1159, _S1160);
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1161, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1162, float3  _S1163)
{
    _d_mul_0(_S1161, _S1162, _S1163);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1164, float _S1165)
{
    _d_abs_0(_S1164, _S1165);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1166, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1167, Matrix<float, 3, 3>  _S1168)
{
    mul_2(_S1166, _S1167, _S1168);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1169, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1170, float2  _S1171)
{
    _d_mul_1(_S1169, _S1170, _S1171);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1172, float _S1173)
{
    _d_exp2_0(_S1172, _S1173);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1174 = *dprgb_in_0;
    float3  _S1175 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1176 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1177 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1178 = {
        _S1177, _S1177, _S1177
    };
    float2  _S1179 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1180 = { _S1179, _S1179, _S1179, _S1179 };
    CRFPPISPChannelParams_0 _S1181 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1182 = {
        _S1181, _S1181, _S1181
    };
    PPISPParams_0 _S1183;
    (&_S1183)->exposure_0 = dpparams_0->primal_0[int(0)];
    (&_S1183)->vignette_params_0 = _S1178;
    (&_S1183)->color_params_0 = _S1180;
    (&_S1183)->crf_params_0 = _S1182;
    (&(&_S1183)->vignette_params_0[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1183)->vignette_params_0[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1184 = dpparams_0->primal_0[int(3)];
    (&(&_S1183)->vignette_params_0[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1185 = dpparams_0->primal_0[int(4)];
    (&(&_S1183)->vignette_params_0[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1186 = dpparams_0->primal_0[int(5)];
    (&(&_S1183)->vignette_params_0[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1183)->vignette_params_0[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1183)->vignette_params_0[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1187 = dpparams_0->primal_0[int(8)];
    (&(&_S1183)->vignette_params_0[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1188 = dpparams_0->primal_0[int(9)];
    (&(&_S1183)->vignette_params_0[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1189 = dpparams_0->primal_0[int(10)];
    (&(&_S1183)->vignette_params_0[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1183)->vignette_params_0[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1183)->vignette_params_0[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1190 = dpparams_0->primal_0[int(13)];
    (&(&_S1183)->vignette_params_0[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1191 = dpparams_0->primal_0[int(14)];
    (&(&_S1183)->vignette_params_0[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1192 = dpparams_0->primal_0[int(15)];
    (&(&_S1183)->vignette_params_0[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1183)->color_params_0)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1183)->color_params_0)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1183)->color_params_0)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1183)->color_params_0)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1183)->color_params_0)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1183)->color_params_0)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1183)->color_params_0)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1183)->color_params_0)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1193 = dpparams_0->primal_0[int(24)];
    (&(&_S1183)->crf_params_0[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1194 = dpparams_0->primal_0[int(25)];
    (&(&_S1183)->crf_params_0[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1195 = dpparams_0->primal_0[int(26)];
    (&(&_S1183)->crf_params_0[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1196 = dpparams_0->primal_0[int(27)];
    (&(&_S1183)->crf_params_0[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1197 = dpparams_0->primal_0[int(28)];
    (&(&_S1183)->crf_params_0[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1198 = dpparams_0->primal_0[int(29)];
    (&(&_S1183)->crf_params_0[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1199 = dpparams_0->primal_0[int(30)];
    (&(&_S1183)->crf_params_0[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1200 = dpparams_0->primal_0[int(31)];
    (&(&_S1183)->crf_params_0[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1201 = dpparams_0->primal_0[int(32)];
    (&(&_S1183)->crf_params_0[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1202 = dpparams_0->primal_0[int(33)];
    (&(&_S1183)->crf_params_0[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1203 = dpparams_0->primal_0[int(34)];
    (&(&_S1183)->crf_params_0[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1204 = dpparams_0->primal_0[int(35)];
    (&(&_S1183)->crf_params_0[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1205 = _S1183;
    float _S1206 = s_primal_ctx_exp2_0(_S1183.exposure_0);
    float3  _S1207 = make_float3 (_S1206);
    float3  rgb_out_2 = (*dprgb_in_0).primal_0 * make_float3 (_S1206);
    float _S1208 = s_primal_ctx_max_0(img_size_1.x, img_size_1.y);
    float _S1209 = (pix_coord_1.x - image_center_1.x) / _S1208;
    float _S1210 = (pix_coord_1.y - image_center_1.y) / _S1208;
    float dx_5 = _S1209 - dpparams_0->primal_0[int(1)];
    float dy_3 = _S1210 - dpparams_0->primal_0[int(2)];
    float r2_8 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_8 * r2_8;
    float r6_0 = r4_3 * r2_8;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_3 + dpparams_0->primal_0[int(3)] * r2_8 + 1.0f;
    float _S1211 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1212 = rgb_out_2.x * _S1211;
    float3  _S1213 = rgb_out_2;
    *&((&_S1213)->x) = _S1212;
    float dx_6 = _S1209 - dpparams_0->primal_0[int(6)];
    float dy_4 = _S1210 - dpparams_0->primal_0[int(7)];
    float r2_9 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_9 * r2_9;
    float r6_1 = r4_4 * r2_9;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_4 + dpparams_0->primal_0[int(8)] * r2_9 + 1.0f;
    float _S1214 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1213)->y) = rgb_out_2.y * _S1214;
    float dx_7 = _S1209 - dpparams_0->primal_0[int(11)];
    float dy_5 = _S1210 - dpparams_0->primal_0[int(12)];
    float r2_10 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_10 * r2_10;
    float r6_2 = r4_5 * r2_10;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_5 + dpparams_0->primal_0[int(13)] * r2_10 + 1.0f;
    float _S1215 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1213)->z) = rgb_out_2.z * _S1215;
    PPISPParams_0 _S1216 = _S1183;
    float2  _S1217 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1183.color_params_0.b_0);
    float2  _S1218 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1183.color_params_0.r_0);
    float2  _S1219 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1183.color_params_0.g_0);
    float2  _S1220 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1183.color_params_0.n_0);
    float _S1221 = 0.3333333432674408f + _S1220.x;
    float _S1222 = 0.3333333432674408f + _S1220.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (_S1217.x, 1.0f + _S1218.x, _S1219.x, _S1217.y, _S1218.y, 1.0f + _S1219.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1222, 1.0f, 0.0f, - _S1221, - _S1222, _S1221, 0.0f);
    Matrix<float, 3, 3>  _S1223 = s_primal_ctx_mul_1(skew_0, T_1);
    float3  r0_1 = make_float3 (_S1223.rows[int(0)].x, _S1223.rows[int(0)].y, _S1223.rows[int(0)].z);
    float3  r1_1 = make_float3 (_S1223.rows[int(1)].x, _S1223.rows[int(1)].y, _S1223.rows[int(1)].z);
    float3  r2_11 = make_float3 (_S1223.rows[int(2)].x, _S1223.rows[int(2)].y, _S1223.rows[int(2)].z);
    float3  _S1224 = s_primal_ctx_cross_0(r0_1, r1_1);
    bool _S1225 = (s_primal_ctx_dot_0(_S1224, _S1224)) < 9.99999968265522539e-21f;
    float3  lambda_v_3;
    float3  _S1226;
    bool _S1227;
    if(_S1225)
    {
        float3  _S1228 = s_primal_ctx_cross_0(r0_1, r2_11);
        bool _S1229 = (s_primal_ctx_dot_0(_S1228, _S1228)) < 9.99999968265522539e-21f;
        if(_S1229)
        {
            lambda_v_3 = s_primal_ctx_cross_0(r1_1, r2_11);
        }
        else
        {
            lambda_v_3 = _S1228;
        }
        _S1227 = _S1229;
        _S1226 = _S1228;
    }
    else
    {
        lambda_v_3 = _S1224;
        _S1227 = false;
        _S1226 = _S1175;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_3.x, 0.0f, 0.0f, 0.0f, lambda_v_3.y, 0.0f, 0.0f, 0.0f, lambda_v_3.z);
    Matrix<float, 3, 3>  _S1230 = s_primal_ctx_mul_1(T_1, D_0);
    Matrix<float, 3, 3>  _S1231 = s_primal_ctx_mul_1(_S1230, S_inv_0);
    bool _S1232 = (s_primal_ctx_abs_0(_S1231.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_2;
    Matrix<float, 3, 3>  _S1233;
    float _S1234;
    if(_S1232)
    {
        float inv_s_0 = 1.0f / _S1231.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1235 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1236 = _S1231.rows[int(2)].z * _S1231.rows[int(2)].z;
        H_2 = _S1231 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1233 = _S1235;
        _S1234 = _S1236;
    }
    else
    {
        H_2 = _S1231;
        _S1233 = _S1176;
        _S1234 = 0.0f;
    }
    float _S1237 = _S1213.x;
    float _S1238 = _S1213.y;
    float intensity_1 = _S1237 + _S1238 + _S1213.z;
    float3  rgi_in_0 = make_float3 (_S1237, _S1238, intensity_1);
    float3  _S1239 = s_primal_ctx_mul_2(H_2, rgi_in_0);
    float _S1240 = _S1239.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_1 / _S1240;
    float3  _S1241 = make_float3 (norm_factor_0);
    float _S1242 = _S1240 * _S1240;
    float3  rgi_out_2 = _S1239 * make_float3 (norm_factor_0);
    float _S1243 = rgi_out_2.x;
    float _S1244 = rgi_out_2.y;
    float3  _S1245 = make_float3 (_S1243, _S1244, rgi_out_2.z - _S1243 - _S1244);
    float3  _S1246 = make_float3 (0.0f);
    float3  _S1247 = make_float3 (1.0f);
    float3  _S1248 = s_primal_ctx_clamp_1(_S1245, _S1246, _S1247);
    float _S1249 = _S1248.x;
    float _S1250 = 1.0f + s_primal_ctx_exp_0(_S1193);
    float _S1251 = 0.30000001192092896f + s_primal_ctx_log_0(_S1250);
    float _S1252 = 1.0f + s_primal_ctx_exp_0(_S1194);
    float _S1253 = 0.30000001192092896f + s_primal_ctx_log_0(_S1252);
    float _S1254 = 1.0f + s_primal_ctx_exp_0(_S1195);
    float _S1255 = 0.10000000149011612f + s_primal_ctx_log_0(_S1254);
    float _S1256 = - _S1196;
    float _S1257 = 1.0f + s_primal_ctx_exp_0(_S1256);
    float _S1258 = 1.0f / _S1257;
    float _S1259 = _S1257 * _S1257;
    float _S1260 = s_primal_ctx_lerp_0(_S1251, _S1253, _S1258);
    float _S1261 = _S1253 * _S1258;
    float a_4 = _S1261 / _S1260;
    float _S1262 = _S1260 * _S1260;
    float b_5 = 1.0f - a_4;
    bool _S1263 = _S1249 <= _S1258;
    float y_11;
    float _S1264;
    float _S1265;
    float _S1266;
    float _S1267;
    float _S1268;
    float _S1269;
    float _S1270;
    float _S1271;
    if(_S1263)
    {
        float _S1272 = _S1249 / _S1258;
        float _S1273 = _S1258 * _S1258;
        float _S1274 = s_primal_ctx_pow_0(_S1272, _S1251);
        y_11 = a_4 * _S1274;
        _S1264 = 0.0f;
        _S1265 = 0.0f;
        _S1266 = 0.0f;
        _S1267 = 0.0f;
        _S1268 = 0.0f;
        _S1269 = _S1274;
        _S1270 = _S1272;
        _S1271 = _S1273;
    }
    else
    {
        float _S1275 = 1.0f - _S1249;
        float _S1276 = 1.0f - _S1258;
        float _S1277 = _S1275 / _S1276;
        float _S1278 = _S1276 * _S1276;
        float _S1279 = s_primal_ctx_pow_0(_S1277, _S1253);
        y_11 = 1.0f - b_5 * _S1279;
        _S1264 = _S1279;
        _S1265 = _S1277;
        _S1266 = _S1278;
        _S1267 = _S1275;
        _S1268 = _S1276;
        _S1269 = 0.0f;
        _S1270 = 0.0f;
        _S1271 = 0.0f;
    }
    float _S1280 = s_primal_ctx_max_0(0.0f, y_11);
    float _S1281 = _S1248.y;
    float _S1282 = 1.0f + s_primal_ctx_exp_0(_S1197);
    float _S1283 = 0.30000001192092896f + s_primal_ctx_log_0(_S1282);
    float _S1284 = 1.0f + s_primal_ctx_exp_0(_S1198);
    float _S1285 = 0.30000001192092896f + s_primal_ctx_log_0(_S1284);
    float _S1286 = 1.0f + s_primal_ctx_exp_0(_S1199);
    float _S1287 = 0.10000000149011612f + s_primal_ctx_log_0(_S1286);
    float _S1288 = - _S1200;
    float _S1289 = 1.0f + s_primal_ctx_exp_0(_S1288);
    float _S1290 = 1.0f / _S1289;
    float _S1291 = _S1289 * _S1289;
    float _S1292 = s_primal_ctx_lerp_0(_S1283, _S1285, _S1290);
    float _S1293 = _S1285 * _S1290;
    float a_5 = _S1293 / _S1292;
    float _S1294 = _S1292 * _S1292;
    float b_6 = 1.0f - a_5;
    bool _S1295 = _S1281 <= _S1290;
    float y_12;
    float _S1296;
    float _S1297;
    float _S1298;
    float _S1299;
    float _S1300;
    float _S1301;
    float _S1302;
    float _S1303;
    if(_S1295)
    {
        float _S1304 = _S1281 / _S1290;
        float _S1305 = _S1290 * _S1290;
        float _S1306 = s_primal_ctx_pow_0(_S1304, _S1283);
        y_12 = a_5 * _S1306;
        _S1296 = 0.0f;
        _S1297 = 0.0f;
        _S1298 = 0.0f;
        _S1299 = 0.0f;
        _S1300 = 0.0f;
        _S1301 = _S1306;
        _S1302 = _S1304;
        _S1303 = _S1305;
    }
    else
    {
        float _S1307 = 1.0f - _S1281;
        float _S1308 = 1.0f - _S1290;
        float _S1309 = _S1307 / _S1308;
        float _S1310 = _S1308 * _S1308;
        float _S1311 = s_primal_ctx_pow_0(_S1309, _S1285);
        y_12 = 1.0f - b_6 * _S1311;
        _S1296 = _S1311;
        _S1297 = _S1309;
        _S1298 = _S1310;
        _S1299 = _S1307;
        _S1300 = _S1308;
        _S1301 = 0.0f;
        _S1302 = 0.0f;
        _S1303 = 0.0f;
    }
    float _S1312 = s_primal_ctx_max_0(0.0f, y_12);
    float _S1313 = _S1248.z;
    float _S1314 = 1.0f + s_primal_ctx_exp_0(_S1201);
    float _S1315 = 0.30000001192092896f + s_primal_ctx_log_0(_S1314);
    float _S1316 = 1.0f + s_primal_ctx_exp_0(_S1202);
    float _S1317 = 0.30000001192092896f + s_primal_ctx_log_0(_S1316);
    float _S1318 = 1.0f + s_primal_ctx_exp_0(_S1203);
    float _S1319 = 0.10000000149011612f + s_primal_ctx_log_0(_S1318);
    float _S1320 = - _S1204;
    float _S1321 = 1.0f + s_primal_ctx_exp_0(_S1320);
    float _S1322 = 1.0f / _S1321;
    float _S1323 = _S1321 * _S1321;
    float _S1324 = s_primal_ctx_lerp_0(_S1315, _S1317, _S1322);
    float _S1325 = _S1317 * _S1322;
    float a_6 = _S1325 / _S1324;
    float _S1326 = _S1324 * _S1324;
    float b_7 = 1.0f - a_6;
    bool _S1327 = _S1313 <= _S1322;
    float y_13;
    float _S1328;
    float _S1329;
    float _S1330;
    float _S1331;
    float _S1332;
    float _S1333;
    float _S1334;
    float _S1335;
    if(_S1327)
    {
        float _S1336 = _S1313 / _S1322;
        float _S1337 = _S1322 * _S1322;
        float _S1338 = s_primal_ctx_pow_0(_S1336, _S1315);
        y_13 = a_6 * _S1338;
        _S1328 = 0.0f;
        _S1329 = 0.0f;
        _S1330 = 0.0f;
        _S1331 = 0.0f;
        _S1332 = 0.0f;
        _S1333 = _S1338;
        _S1334 = _S1336;
        _S1335 = _S1337;
    }
    else
    {
        float _S1339 = 1.0f - _S1313;
        float _S1340 = 1.0f - _S1322;
        float _S1341 = _S1339 / _S1340;
        float _S1342 = _S1340 * _S1340;
        float _S1343 = s_primal_ctx_pow_0(_S1341, _S1317);
        y_13 = 1.0f - b_7 * _S1343;
        _S1328 = _S1343;
        _S1329 = _S1341;
        _S1330 = _S1342;
        _S1331 = _S1339;
        _S1332 = _S1340;
        _S1333 = 0.0f;
        _S1334 = 0.0f;
        _S1335 = 0.0f;
    }
    float _S1344 = s_primal_ctx_max_0(0.0f, y_13);
    DiffPair_float_0 _S1345;
    (&_S1345)->primal_0 = _S1344;
    (&_S1345)->differential_0 = 0.0f;
    DiffPair_float_0 _S1346;
    (&_S1346)->primal_0 = _S1319;
    (&_S1346)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1345, &_S1346, _s_dOut_11.z);
    DiffPair_float_0 _S1347 = _S1346;
    DiffPair_float_0 _S1348;
    (&_S1348)->primal_0 = 0.0f;
    (&_S1348)->differential_0 = 0.0f;
    DiffPair_float_0 _S1349;
    (&_S1349)->primal_0 = y_13;
    (&_S1349)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1348, &_S1349, _S1345.differential_0);
    DiffPair_float_0 _S1350 = _S1349;
    if(_S1327)
    {
        float _S1351 = a_6 * _S1350.differential_0;
        float _S1352 = _S1333 * _S1350.differential_0;
        DiffPair_float_0 _S1353;
        (&_S1353)->primal_0 = _S1334;
        (&_S1353)->differential_0 = 0.0f;
        DiffPair_float_0 _S1354;
        (&_S1354)->primal_0 = _S1315;
        (&_S1354)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1353, &_S1354, _S1351);
        float _S1355 = _S1353.differential_0 / _S1335;
        float _S1356 = _S1313 * - _S1355;
        float _S1357 = _S1322 * _S1355;
        y_13 = 0.0f;
        _S1328 = _S1352;
        _S1329 = _S1356;
        _S1330 = 0.0f;
        _S1331 = _S1354.differential_0;
        _S1332 = _S1357;
    }
    else
    {
        float _S1358 = - _S1350.differential_0;
        float _S1359 = b_7 * _S1358;
        float _S1360 = _S1328 * _S1358;
        DiffPair_float_0 _S1361;
        (&_S1361)->primal_0 = _S1329;
        (&_S1361)->differential_0 = 0.0f;
        DiffPair_float_0 _S1362;
        (&_S1362)->primal_0 = _S1317;
        (&_S1362)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1361, &_S1362, _S1359);
        float _S1363 = _S1361.differential_0 / _S1330;
        float _S1364 = - (_S1331 * - _S1363);
        float _S1365 = - (_S1332 * _S1363);
        y_13 = _S1360;
        _S1328 = 0.0f;
        _S1329 = _S1364;
        _S1330 = _S1362.differential_0;
        _S1331 = 0.0f;
        _S1332 = _S1365;
    }
    float _S1366 = (- y_13 + _S1328) / _S1326;
    float _S1367 = _S1325 * - _S1366;
    float _S1368 = _S1324 * _S1366;
    float _S1369 = _S1317 * _S1368;
    float _S1370 = _S1322 * _S1368;
    DiffPair_float_0 _S1371;
    (&_S1371)->primal_0 = _S1315;
    (&_S1371)->differential_0 = 0.0f;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1317;
    (&_S1372)->differential_0 = 0.0f;
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S1322;
    (&_S1373)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1371, &_S1372, &_S1373, _S1367);
    float _S1374 = - ((_S1369 + _S1373.differential_0 + _S1329) / _S1323);
    DiffPair_float_0 _S1375;
    (&_S1375)->primal_0 = _S1320;
    (&_S1375)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1375, _S1374);
    float _S1376 = - _S1375.differential_0;
    DiffPair_float_0 _S1377;
    (&_S1377)->primal_0 = _S1318;
    (&_S1377)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1377, _S1347.differential_0);
    DiffPair_float_0 _S1378;
    (&_S1378)->primal_0 = _S1203;
    (&_S1378)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1378, _S1377.differential_0);
    DiffPair_float_0 _S1379 = _S1378;
    float _S1380 = _S1370 + _S1372.differential_0 + _S1330;
    DiffPair_float_0 _S1381;
    (&_S1381)->primal_0 = _S1316;
    (&_S1381)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1381, _S1380);
    DiffPair_float_0 _S1382;
    (&_S1382)->primal_0 = _S1202;
    (&_S1382)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1382, _S1381.differential_0);
    DiffPair_float_0 _S1383 = _S1382;
    float _S1384 = _S1371.differential_0 + _S1331;
    DiffPair_float_0 _S1385;
    (&_S1385)->primal_0 = _S1314;
    (&_S1385)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1385, _S1384);
    DiffPair_float_0 _S1386;
    (&_S1386)->primal_0 = _S1201;
    (&_S1386)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1386, _S1385.differential_0);
    DiffPair_float_0 _S1387 = _S1386;
    float3  _S1388 = make_float3 (0.0f, 0.0f, _S1332);
    DiffPair_float_0 _S1389;
    (&_S1389)->primal_0 = _S1312;
    (&_S1389)->differential_0 = 0.0f;
    DiffPair_float_0 _S1390;
    (&_S1390)->primal_0 = _S1287;
    (&_S1390)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1389, &_S1390, _s_dOut_11.y);
    DiffPair_float_0 _S1391 = _S1390;
    DiffPair_float_0 _S1392;
    (&_S1392)->primal_0 = 0.0f;
    (&_S1392)->differential_0 = 0.0f;
    DiffPair_float_0 _S1393;
    (&_S1393)->primal_0 = y_12;
    (&_S1393)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1392, &_S1393, _S1389.differential_0);
    DiffPair_float_0 _S1394 = _S1393;
    if(_S1295)
    {
        float _S1395 = a_5 * _S1394.differential_0;
        float _S1396 = _S1301 * _S1394.differential_0;
        DiffPair_float_0 _S1397;
        (&_S1397)->primal_0 = _S1302;
        (&_S1397)->differential_0 = 0.0f;
        DiffPair_float_0 _S1398;
        (&_S1398)->primal_0 = _S1283;
        (&_S1398)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1397, &_S1398, _S1395);
        float _S1399 = _S1397.differential_0 / _S1303;
        float _S1400 = _S1281 * - _S1399;
        float _S1401 = _S1290 * _S1399;
        y_12 = 0.0f;
        _S1296 = _S1396;
        _S1297 = _S1400;
        _S1298 = 0.0f;
        _S1299 = _S1398.differential_0;
        _S1300 = _S1401;
    }
    else
    {
        float _S1402 = - _S1394.differential_0;
        float _S1403 = b_6 * _S1402;
        float _S1404 = _S1296 * _S1402;
        DiffPair_float_0 _S1405;
        (&_S1405)->primal_0 = _S1297;
        (&_S1405)->differential_0 = 0.0f;
        DiffPair_float_0 _S1406;
        (&_S1406)->primal_0 = _S1285;
        (&_S1406)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1405, &_S1406, _S1403);
        float _S1407 = _S1405.differential_0 / _S1298;
        float _S1408 = - (_S1299 * - _S1407);
        float _S1409 = - (_S1300 * _S1407);
        y_12 = _S1404;
        _S1296 = 0.0f;
        _S1297 = _S1408;
        _S1298 = _S1406.differential_0;
        _S1299 = 0.0f;
        _S1300 = _S1409;
    }
    float _S1410 = (- y_12 + _S1296) / _S1294;
    float _S1411 = _S1293 * - _S1410;
    float _S1412 = _S1292 * _S1410;
    float _S1413 = _S1285 * _S1412;
    float _S1414 = _S1290 * _S1412;
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1283;
    (&_S1415)->differential_0 = 0.0f;
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1285;
    (&_S1416)->differential_0 = 0.0f;
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1290;
    (&_S1417)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1415, &_S1416, &_S1417, _S1411);
    float _S1418 = - ((_S1413 + _S1417.differential_0 + _S1297) / _S1291);
    DiffPair_float_0 _S1419;
    (&_S1419)->primal_0 = _S1288;
    (&_S1419)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1419, _S1418);
    float _S1420 = - _S1419.differential_0;
    DiffPair_float_0 _S1421;
    (&_S1421)->primal_0 = _S1286;
    (&_S1421)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1421, _S1391.differential_0);
    DiffPair_float_0 _S1422;
    (&_S1422)->primal_0 = _S1199;
    (&_S1422)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1422, _S1421.differential_0);
    DiffPair_float_0 _S1423 = _S1422;
    float _S1424 = _S1414 + _S1416.differential_0 + _S1298;
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = _S1284;
    (&_S1425)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1425, _S1424);
    DiffPair_float_0 _S1426;
    (&_S1426)->primal_0 = _S1198;
    (&_S1426)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1426, _S1425.differential_0);
    DiffPair_float_0 _S1427 = _S1426;
    float _S1428 = _S1415.differential_0 + _S1299;
    DiffPair_float_0 _S1429;
    (&_S1429)->primal_0 = _S1282;
    (&_S1429)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1429, _S1428);
    DiffPair_float_0 _S1430;
    (&_S1430)->primal_0 = _S1197;
    (&_S1430)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1430, _S1429.differential_0);
    DiffPair_float_0 _S1431 = _S1430;
    float3  _S1432 = _S1388 + make_float3 (0.0f, _S1300, 0.0f);
    DiffPair_float_0 _S1433;
    (&_S1433)->primal_0 = _S1280;
    (&_S1433)->differential_0 = 0.0f;
    DiffPair_float_0 _S1434;
    (&_S1434)->primal_0 = _S1255;
    (&_S1434)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1433, &_S1434, _s_dOut_11.x);
    DiffPair_float_0 _S1435 = _S1434;
    DiffPair_float_0 _S1436;
    (&_S1436)->primal_0 = 0.0f;
    (&_S1436)->differential_0 = 0.0f;
    DiffPair_float_0 _S1437;
    (&_S1437)->primal_0 = y_11;
    (&_S1437)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1436, &_S1437, _S1433.differential_0);
    DiffPair_float_0 _S1438 = _S1437;
    if(_S1263)
    {
        float _S1439 = a_4 * _S1438.differential_0;
        float _S1440 = _S1269 * _S1438.differential_0;
        DiffPair_float_0 _S1441;
        (&_S1441)->primal_0 = _S1270;
        (&_S1441)->differential_0 = 0.0f;
        DiffPair_float_0 _S1442;
        (&_S1442)->primal_0 = _S1251;
        (&_S1442)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1441, &_S1442, _S1439);
        float _S1443 = _S1441.differential_0 / _S1271;
        float _S1444 = _S1249 * - _S1443;
        float _S1445 = _S1258 * _S1443;
        y_11 = 0.0f;
        _S1264 = _S1440;
        _S1265 = _S1444;
        _S1266 = 0.0f;
        _S1267 = _S1442.differential_0;
        _S1268 = _S1445;
    }
    else
    {
        float _S1446 = - _S1438.differential_0;
        float _S1447 = b_5 * _S1446;
        float _S1448 = _S1264 * _S1446;
        DiffPair_float_0 _S1449;
        (&_S1449)->primal_0 = _S1265;
        (&_S1449)->differential_0 = 0.0f;
        DiffPair_float_0 _S1450;
        (&_S1450)->primal_0 = _S1253;
        (&_S1450)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1449, &_S1450, _S1447);
        float _S1451 = _S1449.differential_0 / _S1266;
        float _S1452 = - (_S1267 * - _S1451);
        float _S1453 = - (_S1268 * _S1451);
        y_11 = _S1448;
        _S1264 = 0.0f;
        _S1265 = _S1452;
        _S1266 = _S1450.differential_0;
        _S1267 = 0.0f;
        _S1268 = _S1453;
    }
    float _S1454 = (- y_11 + _S1264) / _S1262;
    float _S1455 = _S1261 * - _S1454;
    float _S1456 = _S1260 * _S1454;
    float _S1457 = _S1253 * _S1456;
    float _S1458 = _S1258 * _S1456;
    DiffPair_float_0 _S1459;
    (&_S1459)->primal_0 = _S1251;
    (&_S1459)->differential_0 = 0.0f;
    DiffPair_float_0 _S1460;
    (&_S1460)->primal_0 = _S1253;
    (&_S1460)->differential_0 = 0.0f;
    DiffPair_float_0 _S1461;
    (&_S1461)->primal_0 = _S1258;
    (&_S1461)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1459, &_S1460, &_S1461, _S1455);
    float _S1462 = - ((_S1457 + _S1461.differential_0 + _S1265) / _S1259);
    DiffPair_float_0 _S1463;
    (&_S1463)->primal_0 = _S1256;
    (&_S1463)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1463, _S1462);
    float _S1464 = - _S1463.differential_0;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = _S1254;
    (&_S1465)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1465, _S1435.differential_0);
    DiffPair_float_0 _S1466;
    (&_S1466)->primal_0 = _S1195;
    (&_S1466)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1466, _S1465.differential_0);
    DiffPair_float_0 _S1467 = _S1466;
    float _S1468 = _S1458 + _S1460.differential_0 + _S1266;
    DiffPair_float_0 _S1469;
    (&_S1469)->primal_0 = _S1252;
    (&_S1469)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1469, _S1468);
    DiffPair_float_0 _S1470;
    (&_S1470)->primal_0 = _S1194;
    (&_S1470)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1470, _S1469.differential_0);
    DiffPair_float_0 _S1471 = _S1470;
    float _S1472 = _S1459.differential_0 + _S1267;
    DiffPair_float_0 _S1473;
    (&_S1473)->primal_0 = _S1250;
    (&_S1473)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1473, _S1472);
    DiffPair_float_0 _S1474;
    (&_S1474)->primal_0 = _S1193;
    (&_S1474)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1474, _S1473.differential_0);
    DiffPair_float_0 _S1475 = _S1474;
    float3  _S1476 = _S1432 + make_float3 (_S1268, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1477;
    (&_S1477)->primal_0 = _S1245;
    (&_S1477)->differential_0 = _S1175;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1478;
    (&_S1478)->primal_0 = _S1246;
    (&_S1478)->differential_0 = _S1175;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1479;
    (&_S1479)->primal_0 = _S1247;
    (&_S1479)->differential_0 = _S1175;
    s_bwd_prop_clamp_1(&_S1477, &_S1478, &_S1479, _S1476);
    float _S1480 = - _S1477.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1477.differential_0.x + _S1480, _S1477.differential_0.y + _S1480, _S1477.differential_0.z);
    float3  _S1481 = _S1239 * s_diff_rgi_out_T_0;
    float _S1482 = (_S1481.x + _S1481.y + _S1481.z) / _S1242;
    float _S1483 = _S1240 * _S1482;
    float3  _S1484 = _S1241 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_1 * - _S1482);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1485;
    (&_S1485)->primal_0 = H_2;
    (&_S1485)->differential_0 = _S1176;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1486;
    (&_S1486)->primal_0 = rgi_in_0;
    (&_S1486)->differential_0 = _S1175;
    s_bwd_prop_mul_1(&_S1485, &_S1486, _S1484);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1487 = _S1485;
    float _S1488 = _S1483 + _S1486.differential_0.z;
    float _S1489 = _S1486.differential_0.y + _S1488;
    float _S1490 = _S1486.differential_0.x + _S1488;
    float3  _S1491 = make_float3 (_S1490, _S1489, _S1488);
    if(_S1232)
    {
        Matrix<float, 3, 3>  _S1492 = _S1231 * _S1487.differential_0;
        Matrix<float, 3, 3>  _S1493 = _S1233 * _S1487.differential_0;
        _S1234 = - ((_S1492.rows[int(0)].x + _S1492.rows[int(0)].y + _S1492.rows[int(0)].z + _S1492.rows[int(1)].x + _S1492.rows[int(1)].y + _S1492.rows[int(1)].z + _S1492.rows[int(2)].x + _S1492.rows[int(2)].y + _S1492.rows[int(2)].z) / _S1234);
        H_2 = _S1493;
    }
    else
    {
        _S1234 = 0.0f;
        H_2 = _S1487.differential_0;
    }
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = _S1231.rows[int(2)].z;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1494, 0.0f);
    float _S1495 = _S1494.differential_0 + _S1234;
    float3  _S1496 = _S1175;
    *&((&_S1496)->z) = _S1495;
    Matrix<float, 3, 3>  _S1497 = _S1176;
    _S1497[int(2)] = _S1496;
    Matrix<float, 3, 3>  _S1498 = H_2 + _S1497;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1499;
    (&_S1499)->primal_0 = _S1230;
    (&_S1499)->differential_0 = _S1176;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1500;
    (&_S1500)->primal_0 = S_inv_0;
    (&_S1500)->differential_0 = _S1176;
    s_bwd_prop_mul_2(&_S1499, &_S1500, _S1498);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1501;
    (&_S1501)->primal_0 = T_1;
    (&_S1501)->differential_0 = _S1176;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1502;
    (&_S1502)->primal_0 = D_0;
    (&_S1502)->differential_0 = _S1176;
    s_bwd_prop_mul_2(&_S1501, &_S1502, _S1499.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1503 = _S1501;
    float3  _S1504 = make_float3 (_S1502.differential_0.rows[int(0)].x, _S1502.differential_0.rows[int(1)].y, _S1502.differential_0.rows[int(2)].z);
    float3  _S1505;
    if(_S1225)
    {
        if(_S1227)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1506;
            (&_S1506)->primal_0 = r1_1;
            (&_S1506)->differential_0 = _S1175;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1507;
            (&_S1507)->primal_0 = r2_11;
            (&_S1507)->differential_0 = _S1175;
            s_bwd_prop_cross_0(&_S1506, &_S1507, _S1504);
            _S1213 = _S1175;
            lambda_v_3 = _S1507.differential_0;
            _S1505 = _S1506.differential_0;
        }
        else
        {
            _S1213 = _S1504;
            lambda_v_3 = _S1175;
            _S1505 = _S1175;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1508;
        (&_S1508)->primal_0 = _S1226;
        (&_S1508)->differential_0 = _S1175;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1509;
        (&_S1509)->primal_0 = _S1226;
        (&_S1509)->differential_0 = _S1175;
        s_bwd_prop_dot_0(&_S1508, &_S1509, 0.0f);
        float3  _S1510 = _S1509.differential_0 + _S1508.differential_0 + _S1213;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1511;
        (&_S1511)->primal_0 = r0_1;
        (&_S1511)->differential_0 = _S1175;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1512;
        (&_S1512)->primal_0 = r2_11;
        (&_S1512)->differential_0 = _S1175;
        s_bwd_prop_cross_0(&_S1511, &_S1512, _S1510);
        float3  _S1513 = _S1512.differential_0 + lambda_v_3;
        _S1213 = _S1175;
        lambda_v_3 = _S1513;
        _S1226 = _S1505;
        _S1505 = _S1511.differential_0;
    }
    else
    {
        _S1213 = _S1504;
        lambda_v_3 = _S1175;
        _S1226 = _S1175;
        _S1505 = _S1175;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1514;
    (&_S1514)->primal_0 = _S1224;
    (&_S1514)->differential_0 = _S1175;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1515;
    (&_S1515)->primal_0 = _S1224;
    (&_S1515)->differential_0 = _S1175;
    s_bwd_prop_dot_0(&_S1514, &_S1515, 0.0f);
    float3  _S1516 = _S1515.differential_0 + _S1514.differential_0 + _S1213;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1517;
    (&_S1517)->primal_0 = r0_1;
    (&_S1517)->differential_0 = _S1175;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1518;
    (&_S1518)->primal_0 = r1_1;
    (&_S1518)->differential_0 = _S1175;
    s_bwd_prop_cross_0(&_S1517, &_S1518, _S1516);
    float3  _S1519 = _S1175;
    *&((&_S1519)->z) = lambda_v_3.z;
    *&((&_S1519)->y) = lambda_v_3.y;
    *&((&_S1519)->x) = lambda_v_3.x;
    float3  _S1520 = _S1518.differential_0 + _S1226;
    float3  _S1521 = _S1175;
    *&((&_S1521)->z) = _S1520.z;
    *&((&_S1521)->y) = _S1520.y;
    *&((&_S1521)->x) = _S1520.x;
    float3  _S1522 = _S1517.differential_0 + _S1505;
    float3  _S1523 = _S1175;
    *&((&_S1523)->z) = _S1522.z;
    *&((&_S1523)->y) = _S1522.y;
    *&((&_S1523)->x) = _S1522.x;
    Matrix<float, 3, 3>  _S1524 = _S1176;
    _S1524[int(2)] = _S1519;
    _S1524[int(1)] = _S1521;
    _S1524[int(0)] = _S1523;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1525;
    (&_S1525)->primal_0 = skew_0;
    (&_S1525)->differential_0 = _S1176;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1526;
    (&_S1526)->primal_0 = T_1;
    (&_S1526)->differential_0 = _S1176;
    s_bwd_prop_mul_2(&_S1525, &_S1526, _S1524);
    Matrix<float, 3, 3>  _S1527 = _S1526.differential_0 + _S1503.differential_0;
    float2  _S1528 = make_float2 (_S1525.differential_0.rows[int(2)].y + - _S1525.differential_0.rows[int(1)].z, _S1525.differential_0.rows[int(0)].z + - _S1525.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1529 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1530;
    (&_S1530)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1530)->differential_0 = _S1529;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1531;
    (&_S1531)->primal_0 = _S1216.color_params_0.n_0;
    (&_S1531)->differential_0 = _S1179;
    s_bwd_prop_mul_3(&_S1530, &_S1531, _S1528);
    float2  _S1532 = make_float2 (_S1527.rows[int(0)].z, _S1527.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1533;
    (&_S1533)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1533)->differential_0 = _S1529;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1534;
    (&_S1534)->primal_0 = _S1216.color_params_0.g_0;
    (&_S1534)->differential_0 = _S1179;
    s_bwd_prop_mul_3(&_S1533, &_S1534, _S1532);
    float2  _S1535 = make_float2 (_S1527.rows[int(0)].y, _S1527.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1536;
    (&_S1536)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1536)->differential_0 = _S1529;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1537;
    (&_S1537)->primal_0 = _S1216.color_params_0.r_0;
    (&_S1537)->differential_0 = _S1179;
    s_bwd_prop_mul_3(&_S1536, &_S1537, _S1535);
    float2  _S1538 = make_float2 (_S1527.rows[int(0)].x, _S1527.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1539;
    (&_S1539)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1539)->differential_0 = _S1529;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1540;
    (&_S1540)->primal_0 = _S1216.color_params_0.b_0;
    (&_S1540)->differential_0 = _S1179;
    s_bwd_prop_mul_3(&_S1539, &_S1540, _S1538);
    ColorPPISPParams_0 _S1541 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1541)->n_0 = _S1531.differential_0;
    (&_S1541)->g_0 = _S1534.differential_0;
    (&_S1541)->r_0 = _S1537.differential_0;
    (&_S1541)->b_0 = _S1540.differential_0;
    _S1213 = _S1491;
    *&((&_S1213)->z) = 0.0f;
    float _S1542 = rgb_out_2.z * _S1488;
    float _S1543 = _S1215 * _S1488;
    DiffPair_float_0 _S1544;
    (&_S1544)->primal_0 = falloff_2;
    (&_S1544)->differential_0 = 0.0f;
    DiffPair_float_0 _S1545;
    (&_S1545)->primal_0 = 0.0f;
    (&_S1545)->differential_0 = 0.0f;
    DiffPair_float_0 _S1546;
    (&_S1546)->primal_0 = 1.0f;
    (&_S1546)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1544, &_S1545, &_S1546, _S1542);
    float _S1547 = r2_10 * _S1544.differential_0;
    float _S1548 = r4_5 * _S1544.differential_0;
    float s_diff_r6_T_0 = _S1192 * _S1544.differential_0;
    float _S1549 = r6_2 * _S1544.differential_0;
    float _S1550 = r2_10 * (_S1191 * _S1544.differential_0 + r2_10 * s_diff_r6_T_0);
    float _S1551 = _S1190 * _S1544.differential_0 + r4_5 * s_diff_r6_T_0 + _S1550 + _S1550;
    float _S1552 = dy_5 * _S1551;
    float _S1553 = dx_7 * _S1551;
    float _S1554 = - (_S1552 + _S1552);
    float _S1555 = - (_S1553 + _S1553);
    *&((&_S1213)->y) = 0.0f;
    float _S1556 = rgb_out_2.y * _S1489;
    float _S1557 = _S1214 * _S1489;
    DiffPair_float_0 _S1558;
    (&_S1558)->primal_0 = falloff_1;
    (&_S1558)->differential_0 = 0.0f;
    DiffPair_float_0 _S1559;
    (&_S1559)->primal_0 = 0.0f;
    (&_S1559)->differential_0 = 0.0f;
    DiffPair_float_0 _S1560;
    (&_S1560)->primal_0 = 1.0f;
    (&_S1560)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1558, &_S1559, &_S1560, _S1556);
    float _S1561 = r2_9 * _S1558.differential_0;
    float _S1562 = r4_4 * _S1558.differential_0;
    float s_diff_r6_T_1 = _S1189 * _S1558.differential_0;
    float _S1563 = r6_1 * _S1558.differential_0;
    float _S1564 = r2_9 * (_S1188 * _S1558.differential_0 + r2_9 * s_diff_r6_T_1);
    float _S1565 = _S1187 * _S1558.differential_0 + r4_4 * s_diff_r6_T_1 + _S1564 + _S1564;
    float _S1566 = dy_4 * _S1565;
    float _S1567 = dx_6 * _S1565;
    float _S1568 = - (_S1566 + _S1566);
    float _S1569 = - (_S1567 + _S1567);
    *&((&_S1213)->x) = 0.0f;
    float _S1570 = rgb_out_2.x * _S1490;
    float _S1571 = _S1211 * _S1490;
    DiffPair_float_0 _S1572;
    (&_S1572)->primal_0 = falloff_0;
    (&_S1572)->differential_0 = 0.0f;
    DiffPair_float_0 _S1573;
    (&_S1573)->primal_0 = 0.0f;
    (&_S1573)->differential_0 = 0.0f;
    DiffPair_float_0 _S1574;
    (&_S1574)->primal_0 = 1.0f;
    (&_S1574)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1572, &_S1573, &_S1574, _S1570);
    float _S1575 = r2_8 * _S1572.differential_0;
    float _S1576 = r4_3 * _S1572.differential_0;
    float s_diff_r6_T_2 = _S1186 * _S1572.differential_0;
    float _S1577 = r6_0 * _S1572.differential_0;
    float _S1578 = r2_8 * (_S1185 * _S1572.differential_0 + r2_8 * s_diff_r6_T_2);
    float _S1579 = _S1184 * _S1572.differential_0 + r4_3 * s_diff_r6_T_2 + _S1578 + _S1578;
    float _S1580 = dy_3 * _S1579;
    float _S1581 = dx_5 * _S1579;
    float _S1582 = - (_S1580 + _S1580);
    float _S1583 = - (_S1581 + _S1581);
    float3  _S1584 = _S1175;
    *&((&_S1584)->z) = _S1543;
    *&((&_S1584)->y) = _S1557;
    *&((&_S1584)->x) = _S1571;
    float3  _S1585 = _S1213 + _S1584;
    float3  _S1586 = _S1174.primal_0 * _S1585;
    float3  _S1587 = _S1207 * _S1585;
    float _S1588 = _S1586.x + _S1586.y + _S1586.z;
    DiffPair_float_0 _S1589;
    (&_S1589)->primal_0 = _S1205.exposure_0;
    (&_S1589)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1589, _S1588);
    PPISPParams_0 _S1590 = PPISPParams_x24_syn_dzero_0();
    (&_S1590)->color_params_0 = _S1541;
    (&_S1590)->exposure_0 = _S1589.differential_0;
    _S1183 = _S1590;
    (&(&_S1183)->crf_params_0[int(2)])->center_0 = 0.0f;
    float _S1591 = _S1590.crf_params_0[int(2)].center_0 + _S1376;
    (&(&_S1183)->crf_params_0[int(2)])->gamma_0 = 0.0f;
    float _S1592 = _S1590.crf_params_0[int(2)].gamma_0 + _S1379.differential_0;
    (&(&_S1183)->crf_params_0[int(2)])->shoulder_0 = 0.0f;
    float _S1593 = _S1590.crf_params_0[int(2)].shoulder_0 + _S1383.differential_0;
    (&(&_S1183)->crf_params_0[int(2)])->toe_0 = 0.0f;
    float _S1594 = _S1590.crf_params_0[int(2)].toe_0 + _S1387.differential_0;
    (&(&_S1183)->crf_params_0[int(1)])->center_0 = 0.0f;
    float _S1595 = _S1590.crf_params_0[int(1)].center_0 + _S1420;
    (&(&_S1183)->crf_params_0[int(1)])->gamma_0 = 0.0f;
    float _S1596 = _S1590.crf_params_0[int(1)].gamma_0 + _S1423.differential_0;
    (&(&_S1183)->crf_params_0[int(1)])->shoulder_0 = 0.0f;
    float _S1597 = _S1590.crf_params_0[int(1)].shoulder_0 + _S1427.differential_0;
    (&(&_S1183)->crf_params_0[int(1)])->toe_0 = 0.0f;
    float _S1598 = _S1590.crf_params_0[int(1)].toe_0 + _S1431.differential_0;
    (&(&_S1183)->crf_params_0[int(0)])->center_0 = 0.0f;
    float _S1599 = _S1590.crf_params_0[int(0)].center_0 + _S1464;
    (&(&_S1183)->crf_params_0[int(0)])->gamma_0 = 0.0f;
    float _S1600 = _S1590.crf_params_0[int(0)].gamma_0 + _S1467.differential_0;
    (&(&_S1183)->crf_params_0[int(0)])->shoulder_0 = 0.0f;
    float _S1601 = _S1590.crf_params_0[int(0)].shoulder_0 + _S1471.differential_0;
    (&(&_S1183)->crf_params_0[int(0)])->toe_0 = 0.0f;
    float _S1602 = _S1590.crf_params_0[int(0)].toe_0 + _S1475.differential_0;
    *&((&(&(&_S1183)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1183)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1183)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1603 = _S1549 + _S1590.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1183)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1604 = _S1548 + _S1590.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1183)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1605 = _S1547 + _S1590.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1183)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1606 = _S1554 + _S1590.vignette_params_0[int(2)].cy_0;
    (&(&_S1183)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1607 = _S1555 + _S1590.vignette_params_0[int(2)].cx_0;
    (&(&_S1183)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1608 = _S1563 + _S1590.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1183)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1609 = _S1562 + _S1590.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1183)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1610 = _S1561 + _S1590.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1183)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1611 = _S1568 + _S1590.vignette_params_0[int(1)].cy_0;
    (&(&_S1183)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1612 = _S1569 + _S1590.vignette_params_0[int(1)].cx_0;
    (&(&_S1183)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1613 = _S1577 + _S1590.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1183)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1614 = _S1576 + _S1590.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1183)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1615 = _S1575 + _S1590.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1183)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1616 = _S1582 + _S1590.vignette_params_0[int(0)].cy_0;
    (&(&_S1183)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1617 = _S1583 + _S1590.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 36>  _S1618;
    _S1618[int(0)] = 0.0f;
    _S1618[int(1)] = 0.0f;
    _S1618[int(2)] = 0.0f;
    _S1618[int(3)] = 0.0f;
    _S1618[int(4)] = 0.0f;
    _S1618[int(5)] = 0.0f;
    _S1618[int(6)] = 0.0f;
    _S1618[int(7)] = 0.0f;
    _S1618[int(8)] = 0.0f;
    _S1618[int(9)] = 0.0f;
    _S1618[int(10)] = 0.0f;
    _S1618[int(11)] = 0.0f;
    _S1618[int(12)] = 0.0f;
    _S1618[int(13)] = 0.0f;
    _S1618[int(14)] = 0.0f;
    _S1618[int(15)] = 0.0f;
    _S1618[int(16)] = 0.0f;
    _S1618[int(17)] = 0.0f;
    _S1618[int(18)] = 0.0f;
    _S1618[int(19)] = 0.0f;
    _S1618[int(20)] = 0.0f;
    _S1618[int(21)] = 0.0f;
    _S1618[int(22)] = 0.0f;
    _S1618[int(23)] = 0.0f;
    _S1618[int(24)] = 0.0f;
    _S1618[int(25)] = 0.0f;
    _S1618[int(26)] = 0.0f;
    _S1618[int(27)] = 0.0f;
    _S1618[int(28)] = 0.0f;
    _S1618[int(29)] = 0.0f;
    _S1618[int(30)] = 0.0f;
    _S1618[int(31)] = 0.0f;
    _S1618[int(32)] = 0.0f;
    _S1618[int(33)] = 0.0f;
    _S1618[int(34)] = 0.0f;
    _S1618[int(35)] = 0.0f;
    _S1618[int(8)] = _S1610;
    _S1618[int(16)] = _S1590.color_params_0.b_0.x;
    _S1618[int(15)] = _S1603;
    _S1618[int(14)] = _S1604;
    _S1618[int(13)] = _S1605;
    _S1618[int(12)] = _S1606;
    _S1618[int(11)] = _S1607;
    _S1618[int(10)] = _S1608;
    _S1618[int(9)] = _S1609;
    _S1618[int(17)] = _S1590.color_params_0.b_0.y;
    _S1618[int(7)] = _S1611;
    _S1618[int(6)] = _S1612;
    _S1618[int(5)] = _S1613;
    _S1618[int(4)] = _S1614;
    _S1618[int(3)] = _S1615;
    _S1618[int(2)] = _S1616;
    _S1618[int(1)] = _S1617;
    _S1618[int(0)] = _S1183.exposure_0;
    _S1618[int(26)] = _S1600;
    _S1618[int(34)] = _S1592;
    _S1618[int(33)] = _S1593;
    _S1618[int(32)] = _S1594;
    _S1618[int(31)] = _S1595;
    _S1618[int(30)] = _S1596;
    _S1618[int(29)] = _S1597;
    _S1618[int(28)] = _S1598;
    _S1618[int(27)] = _S1599;
    _S1618[int(35)] = _S1591;
    _S1618[int(25)] = _S1601;
    _S1618[int(24)] = _S1602;
    _S1618[int(23)] = _S1590.color_params_0.n_0.y;
    _S1618[int(22)] = _S1590.color_params_0.n_0.x;
    _S1618[int(21)] = _S1590.color_params_0.g_0.y;
    _S1618[int(20)] = _S1590.color_params_0.g_0.x;
    _S1618[int(19)] = _S1590.color_params_0.r_0.y;
    _S1618[int(18)] = _S1590.color_params_0.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1618;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1587;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1619, float2  _S1620, float2  _S1621, float2  _S1622, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1623, float3  _S1624)
{
    s_bwd_prop_apply_ppisp_0(_S1619, _S1620, _S1621, _S1622, _S1623, _S1624);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_1, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, FixedArray<float, 36>  * params_1, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1625 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_1;
    (&dp_rgb_in_0)->differential_0 = _S1625;
    FixedArray<float, 36>  _S1626 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = *params_1;
    (&dp_params_0)->differential_0 = _S1626;
    s_bwd_apply_ppisp_0(&dp_rgb_in_0, pix_coord_2, image_center_2, img_size_2, &dp_params_0, grad_out_0);
    *grad_rgb_in_0 = dp_rgb_in_0.differential_0;
    *grad_params_0 = (&dp_params_0)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  * params_2, FixedArray<float, 22>  * _S1627)
{
    PPISPParams_0 p_1;
    (&p_1)->exposure_0 = (*params_2)[int(0)];
    (&(&p_1)->vignette_params_0[int(0)])->cx_0 = (*params_2)[int(1)];
    (&(&p_1)->vignette_params_0[int(0)])->cy_0 = (*params_2)[int(2)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha0_0 = (*params_2)[int(3)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha1_0 = (*params_2)[int(4)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha2_0 = (*params_2)[int(5)];
    (&(&p_1)->vignette_params_0[int(1)])->cx_0 = (*params_2)[int(6)];
    (&(&p_1)->vignette_params_0[int(1)])->cy_0 = (*params_2)[int(7)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha0_0 = (*params_2)[int(8)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha1_0 = (*params_2)[int(9)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha2_0 = (*params_2)[int(10)];
    (&(&p_1)->vignette_params_0[int(2)])->cx_0 = (*params_2)[int(11)];
    (&(&p_1)->vignette_params_0[int(2)])->cy_0 = (*params_2)[int(12)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha0_0 = (*params_2)[int(13)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha1_0 = (*params_2)[int(14)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha2_0 = (*params_2)[int(15)];
    *&((&(&(&p_1)->color_params_0)->b_0)->x) = (*params_2)[int(16)];
    *&((&(&(&p_1)->color_params_0)->b_0)->y) = (*params_2)[int(17)];
    *&((&(&(&p_1)->color_params_0)->r_0)->x) = (*params_2)[int(18)];
    *&((&(&(&p_1)->color_params_0)->r_0)->y) = (*params_2)[int(19)];
    *&((&(&(&p_1)->color_params_0)->g_0)->x) = (*params_2)[int(20)];
    *&((&(&(&p_1)->color_params_0)->g_0)->y) = (*params_2)[int(21)];
    *&((&(&(&p_1)->color_params_0)->n_0)->x) = (*params_2)[int(22)];
    *&((&(&(&p_1)->color_params_0)->n_0)->y) = (*params_2)[int(23)];
    (&(&p_1)->crf_params_0[int(0)])->toe_0 = (*params_2)[int(24)];
    (&(&p_1)->crf_params_0[int(0)])->shoulder_0 = (*params_2)[int(25)];
    (&(&p_1)->crf_params_0[int(0)])->gamma_0 = (*params_2)[int(26)];
    (&(&p_1)->crf_params_0[int(0)])->center_0 = (*params_2)[int(27)];
    (&(&p_1)->crf_params_0[int(1)])->toe_0 = (*params_2)[int(28)];
    (&(&p_1)->crf_params_0[int(1)])->shoulder_0 = (*params_2)[int(29)];
    (&(&p_1)->crf_params_0[int(1)])->gamma_0 = (*params_2)[int(30)];
    (&(&p_1)->crf_params_0[int(1)])->center_0 = (*params_2)[int(31)];
    (&(&p_1)->crf_params_0[int(2)])->toe_0 = (*params_2)[int(32)];
    (&(&p_1)->crf_params_0[int(2)])->shoulder_0 = (*params_2)[int(33)];
    (&(&p_1)->crf_params_0[int(2)])->gamma_0 = (*params_2)[int(34)];
    (&(&p_1)->crf_params_0[int(2)])->center_0 = (*params_2)[int(35)];
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
    losses_3[int(0)] = p_1.exposure_0;
    float _S1628 = p_1.vignette_params_0[int(0)].cx_0;
    float _S1629 = p_1.vignette_params_0[int(0)].cy_0;
    float _S1630 = p_1.vignette_params_0[int(1)].cx_0;
    float _S1631 = p_1.vignette_params_0[int(1)].cy_0;
    float _S1632 = p_1.vignette_params_0[int(2)].cx_0;
    float _S1633 = p_1.vignette_params_0[int(2)].cy_0;
    losses_3[int(1)] = _S1628 * _S1628 + _S1629 * _S1629 + _S1630 * _S1630 + _S1631 * _S1631 + _S1632 * _S1632 + _S1633 * _S1633;
    losses_3[int(2)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha0_0)));
    losses_3[int(3)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha1_0)));
    losses_3[int(4)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha2_0)));
    float mean_3 = (p_1.vignette_params_0[int(0)].cx_0 + p_1.vignette_params_0[int(1)].cx_0 + p_1.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S1634 = p_1.vignette_params_0[int(0)].cx_0 - mean_3;
    float _S1635 = p_1.vignette_params_0[int(1)].cx_0 - mean_3;
    float _S1636 = p_1.vignette_params_0[int(2)].cx_0 - mean_3;
    losses_3[int(5)] = (_S1634 * _S1634 + _S1635 * _S1635 + _S1636 * _S1636) / 3.0f;
    float mean_4 = (p_1.vignette_params_0[int(0)].cy_0 + p_1.vignette_params_0[int(1)].cy_0 + p_1.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S1637 = p_1.vignette_params_0[int(0)].cy_0 - mean_4;
    float _S1638 = p_1.vignette_params_0[int(1)].cy_0 - mean_4;
    float _S1639 = p_1.vignette_params_0[int(2)].cy_0 - mean_4;
    losses_3[int(6)] = (_S1637 * _S1637 + _S1638 * _S1638 + _S1639 * _S1639) / 3.0f;
    float mean_5 = (p_1.vignette_params_0[int(0)].alpha0_0 + p_1.vignette_params_0[int(1)].alpha0_0 + p_1.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S1640 = p_1.vignette_params_0[int(0)].alpha0_0 - mean_5;
    float _S1641 = p_1.vignette_params_0[int(1)].alpha0_0 - mean_5;
    float _S1642 = p_1.vignette_params_0[int(2)].alpha0_0 - mean_5;
    losses_3[int(7)] = (_S1640 * _S1640 + _S1641 * _S1641 + _S1642 * _S1642) / 3.0f;
    float mean_6 = (p_1.vignette_params_0[int(0)].alpha1_0 + p_1.vignette_params_0[int(1)].alpha1_0 + p_1.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S1643 = p_1.vignette_params_0[int(0)].alpha1_0 - mean_6;
    float _S1644 = p_1.vignette_params_0[int(1)].alpha1_0 - mean_6;
    float _S1645 = p_1.vignette_params_0[int(2)].alpha1_0 - mean_6;
    losses_3[int(8)] = (_S1643 * _S1643 + _S1644 * _S1644 + _S1645 * _S1645) / 3.0f;
    float mean_7 = (p_1.vignette_params_0[int(0)].alpha2_0 + p_1.vignette_params_0[int(1)].alpha2_0 + p_1.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S1646 = p_1.vignette_params_0[int(0)].alpha2_0 - mean_7;
    float _S1647 = p_1.vignette_params_0[int(1)].alpha2_0 - mean_7;
    float _S1648 = p_1.vignette_params_0[int(2)].alpha2_0 - mean_7;
    losses_3[int(9)] = (_S1646 * _S1646 + _S1647 * _S1647 + _S1648 * _S1648) / 3.0f;
    float2  bd_1 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    losses_3[int(10)] = bd_1.x;
    losses_3[int(11)] = bd_1.y;
    losses_3[int(12)] = rd_1.x;
    losses_3[int(13)] = rd_1.y;
    losses_3[int(14)] = gd_1.x;
    losses_3[int(15)] = gd_1.y;
    losses_3[int(16)] = nd_1.x;
    losses_3[int(17)] = nd_1.y;
    float mean_8 = (p_1.crf_params_0[int(0)].toe_0 + p_1.crf_params_0[int(1)].toe_0 + p_1.crf_params_0[int(2)].toe_0) / 3.0f;
    float _S1649 = p_1.crf_params_0[int(0)].toe_0 - mean_8;
    float _S1650 = p_1.crf_params_0[int(1)].toe_0 - mean_8;
    float _S1651 = p_1.crf_params_0[int(2)].toe_0 - mean_8;
    losses_3[int(18)] = (_S1649 * _S1649 + _S1650 * _S1650 + _S1651 * _S1651) / 3.0f;
    float mean_9 = (p_1.crf_params_0[int(0)].shoulder_0 + p_1.crf_params_0[int(1)].shoulder_0 + p_1.crf_params_0[int(2)].shoulder_0) / 3.0f;
    float _S1652 = p_1.crf_params_0[int(0)].shoulder_0 - mean_9;
    float _S1653 = p_1.crf_params_0[int(1)].shoulder_0 - mean_9;
    float _S1654 = p_1.crf_params_0[int(2)].shoulder_0 - mean_9;
    losses_3[int(19)] = (_S1652 * _S1652 + _S1653 * _S1653 + _S1654 * _S1654) / 3.0f;
    float mean_10 = (p_1.crf_params_0[int(0)].gamma_0 + p_1.crf_params_0[int(1)].gamma_0 + p_1.crf_params_0[int(2)].gamma_0) / 3.0f;
    float _S1655 = p_1.crf_params_0[int(0)].gamma_0 - mean_10;
    float _S1656 = p_1.crf_params_0[int(1)].gamma_0 - mean_10;
    float _S1657 = p_1.crf_params_0[int(2)].gamma_0 - mean_10;
    losses_3[int(20)] = (_S1655 * _S1655 + _S1656 * _S1656 + _S1657 * _S1657) / 3.0f;
    float mean_11 = (p_1.crf_params_0[int(0)].center_0 + p_1.crf_params_0[int(1)].center_0 + p_1.crf_params_0[int(2)].center_0) / 3.0f;
    float _S1658 = p_1.crf_params_0[int(0)].center_0 - mean_11;
    float _S1659 = p_1.crf_params_0[int(1)].center_0 - mean_11;
    float _S1660 = p_1.crf_params_0[int(2)].center_0 - mean_11;
    losses_3[int(21)] = (_S1658 * _S1658 + _S1659 * _S1659 + _S1660 * _S1660) / 3.0f;
    *_S1627 = losses_3;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_1, FixedArray<float, 22>  * _s_dOut_12)
{
    VignettingChannelParams_0 _S1661 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1662 = {
        _S1661, _S1661, _S1661
    };
    float2  _S1663 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1664 = { _S1663, _S1663, _S1663, _S1663 };
    CRFPPISPChannelParams_0 _S1665 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1666 = {
        _S1665, _S1665, _S1665
    };
    PPISPParams_0 _S1667;
    (&_S1667)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1667)->vignette_params_0 = _S1662;
    (&_S1667)->color_params_0 = _S1664;
    (&_S1667)->crf_params_0 = _S1666;
    (&(&_S1667)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1667)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    (&(&_S1667)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    (&(&_S1667)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    (&(&_S1667)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1667)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1667)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    (&(&_S1667)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    (&(&_S1667)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    (&(&_S1667)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1667)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1667)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    (&(&_S1667)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    (&(&_S1667)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    (&(&_S1667)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1667)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1667)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1667)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1667)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1667)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1667)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1667)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1667)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    (&(&_S1667)->crf_params_0[int(0)])->toe_0 = dpparams_1->primal_0[int(24)];
    (&(&_S1667)->crf_params_0[int(0)])->shoulder_0 = dpparams_1->primal_0[int(25)];
    (&(&_S1667)->crf_params_0[int(0)])->gamma_0 = dpparams_1->primal_0[int(26)];
    (&(&_S1667)->crf_params_0[int(0)])->center_0 = dpparams_1->primal_0[int(27)];
    (&(&_S1667)->crf_params_0[int(1)])->toe_0 = dpparams_1->primal_0[int(28)];
    (&(&_S1667)->crf_params_0[int(1)])->shoulder_0 = dpparams_1->primal_0[int(29)];
    (&(&_S1667)->crf_params_0[int(1)])->gamma_0 = dpparams_1->primal_0[int(30)];
    (&(&_S1667)->crf_params_0[int(1)])->center_0 = dpparams_1->primal_0[int(31)];
    (&(&_S1667)->crf_params_0[int(2)])->toe_0 = dpparams_1->primal_0[int(32)];
    (&(&_S1667)->crf_params_0[int(2)])->shoulder_0 = dpparams_1->primal_0[int(33)];
    (&(&_S1667)->crf_params_0[int(2)])->gamma_0 = dpparams_1->primal_0[int(34)];
    (&(&_S1667)->crf_params_0[int(2)])->center_0 = dpparams_1->primal_0[int(35)];
    float mean_12 = (dpparams_1->primal_0[int(1)] + dpparams_1->primal_0[int(6)] + dpparams_1->primal_0[int(11)]) / 3.0f;
    float _S1668 = dpparams_1->primal_0[int(1)] - mean_12;
    float _S1669 = dpparams_1->primal_0[int(6)] - mean_12;
    float _S1670 = dpparams_1->primal_0[int(11)] - mean_12;
    float mean_13 = (dpparams_1->primal_0[int(2)] + dpparams_1->primal_0[int(7)] + dpparams_1->primal_0[int(12)]) / 3.0f;
    float _S1671 = dpparams_1->primal_0[int(2)] - mean_13;
    float _S1672 = dpparams_1->primal_0[int(7)] - mean_13;
    float _S1673 = dpparams_1->primal_0[int(12)] - mean_13;
    float mean_14 = (dpparams_1->primal_0[int(3)] + dpparams_1->primal_0[int(8)] + dpparams_1->primal_0[int(13)]) / 3.0f;
    float _S1674 = dpparams_1->primal_0[int(3)] - mean_14;
    float _S1675 = dpparams_1->primal_0[int(8)] - mean_14;
    float _S1676 = dpparams_1->primal_0[int(13)] - mean_14;
    float mean_15 = (dpparams_1->primal_0[int(4)] + dpparams_1->primal_0[int(9)] + dpparams_1->primal_0[int(14)]) / 3.0f;
    float _S1677 = dpparams_1->primal_0[int(4)] - mean_15;
    float _S1678 = dpparams_1->primal_0[int(9)] - mean_15;
    float _S1679 = dpparams_1->primal_0[int(14)] - mean_15;
    float mean_16 = (dpparams_1->primal_0[int(5)] + dpparams_1->primal_0[int(10)] + dpparams_1->primal_0[int(15)]) / 3.0f;
    float _S1680 = dpparams_1->primal_0[int(5)] - mean_16;
    float _S1681 = dpparams_1->primal_0[int(10)] - mean_16;
    float _S1682 = dpparams_1->primal_0[int(15)] - mean_16;
    float mean_17 = (dpparams_1->primal_0[int(24)] + dpparams_1->primal_0[int(28)] + dpparams_1->primal_0[int(32)]) / 3.0f;
    float mean_18 = (dpparams_1->primal_0[int(25)] + dpparams_1->primal_0[int(29)] + dpparams_1->primal_0[int(33)]) / 3.0f;
    float mean_19 = (dpparams_1->primal_0[int(26)] + dpparams_1->primal_0[int(30)] + dpparams_1->primal_0[int(34)]) / 3.0f;
    float mean_20 = (dpparams_1->primal_0[int(27)] + dpparams_1->primal_0[int(31)] + dpparams_1->primal_0[int(35)]) / 3.0f;
    float _S1683 = 0.3333333432674408f * (*_s_dOut_12)[int(21)];
    float _S1684 = (dpparams_1->primal_0[int(35)] - mean_20) * _S1683;
    float _S1685 = _S1684 + _S1684;
    float _S1686 = (dpparams_1->primal_0[int(31)] - mean_20) * _S1683;
    float _S1687 = _S1686 + _S1686;
    float _S1688 = (dpparams_1->primal_0[int(27)] - mean_20) * _S1683;
    float _S1689 = _S1688 + _S1688;
    float _S1690 = 0.3333333432674408f * (- _S1685 + - _S1687 + - _S1689);
    float _S1691 = 0.3333333432674408f * (*_s_dOut_12)[int(20)];
    float _S1692 = (dpparams_1->primal_0[int(34)] - mean_19) * _S1691;
    float _S1693 = _S1692 + _S1692;
    float _S1694 = (dpparams_1->primal_0[int(30)] - mean_19) * _S1691;
    float _S1695 = _S1694 + _S1694;
    float _S1696 = (dpparams_1->primal_0[int(26)] - mean_19) * _S1691;
    float _S1697 = _S1696 + _S1696;
    float _S1698 = 0.3333333432674408f * (- _S1693 + - _S1695 + - _S1697);
    float _S1699 = 0.3333333432674408f * (*_s_dOut_12)[int(19)];
    float _S1700 = (dpparams_1->primal_0[int(33)] - mean_18) * _S1699;
    float _S1701 = _S1700 + _S1700;
    float _S1702 = (dpparams_1->primal_0[int(29)] - mean_18) * _S1699;
    float _S1703 = _S1702 + _S1702;
    float _S1704 = (dpparams_1->primal_0[int(25)] - mean_18) * _S1699;
    float _S1705 = _S1704 + _S1704;
    float _S1706 = 0.3333333432674408f * (- _S1701 + - _S1703 + - _S1705);
    float _S1707 = 0.3333333432674408f * (*_s_dOut_12)[int(18)];
    float _S1708 = (dpparams_1->primal_0[int(32)] - mean_17) * _S1707;
    float _S1709 = _S1708 + _S1708;
    float _S1710 = (dpparams_1->primal_0[int(28)] - mean_17) * _S1707;
    float _S1711 = _S1710 + _S1710;
    float _S1712 = (dpparams_1->primal_0[int(24)] - mean_17) * _S1707;
    float _S1713 = _S1712 + _S1712;
    float _S1714 = 0.3333333432674408f * (- _S1709 + - _S1711 + - _S1713);
    float2  _S1715 = make_float2 ((*_s_dOut_12)[int(16)], (*_s_dOut_12)[int(17)]);
    Matrix<float, 2, 2>  _S1716 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1717;
    (&_S1717)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1717)->differential_0 = _S1716;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1718;
    (&_S1718)->primal_0 = _S1667.color_params_0.n_0;
    (&_S1718)->differential_0 = _S1663;
    s_bwd_prop_mul_3(&_S1717, &_S1718, _S1715);
    float2  _S1719 = make_float2 ((*_s_dOut_12)[int(14)], (*_s_dOut_12)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1720;
    (&_S1720)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1720)->differential_0 = _S1716;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1721;
    (&_S1721)->primal_0 = _S1667.color_params_0.g_0;
    (&_S1721)->differential_0 = _S1663;
    s_bwd_prop_mul_3(&_S1720, &_S1721, _S1719);
    float2  _S1722 = make_float2 ((*_s_dOut_12)[int(12)], (*_s_dOut_12)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1723;
    (&_S1723)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1723)->differential_0 = _S1716;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1724;
    (&_S1724)->primal_0 = _S1667.color_params_0.r_0;
    (&_S1724)->differential_0 = _S1663;
    s_bwd_prop_mul_3(&_S1723, &_S1724, _S1722);
    float2  _S1725 = make_float2 ((*_s_dOut_12)[int(10)], (*_s_dOut_12)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1726;
    (&_S1726)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1726)->differential_0 = _S1716;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1727;
    (&_S1727)->primal_0 = _S1667.color_params_0.b_0;
    (&_S1727)->differential_0 = _S1663;
    s_bwd_prop_mul_3(&_S1726, &_S1727, _S1725);
    ColorPPISPParams_0 _S1728 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1728)->n_0 = _S1718.differential_0;
    (&_S1728)->g_0 = _S1721.differential_0;
    (&_S1728)->r_0 = _S1724.differential_0;
    (&_S1728)->b_0 = _S1727.differential_0;
    float _S1729 = 0.3333333432674408f * (*_s_dOut_12)[int(9)];
    float _S1730 = _S1682 * _S1729;
    float _S1731 = _S1730 + _S1730;
    float _S1732 = _S1681 * _S1729;
    float _S1733 = _S1732 + _S1732;
    float _S1734 = _S1680 * _S1729;
    float _S1735 = _S1734 + _S1734;
    float _S1736 = 0.3333333432674408f * (- _S1731 + - _S1733 + - _S1735);
    float _S1737 = 0.3333333432674408f * (*_s_dOut_12)[int(8)];
    float _S1738 = _S1679 * _S1737;
    float _S1739 = _S1738 + _S1738;
    float _S1740 = _S1678 * _S1737;
    float _S1741 = _S1740 + _S1740;
    float _S1742 = _S1677 * _S1737;
    float _S1743 = _S1742 + _S1742;
    float _S1744 = 0.3333333432674408f * (- _S1739 + - _S1741 + - _S1743);
    float _S1745 = 0.3333333432674408f * (*_s_dOut_12)[int(7)];
    float _S1746 = _S1676 * _S1745;
    float _S1747 = _S1746 + _S1746;
    float _S1748 = _S1675 * _S1745;
    float _S1749 = _S1748 + _S1748;
    float _S1750 = _S1674 * _S1745;
    float _S1751 = _S1750 + _S1750;
    float _S1752 = 0.3333333432674408f * (- _S1747 + - _S1749 + - _S1751);
    float _S1753 = 0.3333333432674408f * (*_s_dOut_12)[int(6)];
    float _S1754 = _S1673 * _S1753;
    float _S1755 = _S1754 + _S1754;
    float _S1756 = _S1672 * _S1753;
    float _S1757 = _S1756 + _S1756;
    float _S1758 = _S1671 * _S1753;
    float _S1759 = _S1758 + _S1758;
    float _S1760 = 0.3333333432674408f * (- _S1755 + - _S1757 + - _S1759);
    float _S1761 = 0.3333333432674408f * (*_s_dOut_12)[int(5)];
    float _S1762 = _S1670 * _S1761;
    float _S1763 = _S1762 + _S1762;
    float _S1764 = _S1669 * _S1761;
    float _S1765 = _S1764 + _S1764;
    float _S1766 = _S1668 * _S1761;
    float _S1767 = _S1766 + _S1766;
    float _S1768 = 0.3333333432674408f * (- _S1763 + - _S1765 + - _S1767);
    DiffPair_float_0 _S1769;
    (&_S1769)->primal_0 = 0.0f;
    (&_S1769)->differential_0 = 0.0f;
    DiffPair_float_0 _S1770;
    (&_S1770)->primal_0 = dpparams_1->primal_0[int(15)];
    (&_S1770)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1769, &_S1770, (*_s_dOut_12)[int(4)]);
    DiffPair_float_0 _S1771;
    (&_S1771)->primal_0 = 0.0f;
    (&_S1771)->differential_0 = 0.0f;
    DiffPair_float_0 _S1772;
    (&_S1772)->primal_0 = dpparams_1->primal_0[int(10)];
    (&_S1772)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1771, &_S1772, (*_s_dOut_12)[int(4)]);
    DiffPair_float_0 _S1773;
    (&_S1773)->primal_0 = 0.0f;
    (&_S1773)->differential_0 = 0.0f;
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = dpparams_1->primal_0[int(5)];
    (&_S1774)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1773, &_S1774, (*_s_dOut_12)[int(4)]);
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = 0.0f;
    (&_S1775)->differential_0 = 0.0f;
    DiffPair_float_0 _S1776;
    (&_S1776)->primal_0 = dpparams_1->primal_0[int(14)];
    (&_S1776)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1775, &_S1776, (*_s_dOut_12)[int(3)]);
    DiffPair_float_0 _S1777;
    (&_S1777)->primal_0 = 0.0f;
    (&_S1777)->differential_0 = 0.0f;
    DiffPair_float_0 _S1778;
    (&_S1778)->primal_0 = dpparams_1->primal_0[int(9)];
    (&_S1778)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1777, &_S1778, (*_s_dOut_12)[int(3)]);
    DiffPair_float_0 _S1779;
    (&_S1779)->primal_0 = 0.0f;
    (&_S1779)->differential_0 = 0.0f;
    DiffPair_float_0 _S1780;
    (&_S1780)->primal_0 = dpparams_1->primal_0[int(4)];
    (&_S1780)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1779, &_S1780, (*_s_dOut_12)[int(3)]);
    DiffPair_float_0 _S1781;
    (&_S1781)->primal_0 = 0.0f;
    (&_S1781)->differential_0 = 0.0f;
    DiffPair_float_0 _S1782;
    (&_S1782)->primal_0 = dpparams_1->primal_0[int(13)];
    (&_S1782)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1781, &_S1782, (*_s_dOut_12)[int(2)]);
    DiffPair_float_0 _S1783;
    (&_S1783)->primal_0 = 0.0f;
    (&_S1783)->differential_0 = 0.0f;
    DiffPair_float_0 _S1784;
    (&_S1784)->primal_0 = dpparams_1->primal_0[int(8)];
    (&_S1784)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1783, &_S1784, (*_s_dOut_12)[int(2)]);
    DiffPair_float_0 _S1785;
    (&_S1785)->primal_0 = 0.0f;
    (&_S1785)->differential_0 = 0.0f;
    DiffPair_float_0 _S1786;
    (&_S1786)->primal_0 = dpparams_1->primal_0[int(3)];
    (&_S1786)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1785, &_S1786, (*_s_dOut_12)[int(2)]);
    float _S1787 = dpparams_1->primal_0[int(12)] * (*_s_dOut_12)[int(1)];
    float _S1788 = dpparams_1->primal_0[int(11)] * (*_s_dOut_12)[int(1)];
    float _S1789 = dpparams_1->primal_0[int(7)] * (*_s_dOut_12)[int(1)];
    float _S1790 = dpparams_1->primal_0[int(6)] * (*_s_dOut_12)[int(1)];
    float _S1791 = dpparams_1->primal_0[int(2)] * (*_s_dOut_12)[int(1)];
    float _S1792 = dpparams_1->primal_0[int(1)] * (*_s_dOut_12)[int(1)];
    PPISPParams_0 _S1793 = PPISPParams_x24_syn_dzero_0();
    (&_S1793)->color_params_0 = _S1728;
    (&_S1793)->exposure_0 = (*_s_dOut_12)[int(0)];
    _S1667 = _S1793;
    (&(&_S1667)->crf_params_0[int(2)])->center_0 = 0.0f;
    float _S1794 = _S1685 + _S1690 + _S1793.crf_params_0[int(2)].center_0;
    (&(&_S1667)->crf_params_0[int(2)])->gamma_0 = 0.0f;
    float _S1795 = _S1693 + _S1698 + _S1793.crf_params_0[int(2)].gamma_0;
    (&(&_S1667)->crf_params_0[int(2)])->shoulder_0 = 0.0f;
    float _S1796 = _S1701 + _S1706 + _S1793.crf_params_0[int(2)].shoulder_0;
    (&(&_S1667)->crf_params_0[int(2)])->toe_0 = 0.0f;
    float _S1797 = _S1709 + _S1714 + _S1793.crf_params_0[int(2)].toe_0;
    (&(&_S1667)->crf_params_0[int(1)])->center_0 = 0.0f;
    float _S1798 = _S1687 + _S1690 + _S1793.crf_params_0[int(1)].center_0;
    (&(&_S1667)->crf_params_0[int(1)])->gamma_0 = 0.0f;
    float _S1799 = _S1695 + _S1698 + _S1793.crf_params_0[int(1)].gamma_0;
    (&(&_S1667)->crf_params_0[int(1)])->shoulder_0 = 0.0f;
    float _S1800 = _S1703 + _S1706 + _S1793.crf_params_0[int(1)].shoulder_0;
    (&(&_S1667)->crf_params_0[int(1)])->toe_0 = 0.0f;
    float _S1801 = _S1711 + _S1714 + _S1793.crf_params_0[int(1)].toe_0;
    (&(&_S1667)->crf_params_0[int(0)])->center_0 = 0.0f;
    float _S1802 = _S1689 + _S1690 + _S1793.crf_params_0[int(0)].center_0;
    (&(&_S1667)->crf_params_0[int(0)])->gamma_0 = 0.0f;
    float _S1803 = _S1697 + _S1698 + _S1793.crf_params_0[int(0)].gamma_0;
    (&(&_S1667)->crf_params_0[int(0)])->shoulder_0 = 0.0f;
    float _S1804 = _S1705 + _S1706 + _S1793.crf_params_0[int(0)].shoulder_0;
    (&(&_S1667)->crf_params_0[int(0)])->toe_0 = 0.0f;
    float _S1805 = _S1713 + _S1714 + _S1793.crf_params_0[int(0)].toe_0;
    *&((&(&(&_S1667)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1667)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1667)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1806 = _S1731 + _S1736 + _S1770.differential_0 + _S1793.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1667)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1807 = _S1739 + _S1744 + _S1776.differential_0 + _S1793.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1667)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1808 = _S1747 + _S1752 + _S1782.differential_0 + _S1793.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1667)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1809 = _S1755 + _S1760 + _S1787 + _S1787 + _S1793.vignette_params_0[int(2)].cy_0;
    (&(&_S1667)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1810 = _S1763 + _S1768 + _S1788 + _S1788 + _S1793.vignette_params_0[int(2)].cx_0;
    (&(&_S1667)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1811 = _S1733 + _S1736 + _S1772.differential_0 + _S1793.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1667)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1812 = _S1741 + _S1744 + _S1778.differential_0 + _S1793.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1667)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1813 = _S1749 + _S1752 + _S1784.differential_0 + _S1793.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1667)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1814 = _S1757 + _S1760 + _S1789 + _S1789 + _S1793.vignette_params_0[int(1)].cy_0;
    (&(&_S1667)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1815 = _S1765 + _S1768 + _S1790 + _S1790 + _S1793.vignette_params_0[int(1)].cx_0;
    (&(&_S1667)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1816 = _S1735 + _S1736 + _S1774.differential_0 + _S1793.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1667)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1817 = _S1743 + _S1744 + _S1780.differential_0 + _S1793.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1667)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1818 = _S1751 + _S1752 + _S1786.differential_0 + _S1793.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1667)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1819 = _S1759 + _S1760 + _S1791 + _S1791 + _S1793.vignette_params_0[int(0)].cy_0;
    (&(&_S1667)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1820 = _S1767 + _S1768 + _S1792 + _S1792 + _S1793.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 36>  _S1821;
    _S1821[int(0)] = 0.0f;
    _S1821[int(1)] = 0.0f;
    _S1821[int(2)] = 0.0f;
    _S1821[int(3)] = 0.0f;
    _S1821[int(4)] = 0.0f;
    _S1821[int(5)] = 0.0f;
    _S1821[int(6)] = 0.0f;
    _S1821[int(7)] = 0.0f;
    _S1821[int(8)] = 0.0f;
    _S1821[int(9)] = 0.0f;
    _S1821[int(10)] = 0.0f;
    _S1821[int(11)] = 0.0f;
    _S1821[int(12)] = 0.0f;
    _S1821[int(13)] = 0.0f;
    _S1821[int(14)] = 0.0f;
    _S1821[int(15)] = 0.0f;
    _S1821[int(16)] = 0.0f;
    _S1821[int(17)] = 0.0f;
    _S1821[int(18)] = 0.0f;
    _S1821[int(19)] = 0.0f;
    _S1821[int(20)] = 0.0f;
    _S1821[int(21)] = 0.0f;
    _S1821[int(22)] = 0.0f;
    _S1821[int(23)] = 0.0f;
    _S1821[int(24)] = 0.0f;
    _S1821[int(25)] = 0.0f;
    _S1821[int(26)] = 0.0f;
    _S1821[int(27)] = 0.0f;
    _S1821[int(28)] = 0.0f;
    _S1821[int(29)] = 0.0f;
    _S1821[int(30)] = 0.0f;
    _S1821[int(31)] = 0.0f;
    _S1821[int(32)] = 0.0f;
    _S1821[int(33)] = 0.0f;
    _S1821[int(34)] = 0.0f;
    _S1821[int(35)] = 0.0f;
    _S1821[int(8)] = _S1813;
    _S1821[int(16)] = _S1793.color_params_0.b_0.x;
    _S1821[int(15)] = _S1806;
    _S1821[int(14)] = _S1807;
    _S1821[int(13)] = _S1808;
    _S1821[int(12)] = _S1809;
    _S1821[int(11)] = _S1810;
    _S1821[int(10)] = _S1811;
    _S1821[int(9)] = _S1812;
    _S1821[int(17)] = _S1793.color_params_0.b_0.y;
    _S1821[int(7)] = _S1814;
    _S1821[int(6)] = _S1815;
    _S1821[int(5)] = _S1816;
    _S1821[int(4)] = _S1817;
    _S1821[int(3)] = _S1818;
    _S1821[int(2)] = _S1819;
    _S1821[int(1)] = _S1820;
    _S1821[int(0)] = _S1667.exposure_0;
    _S1821[int(26)] = _S1803;
    _S1821[int(34)] = _S1795;
    _S1821[int(33)] = _S1796;
    _S1821[int(32)] = _S1797;
    _S1821[int(31)] = _S1798;
    _S1821[int(30)] = _S1799;
    _S1821[int(29)] = _S1800;
    _S1821[int(28)] = _S1801;
    _S1821[int(27)] = _S1802;
    _S1821[int(35)] = _S1794;
    _S1821[int(25)] = _S1804;
    _S1821[int(24)] = _S1805;
    _S1821[int(23)] = _S1793.color_params_0.n_0.y;
    _S1821[int(22)] = _S1793.color_params_0.n_0.x;
    _S1821[int(21)] = _S1793.color_params_0.g_0.y;
    _S1821[int(20)] = _S1793.color_params_0.g_0.x;
    _S1821[int(19)] = _S1793.color_params_0.r_0.y;
    _S1821[int(18)] = _S1793.color_params_0.r_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S1821;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1822, FixedArray<float, 22>  * _S1823)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S1822, _S1823);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  * params_3, FixedArray<float, 22>  * grad_out_1, FixedArray<float, 36>  * _S1824)
{
    FixedArray<float, 36>  _S1825 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = *params_3;
    (&dp_params_1)->differential_0 = _S1825;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_1, grad_out_1);
    *_S1824 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  * raw_losses_2, int num_cameras_0, FixedArray<float, 6>  * loss_weights_0, FixedArray<float, 6>  * _S1826)
{
    float _S1827;
    FixedArray<float, 6>  losses_4;
    float _S1828 = float(num_cameras_0);
    float _S1829 = (*raw_losses_2)[int(0)] / _S1828;
    for(;;)
    {
        float _S1830 = (F32_abs((_S1829)));
        if(_S1830 < 0.10000000149011612f)
        {
            _S1827 = 0.5f * _S1829 * _S1829 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1827 = _S1830 - 0.05000000074505806f;
            break;
        }
    }
    losses_4[int(0)] = _S1827;
    losses_4[int(1)] = (*raw_losses_2)[int(1)] / (3.0f * _S1828);
    losses_4[int(2)] = ((*raw_losses_2)[int(2)] + (*raw_losses_2)[int(3)] + (*raw_losses_2)[int(4)]) / (9.0f * _S1828);
    losses_4[int(3)] = ((*raw_losses_2)[int(5)] + (*raw_losses_2)[int(6)] + (*raw_losses_2)[int(7)] + (*raw_losses_2)[int(8)] + (*raw_losses_2)[int(9)]) / (5.0f * _S1828);
    float _S1831 = (*raw_losses_2)[int(10)] / _S1828;
    for(;;)
    {
        float _S1832 = (F32_abs((_S1831)));
        if(_S1832 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1831 * _S1831 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1832 - 0.00249999994412065f;
            break;
        }
    }
    float _S1833;
    float _S1834 = (*raw_losses_2)[int(11)] / _S1828;
    for(;;)
    {
        float _S1835 = (F32_abs((_S1834)));
        if(_S1835 < 0.00499999988824129f)
        {
            _S1833 = 0.5f * _S1834 * _S1834 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1833 = _S1835 - 0.00249999994412065f;
            break;
        }
    }
    float _S1836 = _S1827 + _S1833;
    float _S1837 = (*raw_losses_2)[int(12)] / _S1828;
    for(;;)
    {
        float _S1838 = (F32_abs((_S1837)));
        if(_S1838 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1837 * _S1837 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1838 - 0.00249999994412065f;
            break;
        }
    }
    float _S1839 = _S1836 + _S1827;
    float _S1840 = (*raw_losses_2)[int(13)] / _S1828;
    for(;;)
    {
        float _S1841 = (F32_abs((_S1840)));
        if(_S1841 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1840 * _S1840 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1841 - 0.00249999994412065f;
            break;
        }
    }
    float _S1842 = _S1839 + _S1827;
    float _S1843 = (*raw_losses_2)[int(14)] / _S1828;
    for(;;)
    {
        float _S1844 = (F32_abs((_S1843)));
        if(_S1844 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1843 * _S1843 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1844 - 0.00249999994412065f;
            break;
        }
    }
    float _S1845 = _S1842 + _S1827;
    float _S1846 = (*raw_losses_2)[int(15)] / _S1828;
    for(;;)
    {
        float _S1847 = (F32_abs((_S1846)));
        if(_S1847 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1846 * _S1846 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1847 - 0.00249999994412065f;
            break;
        }
    }
    float _S1848 = _S1845 + _S1827;
    float _S1849 = (*raw_losses_2)[int(16)] / _S1828;
    for(;;)
    {
        float _S1850 = (F32_abs((_S1849)));
        if(_S1850 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1849 * _S1849 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1850 - 0.00249999994412065f;
            break;
        }
    }
    float _S1851 = _S1848 + _S1827;
    float _S1852 = (*raw_losses_2)[int(17)] / _S1828;
    for(;;)
    {
        float _S1853 = (F32_abs((_S1852)));
        if(_S1853 < 0.00499999988824129f)
        {
            _S1827 = 0.5f * _S1852 * _S1852 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1827 = _S1853 - 0.00249999994412065f;
            break;
        }
    }
    float _S1854 = (_S1851 + _S1827) / 8.0f;
    float _S1855 = ((*raw_losses_2)[int(18)] + (*raw_losses_2)[int(19)] + (*raw_losses_2)[int(20)] + (*raw_losses_2)[int(21)]) / (4.0f * _S1828);
    losses_4[int(0)] = losses_4[int(0)] * (*loss_weights_0)[int(0)];
    losses_4[int(1)] = losses_4[int(1)] * (*loss_weights_0)[int(1)];
    losses_4[int(2)] = losses_4[int(2)] * (*loss_weights_0)[int(2)];
    losses_4[int(3)] = losses_4[int(3)] * (*loss_weights_0)[int(3)];
    losses_4[int(4)] = _S1854 * (*loss_weights_0)[int(4)];
    losses_4[int(5)] = _S1855 * (*loss_weights_0)[int(5)];
    *_S1826 = losses_4;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_1, FixedArray<float, 6>  * loss_weights_1, FixedArray<float, 6>  * _s_dOut_13)
{
    FixedArray<float, 22>  _S1856 = dpraw_losses_1->primal_0;
    float _S1857 = float(num_cameras_1);
    float _S1858 = dpraw_losses_1->primal_0[int(0)] / _S1857;
    bool _S1859 = (s_primal_ctx_abs_0(_S1858)) < 0.10000000149011612f;
    float _S1860;
    if(_S1859)
    {
        _S1860 = 0.5f * _S1858;
    }
    else
    {
        _S1860 = 0.0f;
    }
    float _S1861 = 3.0f * _S1857;
    float _S1862 = 9.0f * _S1857;
    float _S1863 = 5.0f * _S1857;
    float _S1864 = _S1856[int(10)] / _S1857;
    bool _S1865 = (s_primal_ctx_abs_0(_S1864)) < 0.00499999988824129f;
    float _S1866;
    if(_S1865)
    {
        _S1866 = 0.5f * _S1864;
    }
    else
    {
        _S1866 = 0.0f;
    }
    float _S1867 = _S1856[int(11)] / _S1857;
    bool _S1868 = (s_primal_ctx_abs_0(_S1867)) < 0.00499999988824129f;
    float _S1869;
    if(_S1868)
    {
        _S1869 = 0.5f * _S1867;
    }
    else
    {
        _S1869 = 0.0f;
    }
    float _S1870 = _S1856[int(12)] / _S1857;
    bool _S1871 = (s_primal_ctx_abs_0(_S1870)) < 0.00499999988824129f;
    float _S1872;
    if(_S1871)
    {
        _S1872 = 0.5f * _S1870;
    }
    else
    {
        _S1872 = 0.0f;
    }
    float _S1873 = _S1856[int(13)] / _S1857;
    bool _S1874 = (s_primal_ctx_abs_0(_S1873)) < 0.00499999988824129f;
    float _S1875;
    if(_S1874)
    {
        _S1875 = 0.5f * _S1873;
    }
    else
    {
        _S1875 = 0.0f;
    }
    float _S1876 = _S1856[int(14)] / _S1857;
    bool _S1877 = (s_primal_ctx_abs_0(_S1876)) < 0.00499999988824129f;
    float _S1878;
    if(_S1877)
    {
        _S1878 = 0.5f * _S1876;
    }
    else
    {
        _S1878 = 0.0f;
    }
    float _S1879 = _S1856[int(15)] / _S1857;
    bool _S1880 = (s_primal_ctx_abs_0(_S1879)) < 0.00499999988824129f;
    float _S1881;
    if(_S1880)
    {
        _S1881 = 0.5f * _S1879;
    }
    else
    {
        _S1881 = 0.0f;
    }
    float _S1882 = _S1856[int(16)] / _S1857;
    bool _S1883 = (s_primal_ctx_abs_0(_S1882)) < 0.00499999988824129f;
    float _S1884;
    if(_S1883)
    {
        _S1884 = 0.5f * _S1882;
    }
    else
    {
        _S1884 = 0.0f;
    }
    float _S1885 = _S1856[int(17)] / _S1857;
    bool _S1886 = (s_primal_ctx_abs_0(_S1885)) < 0.00499999988824129f;
    float _S1887;
    if(_S1886)
    {
        _S1887 = 0.5f * _S1885;
    }
    else
    {
        _S1887 = 0.0f;
    }
    float _S1888 = (*loss_weights_1)[int(3)] * (*_s_dOut_13)[int(3)];
    float _S1889 = (*loss_weights_1)[int(2)] * (*_s_dOut_13)[int(2)];
    float _S1890 = (*loss_weights_1)[int(1)] * (*_s_dOut_13)[int(1)];
    float _S1891 = (*loss_weights_1)[int(0)] * (*_s_dOut_13)[int(0)];
    float _S1892 = (*loss_weights_1)[int(5)] * (*_s_dOut_13)[int(5)] / (4.0f * _S1857);
    float _S1893 = 0.125f * ((*loss_weights_1)[int(4)] * (*_s_dOut_13)[int(4)]);
    FixedArray<float, 22>  _S1894;
    _S1894[int(0)] = 0.0f;
    _S1894[int(1)] = 0.0f;
    _S1894[int(2)] = 0.0f;
    _S1894[int(3)] = 0.0f;
    _S1894[int(4)] = 0.0f;
    _S1894[int(5)] = 0.0f;
    _S1894[int(6)] = 0.0f;
    _S1894[int(7)] = 0.0f;
    _S1894[int(8)] = 0.0f;
    _S1894[int(9)] = 0.0f;
    _S1894[int(10)] = 0.0f;
    _S1894[int(11)] = 0.0f;
    _S1894[int(12)] = 0.0f;
    _S1894[int(13)] = 0.0f;
    _S1894[int(14)] = 0.0f;
    _S1894[int(15)] = 0.0f;
    _S1894[int(16)] = 0.0f;
    _S1894[int(17)] = 0.0f;
    _S1894[int(18)] = 0.0f;
    _S1894[int(19)] = 0.0f;
    _S1894[int(20)] = 0.0f;
    _S1894[int(21)] = 0.0f;
    _S1894[int(21)] = _S1892;
    _S1894[int(20)] = _S1892;
    _S1894[int(19)] = _S1892;
    _S1894[int(18)] = _S1892;
    float _S1895 = _S1894[int(0)];
    float _S1896 = _S1894[int(1)];
    float _S1897 = _S1894[int(2)];
    float _S1898 = _S1894[int(3)];
    float _S1899 = _S1894[int(4)];
    float _S1900 = _S1894[int(5)];
    float _S1901 = _S1894[int(6)];
    float _S1902 = _S1894[int(7)];
    float _S1903 = _S1894[int(8)];
    float _S1904 = _S1894[int(9)];
    float _S1905 = _S1894[int(10)];
    float _S1906 = _S1894[int(11)];
    float _S1907 = _S1894[int(12)];
    float _S1908 = _S1894[int(13)];
    float _S1909 = _S1894[int(14)];
    float _S1910 = _S1894[int(15)];
    float _S1911 = _S1894[int(16)];
    float _S1912 = _S1894[int(17)];
    float _S1913 = _S1894[int(18)];
    float _S1914 = _S1894[int(19)];
    float _S1915 = _S1894[int(20)];
    float _S1916 = _S1894[int(21)];
    float _S1917;
    if(_S1886)
    {
        float _S1918 = 200.0f * _S1893;
        float _S1919 = _S1887 * _S1918 + 0.5f * (_S1885 * _S1918);
        _S1887 = 0.0f;
        _S1917 = _S1919;
    }
    else
    {
        _S1887 = _S1893;
        _S1917 = 0.0f;
    }
    DiffPair_float_0 _S1920;
    (&_S1920)->primal_0 = _S1885;
    (&_S1920)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1920, _S1887);
    float _S1921 = (_S1920.differential_0 + _S1917) / _S1857;
    FixedArray<float, 22>  _S1922;
    _S1922[int(0)] = 0.0f;
    _S1922[int(1)] = 0.0f;
    _S1922[int(2)] = 0.0f;
    _S1922[int(3)] = 0.0f;
    _S1922[int(4)] = 0.0f;
    _S1922[int(5)] = 0.0f;
    _S1922[int(6)] = 0.0f;
    _S1922[int(7)] = 0.0f;
    _S1922[int(8)] = 0.0f;
    _S1922[int(9)] = 0.0f;
    _S1922[int(10)] = 0.0f;
    _S1922[int(11)] = 0.0f;
    _S1922[int(12)] = 0.0f;
    _S1922[int(13)] = 0.0f;
    _S1922[int(14)] = 0.0f;
    _S1922[int(15)] = 0.0f;
    _S1922[int(16)] = 0.0f;
    _S1922[int(17)] = 0.0f;
    _S1922[int(18)] = 0.0f;
    _S1922[int(19)] = 0.0f;
    _S1922[int(20)] = 0.0f;
    _S1922[int(21)] = 0.0f;
    _S1922[int(17)] = _S1921;
    float _S1923 = _S1895 + _S1922[int(0)];
    float _S1924 = _S1896 + _S1922[int(1)];
    float _S1925 = _S1897 + _S1922[int(2)];
    float _S1926 = _S1898 + _S1922[int(3)];
    float _S1927 = _S1899 + _S1922[int(4)];
    float _S1928 = _S1900 + _S1922[int(5)];
    float _S1929 = _S1901 + _S1922[int(6)];
    float _S1930 = _S1902 + _S1922[int(7)];
    float _S1931 = _S1903 + _S1922[int(8)];
    float _S1932 = _S1904 + _S1922[int(9)];
    float _S1933 = _S1905 + _S1922[int(10)];
    float _S1934 = _S1906 + _S1922[int(11)];
    float _S1935 = _S1907 + _S1922[int(12)];
    float _S1936 = _S1908 + _S1922[int(13)];
    float _S1937 = _S1909 + _S1922[int(14)];
    float _S1938 = _S1910 + _S1922[int(15)];
    float _S1939 = _S1911 + _S1922[int(16)];
    float _S1940 = _S1912 + _S1922[int(17)];
    float _S1941 = _S1913 + _S1922[int(18)];
    float _S1942 = _S1914 + _S1922[int(19)];
    float _S1943 = _S1915 + _S1922[int(20)];
    float _S1944 = _S1916 + _S1922[int(21)];
    if(_S1883)
    {
        float _S1945 = 200.0f * _S1893;
        float _S1946 = _S1884 * _S1945 + 0.5f * (_S1882 * _S1945);
        _S1884 = 0.0f;
        _S1887 = _S1946;
    }
    else
    {
        _S1884 = _S1893;
        _S1887 = 0.0f;
    }
    DiffPair_float_0 _S1947;
    (&_S1947)->primal_0 = _S1882;
    (&_S1947)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1947, _S1884);
    float _S1948 = (_S1947.differential_0 + _S1887) / _S1857;
    FixedArray<float, 22>  _S1949;
    _S1949[int(0)] = 0.0f;
    _S1949[int(1)] = 0.0f;
    _S1949[int(2)] = 0.0f;
    _S1949[int(3)] = 0.0f;
    _S1949[int(4)] = 0.0f;
    _S1949[int(5)] = 0.0f;
    _S1949[int(6)] = 0.0f;
    _S1949[int(7)] = 0.0f;
    _S1949[int(8)] = 0.0f;
    _S1949[int(9)] = 0.0f;
    _S1949[int(10)] = 0.0f;
    _S1949[int(11)] = 0.0f;
    _S1949[int(12)] = 0.0f;
    _S1949[int(13)] = 0.0f;
    _S1949[int(14)] = 0.0f;
    _S1949[int(15)] = 0.0f;
    _S1949[int(16)] = 0.0f;
    _S1949[int(17)] = 0.0f;
    _S1949[int(18)] = 0.0f;
    _S1949[int(19)] = 0.0f;
    _S1949[int(20)] = 0.0f;
    _S1949[int(21)] = 0.0f;
    _S1949[int(16)] = _S1948;
    float _S1950 = _S1923 + _S1949[int(0)];
    float _S1951 = _S1924 + _S1949[int(1)];
    float _S1952 = _S1925 + _S1949[int(2)];
    float _S1953 = _S1926 + _S1949[int(3)];
    float _S1954 = _S1927 + _S1949[int(4)];
    float _S1955 = _S1928 + _S1949[int(5)];
    float _S1956 = _S1929 + _S1949[int(6)];
    float _S1957 = _S1930 + _S1949[int(7)];
    float _S1958 = _S1931 + _S1949[int(8)];
    float _S1959 = _S1932 + _S1949[int(9)];
    float _S1960 = _S1933 + _S1949[int(10)];
    float _S1961 = _S1934 + _S1949[int(11)];
    float _S1962 = _S1935 + _S1949[int(12)];
    float _S1963 = _S1936 + _S1949[int(13)];
    float _S1964 = _S1937 + _S1949[int(14)];
    float _S1965 = _S1938 + _S1949[int(15)];
    float _S1966 = _S1939 + _S1949[int(16)];
    float _S1967 = _S1940 + _S1949[int(17)];
    float _S1968 = _S1941 + _S1949[int(18)];
    float _S1969 = _S1942 + _S1949[int(19)];
    float _S1970 = _S1943 + _S1949[int(20)];
    float _S1971 = _S1944 + _S1949[int(21)];
    if(_S1880)
    {
        float _S1972 = 200.0f * _S1893;
        float _S1973 = _S1881 * _S1972 + 0.5f * (_S1879 * _S1972);
        _S1881 = 0.0f;
        _S1884 = _S1973;
    }
    else
    {
        _S1881 = _S1893;
        _S1884 = 0.0f;
    }
    DiffPair_float_0 _S1974;
    (&_S1974)->primal_0 = _S1879;
    (&_S1974)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1974, _S1881);
    float _S1975 = (_S1974.differential_0 + _S1884) / _S1857;
    FixedArray<float, 22>  _S1976;
    _S1976[int(0)] = 0.0f;
    _S1976[int(1)] = 0.0f;
    _S1976[int(2)] = 0.0f;
    _S1976[int(3)] = 0.0f;
    _S1976[int(4)] = 0.0f;
    _S1976[int(5)] = 0.0f;
    _S1976[int(6)] = 0.0f;
    _S1976[int(7)] = 0.0f;
    _S1976[int(8)] = 0.0f;
    _S1976[int(9)] = 0.0f;
    _S1976[int(10)] = 0.0f;
    _S1976[int(11)] = 0.0f;
    _S1976[int(12)] = 0.0f;
    _S1976[int(13)] = 0.0f;
    _S1976[int(14)] = 0.0f;
    _S1976[int(15)] = 0.0f;
    _S1976[int(16)] = 0.0f;
    _S1976[int(17)] = 0.0f;
    _S1976[int(18)] = 0.0f;
    _S1976[int(19)] = 0.0f;
    _S1976[int(20)] = 0.0f;
    _S1976[int(21)] = 0.0f;
    _S1976[int(15)] = _S1975;
    float _S1977 = _S1950 + _S1976[int(0)];
    float _S1978 = _S1951 + _S1976[int(1)];
    float _S1979 = _S1952 + _S1976[int(2)];
    float _S1980 = _S1953 + _S1976[int(3)];
    float _S1981 = _S1954 + _S1976[int(4)];
    float _S1982 = _S1955 + _S1976[int(5)];
    float _S1983 = _S1956 + _S1976[int(6)];
    float _S1984 = _S1957 + _S1976[int(7)];
    float _S1985 = _S1958 + _S1976[int(8)];
    float _S1986 = _S1959 + _S1976[int(9)];
    float _S1987 = _S1960 + _S1976[int(10)];
    float _S1988 = _S1961 + _S1976[int(11)];
    float _S1989 = _S1962 + _S1976[int(12)];
    float _S1990 = _S1963 + _S1976[int(13)];
    float _S1991 = _S1964 + _S1976[int(14)];
    float _S1992 = _S1965 + _S1976[int(15)];
    float _S1993 = _S1966 + _S1976[int(16)];
    float _S1994 = _S1967 + _S1976[int(17)];
    float _S1995 = _S1968 + _S1976[int(18)];
    float _S1996 = _S1969 + _S1976[int(19)];
    float _S1997 = _S1970 + _S1976[int(20)];
    float _S1998 = _S1971 + _S1976[int(21)];
    if(_S1877)
    {
        float _S1999 = 200.0f * _S1893;
        float _S2000 = _S1878 * _S1999 + 0.5f * (_S1876 * _S1999);
        _S1878 = 0.0f;
        _S1881 = _S2000;
    }
    else
    {
        _S1878 = _S1893;
        _S1881 = 0.0f;
    }
    DiffPair_float_0 _S2001;
    (&_S2001)->primal_0 = _S1876;
    (&_S2001)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2001, _S1878);
    float _S2002 = (_S2001.differential_0 + _S1881) / _S1857;
    FixedArray<float, 22>  _S2003;
    _S2003[int(0)] = 0.0f;
    _S2003[int(1)] = 0.0f;
    _S2003[int(2)] = 0.0f;
    _S2003[int(3)] = 0.0f;
    _S2003[int(4)] = 0.0f;
    _S2003[int(5)] = 0.0f;
    _S2003[int(6)] = 0.0f;
    _S2003[int(7)] = 0.0f;
    _S2003[int(8)] = 0.0f;
    _S2003[int(9)] = 0.0f;
    _S2003[int(10)] = 0.0f;
    _S2003[int(11)] = 0.0f;
    _S2003[int(12)] = 0.0f;
    _S2003[int(13)] = 0.0f;
    _S2003[int(14)] = 0.0f;
    _S2003[int(15)] = 0.0f;
    _S2003[int(16)] = 0.0f;
    _S2003[int(17)] = 0.0f;
    _S2003[int(18)] = 0.0f;
    _S2003[int(19)] = 0.0f;
    _S2003[int(20)] = 0.0f;
    _S2003[int(21)] = 0.0f;
    _S2003[int(14)] = _S2002;
    float _S2004 = _S1977 + _S2003[int(0)];
    float _S2005 = _S1978 + _S2003[int(1)];
    float _S2006 = _S1979 + _S2003[int(2)];
    float _S2007 = _S1980 + _S2003[int(3)];
    float _S2008 = _S1981 + _S2003[int(4)];
    float _S2009 = _S1982 + _S2003[int(5)];
    float _S2010 = _S1983 + _S2003[int(6)];
    float _S2011 = _S1984 + _S2003[int(7)];
    float _S2012 = _S1985 + _S2003[int(8)];
    float _S2013 = _S1986 + _S2003[int(9)];
    float _S2014 = _S1987 + _S2003[int(10)];
    float _S2015 = _S1988 + _S2003[int(11)];
    float _S2016 = _S1989 + _S2003[int(12)];
    float _S2017 = _S1990 + _S2003[int(13)];
    float _S2018 = _S1991 + _S2003[int(14)];
    float _S2019 = _S1992 + _S2003[int(15)];
    float _S2020 = _S1993 + _S2003[int(16)];
    float _S2021 = _S1994 + _S2003[int(17)];
    float _S2022 = _S1995 + _S2003[int(18)];
    float _S2023 = _S1996 + _S2003[int(19)];
    float _S2024 = _S1997 + _S2003[int(20)];
    float _S2025 = _S1998 + _S2003[int(21)];
    if(_S1874)
    {
        float _S2026 = 200.0f * _S1893;
        float _S2027 = _S1875 * _S2026 + 0.5f * (_S1873 * _S2026);
        _S1875 = 0.0f;
        _S1878 = _S2027;
    }
    else
    {
        _S1875 = _S1893;
        _S1878 = 0.0f;
    }
    DiffPair_float_0 _S2028;
    (&_S2028)->primal_0 = _S1873;
    (&_S2028)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2028, _S1875);
    float _S2029 = (_S2028.differential_0 + _S1878) / _S1857;
    FixedArray<float, 22>  _S2030;
    _S2030[int(0)] = 0.0f;
    _S2030[int(1)] = 0.0f;
    _S2030[int(2)] = 0.0f;
    _S2030[int(3)] = 0.0f;
    _S2030[int(4)] = 0.0f;
    _S2030[int(5)] = 0.0f;
    _S2030[int(6)] = 0.0f;
    _S2030[int(7)] = 0.0f;
    _S2030[int(8)] = 0.0f;
    _S2030[int(9)] = 0.0f;
    _S2030[int(10)] = 0.0f;
    _S2030[int(11)] = 0.0f;
    _S2030[int(12)] = 0.0f;
    _S2030[int(13)] = 0.0f;
    _S2030[int(14)] = 0.0f;
    _S2030[int(15)] = 0.0f;
    _S2030[int(16)] = 0.0f;
    _S2030[int(17)] = 0.0f;
    _S2030[int(18)] = 0.0f;
    _S2030[int(19)] = 0.0f;
    _S2030[int(20)] = 0.0f;
    _S2030[int(21)] = 0.0f;
    _S2030[int(13)] = _S2029;
    float _S2031 = _S2004 + _S2030[int(0)];
    float _S2032 = _S2005 + _S2030[int(1)];
    float _S2033 = _S2006 + _S2030[int(2)];
    float _S2034 = _S2007 + _S2030[int(3)];
    float _S2035 = _S2008 + _S2030[int(4)];
    float _S2036 = _S2009 + _S2030[int(5)];
    float _S2037 = _S2010 + _S2030[int(6)];
    float _S2038 = _S2011 + _S2030[int(7)];
    float _S2039 = _S2012 + _S2030[int(8)];
    float _S2040 = _S2013 + _S2030[int(9)];
    float _S2041 = _S2014 + _S2030[int(10)];
    float _S2042 = _S2015 + _S2030[int(11)];
    float _S2043 = _S2016 + _S2030[int(12)];
    float _S2044 = _S2017 + _S2030[int(13)];
    float _S2045 = _S2018 + _S2030[int(14)];
    float _S2046 = _S2019 + _S2030[int(15)];
    float _S2047 = _S2020 + _S2030[int(16)];
    float _S2048 = _S2021 + _S2030[int(17)];
    float _S2049 = _S2022 + _S2030[int(18)];
    float _S2050 = _S2023 + _S2030[int(19)];
    float _S2051 = _S2024 + _S2030[int(20)];
    float _S2052 = _S2025 + _S2030[int(21)];
    if(_S1871)
    {
        float _S2053 = 200.0f * _S1893;
        float _S2054 = _S1872 * _S2053 + 0.5f * (_S1870 * _S2053);
        _S1872 = 0.0f;
        _S1875 = _S2054;
    }
    else
    {
        _S1872 = _S1893;
        _S1875 = 0.0f;
    }
    DiffPair_float_0 _S2055;
    (&_S2055)->primal_0 = _S1870;
    (&_S2055)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2055, _S1872);
    float _S2056 = (_S2055.differential_0 + _S1875) / _S1857;
    FixedArray<float, 22>  _S2057;
    _S2057[int(0)] = 0.0f;
    _S2057[int(1)] = 0.0f;
    _S2057[int(2)] = 0.0f;
    _S2057[int(3)] = 0.0f;
    _S2057[int(4)] = 0.0f;
    _S2057[int(5)] = 0.0f;
    _S2057[int(6)] = 0.0f;
    _S2057[int(7)] = 0.0f;
    _S2057[int(8)] = 0.0f;
    _S2057[int(9)] = 0.0f;
    _S2057[int(10)] = 0.0f;
    _S2057[int(11)] = 0.0f;
    _S2057[int(12)] = 0.0f;
    _S2057[int(13)] = 0.0f;
    _S2057[int(14)] = 0.0f;
    _S2057[int(15)] = 0.0f;
    _S2057[int(16)] = 0.0f;
    _S2057[int(17)] = 0.0f;
    _S2057[int(18)] = 0.0f;
    _S2057[int(19)] = 0.0f;
    _S2057[int(20)] = 0.0f;
    _S2057[int(21)] = 0.0f;
    _S2057[int(12)] = _S2056;
    float _S2058 = _S2031 + _S2057[int(0)];
    float _S2059 = _S2032 + _S2057[int(1)];
    float _S2060 = _S2033 + _S2057[int(2)];
    float _S2061 = _S2034 + _S2057[int(3)];
    float _S2062 = _S2035 + _S2057[int(4)];
    float _S2063 = _S2036 + _S2057[int(5)];
    float _S2064 = _S2037 + _S2057[int(6)];
    float _S2065 = _S2038 + _S2057[int(7)];
    float _S2066 = _S2039 + _S2057[int(8)];
    float _S2067 = _S2040 + _S2057[int(9)];
    float _S2068 = _S2041 + _S2057[int(10)];
    float _S2069 = _S2042 + _S2057[int(11)];
    float _S2070 = _S2043 + _S2057[int(12)];
    float _S2071 = _S2044 + _S2057[int(13)];
    float _S2072 = _S2045 + _S2057[int(14)];
    float _S2073 = _S2046 + _S2057[int(15)];
    float _S2074 = _S2047 + _S2057[int(16)];
    float _S2075 = _S2048 + _S2057[int(17)];
    float _S2076 = _S2049 + _S2057[int(18)];
    float _S2077 = _S2050 + _S2057[int(19)];
    float _S2078 = _S2051 + _S2057[int(20)];
    float _S2079 = _S2052 + _S2057[int(21)];
    if(_S1868)
    {
        float _S2080 = 200.0f * _S1893;
        float _S2081 = _S1869 * _S2080 + 0.5f * (_S1867 * _S2080);
        _S1869 = 0.0f;
        _S1872 = _S2081;
    }
    else
    {
        _S1869 = _S1893;
        _S1872 = 0.0f;
    }
    DiffPair_float_0 _S2082;
    (&_S2082)->primal_0 = _S1867;
    (&_S2082)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2082, _S1869);
    float _S2083 = (_S2082.differential_0 + _S1872) / _S1857;
    FixedArray<float, 22>  _S2084;
    _S2084[int(0)] = 0.0f;
    _S2084[int(1)] = 0.0f;
    _S2084[int(2)] = 0.0f;
    _S2084[int(3)] = 0.0f;
    _S2084[int(4)] = 0.0f;
    _S2084[int(5)] = 0.0f;
    _S2084[int(6)] = 0.0f;
    _S2084[int(7)] = 0.0f;
    _S2084[int(8)] = 0.0f;
    _S2084[int(9)] = 0.0f;
    _S2084[int(10)] = 0.0f;
    _S2084[int(11)] = 0.0f;
    _S2084[int(12)] = 0.0f;
    _S2084[int(13)] = 0.0f;
    _S2084[int(14)] = 0.0f;
    _S2084[int(15)] = 0.0f;
    _S2084[int(16)] = 0.0f;
    _S2084[int(17)] = 0.0f;
    _S2084[int(18)] = 0.0f;
    _S2084[int(19)] = 0.0f;
    _S2084[int(20)] = 0.0f;
    _S2084[int(21)] = 0.0f;
    _S2084[int(11)] = _S2083;
    float _S2085 = _S2058 + _S2084[int(0)];
    float _S2086 = _S2059 + _S2084[int(1)];
    float _S2087 = _S2060 + _S2084[int(2)];
    float _S2088 = _S2061 + _S2084[int(3)];
    float _S2089 = _S2062 + _S2084[int(4)];
    float _S2090 = _S2063 + _S2084[int(5)];
    float _S2091 = _S2064 + _S2084[int(6)];
    float _S2092 = _S2065 + _S2084[int(7)];
    float _S2093 = _S2066 + _S2084[int(8)];
    float _S2094 = _S2067 + _S2084[int(9)];
    float _S2095 = _S2068 + _S2084[int(10)];
    float _S2096 = _S2069 + _S2084[int(11)];
    float _S2097 = _S2070 + _S2084[int(12)];
    float _S2098 = _S2071 + _S2084[int(13)];
    float _S2099 = _S2072 + _S2084[int(14)];
    float _S2100 = _S2073 + _S2084[int(15)];
    float _S2101 = _S2074 + _S2084[int(16)];
    float _S2102 = _S2075 + _S2084[int(17)];
    float _S2103 = _S2076 + _S2084[int(18)];
    float _S2104 = _S2077 + _S2084[int(19)];
    float _S2105 = _S2078 + _S2084[int(20)];
    float _S2106 = _S2079 + _S2084[int(21)];
    if(_S1865)
    {
        float _S2107 = 200.0f * _S1893;
        float _S2108 = _S1866 * _S2107 + 0.5f * (_S1864 * _S2107);
        _S1866 = 0.0f;
        _S1869 = _S2108;
    }
    else
    {
        _S1866 = _S1893;
        _S1869 = 0.0f;
    }
    DiffPair_float_0 _S2109;
    (&_S2109)->primal_0 = _S1864;
    (&_S2109)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2109, _S1866);
    float _S2110 = (_S2109.differential_0 + _S1869) / _S1857;
    float _S2111 = _S1888 / _S1863;
    float _S2112 = _S1889 / _S1862;
    float _S2113 = _S1890 / _S1861;
    FixedArray<float, 22>  _S2114;
    _S2114[int(0)] = 0.0f;
    _S2114[int(1)] = 0.0f;
    _S2114[int(2)] = 0.0f;
    _S2114[int(3)] = 0.0f;
    _S2114[int(4)] = 0.0f;
    _S2114[int(5)] = 0.0f;
    _S2114[int(6)] = 0.0f;
    _S2114[int(7)] = 0.0f;
    _S2114[int(8)] = 0.0f;
    _S2114[int(9)] = 0.0f;
    _S2114[int(10)] = 0.0f;
    _S2114[int(11)] = 0.0f;
    _S2114[int(12)] = 0.0f;
    _S2114[int(13)] = 0.0f;
    _S2114[int(14)] = 0.0f;
    _S2114[int(15)] = 0.0f;
    _S2114[int(16)] = 0.0f;
    _S2114[int(17)] = 0.0f;
    _S2114[int(18)] = 0.0f;
    _S2114[int(19)] = 0.0f;
    _S2114[int(20)] = 0.0f;
    _S2114[int(21)] = 0.0f;
    _S2114[int(10)] = _S2110;
    _S2114[int(9)] = _S2111;
    _S2114[int(8)] = _S2111;
    _S2114[int(7)] = _S2111;
    _S2114[int(6)] = _S2111;
    _S2114[int(5)] = _S2111;
    _S2114[int(4)] = _S2112;
    _S2114[int(3)] = _S2112;
    _S2114[int(2)] = _S2112;
    _S2114[int(1)] = _S2113;
    float _S2115 = _S2085 + _S2114[int(0)];
    float _S2116 = _S2086 + _S2114[int(1)];
    float _S2117 = _S2087 + _S2114[int(2)];
    float _S2118 = _S2088 + _S2114[int(3)];
    float _S2119 = _S2089 + _S2114[int(4)];
    float _S2120 = _S2090 + _S2114[int(5)];
    float _S2121 = _S2091 + _S2114[int(6)];
    float _S2122 = _S2092 + _S2114[int(7)];
    float _S2123 = _S2093 + _S2114[int(8)];
    float _S2124 = _S2094 + _S2114[int(9)];
    float _S2125 = _S2095 + _S2114[int(10)];
    float _S2126 = _S2096 + _S2114[int(11)];
    float _S2127 = _S2097 + _S2114[int(12)];
    float _S2128 = _S2098 + _S2114[int(13)];
    float _S2129 = _S2099 + _S2114[int(14)];
    float _S2130 = _S2100 + _S2114[int(15)];
    float _S2131 = _S2101 + _S2114[int(16)];
    float _S2132 = _S2102 + _S2114[int(17)];
    float _S2133 = _S2103 + _S2114[int(18)];
    float _S2134 = _S2104 + _S2114[int(19)];
    float _S2135 = _S2105 + _S2114[int(20)];
    float _S2136 = _S2106 + _S2114[int(21)];
    if(_S1859)
    {
        float _S2137 = 10.0f * _S1891;
        float _S2138 = _S1860 * _S2137 + 0.5f * (_S1858 * _S2137);
        _S1860 = 0.0f;
        _S1866 = _S2138;
    }
    else
    {
        _S1860 = _S1891;
        _S1866 = 0.0f;
    }
    DiffPair_float_0 _S2139;
    (&_S2139)->primal_0 = _S1858;
    (&_S2139)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2139, _S1860);
    float _S2140 = (_S2139.differential_0 + _S1866) / _S1857;
    FixedArray<float, 22>  _S2141;
    _S2141[int(0)] = 0.0f;
    _S2141[int(1)] = 0.0f;
    _S2141[int(2)] = 0.0f;
    _S2141[int(3)] = 0.0f;
    _S2141[int(4)] = 0.0f;
    _S2141[int(5)] = 0.0f;
    _S2141[int(6)] = 0.0f;
    _S2141[int(7)] = 0.0f;
    _S2141[int(8)] = 0.0f;
    _S2141[int(9)] = 0.0f;
    _S2141[int(10)] = 0.0f;
    _S2141[int(11)] = 0.0f;
    _S2141[int(12)] = 0.0f;
    _S2141[int(13)] = 0.0f;
    _S2141[int(14)] = 0.0f;
    _S2141[int(15)] = 0.0f;
    _S2141[int(16)] = 0.0f;
    _S2141[int(17)] = 0.0f;
    _S2141[int(18)] = 0.0f;
    _S2141[int(19)] = 0.0f;
    _S2141[int(20)] = 0.0f;
    _S2141[int(21)] = 0.0f;
    _S2141[int(0)] = _S2140;
    FixedArray<float, 22>  _S2142 = {
        _S2115 + _S2141[int(0)], _S2116 + _S2141[int(1)], _S2117 + _S2141[int(2)], _S2118 + _S2141[int(3)], _S2119 + _S2141[int(4)], _S2120 + _S2141[int(5)], _S2121 + _S2141[int(6)], _S2122 + _S2141[int(7)], _S2123 + _S2141[int(8)], _S2124 + _S2141[int(9)], _S2125 + _S2141[int(10)], _S2126 + _S2141[int(11)], _S2127 + _S2141[int(12)], _S2128 + _S2141[int(13)], _S2129 + _S2141[int(14)], _S2130 + _S2141[int(15)], _S2131 + _S2141[int(16)], _S2132 + _S2141[int(17)], _S2133 + _S2141[int(18)], _S2134 + _S2141[int(19)], _S2135 + _S2141[int(20)], _S2136 + _S2141[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2142;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S2143, int _S2144, FixedArray<float, 6>  * _S2145, FixedArray<float, 6>  * _S2146)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2143, _S2144, _S2145, _S2146);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  * raw_losses_3, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * grad_out_2, FixedArray<float, 22>  * _S2147)
{
    FixedArray<float, 22>  _S2148 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = *raw_losses_3;
    (&dp_raw_losses_1)->differential_0 = _S2148;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_2, loss_weights_2, grad_out_2);
    *_S2147 = (&dp_raw_losses_1)->differential_0;
    return;
}

