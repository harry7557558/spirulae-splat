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

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_3, float3  dOut_13)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_9).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_3).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_13.x);
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = left_dp_0.differential_0;
    float3  right_d_result_4;
    *&((&right_d_result_4)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_9).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_3).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_13.y);
    *&((&left_d_result_4)->y) = left_dp_1.differential_0;
    *&((&right_d_result_4)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_9).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_3).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_13.z);
    *&((&left_d_result_4)->z) = left_dp_2.differential_0;
    *&((&right_d_result_4)->z) = right_dp_2.differential_0;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = left_d_result_4;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  max_0(float3  x_22, float3  y_7)
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
        *_slang_vector_get_element_ptr(&result_18, i_11) = (F32_max((_slang_vector_get_element(x_22, i_11)), (_slang_vector_get_element(y_7, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_18;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_14)
{
    float3  _S418 = make_float3 (1.0f) / (*dpx_10).primal_0 * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S418;
    return;
}

inline __device__ float3  log_0(float3  x_23)
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
        *_slang_vector_get_element_ptr(&result_19, i_12) = (F32_log((_slang_vector_get_element(x_23, i_12))));
        i_12 = i_12 + int(1);
    }
    return result_19;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_11, float dOut_15)
{
    float _S419 = _slang_select(((*dpx_11).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_11).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_15;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S419;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, float3  dOut_16)
{
    float3  _S420 = _slang_select(((*dpx_12).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_12).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_16;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S420;
    return;
}

inline __device__ float3  abs_0(float3  x_24)
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
        *_slang_vector_get_element_ptr(&result_20, i_13) = (F32_abs((_slang_vector_get_element(x_24, i_13))));
        i_13 = i_13 + int(1);
    }
    return result_20;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_13, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_17)
{
    DiffPair_float_0 _S421 = *dpx_13;
    bool _S422;
    if(((*dpx_13).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S422 = ((*dpx_13).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S422 = false;
    }
    float _S423;
    if(_S422)
    {
        _S423 = dOut_17;
    }
    else
    {
        _S423 = 0.0f;
    }
    dpx_13->primal_0 = _S421.primal_0;
    dpx_13->differential_0 = _S423;
    DiffPair_float_0 _S424 = *dpMin_0;
    if((_S421.primal_0) < ((*dpMin_0).primal_0))
    {
        _S423 = dOut_17;
    }
    else
    {
        _S423 = 0.0f;
    }
    dpMin_0->primal_0 = _S424.primal_0;
    dpMin_0->differential_0 = _S423;
    DiffPair_float_0 _S425 = *dpMax_0;
    if(((*dpx_13).primal_0) > ((*dpMax_0).primal_0))
    {
        _S423 = dOut_17;
    }
    else
    {
        _S423 = 0.0f;
    }
    dpMax_0->primal_0 = _S425.primal_0;
    dpMax_0->differential_0 = _S423;
    return;
}

inline __device__ float clamp_0(float x_25, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_25), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_14, float dOut_18)
{
    float _S426 = -0.5f / ((*dpx_14).primal_0 * (F32_sqrt(((*dpx_14).primal_0)))) * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S426;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_15, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_19)
{
    float _S427 = (1.0f - (*dps_0).primal_0) * dOut_19;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S427;
    DiffPair_float_0 _S428 = *dpy_4;
    float _S429 = (*dps_0).primal_0 * dOut_19;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S429;
    float _S430 = (_S428.primal_0 - (*dpx_15).primal_0) * dOut_19;
    dps_0->primal_0 = _S428.primal_0;
    dps_0->differential_0 = _S430;
    return;
}

inline __device__ float lerp_0(float x_26, float y_8, float s_4)
{
    return x_26 + (y_8 - x_26) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 11>  * weights_0, FixedArray<float, 23>  * _S431)
{
    float3  _S432;
    bool _S433;
    bool _S434;
    FixedArray<float, 23>  losses_1;
    float _S435 = float(mask_0);
    float _S436 = 1.0f - (*weights_0)[int(0)];
    float3  _S437 = make_float3 (0.0f);
    float3  _S438 = abs_0(ref_rgb_0 * make_float3 (_S436) + (log_0(max_0(ref_rgb_0, _S437) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)]) - (render_rgb_0 * make_float3 (_S436) + (log_0(max_0(render_rgb_0, _S437) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)])));
    losses_1[int(0)] = (*weights_0)[int(1)] * _S435 * ((_S438.x + _S438.y + _S438.z) * 0.3333333432674408f);
    float3  _S439 = ref_rgb_0 - render_rgb_0;
    losses_1[int(1)] = _S435 * clamp_0(dot_0(_S439, _S439) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S440 = float(depth_mask_0 & mask_0);
    float _S441 = _S440 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S442 = _S440 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S441;
    losses_1[int(3)] = _S442;
    losses_1[int(4)] = _S441 * _S441;
    losses_1[int(5)] = _S442 * _S442;
    losses_1[int(6)] = _S441 * _S442;
    bool _S443 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S444 = norm2_0 == 0.0f;
        _S433 = _S444;
        if(_S444)
        {
            _S432 = _S437;
            break;
        }
        _S432 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S445;
    bool _S446 = !_S433;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S447 = norm2_1 == 0.0f;
        _S434 = _S447;
        if(_S447)
        {
            _S445 = _S437;
            break;
        }
        _S445 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S448;
    float3  _S449;
    bool _S450 = !_S434;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S449 = _S437;
            _S448 = false;
            break;
        }
        _S449 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S448 = _S443;
        break;
    }
    float _S451 = float(_S446 & _S448);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S432, _S449);
    losses_1[int(7)] = (*weights_0)[int(3)] * _S451 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S452 = float(_S450 & _S448);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S445, _S449);
    losses_1[int(8)] = (*weights_0)[int(3)] * _S452 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S453 = float(_S446 & _S450);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S432, _S445);
    losses_1[int(11)] = (*weights_0)[int(6)] * _S453 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S454 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S455 = float(alpha_mask_0);
    float _S456 = float(ref_alpha_0);
    float _S457 = (F32_max((_S454), (_S456)));
    losses_1[int(9)] = (*weights_0)[int(4)] * _S455 * - lerp_0((F32_log(((F32_max((1.0f - _S457), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S457), (9.99999997475242708e-07f)))))), _S456);
    float _S458 = 1.0f - _S454;
    float _S459 = 1.0f - _S456;
    float _S460 = (F32_max((_S458), (_S459)));
    losses_1[int(10)] = (*weights_0)[int(5)] * _S455 * - lerp_0((F32_log(((F32_max((1.0f - _S460), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S460), (9.99999997475242708e-07f)))))), _S459);
    losses_1[int(12)] = (*weights_0)[int(7)] * 4.0f * _S454 * _S458;
    float _S461 = (F32_max((_S454), (9.999999960041972e-13f)));
    losses_1[int(13)] = (*weights_0)[int(8)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S461;
    losses_1[int(14)] = (*weights_0)[int(9)] * depth_dist_0 / _S461;
    losses_1[int(15)] = (*weights_0)[int(10)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S461;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S435;
    losses_1[int(18)] = _S440;
    losses_1[int(19)] = _S451;
    losses_1[int(20)] = _S452;
    losses_1[int(21)] = _S453;
    losses_1[int(22)] = _S455;
    *_S431 = losses_1;
    return;
}

inline __device__ float3  s_primal_ctx_max_1(float3  _S462, float3  _S463)
{
    return max_0(_S462, _S463);
}

inline __device__ float3  s_primal_ctx_log_1(float3  _S464)
{
    return log_0(_S464);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S465, float3  _S466)
{
    return dot_0(_S465, _S466);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S467)
{
    return (F32_rsqrt((_S467)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S468, float _S469, float _S470)
{
    return clamp_0(_S468, _S469, _S470);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S471, DiffPair_float_0 * _S472, DiffPair_float_0 * _S473, float _S474)
{
    _d_lerp_0(_S471, _S472, _S473, _S474);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S475, DiffPair_float_0 * _S476, DiffPair_float_0 * _S477, float _S478)
{
    _d_clamp_0(_S475, _S476, _S477, _S478);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S479, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S480, float _S481)
{
    _d_dot_0(_S479, _S480, _S481);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S482, float _S483)
{
    _d_rsqrt_0(_S482, _S483);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S484, float3  _S485)
{
    _d_abs_vector_0(_S484, _S485);
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S486, float3  _S487)
{
    _d_log_vector_0(_S486, _S487);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S488, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S489, float3  _S490)
{
    _d_max_vector_0(_S488, _S489, _S490);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 11>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S491 = *dprender_rgb_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S492 = *dpref_rgb_0;
    DiffPair_float_0 _S493 = *dprender_depth_0;
    DiffPair_float_0 _S494 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S495 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S496 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S497 = *dpref_normal_0;
    DiffPair_float_0 _S498 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S499 = *dprgb_dist_0;
    DiffPair_float_0 _S500 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S501 = *dpnormal_dist_0;
    float3  _S502 = make_float3 (0.0f);
    float _S503 = float(mask_1);
    float _S504 = (*weights_1)[int(1)] * _S503;
    float3  _S505 = make_float3 ((*weights_1)[int(0)]);
    float _S506 = 1.0f - (*weights_1)[int(0)];
    float3  _S507 = make_float3 (_S506);
    float3  _S508 = make_float3 (0.0f);
    float3  _S509 = s_primal_ctx_max_1((*dprender_rgb_0).primal_0, _S508) + make_float3 (0.00392156885936856f);
    float3  _S510 = s_primal_ctx_max_1((*dpref_rgb_0).primal_0, _S508) + make_float3 (0.00392156885936856f);
    float3  _S511 = (*dpref_rgb_0).primal_0 * make_float3 (_S506) + (s_primal_ctx_log_1(_S510) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]) - ((*dprender_rgb_0).primal_0 * make_float3 (_S506) + (s_primal_ctx_log_1(_S509) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]));
    float3  _S512 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S513 = s_primal_ctx_dot_0(_S512, _S512) * 0.3333333432674408f;
    float _S514 = float(depth_mask_1 & mask_1);
    float _S515 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S516 = _S514 * s_primal_ctx_log_0(_S515);
    float _S517 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S518 = _S514 * s_primal_ctx_log_0(_S517);
    bool _S519 = normal_mask_1 & mask_1;
    float _S520 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S521 = !(_S520 == 0.0f);
    float3  _S522;
    float3  _S523;
    if(_S521)
    {
        float _S524 = s_primal_ctx_rsqrt_0(_S520);
        float3  _S525 = make_float3 (_S524);
        _S522 = _S495.primal_0 * make_float3 (_S524);
        _S523 = _S525;
    }
    else
    {
        _S522 = _S508;
        _S523 = _S502;
    }
    float _S526 = s_primal_ctx_dot_0(_S496.primal_0, _S496.primal_0);
    bool _S527 = !(_S526 == 0.0f);
    float3  _S528;
    float3  _S529;
    if(_S527)
    {
        float _S530 = s_primal_ctx_rsqrt_0(_S526);
        float3  _S531 = make_float3 (_S530);
        _S528 = _S496.primal_0 * make_float3 (_S530);
        _S529 = _S531;
    }
    else
    {
        _S528 = _S508;
        _S529 = _S502;
    }
    float _S532 = s_primal_ctx_dot_0(_S497.primal_0, _S497.primal_0);
    bool _S533 = _S532 == 0.0f;
    bool _S534;
    if(_S533)
    {
        _S534 = false;
    }
    else
    {
        _S534 = _S519;
    }
    bool _S535 = !_S533;
    float3  _S536;
    float3  _S537;
    if(_S535)
    {
        float _S538 = s_primal_ctx_rsqrt_0(_S532);
        float3  _S539 = make_float3 (_S538);
        _S536 = _S497.primal_0 * make_float3 (_S538);
        _S537 = _S539;
    }
    else
    {
        _S536 = _S508;
        _S537 = _S502;
    }
    float _S540 = (*weights_1)[int(3)] * float(_S521 & _S534);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S522, _S536);
    float _S541 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S542 = (*weights_1)[int(3)] * float(_S527 & _S534);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S528, _S536);
    float _S543 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S544 = (*weights_1)[int(6)] * float(_S521 & _S527);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S522, _S528);
    float _S545 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S546 = s_primal_ctx_clamp_0(_S498.primal_0, 0.0f, 1.0f);
    float _S547 = float(alpha_mask_1);
    float _S548 = (*weights_1)[int(4)] * _S547;
    float _S549 = float(ref_alpha_1);
    float _S550 = s_primal_ctx_max_0(_S546, _S549);
    float _S551 = 1.0f - _S550;
    float _S552 = s_primal_ctx_max_0(_S551, 9.99999997475242708e-07f);
    float _S553 = s_primal_ctx_log_0(_S552);
    float _S554 = s_primal_ctx_max_0(_S550, 9.99999997475242708e-07f);
    float _S555 = s_primal_ctx_log_0(_S554);
    float _S556 = (*weights_1)[int(5)] * _S547;
    float _S557 = 1.0f - _S546;
    float _S558 = 1.0f - _S549;
    float _S559 = s_primal_ctx_max_0(_S557, _S558);
    float _S560 = 1.0f - _S559;
    float _S561 = s_primal_ctx_max_0(_S560, 9.99999997475242708e-07f);
    float _S562 = s_primal_ctx_log_0(_S561);
    float _S563 = s_primal_ctx_max_0(_S559, 9.99999997475242708e-07f);
    float _S564 = s_primal_ctx_log_0(_S563);
    float _S565 = (*weights_1)[int(7)] * 4.0f;
    float _S566 = _S565 * _S546;
    float _S567 = s_primal_ctx_max_0(_S546, 9.999999960041972e-13f);
    float _S568 = _S567 * _S567;
    float _S569 = (*_s_dOut_4)[int(15)] / _S568;
    float _S570 = 0.3333333432674408f * ((*weights_1)[int(10)] * (_S567 * _S569));
    float _S571 = (*_s_dOut_4)[int(14)] / _S568;
    float _S572 = (*weights_1)[int(9)] * (_S567 * _S571);
    float _S573 = (*_s_dOut_4)[int(13)] / _S568;
    float _S574 = _S567 * _S573;
    float _S575 = (*weights_1)[int(10)] * ((_S501.primal_0.x + _S501.primal_0.y + _S501.primal_0.z) * 0.3333333432674408f) * - _S569 + (*weights_1)[int(9)] * _S500.primal_0 * - _S571 + (*weights_1)[int(8)] * ((_S499.primal_0.x + _S499.primal_0.y + _S499.primal_0.z) * 0.3333333432674408f) * - _S573;
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S546;
    (&_S576)->differential_0 = 0.0f;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = 9.999999960041972e-13f;
    (&_S577)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S576, &_S577, _S575);
    float _S578 = 0.3333333432674408f * ((*weights_1)[int(8)] * _S574);
    float _S579 = _S566 * (*_s_dOut_4)[int(12)];
    float _S580 = _S565 * (_S557 * (*_s_dOut_4)[int(12)]);
    float _S581 = - (_S556 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S582;
    (&_S582)->primal_0 = _S562;
    (&_S582)->differential_0 = 0.0f;
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S564;
    (&_S583)->differential_0 = 0.0f;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S558;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S582, &_S583, &_S584, _S581);
    DiffPair_float_0 _S585;
    (&_S585)->primal_0 = _S563;
    (&_S585)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S585, _S583.differential_0);
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = _S559;
    (&_S586)->differential_0 = 0.0f;
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = 9.99999997475242708e-07f;
    (&_S587)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S586, &_S587, _S585.differential_0);
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = _S561;
    (&_S588)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S588, _S582.differential_0);
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = _S560;
    (&_S589)->differential_0 = 0.0f;
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = 9.99999997475242708e-07f;
    (&_S590)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S589, &_S590, _S588.differential_0);
    float _S591 = _S586.differential_0 + - _S589.differential_0;
    DiffPair_float_0 _S592;
    (&_S592)->primal_0 = _S557;
    (&_S592)->differential_0 = 0.0f;
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = _S558;
    (&_S593)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S592, &_S593, _S591);
    float _S594 = - (_S579 + _S592.differential_0);
    float _S595 = - (_S548 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S596;
    (&_S596)->primal_0 = _S553;
    (&_S596)->differential_0 = 0.0f;
    DiffPair_float_0 _S597;
    (&_S597)->primal_0 = _S555;
    (&_S597)->differential_0 = 0.0f;
    DiffPair_float_0 _S598;
    (&_S598)->primal_0 = _S549;
    (&_S598)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S596, &_S597, &_S598, _S595);
    DiffPair_float_0 _S599;
    (&_S599)->primal_0 = _S554;
    (&_S599)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S599, _S597.differential_0);
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = _S550;
    (&_S600)->differential_0 = 0.0f;
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = 9.99999997475242708e-07f;
    (&_S601)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S600, &_S601, _S599.differential_0);
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = _S552;
    (&_S602)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S602, _S596.differential_0);
    DiffPair_float_0 _S603;
    (&_S603)->primal_0 = _S551;
    (&_S603)->differential_0 = 0.0f;
    DiffPair_float_0 _S604;
    (&_S604)->primal_0 = 9.99999997475242708e-07f;
    (&_S604)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S603, &_S604, _S602.differential_0);
    float _S605 = _S600.differential_0 + - _S603.differential_0;
    DiffPair_float_0 _S606;
    (&_S606)->primal_0 = _S546;
    (&_S606)->differential_0 = 0.0f;
    DiffPair_float_0 _S607;
    (&_S607)->primal_0 = _S549;
    (&_S607)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S606, &_S607, _S605);
    float _S608 = _S576.differential_0 + _S580 + _S594 + _S606.differential_0;
    DiffPair_float_0 _S609;
    (&_S609)->primal_0 = _S498.primal_0;
    (&_S609)->differential_0 = 0.0f;
    DiffPair_float_0 _S610;
    (&_S610)->primal_0 = 0.0f;
    (&_S610)->differential_0 = 0.0f;
    DiffPair_float_0 _S611;
    (&_S611)->primal_0 = 1.0f;
    (&_S611)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S609, &_S610, &_S611, _S608);
    DiffPair_float_0 _S612 = _S609;
    float _S613 = _S544 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S545;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S614, _S613);
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = cos_sim_loss_5;
    (&_S615)->differential_0 = 0.0f;
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = 9.999999960041972e-13f;
    (&_S616)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S615, &_S616, _S614.differential_0);
    float _S617 = 0.5f * - (_S613 + _S615.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S618;
    (&_S618)->primal_0 = _S522;
    (&_S618)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S619;
    (&_S619)->primal_0 = _S528;
    (&_S619)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S618, &_S619, _S617);
    float _S620 = _S542 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S543;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S621, _S620);
    DiffPair_float_0 _S622;
    (&_S622)->primal_0 = cos_sim_loss_4;
    (&_S622)->differential_0 = 0.0f;
    DiffPair_float_0 _S623;
    (&_S623)->primal_0 = 9.999999960041972e-13f;
    (&_S623)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S622, &_S623, _S621.differential_0);
    float _S624 = 0.5f * - (_S620 + _S622.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S625;
    (&_S625)->primal_0 = _S528;
    (&_S625)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S626;
    (&_S626)->primal_0 = _S536;
    (&_S626)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S625, &_S626, _S624);
    float _S627 = _S540 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = _S541;
    (&_S628)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S628, _S627);
    DiffPair_float_0 _S629;
    (&_S629)->primal_0 = cos_sim_loss_3;
    (&_S629)->differential_0 = 0.0f;
    DiffPair_float_0 _S630;
    (&_S630)->primal_0 = 9.999999960041972e-13f;
    (&_S630)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S629, &_S630, _S628.differential_0);
    float _S631 = 0.5f * - (_S627 + _S629.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
    (&_S632)->primal_0 = _S522;
    (&_S632)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S633;
    (&_S633)->primal_0 = _S536;
    (&_S633)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S632, &_S633, _S631);
    float3  _S634 = _S626.differential_0 + _S633.differential_0;
    float3  _S635 = _S618.differential_0 + _S632.differential_0;
    float3  _S636 = make_float3 (_S570, _S570, _S570);
    float3  _S637 = make_float3 (_S578, _S578, _S578);
    float3  _S638 = _S619.differential_0 + _S625.differential_0;
    float _S639;
    if(_S535)
    {
        float3  _S640 = _S497.primal_0 * _S634;
        float3  _S641 = _S537 * _S634;
        float _S642 = _S640.x + _S640.y + _S640.z;
        DiffPair_float_0 _S643;
        (&_S643)->primal_0 = _S532;
        (&_S643)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S643, _S642);
        _S639 = _S643.differential_0;
        _S522 = _S641;
    }
    else
    {
        _S639 = 0.0f;
        _S522 = _S502;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S644;
    (&_S644)->primal_0 = _S497.primal_0;
    (&_S644)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S645;
    (&_S645)->primal_0 = _S497.primal_0;
    (&_S645)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S644, &_S645, _S639);
    float3  _S646 = _S645.differential_0 + _S644.differential_0 + _S522;
    if(_S527)
    {
        float3  _S647 = _S496.primal_0 * _S638;
        float3  _S648 = _S529 * _S638;
        float _S649 = _S647.x + _S647.y + _S647.z;
        DiffPair_float_0 _S650;
        (&_S650)->primal_0 = _S526;
        (&_S650)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S650, _S649);
        _S639 = _S650.differential_0;
        _S522 = _S648;
    }
    else
    {
        _S639 = 0.0f;
        _S522 = _S502;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S651;
    (&_S651)->primal_0 = _S496.primal_0;
    (&_S651)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S652;
    (&_S652)->primal_0 = _S496.primal_0;
    (&_S652)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S651, &_S652, _S639);
    float3  _S653 = _S652.differential_0 + _S651.differential_0 + _S522;
    if(_S521)
    {
        float3  _S654 = _S495.primal_0 * _S635;
        float3  _S655 = _S523 * _S635;
        float _S656 = _S654.x + _S654.y + _S654.z;
        DiffPair_float_0 _S657;
        (&_S657)->primal_0 = _S520;
        (&_S657)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S657, _S656);
        _S639 = _S657.differential_0;
        _S522 = _S655;
    }
    else
    {
        _S639 = 0.0f;
        _S522 = _S502;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S658;
    (&_S658)->primal_0 = _S495.primal_0;
    (&_S658)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S659;
    (&_S659)->primal_0 = _S495.primal_0;
    (&_S659)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S658, &_S659, _S639);
    float _S660 = _S518 * (*_s_dOut_4)[int(6)];
    float _S661 = _S518 * (*_s_dOut_4)[int(5)];
    float _S662 = _S516 * (*_s_dOut_4)[int(4)];
    float _S663 = _S514 * (_S516 * (*_s_dOut_4)[int(6)] + _S661 + _S661 + (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S664;
    (&_S664)->primal_0 = _S517;
    (&_S664)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S664, _S663);
    DiffPair_float_0 _S665;
    (&_S665)->primal_0 = _S494.primal_0;
    (&_S665)->differential_0 = 0.0f;
    DiffPair_float_0 _S666;
    (&_S666)->primal_0 = 0.00009999999747379f;
    (&_S666)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S665, &_S666, _S664.differential_0);
    float _S667 = _S514 * (_S660 + _S662 + _S662 + (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = _S515;
    (&_S668)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S668, _S667);
    DiffPair_float_0 _S669;
    (&_S669)->primal_0 = _S493.primal_0;
    (&_S669)->differential_0 = 0.0f;
    DiffPair_float_0 _S670;
    (&_S670)->primal_0 = 0.00009999999747379f;
    (&_S670)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S669, &_S670, _S668.differential_0);
    float _S671 = _S503 * (*_s_dOut_4)[int(1)];
    DiffPair_float_0 _S672;
    (&_S672)->primal_0 = _S513;
    (&_S672)->differential_0 = 0.0f;
    DiffPair_float_0 _S673;
    (&_S673)->primal_0 = 0.0f;
    (&_S673)->differential_0 = 0.0f;
    DiffPair_float_0 _S674;
    (&_S674)->primal_0 = 1.0f;
    (&_S674)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S672, &_S673, &_S674, _S671);
    float _S675 = 0.3333333432674408f * _S672.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S676;
    (&_S676)->primal_0 = _S512;
    (&_S676)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S677;
    (&_S677)->primal_0 = _S512;
    (&_S677)->differential_0 = _S502;
    s_bwd_prop_dot_0(&_S676, &_S677, _S675);
    float3  _S678 = _S677.differential_0 + _S676.differential_0;
    float3  _S679 = - _S678;
    float _S680 = 0.3333333432674408f * (_S504 * (*_s_dOut_4)[int(0)]);
    float3  _S681 = make_float3 (_S680, _S680, _S680);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S682;
    (&_S682)->primal_0 = _S511;
    (&_S682)->differential_0 = _S502;
    s_bwd_prop_abs_0(&_S682, _S681);
    float3  _S683 = - _S682.differential_0;
    float3  _S684 = _S505 * _S682.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S685;
    (&_S685)->primal_0 = _S510;
    (&_S685)->differential_0 = _S502;
    s_bwd_prop_log_1(&_S685, _S684);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S686;
    (&_S686)->primal_0 = _S492.primal_0;
    (&_S686)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S687;
    (&_S687)->primal_0 = _S508;
    (&_S687)->differential_0 = _S502;
    s_bwd_prop_max_1(&_S686, &_S687, _S685.differential_0);
    float3  _S688 = _S507 * _S682.differential_0;
    float3  _S689 = _S505 * _S683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S690;
    (&_S690)->primal_0 = _S509;
    (&_S690)->differential_0 = _S502;
    s_bwd_prop_log_1(&_S690, _S689);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S691;
    (&_S691)->primal_0 = _S491.primal_0;
    (&_S691)->differential_0 = _S502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S692;
    (&_S692)->primal_0 = _S508;
    (&_S692)->differential_0 = _S502;
    s_bwd_prop_max_1(&_S691, &_S692, _S690.differential_0);
    float3  _S693 = _S507 * _S683;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S636;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S572;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S637;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S612.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S646;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S653;
    float3  _S694 = _S659.differential_0 + _S658.differential_0 + _S522;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S694;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S665.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S669.differential_0;
    float3  _S695 = _S678 + _S686.differential_0 + _S688;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S695;
    float3  _S696 = _S679 + _S691.differential_0 + _S693;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S696;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S697, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S698, DiffPair_float_0 * _S699, DiffPair_float_0 * _S700, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S701, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S702, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S703, DiffPair_float_0 * _S704, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S705, DiffPair_float_0 * _S706, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S707, bool _S708, bool _S709, bool _S710, bool _S711, bool _S712, FixedArray<float, 11>  * _S713, FixedArray<float, 23>  * _S714)
{
    s_bwd_prop_per_pixel_losses_0(_S697, _S698, _S699, _S700, _S701, _S702, _S703, _S704, _S705, _S706, _S707, _S708, _S709, _S710, _S711, _S712, _S713, _S714);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 11>  * weights_2, FixedArray<float, 23>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S715 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S715;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S715;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S715;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S715;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S715;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S715;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S715;
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_16, float dOut_20)
{
    float _S716 = 1.0f / ((*dpx_16).primal_0 * 52.30258560180664062f) * dOut_20;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S716;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  * raw_losses_0, FixedArray<float, 11>  * weights_3, FixedArray<float, 10>  * _S717)
{
    FixedArray<float, 10>  losses_2;
    float _S718 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*raw_losses_0)[int(0)] / _S718;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S718)));
    bool _S719;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S719 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S719 = false;
    }
    float _S720;
    if(_S719)
    {
        _S720 = (*weights_3)[int(2)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S720 = 0.0f;
    }
    losses_2[int(2)] = _S720;
    losses_2[int(3)] = ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = ((*raw_losses_0)[int(9)] + (*raw_losses_0)[int(10)]) / (F32_max(((*raw_losses_0)[int(22)]), (1.0f)));
    losses_2[int(5)] = (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    float _S721 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(6)] = (*raw_losses_0)[int(12)] / _S721;
    losses_2[int(7)] = (*raw_losses_0)[int(13)] / _S721;
    losses_2[int(8)] = (*raw_losses_0)[int(14)] / _S721;
    losses_2[int(9)] = (*raw_losses_0)[int(15)] / _S721;
    *_S717 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S722)
{
    return (F32_sqrt((_S722)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S723, float _S724)
{
    _d_log10_0(_S723, _S724);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 11>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S725 = dpraw_losses_0->primal_0;
    float _S726 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S727 = _S726 * _S726;
    float _S728 = dpraw_losses_0->primal_0[int(1)] / _S726;
    bool _S729 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S730;
    if(_S729)
    {
        _S730 = (_S725[int(3)]) != 0.0f;
    }
    else
    {
        _S730 = false;
    }
    float _S731;
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
    if(_S730)
    {
        float _S746 = _S725[int(2)] * _S725[int(3)];
        float _S747 = _S725[int(18)] * _S725[int(18)];
        float _S748 = _S725[int(6)] - _S746 / _S725[int(18)];
        float _S749 = _S725[int(2)] * _S725[int(2)];
        float _S750 = _S725[int(4)] - _S749 / _S725[int(18)];
        float _S751 = _S725[int(3)] * _S725[int(3)];
        float _S752 = _S725[int(5)] - _S751 / _S725[int(18)];
        float _S753 = _S750 * _S752 + 1.0f;
        float _S754 = s_primal_ctx_max_0(9.999999960041972e-13f, _S753);
        float _S755 = s_primal_ctx_sqrt_0(_S754);
        float _S756 = _S755 * _S755;
        float _S757 = 1.0f - _S748 / _S755;
        _S731 = (*weights_4)[int(2)];
        _S732 = _S757;
        _S733 = _S756;
        _S734 = _S748;
        _S735 = _S755;
        _S736 = _S754;
        _S737 = _S753;
        _S738 = _S750;
        _S739 = _S752;
        _S740 = _S747;
        _S741 = _S751;
        _S742 = _S725[int(3)];
        _S743 = _S749;
        _S744 = _S725[int(2)];
        _S745 = _S746;
    }
    else
    {
        _S731 = 0.0f;
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
    }
    float _S758 = s_primal_ctx_max_0(_S725[int(19)], 1.0f);
    float _S759 = _S758 * _S758;
    float _S760 = s_primal_ctx_max_0(_S725[int(20)], 1.0f);
    float _S761 = _S760 * _S760;
    float _S762 = float((I32_max((int((_S725[int(19)]) > 0.5f) + int((_S725[int(20)]) > 0.5f)), (int(1)))));
    float _S763 = _S725[int(9)] + _S725[int(10)];
    float _S764 = s_primal_ctx_max_0(_S725[int(22)], 1.0f);
    float _S765 = _S764 * _S764;
    float _S766 = s_primal_ctx_max_0(_S725[int(21)], 1.0f);
    float _S767 = _S766 * _S766;
    float _S768 = s_primal_ctx_max_0(_S725[int(16)], 1.0f);
    float _S769 = _S768 * _S768;
    float _S770 = (*_s_dOut_5)[int(9)] / _S769;
    float _S771 = _S768 * _S770;
    float _S772 = (*_s_dOut_5)[int(8)] / _S769;
    float _S773 = _S768 * _S772;
    float _S774 = (*_s_dOut_5)[int(7)] / _S769;
    float _S775 = _S768 * _S774;
    float _S776 = (*_s_dOut_5)[int(6)] / _S769;
    float _S777 = _S768 * _S776;
    float _S778 = _S725[int(15)] * - _S770 + _S725[int(14)] * - _S772 + _S725[int(13)] * - _S774 + _S725[int(12)] * - _S776;
    DiffPair_float_0 _S779;
    (&_S779)->primal_0 = _S725[int(16)];
    (&_S779)->differential_0 = 0.0f;
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = 1.0f;
    (&_S780)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S779, &_S780, _S778);
    float _S781 = (*_s_dOut_5)[int(5)] / _S767;
    float _S782 = _S725[int(11)] * - _S781;
    float _S783 = _S766 * _S781;
    DiffPair_float_0 _S784;
    (&_S784)->primal_0 = _S725[int(21)];
    (&_S784)->differential_0 = 0.0f;
    DiffPair_float_0 _S785;
    (&_S785)->primal_0 = 1.0f;
    (&_S785)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S784, &_S785, _S782);
    float _S786 = (*_s_dOut_5)[int(4)] / _S765;
    float _S787 = _S763 * - _S786;
    float _S788 = _S764 * _S786;
    DiffPair_float_0 _S789;
    (&_S789)->primal_0 = _S725[int(22)];
    (&_S789)->differential_0 = 0.0f;
    DiffPair_float_0 _S790;
    (&_S790)->primal_0 = 1.0f;
    (&_S790)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S789, &_S790, _S787);
    float _S791 = (*_s_dOut_5)[int(3)] / _S762;
    float _S792 = _S791 / _S761;
    float _S793 = _S725[int(8)] * - _S792;
    float _S794 = _S760 * _S792;
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = _S725[int(20)];
    (&_S795)->differential_0 = 0.0f;
    DiffPair_float_0 _S796;
    (&_S796)->primal_0 = 1.0f;
    (&_S796)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S795, &_S796, _S793);
    float _S797 = _S791 / _S759;
    float _S798 = _S725[int(7)] * - _S797;
    float _S799 = _S758 * _S797;
    DiffPair_float_0 _S800;
    (&_S800)->primal_0 = _S725[int(19)];
    (&_S800)->differential_0 = 0.0f;
    DiffPair_float_0 _S801;
    (&_S801)->primal_0 = 1.0f;
    (&_S801)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S800, &_S801, _S798);
    FixedArray<float, 23>  _S802;
    _S802[int(0)] = 0.0f;
    _S802[int(1)] = 0.0f;
    _S802[int(2)] = 0.0f;
    _S802[int(3)] = 0.0f;
    _S802[int(4)] = 0.0f;
    _S802[int(5)] = 0.0f;
    _S802[int(6)] = 0.0f;
    _S802[int(7)] = 0.0f;
    _S802[int(8)] = 0.0f;
    _S802[int(9)] = 0.0f;
    _S802[int(10)] = 0.0f;
    _S802[int(11)] = 0.0f;
    _S802[int(12)] = 0.0f;
    _S802[int(13)] = 0.0f;
    _S802[int(14)] = 0.0f;
    _S802[int(15)] = 0.0f;
    _S802[int(16)] = 0.0f;
    _S802[int(17)] = 0.0f;
    _S802[int(18)] = 0.0f;
    _S802[int(19)] = 0.0f;
    _S802[int(20)] = 0.0f;
    _S802[int(21)] = 0.0f;
    _S802[int(22)] = 0.0f;
    _S802[int(15)] = _S771;
    _S802[int(14)] = _S773;
    _S802[int(13)] = _S775;
    _S802[int(16)] = _S779.differential_0;
    _S802[int(12)] = _S777;
    _S802[int(21)] = _S784.differential_0;
    _S802[int(11)] = _S783;
    _S802[int(22)] = _S789.differential_0;
    _S802[int(10)] = _S788;
    _S802[int(9)] = _S788;
    _S802[int(20)] = _S795.differential_0;
    _S802[int(8)] = _S794;
    _S802[int(19)] = _S800.differential_0;
    _S802[int(7)] = _S799;
    float _S803 = _S802[int(0)];
    float _S804 = _S802[int(1)];
    float _S805 = _S802[int(2)];
    float _S806 = _S802[int(3)];
    float _S807 = _S802[int(4)];
    float _S808 = _S802[int(5)];
    float _S809 = _S802[int(6)];
    float _S810 = _S802[int(7)];
    float _S811 = _S802[int(8)];
    float _S812 = _S802[int(9)];
    float _S813 = _S802[int(10)];
    float _S814 = _S802[int(11)];
    float _S815 = _S802[int(12)];
    float _S816 = _S802[int(13)];
    float _S817 = _S802[int(14)];
    float _S818 = _S802[int(15)];
    float _S819 = _S802[int(16)];
    float _S820 = _S802[int(17)];
    float _S821 = _S802[int(18)];
    float _S822 = _S802[int(19)];
    float _S823 = _S802[int(20)];
    float _S824 = _S802[int(21)];
    float _S825 = _S802[int(22)];
    FixedArray<float, 23>  _S826;
    if(_S730)
    {
        float _S827 = _S731 * (*_s_dOut_5)[int(2)];
        DiffPair_float_0 _S828;
        (&_S828)->primal_0 = _S732;
        (&_S828)->differential_0 = 0.0f;
        DiffPair_float_0 _S829;
        (&_S829)->primal_0 = 0.0f;
        (&_S829)->differential_0 = 0.0f;
        DiffPair_float_0 _S830;
        (&_S830)->primal_0 = 2.0f;
        (&_S830)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S828, &_S829, &_S830, _S827);
        float _S831 = - _S828.differential_0 / _S733;
        float _S832 = _S734 * - _S831;
        float _S833 = _S735 * _S831;
        DiffPair_float_0 _S834;
        (&_S834)->primal_0 = _S736;
        (&_S834)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S834, _S832);
        DiffPair_float_0 _S835;
        (&_S835)->primal_0 = 9.999999960041972e-13f;
        (&_S835)->differential_0 = 0.0f;
        DiffPair_float_0 _S836;
        (&_S836)->primal_0 = _S737;
        (&_S836)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S835, &_S836, _S834.differential_0);
        float _S837 = _S738 * _S836.differential_0;
        float _S838 = _S739 * _S836.differential_0;
        float _S839 = - _S837 / _S740;
        float _S840 = _S742 * (_S725[int(18)] * _S839);
        float _S841 = - _S838 / _S740;
        float _S842 = _S744 * (_S725[int(18)] * _S841);
        float _S843 = - _S833 / _S740;
        float _S844 = _S725[int(18)] * _S843;
        float _S845 = _S840 + _S840 + _S744 * _S844;
        float _S846 = _S842 + _S842 + _S742 * _S844;
        float _S847 = _S741 * - _S839 + _S743 * - _S841 + _S745 * - _S843;
        FixedArray<float, 23>  _S848;
        _S848[int(0)] = 0.0f;
        _S848[int(1)] = 0.0f;
        _S848[int(2)] = 0.0f;
        _S848[int(3)] = 0.0f;
        _S848[int(4)] = 0.0f;
        _S848[int(5)] = 0.0f;
        _S848[int(6)] = 0.0f;
        _S848[int(7)] = 0.0f;
        _S848[int(8)] = 0.0f;
        _S848[int(9)] = 0.0f;
        _S848[int(10)] = 0.0f;
        _S848[int(11)] = 0.0f;
        _S848[int(12)] = 0.0f;
        _S848[int(13)] = 0.0f;
        _S848[int(14)] = 0.0f;
        _S848[int(15)] = 0.0f;
        _S848[int(16)] = 0.0f;
        _S848[int(17)] = 0.0f;
        _S848[int(18)] = 0.0f;
        _S848[int(19)] = 0.0f;
        _S848[int(20)] = 0.0f;
        _S848[int(21)] = 0.0f;
        _S848[int(22)] = 0.0f;
        _S848[int(5)] = _S837;
        _S848[int(4)] = _S838;
        _S848[int(3)] = _S845;
        _S848[int(2)] = _S846;
        _S848[int(6)] = _S833;
        float _S849 = _S804 + _S848[int(1)];
        float _S850 = _S805 + _S848[int(2)];
        float _S851 = _S806 + _S848[int(3)];
        float _S852 = _S807 + _S848[int(4)];
        float _S853 = _S808 + _S848[int(5)];
        float _S854 = _S809 + _S848[int(6)];
        float _S855 = _S810 + _S848[int(7)];
        float _S856 = _S811 + _S848[int(8)];
        float _S857 = _S812 + _S848[int(9)];
        float _S858 = _S813 + _S848[int(10)];
        float _S859 = _S814 + _S848[int(11)];
        float _S860 = _S815 + _S848[int(12)];
        float _S861 = _S816 + _S848[int(13)];
        float _S862 = _S817 + _S848[int(14)];
        float _S863 = _S818 + _S848[int(15)];
        float _S864 = _S819 + _S848[int(16)];
        float _S865 = _S820 + _S848[int(17)];
        float _S866 = _S821 + _S848[int(18)];
        float _S867 = _S822 + _S848[int(19)];
        float _S868 = _S823 + _S848[int(20)];
        float _S869 = _S824 + _S848[int(21)];
        float _S870 = _S825 + _S848[int(22)];
        _S826[int(0)] = _S803 + _S848[int(0)];
        _S826[int(1)] = _S849;
        _S826[int(2)] = _S850;
        _S826[int(3)] = _S851;
        _S826[int(4)] = _S852;
        _S826[int(5)] = _S853;
        _S826[int(6)] = _S854;
        _S826[int(7)] = _S855;
        _S826[int(8)] = _S856;
        _S826[int(9)] = _S857;
        _S826[int(10)] = _S858;
        _S826[int(11)] = _S859;
        _S826[int(12)] = _S860;
        _S826[int(13)] = _S861;
        _S826[int(14)] = _S862;
        _S826[int(15)] = _S863;
        _S826[int(16)] = _S864;
        _S826[int(17)] = _S865;
        _S826[int(18)] = _S866;
        _S826[int(19)] = _S867;
        _S826[int(20)] = _S868;
        _S826[int(21)] = _S869;
        _S826[int(22)] = _S870;
        _S731 = _S847;
    }
    else
    {
        _S826[int(0)] = _S803;
        _S826[int(1)] = _S804;
        _S826[int(2)] = _S805;
        _S826[int(3)] = _S806;
        _S826[int(4)] = _S807;
        _S826[int(5)] = _S808;
        _S826[int(6)] = _S809;
        _S826[int(7)] = _S810;
        _S826[int(8)] = _S811;
        _S826[int(9)] = _S812;
        _S826[int(10)] = _S813;
        _S826[int(11)] = _S814;
        _S826[int(12)] = _S815;
        _S826[int(13)] = _S816;
        _S826[int(14)] = _S817;
        _S826[int(15)] = _S818;
        _S826[int(16)] = _S819;
        _S826[int(17)] = _S820;
        _S826[int(18)] = _S821;
        _S826[int(19)] = _S822;
        _S826[int(20)] = _S823;
        _S826[int(21)] = _S824;
        _S826[int(22)] = _S825;
        _S731 = 0.0f;
    }
    if(_S729)
    {
        FixedArray<float, 23>  _S871;
        _S871[int(0)] = 0.0f;
        _S871[int(1)] = 0.0f;
        _S871[int(2)] = 0.0f;
        _S871[int(3)] = 0.0f;
        _S871[int(4)] = 0.0f;
        _S871[int(5)] = 0.0f;
        _S871[int(6)] = 0.0f;
        _S871[int(7)] = 0.0f;
        _S871[int(8)] = 0.0f;
        _S871[int(9)] = 0.0f;
        _S871[int(10)] = 0.0f;
        _S871[int(11)] = 0.0f;
        _S871[int(12)] = 0.0f;
        _S871[int(13)] = 0.0f;
        _S871[int(14)] = 0.0f;
        _S871[int(15)] = 0.0f;
        _S871[int(16)] = 0.0f;
        _S871[int(17)] = 0.0f;
        _S871[int(18)] = 0.0f;
        _S871[int(19)] = 0.0f;
        _S871[int(20)] = 0.0f;
        _S871[int(21)] = 0.0f;
        _S871[int(22)] = 0.0f;
        _S871[int(3)] = 0.0f;
        float _S872 = _S826[int(1)] + _S871[int(1)];
        float _S873 = _S826[int(2)] + _S871[int(2)];
        float _S874 = _S826[int(3)] + _S871[int(3)];
        float _S875 = _S826[int(4)] + _S871[int(4)];
        float _S876 = _S826[int(5)] + _S871[int(5)];
        float _S877 = _S826[int(6)] + _S871[int(6)];
        float _S878 = _S826[int(7)] + _S871[int(7)];
        float _S879 = _S826[int(8)] + _S871[int(8)];
        float _S880 = _S826[int(9)] + _S871[int(9)];
        float _S881 = _S826[int(10)] + _S871[int(10)];
        float _S882 = _S826[int(11)] + _S871[int(11)];
        float _S883 = _S826[int(12)] + _S871[int(12)];
        float _S884 = _S826[int(13)] + _S871[int(13)];
        float _S885 = _S826[int(14)] + _S871[int(14)];
        float _S886 = _S826[int(15)] + _S871[int(15)];
        float _S887 = _S826[int(16)] + _S871[int(16)];
        float _S888 = _S826[int(17)] + _S871[int(17)];
        float _S889 = _S826[int(18)] + _S871[int(18)];
        float _S890 = _S826[int(19)] + _S871[int(19)];
        float _S891 = _S826[int(20)] + _S871[int(20)];
        float _S892 = _S826[int(21)] + _S871[int(21)];
        float _S893 = _S826[int(22)] + _S871[int(22)];
        _S826[int(0)] = _S826[int(0)] + _S871[int(0)];
        _S826[int(1)] = _S872;
        _S826[int(2)] = _S873;
        _S826[int(3)] = _S874;
        _S826[int(4)] = _S875;
        _S826[int(5)] = _S876;
        _S826[int(6)] = _S877;
        _S826[int(7)] = _S878;
        _S826[int(8)] = _S879;
        _S826[int(9)] = _S880;
        _S826[int(10)] = _S881;
        _S826[int(11)] = _S882;
        _S826[int(12)] = _S883;
        _S826[int(13)] = _S884;
        _S826[int(14)] = _S885;
        _S826[int(15)] = _S886;
        _S826[int(16)] = _S887;
        _S826[int(17)] = _S888;
        _S826[int(18)] = _S889;
        _S826[int(19)] = _S890;
        _S826[int(20)] = _S891;
        _S826[int(21)] = _S892;
        _S826[int(22)] = _S893;
    }
    float _S894 = -10.0f * (*_s_dOut_5)[int(1)];
    DiffPair_float_0 _S895;
    (&_S895)->primal_0 = _S728;
    (&_S895)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S895, _S894);
    float _S896 = _S895.differential_0 / _S727;
    float _S897 = _S726 * _S896;
    float _S898 = (*_s_dOut_5)[int(0)] / _S727;
    float _S899 = _S726 * _S898;
    float _S900 = _S725[int(1)] * - _S896 + _S725[int(0)] * - _S898;
    DiffPair_float_0 _S901;
    (&_S901)->primal_0 = _S725[int(17)];
    (&_S901)->differential_0 = 0.0f;
    DiffPair_float_0 _S902;
    (&_S902)->primal_0 = 1.0f;
    (&_S902)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S901, &_S902, _S900);
    FixedArray<float, 23>  _S903;
    _S903[int(0)] = 0.0f;
    _S903[int(1)] = 0.0f;
    _S903[int(2)] = 0.0f;
    _S903[int(3)] = 0.0f;
    _S903[int(4)] = 0.0f;
    _S903[int(5)] = 0.0f;
    _S903[int(6)] = 0.0f;
    _S903[int(7)] = 0.0f;
    _S903[int(8)] = 0.0f;
    _S903[int(9)] = 0.0f;
    _S903[int(10)] = 0.0f;
    _S903[int(11)] = 0.0f;
    _S903[int(12)] = 0.0f;
    _S903[int(13)] = 0.0f;
    _S903[int(14)] = 0.0f;
    _S903[int(15)] = 0.0f;
    _S903[int(16)] = 0.0f;
    _S903[int(17)] = 0.0f;
    _S903[int(18)] = 0.0f;
    _S903[int(19)] = 0.0f;
    _S903[int(20)] = 0.0f;
    _S903[int(21)] = 0.0f;
    _S903[int(22)] = 0.0f;
    _S903[int(18)] = _S731;
    _S903[int(1)] = _S897;
    _S903[int(17)] = _S901.differential_0;
    _S903[int(0)] = _S899;
    FixedArray<float, 23>  _S904 = {
        _S826[int(0)] + _S903[int(0)], _S826[int(1)] + _S903[int(1)], _S826[int(2)] + _S903[int(2)], _S826[int(3)] + _S903[int(3)], _S826[int(4)] + _S903[int(4)], _S826[int(5)] + _S903[int(5)], _S826[int(6)] + _S903[int(6)], _S826[int(7)] + _S903[int(7)], _S826[int(8)] + _S903[int(8)], _S826[int(9)] + _S903[int(9)], _S826[int(10)] + _S903[int(10)], _S826[int(11)] + _S903[int(11)], _S826[int(12)] + _S903[int(12)], _S826[int(13)] + _S903[int(13)], _S826[int(14)] + _S903[int(14)], _S826[int(15)] + _S903[int(15)], _S826[int(16)] + _S903[int(16)], _S826[int(17)] + _S903[int(17)], _S826[int(18)] + _S903[int(18)], _S826[int(19)] + _S903[int(19)], _S826[int(20)] + _S903[int(20)], _S826[int(21)] + _S903[int(21)], _S826[int(22)] + _S903[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S904;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S905, FixedArray<float, 11>  * _S906, FixedArray<float, 10>  * _S907)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S905, _S906, _S907);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  * raw_losses_1, FixedArray<float, 11>  * weights_5, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 23>  * _S908)
{
    FixedArray<float, 23>  _S909 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S909;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_5, v_losses_1);
    *_S908 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_27, float3  y_9)
{
    float3  result_21;
    int i_14 = int(0);
    for(;;)
    {
        if(i_14 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_21, i_14) = (F32_min((_slang_vector_get_element(x_27, i_14)), (_slang_vector_get_element(y_9, i_14))));
        i_14 = i_14 + int(1);
    }
    return result_21;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_21)
{
    DiffPair_float_0 left_dp_3;
    (&left_dp_3)->primal_0 = (*dpx_17).primal_0.x;
    (&left_dp_3)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_5).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_3;
    (&right_dp_3)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_3)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_3, &middle_dp_0, &right_dp_3, dOut_21.x);
    float3  left_d_result_5;
    *&((&left_d_result_5)->x) = left_dp_3.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_5;
    *&((&right_d_result_5)->x) = right_dp_3.differential_0;
    DiffPair_float_0 left_dp_4;
    (&left_dp_4)->primal_0 = (*dpx_17).primal_0.y;
    (&left_dp_4)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_4;
    (&right_dp_4)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_4)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_4, &middle_dp_1, &right_dp_4, dOut_21.y);
    *&((&left_d_result_5)->y) = left_dp_4.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_5)->y) = right_dp_4.differential_0;
    DiffPair_float_0 left_dp_5;
    (&left_dp_5)->primal_0 = (*dpx_17).primal_0.z;
    (&left_dp_5)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_5;
    (&right_dp_5)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_5)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_5, &middle_dp_2, &right_dp_5, dOut_21.z);
    *&((&left_d_result_5)->z) = left_dp_5.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_5)->z) = right_dp_5.differential_0;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = left_d_result_5;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_5;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S910, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S911, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S912, float3  _S913)
{
    _d_clamp_vector_0(_S910, _S911, _S912, _S913);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S914 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S915 = make_float3 (_S914);
    float3  _S916 = make_float3 (0.0f);
    float3  _S917 = make_float3 (1.0f);
    float3  _S918 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S919;
    (&_S919)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S914) * (*dpbackground_0).primal_0;
    (&_S919)->differential_0 = _S918;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S920;
    (&_S920)->primal_0 = _S916;
    (&_S920)->differential_0 = _S918;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S921;
    (&_S921)->primal_0 = _S917;
    (&_S921)->differential_0 = _S918;
    s_bwd_prop_clamp_1(&_S919, &_S920, &_S921, _s_dOut_6);
    float3  _S922 = _S915 * _S919.differential_0;
    float3  _S923 = (*dpbackground_0).primal_0 * _S919.differential_0;
    float _S924 = - (_S923.x + _S923.y + _S923.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S922;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S924;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S919.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S925, DiffPair_float_0 * _S926, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S927, float3  _S928)
{
    s_bwd_prop_blend_background_0(_S925, _S926, _S927, _S928);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S929 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S929;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S929;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ float3  log_map_image(float3  rgb_2, float t_3)
{
    return rgb_2 * make_float3 (1.0f - t_3) + (log_0(max_0(rgb_2, make_float3 (0.0f)) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 (t_3);
}

inline __device__ void s_bwd_prop_log_map_image_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float t_4, float3  _s_dOut_7)
{
    float3  _S930 = make_float3 (1.0f - t_4);
    float3  _S931 = make_float3 (0.0f);
    float3  _S932 = make_float3 (t_4) * _s_dOut_7;
    float3  _S933 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S934;
    (&_S934)->primal_0 = s_primal_ctx_max_1((*dprgb_1).primal_0, _S931) + make_float3 (0.00392156885936856f);
    (&_S934)->differential_0 = _S933;
    s_bwd_prop_log_1(&_S934, _S932);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S935;
    (&_S935)->primal_0 = (*dprgb_1).primal_0;
    (&_S935)->differential_0 = _S933;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S936;
    (&_S936)->primal_0 = _S931;
    (&_S936)->differential_0 = _S933;
    s_bwd_prop_max_1(&_S935, &_S936, _S934.differential_0);
    float3  _S937 = _S935.differential_0 + _S930 * _s_dOut_7;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S937;
    return;
}

inline __device__ void s_bwd_log_map_image_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S938, float _S939, float3  _S940)
{
    s_bwd_prop_log_map_image_0(_S938, _S939, _S940);
    return;
}

inline __device__ float3  log_map_image_bwd(float3  rgb_3, float t_5, float3  v_out_rgb_1)
{
    float3  _S941 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S941;
    s_bwd_log_map_image_0(&p_rgb_1, t_5, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_7, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S942 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = _S942;
    bool _S943 = undistort_point_0(_S942, dist_coeffs_7, int(12), &uv_7);
    if(!_S943)
    {
        int3  _S944 = make_int3 (int(0));
        float3  _S945 = make_float3 ((float)_S944.x, (float)_S944.y, (float)_S944.z);
        return _S945;
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S946, float3  _S947)
{
    return cross_0(_S946, _S947);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_18, float _s_dOut_8)
{
    float _S948 = (*dpx_18).primal_0.x;
    float _S949 = (*dpx_18).primal_0.y;
    float _S950 = (*dpx_18).primal_0.z;
    DiffPair_float_0 _S951;
    (&_S951)->primal_0 = _S948 * _S948 + _S949 * _S949 + _S950 * _S950;
    (&_S951)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S951, _s_dOut_8);
    float _S952 = (*dpx_18).primal_0.z * _S951.differential_0;
    float _S953 = _S952 + _S952;
    float _S954 = (*dpx_18).primal_0.y * _S951.differential_0;
    float _S955 = _S954 + _S954;
    float _S956 = (*dpx_18).primal_0.x * _S951.differential_0;
    float _S957 = _S956 + _S956;
    float3  _S958 = make_float3 (0.0f);
    *&((&_S958)->z) = _S953;
    *&((&_S958)->y) = _S955;
    *&((&_S958)->x) = _S957;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S958;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S959, float _S960)
{
    s_bwd_prop_length_impl_1(_S959, _S960);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S961, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S962, float3  _S963)
{
    _d_cross_0(_S961, _S962, _S963);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S964 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S965 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S966 = s_primal_ctx_cross_0(dx_0, _S965);
    bool _S967 = (s_primal_ctx_dot_0(_S966, _S966)) != 0.0f;
    float3  _S968;
    float3  _S969;
    if(_S967)
    {
        float _S970 = length_2(_S966);
        float3  _S971 = make_float3 (_S970);
        _S968 = make_float3 (_S970 * _S970);
        _S969 = _S971;
    }
    else
    {
        _S968 = _S964;
        _S969 = _S964;
    }
    if(_S967)
    {
        float3  _S972 = _s_dOut_9 / _S968;
        float3  _S973 = _S966 * - _S972;
        float3  _S974 = _S969 * _S972;
        float _S975 = _S973.x + _S973.y + _S973.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S976;
        (&_S976)->primal_0 = _S966;
        (&_S976)->differential_0 = _S964;
        s_bwd_length_impl_1(&_S976, _S975);
        _S968 = _S974 + _S976.differential_0;
    }
    else
    {
        _S968 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S977;
    (&_S977)->primal_0 = _S966;
    (&_S977)->differential_0 = _S964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S978;
    (&_S978)->primal_0 = _S966;
    (&_S978)->differential_0 = _S964;
    s_bwd_prop_dot_0(&_S977, &_S978, 0.0f);
    float3  _S979 = _S978.differential_0 + _S977.differential_0 + _S968;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S980;
    (&_S980)->primal_0 = dx_0;
    (&_S980)->differential_0 = _S964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S981;
    (&_S981)->primal_0 = _S965;
    (&_S981)->differential_0 = _S964;
    s_bwd_prop_cross_0(&_S980, &_S981, _S979);
    float3  s_diff_dy_T_0 = - _S981.differential_0;
    float3  _S982 = - s_diff_dy_T_0;
    float3  _S983 = - _S980.differential_0;
    FixedArray<float3 , 4>  _S984;
    _S984[int(0)] = _S964;
    _S984[int(1)] = _S964;
    _S984[int(2)] = _S964;
    _S984[int(3)] = _S964;
    _S984[int(2)] = _S982;
    _S984[int(3)] = s_diff_dy_T_0;
    _S984[int(0)] = _S983;
    _S984[int(1)] = _S980.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S984;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S985, float3  _S986)
{
    s_bwd_prop_points_to_normal_0(_S985, _S986);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  * points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S987 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = *points_1;
    (&dp_points_0)->differential_0 = _S987;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S988 = float2 {intrins_1.z, intrins_1.w};
    float2  _S989 = float2 {intrins_1.x, intrins_1.y};
    float2  _S990 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S988) / _S989;
    float2  uv_8 = _S990;
    bool _S991 = undistort_point_0(_S990, dist_coeffs_8, int(12), &uv_8);
    if(!_S991)
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
    float2  _S992 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S988) / _S989;
    float2  uv_9 = _S992;
    bool _S993 = undistort_point_0(_S992, dist_coeffs_8, int(12), &uv_9);
    if(!_S993)
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
    float2  _S994 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S988) / _S989;
    float2  uv_10 = _S994;
    bool _S995 = undistort_point_0(_S994, dist_coeffs_8, int(12), &uv_10);
    if(!_S995)
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
    float2  _S996 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S988) / _S989;
    float2  uv_11 = _S996;
    bool _S997 = undistort_point_0(_S996, dist_coeffs_8, int(12), &uv_11);
    if(!_S997)
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
    float2  _S998;
    bool _S999;
    float2  _S1000;
    bool _S1001;
    float2  _S1002;
    bool _S1003;
    float2  _S1004;
    bool _S1005;
};

inline __device__ float s_primal_ctx_sin_0(float _S1006)
{
    return (F32_sin((_S1006)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1007)
{
    return (F32_cos((_S1007)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1008 = make_float2 (0.0f);
    _s_diff_ctx_0->_S998 = _S1008;
    _s_diff_ctx_0->_S999 = false;
    _s_diff_ctx_0->_S1000 = _S1008;
    _s_diff_ctx_0->_S1001 = false;
    _s_diff_ctx_0->_S1002 = _S1008;
    _s_diff_ctx_0->_S1003 = false;
    _s_diff_ctx_0->_S1004 = _S1008;
    _s_diff_ctx_0->_S1005 = false;
    _s_diff_ctx_0->_S1000 = _S1008;
    _s_diff_ctx_0->_S1001 = false;
    _s_diff_ctx_0->_S1002 = _S1008;
    _s_diff_ctx_0->_S1003 = false;
    _s_diff_ctx_0->_S1004 = _S1008;
    _s_diff_ctx_0->_S1005 = false;
    float3  _S1009 = make_float3 (0.0f);
    float2  _S1010 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1011 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1012 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1010) / _S1011;
    float2  _S1013 = _S1012;
    bool _S1014 = undistort_point_0(_S1012, dist_coeffs_9, int(12), &_S1013);
    _s_diff_ctx_0->_S998 = _S1013;
    _s_diff_ctx_0->_S999 = _S1014;
    float2  uv_12 = _S1013;
    bool _S1015 = !_S1014;
    float3  normal_4;
    if(_S1015)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1016 = !_S1015;
    int _S1017;
    FixedArray<float3 , 4>  points_3;
    if(_S1016)
    {
        float3  raydir_18;
        if(is_fisheye_6)
        {
            float _S1018 = length_1(uv_12);
            float3  raydir_19 = make_float3 ((uv_12 / make_float2 (s_primal_ctx_max_0(_S1018, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1018))).x, (uv_12 / make_float2 (s_primal_ctx_max_0(_S1018, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1018))).y, s_primal_ctx_cos_0(_S1018));
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
        float3  _S1019 = make_float3 (dpdepths_0.x) * raydir_18;
        float2  _S1020 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1010) / _S1011;
        float2  _S1021 = _S1020;
        bool _S1022 = undistort_point_0(_S1020, dist_coeffs_9, int(12), &_S1021);
        _s_diff_ctx_0->_S1000 = _S1021;
        _s_diff_ctx_0->_S1001 = _S1022;
        float2  uv_13 = _S1021;
        bool _S1023 = !_S1022;
        if(_S1023)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1024 = !_S1023;
        if(_S1024)
        {
            if(is_fisheye_6)
            {
                float _S1025 = length_1(uv_13);
                float3  raydir_21 = make_float3 ((uv_13 / make_float2 (s_primal_ctx_max_0(_S1025, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1025))).x, (uv_13 / make_float2 (s_primal_ctx_max_0(_S1025, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1025))).y, s_primal_ctx_cos_0(_S1025));
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
            float3  _S1026 = make_float3 (dpdepths_0.y) * raydir_18;
            _S1017 = int(2);
            points_3[int(0)] = _S1019;
            points_3[int(1)] = _S1026;
            points_3[int(2)] = _S1009;
            points_3[int(3)] = _S1009;
        }
        else
        {
            _S1017 = int(0);
            points_3[int(0)] = _S1019;
            points_3[int(1)] = _S1009;
            points_3[int(2)] = _S1009;
            points_3[int(3)] = _S1009;
        }
        bool _runFlag_0;
        if(_S1017 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1016;
            _S1017 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1027 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1010) / _S1011;
            float2  _S1028 = _S1027;
            bool _S1029 = undistort_point_0(_S1027, dist_coeffs_9, int(12), &_S1028);
            _s_diff_ctx_0->_S1002 = _S1028;
            _s_diff_ctx_0->_S1003 = _S1029;
            float2  uv_14 = _S1028;
            if(!_S1029)
            {
                float3  _S1030 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1017 = int(0);
                normal_4 = _S1030;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1031 = length_1(uv_14);
                    float3  raydir_23 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S1031, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1031))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S1031, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1031))).y, s_primal_ctx_cos_0(_S1031));
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
                float2  _S1032 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1010) / _S1011;
                float2  _S1033 = _S1032;
                bool _S1034 = undistort_point_0(_S1032, dist_coeffs_9, int(12), &_S1033);
                _s_diff_ctx_0->_S1004 = _S1033;
                _s_diff_ctx_0->_S1005 = _S1034;
                float2  uv_15 = _S1033;
                bool _S1035 = !_S1034;
                if(_S1035)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1036 = !_S1035;
                int _S1037;
                if(_S1036)
                {
                    if(is_fisheye_6)
                    {
                        float _S1038 = length_1(uv_15);
                        float3  raydir_25 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S1038, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1038))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S1038, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1038))).y, s_primal_ctx_cos_0(_S1038));
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
                    _S1037 = int(2);
                }
                else
                {
                    _S1037 = int(0);
                }
                if(_S1037 != int(2))
                {
                    _runFlag_0 = false;
                    _S1017 = _S1037;
                }
                if(_runFlag_0)
                {
                    _S1017 = int(1);
                }
            }
        }
    }
    else
    {
        _S1017 = int(0);
        points_3[int(0)] = _S1009;
        points_3[int(1)] = _S1009;
        points_3[int(2)] = _S1009;
        points_3[int(3)] = _S1009;
    }
    if(!(_S1017 != int(1)))
    {
        float3  _S1039 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1039, _S1039)) != 0.0f)
        {
            normal_4 = _S1039 / make_float3 (length_2(_S1039));
        }
        else
        {
            normal_4 = _S1039;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1040 = *dpdepths_1;
    float3  _S1041 = make_float3 (0.0f);
    bool _S1042 = !!_s_diff_ctx_1->_S999;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    float3  raydir_30;
    int _S1043;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1044;
    if(_S1042)
    {
        if(is_fisheye_7)
        {
            float _S1045 = length_1(_s_diff_ctx_1->_S998);
            float3  raydir_31 = make_float3 ((_s_diff_ctx_1->_S998 / make_float2 (s_primal_ctx_max_0(_S1045, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1045))).x, (_s_diff_ctx_1->_S998 / make_float2 (s_primal_ctx_max_0(_S1045, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1045))).y, s_primal_ctx_cos_0(_S1045));
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
            float3  raydir_32 = make_float3 (_s_diff_ctx_1->_S998.x, _s_diff_ctx_1->_S998.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_27 = normalize_0(raydir_32);
            }
            else
            {
                raydir_27 = raydir_32;
            }
        }
        float3  _S1046 = make_float3 (_S1040.primal_0.x) * raydir_27;
        bool _S1047 = !!_s_diff_ctx_1->_S1001;
        if(_S1047)
        {
            if(is_fisheye_7)
            {
                float _S1048 = length_1(_s_diff_ctx_1->_S1000);
                float3  raydir_33 = make_float3 ((_s_diff_ctx_1->_S1000 / make_float2 (s_primal_ctx_max_0(_S1048, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1048))).x, (_s_diff_ctx_1->_S1000 / make_float2 (s_primal_ctx_max_0(_S1048, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1048))).y, s_primal_ctx_cos_0(_S1048));
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
                float3  raydir_34 = make_float3 (_s_diff_ctx_1->_S1000.x, _s_diff_ctx_1->_S1000.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_28 = normalize_0(raydir_34);
                }
                else
                {
                    raydir_28 = raydir_34;
                }
            }
            float3  _S1049 = make_float3 (_S1040.primal_0.y) * raydir_28;
            _S1043 = int(2);
            points_4[int(0)] = _S1046;
            points_4[int(1)] = _S1049;
            points_4[int(2)] = _S1041;
            points_4[int(3)] = _S1041;
        }
        else
        {
            _S1043 = int(0);
            points_4[int(0)] = _S1046;
            points_4[int(1)] = _S1041;
            points_4[int(2)] = _S1041;
            points_4[int(3)] = _S1041;
            raydir_28 = _S1041;
        }
        if(_S1043 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1042;
            _S1043 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_1->_S1003)
            {
                _runFlag_2 = false;
                _S1043 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1050 = length_1(_s_diff_ctx_1->_S1002);
                    float3  raydir_35 = make_float3 ((_s_diff_ctx_1->_S1002 / make_float2 (s_primal_ctx_max_0(_S1050, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1050))).x, (_s_diff_ctx_1->_S1002 / make_float2 (s_primal_ctx_max_0(_S1050, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1050))).y, s_primal_ctx_cos_0(_S1050));
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
                    float3  raydir_36 = make_float3 (_s_diff_ctx_1->_S1002.x, _s_diff_ctx_1->_S1002.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_29 = normalize_0(raydir_36);
                    }
                    else
                    {
                        raydir_29 = raydir_36;
                    }
                }
                points_4[int(2)] = make_float3 (_S1040.primal_0.z) * raydir_29;
                bool _S1051 = !!_s_diff_ctx_1->_S1005;
                int _S1052;
                if(_S1051)
                {
                    if(is_fisheye_7)
                    {
                        float _S1053 = length_1(_s_diff_ctx_1->_S1004);
                        float3  raydir_37 = make_float3 ((_s_diff_ctx_1->_S1004 / make_float2 (s_primal_ctx_max_0(_S1053, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1053))).x, (_s_diff_ctx_1->_S1004 / make_float2 (s_primal_ctx_max_0(_S1053, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1053))).y, s_primal_ctx_cos_0(_S1053));
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
                        float3  raydir_38 = make_float3 (_s_diff_ctx_1->_S1004.x, _s_diff_ctx_1->_S1004.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_30 = normalize_0(raydir_38);
                        }
                        else
                        {
                            raydir_30 = raydir_38;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1040.primal_0.w) * raydir_30;
                    _S1052 = int(2);
                }
                else
                {
                    _S1052 = int(0);
                    raydir_30 = _S1041;
                }
                if(_S1052 != int(2))
                {
                    _runFlag_3 = false;
                    _S1043 = _S1052;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1043 = int(1);
                }
                float3  _S1054 = raydir_29;
                _runFlag_3 = _S1051;
                raydir_29 = raydir_30;
                raydir_30 = _S1054;
            }
            else
            {
                _runFlag_3 = false;
                raydir_29 = _S1041;
                raydir_30 = _S1041;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_29 = _S1041;
            raydir_30 = _S1041;
        }
        float3  _S1055 = raydir_27;
        float3  _S1056 = raydir_28;
        raydir_27 = raydir_29;
        raydir_28 = raydir_30;
        _S1044 = _S1047;
        raydir_29 = _S1056;
        raydir_30 = _S1055;
    }
    else
    {
        _S1043 = int(0);
        points_4[int(0)] = _S1041;
        points_4[int(1)] = _S1041;
        points_4[int(2)] = _S1041;
        points_4[int(3)] = _S1041;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_27 = _S1041;
        raydir_28 = _S1041;
        _S1044 = false;
        raydir_29 = _S1041;
        raydir_30 = _S1041;
    }
    bool _S1057 = !(_S1043 != int(1));
    float3  _S1058;
    float3  _S1059;
    float3  _S1060;
    float3  _S1061;
    float3  _S1062;
    bool _S1063;
    if(_S1057)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1064 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1065 = s_primal_ctx_cross_0(dx_1, _S1064);
        bool _S1066 = (s_primal_ctx_dot_0(_S1065, _S1065)) != 0.0f;
        if(_S1066)
        {
            float _S1067 = length_2(_S1065);
            float3  _S1068 = make_float3 (_S1067);
            _S1058 = make_float3 (_S1067 * _S1067);
            _S1059 = _S1068;
        }
        else
        {
            _S1058 = _S1041;
            _S1059 = _S1041;
        }
        float3  _S1069 = _S1059;
        _S1063 = _S1066;
        _S1059 = _S1065;
        _S1060 = _S1069;
        _S1061 = dx_1;
        _S1062 = _S1064;
    }
    else
    {
        _S1063 = false;
        _S1058 = _S1041;
        _S1059 = _S1041;
        _S1060 = _S1041;
        _S1061 = _S1041;
        _S1062 = _S1041;
    }
    float4  _S1070 = make_float4 (0.0f);
    if(_S1057)
    {
        if(_S1063)
        {
            float3  _S1071 = _s_dOut_10 / _S1058;
            float3  _S1072 = _S1059 * - _S1071;
            float3  _S1073 = _S1060 * _S1071;
            float _S1074 = _S1072.x + _S1072.y + _S1072.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1075;
            (&_S1075)->primal_0 = _S1059;
            (&_S1075)->differential_0 = _S1041;
            s_bwd_length_impl_1(&_S1075, _S1074);
            _S1058 = _S1073 + _S1075.differential_0;
        }
        else
        {
            _S1058 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1076;
        (&_S1076)->primal_0 = _S1059;
        (&_S1076)->differential_0 = _S1041;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1077;
        (&_S1077)->primal_0 = _S1059;
        (&_S1077)->differential_0 = _S1041;
        s_bwd_prop_dot_0(&_S1076, &_S1077, 0.0f);
        float3  _S1078 = _S1077.differential_0 + _S1076.differential_0 + _S1058;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1079;
        (&_S1079)->primal_0 = _S1061;
        (&_S1079)->differential_0 = _S1041;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1080;
        (&_S1080)->primal_0 = _S1062;
        (&_S1080)->differential_0 = _S1041;
        s_bwd_prop_cross_0(&_S1079, &_S1080, _S1078);
        float3  s_diff_dy_T_1 = - _S1080.differential_0;
        float3  _S1081 = - s_diff_dy_T_1;
        float3  _S1082 = - _S1079.differential_0;
        FixedArray<float3 , 4>  _S1083;
        _S1083[int(0)] = _S1041;
        _S1083[int(1)] = _S1041;
        _S1083[int(2)] = _S1041;
        _S1083[int(3)] = _S1041;
        _S1083[int(2)] = _S1081;
        _S1083[int(3)] = s_diff_dy_T_1;
        _S1083[int(0)] = _S1082;
        _S1083[int(1)] = _S1079.differential_0;
        points_4[int(0)] = _S1083[int(0)];
        points_4[int(1)] = _S1083[int(1)];
        points_4[int(2)] = _S1083[int(2)];
        points_4[int(3)] = _S1083[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1041;
        points_4[int(1)] = _S1041;
        points_4[int(2)] = _S1041;
        points_4[int(3)] = _S1041;
    }
    float4  _S1084;
    if(_S1042)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1085 = points_4;
                FixedArray<float3 , 4>  _S1086 = points_4;
                FixedArray<float3 , 4>  _S1087 = points_4;
                FixedArray<float3 , 4>  _S1088 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1089 = raydir_27 * _S1088[int(3)];
                    float _S1090 = _S1089.x + _S1089.y + _S1089.z;
                    float4  _S1091 = _S1070;
                    *&((&_S1091)->w) = _S1090;
                    points_4[int(0)] = _S1085[int(0)];
                    points_4[int(1)] = _S1086[int(1)];
                    points_4[int(2)] = _S1087[int(2)];
                    points_4[int(3)] = _S1041;
                    _S1084 = _S1091;
                }
                else
                {
                    points_4[int(0)] = _S1085[int(0)];
                    points_4[int(1)] = _S1086[int(1)];
                    points_4[int(2)] = _S1087[int(2)];
                    points_4[int(3)] = _S1088[int(3)];
                    _S1084 = _S1070;
                }
                float3  _S1092 = raydir_28 * points_4[int(2)];
                float _S1093 = _S1092.x + _S1092.y + _S1092.z;
                FixedArray<float3 , 4>  _S1094 = points_4;
                FixedArray<float3 , 4>  _S1095 = points_4;
                float4  _S1096 = _S1070;
                *&((&_S1096)->z) = _S1093;
                float4  _S1097 = _S1084 + _S1096;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1094[int(1)];
                points_4[int(2)] = _S1041;
                points_4[int(3)] = _S1095[int(3)];
                _S1084 = _S1097;
            }
            else
            {
                FixedArray<float3 , 4>  _S1098 = points_4;
                FixedArray<float3 , 4>  _S1099 = points_4;
                FixedArray<float3 , 4>  _S1100 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1098[int(1)];
                points_4[int(2)] = _S1099[int(2)];
                points_4[int(3)] = _S1100[int(3)];
                _S1084 = _S1070;
            }
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
            _S1084 = _S1070;
        }
        if(_S1044)
        {
            FixedArray<float3 , 4>  _S1104 = points_4;
            float3  _S1105 = raydir_29 * points_4[int(1)];
            float _S1106 = _S1105.x + _S1105.y + _S1105.z;
            float4  _S1107 = _S1070;
            *&((&_S1107)->y) = _S1106;
            float4  _S1108 = _S1084 + _S1107;
            points_4[int(0)] = _S1041;
            points_4[int(1)] = _S1041;
            points_4[int(2)] = _S1041;
            points_4[int(3)] = _S1041;
            raydir_27 = _S1104[int(0)];
            _S1084 = _S1108;
        }
        else
        {
            FixedArray<float3 , 4>  _S1109 = points_4;
            FixedArray<float3 , 4>  _S1110 = points_4;
            FixedArray<float3 , 4>  _S1111 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1109[int(1)];
            points_4[int(2)] = _S1110[int(2)];
            points_4[int(3)] = _S1111[int(3)];
            raydir_27 = _S1041;
        }
        float3  _S1112 = raydir_30 * (points_4[int(0)] + raydir_27);
        float _S1113 = _S1112.x + _S1112.y + _S1112.z;
        float4  _S1114 = _S1070;
        *&((&_S1114)->x) = _S1113;
        _S1084 = _S1084 + _S1114;
    }
    else
    {
        _S1084 = _S1070;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1084;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1115, float4  _S1116, FixedArray<float, 10>  * _S1117, bool _S1118, bool _S1119, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1120, float3  _S1121)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1122;
    float3  _S1123 = s_primal_ctx_depth_to_normal_0(_S1115, _S1116, _S1117, _S1118, _S1119, (*_S1120).primal_0, &_S1122);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1124 = _S1122;
    s_bwd_prop_depth_to_normal_0(_S1115, _S1116, _S1117, _S1118, _S1119, _S1120, _S1121, &_S1124);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1125 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1125;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, dist_coeffs_11, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9)
{
    float2  _S1126 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_16 = _S1126;
    bool _S1127 = undistort_point_0(_S1126, dist_coeffs_12, int(12), &uv_16);
    if(!_S1127)
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

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_19, float dOut_22)
{
    float _S1128 = (F32_exp2(((*dpx_19).primal_0))) * 50.693145751953125f * dOut_22;
    dpx_19->primal_0 = (*dpx_19).primal_0;
    dpx_19->differential_0 = _S1128;
    return;
}

inline __device__ void _d_fma_0(DiffPair_float_0 * dpx_20, DiffPair_float_0 * dpy_6, DiffPair_float_0 * dpz_1, float dOut_23)
{
    DiffPair_float_0 _S1129 = *dpx_20;
    float _S1130 = (*dpy_6).primal_0 * dOut_23;
    dpx_20->primal_0 = (*dpx_20).primal_0;
    dpx_20->differential_0 = _S1130;
    float _S1131 = _S1129.primal_0 * dOut_23;
    dpy_6->primal_0 = (*dpy_6).primal_0;
    dpy_6->differential_0 = _S1131;
    dpz_1->primal_0 = (*dpz_1).primal_0;
    dpz_1->differential_0 = dOut_23;
    return;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_21, DiffPair_float_0 * dpy_7, float dOut_24)
{
    if(((*dpx_21).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_21->primal_0 = (*dpx_21).primal_0;
        dpx_21->differential_0 = 0.0f;
        dpy_7->primal_0 = (*dpy_7).primal_0;
        dpy_7->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_21).primal_0), ((*dpy_7).primal_0)));
        DiffPair_float_0 _S1132 = *dpx_21;
        float _S1133 = val_0 * (*dpy_7).primal_0 / (*dpx_21).primal_0 * dOut_24;
        dpx_21->primal_0 = (*dpx_21).primal_0;
        dpx_21->differential_0 = _S1133;
        float _S1134 = val_0 * (F32_log((_S1132.primal_0))) * dOut_24;
        dpy_7->primal_0 = (*dpy_7).primal_0;
        dpy_7->differential_0 = _S1134;
    }
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
    PPISPParams_0 _S1135 = p_0;
    float max_res_0 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1136 = (pix_coord_0.x - image_center_0.x) / max_res_0;
    float _S1137 = (pix_coord_0.y - image_center_0.y) / max_res_0;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_0))));
    float dx_2 = _S1136 - p_0.vignette_params_0[int(0)].cx_0;
    float dy_0 = _S1137 - p_0.vignette_params_0[int(0)].cy_0;
    float r2_4 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_4 * r2_4;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_0[int(0)].alpha2_0 * (r4_0 * r2_4) + p_0.vignette_params_0[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_0[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1136 - p_0.vignette_params_0[int(1)].cx_0;
    float dy_1 = _S1137 - p_0.vignette_params_0[int(1)].cy_0;
    float r2_5 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_5 * r2_5;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_0[int(1)].alpha2_0 * (r4_1 * r2_5) + p_0.vignette_params_0[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_0[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1136 - p_0.vignette_params_0[int(2)].cx_0;
    float dy_2 = _S1137 - p_0.vignette_params_0[int(2)].cy_0;
    float r2_6 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_6 * r2_6;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_0[int(2)].alpha2_0 * (r4_2 * r2_6) + p_0.vignette_params_0[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_0[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S1138 = rgb_out_0;
    float2  bd_0 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_0.b_0);
    float2  rd_0 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_0.r_0);
    float2  gd_0 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_0.g_0);
    float2  nd_0 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_0.n_0);
    float _S1139 = 0.3333333432674408f + nd_0.x;
    float _S1140 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_3 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1140, 1.0f, 0.0f, - _S1139, - _S1140, _S1139, 0.0f), T_0);
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
    float _S1141 = _S1138.x;
    float _S1142 = _S1138.y;
    float intensity_0 = _S1141 + _S1142 + _S1138.z;
    float3  rgi_out_0 = mul_0(H_1, make_float3 (_S1141, _S1142, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1143 = rgi_out_1.x;
    float _S1144 = rgi_out_1.y;
    float3  _S1145 = clamp_1(make_float3 (_S1143, _S1144, rgi_out_1.z - _S1143 - _S1144), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1146 = _S1145.x;
    float _S1147 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(0)].toe_0))))));
    float _S1148 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(0)].shoulder_0))))));
    float _S1149 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(0)].gamma_0))))));
    float _S1150 = 1.0f / (1.0f + (F32_exp((- _S1135.crf_params_0[int(0)].center_0))));
    float a_1 = _S1148 * _S1150 / (F32_fma((_S1148 - _S1147), (_S1150), (_S1147)));
    float b_2 = 1.0f - a_1;
    float y_10;
    if(_S1146 <= _S1150)
    {
        y_10 = a_1 * (F32_pow((_S1146 / _S1150), (_S1147)));
    }
    else
    {
        y_10 = 1.0f - b_2 * (F32_pow(((1.0f - _S1146) / (1.0f - _S1150)), (_S1148)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1149)));
    float _S1151 = _S1145.y;
    float _S1152 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(1)].toe_0))))));
    float _S1153 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(1)].shoulder_0))))));
    float _S1154 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(1)].gamma_0))))));
    float _S1155 = 1.0f / (1.0f + (F32_exp((- _S1135.crf_params_0[int(1)].center_0))));
    float a_2 = _S1153 * _S1155 / (F32_fma((_S1153 - _S1152), (_S1155), (_S1152)));
    float b_3 = 1.0f - a_2;
    if(_S1151 <= _S1155)
    {
        y_10 = a_2 * (F32_pow((_S1151 / _S1155), (_S1152)));
    }
    else
    {
        y_10 = 1.0f - b_3 * (F32_pow(((1.0f - _S1151) / (1.0f - _S1155)), (_S1153)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1154)));
    float _S1156 = _S1145.z;
    float _S1157 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(2)].toe_0))))));
    float _S1158 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(2)].shoulder_0))))));
    float _S1159 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1135.crf_params_0[int(2)].gamma_0))))));
    float _S1160 = 1.0f / (1.0f + (F32_exp((- _S1135.crf_params_0[int(2)].center_0))));
    float a_3 = _S1158 * _S1160 / (F32_fma((_S1158 - _S1157), (_S1160), (_S1157)));
    float b_4 = 1.0f - a_3;
    if(_S1156 <= _S1160)
    {
        y_10 = a_3 * (F32_pow((_S1156 / _S1160), (_S1157)));
    }
    else
    {
        y_10 = 1.0f - b_4 * (F32_pow(((1.0f - _S1156) / (1.0f - _S1160)), (_S1158)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_10)))), (_S1159)));
    return rgb_out_1;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1161)
{
    return (F32_exp2((_S1161)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1162, float2  _S1163)
{
    return mul_1(_S1162, _S1163);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1164, Matrix<float, 3, 3>  _S1165)
{
    return mul_3(_S1164, _S1165);
}

inline __device__ float s_primal_ctx_abs_0(float _S1166)
{
    return (F32_abs((_S1166)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1167, float3  _S1168)
{
    return mul_0(_S1167, _S1168);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1169, float3  _S1170, float3  _S1171)
{
    return clamp_1(_S1169, _S1170, _S1171);
}

inline __device__ float s_primal_ctx_fma_0(float _S1172, float _S1173, float _S1174)
{
    return (F32_fma((_S1172), (_S1173), (_S1174)));
}

inline __device__ float s_primal_ctx_pow_0(float _S1175, float _S1176)
{
    return (F32_pow((_S1175), (_S1176)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1177, DiffPair_float_0 * _S1178, float _S1179)
{
    _d_pow_0(_S1177, _S1178, _S1179);
    return;
}

inline __device__ void s_bwd_prop_fma_0(DiffPair_float_0 * _S1180, DiffPair_float_0 * _S1181, DiffPair_float_0 * _S1182, float _S1183)
{
    _d_fma_0(_S1180, _S1181, _S1182, _S1183);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1184, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1185, float3  _S1186)
{
    _d_mul_0(_S1184, _S1185, _S1186);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1187, float _S1188)
{
    _d_abs_0(_S1187, _S1188);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1189, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1190, Matrix<float, 3, 3>  _S1191)
{
    mul_2(_S1189, _S1190, _S1191);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1192, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1193, float2  _S1194)
{
    _d_mul_1(_S1192, _S1193, _S1194);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1195, float _S1196)
{
    _d_exp2_0(_S1195, _S1196);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1197 = *dprgb_in_0;
    float3  _S1198 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1199 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1200 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1201 = {
        _S1200, _S1200, _S1200
    };
    float2  _S1202 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1203 = { _S1202, _S1202, _S1202, _S1202 };
    CRFPPISPChannelParams_0 _S1204 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1205 = {
        _S1204, _S1204, _S1204
    };
    PPISPParams_0 _S1206;
    (&_S1206)->exposure_0 = dpparams_0->primal_0[int(0)];
    (&_S1206)->vignette_params_0 = _S1201;
    (&_S1206)->color_params_0 = _S1203;
    (&_S1206)->crf_params_0 = _S1205;
    (&(&_S1206)->vignette_params_0[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1206)->vignette_params_0[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1207 = dpparams_0->primal_0[int(3)];
    (&(&_S1206)->vignette_params_0[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1208 = dpparams_0->primal_0[int(4)];
    (&(&_S1206)->vignette_params_0[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1209 = dpparams_0->primal_0[int(5)];
    (&(&_S1206)->vignette_params_0[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1206)->vignette_params_0[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1206)->vignette_params_0[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1210 = dpparams_0->primal_0[int(8)];
    (&(&_S1206)->vignette_params_0[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1211 = dpparams_0->primal_0[int(9)];
    (&(&_S1206)->vignette_params_0[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1212 = dpparams_0->primal_0[int(10)];
    (&(&_S1206)->vignette_params_0[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1206)->vignette_params_0[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1206)->vignette_params_0[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1213 = dpparams_0->primal_0[int(13)];
    (&(&_S1206)->vignette_params_0[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1214 = dpparams_0->primal_0[int(14)];
    (&(&_S1206)->vignette_params_0[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1215 = dpparams_0->primal_0[int(15)];
    (&(&_S1206)->vignette_params_0[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1206)->color_params_0)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1206)->color_params_0)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1206)->color_params_0)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1206)->color_params_0)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1206)->color_params_0)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1206)->color_params_0)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1206)->color_params_0)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1206)->color_params_0)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1216 = dpparams_0->primal_0[int(24)];
    (&(&_S1206)->crf_params_0[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1217 = dpparams_0->primal_0[int(25)];
    (&(&_S1206)->crf_params_0[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1218 = dpparams_0->primal_0[int(26)];
    (&(&_S1206)->crf_params_0[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1219 = dpparams_0->primal_0[int(27)];
    (&(&_S1206)->crf_params_0[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1220 = dpparams_0->primal_0[int(28)];
    (&(&_S1206)->crf_params_0[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1221 = dpparams_0->primal_0[int(29)];
    (&(&_S1206)->crf_params_0[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1222 = dpparams_0->primal_0[int(30)];
    (&(&_S1206)->crf_params_0[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1223 = dpparams_0->primal_0[int(31)];
    (&(&_S1206)->crf_params_0[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1224 = dpparams_0->primal_0[int(32)];
    (&(&_S1206)->crf_params_0[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1225 = dpparams_0->primal_0[int(33)];
    (&(&_S1206)->crf_params_0[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1226 = dpparams_0->primal_0[int(34)];
    (&(&_S1206)->crf_params_0[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1227 = dpparams_0->primal_0[int(35)];
    (&(&_S1206)->crf_params_0[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1228 = _S1206;
    float _S1229 = s_primal_ctx_exp2_0(_S1206.exposure_0);
    float3  _S1230 = make_float3 (_S1229);
    float3  rgb_out_2 = (*dprgb_in_0).primal_0 * make_float3 (_S1229);
    float _S1231 = s_primal_ctx_max_0(img_size_1.x, img_size_1.y);
    float _S1232 = (pix_coord_1.x - image_center_1.x) / _S1231;
    float _S1233 = (pix_coord_1.y - image_center_1.y) / _S1231;
    float dx_5 = _S1232 - dpparams_0->primal_0[int(1)];
    float dy_3 = _S1233 - dpparams_0->primal_0[int(2)];
    float r2_8 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_8 * r2_8;
    float r6_0 = r4_3 * r2_8;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_3 + dpparams_0->primal_0[int(3)] * r2_8 + 1.0f;
    float _S1234 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1235 = rgb_out_2.x * _S1234;
    float3  _S1236 = rgb_out_2;
    *&((&_S1236)->x) = _S1235;
    float dx_6 = _S1232 - dpparams_0->primal_0[int(6)];
    float dy_4 = _S1233 - dpparams_0->primal_0[int(7)];
    float r2_9 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_9 * r2_9;
    float r6_1 = r4_4 * r2_9;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_4 + dpparams_0->primal_0[int(8)] * r2_9 + 1.0f;
    float _S1237 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1236)->y) = rgb_out_2.y * _S1237;
    float dx_7 = _S1232 - dpparams_0->primal_0[int(11)];
    float dy_5 = _S1233 - dpparams_0->primal_0[int(12)];
    float r2_10 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_10 * r2_10;
    float r6_2 = r4_5 * r2_10;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_5 + dpparams_0->primal_0[int(13)] * r2_10 + 1.0f;
    float _S1238 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1236)->z) = rgb_out_2.z * _S1238;
    PPISPParams_0 _S1239 = _S1206;
    float2  _S1240 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1206.color_params_0.b_0);
    float2  _S1241 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1206.color_params_0.r_0);
    float2  _S1242 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1206.color_params_0.g_0);
    float2  _S1243 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1206.color_params_0.n_0);
    float _S1244 = 0.3333333432674408f + _S1243.x;
    float _S1245 = 0.3333333432674408f + _S1243.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (_S1240.x, 1.0f + _S1241.x, _S1242.x, _S1240.y, _S1241.y, 1.0f + _S1242.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1245, 1.0f, 0.0f, - _S1244, - _S1245, _S1244, 0.0f);
    Matrix<float, 3, 3>  _S1246 = s_primal_ctx_mul_1(skew_0, T_1);
    float3  r0_1 = make_float3 (_S1246.rows[int(0)].x, _S1246.rows[int(0)].y, _S1246.rows[int(0)].z);
    float3  r1_1 = make_float3 (_S1246.rows[int(1)].x, _S1246.rows[int(1)].y, _S1246.rows[int(1)].z);
    float3  r2_11 = make_float3 (_S1246.rows[int(2)].x, _S1246.rows[int(2)].y, _S1246.rows[int(2)].z);
    float3  _S1247 = s_primal_ctx_cross_0(r0_1, r1_1);
    bool _S1248 = (s_primal_ctx_dot_0(_S1247, _S1247)) < 9.99999968265522539e-21f;
    float3  lambda_v_3;
    float3  _S1249;
    bool _S1250;
    if(_S1248)
    {
        float3  _S1251 = s_primal_ctx_cross_0(r0_1, r2_11);
        bool _S1252 = (s_primal_ctx_dot_0(_S1251, _S1251)) < 9.99999968265522539e-21f;
        if(_S1252)
        {
            lambda_v_3 = s_primal_ctx_cross_0(r1_1, r2_11);
        }
        else
        {
            lambda_v_3 = _S1251;
        }
        _S1250 = _S1252;
        _S1249 = _S1251;
    }
    else
    {
        lambda_v_3 = _S1247;
        _S1250 = false;
        _S1249 = _S1198;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_3.x, 0.0f, 0.0f, 0.0f, lambda_v_3.y, 0.0f, 0.0f, 0.0f, lambda_v_3.z);
    Matrix<float, 3, 3>  _S1253 = s_primal_ctx_mul_1(T_1, D_0);
    Matrix<float, 3, 3>  _S1254 = s_primal_ctx_mul_1(_S1253, S_inv_0);
    bool _S1255 = (s_primal_ctx_abs_0(_S1254.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_2;
    Matrix<float, 3, 3>  _S1256;
    float _S1257;
    if(_S1255)
    {
        float inv_s_0 = 1.0f / _S1254.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1258 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1259 = _S1254.rows[int(2)].z * _S1254.rows[int(2)].z;
        H_2 = _S1254 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1256 = _S1258;
        _S1257 = _S1259;
    }
    else
    {
        H_2 = _S1254;
        _S1256 = _S1199;
        _S1257 = 0.0f;
    }
    float _S1260 = _S1236.x;
    float _S1261 = _S1236.y;
    float intensity_1 = _S1260 + _S1261 + _S1236.z;
    float3  rgi_in_0 = make_float3 (_S1260, _S1261, intensity_1);
    float3  _S1262 = s_primal_ctx_mul_2(H_2, rgi_in_0);
    float _S1263 = _S1262.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_1 / _S1263;
    float3  _S1264 = make_float3 (norm_factor_0);
    float _S1265 = _S1263 * _S1263;
    float3  rgi_out_2 = _S1262 * make_float3 (norm_factor_0);
    float _S1266 = rgi_out_2.x;
    float _S1267 = rgi_out_2.y;
    float3  _S1268 = make_float3 (_S1266, _S1267, rgi_out_2.z - _S1266 - _S1267);
    float3  _S1269 = make_float3 (0.0f);
    float3  _S1270 = make_float3 (1.0f);
    float3  _S1271 = s_primal_ctx_clamp_1(_S1268, _S1269, _S1270);
    float _S1272 = _S1271.x;
    float _S1273 = 1.0f + s_primal_ctx_exp_0(_S1216);
    float _S1274 = 0.30000001192092896f + s_primal_ctx_log_0(_S1273);
    float _S1275 = 1.0f + s_primal_ctx_exp_0(_S1217);
    float _S1276 = 0.30000001192092896f + s_primal_ctx_log_0(_S1275);
    float _S1277 = 1.0f + s_primal_ctx_exp_0(_S1218);
    float _S1278 = 0.10000000149011612f + s_primal_ctx_log_0(_S1277);
    float _S1279 = - _S1219;
    float _S1280 = 1.0f + s_primal_ctx_exp_0(_S1279);
    float _S1281 = 1.0f / _S1280;
    float _S1282 = _S1280 * _S1280;
    float _S1283 = _S1276 - _S1274;
    float _S1284 = s_primal_ctx_fma_0(_S1283, _S1281, _S1274);
    float _S1285 = _S1276 * _S1281;
    float a_4 = _S1285 / _S1284;
    float _S1286 = _S1284 * _S1284;
    float b_5 = 1.0f - a_4;
    bool _S1287 = _S1272 <= _S1281;
    float y_11;
    float _S1288;
    float _S1289;
    float _S1290;
    float _S1291;
    float _S1292;
    float _S1293;
    float _S1294;
    float _S1295;
    if(_S1287)
    {
        float _S1296 = _S1272 / _S1281;
        float _S1297 = _S1281 * _S1281;
        float _S1298 = s_primal_ctx_pow_0(_S1296, _S1274);
        y_11 = a_4 * _S1298;
        _S1288 = 0.0f;
        _S1289 = 0.0f;
        _S1290 = 0.0f;
        _S1291 = 0.0f;
        _S1292 = 0.0f;
        _S1293 = _S1298;
        _S1294 = _S1296;
        _S1295 = _S1297;
    }
    else
    {
        float _S1299 = 1.0f - _S1272;
        float _S1300 = 1.0f - _S1281;
        float _S1301 = _S1299 / _S1300;
        float _S1302 = _S1300 * _S1300;
        float _S1303 = s_primal_ctx_pow_0(_S1301, _S1276);
        y_11 = 1.0f - b_5 * _S1303;
        _S1288 = _S1303;
        _S1289 = _S1301;
        _S1290 = _S1302;
        _S1291 = _S1299;
        _S1292 = _S1300;
        _S1293 = 0.0f;
        _S1294 = 0.0f;
        _S1295 = 0.0f;
    }
    float _S1304 = s_primal_ctx_max_0(0.0f, y_11);
    float _S1305 = _S1271.y;
    float _S1306 = 1.0f + s_primal_ctx_exp_0(_S1220);
    float _S1307 = 0.30000001192092896f + s_primal_ctx_log_0(_S1306);
    float _S1308 = 1.0f + s_primal_ctx_exp_0(_S1221);
    float _S1309 = 0.30000001192092896f + s_primal_ctx_log_0(_S1308);
    float _S1310 = 1.0f + s_primal_ctx_exp_0(_S1222);
    float _S1311 = 0.10000000149011612f + s_primal_ctx_log_0(_S1310);
    float _S1312 = - _S1223;
    float _S1313 = 1.0f + s_primal_ctx_exp_0(_S1312);
    float _S1314 = 1.0f / _S1313;
    float _S1315 = _S1313 * _S1313;
    float _S1316 = _S1309 - _S1307;
    float _S1317 = s_primal_ctx_fma_0(_S1316, _S1314, _S1307);
    float _S1318 = _S1309 * _S1314;
    float a_5 = _S1318 / _S1317;
    float _S1319 = _S1317 * _S1317;
    float b_6 = 1.0f - a_5;
    bool _S1320 = _S1305 <= _S1314;
    float y_12;
    float _S1321;
    float _S1322;
    float _S1323;
    float _S1324;
    float _S1325;
    float _S1326;
    float _S1327;
    float _S1328;
    if(_S1320)
    {
        float _S1329 = _S1305 / _S1314;
        float _S1330 = _S1314 * _S1314;
        float _S1331 = s_primal_ctx_pow_0(_S1329, _S1307);
        y_12 = a_5 * _S1331;
        _S1321 = 0.0f;
        _S1322 = 0.0f;
        _S1323 = 0.0f;
        _S1324 = 0.0f;
        _S1325 = 0.0f;
        _S1326 = _S1331;
        _S1327 = _S1329;
        _S1328 = _S1330;
    }
    else
    {
        float _S1332 = 1.0f - _S1305;
        float _S1333 = 1.0f - _S1314;
        float _S1334 = _S1332 / _S1333;
        float _S1335 = _S1333 * _S1333;
        float _S1336 = s_primal_ctx_pow_0(_S1334, _S1309);
        y_12 = 1.0f - b_6 * _S1336;
        _S1321 = _S1336;
        _S1322 = _S1334;
        _S1323 = _S1335;
        _S1324 = _S1332;
        _S1325 = _S1333;
        _S1326 = 0.0f;
        _S1327 = 0.0f;
        _S1328 = 0.0f;
    }
    float _S1337 = s_primal_ctx_max_0(0.0f, y_12);
    float _S1338 = _S1271.z;
    float _S1339 = 1.0f + s_primal_ctx_exp_0(_S1224);
    float _S1340 = 0.30000001192092896f + s_primal_ctx_log_0(_S1339);
    float _S1341 = 1.0f + s_primal_ctx_exp_0(_S1225);
    float _S1342 = 0.30000001192092896f + s_primal_ctx_log_0(_S1341);
    float _S1343 = 1.0f + s_primal_ctx_exp_0(_S1226);
    float _S1344 = 0.10000000149011612f + s_primal_ctx_log_0(_S1343);
    float _S1345 = - _S1227;
    float _S1346 = 1.0f + s_primal_ctx_exp_0(_S1345);
    float _S1347 = 1.0f / _S1346;
    float _S1348 = _S1346 * _S1346;
    float _S1349 = _S1342 - _S1340;
    float _S1350 = s_primal_ctx_fma_0(_S1349, _S1347, _S1340);
    float _S1351 = _S1342 * _S1347;
    float a_6 = _S1351 / _S1350;
    float _S1352 = _S1350 * _S1350;
    float b_7 = 1.0f - a_6;
    bool _S1353 = _S1338 <= _S1347;
    float y_13;
    float _S1354;
    float _S1355;
    float _S1356;
    float _S1357;
    float _S1358;
    float _S1359;
    float _S1360;
    float _S1361;
    if(_S1353)
    {
        float _S1362 = _S1338 / _S1347;
        float _S1363 = _S1347 * _S1347;
        float _S1364 = s_primal_ctx_pow_0(_S1362, _S1340);
        y_13 = a_6 * _S1364;
        _S1354 = 0.0f;
        _S1355 = 0.0f;
        _S1356 = 0.0f;
        _S1357 = 0.0f;
        _S1358 = 0.0f;
        _S1359 = _S1364;
        _S1360 = _S1362;
        _S1361 = _S1363;
    }
    else
    {
        float _S1365 = 1.0f - _S1338;
        float _S1366 = 1.0f - _S1347;
        float _S1367 = _S1365 / _S1366;
        float _S1368 = _S1366 * _S1366;
        float _S1369 = s_primal_ctx_pow_0(_S1367, _S1342);
        y_13 = 1.0f - b_7 * _S1369;
        _S1354 = _S1369;
        _S1355 = _S1367;
        _S1356 = _S1368;
        _S1357 = _S1365;
        _S1358 = _S1366;
        _S1359 = 0.0f;
        _S1360 = 0.0f;
        _S1361 = 0.0f;
    }
    float _S1370 = s_primal_ctx_max_0(0.0f, y_13);
    DiffPair_float_0 _S1371;
    (&_S1371)->primal_0 = _S1370;
    (&_S1371)->differential_0 = 0.0f;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1344;
    (&_S1372)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1371, &_S1372, _s_dOut_11.z);
    DiffPair_float_0 _S1373 = _S1372;
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = 0.0f;
    (&_S1374)->differential_0 = 0.0f;
    DiffPair_float_0 _S1375;
    (&_S1375)->primal_0 = y_13;
    (&_S1375)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1374, &_S1375, _S1371.differential_0);
    DiffPair_float_0 _S1376 = _S1375;
    if(_S1353)
    {
        float _S1377 = a_6 * _S1376.differential_0;
        float _S1378 = _S1359 * _S1376.differential_0;
        DiffPair_float_0 _S1379;
        (&_S1379)->primal_0 = _S1360;
        (&_S1379)->differential_0 = 0.0f;
        DiffPair_float_0 _S1380;
        (&_S1380)->primal_0 = _S1340;
        (&_S1380)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1379, &_S1380, _S1377);
        float _S1381 = _S1379.differential_0 / _S1361;
        float _S1382 = _S1338 * - _S1381;
        float _S1383 = _S1347 * _S1381;
        y_13 = 0.0f;
        _S1354 = _S1378;
        _S1355 = _S1382;
        _S1356 = 0.0f;
        _S1357 = _S1380.differential_0;
        _S1358 = _S1383;
    }
    else
    {
        float _S1384 = - _S1376.differential_0;
        float _S1385 = b_7 * _S1384;
        float _S1386 = _S1354 * _S1384;
        DiffPair_float_0 _S1387;
        (&_S1387)->primal_0 = _S1355;
        (&_S1387)->differential_0 = 0.0f;
        DiffPair_float_0 _S1388;
        (&_S1388)->primal_0 = _S1342;
        (&_S1388)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1387, &_S1388, _S1385);
        float _S1389 = _S1387.differential_0 / _S1356;
        float _S1390 = - (_S1357 * - _S1389);
        float _S1391 = - (_S1358 * _S1389);
        y_13 = _S1386;
        _S1354 = 0.0f;
        _S1355 = _S1390;
        _S1356 = _S1388.differential_0;
        _S1357 = 0.0f;
        _S1358 = _S1391;
    }
    float _S1392 = (- y_13 + _S1354) / _S1352;
    float _S1393 = _S1351 * - _S1392;
    float _S1394 = _S1350 * _S1392;
    float _S1395 = _S1342 * _S1394;
    float _S1396 = _S1347 * _S1394;
    DiffPair_float_0 _S1397;
    (&_S1397)->primal_0 = _S1349;
    (&_S1397)->differential_0 = 0.0f;
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1347;
    (&_S1398)->differential_0 = 0.0f;
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = _S1340;
    (&_S1399)->differential_0 = 0.0f;
    s_bwd_prop_fma_0(&_S1397, &_S1398, &_S1399, _S1393);
    float _S1400 = - _S1397.differential_0;
    float _S1401 = - ((_S1395 + _S1398.differential_0 + _S1355) / _S1348);
    DiffPair_float_0 _S1402;
    (&_S1402)->primal_0 = _S1345;
    (&_S1402)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1402, _S1401);
    float _S1403 = - _S1402.differential_0;
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = _S1343;
    (&_S1404)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1404, _S1373.differential_0);
    DiffPair_float_0 _S1405;
    (&_S1405)->primal_0 = _S1226;
    (&_S1405)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1405, _S1404.differential_0);
    DiffPair_float_0 _S1406 = _S1405;
    float _S1407 = _S1396 + _S1397.differential_0 + _S1356;
    DiffPair_float_0 _S1408;
    (&_S1408)->primal_0 = _S1341;
    (&_S1408)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1408, _S1407);
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1225;
    (&_S1409)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1409, _S1408.differential_0);
    DiffPair_float_0 _S1410 = _S1409;
    float _S1411 = _S1399.differential_0 + _S1400 + _S1357;
    DiffPair_float_0 _S1412;
    (&_S1412)->primal_0 = _S1339;
    (&_S1412)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1412, _S1411);
    DiffPair_float_0 _S1413;
    (&_S1413)->primal_0 = _S1224;
    (&_S1413)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1413, _S1412.differential_0);
    DiffPair_float_0 _S1414 = _S1413;
    float3  _S1415 = make_float3 (0.0f, 0.0f, _S1358);
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = _S1337;
    (&_S1416)->differential_0 = 0.0f;
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1311;
    (&_S1417)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1416, &_S1417, _s_dOut_11.y);
    DiffPair_float_0 _S1418 = _S1417;
    DiffPair_float_0 _S1419;
    (&_S1419)->primal_0 = 0.0f;
    (&_S1419)->differential_0 = 0.0f;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = y_12;
    (&_S1420)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1419, &_S1420, _S1416.differential_0);
    DiffPair_float_0 _S1421 = _S1420;
    if(_S1320)
    {
        float _S1422 = a_5 * _S1421.differential_0;
        float _S1423 = _S1326 * _S1421.differential_0;
        DiffPair_float_0 _S1424;
        (&_S1424)->primal_0 = _S1327;
        (&_S1424)->differential_0 = 0.0f;
        DiffPair_float_0 _S1425;
        (&_S1425)->primal_0 = _S1307;
        (&_S1425)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1424, &_S1425, _S1422);
        float _S1426 = _S1424.differential_0 / _S1328;
        float _S1427 = _S1305 * - _S1426;
        float _S1428 = _S1314 * _S1426;
        y_12 = 0.0f;
        _S1321 = _S1423;
        _S1322 = _S1427;
        _S1323 = 0.0f;
        _S1324 = _S1425.differential_0;
        _S1325 = _S1428;
    }
    else
    {
        float _S1429 = - _S1421.differential_0;
        float _S1430 = b_6 * _S1429;
        float _S1431 = _S1321 * _S1429;
        DiffPair_float_0 _S1432;
        (&_S1432)->primal_0 = _S1322;
        (&_S1432)->differential_0 = 0.0f;
        DiffPair_float_0 _S1433;
        (&_S1433)->primal_0 = _S1309;
        (&_S1433)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1432, &_S1433, _S1430);
        float _S1434 = _S1432.differential_0 / _S1323;
        float _S1435 = - (_S1324 * - _S1434);
        float _S1436 = - (_S1325 * _S1434);
        y_12 = _S1431;
        _S1321 = 0.0f;
        _S1322 = _S1435;
        _S1323 = _S1433.differential_0;
        _S1324 = 0.0f;
        _S1325 = _S1436;
    }
    float _S1437 = (- y_12 + _S1321) / _S1319;
    float _S1438 = _S1318 * - _S1437;
    float _S1439 = _S1317 * _S1437;
    float _S1440 = _S1309 * _S1439;
    float _S1441 = _S1314 * _S1439;
    DiffPair_float_0 _S1442;
    (&_S1442)->primal_0 = _S1316;
    (&_S1442)->differential_0 = 0.0f;
    DiffPair_float_0 _S1443;
    (&_S1443)->primal_0 = _S1314;
    (&_S1443)->differential_0 = 0.0f;
    DiffPair_float_0 _S1444;
    (&_S1444)->primal_0 = _S1307;
    (&_S1444)->differential_0 = 0.0f;
    s_bwd_prop_fma_0(&_S1442, &_S1443, &_S1444, _S1438);
    float _S1445 = - _S1442.differential_0;
    float _S1446 = - ((_S1440 + _S1443.differential_0 + _S1322) / _S1315);
    DiffPair_float_0 _S1447;
    (&_S1447)->primal_0 = _S1312;
    (&_S1447)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1447, _S1446);
    float _S1448 = - _S1447.differential_0;
    DiffPair_float_0 _S1449;
    (&_S1449)->primal_0 = _S1310;
    (&_S1449)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1449, _S1418.differential_0);
    DiffPair_float_0 _S1450;
    (&_S1450)->primal_0 = _S1222;
    (&_S1450)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1450, _S1449.differential_0);
    DiffPair_float_0 _S1451 = _S1450;
    float _S1452 = _S1441 + _S1442.differential_0 + _S1323;
    DiffPair_float_0 _S1453;
    (&_S1453)->primal_0 = _S1308;
    (&_S1453)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1453, _S1452);
    DiffPair_float_0 _S1454;
    (&_S1454)->primal_0 = _S1221;
    (&_S1454)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1454, _S1453.differential_0);
    DiffPair_float_0 _S1455 = _S1454;
    float _S1456 = _S1444.differential_0 + _S1445 + _S1324;
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1306;
    (&_S1457)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1457, _S1456);
    DiffPair_float_0 _S1458;
    (&_S1458)->primal_0 = _S1220;
    (&_S1458)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1458, _S1457.differential_0);
    DiffPair_float_0 _S1459 = _S1458;
    float3  _S1460 = _S1415 + make_float3 (0.0f, _S1325, 0.0f);
    DiffPair_float_0 _S1461;
    (&_S1461)->primal_0 = _S1304;
    (&_S1461)->differential_0 = 0.0f;
    DiffPair_float_0 _S1462;
    (&_S1462)->primal_0 = _S1278;
    (&_S1462)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1461, &_S1462, _s_dOut_11.x);
    DiffPair_float_0 _S1463 = _S1462;
    DiffPair_float_0 _S1464;
    (&_S1464)->primal_0 = 0.0f;
    (&_S1464)->differential_0 = 0.0f;
    DiffPair_float_0 _S1465;
    (&_S1465)->primal_0 = y_11;
    (&_S1465)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1464, &_S1465, _S1461.differential_0);
    DiffPair_float_0 _S1466 = _S1465;
    if(_S1287)
    {
        float _S1467 = a_4 * _S1466.differential_0;
        float _S1468 = _S1293 * _S1466.differential_0;
        DiffPair_float_0 _S1469;
        (&_S1469)->primal_0 = _S1294;
        (&_S1469)->differential_0 = 0.0f;
        DiffPair_float_0 _S1470;
        (&_S1470)->primal_0 = _S1274;
        (&_S1470)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1469, &_S1470, _S1467);
        float _S1471 = _S1469.differential_0 / _S1295;
        float _S1472 = _S1272 * - _S1471;
        float _S1473 = _S1281 * _S1471;
        y_11 = 0.0f;
        _S1288 = _S1468;
        _S1289 = _S1472;
        _S1290 = 0.0f;
        _S1291 = _S1470.differential_0;
        _S1292 = _S1473;
    }
    else
    {
        float _S1474 = - _S1466.differential_0;
        float _S1475 = b_5 * _S1474;
        float _S1476 = _S1288 * _S1474;
        DiffPair_float_0 _S1477;
        (&_S1477)->primal_0 = _S1289;
        (&_S1477)->differential_0 = 0.0f;
        DiffPair_float_0 _S1478;
        (&_S1478)->primal_0 = _S1276;
        (&_S1478)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1477, &_S1478, _S1475);
        float _S1479 = _S1477.differential_0 / _S1290;
        float _S1480 = - (_S1291 * - _S1479);
        float _S1481 = - (_S1292 * _S1479);
        y_11 = _S1476;
        _S1288 = 0.0f;
        _S1289 = _S1480;
        _S1290 = _S1478.differential_0;
        _S1291 = 0.0f;
        _S1292 = _S1481;
    }
    float _S1482 = (- y_11 + _S1288) / _S1286;
    float _S1483 = _S1285 * - _S1482;
    float _S1484 = _S1284 * _S1482;
    float _S1485 = _S1276 * _S1484;
    float _S1486 = _S1281 * _S1484;
    DiffPair_float_0 _S1487;
    (&_S1487)->primal_0 = _S1283;
    (&_S1487)->differential_0 = 0.0f;
    DiffPair_float_0 _S1488;
    (&_S1488)->primal_0 = _S1281;
    (&_S1488)->differential_0 = 0.0f;
    DiffPair_float_0 _S1489;
    (&_S1489)->primal_0 = _S1274;
    (&_S1489)->differential_0 = 0.0f;
    s_bwd_prop_fma_0(&_S1487, &_S1488, &_S1489, _S1483);
    float _S1490 = - _S1487.differential_0;
    float _S1491 = - ((_S1485 + _S1488.differential_0 + _S1289) / _S1282);
    DiffPair_float_0 _S1492;
    (&_S1492)->primal_0 = _S1279;
    (&_S1492)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1492, _S1491);
    float _S1493 = - _S1492.differential_0;
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = _S1277;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1494, _S1463.differential_0);
    DiffPair_float_0 _S1495;
    (&_S1495)->primal_0 = _S1218;
    (&_S1495)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1495, _S1494.differential_0);
    DiffPair_float_0 _S1496 = _S1495;
    float _S1497 = _S1486 + _S1487.differential_0 + _S1290;
    DiffPair_float_0 _S1498;
    (&_S1498)->primal_0 = _S1275;
    (&_S1498)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1498, _S1497);
    DiffPair_float_0 _S1499;
    (&_S1499)->primal_0 = _S1217;
    (&_S1499)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1499, _S1498.differential_0);
    DiffPair_float_0 _S1500 = _S1499;
    float _S1501 = _S1489.differential_0 + _S1490 + _S1291;
    DiffPair_float_0 _S1502;
    (&_S1502)->primal_0 = _S1273;
    (&_S1502)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1502, _S1501);
    DiffPair_float_0 _S1503;
    (&_S1503)->primal_0 = _S1216;
    (&_S1503)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1503, _S1502.differential_0);
    DiffPair_float_0 _S1504 = _S1503;
    float3  _S1505 = _S1460 + make_float3 (_S1292, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1506;
    (&_S1506)->primal_0 = _S1268;
    (&_S1506)->differential_0 = _S1198;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1507;
    (&_S1507)->primal_0 = _S1269;
    (&_S1507)->differential_0 = _S1198;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1508;
    (&_S1508)->primal_0 = _S1270;
    (&_S1508)->differential_0 = _S1198;
    s_bwd_prop_clamp_1(&_S1506, &_S1507, &_S1508, _S1505);
    float _S1509 = - _S1506.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1506.differential_0.x + _S1509, _S1506.differential_0.y + _S1509, _S1506.differential_0.z);
    float3  _S1510 = _S1262 * s_diff_rgi_out_T_0;
    float _S1511 = (_S1510.x + _S1510.y + _S1510.z) / _S1265;
    float _S1512 = _S1263 * _S1511;
    float3  _S1513 = _S1264 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_1 * - _S1511);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1514;
    (&_S1514)->primal_0 = H_2;
    (&_S1514)->differential_0 = _S1199;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1515;
    (&_S1515)->primal_0 = rgi_in_0;
    (&_S1515)->differential_0 = _S1198;
    s_bwd_prop_mul_1(&_S1514, &_S1515, _S1513);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1516 = _S1514;
    float _S1517 = _S1512 + _S1515.differential_0.z;
    float _S1518 = _S1515.differential_0.y + _S1517;
    float _S1519 = _S1515.differential_0.x + _S1517;
    float3  _S1520 = make_float3 (_S1519, _S1518, _S1517);
    if(_S1255)
    {
        Matrix<float, 3, 3>  _S1521 = _S1254 * _S1516.differential_0;
        Matrix<float, 3, 3>  _S1522 = _S1256 * _S1516.differential_0;
        _S1257 = - ((_S1521.rows[int(0)].x + _S1521.rows[int(0)].y + _S1521.rows[int(0)].z + _S1521.rows[int(1)].x + _S1521.rows[int(1)].y + _S1521.rows[int(1)].z + _S1521.rows[int(2)].x + _S1521.rows[int(2)].y + _S1521.rows[int(2)].z) / _S1257);
        H_2 = _S1522;
    }
    else
    {
        _S1257 = 0.0f;
        H_2 = _S1516.differential_0;
    }
    DiffPair_float_0 _S1523;
    (&_S1523)->primal_0 = _S1254.rows[int(2)].z;
    (&_S1523)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1523, 0.0f);
    float _S1524 = _S1523.differential_0 + _S1257;
    float3  _S1525 = _S1198;
    *&((&_S1525)->z) = _S1524;
    Matrix<float, 3, 3>  _S1526 = _S1199;
    _S1526[int(2)] = _S1525;
    Matrix<float, 3, 3>  _S1527 = H_2 + _S1526;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1528;
    (&_S1528)->primal_0 = _S1253;
    (&_S1528)->differential_0 = _S1199;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1529;
    (&_S1529)->primal_0 = S_inv_0;
    (&_S1529)->differential_0 = _S1199;
    s_bwd_prop_mul_2(&_S1528, &_S1529, _S1527);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1530;
    (&_S1530)->primal_0 = T_1;
    (&_S1530)->differential_0 = _S1199;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1531;
    (&_S1531)->primal_0 = D_0;
    (&_S1531)->differential_0 = _S1199;
    s_bwd_prop_mul_2(&_S1530, &_S1531, _S1528.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1532 = _S1530;
    float3  _S1533 = make_float3 (_S1531.differential_0.rows[int(0)].x, _S1531.differential_0.rows[int(1)].y, _S1531.differential_0.rows[int(2)].z);
    float3  _S1534;
    if(_S1248)
    {
        if(_S1250)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1535;
            (&_S1535)->primal_0 = r1_1;
            (&_S1535)->differential_0 = _S1198;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1536;
            (&_S1536)->primal_0 = r2_11;
            (&_S1536)->differential_0 = _S1198;
            s_bwd_prop_cross_0(&_S1535, &_S1536, _S1533);
            _S1236 = _S1198;
            lambda_v_3 = _S1536.differential_0;
            _S1534 = _S1535.differential_0;
        }
        else
        {
            _S1236 = _S1533;
            lambda_v_3 = _S1198;
            _S1534 = _S1198;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1537;
        (&_S1537)->primal_0 = _S1249;
        (&_S1537)->differential_0 = _S1198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1538;
        (&_S1538)->primal_0 = _S1249;
        (&_S1538)->differential_0 = _S1198;
        s_bwd_prop_dot_0(&_S1537, &_S1538, 0.0f);
        float3  _S1539 = _S1538.differential_0 + _S1537.differential_0 + _S1236;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1540;
        (&_S1540)->primal_0 = r0_1;
        (&_S1540)->differential_0 = _S1198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1541;
        (&_S1541)->primal_0 = r2_11;
        (&_S1541)->differential_0 = _S1198;
        s_bwd_prop_cross_0(&_S1540, &_S1541, _S1539);
        float3  _S1542 = _S1541.differential_0 + lambda_v_3;
        _S1236 = _S1198;
        lambda_v_3 = _S1542;
        _S1249 = _S1534;
        _S1534 = _S1540.differential_0;
    }
    else
    {
        _S1236 = _S1533;
        lambda_v_3 = _S1198;
        _S1249 = _S1198;
        _S1534 = _S1198;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1543;
    (&_S1543)->primal_0 = _S1247;
    (&_S1543)->differential_0 = _S1198;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1544;
    (&_S1544)->primal_0 = _S1247;
    (&_S1544)->differential_0 = _S1198;
    s_bwd_prop_dot_0(&_S1543, &_S1544, 0.0f);
    float3  _S1545 = _S1544.differential_0 + _S1543.differential_0 + _S1236;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1546;
    (&_S1546)->primal_0 = r0_1;
    (&_S1546)->differential_0 = _S1198;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1547;
    (&_S1547)->primal_0 = r1_1;
    (&_S1547)->differential_0 = _S1198;
    s_bwd_prop_cross_0(&_S1546, &_S1547, _S1545);
    float3  _S1548 = _S1198;
    *&((&_S1548)->z) = lambda_v_3.z;
    *&((&_S1548)->y) = lambda_v_3.y;
    *&((&_S1548)->x) = lambda_v_3.x;
    float3  _S1549 = _S1547.differential_0 + _S1249;
    float3  _S1550 = _S1198;
    *&((&_S1550)->z) = _S1549.z;
    *&((&_S1550)->y) = _S1549.y;
    *&((&_S1550)->x) = _S1549.x;
    float3  _S1551 = _S1546.differential_0 + _S1534;
    float3  _S1552 = _S1198;
    *&((&_S1552)->z) = _S1551.z;
    *&((&_S1552)->y) = _S1551.y;
    *&((&_S1552)->x) = _S1551.x;
    Matrix<float, 3, 3>  _S1553 = _S1199;
    _S1553[int(2)] = _S1548;
    _S1553[int(1)] = _S1550;
    _S1553[int(0)] = _S1552;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1554;
    (&_S1554)->primal_0 = skew_0;
    (&_S1554)->differential_0 = _S1199;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1555;
    (&_S1555)->primal_0 = T_1;
    (&_S1555)->differential_0 = _S1199;
    s_bwd_prop_mul_2(&_S1554, &_S1555, _S1553);
    Matrix<float, 3, 3>  _S1556 = _S1555.differential_0 + _S1532.differential_0;
    float2  _S1557 = make_float2 (_S1554.differential_0.rows[int(2)].y + - _S1554.differential_0.rows[int(1)].z, _S1554.differential_0.rows[int(0)].z + - _S1554.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1558 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1559;
    (&_S1559)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1559)->differential_0 = _S1558;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1560;
    (&_S1560)->primal_0 = _S1239.color_params_0.n_0;
    (&_S1560)->differential_0 = _S1202;
    s_bwd_prop_mul_3(&_S1559, &_S1560, _S1557);
    float2  _S1561 = make_float2 (_S1556.rows[int(0)].z, _S1556.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1562;
    (&_S1562)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1562)->differential_0 = _S1558;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1563;
    (&_S1563)->primal_0 = _S1239.color_params_0.g_0;
    (&_S1563)->differential_0 = _S1202;
    s_bwd_prop_mul_3(&_S1562, &_S1563, _S1561);
    float2  _S1564 = make_float2 (_S1556.rows[int(0)].y, _S1556.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1565;
    (&_S1565)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1565)->differential_0 = _S1558;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1566;
    (&_S1566)->primal_0 = _S1239.color_params_0.r_0;
    (&_S1566)->differential_0 = _S1202;
    s_bwd_prop_mul_3(&_S1565, &_S1566, _S1564);
    float2  _S1567 = make_float2 (_S1556.rows[int(0)].x, _S1556.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1568;
    (&_S1568)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1568)->differential_0 = _S1558;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1569;
    (&_S1569)->primal_0 = _S1239.color_params_0.b_0;
    (&_S1569)->differential_0 = _S1202;
    s_bwd_prop_mul_3(&_S1568, &_S1569, _S1567);
    ColorPPISPParams_0 _S1570 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1570)->n_0 = _S1560.differential_0;
    (&_S1570)->g_0 = _S1563.differential_0;
    (&_S1570)->r_0 = _S1566.differential_0;
    (&_S1570)->b_0 = _S1569.differential_0;
    _S1236 = _S1520;
    *&((&_S1236)->z) = 0.0f;
    float _S1571 = rgb_out_2.z * _S1517;
    float _S1572 = _S1238 * _S1517;
    DiffPair_float_0 _S1573;
    (&_S1573)->primal_0 = falloff_2;
    (&_S1573)->differential_0 = 0.0f;
    DiffPair_float_0 _S1574;
    (&_S1574)->primal_0 = 0.0f;
    (&_S1574)->differential_0 = 0.0f;
    DiffPair_float_0 _S1575;
    (&_S1575)->primal_0 = 1.0f;
    (&_S1575)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1573, &_S1574, &_S1575, _S1571);
    float _S1576 = r2_10 * _S1573.differential_0;
    float _S1577 = r4_5 * _S1573.differential_0;
    float s_diff_r6_T_0 = _S1215 * _S1573.differential_0;
    float _S1578 = r6_2 * _S1573.differential_0;
    float _S1579 = r2_10 * (_S1214 * _S1573.differential_0 + r2_10 * s_diff_r6_T_0);
    float _S1580 = _S1213 * _S1573.differential_0 + r4_5 * s_diff_r6_T_0 + _S1579 + _S1579;
    float _S1581 = dy_5 * _S1580;
    float _S1582 = dx_7 * _S1580;
    float _S1583 = - (_S1581 + _S1581);
    float _S1584 = - (_S1582 + _S1582);
    *&((&_S1236)->y) = 0.0f;
    float _S1585 = rgb_out_2.y * _S1518;
    float _S1586 = _S1237 * _S1518;
    DiffPair_float_0 _S1587;
    (&_S1587)->primal_0 = falloff_1;
    (&_S1587)->differential_0 = 0.0f;
    DiffPair_float_0 _S1588;
    (&_S1588)->primal_0 = 0.0f;
    (&_S1588)->differential_0 = 0.0f;
    DiffPair_float_0 _S1589;
    (&_S1589)->primal_0 = 1.0f;
    (&_S1589)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1587, &_S1588, &_S1589, _S1585);
    float _S1590 = r2_9 * _S1587.differential_0;
    float _S1591 = r4_4 * _S1587.differential_0;
    float s_diff_r6_T_1 = _S1212 * _S1587.differential_0;
    float _S1592 = r6_1 * _S1587.differential_0;
    float _S1593 = r2_9 * (_S1211 * _S1587.differential_0 + r2_9 * s_diff_r6_T_1);
    float _S1594 = _S1210 * _S1587.differential_0 + r4_4 * s_diff_r6_T_1 + _S1593 + _S1593;
    float _S1595 = dy_4 * _S1594;
    float _S1596 = dx_6 * _S1594;
    float _S1597 = - (_S1595 + _S1595);
    float _S1598 = - (_S1596 + _S1596);
    *&((&_S1236)->x) = 0.0f;
    float _S1599 = rgb_out_2.x * _S1519;
    float _S1600 = _S1234 * _S1519;
    DiffPair_float_0 _S1601;
    (&_S1601)->primal_0 = falloff_0;
    (&_S1601)->differential_0 = 0.0f;
    DiffPair_float_0 _S1602;
    (&_S1602)->primal_0 = 0.0f;
    (&_S1602)->differential_0 = 0.0f;
    DiffPair_float_0 _S1603;
    (&_S1603)->primal_0 = 1.0f;
    (&_S1603)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1601, &_S1602, &_S1603, _S1599);
    float _S1604 = r2_8 * _S1601.differential_0;
    float _S1605 = r4_3 * _S1601.differential_0;
    float s_diff_r6_T_2 = _S1209 * _S1601.differential_0;
    float _S1606 = r6_0 * _S1601.differential_0;
    float _S1607 = r2_8 * (_S1208 * _S1601.differential_0 + r2_8 * s_diff_r6_T_2);
    float _S1608 = _S1207 * _S1601.differential_0 + r4_3 * s_diff_r6_T_2 + _S1607 + _S1607;
    float _S1609 = dy_3 * _S1608;
    float _S1610 = dx_5 * _S1608;
    float _S1611 = - (_S1609 + _S1609);
    float _S1612 = - (_S1610 + _S1610);
    float3  _S1613 = _S1198;
    *&((&_S1613)->z) = _S1572;
    *&((&_S1613)->y) = _S1586;
    *&((&_S1613)->x) = _S1600;
    float3  _S1614 = _S1236 + _S1613;
    float3  _S1615 = _S1197.primal_0 * _S1614;
    float3  _S1616 = _S1230 * _S1614;
    float _S1617 = _S1615.x + _S1615.y + _S1615.z;
    DiffPair_float_0 _S1618;
    (&_S1618)->primal_0 = _S1228.exposure_0;
    (&_S1618)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1618, _S1617);
    PPISPParams_0 _S1619 = PPISPParams_x24_syn_dzero_0();
    (&_S1619)->color_params_0 = _S1570;
    (&_S1619)->exposure_0 = _S1618.differential_0;
    _S1206 = _S1619;
    (&(&_S1206)->crf_params_0[int(2)])->center_0 = 0.0f;
    float _S1620 = _S1619.crf_params_0[int(2)].center_0 + _S1403;
    (&(&_S1206)->crf_params_0[int(2)])->gamma_0 = 0.0f;
    float _S1621 = _S1619.crf_params_0[int(2)].gamma_0 + _S1406.differential_0;
    (&(&_S1206)->crf_params_0[int(2)])->shoulder_0 = 0.0f;
    float _S1622 = _S1619.crf_params_0[int(2)].shoulder_0 + _S1410.differential_0;
    (&(&_S1206)->crf_params_0[int(2)])->toe_0 = 0.0f;
    float _S1623 = _S1619.crf_params_0[int(2)].toe_0 + _S1414.differential_0;
    (&(&_S1206)->crf_params_0[int(1)])->center_0 = 0.0f;
    float _S1624 = _S1619.crf_params_0[int(1)].center_0 + _S1448;
    (&(&_S1206)->crf_params_0[int(1)])->gamma_0 = 0.0f;
    float _S1625 = _S1619.crf_params_0[int(1)].gamma_0 + _S1451.differential_0;
    (&(&_S1206)->crf_params_0[int(1)])->shoulder_0 = 0.0f;
    float _S1626 = _S1619.crf_params_0[int(1)].shoulder_0 + _S1455.differential_0;
    (&(&_S1206)->crf_params_0[int(1)])->toe_0 = 0.0f;
    float _S1627 = _S1619.crf_params_0[int(1)].toe_0 + _S1459.differential_0;
    (&(&_S1206)->crf_params_0[int(0)])->center_0 = 0.0f;
    float _S1628 = _S1619.crf_params_0[int(0)].center_0 + _S1493;
    (&(&_S1206)->crf_params_0[int(0)])->gamma_0 = 0.0f;
    float _S1629 = _S1619.crf_params_0[int(0)].gamma_0 + _S1496.differential_0;
    (&(&_S1206)->crf_params_0[int(0)])->shoulder_0 = 0.0f;
    float _S1630 = _S1619.crf_params_0[int(0)].shoulder_0 + _S1500.differential_0;
    (&(&_S1206)->crf_params_0[int(0)])->toe_0 = 0.0f;
    float _S1631 = _S1619.crf_params_0[int(0)].toe_0 + _S1504.differential_0;
    *&((&(&(&_S1206)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1206)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1206)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1632 = _S1578 + _S1619.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1206)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1633 = _S1577 + _S1619.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1206)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1634 = _S1576 + _S1619.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1206)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1635 = _S1583 + _S1619.vignette_params_0[int(2)].cy_0;
    (&(&_S1206)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1636 = _S1584 + _S1619.vignette_params_0[int(2)].cx_0;
    (&(&_S1206)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1637 = _S1592 + _S1619.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1206)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1638 = _S1591 + _S1619.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1206)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1639 = _S1590 + _S1619.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1206)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1640 = _S1597 + _S1619.vignette_params_0[int(1)].cy_0;
    (&(&_S1206)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1641 = _S1598 + _S1619.vignette_params_0[int(1)].cx_0;
    (&(&_S1206)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1642 = _S1606 + _S1619.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1206)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1643 = _S1605 + _S1619.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1206)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1644 = _S1604 + _S1619.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1206)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1645 = _S1611 + _S1619.vignette_params_0[int(0)].cy_0;
    (&(&_S1206)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1646 = _S1612 + _S1619.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 36>  _S1647;
    _S1647[int(0)] = 0.0f;
    _S1647[int(1)] = 0.0f;
    _S1647[int(2)] = 0.0f;
    _S1647[int(3)] = 0.0f;
    _S1647[int(4)] = 0.0f;
    _S1647[int(5)] = 0.0f;
    _S1647[int(6)] = 0.0f;
    _S1647[int(7)] = 0.0f;
    _S1647[int(8)] = 0.0f;
    _S1647[int(9)] = 0.0f;
    _S1647[int(10)] = 0.0f;
    _S1647[int(11)] = 0.0f;
    _S1647[int(12)] = 0.0f;
    _S1647[int(13)] = 0.0f;
    _S1647[int(14)] = 0.0f;
    _S1647[int(15)] = 0.0f;
    _S1647[int(16)] = 0.0f;
    _S1647[int(17)] = 0.0f;
    _S1647[int(18)] = 0.0f;
    _S1647[int(19)] = 0.0f;
    _S1647[int(20)] = 0.0f;
    _S1647[int(21)] = 0.0f;
    _S1647[int(22)] = 0.0f;
    _S1647[int(23)] = 0.0f;
    _S1647[int(24)] = 0.0f;
    _S1647[int(25)] = 0.0f;
    _S1647[int(26)] = 0.0f;
    _S1647[int(27)] = 0.0f;
    _S1647[int(28)] = 0.0f;
    _S1647[int(29)] = 0.0f;
    _S1647[int(30)] = 0.0f;
    _S1647[int(31)] = 0.0f;
    _S1647[int(32)] = 0.0f;
    _S1647[int(33)] = 0.0f;
    _S1647[int(34)] = 0.0f;
    _S1647[int(35)] = 0.0f;
    _S1647[int(8)] = _S1639;
    _S1647[int(16)] = _S1619.color_params_0.b_0.x;
    _S1647[int(15)] = _S1632;
    _S1647[int(14)] = _S1633;
    _S1647[int(13)] = _S1634;
    _S1647[int(12)] = _S1635;
    _S1647[int(11)] = _S1636;
    _S1647[int(10)] = _S1637;
    _S1647[int(9)] = _S1638;
    _S1647[int(17)] = _S1619.color_params_0.b_0.y;
    _S1647[int(7)] = _S1640;
    _S1647[int(6)] = _S1641;
    _S1647[int(5)] = _S1642;
    _S1647[int(4)] = _S1643;
    _S1647[int(3)] = _S1644;
    _S1647[int(2)] = _S1645;
    _S1647[int(1)] = _S1646;
    _S1647[int(0)] = _S1206.exposure_0;
    _S1647[int(26)] = _S1629;
    _S1647[int(34)] = _S1621;
    _S1647[int(33)] = _S1622;
    _S1647[int(32)] = _S1623;
    _S1647[int(31)] = _S1624;
    _S1647[int(30)] = _S1625;
    _S1647[int(29)] = _S1626;
    _S1647[int(28)] = _S1627;
    _S1647[int(27)] = _S1628;
    _S1647[int(35)] = _S1620;
    _S1647[int(25)] = _S1630;
    _S1647[int(24)] = _S1631;
    _S1647[int(23)] = _S1619.color_params_0.n_0.y;
    _S1647[int(22)] = _S1619.color_params_0.n_0.x;
    _S1647[int(21)] = _S1619.color_params_0.g_0.y;
    _S1647[int(20)] = _S1619.color_params_0.g_0.x;
    _S1647[int(19)] = _S1619.color_params_0.r_0.y;
    _S1647[int(18)] = _S1619.color_params_0.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S1647;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S1616;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1648, float2  _S1649, float2  _S1650, float2  _S1651, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1652, float3  _S1653)
{
    s_bwd_prop_apply_ppisp_0(_S1648, _S1649, _S1650, _S1651, _S1652, _S1653);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_1, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, FixedArray<float, 36>  * params_1, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S1654 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_1;
    (&dp_rgb_in_0)->differential_0 = _S1654;
    FixedArray<float, 36>  _S1655 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = *params_1;
    (&dp_params_0)->differential_0 = _S1655;
    s_bwd_apply_ppisp_0(&dp_rgb_in_0, pix_coord_2, image_center_2, img_size_2, &dp_params_0, grad_out_0);
    *grad_rgb_in_0 = dp_rgb_in_0.differential_0;
    *grad_params_0 = (&dp_params_0)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  * params_2, FixedArray<float, 27>  * _S1656)
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
    FixedArray<float, 27>  losses_3;
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
    losses_3[int(22)] = 0.0f;
    losses_3[int(23)] = 0.0f;
    losses_3[int(24)] = 0.0f;
    losses_3[int(25)] = 0.0f;
    losses_3[int(26)] = 0.0f;
    losses_3[int(0)] = p_1.exposure_0;
    losses_3[int(1)] = p_1.vignette_params_0[int(0)].cx_0;
    losses_3[int(2)] = p_1.vignette_params_0[int(0)].cy_0;
    losses_3[int(3)] = p_1.vignette_params_0[int(1)].cx_0;
    losses_3[int(4)] = p_1.vignette_params_0[int(1)].cy_0;
    losses_3[int(5)] = p_1.vignette_params_0[int(2)].cx_0;
    losses_3[int(6)] = p_1.vignette_params_0[int(2)].cy_0;
    losses_3[int(7)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha0_0)));
    losses_3[int(8)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha1_0)));
    losses_3[int(9)] = (F32_max((0.0f), (p_1.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_1.vignette_params_0[int(2)].alpha2_0)));
    float mean_3 = (p_1.vignette_params_0[int(0)].cx_0 + p_1.vignette_params_0[int(1)].cx_0 + p_1.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S1657 = p_1.vignette_params_0[int(0)].cx_0 - mean_3;
    float _S1658 = p_1.vignette_params_0[int(1)].cx_0 - mean_3;
    float _S1659 = p_1.vignette_params_0[int(2)].cx_0 - mean_3;
    losses_3[int(10)] = (_S1657 * _S1657 + _S1658 * _S1658 + _S1659 * _S1659) / 3.0f;
    float mean_4 = (p_1.vignette_params_0[int(0)].cy_0 + p_1.vignette_params_0[int(1)].cy_0 + p_1.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S1660 = p_1.vignette_params_0[int(0)].cy_0 - mean_4;
    float _S1661 = p_1.vignette_params_0[int(1)].cy_0 - mean_4;
    float _S1662 = p_1.vignette_params_0[int(2)].cy_0 - mean_4;
    losses_3[int(11)] = (_S1660 * _S1660 + _S1661 * _S1661 + _S1662 * _S1662) / 3.0f;
    float mean_5 = (p_1.vignette_params_0[int(0)].alpha0_0 + p_1.vignette_params_0[int(1)].alpha0_0 + p_1.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S1663 = p_1.vignette_params_0[int(0)].alpha0_0 - mean_5;
    float _S1664 = p_1.vignette_params_0[int(1)].alpha0_0 - mean_5;
    float _S1665 = p_1.vignette_params_0[int(2)].alpha0_0 - mean_5;
    losses_3[int(12)] = (_S1663 * _S1663 + _S1664 * _S1664 + _S1665 * _S1665) / 3.0f;
    float mean_6 = (p_1.vignette_params_0[int(0)].alpha1_0 + p_1.vignette_params_0[int(1)].alpha1_0 + p_1.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S1666 = p_1.vignette_params_0[int(0)].alpha1_0 - mean_6;
    float _S1667 = p_1.vignette_params_0[int(1)].alpha1_0 - mean_6;
    float _S1668 = p_1.vignette_params_0[int(2)].alpha1_0 - mean_6;
    losses_3[int(13)] = (_S1666 * _S1666 + _S1667 * _S1667 + _S1668 * _S1668) / 3.0f;
    float mean_7 = (p_1.vignette_params_0[int(0)].alpha2_0 + p_1.vignette_params_0[int(1)].alpha2_0 + p_1.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S1669 = p_1.vignette_params_0[int(0)].alpha2_0 - mean_7;
    float _S1670 = p_1.vignette_params_0[int(1)].alpha2_0 - mean_7;
    float _S1671 = p_1.vignette_params_0[int(2)].alpha2_0 - mean_7;
    losses_3[int(14)] = (_S1669 * _S1669 + _S1670 * _S1670 + _S1671 * _S1671) / 3.0f;
    float2  bd_1 = mul_1(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_1(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_1(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_1(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    losses_3[int(15)] = bd_1.x;
    losses_3[int(16)] = bd_1.y;
    losses_3[int(17)] = rd_1.x;
    losses_3[int(18)] = rd_1.y;
    losses_3[int(19)] = gd_1.x;
    losses_3[int(20)] = gd_1.y;
    losses_3[int(21)] = nd_1.x;
    losses_3[int(22)] = nd_1.y;
    float mean_8 = (p_1.crf_params_0[int(0)].toe_0 + p_1.crf_params_0[int(1)].toe_0 + p_1.crf_params_0[int(2)].toe_0) / 3.0f;
    float _S1672 = p_1.crf_params_0[int(0)].toe_0 - mean_8;
    float _S1673 = p_1.crf_params_0[int(1)].toe_0 - mean_8;
    float _S1674 = p_1.crf_params_0[int(2)].toe_0 - mean_8;
    losses_3[int(23)] = (_S1672 * _S1672 + _S1673 * _S1673 + _S1674 * _S1674) / 3.0f;
    float mean_9 = (p_1.crf_params_0[int(0)].shoulder_0 + p_1.crf_params_0[int(1)].shoulder_0 + p_1.crf_params_0[int(2)].shoulder_0) / 3.0f;
    float _S1675 = p_1.crf_params_0[int(0)].shoulder_0 - mean_9;
    float _S1676 = p_1.crf_params_0[int(1)].shoulder_0 - mean_9;
    float _S1677 = p_1.crf_params_0[int(2)].shoulder_0 - mean_9;
    losses_3[int(24)] = (_S1675 * _S1675 + _S1676 * _S1676 + _S1677 * _S1677) / 3.0f;
    float mean_10 = (p_1.crf_params_0[int(0)].gamma_0 + p_1.crf_params_0[int(1)].gamma_0 + p_1.crf_params_0[int(2)].gamma_0) / 3.0f;
    float _S1678 = p_1.crf_params_0[int(0)].gamma_0 - mean_10;
    float _S1679 = p_1.crf_params_0[int(1)].gamma_0 - mean_10;
    float _S1680 = p_1.crf_params_0[int(2)].gamma_0 - mean_10;
    losses_3[int(25)] = (_S1678 * _S1678 + _S1679 * _S1679 + _S1680 * _S1680) / 3.0f;
    float mean_11 = (p_1.crf_params_0[int(0)].center_0 + p_1.crf_params_0[int(1)].center_0 + p_1.crf_params_0[int(2)].center_0) / 3.0f;
    float _S1681 = p_1.crf_params_0[int(0)].center_0 - mean_11;
    float _S1682 = p_1.crf_params_0[int(1)].center_0 - mean_11;
    float _S1683 = p_1.crf_params_0[int(2)].center_0 - mean_11;
    losses_3[int(26)] = (_S1681 * _S1681 + _S1682 * _S1682 + _S1683 * _S1683) / 3.0f;
    *_S1656 = losses_3;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_1, FixedArray<float, 27>  * _s_dOut_12)
{
    VignettingChannelParams_0 _S1684 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1685 = {
        _S1684, _S1684, _S1684
    };
    float2  _S1686 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1687 = { _S1686, _S1686, _S1686, _S1686 };
    CRFPPISPChannelParams_0 _S1688 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1689 = {
        _S1688, _S1688, _S1688
    };
    PPISPParams_0 _S1690;
    (&_S1690)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S1690)->vignette_params_0 = _S1685;
    (&_S1690)->color_params_0 = _S1687;
    (&_S1690)->crf_params_0 = _S1689;
    (&(&_S1690)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S1690)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    (&(&_S1690)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    (&(&_S1690)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    (&(&_S1690)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S1690)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S1690)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    (&(&_S1690)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    (&(&_S1690)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    (&(&_S1690)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S1690)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S1690)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    (&(&_S1690)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    (&(&_S1690)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    (&(&_S1690)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S1690)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S1690)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S1690)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S1690)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S1690)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S1690)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S1690)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S1690)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    (&(&_S1690)->crf_params_0[int(0)])->toe_0 = dpparams_1->primal_0[int(24)];
    (&(&_S1690)->crf_params_0[int(0)])->shoulder_0 = dpparams_1->primal_0[int(25)];
    (&(&_S1690)->crf_params_0[int(0)])->gamma_0 = dpparams_1->primal_0[int(26)];
    (&(&_S1690)->crf_params_0[int(0)])->center_0 = dpparams_1->primal_0[int(27)];
    (&(&_S1690)->crf_params_0[int(1)])->toe_0 = dpparams_1->primal_0[int(28)];
    (&(&_S1690)->crf_params_0[int(1)])->shoulder_0 = dpparams_1->primal_0[int(29)];
    (&(&_S1690)->crf_params_0[int(1)])->gamma_0 = dpparams_1->primal_0[int(30)];
    (&(&_S1690)->crf_params_0[int(1)])->center_0 = dpparams_1->primal_0[int(31)];
    (&(&_S1690)->crf_params_0[int(2)])->toe_0 = dpparams_1->primal_0[int(32)];
    (&(&_S1690)->crf_params_0[int(2)])->shoulder_0 = dpparams_1->primal_0[int(33)];
    (&(&_S1690)->crf_params_0[int(2)])->gamma_0 = dpparams_1->primal_0[int(34)];
    (&(&_S1690)->crf_params_0[int(2)])->center_0 = dpparams_1->primal_0[int(35)];
    float mean_12 = (dpparams_1->primal_0[int(1)] + dpparams_1->primal_0[int(6)] + dpparams_1->primal_0[int(11)]) / 3.0f;
    float _S1691 = dpparams_1->primal_0[int(1)] - mean_12;
    float _S1692 = dpparams_1->primal_0[int(6)] - mean_12;
    float _S1693 = dpparams_1->primal_0[int(11)] - mean_12;
    float mean_13 = (dpparams_1->primal_0[int(2)] + dpparams_1->primal_0[int(7)] + dpparams_1->primal_0[int(12)]) / 3.0f;
    float _S1694 = dpparams_1->primal_0[int(2)] - mean_13;
    float _S1695 = dpparams_1->primal_0[int(7)] - mean_13;
    float _S1696 = dpparams_1->primal_0[int(12)] - mean_13;
    float mean_14 = (dpparams_1->primal_0[int(3)] + dpparams_1->primal_0[int(8)] + dpparams_1->primal_0[int(13)]) / 3.0f;
    float _S1697 = dpparams_1->primal_0[int(3)] - mean_14;
    float _S1698 = dpparams_1->primal_0[int(8)] - mean_14;
    float _S1699 = dpparams_1->primal_0[int(13)] - mean_14;
    float mean_15 = (dpparams_1->primal_0[int(4)] + dpparams_1->primal_0[int(9)] + dpparams_1->primal_0[int(14)]) / 3.0f;
    float _S1700 = dpparams_1->primal_0[int(4)] - mean_15;
    float _S1701 = dpparams_1->primal_0[int(9)] - mean_15;
    float _S1702 = dpparams_1->primal_0[int(14)] - mean_15;
    float mean_16 = (dpparams_1->primal_0[int(5)] + dpparams_1->primal_0[int(10)] + dpparams_1->primal_0[int(15)]) / 3.0f;
    float _S1703 = dpparams_1->primal_0[int(5)] - mean_16;
    float _S1704 = dpparams_1->primal_0[int(10)] - mean_16;
    float _S1705 = dpparams_1->primal_0[int(15)] - mean_16;
    float mean_17 = (dpparams_1->primal_0[int(24)] + dpparams_1->primal_0[int(28)] + dpparams_1->primal_0[int(32)]) / 3.0f;
    float mean_18 = (dpparams_1->primal_0[int(25)] + dpparams_1->primal_0[int(29)] + dpparams_1->primal_0[int(33)]) / 3.0f;
    float mean_19 = (dpparams_1->primal_0[int(26)] + dpparams_1->primal_0[int(30)] + dpparams_1->primal_0[int(34)]) / 3.0f;
    float mean_20 = (dpparams_1->primal_0[int(27)] + dpparams_1->primal_0[int(31)] + dpparams_1->primal_0[int(35)]) / 3.0f;
    float _S1706 = 0.3333333432674408f * (*_s_dOut_12)[int(26)];
    float _S1707 = (dpparams_1->primal_0[int(35)] - mean_20) * _S1706;
    float _S1708 = _S1707 + _S1707;
    float _S1709 = (dpparams_1->primal_0[int(31)] - mean_20) * _S1706;
    float _S1710 = _S1709 + _S1709;
    float _S1711 = (dpparams_1->primal_0[int(27)] - mean_20) * _S1706;
    float _S1712 = _S1711 + _S1711;
    float _S1713 = 0.3333333432674408f * (- _S1708 + - _S1710 + - _S1712);
    float _S1714 = 0.3333333432674408f * (*_s_dOut_12)[int(25)];
    float _S1715 = (dpparams_1->primal_0[int(34)] - mean_19) * _S1714;
    float _S1716 = _S1715 + _S1715;
    float _S1717 = (dpparams_1->primal_0[int(30)] - mean_19) * _S1714;
    float _S1718 = _S1717 + _S1717;
    float _S1719 = (dpparams_1->primal_0[int(26)] - mean_19) * _S1714;
    float _S1720 = _S1719 + _S1719;
    float _S1721 = 0.3333333432674408f * (- _S1716 + - _S1718 + - _S1720);
    float _S1722 = 0.3333333432674408f * (*_s_dOut_12)[int(24)];
    float _S1723 = (dpparams_1->primal_0[int(33)] - mean_18) * _S1722;
    float _S1724 = _S1723 + _S1723;
    float _S1725 = (dpparams_1->primal_0[int(29)] - mean_18) * _S1722;
    float _S1726 = _S1725 + _S1725;
    float _S1727 = (dpparams_1->primal_0[int(25)] - mean_18) * _S1722;
    float _S1728 = _S1727 + _S1727;
    float _S1729 = 0.3333333432674408f * (- _S1724 + - _S1726 + - _S1728);
    float _S1730 = 0.3333333432674408f * (*_s_dOut_12)[int(23)];
    float _S1731 = (dpparams_1->primal_0[int(32)] - mean_17) * _S1730;
    float _S1732 = _S1731 + _S1731;
    float _S1733 = (dpparams_1->primal_0[int(28)] - mean_17) * _S1730;
    float _S1734 = _S1733 + _S1733;
    float _S1735 = (dpparams_1->primal_0[int(24)] - mean_17) * _S1730;
    float _S1736 = _S1735 + _S1735;
    float _S1737 = 0.3333333432674408f * (- _S1732 + - _S1734 + - _S1736);
    float2  _S1738 = make_float2 ((*_s_dOut_12)[int(21)], (*_s_dOut_12)[int(22)]);
    Matrix<float, 2, 2>  _S1739 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1740;
    (&_S1740)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1740)->differential_0 = _S1739;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1741;
    (&_S1741)->primal_0 = _S1690.color_params_0.n_0;
    (&_S1741)->differential_0 = _S1686;
    s_bwd_prop_mul_3(&_S1740, &_S1741, _S1738);
    float2  _S1742 = make_float2 ((*_s_dOut_12)[int(19)], (*_s_dOut_12)[int(20)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1743;
    (&_S1743)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1743)->differential_0 = _S1739;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1744;
    (&_S1744)->primal_0 = _S1690.color_params_0.g_0;
    (&_S1744)->differential_0 = _S1686;
    s_bwd_prop_mul_3(&_S1743, &_S1744, _S1742);
    float2  _S1745 = make_float2 ((*_s_dOut_12)[int(17)], (*_s_dOut_12)[int(18)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1746;
    (&_S1746)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1746)->differential_0 = _S1739;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1747;
    (&_S1747)->primal_0 = _S1690.color_params_0.r_0;
    (&_S1747)->differential_0 = _S1686;
    s_bwd_prop_mul_3(&_S1746, &_S1747, _S1745);
    float2  _S1748 = make_float2 ((*_s_dOut_12)[int(15)], (*_s_dOut_12)[int(16)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1749;
    (&_S1749)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1749)->differential_0 = _S1739;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1750;
    (&_S1750)->primal_0 = _S1690.color_params_0.b_0;
    (&_S1750)->differential_0 = _S1686;
    s_bwd_prop_mul_3(&_S1749, &_S1750, _S1748);
    ColorPPISPParams_0 _S1751 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1751)->n_0 = _S1741.differential_0;
    (&_S1751)->g_0 = _S1744.differential_0;
    (&_S1751)->r_0 = _S1747.differential_0;
    (&_S1751)->b_0 = _S1750.differential_0;
    float _S1752 = 0.3333333432674408f * (*_s_dOut_12)[int(14)];
    float _S1753 = _S1705 * _S1752;
    float _S1754 = _S1753 + _S1753;
    float _S1755 = _S1704 * _S1752;
    float _S1756 = _S1755 + _S1755;
    float _S1757 = _S1703 * _S1752;
    float _S1758 = _S1757 + _S1757;
    float _S1759 = 0.3333333432674408f * (- _S1754 + - _S1756 + - _S1758);
    float _S1760 = 0.3333333432674408f * (*_s_dOut_12)[int(13)];
    float _S1761 = _S1702 * _S1760;
    float _S1762 = _S1761 + _S1761;
    float _S1763 = _S1701 * _S1760;
    float _S1764 = _S1763 + _S1763;
    float _S1765 = _S1700 * _S1760;
    float _S1766 = _S1765 + _S1765;
    float _S1767 = 0.3333333432674408f * (- _S1762 + - _S1764 + - _S1766);
    float _S1768 = 0.3333333432674408f * (*_s_dOut_12)[int(12)];
    float _S1769 = _S1699 * _S1768;
    float _S1770 = _S1769 + _S1769;
    float _S1771 = _S1698 * _S1768;
    float _S1772 = _S1771 + _S1771;
    float _S1773 = _S1697 * _S1768;
    float _S1774 = _S1773 + _S1773;
    float _S1775 = 0.3333333432674408f * (- _S1770 + - _S1772 + - _S1774);
    float _S1776 = 0.3333333432674408f * (*_s_dOut_12)[int(11)];
    float _S1777 = _S1696 * _S1776;
    float _S1778 = _S1777 + _S1777;
    float _S1779 = _S1695 * _S1776;
    float _S1780 = _S1779 + _S1779;
    float _S1781 = _S1694 * _S1776;
    float _S1782 = _S1781 + _S1781;
    float _S1783 = 0.3333333432674408f * (- _S1778 + - _S1780 + - _S1782);
    float _S1784 = 0.3333333432674408f * (*_s_dOut_12)[int(10)];
    float _S1785 = _S1693 * _S1784;
    float _S1786 = _S1785 + _S1785;
    float _S1787 = _S1692 * _S1784;
    float _S1788 = _S1787 + _S1787;
    float _S1789 = _S1691 * _S1784;
    float _S1790 = _S1789 + _S1789;
    float _S1791 = 0.3333333432674408f * (- _S1786 + - _S1788 + - _S1790);
    DiffPair_float_0 _S1792;
    (&_S1792)->primal_0 = 0.0f;
    (&_S1792)->differential_0 = 0.0f;
    DiffPair_float_0 _S1793;
    (&_S1793)->primal_0 = dpparams_1->primal_0[int(15)];
    (&_S1793)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1792, &_S1793, (*_s_dOut_12)[int(9)]);
    DiffPair_float_0 _S1794;
    (&_S1794)->primal_0 = 0.0f;
    (&_S1794)->differential_0 = 0.0f;
    DiffPair_float_0 _S1795;
    (&_S1795)->primal_0 = dpparams_1->primal_0[int(10)];
    (&_S1795)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1794, &_S1795, (*_s_dOut_12)[int(9)]);
    DiffPair_float_0 _S1796;
    (&_S1796)->primal_0 = 0.0f;
    (&_S1796)->differential_0 = 0.0f;
    DiffPair_float_0 _S1797;
    (&_S1797)->primal_0 = dpparams_1->primal_0[int(5)];
    (&_S1797)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1796, &_S1797, (*_s_dOut_12)[int(9)]);
    DiffPair_float_0 _S1798;
    (&_S1798)->primal_0 = 0.0f;
    (&_S1798)->differential_0 = 0.0f;
    DiffPair_float_0 _S1799;
    (&_S1799)->primal_0 = dpparams_1->primal_0[int(14)];
    (&_S1799)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1798, &_S1799, (*_s_dOut_12)[int(8)]);
    DiffPair_float_0 _S1800;
    (&_S1800)->primal_0 = 0.0f;
    (&_S1800)->differential_0 = 0.0f;
    DiffPair_float_0 _S1801;
    (&_S1801)->primal_0 = dpparams_1->primal_0[int(9)];
    (&_S1801)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1800, &_S1801, (*_s_dOut_12)[int(8)]);
    DiffPair_float_0 _S1802;
    (&_S1802)->primal_0 = 0.0f;
    (&_S1802)->differential_0 = 0.0f;
    DiffPair_float_0 _S1803;
    (&_S1803)->primal_0 = dpparams_1->primal_0[int(4)];
    (&_S1803)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1802, &_S1803, (*_s_dOut_12)[int(8)]);
    DiffPair_float_0 _S1804;
    (&_S1804)->primal_0 = 0.0f;
    (&_S1804)->differential_0 = 0.0f;
    DiffPair_float_0 _S1805;
    (&_S1805)->primal_0 = dpparams_1->primal_0[int(13)];
    (&_S1805)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1804, &_S1805, (*_s_dOut_12)[int(7)]);
    DiffPair_float_0 _S1806;
    (&_S1806)->primal_0 = 0.0f;
    (&_S1806)->differential_0 = 0.0f;
    DiffPair_float_0 _S1807;
    (&_S1807)->primal_0 = dpparams_1->primal_0[int(8)];
    (&_S1807)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1806, &_S1807, (*_s_dOut_12)[int(7)]);
    DiffPair_float_0 _S1808;
    (&_S1808)->primal_0 = 0.0f;
    (&_S1808)->differential_0 = 0.0f;
    DiffPair_float_0 _S1809;
    (&_S1809)->primal_0 = dpparams_1->primal_0[int(3)];
    (&_S1809)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1808, &_S1809, (*_s_dOut_12)[int(7)]);
    PPISPParams_0 _S1810 = PPISPParams_x24_syn_dzero_0();
    (&_S1810)->color_params_0 = _S1751;
    (&_S1810)->exposure_0 = (*_s_dOut_12)[int(0)];
    _S1690 = _S1810;
    (&(&_S1690)->crf_params_0[int(2)])->center_0 = 0.0f;
    float _S1811 = _S1708 + _S1713 + _S1810.crf_params_0[int(2)].center_0;
    (&(&_S1690)->crf_params_0[int(2)])->gamma_0 = 0.0f;
    float _S1812 = _S1716 + _S1721 + _S1810.crf_params_0[int(2)].gamma_0;
    (&(&_S1690)->crf_params_0[int(2)])->shoulder_0 = 0.0f;
    float _S1813 = _S1724 + _S1729 + _S1810.crf_params_0[int(2)].shoulder_0;
    (&(&_S1690)->crf_params_0[int(2)])->toe_0 = 0.0f;
    float _S1814 = _S1732 + _S1737 + _S1810.crf_params_0[int(2)].toe_0;
    (&(&_S1690)->crf_params_0[int(1)])->center_0 = 0.0f;
    float _S1815 = _S1710 + _S1713 + _S1810.crf_params_0[int(1)].center_0;
    (&(&_S1690)->crf_params_0[int(1)])->gamma_0 = 0.0f;
    float _S1816 = _S1718 + _S1721 + _S1810.crf_params_0[int(1)].gamma_0;
    (&(&_S1690)->crf_params_0[int(1)])->shoulder_0 = 0.0f;
    float _S1817 = _S1726 + _S1729 + _S1810.crf_params_0[int(1)].shoulder_0;
    (&(&_S1690)->crf_params_0[int(1)])->toe_0 = 0.0f;
    float _S1818 = _S1734 + _S1737 + _S1810.crf_params_0[int(1)].toe_0;
    (&(&_S1690)->crf_params_0[int(0)])->center_0 = 0.0f;
    float _S1819 = _S1712 + _S1713 + _S1810.crf_params_0[int(0)].center_0;
    (&(&_S1690)->crf_params_0[int(0)])->gamma_0 = 0.0f;
    float _S1820 = _S1720 + _S1721 + _S1810.crf_params_0[int(0)].gamma_0;
    (&(&_S1690)->crf_params_0[int(0)])->shoulder_0 = 0.0f;
    float _S1821 = _S1728 + _S1729 + _S1810.crf_params_0[int(0)].shoulder_0;
    (&(&_S1690)->crf_params_0[int(0)])->toe_0 = 0.0f;
    float _S1822 = _S1736 + _S1737 + _S1810.crf_params_0[int(0)].toe_0;
    *&((&(&(&_S1690)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1690)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1690)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1823 = _S1754 + _S1759 + _S1793.differential_0 + _S1810.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1690)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1824 = _S1762 + _S1767 + _S1799.differential_0 + _S1810.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1690)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1825 = _S1770 + _S1775 + _S1805.differential_0 + _S1810.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1690)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1826 = (*_s_dOut_12)[int(6)] + _S1778 + _S1783 + _S1810.vignette_params_0[int(2)].cy_0;
    (&(&_S1690)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1827 = (*_s_dOut_12)[int(5)] + _S1786 + _S1791 + _S1810.vignette_params_0[int(2)].cx_0;
    (&(&_S1690)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1828 = _S1756 + _S1759 + _S1795.differential_0 + _S1810.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1690)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1829 = _S1764 + _S1767 + _S1801.differential_0 + _S1810.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1690)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1830 = _S1772 + _S1775 + _S1807.differential_0 + _S1810.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1690)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1831 = (*_s_dOut_12)[int(4)] + _S1780 + _S1783 + _S1810.vignette_params_0[int(1)].cy_0;
    (&(&_S1690)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1832 = (*_s_dOut_12)[int(3)] + _S1788 + _S1791 + _S1810.vignette_params_0[int(1)].cx_0;
    (&(&_S1690)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1833 = _S1758 + _S1759 + _S1797.differential_0 + _S1810.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1690)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1834 = _S1766 + _S1767 + _S1803.differential_0 + _S1810.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1690)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1835 = _S1774 + _S1775 + _S1809.differential_0 + _S1810.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1690)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1836 = (*_s_dOut_12)[int(2)] + _S1782 + _S1783 + _S1810.vignette_params_0[int(0)].cy_0;
    (&(&_S1690)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1837 = (*_s_dOut_12)[int(1)] + _S1790 + _S1791 + _S1810.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 36>  _S1838;
    _S1838[int(0)] = 0.0f;
    _S1838[int(1)] = 0.0f;
    _S1838[int(2)] = 0.0f;
    _S1838[int(3)] = 0.0f;
    _S1838[int(4)] = 0.0f;
    _S1838[int(5)] = 0.0f;
    _S1838[int(6)] = 0.0f;
    _S1838[int(7)] = 0.0f;
    _S1838[int(8)] = 0.0f;
    _S1838[int(9)] = 0.0f;
    _S1838[int(10)] = 0.0f;
    _S1838[int(11)] = 0.0f;
    _S1838[int(12)] = 0.0f;
    _S1838[int(13)] = 0.0f;
    _S1838[int(14)] = 0.0f;
    _S1838[int(15)] = 0.0f;
    _S1838[int(16)] = 0.0f;
    _S1838[int(17)] = 0.0f;
    _S1838[int(18)] = 0.0f;
    _S1838[int(19)] = 0.0f;
    _S1838[int(20)] = 0.0f;
    _S1838[int(21)] = 0.0f;
    _S1838[int(22)] = 0.0f;
    _S1838[int(23)] = 0.0f;
    _S1838[int(24)] = 0.0f;
    _S1838[int(25)] = 0.0f;
    _S1838[int(26)] = 0.0f;
    _S1838[int(27)] = 0.0f;
    _S1838[int(28)] = 0.0f;
    _S1838[int(29)] = 0.0f;
    _S1838[int(30)] = 0.0f;
    _S1838[int(31)] = 0.0f;
    _S1838[int(32)] = 0.0f;
    _S1838[int(33)] = 0.0f;
    _S1838[int(34)] = 0.0f;
    _S1838[int(35)] = 0.0f;
    _S1838[int(8)] = _S1830;
    _S1838[int(16)] = _S1810.color_params_0.b_0.x;
    _S1838[int(15)] = _S1823;
    _S1838[int(14)] = _S1824;
    _S1838[int(13)] = _S1825;
    _S1838[int(12)] = _S1826;
    _S1838[int(11)] = _S1827;
    _S1838[int(10)] = _S1828;
    _S1838[int(9)] = _S1829;
    _S1838[int(17)] = _S1810.color_params_0.b_0.y;
    _S1838[int(7)] = _S1831;
    _S1838[int(6)] = _S1832;
    _S1838[int(5)] = _S1833;
    _S1838[int(4)] = _S1834;
    _S1838[int(3)] = _S1835;
    _S1838[int(2)] = _S1836;
    _S1838[int(1)] = _S1837;
    _S1838[int(0)] = _S1690.exposure_0;
    _S1838[int(26)] = _S1820;
    _S1838[int(34)] = _S1812;
    _S1838[int(33)] = _S1813;
    _S1838[int(32)] = _S1814;
    _S1838[int(31)] = _S1815;
    _S1838[int(30)] = _S1816;
    _S1838[int(29)] = _S1817;
    _S1838[int(28)] = _S1818;
    _S1838[int(27)] = _S1819;
    _S1838[int(35)] = _S1811;
    _S1838[int(25)] = _S1821;
    _S1838[int(24)] = _S1822;
    _S1838[int(23)] = _S1810.color_params_0.n_0.y;
    _S1838[int(22)] = _S1810.color_params_0.n_0.x;
    _S1838[int(21)] = _S1810.color_params_0.g_0.y;
    _S1838[int(20)] = _S1810.color_params_0.g_0.x;
    _S1838[int(19)] = _S1810.color_params_0.r_0.y;
    _S1838[int(18)] = _S1810.color_params_0.r_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S1838;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1839, FixedArray<float, 27>  * _S1840)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S1839, _S1840);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  * params_3, FixedArray<float, 27>  * grad_out_1, FixedArray<float, 36>  * _S1841)
{
    FixedArray<float, 36>  _S1842 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = *params_3;
    (&dp_params_1)->differential_0 = _S1842;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_1, grad_out_1);
    *_S1841 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 27>  * raw_losses_2, int num_cameras_0, FixedArray<float, 6>  * loss_weights_0, FixedArray<float, 6>  * _S1843)
{
    float _S1844;
    FixedArray<float, 6>  losses_4;
    float _S1845 = float(num_cameras_0);
    float _S1846 = (*raw_losses_2)[int(0)] / _S1845;
    for(;;)
    {
        float _S1847 = (F32_abs((_S1846)));
        if(_S1847 < 0.10000000149011612f)
        {
            _S1844 = 0.5f * _S1846 * _S1846 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1844 = _S1847 - 0.05000000074505806f;
            break;
        }
    }
    losses_4[int(0)] = _S1844;
    float3  mean_cx_0 = make_float3 ((*raw_losses_2)[int(1)], (*raw_losses_2)[int(3)], (*raw_losses_2)[int(5)]) / make_float3 (_S1845);
    float3  mean_cy_0 = make_float3 ((*raw_losses_2)[int(2)], (*raw_losses_2)[int(4)], (*raw_losses_2)[int(6)]) / make_float3 (_S1845);
    losses_4[int(1)] = (dot_0(mean_cx_0, mean_cx_0) + dot_0(mean_cy_0, mean_cy_0)) / 3.0f;
    losses_4[int(2)] = ((*raw_losses_2)[int(7)] + (*raw_losses_2)[int(8)] + (*raw_losses_2)[int(9)]) / (3.0f * _S1845);
    losses_4[int(3)] = ((*raw_losses_2)[int(10)] + (*raw_losses_2)[int(11)] + (*raw_losses_2)[int(12)] + (*raw_losses_2)[int(13)] + (*raw_losses_2)[int(14)]) / (5.0f * _S1845);
    float _S1848 = (*raw_losses_2)[int(15)] / _S1845;
    for(;;)
    {
        float _S1849 = (F32_abs((_S1848)));
        if(_S1849 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1848 * _S1848 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1849 - 0.00249999994412065f;
            break;
        }
    }
    float _S1850;
    float _S1851 = (*raw_losses_2)[int(16)] / _S1845;
    for(;;)
    {
        float _S1852 = (F32_abs((_S1851)));
        if(_S1852 < 0.00499999988824129f)
        {
            _S1850 = 0.5f * _S1851 * _S1851 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1850 = _S1852 - 0.00249999994412065f;
            break;
        }
    }
    float _S1853 = _S1844 + _S1850;
    float _S1854 = (*raw_losses_2)[int(17)] / _S1845;
    for(;;)
    {
        float _S1855 = (F32_abs((_S1854)));
        if(_S1855 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1854 * _S1854 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1855 - 0.00249999994412065f;
            break;
        }
    }
    float _S1856 = _S1853 + _S1844;
    float _S1857 = (*raw_losses_2)[int(18)] / _S1845;
    for(;;)
    {
        float _S1858 = (F32_abs((_S1857)));
        if(_S1858 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1857 * _S1857 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1858 - 0.00249999994412065f;
            break;
        }
    }
    float _S1859 = _S1856 + _S1844;
    float _S1860 = (*raw_losses_2)[int(19)] / _S1845;
    for(;;)
    {
        float _S1861 = (F32_abs((_S1860)));
        if(_S1861 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1860 * _S1860 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1861 - 0.00249999994412065f;
            break;
        }
    }
    float _S1862 = _S1859 + _S1844;
    float _S1863 = (*raw_losses_2)[int(20)] / _S1845;
    for(;;)
    {
        float _S1864 = (F32_abs((_S1863)));
        if(_S1864 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1863 * _S1863 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1864 - 0.00249999994412065f;
            break;
        }
    }
    float _S1865 = _S1862 + _S1844;
    float _S1866 = (*raw_losses_2)[int(21)] / _S1845;
    for(;;)
    {
        float _S1867 = (F32_abs((_S1866)));
        if(_S1867 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1866 * _S1866 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1867 - 0.00249999994412065f;
            break;
        }
    }
    float _S1868 = _S1865 + _S1844;
    float _S1869 = (*raw_losses_2)[int(22)] / _S1845;
    for(;;)
    {
        float _S1870 = (F32_abs((_S1869)));
        if(_S1870 < 0.00499999988824129f)
        {
            _S1844 = 0.5f * _S1869 * _S1869 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1844 = _S1870 - 0.00249999994412065f;
            break;
        }
    }
    float _S1871 = (_S1868 + _S1844) / 8.0f;
    float _S1872 = ((*raw_losses_2)[int(23)] + (*raw_losses_2)[int(24)] + (*raw_losses_2)[int(25)] + (*raw_losses_2)[int(26)]) / (4.0f * _S1845);
    losses_4[int(0)] = losses_4[int(0)] * (*loss_weights_0)[int(0)];
    losses_4[int(1)] = losses_4[int(1)] * (*loss_weights_0)[int(1)];
    losses_4[int(2)] = losses_4[int(2)] * (*loss_weights_0)[int(2)];
    losses_4[int(3)] = losses_4[int(3)] * (*loss_weights_0)[int(3)];
    losses_4[int(4)] = _S1871 * (*loss_weights_0)[int(4)];
    losses_4[int(5)] = _S1872 * (*loss_weights_0)[int(5)];
    *_S1843 = losses_4;
    return;
}

struct DiffPair_arrayx3Cfloatx2C27x3E_0
{
    FixedArray<float, 27>  primal_0;
    FixedArray<float, 27>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C27x3E_0 * dpraw_losses_1, int num_cameras_1, FixedArray<float, 6>  * loss_weights_1, FixedArray<float, 6>  * _s_dOut_13)
{
    FixedArray<float, 27>  _S1873 = dpraw_losses_1->primal_0;
    float _S1874 = float(num_cameras_1);
    float3  _S1875 = make_float3 (_S1874);
    float _S1876 = dpraw_losses_1->primal_0[int(0)] / _S1874;
    bool _S1877 = (s_primal_ctx_abs_0(_S1876)) < 0.10000000149011612f;
    float _S1878;
    if(_S1877)
    {
        _S1878 = 0.5f * _S1876;
    }
    else
    {
        _S1878 = 0.0f;
    }
    float3  mean_cx_1 = make_float3 (_S1873[int(1)], _S1873[int(3)], _S1873[int(5)]) / make_float3 (_S1874);
    float3  mean_cy_1 = make_float3 (_S1873[int(2)], _S1873[int(4)], _S1873[int(6)]) / make_float3 (_S1874);
    float _S1879 = 3.0f * _S1874;
    float _S1880 = 5.0f * _S1874;
    float _S1881 = _S1873[int(15)] / _S1874;
    bool _S1882 = (s_primal_ctx_abs_0(_S1881)) < 0.00499999988824129f;
    float _S1883;
    if(_S1882)
    {
        _S1883 = 0.5f * _S1881;
    }
    else
    {
        _S1883 = 0.0f;
    }
    float _S1884 = _S1873[int(16)] / _S1874;
    bool _S1885 = (s_primal_ctx_abs_0(_S1884)) < 0.00499999988824129f;
    float _S1886;
    if(_S1885)
    {
        _S1886 = 0.5f * _S1884;
    }
    else
    {
        _S1886 = 0.0f;
    }
    float _S1887 = _S1873[int(17)] / _S1874;
    bool _S1888 = (s_primal_ctx_abs_0(_S1887)) < 0.00499999988824129f;
    float _S1889;
    if(_S1888)
    {
        _S1889 = 0.5f * _S1887;
    }
    else
    {
        _S1889 = 0.0f;
    }
    float _S1890 = _S1873[int(18)] / _S1874;
    bool _S1891 = (s_primal_ctx_abs_0(_S1890)) < 0.00499999988824129f;
    float _S1892;
    if(_S1891)
    {
        _S1892 = 0.5f * _S1890;
    }
    else
    {
        _S1892 = 0.0f;
    }
    float _S1893 = _S1873[int(19)] / _S1874;
    bool _S1894 = (s_primal_ctx_abs_0(_S1893)) < 0.00499999988824129f;
    float _S1895;
    if(_S1894)
    {
        _S1895 = 0.5f * _S1893;
    }
    else
    {
        _S1895 = 0.0f;
    }
    float _S1896 = _S1873[int(20)] / _S1874;
    bool _S1897 = (s_primal_ctx_abs_0(_S1896)) < 0.00499999988824129f;
    float _S1898;
    if(_S1897)
    {
        _S1898 = 0.5f * _S1896;
    }
    else
    {
        _S1898 = 0.0f;
    }
    float _S1899 = _S1873[int(21)] / _S1874;
    bool _S1900 = (s_primal_ctx_abs_0(_S1899)) < 0.00499999988824129f;
    float _S1901;
    if(_S1900)
    {
        _S1901 = 0.5f * _S1899;
    }
    else
    {
        _S1901 = 0.0f;
    }
    float _S1902 = _S1873[int(22)] / _S1874;
    bool _S1903 = (s_primal_ctx_abs_0(_S1902)) < 0.00499999988824129f;
    float _S1904;
    if(_S1903)
    {
        _S1904 = 0.5f * _S1902;
    }
    else
    {
        _S1904 = 0.0f;
    }
    float _S1905 = (*loss_weights_1)[int(3)] * (*_s_dOut_13)[int(3)];
    float _S1906 = (*loss_weights_1)[int(2)] * (*_s_dOut_13)[int(2)];
    float _S1907 = (*loss_weights_1)[int(1)] * (*_s_dOut_13)[int(1)];
    float _S1908 = (*loss_weights_1)[int(0)] * (*_s_dOut_13)[int(0)];
    float _S1909 = (*loss_weights_1)[int(5)] * (*_s_dOut_13)[int(5)] / (4.0f * _S1874);
    float _S1910 = 0.125f * ((*loss_weights_1)[int(4)] * (*_s_dOut_13)[int(4)]);
    FixedArray<float, 27>  _S1911;
    _S1911[int(0)] = 0.0f;
    _S1911[int(1)] = 0.0f;
    _S1911[int(2)] = 0.0f;
    _S1911[int(3)] = 0.0f;
    _S1911[int(4)] = 0.0f;
    _S1911[int(5)] = 0.0f;
    _S1911[int(6)] = 0.0f;
    _S1911[int(7)] = 0.0f;
    _S1911[int(8)] = 0.0f;
    _S1911[int(9)] = 0.0f;
    _S1911[int(10)] = 0.0f;
    _S1911[int(11)] = 0.0f;
    _S1911[int(12)] = 0.0f;
    _S1911[int(13)] = 0.0f;
    _S1911[int(14)] = 0.0f;
    _S1911[int(15)] = 0.0f;
    _S1911[int(16)] = 0.0f;
    _S1911[int(17)] = 0.0f;
    _S1911[int(18)] = 0.0f;
    _S1911[int(19)] = 0.0f;
    _S1911[int(20)] = 0.0f;
    _S1911[int(21)] = 0.0f;
    _S1911[int(22)] = 0.0f;
    _S1911[int(23)] = 0.0f;
    _S1911[int(24)] = 0.0f;
    _S1911[int(25)] = 0.0f;
    _S1911[int(26)] = 0.0f;
    _S1911[int(26)] = _S1909;
    _S1911[int(25)] = _S1909;
    _S1911[int(24)] = _S1909;
    _S1911[int(23)] = _S1909;
    float _S1912 = _S1911[int(0)];
    float _S1913 = _S1911[int(1)];
    float _S1914 = _S1911[int(2)];
    float _S1915 = _S1911[int(3)];
    float _S1916 = _S1911[int(4)];
    float _S1917 = _S1911[int(5)];
    float _S1918 = _S1911[int(6)];
    float _S1919 = _S1911[int(7)];
    float _S1920 = _S1911[int(8)];
    float _S1921 = _S1911[int(9)];
    float _S1922 = _S1911[int(10)];
    float _S1923 = _S1911[int(11)];
    float _S1924 = _S1911[int(12)];
    float _S1925 = _S1911[int(13)];
    float _S1926 = _S1911[int(14)];
    float _S1927 = _S1911[int(15)];
    float _S1928 = _S1911[int(16)];
    float _S1929 = _S1911[int(17)];
    float _S1930 = _S1911[int(18)];
    float _S1931 = _S1911[int(19)];
    float _S1932 = _S1911[int(20)];
    float _S1933 = _S1911[int(21)];
    float _S1934 = _S1911[int(22)];
    float _S1935 = _S1911[int(23)];
    float _S1936 = _S1911[int(24)];
    float _S1937 = _S1911[int(25)];
    float _S1938 = _S1911[int(26)];
    float _S1939;
    if(_S1903)
    {
        float _S1940 = 200.0f * _S1910;
        float _S1941 = _S1904 * _S1940 + 0.5f * (_S1902 * _S1940);
        _S1904 = 0.0f;
        _S1939 = _S1941;
    }
    else
    {
        _S1904 = _S1910;
        _S1939 = 0.0f;
    }
    DiffPair_float_0 _S1942;
    (&_S1942)->primal_0 = _S1902;
    (&_S1942)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1942, _S1904);
    float _S1943 = (_S1942.differential_0 + _S1939) / _S1874;
    FixedArray<float, 27>  _S1944;
    _S1944[int(0)] = 0.0f;
    _S1944[int(1)] = 0.0f;
    _S1944[int(2)] = 0.0f;
    _S1944[int(3)] = 0.0f;
    _S1944[int(4)] = 0.0f;
    _S1944[int(5)] = 0.0f;
    _S1944[int(6)] = 0.0f;
    _S1944[int(7)] = 0.0f;
    _S1944[int(8)] = 0.0f;
    _S1944[int(9)] = 0.0f;
    _S1944[int(10)] = 0.0f;
    _S1944[int(11)] = 0.0f;
    _S1944[int(12)] = 0.0f;
    _S1944[int(13)] = 0.0f;
    _S1944[int(14)] = 0.0f;
    _S1944[int(15)] = 0.0f;
    _S1944[int(16)] = 0.0f;
    _S1944[int(17)] = 0.0f;
    _S1944[int(18)] = 0.0f;
    _S1944[int(19)] = 0.0f;
    _S1944[int(20)] = 0.0f;
    _S1944[int(21)] = 0.0f;
    _S1944[int(22)] = 0.0f;
    _S1944[int(23)] = 0.0f;
    _S1944[int(24)] = 0.0f;
    _S1944[int(25)] = 0.0f;
    _S1944[int(26)] = 0.0f;
    _S1944[int(22)] = _S1943;
    float _S1945 = _S1912 + _S1944[int(0)];
    float _S1946 = _S1913 + _S1944[int(1)];
    float _S1947 = _S1914 + _S1944[int(2)];
    float _S1948 = _S1915 + _S1944[int(3)];
    float _S1949 = _S1916 + _S1944[int(4)];
    float _S1950 = _S1917 + _S1944[int(5)];
    float _S1951 = _S1918 + _S1944[int(6)];
    float _S1952 = _S1919 + _S1944[int(7)];
    float _S1953 = _S1920 + _S1944[int(8)];
    float _S1954 = _S1921 + _S1944[int(9)];
    float _S1955 = _S1922 + _S1944[int(10)];
    float _S1956 = _S1923 + _S1944[int(11)];
    float _S1957 = _S1924 + _S1944[int(12)];
    float _S1958 = _S1925 + _S1944[int(13)];
    float _S1959 = _S1926 + _S1944[int(14)];
    float _S1960 = _S1927 + _S1944[int(15)];
    float _S1961 = _S1928 + _S1944[int(16)];
    float _S1962 = _S1929 + _S1944[int(17)];
    float _S1963 = _S1930 + _S1944[int(18)];
    float _S1964 = _S1931 + _S1944[int(19)];
    float _S1965 = _S1932 + _S1944[int(20)];
    float _S1966 = _S1933 + _S1944[int(21)];
    float _S1967 = _S1934 + _S1944[int(22)];
    float _S1968 = _S1935 + _S1944[int(23)];
    float _S1969 = _S1936 + _S1944[int(24)];
    float _S1970 = _S1937 + _S1944[int(25)];
    float _S1971 = _S1938 + _S1944[int(26)];
    if(_S1900)
    {
        float _S1972 = 200.0f * _S1910;
        float _S1973 = _S1901 * _S1972 + 0.5f * (_S1899 * _S1972);
        _S1901 = 0.0f;
        _S1904 = _S1973;
    }
    else
    {
        _S1901 = _S1910;
        _S1904 = 0.0f;
    }
    DiffPair_float_0 _S1974;
    (&_S1974)->primal_0 = _S1899;
    (&_S1974)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1974, _S1901);
    float _S1975 = (_S1974.differential_0 + _S1904) / _S1874;
    FixedArray<float, 27>  _S1976;
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
    _S1976[int(22)] = 0.0f;
    _S1976[int(23)] = 0.0f;
    _S1976[int(24)] = 0.0f;
    _S1976[int(25)] = 0.0f;
    _S1976[int(26)] = 0.0f;
    _S1976[int(21)] = _S1975;
    float _S1977 = _S1945 + _S1976[int(0)];
    float _S1978 = _S1946 + _S1976[int(1)];
    float _S1979 = _S1947 + _S1976[int(2)];
    float _S1980 = _S1948 + _S1976[int(3)];
    float _S1981 = _S1949 + _S1976[int(4)];
    float _S1982 = _S1950 + _S1976[int(5)];
    float _S1983 = _S1951 + _S1976[int(6)];
    float _S1984 = _S1952 + _S1976[int(7)];
    float _S1985 = _S1953 + _S1976[int(8)];
    float _S1986 = _S1954 + _S1976[int(9)];
    float _S1987 = _S1955 + _S1976[int(10)];
    float _S1988 = _S1956 + _S1976[int(11)];
    float _S1989 = _S1957 + _S1976[int(12)];
    float _S1990 = _S1958 + _S1976[int(13)];
    float _S1991 = _S1959 + _S1976[int(14)];
    float _S1992 = _S1960 + _S1976[int(15)];
    float _S1993 = _S1961 + _S1976[int(16)];
    float _S1994 = _S1962 + _S1976[int(17)];
    float _S1995 = _S1963 + _S1976[int(18)];
    float _S1996 = _S1964 + _S1976[int(19)];
    float _S1997 = _S1965 + _S1976[int(20)];
    float _S1998 = _S1966 + _S1976[int(21)];
    float _S1999 = _S1967 + _S1976[int(22)];
    float _S2000 = _S1968 + _S1976[int(23)];
    float _S2001 = _S1969 + _S1976[int(24)];
    float _S2002 = _S1970 + _S1976[int(25)];
    float _S2003 = _S1971 + _S1976[int(26)];
    if(_S1897)
    {
        float _S2004 = 200.0f * _S1910;
        float _S2005 = _S1898 * _S2004 + 0.5f * (_S1896 * _S2004);
        _S1898 = 0.0f;
        _S1901 = _S2005;
    }
    else
    {
        _S1898 = _S1910;
        _S1901 = 0.0f;
    }
    DiffPair_float_0 _S2006;
    (&_S2006)->primal_0 = _S1896;
    (&_S2006)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2006, _S1898);
    float _S2007 = (_S2006.differential_0 + _S1901) / _S1874;
    FixedArray<float, 27>  _S2008;
    _S2008[int(0)] = 0.0f;
    _S2008[int(1)] = 0.0f;
    _S2008[int(2)] = 0.0f;
    _S2008[int(3)] = 0.0f;
    _S2008[int(4)] = 0.0f;
    _S2008[int(5)] = 0.0f;
    _S2008[int(6)] = 0.0f;
    _S2008[int(7)] = 0.0f;
    _S2008[int(8)] = 0.0f;
    _S2008[int(9)] = 0.0f;
    _S2008[int(10)] = 0.0f;
    _S2008[int(11)] = 0.0f;
    _S2008[int(12)] = 0.0f;
    _S2008[int(13)] = 0.0f;
    _S2008[int(14)] = 0.0f;
    _S2008[int(15)] = 0.0f;
    _S2008[int(16)] = 0.0f;
    _S2008[int(17)] = 0.0f;
    _S2008[int(18)] = 0.0f;
    _S2008[int(19)] = 0.0f;
    _S2008[int(20)] = 0.0f;
    _S2008[int(21)] = 0.0f;
    _S2008[int(22)] = 0.0f;
    _S2008[int(23)] = 0.0f;
    _S2008[int(24)] = 0.0f;
    _S2008[int(25)] = 0.0f;
    _S2008[int(26)] = 0.0f;
    _S2008[int(20)] = _S2007;
    float _S2009 = _S1977 + _S2008[int(0)];
    float _S2010 = _S1978 + _S2008[int(1)];
    float _S2011 = _S1979 + _S2008[int(2)];
    float _S2012 = _S1980 + _S2008[int(3)];
    float _S2013 = _S1981 + _S2008[int(4)];
    float _S2014 = _S1982 + _S2008[int(5)];
    float _S2015 = _S1983 + _S2008[int(6)];
    float _S2016 = _S1984 + _S2008[int(7)];
    float _S2017 = _S1985 + _S2008[int(8)];
    float _S2018 = _S1986 + _S2008[int(9)];
    float _S2019 = _S1987 + _S2008[int(10)];
    float _S2020 = _S1988 + _S2008[int(11)];
    float _S2021 = _S1989 + _S2008[int(12)];
    float _S2022 = _S1990 + _S2008[int(13)];
    float _S2023 = _S1991 + _S2008[int(14)];
    float _S2024 = _S1992 + _S2008[int(15)];
    float _S2025 = _S1993 + _S2008[int(16)];
    float _S2026 = _S1994 + _S2008[int(17)];
    float _S2027 = _S1995 + _S2008[int(18)];
    float _S2028 = _S1996 + _S2008[int(19)];
    float _S2029 = _S1997 + _S2008[int(20)];
    float _S2030 = _S1998 + _S2008[int(21)];
    float _S2031 = _S1999 + _S2008[int(22)];
    float _S2032 = _S2000 + _S2008[int(23)];
    float _S2033 = _S2001 + _S2008[int(24)];
    float _S2034 = _S2002 + _S2008[int(25)];
    float _S2035 = _S2003 + _S2008[int(26)];
    if(_S1894)
    {
        float _S2036 = 200.0f * _S1910;
        float _S2037 = _S1895 * _S2036 + 0.5f * (_S1893 * _S2036);
        _S1895 = 0.0f;
        _S1898 = _S2037;
    }
    else
    {
        _S1895 = _S1910;
        _S1898 = 0.0f;
    }
    DiffPair_float_0 _S2038;
    (&_S2038)->primal_0 = _S1893;
    (&_S2038)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2038, _S1895);
    float _S2039 = (_S2038.differential_0 + _S1898) / _S1874;
    FixedArray<float, 27>  _S2040;
    _S2040[int(0)] = 0.0f;
    _S2040[int(1)] = 0.0f;
    _S2040[int(2)] = 0.0f;
    _S2040[int(3)] = 0.0f;
    _S2040[int(4)] = 0.0f;
    _S2040[int(5)] = 0.0f;
    _S2040[int(6)] = 0.0f;
    _S2040[int(7)] = 0.0f;
    _S2040[int(8)] = 0.0f;
    _S2040[int(9)] = 0.0f;
    _S2040[int(10)] = 0.0f;
    _S2040[int(11)] = 0.0f;
    _S2040[int(12)] = 0.0f;
    _S2040[int(13)] = 0.0f;
    _S2040[int(14)] = 0.0f;
    _S2040[int(15)] = 0.0f;
    _S2040[int(16)] = 0.0f;
    _S2040[int(17)] = 0.0f;
    _S2040[int(18)] = 0.0f;
    _S2040[int(19)] = 0.0f;
    _S2040[int(20)] = 0.0f;
    _S2040[int(21)] = 0.0f;
    _S2040[int(22)] = 0.0f;
    _S2040[int(23)] = 0.0f;
    _S2040[int(24)] = 0.0f;
    _S2040[int(25)] = 0.0f;
    _S2040[int(26)] = 0.0f;
    _S2040[int(19)] = _S2039;
    float _S2041 = _S2009 + _S2040[int(0)];
    float _S2042 = _S2010 + _S2040[int(1)];
    float _S2043 = _S2011 + _S2040[int(2)];
    float _S2044 = _S2012 + _S2040[int(3)];
    float _S2045 = _S2013 + _S2040[int(4)];
    float _S2046 = _S2014 + _S2040[int(5)];
    float _S2047 = _S2015 + _S2040[int(6)];
    float _S2048 = _S2016 + _S2040[int(7)];
    float _S2049 = _S2017 + _S2040[int(8)];
    float _S2050 = _S2018 + _S2040[int(9)];
    float _S2051 = _S2019 + _S2040[int(10)];
    float _S2052 = _S2020 + _S2040[int(11)];
    float _S2053 = _S2021 + _S2040[int(12)];
    float _S2054 = _S2022 + _S2040[int(13)];
    float _S2055 = _S2023 + _S2040[int(14)];
    float _S2056 = _S2024 + _S2040[int(15)];
    float _S2057 = _S2025 + _S2040[int(16)];
    float _S2058 = _S2026 + _S2040[int(17)];
    float _S2059 = _S2027 + _S2040[int(18)];
    float _S2060 = _S2028 + _S2040[int(19)];
    float _S2061 = _S2029 + _S2040[int(20)];
    float _S2062 = _S2030 + _S2040[int(21)];
    float _S2063 = _S2031 + _S2040[int(22)];
    float _S2064 = _S2032 + _S2040[int(23)];
    float _S2065 = _S2033 + _S2040[int(24)];
    float _S2066 = _S2034 + _S2040[int(25)];
    float _S2067 = _S2035 + _S2040[int(26)];
    if(_S1891)
    {
        float _S2068 = 200.0f * _S1910;
        float _S2069 = _S1892 * _S2068 + 0.5f * (_S1890 * _S2068);
        _S1892 = 0.0f;
        _S1895 = _S2069;
    }
    else
    {
        _S1892 = _S1910;
        _S1895 = 0.0f;
    }
    DiffPair_float_0 _S2070;
    (&_S2070)->primal_0 = _S1890;
    (&_S2070)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2070, _S1892);
    float _S2071 = (_S2070.differential_0 + _S1895) / _S1874;
    FixedArray<float, 27>  _S2072;
    _S2072[int(0)] = 0.0f;
    _S2072[int(1)] = 0.0f;
    _S2072[int(2)] = 0.0f;
    _S2072[int(3)] = 0.0f;
    _S2072[int(4)] = 0.0f;
    _S2072[int(5)] = 0.0f;
    _S2072[int(6)] = 0.0f;
    _S2072[int(7)] = 0.0f;
    _S2072[int(8)] = 0.0f;
    _S2072[int(9)] = 0.0f;
    _S2072[int(10)] = 0.0f;
    _S2072[int(11)] = 0.0f;
    _S2072[int(12)] = 0.0f;
    _S2072[int(13)] = 0.0f;
    _S2072[int(14)] = 0.0f;
    _S2072[int(15)] = 0.0f;
    _S2072[int(16)] = 0.0f;
    _S2072[int(17)] = 0.0f;
    _S2072[int(18)] = 0.0f;
    _S2072[int(19)] = 0.0f;
    _S2072[int(20)] = 0.0f;
    _S2072[int(21)] = 0.0f;
    _S2072[int(22)] = 0.0f;
    _S2072[int(23)] = 0.0f;
    _S2072[int(24)] = 0.0f;
    _S2072[int(25)] = 0.0f;
    _S2072[int(26)] = 0.0f;
    _S2072[int(18)] = _S2071;
    float _S2073 = _S2041 + _S2072[int(0)];
    float _S2074 = _S2042 + _S2072[int(1)];
    float _S2075 = _S2043 + _S2072[int(2)];
    float _S2076 = _S2044 + _S2072[int(3)];
    float _S2077 = _S2045 + _S2072[int(4)];
    float _S2078 = _S2046 + _S2072[int(5)];
    float _S2079 = _S2047 + _S2072[int(6)];
    float _S2080 = _S2048 + _S2072[int(7)];
    float _S2081 = _S2049 + _S2072[int(8)];
    float _S2082 = _S2050 + _S2072[int(9)];
    float _S2083 = _S2051 + _S2072[int(10)];
    float _S2084 = _S2052 + _S2072[int(11)];
    float _S2085 = _S2053 + _S2072[int(12)];
    float _S2086 = _S2054 + _S2072[int(13)];
    float _S2087 = _S2055 + _S2072[int(14)];
    float _S2088 = _S2056 + _S2072[int(15)];
    float _S2089 = _S2057 + _S2072[int(16)];
    float _S2090 = _S2058 + _S2072[int(17)];
    float _S2091 = _S2059 + _S2072[int(18)];
    float _S2092 = _S2060 + _S2072[int(19)];
    float _S2093 = _S2061 + _S2072[int(20)];
    float _S2094 = _S2062 + _S2072[int(21)];
    float _S2095 = _S2063 + _S2072[int(22)];
    float _S2096 = _S2064 + _S2072[int(23)];
    float _S2097 = _S2065 + _S2072[int(24)];
    float _S2098 = _S2066 + _S2072[int(25)];
    float _S2099 = _S2067 + _S2072[int(26)];
    if(_S1888)
    {
        float _S2100 = 200.0f * _S1910;
        float _S2101 = _S1889 * _S2100 + 0.5f * (_S1887 * _S2100);
        _S1889 = 0.0f;
        _S1892 = _S2101;
    }
    else
    {
        _S1889 = _S1910;
        _S1892 = 0.0f;
    }
    DiffPair_float_0 _S2102;
    (&_S2102)->primal_0 = _S1887;
    (&_S2102)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2102, _S1889);
    float _S2103 = (_S2102.differential_0 + _S1892) / _S1874;
    FixedArray<float, 27>  _S2104;
    _S2104[int(0)] = 0.0f;
    _S2104[int(1)] = 0.0f;
    _S2104[int(2)] = 0.0f;
    _S2104[int(3)] = 0.0f;
    _S2104[int(4)] = 0.0f;
    _S2104[int(5)] = 0.0f;
    _S2104[int(6)] = 0.0f;
    _S2104[int(7)] = 0.0f;
    _S2104[int(8)] = 0.0f;
    _S2104[int(9)] = 0.0f;
    _S2104[int(10)] = 0.0f;
    _S2104[int(11)] = 0.0f;
    _S2104[int(12)] = 0.0f;
    _S2104[int(13)] = 0.0f;
    _S2104[int(14)] = 0.0f;
    _S2104[int(15)] = 0.0f;
    _S2104[int(16)] = 0.0f;
    _S2104[int(17)] = 0.0f;
    _S2104[int(18)] = 0.0f;
    _S2104[int(19)] = 0.0f;
    _S2104[int(20)] = 0.0f;
    _S2104[int(21)] = 0.0f;
    _S2104[int(22)] = 0.0f;
    _S2104[int(23)] = 0.0f;
    _S2104[int(24)] = 0.0f;
    _S2104[int(25)] = 0.0f;
    _S2104[int(26)] = 0.0f;
    _S2104[int(17)] = _S2103;
    float _S2105 = _S2073 + _S2104[int(0)];
    float _S2106 = _S2074 + _S2104[int(1)];
    float _S2107 = _S2075 + _S2104[int(2)];
    float _S2108 = _S2076 + _S2104[int(3)];
    float _S2109 = _S2077 + _S2104[int(4)];
    float _S2110 = _S2078 + _S2104[int(5)];
    float _S2111 = _S2079 + _S2104[int(6)];
    float _S2112 = _S2080 + _S2104[int(7)];
    float _S2113 = _S2081 + _S2104[int(8)];
    float _S2114 = _S2082 + _S2104[int(9)];
    float _S2115 = _S2083 + _S2104[int(10)];
    float _S2116 = _S2084 + _S2104[int(11)];
    float _S2117 = _S2085 + _S2104[int(12)];
    float _S2118 = _S2086 + _S2104[int(13)];
    float _S2119 = _S2087 + _S2104[int(14)];
    float _S2120 = _S2088 + _S2104[int(15)];
    float _S2121 = _S2089 + _S2104[int(16)];
    float _S2122 = _S2090 + _S2104[int(17)];
    float _S2123 = _S2091 + _S2104[int(18)];
    float _S2124 = _S2092 + _S2104[int(19)];
    float _S2125 = _S2093 + _S2104[int(20)];
    float _S2126 = _S2094 + _S2104[int(21)];
    float _S2127 = _S2095 + _S2104[int(22)];
    float _S2128 = _S2096 + _S2104[int(23)];
    float _S2129 = _S2097 + _S2104[int(24)];
    float _S2130 = _S2098 + _S2104[int(25)];
    float _S2131 = _S2099 + _S2104[int(26)];
    if(_S1885)
    {
        float _S2132 = 200.0f * _S1910;
        float _S2133 = _S1886 * _S2132 + 0.5f * (_S1884 * _S2132);
        _S1886 = 0.0f;
        _S1889 = _S2133;
    }
    else
    {
        _S1886 = _S1910;
        _S1889 = 0.0f;
    }
    DiffPair_float_0 _S2134;
    (&_S2134)->primal_0 = _S1884;
    (&_S2134)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2134, _S1886);
    float _S2135 = (_S2134.differential_0 + _S1889) / _S1874;
    FixedArray<float, 27>  _S2136;
    _S2136[int(0)] = 0.0f;
    _S2136[int(1)] = 0.0f;
    _S2136[int(2)] = 0.0f;
    _S2136[int(3)] = 0.0f;
    _S2136[int(4)] = 0.0f;
    _S2136[int(5)] = 0.0f;
    _S2136[int(6)] = 0.0f;
    _S2136[int(7)] = 0.0f;
    _S2136[int(8)] = 0.0f;
    _S2136[int(9)] = 0.0f;
    _S2136[int(10)] = 0.0f;
    _S2136[int(11)] = 0.0f;
    _S2136[int(12)] = 0.0f;
    _S2136[int(13)] = 0.0f;
    _S2136[int(14)] = 0.0f;
    _S2136[int(15)] = 0.0f;
    _S2136[int(16)] = 0.0f;
    _S2136[int(17)] = 0.0f;
    _S2136[int(18)] = 0.0f;
    _S2136[int(19)] = 0.0f;
    _S2136[int(20)] = 0.0f;
    _S2136[int(21)] = 0.0f;
    _S2136[int(22)] = 0.0f;
    _S2136[int(23)] = 0.0f;
    _S2136[int(24)] = 0.0f;
    _S2136[int(25)] = 0.0f;
    _S2136[int(26)] = 0.0f;
    _S2136[int(16)] = _S2135;
    float _S2137 = _S2105 + _S2136[int(0)];
    float _S2138 = _S2106 + _S2136[int(1)];
    float _S2139 = _S2107 + _S2136[int(2)];
    float _S2140 = _S2108 + _S2136[int(3)];
    float _S2141 = _S2109 + _S2136[int(4)];
    float _S2142 = _S2110 + _S2136[int(5)];
    float _S2143 = _S2111 + _S2136[int(6)];
    float _S2144 = _S2112 + _S2136[int(7)];
    float _S2145 = _S2113 + _S2136[int(8)];
    float _S2146 = _S2114 + _S2136[int(9)];
    float _S2147 = _S2115 + _S2136[int(10)];
    float _S2148 = _S2116 + _S2136[int(11)];
    float _S2149 = _S2117 + _S2136[int(12)];
    float _S2150 = _S2118 + _S2136[int(13)];
    float _S2151 = _S2119 + _S2136[int(14)];
    float _S2152 = _S2120 + _S2136[int(15)];
    float _S2153 = _S2121 + _S2136[int(16)];
    float _S2154 = _S2122 + _S2136[int(17)];
    float _S2155 = _S2123 + _S2136[int(18)];
    float _S2156 = _S2124 + _S2136[int(19)];
    float _S2157 = _S2125 + _S2136[int(20)];
    float _S2158 = _S2126 + _S2136[int(21)];
    float _S2159 = _S2127 + _S2136[int(22)];
    float _S2160 = _S2128 + _S2136[int(23)];
    float _S2161 = _S2129 + _S2136[int(24)];
    float _S2162 = _S2130 + _S2136[int(25)];
    float _S2163 = _S2131 + _S2136[int(26)];
    if(_S1882)
    {
        float _S2164 = 200.0f * _S1910;
        float _S2165 = _S1883 * _S2164 + 0.5f * (_S1881 * _S2164);
        _S1883 = 0.0f;
        _S1886 = _S2165;
    }
    else
    {
        _S1883 = _S1910;
        _S1886 = 0.0f;
    }
    DiffPair_float_0 _S2166;
    (&_S2166)->primal_0 = _S1881;
    (&_S2166)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2166, _S1883);
    float _S2167 = (_S2166.differential_0 + _S1886) / _S1874;
    float _S2168 = _S1905 / _S1880;
    float _S2169 = _S1906 / _S1879;
    float _S2170 = 0.3333333432674408f * _S1907;
    float3  _S2171 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2172;
    (&_S2172)->primal_0 = mean_cy_1;
    (&_S2172)->differential_0 = _S2171;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2173;
    (&_S2173)->primal_0 = mean_cy_1;
    (&_S2173)->differential_0 = _S2171;
    s_bwd_prop_dot_0(&_S2172, &_S2173, _S2170);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2174;
    (&_S2174)->primal_0 = mean_cx_1;
    (&_S2174)->differential_0 = _S2171;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2175;
    (&_S2175)->primal_0 = mean_cx_1;
    (&_S2175)->differential_0 = _S2171;
    s_bwd_prop_dot_0(&_S2174, &_S2175, _S2170);
    float3  _S2176 = (_S2173.differential_0 + _S2172.differential_0) / _S1875;
    float3  _S2177 = (_S2175.differential_0 + _S2174.differential_0) / _S1875;
    FixedArray<float, 27>  _S2178;
    _S2178[int(0)] = 0.0f;
    _S2178[int(1)] = 0.0f;
    _S2178[int(2)] = 0.0f;
    _S2178[int(3)] = 0.0f;
    _S2178[int(4)] = 0.0f;
    _S2178[int(5)] = 0.0f;
    _S2178[int(6)] = 0.0f;
    _S2178[int(7)] = 0.0f;
    _S2178[int(8)] = 0.0f;
    _S2178[int(9)] = 0.0f;
    _S2178[int(10)] = 0.0f;
    _S2178[int(11)] = 0.0f;
    _S2178[int(12)] = 0.0f;
    _S2178[int(13)] = 0.0f;
    _S2178[int(14)] = 0.0f;
    _S2178[int(15)] = 0.0f;
    _S2178[int(16)] = 0.0f;
    _S2178[int(17)] = 0.0f;
    _S2178[int(18)] = 0.0f;
    _S2178[int(19)] = 0.0f;
    _S2178[int(20)] = 0.0f;
    _S2178[int(21)] = 0.0f;
    _S2178[int(22)] = 0.0f;
    _S2178[int(23)] = 0.0f;
    _S2178[int(24)] = 0.0f;
    _S2178[int(25)] = 0.0f;
    _S2178[int(26)] = 0.0f;
    _S2178[int(15)] = _S2167;
    _S2178[int(14)] = _S2168;
    _S2178[int(13)] = _S2168;
    _S2178[int(12)] = _S2168;
    _S2178[int(11)] = _S2168;
    _S2178[int(10)] = _S2168;
    _S2178[int(9)] = _S2169;
    _S2178[int(8)] = _S2169;
    _S2178[int(7)] = _S2169;
    _S2178[int(6)] = _S2176.z;
    _S2178[int(4)] = _S2176.y;
    _S2178[int(2)] = _S2176.x;
    _S2178[int(5)] = _S2177.z;
    _S2178[int(3)] = _S2177.y;
    _S2178[int(1)] = _S2177.x;
    float _S2179 = _S2137 + _S2178[int(0)];
    float _S2180 = _S2138 + _S2178[int(1)];
    float _S2181 = _S2139 + _S2178[int(2)];
    float _S2182 = _S2140 + _S2178[int(3)];
    float _S2183 = _S2141 + _S2178[int(4)];
    float _S2184 = _S2142 + _S2178[int(5)];
    float _S2185 = _S2143 + _S2178[int(6)];
    float _S2186 = _S2144 + _S2178[int(7)];
    float _S2187 = _S2145 + _S2178[int(8)];
    float _S2188 = _S2146 + _S2178[int(9)];
    float _S2189 = _S2147 + _S2178[int(10)];
    float _S2190 = _S2148 + _S2178[int(11)];
    float _S2191 = _S2149 + _S2178[int(12)];
    float _S2192 = _S2150 + _S2178[int(13)];
    float _S2193 = _S2151 + _S2178[int(14)];
    float _S2194 = _S2152 + _S2178[int(15)];
    float _S2195 = _S2153 + _S2178[int(16)];
    float _S2196 = _S2154 + _S2178[int(17)];
    float _S2197 = _S2155 + _S2178[int(18)];
    float _S2198 = _S2156 + _S2178[int(19)];
    float _S2199 = _S2157 + _S2178[int(20)];
    float _S2200 = _S2158 + _S2178[int(21)];
    float _S2201 = _S2159 + _S2178[int(22)];
    float _S2202 = _S2160 + _S2178[int(23)];
    float _S2203 = _S2161 + _S2178[int(24)];
    float _S2204 = _S2162 + _S2178[int(25)];
    float _S2205 = _S2163 + _S2178[int(26)];
    if(_S1877)
    {
        float _S2206 = 10.0f * _S1908;
        float _S2207 = _S1878 * _S2206 + 0.5f * (_S1876 * _S2206);
        _S1878 = 0.0f;
        _S1883 = _S2207;
    }
    else
    {
        _S1878 = _S1908;
        _S1883 = 0.0f;
    }
    DiffPair_float_0 _S2208;
    (&_S2208)->primal_0 = _S1876;
    (&_S2208)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2208, _S1878);
    float _S2209 = (_S2208.differential_0 + _S1883) / _S1874;
    FixedArray<float, 27>  _S2210;
    _S2210[int(0)] = 0.0f;
    _S2210[int(1)] = 0.0f;
    _S2210[int(2)] = 0.0f;
    _S2210[int(3)] = 0.0f;
    _S2210[int(4)] = 0.0f;
    _S2210[int(5)] = 0.0f;
    _S2210[int(6)] = 0.0f;
    _S2210[int(7)] = 0.0f;
    _S2210[int(8)] = 0.0f;
    _S2210[int(9)] = 0.0f;
    _S2210[int(10)] = 0.0f;
    _S2210[int(11)] = 0.0f;
    _S2210[int(12)] = 0.0f;
    _S2210[int(13)] = 0.0f;
    _S2210[int(14)] = 0.0f;
    _S2210[int(15)] = 0.0f;
    _S2210[int(16)] = 0.0f;
    _S2210[int(17)] = 0.0f;
    _S2210[int(18)] = 0.0f;
    _S2210[int(19)] = 0.0f;
    _S2210[int(20)] = 0.0f;
    _S2210[int(21)] = 0.0f;
    _S2210[int(22)] = 0.0f;
    _S2210[int(23)] = 0.0f;
    _S2210[int(24)] = 0.0f;
    _S2210[int(25)] = 0.0f;
    _S2210[int(26)] = 0.0f;
    _S2210[int(0)] = _S2209;
    FixedArray<float, 27>  _S2211 = {
        _S2179 + _S2210[int(0)], _S2180 + _S2210[int(1)], _S2181 + _S2210[int(2)], _S2182 + _S2210[int(3)], _S2183 + _S2210[int(4)], _S2184 + _S2210[int(5)], _S2185 + _S2210[int(6)], _S2186 + _S2210[int(7)], _S2187 + _S2210[int(8)], _S2188 + _S2210[int(9)], _S2189 + _S2210[int(10)], _S2190 + _S2210[int(11)], _S2191 + _S2210[int(12)], _S2192 + _S2210[int(13)], _S2193 + _S2210[int(14)], _S2194 + _S2210[int(15)], _S2195 + _S2210[int(16)], _S2196 + _S2210[int(17)], _S2197 + _S2210[int(18)], _S2198 + _S2210[int(19)], _S2199 + _S2210[int(20)], _S2200 + _S2210[int(21)], _S2201 + _S2210[int(22)], _S2202 + _S2210[int(23)], _S2203 + _S2210[int(24)], _S2204 + _S2210[int(25)], _S2205 + _S2210[int(26)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2211;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C27x3E_0 * _S2212, int _S2213, FixedArray<float, 6>  * _S2214, FixedArray<float, 6>  * _S2215)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2212, _S2213, _S2214, _S2215);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 27>  * raw_losses_3, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * grad_out_2, FixedArray<float, 27>  * _S2216)
{
    FixedArray<float, 27>  _S2217 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C27x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = *raw_losses_3;
    (&dp_raw_losses_1)->differential_0 = _S2217;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_2, loss_weights_2, grad_out_2);
    *_S2216 = (&dp_raw_losses_1)->differential_0;
    return;
}

