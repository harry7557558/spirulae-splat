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

inline __device__ DiffPair_float_0 _d_exp_1(DiffPair_float_0 * dpx_1)
{
    float _S7 = (F32_exp((dpx_1->primal_0)));
    DiffPair_float_0 _S8 = { _S7, _S7 * dpx_1->differential_0 };
    return _S8;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S9 = *dpx_2;
    float _S10;
    if(((*dpx_2).primal_0) > ((*dpy_0).primal_0))
    {
        _S10 = dOut_1;
    }
    else
    {
        if(((*dpx_2).primal_0) < ((*dpy_0).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_1;
        }
    }
    dpx_2->primal_0 = _S9.primal_0;
    dpx_2->differential_0 = _S10;
    DiffPair_float_0 _S11 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S9.primal_0))
    {
        _S10 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_2).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S11.primal_0;
    dpy_0->differential_0 = _S10;
    return;
}

inline __device__ DiffPair_float_0 _d_max_1(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_1)
{
    float _S12 = dpx_3->primal_0;
    float _S13 = dpy_1->primal_0;
    float _S14 = (F32_max((dpx_3->primal_0), (dpy_1->primal_0)));
    float _S15;
    if((dpx_3->primal_0) > (dpy_1->primal_0))
    {
        _S15 = dpx_3->differential_0;
    }
    else
    {
        if(_S12 < _S13)
        {
            _S15 = dpy_1->differential_0;
        }
        else
        {
            _S15 = 0.5f * (dpx_3->differential_0 + dpy_1->differential_0);
        }
    }
    DiffPair_float_0 _S16 = { _S14, _S15 };
    return _S16;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_4, float dOut_2)
{
    float _S17 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_4).primal_0)))))) * dOut_2;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S17;
    return;
}

inline __device__ DiffPair_float_0 _d_sqrt_1(DiffPair_float_0 * dpx_5)
{
    DiffPair_float_0 _S18 = { (F32_sqrt((dpx_5->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_5->primal_0)))))) * dpx_5->differential_0 };
    return _S18;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, float dOut_3)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_3;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_6).primal_0.x * dOut_3;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_3;
    *&((&y_d_result_0)->y) = (*dpx_6).primal_0.y * dOut_3;
    *&((&x_d_result_0)->z) = (*dpy_2).primal_0.z * dOut_3;
    *&((&y_d_result_0)->z) = (*dpx_6).primal_0.z * dOut_3;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
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

inline __device__ float length_1(float3  x_4)
{
    return (F32_sqrt((dot_0(x_4, x_4))));
}

inline __device__ float length_2(float2  x_5)
{
    return (F32_sqrt((dot_2(x_5, x_5))));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_7, float dOut_4)
{
    float _S19 = 1.0f / (*dpx_7).primal_0 * dOut_4;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S19;
    return;
}

inline __device__ DiffPair_float_0 _d_log_1(DiffPair_float_0 * dpx_8)
{
    DiffPair_float_0 _S20 = { (F32_log((dpx_8->primal_0))), 1.0f / dpx_8->primal_0 * dpx_8->differential_0 };
    return _S20;
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

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float3  dOut_5)
{
    float3  _S21 = exp_0((*dpx_9).primal_0) * dOut_5;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S21;
    return;
}

inline __device__ DiffPair_vectorx3Cfloatx2C3x3E_0 _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10)
{
    float3  _S22 = exp_0(dpx_10->primal_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S23 = { _S22, _S22 * dpx_10->differential_0 };
    return _S23;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpy_3, float dOut_6)
{
    DiffPair_float_0 _S24 = *dpx_11;
    float _S25;
    if(((*dpx_11).primal_0) < ((*dpy_3).primal_0))
    {
        _S25 = dOut_6;
    }
    else
    {
        if(((*dpx_11).primal_0) > ((*dpy_3).primal_0))
        {
            _S25 = 0.0f;
        }
        else
        {
            _S25 = 0.5f * dOut_6;
        }
    }
    dpx_11->primal_0 = _S24.primal_0;
    dpx_11->differential_0 = _S25;
    DiffPair_float_0 _S26 = *dpy_3;
    if(((*dpy_3).primal_0) < (_S24.primal_0))
    {
        _S25 = dOut_6;
    }
    else
    {
        if(((*dpy_3).primal_0) > ((*dpx_11).primal_0))
        {
            _S25 = 0.0f;
        }
        else
        {
            _S25 = 0.5f * dOut_6;
        }
    }
    dpy_3->primal_0 = _S26.primal_0;
    dpy_3->differential_0 = _S25;
    return;
}

inline __device__ DiffPair_float_0 _d_min_1(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpy_4)
{
    float _S27 = dpx_12->primal_0;
    float _S28 = dpy_4->primal_0;
    float _S29 = (F32_min((dpx_12->primal_0), (dpy_4->primal_0)));
    float _S30;
    if((dpx_12->primal_0) < (dpy_4->primal_0))
    {
        _S30 = dpx_12->differential_0;
    }
    else
    {
        if(_S27 > _S28)
        {
            _S30 = dpy_4->differential_0;
        }
        else
        {
            _S30 = 0.5f * (dpx_12->differential_0 + dpy_4->differential_0);
        }
    }
    DiffPair_float_0 _S31 = { _S29, _S30 };
    return _S31;
}

inline __device__ void per_splat_losses(float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S32)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    float3  _S33 = exp_0(scales_0);
    float _S34 = _S33.x;
    float _S35 = _S33.y;
    float _S36 = _S33.z;
    losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S34 + _S35 + _S36) / 3.0f;
    losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S34), (_S35)))), (_S36))) / (F32_min(((F32_min((_S34), (_S35)))), (_S36)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
    float3  _S37 = exp_0(make_float3 (2.0f) * scales_0);
    float x_7 = _S37.x;
    float y_3 = _S37.y;
    float z_0 = _S37.z;
    float s_0 = x_7 + y_3 + z_0;
    float s1_0 = (F32_max(((F32_max((x_7), (y_3)))), (z_0))) / s_0;
    float s3_0 = (F32_min(((F32_min((x_7), (y_3)))), (z_0))) / s_0;
    float s2_0 = 1.0f - s1_0 - s3_0;
    losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    *_S32 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_exp_0(float _S38)
{
    return (F32_exp((_S38)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S39)
{
    return exp_0(_S39);
}

inline __device__ float s_primal_ctx_log_0(float _S40)
{
    return (F32_log((_S40)));
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S41, float _S42)
{
    _d_log_0(_S41, _S42);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S43, float _S44)
{
    _d_exp_0(_S43, _S44);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S45, float3  _S46)
{
    _d_exp_vector_0(_S45, _S46);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S47, float _S48)
{
    _d_sqrt_0(_S47, _S48);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_13, float _s_dOut_0)
{
    float _S49 = (*dpx_13).primal_0.x;
    float _S50 = (*dpx_13).primal_0.y;
    float _S51 = (*dpx_13).primal_0.z;
    float _S52 = (*dpx_13).primal_0.w;
    DiffPair_float_0 _S53;
    (&_S53)->primal_0 = _S49 * _S49 + _S50 * _S50 + _S51 * _S51 + _S52 * _S52;
    (&_S53)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S53, _s_dOut_0);
    float _S54 = (*dpx_13).primal_0.w * _S53.differential_0;
    float _S55 = _S54 + _S54;
    float _S56 = (*dpx_13).primal_0.z * _S53.differential_0;
    float _S57 = _S56 + _S56;
    float _S58 = (*dpx_13).primal_0.y * _S53.differential_0;
    float _S59 = _S58 + _S58;
    float _S60 = (*dpx_13).primal_0.x * _S53.differential_0;
    float _S61 = _S60 + _S60;
    float4  _S62 = make_float4 (0.0f);
    *&((&_S62)->w) = _S55;
    *&((&_S62)->z) = _S57;
    *&((&_S62)->y) = _S59;
    *&((&_S62)->x) = _S61;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S62;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S63, float _S64)
{
    s_bwd_prop_length_impl_0(_S63, _S64);
    return;
}

inline __device__ void per_splat_losses_bwd(float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1)
{
    float _S65 = - opacity_1;
    float _S66 = 1.0f + s_primal_ctx_exp_0(_S65);
    float _S67 = _S66 * _S66;
    float _S68 = length_0(quat_1);
    float3  _S69 = s_primal_ctx_exp_1(scales_1);
    float _S70 = _S69.x;
    float _S71 = _S69.y;
    float _S72 = _S69.z;
    float _S73 = (F32_max((_S70), (_S71)));
    float _S74 = (F32_max((_S73), (_S72)));
    float _S75 = (F32_min((_S70), (_S71)));
    float _S76 = (F32_min((_S75), (_S72)));
    float _S77 = _S74 / _S76;
    float _S78 = _S76 * _S76;
    float3  _S79 = make_float3 (2.0f) * scales_1;
    float3  _S80 = s_primal_ctx_exp_1(_S79);
    float x_8 = _S80.x;
    float y_4 = _S80.y;
    float z_1 = _S80.z;
    float s_1 = x_8 + y_4 + z_1;
    float _S81 = (F32_max((x_8), (y_4)));
    float _S82 = (F32_max((_S81), (z_1)));
    float s1_1 = _S82 / s_1;
    float _S83 = s_1 * s_1;
    float _S84 = (F32_min((x_8), (y_4)));
    float _S85 = (F32_min((_S84), (z_1)));
    float s3_1 = _S85 / s_1;
    float s2_1 = 1.0f - s1_1 - s3_1;
    float _S86 = - s1_1;
    float _S87 = s_primal_ctx_log_0(s1_1);
    float _S88 = s_primal_ctx_log_0(s2_1);
    float _S89 = s_primal_ctx_log_0(s3_1);
    float _S90 = _S86 * _S87 - s2_1 * _S88 - s3_1 * _S89;
    float _S91 = s_primal_ctx_exp_0(_S90) - 0.99998998641967773f;
    float _S92 = erank_reg_weight_s3_1 * v_loss_0[int(3)];
    float _S93 = erank_reg_weight_1 * v_loss_0[int(3)];
    DiffPair_float_0 _S94;
    (&_S94)->primal_0 = - s_primal_ctx_log_0(_S91);
    (&_S94)->differential_0 = 0.0f;
    DiffPair_float_0 _S95;
    (&_S95)->primal_0 = 0.0f;
    (&_S95)->differential_0 = 0.0f;
    _d_max_0(&_S94, &_S95, _S93);
    float _S96 = - _S94.differential_0;
    DiffPair_float_0 _S97;
    (&_S97)->primal_0 = _S91;
    (&_S97)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S97, _S96);
    DiffPair_float_0 _S98;
    (&_S98)->primal_0 = _S90;
    (&_S98)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S98, _S97.differential_0);
    float _S99 = - _S98.differential_0;
    float _S100 = s3_1 * _S99;
    float _S101 = _S89 * _S99;
    DiffPair_float_0 _S102;
    (&_S102)->primal_0 = s3_1;
    (&_S102)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S102, _S100);
    float _S103 = s2_1 * _S99;
    float _S104 = _S88 * _S99;
    DiffPair_float_0 _S105;
    (&_S105)->primal_0 = s2_1;
    (&_S105)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S105, _S103);
    float _S106 = _S86 * _S98.differential_0;
    float _S107 = _S87 * _S98.differential_0;
    DiffPair_float_0 _S108;
    (&_S108)->primal_0 = s1_1;
    (&_S108)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S108, _S106);
    float _S109 = - _S107;
    float _S110 = - (_S104 + _S105.differential_0);
    float _S111 = (_S92 + _S101 + _S102.differential_0 + _S110) / _S83;
    float _S112 = _S85 * - _S111;
    float _S113 = s_1 * _S111;
    DiffPair_float_0 _S114;
    (&_S114)->primal_0 = _S84;
    (&_S114)->differential_0 = 0.0f;
    DiffPair_float_0 _S115;
    (&_S115)->primal_0 = z_1;
    (&_S115)->differential_0 = 0.0f;
    _d_min_0(&_S114, &_S115, _S113);
    DiffPair_float_0 _S116;
    (&_S116)->primal_0 = x_8;
    (&_S116)->differential_0 = 0.0f;
    DiffPair_float_0 _S117;
    (&_S117)->primal_0 = y_4;
    (&_S117)->differential_0 = 0.0f;
    _d_min_0(&_S116, &_S117, _S114.differential_0);
    float _S118 = (_S108.differential_0 + _S109 + _S110) / _S83;
    float _S119 = _S82 * - _S118;
    float _S120 = s_1 * _S118;
    DiffPair_float_0 _S121;
    (&_S121)->primal_0 = _S81;
    (&_S121)->differential_0 = 0.0f;
    DiffPair_float_0 _S122;
    (&_S122)->primal_0 = z_1;
    (&_S122)->differential_0 = 0.0f;
    _d_max_0(&_S121, &_S122, _S120);
    DiffPair_float_0 _S123;
    (&_S123)->primal_0 = x_8;
    (&_S123)->differential_0 = 0.0f;
    DiffPair_float_0 _S124;
    (&_S124)->primal_0 = y_4;
    (&_S124)->differential_0 = 0.0f;
    _d_max_0(&_S123, &_S124, _S121.differential_0);
    float _S125 = _S112 + _S119;
    float3  _S126 = make_float3 (_S116.differential_0 + _S123.differential_0 + _S125, _S117.differential_0 + _S124.differential_0 + _S125, _S115.differential_0 + _S122.differential_0 + _S125);
    float3  _S127 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S128;
    (&_S128)->primal_0 = _S79;
    (&_S128)->differential_0 = _S127;
    s_bwd_prop_exp_1(&_S128, _S126);
    float3  _S129 = make_float3 (2.0f) * _S128.differential_0;
    float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * v_loss_0[int(2)];
    DiffPair_float_0 _S130;
    (&_S130)->primal_0 = _S77;
    (&_S130)->differential_0 = 0.0f;
    DiffPair_float_0 _S131;
    (&_S131)->primal_0 = max_gauss_ratio_1;
    (&_S131)->differential_0 = 0.0f;
    _d_max_0(&_S130, &_S131, s_diff_scale_reg_T_0);
    float _S132 = _S130.differential_0 / _S78;
    float _S133 = _S74 * - _S132;
    float _S134 = _S76 * _S132;
    DiffPair_float_0 _S135;
    (&_S135)->primal_0 = _S75;
    (&_S135)->differential_0 = 0.0f;
    DiffPair_float_0 _S136;
    (&_S136)->primal_0 = _S72;
    (&_S136)->differential_0 = 0.0f;
    _d_min_0(&_S135, &_S136, _S133);
    DiffPair_float_0 _S137;
    (&_S137)->primal_0 = _S70;
    (&_S137)->differential_0 = 0.0f;
    DiffPair_float_0 _S138;
    (&_S138)->primal_0 = _S71;
    (&_S138)->differential_0 = 0.0f;
    _d_min_0(&_S137, &_S138, _S135.differential_0);
    DiffPair_float_0 _S139;
    (&_S139)->primal_0 = _S73;
    (&_S139)->differential_0 = 0.0f;
    DiffPair_float_0 _S140;
    (&_S140)->primal_0 = _S72;
    (&_S140)->differential_0 = 0.0f;
    _d_max_0(&_S139, &_S140, _S134);
    DiffPair_float_0 _S141;
    (&_S141)->primal_0 = _S70;
    (&_S141)->differential_0 = 0.0f;
    DiffPair_float_0 _S142;
    (&_S142)->primal_0 = _S71;
    (&_S142)->differential_0 = 0.0f;
    _d_max_0(&_S141, &_S142, _S139.differential_0);
    float _S143 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * v_loss_0[int(1)]);
    float3  _S144 = make_float3 (_S137.differential_0 + _S141.differential_0 + _S143, _S138.differential_0 + _S142.differential_0 + _S143, _S136.differential_0 + _S140.differential_0 + _S143);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S145;
    (&_S145)->primal_0 = scales_1;
    (&_S145)->differential_0 = _S127;
    s_bwd_prop_exp_1(&_S145, _S144);
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * v_loss_0[int(4)];
    float _S146 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S147;
    (&_S147)->primal_0 = _S68;
    (&_S147)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S147, _S146);
    float _S148 = _S147.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S149 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S150;
    (&_S150)->primal_0 = quat_1;
    (&_S150)->differential_0 = _S149;
    s_bwd_length_impl_0(&_S150, _S148);
    float _S151 = - (mcmc_opacity_reg_weight_1 * v_loss_0[int(0)] / _S67);
    DiffPair_float_0 _S152;
    (&_S152)->primal_0 = _S65;
    (&_S152)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S152, _S151);
    float _S153 = - _S152.differential_0;
    *v_scales_0 = _S129 + _S145.differential_0;
    *v_opacity_0 = _S153;
    *v_quat_0 = _S150.differential_0;
    return;
}

inline __device__ DiffPair_float_0 s_fwd_s_primal_ctx_exp_0(DiffPair_float_0 * _S154)
{
    DiffPair_float_0 _S155;
    (&_S155)->primal_0 = _S154->primal_0;
    (&_S155)->differential_0 = _S154->differential_0;
    DiffPair_float_0 _S156 = _d_exp_1(&_S155);
    DiffPair_float_0 _S157 = { _S156.primal_0, _S156.differential_0 };
    return _S157;
}

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_14)
{
    float _S158 = *&((&dpx_14->differential_0)->x) * *&((&dpx_14->primal_0)->x);
    float _S159 = *&((&dpx_14->differential_0)->y) * *&((&dpx_14->primal_0)->y);
    float _S160 = *&((&dpx_14->differential_0)->z) * *&((&dpx_14->primal_0)->z);
    float _S161 = *&((&dpx_14->differential_0)->w) * *&((&dpx_14->primal_0)->w);
    float s_diff_len_0 = _S158 + _S158 + (_S159 + _S159) + (_S160 + _S160) + (_S161 + _S161);
    DiffPair_float_0 _S162;
    (&_S162)->primal_0 = *&((&dpx_14->primal_0)->x) * *&((&dpx_14->primal_0)->x) + *&((&dpx_14->primal_0)->y) * *&((&dpx_14->primal_0)->y) + *&((&dpx_14->primal_0)->z) * *&((&dpx_14->primal_0)->z) + *&((&dpx_14->primal_0)->w) * *&((&dpx_14->primal_0)->w);
    (&_S162)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S163 = _d_sqrt_1(&_S162);
    DiffPair_float_0 _S164 = { _S163.primal_0, _S163.differential_0 };
    return _S164;
}

inline __device__ DiffPair_vectorx3Cfloatx2C3x3E_0 s_fwd_s_primal_ctx_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S165)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S166;
    (&_S166)->primal_0 = _S165->primal_0;
    (&_S166)->differential_0 = _S165->differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S167 = _d_exp_vector_1(&_S166);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S168 = { _S167.primal_0, _S167.differential_0 };
    return _S168;
}

inline __device__ DiffPair_float_0 s_fwd_s_primal_ctx_log_0(DiffPair_float_0 * _S169)
{
    DiffPair_float_0 _S170;
    (&_S170)->primal_0 = _S169->primal_0;
    (&_S170)->differential_0 = _S169->differential_0;
    DiffPair_float_0 _S171 = _d_log_1(&_S170);
    DiffPair_float_0 _S172 = { _S171.primal_0, _S171.differential_0 };
    return _S172;
}

struct DiffPair_0
{
    DiffPair_float_0 primal_0;
    DiffPair_float_0 differential_0;
};

inline __device__ void s_fwd_d_max_0(DiffPair_0 * dpdpx_0, DiffPair_0 * dpdpy_0, DiffPair_float_0 * dpdOut_0)
{
    DiffPair_0 _S173 = *dpdpx_0;
    DiffPair_0 _S174 = *dpdpy_0;
    float _S175 = dpdOut_0->differential_0;
    float _S176 = dpdOut_0->primal_0;
    float _S177;
    float _S178;
    if(((*dpdpx_0).primal_0.primal_0) > ((*dpdpy_0).primal_0.primal_0))
    {
        _S177 = _S176;
        _S178 = _S175;
    }
    else
    {
        if((_S173.primal_0.primal_0) < (_S174.primal_0.primal_0))
        {
            _S177 = 0.0f;
            _S178 = 0.0f;
        }
        else
        {
            float _S179 = _S175 * 0.5f;
            _S177 = 0.5f * _S176;
            _S178 = _S179;
        }
    }
    DiffPair_float_0 _S180 = { _S173.primal_0.primal_0, _S177 };
    DiffPair_float_0 _S181 = { _S173.differential_0.primal_0, _S178 };
    if((_S174.primal_0.primal_0) > (_S173.primal_0.primal_0))
    {
        _S177 = _S176;
        _S178 = _S175;
    }
    else
    {
        if((_S174.primal_0.primal_0) < (_S173.primal_0.primal_0))
        {
            _S177 = 0.0f;
            _S178 = 0.0f;
        }
        else
        {
            float _S182 = _S175 * 0.5f;
            _S177 = 0.5f * _S176;
            _S178 = _S182;
        }
    }
    DiffPair_float_0 _S183 = { _S174.primal_0.primal_0, _S177 };
    DiffPair_float_0 _S184 = { _S174.differential_0.primal_0, _S178 };
    dpdpx_0->primal_0 = _S180;
    dpdpx_0->differential_0 = _S181;
    dpdpy_0->primal_0 = _S183;
    dpdpy_0->differential_0 = _S184;
    return;
}

inline __device__ void s_fwd_d_log_0(DiffPair_0 * dpdpx_1, DiffPair_float_0 * dpdOut_1)
{
    float _S185 = 1.0f / (*dpdpx_1).primal_0.primal_0;
    DiffPair_float_0 _S186 = { (*dpdpx_1).primal_0.primal_0, _S185 * dpdOut_1->primal_0 };
    DiffPair_float_0 _S187 = { (*dpdpx_1).differential_0.primal_0, (0.0f - (*dpdpx_1).differential_0.primal_0) / ((*dpdpx_1).primal_0.primal_0 * (*dpdpx_1).primal_0.primal_0) * dpdOut_1->primal_0 + dpdOut_1->differential_0 * _S185 };
    dpdpx_1->primal_0 = _S186;
    dpdpx_1->differential_0 = _S187;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_log_0(DiffPair_0 * _S188, DiffPair_float_0 * _S189)
{
    DiffPair_0 _S190;
    (&_S190)->primal_0 = (*_S188).primal_0;
    (&_S190)->differential_0 = (*_S188).differential_0;
    DiffPair_float_0 _S191;
    (&_S191)->primal_0 = _S189->primal_0;
    (&_S191)->differential_0 = _S189->differential_0;
    s_fwd_d_log_0(&_S190, &_S191);
    _S188->primal_0 = _S190.primal_0;
    _S188->differential_0 = _S190.differential_0;
    return;
}

inline __device__ void s_fwd_d_exp_0(DiffPair_0 * dpdpx_2, DiffPair_float_0 * dpdOut_2)
{
    DiffPair_float_0 _S192;
    (&_S192)->primal_0 = (*dpdpx_2).primal_0.primal_0;
    (&_S192)->differential_0 = (*dpdpx_2).differential_0.primal_0;
    DiffPair_float_0 _S193 = _d_exp_1(&_S192);
    DiffPair_float_0 _S194 = { (*dpdpx_2).primal_0.primal_0, _S193.primal_0 * dpdOut_2->primal_0 };
    DiffPair_float_0 _S195 = { (*dpdpx_2).differential_0.primal_0, _S193.differential_0 * dpdOut_2->primal_0 + dpdOut_2->differential_0 * _S193.primal_0 };
    dpdpx_2->primal_0 = _S194;
    dpdpx_2->differential_0 = _S195;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_exp_0(DiffPair_0 * _S196, DiffPair_float_0 * _S197)
{
    DiffPair_0 _S198;
    (&_S198)->primal_0 = (*_S196).primal_0;
    (&_S198)->differential_0 = (*_S196).differential_0;
    DiffPair_float_0 _S199;
    (&_S199)->primal_0 = _S197->primal_0;
    (&_S199)->differential_0 = _S197->differential_0;
    s_fwd_d_exp_0(&_S198, &_S199);
    _S196->primal_0 = _S198.primal_0;
    _S196->differential_0 = _S198.differential_0;
    return;
}

inline __device__ void s_fwd_d_min_0(DiffPair_0 * dpdpx_3, DiffPair_0 * dpdpy_1, DiffPair_float_0 * dpdOut_3)
{
    DiffPair_0 _S200 = *dpdpx_3;
    DiffPair_0 _S201 = *dpdpy_1;
    float _S202 = dpdOut_3->differential_0;
    float _S203 = dpdOut_3->primal_0;
    float _S204;
    float _S205;
    if(((*dpdpx_3).primal_0.primal_0) < ((*dpdpy_1).primal_0.primal_0))
    {
        _S204 = _S203;
        _S205 = _S202;
    }
    else
    {
        if((_S200.primal_0.primal_0) > (_S201.primal_0.primal_0))
        {
            _S204 = 0.0f;
            _S205 = 0.0f;
        }
        else
        {
            float _S206 = _S202 * 0.5f;
            _S204 = 0.5f * _S203;
            _S205 = _S206;
        }
    }
    DiffPair_float_0 _S207 = { _S200.primal_0.primal_0, _S204 };
    DiffPair_float_0 _S208 = { _S200.differential_0.primal_0, _S205 };
    if((_S201.primal_0.primal_0) < (_S200.primal_0.primal_0))
    {
        _S204 = _S203;
        _S205 = _S202;
    }
    else
    {
        if((_S201.primal_0.primal_0) > (_S200.primal_0.primal_0))
        {
            _S204 = 0.0f;
            _S205 = 0.0f;
        }
        else
        {
            float _S209 = _S202 * 0.5f;
            _S204 = 0.5f * _S203;
            _S205 = _S209;
        }
    }
    DiffPair_float_0 _S210 = { _S201.primal_0.primal_0, _S204 };
    DiffPair_float_0 _S211 = { _S201.differential_0.primal_0, _S205 };
    dpdpx_3->primal_0 = _S207;
    dpdpx_3->differential_0 = _S208;
    dpdpy_1->primal_0 = _S210;
    dpdpy_1->differential_0 = _S211;
    return;
}

struct DiffPair_1
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 differential_0;
};

inline __device__ void s_fwd_d_exp_vector_0(DiffPair_1 * dpdpx_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S212;
    (&_S212)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S212)->differential_0 = (*dpdpx_4).differential_0.primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S213 = _d_exp_vector_1(&_S212);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S214 = { (*dpdpx_4).primal_0.primal_0, _S213.primal_0 * dpdOut_4->primal_0 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S215 = { (*dpdpx_4).differential_0.primal_0, _S213.differential_0 * dpdOut_4->primal_0 + dpdOut_4->differential_0 * _S213.primal_0 };
    dpdpx_4->primal_0 = _S214;
    dpdpx_4->differential_0 = _S215;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_exp_1(DiffPair_1 * _S216, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S217)
{
    DiffPair_1 _S218;
    (&_S218)->primal_0 = (*_S216).primal_0;
    (&_S218)->differential_0 = (*_S216).differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S219;
    (&_S219)->primal_0 = _S217->primal_0;
    (&_S219)->differential_0 = _S217->differential_0;
    s_fwd_d_exp_vector_0(&_S218, &_S219);
    _S216->primal_0 = _S218.primal_0;
    _S216->differential_0 = _S218.differential_0;
    return;
}

struct DiffPair_2
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 differential_0;
};

inline __device__ void s_fwd_d_sqrt_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dpdOut_5)
{
    DiffPair_float_0 _S220;
    (&_S220)->primal_0 = 1.00000001168609742e-07f;
    (&_S220)->differential_0 = 0.0f;
    DiffPair_float_0 _S221;
    (&_S221)->primal_0 = (*dpdpx_5).primal_0.primal_0;
    (&_S221)->differential_0 = (*dpdpx_5).differential_0.primal_0;
    DiffPair_float_0 _S222 = _d_max_1(&_S220, &_S221);
    DiffPair_float_0 _S223;
    (&_S223)->primal_0 = _S222.primal_0;
    (&_S223)->differential_0 = _S222.differential_0;
    DiffPair_float_0 _S224 = _d_sqrt_1(&_S223);
    float _S225 = 0.5f / _S224.primal_0;
    DiffPair_float_0 _S226 = { (*dpdpx_5).primal_0.primal_0, _S225 * dpdOut_5->primal_0 };
    DiffPair_float_0 _S227 = { (*dpdpx_5).differential_0.primal_0, (0.0f - 0.5f * _S224.differential_0) / (_S224.primal_0 * _S224.primal_0) * dpdOut_5->primal_0 + dpdOut_5->differential_0 * _S225 };
    dpdpx_5->primal_0 = _S226;
    dpdpx_5->differential_0 = _S227;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_sqrt_0(DiffPair_0 * _S228, DiffPair_float_0 * _S229)
{
    DiffPair_0 _S230;
    (&_S230)->primal_0 = (*_S228).primal_0;
    (&_S230)->differential_0 = (*_S228).differential_0;
    DiffPair_float_0 _S231;
    (&_S231)->primal_0 = _S229->primal_0;
    (&_S231)->differential_0 = _S229->differential_0;
    s_fwd_d_sqrt_0(&_S230, &_S231);
    _S228->primal_0 = _S230.primal_0;
    _S228->differential_0 = _S230.differential_0;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_length_impl_0(DiffPair_2 * dpdpx_6, DiffPair_float_0 * dp_s_dOut_0)
{
    float _S232 = (*dpdpx_6).primal_0.primal_0.x;
    float _S233 = (*dpdpx_6).differential_0.primal_0.x * (*dpdpx_6).primal_0.primal_0.x;
    float _S234 = (*dpdpx_6).primal_0.primal_0.y;
    float _S235 = (*dpdpx_6).differential_0.primal_0.y * (*dpdpx_6).primal_0.primal_0.y;
    float _S236 = (*dpdpx_6).primal_0.primal_0.z;
    float _S237 = (*dpdpx_6).differential_0.primal_0.z * (*dpdpx_6).primal_0.primal_0.z;
    float _S238 = (*dpdpx_6).primal_0.primal_0.w;
    float _S239 = (*dpdpx_6).differential_0.primal_0.w * (*dpdpx_6).primal_0.primal_0.w;
    DiffPair_float_0 _S240 = { _S232 * _S232 + _S234 * _S234 + _S236 * _S236 + _S238 * _S238, 0.0f };
    DiffPair_float_0 _S241 = { _S233 + _S233 + (_S235 + _S235) + (_S237 + _S237) + (_S239 + _S239), 0.0f };
    DiffPair_0 _S242;
    (&_S242)->primal_0 = _S240;
    (&_S242)->differential_0 = _S241;
    DiffPair_float_0 _S243;
    (&_S243)->primal_0 = dp_s_dOut_0->primal_0;
    (&_S243)->differential_0 = dp_s_dOut_0->differential_0;
    s_fwd_s_bwd_prop_sqrt_0(&_S242, &_S243);
    float _S244 = (*dpdpx_6).primal_0.primal_0.w * _S242.primal_0.differential_0;
    float _S245 = (*dpdpx_6).differential_0.primal_0.w * _S242.primal_0.differential_0 + _S242.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.w;
    float _S246 = _S244 + _S244;
    float _S247 = _S245 + _S245;
    float _S248 = (*dpdpx_6).primal_0.primal_0.z * _S242.primal_0.differential_0;
    float _S249 = (*dpdpx_6).differential_0.primal_0.z * _S242.primal_0.differential_0 + _S242.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.z;
    float _S250 = _S248 + _S248;
    float _S251 = _S249 + _S249;
    float _S252 = (*dpdpx_6).primal_0.primal_0.y * _S242.primal_0.differential_0;
    float _S253 = (*dpdpx_6).differential_0.primal_0.y * _S242.primal_0.differential_0 + _S242.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.y;
    float _S254 = _S252 + _S252;
    float _S255 = _S253 + _S253;
    float _S256 = (*dpdpx_6).primal_0.primal_0.x * _S242.primal_0.differential_0;
    float _S257 = (*dpdpx_6).differential_0.primal_0.x * _S242.primal_0.differential_0 + _S242.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.x;
    float _S258 = _S256 + _S256;
    float _S259 = _S257 + _S257;
    float4  _S260 = make_float4 (0.0f);
    float4  _S261 = _S260;
    *&((&_S261)->w) = _S246;
    float4  _S262 = _S260;
    *&((&_S262)->w) = _S247;
    *&((&_S261)->z) = _S250;
    *&((&_S262)->z) = _S251;
    *&((&_S261)->y) = _S254;
    *&((&_S262)->y) = _S255;
    *&((&_S261)->x) = _S258;
    *&((&_S262)->x) = _S259;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S263 = { (*dpdpx_6).primal_0.primal_0, _S261 };
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S264 = { (*dpdpx_6).differential_0.primal_0, _S262 };
    dpdpx_6->primal_0 = _S263;
    dpdpx_6->differential_0 = _S264;
    return;
}

inline __device__ void s_fwd_s_bwd_length_impl_0(DiffPair_2 * _S265, DiffPair_float_0 * _S266)
{
    DiffPair_2 _S267;
    (&_S267)->primal_0 = (*_S265).primal_0;
    (&_S267)->differential_0 = (*_S265).differential_0;
    DiffPair_float_0 _S268;
    (&_S268)->primal_0 = _S266->primal_0;
    (&_S268)->differential_0 = _S266->differential_0;
    s_fwd_s_bwd_prop_length_impl_0(&_S267, &_S268);
    _S265->primal_0 = _S267.primal_0;
    _S265->differential_0 = _S267.differential_0;
    return;
}

inline __device__ void s_fwd_per_splat_losses_bwd_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, FixedArray<float, 5>  * v_loss_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpv_scales_0, DiffPair_float_0 * dpv_opacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpv_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float _S269 = - dpopacity_0->primal_0;
    float _S270 = - dpopacity_0->differential_0;
    DiffPair_float_0 _S271;
    (&_S271)->primal_0 = _S269;
    (&_S271)->differential_0 = _S270;
    DiffPair_float_0 _S272 = s_fwd_s_primal_ctx_exp_0(&_S271);
    float _S273 = 1.0f + _S272.primal_0;
    float _S274 = _S273 * _S273;
    float _S275 = _S272.differential_0 * _S273;
    float _S276 = _S275 + _S275;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S277;
    (&_S277)->primal_0 = dpquat_0->primal_0;
    (&_S277)->differential_0 = dpquat_0->differential_0;
    DiffPair_float_0 _S278 = s_fwd_length_impl_0(&_S277);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S279;
    (&_S279)->primal_0 = dpscales_0->primal_0;
    (&_S279)->differential_0 = dpscales_0->differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S280 = s_fwd_s_primal_ctx_exp_1(&_S279);
    float _S281 = _S280.primal_0.x;
    float _S282 = _S280.differential_0.x;
    float _S283 = _S280.primal_0.y;
    float _S284 = _S280.differential_0.y;
    float _S285 = _S280.primal_0.z;
    float _S286 = _S280.differential_0.z;
    DiffPair_float_0 _S287;
    (&_S287)->primal_0 = _S281;
    (&_S287)->differential_0 = _S282;
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = _S283;
    (&_S288)->differential_0 = _S284;
    DiffPair_float_0 _S289 = _d_max_1(&_S287, &_S288);
    DiffPair_float_0 _S290;
    (&_S290)->primal_0 = _S289.primal_0;
    (&_S290)->differential_0 = _S289.differential_0;
    DiffPair_float_0 _S291;
    (&_S291)->primal_0 = _S285;
    (&_S291)->differential_0 = _S286;
    DiffPair_float_0 _S292 = _d_max_1(&_S290, &_S291);
    DiffPair_float_0 _S293;
    (&_S293)->primal_0 = _S281;
    (&_S293)->differential_0 = _S282;
    DiffPair_float_0 _S294;
    (&_S294)->primal_0 = _S283;
    (&_S294)->differential_0 = _S284;
    DiffPair_float_0 _S295 = _d_min_1(&_S293, &_S294);
    DiffPair_float_0 _S296;
    (&_S296)->primal_0 = _S295.primal_0;
    (&_S296)->differential_0 = _S295.differential_0;
    DiffPair_float_0 _S297;
    (&_S297)->primal_0 = _S285;
    (&_S297)->differential_0 = _S286;
    DiffPair_float_0 _S298 = _d_min_1(&_S296, &_S297);
    float _S299 = _S292.primal_0 / _S298.primal_0;
    float _S300 = _S298.primal_0 * _S298.primal_0;
    float _S301 = (_S292.differential_0 * _S298.primal_0 - _S292.primal_0 * _S298.differential_0) / _S300;
    float _S302 = _S298.differential_0 * _S298.primal_0;
    float _S303 = _S302 + _S302;
    float3  _S304 = make_float3 (2.0f) * dpscales_0->primal_0;
    float3  _S305 = dpscales_0->differential_0 * make_float3 (2.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S306;
    (&_S306)->primal_0 = _S304;
    (&_S306)->differential_0 = _S305;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S307 = s_fwd_s_primal_ctx_exp_1(&_S306);
    float x_9 = _S307.primal_0.x;
    float s_diff_x_0 = _S307.differential_0.x;
    float y_5 = _S307.primal_0.y;
    float s_diff_y_0 = _S307.differential_0.y;
    float z_2 = _S307.primal_0.z;
    float s_diff_z_0 = _S307.differential_0.z;
    float s_2 = x_9 + y_5 + z_2;
    float s_diff_s_0 = s_diff_x_0 + s_diff_y_0 + s_diff_z_0;
    DiffPair_float_0 _S308;
    (&_S308)->primal_0 = x_9;
    (&_S308)->differential_0 = s_diff_x_0;
    DiffPair_float_0 _S309;
    (&_S309)->primal_0 = y_5;
    (&_S309)->differential_0 = s_diff_y_0;
    DiffPair_float_0 _S310 = _d_max_1(&_S308, &_S309);
    DiffPair_float_0 _S311;
    (&_S311)->primal_0 = _S310.primal_0;
    (&_S311)->differential_0 = _S310.differential_0;
    DiffPair_float_0 _S312;
    (&_S312)->primal_0 = z_2;
    (&_S312)->differential_0 = s_diff_z_0;
    DiffPair_float_0 _S313 = _d_max_1(&_S311, &_S312);
    float s1_2 = _S313.primal_0 / s_2;
    float _S314 = s_2 * s_2;
    float s_diff_s1_0 = (_S313.differential_0 * s_2 - _S313.primal_0 * s_diff_s_0) / _S314;
    float _S315 = s_diff_s_0 * s_2;
    float _S316 = _S315 + _S315;
    DiffPair_float_0 _S317;
    (&_S317)->primal_0 = x_9;
    (&_S317)->differential_0 = s_diff_x_0;
    DiffPair_float_0 _S318;
    (&_S318)->primal_0 = y_5;
    (&_S318)->differential_0 = s_diff_y_0;
    DiffPair_float_0 _S319 = _d_min_1(&_S317, &_S318);
    DiffPair_float_0 _S320;
    (&_S320)->primal_0 = _S319.primal_0;
    (&_S320)->differential_0 = _S319.differential_0;
    DiffPair_float_0 _S321;
    (&_S321)->primal_0 = z_2;
    (&_S321)->differential_0 = s_diff_z_0;
    DiffPair_float_0 _S322 = _d_min_1(&_S320, &_S321);
    float s3_2 = _S322.primal_0 / s_2;
    float s_diff_s3_0 = (_S322.differential_0 * s_2 - _S322.primal_0 * s_diff_s_0) / _S314;
    float s2_2 = 1.0f - s1_2 - s3_2;
    float s_diff_s2_0 = 0.0f - s_diff_s1_0 - s_diff_s3_0;
    float _S323 = - s1_2;
    float _S324 = - s_diff_s1_0;
    DiffPair_float_0 _S325;
    (&_S325)->primal_0 = s1_2;
    (&_S325)->differential_0 = s_diff_s1_0;
    DiffPair_float_0 _S326 = s_fwd_s_primal_ctx_log_0(&_S325);
    float _S327 = _S323 * _S326.primal_0;
    float _S328 = _S324 * _S326.primal_0 + _S326.differential_0 * _S323;
    DiffPair_float_0 _S329;
    (&_S329)->primal_0 = s2_2;
    (&_S329)->differential_0 = s_diff_s2_0;
    DiffPair_float_0 _S330 = s_fwd_s_primal_ctx_log_0(&_S329);
    float _S331 = _S327 - s2_2 * _S330.primal_0;
    float _S332 = _S328 - (s_diff_s2_0 * _S330.primal_0 + _S330.differential_0 * s2_2);
    DiffPair_float_0 _S333;
    (&_S333)->primal_0 = s3_2;
    (&_S333)->differential_0 = s_diff_s3_0;
    DiffPair_float_0 _S334 = s_fwd_s_primal_ctx_log_0(&_S333);
    float _S335 = _S331 - s3_2 * _S334.primal_0;
    float _S336 = _S332 - (s_diff_s3_0 * _S334.primal_0 + _S334.differential_0 * s3_2);
    DiffPair_float_0 _S337;
    (&_S337)->primal_0 = _S335;
    (&_S337)->differential_0 = _S336;
    DiffPair_float_0 _S338 = s_fwd_s_primal_ctx_exp_0(&_S337);
    float _S339 = _S338.primal_0 - 0.99998998641967773f;
    DiffPair_float_0 _S340;
    (&_S340)->primal_0 = _S339;
    (&_S340)->differential_0 = _S338.differential_0;
    DiffPair_float_0 _S341 = s_fwd_s_primal_ctx_log_0(&_S340);
    float _S342 = erank_reg_weight_s3_2 * (*v_loss_1)[int(3)];
    float _S343 = erank_reg_weight_2 * (*v_loss_1)[int(3)];
    DiffPair_float_0 _S344 = { - _S341.primal_0, 0.0f };
    DiffPair_float_0 _S345 = { - _S341.differential_0, 0.0f };
    DiffPair_float_0 _S346 = { 0.0f, 0.0f };
    DiffPair_0 _S347;
    (&_S347)->primal_0 = _S344;
    (&_S347)->differential_0 = _S345;
    DiffPair_0 _S348;
    (&_S348)->primal_0 = _S346;
    (&_S348)->differential_0 = _S346;
    DiffPair_float_0 _S349;
    (&_S349)->primal_0 = _S343;
    (&_S349)->differential_0 = 0.0f;
    s_fwd_d_max_0(&_S347, &_S348, &_S349);
    float _S350 = - _S347.primal_0.differential_0;
    float _S351 = - _S347.differential_0.differential_0;
    DiffPair_float_0 _S352 = { _S339, 0.0f };
    DiffPair_float_0 _S353 = { _S338.differential_0, 0.0f };
    DiffPair_0 _S354;
    (&_S354)->primal_0 = _S352;
    (&_S354)->differential_0 = _S353;
    DiffPair_float_0 _S355;
    (&_S355)->primal_0 = _S350;
    (&_S355)->differential_0 = _S351;
    s_fwd_s_bwd_prop_log_0(&_S354, &_S355);
    DiffPair_float_0 _S356 = { _S335, 0.0f };
    DiffPair_float_0 _S357 = { _S336, 0.0f };
    DiffPair_0 _S358;
    (&_S358)->primal_0 = _S356;
    (&_S358)->differential_0 = _S357;
    DiffPair_float_0 _S359;
    (&_S359)->primal_0 = _S354.primal_0.differential_0;
    (&_S359)->differential_0 = _S354.differential_0.differential_0;
    s_fwd_s_bwd_prop_exp_0(&_S358, &_S359);
    float _S360 = - _S358.primal_0.differential_0;
    float _S361 = - _S358.differential_0.differential_0;
    float _S362 = s3_2 * _S360;
    float _S363 = s_diff_s3_0 * _S360 + _S361 * s3_2;
    float _S364 = _S334.primal_0 * _S360;
    float _S365 = _S334.differential_0 * _S360 + _S361 * _S334.primal_0;
    DiffPair_float_0 _S366 = { s3_2, 0.0f };
    DiffPair_float_0 _S367 = { s_diff_s3_0, 0.0f };
    DiffPair_0 _S368;
    (&_S368)->primal_0 = _S366;
    (&_S368)->differential_0 = _S367;
    DiffPair_float_0 _S369;
    (&_S369)->primal_0 = _S362;
    (&_S369)->differential_0 = _S363;
    s_fwd_s_bwd_prop_log_0(&_S368, &_S369);
    float _S370 = s2_2 * _S360;
    float _S371 = s_diff_s2_0 * _S360 + _S361 * s2_2;
    float _S372 = _S330.primal_0 * _S360;
    float _S373 = _S330.differential_0 * _S360 + _S361 * _S330.primal_0;
    DiffPair_float_0 _S374 = { s2_2, 0.0f };
    DiffPair_float_0 _S375 = { s_diff_s2_0, 0.0f };
    DiffPair_0 _S376;
    (&_S376)->primal_0 = _S374;
    (&_S376)->differential_0 = _S375;
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = _S370;
    (&_S377)->differential_0 = _S371;
    s_fwd_s_bwd_prop_log_0(&_S376, &_S377);
    float _S378 = _S323 * _S358.primal_0.differential_0;
    float _S379 = _S324 * _S358.primal_0.differential_0 + _S358.differential_0.differential_0 * _S323;
    float _S380 = _S326.primal_0 * _S358.primal_0.differential_0;
    float _S381 = _S326.differential_0 * _S358.primal_0.differential_0 + _S358.differential_0.differential_0 * _S326.primal_0;
    DiffPair_float_0 _S382 = { s1_2, 0.0f };
    DiffPair_float_0 _S383 = { s_diff_s1_0, 0.0f };
    DiffPair_0 _S384;
    (&_S384)->primal_0 = _S382;
    (&_S384)->differential_0 = _S383;
    DiffPair_float_0 _S385;
    (&_S385)->primal_0 = _S378;
    (&_S385)->differential_0 = _S379;
    s_fwd_s_bwd_prop_log_0(&_S384, &_S385);
    float _S386 = - _S380;
    float _S387 = - _S381;
    float _S388 = - (_S372 + _S376.primal_0.differential_0);
    float _S389 = - (_S373 + _S376.differential_0.differential_0);
    float _S390 = _S342 + _S364 + _S368.primal_0.differential_0 + _S388;
    float _S391 = _S390 / _S314;
    float _S392 = _S314 * _S314;
    float _S393 = ((_S365 + _S368.differential_0.differential_0 + _S389) * _S314 - _S390 * _S316) / _S392;
    float _S394 = - _S391;
    float _S395 = _S322.primal_0 * _S394;
    float _S396 = _S322.differential_0 * _S394 + - _S393 * _S322.primal_0;
    float _S397 = s_2 * _S391;
    float _S398 = s_diff_s_0 * _S391 + _S393 * s_2;
    DiffPair_float_0 _S399 = { _S319.primal_0, 0.0f };
    DiffPair_float_0 _S400 = { _S319.differential_0, 0.0f };
    DiffPair_float_0 _S401 = { z_2, 0.0f };
    DiffPair_float_0 _S402 = { s_diff_z_0, 0.0f };
    DiffPair_0 _S403;
    (&_S403)->primal_0 = _S399;
    (&_S403)->differential_0 = _S400;
    DiffPair_0 _S404;
    (&_S404)->primal_0 = _S401;
    (&_S404)->differential_0 = _S402;
    DiffPair_float_0 _S405;
    (&_S405)->primal_0 = _S397;
    (&_S405)->differential_0 = _S398;
    s_fwd_d_min_0(&_S403, &_S404, &_S405);
    DiffPair_float_0 _S406 = { x_9, 0.0f };
    DiffPair_float_0 _S407 = { s_diff_x_0, 0.0f };
    DiffPair_float_0 _S408 = { y_5, 0.0f };
    DiffPair_float_0 _S409 = { s_diff_y_0, 0.0f };
    DiffPair_0 _S410;
    (&_S410)->primal_0 = _S406;
    (&_S410)->differential_0 = _S407;
    DiffPair_0 _S411;
    (&_S411)->primal_0 = _S408;
    (&_S411)->differential_0 = _S409;
    DiffPair_float_0 _S412;
    (&_S412)->primal_0 = _S403.primal_0.differential_0;
    (&_S412)->differential_0 = _S403.differential_0.differential_0;
    s_fwd_d_min_0(&_S410, &_S411, &_S412);
    float _S413 = _S384.primal_0.differential_0 + _S386 + _S388;
    float _S414 = _S413 / _S314;
    float _S415 = ((_S384.differential_0.differential_0 + _S387 + _S389) * _S314 - _S413 * _S316) / _S392;
    float _S416 = - _S414;
    float _S417 = _S313.primal_0 * _S416;
    float _S418 = _S313.differential_0 * _S416 + - _S415 * _S313.primal_0;
    float _S419 = s_2 * _S414;
    float _S420 = s_diff_s_0 * _S414 + _S415 * s_2;
    DiffPair_float_0 _S421 = { _S310.primal_0, 0.0f };
    DiffPair_float_0 _S422 = { _S310.differential_0, 0.0f };
    DiffPair_0 _S423;
    (&_S423)->primal_0 = _S421;
    (&_S423)->differential_0 = _S422;
    DiffPair_0 _S424;
    (&_S424)->primal_0 = _S401;
    (&_S424)->differential_0 = _S402;
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = _S419;
    (&_S425)->differential_0 = _S420;
    s_fwd_d_max_0(&_S423, &_S424, &_S425);
    DiffPair_0 _S426;
    (&_S426)->primal_0 = _S406;
    (&_S426)->differential_0 = _S407;
    DiffPair_0 _S427;
    (&_S427)->primal_0 = _S408;
    (&_S427)->differential_0 = _S409;
    DiffPair_float_0 _S428;
    (&_S428)->primal_0 = _S423.primal_0.differential_0;
    (&_S428)->differential_0 = _S423.differential_0.differential_0;
    s_fwd_d_max_0(&_S426, &_S427, &_S428);
    float _S429 = _S395 + _S417;
    float _S430 = _S396 + _S418;
    float3  _S431 = make_float3 (_S410.primal_0.differential_0 + _S426.primal_0.differential_0 + _S429, _S411.primal_0.differential_0 + _S427.primal_0.differential_0 + _S429, _S404.primal_0.differential_0 + _S424.primal_0.differential_0 + _S429);
    float3  _S432 = make_float3 (_S410.differential_0.differential_0 + _S426.differential_0.differential_0 + _S430, _S411.differential_0.differential_0 + _S427.differential_0.differential_0 + _S430, _S404.differential_0.differential_0 + _S424.differential_0.differential_0 + _S430);
    float3  _S433 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S434 = { _S304, _S433 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S435 = { _S305, _S433 };
    DiffPair_1 _S436;
    (&_S436)->primal_0 = _S434;
    (&_S436)->differential_0 = _S435;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S437;
    (&_S437)->primal_0 = _S431;
    (&_S437)->differential_0 = _S432;
    s_fwd_s_bwd_prop_exp_1(&_S436, &_S437);
    float3  _S438 = make_float3 (2.0f) * _S436.primal_0.differential_0;
    float3  _S439 = _S436.differential_0.differential_0 * make_float3 (2.0f);
    float s_diff_scale_reg_T_1 = scale_regularization_weight_2 * (*v_loss_1)[int(2)];
    DiffPair_float_0 _S440 = { _S299, 0.0f };
    DiffPair_float_0 _S441 = { _S301, 0.0f };
    DiffPair_float_0 _S442 = { max_gauss_ratio_2, 0.0f };
    DiffPair_0 _S443;
    (&_S443)->primal_0 = _S440;
    (&_S443)->differential_0 = _S441;
    DiffPair_0 _S444;
    (&_S444)->primal_0 = _S442;
    (&_S444)->differential_0 = _S346;
    DiffPair_float_0 _S445;
    (&_S445)->primal_0 = s_diff_scale_reg_T_1;
    (&_S445)->differential_0 = 0.0f;
    s_fwd_d_max_0(&_S443, &_S444, &_S445);
    float _S446 = _S443.primal_0.differential_0 / _S300;
    float _S447 = (_S443.differential_0.differential_0 * _S300 - _S443.primal_0.differential_0 * _S303) / (_S300 * _S300);
    float _S448 = - _S446;
    float _S449 = _S292.primal_0 * _S448;
    float _S450 = _S292.differential_0 * _S448 + - _S447 * _S292.primal_0;
    float _S451 = _S298.primal_0 * _S446;
    float _S452 = _S298.differential_0 * _S446 + _S447 * _S298.primal_0;
    DiffPair_float_0 _S453 = { _S295.primal_0, 0.0f };
    DiffPair_float_0 _S454 = { _S295.differential_0, 0.0f };
    DiffPair_float_0 _S455 = { _S285, 0.0f };
    DiffPair_float_0 _S456 = { _S286, 0.0f };
    DiffPair_0 _S457;
    (&_S457)->primal_0 = _S453;
    (&_S457)->differential_0 = _S454;
    DiffPair_0 _S458;
    (&_S458)->primal_0 = _S455;
    (&_S458)->differential_0 = _S456;
    DiffPair_float_0 _S459;
    (&_S459)->primal_0 = _S449;
    (&_S459)->differential_0 = _S450;
    s_fwd_d_min_0(&_S457, &_S458, &_S459);
    DiffPair_float_0 _S460 = { _S281, 0.0f };
    DiffPair_float_0 _S461 = { _S282, 0.0f };
    DiffPair_float_0 _S462 = { _S283, 0.0f };
    DiffPair_float_0 _S463 = { _S284, 0.0f };
    DiffPair_0 _S464;
    (&_S464)->primal_0 = _S460;
    (&_S464)->differential_0 = _S461;
    DiffPair_0 _S465;
    (&_S465)->primal_0 = _S462;
    (&_S465)->differential_0 = _S463;
    DiffPair_float_0 _S466;
    (&_S466)->primal_0 = _S457.primal_0.differential_0;
    (&_S466)->differential_0 = _S457.differential_0.differential_0;
    s_fwd_d_min_0(&_S464, &_S465, &_S466);
    DiffPair_float_0 _S467 = { _S289.primal_0, 0.0f };
    DiffPair_float_0 _S468 = { _S289.differential_0, 0.0f };
    DiffPair_0 _S469;
    (&_S469)->primal_0 = _S467;
    (&_S469)->differential_0 = _S468;
    DiffPair_0 _S470;
    (&_S470)->primal_0 = _S455;
    (&_S470)->differential_0 = _S456;
    DiffPair_float_0 _S471;
    (&_S471)->primal_0 = _S451;
    (&_S471)->differential_0 = _S452;
    s_fwd_d_max_0(&_S469, &_S470, &_S471);
    DiffPair_0 _S472;
    (&_S472)->primal_0 = _S460;
    (&_S472)->differential_0 = _S461;
    DiffPair_0 _S473;
    (&_S473)->primal_0 = _S462;
    (&_S473)->differential_0 = _S463;
    DiffPair_float_0 _S474;
    (&_S474)->primal_0 = _S469.primal_0.differential_0;
    (&_S474)->differential_0 = _S469.differential_0.differential_0;
    s_fwd_d_max_0(&_S472, &_S473, &_S474);
    float _S475 = mcmc_scale_reg_weight_2 * (0.3333333432674408f * (*v_loss_1)[int(1)]);
    float3  _S476 = make_float3 (_S464.primal_0.differential_0 + _S472.primal_0.differential_0 + _S475, _S465.primal_0.differential_0 + _S473.primal_0.differential_0 + _S475, _S458.primal_0.differential_0 + _S470.primal_0.differential_0 + _S475);
    float3  _S477 = make_float3 (_S464.differential_0.differential_0 + _S472.differential_0.differential_0, _S465.differential_0.differential_0 + _S473.differential_0.differential_0, _S458.differential_0.differential_0 + _S470.differential_0.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S478 = { dpscales_0->primal_0, _S433 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S479 = { dpscales_0->differential_0, _S433 };
    DiffPair_1 _S480;
    (&_S480)->primal_0 = _S478;
    (&_S480)->differential_0 = _S479;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S481;
    (&_S481)->primal_0 = _S476;
    (&_S481)->differential_0 = _S477;
    s_fwd_s_bwd_prop_exp_1(&_S480, &_S481);
    float s_diff_quat_norm_reg_T_1 = quat_norm_reg_weight_2 * (*v_loss_1)[int(4)];
    float _S482 = - s_diff_quat_norm_reg_T_1;
    DiffPair_float_0 _S483 = { _S278.primal_0, 0.0f };
    DiffPair_float_0 _S484 = { _S278.differential_0, 0.0f };
    DiffPair_0 _S485;
    (&_S485)->primal_0 = _S483;
    (&_S485)->differential_0 = _S484;
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = _S482;
    (&_S486)->differential_0 = -0.0f;
    s_fwd_s_bwd_prop_log_0(&_S485, &_S486);
    float _S487 = _S485.primal_0.differential_0 + s_diff_quat_norm_reg_T_1;
    float4  _S488 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S489 = { dpquat_0->primal_0, _S488 };
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S490 = { dpquat_0->differential_0, _S488 };
    DiffPair_2 _S491;
    (&_S491)->primal_0 = _S489;
    (&_S491)->differential_0 = _S490;
    DiffPair_float_0 _S492;
    (&_S492)->primal_0 = _S487;
    (&_S492)->differential_0 = _S485.differential_0.differential_0;
    s_fwd_s_bwd_length_impl_0(&_S491, &_S492);
    float s_diff_reg_T_0 = mcmc_opacity_reg_weight_2 * (*v_loss_1)[int(0)];
    float _S493 = - (s_diff_reg_T_0 / _S274);
    float _S494 = - ((0.0f - s_diff_reg_T_0 * _S276) / (_S274 * _S274));
    DiffPair_float_0 _S495 = { _S269, 0.0f };
    DiffPair_float_0 _S496 = { _S270, 0.0f };
    DiffPair_0 _S497;
    (&_S497)->primal_0 = _S495;
    (&_S497)->differential_0 = _S496;
    DiffPair_float_0 _S498;
    (&_S498)->primal_0 = _S493;
    (&_S498)->differential_0 = _S494;
    s_fwd_s_bwd_prop_exp_0(&_S497, &_S498);
    float _S499 = - _S497.primal_0.differential_0;
    float _S500 = - _S497.differential_0.differential_0;
    float3  _S501 = _S439 + _S480.differential_0.differential_0;
    dpv_scales_0->primal_0 = _S438 + _S480.primal_0.differential_0;
    dpv_scales_0->differential_0 = _S501;
    dpv_opacity_0->primal_0 = _S499;
    dpv_opacity_0->differential_0 = _S500;
    dpv_quat_0->primal_0 = _S491.primal_0.differential_0;
    dpv_quat_0->differential_0 = _S491.differential_0.differential_0;
    return;
}

inline __device__ void per_splat_losses_bwd(float3  scales_2, float opacity_2, float4  quat_2, FixedArray<float, 5>  v_loss_2, float3  * v_scales_1, float * v_opacity_1, float4  * v_quat_1, float3  * vr_scales_0, float * vr_opacity_0, float4  * vr_quat_0, float3  * h_scales_0, float * h_opacity_0, float4  * h_quat_0, float mcmc_opacity_reg_weight_3, float mcmc_scale_reg_weight_3, float max_gauss_ratio_3, float scale_regularization_weight_3, float erank_reg_weight_3, float erank_reg_weight_s3_3, float quat_norm_reg_weight_3)
{
    float _S502 = - opacity_2;
    float _S503 = 1.0f + s_primal_ctx_exp_0(_S502);
    float _S504 = _S503 * _S503;
    float _S505 = length_0(quat_2);
    float3  _S506 = s_primal_ctx_exp_1(scales_2);
    float _S507 = _S506.x;
    float _S508 = _S506.y;
    float _S509 = _S506.z;
    float _S510 = (F32_max((_S507), (_S508)));
    float _S511 = (F32_max((_S510), (_S509)));
    float _S512 = (F32_min((_S507), (_S508)));
    float _S513 = (F32_min((_S512), (_S509)));
    float _S514 = _S511 / _S513;
    float _S515 = _S513 * _S513;
    float3  _S516 = make_float3 (2.0f) * scales_2;
    float3  _S517 = s_primal_ctx_exp_1(_S516);
    float x_10 = _S517.x;
    float y_6 = _S517.y;
    float z_3 = _S517.z;
    float s_3 = x_10 + y_6 + z_3;
    float _S518 = (F32_max((x_10), (y_6)));
    float _S519 = (F32_max((_S518), (z_3)));
    float s1_3 = _S519 / s_3;
    float _S520 = s_3 * s_3;
    float _S521 = (F32_min((x_10), (y_6)));
    float _S522 = (F32_min((_S521), (z_3)));
    float s3_3 = _S522 / s_3;
    float s2_3 = 1.0f - s1_3 - s3_3;
    float _S523 = - s1_3;
    float _S524 = s_primal_ctx_log_0(s1_3);
    float _S525 = s_primal_ctx_log_0(s2_3);
    float _S526 = s_primal_ctx_log_0(s3_3);
    float _S527 = _S523 * _S524 - s2_3 * _S525 - s3_3 * _S526;
    float _S528 = s_primal_ctx_exp_0(_S527) - 0.99998998641967773f;
    float _S529 = erank_reg_weight_s3_3 * v_loss_2[int(3)];
    float _S530 = erank_reg_weight_3 * v_loss_2[int(3)];
    DiffPair_float_0 _S531;
    (&_S531)->primal_0 = - s_primal_ctx_log_0(_S528);
    (&_S531)->differential_0 = 0.0f;
    DiffPair_float_0 _S532;
    (&_S532)->primal_0 = 0.0f;
    (&_S532)->differential_0 = 0.0f;
    _d_max_0(&_S531, &_S532, _S530);
    float _S533 = - _S531.differential_0;
    DiffPair_float_0 _S534;
    (&_S534)->primal_0 = _S528;
    (&_S534)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S534, _S533);
    DiffPair_float_0 _S535;
    (&_S535)->primal_0 = _S527;
    (&_S535)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S535, _S534.differential_0);
    float _S536 = - _S535.differential_0;
    float _S537 = s3_3 * _S536;
    float _S538 = _S526 * _S536;
    DiffPair_float_0 _S539;
    (&_S539)->primal_0 = s3_3;
    (&_S539)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S539, _S537);
    float _S540 = s2_3 * _S536;
    float _S541 = _S525 * _S536;
    DiffPair_float_0 _S542;
    (&_S542)->primal_0 = s2_3;
    (&_S542)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S542, _S540);
    float _S543 = _S523 * _S535.differential_0;
    float _S544 = _S524 * _S535.differential_0;
    DiffPair_float_0 _S545;
    (&_S545)->primal_0 = s1_3;
    (&_S545)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S545, _S543);
    float _S546 = - _S544;
    float _S547 = - (_S541 + _S542.differential_0);
    float _S548 = (_S529 + _S538 + _S539.differential_0 + _S547) / _S520;
    float _S549 = _S522 * - _S548;
    float _S550 = s_3 * _S548;
    DiffPair_float_0 _S551;
    (&_S551)->primal_0 = _S521;
    (&_S551)->differential_0 = 0.0f;
    DiffPair_float_0 _S552;
    (&_S552)->primal_0 = z_3;
    (&_S552)->differential_0 = 0.0f;
    _d_min_0(&_S551, &_S552, _S550);
    DiffPair_float_0 _S553;
    (&_S553)->primal_0 = x_10;
    (&_S553)->differential_0 = 0.0f;
    DiffPair_float_0 _S554;
    (&_S554)->primal_0 = y_6;
    (&_S554)->differential_0 = 0.0f;
    _d_min_0(&_S553, &_S554, _S551.differential_0);
    float _S555 = (_S545.differential_0 + _S546 + _S547) / _S520;
    float _S556 = _S519 * - _S555;
    float _S557 = s_3 * _S555;
    DiffPair_float_0 _S558;
    (&_S558)->primal_0 = _S518;
    (&_S558)->differential_0 = 0.0f;
    DiffPair_float_0 _S559;
    (&_S559)->primal_0 = z_3;
    (&_S559)->differential_0 = 0.0f;
    _d_max_0(&_S558, &_S559, _S557);
    DiffPair_float_0 _S560;
    (&_S560)->primal_0 = x_10;
    (&_S560)->differential_0 = 0.0f;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = y_6;
    (&_S561)->differential_0 = 0.0f;
    _d_max_0(&_S560, &_S561, _S558.differential_0);
    float _S562 = _S549 + _S556;
    float3  _S563 = make_float3 (_S553.differential_0 + _S560.differential_0 + _S562, _S554.differential_0 + _S561.differential_0 + _S562, _S552.differential_0 + _S559.differential_0 + _S562);
    float3  _S564 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S565;
    (&_S565)->primal_0 = _S516;
    (&_S565)->differential_0 = _S564;
    s_bwd_prop_exp_1(&_S565, _S563);
    float3  _S566 = make_float3 (2.0f) * _S565.differential_0;
    float s_diff_scale_reg_T_2 = scale_regularization_weight_3 * v_loss_2[int(2)];
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S514;
    (&_S567)->differential_0 = 0.0f;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = max_gauss_ratio_3;
    (&_S568)->differential_0 = 0.0f;
    _d_max_0(&_S567, &_S568, s_diff_scale_reg_T_2);
    float _S569 = _S567.differential_0 / _S515;
    float _S570 = _S511 * - _S569;
    float _S571 = _S513 * _S569;
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S512;
    (&_S572)->differential_0 = 0.0f;
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S509;
    (&_S573)->differential_0 = 0.0f;
    _d_min_0(&_S572, &_S573, _S570);
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = _S507;
    (&_S574)->differential_0 = 0.0f;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = _S508;
    (&_S575)->differential_0 = 0.0f;
    _d_min_0(&_S574, &_S575, _S572.differential_0);
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S510;
    (&_S576)->differential_0 = 0.0f;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = _S509;
    (&_S577)->differential_0 = 0.0f;
    _d_max_0(&_S576, &_S577, _S571);
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = _S507;
    (&_S578)->differential_0 = 0.0f;
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S508;
    (&_S579)->differential_0 = 0.0f;
    _d_max_0(&_S578, &_S579, _S576.differential_0);
    float _S580 = mcmc_scale_reg_weight_3 * (0.3333333432674408f * v_loss_2[int(1)]);
    float3  _S581 = make_float3 (_S574.differential_0 + _S578.differential_0 + _S580, _S575.differential_0 + _S579.differential_0 + _S580, _S573.differential_0 + _S577.differential_0 + _S580);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S582;
    (&_S582)->primal_0 = scales_2;
    (&_S582)->differential_0 = _S564;
    s_bwd_prop_exp_1(&_S582, _S581);
    float s_diff_quat_norm_reg_T_2 = quat_norm_reg_weight_3 * v_loss_2[int(4)];
    float _S583 = - s_diff_quat_norm_reg_T_2;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S505;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S584, _S583);
    float _S585 = _S584.differential_0 + s_diff_quat_norm_reg_T_2;
    float4  _S586 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S587;
    (&_S587)->primal_0 = quat_2;
    (&_S587)->differential_0 = _S586;
    s_bwd_length_impl_0(&_S587, _S585);
    float _S588 = - (mcmc_opacity_reg_weight_3 * v_loss_2[int(0)] / _S504);
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = _S502;
    (&_S589)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S589, _S588);
    float _S590 = - _S589.differential_0;
    *v_scales_1 = _S566 + _S582.differential_0;
    *v_opacity_1 = _S590;
    *v_quat_1 = _S587.differential_0;
    FixedArray<float, 5>  losses_1;
    losses_1[int(0)] = mcmc_opacity_reg_weight_3 * (1.0f / (1.0f + (F32_exp((_S502)))));
    losses_1[int(4)] = quat_norm_reg_weight_3 * (_S505 - 1.0f - (F32_log((_S505))));
    float3  _S591 = exp_0(scales_2);
    float _S592 = _S591.x;
    float _S593 = _S591.y;
    float _S594 = _S591.z;
    losses_1[int(1)] = mcmc_scale_reg_weight_3 * (_S592 + _S593 + _S594) / 3.0f;
    losses_1[int(2)] = scale_regularization_weight_3 * ((F32_max(((F32_max(((F32_max((_S592), (_S593)))), (_S594))) / (F32_min(((F32_min((_S592), (_S593)))), (_S594)))), (max_gauss_ratio_3))) - max_gauss_ratio_3);
    float3  _S595 = exp_0(_S516);
    float x_11 = _S595.x;
    float y_7 = _S595.y;
    float z_4 = _S595.z;
    float s_4 = x_11 + y_7 + z_4;
    float s1_4 = (F32_max(((F32_max((x_11), (y_7)))), (z_4))) / s_4;
    float s3_4 = (F32_min(((F32_min((x_11), (y_7)))), (z_4))) / s_4;
    float s2_4 = 1.0f - s1_4 - s3_4;
    losses_1[int(3)] = erank_reg_weight_3 * (F32_max((- (F32_log(((F32_exp((- s1_4 * (F32_log((s1_4))) - s2_4 * (F32_log((s2_4))) - s3_4 * (F32_log((s3_4)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_3 * s3_4;
    FixedArray<float, 5>  residuals_0;
    float _S596 = (F32_sqrt(((F32_max((losses_1[int(0)]), (0.0f))))));
    residuals_0[int(0)] = _S596;
    float _S597 = (F32_sqrt(((F32_max((losses_1[int(1)]), (0.0f))))));
    residuals_0[int(1)] = _S597;
    float _S598 = (F32_sqrt(((F32_max((losses_1[int(2)]), (0.0f))))));
    residuals_0[int(2)] = _S598;
    float _S599 = (F32_sqrt(((F32_max((losses_1[int(3)]), (0.0f))))));
    residuals_0[int(3)] = _S599;
    float _S600 = (F32_sqrt(((F32_max((losses_1[int(4)]), (0.0f))))));
    residuals_0[int(4)] = _S600;
    float _S601 = _S511 / _S513;
    float s1_5 = _S519 / s_3;
    float s3_5 = _S522 / s_3;
    float s2_5 = 1.0f - s1_5 - s3_5;
    float _S602 = - s1_5;
    float _S603 = s_primal_ctx_log_0(s1_5);
    float _S604 = s_primal_ctx_log_0(s2_5);
    float _S605 = s_primal_ctx_log_0(s3_5);
    float _S606 = _S602 * _S603 - s2_5 * _S604 - s3_5 * _S605;
    float _S607 = s_primal_ctx_exp_0(_S606) - 0.99998998641967773f;
    float _S608 = erank_reg_weight_s3_3 * _S599;
    float _S609 = erank_reg_weight_3 * _S599;
    DiffPair_float_0 _S610;
    (&_S610)->primal_0 = - s_primal_ctx_log_0(_S607);
    (&_S610)->differential_0 = 0.0f;
    DiffPair_float_0 _S611;
    (&_S611)->primal_0 = 0.0f;
    (&_S611)->differential_0 = 0.0f;
    _d_max_0(&_S610, &_S611, _S609);
    float _S612 = - _S610.differential_0;
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = _S607;
    (&_S613)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S613, _S612);
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S606;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S614, _S613.differential_0);
    float _S615 = - _S614.differential_0;
    float _S616 = s3_5 * _S615;
    float _S617 = _S605 * _S615;
    DiffPair_float_0 _S618;
    (&_S618)->primal_0 = s3_5;
    (&_S618)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S618, _S616);
    float _S619 = s2_5 * _S615;
    float _S620 = _S604 * _S615;
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = s2_5;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S621, _S619);
    float _S622 = _S602 * _S614.differential_0;
    float _S623 = _S603 * _S614.differential_0;
    DiffPair_float_0 _S624;
    (&_S624)->primal_0 = s1_5;
    (&_S624)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S624, _S622);
    float _S625 = - _S623;
    float _S626 = - (_S620 + _S621.differential_0);
    float _S627 = (_S608 + _S617 + _S618.differential_0 + _S626) / _S520;
    float _S628 = _S522 * - _S627;
    float _S629 = s_3 * _S627;
    DiffPair_float_0 _S630;
    (&_S630)->primal_0 = _S521;
    (&_S630)->differential_0 = 0.0f;
    DiffPair_float_0 _S631;
    (&_S631)->primal_0 = z_3;
    (&_S631)->differential_0 = 0.0f;
    _d_min_0(&_S630, &_S631, _S629);
    DiffPair_float_0 _S632;
    (&_S632)->primal_0 = x_10;
    (&_S632)->differential_0 = 0.0f;
    DiffPair_float_0 _S633;
    (&_S633)->primal_0 = y_6;
    (&_S633)->differential_0 = 0.0f;
    _d_min_0(&_S632, &_S633, _S630.differential_0);
    float _S634 = (_S624.differential_0 + _S625 + _S626) / _S520;
    float _S635 = _S519 * - _S634;
    float _S636 = s_3 * _S634;
    DiffPair_float_0 _S637;
    (&_S637)->primal_0 = _S518;
    (&_S637)->differential_0 = 0.0f;
    DiffPair_float_0 _S638;
    (&_S638)->primal_0 = z_3;
    (&_S638)->differential_0 = 0.0f;
    _d_max_0(&_S637, &_S638, _S636);
    DiffPair_float_0 _S639;
    (&_S639)->primal_0 = x_10;
    (&_S639)->differential_0 = 0.0f;
    DiffPair_float_0 _S640;
    (&_S640)->primal_0 = y_6;
    (&_S640)->differential_0 = 0.0f;
    _d_max_0(&_S639, &_S640, _S637.differential_0);
    float _S641 = _S628 + _S635;
    float3  _S642 = make_float3 (_S632.differential_0 + _S639.differential_0 + _S641, _S633.differential_0 + _S640.differential_0 + _S641, _S631.differential_0 + _S638.differential_0 + _S641);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S643;
    (&_S643)->primal_0 = _S516;
    (&_S643)->differential_0 = _S564;
    s_bwd_prop_exp_1(&_S643, _S642);
    float3  _S644 = make_float3 (2.0f) * _S643.differential_0;
    float s_diff_scale_reg_T_3 = scale_regularization_weight_3 * _S598;
    DiffPair_float_0 _S645;
    (&_S645)->primal_0 = _S601;
    (&_S645)->differential_0 = 0.0f;
    DiffPair_float_0 _S646;
    (&_S646)->primal_0 = max_gauss_ratio_3;
    (&_S646)->differential_0 = 0.0f;
    _d_max_0(&_S645, &_S646, s_diff_scale_reg_T_3);
    float _S647 = _S645.differential_0 / _S515;
    float _S648 = _S511 * - _S647;
    float _S649 = _S513 * _S647;
    DiffPair_float_0 _S650;
    (&_S650)->primal_0 = _S512;
    (&_S650)->differential_0 = 0.0f;
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = _S509;
    (&_S651)->differential_0 = 0.0f;
    _d_min_0(&_S650, &_S651, _S648);
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = _S507;
    (&_S652)->differential_0 = 0.0f;
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = _S508;
    (&_S653)->differential_0 = 0.0f;
    _d_min_0(&_S652, &_S653, _S650.differential_0);
    DiffPair_float_0 _S654;
    (&_S654)->primal_0 = _S510;
    (&_S654)->differential_0 = 0.0f;
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = _S509;
    (&_S655)->differential_0 = 0.0f;
    _d_max_0(&_S654, &_S655, _S649);
    DiffPair_float_0 _S656;
    (&_S656)->primal_0 = _S507;
    (&_S656)->differential_0 = 0.0f;
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = _S508;
    (&_S657)->differential_0 = 0.0f;
    _d_max_0(&_S656, &_S657, _S654.differential_0);
    float _S658 = mcmc_scale_reg_weight_3 * (0.3333333432674408f * _S597);
    float3  _S659 = make_float3 (_S652.differential_0 + _S656.differential_0 + _S658, _S653.differential_0 + _S657.differential_0 + _S658, _S651.differential_0 + _S655.differential_0 + _S658);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S660;
    (&_S660)->primal_0 = scales_2;
    (&_S660)->differential_0 = _S564;
    s_bwd_prop_exp_1(&_S660, _S659);
    float s_diff_quat_norm_reg_T_3 = quat_norm_reg_weight_3 * _S600;
    float _S661 = - s_diff_quat_norm_reg_T_3;
    DiffPair_float_0 _S662;
    (&_S662)->primal_0 = _S505;
    (&_S662)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S662, _S661);
    float _S663 = _S662.differential_0 + s_diff_quat_norm_reg_T_3;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S664;
    (&_S664)->primal_0 = quat_2;
    (&_S664)->differential_0 = _S586;
    s_bwd_length_impl_0(&_S664, _S663);
    float _S665 = - (mcmc_opacity_reg_weight_3 * _S596 / _S504);
    DiffPair_float_0 _S666;
    (&_S666)->primal_0 = _S502;
    (&_S666)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S666, _S665);
    float _S667 = - _S666.differential_0;
    *vr_scales_0 = _S644 + _S660.differential_0;
    *vr_opacity_0 = _S667;
    *vr_quat_0 = _S664.differential_0;
    float3  _S668 = make_float3 (1.0f, 0.0f, 0.0f);
    int4  _S669 = make_int4 (int(0));
    float4  _S670 = make_float4 ((float)_S669.x, (float)_S669.y, (float)_S669.z, (float)_S669.w);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S671;
    (&_S671)->primal_0 = scales_2;
    (&_S671)->differential_0 = _S668;
    DiffPair_float_0 _S672;
    (&_S672)->primal_0 = opacity_2;
    (&_S672)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S673;
    (&_S673)->primal_0 = quat_2;
    (&_S673)->differential_0 = _S670;
    FixedArray<float, 5>  _S674;
    _S674[int(0)] = 1.0f;
    _S674[int(1)] = 1.0f;
    _S674[int(2)] = 1.0f;
    _S674[int(3)] = 1.0f;
    _S674[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_0;
    DiffPair_float_0 dp_v_opacity_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_0;
    s_fwd_per_splat_losses_bwd_0(&_S671, &_S672, &_S673, &_S674, &dp_v_scales_0, &dp_v_opacity_0, &dp_v_quat_0, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->x) = dp_v_scales_0.differential_0.x;
    float3  _S675 = make_float3 (0.0f, 1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S676;
    (&_S676)->primal_0 = scales_2;
    (&_S676)->differential_0 = _S675;
    DiffPair_float_0 _S677;
    (&_S677)->primal_0 = opacity_2;
    (&_S677)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S678;
    (&_S678)->primal_0 = quat_2;
    (&_S678)->differential_0 = _S670;
    FixedArray<float, 5>  _S679;
    _S679[int(0)] = 1.0f;
    _S679[int(1)] = 1.0f;
    _S679[int(2)] = 1.0f;
    _S679[int(3)] = 1.0f;
    _S679[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_1;
    DiffPair_float_0 dp_v_opacity_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_1;
    s_fwd_per_splat_losses_bwd_0(&_S676, &_S677, &_S678, &_S679, &dp_v_scales_1, &dp_v_opacity_1, &dp_v_quat_1, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->y) = dp_v_scales_1.differential_0.y;
    float3  _S680 = make_float3 (0.0f, 0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S681;
    (&_S681)->primal_0 = scales_2;
    (&_S681)->differential_0 = _S680;
    DiffPair_float_0 _S682;
    (&_S682)->primal_0 = opacity_2;
    (&_S682)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S683;
    (&_S683)->primal_0 = quat_2;
    (&_S683)->differential_0 = _S670;
    FixedArray<float, 5>  _S684;
    _S684[int(0)] = 1.0f;
    _S684[int(1)] = 1.0f;
    _S684[int(2)] = 1.0f;
    _S684[int(3)] = 1.0f;
    _S684[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_2;
    DiffPair_float_0 dp_v_opacity_2;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_2;
    s_fwd_per_splat_losses_bwd_0(&_S681, &_S682, &_S683, &_S684, &dp_v_scales_2, &dp_v_opacity_2, &dp_v_quat_2, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->z) = dp_v_scales_2.differential_0.z;
    int3  _S685 = make_int3 (int(0));
    float3  _S686 = make_float3 ((float)_S685.x, (float)_S685.y, (float)_S685.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S687;
    (&_S687)->primal_0 = scales_2;
    (&_S687)->differential_0 = _S686;
    DiffPair_float_0 _S688;
    (&_S688)->primal_0 = opacity_2;
    (&_S688)->differential_0 = 1.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S689;
    (&_S689)->primal_0 = quat_2;
    (&_S689)->differential_0 = _S670;
    FixedArray<float, 5>  _S690;
    _S690[int(0)] = 1.0f;
    _S690[int(1)] = 1.0f;
    _S690[int(2)] = 1.0f;
    _S690[int(3)] = 1.0f;
    _S690[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_3;
    DiffPair_float_0 dp_v_opacity_3;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_3;
    s_fwd_per_splat_losses_bwd_0(&_S687, &_S688, &_S689, &_S690, &dp_v_scales_3, &dp_v_opacity_3, &dp_v_quat_3, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *h_opacity_0 = dp_v_opacity_3.differential_0;
    float4  _S691 = make_float4 (1.0f, 0.0f, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S692;
    (&_S692)->primal_0 = scales_2;
    (&_S692)->differential_0 = _S686;
    DiffPair_float_0 _S693;
    (&_S693)->primal_0 = opacity_2;
    (&_S693)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S694;
    (&_S694)->primal_0 = quat_2;
    (&_S694)->differential_0 = _S691;
    FixedArray<float, 5>  _S695;
    _S695[int(0)] = 1.0f;
    _S695[int(1)] = 1.0f;
    _S695[int(2)] = 1.0f;
    _S695[int(3)] = 1.0f;
    _S695[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_4;
    DiffPair_float_0 dp_v_opacity_4;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_4;
    s_fwd_per_splat_losses_bwd_0(&_S692, &_S693, &_S694, &_S695, &dp_v_scales_4, &dp_v_opacity_4, &dp_v_quat_4, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->x) = dp_v_quat_4.differential_0.x;
    float4  _S696 = make_float4 (0.0f, 1.0f, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S697;
    (&_S697)->primal_0 = scales_2;
    (&_S697)->differential_0 = _S686;
    DiffPair_float_0 _S698;
    (&_S698)->primal_0 = opacity_2;
    (&_S698)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S699;
    (&_S699)->primal_0 = quat_2;
    (&_S699)->differential_0 = _S696;
    FixedArray<float, 5>  _S700;
    _S700[int(0)] = 1.0f;
    _S700[int(1)] = 1.0f;
    _S700[int(2)] = 1.0f;
    _S700[int(3)] = 1.0f;
    _S700[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_5;
    DiffPair_float_0 dp_v_opacity_5;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_5;
    s_fwd_per_splat_losses_bwd_0(&_S697, &_S698, &_S699, &_S700, &dp_v_scales_5, &dp_v_opacity_5, &dp_v_quat_5, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->y) = dp_v_quat_5.differential_0.y;
    float4  _S701 = make_float4 (0.0f, 0.0f, 1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S702;
    (&_S702)->primal_0 = scales_2;
    (&_S702)->differential_0 = _S686;
    DiffPair_float_0 _S703;
    (&_S703)->primal_0 = opacity_2;
    (&_S703)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S704;
    (&_S704)->primal_0 = quat_2;
    (&_S704)->differential_0 = _S701;
    FixedArray<float, 5>  _S705;
    _S705[int(0)] = 1.0f;
    _S705[int(1)] = 1.0f;
    _S705[int(2)] = 1.0f;
    _S705[int(3)] = 1.0f;
    _S705[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_6;
    DiffPair_float_0 dp_v_opacity_6;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_6;
    s_fwd_per_splat_losses_bwd_0(&_S702, &_S703, &_S704, &_S705, &dp_v_scales_6, &dp_v_opacity_6, &dp_v_quat_6, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->z) = dp_v_quat_6.differential_0.z;
    float4  _S706 = make_float4 (0.0f, 0.0f, 0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S707;
    (&_S707)->primal_0 = scales_2;
    (&_S707)->differential_0 = _S686;
    DiffPair_float_0 _S708;
    (&_S708)->primal_0 = opacity_2;
    (&_S708)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S709;
    (&_S709)->primal_0 = quat_2;
    (&_S709)->differential_0 = _S706;
    FixedArray<float, 5>  _S710;
    _S710[int(0)] = 1.0f;
    _S710[int(1)] = 1.0f;
    _S710[int(2)] = 1.0f;
    _S710[int(3)] = 1.0f;
    _S710[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_7;
    DiffPair_float_0 dp_v_opacity_7;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_7;
    s_fwd_per_splat_losses_bwd_0(&_S707, &_S708, &_S709, &_S710, &dp_v_scales_7, &dp_v_opacity_7, &dp_v_quat_7, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->w) = dp_v_quat_7.differential_0.w;
    return;
}

inline __device__ float4  normalize_0(float4  x_12)
{
    return x_12 / make_float4 (length_0(x_12));
}

inline __device__ float3  normalize_1(float3  x_13)
{
    return x_13 / make_float3 (length_1(x_13));
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_14)
{
    Matrix<float, 3, 3>  result_13;
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
            *_slang_vector_get_element_ptr(((&result_13)->rows + (r_1)), c_0) = _slang_vector_get_element(x_14.rows[c_0], r_1);
            c_0 = c_0 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_13;
}

inline __device__ Matrix<float, 2, 2>  transpose_1(Matrix<float, 2, 2>  x_15)
{
    Matrix<float, 2, 2>  result_14;
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
            *_slang_vector_get_element_ptr(((&result_14)->rows + (r_2)), c_1) = _slang_vector_get_element(x_15.rows[c_1], r_2);
            c_1 = c_1 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_14;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_0, Matrix<float, 3, 3>  dOut_7)
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
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_7.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(0)].x * dOut_7.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_7.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(0)].y * dOut_7.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_7.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(0)].z * dOut_7.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_7.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(0)].x * dOut_7.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_7.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(0)].y * dOut_7.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_7.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(0)].z * dOut_7.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_7.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(0)].x * dOut_7.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_7.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(0)].y * dOut_7.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_7.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(0)].z * dOut_7.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_7.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(1)].x * dOut_7.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_7.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(1)].y * dOut_7.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_7.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(1)].z * dOut_7.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_7.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(1)].x * dOut_7.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_7.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(1)].y * dOut_7.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_7.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(1)].z * dOut_7.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_7.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(1)].x * dOut_7.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_7.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(1)].y * dOut_7.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_7.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(1)].z * dOut_7.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_7.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(2)].x * dOut_7.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_7.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(2)].y * dOut_7.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_7.rows[int(2)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(2)].z * dOut_7.rows[int(2)].x;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_7.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(2)].x * dOut_7.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_7.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(2)].y * dOut_7.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_7.rows[int(2)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(2)].z * dOut_7.rows[int(2)].y;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = *&(((&left_d_result_0)->rows + (int(2)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_7.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(2)].x * dOut_7.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = *&(((&left_d_result_0)->rows + (int(2)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_7.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(2)].y * dOut_7.rows[int(2)].z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = *&(((&left_d_result_0)->rows + (int(2)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_7.rows[int(2)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(2)].z * dOut_7.rows[int(2)].z;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
{
    Matrix<float, 3, 3>  result_15;
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
            int i_4 = int(0);
            float sum_0 = 0.0f;
            for(;;)
            {
                if(i_4 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_1 = sum_0 + _slang_vector_get_element(left_1.rows[r_3], i_4) * _slang_vector_get_element(right_1.rows[i_4], c_2);
                i_4 = i_4 + int(1);
                sum_0 = sum_1;
            }
            *_slang_vector_get_element_ptr(((&result_15)->rows + (r_3)), c_2) = sum_0;
            c_2 = c_2 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_15;
}

inline __device__ float4  floor_0(float4  x_16)
{
    float4  result_16;
    int i_5 = int(0);
    for(;;)
    {
        if(i_5 < int(4))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_16, i_5) = (F32_floor((_slang_vector_get_element(x_16, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_16;
}

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_2, float3  dOut_8)
{
    float _S711 = (*left_2).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_2).primal_0.x * dOut_8.x;
    float sum_2 = _S711 + (*left_2).primal_0.rows[int(1)].x * dOut_8.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_2).primal_0.x * dOut_8.y;
    float sum_3 = sum_2 + (*left_2).primal_0.rows[int(2)].x * dOut_8.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_2).primal_0.x * dOut_8.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_3;
    float _S712 = (*left_2).primal_0.rows[int(0)].y * dOut_8.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_2).primal_0.y * dOut_8.x;
    float sum_4 = _S712 + (*left_2).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_2).primal_0.y * dOut_8.y;
    float sum_5 = sum_4 + (*left_2).primal_0.rows[int(2)].y * dOut_8.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_2).primal_0.y * dOut_8.z;
    *&((&right_d_result_1)->y) = sum_5;
    float _S713 = (*left_2).primal_0.rows[int(0)].z * dOut_8.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_2).primal_0.z * dOut_8.x;
    float sum_6 = _S713 + (*left_2).primal_0.rows[int(1)].z * dOut_8.y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = (*right_2).primal_0.z * dOut_8.y;
    float sum_7 = sum_6 + (*left_2).primal_0.rows[int(2)].z * dOut_8.z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = (*right_2).primal_0.z * dOut_8.z;
    *&((&right_d_result_1)->z) = sum_7;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

struct DiffPair_matrixx3Cfloatx2C2x2C2x3E_0
{
    Matrix<float, 2, 2>  primal_0;
    Matrix<float, 2, 2>  differential_0;
};

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * left_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * right_3, float2  dOut_9)
{
    float _S714 = (*left_3).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 2, 2>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = (*right_3).primal_0.x * dOut_9.x;
    float sum_8 = _S714 + (*left_3).primal_0.rows[int(1)].x * dOut_9.y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = (*right_3).primal_0.x * dOut_9.y;
    float2  right_d_result_2;
    *&((&right_d_result_2)->x) = sum_8;
    float _S715 = (*left_3).primal_0.rows[int(0)].y * dOut_9.x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = (*right_3).primal_0.y * dOut_9.x;
    float sum_9 = _S715 + (*left_3).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = (*right_3).primal_0.y * dOut_9.y;
    *&((&right_d_result_2)->y) = sum_9;
    left_3->primal_0 = (*left_3).primal_0;
    left_3->differential_0 = left_d_result_2;
    right_3->primal_0 = (*right_3).primal_0;
    right_3->differential_0 = right_d_result_2;
    return;
}

inline __device__ float3  mul_2(Matrix<float, 3, 3>  left_4, float3  right_4)
{
    float3  result_17;
    int i_6 = int(0);
    for(;;)
    {
        if(i_6 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_4.rows[i_6], j_0) * _slang_vector_get_element(right_4, j_0);
            j_0 = j_0 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_17, i_6) = sum_10;
        i_6 = i_6 + int(1);
    }
    return result_17;
}

inline __device__ float2  mul_3(Matrix<float, 2, 2>  left_5, float2  right_5)
{
    float2  result_18;
    int i_7 = int(0);
    for(;;)
    {
        if(i_7 < int(2))
        {
        }
        else
        {
            break;
        }
        int j_1 = int(0);
        float sum_12 = 0.0f;
        for(;;)
        {
            if(j_1 < int(2))
            {
            }
            else
            {
                break;
            }
            float sum_13 = sum_12 + _slang_vector_get_element(left_5.rows[i_7], j_1) * _slang_vector_get_element(right_5, j_1);
            j_1 = j_1 + int(1);
            sum_12 = sum_13;
        }
        *_slang_vector_get_element_ptr(&result_18, i_7) = sum_12;
        i_7 = i_7 + int(1);
    }
    return result_18;
}

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_0, float3  scale_0, float4  quat_3, float opac_0)
{
    float4  _S716 = normalize_0(quat_3);
    float3  _S717 = exp_0(scale_0);
    float x_17 = _S716.y;
    float x2_0 = x_17 * x_17;
    float y2_0 = _S716.z * _S716.z;
    float z2_0 = _S716.w * _S716.w;
    float xy_0 = _S716.y * _S716.z;
    float xz_0 = _S716.y * _S716.w;
    float yz_0 = _S716.z * _S716.w;
    float wx_0 = _S716.x * _S716.y;
    float wy_0 = _S716.x * _S716.z;
    float wz_0 = _S716.x * _S716.w;
    Matrix<float, 3, 3>  M_0 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0))), makeMatrix<float, 3, 3> (_S717.x, 0.0f, 0.0f, 0.0f, _S717.y, 0.0f, 0.0f, 0.0f, _S717.z));
    float4  _S718 = make_float4 (dot_0(*mean_0, *mean_0), dot_0(*mean_0, scale_0), dot_0(scale_0, scale_0), dot_1(quat_3, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S719 = _S718 - floor_0(_S718);
    float4  _S720 = _S719 + make_float4 (dot_1(_S719, float4 {_S719.w, _S719.z, _S719.x, _S719.y} + make_float4 (33.3300018310546875f)));
    float4  _S721 = (float4 {_S720.x, _S720.x, _S720.y, _S720.z} + float4 {_S720.y, _S720.z, _S720.z, _S720.w}) * float4 {_S720.z, _S720.y, _S720.w, _S720.x};
    float4  _S722 = _S721 - floor_0(_S721);
    float2  _S723 = float2 {_S722.x, _S722.z};
    float _S724 = 6.28318548202514648f * _S723.y;
    float2  _S725 = float2 {_S722.y, _S722.w};
    float _S726 = 6.28318548202514648f * _S725.y;
    *mean_0 = *mean_0 + mul_2(mul_1(M_0, transpose_0(M_0)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S723.x))))))) * make_float2 ((F32_cos((_S724))), (F32_sin((_S724))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S723.x))))))) * make_float2 ((F32_cos((_S724))), (F32_sin((_S724))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S725.x))))))) * make_float2 ((F32_cos((_S726))), (F32_sin((_S726))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_10)
{
    float _S727 = dOut_10.y;
    float _S728 = dOut_10.z;
    float _S729 = dOut_10.x;
    float _S730 = (*a_0).primal_0.z * _S727 + - (*a_0).primal_0.y * _S728;
    float _S731 = - (*a_0).primal_0.z * _S729 + (*a_0).primal_0.x * _S728;
    float _S732 = (*a_0).primal_0.y * _S729 + - (*a_0).primal_0.x * _S727;
    float3  _S733 = make_float3 (- (*b_1).primal_0.z * _S727 + (*b_1).primal_0.y * _S728, (*b_1).primal_0.z * _S729 + - (*b_1).primal_0.x * _S728, - (*b_1).primal_0.y * _S729 + (*b_1).primal_0.x * _S727);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S733;
    float3  _S734 = make_float3 (_S730, _S731, _S732);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S734;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S735 = left_6.y;
    float _S736 = right_6.z;
    float _S737 = left_6.z;
    float _S738 = right_6.y;
    float _S739 = right_6.x;
    float _S740 = left_6.x;
    return make_float3 (_S735 * _S736 - _S737 * _S738, _S737 * _S739 - _S740 * _S736, _S740 * _S738 - _S735 * _S739);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_1, float3  scale_1, float4  quat_4, float opac_1)
{
    float4  _S741 = normalize_0(quat_4);
    float _S742 = scale_1.x;
    float sx_0 = (F32_exp((_S742)));
    float _S743 = scale_1.y;
    float sy_0 = (F32_exp((_S743)));
    float sz_0 = scale_1.z - 0.5f * (_S742 + _S743);
    float x_18 = _S741.y;
    float x2_1 = x_18 * x_18;
    float y2_1 = _S741.z * _S741.z;
    float z2_1 = _S741.w * _S741.w;
    float xy_1 = _S741.y * _S741.z;
    float xz_1 = _S741.y * _S741.w;
    float yz_1 = _S741.z * _S741.w;
    float wx_1 = _S741.x * _S741.y;
    float wy_1 = _S741.x * _S741.z;
    float wz_1 = _S741.x * _S741.w;
    Matrix<float, 3, 3>  _S744 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1)));
    float3  vert0_0 = mul_2(_S744, make_float3 (sx_0, 0.0f, 0.0f)) + *mean_1;
    float3  vert1_0 = mul_2(_S744, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + *mean_1;
    float3  vert2_0 = mul_2(_S744, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + *mean_1;
    float3  vertc_0 = (vert0_0 + vert1_0 + vert2_0) / make_float3 (3.0f);
    float3  d0_0 = vert0_0 - vertc_0;
    float3  d1_0 = vert1_0 - vertc_0;
    float3  d2_0 = vert2_0 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_1(d0_0)), (length_1(d1_0))))), (length_1(d2_0))))) * normalize_1(cross_0(d0_0, d1_0));
    float4  _S745 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_1), dot_0(scale_1, scale_1), dot_1(quat_4, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S746 = _S745 - floor_0(_S745);
    float4  _S747 = _S746 + make_float4 (dot_1(_S746, float4 {_S746.w, _S746.z, _S746.x, _S746.y} + make_float4 (33.3300018310546875f)));
    float4  _S748 = (float4 {_S747.x, _S747.x, _S747.y, _S747.z} + float4 {_S747.y, _S747.z, _S747.z, _S747.w}) * float4 {_S747.z, _S747.y, _S747.w, _S747.x};
    float4  _S749 = _S748 - floor_0(_S748);
    float2  _S750 = float2 {_S749.x, _S749.z};
    float _S751 = 6.28318548202514648f * _S750.y;
    float2  _S752 = float2 {_S749.y, _S749.w};
    float _S753 = 6.28318548202514648f * _S752.y;
    *mean_1 = *mean_1 + mul_2(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S750.x))))))) * make_float2 ((F32_cos((_S751))), (F32_sin((_S751))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S750.x))))))) * make_float2 ((F32_cos((_S751))), (F32_sin((_S751))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S752.x))))))) * make_float2 ((F32_cos((_S753))), (F32_sin((_S753))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_15, float dOut_11)
{
    float _S754 = _slang_select(((*dpx_15).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_15).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_11;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S754;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  dOut_12)
{
    float3  _S755 = _slang_select(((*dpx_16).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_16).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_12;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S755;
    return;
}

inline __device__ float3  abs_0(float3  x_19)
{
    float3  result_19;
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
        *_slang_vector_get_element_ptr(&result_19, i_8) = (F32_abs((_slang_vector_get_element(x_19, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_19;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_17, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_13)
{
    DiffPair_float_0 _S756 = *dpx_17;
    bool _S757;
    if(((*dpx_17).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S757 = ((*dpx_17).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S757 = false;
    }
    float _S758;
    if(_S757)
    {
        _S758 = dOut_13;
    }
    else
    {
        _S758 = 0.0f;
    }
    dpx_17->primal_0 = _S756.primal_0;
    dpx_17->differential_0 = _S758;
    DiffPair_float_0 _S759 = *dpMin_0;
    if((_S756.primal_0) < ((*dpMin_0).primal_0))
    {
        _S758 = dOut_13;
    }
    else
    {
        _S758 = 0.0f;
    }
    dpMin_0->primal_0 = _S759.primal_0;
    dpMin_0->differential_0 = _S758;
    DiffPair_float_0 _S760 = *dpMax_0;
    if(((*dpx_17).primal_0) > ((*dpMax_0).primal_0))
    {
        _S758 = dOut_13;
    }
    else
    {
        _S758 = 0.0f;
    }
    dpMax_0->primal_0 = _S760.primal_0;
    dpMax_0->differential_0 = _S758;
    return;
}

inline __device__ float clamp_0(float x_20, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_20), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_18, float dOut_14)
{
    float _S761 = -0.5f / ((*dpx_18).primal_0 * (F32_sqrt(((*dpx_18).primal_0)))) * dOut_14;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S761;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_19, DiffPair_float_0 * dpy_5, DiffPair_float_0 * dps_0, float dOut_15)
{
    float _S762 = (1.0f - (*dps_0).primal_0) * dOut_15;
    dpx_19->primal_0 = (*dpx_19).primal_0;
    dpx_19->differential_0 = _S762;
    DiffPair_float_0 _S763 = *dpy_5;
    float _S764 = (*dps_0).primal_0 * dOut_15;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = _S764;
    float _S765 = (_S763.primal_0 - (*dpx_19).primal_0) * dOut_15;
    dps_0->primal_0 = _S763.primal_0;
    dps_0->differential_0 = _S765;
    return;
}

inline __device__ float lerp_0(float x_21, float y_8, float s_5)
{
    return x_21 + (y_8 - x_21) * s_5;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 10>  weights_0, FixedArray<float, 23>  * _S766)
{
    float3  _S767;
    bool _S768;
    bool _S769;
    FixedArray<float, 23>  losses_2;
    float _S770 = float(mask_0);
    float3  _S771 = ref_rgb_0 - render_rgb_0;
    float3  _S772 = abs_0(_S771);
    losses_2[int(0)] = weights_0[int(0)] * _S770 * ((_S772.x + _S772.y + _S772.z) * 0.3333333432674408f);
    losses_2[int(1)] = _S770 * clamp_0(dot_0(_S771, _S771) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S773 = float(depth_mask_0 & mask_0);
    float _S774 = _S773 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S775 = _S773 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_2[int(2)] = _S774;
    losses_2[int(3)] = _S775;
    losses_2[int(4)] = _S774 * _S774;
    losses_2[int(5)] = _S775 * _S775;
    losses_2[int(6)] = _S774 * _S775;
    bool _S776 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S777 = norm2_0 == 0.0f;
        _S768 = _S777;
        if(_S777)
        {
            _S767 = make_float3 (0.0f);
            break;
        }
        _S767 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S778;
    bool _S779 = !_S768;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S780 = norm2_1 == 0.0f;
        _S769 = _S780;
        if(_S780)
        {
            _S778 = make_float3 (0.0f);
            break;
        }
        _S778 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S781;
    float3  _S782;
    bool _S783 = !_S769;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S782 = make_float3 (0.0f);
            _S781 = false;
            break;
        }
        _S782 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S781 = _S776;
        break;
    }
    float _S784 = float(_S779 & _S781);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S767, _S782);
    losses_2[int(7)] = weights_0[int(2)] * _S784 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S785 = float(_S783 & _S781);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S778, _S782);
    losses_2[int(8)] = weights_0[int(2)] * _S785 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S786 = float(_S779 & _S783);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S767, _S778);
    losses_2[int(11)] = weights_0[int(5)] * _S786 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S787 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S788 = float(alpha_mask_0);
    float _S789 = float(ref_alpha_0);
    float _S790 = (F32_max((_S787), (_S789)));
    losses_2[int(9)] = weights_0[int(3)] * _S788 * - lerp_0((F32_log(((F32_max((1.0f - _S790), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S790), (9.99999997475242708e-07f)))))), _S789);
    float _S791 = 1.0f - _S787;
    float _S792 = 1.0f - _S789;
    float _S793 = (F32_max((_S791), (_S792)));
    losses_2[int(10)] = weights_0[int(4)] * _S788 * - lerp_0((F32_log(((F32_max((1.0f - _S793), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S793), (9.99999997475242708e-07f)))))), _S792);
    losses_2[int(12)] = weights_0[int(6)] * 4.0f * _S787 * _S791;
    float _S794 = (F32_max((_S787), (9.999999960041972e-13f)));
    losses_2[int(13)] = weights_0[int(7)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S794;
    losses_2[int(14)] = weights_0[int(8)] * depth_dist_0 / _S794;
    losses_2[int(15)] = weights_0[int(9)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S794;
    losses_2[int(16)] = 1.0f;
    losses_2[int(17)] = _S770;
    losses_2[int(18)] = _S773;
    losses_2[int(19)] = _S784;
    losses_2[int(20)] = _S785;
    losses_2[int(21)] = _S786;
    losses_2[int(22)] = _S788;
    *_S766 = losses_2;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S795, float3  _S796)
{
    return dot_0(_S795, _S796);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S797)
{
    return (F32_rsqrt((_S797)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S798, float _S799, float _S800)
{
    return clamp_0(_S798, _S799, _S800);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S801, DiffPair_float_0 * _S802, DiffPair_float_0 * _S803, float _S804)
{
    _d_lerp_0(_S801, _S802, _S803, _S804);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S805, DiffPair_float_0 * _S806, DiffPair_float_0 * _S807, float _S808)
{
    _d_clamp_0(_S805, _S806, _S807, _S808);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S809, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S810, float _S811)
{
    _d_dot_0(_S809, _S810, _S811);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S812, float _S813)
{
    _d_rsqrt_0(_S812, _S813);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S814, float3  _S815)
{
    _d_abs_vector_0(_S814, _S815);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 23>  * _s_dOut_1)
{
    DiffPair_float_0 _S816 = *dprender_depth_0;
    DiffPair_float_0 _S817 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S818 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S819 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S820 = *dpref_normal_0;
    DiffPair_float_0 _S821 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S822 = *dprgb_dist_0;
    DiffPair_float_0 _S823 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S824 = *dpnormal_dist_0;
    float3  _S825 = make_float3 (0.0f);
    float _S826 = float(mask_1);
    float _S827 = (*weights_1)[int(0)] * _S826;
    float3  _S828 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S829 = s_primal_ctx_dot_0(_S828, _S828) * 0.3333333432674408f;
    float _S830 = float(depth_mask_1 & mask_1);
    float _S831 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S832 = _S830 * s_primal_ctx_log_0(_S831);
    float _S833 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S834 = _S830 * s_primal_ctx_log_0(_S833);
    bool _S835 = normal_mask_1 & mask_1;
    float _S836 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S837 = _S836 == 0.0f;
    float3  _S838;
    if(_S837)
    {
        _S838 = make_float3 (0.0f);
    }
    bool _S839 = !_S837;
    float3  _S840;
    if(_S839)
    {
        float _S841 = s_primal_ctx_rsqrt_0(_S836);
        float3  _S842 = make_float3 (_S841);
        _S838 = _S818.primal_0 * make_float3 (_S841);
        _S840 = _S842;
    }
    else
    {
        _S840 = _S825;
    }
    float _S843 = s_primal_ctx_dot_0(_S819.primal_0, _S819.primal_0);
    bool _S844 = _S843 == 0.0f;
    float3  _S845;
    if(_S844)
    {
        _S845 = make_float3 (0.0f);
    }
    bool _S846 = !_S844;
    float3  _S847;
    if(_S846)
    {
        float _S848 = s_primal_ctx_rsqrt_0(_S843);
        float3  _S849 = make_float3 (_S848);
        _S845 = _S819.primal_0 * make_float3 (_S848);
        _S847 = _S849;
    }
    else
    {
        _S847 = _S825;
    }
    float _S850 = s_primal_ctx_dot_0(_S820.primal_0, _S820.primal_0);
    bool _S851 = _S850 == 0.0f;
    float3  _S852;
    bool _S853;
    if(_S851)
    {
        float3  _S854 = make_float3 (0.0f);
        _S853 = false;
        _S852 = _S854;
    }
    else
    {
        _S853 = _S835;
    }
    bool _S855 = !_S851;
    float3  _S856;
    if(_S855)
    {
        float _S857 = s_primal_ctx_rsqrt_0(_S850);
        float3  _S858 = make_float3 (_S857);
        _S852 = _S820.primal_0 * make_float3 (_S857);
        _S856 = _S858;
    }
    else
    {
        _S856 = _S825;
    }
    float _S859 = (*weights_1)[int(2)] * float(_S839 & _S853);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S838, _S852);
    float _S860 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S861 = (*weights_1)[int(2)] * float(_S846 & _S853);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S845, _S852);
    float _S862 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S863 = (*weights_1)[int(5)] * float(_S839 & _S846);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S838, _S845);
    float _S864 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S865 = s_primal_ctx_clamp_0(_S821.primal_0, 0.0f, 1.0f);
    float _S866 = float(alpha_mask_1);
    float _S867 = (*weights_1)[int(3)] * _S866;
    float _S868 = float(ref_alpha_1);
    float _S869 = (F32_max((_S865), (_S868)));
    float _S870 = 1.0f - _S869;
    float _S871 = (F32_max((_S870), (9.99999997475242708e-07f)));
    float _S872 = s_primal_ctx_log_0(_S871);
    float _S873 = (F32_max((_S869), (9.99999997475242708e-07f)));
    float _S874 = s_primal_ctx_log_0(_S873);
    float _S875 = (*weights_1)[int(4)] * _S866;
    float _S876 = 1.0f - _S865;
    float _S877 = 1.0f - _S868;
    float _S878 = (F32_max((_S876), (_S877)));
    float _S879 = 1.0f - _S878;
    float _S880 = (F32_max((_S879), (9.99999997475242708e-07f)));
    float _S881 = s_primal_ctx_log_0(_S880);
    float _S882 = (F32_max((_S878), (9.99999997475242708e-07f)));
    float _S883 = s_primal_ctx_log_0(_S882);
    float _S884 = (*weights_1)[int(6)] * 4.0f;
    float _S885 = _S884 * _S865;
    float _S886 = (F32_max((_S865), (9.999999960041972e-13f)));
    float _S887 = _S886 * _S886;
    float _S888 = (*_s_dOut_1)[int(0)];
    float _S889 = (*_s_dOut_1)[int(1)];
    float _S890 = (*_s_dOut_1)[int(2)];
    float _S891 = (*_s_dOut_1)[int(3)];
    float _S892 = (*_s_dOut_1)[int(4)];
    float _S893 = (*_s_dOut_1)[int(5)];
    float _S894 = (*_s_dOut_1)[int(6)];
    float _S895 = (*_s_dOut_1)[int(15)] / _S887;
    float _S896 = 0.3333333432674408f * ((*weights_1)[int(9)] * (_S886 * _S895));
    float _S897 = (*_s_dOut_1)[int(14)] / _S887;
    float _S898 = (*weights_1)[int(8)] * (_S886 * _S897);
    float _S899 = (*_s_dOut_1)[int(13)] / _S887;
    float _S900 = _S886 * _S899;
    float _S901 = (*weights_1)[int(9)] * ((_S824.primal_0.x + _S824.primal_0.y + _S824.primal_0.z) * 0.3333333432674408f) * - _S895 + (*weights_1)[int(8)] * _S823.primal_0 * - _S897 + (*weights_1)[int(7)] * ((_S822.primal_0.x + _S822.primal_0.y + _S822.primal_0.z) * 0.3333333432674408f) * - _S899;
    DiffPair_float_0 _S902;
    (&_S902)->primal_0 = _S865;
    (&_S902)->differential_0 = 0.0f;
    DiffPair_float_0 _S903;
    (&_S903)->primal_0 = 9.999999960041972e-13f;
    (&_S903)->differential_0 = 0.0f;
    _d_max_0(&_S902, &_S903, _S901);
    float _S904 = 0.3333333432674408f * ((*weights_1)[int(7)] * _S900);
    float _S905 = _S885 * (*_s_dOut_1)[int(12)];
    float _S906 = _S884 * (_S876 * (*_s_dOut_1)[int(12)]);
    float _S907 = - (_S875 * (*_s_dOut_1)[int(10)]);
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S881;
    (&_S908)->differential_0 = 0.0f;
    DiffPair_float_0 _S909;
    (&_S909)->primal_0 = _S883;
    (&_S909)->differential_0 = 0.0f;
    DiffPair_float_0 _S910;
    (&_S910)->primal_0 = _S877;
    (&_S910)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S908, &_S909, &_S910, _S907);
    DiffPair_float_0 _S911;
    (&_S911)->primal_0 = _S882;
    (&_S911)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S911, _S909.differential_0);
    DiffPair_float_0 _S912;
    (&_S912)->primal_0 = _S878;
    (&_S912)->differential_0 = 0.0f;
    DiffPair_float_0 _S913;
    (&_S913)->primal_0 = 9.99999997475242708e-07f;
    (&_S913)->differential_0 = 0.0f;
    _d_max_0(&_S912, &_S913, _S911.differential_0);
    DiffPair_float_0 _S914;
    (&_S914)->primal_0 = _S880;
    (&_S914)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S914, _S908.differential_0);
    DiffPair_float_0 _S915;
    (&_S915)->primal_0 = _S879;
    (&_S915)->differential_0 = 0.0f;
    DiffPair_float_0 _S916;
    (&_S916)->primal_0 = 9.99999997475242708e-07f;
    (&_S916)->differential_0 = 0.0f;
    _d_max_0(&_S915, &_S916, _S914.differential_0);
    float _S917 = _S912.differential_0 + - _S915.differential_0;
    DiffPair_float_0 _S918;
    (&_S918)->primal_0 = _S876;
    (&_S918)->differential_0 = 0.0f;
    DiffPair_float_0 _S919;
    (&_S919)->primal_0 = _S877;
    (&_S919)->differential_0 = 0.0f;
    _d_max_0(&_S918, &_S919, _S917);
    float _S920 = - (_S905 + _S918.differential_0);
    float _S921 = - (_S867 * (*_s_dOut_1)[int(9)]);
    DiffPair_float_0 _S922;
    (&_S922)->primal_0 = _S872;
    (&_S922)->differential_0 = 0.0f;
    DiffPair_float_0 _S923;
    (&_S923)->primal_0 = _S874;
    (&_S923)->differential_0 = 0.0f;
    DiffPair_float_0 _S924;
    (&_S924)->primal_0 = _S868;
    (&_S924)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S922, &_S923, &_S924, _S921);
    DiffPair_float_0 _S925;
    (&_S925)->primal_0 = _S873;
    (&_S925)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S925, _S923.differential_0);
    DiffPair_float_0 _S926;
    (&_S926)->primal_0 = _S869;
    (&_S926)->differential_0 = 0.0f;
    DiffPair_float_0 _S927;
    (&_S927)->primal_0 = 9.99999997475242708e-07f;
    (&_S927)->differential_0 = 0.0f;
    _d_max_0(&_S926, &_S927, _S925.differential_0);
    DiffPair_float_0 _S928;
    (&_S928)->primal_0 = _S871;
    (&_S928)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S928, _S922.differential_0);
    DiffPair_float_0 _S929;
    (&_S929)->primal_0 = _S870;
    (&_S929)->differential_0 = 0.0f;
    DiffPair_float_0 _S930;
    (&_S930)->primal_0 = 9.99999997475242708e-07f;
    (&_S930)->differential_0 = 0.0f;
    _d_max_0(&_S929, &_S930, _S928.differential_0);
    float _S931 = _S926.differential_0 + - _S929.differential_0;
    DiffPair_float_0 _S932;
    (&_S932)->primal_0 = _S865;
    (&_S932)->differential_0 = 0.0f;
    DiffPair_float_0 _S933;
    (&_S933)->primal_0 = _S868;
    (&_S933)->differential_0 = 0.0f;
    _d_max_0(&_S932, &_S933, _S931);
    float _S934 = _S902.differential_0 + _S906 + _S920 + _S932.differential_0;
    DiffPair_float_0 _S935;
    (&_S935)->primal_0 = _S821.primal_0;
    (&_S935)->differential_0 = 0.0f;
    DiffPair_float_0 _S936;
    (&_S936)->primal_0 = 0.0f;
    (&_S936)->differential_0 = 0.0f;
    DiffPair_float_0 _S937;
    (&_S937)->primal_0 = 1.0f;
    (&_S937)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S935, &_S936, &_S937, _S934);
    DiffPair_float_0 _S938 = _S935;
    float _S939 = _S863 * (*_s_dOut_1)[int(11)];
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S864;
    (&_S940)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S940, _S939);
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = cos_sim_loss_5;
    (&_S941)->differential_0 = 0.0f;
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = 9.999999960041972e-13f;
    (&_S942)->differential_0 = 0.0f;
    _d_max_0(&_S941, &_S942, _S940.differential_0);
    float _S943 = 0.5f * - (_S939 + _S941.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S944;
    (&_S944)->primal_0 = _S838;
    (&_S944)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S945;
    (&_S945)->primal_0 = _S845;
    (&_S945)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S944, &_S945, _S943);
    float _S946 = _S861 * (*_s_dOut_1)[int(8)];
    DiffPair_float_0 _S947;
    (&_S947)->primal_0 = _S862;
    (&_S947)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S947, _S946);
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = cos_sim_loss_4;
    (&_S948)->differential_0 = 0.0f;
    DiffPair_float_0 _S949;
    (&_S949)->primal_0 = 9.999999960041972e-13f;
    (&_S949)->differential_0 = 0.0f;
    _d_max_0(&_S948, &_S949, _S947.differential_0);
    float _S950 = 0.5f * - (_S946 + _S948.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S951;
    (&_S951)->primal_0 = _S845;
    (&_S951)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S952;
    (&_S952)->primal_0 = _S852;
    (&_S952)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S951, &_S952, _S950);
    float _S953 = _S859 * (*_s_dOut_1)[int(7)];
    DiffPair_float_0 _S954;
    (&_S954)->primal_0 = _S860;
    (&_S954)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S954, _S953);
    DiffPair_float_0 _S955;
    (&_S955)->primal_0 = cos_sim_loss_3;
    (&_S955)->differential_0 = 0.0f;
    DiffPair_float_0 _S956;
    (&_S956)->primal_0 = 9.999999960041972e-13f;
    (&_S956)->differential_0 = 0.0f;
    _d_max_0(&_S955, &_S956, _S954.differential_0);
    float _S957 = 0.5f * - (_S953 + _S955.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S958;
    (&_S958)->primal_0 = _S838;
    (&_S958)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S959;
    (&_S959)->primal_0 = _S852;
    (&_S959)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S958, &_S959, _S957);
    float3  _S960 = _S952.differential_0 + _S959.differential_0;
    float3  _S961 = _S944.differential_0 + _S958.differential_0;
    float3  _S962 = make_float3 (_S896, _S896, _S896);
    float3  _S963 = make_float3 (_S904, _S904, _S904);
    float3  _S964 = _S945.differential_0 + _S951.differential_0;
    float _S965;
    if(_S855)
    {
        float3  _S966 = _S820.primal_0 * _S960;
        float3  _S967 = _S856 * _S960;
        float _S968 = _S966.x + _S966.y + _S966.z;
        DiffPair_float_0 _S969;
        (&_S969)->primal_0 = _S850;
        (&_S969)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S969, _S968);
        _S965 = _S969.differential_0;
        _S838 = _S967;
    }
    else
    {
        _S965 = 0.0f;
        _S838 = _S825;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S970;
    (&_S970)->primal_0 = _S820.primal_0;
    (&_S970)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S971;
    (&_S971)->primal_0 = _S820.primal_0;
    (&_S971)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S970, &_S971, _S965);
    float3  _S972 = _S971.differential_0 + _S970.differential_0 + _S838;
    if(_S846)
    {
        float3  _S973 = _S819.primal_0 * _S964;
        float3  _S974 = _S847 * _S964;
        float _S975 = _S973.x + _S973.y + _S973.z;
        DiffPair_float_0 _S976;
        (&_S976)->primal_0 = _S843;
        (&_S976)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S976, _S975);
        _S965 = _S976.differential_0;
        _S838 = _S974;
    }
    else
    {
        _S965 = 0.0f;
        _S838 = _S825;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S977;
    (&_S977)->primal_0 = _S819.primal_0;
    (&_S977)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S978;
    (&_S978)->primal_0 = _S819.primal_0;
    (&_S978)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S977, &_S978, _S965);
    float3  _S979 = _S978.differential_0 + _S977.differential_0 + _S838;
    if(_S839)
    {
        float3  _S980 = _S818.primal_0 * _S961;
        float3  _S981 = _S840 * _S961;
        float _S982 = _S980.x + _S980.y + _S980.z;
        DiffPair_float_0 _S983;
        (&_S983)->primal_0 = _S836;
        (&_S983)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S983, _S982);
        _S965 = _S983.differential_0;
        _S838 = _S981;
    }
    else
    {
        _S965 = 0.0f;
        _S838 = _S825;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S984;
    (&_S984)->primal_0 = _S818.primal_0;
    (&_S984)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S985;
    (&_S985)->primal_0 = _S818.primal_0;
    (&_S985)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S984, &_S985, _S965);
    float _S986 = _S834 * _S894;
    float _S987 = _S834 * _S893;
    float _S988 = _S832 * _S892;
    float _S989 = _S830 * (_S832 * _S894 + _S987 + _S987 + _S891);
    DiffPair_float_0 _S990;
    (&_S990)->primal_0 = _S833;
    (&_S990)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S990, _S989);
    DiffPair_float_0 _S991;
    (&_S991)->primal_0 = _S817.primal_0;
    (&_S991)->differential_0 = 0.0f;
    DiffPair_float_0 _S992;
    (&_S992)->primal_0 = 0.00009999999747379f;
    (&_S992)->differential_0 = 0.0f;
    _d_max_0(&_S991, &_S992, _S990.differential_0);
    float _S993 = _S830 * (_S986 + _S988 + _S988 + _S890);
    DiffPair_float_0 _S994;
    (&_S994)->primal_0 = _S831;
    (&_S994)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S994, _S993);
    DiffPair_float_0 _S995;
    (&_S995)->primal_0 = _S816.primal_0;
    (&_S995)->differential_0 = 0.0f;
    DiffPair_float_0 _S996;
    (&_S996)->primal_0 = 0.00009999999747379f;
    (&_S996)->differential_0 = 0.0f;
    _d_max_0(&_S995, &_S996, _S994.differential_0);
    float _S997 = _S826 * _S889;
    DiffPair_float_0 _S998;
    (&_S998)->primal_0 = _S829;
    (&_S998)->differential_0 = 0.0f;
    DiffPair_float_0 _S999;
    (&_S999)->primal_0 = 0.0f;
    (&_S999)->differential_0 = 0.0f;
    DiffPair_float_0 _S1000;
    (&_S1000)->primal_0 = 1.0f;
    (&_S1000)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S998, &_S999, &_S1000, _S997);
    float _S1001 = 0.3333333432674408f * _S998.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1002;
    (&_S1002)->primal_0 = _S828;
    (&_S1002)->differential_0 = _S825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1003;
    (&_S1003)->primal_0 = _S828;
    (&_S1003)->differential_0 = _S825;
    s_bwd_prop_dot_0(&_S1002, &_S1003, _S1001);
    float _S1004 = 0.3333333432674408f * (_S827 * _S888);
    float3  _S1005 = make_float3 (_S1004, _S1004, _S1004);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1006;
    (&_S1006)->primal_0 = _S828;
    (&_S1006)->differential_0 = _S825;
    s_bwd_prop_abs_0(&_S1006, _S1005);
    float3  _S1007 = _S1003.differential_0 + _S1002.differential_0 + _S1006.differential_0;
    float3  _S1008 = - _S1007;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S962;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S898;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S963;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S938.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S972;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S979;
    float3  _S1009 = _S985.differential_0 + _S984.differential_0 + _S838;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S1009;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S991.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S995.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S1007;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S1008;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1010, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1011, DiffPair_float_0 * _S1012, DiffPair_float_0 * _S1013, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1014, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1015, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1016, DiffPair_float_0 * _S1017, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1018, DiffPair_float_0 * _S1019, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1020, bool _S1021, bool _S1022, bool _S1023, bool _S1024, bool _S1025, FixedArray<float, 10>  * _S1026, FixedArray<float, 23>  * _S1027)
{
    s_bwd_prop_per_pixel_losses_0(_S1010, _S1011, _S1012, _S1013, _S1014, _S1015, _S1016, _S1017, _S1018, _S1019, _S1020, _S1021, _S1022, _S1023, _S1024, _S1025, _S1026, _S1027);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 10>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S1028 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S1028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S1028;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S1028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S1028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S1028;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S1028;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S1028;
    FixedArray<float, 10>  _S1029 = weights_2;
    FixedArray<float, 23>  _S1030 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S1029, &_S1030);
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_20, float dOut_16)
{
    float _S1031 = 1.0f / ((*dpx_20).primal_0 * 2.30258512496948242f) * dOut_16;
    dpx_20->primal_0 = (*dpx_20).primal_0;
    dpx_20->differential_0 = _S1031;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 10>  weights_3, FixedArray<float, 10>  * _S1032)
{
    FixedArray<float, 10>  losses_3;
    float _S1033 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_3[int(0)] = raw_losses_0[int(0)] / _S1033;
    losses_3[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S1033)));
    bool _S1034;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S1034 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S1034 = false;
    }
    float _S1035;
    if(_S1034)
    {
        _S1035 = weights_3[int(1)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S1035 = 0.0f;
    }
    losses_3[int(2)] = _S1035;
    losses_3[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_3[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_3[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S1036 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_3[int(6)] = raw_losses_0[int(12)] / _S1036;
    losses_3[int(7)] = raw_losses_0[int(13)] / _S1036;
    losses_3[int(8)] = raw_losses_0[int(14)] / _S1036;
    losses_3[int(9)] = raw_losses_0[int(15)] / _S1036;
    *_S1032 = losses_3;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S1037)
{
    return (F32_sqrt((_S1037)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S1038, float _S1039)
{
    _d_log10_0(_S1038, _S1039);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_4, FixedArray<float, 10>  * _s_dOut_2)
{
    FixedArray<float, 23>  _S1040 = dpraw_losses_0->primal_0;
    float _S1041 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S1042 = _S1041 * _S1041;
    float _S1043 = dpraw_losses_0->primal_0[int(1)] / _S1041;
    bool _S1044 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S1045;
    if(_S1044)
    {
        _S1045 = (_S1040[int(3)]) != 0.0f;
    }
    else
    {
        _S1045 = false;
    }
    float _S1046;
    float _S1047;
    float _S1048;
    float _S1049;
    float _S1050;
    float _S1051;
    float _S1052;
    float _S1053;
    float _S1054;
    float _S1055;
    float _S1056;
    float _S1057;
    float _S1058;
    float _S1059;
    float _S1060;
    if(_S1045)
    {
        float _S1061 = _S1040[int(2)] * _S1040[int(3)];
        float _S1062 = _S1040[int(18)] * _S1040[int(18)];
        float _S1063 = _S1040[int(6)] - _S1061 / _S1040[int(18)];
        float _S1064 = _S1040[int(2)] * _S1040[int(2)];
        float _S1065 = _S1040[int(4)] - _S1064 / _S1040[int(18)];
        float _S1066 = _S1040[int(3)] * _S1040[int(3)];
        float _S1067 = _S1040[int(5)] - _S1066 / _S1040[int(18)];
        float _S1068 = _S1065 * _S1067 + 1.0f;
        float _S1069 = (F32_max((9.999999960041972e-13f), (_S1068)));
        float _S1070 = s_primal_ctx_sqrt_0(_S1069);
        float _S1071 = _S1070 * _S1070;
        float _S1072 = 1.0f - _S1063 / _S1070;
        _S1046 = (*weights_4)[int(1)];
        _S1047 = _S1072;
        _S1048 = _S1071;
        _S1049 = _S1063;
        _S1050 = _S1070;
        _S1051 = _S1069;
        _S1052 = _S1068;
        _S1053 = _S1065;
        _S1054 = _S1067;
        _S1055 = _S1062;
        _S1056 = _S1066;
        _S1057 = _S1040[int(3)];
        _S1058 = _S1064;
        _S1059 = _S1040[int(2)];
        _S1060 = _S1061;
    }
    else
    {
        _S1046 = 0.0f;
        _S1047 = 0.0f;
        _S1048 = 0.0f;
        _S1049 = 0.0f;
        _S1050 = 0.0f;
        _S1051 = 0.0f;
        _S1052 = 0.0f;
        _S1053 = 0.0f;
        _S1054 = 0.0f;
        _S1055 = 0.0f;
        _S1056 = 0.0f;
        _S1057 = 0.0f;
        _S1058 = 0.0f;
        _S1059 = 0.0f;
        _S1060 = 0.0f;
    }
    float _S1073 = (F32_max((_S1040[int(19)]), (1.0f)));
    float _S1074 = _S1073 * _S1073;
    float _S1075 = (F32_max((_S1040[int(20)]), (1.0f)));
    float _S1076 = _S1075 * _S1075;
    float _S1077 = float((I32_max((int((_S1040[int(19)]) > 0.5f) + int((_S1040[int(20)]) > 0.5f)), (int(1)))));
    float _S1078 = _S1040[int(9)] + _S1040[int(10)];
    float _S1079 = (F32_max((_S1040[int(22)]), (1.0f)));
    float _S1080 = _S1079 * _S1079;
    float _S1081 = (F32_max((_S1040[int(21)]), (1.0f)));
    float _S1082 = _S1081 * _S1081;
    float _S1083 = (F32_max((_S1040[int(16)]), (1.0f)));
    float _S1084 = _S1083 * _S1083;
    float _S1085 = (*_s_dOut_2)[int(0)];
    float _S1086 = (*_s_dOut_2)[int(1)];
    float _S1087 = (*_s_dOut_2)[int(2)];
    float _S1088 = (*_s_dOut_2)[int(9)] / _S1084;
    float _S1089 = _S1083 * _S1088;
    float _S1090 = (*_s_dOut_2)[int(8)] / _S1084;
    float _S1091 = _S1083 * _S1090;
    float _S1092 = (*_s_dOut_2)[int(7)] / _S1084;
    float _S1093 = _S1083 * _S1092;
    float _S1094 = (*_s_dOut_2)[int(6)] / _S1084;
    float _S1095 = _S1083 * _S1094;
    float _S1096 = _S1040[int(15)] * - _S1088 + _S1040[int(14)] * - _S1090 + _S1040[int(13)] * - _S1092 + _S1040[int(12)] * - _S1094;
    DiffPair_float_0 _S1097;
    (&_S1097)->primal_0 = _S1040[int(16)];
    (&_S1097)->differential_0 = 0.0f;
    DiffPair_float_0 _S1098;
    (&_S1098)->primal_0 = 1.0f;
    (&_S1098)->differential_0 = 0.0f;
    _d_max_0(&_S1097, &_S1098, _S1096);
    float _S1099 = (*_s_dOut_2)[int(5)] / _S1082;
    float _S1100 = _S1040[int(11)] * - _S1099;
    float _S1101 = _S1081 * _S1099;
    DiffPair_float_0 _S1102;
    (&_S1102)->primal_0 = _S1040[int(21)];
    (&_S1102)->differential_0 = 0.0f;
    DiffPair_float_0 _S1103;
    (&_S1103)->primal_0 = 1.0f;
    (&_S1103)->differential_0 = 0.0f;
    _d_max_0(&_S1102, &_S1103, _S1100);
    float _S1104 = (*_s_dOut_2)[int(4)] / _S1080;
    float _S1105 = _S1078 * - _S1104;
    float _S1106 = _S1079 * _S1104;
    DiffPair_float_0 _S1107;
    (&_S1107)->primal_0 = _S1040[int(22)];
    (&_S1107)->differential_0 = 0.0f;
    DiffPair_float_0 _S1108;
    (&_S1108)->primal_0 = 1.0f;
    (&_S1108)->differential_0 = 0.0f;
    _d_max_0(&_S1107, &_S1108, _S1105);
    float _S1109 = (*_s_dOut_2)[int(3)] / _S1077;
    float _S1110 = _S1109 / _S1076;
    float _S1111 = _S1040[int(8)] * - _S1110;
    float _S1112 = _S1075 * _S1110;
    DiffPair_float_0 _S1113;
    (&_S1113)->primal_0 = _S1040[int(20)];
    (&_S1113)->differential_0 = 0.0f;
    DiffPair_float_0 _S1114;
    (&_S1114)->primal_0 = 1.0f;
    (&_S1114)->differential_0 = 0.0f;
    _d_max_0(&_S1113, &_S1114, _S1111);
    float _S1115 = _S1109 / _S1074;
    float _S1116 = _S1040[int(7)] * - _S1115;
    float _S1117 = _S1073 * _S1115;
    DiffPair_float_0 _S1118;
    (&_S1118)->primal_0 = _S1040[int(19)];
    (&_S1118)->differential_0 = 0.0f;
    DiffPair_float_0 _S1119;
    (&_S1119)->primal_0 = 1.0f;
    (&_S1119)->differential_0 = 0.0f;
    _d_max_0(&_S1118, &_S1119, _S1116);
    FixedArray<float, 23>  _S1120;
    _S1120[int(0)] = 0.0f;
    _S1120[int(1)] = 0.0f;
    _S1120[int(2)] = 0.0f;
    _S1120[int(3)] = 0.0f;
    _S1120[int(4)] = 0.0f;
    _S1120[int(5)] = 0.0f;
    _S1120[int(6)] = 0.0f;
    _S1120[int(7)] = 0.0f;
    _S1120[int(8)] = 0.0f;
    _S1120[int(9)] = 0.0f;
    _S1120[int(10)] = 0.0f;
    _S1120[int(11)] = 0.0f;
    _S1120[int(12)] = 0.0f;
    _S1120[int(13)] = 0.0f;
    _S1120[int(14)] = 0.0f;
    _S1120[int(15)] = 0.0f;
    _S1120[int(16)] = 0.0f;
    _S1120[int(17)] = 0.0f;
    _S1120[int(18)] = 0.0f;
    _S1120[int(19)] = 0.0f;
    _S1120[int(20)] = 0.0f;
    _S1120[int(21)] = 0.0f;
    _S1120[int(22)] = 0.0f;
    _S1120[int(15)] = _S1089;
    _S1120[int(14)] = _S1091;
    _S1120[int(13)] = _S1093;
    _S1120[int(16)] = _S1097.differential_0;
    _S1120[int(12)] = _S1095;
    _S1120[int(21)] = _S1102.differential_0;
    _S1120[int(11)] = _S1101;
    _S1120[int(22)] = _S1107.differential_0;
    _S1120[int(10)] = _S1106;
    _S1120[int(9)] = _S1106;
    _S1120[int(20)] = _S1113.differential_0;
    _S1120[int(8)] = _S1112;
    _S1120[int(19)] = _S1118.differential_0;
    _S1120[int(7)] = _S1117;
    float _S1121 = _S1120[int(0)];
    float _S1122 = _S1120[int(1)];
    float _S1123 = _S1120[int(2)];
    float _S1124 = _S1120[int(3)];
    float _S1125 = _S1120[int(4)];
    float _S1126 = _S1120[int(5)];
    float _S1127 = _S1120[int(6)];
    float _S1128 = _S1120[int(7)];
    float _S1129 = _S1120[int(8)];
    float _S1130 = _S1120[int(9)];
    float _S1131 = _S1120[int(10)];
    float _S1132 = _S1120[int(11)];
    float _S1133 = _S1120[int(12)];
    float _S1134 = _S1120[int(13)];
    float _S1135 = _S1120[int(14)];
    float _S1136 = _S1120[int(15)];
    float _S1137 = _S1120[int(16)];
    float _S1138 = _S1120[int(17)];
    float _S1139 = _S1120[int(18)];
    float _S1140 = _S1120[int(19)];
    float _S1141 = _S1120[int(20)];
    float _S1142 = _S1120[int(21)];
    float _S1143 = _S1120[int(22)];
    FixedArray<float, 23>  _S1144;
    if(_S1045)
    {
        float _S1145 = _S1046 * _S1087;
        DiffPair_float_0 _S1146;
        (&_S1146)->primal_0 = _S1047;
        (&_S1146)->differential_0 = 0.0f;
        DiffPair_float_0 _S1147;
        (&_S1147)->primal_0 = 0.0f;
        (&_S1147)->differential_0 = 0.0f;
        DiffPair_float_0 _S1148;
        (&_S1148)->primal_0 = 2.0f;
        (&_S1148)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S1146, &_S1147, &_S1148, _S1145);
        float _S1149 = - _S1146.differential_0 / _S1048;
        float _S1150 = _S1049 * - _S1149;
        float _S1151 = _S1050 * _S1149;
        DiffPair_float_0 _S1152;
        (&_S1152)->primal_0 = _S1051;
        (&_S1152)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S1152, _S1150);
        DiffPair_float_0 _S1153;
        (&_S1153)->primal_0 = 9.999999960041972e-13f;
        (&_S1153)->differential_0 = 0.0f;
        DiffPair_float_0 _S1154;
        (&_S1154)->primal_0 = _S1052;
        (&_S1154)->differential_0 = 0.0f;
        _d_max_0(&_S1153, &_S1154, _S1152.differential_0);
        float _S1155 = _S1053 * _S1154.differential_0;
        float _S1156 = _S1054 * _S1154.differential_0;
        float _S1157 = - _S1155 / _S1055;
        float _S1158 = _S1057 * (_S1040[int(18)] * _S1157);
        float _S1159 = - _S1156 / _S1055;
        float _S1160 = _S1059 * (_S1040[int(18)] * _S1159);
        float _S1161 = - _S1151 / _S1055;
        float _S1162 = _S1040[int(18)] * _S1161;
        float _S1163 = _S1158 + _S1158 + _S1059 * _S1162;
        float _S1164 = _S1160 + _S1160 + _S1057 * _S1162;
        float _S1165 = _S1056 * - _S1157 + _S1058 * - _S1159 + _S1060 * - _S1161;
        FixedArray<float, 23>  _S1166;
        _S1166[int(0)] = 0.0f;
        _S1166[int(1)] = 0.0f;
        _S1166[int(2)] = 0.0f;
        _S1166[int(3)] = 0.0f;
        _S1166[int(4)] = 0.0f;
        _S1166[int(5)] = 0.0f;
        _S1166[int(6)] = 0.0f;
        _S1166[int(7)] = 0.0f;
        _S1166[int(8)] = 0.0f;
        _S1166[int(9)] = 0.0f;
        _S1166[int(10)] = 0.0f;
        _S1166[int(11)] = 0.0f;
        _S1166[int(12)] = 0.0f;
        _S1166[int(13)] = 0.0f;
        _S1166[int(14)] = 0.0f;
        _S1166[int(15)] = 0.0f;
        _S1166[int(16)] = 0.0f;
        _S1166[int(17)] = 0.0f;
        _S1166[int(18)] = 0.0f;
        _S1166[int(19)] = 0.0f;
        _S1166[int(20)] = 0.0f;
        _S1166[int(21)] = 0.0f;
        _S1166[int(22)] = 0.0f;
        _S1166[int(5)] = _S1155;
        _S1166[int(4)] = _S1156;
        _S1166[int(3)] = _S1163;
        _S1166[int(2)] = _S1164;
        _S1166[int(6)] = _S1151;
        float _S1167 = _S1122 + _S1166[int(1)];
        float _S1168 = _S1123 + _S1166[int(2)];
        float _S1169 = _S1124 + _S1166[int(3)];
        float _S1170 = _S1125 + _S1166[int(4)];
        float _S1171 = _S1126 + _S1166[int(5)];
        float _S1172 = _S1127 + _S1166[int(6)];
        float _S1173 = _S1128 + _S1166[int(7)];
        float _S1174 = _S1129 + _S1166[int(8)];
        float _S1175 = _S1130 + _S1166[int(9)];
        float _S1176 = _S1131 + _S1166[int(10)];
        float _S1177 = _S1132 + _S1166[int(11)];
        float _S1178 = _S1133 + _S1166[int(12)];
        float _S1179 = _S1134 + _S1166[int(13)];
        float _S1180 = _S1135 + _S1166[int(14)];
        float _S1181 = _S1136 + _S1166[int(15)];
        float _S1182 = _S1137 + _S1166[int(16)];
        float _S1183 = _S1138 + _S1166[int(17)];
        float _S1184 = _S1139 + _S1166[int(18)];
        float _S1185 = _S1140 + _S1166[int(19)];
        float _S1186 = _S1141 + _S1166[int(20)];
        float _S1187 = _S1142 + _S1166[int(21)];
        float _S1188 = _S1143 + _S1166[int(22)];
        _S1144[int(0)] = _S1121 + _S1166[int(0)];
        _S1144[int(1)] = _S1167;
        _S1144[int(2)] = _S1168;
        _S1144[int(3)] = _S1169;
        _S1144[int(4)] = _S1170;
        _S1144[int(5)] = _S1171;
        _S1144[int(6)] = _S1172;
        _S1144[int(7)] = _S1173;
        _S1144[int(8)] = _S1174;
        _S1144[int(9)] = _S1175;
        _S1144[int(10)] = _S1176;
        _S1144[int(11)] = _S1177;
        _S1144[int(12)] = _S1178;
        _S1144[int(13)] = _S1179;
        _S1144[int(14)] = _S1180;
        _S1144[int(15)] = _S1181;
        _S1144[int(16)] = _S1182;
        _S1144[int(17)] = _S1183;
        _S1144[int(18)] = _S1184;
        _S1144[int(19)] = _S1185;
        _S1144[int(20)] = _S1186;
        _S1144[int(21)] = _S1187;
        _S1144[int(22)] = _S1188;
        _S1046 = _S1165;
    }
    else
    {
        _S1144[int(0)] = _S1121;
        _S1144[int(1)] = _S1122;
        _S1144[int(2)] = _S1123;
        _S1144[int(3)] = _S1124;
        _S1144[int(4)] = _S1125;
        _S1144[int(5)] = _S1126;
        _S1144[int(6)] = _S1127;
        _S1144[int(7)] = _S1128;
        _S1144[int(8)] = _S1129;
        _S1144[int(9)] = _S1130;
        _S1144[int(10)] = _S1131;
        _S1144[int(11)] = _S1132;
        _S1144[int(12)] = _S1133;
        _S1144[int(13)] = _S1134;
        _S1144[int(14)] = _S1135;
        _S1144[int(15)] = _S1136;
        _S1144[int(16)] = _S1137;
        _S1144[int(17)] = _S1138;
        _S1144[int(18)] = _S1139;
        _S1144[int(19)] = _S1140;
        _S1144[int(20)] = _S1141;
        _S1144[int(21)] = _S1142;
        _S1144[int(22)] = _S1143;
        _S1046 = 0.0f;
    }
    if(_S1044)
    {
        FixedArray<float, 23>  _S1189;
        _S1189[int(0)] = 0.0f;
        _S1189[int(1)] = 0.0f;
        _S1189[int(2)] = 0.0f;
        _S1189[int(3)] = 0.0f;
        _S1189[int(4)] = 0.0f;
        _S1189[int(5)] = 0.0f;
        _S1189[int(6)] = 0.0f;
        _S1189[int(7)] = 0.0f;
        _S1189[int(8)] = 0.0f;
        _S1189[int(9)] = 0.0f;
        _S1189[int(10)] = 0.0f;
        _S1189[int(11)] = 0.0f;
        _S1189[int(12)] = 0.0f;
        _S1189[int(13)] = 0.0f;
        _S1189[int(14)] = 0.0f;
        _S1189[int(15)] = 0.0f;
        _S1189[int(16)] = 0.0f;
        _S1189[int(17)] = 0.0f;
        _S1189[int(18)] = 0.0f;
        _S1189[int(19)] = 0.0f;
        _S1189[int(20)] = 0.0f;
        _S1189[int(21)] = 0.0f;
        _S1189[int(22)] = 0.0f;
        _S1189[int(3)] = 0.0f;
        float _S1190 = _S1144[int(1)] + _S1189[int(1)];
        float _S1191 = _S1144[int(2)] + _S1189[int(2)];
        float _S1192 = _S1144[int(3)] + _S1189[int(3)];
        float _S1193 = _S1144[int(4)] + _S1189[int(4)];
        float _S1194 = _S1144[int(5)] + _S1189[int(5)];
        float _S1195 = _S1144[int(6)] + _S1189[int(6)];
        float _S1196 = _S1144[int(7)] + _S1189[int(7)];
        float _S1197 = _S1144[int(8)] + _S1189[int(8)];
        float _S1198 = _S1144[int(9)] + _S1189[int(9)];
        float _S1199 = _S1144[int(10)] + _S1189[int(10)];
        float _S1200 = _S1144[int(11)] + _S1189[int(11)];
        float _S1201 = _S1144[int(12)] + _S1189[int(12)];
        float _S1202 = _S1144[int(13)] + _S1189[int(13)];
        float _S1203 = _S1144[int(14)] + _S1189[int(14)];
        float _S1204 = _S1144[int(15)] + _S1189[int(15)];
        float _S1205 = _S1144[int(16)] + _S1189[int(16)];
        float _S1206 = _S1144[int(17)] + _S1189[int(17)];
        float _S1207 = _S1144[int(18)] + _S1189[int(18)];
        float _S1208 = _S1144[int(19)] + _S1189[int(19)];
        float _S1209 = _S1144[int(20)] + _S1189[int(20)];
        float _S1210 = _S1144[int(21)] + _S1189[int(21)];
        float _S1211 = _S1144[int(22)] + _S1189[int(22)];
        _S1144[int(0)] = _S1144[int(0)] + _S1189[int(0)];
        _S1144[int(1)] = _S1190;
        _S1144[int(2)] = _S1191;
        _S1144[int(3)] = _S1192;
        _S1144[int(4)] = _S1193;
        _S1144[int(5)] = _S1194;
        _S1144[int(6)] = _S1195;
        _S1144[int(7)] = _S1196;
        _S1144[int(8)] = _S1197;
        _S1144[int(9)] = _S1198;
        _S1144[int(10)] = _S1199;
        _S1144[int(11)] = _S1200;
        _S1144[int(12)] = _S1201;
        _S1144[int(13)] = _S1202;
        _S1144[int(14)] = _S1203;
        _S1144[int(15)] = _S1204;
        _S1144[int(16)] = _S1205;
        _S1144[int(17)] = _S1206;
        _S1144[int(18)] = _S1207;
        _S1144[int(19)] = _S1208;
        _S1144[int(20)] = _S1209;
        _S1144[int(21)] = _S1210;
        _S1144[int(22)] = _S1211;
    }
    float _S1212 = -10.0f * _S1086;
    DiffPair_float_0 _S1213;
    (&_S1213)->primal_0 = _S1043;
    (&_S1213)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S1213, _S1212);
    float _S1214 = _S1213.differential_0 / _S1042;
    float _S1215 = _S1041 * _S1214;
    float _S1216 = _S1085 / _S1042;
    float _S1217 = _S1041 * _S1216;
    float _S1218 = _S1040[int(1)] * - _S1214 + _S1040[int(0)] * - _S1216;
    DiffPair_float_0 _S1219;
    (&_S1219)->primal_0 = _S1040[int(17)];
    (&_S1219)->differential_0 = 0.0f;
    DiffPair_float_0 _S1220;
    (&_S1220)->primal_0 = 1.0f;
    (&_S1220)->differential_0 = 0.0f;
    _d_max_0(&_S1219, &_S1220, _S1218);
    FixedArray<float, 23>  _S1221;
    _S1221[int(0)] = 0.0f;
    _S1221[int(1)] = 0.0f;
    _S1221[int(2)] = 0.0f;
    _S1221[int(3)] = 0.0f;
    _S1221[int(4)] = 0.0f;
    _S1221[int(5)] = 0.0f;
    _S1221[int(6)] = 0.0f;
    _S1221[int(7)] = 0.0f;
    _S1221[int(8)] = 0.0f;
    _S1221[int(9)] = 0.0f;
    _S1221[int(10)] = 0.0f;
    _S1221[int(11)] = 0.0f;
    _S1221[int(12)] = 0.0f;
    _S1221[int(13)] = 0.0f;
    _S1221[int(14)] = 0.0f;
    _S1221[int(15)] = 0.0f;
    _S1221[int(16)] = 0.0f;
    _S1221[int(17)] = 0.0f;
    _S1221[int(18)] = 0.0f;
    _S1221[int(19)] = 0.0f;
    _S1221[int(20)] = 0.0f;
    _S1221[int(21)] = 0.0f;
    _S1221[int(22)] = 0.0f;
    _S1221[int(18)] = _S1046;
    _S1221[int(1)] = _S1215;
    _S1221[int(17)] = _S1219.differential_0;
    _S1221[int(0)] = _S1217;
    FixedArray<float, 23>  _S1222 = {
        _S1144[int(0)] + _S1221[int(0)], _S1144[int(1)] + _S1221[int(1)], _S1144[int(2)] + _S1221[int(2)], _S1144[int(3)] + _S1221[int(3)], _S1144[int(4)] + _S1221[int(4)], _S1144[int(5)] + _S1221[int(5)], _S1144[int(6)] + _S1221[int(6)], _S1144[int(7)] + _S1221[int(7)], _S1144[int(8)] + _S1221[int(8)], _S1144[int(9)] + _S1221[int(9)], _S1144[int(10)] + _S1221[int(10)], _S1144[int(11)] + _S1221[int(11)], _S1144[int(12)] + _S1221[int(12)], _S1144[int(13)] + _S1221[int(13)], _S1144[int(14)] + _S1221[int(14)], _S1144[int(15)] + _S1221[int(15)], _S1144[int(16)] + _S1221[int(16)], _S1144[int(17)] + _S1221[int(17)], _S1144[int(18)] + _S1221[int(18)], _S1144[int(19)] + _S1221[int(19)], _S1144[int(20)] + _S1221[int(20)], _S1144[int(21)] + _S1221[int(21)], _S1144[int(22)] + _S1221[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S1222;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S1223, FixedArray<float, 10>  * _S1224, FixedArray<float, 10>  * _S1225)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S1223, _S1224, _S1225);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 10>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S1226)
{
    FixedArray<float, 23>  _S1227 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S1227;
    FixedArray<float, 10>  _S1228 = weights_5;
    FixedArray<float, 10>  _S1229 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S1228, &_S1229);
    *_S1226 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_22, float3  y_9)
{
    float3  result_20;
    int i_9 = int(0);
    for(;;)
    {
        if(i_9 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_20, i_9) = (F32_min((_slang_vector_get_element(x_22, i_9)), (_slang_vector_get_element(y_9, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_20;
}

inline __device__ float3  max_0(float3  x_23, float3  y_10)
{
    float3  result_21;
    int i_10 = int(0);
    for(;;)
    {
        if(i_10 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_21, i_10) = (F32_max((_slang_vector_get_element(x_23, i_10)), (_slang_vector_get_element(y_10, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_21;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_21, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_17)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_21).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_6).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_17.x);
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_3;
    *&((&right_d_result_3)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_21).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_6).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_17.y);
    *&((&left_d_result_3)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_3)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_21).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_6).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_17.z);
    *&((&left_d_result_3)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_3)->z) = right_dp_2.differential_0;
    dpx_21->primal_0 = (*dpx_21).primal_0;
    dpx_21->differential_0 = left_d_result_3;
    dpy_6->primal_0 = (*dpy_6).primal_0;
    dpy_6->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_3;
    return;
}

inline __device__ float3  clamp_1(float3  x_24, float3  minBound_1, float3  maxBound_1)
{
    return min_0(max_0(x_24, minBound_1), maxBound_1);
}

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_1(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1230, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1231, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1232, float3  _S1233)
{
    _d_clamp_vector_0(_S1230, _S1231, _S1232, _S1233);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_3)
{
    float _S1234 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S1235 = make_float3 (_S1234);
    float3  _S1236 = make_float3 (0.0f);
    float3  _S1237 = make_float3 (1.0f);
    float3  _S1238 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239;
    (&_S1239)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S1234) * (*dpbackground_0).primal_0;
    (&_S1239)->differential_0 = _S1238;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1240;
    (&_S1240)->primal_0 = _S1236;
    (&_S1240)->differential_0 = _S1238;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1241;
    (&_S1241)->primal_0 = _S1237;
    (&_S1241)->differential_0 = _S1238;
    s_bwd_prop_clamp_1(&_S1239, &_S1240, &_S1241, _s_dOut_3);
    float3  _S1242 = _S1235 * _S1239.differential_0;
    float3  _S1243 = (*dpbackground_0).primal_0 * _S1239.differential_0;
    float _S1244 = - (_S1243.x + _S1243.y + _S1243.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S1242;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S1244;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S1239.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1245, DiffPair_float_0 * _S1246, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1247, float3  _S1248)
{
    s_bwd_prop_blend_background_0(_S1245, _S1246, _S1247, _S1248);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S1249 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S1249;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S1249;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_22, DiffPair_float_0 * dpy_7, float dOut_18)
{
    if(((*dpx_22).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_22->primal_0 = (*dpx_22).primal_0;
        dpx_22->differential_0 = 0.0f;
        dpy_7->primal_0 = (*dpy_7).primal_0;
        dpy_7->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_22).primal_0), ((*dpy_7).primal_0)));
        DiffPair_float_0 _S1250 = *dpx_22;
        float _S1251 = val_0 * (*dpy_7).primal_0 / (*dpx_22).primal_0 * dOut_18;
        dpx_22->primal_0 = (*dpx_22).primal_0;
        dpx_22->differential_0 = _S1251;
        float _S1252 = val_0 * (F32_log((_S1250.primal_0))) * dOut_18;
        dpy_7->primal_0 = (*dpy_7).primal_0;
        dpy_7->differential_0 = _S1252;
    }
    return;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S1253 = rgb_2;
    float _S1254;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S1254 = _S1253.x * 12.92000007629394531f;
    }
    else
    {
        _S1254 = 1.0549999475479126f * (F32_pow((_S1253.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1253)->x) = _S1254;
    if((_S1253.y) < 0.00313080009073019f)
    {
        _S1254 = _S1253.y * 12.92000007629394531f;
    }
    else
    {
        _S1254 = 1.0549999475479126f * (F32_pow((_S1253.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1253)->y) = _S1254;
    if((_S1253.z) < 0.00313080009073019f)
    {
        _S1254 = _S1253.z * 12.92000007629394531f;
    }
    else
    {
        _S1254 = 1.0549999475479126f * (F32_pow((_S1253.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S1253)->z) = _S1254;
    return _S1253;
}

inline __device__ float s_primal_ctx_pow_0(float _S1255, float _S1256)
{
    return (F32_pow((_S1255), (_S1256)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S1257, DiffPair_float_0 * _S1258, float _S1259)
{
    _d_pow_0(_S1257, _S1258, _S1259);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1260 = *dprgb_1;
    float _S1261 = (*dprgb_1).primal_0.x;
    bool _S1262 = _S1261 < 0.00313080009073019f;
    float _S1263;
    if(_S1262)
    {
        _S1263 = _S1261 * 12.92000007629394531f;
    }
    else
    {
        _S1263 = 1.0549999475479126f * s_primal_ctx_pow_0(_S1261, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S1264 = _S1260.primal_0;
    *&((&_S1264)->x) = _S1263;
    float _S1265 = _S1264.y;
    bool _S1266 = _S1265 < 0.00313080009073019f;
    if(_S1266)
    {
        _S1263 = _S1265 * 12.92000007629394531f;
    }
    else
    {
        _S1263 = 1.0549999475479126f * s_primal_ctx_pow_0(_S1265, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S1264)->y) = _S1263;
    float _S1267 = _S1264.z;
    bool _S1268 = _S1267 < 0.00313080009073019f;
    _S1264 = _s_dOut_4;
    *&((&_S1264)->z) = 0.0f;
    if(_S1268)
    {
        _S1263 = 12.92000007629394531f * _s_dOut_4.z;
    }
    else
    {
        float _S1269 = 1.0549999475479126f * _s_dOut_4.z;
        DiffPair_float_0 _S1270;
        (&_S1270)->primal_0 = _S1267;
        (&_S1270)->differential_0 = 0.0f;
        DiffPair_float_0 _S1271;
        (&_S1271)->primal_0 = 0.4166666567325592f;
        (&_S1271)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1270, &_S1271, _S1269);
        _S1263 = _S1270.differential_0;
    }
    float3  _S1272 = _S1264 + make_float3 (0.0f, 0.0f, _S1263);
    _S1264 = _S1272;
    *&((&_S1264)->y) = 0.0f;
    if(_S1266)
    {
        _S1263 = 12.92000007629394531f * _S1272.y;
    }
    else
    {
        float _S1273 = 1.0549999475479126f * _S1272.y;
        DiffPair_float_0 _S1274;
        (&_S1274)->primal_0 = _S1265;
        (&_S1274)->differential_0 = 0.0f;
        DiffPair_float_0 _S1275;
        (&_S1275)->primal_0 = 0.4166666567325592f;
        (&_S1275)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1274, &_S1275, _S1273);
        _S1263 = _S1274.differential_0;
    }
    float3  _S1276 = _S1264 + make_float3 (0.0f, _S1263, 0.0f);
    _S1264 = _S1276;
    *&((&_S1264)->x) = 0.0f;
    if(_S1262)
    {
        _S1263 = 12.92000007629394531f * _S1276.x;
    }
    else
    {
        float _S1277 = 1.0549999475479126f * _S1276.x;
        DiffPair_float_0 _S1278;
        (&_S1278)->primal_0 = _S1261;
        (&_S1278)->differential_0 = 0.0f;
        DiffPair_float_0 _S1279;
        (&_S1279)->primal_0 = 0.4166666567325592f;
        (&_S1279)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1278, &_S1279, _S1277);
        _S1263 = _S1278.differential_0;
    }
    float3  _S1280 = _S1264 + make_float3 (_S1263, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S1280;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1281, float3  _S1282)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S1281, _S1282);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S1283 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S1283;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool undistort_point_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0, int maxiter_0, float2  * uv_undist_0)
{
    int i_11 = int(0);
    float2  q_0 = uv_0;
    for(;;)
    {
        if(i_11 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S1284 = (*dist_coeffs_0)[int(3)];
        float _S1285 = (*dist_coeffs_0)[int(4)];
        float _S1286 = (*dist_coeffs_0)[int(5)];
        float _S1287 = (*dist_coeffs_0)[int(6)];
        float _S1288 = (*dist_coeffs_0)[int(7)];
        float _S1289 = (*dist_coeffs_0)[int(8)];
        float _S1290 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S1291 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S1292 = (*dist_coeffs_0)[int(1)] + r2_0 * _S1291;
        float _S1293 = (*dist_coeffs_0)[int(0)] + r2_0 * _S1292;
        float radial_0 = 1.0f + r2_0 * _S1293;
        float _S1294 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S1295 = _S1294 * u_0;
        float _S1296 = 2.0f * u_0;
        float _S1297 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S1298 = _S1297 * u_0;
        float _S1299 = 2.0f * v_0;
        float2  _S1300 = q_0 * make_float2 (radial_0) + make_float2 (_S1295 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S1296 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S1298 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S1299 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_4 = _S1300 + make_float2 ((*dist_coeffs_0)[int(8)] * _S1300.x + (*dist_coeffs_0)[int(9)] * _S1300.y, 0.0f) - uv_0;
        float _S1301 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S1301 + _S1301);
        float2  _S1302 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S1293 + (s_diff_r2_0 * _S1292 + (s_diff_r2_0 * _S1291 + s_diff_r2_0 * _S1284 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S1294 * v_0 + 0.0f * _S1295 + (s_diff_r2_0 + (_S1296 + _S1296)) * _S1286 + s_diff_r2_0 * _S1287, _S1297 * v_0 + 0.0f * _S1298 + (s_diff_r2_0 + (_S1301 + 0.0f * _S1299)) * _S1285 + s_diff_r2_0 * _S1288);
        float _S1303 = 0.0f * u_0;
        float s_diff_r2_1 = _S1303 + _S1303 + (v_0 + v_0);
        float2  _S1304 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S1293 + (s_diff_r2_1 * _S1292 + (s_diff_r2_1 * _S1291 + s_diff_r2_1 * _S1284 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S1294 * v_0 + _S1295 + (s_diff_r2_1 + (_S1303 + 0.0f * _S1296)) * _S1286 + s_diff_r2_1 * _S1287, 0.0f * _S1297 * v_0 + _S1298 + (s_diff_r2_1 + (_S1299 + _S1299)) * _S1285 + s_diff_r2_1 * _S1288);
        Matrix<float, 2, 2>  _S1305 = transpose_1(makeMatrix<float, 2, 2> (_S1302 + make_float2 (_S1302.x * _S1289 + _S1302.y * _S1290, 0.0f), _S1304 + make_float2 (_S1304.x * _S1289 + _S1304.y * _S1290, 0.0f)));
        float inv_det_0 = 1.0f / (_S1305.rows[int(0)].x * _S1305.rows[int(1)].y - _S1305.rows[int(0)].y * _S1305.rows[int(1)].x);
        float _S1306 = r_4.x;
        float _S1307 = r_4.y;
        float2  q_1 = q_0 - make_float2 ((_S1306 * _S1305.rows[int(1)].y - _S1307 * _S1305.rows[int(0)].y) * inv_det_0, (- _S1306 * _S1305.rows[int(1)].x + _S1307 * _S1305.rows[int(0)].x) * inv_det_0);
        i_11 = i_11 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S1308 = (*dist_coeffs_0)[int(0)];
    float _S1309 = (*dist_coeffs_0)[int(1)];
    float _S1310 = (*dist_coeffs_0)[int(2)];
    float _S1311 = (*dist_coeffs_0)[int(3)];
    float _S1312 = (*dist_coeffs_0)[int(4)];
    float _S1313 = (*dist_coeffs_0)[int(5)];
    float _S1314 = (*dist_coeffs_0)[int(6)];
    float _S1315 = (*dist_coeffs_0)[int(7)];
    float _S1316 = (*dist_coeffs_0)[int(8)];
    float _S1317 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S1318 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S1318 + _S1318);
    float _S1319 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S1320 = (*dist_coeffs_0)[int(1)] + r2_1 * _S1319;
    float _S1321 = (*dist_coeffs_0)[int(0)] + r2_1 * _S1320;
    float radial_1 = 1.0f + r2_1 * _S1321;
    float _S1322 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S1323 = _S1322 * u_1;
    float _S1324 = 2.0f * u_1;
    float _S1325 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S1326 = _S1325 * u_1;
    float _S1327 = 2.0f * v_1;
    float2  _S1328 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S1321 + (s_diff_r2_2 * _S1320 + (s_diff_r2_2 * _S1319 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S1322 * v_1 + 0.0f * _S1323 + (s_diff_r2_2 + (_S1324 + _S1324)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S1325 * v_1 + 0.0f * _S1326 + (s_diff_r2_2 + (_S1318 + 0.0f * _S1327)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S1329 = 0.0f * u_1;
    float s_diff_r2_3 = _S1329 + _S1329 + (v_1 + v_1);
    float2  _S1330 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S1321 + (s_diff_r2_3 * _S1320 + (s_diff_r2_3 * _S1319 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S1322 * v_1 + _S1323 + (s_diff_r2_3 + (_S1329 + 0.0f * _S1324)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S1325 * v_1 + _S1326 + (s_diff_r2_3 + (_S1327 + _S1327)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S1331 = transpose_1(makeMatrix<float, 2, 2> (_S1328 + make_float2 (_S1328.x * (*dist_coeffs_0)[int(8)] + _S1328.y * (*dist_coeffs_0)[int(9)], 0.0f), _S1330 + make_float2 (_S1330.x * (*dist_coeffs_0)[int(8)] + _S1330.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S1332;
    if((F32_min((determinant_0(_S1331)), ((F32_min((_S1331.rows[int(0)].x), (_S1331.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S1333 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S1308 + r2_2 * (_S1309 + r2_2 * (_S1310 + r2_2 * _S1311)))) + make_float2 (_S1322 * u_2 * v_2 + _S1313 * (r2_2 + 2.0f * u_2 * u_2) + _S1314 * r2_2, _S1325 * u_2 * v_2 + _S1312 * (r2_2 + 2.0f * v_2 * v_2) + _S1315 * r2_2);
        _S1332 = (length_2(_S1333 + make_float2 (_S1316 * _S1333.x + _S1317 * _S1333.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S1332 = false;
    }
    return _S1332;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S1334 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S1334;
    FixedArray<float, 10>  _S1335 = dist_coeffs_1;
    bool _S1336 = undistort_point_0(_S1334, &_S1335, int(12), &uv_1);
    if(!_S1336)
    {
        int3  _S1337 = make_int3 (int(0));
        float3  _S1338 = make_float3 ((float)_S1337.x, (float)_S1337.y, (float)_S1337.z);
        return _S1338;
    }
    float3  raydir_0;
    if(is_fisheye_0)
    {
        float theta_0 = length_2(uv_1);
        float3  raydir_1 = make_float3 ((uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).x, (uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).y, (F32_cos((theta_0))));
        if(!is_ray_depth_0)
        {
            raydir_0 = raydir_1 / make_float3 (raydir_1.z);
        }
        else
        {
            raydir_0 = raydir_1;
        }
    }
    else
    {
        float3  raydir_2 = make_float3 (uv_1.x, uv_1.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_0 = normalize_1(raydir_2);
        }
        else
        {
            raydir_0 = raydir_2;
        }
    }
    return raydir_0;
}

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  points_0)
{
    float3  normal_0 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    float3  normal_1;
    if((dot_0(normal_0, normal_0)) != 0.0f)
    {
        normal_1 = normal_0 / make_float3 (length_1(normal_0));
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1339, float3  _S1340)
{
    return cross_0(_S1339, _S1340);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_23, float _s_dOut_5)
{
    float _S1341 = (*dpx_23).primal_0.x;
    float _S1342 = (*dpx_23).primal_0.y;
    float _S1343 = (*dpx_23).primal_0.z;
    DiffPair_float_0 _S1344;
    (&_S1344)->primal_0 = _S1341 * _S1341 + _S1342 * _S1342 + _S1343 * _S1343;
    (&_S1344)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1344, _s_dOut_5);
    float _S1345 = (*dpx_23).primal_0.z * _S1344.differential_0;
    float _S1346 = _S1345 + _S1345;
    float _S1347 = (*dpx_23).primal_0.y * _S1344.differential_0;
    float _S1348 = _S1347 + _S1347;
    float _S1349 = (*dpx_23).primal_0.x * _S1344.differential_0;
    float _S1350 = _S1349 + _S1349;
    float3  _S1351 = make_float3 (0.0f);
    *&((&_S1351)->z) = _S1346;
    *&((&_S1351)->y) = _S1348;
    *&((&_S1351)->x) = _S1350;
    dpx_23->primal_0 = (*dpx_23).primal_0;
    dpx_23->differential_0 = _S1351;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1352, float _S1353)
{
    s_bwd_prop_length_impl_1(_S1352, _S1353);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1354, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1355, float3  _S1356)
{
    _d_cross_0(_S1354, _S1355, _S1356);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_6)
{
    float3  _S1357 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S1358 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S1359 = s_primal_ctx_cross_0(dx_0, _S1358);
    bool _S1360 = (s_primal_ctx_dot_0(_S1359, _S1359)) != 0.0f;
    float3  _S1361;
    float3  _S1362;
    if(_S1360)
    {
        float _S1363 = length_1(_S1359);
        float3  _S1364 = make_float3 (_S1363);
        _S1361 = make_float3 (_S1363 * _S1363);
        _S1362 = _S1364;
    }
    else
    {
        _S1361 = _S1357;
        _S1362 = _S1357;
    }
    if(_S1360)
    {
        float3  _S1365 = _s_dOut_6 / _S1361;
        float3  _S1366 = _S1359 * - _S1365;
        float3  _S1367 = _S1362 * _S1365;
        float _S1368 = _S1366.x + _S1366.y + _S1366.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1369;
        (&_S1369)->primal_0 = _S1359;
        (&_S1369)->differential_0 = _S1357;
        s_bwd_length_impl_1(&_S1369, _S1368);
        _S1361 = _S1367 + _S1369.differential_0;
    }
    else
    {
        _S1361 = _s_dOut_6;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1370;
    (&_S1370)->primal_0 = _S1359;
    (&_S1370)->differential_0 = _S1357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1371;
    (&_S1371)->primal_0 = _S1359;
    (&_S1371)->differential_0 = _S1357;
    s_bwd_prop_dot_0(&_S1370, &_S1371, 0.0f);
    float3  _S1372 = _S1371.differential_0 + _S1370.differential_0 + _S1361;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1373;
    (&_S1373)->primal_0 = dx_0;
    (&_S1373)->differential_0 = _S1357;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1374;
    (&_S1374)->primal_0 = _S1358;
    (&_S1374)->differential_0 = _S1357;
    s_bwd_prop_cross_0(&_S1373, &_S1374, _S1372);
    float3  s_diff_dy_T_0 = - _S1374.differential_0;
    float3  _S1375 = - s_diff_dy_T_0;
    float3  _S1376 = - _S1373.differential_0;
    FixedArray<float3 , 4>  _S1377;
    _S1377[int(0)] = _S1357;
    _S1377[int(1)] = _S1357;
    _S1377[int(2)] = _S1357;
    _S1377[int(3)] = _S1357;
    _S1377[int(2)] = _S1375;
    _S1377[int(3)] = s_diff_dy_T_0;
    _S1377[int(0)] = _S1376;
    _S1377[int(1)] = _S1373.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S1377;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S1378, float3  _S1379)
{
    s_bwd_prop_points_to_normal_0(_S1378, _S1379);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S1380 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S1380;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S1381 = float2 {intrins_1.z, intrins_1.w};
    float2  _S1382 = float2 {intrins_1.x, intrins_1.y};
    float2  _S1383 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S1381) / _S1382;
    float2  uv_2 = _S1383;
    FixedArray<float, 10>  _S1384 = dist_coeffs_2;
    bool _S1385 = undistort_point_0(_S1383, &_S1384, int(12), &uv_2);
    if(!_S1385)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_3;
    if(is_fisheye_1)
    {
        float theta_1 = length_2(uv_2);
        float3  raydir_4 = make_float3 ((uv_2 / make_float2 ((F32_max((theta_1), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_1))))).x, (uv_2 / make_float2 ((F32_max((theta_1), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_1))))).y, (F32_cos((theta_1))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_4 / make_float3 (raydir_4.z);
        }
        else
        {
            raydir_3 = raydir_4;
        }
    }
    else
    {
        float3  raydir_5 = make_float3 (uv_2.x, uv_2.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_1(raydir_5);
        }
        else
        {
            raydir_3 = raydir_5;
        }
    }
    points_2[int(0)] = make_float3 (depths_0.x) * raydir_3;
    float2  _S1386 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S1381) / _S1382;
    float2  uv_3 = _S1386;
    FixedArray<float, 10>  _S1387 = dist_coeffs_2;
    bool _S1388 = undistort_point_0(_S1386, &_S1387, int(12), &uv_3);
    if(!_S1388)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_2 = length_2(uv_3);
        float3  raydir_6 = make_float3 ((uv_3 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).x, (uv_3 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).y, (F32_cos((theta_2))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_6 / make_float3 (raydir_6.z);
        }
        else
        {
            raydir_3 = raydir_6;
        }
    }
    else
    {
        float3  raydir_7 = make_float3 (uv_3.x, uv_3.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_1(raydir_7);
        }
        else
        {
            raydir_3 = raydir_7;
        }
    }
    points_2[int(1)] = make_float3 (depths_0.y) * raydir_3;
    float2  _S1389 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S1381) / _S1382;
    float2  uv_4 = _S1389;
    FixedArray<float, 10>  _S1390 = dist_coeffs_2;
    bool _S1391 = undistort_point_0(_S1389, &_S1390, int(12), &uv_4);
    if(!_S1391)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_3 = length_2(uv_4);
        float3  raydir_8 = make_float3 ((uv_4 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).x, (uv_4 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).y, (F32_cos((theta_3))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_8 / make_float3 (raydir_8.z);
        }
        else
        {
            raydir_3 = raydir_8;
        }
    }
    else
    {
        float3  raydir_9 = make_float3 (uv_4.x, uv_4.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_1(raydir_9);
        }
        else
        {
            raydir_3 = raydir_9;
        }
    }
    points_2[int(2)] = make_float3 (depths_0.z) * raydir_3;
    float2  _S1392 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S1381) / _S1382;
    float2  uv_5 = _S1392;
    FixedArray<float, 10>  _S1393 = dist_coeffs_2;
    bool _S1394 = undistort_point_0(_S1392, &_S1393, int(12), &uv_5);
    if(!_S1394)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_4 = length_2(uv_5);
        float3  raydir_10 = make_float3 ((uv_5 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_5 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_10 / make_float3 (raydir_10.z);
        }
        else
        {
            raydir_3 = raydir_10;
        }
    }
    else
    {
        float3  raydir_11 = make_float3 (uv_5.x, uv_5.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_1(raydir_11);
        }
        else
        {
            raydir_3 = raydir_11;
        }
    }
    points_2[int(3)] = make_float3 (depths_0.w) * raydir_3;
    float3  normal_2 = cross_0(points_2[int(1)] - points_2[int(0)], - (points_2[int(3)] - points_2[int(2)]));
    float3  normal_3;
    if((dot_0(normal_2, normal_2)) != 0.0f)
    {
        normal_3 = normal_2 / make_float3 (length_1(normal_2));
    }
    else
    {
        normal_3 = normal_2;
    }
    return normal_3;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S1395;
    bool _S1396;
    float2  _S1397;
    bool _S1398;
    float2  _S1399;
    bool _S1400;
    float2  _S1401;
    bool _S1402;
};

inline __device__ float s_primal_ctx_sin_0(float _S1403)
{
    return (F32_sin((_S1403)));
}

inline __device__ float s_primal_ctx_cos_0(float _S1404)
{
    return (F32_cos((_S1404)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S1405 = make_float2 (0.0f);
    _s_diff_ctx_0->_S1395 = _S1405;
    _s_diff_ctx_0->_S1396 = false;
    _s_diff_ctx_0->_S1397 = _S1405;
    _s_diff_ctx_0->_S1398 = false;
    _s_diff_ctx_0->_S1399 = _S1405;
    _s_diff_ctx_0->_S1400 = false;
    _s_diff_ctx_0->_S1401 = _S1405;
    _s_diff_ctx_0->_S1402 = false;
    _s_diff_ctx_0->_S1397 = _S1405;
    _s_diff_ctx_0->_S1398 = false;
    _s_diff_ctx_0->_S1399 = _S1405;
    _s_diff_ctx_0->_S1400 = false;
    _s_diff_ctx_0->_S1401 = _S1405;
    _s_diff_ctx_0->_S1402 = false;
    float3  _S1406 = make_float3 (0.0f);
    float2  _S1407 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1408 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1409 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1407) / _S1408;
    float2  _S1410 = _S1409;
    bool _S1411 = undistort_point_0(_S1409, dist_coeffs_3, int(12), &_S1410);
    _s_diff_ctx_0->_S1395 = _S1410;
    _s_diff_ctx_0->_S1396 = _S1411;
    float2  uv_6 = _S1410;
    bool _S1412 = !_S1411;
    float3  normal_4;
    if(_S1412)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1413 = !_S1412;
    int _S1414;
    FixedArray<float3 , 4>  points_3;
    if(_S1413)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S1415 = length_2(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S1415), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1415))).x, (uv_6 / make_float2 ((F32_max((_S1415), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1415))).y, s_primal_ctx_cos_0(_S1415));
            if(!is_ray_depth_2)
            {
                raydir_12 = raydir_13 / make_float3 (raydir_13.z);
            }
            else
            {
                raydir_12 = raydir_13;
            }
        }
        else
        {
            float3  raydir_14 = make_float3 (uv_6.x, uv_6.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_12 = normalize_1(raydir_14);
            }
            else
            {
                raydir_12 = raydir_14;
            }
        }
        float3  _S1416 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S1417 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1407) / _S1408;
        float2  _S1418 = _S1417;
        bool _S1419 = undistort_point_0(_S1417, dist_coeffs_3, int(12), &_S1418);
        _s_diff_ctx_0->_S1397 = _S1418;
        _s_diff_ctx_0->_S1398 = _S1419;
        float2  uv_7 = _S1418;
        bool _S1420 = !_S1419;
        if(_S1420)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1421 = !_S1420;
        if(_S1421)
        {
            if(is_fisheye_2)
            {
                float _S1422 = length_2(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S1422), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1422))).x, (uv_7 / make_float2 ((F32_max((_S1422), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1422))).y, s_primal_ctx_cos_0(_S1422));
                if(!is_ray_depth_2)
                {
                    raydir_12 = raydir_15 / make_float3 (raydir_15.z);
                }
                else
                {
                    raydir_12 = raydir_15;
                }
            }
            else
            {
                float3  raydir_16 = make_float3 (uv_7.x, uv_7.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_12 = normalize_1(raydir_16);
                }
                else
                {
                    raydir_12 = raydir_16;
                }
            }
            float3  _S1423 = make_float3 (dpdepths_0.y) * raydir_12;
            _S1414 = int(2);
            points_3[int(0)] = _S1416;
            points_3[int(1)] = _S1423;
            points_3[int(2)] = _S1406;
            points_3[int(3)] = _S1406;
        }
        else
        {
            _S1414 = int(0);
            points_3[int(0)] = _S1416;
            points_3[int(1)] = _S1406;
            points_3[int(2)] = _S1406;
            points_3[int(3)] = _S1406;
        }
        bool _runFlag_0;
        if(_S1414 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1413;
            _S1414 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1424 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1407) / _S1408;
            float2  _S1425 = _S1424;
            bool _S1426 = undistort_point_0(_S1424, dist_coeffs_3, int(12), &_S1425);
            _s_diff_ctx_0->_S1399 = _S1425;
            _s_diff_ctx_0->_S1400 = _S1426;
            float2  uv_8 = _S1425;
            if(!_S1426)
            {
                float3  _S1427 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1414 = int(0);
                normal_4 = _S1427;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S1428 = length_2(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S1428), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1428))).x, (uv_8 / make_float2 ((F32_max((_S1428), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1428))).y, s_primal_ctx_cos_0(_S1428));
                    if(!is_ray_depth_2)
                    {
                        raydir_12 = raydir_17 / make_float3 (raydir_17.z);
                    }
                    else
                    {
                        raydir_12 = raydir_17;
                    }
                }
                else
                {
                    float3  raydir_18 = make_float3 (uv_8.x, uv_8.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_12 = normalize_1(raydir_18);
                    }
                    else
                    {
                        raydir_12 = raydir_18;
                    }
                }
                points_3[int(2)] = make_float3 (dpdepths_0.z) * raydir_12;
                float2  _S1429 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1407) / _S1408;
                float2  _S1430 = _S1429;
                bool _S1431 = undistort_point_0(_S1429, dist_coeffs_3, int(12), &_S1430);
                _s_diff_ctx_0->_S1401 = _S1430;
                _s_diff_ctx_0->_S1402 = _S1431;
                float2  uv_9 = _S1430;
                bool _S1432 = !_S1431;
                if(_S1432)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1433 = !_S1432;
                int _S1434;
                if(_S1433)
                {
                    if(is_fisheye_2)
                    {
                        float _S1435 = length_2(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S1435), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1435))).x, (uv_9 / make_float2 ((F32_max((_S1435), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1435))).y, s_primal_ctx_cos_0(_S1435));
                        if(!is_ray_depth_2)
                        {
                            raydir_12 = raydir_19 / make_float3 (raydir_19.z);
                        }
                        else
                        {
                            raydir_12 = raydir_19;
                        }
                    }
                    else
                    {
                        float3  raydir_20 = make_float3 (uv_9.x, uv_9.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_12 = normalize_1(raydir_20);
                        }
                        else
                        {
                            raydir_12 = raydir_20;
                        }
                    }
                    points_3[int(3)] = make_float3 (dpdepths_0.w) * raydir_12;
                    _S1434 = int(2);
                }
                else
                {
                    _S1434 = int(0);
                }
                if(_S1434 != int(2))
                {
                    _runFlag_0 = false;
                    _S1414 = _S1434;
                }
                if(_runFlag_0)
                {
                    _S1414 = int(1);
                }
            }
        }
    }
    else
    {
        _S1414 = int(0);
        points_3[int(0)] = _S1406;
        points_3[int(1)] = _S1406;
        points_3[int(2)] = _S1406;
        points_3[int(3)] = _S1406;
    }
    if(!(_S1414 != int(1)))
    {
        float3  _S1436 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1436, _S1436)) != 0.0f)
        {
            normal_4 = _S1436 / make_float3 (length_1(_S1436));
        }
        else
        {
            normal_4 = _S1436;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_7, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1437 = *dpdepths_1;
    float3  _S1438 = make_float3 (0.0f);
    float2  _S1439 = _s_diff_ctx_1->_S1395;
    bool _S1440 = !!_s_diff_ctx_1->_S1396;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S1441;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1442;
    if(_S1440)
    {
        if(is_fisheye_3)
        {
            float _S1443 = length_2(_S1439);
            float3  raydir_25 = make_float3 ((_S1439 / make_float2 ((F32_max((_S1443), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1443))).x, (_S1439 / make_float2 ((F32_max((_S1443), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1443))).y, s_primal_ctx_cos_0(_S1443));
            if(!is_ray_depth_3)
            {
                raydir_21 = raydir_25 / make_float3 (raydir_25.z);
            }
            else
            {
                raydir_21 = raydir_25;
            }
        }
        else
        {
            float3  raydir_26 = make_float3 (_S1439.x, _S1439.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_1(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S1444 = make_float3 (_S1437.primal_0.x) * raydir_21;
        float2  _S1445 = _s_diff_ctx_1->_S1397;
        bool _S1446 = !!_s_diff_ctx_1->_S1398;
        if(_S1446)
        {
            if(is_fisheye_3)
            {
                float _S1447 = length_2(_S1445);
                float3  raydir_27 = make_float3 ((_S1445 / make_float2 ((F32_max((_S1447), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1447))).x, (_S1445 / make_float2 ((F32_max((_S1447), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1447))).y, s_primal_ctx_cos_0(_S1447));
                if(!is_ray_depth_3)
                {
                    raydir_22 = raydir_27 / make_float3 (raydir_27.z);
                }
                else
                {
                    raydir_22 = raydir_27;
                }
            }
            else
            {
                float3  raydir_28 = make_float3 (_S1445.x, _S1445.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_1(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S1448 = make_float3 (_S1437.primal_0.y) * raydir_22;
            _S1441 = int(2);
            points_4[int(0)] = _S1444;
            points_4[int(1)] = _S1448;
            points_4[int(2)] = _S1438;
            points_4[int(3)] = _S1438;
        }
        else
        {
            _S1441 = int(0);
            points_4[int(0)] = _S1444;
            points_4[int(1)] = _S1438;
            points_4[int(2)] = _S1438;
            points_4[int(3)] = _S1438;
            raydir_22 = _S1438;
        }
        if(_S1441 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1440;
            _S1441 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S1449 = _s_diff_ctx_1->_S1399;
            if(!_s_diff_ctx_1->_S1400)
            {
                _runFlag_2 = false;
                _S1441 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S1450 = length_2(_S1449);
                    float3  raydir_29 = make_float3 ((_S1449 / make_float2 ((F32_max((_S1450), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1450))).x, (_S1449 / make_float2 ((F32_max((_S1450), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1450))).y, s_primal_ctx_cos_0(_S1450));
                    if(!is_ray_depth_3)
                    {
                        raydir_23 = raydir_29 / make_float3 (raydir_29.z);
                    }
                    else
                    {
                        raydir_23 = raydir_29;
                    }
                }
                else
                {
                    float3  raydir_30 = make_float3 (_S1449.x, _S1449.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_1(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S1437.primal_0.z) * raydir_23;
                float2  _S1451 = _s_diff_ctx_1->_S1401;
                bool _S1452 = !!_s_diff_ctx_1->_S1402;
                int _S1453;
                if(_S1452)
                {
                    if(is_fisheye_3)
                    {
                        float _S1454 = length_2(_S1451);
                        float3  raydir_31 = make_float3 ((_S1451 / make_float2 ((F32_max((_S1454), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1454))).x, (_S1451 / make_float2 ((F32_max((_S1454), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S1454))).y, s_primal_ctx_cos_0(_S1454));
                        if(!is_ray_depth_3)
                        {
                            raydir_24 = raydir_31 / make_float3 (raydir_31.z);
                        }
                        else
                        {
                            raydir_24 = raydir_31;
                        }
                    }
                    else
                    {
                        float3  raydir_32 = make_float3 (_S1451.x, _S1451.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_1(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1437.primal_0.w) * raydir_24;
                    _S1453 = int(2);
                }
                else
                {
                    _S1453 = int(0);
                    raydir_24 = _S1438;
                }
                if(_S1453 != int(2))
                {
                    _runFlag_3 = false;
                    _S1441 = _S1453;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1441 = int(1);
                }
                float3  _S1455 = raydir_23;
                _runFlag_3 = _S1452;
                raydir_23 = raydir_24;
                raydir_24 = _S1455;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S1438;
                raydir_24 = _S1438;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S1438;
            raydir_24 = _S1438;
        }
        float3  _S1456 = raydir_21;
        float3  _S1457 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S1442 = _S1446;
        raydir_23 = _S1457;
        raydir_24 = _S1456;
    }
    else
    {
        _S1441 = int(0);
        points_4[int(0)] = _S1438;
        points_4[int(1)] = _S1438;
        points_4[int(2)] = _S1438;
        points_4[int(3)] = _S1438;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S1438;
        raydir_22 = _S1438;
        _S1442 = false;
        raydir_23 = _S1438;
        raydir_24 = _S1438;
    }
    bool _S1458 = !(_S1441 != int(1));
    float3  _S1459;
    float3  _S1460;
    float3  _S1461;
    float3  _S1462;
    float3  _S1463;
    bool _S1464;
    if(_S1458)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1465 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1466 = s_primal_ctx_cross_0(dx_1, _S1465);
        bool _S1467 = (s_primal_ctx_dot_0(_S1466, _S1466)) != 0.0f;
        if(_S1467)
        {
            float _S1468 = length_1(_S1466);
            float3  _S1469 = make_float3 (_S1468);
            _S1459 = make_float3 (_S1468 * _S1468);
            _S1460 = _S1469;
        }
        else
        {
            _S1459 = _S1438;
            _S1460 = _S1438;
        }
        float3  _S1470 = _S1460;
        _S1464 = _S1467;
        _S1460 = _S1466;
        _S1461 = _S1470;
        _S1462 = dx_1;
        _S1463 = _S1465;
    }
    else
    {
        _S1464 = false;
        _S1459 = _S1438;
        _S1460 = _S1438;
        _S1461 = _S1438;
        _S1462 = _S1438;
        _S1463 = _S1438;
    }
    float4  _S1471 = make_float4 (0.0f);
    if(_S1458)
    {
        if(_S1464)
        {
            float3  _S1472 = _s_dOut_7 / _S1459;
            float3  _S1473 = _S1460 * - _S1472;
            float3  _S1474 = _S1461 * _S1472;
            float _S1475 = _S1473.x + _S1473.y + _S1473.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1476;
            (&_S1476)->primal_0 = _S1460;
            (&_S1476)->differential_0 = _S1438;
            s_bwd_length_impl_1(&_S1476, _S1475);
            _S1459 = _S1474 + _S1476.differential_0;
        }
        else
        {
            _S1459 = _s_dOut_7;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1477;
        (&_S1477)->primal_0 = _S1460;
        (&_S1477)->differential_0 = _S1438;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1478;
        (&_S1478)->primal_0 = _S1460;
        (&_S1478)->differential_0 = _S1438;
        s_bwd_prop_dot_0(&_S1477, &_S1478, 0.0f);
        float3  _S1479 = _S1478.differential_0 + _S1477.differential_0 + _S1459;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1480;
        (&_S1480)->primal_0 = _S1462;
        (&_S1480)->differential_0 = _S1438;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1481;
        (&_S1481)->primal_0 = _S1463;
        (&_S1481)->differential_0 = _S1438;
        s_bwd_prop_cross_0(&_S1480, &_S1481, _S1479);
        float3  s_diff_dy_T_1 = - _S1481.differential_0;
        float3  _S1482 = - s_diff_dy_T_1;
        float3  _S1483 = - _S1480.differential_0;
        FixedArray<float3 , 4>  _S1484;
        _S1484[int(0)] = _S1438;
        _S1484[int(1)] = _S1438;
        _S1484[int(2)] = _S1438;
        _S1484[int(3)] = _S1438;
        _S1484[int(2)] = _S1482;
        _S1484[int(3)] = s_diff_dy_T_1;
        _S1484[int(0)] = _S1483;
        _S1484[int(1)] = _S1480.differential_0;
        points_4[int(0)] = _S1484[int(0)];
        points_4[int(1)] = _S1484[int(1)];
        points_4[int(2)] = _S1484[int(2)];
        points_4[int(3)] = _S1484[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1438;
        points_4[int(1)] = _S1438;
        points_4[int(2)] = _S1438;
        points_4[int(3)] = _S1438;
    }
    float4  _S1485;
    if(_S1440)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1486 = points_4;
                FixedArray<float3 , 4>  _S1487 = points_4;
                FixedArray<float3 , 4>  _S1488 = points_4;
                FixedArray<float3 , 4>  _S1489 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1490 = raydir_21 * _S1489[int(3)];
                    float _S1491 = _S1490.x + _S1490.y + _S1490.z;
                    float4  _S1492 = _S1471;
                    *&((&_S1492)->w) = _S1491;
                    points_4[int(0)] = _S1486[int(0)];
                    points_4[int(1)] = _S1487[int(1)];
                    points_4[int(2)] = _S1488[int(2)];
                    points_4[int(3)] = _S1438;
                    _S1485 = _S1492;
                }
                else
                {
                    points_4[int(0)] = _S1486[int(0)];
                    points_4[int(1)] = _S1487[int(1)];
                    points_4[int(2)] = _S1488[int(2)];
                    points_4[int(3)] = _S1489[int(3)];
                    _S1485 = _S1471;
                }
                float3  _S1493 = raydir_22 * points_4[int(2)];
                float _S1494 = _S1493.x + _S1493.y + _S1493.z;
                FixedArray<float3 , 4>  _S1495 = points_4;
                FixedArray<float3 , 4>  _S1496 = points_4;
                float4  _S1497 = _S1471;
                *&((&_S1497)->z) = _S1494;
                float4  _S1498 = _S1485 + _S1497;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1495[int(1)];
                points_4[int(2)] = _S1438;
                points_4[int(3)] = _S1496[int(3)];
                _S1485 = _S1498;
            }
            else
            {
                FixedArray<float3 , 4>  _S1499 = points_4;
                FixedArray<float3 , 4>  _S1500 = points_4;
                FixedArray<float3 , 4>  _S1501 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1499[int(1)];
                points_4[int(2)] = _S1500[int(2)];
                points_4[int(3)] = _S1501[int(3)];
                _S1485 = _S1471;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1502 = points_4;
            FixedArray<float3 , 4>  _S1503 = points_4;
            FixedArray<float3 , 4>  _S1504 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1502[int(1)];
            points_4[int(2)] = _S1503[int(2)];
            points_4[int(3)] = _S1504[int(3)];
            _S1485 = _S1471;
        }
        if(_S1442)
        {
            FixedArray<float3 , 4>  _S1505 = points_4;
            float3  _S1506 = raydir_23 * points_4[int(1)];
            float _S1507 = _S1506.x + _S1506.y + _S1506.z;
            float4  _S1508 = _S1471;
            *&((&_S1508)->y) = _S1507;
            float4  _S1509 = _S1485 + _S1508;
            points_4[int(0)] = _S1438;
            points_4[int(1)] = _S1438;
            points_4[int(2)] = _S1438;
            points_4[int(3)] = _S1438;
            raydir_21 = _S1505[int(0)];
            _S1485 = _S1509;
        }
        else
        {
            FixedArray<float3 , 4>  _S1510 = points_4;
            FixedArray<float3 , 4>  _S1511 = points_4;
            FixedArray<float3 , 4>  _S1512 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1510[int(1)];
            points_4[int(2)] = _S1511[int(2)];
            points_4[int(3)] = _S1512[int(3)];
            raydir_21 = _S1438;
        }
        float3  _S1513 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S1514 = _S1513.x + _S1513.y + _S1513.z;
        float4  _S1515 = _S1471;
        *&((&_S1515)->x) = _S1514;
        _S1485 = _S1485 + _S1515;
    }
    else
    {
        _S1485 = _S1471;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1485;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1516, float4  _S1517, FixedArray<float, 10>  * _S1518, bool _S1519, bool _S1520, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1521, float3  _S1522)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1523;
    float3  _S1524 = s_primal_ctx_depth_to_normal_0(_S1516, _S1517, _S1518, _S1519, _S1520, (*_S1521).primal_0, &_S1523);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1525 = _S1523;
    s_bwd_prop_depth_to_normal_0(_S1516, _S1517, _S1518, _S1519, _S1520, _S1521, _S1522, &_S1525);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1526 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1526;
    FixedArray<float, 10>  _S1527 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S1527, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S1528 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S1528;
    FixedArray<float, 10>  _S1529 = dist_coeffs_6;
    bool _S1530 = undistort_point_0(_S1528, &_S1529, int(12), &uv_10);
    if(!_S1530)
    {
        return 0.0f;
    }
    float3  raydir_33;
    if(is_fisheye_5)
    {
        float theta_5 = length_2(uv_10);
        float3  raydir_34 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_10 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        raydir_33 = raydir_34 / make_float3 (raydir_34.z);
    }
    else
    {
        raydir_33 = make_float3 (uv_10.x, uv_10.y, 1.0f);
    }
    return float((F32_sign((raydir_33.z)))) / length_1(raydir_33);
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_24, float dOut_19)
{
    float _S1531 = (F32_exp2(((*dpx_24).primal_0))) * 0.69314718246459961f * dOut_19;
    dpx_24->primal_0 = (*dpx_24).primal_0;
    dpx_24->differential_0 = _S1531;
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
    PPISPParams_0 _S1532 = p_0;
    float _S1533 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S1534 = (pix_coord_0.x - image_center_0.x) / _S1533;
    float _S1535 = (pix_coord_0.y - image_center_0.y) / _S1533;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_2 = _S1534 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S1535 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_3 = dx_2 * dx_2 + dy_0 * dy_0;
    float r4_0 = r2_3 * r2_3;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_3) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_3 + 1.0f, 0.0f, 1.0f);
    float dx_3 = _S1534 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S1535 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_4 = dx_3 * dx_3 + dy_1 * dy_1;
    float r4_1 = r2_4 * r2_4;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_4) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S1534 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S1535 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_5 = dx_4 * dx_4 + dy_2 * dy_2;
    float r4_2 = r2_5 * r2_5;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_5) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float3  _S1536 = rgb_out_0;
    float2  bd_0 = mul_3(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_3(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_3(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_3(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S1537 = 0.3333333432674408f + nd_0.x;
    float _S1538 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_1 = mul_1(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1538, 1.0f, 0.0f, - _S1537, - _S1538, _S1537, 0.0f), T_0);
    float3  r0_0 = make_float3 (M_1.rows[int(0)].x, M_1.rows[int(0)].y, M_1.rows[int(0)].z);
    float3  r1_0 = make_float3 (M_1.rows[int(1)].x, M_1.rows[int(1)].y, M_1.rows[int(1)].z);
    float3  r2_6 = make_float3 (M_1.rows[int(2)].x, M_1.rows[int(2)].y, M_1.rows[int(2)].z);
    float3  lambda_v_0 = cross_0(r0_0, r1_0);
    float3  lambda_v_1;
    if((dot_0(lambda_v_0, lambda_v_0)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_2 = cross_0(r0_0, r2_6);
        if((dot_0(lambda_v_2, lambda_v_2)) < 9.99999968265522539e-21f)
        {
            lambda_v_1 = cross_0(r1_0, r2_6);
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
    float _S1539 = _S1536.x;
    float _S1540 = _S1536.y;
    float intensity_0 = _S1539 + _S1540 + _S1536.z;
    float3  rgi_out_0 = mul_2(H_1, make_float3 (_S1539, _S1540, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / (rgi_out_0.z + 0.00000999999974738f));
    float _S1541 = rgi_out_1.x;
    float _S1542 = rgi_out_1.y;
    float3  _S1543 = clamp_1(make_float3 (_S1541, _S1542, rgi_out_1.z - _S1541 - _S1542), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S1544 = _S1543.x;
    float _S1545 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(0)].toe_0))))));
    float _S1546 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(0)].shoulder_0))))));
    float _S1547 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(0)].gamma_0))))));
    float _S1548 = 1.0f / (1.0f + (F32_exp((- _S1532.crf_params_1[int(0)].center_0))));
    float a_1 = _S1546 * _S1548 / lerp_0(_S1545, _S1546, _S1548);
    float b_2 = 1.0f - a_1;
    float y_11;
    if(_S1544 <= _S1548)
    {
        y_11 = a_1 * (F32_pow((_S1544 / _S1548), (_S1545)));
    }
    else
    {
        y_11 = 1.0f - b_2 * (F32_pow(((1.0f - _S1544) / (1.0f - _S1548)), (_S1546)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_11)))), (_S1547)));
    float _S1549 = _S1543.y;
    float _S1550 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(1)].toe_0))))));
    float _S1551 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(1)].shoulder_0))))));
    float _S1552 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(1)].gamma_0))))));
    float _S1553 = 1.0f / (1.0f + (F32_exp((- _S1532.crf_params_1[int(1)].center_0))));
    float a_2 = _S1551 * _S1553 / lerp_0(_S1550, _S1551, _S1553);
    float b_3 = 1.0f - a_2;
    if(_S1549 <= _S1553)
    {
        y_11 = a_2 * (F32_pow((_S1549 / _S1553), (_S1550)));
    }
    else
    {
        y_11 = 1.0f - b_3 * (F32_pow(((1.0f - _S1549) / (1.0f - _S1553)), (_S1551)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_11)))), (_S1552)));
    float _S1554 = _S1543.z;
    float _S1555 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(2)].toe_0))))));
    float _S1556 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(2)].shoulder_0))))));
    float _S1557 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S1532.crf_params_1[int(2)].gamma_0))))));
    float _S1558 = 1.0f / (1.0f + (F32_exp((- _S1532.crf_params_1[int(2)].center_0))));
    float a_3 = _S1556 * _S1558 / lerp_0(_S1555, _S1556, _S1558);
    float b_4 = 1.0f - a_3;
    if(_S1554 <= _S1558)
    {
        y_11 = a_3 * (F32_pow((_S1554 / _S1558), (_S1555)));
    }
    else
    {
        y_11 = 1.0f - b_4 * (F32_pow(((1.0f - _S1554) / (1.0f - _S1558)), (_S1556)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_11)))), (_S1557)));
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
    PPISPParamsRQS_0 _S1559 = p_1;
    float _S1560 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S1561 = (pix_coord_1.x - image_center_1.x) / _S1560;
    float _S1562 = (pix_coord_1.y - image_center_1.y) / _S1560;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_5 = _S1561 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S1562 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_7 = dx_5 * dx_5 + dy_3 * dy_3;
    float r4_3 = r2_7 * r2_7;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_7) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_7 + 1.0f, 0.0f, 1.0f);
    float dx_6 = _S1561 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S1562 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_8 = dx_6 * dx_6 + dy_4 * dy_4;
    float r4_4 = r2_8 * r2_8;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_8) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_8 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S1561 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S1562 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_9 = dx_7 * dx_7 + dy_5 * dy_5;
    float r4_5 = r2_9 * r2_9;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_9) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float3  _S1563 = rgb_out_2;
    float2  bd_1 = mul_3(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_3(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_3(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_3(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S1564 = 0.3333333432674408f + nd_1.x;
    float _S1565 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_2 = mul_1(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1565, 1.0f, 0.0f, - _S1564, - _S1565, _S1564, 0.0f), T_1);
    float3  r0_1 = make_float3 (M_2.rows[int(0)].x, M_2.rows[int(0)].y, M_2.rows[int(0)].z);
    float3  r1_1 = make_float3 (M_2.rows[int(1)].x, M_2.rows[int(1)].y, M_2.rows[int(1)].z);
    float3  r2_10 = make_float3 (M_2.rows[int(2)].x, M_2.rows[int(2)].y, M_2.rows[int(2)].z);
    float3  lambda_v_3 = cross_0(r0_1, r1_1);
    float3  lambda_v_4;
    if((dot_0(lambda_v_3, lambda_v_3)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_5 = cross_0(r0_1, r2_10);
        if((dot_0(lambda_v_5, lambda_v_5)) < 9.99999968265522539e-21f)
        {
            lambda_v_4 = cross_0(r1_1, r2_10);
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
    float _S1566 = _S1563.x;
    float _S1567 = _S1563.y;
    float intensity_1 = _S1566 + _S1567 + _S1563.z;
    float3  rgi_out_2 = mul_2(H_3, make_float3 (_S1566, _S1567, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / (rgi_out_2.z + 0.00000999999974738f));
    float _S1568 = rgi_out_3.x;
    float _S1569 = rgi_out_3.y;
    float3  _S1570 = clamp_1(make_float3 (_S1568, _S1569, rgi_out_3.z - _S1568 - _S1569), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S1571 = _S1570.x;
    float g0_1 = (F32_exp((_S1559.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S1559.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S1559.crf_params_0[int(0)].gc_0)));
    float y_12;
    if(_S1571 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S1571 / x0_1;
        float _S1572 = 1.0f - t0_0;
        y_12 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S1572) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S1572);
    }
    else
    {
        float _S1573 = 1.0f - y0_1;
        float _S1574 = 1.0f - x0_1;
        float s1_6 = _S1573 / _S1574;
        float t1_0 = (_S1571 - x0_1) / _S1574;
        float _S1575 = 1.0f - t1_0;
        y_12 = y0_1 + _S1573 * (s1_6 * t1_0 * t1_0 + gc_1 * t1_0 * _S1575) / (s1_6 + (gc_1 + g1_1 - 2.0f * s1_6) * t1_0 * _S1575);
    }
    *&((&rgb_out_3)->x) = y_12;
    float _S1576 = _S1570.y;
    float g0_2 = (F32_exp((_S1559.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S1559.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S1559.crf_params_0[int(1)].gc_0)));
    if(_S1576 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S1576 / x0_2;
        float _S1577 = 1.0f - t0_1;
        y_12 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S1577) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S1577);
    }
    else
    {
        float _S1578 = 1.0f - y0_2;
        float _S1579 = 1.0f - x0_2;
        float s1_7 = _S1578 / _S1579;
        float t1_1 = (_S1576 - x0_2) / _S1579;
        float _S1580 = 1.0f - t1_1;
        y_12 = y0_2 + _S1578 * (s1_7 * t1_1 * t1_1 + gc_2 * t1_1 * _S1580) / (s1_7 + (gc_2 + g1_2 - 2.0f * s1_7) * t1_1 * _S1580);
    }
    *&((&rgb_out_3)->y) = y_12;
    float _S1581 = _S1570.z;
    float g0_3 = (F32_exp((_S1559.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S1559.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S1559.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S1559.crf_params_0[int(2)].gc_0)));
    if(_S1581 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S1581 / x0_3;
        float _S1582 = 1.0f - t0_2;
        y_12 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S1582) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S1582);
    }
    else
    {
        float _S1583 = 1.0f - y0_3;
        float _S1584 = 1.0f - x0_3;
        float s1_8 = _S1583 / _S1584;
        float t1_2 = (_S1581 - x0_3) / _S1584;
        float _S1585 = 1.0f - t1_2;
        y_12 = y0_3 + _S1583 * (s1_8 * t1_2 * t1_2 + gc_3 * t1_2 * _S1585) / (s1_8 + (gc_3 + g1_3 - 2.0f * s1_8) * t1_2 * _S1585);
    }
    *&((&rgb_out_3)->z) = y_12;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S1586)
{
    return (F32_exp2((_S1586)));
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S1587, float2  _S1588)
{
    return mul_3(_S1587, _S1588);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1589, Matrix<float, 3, 3>  _S1590)
{
    return mul_1(_S1589, _S1590);
}

inline __device__ float s_primal_ctx_abs_0(float _S1591)
{
    return (F32_abs((_S1591)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S1592, float3  _S1593)
{
    return mul_2(_S1592, _S1593);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S1594, float3  _S1595, float3  _S1596)
{
    return clamp_1(_S1594, _S1595, _S1596);
}

inline __device__ float s_primal_ctx_lerp_0(float _S1597, float _S1598, float _S1599)
{
    return lerp_0(_S1597, _S1598, _S1599);
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1600, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1601, float3  _S1602)
{
    _d_mul_0(_S1600, _S1601, _S1602);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_float_0 * _S1603, float _S1604)
{
    _d_abs_0(_S1603, _S1604);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1605, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1606, Matrix<float, 3, 3>  _S1607)
{
    mul_0(_S1605, _S1606, _S1607);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S1608, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1609, float2  _S1610)
{
    _d_mul_1(_S1608, _S1609, _S1610);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1611, float _S1612)
{
    _d_exp2_0(_S1611, _S1612);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_8)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1613 = *dprgb_in_0;
    float3  _S1614 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1615 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1616 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1617 = {
        _S1616, _S1616, _S1616
    };
    float2  _S1618 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1619 = { _S1618, _S1618, _S1618, _S1618 };
    CRFPPISPChannelParams_0 _S1620 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1621 = {
        _S1620, _S1620, _S1620
    };
    PPISPParams_0 _S1622;
    (&_S1622)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S1622)->vignette_params_1 = _S1617;
    (&_S1622)->color_params_1 = _S1619;
    (&_S1622)->crf_params_1 = _S1621;
    (&(&_S1622)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S1622)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S1623 = dpparams_0->primal_0[int(3)];
    (&(&_S1622)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S1624 = dpparams_0->primal_0[int(4)];
    (&(&_S1622)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S1625 = dpparams_0->primal_0[int(5)];
    (&(&_S1622)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S1622)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S1622)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S1626 = dpparams_0->primal_0[int(8)];
    (&(&_S1622)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S1627 = dpparams_0->primal_0[int(9)];
    (&(&_S1622)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S1628 = dpparams_0->primal_0[int(10)];
    (&(&_S1622)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S1622)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S1622)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S1629 = dpparams_0->primal_0[int(13)];
    (&(&_S1622)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S1630 = dpparams_0->primal_0[int(14)];
    (&(&_S1622)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S1631 = dpparams_0->primal_0[int(15)];
    (&(&_S1622)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S1622)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S1622)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S1622)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S1622)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S1622)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S1622)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S1622)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S1622)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S1632 = dpparams_0->primal_0[int(24)];
    (&(&_S1622)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S1633 = dpparams_0->primal_0[int(25)];
    (&(&_S1622)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S1634 = dpparams_0->primal_0[int(26)];
    (&(&_S1622)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S1635 = dpparams_0->primal_0[int(27)];
    (&(&_S1622)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S1636 = dpparams_0->primal_0[int(28)];
    (&(&_S1622)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S1637 = dpparams_0->primal_0[int(29)];
    (&(&_S1622)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S1638 = dpparams_0->primal_0[int(30)];
    (&(&_S1622)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S1639 = dpparams_0->primal_0[int(31)];
    (&(&_S1622)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S1640 = dpparams_0->primal_0[int(32)];
    (&(&_S1622)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S1641 = dpparams_0->primal_0[int(33)];
    (&(&_S1622)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S1642 = dpparams_0->primal_0[int(34)];
    (&(&_S1622)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S1643 = dpparams_0->primal_0[int(35)];
    (&(&_S1622)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S1644 = _S1622;
    float _S1645 = s_primal_ctx_exp2_0(_S1622.exposure_1);
    float3  _S1646 = make_float3 (_S1645);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S1645);
    float _S1647 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S1648 = (pix_coord_2.x - image_center_2.x) / _S1647;
    float _S1649 = (pix_coord_2.y - image_center_2.y) / _S1647;
    float dx_8 = _S1648 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S1649 - dpparams_0->primal_0[int(2)];
    float r2_11 = dx_8 * dx_8 + dy_6 * dy_6;
    float r4_6 = r2_11 * r2_11;
    float r6_0 = r4_6 * r2_11;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_11 + 1.0f;
    float _S1650 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S1651 = rgb_out_4.x * _S1650;
    float3  _S1652 = rgb_out_4;
    *&((&_S1652)->x) = _S1651;
    float dx_9 = _S1648 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S1649 - dpparams_0->primal_0[int(7)];
    float r2_12 = dx_9 * dx_9 + dy_7 * dy_7;
    float r4_7 = r2_12 * r2_12;
    float r6_1 = r4_7 * r2_12;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_12 + 1.0f;
    float _S1653 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S1652)->y) = rgb_out_4.y * _S1653;
    float dx_10 = _S1648 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S1649 - dpparams_0->primal_0[int(12)];
    float r2_13 = dx_10 * dx_10 + dy_8 * dy_8;
    float r4_8 = r2_13 * r2_13;
    float r6_2 = r4_8 * r2_13;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_13 + 1.0f;
    float _S1654 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S1652)->z) = rgb_out_4.z * _S1654;
    PPISPParams_0 _S1655 = _S1622;
    float2  _S1656 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1622.color_params_1.b_0);
    float2  _S1657 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1622.color_params_1.r_0);
    float2  _S1658 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1622.color_params_1.g_0);
    float2  _S1659 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1622.color_params_1.n_0);
    float _S1660 = 0.3333333432674408f + _S1659.x;
    float _S1661 = 0.3333333432674408f + _S1659.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S1656.x, 1.0f + _S1657.x, _S1658.x, _S1656.y, _S1657.y, 1.0f + _S1658.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1661, 1.0f, 0.0f, - _S1660, - _S1661, _S1660, 0.0f);
    Matrix<float, 3, 3>  _S1662 = s_primal_ctx_mul_1(skew_0, T_2);
    float3  r0_2 = make_float3 (_S1662.rows[int(0)].x, _S1662.rows[int(0)].y, _S1662.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S1662.rows[int(1)].x, _S1662.rows[int(1)].y, _S1662.rows[int(1)].z);
    float3  r2_14 = make_float3 (_S1662.rows[int(2)].x, _S1662.rows[int(2)].y, _S1662.rows[int(2)].z);
    float3  _S1663 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S1664 = (s_primal_ctx_dot_0(_S1663, _S1663)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S1665;
    bool _S1666;
    if(_S1664)
    {
        float3  _S1667 = s_primal_ctx_cross_0(r0_2, r2_14);
        bool _S1668 = (s_primal_ctx_dot_0(_S1667, _S1667)) < 9.99999968265522539e-21f;
        if(_S1668)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_14);
        }
        else
        {
            lambda_v_6 = _S1667;
        }
        _S1666 = _S1668;
        _S1665 = _S1667;
    }
    else
    {
        lambda_v_6 = _S1663;
        _S1666 = false;
        _S1665 = _S1614;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S1669 = s_primal_ctx_mul_1(T_2, D_0);
    Matrix<float, 3, 3>  _S1670 = s_primal_ctx_mul_1(_S1669, S_inv_0);
    bool _S1671 = (s_primal_ctx_abs_0(_S1670.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S1672;
    float _S1673;
    if(_S1671)
    {
        float inv_s_0 = 1.0f / _S1670.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1674 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S1675 = _S1670.rows[int(2)].z * _S1670.rows[int(2)].z;
        H_4 = _S1670 * makeMatrix<float, 3, 3> (inv_s_0);
        _S1672 = _S1674;
        _S1673 = _S1675;
    }
    else
    {
        H_4 = _S1670;
        _S1672 = _S1615;
        _S1673 = 0.0f;
    }
    float _S1676 = _S1652.x;
    float _S1677 = _S1652.y;
    float intensity_2 = _S1676 + _S1677 + _S1652.z;
    float3  rgi_in_0 = make_float3 (_S1676, _S1677, intensity_2);
    float3  _S1678 = s_primal_ctx_mul_2(H_4, rgi_in_0);
    float _S1679 = _S1678.z + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S1679;
    float3  _S1680 = make_float3 (norm_factor_0);
    float _S1681 = _S1679 * _S1679;
    float3  rgi_out_4 = _S1678 * make_float3 (norm_factor_0);
    float _S1682 = rgi_out_4.x;
    float _S1683 = rgi_out_4.y;
    float3  _S1684 = make_float3 (_S1682, _S1683, rgi_out_4.z - _S1682 - _S1683);
    float3  _S1685 = make_float3 (0.0f);
    float3  _S1686 = make_float3 (1.0f);
    float3  _S1687 = s_primal_ctx_clamp_1(_S1684, _S1685, _S1686);
    float _S1688 = _S1687.x;
    float _S1689 = 1.0f + s_primal_ctx_exp_0(_S1632);
    float _S1690 = 0.30000001192092896f + s_primal_ctx_log_0(_S1689);
    float _S1691 = 1.0f + s_primal_ctx_exp_0(_S1633);
    float _S1692 = 0.30000001192092896f + s_primal_ctx_log_0(_S1691);
    float _S1693 = 1.0f + s_primal_ctx_exp_0(_S1634);
    float _S1694 = 0.10000000149011612f + s_primal_ctx_log_0(_S1693);
    float _S1695 = - _S1635;
    float _S1696 = 1.0f + s_primal_ctx_exp_0(_S1695);
    float _S1697 = 1.0f / _S1696;
    float _S1698 = _S1696 * _S1696;
    float _S1699 = s_primal_ctx_lerp_0(_S1690, _S1692, _S1697);
    float _S1700 = _S1692 * _S1697;
    float a_4 = _S1700 / _S1699;
    float _S1701 = _S1699 * _S1699;
    float b_5 = 1.0f - a_4;
    bool _S1702 = _S1688 <= _S1697;
    float y_13;
    float _S1703;
    float _S1704;
    float _S1705;
    float _S1706;
    float _S1707;
    float _S1708;
    float _S1709;
    float _S1710;
    if(_S1702)
    {
        float _S1711 = _S1688 / _S1697;
        float _S1712 = _S1697 * _S1697;
        float _S1713 = s_primal_ctx_pow_0(_S1711, _S1690);
        y_13 = a_4 * _S1713;
        _S1703 = _S1713;
        _S1704 = _S1711;
        _S1705 = _S1712;
        _S1706 = 0.0f;
        _S1707 = 0.0f;
        _S1708 = 0.0f;
        _S1709 = 0.0f;
        _S1710 = 0.0f;
    }
    else
    {
        float _S1714 = 1.0f - _S1688;
        float _S1715 = 1.0f - _S1697;
        float _S1716 = _S1714 / _S1715;
        float _S1717 = _S1715 * _S1715;
        float _S1718 = s_primal_ctx_pow_0(_S1716, _S1692);
        y_13 = 1.0f - b_5 * _S1718;
        _S1703 = 0.0f;
        _S1704 = 0.0f;
        _S1705 = 0.0f;
        _S1706 = _S1718;
        _S1707 = _S1716;
        _S1708 = _S1717;
        _S1709 = _S1714;
        _S1710 = _S1715;
    }
    float _S1719 = (F32_max((0.0f), (y_13)));
    float _S1720 = _S1687.y;
    float _S1721 = 1.0f + s_primal_ctx_exp_0(_S1636);
    float _S1722 = 0.30000001192092896f + s_primal_ctx_log_0(_S1721);
    float _S1723 = 1.0f + s_primal_ctx_exp_0(_S1637);
    float _S1724 = 0.30000001192092896f + s_primal_ctx_log_0(_S1723);
    float _S1725 = 1.0f + s_primal_ctx_exp_0(_S1638);
    float _S1726 = 0.10000000149011612f + s_primal_ctx_log_0(_S1725);
    float _S1727 = - _S1639;
    float _S1728 = 1.0f + s_primal_ctx_exp_0(_S1727);
    float _S1729 = 1.0f / _S1728;
    float _S1730 = _S1728 * _S1728;
    float _S1731 = s_primal_ctx_lerp_0(_S1722, _S1724, _S1729);
    float _S1732 = _S1724 * _S1729;
    float a_5 = _S1732 / _S1731;
    float _S1733 = _S1731 * _S1731;
    float b_6 = 1.0f - a_5;
    bool _S1734 = _S1720 <= _S1729;
    float y_14;
    float _S1735;
    float _S1736;
    float _S1737;
    float _S1738;
    float _S1739;
    float _S1740;
    float _S1741;
    float _S1742;
    if(_S1734)
    {
        float _S1743 = _S1720 / _S1729;
        float _S1744 = _S1729 * _S1729;
        float _S1745 = s_primal_ctx_pow_0(_S1743, _S1722);
        y_14 = a_5 * _S1745;
        _S1735 = _S1745;
        _S1736 = _S1743;
        _S1737 = _S1744;
        _S1738 = 0.0f;
        _S1739 = 0.0f;
        _S1740 = 0.0f;
        _S1741 = 0.0f;
        _S1742 = 0.0f;
    }
    else
    {
        float _S1746 = 1.0f - _S1720;
        float _S1747 = 1.0f - _S1729;
        float _S1748 = _S1746 / _S1747;
        float _S1749 = _S1747 * _S1747;
        float _S1750 = s_primal_ctx_pow_0(_S1748, _S1724);
        y_14 = 1.0f - b_6 * _S1750;
        _S1735 = 0.0f;
        _S1736 = 0.0f;
        _S1737 = 0.0f;
        _S1738 = _S1750;
        _S1739 = _S1748;
        _S1740 = _S1749;
        _S1741 = _S1746;
        _S1742 = _S1747;
    }
    float _S1751 = (F32_max((0.0f), (y_14)));
    float _S1752 = _S1687.z;
    float _S1753 = 1.0f + s_primal_ctx_exp_0(_S1640);
    float _S1754 = 0.30000001192092896f + s_primal_ctx_log_0(_S1753);
    float _S1755 = 1.0f + s_primal_ctx_exp_0(_S1641);
    float _S1756 = 0.30000001192092896f + s_primal_ctx_log_0(_S1755);
    float _S1757 = 1.0f + s_primal_ctx_exp_0(_S1642);
    float _S1758 = 0.10000000149011612f + s_primal_ctx_log_0(_S1757);
    float _S1759 = - _S1643;
    float _S1760 = 1.0f + s_primal_ctx_exp_0(_S1759);
    float _S1761 = 1.0f / _S1760;
    float _S1762 = _S1760 * _S1760;
    float _S1763 = s_primal_ctx_lerp_0(_S1754, _S1756, _S1761);
    float _S1764 = _S1756 * _S1761;
    float a_6 = _S1764 / _S1763;
    float _S1765 = _S1763 * _S1763;
    float b_7 = 1.0f - a_6;
    bool _S1766 = _S1752 <= _S1761;
    float y_15;
    float _S1767;
    float _S1768;
    float _S1769;
    float _S1770;
    float _S1771;
    float _S1772;
    float _S1773;
    float _S1774;
    if(_S1766)
    {
        float _S1775 = _S1752 / _S1761;
        float _S1776 = _S1761 * _S1761;
        float _S1777 = s_primal_ctx_pow_0(_S1775, _S1754);
        y_15 = a_6 * _S1777;
        _S1767 = _S1777;
        _S1768 = _S1775;
        _S1769 = _S1776;
        _S1770 = 0.0f;
        _S1771 = 0.0f;
        _S1772 = 0.0f;
        _S1773 = 0.0f;
        _S1774 = 0.0f;
    }
    else
    {
        float _S1778 = 1.0f - _S1752;
        float _S1779 = 1.0f - _S1761;
        float _S1780 = _S1778 / _S1779;
        float _S1781 = _S1779 * _S1779;
        float _S1782 = s_primal_ctx_pow_0(_S1780, _S1756);
        y_15 = 1.0f - b_7 * _S1782;
        _S1767 = 0.0f;
        _S1768 = 0.0f;
        _S1769 = 0.0f;
        _S1770 = _S1782;
        _S1771 = _S1780;
        _S1772 = _S1781;
        _S1773 = _S1778;
        _S1774 = _S1779;
    }
    float _S1783 = (F32_max((0.0f), (y_15)));
    DiffPair_float_0 _S1784;
    (&_S1784)->primal_0 = _S1783;
    (&_S1784)->differential_0 = 0.0f;
    DiffPair_float_0 _S1785;
    (&_S1785)->primal_0 = _S1758;
    (&_S1785)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1784, &_S1785, _s_dOut_8.z);
    DiffPair_float_0 _S1786 = _S1785;
    DiffPair_float_0 _S1787;
    (&_S1787)->primal_0 = 0.0f;
    (&_S1787)->differential_0 = 0.0f;
    DiffPair_float_0 _S1788;
    (&_S1788)->primal_0 = y_15;
    (&_S1788)->differential_0 = 0.0f;
    _d_max_0(&_S1787, &_S1788, _S1784.differential_0);
    DiffPair_float_0 _S1789 = _S1788;
    if(_S1766)
    {
        float _S1790 = a_6 * _S1789.differential_0;
        float _S1791 = _S1767 * _S1789.differential_0;
        DiffPair_float_0 _S1792;
        (&_S1792)->primal_0 = _S1768;
        (&_S1792)->differential_0 = 0.0f;
        DiffPair_float_0 _S1793;
        (&_S1793)->primal_0 = _S1754;
        (&_S1793)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1792, &_S1793, _S1790);
        float _S1794 = _S1792.differential_0 / _S1769;
        float _S1795 = _S1752 * - _S1794;
        float _S1796 = _S1761 * _S1794;
        y_15 = 0.0f;
        _S1767 = _S1791;
        _S1768 = _S1795;
        _S1769 = 0.0f;
        _S1770 = _S1793.differential_0;
        _S1771 = _S1796;
    }
    else
    {
        float _S1797 = - _S1789.differential_0;
        float _S1798 = b_7 * _S1797;
        float _S1799 = _S1770 * _S1797;
        DiffPair_float_0 _S1800;
        (&_S1800)->primal_0 = _S1771;
        (&_S1800)->differential_0 = 0.0f;
        DiffPair_float_0 _S1801;
        (&_S1801)->primal_0 = _S1756;
        (&_S1801)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1800, &_S1801, _S1798);
        float _S1802 = _S1800.differential_0 / _S1772;
        float _S1803 = - (_S1773 * - _S1802);
        float _S1804 = - (_S1774 * _S1802);
        y_15 = _S1799;
        _S1767 = 0.0f;
        _S1768 = _S1803;
        _S1769 = _S1801.differential_0;
        _S1770 = 0.0f;
        _S1771 = _S1804;
    }
    float _S1805 = (- y_15 + _S1767) / _S1765;
    float _S1806 = _S1764 * - _S1805;
    float _S1807 = _S1763 * _S1805;
    float _S1808 = _S1756 * _S1807;
    float _S1809 = _S1761 * _S1807;
    DiffPair_float_0 _S1810;
    (&_S1810)->primal_0 = _S1754;
    (&_S1810)->differential_0 = 0.0f;
    DiffPair_float_0 _S1811;
    (&_S1811)->primal_0 = _S1756;
    (&_S1811)->differential_0 = 0.0f;
    DiffPair_float_0 _S1812;
    (&_S1812)->primal_0 = _S1761;
    (&_S1812)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1810, &_S1811, &_S1812, _S1806);
    float _S1813 = - ((_S1808 + _S1812.differential_0 + _S1768) / _S1762);
    DiffPair_float_0 _S1814;
    (&_S1814)->primal_0 = _S1759;
    (&_S1814)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1814, _S1813);
    float _S1815 = - _S1814.differential_0;
    DiffPair_float_0 _S1816;
    (&_S1816)->primal_0 = _S1757;
    (&_S1816)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1816, _S1786.differential_0);
    DiffPair_float_0 _S1817;
    (&_S1817)->primal_0 = _S1642;
    (&_S1817)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1817, _S1816.differential_0);
    DiffPair_float_0 _S1818 = _S1817;
    float _S1819 = _S1809 + _S1811.differential_0 + _S1769;
    DiffPair_float_0 _S1820;
    (&_S1820)->primal_0 = _S1755;
    (&_S1820)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1820, _S1819);
    DiffPair_float_0 _S1821;
    (&_S1821)->primal_0 = _S1641;
    (&_S1821)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1821, _S1820.differential_0);
    DiffPair_float_0 _S1822 = _S1821;
    float _S1823 = _S1810.differential_0 + _S1770;
    DiffPair_float_0 _S1824;
    (&_S1824)->primal_0 = _S1753;
    (&_S1824)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1824, _S1823);
    DiffPair_float_0 _S1825;
    (&_S1825)->primal_0 = _S1640;
    (&_S1825)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1825, _S1824.differential_0);
    DiffPair_float_0 _S1826 = _S1825;
    float3  _S1827 = make_float3 (0.0f, 0.0f, _S1771);
    DiffPair_float_0 _S1828;
    (&_S1828)->primal_0 = _S1751;
    (&_S1828)->differential_0 = 0.0f;
    DiffPair_float_0 _S1829;
    (&_S1829)->primal_0 = _S1726;
    (&_S1829)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1828, &_S1829, _s_dOut_8.y);
    DiffPair_float_0 _S1830 = _S1829;
    DiffPair_float_0 _S1831;
    (&_S1831)->primal_0 = 0.0f;
    (&_S1831)->differential_0 = 0.0f;
    DiffPair_float_0 _S1832;
    (&_S1832)->primal_0 = y_14;
    (&_S1832)->differential_0 = 0.0f;
    _d_max_0(&_S1831, &_S1832, _S1828.differential_0);
    DiffPair_float_0 _S1833 = _S1832;
    if(_S1734)
    {
        float _S1834 = a_5 * _S1833.differential_0;
        float _S1835 = _S1735 * _S1833.differential_0;
        DiffPair_float_0 _S1836;
        (&_S1836)->primal_0 = _S1736;
        (&_S1836)->differential_0 = 0.0f;
        DiffPair_float_0 _S1837;
        (&_S1837)->primal_0 = _S1722;
        (&_S1837)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1836, &_S1837, _S1834);
        float _S1838 = _S1836.differential_0 / _S1737;
        float _S1839 = _S1720 * - _S1838;
        float _S1840 = _S1729 * _S1838;
        y_14 = 0.0f;
        _S1735 = _S1835;
        _S1736 = _S1839;
        _S1737 = 0.0f;
        _S1738 = _S1837.differential_0;
        _S1739 = _S1840;
    }
    else
    {
        float _S1841 = - _S1833.differential_0;
        float _S1842 = b_6 * _S1841;
        float _S1843 = _S1738 * _S1841;
        DiffPair_float_0 _S1844;
        (&_S1844)->primal_0 = _S1739;
        (&_S1844)->differential_0 = 0.0f;
        DiffPair_float_0 _S1845;
        (&_S1845)->primal_0 = _S1724;
        (&_S1845)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1844, &_S1845, _S1842);
        float _S1846 = _S1844.differential_0 / _S1740;
        float _S1847 = - (_S1741 * - _S1846);
        float _S1848 = - (_S1742 * _S1846);
        y_14 = _S1843;
        _S1735 = 0.0f;
        _S1736 = _S1847;
        _S1737 = _S1845.differential_0;
        _S1738 = 0.0f;
        _S1739 = _S1848;
    }
    float _S1849 = (- y_14 + _S1735) / _S1733;
    float _S1850 = _S1732 * - _S1849;
    float _S1851 = _S1731 * _S1849;
    float _S1852 = _S1724 * _S1851;
    float _S1853 = _S1729 * _S1851;
    DiffPair_float_0 _S1854;
    (&_S1854)->primal_0 = _S1722;
    (&_S1854)->differential_0 = 0.0f;
    DiffPair_float_0 _S1855;
    (&_S1855)->primal_0 = _S1724;
    (&_S1855)->differential_0 = 0.0f;
    DiffPair_float_0 _S1856;
    (&_S1856)->primal_0 = _S1729;
    (&_S1856)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1854, &_S1855, &_S1856, _S1850);
    float _S1857 = - ((_S1852 + _S1856.differential_0 + _S1736) / _S1730);
    DiffPair_float_0 _S1858;
    (&_S1858)->primal_0 = _S1727;
    (&_S1858)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1858, _S1857);
    float _S1859 = - _S1858.differential_0;
    DiffPair_float_0 _S1860;
    (&_S1860)->primal_0 = _S1725;
    (&_S1860)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1860, _S1830.differential_0);
    DiffPair_float_0 _S1861;
    (&_S1861)->primal_0 = _S1638;
    (&_S1861)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1861, _S1860.differential_0);
    DiffPair_float_0 _S1862 = _S1861;
    float _S1863 = _S1853 + _S1855.differential_0 + _S1737;
    DiffPair_float_0 _S1864;
    (&_S1864)->primal_0 = _S1723;
    (&_S1864)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1864, _S1863);
    DiffPair_float_0 _S1865;
    (&_S1865)->primal_0 = _S1637;
    (&_S1865)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1865, _S1864.differential_0);
    DiffPair_float_0 _S1866 = _S1865;
    float _S1867 = _S1854.differential_0 + _S1738;
    DiffPair_float_0 _S1868;
    (&_S1868)->primal_0 = _S1721;
    (&_S1868)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1868, _S1867);
    DiffPair_float_0 _S1869;
    (&_S1869)->primal_0 = _S1636;
    (&_S1869)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1869, _S1868.differential_0);
    DiffPair_float_0 _S1870 = _S1869;
    float3  _S1871 = _S1827 + make_float3 (0.0f, _S1739, 0.0f);
    DiffPair_float_0 _S1872;
    (&_S1872)->primal_0 = _S1719;
    (&_S1872)->differential_0 = 0.0f;
    DiffPair_float_0 _S1873;
    (&_S1873)->primal_0 = _S1694;
    (&_S1873)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S1872, &_S1873, _s_dOut_8.x);
    DiffPair_float_0 _S1874 = _S1873;
    DiffPair_float_0 _S1875;
    (&_S1875)->primal_0 = 0.0f;
    (&_S1875)->differential_0 = 0.0f;
    DiffPair_float_0 _S1876;
    (&_S1876)->primal_0 = y_13;
    (&_S1876)->differential_0 = 0.0f;
    _d_max_0(&_S1875, &_S1876, _S1872.differential_0);
    DiffPair_float_0 _S1877 = _S1876;
    if(_S1702)
    {
        float _S1878 = a_4 * _S1877.differential_0;
        float _S1879 = _S1703 * _S1877.differential_0;
        DiffPair_float_0 _S1880;
        (&_S1880)->primal_0 = _S1704;
        (&_S1880)->differential_0 = 0.0f;
        DiffPair_float_0 _S1881;
        (&_S1881)->primal_0 = _S1690;
        (&_S1881)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1880, &_S1881, _S1878);
        float _S1882 = _S1880.differential_0 / _S1705;
        float _S1883 = _S1688 * - _S1882;
        float _S1884 = _S1697 * _S1882;
        y_13 = 0.0f;
        _S1703 = _S1879;
        _S1704 = _S1883;
        _S1705 = 0.0f;
        _S1706 = _S1881.differential_0;
        _S1707 = _S1884;
    }
    else
    {
        float _S1885 = - _S1877.differential_0;
        float _S1886 = b_5 * _S1885;
        float _S1887 = _S1706 * _S1885;
        DiffPair_float_0 _S1888;
        (&_S1888)->primal_0 = _S1707;
        (&_S1888)->differential_0 = 0.0f;
        DiffPair_float_0 _S1889;
        (&_S1889)->primal_0 = _S1692;
        (&_S1889)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S1888, &_S1889, _S1886);
        float _S1890 = _S1888.differential_0 / _S1708;
        float _S1891 = - (_S1709 * - _S1890);
        float _S1892 = - (_S1710 * _S1890);
        y_13 = _S1887;
        _S1703 = 0.0f;
        _S1704 = _S1891;
        _S1705 = _S1889.differential_0;
        _S1706 = 0.0f;
        _S1707 = _S1892;
    }
    float _S1893 = (- y_13 + _S1703) / _S1701;
    float _S1894 = _S1700 * - _S1893;
    float _S1895 = _S1699 * _S1893;
    float _S1896 = _S1692 * _S1895;
    float _S1897 = _S1697 * _S1895;
    DiffPair_float_0 _S1898;
    (&_S1898)->primal_0 = _S1690;
    (&_S1898)->differential_0 = 0.0f;
    DiffPair_float_0 _S1899;
    (&_S1899)->primal_0 = _S1692;
    (&_S1899)->differential_0 = 0.0f;
    DiffPair_float_0 _S1900;
    (&_S1900)->primal_0 = _S1697;
    (&_S1900)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S1898, &_S1899, &_S1900, _S1894);
    float _S1901 = - ((_S1896 + _S1900.differential_0 + _S1704) / _S1698);
    DiffPair_float_0 _S1902;
    (&_S1902)->primal_0 = _S1695;
    (&_S1902)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1902, _S1901);
    float _S1903 = - _S1902.differential_0;
    DiffPair_float_0 _S1904;
    (&_S1904)->primal_0 = _S1693;
    (&_S1904)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1904, _S1874.differential_0);
    DiffPair_float_0 _S1905;
    (&_S1905)->primal_0 = _S1634;
    (&_S1905)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1905, _S1904.differential_0);
    DiffPair_float_0 _S1906 = _S1905;
    float _S1907 = _S1897 + _S1899.differential_0 + _S1705;
    DiffPair_float_0 _S1908;
    (&_S1908)->primal_0 = _S1691;
    (&_S1908)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1908, _S1907);
    DiffPair_float_0 _S1909;
    (&_S1909)->primal_0 = _S1633;
    (&_S1909)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1909, _S1908.differential_0);
    DiffPair_float_0 _S1910 = _S1909;
    float _S1911 = _S1898.differential_0 + _S1706;
    DiffPair_float_0 _S1912;
    (&_S1912)->primal_0 = _S1689;
    (&_S1912)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1912, _S1911);
    DiffPair_float_0 _S1913;
    (&_S1913)->primal_0 = _S1632;
    (&_S1913)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1913, _S1912.differential_0);
    DiffPair_float_0 _S1914 = _S1913;
    float3  _S1915 = _S1871 + make_float3 (_S1707, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1916;
    (&_S1916)->primal_0 = _S1684;
    (&_S1916)->differential_0 = _S1614;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1917;
    (&_S1917)->primal_0 = _S1685;
    (&_S1917)->differential_0 = _S1614;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1918;
    (&_S1918)->primal_0 = _S1686;
    (&_S1918)->differential_0 = _S1614;
    s_bwd_prop_clamp_1(&_S1916, &_S1917, &_S1918, _S1915);
    float _S1919 = - _S1916.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S1916.differential_0.x + _S1919, _S1916.differential_0.y + _S1919, _S1916.differential_0.z);
    float3  _S1920 = _S1678 * s_diff_rgi_out_T_0;
    float _S1921 = (_S1920.x + _S1920.y + _S1920.z) / _S1681;
    float _S1922 = _S1679 * _S1921;
    float3  _S1923 = _S1680 * s_diff_rgi_out_T_0 + make_float3 (0.0f, 0.0f, intensity_2 * - _S1921);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1924;
    (&_S1924)->primal_0 = H_4;
    (&_S1924)->differential_0 = _S1615;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1925;
    (&_S1925)->primal_0 = rgi_in_0;
    (&_S1925)->differential_0 = _S1614;
    s_bwd_prop_mul_0(&_S1924, &_S1925, _S1923);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1926 = _S1924;
    float _S1927 = _S1922 + _S1925.differential_0.z;
    float _S1928 = _S1925.differential_0.y + _S1927;
    float _S1929 = _S1925.differential_0.x + _S1927;
    float3  _S1930 = make_float3 (_S1929, _S1928, _S1927);
    if(_S1671)
    {
        Matrix<float, 3, 3>  _S1931 = _S1670 * _S1926.differential_0;
        Matrix<float, 3, 3>  _S1932 = _S1672 * _S1926.differential_0;
        _S1673 = - ((_S1931.rows[int(0)].x + _S1931.rows[int(0)].y + _S1931.rows[int(0)].z + _S1931.rows[int(1)].x + _S1931.rows[int(1)].y + _S1931.rows[int(1)].z + _S1931.rows[int(2)].x + _S1931.rows[int(2)].y + _S1931.rows[int(2)].z) / _S1673);
        H_4 = _S1932;
    }
    else
    {
        _S1673 = 0.0f;
        H_4 = _S1926.differential_0;
    }
    DiffPair_float_0 _S1933;
    (&_S1933)->primal_0 = _S1670.rows[int(2)].z;
    (&_S1933)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S1933, 0.0f);
    float _S1934 = _S1933.differential_0 + _S1673;
    float3  _S1935 = _S1614;
    *&((&_S1935)->z) = _S1934;
    Matrix<float, 3, 3>  _S1936 = _S1615;
    _S1936[int(2)] = _S1935;
    Matrix<float, 3, 3>  _S1937 = H_4 + _S1936;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1938;
    (&_S1938)->primal_0 = _S1669;
    (&_S1938)->differential_0 = _S1615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1939;
    (&_S1939)->primal_0 = S_inv_0;
    (&_S1939)->differential_0 = _S1615;
    s_bwd_prop_mul_1(&_S1938, &_S1939, _S1937);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1940;
    (&_S1940)->primal_0 = T_2;
    (&_S1940)->differential_0 = _S1615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1941;
    (&_S1941)->primal_0 = D_0;
    (&_S1941)->differential_0 = _S1615;
    s_bwd_prop_mul_1(&_S1940, &_S1941, _S1938.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1942 = _S1940;
    float3  _S1943 = make_float3 (_S1941.differential_0.rows[int(0)].x, _S1941.differential_0.rows[int(1)].y, _S1941.differential_0.rows[int(2)].z);
    float3  _S1944;
    if(_S1664)
    {
        if(_S1666)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1945;
            (&_S1945)->primal_0 = r1_2;
            (&_S1945)->differential_0 = _S1614;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1946;
            (&_S1946)->primal_0 = r2_14;
            (&_S1946)->differential_0 = _S1614;
            s_bwd_prop_cross_0(&_S1945, &_S1946, _S1943);
            _S1652 = _S1614;
            lambda_v_6 = _S1946.differential_0;
            _S1944 = _S1945.differential_0;
        }
        else
        {
            _S1652 = _S1943;
            lambda_v_6 = _S1614;
            _S1944 = _S1614;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1947;
        (&_S1947)->primal_0 = _S1665;
        (&_S1947)->differential_0 = _S1614;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1948;
        (&_S1948)->primal_0 = _S1665;
        (&_S1948)->differential_0 = _S1614;
        s_bwd_prop_dot_0(&_S1947, &_S1948, 0.0f);
        float3  _S1949 = _S1948.differential_0 + _S1947.differential_0 + _S1652;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1950;
        (&_S1950)->primal_0 = r0_2;
        (&_S1950)->differential_0 = _S1614;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1951;
        (&_S1951)->primal_0 = r2_14;
        (&_S1951)->differential_0 = _S1614;
        s_bwd_prop_cross_0(&_S1950, &_S1951, _S1949);
        float3  _S1952 = _S1951.differential_0 + lambda_v_6;
        _S1652 = _S1614;
        lambda_v_6 = _S1952;
        _S1665 = _S1944;
        _S1944 = _S1950.differential_0;
    }
    else
    {
        _S1652 = _S1943;
        lambda_v_6 = _S1614;
        _S1665 = _S1614;
        _S1944 = _S1614;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1953;
    (&_S1953)->primal_0 = _S1663;
    (&_S1953)->differential_0 = _S1614;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1954;
    (&_S1954)->primal_0 = _S1663;
    (&_S1954)->differential_0 = _S1614;
    s_bwd_prop_dot_0(&_S1953, &_S1954, 0.0f);
    float3  _S1955 = _S1954.differential_0 + _S1953.differential_0 + _S1652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1956;
    (&_S1956)->primal_0 = r0_2;
    (&_S1956)->differential_0 = _S1614;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1957;
    (&_S1957)->primal_0 = r1_2;
    (&_S1957)->differential_0 = _S1614;
    s_bwd_prop_cross_0(&_S1956, &_S1957, _S1955);
    float3  _S1958 = _S1614;
    *&((&_S1958)->z) = lambda_v_6.z;
    *&((&_S1958)->y) = lambda_v_6.y;
    *&((&_S1958)->x) = lambda_v_6.x;
    float3  _S1959 = _S1957.differential_0 + _S1665;
    float3  _S1960 = _S1614;
    *&((&_S1960)->z) = _S1959.z;
    *&((&_S1960)->y) = _S1959.y;
    *&((&_S1960)->x) = _S1959.x;
    float3  _S1961 = _S1956.differential_0 + _S1944;
    float3  _S1962 = _S1614;
    *&((&_S1962)->z) = _S1961.z;
    *&((&_S1962)->y) = _S1961.y;
    *&((&_S1962)->x) = _S1961.x;
    Matrix<float, 3, 3>  _S1963 = _S1615;
    _S1963[int(2)] = _S1958;
    _S1963[int(1)] = _S1960;
    _S1963[int(0)] = _S1962;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1964;
    (&_S1964)->primal_0 = skew_0;
    (&_S1964)->differential_0 = _S1615;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1965;
    (&_S1965)->primal_0 = T_2;
    (&_S1965)->differential_0 = _S1615;
    s_bwd_prop_mul_1(&_S1964, &_S1965, _S1963);
    Matrix<float, 3, 3>  _S1966 = _S1965.differential_0 + _S1942.differential_0;
    float2  _S1967 = make_float2 (_S1964.differential_0.rows[int(2)].y + - _S1964.differential_0.rows[int(1)].z, _S1964.differential_0.rows[int(0)].z + - _S1964.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1968 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1969;
    (&_S1969)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1969)->differential_0 = _S1968;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1970;
    (&_S1970)->primal_0 = _S1655.color_params_1.n_0;
    (&_S1970)->differential_0 = _S1618;
    s_bwd_prop_mul_2(&_S1969, &_S1970, _S1967);
    float2  _S1971 = make_float2 (_S1966.rows[int(0)].z, _S1966.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1972;
    (&_S1972)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1972)->differential_0 = _S1968;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1973;
    (&_S1973)->primal_0 = _S1655.color_params_1.g_0;
    (&_S1973)->differential_0 = _S1618;
    s_bwd_prop_mul_2(&_S1972, &_S1973, _S1971);
    float2  _S1974 = make_float2 (_S1966.rows[int(0)].y, _S1966.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1975;
    (&_S1975)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1975)->differential_0 = _S1968;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1976;
    (&_S1976)->primal_0 = _S1655.color_params_1.r_0;
    (&_S1976)->differential_0 = _S1618;
    s_bwd_prop_mul_2(&_S1975, &_S1976, _S1974);
    float2  _S1977 = make_float2 (_S1966.rows[int(0)].x, _S1966.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1978;
    (&_S1978)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1978)->differential_0 = _S1968;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1979;
    (&_S1979)->primal_0 = _S1655.color_params_1.b_0;
    (&_S1979)->differential_0 = _S1618;
    s_bwd_prop_mul_2(&_S1978, &_S1979, _S1977);
    ColorPPISPParams_0 _S1980 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1980)->n_0 = _S1970.differential_0;
    (&_S1980)->g_0 = _S1973.differential_0;
    (&_S1980)->r_0 = _S1976.differential_0;
    (&_S1980)->b_0 = _S1979.differential_0;
    _S1652 = _S1930;
    *&((&_S1652)->z) = 0.0f;
    float _S1981 = rgb_out_4.z * _S1927;
    float _S1982 = _S1654 * _S1927;
    DiffPair_float_0 _S1983;
    (&_S1983)->primal_0 = falloff_2;
    (&_S1983)->differential_0 = 0.0f;
    DiffPair_float_0 _S1984;
    (&_S1984)->primal_0 = 0.0f;
    (&_S1984)->differential_0 = 0.0f;
    DiffPair_float_0 _S1985;
    (&_S1985)->primal_0 = 1.0f;
    (&_S1985)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1983, &_S1984, &_S1985, _S1981);
    float _S1986 = r2_13 * _S1983.differential_0;
    float _S1987 = r4_8 * _S1983.differential_0;
    float s_diff_r6_T_0 = _S1631 * _S1983.differential_0;
    float _S1988 = r6_2 * _S1983.differential_0;
    float _S1989 = r2_13 * (_S1630 * _S1983.differential_0 + r2_13 * s_diff_r6_T_0);
    float _S1990 = _S1629 * _S1983.differential_0 + r4_8 * s_diff_r6_T_0 + _S1989 + _S1989;
    float _S1991 = dy_8 * _S1990;
    float _S1992 = dx_10 * _S1990;
    float _S1993 = - (_S1991 + _S1991);
    float _S1994 = - (_S1992 + _S1992);
    *&((&_S1652)->y) = 0.0f;
    float _S1995 = rgb_out_4.y * _S1928;
    float _S1996 = _S1653 * _S1928;
    DiffPair_float_0 _S1997;
    (&_S1997)->primal_0 = falloff_1;
    (&_S1997)->differential_0 = 0.0f;
    DiffPair_float_0 _S1998;
    (&_S1998)->primal_0 = 0.0f;
    (&_S1998)->differential_0 = 0.0f;
    DiffPair_float_0 _S1999;
    (&_S1999)->primal_0 = 1.0f;
    (&_S1999)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1997, &_S1998, &_S1999, _S1995);
    float _S2000 = r2_12 * _S1997.differential_0;
    float _S2001 = r4_7 * _S1997.differential_0;
    float s_diff_r6_T_1 = _S1628 * _S1997.differential_0;
    float _S2002 = r6_1 * _S1997.differential_0;
    float _S2003 = r2_12 * (_S1627 * _S1997.differential_0 + r2_12 * s_diff_r6_T_1);
    float _S2004 = _S1626 * _S1997.differential_0 + r4_7 * s_diff_r6_T_1 + _S2003 + _S2003;
    float _S2005 = dy_7 * _S2004;
    float _S2006 = dx_9 * _S2004;
    float _S2007 = - (_S2005 + _S2005);
    float _S2008 = - (_S2006 + _S2006);
    *&((&_S1652)->x) = 0.0f;
    float _S2009 = rgb_out_4.x * _S1929;
    float _S2010 = _S1650 * _S1929;
    DiffPair_float_0 _S2011;
    (&_S2011)->primal_0 = falloff_0;
    (&_S2011)->differential_0 = 0.0f;
    DiffPair_float_0 _S2012;
    (&_S2012)->primal_0 = 0.0f;
    (&_S2012)->differential_0 = 0.0f;
    DiffPair_float_0 _S2013;
    (&_S2013)->primal_0 = 1.0f;
    (&_S2013)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2011, &_S2012, &_S2013, _S2009);
    float _S2014 = r2_11 * _S2011.differential_0;
    float _S2015 = r4_6 * _S2011.differential_0;
    float s_diff_r6_T_2 = _S1625 * _S2011.differential_0;
    float _S2016 = r6_0 * _S2011.differential_0;
    float _S2017 = r2_11 * (_S1624 * _S2011.differential_0 + r2_11 * s_diff_r6_T_2);
    float _S2018 = _S1623 * _S2011.differential_0 + r4_6 * s_diff_r6_T_2 + _S2017 + _S2017;
    float _S2019 = dy_6 * _S2018;
    float _S2020 = dx_8 * _S2018;
    float _S2021 = - (_S2019 + _S2019);
    float _S2022 = - (_S2020 + _S2020);
    float3  _S2023 = _S1614;
    *&((&_S2023)->z) = _S1982;
    *&((&_S2023)->y) = _S1996;
    *&((&_S2023)->x) = _S2010;
    float3  _S2024 = _S1652 + _S2023;
    float3  _S2025 = _S1613.primal_0 * _S2024;
    float3  _S2026 = _S1646 * _S2024;
    float _S2027 = _S2025.x + _S2025.y + _S2025.z;
    DiffPair_float_0 _S2028;
    (&_S2028)->primal_0 = _S1644.exposure_1;
    (&_S2028)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2028, _S2027);
    PPISPParams_0 _S2029 = PPISPParams_x24_syn_dzero_0();
    (&_S2029)->color_params_1 = _S1980;
    (&_S2029)->exposure_1 = _S2028.differential_0;
    _S1622 = _S2029;
    (&(&_S1622)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2030 = _S2029.crf_params_1[int(2)].center_0 + _S1815;
    (&(&_S1622)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2031 = _S2029.crf_params_1[int(2)].gamma_0 + _S1818.differential_0;
    (&(&_S1622)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2032 = _S2029.crf_params_1[int(2)].shoulder_0 + _S1822.differential_0;
    (&(&_S1622)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2033 = _S2029.crf_params_1[int(2)].toe_0 + _S1826.differential_0;
    (&(&_S1622)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2034 = _S2029.crf_params_1[int(1)].center_0 + _S1859;
    (&(&_S1622)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2035 = _S2029.crf_params_1[int(1)].gamma_0 + _S1862.differential_0;
    (&(&_S1622)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2036 = _S2029.crf_params_1[int(1)].shoulder_0 + _S1866.differential_0;
    (&(&_S1622)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2037 = _S2029.crf_params_1[int(1)].toe_0 + _S1870.differential_0;
    (&(&_S1622)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2038 = _S2029.crf_params_1[int(0)].center_0 + _S1903;
    (&(&_S1622)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2039 = _S2029.crf_params_1[int(0)].gamma_0 + _S1906.differential_0;
    (&(&_S1622)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2040 = _S2029.crf_params_1[int(0)].shoulder_0 + _S1910.differential_0;
    (&(&_S1622)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2041 = _S2029.crf_params_1[int(0)].toe_0 + _S1914.differential_0;
    *&((&(&(&_S1622)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1622)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1622)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2042 = _S1988 + _S2029.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1622)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2043 = _S1987 + _S2029.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1622)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2044 = _S1986 + _S2029.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1622)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2045 = _S1993 + _S2029.vignette_params_1[int(2)].cy_0;
    (&(&_S1622)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2046 = _S1994 + _S2029.vignette_params_1[int(2)].cx_0;
    (&(&_S1622)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2047 = _S2002 + _S2029.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1622)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2048 = _S2001 + _S2029.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1622)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2049 = _S2000 + _S2029.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1622)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2050 = _S2007 + _S2029.vignette_params_1[int(1)].cy_0;
    (&(&_S1622)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2051 = _S2008 + _S2029.vignette_params_1[int(1)].cx_0;
    (&(&_S1622)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2052 = _S2016 + _S2029.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1622)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2053 = _S2015 + _S2029.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1622)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2054 = _S2014 + _S2029.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1622)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2055 = _S2021 + _S2029.vignette_params_1[int(0)].cy_0;
    (&(&_S1622)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2056 = _S2022 + _S2029.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2057;
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
    _S2057[int(22)] = 0.0f;
    _S2057[int(23)] = 0.0f;
    _S2057[int(24)] = 0.0f;
    _S2057[int(25)] = 0.0f;
    _S2057[int(26)] = 0.0f;
    _S2057[int(27)] = 0.0f;
    _S2057[int(28)] = 0.0f;
    _S2057[int(29)] = 0.0f;
    _S2057[int(30)] = 0.0f;
    _S2057[int(31)] = 0.0f;
    _S2057[int(32)] = 0.0f;
    _S2057[int(33)] = 0.0f;
    _S2057[int(34)] = 0.0f;
    _S2057[int(35)] = 0.0f;
    _S2057[int(8)] = _S2049;
    _S2057[int(16)] = _S2029.color_params_1.b_0.x;
    _S2057[int(15)] = _S2042;
    _S2057[int(14)] = _S2043;
    _S2057[int(13)] = _S2044;
    _S2057[int(12)] = _S2045;
    _S2057[int(11)] = _S2046;
    _S2057[int(10)] = _S2047;
    _S2057[int(9)] = _S2048;
    _S2057[int(17)] = _S2029.color_params_1.b_0.y;
    _S2057[int(7)] = _S2050;
    _S2057[int(6)] = _S2051;
    _S2057[int(5)] = _S2052;
    _S2057[int(4)] = _S2053;
    _S2057[int(3)] = _S2054;
    _S2057[int(2)] = _S2055;
    _S2057[int(1)] = _S2056;
    _S2057[int(0)] = _S1622.exposure_1;
    _S2057[int(26)] = _S2039;
    _S2057[int(34)] = _S2031;
    _S2057[int(33)] = _S2032;
    _S2057[int(32)] = _S2033;
    _S2057[int(31)] = _S2034;
    _S2057[int(30)] = _S2035;
    _S2057[int(29)] = _S2036;
    _S2057[int(28)] = _S2037;
    _S2057[int(27)] = _S2038;
    _S2057[int(35)] = _S2030;
    _S2057[int(25)] = _S2040;
    _S2057[int(24)] = _S2041;
    _S2057[int(23)] = _S2029.color_params_1.n_0.y;
    _S2057[int(22)] = _S2029.color_params_1.n_0.x;
    _S2057[int(21)] = _S2029.color_params_1.g_0.y;
    _S2057[int(20)] = _S2029.color_params_1.g_0.x;
    _S2057[int(19)] = _S2029.color_params_1.r_0.y;
    _S2057[int(18)] = _S2029.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S2057;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S2026;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2058, float2  _S2059, float2  _S2060, float2  _S2061, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2062, float3  _S2063)
{
    s_bwd_prop_apply_ppisp_0(_S2058, _S2059, _S2060, _S2061, _S2062, _S2063);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S2064 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S2064;
    FixedArray<float, 36>  _S2065 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_2;
    (&dp_params_0)->differential_0 = _S2065;
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

inline __device__ void s_bwd_prop_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_1, float2  pix_coord_4, float2  image_center_4, float2  img_size_4, DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_1, float3  _s_dOut_9)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2066 = *dprgb_in_1;
    float3  _S2067 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S2068 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S2069 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2070 = {
        _S2069, _S2069, _S2069
    };
    float2  _S2071 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2072 = { _S2071, _S2071, _S2071, _S2071 };
    RQSCRFPPISPChannelParams_0 _S2073 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2074 = {
        _S2073, _S2073, _S2073
    };
    PPISPParamsRQS_0 _S2075;
    (&_S2075)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S2075)->vignette_params_0 = _S2070;
    (&_S2075)->color_params_0 = _S2072;
    (&_S2075)->crf_params_0 = _S2074;
    (&(&_S2075)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S2075)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S2076 = dpparams_1->primal_0[int(3)];
    (&(&_S2075)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S2077 = dpparams_1->primal_0[int(4)];
    (&(&_S2075)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S2078 = dpparams_1->primal_0[int(5)];
    (&(&_S2075)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S2075)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S2075)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S2079 = dpparams_1->primal_0[int(8)];
    (&(&_S2075)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S2080 = dpparams_1->primal_0[int(9)];
    (&(&_S2075)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S2081 = dpparams_1->primal_0[int(10)];
    (&(&_S2075)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S2075)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S2075)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S2082 = dpparams_1->primal_0[int(13)];
    (&(&_S2075)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S2083 = dpparams_1->primal_0[int(14)];
    (&(&_S2075)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S2084 = dpparams_1->primal_0[int(15)];
    (&(&_S2075)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S2075)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S2075)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S2075)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S2075)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S2075)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S2075)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S2075)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S2075)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S2085 = dpparams_1->primal_0[int(24)];
    (&(&_S2075)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S2086 = dpparams_1->primal_0[int(25)];
    (&(&_S2075)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S2087 = dpparams_1->primal_0[int(26)];
    (&(&_S2075)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S2088 = dpparams_1->primal_0[int(27)];
    (&(&_S2075)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S2089 = dpparams_1->primal_0[int(28)];
    (&(&_S2075)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S2090 = dpparams_1->primal_0[int(29)];
    (&(&_S2075)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S2091 = dpparams_1->primal_0[int(30)];
    (&(&_S2075)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S2092 = dpparams_1->primal_0[int(31)];
    (&(&_S2075)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S2093 = dpparams_1->primal_0[int(32)];
    (&(&_S2075)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S2094 = dpparams_1->primal_0[int(33)];
    (&(&_S2075)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S2095 = dpparams_1->primal_0[int(34)];
    (&(&_S2075)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S2096 = dpparams_1->primal_0[int(35)];
    (&(&_S2075)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S2097 = dpparams_1->primal_0[int(36)];
    (&(&_S2075)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S2098 = dpparams_1->primal_0[int(37)];
    (&(&_S2075)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S2099 = dpparams_1->primal_0[int(38)];
    (&(&_S2075)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S2100 = _S2075;
    float _S2101 = s_primal_ctx_exp2_0(_S2075.exposure_0);
    float3  _S2102 = make_float3 (_S2101);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S2101);
    float _S2103 = (F32_max((img_size_4.x), (img_size_4.y)));
    float _S2104 = (pix_coord_4.x - image_center_4.x) / _S2103;
    float _S2105 = (pix_coord_4.y - image_center_4.y) / _S2103;
    float dx_11 = _S2104 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S2105 - dpparams_1->primal_0[int(2)];
    float r2_15 = dx_11 * dx_11 + dy_9 * dy_9;
    float r4_9 = r2_15 * r2_15;
    float r6_3 = r4_9 * r2_15;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_15 + 1.0f;
    float _S2106 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S2107 = rgb_out_5.x * _S2106;
    float3  _S2108 = rgb_out_5;
    *&((&_S2108)->x) = _S2107;
    float dx_12 = _S2104 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S2105 - dpparams_1->primal_0[int(7)];
    float r2_16 = dx_12 * dx_12 + dy_10 * dy_10;
    float r4_10 = r2_16 * r2_16;
    float r6_4 = r4_10 * r2_16;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_16 + 1.0f;
    float _S2109 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S2108)->y) = rgb_out_5.y * _S2109;
    float dx_13 = _S2104 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S2105 - dpparams_1->primal_0[int(12)];
    float r2_17 = dx_13 * dx_13 + dy_11 * dy_11;
    float r4_11 = r2_17 * r2_17;
    float r6_5 = r4_11 * r2_17;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_17 + 1.0f;
    float _S2110 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S2108)->z) = rgb_out_5.z * _S2110;
    PPISPParamsRQS_0 _S2111 = _S2075;
    float2  _S2112 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S2075.color_params_0.b_0);
    float2  _S2113 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S2075.color_params_0.r_0);
    float2  _S2114 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S2075.color_params_0.g_0);
    float2  _S2115 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S2075.color_params_0.n_0);
    float _S2116 = 0.3333333432674408f + _S2115.x;
    float _S2117 = 0.3333333432674408f + _S2115.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S2112.x, 1.0f + _S2113.x, _S2114.x, _S2112.y, _S2113.y, 1.0f + _S2114.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S2117, 1.0f, 0.0f, - _S2116, - _S2117, _S2116, 0.0f);
    Matrix<float, 3, 3>  _S2118 = s_primal_ctx_mul_1(skew_1, T_3);
    float3  r0_3 = make_float3 (_S2118.rows[int(0)].x, _S2118.rows[int(0)].y, _S2118.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S2118.rows[int(1)].x, _S2118.rows[int(1)].y, _S2118.rows[int(1)].z);
    float3  r2_18 = make_float3 (_S2118.rows[int(2)].x, _S2118.rows[int(2)].y, _S2118.rows[int(2)].z);
    float3  _S2119 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S2120 = (s_primal_ctx_dot_0(_S2119, _S2119)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S2121;
    bool _S2122;
    if(_S2120)
    {
        float3  _S2123 = s_primal_ctx_cross_0(r0_3, r2_18);
        bool _S2124 = (s_primal_ctx_dot_0(_S2123, _S2123)) < 9.99999968265522539e-21f;
        if(_S2124)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_18);
        }
        else
        {
            lambda_v_7 = _S2123;
        }
        _S2122 = _S2124;
        _S2121 = _S2123;
    }
    else
    {
        lambda_v_7 = _S2119;
        _S2122 = false;
        _S2121 = _S2067;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S2125 = s_primal_ctx_mul_1(T_3, D_1);
    Matrix<float, 3, 3>  _S2126 = s_primal_ctx_mul_1(_S2125, S_inv_1);
    bool _S2127 = (s_primal_ctx_abs_0(_S2126.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S2128;
    float _S2129;
    if(_S2127)
    {
        float inv_s_1 = 1.0f / _S2126.rows[int(2)].z;
        Matrix<float, 3, 3>  _S2130 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S2131 = _S2126.rows[int(2)].z * _S2126.rows[int(2)].z;
        H_5 = _S2126 * makeMatrix<float, 3, 3> (inv_s_1);
        _S2128 = _S2130;
        _S2129 = _S2131;
    }
    else
    {
        H_5 = _S2126;
        _S2128 = _S2068;
        _S2129 = 0.0f;
    }
    float _S2132 = _S2108.x;
    float _S2133 = _S2108.y;
    float intensity_3 = _S2132 + _S2133 + _S2108.z;
    float3  rgi_in_1 = make_float3 (_S2132, _S2133, intensity_3);
    float3  _S2134 = s_primal_ctx_mul_2(H_5, rgi_in_1);
    float _S2135 = _S2134.z + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S2135;
    float3  _S2136 = make_float3 (norm_factor_1);
    float _S2137 = _S2135 * _S2135;
    float3  rgi_out_5 = _S2134 * make_float3 (norm_factor_1);
    float _S2138 = rgi_out_5.x;
    float _S2139 = rgi_out_5.y;
    float3  _S2140 = make_float3 (_S2138, _S2139, rgi_out_5.z - _S2138 - _S2139);
    float3  _S2141 = make_float3 (0.0f);
    float3  _S2142 = make_float3 (1.0f);
    float3  _S2143 = s_primal_ctx_clamp_1(_S2140, _S2141, _S2142);
    float _S2144 = _S2143.x;
    float _S2145 = s_primal_ctx_exp_0(_S2085);
    float _S2146 = s_primal_ctx_exp_0(_S2086);
    float _S2147 = - _S2087;
    float _S2148 = 1.0f + s_primal_ctx_exp_0(_S2147);
    float x0_4 = 1.0f / _S2148;
    float _S2149 = _S2148 * _S2148;
    float _S2150 = - _S2088;
    float _S2151 = 1.0f + s_primal_ctx_exp_0(_S2150);
    float y0_4 = 1.0f / _S2151;
    float _S2152 = _S2151 * _S2151;
    float _S2153 = s_primal_ctx_exp_0(_S2089);
    bool _S2154 = _S2144 < x0_4;
    float _S2155;
    float _S2156;
    float _S2157;
    float _S2158;
    float _S2159;
    float _S2160;
    float _S2161;
    float _S2162;
    float _S2163;
    float _S2164;
    float _S2165;
    float _S2166;
    float _S2167;
    float _S2168;
    float _S2169;
    float _S2170;
    float _S2171;
    float _S2172;
    float _S2173;
    float _S2174;
    float _S2175;
    float _S2176;
    float _S2177;
    float _S2178;
    float _S2179;
    float _S2180;
    float _S2181;
    if(_S2154)
    {
        float s0_3 = y0_4 / x0_4;
        float _S2182 = x0_4 * x0_4;
        float t0_3 = _S2144 / x0_4;
        float _S2183 = s0_3 * t0_3;
        float _S2184 = _S2145 * t0_3;
        float _S2185 = 1.0f - t0_3;
        float _S2186 = _S2183 * t0_3 + _S2184 * _S2185;
        float _S2187 = y0_4 * _S2186;
        float _S2188 = _S2145 + _S2153 - 2.0f * s0_3;
        float _S2189 = _S2188 * t0_3;
        float _S2190 = s0_3 + _S2189 * _S2185;
        _S2155 = _S2190 * _S2190;
        _S2156 = _S2187;
        _S2157 = _S2190;
        _S2158 = _S2189;
        _S2159 = _S2185;
        _S2160 = _S2188;
        _S2161 = t0_3;
        _S2162 = _S2186;
        _S2163 = _S2184;
        _S2164 = _S2183;
        _S2165 = s0_3;
        _S2166 = _S2182;
        _S2167 = 0.0f;
        _S2168 = 0.0f;
        _S2169 = 0.0f;
        _S2170 = 0.0f;
        _S2171 = 0.0f;
        _S2172 = 0.0f;
        _S2173 = 0.0f;
        _S2174 = 0.0f;
        _S2175 = 0.0f;
        _S2176 = 0.0f;
        _S2177 = 0.0f;
        _S2178 = 0.0f;
        _S2179 = 0.0f;
        _S2180 = 0.0f;
        _S2181 = 0.0f;
    }
    else
    {
        float _S2191 = 1.0f - y0_4;
        float _S2192 = 1.0f - x0_4;
        float s1_9 = _S2191 / _S2192;
        float _S2193 = _S2192 * _S2192;
        float _S2194 = _S2144 - x0_4;
        float t1_3 = _S2194 / _S2192;
        float _S2195 = s1_9 * t1_3;
        float _S2196 = _S2153 * t1_3;
        float _S2197 = 1.0f - t1_3;
        float _S2198 = _S2195 * t1_3 + _S2196 * _S2197;
        float _S2199 = _S2191 * _S2198;
        float _S2200 = _S2153 + _S2146 - 2.0f * s1_9;
        float _S2201 = _S2200 * t1_3;
        float _S2202 = s1_9 + _S2201 * _S2197;
        float _S2203 = _S2202 * _S2202;
        _S2155 = 0.0f;
        _S2156 = 0.0f;
        _S2157 = 0.0f;
        _S2158 = 0.0f;
        _S2159 = 0.0f;
        _S2160 = 0.0f;
        _S2161 = 0.0f;
        _S2162 = 0.0f;
        _S2163 = 0.0f;
        _S2164 = 0.0f;
        _S2165 = 0.0f;
        _S2166 = 0.0f;
        _S2167 = _S2203;
        _S2168 = _S2199;
        _S2169 = _S2202;
        _S2170 = _S2201;
        _S2171 = _S2197;
        _S2172 = _S2200;
        _S2173 = t1_3;
        _S2174 = _S2191;
        _S2175 = _S2198;
        _S2176 = _S2196;
        _S2177 = _S2195;
        _S2178 = s1_9;
        _S2179 = _S2193;
        _S2180 = _S2194;
        _S2181 = _S2192;
    }
    float _S2204 = _S2143.y;
    float _S2205 = s_primal_ctx_exp_0(_S2090);
    float _S2206 = s_primal_ctx_exp_0(_S2091);
    float _S2207 = - _S2092;
    float _S2208 = 1.0f + s_primal_ctx_exp_0(_S2207);
    float x0_5 = 1.0f / _S2208;
    float _S2209 = _S2208 * _S2208;
    float _S2210 = - _S2093;
    float _S2211 = 1.0f + s_primal_ctx_exp_0(_S2210);
    float y0_5 = 1.0f / _S2211;
    float _S2212 = _S2211 * _S2211;
    float _S2213 = s_primal_ctx_exp_0(_S2094);
    bool _S2214 = _S2204 < x0_5;
    float _S2215;
    float _S2216;
    float _S2217;
    float _S2218;
    float _S2219;
    float _S2220;
    float _S2221;
    float _S2222;
    float _S2223;
    float _S2224;
    float _S2225;
    float _S2226;
    float _S2227;
    float _S2228;
    float _S2229;
    float _S2230;
    float _S2231;
    float _S2232;
    float _S2233;
    float _S2234;
    float _S2235;
    float _S2236;
    float _S2237;
    float _S2238;
    float _S2239;
    float _S2240;
    float _S2241;
    if(_S2214)
    {
        float s0_4 = y0_5 / x0_5;
        float _S2242 = x0_5 * x0_5;
        float t0_4 = _S2204 / x0_5;
        float _S2243 = s0_4 * t0_4;
        float _S2244 = _S2205 * t0_4;
        float _S2245 = 1.0f - t0_4;
        float _S2246 = _S2243 * t0_4 + _S2244 * _S2245;
        float _S2247 = y0_5 * _S2246;
        float _S2248 = _S2205 + _S2213 - 2.0f * s0_4;
        float _S2249 = _S2248 * t0_4;
        float _S2250 = s0_4 + _S2249 * _S2245;
        _S2215 = _S2250 * _S2250;
        _S2216 = _S2247;
        _S2217 = _S2250;
        _S2218 = _S2249;
        _S2219 = _S2245;
        _S2220 = _S2248;
        _S2221 = t0_4;
        _S2222 = _S2246;
        _S2223 = _S2244;
        _S2224 = _S2243;
        _S2225 = s0_4;
        _S2226 = _S2242;
        _S2227 = 0.0f;
        _S2228 = 0.0f;
        _S2229 = 0.0f;
        _S2230 = 0.0f;
        _S2231 = 0.0f;
        _S2232 = 0.0f;
        _S2233 = 0.0f;
        _S2234 = 0.0f;
        _S2235 = 0.0f;
        _S2236 = 0.0f;
        _S2237 = 0.0f;
        _S2238 = 0.0f;
        _S2239 = 0.0f;
        _S2240 = 0.0f;
        _S2241 = 0.0f;
    }
    else
    {
        float _S2251 = 1.0f - y0_5;
        float _S2252 = 1.0f - x0_5;
        float s1_10 = _S2251 / _S2252;
        float _S2253 = _S2252 * _S2252;
        float _S2254 = _S2204 - x0_5;
        float t1_4 = _S2254 / _S2252;
        float _S2255 = s1_10 * t1_4;
        float _S2256 = _S2213 * t1_4;
        float _S2257 = 1.0f - t1_4;
        float _S2258 = _S2255 * t1_4 + _S2256 * _S2257;
        float _S2259 = _S2251 * _S2258;
        float _S2260 = _S2213 + _S2206 - 2.0f * s1_10;
        float _S2261 = _S2260 * t1_4;
        float _S2262 = s1_10 + _S2261 * _S2257;
        float _S2263 = _S2262 * _S2262;
        _S2215 = 0.0f;
        _S2216 = 0.0f;
        _S2217 = 0.0f;
        _S2218 = 0.0f;
        _S2219 = 0.0f;
        _S2220 = 0.0f;
        _S2221 = 0.0f;
        _S2222 = 0.0f;
        _S2223 = 0.0f;
        _S2224 = 0.0f;
        _S2225 = 0.0f;
        _S2226 = 0.0f;
        _S2227 = _S2263;
        _S2228 = _S2259;
        _S2229 = _S2262;
        _S2230 = _S2261;
        _S2231 = _S2257;
        _S2232 = _S2260;
        _S2233 = t1_4;
        _S2234 = _S2251;
        _S2235 = _S2258;
        _S2236 = _S2256;
        _S2237 = _S2255;
        _S2238 = s1_10;
        _S2239 = _S2253;
        _S2240 = _S2254;
        _S2241 = _S2252;
    }
    float _S2264 = _S2143.z;
    float _S2265 = s_primal_ctx_exp_0(_S2095);
    float _S2266 = s_primal_ctx_exp_0(_S2096);
    float _S2267 = - _S2097;
    float _S2268 = 1.0f + s_primal_ctx_exp_0(_S2267);
    float x0_6 = 1.0f / _S2268;
    float _S2269 = _S2268 * _S2268;
    float _S2270 = - _S2098;
    float _S2271 = 1.0f + s_primal_ctx_exp_0(_S2270);
    float y0_6 = 1.0f / _S2271;
    float _S2272 = _S2271 * _S2271;
    float _S2273 = s_primal_ctx_exp_0(_S2099);
    bool _S2274 = _S2264 < x0_6;
    float _S2275;
    float _S2276;
    float _S2277;
    float _S2278;
    float _S2279;
    float _S2280;
    float _S2281;
    float _S2282;
    float _S2283;
    float _S2284;
    float _S2285;
    float _S2286;
    float _S2287;
    float _S2288;
    float _S2289;
    float _S2290;
    float _S2291;
    float _S2292;
    float _S2293;
    float _S2294;
    float _S2295;
    float _S2296;
    float _S2297;
    float _S2298;
    float _S2299;
    float _S2300;
    float _S2301;
    if(_S2274)
    {
        float s0_5 = y0_6 / x0_6;
        float _S2302 = x0_6 * x0_6;
        float t0_5 = _S2264 / x0_6;
        float _S2303 = s0_5 * t0_5;
        float _S2304 = _S2265 * t0_5;
        float _S2305 = 1.0f - t0_5;
        float _S2306 = _S2303 * t0_5 + _S2304 * _S2305;
        float _S2307 = y0_6 * _S2306;
        float _S2308 = _S2265 + _S2273 - 2.0f * s0_5;
        float _S2309 = _S2308 * t0_5;
        float _S2310 = s0_5 + _S2309 * _S2305;
        _S2275 = _S2310 * _S2310;
        _S2276 = _S2307;
        _S2277 = _S2310;
        _S2278 = _S2309;
        _S2279 = _S2305;
        _S2280 = _S2308;
        _S2281 = t0_5;
        _S2282 = _S2306;
        _S2283 = _S2304;
        _S2284 = _S2303;
        _S2285 = s0_5;
        _S2286 = _S2302;
        _S2287 = 0.0f;
        _S2288 = 0.0f;
        _S2289 = 0.0f;
        _S2290 = 0.0f;
        _S2291 = 0.0f;
        _S2292 = 0.0f;
        _S2293 = 0.0f;
        _S2294 = 0.0f;
        _S2295 = 0.0f;
        _S2296 = 0.0f;
        _S2297 = 0.0f;
        _S2298 = 0.0f;
        _S2299 = 0.0f;
        _S2300 = 0.0f;
        _S2301 = 0.0f;
    }
    else
    {
        float _S2311 = 1.0f - y0_6;
        float _S2312 = 1.0f - x0_6;
        float s1_11 = _S2311 / _S2312;
        float _S2313 = _S2312 * _S2312;
        float _S2314 = _S2264 - x0_6;
        float t1_5 = _S2314 / _S2312;
        float _S2315 = s1_11 * t1_5;
        float _S2316 = _S2273 * t1_5;
        float _S2317 = 1.0f - t1_5;
        float _S2318 = _S2315 * t1_5 + _S2316 * _S2317;
        float _S2319 = _S2311 * _S2318;
        float _S2320 = _S2273 + _S2266 - 2.0f * s1_11;
        float _S2321 = _S2320 * t1_5;
        float _S2322 = s1_11 + _S2321 * _S2317;
        float _S2323 = _S2322 * _S2322;
        _S2275 = 0.0f;
        _S2276 = 0.0f;
        _S2277 = 0.0f;
        _S2278 = 0.0f;
        _S2279 = 0.0f;
        _S2280 = 0.0f;
        _S2281 = 0.0f;
        _S2282 = 0.0f;
        _S2283 = 0.0f;
        _S2284 = 0.0f;
        _S2285 = 0.0f;
        _S2286 = 0.0f;
        _S2287 = _S2323;
        _S2288 = _S2319;
        _S2289 = _S2322;
        _S2290 = _S2321;
        _S2291 = _S2317;
        _S2292 = _S2320;
        _S2293 = t1_5;
        _S2294 = _S2311;
        _S2295 = _S2318;
        _S2296 = _S2316;
        _S2297 = _S2315;
        _S2298 = s1_11;
        _S2299 = _S2313;
        _S2300 = _S2314;
        _S2301 = _S2312;
    }
    if(_S2274)
    {
        float _S2324 = _s_dOut_9.z / _S2275;
        float _S2325 = _S2276 * - _S2324;
        float _S2326 = _S2277 * _S2324;
        float _S2327 = _S2279 * _S2325;
        float _S2328 = _S2281 * _S2327;
        float _S2329 = y0_6 * _S2326;
        float _S2330 = _S2279 * _S2329;
        float _S2331 = _S2281 * _S2329;
        float _S2332 = (_S2280 * _S2327 + - (_S2278 * _S2325 + _S2283 * _S2329) + _S2265 * _S2330 + _S2284 * _S2329 + _S2285 * _S2331) / _S2286;
        float _S2333 = x0_6 * _S2332;
        float _S2334 = (_S2325 + 2.0f * - _S2328 + _S2281 * _S2331) / _S2286;
        float _S2335 = _S2282 * _S2326 + x0_6 * _S2334;
        float _S2336 = _S2328 + _S2281 * _S2330;
        float _S2337 = _S2264 * - _S2332 + y0_6 * - _S2334;
        _S2275 = _S2328;
        _S2276 = _S2335;
        _S2277 = _S2337;
        _S2278 = 0.0f;
        _S2279 = _S2336;
        _S2280 = _S2333;
    }
    else
    {
        float _S2338 = _s_dOut_9.z / _S2287;
        float _S2339 = _S2288 * - _S2338;
        float _S2340 = _S2289 * _S2338;
        float _S2341 = _S2291 * _S2339;
        float _S2342 = _S2293 * _S2341;
        float _S2343 = _S2294 * _S2340;
        float _S2344 = _S2291 * _S2343;
        float _S2345 = _S2293 * _S2343;
        float _S2346 = (_S2292 * _S2341 + - (_S2290 * _S2339 + _S2296 * _S2343) + _S2273 * _S2344 + _S2297 * _S2343 + _S2298 * _S2345) / _S2299;
        float _S2347 = _S2301 * _S2346;
        float _S2348 = (_S2339 + 2.0f * - _S2342 + _S2293 * _S2345) / _S2299;
        float _S2349 = _s_dOut_9.z + - (_S2295 * _S2340 + _S2301 * _S2348);
        float _S2350 = - _S2347 + - (_S2300 * - _S2346 + _S2294 * - _S2348);
        _S2275 = _S2342 + _S2293 * _S2344;
        _S2276 = _S2349;
        _S2277 = _S2350;
        _S2278 = _S2342;
        _S2279 = 0.0f;
        _S2280 = _S2347;
    }
    DiffPair_float_0 _S2351;
    (&_S2351)->primal_0 = _S2099;
    (&_S2351)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2351, _S2275);
    DiffPair_float_0 _S2352 = _S2351;
    float _S2353 = - (_S2276 / _S2272);
    DiffPair_float_0 _S2354;
    (&_S2354)->primal_0 = _S2270;
    (&_S2354)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2354, _S2353);
    float _S2355 = - _S2354.differential_0;
    float _S2356 = - (_S2277 / _S2269);
    DiffPair_float_0 _S2357;
    (&_S2357)->primal_0 = _S2267;
    (&_S2357)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2357, _S2356);
    float _S2358 = - _S2357.differential_0;
    DiffPair_float_0 _S2359;
    (&_S2359)->primal_0 = _S2096;
    (&_S2359)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2359, _S2278);
    DiffPair_float_0 _S2360 = _S2359;
    DiffPair_float_0 _S2361;
    (&_S2361)->primal_0 = _S2095;
    (&_S2361)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2361, _S2279);
    DiffPair_float_0 _S2362 = _S2361;
    float3  _S2363 = make_float3 (0.0f, 0.0f, _S2280);
    if(_S2214)
    {
        float _S2364 = _s_dOut_9.y / _S2215;
        float _S2365 = _S2216 * - _S2364;
        float _S2366 = _S2217 * _S2364;
        float _S2367 = _S2219 * _S2365;
        float _S2368 = _S2221 * _S2367;
        float _S2369 = y0_5 * _S2366;
        float _S2370 = _S2219 * _S2369;
        float _S2371 = _S2221 * _S2369;
        float _S2372 = (_S2220 * _S2367 + - (_S2218 * _S2365 + _S2223 * _S2369) + _S2205 * _S2370 + _S2224 * _S2369 + _S2225 * _S2371) / _S2226;
        float _S2373 = x0_5 * _S2372;
        float _S2374 = (_S2365 + 2.0f * - _S2368 + _S2221 * _S2371) / _S2226;
        float _S2375 = _S2222 * _S2366 + x0_5 * _S2374;
        float _S2376 = _S2368 + _S2221 * _S2370;
        float _S2377 = _S2204 * - _S2372 + y0_5 * - _S2374;
        _S2215 = _S2368;
        _S2216 = _S2375;
        _S2217 = _S2377;
        _S2218 = 0.0f;
        _S2219 = _S2376;
        _S2220 = _S2373;
    }
    else
    {
        float _S2378 = _s_dOut_9.y / _S2227;
        float _S2379 = _S2228 * - _S2378;
        float _S2380 = _S2229 * _S2378;
        float _S2381 = _S2231 * _S2379;
        float _S2382 = _S2233 * _S2381;
        float _S2383 = _S2234 * _S2380;
        float _S2384 = _S2231 * _S2383;
        float _S2385 = _S2233 * _S2383;
        float _S2386 = (_S2232 * _S2381 + - (_S2230 * _S2379 + _S2236 * _S2383) + _S2213 * _S2384 + _S2237 * _S2383 + _S2238 * _S2385) / _S2239;
        float _S2387 = _S2241 * _S2386;
        float _S2388 = (_S2379 + 2.0f * - _S2382 + _S2233 * _S2385) / _S2239;
        float _S2389 = _s_dOut_9.y + - (_S2235 * _S2380 + _S2241 * _S2388);
        float _S2390 = - _S2387 + - (_S2240 * - _S2386 + _S2234 * - _S2388);
        _S2215 = _S2382 + _S2233 * _S2384;
        _S2216 = _S2389;
        _S2217 = _S2390;
        _S2218 = _S2382;
        _S2219 = 0.0f;
        _S2220 = _S2387;
    }
    DiffPair_float_0 _S2391;
    (&_S2391)->primal_0 = _S2094;
    (&_S2391)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2391, _S2215);
    DiffPair_float_0 _S2392 = _S2391;
    float _S2393 = - (_S2216 / _S2212);
    DiffPair_float_0 _S2394;
    (&_S2394)->primal_0 = _S2210;
    (&_S2394)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2394, _S2393);
    float _S2395 = - _S2394.differential_0;
    float _S2396 = - (_S2217 / _S2209);
    DiffPair_float_0 _S2397;
    (&_S2397)->primal_0 = _S2207;
    (&_S2397)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2397, _S2396);
    float _S2398 = - _S2397.differential_0;
    DiffPair_float_0 _S2399;
    (&_S2399)->primal_0 = _S2091;
    (&_S2399)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2399, _S2218);
    DiffPair_float_0 _S2400 = _S2399;
    DiffPair_float_0 _S2401;
    (&_S2401)->primal_0 = _S2090;
    (&_S2401)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2401, _S2219);
    DiffPair_float_0 _S2402 = _S2401;
    float3  _S2403 = _S2363 + make_float3 (0.0f, _S2220, 0.0f);
    if(_S2154)
    {
        float _S2404 = _s_dOut_9.x / _S2155;
        float _S2405 = _S2156 * - _S2404;
        float _S2406 = _S2157 * _S2404;
        float _S2407 = _S2159 * _S2405;
        float _S2408 = _S2161 * _S2407;
        float _S2409 = y0_4 * _S2406;
        float _S2410 = _S2159 * _S2409;
        float _S2411 = _S2161 * _S2409;
        float _S2412 = (_S2160 * _S2407 + - (_S2158 * _S2405 + _S2163 * _S2409) + _S2145 * _S2410 + _S2164 * _S2409 + _S2165 * _S2411) / _S2166;
        float _S2413 = x0_4 * _S2412;
        float _S2414 = (_S2405 + 2.0f * - _S2408 + _S2161 * _S2411) / _S2166;
        float _S2415 = _S2162 * _S2406 + x0_4 * _S2414;
        float _S2416 = _S2408 + _S2161 * _S2410;
        float _S2417 = _S2144 * - _S2412 + y0_4 * - _S2414;
        _S2155 = _S2408;
        _S2156 = _S2415;
        _S2157 = _S2417;
        _S2158 = 0.0f;
        _S2159 = _S2416;
        _S2160 = _S2413;
    }
    else
    {
        float _S2418 = _s_dOut_9.x / _S2167;
        float _S2419 = _S2168 * - _S2418;
        float _S2420 = _S2169 * _S2418;
        float _S2421 = _S2171 * _S2419;
        float _S2422 = _S2173 * _S2421;
        float _S2423 = _S2174 * _S2420;
        float _S2424 = _S2171 * _S2423;
        float _S2425 = _S2173 * _S2423;
        float _S2426 = (_S2172 * _S2421 + - (_S2170 * _S2419 + _S2176 * _S2423) + _S2153 * _S2424 + _S2177 * _S2423 + _S2178 * _S2425) / _S2179;
        float _S2427 = _S2181 * _S2426;
        float _S2428 = (_S2419 + 2.0f * - _S2422 + _S2173 * _S2425) / _S2179;
        float _S2429 = _s_dOut_9.x + - (_S2175 * _S2420 + _S2181 * _S2428);
        float _S2430 = - _S2427 + - (_S2180 * - _S2426 + _S2174 * - _S2428);
        _S2155 = _S2422 + _S2173 * _S2424;
        _S2156 = _S2429;
        _S2157 = _S2430;
        _S2158 = _S2422;
        _S2159 = 0.0f;
        _S2160 = _S2427;
    }
    DiffPair_float_0 _S2431;
    (&_S2431)->primal_0 = _S2089;
    (&_S2431)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2431, _S2155);
    DiffPair_float_0 _S2432 = _S2431;
    float _S2433 = - (_S2156 / _S2152);
    DiffPair_float_0 _S2434;
    (&_S2434)->primal_0 = _S2150;
    (&_S2434)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2434, _S2433);
    float _S2435 = - _S2434.differential_0;
    float _S2436 = - (_S2157 / _S2149);
    DiffPair_float_0 _S2437;
    (&_S2437)->primal_0 = _S2147;
    (&_S2437)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2437, _S2436);
    float _S2438 = - _S2437.differential_0;
    DiffPair_float_0 _S2439;
    (&_S2439)->primal_0 = _S2086;
    (&_S2439)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2439, _S2158);
    DiffPair_float_0 _S2440 = _S2439;
    DiffPair_float_0 _S2441;
    (&_S2441)->primal_0 = _S2085;
    (&_S2441)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2441, _S2159);
    DiffPair_float_0 _S2442 = _S2441;
    float3  _S2443 = _S2403 + make_float3 (_S2160, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2444;
    (&_S2444)->primal_0 = _S2140;
    (&_S2444)->differential_0 = _S2067;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2445;
    (&_S2445)->primal_0 = _S2141;
    (&_S2445)->differential_0 = _S2067;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2446;
    (&_S2446)->primal_0 = _S2142;
    (&_S2446)->differential_0 = _S2067;
    s_bwd_prop_clamp_1(&_S2444, &_S2445, &_S2446, _S2443);
    float _S2447 = - _S2444.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S2444.differential_0.x + _S2447, _S2444.differential_0.y + _S2447, _S2444.differential_0.z);
    float3  _S2448 = _S2134 * s_diff_rgi_out_T_1;
    float _S2449 = (_S2448.x + _S2448.y + _S2448.z) / _S2137;
    float _S2450 = _S2135 * _S2449;
    float3  _S2451 = _S2136 * s_diff_rgi_out_T_1 + make_float3 (0.0f, 0.0f, intensity_3 * - _S2449);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2452;
    (&_S2452)->primal_0 = H_5;
    (&_S2452)->differential_0 = _S2068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2453;
    (&_S2453)->primal_0 = rgi_in_1;
    (&_S2453)->differential_0 = _S2067;
    s_bwd_prop_mul_0(&_S2452, &_S2453, _S2451);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2454 = _S2452;
    float _S2455 = _S2450 + _S2453.differential_0.z;
    float _S2456 = _S2453.differential_0.y + _S2455;
    float _S2457 = _S2453.differential_0.x + _S2455;
    float3  _S2458 = make_float3 (_S2457, _S2456, _S2455);
    if(_S2127)
    {
        Matrix<float, 3, 3>  _S2459 = _S2126 * _S2454.differential_0;
        Matrix<float, 3, 3>  _S2460 = _S2128 * _S2454.differential_0;
        _S2129 = - ((_S2459.rows[int(0)].x + _S2459.rows[int(0)].y + _S2459.rows[int(0)].z + _S2459.rows[int(1)].x + _S2459.rows[int(1)].y + _S2459.rows[int(1)].z + _S2459.rows[int(2)].x + _S2459.rows[int(2)].y + _S2459.rows[int(2)].z) / _S2129);
        H_5 = _S2460;
    }
    else
    {
        _S2129 = 0.0f;
        H_5 = _S2454.differential_0;
    }
    DiffPair_float_0 _S2461;
    (&_S2461)->primal_0 = _S2126.rows[int(2)].z;
    (&_S2461)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S2461, 0.0f);
    float _S2462 = _S2461.differential_0 + _S2129;
    float3  _S2463 = _S2067;
    *&((&_S2463)->z) = _S2462;
    Matrix<float, 3, 3>  _S2464 = _S2068;
    _S2464[int(2)] = _S2463;
    Matrix<float, 3, 3>  _S2465 = H_5 + _S2464;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2466;
    (&_S2466)->primal_0 = _S2125;
    (&_S2466)->differential_0 = _S2068;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2467;
    (&_S2467)->primal_0 = S_inv_1;
    (&_S2467)->differential_0 = _S2068;
    s_bwd_prop_mul_1(&_S2466, &_S2467, _S2465);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2468;
    (&_S2468)->primal_0 = T_3;
    (&_S2468)->differential_0 = _S2068;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2469;
    (&_S2469)->primal_0 = D_1;
    (&_S2469)->differential_0 = _S2068;
    s_bwd_prop_mul_1(&_S2468, &_S2469, _S2466.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2470 = _S2468;
    float3  _S2471 = make_float3 (_S2469.differential_0.rows[int(0)].x, _S2469.differential_0.rows[int(1)].y, _S2469.differential_0.rows[int(2)].z);
    float3  _S2472;
    if(_S2120)
    {
        if(_S2122)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2473;
            (&_S2473)->primal_0 = r1_3;
            (&_S2473)->differential_0 = _S2067;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S2474;
            (&_S2474)->primal_0 = r2_18;
            (&_S2474)->differential_0 = _S2067;
            s_bwd_prop_cross_0(&_S2473, &_S2474, _S2471);
            _S2108 = _S2067;
            lambda_v_7 = _S2474.differential_0;
            _S2472 = _S2473.differential_0;
        }
        else
        {
            _S2108 = _S2471;
            lambda_v_7 = _S2067;
            _S2472 = _S2067;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2475;
        (&_S2475)->primal_0 = _S2121;
        (&_S2475)->differential_0 = _S2067;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2476;
        (&_S2476)->primal_0 = _S2121;
        (&_S2476)->differential_0 = _S2067;
        s_bwd_prop_dot_0(&_S2475, &_S2476, 0.0f);
        float3  _S2477 = _S2476.differential_0 + _S2475.differential_0 + _S2108;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2478;
        (&_S2478)->primal_0 = r0_3;
        (&_S2478)->differential_0 = _S2067;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S2479;
        (&_S2479)->primal_0 = r2_18;
        (&_S2479)->differential_0 = _S2067;
        s_bwd_prop_cross_0(&_S2478, &_S2479, _S2477);
        float3  _S2480 = _S2479.differential_0 + lambda_v_7;
        _S2108 = _S2067;
        lambda_v_7 = _S2480;
        _S2121 = _S2472;
        _S2472 = _S2478.differential_0;
    }
    else
    {
        _S2108 = _S2471;
        lambda_v_7 = _S2067;
        _S2121 = _S2067;
        _S2472 = _S2067;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2481;
    (&_S2481)->primal_0 = _S2119;
    (&_S2481)->differential_0 = _S2067;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2482;
    (&_S2482)->primal_0 = _S2119;
    (&_S2482)->differential_0 = _S2067;
    s_bwd_prop_dot_0(&_S2481, &_S2482, 0.0f);
    float3  _S2483 = _S2482.differential_0 + _S2481.differential_0 + _S2108;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2484;
    (&_S2484)->primal_0 = r0_3;
    (&_S2484)->differential_0 = _S2067;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2485;
    (&_S2485)->primal_0 = r1_3;
    (&_S2485)->differential_0 = _S2067;
    s_bwd_prop_cross_0(&_S2484, &_S2485, _S2483);
    float3  _S2486 = _S2067;
    *&((&_S2486)->z) = lambda_v_7.z;
    *&((&_S2486)->y) = lambda_v_7.y;
    *&((&_S2486)->x) = lambda_v_7.x;
    float3  _S2487 = _S2485.differential_0 + _S2121;
    float3  _S2488 = _S2067;
    *&((&_S2488)->z) = _S2487.z;
    *&((&_S2488)->y) = _S2487.y;
    *&((&_S2488)->x) = _S2487.x;
    float3  _S2489 = _S2484.differential_0 + _S2472;
    float3  _S2490 = _S2067;
    *&((&_S2490)->z) = _S2489.z;
    *&((&_S2490)->y) = _S2489.y;
    *&((&_S2490)->x) = _S2489.x;
    Matrix<float, 3, 3>  _S2491 = _S2068;
    _S2491[int(2)] = _S2486;
    _S2491[int(1)] = _S2488;
    _S2491[int(0)] = _S2490;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2492;
    (&_S2492)->primal_0 = skew_1;
    (&_S2492)->differential_0 = _S2068;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2493;
    (&_S2493)->primal_0 = T_3;
    (&_S2493)->differential_0 = _S2068;
    s_bwd_prop_mul_1(&_S2492, &_S2493, _S2491);
    Matrix<float, 3, 3>  _S2494 = _S2493.differential_0 + _S2470.differential_0;
    float2  _S2495 = make_float2 (_S2492.differential_0.rows[int(2)].y + - _S2492.differential_0.rows[int(1)].z, _S2492.differential_0.rows[int(0)].z + - _S2492.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S2496 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2497;
    (&_S2497)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2497)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2498;
    (&_S2498)->primal_0 = _S2111.color_params_0.n_0;
    (&_S2498)->differential_0 = _S2071;
    s_bwd_prop_mul_2(&_S2497, &_S2498, _S2495);
    float2  _S2499 = make_float2 (_S2494.rows[int(0)].z, _S2494.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2500;
    (&_S2500)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2500)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2501;
    (&_S2501)->primal_0 = _S2111.color_params_0.g_0;
    (&_S2501)->differential_0 = _S2071;
    s_bwd_prop_mul_2(&_S2500, &_S2501, _S2499);
    float2  _S2502 = make_float2 (_S2494.rows[int(0)].y, _S2494.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2503;
    (&_S2503)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2503)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2504;
    (&_S2504)->primal_0 = _S2111.color_params_0.r_0;
    (&_S2504)->differential_0 = _S2071;
    s_bwd_prop_mul_2(&_S2503, &_S2504, _S2502);
    float2  _S2505 = make_float2 (_S2494.rows[int(0)].x, _S2494.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2506;
    (&_S2506)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2506)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2507;
    (&_S2507)->primal_0 = _S2111.color_params_0.b_0;
    (&_S2507)->differential_0 = _S2071;
    s_bwd_prop_mul_2(&_S2506, &_S2507, _S2505);
    ColorPPISPParams_0 _S2508 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2508)->n_0 = _S2498.differential_0;
    (&_S2508)->g_0 = _S2501.differential_0;
    (&_S2508)->r_0 = _S2504.differential_0;
    (&_S2508)->b_0 = _S2507.differential_0;
    _S2108 = _S2458;
    *&((&_S2108)->z) = 0.0f;
    float _S2509 = rgb_out_5.z * _S2455;
    float _S2510 = _S2110 * _S2455;
    DiffPair_float_0 _S2511;
    (&_S2511)->primal_0 = falloff_5;
    (&_S2511)->differential_0 = 0.0f;
    DiffPair_float_0 _S2512;
    (&_S2512)->primal_0 = 0.0f;
    (&_S2512)->differential_0 = 0.0f;
    DiffPair_float_0 _S2513;
    (&_S2513)->primal_0 = 1.0f;
    (&_S2513)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2511, &_S2512, &_S2513, _S2509);
    float _S2514 = r2_17 * _S2511.differential_0;
    float _S2515 = r4_11 * _S2511.differential_0;
    float s_diff_r6_T_3 = _S2084 * _S2511.differential_0;
    float _S2516 = r6_5 * _S2511.differential_0;
    float _S2517 = r2_17 * (_S2083 * _S2511.differential_0 + r2_17 * s_diff_r6_T_3);
    float _S2518 = _S2082 * _S2511.differential_0 + r4_11 * s_diff_r6_T_3 + _S2517 + _S2517;
    float _S2519 = dy_11 * _S2518;
    float _S2520 = dx_13 * _S2518;
    float _S2521 = - (_S2519 + _S2519);
    float _S2522 = - (_S2520 + _S2520);
    *&((&_S2108)->y) = 0.0f;
    float _S2523 = rgb_out_5.y * _S2456;
    float _S2524 = _S2109 * _S2456;
    DiffPair_float_0 _S2525;
    (&_S2525)->primal_0 = falloff_4;
    (&_S2525)->differential_0 = 0.0f;
    DiffPair_float_0 _S2526;
    (&_S2526)->primal_0 = 0.0f;
    (&_S2526)->differential_0 = 0.0f;
    DiffPair_float_0 _S2527;
    (&_S2527)->primal_0 = 1.0f;
    (&_S2527)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2525, &_S2526, &_S2527, _S2523);
    float _S2528 = r2_16 * _S2525.differential_0;
    float _S2529 = r4_10 * _S2525.differential_0;
    float s_diff_r6_T_4 = _S2081 * _S2525.differential_0;
    float _S2530 = r6_4 * _S2525.differential_0;
    float _S2531 = r2_16 * (_S2080 * _S2525.differential_0 + r2_16 * s_diff_r6_T_4);
    float _S2532 = _S2079 * _S2525.differential_0 + r4_10 * s_diff_r6_T_4 + _S2531 + _S2531;
    float _S2533 = dy_10 * _S2532;
    float _S2534 = dx_12 * _S2532;
    float _S2535 = - (_S2533 + _S2533);
    float _S2536 = - (_S2534 + _S2534);
    *&((&_S2108)->x) = 0.0f;
    float _S2537 = rgb_out_5.x * _S2457;
    float _S2538 = _S2106 * _S2457;
    DiffPair_float_0 _S2539;
    (&_S2539)->primal_0 = falloff_3;
    (&_S2539)->differential_0 = 0.0f;
    DiffPair_float_0 _S2540;
    (&_S2540)->primal_0 = 0.0f;
    (&_S2540)->differential_0 = 0.0f;
    DiffPair_float_0 _S2541;
    (&_S2541)->primal_0 = 1.0f;
    (&_S2541)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2539, &_S2540, &_S2541, _S2537);
    float _S2542 = r2_15 * _S2539.differential_0;
    float _S2543 = r4_9 * _S2539.differential_0;
    float s_diff_r6_T_5 = _S2078 * _S2539.differential_0;
    float _S2544 = r6_3 * _S2539.differential_0;
    float _S2545 = r2_15 * (_S2077 * _S2539.differential_0 + r2_15 * s_diff_r6_T_5);
    float _S2546 = _S2076 * _S2539.differential_0 + r4_9 * s_diff_r6_T_5 + _S2545 + _S2545;
    float _S2547 = dy_9 * _S2546;
    float _S2548 = dx_11 * _S2546;
    float _S2549 = - (_S2547 + _S2547);
    float _S2550 = - (_S2548 + _S2548);
    float3  _S2551 = _S2067;
    *&((&_S2551)->z) = _S2510;
    *&((&_S2551)->y) = _S2524;
    *&((&_S2551)->x) = _S2538;
    float3  _S2552 = _S2108 + _S2551;
    float3  _S2553 = _S2066.primal_0 * _S2552;
    float3  _S2554 = _S2102 * _S2552;
    float _S2555 = _S2553.x + _S2553.y + _S2553.z;
    DiffPair_float_0 _S2556;
    (&_S2556)->primal_0 = _S2100.exposure_0;
    (&_S2556)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2556, _S2555);
    PPISPParamsRQS_0 _S2557 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2557)->color_params_0 = _S2508;
    (&_S2557)->exposure_0 = _S2556.differential_0;
    _S2075 = _S2557;
    (&(&_S2075)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2558 = _S2557.crf_params_0[int(2)].gc_0 + _S2352.differential_0;
    (&(&_S2075)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2559 = _S2557.crf_params_0[int(2)].y0_0 + _S2355;
    (&(&_S2075)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2560 = _S2557.crf_params_0[int(2)].x0_0 + _S2358;
    (&(&_S2075)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2561 = _S2557.crf_params_0[int(2)].g1_0 + _S2360.differential_0;
    (&(&_S2075)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2562 = _S2557.crf_params_0[int(2)].g0_0 + _S2362.differential_0;
    (&(&_S2075)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2563 = _S2557.crf_params_0[int(1)].gc_0 + _S2392.differential_0;
    (&(&_S2075)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2564 = _S2557.crf_params_0[int(1)].y0_0 + _S2395;
    (&(&_S2075)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2565 = _S2557.crf_params_0[int(1)].x0_0 + _S2398;
    (&(&_S2075)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2566 = _S2557.crf_params_0[int(1)].g1_0 + _S2400.differential_0;
    (&(&_S2075)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2567 = _S2557.crf_params_0[int(1)].g0_0 + _S2402.differential_0;
    (&(&_S2075)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2568 = _S2557.crf_params_0[int(0)].gc_0 + _S2432.differential_0;
    (&(&_S2075)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2569 = _S2557.crf_params_0[int(0)].y0_0 + _S2435;
    (&(&_S2075)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2570 = _S2557.crf_params_0[int(0)].x0_0 + _S2438;
    (&(&_S2075)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2571 = _S2557.crf_params_0[int(0)].g1_0 + _S2440.differential_0;
    (&(&_S2075)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2572 = _S2557.crf_params_0[int(0)].g0_0 + _S2442.differential_0;
    *&((&(&(&_S2075)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2075)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2075)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2573 = _S2516 + _S2557.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2075)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2574 = _S2515 + _S2557.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2075)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2575 = _S2514 + _S2557.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2075)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2576 = _S2521 + _S2557.vignette_params_0[int(2)].cy_0;
    (&(&_S2075)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2577 = _S2522 + _S2557.vignette_params_0[int(2)].cx_0;
    (&(&_S2075)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2578 = _S2530 + _S2557.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2075)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2579 = _S2529 + _S2557.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2075)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2580 = _S2528 + _S2557.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2075)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2581 = _S2535 + _S2557.vignette_params_0[int(1)].cy_0;
    (&(&_S2075)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2582 = _S2536 + _S2557.vignette_params_0[int(1)].cx_0;
    (&(&_S2075)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S2583 = _S2544 + _S2557.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2075)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S2584 = _S2543 + _S2557.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2075)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S2585 = _S2542 + _S2557.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2075)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S2586 = _S2549 + _S2557.vignette_params_0[int(0)].cy_0;
    (&(&_S2075)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S2587 = _S2550 + _S2557.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S2588;
    _S2588[int(0)] = 0.0f;
    _S2588[int(1)] = 0.0f;
    _S2588[int(2)] = 0.0f;
    _S2588[int(3)] = 0.0f;
    _S2588[int(4)] = 0.0f;
    _S2588[int(5)] = 0.0f;
    _S2588[int(6)] = 0.0f;
    _S2588[int(7)] = 0.0f;
    _S2588[int(8)] = 0.0f;
    _S2588[int(9)] = 0.0f;
    _S2588[int(10)] = 0.0f;
    _S2588[int(11)] = 0.0f;
    _S2588[int(12)] = 0.0f;
    _S2588[int(13)] = 0.0f;
    _S2588[int(14)] = 0.0f;
    _S2588[int(15)] = 0.0f;
    _S2588[int(16)] = 0.0f;
    _S2588[int(17)] = 0.0f;
    _S2588[int(18)] = 0.0f;
    _S2588[int(19)] = 0.0f;
    _S2588[int(20)] = 0.0f;
    _S2588[int(21)] = 0.0f;
    _S2588[int(22)] = 0.0f;
    _S2588[int(23)] = 0.0f;
    _S2588[int(24)] = 0.0f;
    _S2588[int(25)] = 0.0f;
    _S2588[int(26)] = 0.0f;
    _S2588[int(27)] = 0.0f;
    _S2588[int(28)] = 0.0f;
    _S2588[int(29)] = 0.0f;
    _S2588[int(30)] = 0.0f;
    _S2588[int(31)] = 0.0f;
    _S2588[int(32)] = 0.0f;
    _S2588[int(33)] = 0.0f;
    _S2588[int(34)] = 0.0f;
    _S2588[int(35)] = 0.0f;
    _S2588[int(36)] = 0.0f;
    _S2588[int(37)] = 0.0f;
    _S2588[int(38)] = 0.0f;
    _S2588[int(9)] = _S2579;
    _S2588[int(18)] = _S2557.color_params_0.r_0.x;
    _S2588[int(17)] = _S2557.color_params_0.b_0.y;
    _S2588[int(16)] = _S2557.color_params_0.b_0.x;
    _S2588[int(15)] = _S2573;
    _S2588[int(14)] = _S2574;
    _S2588[int(13)] = _S2575;
    _S2588[int(12)] = _S2576;
    _S2588[int(11)] = _S2577;
    _S2588[int(10)] = _S2578;
    _S2588[int(19)] = _S2557.color_params_0.r_0.y;
    _S2588[int(8)] = _S2580;
    _S2588[int(7)] = _S2581;
    _S2588[int(6)] = _S2582;
    _S2588[int(5)] = _S2583;
    _S2588[int(4)] = _S2584;
    _S2588[int(3)] = _S2585;
    _S2588[int(2)] = _S2586;
    _S2588[int(1)] = _S2587;
    _S2588[int(0)] = _S2075.exposure_0;
    _S2588[int(28)] = _S2568;
    _S2588[int(37)] = _S2559;
    _S2588[int(36)] = _S2560;
    _S2588[int(35)] = _S2561;
    _S2588[int(34)] = _S2562;
    _S2588[int(33)] = _S2563;
    _S2588[int(32)] = _S2564;
    _S2588[int(31)] = _S2565;
    _S2588[int(30)] = _S2566;
    _S2588[int(29)] = _S2567;
    _S2588[int(38)] = _S2558;
    _S2588[int(27)] = _S2569;
    _S2588[int(26)] = _S2570;
    _S2588[int(25)] = _S2571;
    _S2588[int(24)] = _S2572;
    _S2588[int(23)] = _S2557.color_params_0.n_0.y;
    _S2588[int(22)] = _S2557.color_params_0.n_0.x;
    _S2588[int(21)] = _S2557.color_params_0.g_0.y;
    _S2588[int(20)] = _S2557.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S2588;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S2554;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2589, float2  _S2590, float2  _S2591, float2  _S2592, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S2593, float3  _S2594)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S2589, _S2590, _S2591, _S2592, _S2593, _S2594);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S2595 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S2595;
    FixedArray<float, 39>  _S2596 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_3;
    (&dp_params_1)->differential_0 = _S2596;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_4, FixedArray<float, 22>  * _S2597)
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
    FixedArray<float, 22>  losses_4;
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
    losses_4[int(0)] = p_2.exposure_1;
    float _S2598 = p_2.vignette_params_1[int(0)].cx_0;
    float _S2599 = p_2.vignette_params_1[int(0)].cy_0;
    float _S2600 = p_2.vignette_params_1[int(1)].cx_0;
    float _S2601 = p_2.vignette_params_1[int(1)].cy_0;
    float _S2602 = p_2.vignette_params_1[int(2)].cx_0;
    float _S2603 = p_2.vignette_params_1[int(2)].cy_0;
    losses_4[int(1)] = _S2598 * _S2598 + _S2599 * _S2599 + _S2600 * _S2600 + _S2601 * _S2601 + _S2602 * _S2602 + _S2603 * _S2603;
    losses_4[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_4[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_4[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_2 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S2604 = p_2.vignette_params_1[int(0)].cx_0 - mean_2;
    float _S2605 = p_2.vignette_params_1[int(1)].cx_0 - mean_2;
    float _S2606 = p_2.vignette_params_1[int(2)].cx_0 - mean_2;
    losses_4[int(5)] = (_S2604 * _S2604 + _S2605 * _S2605 + _S2606 * _S2606) / 3.0f;
    float mean_3 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S2607 = p_2.vignette_params_1[int(0)].cy_0 - mean_3;
    float _S2608 = p_2.vignette_params_1[int(1)].cy_0 - mean_3;
    float _S2609 = p_2.vignette_params_1[int(2)].cy_0 - mean_3;
    losses_4[int(6)] = (_S2607 * _S2607 + _S2608 * _S2608 + _S2609 * _S2609) / 3.0f;
    float mean_4 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S2610 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_4;
    float _S2611 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_4;
    float _S2612 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_4;
    losses_4[int(7)] = (_S2610 * _S2610 + _S2611 * _S2611 + _S2612 * _S2612) / 3.0f;
    float mean_5 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S2613 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_5;
    float _S2614 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_5;
    float _S2615 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_5;
    losses_4[int(8)] = (_S2613 * _S2613 + _S2614 * _S2614 + _S2615 * _S2615) / 3.0f;
    float mean_6 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S2616 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_6;
    float _S2617 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_6;
    float _S2618 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_6;
    losses_4[int(9)] = (_S2616 * _S2616 + _S2617 * _S2617 + _S2618 * _S2618) / 3.0f;
    float2  bd_2 = mul_3(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_2.color_params_1.b_0);
    float2  rd_2 = mul_3(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_2.color_params_1.r_0);
    float2  gd_2 = mul_3(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_2.color_params_1.g_0);
    float2  nd_2 = mul_3(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_2.color_params_1.n_0);
    losses_4[int(10)] = bd_2.x;
    losses_4[int(11)] = bd_2.y;
    losses_4[int(12)] = rd_2.x;
    losses_4[int(13)] = rd_2.y;
    losses_4[int(14)] = gd_2.x;
    losses_4[int(15)] = gd_2.y;
    losses_4[int(16)] = nd_2.x;
    losses_4[int(17)] = nd_2.y;
    float mean_7 = (p_2.crf_params_1[int(0)].toe_0 + p_2.crf_params_1[int(1)].toe_0 + p_2.crf_params_1[int(2)].toe_0) / 3.0f;
    float _S2619 = p_2.crf_params_1[int(0)].toe_0 - mean_7;
    float _S2620 = p_2.crf_params_1[int(1)].toe_0 - mean_7;
    float _S2621 = p_2.crf_params_1[int(2)].toe_0 - mean_7;
    losses_4[int(18)] = (_S2619 * _S2619 + _S2620 * _S2620 + _S2621 * _S2621) / 3.0f;
    float mean_8 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S2622 = p_2.crf_params_1[int(0)].shoulder_0 - mean_8;
    float _S2623 = p_2.crf_params_1[int(1)].shoulder_0 - mean_8;
    float _S2624 = p_2.crf_params_1[int(2)].shoulder_0 - mean_8;
    losses_4[int(19)] = (_S2622 * _S2622 + _S2623 * _S2623 + _S2624 * _S2624) / 3.0f;
    float mean_9 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S2625 = p_2.crf_params_1[int(0)].gamma_0 - mean_9;
    float _S2626 = p_2.crf_params_1[int(1)].gamma_0 - mean_9;
    float _S2627 = p_2.crf_params_1[int(2)].gamma_0 - mean_9;
    losses_4[int(20)] = (_S2625 * _S2625 + _S2626 * _S2626 + _S2627 * _S2627) / 3.0f;
    float mean_10 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S2628 = p_2.crf_params_1[int(0)].center_0 - mean_10;
    float _S2629 = p_2.crf_params_1[int(1)].center_0 - mean_10;
    float _S2630 = p_2.crf_params_1[int(2)].center_0 - mean_10;
    losses_4[int(21)] = (_S2628 * _S2628 + _S2629 * _S2629 + _S2630 * _S2630) / 3.0f;
    *_S2597 = losses_4;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_5, FixedArray<float, 23>  * _S2631)
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
    FixedArray<float, 23>  losses_5;
    losses_5[int(0)] = 0.0f;
    losses_5[int(1)] = 0.0f;
    losses_5[int(2)] = 0.0f;
    losses_5[int(3)] = 0.0f;
    losses_5[int(4)] = 0.0f;
    losses_5[int(5)] = 0.0f;
    losses_5[int(6)] = 0.0f;
    losses_5[int(7)] = 0.0f;
    losses_5[int(8)] = 0.0f;
    losses_5[int(9)] = 0.0f;
    losses_5[int(10)] = 0.0f;
    losses_5[int(11)] = 0.0f;
    losses_5[int(12)] = 0.0f;
    losses_5[int(13)] = 0.0f;
    losses_5[int(14)] = 0.0f;
    losses_5[int(15)] = 0.0f;
    losses_5[int(16)] = 0.0f;
    losses_5[int(17)] = 0.0f;
    losses_5[int(18)] = 0.0f;
    losses_5[int(19)] = 0.0f;
    losses_5[int(20)] = 0.0f;
    losses_5[int(21)] = 0.0f;
    losses_5[int(22)] = 0.0f;
    losses_5[int(0)] = p_3.exposure_0;
    float _S2632 = p_3.vignette_params_0[int(0)].cx_0;
    float _S2633 = p_3.vignette_params_0[int(0)].cy_0;
    float _S2634 = p_3.vignette_params_0[int(1)].cx_0;
    float _S2635 = p_3.vignette_params_0[int(1)].cy_0;
    float _S2636 = p_3.vignette_params_0[int(2)].cx_0;
    float _S2637 = p_3.vignette_params_0[int(2)].cy_0;
    losses_5[int(1)] = _S2632 * _S2632 + _S2633 * _S2633 + _S2634 * _S2634 + _S2635 * _S2635 + _S2636 * _S2636 + _S2637 * _S2637;
    losses_5[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_5[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_5[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_11 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S2638 = p_3.vignette_params_0[int(0)].cx_0 - mean_11;
    float _S2639 = p_3.vignette_params_0[int(1)].cx_0 - mean_11;
    float _S2640 = p_3.vignette_params_0[int(2)].cx_0 - mean_11;
    losses_5[int(5)] = (_S2638 * _S2638 + _S2639 * _S2639 + _S2640 * _S2640) / 3.0f;
    float mean_12 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S2641 = p_3.vignette_params_0[int(0)].cy_0 - mean_12;
    float _S2642 = p_3.vignette_params_0[int(1)].cy_0 - mean_12;
    float _S2643 = p_3.vignette_params_0[int(2)].cy_0 - mean_12;
    losses_5[int(6)] = (_S2641 * _S2641 + _S2642 * _S2642 + _S2643 * _S2643) / 3.0f;
    float mean_13 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S2644 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_13;
    float _S2645 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_13;
    float _S2646 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_13;
    losses_5[int(7)] = (_S2644 * _S2644 + _S2645 * _S2645 + _S2646 * _S2646) / 3.0f;
    float mean_14 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S2647 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_14;
    float _S2648 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_14;
    float _S2649 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_14;
    losses_5[int(8)] = (_S2647 * _S2647 + _S2648 * _S2648 + _S2649 * _S2649) / 3.0f;
    float mean_15 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S2650 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_15;
    float _S2651 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_15;
    float _S2652 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_15;
    losses_5[int(9)] = (_S2650 * _S2650 + _S2651 * _S2651 + _S2652 * _S2652) / 3.0f;
    float2  bd_3 = mul_3(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_3.color_params_0.b_0);
    float2  rd_3 = mul_3(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_3.color_params_0.r_0);
    float2  gd_3 = mul_3(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_3.color_params_0.g_0);
    float2  nd_3 = mul_3(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_3.color_params_0.n_0);
    losses_5[int(10)] = bd_3.x;
    losses_5[int(11)] = bd_3.y;
    losses_5[int(12)] = rd_3.x;
    losses_5[int(13)] = rd_3.y;
    losses_5[int(14)] = gd_3.x;
    losses_5[int(15)] = gd_3.y;
    losses_5[int(16)] = nd_3.x;
    losses_5[int(17)] = nd_3.y;
    float mean_16 = (p_3.crf_params_0[int(0)].g0_0 + p_3.crf_params_0[int(1)].g0_0 + p_3.crf_params_0[int(2)].g0_0) / 3.0f;
    float _S2653 = p_3.crf_params_0[int(0)].g0_0 - mean_16;
    float _S2654 = p_3.crf_params_0[int(1)].g0_0 - mean_16;
    float _S2655 = p_3.crf_params_0[int(2)].g0_0 - mean_16;
    losses_5[int(18)] = (_S2653 * _S2653 + _S2654 * _S2654 + _S2655 * _S2655) / 3.0f;
    float mean_17 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S2656 = p_3.crf_params_0[int(0)].g1_0 - mean_17;
    float _S2657 = p_3.crf_params_0[int(1)].g1_0 - mean_17;
    float _S2658 = p_3.crf_params_0[int(2)].g1_0 - mean_17;
    losses_5[int(19)] = (_S2656 * _S2656 + _S2657 * _S2657 + _S2658 * _S2658) / 3.0f;
    float mean_18 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S2659 = p_3.crf_params_0[int(0)].x0_0 - mean_18;
    float _S2660 = p_3.crf_params_0[int(1)].x0_0 - mean_18;
    float _S2661 = p_3.crf_params_0[int(2)].x0_0 - mean_18;
    losses_5[int(20)] = (_S2659 * _S2659 + _S2660 * _S2660 + _S2661 * _S2661) / 3.0f;
    float mean_19 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S2662 = p_3.crf_params_0[int(0)].y0_0 - mean_19;
    float _S2663 = p_3.crf_params_0[int(1)].y0_0 - mean_19;
    float _S2664 = p_3.crf_params_0[int(2)].y0_0 - mean_19;
    losses_5[int(21)] = (_S2662 * _S2662 + _S2663 * _S2663 + _S2664 * _S2664) / 3.0f;
    float mean_20 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S2665 = p_3.crf_params_0[int(0)].gc_0 - mean_20;
    float _S2666 = p_3.crf_params_0[int(1)].gc_0 - mean_20;
    float _S2667 = p_3.crf_params_0[int(2)].gc_0 - mean_20;
    losses_5[int(22)] = (_S2665 * _S2665 + _S2666 * _S2666 + _S2667 * _S2667) / 3.0f;
    *_S2631 = losses_5;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_10)
{
    VignettingChannelParams_0 _S2668 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2669 = {
        _S2668, _S2668, _S2668
    };
    float2  _S2670 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2671 = { _S2670, _S2670, _S2670, _S2670 };
    CRFPPISPChannelParams_0 _S2672 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S2673 = {
        _S2672, _S2672, _S2672
    };
    PPISPParams_0 _S2674;
    (&_S2674)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S2674)->vignette_params_1 = _S2669;
    (&_S2674)->color_params_1 = _S2671;
    (&_S2674)->crf_params_1 = _S2673;
    (&(&_S2674)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S2674)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S2674)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S2674)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S2674)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S2674)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S2674)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S2674)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S2674)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S2674)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S2674)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S2674)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S2674)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S2674)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S2674)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S2674)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S2674)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S2674)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S2674)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S2674)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S2674)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S2674)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S2674)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S2674)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S2674)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S2674)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S2674)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S2674)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S2674)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S2674)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S2674)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S2674)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S2674)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S2674)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S2674)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_21 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S2675 = dpparams_2->primal_0[int(1)] - mean_21;
    float _S2676 = dpparams_2->primal_0[int(6)] - mean_21;
    float _S2677 = dpparams_2->primal_0[int(11)] - mean_21;
    float mean_22 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S2678 = dpparams_2->primal_0[int(2)] - mean_22;
    float _S2679 = dpparams_2->primal_0[int(7)] - mean_22;
    float _S2680 = dpparams_2->primal_0[int(12)] - mean_22;
    float mean_23 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S2681 = dpparams_2->primal_0[int(3)] - mean_23;
    float _S2682 = dpparams_2->primal_0[int(8)] - mean_23;
    float _S2683 = dpparams_2->primal_0[int(13)] - mean_23;
    float mean_24 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S2684 = dpparams_2->primal_0[int(4)] - mean_24;
    float _S2685 = dpparams_2->primal_0[int(9)] - mean_24;
    float _S2686 = dpparams_2->primal_0[int(14)] - mean_24;
    float mean_25 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S2687 = dpparams_2->primal_0[int(5)] - mean_25;
    float _S2688 = dpparams_2->primal_0[int(10)] - mean_25;
    float _S2689 = dpparams_2->primal_0[int(15)] - mean_25;
    float mean_26 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_27 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_28 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_29 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S2690 = 0.3333333432674408f * (*_s_dOut_10)[int(21)];
    float _S2691 = (dpparams_2->primal_0[int(35)] - mean_29) * _S2690;
    float _S2692 = _S2691 + _S2691;
    float _S2693 = (dpparams_2->primal_0[int(31)] - mean_29) * _S2690;
    float _S2694 = _S2693 + _S2693;
    float _S2695 = (dpparams_2->primal_0[int(27)] - mean_29) * _S2690;
    float _S2696 = _S2695 + _S2695;
    float _S2697 = 0.3333333432674408f * (- _S2692 + - _S2694 + - _S2696);
    float _S2698 = 0.3333333432674408f * (*_s_dOut_10)[int(20)];
    float _S2699 = (dpparams_2->primal_0[int(34)] - mean_28) * _S2698;
    float _S2700 = _S2699 + _S2699;
    float _S2701 = (dpparams_2->primal_0[int(30)] - mean_28) * _S2698;
    float _S2702 = _S2701 + _S2701;
    float _S2703 = (dpparams_2->primal_0[int(26)] - mean_28) * _S2698;
    float _S2704 = _S2703 + _S2703;
    float _S2705 = 0.3333333432674408f * (- _S2700 + - _S2702 + - _S2704);
    float _S2706 = 0.3333333432674408f * (*_s_dOut_10)[int(19)];
    float _S2707 = (dpparams_2->primal_0[int(33)] - mean_27) * _S2706;
    float _S2708 = _S2707 + _S2707;
    float _S2709 = (dpparams_2->primal_0[int(29)] - mean_27) * _S2706;
    float _S2710 = _S2709 + _S2709;
    float _S2711 = (dpparams_2->primal_0[int(25)] - mean_27) * _S2706;
    float _S2712 = _S2711 + _S2711;
    float _S2713 = 0.3333333432674408f * (- _S2708 + - _S2710 + - _S2712);
    float _S2714 = 0.3333333432674408f * (*_s_dOut_10)[int(18)];
    float _S2715 = (dpparams_2->primal_0[int(32)] - mean_26) * _S2714;
    float _S2716 = _S2715 + _S2715;
    float _S2717 = (dpparams_2->primal_0[int(28)] - mean_26) * _S2714;
    float _S2718 = _S2717 + _S2717;
    float _S2719 = (dpparams_2->primal_0[int(24)] - mean_26) * _S2714;
    float _S2720 = _S2719 + _S2719;
    float _S2721 = 0.3333333432674408f * (- _S2716 + - _S2718 + - _S2720);
    float2  _S2722 = make_float2 ((*_s_dOut_10)[int(16)], (*_s_dOut_10)[int(17)]);
    Matrix<float, 2, 2>  _S2723 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2724;
    (&_S2724)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2724)->differential_0 = _S2723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2725;
    (&_S2725)->primal_0 = _S2674.color_params_1.n_0;
    (&_S2725)->differential_0 = _S2670;
    s_bwd_prop_mul_2(&_S2724, &_S2725, _S2722);
    float2  _S2726 = make_float2 ((*_s_dOut_10)[int(14)], (*_s_dOut_10)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2727;
    (&_S2727)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2727)->differential_0 = _S2723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2728;
    (&_S2728)->primal_0 = _S2674.color_params_1.g_0;
    (&_S2728)->differential_0 = _S2670;
    s_bwd_prop_mul_2(&_S2727, &_S2728, _S2726);
    float2  _S2729 = make_float2 ((*_s_dOut_10)[int(12)], (*_s_dOut_10)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2730;
    (&_S2730)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2730)->differential_0 = _S2723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2731;
    (&_S2731)->primal_0 = _S2674.color_params_1.r_0;
    (&_S2731)->differential_0 = _S2670;
    s_bwd_prop_mul_2(&_S2730, &_S2731, _S2729);
    float2  _S2732 = make_float2 ((*_s_dOut_10)[int(10)], (*_s_dOut_10)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2733;
    (&_S2733)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2733)->differential_0 = _S2723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2734;
    (&_S2734)->primal_0 = _S2674.color_params_1.b_0;
    (&_S2734)->differential_0 = _S2670;
    s_bwd_prop_mul_2(&_S2733, &_S2734, _S2732);
    ColorPPISPParams_0 _S2735 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2735)->n_0 = _S2725.differential_0;
    (&_S2735)->g_0 = _S2728.differential_0;
    (&_S2735)->r_0 = _S2731.differential_0;
    (&_S2735)->b_0 = _S2734.differential_0;
    float _S2736 = 0.3333333432674408f * (*_s_dOut_10)[int(9)];
    float _S2737 = _S2689 * _S2736;
    float _S2738 = _S2737 + _S2737;
    float _S2739 = _S2688 * _S2736;
    float _S2740 = _S2739 + _S2739;
    float _S2741 = _S2687 * _S2736;
    float _S2742 = _S2741 + _S2741;
    float _S2743 = 0.3333333432674408f * (- _S2738 + - _S2740 + - _S2742);
    float _S2744 = 0.3333333432674408f * (*_s_dOut_10)[int(8)];
    float _S2745 = _S2686 * _S2744;
    float _S2746 = _S2745 + _S2745;
    float _S2747 = _S2685 * _S2744;
    float _S2748 = _S2747 + _S2747;
    float _S2749 = _S2684 * _S2744;
    float _S2750 = _S2749 + _S2749;
    float _S2751 = 0.3333333432674408f * (- _S2746 + - _S2748 + - _S2750);
    float _S2752 = 0.3333333432674408f * (*_s_dOut_10)[int(7)];
    float _S2753 = _S2683 * _S2752;
    float _S2754 = _S2753 + _S2753;
    float _S2755 = _S2682 * _S2752;
    float _S2756 = _S2755 + _S2755;
    float _S2757 = _S2681 * _S2752;
    float _S2758 = _S2757 + _S2757;
    float _S2759 = 0.3333333432674408f * (- _S2754 + - _S2756 + - _S2758);
    float _S2760 = 0.3333333432674408f * (*_s_dOut_10)[int(6)];
    float _S2761 = _S2680 * _S2760;
    float _S2762 = _S2761 + _S2761;
    float _S2763 = _S2679 * _S2760;
    float _S2764 = _S2763 + _S2763;
    float _S2765 = _S2678 * _S2760;
    float _S2766 = _S2765 + _S2765;
    float _S2767 = 0.3333333432674408f * (- _S2762 + - _S2764 + - _S2766);
    float _S2768 = 0.3333333432674408f * (*_s_dOut_10)[int(5)];
    float _S2769 = _S2677 * _S2768;
    float _S2770 = _S2769 + _S2769;
    float _S2771 = _S2676 * _S2768;
    float _S2772 = _S2771 + _S2771;
    float _S2773 = _S2675 * _S2768;
    float _S2774 = _S2773 + _S2773;
    float _S2775 = 0.3333333432674408f * (- _S2770 + - _S2772 + - _S2774);
    DiffPair_float_0 _S2776;
    (&_S2776)->primal_0 = 0.0f;
    (&_S2776)->differential_0 = 0.0f;
    DiffPair_float_0 _S2777;
    (&_S2777)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S2777)->differential_0 = 0.0f;
    _d_max_0(&_S2776, &_S2777, (*_s_dOut_10)[int(4)]);
    DiffPair_float_0 _S2778;
    (&_S2778)->primal_0 = 0.0f;
    (&_S2778)->differential_0 = 0.0f;
    DiffPair_float_0 _S2779;
    (&_S2779)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S2779)->differential_0 = 0.0f;
    _d_max_0(&_S2778, &_S2779, (*_s_dOut_10)[int(4)]);
    DiffPair_float_0 _S2780;
    (&_S2780)->primal_0 = 0.0f;
    (&_S2780)->differential_0 = 0.0f;
    DiffPair_float_0 _S2781;
    (&_S2781)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S2781)->differential_0 = 0.0f;
    _d_max_0(&_S2780, &_S2781, (*_s_dOut_10)[int(4)]);
    DiffPair_float_0 _S2782;
    (&_S2782)->primal_0 = 0.0f;
    (&_S2782)->differential_0 = 0.0f;
    DiffPair_float_0 _S2783;
    (&_S2783)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S2783)->differential_0 = 0.0f;
    _d_max_0(&_S2782, &_S2783, (*_s_dOut_10)[int(3)]);
    DiffPair_float_0 _S2784;
    (&_S2784)->primal_0 = 0.0f;
    (&_S2784)->differential_0 = 0.0f;
    DiffPair_float_0 _S2785;
    (&_S2785)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S2785)->differential_0 = 0.0f;
    _d_max_0(&_S2784, &_S2785, (*_s_dOut_10)[int(3)]);
    DiffPair_float_0 _S2786;
    (&_S2786)->primal_0 = 0.0f;
    (&_S2786)->differential_0 = 0.0f;
    DiffPair_float_0 _S2787;
    (&_S2787)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S2787)->differential_0 = 0.0f;
    _d_max_0(&_S2786, &_S2787, (*_s_dOut_10)[int(3)]);
    DiffPair_float_0 _S2788;
    (&_S2788)->primal_0 = 0.0f;
    (&_S2788)->differential_0 = 0.0f;
    DiffPair_float_0 _S2789;
    (&_S2789)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S2789)->differential_0 = 0.0f;
    _d_max_0(&_S2788, &_S2789, (*_s_dOut_10)[int(2)]);
    DiffPair_float_0 _S2790;
    (&_S2790)->primal_0 = 0.0f;
    (&_S2790)->differential_0 = 0.0f;
    DiffPair_float_0 _S2791;
    (&_S2791)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S2791)->differential_0 = 0.0f;
    _d_max_0(&_S2790, &_S2791, (*_s_dOut_10)[int(2)]);
    DiffPair_float_0 _S2792;
    (&_S2792)->primal_0 = 0.0f;
    (&_S2792)->differential_0 = 0.0f;
    DiffPair_float_0 _S2793;
    (&_S2793)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S2793)->differential_0 = 0.0f;
    _d_max_0(&_S2792, &_S2793, (*_s_dOut_10)[int(2)]);
    float _S2794 = dpparams_2->primal_0[int(12)] * (*_s_dOut_10)[int(1)];
    float _S2795 = dpparams_2->primal_0[int(11)] * (*_s_dOut_10)[int(1)];
    float _S2796 = dpparams_2->primal_0[int(7)] * (*_s_dOut_10)[int(1)];
    float _S2797 = dpparams_2->primal_0[int(6)] * (*_s_dOut_10)[int(1)];
    float _S2798 = dpparams_2->primal_0[int(2)] * (*_s_dOut_10)[int(1)];
    float _S2799 = dpparams_2->primal_0[int(1)] * (*_s_dOut_10)[int(1)];
    PPISPParams_0 _S2800 = PPISPParams_x24_syn_dzero_0();
    (&_S2800)->color_params_1 = _S2735;
    (&_S2800)->exposure_1 = (*_s_dOut_10)[int(0)];
    _S2674 = _S2800;
    (&(&_S2674)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S2801 = _S2692 + _S2697 + _S2800.crf_params_1[int(2)].center_0;
    (&(&_S2674)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S2802 = _S2700 + _S2705 + _S2800.crf_params_1[int(2)].gamma_0;
    (&(&_S2674)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S2803 = _S2708 + _S2713 + _S2800.crf_params_1[int(2)].shoulder_0;
    (&(&_S2674)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S2804 = _S2716 + _S2721 + _S2800.crf_params_1[int(2)].toe_0;
    (&(&_S2674)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S2805 = _S2694 + _S2697 + _S2800.crf_params_1[int(1)].center_0;
    (&(&_S2674)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S2806 = _S2702 + _S2705 + _S2800.crf_params_1[int(1)].gamma_0;
    (&(&_S2674)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S2807 = _S2710 + _S2713 + _S2800.crf_params_1[int(1)].shoulder_0;
    (&(&_S2674)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S2808 = _S2718 + _S2721 + _S2800.crf_params_1[int(1)].toe_0;
    (&(&_S2674)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S2809 = _S2696 + _S2697 + _S2800.crf_params_1[int(0)].center_0;
    (&(&_S2674)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S2810 = _S2704 + _S2705 + _S2800.crf_params_1[int(0)].gamma_0;
    (&(&_S2674)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S2811 = _S2712 + _S2713 + _S2800.crf_params_1[int(0)].shoulder_0;
    (&(&_S2674)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S2812 = _S2720 + _S2721 + _S2800.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S2674)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S2674)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S2674)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S2813 = _S2738 + _S2743 + _S2777.differential_0 + _S2800.vignette_params_1[int(2)].alpha2_0;
    (&(&_S2674)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S2814 = _S2746 + _S2751 + _S2783.differential_0 + _S2800.vignette_params_1[int(2)].alpha1_0;
    (&(&_S2674)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S2815 = _S2754 + _S2759 + _S2789.differential_0 + _S2800.vignette_params_1[int(2)].alpha0_0;
    (&(&_S2674)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S2816 = _S2762 + _S2767 + _S2794 + _S2794 + _S2800.vignette_params_1[int(2)].cy_0;
    (&(&_S2674)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S2817 = _S2770 + _S2775 + _S2795 + _S2795 + _S2800.vignette_params_1[int(2)].cx_0;
    (&(&_S2674)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S2818 = _S2740 + _S2743 + _S2779.differential_0 + _S2800.vignette_params_1[int(1)].alpha2_0;
    (&(&_S2674)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S2819 = _S2748 + _S2751 + _S2785.differential_0 + _S2800.vignette_params_1[int(1)].alpha1_0;
    (&(&_S2674)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S2820 = _S2756 + _S2759 + _S2791.differential_0 + _S2800.vignette_params_1[int(1)].alpha0_0;
    (&(&_S2674)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S2821 = _S2764 + _S2767 + _S2796 + _S2796 + _S2800.vignette_params_1[int(1)].cy_0;
    (&(&_S2674)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S2822 = _S2772 + _S2775 + _S2797 + _S2797 + _S2800.vignette_params_1[int(1)].cx_0;
    (&(&_S2674)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S2823 = _S2742 + _S2743 + _S2781.differential_0 + _S2800.vignette_params_1[int(0)].alpha2_0;
    (&(&_S2674)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S2824 = _S2750 + _S2751 + _S2787.differential_0 + _S2800.vignette_params_1[int(0)].alpha1_0;
    (&(&_S2674)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S2825 = _S2758 + _S2759 + _S2793.differential_0 + _S2800.vignette_params_1[int(0)].alpha0_0;
    (&(&_S2674)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S2826 = _S2766 + _S2767 + _S2798 + _S2798 + _S2800.vignette_params_1[int(0)].cy_0;
    (&(&_S2674)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S2827 = _S2774 + _S2775 + _S2799 + _S2799 + _S2800.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S2828;
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
    _S2828[int(22)] = 0.0f;
    _S2828[int(23)] = 0.0f;
    _S2828[int(24)] = 0.0f;
    _S2828[int(25)] = 0.0f;
    _S2828[int(26)] = 0.0f;
    _S2828[int(27)] = 0.0f;
    _S2828[int(28)] = 0.0f;
    _S2828[int(29)] = 0.0f;
    _S2828[int(30)] = 0.0f;
    _S2828[int(31)] = 0.0f;
    _S2828[int(32)] = 0.0f;
    _S2828[int(33)] = 0.0f;
    _S2828[int(34)] = 0.0f;
    _S2828[int(35)] = 0.0f;
    _S2828[int(8)] = _S2820;
    _S2828[int(16)] = _S2800.color_params_1.b_0.x;
    _S2828[int(15)] = _S2813;
    _S2828[int(14)] = _S2814;
    _S2828[int(13)] = _S2815;
    _S2828[int(12)] = _S2816;
    _S2828[int(11)] = _S2817;
    _S2828[int(10)] = _S2818;
    _S2828[int(9)] = _S2819;
    _S2828[int(17)] = _S2800.color_params_1.b_0.y;
    _S2828[int(7)] = _S2821;
    _S2828[int(6)] = _S2822;
    _S2828[int(5)] = _S2823;
    _S2828[int(4)] = _S2824;
    _S2828[int(3)] = _S2825;
    _S2828[int(2)] = _S2826;
    _S2828[int(1)] = _S2827;
    _S2828[int(0)] = _S2674.exposure_1;
    _S2828[int(26)] = _S2810;
    _S2828[int(34)] = _S2802;
    _S2828[int(33)] = _S2803;
    _S2828[int(32)] = _S2804;
    _S2828[int(31)] = _S2805;
    _S2828[int(30)] = _S2806;
    _S2828[int(29)] = _S2807;
    _S2828[int(28)] = _S2808;
    _S2828[int(27)] = _S2809;
    _S2828[int(35)] = _S2801;
    _S2828[int(25)] = _S2811;
    _S2828[int(24)] = _S2812;
    _S2828[int(23)] = _S2800.color_params_1.n_0.y;
    _S2828[int(22)] = _S2800.color_params_1.n_0.x;
    _S2828[int(21)] = _S2800.color_params_1.g_0.y;
    _S2828[int(20)] = _S2800.color_params_1.g_0.x;
    _S2828[int(19)] = _S2800.color_params_1.r_0.y;
    _S2828[int(18)] = _S2800.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S2828;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S2829, FixedArray<float, 22>  * _S2830)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S2829, _S2830);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_6, FixedArray<float, 22>  grad_out_2, FixedArray<float, 36>  * _S2831)
{
    FixedArray<float, 36>  _S2832 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_6;
    (&dp_params_2)->differential_0 = _S2832;
    FixedArray<float, 22>  _S2833 = grad_out_2;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, &_S2833);
    *_S2831 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_11)
{
    VignettingChannelParams_0 _S2834 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S2835 = {
        _S2834, _S2834, _S2834
    };
    float2  _S2836 = make_float2 (0.0f);
    ColorPPISPParams_0 _S2837 = { _S2836, _S2836, _S2836, _S2836 };
    RQSCRFPPISPChannelParams_0 _S2838 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S2839 = {
        _S2838, _S2838, _S2838
    };
    PPISPParamsRQS_0 _S2840;
    (&_S2840)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S2840)->vignette_params_0 = _S2835;
    (&_S2840)->color_params_0 = _S2837;
    (&_S2840)->crf_params_0 = _S2839;
    (&(&_S2840)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S2840)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S2840)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S2840)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S2840)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S2840)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S2840)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S2840)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S2840)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S2840)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S2840)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S2840)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S2840)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S2840)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S2840)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S2840)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S2840)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S2840)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S2840)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S2840)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S2840)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S2840)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S2840)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S2840)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S2840)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S2840)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S2840)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S2840)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S2840)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S2840)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S2840)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S2840)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S2840)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S2840)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S2840)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S2840)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S2840)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S2840)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_30 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S2841 = dpparams_3->primal_0[int(1)] - mean_30;
    float _S2842 = dpparams_3->primal_0[int(6)] - mean_30;
    float _S2843 = dpparams_3->primal_0[int(11)] - mean_30;
    float mean_31 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S2844 = dpparams_3->primal_0[int(2)] - mean_31;
    float _S2845 = dpparams_3->primal_0[int(7)] - mean_31;
    float _S2846 = dpparams_3->primal_0[int(12)] - mean_31;
    float mean_32 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S2847 = dpparams_3->primal_0[int(3)] - mean_32;
    float _S2848 = dpparams_3->primal_0[int(8)] - mean_32;
    float _S2849 = dpparams_3->primal_0[int(13)] - mean_32;
    float mean_33 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S2850 = dpparams_3->primal_0[int(4)] - mean_33;
    float _S2851 = dpparams_3->primal_0[int(9)] - mean_33;
    float _S2852 = dpparams_3->primal_0[int(14)] - mean_33;
    float mean_34 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S2853 = dpparams_3->primal_0[int(5)] - mean_34;
    float _S2854 = dpparams_3->primal_0[int(10)] - mean_34;
    float _S2855 = dpparams_3->primal_0[int(15)] - mean_34;
    float mean_35 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_36 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_37 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_38 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_39 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S2856 = 0.3333333432674408f * (*_s_dOut_11)[int(22)];
    float _S2857 = (dpparams_3->primal_0[int(38)] - mean_39) * _S2856;
    float _S2858 = _S2857 + _S2857;
    float _S2859 = (dpparams_3->primal_0[int(33)] - mean_39) * _S2856;
    float _S2860 = _S2859 + _S2859;
    float _S2861 = (dpparams_3->primal_0[int(28)] - mean_39) * _S2856;
    float _S2862 = _S2861 + _S2861;
    float _S2863 = 0.3333333432674408f * (- _S2858 + - _S2860 + - _S2862);
    float _S2864 = 0.3333333432674408f * (*_s_dOut_11)[int(21)];
    float _S2865 = (dpparams_3->primal_0[int(37)] - mean_38) * _S2864;
    float _S2866 = _S2865 + _S2865;
    float _S2867 = (dpparams_3->primal_0[int(32)] - mean_38) * _S2864;
    float _S2868 = _S2867 + _S2867;
    float _S2869 = (dpparams_3->primal_0[int(27)] - mean_38) * _S2864;
    float _S2870 = _S2869 + _S2869;
    float _S2871 = 0.3333333432674408f * (- _S2866 + - _S2868 + - _S2870);
    float _S2872 = 0.3333333432674408f * (*_s_dOut_11)[int(20)];
    float _S2873 = (dpparams_3->primal_0[int(36)] - mean_37) * _S2872;
    float _S2874 = _S2873 + _S2873;
    float _S2875 = (dpparams_3->primal_0[int(31)] - mean_37) * _S2872;
    float _S2876 = _S2875 + _S2875;
    float _S2877 = (dpparams_3->primal_0[int(26)] - mean_37) * _S2872;
    float _S2878 = _S2877 + _S2877;
    float _S2879 = 0.3333333432674408f * (- _S2874 + - _S2876 + - _S2878);
    float _S2880 = 0.3333333432674408f * (*_s_dOut_11)[int(19)];
    float _S2881 = (dpparams_3->primal_0[int(35)] - mean_36) * _S2880;
    float _S2882 = _S2881 + _S2881;
    float _S2883 = (dpparams_3->primal_0[int(30)] - mean_36) * _S2880;
    float _S2884 = _S2883 + _S2883;
    float _S2885 = (dpparams_3->primal_0[int(25)] - mean_36) * _S2880;
    float _S2886 = _S2885 + _S2885;
    float _S2887 = 0.3333333432674408f * (- _S2882 + - _S2884 + - _S2886);
    float _S2888 = 0.3333333432674408f * (*_s_dOut_11)[int(18)];
    float _S2889 = (dpparams_3->primal_0[int(34)] - mean_35) * _S2888;
    float _S2890 = _S2889 + _S2889;
    float _S2891 = (dpparams_3->primal_0[int(29)] - mean_35) * _S2888;
    float _S2892 = _S2891 + _S2891;
    float _S2893 = (dpparams_3->primal_0[int(24)] - mean_35) * _S2888;
    float _S2894 = _S2893 + _S2893;
    float _S2895 = 0.3333333432674408f * (- _S2890 + - _S2892 + - _S2894);
    float2  _S2896 = make_float2 ((*_s_dOut_11)[int(16)], (*_s_dOut_11)[int(17)]);
    Matrix<float, 2, 2>  _S2897 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2898;
    (&_S2898)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S2898)->differential_0 = _S2897;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2899;
    (&_S2899)->primal_0 = _S2840.color_params_0.n_0;
    (&_S2899)->differential_0 = _S2836;
    s_bwd_prop_mul_2(&_S2898, &_S2899, _S2896);
    float2  _S2900 = make_float2 ((*_s_dOut_11)[int(14)], (*_s_dOut_11)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2901;
    (&_S2901)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S2901)->differential_0 = _S2897;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2902;
    (&_S2902)->primal_0 = _S2840.color_params_0.g_0;
    (&_S2902)->differential_0 = _S2836;
    s_bwd_prop_mul_2(&_S2901, &_S2902, _S2900);
    float2  _S2903 = make_float2 ((*_s_dOut_11)[int(12)], (*_s_dOut_11)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2904;
    (&_S2904)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S2904)->differential_0 = _S2897;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2905;
    (&_S2905)->primal_0 = _S2840.color_params_0.r_0;
    (&_S2905)->differential_0 = _S2836;
    s_bwd_prop_mul_2(&_S2904, &_S2905, _S2903);
    float2  _S2906 = make_float2 ((*_s_dOut_11)[int(10)], (*_s_dOut_11)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S2907;
    (&_S2907)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S2907)->differential_0 = _S2897;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2908;
    (&_S2908)->primal_0 = _S2840.color_params_0.b_0;
    (&_S2908)->differential_0 = _S2836;
    s_bwd_prop_mul_2(&_S2907, &_S2908, _S2906);
    ColorPPISPParams_0 _S2909 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S2909)->n_0 = _S2899.differential_0;
    (&_S2909)->g_0 = _S2902.differential_0;
    (&_S2909)->r_0 = _S2905.differential_0;
    (&_S2909)->b_0 = _S2908.differential_0;
    float _S2910 = 0.3333333432674408f * (*_s_dOut_11)[int(9)];
    float _S2911 = _S2855 * _S2910;
    float _S2912 = _S2911 + _S2911;
    float _S2913 = _S2854 * _S2910;
    float _S2914 = _S2913 + _S2913;
    float _S2915 = _S2853 * _S2910;
    float _S2916 = _S2915 + _S2915;
    float _S2917 = 0.3333333432674408f * (- _S2912 + - _S2914 + - _S2916);
    float _S2918 = 0.3333333432674408f * (*_s_dOut_11)[int(8)];
    float _S2919 = _S2852 * _S2918;
    float _S2920 = _S2919 + _S2919;
    float _S2921 = _S2851 * _S2918;
    float _S2922 = _S2921 + _S2921;
    float _S2923 = _S2850 * _S2918;
    float _S2924 = _S2923 + _S2923;
    float _S2925 = 0.3333333432674408f * (- _S2920 + - _S2922 + - _S2924);
    float _S2926 = 0.3333333432674408f * (*_s_dOut_11)[int(7)];
    float _S2927 = _S2849 * _S2926;
    float _S2928 = _S2927 + _S2927;
    float _S2929 = _S2848 * _S2926;
    float _S2930 = _S2929 + _S2929;
    float _S2931 = _S2847 * _S2926;
    float _S2932 = _S2931 + _S2931;
    float _S2933 = 0.3333333432674408f * (- _S2928 + - _S2930 + - _S2932);
    float _S2934 = 0.3333333432674408f * (*_s_dOut_11)[int(6)];
    float _S2935 = _S2846 * _S2934;
    float _S2936 = _S2935 + _S2935;
    float _S2937 = _S2845 * _S2934;
    float _S2938 = _S2937 + _S2937;
    float _S2939 = _S2844 * _S2934;
    float _S2940 = _S2939 + _S2939;
    float _S2941 = 0.3333333432674408f * (- _S2936 + - _S2938 + - _S2940);
    float _S2942 = 0.3333333432674408f * (*_s_dOut_11)[int(5)];
    float _S2943 = _S2843 * _S2942;
    float _S2944 = _S2943 + _S2943;
    float _S2945 = _S2842 * _S2942;
    float _S2946 = _S2945 + _S2945;
    float _S2947 = _S2841 * _S2942;
    float _S2948 = _S2947 + _S2947;
    float _S2949 = 0.3333333432674408f * (- _S2944 + - _S2946 + - _S2948);
    DiffPair_float_0 _S2950;
    (&_S2950)->primal_0 = 0.0f;
    (&_S2950)->differential_0 = 0.0f;
    DiffPair_float_0 _S2951;
    (&_S2951)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S2951)->differential_0 = 0.0f;
    _d_max_0(&_S2950, &_S2951, (*_s_dOut_11)[int(4)]);
    DiffPair_float_0 _S2952;
    (&_S2952)->primal_0 = 0.0f;
    (&_S2952)->differential_0 = 0.0f;
    DiffPair_float_0 _S2953;
    (&_S2953)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S2953)->differential_0 = 0.0f;
    _d_max_0(&_S2952, &_S2953, (*_s_dOut_11)[int(4)]);
    DiffPair_float_0 _S2954;
    (&_S2954)->primal_0 = 0.0f;
    (&_S2954)->differential_0 = 0.0f;
    DiffPair_float_0 _S2955;
    (&_S2955)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S2955)->differential_0 = 0.0f;
    _d_max_0(&_S2954, &_S2955, (*_s_dOut_11)[int(4)]);
    DiffPair_float_0 _S2956;
    (&_S2956)->primal_0 = 0.0f;
    (&_S2956)->differential_0 = 0.0f;
    DiffPair_float_0 _S2957;
    (&_S2957)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S2957)->differential_0 = 0.0f;
    _d_max_0(&_S2956, &_S2957, (*_s_dOut_11)[int(3)]);
    DiffPair_float_0 _S2958;
    (&_S2958)->primal_0 = 0.0f;
    (&_S2958)->differential_0 = 0.0f;
    DiffPair_float_0 _S2959;
    (&_S2959)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S2959)->differential_0 = 0.0f;
    _d_max_0(&_S2958, &_S2959, (*_s_dOut_11)[int(3)]);
    DiffPair_float_0 _S2960;
    (&_S2960)->primal_0 = 0.0f;
    (&_S2960)->differential_0 = 0.0f;
    DiffPair_float_0 _S2961;
    (&_S2961)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S2961)->differential_0 = 0.0f;
    _d_max_0(&_S2960, &_S2961, (*_s_dOut_11)[int(3)]);
    DiffPair_float_0 _S2962;
    (&_S2962)->primal_0 = 0.0f;
    (&_S2962)->differential_0 = 0.0f;
    DiffPair_float_0 _S2963;
    (&_S2963)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S2963)->differential_0 = 0.0f;
    _d_max_0(&_S2962, &_S2963, (*_s_dOut_11)[int(2)]);
    DiffPair_float_0 _S2964;
    (&_S2964)->primal_0 = 0.0f;
    (&_S2964)->differential_0 = 0.0f;
    DiffPair_float_0 _S2965;
    (&_S2965)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S2965)->differential_0 = 0.0f;
    _d_max_0(&_S2964, &_S2965, (*_s_dOut_11)[int(2)]);
    DiffPair_float_0 _S2966;
    (&_S2966)->primal_0 = 0.0f;
    (&_S2966)->differential_0 = 0.0f;
    DiffPair_float_0 _S2967;
    (&_S2967)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S2967)->differential_0 = 0.0f;
    _d_max_0(&_S2966, &_S2967, (*_s_dOut_11)[int(2)]);
    float _S2968 = dpparams_3->primal_0[int(12)] * (*_s_dOut_11)[int(1)];
    float _S2969 = dpparams_3->primal_0[int(11)] * (*_s_dOut_11)[int(1)];
    float _S2970 = dpparams_3->primal_0[int(7)] * (*_s_dOut_11)[int(1)];
    float _S2971 = dpparams_3->primal_0[int(6)] * (*_s_dOut_11)[int(1)];
    float _S2972 = dpparams_3->primal_0[int(2)] * (*_s_dOut_11)[int(1)];
    float _S2973 = dpparams_3->primal_0[int(1)] * (*_s_dOut_11)[int(1)];
    PPISPParamsRQS_0 _S2974 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S2974)->color_params_0 = _S2909;
    (&_S2974)->exposure_0 = (*_s_dOut_11)[int(0)];
    _S2840 = _S2974;
    (&(&_S2840)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S2975 = _S2858 + _S2863 + _S2974.crf_params_0[int(2)].gc_0;
    (&(&_S2840)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S2976 = _S2866 + _S2871 + _S2974.crf_params_0[int(2)].y0_0;
    (&(&_S2840)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S2977 = _S2874 + _S2879 + _S2974.crf_params_0[int(2)].x0_0;
    (&(&_S2840)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S2978 = _S2882 + _S2887 + _S2974.crf_params_0[int(2)].g1_0;
    (&(&_S2840)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S2979 = _S2890 + _S2895 + _S2974.crf_params_0[int(2)].g0_0;
    (&(&_S2840)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S2980 = _S2860 + _S2863 + _S2974.crf_params_0[int(1)].gc_0;
    (&(&_S2840)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S2981 = _S2868 + _S2871 + _S2974.crf_params_0[int(1)].y0_0;
    (&(&_S2840)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S2982 = _S2876 + _S2879 + _S2974.crf_params_0[int(1)].x0_0;
    (&(&_S2840)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S2983 = _S2884 + _S2887 + _S2974.crf_params_0[int(1)].g1_0;
    (&(&_S2840)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S2984 = _S2892 + _S2895 + _S2974.crf_params_0[int(1)].g0_0;
    (&(&_S2840)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S2985 = _S2862 + _S2863 + _S2974.crf_params_0[int(0)].gc_0;
    (&(&_S2840)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S2986 = _S2870 + _S2871 + _S2974.crf_params_0[int(0)].y0_0;
    (&(&_S2840)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S2987 = _S2878 + _S2879 + _S2974.crf_params_0[int(0)].x0_0;
    (&(&_S2840)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S2988 = _S2886 + _S2887 + _S2974.crf_params_0[int(0)].g1_0;
    (&(&_S2840)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S2989 = _S2894 + _S2895 + _S2974.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S2840)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S2840)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S2840)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S2990 = _S2912 + _S2917 + _S2951.differential_0 + _S2974.vignette_params_0[int(2)].alpha2_0;
    (&(&_S2840)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S2991 = _S2920 + _S2925 + _S2957.differential_0 + _S2974.vignette_params_0[int(2)].alpha1_0;
    (&(&_S2840)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S2992 = _S2928 + _S2933 + _S2963.differential_0 + _S2974.vignette_params_0[int(2)].alpha0_0;
    (&(&_S2840)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S2993 = _S2936 + _S2941 + _S2968 + _S2968 + _S2974.vignette_params_0[int(2)].cy_0;
    (&(&_S2840)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S2994 = _S2944 + _S2949 + _S2969 + _S2969 + _S2974.vignette_params_0[int(2)].cx_0;
    (&(&_S2840)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S2995 = _S2914 + _S2917 + _S2953.differential_0 + _S2974.vignette_params_0[int(1)].alpha2_0;
    (&(&_S2840)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S2996 = _S2922 + _S2925 + _S2959.differential_0 + _S2974.vignette_params_0[int(1)].alpha1_0;
    (&(&_S2840)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S2997 = _S2930 + _S2933 + _S2965.differential_0 + _S2974.vignette_params_0[int(1)].alpha0_0;
    (&(&_S2840)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S2998 = _S2938 + _S2941 + _S2970 + _S2970 + _S2974.vignette_params_0[int(1)].cy_0;
    (&(&_S2840)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S2999 = _S2946 + _S2949 + _S2971 + _S2971 + _S2974.vignette_params_0[int(1)].cx_0;
    (&(&_S2840)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S3000 = _S2916 + _S2917 + _S2955.differential_0 + _S2974.vignette_params_0[int(0)].alpha2_0;
    (&(&_S2840)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S3001 = _S2924 + _S2925 + _S2961.differential_0 + _S2974.vignette_params_0[int(0)].alpha1_0;
    (&(&_S2840)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S3002 = _S2932 + _S2933 + _S2967.differential_0 + _S2974.vignette_params_0[int(0)].alpha0_0;
    (&(&_S2840)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S3003 = _S2940 + _S2941 + _S2972 + _S2972 + _S2974.vignette_params_0[int(0)].cy_0;
    (&(&_S2840)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S3004 = _S2948 + _S2949 + _S2973 + _S2973 + _S2974.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S3005;
    _S3005[int(0)] = 0.0f;
    _S3005[int(1)] = 0.0f;
    _S3005[int(2)] = 0.0f;
    _S3005[int(3)] = 0.0f;
    _S3005[int(4)] = 0.0f;
    _S3005[int(5)] = 0.0f;
    _S3005[int(6)] = 0.0f;
    _S3005[int(7)] = 0.0f;
    _S3005[int(8)] = 0.0f;
    _S3005[int(9)] = 0.0f;
    _S3005[int(10)] = 0.0f;
    _S3005[int(11)] = 0.0f;
    _S3005[int(12)] = 0.0f;
    _S3005[int(13)] = 0.0f;
    _S3005[int(14)] = 0.0f;
    _S3005[int(15)] = 0.0f;
    _S3005[int(16)] = 0.0f;
    _S3005[int(17)] = 0.0f;
    _S3005[int(18)] = 0.0f;
    _S3005[int(19)] = 0.0f;
    _S3005[int(20)] = 0.0f;
    _S3005[int(21)] = 0.0f;
    _S3005[int(22)] = 0.0f;
    _S3005[int(23)] = 0.0f;
    _S3005[int(24)] = 0.0f;
    _S3005[int(25)] = 0.0f;
    _S3005[int(26)] = 0.0f;
    _S3005[int(27)] = 0.0f;
    _S3005[int(28)] = 0.0f;
    _S3005[int(29)] = 0.0f;
    _S3005[int(30)] = 0.0f;
    _S3005[int(31)] = 0.0f;
    _S3005[int(32)] = 0.0f;
    _S3005[int(33)] = 0.0f;
    _S3005[int(34)] = 0.0f;
    _S3005[int(35)] = 0.0f;
    _S3005[int(36)] = 0.0f;
    _S3005[int(37)] = 0.0f;
    _S3005[int(38)] = 0.0f;
    _S3005[int(9)] = _S2996;
    _S3005[int(18)] = _S2974.color_params_0.r_0.x;
    _S3005[int(17)] = _S2974.color_params_0.b_0.y;
    _S3005[int(16)] = _S2974.color_params_0.b_0.x;
    _S3005[int(15)] = _S2990;
    _S3005[int(14)] = _S2991;
    _S3005[int(13)] = _S2992;
    _S3005[int(12)] = _S2993;
    _S3005[int(11)] = _S2994;
    _S3005[int(10)] = _S2995;
    _S3005[int(19)] = _S2974.color_params_0.r_0.y;
    _S3005[int(8)] = _S2997;
    _S3005[int(7)] = _S2998;
    _S3005[int(6)] = _S2999;
    _S3005[int(5)] = _S3000;
    _S3005[int(4)] = _S3001;
    _S3005[int(3)] = _S3002;
    _S3005[int(2)] = _S3003;
    _S3005[int(1)] = _S3004;
    _S3005[int(0)] = _S2840.exposure_0;
    _S3005[int(28)] = _S2985;
    _S3005[int(37)] = _S2976;
    _S3005[int(36)] = _S2977;
    _S3005[int(35)] = _S2978;
    _S3005[int(34)] = _S2979;
    _S3005[int(33)] = _S2980;
    _S3005[int(32)] = _S2981;
    _S3005[int(31)] = _S2982;
    _S3005[int(30)] = _S2983;
    _S3005[int(29)] = _S2984;
    _S3005[int(38)] = _S2975;
    _S3005[int(27)] = _S2986;
    _S3005[int(26)] = _S2987;
    _S3005[int(25)] = _S2988;
    _S3005[int(24)] = _S2989;
    _S3005[int(23)] = _S2974.color_params_0.n_0.y;
    _S3005[int(22)] = _S2974.color_params_0.n_0.x;
    _S3005[int(21)] = _S2974.color_params_0.g_0.y;
    _S3005[int(20)] = _S2974.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S3005;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S3006, FixedArray<float, 23>  * _S3007)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S3006, _S3007);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_7, FixedArray<float, 23>  grad_out_3, FixedArray<float, 39>  * _S3008)
{
    FixedArray<float, 39>  _S3009 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_7;
    (&dp_params_3)->differential_0 = _S3009;
    FixedArray<float, 23>  _S3010 = grad_out_3;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, &_S3010);
    *_S3008 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_2, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S3011)
{
    float _S3012;
    FixedArray<float, 6>  losses_6;
    float _S3013 = float(num_cameras_0);
    float _S3014 = raw_losses_2[int(0)] / _S3013;
    for(;;)
    {
        float _S3015 = (F32_abs((_S3014)));
        if(_S3015 < 0.10000000149011612f)
        {
            _S3012 = 0.5f * _S3014 * _S3014 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S3012 = _S3015 - 0.05000000074505806f;
            break;
        }
    }
    losses_6[int(0)] = _S3012;
    losses_6[int(1)] = raw_losses_2[int(1)] / (3.0f * _S3013);
    losses_6[int(2)] = (raw_losses_2[int(2)] + raw_losses_2[int(3)] + raw_losses_2[int(4)]) / (9.0f * _S3013);
    losses_6[int(3)] = (raw_losses_2[int(5)] + raw_losses_2[int(6)] + raw_losses_2[int(7)] + raw_losses_2[int(8)] + raw_losses_2[int(9)]) / (5.0f * _S3013);
    float _S3016 = raw_losses_2[int(10)] / _S3013;
    for(;;)
    {
        float _S3017 = (F32_abs((_S3016)));
        if(_S3017 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3016 * _S3016 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3017 - 0.00249999994412065f;
            break;
        }
    }
    float _S3018;
    float _S3019 = raw_losses_2[int(11)] / _S3013;
    for(;;)
    {
        float _S3020 = (F32_abs((_S3019)));
        if(_S3020 < 0.00499999988824129f)
        {
            _S3018 = 0.5f * _S3019 * _S3019 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3018 = _S3020 - 0.00249999994412065f;
            break;
        }
    }
    float _S3021 = _S3012 + _S3018;
    float _S3022 = raw_losses_2[int(12)] / _S3013;
    for(;;)
    {
        float _S3023 = (F32_abs((_S3022)));
        if(_S3023 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3022 * _S3022 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3023 - 0.00249999994412065f;
            break;
        }
    }
    float _S3024 = _S3021 + _S3012;
    float _S3025 = raw_losses_2[int(13)] / _S3013;
    for(;;)
    {
        float _S3026 = (F32_abs((_S3025)));
        if(_S3026 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3025 * _S3025 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3026 - 0.00249999994412065f;
            break;
        }
    }
    float _S3027 = _S3024 + _S3012;
    float _S3028 = raw_losses_2[int(14)] / _S3013;
    for(;;)
    {
        float _S3029 = (F32_abs((_S3028)));
        if(_S3029 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3028 * _S3028 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3029 - 0.00249999994412065f;
            break;
        }
    }
    float _S3030 = _S3027 + _S3012;
    float _S3031 = raw_losses_2[int(15)] / _S3013;
    for(;;)
    {
        float _S3032 = (F32_abs((_S3031)));
        if(_S3032 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3031 * _S3031 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3032 - 0.00249999994412065f;
            break;
        }
    }
    float _S3033 = _S3030 + _S3012;
    float _S3034 = raw_losses_2[int(16)] / _S3013;
    for(;;)
    {
        float _S3035 = (F32_abs((_S3034)));
        if(_S3035 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3034 * _S3034 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3035 - 0.00249999994412065f;
            break;
        }
    }
    float _S3036 = _S3033 + _S3012;
    float _S3037 = raw_losses_2[int(17)] / _S3013;
    for(;;)
    {
        float _S3038 = (F32_abs((_S3037)));
        if(_S3038 < 0.00499999988824129f)
        {
            _S3012 = 0.5f * _S3037 * _S3037 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3012 = _S3038 - 0.00249999994412065f;
            break;
        }
    }
    float _S3039 = (_S3036 + _S3012) / 8.0f;
    float _S3040 = (raw_losses_2[int(18)] + raw_losses_2[int(19)] + raw_losses_2[int(20)] + raw_losses_2[int(21)]) / (4.0f * _S3013);
    losses_6[int(0)] = losses_6[int(0)] * loss_weights_0[int(0)];
    losses_6[int(1)] = losses_6[int(1)] * loss_weights_0[int(1)];
    losses_6[int(2)] = losses_6[int(2)] * loss_weights_0[int(2)];
    losses_6[int(3)] = losses_6[int(3)] * loss_weights_0[int(3)];
    losses_6[int(4)] = _S3039 * loss_weights_0[int(4)];
    losses_6[int(5)] = _S3040 * loss_weights_0[int(5)];
    *_S3011 = losses_6;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_3, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S3041)
{
    float _S3042;
    FixedArray<float, 6>  losses_7;
    float _S3043 = float(num_cameras_1);
    float _S3044 = raw_losses_3[int(0)] / _S3043;
    for(;;)
    {
        float _S3045 = (F32_abs((_S3044)));
        if(_S3045 < 0.10000000149011612f)
        {
            _S3042 = 0.5f * _S3044 * _S3044 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S3042 = _S3045 - 0.05000000074505806f;
            break;
        }
    }
    losses_7[int(0)] = _S3042;
    losses_7[int(1)] = raw_losses_3[int(1)] / (3.0f * _S3043);
    losses_7[int(2)] = (raw_losses_3[int(2)] + raw_losses_3[int(3)] + raw_losses_3[int(4)]) / (9.0f * _S3043);
    float _S3046 = 5.0f * _S3043;
    losses_7[int(3)] = (raw_losses_3[int(5)] + raw_losses_3[int(6)] + raw_losses_3[int(7)] + raw_losses_3[int(8)] + raw_losses_3[int(9)]) / _S3046;
    float _S3047 = raw_losses_3[int(10)] / _S3043;
    for(;;)
    {
        float _S3048 = (F32_abs((_S3047)));
        if(_S3048 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3047 * _S3047 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3048 - 0.00249999994412065f;
            break;
        }
    }
    float _S3049;
    float _S3050 = raw_losses_3[int(11)] / _S3043;
    for(;;)
    {
        float _S3051 = (F32_abs((_S3050)));
        if(_S3051 < 0.00499999988824129f)
        {
            _S3049 = 0.5f * _S3050 * _S3050 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3049 = _S3051 - 0.00249999994412065f;
            break;
        }
    }
    float _S3052 = _S3042 + _S3049;
    float _S3053 = raw_losses_3[int(12)] / _S3043;
    for(;;)
    {
        float _S3054 = (F32_abs((_S3053)));
        if(_S3054 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3053 * _S3053 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3054 - 0.00249999994412065f;
            break;
        }
    }
    float _S3055 = _S3052 + _S3042;
    float _S3056 = raw_losses_3[int(13)] / _S3043;
    for(;;)
    {
        float _S3057 = (F32_abs((_S3056)));
        if(_S3057 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3056 * _S3056 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3057 - 0.00249999994412065f;
            break;
        }
    }
    float _S3058 = _S3055 + _S3042;
    float _S3059 = raw_losses_3[int(14)] / _S3043;
    for(;;)
    {
        float _S3060 = (F32_abs((_S3059)));
        if(_S3060 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3059 * _S3059 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3060 - 0.00249999994412065f;
            break;
        }
    }
    float _S3061 = _S3058 + _S3042;
    float _S3062 = raw_losses_3[int(15)] / _S3043;
    for(;;)
    {
        float _S3063 = (F32_abs((_S3062)));
        if(_S3063 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3062 * _S3062 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3063 - 0.00249999994412065f;
            break;
        }
    }
    float _S3064 = _S3061 + _S3042;
    float _S3065 = raw_losses_3[int(16)] / _S3043;
    for(;;)
    {
        float _S3066 = (F32_abs((_S3065)));
        if(_S3066 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3065 * _S3065 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3066 - 0.00249999994412065f;
            break;
        }
    }
    float _S3067 = _S3064 + _S3042;
    float _S3068 = raw_losses_3[int(17)] / _S3043;
    for(;;)
    {
        float _S3069 = (F32_abs((_S3068)));
        if(_S3069 < 0.00499999988824129f)
        {
            _S3042 = 0.5f * _S3068 * _S3068 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S3042 = _S3069 - 0.00249999994412065f;
            break;
        }
    }
    float _S3070 = (_S3067 + _S3042) / 8.0f;
    float _S3071 = (raw_losses_3[int(18)] + raw_losses_3[int(19)] + raw_losses_3[int(20)] + raw_losses_3[int(21)] + raw_losses_3[int(22)]) / _S3046;
    losses_7[int(0)] = losses_7[int(0)] * loss_weights_1[int(0)];
    losses_7[int(1)] = losses_7[int(1)] * loss_weights_1[int(1)];
    losses_7[int(2)] = losses_7[int(2)] * loss_weights_1[int(2)];
    losses_7[int(3)] = losses_7[int(3)] * loss_weights_1[int(3)];
    losses_7[int(4)] = _S3070 * loss_weights_1[int(4)];
    losses_7[int(5)] = _S3071 * loss_weights_1[int(5)];
    *_S3041 = losses_7;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_1, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_12)
{
    FixedArray<float, 22>  _S3072 = dpraw_losses_1->primal_0;
    float _S3073 = float(num_cameras_2);
    float _S3074 = dpraw_losses_1->primal_0[int(0)] / _S3073;
    bool _S3075 = (s_primal_ctx_abs_0(_S3074)) < 0.10000000149011612f;
    float _S3076;
    if(_S3075)
    {
        _S3076 = 0.5f * _S3074;
    }
    else
    {
        _S3076 = 0.0f;
    }
    float _S3077 = 3.0f * _S3073;
    float _S3078 = 9.0f * _S3073;
    float _S3079 = 5.0f * _S3073;
    float _S3080 = _S3072[int(10)] / _S3073;
    bool _S3081 = (s_primal_ctx_abs_0(_S3080)) < 0.00499999988824129f;
    float _S3082;
    if(_S3081)
    {
        _S3082 = 0.5f * _S3080;
    }
    else
    {
        _S3082 = 0.0f;
    }
    float _S3083 = _S3072[int(11)] / _S3073;
    bool _S3084 = (s_primal_ctx_abs_0(_S3083)) < 0.00499999988824129f;
    float _S3085;
    if(_S3084)
    {
        _S3085 = 0.5f * _S3083;
    }
    else
    {
        _S3085 = 0.0f;
    }
    float _S3086 = _S3072[int(12)] / _S3073;
    bool _S3087 = (s_primal_ctx_abs_0(_S3086)) < 0.00499999988824129f;
    float _S3088;
    if(_S3087)
    {
        _S3088 = 0.5f * _S3086;
    }
    else
    {
        _S3088 = 0.0f;
    }
    float _S3089 = _S3072[int(13)] / _S3073;
    bool _S3090 = (s_primal_ctx_abs_0(_S3089)) < 0.00499999988824129f;
    float _S3091;
    if(_S3090)
    {
        _S3091 = 0.5f * _S3089;
    }
    else
    {
        _S3091 = 0.0f;
    }
    float _S3092 = _S3072[int(14)] / _S3073;
    bool _S3093 = (s_primal_ctx_abs_0(_S3092)) < 0.00499999988824129f;
    float _S3094;
    if(_S3093)
    {
        _S3094 = 0.5f * _S3092;
    }
    else
    {
        _S3094 = 0.0f;
    }
    float _S3095 = _S3072[int(15)] / _S3073;
    bool _S3096 = (s_primal_ctx_abs_0(_S3095)) < 0.00499999988824129f;
    float _S3097;
    if(_S3096)
    {
        _S3097 = 0.5f * _S3095;
    }
    else
    {
        _S3097 = 0.0f;
    }
    float _S3098 = _S3072[int(16)] / _S3073;
    bool _S3099 = (s_primal_ctx_abs_0(_S3098)) < 0.00499999988824129f;
    float _S3100;
    if(_S3099)
    {
        _S3100 = 0.5f * _S3098;
    }
    else
    {
        _S3100 = 0.0f;
    }
    float _S3101 = _S3072[int(17)] / _S3073;
    bool _S3102 = (s_primal_ctx_abs_0(_S3101)) < 0.00499999988824129f;
    float _S3103;
    if(_S3102)
    {
        _S3103 = 0.5f * _S3101;
    }
    else
    {
        _S3103 = 0.0f;
    }
    float _S3104 = (*loss_weights_2)[int(3)] * (*_s_dOut_12)[int(3)];
    float _S3105 = (*loss_weights_2)[int(2)] * (*_s_dOut_12)[int(2)];
    float _S3106 = (*loss_weights_2)[int(1)] * (*_s_dOut_12)[int(1)];
    float _S3107 = (*loss_weights_2)[int(0)] * (*_s_dOut_12)[int(0)];
    float _S3108 = (*loss_weights_2)[int(5)] * (*_s_dOut_12)[int(5)] / (4.0f * _S3073);
    float _S3109 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_12)[int(4)]);
    FixedArray<float, 22>  _S3110;
    _S3110[int(0)] = 0.0f;
    _S3110[int(1)] = 0.0f;
    _S3110[int(2)] = 0.0f;
    _S3110[int(3)] = 0.0f;
    _S3110[int(4)] = 0.0f;
    _S3110[int(5)] = 0.0f;
    _S3110[int(6)] = 0.0f;
    _S3110[int(7)] = 0.0f;
    _S3110[int(8)] = 0.0f;
    _S3110[int(9)] = 0.0f;
    _S3110[int(10)] = 0.0f;
    _S3110[int(11)] = 0.0f;
    _S3110[int(12)] = 0.0f;
    _S3110[int(13)] = 0.0f;
    _S3110[int(14)] = 0.0f;
    _S3110[int(15)] = 0.0f;
    _S3110[int(16)] = 0.0f;
    _S3110[int(17)] = 0.0f;
    _S3110[int(18)] = 0.0f;
    _S3110[int(19)] = 0.0f;
    _S3110[int(20)] = 0.0f;
    _S3110[int(21)] = 0.0f;
    _S3110[int(21)] = _S3108;
    _S3110[int(20)] = _S3108;
    _S3110[int(19)] = _S3108;
    _S3110[int(18)] = _S3108;
    float _S3111 = _S3110[int(0)];
    float _S3112 = _S3110[int(1)];
    float _S3113 = _S3110[int(2)];
    float _S3114 = _S3110[int(3)];
    float _S3115 = _S3110[int(4)];
    float _S3116 = _S3110[int(5)];
    float _S3117 = _S3110[int(6)];
    float _S3118 = _S3110[int(7)];
    float _S3119 = _S3110[int(8)];
    float _S3120 = _S3110[int(9)];
    float _S3121 = _S3110[int(10)];
    float _S3122 = _S3110[int(11)];
    float _S3123 = _S3110[int(12)];
    float _S3124 = _S3110[int(13)];
    float _S3125 = _S3110[int(14)];
    float _S3126 = _S3110[int(15)];
    float _S3127 = _S3110[int(16)];
    float _S3128 = _S3110[int(17)];
    float _S3129 = _S3110[int(18)];
    float _S3130 = _S3110[int(19)];
    float _S3131 = _S3110[int(20)];
    float _S3132 = _S3110[int(21)];
    float _S3133;
    if(_S3102)
    {
        float _S3134 = 200.0f * _S3109;
        float _S3135 = _S3103 * _S3134 + 0.5f * (_S3101 * _S3134);
        _S3103 = 0.0f;
        _S3133 = _S3135;
    }
    else
    {
        _S3103 = _S3109;
        _S3133 = 0.0f;
    }
    DiffPair_float_0 _S3136;
    (&_S3136)->primal_0 = _S3101;
    (&_S3136)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3136, _S3103);
    float _S3137 = (_S3136.differential_0 + _S3133) / _S3073;
    FixedArray<float, 22>  _S3138;
    _S3138[int(0)] = 0.0f;
    _S3138[int(1)] = 0.0f;
    _S3138[int(2)] = 0.0f;
    _S3138[int(3)] = 0.0f;
    _S3138[int(4)] = 0.0f;
    _S3138[int(5)] = 0.0f;
    _S3138[int(6)] = 0.0f;
    _S3138[int(7)] = 0.0f;
    _S3138[int(8)] = 0.0f;
    _S3138[int(9)] = 0.0f;
    _S3138[int(10)] = 0.0f;
    _S3138[int(11)] = 0.0f;
    _S3138[int(12)] = 0.0f;
    _S3138[int(13)] = 0.0f;
    _S3138[int(14)] = 0.0f;
    _S3138[int(15)] = 0.0f;
    _S3138[int(16)] = 0.0f;
    _S3138[int(17)] = 0.0f;
    _S3138[int(18)] = 0.0f;
    _S3138[int(19)] = 0.0f;
    _S3138[int(20)] = 0.0f;
    _S3138[int(21)] = 0.0f;
    _S3138[int(17)] = _S3137;
    float _S3139 = _S3111 + _S3138[int(0)];
    float _S3140 = _S3112 + _S3138[int(1)];
    float _S3141 = _S3113 + _S3138[int(2)];
    float _S3142 = _S3114 + _S3138[int(3)];
    float _S3143 = _S3115 + _S3138[int(4)];
    float _S3144 = _S3116 + _S3138[int(5)];
    float _S3145 = _S3117 + _S3138[int(6)];
    float _S3146 = _S3118 + _S3138[int(7)];
    float _S3147 = _S3119 + _S3138[int(8)];
    float _S3148 = _S3120 + _S3138[int(9)];
    float _S3149 = _S3121 + _S3138[int(10)];
    float _S3150 = _S3122 + _S3138[int(11)];
    float _S3151 = _S3123 + _S3138[int(12)];
    float _S3152 = _S3124 + _S3138[int(13)];
    float _S3153 = _S3125 + _S3138[int(14)];
    float _S3154 = _S3126 + _S3138[int(15)];
    float _S3155 = _S3127 + _S3138[int(16)];
    float _S3156 = _S3128 + _S3138[int(17)];
    float _S3157 = _S3129 + _S3138[int(18)];
    float _S3158 = _S3130 + _S3138[int(19)];
    float _S3159 = _S3131 + _S3138[int(20)];
    float _S3160 = _S3132 + _S3138[int(21)];
    if(_S3099)
    {
        float _S3161 = 200.0f * _S3109;
        float _S3162 = _S3100 * _S3161 + 0.5f * (_S3098 * _S3161);
        _S3100 = 0.0f;
        _S3103 = _S3162;
    }
    else
    {
        _S3100 = _S3109;
        _S3103 = 0.0f;
    }
    DiffPair_float_0 _S3163;
    (&_S3163)->primal_0 = _S3098;
    (&_S3163)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3163, _S3100);
    float _S3164 = (_S3163.differential_0 + _S3103) / _S3073;
    FixedArray<float, 22>  _S3165;
    _S3165[int(0)] = 0.0f;
    _S3165[int(1)] = 0.0f;
    _S3165[int(2)] = 0.0f;
    _S3165[int(3)] = 0.0f;
    _S3165[int(4)] = 0.0f;
    _S3165[int(5)] = 0.0f;
    _S3165[int(6)] = 0.0f;
    _S3165[int(7)] = 0.0f;
    _S3165[int(8)] = 0.0f;
    _S3165[int(9)] = 0.0f;
    _S3165[int(10)] = 0.0f;
    _S3165[int(11)] = 0.0f;
    _S3165[int(12)] = 0.0f;
    _S3165[int(13)] = 0.0f;
    _S3165[int(14)] = 0.0f;
    _S3165[int(15)] = 0.0f;
    _S3165[int(16)] = 0.0f;
    _S3165[int(17)] = 0.0f;
    _S3165[int(18)] = 0.0f;
    _S3165[int(19)] = 0.0f;
    _S3165[int(20)] = 0.0f;
    _S3165[int(21)] = 0.0f;
    _S3165[int(16)] = _S3164;
    float _S3166 = _S3139 + _S3165[int(0)];
    float _S3167 = _S3140 + _S3165[int(1)];
    float _S3168 = _S3141 + _S3165[int(2)];
    float _S3169 = _S3142 + _S3165[int(3)];
    float _S3170 = _S3143 + _S3165[int(4)];
    float _S3171 = _S3144 + _S3165[int(5)];
    float _S3172 = _S3145 + _S3165[int(6)];
    float _S3173 = _S3146 + _S3165[int(7)];
    float _S3174 = _S3147 + _S3165[int(8)];
    float _S3175 = _S3148 + _S3165[int(9)];
    float _S3176 = _S3149 + _S3165[int(10)];
    float _S3177 = _S3150 + _S3165[int(11)];
    float _S3178 = _S3151 + _S3165[int(12)];
    float _S3179 = _S3152 + _S3165[int(13)];
    float _S3180 = _S3153 + _S3165[int(14)];
    float _S3181 = _S3154 + _S3165[int(15)];
    float _S3182 = _S3155 + _S3165[int(16)];
    float _S3183 = _S3156 + _S3165[int(17)];
    float _S3184 = _S3157 + _S3165[int(18)];
    float _S3185 = _S3158 + _S3165[int(19)];
    float _S3186 = _S3159 + _S3165[int(20)];
    float _S3187 = _S3160 + _S3165[int(21)];
    if(_S3096)
    {
        float _S3188 = 200.0f * _S3109;
        float _S3189 = _S3097 * _S3188 + 0.5f * (_S3095 * _S3188);
        _S3097 = 0.0f;
        _S3100 = _S3189;
    }
    else
    {
        _S3097 = _S3109;
        _S3100 = 0.0f;
    }
    DiffPair_float_0 _S3190;
    (&_S3190)->primal_0 = _S3095;
    (&_S3190)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3190, _S3097);
    float _S3191 = (_S3190.differential_0 + _S3100) / _S3073;
    FixedArray<float, 22>  _S3192;
    _S3192[int(0)] = 0.0f;
    _S3192[int(1)] = 0.0f;
    _S3192[int(2)] = 0.0f;
    _S3192[int(3)] = 0.0f;
    _S3192[int(4)] = 0.0f;
    _S3192[int(5)] = 0.0f;
    _S3192[int(6)] = 0.0f;
    _S3192[int(7)] = 0.0f;
    _S3192[int(8)] = 0.0f;
    _S3192[int(9)] = 0.0f;
    _S3192[int(10)] = 0.0f;
    _S3192[int(11)] = 0.0f;
    _S3192[int(12)] = 0.0f;
    _S3192[int(13)] = 0.0f;
    _S3192[int(14)] = 0.0f;
    _S3192[int(15)] = 0.0f;
    _S3192[int(16)] = 0.0f;
    _S3192[int(17)] = 0.0f;
    _S3192[int(18)] = 0.0f;
    _S3192[int(19)] = 0.0f;
    _S3192[int(20)] = 0.0f;
    _S3192[int(21)] = 0.0f;
    _S3192[int(15)] = _S3191;
    float _S3193 = _S3166 + _S3192[int(0)];
    float _S3194 = _S3167 + _S3192[int(1)];
    float _S3195 = _S3168 + _S3192[int(2)];
    float _S3196 = _S3169 + _S3192[int(3)];
    float _S3197 = _S3170 + _S3192[int(4)];
    float _S3198 = _S3171 + _S3192[int(5)];
    float _S3199 = _S3172 + _S3192[int(6)];
    float _S3200 = _S3173 + _S3192[int(7)];
    float _S3201 = _S3174 + _S3192[int(8)];
    float _S3202 = _S3175 + _S3192[int(9)];
    float _S3203 = _S3176 + _S3192[int(10)];
    float _S3204 = _S3177 + _S3192[int(11)];
    float _S3205 = _S3178 + _S3192[int(12)];
    float _S3206 = _S3179 + _S3192[int(13)];
    float _S3207 = _S3180 + _S3192[int(14)];
    float _S3208 = _S3181 + _S3192[int(15)];
    float _S3209 = _S3182 + _S3192[int(16)];
    float _S3210 = _S3183 + _S3192[int(17)];
    float _S3211 = _S3184 + _S3192[int(18)];
    float _S3212 = _S3185 + _S3192[int(19)];
    float _S3213 = _S3186 + _S3192[int(20)];
    float _S3214 = _S3187 + _S3192[int(21)];
    if(_S3093)
    {
        float _S3215 = 200.0f * _S3109;
        float _S3216 = _S3094 * _S3215 + 0.5f * (_S3092 * _S3215);
        _S3094 = 0.0f;
        _S3097 = _S3216;
    }
    else
    {
        _S3094 = _S3109;
        _S3097 = 0.0f;
    }
    DiffPair_float_0 _S3217;
    (&_S3217)->primal_0 = _S3092;
    (&_S3217)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3217, _S3094);
    float _S3218 = (_S3217.differential_0 + _S3097) / _S3073;
    FixedArray<float, 22>  _S3219;
    _S3219[int(0)] = 0.0f;
    _S3219[int(1)] = 0.0f;
    _S3219[int(2)] = 0.0f;
    _S3219[int(3)] = 0.0f;
    _S3219[int(4)] = 0.0f;
    _S3219[int(5)] = 0.0f;
    _S3219[int(6)] = 0.0f;
    _S3219[int(7)] = 0.0f;
    _S3219[int(8)] = 0.0f;
    _S3219[int(9)] = 0.0f;
    _S3219[int(10)] = 0.0f;
    _S3219[int(11)] = 0.0f;
    _S3219[int(12)] = 0.0f;
    _S3219[int(13)] = 0.0f;
    _S3219[int(14)] = 0.0f;
    _S3219[int(15)] = 0.0f;
    _S3219[int(16)] = 0.0f;
    _S3219[int(17)] = 0.0f;
    _S3219[int(18)] = 0.0f;
    _S3219[int(19)] = 0.0f;
    _S3219[int(20)] = 0.0f;
    _S3219[int(21)] = 0.0f;
    _S3219[int(14)] = _S3218;
    float _S3220 = _S3193 + _S3219[int(0)];
    float _S3221 = _S3194 + _S3219[int(1)];
    float _S3222 = _S3195 + _S3219[int(2)];
    float _S3223 = _S3196 + _S3219[int(3)];
    float _S3224 = _S3197 + _S3219[int(4)];
    float _S3225 = _S3198 + _S3219[int(5)];
    float _S3226 = _S3199 + _S3219[int(6)];
    float _S3227 = _S3200 + _S3219[int(7)];
    float _S3228 = _S3201 + _S3219[int(8)];
    float _S3229 = _S3202 + _S3219[int(9)];
    float _S3230 = _S3203 + _S3219[int(10)];
    float _S3231 = _S3204 + _S3219[int(11)];
    float _S3232 = _S3205 + _S3219[int(12)];
    float _S3233 = _S3206 + _S3219[int(13)];
    float _S3234 = _S3207 + _S3219[int(14)];
    float _S3235 = _S3208 + _S3219[int(15)];
    float _S3236 = _S3209 + _S3219[int(16)];
    float _S3237 = _S3210 + _S3219[int(17)];
    float _S3238 = _S3211 + _S3219[int(18)];
    float _S3239 = _S3212 + _S3219[int(19)];
    float _S3240 = _S3213 + _S3219[int(20)];
    float _S3241 = _S3214 + _S3219[int(21)];
    if(_S3090)
    {
        float _S3242 = 200.0f * _S3109;
        float _S3243 = _S3091 * _S3242 + 0.5f * (_S3089 * _S3242);
        _S3091 = 0.0f;
        _S3094 = _S3243;
    }
    else
    {
        _S3091 = _S3109;
        _S3094 = 0.0f;
    }
    DiffPair_float_0 _S3244;
    (&_S3244)->primal_0 = _S3089;
    (&_S3244)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3244, _S3091);
    float _S3245 = (_S3244.differential_0 + _S3094) / _S3073;
    FixedArray<float, 22>  _S3246;
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
    _S3246[int(13)] = _S3245;
    float _S3247 = _S3220 + _S3246[int(0)];
    float _S3248 = _S3221 + _S3246[int(1)];
    float _S3249 = _S3222 + _S3246[int(2)];
    float _S3250 = _S3223 + _S3246[int(3)];
    float _S3251 = _S3224 + _S3246[int(4)];
    float _S3252 = _S3225 + _S3246[int(5)];
    float _S3253 = _S3226 + _S3246[int(6)];
    float _S3254 = _S3227 + _S3246[int(7)];
    float _S3255 = _S3228 + _S3246[int(8)];
    float _S3256 = _S3229 + _S3246[int(9)];
    float _S3257 = _S3230 + _S3246[int(10)];
    float _S3258 = _S3231 + _S3246[int(11)];
    float _S3259 = _S3232 + _S3246[int(12)];
    float _S3260 = _S3233 + _S3246[int(13)];
    float _S3261 = _S3234 + _S3246[int(14)];
    float _S3262 = _S3235 + _S3246[int(15)];
    float _S3263 = _S3236 + _S3246[int(16)];
    float _S3264 = _S3237 + _S3246[int(17)];
    float _S3265 = _S3238 + _S3246[int(18)];
    float _S3266 = _S3239 + _S3246[int(19)];
    float _S3267 = _S3240 + _S3246[int(20)];
    float _S3268 = _S3241 + _S3246[int(21)];
    if(_S3087)
    {
        float _S3269 = 200.0f * _S3109;
        float _S3270 = _S3088 * _S3269 + 0.5f * (_S3086 * _S3269);
        _S3088 = 0.0f;
        _S3091 = _S3270;
    }
    else
    {
        _S3088 = _S3109;
        _S3091 = 0.0f;
    }
    DiffPair_float_0 _S3271;
    (&_S3271)->primal_0 = _S3086;
    (&_S3271)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3271, _S3088);
    float _S3272 = (_S3271.differential_0 + _S3091) / _S3073;
    FixedArray<float, 22>  _S3273;
    _S3273[int(0)] = 0.0f;
    _S3273[int(1)] = 0.0f;
    _S3273[int(2)] = 0.0f;
    _S3273[int(3)] = 0.0f;
    _S3273[int(4)] = 0.0f;
    _S3273[int(5)] = 0.0f;
    _S3273[int(6)] = 0.0f;
    _S3273[int(7)] = 0.0f;
    _S3273[int(8)] = 0.0f;
    _S3273[int(9)] = 0.0f;
    _S3273[int(10)] = 0.0f;
    _S3273[int(11)] = 0.0f;
    _S3273[int(12)] = 0.0f;
    _S3273[int(13)] = 0.0f;
    _S3273[int(14)] = 0.0f;
    _S3273[int(15)] = 0.0f;
    _S3273[int(16)] = 0.0f;
    _S3273[int(17)] = 0.0f;
    _S3273[int(18)] = 0.0f;
    _S3273[int(19)] = 0.0f;
    _S3273[int(20)] = 0.0f;
    _S3273[int(21)] = 0.0f;
    _S3273[int(12)] = _S3272;
    float _S3274 = _S3247 + _S3273[int(0)];
    float _S3275 = _S3248 + _S3273[int(1)];
    float _S3276 = _S3249 + _S3273[int(2)];
    float _S3277 = _S3250 + _S3273[int(3)];
    float _S3278 = _S3251 + _S3273[int(4)];
    float _S3279 = _S3252 + _S3273[int(5)];
    float _S3280 = _S3253 + _S3273[int(6)];
    float _S3281 = _S3254 + _S3273[int(7)];
    float _S3282 = _S3255 + _S3273[int(8)];
    float _S3283 = _S3256 + _S3273[int(9)];
    float _S3284 = _S3257 + _S3273[int(10)];
    float _S3285 = _S3258 + _S3273[int(11)];
    float _S3286 = _S3259 + _S3273[int(12)];
    float _S3287 = _S3260 + _S3273[int(13)];
    float _S3288 = _S3261 + _S3273[int(14)];
    float _S3289 = _S3262 + _S3273[int(15)];
    float _S3290 = _S3263 + _S3273[int(16)];
    float _S3291 = _S3264 + _S3273[int(17)];
    float _S3292 = _S3265 + _S3273[int(18)];
    float _S3293 = _S3266 + _S3273[int(19)];
    float _S3294 = _S3267 + _S3273[int(20)];
    float _S3295 = _S3268 + _S3273[int(21)];
    if(_S3084)
    {
        float _S3296 = 200.0f * _S3109;
        float _S3297 = _S3085 * _S3296 + 0.5f * (_S3083 * _S3296);
        _S3085 = 0.0f;
        _S3088 = _S3297;
    }
    else
    {
        _S3085 = _S3109;
        _S3088 = 0.0f;
    }
    DiffPair_float_0 _S3298;
    (&_S3298)->primal_0 = _S3083;
    (&_S3298)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3298, _S3085);
    float _S3299 = (_S3298.differential_0 + _S3088) / _S3073;
    FixedArray<float, 22>  _S3300;
    _S3300[int(0)] = 0.0f;
    _S3300[int(1)] = 0.0f;
    _S3300[int(2)] = 0.0f;
    _S3300[int(3)] = 0.0f;
    _S3300[int(4)] = 0.0f;
    _S3300[int(5)] = 0.0f;
    _S3300[int(6)] = 0.0f;
    _S3300[int(7)] = 0.0f;
    _S3300[int(8)] = 0.0f;
    _S3300[int(9)] = 0.0f;
    _S3300[int(10)] = 0.0f;
    _S3300[int(11)] = 0.0f;
    _S3300[int(12)] = 0.0f;
    _S3300[int(13)] = 0.0f;
    _S3300[int(14)] = 0.0f;
    _S3300[int(15)] = 0.0f;
    _S3300[int(16)] = 0.0f;
    _S3300[int(17)] = 0.0f;
    _S3300[int(18)] = 0.0f;
    _S3300[int(19)] = 0.0f;
    _S3300[int(20)] = 0.0f;
    _S3300[int(21)] = 0.0f;
    _S3300[int(11)] = _S3299;
    float _S3301 = _S3274 + _S3300[int(0)];
    float _S3302 = _S3275 + _S3300[int(1)];
    float _S3303 = _S3276 + _S3300[int(2)];
    float _S3304 = _S3277 + _S3300[int(3)];
    float _S3305 = _S3278 + _S3300[int(4)];
    float _S3306 = _S3279 + _S3300[int(5)];
    float _S3307 = _S3280 + _S3300[int(6)];
    float _S3308 = _S3281 + _S3300[int(7)];
    float _S3309 = _S3282 + _S3300[int(8)];
    float _S3310 = _S3283 + _S3300[int(9)];
    float _S3311 = _S3284 + _S3300[int(10)];
    float _S3312 = _S3285 + _S3300[int(11)];
    float _S3313 = _S3286 + _S3300[int(12)];
    float _S3314 = _S3287 + _S3300[int(13)];
    float _S3315 = _S3288 + _S3300[int(14)];
    float _S3316 = _S3289 + _S3300[int(15)];
    float _S3317 = _S3290 + _S3300[int(16)];
    float _S3318 = _S3291 + _S3300[int(17)];
    float _S3319 = _S3292 + _S3300[int(18)];
    float _S3320 = _S3293 + _S3300[int(19)];
    float _S3321 = _S3294 + _S3300[int(20)];
    float _S3322 = _S3295 + _S3300[int(21)];
    if(_S3081)
    {
        float _S3323 = 200.0f * _S3109;
        float _S3324 = _S3082 * _S3323 + 0.5f * (_S3080 * _S3323);
        _S3082 = 0.0f;
        _S3085 = _S3324;
    }
    else
    {
        _S3082 = _S3109;
        _S3085 = 0.0f;
    }
    DiffPair_float_0 _S3325;
    (&_S3325)->primal_0 = _S3080;
    (&_S3325)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3325, _S3082);
    float _S3326 = (_S3325.differential_0 + _S3085) / _S3073;
    float _S3327 = _S3104 / _S3079;
    float _S3328 = _S3105 / _S3078;
    float _S3329 = _S3106 / _S3077;
    FixedArray<float, 22>  _S3330;
    _S3330[int(0)] = 0.0f;
    _S3330[int(1)] = 0.0f;
    _S3330[int(2)] = 0.0f;
    _S3330[int(3)] = 0.0f;
    _S3330[int(4)] = 0.0f;
    _S3330[int(5)] = 0.0f;
    _S3330[int(6)] = 0.0f;
    _S3330[int(7)] = 0.0f;
    _S3330[int(8)] = 0.0f;
    _S3330[int(9)] = 0.0f;
    _S3330[int(10)] = 0.0f;
    _S3330[int(11)] = 0.0f;
    _S3330[int(12)] = 0.0f;
    _S3330[int(13)] = 0.0f;
    _S3330[int(14)] = 0.0f;
    _S3330[int(15)] = 0.0f;
    _S3330[int(16)] = 0.0f;
    _S3330[int(17)] = 0.0f;
    _S3330[int(18)] = 0.0f;
    _S3330[int(19)] = 0.0f;
    _S3330[int(20)] = 0.0f;
    _S3330[int(21)] = 0.0f;
    _S3330[int(10)] = _S3326;
    _S3330[int(9)] = _S3327;
    _S3330[int(8)] = _S3327;
    _S3330[int(7)] = _S3327;
    _S3330[int(6)] = _S3327;
    _S3330[int(5)] = _S3327;
    _S3330[int(4)] = _S3328;
    _S3330[int(3)] = _S3328;
    _S3330[int(2)] = _S3328;
    _S3330[int(1)] = _S3329;
    float _S3331 = _S3301 + _S3330[int(0)];
    float _S3332 = _S3302 + _S3330[int(1)];
    float _S3333 = _S3303 + _S3330[int(2)];
    float _S3334 = _S3304 + _S3330[int(3)];
    float _S3335 = _S3305 + _S3330[int(4)];
    float _S3336 = _S3306 + _S3330[int(5)];
    float _S3337 = _S3307 + _S3330[int(6)];
    float _S3338 = _S3308 + _S3330[int(7)];
    float _S3339 = _S3309 + _S3330[int(8)];
    float _S3340 = _S3310 + _S3330[int(9)];
    float _S3341 = _S3311 + _S3330[int(10)];
    float _S3342 = _S3312 + _S3330[int(11)];
    float _S3343 = _S3313 + _S3330[int(12)];
    float _S3344 = _S3314 + _S3330[int(13)];
    float _S3345 = _S3315 + _S3330[int(14)];
    float _S3346 = _S3316 + _S3330[int(15)];
    float _S3347 = _S3317 + _S3330[int(16)];
    float _S3348 = _S3318 + _S3330[int(17)];
    float _S3349 = _S3319 + _S3330[int(18)];
    float _S3350 = _S3320 + _S3330[int(19)];
    float _S3351 = _S3321 + _S3330[int(20)];
    float _S3352 = _S3322 + _S3330[int(21)];
    if(_S3075)
    {
        float _S3353 = 10.0f * _S3107;
        float _S3354 = _S3076 * _S3353 + 0.5f * (_S3074 * _S3353);
        _S3076 = 0.0f;
        _S3082 = _S3354;
    }
    else
    {
        _S3076 = _S3107;
        _S3082 = 0.0f;
    }
    DiffPair_float_0 _S3355;
    (&_S3355)->primal_0 = _S3074;
    (&_S3355)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3355, _S3076);
    float _S3356 = (_S3355.differential_0 + _S3082) / _S3073;
    FixedArray<float, 22>  _S3357;
    _S3357[int(0)] = 0.0f;
    _S3357[int(1)] = 0.0f;
    _S3357[int(2)] = 0.0f;
    _S3357[int(3)] = 0.0f;
    _S3357[int(4)] = 0.0f;
    _S3357[int(5)] = 0.0f;
    _S3357[int(6)] = 0.0f;
    _S3357[int(7)] = 0.0f;
    _S3357[int(8)] = 0.0f;
    _S3357[int(9)] = 0.0f;
    _S3357[int(10)] = 0.0f;
    _S3357[int(11)] = 0.0f;
    _S3357[int(12)] = 0.0f;
    _S3357[int(13)] = 0.0f;
    _S3357[int(14)] = 0.0f;
    _S3357[int(15)] = 0.0f;
    _S3357[int(16)] = 0.0f;
    _S3357[int(17)] = 0.0f;
    _S3357[int(18)] = 0.0f;
    _S3357[int(19)] = 0.0f;
    _S3357[int(20)] = 0.0f;
    _S3357[int(21)] = 0.0f;
    _S3357[int(0)] = _S3356;
    FixedArray<float, 22>  _S3358 = {
        _S3331 + _S3357[int(0)], _S3332 + _S3357[int(1)], _S3333 + _S3357[int(2)], _S3334 + _S3357[int(3)], _S3335 + _S3357[int(4)], _S3336 + _S3357[int(5)], _S3337 + _S3357[int(6)], _S3338 + _S3357[int(7)], _S3339 + _S3357[int(8)], _S3340 + _S3357[int(9)], _S3341 + _S3357[int(10)], _S3342 + _S3357[int(11)], _S3343 + _S3357[int(12)], _S3344 + _S3357[int(13)], _S3345 + _S3357[int(14)], _S3346 + _S3357[int(15)], _S3347 + _S3357[int(16)], _S3348 + _S3357[int(17)], _S3349 + _S3357[int(18)], _S3350 + _S3357[int(19)], _S3351 + _S3357[int(20)], _S3352 + _S3357[int(21)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S3358;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S3359, int _S3360, FixedArray<float, 6>  * _S3361, FixedArray<float, 6>  * _S3362)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S3359, _S3360, _S3361, _S3362);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_4, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_4, FixedArray<float, 22>  * _S3363)
{
    FixedArray<float, 22>  _S3364 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_4;
    (&dp_raw_losses_1)->differential_0 = _S3364;
    FixedArray<float, 6>  _S3365 = loss_weights_3;
    FixedArray<float, 6>  _S3366 = grad_out_4;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_1, num_cameras_3, &_S3365, &_S3366);
    *_S3363 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_2, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_13)
{
    FixedArray<float, 23>  _S3367 = dpraw_losses_2->primal_0;
    float _S3368 = float(num_cameras_4);
    float _S3369 = dpraw_losses_2->primal_0[int(0)] / _S3368;
    bool _S3370 = (s_primal_ctx_abs_0(_S3369)) < 0.10000000149011612f;
    float _S3371;
    if(_S3370)
    {
        _S3371 = 0.5f * _S3369;
    }
    else
    {
        _S3371 = 0.0f;
    }
    float _S3372 = 3.0f * _S3368;
    float _S3373 = 9.0f * _S3368;
    float _S3374 = 5.0f * _S3368;
    float _S3375 = _S3367[int(10)] / _S3368;
    bool _S3376 = (s_primal_ctx_abs_0(_S3375)) < 0.00499999988824129f;
    float _S3377;
    if(_S3376)
    {
        _S3377 = 0.5f * _S3375;
    }
    else
    {
        _S3377 = 0.0f;
    }
    float _S3378 = _S3367[int(11)] / _S3368;
    bool _S3379 = (s_primal_ctx_abs_0(_S3378)) < 0.00499999988824129f;
    float _S3380;
    if(_S3379)
    {
        _S3380 = 0.5f * _S3378;
    }
    else
    {
        _S3380 = 0.0f;
    }
    float _S3381 = _S3367[int(12)] / _S3368;
    bool _S3382 = (s_primal_ctx_abs_0(_S3381)) < 0.00499999988824129f;
    float _S3383;
    if(_S3382)
    {
        _S3383 = 0.5f * _S3381;
    }
    else
    {
        _S3383 = 0.0f;
    }
    float _S3384 = _S3367[int(13)] / _S3368;
    bool _S3385 = (s_primal_ctx_abs_0(_S3384)) < 0.00499999988824129f;
    float _S3386;
    if(_S3385)
    {
        _S3386 = 0.5f * _S3384;
    }
    else
    {
        _S3386 = 0.0f;
    }
    float _S3387 = _S3367[int(14)] / _S3368;
    bool _S3388 = (s_primal_ctx_abs_0(_S3387)) < 0.00499999988824129f;
    float _S3389;
    if(_S3388)
    {
        _S3389 = 0.5f * _S3387;
    }
    else
    {
        _S3389 = 0.0f;
    }
    float _S3390 = _S3367[int(15)] / _S3368;
    bool _S3391 = (s_primal_ctx_abs_0(_S3390)) < 0.00499999988824129f;
    float _S3392;
    if(_S3391)
    {
        _S3392 = 0.5f * _S3390;
    }
    else
    {
        _S3392 = 0.0f;
    }
    float _S3393 = _S3367[int(16)] / _S3368;
    bool _S3394 = (s_primal_ctx_abs_0(_S3393)) < 0.00499999988824129f;
    float _S3395;
    if(_S3394)
    {
        _S3395 = 0.5f * _S3393;
    }
    else
    {
        _S3395 = 0.0f;
    }
    float _S3396 = _S3367[int(17)] / _S3368;
    bool _S3397 = (s_primal_ctx_abs_0(_S3396)) < 0.00499999988824129f;
    float _S3398;
    if(_S3397)
    {
        _S3398 = 0.5f * _S3396;
    }
    else
    {
        _S3398 = 0.0f;
    }
    float _S3399 = (*loss_weights_4)[int(3)] * (*_s_dOut_13)[int(3)];
    float _S3400 = (*loss_weights_4)[int(2)] * (*_s_dOut_13)[int(2)];
    float _S3401 = (*loss_weights_4)[int(1)] * (*_s_dOut_13)[int(1)];
    float _S3402 = (*loss_weights_4)[int(0)] * (*_s_dOut_13)[int(0)];
    float _S3403 = (*loss_weights_4)[int(5)] * (*_s_dOut_13)[int(5)] / _S3374;
    float _S3404 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_13)[int(4)]);
    FixedArray<float, 23>  _S3405;
    _S3405[int(0)] = 0.0f;
    _S3405[int(1)] = 0.0f;
    _S3405[int(2)] = 0.0f;
    _S3405[int(3)] = 0.0f;
    _S3405[int(4)] = 0.0f;
    _S3405[int(5)] = 0.0f;
    _S3405[int(6)] = 0.0f;
    _S3405[int(7)] = 0.0f;
    _S3405[int(8)] = 0.0f;
    _S3405[int(9)] = 0.0f;
    _S3405[int(10)] = 0.0f;
    _S3405[int(11)] = 0.0f;
    _S3405[int(12)] = 0.0f;
    _S3405[int(13)] = 0.0f;
    _S3405[int(14)] = 0.0f;
    _S3405[int(15)] = 0.0f;
    _S3405[int(16)] = 0.0f;
    _S3405[int(17)] = 0.0f;
    _S3405[int(18)] = 0.0f;
    _S3405[int(19)] = 0.0f;
    _S3405[int(20)] = 0.0f;
    _S3405[int(21)] = 0.0f;
    _S3405[int(22)] = 0.0f;
    _S3405[int(22)] = _S3403;
    _S3405[int(21)] = _S3403;
    _S3405[int(20)] = _S3403;
    _S3405[int(19)] = _S3403;
    _S3405[int(18)] = _S3403;
    float _S3406 = _S3405[int(0)];
    float _S3407 = _S3405[int(1)];
    float _S3408 = _S3405[int(2)];
    float _S3409 = _S3405[int(3)];
    float _S3410 = _S3405[int(4)];
    float _S3411 = _S3405[int(5)];
    float _S3412 = _S3405[int(6)];
    float _S3413 = _S3405[int(7)];
    float _S3414 = _S3405[int(8)];
    float _S3415 = _S3405[int(9)];
    float _S3416 = _S3405[int(10)];
    float _S3417 = _S3405[int(11)];
    float _S3418 = _S3405[int(12)];
    float _S3419 = _S3405[int(13)];
    float _S3420 = _S3405[int(14)];
    float _S3421 = _S3405[int(15)];
    float _S3422 = _S3405[int(16)];
    float _S3423 = _S3405[int(17)];
    float _S3424 = _S3405[int(18)];
    float _S3425 = _S3405[int(19)];
    float _S3426 = _S3405[int(20)];
    float _S3427 = _S3405[int(21)];
    float _S3428 = _S3405[int(22)];
    float _S3429;
    if(_S3397)
    {
        float _S3430 = 200.0f * _S3404;
        float _S3431 = _S3398 * _S3430 + 0.5f * (_S3396 * _S3430);
        _S3398 = 0.0f;
        _S3429 = _S3431;
    }
    else
    {
        _S3398 = _S3404;
        _S3429 = 0.0f;
    }
    DiffPair_float_0 _S3432;
    (&_S3432)->primal_0 = _S3396;
    (&_S3432)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3432, _S3398);
    float _S3433 = (_S3432.differential_0 + _S3429) / _S3368;
    FixedArray<float, 23>  _S3434;
    _S3434[int(0)] = 0.0f;
    _S3434[int(1)] = 0.0f;
    _S3434[int(2)] = 0.0f;
    _S3434[int(3)] = 0.0f;
    _S3434[int(4)] = 0.0f;
    _S3434[int(5)] = 0.0f;
    _S3434[int(6)] = 0.0f;
    _S3434[int(7)] = 0.0f;
    _S3434[int(8)] = 0.0f;
    _S3434[int(9)] = 0.0f;
    _S3434[int(10)] = 0.0f;
    _S3434[int(11)] = 0.0f;
    _S3434[int(12)] = 0.0f;
    _S3434[int(13)] = 0.0f;
    _S3434[int(14)] = 0.0f;
    _S3434[int(15)] = 0.0f;
    _S3434[int(16)] = 0.0f;
    _S3434[int(17)] = 0.0f;
    _S3434[int(18)] = 0.0f;
    _S3434[int(19)] = 0.0f;
    _S3434[int(20)] = 0.0f;
    _S3434[int(21)] = 0.0f;
    _S3434[int(22)] = 0.0f;
    _S3434[int(17)] = _S3433;
    float _S3435 = _S3406 + _S3434[int(0)];
    float _S3436 = _S3407 + _S3434[int(1)];
    float _S3437 = _S3408 + _S3434[int(2)];
    float _S3438 = _S3409 + _S3434[int(3)];
    float _S3439 = _S3410 + _S3434[int(4)];
    float _S3440 = _S3411 + _S3434[int(5)];
    float _S3441 = _S3412 + _S3434[int(6)];
    float _S3442 = _S3413 + _S3434[int(7)];
    float _S3443 = _S3414 + _S3434[int(8)];
    float _S3444 = _S3415 + _S3434[int(9)];
    float _S3445 = _S3416 + _S3434[int(10)];
    float _S3446 = _S3417 + _S3434[int(11)];
    float _S3447 = _S3418 + _S3434[int(12)];
    float _S3448 = _S3419 + _S3434[int(13)];
    float _S3449 = _S3420 + _S3434[int(14)];
    float _S3450 = _S3421 + _S3434[int(15)];
    float _S3451 = _S3422 + _S3434[int(16)];
    float _S3452 = _S3423 + _S3434[int(17)];
    float _S3453 = _S3424 + _S3434[int(18)];
    float _S3454 = _S3425 + _S3434[int(19)];
    float _S3455 = _S3426 + _S3434[int(20)];
    float _S3456 = _S3427 + _S3434[int(21)];
    float _S3457 = _S3428 + _S3434[int(22)];
    if(_S3394)
    {
        float _S3458 = 200.0f * _S3404;
        float _S3459 = _S3395 * _S3458 + 0.5f * (_S3393 * _S3458);
        _S3395 = 0.0f;
        _S3398 = _S3459;
    }
    else
    {
        _S3395 = _S3404;
        _S3398 = 0.0f;
    }
    DiffPair_float_0 _S3460;
    (&_S3460)->primal_0 = _S3393;
    (&_S3460)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3460, _S3395);
    float _S3461 = (_S3460.differential_0 + _S3398) / _S3368;
    FixedArray<float, 23>  _S3462;
    _S3462[int(0)] = 0.0f;
    _S3462[int(1)] = 0.0f;
    _S3462[int(2)] = 0.0f;
    _S3462[int(3)] = 0.0f;
    _S3462[int(4)] = 0.0f;
    _S3462[int(5)] = 0.0f;
    _S3462[int(6)] = 0.0f;
    _S3462[int(7)] = 0.0f;
    _S3462[int(8)] = 0.0f;
    _S3462[int(9)] = 0.0f;
    _S3462[int(10)] = 0.0f;
    _S3462[int(11)] = 0.0f;
    _S3462[int(12)] = 0.0f;
    _S3462[int(13)] = 0.0f;
    _S3462[int(14)] = 0.0f;
    _S3462[int(15)] = 0.0f;
    _S3462[int(16)] = 0.0f;
    _S3462[int(17)] = 0.0f;
    _S3462[int(18)] = 0.0f;
    _S3462[int(19)] = 0.0f;
    _S3462[int(20)] = 0.0f;
    _S3462[int(21)] = 0.0f;
    _S3462[int(22)] = 0.0f;
    _S3462[int(16)] = _S3461;
    float _S3463 = _S3435 + _S3462[int(0)];
    float _S3464 = _S3436 + _S3462[int(1)];
    float _S3465 = _S3437 + _S3462[int(2)];
    float _S3466 = _S3438 + _S3462[int(3)];
    float _S3467 = _S3439 + _S3462[int(4)];
    float _S3468 = _S3440 + _S3462[int(5)];
    float _S3469 = _S3441 + _S3462[int(6)];
    float _S3470 = _S3442 + _S3462[int(7)];
    float _S3471 = _S3443 + _S3462[int(8)];
    float _S3472 = _S3444 + _S3462[int(9)];
    float _S3473 = _S3445 + _S3462[int(10)];
    float _S3474 = _S3446 + _S3462[int(11)];
    float _S3475 = _S3447 + _S3462[int(12)];
    float _S3476 = _S3448 + _S3462[int(13)];
    float _S3477 = _S3449 + _S3462[int(14)];
    float _S3478 = _S3450 + _S3462[int(15)];
    float _S3479 = _S3451 + _S3462[int(16)];
    float _S3480 = _S3452 + _S3462[int(17)];
    float _S3481 = _S3453 + _S3462[int(18)];
    float _S3482 = _S3454 + _S3462[int(19)];
    float _S3483 = _S3455 + _S3462[int(20)];
    float _S3484 = _S3456 + _S3462[int(21)];
    float _S3485 = _S3457 + _S3462[int(22)];
    if(_S3391)
    {
        float _S3486 = 200.0f * _S3404;
        float _S3487 = _S3392 * _S3486 + 0.5f * (_S3390 * _S3486);
        _S3392 = 0.0f;
        _S3395 = _S3487;
    }
    else
    {
        _S3392 = _S3404;
        _S3395 = 0.0f;
    }
    DiffPair_float_0 _S3488;
    (&_S3488)->primal_0 = _S3390;
    (&_S3488)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3488, _S3392);
    float _S3489 = (_S3488.differential_0 + _S3395) / _S3368;
    FixedArray<float, 23>  _S3490;
    _S3490[int(0)] = 0.0f;
    _S3490[int(1)] = 0.0f;
    _S3490[int(2)] = 0.0f;
    _S3490[int(3)] = 0.0f;
    _S3490[int(4)] = 0.0f;
    _S3490[int(5)] = 0.0f;
    _S3490[int(6)] = 0.0f;
    _S3490[int(7)] = 0.0f;
    _S3490[int(8)] = 0.0f;
    _S3490[int(9)] = 0.0f;
    _S3490[int(10)] = 0.0f;
    _S3490[int(11)] = 0.0f;
    _S3490[int(12)] = 0.0f;
    _S3490[int(13)] = 0.0f;
    _S3490[int(14)] = 0.0f;
    _S3490[int(15)] = 0.0f;
    _S3490[int(16)] = 0.0f;
    _S3490[int(17)] = 0.0f;
    _S3490[int(18)] = 0.0f;
    _S3490[int(19)] = 0.0f;
    _S3490[int(20)] = 0.0f;
    _S3490[int(21)] = 0.0f;
    _S3490[int(22)] = 0.0f;
    _S3490[int(15)] = _S3489;
    float _S3491 = _S3463 + _S3490[int(0)];
    float _S3492 = _S3464 + _S3490[int(1)];
    float _S3493 = _S3465 + _S3490[int(2)];
    float _S3494 = _S3466 + _S3490[int(3)];
    float _S3495 = _S3467 + _S3490[int(4)];
    float _S3496 = _S3468 + _S3490[int(5)];
    float _S3497 = _S3469 + _S3490[int(6)];
    float _S3498 = _S3470 + _S3490[int(7)];
    float _S3499 = _S3471 + _S3490[int(8)];
    float _S3500 = _S3472 + _S3490[int(9)];
    float _S3501 = _S3473 + _S3490[int(10)];
    float _S3502 = _S3474 + _S3490[int(11)];
    float _S3503 = _S3475 + _S3490[int(12)];
    float _S3504 = _S3476 + _S3490[int(13)];
    float _S3505 = _S3477 + _S3490[int(14)];
    float _S3506 = _S3478 + _S3490[int(15)];
    float _S3507 = _S3479 + _S3490[int(16)];
    float _S3508 = _S3480 + _S3490[int(17)];
    float _S3509 = _S3481 + _S3490[int(18)];
    float _S3510 = _S3482 + _S3490[int(19)];
    float _S3511 = _S3483 + _S3490[int(20)];
    float _S3512 = _S3484 + _S3490[int(21)];
    float _S3513 = _S3485 + _S3490[int(22)];
    if(_S3388)
    {
        float _S3514 = 200.0f * _S3404;
        float _S3515 = _S3389 * _S3514 + 0.5f * (_S3387 * _S3514);
        _S3389 = 0.0f;
        _S3392 = _S3515;
    }
    else
    {
        _S3389 = _S3404;
        _S3392 = 0.0f;
    }
    DiffPair_float_0 _S3516;
    (&_S3516)->primal_0 = _S3387;
    (&_S3516)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3516, _S3389);
    float _S3517 = (_S3516.differential_0 + _S3392) / _S3368;
    FixedArray<float, 23>  _S3518;
    _S3518[int(0)] = 0.0f;
    _S3518[int(1)] = 0.0f;
    _S3518[int(2)] = 0.0f;
    _S3518[int(3)] = 0.0f;
    _S3518[int(4)] = 0.0f;
    _S3518[int(5)] = 0.0f;
    _S3518[int(6)] = 0.0f;
    _S3518[int(7)] = 0.0f;
    _S3518[int(8)] = 0.0f;
    _S3518[int(9)] = 0.0f;
    _S3518[int(10)] = 0.0f;
    _S3518[int(11)] = 0.0f;
    _S3518[int(12)] = 0.0f;
    _S3518[int(13)] = 0.0f;
    _S3518[int(14)] = 0.0f;
    _S3518[int(15)] = 0.0f;
    _S3518[int(16)] = 0.0f;
    _S3518[int(17)] = 0.0f;
    _S3518[int(18)] = 0.0f;
    _S3518[int(19)] = 0.0f;
    _S3518[int(20)] = 0.0f;
    _S3518[int(21)] = 0.0f;
    _S3518[int(22)] = 0.0f;
    _S3518[int(14)] = _S3517;
    float _S3519 = _S3491 + _S3518[int(0)];
    float _S3520 = _S3492 + _S3518[int(1)];
    float _S3521 = _S3493 + _S3518[int(2)];
    float _S3522 = _S3494 + _S3518[int(3)];
    float _S3523 = _S3495 + _S3518[int(4)];
    float _S3524 = _S3496 + _S3518[int(5)];
    float _S3525 = _S3497 + _S3518[int(6)];
    float _S3526 = _S3498 + _S3518[int(7)];
    float _S3527 = _S3499 + _S3518[int(8)];
    float _S3528 = _S3500 + _S3518[int(9)];
    float _S3529 = _S3501 + _S3518[int(10)];
    float _S3530 = _S3502 + _S3518[int(11)];
    float _S3531 = _S3503 + _S3518[int(12)];
    float _S3532 = _S3504 + _S3518[int(13)];
    float _S3533 = _S3505 + _S3518[int(14)];
    float _S3534 = _S3506 + _S3518[int(15)];
    float _S3535 = _S3507 + _S3518[int(16)];
    float _S3536 = _S3508 + _S3518[int(17)];
    float _S3537 = _S3509 + _S3518[int(18)];
    float _S3538 = _S3510 + _S3518[int(19)];
    float _S3539 = _S3511 + _S3518[int(20)];
    float _S3540 = _S3512 + _S3518[int(21)];
    float _S3541 = _S3513 + _S3518[int(22)];
    if(_S3385)
    {
        float _S3542 = 200.0f * _S3404;
        float _S3543 = _S3386 * _S3542 + 0.5f * (_S3384 * _S3542);
        _S3386 = 0.0f;
        _S3389 = _S3543;
    }
    else
    {
        _S3386 = _S3404;
        _S3389 = 0.0f;
    }
    DiffPair_float_0 _S3544;
    (&_S3544)->primal_0 = _S3384;
    (&_S3544)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3544, _S3386);
    float _S3545 = (_S3544.differential_0 + _S3389) / _S3368;
    FixedArray<float, 23>  _S3546;
    _S3546[int(0)] = 0.0f;
    _S3546[int(1)] = 0.0f;
    _S3546[int(2)] = 0.0f;
    _S3546[int(3)] = 0.0f;
    _S3546[int(4)] = 0.0f;
    _S3546[int(5)] = 0.0f;
    _S3546[int(6)] = 0.0f;
    _S3546[int(7)] = 0.0f;
    _S3546[int(8)] = 0.0f;
    _S3546[int(9)] = 0.0f;
    _S3546[int(10)] = 0.0f;
    _S3546[int(11)] = 0.0f;
    _S3546[int(12)] = 0.0f;
    _S3546[int(13)] = 0.0f;
    _S3546[int(14)] = 0.0f;
    _S3546[int(15)] = 0.0f;
    _S3546[int(16)] = 0.0f;
    _S3546[int(17)] = 0.0f;
    _S3546[int(18)] = 0.0f;
    _S3546[int(19)] = 0.0f;
    _S3546[int(20)] = 0.0f;
    _S3546[int(21)] = 0.0f;
    _S3546[int(22)] = 0.0f;
    _S3546[int(13)] = _S3545;
    float _S3547 = _S3519 + _S3546[int(0)];
    float _S3548 = _S3520 + _S3546[int(1)];
    float _S3549 = _S3521 + _S3546[int(2)];
    float _S3550 = _S3522 + _S3546[int(3)];
    float _S3551 = _S3523 + _S3546[int(4)];
    float _S3552 = _S3524 + _S3546[int(5)];
    float _S3553 = _S3525 + _S3546[int(6)];
    float _S3554 = _S3526 + _S3546[int(7)];
    float _S3555 = _S3527 + _S3546[int(8)];
    float _S3556 = _S3528 + _S3546[int(9)];
    float _S3557 = _S3529 + _S3546[int(10)];
    float _S3558 = _S3530 + _S3546[int(11)];
    float _S3559 = _S3531 + _S3546[int(12)];
    float _S3560 = _S3532 + _S3546[int(13)];
    float _S3561 = _S3533 + _S3546[int(14)];
    float _S3562 = _S3534 + _S3546[int(15)];
    float _S3563 = _S3535 + _S3546[int(16)];
    float _S3564 = _S3536 + _S3546[int(17)];
    float _S3565 = _S3537 + _S3546[int(18)];
    float _S3566 = _S3538 + _S3546[int(19)];
    float _S3567 = _S3539 + _S3546[int(20)];
    float _S3568 = _S3540 + _S3546[int(21)];
    float _S3569 = _S3541 + _S3546[int(22)];
    if(_S3382)
    {
        float _S3570 = 200.0f * _S3404;
        float _S3571 = _S3383 * _S3570 + 0.5f * (_S3381 * _S3570);
        _S3383 = 0.0f;
        _S3386 = _S3571;
    }
    else
    {
        _S3383 = _S3404;
        _S3386 = 0.0f;
    }
    DiffPair_float_0 _S3572;
    (&_S3572)->primal_0 = _S3381;
    (&_S3572)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3572, _S3383);
    float _S3573 = (_S3572.differential_0 + _S3386) / _S3368;
    FixedArray<float, 23>  _S3574;
    _S3574[int(0)] = 0.0f;
    _S3574[int(1)] = 0.0f;
    _S3574[int(2)] = 0.0f;
    _S3574[int(3)] = 0.0f;
    _S3574[int(4)] = 0.0f;
    _S3574[int(5)] = 0.0f;
    _S3574[int(6)] = 0.0f;
    _S3574[int(7)] = 0.0f;
    _S3574[int(8)] = 0.0f;
    _S3574[int(9)] = 0.0f;
    _S3574[int(10)] = 0.0f;
    _S3574[int(11)] = 0.0f;
    _S3574[int(12)] = 0.0f;
    _S3574[int(13)] = 0.0f;
    _S3574[int(14)] = 0.0f;
    _S3574[int(15)] = 0.0f;
    _S3574[int(16)] = 0.0f;
    _S3574[int(17)] = 0.0f;
    _S3574[int(18)] = 0.0f;
    _S3574[int(19)] = 0.0f;
    _S3574[int(20)] = 0.0f;
    _S3574[int(21)] = 0.0f;
    _S3574[int(22)] = 0.0f;
    _S3574[int(12)] = _S3573;
    float _S3575 = _S3547 + _S3574[int(0)];
    float _S3576 = _S3548 + _S3574[int(1)];
    float _S3577 = _S3549 + _S3574[int(2)];
    float _S3578 = _S3550 + _S3574[int(3)];
    float _S3579 = _S3551 + _S3574[int(4)];
    float _S3580 = _S3552 + _S3574[int(5)];
    float _S3581 = _S3553 + _S3574[int(6)];
    float _S3582 = _S3554 + _S3574[int(7)];
    float _S3583 = _S3555 + _S3574[int(8)];
    float _S3584 = _S3556 + _S3574[int(9)];
    float _S3585 = _S3557 + _S3574[int(10)];
    float _S3586 = _S3558 + _S3574[int(11)];
    float _S3587 = _S3559 + _S3574[int(12)];
    float _S3588 = _S3560 + _S3574[int(13)];
    float _S3589 = _S3561 + _S3574[int(14)];
    float _S3590 = _S3562 + _S3574[int(15)];
    float _S3591 = _S3563 + _S3574[int(16)];
    float _S3592 = _S3564 + _S3574[int(17)];
    float _S3593 = _S3565 + _S3574[int(18)];
    float _S3594 = _S3566 + _S3574[int(19)];
    float _S3595 = _S3567 + _S3574[int(20)];
    float _S3596 = _S3568 + _S3574[int(21)];
    float _S3597 = _S3569 + _S3574[int(22)];
    if(_S3379)
    {
        float _S3598 = 200.0f * _S3404;
        float _S3599 = _S3380 * _S3598 + 0.5f * (_S3378 * _S3598);
        _S3380 = 0.0f;
        _S3383 = _S3599;
    }
    else
    {
        _S3380 = _S3404;
        _S3383 = 0.0f;
    }
    DiffPair_float_0 _S3600;
    (&_S3600)->primal_0 = _S3378;
    (&_S3600)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3600, _S3380);
    float _S3601 = (_S3600.differential_0 + _S3383) / _S3368;
    FixedArray<float, 23>  _S3602;
    _S3602[int(0)] = 0.0f;
    _S3602[int(1)] = 0.0f;
    _S3602[int(2)] = 0.0f;
    _S3602[int(3)] = 0.0f;
    _S3602[int(4)] = 0.0f;
    _S3602[int(5)] = 0.0f;
    _S3602[int(6)] = 0.0f;
    _S3602[int(7)] = 0.0f;
    _S3602[int(8)] = 0.0f;
    _S3602[int(9)] = 0.0f;
    _S3602[int(10)] = 0.0f;
    _S3602[int(11)] = 0.0f;
    _S3602[int(12)] = 0.0f;
    _S3602[int(13)] = 0.0f;
    _S3602[int(14)] = 0.0f;
    _S3602[int(15)] = 0.0f;
    _S3602[int(16)] = 0.0f;
    _S3602[int(17)] = 0.0f;
    _S3602[int(18)] = 0.0f;
    _S3602[int(19)] = 0.0f;
    _S3602[int(20)] = 0.0f;
    _S3602[int(21)] = 0.0f;
    _S3602[int(22)] = 0.0f;
    _S3602[int(11)] = _S3601;
    float _S3603 = _S3575 + _S3602[int(0)];
    float _S3604 = _S3576 + _S3602[int(1)];
    float _S3605 = _S3577 + _S3602[int(2)];
    float _S3606 = _S3578 + _S3602[int(3)];
    float _S3607 = _S3579 + _S3602[int(4)];
    float _S3608 = _S3580 + _S3602[int(5)];
    float _S3609 = _S3581 + _S3602[int(6)];
    float _S3610 = _S3582 + _S3602[int(7)];
    float _S3611 = _S3583 + _S3602[int(8)];
    float _S3612 = _S3584 + _S3602[int(9)];
    float _S3613 = _S3585 + _S3602[int(10)];
    float _S3614 = _S3586 + _S3602[int(11)];
    float _S3615 = _S3587 + _S3602[int(12)];
    float _S3616 = _S3588 + _S3602[int(13)];
    float _S3617 = _S3589 + _S3602[int(14)];
    float _S3618 = _S3590 + _S3602[int(15)];
    float _S3619 = _S3591 + _S3602[int(16)];
    float _S3620 = _S3592 + _S3602[int(17)];
    float _S3621 = _S3593 + _S3602[int(18)];
    float _S3622 = _S3594 + _S3602[int(19)];
    float _S3623 = _S3595 + _S3602[int(20)];
    float _S3624 = _S3596 + _S3602[int(21)];
    float _S3625 = _S3597 + _S3602[int(22)];
    if(_S3376)
    {
        float _S3626 = 200.0f * _S3404;
        float _S3627 = _S3377 * _S3626 + 0.5f * (_S3375 * _S3626);
        _S3377 = 0.0f;
        _S3380 = _S3627;
    }
    else
    {
        _S3377 = _S3404;
        _S3380 = 0.0f;
    }
    DiffPair_float_0 _S3628;
    (&_S3628)->primal_0 = _S3375;
    (&_S3628)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3628, _S3377);
    float _S3629 = (_S3628.differential_0 + _S3380) / _S3368;
    float _S3630 = _S3399 / _S3374;
    float _S3631 = _S3400 / _S3373;
    float _S3632 = _S3401 / _S3372;
    FixedArray<float, 23>  _S3633;
    _S3633[int(0)] = 0.0f;
    _S3633[int(1)] = 0.0f;
    _S3633[int(2)] = 0.0f;
    _S3633[int(3)] = 0.0f;
    _S3633[int(4)] = 0.0f;
    _S3633[int(5)] = 0.0f;
    _S3633[int(6)] = 0.0f;
    _S3633[int(7)] = 0.0f;
    _S3633[int(8)] = 0.0f;
    _S3633[int(9)] = 0.0f;
    _S3633[int(10)] = 0.0f;
    _S3633[int(11)] = 0.0f;
    _S3633[int(12)] = 0.0f;
    _S3633[int(13)] = 0.0f;
    _S3633[int(14)] = 0.0f;
    _S3633[int(15)] = 0.0f;
    _S3633[int(16)] = 0.0f;
    _S3633[int(17)] = 0.0f;
    _S3633[int(18)] = 0.0f;
    _S3633[int(19)] = 0.0f;
    _S3633[int(20)] = 0.0f;
    _S3633[int(21)] = 0.0f;
    _S3633[int(22)] = 0.0f;
    _S3633[int(10)] = _S3629;
    _S3633[int(9)] = _S3630;
    _S3633[int(8)] = _S3630;
    _S3633[int(7)] = _S3630;
    _S3633[int(6)] = _S3630;
    _S3633[int(5)] = _S3630;
    _S3633[int(4)] = _S3631;
    _S3633[int(3)] = _S3631;
    _S3633[int(2)] = _S3631;
    _S3633[int(1)] = _S3632;
    float _S3634 = _S3603 + _S3633[int(0)];
    float _S3635 = _S3604 + _S3633[int(1)];
    float _S3636 = _S3605 + _S3633[int(2)];
    float _S3637 = _S3606 + _S3633[int(3)];
    float _S3638 = _S3607 + _S3633[int(4)];
    float _S3639 = _S3608 + _S3633[int(5)];
    float _S3640 = _S3609 + _S3633[int(6)];
    float _S3641 = _S3610 + _S3633[int(7)];
    float _S3642 = _S3611 + _S3633[int(8)];
    float _S3643 = _S3612 + _S3633[int(9)];
    float _S3644 = _S3613 + _S3633[int(10)];
    float _S3645 = _S3614 + _S3633[int(11)];
    float _S3646 = _S3615 + _S3633[int(12)];
    float _S3647 = _S3616 + _S3633[int(13)];
    float _S3648 = _S3617 + _S3633[int(14)];
    float _S3649 = _S3618 + _S3633[int(15)];
    float _S3650 = _S3619 + _S3633[int(16)];
    float _S3651 = _S3620 + _S3633[int(17)];
    float _S3652 = _S3621 + _S3633[int(18)];
    float _S3653 = _S3622 + _S3633[int(19)];
    float _S3654 = _S3623 + _S3633[int(20)];
    float _S3655 = _S3624 + _S3633[int(21)];
    float _S3656 = _S3625 + _S3633[int(22)];
    if(_S3370)
    {
        float _S3657 = 10.0f * _S3402;
        float _S3658 = _S3371 * _S3657 + 0.5f * (_S3369 * _S3657);
        _S3371 = 0.0f;
        _S3377 = _S3658;
    }
    else
    {
        _S3371 = _S3402;
        _S3377 = 0.0f;
    }
    DiffPair_float_0 _S3659;
    (&_S3659)->primal_0 = _S3369;
    (&_S3659)->differential_0 = 0.0f;
    s_bwd_prop_abs_1(&_S3659, _S3371);
    float _S3660 = (_S3659.differential_0 + _S3377) / _S3368;
    FixedArray<float, 23>  _S3661;
    _S3661[int(0)] = 0.0f;
    _S3661[int(1)] = 0.0f;
    _S3661[int(2)] = 0.0f;
    _S3661[int(3)] = 0.0f;
    _S3661[int(4)] = 0.0f;
    _S3661[int(5)] = 0.0f;
    _S3661[int(6)] = 0.0f;
    _S3661[int(7)] = 0.0f;
    _S3661[int(8)] = 0.0f;
    _S3661[int(9)] = 0.0f;
    _S3661[int(10)] = 0.0f;
    _S3661[int(11)] = 0.0f;
    _S3661[int(12)] = 0.0f;
    _S3661[int(13)] = 0.0f;
    _S3661[int(14)] = 0.0f;
    _S3661[int(15)] = 0.0f;
    _S3661[int(16)] = 0.0f;
    _S3661[int(17)] = 0.0f;
    _S3661[int(18)] = 0.0f;
    _S3661[int(19)] = 0.0f;
    _S3661[int(20)] = 0.0f;
    _S3661[int(21)] = 0.0f;
    _S3661[int(22)] = 0.0f;
    _S3661[int(0)] = _S3660;
    FixedArray<float, 23>  _S3662 = {
        _S3634 + _S3661[int(0)], _S3635 + _S3661[int(1)], _S3636 + _S3661[int(2)], _S3637 + _S3661[int(3)], _S3638 + _S3661[int(4)], _S3639 + _S3661[int(5)], _S3640 + _S3661[int(6)], _S3641 + _S3661[int(7)], _S3642 + _S3661[int(8)], _S3643 + _S3661[int(9)], _S3644 + _S3661[int(10)], _S3645 + _S3661[int(11)], _S3646 + _S3661[int(12)], _S3647 + _S3661[int(13)], _S3648 + _S3661[int(14)], _S3649 + _S3661[int(15)], _S3650 + _S3661[int(16)], _S3651 + _S3661[int(17)], _S3652 + _S3661[int(18)], _S3653 + _S3661[int(19)], _S3654 + _S3661[int(20)], _S3655 + _S3661[int(21)], _S3656 + _S3661[int(22)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S3662;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S3663, int _S3664, FixedArray<float, 6>  * _S3665, FixedArray<float, 6>  * _S3666)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S3663, _S3664, _S3665, _S3666);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_5, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_5, FixedArray<float, 23>  * _S3667)
{
    FixedArray<float, 23>  _S3668 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S3668;
    FixedArray<float, 6>  _S3669 = loss_weights_5;
    FixedArray<float, 6>  _S3670 = grad_out_5;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_2, num_cameras_5, &_S3669, &_S3670);
    *_S3667 = (&dp_raw_losses_2)->differential_0;
    return;
}

