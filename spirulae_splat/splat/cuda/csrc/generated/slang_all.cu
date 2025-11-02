#pragma once

#include "slang.cuh"

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S1 = (F32_exp(((*dpx_0).primal_0))) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S1;
    return;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S2 = *dpx_1;
    float _S3;
    if(((*dpx_1).primal_0) > ((*dpy_0).primal_0))
    {
        _S3 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_0).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S2.primal_0;
    dpx_1->differential_0 = _S3;
    DiffPair_float_0 _S4 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S2.primal_0))
    {
        _S3 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_1).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S4.primal_0;
    dpy_0->differential_0 = _S3;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S5 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_2).primal_0)))))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S5;
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
    float result_0 = 0.0f;
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_1 = result_0 + _slang_vector_get_element(x_0, i_0) * _slang_vector_get_element(y_0, i_0);
        i_0 = i_0 + int(1);
        result_0 = result_1;
    }
    return result_0;
}

inline __device__ float dot_1(float4  x_1, float4  y_1)
{
    int i_1 = int(0);
    float result_2 = 0.0f;
    for(;;)
    {
        if(i_1 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_3 = result_2 + _slang_vector_get_element(x_1, i_1) * _slang_vector_get_element(y_1, i_1);
        i_1 = i_1 + int(1);
        result_2 = result_3;
    }
    return result_2;
}

inline __device__ float dot_2(float2  x_2, float2  y_2)
{
    int i_2 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_2 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_5 = result_4 + _slang_vector_get_element(x_2, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_4 = result_5;
    }
    return result_4;
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
    float _S6 = 1.0f / (*dpx_4).primal_0 * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S6;
    return;
}

inline __device__ float3  exp_0(float3  x_6)
{
    float3  result_6;
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
        *_slang_vector_get_element_ptr(&result_6, i_3) = (F32_exp((_slang_vector_get_element(x_6, i_3))));
        i_3 = i_3 + int(1);
    }
    return result_6;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float3  dOut_5)
{
    float3  _S7 = exp_0((*dpx_5).primal_0) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S7;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ float2  exp_1(float2  x_7)
{
    float2  result_7;
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
        *_slang_vector_get_element_ptr(&result_7, i_4) = (F32_exp((_slang_vector_get_element(x_7, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_7;
}

inline __device__ void _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_6, float2  dOut_6)
{
    float2  _S8 = exp_1((*dpx_6).primal_0) * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S8;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpy_2, float dOut_7)
{
    DiffPair_float_0 _S9 = *dpx_7;
    float _S10;
    if(((*dpx_7).primal_0) < ((*dpy_2).primal_0))
    {
        _S10 = dOut_7;
    }
    else
    {
        if(((*dpx_7).primal_0) > ((*dpy_2).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_7;
        }
    }
    dpx_7->primal_0 = _S9.primal_0;
    dpx_7->differential_0 = _S10;
    DiffPair_float_0 _S11 = *dpy_2;
    if(((*dpy_2).primal_0) < (_S9.primal_0))
    {
        _S10 = dOut_7;
    }
    else
    {
        if(((*dpy_2).primal_0) > ((*dpx_7).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_7;
        }
    }
    dpy_2->primal_0 = _S11.primal_0;
    dpy_2->differential_0 = _S10;
    return;
}

inline __device__ void per_splat_losses(bool is_3dgs_0, float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S12)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    if(is_3dgs_0)
    {
        float3  _S13 = exp_0(scales_0);
        float _S14 = _S13.x;
        float _S15 = _S13.y;
        float _S16 = _S13.z;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S14 + _S15 + _S16) / 3.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S14), (_S15)))), (_S16))) / (F32_min(((F32_min((_S14), (_S15)))), (_S16)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float3  _S17 = exp_0(make_float3 (2.0f) * scales_0);
        float x_8 = _S17.x;
        float y_3 = _S17.y;
        float z_0 = _S17.z;
        float s_0 = x_8 + y_3 + z_0;
        float s1_0 = (F32_max(((F32_max((x_8), (y_3)))), (z_0))) / s_0;
        float s3_0 = (F32_min(((F32_min((x_8), (y_3)))), (z_0))) / s_0;
        float s2_0 = 1.0f - s1_0 - s3_0;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    }
    else
    {
        float2  _S18 = float2 {scales_0.x, scales_0.y};
        float2  _S19 = exp_1(_S18);
        float _S20 = _S19.x;
        float _S21 = _S19.y;
        losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S20 + _S21) / 2.0f;
        losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max((_S20), (_S21))) / (F32_min((_S20), (_S21)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
        float2  _S22 = exp_1(make_float2 (2.0f) * _S18);
        float x_9 = _S22.x;
        float y_4 = _S22.y;
        float s_1 = x_9 + y_4;
        float s1_1 = (F32_max((x_9), (y_4))) / s_1;
        float s2_1 = (F32_min((x_9), (y_4))) / s_1;
        losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_1 * (F32_log((s1_1))) - s2_1 * (F32_log((s2_1)))))) - 0.99998998641967773f)))), (0.0f)));
    }
    *_S12 = losses_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ float s_primal_ctx_exp_0(float _S23)
{
    return (F32_exp((_S23)));
}

inline __device__ float2  s_primal_ctx_exp_1(float2  _S24)
{
    return exp_1(_S24);
}

inline __device__ float s_primal_ctx_max_0(float _S25, float _S26)
{
    return (F32_max((_S25), (_S26)));
}

inline __device__ float s_primal_ctx_min_0(float _S27, float _S28)
{
    return (F32_min((_S27), (_S28)));
}

inline __device__ float s_primal_ctx_log_0(float _S29)
{
    return (F32_log((_S29)));
}

inline __device__ float3  s_primal_ctx_exp_2(float3  _S30)
{
    return exp_0(_S30);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S31, DiffPair_float_0 * _S32, float _S33)
{
    _d_max_0(_S31, _S32, _S33);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S34, float _S35)
{
    _d_log_0(_S34, _S35);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S36, float _S37)
{
    _d_exp_0(_S36, _S37);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S38, DiffPair_float_0 * _S39, float _S40)
{
    _d_min_0(_S38, _S39, _S40);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S41, float2  _S42)
{
    _d_exp_vector_1(_S41, _S42);
    return;
}

inline __device__ void s_bwd_prop_exp_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S43, float3  _S44)
{
    _d_exp_vector_0(_S43, _S44);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S45, float _S46)
{
    _d_sqrt_0(_S45, _S46);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_8, float _s_dOut_0)
{
    float _S47 = (*dpx_8).primal_0.x;
    float _S48 = (*dpx_8).primal_0.y;
    float _S49 = (*dpx_8).primal_0.z;
    float _S50 = (*dpx_8).primal_0.w;
    DiffPair_float_0 _S51;
    (&_S51)->primal_0 = _S47 * _S47 + _S48 * _S48 + _S49 * _S49 + _S50 * _S50;
    (&_S51)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S51, _s_dOut_0);
    float _S52 = (*dpx_8).primal_0.w * _S51.differential_0;
    float _S53 = _S52 + _S52;
    float _S54 = (*dpx_8).primal_0.z * _S51.differential_0;
    float _S55 = _S54 + _S54;
    float _S56 = (*dpx_8).primal_0.y * _S51.differential_0;
    float _S57 = _S56 + _S56;
    float _S58 = (*dpx_8).primal_0.x * _S51.differential_0;
    float _S59 = _S58 + _S58;
    float4  _S60 = make_float4 (0.0f);
    *&((&_S60)->w) = _S53;
    *&((&_S60)->z) = _S55;
    *&((&_S60)->y) = _S57;
    *&((&_S60)->x) = _S59;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S60;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S61, float _S62)
{
    s_bwd_prop_length_impl_0(_S61, _S62);
    return;
}

inline __device__ void s_bwd_prop_per_splat_losses_0(bool is_3dgs_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1, FixedArray<float, 5>  * _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S63 = *dpscales_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S64 = *dpquat_0;
    float2  _S65 = make_float2 (0.0f);
    float3  _S66 = make_float3 (0.0f);
    float _S67 = - (*dpopacity_0).primal_0;
    float _S68 = 1.0f + s_primal_ctx_exp_0(_S67);
    float _S69 = _S68 * _S68;
    float _S70 = length_0((*dpquat_0).primal_0);
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
    float _S114;
    float _S115;
    float _S116;
    float _S117;
    float _S118;
    float2  _S119;
    float2  _S120;
    float3  _S121;
    if(is_3dgs_1)
    {
        float3  _S122 = s_primal_ctx_exp_2(_S63.primal_0);
        float _S123 = _S122.x;
        float _S124 = _S122.y;
        float _S125 = _S122.z;
        float _S126 = s_primal_ctx_max_0(_S123, _S124);
        float _S127 = s_primal_ctx_max_0(_S126, _S125);
        float _S128 = s_primal_ctx_min_0(_S123, _S124);
        float _S129 = s_primal_ctx_min_0(_S128, _S125);
        float _S130 = _S127 / _S129;
        float _S131 = _S129 * _S129;
        float3  _S132 = make_float3 (2.0f) * _S63.primal_0;
        float3  _S133 = s_primal_ctx_exp_2(_S132);
        float x_10 = _S133.x;
        float y_5 = _S133.y;
        float z_1 = _S133.z;
        float s_2 = x_10 + y_5 + z_1;
        float _S134 = s_primal_ctx_max_0(x_10, y_5);
        float _S135 = s_primal_ctx_max_0(_S134, z_1);
        float s1_2 = _S135 / s_2;
        float _S136 = s_2 * s_2;
        float _S137 = s_primal_ctx_min_0(x_10, y_5);
        float _S138 = s_primal_ctx_min_0(_S137, z_1);
        float s3_1 = _S138 / s_2;
        float s2_2 = 1.0f - s1_2 - s3_1;
        float _S139 = - s1_2;
        float _S140 = s_primal_ctx_log_0(s1_2);
        float _S141 = s_primal_ctx_log_0(s2_2);
        float _S142 = s_primal_ctx_log_0(s3_1);
        float _S143 = _S139 * _S140 - s2_2 * _S141 - s3_1 * _S142;
        float _S144 = s_primal_ctx_exp_0(_S143) - 0.99998998641967773f;
        float _S145 = - s_primal_ctx_log_0(_S144);
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
        _S119 = _S65;
        _S85 = 0.0f;
        _S86 = 0.0f;
        _S87 = 0.0f;
        _S88 = 0.0f;
        _S89 = 0.0f;
        _S90 = 0.0f;
        _S120 = _S65;
        _S91 = _S145;
        _S92 = _S144;
        _S93 = _S143;
        _S94 = s3_1;
        _S95 = _S142;
        _S96 = s2_2;
        _S97 = _S141;
        _S98 = _S139;
        _S99 = _S140;
        _S100 = s1_2;
        _S101 = _S136;
        _S102 = _S138;
        _S103 = s_2;
        _S104 = _S137;
        _S105 = z_1;
        _S106 = x_10;
        _S107 = y_5;
        _S108 = _S135;
        _S109 = _S134;
        _S121 = _S132;
        _S110 = _S130;
        _S111 = _S131;
        _S112 = _S127;
        _S113 = _S129;
        _S114 = _S128;
        _S115 = _S125;
        _S116 = _S123;
        _S117 = _S124;
        _S118 = _S126;
    }
    else
    {
        float2  _S146 = float2 {_S63.primal_0.x, _S63.primal_0.y};
        float2  _S147 = s_primal_ctx_exp_1(_S146);
        float _S148 = _S147.x;
        float _S149 = _S147.y;
        float _S150 = s_primal_ctx_max_0(_S148, _S149);
        float _S151 = s_primal_ctx_min_0(_S148, _S149);
        float _S152 = _S150 / _S151;
        float _S153 = _S151 * _S151;
        float2  _S154 = make_float2 (2.0f) * _S146;
        float2  _S155 = s_primal_ctx_exp_1(_S154);
        float x_11 = _S155.x;
        float y_6 = _S155.y;
        float s_3 = x_11 + y_6;
        float _S156 = s_primal_ctx_max_0(x_11, y_6);
        float s1_3 = _S156 / s_3;
        float _S157 = s_3 * s_3;
        float _S158 = s_primal_ctx_min_0(x_11, y_6);
        float s2_3 = _S158 / s_3;
        float _S159 = - s1_3;
        float _S160 = s_primal_ctx_log_0(s1_3);
        float _S161 = s_primal_ctx_log_0(s2_3);
        float _S162 = _S159 * _S160 - s2_3 * _S161;
        float _S163 = s_primal_ctx_exp_0(_S162) - 0.99998998641967773f;
        _S71 = - s_primal_ctx_log_0(_S163);
        _S72 = _S163;
        _S73 = _S162;
        _S74 = s2_3;
        _S75 = _S161;
        _S76 = _S159;
        _S77 = _S160;
        _S78 = s1_3;
        _S79 = _S157;
        _S80 = _S158;
        _S81 = s_3;
        _S82 = x_11;
        _S83 = y_6;
        _S84 = _S156;
        _S119 = _S154;
        _S85 = _S152;
        _S86 = _S153;
        _S87 = _S150;
        _S88 = _S151;
        _S89 = _S148;
        _S90 = _S149;
        _S120 = _S146;
        _S91 = 0.0f;
        _S92 = 0.0f;
        _S93 = 0.0f;
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
        _S121 = _S66;
        _S110 = 0.0f;
        _S111 = 0.0f;
        _S112 = 0.0f;
        _S113 = 0.0f;
        _S114 = 0.0f;
        _S115 = 0.0f;
        _S116 = 0.0f;
        _S117 = 0.0f;
        _S118 = 0.0f;
    }
    if(is_3dgs_1)
    {
        float _S164 = erank_reg_weight_s3_1 * (*_s_dOut_1)[int(3)];
        float _S165 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S166;
        (&_S166)->primal_0 = _S91;
        (&_S166)->differential_0 = 0.0f;
        DiffPair_float_0 _S167;
        (&_S167)->primal_0 = 0.0f;
        (&_S167)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S166, &_S167, _S165);
        float _S168 = - _S166.differential_0;
        DiffPair_float_0 _S169;
        (&_S169)->primal_0 = _S92;
        (&_S169)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S169, _S168);
        DiffPair_float_0 _S170;
        (&_S170)->primal_0 = _S93;
        (&_S170)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S170, _S169.differential_0);
        float _S171 = - _S170.differential_0;
        float _S172 = _S94 * _S171;
        float _S173 = _S95 * _S171;
        DiffPair_float_0 _S174;
        (&_S174)->primal_0 = _S94;
        (&_S174)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S174, _S172);
        float _S175 = _S96 * _S171;
        float _S176 = _S97 * _S171;
        DiffPair_float_0 _S177;
        (&_S177)->primal_0 = _S96;
        (&_S177)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S177, _S175);
        float _S178 = _S98 * _S170.differential_0;
        float _S179 = _S99 * _S170.differential_0;
        DiffPair_float_0 _S180;
        (&_S180)->primal_0 = _S100;
        (&_S180)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S180, _S178);
        float _S181 = - _S179;
        float _S182 = - (_S176 + _S177.differential_0);
        float _S183 = (_S164 + _S173 + _S174.differential_0 + _S182) / _S101;
        float _S184 = _S102 * - _S183;
        float _S185 = _S103 * _S183;
        DiffPair_float_0 _S186;
        (&_S186)->primal_0 = _S104;
        (&_S186)->differential_0 = 0.0f;
        DiffPair_float_0 _S187;
        (&_S187)->primal_0 = _S105;
        (&_S187)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S186, &_S187, _S185);
        DiffPair_float_0 _S188;
        (&_S188)->primal_0 = _S106;
        (&_S188)->differential_0 = 0.0f;
        DiffPair_float_0 _S189;
        (&_S189)->primal_0 = _S107;
        (&_S189)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S188, &_S189, _S186.differential_0);
        float _S190 = (_S180.differential_0 + _S181 + _S182) / _S101;
        float _S191 = _S108 * - _S190;
        float _S192 = _S103 * _S190;
        DiffPair_float_0 _S193;
        (&_S193)->primal_0 = _S109;
        (&_S193)->differential_0 = 0.0f;
        DiffPair_float_0 _S194;
        (&_S194)->primal_0 = _S105;
        (&_S194)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S193, &_S194, _S192);
        DiffPair_float_0 _S195;
        (&_S195)->primal_0 = _S106;
        (&_S195)->differential_0 = 0.0f;
        DiffPair_float_0 _S196;
        (&_S196)->primal_0 = _S107;
        (&_S196)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S195, &_S196, _S193.differential_0);
        float _S197 = _S184 + _S191;
        float3  _S198 = make_float3 (_S188.differential_0 + _S195.differential_0 + _S197, _S189.differential_0 + _S196.differential_0 + _S197, _S187.differential_0 + _S194.differential_0 + _S197);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S199;
        (&_S199)->primal_0 = _S121;
        (&_S199)->differential_0 = _S66;
        s_bwd_prop_exp_2(&_S199, _S198);
        float3  _S200 = make_float3 (2.0f) * _S199.differential_0;
        float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S201;
        (&_S201)->primal_0 = _S110;
        (&_S201)->differential_0 = 0.0f;
        DiffPair_float_0 _S202;
        (&_S202)->primal_0 = max_gauss_ratio_1;
        (&_S202)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S201, &_S202, s_diff_scale_reg_T_0);
        float _S203 = _S201.differential_0 / _S111;
        float _S204 = _S112 * - _S203;
        float _S205 = _S113 * _S203;
        DiffPair_float_0 _S206;
        (&_S206)->primal_0 = _S114;
        (&_S206)->differential_0 = 0.0f;
        DiffPair_float_0 _S207;
        (&_S207)->primal_0 = _S115;
        (&_S207)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S206, &_S207, _S204);
        DiffPair_float_0 _S208;
        (&_S208)->primal_0 = _S116;
        (&_S208)->differential_0 = 0.0f;
        DiffPair_float_0 _S209;
        (&_S209)->primal_0 = _S117;
        (&_S209)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S208, &_S209, _S206.differential_0);
        DiffPair_float_0 _S210;
        (&_S210)->primal_0 = _S118;
        (&_S210)->differential_0 = 0.0f;
        DiffPair_float_0 _S211;
        (&_S211)->primal_0 = _S115;
        (&_S211)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S210, &_S211, _S205);
        DiffPair_float_0 _S212;
        (&_S212)->primal_0 = _S116;
        (&_S212)->differential_0 = 0.0f;
        DiffPair_float_0 _S213;
        (&_S213)->primal_0 = _S117;
        (&_S213)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S212, &_S213, _S210.differential_0);
        float _S214 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * (*_s_dOut_1)[int(1)]);
        float3  _S215 = make_float3 (_S208.differential_0 + _S212.differential_0 + _S214, _S209.differential_0 + _S213.differential_0 + _S214, _S207.differential_0 + _S211.differential_0 + _S214);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S216;
        (&_S216)->primal_0 = _S63.primal_0;
        (&_S216)->differential_0 = _S66;
        s_bwd_prop_exp_2(&_S216, _S215);
        float3  _S217 = _S200 + _S216.differential_0;
        _S71 = (*_s_dOut_1)[int(4)];
        _S72 = (*_s_dOut_1)[int(0)];
        _S121 = _S217;
    }
    else
    {
        float _S218 = erank_reg_weight_1 * (*_s_dOut_1)[int(3)];
        DiffPair_float_0 _S219;
        (&_S219)->primal_0 = _S71;
        (&_S219)->differential_0 = 0.0f;
        DiffPair_float_0 _S220;
        (&_S220)->primal_0 = 0.0f;
        (&_S220)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S219, &_S220, _S218);
        float _S221 = - _S219.differential_0;
        DiffPair_float_0 _S222;
        (&_S222)->primal_0 = _S72;
        (&_S222)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S222, _S221);
        DiffPair_float_0 _S223;
        (&_S223)->primal_0 = _S73;
        (&_S223)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S223, _S222.differential_0);
        float _S224 = - _S223.differential_0;
        float _S225 = _S74 * _S224;
        float _S226 = _S75 * _S224;
        DiffPair_float_0 _S227;
        (&_S227)->primal_0 = _S74;
        (&_S227)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S227, _S225);
        float _S228 = _S76 * _S223.differential_0;
        float _S229 = _S77 * _S223.differential_0;
        DiffPair_float_0 _S230;
        (&_S230)->primal_0 = _S78;
        (&_S230)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S230, _S228);
        float _S231 = - _S229;
        float _S232 = (_S226 + _S227.differential_0) / _S79;
        float _S233 = _S80 * - _S232;
        float _S234 = _S81 * _S232;
        DiffPair_float_0 _S235;
        (&_S235)->primal_0 = _S82;
        (&_S235)->differential_0 = 0.0f;
        DiffPair_float_0 _S236;
        (&_S236)->primal_0 = _S83;
        (&_S236)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S235, &_S236, _S234);
        float _S237 = (_S230.differential_0 + _S231) / _S79;
        float _S238 = _S84 * - _S237;
        float _S239 = _S81 * _S237;
        DiffPair_float_0 _S240;
        (&_S240)->primal_0 = _S82;
        (&_S240)->differential_0 = 0.0f;
        DiffPair_float_0 _S241;
        (&_S241)->primal_0 = _S83;
        (&_S241)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S240, &_S241, _S239);
        float _S242 = _S233 + _S238;
        float2  _S243 = make_float2 (_S235.differential_0 + _S240.differential_0 + _S242, _S236.differential_0 + _S241.differential_0 + _S242);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S244;
        (&_S244)->primal_0 = _S119;
        (&_S244)->differential_0 = _S65;
        s_bwd_prop_exp_1(&_S244, _S243);
        float2  _S245 = make_float2 (2.0f) * _S244.differential_0;
        float s_diff_scale_reg_T_1 = scale_regularization_weight_1 * (*_s_dOut_1)[int(2)];
        DiffPair_float_0 _S246;
        (&_S246)->primal_0 = _S85;
        (&_S246)->differential_0 = 0.0f;
        DiffPair_float_0 _S247;
        (&_S247)->primal_0 = max_gauss_ratio_1;
        (&_S247)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S246, &_S247, s_diff_scale_reg_T_1);
        float _S248 = _S246.differential_0 / _S86;
        float _S249 = _S87 * - _S248;
        float _S250 = _S88 * _S248;
        DiffPair_float_0 _S251;
        (&_S251)->primal_0 = _S89;
        (&_S251)->differential_0 = 0.0f;
        DiffPair_float_0 _S252;
        (&_S252)->primal_0 = _S90;
        (&_S252)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S251, &_S252, _S249);
        DiffPair_float_0 _S253;
        (&_S253)->primal_0 = _S89;
        (&_S253)->differential_0 = 0.0f;
        DiffPair_float_0 _S254;
        (&_S254)->primal_0 = _S90;
        (&_S254)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S253, &_S254, _S250);
        float _S255 = mcmc_scale_reg_weight_1 * (0.5f * (*_s_dOut_1)[int(1)]);
        float2  _S256 = make_float2 (_S251.differential_0 + _S253.differential_0 + _S255, _S252.differential_0 + _S254.differential_0 + _S255);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S257;
        (&_S257)->primal_0 = _S120;
        (&_S257)->differential_0 = _S65;
        s_bwd_prop_exp_1(&_S257, _S256);
        float2  _S258 = _S245 + _S257.differential_0;
        float3  _S259 = make_float3 (_S258.x, _S258.y, 0.0f);
        _S71 = (*_s_dOut_1)[int(4)];
        _S72 = (*_s_dOut_1)[int(0)];
        _S121 = _S259;
    }
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * _S71;
    float _S260 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S261;
    (&_S261)->primal_0 = _S70;
    (&_S261)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S261, _S260);
    float _S262 = _S261.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S263 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S264;
    (&_S264)->primal_0 = _S64.primal_0;
    (&_S264)->differential_0 = _S263;
    s_bwd_length_impl_0(&_S264, _S262);
    float _S265 = - (mcmc_opacity_reg_weight_1 * _S72 / _S69);
    DiffPair_float_0 _S266;
    (&_S266)->primal_0 = _S67;
    (&_S266)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S266, _S265);
    float _S267 = - _S266.differential_0;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S264.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S267;
    dpscales_0->primal_0 = (*dpscales_0).primal_0;
    dpscales_0->differential_0 = _S121;
    return;
}

inline __device__ void s_bwd_per_splat_losses_0(bool _S268, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S269, DiffPair_float_0 * _S270, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S271, float _S272, float _S273, float _S274, float _S275, float _S276, float _S277, float _S278, FixedArray<float, 5>  * _S279)
{
    s_bwd_prop_per_splat_losses_0(_S268, _S269, _S270, _S271, _S272, _S273, _S274, _S275, _S276, _S277, _S278, _S279);
    return;
}

inline __device__ void per_splat_losses_bwd(bool is_3dgs_2, float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  * v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float3  _S280 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_scales_0;
    (&p_scales_0)->primal_0 = scales_1;
    (&p_scales_0)->differential_0 = _S280;
    DiffPair_float_0 p_opacity_0;
    (&p_opacity_0)->primal_0 = opacity_1;
    (&p_opacity_0)->differential_0 = 0.0f;
    float4  _S281 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 p_quat_0;
    (&p_quat_0)->primal_0 = quat_1;
    (&p_quat_0)->differential_0 = _S281;
    s_bwd_per_splat_losses_0(is_3dgs_2, &p_scales_0, &p_opacity_0, &p_quat_0, mcmc_opacity_reg_weight_2, mcmc_scale_reg_weight_2, max_gauss_ratio_2, scale_regularization_weight_2, erank_reg_weight_2, erank_reg_weight_s3_2, quat_norm_reg_weight_2, v_loss_0);
    *v_scales_0 = p_scales_0.differential_0;
    *v_opacity_0 = p_opacity_0.differential_0;
    *v_quat_0 = p_quat_0.differential_0;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float3  dOut_8)
{
    float3  _S282 = _slang_select(((*dpx_9).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_9).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_8;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S282;
    return;
}

inline __device__ float3  abs_0(float3  x_12)
{
    float3  result_8;
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
        *_slang_vector_get_element_ptr(&result_8, i_5) = (F32_abs((_slang_vector_get_element(x_12, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_8;
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_10, float dOut_9)
{
    float _S283 = -0.5f / ((*dpx_10).primal_0 * (F32_sqrt(((*dpx_10).primal_0)))) * dOut_9;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S283;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_10)
{
    DiffPair_float_0 _S284 = *dpx_11;
    bool _S285;
    if(((*dpx_11).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S285 = ((*dpx_11).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S285 = false;
    }
    float _S286;
    if(_S285)
    {
        _S286 = dOut_10;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpx_11->primal_0 = _S284.primal_0;
    dpx_11->differential_0 = _S286;
    DiffPair_float_0 _S287 = *dpMin_0;
    if((_S284.primal_0) < ((*dpMin_0).primal_0))
    {
        _S286 = dOut_10;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpMin_0->primal_0 = _S287.primal_0;
    dpMin_0->differential_0 = _S286;
    DiffPair_float_0 _S288 = *dpMax_0;
    if(((*dpx_11).primal_0) > ((*dpMax_0).primal_0))
    {
        _S286 = dOut_10;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpMax_0->primal_0 = _S288.primal_0;
    dpMax_0->differential_0 = _S286;
    return;
}

inline __device__ float clamp_0(float x_13, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_13), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_11)
{
    float _S289 = (1.0f - (*dps_0).primal_0) * dOut_11;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S289;
    DiffPair_float_0 _S290 = *dpy_3;
    float _S291 = (*dps_0).primal_0 * dOut_11;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S291;
    float _S292 = (_S290.primal_0 - (*dpx_12).primal_0) * dOut_11;
    dps_0->primal_0 = _S290.primal_0;
    dps_0->differential_0 = _S292;
    return;
}

inline __device__ float lerp_0(float x_14, float y_7, float s_4)
{
    return x_14 + (y_7 - x_14) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, FixedArray<float, 21>  * _S293)
{
    float3  _S294;
    bool _S295;
    bool _S296;
    FixedArray<float, 21>  losses_1;
    float _S297 = float(mask_0);
    float3  _S298 = ref_rgb_0 - render_rgb_0;
    float3  _S299 = abs_0(_S298);
    losses_1[int(0)] = _S297 * ((_S299.x + _S299.y + _S299.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S297 * (dot_0(_S298, _S298) * 0.3333333432674408f);
    float _S300 = float(depth_mask_0 & mask_0);
    float _S301 = _S300 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S302 = _S300 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S301;
    losses_1[int(3)] = _S302;
    losses_1[int(4)] = _S301 * _S301;
    losses_1[int(5)] = _S302 * _S302;
    losses_1[int(6)] = _S301 * _S302;
    bool _S303 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S304 = norm2_0 == 0.0f;
        _S295 = _S304;
        if(_S304)
        {
            _S294 = make_float3 (0.0f);
            break;
        }
        _S294 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S305;
    bool _S306 = !_S295;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S307 = norm2_1 == 0.0f;
        _S296 = _S307;
        if(_S307)
        {
            _S305 = make_float3 (0.0f);
            break;
        }
        _S305 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S308;
    float3  _S309;
    bool _S310 = !_S296;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S309 = make_float3 (0.0f);
            _S308 = false;
            break;
        }
        _S309 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S308 = _S303;
        break;
    }
    float _S311 = float(_S306 & _S308);
    float3  _S312 = abs_0(_S309 - _S294);
    losses_1[int(7)] = _S311 * ((_S312.x + _S312.y + _S312.z) * 0.3333333432674408f);
    float _S313 = float(_S310 & _S308);
    float3  _S314 = abs_0(_S309 - _S305);
    losses_1[int(8)] = _S313 * ((_S314.x + _S314.y + _S314.z) * 0.3333333432674408f);
    float _S315 = float(_S306 & _S310);
    float3  _S316 = abs_0(_S305 - _S294);
    losses_1[int(11)] = _S315 * ((_S316.x + _S316.y + _S316.z) * 0.3333333432674408f);
    float _S317 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S318 = float(ref_alpha_0);
    float _S319 = (F32_max((_S317), (_S318)));
    losses_1[int(9)] = - lerp_0((F32_log(((F32_max((1.0f - _S319), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S319), (9.99999997475242708e-07f)))))), _S318);
    float _S320 = 1.0f - _S317;
    float _S321 = 1.0f - _S318;
    float _S322 = (F32_max((_S320), (_S321)));
    losses_1[int(10)] = - lerp_0((F32_log(((F32_max((1.0f - _S322), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S322), (9.99999997475242708e-07f)))))), _S321);
    losses_1[int(12)] = 4.0f * _S317 * _S320;
    float _S323 = (F32_max((_S317), (9.999999960041972e-13f)));
    losses_1[int(13)] = (rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f / _S323;
    losses_1[int(14)] = depth_dist_0 / _S323;
    losses_1[int(15)] = (normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f / _S323;
    losses_1[int(16)] = _S297;
    losses_1[int(17)] = _S300;
    losses_1[int(18)] = _S311;
    losses_1[int(19)] = _S313;
    losses_1[int(20)] = _S315;
    *_S293 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S324, float3  _S325)
{
    return dot_0(_S324, _S325);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S326)
{
    return (F32_rsqrt((_S326)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S327, float _S328, float _S329)
{
    return clamp_0(_S327, _S328, _S329);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S330, DiffPair_float_0 * _S331, DiffPair_float_0 * _S332, float _S333)
{
    _d_lerp_0(_S330, _S331, _S332, _S333);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S334, DiffPair_float_0 * _S335, DiffPair_float_0 * _S336, float _S337)
{
    _d_clamp_0(_S334, _S335, _S336, _S337);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S338, float3  _S339)
{
    _d_abs_vector_0(_S338, _S339);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S340, float _S341)
{
    _d_rsqrt_0(_S340, _S341);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S342, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S343, float _S344)
{
    _d_dot_0(_S342, _S343, _S344);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, FixedArray<float, 21>  * _s_dOut_2)
{
    DiffPair_float_0 _S345 = *dprender_depth_0;
    DiffPair_float_0 _S346 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S347 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S348 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S349 = *dpref_normal_0;
    DiffPair_float_0 _S350 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S351 = *dprgb_dist_0;
    DiffPair_float_0 _S352 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S353 = *dpnormal_dist_0;
    float3  _S354 = make_float3 (0.0f);
    float _S355 = float(mask_1);
    float3  _S356 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S357 = float(depth_mask_1 & mask_1);
    float _S358 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S359 = _S357 * s_primal_ctx_log_0(_S358);
    float _S360 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S361 = _S357 * s_primal_ctx_log_0(_S360);
    bool _S362 = normal_mask_1 & mask_1;
    float _S363 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S364 = _S363 == 0.0f;
    float3  _S365;
    if(_S364)
    {
        _S365 = make_float3 (0.0f);
    }
    bool _S366 = !_S364;
    float3  _S367;
    if(_S366)
    {
        float _S368 = s_primal_ctx_rsqrt_0(_S363);
        float3  _S369 = make_float3 (_S368);
        _S365 = _S347.primal_0 * make_float3 (_S368);
        _S367 = _S369;
    }
    else
    {
        _S367 = _S354;
    }
    float _S370 = s_primal_ctx_dot_0(_S348.primal_0, _S348.primal_0);
    bool _S371 = _S370 == 0.0f;
    float3  _S372;
    if(_S371)
    {
        _S372 = make_float3 (0.0f);
    }
    bool _S373 = !_S371;
    float3  _S374;
    if(_S373)
    {
        float _S375 = s_primal_ctx_rsqrt_0(_S370);
        float3  _S376 = make_float3 (_S375);
        _S372 = _S348.primal_0 * make_float3 (_S375);
        _S374 = _S376;
    }
    else
    {
        _S374 = _S354;
    }
    float _S377 = s_primal_ctx_dot_0(_S349.primal_0, _S349.primal_0);
    bool _S378 = _S377 == 0.0f;
    float3  _S379;
    bool _S380;
    if(_S378)
    {
        float3  _S381 = make_float3 (0.0f);
        _S380 = false;
        _S379 = _S381;
    }
    else
    {
        _S380 = _S362;
    }
    bool _S382 = !_S378;
    float3  _S383;
    if(_S382)
    {
        float _S384 = s_primal_ctx_rsqrt_0(_S377);
        float3  _S385 = make_float3 (_S384);
        _S379 = _S349.primal_0 * make_float3 (_S384);
        _S383 = _S385;
    }
    else
    {
        _S383 = _S354;
    }
    float _S386 = float(_S366 & _S380);
    float3  _S387 = _S379 - _S365;
    float _S388 = float(_S373 & _S380);
    float3  _S389 = _S379 - _S372;
    float _S390 = float(_S366 & _S373);
    float3  _S391 = _S372 - _S365;
    float _S392 = s_primal_ctx_clamp_0(_S350.primal_0, 0.0f, 1.0f);
    float _S393 = float(ref_alpha_1);
    float _S394 = s_primal_ctx_max_0(_S392, _S393);
    float _S395 = 1.0f - _S394;
    float _S396 = s_primal_ctx_max_0(_S395, 9.99999997475242708e-07f);
    float _S397 = s_primal_ctx_log_0(_S396);
    float _S398 = s_primal_ctx_max_0(_S394, 9.99999997475242708e-07f);
    float _S399 = s_primal_ctx_log_0(_S398);
    float _S400 = 1.0f - _S392;
    float _S401 = 1.0f - _S393;
    float _S402 = s_primal_ctx_max_0(_S400, _S401);
    float _S403 = 1.0f - _S402;
    float _S404 = s_primal_ctx_max_0(_S403, 9.99999997475242708e-07f);
    float _S405 = s_primal_ctx_log_0(_S404);
    float _S406 = s_primal_ctx_max_0(_S402, 9.99999997475242708e-07f);
    float _S407 = s_primal_ctx_log_0(_S406);
    float _S408 = 4.0f * _S392;
    float _S409 = s_primal_ctx_max_0(_S392, 9.999999960041972e-13f);
    float _S410 = _S409 * _S409;
    float _S411 = (*_s_dOut_2)[int(15)] / _S410;
    float _S412 = 0.3333333432674408f * (_S409 * _S411);
    float _S413 = (*_s_dOut_2)[int(14)] / _S410;
    float _S414 = _S409 * _S413;
    float _S415 = (*_s_dOut_2)[int(13)] / _S410;
    float _S416 = _S409 * _S415;
    float _S417 = (_S353.primal_0.x + _S353.primal_0.y + _S353.primal_0.z) * 0.3333333432674408f * - _S411 + _S352.primal_0 * - _S413 + (_S351.primal_0.x + _S351.primal_0.y + _S351.primal_0.z) * 0.3333333432674408f * - _S415;
    DiffPair_float_0 _S418;
    (&_S418)->primal_0 = _S392;
    (&_S418)->differential_0 = 0.0f;
    DiffPair_float_0 _S419;
    (&_S419)->primal_0 = 9.999999960041972e-13f;
    (&_S419)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S418, &_S419, _S417);
    float _S420 = 0.3333333432674408f * _S416;
    float _S421 = _S408 * (*_s_dOut_2)[int(12)];
    float _S422 = 4.0f * (_S400 * (*_s_dOut_2)[int(12)]);
    float _S423 = - (*_s_dOut_2)[int(10)];
    DiffPair_float_0 _S424;
    (&_S424)->primal_0 = _S405;
    (&_S424)->differential_0 = 0.0f;
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = _S407;
    (&_S425)->differential_0 = 0.0f;
    DiffPair_float_0 _S426;
    (&_S426)->primal_0 = _S401;
    (&_S426)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S424, &_S425, &_S426, _S423);
    DiffPair_float_0 _S427;
    (&_S427)->primal_0 = _S406;
    (&_S427)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S427, _S425.differential_0);
    DiffPair_float_0 _S428;
    (&_S428)->primal_0 = _S402;
    (&_S428)->differential_0 = 0.0f;
    DiffPair_float_0 _S429;
    (&_S429)->primal_0 = 9.99999997475242708e-07f;
    (&_S429)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S428, &_S429, _S427.differential_0);
    DiffPair_float_0 _S430;
    (&_S430)->primal_0 = _S404;
    (&_S430)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S430, _S424.differential_0);
    DiffPair_float_0 _S431;
    (&_S431)->primal_0 = _S403;
    (&_S431)->differential_0 = 0.0f;
    DiffPair_float_0 _S432;
    (&_S432)->primal_0 = 9.99999997475242708e-07f;
    (&_S432)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S431, &_S432, _S430.differential_0);
    float _S433 = _S428.differential_0 + - _S431.differential_0;
    DiffPair_float_0 _S434;
    (&_S434)->primal_0 = _S400;
    (&_S434)->differential_0 = 0.0f;
    DiffPair_float_0 _S435;
    (&_S435)->primal_0 = _S401;
    (&_S435)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S434, &_S435, _S433);
    float _S436 = - (_S421 + _S434.differential_0);
    float _S437 = - (*_s_dOut_2)[int(9)];
    DiffPair_float_0 _S438;
    (&_S438)->primal_0 = _S397;
    (&_S438)->differential_0 = 0.0f;
    DiffPair_float_0 _S439;
    (&_S439)->primal_0 = _S399;
    (&_S439)->differential_0 = 0.0f;
    DiffPair_float_0 _S440;
    (&_S440)->primal_0 = _S393;
    (&_S440)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S438, &_S439, &_S440, _S437);
    DiffPair_float_0 _S441;
    (&_S441)->primal_0 = _S398;
    (&_S441)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S441, _S439.differential_0);
    DiffPair_float_0 _S442;
    (&_S442)->primal_0 = _S394;
    (&_S442)->differential_0 = 0.0f;
    DiffPair_float_0 _S443;
    (&_S443)->primal_0 = 9.99999997475242708e-07f;
    (&_S443)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S442, &_S443, _S441.differential_0);
    DiffPair_float_0 _S444;
    (&_S444)->primal_0 = _S396;
    (&_S444)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S444, _S438.differential_0);
    DiffPair_float_0 _S445;
    (&_S445)->primal_0 = _S395;
    (&_S445)->differential_0 = 0.0f;
    DiffPair_float_0 _S446;
    (&_S446)->primal_0 = 9.99999997475242708e-07f;
    (&_S446)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S445, &_S446, _S444.differential_0);
    float _S447 = _S442.differential_0 + - _S445.differential_0;
    DiffPair_float_0 _S448;
    (&_S448)->primal_0 = _S392;
    (&_S448)->differential_0 = 0.0f;
    DiffPair_float_0 _S449;
    (&_S449)->primal_0 = _S393;
    (&_S449)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S448, &_S449, _S447);
    float _S450 = _S418.differential_0 + _S422 + _S436 + _S448.differential_0;
    DiffPair_float_0 _S451;
    (&_S451)->primal_0 = _S350.primal_0;
    (&_S451)->differential_0 = 0.0f;
    DiffPair_float_0 _S452;
    (&_S452)->primal_0 = 0.0f;
    (&_S452)->differential_0 = 0.0f;
    DiffPair_float_0 _S453;
    (&_S453)->primal_0 = 1.0f;
    (&_S453)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S451, &_S452, &_S453, _S450);
    DiffPair_float_0 _S454 = _S451;
    float _S455 = 0.3333333432674408f * (_S390 * (*_s_dOut_2)[int(11)]);
    float3  _S456 = make_float3 (_S455, _S455, _S455);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S457;
    (&_S457)->primal_0 = _S391;
    (&_S457)->differential_0 = _S354;
    s_bwd_prop_abs_0(&_S457, _S456);
    float3  _S458 = - _S457.differential_0;
    float _S459 = 0.3333333432674408f * (_S388 * (*_s_dOut_2)[int(8)]);
    float3  _S460 = make_float3 (_S459, _S459, _S459);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S461;
    (&_S461)->primal_0 = _S389;
    (&_S461)->differential_0 = _S354;
    s_bwd_prop_abs_0(&_S461, _S460);
    float3  _S462 = - _S461.differential_0;
    float _S463 = 0.3333333432674408f * (_S386 * (*_s_dOut_2)[int(7)]);
    float3  _S464 = make_float3 (_S463, _S463, _S463);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S465;
    (&_S465)->primal_0 = _S387;
    (&_S465)->differential_0 = _S354;
    s_bwd_prop_abs_0(&_S465, _S464);
    float3  _S466 = _S461.differential_0 + _S465.differential_0;
    float3  _S467 = _S458 + - _S465.differential_0;
    float3  _S468 = make_float3 (_S412, _S412, _S412);
    float3  _S469 = make_float3 (_S420, _S420, _S420);
    float3  _S470 = _S457.differential_0 + _S462;
    float _S471;
    if(_S382)
    {
        float3  _S472 = _S349.primal_0 * _S466;
        float3  _S473 = _S383 * _S466;
        float _S474 = _S472.x + _S472.y + _S472.z;
        DiffPair_float_0 _S475;
        (&_S475)->primal_0 = _S377;
        (&_S475)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S475, _S474);
        _S471 = _S475.differential_0;
        _S365 = _S473;
    }
    else
    {
        _S471 = 0.0f;
        _S365 = _S354;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S476;
    (&_S476)->primal_0 = _S349.primal_0;
    (&_S476)->differential_0 = _S354;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S477;
    (&_S477)->primal_0 = _S349.primal_0;
    (&_S477)->differential_0 = _S354;
    s_bwd_prop_dot_0(&_S476, &_S477, _S471);
    float3  _S478 = _S477.differential_0 + _S476.differential_0 + _S365;
    if(_S373)
    {
        float3  _S479 = _S348.primal_0 * _S470;
        float3  _S480 = _S374 * _S470;
        float _S481 = _S479.x + _S479.y + _S479.z;
        DiffPair_float_0 _S482;
        (&_S482)->primal_0 = _S370;
        (&_S482)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S482, _S481);
        _S471 = _S482.differential_0;
        _S365 = _S480;
    }
    else
    {
        _S471 = 0.0f;
        _S365 = _S354;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S483;
    (&_S483)->primal_0 = _S348.primal_0;
    (&_S483)->differential_0 = _S354;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S348.primal_0;
    (&_S484)->differential_0 = _S354;
    s_bwd_prop_dot_0(&_S483, &_S484, _S471);
    float3  _S485 = _S484.differential_0 + _S483.differential_0 + _S365;
    if(_S366)
    {
        float3  _S486 = _S347.primal_0 * _S467;
        float3  _S487 = _S367 * _S467;
        float _S488 = _S486.x + _S486.y + _S486.z;
        DiffPair_float_0 _S489;
        (&_S489)->primal_0 = _S363;
        (&_S489)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S489, _S488);
        _S471 = _S489.differential_0;
        _S365 = _S487;
    }
    else
    {
        _S471 = 0.0f;
        _S365 = _S354;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S490;
    (&_S490)->primal_0 = _S347.primal_0;
    (&_S490)->differential_0 = _S354;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S491;
    (&_S491)->primal_0 = _S347.primal_0;
    (&_S491)->differential_0 = _S354;
    s_bwd_prop_dot_0(&_S490, &_S491, _S471);
    float _S492 = _S361 * (*_s_dOut_2)[int(6)];
    float _S493 = _S361 * (*_s_dOut_2)[int(5)];
    float _S494 = _S359 * (*_s_dOut_2)[int(4)];
    float _S495 = _S357 * (_S359 * (*_s_dOut_2)[int(6)] + _S493 + _S493 + (*_s_dOut_2)[int(3)]);
    DiffPair_float_0 _S496;
    (&_S496)->primal_0 = _S360;
    (&_S496)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S496, _S495);
    DiffPair_float_0 _S497;
    (&_S497)->primal_0 = _S346.primal_0;
    (&_S497)->differential_0 = 0.0f;
    DiffPair_float_0 _S498;
    (&_S498)->primal_0 = 0.00009999999747379f;
    (&_S498)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S497, &_S498, _S496.differential_0);
    float _S499 = _S357 * (_S492 + _S494 + _S494 + (*_s_dOut_2)[int(2)]);
    DiffPair_float_0 _S500;
    (&_S500)->primal_0 = _S358;
    (&_S500)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S500, _S499);
    DiffPair_float_0 _S501;
    (&_S501)->primal_0 = _S345.primal_0;
    (&_S501)->differential_0 = 0.0f;
    DiffPair_float_0 _S502;
    (&_S502)->primal_0 = 0.00009999999747379f;
    (&_S502)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S501, &_S502, _S500.differential_0);
    float _S503 = 0.3333333432674408f * (_S355 * (*_s_dOut_2)[int(1)]);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S504;
    (&_S504)->primal_0 = _S356;
    (&_S504)->differential_0 = _S354;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S505;
    (&_S505)->primal_0 = _S356;
    (&_S505)->differential_0 = _S354;
    s_bwd_prop_dot_0(&_S504, &_S505, _S503);
    float _S506 = 0.3333333432674408f * (_S355 * (*_s_dOut_2)[int(0)]);
    float3  _S507 = make_float3 (_S506, _S506, _S506);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S508;
    (&_S508)->primal_0 = _S356;
    (&_S508)->differential_0 = _S354;
    s_bwd_prop_abs_0(&_S508, _S507);
    float3  _S509 = _S505.differential_0 + _S504.differential_0 + _S508.differential_0;
    float3  _S510 = - _S509;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S468;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S414;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S469;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S454.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S478;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S485;
    float3  _S511 = _S491.differential_0 + _S490.differential_0 + _S365;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S511;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S497.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S501.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S509;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S510;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S512, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S513, DiffPair_float_0 * _S514, DiffPair_float_0 * _S515, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S516, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S517, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S518, DiffPair_float_0 * _S519, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S520, DiffPair_float_0 * _S521, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S522, bool _S523, bool _S524, bool _S525, bool _S526, FixedArray<float, 21>  * _S527)
{
    s_bwd_prop_per_pixel_losses_0(_S512, _S513, _S514, _S515, _S516, _S517, _S518, _S519, _S520, _S521, _S522, _S523, _S524, _S525, _S526, _S527);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_0, FixedArray<float, 21>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S528 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S528;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S528;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S528;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S528;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_alpha_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, mask_2, depth_mask_2, normal_mask_2, alpha_mask_0, v_losses_0);
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_13, float dOut_12)
{
    float _S529 = 1.0f / ((*dpx_13).primal_0 * 52.30258560180664062f) * dOut_12;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S529;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 21>  * raw_losses_0, float num_pixels_0, FixedArray<float, 10>  * weights_0, FixedArray<float, 10>  * _S530)
{
    FixedArray<float, 10>  losses_2;
    float _S531 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(0)] = (*weights_0)[int(0)] * (*raw_losses_0)[int(0)] / _S531;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S531)));
    float _S532;
    if(((*raw_losses_0)[int(17)]) > 0.0f)
    {
        _S532 = (*weights_0)[int(1)] * (1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(17)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(17)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(17)]) + 1.0f)))))));
    }
    else
    {
        _S532 = 0.0f;
    }
    losses_2[int(2)] = _S532;
    losses_2[int(3)] = (*weights_0)[int(2)] * ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(18)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(18)]) > 0.5f) + int(((*raw_losses_0)[int(19)]) > 0.5f)), (int(1)))));
    float _S533 = (F32_max((num_pixels_0), (1.0f)));
    losses_2[int(4)] = ((*weights_0)[int(3)] * (*raw_losses_0)[int(9)] + (*weights_0)[int(4)] * (*raw_losses_0)[int(10)]) / _S533;
    losses_2[int(5)] = (*weights_0)[int(5)] * (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)));
    losses_2[int(6)] = (*weights_0)[int(6)] * (*raw_losses_0)[int(12)] / _S533;
    losses_2[int(7)] = (*weights_0)[int(7)] * (*raw_losses_0)[int(13)] / _S533;
    losses_2[int(8)] = (*weights_0)[int(8)] * (*raw_losses_0)[int(14)] / _S533;
    losses_2[int(9)] = (*weights_0)[int(9)] * (*raw_losses_0)[int(15)] / _S533;
    *_S530 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C21x3E_0
{
    FixedArray<float, 21>  primal_0;
    FixedArray<float, 21>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S534)
{
    return (F32_sqrt((_S534)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S535, float _S536)
{
    _d_log10_0(_S535, _S536);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C21x3E_0 * dpraw_losses_0, float num_pixels_1, FixedArray<float, 10>  * weights_1, FixedArray<float, 10>  * _s_dOut_3)
{
    FixedArray<float, 21>  _S537 = dpraw_losses_0->primal_0;
    float _S538 = (*weights_1)[int(0)] * dpraw_losses_0->primal_0[int(0)];
    float _S539 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(16)], 1.0f);
    float _S540 = _S539 * _S539;
    float _S541 = dpraw_losses_0->primal_0[int(1)] / _S539;
    bool _S542 = (dpraw_losses_0->primal_0[int(17)]) > 0.0f;
    float _S543;
    float _S544;
    float _S545;
    float _S546;
    float _S547;
    float _S548;
    float _S549;
    float _S550;
    float _S551;
    float _S552;
    float _S553;
    float _S554;
    float _S555;
    float _S556;
    if(_S542)
    {
        float _S557 = _S537[int(2)] * _S537[int(3)];
        float _S558 = _S537[int(17)] * _S537[int(17)];
        float _S559 = _S537[int(6)] - _S557 / _S537[int(17)];
        float _S560 = _S537[int(2)] * _S537[int(2)];
        float _S561 = _S537[int(4)] - _S560 / _S537[int(17)];
        float _S562 = _S537[int(3)] * _S537[int(3)];
        float _S563 = _S537[int(5)] - _S562 / _S537[int(17)];
        float _S564 = _S561 * _S563 + 1.0f;
        float _S565 = s_primal_ctx_max_0(9.999999960041972e-13f, _S564);
        float _S566 = s_primal_ctx_sqrt_0(_S565);
        float _S567 = _S566 * _S566;
        _S543 = (*weights_1)[int(1)];
        _S544 = _S567;
        _S545 = _S559;
        _S546 = _S566;
        _S547 = _S565;
        _S548 = _S564;
        _S549 = _S561;
        _S550 = _S563;
        _S551 = _S558;
        _S552 = _S562;
        _S553 = _S537[int(3)];
        _S554 = _S560;
        _S555 = _S537[int(2)];
        _S556 = _S557;
    }
    else
    {
        _S543 = 0.0f;
        _S544 = 0.0f;
        _S545 = 0.0f;
        _S546 = 0.0f;
        _S547 = 0.0f;
        _S548 = 0.0f;
        _S549 = 0.0f;
        _S550 = 0.0f;
        _S551 = 0.0f;
        _S552 = 0.0f;
        _S553 = 0.0f;
        _S554 = 0.0f;
        _S555 = 0.0f;
        _S556 = 0.0f;
    }
    float _S568 = s_primal_ctx_max_0(_S537[int(18)], 1.0f);
    float _S569 = _S568 * _S568;
    float _S570 = s_primal_ctx_max_0(_S537[int(19)], 1.0f);
    float _S571 = _S570 * _S570;
    float _S572 = float((I32_max((int((_S537[int(18)]) > 0.5f) + int((_S537[int(19)]) > 0.5f)), (int(1)))));
    float _S573 = s_primal_ctx_max_0(num_pixels_1, 1.0f);
    float _S574 = _S573 * _S573;
    float _S575 = (*weights_1)[int(5)] * _S537[int(11)];
    float _S576 = s_primal_ctx_max_0(_S537[int(20)], 1.0f);
    float _S577 = _S576 * _S576;
    float _S578 = (*weights_1)[int(9)] * (_S573 * ((*_s_dOut_3)[int(9)] / _S574));
    float _S579 = (*weights_1)[int(8)] * (_S573 * ((*_s_dOut_3)[int(8)] / _S574));
    float _S580 = (*weights_1)[int(7)] * (_S573 * ((*_s_dOut_3)[int(7)] / _S574));
    float _S581 = (*weights_1)[int(6)] * (_S573 * ((*_s_dOut_3)[int(6)] / _S574));
    float _S582 = (*_s_dOut_3)[int(5)] / _S577;
    float _S583 = _S575 * - _S582;
    float _S584 = _S576 * _S582;
    DiffPair_float_0 _S585;
    (&_S585)->primal_0 = _S537[int(20)];
    (&_S585)->differential_0 = 0.0f;
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = 1.0f;
    (&_S586)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S585, &_S586, _S583);
    float _S587 = (*weights_1)[int(5)] * _S584;
    float _S588 = _S573 * ((*_s_dOut_3)[int(4)] / _S574);
    float _S589 = (*weights_1)[int(4)] * _S588;
    float _S590 = (*weights_1)[int(3)] * _S588;
    float _S591 = (*weights_1)[int(2)] * ((*_s_dOut_3)[int(3)] / _S572);
    float _S592 = _S591 / _S571;
    float _S593 = _S537[int(8)] * - _S592;
    float _S594 = _S570 * _S592;
    DiffPair_float_0 _S595;
    (&_S595)->primal_0 = _S537[int(19)];
    (&_S595)->differential_0 = 0.0f;
    DiffPair_float_0 _S596;
    (&_S596)->primal_0 = 1.0f;
    (&_S596)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S595, &_S596, _S593);
    float _S597 = _S591 / _S569;
    float _S598 = _S537[int(7)] * - _S597;
    float _S599 = _S568 * _S597;
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = _S537[int(18)];
    (&_S600)->differential_0 = 0.0f;
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = 1.0f;
    (&_S601)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S600, &_S601, _S598);
    FixedArray<float, 21>  _S602;
    _S602[int(0)] = 0.0f;
    _S602[int(1)] = 0.0f;
    _S602[int(2)] = 0.0f;
    _S602[int(3)] = 0.0f;
    _S602[int(4)] = 0.0f;
    _S602[int(5)] = 0.0f;
    _S602[int(6)] = 0.0f;
    _S602[int(7)] = 0.0f;
    _S602[int(8)] = 0.0f;
    _S602[int(9)] = 0.0f;
    _S602[int(10)] = 0.0f;
    _S602[int(11)] = 0.0f;
    _S602[int(12)] = 0.0f;
    _S602[int(13)] = 0.0f;
    _S602[int(14)] = 0.0f;
    _S602[int(15)] = 0.0f;
    _S602[int(16)] = 0.0f;
    _S602[int(17)] = 0.0f;
    _S602[int(18)] = 0.0f;
    _S602[int(19)] = 0.0f;
    _S602[int(20)] = 0.0f;
    _S602[int(15)] = _S578;
    _S602[int(14)] = _S579;
    _S602[int(13)] = _S580;
    _S602[int(12)] = _S581;
    _S602[int(20)] = _S585.differential_0;
    _S602[int(11)] = _S587;
    _S602[int(10)] = _S589;
    _S602[int(9)] = _S590;
    _S602[int(19)] = _S595.differential_0;
    _S602[int(8)] = _S594;
    _S602[int(18)] = _S600.differential_0;
    _S602[int(7)] = _S599;
    float _S603 = _S602[int(0)];
    float _S604 = _S602[int(1)];
    float _S605 = _S602[int(2)];
    float _S606 = _S602[int(3)];
    float _S607 = _S602[int(4)];
    float _S608 = _S602[int(5)];
    float _S609 = _S602[int(6)];
    float _S610 = _S602[int(7)];
    float _S611 = _S602[int(8)];
    float _S612 = _S602[int(9)];
    float _S613 = _S602[int(10)];
    float _S614 = _S602[int(11)];
    float _S615 = _S602[int(12)];
    float _S616 = _S602[int(13)];
    float _S617 = _S602[int(14)];
    float _S618 = _S602[int(15)];
    float _S619 = _S602[int(16)];
    float _S620 = _S602[int(17)];
    float _S621 = _S602[int(18)];
    float _S622 = _S602[int(19)];
    float _S623 = _S602[int(20)];
    FixedArray<float, 21>  _S624;
    if(_S542)
    {
        float _S625 = - (_S543 * (*_s_dOut_3)[int(2)]) / _S544;
        float _S626 = _S545 * - _S625;
        float _S627 = _S546 * _S625;
        DiffPair_float_0 _S628;
        (&_S628)->primal_0 = _S547;
        (&_S628)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S628, _S626);
        DiffPair_float_0 _S629;
        (&_S629)->primal_0 = 9.999999960041972e-13f;
        (&_S629)->differential_0 = 0.0f;
        DiffPair_float_0 _S630;
        (&_S630)->primal_0 = _S548;
        (&_S630)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S629, &_S630, _S628.differential_0);
        float _S631 = _S549 * _S630.differential_0;
        float _S632 = _S550 * _S630.differential_0;
        float _S633 = - _S631 / _S551;
        float _S634 = _S553 * (_S537[int(17)] * _S633);
        float _S635 = - _S632 / _S551;
        float _S636 = _S555 * (_S537[int(17)] * _S635);
        float _S637 = - _S627 / _S551;
        float _S638 = _S537[int(17)] * _S637;
        float _S639 = _S634 + _S634 + _S555 * _S638;
        float _S640 = _S636 + _S636 + _S553 * _S638;
        float _S641 = _S552 * - _S633 + _S554 * - _S635 + _S556 * - _S637;
        FixedArray<float, 21>  _S642;
        _S642[int(0)] = 0.0f;
        _S642[int(1)] = 0.0f;
        _S642[int(2)] = 0.0f;
        _S642[int(3)] = 0.0f;
        _S642[int(4)] = 0.0f;
        _S642[int(5)] = 0.0f;
        _S642[int(6)] = 0.0f;
        _S642[int(7)] = 0.0f;
        _S642[int(8)] = 0.0f;
        _S642[int(9)] = 0.0f;
        _S642[int(10)] = 0.0f;
        _S642[int(11)] = 0.0f;
        _S642[int(12)] = 0.0f;
        _S642[int(13)] = 0.0f;
        _S642[int(14)] = 0.0f;
        _S642[int(15)] = 0.0f;
        _S642[int(16)] = 0.0f;
        _S642[int(17)] = 0.0f;
        _S642[int(18)] = 0.0f;
        _S642[int(19)] = 0.0f;
        _S642[int(20)] = 0.0f;
        _S642[int(5)] = _S631;
        _S642[int(4)] = _S632;
        _S642[int(3)] = _S639;
        _S642[int(2)] = _S640;
        _S642[int(6)] = _S627;
        float _S643 = _S603 + _S642[int(0)];
        float _S644 = _S604 + _S642[int(1)];
        float _S645 = _S605 + _S642[int(2)];
        float _S646 = _S606 + _S642[int(3)];
        float _S647 = _S607 + _S642[int(4)];
        float _S648 = _S608 + _S642[int(5)];
        float _S649 = _S609 + _S642[int(6)];
        float _S650 = _S610 + _S642[int(7)];
        float _S651 = _S611 + _S642[int(8)];
        float _S652 = _S612 + _S642[int(9)];
        float _S653 = _S613 + _S642[int(10)];
        float _S654 = _S614 + _S642[int(11)];
        float _S655 = _S615 + _S642[int(12)];
        float _S656 = _S616 + _S642[int(13)];
        float _S657 = _S617 + _S642[int(14)];
        float _S658 = _S618 + _S642[int(15)];
        float _S659 = _S619 + _S642[int(16)];
        float _S660 = _S620 + _S642[int(17)];
        float _S661 = _S621 + _S642[int(18)];
        float _S662 = _S622 + _S642[int(19)];
        float _S663 = _S623 + _S642[int(20)];
        _S543 = _S641;
        _S624[int(0)] = _S643;
        _S624[int(1)] = _S644;
        _S624[int(2)] = _S645;
        _S624[int(3)] = _S646;
        _S624[int(4)] = _S647;
        _S624[int(5)] = _S648;
        _S624[int(6)] = _S649;
        _S624[int(7)] = _S650;
        _S624[int(8)] = _S651;
        _S624[int(9)] = _S652;
        _S624[int(10)] = _S653;
        _S624[int(11)] = _S654;
        _S624[int(12)] = _S655;
        _S624[int(13)] = _S656;
        _S624[int(14)] = _S657;
        _S624[int(15)] = _S658;
        _S624[int(16)] = _S659;
        _S624[int(17)] = _S660;
        _S624[int(18)] = _S661;
        _S624[int(19)] = _S662;
        _S624[int(20)] = _S663;
    }
    else
    {
        _S543 = 0.0f;
        _S624[int(0)] = _S603;
        _S624[int(1)] = _S604;
        _S624[int(2)] = _S605;
        _S624[int(3)] = _S606;
        _S624[int(4)] = _S607;
        _S624[int(5)] = _S608;
        _S624[int(6)] = _S609;
        _S624[int(7)] = _S610;
        _S624[int(8)] = _S611;
        _S624[int(9)] = _S612;
        _S624[int(10)] = _S613;
        _S624[int(11)] = _S614;
        _S624[int(12)] = _S615;
        _S624[int(13)] = _S616;
        _S624[int(14)] = _S617;
        _S624[int(15)] = _S618;
        _S624[int(16)] = _S619;
        _S624[int(17)] = _S620;
        _S624[int(18)] = _S621;
        _S624[int(19)] = _S622;
        _S624[int(20)] = _S623;
    }
    float _S664 = -10.0f * (*_s_dOut_3)[int(1)];
    DiffPair_float_0 _S665;
    (&_S665)->primal_0 = _S541;
    (&_S665)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S665, _S664);
    float _S666 = _S665.differential_0 / _S540;
    float _S667 = _S539 * _S666;
    float _S668 = (*_s_dOut_3)[int(0)] / _S540;
    float _S669 = _S539 * _S668;
    float _S670 = _S537[int(1)] * - _S666 + _S538 * - _S668;
    DiffPair_float_0 _S671;
    (&_S671)->primal_0 = _S537[int(16)];
    (&_S671)->differential_0 = 0.0f;
    DiffPair_float_0 _S672;
    (&_S672)->primal_0 = 1.0f;
    (&_S672)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S671, &_S672, _S670);
    float _S673 = (*weights_1)[int(0)] * _S669;
    FixedArray<float, 21>  _S674;
    _S674[int(0)] = 0.0f;
    _S674[int(1)] = 0.0f;
    _S674[int(2)] = 0.0f;
    _S674[int(3)] = 0.0f;
    _S674[int(4)] = 0.0f;
    _S674[int(5)] = 0.0f;
    _S674[int(6)] = 0.0f;
    _S674[int(7)] = 0.0f;
    _S674[int(8)] = 0.0f;
    _S674[int(9)] = 0.0f;
    _S674[int(10)] = 0.0f;
    _S674[int(11)] = 0.0f;
    _S674[int(12)] = 0.0f;
    _S674[int(13)] = 0.0f;
    _S674[int(14)] = 0.0f;
    _S674[int(15)] = 0.0f;
    _S674[int(16)] = 0.0f;
    _S674[int(17)] = 0.0f;
    _S674[int(18)] = 0.0f;
    _S674[int(19)] = 0.0f;
    _S674[int(20)] = 0.0f;
    _S674[int(17)] = _S543;
    _S674[int(1)] = _S667;
    _S674[int(16)] = _S671.differential_0;
    _S674[int(0)] = _S673;
    FixedArray<float, 21>  _S675 = {
        _S624[int(0)] + _S674[int(0)], _S624[int(1)] + _S674[int(1)], _S624[int(2)] + _S674[int(2)], _S624[int(3)] + _S674[int(3)], _S624[int(4)] + _S674[int(4)], _S624[int(5)] + _S674[int(5)], _S624[int(6)] + _S674[int(6)], _S624[int(7)] + _S674[int(7)], _S624[int(8)] + _S674[int(8)], _S624[int(9)] + _S674[int(9)], _S624[int(10)] + _S674[int(10)], _S624[int(11)] + _S674[int(11)], _S624[int(12)] + _S674[int(12)], _S624[int(13)] + _S674[int(13)], _S624[int(14)] + _S674[int(14)], _S624[int(15)] + _S674[int(15)], _S624[int(16)] + _S674[int(16)], _S624[int(17)] + _S674[int(17)], _S624[int(18)] + _S674[int(18)], _S624[int(19)] + _S674[int(19)], _S624[int(20)] + _S674[int(20)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S675;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C21x3E_0 * _S676, float _S677, FixedArray<float, 10>  * _S678, FixedArray<float, 10>  * _S679)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S676, _S677, _S678, _S679);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 21>  * raw_losses_1, float num_pixels_2, FixedArray<float, 10>  * weights_2, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 21>  * _S680)
{
    FixedArray<float, 21>  _S681 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C21x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S681;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, num_pixels_2, weights_2, v_losses_1);
    *_S680 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_15, float3  y_8)
{
    float3  result_9;
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
        *_slang_vector_get_element_ptr(&result_9, i_6) = (F32_min((_slang_vector_get_element(x_15, i_6)), (_slang_vector_get_element(y_8, i_6))));
        i_6 = i_6 + int(1);
    }
    return result_9;
}

inline __device__ float3  max_0(float3  x_16, float3  y_9)
{
    float3  result_10;
    int i_7 = int(0);
    for(;;)
    {
        if(i_7 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_10, i_7) = (F32_max((_slang_vector_get_element(x_16, i_7)), (_slang_vector_get_element(y_9, i_7))));
        i_7 = i_7 + int(1);
    }
    return result_10;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_14, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_13)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_14).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_4).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_13.x);
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_14).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_4).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_13.y);
    *&((&left_d_result_0)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_0)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_14).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_4).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_13.z);
    *&((&left_d_result_0)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_0)->z) = right_dp_2.differential_0;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = left_d_result_0;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  clamp_1(float3  x_17, float3  minBound_1, float3  maxBound_1)
{
    return min_0(max_0(x_17, minBound_1), maxBound_1);
}

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_1(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S682, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S683, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S684, float3  _S685)
{
    _d_clamp_vector_0(_S682, _S683, _S684, _S685);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_4)
{
    float _S686 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S687 = make_float3 (_S686);
    float3  _S688 = make_float3 (0.0f);
    float3  _S689 = make_float3 (1.0f);
    float3  _S690 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S691;
    (&_S691)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S686) * (*dpbackground_0).primal_0;
    (&_S691)->differential_0 = _S690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S692;
    (&_S692)->primal_0 = _S688;
    (&_S692)->differential_0 = _S690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S693;
    (&_S693)->primal_0 = _S689;
    (&_S693)->differential_0 = _S690;
    s_bwd_prop_clamp_1(&_S691, &_S692, &_S693, _s_dOut_4);
    float3  _S694 = _S687 * _S691.differential_0;
    float3  _S695 = (*dpbackground_0).primal_0 * _S691.differential_0;
    float _S696 = - (_S695.x + _S695.y + _S695.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S694;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S696;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S691.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S697, DiffPair_float_0 * _S698, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S699, float3  _S700)
{
    s_bwd_prop_blend_background_0(_S697, _S698, _S699, _S700);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S701 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S701;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S701;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_18)
{
    Matrix<float, 3, 3>  result_11;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(3))
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
            *_slang_vector_get_element_ptr(((&result_11)->rows + (r_0)), c_0) = _slang_vector_get_element(x_18.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_11;
}

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_2)
{
    float x_19 = quat_2.y;
    float x2_0 = x_19 * x_19;
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

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_0, float3  right_0)
{
    float3  result_12;
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
        int j_0 = int(0);
        float sum_0 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_1 = sum_0 + _slang_vector_get_element(left_0.rows[i_8], j_0) * _slang_vector_get_element(right_0, j_0);
            j_0 = j_0 + int(1);
            sum_0 = sum_1;
        }
        *_slang_vector_get_element_ptr(&result_12, i_8) = sum_0;
        i_8 = i_8 + int(1);
    }
    return result_12;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
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
            int i_9 = int(0);
            float sum_2 = 0.0f;
            for(;;)
            {
                if(i_9 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_3 = sum_2 + _slang_vector_get_element(left_1.rows[r_1], i_9) * _slang_vector_get_element(right_1.rows[i_9], c_1);
                i_9 = i_9 + int(1);
                sum_2 = sum_3;
            }
            *_slang_vector_get_element_ptr(((&result_13)->rows + (r_1)), c_1) = sum_2;
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_13;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_1(mul_1(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_3, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_20 = quat_3.y;
    float x2_1 = x_20 * x_20;
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

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_4, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_21 = quat_4.y;
    float x2_2 = x_21 * x_21;
    float y2_2 = quat_4.z * quat_4.z;
    float z2_2 = quat_4.w * quat_4.w;
    float xy_2 = quat_4.y * quat_4.z;
    float xz_2 = quat_4.y * quat_4.w;
    float yz_2 = quat_4.z * quat_4.w;
    float wx_2 = quat_4.x * quat_4.y;
    float wy_2 = quat_4.x * quat_4.z;
    float wz_2 = quat_4.x * quat_4.w;
    *M_1 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

struct CameraDistortion_0
{
    float4  radial_coeffs_0;
    float2  tangential_coeffs_0;
    float2  thin_prism_coeffs_0;
};

inline __device__ CameraDistortion_0 CameraDistortion_x24init_0(float4  radial_coeffs_1, float2  tangential_coeffs_1, float2  thin_prism_coeffs_1)
{
    CameraDistortion_0 _S702;
    (&_S702)->radial_coeffs_0 = radial_coeffs_1;
    (&_S702)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S702)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    return _S702;
}

inline __device__ float2  distort_point(float2  uv_0, bool is_fisheye_0, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2)
{
    float2  _S703;
    if(is_fisheye_0)
    {
        float r_2 = length_1(uv_0);
        float theta_0 = (F32_atan((r_2)));
        float _S704;
        if(r_2 < 0.00100000004749745f)
        {
            _S704 = 1.0f - r_2 * r_2 / 3.0f;
        }
        else
        {
            _S704 = theta_0 / r_2;
        }
        _S703 = uv_0 * make_float2 (_S704);
    }
    else
    {
        _S703 = uv_0;
    }
    CameraDistortion_0 _S705 = CameraDistortion_x24init_0(radial_coeffs_2, tangential_coeffs_2, thin_prism_coeffs_2);
    float p1_0 = _S705.tangential_coeffs_0.x;
    float p2_0 = _S705.tangential_coeffs_0.y;
    float u_0 = _S703.x;
    float v_0 = _S703.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    return _S703 * make_float2 (1.0f + r2_0 * (_S705.radial_coeffs_0.x + r2_0 * (_S705.radial_coeffs_0.y + r2_0 * (_S705.radial_coeffs_0.z + r2_0 * _S705.radial_coeffs_0.w)))) + make_float2 (2.0f * p1_0 * u_0 * v_0 + p2_0 * (r2_0 + 2.0f * u_0 * u_0) + _S705.thin_prism_coeffs_0.x * r2_0, 2.0f * p2_0 * u_0 * v_0 + p1_0 * (r2_0 + 2.0f * v_0 * v_0) + _S705.thin_prism_coeffs_0.y * r2_0);
}

inline __device__ float2  undistort_point_0(float2  uv_1, CameraDistortion_0 * dist_coeffs_0, int maxiter_0)
{
    int i_10 = int(0);
    float2  q_0 = uv_1;
    for(;;)
    {
        if(i_10 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float k4_0 = dist_coeffs_0->radial_coeffs_0.w;
        float p1_1 = dist_coeffs_0->tangential_coeffs_0.x;
        float p2_1 = dist_coeffs_0->tangential_coeffs_0.y;
        float sx1_0 = dist_coeffs_0->thin_prism_coeffs_0.x;
        float sy1_0 = dist_coeffs_0->thin_prism_coeffs_0.y;
        float u_1 = q_0.x;
        float v_1 = q_0.y;
        float r2_1 = u_1 * u_1 + v_1 * v_1;
        float _S706 = dist_coeffs_0->radial_coeffs_0.z + r2_1 * k4_0;
        float _S707 = dist_coeffs_0->radial_coeffs_0.y + r2_1 * _S706;
        float _S708 = dist_coeffs_0->radial_coeffs_0.x + r2_1 * _S707;
        float radial_0 = 1.0f + r2_1 * _S708;
        float _S709 = 2.0f * p1_1;
        float _S710 = _S709 * u_1;
        float _S711 = 2.0f * u_1;
        float _S712 = 2.0f * p2_1;
        float _S713 = _S712 * u_1;
        float _S714 = 2.0f * v_1;
        float2  _S715 = q_0 * make_float2 (radial_0) + make_float2 (_S710 * v_1 + p2_1 * (r2_1 + _S711 * u_1) + sx1_0 * r2_1, _S713 * v_1 + p1_1 * (r2_1 + _S714 * v_1) + sy1_0 * r2_1);
        float2  _S716 = make_float2 (0.0f);
        float2  seed_0 = _S716;
        *&((&seed_0)->x) = 1.0f;
        float2  _S717 = make_float2 (radial_0);
        float2  _S718 = q_0 * seed_0;
        float _S719 = p1_1 * seed_0.y;
        float _S720 = p2_1 * seed_0.x;
        float _S721 = _S718.x + _S718.y;
        float _S722 = r2_1 * _S721;
        float _S723 = r2_1 * _S722;
        float _S724 = sy1_0 * seed_0.y + _S719 + sx1_0 * seed_0.x + _S720 + _S708 * _S721 + _S707 * _S722 + _S706 * _S723 + k4_0 * (r2_1 * _S723);
        float _S725 = v_1 * _S724;
        float _S726 = u_1 * _S724;
        Matrix<float, 2, 2>  J_0;
        J_0[int(0)] = _S717 * seed_0 + make_float2 (_S712 * (v_1 * seed_0.y) + _S711 * _S720 + 2.0f * (u_1 * _S720) + _S709 * (v_1 * seed_0.x) + _S726 + _S726, _S714 * _S719 + 2.0f * (v_1 * _S719) + _S713 * seed_0.y + _S710 * seed_0.x + _S725 + _S725);
        float2  seed_1 = _S716;
        *&((&seed_1)->y) = 1.0f;
        float2  _S727 = q_0 * seed_1;
        float _S728 = p1_1 * seed_1.y;
        float _S729 = p2_1 * seed_1.x;
        float _S730 = _S727.x + _S727.y;
        float _S731 = r2_1 * _S730;
        float _S732 = r2_1 * _S731;
        float _S733 = sy1_0 * seed_1.y + _S728 + sx1_0 * seed_1.x + _S729 + _S708 * _S730 + _S707 * _S731 + _S706 * _S732 + k4_0 * (r2_1 * _S732);
        float _S734 = v_1 * _S733;
        float _S735 = u_1 * _S733;
        J_0[int(1)] = _S717 * seed_1 + make_float2 (_S712 * (v_1 * seed_1.y) + _S711 * _S729 + 2.0f * (u_1 * _S729) + _S709 * (v_1 * seed_1.x) + _S735 + _S735, _S714 * _S728 + 2.0f * (v_1 * _S728) + _S713 * seed_1.y + _S710 * seed_1.x + _S734 + _S734);
        float2  _S736 = _S715 - uv_1;
        float inv_det_0 = 1.0f / (J_0.rows[int(0)].x * J_0.rows[int(1)].y - J_0.rows[int(0)].y * J_0.rows[int(1)].x);
        float _S737 = _S736.x;
        float _S738 = _S736.y;
        float2  q_1 = q_0 - make_float2 ((_S737 * J_0.rows[int(1)].y - _S738 * J_0.rows[int(0)].y) * inv_det_0, (- _S737 * J_0.rows[int(1)].x + _S738 * J_0.rows[int(0)].x) * inv_det_0);
        i_10 = i_10 + int(1);
        q_0 = q_1;
    }
    return q_0;
}

inline __device__ float2  undistort_point(float2  uv_2, bool is_fisheye_1, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3)
{
    CameraDistortion_0 _S739;
    (&_S739)->radial_coeffs_0 = radial_coeffs_3;
    (&_S739)->tangential_coeffs_0 = tangential_coeffs_3;
    (&_S739)->thin_prism_coeffs_0 = thin_prism_coeffs_3;
    float2  _S740 = undistort_point_0(uv_2, &_S739, int(8));
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float theta_1 = length_1(_S740);
        float _S741;
        if(theta_1 < 0.00100000004749745f)
        {
            _S741 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S741 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S742 = make_float3 ((_S740 * make_float2 (_S741)).x, (_S740 * make_float2 (_S741)).y, (F32_cos((theta_1))));
        raydir_0 = _S742;
    }
    else
    {
        raydir_0 = make_float3 (_S740.x, _S740.y, 1.0f);
    }
    return float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_14)
{
    float _S743 = (*right_2).primal_0.rows[int(0)].x * dOut_14.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_14.x;
    float sum_4 = _S743 + (*right_2).primal_0.rows[int(0)].y * dOut_14.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_14.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(0)].z * dOut_14.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_14.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_5;
    float _S744 = (*right_2).primal_0.rows[int(1)].x * dOut_14.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_14.x;
    float sum_6 = _S744 + (*right_2).primal_0.rows[int(1)].y * dOut_14.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_14.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(1)].z * dOut_14.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_14.z;
    *&((&left_d_result_1)->y) = sum_7;
    float _S745 = (*right_2).primal_0.rows[int(2)].x * dOut_14.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_14.x;
    float sum_8 = _S745 + (*right_2).primal_0.rows[int(2)].y * dOut_14.y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_14.y;
    float sum_9 = sum_8 + (*right_2).primal_0.rows[int(2)].z * dOut_14.z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_14.z;
    *&((&left_d_result_1)->z) = sum_9;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_2(float3  left_3, Matrix<float, 3, 3>  right_3)
{
    float3  result_14;
    int j_1 = int(0);
    for(;;)
    {
        if(j_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int i_11 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(i_11 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_3, i_11) * _slang_vector_get_element(right_3.rows[i_11], j_1);
            i_11 = i_11 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_14, j_1) = sum_10;
        j_1 = j_1 + int(1);
    }
    return result_14;
}

inline __device__ float3  normalize_0(float3  x_22)
{
    return x_22 / make_float3 (length_2(x_22));
}

inline __device__ void generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_3, bool is_fisheye_2, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, float3  * ray_o_0, float3  * ray_d_0)
{
    *ray_o_0 = - mul_2(t_1, R_2);
    CameraDistortion_0 _S746;
    (&_S746)->radial_coeffs_0 = radial_coeffs_4;
    (&_S746)->tangential_coeffs_0 = tangential_coeffs_4;
    (&_S746)->thin_prism_coeffs_0 = thin_prism_coeffs_4;
    float2  _S747 = undistort_point_0(uv_3, &_S746, int(8));
    float3  raydir_1;
    if(is_fisheye_2)
    {
        float theta_2 = length_1(_S747);
        float _S748;
        if(theta_2 < 0.00100000004749745f)
        {
            _S748 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S748 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S749 = make_float3 ((_S747 * make_float2 (_S748)).x, (_S747 * make_float2 (_S748)).y, (F32_cos((theta_2))));
        raydir_1 = _S749;
    }
    else
    {
        raydir_1 = make_float3 (_S747.x, _S747.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_2(raydir_1, R_2));
    return;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S750;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S751, Matrix<float, 3, 3>  _S752)
{
    return mul_2(_S751, _S752);
}

inline __device__ float s_primal_ctx_sin_0(float _S753)
{
    return (F32_sin((_S753)));
}

inline __device__ float s_primal_ctx_cos_0(float _S754)
{
    return (F32_cos((_S754)));
}

inline __device__ void s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_4, bool is_fisheye_3, float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S750 = make_float2 (0.0f);
    float3  _S755 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    CameraDistortion_0 _S756;
    (&_S756)->radial_coeffs_0 = radial_coeffs_5;
    (&_S756)->tangential_coeffs_0 = tangential_coeffs_5;
    (&_S756)->thin_prism_coeffs_0 = thin_prism_coeffs_5;
    float2  _S757 = undistort_point_0(uv_4, &_S756, int(8));
    _s_diff_ctx_0->_S750 = _S757;
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float _S758 = length_1(_S757);
        float _S759;
        if(_S758 < 0.00100000004749745f)
        {
            _S759 = 1.0f - _S758 * _S758 / 6.0f;
        }
        else
        {
            _S759 = s_primal_ctx_sin_0(_S758) / _S758;
        }
        float3  _S760 = make_float3 ((_S757 * make_float2 (_S759)).x, (_S757 * make_float2 (_S759)).y, s_primal_ctx_cos_0(_S758));
        raydir_2 = _S760;
    }
    else
    {
        raydir_2 = make_float3 (_S757.x, _S757.y, 1.0f);
    }
    float3  _S761 = normalize_0(s_primal_ctx_mul_0(raydir_2, dpR_0));
    *dpray_o_0 = _S755;
    *dpray_d_0 = _S761;
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float _s_dOut_5)
{
    float _S762 = (*dpx_15).primal_0.x;
    float _S763 = (*dpx_15).primal_0.y;
    float _S764 = (*dpx_15).primal_0.z;
    DiffPair_float_0 _S765;
    (&_S765)->primal_0 = _S762 * _S762 + _S763 * _S763 + _S764 * _S764;
    (&_S765)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S765, _s_dOut_5);
    float _S766 = (*dpx_15).primal_0.z * _S765.differential_0;
    float _S767 = _S766 + _S766;
    float _S768 = (*dpx_15).primal_0.y * _S765.differential_0;
    float _S769 = _S768 + _S768;
    float _S770 = (*dpx_15).primal_0.x * _S765.differential_0;
    float _S771 = _S770 + _S770;
    float3  _S772 = make_float3 (0.0f);
    *&((&_S772)->z) = _S767;
    *&((&_S772)->y) = _S769;
    *&((&_S772)->x) = _S771;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S772;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S773, float _S774)
{
    s_bwd_prop_length_impl_1(_S773, _S774);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  _s_dOut_6)
{
    float _S775 = length_2((*dpx_16).primal_0);
    float3  _S776 = (*dpx_16).primal_0 * _s_dOut_6;
    float3  _S777 = make_float3 (1.0f / _S775) * _s_dOut_6;
    float _S778 = - ((_S776.x + _S776.y + _S776.z) / (_S775 * _S775));
    float3  _S779 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S780;
    (&_S780)->primal_0 = (*dpx_16).primal_0;
    (&_S780)->differential_0 = _S779;
    s_bwd_length_impl_1(&_S780, _S778);
    float3  _S781 = _S777 + _S780.differential_0;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S781;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S782, float3  _S783)
{
    s_bwd_prop_normalize_impl_0(_S782, _S783);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S784, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S785, float3  _S786)
{
    _d_mul_0(_S784, _S785, _S786);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_5, bool is_fisheye_4, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S787 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S788 = *dpt_1;
    float3  raydir_3;
    if(is_fisheye_4)
    {
        float _S789 = length_1(_s_diff_ctx_1->_S750);
        float _S790;
        if(_S789 < 0.00100000004749745f)
        {
            _S790 = 1.0f - _S789 * _S789 / 6.0f;
        }
        else
        {
            _S790 = s_primal_ctx_sin_0(_S789) / _S789;
        }
        float3  _S791 = make_float3 ((_s_diff_ctx_1->_S750 * make_float2 (_S790)).x, (_s_diff_ctx_1->_S750 * make_float2 (_S790)).y, s_primal_ctx_cos_0(_S789));
        raydir_3 = _S791;
    }
    else
    {
        raydir_3 = make_float3 (_s_diff_ctx_1->_S750.x, _s_diff_ctx_1->_S750.y, 1.0f);
    }
    float3  _S792 = s_primal_ctx_mul_0(raydir_3, _S787.primal_0);
    float3  _S793 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S794;
    (&_S794)->primal_0 = _S792;
    (&_S794)->differential_0 = _S793;
    s_bwd_normalize_impl_0(&_S794, dpray_d_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S795;
    (&_S795)->primal_0 = raydir_3;
    (&_S795)->differential_0 = _S793;
    Matrix<float, 3, 3>  _S796 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S797;
    (&_S797)->primal_0 = _S787.primal_0;
    (&_S797)->differential_0 = _S796;
    s_bwd_prop_mul_0(&_S795, &_S797, _S794.differential_0);
    float3  _S798 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S799;
    (&_S799)->primal_0 = _S788.primal_0;
    (&_S799)->differential_0 = _S793;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S800;
    (&_S800)->primal_0 = _S787.primal_0;
    (&_S800)->differential_0 = _S796;
    s_bwd_prop_mul_0(&_S799, &_S800, _S798);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S799.differential_0;
    Matrix<float, 3, 3>  _S801 = _S800.differential_0 + _S797.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S801;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S802, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S803, float2  _S804, bool _S805, float4  _S806, float2  _S807, float2  _S808, float3  _S809, float3  _S810)
{
    float3  _S811;
    float3  _S812;
    s_bwd_prop_generate_ray_Intermediates_0 _S813;
    s_primal_ctx_generate_ray_0((*_S802).primal_0, (*_S803).primal_0, _S804, _S805, _S806, _S807, _S808, &_S811, &_S812, &_S813);
    s_bwd_prop_generate_ray_Intermediates_0 _S814 = _S813;
    s_bwd_prop_generate_ray_0(_S802, _S803, _S804, _S805, _S806, _S807, _S808, _S809, _S810, &_S814);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_6, bool is_fisheye_5, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S815 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S815;
    float3  _S816 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S816;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_6, is_fisheye_5, radial_coeffs_7, tangential_coeffs_7, thin_prism_coeffs_7, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S817 = dOut_15.y;
    float _S818 = dOut_15.z;
    float _S819 = dOut_15.x;
    float _S820 = (*a_0).primal_0.z * _S817 + - (*a_0).primal_0.y * _S818;
    float _S821 = - (*a_0).primal_0.z * _S819 + (*a_0).primal_0.x * _S818;
    float _S822 = (*a_0).primal_0.y * _S819 + - (*a_0).primal_0.x * _S817;
    float3  _S823 = make_float3 (- (*b_0).primal_0.z * _S817 + (*b_0).primal_0.y * _S818, (*b_0).primal_0.z * _S819 + - (*b_0).primal_0.x * _S818, - (*b_0).primal_0.y * _S819 + (*b_0).primal_0.x * _S817);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S823;
    float3  _S824 = make_float3 (_S820, _S821, _S822);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S824;
    return;
}

inline __device__ float3  cross_0(float3  left_4, float3  right_4)
{
    float _S825 = left_4.y;
    float _S826 = right_4.z;
    float _S827 = left_4.z;
    float _S828 = right_4.y;
    float _S829 = right_4.x;
    float _S830 = left_4.x;
    return make_float3 (_S825 * _S826 - _S827 * _S828, _S827 * _S829 - _S830 * _S826, _S830 * _S828 - _S825 * _S829);
}

inline __device__ void depth_to_normal(uint width_0, uint height_0, float2  pix_center_0, float4  intrins_0, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, bool is_fisheye_6, bool is_ray_depth_0, float4  depths_0, float3  * normal_0)
{
    FixedArray<float3 , 4>  points_0;
    float2  _S831 = float2 {intrins_0.z, intrins_0.w};
    float2  _S832 = float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S831) / _S832;
    float3  raydir_4;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S833;
        (&_S833)->radial_coeffs_0 = radial_coeffs_8;
        (&_S833)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S833)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S834 = undistort_point_0(uv_7, &_S833, int(12));
        float theta_3 = length_1(_S834);
        float3  raydir_5 = make_float3 ((_S834 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).x, (_S834 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).y, (F32_cos((theta_3))));
        if(!is_ray_depth_0)
        {
            raydir_4 = raydir_5 / make_float3 (raydir_5.z);
        }
        else
        {
            raydir_4 = raydir_5;
        }
    }
    else
    {
        float3  raydir_6 = make_float3 (uv_7.x, uv_7.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_4 = normalize_0(raydir_6);
        }
        else
        {
            raydir_4 = raydir_6;
        }
    }
    points_0[int(0)] = make_float3 (depths_0.x) * raydir_4;
    float2  uv_8 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S831) / _S832;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S835;
        (&_S835)->radial_coeffs_0 = radial_coeffs_8;
        (&_S835)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S835)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S836 = undistort_point_0(uv_8, &_S835, int(12));
        float theta_4 = length_1(_S836);
        float3  raydir_7 = make_float3 ((_S836 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (_S836 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_0)
        {
            raydir_4 = raydir_7 / make_float3 (raydir_7.z);
        }
        else
        {
            raydir_4 = raydir_7;
        }
    }
    else
    {
        float3  raydir_8 = make_float3 (uv_8.x, uv_8.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_4 = normalize_0(raydir_8);
        }
        else
        {
            raydir_4 = raydir_8;
        }
    }
    points_0[int(1)] = make_float3 (depths_0.y) * raydir_4;
    float2  uv_9 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S831) / _S832;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S837;
        (&_S837)->radial_coeffs_0 = radial_coeffs_8;
        (&_S837)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S837)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S838 = undistort_point_0(uv_9, &_S837, int(12));
        float theta_5 = length_1(_S838);
        float3  raydir_9 = make_float3 ((_S838 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (_S838 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        if(!is_ray_depth_0)
        {
            raydir_4 = raydir_9 / make_float3 (raydir_9.z);
        }
        else
        {
            raydir_4 = raydir_9;
        }
    }
    else
    {
        float3  raydir_10 = make_float3 (uv_9.x, uv_9.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_4 = normalize_0(raydir_10);
        }
        else
        {
            raydir_4 = raydir_10;
        }
    }
    points_0[int(2)] = make_float3 (depths_0.z) * raydir_4;
    float2  uv_10 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S831) / _S832;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S839;
        (&_S839)->radial_coeffs_0 = radial_coeffs_8;
        (&_S839)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S839)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S840 = undistort_point_0(uv_10, &_S839, int(12));
        float theta_6 = length_1(_S840);
        float3  raydir_11 = make_float3 ((_S840 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (_S840 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
        if(!is_ray_depth_0)
        {
            raydir_4 = raydir_11 / make_float3 (raydir_11.z);
        }
        else
        {
            raydir_4 = raydir_11;
        }
    }
    else
    {
        float3  raydir_12 = make_float3 (uv_10.x, uv_10.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_4 = normalize_0(raydir_12);
        }
        else
        {
            raydir_4 = raydir_12;
        }
    }
    points_0[int(3)] = make_float3 (depths_0.w) * raydir_4;
    float3  _S841 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    *normal_0 = _S841;
    if((length_2(_S841)) != 0.0f)
    {
        *normal_0 = *normal_0 / make_float3 (length_2(*normal_0));
    }
    return;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S842;
    float2  _S843;
    float2  _S844;
    float2  _S845;
};

inline __device__ float3  s_primal_ctx_cross_0(float3  _S846, float3  _S847)
{
    return cross_0(_S846, _S847);
}

inline __device__ void s_primal_ctx_depth_to_normal_0(uint width_1, uint height_1, float2  pix_center_1, float4  intrins_1, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, bool is_fisheye_7, bool is_ray_depth_1, float4  dpdepths_0, float3  * dpnormal_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S848 = make_float2 (0.0f);
    _s_diff_ctx_2->_S842 = _S848;
    _s_diff_ctx_2->_S843 = _S848;
    _s_diff_ctx_2->_S844 = _S848;
    _s_diff_ctx_2->_S845 = _S848;
    _s_diff_ctx_2->_S842 = _S848;
    _s_diff_ctx_2->_S843 = _S848;
    _s_diff_ctx_2->_S844 = _S848;
    _s_diff_ctx_2->_S845 = _S848;
    float2  _S849 = float2 {intrins_1.z, intrins_1.w};
    float2  _S850 = float2 {intrins_1.x, intrins_1.y};
    float2  uv_11 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S849) / _S850;
    float3  raydir_13;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S851;
        (&_S851)->radial_coeffs_0 = radial_coeffs_9;
        (&_S851)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S851)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S852 = undistort_point_0(uv_11, &_S851, int(12));
        _s_diff_ctx_2->_S842 = _S852;
        float _S853 = length_1(_S852);
        float3  raydir_14 = make_float3 ((_S852 / make_float2 (s_primal_ctx_max_0(_S853, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S853))).x, (_S852 / make_float2 (s_primal_ctx_max_0(_S853, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S853))).y, s_primal_ctx_cos_0(_S853));
        if(!is_ray_depth_1)
        {
            raydir_13 = raydir_14 / make_float3 (raydir_14.z);
        }
        else
        {
            raydir_13 = raydir_14;
        }
    }
    else
    {
        float3  raydir_15 = make_float3 (uv_11.x, uv_11.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_13 = normalize_0(raydir_15);
        }
        else
        {
            raydir_13 = raydir_15;
        }
    }
    float3  _S854 = make_float3 (dpdepths_0.x) * raydir_13;
    float2  uv_12 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S849) / _S850;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S855;
        (&_S855)->radial_coeffs_0 = radial_coeffs_9;
        (&_S855)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S855)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S856 = undistort_point_0(uv_12, &_S855, int(12));
        _s_diff_ctx_2->_S843 = _S856;
        float _S857 = length_1(_S856);
        float3  raydir_16 = make_float3 ((_S856 / make_float2 (s_primal_ctx_max_0(_S857, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S857))).x, (_S856 / make_float2 (s_primal_ctx_max_0(_S857, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S857))).y, s_primal_ctx_cos_0(_S857));
        if(!is_ray_depth_1)
        {
            raydir_13 = raydir_16 / make_float3 (raydir_16.z);
        }
        else
        {
            raydir_13 = raydir_16;
        }
    }
    else
    {
        float3  raydir_17 = make_float3 (uv_12.x, uv_12.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_13 = normalize_0(raydir_17);
        }
        else
        {
            raydir_13 = raydir_17;
        }
    }
    float3  _S858 = make_float3 (dpdepths_0.y) * raydir_13;
    float2  uv_13 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S849) / _S850;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S859;
        (&_S859)->radial_coeffs_0 = radial_coeffs_9;
        (&_S859)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S859)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S860 = undistort_point_0(uv_13, &_S859, int(12));
        _s_diff_ctx_2->_S844 = _S860;
        float _S861 = length_1(_S860);
        float3  raydir_18 = make_float3 ((_S860 / make_float2 (s_primal_ctx_max_0(_S861, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S861))).x, (_S860 / make_float2 (s_primal_ctx_max_0(_S861, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S861))).y, s_primal_ctx_cos_0(_S861));
        if(!is_ray_depth_1)
        {
            raydir_13 = raydir_18 / make_float3 (raydir_18.z);
        }
        else
        {
            raydir_13 = raydir_18;
        }
    }
    else
    {
        float3  raydir_19 = make_float3 (uv_13.x, uv_13.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_13 = normalize_0(raydir_19);
        }
        else
        {
            raydir_13 = raydir_19;
        }
    }
    float3  _S862 = make_float3 (dpdepths_0.z) * raydir_13;
    float2  uv_14 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S849) / _S850;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S863;
        (&_S863)->radial_coeffs_0 = radial_coeffs_9;
        (&_S863)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S863)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S864 = undistort_point_0(uv_14, &_S863, int(12));
        _s_diff_ctx_2->_S845 = _S864;
        float _S865 = length_1(_S864);
        float3  raydir_20 = make_float3 ((_S864 / make_float2 (s_primal_ctx_max_0(_S865, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S865))).x, (_S864 / make_float2 (s_primal_ctx_max_0(_S865, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S865))).y, s_primal_ctx_cos_0(_S865));
        if(!is_ray_depth_1)
        {
            raydir_13 = raydir_20 / make_float3 (raydir_20.z);
        }
        else
        {
            raydir_13 = raydir_20;
        }
    }
    else
    {
        float3  raydir_21 = make_float3 (uv_14.x, uv_14.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_13 = normalize_0(raydir_21);
        }
        else
        {
            raydir_13 = raydir_21;
        }
    }
    float3  _S866 = make_float3 (dpdepths_0.w) * raydir_13;
    float3  _S867 = s_primal_ctx_cross_0(_S858 - _S854, - (_S866 - _S862));
    float _S868 = length_2(_S867);
    if(_S868 != 0.0f)
    {
        raydir_13 = _S867 / make_float3 (_S868);
    }
    else
    {
        raydir_13 = _S867;
    }
    *dpnormal_0 = raydir_13;
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S869, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S870, float3  _S871)
{
    _d_cross_0(_S869, _S870, _S871);
    return;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(uint width_2, uint height_2, float2  pix_center_2, float4  intrins_2, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, bool is_fisheye_8, bool is_ray_depth_2, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  dpnormal_1, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S872 = *dpdepths_1;
    float3  _S873 = make_float3 (0.0f);
    float2  _S874 = float2 {intrins_2.z, intrins_2.w};
    float2  _S875 = float2 {intrins_2.x, intrins_2.y};
    float2  uv_15 = (pix_center_2 + make_float2 (-1.0f, -0.0f) - _S874) / _S875;
    float3  raydir_22;
    if(is_fisheye_8)
    {
        float _S876 = length_1(_s_diff_ctx_3->_S842);
        float3  raydir_23 = make_float3 ((_s_diff_ctx_3->_S842 / make_float2 (s_primal_ctx_max_0(_S876, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S876))).x, (_s_diff_ctx_3->_S842 / make_float2 (s_primal_ctx_max_0(_S876, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S876))).y, s_primal_ctx_cos_0(_S876));
        if(!is_ray_depth_2)
        {
            raydir_22 = raydir_23 / make_float3 (raydir_23.z);
        }
        else
        {
            raydir_22 = raydir_23;
        }
    }
    else
    {
        float3  raydir_24 = make_float3 (uv_15.x, uv_15.y, 1.0f);
        if(is_ray_depth_2)
        {
            raydir_22 = normalize_0(raydir_24);
        }
        else
        {
            raydir_22 = raydir_24;
        }
    }
    float3  _S877 = make_float3 (_S872.primal_0.x) * raydir_22;
    float2  uv_16 = (pix_center_2 + make_float2 (1.0f, -0.0f) - _S874) / _S875;
    float3  raydir_25;
    if(is_fisheye_8)
    {
        float _S878 = length_1(_s_diff_ctx_3->_S843);
        float3  raydir_26 = make_float3 ((_s_diff_ctx_3->_S843 / make_float2 (s_primal_ctx_max_0(_S878, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S878))).x, (_s_diff_ctx_3->_S843 / make_float2 (s_primal_ctx_max_0(_S878, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S878))).y, s_primal_ctx_cos_0(_S878));
        if(!is_ray_depth_2)
        {
            raydir_25 = raydir_26 / make_float3 (raydir_26.z);
        }
        else
        {
            raydir_25 = raydir_26;
        }
    }
    else
    {
        float3  raydir_27 = make_float3 (uv_16.x, uv_16.y, 1.0f);
        if(is_ray_depth_2)
        {
            raydir_25 = normalize_0(raydir_27);
        }
        else
        {
            raydir_25 = raydir_27;
        }
    }
    float3  _S879 = make_float3 (_S872.primal_0.y) * raydir_25;
    float2  uv_17 = (pix_center_2 + make_float2 (0.0f, -1.0f) - _S874) / _S875;
    float3  raydir_28;
    if(is_fisheye_8)
    {
        float _S880 = length_1(_s_diff_ctx_3->_S844);
        float3  raydir_29 = make_float3 ((_s_diff_ctx_3->_S844 / make_float2 (s_primal_ctx_max_0(_S880, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S880))).x, (_s_diff_ctx_3->_S844 / make_float2 (s_primal_ctx_max_0(_S880, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S880))).y, s_primal_ctx_cos_0(_S880));
        if(!is_ray_depth_2)
        {
            raydir_28 = raydir_29 / make_float3 (raydir_29.z);
        }
        else
        {
            raydir_28 = raydir_29;
        }
    }
    else
    {
        float3  raydir_30 = make_float3 (uv_17.x, uv_17.y, 1.0f);
        if(is_ray_depth_2)
        {
            raydir_28 = normalize_0(raydir_30);
        }
        else
        {
            raydir_28 = raydir_30;
        }
    }
    float3  _S881 = make_float3 (_S872.primal_0.z) * raydir_28;
    float2  uv_18 = (pix_center_2 + make_float2 (0.0f, 1.0f) - _S874) / _S875;
    float3  raydir_31;
    if(is_fisheye_8)
    {
        float _S882 = length_1(_s_diff_ctx_3->_S845);
        float3  raydir_32 = make_float3 ((_s_diff_ctx_3->_S845 / make_float2 (s_primal_ctx_max_0(_S882, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S882))).x, (_s_diff_ctx_3->_S845 / make_float2 (s_primal_ctx_max_0(_S882, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S882))).y, s_primal_ctx_cos_0(_S882));
        if(!is_ray_depth_2)
        {
            raydir_31 = raydir_32 / make_float3 (raydir_32.z);
        }
        else
        {
            raydir_31 = raydir_32;
        }
    }
    else
    {
        float3  raydir_33 = make_float3 (uv_18.x, uv_18.y, 1.0f);
        if(is_ray_depth_2)
        {
            raydir_31 = normalize_0(raydir_33);
        }
        else
        {
            raydir_31 = raydir_33;
        }
    }
    float3  _S883 = make_float3 (_S872.primal_0.w) * raydir_31;
    float3  dx_0 = _S879 - _S877;
    float3  _S884 = - (_S883 - _S881);
    float3  _S885 = s_primal_ctx_cross_0(dx_0, _S884);
    float _S886 = length_2(_S885);
    float3  _S887 = make_float3 (_S886);
    bool _S888 = _S886 != 0.0f;
    float3  _S889;
    if(_S888)
    {
        _S889 = make_float3 (_S886 * _S886);
    }
    else
    {
        _S889 = _S873;
    }
    float3  _S890;
    if(_S888)
    {
        float3  _S891 = dpnormal_1 / _S889;
        float3  _S892 = _S887 * _S891;
        _S889 = _S885 * - _S891;
        _S890 = _S892;
    }
    else
    {
        _S889 = _S873;
        _S890 = dpnormal_1;
    }
    float _S893 = _S889.x + _S889.y + _S889.z;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S894;
    (&_S894)->primal_0 = _S885;
    (&_S894)->differential_0 = _S873;
    s_bwd_length_impl_1(&_S894, _S893);
    float3  _S895 = _S894.differential_0 + _S890;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S896;
    (&_S896)->primal_0 = dx_0;
    (&_S896)->differential_0 = _S873;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S897;
    (&_S897)->primal_0 = _S884;
    (&_S897)->differential_0 = _S873;
    s_bwd_prop_cross_0(&_S896, &_S897, _S895);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S898 = _S896;
    float3  s_diff_dy_T_0 = - _S897.differential_0;
    float3  _S899 = - s_diff_dy_T_0;
    float3  _S900 = - _S896.differential_0;
    float3  _S901 = raydir_31 * s_diff_dy_T_0;
    float _S902 = _S901.x + _S901.y + _S901.z;
    float4  _S903 = make_float4 (0.0f);
    float4  _S904 = _S903;
    *&((&_S904)->w) = _S902;
    float3  _S905 = raydir_28 * _S899;
    float _S906 = _S905.x + _S905.y + _S905.z;
    float4  _S907 = _S903;
    *&((&_S907)->z) = _S906;
    float4  _S908 = _S904 + _S907;
    float3  _S909 = raydir_25 * _S898.differential_0;
    float _S910 = _S909.x + _S909.y + _S909.z;
    float4  _S911 = _S903;
    *&((&_S911)->y) = _S910;
    float4  _S912 = _S908 + _S911;
    float3  _S913 = raydir_22 * _S900;
    float _S914 = _S913.x + _S913.y + _S913.z;
    float4  _S915 = _S903;
    *&((&_S915)->x) = _S914;
    float4  _S916 = _S912 + _S915;
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S916;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(uint _S917, uint _S918, float2  _S919, float4  _S920, float4  _S921, float2  _S922, float2  _S923, bool _S924, bool _S925, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S926, float3  _S927)
{
    float3  _S928;
    s_bwd_prop_depth_to_normal_Intermediates_0 _S929;
    s_primal_ctx_depth_to_normal_0(_S917, _S918, _S919, _S920, _S921, _S922, _S923, _S924, _S925, (*_S926).primal_0, &_S928, &_S929);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S930 = _S929;
    s_bwd_prop_depth_to_normal_0(_S917, _S918, _S919, _S920, _S921, _S922, _S923, _S924, _S925, _S926, _S927, &_S930);
    return;
}

inline __device__ void depth_to_normal_vjp(uint width_3, uint height_3, float2  pix_center_3, float4  intrins_3, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, bool is_fisheye_9, bool is_ray_depth_3, float4  depths_1, float3  v_normal_0, float4  * v_depths_0)
{
    float4  _S931 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S931;
    s_bwd_depth_to_normal_0(width_3, height_3, pix_center_3, intrins_3, radial_coeffs_11, tangential_coeffs_11, thin_prism_coeffs_11, is_fisheye_9, is_ray_depth_3, &dp_depths_0, v_normal_0);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(uint width_4, uint height_4, float2  pix_center_4, float4  intrins_4, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, bool is_fisheye_10)
{
    float2  uv_19 = (pix_center_4 - float2 {intrins_4.z, intrins_4.w}) / float2 {intrins_4.x, intrins_4.y};
    float3  raydir_34;
    if(is_fisheye_10)
    {
        CameraDistortion_0 _S932;
        (&_S932)->radial_coeffs_0 = radial_coeffs_12;
        (&_S932)->tangential_coeffs_0 = tangential_coeffs_12;
        (&_S932)->thin_prism_coeffs_0 = thin_prism_coeffs_12;
        float2  _S933 = undistort_point_0(uv_19, &_S932, int(12));
        float theta_7 = length_1(_S933);
        float3  raydir_35 = make_float3 ((_S933 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (_S933 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        raydir_34 = raydir_35 / make_float3 (raydir_35.z);
    }
    else
    {
        raydir_34 = make_float3 (uv_19.x, uv_19.y, 1.0f);
    }
    return float((F32_sign((raydir_34.z)))) / length_2(raydir_34);
}

