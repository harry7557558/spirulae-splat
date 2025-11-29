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

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_12)
{
    Matrix<float, 3, 3>  result_8;
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
            *_slang_vector_get_element_ptr(((&result_8)->rows + (r_0)), c_0) = _slang_vector_get_element(x_12.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_8;
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

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_0, float3  right_0)
{
    float3  result_9;
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
            float sum_1 = sum_0 + _slang_vector_get_element(left_0.rows[i_5], j_0) * _slang_vector_get_element(right_0, j_0);
            j_0 = j_0 + int(1);
            sum_0 = sum_1;
        }
        *_slang_vector_get_element_ptr(&result_9, i_5) = sum_0;
        i_5 = i_5 + int(1);
    }
    return result_9;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
{
    Matrix<float, 3, 3>  result_10;
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
            int i_6 = int(0);
            float sum_2 = 0.0f;
            for(;;)
            {
                if(i_6 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_3 = sum_2 + _slang_vector_get_element(left_1.rows[r_1], i_6) * _slang_vector_get_element(right_1.rows[i_6], c_1);
                i_6 = i_6 + int(1);
                sum_2 = sum_3;
            }
            *_slang_vector_get_element_ptr(((&result_10)->rows + (r_1)), c_1) = sum_2;
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_10;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_1(mul_1(R_1, covarW_0), transpose_0(R_1));
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
    Matrix<float, 3, 3>  M_0 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_1(M_0, transpose_0(M_0));
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
    *M_1 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
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
    float _S282 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S283 = (*dist_coeffs_0)[int(1)] + r2_0 * _S282;
    float _S284 = (*dist_coeffs_0)[int(0)] + r2_0 * _S283;
    float2  _S285 = make_float2 (1.0f + r2_0 * _S284);
    float _S286 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S287 = _S286 * u_0;
    float _S288 = 2.0f * u_0;
    float _S289 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S290 = _S289 * u_0;
    float _S291 = 2.0f * v_0;
    float2  _S292 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S293 = uv_0 * _S292;
    float _S294 = (*dist_coeffs_0)[int(4)] * _S292.y;
    float _S295 = (*dist_coeffs_0)[int(5)] * _S292.x;
    float _S296 = _S293.x + _S293.y;
    float _S297 = r2_0 * _S296;
    float _S298 = r2_0 * _S297;
    float _S299 = (*dist_coeffs_0)[int(7)] * _S292.y + _S294 + (*dist_coeffs_0)[int(6)] * _S292.x + _S295 + _S284 * _S296 + _S283 * _S297 + _S282 * _S298 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S298);
    float _S300 = v_0 * _S299;
    float _S301 = u_0 * _S299;
    float2  _S302 = make_float2 (0.0f, 1.0f) + make_float2 ((*dist_coeffs_0)[int(8)] * 0.0f, (*dist_coeffs_0)[int(9)] * 0.0f);
    float2  _S303 = uv_0 * _S302;
    float _S304 = (*dist_coeffs_0)[int(4)] * _S302.y;
    float _S305 = (*dist_coeffs_0)[int(5)] * _S302.x;
    float _S306 = _S303.x + _S303.y;
    float _S307 = r2_0 * _S306;
    float _S308 = r2_0 * _S307;
    float _S309 = (*dist_coeffs_0)[int(7)] * _S302.y + _S304 + (*dist_coeffs_0)[int(6)] * _S302.x + _S305 + _S284 * _S306 + _S283 * _S307 + _S282 * _S308 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S308);
    float _S310 = v_0 * _S309;
    float _S311 = u_0 * _S309;
    return makeMatrix<float, 2, 2> (_S285 * _S292 + make_float2 (_S289 * (v_0 * _S292.y) + _S288 * _S295 + 2.0f * (u_0 * _S295) + _S286 * (v_0 * _S292.x) + _S301 + _S301, _S291 * _S294 + 2.0f * (v_0 * _S294) + _S290 * _S292.y + _S287 * _S292.x + _S300 + _S300), _S285 * _S302 + make_float2 (_S289 * (v_0 * _S302.y) + _S288 * _S305 + 2.0f * (u_0 * _S305) + _S286 * (v_0 * _S302.x) + _S311 + _S311, _S291 * _S304 + 2.0f * (v_0 * _S304) + _S290 * _S302.y + _S287 * _S302.x + _S310 + _S310));
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  * dist_coeffs_1)
{
    Matrix<float, 2, 2>  _S312 = camera_distortion_jac_0(uv_1, dist_coeffs_1);
    return (F32_min((determinant_0(_S312)), ((F32_min((_S312.rows[int(0)].x), (_S312.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_2)
{
    float2  _S313;
    if(is_fisheye_0)
    {
        float r_2 = length_1(uv_2);
        float theta_0 = (F32_atan((r_2)));
        float _S314;
        if(r_2 < 0.00100000004749745f)
        {
            _S314 = 1.0f - r_2 * r_2 / 3.0f;
        }
        else
        {
            _S314 = theta_0 / r_2;
        }
        _S313 = uv_2 * make_float2 (_S314);
    }
    else
    {
        _S313 = uv_2;
    }
    float u_1 = _S313.x;
    float v_1 = _S313.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S315 = _S313 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_2)[int(0)] + r2_1 * ((*dist_coeffs_2)[int(1)] + r2_1 * ((*dist_coeffs_2)[int(2)] + r2_1 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_1 * v_1 + (*dist_coeffs_2)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_2)[int(6)] * r2_1, 2.0f * (*dist_coeffs_2)[int(5)] * u_1 * v_1 + (*dist_coeffs_2)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_2)[int(7)] * r2_1);
    return _S315 + make_float2 ((*dist_coeffs_2)[int(8)] * _S315.x + (*dist_coeffs_2)[int(9)] * _S315.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_3, FixedArray<float, 10>  * dist_coeffs_3, int maxiter_0, float2  * uv_undist_0)
{
    int i_7 = int(0);
    float2  q_0 = uv_3;
    for(;;)
    {
        if(i_7 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S316 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_3 = _S316 + make_float2 ((*dist_coeffs_3)[int(8)] * _S316.x + (*dist_coeffs_3)[int(9)] * _S316.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S317 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S317.rows[int(0)].x * _S317.rows[int(1)].y - _S317.rows[int(0)].y * _S317.rows[int(1)].x);
        float _S318 = r_3.x;
        float _S319 = r_3.y;
        float2  q_1 = q_0 - make_float2 ((_S318 * _S317.rows[int(1)].y - _S319 * _S317.rows[int(0)].y) * inv_det_0, (- _S318 * _S317.rows[int(1)].x + _S319 * _S317.rows[int(0)].x) * inv_det_0);
        i_7 = i_7 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S320 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S321;
    if((F32_min((determinant_0(_S320)), ((F32_min((_S320.rows[int(0)].x), (_S320.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S322 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S321 = (length_1(_S322 + make_float2 ((*dist_coeffs_3)[int(8)] * _S322.x + (*dist_coeffs_3)[int(9)] * _S322.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S321 = false;
    }
    return _S321;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S323 = uv_4;
    bool _S324 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S323);
    if(!_S324)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S325 = _S323;
        float theta_1 = length_1(_S323);
        float _S326;
        if(theta_1 < 0.00100000004749745f)
        {
            _S326 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S326 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S327 = make_float3 ((_S325 * make_float2 (_S326)).x, (_S325 * make_float2 (_S326)).y, (F32_cos((theta_1))));
        raydir_0 = _S327;
    }
    else
    {
        raydir_0 = make_float3 (_S323.x, _S323.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_5, float3  * raydir_1)
{
    float2  _S328 = uv_5;
    int3  _S329 = make_int3 (int(0));
    float3  _S330 = make_float3 ((float)_S329.x, (float)_S329.y, (float)_S329.z);
    *raydir_1 = _S330;
    bool _S331 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S328);
    if(!_S331)
    {
        return false;
    }
    float3  _S332;
    if(is_fisheye_2)
    {
        float2  _S333 = _S328;
        float theta_2 = length_1(_S328);
        float _S334;
        if(theta_2 < 0.00100000004749745f)
        {
            _S334 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S334 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S335 = make_float3 ((_S333 * make_float2 (_S334)).x, (_S333 * make_float2 (_S334)).y, (F32_cos((theta_2))));
        _S332 = _S335;
    }
    else
    {
        _S332 = make_float3 (_S328.x, _S328.y, 1.0f);
    }
    *raydir_1 = _S332;
    return true;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_8)
{
    float _S336 = (*right_2).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  right_d_result_0;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_8.x;
    float sum_4 = _S336 + (*right_2).primal_0.rows[int(0)].y * dOut_8.y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_8.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(0)].z * dOut_8.z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_8.z;
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = sum_5;
    float _S337 = (*right_2).primal_0.rows[int(1)].x * dOut_8.x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_8.x;
    float sum_6 = _S337 + (*right_2).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_8.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(1)].z * dOut_8.z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_8.z;
    *&((&left_d_result_0)->y) = sum_7;
    float _S338 = (*right_2).primal_0.rows[int(2)].x * dOut_8.x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_8.x;
    float sum_8 = _S338 + (*right_2).primal_0.rows[int(2)].y * dOut_8.y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_8.y;
    float sum_9 = sum_8 + (*right_2).primal_0.rows[int(2)].z * dOut_8.z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_8.z;
    *&((&left_d_result_0)->z) = sum_9;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_0;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_2(float3  left_3, Matrix<float, 3, 3>  right_3)
{
    float3  result_11;
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
        int i_8 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(i_8 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_3, i_8) * _slang_vector_get_element(right_3.rows[i_8], j_1);
            i_8 = i_8 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_11, j_1) = sum_10;
        j_1 = j_1 + int(1);
    }
    return result_11;
}

inline __device__ float3  normalize_0(float3  x_16)
{
    return x_16 / make_float3 (length_2(x_16));
}

inline __device__ float4  normalize_1(float4  x_17)
{
    return x_17 / make_float4 (length_0(x_17));
}

inline __device__ bool generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_6, float3  * ray_o_0, float3  * ray_d_0)
{
    float2  _S339 = uv_6;
    *ray_o_0 = - mul_2(t_1, R_2);
    bool _S340 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S339);
    if(!_S340)
    {
        return false;
    }
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float2  _S341 = _S339;
        float theta_3 = length_1(_S339);
        float _S342;
        if(theta_3 < 0.00100000004749745f)
        {
            _S342 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S342 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S343 = make_float3 ((_S341 * make_float2 (_S342)).x, (_S341 * make_float2 (_S342)).y, (F32_cos((theta_3))));
        raydir_2 = _S343;
    }
    else
    {
        raydir_2 = make_float3 (_S339.x, _S339.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_2(raydir_2, R_2));
    return true;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S344;
    bool _S345;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S346, Matrix<float, 3, 3>  _S347)
{
    return mul_2(_S346, _S347);
}

inline __device__ float s_primal_ctx_sin_0(float _S348)
{
    return (F32_sin((_S348)));
}

inline __device__ float s_primal_ctx_cos_0(float _S349)
{
    return (F32_cos((_S349)));
}

inline __device__ bool s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_7, bool is_fisheye_4, FixedArray<float, 10>  * dist_coeffs_7, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S344 = make_float2 (0.0f);
    _s_diff_ctx_0->_S345 = false;
    float3  _S350 = make_float3 (0.0f);
    float3  _S351 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    float2  _S352 = uv_7;
    bool _S353 = undistort_point_0(uv_7, dist_coeffs_7, int(8), &_S352);
    _s_diff_ctx_0->_S344 = _S352;
    _s_diff_ctx_0->_S345 = _S353;
    float2  _S354 = _S352;
    float3  raydir_3;
    bool _S355;
    if(!!_S353)
    {
        if(is_fisheye_4)
        {
            float _S356 = length_1(_S354);
            float _S357;
            if(_S356 < 0.00100000004749745f)
            {
                _S357 = 1.0f - _S356 * _S356 / 6.0f;
            }
            else
            {
                _S357 = s_primal_ctx_sin_0(_S356) / _S356;
            }
            float3  _S358 = make_float3 ((_S354 * make_float2 (_S357)).x, (_S354 * make_float2 (_S357)).y, s_primal_ctx_cos_0(_S356));
            raydir_3 = _S358;
        }
        else
        {
            raydir_3 = make_float3 (_S354.x, _S354.y, 1.0f);
        }
        float3  _S359 = normalize_0(s_primal_ctx_mul_0(raydir_3, dpR_0));
        _S355 = true;
        raydir_3 = _S359;
    }
    else
    {
        _S355 = false;
        raydir_3 = _S350;
    }
    *dpray_o_0 = _S351;
    *dpray_d_0 = raydir_3;
    return _S355;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float _s_dOut_2)
{
    float _S360 = (*dpx_9).primal_0.x;
    float _S361 = (*dpx_9).primal_0.y;
    float _S362 = (*dpx_9).primal_0.z;
    DiffPair_float_0 _S363;
    (&_S363)->primal_0 = _S360 * _S360 + _S361 * _S361 + _S362 * _S362;
    (&_S363)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S363, _s_dOut_2);
    float _S364 = (*dpx_9).primal_0.z * _S363.differential_0;
    float _S365 = _S364 + _S364;
    float _S366 = (*dpx_9).primal_0.y * _S363.differential_0;
    float _S367 = _S366 + _S366;
    float _S368 = (*dpx_9).primal_0.x * _S363.differential_0;
    float _S369 = _S368 + _S368;
    float3  _S370 = make_float3 (0.0f);
    *&((&_S370)->z) = _S365;
    *&((&_S370)->y) = _S367;
    *&((&_S370)->x) = _S369;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S370;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S371, float _S372)
{
    s_bwd_prop_length_impl_1(_S371, _S372);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  _s_dOut_3)
{
    float _S373 = length_2((*dpx_10).primal_0);
    float3  _S374 = (*dpx_10).primal_0 * _s_dOut_3;
    float3  _S375 = make_float3 (1.0f / _S373) * _s_dOut_3;
    float _S376 = - ((_S374.x + _S374.y + _S374.z) / (_S373 * _S373));
    float3  _S377 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S378;
    (&_S378)->primal_0 = (*dpx_10).primal_0;
    (&_S378)->differential_0 = _S377;
    s_bwd_length_impl_1(&_S378, _S376);
    float3  _S379 = _S375 + _S378.differential_0;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S379;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S380, float3  _S381)
{
    s_bwd_prop_normalize_impl_0(_S380, _S381);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S382, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S383, float3  _S384)
{
    _d_mul_0(_S382, _S383, _S384);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_8, bool is_fisheye_5, FixedArray<float, 10>  * dist_coeffs_8, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S385 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S386 = *dpt_1;
    float3  _S387 = make_float3 (0.0f);
    bool _S388 = !!_s_diff_ctx_1->_S345;
    float3  raydir_4;
    float3  _S389;
    if(_S388)
    {
        if(is_fisheye_5)
        {
            float _S390 = length_1(_s_diff_ctx_1->_S344);
            float _S391;
            if(_S390 < 0.00100000004749745f)
            {
                _S391 = 1.0f - _S390 * _S390 / 6.0f;
            }
            else
            {
                _S391 = s_primal_ctx_sin_0(_S390) / _S390;
            }
            float3  _S392 = make_float3 ((_s_diff_ctx_1->_S344 * make_float2 (_S391)).x, (_s_diff_ctx_1->_S344 * make_float2 (_S391)).y, s_primal_ctx_cos_0(_S390));
            raydir_4 = _S392;
        }
        else
        {
            raydir_4 = make_float3 (_s_diff_ctx_1->_S344.x, _s_diff_ctx_1->_S344.y, 1.0f);
        }
        float3  _S393 = raydir_4;
        raydir_4 = s_primal_ctx_mul_0(raydir_4, _S385.primal_0);
        _S389 = _S393;
    }
    else
    {
        raydir_4 = _S387;
        _S389 = _S387;
    }
    Matrix<float, 3, 3>  _S394 = makeMatrix<float, 3, 3> (0.0f);
    Matrix<float, 3, 3>  _S395;
    if(_S388)
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S396;
        (&_S396)->primal_0 = raydir_4;
        (&_S396)->differential_0 = _S387;
        s_bwd_normalize_impl_0(&_S396, dpray_d_1);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S397;
        (&_S397)->primal_0 = _S389;
        (&_S397)->differential_0 = _S387;
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S398;
        (&_S398)->primal_0 = _S385.primal_0;
        (&_S398)->differential_0 = _S394;
        s_bwd_prop_mul_0(&_S397, &_S398, _S396.differential_0);
        _S395 = _S398.differential_0;
    }
    else
    {
        _S395 = _S394;
    }
    float3  _S399 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S400;
    (&_S400)->primal_0 = _S386.primal_0;
    (&_S400)->differential_0 = _S387;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S401;
    (&_S401)->primal_0 = _S385.primal_0;
    (&_S401)->differential_0 = _S394;
    s_bwd_prop_mul_0(&_S400, &_S401, _S399);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S400.differential_0;
    Matrix<float, 3, 3>  _S402 = _S401.differential_0 + _S395;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S402;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S403, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S404, float2  _S405, bool _S406, FixedArray<float, 10>  * _S407, float3  _S408, float3  _S409)
{
    float3  _S410;
    float3  _S411;
    s_bwd_prop_generate_ray_Intermediates_0 _S412;
    bool _S413 = s_primal_ctx_generate_ray_0((*_S403).primal_0, (*_S404).primal_0, _S405, _S406, _S407, &_S410, &_S411, &_S412);
    s_bwd_prop_generate_ray_Intermediates_0 _S414 = _S412;
    s_bwd_prop_generate_ray_0(_S403, _S404, _S405, _S406, _S407, _S408, _S409, &_S414);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_9, bool is_fisheye_6, FixedArray<float, 10>  * dist_coeffs_9, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S415 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S415;
    float3  _S416 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S416;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_9, is_fisheye_6, dist_coeffs_9, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S417 = scale_2.x;
    float sx_0 = (F32_exp((_S417)));
    float _S418 = scale_2.y;
    float sy_0 = (F32_exp((_S418)));
    float sz_0 = scale_2.z - 0.5f * (_S417 + _S418);
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
    Matrix<float, 3, 3>  _S419 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S419, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S419, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S419, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
    return;
}

inline __device__ float4  floor_0(float4  x_19)
{
    float4  result_12;
    int i_9 = int(0);
    for(;;)
    {
        if(i_9 < int(4))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_12, i_9) = (F32_floor((_slang_vector_get_element(x_19, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_12;
}

inline __device__ void mcmc_add_noise_3dgs(float scaler_0, float min_opacity_0, float3  * mean_1, float3  scale_3, float4  quat_6, float opac_0)
{
    float4  _S420 = normalize_1(quat_6);
    float3  _S421 = exp_0(scale_3);
    float x_20 = _S420.y;
    float x2_4 = x_20 * x_20;
    float y2_4 = _S420.z * _S420.z;
    float z2_4 = _S420.w * _S420.w;
    float xy_4 = _S420.y * _S420.z;
    float xz_4 = _S420.y * _S420.w;
    float yz_4 = _S420.z * _S420.w;
    float wx_4 = _S420.x * _S420.y;
    float wy_4 = _S420.x * _S420.z;
    float wz_4 = _S420.x * _S420.w;
    Matrix<float, 3, 3>  M_2 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S421.x, 0.0f, 0.0f, 0.0f, _S421.y, 0.0f, 0.0f, 0.0f, _S421.z));
    float4  _S422 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S423 = _S422 - floor_0(_S422);
    float4  _S424 = _S423 + make_float4 (dot_1(_S423, float4 {_S423.w, _S423.z, _S423.x, _S423.y} + make_float4 (33.3300018310546875f)));
    float4  _S425 = (float4 {_S424.x, _S424.x, _S424.y, _S424.z} + float4 {_S424.y, _S424.z, _S424.z, _S424.w}) * float4 {_S424.z, _S424.y, _S424.w, _S424.x};
    float4  _S426 = _S425 - floor_0(_S425);
    float2  _S427 = float2 {_S426.x, _S426.z};
    float _S428 = 6.28318548202514648f * _S427.y;
    float2  _S429 = float2 {_S426.y, _S426.w};
    float _S430 = 6.28318548202514648f * _S429.y;
    *mean_1 = *mean_1 + mul_0(mul_1(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S427.x))))))) * make_float2 ((F32_cos((_S428))), (F32_sin((_S428))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S427.x))))))) * make_float2 ((F32_cos((_S428))), (F32_sin((_S428))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S429.x))))))) * make_float2 ((F32_cos((_S430))), (F32_sin((_S430))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_9)
{
    float _S431 = dOut_9.y;
    float _S432 = dOut_9.z;
    float _S433 = dOut_9.x;
    float _S434 = (*a_0).primal_0.z * _S431 + - (*a_0).primal_0.y * _S432;
    float _S435 = - (*a_0).primal_0.z * _S433 + (*a_0).primal_0.x * _S432;
    float _S436 = (*a_0).primal_0.y * _S433 + - (*a_0).primal_0.x * _S431;
    float3  _S437 = make_float3 (- (*b_0).primal_0.z * _S431 + (*b_0).primal_0.y * _S432, (*b_0).primal_0.z * _S433 + - (*b_0).primal_0.x * _S432, - (*b_0).primal_0.y * _S433 + (*b_0).primal_0.x * _S431);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S437;
    float3  _S438 = make_float3 (_S434, _S435, _S436);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S438;
    return;
}

inline __device__ float3  cross_0(float3  left_4, float3  right_4)
{
    float _S439 = left_4.y;
    float _S440 = right_4.z;
    float _S441 = left_4.z;
    float _S442 = right_4.y;
    float _S443 = right_4.x;
    float _S444 = left_4.x;
    return make_float3 (_S439 * _S440 - _S441 * _S442, _S441 * _S443 - _S444 * _S440, _S444 * _S442 - _S439 * _S443);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S445 = normalize_1(quat_7);
    float _S446 = scale_4.x;
    float sx_1 = (F32_exp((_S446)));
    float _S447 = scale_4.y;
    float sy_1 = (F32_exp((_S447)));
    float sz_1 = scale_4.z - 0.5f * (_S446 + _S447);
    float x_21 = _S445.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S445.z * _S445.z;
    float z2_5 = _S445.w * _S445.w;
    float xy_5 = _S445.y * _S445.z;
    float xz_5 = _S445.y * _S445.w;
    float yz_5 = _S445.z * _S445.w;
    float wx_5 = _S445.x * _S445.y;
    float wy_5 = _S445.x * _S445.z;
    float wz_5 = _S445.x * _S445.w;
    Matrix<float, 3, 3>  _S448 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S448, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S448, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S448, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S449 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S450 = _S449 - floor_0(_S449);
    float4  _S451 = _S450 + make_float4 (dot_1(_S450, float4 {_S450.w, _S450.z, _S450.x, _S450.y} + make_float4 (33.3300018310546875f)));
    float4  _S452 = (float4 {_S451.x, _S451.x, _S451.y, _S451.z} + float4 {_S451.y, _S451.z, _S451.z, _S451.w}) * float4 {_S451.z, _S451.y, _S451.w, _S451.x};
    float4  _S453 = _S452 - floor_0(_S452);
    float2  _S454 = float2 {_S453.x, _S453.z};
    float _S455 = 6.28318548202514648f * _S454.y;
    float2  _S456 = float2 {_S453.y, _S453.w};
    float _S457 = 6.28318548202514648f * _S456.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S454.x))))))) * make_float2 ((F32_cos((_S455))), (F32_sin((_S455))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S454.x))))))) * make_float2 ((F32_cos((_S455))), (F32_sin((_S455))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S456.x))))))) * make_float2 ((F32_cos((_S457))), (F32_sin((_S457))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_3, float3  dOut_10)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_11).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_3).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_10.x);
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = left_dp_0.differential_0;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_11).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_3).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_10.y);
    *&((&left_d_result_1)->y) = left_dp_1.differential_0;
    *&((&right_d_result_1)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_11).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_3).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_10.z);
    *&((&left_d_result_1)->z) = left_dp_2.differential_0;
    *&((&right_d_result_1)->z) = right_dp_2.differential_0;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = left_d_result_1;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  max_0(float3  x_22, float3  y_7)
{
    float3  result_13;
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
        *_slang_vector_get_element_ptr(&result_13, i_10) = (F32_max((_slang_vector_get_element(x_22, i_10)), (_slang_vector_get_element(y_7, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_13;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, float3  dOut_11)
{
    float3  _S458 = make_float3 (1.0f) / (*dpx_12).primal_0 * dOut_11;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S458;
    return;
}

inline __device__ float3  log_0(float3  x_23)
{
    float3  result_14;
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
        *_slang_vector_get_element_ptr(&result_14, i_11) = (F32_log((_slang_vector_get_element(x_23, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_14;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_13, float3  dOut_12)
{
    float3  _S459 = _slang_select(((*dpx_13).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_13).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_12;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S459;
    return;
}

inline __device__ float3  abs_0(float3  x_24)
{
    float3  result_15;
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
        *_slang_vector_get_element_ptr(&result_15, i_12) = (F32_abs((_slang_vector_get_element(x_24, i_12))));
        i_12 = i_12 + int(1);
    }
    return result_15;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_14, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_13)
{
    DiffPair_float_0 _S460 = *dpx_14;
    bool _S461;
    if(((*dpx_14).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S461 = ((*dpx_14).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S461 = false;
    }
    float _S462;
    if(_S461)
    {
        _S462 = dOut_13;
    }
    else
    {
        _S462 = 0.0f;
    }
    dpx_14->primal_0 = _S460.primal_0;
    dpx_14->differential_0 = _S462;
    DiffPair_float_0 _S463 = *dpMin_0;
    if((_S460.primal_0) < ((*dpMin_0).primal_0))
    {
        _S462 = dOut_13;
    }
    else
    {
        _S462 = 0.0f;
    }
    dpMin_0->primal_0 = _S463.primal_0;
    dpMin_0->differential_0 = _S462;
    DiffPair_float_0 _S464 = *dpMax_0;
    if(((*dpx_14).primal_0) > ((*dpMax_0).primal_0))
    {
        _S462 = dOut_13;
    }
    else
    {
        _S462 = 0.0f;
    }
    dpMax_0->primal_0 = _S464.primal_0;
    dpMax_0->differential_0 = _S462;
    return;
}

inline __device__ float clamp_0(float x_25, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_25), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_15, float dOut_14)
{
    float _S465 = -0.5f / ((*dpx_15).primal_0 * (F32_sqrt(((*dpx_15).primal_0)))) * dOut_14;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S465;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_15)
{
    float _S466 = (1.0f - (*dps_0).primal_0) * dOut_15;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S466;
    DiffPair_float_0 _S467 = *dpy_4;
    float _S468 = (*dps_0).primal_0 * dOut_15;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S468;
    float _S469 = (_S467.primal_0 - (*dpx_16).primal_0) * dOut_15;
    dps_0->primal_0 = _S467.primal_0;
    dps_0->differential_0 = _S469;
    return;
}

inline __device__ float lerp_0(float x_26, float y_8, float s_4)
{
    return x_26 + (y_8 - x_26) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 11>  * weights_0, FixedArray<float, 23>  * _S470)
{
    float3  _S471;
    bool _S472;
    bool _S473;
    FixedArray<float, 23>  losses_1;
    float _S474 = float(mask_0);
    float _S475 = 1.0f - (*weights_0)[int(0)];
    float3  _S476 = make_float3 (0.0f);
    float3  _S477 = abs_0(ref_rgb_0 * make_float3 (_S475) + (log_0(max_0(ref_rgb_0, _S476) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)]) - (render_rgb_0 * make_float3 (_S475) + (log_0(max_0(render_rgb_0, _S476) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)])));
    losses_1[int(0)] = (*weights_0)[int(1)] * _S474 * ((_S477.x + _S477.y + _S477.z) * 0.3333333432674408f);
    float3  _S478 = ref_rgb_0 - render_rgb_0;
    losses_1[int(1)] = _S474 * clamp_0(dot_0(_S478, _S478) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S479 = float(depth_mask_0 & mask_0);
    float _S480 = _S479 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S481 = _S479 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S480;
    losses_1[int(3)] = _S481;
    losses_1[int(4)] = _S480 * _S480;
    losses_1[int(5)] = _S481 * _S481;
    losses_1[int(6)] = _S480 * _S481;
    bool _S482 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S483 = norm2_0 == 0.0f;
        _S472 = _S483;
        if(_S483)
        {
            _S471 = _S476;
            break;
        }
        _S471 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S484;
    bool _S485 = !_S472;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S486 = norm2_1 == 0.0f;
        _S473 = _S486;
        if(_S486)
        {
            _S484 = _S476;
            break;
        }
        _S484 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S487;
    float3  _S488;
    bool _S489 = !_S473;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S488 = _S476;
            _S487 = false;
            break;
        }
        _S488 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S487 = _S482;
        break;
    }
    float _S490 = float(_S485 & _S487);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S471, _S488);
    losses_1[int(7)] = (*weights_0)[int(3)] * _S490 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S491 = float(_S489 & _S487);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S484, _S488);
    losses_1[int(8)] = (*weights_0)[int(3)] * _S491 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S492 = float(_S485 & _S489);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S471, _S484);
    losses_1[int(11)] = (*weights_0)[int(6)] * _S492 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S493 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S494 = float(alpha_mask_0);
    float _S495 = float(ref_alpha_0);
    float _S496 = (F32_max((_S493), (_S495)));
    losses_1[int(9)] = (*weights_0)[int(4)] * _S494 * - lerp_0((F32_log(((F32_max((1.0f - _S496), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S496), (9.99999997475242708e-07f)))))), _S495);
    float _S497 = 1.0f - _S493;
    float _S498 = 1.0f - _S495;
    float _S499 = (F32_max((_S497), (_S498)));
    losses_1[int(10)] = (*weights_0)[int(5)] * _S494 * - lerp_0((F32_log(((F32_max((1.0f - _S499), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S499), (9.99999997475242708e-07f)))))), _S498);
    losses_1[int(12)] = (*weights_0)[int(7)] * 4.0f * _S493 * _S497;
    float _S500 = (F32_max((_S493), (9.999999960041972e-13f)));
    losses_1[int(13)] = (*weights_0)[int(8)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S500;
    losses_1[int(14)] = (*weights_0)[int(9)] * depth_dist_0 / _S500;
    losses_1[int(15)] = (*weights_0)[int(10)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S500;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S474;
    losses_1[int(18)] = _S479;
    losses_1[int(19)] = _S490;
    losses_1[int(20)] = _S491;
    losses_1[int(21)] = _S492;
    losses_1[int(22)] = _S494;
    *_S470 = losses_1;
    return;
}

inline __device__ float3  s_primal_ctx_max_1(float3  _S501, float3  _S502)
{
    return max_0(_S501, _S502);
}

inline __device__ float3  s_primal_ctx_log_1(float3  _S503)
{
    return log_0(_S503);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S504, float3  _S505)
{
    return dot_0(_S504, _S505);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S506)
{
    return (F32_rsqrt((_S506)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S507, float _S508, float _S509)
{
    return clamp_0(_S507, _S508, _S509);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S510, DiffPair_float_0 * _S511, DiffPair_float_0 * _S512, float _S513)
{
    _d_lerp_0(_S510, _S511, _S512, _S513);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S514, DiffPair_float_0 * _S515, DiffPair_float_0 * _S516, float _S517)
{
    _d_clamp_0(_S514, _S515, _S516, _S517);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S518, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S519, float _S520)
{
    _d_dot_0(_S518, _S519, _S520);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S521, float _S522)
{
    _d_rsqrt_0(_S521, _S522);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S523, float3  _S524)
{
    _d_abs_vector_0(_S523, _S524);
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S525, float3  _S526)
{
    _d_log_vector_0(_S525, _S526);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S527, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S528, float3  _S529)
{
    _d_max_vector_0(_S527, _S528, _S529);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 11>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S530 = *dprender_rgb_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S531 = *dpref_rgb_0;
    DiffPair_float_0 _S532 = *dprender_depth_0;
    DiffPair_float_0 _S533 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S534 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S535 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S536 = *dpref_normal_0;
    DiffPair_float_0 _S537 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S538 = *dprgb_dist_0;
    DiffPair_float_0 _S539 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S540 = *dpnormal_dist_0;
    float3  _S541 = make_float3 (0.0f);
    float _S542 = float(mask_1);
    float _S543 = (*weights_1)[int(1)] * _S542;
    float3  _S544 = make_float3 ((*weights_1)[int(0)]);
    float _S545 = 1.0f - (*weights_1)[int(0)];
    float3  _S546 = make_float3 (_S545);
    float3  _S547 = make_float3 (0.0f);
    float3  _S548 = s_primal_ctx_max_1((*dprender_rgb_0).primal_0, _S547) + make_float3 (0.00392156885936856f);
    float3  _S549 = s_primal_ctx_max_1((*dpref_rgb_0).primal_0, _S547) + make_float3 (0.00392156885936856f);
    float3  _S550 = (*dpref_rgb_0).primal_0 * make_float3 (_S545) + (s_primal_ctx_log_1(_S549) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]) - ((*dprender_rgb_0).primal_0 * make_float3 (_S545) + (s_primal_ctx_log_1(_S548) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]));
    float3  _S551 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S552 = s_primal_ctx_dot_0(_S551, _S551) * 0.3333333432674408f;
    float _S553 = float(depth_mask_1 & mask_1);
    float _S554 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S555 = _S553 * s_primal_ctx_log_0(_S554);
    float _S556 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S557 = _S553 * s_primal_ctx_log_0(_S556);
    bool _S558 = normal_mask_1 & mask_1;
    float _S559 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S560 = !(_S559 == 0.0f);
    float3  _S561;
    float3  _S562;
    if(_S560)
    {
        float _S563 = s_primal_ctx_rsqrt_0(_S559);
        float3  _S564 = make_float3 (_S563);
        _S561 = _S534.primal_0 * make_float3 (_S563);
        _S562 = _S564;
    }
    else
    {
        _S561 = _S547;
        _S562 = _S541;
    }
    float _S565 = s_primal_ctx_dot_0(_S535.primal_0, _S535.primal_0);
    bool _S566 = !(_S565 == 0.0f);
    float3  _S567;
    float3  _S568;
    if(_S566)
    {
        float _S569 = s_primal_ctx_rsqrt_0(_S565);
        float3  _S570 = make_float3 (_S569);
        _S567 = _S535.primal_0 * make_float3 (_S569);
        _S568 = _S570;
    }
    else
    {
        _S567 = _S547;
        _S568 = _S541;
    }
    float _S571 = s_primal_ctx_dot_0(_S536.primal_0, _S536.primal_0);
    bool _S572 = _S571 == 0.0f;
    bool _S573;
    if(_S572)
    {
        _S573 = false;
    }
    else
    {
        _S573 = _S558;
    }
    bool _S574 = !_S572;
    float3  _S575;
    float3  _S576;
    if(_S574)
    {
        float _S577 = s_primal_ctx_rsqrt_0(_S571);
        float3  _S578 = make_float3 (_S577);
        _S575 = _S536.primal_0 * make_float3 (_S577);
        _S576 = _S578;
    }
    else
    {
        _S575 = _S547;
        _S576 = _S541;
    }
    float _S579 = (*weights_1)[int(3)] * float(_S560 & _S573);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S561, _S575);
    float _S580 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S581 = (*weights_1)[int(3)] * float(_S566 & _S573);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S567, _S575);
    float _S582 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S583 = (*weights_1)[int(6)] * float(_S560 & _S566);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S561, _S567);
    float _S584 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S585 = s_primal_ctx_clamp_0(_S537.primal_0, 0.0f, 1.0f);
    float _S586 = float(alpha_mask_1);
    float _S587 = (*weights_1)[int(4)] * _S586;
    float _S588 = float(ref_alpha_1);
    float _S589 = s_primal_ctx_max_0(_S585, _S588);
    float _S590 = 1.0f - _S589;
    float _S591 = s_primal_ctx_max_0(_S590, 9.99999997475242708e-07f);
    float _S592 = s_primal_ctx_log_0(_S591);
    float _S593 = s_primal_ctx_max_0(_S589, 9.99999997475242708e-07f);
    float _S594 = s_primal_ctx_log_0(_S593);
    float _S595 = (*weights_1)[int(5)] * _S586;
    float _S596 = 1.0f - _S585;
    float _S597 = 1.0f - _S588;
    float _S598 = s_primal_ctx_max_0(_S596, _S597);
    float _S599 = 1.0f - _S598;
    float _S600 = s_primal_ctx_max_0(_S599, 9.99999997475242708e-07f);
    float _S601 = s_primal_ctx_log_0(_S600);
    float _S602 = s_primal_ctx_max_0(_S598, 9.99999997475242708e-07f);
    float _S603 = s_primal_ctx_log_0(_S602);
    float _S604 = (*weights_1)[int(7)] * 4.0f;
    float _S605 = _S604 * _S585;
    float _S606 = s_primal_ctx_max_0(_S585, 9.999999960041972e-13f);
    float _S607 = _S606 * _S606;
    float _S608 = (*_s_dOut_4)[int(15)] / _S607;
    float _S609 = 0.3333333432674408f * ((*weights_1)[int(10)] * (_S606 * _S608));
    float _S610 = (*_s_dOut_4)[int(14)] / _S607;
    float _S611 = (*weights_1)[int(9)] * (_S606 * _S610);
    float _S612 = (*_s_dOut_4)[int(13)] / _S607;
    float _S613 = _S606 * _S612;
    float _S614 = (*weights_1)[int(10)] * ((_S540.primal_0.x + _S540.primal_0.y + _S540.primal_0.z) * 0.3333333432674408f) * - _S608 + (*weights_1)[int(9)] * _S539.primal_0 * - _S610 + (*weights_1)[int(8)] * ((_S538.primal_0.x + _S538.primal_0.y + _S538.primal_0.z) * 0.3333333432674408f) * - _S612;
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = _S585;
    (&_S615)->differential_0 = 0.0f;
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = 9.999999960041972e-13f;
    (&_S616)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S615, &_S616, _S614);
    float _S617 = 0.3333333432674408f * ((*weights_1)[int(8)] * _S613);
    float _S618 = _S605 * (*_s_dOut_4)[int(12)];
    float _S619 = _S604 * (_S596 * (*_s_dOut_4)[int(12)]);
    float _S620 = - (_S595 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S601;
    (&_S621)->differential_0 = 0.0f;
    DiffPair_float_0 _S622;
    (&_S622)->primal_0 = _S603;
    (&_S622)->differential_0 = 0.0f;
    DiffPair_float_0 _S623;
    (&_S623)->primal_0 = _S597;
    (&_S623)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S621, &_S622, &_S623, _S620);
    DiffPair_float_0 _S624;
    (&_S624)->primal_0 = _S602;
    (&_S624)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S624, _S622.differential_0);
    DiffPair_float_0 _S625;
    (&_S625)->primal_0 = _S598;
    (&_S625)->differential_0 = 0.0f;
    DiffPair_float_0 _S626;
    (&_S626)->primal_0 = 9.99999997475242708e-07f;
    (&_S626)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S625, &_S626, _S624.differential_0);
    DiffPair_float_0 _S627;
    (&_S627)->primal_0 = _S600;
    (&_S627)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S627, _S621.differential_0);
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = _S599;
    (&_S628)->differential_0 = 0.0f;
    DiffPair_float_0 _S629;
    (&_S629)->primal_0 = 9.99999997475242708e-07f;
    (&_S629)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S628, &_S629, _S627.differential_0);
    float _S630 = _S625.differential_0 + - _S628.differential_0;
    DiffPair_float_0 _S631;
    (&_S631)->primal_0 = _S596;
    (&_S631)->differential_0 = 0.0f;
    DiffPair_float_0 _S632;
    (&_S632)->primal_0 = _S597;
    (&_S632)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S631, &_S632, _S630);
    float _S633 = - (_S618 + _S631.differential_0);
    float _S634 = - (_S587 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S635;
    (&_S635)->primal_0 = _S592;
    (&_S635)->differential_0 = 0.0f;
    DiffPair_float_0 _S636;
    (&_S636)->primal_0 = _S594;
    (&_S636)->differential_0 = 0.0f;
    DiffPair_float_0 _S637;
    (&_S637)->primal_0 = _S588;
    (&_S637)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S635, &_S636, &_S637, _S634);
    DiffPair_float_0 _S638;
    (&_S638)->primal_0 = _S593;
    (&_S638)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S638, _S636.differential_0);
    DiffPair_float_0 _S639;
    (&_S639)->primal_0 = _S589;
    (&_S639)->differential_0 = 0.0f;
    DiffPair_float_0 _S640;
    (&_S640)->primal_0 = 9.99999997475242708e-07f;
    (&_S640)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S639, &_S640, _S638.differential_0);
    DiffPair_float_0 _S641;
    (&_S641)->primal_0 = _S591;
    (&_S641)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S641, _S635.differential_0);
    DiffPair_float_0 _S642;
    (&_S642)->primal_0 = _S590;
    (&_S642)->differential_0 = 0.0f;
    DiffPair_float_0 _S643;
    (&_S643)->primal_0 = 9.99999997475242708e-07f;
    (&_S643)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S642, &_S643, _S641.differential_0);
    float _S644 = _S639.differential_0 + - _S642.differential_0;
    DiffPair_float_0 _S645;
    (&_S645)->primal_0 = _S585;
    (&_S645)->differential_0 = 0.0f;
    DiffPair_float_0 _S646;
    (&_S646)->primal_0 = _S588;
    (&_S646)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S645, &_S646, _S644);
    float _S647 = _S615.differential_0 + _S619 + _S633 + _S645.differential_0;
    DiffPair_float_0 _S648;
    (&_S648)->primal_0 = _S537.primal_0;
    (&_S648)->differential_0 = 0.0f;
    DiffPair_float_0 _S649;
    (&_S649)->primal_0 = 0.0f;
    (&_S649)->differential_0 = 0.0f;
    DiffPair_float_0 _S650;
    (&_S650)->primal_0 = 1.0f;
    (&_S650)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S648, &_S649, &_S650, _S647);
    DiffPair_float_0 _S651 = _S648;
    float _S652 = _S583 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = _S584;
    (&_S653)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S653, _S652);
    DiffPair_float_0 _S654;
    (&_S654)->primal_0 = cos_sim_loss_5;
    (&_S654)->differential_0 = 0.0f;
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = 9.999999960041972e-13f;
    (&_S655)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S654, &_S655, _S653.differential_0);
    float _S656 = 0.5f * - (_S652 + _S654.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S657;
    (&_S657)->primal_0 = _S561;
    (&_S657)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S658;
    (&_S658)->primal_0 = _S567;
    (&_S658)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S657, &_S658, _S656);
    float _S659 = _S581 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S660;
    (&_S660)->primal_0 = _S582;
    (&_S660)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S660, _S659);
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = cos_sim_loss_4;
    (&_S661)->differential_0 = 0.0f;
    DiffPair_float_0 _S662;
    (&_S662)->primal_0 = 9.999999960041972e-13f;
    (&_S662)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S661, &_S662, _S660.differential_0);
    float _S663 = 0.5f * - (_S659 + _S661.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S664;
    (&_S664)->primal_0 = _S567;
    (&_S664)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S665;
    (&_S665)->primal_0 = _S575;
    (&_S665)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S664, &_S665, _S663);
    float _S666 = _S579 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S667;
    (&_S667)->primal_0 = _S580;
    (&_S667)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S667, _S666);
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = cos_sim_loss_3;
    (&_S668)->differential_0 = 0.0f;
    DiffPair_float_0 _S669;
    (&_S669)->primal_0 = 9.999999960041972e-13f;
    (&_S669)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S668, &_S669, _S667.differential_0);
    float _S670 = 0.5f * - (_S666 + _S668.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S671;
    (&_S671)->primal_0 = _S561;
    (&_S671)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S672;
    (&_S672)->primal_0 = _S575;
    (&_S672)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S671, &_S672, _S670);
    float3  _S673 = _S665.differential_0 + _S672.differential_0;
    float3  _S674 = _S657.differential_0 + _S671.differential_0;
    float3  _S675 = make_float3 (_S609, _S609, _S609);
    float3  _S676 = make_float3 (_S617, _S617, _S617);
    float3  _S677 = _S658.differential_0 + _S664.differential_0;
    float _S678;
    if(_S574)
    {
        float3  _S679 = _S536.primal_0 * _S673;
        float3  _S680 = _S576 * _S673;
        float _S681 = _S679.x + _S679.y + _S679.z;
        DiffPair_float_0 _S682;
        (&_S682)->primal_0 = _S571;
        (&_S682)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S682, _S681);
        _S678 = _S682.differential_0;
        _S561 = _S680;
    }
    else
    {
        _S678 = 0.0f;
        _S561 = _S541;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S683;
    (&_S683)->primal_0 = _S536.primal_0;
    (&_S683)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S684;
    (&_S684)->primal_0 = _S536.primal_0;
    (&_S684)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S683, &_S684, _S678);
    float3  _S685 = _S684.differential_0 + _S683.differential_0 + _S561;
    if(_S566)
    {
        float3  _S686 = _S535.primal_0 * _S677;
        float3  _S687 = _S568 * _S677;
        float _S688 = _S686.x + _S686.y + _S686.z;
        DiffPair_float_0 _S689;
        (&_S689)->primal_0 = _S565;
        (&_S689)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S689, _S688);
        _S678 = _S689.differential_0;
        _S561 = _S687;
    }
    else
    {
        _S678 = 0.0f;
        _S561 = _S541;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S690;
    (&_S690)->primal_0 = _S535.primal_0;
    (&_S690)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S691;
    (&_S691)->primal_0 = _S535.primal_0;
    (&_S691)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S690, &_S691, _S678);
    float3  _S692 = _S691.differential_0 + _S690.differential_0 + _S561;
    if(_S560)
    {
        float3  _S693 = _S534.primal_0 * _S674;
        float3  _S694 = _S562 * _S674;
        float _S695 = _S693.x + _S693.y + _S693.z;
        DiffPair_float_0 _S696;
        (&_S696)->primal_0 = _S559;
        (&_S696)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S696, _S695);
        _S678 = _S696.differential_0;
        _S561 = _S694;
    }
    else
    {
        _S678 = 0.0f;
        _S561 = _S541;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S697;
    (&_S697)->primal_0 = _S534.primal_0;
    (&_S697)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S698;
    (&_S698)->primal_0 = _S534.primal_0;
    (&_S698)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S697, &_S698, _S678);
    float _S699 = _S557 * (*_s_dOut_4)[int(6)];
    float _S700 = _S557 * (*_s_dOut_4)[int(5)];
    float _S701 = _S555 * (*_s_dOut_4)[int(4)];
    float _S702 = _S553 * (_S555 * (*_s_dOut_4)[int(6)] + _S700 + _S700 + (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S703;
    (&_S703)->primal_0 = _S556;
    (&_S703)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S703, _S702);
    DiffPair_float_0 _S704;
    (&_S704)->primal_0 = _S533.primal_0;
    (&_S704)->differential_0 = 0.0f;
    DiffPair_float_0 _S705;
    (&_S705)->primal_0 = 0.00009999999747379f;
    (&_S705)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S704, &_S705, _S703.differential_0);
    float _S706 = _S553 * (_S699 + _S701 + _S701 + (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S707;
    (&_S707)->primal_0 = _S554;
    (&_S707)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S707, _S706);
    DiffPair_float_0 _S708;
    (&_S708)->primal_0 = _S532.primal_0;
    (&_S708)->differential_0 = 0.0f;
    DiffPair_float_0 _S709;
    (&_S709)->primal_0 = 0.00009999999747379f;
    (&_S709)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S708, &_S709, _S707.differential_0);
    float _S710 = _S542 * (*_s_dOut_4)[int(1)];
    DiffPair_float_0 _S711;
    (&_S711)->primal_0 = _S552;
    (&_S711)->differential_0 = 0.0f;
    DiffPair_float_0 _S712;
    (&_S712)->primal_0 = 0.0f;
    (&_S712)->differential_0 = 0.0f;
    DiffPair_float_0 _S713;
    (&_S713)->primal_0 = 1.0f;
    (&_S713)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S711, &_S712, &_S713, _S710);
    float _S714 = 0.3333333432674408f * _S711.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S715;
    (&_S715)->primal_0 = _S551;
    (&_S715)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S716;
    (&_S716)->primal_0 = _S551;
    (&_S716)->differential_0 = _S541;
    s_bwd_prop_dot_0(&_S715, &_S716, _S714);
    float3  _S717 = _S716.differential_0 + _S715.differential_0;
    float3  _S718 = - _S717;
    float _S719 = 0.3333333432674408f * (_S543 * (*_s_dOut_4)[int(0)]);
    float3  _S720 = make_float3 (_S719, _S719, _S719);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S721;
    (&_S721)->primal_0 = _S550;
    (&_S721)->differential_0 = _S541;
    s_bwd_prop_abs_0(&_S721, _S720);
    float3  _S722 = - _S721.differential_0;
    float3  _S723 = _S544 * _S721.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S724;
    (&_S724)->primal_0 = _S549;
    (&_S724)->differential_0 = _S541;
    s_bwd_prop_log_1(&_S724, _S723);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S725;
    (&_S725)->primal_0 = _S531.primal_0;
    (&_S725)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S726;
    (&_S726)->primal_0 = _S547;
    (&_S726)->differential_0 = _S541;
    s_bwd_prop_max_1(&_S725, &_S726, _S724.differential_0);
    float3  _S727 = _S546 * _S721.differential_0;
    float3  _S728 = _S544 * _S722;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S729;
    (&_S729)->primal_0 = _S548;
    (&_S729)->differential_0 = _S541;
    s_bwd_prop_log_1(&_S729, _S728);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S730;
    (&_S730)->primal_0 = _S530.primal_0;
    (&_S730)->differential_0 = _S541;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S731;
    (&_S731)->primal_0 = _S547;
    (&_S731)->differential_0 = _S541;
    s_bwd_prop_max_1(&_S730, &_S731, _S729.differential_0);
    float3  _S732 = _S546 * _S722;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S675;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S611;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S676;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S651.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S685;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S692;
    float3  _S733 = _S698.differential_0 + _S697.differential_0 + _S561;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S733;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S704.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S708.differential_0;
    float3  _S734 = _S717 + _S725.differential_0 + _S727;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S734;
    float3  _S735 = _S718 + _S730.differential_0 + _S732;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S735;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S736, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S737, DiffPair_float_0 * _S738, DiffPair_float_0 * _S739, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S740, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S741, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S742, DiffPair_float_0 * _S743, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S744, DiffPair_float_0 * _S745, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S746, bool _S747, bool _S748, bool _S749, bool _S750, bool _S751, FixedArray<float, 11>  * _S752, FixedArray<float, 23>  * _S753)
{
    s_bwd_prop_per_pixel_losses_0(_S736, _S737, _S738, _S739, _S740, _S741, _S742, _S743, _S744, _S745, _S746, _S747, _S748, _S749, _S750, _S751, _S752, _S753);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 11>  * weights_2, FixedArray<float, 23>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S754 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S754;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S754;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S754;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S754;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S754;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S754;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S754;
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_17, float dOut_16)
{
    float _S755 = 1.0f / ((*dpx_17).primal_0 * 52.30258560180664062f) * dOut_16;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S755;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  * raw_losses_0, FixedArray<float, 11>  * weights_3, FixedArray<float, 10>  * _S756)
{
    FixedArray<float, 10>  losses_2;
    float _S757 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*raw_losses_0)[int(0)] / _S757;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S757)));
    bool _S758;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S758 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S758 = false;
    }
    float _S759;
    if(_S758)
    {
        _S759 = (*weights_3)[int(2)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S759 = 0.0f;
    }
    losses_2[int(2)] = _S759;
    losses_2[int(3)] = ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = ((*raw_losses_0)[int(9)] + (*raw_losses_0)[int(10)]) / (F32_max(((*raw_losses_0)[int(22)]), (1.0f)));
    losses_2[int(5)] = (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    float _S760 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(6)] = (*raw_losses_0)[int(12)] / _S760;
    losses_2[int(7)] = (*raw_losses_0)[int(13)] / _S760;
    losses_2[int(8)] = (*raw_losses_0)[int(14)] / _S760;
    losses_2[int(9)] = (*raw_losses_0)[int(15)] / _S760;
    *_S756 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S761)
{
    return (F32_sqrt((_S761)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S762, float _S763)
{
    _d_log10_0(_S762, _S763);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 11>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S764 = dpraw_losses_0->primal_0;
    float _S765 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S766 = _S765 * _S765;
    float _S767 = dpraw_losses_0->primal_0[int(1)] / _S765;
    bool _S768 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S769;
    if(_S768)
    {
        _S769 = (_S764[int(3)]) != 0.0f;
    }
    else
    {
        _S769 = false;
    }
    float _S770;
    float _S771;
    float _S772;
    float _S773;
    float _S774;
    float _S775;
    float _S776;
    float _S777;
    float _S778;
    float _S779;
    float _S780;
    float _S781;
    float _S782;
    float _S783;
    float _S784;
    if(_S769)
    {
        float _S785 = _S764[int(2)] * _S764[int(3)];
        float _S786 = _S764[int(18)] * _S764[int(18)];
        float _S787 = _S764[int(6)] - _S785 / _S764[int(18)];
        float _S788 = _S764[int(2)] * _S764[int(2)];
        float _S789 = _S764[int(4)] - _S788 / _S764[int(18)];
        float _S790 = _S764[int(3)] * _S764[int(3)];
        float _S791 = _S764[int(5)] - _S790 / _S764[int(18)];
        float _S792 = _S789 * _S791 + 1.0f;
        float _S793 = s_primal_ctx_max_0(9.999999960041972e-13f, _S792);
        float _S794 = s_primal_ctx_sqrt_0(_S793);
        float _S795 = _S794 * _S794;
        float _S796 = 1.0f - _S787 / _S794;
        _S770 = (*weights_4)[int(2)];
        _S771 = _S796;
        _S772 = _S795;
        _S773 = _S787;
        _S774 = _S794;
        _S775 = _S793;
        _S776 = _S792;
        _S777 = _S789;
        _S778 = _S791;
        _S779 = _S786;
        _S780 = _S790;
        _S781 = _S764[int(3)];
        _S782 = _S788;
        _S783 = _S764[int(2)];
        _S784 = _S785;
    }
    else
    {
        _S770 = 0.0f;
        _S771 = 0.0f;
        _S772 = 0.0f;
        _S773 = 0.0f;
        _S774 = 0.0f;
        _S775 = 0.0f;
        _S776 = 0.0f;
        _S777 = 0.0f;
        _S778 = 0.0f;
        _S779 = 0.0f;
        _S780 = 0.0f;
        _S781 = 0.0f;
        _S782 = 0.0f;
        _S783 = 0.0f;
        _S784 = 0.0f;
    }
    float _S797 = s_primal_ctx_max_0(_S764[int(19)], 1.0f);
    float _S798 = _S797 * _S797;
    float _S799 = s_primal_ctx_max_0(_S764[int(20)], 1.0f);
    float _S800 = _S799 * _S799;
    float _S801 = float((I32_max((int((_S764[int(19)]) > 0.5f) + int((_S764[int(20)]) > 0.5f)), (int(1)))));
    float _S802 = _S764[int(9)] + _S764[int(10)];
    float _S803 = s_primal_ctx_max_0(_S764[int(22)], 1.0f);
    float _S804 = _S803 * _S803;
    float _S805 = s_primal_ctx_max_0(_S764[int(21)], 1.0f);
    float _S806 = _S805 * _S805;
    float _S807 = s_primal_ctx_max_0(_S764[int(16)], 1.0f);
    float _S808 = _S807 * _S807;
    float _S809 = (*_s_dOut_5)[int(9)] / _S808;
    float _S810 = _S807 * _S809;
    float _S811 = (*_s_dOut_5)[int(8)] / _S808;
    float _S812 = _S807 * _S811;
    float _S813 = (*_s_dOut_5)[int(7)] / _S808;
    float _S814 = _S807 * _S813;
    float _S815 = (*_s_dOut_5)[int(6)] / _S808;
    float _S816 = _S807 * _S815;
    float _S817 = _S764[int(15)] * - _S809 + _S764[int(14)] * - _S811 + _S764[int(13)] * - _S813 + _S764[int(12)] * - _S815;
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = _S764[int(16)];
    (&_S818)->differential_0 = 0.0f;
    DiffPair_float_0 _S819;
    (&_S819)->primal_0 = 1.0f;
    (&_S819)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S818, &_S819, _S817);
    float _S820 = (*_s_dOut_5)[int(5)] / _S806;
    float _S821 = _S764[int(11)] * - _S820;
    float _S822 = _S805 * _S820;
    DiffPair_float_0 _S823;
    (&_S823)->primal_0 = _S764[int(21)];
    (&_S823)->differential_0 = 0.0f;
    DiffPair_float_0 _S824;
    (&_S824)->primal_0 = 1.0f;
    (&_S824)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S823, &_S824, _S821);
    float _S825 = (*_s_dOut_5)[int(4)] / _S804;
    float _S826 = _S802 * - _S825;
    float _S827 = _S803 * _S825;
    DiffPair_float_0 _S828;
    (&_S828)->primal_0 = _S764[int(22)];
    (&_S828)->differential_0 = 0.0f;
    DiffPair_float_0 _S829;
    (&_S829)->primal_0 = 1.0f;
    (&_S829)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S828, &_S829, _S826);
    float _S830 = (*_s_dOut_5)[int(3)] / _S801;
    float _S831 = _S830 / _S800;
    float _S832 = _S764[int(8)] * - _S831;
    float _S833 = _S799 * _S831;
    DiffPair_float_0 _S834;
    (&_S834)->primal_0 = _S764[int(20)];
    (&_S834)->differential_0 = 0.0f;
    DiffPair_float_0 _S835;
    (&_S835)->primal_0 = 1.0f;
    (&_S835)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S834, &_S835, _S832);
    float _S836 = _S830 / _S798;
    float _S837 = _S764[int(7)] * - _S836;
    float _S838 = _S797 * _S836;
    DiffPair_float_0 _S839;
    (&_S839)->primal_0 = _S764[int(19)];
    (&_S839)->differential_0 = 0.0f;
    DiffPair_float_0 _S840;
    (&_S840)->primal_0 = 1.0f;
    (&_S840)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S839, &_S840, _S837);
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
    _S841[int(15)] = _S810;
    _S841[int(14)] = _S812;
    _S841[int(13)] = _S814;
    _S841[int(16)] = _S818.differential_0;
    _S841[int(12)] = _S816;
    _S841[int(21)] = _S823.differential_0;
    _S841[int(11)] = _S822;
    _S841[int(22)] = _S828.differential_0;
    _S841[int(10)] = _S827;
    _S841[int(9)] = _S827;
    _S841[int(20)] = _S834.differential_0;
    _S841[int(8)] = _S833;
    _S841[int(19)] = _S839.differential_0;
    _S841[int(7)] = _S838;
    float _S842 = _S841[int(0)];
    float _S843 = _S841[int(1)];
    float _S844 = _S841[int(2)];
    float _S845 = _S841[int(3)];
    float _S846 = _S841[int(4)];
    float _S847 = _S841[int(5)];
    float _S848 = _S841[int(6)];
    float _S849 = _S841[int(7)];
    float _S850 = _S841[int(8)];
    float _S851 = _S841[int(9)];
    float _S852 = _S841[int(10)];
    float _S853 = _S841[int(11)];
    float _S854 = _S841[int(12)];
    float _S855 = _S841[int(13)];
    float _S856 = _S841[int(14)];
    float _S857 = _S841[int(15)];
    float _S858 = _S841[int(16)];
    float _S859 = _S841[int(17)];
    float _S860 = _S841[int(18)];
    float _S861 = _S841[int(19)];
    float _S862 = _S841[int(20)];
    float _S863 = _S841[int(21)];
    float _S864 = _S841[int(22)];
    FixedArray<float, 23>  _S865;
    if(_S769)
    {
        float _S866 = _S770 * (*_s_dOut_5)[int(2)];
        DiffPair_float_0 _S867;
        (&_S867)->primal_0 = _S771;
        (&_S867)->differential_0 = 0.0f;
        DiffPair_float_0 _S868;
        (&_S868)->primal_0 = 0.0f;
        (&_S868)->differential_0 = 0.0f;
        DiffPair_float_0 _S869;
        (&_S869)->primal_0 = 2.0f;
        (&_S869)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S867, &_S868, &_S869, _S866);
        float _S870 = - _S867.differential_0 / _S772;
        float _S871 = _S773 * - _S870;
        float _S872 = _S774 * _S870;
        DiffPair_float_0 _S873;
        (&_S873)->primal_0 = _S775;
        (&_S873)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S873, _S871);
        DiffPair_float_0 _S874;
        (&_S874)->primal_0 = 9.999999960041972e-13f;
        (&_S874)->differential_0 = 0.0f;
        DiffPair_float_0 _S875;
        (&_S875)->primal_0 = _S776;
        (&_S875)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S874, &_S875, _S873.differential_0);
        float _S876 = _S777 * _S875.differential_0;
        float _S877 = _S778 * _S875.differential_0;
        float _S878 = - _S876 / _S779;
        float _S879 = _S781 * (_S764[int(18)] * _S878);
        float _S880 = - _S877 / _S779;
        float _S881 = _S783 * (_S764[int(18)] * _S880);
        float _S882 = - _S872 / _S779;
        float _S883 = _S764[int(18)] * _S882;
        float _S884 = _S879 + _S879 + _S783 * _S883;
        float _S885 = _S881 + _S881 + _S781 * _S883;
        float _S886 = _S780 * - _S878 + _S782 * - _S880 + _S784 * - _S882;
        FixedArray<float, 23>  _S887;
        _S887[int(0)] = 0.0f;
        _S887[int(1)] = 0.0f;
        _S887[int(2)] = 0.0f;
        _S887[int(3)] = 0.0f;
        _S887[int(4)] = 0.0f;
        _S887[int(5)] = 0.0f;
        _S887[int(6)] = 0.0f;
        _S887[int(7)] = 0.0f;
        _S887[int(8)] = 0.0f;
        _S887[int(9)] = 0.0f;
        _S887[int(10)] = 0.0f;
        _S887[int(11)] = 0.0f;
        _S887[int(12)] = 0.0f;
        _S887[int(13)] = 0.0f;
        _S887[int(14)] = 0.0f;
        _S887[int(15)] = 0.0f;
        _S887[int(16)] = 0.0f;
        _S887[int(17)] = 0.0f;
        _S887[int(18)] = 0.0f;
        _S887[int(19)] = 0.0f;
        _S887[int(20)] = 0.0f;
        _S887[int(21)] = 0.0f;
        _S887[int(22)] = 0.0f;
        _S887[int(5)] = _S876;
        _S887[int(4)] = _S877;
        _S887[int(3)] = _S884;
        _S887[int(2)] = _S885;
        _S887[int(6)] = _S872;
        float _S888 = _S843 + _S887[int(1)];
        float _S889 = _S844 + _S887[int(2)];
        float _S890 = _S845 + _S887[int(3)];
        float _S891 = _S846 + _S887[int(4)];
        float _S892 = _S847 + _S887[int(5)];
        float _S893 = _S848 + _S887[int(6)];
        float _S894 = _S849 + _S887[int(7)];
        float _S895 = _S850 + _S887[int(8)];
        float _S896 = _S851 + _S887[int(9)];
        float _S897 = _S852 + _S887[int(10)];
        float _S898 = _S853 + _S887[int(11)];
        float _S899 = _S854 + _S887[int(12)];
        float _S900 = _S855 + _S887[int(13)];
        float _S901 = _S856 + _S887[int(14)];
        float _S902 = _S857 + _S887[int(15)];
        float _S903 = _S858 + _S887[int(16)];
        float _S904 = _S859 + _S887[int(17)];
        float _S905 = _S860 + _S887[int(18)];
        float _S906 = _S861 + _S887[int(19)];
        float _S907 = _S862 + _S887[int(20)];
        float _S908 = _S863 + _S887[int(21)];
        float _S909 = _S864 + _S887[int(22)];
        _S865[int(0)] = _S842 + _S887[int(0)];
        _S865[int(1)] = _S888;
        _S865[int(2)] = _S889;
        _S865[int(3)] = _S890;
        _S865[int(4)] = _S891;
        _S865[int(5)] = _S892;
        _S865[int(6)] = _S893;
        _S865[int(7)] = _S894;
        _S865[int(8)] = _S895;
        _S865[int(9)] = _S896;
        _S865[int(10)] = _S897;
        _S865[int(11)] = _S898;
        _S865[int(12)] = _S899;
        _S865[int(13)] = _S900;
        _S865[int(14)] = _S901;
        _S865[int(15)] = _S902;
        _S865[int(16)] = _S903;
        _S865[int(17)] = _S904;
        _S865[int(18)] = _S905;
        _S865[int(19)] = _S906;
        _S865[int(20)] = _S907;
        _S865[int(21)] = _S908;
        _S865[int(22)] = _S909;
        _S770 = _S886;
    }
    else
    {
        _S865[int(0)] = _S842;
        _S865[int(1)] = _S843;
        _S865[int(2)] = _S844;
        _S865[int(3)] = _S845;
        _S865[int(4)] = _S846;
        _S865[int(5)] = _S847;
        _S865[int(6)] = _S848;
        _S865[int(7)] = _S849;
        _S865[int(8)] = _S850;
        _S865[int(9)] = _S851;
        _S865[int(10)] = _S852;
        _S865[int(11)] = _S853;
        _S865[int(12)] = _S854;
        _S865[int(13)] = _S855;
        _S865[int(14)] = _S856;
        _S865[int(15)] = _S857;
        _S865[int(16)] = _S858;
        _S865[int(17)] = _S859;
        _S865[int(18)] = _S860;
        _S865[int(19)] = _S861;
        _S865[int(20)] = _S862;
        _S865[int(21)] = _S863;
        _S865[int(22)] = _S864;
        _S770 = 0.0f;
    }
    if(_S768)
    {
        FixedArray<float, 23>  _S910;
        _S910[int(0)] = 0.0f;
        _S910[int(1)] = 0.0f;
        _S910[int(2)] = 0.0f;
        _S910[int(3)] = 0.0f;
        _S910[int(4)] = 0.0f;
        _S910[int(5)] = 0.0f;
        _S910[int(6)] = 0.0f;
        _S910[int(7)] = 0.0f;
        _S910[int(8)] = 0.0f;
        _S910[int(9)] = 0.0f;
        _S910[int(10)] = 0.0f;
        _S910[int(11)] = 0.0f;
        _S910[int(12)] = 0.0f;
        _S910[int(13)] = 0.0f;
        _S910[int(14)] = 0.0f;
        _S910[int(15)] = 0.0f;
        _S910[int(16)] = 0.0f;
        _S910[int(17)] = 0.0f;
        _S910[int(18)] = 0.0f;
        _S910[int(19)] = 0.0f;
        _S910[int(20)] = 0.0f;
        _S910[int(21)] = 0.0f;
        _S910[int(22)] = 0.0f;
        _S910[int(3)] = 0.0f;
        float _S911 = _S865[int(1)] + _S910[int(1)];
        float _S912 = _S865[int(2)] + _S910[int(2)];
        float _S913 = _S865[int(3)] + _S910[int(3)];
        float _S914 = _S865[int(4)] + _S910[int(4)];
        float _S915 = _S865[int(5)] + _S910[int(5)];
        float _S916 = _S865[int(6)] + _S910[int(6)];
        float _S917 = _S865[int(7)] + _S910[int(7)];
        float _S918 = _S865[int(8)] + _S910[int(8)];
        float _S919 = _S865[int(9)] + _S910[int(9)];
        float _S920 = _S865[int(10)] + _S910[int(10)];
        float _S921 = _S865[int(11)] + _S910[int(11)];
        float _S922 = _S865[int(12)] + _S910[int(12)];
        float _S923 = _S865[int(13)] + _S910[int(13)];
        float _S924 = _S865[int(14)] + _S910[int(14)];
        float _S925 = _S865[int(15)] + _S910[int(15)];
        float _S926 = _S865[int(16)] + _S910[int(16)];
        float _S927 = _S865[int(17)] + _S910[int(17)];
        float _S928 = _S865[int(18)] + _S910[int(18)];
        float _S929 = _S865[int(19)] + _S910[int(19)];
        float _S930 = _S865[int(20)] + _S910[int(20)];
        float _S931 = _S865[int(21)] + _S910[int(21)];
        float _S932 = _S865[int(22)] + _S910[int(22)];
        _S865[int(0)] = _S865[int(0)] + _S910[int(0)];
        _S865[int(1)] = _S911;
        _S865[int(2)] = _S912;
        _S865[int(3)] = _S913;
        _S865[int(4)] = _S914;
        _S865[int(5)] = _S915;
        _S865[int(6)] = _S916;
        _S865[int(7)] = _S917;
        _S865[int(8)] = _S918;
        _S865[int(9)] = _S919;
        _S865[int(10)] = _S920;
        _S865[int(11)] = _S921;
        _S865[int(12)] = _S922;
        _S865[int(13)] = _S923;
        _S865[int(14)] = _S924;
        _S865[int(15)] = _S925;
        _S865[int(16)] = _S926;
        _S865[int(17)] = _S927;
        _S865[int(18)] = _S928;
        _S865[int(19)] = _S929;
        _S865[int(20)] = _S930;
        _S865[int(21)] = _S931;
        _S865[int(22)] = _S932;
    }
    float _S933 = -10.0f * (*_s_dOut_5)[int(1)];
    DiffPair_float_0 _S934;
    (&_S934)->primal_0 = _S767;
    (&_S934)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S934, _S933);
    float _S935 = _S934.differential_0 / _S766;
    float _S936 = _S765 * _S935;
    float _S937 = (*_s_dOut_5)[int(0)] / _S766;
    float _S938 = _S765 * _S937;
    float _S939 = _S764[int(1)] * - _S935 + _S764[int(0)] * - _S937;
    DiffPair_float_0 _S940;
    (&_S940)->primal_0 = _S764[int(17)];
    (&_S940)->differential_0 = 0.0f;
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = 1.0f;
    (&_S941)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S940, &_S941, _S939);
    FixedArray<float, 23>  _S942;
    _S942[int(0)] = 0.0f;
    _S942[int(1)] = 0.0f;
    _S942[int(2)] = 0.0f;
    _S942[int(3)] = 0.0f;
    _S942[int(4)] = 0.0f;
    _S942[int(5)] = 0.0f;
    _S942[int(6)] = 0.0f;
    _S942[int(7)] = 0.0f;
    _S942[int(8)] = 0.0f;
    _S942[int(9)] = 0.0f;
    _S942[int(10)] = 0.0f;
    _S942[int(11)] = 0.0f;
    _S942[int(12)] = 0.0f;
    _S942[int(13)] = 0.0f;
    _S942[int(14)] = 0.0f;
    _S942[int(15)] = 0.0f;
    _S942[int(16)] = 0.0f;
    _S942[int(17)] = 0.0f;
    _S942[int(18)] = 0.0f;
    _S942[int(19)] = 0.0f;
    _S942[int(20)] = 0.0f;
    _S942[int(21)] = 0.0f;
    _S942[int(22)] = 0.0f;
    _S942[int(18)] = _S770;
    _S942[int(1)] = _S936;
    _S942[int(17)] = _S940.differential_0;
    _S942[int(0)] = _S938;
    FixedArray<float, 23>  _S943 = {
        _S865[int(0)] + _S942[int(0)], _S865[int(1)] + _S942[int(1)], _S865[int(2)] + _S942[int(2)], _S865[int(3)] + _S942[int(3)], _S865[int(4)] + _S942[int(4)], _S865[int(5)] + _S942[int(5)], _S865[int(6)] + _S942[int(6)], _S865[int(7)] + _S942[int(7)], _S865[int(8)] + _S942[int(8)], _S865[int(9)] + _S942[int(9)], _S865[int(10)] + _S942[int(10)], _S865[int(11)] + _S942[int(11)], _S865[int(12)] + _S942[int(12)], _S865[int(13)] + _S942[int(13)], _S865[int(14)] + _S942[int(14)], _S865[int(15)] + _S942[int(15)], _S865[int(16)] + _S942[int(16)], _S865[int(17)] + _S942[int(17)], _S865[int(18)] + _S942[int(18)], _S865[int(19)] + _S942[int(19)], _S865[int(20)] + _S942[int(20)], _S865[int(21)] + _S942[int(21)], _S865[int(22)] + _S942[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S943;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S944, FixedArray<float, 11>  * _S945, FixedArray<float, 10>  * _S946)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S944, _S945, _S946);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  * raw_losses_1, FixedArray<float, 11>  * weights_5, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 23>  * _S947)
{
    FixedArray<float, 23>  _S948 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S948;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_5, v_losses_1);
    *_S947 = (&dp_raw_losses_0)->differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_27, float3  y_9)
{
    float3  result_16;
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
        *_slang_vector_get_element_ptr(&result_16, i_13) = (F32_min((_slang_vector_get_element(x_27, i_13)), (_slang_vector_get_element(y_9, i_13))));
        i_13 = i_13 + int(1);
    }
    return result_16;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_18, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_17)
{
    DiffPair_float_0 left_dp_3;
    (&left_dp_3)->primal_0 = (*dpx_18).primal_0.x;
    (&left_dp_3)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_5).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_3;
    (&right_dp_3)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_3)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_3, &middle_dp_0, &right_dp_3, dOut_17.x);
    float3  left_d_result_2;
    *&((&left_d_result_2)->x) = left_dp_3.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_2;
    *&((&right_d_result_2)->x) = right_dp_3.differential_0;
    DiffPair_float_0 left_dp_4;
    (&left_dp_4)->primal_0 = (*dpx_18).primal_0.y;
    (&left_dp_4)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_4;
    (&right_dp_4)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_4)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_4, &middle_dp_1, &right_dp_4, dOut_17.y);
    *&((&left_d_result_2)->y) = left_dp_4.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_2)->y) = right_dp_4.differential_0;
    DiffPair_float_0 left_dp_5;
    (&left_dp_5)->primal_0 = (*dpx_18).primal_0.z;
    (&left_dp_5)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_5;
    (&right_dp_5)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_5)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_5, &middle_dp_2, &right_dp_5, dOut_17.z);
    *&((&left_d_result_2)->z) = left_dp_5.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_2)->z) = right_dp_5.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = left_d_result_2;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_2;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S949, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S950, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S951, float3  _S952)
{
    _d_clamp_vector_0(_S949, _S950, _S951, _S952);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S953 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S954 = make_float3 (_S953);
    float3  _S955 = make_float3 (0.0f);
    float3  _S956 = make_float3 (1.0f);
    float3  _S957 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S958;
    (&_S958)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S953) * (*dpbackground_0).primal_0;
    (&_S958)->differential_0 = _S957;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S959;
    (&_S959)->primal_0 = _S955;
    (&_S959)->differential_0 = _S957;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S960;
    (&_S960)->primal_0 = _S956;
    (&_S960)->differential_0 = _S957;
    s_bwd_prop_clamp_1(&_S958, &_S959, &_S960, _s_dOut_6);
    float3  _S961 = _S954 * _S958.differential_0;
    float3  _S962 = (*dpbackground_0).primal_0 * _S958.differential_0;
    float _S963 = - (_S962.x + _S962.y + _S962.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S961;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S963;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S958.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S964, DiffPair_float_0 * _S965, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S966, float3  _S967)
{
    s_bwd_prop_blend_background_0(_S964, _S965, _S966, _S967);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S968 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S968;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S968;
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
    float3  _S969 = make_float3 (1.0f - t_4);
    float3  _S970 = make_float3 (0.0f);
    float3  _S971 = make_float3 (t_4) * _s_dOut_7;
    float3  _S972 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S973;
    (&_S973)->primal_0 = s_primal_ctx_max_1((*dprgb_1).primal_0, _S970) + make_float3 (0.00392156885936856f);
    (&_S973)->differential_0 = _S972;
    s_bwd_prop_log_1(&_S973, _S971);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S974;
    (&_S974)->primal_0 = (*dprgb_1).primal_0;
    (&_S974)->differential_0 = _S972;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S975;
    (&_S975)->primal_0 = _S970;
    (&_S975)->differential_0 = _S972;
    s_bwd_prop_max_1(&_S974, &_S975, _S973.differential_0);
    float3  _S976 = _S974.differential_0 + _S969 * _s_dOut_7;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S976;
    return;
}

inline __device__ void s_bwd_log_map_image_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S977, float _S978, float3  _S979)
{
    s_bwd_prop_log_map_image_0(_S977, _S978, _S979);
    return;
}

inline __device__ float3  log_map_image_bwd(float3  rgb_3, float t_5, float3  v_out_rgb_1)
{
    float3  _S980 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S980;
    s_bwd_log_map_image_0(&p_rgb_1, t_5, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ void depth_to_normal(uint width_0, uint height_0, float2  pix_center_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_0, float4  depths_0, float3  * normal_0)
{
    FixedArray<float3 , 4>  points_0;
    float2  _S981 = float2 {intrins_0.z, intrins_0.w};
    float2  _S982 = float2 {intrins_0.x, intrins_0.y};
    float2  _S983 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S981) / _S982;
    float2  uv_10 = _S983;
    bool _S984 = undistort_point_0(_S983, dist_coeffs_10, int(12), &uv_10);
    if(!_S984)
    {
        *normal_0 = make_float3 (0.0f);
        return;
    }
    float3  raydir_5;
    if(is_fisheye_7)
    {
        float theta_4 = length_1(uv_10);
        float3  raydir_6 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_10 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_0)
        {
            raydir_5 = raydir_6 / make_float3 (raydir_6.z);
        }
        else
        {
            raydir_5 = raydir_6;
        }
    }
    else
    {
        float3  raydir_7 = make_float3 (uv_10.x, uv_10.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_5 = normalize_0(raydir_7);
        }
        else
        {
            raydir_5 = raydir_7;
        }
    }
    points_0[int(0)] = make_float3 (depths_0.x) * raydir_5;
    float2  _S985 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S981) / _S982;
    float2  uv_11 = _S985;
    bool _S986 = undistort_point_0(_S985, dist_coeffs_10, int(12), &uv_11);
    if(!_S986)
    {
        *normal_0 = make_float3 (0.0f);
        return;
    }
    if(is_fisheye_7)
    {
        float theta_5 = length_1(uv_11);
        float3  raydir_8 = make_float3 ((uv_11 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_11 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        if(!is_ray_depth_0)
        {
            raydir_5 = raydir_8 / make_float3 (raydir_8.z);
        }
        else
        {
            raydir_5 = raydir_8;
        }
    }
    else
    {
        float3  raydir_9 = make_float3 (uv_11.x, uv_11.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_5 = normalize_0(raydir_9);
        }
        else
        {
            raydir_5 = raydir_9;
        }
    }
    points_0[int(1)] = make_float3 (depths_0.y) * raydir_5;
    float2  _S987 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S981) / _S982;
    float2  uv_12 = _S987;
    bool _S988 = undistort_point_0(_S987, dist_coeffs_10, int(12), &uv_12);
    if(!_S988)
    {
        *normal_0 = make_float3 (0.0f);
        return;
    }
    if(is_fisheye_7)
    {
        float theta_6 = length_1(uv_12);
        float3  raydir_10 = make_float3 ((uv_12 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (uv_12 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
        if(!is_ray_depth_0)
        {
            raydir_5 = raydir_10 / make_float3 (raydir_10.z);
        }
        else
        {
            raydir_5 = raydir_10;
        }
    }
    else
    {
        float3  raydir_11 = make_float3 (uv_12.x, uv_12.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_5 = normalize_0(raydir_11);
        }
        else
        {
            raydir_5 = raydir_11;
        }
    }
    points_0[int(2)] = make_float3 (depths_0.z) * raydir_5;
    float2  _S989 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S981) / _S982;
    float2  uv_13 = _S989;
    bool _S990 = undistort_point_0(_S989, dist_coeffs_10, int(12), &uv_13);
    if(!_S990)
    {
        *normal_0 = make_float3 (0.0f);
        return;
    }
    if(is_fisheye_7)
    {
        float theta_7 = length_1(uv_13);
        float3  raydir_12 = make_float3 ((uv_13 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (uv_13 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        if(!is_ray_depth_0)
        {
            raydir_5 = raydir_12 / make_float3 (raydir_12.z);
        }
        else
        {
            raydir_5 = raydir_12;
        }
    }
    else
    {
        float3  raydir_13 = make_float3 (uv_13.x, uv_13.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_5 = normalize_0(raydir_13);
        }
        else
        {
            raydir_5 = raydir_13;
        }
    }
    points_0[int(3)] = make_float3 (depths_0.w) * raydir_5;
    float3  _S991 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    *normal_0 = _S991;
    if((length_2(_S991)) != 0.0f)
    {
        *normal_0 = *normal_0 / make_float3 (length_2(*normal_0));
    }
    return;
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1000, float3  _S1001)
{
    return cross_0(_S1000, _S1001);
}

inline __device__ void s_primal_ctx_depth_to_normal_0(uint width_1, uint height_1, float2  pix_center_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_1, float4  dpdepths_0, float3  * dpnormal_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S1002 = make_float2 (0.0f);
    _s_diff_ctx_2->_S992 = _S1002;
    _s_diff_ctx_2->_S993 = false;
    _s_diff_ctx_2->_S994 = _S1002;
    _s_diff_ctx_2->_S995 = false;
    _s_diff_ctx_2->_S996 = _S1002;
    _s_diff_ctx_2->_S997 = false;
    _s_diff_ctx_2->_S998 = _S1002;
    _s_diff_ctx_2->_S999 = false;
    _s_diff_ctx_2->_S994 = _S1002;
    _s_diff_ctx_2->_S995 = false;
    _s_diff_ctx_2->_S996 = _S1002;
    _s_diff_ctx_2->_S997 = false;
    _s_diff_ctx_2->_S998 = _S1002;
    _s_diff_ctx_2->_S999 = false;
    float3  _S1003 = make_float3 (0.0f);
    float2  _S1004 = float2 {intrins_1.z, intrins_1.w};
    float2  _S1005 = float2 {intrins_1.x, intrins_1.y};
    float2  _S1006 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1004) / _S1005;
    float2  _S1007 = _S1006;
    bool _S1008 = undistort_point_0(_S1006, dist_coeffs_11, int(12), &_S1007);
    _s_diff_ctx_2->_S992 = _S1007;
    _s_diff_ctx_2->_S993 = _S1008;
    float2  uv_14 = _S1007;
    bool _S1009 = !_S1008;
    float3  _S1010;
    if(_S1009)
    {
        _S1010 = make_float3 (0.0f);
    }
    else
    {
        _S1010 = _S1003;
    }
    bool _S1011 = !_S1009;
    int _S1012;
    FixedArray<float3 , 4>  points_1;
    if(_S1011)
    {
        float3  raydir_14;
        if(is_fisheye_8)
        {
            float _S1013 = length_1(uv_14);
            float3  raydir_15 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S1013, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1013))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S1013, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1013))).y, s_primal_ctx_cos_0(_S1013));
            if(!is_ray_depth_1)
            {
                raydir_14 = raydir_15 / make_float3 (raydir_15.z);
            }
            else
            {
                raydir_14 = raydir_15;
            }
        }
        else
        {
            float3  raydir_16 = make_float3 (uv_14.x, uv_14.y, 1.0f);
            if(is_ray_depth_1)
            {
                raydir_14 = normalize_0(raydir_16);
            }
            else
            {
                raydir_14 = raydir_16;
            }
        }
        float3  _S1014 = make_float3 (dpdepths_0.x) * raydir_14;
        float2  _S1015 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1004) / _S1005;
        float2  _S1016 = _S1015;
        bool _S1017 = undistort_point_0(_S1015, dist_coeffs_11, int(12), &_S1016);
        _s_diff_ctx_2->_S994 = _S1016;
        _s_diff_ctx_2->_S995 = _S1017;
        float2  uv_15 = _S1016;
        bool _S1018 = !_S1017;
        if(_S1018)
        {
            _S1010 = make_float3 (0.0f);
        }
        bool _S1019 = !_S1018;
        if(_S1019)
        {
            if(is_fisheye_8)
            {
                float _S1020 = length_1(uv_15);
                float3  raydir_17 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S1020, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1020))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S1020, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1020))).y, s_primal_ctx_cos_0(_S1020));
                if(!is_ray_depth_1)
                {
                    raydir_14 = raydir_17 / make_float3 (raydir_17.z);
                }
                else
                {
                    raydir_14 = raydir_17;
                }
            }
            else
            {
                float3  raydir_18 = make_float3 (uv_15.x, uv_15.y, 1.0f);
                if(is_ray_depth_1)
                {
                    raydir_14 = normalize_0(raydir_18);
                }
                else
                {
                    raydir_14 = raydir_18;
                }
            }
            float3  _S1021 = make_float3 (dpdepths_0.y) * raydir_14;
            _S1012 = int(2);
            points_1[int(0)] = _S1014;
            points_1[int(1)] = _S1021;
            points_1[int(2)] = _S1003;
            points_1[int(3)] = _S1003;
        }
        else
        {
            _S1012 = int(0);
            points_1[int(0)] = _S1014;
            points_1[int(1)] = _S1003;
            points_1[int(2)] = _S1003;
            points_1[int(3)] = _S1003;
        }
        bool _runFlag_0;
        if(_S1012 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1011;
            _S1012 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1022 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1004) / _S1005;
            float2  _S1023 = _S1022;
            bool _S1024 = undistort_point_0(_S1022, dist_coeffs_11, int(12), &_S1023);
            _s_diff_ctx_2->_S996 = _S1023;
            _s_diff_ctx_2->_S997 = _S1024;
            float2  uv_16 = _S1023;
            if(!_S1024)
            {
                float3  _S1025 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1012 = int(0);
                _S1010 = _S1025;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_8)
                {
                    float _S1026 = length_1(uv_16);
                    float3  raydir_19 = make_float3 ((uv_16 / make_float2 (s_primal_ctx_max_0(_S1026, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1026))).x, (uv_16 / make_float2 (s_primal_ctx_max_0(_S1026, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1026))).y, s_primal_ctx_cos_0(_S1026));
                    if(!is_ray_depth_1)
                    {
                        raydir_14 = raydir_19 / make_float3 (raydir_19.z);
                    }
                    else
                    {
                        raydir_14 = raydir_19;
                    }
                }
                else
                {
                    float3  raydir_20 = make_float3 (uv_16.x, uv_16.y, 1.0f);
                    if(is_ray_depth_1)
                    {
                        raydir_14 = normalize_0(raydir_20);
                    }
                    else
                    {
                        raydir_14 = raydir_20;
                    }
                }
                points_1[int(2)] = make_float3 (dpdepths_0.z) * raydir_14;
                float2  _S1027 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1004) / _S1005;
                float2  _S1028 = _S1027;
                bool _S1029 = undistort_point_0(_S1027, dist_coeffs_11, int(12), &_S1028);
                _s_diff_ctx_2->_S998 = _S1028;
                _s_diff_ctx_2->_S999 = _S1029;
                float2  uv_17 = _S1028;
                bool _S1030 = !_S1029;
                if(_S1030)
                {
                    _S1010 = make_float3 (0.0f);
                }
                bool _S1031 = !_S1030;
                int _S1032;
                if(_S1031)
                {
                    if(is_fisheye_8)
                    {
                        float _S1033 = length_1(uv_17);
                        float3  raydir_21 = make_float3 ((uv_17 / make_float2 (s_primal_ctx_max_0(_S1033, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1033))).x, (uv_17 / make_float2 (s_primal_ctx_max_0(_S1033, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1033))).y, s_primal_ctx_cos_0(_S1033));
                        if(!is_ray_depth_1)
                        {
                            raydir_14 = raydir_21 / make_float3 (raydir_21.z);
                        }
                        else
                        {
                            raydir_14 = raydir_21;
                        }
                    }
                    else
                    {
                        float3  raydir_22 = make_float3 (uv_17.x, uv_17.y, 1.0f);
                        if(is_ray_depth_1)
                        {
                            raydir_14 = normalize_0(raydir_22);
                        }
                        else
                        {
                            raydir_14 = raydir_22;
                        }
                    }
                    points_1[int(3)] = make_float3 (dpdepths_0.w) * raydir_14;
                    _S1032 = int(2);
                }
                else
                {
                    _S1032 = int(0);
                }
                if(_S1032 != int(2))
                {
                    _runFlag_0 = false;
                    _S1012 = _S1032;
                }
                if(_runFlag_0)
                {
                    _S1012 = int(1);
                }
            }
        }
    }
    else
    {
        _S1012 = int(0);
        points_1[int(0)] = _S1003;
        points_1[int(1)] = _S1003;
        points_1[int(2)] = _S1003;
        points_1[int(3)] = _S1003;
    }
    if(!(_S1012 != int(1)))
    {
        float3  _S1034 = s_primal_ctx_cross_0(points_1[int(1)] - points_1[int(0)], - (points_1[int(3)] - points_1[int(2)]));
        float _S1035 = length_2(_S1034);
        if(_S1035 != 0.0f)
        {
            _S1010 = _S1034 / make_float3 (_S1035);
        }
        else
        {
            _S1010 = _S1034;
        }
    }
    *dpnormal_0 = _S1010;
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1036, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1037, float3  _S1038)
{
    _d_cross_0(_S1036, _S1037, _S1038);
    return;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(uint width_2, uint height_2, float2  pix_center_2, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9, bool is_ray_depth_2, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  dpnormal_1, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1039 = *dpdepths_1;
    float3  _S1040 = make_float3 (0.0f);
    bool _S1041 = !!_s_diff_ctx_3->_S993;
    float3  raydir_23;
    float3  raydir_24;
    float3  raydir_25;
    float3  raydir_26;
    int _S1042;
    FixedArray<float3 , 4>  points_2;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1043;
    if(_S1041)
    {
        if(is_fisheye_9)
        {
            float _S1044 = length_1(_s_diff_ctx_3->_S992);
            float3  raydir_27 = make_float3 ((_s_diff_ctx_3->_S992 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).x, (_s_diff_ctx_3->_S992 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).y, s_primal_ctx_cos_0(_S1044));
            if(!is_ray_depth_2)
            {
                raydir_23 = raydir_27 / make_float3 (raydir_27.z);
            }
            else
            {
                raydir_23 = raydir_27;
            }
        }
        else
        {
            float3  raydir_28 = make_float3 (_s_diff_ctx_3->_S992.x, _s_diff_ctx_3->_S992.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_23 = normalize_0(raydir_28);
            }
            else
            {
                raydir_23 = raydir_28;
            }
        }
        float3  _S1045 = make_float3 (_S1039.primal_0.x) * raydir_23;
        bool _S1046 = !!_s_diff_ctx_3->_S995;
        if(_S1046)
        {
            if(is_fisheye_9)
            {
                float _S1047 = length_1(_s_diff_ctx_3->_S994);
                float3  raydir_29 = make_float3 ((_s_diff_ctx_3->_S994 / make_float2 (s_primal_ctx_max_0(_S1047, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1047))).x, (_s_diff_ctx_3->_S994 / make_float2 (s_primal_ctx_max_0(_S1047, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1047))).y, s_primal_ctx_cos_0(_S1047));
                if(!is_ray_depth_2)
                {
                    raydir_24 = raydir_29 / make_float3 (raydir_29.z);
                }
                else
                {
                    raydir_24 = raydir_29;
                }
            }
            else
            {
                float3  raydir_30 = make_float3 (_s_diff_ctx_3->_S994.x, _s_diff_ctx_3->_S994.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_24 = normalize_0(raydir_30);
                }
                else
                {
                    raydir_24 = raydir_30;
                }
            }
            float3  _S1048 = make_float3 (_S1039.primal_0.y) * raydir_24;
            _S1042 = int(2);
            points_2[int(0)] = _S1045;
            points_2[int(1)] = _S1048;
            points_2[int(2)] = _S1040;
            points_2[int(3)] = _S1040;
        }
        else
        {
            _S1042 = int(0);
            points_2[int(0)] = _S1045;
            points_2[int(1)] = _S1040;
            points_2[int(2)] = _S1040;
            points_2[int(3)] = _S1040;
            raydir_24 = _S1040;
        }
        if(_S1042 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1041;
            _S1042 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_3->_S997)
            {
                _runFlag_2 = false;
                _S1042 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_9)
                {
                    float _S1049 = length_1(_s_diff_ctx_3->_S996);
                    float3  raydir_31 = make_float3 ((_s_diff_ctx_3->_S996 / make_float2 (s_primal_ctx_max_0(_S1049, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1049))).x, (_s_diff_ctx_3->_S996 / make_float2 (s_primal_ctx_max_0(_S1049, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1049))).y, s_primal_ctx_cos_0(_S1049));
                    if(!is_ray_depth_2)
                    {
                        raydir_25 = raydir_31 / make_float3 (raydir_31.z);
                    }
                    else
                    {
                        raydir_25 = raydir_31;
                    }
                }
                else
                {
                    float3  raydir_32 = make_float3 (_s_diff_ctx_3->_S996.x, _s_diff_ctx_3->_S996.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_25 = normalize_0(raydir_32);
                    }
                    else
                    {
                        raydir_25 = raydir_32;
                    }
                }
                points_2[int(2)] = make_float3 (_S1039.primal_0.z) * raydir_25;
                bool _S1050 = !!_s_diff_ctx_3->_S999;
                int _S1051;
                if(_S1050)
                {
                    if(is_fisheye_9)
                    {
                        float _S1052 = length_1(_s_diff_ctx_3->_S998);
                        float3  raydir_33 = make_float3 ((_s_diff_ctx_3->_S998 / make_float2 (s_primal_ctx_max_0(_S1052, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1052))).x, (_s_diff_ctx_3->_S998 / make_float2 (s_primal_ctx_max_0(_S1052, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1052))).y, s_primal_ctx_cos_0(_S1052));
                        if(!is_ray_depth_2)
                        {
                            raydir_26 = raydir_33 / make_float3 (raydir_33.z);
                        }
                        else
                        {
                            raydir_26 = raydir_33;
                        }
                    }
                    else
                    {
                        float3  raydir_34 = make_float3 (_s_diff_ctx_3->_S998.x, _s_diff_ctx_3->_S998.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_26 = normalize_0(raydir_34);
                        }
                        else
                        {
                            raydir_26 = raydir_34;
                        }
                    }
                    points_2[int(3)] = make_float3 (_S1039.primal_0.w) * raydir_26;
                    _S1051 = int(2);
                }
                else
                {
                    _S1051 = int(0);
                    raydir_26 = _S1040;
                }
                if(_S1051 != int(2))
                {
                    _runFlag_3 = false;
                    _S1042 = _S1051;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1042 = int(1);
                }
                float3  _S1053 = raydir_25;
                _runFlag_3 = _S1050;
                raydir_25 = raydir_26;
                raydir_26 = _S1053;
            }
            else
            {
                _runFlag_3 = false;
                raydir_25 = _S1040;
                raydir_26 = _S1040;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_25 = _S1040;
            raydir_26 = _S1040;
        }
        float3  _S1054 = raydir_23;
        float3  _S1055 = raydir_24;
        raydir_23 = raydir_25;
        raydir_24 = raydir_26;
        _S1043 = _S1046;
        raydir_25 = _S1055;
        raydir_26 = _S1054;
    }
    else
    {
        _S1042 = int(0);
        points_2[int(0)] = _S1040;
        points_2[int(1)] = _S1040;
        points_2[int(2)] = _S1040;
        points_2[int(3)] = _S1040;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_23 = _S1040;
        raydir_24 = _S1040;
        _S1043 = false;
        raydir_25 = _S1040;
        raydir_26 = _S1040;
    }
    bool _S1056 = !(_S1042 != int(1));
    float3  _S1057;
    float3  _S1058;
    float3  _S1059;
    float3  _S1060;
    float3  _S1061;
    bool _S1062;
    if(_S1056)
    {
        float3  dx_0 = points_2[int(1)] - points_2[int(0)];
        float3  _S1063 = - (points_2[int(3)] - points_2[int(2)]);
        float3  _S1064 = s_primal_ctx_cross_0(dx_0, _S1063);
        float _S1065 = length_2(_S1064);
        float3  _S1066 = make_float3 (_S1065);
        bool _S1067 = _S1065 != 0.0f;
        if(_S1067)
        {
            _S1057 = make_float3 (_S1065 * _S1065);
        }
        else
        {
            _S1057 = _S1040;
        }
        _S1062 = _S1067;
        _S1058 = _S1064;
        _S1059 = _S1066;
        _S1060 = dx_0;
        _S1061 = _S1063;
    }
    else
    {
        _S1062 = false;
        _S1057 = _S1040;
        _S1058 = _S1040;
        _S1059 = _S1040;
        _S1060 = _S1040;
        _S1061 = _S1040;
    }
    float4  _S1068 = make_float4 (0.0f);
    if(_S1056)
    {
        if(_S1062)
        {
            float3  _S1069 = dpnormal_1 / _S1057;
            float3  _S1070 = _S1059 * _S1069;
            _S1057 = _S1058 * - _S1069;
            _S1059 = _S1070;
        }
        else
        {
            _S1057 = _S1040;
            _S1059 = dpnormal_1;
        }
        float _S1071 = _S1057.x + _S1057.y + _S1057.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1072;
        (&_S1072)->primal_0 = _S1058;
        (&_S1072)->differential_0 = _S1040;
        s_bwd_length_impl_1(&_S1072, _S1071);
        float3  _S1073 = _S1072.differential_0 + _S1059;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1074;
        (&_S1074)->primal_0 = _S1060;
        (&_S1074)->differential_0 = _S1040;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1075;
        (&_S1075)->primal_0 = _S1061;
        (&_S1075)->differential_0 = _S1040;
        s_bwd_prop_cross_0(&_S1074, &_S1075, _S1073);
        float3  s_diff_dy_T_0 = - _S1075.differential_0;
        float3  _S1076 = - s_diff_dy_T_0;
        float3  _S1077 = - _S1074.differential_0;
        FixedArray<float3 , 4>  _S1078;
        _S1078[int(0)] = _S1040;
        _S1078[int(1)] = _S1040;
        _S1078[int(2)] = _S1040;
        _S1078[int(3)] = _S1040;
        _S1078[int(2)] = _S1076;
        _S1078[int(3)] = s_diff_dy_T_0;
        _S1078[int(0)] = _S1077;
        _S1078[int(1)] = _S1074.differential_0;
        points_2[int(0)] = _S1078[int(0)];
        points_2[int(1)] = _S1078[int(1)];
        points_2[int(2)] = _S1078[int(2)];
        points_2[int(3)] = _S1078[int(3)];
    }
    else
    {
        points_2[int(0)] = _S1040;
        points_2[int(1)] = _S1040;
        points_2[int(2)] = _S1040;
        points_2[int(3)] = _S1040;
    }
    float4  _S1079;
    if(_S1041)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1080 = points_2;
                FixedArray<float3 , 4>  _S1081 = points_2;
                FixedArray<float3 , 4>  _S1082 = points_2;
                FixedArray<float3 , 4>  _S1083 = points_2;
                if(_runFlag_3)
                {
                    float3  _S1084 = raydir_23 * _S1083[int(3)];
                    float _S1085 = _S1084.x + _S1084.y + _S1084.z;
                    float4  _S1086 = _S1068;
                    *&((&_S1086)->w) = _S1085;
                    points_2[int(0)] = _S1080[int(0)];
                    points_2[int(1)] = _S1081[int(1)];
                    points_2[int(2)] = _S1082[int(2)];
                    points_2[int(3)] = _S1040;
                    _S1079 = _S1086;
                }
                else
                {
                    points_2[int(0)] = _S1080[int(0)];
                    points_2[int(1)] = _S1081[int(1)];
                    points_2[int(2)] = _S1082[int(2)];
                    points_2[int(3)] = _S1083[int(3)];
                    _S1079 = _S1068;
                }
                float3  _S1087 = raydir_24 * points_2[int(2)];
                float _S1088 = _S1087.x + _S1087.y + _S1087.z;
                FixedArray<float3 , 4>  _S1089 = points_2;
                FixedArray<float3 , 4>  _S1090 = points_2;
                float4  _S1091 = _S1068;
                *&((&_S1091)->z) = _S1088;
                float4  _S1092 = _S1079 + _S1091;
                points_2[int(0)] = points_2[int(0)];
                points_2[int(1)] = _S1089[int(1)];
                points_2[int(2)] = _S1040;
                points_2[int(3)] = _S1090[int(3)];
                _S1079 = _S1092;
            }
            else
            {
                FixedArray<float3 , 4>  _S1093 = points_2;
                FixedArray<float3 , 4>  _S1094 = points_2;
                FixedArray<float3 , 4>  _S1095 = points_2;
                points_2[int(0)] = points_2[int(0)];
                points_2[int(1)] = _S1093[int(1)];
                points_2[int(2)] = _S1094[int(2)];
                points_2[int(3)] = _S1095[int(3)];
                _S1079 = _S1068;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1096 = points_2;
            FixedArray<float3 , 4>  _S1097 = points_2;
            FixedArray<float3 , 4>  _S1098 = points_2;
            points_2[int(0)] = points_2[int(0)];
            points_2[int(1)] = _S1096[int(1)];
            points_2[int(2)] = _S1097[int(2)];
            points_2[int(3)] = _S1098[int(3)];
            _S1079 = _S1068;
        }
        if(_S1043)
        {
            FixedArray<float3 , 4>  _S1099 = points_2;
            float3  _S1100 = raydir_25 * points_2[int(1)];
            float _S1101 = _S1100.x + _S1100.y + _S1100.z;
            float4  _S1102 = _S1068;
            *&((&_S1102)->y) = _S1101;
            float4  _S1103 = _S1079 + _S1102;
            points_2[int(0)] = _S1040;
            points_2[int(1)] = _S1040;
            points_2[int(2)] = _S1040;
            points_2[int(3)] = _S1040;
            raydir_23 = _S1099[int(0)];
            _S1079 = _S1103;
        }
        else
        {
            FixedArray<float3 , 4>  _S1104 = points_2;
            FixedArray<float3 , 4>  _S1105 = points_2;
            FixedArray<float3 , 4>  _S1106 = points_2;
            points_2[int(0)] = points_2[int(0)];
            points_2[int(1)] = _S1104[int(1)];
            points_2[int(2)] = _S1105[int(2)];
            points_2[int(3)] = _S1106[int(3)];
            raydir_23 = _S1040;
        }
        float3  _S1107 = raydir_26 * (points_2[int(0)] + raydir_23);
        float _S1108 = _S1107.x + _S1107.y + _S1107.z;
        float4  _S1109 = _S1068;
        *&((&_S1109)->x) = _S1108;
        _S1079 = _S1079 + _S1109;
    }
    else
    {
        _S1079 = _S1068;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1079;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(uint _S1110, uint _S1111, float2  _S1112, float4  _S1113, FixedArray<float, 10>  * _S1114, bool _S1115, bool _S1116, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1117, float3  _S1118)
{
    float3  _S1119;
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1120;
    s_primal_ctx_depth_to_normal_0(_S1110, _S1111, _S1112, _S1113, _S1114, _S1115, _S1116, (*_S1117).primal_0, &_S1119, &_S1120);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1121 = _S1120;
    s_bwd_prop_depth_to_normal_0(_S1110, _S1111, _S1112, _S1113, _S1114, _S1115, _S1116, _S1117, _S1118, &_S1121);
    return;
}

inline __device__ void depth_to_normal_vjp(uint width_3, uint height_3, float2  pix_center_3, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_13, bool is_fisheye_10, bool is_ray_depth_3, float4  depths_1, float3  v_normal_0, float4  * v_depths_0)
{
    float4  _S1122 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1122;
    s_bwd_depth_to_normal_0(width_3, height_3, pix_center_3, intrins_3, dist_coeffs_13, is_fisheye_10, is_ray_depth_3, &dp_depths_0, v_normal_0);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(uint width_4, uint height_4, float2  pix_center_4, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_14, bool is_fisheye_11)
{
    float2  _S1123 = (pix_center_4 - float2 {intrins_4.z, intrins_4.w}) / float2 {intrins_4.x, intrins_4.y};
    float2  uv_18 = _S1123;
    bool _S1124 = undistort_point_0(_S1123, dist_coeffs_14, int(12), &uv_18);
    if(!_S1124)
    {
        return 0.0f;
    }
    float3  raydir_35;
    if(is_fisheye_11)
    {
        float theta_8 = length_1(uv_18);
        float3  raydir_36 = make_float3 ((uv_18 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).x, (uv_18 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).y, (F32_cos((theta_8))));
        raydir_35 = raydir_36 / make_float3 (raydir_36.z);
    }
    else
    {
        raydir_35 = make_float3 (uv_18.x, uv_18.y, 1.0f);
    }
    return float((F32_sign((raydir_35.z)))) / length_2(raydir_35);
}

