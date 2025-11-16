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

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_10, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_9)
{
    DiffPair_float_0 _S283 = *dpx_10;
    bool _S284;
    if(((*dpx_10).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S284 = ((*dpx_10).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S284 = false;
    }
    float _S285;
    if(_S284)
    {
        _S285 = dOut_9;
    }
    else
    {
        _S285 = 0.0f;
    }
    dpx_10->primal_0 = _S283.primal_0;
    dpx_10->differential_0 = _S285;
    DiffPair_float_0 _S286 = *dpMin_0;
    if((_S283.primal_0) < ((*dpMin_0).primal_0))
    {
        _S285 = dOut_9;
    }
    else
    {
        _S285 = 0.0f;
    }
    dpMin_0->primal_0 = _S286.primal_0;
    dpMin_0->differential_0 = _S285;
    DiffPair_float_0 _S287 = *dpMax_0;
    if(((*dpx_10).primal_0) > ((*dpMax_0).primal_0))
    {
        _S285 = dOut_9;
    }
    else
    {
        _S285 = 0.0f;
    }
    dpMax_0->primal_0 = _S287.primal_0;
    dpMax_0->differential_0 = _S285;
    return;
}

inline __device__ float clamp_0(float x_13, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_13), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_11, float dOut_10)
{
    float _S288 = -0.5f / ((*dpx_11).primal_0 * (F32_sqrt(((*dpx_11).primal_0)))) * dOut_10;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S288;
    return;
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

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, FixedArray<float, 22>  * _S293)
{
    float3  _S294;
    bool _S295;
    bool _S296;
    FixedArray<float, 22>  losses_1;
    float _S297 = float(mask_0);
    float3  _S298 = ref_rgb_0 - render_rgb_0;
    float3  _S299 = abs_0(_S298);
    losses_1[int(0)] = _S297 * ((_S299.x + _S299.y + _S299.z) * 0.3333333432674408f);
    losses_1[int(1)] = _S297 * clamp_0(dot_0(_S298, _S298) * 0.3333333432674408f, 0.0f, 1.0f);
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
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S294, _S309);
    losses_1[int(7)] = _S311 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S312 = float(_S310 & _S308);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S305, _S309);
    losses_1[int(8)] = _S312 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S313 = float(_S306 & _S310);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S294, _S305);
    losses_1[int(11)] = _S313 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S314 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S315 = float(ref_alpha_0);
    float _S316 = (F32_max((_S314), (_S315)));
    losses_1[int(9)] = - lerp_0((F32_log(((F32_max((1.0f - _S316), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S316), (9.99999997475242708e-07f)))))), _S315);
    float _S317 = 1.0f - _S314;
    float _S318 = 1.0f - _S315;
    float _S319 = (F32_max((_S317), (_S318)));
    losses_1[int(10)] = - lerp_0((F32_log(((F32_max((1.0f - _S319), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S319), (9.99999997475242708e-07f)))))), _S318);
    losses_1[int(12)] = 4.0f * _S314 * _S317;
    float _S320 = (F32_max((_S314), (9.999999960041972e-13f)));
    losses_1[int(13)] = (rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f / _S320;
    losses_1[int(14)] = depth_dist_0 / _S320;
    losses_1[int(15)] = (normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f / _S320;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S297;
    losses_1[int(18)] = _S300;
    losses_1[int(19)] = _S311;
    losses_1[int(20)] = _S312;
    losses_1[int(21)] = _S313;
    *_S293 = losses_1;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S321, float3  _S322)
{
    return dot_0(_S321, _S322);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S323)
{
    return (F32_rsqrt((_S323)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S324, float _S325, float _S326)
{
    return clamp_0(_S324, _S325, _S326);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S327, DiffPair_float_0 * _S328, DiffPair_float_0 * _S329, float _S330)
{
    _d_lerp_0(_S327, _S328, _S329, _S330);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S331, DiffPair_float_0 * _S332, DiffPair_float_0 * _S333, float _S334)
{
    _d_clamp_0(_S331, _S332, _S333, _S334);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S335, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S336, float _S337)
{
    _d_dot_0(_S335, _S336, _S337);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S338, float _S339)
{
    _d_rsqrt_0(_S338, _S339);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S340, float3  _S341)
{
    _d_abs_vector_0(_S340, _S341);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, FixedArray<float, 22>  * _s_dOut_2)
{
    DiffPair_float_0 _S342 = *dprender_depth_0;
    DiffPair_float_0 _S343 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S344 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S345 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S346 = *dpref_normal_0;
    DiffPair_float_0 _S347 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S348 = *dprgb_dist_0;
    DiffPair_float_0 _S349 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S350 = *dpnormal_dist_0;
    float3  _S351 = make_float3 (0.0f);
    float _S352 = float(mask_1);
    float3  _S353 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S354 = s_primal_ctx_dot_0(_S353, _S353) * 0.3333333432674408f;
    float _S355 = float(depth_mask_1 & mask_1);
    float _S356 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S357 = _S355 * s_primal_ctx_log_0(_S356);
    float _S358 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S359 = _S355 * s_primal_ctx_log_0(_S358);
    bool _S360 = normal_mask_1 & mask_1;
    float _S361 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S362 = _S361 == 0.0f;
    float3  _S363;
    if(_S362)
    {
        _S363 = make_float3 (0.0f);
    }
    bool _S364 = !_S362;
    float3  _S365;
    if(_S364)
    {
        float _S366 = s_primal_ctx_rsqrt_0(_S361);
        float3  _S367 = make_float3 (_S366);
        _S363 = _S344.primal_0 * make_float3 (_S366);
        _S365 = _S367;
    }
    else
    {
        _S365 = _S351;
    }
    float _S368 = s_primal_ctx_dot_0(_S345.primal_0, _S345.primal_0);
    bool _S369 = _S368 == 0.0f;
    float3  _S370;
    if(_S369)
    {
        _S370 = make_float3 (0.0f);
    }
    bool _S371 = !_S369;
    float3  _S372;
    if(_S371)
    {
        float _S373 = s_primal_ctx_rsqrt_0(_S368);
        float3  _S374 = make_float3 (_S373);
        _S370 = _S345.primal_0 * make_float3 (_S373);
        _S372 = _S374;
    }
    else
    {
        _S372 = _S351;
    }
    float _S375 = s_primal_ctx_dot_0(_S346.primal_0, _S346.primal_0);
    bool _S376 = _S375 == 0.0f;
    float3  _S377;
    bool _S378;
    if(_S376)
    {
        float3  _S379 = make_float3 (0.0f);
        _S378 = false;
        _S377 = _S379;
    }
    else
    {
        _S378 = _S360;
    }
    bool _S380 = !_S376;
    float3  _S381;
    if(_S380)
    {
        float _S382 = s_primal_ctx_rsqrt_0(_S375);
        float3  _S383 = make_float3 (_S382);
        _S377 = _S346.primal_0 * make_float3 (_S382);
        _S381 = _S383;
    }
    else
    {
        _S381 = _S351;
    }
    float _S384 = float(_S364 & _S378);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S363, _S377);
    float _S385 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S386 = float(_S371 & _S378);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S370, _S377);
    float _S387 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S388 = float(_S364 & _S371);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S363, _S370);
    float _S389 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S390 = s_primal_ctx_clamp_0(_S347.primal_0, 0.0f, 1.0f);
    float _S391 = float(ref_alpha_1);
    float _S392 = s_primal_ctx_max_0(_S390, _S391);
    float _S393 = 1.0f - _S392;
    float _S394 = s_primal_ctx_max_0(_S393, 9.99999997475242708e-07f);
    float _S395 = s_primal_ctx_log_0(_S394);
    float _S396 = s_primal_ctx_max_0(_S392, 9.99999997475242708e-07f);
    float _S397 = s_primal_ctx_log_0(_S396);
    float _S398 = 1.0f - _S390;
    float _S399 = 1.0f - _S391;
    float _S400 = s_primal_ctx_max_0(_S398, _S399);
    float _S401 = 1.0f - _S400;
    float _S402 = s_primal_ctx_max_0(_S401, 9.99999997475242708e-07f);
    float _S403 = s_primal_ctx_log_0(_S402);
    float _S404 = s_primal_ctx_max_0(_S400, 9.99999997475242708e-07f);
    float _S405 = s_primal_ctx_log_0(_S404);
    float _S406 = 4.0f * _S390;
    float _S407 = s_primal_ctx_max_0(_S390, 9.999999960041972e-13f);
    float _S408 = _S407 * _S407;
    float _S409 = (*_s_dOut_2)[int(15)] / _S408;
    float _S410 = 0.3333333432674408f * (_S407 * _S409);
    float _S411 = (*_s_dOut_2)[int(14)] / _S408;
    float _S412 = _S407 * _S411;
    float _S413 = (*_s_dOut_2)[int(13)] / _S408;
    float _S414 = _S407 * _S413;
    float _S415 = (_S350.primal_0.x + _S350.primal_0.y + _S350.primal_0.z) * 0.3333333432674408f * - _S409 + _S349.primal_0 * - _S411 + (_S348.primal_0.x + _S348.primal_0.y + _S348.primal_0.z) * 0.3333333432674408f * - _S413;
    DiffPair_float_0 _S416;
    (&_S416)->primal_0 = _S390;
    (&_S416)->differential_0 = 0.0f;
    DiffPair_float_0 _S417;
    (&_S417)->primal_0 = 9.999999960041972e-13f;
    (&_S417)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S416, &_S417, _S415);
    float _S418 = 0.3333333432674408f * _S414;
    float _S419 = _S406 * (*_s_dOut_2)[int(12)];
    float _S420 = 4.0f * (_S398 * (*_s_dOut_2)[int(12)]);
    float _S421 = - (*_s_dOut_2)[int(10)];
    DiffPair_float_0 _S422;
    (&_S422)->primal_0 = _S403;
    (&_S422)->differential_0 = 0.0f;
    DiffPair_float_0 _S423;
    (&_S423)->primal_0 = _S405;
    (&_S423)->differential_0 = 0.0f;
    DiffPair_float_0 _S424;
    (&_S424)->primal_0 = _S399;
    (&_S424)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S422, &_S423, &_S424, _S421);
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = _S404;
    (&_S425)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S425, _S423.differential_0);
    DiffPair_float_0 _S426;
    (&_S426)->primal_0 = _S400;
    (&_S426)->differential_0 = 0.0f;
    DiffPair_float_0 _S427;
    (&_S427)->primal_0 = 9.99999997475242708e-07f;
    (&_S427)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S426, &_S427, _S425.differential_0);
    DiffPair_float_0 _S428;
    (&_S428)->primal_0 = _S402;
    (&_S428)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S428, _S422.differential_0);
    DiffPair_float_0 _S429;
    (&_S429)->primal_0 = _S401;
    (&_S429)->differential_0 = 0.0f;
    DiffPair_float_0 _S430;
    (&_S430)->primal_0 = 9.99999997475242708e-07f;
    (&_S430)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S429, &_S430, _S428.differential_0);
    float _S431 = _S426.differential_0 + - _S429.differential_0;
    DiffPair_float_0 _S432;
    (&_S432)->primal_0 = _S398;
    (&_S432)->differential_0 = 0.0f;
    DiffPair_float_0 _S433;
    (&_S433)->primal_0 = _S399;
    (&_S433)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S432, &_S433, _S431);
    float _S434 = - (_S419 + _S432.differential_0);
    float _S435 = - (*_s_dOut_2)[int(9)];
    DiffPair_float_0 _S436;
    (&_S436)->primal_0 = _S395;
    (&_S436)->differential_0 = 0.0f;
    DiffPair_float_0 _S437;
    (&_S437)->primal_0 = _S397;
    (&_S437)->differential_0 = 0.0f;
    DiffPair_float_0 _S438;
    (&_S438)->primal_0 = _S391;
    (&_S438)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S436, &_S437, &_S438, _S435);
    DiffPair_float_0 _S439;
    (&_S439)->primal_0 = _S396;
    (&_S439)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S439, _S437.differential_0);
    DiffPair_float_0 _S440;
    (&_S440)->primal_0 = _S392;
    (&_S440)->differential_0 = 0.0f;
    DiffPair_float_0 _S441;
    (&_S441)->primal_0 = 9.99999997475242708e-07f;
    (&_S441)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S440, &_S441, _S439.differential_0);
    DiffPair_float_0 _S442;
    (&_S442)->primal_0 = _S394;
    (&_S442)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S442, _S436.differential_0);
    DiffPair_float_0 _S443;
    (&_S443)->primal_0 = _S393;
    (&_S443)->differential_0 = 0.0f;
    DiffPair_float_0 _S444;
    (&_S444)->primal_0 = 9.99999997475242708e-07f;
    (&_S444)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S443, &_S444, _S442.differential_0);
    float _S445 = _S440.differential_0 + - _S443.differential_0;
    DiffPair_float_0 _S446;
    (&_S446)->primal_0 = _S390;
    (&_S446)->differential_0 = 0.0f;
    DiffPair_float_0 _S447;
    (&_S447)->primal_0 = _S391;
    (&_S447)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S446, &_S447, _S445);
    float _S448 = _S416.differential_0 + _S420 + _S434 + _S446.differential_0;
    DiffPair_float_0 _S449;
    (&_S449)->primal_0 = _S347.primal_0;
    (&_S449)->differential_0 = 0.0f;
    DiffPair_float_0 _S450;
    (&_S450)->primal_0 = 0.0f;
    (&_S450)->differential_0 = 0.0f;
    DiffPair_float_0 _S451;
    (&_S451)->primal_0 = 1.0f;
    (&_S451)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S449, &_S450, &_S451, _S448);
    DiffPair_float_0 _S452 = _S449;
    float _S453 = _S388 * (*_s_dOut_2)[int(11)];
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S389;
    (&_S454)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S454, _S453);
    DiffPair_float_0 _S455;
    (&_S455)->primal_0 = cos_sim_loss_5;
    (&_S455)->differential_0 = 0.0f;
    DiffPair_float_0 _S456;
    (&_S456)->primal_0 = 9.999999960041972e-13f;
    (&_S456)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S455, &_S456, _S454.differential_0);
    float _S457 = 0.5f * - (_S453 + _S455.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S458;
    (&_S458)->primal_0 = _S363;
    (&_S458)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S459;
    (&_S459)->primal_0 = _S370;
    (&_S459)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S458, &_S459, _S457);
    float _S460 = _S386 * (*_s_dOut_2)[int(8)];
    DiffPair_float_0 _S461;
    (&_S461)->primal_0 = _S387;
    (&_S461)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S461, _S460);
    DiffPair_float_0 _S462;
    (&_S462)->primal_0 = cos_sim_loss_4;
    (&_S462)->differential_0 = 0.0f;
    DiffPair_float_0 _S463;
    (&_S463)->primal_0 = 9.999999960041972e-13f;
    (&_S463)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S462, &_S463, _S461.differential_0);
    float _S464 = 0.5f * - (_S460 + _S462.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S465;
    (&_S465)->primal_0 = _S370;
    (&_S465)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S466;
    (&_S466)->primal_0 = _S377;
    (&_S466)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S465, &_S466, _S464);
    float _S467 = _S384 * (*_s_dOut_2)[int(7)];
    DiffPair_float_0 _S468;
    (&_S468)->primal_0 = _S385;
    (&_S468)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S468, _S467);
    DiffPair_float_0 _S469;
    (&_S469)->primal_0 = cos_sim_loss_3;
    (&_S469)->differential_0 = 0.0f;
    DiffPair_float_0 _S470;
    (&_S470)->primal_0 = 9.999999960041972e-13f;
    (&_S470)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S469, &_S470, _S468.differential_0);
    float _S471 = 0.5f * - (_S467 + _S469.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S472;
    (&_S472)->primal_0 = _S363;
    (&_S472)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S473;
    (&_S473)->primal_0 = _S377;
    (&_S473)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S472, &_S473, _S471);
    float3  _S474 = _S466.differential_0 + _S473.differential_0;
    float3  _S475 = _S458.differential_0 + _S472.differential_0;
    float3  _S476 = make_float3 (_S410, _S410, _S410);
    float3  _S477 = make_float3 (_S418, _S418, _S418);
    float3  _S478 = _S459.differential_0 + _S465.differential_0;
    float _S479;
    if(_S380)
    {
        float3  _S480 = _S346.primal_0 * _S474;
        float3  _S481 = _S381 * _S474;
        float _S482 = _S480.x + _S480.y + _S480.z;
        DiffPair_float_0 _S483;
        (&_S483)->primal_0 = _S375;
        (&_S483)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S483, _S482);
        _S479 = _S483.differential_0;
        _S363 = _S481;
    }
    else
    {
        _S479 = 0.0f;
        _S363 = _S351;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S346.primal_0;
    (&_S484)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S485;
    (&_S485)->primal_0 = _S346.primal_0;
    (&_S485)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S484, &_S485, _S479);
    float3  _S486 = _S485.differential_0 + _S484.differential_0 + _S363;
    if(_S371)
    {
        float3  _S487 = _S345.primal_0 * _S478;
        float3  _S488 = _S372 * _S478;
        float _S489 = _S487.x + _S487.y + _S487.z;
        DiffPair_float_0 _S490;
        (&_S490)->primal_0 = _S368;
        (&_S490)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S490, _S489);
        _S479 = _S490.differential_0;
        _S363 = _S488;
    }
    else
    {
        _S479 = 0.0f;
        _S363 = _S351;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S491;
    (&_S491)->primal_0 = _S345.primal_0;
    (&_S491)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S492;
    (&_S492)->primal_0 = _S345.primal_0;
    (&_S492)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S491, &_S492, _S479);
    float3  _S493 = _S492.differential_0 + _S491.differential_0 + _S363;
    if(_S364)
    {
        float3  _S494 = _S344.primal_0 * _S475;
        float3  _S495 = _S365 * _S475;
        float _S496 = _S494.x + _S494.y + _S494.z;
        DiffPair_float_0 _S497;
        (&_S497)->primal_0 = _S361;
        (&_S497)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S497, _S496);
        _S479 = _S497.differential_0;
        _S363 = _S495;
    }
    else
    {
        _S479 = 0.0f;
        _S363 = _S351;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S498;
    (&_S498)->primal_0 = _S344.primal_0;
    (&_S498)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S499;
    (&_S499)->primal_0 = _S344.primal_0;
    (&_S499)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S498, &_S499, _S479);
    float _S500 = _S359 * (*_s_dOut_2)[int(6)];
    float _S501 = _S359 * (*_s_dOut_2)[int(5)];
    float _S502 = _S357 * (*_s_dOut_2)[int(4)];
    float _S503 = _S355 * (_S357 * (*_s_dOut_2)[int(6)] + _S501 + _S501 + (*_s_dOut_2)[int(3)]);
    DiffPair_float_0 _S504;
    (&_S504)->primal_0 = _S358;
    (&_S504)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S504, _S503);
    DiffPair_float_0 _S505;
    (&_S505)->primal_0 = _S343.primal_0;
    (&_S505)->differential_0 = 0.0f;
    DiffPair_float_0 _S506;
    (&_S506)->primal_0 = 0.00009999999747379f;
    (&_S506)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S505, &_S506, _S504.differential_0);
    float _S507 = _S355 * (_S500 + _S502 + _S502 + (*_s_dOut_2)[int(2)]);
    DiffPair_float_0 _S508;
    (&_S508)->primal_0 = _S356;
    (&_S508)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S508, _S507);
    DiffPair_float_0 _S509;
    (&_S509)->primal_0 = _S342.primal_0;
    (&_S509)->differential_0 = 0.0f;
    DiffPair_float_0 _S510;
    (&_S510)->primal_0 = 0.00009999999747379f;
    (&_S510)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S509, &_S510, _S508.differential_0);
    float _S511 = _S352 * (*_s_dOut_2)[int(1)];
    DiffPair_float_0 _S512;
    (&_S512)->primal_0 = _S354;
    (&_S512)->differential_0 = 0.0f;
    DiffPair_float_0 _S513;
    (&_S513)->primal_0 = 0.0f;
    (&_S513)->differential_0 = 0.0f;
    DiffPair_float_0 _S514;
    (&_S514)->primal_0 = 1.0f;
    (&_S514)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S512, &_S513, &_S514, _S511);
    float _S515 = 0.3333333432674408f * _S512.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S516;
    (&_S516)->primal_0 = _S353;
    (&_S516)->differential_0 = _S351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S517;
    (&_S517)->primal_0 = _S353;
    (&_S517)->differential_0 = _S351;
    s_bwd_prop_dot_0(&_S516, &_S517, _S515);
    float _S518 = 0.3333333432674408f * (_S352 * (*_s_dOut_2)[int(0)]);
    float3  _S519 = make_float3 (_S518, _S518, _S518);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S520;
    (&_S520)->primal_0 = _S353;
    (&_S520)->differential_0 = _S351;
    s_bwd_prop_abs_0(&_S520, _S519);
    float3  _S521 = _S517.differential_0 + _S516.differential_0 + _S520.differential_0;
    float3  _S522 = - _S521;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S476;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S412;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S477;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S452.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S486;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S493;
    float3  _S523 = _S499.differential_0 + _S498.differential_0 + _S363;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S523;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S505.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S509.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S521;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S522;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S524, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S525, DiffPair_float_0 * _S526, DiffPair_float_0 * _S527, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S528, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S529, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S530, DiffPair_float_0 * _S531, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S532, DiffPair_float_0 * _S533, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S534, bool _S535, bool _S536, bool _S537, bool _S538, FixedArray<float, 22>  * _S539)
{
    s_bwd_prop_per_pixel_losses_0(_S524, _S525, _S526, _S527, _S528, _S529, _S530, _S531, _S532, _S533, _S534, _S535, _S536, _S537, _S538, _S539);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_0, FixedArray<float, 22>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S540 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S540;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S540;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S540;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S540;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S540;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S540;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S540;
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
    float _S541 = 1.0f / ((*dpx_13).primal_0 * 52.30258560180664062f) * dOut_12;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S541;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 22>  * raw_losses_0, FixedArray<float, 10>  * weights_0, FixedArray<float, 10>  * _S542)
{
    FixedArray<float, 10>  losses_2;
    float _S543 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*weights_0)[int(0)] * (*raw_losses_0)[int(0)] / _S543;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S543)));
    bool _S544;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S544 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S544 = false;
    }
    float _S545;
    if(_S544)
    {
        _S545 = (*weights_0)[int(1)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S545 = 0.0f;
    }
    losses_2[int(2)] = _S545;
    losses_2[int(3)] = (*weights_0)[int(2)] * ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    float _S546 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(4)] = ((*weights_0)[int(3)] * (*raw_losses_0)[int(9)] + (*weights_0)[int(4)] * (*raw_losses_0)[int(10)]) / _S546;
    losses_2[int(5)] = (*weights_0)[int(5)] * (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    losses_2[int(6)] = (*weights_0)[int(6)] * (*raw_losses_0)[int(12)] / _S546;
    losses_2[int(7)] = (*weights_0)[int(7)] * (*raw_losses_0)[int(13)] / _S546;
    losses_2[int(8)] = (*weights_0)[int(8)] * (*raw_losses_0)[int(14)] / _S546;
    losses_2[int(9)] = (*weights_0)[int(9)] * (*raw_losses_0)[int(15)] / _S546;
    *_S542 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S547)
{
    return (F32_sqrt((_S547)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S548, float _S549)
{
    _d_log10_0(_S548, _S549);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_0, FixedArray<float, 10>  * weights_1, FixedArray<float, 10>  * _s_dOut_3)
{
    FixedArray<float, 22>  _S550 = dpraw_losses_0->primal_0;
    float _S551 = (*weights_1)[int(0)] * dpraw_losses_0->primal_0[int(0)];
    float _S552 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S553 = _S552 * _S552;
    float _S554 = dpraw_losses_0->primal_0[int(1)] / _S552;
    bool _S555 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S556;
    if(_S555)
    {
        _S556 = (_S550[int(3)]) != 0.0f;
    }
    else
    {
        _S556 = false;
    }
    float _S557;
    float _S558;
    float _S559;
    float _S560;
    float _S561;
    float _S562;
    float _S563;
    float _S564;
    float _S565;
    float _S566;
    float _S567;
    float _S568;
    float _S569;
    float _S570;
    float _S571;
    if(_S556)
    {
        float _S572 = _S550[int(2)] * _S550[int(3)];
        float _S573 = _S550[int(18)] * _S550[int(18)];
        float _S574 = _S550[int(6)] - _S572 / _S550[int(18)];
        float _S575 = _S550[int(2)] * _S550[int(2)];
        float _S576 = _S550[int(4)] - _S575 / _S550[int(18)];
        float _S577 = _S550[int(3)] * _S550[int(3)];
        float _S578 = _S550[int(5)] - _S577 / _S550[int(18)];
        float _S579 = _S576 * _S578 + 1.0f;
        float _S580 = s_primal_ctx_max_0(9.999999960041972e-13f, _S579);
        float _S581 = s_primal_ctx_sqrt_0(_S580);
        float _S582 = _S581 * _S581;
        float _S583 = 1.0f - _S574 / _S581;
        _S557 = (*weights_1)[int(1)];
        _S558 = _S583;
        _S559 = _S582;
        _S560 = _S574;
        _S561 = _S581;
        _S562 = _S580;
        _S563 = _S579;
        _S564 = _S576;
        _S565 = _S578;
        _S566 = _S573;
        _S567 = _S577;
        _S568 = _S550[int(3)];
        _S569 = _S575;
        _S570 = _S550[int(2)];
        _S571 = _S572;
    }
    else
    {
        _S557 = 0.0f;
        _S558 = 0.0f;
        _S559 = 0.0f;
        _S560 = 0.0f;
        _S561 = 0.0f;
        _S562 = 0.0f;
        _S563 = 0.0f;
        _S564 = 0.0f;
        _S565 = 0.0f;
        _S566 = 0.0f;
        _S567 = 0.0f;
        _S568 = 0.0f;
        _S569 = 0.0f;
        _S570 = 0.0f;
        _S571 = 0.0f;
    }
    float _S584 = s_primal_ctx_max_0(_S550[int(19)], 1.0f);
    float _S585 = _S584 * _S584;
    float _S586 = s_primal_ctx_max_0(_S550[int(20)], 1.0f);
    float _S587 = _S586 * _S586;
    float _S588 = float((I32_max((int((_S550[int(19)]) > 0.5f) + int((_S550[int(20)]) > 0.5f)), (int(1)))));
    float _S589 = (*weights_1)[int(3)] * _S550[int(9)] + (*weights_1)[int(4)] * _S550[int(10)];
    float _S590 = s_primal_ctx_max_0(_S550[int(16)], 1.0f);
    float _S591 = _S590 * _S590;
    float _S592 = (*weights_1)[int(5)] * _S550[int(11)];
    float _S593 = s_primal_ctx_max_0(_S550[int(21)], 1.0f);
    float _S594 = _S593 * _S593;
    float _S595 = (*weights_1)[int(6)] * _S550[int(12)];
    float _S596 = (*weights_1)[int(7)] * _S550[int(13)];
    float _S597 = (*weights_1)[int(8)] * _S550[int(14)];
    float _S598 = (*weights_1)[int(9)] * _S550[int(15)];
    float _S599 = (*_s_dOut_3)[int(9)] / _S591;
    float _S600 = _S598 * - _S599;
    float _S601 = (*weights_1)[int(9)] * (_S590 * _S599);
    float _S602 = (*_s_dOut_3)[int(8)] / _S591;
    float _S603 = _S597 * - _S602;
    float _S604 = (*weights_1)[int(8)] * (_S590 * _S602);
    float _S605 = (*_s_dOut_3)[int(7)] / _S591;
    float _S606 = _S596 * - _S605;
    float _S607 = (*weights_1)[int(7)] * (_S590 * _S605);
    float _S608 = (*_s_dOut_3)[int(6)] / _S591;
    float _S609 = _S595 * - _S608;
    float _S610 = (*weights_1)[int(6)] * (_S590 * _S608);
    float _S611 = (*_s_dOut_3)[int(5)] / _S594;
    float _S612 = _S592 * - _S611;
    float _S613 = _S593 * _S611;
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S550[int(21)];
    (&_S614)->differential_0 = 0.0f;
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = 1.0f;
    (&_S615)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S614, &_S615, _S612);
    float _S616 = (*weights_1)[int(5)] * _S613;
    float _S617 = (*_s_dOut_3)[int(4)] / _S591;
    float _S618 = _S590 * _S617;
    float _S619 = _S600 + _S603 + _S606 + _S609 + _S589 * - _S617;
    DiffPair_float_0 _S620;
    (&_S620)->primal_0 = _S550[int(16)];
    (&_S620)->differential_0 = 0.0f;
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = 1.0f;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S620, &_S621, _S619);
    float _S622 = (*weights_1)[int(4)] * _S618;
    float _S623 = (*weights_1)[int(3)] * _S618;
    float _S624 = (*weights_1)[int(2)] * ((*_s_dOut_3)[int(3)] / _S588);
    float _S625 = _S624 / _S587;
    float _S626 = _S550[int(8)] * - _S625;
    float _S627 = _S586 * _S625;
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = _S550[int(20)];
    (&_S628)->differential_0 = 0.0f;
    DiffPair_float_0 _S629;
    (&_S629)->primal_0 = 1.0f;
    (&_S629)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S628, &_S629, _S626);
    float _S630 = _S624 / _S585;
    float _S631 = _S550[int(7)] * - _S630;
    float _S632 = _S584 * _S630;
    DiffPair_float_0 _S633;
    (&_S633)->primal_0 = _S550[int(19)];
    (&_S633)->differential_0 = 0.0f;
    DiffPair_float_0 _S634;
    (&_S634)->primal_0 = 1.0f;
    (&_S634)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S633, &_S634, _S631);
    FixedArray<float, 22>  _S635;
    _S635[int(0)] = 0.0f;
    _S635[int(1)] = 0.0f;
    _S635[int(2)] = 0.0f;
    _S635[int(3)] = 0.0f;
    _S635[int(4)] = 0.0f;
    _S635[int(5)] = 0.0f;
    _S635[int(6)] = 0.0f;
    _S635[int(7)] = 0.0f;
    _S635[int(8)] = 0.0f;
    _S635[int(9)] = 0.0f;
    _S635[int(10)] = 0.0f;
    _S635[int(11)] = 0.0f;
    _S635[int(12)] = 0.0f;
    _S635[int(13)] = 0.0f;
    _S635[int(14)] = 0.0f;
    _S635[int(15)] = 0.0f;
    _S635[int(16)] = 0.0f;
    _S635[int(17)] = 0.0f;
    _S635[int(18)] = 0.0f;
    _S635[int(19)] = 0.0f;
    _S635[int(20)] = 0.0f;
    _S635[int(21)] = 0.0f;
    _S635[int(15)] = _S601;
    _S635[int(14)] = _S604;
    _S635[int(13)] = _S607;
    _S635[int(12)] = _S610;
    _S635[int(21)] = _S614.differential_0;
    _S635[int(11)] = _S616;
    _S635[int(16)] = _S620.differential_0;
    _S635[int(10)] = _S622;
    _S635[int(9)] = _S623;
    _S635[int(20)] = _S628.differential_0;
    _S635[int(8)] = _S627;
    _S635[int(19)] = _S633.differential_0;
    _S635[int(7)] = _S632;
    float _S636 = _S635[int(0)];
    float _S637 = _S635[int(1)];
    float _S638 = _S635[int(2)];
    float _S639 = _S635[int(3)];
    float _S640 = _S635[int(4)];
    float _S641 = _S635[int(5)];
    float _S642 = _S635[int(6)];
    float _S643 = _S635[int(7)];
    float _S644 = _S635[int(8)];
    float _S645 = _S635[int(9)];
    float _S646 = _S635[int(10)];
    float _S647 = _S635[int(11)];
    float _S648 = _S635[int(12)];
    float _S649 = _S635[int(13)];
    float _S650 = _S635[int(14)];
    float _S651 = _S635[int(15)];
    float _S652 = _S635[int(16)];
    float _S653 = _S635[int(17)];
    float _S654 = _S635[int(18)];
    float _S655 = _S635[int(19)];
    float _S656 = _S635[int(20)];
    float _S657 = _S635[int(21)];
    FixedArray<float, 22>  _S658;
    if(_S556)
    {
        float _S659 = _S557 * (*_s_dOut_3)[int(2)];
        DiffPair_float_0 _S660;
        (&_S660)->primal_0 = _S558;
        (&_S660)->differential_0 = 0.0f;
        DiffPair_float_0 _S661;
        (&_S661)->primal_0 = 0.0f;
        (&_S661)->differential_0 = 0.0f;
        DiffPair_float_0 _S662;
        (&_S662)->primal_0 = 2.0f;
        (&_S662)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S660, &_S661, &_S662, _S659);
        float _S663 = - _S660.differential_0 / _S559;
        float _S664 = _S560 * - _S663;
        float _S665 = _S561 * _S663;
        DiffPair_float_0 _S666;
        (&_S666)->primal_0 = _S562;
        (&_S666)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S666, _S664);
        DiffPair_float_0 _S667;
        (&_S667)->primal_0 = 9.999999960041972e-13f;
        (&_S667)->differential_0 = 0.0f;
        DiffPair_float_0 _S668;
        (&_S668)->primal_0 = _S563;
        (&_S668)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S667, &_S668, _S666.differential_0);
        float _S669 = _S564 * _S668.differential_0;
        float _S670 = _S565 * _S668.differential_0;
        float _S671 = - _S669 / _S566;
        float _S672 = _S568 * (_S550[int(18)] * _S671);
        float _S673 = - _S670 / _S566;
        float _S674 = _S570 * (_S550[int(18)] * _S673);
        float _S675 = - _S665 / _S566;
        float _S676 = _S550[int(18)] * _S675;
        float _S677 = _S672 + _S672 + _S570 * _S676;
        float _S678 = _S674 + _S674 + _S568 * _S676;
        float _S679 = _S567 * - _S671 + _S569 * - _S673 + _S571 * - _S675;
        FixedArray<float, 22>  _S680;
        _S680[int(0)] = 0.0f;
        _S680[int(1)] = 0.0f;
        _S680[int(2)] = 0.0f;
        _S680[int(3)] = 0.0f;
        _S680[int(4)] = 0.0f;
        _S680[int(5)] = 0.0f;
        _S680[int(6)] = 0.0f;
        _S680[int(7)] = 0.0f;
        _S680[int(8)] = 0.0f;
        _S680[int(9)] = 0.0f;
        _S680[int(10)] = 0.0f;
        _S680[int(11)] = 0.0f;
        _S680[int(12)] = 0.0f;
        _S680[int(13)] = 0.0f;
        _S680[int(14)] = 0.0f;
        _S680[int(15)] = 0.0f;
        _S680[int(16)] = 0.0f;
        _S680[int(17)] = 0.0f;
        _S680[int(18)] = 0.0f;
        _S680[int(19)] = 0.0f;
        _S680[int(20)] = 0.0f;
        _S680[int(21)] = 0.0f;
        _S680[int(5)] = _S669;
        _S680[int(4)] = _S670;
        _S680[int(3)] = _S677;
        _S680[int(2)] = _S678;
        _S680[int(6)] = _S665;
        float _S681 = _S637 + _S680[int(1)];
        float _S682 = _S638 + _S680[int(2)];
        float _S683 = _S639 + _S680[int(3)];
        float _S684 = _S640 + _S680[int(4)];
        float _S685 = _S641 + _S680[int(5)];
        float _S686 = _S642 + _S680[int(6)];
        float _S687 = _S643 + _S680[int(7)];
        float _S688 = _S644 + _S680[int(8)];
        float _S689 = _S645 + _S680[int(9)];
        float _S690 = _S646 + _S680[int(10)];
        float _S691 = _S647 + _S680[int(11)];
        float _S692 = _S648 + _S680[int(12)];
        float _S693 = _S649 + _S680[int(13)];
        float _S694 = _S650 + _S680[int(14)];
        float _S695 = _S651 + _S680[int(15)];
        float _S696 = _S652 + _S680[int(16)];
        float _S697 = _S653 + _S680[int(17)];
        float _S698 = _S654 + _S680[int(18)];
        float _S699 = _S655 + _S680[int(19)];
        float _S700 = _S656 + _S680[int(20)];
        float _S701 = _S657 + _S680[int(21)];
        _S658[int(0)] = _S636 + _S680[int(0)];
        _S658[int(1)] = _S681;
        _S658[int(2)] = _S682;
        _S658[int(3)] = _S683;
        _S658[int(4)] = _S684;
        _S658[int(5)] = _S685;
        _S658[int(6)] = _S686;
        _S658[int(7)] = _S687;
        _S658[int(8)] = _S688;
        _S658[int(9)] = _S689;
        _S658[int(10)] = _S690;
        _S658[int(11)] = _S691;
        _S658[int(12)] = _S692;
        _S658[int(13)] = _S693;
        _S658[int(14)] = _S694;
        _S658[int(15)] = _S695;
        _S658[int(16)] = _S696;
        _S658[int(17)] = _S697;
        _S658[int(18)] = _S698;
        _S658[int(19)] = _S699;
        _S658[int(20)] = _S700;
        _S658[int(21)] = _S701;
        _S557 = _S679;
    }
    else
    {
        _S658[int(0)] = _S636;
        _S658[int(1)] = _S637;
        _S658[int(2)] = _S638;
        _S658[int(3)] = _S639;
        _S658[int(4)] = _S640;
        _S658[int(5)] = _S641;
        _S658[int(6)] = _S642;
        _S658[int(7)] = _S643;
        _S658[int(8)] = _S644;
        _S658[int(9)] = _S645;
        _S658[int(10)] = _S646;
        _S658[int(11)] = _S647;
        _S658[int(12)] = _S648;
        _S658[int(13)] = _S649;
        _S658[int(14)] = _S650;
        _S658[int(15)] = _S651;
        _S658[int(16)] = _S652;
        _S658[int(17)] = _S653;
        _S658[int(18)] = _S654;
        _S658[int(19)] = _S655;
        _S658[int(20)] = _S656;
        _S658[int(21)] = _S657;
        _S557 = 0.0f;
    }
    if(_S555)
    {
        FixedArray<float, 22>  _S702;
        _S702[int(0)] = 0.0f;
        _S702[int(1)] = 0.0f;
        _S702[int(2)] = 0.0f;
        _S702[int(3)] = 0.0f;
        _S702[int(4)] = 0.0f;
        _S702[int(5)] = 0.0f;
        _S702[int(6)] = 0.0f;
        _S702[int(7)] = 0.0f;
        _S702[int(8)] = 0.0f;
        _S702[int(9)] = 0.0f;
        _S702[int(10)] = 0.0f;
        _S702[int(11)] = 0.0f;
        _S702[int(12)] = 0.0f;
        _S702[int(13)] = 0.0f;
        _S702[int(14)] = 0.0f;
        _S702[int(15)] = 0.0f;
        _S702[int(16)] = 0.0f;
        _S702[int(17)] = 0.0f;
        _S702[int(18)] = 0.0f;
        _S702[int(19)] = 0.0f;
        _S702[int(20)] = 0.0f;
        _S702[int(21)] = 0.0f;
        _S702[int(3)] = 0.0f;
        float _S703 = _S658[int(1)] + _S702[int(1)];
        float _S704 = _S658[int(2)] + _S702[int(2)];
        float _S705 = _S658[int(3)] + _S702[int(3)];
        float _S706 = _S658[int(4)] + _S702[int(4)];
        float _S707 = _S658[int(5)] + _S702[int(5)];
        float _S708 = _S658[int(6)] + _S702[int(6)];
        float _S709 = _S658[int(7)] + _S702[int(7)];
        float _S710 = _S658[int(8)] + _S702[int(8)];
        float _S711 = _S658[int(9)] + _S702[int(9)];
        float _S712 = _S658[int(10)] + _S702[int(10)];
        float _S713 = _S658[int(11)] + _S702[int(11)];
        float _S714 = _S658[int(12)] + _S702[int(12)];
        float _S715 = _S658[int(13)] + _S702[int(13)];
        float _S716 = _S658[int(14)] + _S702[int(14)];
        float _S717 = _S658[int(15)] + _S702[int(15)];
        float _S718 = _S658[int(16)] + _S702[int(16)];
        float _S719 = _S658[int(17)] + _S702[int(17)];
        float _S720 = _S658[int(18)] + _S702[int(18)];
        float _S721 = _S658[int(19)] + _S702[int(19)];
        float _S722 = _S658[int(20)] + _S702[int(20)];
        float _S723 = _S658[int(21)] + _S702[int(21)];
        _S658[int(0)] = _S658[int(0)] + _S702[int(0)];
        _S658[int(1)] = _S703;
        _S658[int(2)] = _S704;
        _S658[int(3)] = _S705;
        _S658[int(4)] = _S706;
        _S658[int(5)] = _S707;
        _S658[int(6)] = _S708;
        _S658[int(7)] = _S709;
        _S658[int(8)] = _S710;
        _S658[int(9)] = _S711;
        _S658[int(10)] = _S712;
        _S658[int(11)] = _S713;
        _S658[int(12)] = _S714;
        _S658[int(13)] = _S715;
        _S658[int(14)] = _S716;
        _S658[int(15)] = _S717;
        _S658[int(16)] = _S718;
        _S658[int(17)] = _S719;
        _S658[int(18)] = _S720;
        _S658[int(19)] = _S721;
        _S658[int(20)] = _S722;
        _S658[int(21)] = _S723;
    }
    float _S724 = -10.0f * (*_s_dOut_3)[int(1)];
    DiffPair_float_0 _S725;
    (&_S725)->primal_0 = _S554;
    (&_S725)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S725, _S724);
    float _S726 = _S725.differential_0 / _S553;
    float _S727 = _S552 * _S726;
    float _S728 = (*_s_dOut_3)[int(0)] / _S553;
    float _S729 = _S552 * _S728;
    float _S730 = _S550[int(1)] * - _S726 + _S551 * - _S728;
    DiffPair_float_0 _S731;
    (&_S731)->primal_0 = _S550[int(17)];
    (&_S731)->differential_0 = 0.0f;
    DiffPair_float_0 _S732;
    (&_S732)->primal_0 = 1.0f;
    (&_S732)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S731, &_S732, _S730);
    float _S733 = (*weights_1)[int(0)] * _S729;
    FixedArray<float, 22>  _S734;
    _S734[int(0)] = 0.0f;
    _S734[int(1)] = 0.0f;
    _S734[int(2)] = 0.0f;
    _S734[int(3)] = 0.0f;
    _S734[int(4)] = 0.0f;
    _S734[int(5)] = 0.0f;
    _S734[int(6)] = 0.0f;
    _S734[int(7)] = 0.0f;
    _S734[int(8)] = 0.0f;
    _S734[int(9)] = 0.0f;
    _S734[int(10)] = 0.0f;
    _S734[int(11)] = 0.0f;
    _S734[int(12)] = 0.0f;
    _S734[int(13)] = 0.0f;
    _S734[int(14)] = 0.0f;
    _S734[int(15)] = 0.0f;
    _S734[int(16)] = 0.0f;
    _S734[int(17)] = 0.0f;
    _S734[int(18)] = 0.0f;
    _S734[int(19)] = 0.0f;
    _S734[int(20)] = 0.0f;
    _S734[int(21)] = 0.0f;
    _S734[int(18)] = _S557;
    _S734[int(1)] = _S727;
    _S734[int(17)] = _S731.differential_0;
    _S734[int(0)] = _S733;
    FixedArray<float, 22>  _S735 = {
        _S658[int(0)] + _S734[int(0)], _S658[int(1)] + _S734[int(1)], _S658[int(2)] + _S734[int(2)], _S658[int(3)] + _S734[int(3)], _S658[int(4)] + _S734[int(4)], _S658[int(5)] + _S734[int(5)], _S658[int(6)] + _S734[int(6)], _S658[int(7)] + _S734[int(7)], _S658[int(8)] + _S734[int(8)], _S658[int(9)] + _S734[int(9)], _S658[int(10)] + _S734[int(10)], _S658[int(11)] + _S734[int(11)], _S658[int(12)] + _S734[int(12)], _S658[int(13)] + _S734[int(13)], _S658[int(14)] + _S734[int(14)], _S658[int(15)] + _S734[int(15)], _S658[int(16)] + _S734[int(16)], _S658[int(17)] + _S734[int(17)], _S658[int(18)] + _S734[int(18)], _S658[int(19)] + _S734[int(19)], _S658[int(20)] + _S734[int(20)], _S658[int(21)] + _S734[int(21)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S735;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S736, FixedArray<float, 10>  * _S737, FixedArray<float, 10>  * _S738)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S736, _S737, _S738);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 22>  * raw_losses_1, FixedArray<float, 10>  * weights_2, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 22>  * _S739)
{
    FixedArray<float, 22>  _S740 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S740;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_2, v_losses_1);
    *_S739 = (&dp_raw_losses_0)->differential_0;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S741, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S742, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S743, float3  _S744)
{
    _d_clamp_vector_0(_S741, _S742, _S743, _S744);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_4)
{
    float _S745 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S746 = make_float3 (_S745);
    float3  _S747 = make_float3 (0.0f);
    float3  _S748 = make_float3 (1.0f);
    float3  _S749 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S750;
    (&_S750)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S745) * (*dpbackground_0).primal_0;
    (&_S750)->differential_0 = _S749;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S751;
    (&_S751)->primal_0 = _S747;
    (&_S751)->differential_0 = _S749;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S752;
    (&_S752)->primal_0 = _S748;
    (&_S752)->differential_0 = _S749;
    s_bwd_prop_clamp_1(&_S750, &_S751, &_S752, _s_dOut_4);
    float3  _S753 = _S746 * _S750.differential_0;
    float3  _S754 = (*dpbackground_0).primal_0 * _S750.differential_0;
    float _S755 = - (_S754.x + _S754.y + _S754.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S753;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S755;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S750.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S756, DiffPair_float_0 * _S757, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S758, float3  _S759)
{
    s_bwd_prop_blend_background_0(_S756, _S757, _S758, _S759);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S760 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S760;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S760;
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

inline __device__ Matrix<float, 2, 2>  camera_distortion_jac_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float _S761 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S762 = (*dist_coeffs_0)[int(1)] + r2_0 * _S761;
    float _S763 = (*dist_coeffs_0)[int(0)] + r2_0 * _S762;
    float2  _S764 = make_float2 (1.0f + r2_0 * _S763);
    float _S765 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S766 = _S765 * u_0;
    float _S767 = 2.0f * u_0;
    float _S768 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S769 = _S768 * u_0;
    float _S770 = 2.0f * v_0;
    float2  _S771 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S772 = uv_0 * _S771;
    float _S773 = (*dist_coeffs_0)[int(4)] * _S771.y;
    float _S774 = (*dist_coeffs_0)[int(5)] * _S771.x;
    float _S775 = _S772.x + _S772.y;
    float _S776 = r2_0 * _S775;
    float _S777 = r2_0 * _S776;
    float _S778 = (*dist_coeffs_0)[int(7)] * _S771.y + _S773 + (*dist_coeffs_0)[int(6)] * _S771.x + _S774 + _S763 * _S775 + _S762 * _S776 + _S761 * _S777 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S777);
    float _S779 = v_0 * _S778;
    float _S780 = u_0 * _S778;
    float2  _S781 = make_float2 (0.0f, 1.0f) + make_float2 ((*dist_coeffs_0)[int(8)] * 0.0f, (*dist_coeffs_0)[int(9)] * 0.0f);
    float2  _S782 = uv_0 * _S781;
    float _S783 = (*dist_coeffs_0)[int(4)] * _S781.y;
    float _S784 = (*dist_coeffs_0)[int(5)] * _S781.x;
    float _S785 = _S782.x + _S782.y;
    float _S786 = r2_0 * _S785;
    float _S787 = r2_0 * _S786;
    float _S788 = (*dist_coeffs_0)[int(7)] * _S781.y + _S783 + (*dist_coeffs_0)[int(6)] * _S781.x + _S784 + _S763 * _S785 + _S762 * _S786 + _S761 * _S787 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S787);
    float _S789 = v_0 * _S788;
    float _S790 = u_0 * _S788;
    return makeMatrix<float, 2, 2> (_S764 * _S771 + make_float2 (_S768 * (v_0 * _S771.y) + _S767 * _S774 + 2.0f * (u_0 * _S774) + _S765 * (v_0 * _S771.x) + _S780 + _S780, _S770 * _S773 + 2.0f * (v_0 * _S773) + _S769 * _S771.y + _S766 * _S771.x + _S779 + _S779), _S764 * _S781 + make_float2 (_S768 * (v_0 * _S781.y) + _S767 * _S784 + 2.0f * (u_0 * _S784) + _S765 * (v_0 * _S781.x) + _S790 + _S790, _S770 * _S783 + 2.0f * (v_0 * _S783) + _S769 * _S781.y + _S766 * _S781.x + _S789 + _S789));
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  * dist_coeffs_1)
{
    Matrix<float, 2, 2>  _S791 = camera_distortion_jac_0(uv_1, dist_coeffs_1);
    return (F32_min((determinant_0(_S791)), ((F32_min((_S791.rows[int(0)].x), (_S791.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_2)
{
    float2  _S792;
    if(is_fisheye_0)
    {
        float r_2 = length_1(uv_2);
        float theta_0 = (F32_atan((r_2)));
        float _S793;
        if(r_2 < 0.00100000004749745f)
        {
            _S793 = 1.0f - r_2 * r_2 / 3.0f;
        }
        else
        {
            _S793 = theta_0 / r_2;
        }
        _S792 = uv_2 * make_float2 (_S793);
    }
    else
    {
        _S792 = uv_2;
    }
    float u_1 = _S792.x;
    float v_1 = _S792.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S794 = _S792 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_2)[int(0)] + r2_1 * ((*dist_coeffs_2)[int(1)] + r2_1 * ((*dist_coeffs_2)[int(2)] + r2_1 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_1 * v_1 + (*dist_coeffs_2)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_2)[int(6)] * r2_1, 2.0f * (*dist_coeffs_2)[int(5)] * u_1 * v_1 + (*dist_coeffs_2)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_2)[int(7)] * r2_1);
    return _S794 + make_float2 ((*dist_coeffs_2)[int(8)] * _S794.x + (*dist_coeffs_2)[int(9)] * _S794.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_3, FixedArray<float, 10>  * dist_coeffs_3, int maxiter_0, float2  * uv_undist_0)
{
    int i_10 = int(0);
    float2  q_0 = uv_3;
    for(;;)
    {
        if(i_10 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S795 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_3 = _S795 + make_float2 ((*dist_coeffs_3)[int(8)] * _S795.x + (*dist_coeffs_3)[int(9)] * _S795.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S796 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S796.rows[int(0)].x * _S796.rows[int(1)].y - _S796.rows[int(0)].y * _S796.rows[int(1)].x);
        float _S797 = r_3.x;
        float _S798 = r_3.y;
        float2  q_1 = q_0 - make_float2 ((_S797 * _S796.rows[int(1)].y - _S798 * _S796.rows[int(0)].y) * inv_det_0, (- _S797 * _S796.rows[int(1)].x + _S798 * _S796.rows[int(0)].x) * inv_det_0);
        i_10 = i_10 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S799 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S800;
    if((F32_min((determinant_0(_S799)), ((F32_min((_S799.rows[int(0)].x), (_S799.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S801 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S800 = (length_1(_S801 + make_float2 ((*dist_coeffs_3)[int(8)] * _S801.x + (*dist_coeffs_3)[int(9)] * _S801.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S800 = false;
    }
    return _S800;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S802 = uv_4;
    bool _S803 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S802);
    if(!_S803)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S804 = _S802;
        float theta_1 = length_1(_S802);
        float _S805;
        if(theta_1 < 0.00100000004749745f)
        {
            _S805 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S805 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S806 = make_float3 ((_S804 * make_float2 (_S805)).x, (_S804 * make_float2 (_S805)).y, (F32_cos((theta_1))));
        raydir_0 = _S806;
    }
    else
    {
        raydir_0 = make_float3 (_S802.x, _S802.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_5, float3  * raydir_1)
{
    float2  _S807 = uv_5;
    int3  _S808 = make_int3 (int(0));
    float3  _S809 = make_float3 ((float)_S808.x, (float)_S808.y, (float)_S808.z);
    *raydir_1 = _S809;
    bool _S810 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S807);
    if(!_S810)
    {
        return false;
    }
    float3  _S811;
    if(is_fisheye_2)
    {
        float2  _S812 = _S807;
        float theta_2 = length_1(_S807);
        float _S813;
        if(theta_2 < 0.00100000004749745f)
        {
            _S813 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S813 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S814 = make_float3 ((_S812 * make_float2 (_S813)).x, (_S812 * make_float2 (_S813)).y, (F32_cos((theta_2))));
        _S811 = _S814;
    }
    else
    {
        _S811 = make_float3 (_S807.x, _S807.y, 1.0f);
    }
    *raydir_1 = _S811;
    return true;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_14)
{
    float _S815 = (*right_2).primal_0.rows[int(0)].x * dOut_14.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_14.x;
    float sum_4 = _S815 + (*right_2).primal_0.rows[int(0)].y * dOut_14.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_14.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(0)].z * dOut_14.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_14.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_5;
    float _S816 = (*right_2).primal_0.rows[int(1)].x * dOut_14.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_14.x;
    float sum_6 = _S816 + (*right_2).primal_0.rows[int(1)].y * dOut_14.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_14.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(1)].z * dOut_14.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_14.z;
    *&((&left_d_result_1)->y) = sum_7;
    float _S817 = (*right_2).primal_0.rows[int(2)].x * dOut_14.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_14.x;
    float sum_8 = _S817 + (*right_2).primal_0.rows[int(2)].y * dOut_14.y;
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

inline __device__ bool generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_6, float3  * ray_o_0, float3  * ray_d_0)
{
    float2  _S818 = uv_6;
    *ray_o_0 = - mul_2(t_1, R_2);
    bool _S819 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S818);
    if(!_S819)
    {
        return false;
    }
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float2  _S820 = _S818;
        float theta_3 = length_1(_S818);
        float _S821;
        if(theta_3 < 0.00100000004749745f)
        {
            _S821 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S821 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S822 = make_float3 ((_S820 * make_float2 (_S821)).x, (_S820 * make_float2 (_S821)).y, (F32_cos((theta_3))));
        raydir_2 = _S822;
    }
    else
    {
        raydir_2 = make_float3 (_S818.x, _S818.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_2(raydir_2, R_2));
    return true;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S823;
    bool _S824;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S825, Matrix<float, 3, 3>  _S826)
{
    return mul_2(_S825, _S826);
}

inline __device__ float s_primal_ctx_sin_0(float _S827)
{
    return (F32_sin((_S827)));
}

inline __device__ float s_primal_ctx_cos_0(float _S828)
{
    return (F32_cos((_S828)));
}

inline __device__ bool s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_7, bool is_fisheye_4, FixedArray<float, 10>  * dist_coeffs_7, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S823 = make_float2 (0.0f);
    _s_diff_ctx_0->_S824 = false;
    float3  _S829 = make_float3 (0.0f);
    float3  _S830 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    float2  _S831 = uv_7;
    bool _S832 = undistort_point_0(uv_7, dist_coeffs_7, int(8), &_S831);
    _s_diff_ctx_0->_S823 = _S831;
    _s_diff_ctx_0->_S824 = _S832;
    float2  _S833 = _S831;
    float3  raydir_3;
    bool _S834;
    if(!!_S832)
    {
        if(is_fisheye_4)
        {
            float _S835 = length_1(_S833);
            float _S836;
            if(_S835 < 0.00100000004749745f)
            {
                _S836 = 1.0f - _S835 * _S835 / 6.0f;
            }
            else
            {
                _S836 = s_primal_ctx_sin_0(_S835) / _S835;
            }
            float3  _S837 = make_float3 ((_S833 * make_float2 (_S836)).x, (_S833 * make_float2 (_S836)).y, s_primal_ctx_cos_0(_S835));
            raydir_3 = _S837;
        }
        else
        {
            raydir_3 = make_float3 (_S833.x, _S833.y, 1.0f);
        }
        float3  _S838 = normalize_0(s_primal_ctx_mul_0(raydir_3, dpR_0));
        _S834 = true;
        raydir_3 = _S838;
    }
    else
    {
        _S834 = false;
        raydir_3 = _S829;
    }
    *dpray_o_0 = _S830;
    *dpray_d_0 = raydir_3;
    return _S834;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float _s_dOut_5)
{
    float _S839 = (*dpx_15).primal_0.x;
    float _S840 = (*dpx_15).primal_0.y;
    float _S841 = (*dpx_15).primal_0.z;
    DiffPair_float_0 _S842;
    (&_S842)->primal_0 = _S839 * _S839 + _S840 * _S840 + _S841 * _S841;
    (&_S842)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S842, _s_dOut_5);
    float _S843 = (*dpx_15).primal_0.z * _S842.differential_0;
    float _S844 = _S843 + _S843;
    float _S845 = (*dpx_15).primal_0.y * _S842.differential_0;
    float _S846 = _S845 + _S845;
    float _S847 = (*dpx_15).primal_0.x * _S842.differential_0;
    float _S848 = _S847 + _S847;
    float3  _S849 = make_float3 (0.0f);
    *&((&_S849)->z) = _S844;
    *&((&_S849)->y) = _S846;
    *&((&_S849)->x) = _S848;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S849;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S850, float _S851)
{
    s_bwd_prop_length_impl_1(_S850, _S851);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  _s_dOut_6)
{
    float _S852 = length_2((*dpx_16).primal_0);
    float3  _S853 = (*dpx_16).primal_0 * _s_dOut_6;
    float3  _S854 = make_float3 (1.0f / _S852) * _s_dOut_6;
    float _S855 = - ((_S853.x + _S853.y + _S853.z) / (_S852 * _S852));
    float3  _S856 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S857;
    (&_S857)->primal_0 = (*dpx_16).primal_0;
    (&_S857)->differential_0 = _S856;
    s_bwd_length_impl_1(&_S857, _S855);
    float3  _S858 = _S854 + _S857.differential_0;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S858;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S859, float3  _S860)
{
    s_bwd_prop_normalize_impl_0(_S859, _S860);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S861, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S862, float3  _S863)
{
    _d_mul_0(_S861, _S862, _S863);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_8, bool is_fisheye_5, FixedArray<float, 10>  * dist_coeffs_8, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S864 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S865 = *dpt_1;
    float3  _S866 = make_float3 (0.0f);
    bool _S867 = !!_s_diff_ctx_1->_S824;
    float3  raydir_4;
    float3  _S868;
    if(_S867)
    {
        if(is_fisheye_5)
        {
            float _S869 = length_1(_s_diff_ctx_1->_S823);
            float _S870;
            if(_S869 < 0.00100000004749745f)
            {
                _S870 = 1.0f - _S869 * _S869 / 6.0f;
            }
            else
            {
                _S870 = s_primal_ctx_sin_0(_S869) / _S869;
            }
            float3  _S871 = make_float3 ((_s_diff_ctx_1->_S823 * make_float2 (_S870)).x, (_s_diff_ctx_1->_S823 * make_float2 (_S870)).y, s_primal_ctx_cos_0(_S869));
            raydir_4 = _S871;
        }
        else
        {
            raydir_4 = make_float3 (_s_diff_ctx_1->_S823.x, _s_diff_ctx_1->_S823.y, 1.0f);
        }
        float3  _S872 = raydir_4;
        raydir_4 = s_primal_ctx_mul_0(raydir_4, _S864.primal_0);
        _S868 = _S872;
    }
    else
    {
        raydir_4 = _S866;
        _S868 = _S866;
    }
    Matrix<float, 3, 3>  _S873 = makeMatrix<float, 3, 3> (0.0f);
    Matrix<float, 3, 3>  _S874;
    if(_S867)
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S875;
        (&_S875)->primal_0 = raydir_4;
        (&_S875)->differential_0 = _S866;
        s_bwd_normalize_impl_0(&_S875, dpray_d_1);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S876;
        (&_S876)->primal_0 = _S868;
        (&_S876)->differential_0 = _S866;
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S877;
        (&_S877)->primal_0 = _S864.primal_0;
        (&_S877)->differential_0 = _S873;
        s_bwd_prop_mul_0(&_S876, &_S877, _S875.differential_0);
        _S874 = _S877.differential_0;
    }
    else
    {
        _S874 = _S873;
    }
    float3  _S878 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S879;
    (&_S879)->primal_0 = _S865.primal_0;
    (&_S879)->differential_0 = _S866;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S880;
    (&_S880)->primal_0 = _S864.primal_0;
    (&_S880)->differential_0 = _S873;
    s_bwd_prop_mul_0(&_S879, &_S880, _S878);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S879.differential_0;
    Matrix<float, 3, 3>  _S881 = _S880.differential_0 + _S874;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S881;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S882, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S883, float2  _S884, bool _S885, FixedArray<float, 10>  * _S886, float3  _S887, float3  _S888)
{
    float3  _S889;
    float3  _S890;
    s_bwd_prop_generate_ray_Intermediates_0 _S891;
    bool _S892 = s_primal_ctx_generate_ray_0((*_S882).primal_0, (*_S883).primal_0, _S884, _S885, _S886, &_S889, &_S890, &_S891);
    s_bwd_prop_generate_ray_Intermediates_0 _S893 = _S891;
    s_bwd_prop_generate_ray_0(_S882, _S883, _S884, _S885, _S886, _S887, _S888, &_S893);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_9, bool is_fisheye_6, FixedArray<float, 10>  * dist_coeffs_9, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S894 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S894;
    float3  _S895 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S895;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_9, is_fisheye_6, dist_coeffs_9, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S896 = dOut_15.y;
    float _S897 = dOut_15.z;
    float _S898 = dOut_15.x;
    float _S899 = (*a_0).primal_0.z * _S896 + - (*a_0).primal_0.y * _S897;
    float _S900 = - (*a_0).primal_0.z * _S898 + (*a_0).primal_0.x * _S897;
    float _S901 = (*a_0).primal_0.y * _S898 + - (*a_0).primal_0.x * _S896;
    float3  _S902 = make_float3 (- (*b_0).primal_0.z * _S896 + (*b_0).primal_0.y * _S897, (*b_0).primal_0.z * _S898 + - (*b_0).primal_0.x * _S897, - (*b_0).primal_0.y * _S898 + (*b_0).primal_0.x * _S896);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S902;
    float3  _S903 = make_float3 (_S899, _S900, _S901);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S903;
    return;
}

inline __device__ float3  cross_0(float3  left_4, float3  right_4)
{
    float _S904 = left_4.y;
    float _S905 = right_4.z;
    float _S906 = left_4.z;
    float _S907 = right_4.y;
    float _S908 = right_4.x;
    float _S909 = left_4.x;
    return make_float3 (_S904 * _S905 - _S906 * _S907, _S906 * _S908 - _S909 * _S905, _S909 * _S907 - _S904 * _S908);
}

inline __device__ void depth_to_normal(uint width_0, uint height_0, float2  pix_center_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_0, float4  depths_0, float3  * normal_0)
{
    FixedArray<float3 , 4>  points_0;
    float2  _S910 = float2 {intrins_0.z, intrins_0.w};
    float2  _S911 = float2 {intrins_0.x, intrins_0.y};
    float2  _S912 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S910) / _S911;
    float2  uv_10 = _S912;
    bool _S913 = undistort_point_0(_S912, dist_coeffs_10, int(12), &uv_10);
    if(!_S913)
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
    float2  _S914 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S910) / _S911;
    float2  uv_11 = _S914;
    bool _S915 = undistort_point_0(_S914, dist_coeffs_10, int(12), &uv_11);
    if(!_S915)
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
    float2  _S916 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S910) / _S911;
    float2  uv_12 = _S916;
    bool _S917 = undistort_point_0(_S916, dist_coeffs_10, int(12), &uv_12);
    if(!_S917)
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
    float2  _S918 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S910) / _S911;
    float2  uv_13 = _S918;
    bool _S919 = undistort_point_0(_S918, dist_coeffs_10, int(12), &uv_13);
    if(!_S919)
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
    float3  _S920 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    *normal_0 = _S920;
    if((length_2(_S920)) != 0.0f)
    {
        *normal_0 = *normal_0 / make_float3 (length_2(*normal_0));
    }
    return;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S921;
    bool _S922;
    float2  _S923;
    bool _S924;
    float2  _S925;
    bool _S926;
    float2  _S927;
    bool _S928;
};

inline __device__ float3  s_primal_ctx_cross_0(float3  _S929, float3  _S930)
{
    return cross_0(_S929, _S930);
}

inline __device__ void s_primal_ctx_depth_to_normal_0(uint width_1, uint height_1, float2  pix_center_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_1, float4  dpdepths_0, float3  * dpnormal_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S931 = make_float2 (0.0f);
    _s_diff_ctx_2->_S921 = _S931;
    _s_diff_ctx_2->_S922 = false;
    _s_diff_ctx_2->_S923 = _S931;
    _s_diff_ctx_2->_S924 = false;
    _s_diff_ctx_2->_S925 = _S931;
    _s_diff_ctx_2->_S926 = false;
    _s_diff_ctx_2->_S927 = _S931;
    _s_diff_ctx_2->_S928 = false;
    _s_diff_ctx_2->_S923 = _S931;
    _s_diff_ctx_2->_S924 = false;
    _s_diff_ctx_2->_S925 = _S931;
    _s_diff_ctx_2->_S926 = false;
    _s_diff_ctx_2->_S927 = _S931;
    _s_diff_ctx_2->_S928 = false;
    float3  _S932 = make_float3 (0.0f);
    float2  _S933 = float2 {intrins_1.z, intrins_1.w};
    float2  _S934 = float2 {intrins_1.x, intrins_1.y};
    float2  _S935 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S933) / _S934;
    float2  _S936 = _S935;
    bool _S937 = undistort_point_0(_S935, dist_coeffs_11, int(12), &_S936);
    _s_diff_ctx_2->_S921 = _S936;
    _s_diff_ctx_2->_S922 = _S937;
    float2  uv_14 = _S936;
    bool _S938 = !_S937;
    float3  _S939;
    if(_S938)
    {
        _S939 = make_float3 (0.0f);
    }
    else
    {
        _S939 = _S932;
    }
    bool _S940 = !_S938;
    int _S941;
    FixedArray<float3 , 4>  points_1;
    if(_S940)
    {
        float3  raydir_14;
        if(is_fisheye_8)
        {
            float _S942 = length_1(uv_14);
            float3  raydir_15 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S942, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S942))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S942, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S942))).y, s_primal_ctx_cos_0(_S942));
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
        float3  _S943 = make_float3 (dpdepths_0.x) * raydir_14;
        float2  _S944 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S933) / _S934;
        float2  _S945 = _S944;
        bool _S946 = undistort_point_0(_S944, dist_coeffs_11, int(12), &_S945);
        _s_diff_ctx_2->_S923 = _S945;
        _s_diff_ctx_2->_S924 = _S946;
        float2  uv_15 = _S945;
        bool _S947 = !_S946;
        if(_S947)
        {
            _S939 = make_float3 (0.0f);
        }
        bool _S948 = !_S947;
        if(_S948)
        {
            if(is_fisheye_8)
            {
                float _S949 = length_1(uv_15);
                float3  raydir_17 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S949, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S949))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S949, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S949))).y, s_primal_ctx_cos_0(_S949));
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
            float3  _S950 = make_float3 (dpdepths_0.y) * raydir_14;
            _S941 = int(2);
            points_1[int(0)] = _S943;
            points_1[int(1)] = _S950;
            points_1[int(2)] = _S932;
            points_1[int(3)] = _S932;
        }
        else
        {
            _S941 = int(0);
            points_1[int(0)] = _S943;
            points_1[int(1)] = _S932;
            points_1[int(2)] = _S932;
            points_1[int(3)] = _S932;
        }
        bool _runFlag_0;
        if(_S941 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S940;
            _S941 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S951 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S933) / _S934;
            float2  _S952 = _S951;
            bool _S953 = undistort_point_0(_S951, dist_coeffs_11, int(12), &_S952);
            _s_diff_ctx_2->_S925 = _S952;
            _s_diff_ctx_2->_S926 = _S953;
            float2  uv_16 = _S952;
            if(!_S953)
            {
                float3  _S954 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S941 = int(0);
                _S939 = _S954;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_8)
                {
                    float _S955 = length_1(uv_16);
                    float3  raydir_19 = make_float3 ((uv_16 / make_float2 (s_primal_ctx_max_0(_S955, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S955))).x, (uv_16 / make_float2 (s_primal_ctx_max_0(_S955, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S955))).y, s_primal_ctx_cos_0(_S955));
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
                float2  _S956 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S933) / _S934;
                float2  _S957 = _S956;
                bool _S958 = undistort_point_0(_S956, dist_coeffs_11, int(12), &_S957);
                _s_diff_ctx_2->_S927 = _S957;
                _s_diff_ctx_2->_S928 = _S958;
                float2  uv_17 = _S957;
                bool _S959 = !_S958;
                if(_S959)
                {
                    _S939 = make_float3 (0.0f);
                }
                bool _S960 = !_S959;
                int _S961;
                if(_S960)
                {
                    if(is_fisheye_8)
                    {
                        float _S962 = length_1(uv_17);
                        float3  raydir_21 = make_float3 ((uv_17 / make_float2 (s_primal_ctx_max_0(_S962, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S962))).x, (uv_17 / make_float2 (s_primal_ctx_max_0(_S962, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S962))).y, s_primal_ctx_cos_0(_S962));
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
                    _S961 = int(2);
                }
                else
                {
                    _S961 = int(0);
                }
                if(_S961 != int(2))
                {
                    _runFlag_0 = false;
                    _S941 = _S961;
                }
                if(_runFlag_0)
                {
                    _S941 = int(1);
                }
            }
        }
    }
    else
    {
        _S941 = int(0);
        points_1[int(0)] = _S932;
        points_1[int(1)] = _S932;
        points_1[int(2)] = _S932;
        points_1[int(3)] = _S932;
    }
    if(!(_S941 != int(1)))
    {
        float3  _S963 = s_primal_ctx_cross_0(points_1[int(1)] - points_1[int(0)], - (points_1[int(3)] - points_1[int(2)]));
        float _S964 = length_2(_S963);
        if(_S964 != 0.0f)
        {
            _S939 = _S963 / make_float3 (_S964);
        }
        else
        {
            _S939 = _S963;
        }
    }
    *dpnormal_0 = _S939;
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S965, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S966, float3  _S967)
{
    _d_cross_0(_S965, _S966, _S967);
    return;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(uint width_2, uint height_2, float2  pix_center_2, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9, bool is_ray_depth_2, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  dpnormal_1, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S968 = *dpdepths_1;
    float3  _S969 = make_float3 (0.0f);
    bool _S970 = !!_s_diff_ctx_3->_S922;
    float3  raydir_23;
    float3  raydir_24;
    float3  raydir_25;
    float3  raydir_26;
    int _S971;
    FixedArray<float3 , 4>  points_2;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S972;
    if(_S970)
    {
        if(is_fisheye_9)
        {
            float _S973 = length_1(_s_diff_ctx_3->_S921);
            float3  raydir_27 = make_float3 ((_s_diff_ctx_3->_S921 / make_float2 (s_primal_ctx_max_0(_S973, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S973))).x, (_s_diff_ctx_3->_S921 / make_float2 (s_primal_ctx_max_0(_S973, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S973))).y, s_primal_ctx_cos_0(_S973));
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
            float3  raydir_28 = make_float3 (_s_diff_ctx_3->_S921.x, _s_diff_ctx_3->_S921.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_23 = normalize_0(raydir_28);
            }
            else
            {
                raydir_23 = raydir_28;
            }
        }
        float3  _S974 = make_float3 (_S968.primal_0.x) * raydir_23;
        bool _S975 = !!_s_diff_ctx_3->_S924;
        if(_S975)
        {
            if(is_fisheye_9)
            {
                float _S976 = length_1(_s_diff_ctx_3->_S923);
                float3  raydir_29 = make_float3 ((_s_diff_ctx_3->_S923 / make_float2 (s_primal_ctx_max_0(_S976, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S976))).x, (_s_diff_ctx_3->_S923 / make_float2 (s_primal_ctx_max_0(_S976, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S976))).y, s_primal_ctx_cos_0(_S976));
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
                float3  raydir_30 = make_float3 (_s_diff_ctx_3->_S923.x, _s_diff_ctx_3->_S923.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_24 = normalize_0(raydir_30);
                }
                else
                {
                    raydir_24 = raydir_30;
                }
            }
            float3  _S977 = make_float3 (_S968.primal_0.y) * raydir_24;
            _S971 = int(2);
            points_2[int(0)] = _S974;
            points_2[int(1)] = _S977;
            points_2[int(2)] = _S969;
            points_2[int(3)] = _S969;
        }
        else
        {
            _S971 = int(0);
            points_2[int(0)] = _S974;
            points_2[int(1)] = _S969;
            points_2[int(2)] = _S969;
            points_2[int(3)] = _S969;
            raydir_24 = _S969;
        }
        if(_S971 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S970;
            _S971 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_3->_S926)
            {
                _runFlag_2 = false;
                _S971 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_9)
                {
                    float _S978 = length_1(_s_diff_ctx_3->_S925);
                    float3  raydir_31 = make_float3 ((_s_diff_ctx_3->_S925 / make_float2 (s_primal_ctx_max_0(_S978, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S978))).x, (_s_diff_ctx_3->_S925 / make_float2 (s_primal_ctx_max_0(_S978, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S978))).y, s_primal_ctx_cos_0(_S978));
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
                    float3  raydir_32 = make_float3 (_s_diff_ctx_3->_S925.x, _s_diff_ctx_3->_S925.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_25 = normalize_0(raydir_32);
                    }
                    else
                    {
                        raydir_25 = raydir_32;
                    }
                }
                points_2[int(2)] = make_float3 (_S968.primal_0.z) * raydir_25;
                bool _S979 = !!_s_diff_ctx_3->_S928;
                int _S980;
                if(_S979)
                {
                    if(is_fisheye_9)
                    {
                        float _S981 = length_1(_s_diff_ctx_3->_S927);
                        float3  raydir_33 = make_float3 ((_s_diff_ctx_3->_S927 / make_float2 (s_primal_ctx_max_0(_S981, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S981))).x, (_s_diff_ctx_3->_S927 / make_float2 (s_primal_ctx_max_0(_S981, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S981))).y, s_primal_ctx_cos_0(_S981));
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
                        float3  raydir_34 = make_float3 (_s_diff_ctx_3->_S927.x, _s_diff_ctx_3->_S927.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_26 = normalize_0(raydir_34);
                        }
                        else
                        {
                            raydir_26 = raydir_34;
                        }
                    }
                    points_2[int(3)] = make_float3 (_S968.primal_0.w) * raydir_26;
                    _S980 = int(2);
                }
                else
                {
                    _S980 = int(0);
                    raydir_26 = _S969;
                }
                if(_S980 != int(2))
                {
                    _runFlag_3 = false;
                    _S971 = _S980;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S971 = int(1);
                }
                float3  _S982 = raydir_25;
                _runFlag_3 = _S979;
                raydir_25 = raydir_26;
                raydir_26 = _S982;
            }
            else
            {
                _runFlag_3 = false;
                raydir_25 = _S969;
                raydir_26 = _S969;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_25 = _S969;
            raydir_26 = _S969;
        }
        float3  _S983 = raydir_23;
        float3  _S984 = raydir_24;
        raydir_23 = raydir_25;
        raydir_24 = raydir_26;
        _S972 = _S975;
        raydir_25 = _S984;
        raydir_26 = _S983;
    }
    else
    {
        _S971 = int(0);
        points_2[int(0)] = _S969;
        points_2[int(1)] = _S969;
        points_2[int(2)] = _S969;
        points_2[int(3)] = _S969;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_23 = _S969;
        raydir_24 = _S969;
        _S972 = false;
        raydir_25 = _S969;
        raydir_26 = _S969;
    }
    bool _S985 = !(_S971 != int(1));
    float3  _S986;
    float3  _S987;
    float3  _S988;
    float3  _S989;
    float3  _S990;
    bool _S991;
    if(_S985)
    {
        float3  dx_0 = points_2[int(1)] - points_2[int(0)];
        float3  _S992 = - (points_2[int(3)] - points_2[int(2)]);
        float3  _S993 = s_primal_ctx_cross_0(dx_0, _S992);
        float _S994 = length_2(_S993);
        float3  _S995 = make_float3 (_S994);
        bool _S996 = _S994 != 0.0f;
        if(_S996)
        {
            _S986 = make_float3 (_S994 * _S994);
        }
        else
        {
            _S986 = _S969;
        }
        _S991 = _S996;
        _S987 = _S993;
        _S988 = _S995;
        _S989 = dx_0;
        _S990 = _S992;
    }
    else
    {
        _S991 = false;
        _S986 = _S969;
        _S987 = _S969;
        _S988 = _S969;
        _S989 = _S969;
        _S990 = _S969;
    }
    float4  _S997 = make_float4 (0.0f);
    if(_S985)
    {
        if(_S991)
        {
            float3  _S998 = dpnormal_1 / _S986;
            float3  _S999 = _S988 * _S998;
            _S986 = _S987 * - _S998;
            _S988 = _S999;
        }
        else
        {
            _S986 = _S969;
            _S988 = dpnormal_1;
        }
        float _S1000 = _S986.x + _S986.y + _S986.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1001;
        (&_S1001)->primal_0 = _S987;
        (&_S1001)->differential_0 = _S969;
        s_bwd_length_impl_1(&_S1001, _S1000);
        float3  _S1002 = _S1001.differential_0 + _S988;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1003;
        (&_S1003)->primal_0 = _S989;
        (&_S1003)->differential_0 = _S969;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1004;
        (&_S1004)->primal_0 = _S990;
        (&_S1004)->differential_0 = _S969;
        s_bwd_prop_cross_0(&_S1003, &_S1004, _S1002);
        float3  s_diff_dy_T_0 = - _S1004.differential_0;
        float3  _S1005 = - s_diff_dy_T_0;
        float3  _S1006 = - _S1003.differential_0;
        FixedArray<float3 , 4>  _S1007;
        _S1007[int(0)] = _S969;
        _S1007[int(1)] = _S969;
        _S1007[int(2)] = _S969;
        _S1007[int(3)] = _S969;
        _S1007[int(2)] = _S1005;
        _S1007[int(3)] = s_diff_dy_T_0;
        _S1007[int(0)] = _S1006;
        _S1007[int(1)] = _S1003.differential_0;
        points_2[int(0)] = _S1007[int(0)];
        points_2[int(1)] = _S1007[int(1)];
        points_2[int(2)] = _S1007[int(2)];
        points_2[int(3)] = _S1007[int(3)];
    }
    else
    {
        points_2[int(0)] = _S969;
        points_2[int(1)] = _S969;
        points_2[int(2)] = _S969;
        points_2[int(3)] = _S969;
    }
    float4  _S1008;
    if(_S970)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1009 = points_2;
                FixedArray<float3 , 4>  _S1010 = points_2;
                FixedArray<float3 , 4>  _S1011 = points_2;
                FixedArray<float3 , 4>  _S1012 = points_2;
                if(_runFlag_3)
                {
                    float3  _S1013 = raydir_23 * _S1012[int(3)];
                    float _S1014 = _S1013.x + _S1013.y + _S1013.z;
                    float4  _S1015 = _S997;
                    *&((&_S1015)->w) = _S1014;
                    points_2[int(0)] = _S1009[int(0)];
                    points_2[int(1)] = _S1010[int(1)];
                    points_2[int(2)] = _S1011[int(2)];
                    points_2[int(3)] = _S969;
                    _S1008 = _S1015;
                }
                else
                {
                    points_2[int(0)] = _S1009[int(0)];
                    points_2[int(1)] = _S1010[int(1)];
                    points_2[int(2)] = _S1011[int(2)];
                    points_2[int(3)] = _S1012[int(3)];
                    _S1008 = _S997;
                }
                float3  _S1016 = raydir_24 * points_2[int(2)];
                float _S1017 = _S1016.x + _S1016.y + _S1016.z;
                FixedArray<float3 , 4>  _S1018 = points_2;
                FixedArray<float3 , 4>  _S1019 = points_2;
                float4  _S1020 = _S997;
                *&((&_S1020)->z) = _S1017;
                float4  _S1021 = _S1008 + _S1020;
                points_2[int(0)] = points_2[int(0)];
                points_2[int(1)] = _S1018[int(1)];
                points_2[int(2)] = _S969;
                points_2[int(3)] = _S1019[int(3)];
                _S1008 = _S1021;
            }
            else
            {
                FixedArray<float3 , 4>  _S1022 = points_2;
                FixedArray<float3 , 4>  _S1023 = points_2;
                FixedArray<float3 , 4>  _S1024 = points_2;
                points_2[int(0)] = points_2[int(0)];
                points_2[int(1)] = _S1022[int(1)];
                points_2[int(2)] = _S1023[int(2)];
                points_2[int(3)] = _S1024[int(3)];
                _S1008 = _S997;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S1025 = points_2;
            FixedArray<float3 , 4>  _S1026 = points_2;
            FixedArray<float3 , 4>  _S1027 = points_2;
            points_2[int(0)] = points_2[int(0)];
            points_2[int(1)] = _S1025[int(1)];
            points_2[int(2)] = _S1026[int(2)];
            points_2[int(3)] = _S1027[int(3)];
            _S1008 = _S997;
        }
        if(_S972)
        {
            FixedArray<float3 , 4>  _S1028 = points_2;
            float3  _S1029 = raydir_25 * points_2[int(1)];
            float _S1030 = _S1029.x + _S1029.y + _S1029.z;
            float4  _S1031 = _S997;
            *&((&_S1031)->y) = _S1030;
            float4  _S1032 = _S1008 + _S1031;
            points_2[int(0)] = _S969;
            points_2[int(1)] = _S969;
            points_2[int(2)] = _S969;
            points_2[int(3)] = _S969;
            raydir_23 = _S1028[int(0)];
            _S1008 = _S1032;
        }
        else
        {
            FixedArray<float3 , 4>  _S1033 = points_2;
            FixedArray<float3 , 4>  _S1034 = points_2;
            FixedArray<float3 , 4>  _S1035 = points_2;
            points_2[int(0)] = points_2[int(0)];
            points_2[int(1)] = _S1033[int(1)];
            points_2[int(2)] = _S1034[int(2)];
            points_2[int(3)] = _S1035[int(3)];
            raydir_23 = _S969;
        }
        float3  _S1036 = raydir_26 * (points_2[int(0)] + raydir_23);
        float _S1037 = _S1036.x + _S1036.y + _S1036.z;
        float4  _S1038 = _S997;
        *&((&_S1038)->x) = _S1037;
        _S1008 = _S1008 + _S1038;
    }
    else
    {
        _S1008 = _S997;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1008;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(uint _S1039, uint _S1040, float2  _S1041, float4  _S1042, FixedArray<float, 10>  * _S1043, bool _S1044, bool _S1045, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1046, float3  _S1047)
{
    float3  _S1048;
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1049;
    s_primal_ctx_depth_to_normal_0(_S1039, _S1040, _S1041, _S1042, _S1043, _S1044, _S1045, (*_S1046).primal_0, &_S1048, &_S1049);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1050 = _S1049;
    s_bwd_prop_depth_to_normal_0(_S1039, _S1040, _S1041, _S1042, _S1043, _S1044, _S1045, _S1046, _S1047, &_S1050);
    return;
}

inline __device__ void depth_to_normal_vjp(uint width_3, uint height_3, float2  pix_center_3, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_13, bool is_fisheye_10, bool is_ray_depth_3, float4  depths_1, float3  v_normal_0, float4  * v_depths_0)
{
    float4  _S1051 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1051;
    s_bwd_depth_to_normal_0(width_3, height_3, pix_center_3, intrins_3, dist_coeffs_13, is_fisheye_10, is_ray_depth_3, &dp_depths_0, v_normal_0);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(uint width_4, uint height_4, float2  pix_center_4, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_14, bool is_fisheye_11)
{
    float2  _S1052 = (pix_center_4 - float2 {intrins_4.z, intrins_4.w}) / float2 {intrins_4.x, intrins_4.y};
    float2  uv_18 = _S1052;
    bool _S1053 = undistort_point_0(_S1052, dist_coeffs_14, int(12), &uv_18);
    if(!_S1053)
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

