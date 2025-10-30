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

inline __device__ float dot_0(float4  x_0, float4  y_0)
{
    int i_0 = int(0);
    float result_0 = 0.0f;
    for(;;)
    {
        if(i_0 < int(4))
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

inline __device__ float dot_1(float2  x_1, float2  y_1)
{
    int i_1 = int(0);
    float result_2 = 0.0f;
    for(;;)
    {
        if(i_1 < int(2))
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

inline __device__ float dot_2(float3  x_2, float3  y_2)
{
    int i_2 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_2 < int(3))
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
    return (F32_sqrt((dot_0(x_3, x_3))));
}

inline __device__ float length_1(float2  x_4)
{
    return (F32_sqrt((dot_1(x_4, x_4))));
}

inline __device__ float length_2(float3  x_5)
{
    return (F32_sqrt((dot_2(x_5, x_5))));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_3, float dOut_3)
{
    float _S6 = 1.0f / (*dpx_3).primal_0 * dOut_3;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S6;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

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

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, float3  dOut_4)
{
    float3  _S7 = exp_0((*dpx_4).primal_0) * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S7;
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

inline __device__ void _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_5, float2  dOut_5)
{
    float2  _S8 = exp_1((*dpx_5).primal_0) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S8;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_6, DiffPair_float_0 * dpy_1, float dOut_6)
{
    DiffPair_float_0 _S9 = *dpx_6;
    float _S10;
    if(((*dpx_6).primal_0) < ((*dpy_1).primal_0))
    {
        _S10 = dOut_6;
    }
    else
    {
        if(((*dpx_6).primal_0) > ((*dpy_1).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_6;
        }
    }
    dpx_6->primal_0 = _S9.primal_0;
    dpx_6->differential_0 = _S10;
    DiffPair_float_0 _S11 = *dpy_1;
    if(((*dpy_1).primal_0) < (_S9.primal_0))
    {
        _S10 = dOut_6;
    }
    else
    {
        if(((*dpy_1).primal_0) > ((*dpx_6).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_6;
        }
    }
    dpy_1->primal_0 = _S11.primal_0;
    dpy_1->differential_0 = _S10;
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

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_7, float _s_dOut_0)
{
    float _S47 = (*dpx_7).primal_0.x;
    float _S48 = (*dpx_7).primal_0.y;
    float _S49 = (*dpx_7).primal_0.z;
    float _S50 = (*dpx_7).primal_0.w;
    DiffPair_float_0 _S51;
    (&_S51)->primal_0 = _S47 * _S47 + _S48 * _S48 + _S49 * _S49 + _S50 * _S50;
    (&_S51)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S51, _s_dOut_0);
    float _S52 = (*dpx_7).primal_0.w * _S51.differential_0;
    float _S53 = _S52 + _S52;
    float _S54 = (*dpx_7).primal_0.z * _S51.differential_0;
    float _S55 = _S54 + _S54;
    float _S56 = (*dpx_7).primal_0.y * _S51.differential_0;
    float _S57 = _S56 + _S56;
    float _S58 = (*dpx_7).primal_0.x * _S51.differential_0;
    float _S59 = _S58 + _S58;
    float4  _S60 = make_float4 (0.0f);
    *&((&_S60)->w) = _S53;
    *&((&_S60)->z) = _S55;
    *&((&_S60)->y) = _S57;
    *&((&_S60)->x) = _S59;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S60;
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

inline __device__ float3  min_0(float3  x_12, float3  y_7)
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
        *_slang_vector_get_element_ptr(&result_8, i_5) = (F32_min((_slang_vector_get_element(x_12, i_5)), (_slang_vector_get_element(y_7, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_8;
}

inline __device__ float3  max_0(float3  x_13, float3  y_8)
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
        *_slang_vector_get_element_ptr(&result_9, i_6) = (F32_max((_slang_vector_get_element(x_13, i_6)), (_slang_vector_get_element(y_8, i_6))));
        i_6 = i_6 + int(1);
    }
    return result_9;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_7)
{
    DiffPair_float_0 _S282 = *dpx_8;
    bool _S283;
    if(((*dpx_8).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S283 = ((*dpx_8).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S283 = false;
    }
    float _S284;
    if(_S283)
    {
        _S284 = dOut_7;
    }
    else
    {
        _S284 = 0.0f;
    }
    dpx_8->primal_0 = _S282.primal_0;
    dpx_8->differential_0 = _S284;
    DiffPair_float_0 _S285 = *dpMin_0;
    if((_S282.primal_0) < ((*dpMin_0).primal_0))
    {
        _S284 = dOut_7;
    }
    else
    {
        _S284 = 0.0f;
    }
    dpMin_0->primal_0 = _S285.primal_0;
    dpMin_0->differential_0 = _S284;
    DiffPair_float_0 _S286 = *dpMax_0;
    if(((*dpx_8).primal_0) > ((*dpMax_0).primal_0))
    {
        _S284 = dOut_7;
    }
    else
    {
        _S284 = 0.0f;
    }
    dpMax_0->primal_0 = _S286.primal_0;
    dpMax_0->differential_0 = _S284;
    return;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_8)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_9).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_2).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_8.x);
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_9).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_2).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_8.y);
    *&((&left_d_result_0)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_0)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_9).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_2).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_8.z);
    *&((&left_d_result_0)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_0)->z) = right_dp_2.differential_0;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = left_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  clamp_0(float3  x_14, float3  minBound_0, float3  maxBound_0)
{
    return min_0(max_0(x_14, minBound_0), maxBound_0);
}

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_0(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S287, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S288, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S289, float3  _S290)
{
    _d_clamp_vector_0(_S287, _S288, _S289, _S290);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_2)
{
    float _S291 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S292 = make_float3 (_S291);
    float3  _S293 = make_float3 (0.0f);
    float3  _S294 = make_float3 (1.0f);
    float3  _S295 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S296;
    (&_S296)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S291) * (*dpbackground_0).primal_0;
    (&_S296)->differential_0 = _S295;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S297;
    (&_S297)->primal_0 = _S293;
    (&_S297)->differential_0 = _S295;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S298;
    (&_S298)->primal_0 = _S294;
    (&_S298)->differential_0 = _S295;
    s_bwd_prop_clamp_0(&_S296, &_S297, &_S298, _s_dOut_2);
    float3  _S299 = _S292 * _S296.differential_0;
    float3  _S300 = (*dpbackground_0).primal_0 * _S296.differential_0;
    float _S301 = - (_S300.x + _S300.y + _S300.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S299;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S301;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S296.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S302, DiffPair_float_0 * _S303, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S304, float3  _S305)
{
    s_bwd_prop_blend_background_0(_S302, _S303, _S304, _S305);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S306 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S306;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S306;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_15)
{
    Matrix<float, 3, 3>  result_10;
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
            *_slang_vector_get_element_ptr(((&result_10)->rows + (r_0)), c_0) = _slang_vector_get_element(x_15.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_10;
}

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_2)
{
    float x_16 = quat_2.y;
    float x2_0 = x_16 * x_16;
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
    float3  result_11;
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
            float sum_1 = sum_0 + _slang_vector_get_element(left_0.rows[i_7], j_0) * _slang_vector_get_element(right_0, j_0);
            j_0 = j_0 + int(1);
            sum_0 = sum_1;
        }
        *_slang_vector_get_element_ptr(&result_11, i_7) = sum_0;
        i_7 = i_7 + int(1);
    }
    return result_11;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_1, Matrix<float, 3, 3>  right_1)
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
            int i_8 = int(0);
            float sum_2 = 0.0f;
            for(;;)
            {
                if(i_8 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_3 = sum_2 + _slang_vector_get_element(left_1.rows[r_1], i_8) * _slang_vector_get_element(right_1.rows[i_8], c_1);
                i_8 = i_8 + int(1);
                sum_2 = sum_3;
            }
            *_slang_vector_get_element_ptr(((&result_12)->rows + (r_1)), c_1) = sum_2;
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_12;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_1(mul_1(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_3, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_17 = quat_3.y;
    float x2_1 = x_17 * x_17;
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
    CameraDistortion_0 _S307;
    (&_S307)->radial_coeffs_0 = radial_coeffs_1;
    (&_S307)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S307)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    return _S307;
}

inline __device__ float2  distort_point(float2  uv_0, bool is_fisheye_0, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2)
{
    float2  _S308;
    if(is_fisheye_0)
    {
        float r_2 = length_1(uv_0);
        float theta_0 = (F32_atan((r_2)));
        float _S309;
        if(r_2 < 0.00100000004749745f)
        {
            _S309 = 1.0f - r_2 * r_2 / 3.0f;
        }
        else
        {
            _S309 = theta_0 / r_2;
        }
        _S308 = uv_0 * make_float2 (_S309);
    }
    else
    {
        _S308 = uv_0;
    }
    CameraDistortion_0 _S310 = CameraDistortion_x24init_0(radial_coeffs_2, tangential_coeffs_2, thin_prism_coeffs_2);
    float p1_0 = _S310.tangential_coeffs_0.x;
    float p2_0 = _S310.tangential_coeffs_0.y;
    float u_0 = _S308.x;
    float v_0 = _S308.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    return _S308 * make_float2 (1.0f + r2_0 * (_S310.radial_coeffs_0.x + r2_0 * (_S310.radial_coeffs_0.y + r2_0 * (_S310.radial_coeffs_0.z + r2_0 * _S310.radial_coeffs_0.w)))) + make_float2 (2.0f * p1_0 * u_0 * v_0 + p2_0 * (r2_0 + 2.0f * u_0 * u_0) + _S310.thin_prism_coeffs_0.x * r2_0, 2.0f * p2_0 * u_0 * v_0 + p1_0 * (r2_0 + 2.0f * v_0 * v_0) + _S310.thin_prism_coeffs_0.y * r2_0);
}

inline __device__ float2  undistort_point_0(float2  uv_1, CameraDistortion_0 * dist_coeffs_0, int maxiter_0)
{
    int i_9 = int(0);
    float2  q_0 = uv_1;
    for(;;)
    {
        if(i_9 < maxiter_0)
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
        float _S311 = dist_coeffs_0->radial_coeffs_0.z + r2_1 * k4_0;
        float _S312 = dist_coeffs_0->radial_coeffs_0.y + r2_1 * _S311;
        float _S313 = dist_coeffs_0->radial_coeffs_0.x + r2_1 * _S312;
        float radial_0 = 1.0f + r2_1 * _S313;
        float _S314 = 2.0f * p1_1;
        float _S315 = _S314 * u_1;
        float _S316 = 2.0f * u_1;
        float _S317 = 2.0f * p2_1;
        float _S318 = _S317 * u_1;
        float _S319 = 2.0f * v_1;
        float2  _S320 = q_0 * make_float2 (radial_0) + make_float2 (_S315 * v_1 + p2_1 * (r2_1 + _S316 * u_1) + sx1_0 * r2_1, _S318 * v_1 + p1_1 * (r2_1 + _S319 * v_1) + sy1_0 * r2_1);
        float2  _S321 = make_float2 (0.0f);
        float2  seed_0 = _S321;
        *&((&seed_0)->x) = 1.0f;
        float2  _S322 = make_float2 (radial_0);
        float2  _S323 = q_0 * seed_0;
        float _S324 = p1_1 * seed_0.y;
        float _S325 = p2_1 * seed_0.x;
        float _S326 = _S323.x + _S323.y;
        float _S327 = r2_1 * _S326;
        float _S328 = r2_1 * _S327;
        float _S329 = sy1_0 * seed_0.y + _S324 + sx1_0 * seed_0.x + _S325 + _S313 * _S326 + _S312 * _S327 + _S311 * _S328 + k4_0 * (r2_1 * _S328);
        float _S330 = v_1 * _S329;
        float _S331 = u_1 * _S329;
        Matrix<float, 2, 2>  J_0;
        J_0[int(0)] = _S322 * seed_0 + make_float2 (_S317 * (v_1 * seed_0.y) + _S316 * _S325 + 2.0f * (u_1 * _S325) + _S314 * (v_1 * seed_0.x) + _S331 + _S331, _S319 * _S324 + 2.0f * (v_1 * _S324) + _S318 * seed_0.y + _S315 * seed_0.x + _S330 + _S330);
        float2  seed_1 = _S321;
        *&((&seed_1)->y) = 1.0f;
        float2  _S332 = q_0 * seed_1;
        float _S333 = p1_1 * seed_1.y;
        float _S334 = p2_1 * seed_1.x;
        float _S335 = _S332.x + _S332.y;
        float _S336 = r2_1 * _S335;
        float _S337 = r2_1 * _S336;
        float _S338 = sy1_0 * seed_1.y + _S333 + sx1_0 * seed_1.x + _S334 + _S313 * _S335 + _S312 * _S336 + _S311 * _S337 + k4_0 * (r2_1 * _S337);
        float _S339 = v_1 * _S338;
        float _S340 = u_1 * _S338;
        J_0[int(1)] = _S322 * seed_1 + make_float2 (_S317 * (v_1 * seed_1.y) + _S316 * _S334 + 2.0f * (u_1 * _S334) + _S314 * (v_1 * seed_1.x) + _S340 + _S340, _S319 * _S333 + 2.0f * (v_1 * _S333) + _S318 * seed_1.y + _S315 * seed_1.x + _S339 + _S339);
        float2  _S341 = _S320 - uv_1;
        float inv_det_0 = 1.0f / (J_0.rows[int(0)].x * J_0.rows[int(1)].y - J_0.rows[int(0)].y * J_0.rows[int(1)].x);
        float _S342 = _S341.x;
        float _S343 = _S341.y;
        float2  q_1 = q_0 - make_float2 ((_S342 * J_0.rows[int(1)].y - _S343 * J_0.rows[int(0)].y) * inv_det_0, (- _S342 * J_0.rows[int(1)].x + _S343 * J_0.rows[int(0)].x) * inv_det_0);
        i_9 = i_9 + int(1);
        q_0 = q_1;
    }
    return q_0;
}

inline __device__ float2  undistort_point(float2  uv_2, bool is_fisheye_1, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3)
{
    CameraDistortion_0 _S344;
    (&_S344)->radial_coeffs_0 = radial_coeffs_3;
    (&_S344)->tangential_coeffs_0 = tangential_coeffs_3;
    (&_S344)->thin_prism_coeffs_0 = thin_prism_coeffs_3;
    float2  _S345 = undistort_point_0(uv_2, &_S344, int(8));
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float theta_1 = length_1(_S345);
        float _S346;
        if(theta_1 < 0.00100000004749745f)
        {
            _S346 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S346 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S347 = make_float3 ((_S345 * make_float2 (_S346)).x, (_S345 * make_float2 (_S346)).y, (F32_cos((theta_1))));
        raydir_0 = _S347;
    }
    else
    {
        raydir_0 = make_float3 (_S345.x, _S345.y, 1.0f);
    }
    return float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_9)
{
    float _S348 = (*right_2).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_9.x;
    float sum_4 = _S348 + (*right_2).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_9.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_9.z;
    float3  left_d_result_1;
    *&((&left_d_result_1)->x) = sum_5;
    float _S349 = (*right_2).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_9.x;
    float sum_6 = _S349 + (*right_2).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_9.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_9.z;
    *&((&left_d_result_1)->y) = sum_7;
    float _S350 = (*right_2).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_9.x;
    float sum_8 = _S350 + (*right_2).primal_0.rows[int(2)].y * dOut_9.y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = (*left_2).primal_0.z * dOut_9.y;
    float sum_9 = sum_8 + (*right_2).primal_0.rows[int(2)].z * dOut_9.z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = (*left_2).primal_0.z * dOut_9.z;
    *&((&left_d_result_1)->z) = sum_9;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_2(float3  left_3, Matrix<float, 3, 3>  right_3)
{
    float3  result_13;
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
        int i_10 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(i_10 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_3, i_10) * _slang_vector_get_element(right_3.rows[i_10], j_1);
            i_10 = i_10 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_13, j_1) = sum_10;
        j_1 = j_1 + int(1);
    }
    return result_13;
}

inline __device__ float3  normalize_0(float3  x_19)
{
    return x_19 / make_float3 (length_2(x_19));
}

inline __device__ void generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_3, bool is_fisheye_2, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, float3  * ray_o_0, float3  * ray_d_0)
{
    *ray_o_0 = - mul_2(t_1, R_2);
    CameraDistortion_0 _S351;
    (&_S351)->radial_coeffs_0 = radial_coeffs_4;
    (&_S351)->tangential_coeffs_0 = tangential_coeffs_4;
    (&_S351)->thin_prism_coeffs_0 = thin_prism_coeffs_4;
    float2  _S352 = undistort_point_0(uv_3, &_S351, int(8));
    float3  raydir_1;
    if(is_fisheye_2)
    {
        float theta_2 = length_1(_S352);
        float _S353;
        if(theta_2 < 0.00100000004749745f)
        {
            _S353 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S353 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S354 = make_float3 ((_S352 * make_float2 (_S353)).x, (_S352 * make_float2 (_S353)).y, (F32_cos((theta_2))));
        raydir_1 = _S354;
    }
    else
    {
        raydir_1 = make_float3 (_S352.x, _S352.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_2(raydir_1, R_2));
    return;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S355;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S356, Matrix<float, 3, 3>  _S357)
{
    return mul_2(_S356, _S357);
}

inline __device__ float s_primal_ctx_sin_0(float _S358)
{
    return (F32_sin((_S358)));
}

inline __device__ float s_primal_ctx_cos_0(float _S359)
{
    return (F32_cos((_S359)));
}

inline __device__ void s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_4, bool is_fisheye_3, float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S355 = make_float2 (0.0f);
    float3  _S360 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    CameraDistortion_0 _S361;
    (&_S361)->radial_coeffs_0 = radial_coeffs_5;
    (&_S361)->tangential_coeffs_0 = tangential_coeffs_5;
    (&_S361)->thin_prism_coeffs_0 = thin_prism_coeffs_5;
    float2  _S362 = undistort_point_0(uv_4, &_S361, int(8));
    _s_diff_ctx_0->_S355 = _S362;
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float _S363 = length_1(_S362);
        float _S364;
        if(_S363 < 0.00100000004749745f)
        {
            _S364 = 1.0f - _S363 * _S363 / 6.0f;
        }
        else
        {
            _S364 = s_primal_ctx_sin_0(_S363) / _S363;
        }
        float3  _S365 = make_float3 ((_S362 * make_float2 (_S364)).x, (_S362 * make_float2 (_S364)).y, s_primal_ctx_cos_0(_S363));
        raydir_2 = _S365;
    }
    else
    {
        raydir_2 = make_float3 (_S362.x, _S362.y, 1.0f);
    }
    float3  _S366 = normalize_0(s_primal_ctx_mul_0(raydir_2, dpR_0));
    *dpray_o_0 = _S360;
    *dpray_d_0 = _S366;
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float _s_dOut_3)
{
    float _S367 = (*dpx_10).primal_0.x;
    float _S368 = (*dpx_10).primal_0.y;
    float _S369 = (*dpx_10).primal_0.z;
    DiffPair_float_0 _S370;
    (&_S370)->primal_0 = _S367 * _S367 + _S368 * _S368 + _S369 * _S369;
    (&_S370)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S370, _s_dOut_3);
    float _S371 = (*dpx_10).primal_0.z * _S370.differential_0;
    float _S372 = _S371 + _S371;
    float _S373 = (*dpx_10).primal_0.y * _S370.differential_0;
    float _S374 = _S373 + _S373;
    float _S375 = (*dpx_10).primal_0.x * _S370.differential_0;
    float _S376 = _S375 + _S375;
    float3  _S377 = make_float3 (0.0f);
    *&((&_S377)->z) = _S372;
    *&((&_S377)->y) = _S374;
    *&((&_S377)->x) = _S376;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S377;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S378, float _S379)
{
    s_bwd_prop_length_impl_1(_S378, _S379);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, float3  _s_dOut_4)
{
    float _S380 = length_2((*dpx_11).primal_0);
    float3  _S381 = (*dpx_11).primal_0 * _s_dOut_4;
    float3  _S382 = make_float3 (1.0f / _S380) * _s_dOut_4;
    float _S383 = - ((_S381.x + _S381.y + _S381.z) / (_S380 * _S380));
    float3  _S384 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S385;
    (&_S385)->primal_0 = (*dpx_11).primal_0;
    (&_S385)->differential_0 = _S384;
    s_bwd_length_impl_1(&_S385, _S383);
    float3  _S386 = _S382 + _S385.differential_0;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S386;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S387, float3  _S388)
{
    s_bwd_prop_normalize_impl_0(_S387, _S388);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S389, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S390, float3  _S391)
{
    _d_mul_0(_S389, _S390, _S391);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_5, bool is_fisheye_4, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S392 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S393 = *dpt_1;
    float3  raydir_3;
    if(is_fisheye_4)
    {
        float _S394 = length_1(_s_diff_ctx_1->_S355);
        float _S395;
        if(_S394 < 0.00100000004749745f)
        {
            _S395 = 1.0f - _S394 * _S394 / 6.0f;
        }
        else
        {
            _S395 = s_primal_ctx_sin_0(_S394) / _S394;
        }
        float3  _S396 = make_float3 ((_s_diff_ctx_1->_S355 * make_float2 (_S395)).x, (_s_diff_ctx_1->_S355 * make_float2 (_S395)).y, s_primal_ctx_cos_0(_S394));
        raydir_3 = _S396;
    }
    else
    {
        raydir_3 = make_float3 (_s_diff_ctx_1->_S355.x, _s_diff_ctx_1->_S355.y, 1.0f);
    }
    float3  _S397 = s_primal_ctx_mul_0(raydir_3, _S392.primal_0);
    float3  _S398 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S399;
    (&_S399)->primal_0 = _S397;
    (&_S399)->differential_0 = _S398;
    s_bwd_normalize_impl_0(&_S399, dpray_d_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S400;
    (&_S400)->primal_0 = raydir_3;
    (&_S400)->differential_0 = _S398;
    Matrix<float, 3, 3>  _S401 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S402;
    (&_S402)->primal_0 = _S392.primal_0;
    (&_S402)->differential_0 = _S401;
    s_bwd_prop_mul_0(&_S400, &_S402, _S399.differential_0);
    float3  _S403 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S404;
    (&_S404)->primal_0 = _S393.primal_0;
    (&_S404)->differential_0 = _S398;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S405;
    (&_S405)->primal_0 = _S392.primal_0;
    (&_S405)->differential_0 = _S401;
    s_bwd_prop_mul_0(&_S404, &_S405, _S403);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S404.differential_0;
    Matrix<float, 3, 3>  _S406 = _S405.differential_0 + _S402.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S406;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S407, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S408, float2  _S409, bool _S410, float4  _S411, float2  _S412, float2  _S413, float3  _S414, float3  _S415)
{
    float3  _S416;
    float3  _S417;
    s_bwd_prop_generate_ray_Intermediates_0 _S418;
    s_primal_ctx_generate_ray_0((*_S407).primal_0, (*_S408).primal_0, _S409, _S410, _S411, _S412, _S413, &_S416, &_S417, &_S418);
    s_bwd_prop_generate_ray_Intermediates_0 _S419 = _S418;
    s_bwd_prop_generate_ray_0(_S407, _S408, _S409, _S410, _S411, _S412, _S413, _S414, _S415, &_S419);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_6, bool is_fisheye_5, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S420 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S420;
    float3  _S421 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S421;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_6, is_fisheye_5, radial_coeffs_7, tangential_coeffs_7, thin_prism_coeffs_7, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_10)
{
    float _S422 = dOut_10.y;
    float _S423 = dOut_10.z;
    float _S424 = dOut_10.x;
    float _S425 = (*a_0).primal_0.z * _S422 + - (*a_0).primal_0.y * _S423;
    float _S426 = - (*a_0).primal_0.z * _S424 + (*a_0).primal_0.x * _S423;
    float _S427 = (*a_0).primal_0.y * _S424 + - (*a_0).primal_0.x * _S422;
    float3  _S428 = make_float3 (- (*b_0).primal_0.z * _S422 + (*b_0).primal_0.y * _S423, (*b_0).primal_0.z * _S424 + - (*b_0).primal_0.x * _S423, - (*b_0).primal_0.y * _S424 + (*b_0).primal_0.x * _S422);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S428;
    float3  _S429 = make_float3 (_S425, _S426, _S427);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S429;
    return;
}

inline __device__ float3  cross_0(float3  left_4, float3  right_4)
{
    float _S430 = left_4.y;
    float _S431 = right_4.z;
    float _S432 = left_4.z;
    float _S433 = right_4.y;
    float _S434 = right_4.x;
    float _S435 = left_4.x;
    return make_float3 (_S430 * _S431 - _S432 * _S433, _S432 * _S434 - _S435 * _S431, _S435 * _S433 - _S430 * _S434);
}

inline __device__ void depth_to_normal(uint width_0, uint height_0, float2  pix_center_0, float4  intrins_0, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, bool is_fisheye_6, bool is_ray_depth_0, float4  depths_0, float3  * normal_0)
{
    FixedArray<float3 , 4>  points_0;
    float2  _S436 = float2 {intrins_0.z, intrins_0.w};
    float2  _S437 = float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S436) / _S437;
    float3  raydir_4;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S438;
        (&_S438)->radial_coeffs_0 = radial_coeffs_8;
        (&_S438)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S438)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S439 = undistort_point_0(uv_7, &_S438, int(12));
        float theta_3 = length_1(_S439);
        float3  raydir_5 = make_float3 ((_S439 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).x, (_S439 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).y, (F32_cos((theta_3))));
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
    float2  uv_8 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S436) / _S437;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S440;
        (&_S440)->radial_coeffs_0 = radial_coeffs_8;
        (&_S440)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S440)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S441 = undistort_point_0(uv_8, &_S440, int(12));
        float theta_4 = length_1(_S441);
        float3  raydir_7 = make_float3 ((_S441 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (_S441 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
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
    float2  uv_9 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S436) / _S437;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S442;
        (&_S442)->radial_coeffs_0 = radial_coeffs_8;
        (&_S442)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S442)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S443 = undistort_point_0(uv_9, &_S442, int(12));
        float theta_5 = length_1(_S443);
        float3  raydir_9 = make_float3 ((_S443 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (_S443 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
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
    float2  uv_10 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S436) / _S437;
    if(is_fisheye_6)
    {
        CameraDistortion_0 _S444;
        (&_S444)->radial_coeffs_0 = radial_coeffs_8;
        (&_S444)->tangential_coeffs_0 = tangential_coeffs_8;
        (&_S444)->thin_prism_coeffs_0 = thin_prism_coeffs_8;
        float2  _S445 = undistort_point_0(uv_10, &_S444, int(12));
        float theta_6 = length_1(_S445);
        float3  raydir_11 = make_float3 ((_S445 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (_S445 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
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
    float3  _S446 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    *normal_0 = _S446;
    if((length_2(_S446)) != 0.0f)
    {
        *normal_0 = *normal_0 / make_float3 (length_2(*normal_0));
    }
    return;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S447;
    float2  _S448;
    float2  _S449;
    float2  _S450;
};

inline __device__ float3  s_primal_ctx_cross_0(float3  _S451, float3  _S452)
{
    return cross_0(_S451, _S452);
}

inline __device__ void s_primal_ctx_depth_to_normal_0(uint width_1, uint height_1, float2  pix_center_1, float4  intrins_1, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, bool is_fisheye_7, bool is_ray_depth_1, float4  dpdepths_0, float3  * dpnormal_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S453 = make_float2 (0.0f);
    _s_diff_ctx_2->_S447 = _S453;
    _s_diff_ctx_2->_S448 = _S453;
    _s_diff_ctx_2->_S449 = _S453;
    _s_diff_ctx_2->_S450 = _S453;
    _s_diff_ctx_2->_S447 = _S453;
    _s_diff_ctx_2->_S448 = _S453;
    _s_diff_ctx_2->_S449 = _S453;
    _s_diff_ctx_2->_S450 = _S453;
    float2  _S454 = float2 {intrins_1.z, intrins_1.w};
    float2  _S455 = float2 {intrins_1.x, intrins_1.y};
    float2  uv_11 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S454) / _S455;
    float3  raydir_13;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S456;
        (&_S456)->radial_coeffs_0 = radial_coeffs_9;
        (&_S456)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S456)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S457 = undistort_point_0(uv_11, &_S456, int(12));
        _s_diff_ctx_2->_S447 = _S457;
        float _S458 = length_1(_S457);
        float3  raydir_14 = make_float3 ((_S457 / make_float2 (s_primal_ctx_max_0(_S458, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S458))).x, (_S457 / make_float2 (s_primal_ctx_max_0(_S458, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S458))).y, s_primal_ctx_cos_0(_S458));
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
    float3  _S459 = make_float3 (dpdepths_0.x) * raydir_13;
    float2  uv_12 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S454) / _S455;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S460;
        (&_S460)->radial_coeffs_0 = radial_coeffs_9;
        (&_S460)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S460)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S461 = undistort_point_0(uv_12, &_S460, int(12));
        _s_diff_ctx_2->_S448 = _S461;
        float _S462 = length_1(_S461);
        float3  raydir_16 = make_float3 ((_S461 / make_float2 (s_primal_ctx_max_0(_S462, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S462))).x, (_S461 / make_float2 (s_primal_ctx_max_0(_S462, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S462))).y, s_primal_ctx_cos_0(_S462));
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
    float3  _S463 = make_float3 (dpdepths_0.y) * raydir_13;
    float2  uv_13 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S454) / _S455;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S464;
        (&_S464)->radial_coeffs_0 = radial_coeffs_9;
        (&_S464)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S464)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S465 = undistort_point_0(uv_13, &_S464, int(12));
        _s_diff_ctx_2->_S449 = _S465;
        float _S466 = length_1(_S465);
        float3  raydir_18 = make_float3 ((_S465 / make_float2 (s_primal_ctx_max_0(_S466, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S466))).x, (_S465 / make_float2 (s_primal_ctx_max_0(_S466, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S466))).y, s_primal_ctx_cos_0(_S466));
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
    float3  _S467 = make_float3 (dpdepths_0.z) * raydir_13;
    float2  uv_14 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S454) / _S455;
    if(is_fisheye_7)
    {
        CameraDistortion_0 _S468;
        (&_S468)->radial_coeffs_0 = radial_coeffs_9;
        (&_S468)->tangential_coeffs_0 = tangential_coeffs_9;
        (&_S468)->thin_prism_coeffs_0 = thin_prism_coeffs_9;
        float2  _S469 = undistort_point_0(uv_14, &_S468, int(12));
        _s_diff_ctx_2->_S450 = _S469;
        float _S470 = length_1(_S469);
        float3  raydir_20 = make_float3 ((_S469 / make_float2 (s_primal_ctx_max_0(_S470, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S470))).x, (_S469 / make_float2 (s_primal_ctx_max_0(_S470, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S470))).y, s_primal_ctx_cos_0(_S470));
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
    float3  _S471 = make_float3 (dpdepths_0.w) * raydir_13;
    float3  _S472 = s_primal_ctx_cross_0(_S463 - _S459, - (_S471 - _S467));
    float _S473 = length_2(_S472);
    if(_S473 != 0.0f)
    {
        raydir_13 = _S472 / make_float3 (_S473);
    }
    else
    {
        raydir_13 = _S472;
    }
    *dpnormal_0 = raydir_13;
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S474, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S475, float3  _S476)
{
    _d_cross_0(_S474, _S475, _S476);
    return;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(uint width_2, uint height_2, float2  pix_center_2, float4  intrins_2, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, bool is_fisheye_8, bool is_ray_depth_2, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  dpnormal_1, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S477 = *dpdepths_1;
    float3  _S478 = make_float3 (0.0f);
    float2  _S479 = float2 {intrins_2.z, intrins_2.w};
    float2  _S480 = float2 {intrins_2.x, intrins_2.y};
    float2  uv_15 = (pix_center_2 + make_float2 (-1.0f, -0.0f) - _S479) / _S480;
    float3  raydir_22;
    if(is_fisheye_8)
    {
        float _S481 = length_1(_s_diff_ctx_3->_S447);
        float3  raydir_23 = make_float3 ((_s_diff_ctx_3->_S447 / make_float2 (s_primal_ctx_max_0(_S481, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S481))).x, (_s_diff_ctx_3->_S447 / make_float2 (s_primal_ctx_max_0(_S481, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S481))).y, s_primal_ctx_cos_0(_S481));
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
    float3  _S482 = make_float3 (_S477.primal_0.x) * raydir_22;
    float2  uv_16 = (pix_center_2 + make_float2 (1.0f, -0.0f) - _S479) / _S480;
    float3  raydir_25;
    if(is_fisheye_8)
    {
        float _S483 = length_1(_s_diff_ctx_3->_S448);
        float3  raydir_26 = make_float3 ((_s_diff_ctx_3->_S448 / make_float2 (s_primal_ctx_max_0(_S483, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S483))).x, (_s_diff_ctx_3->_S448 / make_float2 (s_primal_ctx_max_0(_S483, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S483))).y, s_primal_ctx_cos_0(_S483));
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
    float3  _S484 = make_float3 (_S477.primal_0.y) * raydir_25;
    float2  uv_17 = (pix_center_2 + make_float2 (0.0f, -1.0f) - _S479) / _S480;
    float3  raydir_28;
    if(is_fisheye_8)
    {
        float _S485 = length_1(_s_diff_ctx_3->_S449);
        float3  raydir_29 = make_float3 ((_s_diff_ctx_3->_S449 / make_float2 (s_primal_ctx_max_0(_S485, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S485))).x, (_s_diff_ctx_3->_S449 / make_float2 (s_primal_ctx_max_0(_S485, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S485))).y, s_primal_ctx_cos_0(_S485));
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
    float3  _S486 = make_float3 (_S477.primal_0.z) * raydir_28;
    float2  uv_18 = (pix_center_2 + make_float2 (0.0f, 1.0f) - _S479) / _S480;
    float3  raydir_31;
    if(is_fisheye_8)
    {
        float _S487 = length_1(_s_diff_ctx_3->_S450);
        float3  raydir_32 = make_float3 ((_s_diff_ctx_3->_S450 / make_float2 (s_primal_ctx_max_0(_S487, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S487))).x, (_s_diff_ctx_3->_S450 / make_float2 (s_primal_ctx_max_0(_S487, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S487))).y, s_primal_ctx_cos_0(_S487));
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
    float3  _S488 = make_float3 (_S477.primal_0.w) * raydir_31;
    float3  dx_0 = _S484 - _S482;
    float3  _S489 = - (_S488 - _S486);
    float3  _S490 = s_primal_ctx_cross_0(dx_0, _S489);
    float _S491 = length_2(_S490);
    float3  _S492 = make_float3 (_S491);
    bool _S493 = _S491 != 0.0f;
    float3  _S494;
    if(_S493)
    {
        _S494 = make_float3 (_S491 * _S491);
    }
    else
    {
        _S494 = _S478;
    }
    float3  _S495;
    if(_S493)
    {
        float3  _S496 = dpnormal_1 / _S494;
        float3  _S497 = _S492 * _S496;
        _S494 = _S490 * - _S496;
        _S495 = _S497;
    }
    else
    {
        _S494 = _S478;
        _S495 = dpnormal_1;
    }
    float _S498 = _S494.x + _S494.y + _S494.z;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S499;
    (&_S499)->primal_0 = _S490;
    (&_S499)->differential_0 = _S478;
    s_bwd_length_impl_1(&_S499, _S498);
    float3  _S500 = _S499.differential_0 + _S495;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S501;
    (&_S501)->primal_0 = dx_0;
    (&_S501)->differential_0 = _S478;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S502;
    (&_S502)->primal_0 = _S489;
    (&_S502)->differential_0 = _S478;
    s_bwd_prop_cross_0(&_S501, &_S502, _S500);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S503 = _S501;
    float3  s_diff_dy_T_0 = - _S502.differential_0;
    float3  _S504 = - s_diff_dy_T_0;
    float3  _S505 = - _S501.differential_0;
    float3  _S506 = raydir_31 * s_diff_dy_T_0;
    float _S507 = _S506.x + _S506.y + _S506.z;
    float4  _S508 = make_float4 (0.0f);
    float4  _S509 = _S508;
    *&((&_S509)->w) = _S507;
    float3  _S510 = raydir_28 * _S504;
    float _S511 = _S510.x + _S510.y + _S510.z;
    float4  _S512 = _S508;
    *&((&_S512)->z) = _S511;
    float4  _S513 = _S509 + _S512;
    float3  _S514 = raydir_25 * _S503.differential_0;
    float _S515 = _S514.x + _S514.y + _S514.z;
    float4  _S516 = _S508;
    *&((&_S516)->y) = _S515;
    float4  _S517 = _S513 + _S516;
    float3  _S518 = raydir_22 * _S505;
    float _S519 = _S518.x + _S518.y + _S518.z;
    float4  _S520 = _S508;
    *&((&_S520)->x) = _S519;
    float4  _S521 = _S517 + _S520;
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S521;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(uint _S522, uint _S523, float2  _S524, float4  _S525, float4  _S526, float2  _S527, float2  _S528, bool _S529, bool _S530, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S531, float3  _S532)
{
    float3  _S533;
    s_bwd_prop_depth_to_normal_Intermediates_0 _S534;
    s_primal_ctx_depth_to_normal_0(_S522, _S523, _S524, _S525, _S526, _S527, _S528, _S529, _S530, (*_S531).primal_0, &_S533, &_S534);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S535 = _S534;
    s_bwd_prop_depth_to_normal_0(_S522, _S523, _S524, _S525, _S526, _S527, _S528, _S529, _S530, _S531, _S532, &_S535);
    return;
}

inline __device__ void depth_to_normal_vjp(uint width_3, uint height_3, float2  pix_center_3, float4  intrins_3, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, bool is_fisheye_9, bool is_ray_depth_3, float4  depths_1, float3  v_normal_0, float4  * v_depths_0)
{
    float4  _S536 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S536;
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
        CameraDistortion_0 _S537;
        (&_S537)->radial_coeffs_0 = radial_coeffs_12;
        (&_S537)->tangential_coeffs_0 = tangential_coeffs_12;
        (&_S537)->thin_prism_coeffs_0 = thin_prism_coeffs_12;
        float2  _S538 = undistort_point_0(uv_19, &_S537, int(12));
        float theta_7 = length_1(_S538);
        float3  raydir_35 = make_float3 ((_S538 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (_S538 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        raydir_34 = raydir_35 / make_float3 (raydir_35.z);
    }
    else
    {
        raydir_34 = make_float3 (uv_19.x, uv_19.y, 1.0f);
    }
    return float((F32_sign((raydir_34.z)))) / length_2(raydir_34);
}

