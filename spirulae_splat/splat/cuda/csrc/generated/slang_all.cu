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
    float2  _S336 = uv_6;
    bool _S337 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S336);
    if(!_S337)
    {
        int3  _S338 = make_int3 (int(0));
        float3  _S339 = make_float3 ((float)_S338.x, (float)_S338.y, (float)_S338.z);
        *raydir_2 = _S339;
        return false;
    }
    float3  _S340;
    if(is_fisheye_3)
    {
        float2  _S341 = _S336;
        float theta_3 = length_1(_S336);
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
        _S340 = _S343;
    }
    else
    {
        _S340 = make_float3 (_S336.x, _S336.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S340);
    return true;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, float3  dOut_8)
{
    float _S344 = (*right_2).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  right_d_result_0;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = (*left_2).primal_0.x * dOut_8.x;
    float sum_4 = _S344 + (*right_2).primal_0.rows[int(0)].y * dOut_8.y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = (*left_2).primal_0.x * dOut_8.y;
    float sum_5 = sum_4 + (*right_2).primal_0.rows[int(0)].z * dOut_8.z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = (*left_2).primal_0.x * dOut_8.z;
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = sum_5;
    float _S345 = (*right_2).primal_0.rows[int(1)].x * dOut_8.x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = (*left_2).primal_0.y * dOut_8.x;
    float sum_6 = _S345 + (*right_2).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = (*left_2).primal_0.y * dOut_8.y;
    float sum_7 = sum_6 + (*right_2).primal_0.rows[int(1)].z * dOut_8.z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = (*left_2).primal_0.y * dOut_8.z;
    *&((&left_d_result_0)->y) = sum_7;
    float _S346 = (*right_2).primal_0.rows[int(2)].x * dOut_8.x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = (*left_2).primal_0.z * dOut_8.x;
    float sum_8 = _S346 + (*right_2).primal_0.rows[int(2)].y * dOut_8.y;
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

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_2, float3  t_1)
{
    return - mul_2(t_1, R_2);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_3, float3  raydir_3)
{
    return mul_2(raydir_3, R_3);
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S347, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S348, float3  _S349)
{
    _d_mul_0(_S347, _S348, _S349);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_2)
{
    float3  _S350 = - _s_dOut_2;
    float3  _S351 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S352;
    (&_S352)->primal_0 = (*dpt_0).primal_0;
    (&_S352)->differential_0 = _S351;
    Matrix<float, 3, 3>  _S353 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S354;
    (&_S354)->primal_0 = (*dpR_0).primal_0;
    (&_S354)->differential_0 = _S353;
    s_bwd_prop_mul_0(&_S352, &_S354, _S350);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S352.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S354.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S355, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S356, float3  _S357)
{
    s_bwd_prop_transform_ray_o_0(_S355, _S356, _S357);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_4, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S358 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_4;
    (&dp_R_0)->differential_0 = _S358;
    float3  _S359 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S359;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_3)
{
    float3  _S360 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S361;
    (&_S361)->primal_0 = (*dpraydir_0).primal_0;
    (&_S361)->differential_0 = _S360;
    Matrix<float, 3, 3>  _S362 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S363;
    (&_S363)->primal_0 = (*dpR_1).primal_0;
    (&_S363)->differential_0 = _S362;
    s_bwd_prop_mul_0(&_S361, &_S363, _s_dOut_3);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S361.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S363.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S364, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S365, float3  _S366)
{
    s_bwd_prop_transform_ray_d_0(_S364, _S365, _S366);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_5, float3  raydir_4, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S367 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_5;
    (&dp_R_1)->differential_0 = _S367;
    float3  _S368 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_4;
    (&dp_raydir_0)->differential_0 = _S368;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_5, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S369 = scale_2.x;
    float sx_0 = (F32_exp((_S369)));
    float _S370 = scale_2.y;
    float sy_0 = (F32_exp((_S370)));
    float sz_0 = scale_2.z - 0.5f * (_S369 + _S370);
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
    Matrix<float, 3, 3>  _S371 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S371, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S371, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S371, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
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
    float4  _S372 = normalize_1(quat_6);
    float3  _S373 = exp_0(scale_3);
    float x_20 = _S372.y;
    float x2_4 = x_20 * x_20;
    float y2_4 = _S372.z * _S372.z;
    float z2_4 = _S372.w * _S372.w;
    float xy_4 = _S372.y * _S372.z;
    float xz_4 = _S372.y * _S372.w;
    float yz_4 = _S372.z * _S372.w;
    float wx_4 = _S372.x * _S372.y;
    float wy_4 = _S372.x * _S372.z;
    float wz_4 = _S372.x * _S372.w;
    Matrix<float, 3, 3>  M_2 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S373.x, 0.0f, 0.0f, 0.0f, _S373.y, 0.0f, 0.0f, 0.0f, _S373.z));
    float4  _S374 = make_float4 (dot_0(*mean_1, *mean_1), dot_0(*mean_1, scale_3), dot_0(scale_3, scale_3), dot_1(quat_6, make_float4 (opac_0))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S375 = _S374 - floor_0(_S374);
    float4  _S376 = _S375 + make_float4 (dot_1(_S375, float4 {_S375.w, _S375.z, _S375.x, _S375.y} + make_float4 (33.3300018310546875f)));
    float4  _S377 = (float4 {_S376.x, _S376.x, _S376.y, _S376.z} + float4 {_S376.y, _S376.z, _S376.z, _S376.w}) * float4 {_S376.z, _S376.y, _S376.w, _S376.x};
    float4  _S378 = _S377 - floor_0(_S377);
    float2  _S379 = float2 {_S378.x, _S378.z};
    float _S380 = 6.28318548202514648f * _S379.y;
    float2  _S381 = float2 {_S378.y, _S378.w};
    float _S382 = 6.28318548202514648f * _S381.y;
    *mean_1 = *mean_1 + mul_0(mul_1(M_2, transpose_0(M_2)), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S379.x))))))) * make_float2 ((F32_cos((_S380))), (F32_sin((_S380))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S379.x))))))) * make_float2 ((F32_cos((_S380))), (F32_sin((_S380))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S381.x))))))) * make_float2 ((F32_cos((_S382))), (F32_sin((_S382))))).x) * make_float3 (scaler_0) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_0) * (1.0f - opac_0 - (1.0f - min_opacity_0))))))));
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_9)
{
    float _S383 = dOut_9.y;
    float _S384 = dOut_9.z;
    float _S385 = dOut_9.x;
    float _S386 = (*a_0).primal_0.z * _S383 + - (*a_0).primal_0.y * _S384;
    float _S387 = - (*a_0).primal_0.z * _S385 + (*a_0).primal_0.x * _S384;
    float _S388 = (*a_0).primal_0.y * _S385 + - (*a_0).primal_0.x * _S383;
    float3  _S389 = make_float3 (- (*b_0).primal_0.z * _S383 + (*b_0).primal_0.y * _S384, (*b_0).primal_0.z * _S385 + - (*b_0).primal_0.x * _S384, - (*b_0).primal_0.y * _S385 + (*b_0).primal_0.x * _S383);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S389;
    float3  _S390 = make_float3 (_S386, _S387, _S388);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S390;
    return;
}

inline __device__ float3  cross_0(float3  left_4, float3  right_4)
{
    float _S391 = left_4.y;
    float _S392 = right_4.z;
    float _S393 = left_4.z;
    float _S394 = right_4.y;
    float _S395 = right_4.x;
    float _S396 = left_4.x;
    return make_float3 (_S391 * _S392 - _S393 * _S394, _S393 * _S395 - _S396 * _S392, _S396 * _S394 - _S391 * _S395);
}

inline __device__ void mcmc_add_noise_triangle(float scaler_1, float min_opacity_1, float3  * mean_2, float3  scale_4, float4  quat_7, float opac_1)
{
    float4  _S397 = normalize_1(quat_7);
    float _S398 = scale_4.x;
    float sx_1 = (F32_exp((_S398)));
    float _S399 = scale_4.y;
    float sy_1 = (F32_exp((_S399)));
    float sz_1 = scale_4.z - 0.5f * (_S398 + _S399);
    float x_21 = _S397.y;
    float x2_5 = x_21 * x_21;
    float y2_5 = _S397.z * _S397.z;
    float z2_5 = _S397.w * _S397.w;
    float xy_5 = _S397.y * _S397.z;
    float xz_5 = _S397.y * _S397.w;
    float yz_5 = _S397.z * _S397.w;
    float wx_5 = _S397.x * _S397.y;
    float wy_5 = _S397.x * _S397.z;
    float wz_5 = _S397.x * _S397.w;
    Matrix<float, 3, 3>  _S400 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5)));
    float3  vert0_1 = mul_0(_S400, make_float3 (sx_1, 0.0f, 0.0f)) + *mean_2;
    float3  vert1_1 = mul_0(_S400, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + *mean_2;
    float3  vert2_1 = mul_0(_S400, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + *mean_2;
    float3  vertc_0 = (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f);
    float3  d0_0 = vert0_1 - vertc_0;
    float3  d1_0 = vert1_1 - vertc_0;
    float3  d2_0 = vert2_1 - vertc_0;
    float3  dn_0 = make_float3 (0.5f * (F32_min(((F32_min((length_2(d0_0)), (length_2(d1_0))))), (length_2(d2_0))))) * normalize_0(cross_0(d0_0, d1_0));
    float4  _S401 = make_float4 (dot_0(*mean_2, *mean_2), dot_0(*mean_2, scale_4), dot_0(scale_4, scale_4), dot_1(quat_7, make_float4 (opac_1))) * make_float4 (0.1031000018119812f, 0.10300000011920929f, 0.09730000048875809f, 0.10989999771118164f);
    float4  _S402 = _S401 - floor_0(_S401);
    float4  _S403 = _S402 + make_float4 (dot_1(_S402, float4 {_S402.w, _S402.z, _S402.x, _S402.y} + make_float4 (33.3300018310546875f)));
    float4  _S404 = (float4 {_S403.x, _S403.x, _S403.y, _S403.z} + float4 {_S403.y, _S403.z, _S403.z, _S403.w}) * float4 {_S403.z, _S403.y, _S403.w, _S403.x};
    float4  _S405 = _S404 - floor_0(_S404);
    float2  _S406 = float2 {_S405.x, _S405.z};
    float _S407 = 6.28318548202514648f * _S406.y;
    float2  _S408 = float2 {_S405.y, _S405.w};
    float _S409 = 6.28318548202514648f * _S408.y;
    *mean_2 = *mean_2 + mul_0(makeMatrix<float, 3, 3> (0.5f) * (makeMatrix<float, 3, 3> (make_float3 (d0_0.x) * d0_0, make_float3 (d0_0.y) * d0_0, make_float3 (d0_0.z) * d0_0) + makeMatrix<float, 3, 3> (make_float3 (d1_0.x) * d1_0, make_float3 (d1_0.y) * d1_0, make_float3 (d1_0.z) * d1_0) + makeMatrix<float, 3, 3> (make_float3 (d2_0.x) * d2_0, make_float3 (d2_0.y) * d2_0, make_float3 (d2_0.z) * d2_0) + makeMatrix<float, 3, 3> (make_float3 (dn_0.x) * dn_0, make_float3 (dn_0.y) * dn_0, make_float3 (dn_0.z) * dn_0)) / makeMatrix<float, 3, 3> (3.5f), make_float3 ((make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S406.x))))))) * make_float2 ((F32_cos((_S407))), (F32_sin((_S407))))).x, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S406.x))))))) * make_float2 ((F32_cos((_S407))), (F32_sin((_S407))))).y, (make_float2 ((F32_sqrt((-2.0f * (F32_log((1.0f - _S408.x))))))) * make_float2 ((F32_cos((_S409))), (F32_sin((_S409))))).x) * make_float3 (scaler_1) * make_float3 (1.0f / (1.0f + (F32_exp((- (0.5f / min_opacity_1) * (1.0f - opac_1 - (1.0f - min_opacity_1))))))));
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_3, float3  dOut_10)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_9).primal_0.x;
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
    (&left_dp_1)->primal_0 = (*dpx_9).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_3).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_10.y);
    *&((&left_d_result_1)->y) = left_dp_1.differential_0;
    *&((&right_d_result_1)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_9).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_3).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_10.z);
    *&((&left_d_result_1)->z) = left_dp_2.differential_0;
    *&((&right_d_result_1)->z) = right_dp_2.differential_0;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = left_d_result_1;
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

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_11)
{
    float3  _S410 = make_float3 (1.0f) / (*dpx_10).primal_0 * dOut_11;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S410;
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

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, float3  dOut_12)
{
    float3  _S411 = _slang_select(((*dpx_11).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_11).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_12;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S411;
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

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_13)
{
    DiffPair_float_0 _S412 = *dpx_12;
    bool _S413;
    if(((*dpx_12).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S413 = ((*dpx_12).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S413 = false;
    }
    float _S414;
    if(_S413)
    {
        _S414 = dOut_13;
    }
    else
    {
        _S414 = 0.0f;
    }
    dpx_12->primal_0 = _S412.primal_0;
    dpx_12->differential_0 = _S414;
    DiffPair_float_0 _S415 = *dpMin_0;
    if((_S412.primal_0) < ((*dpMin_0).primal_0))
    {
        _S414 = dOut_13;
    }
    else
    {
        _S414 = 0.0f;
    }
    dpMin_0->primal_0 = _S415.primal_0;
    dpMin_0->differential_0 = _S414;
    DiffPair_float_0 _S416 = *dpMax_0;
    if(((*dpx_12).primal_0) > ((*dpMax_0).primal_0))
    {
        _S414 = dOut_13;
    }
    else
    {
        _S414 = 0.0f;
    }
    dpMax_0->primal_0 = _S416.primal_0;
    dpMax_0->differential_0 = _S414;
    return;
}

inline __device__ float clamp_0(float x_25, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_25), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_13, float dOut_14)
{
    float _S417 = -0.5f / ((*dpx_13).primal_0 * (F32_sqrt(((*dpx_13).primal_0)))) * dOut_14;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S417;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_14, DiffPair_float_0 * dpy_4, DiffPair_float_0 * dps_0, float dOut_15)
{
    float _S418 = (1.0f - (*dps_0).primal_0) * dOut_15;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S418;
    DiffPair_float_0 _S419 = *dpy_4;
    float _S420 = (*dps_0).primal_0 * dOut_15;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S420;
    float _S421 = (_S419.primal_0 - (*dpx_14).primal_0) * dOut_15;
    dps_0->primal_0 = _S419.primal_0;
    dps_0->differential_0 = _S421;
    return;
}

inline __device__ float lerp_0(float x_26, float y_8, float s_4)
{
    return x_26 + (y_8 - x_26) * s_4;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_alpha_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 11>  * weights_0, FixedArray<float, 23>  * _S422)
{
    float3  _S423;
    bool _S424;
    bool _S425;
    FixedArray<float, 23>  losses_1;
    float _S426 = float(mask_0);
    float _S427 = 1.0f - (*weights_0)[int(0)];
    float3  _S428 = make_float3 (0.0f);
    float3  _S429 = abs_0(ref_rgb_0 * make_float3 (_S427) + (log_0(max_0(ref_rgb_0, _S428) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)]) - (render_rgb_0 * make_float3 (_S427) + (log_0(max_0(render_rgb_0, _S428) + make_float3 (0.00392156885936856f)) + make_float3 (1.0f)) * make_float3 ((*weights_0)[int(0)])));
    losses_1[int(0)] = (*weights_0)[int(1)] * _S426 * ((_S429.x + _S429.y + _S429.z) * 0.3333333432674408f);
    float3  _S430 = ref_rgb_0 - render_rgb_0;
    losses_1[int(1)] = _S426 * clamp_0(dot_0(_S430, _S430) * 0.3333333432674408f, 0.0f, 1.0f);
    float _S431 = float(depth_mask_0 & mask_0);
    float _S432 = _S431 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S433 = _S431 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_1[int(2)] = _S432;
    losses_1[int(3)] = _S433;
    losses_1[int(4)] = _S432 * _S432;
    losses_1[int(5)] = _S433 * _S433;
    losses_1[int(6)] = _S432 * _S433;
    bool _S434 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S435 = norm2_0 == 0.0f;
        _S424 = _S435;
        if(_S435)
        {
            _S423 = _S428;
            break;
        }
        _S423 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S436;
    bool _S437 = !_S424;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S438 = norm2_1 == 0.0f;
        _S425 = _S438;
        if(_S438)
        {
            _S436 = _S428;
            break;
        }
        _S436 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S439;
    float3  _S440;
    bool _S441 = !_S425;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S440 = _S428;
            _S439 = false;
            break;
        }
        _S440 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S439 = _S434;
        break;
    }
    float _S442 = float(_S437 & _S439);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S423, _S440);
    losses_1[int(7)] = (*weights_0)[int(3)] * _S442 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S443 = float(_S441 & _S439);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S436, _S440);
    losses_1[int(8)] = (*weights_0)[int(3)] * _S443 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S444 = float(_S437 & _S441);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S423, _S436);
    losses_1[int(11)] = (*weights_0)[int(6)] * _S444 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float _S445 = clamp_0(render_alpha_0, 0.0f, 1.0f);
    float _S446 = float(alpha_mask_0);
    float _S447 = float(ref_alpha_0);
    float _S448 = (F32_max((_S445), (_S447)));
    losses_1[int(9)] = (*weights_0)[int(4)] * _S446 * - lerp_0((F32_log(((F32_max((1.0f - _S448), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S448), (9.99999997475242708e-07f)))))), _S447);
    float _S449 = 1.0f - _S445;
    float _S450 = 1.0f - _S447;
    float _S451 = (F32_max((_S449), (_S450)));
    losses_1[int(10)] = (*weights_0)[int(5)] * _S446 * - lerp_0((F32_log(((F32_max((1.0f - _S451), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S451), (9.99999997475242708e-07f)))))), _S450);
    losses_1[int(12)] = (*weights_0)[int(7)] * 4.0f * _S445 * _S449;
    float _S452 = (F32_max((_S445), (9.999999960041972e-13f)));
    losses_1[int(13)] = (*weights_0)[int(8)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S452;
    losses_1[int(14)] = (*weights_0)[int(9)] * depth_dist_0 / _S452;
    losses_1[int(15)] = (*weights_0)[int(10)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S452;
    losses_1[int(16)] = 1.0f;
    losses_1[int(17)] = _S426;
    losses_1[int(18)] = _S431;
    losses_1[int(19)] = _S442;
    losses_1[int(20)] = _S443;
    losses_1[int(21)] = _S444;
    losses_1[int(22)] = _S446;
    *_S422 = losses_1;
    return;
}

inline __device__ float3  s_primal_ctx_max_1(float3  _S453, float3  _S454)
{
    return max_0(_S453, _S454);
}

inline __device__ float3  s_primal_ctx_log_1(float3  _S455)
{
    return log_0(_S455);
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

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S477, float3  _S478)
{
    _d_log_vector_0(_S477, _S478);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S479, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S480, float3  _S481)
{
    _d_max_vector_0(_S479, _S480, _S481);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_alpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 11>  * weights_1, FixedArray<float, 23>  * _s_dOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S482 = *dprender_rgb_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S483 = *dpref_rgb_0;
    DiffPair_float_0 _S484 = *dprender_depth_0;
    DiffPair_float_0 _S485 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S486 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S487 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S488 = *dpref_normal_0;
    DiffPair_float_0 _S489 = *dprender_alpha_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S490 = *dprgb_dist_0;
    DiffPair_float_0 _S491 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S492 = *dpnormal_dist_0;
    float3  _S493 = make_float3 (0.0f);
    float _S494 = float(mask_1);
    float _S495 = (*weights_1)[int(1)] * _S494;
    float3  _S496 = make_float3 ((*weights_1)[int(0)]);
    float _S497 = 1.0f - (*weights_1)[int(0)];
    float3  _S498 = make_float3 (_S497);
    float3  _S499 = make_float3 (0.0f);
    float3  _S500 = s_primal_ctx_max_1((*dprender_rgb_0).primal_0, _S499) + make_float3 (0.00392156885936856f);
    float3  _S501 = s_primal_ctx_max_1((*dpref_rgb_0).primal_0, _S499) + make_float3 (0.00392156885936856f);
    float3  _S502 = (*dpref_rgb_0).primal_0 * make_float3 (_S497) + (s_primal_ctx_log_1(_S501) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]) - ((*dprender_rgb_0).primal_0 * make_float3 (_S497) + (s_primal_ctx_log_1(_S500) + make_float3 (1.0f)) * make_float3 ((*weights_1)[int(0)]));
    float3  _S503 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S504 = s_primal_ctx_dot_0(_S503, _S503) * 0.3333333432674408f;
    float _S505 = float(depth_mask_1 & mask_1);
    float _S506 = s_primal_ctx_max_0((*dprender_depth_0).primal_0, 0.00009999999747379f);
    float _S507 = _S505 * s_primal_ctx_log_0(_S506);
    float _S508 = s_primal_ctx_max_0((*dpref_depth_0).primal_0, 0.00009999999747379f);
    float _S509 = _S505 * s_primal_ctx_log_0(_S508);
    bool _S510 = normal_mask_1 & mask_1;
    float _S511 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S512 = !(_S511 == 0.0f);
    float3  _S513;
    float3  _S514;
    if(_S512)
    {
        float _S515 = s_primal_ctx_rsqrt_0(_S511);
        float3  _S516 = make_float3 (_S515);
        _S513 = _S486.primal_0 * make_float3 (_S515);
        _S514 = _S516;
    }
    else
    {
        _S513 = _S499;
        _S514 = _S493;
    }
    float _S517 = s_primal_ctx_dot_0(_S487.primal_0, _S487.primal_0);
    bool _S518 = !(_S517 == 0.0f);
    float3  _S519;
    float3  _S520;
    if(_S518)
    {
        float _S521 = s_primal_ctx_rsqrt_0(_S517);
        float3  _S522 = make_float3 (_S521);
        _S519 = _S487.primal_0 * make_float3 (_S521);
        _S520 = _S522;
    }
    else
    {
        _S519 = _S499;
        _S520 = _S493;
    }
    float _S523 = s_primal_ctx_dot_0(_S488.primal_0, _S488.primal_0);
    bool _S524 = _S523 == 0.0f;
    bool _S525;
    if(_S524)
    {
        _S525 = false;
    }
    else
    {
        _S525 = _S510;
    }
    bool _S526 = !_S524;
    float3  _S527;
    float3  _S528;
    if(_S526)
    {
        float _S529 = s_primal_ctx_rsqrt_0(_S523);
        float3  _S530 = make_float3 (_S529);
        _S527 = _S488.primal_0 * make_float3 (_S529);
        _S528 = _S530;
    }
    else
    {
        _S527 = _S499;
        _S528 = _S493;
    }
    float _S531 = (*weights_1)[int(3)] * float(_S512 & _S525);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S513, _S527);
    float _S532 = s_primal_ctx_max_0(cos_sim_loss_3, 9.999999960041972e-13f);
    float _S533 = (*weights_1)[int(3)] * float(_S518 & _S525);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S519, _S527);
    float _S534 = s_primal_ctx_max_0(cos_sim_loss_4, 9.999999960041972e-13f);
    float _S535 = (*weights_1)[int(6)] * float(_S512 & _S518);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S513, _S519);
    float _S536 = s_primal_ctx_max_0(cos_sim_loss_5, 9.999999960041972e-13f);
    float _S537 = s_primal_ctx_clamp_0(_S489.primal_0, 0.0f, 1.0f);
    float _S538 = float(alpha_mask_1);
    float _S539 = (*weights_1)[int(4)] * _S538;
    float _S540 = float(ref_alpha_1);
    float _S541 = s_primal_ctx_max_0(_S537, _S540);
    float _S542 = 1.0f - _S541;
    float _S543 = s_primal_ctx_max_0(_S542, 9.99999997475242708e-07f);
    float _S544 = s_primal_ctx_log_0(_S543);
    float _S545 = s_primal_ctx_max_0(_S541, 9.99999997475242708e-07f);
    float _S546 = s_primal_ctx_log_0(_S545);
    float _S547 = (*weights_1)[int(5)] * _S538;
    float _S548 = 1.0f - _S537;
    float _S549 = 1.0f - _S540;
    float _S550 = s_primal_ctx_max_0(_S548, _S549);
    float _S551 = 1.0f - _S550;
    float _S552 = s_primal_ctx_max_0(_S551, 9.99999997475242708e-07f);
    float _S553 = s_primal_ctx_log_0(_S552);
    float _S554 = s_primal_ctx_max_0(_S550, 9.99999997475242708e-07f);
    float _S555 = s_primal_ctx_log_0(_S554);
    float _S556 = (*weights_1)[int(7)] * 4.0f;
    float _S557 = _S556 * _S537;
    float _S558 = s_primal_ctx_max_0(_S537, 9.999999960041972e-13f);
    float _S559 = _S558 * _S558;
    float _S560 = (*_s_dOut_4)[int(15)] / _S559;
    float _S561 = 0.3333333432674408f * ((*weights_1)[int(10)] * (_S558 * _S560));
    float _S562 = (*_s_dOut_4)[int(14)] / _S559;
    float _S563 = (*weights_1)[int(9)] * (_S558 * _S562);
    float _S564 = (*_s_dOut_4)[int(13)] / _S559;
    float _S565 = _S558 * _S564;
    float _S566 = (*weights_1)[int(10)] * ((_S492.primal_0.x + _S492.primal_0.y + _S492.primal_0.z) * 0.3333333432674408f) * - _S560 + (*weights_1)[int(9)] * _S491.primal_0 * - _S562 + (*weights_1)[int(8)] * ((_S490.primal_0.x + _S490.primal_0.y + _S490.primal_0.z) * 0.3333333432674408f) * - _S564;
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S537;
    (&_S567)->differential_0 = 0.0f;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = 9.999999960041972e-13f;
    (&_S568)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S567, &_S568, _S566);
    float _S569 = 0.3333333432674408f * ((*weights_1)[int(8)] * _S565);
    float _S570 = _S557 * (*_s_dOut_4)[int(12)];
    float _S571 = _S556 * (_S548 * (*_s_dOut_4)[int(12)]);
    float _S572 = - (_S547 * (*_s_dOut_4)[int(10)]);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S553;
    (&_S573)->differential_0 = 0.0f;
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = _S555;
    (&_S574)->differential_0 = 0.0f;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = _S549;
    (&_S575)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S573, &_S574, &_S575, _S572);
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S554;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S576, _S574.differential_0);
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = _S550;
    (&_S577)->differential_0 = 0.0f;
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = 9.99999997475242708e-07f;
    (&_S578)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S577, &_S578, _S576.differential_0);
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S552;
    (&_S579)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S579, _S573.differential_0);
    DiffPair_float_0 _S580;
    (&_S580)->primal_0 = _S551;
    (&_S580)->differential_0 = 0.0f;
    DiffPair_float_0 _S581;
    (&_S581)->primal_0 = 9.99999997475242708e-07f;
    (&_S581)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S580, &_S581, _S579.differential_0);
    float _S582 = _S577.differential_0 + - _S580.differential_0;
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S548;
    (&_S583)->differential_0 = 0.0f;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S549;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S583, &_S584, _S582);
    float _S585 = - (_S570 + _S583.differential_0);
    float _S586 = - (_S539 * (*_s_dOut_4)[int(9)]);
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = _S544;
    (&_S587)->differential_0 = 0.0f;
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = _S546;
    (&_S588)->differential_0 = 0.0f;
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = _S540;
    (&_S589)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S587, &_S588, &_S589, _S586);
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = _S545;
    (&_S590)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S590, _S588.differential_0);
    DiffPair_float_0 _S591;
    (&_S591)->primal_0 = _S541;
    (&_S591)->differential_0 = 0.0f;
    DiffPair_float_0 _S592;
    (&_S592)->primal_0 = 9.99999997475242708e-07f;
    (&_S592)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S591, &_S592, _S590.differential_0);
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = _S543;
    (&_S593)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S593, _S587.differential_0);
    DiffPair_float_0 _S594;
    (&_S594)->primal_0 = _S542;
    (&_S594)->differential_0 = 0.0f;
    DiffPair_float_0 _S595;
    (&_S595)->primal_0 = 9.99999997475242708e-07f;
    (&_S595)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S594, &_S595, _S593.differential_0);
    float _S596 = _S591.differential_0 + - _S594.differential_0;
    DiffPair_float_0 _S597;
    (&_S597)->primal_0 = _S537;
    (&_S597)->differential_0 = 0.0f;
    DiffPair_float_0 _S598;
    (&_S598)->primal_0 = _S540;
    (&_S598)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S597, &_S598, _S596);
    float _S599 = _S567.differential_0 + _S571 + _S585 + _S597.differential_0;
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = _S489.primal_0;
    (&_S600)->differential_0 = 0.0f;
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = 0.0f;
    (&_S601)->differential_0 = 0.0f;
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = 1.0f;
    (&_S602)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S600, &_S601, &_S602, _S599);
    DiffPair_float_0 _S603 = _S600;
    float _S604 = _S535 * (*_s_dOut_4)[int(11)];
    DiffPair_float_0 _S605;
    (&_S605)->primal_0 = _S536;
    (&_S605)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S605, _S604);
    DiffPair_float_0 _S606;
    (&_S606)->primal_0 = cos_sim_loss_5;
    (&_S606)->differential_0 = 0.0f;
    DiffPair_float_0 _S607;
    (&_S607)->primal_0 = 9.999999960041972e-13f;
    (&_S607)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S606, &_S607, _S605.differential_0);
    float _S608 = 0.5f * - (_S604 + _S606.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S609;
    (&_S609)->primal_0 = _S513;
    (&_S609)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S610;
    (&_S610)->primal_0 = _S519;
    (&_S610)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S609, &_S610, _S608);
    float _S611 = _S533 * (*_s_dOut_4)[int(8)];
    DiffPair_float_0 _S612;
    (&_S612)->primal_0 = _S534;
    (&_S612)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S612, _S611);
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = cos_sim_loss_4;
    (&_S613)->differential_0 = 0.0f;
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = 9.999999960041972e-13f;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S613, &_S614, _S612.differential_0);
    float _S615 = 0.5f * - (_S611 + _S613.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S616;
    (&_S616)->primal_0 = _S519;
    (&_S616)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S617;
    (&_S617)->primal_0 = _S527;
    (&_S617)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S616, &_S617, _S615);
    float _S618 = _S531 * (*_s_dOut_4)[int(7)];
    DiffPair_float_0 _S619;
    (&_S619)->primal_0 = _S532;
    (&_S619)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S619, _S618);
    DiffPair_float_0 _S620;
    (&_S620)->primal_0 = cos_sim_loss_3;
    (&_S620)->differential_0 = 0.0f;
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = 9.999999960041972e-13f;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S620, &_S621, _S619.differential_0);
    float _S622 = 0.5f * - (_S618 + _S620.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S623;
    (&_S623)->primal_0 = _S513;
    (&_S623)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S624;
    (&_S624)->primal_0 = _S527;
    (&_S624)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S623, &_S624, _S622);
    float3  _S625 = _S617.differential_0 + _S624.differential_0;
    float3  _S626 = _S609.differential_0 + _S623.differential_0;
    float3  _S627 = make_float3 (_S561, _S561, _S561);
    float3  _S628 = make_float3 (_S569, _S569, _S569);
    float3  _S629 = _S610.differential_0 + _S616.differential_0;
    float _S630;
    if(_S526)
    {
        float3  _S631 = _S488.primal_0 * _S625;
        float3  _S632 = _S528 * _S625;
        float _S633 = _S631.x + _S631.y + _S631.z;
        DiffPair_float_0 _S634;
        (&_S634)->primal_0 = _S523;
        (&_S634)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S634, _S633);
        _S630 = _S634.differential_0;
        _S513 = _S632;
    }
    else
    {
        _S630 = 0.0f;
        _S513 = _S493;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S635;
    (&_S635)->primal_0 = _S488.primal_0;
    (&_S635)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S636;
    (&_S636)->primal_0 = _S488.primal_0;
    (&_S636)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S635, &_S636, _S630);
    float3  _S637 = _S636.differential_0 + _S635.differential_0 + _S513;
    if(_S518)
    {
        float3  _S638 = _S487.primal_0 * _S629;
        float3  _S639 = _S520 * _S629;
        float _S640 = _S638.x + _S638.y + _S638.z;
        DiffPair_float_0 _S641;
        (&_S641)->primal_0 = _S517;
        (&_S641)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S641, _S640);
        _S630 = _S641.differential_0;
        _S513 = _S639;
    }
    else
    {
        _S630 = 0.0f;
        _S513 = _S493;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S642;
    (&_S642)->primal_0 = _S487.primal_0;
    (&_S642)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S643;
    (&_S643)->primal_0 = _S487.primal_0;
    (&_S643)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S642, &_S643, _S630);
    float3  _S644 = _S643.differential_0 + _S642.differential_0 + _S513;
    if(_S512)
    {
        float3  _S645 = _S486.primal_0 * _S626;
        float3  _S646 = _S514 * _S626;
        float _S647 = _S645.x + _S645.y + _S645.z;
        DiffPair_float_0 _S648;
        (&_S648)->primal_0 = _S511;
        (&_S648)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S648, _S647);
        _S630 = _S648.differential_0;
        _S513 = _S646;
    }
    else
    {
        _S630 = 0.0f;
        _S513 = _S493;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S649;
    (&_S649)->primal_0 = _S486.primal_0;
    (&_S649)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S650;
    (&_S650)->primal_0 = _S486.primal_0;
    (&_S650)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S649, &_S650, _S630);
    float _S651 = _S509 * (*_s_dOut_4)[int(6)];
    float _S652 = _S509 * (*_s_dOut_4)[int(5)];
    float _S653 = _S507 * (*_s_dOut_4)[int(4)];
    float _S654 = _S505 * (_S507 * (*_s_dOut_4)[int(6)] + _S652 + _S652 + (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = _S508;
    (&_S655)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S655, _S654);
    DiffPair_float_0 _S656;
    (&_S656)->primal_0 = _S485.primal_0;
    (&_S656)->differential_0 = 0.0f;
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = 0.00009999999747379f;
    (&_S657)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S656, &_S657, _S655.differential_0);
    float _S658 = _S505 * (_S651 + _S653 + _S653 + (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S659;
    (&_S659)->primal_0 = _S506;
    (&_S659)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S659, _S658);
    DiffPair_float_0 _S660;
    (&_S660)->primal_0 = _S484.primal_0;
    (&_S660)->differential_0 = 0.0f;
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = 0.00009999999747379f;
    (&_S661)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S660, &_S661, _S659.differential_0);
    float _S662 = _S494 * (*_s_dOut_4)[int(1)];
    DiffPair_float_0 _S663;
    (&_S663)->primal_0 = _S504;
    (&_S663)->differential_0 = 0.0f;
    DiffPair_float_0 _S664;
    (&_S664)->primal_0 = 0.0f;
    (&_S664)->differential_0 = 0.0f;
    DiffPair_float_0 _S665;
    (&_S665)->primal_0 = 1.0f;
    (&_S665)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S663, &_S664, &_S665, _S662);
    float _S666 = 0.3333333432674408f * _S663.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S667;
    (&_S667)->primal_0 = _S503;
    (&_S667)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S668;
    (&_S668)->primal_0 = _S503;
    (&_S668)->differential_0 = _S493;
    s_bwd_prop_dot_0(&_S667, &_S668, _S666);
    float3  _S669 = _S668.differential_0 + _S667.differential_0;
    float3  _S670 = - _S669;
    float _S671 = 0.3333333432674408f * (_S495 * (*_s_dOut_4)[int(0)]);
    float3  _S672 = make_float3 (_S671, _S671, _S671);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S673;
    (&_S673)->primal_0 = _S502;
    (&_S673)->differential_0 = _S493;
    s_bwd_prop_abs_0(&_S673, _S672);
    float3  _S674 = - _S673.differential_0;
    float3  _S675 = _S496 * _S673.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S676;
    (&_S676)->primal_0 = _S501;
    (&_S676)->differential_0 = _S493;
    s_bwd_prop_log_1(&_S676, _S675);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S677;
    (&_S677)->primal_0 = _S483.primal_0;
    (&_S677)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S678;
    (&_S678)->primal_0 = _S499;
    (&_S678)->differential_0 = _S493;
    s_bwd_prop_max_1(&_S677, &_S678, _S676.differential_0);
    float3  _S679 = _S498 * _S673.differential_0;
    float3  _S680 = _S496 * _S674;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S681;
    (&_S681)->primal_0 = _S500;
    (&_S681)->differential_0 = _S493;
    s_bwd_prop_log_1(&_S681, _S680);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S682;
    (&_S682)->primal_0 = _S482.primal_0;
    (&_S682)->differential_0 = _S493;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S683;
    (&_S683)->primal_0 = _S499;
    (&_S683)->differential_0 = _S493;
    s_bwd_prop_max_1(&_S682, &_S683, _S681.differential_0);
    float3  _S684 = _S498 * _S674;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S627;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S563;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S628;
    dprender_alpha_0->primal_0 = (*dprender_alpha_0).primal_0;
    dprender_alpha_0->differential_0 = _S603.differential_0;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S637;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S644;
    float3  _S685 = _S650.differential_0 + _S649.differential_0 + _S513;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S685;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S656.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S660.differential_0;
    float3  _S686 = _S669 + _S677.differential_0 + _S679;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S686;
    float3  _S687 = _S670 + _S682.differential_0 + _S684;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S687;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S688, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S689, DiffPair_float_0 * _S690, DiffPair_float_0 * _S691, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S692, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S693, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S694, DiffPair_float_0 * _S695, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S696, DiffPair_float_0 * _S697, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S698, bool _S699, bool _S700, bool _S701, bool _S702, bool _S703, FixedArray<float, 11>  * _S704, FixedArray<float, 23>  * _S705)
{
    s_bwd_prop_per_pixel_losses_0(_S688, _S689, _S690, _S691, _S692, _S693, _S694, _S695, _S696, _S697, _S698, _S699, _S700, _S701, _S702, _S703, _S704, _S705);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_alpha_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 11>  * weights_2, FixedArray<float, 23>  * v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_alpha_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S706 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S706;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S706;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S706;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S706;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S706;
    DiffPair_float_0 dp_render_alpha_0;
    (&dp_render_alpha_0)->primal_0 = render_alpha_1;
    (&dp_render_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S706;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S706;
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_15, float dOut_16)
{
    float _S707 = 1.0f / ((*dpx_15).primal_0 * 52.30258560180664062f) * dOut_16;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S707;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  * raw_losses_0, FixedArray<float, 11>  * weights_3, FixedArray<float, 10>  * _S708)
{
    FixedArray<float, 10>  losses_2;
    float _S709 = (F32_max(((*raw_losses_0)[int(17)]), (1.0f)));
    losses_2[int(0)] = (*raw_losses_0)[int(0)] / _S709;
    losses_2[int(1)] = -10.0f * (F32_log10(((*raw_losses_0)[int(1)] / _S709)));
    bool _S710;
    if(((*raw_losses_0)[int(18)]) > 0.0f)
    {
        _S710 = ((*raw_losses_0)[int(3)]) != 0.0f;
    }
    else
    {
        _S710 = false;
    }
    float _S711;
    if(_S710)
    {
        _S711 = (*weights_3)[int(2)] * clamp_0(1.0f - ((*raw_losses_0)[int(6)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), (((*raw_losses_0)[int(4)] - (*raw_losses_0)[int(2)] * (*raw_losses_0)[int(2)] / (*raw_losses_0)[int(18)]) * ((*raw_losses_0)[int(5)] - (*raw_losses_0)[int(3)] * (*raw_losses_0)[int(3)] / (*raw_losses_0)[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S711 = 0.0f;
    }
    losses_2[int(2)] = _S711;
    losses_2[int(3)] = ((*raw_losses_0)[int(7)] / (F32_max(((*raw_losses_0)[int(19)]), (1.0f))) + (*raw_losses_0)[int(8)] / (F32_max(((*raw_losses_0)[int(20)]), (1.0f)))) / float((I32_max((int(((*raw_losses_0)[int(19)]) > 0.5f) + int(((*raw_losses_0)[int(20)]) > 0.5f)), (int(1)))));
    losses_2[int(4)] = ((*raw_losses_0)[int(9)] + (*raw_losses_0)[int(10)]) / (F32_max(((*raw_losses_0)[int(22)]), (1.0f)));
    losses_2[int(5)] = (*raw_losses_0)[int(11)] / (F32_max(((*raw_losses_0)[int(21)]), (1.0f)));
    float _S712 = (F32_max(((*raw_losses_0)[int(16)]), (1.0f)));
    losses_2[int(6)] = (*raw_losses_0)[int(12)] / _S712;
    losses_2[int(7)] = (*raw_losses_0)[int(13)] / _S712;
    losses_2[int(8)] = (*raw_losses_0)[int(14)] / _S712;
    losses_2[int(9)] = (*raw_losses_0)[int(15)] / _S712;
    *_S708 = losses_2;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S713)
{
    return (F32_sqrt((_S713)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S714, float _S715)
{
    _d_log10_0(_S714, _S715);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 11>  * weights_4, FixedArray<float, 10>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S716 = dpraw_losses_0->primal_0;
    float _S717 = s_primal_ctx_max_0(dpraw_losses_0->primal_0[int(17)], 1.0f);
    float _S718 = _S717 * _S717;
    float _S719 = dpraw_losses_0->primal_0[int(1)] / _S717;
    bool _S720 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S721;
    if(_S720)
    {
        _S721 = (_S716[int(3)]) != 0.0f;
    }
    else
    {
        _S721 = false;
    }
    float _S722;
    float _S723;
    float _S724;
    float _S725;
    float _S726;
    float _S727;
    float _S728;
    float _S729;
    float _S730;
    float _S731;
    float _S732;
    float _S733;
    float _S734;
    float _S735;
    float _S736;
    if(_S721)
    {
        float _S737 = _S716[int(2)] * _S716[int(3)];
        float _S738 = _S716[int(18)] * _S716[int(18)];
        float _S739 = _S716[int(6)] - _S737 / _S716[int(18)];
        float _S740 = _S716[int(2)] * _S716[int(2)];
        float _S741 = _S716[int(4)] - _S740 / _S716[int(18)];
        float _S742 = _S716[int(3)] * _S716[int(3)];
        float _S743 = _S716[int(5)] - _S742 / _S716[int(18)];
        float _S744 = _S741 * _S743 + 1.0f;
        float _S745 = s_primal_ctx_max_0(9.999999960041972e-13f, _S744);
        float _S746 = s_primal_ctx_sqrt_0(_S745);
        float _S747 = _S746 * _S746;
        float _S748 = 1.0f - _S739 / _S746;
        _S722 = (*weights_4)[int(2)];
        _S723 = _S748;
        _S724 = _S747;
        _S725 = _S739;
        _S726 = _S746;
        _S727 = _S745;
        _S728 = _S744;
        _S729 = _S741;
        _S730 = _S743;
        _S731 = _S738;
        _S732 = _S742;
        _S733 = _S716[int(3)];
        _S734 = _S740;
        _S735 = _S716[int(2)];
        _S736 = _S737;
    }
    else
    {
        _S722 = 0.0f;
        _S723 = 0.0f;
        _S724 = 0.0f;
        _S725 = 0.0f;
        _S726 = 0.0f;
        _S727 = 0.0f;
        _S728 = 0.0f;
        _S729 = 0.0f;
        _S730 = 0.0f;
        _S731 = 0.0f;
        _S732 = 0.0f;
        _S733 = 0.0f;
        _S734 = 0.0f;
        _S735 = 0.0f;
        _S736 = 0.0f;
    }
    float _S749 = s_primal_ctx_max_0(_S716[int(19)], 1.0f);
    float _S750 = _S749 * _S749;
    float _S751 = s_primal_ctx_max_0(_S716[int(20)], 1.0f);
    float _S752 = _S751 * _S751;
    float _S753 = float((I32_max((int((_S716[int(19)]) > 0.5f) + int((_S716[int(20)]) > 0.5f)), (int(1)))));
    float _S754 = _S716[int(9)] + _S716[int(10)];
    float _S755 = s_primal_ctx_max_0(_S716[int(22)], 1.0f);
    float _S756 = _S755 * _S755;
    float _S757 = s_primal_ctx_max_0(_S716[int(21)], 1.0f);
    float _S758 = _S757 * _S757;
    float _S759 = s_primal_ctx_max_0(_S716[int(16)], 1.0f);
    float _S760 = _S759 * _S759;
    float _S761 = (*_s_dOut_5)[int(9)] / _S760;
    float _S762 = _S759 * _S761;
    float _S763 = (*_s_dOut_5)[int(8)] / _S760;
    float _S764 = _S759 * _S763;
    float _S765 = (*_s_dOut_5)[int(7)] / _S760;
    float _S766 = _S759 * _S765;
    float _S767 = (*_s_dOut_5)[int(6)] / _S760;
    float _S768 = _S759 * _S767;
    float _S769 = _S716[int(15)] * - _S761 + _S716[int(14)] * - _S763 + _S716[int(13)] * - _S765 + _S716[int(12)] * - _S767;
    DiffPair_float_0 _S770;
    (&_S770)->primal_0 = _S716[int(16)];
    (&_S770)->differential_0 = 0.0f;
    DiffPair_float_0 _S771;
    (&_S771)->primal_0 = 1.0f;
    (&_S771)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S770, &_S771, _S769);
    float _S772 = (*_s_dOut_5)[int(5)] / _S758;
    float _S773 = _S716[int(11)] * - _S772;
    float _S774 = _S757 * _S772;
    DiffPair_float_0 _S775;
    (&_S775)->primal_0 = _S716[int(21)];
    (&_S775)->differential_0 = 0.0f;
    DiffPair_float_0 _S776;
    (&_S776)->primal_0 = 1.0f;
    (&_S776)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S775, &_S776, _S773);
    float _S777 = (*_s_dOut_5)[int(4)] / _S756;
    float _S778 = _S754 * - _S777;
    float _S779 = _S755 * _S777;
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = _S716[int(22)];
    (&_S780)->differential_0 = 0.0f;
    DiffPair_float_0 _S781;
    (&_S781)->primal_0 = 1.0f;
    (&_S781)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S780, &_S781, _S778);
    float _S782 = (*_s_dOut_5)[int(3)] / _S753;
    float _S783 = _S782 / _S752;
    float _S784 = _S716[int(8)] * - _S783;
    float _S785 = _S751 * _S783;
    DiffPair_float_0 _S786;
    (&_S786)->primal_0 = _S716[int(20)];
    (&_S786)->differential_0 = 0.0f;
    DiffPair_float_0 _S787;
    (&_S787)->primal_0 = 1.0f;
    (&_S787)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S786, &_S787, _S784);
    float _S788 = _S782 / _S750;
    float _S789 = _S716[int(7)] * - _S788;
    float _S790 = _S749 * _S788;
    DiffPair_float_0 _S791;
    (&_S791)->primal_0 = _S716[int(19)];
    (&_S791)->differential_0 = 0.0f;
    DiffPair_float_0 _S792;
    (&_S792)->primal_0 = 1.0f;
    (&_S792)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S791, &_S792, _S789);
    FixedArray<float, 23>  _S793;
    _S793[int(0)] = 0.0f;
    _S793[int(1)] = 0.0f;
    _S793[int(2)] = 0.0f;
    _S793[int(3)] = 0.0f;
    _S793[int(4)] = 0.0f;
    _S793[int(5)] = 0.0f;
    _S793[int(6)] = 0.0f;
    _S793[int(7)] = 0.0f;
    _S793[int(8)] = 0.0f;
    _S793[int(9)] = 0.0f;
    _S793[int(10)] = 0.0f;
    _S793[int(11)] = 0.0f;
    _S793[int(12)] = 0.0f;
    _S793[int(13)] = 0.0f;
    _S793[int(14)] = 0.0f;
    _S793[int(15)] = 0.0f;
    _S793[int(16)] = 0.0f;
    _S793[int(17)] = 0.0f;
    _S793[int(18)] = 0.0f;
    _S793[int(19)] = 0.0f;
    _S793[int(20)] = 0.0f;
    _S793[int(21)] = 0.0f;
    _S793[int(22)] = 0.0f;
    _S793[int(15)] = _S762;
    _S793[int(14)] = _S764;
    _S793[int(13)] = _S766;
    _S793[int(16)] = _S770.differential_0;
    _S793[int(12)] = _S768;
    _S793[int(21)] = _S775.differential_0;
    _S793[int(11)] = _S774;
    _S793[int(22)] = _S780.differential_0;
    _S793[int(10)] = _S779;
    _S793[int(9)] = _S779;
    _S793[int(20)] = _S786.differential_0;
    _S793[int(8)] = _S785;
    _S793[int(19)] = _S791.differential_0;
    _S793[int(7)] = _S790;
    float _S794 = _S793[int(0)];
    float _S795 = _S793[int(1)];
    float _S796 = _S793[int(2)];
    float _S797 = _S793[int(3)];
    float _S798 = _S793[int(4)];
    float _S799 = _S793[int(5)];
    float _S800 = _S793[int(6)];
    float _S801 = _S793[int(7)];
    float _S802 = _S793[int(8)];
    float _S803 = _S793[int(9)];
    float _S804 = _S793[int(10)];
    float _S805 = _S793[int(11)];
    float _S806 = _S793[int(12)];
    float _S807 = _S793[int(13)];
    float _S808 = _S793[int(14)];
    float _S809 = _S793[int(15)];
    float _S810 = _S793[int(16)];
    float _S811 = _S793[int(17)];
    float _S812 = _S793[int(18)];
    float _S813 = _S793[int(19)];
    float _S814 = _S793[int(20)];
    float _S815 = _S793[int(21)];
    float _S816 = _S793[int(22)];
    FixedArray<float, 23>  _S817;
    if(_S721)
    {
        float _S818 = _S722 * (*_s_dOut_5)[int(2)];
        DiffPair_float_0 _S819;
        (&_S819)->primal_0 = _S723;
        (&_S819)->differential_0 = 0.0f;
        DiffPair_float_0 _S820;
        (&_S820)->primal_0 = 0.0f;
        (&_S820)->differential_0 = 0.0f;
        DiffPair_float_0 _S821;
        (&_S821)->primal_0 = 2.0f;
        (&_S821)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S819, &_S820, &_S821, _S818);
        float _S822 = - _S819.differential_0 / _S724;
        float _S823 = _S725 * - _S822;
        float _S824 = _S726 * _S822;
        DiffPair_float_0 _S825;
        (&_S825)->primal_0 = _S727;
        (&_S825)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S825, _S823);
        DiffPair_float_0 _S826;
        (&_S826)->primal_0 = 9.999999960041972e-13f;
        (&_S826)->differential_0 = 0.0f;
        DiffPair_float_0 _S827;
        (&_S827)->primal_0 = _S728;
        (&_S827)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S826, &_S827, _S825.differential_0);
        float _S828 = _S729 * _S827.differential_0;
        float _S829 = _S730 * _S827.differential_0;
        float _S830 = - _S828 / _S731;
        float _S831 = _S733 * (_S716[int(18)] * _S830);
        float _S832 = - _S829 / _S731;
        float _S833 = _S735 * (_S716[int(18)] * _S832);
        float _S834 = - _S824 / _S731;
        float _S835 = _S716[int(18)] * _S834;
        float _S836 = _S831 + _S831 + _S735 * _S835;
        float _S837 = _S833 + _S833 + _S733 * _S835;
        float _S838 = _S732 * - _S830 + _S734 * - _S832 + _S736 * - _S834;
        FixedArray<float, 23>  _S839;
        _S839[int(0)] = 0.0f;
        _S839[int(1)] = 0.0f;
        _S839[int(2)] = 0.0f;
        _S839[int(3)] = 0.0f;
        _S839[int(4)] = 0.0f;
        _S839[int(5)] = 0.0f;
        _S839[int(6)] = 0.0f;
        _S839[int(7)] = 0.0f;
        _S839[int(8)] = 0.0f;
        _S839[int(9)] = 0.0f;
        _S839[int(10)] = 0.0f;
        _S839[int(11)] = 0.0f;
        _S839[int(12)] = 0.0f;
        _S839[int(13)] = 0.0f;
        _S839[int(14)] = 0.0f;
        _S839[int(15)] = 0.0f;
        _S839[int(16)] = 0.0f;
        _S839[int(17)] = 0.0f;
        _S839[int(18)] = 0.0f;
        _S839[int(19)] = 0.0f;
        _S839[int(20)] = 0.0f;
        _S839[int(21)] = 0.0f;
        _S839[int(22)] = 0.0f;
        _S839[int(5)] = _S828;
        _S839[int(4)] = _S829;
        _S839[int(3)] = _S836;
        _S839[int(2)] = _S837;
        _S839[int(6)] = _S824;
        float _S840 = _S795 + _S839[int(1)];
        float _S841 = _S796 + _S839[int(2)];
        float _S842 = _S797 + _S839[int(3)];
        float _S843 = _S798 + _S839[int(4)];
        float _S844 = _S799 + _S839[int(5)];
        float _S845 = _S800 + _S839[int(6)];
        float _S846 = _S801 + _S839[int(7)];
        float _S847 = _S802 + _S839[int(8)];
        float _S848 = _S803 + _S839[int(9)];
        float _S849 = _S804 + _S839[int(10)];
        float _S850 = _S805 + _S839[int(11)];
        float _S851 = _S806 + _S839[int(12)];
        float _S852 = _S807 + _S839[int(13)];
        float _S853 = _S808 + _S839[int(14)];
        float _S854 = _S809 + _S839[int(15)];
        float _S855 = _S810 + _S839[int(16)];
        float _S856 = _S811 + _S839[int(17)];
        float _S857 = _S812 + _S839[int(18)];
        float _S858 = _S813 + _S839[int(19)];
        float _S859 = _S814 + _S839[int(20)];
        float _S860 = _S815 + _S839[int(21)];
        float _S861 = _S816 + _S839[int(22)];
        _S817[int(0)] = _S794 + _S839[int(0)];
        _S817[int(1)] = _S840;
        _S817[int(2)] = _S841;
        _S817[int(3)] = _S842;
        _S817[int(4)] = _S843;
        _S817[int(5)] = _S844;
        _S817[int(6)] = _S845;
        _S817[int(7)] = _S846;
        _S817[int(8)] = _S847;
        _S817[int(9)] = _S848;
        _S817[int(10)] = _S849;
        _S817[int(11)] = _S850;
        _S817[int(12)] = _S851;
        _S817[int(13)] = _S852;
        _S817[int(14)] = _S853;
        _S817[int(15)] = _S854;
        _S817[int(16)] = _S855;
        _S817[int(17)] = _S856;
        _S817[int(18)] = _S857;
        _S817[int(19)] = _S858;
        _S817[int(20)] = _S859;
        _S817[int(21)] = _S860;
        _S817[int(22)] = _S861;
        _S722 = _S838;
    }
    else
    {
        _S817[int(0)] = _S794;
        _S817[int(1)] = _S795;
        _S817[int(2)] = _S796;
        _S817[int(3)] = _S797;
        _S817[int(4)] = _S798;
        _S817[int(5)] = _S799;
        _S817[int(6)] = _S800;
        _S817[int(7)] = _S801;
        _S817[int(8)] = _S802;
        _S817[int(9)] = _S803;
        _S817[int(10)] = _S804;
        _S817[int(11)] = _S805;
        _S817[int(12)] = _S806;
        _S817[int(13)] = _S807;
        _S817[int(14)] = _S808;
        _S817[int(15)] = _S809;
        _S817[int(16)] = _S810;
        _S817[int(17)] = _S811;
        _S817[int(18)] = _S812;
        _S817[int(19)] = _S813;
        _S817[int(20)] = _S814;
        _S817[int(21)] = _S815;
        _S817[int(22)] = _S816;
        _S722 = 0.0f;
    }
    if(_S720)
    {
        FixedArray<float, 23>  _S862;
        _S862[int(0)] = 0.0f;
        _S862[int(1)] = 0.0f;
        _S862[int(2)] = 0.0f;
        _S862[int(3)] = 0.0f;
        _S862[int(4)] = 0.0f;
        _S862[int(5)] = 0.0f;
        _S862[int(6)] = 0.0f;
        _S862[int(7)] = 0.0f;
        _S862[int(8)] = 0.0f;
        _S862[int(9)] = 0.0f;
        _S862[int(10)] = 0.0f;
        _S862[int(11)] = 0.0f;
        _S862[int(12)] = 0.0f;
        _S862[int(13)] = 0.0f;
        _S862[int(14)] = 0.0f;
        _S862[int(15)] = 0.0f;
        _S862[int(16)] = 0.0f;
        _S862[int(17)] = 0.0f;
        _S862[int(18)] = 0.0f;
        _S862[int(19)] = 0.0f;
        _S862[int(20)] = 0.0f;
        _S862[int(21)] = 0.0f;
        _S862[int(22)] = 0.0f;
        _S862[int(3)] = 0.0f;
        float _S863 = _S817[int(1)] + _S862[int(1)];
        float _S864 = _S817[int(2)] + _S862[int(2)];
        float _S865 = _S817[int(3)] + _S862[int(3)];
        float _S866 = _S817[int(4)] + _S862[int(4)];
        float _S867 = _S817[int(5)] + _S862[int(5)];
        float _S868 = _S817[int(6)] + _S862[int(6)];
        float _S869 = _S817[int(7)] + _S862[int(7)];
        float _S870 = _S817[int(8)] + _S862[int(8)];
        float _S871 = _S817[int(9)] + _S862[int(9)];
        float _S872 = _S817[int(10)] + _S862[int(10)];
        float _S873 = _S817[int(11)] + _S862[int(11)];
        float _S874 = _S817[int(12)] + _S862[int(12)];
        float _S875 = _S817[int(13)] + _S862[int(13)];
        float _S876 = _S817[int(14)] + _S862[int(14)];
        float _S877 = _S817[int(15)] + _S862[int(15)];
        float _S878 = _S817[int(16)] + _S862[int(16)];
        float _S879 = _S817[int(17)] + _S862[int(17)];
        float _S880 = _S817[int(18)] + _S862[int(18)];
        float _S881 = _S817[int(19)] + _S862[int(19)];
        float _S882 = _S817[int(20)] + _S862[int(20)];
        float _S883 = _S817[int(21)] + _S862[int(21)];
        float _S884 = _S817[int(22)] + _S862[int(22)];
        _S817[int(0)] = _S817[int(0)] + _S862[int(0)];
        _S817[int(1)] = _S863;
        _S817[int(2)] = _S864;
        _S817[int(3)] = _S865;
        _S817[int(4)] = _S866;
        _S817[int(5)] = _S867;
        _S817[int(6)] = _S868;
        _S817[int(7)] = _S869;
        _S817[int(8)] = _S870;
        _S817[int(9)] = _S871;
        _S817[int(10)] = _S872;
        _S817[int(11)] = _S873;
        _S817[int(12)] = _S874;
        _S817[int(13)] = _S875;
        _S817[int(14)] = _S876;
        _S817[int(15)] = _S877;
        _S817[int(16)] = _S878;
        _S817[int(17)] = _S879;
        _S817[int(18)] = _S880;
        _S817[int(19)] = _S881;
        _S817[int(20)] = _S882;
        _S817[int(21)] = _S883;
        _S817[int(22)] = _S884;
    }
    float _S885 = -10.0f * (*_s_dOut_5)[int(1)];
    DiffPair_float_0 _S886;
    (&_S886)->primal_0 = _S719;
    (&_S886)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S886, _S885);
    float _S887 = _S886.differential_0 / _S718;
    float _S888 = _S717 * _S887;
    float _S889 = (*_s_dOut_5)[int(0)] / _S718;
    float _S890 = _S717 * _S889;
    float _S891 = _S716[int(1)] * - _S887 + _S716[int(0)] * - _S889;
    DiffPair_float_0 _S892;
    (&_S892)->primal_0 = _S716[int(17)];
    (&_S892)->differential_0 = 0.0f;
    DiffPair_float_0 _S893;
    (&_S893)->primal_0 = 1.0f;
    (&_S893)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S892, &_S893, _S891);
    FixedArray<float, 23>  _S894;
    _S894[int(0)] = 0.0f;
    _S894[int(1)] = 0.0f;
    _S894[int(2)] = 0.0f;
    _S894[int(3)] = 0.0f;
    _S894[int(4)] = 0.0f;
    _S894[int(5)] = 0.0f;
    _S894[int(6)] = 0.0f;
    _S894[int(7)] = 0.0f;
    _S894[int(8)] = 0.0f;
    _S894[int(9)] = 0.0f;
    _S894[int(10)] = 0.0f;
    _S894[int(11)] = 0.0f;
    _S894[int(12)] = 0.0f;
    _S894[int(13)] = 0.0f;
    _S894[int(14)] = 0.0f;
    _S894[int(15)] = 0.0f;
    _S894[int(16)] = 0.0f;
    _S894[int(17)] = 0.0f;
    _S894[int(18)] = 0.0f;
    _S894[int(19)] = 0.0f;
    _S894[int(20)] = 0.0f;
    _S894[int(21)] = 0.0f;
    _S894[int(22)] = 0.0f;
    _S894[int(18)] = _S722;
    _S894[int(1)] = _S888;
    _S894[int(17)] = _S892.differential_0;
    _S894[int(0)] = _S890;
    FixedArray<float, 23>  _S895 = {
        _S817[int(0)] + _S894[int(0)], _S817[int(1)] + _S894[int(1)], _S817[int(2)] + _S894[int(2)], _S817[int(3)] + _S894[int(3)], _S817[int(4)] + _S894[int(4)], _S817[int(5)] + _S894[int(5)], _S817[int(6)] + _S894[int(6)], _S817[int(7)] + _S894[int(7)], _S817[int(8)] + _S894[int(8)], _S817[int(9)] + _S894[int(9)], _S817[int(10)] + _S894[int(10)], _S817[int(11)] + _S894[int(11)], _S817[int(12)] + _S894[int(12)], _S817[int(13)] + _S894[int(13)], _S817[int(14)] + _S894[int(14)], _S817[int(15)] + _S894[int(15)], _S817[int(16)] + _S894[int(16)], _S817[int(17)] + _S894[int(17)], _S817[int(18)] + _S894[int(18)], _S817[int(19)] + _S894[int(19)], _S817[int(20)] + _S894[int(20)], _S817[int(21)] + _S894[int(21)], _S817[int(22)] + _S894[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S895;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S896, FixedArray<float, 11>  * _S897, FixedArray<float, 10>  * _S898)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S896, _S897, _S898);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  * raw_losses_1, FixedArray<float, 11>  * weights_5, FixedArray<float, 10>  * v_losses_1, FixedArray<float, 23>  * _S899)
{
    FixedArray<float, 23>  _S900 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = *raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S900;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, weights_5, v_losses_1);
    *_S899 = (&dp_raw_losses_0)->differential_0;
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

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_17)
{
    DiffPair_float_0 left_dp_3;
    (&left_dp_3)->primal_0 = (*dpx_16).primal_0.x;
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
    (&left_dp_4)->primal_0 = (*dpx_16).primal_0.y;
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
    (&left_dp_5)->primal_0 = (*dpx_16).primal_0.z;
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
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = left_d_result_2;
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

inline __device__ void s_bwd_prop_clamp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S901, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S902, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S903, float3  _S904)
{
    _d_clamp_vector_0(_S901, _S902, _S903, _S904);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_6)
{
    float _S905 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S906 = make_float3 (_S905);
    float3  _S907 = make_float3 (0.0f);
    float3  _S908 = make_float3 (1.0f);
    float3  _S909 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S910;
    (&_S910)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S905) * (*dpbackground_0).primal_0;
    (&_S910)->differential_0 = _S909;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S911;
    (&_S911)->primal_0 = _S907;
    (&_S911)->differential_0 = _S909;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S912;
    (&_S912)->primal_0 = _S908;
    (&_S912)->differential_0 = _S909;
    s_bwd_prop_clamp_1(&_S910, &_S911, &_S912, _s_dOut_6);
    float3  _S913 = _S906 * _S910.differential_0;
    float3  _S914 = (*dpbackground_0).primal_0 * _S910.differential_0;
    float _S915 = - (_S914.x + _S914.y + _S914.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S913;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S915;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S910.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S916, DiffPair_float_0 * _S917, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S918, float3  _S919)
{
    s_bwd_prop_blend_background_0(_S916, _S917, _S918, _S919);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S920 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S920;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S920;
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
    float3  _S921 = make_float3 (1.0f - t_4);
    float3  _S922 = make_float3 (0.0f);
    float3  _S923 = make_float3 (t_4) * _s_dOut_7;
    float3  _S924 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S925;
    (&_S925)->primal_0 = s_primal_ctx_max_1((*dprgb_1).primal_0, _S922) + make_float3 (0.00392156885936856f);
    (&_S925)->differential_0 = _S924;
    s_bwd_prop_log_1(&_S925, _S923);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S926;
    (&_S926)->primal_0 = (*dprgb_1).primal_0;
    (&_S926)->differential_0 = _S924;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S927;
    (&_S927)->primal_0 = _S922;
    (&_S927)->differential_0 = _S924;
    s_bwd_prop_max_1(&_S926, &_S927, _S925.differential_0);
    float3  _S928 = _S926.differential_0 + _S921 * _s_dOut_7;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S928;
    return;
}

inline __device__ void s_bwd_log_map_image_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S929, float _S930, float3  _S931)
{
    s_bwd_prop_log_map_image_0(_S929, _S930, _S931);
    return;
}

inline __device__ float3  log_map_image_bwd(float3  rgb_3, float t_5, float3  v_out_rgb_1)
{
    float3  _S932 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S932;
    s_bwd_log_map_image_0(&p_rgb_1, t_5, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_7, bool is_fisheye_4, bool is_ray_depth_0)
{
    float2  _S933 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_7 = _S933;
    bool _S934 = undistort_point_0(_S933, dist_coeffs_7, int(12), &uv_7);
    if(!_S934)
    {
        int3  _S935 = make_int3 (int(0));
        float3  _S936 = make_float3 ((float)_S935.x, (float)_S935.y, (float)_S935.z);
        return _S936;
    }
    float3  raydir_5;
    if(is_fisheye_4)
    {
        float theta_4 = length_1(uv_7);
        float3  raydir_6 = make_float3 ((uv_7 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_7 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
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
        float3  raydir_7 = make_float3 (uv_7.x, uv_7.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_5 = normalize_0(raydir_7);
        }
        else
        {
            raydir_5 = raydir_7;
        }
    }
    return raydir_5;
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S937, float3  _S938)
{
    return cross_0(_S937, _S938);
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_17, float _s_dOut_8)
{
    float _S939 = (*dpx_17).primal_0.x;
    float _S940 = (*dpx_17).primal_0.y;
    float _S941 = (*dpx_17).primal_0.z;
    DiffPair_float_0 _S942;
    (&_S942)->primal_0 = _S939 * _S939 + _S940 * _S940 + _S941 * _S941;
    (&_S942)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S942, _s_dOut_8);
    float _S943 = (*dpx_17).primal_0.z * _S942.differential_0;
    float _S944 = _S943 + _S943;
    float _S945 = (*dpx_17).primal_0.y * _S942.differential_0;
    float _S946 = _S945 + _S945;
    float _S947 = (*dpx_17).primal_0.x * _S942.differential_0;
    float _S948 = _S947 + _S947;
    float3  _S949 = make_float3 (0.0f);
    *&((&_S949)->z) = _S944;
    *&((&_S949)->y) = _S946;
    *&((&_S949)->x) = _S948;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S949;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S950, float _S951)
{
    s_bwd_prop_length_impl_1(_S950, _S951);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S952, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S953, float3  _S954)
{
    _d_cross_0(_S952, _S953, _S954);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_9)
{
    float3  _S955 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S956 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S957 = s_primal_ctx_cross_0(dx_0, _S956);
    bool _S958 = (s_primal_ctx_dot_0(_S957, _S957)) != 0.0f;
    float3  _S959;
    float3  _S960;
    if(_S958)
    {
        float _S961 = length_2(_S957);
        float3  _S962 = make_float3 (_S961);
        _S959 = make_float3 (_S961 * _S961);
        _S960 = _S962;
    }
    else
    {
        _S959 = _S955;
        _S960 = _S955;
    }
    if(_S958)
    {
        float3  _S963 = _s_dOut_9 / _S959;
        float3  _S964 = _S957 * - _S963;
        float3  _S965 = _S960 * _S963;
        float _S966 = _S964.x + _S964.y + _S964.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S967;
        (&_S967)->primal_0 = _S957;
        (&_S967)->differential_0 = _S955;
        s_bwd_length_impl_1(&_S967, _S966);
        _S959 = _S965 + _S967.differential_0;
    }
    else
    {
        _S959 = _s_dOut_9;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S968;
    (&_S968)->primal_0 = _S957;
    (&_S968)->differential_0 = _S955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S969;
    (&_S969)->primal_0 = _S957;
    (&_S969)->differential_0 = _S955;
    s_bwd_prop_dot_0(&_S968, &_S969, 0.0f);
    float3  _S970 = _S969.differential_0 + _S968.differential_0 + _S959;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S971;
    (&_S971)->primal_0 = dx_0;
    (&_S971)->differential_0 = _S955;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S972;
    (&_S972)->primal_0 = _S956;
    (&_S972)->differential_0 = _S955;
    s_bwd_prop_cross_0(&_S971, &_S972, _S970);
    float3  s_diff_dy_T_0 = - _S972.differential_0;
    float3  _S973 = - s_diff_dy_T_0;
    float3  _S974 = - _S971.differential_0;
    FixedArray<float3 , 4>  _S975;
    _S975[int(0)] = _S955;
    _S975[int(1)] = _S955;
    _S975[int(2)] = _S955;
    _S975[int(3)] = _S955;
    _S975[int(2)] = _S973;
    _S975[int(3)] = s_diff_dy_T_0;
    _S975[int(0)] = _S974;
    _S975[int(1)] = _S971.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S975;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S976, float3  _S977)
{
    s_bwd_prop_points_to_normal_0(_S976, _S977);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  * points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S978 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = *points_1;
    (&dp_points_0)->differential_0 = _S978;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_5, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S979 = float2 {intrins_1.z, intrins_1.w};
    float2  _S980 = float2 {intrins_1.x, intrins_1.y};
    float2  _S981 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S979) / _S980;
    float2  uv_8 = _S981;
    bool _S982 = undistort_point_0(_S981, dist_coeffs_8, int(12), &uv_8);
    if(!_S982)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_8;
    if(is_fisheye_5)
    {
        float theta_5 = length_1(uv_8);
        float3  raydir_9 = make_float3 ((uv_8 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_8 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        if(!is_ray_depth_1)
        {
            raydir_8 = raydir_9 / make_float3 (raydir_9.z);
        }
        else
        {
            raydir_8 = raydir_9;
        }
    }
    else
    {
        float3  raydir_10 = make_float3 (uv_8.x, uv_8.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_8 = normalize_0(raydir_10);
        }
        else
        {
            raydir_8 = raydir_10;
        }
    }
    points_2[int(0)] = make_float3 (depths_0.x) * raydir_8;
    float2  _S983 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S979) / _S980;
    float2  uv_9 = _S983;
    bool _S984 = undistort_point_0(_S983, dist_coeffs_8, int(12), &uv_9);
    if(!_S984)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_6 = length_1(uv_9);
        float3  raydir_11 = make_float3 ((uv_9 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (uv_9 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
        if(!is_ray_depth_1)
        {
            raydir_8 = raydir_11 / make_float3 (raydir_11.z);
        }
        else
        {
            raydir_8 = raydir_11;
        }
    }
    else
    {
        float3  raydir_12 = make_float3 (uv_9.x, uv_9.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_8 = normalize_0(raydir_12);
        }
        else
        {
            raydir_8 = raydir_12;
        }
    }
    points_2[int(1)] = make_float3 (depths_0.y) * raydir_8;
    float2  _S985 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S979) / _S980;
    float2  uv_10 = _S985;
    bool _S986 = undistort_point_0(_S985, dist_coeffs_8, int(12), &uv_10);
    if(!_S986)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_7 = length_1(uv_10);
        float3  raydir_13 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (uv_10 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        if(!is_ray_depth_1)
        {
            raydir_8 = raydir_13 / make_float3 (raydir_13.z);
        }
        else
        {
            raydir_8 = raydir_13;
        }
    }
    else
    {
        float3  raydir_14 = make_float3 (uv_10.x, uv_10.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_8 = normalize_0(raydir_14);
        }
        else
        {
            raydir_8 = raydir_14;
        }
    }
    points_2[int(2)] = make_float3 (depths_0.z) * raydir_8;
    float2  _S987 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S979) / _S980;
    float2  uv_11 = _S987;
    bool _S988 = undistort_point_0(_S987, dist_coeffs_8, int(12), &uv_11);
    if(!_S988)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_8 = length_1(uv_11);
        float3  raydir_15 = make_float3 ((uv_11 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).x, (uv_11 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).y, (F32_cos((theta_8))));
        if(!is_ray_depth_1)
        {
            raydir_8 = raydir_15 / make_float3 (raydir_15.z);
        }
        else
        {
            raydir_8 = raydir_15;
        }
    }
    else
    {
        float3  raydir_16 = make_float3 (uv_11.x, uv_11.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_8 = normalize_0(raydir_16);
        }
        else
        {
            raydir_8 = raydir_16;
        }
    }
    points_2[int(3)] = make_float3 (depths_0.w) * raydir_8;
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
    float2  _S989;
    bool _S990;
    float2  _S991;
    bool _S992;
    float2  _S993;
    bool _S994;
    float2  _S995;
    bool _S996;
};

inline __device__ float s_primal_ctx_sin_0(float _S997)
{
    return (F32_sin((_S997)));
}

inline __device__ float s_primal_ctx_cos_0(float _S998)
{
    return (F32_cos((_S998)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, bool is_fisheye_6, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S999 = make_float2 (0.0f);
    _s_diff_ctx_0->_S989 = _S999;
    _s_diff_ctx_0->_S990 = false;
    _s_diff_ctx_0->_S991 = _S999;
    _s_diff_ctx_0->_S992 = false;
    _s_diff_ctx_0->_S993 = _S999;
    _s_diff_ctx_0->_S994 = false;
    _s_diff_ctx_0->_S995 = _S999;
    _s_diff_ctx_0->_S996 = false;
    _s_diff_ctx_0->_S991 = _S999;
    _s_diff_ctx_0->_S992 = false;
    _s_diff_ctx_0->_S993 = _S999;
    _s_diff_ctx_0->_S994 = false;
    _s_diff_ctx_0->_S995 = _S999;
    _s_diff_ctx_0->_S996 = false;
    float3  _S1000 = make_float3 (0.0f);
    float2  _S1001 = float2 {intrins_2.z, intrins_2.w};
    float2  _S1002 = float2 {intrins_2.x, intrins_2.y};
    float2  _S1003 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S1001) / _S1002;
    float2  _S1004 = _S1003;
    bool _S1005 = undistort_point_0(_S1003, dist_coeffs_9, int(12), &_S1004);
    _s_diff_ctx_0->_S989 = _S1004;
    _s_diff_ctx_0->_S990 = _S1005;
    float2  uv_12 = _S1004;
    bool _S1006 = !_S1005;
    float3  normal_4;
    if(_S1006)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S1007 = !_S1006;
    int _S1008;
    FixedArray<float3 , 4>  points_3;
    if(_S1007)
    {
        float3  raydir_17;
        if(is_fisheye_6)
        {
            float _S1009 = length_1(uv_12);
            float3  raydir_18 = make_float3 ((uv_12 / make_float2 (s_primal_ctx_max_0(_S1009, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1009))).x, (uv_12 / make_float2 (s_primal_ctx_max_0(_S1009, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1009))).y, s_primal_ctx_cos_0(_S1009));
            if(!is_ray_depth_2)
            {
                raydir_17 = raydir_18 / make_float3 (raydir_18.z);
            }
            else
            {
                raydir_17 = raydir_18;
            }
        }
        else
        {
            float3  raydir_19 = make_float3 (uv_12.x, uv_12.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_17 = normalize_0(raydir_19);
            }
            else
            {
                raydir_17 = raydir_19;
            }
        }
        float3  _S1010 = make_float3 (dpdepths_0.x) * raydir_17;
        float2  _S1011 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S1001) / _S1002;
        float2  _S1012 = _S1011;
        bool _S1013 = undistort_point_0(_S1011, dist_coeffs_9, int(12), &_S1012);
        _s_diff_ctx_0->_S991 = _S1012;
        _s_diff_ctx_0->_S992 = _S1013;
        float2  uv_13 = _S1012;
        bool _S1014 = !_S1013;
        if(_S1014)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S1015 = !_S1014;
        if(_S1015)
        {
            if(is_fisheye_6)
            {
                float _S1016 = length_1(uv_13);
                float3  raydir_20 = make_float3 ((uv_13 / make_float2 (s_primal_ctx_max_0(_S1016, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1016))).x, (uv_13 / make_float2 (s_primal_ctx_max_0(_S1016, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1016))).y, s_primal_ctx_cos_0(_S1016));
                if(!is_ray_depth_2)
                {
                    raydir_17 = raydir_20 / make_float3 (raydir_20.z);
                }
                else
                {
                    raydir_17 = raydir_20;
                }
            }
            else
            {
                float3  raydir_21 = make_float3 (uv_13.x, uv_13.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_17 = normalize_0(raydir_21);
                }
                else
                {
                    raydir_17 = raydir_21;
                }
            }
            float3  _S1017 = make_float3 (dpdepths_0.y) * raydir_17;
            _S1008 = int(2);
            points_3[int(0)] = _S1010;
            points_3[int(1)] = _S1017;
            points_3[int(2)] = _S1000;
            points_3[int(3)] = _S1000;
        }
        else
        {
            _S1008 = int(0);
            points_3[int(0)] = _S1010;
            points_3[int(1)] = _S1000;
            points_3[int(2)] = _S1000;
            points_3[int(3)] = _S1000;
        }
        bool _runFlag_0;
        if(_S1008 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S1007;
            _S1008 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S1018 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S1001) / _S1002;
            float2  _S1019 = _S1018;
            bool _S1020 = undistort_point_0(_S1018, dist_coeffs_9, int(12), &_S1019);
            _s_diff_ctx_0->_S993 = _S1019;
            _s_diff_ctx_0->_S994 = _S1020;
            float2  uv_14 = _S1019;
            if(!_S1020)
            {
                float3  _S1021 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S1008 = int(0);
                normal_4 = _S1021;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S1022 = length_1(uv_14);
                    float3  raydir_22 = make_float3 ((uv_14 / make_float2 (s_primal_ctx_max_0(_S1022, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1022))).x, (uv_14 / make_float2 (s_primal_ctx_max_0(_S1022, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1022))).y, s_primal_ctx_cos_0(_S1022));
                    if(!is_ray_depth_2)
                    {
                        raydir_17 = raydir_22 / make_float3 (raydir_22.z);
                    }
                    else
                    {
                        raydir_17 = raydir_22;
                    }
                }
                else
                {
                    float3  raydir_23 = make_float3 (uv_14.x, uv_14.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_17 = normalize_0(raydir_23);
                    }
                    else
                    {
                        raydir_17 = raydir_23;
                    }
                }
                points_3[int(2)] = make_float3 (dpdepths_0.z) * raydir_17;
                float2  _S1023 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S1001) / _S1002;
                float2  _S1024 = _S1023;
                bool _S1025 = undistort_point_0(_S1023, dist_coeffs_9, int(12), &_S1024);
                _s_diff_ctx_0->_S995 = _S1024;
                _s_diff_ctx_0->_S996 = _S1025;
                float2  uv_15 = _S1024;
                bool _S1026 = !_S1025;
                if(_S1026)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S1027 = !_S1026;
                int _S1028;
                if(_S1027)
                {
                    if(is_fisheye_6)
                    {
                        float _S1029 = length_1(uv_15);
                        float3  raydir_24 = make_float3 ((uv_15 / make_float2 (s_primal_ctx_max_0(_S1029, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1029))).x, (uv_15 / make_float2 (s_primal_ctx_max_0(_S1029, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1029))).y, s_primal_ctx_cos_0(_S1029));
                        if(!is_ray_depth_2)
                        {
                            raydir_17 = raydir_24 / make_float3 (raydir_24.z);
                        }
                        else
                        {
                            raydir_17 = raydir_24;
                        }
                    }
                    else
                    {
                        float3  raydir_25 = make_float3 (uv_15.x, uv_15.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_17 = normalize_0(raydir_25);
                        }
                        else
                        {
                            raydir_17 = raydir_25;
                        }
                    }
                    points_3[int(3)] = make_float3 (dpdepths_0.w) * raydir_17;
                    _S1028 = int(2);
                }
                else
                {
                    _S1028 = int(0);
                }
                if(_S1028 != int(2))
                {
                    _runFlag_0 = false;
                    _S1008 = _S1028;
                }
                if(_runFlag_0)
                {
                    _S1008 = int(1);
                }
            }
        }
    }
    else
    {
        _S1008 = int(0);
        points_3[int(0)] = _S1000;
        points_3[int(1)] = _S1000;
        points_3[int(2)] = _S1000;
        points_3[int(3)] = _S1000;
    }
    if(!(_S1008 != int(1)))
    {
        float3  _S1030 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S1030, _S1030)) != 0.0f)
        {
            normal_4 = _S1030 / make_float3 (length_2(_S1030));
        }
        else
        {
            normal_4 = _S1030;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, bool is_fisheye_7, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_10, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1031 = *dpdepths_1;
    float3  _S1032 = make_float3 (0.0f);
    bool _S1033 = !!_s_diff_ctx_1->_S990;
    float3  raydir_26;
    float3  raydir_27;
    float3  raydir_28;
    float3  raydir_29;
    int _S1034;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S1035;
    if(_S1033)
    {
        if(is_fisheye_7)
        {
            float _S1036 = length_1(_s_diff_ctx_1->_S989);
            float3  raydir_30 = make_float3 ((_s_diff_ctx_1->_S989 / make_float2 (s_primal_ctx_max_0(_S1036, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1036))).x, (_s_diff_ctx_1->_S989 / make_float2 (s_primal_ctx_max_0(_S1036, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1036))).y, s_primal_ctx_cos_0(_S1036));
            if(!is_ray_depth_3)
            {
                raydir_26 = raydir_30 / make_float3 (raydir_30.z);
            }
            else
            {
                raydir_26 = raydir_30;
            }
        }
        else
        {
            float3  raydir_31 = make_float3 (_s_diff_ctx_1->_S989.x, _s_diff_ctx_1->_S989.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_26 = normalize_0(raydir_31);
            }
            else
            {
                raydir_26 = raydir_31;
            }
        }
        float3  _S1037 = make_float3 (_S1031.primal_0.x) * raydir_26;
        bool _S1038 = !!_s_diff_ctx_1->_S992;
        if(_S1038)
        {
            if(is_fisheye_7)
            {
                float _S1039 = length_1(_s_diff_ctx_1->_S991);
                float3  raydir_32 = make_float3 ((_s_diff_ctx_1->_S991 / make_float2 (s_primal_ctx_max_0(_S1039, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1039))).x, (_s_diff_ctx_1->_S991 / make_float2 (s_primal_ctx_max_0(_S1039, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1039))).y, s_primal_ctx_cos_0(_S1039));
                if(!is_ray_depth_3)
                {
                    raydir_27 = raydir_32 / make_float3 (raydir_32.z);
                }
                else
                {
                    raydir_27 = raydir_32;
                }
            }
            else
            {
                float3  raydir_33 = make_float3 (_s_diff_ctx_1->_S991.x, _s_diff_ctx_1->_S991.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_27 = normalize_0(raydir_33);
                }
                else
                {
                    raydir_27 = raydir_33;
                }
            }
            float3  _S1040 = make_float3 (_S1031.primal_0.y) * raydir_27;
            _S1034 = int(2);
            points_4[int(0)] = _S1037;
            points_4[int(1)] = _S1040;
            points_4[int(2)] = _S1032;
            points_4[int(3)] = _S1032;
        }
        else
        {
            _S1034 = int(0);
            points_4[int(0)] = _S1037;
            points_4[int(1)] = _S1032;
            points_4[int(2)] = _S1032;
            points_4[int(3)] = _S1032;
            raydir_27 = _S1032;
        }
        if(_S1034 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S1033;
            _S1034 = int(0);
        }
        if(_runFlag_1)
        {
            if(!_s_diff_ctx_1->_S994)
            {
                _runFlag_2 = false;
                _S1034 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S1041 = length_1(_s_diff_ctx_1->_S993);
                    float3  raydir_34 = make_float3 ((_s_diff_ctx_1->_S993 / make_float2 (s_primal_ctx_max_0(_S1041, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1041))).x, (_s_diff_ctx_1->_S993 / make_float2 (s_primal_ctx_max_0(_S1041, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1041))).y, s_primal_ctx_cos_0(_S1041));
                    if(!is_ray_depth_3)
                    {
                        raydir_28 = raydir_34 / make_float3 (raydir_34.z);
                    }
                    else
                    {
                        raydir_28 = raydir_34;
                    }
                }
                else
                {
                    float3  raydir_35 = make_float3 (_s_diff_ctx_1->_S993.x, _s_diff_ctx_1->_S993.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_28 = normalize_0(raydir_35);
                    }
                    else
                    {
                        raydir_28 = raydir_35;
                    }
                }
                points_4[int(2)] = make_float3 (_S1031.primal_0.z) * raydir_28;
                bool _S1042 = !!_s_diff_ctx_1->_S996;
                int _S1043;
                if(_S1042)
                {
                    if(is_fisheye_7)
                    {
                        float _S1044 = length_1(_s_diff_ctx_1->_S995);
                        float3  raydir_36 = make_float3 ((_s_diff_ctx_1->_S995 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).x, (_s_diff_ctx_1->_S995 / make_float2 (s_primal_ctx_max_0(_S1044, 1.00000001168609742e-07f)) * make_float2 (s_primal_ctx_sin_0(_S1044))).y, s_primal_ctx_cos_0(_S1044));
                        if(!is_ray_depth_3)
                        {
                            raydir_29 = raydir_36 / make_float3 (raydir_36.z);
                        }
                        else
                        {
                            raydir_29 = raydir_36;
                        }
                    }
                    else
                    {
                        float3  raydir_37 = make_float3 (_s_diff_ctx_1->_S995.x, _s_diff_ctx_1->_S995.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_29 = normalize_0(raydir_37);
                        }
                        else
                        {
                            raydir_29 = raydir_37;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S1031.primal_0.w) * raydir_29;
                    _S1043 = int(2);
                }
                else
                {
                    _S1043 = int(0);
                    raydir_29 = _S1032;
                }
                if(_S1043 != int(2))
                {
                    _runFlag_3 = false;
                    _S1034 = _S1043;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S1034 = int(1);
                }
                float3  _S1045 = raydir_28;
                _runFlag_3 = _S1042;
                raydir_28 = raydir_29;
                raydir_29 = _S1045;
            }
            else
            {
                _runFlag_3 = false;
                raydir_28 = _S1032;
                raydir_29 = _S1032;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_28 = _S1032;
            raydir_29 = _S1032;
        }
        float3  _S1046 = raydir_26;
        float3  _S1047 = raydir_27;
        raydir_26 = raydir_28;
        raydir_27 = raydir_29;
        _S1035 = _S1038;
        raydir_28 = _S1047;
        raydir_29 = _S1046;
    }
    else
    {
        _S1034 = int(0);
        points_4[int(0)] = _S1032;
        points_4[int(1)] = _S1032;
        points_4[int(2)] = _S1032;
        points_4[int(3)] = _S1032;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_26 = _S1032;
        raydir_27 = _S1032;
        _S1035 = false;
        raydir_28 = _S1032;
        raydir_29 = _S1032;
    }
    bool _S1048 = !(_S1034 != int(1));
    float3  _S1049;
    float3  _S1050;
    float3  _S1051;
    float3  _S1052;
    float3  _S1053;
    bool _S1054;
    if(_S1048)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S1055 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S1056 = s_primal_ctx_cross_0(dx_1, _S1055);
        bool _S1057 = (s_primal_ctx_dot_0(_S1056, _S1056)) != 0.0f;
        if(_S1057)
        {
            float _S1058 = length_2(_S1056);
            float3  _S1059 = make_float3 (_S1058);
            _S1049 = make_float3 (_S1058 * _S1058);
            _S1050 = _S1059;
        }
        else
        {
            _S1049 = _S1032;
            _S1050 = _S1032;
        }
        float3  _S1060 = _S1050;
        _S1054 = _S1057;
        _S1050 = _S1056;
        _S1051 = _S1060;
        _S1052 = dx_1;
        _S1053 = _S1055;
    }
    else
    {
        _S1054 = false;
        _S1049 = _S1032;
        _S1050 = _S1032;
        _S1051 = _S1032;
        _S1052 = _S1032;
        _S1053 = _S1032;
    }
    float4  _S1061 = make_float4 (0.0f);
    if(_S1048)
    {
        if(_S1054)
        {
            float3  _S1062 = _s_dOut_10 / _S1049;
            float3  _S1063 = _S1050 * - _S1062;
            float3  _S1064 = _S1051 * _S1062;
            float _S1065 = _S1063.x + _S1063.y + _S1063.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1066;
            (&_S1066)->primal_0 = _S1050;
            (&_S1066)->differential_0 = _S1032;
            s_bwd_length_impl_1(&_S1066, _S1065);
            _S1049 = _S1064 + _S1066.differential_0;
        }
        else
        {
            _S1049 = _s_dOut_10;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1067;
        (&_S1067)->primal_0 = _S1050;
        (&_S1067)->differential_0 = _S1032;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1068;
        (&_S1068)->primal_0 = _S1050;
        (&_S1068)->differential_0 = _S1032;
        s_bwd_prop_dot_0(&_S1067, &_S1068, 0.0f);
        float3  _S1069 = _S1068.differential_0 + _S1067.differential_0 + _S1049;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1070;
        (&_S1070)->primal_0 = _S1052;
        (&_S1070)->differential_0 = _S1032;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1071;
        (&_S1071)->primal_0 = _S1053;
        (&_S1071)->differential_0 = _S1032;
        s_bwd_prop_cross_0(&_S1070, &_S1071, _S1069);
        float3  s_diff_dy_T_1 = - _S1071.differential_0;
        float3  _S1072 = - s_diff_dy_T_1;
        float3  _S1073 = - _S1070.differential_0;
        FixedArray<float3 , 4>  _S1074;
        _S1074[int(0)] = _S1032;
        _S1074[int(1)] = _S1032;
        _S1074[int(2)] = _S1032;
        _S1074[int(3)] = _S1032;
        _S1074[int(2)] = _S1072;
        _S1074[int(3)] = s_diff_dy_T_1;
        _S1074[int(0)] = _S1073;
        _S1074[int(1)] = _S1070.differential_0;
        points_4[int(0)] = _S1074[int(0)];
        points_4[int(1)] = _S1074[int(1)];
        points_4[int(2)] = _S1074[int(2)];
        points_4[int(3)] = _S1074[int(3)];
    }
    else
    {
        points_4[int(0)] = _S1032;
        points_4[int(1)] = _S1032;
        points_4[int(2)] = _S1032;
        points_4[int(3)] = _S1032;
    }
    float4  _S1075;
    if(_S1033)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S1076 = points_4;
                FixedArray<float3 , 4>  _S1077 = points_4;
                FixedArray<float3 , 4>  _S1078 = points_4;
                FixedArray<float3 , 4>  _S1079 = points_4;
                if(_runFlag_3)
                {
                    float3  _S1080 = raydir_26 * _S1079[int(3)];
                    float _S1081 = _S1080.x + _S1080.y + _S1080.z;
                    float4  _S1082 = _S1061;
                    *&((&_S1082)->w) = _S1081;
                    points_4[int(0)] = _S1076[int(0)];
                    points_4[int(1)] = _S1077[int(1)];
                    points_4[int(2)] = _S1078[int(2)];
                    points_4[int(3)] = _S1032;
                    _S1075 = _S1082;
                }
                else
                {
                    points_4[int(0)] = _S1076[int(0)];
                    points_4[int(1)] = _S1077[int(1)];
                    points_4[int(2)] = _S1078[int(2)];
                    points_4[int(3)] = _S1079[int(3)];
                    _S1075 = _S1061;
                }
                float3  _S1083 = raydir_27 * points_4[int(2)];
                float _S1084 = _S1083.x + _S1083.y + _S1083.z;
                FixedArray<float3 , 4>  _S1085 = points_4;
                FixedArray<float3 , 4>  _S1086 = points_4;
                float4  _S1087 = _S1061;
                *&((&_S1087)->z) = _S1084;
                float4  _S1088 = _S1075 + _S1087;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1085[int(1)];
                points_4[int(2)] = _S1032;
                points_4[int(3)] = _S1086[int(3)];
                _S1075 = _S1088;
            }
            else
            {
                FixedArray<float3 , 4>  _S1089 = points_4;
                FixedArray<float3 , 4>  _S1090 = points_4;
                FixedArray<float3 , 4>  _S1091 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S1089[int(1)];
                points_4[int(2)] = _S1090[int(2)];
                points_4[int(3)] = _S1091[int(3)];
                _S1075 = _S1061;
            }
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
            _S1075 = _S1061;
        }
        if(_S1035)
        {
            FixedArray<float3 , 4>  _S1095 = points_4;
            float3  _S1096 = raydir_28 * points_4[int(1)];
            float _S1097 = _S1096.x + _S1096.y + _S1096.z;
            float4  _S1098 = _S1061;
            *&((&_S1098)->y) = _S1097;
            float4  _S1099 = _S1075 + _S1098;
            points_4[int(0)] = _S1032;
            points_4[int(1)] = _S1032;
            points_4[int(2)] = _S1032;
            points_4[int(3)] = _S1032;
            raydir_26 = _S1095[int(0)];
            _S1075 = _S1099;
        }
        else
        {
            FixedArray<float3 , 4>  _S1100 = points_4;
            FixedArray<float3 , 4>  _S1101 = points_4;
            FixedArray<float3 , 4>  _S1102 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S1100[int(1)];
            points_4[int(2)] = _S1101[int(2)];
            points_4[int(3)] = _S1102[int(3)];
            raydir_26 = _S1032;
        }
        float3  _S1103 = raydir_29 * (points_4[int(0)] + raydir_26);
        float _S1104 = _S1103.x + _S1103.y + _S1103.z;
        float4  _S1105 = _S1061;
        *&((&_S1105)->x) = _S1104;
        _S1075 = _S1075 + _S1105;
    }
    else
    {
        _S1075 = _S1061;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S1075;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S1106, float4  _S1107, FixedArray<float, 10>  * _S1108, bool _S1109, bool _S1110, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1111, float3  _S1112)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1113;
    float3  _S1114 = s_primal_ctx_depth_to_normal_0(_S1106, _S1107, _S1108, _S1109, _S1110, (*_S1111).primal_0, &_S1113);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S1115 = _S1113;
    s_bwd_prop_depth_to_normal_0(_S1106, _S1107, _S1108, _S1109, _S1110, _S1111, _S1112, &_S1115);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_11, bool is_fisheye_8, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S1116 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S1116;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, dist_coeffs_11, is_fisheye_8, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_9)
{
    float2  _S1117 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_16 = _S1117;
    bool _S1118 = undistort_point_0(_S1117, dist_coeffs_12, int(12), &uv_16);
    if(!_S1118)
    {
        return 0.0f;
    }
    float3  raydir_38;
    if(is_fisheye_9)
    {
        float theta_9 = length_1(uv_16);
        float3  raydir_39 = make_float3 ((uv_16 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).x, (uv_16 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).y, (F32_cos((theta_9))));
        raydir_38 = raydir_39 / make_float3 (raydir_39.z);
    }
    else
    {
        raydir_38 = make_float3 (uv_16.x, uv_16.y, 1.0f);
    }
    return float((F32_sign((raydir_38.z)))) / length_2(raydir_38);
}

