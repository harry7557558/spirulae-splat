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

inline __device__ DiffPair_float_0 _d_exp_1(DiffPair_float_0 * dpx_1)
{
    float _S2 = (F32_exp((dpx_1->primal_0)));
    DiffPair_float_0 _S3 = { _S2, _S2 * dpx_1->differential_0 };
    return _S3;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S4 = *dpx_2;
    float _S5;
    if(((*dpx_2).primal_0) > ((*dpy_0).primal_0))
    {
        _S5 = dOut_1;
    }
    else
    {
        if(((*dpx_2).primal_0) < ((*dpy_0).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_1;
        }
    }
    dpx_2->primal_0 = _S4.primal_0;
    dpx_2->differential_0 = _S5;
    DiffPair_float_0 _S6 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S4.primal_0))
    {
        _S5 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_2).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S6.primal_0;
    dpy_0->differential_0 = _S5;
    return;
}

inline __device__ DiffPair_float_0 _d_max_1(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_1)
{
    float _S7 = dpx_3->primal_0;
    float _S8 = dpy_1->primal_0;
    float _S9 = (F32_max((dpx_3->primal_0), (dpy_1->primal_0)));
    float _S10;
    if((dpx_3->primal_0) > (dpy_1->primal_0))
    {
        _S10 = dpx_3->differential_0;
    }
    else
    {
        if(_S7 < _S8)
        {
            _S10 = dpy_1->differential_0;
        }
        else
        {
            _S10 = 0.5f * (dpx_3->differential_0 + dpy_1->differential_0);
        }
    }
    DiffPair_float_0 _S11 = { _S9, _S10 };
    return _S11;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_4, float dOut_2)
{
    float _S12 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_4).primal_0)))))) * dOut_2;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S12;
    return;
}

inline __device__ DiffPair_float_0 _d_sqrt_1(DiffPair_float_0 * dpx_5)
{
    DiffPair_float_0 _S13 = { (F32_sqrt((dpx_5->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_5->primal_0)))))) * dpx_5->differential_0 };
    return _S13;
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

inline __device__ float length_0(float4  x_1)
{
    return (F32_sqrt((dot_0(x_1, x_1))));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_6, float dOut_3)
{
    float _S14 = 1.0f / (*dpx_6).primal_0 * dOut_3;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S14;
    return;
}

inline __device__ DiffPair_float_0 _d_log_1(DiffPair_float_0 * dpx_7)
{
    DiffPair_float_0 _S15 = { (F32_log((dpx_7->primal_0))), 1.0f / dpx_7->primal_0 * dpx_7->differential_0 };
    return _S15;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ float3  exp_0(float3  x_2)
{
    float3  result_2;
    int i_1 = int(0);
    for(;;)
    {
        if(i_1 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_2, i_1) = (F32_exp((_slang_vector_get_element(x_2, i_1))));
        i_1 = i_1 + int(1);
    }
    return result_2;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  dOut_4)
{
    float3  _S16 = exp_0((*dpx_8).primal_0) * dOut_4;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S16;
    return;
}

inline __device__ DiffPair_vectorx3Cfloatx2C3x3E_0 _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9)
{
    float3  _S17 = exp_0(dpx_9->primal_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S18 = { _S17, _S17 * dpx_9->differential_0 };
    return _S18;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_10, DiffPair_float_0 * dpy_2, float dOut_5)
{
    DiffPair_float_0 _S19 = *dpx_10;
    float _S20;
    if(((*dpx_10).primal_0) < ((*dpy_2).primal_0))
    {
        _S20 = dOut_5;
    }
    else
    {
        if(((*dpx_10).primal_0) > ((*dpy_2).primal_0))
        {
            _S20 = 0.0f;
        }
        else
        {
            _S20 = 0.5f * dOut_5;
        }
    }
    dpx_10->primal_0 = _S19.primal_0;
    dpx_10->differential_0 = _S20;
    DiffPair_float_0 _S21 = *dpy_2;
    if(((*dpy_2).primal_0) < (_S19.primal_0))
    {
        _S20 = dOut_5;
    }
    else
    {
        if(((*dpy_2).primal_0) > ((*dpx_10).primal_0))
        {
            _S20 = 0.0f;
        }
        else
        {
            _S20 = 0.5f * dOut_5;
        }
    }
    dpy_2->primal_0 = _S21.primal_0;
    dpy_2->differential_0 = _S20;
    return;
}

inline __device__ DiffPair_float_0 _d_min_1(DiffPair_float_0 * dpx_11, DiffPair_float_0 * dpy_3)
{
    float _S22 = dpx_11->primal_0;
    float _S23 = dpy_3->primal_0;
    float _S24 = (F32_min((dpx_11->primal_0), (dpy_3->primal_0)));
    float _S25;
    if((dpx_11->primal_0) < (dpy_3->primal_0))
    {
        _S25 = dpx_11->differential_0;
    }
    else
    {
        if(_S22 > _S23)
        {
            _S25 = dpy_3->differential_0;
        }
        else
        {
            _S25 = 0.5f * (dpx_11->differential_0 + dpy_3->differential_0);
        }
    }
    DiffPair_float_0 _S26 = { _S24, _S25 };
    return _S26;
}

inline __device__ void per_splat_losses(float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S27)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    float3  _S28 = exp_0(scales_0);
    float _S29 = _S28.x;
    float _S30 = _S28.y;
    float _S31 = _S28.z;
    losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S29 + _S30 + _S31) / 3.0f;
    losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S29), (_S30)))), (_S31))) / (F32_min(((F32_min((_S29), (_S30)))), (_S31)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
    float3  _S32 = exp_0(make_float3 (2.0f) * scales_0);
    float x_3 = _S32.x;
    float y_1 = _S32.y;
    float z_0 = _S32.z;
    float s_0 = x_3 + y_1 + z_0;
    float s1_0 = (F32_max(((F32_max((x_3), (y_1)))), (z_0))) / s_0;
    float s3_0 = (F32_min(((F32_min((x_3), (y_1)))), (z_0))) / s_0;
    float s2_0 = 1.0f - s1_0 - s3_0;
    losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    *_S27 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_exp_0(float _S33)
{
    return (F32_exp((_S33)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S34)
{
    return exp_0(_S34);
}

inline __device__ float s_primal_ctx_log_0(float _S35)
{
    return (F32_log((_S35)));
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S36, float _S37)
{
    _d_log_0(_S36, _S37);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S38, float _S39)
{
    _d_exp_0(_S38, _S39);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S40, float3  _S41)
{
    _d_exp_vector_0(_S40, _S41);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S42, float _S43)
{
    _d_sqrt_0(_S42, _S43);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_12, float _s_dOut_0)
{
    float _S44 = (*dpx_12).primal_0.x;
    float _S45 = (*dpx_12).primal_0.y;
    float _S46 = (*dpx_12).primal_0.z;
    float _S47 = (*dpx_12).primal_0.w;
    DiffPair_float_0 _S48;
    (&_S48)->primal_0 = _S44 * _S44 + _S45 * _S45 + _S46 * _S46 + _S47 * _S47;
    (&_S48)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S48, _s_dOut_0);
    float _S49 = (*dpx_12).primal_0.w * _S48.differential_0;
    float _S50 = _S49 + _S49;
    float _S51 = (*dpx_12).primal_0.z * _S48.differential_0;
    float _S52 = _S51 + _S51;
    float _S53 = (*dpx_12).primal_0.y * _S48.differential_0;
    float _S54 = _S53 + _S53;
    float _S55 = (*dpx_12).primal_0.x * _S48.differential_0;
    float _S56 = _S55 + _S55;
    float4  _S57 = make_float4 (0.0f);
    *&((&_S57)->w) = _S50;
    *&((&_S57)->z) = _S52;
    *&((&_S57)->y) = _S54;
    *&((&_S57)->x) = _S56;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S57;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S58, float _S59)
{
    s_bwd_prop_length_impl_0(_S58, _S59);
    return;
}

inline __device__ void per_splat_losses_bwd(float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1)
{
    float _S60 = - opacity_1;
    float _S61 = 1.0f + s_primal_ctx_exp_0(_S60);
    float _S62 = _S61 * _S61;
    float _S63 = length_0(quat_1);
    float3  _S64 = s_primal_ctx_exp_1(scales_1);
    float _S65 = _S64.x;
    float _S66 = _S64.y;
    float _S67 = _S64.z;
    float _S68 = (F32_max((_S65), (_S66)));
    float _S69 = (F32_max((_S68), (_S67)));
    float _S70 = (F32_min((_S65), (_S66)));
    float _S71 = (F32_min((_S70), (_S67)));
    float _S72 = _S69 / _S71;
    float _S73 = _S71 * _S71;
    float3  _S74 = make_float3 (2.0f) * scales_1;
    float3  _S75 = s_primal_ctx_exp_1(_S74);
    float x_4 = _S75.x;
    float y_2 = _S75.y;
    float z_1 = _S75.z;
    float s_1 = x_4 + y_2 + z_1;
    float _S76 = (F32_max((x_4), (y_2)));
    float _S77 = (F32_max((_S76), (z_1)));
    float s1_1 = _S77 / s_1;
    float _S78 = s_1 * s_1;
    float _S79 = (F32_min((x_4), (y_2)));
    float _S80 = (F32_min((_S79), (z_1)));
    float s3_1 = _S80 / s_1;
    float s2_1 = 1.0f - s1_1 - s3_1;
    float _S81 = - s1_1;
    float _S82 = s_primal_ctx_log_0(s1_1);
    float _S83 = s_primal_ctx_log_0(s2_1);
    float _S84 = s_primal_ctx_log_0(s3_1);
    float _S85 = _S81 * _S82 - s2_1 * _S83 - s3_1 * _S84;
    float _S86 = s_primal_ctx_exp_0(_S85) - 0.99998998641967773f;
    float _S87 = erank_reg_weight_s3_1 * v_loss_0[int(3)];
    float _S88 = erank_reg_weight_1 * v_loss_0[int(3)];
    DiffPair_float_0 _S89;
    (&_S89)->primal_0 = - s_primal_ctx_log_0(_S86);
    (&_S89)->differential_0 = 0.0f;
    DiffPair_float_0 _S90;
    (&_S90)->primal_0 = 0.0f;
    (&_S90)->differential_0 = 0.0f;
    _d_max_0(&_S89, &_S90, _S88);
    float _S91 = - _S89.differential_0;
    DiffPair_float_0 _S92;
    (&_S92)->primal_0 = _S86;
    (&_S92)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S92, _S91);
    DiffPair_float_0 _S93;
    (&_S93)->primal_0 = _S85;
    (&_S93)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S93, _S92.differential_0);
    float _S94 = - _S93.differential_0;
    float _S95 = s3_1 * _S94;
    float _S96 = _S84 * _S94;
    DiffPair_float_0 _S97;
    (&_S97)->primal_0 = s3_1;
    (&_S97)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S97, _S95);
    float _S98 = s2_1 * _S94;
    float _S99 = _S83 * _S94;
    DiffPair_float_0 _S100;
    (&_S100)->primal_0 = s2_1;
    (&_S100)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S100, _S98);
    float _S101 = _S81 * _S93.differential_0;
    float _S102 = _S82 * _S93.differential_0;
    DiffPair_float_0 _S103;
    (&_S103)->primal_0 = s1_1;
    (&_S103)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S103, _S101);
    float _S104 = - _S102;
    float _S105 = - (_S99 + _S100.differential_0);
    float _S106 = (_S87 + _S96 + _S97.differential_0 + _S105) / _S78;
    float _S107 = _S80 * - _S106;
    float _S108 = s_1 * _S106;
    DiffPair_float_0 _S109;
    (&_S109)->primal_0 = _S79;
    (&_S109)->differential_0 = 0.0f;
    DiffPair_float_0 _S110;
    (&_S110)->primal_0 = z_1;
    (&_S110)->differential_0 = 0.0f;
    _d_min_0(&_S109, &_S110, _S108);
    DiffPair_float_0 _S111;
    (&_S111)->primal_0 = x_4;
    (&_S111)->differential_0 = 0.0f;
    DiffPair_float_0 _S112;
    (&_S112)->primal_0 = y_2;
    (&_S112)->differential_0 = 0.0f;
    _d_min_0(&_S111, &_S112, _S109.differential_0);
    float _S113 = (_S103.differential_0 + _S104 + _S105) / _S78;
    float _S114 = _S77 * - _S113;
    float _S115 = s_1 * _S113;
    DiffPair_float_0 _S116;
    (&_S116)->primal_0 = _S76;
    (&_S116)->differential_0 = 0.0f;
    DiffPair_float_0 _S117;
    (&_S117)->primal_0 = z_1;
    (&_S117)->differential_0 = 0.0f;
    _d_max_0(&_S116, &_S117, _S115);
    DiffPair_float_0 _S118;
    (&_S118)->primal_0 = x_4;
    (&_S118)->differential_0 = 0.0f;
    DiffPair_float_0 _S119;
    (&_S119)->primal_0 = y_2;
    (&_S119)->differential_0 = 0.0f;
    _d_max_0(&_S118, &_S119, _S116.differential_0);
    float _S120 = _S107 + _S114;
    float3  _S121 = make_float3 (_S111.differential_0 + _S118.differential_0 + _S120, _S112.differential_0 + _S119.differential_0 + _S120, _S110.differential_0 + _S117.differential_0 + _S120);
    float3  _S122 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S123;
    (&_S123)->primal_0 = _S74;
    (&_S123)->differential_0 = _S122;
    s_bwd_prop_exp_1(&_S123, _S121);
    float3  _S124 = make_float3 (2.0f) * _S123.differential_0;
    float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * v_loss_0[int(2)];
    DiffPair_float_0 _S125;
    (&_S125)->primal_0 = _S72;
    (&_S125)->differential_0 = 0.0f;
    DiffPair_float_0 _S126;
    (&_S126)->primal_0 = max_gauss_ratio_1;
    (&_S126)->differential_0 = 0.0f;
    _d_max_0(&_S125, &_S126, s_diff_scale_reg_T_0);
    float _S127 = _S125.differential_0 / _S73;
    float _S128 = _S69 * - _S127;
    float _S129 = _S71 * _S127;
    DiffPair_float_0 _S130;
    (&_S130)->primal_0 = _S70;
    (&_S130)->differential_0 = 0.0f;
    DiffPair_float_0 _S131;
    (&_S131)->primal_0 = _S67;
    (&_S131)->differential_0 = 0.0f;
    _d_min_0(&_S130, &_S131, _S128);
    DiffPair_float_0 _S132;
    (&_S132)->primal_0 = _S65;
    (&_S132)->differential_0 = 0.0f;
    DiffPair_float_0 _S133;
    (&_S133)->primal_0 = _S66;
    (&_S133)->differential_0 = 0.0f;
    _d_min_0(&_S132, &_S133, _S130.differential_0);
    DiffPair_float_0 _S134;
    (&_S134)->primal_0 = _S68;
    (&_S134)->differential_0 = 0.0f;
    DiffPair_float_0 _S135;
    (&_S135)->primal_0 = _S67;
    (&_S135)->differential_0 = 0.0f;
    _d_max_0(&_S134, &_S135, _S129);
    DiffPair_float_0 _S136;
    (&_S136)->primal_0 = _S65;
    (&_S136)->differential_0 = 0.0f;
    DiffPair_float_0 _S137;
    (&_S137)->primal_0 = _S66;
    (&_S137)->differential_0 = 0.0f;
    _d_max_0(&_S136, &_S137, _S134.differential_0);
    float _S138 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * v_loss_0[int(1)]);
    float3  _S139 = make_float3 (_S132.differential_0 + _S136.differential_0 + _S138, _S133.differential_0 + _S137.differential_0 + _S138, _S131.differential_0 + _S135.differential_0 + _S138);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S140;
    (&_S140)->primal_0 = scales_1;
    (&_S140)->differential_0 = _S122;
    s_bwd_prop_exp_1(&_S140, _S139);
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * v_loss_0[int(4)];
    float _S141 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S142;
    (&_S142)->primal_0 = _S63;
    (&_S142)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S142, _S141);
    float _S143 = _S142.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S144 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S145;
    (&_S145)->primal_0 = quat_1;
    (&_S145)->differential_0 = _S144;
    s_bwd_length_impl_0(&_S145, _S143);
    float _S146 = - (mcmc_opacity_reg_weight_1 * v_loss_0[int(0)] / _S62);
    DiffPair_float_0 _S147;
    (&_S147)->primal_0 = _S60;
    (&_S147)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S147, _S146);
    float _S148 = - _S147.differential_0;
    *v_scales_0 = _S124 + _S140.differential_0;
    *v_opacity_0 = _S148;
    *v_quat_0 = _S145.differential_0;
    return;
}

inline __device__ DiffPair_float_0 s_fwd_s_primal_ctx_exp_0(DiffPair_float_0 * _S149)
{
    DiffPair_float_0 _S150;
    (&_S150)->primal_0 = _S149->primal_0;
    (&_S150)->differential_0 = _S149->differential_0;
    DiffPair_float_0 _S151 = _d_exp_1(&_S150);
    DiffPair_float_0 _S152 = { _S151.primal_0, _S151.differential_0 };
    return _S152;
}

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_13)
{
    float _S153 = *&((&dpx_13->differential_0)->x) * *&((&dpx_13->primal_0)->x);
    float _S154 = *&((&dpx_13->differential_0)->y) * *&((&dpx_13->primal_0)->y);
    float _S155 = *&((&dpx_13->differential_0)->z) * *&((&dpx_13->primal_0)->z);
    float _S156 = *&((&dpx_13->differential_0)->w) * *&((&dpx_13->primal_0)->w);
    float s_diff_len_0 = _S153 + _S153 + (_S154 + _S154) + (_S155 + _S155) + (_S156 + _S156);
    DiffPair_float_0 _S157;
    (&_S157)->primal_0 = *&((&dpx_13->primal_0)->x) * *&((&dpx_13->primal_0)->x) + *&((&dpx_13->primal_0)->y) * *&((&dpx_13->primal_0)->y) + *&((&dpx_13->primal_0)->z) * *&((&dpx_13->primal_0)->z) + *&((&dpx_13->primal_0)->w) * *&((&dpx_13->primal_0)->w);
    (&_S157)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S158 = _d_sqrt_1(&_S157);
    DiffPair_float_0 _S159 = { _S158.primal_0, _S158.differential_0 };
    return _S159;
}

inline __device__ DiffPair_vectorx3Cfloatx2C3x3E_0 s_fwd_s_primal_ctx_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S160)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S161;
    (&_S161)->primal_0 = _S160->primal_0;
    (&_S161)->differential_0 = _S160->differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S162 = _d_exp_vector_1(&_S161);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S163 = { _S162.primal_0, _S162.differential_0 };
    return _S163;
}

inline __device__ DiffPair_float_0 s_fwd_s_primal_ctx_log_0(DiffPair_float_0 * _S164)
{
    DiffPair_float_0 _S165;
    (&_S165)->primal_0 = _S164->primal_0;
    (&_S165)->differential_0 = _S164->differential_0;
    DiffPair_float_0 _S166 = _d_log_1(&_S165);
    DiffPair_float_0 _S167 = { _S166.primal_0, _S166.differential_0 };
    return _S167;
}

struct DiffPair_0
{
    DiffPair_float_0 primal_0;
    DiffPair_float_0 differential_0;
};

inline __device__ void s_fwd_d_max_0(DiffPair_0 * dpdpx_0, DiffPair_0 * dpdpy_0, DiffPair_float_0 * dpdOut_0)
{
    DiffPair_0 _S168 = *dpdpx_0;
    DiffPair_0 _S169 = *dpdpy_0;
    float _S170 = dpdOut_0->differential_0;
    float _S171 = dpdOut_0->primal_0;
    float _S172;
    float _S173;
    if(((*dpdpx_0).primal_0.primal_0) > ((*dpdpy_0).primal_0.primal_0))
    {
        _S172 = _S171;
        _S173 = _S170;
    }
    else
    {
        if((_S168.primal_0.primal_0) < (_S169.primal_0.primal_0))
        {
            _S172 = 0.0f;
            _S173 = 0.0f;
        }
        else
        {
            float _S174 = _S170 * 0.5f;
            _S172 = 0.5f * _S171;
            _S173 = _S174;
        }
    }
    DiffPair_float_0 _S175 = { _S168.primal_0.primal_0, _S172 };
    DiffPair_float_0 _S176 = { _S168.differential_0.primal_0, _S173 };
    if((_S169.primal_0.primal_0) > (_S168.primal_0.primal_0))
    {
        _S172 = _S171;
        _S173 = _S170;
    }
    else
    {
        if((_S169.primal_0.primal_0) < (_S168.primal_0.primal_0))
        {
            _S172 = 0.0f;
            _S173 = 0.0f;
        }
        else
        {
            float _S177 = _S170 * 0.5f;
            _S172 = 0.5f * _S171;
            _S173 = _S177;
        }
    }
    DiffPair_float_0 _S178 = { _S169.primal_0.primal_0, _S172 };
    DiffPair_float_0 _S179 = { _S169.differential_0.primal_0, _S173 };
    dpdpx_0->primal_0 = _S175;
    dpdpx_0->differential_0 = _S176;
    dpdpy_0->primal_0 = _S178;
    dpdpy_0->differential_0 = _S179;
    return;
}

inline __device__ void s_fwd_d_log_0(DiffPair_0 * dpdpx_1, DiffPair_float_0 * dpdOut_1)
{
    float _S180 = 1.0f / (*dpdpx_1).primal_0.primal_0;
    DiffPair_float_0 _S181 = { (*dpdpx_1).primal_0.primal_0, _S180 * dpdOut_1->primal_0 };
    DiffPair_float_0 _S182 = { (*dpdpx_1).differential_0.primal_0, (0.0f - (*dpdpx_1).differential_0.primal_0) / ((*dpdpx_1).primal_0.primal_0 * (*dpdpx_1).primal_0.primal_0) * dpdOut_1->primal_0 + dpdOut_1->differential_0 * _S180 };
    dpdpx_1->primal_0 = _S181;
    dpdpx_1->differential_0 = _S182;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_log_0(DiffPair_0 * _S183, DiffPair_float_0 * _S184)
{
    DiffPair_0 _S185;
    (&_S185)->primal_0 = (*_S183).primal_0;
    (&_S185)->differential_0 = (*_S183).differential_0;
    DiffPair_float_0 _S186;
    (&_S186)->primal_0 = _S184->primal_0;
    (&_S186)->differential_0 = _S184->differential_0;
    s_fwd_d_log_0(&_S185, &_S186);
    _S183->primal_0 = _S185.primal_0;
    _S183->differential_0 = _S185.differential_0;
    return;
}

inline __device__ void s_fwd_d_exp_0(DiffPair_0 * dpdpx_2, DiffPair_float_0 * dpdOut_2)
{
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = (*dpdpx_2).primal_0.primal_0;
    (&_S187)->differential_0 = (*dpdpx_2).differential_0.primal_0;
    DiffPair_float_0 _S188 = _d_exp_1(&_S187);
    DiffPair_float_0 _S189 = { (*dpdpx_2).primal_0.primal_0, _S188.primal_0 * dpdOut_2->primal_0 };
    DiffPair_float_0 _S190 = { (*dpdpx_2).differential_0.primal_0, _S188.differential_0 * dpdOut_2->primal_0 + dpdOut_2->differential_0 * _S188.primal_0 };
    dpdpx_2->primal_0 = _S189;
    dpdpx_2->differential_0 = _S190;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_exp_0(DiffPair_0 * _S191, DiffPair_float_0 * _S192)
{
    DiffPair_0 _S193;
    (&_S193)->primal_0 = (*_S191).primal_0;
    (&_S193)->differential_0 = (*_S191).differential_0;
    DiffPair_float_0 _S194;
    (&_S194)->primal_0 = _S192->primal_0;
    (&_S194)->differential_0 = _S192->differential_0;
    s_fwd_d_exp_0(&_S193, &_S194);
    _S191->primal_0 = _S193.primal_0;
    _S191->differential_0 = _S193.differential_0;
    return;
}

inline __device__ void s_fwd_d_min_0(DiffPair_0 * dpdpx_3, DiffPair_0 * dpdpy_1, DiffPair_float_0 * dpdOut_3)
{
    DiffPair_0 _S195 = *dpdpx_3;
    DiffPair_0 _S196 = *dpdpy_1;
    float _S197 = dpdOut_3->differential_0;
    float _S198 = dpdOut_3->primal_0;
    float _S199;
    float _S200;
    if(((*dpdpx_3).primal_0.primal_0) < ((*dpdpy_1).primal_0.primal_0))
    {
        _S199 = _S198;
        _S200 = _S197;
    }
    else
    {
        if((_S195.primal_0.primal_0) > (_S196.primal_0.primal_0))
        {
            _S199 = 0.0f;
            _S200 = 0.0f;
        }
        else
        {
            float _S201 = _S197 * 0.5f;
            _S199 = 0.5f * _S198;
            _S200 = _S201;
        }
    }
    DiffPair_float_0 _S202 = { _S195.primal_0.primal_0, _S199 };
    DiffPair_float_0 _S203 = { _S195.differential_0.primal_0, _S200 };
    if((_S196.primal_0.primal_0) < (_S195.primal_0.primal_0))
    {
        _S199 = _S198;
        _S200 = _S197;
    }
    else
    {
        if((_S196.primal_0.primal_0) > (_S195.primal_0.primal_0))
        {
            _S199 = 0.0f;
            _S200 = 0.0f;
        }
        else
        {
            float _S204 = _S197 * 0.5f;
            _S199 = 0.5f * _S198;
            _S200 = _S204;
        }
    }
    DiffPair_float_0 _S205 = { _S196.primal_0.primal_0, _S199 };
    DiffPair_float_0 _S206 = { _S196.differential_0.primal_0, _S200 };
    dpdpx_3->primal_0 = _S202;
    dpdpx_3->differential_0 = _S203;
    dpdpy_1->primal_0 = _S205;
    dpdpy_1->differential_0 = _S206;
    return;
}

struct DiffPair_1
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 differential_0;
};

inline __device__ void s_fwd_d_exp_vector_0(DiffPair_1 * dpdpx_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdOut_4)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S207;
    (&_S207)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S207)->differential_0 = (*dpdpx_4).differential_0.primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S208 = _d_exp_vector_1(&_S207);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S209 = { (*dpdpx_4).primal_0.primal_0, _S208.primal_0 * dpdOut_4->primal_0 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S210 = { (*dpdpx_4).differential_0.primal_0, _S208.differential_0 * dpdOut_4->primal_0 + dpdOut_4->differential_0 * _S208.primal_0 };
    dpdpx_4->primal_0 = _S209;
    dpdpx_4->differential_0 = _S210;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_exp_1(DiffPair_1 * _S211, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S212)
{
    DiffPair_1 _S213;
    (&_S213)->primal_0 = (*_S211).primal_0;
    (&_S213)->differential_0 = (*_S211).differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S214;
    (&_S214)->primal_0 = _S212->primal_0;
    (&_S214)->differential_0 = _S212->differential_0;
    s_fwd_d_exp_vector_0(&_S213, &_S214);
    _S211->primal_0 = _S213.primal_0;
    _S211->differential_0 = _S213.differential_0;
    return;
}

struct DiffPair_2
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 differential_0;
};

inline __device__ void s_fwd_d_sqrt_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dpdOut_5)
{
    DiffPair_float_0 _S215;
    (&_S215)->primal_0 = 1.00000001168609742e-07f;
    (&_S215)->differential_0 = 0.0f;
    DiffPair_float_0 _S216;
    (&_S216)->primal_0 = (*dpdpx_5).primal_0.primal_0;
    (&_S216)->differential_0 = (*dpdpx_5).differential_0.primal_0;
    DiffPair_float_0 _S217 = _d_max_1(&_S215, &_S216);
    DiffPair_float_0 _S218;
    (&_S218)->primal_0 = _S217.primal_0;
    (&_S218)->differential_0 = _S217.differential_0;
    DiffPair_float_0 _S219 = _d_sqrt_1(&_S218);
    float _S220 = 0.5f / _S219.primal_0;
    DiffPair_float_0 _S221 = { (*dpdpx_5).primal_0.primal_0, _S220 * dpdOut_5->primal_0 };
    DiffPair_float_0 _S222 = { (*dpdpx_5).differential_0.primal_0, (0.0f - 0.5f * _S219.differential_0) / (_S219.primal_0 * _S219.primal_0) * dpdOut_5->primal_0 + dpdOut_5->differential_0 * _S220 };
    dpdpx_5->primal_0 = _S221;
    dpdpx_5->differential_0 = _S222;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_sqrt_0(DiffPair_0 * _S223, DiffPair_float_0 * _S224)
{
    DiffPair_0 _S225;
    (&_S225)->primal_0 = (*_S223).primal_0;
    (&_S225)->differential_0 = (*_S223).differential_0;
    DiffPair_float_0 _S226;
    (&_S226)->primal_0 = _S224->primal_0;
    (&_S226)->differential_0 = _S224->differential_0;
    s_fwd_d_sqrt_0(&_S225, &_S226);
    _S223->primal_0 = _S225.primal_0;
    _S223->differential_0 = _S225.differential_0;
    return;
}

inline __device__ void s_fwd_s_bwd_prop_length_impl_0(DiffPair_2 * dpdpx_6, DiffPair_float_0 * dp_s_dOut_0)
{
    float _S227 = (*dpdpx_6).primal_0.primal_0.x;
    float _S228 = (*dpdpx_6).differential_0.primal_0.x * (*dpdpx_6).primal_0.primal_0.x;
    float _S229 = (*dpdpx_6).primal_0.primal_0.y;
    float _S230 = (*dpdpx_6).differential_0.primal_0.y * (*dpdpx_6).primal_0.primal_0.y;
    float _S231 = (*dpdpx_6).primal_0.primal_0.z;
    float _S232 = (*dpdpx_6).differential_0.primal_0.z * (*dpdpx_6).primal_0.primal_0.z;
    float _S233 = (*dpdpx_6).primal_0.primal_0.w;
    float _S234 = (*dpdpx_6).differential_0.primal_0.w * (*dpdpx_6).primal_0.primal_0.w;
    DiffPair_float_0 _S235 = { _S227 * _S227 + _S229 * _S229 + _S231 * _S231 + _S233 * _S233, 0.0f };
    DiffPair_float_0 _S236 = { _S228 + _S228 + (_S230 + _S230) + (_S232 + _S232) + (_S234 + _S234), 0.0f };
    DiffPair_0 _S237;
    (&_S237)->primal_0 = _S235;
    (&_S237)->differential_0 = _S236;
    DiffPair_float_0 _S238;
    (&_S238)->primal_0 = dp_s_dOut_0->primal_0;
    (&_S238)->differential_0 = dp_s_dOut_0->differential_0;
    s_fwd_s_bwd_prop_sqrt_0(&_S237, &_S238);
    float _S239 = (*dpdpx_6).primal_0.primal_0.w * _S237.primal_0.differential_0;
    float _S240 = (*dpdpx_6).differential_0.primal_0.w * _S237.primal_0.differential_0 + _S237.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.w;
    float _S241 = _S239 + _S239;
    float _S242 = _S240 + _S240;
    float _S243 = (*dpdpx_6).primal_0.primal_0.z * _S237.primal_0.differential_0;
    float _S244 = (*dpdpx_6).differential_0.primal_0.z * _S237.primal_0.differential_0 + _S237.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.z;
    float _S245 = _S243 + _S243;
    float _S246 = _S244 + _S244;
    float _S247 = (*dpdpx_6).primal_0.primal_0.y * _S237.primal_0.differential_0;
    float _S248 = (*dpdpx_6).differential_0.primal_0.y * _S237.primal_0.differential_0 + _S237.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.y;
    float _S249 = _S247 + _S247;
    float _S250 = _S248 + _S248;
    float _S251 = (*dpdpx_6).primal_0.primal_0.x * _S237.primal_0.differential_0;
    float _S252 = (*dpdpx_6).differential_0.primal_0.x * _S237.primal_0.differential_0 + _S237.differential_0.differential_0 * (*dpdpx_6).primal_0.primal_0.x;
    float _S253 = _S251 + _S251;
    float _S254 = _S252 + _S252;
    float4  _S255 = make_float4 (0.0f);
    float4  _S256 = _S255;
    *&((&_S256)->w) = _S241;
    float4  _S257 = _S255;
    *&((&_S257)->w) = _S242;
    *&((&_S256)->z) = _S245;
    *&((&_S257)->z) = _S246;
    *&((&_S256)->y) = _S249;
    *&((&_S257)->y) = _S250;
    *&((&_S256)->x) = _S253;
    *&((&_S257)->x) = _S254;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S258 = { (*dpdpx_6).primal_0.primal_0, _S256 };
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S259 = { (*dpdpx_6).differential_0.primal_0, _S257 };
    dpdpx_6->primal_0 = _S258;
    dpdpx_6->differential_0 = _S259;
    return;
}

inline __device__ void s_fwd_s_bwd_length_impl_0(DiffPair_2 * _S260, DiffPair_float_0 * _S261)
{
    DiffPair_2 _S262;
    (&_S262)->primal_0 = (*_S260).primal_0;
    (&_S262)->differential_0 = (*_S260).differential_0;
    DiffPair_float_0 _S263;
    (&_S263)->primal_0 = _S261->primal_0;
    (&_S263)->differential_0 = _S261->differential_0;
    s_fwd_s_bwd_prop_length_impl_0(&_S262, &_S263);
    _S260->primal_0 = _S262.primal_0;
    _S260->differential_0 = _S262.differential_0;
    return;
}

inline __device__ void s_fwd_per_splat_losses_bwd_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, FixedArray<float, 5>  * v_loss_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpv_scales_0, DiffPair_float_0 * dpv_opacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpv_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float _S264 = - dpopacity_0->primal_0;
    float _S265 = - dpopacity_0->differential_0;
    DiffPair_float_0 _S266;
    (&_S266)->primal_0 = _S264;
    (&_S266)->differential_0 = _S265;
    DiffPair_float_0 _S267 = s_fwd_s_primal_ctx_exp_0(&_S266);
    float _S268 = 1.0f + _S267.primal_0;
    float _S269 = _S268 * _S268;
    float _S270 = _S267.differential_0 * _S268;
    float _S271 = _S270 + _S270;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S272;
    (&_S272)->primal_0 = dpquat_0->primal_0;
    (&_S272)->differential_0 = dpquat_0->differential_0;
    DiffPair_float_0 _S273 = s_fwd_length_impl_0(&_S272);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S274;
    (&_S274)->primal_0 = dpscales_0->primal_0;
    (&_S274)->differential_0 = dpscales_0->differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S275 = s_fwd_s_primal_ctx_exp_1(&_S274);
    float _S276 = _S275.primal_0.x;
    float _S277 = _S275.differential_0.x;
    float _S278 = _S275.primal_0.y;
    float _S279 = _S275.differential_0.y;
    float _S280 = _S275.primal_0.z;
    float _S281 = _S275.differential_0.z;
    DiffPair_float_0 _S282;
    (&_S282)->primal_0 = _S276;
    (&_S282)->differential_0 = _S277;
    DiffPair_float_0 _S283;
    (&_S283)->primal_0 = _S278;
    (&_S283)->differential_0 = _S279;
    DiffPair_float_0 _S284 = _d_max_1(&_S282, &_S283);
    DiffPair_float_0 _S285;
    (&_S285)->primal_0 = _S284.primal_0;
    (&_S285)->differential_0 = _S284.differential_0;
    DiffPair_float_0 _S286;
    (&_S286)->primal_0 = _S280;
    (&_S286)->differential_0 = _S281;
    DiffPair_float_0 _S287 = _d_max_1(&_S285, &_S286);
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = _S276;
    (&_S288)->differential_0 = _S277;
    DiffPair_float_0 _S289;
    (&_S289)->primal_0 = _S278;
    (&_S289)->differential_0 = _S279;
    DiffPair_float_0 _S290 = _d_min_1(&_S288, &_S289);
    DiffPair_float_0 _S291;
    (&_S291)->primal_0 = _S290.primal_0;
    (&_S291)->differential_0 = _S290.differential_0;
    DiffPair_float_0 _S292;
    (&_S292)->primal_0 = _S280;
    (&_S292)->differential_0 = _S281;
    DiffPair_float_0 _S293 = _d_min_1(&_S291, &_S292);
    float _S294 = _S287.primal_0 / _S293.primal_0;
    float _S295 = _S293.primal_0 * _S293.primal_0;
    float _S296 = (_S287.differential_0 * _S293.primal_0 - _S287.primal_0 * _S293.differential_0) / _S295;
    float _S297 = _S293.differential_0 * _S293.primal_0;
    float _S298 = _S297 + _S297;
    float3  _S299 = make_float3 (2.0f) * dpscales_0->primal_0;
    float3  _S300 = dpscales_0->differential_0 * make_float3 (2.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S301;
    (&_S301)->primal_0 = _S299;
    (&_S301)->differential_0 = _S300;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S302 = s_fwd_s_primal_ctx_exp_1(&_S301);
    float x_5 = _S302.primal_0.x;
    float s_diff_x_0 = _S302.differential_0.x;
    float y_3 = _S302.primal_0.y;
    float s_diff_y_0 = _S302.differential_0.y;
    float z_2 = _S302.primal_0.z;
    float s_diff_z_0 = _S302.differential_0.z;
    float s_2 = x_5 + y_3 + z_2;
    float s_diff_s_0 = s_diff_x_0 + s_diff_y_0 + s_diff_z_0;
    DiffPair_float_0 _S303;
    (&_S303)->primal_0 = x_5;
    (&_S303)->differential_0 = s_diff_x_0;
    DiffPair_float_0 _S304;
    (&_S304)->primal_0 = y_3;
    (&_S304)->differential_0 = s_diff_y_0;
    DiffPair_float_0 _S305 = _d_max_1(&_S303, &_S304);
    DiffPair_float_0 _S306;
    (&_S306)->primal_0 = _S305.primal_0;
    (&_S306)->differential_0 = _S305.differential_0;
    DiffPair_float_0 _S307;
    (&_S307)->primal_0 = z_2;
    (&_S307)->differential_0 = s_diff_z_0;
    DiffPair_float_0 _S308 = _d_max_1(&_S306, &_S307);
    float s1_2 = _S308.primal_0 / s_2;
    float _S309 = s_2 * s_2;
    float s_diff_s1_0 = (_S308.differential_0 * s_2 - _S308.primal_0 * s_diff_s_0) / _S309;
    float _S310 = s_diff_s_0 * s_2;
    float _S311 = _S310 + _S310;
    DiffPair_float_0 _S312;
    (&_S312)->primal_0 = x_5;
    (&_S312)->differential_0 = s_diff_x_0;
    DiffPair_float_0 _S313;
    (&_S313)->primal_0 = y_3;
    (&_S313)->differential_0 = s_diff_y_0;
    DiffPair_float_0 _S314 = _d_min_1(&_S312, &_S313);
    DiffPair_float_0 _S315;
    (&_S315)->primal_0 = _S314.primal_0;
    (&_S315)->differential_0 = _S314.differential_0;
    DiffPair_float_0 _S316;
    (&_S316)->primal_0 = z_2;
    (&_S316)->differential_0 = s_diff_z_0;
    DiffPair_float_0 _S317 = _d_min_1(&_S315, &_S316);
    float s3_2 = _S317.primal_0 / s_2;
    float s_diff_s3_0 = (_S317.differential_0 * s_2 - _S317.primal_0 * s_diff_s_0) / _S309;
    float s2_2 = 1.0f - s1_2 - s3_2;
    float s_diff_s2_0 = 0.0f - s_diff_s1_0 - s_diff_s3_0;
    float _S318 = - s1_2;
    float _S319 = - s_diff_s1_0;
    DiffPair_float_0 _S320;
    (&_S320)->primal_0 = s1_2;
    (&_S320)->differential_0 = s_diff_s1_0;
    DiffPair_float_0 _S321 = s_fwd_s_primal_ctx_log_0(&_S320);
    float _S322 = _S318 * _S321.primal_0;
    float _S323 = _S319 * _S321.primal_0 + _S321.differential_0 * _S318;
    DiffPair_float_0 _S324;
    (&_S324)->primal_0 = s2_2;
    (&_S324)->differential_0 = s_diff_s2_0;
    DiffPair_float_0 _S325 = s_fwd_s_primal_ctx_log_0(&_S324);
    float _S326 = _S322 - s2_2 * _S325.primal_0;
    float _S327 = _S323 - (s_diff_s2_0 * _S325.primal_0 + _S325.differential_0 * s2_2);
    DiffPair_float_0 _S328;
    (&_S328)->primal_0 = s3_2;
    (&_S328)->differential_0 = s_diff_s3_0;
    DiffPair_float_0 _S329 = s_fwd_s_primal_ctx_log_0(&_S328);
    float _S330 = _S326 - s3_2 * _S329.primal_0;
    float _S331 = _S327 - (s_diff_s3_0 * _S329.primal_0 + _S329.differential_0 * s3_2);
    DiffPair_float_0 _S332;
    (&_S332)->primal_0 = _S330;
    (&_S332)->differential_0 = _S331;
    DiffPair_float_0 _S333 = s_fwd_s_primal_ctx_exp_0(&_S332);
    float _S334 = _S333.primal_0 - 0.99998998641967773f;
    DiffPair_float_0 _S335;
    (&_S335)->primal_0 = _S334;
    (&_S335)->differential_0 = _S333.differential_0;
    DiffPair_float_0 _S336 = s_fwd_s_primal_ctx_log_0(&_S335);
    float _S337 = erank_reg_weight_s3_2 * (*v_loss_1)[int(3)];
    float _S338 = erank_reg_weight_2 * (*v_loss_1)[int(3)];
    DiffPair_float_0 _S339 = { - _S336.primal_0, 0.0f };
    DiffPair_float_0 _S340 = { - _S336.differential_0, 0.0f };
    DiffPair_float_0 _S341 = { 0.0f, 0.0f };
    DiffPair_0 _S342;
    (&_S342)->primal_0 = _S339;
    (&_S342)->differential_0 = _S340;
    DiffPair_0 _S343;
    (&_S343)->primal_0 = _S341;
    (&_S343)->differential_0 = _S341;
    DiffPair_float_0 _S344;
    (&_S344)->primal_0 = _S338;
    (&_S344)->differential_0 = 0.0f;
    s_fwd_d_max_0(&_S342, &_S343, &_S344);
    float _S345 = - _S342.primal_0.differential_0;
    float _S346 = - _S342.differential_0.differential_0;
    DiffPair_float_0 _S347 = { _S334, 0.0f };
    DiffPair_float_0 _S348 = { _S333.differential_0, 0.0f };
    DiffPair_0 _S349;
    (&_S349)->primal_0 = _S347;
    (&_S349)->differential_0 = _S348;
    DiffPair_float_0 _S350;
    (&_S350)->primal_0 = _S345;
    (&_S350)->differential_0 = _S346;
    s_fwd_s_bwd_prop_log_0(&_S349, &_S350);
    DiffPair_float_0 _S351 = { _S330, 0.0f };
    DiffPair_float_0 _S352 = { _S331, 0.0f };
    DiffPair_0 _S353;
    (&_S353)->primal_0 = _S351;
    (&_S353)->differential_0 = _S352;
    DiffPair_float_0 _S354;
    (&_S354)->primal_0 = _S349.primal_0.differential_0;
    (&_S354)->differential_0 = _S349.differential_0.differential_0;
    s_fwd_s_bwd_prop_exp_0(&_S353, &_S354);
    float _S355 = - _S353.primal_0.differential_0;
    float _S356 = - _S353.differential_0.differential_0;
    float _S357 = s3_2 * _S355;
    float _S358 = s_diff_s3_0 * _S355 + _S356 * s3_2;
    float _S359 = _S329.primal_0 * _S355;
    float _S360 = _S329.differential_0 * _S355 + _S356 * _S329.primal_0;
    DiffPair_float_0 _S361 = { s3_2, 0.0f };
    DiffPair_float_0 _S362 = { s_diff_s3_0, 0.0f };
    DiffPair_0 _S363;
    (&_S363)->primal_0 = _S361;
    (&_S363)->differential_0 = _S362;
    DiffPair_float_0 _S364;
    (&_S364)->primal_0 = _S357;
    (&_S364)->differential_0 = _S358;
    s_fwd_s_bwd_prop_log_0(&_S363, &_S364);
    float _S365 = s2_2 * _S355;
    float _S366 = s_diff_s2_0 * _S355 + _S356 * s2_2;
    float _S367 = _S325.primal_0 * _S355;
    float _S368 = _S325.differential_0 * _S355 + _S356 * _S325.primal_0;
    DiffPair_float_0 _S369 = { s2_2, 0.0f };
    DiffPair_float_0 _S370 = { s_diff_s2_0, 0.0f };
    DiffPair_0 _S371;
    (&_S371)->primal_0 = _S369;
    (&_S371)->differential_0 = _S370;
    DiffPair_float_0 _S372;
    (&_S372)->primal_0 = _S365;
    (&_S372)->differential_0 = _S366;
    s_fwd_s_bwd_prop_log_0(&_S371, &_S372);
    float _S373 = _S318 * _S353.primal_0.differential_0;
    float _S374 = _S319 * _S353.primal_0.differential_0 + _S353.differential_0.differential_0 * _S318;
    float _S375 = _S321.primal_0 * _S353.primal_0.differential_0;
    float _S376 = _S321.differential_0 * _S353.primal_0.differential_0 + _S353.differential_0.differential_0 * _S321.primal_0;
    DiffPair_float_0 _S377 = { s1_2, 0.0f };
    DiffPair_float_0 _S378 = { s_diff_s1_0, 0.0f };
    DiffPair_0 _S379;
    (&_S379)->primal_0 = _S377;
    (&_S379)->differential_0 = _S378;
    DiffPair_float_0 _S380;
    (&_S380)->primal_0 = _S373;
    (&_S380)->differential_0 = _S374;
    s_fwd_s_bwd_prop_log_0(&_S379, &_S380);
    float _S381 = - _S375;
    float _S382 = - _S376;
    float _S383 = - (_S367 + _S371.primal_0.differential_0);
    float _S384 = - (_S368 + _S371.differential_0.differential_0);
    float _S385 = _S337 + _S359 + _S363.primal_0.differential_0 + _S383;
    float _S386 = _S385 / _S309;
    float _S387 = _S309 * _S309;
    float _S388 = ((_S360 + _S363.differential_0.differential_0 + _S384) * _S309 - _S385 * _S311) / _S387;
    float _S389 = - _S386;
    float _S390 = _S317.primal_0 * _S389;
    float _S391 = _S317.differential_0 * _S389 + - _S388 * _S317.primal_0;
    float _S392 = s_2 * _S386;
    float _S393 = s_diff_s_0 * _S386 + _S388 * s_2;
    DiffPair_float_0 _S394 = { _S314.primal_0, 0.0f };
    DiffPair_float_0 _S395 = { _S314.differential_0, 0.0f };
    DiffPair_float_0 _S396 = { z_2, 0.0f };
    DiffPair_float_0 _S397 = { s_diff_z_0, 0.0f };
    DiffPair_0 _S398;
    (&_S398)->primal_0 = _S394;
    (&_S398)->differential_0 = _S395;
    DiffPair_0 _S399;
    (&_S399)->primal_0 = _S396;
    (&_S399)->differential_0 = _S397;
    DiffPair_float_0 _S400;
    (&_S400)->primal_0 = _S392;
    (&_S400)->differential_0 = _S393;
    s_fwd_d_min_0(&_S398, &_S399, &_S400);
    DiffPair_float_0 _S401 = { x_5, 0.0f };
    DiffPair_float_0 _S402 = { s_diff_x_0, 0.0f };
    DiffPair_float_0 _S403 = { y_3, 0.0f };
    DiffPair_float_0 _S404 = { s_diff_y_0, 0.0f };
    DiffPair_0 _S405;
    (&_S405)->primal_0 = _S401;
    (&_S405)->differential_0 = _S402;
    DiffPair_0 _S406;
    (&_S406)->primal_0 = _S403;
    (&_S406)->differential_0 = _S404;
    DiffPair_float_0 _S407;
    (&_S407)->primal_0 = _S398.primal_0.differential_0;
    (&_S407)->differential_0 = _S398.differential_0.differential_0;
    s_fwd_d_min_0(&_S405, &_S406, &_S407);
    float _S408 = _S379.primal_0.differential_0 + _S381 + _S383;
    float _S409 = _S408 / _S309;
    float _S410 = ((_S379.differential_0.differential_0 + _S382 + _S384) * _S309 - _S408 * _S311) / _S387;
    float _S411 = - _S409;
    float _S412 = _S308.primal_0 * _S411;
    float _S413 = _S308.differential_0 * _S411 + - _S410 * _S308.primal_0;
    float _S414 = s_2 * _S409;
    float _S415 = s_diff_s_0 * _S409 + _S410 * s_2;
    DiffPair_float_0 _S416 = { _S305.primal_0, 0.0f };
    DiffPair_float_0 _S417 = { _S305.differential_0, 0.0f };
    DiffPair_0 _S418;
    (&_S418)->primal_0 = _S416;
    (&_S418)->differential_0 = _S417;
    DiffPair_0 _S419;
    (&_S419)->primal_0 = _S396;
    (&_S419)->differential_0 = _S397;
    DiffPair_float_0 _S420;
    (&_S420)->primal_0 = _S414;
    (&_S420)->differential_0 = _S415;
    s_fwd_d_max_0(&_S418, &_S419, &_S420);
    DiffPair_0 _S421;
    (&_S421)->primal_0 = _S401;
    (&_S421)->differential_0 = _S402;
    DiffPair_0 _S422;
    (&_S422)->primal_0 = _S403;
    (&_S422)->differential_0 = _S404;
    DiffPair_float_0 _S423;
    (&_S423)->primal_0 = _S418.primal_0.differential_0;
    (&_S423)->differential_0 = _S418.differential_0.differential_0;
    s_fwd_d_max_0(&_S421, &_S422, &_S423);
    float _S424 = _S390 + _S412;
    float _S425 = _S391 + _S413;
    float3  _S426 = make_float3 (_S405.primal_0.differential_0 + _S421.primal_0.differential_0 + _S424, _S406.primal_0.differential_0 + _S422.primal_0.differential_0 + _S424, _S399.primal_0.differential_0 + _S419.primal_0.differential_0 + _S424);
    float3  _S427 = make_float3 (_S405.differential_0.differential_0 + _S421.differential_0.differential_0 + _S425, _S406.differential_0.differential_0 + _S422.differential_0.differential_0 + _S425, _S399.differential_0.differential_0 + _S419.differential_0.differential_0 + _S425);
    float3  _S428 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S429 = { _S299, _S428 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S430 = { _S300, _S428 };
    DiffPair_1 _S431;
    (&_S431)->primal_0 = _S429;
    (&_S431)->differential_0 = _S430;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S432;
    (&_S432)->primal_0 = _S426;
    (&_S432)->differential_0 = _S427;
    s_fwd_s_bwd_prop_exp_1(&_S431, &_S432);
    float3  _S433 = make_float3 (2.0f) * _S431.primal_0.differential_0;
    float3  _S434 = _S431.differential_0.differential_0 * make_float3 (2.0f);
    float s_diff_scale_reg_T_1 = scale_regularization_weight_2 * (*v_loss_1)[int(2)];
    DiffPair_float_0 _S435 = { _S294, 0.0f };
    DiffPair_float_0 _S436 = { _S296, 0.0f };
    DiffPair_float_0 _S437 = { max_gauss_ratio_2, 0.0f };
    DiffPair_0 _S438;
    (&_S438)->primal_0 = _S435;
    (&_S438)->differential_0 = _S436;
    DiffPair_0 _S439;
    (&_S439)->primal_0 = _S437;
    (&_S439)->differential_0 = _S341;
    DiffPair_float_0 _S440;
    (&_S440)->primal_0 = s_diff_scale_reg_T_1;
    (&_S440)->differential_0 = 0.0f;
    s_fwd_d_max_0(&_S438, &_S439, &_S440);
    float _S441 = _S438.primal_0.differential_0 / _S295;
    float _S442 = (_S438.differential_0.differential_0 * _S295 - _S438.primal_0.differential_0 * _S298) / (_S295 * _S295);
    float _S443 = - _S441;
    float _S444 = _S287.primal_0 * _S443;
    float _S445 = _S287.differential_0 * _S443 + - _S442 * _S287.primal_0;
    float _S446 = _S293.primal_0 * _S441;
    float _S447 = _S293.differential_0 * _S441 + _S442 * _S293.primal_0;
    DiffPair_float_0 _S448 = { _S290.primal_0, 0.0f };
    DiffPair_float_0 _S449 = { _S290.differential_0, 0.0f };
    DiffPair_float_0 _S450 = { _S280, 0.0f };
    DiffPair_float_0 _S451 = { _S281, 0.0f };
    DiffPair_0 _S452;
    (&_S452)->primal_0 = _S448;
    (&_S452)->differential_0 = _S449;
    DiffPair_0 _S453;
    (&_S453)->primal_0 = _S450;
    (&_S453)->differential_0 = _S451;
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S444;
    (&_S454)->differential_0 = _S445;
    s_fwd_d_min_0(&_S452, &_S453, &_S454);
    DiffPair_float_0 _S455 = { _S276, 0.0f };
    DiffPair_float_0 _S456 = { _S277, 0.0f };
    DiffPair_float_0 _S457 = { _S278, 0.0f };
    DiffPair_float_0 _S458 = { _S279, 0.0f };
    DiffPair_0 _S459;
    (&_S459)->primal_0 = _S455;
    (&_S459)->differential_0 = _S456;
    DiffPair_0 _S460;
    (&_S460)->primal_0 = _S457;
    (&_S460)->differential_0 = _S458;
    DiffPair_float_0 _S461;
    (&_S461)->primal_0 = _S452.primal_0.differential_0;
    (&_S461)->differential_0 = _S452.differential_0.differential_0;
    s_fwd_d_min_0(&_S459, &_S460, &_S461);
    DiffPair_float_0 _S462 = { _S284.primal_0, 0.0f };
    DiffPair_float_0 _S463 = { _S284.differential_0, 0.0f };
    DiffPair_0 _S464;
    (&_S464)->primal_0 = _S462;
    (&_S464)->differential_0 = _S463;
    DiffPair_0 _S465;
    (&_S465)->primal_0 = _S450;
    (&_S465)->differential_0 = _S451;
    DiffPair_float_0 _S466;
    (&_S466)->primal_0 = _S446;
    (&_S466)->differential_0 = _S447;
    s_fwd_d_max_0(&_S464, &_S465, &_S466);
    DiffPair_0 _S467;
    (&_S467)->primal_0 = _S455;
    (&_S467)->differential_0 = _S456;
    DiffPair_0 _S468;
    (&_S468)->primal_0 = _S457;
    (&_S468)->differential_0 = _S458;
    DiffPair_float_0 _S469;
    (&_S469)->primal_0 = _S464.primal_0.differential_0;
    (&_S469)->differential_0 = _S464.differential_0.differential_0;
    s_fwd_d_max_0(&_S467, &_S468, &_S469);
    float _S470 = mcmc_scale_reg_weight_2 * (0.3333333432674408f * (*v_loss_1)[int(1)]);
    float3  _S471 = make_float3 (_S459.primal_0.differential_0 + _S467.primal_0.differential_0 + _S470, _S460.primal_0.differential_0 + _S468.primal_0.differential_0 + _S470, _S453.primal_0.differential_0 + _S465.primal_0.differential_0 + _S470);
    float3  _S472 = make_float3 (_S459.differential_0.differential_0 + _S467.differential_0.differential_0, _S460.differential_0.differential_0 + _S468.differential_0.differential_0, _S453.differential_0.differential_0 + _S465.differential_0.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S473 = { dpscales_0->primal_0, _S428 };
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S474 = { dpscales_0->differential_0, _S428 };
    DiffPair_1 _S475;
    (&_S475)->primal_0 = _S473;
    (&_S475)->differential_0 = _S474;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S476;
    (&_S476)->primal_0 = _S471;
    (&_S476)->differential_0 = _S472;
    s_fwd_s_bwd_prop_exp_1(&_S475, &_S476);
    float s_diff_quat_norm_reg_T_1 = quat_norm_reg_weight_2 * (*v_loss_1)[int(4)];
    float _S477 = - s_diff_quat_norm_reg_T_1;
    DiffPair_float_0 _S478 = { _S273.primal_0, 0.0f };
    DiffPair_float_0 _S479 = { _S273.differential_0, 0.0f };
    DiffPair_0 _S480;
    (&_S480)->primal_0 = _S478;
    (&_S480)->differential_0 = _S479;
    DiffPair_float_0 _S481;
    (&_S481)->primal_0 = _S477;
    (&_S481)->differential_0 = -0.0f;
    s_fwd_s_bwd_prop_log_0(&_S480, &_S481);
    float _S482 = _S480.primal_0.differential_0 + s_diff_quat_norm_reg_T_1;
    float4  _S483 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S484 = { dpquat_0->primal_0, _S483 };
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S485 = { dpquat_0->differential_0, _S483 };
    DiffPair_2 _S486;
    (&_S486)->primal_0 = _S484;
    (&_S486)->differential_0 = _S485;
    DiffPair_float_0 _S487;
    (&_S487)->primal_0 = _S482;
    (&_S487)->differential_0 = _S480.differential_0.differential_0;
    s_fwd_s_bwd_length_impl_0(&_S486, &_S487);
    float s_diff_reg_T_0 = mcmc_opacity_reg_weight_2 * (*v_loss_1)[int(0)];
    float _S488 = - (s_diff_reg_T_0 / _S269);
    float _S489 = - ((0.0f - s_diff_reg_T_0 * _S271) / (_S269 * _S269));
    DiffPair_float_0 _S490 = { _S264, 0.0f };
    DiffPair_float_0 _S491 = { _S265, 0.0f };
    DiffPair_0 _S492;
    (&_S492)->primal_0 = _S490;
    (&_S492)->differential_0 = _S491;
    DiffPair_float_0 _S493;
    (&_S493)->primal_0 = _S488;
    (&_S493)->differential_0 = _S489;
    s_fwd_s_bwd_prop_exp_0(&_S492, &_S493);
    float _S494 = - _S492.primal_0.differential_0;
    float _S495 = - _S492.differential_0.differential_0;
    float3  _S496 = _S434 + _S475.differential_0.differential_0;
    dpv_scales_0->primal_0 = _S433 + _S475.primal_0.differential_0;
    dpv_scales_0->differential_0 = _S496;
    dpv_opacity_0->primal_0 = _S494;
    dpv_opacity_0->differential_0 = _S495;
    dpv_quat_0->primal_0 = _S486.primal_0.differential_0;
    dpv_quat_0->differential_0 = _S486.differential_0.differential_0;
    return;
}

inline __device__ void per_splat_losses_bwd(float3  scales_2, float opacity_2, float4  quat_2, FixedArray<float, 5>  v_loss_2, float3  * v_scales_1, float * v_opacity_1, float4  * v_quat_1, float3  * vr_scales_0, float * vr_opacity_0, float4  * vr_quat_0, float3  * h_scales_0, float * h_opacity_0, float4  * h_quat_0, float mcmc_opacity_reg_weight_3, float mcmc_scale_reg_weight_3, float max_gauss_ratio_3, float scale_regularization_weight_3, float erank_reg_weight_3, float erank_reg_weight_s3_3, float quat_norm_reg_weight_3)
{
    float _S497 = - opacity_2;
    float _S498 = 1.0f + s_primal_ctx_exp_0(_S497);
    float _S499 = _S498 * _S498;
    float _S500 = length_0(quat_2);
    float3  _S501 = s_primal_ctx_exp_1(scales_2);
    float _S502 = _S501.x;
    float _S503 = _S501.y;
    float _S504 = _S501.z;
    float _S505 = (F32_max((_S502), (_S503)));
    float _S506 = (F32_max((_S505), (_S504)));
    float _S507 = (F32_min((_S502), (_S503)));
    float _S508 = (F32_min((_S507), (_S504)));
    float _S509 = _S506 / _S508;
    float _S510 = _S508 * _S508;
    float3  _S511 = make_float3 (2.0f) * scales_2;
    float3  _S512 = s_primal_ctx_exp_1(_S511);
    float x_6 = _S512.x;
    float y_4 = _S512.y;
    float z_3 = _S512.z;
    float s_3 = x_6 + y_4 + z_3;
    float _S513 = (F32_max((x_6), (y_4)));
    float _S514 = (F32_max((_S513), (z_3)));
    float s1_3 = _S514 / s_3;
    float _S515 = s_3 * s_3;
    float _S516 = (F32_min((x_6), (y_4)));
    float _S517 = (F32_min((_S516), (z_3)));
    float s3_3 = _S517 / s_3;
    float s2_3 = 1.0f - s1_3 - s3_3;
    float _S518 = - s1_3;
    float _S519 = s_primal_ctx_log_0(s1_3);
    float _S520 = s_primal_ctx_log_0(s2_3);
    float _S521 = s_primal_ctx_log_0(s3_3);
    float _S522 = _S518 * _S519 - s2_3 * _S520 - s3_3 * _S521;
    float _S523 = s_primal_ctx_exp_0(_S522) - 0.99998998641967773f;
    float _S524 = erank_reg_weight_s3_3 * v_loss_2[int(3)];
    float _S525 = erank_reg_weight_3 * v_loss_2[int(3)];
    DiffPair_float_0 _S526;
    (&_S526)->primal_0 = - s_primal_ctx_log_0(_S523);
    (&_S526)->differential_0 = 0.0f;
    DiffPair_float_0 _S527;
    (&_S527)->primal_0 = 0.0f;
    (&_S527)->differential_0 = 0.0f;
    _d_max_0(&_S526, &_S527, _S525);
    float _S528 = - _S526.differential_0;
    DiffPair_float_0 _S529;
    (&_S529)->primal_0 = _S523;
    (&_S529)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S529, _S528);
    DiffPair_float_0 _S530;
    (&_S530)->primal_0 = _S522;
    (&_S530)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S530, _S529.differential_0);
    float _S531 = - _S530.differential_0;
    float _S532 = s3_3 * _S531;
    float _S533 = _S521 * _S531;
    DiffPair_float_0 _S534;
    (&_S534)->primal_0 = s3_3;
    (&_S534)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S534, _S532);
    float _S535 = s2_3 * _S531;
    float _S536 = _S520 * _S531;
    DiffPair_float_0 _S537;
    (&_S537)->primal_0 = s2_3;
    (&_S537)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S537, _S535);
    float _S538 = _S518 * _S530.differential_0;
    float _S539 = _S519 * _S530.differential_0;
    DiffPair_float_0 _S540;
    (&_S540)->primal_0 = s1_3;
    (&_S540)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S540, _S538);
    float _S541 = - _S539;
    float _S542 = - (_S536 + _S537.differential_0);
    float _S543 = (_S524 + _S533 + _S534.differential_0 + _S542) / _S515;
    float _S544 = _S517 * - _S543;
    float _S545 = s_3 * _S543;
    DiffPair_float_0 _S546;
    (&_S546)->primal_0 = _S516;
    (&_S546)->differential_0 = 0.0f;
    DiffPair_float_0 _S547;
    (&_S547)->primal_0 = z_3;
    (&_S547)->differential_0 = 0.0f;
    _d_min_0(&_S546, &_S547, _S545);
    DiffPair_float_0 _S548;
    (&_S548)->primal_0 = x_6;
    (&_S548)->differential_0 = 0.0f;
    DiffPair_float_0 _S549;
    (&_S549)->primal_0 = y_4;
    (&_S549)->differential_0 = 0.0f;
    _d_min_0(&_S548, &_S549, _S546.differential_0);
    float _S550 = (_S540.differential_0 + _S541 + _S542) / _S515;
    float _S551 = _S514 * - _S550;
    float _S552 = s_3 * _S550;
    DiffPair_float_0 _S553;
    (&_S553)->primal_0 = _S513;
    (&_S553)->differential_0 = 0.0f;
    DiffPair_float_0 _S554;
    (&_S554)->primal_0 = z_3;
    (&_S554)->differential_0 = 0.0f;
    _d_max_0(&_S553, &_S554, _S552);
    DiffPair_float_0 _S555;
    (&_S555)->primal_0 = x_6;
    (&_S555)->differential_0 = 0.0f;
    DiffPair_float_0 _S556;
    (&_S556)->primal_0 = y_4;
    (&_S556)->differential_0 = 0.0f;
    _d_max_0(&_S555, &_S556, _S553.differential_0);
    float _S557 = _S544 + _S551;
    float3  _S558 = make_float3 (_S548.differential_0 + _S555.differential_0 + _S557, _S549.differential_0 + _S556.differential_0 + _S557, _S547.differential_0 + _S554.differential_0 + _S557);
    float3  _S559 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S560;
    (&_S560)->primal_0 = _S511;
    (&_S560)->differential_0 = _S559;
    s_bwd_prop_exp_1(&_S560, _S558);
    float3  _S561 = make_float3 (2.0f) * _S560.differential_0;
    float s_diff_scale_reg_T_2 = scale_regularization_weight_3 * v_loss_2[int(2)];
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S509;
    (&_S562)->differential_0 = 0.0f;
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = max_gauss_ratio_3;
    (&_S563)->differential_0 = 0.0f;
    _d_max_0(&_S562, &_S563, s_diff_scale_reg_T_2);
    float _S564 = _S562.differential_0 / _S510;
    float _S565 = _S506 * - _S564;
    float _S566 = _S508 * _S564;
    DiffPair_float_0 _S567;
    (&_S567)->primal_0 = _S507;
    (&_S567)->differential_0 = 0.0f;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S504;
    (&_S568)->differential_0 = 0.0f;
    _d_min_0(&_S567, &_S568, _S565);
    DiffPair_float_0 _S569;
    (&_S569)->primal_0 = _S502;
    (&_S569)->differential_0 = 0.0f;
    DiffPair_float_0 _S570;
    (&_S570)->primal_0 = _S503;
    (&_S570)->differential_0 = 0.0f;
    _d_min_0(&_S569, &_S570, _S567.differential_0);
    DiffPair_float_0 _S571;
    (&_S571)->primal_0 = _S505;
    (&_S571)->differential_0 = 0.0f;
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S504;
    (&_S572)->differential_0 = 0.0f;
    _d_max_0(&_S571, &_S572, _S566);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S502;
    (&_S573)->differential_0 = 0.0f;
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = _S503;
    (&_S574)->differential_0 = 0.0f;
    _d_max_0(&_S573, &_S574, _S571.differential_0);
    float _S575 = mcmc_scale_reg_weight_3 * (0.3333333432674408f * v_loss_2[int(1)]);
    float3  _S576 = make_float3 (_S569.differential_0 + _S573.differential_0 + _S575, _S570.differential_0 + _S574.differential_0 + _S575, _S568.differential_0 + _S572.differential_0 + _S575);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S577;
    (&_S577)->primal_0 = scales_2;
    (&_S577)->differential_0 = _S559;
    s_bwd_prop_exp_1(&_S577, _S576);
    float s_diff_quat_norm_reg_T_2 = quat_norm_reg_weight_3 * v_loss_2[int(4)];
    float _S578 = - s_diff_quat_norm_reg_T_2;
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S500;
    (&_S579)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S579, _S578);
    float _S580 = _S579.differential_0 + s_diff_quat_norm_reg_T_2;
    float4  _S581 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S582;
    (&_S582)->primal_0 = quat_2;
    (&_S582)->differential_0 = _S581;
    s_bwd_length_impl_0(&_S582, _S580);
    float _S583 = - (mcmc_opacity_reg_weight_3 * v_loss_2[int(0)] / _S499);
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S497;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S584, _S583);
    float _S585 = - _S584.differential_0;
    *v_scales_1 = _S561 + _S577.differential_0;
    *v_opacity_1 = _S585;
    *v_quat_1 = _S582.differential_0;
    FixedArray<float, 5>  losses_1;
    losses_1[int(0)] = mcmc_opacity_reg_weight_3 * (1.0f / (1.0f + (F32_exp((_S497)))));
    losses_1[int(4)] = quat_norm_reg_weight_3 * (_S500 - 1.0f - (F32_log((_S500))));
    float3  _S586 = exp_0(scales_2);
    float _S587 = _S586.x;
    float _S588 = _S586.y;
    float _S589 = _S586.z;
    losses_1[int(1)] = mcmc_scale_reg_weight_3 * (_S587 + _S588 + _S589) / 3.0f;
    losses_1[int(2)] = scale_regularization_weight_3 * ((F32_max(((F32_max(((F32_max((_S587), (_S588)))), (_S589))) / (F32_min(((F32_min((_S587), (_S588)))), (_S589)))), (max_gauss_ratio_3))) - max_gauss_ratio_3);
    float3  _S590 = exp_0(_S511);
    float x_7 = _S590.x;
    float y_5 = _S590.y;
    float z_4 = _S590.z;
    float s_4 = x_7 + y_5 + z_4;
    float s1_4 = (F32_max(((F32_max((x_7), (y_5)))), (z_4))) / s_4;
    float s3_4 = (F32_min(((F32_min((x_7), (y_5)))), (z_4))) / s_4;
    float s2_4 = 1.0f - s1_4 - s3_4;
    losses_1[int(3)] = erank_reg_weight_3 * (F32_max((- (F32_log(((F32_exp((- s1_4 * (F32_log((s1_4))) - s2_4 * (F32_log((s2_4))) - s3_4 * (F32_log((s3_4)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_3 * s3_4;
    FixedArray<float, 5>  residuals_0;
    float _S591 = (F32_sqrt(((F32_max((losses_1[int(0)]), (0.0f))))));
    residuals_0[int(0)] = _S591;
    float _S592 = (F32_sqrt(((F32_max((losses_1[int(1)]), (0.0f))))));
    residuals_0[int(1)] = _S592;
    float _S593 = (F32_sqrt(((F32_max((losses_1[int(2)]), (0.0f))))));
    residuals_0[int(2)] = _S593;
    float _S594 = (F32_sqrt(((F32_max((losses_1[int(3)]), (0.0f))))));
    residuals_0[int(3)] = _S594;
    float _S595 = (F32_sqrt(((F32_max((losses_1[int(4)]), (0.0f))))));
    residuals_0[int(4)] = _S595;
    float _S596 = _S506 / _S508;
    float s1_5 = _S514 / s_3;
    float s3_5 = _S517 / s_3;
    float s2_5 = 1.0f - s1_5 - s3_5;
    float _S597 = - s1_5;
    float _S598 = s_primal_ctx_log_0(s1_5);
    float _S599 = s_primal_ctx_log_0(s2_5);
    float _S600 = s_primal_ctx_log_0(s3_5);
    float _S601 = _S597 * _S598 - s2_5 * _S599 - s3_5 * _S600;
    float _S602 = s_primal_ctx_exp_0(_S601) - 0.99998998641967773f;
    float _S603 = erank_reg_weight_s3_3 * _S594;
    float _S604 = erank_reg_weight_3 * _S594;
    DiffPair_float_0 _S605;
    (&_S605)->primal_0 = - s_primal_ctx_log_0(_S602);
    (&_S605)->differential_0 = 0.0f;
    DiffPair_float_0 _S606;
    (&_S606)->primal_0 = 0.0f;
    (&_S606)->differential_0 = 0.0f;
    _d_max_0(&_S605, &_S606, _S604);
    float _S607 = - _S605.differential_0;
    DiffPair_float_0 _S608;
    (&_S608)->primal_0 = _S602;
    (&_S608)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S608, _S607);
    DiffPair_float_0 _S609;
    (&_S609)->primal_0 = _S601;
    (&_S609)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S609, _S608.differential_0);
    float _S610 = - _S609.differential_0;
    float _S611 = s3_5 * _S610;
    float _S612 = _S600 * _S610;
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = s3_5;
    (&_S613)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S613, _S611);
    float _S614 = s2_5 * _S610;
    float _S615 = _S599 * _S610;
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = s2_5;
    (&_S616)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S616, _S614);
    float _S617 = _S597 * _S609.differential_0;
    float _S618 = _S598 * _S609.differential_0;
    DiffPair_float_0 _S619;
    (&_S619)->primal_0 = s1_5;
    (&_S619)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S619, _S617);
    float _S620 = - _S618;
    float _S621 = - (_S615 + _S616.differential_0);
    float _S622 = (_S603 + _S612 + _S613.differential_0 + _S621) / _S515;
    float _S623 = _S517 * - _S622;
    float _S624 = s_3 * _S622;
    DiffPair_float_0 _S625;
    (&_S625)->primal_0 = _S516;
    (&_S625)->differential_0 = 0.0f;
    DiffPair_float_0 _S626;
    (&_S626)->primal_0 = z_3;
    (&_S626)->differential_0 = 0.0f;
    _d_min_0(&_S625, &_S626, _S624);
    DiffPair_float_0 _S627;
    (&_S627)->primal_0 = x_6;
    (&_S627)->differential_0 = 0.0f;
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = y_4;
    (&_S628)->differential_0 = 0.0f;
    _d_min_0(&_S627, &_S628, _S625.differential_0);
    float _S629 = (_S619.differential_0 + _S620 + _S621) / _S515;
    float _S630 = _S514 * - _S629;
    float _S631 = s_3 * _S629;
    DiffPair_float_0 _S632;
    (&_S632)->primal_0 = _S513;
    (&_S632)->differential_0 = 0.0f;
    DiffPair_float_0 _S633;
    (&_S633)->primal_0 = z_3;
    (&_S633)->differential_0 = 0.0f;
    _d_max_0(&_S632, &_S633, _S631);
    DiffPair_float_0 _S634;
    (&_S634)->primal_0 = x_6;
    (&_S634)->differential_0 = 0.0f;
    DiffPair_float_0 _S635;
    (&_S635)->primal_0 = y_4;
    (&_S635)->differential_0 = 0.0f;
    _d_max_0(&_S634, &_S635, _S632.differential_0);
    float _S636 = _S623 + _S630;
    float3  _S637 = make_float3 (_S627.differential_0 + _S634.differential_0 + _S636, _S628.differential_0 + _S635.differential_0 + _S636, _S626.differential_0 + _S633.differential_0 + _S636);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S511;
    (&_S638)->differential_0 = _S559;
    s_bwd_prop_exp_1(&_S638, _S637);
    float3  _S639 = make_float3 (2.0f) * _S638.differential_0;
    float s_diff_scale_reg_T_3 = scale_regularization_weight_3 * _S593;
    DiffPair_float_0 _S640;
    (&_S640)->primal_0 = _S596;
    (&_S640)->differential_0 = 0.0f;
    DiffPair_float_0 _S641;
    (&_S641)->primal_0 = max_gauss_ratio_3;
    (&_S641)->differential_0 = 0.0f;
    _d_max_0(&_S640, &_S641, s_diff_scale_reg_T_3);
    float _S642 = _S640.differential_0 / _S510;
    float _S643 = _S506 * - _S642;
    float _S644 = _S508 * _S642;
    DiffPair_float_0 _S645;
    (&_S645)->primal_0 = _S507;
    (&_S645)->differential_0 = 0.0f;
    DiffPair_float_0 _S646;
    (&_S646)->primal_0 = _S504;
    (&_S646)->differential_0 = 0.0f;
    _d_min_0(&_S645, &_S646, _S643);
    DiffPair_float_0 _S647;
    (&_S647)->primal_0 = _S502;
    (&_S647)->differential_0 = 0.0f;
    DiffPair_float_0 _S648;
    (&_S648)->primal_0 = _S503;
    (&_S648)->differential_0 = 0.0f;
    _d_min_0(&_S647, &_S648, _S645.differential_0);
    DiffPair_float_0 _S649;
    (&_S649)->primal_0 = _S505;
    (&_S649)->differential_0 = 0.0f;
    DiffPair_float_0 _S650;
    (&_S650)->primal_0 = _S504;
    (&_S650)->differential_0 = 0.0f;
    _d_max_0(&_S649, &_S650, _S644);
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = _S502;
    (&_S651)->differential_0 = 0.0f;
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = _S503;
    (&_S652)->differential_0 = 0.0f;
    _d_max_0(&_S651, &_S652, _S649.differential_0);
    float _S653 = mcmc_scale_reg_weight_3 * (0.3333333432674408f * _S592);
    float3  _S654 = make_float3 (_S647.differential_0 + _S651.differential_0 + _S653, _S648.differential_0 + _S652.differential_0 + _S653, _S646.differential_0 + _S650.differential_0 + _S653);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S655;
    (&_S655)->primal_0 = scales_2;
    (&_S655)->differential_0 = _S559;
    s_bwd_prop_exp_1(&_S655, _S654);
    float s_diff_quat_norm_reg_T_3 = quat_norm_reg_weight_3 * _S595;
    float _S656 = - s_diff_quat_norm_reg_T_3;
    DiffPair_float_0 _S657;
    (&_S657)->primal_0 = _S500;
    (&_S657)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S657, _S656);
    float _S658 = _S657.differential_0 + s_diff_quat_norm_reg_T_3;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S659;
    (&_S659)->primal_0 = quat_2;
    (&_S659)->differential_0 = _S581;
    s_bwd_length_impl_0(&_S659, _S658);
    float _S660 = - (mcmc_opacity_reg_weight_3 * _S591 / _S499);
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = _S497;
    (&_S661)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S661, _S660);
    float _S662 = - _S661.differential_0;
    *vr_scales_0 = _S639 + _S655.differential_0;
    *vr_opacity_0 = _S662;
    *vr_quat_0 = _S659.differential_0;
    float3  _S663 = make_float3 (1.0f, 0.0f, 0.0f);
    int4  _S664 = make_int4 (int(0));
    float4  _S665 = make_float4 ((float)_S664.x, (float)_S664.y, (float)_S664.z, (float)_S664.w);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S666;
    (&_S666)->primal_0 = scales_2;
    (&_S666)->differential_0 = _S663;
    DiffPair_float_0 _S667;
    (&_S667)->primal_0 = opacity_2;
    (&_S667)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S668;
    (&_S668)->primal_0 = quat_2;
    (&_S668)->differential_0 = _S665;
    FixedArray<float, 5>  _S669;
    _S669[int(0)] = 1.0f;
    _S669[int(1)] = 1.0f;
    _S669[int(2)] = 1.0f;
    _S669[int(3)] = 1.0f;
    _S669[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_0;
    DiffPair_float_0 dp_v_opacity_0;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_0;
    s_fwd_per_splat_losses_bwd_0(&_S666, &_S667, &_S668, &_S669, &dp_v_scales_0, &dp_v_opacity_0, &dp_v_quat_0, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->x) = dp_v_scales_0.differential_0.x;
    float3  _S670 = make_float3 (0.0f, 1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S671;
    (&_S671)->primal_0 = scales_2;
    (&_S671)->differential_0 = _S670;
    DiffPair_float_0 _S672;
    (&_S672)->primal_0 = opacity_2;
    (&_S672)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S673;
    (&_S673)->primal_0 = quat_2;
    (&_S673)->differential_0 = _S665;
    FixedArray<float, 5>  _S674;
    _S674[int(0)] = 1.0f;
    _S674[int(1)] = 1.0f;
    _S674[int(2)] = 1.0f;
    _S674[int(3)] = 1.0f;
    _S674[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_1;
    DiffPair_float_0 dp_v_opacity_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_1;
    s_fwd_per_splat_losses_bwd_0(&_S671, &_S672, &_S673, &_S674, &dp_v_scales_1, &dp_v_opacity_1, &dp_v_quat_1, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->y) = dp_v_scales_1.differential_0.y;
    float3  _S675 = make_float3 (0.0f, 0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S676;
    (&_S676)->primal_0 = scales_2;
    (&_S676)->differential_0 = _S675;
    DiffPair_float_0 _S677;
    (&_S677)->primal_0 = opacity_2;
    (&_S677)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S678;
    (&_S678)->primal_0 = quat_2;
    (&_S678)->differential_0 = _S665;
    FixedArray<float, 5>  _S679;
    _S679[int(0)] = 1.0f;
    _S679[int(1)] = 1.0f;
    _S679[int(2)] = 1.0f;
    _S679[int(3)] = 1.0f;
    _S679[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_2;
    DiffPair_float_0 dp_v_opacity_2;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_2;
    s_fwd_per_splat_losses_bwd_0(&_S676, &_S677, &_S678, &_S679, &dp_v_scales_2, &dp_v_opacity_2, &dp_v_quat_2, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_scales_0->z) = dp_v_scales_2.differential_0.z;
    int3  _S680 = make_int3 (int(0));
    float3  _S681 = make_float3 ((float)_S680.x, (float)_S680.y, (float)_S680.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S682;
    (&_S682)->primal_0 = scales_2;
    (&_S682)->differential_0 = _S681;
    DiffPair_float_0 _S683;
    (&_S683)->primal_0 = opacity_2;
    (&_S683)->differential_0 = 1.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S684;
    (&_S684)->primal_0 = quat_2;
    (&_S684)->differential_0 = _S665;
    FixedArray<float, 5>  _S685;
    _S685[int(0)] = 1.0f;
    _S685[int(1)] = 1.0f;
    _S685[int(2)] = 1.0f;
    _S685[int(3)] = 1.0f;
    _S685[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_3;
    DiffPair_float_0 dp_v_opacity_3;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_3;
    s_fwd_per_splat_losses_bwd_0(&_S682, &_S683, &_S684, &_S685, &dp_v_scales_3, &dp_v_opacity_3, &dp_v_quat_3, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *h_opacity_0 = dp_v_opacity_3.differential_0;
    float4  _S686 = make_float4 (1.0f, 0.0f, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S687;
    (&_S687)->primal_0 = scales_2;
    (&_S687)->differential_0 = _S681;
    DiffPair_float_0 _S688;
    (&_S688)->primal_0 = opacity_2;
    (&_S688)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S689;
    (&_S689)->primal_0 = quat_2;
    (&_S689)->differential_0 = _S686;
    FixedArray<float, 5>  _S690;
    _S690[int(0)] = 1.0f;
    _S690[int(1)] = 1.0f;
    _S690[int(2)] = 1.0f;
    _S690[int(3)] = 1.0f;
    _S690[int(4)] = 1.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_4;
    DiffPair_float_0 dp_v_opacity_4;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_4;
    s_fwd_per_splat_losses_bwd_0(&_S687, &_S688, &_S689, &_S690, &dp_v_scales_4, &dp_v_opacity_4, &dp_v_quat_4, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->x) = dp_v_quat_4.differential_0.x;
    float4  _S691 = make_float4 (0.0f, 1.0f, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S692;
    (&_S692)->primal_0 = scales_2;
    (&_S692)->differential_0 = _S681;
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
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_5;
    DiffPair_float_0 dp_v_opacity_5;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_5;
    s_fwd_per_splat_losses_bwd_0(&_S692, &_S693, &_S694, &_S695, &dp_v_scales_5, &dp_v_opacity_5, &dp_v_quat_5, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->y) = dp_v_quat_5.differential_0.y;
    float4  _S696 = make_float4 (0.0f, 0.0f, 1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S697;
    (&_S697)->primal_0 = scales_2;
    (&_S697)->differential_0 = _S681;
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
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_6;
    DiffPair_float_0 dp_v_opacity_6;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_6;
    s_fwd_per_splat_losses_bwd_0(&_S697, &_S698, &_S699, &_S700, &dp_v_scales_6, &dp_v_opacity_6, &dp_v_quat_6, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->z) = dp_v_quat_6.differential_0.z;
    float4  _S701 = make_float4 (0.0f, 0.0f, 0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S702;
    (&_S702)->primal_0 = scales_2;
    (&_S702)->differential_0 = _S681;
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
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_v_scales_7;
    DiffPair_float_0 dp_v_opacity_7;
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_v_quat_7;
    s_fwd_per_splat_losses_bwd_0(&_S702, &_S703, &_S704, &_S705, &dp_v_scales_7, &dp_v_opacity_7, &dp_v_quat_7, mcmc_opacity_reg_weight_3, mcmc_scale_reg_weight_3, max_gauss_ratio_3, scale_regularization_weight_3, erank_reg_weight_3, erank_reg_weight_s3_3, quat_norm_reg_weight_3);
    *&(h_quat_0->w) = dp_v_quat_7.differential_0.w;
    return;
}

