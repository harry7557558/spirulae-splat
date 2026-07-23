#pragma once

#include "generated/slang.cuh"

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

inline __device__ float length_0(float4  x_1)
{
    return (F32_sqrt((dot_0(x_1, x_1))));
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

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, float3  dOut_4)
{
    float3  _S7 = exp_0((*dpx_4).primal_0) * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S7;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_5, DiffPair_float_0 * dpy_1, float dOut_5)
{
    DiffPair_float_0 _S8 = *dpx_5;
    float _S9;
    if(((*dpx_5).primal_0) < ((*dpy_1).primal_0))
    {
        _S9 = dOut_5;
    }
    else
    {
        if(((*dpx_5).primal_0) > ((*dpy_1).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_5;
        }
    }
    dpx_5->primal_0 = _S8.primal_0;
    dpx_5->differential_0 = _S9;
    DiffPair_float_0 _S10 = *dpy_1;
    if(((*dpy_1).primal_0) < (_S8.primal_0))
    {
        _S9 = dOut_5;
    }
    else
    {
        if(((*dpy_1).primal_0) > ((*dpx_5).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_5;
        }
    }
    dpy_1->primal_0 = _S10.primal_0;
    dpy_1->differential_0 = _S9;
    return;
}

inline __device__ void per_splat_losses(float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S11)
{
    FixedArray<float, 5>  losses_0;
    losses_0[int(0)] = mcmc_opacity_reg_weight_0 * (1.0f / (1.0f + (F32_exp((- opacity_0)))));
    float quat_norm_0 = length_0(quat_0);
    losses_0[int(4)] = quat_norm_reg_weight_0 * (quat_norm_0 - 1.0f - (F32_log((quat_norm_0))));
    float3  _S12 = exp_0(scales_0);
    float _S13 = _S12.x;
    float _S14 = _S12.y;
    float _S15 = _S12.z;
    losses_0[int(1)] = mcmc_scale_reg_weight_0 * (_S13 + _S14 + _S15) / 3.0f;
    losses_0[int(2)] = scale_regularization_weight_0 * ((F32_max(((F32_max(((F32_max((_S13), (_S14)))), (_S15))) / (F32_min(((F32_min((_S13), (_S14)))), (_S15)))), (max_gauss_ratio_0))) - max_gauss_ratio_0);
    float3  _S16 = exp_0(make_float3 (2.0f) * scales_0);
    float x_3 = _S16.x;
    float y_1 = _S16.y;
    float z_0 = _S16.z;
    float s_0 = x_3 + y_1 + z_0;
    float s1_0 = (F32_max(((F32_max((x_3), (y_1)))), (z_0))) / s_0;
    float s3_0 = (F32_min(((F32_min((x_3), (y_1)))), (z_0))) / s_0;
    float s2_0 = 1.0f - s1_0 - s3_0;
    losses_0[int(3)] = erank_reg_weight_0 * (F32_max((- (F32_log(((F32_exp((- s1_0 * (F32_log((s1_0))) - s2_0 * (F32_log((s2_0))) - s3_0 * (F32_log((s3_0)))))) - 0.99998998641967773f)))), (0.0f))) + erank_reg_weight_s3_0 * s3_0;
    *_S11 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_exp_0(float _S17)
{
    return (F32_exp((_S17)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S18)
{
    return exp_0(_S18);
}

inline __device__ float s_primal_ctx_log_0(float _S19)
{
    return (F32_log((_S19)));
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S20, float _S21)
{
    _d_log_0(_S20, _S21);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S22, float _S23)
{
    _d_exp_0(_S22, _S23);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S24, float3  _S25)
{
    _d_exp_vector_0(_S24, _S25);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S26, float _S27)
{
    _d_sqrt_0(_S26, _S27);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_6, float _s_dOut_0)
{
    float _S28 = (*dpx_6).primal_0.x;
    float _S29 = (*dpx_6).primal_0.y;
    float _S30 = (*dpx_6).primal_0.z;
    float _S31 = (*dpx_6).primal_0.w;
    DiffPair_float_0 _S32;
    (&_S32)->primal_0 = _S28 * _S28 + _S29 * _S29 + _S30 * _S30 + _S31 * _S31;
    (&_S32)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S32, _s_dOut_0);
    float _S33 = (*dpx_6).primal_0.w * _S32.differential_0;
    float _S34 = _S33 + _S33;
    float _S35 = (*dpx_6).primal_0.z * _S32.differential_0;
    float _S36 = _S35 + _S35;
    float _S37 = (*dpx_6).primal_0.y * _S32.differential_0;
    float _S38 = _S37 + _S37;
    float _S39 = (*dpx_6).primal_0.x * _S32.differential_0;
    float _S40 = _S39 + _S39;
    float4  _S41 = make_float4 (0.0f);
    *&((&_S41)->w) = _S34;
    *&((&_S41)->z) = _S36;
    *&((&_S41)->y) = _S38;
    *&((&_S41)->x) = _S40;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S41;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S42, float _S43)
{
    s_bwd_prop_length_impl_0(_S42, _S43);
    return;
}

inline __device__ void per_splat_losses_bwd(float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1)
{
    float _S44 = - opacity_1;
    float _S45 = 1.0f + s_primal_ctx_exp_0(_S44);
    float _S46 = _S45 * _S45;
    float _S47 = length_0(quat_1);
    float3  _S48 = s_primal_ctx_exp_1(scales_1);
    float _S49 = _S48.x;
    float _S50 = _S48.y;
    float _S51 = _S48.z;
    float _S52 = (F32_max((_S49), (_S50)));
    float _S53 = (F32_max((_S52), (_S51)));
    float _S54 = (F32_min((_S49), (_S50)));
    float _S55 = (F32_min((_S54), (_S51)));
    float _S56 = _S53 / _S55;
    float _S57 = _S55 * _S55;
    float3  _S58 = make_float3 (2.0f) * scales_1;
    float3  _S59 = s_primal_ctx_exp_1(_S58);
    float x_4 = _S59.x;
    float y_2 = _S59.y;
    float z_1 = _S59.z;
    float s_1 = x_4 + y_2 + z_1;
    float _S60 = (F32_max((x_4), (y_2)));
    float _S61 = (F32_max((_S60), (z_1)));
    float s1_1 = _S61 / s_1;
    float _S62 = s_1 * s_1;
    float _S63 = (F32_min((x_4), (y_2)));
    float _S64 = (F32_min((_S63), (z_1)));
    float s3_1 = _S64 / s_1;
    float s2_1 = 1.0f - s1_1 - s3_1;
    float _S65 = - s1_1;
    float _S66 = s_primal_ctx_log_0(s1_1);
    float _S67 = s_primal_ctx_log_0(s2_1);
    float _S68 = s_primal_ctx_log_0(s3_1);
    float _S69 = _S65 * _S66 - s2_1 * _S67 - s3_1 * _S68;
    float _S70 = s_primal_ctx_exp_0(_S69) - 0.99998998641967773f;
    float _S71 = erank_reg_weight_s3_1 * v_loss_0[int(3)];
    float _S72 = erank_reg_weight_1 * v_loss_0[int(3)];
    DiffPair_float_0 _S73;
    (&_S73)->primal_0 = - s_primal_ctx_log_0(_S70);
    (&_S73)->differential_0 = 0.0f;
    DiffPair_float_0 _S74;
    (&_S74)->primal_0 = 0.0f;
    (&_S74)->differential_0 = 0.0f;
    _d_max_0(&_S73, &_S74, _S72);
    float _S75 = - _S73.differential_0;
    DiffPair_float_0 _S76;
    (&_S76)->primal_0 = _S70;
    (&_S76)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S76, _S75);
    DiffPair_float_0 _S77;
    (&_S77)->primal_0 = _S69;
    (&_S77)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S77, _S76.differential_0);
    float _S78 = - _S77.differential_0;
    float _S79 = s3_1 * _S78;
    float _S80 = _S68 * _S78;
    DiffPair_float_0 _S81;
    (&_S81)->primal_0 = s3_1;
    (&_S81)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S81, _S79);
    float _S82 = s2_1 * _S78;
    float _S83 = _S67 * _S78;
    DiffPair_float_0 _S84;
    (&_S84)->primal_0 = s2_1;
    (&_S84)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S84, _S82);
    float _S85 = _S65 * _S77.differential_0;
    float _S86 = _S66 * _S77.differential_0;
    DiffPair_float_0 _S87;
    (&_S87)->primal_0 = s1_1;
    (&_S87)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S87, _S85);
    float _S88 = - _S86;
    float _S89 = - (_S83 + _S84.differential_0);
    float _S90 = (_S71 + _S80 + _S81.differential_0 + _S89) / _S62;
    float _S91 = _S64 * - _S90;
    float _S92 = s_1 * _S90;
    DiffPair_float_0 _S93;
    (&_S93)->primal_0 = _S63;
    (&_S93)->differential_0 = 0.0f;
    DiffPair_float_0 _S94;
    (&_S94)->primal_0 = z_1;
    (&_S94)->differential_0 = 0.0f;
    _d_min_0(&_S93, &_S94, _S92);
    DiffPair_float_0 _S95;
    (&_S95)->primal_0 = x_4;
    (&_S95)->differential_0 = 0.0f;
    DiffPair_float_0 _S96;
    (&_S96)->primal_0 = y_2;
    (&_S96)->differential_0 = 0.0f;
    _d_min_0(&_S95, &_S96, _S93.differential_0);
    float _S97 = (_S87.differential_0 + _S88 + _S89) / _S62;
    float _S98 = _S61 * - _S97;
    float _S99 = s_1 * _S97;
    DiffPair_float_0 _S100;
    (&_S100)->primal_0 = _S60;
    (&_S100)->differential_0 = 0.0f;
    DiffPair_float_0 _S101;
    (&_S101)->primal_0 = z_1;
    (&_S101)->differential_0 = 0.0f;
    _d_max_0(&_S100, &_S101, _S99);
    DiffPair_float_0 _S102;
    (&_S102)->primal_0 = x_4;
    (&_S102)->differential_0 = 0.0f;
    DiffPair_float_0 _S103;
    (&_S103)->primal_0 = y_2;
    (&_S103)->differential_0 = 0.0f;
    _d_max_0(&_S102, &_S103, _S100.differential_0);
    float _S104 = _S91 + _S98;
    float3  _S105 = make_float3 (_S95.differential_0 + _S102.differential_0 + _S104, _S96.differential_0 + _S103.differential_0 + _S104, _S94.differential_0 + _S101.differential_0 + _S104);
    float3  _S106 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S107;
    (&_S107)->primal_0 = _S58;
    (&_S107)->differential_0 = _S106;
    s_bwd_prop_exp_1(&_S107, _S105);
    float3  _S108 = make_float3 (2.0f) * _S107.differential_0;
    float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * v_loss_0[int(2)];
    DiffPair_float_0 _S109;
    (&_S109)->primal_0 = _S56;
    (&_S109)->differential_0 = 0.0f;
    DiffPair_float_0 _S110;
    (&_S110)->primal_0 = max_gauss_ratio_1;
    (&_S110)->differential_0 = 0.0f;
    _d_max_0(&_S109, &_S110, s_diff_scale_reg_T_0);
    float _S111 = _S109.differential_0 / _S57;
    float _S112 = _S53 * - _S111;
    float _S113 = _S55 * _S111;
    DiffPair_float_0 _S114;
    (&_S114)->primal_0 = _S54;
    (&_S114)->differential_0 = 0.0f;
    DiffPair_float_0 _S115;
    (&_S115)->primal_0 = _S51;
    (&_S115)->differential_0 = 0.0f;
    _d_min_0(&_S114, &_S115, _S112);
    DiffPair_float_0 _S116;
    (&_S116)->primal_0 = _S49;
    (&_S116)->differential_0 = 0.0f;
    DiffPair_float_0 _S117;
    (&_S117)->primal_0 = _S50;
    (&_S117)->differential_0 = 0.0f;
    _d_min_0(&_S116, &_S117, _S114.differential_0);
    DiffPair_float_0 _S118;
    (&_S118)->primal_0 = _S52;
    (&_S118)->differential_0 = 0.0f;
    DiffPair_float_0 _S119;
    (&_S119)->primal_0 = _S51;
    (&_S119)->differential_0 = 0.0f;
    _d_max_0(&_S118, &_S119, _S113);
    DiffPair_float_0 _S120;
    (&_S120)->primal_0 = _S49;
    (&_S120)->differential_0 = 0.0f;
    DiffPair_float_0 _S121;
    (&_S121)->primal_0 = _S50;
    (&_S121)->differential_0 = 0.0f;
    _d_max_0(&_S120, &_S121, _S118.differential_0);
    float _S122 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * v_loss_0[int(1)]);
    float3  _S123 = make_float3 (_S116.differential_0 + _S120.differential_0 + _S122, _S117.differential_0 + _S121.differential_0 + _S122, _S115.differential_0 + _S119.differential_0 + _S122);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S124;
    (&_S124)->primal_0 = scales_1;
    (&_S124)->differential_0 = _S106;
    s_bwd_prop_exp_1(&_S124, _S123);
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * v_loss_0[int(4)];
    float _S125 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S126;
    (&_S126)->primal_0 = _S47;
    (&_S126)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S126, _S125);
    float _S127 = _S126.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S128 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S129;
    (&_S129)->primal_0 = quat_1;
    (&_S129)->differential_0 = _S128;
    s_bwd_length_impl_0(&_S129, _S127);
    float _S130 = - (mcmc_opacity_reg_weight_1 * v_loss_0[int(0)] / _S46);
    DiffPair_float_0 _S131;
    (&_S131)->primal_0 = _S44;
    (&_S131)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S131, _S130);
    float _S132 = - _S131.differential_0;
    *v_scales_0 = _S108 + _S124.differential_0;
    *v_opacity_0 = _S132;
    *v_quat_0 = _S129.differential_0;
    return;
}

