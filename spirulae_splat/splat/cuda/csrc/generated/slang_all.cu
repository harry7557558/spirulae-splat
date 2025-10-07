#include "slang.cuh"

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

__device__ void _d_exp_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S1 = (F32_exp(((*dpx_0).primal_0))) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S1;
    return;
}

__device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
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

__device__ void _d_sqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S5 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_2).primal_0)))))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S5;
    return;
}

__device__ float dot_0(float4  x_0, float4  y_0)
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

__device__ float length_0(float4  x_1)
{
    return (F32_sqrt((dot_0(x_1, x_1))));
}

__device__ void _d_log_0(DiffPair_float_0 * dpx_3, float dOut_3)
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

__device__ float3  exp_0(float3  x_2)
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

__device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, float3  dOut_4)
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

__device__ float2  exp_1(float2  x_3)
{
    float2  result_3;
    int i_2 = int(0);
    for(;;)
    {
        if(i_2 < int(2))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_3, i_2) = (F32_exp((_slang_vector_get_element(x_3, i_2))));
        i_2 = i_2 + int(1);
    }
    return result_3;
}

__device__ void _d_exp_vector_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_5, float2  dOut_5)
{
    float2  _S8 = exp_1((*dpx_5).primal_0) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S8;
    return;
}

__device__ void _d_min_0(DiffPair_float_0 * dpx_6, DiffPair_float_0 * dpy_1, float dOut_6)
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

__device__ void per_splat_losses(bool is_3dgs_0, float3  scales_0, float opacity_0, float4  quat_0, float mcmc_opacity_reg_weight_0, float mcmc_scale_reg_weight_0, float max_gauss_ratio_0, float scale_regularization_weight_0, float erank_reg_weight_0, float erank_reg_weight_s3_0, float quat_norm_reg_weight_0, FixedArray<float, 5>  * _S12)
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
        float x_4 = _S17.x;
        float y_1 = _S17.y;
        float z_0 = _S17.z;
        float s_0 = x_4 + y_1 + z_0;
        float s1_0 = (F32_max(((F32_max((x_4), (y_1)))), (z_0))) / s_0;
        float s3_0 = (F32_min(((F32_min((x_4), (y_1)))), (z_0))) / s_0;
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
        float x_5 = _S22.x;
        float y_2 = _S22.y;
        float s_1 = x_5 + y_2;
        float s1_1 = (F32_max((x_5), (y_2))) / s_1;
        float s2_1 = (F32_min((x_5), (y_2))) / s_1;
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

__device__ float s_primal_ctx_exp_0(float _S23)
{
    return (F32_exp((_S23)));
}

__device__ float2  s_primal_ctx_exp_1(float2  _S24)
{
    return exp_1(_S24);
}

__device__ float s_primal_ctx_max_0(float _S25, float _S26)
{
    return (F32_max((_S25), (_S26)));
}

__device__ float s_primal_ctx_min_0(float _S27, float _S28)
{
    return (F32_min((_S27), (_S28)));
}

__device__ float s_primal_ctx_log_0(float _S29)
{
    return (F32_log((_S29)));
}

__device__ float3  s_primal_ctx_exp_2(float3  _S30)
{
    return exp_0(_S30);
}

__device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S31, DiffPair_float_0 * _S32, float _S33)
{
    _d_max_0(_S31, _S32, _S33);
    return;
}

__device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S34, float _S35)
{
    _d_log_0(_S34, _S35);
    return;
}

__device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S36, float _S37)
{
    _d_exp_0(_S36, _S37);
    return;
}

__device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S38, DiffPair_float_0 * _S39, float _S40)
{
    _d_min_0(_S38, _S39, _S40);
    return;
}

__device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S41, float2  _S42)
{
    _d_exp_vector_1(_S41, _S42);
    return;
}

__device__ void s_bwd_prop_exp_2(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S43, float3  _S44)
{
    _d_exp_vector_0(_S43, _S44);
    return;
}

__device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S45, float _S46)
{
    _d_sqrt_0(_S45, _S46);
    return;
}

__device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_7, float _s_dOut_0)
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

__device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S61, float _S62)
{
    s_bwd_prop_length_impl_0(_S61, _S62);
    return;
}

__device__ void s_bwd_prop_per_splat_losses_0(bool is_3dgs_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscales_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, float mcmc_opacity_reg_weight_1, float mcmc_scale_reg_weight_1, float max_gauss_ratio_1, float scale_regularization_weight_1, float erank_reg_weight_1, float erank_reg_weight_s3_1, float quat_norm_reg_weight_1, FixedArray<float, 5>  _s_dOut_1)
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
    float _S119;
    float _S120;
    float2  _S121;
    float2  _S122;
    float3  _S123;
    if(is_3dgs_1)
    {
        float3  _S124 = s_primal_ctx_exp_2(_S63.primal_0);
        float _S125 = _S124.x;
        float _S126 = _S124.y;
        float _S127 = _S124.z;
        float _S128 = s_primal_ctx_max_0(_S125, _S126);
        float _S129 = s_primal_ctx_max_0(_S128, _S127);
        float _S130 = s_primal_ctx_min_0(_S125, _S126);
        float _S131 = s_primal_ctx_min_0(_S130, _S127);
        float _S132 = _S129 / _S131;
        float _S133 = _S131 * _S131;
        float3  _S134 = make_float3 (2.0f) * _S63.primal_0;
        float3  _S135 = s_primal_ctx_exp_2(_S134);
        float x_6 = _S135.x;
        float y_3 = _S135.y;
        float z_1 = _S135.z;
        float s_2 = x_6 + y_3 + z_1;
        float _S136 = s_primal_ctx_max_0(x_6, y_3);
        float _S137 = s_primal_ctx_max_0(_S136, z_1);
        float s1_2 = _S137 / s_2;
        float _S138 = s_2 * s_2;
        float _S139 = s_primal_ctx_min_0(x_6, y_3);
        float _S140 = s_primal_ctx_min_0(_S139, z_1);
        float s3_1 = _S140 / s_2;
        float s2_2 = 1.0f - s1_2 - s3_1;
        float _S141 = - s1_2;
        float _S142 = s_primal_ctx_log_0(s1_2);
        float _S143 = s_primal_ctx_log_0(s2_2);
        float _S144 = s_primal_ctx_log_0(s3_1);
        float _S145 = _S141 * _S142 - s2_2 * _S143 - s3_1 * _S144;
        float _S146 = s_primal_ctx_exp_0(_S145) - 0.99998998641967773f;
        float _S147 = - s_primal_ctx_log_0(_S146);
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
        _S85 = 0.0f;
        _S121 = _S65;
        _S86 = 0.0f;
        _S87 = 0.0f;
        _S88 = 0.0f;
        _S89 = 0.0f;
        _S90 = 0.0f;
        _S91 = 0.0f;
        _S122 = _S65;
        _S92 = _S147;
        _S93 = _S146;
        _S94 = _S145;
        _S95 = s3_1;
        _S96 = _S144;
        _S97 = s2_2;
        _S98 = _S143;
        _S99 = _S141;
        _S100 = _S142;
        _S101 = s1_2;
        _S102 = _S138;
        _S103 = _S140;
        _S104 = s_2;
        _S105 = _S139;
        _S106 = z_1;
        _S107 = x_6;
        _S108 = y_3;
        _S109 = _S138;
        _S110 = _S137;
        _S111 = _S136;
        _S123 = _S134;
        _S112 = _S132;
        _S113 = _S133;
        _S114 = _S129;
        _S115 = _S131;
        _S116 = _S130;
        _S117 = _S127;
        _S118 = _S125;
        _S119 = _S126;
        _S120 = _S128;
    }
    else
    {
        float2  _S148 = float2 {_S63.primal_0.x, _S63.primal_0.y};
        float2  _S149 = s_primal_ctx_exp_1(_S148);
        float _S150 = _S149.x;
        float _S151 = _S149.y;
        float _S152 = s_primal_ctx_max_0(_S150, _S151);
        float _S153 = s_primal_ctx_min_0(_S150, _S151);
        float _S154 = _S152 / _S153;
        float _S155 = _S153 * _S153;
        float2  _S156 = make_float2 (2.0f) * _S148;
        float2  _S157 = s_primal_ctx_exp_1(_S156);
        float x_7 = _S157.x;
        float y_4 = _S157.y;
        float s_3 = x_7 + y_4;
        float _S158 = s_primal_ctx_max_0(x_7, y_4);
        float s1_3 = _S158 / s_3;
        float _S159 = s_3 * s_3;
        float _S160 = s_primal_ctx_min_0(x_7, y_4);
        float s2_3 = _S160 / s_3;
        float _S161 = - s1_3;
        float _S162 = s_primal_ctx_log_0(s1_3);
        float _S163 = s_primal_ctx_log_0(s2_3);
        float _S164 = _S161 * _S162 - s2_3 * _S163;
        float _S165 = s_primal_ctx_exp_0(_S164) - 0.99998998641967773f;
        _S71 = - s_primal_ctx_log_0(_S165);
        _S72 = _S165;
        _S73 = _S164;
        _S74 = s2_3;
        _S75 = _S163;
        _S76 = _S161;
        _S77 = _S162;
        _S78 = s1_3;
        _S79 = _S159;
        _S80 = _S160;
        _S81 = s_3;
        _S82 = x_7;
        _S83 = y_4;
        _S84 = _S159;
        _S85 = _S158;
        _S121 = _S156;
        _S86 = _S154;
        _S87 = _S155;
        _S88 = _S152;
        _S89 = _S153;
        _S90 = _S150;
        _S91 = _S151;
        _S122 = _S148;
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
        _S110 = 0.0f;
        _S111 = 0.0f;
        _S123 = _S66;
        _S112 = 0.0f;
        _S113 = 0.0f;
        _S114 = 0.0f;
        _S115 = 0.0f;
        _S116 = 0.0f;
        _S117 = 0.0f;
        _S118 = 0.0f;
        _S119 = 0.0f;
        _S120 = 0.0f;
    }
    if(is_3dgs_1)
    {
        float _S166 = erank_reg_weight_s3_1 * _s_dOut_1[int(3)];
        float _S167 = erank_reg_weight_1 * _s_dOut_1[int(3)];
        DiffPair_float_0 _S168;
        (&_S168)->primal_0 = _S92;
        (&_S168)->differential_0 = 0.0f;
        DiffPair_float_0 _S169;
        (&_S169)->primal_0 = 0.0f;
        (&_S169)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S168, &_S169, _S167);
        float _S170 = - _S168.differential_0;
        DiffPair_float_0 _S171;
        (&_S171)->primal_0 = _S93;
        (&_S171)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S171, _S170);
        DiffPair_float_0 _S172;
        (&_S172)->primal_0 = _S94;
        (&_S172)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S172, _S171.differential_0);
        float _S173 = - _S172.differential_0;
        float _S174 = _S95 * _S173;
        float _S175 = _S96 * _S173;
        DiffPair_float_0 _S176;
        (&_S176)->primal_0 = _S95;
        (&_S176)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S176, _S174);
        float _S177 = _S97 * _S173;
        float _S178 = _S98 * _S173;
        DiffPair_float_0 _S179;
        (&_S179)->primal_0 = _S97;
        (&_S179)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S179, _S177);
        float _S180 = _S99 * _S172.differential_0;
        float _S181 = _S100 * _S172.differential_0;
        DiffPair_float_0 _S182;
        (&_S182)->primal_0 = _S101;
        (&_S182)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S182, _S180);
        float _S183 = - _S181;
        float _S184 = - (_S178 + _S179.differential_0);
        float _S185 = (_S166 + _S175 + _S176.differential_0 + _S184) / _S102;
        float _S186 = _S103 * - _S185;
        float _S187 = _S104 * _S185;
        DiffPair_float_0 _S188;
        (&_S188)->primal_0 = _S105;
        (&_S188)->differential_0 = 0.0f;
        DiffPair_float_0 _S189;
        (&_S189)->primal_0 = _S106;
        (&_S189)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S188, &_S189, _S187);
        DiffPair_float_0 _S190;
        (&_S190)->primal_0 = _S107;
        (&_S190)->differential_0 = 0.0f;
        DiffPair_float_0 _S191;
        (&_S191)->primal_0 = _S108;
        (&_S191)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S190, &_S191, _S188.differential_0);
        float _S192 = (_S182.differential_0 + _S183 + _S184) / _S109;
        float _S193 = _S110 * - _S192;
        float _S194 = _S104 * _S192;
        DiffPair_float_0 _S195;
        (&_S195)->primal_0 = _S111;
        (&_S195)->differential_0 = 0.0f;
        DiffPair_float_0 _S196;
        (&_S196)->primal_0 = _S106;
        (&_S196)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S195, &_S196, _S194);
        DiffPair_float_0 _S197;
        (&_S197)->primal_0 = _S107;
        (&_S197)->differential_0 = 0.0f;
        DiffPair_float_0 _S198;
        (&_S198)->primal_0 = _S108;
        (&_S198)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S197, &_S198, _S195.differential_0);
        float _S199 = _S186 + _S193;
        float3  _S200 = make_float3 (_S190.differential_0 + _S197.differential_0 + _S199, _S191.differential_0 + _S198.differential_0 + _S199, _S189.differential_0 + _S196.differential_0 + _S199);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S201;
        (&_S201)->primal_0 = _S123;
        (&_S201)->differential_0 = _S66;
        s_bwd_prop_exp_2(&_S201, _S200);
        float3  _S202 = make_float3 (2.0f) * _S201.differential_0;
        float s_diff_scale_reg_T_0 = scale_regularization_weight_1 * _s_dOut_1[int(2)];
        DiffPair_float_0 _S203;
        (&_S203)->primal_0 = _S112;
        (&_S203)->differential_0 = 0.0f;
        DiffPair_float_0 _S204;
        (&_S204)->primal_0 = max_gauss_ratio_1;
        (&_S204)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S203, &_S204, s_diff_scale_reg_T_0);
        float _S205 = _S203.differential_0 / _S113;
        float _S206 = _S114 * - _S205;
        float _S207 = _S115 * _S205;
        DiffPair_float_0 _S208;
        (&_S208)->primal_0 = _S116;
        (&_S208)->differential_0 = 0.0f;
        DiffPair_float_0 _S209;
        (&_S209)->primal_0 = _S117;
        (&_S209)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S208, &_S209, _S206);
        DiffPair_float_0 _S210;
        (&_S210)->primal_0 = _S118;
        (&_S210)->differential_0 = 0.0f;
        DiffPair_float_0 _S211;
        (&_S211)->primal_0 = _S119;
        (&_S211)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S210, &_S211, _S208.differential_0);
        DiffPair_float_0 _S212;
        (&_S212)->primal_0 = _S120;
        (&_S212)->differential_0 = 0.0f;
        DiffPair_float_0 _S213;
        (&_S213)->primal_0 = _S117;
        (&_S213)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S212, &_S213, _S207);
        DiffPair_float_0 _S214;
        (&_S214)->primal_0 = _S118;
        (&_S214)->differential_0 = 0.0f;
        DiffPair_float_0 _S215;
        (&_S215)->primal_0 = _S119;
        (&_S215)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S214, &_S215, _S212.differential_0);
        float _S216 = mcmc_scale_reg_weight_1 * (0.3333333432674408f * _s_dOut_1[int(1)]);
        float3  _S217 = make_float3 (_S210.differential_0 + _S214.differential_0 + _S216, _S211.differential_0 + _S215.differential_0 + _S216, _S209.differential_0 + _S213.differential_0 + _S216);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S218;
        (&_S218)->primal_0 = _S63.primal_0;
        (&_S218)->differential_0 = _S66;
        s_bwd_prop_exp_2(&_S218, _S217);
        float3  _S219 = _S202 + _S218.differential_0;
        _S71 = _s_dOut_1[int(4)];
        _S72 = _s_dOut_1[int(0)];
        _S123 = _S219;
    }
    else
    {
        float _S220 = erank_reg_weight_1 * _s_dOut_1[int(3)];
        DiffPair_float_0 _S221;
        (&_S221)->primal_0 = _S71;
        (&_S221)->differential_0 = 0.0f;
        DiffPair_float_0 _S222;
        (&_S222)->primal_0 = 0.0f;
        (&_S222)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S221, &_S222, _S220);
        float _S223 = - _S221.differential_0;
        DiffPair_float_0 _S224;
        (&_S224)->primal_0 = _S72;
        (&_S224)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S224, _S223);
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = _S73;
        (&_S225)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S225, _S224.differential_0);
        float _S226 = - _S225.differential_0;
        float _S227 = _S74 * _S226;
        float _S228 = _S75 * _S226;
        DiffPair_float_0 _S229;
        (&_S229)->primal_0 = _S74;
        (&_S229)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S229, _S227);
        float _S230 = _S76 * _S225.differential_0;
        float _S231 = _S77 * _S225.differential_0;
        DiffPair_float_0 _S232;
        (&_S232)->primal_0 = _S78;
        (&_S232)->differential_0 = 0.0f;
        s_bwd_prop_log_0(&_S232, _S230);
        float _S233 = - _S231;
        float _S234 = (_S228 + _S229.differential_0) / _S79;
        float _S235 = _S80 * - _S234;
        float _S236 = _S81 * _S234;
        DiffPair_float_0 _S237;
        (&_S237)->primal_0 = _S82;
        (&_S237)->differential_0 = 0.0f;
        DiffPair_float_0 _S238;
        (&_S238)->primal_0 = _S83;
        (&_S238)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S237, &_S238, _S236);
        float _S239 = (_S232.differential_0 + _S233) / _S84;
        float _S240 = _S85 * - _S239;
        float _S241 = _S81 * _S239;
        DiffPair_float_0 _S242;
        (&_S242)->primal_0 = _S82;
        (&_S242)->differential_0 = 0.0f;
        DiffPair_float_0 _S243;
        (&_S243)->primal_0 = _S83;
        (&_S243)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S242, &_S243, _S241);
        float _S244 = _S235 + _S240;
        float2  _S245 = make_float2 (_S237.differential_0 + _S242.differential_0 + _S244, _S238.differential_0 + _S243.differential_0 + _S244);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S246;
        (&_S246)->primal_0 = _S121;
        (&_S246)->differential_0 = _S65;
        s_bwd_prop_exp_1(&_S246, _S245);
        float2  _S247 = make_float2 (2.0f) * _S246.differential_0;
        float s_diff_scale_reg_T_1 = scale_regularization_weight_1 * _s_dOut_1[int(2)];
        DiffPair_float_0 _S248;
        (&_S248)->primal_0 = _S86;
        (&_S248)->differential_0 = 0.0f;
        DiffPair_float_0 _S249;
        (&_S249)->primal_0 = max_gauss_ratio_1;
        (&_S249)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S248, &_S249, s_diff_scale_reg_T_1);
        float _S250 = _S248.differential_0 / _S87;
        float _S251 = _S88 * - _S250;
        float _S252 = _S89 * _S250;
        DiffPair_float_0 _S253;
        (&_S253)->primal_0 = _S90;
        (&_S253)->differential_0 = 0.0f;
        DiffPair_float_0 _S254;
        (&_S254)->primal_0 = _S91;
        (&_S254)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S253, &_S254, _S251);
        DiffPair_float_0 _S255;
        (&_S255)->primal_0 = _S90;
        (&_S255)->differential_0 = 0.0f;
        DiffPair_float_0 _S256;
        (&_S256)->primal_0 = _S91;
        (&_S256)->differential_0 = 0.0f;
        s_bwd_prop_max_0(&_S255, &_S256, _S252);
        float _S257 = mcmc_scale_reg_weight_1 * (0.5f * _s_dOut_1[int(1)]);
        float2  _S258 = make_float2 (_S253.differential_0 + _S255.differential_0 + _S257, _S254.differential_0 + _S256.differential_0 + _S257);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S259;
        (&_S259)->primal_0 = _S122;
        (&_S259)->differential_0 = _S65;
        s_bwd_prop_exp_1(&_S259, _S258);
        float2  _S260 = _S247 + _S259.differential_0;
        float3  _S261 = make_float3 (_S260.x, _S260.y, 0.0f);
        _S71 = _s_dOut_1[int(4)];
        _S72 = _s_dOut_1[int(0)];
        _S123 = _S261;
    }
    float s_diff_quat_norm_reg_T_0 = quat_norm_reg_weight_1 * _S71;
    float _S262 = - s_diff_quat_norm_reg_T_0;
    DiffPair_float_0 _S263;
    (&_S263)->primal_0 = _S70;
    (&_S263)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S263, _S262);
    float _S264 = _S263.differential_0 + s_diff_quat_norm_reg_T_0;
    float4  _S265 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S266;
    (&_S266)->primal_0 = _S64.primal_0;
    (&_S266)->differential_0 = _S265;
    s_bwd_length_impl_0(&_S266, _S264);
    float _S267 = - (mcmc_opacity_reg_weight_1 * _S72 / _S69);
    DiffPair_float_0 _S268;
    (&_S268)->primal_0 = _S67;
    (&_S268)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S268, _S267);
    float _S269 = - _S268.differential_0;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S266.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S269;
    dpscales_0->primal_0 = (*dpscales_0).primal_0;
    dpscales_0->differential_0 = _S123;
    return;
}

__device__ void s_bwd_per_splat_losses_0(bool _S270, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S271, DiffPair_float_0 * _S272, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S273, float _S274, float _S275, float _S276, float _S277, float _S278, float _S279, float _S280, FixedArray<float, 5>  _S281)
{
    s_bwd_prop_per_splat_losses_0(_S270, _S271, _S272, _S273, _S274, _S275, _S276, _S277, _S278, _S279, _S280, _S281);
    return;
}

__device__ void per_splat_losses_bwd(bool is_3dgs_2, float3  scales_1, float opacity_1, float4  quat_1, FixedArray<float, 5>  v_loss_0, float3  * v_scales_0, float * v_opacity_0, float4  * v_quat_0, float mcmc_opacity_reg_weight_2, float mcmc_scale_reg_weight_2, float max_gauss_ratio_2, float scale_regularization_weight_2, float erank_reg_weight_2, float erank_reg_weight_s3_2, float quat_norm_reg_weight_2)
{
    float3  _S282 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_scales_0;
    (&p_scales_0)->primal_0 = scales_1;
    (&p_scales_0)->differential_0 = _S282;
    DiffPair_float_0 p_opacity_0;
    (&p_opacity_0)->primal_0 = opacity_1;
    (&p_opacity_0)->differential_0 = 0.0f;
    float4  _S283 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 p_quat_0;
    (&p_quat_0)->primal_0 = quat_1;
    (&p_quat_0)->differential_0 = _S283;
    s_bwd_per_splat_losses_0(is_3dgs_2, &p_scales_0, &p_opacity_0, &p_quat_0, mcmc_opacity_reg_weight_2, mcmc_scale_reg_weight_2, max_gauss_ratio_2, scale_regularization_weight_2, erank_reg_weight_2, erank_reg_weight_s3_2, quat_norm_reg_weight_2, v_loss_0);
    *v_scales_0 = p_scales_0.differential_0;
    *v_opacity_0 = p_opacity_0.differential_0;
    *v_quat_0 = p_quat_0.differential_0;
    return;
}

__device__ float3  min_0(float3  x_8, float3  y_5)
{
    float3  result_4;
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
        *_slang_vector_get_element_ptr(&result_4, i_3) = (F32_min((_slang_vector_get_element(x_8, i_3)), (_slang_vector_get_element(y_5, i_3))));
        i_3 = i_3 + int(1);
    }
    return result_4;
}

__device__ float3  max_0(float3  x_9, float3  y_6)
{
    float3  result_5;
    int i_4 = int(0);
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_5, i_4) = (F32_max((_slang_vector_get_element(x_9, i_4)), (_slang_vector_get_element(y_6, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_5;
}

__device__ void _d_clamp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_7)
{
    DiffPair_float_0 _S284 = *dpx_8;
    bool _S285;
    if(((*dpx_8).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S285 = ((*dpx_8).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S285 = false;
    }
    float _S286;
    if(_S285)
    {
        _S286 = dOut_7;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpx_8->primal_0 = _S284.primal_0;
    dpx_8->differential_0 = _S286;
    DiffPair_float_0 _S287 = *dpMin_0;
    if((_S284.primal_0) < ((*dpMin_0).primal_0))
    {
        _S286 = dOut_7;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpMin_0->primal_0 = _S287.primal_0;
    dpMin_0->differential_0 = _S286;
    DiffPair_float_0 _S288 = *dpMax_0;
    if(((*dpx_8).primal_0) > ((*dpMax_0).primal_0))
    {
        _S286 = dOut_7;
    }
    else
    {
        _S286 = 0.0f;
    }
    dpMax_0->primal_0 = _S288.primal_0;
    dpMax_0->differential_0 = _S286;
    return;
}

__device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_8)
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

__device__ float3  clamp_0(float3  x_10, float3  minBound_0, float3  maxBound_0)
{
    return min_0(max_0(x_10, minBound_0), maxBound_0);
}

__device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_0(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

__device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S289, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S290, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S291, float3  _S292)
{
    _d_clamp_vector_0(_S289, _S290, _S291, _S292);
    return;
}

__device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_2)
{
    float _S293 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S294 = make_float3 (_S293);
    float3  _S295 = make_float3 (0.0f);
    float3  _S296 = make_float3 (1.0f);
    float3  _S297 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S298;
    (&_S298)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S293) * (*dpbackground_0).primal_0;
    (&_S298)->differential_0 = _S297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S299;
    (&_S299)->primal_0 = _S295;
    (&_S299)->differential_0 = _S297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S300;
    (&_S300)->primal_0 = _S296;
    (&_S300)->differential_0 = _S297;
    s_bwd_prop_clamp_0(&_S298, &_S299, &_S300, _s_dOut_2);
    float3  _S301 = _S294 * _S298.differential_0;
    float3  _S302 = (*dpbackground_0).primal_0 * _S298.differential_0;
    float _S303 = - (_S302.x + _S302.y + _S302.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S301;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S303;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S298.differential_0;
    return;
}

__device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S304, DiffPair_float_0 * _S305, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S306, float3  _S307)
{
    s_bwd_prop_blend_background_0(_S304, _S305, _S306, _S307);
    return;
}

__device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S308 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S308;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S308;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

