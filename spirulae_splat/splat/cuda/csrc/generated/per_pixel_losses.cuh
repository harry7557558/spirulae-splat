#pragma once

#include "slang.cuh"

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_0, float3  dOut_0)
{
    float3  _S1 = _slang_select(((*dpx_0).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_0).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S1;
    return;
}

inline __device__ float3  abs_0(float3  x_0)
{
    float3  result_0;
    int i_0 = int(0);
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_0, i_0) = (F32_abs((_slang_vector_get_element(x_0, i_0))));
        i_0 = i_0 + int(1);
    }
    return result_0;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, float dOut_1)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_0).primal_0.x * dOut_1;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_1).primal_0.x * dOut_1;
    *&((&x_d_result_0)->y) = (*dpy_0).primal_0.y * dOut_1;
    *&((&y_d_result_0)->y) = (*dpx_1).primal_0.y * dOut_1;
    *&((&x_d_result_0)->z) = (*dpy_0).primal_0.z * dOut_1;
    *&((&y_d_result_0)->z) = (*dpx_1).primal_0.z * dOut_1;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = x_d_result_0;
    dpy_0->primal_0 = (*dpy_0).primal_0;
    dpy_0->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_1, float3  y_0)
{
    int i_1 = int(0);
    float result_1 = 0.0f;
    for(;;)
    {
        if(i_1 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_2 = result_1 + _slang_vector_get_element(x_1, i_1) * _slang_vector_get_element(y_0, i_1);
        i_1 = i_1 + int(1);
        result_1 = result_2;
    }
    return result_1;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpy_1, float dOut_2)
{
    DiffPair_float_0 _S2 = *dpx_2;
    float _S3;
    if(((*dpx_2).primal_0) > ((*dpy_1).primal_0))
    {
        _S3 = dOut_2;
    }
    else
    {
        if(((*dpx_2).primal_0) < ((*dpy_1).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_2;
        }
    }
    dpx_2->primal_0 = _S2.primal_0;
    dpx_2->differential_0 = _S3;
    DiffPair_float_0 _S4 = *dpy_1;
    if(((*dpy_1).primal_0) > (_S2.primal_0))
    {
        _S3 = dOut_2;
    }
    else
    {
        if(((*dpy_1).primal_0) < ((*dpx_2).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_2;
        }
    }
    dpy_1->primal_0 = _S4.primal_0;
    dpy_1->differential_0 = _S3;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_3)
{
    DiffPair_float_0 _S5 = *dpx_3;
    bool _S6;
    if(((*dpx_3).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S6 = ((*dpx_3).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S6 = false;
    }
    float _S7;
    if(_S6)
    {
        _S7 = dOut_3;
    }
    else
    {
        _S7 = 0.0f;
    }
    dpx_3->primal_0 = _S5.primal_0;
    dpx_3->differential_0 = _S7;
    DiffPair_float_0 _S8 = *dpMin_0;
    if((_S5.primal_0) < ((*dpMin_0).primal_0))
    {
        _S7 = dOut_3;
    }
    else
    {
        _S7 = 0.0f;
    }
    dpMin_0->primal_0 = _S8.primal_0;
    dpMin_0->differential_0 = _S7;
    DiffPair_float_0 _S9 = *dpMax_0;
    if(((*dpx_3).primal_0) > ((*dpMax_0).primal_0))
    {
        _S7 = dOut_3;
    }
    else
    {
        _S7 = 0.0f;
    }
    dpMax_0->primal_0 = _S9.primal_0;
    dpMax_0->differential_0 = _S7;
    return;
}

inline __device__ float clamp_0(float x_2, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_2), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_4, float dOut_4)
{
    float _S10 = 1.0f / (*dpx_4).primal_0 * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S10;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_5, float dOut_5)
{
    float _S11 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_5).primal_0)))))) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S11;
    return;
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_6, float dOut_6)
{
    float _S12 = -0.5f / ((*dpx_6).primal_0 * (F32_sqrt(((*dpx_6).primal_0)))) * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S12;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpy_2, DiffPair_float_0 * dps_0, float dOut_7)
{
    float _S13 = (1.0f - (*dps_0).primal_0) * dOut_7;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S13;
    DiffPair_float_0 _S14 = *dpy_2;
    float _S15 = (*dps_0).primal_0 * dOut_7;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = _S15;
    float _S16 = (_S14.primal_0 - (*dpx_7).primal_0) * dOut_7;
    dps_0->primal_0 = _S14.primal_0;
    dps_0->differential_0 = _S16;
    return;
}

inline __device__ float lerp_0(float x_3, float y_1, float s_0)
{
    return x_3 + (y_1 - x_3) * s_0;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_Ts_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool mask_0, bool depth_mask_0, bool normal_mask_0, bool alpha_mask_0, FixedArray<float, 11>  weights_0, FixedArray<float, 23>  * _S17)
{
    float3  _S18;
    bool _S19;
    bool _S20;
    FixedArray<float, 23>  losses_0;
    float _S21 = float(mask_0);
    float3  _S22 = ref_rgb_0 - render_rgb_0;
    float3  _S23 = abs_0(_S22);
    float _S24 = dot_0(_S22, _S22) * 0.3333333432674408f;
    losses_0[int(0)] = _S21 * (weights_0[int(0)] * ((_S23.x + _S23.y + _S23.z) * 0.3333333432674408f) + weights_0[int(1)] * _S24);
    losses_0[int(1)] = _S21 * clamp_0(_S24, 0.0f, 1.0f);
    float _S25 = float(depth_mask_0 & mask_0);
    float _S26 = _S25 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S27 = _S25 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_0[int(2)] = _S26;
    losses_0[int(3)] = _S27;
    losses_0[int(4)] = _S26 * _S26;
    losses_0[int(5)] = _S27 * _S27;
    losses_0[int(6)] = _S26 * _S27;
    bool _S28 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S29 = norm2_0 == 0.0f;
        _S19 = _S29;
        if(_S29)
        {
            _S18 = make_float3 (0.0f);
            break;
        }
        _S18 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S30;
    bool _S31 = !_S19;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S32 = norm2_1 == 0.0f;
        _S20 = _S32;
        if(_S32)
        {
            _S30 = make_float3 (0.0f);
            break;
        }
        _S30 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    bool _S33;
    float3  _S34;
    bool _S35 = !_S20;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S34 = make_float3 (0.0f);
            _S33 = false;
            break;
        }
        _S34 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        _S33 = _S28;
        break;
    }
    float _S36 = float(_S31 & _S33);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S18, _S34);
    losses_0[int(7)] = weights_0[int(3)] * _S36 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S37 = float(_S35 & _S33);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S30, _S34);
    losses_0[int(8)] = weights_0[int(3)] * _S37 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S38 = float(_S31 & _S35);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S18, _S30);
    losses_0[int(11)] = weights_0[int(6)] * _S38 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float render_alpha_0 = clamp_0(1.0f - render_Ts_0, 0.0f, 1.0f);
    float _S39 = float(alpha_mask_0);
    float _S40 = float(ref_alpha_0);
    float _S41 = (F32_max((render_alpha_0), (_S40)));
    losses_0[int(9)] = weights_0[int(4)] * _S39 * - lerp_0((F32_log(((F32_max((1.0f - _S41), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S41), (9.99999997475242708e-07f)))))), _S40);
    float _S42 = 1.0f - render_alpha_0;
    float _S43 = 1.0f - _S40;
    float _S44 = (F32_max((_S42), (_S43)));
    losses_0[int(10)] = weights_0[int(5)] * _S39 * - lerp_0((F32_log(((F32_max((1.0f - _S44), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S44), (9.99999997475242708e-07f)))))), _S43);
    losses_0[int(12)] = weights_0[int(7)] * 4.0f * render_alpha_0 * _S42;
    float _S45 = (F32_max((render_alpha_0), (9.999999960041972e-13f)));
    losses_0[int(13)] = weights_0[int(8)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f) / _S45;
    losses_0[int(14)] = weights_0[int(9)] * depth_dist_0 / _S45;
    losses_0[int(15)] = weights_0[int(10)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f) / _S45;
    losses_0[int(16)] = 1.0f;
    losses_0[int(17)] = _S21;
    losses_0[int(18)] = _S25;
    losses_0[int(19)] = _S36;
    losses_0[int(20)] = _S37;
    losses_0[int(21)] = _S38;
    losses_0[int(22)] = _S39;
    *_S17 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S46, float3  _S47)
{
    return dot_0(_S46, _S47);
}

inline __device__ float s_primal_ctx_log_0(float _S48)
{
    return (F32_log((_S48)));
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S49)
{
    return (F32_rsqrt((_S49)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S50, float _S51, float _S52)
{
    return clamp_0(_S50, _S51, _S52);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S53, DiffPair_float_0 * _S54, DiffPair_float_0 * _S55, float _S56)
{
    _d_lerp_0(_S53, _S54, _S55, _S56);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S57, float _S58)
{
    _d_log_0(_S57, _S58);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S59, DiffPair_float_0 * _S60, DiffPair_float_0 * _S61, float _S62)
{
    _d_clamp_0(_S59, _S60, _S61, _S62);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S63, float _S64)
{
    _d_sqrt_0(_S63, _S64);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S65, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S66, float _S67)
{
    _d_dot_0(_S65, _S66, _S67);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S68, float _S69)
{
    _d_rsqrt_0(_S68, _S69);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S70, float3  _S71)
{
    _d_abs_vector_0(_S70, _S71);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_Ts_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 11>  * weights_1, FixedArray<float, 23>  * _s_dOut_0)
{
    DiffPair_float_0 _S72 = *dprender_depth_0;
    DiffPair_float_0 _S73 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S74 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S75 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S76 = *dpref_normal_0;
    DiffPair_float_0 _S77 = *dprender_Ts_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S78 = *dprgb_dist_0;
    DiffPair_float_0 _S79 = *dpdepth_dist_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S80 = *dpnormal_dist_0;
    float3  _S81 = make_float3 (0.0f);
    float _S82 = float(mask_1);
    float _S83 = (*weights_1)[int(0)];
    float3  _S84 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S85 = (*weights_1)[int(1)];
    float _S86 = s_primal_ctx_dot_0(_S84, _S84) * 0.3333333432674408f;
    float _S87 = float(depth_mask_1 & mask_1);
    float _S88 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S89 = _S87 * s_primal_ctx_log_0(_S88);
    float _S90 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S91 = _S87 * s_primal_ctx_log_0(_S90);
    bool _S92 = normal_mask_1 & mask_1;
    float _S93 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S94 = _S93 == 0.0f;
    float3  _S95;
    if(_S94)
    {
        _S95 = make_float3 (0.0f);
    }
    bool _S96 = !_S94;
    float3  _S97;
    if(_S96)
    {
        float _S98 = s_primal_ctx_rsqrt_0(_S93);
        float3  _S99 = make_float3 (_S98);
        _S95 = _S74.primal_0 * make_float3 (_S98);
        _S97 = _S99;
    }
    else
    {
        _S97 = _S81;
    }
    float _S100 = s_primal_ctx_dot_0(_S75.primal_0, _S75.primal_0);
    bool _S101 = _S100 == 0.0f;
    float3  _S102;
    if(_S101)
    {
        _S102 = make_float3 (0.0f);
    }
    bool _S103 = !_S101;
    float3  _S104;
    if(_S103)
    {
        float _S105 = s_primal_ctx_rsqrt_0(_S100);
        float3  _S106 = make_float3 (_S105);
        _S102 = _S75.primal_0 * make_float3 (_S105);
        _S104 = _S106;
    }
    else
    {
        _S104 = _S81;
    }
    float _S107 = s_primal_ctx_dot_0(_S76.primal_0, _S76.primal_0);
    bool _S108 = _S107 == 0.0f;
    float3  _S109;
    bool _S110;
    if(_S108)
    {
        float3  _S111 = make_float3 (0.0f);
        _S110 = false;
        _S109 = _S111;
    }
    else
    {
        _S110 = _S92;
    }
    bool _S112 = !_S108;
    float3  _S113;
    if(_S112)
    {
        float _S114 = s_primal_ctx_rsqrt_0(_S107);
        float3  _S115 = make_float3 (_S114);
        _S109 = _S76.primal_0 * make_float3 (_S114);
        _S113 = _S115;
    }
    else
    {
        _S113 = _S81;
    }
    float _S116 = (*weights_1)[int(3)] * float(_S96 & _S110);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S95, _S109);
    float _S117 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S118 = (*weights_1)[int(3)] * float(_S103 & _S110);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S102, _S109);
    float _S119 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S120 = (*weights_1)[int(6)] * float(_S96 & _S103);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S95, _S102);
    float _S121 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S122 = 1.0f - _S77.primal_0;
    float _S123 = s_primal_ctx_clamp_0(_S122, 0.0f, 1.0f);
    float _S124 = float(alpha_mask_1);
    float _S125 = (*weights_1)[int(4)] * _S124;
    float _S126 = float(ref_alpha_1);
    float _S127 = (F32_max((_S123), (_S126)));
    float _S128 = 1.0f - _S127;
    float _S129 = (F32_max((_S128), (9.99999997475242708e-07f)));
    float _S130 = s_primal_ctx_log_0(_S129);
    float _S131 = (F32_max((_S127), (9.99999997475242708e-07f)));
    float _S132 = s_primal_ctx_log_0(_S131);
    float _S133 = (*weights_1)[int(5)] * _S124;
    float _S134 = 1.0f - _S123;
    float _S135 = 1.0f - _S126;
    float _S136 = (F32_max((_S134), (_S135)));
    float _S137 = 1.0f - _S136;
    float _S138 = (F32_max((_S137), (9.99999997475242708e-07f)));
    float _S139 = s_primal_ctx_log_0(_S138);
    float _S140 = (F32_max((_S136), (9.99999997475242708e-07f)));
    float _S141 = s_primal_ctx_log_0(_S140);
    float _S142 = (*weights_1)[int(7)] * 4.0f;
    float _S143 = _S142 * _S123;
    float _S144 = (F32_max((_S123), (9.999999960041972e-13f)));
    float _S145 = _S144 * _S144;
    float _S146 = (*_s_dOut_0)[int(0)];
    float _S147 = (*_s_dOut_0)[int(1)];
    float _S148 = (*_s_dOut_0)[int(2)];
    float _S149 = (*_s_dOut_0)[int(3)];
    float _S150 = (*_s_dOut_0)[int(4)];
    float _S151 = (*_s_dOut_0)[int(5)];
    float _S152 = (*_s_dOut_0)[int(6)];
    float _S153 = (*_s_dOut_0)[int(15)] / _S145;
    float _S154 = 0.3333333432674408f * ((*weights_1)[int(10)] * (_S144 * _S153));
    float _S155 = (*_s_dOut_0)[int(14)] / _S145;
    float _S156 = (*weights_1)[int(9)] * (_S144 * _S155);
    float _S157 = (*_s_dOut_0)[int(13)] / _S145;
    float _S158 = _S144 * _S157;
    float _S159 = (*weights_1)[int(10)] * ((_S80.primal_0.x + _S80.primal_0.y + _S80.primal_0.z) * 0.3333333432674408f) * - _S153 + (*weights_1)[int(9)] * _S79.primal_0 * - _S155 + (*weights_1)[int(8)] * ((_S78.primal_0.x + _S78.primal_0.y + _S78.primal_0.z) * 0.3333333432674408f) * - _S157;
    DiffPair_float_0 _S160;
    (&_S160)->primal_0 = _S123;
    (&_S160)->differential_0 = 0.0f;
    DiffPair_float_0 _S161;
    (&_S161)->primal_0 = 9.999999960041972e-13f;
    (&_S161)->differential_0 = 0.0f;
    _d_max_0(&_S160, &_S161, _S159);
    float _S162 = 0.3333333432674408f * ((*weights_1)[int(8)] * _S158);
    float _S163 = _S143 * (*_s_dOut_0)[int(12)];
    float _S164 = _S142 * (_S134 * (*_s_dOut_0)[int(12)]);
    float _S165 = - (_S133 * (*_s_dOut_0)[int(10)]);
    DiffPair_float_0 _S166;
    (&_S166)->primal_0 = _S139;
    (&_S166)->differential_0 = 0.0f;
    DiffPair_float_0 _S167;
    (&_S167)->primal_0 = _S141;
    (&_S167)->differential_0 = 0.0f;
    DiffPair_float_0 _S168;
    (&_S168)->primal_0 = _S135;
    (&_S168)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S166, &_S167, &_S168, _S165);
    DiffPair_float_0 _S169;
    (&_S169)->primal_0 = _S140;
    (&_S169)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S169, _S167.differential_0);
    DiffPair_float_0 _S170;
    (&_S170)->primal_0 = _S136;
    (&_S170)->differential_0 = 0.0f;
    DiffPair_float_0 _S171;
    (&_S171)->primal_0 = 9.99999997475242708e-07f;
    (&_S171)->differential_0 = 0.0f;
    _d_max_0(&_S170, &_S171, _S169.differential_0);
    DiffPair_float_0 _S172;
    (&_S172)->primal_0 = _S138;
    (&_S172)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S172, _S166.differential_0);
    DiffPair_float_0 _S173;
    (&_S173)->primal_0 = _S137;
    (&_S173)->differential_0 = 0.0f;
    DiffPair_float_0 _S174;
    (&_S174)->primal_0 = 9.99999997475242708e-07f;
    (&_S174)->differential_0 = 0.0f;
    _d_max_0(&_S173, &_S174, _S172.differential_0);
    float _S175 = _S170.differential_0 + - _S173.differential_0;
    DiffPair_float_0 _S176;
    (&_S176)->primal_0 = _S134;
    (&_S176)->differential_0 = 0.0f;
    DiffPair_float_0 _S177;
    (&_S177)->primal_0 = _S135;
    (&_S177)->differential_0 = 0.0f;
    _d_max_0(&_S176, &_S177, _S175);
    float _S178 = - (_S163 + _S176.differential_0);
    float _S179 = - (_S125 * (*_s_dOut_0)[int(9)]);
    DiffPair_float_0 _S180;
    (&_S180)->primal_0 = _S130;
    (&_S180)->differential_0 = 0.0f;
    DiffPair_float_0 _S181;
    (&_S181)->primal_0 = _S132;
    (&_S181)->differential_0 = 0.0f;
    DiffPair_float_0 _S182;
    (&_S182)->primal_0 = _S126;
    (&_S182)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S180, &_S181, &_S182, _S179);
    DiffPair_float_0 _S183;
    (&_S183)->primal_0 = _S131;
    (&_S183)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S183, _S181.differential_0);
    DiffPair_float_0 _S184;
    (&_S184)->primal_0 = _S127;
    (&_S184)->differential_0 = 0.0f;
    DiffPair_float_0 _S185;
    (&_S185)->primal_0 = 9.99999997475242708e-07f;
    (&_S185)->differential_0 = 0.0f;
    _d_max_0(&_S184, &_S185, _S183.differential_0);
    DiffPair_float_0 _S186;
    (&_S186)->primal_0 = _S129;
    (&_S186)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S186, _S180.differential_0);
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = _S128;
    (&_S187)->differential_0 = 0.0f;
    DiffPair_float_0 _S188;
    (&_S188)->primal_0 = 9.99999997475242708e-07f;
    (&_S188)->differential_0 = 0.0f;
    _d_max_0(&_S187, &_S188, _S186.differential_0);
    float _S189 = _S184.differential_0 + - _S187.differential_0;
    DiffPair_float_0 _S190;
    (&_S190)->primal_0 = _S123;
    (&_S190)->differential_0 = 0.0f;
    DiffPair_float_0 _S191;
    (&_S191)->primal_0 = _S126;
    (&_S191)->differential_0 = 0.0f;
    _d_max_0(&_S190, &_S191, _S189);
    float _S192 = _S160.differential_0 + _S164 + _S178 + _S190.differential_0;
    DiffPair_float_0 _S193;
    (&_S193)->primal_0 = _S122;
    (&_S193)->differential_0 = 0.0f;
    DiffPair_float_0 _S194;
    (&_S194)->primal_0 = 0.0f;
    (&_S194)->differential_0 = 0.0f;
    DiffPair_float_0 _S195;
    (&_S195)->primal_0 = 1.0f;
    (&_S195)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S193, &_S194, &_S195, _S192);
    float _S196 = - _S193.differential_0;
    float _S197 = _S120 * (*_s_dOut_0)[int(11)];
    DiffPair_float_0 _S198;
    (&_S198)->primal_0 = _S121;
    (&_S198)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S198, _S197);
    DiffPair_float_0 _S199;
    (&_S199)->primal_0 = cos_sim_loss_5;
    (&_S199)->differential_0 = 0.0f;
    DiffPair_float_0 _S200;
    (&_S200)->primal_0 = 9.999999960041972e-13f;
    (&_S200)->differential_0 = 0.0f;
    _d_max_0(&_S199, &_S200, _S198.differential_0);
    float _S201 = 0.5f * - (_S197 + _S199.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S202;
    (&_S202)->primal_0 = _S95;
    (&_S202)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S203;
    (&_S203)->primal_0 = _S102;
    (&_S203)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S202, &_S203, _S201);
    float _S204 = _S118 * (*_s_dOut_0)[int(8)];
    DiffPair_float_0 _S205;
    (&_S205)->primal_0 = _S119;
    (&_S205)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S205, _S204);
    DiffPair_float_0 _S206;
    (&_S206)->primal_0 = cos_sim_loss_4;
    (&_S206)->differential_0 = 0.0f;
    DiffPair_float_0 _S207;
    (&_S207)->primal_0 = 9.999999960041972e-13f;
    (&_S207)->differential_0 = 0.0f;
    _d_max_0(&_S206, &_S207, _S205.differential_0);
    float _S208 = 0.5f * - (_S204 + _S206.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S209;
    (&_S209)->primal_0 = _S102;
    (&_S209)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S210;
    (&_S210)->primal_0 = _S109;
    (&_S210)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S209, &_S210, _S208);
    float _S211 = _S116 * (*_s_dOut_0)[int(7)];
    DiffPair_float_0 _S212;
    (&_S212)->primal_0 = _S117;
    (&_S212)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S212, _S211);
    DiffPair_float_0 _S213;
    (&_S213)->primal_0 = cos_sim_loss_3;
    (&_S213)->differential_0 = 0.0f;
    DiffPair_float_0 _S214;
    (&_S214)->primal_0 = 9.999999960041972e-13f;
    (&_S214)->differential_0 = 0.0f;
    _d_max_0(&_S213, &_S214, _S212.differential_0);
    float _S215 = 0.5f * - (_S211 + _S213.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S216;
    (&_S216)->primal_0 = _S95;
    (&_S216)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S217;
    (&_S217)->primal_0 = _S109;
    (&_S217)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S216, &_S217, _S215);
    float3  _S218 = _S210.differential_0 + _S217.differential_0;
    float3  _S219 = _S202.differential_0 + _S216.differential_0;
    float3  _S220 = make_float3 (_S154, _S154, _S154);
    float3  _S221 = make_float3 (_S162, _S162, _S162);
    float3  _S222 = _S203.differential_0 + _S209.differential_0;
    float _S223;
    if(_S112)
    {
        float3  _S224 = _S76.primal_0 * _S218;
        float3  _S225 = _S113 * _S218;
        float _S226 = _S224.x + _S224.y + _S224.z;
        DiffPair_float_0 _S227;
        (&_S227)->primal_0 = _S107;
        (&_S227)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S227, _S226);
        _S223 = _S227.differential_0;
        _S95 = _S225;
    }
    else
    {
        _S223 = 0.0f;
        _S95 = _S81;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S228;
    (&_S228)->primal_0 = _S76.primal_0;
    (&_S228)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S229;
    (&_S229)->primal_0 = _S76.primal_0;
    (&_S229)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S228, &_S229, _S223);
    float3  _S230 = _S229.differential_0 + _S228.differential_0 + _S95;
    if(_S103)
    {
        float3  _S231 = _S75.primal_0 * _S222;
        float3  _S232 = _S104 * _S222;
        float _S233 = _S231.x + _S231.y + _S231.z;
        DiffPair_float_0 _S234;
        (&_S234)->primal_0 = _S100;
        (&_S234)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S234, _S233);
        _S223 = _S234.differential_0;
        _S95 = _S232;
    }
    else
    {
        _S223 = 0.0f;
        _S95 = _S81;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S235;
    (&_S235)->primal_0 = _S75.primal_0;
    (&_S235)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S236;
    (&_S236)->primal_0 = _S75.primal_0;
    (&_S236)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S235, &_S236, _S223);
    float3  _S237 = _S236.differential_0 + _S235.differential_0 + _S95;
    if(_S96)
    {
        float3  _S238 = _S74.primal_0 * _S219;
        float3  _S239 = _S97 * _S219;
        float _S240 = _S238.x + _S238.y + _S238.z;
        DiffPair_float_0 _S241;
        (&_S241)->primal_0 = _S93;
        (&_S241)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S241, _S240);
        _S223 = _S241.differential_0;
        _S95 = _S239;
    }
    else
    {
        _S223 = 0.0f;
        _S95 = _S81;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S74.primal_0;
    (&_S242)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S243;
    (&_S243)->primal_0 = _S74.primal_0;
    (&_S243)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S242, &_S243, _S223);
    float _S244 = _S91 * _S152;
    float _S245 = _S91 * _S151;
    float _S246 = _S89 * _S150;
    float _S247 = _S87 * (_S89 * _S152 + _S245 + _S245 + _S149);
    DiffPair_float_0 _S248;
    (&_S248)->primal_0 = _S90;
    (&_S248)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S248, _S247);
    DiffPair_float_0 _S249;
    (&_S249)->primal_0 = _S73.primal_0;
    (&_S249)->differential_0 = 0.0f;
    DiffPair_float_0 _S250;
    (&_S250)->primal_0 = 0.00009999999747379f;
    (&_S250)->differential_0 = 0.0f;
    _d_max_0(&_S249, &_S250, _S248.differential_0);
    float _S251 = _S87 * (_S244 + _S246 + _S246 + _S148);
    DiffPair_float_0 _S252;
    (&_S252)->primal_0 = _S88;
    (&_S252)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S252, _S251);
    DiffPair_float_0 _S253;
    (&_S253)->primal_0 = _S72.primal_0;
    (&_S253)->differential_0 = 0.0f;
    DiffPair_float_0 _S254;
    (&_S254)->primal_0 = 0.00009999999747379f;
    (&_S254)->differential_0 = 0.0f;
    _d_max_0(&_S253, &_S254, _S252.differential_0);
    float _S255 = _S82 * _S147;
    DiffPair_float_0 _S256;
    (&_S256)->primal_0 = _S86;
    (&_S256)->differential_0 = 0.0f;
    DiffPair_float_0 _S257;
    (&_S257)->primal_0 = 0.0f;
    (&_S257)->differential_0 = 0.0f;
    DiffPair_float_0 _S258;
    (&_S258)->primal_0 = 1.0f;
    (&_S258)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S256, &_S257, &_S258, _S255);
    float _S259 = _S82 * _S146;
    float _S260 = 0.3333333432674408f * (_S256.differential_0 + _S85 * _S259);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S261;
    (&_S261)->primal_0 = _S84;
    (&_S261)->differential_0 = _S81;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S262;
    (&_S262)->primal_0 = _S84;
    (&_S262)->differential_0 = _S81;
    s_bwd_prop_dot_0(&_S261, &_S262, _S260);
    float _S263 = 0.3333333432674408f * (_S83 * _S259);
    float3  _S264 = make_float3 (_S263, _S263, _S263);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S265;
    (&_S265)->primal_0 = _S84;
    (&_S265)->differential_0 = _S81;
    s_bwd_prop_abs_0(&_S265, _S264);
    float3  _S266 = _S262.differential_0 + _S261.differential_0 + _S265.differential_0;
    float3  _S267 = - _S266;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S220;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S156;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S221;
    dprender_Ts_0->primal_0 = (*dprender_Ts_0).primal_0;
    dprender_Ts_0->differential_0 = _S196;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S230;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S237;
    float3  _S268 = _S243.differential_0 + _S242.differential_0 + _S95;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S268;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S249.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S253.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S266;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S267;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S269, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S270, DiffPair_float_0 * _S271, DiffPair_float_0 * _S272, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S273, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S274, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S275, DiffPair_float_0 * _S276, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S277, DiffPair_float_0 * _S278, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S279, bool _S280, bool _S281, bool _S282, bool _S283, bool _S284, FixedArray<float, 11>  * _S285, FixedArray<float, 23>  * _S286)
{
    s_bwd_prop_per_pixel_losses_0(_S269, _S270, _S271, _S272, _S273, _S274, _S275, _S276, _S277, _S278, _S279, _S280, _S281, _S282, _S283, _S284, _S285, _S286);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_Ts_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 11>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_Ts_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S287 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S287;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S287;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S287;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S287;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S287;
    DiffPair_float_0 dp_render_Ts_0;
    (&dp_render_Ts_0)->primal_0 = render_Ts_1;
    (&dp_render_Ts_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S287;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S287;
    FixedArray<float, 11>  _S288 = weights_2;
    FixedArray<float, 23>  _S289 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_Ts_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S288, &_S289);
    *v_render_rgb_0 = dp_render_rgb_0.differential_0;
    *v_ref_rgb_0 = dp_ref_rgb_0.differential_0;
    *v_render_depth_0 = dp_render_depth_0.differential_0;
    *v_ref_depth_0 = dp_ref_depth_0.differential_0;
    *v_render_normal_0 = dp_render_normal_0.differential_0;
    *v_depth_normal_0 = dp_depth_normal_0.differential_0;
    *v_ref_normal_0 = dp_ref_normal_0.differential_0;
    *v_render_Ts_0 = dp_render_Ts_0.differential_0;
    *v_rgb_dist_0 = dp_rgb_dist_0.differential_0;
    *v_depth_dist_0 = dp_depth_dist_0.differential_0;
    *v_normal_dist_0 = dp_normal_dist_0.differential_0;
    return;
}

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_8, float dOut_8)
{
    float _S290 = 1.0f / ((*dpx_8).primal_0 * 2.30258512496948242f) * dOut_8;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S290;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 11>  weights_3, FixedArray<float, 10>  * _S291)
{
    FixedArray<float, 10>  losses_1;
    float _S292 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_1[int(0)] = raw_losses_0[int(0)] / _S292;
    losses_1[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S292)));
    bool _S293;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S293 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S293 = false;
    }
    float _S294;
    if(_S293)
    {
        _S294 = weights_3[int(2)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S294 = 0.0f;
    }
    losses_1[int(2)] = _S294;
    losses_1[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_1[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_1[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S295 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_1[int(6)] = raw_losses_0[int(12)] / _S295;
    losses_1[int(7)] = raw_losses_0[int(13)] / _S295;
    losses_1[int(8)] = raw_losses_0[int(14)] / _S295;
    losses_1[int(9)] = raw_losses_0[int(15)] / _S295;
    *_S291 = losses_1;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S296)
{
    return (F32_sqrt((_S296)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S297, float _S298)
{
    _d_log10_0(_S297, _S298);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 11>  * weights_4, FixedArray<float, 10>  * _s_dOut_1)
{
    FixedArray<float, 23>  _S299 = dpraw_losses_0->primal_0;
    float _S300 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S301 = _S300 * _S300;
    float _S302 = dpraw_losses_0->primal_0[int(1)] / _S300;
    bool _S303 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S304;
    if(_S303)
    {
        _S304 = (_S299[int(3)]) != 0.0f;
    }
    else
    {
        _S304 = false;
    }
    float _S305;
    float _S306;
    float _S307;
    float _S308;
    float _S309;
    float _S310;
    float _S311;
    float _S312;
    float _S313;
    float _S314;
    float _S315;
    float _S316;
    float _S317;
    float _S318;
    float _S319;
    if(_S304)
    {
        float _S320 = _S299[int(2)] * _S299[int(3)];
        float _S321 = _S299[int(18)] * _S299[int(18)];
        float _S322 = _S299[int(6)] - _S320 / _S299[int(18)];
        float _S323 = _S299[int(2)] * _S299[int(2)];
        float _S324 = _S299[int(4)] - _S323 / _S299[int(18)];
        float _S325 = _S299[int(3)] * _S299[int(3)];
        float _S326 = _S299[int(5)] - _S325 / _S299[int(18)];
        float _S327 = _S324 * _S326 + 1.0f;
        float _S328 = (F32_max((9.999999960041972e-13f), (_S327)));
        float _S329 = s_primal_ctx_sqrt_0(_S328);
        float _S330 = _S329 * _S329;
        float _S331 = 1.0f - _S322 / _S329;
        _S305 = (*weights_4)[int(2)];
        _S306 = _S331;
        _S307 = _S330;
        _S308 = _S322;
        _S309 = _S329;
        _S310 = _S328;
        _S311 = _S327;
        _S312 = _S324;
        _S313 = _S326;
        _S314 = _S321;
        _S315 = _S325;
        _S316 = _S299[int(3)];
        _S317 = _S323;
        _S318 = _S299[int(2)];
        _S319 = _S320;
    }
    else
    {
        _S305 = 0.0f;
        _S306 = 0.0f;
        _S307 = 0.0f;
        _S308 = 0.0f;
        _S309 = 0.0f;
        _S310 = 0.0f;
        _S311 = 0.0f;
        _S312 = 0.0f;
        _S313 = 0.0f;
        _S314 = 0.0f;
        _S315 = 0.0f;
        _S316 = 0.0f;
        _S317 = 0.0f;
        _S318 = 0.0f;
        _S319 = 0.0f;
    }
    float _S332 = (F32_max((_S299[int(19)]), (1.0f)));
    float _S333 = _S332 * _S332;
    float _S334 = (F32_max((_S299[int(20)]), (1.0f)));
    float _S335 = _S334 * _S334;
    float _S336 = float((I32_max((int((_S299[int(19)]) > 0.5f) + int((_S299[int(20)]) > 0.5f)), (int(1)))));
    float _S337 = _S299[int(9)] + _S299[int(10)];
    float _S338 = (F32_max((_S299[int(22)]), (1.0f)));
    float _S339 = _S338 * _S338;
    float _S340 = (F32_max((_S299[int(21)]), (1.0f)));
    float _S341 = _S340 * _S340;
    float _S342 = (F32_max((_S299[int(16)]), (1.0f)));
    float _S343 = _S342 * _S342;
    float _S344 = (*_s_dOut_1)[int(0)];
    float _S345 = (*_s_dOut_1)[int(1)];
    float _S346 = (*_s_dOut_1)[int(2)];
    float _S347 = (*_s_dOut_1)[int(9)] / _S343;
    float _S348 = _S342 * _S347;
    float _S349 = (*_s_dOut_1)[int(8)] / _S343;
    float _S350 = _S342 * _S349;
    float _S351 = (*_s_dOut_1)[int(7)] / _S343;
    float _S352 = _S342 * _S351;
    float _S353 = (*_s_dOut_1)[int(6)] / _S343;
    float _S354 = _S342 * _S353;
    float _S355 = _S299[int(15)] * - _S347 + _S299[int(14)] * - _S349 + _S299[int(13)] * - _S351 + _S299[int(12)] * - _S353;
    DiffPair_float_0 _S356;
    (&_S356)->primal_0 = _S299[int(16)];
    (&_S356)->differential_0 = 0.0f;
    DiffPair_float_0 _S357;
    (&_S357)->primal_0 = 1.0f;
    (&_S357)->differential_0 = 0.0f;
    _d_max_0(&_S356, &_S357, _S355);
    float _S358 = (*_s_dOut_1)[int(5)] / _S341;
    float _S359 = _S299[int(11)] * - _S358;
    float _S360 = _S340 * _S358;
    DiffPair_float_0 _S361;
    (&_S361)->primal_0 = _S299[int(21)];
    (&_S361)->differential_0 = 0.0f;
    DiffPair_float_0 _S362;
    (&_S362)->primal_0 = 1.0f;
    (&_S362)->differential_0 = 0.0f;
    _d_max_0(&_S361, &_S362, _S359);
    float _S363 = (*_s_dOut_1)[int(4)] / _S339;
    float _S364 = _S337 * - _S363;
    float _S365 = _S338 * _S363;
    DiffPair_float_0 _S366;
    (&_S366)->primal_0 = _S299[int(22)];
    (&_S366)->differential_0 = 0.0f;
    DiffPair_float_0 _S367;
    (&_S367)->primal_0 = 1.0f;
    (&_S367)->differential_0 = 0.0f;
    _d_max_0(&_S366, &_S367, _S364);
    float _S368 = (*_s_dOut_1)[int(3)] / _S336;
    float _S369 = _S368 / _S335;
    float _S370 = _S299[int(8)] * - _S369;
    float _S371 = _S334 * _S369;
    DiffPair_float_0 _S372;
    (&_S372)->primal_0 = _S299[int(20)];
    (&_S372)->differential_0 = 0.0f;
    DiffPair_float_0 _S373;
    (&_S373)->primal_0 = 1.0f;
    (&_S373)->differential_0 = 0.0f;
    _d_max_0(&_S372, &_S373, _S370);
    float _S374 = _S368 / _S333;
    float _S375 = _S299[int(7)] * - _S374;
    float _S376 = _S332 * _S374;
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = _S299[int(19)];
    (&_S377)->differential_0 = 0.0f;
    DiffPair_float_0 _S378;
    (&_S378)->primal_0 = 1.0f;
    (&_S378)->differential_0 = 0.0f;
    _d_max_0(&_S377, &_S378, _S375);
    FixedArray<float, 23>  _S379;
    _S379[int(0)] = 0.0f;
    _S379[int(1)] = 0.0f;
    _S379[int(2)] = 0.0f;
    _S379[int(3)] = 0.0f;
    _S379[int(4)] = 0.0f;
    _S379[int(5)] = 0.0f;
    _S379[int(6)] = 0.0f;
    _S379[int(7)] = 0.0f;
    _S379[int(8)] = 0.0f;
    _S379[int(9)] = 0.0f;
    _S379[int(10)] = 0.0f;
    _S379[int(11)] = 0.0f;
    _S379[int(12)] = 0.0f;
    _S379[int(13)] = 0.0f;
    _S379[int(14)] = 0.0f;
    _S379[int(15)] = 0.0f;
    _S379[int(16)] = 0.0f;
    _S379[int(17)] = 0.0f;
    _S379[int(18)] = 0.0f;
    _S379[int(19)] = 0.0f;
    _S379[int(20)] = 0.0f;
    _S379[int(21)] = 0.0f;
    _S379[int(22)] = 0.0f;
    _S379[int(15)] = _S348;
    _S379[int(14)] = _S350;
    _S379[int(13)] = _S352;
    _S379[int(16)] = _S356.differential_0;
    _S379[int(12)] = _S354;
    _S379[int(21)] = _S361.differential_0;
    _S379[int(11)] = _S360;
    _S379[int(22)] = _S366.differential_0;
    _S379[int(10)] = _S365;
    _S379[int(9)] = _S365;
    _S379[int(20)] = _S372.differential_0;
    _S379[int(8)] = _S371;
    _S379[int(19)] = _S377.differential_0;
    _S379[int(7)] = _S376;
    float _S380 = _S379[int(0)];
    float _S381 = _S379[int(1)];
    float _S382 = _S379[int(2)];
    float _S383 = _S379[int(3)];
    float _S384 = _S379[int(4)];
    float _S385 = _S379[int(5)];
    float _S386 = _S379[int(6)];
    float _S387 = _S379[int(7)];
    float _S388 = _S379[int(8)];
    float _S389 = _S379[int(9)];
    float _S390 = _S379[int(10)];
    float _S391 = _S379[int(11)];
    float _S392 = _S379[int(12)];
    float _S393 = _S379[int(13)];
    float _S394 = _S379[int(14)];
    float _S395 = _S379[int(15)];
    float _S396 = _S379[int(16)];
    float _S397 = _S379[int(17)];
    float _S398 = _S379[int(18)];
    float _S399 = _S379[int(19)];
    float _S400 = _S379[int(20)];
    float _S401 = _S379[int(21)];
    float _S402 = _S379[int(22)];
    FixedArray<float, 23>  _S403;
    if(_S304)
    {
        float _S404 = _S305 * _S346;
        DiffPair_float_0 _S405;
        (&_S405)->primal_0 = _S306;
        (&_S405)->differential_0 = 0.0f;
        DiffPair_float_0 _S406;
        (&_S406)->primal_0 = 0.0f;
        (&_S406)->differential_0 = 0.0f;
        DiffPair_float_0 _S407;
        (&_S407)->primal_0 = 2.0f;
        (&_S407)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S405, &_S406, &_S407, _S404);
        float _S408 = - _S405.differential_0 / _S307;
        float _S409 = _S308 * - _S408;
        float _S410 = _S309 * _S408;
        DiffPair_float_0 _S411;
        (&_S411)->primal_0 = _S310;
        (&_S411)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S411, _S409);
        DiffPair_float_0 _S412;
        (&_S412)->primal_0 = 9.999999960041972e-13f;
        (&_S412)->differential_0 = 0.0f;
        DiffPair_float_0 _S413;
        (&_S413)->primal_0 = _S311;
        (&_S413)->differential_0 = 0.0f;
        _d_max_0(&_S412, &_S413, _S411.differential_0);
        float _S414 = _S312 * _S413.differential_0;
        float _S415 = _S313 * _S413.differential_0;
        float _S416 = - _S414 / _S314;
        float _S417 = _S316 * (_S299[int(18)] * _S416);
        float _S418 = - _S415 / _S314;
        float _S419 = _S318 * (_S299[int(18)] * _S418);
        float _S420 = - _S410 / _S314;
        float _S421 = _S299[int(18)] * _S420;
        float _S422 = _S417 + _S417 + _S318 * _S421;
        float _S423 = _S419 + _S419 + _S316 * _S421;
        float _S424 = _S315 * - _S416 + _S317 * - _S418 + _S319 * - _S420;
        FixedArray<float, 23>  _S425;
        _S425[int(0)] = 0.0f;
        _S425[int(1)] = 0.0f;
        _S425[int(2)] = 0.0f;
        _S425[int(3)] = 0.0f;
        _S425[int(4)] = 0.0f;
        _S425[int(5)] = 0.0f;
        _S425[int(6)] = 0.0f;
        _S425[int(7)] = 0.0f;
        _S425[int(8)] = 0.0f;
        _S425[int(9)] = 0.0f;
        _S425[int(10)] = 0.0f;
        _S425[int(11)] = 0.0f;
        _S425[int(12)] = 0.0f;
        _S425[int(13)] = 0.0f;
        _S425[int(14)] = 0.0f;
        _S425[int(15)] = 0.0f;
        _S425[int(16)] = 0.0f;
        _S425[int(17)] = 0.0f;
        _S425[int(18)] = 0.0f;
        _S425[int(19)] = 0.0f;
        _S425[int(20)] = 0.0f;
        _S425[int(21)] = 0.0f;
        _S425[int(22)] = 0.0f;
        _S425[int(5)] = _S414;
        _S425[int(4)] = _S415;
        _S425[int(3)] = _S422;
        _S425[int(2)] = _S423;
        _S425[int(6)] = _S410;
        float _S426 = _S381 + _S425[int(1)];
        float _S427 = _S382 + _S425[int(2)];
        float _S428 = _S383 + _S425[int(3)];
        float _S429 = _S384 + _S425[int(4)];
        float _S430 = _S385 + _S425[int(5)];
        float _S431 = _S386 + _S425[int(6)];
        float _S432 = _S387 + _S425[int(7)];
        float _S433 = _S388 + _S425[int(8)];
        float _S434 = _S389 + _S425[int(9)];
        float _S435 = _S390 + _S425[int(10)];
        float _S436 = _S391 + _S425[int(11)];
        float _S437 = _S392 + _S425[int(12)];
        float _S438 = _S393 + _S425[int(13)];
        float _S439 = _S394 + _S425[int(14)];
        float _S440 = _S395 + _S425[int(15)];
        float _S441 = _S396 + _S425[int(16)];
        float _S442 = _S397 + _S425[int(17)];
        float _S443 = _S398 + _S425[int(18)];
        float _S444 = _S399 + _S425[int(19)];
        float _S445 = _S400 + _S425[int(20)];
        float _S446 = _S401 + _S425[int(21)];
        float _S447 = _S402 + _S425[int(22)];
        _S403[int(0)] = _S380 + _S425[int(0)];
        _S403[int(1)] = _S426;
        _S403[int(2)] = _S427;
        _S403[int(3)] = _S428;
        _S403[int(4)] = _S429;
        _S403[int(5)] = _S430;
        _S403[int(6)] = _S431;
        _S403[int(7)] = _S432;
        _S403[int(8)] = _S433;
        _S403[int(9)] = _S434;
        _S403[int(10)] = _S435;
        _S403[int(11)] = _S436;
        _S403[int(12)] = _S437;
        _S403[int(13)] = _S438;
        _S403[int(14)] = _S439;
        _S403[int(15)] = _S440;
        _S403[int(16)] = _S441;
        _S403[int(17)] = _S442;
        _S403[int(18)] = _S443;
        _S403[int(19)] = _S444;
        _S403[int(20)] = _S445;
        _S403[int(21)] = _S446;
        _S403[int(22)] = _S447;
        _S305 = _S424;
    }
    else
    {
        _S403[int(0)] = _S380;
        _S403[int(1)] = _S381;
        _S403[int(2)] = _S382;
        _S403[int(3)] = _S383;
        _S403[int(4)] = _S384;
        _S403[int(5)] = _S385;
        _S403[int(6)] = _S386;
        _S403[int(7)] = _S387;
        _S403[int(8)] = _S388;
        _S403[int(9)] = _S389;
        _S403[int(10)] = _S390;
        _S403[int(11)] = _S391;
        _S403[int(12)] = _S392;
        _S403[int(13)] = _S393;
        _S403[int(14)] = _S394;
        _S403[int(15)] = _S395;
        _S403[int(16)] = _S396;
        _S403[int(17)] = _S397;
        _S403[int(18)] = _S398;
        _S403[int(19)] = _S399;
        _S403[int(20)] = _S400;
        _S403[int(21)] = _S401;
        _S403[int(22)] = _S402;
        _S305 = 0.0f;
    }
    if(_S303)
    {
        FixedArray<float, 23>  _S448;
        _S448[int(0)] = 0.0f;
        _S448[int(1)] = 0.0f;
        _S448[int(2)] = 0.0f;
        _S448[int(3)] = 0.0f;
        _S448[int(4)] = 0.0f;
        _S448[int(5)] = 0.0f;
        _S448[int(6)] = 0.0f;
        _S448[int(7)] = 0.0f;
        _S448[int(8)] = 0.0f;
        _S448[int(9)] = 0.0f;
        _S448[int(10)] = 0.0f;
        _S448[int(11)] = 0.0f;
        _S448[int(12)] = 0.0f;
        _S448[int(13)] = 0.0f;
        _S448[int(14)] = 0.0f;
        _S448[int(15)] = 0.0f;
        _S448[int(16)] = 0.0f;
        _S448[int(17)] = 0.0f;
        _S448[int(18)] = 0.0f;
        _S448[int(19)] = 0.0f;
        _S448[int(20)] = 0.0f;
        _S448[int(21)] = 0.0f;
        _S448[int(22)] = 0.0f;
        _S448[int(3)] = 0.0f;
        float _S449 = _S403[int(1)] + _S448[int(1)];
        float _S450 = _S403[int(2)] + _S448[int(2)];
        float _S451 = _S403[int(3)] + _S448[int(3)];
        float _S452 = _S403[int(4)] + _S448[int(4)];
        float _S453 = _S403[int(5)] + _S448[int(5)];
        float _S454 = _S403[int(6)] + _S448[int(6)];
        float _S455 = _S403[int(7)] + _S448[int(7)];
        float _S456 = _S403[int(8)] + _S448[int(8)];
        float _S457 = _S403[int(9)] + _S448[int(9)];
        float _S458 = _S403[int(10)] + _S448[int(10)];
        float _S459 = _S403[int(11)] + _S448[int(11)];
        float _S460 = _S403[int(12)] + _S448[int(12)];
        float _S461 = _S403[int(13)] + _S448[int(13)];
        float _S462 = _S403[int(14)] + _S448[int(14)];
        float _S463 = _S403[int(15)] + _S448[int(15)];
        float _S464 = _S403[int(16)] + _S448[int(16)];
        float _S465 = _S403[int(17)] + _S448[int(17)];
        float _S466 = _S403[int(18)] + _S448[int(18)];
        float _S467 = _S403[int(19)] + _S448[int(19)];
        float _S468 = _S403[int(20)] + _S448[int(20)];
        float _S469 = _S403[int(21)] + _S448[int(21)];
        float _S470 = _S403[int(22)] + _S448[int(22)];
        _S403[int(0)] = _S403[int(0)] + _S448[int(0)];
        _S403[int(1)] = _S449;
        _S403[int(2)] = _S450;
        _S403[int(3)] = _S451;
        _S403[int(4)] = _S452;
        _S403[int(5)] = _S453;
        _S403[int(6)] = _S454;
        _S403[int(7)] = _S455;
        _S403[int(8)] = _S456;
        _S403[int(9)] = _S457;
        _S403[int(10)] = _S458;
        _S403[int(11)] = _S459;
        _S403[int(12)] = _S460;
        _S403[int(13)] = _S461;
        _S403[int(14)] = _S462;
        _S403[int(15)] = _S463;
        _S403[int(16)] = _S464;
        _S403[int(17)] = _S465;
        _S403[int(18)] = _S466;
        _S403[int(19)] = _S467;
        _S403[int(20)] = _S468;
        _S403[int(21)] = _S469;
        _S403[int(22)] = _S470;
    }
    float _S471 = -10.0f * _S345;
    DiffPair_float_0 _S472;
    (&_S472)->primal_0 = _S302;
    (&_S472)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S472, _S471);
    float _S473 = _S472.differential_0 / _S301;
    float _S474 = _S300 * _S473;
    float _S475 = _S344 / _S301;
    float _S476 = _S300 * _S475;
    float _S477 = _S299[int(1)] * - _S473 + _S299[int(0)] * - _S475;
    DiffPair_float_0 _S478;
    (&_S478)->primal_0 = _S299[int(17)];
    (&_S478)->differential_0 = 0.0f;
    DiffPair_float_0 _S479;
    (&_S479)->primal_0 = 1.0f;
    (&_S479)->differential_0 = 0.0f;
    _d_max_0(&_S478, &_S479, _S477);
    FixedArray<float, 23>  _S480;
    _S480[int(0)] = 0.0f;
    _S480[int(1)] = 0.0f;
    _S480[int(2)] = 0.0f;
    _S480[int(3)] = 0.0f;
    _S480[int(4)] = 0.0f;
    _S480[int(5)] = 0.0f;
    _S480[int(6)] = 0.0f;
    _S480[int(7)] = 0.0f;
    _S480[int(8)] = 0.0f;
    _S480[int(9)] = 0.0f;
    _S480[int(10)] = 0.0f;
    _S480[int(11)] = 0.0f;
    _S480[int(12)] = 0.0f;
    _S480[int(13)] = 0.0f;
    _S480[int(14)] = 0.0f;
    _S480[int(15)] = 0.0f;
    _S480[int(16)] = 0.0f;
    _S480[int(17)] = 0.0f;
    _S480[int(18)] = 0.0f;
    _S480[int(19)] = 0.0f;
    _S480[int(20)] = 0.0f;
    _S480[int(21)] = 0.0f;
    _S480[int(22)] = 0.0f;
    _S480[int(18)] = _S305;
    _S480[int(1)] = _S474;
    _S480[int(17)] = _S478.differential_0;
    _S480[int(0)] = _S476;
    FixedArray<float, 23>  _S481 = {
        _S403[int(0)] + _S480[int(0)], _S403[int(1)] + _S480[int(1)], _S403[int(2)] + _S480[int(2)], _S403[int(3)] + _S480[int(3)], _S403[int(4)] + _S480[int(4)], _S403[int(5)] + _S480[int(5)], _S403[int(6)] + _S480[int(6)], _S403[int(7)] + _S480[int(7)], _S403[int(8)] + _S480[int(8)], _S403[int(9)] + _S480[int(9)], _S403[int(10)] + _S480[int(10)], _S403[int(11)] + _S480[int(11)], _S403[int(12)] + _S480[int(12)], _S403[int(13)] + _S480[int(13)], _S403[int(14)] + _S480[int(14)], _S403[int(15)] + _S480[int(15)], _S403[int(16)] + _S480[int(16)], _S403[int(17)] + _S480[int(17)], _S403[int(18)] + _S480[int(18)], _S403[int(19)] + _S480[int(19)], _S403[int(20)] + _S480[int(20)], _S403[int(21)] + _S480[int(21)], _S403[int(22)] + _S480[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S481;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S482, FixedArray<float, 11>  * _S483, FixedArray<float, 10>  * _S484)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S482, _S483, _S484);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 11>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S485)
{
    FixedArray<float, 23>  _S486 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S486;
    FixedArray<float, 11>  _S487 = weights_5;
    FixedArray<float, 10>  _S488 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S487, &_S488);
    *_S485 = (&dp_raw_losses_0)->differential_0;
    return;
}

