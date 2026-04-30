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
    losses_0[int(13)] = weights_0[int(8)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f);
    losses_0[int(14)] = weights_0[int(9)] * depth_dist_0;
    losses_0[int(15)] = weights_0[int(10)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f);
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

inline __device__ float s_primal_ctx_dot_0(float3  _S45, float3  _S46)
{
    return dot_0(_S45, _S46);
}

inline __device__ float s_primal_ctx_log_0(float _S47)
{
    return (F32_log((_S47)));
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S48)
{
    return (F32_rsqrt((_S48)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S49, float _S50, float _S51)
{
    return clamp_0(_S49, _S50, _S51);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S52, DiffPair_float_0 * _S53, DiffPair_float_0 * _S54, float _S55)
{
    _d_lerp_0(_S52, _S53, _S54, _S55);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S56, float _S57)
{
    _d_log_0(_S56, _S57);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S58, DiffPair_float_0 * _S59, DiffPair_float_0 * _S60, float _S61)
{
    _d_clamp_0(_S58, _S59, _S60, _S61);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S62, float _S63)
{
    _d_sqrt_0(_S62, _S63);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S64, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S65, float _S66)
{
    _d_dot_0(_S64, _S65, _S66);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S67, float _S68)
{
    _d_rsqrt_0(_S67, _S68);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S69, float3  _S70)
{
    _d_abs_vector_0(_S69, _S70);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_Ts_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool mask_1, bool depth_mask_1, bool normal_mask_1, bool alpha_mask_1, FixedArray<float, 11>  * weights_1, FixedArray<float, 23>  * _s_dOut_0)
{
    DiffPair_float_0 _S71 = *dprender_depth_0;
    DiffPair_float_0 _S72 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S73 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S74 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S75 = *dpref_normal_0;
    DiffPair_float_0 _S76 = *dprender_Ts_0;
    float3  _S77 = make_float3 (0.0f);
    float _S78 = float(mask_1);
    float _S79 = (*weights_1)[int(0)];
    float3  _S80 = (*dpref_rgb_0).primal_0 - (*dprender_rgb_0).primal_0;
    float _S81 = (*weights_1)[int(1)];
    float _S82 = s_primal_ctx_dot_0(_S80, _S80) * 0.3333333432674408f;
    float _S83 = float(depth_mask_1 & mask_1);
    float _S84 = (F32_max(((*dprender_depth_0).primal_0), (0.00009999999747379f)));
    float _S85 = _S83 * s_primal_ctx_log_0(_S84);
    float _S86 = (F32_max(((*dpref_depth_0).primal_0), (0.00009999999747379f)));
    float _S87 = _S83 * s_primal_ctx_log_0(_S86);
    bool _S88 = normal_mask_1 & mask_1;
    float _S89 = s_primal_ctx_dot_0((*dprender_normal_0).primal_0, (*dprender_normal_0).primal_0);
    bool _S90 = _S89 == 0.0f;
    float3  _S91;
    if(_S90)
    {
        _S91 = make_float3 (0.0f);
    }
    bool _S92 = !_S90;
    float3  _S93;
    if(_S92)
    {
        float _S94 = s_primal_ctx_rsqrt_0(_S89);
        float3  _S95 = make_float3 (_S94);
        _S91 = _S73.primal_0 * make_float3 (_S94);
        _S93 = _S95;
    }
    else
    {
        _S93 = _S77;
    }
    float _S96 = s_primal_ctx_dot_0(_S74.primal_0, _S74.primal_0);
    bool _S97 = _S96 == 0.0f;
    float3  _S98;
    if(_S97)
    {
        _S98 = make_float3 (0.0f);
    }
    bool _S99 = !_S97;
    float3  _S100;
    if(_S99)
    {
        float _S101 = s_primal_ctx_rsqrt_0(_S96);
        float3  _S102 = make_float3 (_S101);
        _S98 = _S74.primal_0 * make_float3 (_S101);
        _S100 = _S102;
    }
    else
    {
        _S100 = _S77;
    }
    float _S103 = s_primal_ctx_dot_0(_S75.primal_0, _S75.primal_0);
    bool _S104 = _S103 == 0.0f;
    float3  _S105;
    bool _S106;
    if(_S104)
    {
        float3  _S107 = make_float3 (0.0f);
        _S106 = false;
        _S105 = _S107;
    }
    else
    {
        _S106 = _S88;
    }
    bool _S108 = !_S104;
    float3  _S109;
    if(_S108)
    {
        float _S110 = s_primal_ctx_rsqrt_0(_S103);
        float3  _S111 = make_float3 (_S110);
        _S105 = _S75.primal_0 * make_float3 (_S110);
        _S109 = _S111;
    }
    else
    {
        _S109 = _S77;
    }
    float _S112 = (*weights_1)[int(3)] * float(_S92 & _S106);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S91, _S105);
    float _S113 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S114 = (*weights_1)[int(3)] * float(_S99 & _S106);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S98, _S105);
    float _S115 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S116 = (*weights_1)[int(6)] * float(_S92 & _S99);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S91, _S98);
    float _S117 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S118 = 1.0f - _S76.primal_0;
    float _S119 = s_primal_ctx_clamp_0(_S118, 0.0f, 1.0f);
    float _S120 = float(alpha_mask_1);
    float _S121 = (*weights_1)[int(4)] * _S120;
    float _S122 = float(ref_alpha_1);
    float _S123 = (F32_max((_S119), (_S122)));
    float _S124 = 1.0f - _S123;
    float _S125 = (F32_max((_S124), (9.99999997475242708e-07f)));
    float _S126 = s_primal_ctx_log_0(_S125);
    float _S127 = (F32_max((_S123), (9.99999997475242708e-07f)));
    float _S128 = s_primal_ctx_log_0(_S127);
    float _S129 = 1.0f - _S119;
    float _S130 = 1.0f - _S122;
    float _S131 = (F32_max((_S129), (_S130)));
    float _S132 = 1.0f - _S131;
    float _S133 = (F32_max((_S132), (9.99999997475242708e-07f)));
    float _S134 = (F32_max((_S131), (9.99999997475242708e-07f)));
    float _S135 = s_primal_ctx_log_0(_S134);
    float _S136 = (*weights_1)[int(7)] * 4.0f;
    float _S137 = (*_s_dOut_0)[int(0)];
    float _S138 = (*_s_dOut_0)[int(1)];
    float _S139 = (*_s_dOut_0)[int(2)];
    float _S140 = (*_s_dOut_0)[int(3)];
    float _S141 = (*_s_dOut_0)[int(4)];
    float _S142 = (*_s_dOut_0)[int(5)];
    float _S143 = (*_s_dOut_0)[int(6)];
    float _S144 = 0.3333333432674408f * ((*weights_1)[int(10)] * (*_s_dOut_0)[int(15)]);
    float _S145 = (*weights_1)[int(9)] * (*_s_dOut_0)[int(14)];
    float _S146 = 0.3333333432674408f * ((*weights_1)[int(8)] * (*_s_dOut_0)[int(13)]);
    float _S147 = _S136 * _S119 * (*_s_dOut_0)[int(12)];
    float _S148 = _S136 * (_S129 * (*_s_dOut_0)[int(12)]);
    float _S149 = - ((*weights_1)[int(5)] * _S120 * (*_s_dOut_0)[int(10)]);
    DiffPair_float_0 _S150;
    (&_S150)->primal_0 = s_primal_ctx_log_0(_S133);
    (&_S150)->differential_0 = 0.0f;
    DiffPair_float_0 _S151;
    (&_S151)->primal_0 = _S135;
    (&_S151)->differential_0 = 0.0f;
    DiffPair_float_0 _S152;
    (&_S152)->primal_0 = _S130;
    (&_S152)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S150, &_S151, &_S152, _S149);
    DiffPair_float_0 _S153;
    (&_S153)->primal_0 = _S134;
    (&_S153)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S153, _S151.differential_0);
    DiffPair_float_0 _S154;
    (&_S154)->primal_0 = _S131;
    (&_S154)->differential_0 = 0.0f;
    DiffPair_float_0 _S155;
    (&_S155)->primal_0 = 9.99999997475242708e-07f;
    (&_S155)->differential_0 = 0.0f;
    _d_max_0(&_S154, &_S155, _S153.differential_0);
    DiffPair_float_0 _S156;
    (&_S156)->primal_0 = _S133;
    (&_S156)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S156, _S150.differential_0);
    DiffPair_float_0 _S157;
    (&_S157)->primal_0 = _S132;
    (&_S157)->differential_0 = 0.0f;
    DiffPair_float_0 _S158;
    (&_S158)->primal_0 = 9.99999997475242708e-07f;
    (&_S158)->differential_0 = 0.0f;
    _d_max_0(&_S157, &_S158, _S156.differential_0);
    float _S159 = _S154.differential_0 + - _S157.differential_0;
    DiffPair_float_0 _S160;
    (&_S160)->primal_0 = _S129;
    (&_S160)->differential_0 = 0.0f;
    DiffPair_float_0 _S161;
    (&_S161)->primal_0 = _S130;
    (&_S161)->differential_0 = 0.0f;
    _d_max_0(&_S160, &_S161, _S159);
    float _S162 = - (_S147 + _S160.differential_0);
    float _S163 = - (_S121 * (*_s_dOut_0)[int(9)]);
    DiffPair_float_0 _S164;
    (&_S164)->primal_0 = _S126;
    (&_S164)->differential_0 = 0.0f;
    DiffPair_float_0 _S165;
    (&_S165)->primal_0 = _S128;
    (&_S165)->differential_0 = 0.0f;
    DiffPair_float_0 _S166;
    (&_S166)->primal_0 = _S122;
    (&_S166)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S164, &_S165, &_S166, _S163);
    DiffPair_float_0 _S167;
    (&_S167)->primal_0 = _S127;
    (&_S167)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S167, _S165.differential_0);
    DiffPair_float_0 _S168;
    (&_S168)->primal_0 = _S123;
    (&_S168)->differential_0 = 0.0f;
    DiffPair_float_0 _S169;
    (&_S169)->primal_0 = 9.99999997475242708e-07f;
    (&_S169)->differential_0 = 0.0f;
    _d_max_0(&_S168, &_S169, _S167.differential_0);
    DiffPair_float_0 _S170;
    (&_S170)->primal_0 = _S125;
    (&_S170)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S170, _S164.differential_0);
    DiffPair_float_0 _S171;
    (&_S171)->primal_0 = _S124;
    (&_S171)->differential_0 = 0.0f;
    DiffPair_float_0 _S172;
    (&_S172)->primal_0 = 9.99999997475242708e-07f;
    (&_S172)->differential_0 = 0.0f;
    _d_max_0(&_S171, &_S172, _S170.differential_0);
    float _S173 = _S168.differential_0 + - _S171.differential_0;
    DiffPair_float_0 _S174;
    (&_S174)->primal_0 = _S119;
    (&_S174)->differential_0 = 0.0f;
    DiffPair_float_0 _S175;
    (&_S175)->primal_0 = _S122;
    (&_S175)->differential_0 = 0.0f;
    _d_max_0(&_S174, &_S175, _S173);
    float _S176 = _S148 + _S162 + _S174.differential_0;
    DiffPair_float_0 _S177;
    (&_S177)->primal_0 = _S118;
    (&_S177)->differential_0 = 0.0f;
    DiffPair_float_0 _S178;
    (&_S178)->primal_0 = 0.0f;
    (&_S178)->differential_0 = 0.0f;
    DiffPair_float_0 _S179;
    (&_S179)->primal_0 = 1.0f;
    (&_S179)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S177, &_S178, &_S179, _S176);
    float _S180 = - _S177.differential_0;
    float _S181 = _S116 * (*_s_dOut_0)[int(11)];
    DiffPair_float_0 _S182;
    (&_S182)->primal_0 = _S117;
    (&_S182)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S182, _S181);
    DiffPair_float_0 _S183;
    (&_S183)->primal_0 = cos_sim_loss_5;
    (&_S183)->differential_0 = 0.0f;
    DiffPair_float_0 _S184;
    (&_S184)->primal_0 = 9.999999960041972e-13f;
    (&_S184)->differential_0 = 0.0f;
    _d_max_0(&_S183, &_S184, _S182.differential_0);
    float _S185 = 0.5f * - (_S181 + _S183.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S186;
    (&_S186)->primal_0 = _S91;
    (&_S186)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S187;
    (&_S187)->primal_0 = _S98;
    (&_S187)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S186, &_S187, _S185);
    float _S188 = _S114 * (*_s_dOut_0)[int(8)];
    DiffPair_float_0 _S189;
    (&_S189)->primal_0 = _S115;
    (&_S189)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S189, _S188);
    DiffPair_float_0 _S190;
    (&_S190)->primal_0 = cos_sim_loss_4;
    (&_S190)->differential_0 = 0.0f;
    DiffPair_float_0 _S191;
    (&_S191)->primal_0 = 9.999999960041972e-13f;
    (&_S191)->differential_0 = 0.0f;
    _d_max_0(&_S190, &_S191, _S189.differential_0);
    float _S192 = 0.5f * - (_S188 + _S190.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S193;
    (&_S193)->primal_0 = _S98;
    (&_S193)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S194;
    (&_S194)->primal_0 = _S105;
    (&_S194)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S193, &_S194, _S192);
    float _S195 = _S112 * (*_s_dOut_0)[int(7)];
    DiffPair_float_0 _S196;
    (&_S196)->primal_0 = _S113;
    (&_S196)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S196, _S195);
    DiffPair_float_0 _S197;
    (&_S197)->primal_0 = cos_sim_loss_3;
    (&_S197)->differential_0 = 0.0f;
    DiffPair_float_0 _S198;
    (&_S198)->primal_0 = 9.999999960041972e-13f;
    (&_S198)->differential_0 = 0.0f;
    _d_max_0(&_S197, &_S198, _S196.differential_0);
    float _S199 = 0.5f * - (_S195 + _S197.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S200;
    (&_S200)->primal_0 = _S91;
    (&_S200)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S201;
    (&_S201)->primal_0 = _S105;
    (&_S201)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S200, &_S201, _S199);
    float3  _S202 = _S194.differential_0 + _S201.differential_0;
    float3  _S203 = _S186.differential_0 + _S200.differential_0;
    float3  _S204 = make_float3 (_S144, _S144, _S144);
    float3  _S205 = make_float3 (_S146, _S146, _S146);
    float3  _S206 = _S187.differential_0 + _S193.differential_0;
    float _S207;
    if(_S108)
    {
        float3  _S208 = _S75.primal_0 * _S202;
        float3  _S209 = _S109 * _S202;
        float _S210 = _S208.x + _S208.y + _S208.z;
        DiffPair_float_0 _S211;
        (&_S211)->primal_0 = _S103;
        (&_S211)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S211, _S210);
        _S207 = _S211.differential_0;
        _S91 = _S209;
    }
    else
    {
        _S207 = 0.0f;
        _S91 = _S77;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S212;
    (&_S212)->primal_0 = _S75.primal_0;
    (&_S212)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S213;
    (&_S213)->primal_0 = _S75.primal_0;
    (&_S213)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S212, &_S213, _S207);
    float3  _S214 = _S213.differential_0 + _S212.differential_0 + _S91;
    if(_S99)
    {
        float3  _S215 = _S74.primal_0 * _S206;
        float3  _S216 = _S100 * _S206;
        float _S217 = _S215.x + _S215.y + _S215.z;
        DiffPair_float_0 _S218;
        (&_S218)->primal_0 = _S96;
        (&_S218)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S218, _S217);
        _S207 = _S218.differential_0;
        _S91 = _S216;
    }
    else
    {
        _S207 = 0.0f;
        _S91 = _S77;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S219;
    (&_S219)->primal_0 = _S74.primal_0;
    (&_S219)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S220;
    (&_S220)->primal_0 = _S74.primal_0;
    (&_S220)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S219, &_S220, _S207);
    float3  _S221 = _S220.differential_0 + _S219.differential_0 + _S91;
    if(_S92)
    {
        float3  _S222 = _S73.primal_0 * _S203;
        float3  _S223 = _S93 * _S203;
        float _S224 = _S222.x + _S222.y + _S222.z;
        DiffPair_float_0 _S225;
        (&_S225)->primal_0 = _S89;
        (&_S225)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S225, _S224);
        _S207 = _S225.differential_0;
        _S91 = _S223;
    }
    else
    {
        _S207 = 0.0f;
        _S91 = _S77;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S226;
    (&_S226)->primal_0 = _S73.primal_0;
    (&_S226)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S227;
    (&_S227)->primal_0 = _S73.primal_0;
    (&_S227)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S226, &_S227, _S207);
    float _S228 = _S87 * _S143;
    float _S229 = _S87 * _S142;
    float _S230 = _S85 * _S141;
    float _S231 = _S83 * (_S85 * _S143 + _S229 + _S229 + _S140);
    DiffPair_float_0 _S232;
    (&_S232)->primal_0 = _S86;
    (&_S232)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S232, _S231);
    DiffPair_float_0 _S233;
    (&_S233)->primal_0 = _S72.primal_0;
    (&_S233)->differential_0 = 0.0f;
    DiffPair_float_0 _S234;
    (&_S234)->primal_0 = 0.00009999999747379f;
    (&_S234)->differential_0 = 0.0f;
    _d_max_0(&_S233, &_S234, _S232.differential_0);
    float _S235 = _S83 * (_S228 + _S230 + _S230 + _S139);
    DiffPair_float_0 _S236;
    (&_S236)->primal_0 = _S84;
    (&_S236)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S236, _S235);
    DiffPair_float_0 _S237;
    (&_S237)->primal_0 = _S71.primal_0;
    (&_S237)->differential_0 = 0.0f;
    DiffPair_float_0 _S238;
    (&_S238)->primal_0 = 0.00009999999747379f;
    (&_S238)->differential_0 = 0.0f;
    _d_max_0(&_S237, &_S238, _S236.differential_0);
    float _S239 = _S78 * _S138;
    DiffPair_float_0 _S240;
    (&_S240)->primal_0 = _S82;
    (&_S240)->differential_0 = 0.0f;
    DiffPair_float_0 _S241;
    (&_S241)->primal_0 = 0.0f;
    (&_S241)->differential_0 = 0.0f;
    DiffPair_float_0 _S242;
    (&_S242)->primal_0 = 1.0f;
    (&_S242)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S240, &_S241, &_S242, _S239);
    float _S243 = _S78 * _S137;
    float _S244 = 0.3333333432674408f * (_S240.differential_0 + _S81 * _S243);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S245;
    (&_S245)->primal_0 = _S80;
    (&_S245)->differential_0 = _S77;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S246;
    (&_S246)->primal_0 = _S80;
    (&_S246)->differential_0 = _S77;
    s_bwd_prop_dot_0(&_S245, &_S246, _S244);
    float _S247 = 0.3333333432674408f * (_S79 * _S243);
    float3  _S248 = make_float3 (_S247, _S247, _S247);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S249;
    (&_S249)->primal_0 = _S80;
    (&_S249)->differential_0 = _S77;
    s_bwd_prop_abs_0(&_S249, _S248);
    float3  _S250 = _S246.differential_0 + _S245.differential_0 + _S249.differential_0;
    float3  _S251 = - _S250;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S204;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S145;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S205;
    dprender_Ts_0->primal_0 = (*dprender_Ts_0).primal_0;
    dprender_Ts_0->differential_0 = _S180;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S214;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S221;
    float3  _S252 = _S227.differential_0 + _S226.differential_0 + _S91;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S252;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S233.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S237.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S250;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S251;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S253, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S254, DiffPair_float_0 * _S255, DiffPair_float_0 * _S256, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S257, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S258, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S259, DiffPair_float_0 * _S260, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S261, DiffPair_float_0 * _S262, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S263, bool _S264, bool _S265, bool _S266, bool _S267, bool _S268, FixedArray<float, 11>  * _S269, FixedArray<float, 23>  * _S270)
{
    s_bwd_prop_per_pixel_losses_0(_S253, _S254, _S255, _S256, _S257, _S258, _S259, _S260, _S261, _S262, _S263, _S264, _S265, _S266, _S267, _S268, _S269, _S270);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_Ts_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool mask_2, bool depth_mask_2, bool normal_mask_2, bool alpha_mask_2, FixedArray<float, 11>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_Ts_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S271 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S271;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S271;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S271;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S271;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S271;
    DiffPair_float_0 dp_render_Ts_0;
    (&dp_render_Ts_0)->primal_0 = render_Ts_1;
    (&dp_render_Ts_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S271;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S271;
    FixedArray<float, 11>  _S272 = weights_2;
    FixedArray<float, 23>  _S273 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_Ts_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, mask_2, depth_mask_2, normal_mask_2, alpha_mask_2, &_S272, &_S273);
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
    float _S274 = 1.0f / ((*dpx_8).primal_0 * 2.30258512496948242f) * dOut_8;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S274;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 11>  weights_3, FixedArray<float, 10>  * _S275)
{
    FixedArray<float, 10>  losses_1;
    float _S276 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_1[int(0)] = raw_losses_0[int(0)] / _S276;
    losses_1[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S276)));
    bool _S277;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S277 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S277 = false;
    }
    float _S278;
    if(_S277)
    {
        _S278 = weights_3[int(2)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S278 = 0.0f;
    }
    losses_1[int(2)] = _S278;
    losses_1[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_1[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_1[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S279 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_1[int(6)] = raw_losses_0[int(12)] / _S279;
    losses_1[int(7)] = raw_losses_0[int(13)] / _S279;
    losses_1[int(8)] = raw_losses_0[int(14)] / _S279;
    losses_1[int(9)] = raw_losses_0[int(15)] / _S279;
    *_S275 = losses_1;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S280)
{
    return (F32_sqrt((_S280)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S281, float _S282)
{
    _d_log10_0(_S281, _S282);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 11>  * weights_4, FixedArray<float, 10>  * _s_dOut_1)
{
    FixedArray<float, 23>  _S283 = dpraw_losses_0->primal_0;
    float _S284 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S285 = _S284 * _S284;
    float _S286 = dpraw_losses_0->primal_0[int(1)] / _S284;
    bool _S287 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S288;
    if(_S287)
    {
        _S288 = (_S283[int(3)]) != 0.0f;
    }
    else
    {
        _S288 = false;
    }
    float _S289;
    float _S290;
    float _S291;
    float _S292;
    float _S293;
    float _S294;
    float _S295;
    float _S296;
    float _S297;
    float _S298;
    float _S299;
    float _S300;
    float _S301;
    float _S302;
    float _S303;
    if(_S288)
    {
        float _S304 = _S283[int(2)] * _S283[int(3)];
        float _S305 = _S283[int(18)] * _S283[int(18)];
        float _S306 = _S283[int(6)] - _S304 / _S283[int(18)];
        float _S307 = _S283[int(2)] * _S283[int(2)];
        float _S308 = _S283[int(4)] - _S307 / _S283[int(18)];
        float _S309 = _S283[int(3)] * _S283[int(3)];
        float _S310 = _S283[int(5)] - _S309 / _S283[int(18)];
        float _S311 = _S308 * _S310 + 1.0f;
        float _S312 = (F32_max((9.999999960041972e-13f), (_S311)));
        float _S313 = s_primal_ctx_sqrt_0(_S312);
        float _S314 = _S313 * _S313;
        float _S315 = 1.0f - _S306 / _S313;
        _S289 = (*weights_4)[int(2)];
        _S290 = _S315;
        _S291 = _S314;
        _S292 = _S306;
        _S293 = _S313;
        _S294 = _S312;
        _S295 = _S311;
        _S296 = _S308;
        _S297 = _S310;
        _S298 = _S305;
        _S299 = _S309;
        _S300 = _S283[int(3)];
        _S301 = _S307;
        _S302 = _S283[int(2)];
        _S303 = _S304;
    }
    else
    {
        _S289 = 0.0f;
        _S290 = 0.0f;
        _S291 = 0.0f;
        _S292 = 0.0f;
        _S293 = 0.0f;
        _S294 = 0.0f;
        _S295 = 0.0f;
        _S296 = 0.0f;
        _S297 = 0.0f;
        _S298 = 0.0f;
        _S299 = 0.0f;
        _S300 = 0.0f;
        _S301 = 0.0f;
        _S302 = 0.0f;
        _S303 = 0.0f;
    }
    float _S316 = (F32_max((_S283[int(19)]), (1.0f)));
    float _S317 = _S316 * _S316;
    float _S318 = (F32_max((_S283[int(20)]), (1.0f)));
    float _S319 = _S318 * _S318;
    float _S320 = float((I32_max((int((_S283[int(19)]) > 0.5f) + int((_S283[int(20)]) > 0.5f)), (int(1)))));
    float _S321 = _S283[int(9)] + _S283[int(10)];
    float _S322 = (F32_max((_S283[int(22)]), (1.0f)));
    float _S323 = _S322 * _S322;
    float _S324 = (F32_max((_S283[int(21)]), (1.0f)));
    float _S325 = _S324 * _S324;
    float _S326 = (F32_max((_S283[int(16)]), (1.0f)));
    float _S327 = _S326 * _S326;
    float _S328 = (*_s_dOut_1)[int(0)];
    float _S329 = (*_s_dOut_1)[int(1)];
    float _S330 = (*_s_dOut_1)[int(2)];
    float _S331 = (*_s_dOut_1)[int(9)] / _S327;
    float _S332 = _S326 * _S331;
    float _S333 = (*_s_dOut_1)[int(8)] / _S327;
    float _S334 = _S326 * _S333;
    float _S335 = (*_s_dOut_1)[int(7)] / _S327;
    float _S336 = _S326 * _S335;
    float _S337 = (*_s_dOut_1)[int(6)] / _S327;
    float _S338 = _S326 * _S337;
    float _S339 = _S283[int(15)] * - _S331 + _S283[int(14)] * - _S333 + _S283[int(13)] * - _S335 + _S283[int(12)] * - _S337;
    DiffPair_float_0 _S340;
    (&_S340)->primal_0 = _S283[int(16)];
    (&_S340)->differential_0 = 0.0f;
    DiffPair_float_0 _S341;
    (&_S341)->primal_0 = 1.0f;
    (&_S341)->differential_0 = 0.0f;
    _d_max_0(&_S340, &_S341, _S339);
    float _S342 = (*_s_dOut_1)[int(5)] / _S325;
    float _S343 = _S283[int(11)] * - _S342;
    float _S344 = _S324 * _S342;
    DiffPair_float_0 _S345;
    (&_S345)->primal_0 = _S283[int(21)];
    (&_S345)->differential_0 = 0.0f;
    DiffPair_float_0 _S346;
    (&_S346)->primal_0 = 1.0f;
    (&_S346)->differential_0 = 0.0f;
    _d_max_0(&_S345, &_S346, _S343);
    float _S347 = (*_s_dOut_1)[int(4)] / _S323;
    float _S348 = _S321 * - _S347;
    float _S349 = _S322 * _S347;
    DiffPair_float_0 _S350;
    (&_S350)->primal_0 = _S283[int(22)];
    (&_S350)->differential_0 = 0.0f;
    DiffPair_float_0 _S351;
    (&_S351)->primal_0 = 1.0f;
    (&_S351)->differential_0 = 0.0f;
    _d_max_0(&_S350, &_S351, _S348);
    float _S352 = (*_s_dOut_1)[int(3)] / _S320;
    float _S353 = _S352 / _S319;
    float _S354 = _S283[int(8)] * - _S353;
    float _S355 = _S318 * _S353;
    DiffPair_float_0 _S356;
    (&_S356)->primal_0 = _S283[int(20)];
    (&_S356)->differential_0 = 0.0f;
    DiffPair_float_0 _S357;
    (&_S357)->primal_0 = 1.0f;
    (&_S357)->differential_0 = 0.0f;
    _d_max_0(&_S356, &_S357, _S354);
    float _S358 = _S352 / _S317;
    float _S359 = _S283[int(7)] * - _S358;
    float _S360 = _S316 * _S358;
    DiffPair_float_0 _S361;
    (&_S361)->primal_0 = _S283[int(19)];
    (&_S361)->differential_0 = 0.0f;
    DiffPair_float_0 _S362;
    (&_S362)->primal_0 = 1.0f;
    (&_S362)->differential_0 = 0.0f;
    _d_max_0(&_S361, &_S362, _S359);
    FixedArray<float, 23>  _S363;
    _S363[int(0)] = 0.0f;
    _S363[int(1)] = 0.0f;
    _S363[int(2)] = 0.0f;
    _S363[int(3)] = 0.0f;
    _S363[int(4)] = 0.0f;
    _S363[int(5)] = 0.0f;
    _S363[int(6)] = 0.0f;
    _S363[int(7)] = 0.0f;
    _S363[int(8)] = 0.0f;
    _S363[int(9)] = 0.0f;
    _S363[int(10)] = 0.0f;
    _S363[int(11)] = 0.0f;
    _S363[int(12)] = 0.0f;
    _S363[int(13)] = 0.0f;
    _S363[int(14)] = 0.0f;
    _S363[int(15)] = 0.0f;
    _S363[int(16)] = 0.0f;
    _S363[int(17)] = 0.0f;
    _S363[int(18)] = 0.0f;
    _S363[int(19)] = 0.0f;
    _S363[int(20)] = 0.0f;
    _S363[int(21)] = 0.0f;
    _S363[int(22)] = 0.0f;
    _S363[int(15)] = _S332;
    _S363[int(14)] = _S334;
    _S363[int(13)] = _S336;
    _S363[int(16)] = _S340.differential_0;
    _S363[int(12)] = _S338;
    _S363[int(21)] = _S345.differential_0;
    _S363[int(11)] = _S344;
    _S363[int(22)] = _S350.differential_0;
    _S363[int(10)] = _S349;
    _S363[int(9)] = _S349;
    _S363[int(20)] = _S356.differential_0;
    _S363[int(8)] = _S355;
    _S363[int(19)] = _S361.differential_0;
    _S363[int(7)] = _S360;
    float _S364 = _S363[int(0)];
    float _S365 = _S363[int(1)];
    float _S366 = _S363[int(2)];
    float _S367 = _S363[int(3)];
    float _S368 = _S363[int(4)];
    float _S369 = _S363[int(5)];
    float _S370 = _S363[int(6)];
    float _S371 = _S363[int(7)];
    float _S372 = _S363[int(8)];
    float _S373 = _S363[int(9)];
    float _S374 = _S363[int(10)];
    float _S375 = _S363[int(11)];
    float _S376 = _S363[int(12)];
    float _S377 = _S363[int(13)];
    float _S378 = _S363[int(14)];
    float _S379 = _S363[int(15)];
    float _S380 = _S363[int(16)];
    float _S381 = _S363[int(17)];
    float _S382 = _S363[int(18)];
    float _S383 = _S363[int(19)];
    float _S384 = _S363[int(20)];
    float _S385 = _S363[int(21)];
    float _S386 = _S363[int(22)];
    FixedArray<float, 23>  _S387;
    if(_S288)
    {
        float _S388 = _S289 * _S330;
        DiffPair_float_0 _S389;
        (&_S389)->primal_0 = _S290;
        (&_S389)->differential_0 = 0.0f;
        DiffPair_float_0 _S390;
        (&_S390)->primal_0 = 0.0f;
        (&_S390)->differential_0 = 0.0f;
        DiffPair_float_0 _S391;
        (&_S391)->primal_0 = 2.0f;
        (&_S391)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S389, &_S390, &_S391, _S388);
        float _S392 = - _S389.differential_0 / _S291;
        float _S393 = _S292 * - _S392;
        float _S394 = _S293 * _S392;
        DiffPair_float_0 _S395;
        (&_S395)->primal_0 = _S294;
        (&_S395)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S395, _S393);
        DiffPair_float_0 _S396;
        (&_S396)->primal_0 = 9.999999960041972e-13f;
        (&_S396)->differential_0 = 0.0f;
        DiffPair_float_0 _S397;
        (&_S397)->primal_0 = _S295;
        (&_S397)->differential_0 = 0.0f;
        _d_max_0(&_S396, &_S397, _S395.differential_0);
        float _S398 = _S296 * _S397.differential_0;
        float _S399 = _S297 * _S397.differential_0;
        float _S400 = - _S398 / _S298;
        float _S401 = _S300 * (_S283[int(18)] * _S400);
        float _S402 = - _S399 / _S298;
        float _S403 = _S302 * (_S283[int(18)] * _S402);
        float _S404 = - _S394 / _S298;
        float _S405 = _S283[int(18)] * _S404;
        float _S406 = _S401 + _S401 + _S302 * _S405;
        float _S407 = _S403 + _S403 + _S300 * _S405;
        float _S408 = _S299 * - _S400 + _S301 * - _S402 + _S303 * - _S404;
        FixedArray<float, 23>  _S409;
        _S409[int(0)] = 0.0f;
        _S409[int(1)] = 0.0f;
        _S409[int(2)] = 0.0f;
        _S409[int(3)] = 0.0f;
        _S409[int(4)] = 0.0f;
        _S409[int(5)] = 0.0f;
        _S409[int(6)] = 0.0f;
        _S409[int(7)] = 0.0f;
        _S409[int(8)] = 0.0f;
        _S409[int(9)] = 0.0f;
        _S409[int(10)] = 0.0f;
        _S409[int(11)] = 0.0f;
        _S409[int(12)] = 0.0f;
        _S409[int(13)] = 0.0f;
        _S409[int(14)] = 0.0f;
        _S409[int(15)] = 0.0f;
        _S409[int(16)] = 0.0f;
        _S409[int(17)] = 0.0f;
        _S409[int(18)] = 0.0f;
        _S409[int(19)] = 0.0f;
        _S409[int(20)] = 0.0f;
        _S409[int(21)] = 0.0f;
        _S409[int(22)] = 0.0f;
        _S409[int(5)] = _S398;
        _S409[int(4)] = _S399;
        _S409[int(3)] = _S406;
        _S409[int(2)] = _S407;
        _S409[int(6)] = _S394;
        float _S410 = _S365 + _S409[int(1)];
        float _S411 = _S366 + _S409[int(2)];
        float _S412 = _S367 + _S409[int(3)];
        float _S413 = _S368 + _S409[int(4)];
        float _S414 = _S369 + _S409[int(5)];
        float _S415 = _S370 + _S409[int(6)];
        float _S416 = _S371 + _S409[int(7)];
        float _S417 = _S372 + _S409[int(8)];
        float _S418 = _S373 + _S409[int(9)];
        float _S419 = _S374 + _S409[int(10)];
        float _S420 = _S375 + _S409[int(11)];
        float _S421 = _S376 + _S409[int(12)];
        float _S422 = _S377 + _S409[int(13)];
        float _S423 = _S378 + _S409[int(14)];
        float _S424 = _S379 + _S409[int(15)];
        float _S425 = _S380 + _S409[int(16)];
        float _S426 = _S381 + _S409[int(17)];
        float _S427 = _S382 + _S409[int(18)];
        float _S428 = _S383 + _S409[int(19)];
        float _S429 = _S384 + _S409[int(20)];
        float _S430 = _S385 + _S409[int(21)];
        float _S431 = _S386 + _S409[int(22)];
        _S387[int(0)] = _S364 + _S409[int(0)];
        _S387[int(1)] = _S410;
        _S387[int(2)] = _S411;
        _S387[int(3)] = _S412;
        _S387[int(4)] = _S413;
        _S387[int(5)] = _S414;
        _S387[int(6)] = _S415;
        _S387[int(7)] = _S416;
        _S387[int(8)] = _S417;
        _S387[int(9)] = _S418;
        _S387[int(10)] = _S419;
        _S387[int(11)] = _S420;
        _S387[int(12)] = _S421;
        _S387[int(13)] = _S422;
        _S387[int(14)] = _S423;
        _S387[int(15)] = _S424;
        _S387[int(16)] = _S425;
        _S387[int(17)] = _S426;
        _S387[int(18)] = _S427;
        _S387[int(19)] = _S428;
        _S387[int(20)] = _S429;
        _S387[int(21)] = _S430;
        _S387[int(22)] = _S431;
        _S289 = _S408;
    }
    else
    {
        _S387[int(0)] = _S364;
        _S387[int(1)] = _S365;
        _S387[int(2)] = _S366;
        _S387[int(3)] = _S367;
        _S387[int(4)] = _S368;
        _S387[int(5)] = _S369;
        _S387[int(6)] = _S370;
        _S387[int(7)] = _S371;
        _S387[int(8)] = _S372;
        _S387[int(9)] = _S373;
        _S387[int(10)] = _S374;
        _S387[int(11)] = _S375;
        _S387[int(12)] = _S376;
        _S387[int(13)] = _S377;
        _S387[int(14)] = _S378;
        _S387[int(15)] = _S379;
        _S387[int(16)] = _S380;
        _S387[int(17)] = _S381;
        _S387[int(18)] = _S382;
        _S387[int(19)] = _S383;
        _S387[int(20)] = _S384;
        _S387[int(21)] = _S385;
        _S387[int(22)] = _S386;
        _S289 = 0.0f;
    }
    if(_S287)
    {
        FixedArray<float, 23>  _S432;
        _S432[int(0)] = 0.0f;
        _S432[int(1)] = 0.0f;
        _S432[int(2)] = 0.0f;
        _S432[int(3)] = 0.0f;
        _S432[int(4)] = 0.0f;
        _S432[int(5)] = 0.0f;
        _S432[int(6)] = 0.0f;
        _S432[int(7)] = 0.0f;
        _S432[int(8)] = 0.0f;
        _S432[int(9)] = 0.0f;
        _S432[int(10)] = 0.0f;
        _S432[int(11)] = 0.0f;
        _S432[int(12)] = 0.0f;
        _S432[int(13)] = 0.0f;
        _S432[int(14)] = 0.0f;
        _S432[int(15)] = 0.0f;
        _S432[int(16)] = 0.0f;
        _S432[int(17)] = 0.0f;
        _S432[int(18)] = 0.0f;
        _S432[int(19)] = 0.0f;
        _S432[int(20)] = 0.0f;
        _S432[int(21)] = 0.0f;
        _S432[int(22)] = 0.0f;
        _S432[int(3)] = 0.0f;
        float _S433 = _S387[int(1)] + _S432[int(1)];
        float _S434 = _S387[int(2)] + _S432[int(2)];
        float _S435 = _S387[int(3)] + _S432[int(3)];
        float _S436 = _S387[int(4)] + _S432[int(4)];
        float _S437 = _S387[int(5)] + _S432[int(5)];
        float _S438 = _S387[int(6)] + _S432[int(6)];
        float _S439 = _S387[int(7)] + _S432[int(7)];
        float _S440 = _S387[int(8)] + _S432[int(8)];
        float _S441 = _S387[int(9)] + _S432[int(9)];
        float _S442 = _S387[int(10)] + _S432[int(10)];
        float _S443 = _S387[int(11)] + _S432[int(11)];
        float _S444 = _S387[int(12)] + _S432[int(12)];
        float _S445 = _S387[int(13)] + _S432[int(13)];
        float _S446 = _S387[int(14)] + _S432[int(14)];
        float _S447 = _S387[int(15)] + _S432[int(15)];
        float _S448 = _S387[int(16)] + _S432[int(16)];
        float _S449 = _S387[int(17)] + _S432[int(17)];
        float _S450 = _S387[int(18)] + _S432[int(18)];
        float _S451 = _S387[int(19)] + _S432[int(19)];
        float _S452 = _S387[int(20)] + _S432[int(20)];
        float _S453 = _S387[int(21)] + _S432[int(21)];
        float _S454 = _S387[int(22)] + _S432[int(22)];
        _S387[int(0)] = _S387[int(0)] + _S432[int(0)];
        _S387[int(1)] = _S433;
        _S387[int(2)] = _S434;
        _S387[int(3)] = _S435;
        _S387[int(4)] = _S436;
        _S387[int(5)] = _S437;
        _S387[int(6)] = _S438;
        _S387[int(7)] = _S439;
        _S387[int(8)] = _S440;
        _S387[int(9)] = _S441;
        _S387[int(10)] = _S442;
        _S387[int(11)] = _S443;
        _S387[int(12)] = _S444;
        _S387[int(13)] = _S445;
        _S387[int(14)] = _S446;
        _S387[int(15)] = _S447;
        _S387[int(16)] = _S448;
        _S387[int(17)] = _S449;
        _S387[int(18)] = _S450;
        _S387[int(19)] = _S451;
        _S387[int(20)] = _S452;
        _S387[int(21)] = _S453;
        _S387[int(22)] = _S454;
    }
    float _S455 = -10.0f * _S329;
    DiffPair_float_0 _S456;
    (&_S456)->primal_0 = _S286;
    (&_S456)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S456, _S455);
    float _S457 = _S456.differential_0 / _S285;
    float _S458 = _S284 * _S457;
    float _S459 = _S328 / _S285;
    float _S460 = _S284 * _S459;
    float _S461 = _S283[int(1)] * - _S457 + _S283[int(0)] * - _S459;
    DiffPair_float_0 _S462;
    (&_S462)->primal_0 = _S283[int(17)];
    (&_S462)->differential_0 = 0.0f;
    DiffPair_float_0 _S463;
    (&_S463)->primal_0 = 1.0f;
    (&_S463)->differential_0 = 0.0f;
    _d_max_0(&_S462, &_S463, _S461);
    FixedArray<float, 23>  _S464;
    _S464[int(0)] = 0.0f;
    _S464[int(1)] = 0.0f;
    _S464[int(2)] = 0.0f;
    _S464[int(3)] = 0.0f;
    _S464[int(4)] = 0.0f;
    _S464[int(5)] = 0.0f;
    _S464[int(6)] = 0.0f;
    _S464[int(7)] = 0.0f;
    _S464[int(8)] = 0.0f;
    _S464[int(9)] = 0.0f;
    _S464[int(10)] = 0.0f;
    _S464[int(11)] = 0.0f;
    _S464[int(12)] = 0.0f;
    _S464[int(13)] = 0.0f;
    _S464[int(14)] = 0.0f;
    _S464[int(15)] = 0.0f;
    _S464[int(16)] = 0.0f;
    _S464[int(17)] = 0.0f;
    _S464[int(18)] = 0.0f;
    _S464[int(19)] = 0.0f;
    _S464[int(20)] = 0.0f;
    _S464[int(21)] = 0.0f;
    _S464[int(22)] = 0.0f;
    _S464[int(18)] = _S289;
    _S464[int(1)] = _S458;
    _S464[int(17)] = _S462.differential_0;
    _S464[int(0)] = _S460;
    FixedArray<float, 23>  _S465 = {
        _S387[int(0)] + _S464[int(0)], _S387[int(1)] + _S464[int(1)], _S387[int(2)] + _S464[int(2)], _S387[int(3)] + _S464[int(3)], _S387[int(4)] + _S464[int(4)], _S387[int(5)] + _S464[int(5)], _S387[int(6)] + _S464[int(6)], _S387[int(7)] + _S464[int(7)], _S387[int(8)] + _S464[int(8)], _S387[int(9)] + _S464[int(9)], _S387[int(10)] + _S464[int(10)], _S387[int(11)] + _S464[int(11)], _S387[int(12)] + _S464[int(12)], _S387[int(13)] + _S464[int(13)], _S387[int(14)] + _S464[int(14)], _S387[int(15)] + _S464[int(15)], _S387[int(16)] + _S464[int(16)], _S387[int(17)] + _S464[int(17)], _S387[int(18)] + _S464[int(18)], _S387[int(19)] + _S464[int(19)], _S387[int(20)] + _S464[int(20)], _S387[int(21)] + _S464[int(21)], _S387[int(22)] + _S464[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S465;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S466, FixedArray<float, 11>  * _S467, FixedArray<float, 10>  * _S468)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S466, _S467, _S468);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 11>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S469)
{
    FixedArray<float, 23>  _S470 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S470;
    FixedArray<float, 11>  _S471 = weights_5;
    FixedArray<float, 10>  _S472 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S471, &_S472);
    *_S469 = (&dp_raw_losses_0)->differential_0;
    return;
}

