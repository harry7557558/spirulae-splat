#pragma once

#include "slang.cuh"

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S1 = _slang_select(((*dpx_0).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_0).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S1;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_1, float3  dOut_1)
{
    float3  _S2 = _slang_select(((*dpx_1).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_1).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_1;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S2;
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

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, float dOut_2)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_0).primal_0.x * dOut_2;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_2).primal_0.x * dOut_2;
    *&((&x_d_result_0)->y) = (*dpy_0).primal_0.y * dOut_2;
    *&((&y_d_result_0)->y) = (*dpx_2).primal_0.y * dOut_2;
    *&((&x_d_result_0)->z) = (*dpy_0).primal_0.z * dOut_2;
    *&((&y_d_result_0)->z) = (*dpx_2).primal_0.z * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = x_d_result_0;
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

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_1, float dOut_3)
{
    DiffPair_float_0 _S3 = *dpx_3;
    float _S4;
    if(((*dpx_3).primal_0) > ((*dpy_1).primal_0))
    {
        _S4 = dOut_3;
    }
    else
    {
        if(((*dpx_3).primal_0) < ((*dpy_1).primal_0))
        {
            _S4 = 0.0f;
        }
        else
        {
            _S4 = 0.5f * dOut_3;
        }
    }
    dpx_3->primal_0 = _S3.primal_0;
    dpx_3->differential_0 = _S4;
    DiffPair_float_0 _S5 = *dpy_1;
    if(((*dpy_1).primal_0) > (_S3.primal_0))
    {
        _S4 = dOut_3;
    }
    else
    {
        if(((*dpy_1).primal_0) < ((*dpx_3).primal_0))
        {
            _S4 = 0.0f;
        }
        else
        {
            _S4 = 0.5f * dOut_3;
        }
    }
    dpy_1->primal_0 = _S5.primal_0;
    dpy_1->differential_0 = _S4;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_4, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_4)
{
    DiffPair_float_0 _S6 = *dpx_4;
    bool _S7;
    if(((*dpx_4).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S7 = ((*dpx_4).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S7 = false;
    }
    float _S8;
    if(_S7)
    {
        _S8 = dOut_4;
    }
    else
    {
        _S8 = 0.0f;
    }
    dpx_4->primal_0 = _S6.primal_0;
    dpx_4->differential_0 = _S8;
    DiffPair_float_0 _S9 = *dpMin_0;
    if((_S6.primal_0) < ((*dpMin_0).primal_0))
    {
        _S8 = dOut_4;
    }
    else
    {
        _S8 = 0.0f;
    }
    dpMin_0->primal_0 = _S9.primal_0;
    dpMin_0->differential_0 = _S8;
    DiffPair_float_0 _S10 = *dpMax_0;
    if(((*dpx_4).primal_0) > ((*dpMax_0).primal_0))
    {
        _S8 = dOut_4;
    }
    else
    {
        _S8 = 0.0f;
    }
    dpMax_0->primal_0 = _S10.primal_0;
    dpMax_0->differential_0 = _S8;
    return;
}

inline __device__ float clamp_0(float x_2, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_2), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_5, float dOut_5)
{
    float _S11 = 1.0f / (*dpx_5).primal_0 * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S11;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_6, float dOut_6)
{
    float _S12 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_6).primal_0)))))) * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S12;
    return;
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_7, float dOut_7)
{
    float _S13 = -0.5f / ((*dpx_7).primal_0 * (F32_sqrt(((*dpx_7).primal_0)))) * dOut_7;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S13;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpy_2, DiffPair_float_0 * dps_0, float dOut_8)
{
    float _S14 = (1.0f - (*dps_0).primal_0) * dOut_8;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S14;
    DiffPair_float_0 _S15 = *dpy_2;
    float _S16 = (*dps_0).primal_0 * dOut_8;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = _S16;
    float _S17 = (_S15.primal_0 - (*dpx_8).primal_0) * dOut_8;
    dps_0->primal_0 = _S15.primal_0;
    dps_0->differential_0 = _S17;
    return;
}

inline __device__ float lerp_0(float x_3, float y_1, float s_0)
{
    return x_3 + (y_1 - x_3) * s_0;
}

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_Ts_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, bool ref_alpha_0, bool has_mask_0, FixedArray<float, 15>  weights_0, FixedArray<float, 23>  * _S18)
{
    bool _S19;
    bool _S20;
    FixedArray<float, 23>  losses_0;
    bool mask_0;
    if(has_mask_0)
    {
        mask_0 = ref_alpha_0;
    }
    else
    {
        mask_0 = true;
    }
    float3  _S21;
    bool depth_mask_0 = ref_depth_0 != 0.0f;
    bool normal_mask_0 = (ref_normal_0.x + ref_normal_0.y + ref_normal_0.z) > -2.36599993705749512f;
    float _S22 = render_rgb_0.x;
    float _S23 = render_rgb_0.y;
    float _S24 = render_rgb_0.z;
    float _S25 = ref_rgb_0.x;
    float _S26 = ref_rgb_0.y;
    float _S27 = ref_rgb_0.z;
    float dY_0 = 0.29899999499320984f * _S22 + 0.58700001239776611f * _S23 + 0.11400000005960464f * _S24 - (0.29899999499320984f * _S25 + 0.58700001239776611f * _S26 + 0.11400000005960464f * _S27);
    float dU_0 = -0.14712999761104584f * _S22 - 0.28885999321937561f * _S23 + 0.43599998950958252f * _S24 - (-0.14712999761104584f * _S25 - 0.28885999321937561f * _S26 + 0.43599998950958252f * _S27);
    float dV_0 = 0.61500000953674316f * _S22 - 0.51498997211456299f * _S23 - 0.10001000016927719f * _S24 - (0.61500000953674316f * _S25 - 0.51498997211456299f * _S26 - 0.10001000016927719f * _S27);
    float _S28 = float(mask_0);
    float3  _S29 = ref_rgb_0 - render_rgb_0;
    float3  _S30 = abs_0(_S29);
    float _S31 = dot_0(_S29, _S29) * 0.3333333432674408f;
    losses_0[int(0)] = _S28 * (weights_0[int(0)] * ((_S30.x + _S30.y + _S30.z) * 0.3333333432674408f) + weights_0[int(1)] * _S31 + weights_0[int(2)] * (F32_abs((dY_0))) + weights_0[int(3)] * dY_0 * dY_0 + weights_0[int(4)] * dU_0 * dU_0 + weights_0[int(5)] * dV_0 * dV_0);
    losses_0[int(1)] = _S28 * clamp_0(_S31, 0.0f, 1.0f);
    float _S32 = float(depth_mask_0 & mask_0);
    float _S33 = _S32 * (F32_log(((F32_max((render_depth_0), (0.00009999999747379f))))));
    float _S34 = _S32 * (F32_log(((F32_max((ref_depth_0), (0.00009999999747379f))))));
    losses_0[int(2)] = _S33;
    losses_0[int(3)] = _S34;
    losses_0[int(4)] = _S33 * _S33;
    losses_0[int(5)] = _S34 * _S34;
    losses_0[int(6)] = _S33 * _S34;
    bool _S35 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S36 = norm2_0 == 0.0f;
        _S19 = _S36;
        if(_S36)
        {
            _S21 = make_float3 (0.0f);
            break;
        }
        _S21 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S37;
    bool _S38 = !_S19;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S39 = norm2_1 == 0.0f;
        _S20 = _S39;
        if(_S39)
        {
            _S37 = make_float3 (0.0f);
            break;
        }
        _S37 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    float3  _S40;
    bool normal_mask_1;
    bool _S41 = !_S20;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S40 = make_float3 (0.0f);
            normal_mask_1 = false;
            break;
        }
        _S40 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        normal_mask_1 = _S35;
        break;
    }
    float _S42 = float(_S38 & normal_mask_1);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S21, _S40);
    losses_0[int(7)] = weights_0[int(7)] * _S42 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S43 = float(_S41 & normal_mask_1);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S37, _S40);
    losses_0[int(8)] = weights_0[int(7)] * _S43 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S44 = float(_S38 & _S41);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S21, _S37);
    losses_0[int(11)] = weights_0[int(10)] * _S44 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    float render_alpha_0 = clamp_0(1.0f - render_Ts_0, 0.0f, 1.0f);
    float _S45 = float(depth_mask_0);
    float _S46 = float(ref_alpha_0);
    float _S47 = (F32_max((render_alpha_0), (_S46)));
    losses_0[int(9)] = weights_0[int(8)] * _S45 * - lerp_0((F32_log(((F32_max((1.0f - _S47), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S47), (9.99999997475242708e-07f)))))), _S46);
    float _S48 = 1.0f - render_alpha_0;
    float _S49 = 1.0f - _S46;
    float _S50 = (F32_max((_S48), (_S49)));
    losses_0[int(10)] = weights_0[int(9)] * _S45 * - lerp_0((F32_log(((F32_max((1.0f - _S50), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S50), (9.99999997475242708e-07f)))))), _S49);
    losses_0[int(12)] = weights_0[int(11)] * 4.0f * render_alpha_0 * _S48;
    losses_0[int(13)] = weights_0[int(12)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f);
    losses_0[int(14)] = weights_0[int(13)] * depth_dist_0;
    losses_0[int(15)] = weights_0[int(14)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f);
    losses_0[int(16)] = 1.0f;
    losses_0[int(17)] = _S28;
    losses_0[int(18)] = _S32;
    losses_0[int(19)] = _S42;
    losses_0[int(20)] = _S43;
    losses_0[int(21)] = _S44;
    losses_0[int(22)] = _S45;
    *_S18 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S51, float3  _S52)
{
    return dot_0(_S51, _S52);
}

inline __device__ float s_primal_ctx_log_0(float _S53)
{
    return (F32_log((_S53)));
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S54)
{
    return (F32_rsqrt((_S54)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S55, float _S56, float _S57)
{
    return clamp_0(_S55, _S56, _S57);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S58, DiffPair_float_0 * _S59, DiffPair_float_0 * _S60, float _S61)
{
    _d_lerp_0(_S58, _S59, _S60, _S61);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S62, float _S63)
{
    _d_log_0(_S62, _S63);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S64, DiffPair_float_0 * _S65, DiffPair_float_0 * _S66, float _S67)
{
    _d_clamp_0(_S64, _S65, _S66, _S67);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S68, float _S69)
{
    _d_sqrt_0(_S68, _S69);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S70, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S71, float _S72)
{
    _d_dot_0(_S70, _S71, _S72);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S73, float _S74)
{
    _d_rsqrt_0(_S73, _S74);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S75, float _S76)
{
    _d_abs_0(_S75, _S76);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S77, float3  _S78)
{
    _d_abs_vector_0(_S77, _S78);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_Ts_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, bool ref_alpha_1, bool has_mask_1, FixedArray<float, 15>  * weights_1, FixedArray<float, 23>  * _s_dOut_0)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S79 = *dprender_rgb_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S80 = *dpref_rgb_0;
    DiffPair_float_0 _S81 = *dprender_depth_0;
    DiffPair_float_0 _S82 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S83 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S84 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S85 = *dpref_normal_0;
    DiffPair_float_0 _S86 = *dprender_Ts_0;
    float3  _S87 = make_float3 (0.0f);
    bool mask_1;
    if(has_mask_1)
    {
        mask_1 = ref_alpha_1;
    }
    else
    {
        mask_1 = true;
    }
    bool depth_mask_1 = (_S82.primal_0) != 0.0f;
    float _S88 = _S79.primal_0.x;
    float _S89 = _S79.primal_0.y;
    float _S90 = _S79.primal_0.z;
    float _S91 = _S80.primal_0.x;
    float _S92 = _S80.primal_0.y;
    float _S93 = _S80.primal_0.z;
    float dY_1 = 0.29899999499320984f * _S88 + 0.58700001239776611f * _S89 + 0.11400000005960464f * _S90 - (0.29899999499320984f * _S91 + 0.58700001239776611f * _S92 + 0.11400000005960464f * _S93);
    float dU_1 = -0.14712999761104584f * _S88 - 0.28885999321937561f * _S89 + 0.43599998950958252f * _S90 - (-0.14712999761104584f * _S91 - 0.28885999321937561f * _S92 + 0.43599998950958252f * _S93);
    float dV_1 = 0.61500000953674316f * _S88 - 0.51498997211456299f * _S89 - 0.10001000016927719f * _S90 - (0.61500000953674316f * _S91 - 0.51498997211456299f * _S92 - 0.10001000016927719f * _S93);
    float _S94 = float(mask_1);
    float _S95 = (*weights_1)[int(0)];
    float3  _S96 = _S80.primal_0 - _S79.primal_0;
    float _S97 = (*weights_1)[int(1)];
    float _S98 = s_primal_ctx_dot_0(_S96, _S96) * 0.3333333432674408f;
    float _S99 = (*weights_1)[int(2)];
    float _S100 = (*weights_1)[int(3)];
    float _S101 = (*weights_1)[int(3)] * dY_1;
    float _S102 = (*weights_1)[int(4)];
    float _S103 = (*weights_1)[int(4)] * dU_1;
    float _S104 = (*weights_1)[int(5)];
    float _S105 = (*weights_1)[int(5)] * dV_1;
    float _S106 = float(depth_mask_1 & mask_1);
    float _S107 = (F32_max((_S81.primal_0), (0.00009999999747379f)));
    float _S108 = _S106 * s_primal_ctx_log_0(_S107);
    float _S109 = (F32_max((_S82.primal_0), (0.00009999999747379f)));
    float _S110 = _S106 * s_primal_ctx_log_0(_S109);
    bool _S111 = ((_S85.primal_0.x + _S85.primal_0.y + _S85.primal_0.z) > -2.36599993705749512f) & mask_1;
    float _S112 = s_primal_ctx_dot_0(_S83.primal_0, _S83.primal_0);
    bool _S113 = _S112 == 0.0f;
    float3  _S114;
    if(_S113)
    {
        _S114 = make_float3 (0.0f);
    }
    bool _S115 = !_S113;
    float3  _S116;
    if(_S115)
    {
        float _S117 = s_primal_ctx_rsqrt_0(_S112);
        float3  _S118 = make_float3 (_S117);
        _S114 = _S83.primal_0 * make_float3 (_S117);
        _S116 = _S118;
    }
    else
    {
        _S116 = _S87;
    }
    float _S119 = s_primal_ctx_dot_0(_S84.primal_0, _S84.primal_0);
    bool _S120 = _S119 == 0.0f;
    float3  _S121;
    if(_S120)
    {
        _S121 = make_float3 (0.0f);
    }
    bool _S122 = !_S120;
    float3  _S123;
    if(_S122)
    {
        float _S124 = s_primal_ctx_rsqrt_0(_S119);
        float3  _S125 = make_float3 (_S124);
        _S121 = _S84.primal_0 * make_float3 (_S124);
        _S123 = _S125;
    }
    else
    {
        _S123 = _S87;
    }
    float _S126 = s_primal_ctx_dot_0(_S85.primal_0, _S85.primal_0);
    bool _S127 = _S126 == 0.0f;
    bool normal_mask_2;
    float3  _S128;
    if(_S127)
    {
        float3  _S129 = make_float3 (0.0f);
        normal_mask_2 = false;
        _S128 = _S129;
    }
    else
    {
        normal_mask_2 = _S111;
    }
    bool _S130 = !_S127;
    float3  _S131;
    if(_S130)
    {
        float _S132 = s_primal_ctx_rsqrt_0(_S126);
        float3  _S133 = make_float3 (_S132);
        _S128 = _S85.primal_0 * make_float3 (_S132);
        _S131 = _S133;
    }
    else
    {
        _S131 = _S87;
    }
    float _S134 = (*weights_1)[int(7)] * float(_S115 & normal_mask_2);
    float cos_sim_loss_3 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S114, _S128);
    float _S135 = (F32_max((cos_sim_loss_3), (9.999999960041972e-13f)));
    float _S136 = (*weights_1)[int(7)] * float(_S122 & normal_mask_2);
    float cos_sim_loss_4 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S121, _S128);
    float _S137 = (F32_max((cos_sim_loss_4), (9.999999960041972e-13f)));
    float _S138 = (*weights_1)[int(10)] * float(_S115 & _S122);
    float cos_sim_loss_5 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S114, _S121);
    float _S139 = (F32_max((cos_sim_loss_5), (9.999999960041972e-13f)));
    float _S140 = 1.0f - _S86.primal_0;
    float _S141 = s_primal_ctx_clamp_0(_S140, 0.0f, 1.0f);
    float _S142 = float(depth_mask_1);
    float _S143 = (*weights_1)[int(8)] * _S142;
    float _S144 = float(ref_alpha_1);
    float _S145 = (F32_max((_S141), (_S144)));
    float _S146 = 1.0f - _S145;
    float _S147 = (F32_max((_S146), (9.99999997475242708e-07f)));
    float _S148 = s_primal_ctx_log_0(_S147);
    float _S149 = (F32_max((_S145), (9.99999997475242708e-07f)));
    float _S150 = s_primal_ctx_log_0(_S149);
    float _S151 = 1.0f - _S141;
    float _S152 = 1.0f - _S144;
    float _S153 = (F32_max((_S151), (_S152)));
    float _S154 = 1.0f - _S153;
    float _S155 = (F32_max((_S154), (9.99999997475242708e-07f)));
    float _S156 = (F32_max((_S153), (9.99999997475242708e-07f)));
    float _S157 = s_primal_ctx_log_0(_S156);
    float _S158 = (*weights_1)[int(11)] * 4.0f;
    float _S159 = (*_s_dOut_0)[int(0)];
    float _S160 = (*_s_dOut_0)[int(1)];
    float _S161 = (*_s_dOut_0)[int(2)];
    float _S162 = (*_s_dOut_0)[int(3)];
    float _S163 = (*_s_dOut_0)[int(4)];
    float _S164 = (*_s_dOut_0)[int(5)];
    float _S165 = (*_s_dOut_0)[int(6)];
    float _S166 = 0.3333333432674408f * ((*weights_1)[int(14)] * (*_s_dOut_0)[int(15)]);
    float _S167 = (*weights_1)[int(13)] * (*_s_dOut_0)[int(14)];
    float _S168 = 0.3333333432674408f * ((*weights_1)[int(12)] * (*_s_dOut_0)[int(13)]);
    float _S169 = _S158 * _S141 * (*_s_dOut_0)[int(12)];
    float _S170 = _S158 * (_S151 * (*_s_dOut_0)[int(12)]);
    float _S171 = - ((*weights_1)[int(9)] * _S142 * (*_s_dOut_0)[int(10)]);
    DiffPair_float_0 _S172;
    (&_S172)->primal_0 = s_primal_ctx_log_0(_S155);
    (&_S172)->differential_0 = 0.0f;
    DiffPair_float_0 _S173;
    (&_S173)->primal_0 = _S157;
    (&_S173)->differential_0 = 0.0f;
    DiffPair_float_0 _S174;
    (&_S174)->primal_0 = _S152;
    (&_S174)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S172, &_S173, &_S174, _S171);
    DiffPair_float_0 _S175;
    (&_S175)->primal_0 = _S156;
    (&_S175)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S175, _S173.differential_0);
    DiffPair_float_0 _S176;
    (&_S176)->primal_0 = _S153;
    (&_S176)->differential_0 = 0.0f;
    DiffPair_float_0 _S177;
    (&_S177)->primal_0 = 9.99999997475242708e-07f;
    (&_S177)->differential_0 = 0.0f;
    _d_max_0(&_S176, &_S177, _S175.differential_0);
    DiffPair_float_0 _S178;
    (&_S178)->primal_0 = _S155;
    (&_S178)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S178, _S172.differential_0);
    DiffPair_float_0 _S179;
    (&_S179)->primal_0 = _S154;
    (&_S179)->differential_0 = 0.0f;
    DiffPair_float_0 _S180;
    (&_S180)->primal_0 = 9.99999997475242708e-07f;
    (&_S180)->differential_0 = 0.0f;
    _d_max_0(&_S179, &_S180, _S178.differential_0);
    float _S181 = _S176.differential_0 + - _S179.differential_0;
    DiffPair_float_0 _S182;
    (&_S182)->primal_0 = _S151;
    (&_S182)->differential_0 = 0.0f;
    DiffPair_float_0 _S183;
    (&_S183)->primal_0 = _S152;
    (&_S183)->differential_0 = 0.0f;
    _d_max_0(&_S182, &_S183, _S181);
    float _S184 = - (_S169 + _S182.differential_0);
    float _S185 = - (_S143 * (*_s_dOut_0)[int(9)]);
    DiffPair_float_0 _S186;
    (&_S186)->primal_0 = _S148;
    (&_S186)->differential_0 = 0.0f;
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = _S150;
    (&_S187)->differential_0 = 0.0f;
    DiffPair_float_0 _S188;
    (&_S188)->primal_0 = _S144;
    (&_S188)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S186, &_S187, &_S188, _S185);
    DiffPair_float_0 _S189;
    (&_S189)->primal_0 = _S149;
    (&_S189)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S189, _S187.differential_0);
    DiffPair_float_0 _S190;
    (&_S190)->primal_0 = _S145;
    (&_S190)->differential_0 = 0.0f;
    DiffPair_float_0 _S191;
    (&_S191)->primal_0 = 9.99999997475242708e-07f;
    (&_S191)->differential_0 = 0.0f;
    _d_max_0(&_S190, &_S191, _S189.differential_0);
    DiffPair_float_0 _S192;
    (&_S192)->primal_0 = _S147;
    (&_S192)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S192, _S186.differential_0);
    DiffPair_float_0 _S193;
    (&_S193)->primal_0 = _S146;
    (&_S193)->differential_0 = 0.0f;
    DiffPair_float_0 _S194;
    (&_S194)->primal_0 = 9.99999997475242708e-07f;
    (&_S194)->differential_0 = 0.0f;
    _d_max_0(&_S193, &_S194, _S192.differential_0);
    float _S195 = _S190.differential_0 + - _S193.differential_0;
    DiffPair_float_0 _S196;
    (&_S196)->primal_0 = _S141;
    (&_S196)->differential_0 = 0.0f;
    DiffPair_float_0 _S197;
    (&_S197)->primal_0 = _S144;
    (&_S197)->differential_0 = 0.0f;
    _d_max_0(&_S196, &_S197, _S195);
    float _S198 = _S170 + _S184 + _S196.differential_0;
    DiffPair_float_0 _S199;
    (&_S199)->primal_0 = _S140;
    (&_S199)->differential_0 = 0.0f;
    DiffPair_float_0 _S200;
    (&_S200)->primal_0 = 0.0f;
    (&_S200)->differential_0 = 0.0f;
    DiffPair_float_0 _S201;
    (&_S201)->primal_0 = 1.0f;
    (&_S201)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S199, &_S200, &_S201, _S198);
    float _S202 = - _S199.differential_0;
    float _S203 = _S138 * (*_s_dOut_0)[int(11)];
    DiffPair_float_0 _S204;
    (&_S204)->primal_0 = _S139;
    (&_S204)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S204, _S203);
    DiffPair_float_0 _S205;
    (&_S205)->primal_0 = cos_sim_loss_5;
    (&_S205)->differential_0 = 0.0f;
    DiffPair_float_0 _S206;
    (&_S206)->primal_0 = 9.999999960041972e-13f;
    (&_S206)->differential_0 = 0.0f;
    _d_max_0(&_S205, &_S206, _S204.differential_0);
    float _S207 = 0.5f * - (_S203 + _S205.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S208;
    (&_S208)->primal_0 = _S114;
    (&_S208)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S209;
    (&_S209)->primal_0 = _S121;
    (&_S209)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S208, &_S209, _S207);
    float _S210 = _S136 * (*_s_dOut_0)[int(8)];
    DiffPair_float_0 _S211;
    (&_S211)->primal_0 = _S137;
    (&_S211)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S211, _S210);
    DiffPair_float_0 _S212;
    (&_S212)->primal_0 = cos_sim_loss_4;
    (&_S212)->differential_0 = 0.0f;
    DiffPair_float_0 _S213;
    (&_S213)->primal_0 = 9.999999960041972e-13f;
    (&_S213)->differential_0 = 0.0f;
    _d_max_0(&_S212, &_S213, _S211.differential_0);
    float _S214 = 0.5f * - (_S210 + _S212.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S215;
    (&_S215)->primal_0 = _S121;
    (&_S215)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S216;
    (&_S216)->primal_0 = _S128;
    (&_S216)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S215, &_S216, _S214);
    float _S217 = _S134 * (*_s_dOut_0)[int(7)];
    DiffPair_float_0 _S218;
    (&_S218)->primal_0 = _S135;
    (&_S218)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S218, _S217);
    DiffPair_float_0 _S219;
    (&_S219)->primal_0 = cos_sim_loss_3;
    (&_S219)->differential_0 = 0.0f;
    DiffPair_float_0 _S220;
    (&_S220)->primal_0 = 9.999999960041972e-13f;
    (&_S220)->differential_0 = 0.0f;
    _d_max_0(&_S219, &_S220, _S218.differential_0);
    float _S221 = 0.5f * - (_S217 + _S219.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S222;
    (&_S222)->primal_0 = _S114;
    (&_S222)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S223;
    (&_S223)->primal_0 = _S128;
    (&_S223)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S222, &_S223, _S221);
    float3  _S224 = _S216.differential_0 + _S223.differential_0;
    float3  _S225 = _S208.differential_0 + _S222.differential_0;
    float3  _S226 = make_float3 (_S166, _S166, _S166);
    float3  _S227 = make_float3 (_S168, _S168, _S168);
    float3  _S228 = _S209.differential_0 + _S215.differential_0;
    float _S229;
    if(_S130)
    {
        float3  _S230 = _S85.primal_0 * _S224;
        float3  _S231 = _S131 * _S224;
        float _S232 = _S230.x + _S230.y + _S230.z;
        DiffPair_float_0 _S233;
        (&_S233)->primal_0 = _S126;
        (&_S233)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S233, _S232);
        _S229 = _S233.differential_0;
        _S114 = _S231;
    }
    else
    {
        _S229 = 0.0f;
        _S114 = _S87;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S234;
    (&_S234)->primal_0 = _S85.primal_0;
    (&_S234)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S235;
    (&_S235)->primal_0 = _S85.primal_0;
    (&_S235)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S234, &_S235, _S229);
    float3  _S236 = _S235.differential_0 + _S234.differential_0 + _S114;
    if(_S122)
    {
        float3  _S237 = _S84.primal_0 * _S228;
        float3  _S238 = _S123 * _S228;
        float _S239 = _S237.x + _S237.y + _S237.z;
        DiffPair_float_0 _S240;
        (&_S240)->primal_0 = _S119;
        (&_S240)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S240, _S239);
        _S229 = _S240.differential_0;
        _S114 = _S238;
    }
    else
    {
        _S229 = 0.0f;
        _S114 = _S87;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S241;
    (&_S241)->primal_0 = _S84.primal_0;
    (&_S241)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S84.primal_0;
    (&_S242)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S241, &_S242, _S229);
    float3  _S243 = _S242.differential_0 + _S241.differential_0 + _S114;
    if(_S115)
    {
        float3  _S244 = _S83.primal_0 * _S225;
        float3  _S245 = _S116 * _S225;
        float _S246 = _S244.x + _S244.y + _S244.z;
        DiffPair_float_0 _S247;
        (&_S247)->primal_0 = _S112;
        (&_S247)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S247, _S246);
        _S229 = _S247.differential_0;
        _S114 = _S245;
    }
    else
    {
        _S229 = 0.0f;
        _S114 = _S87;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S248;
    (&_S248)->primal_0 = _S83.primal_0;
    (&_S248)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S249;
    (&_S249)->primal_0 = _S83.primal_0;
    (&_S249)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S248, &_S249, _S229);
    float3  _S250 = _S249.differential_0 + _S248.differential_0 + _S114;
    float _S251 = _S110 * _S165;
    float _S252 = _S110 * _S164;
    float _S253 = _S108 * _S163;
    float _S254 = _S106 * (_S108 * _S165 + _S252 + _S252 + _S162);
    DiffPair_float_0 _S255;
    (&_S255)->primal_0 = _S109;
    (&_S255)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S255, _S254);
    DiffPair_float_0 _S256;
    (&_S256)->primal_0 = _S82.primal_0;
    (&_S256)->differential_0 = 0.0f;
    DiffPair_float_0 _S257;
    (&_S257)->primal_0 = 0.00009999999747379f;
    (&_S257)->differential_0 = 0.0f;
    _d_max_0(&_S256, &_S257, _S255.differential_0);
    float _S258 = _S106 * (_S251 + _S253 + _S253 + _S161);
    DiffPair_float_0 _S259;
    (&_S259)->primal_0 = _S107;
    (&_S259)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S259, _S258);
    DiffPair_float_0 _S260;
    (&_S260)->primal_0 = _S81.primal_0;
    (&_S260)->differential_0 = 0.0f;
    DiffPair_float_0 _S261;
    (&_S261)->primal_0 = 0.00009999999747379f;
    (&_S261)->differential_0 = 0.0f;
    _d_max_0(&_S260, &_S261, _S259.differential_0);
    float _S262 = _S94 * _S160;
    DiffPair_float_0 _S263;
    (&_S263)->primal_0 = _S98;
    (&_S263)->differential_0 = 0.0f;
    DiffPair_float_0 _S264;
    (&_S264)->primal_0 = 0.0f;
    (&_S264)->differential_0 = 0.0f;
    DiffPair_float_0 _S265;
    (&_S265)->primal_0 = 1.0f;
    (&_S265)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S263, &_S264, &_S265, _S262);
    float _S266 = _S94 * _S159;
    float _S267 = _S105 * _S266;
    float _S268 = _S104 * (dV_1 * _S266);
    float _S269 = _S103 * _S266;
    float _S270 = _S102 * (dU_1 * _S266);
    float _S271 = _S101 * _S266;
    float _S272 = _S100 * (dY_1 * _S266);
    float _S273 = _S99 * _S266;
    DiffPair_float_0 _S274;
    (&_S274)->primal_0 = dY_1;
    (&_S274)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S274, _S273);
    float _S275 = 0.3333333432674408f * (_S263.differential_0 + _S97 * _S266);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S276;
    (&_S276)->primal_0 = _S96;
    (&_S276)->differential_0 = _S87;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S277;
    (&_S277)->primal_0 = _S96;
    (&_S277)->differential_0 = _S87;
    s_bwd_prop_dot_0(&_S276, &_S277, _S275);
    float _S278 = 0.3333333432674408f * (_S95 * _S266);
    float3  _S279 = make_float3 (_S278, _S278, _S278);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S280;
    (&_S280)->primal_0 = _S96;
    (&_S280)->differential_0 = _S87;
    s_bwd_prop_abs_1(&_S280, _S279);
    float3  _S281 = _S277.differential_0 + _S276.differential_0 + _S280.differential_0;
    float _S282 = _S267 + _S268;
    float s_diff_V_T_0 = - _S282;
    float _S283 = _S269 + _S270;
    float s_diff_U_T_0 = - _S283;
    float _S284 = _S271 + _S272 + _S274.differential_0;
    float s_diff_Y_T_0 = - _S284;
    float _S285 = - s_diff_V_T_0;
    float3  _S286 = _S281 + make_float3 (0.61500000953674316f * s_diff_V_T_0 + -0.14712999761104584f * s_diff_U_T_0 + 0.29899999499320984f * s_diff_Y_T_0, 0.51498997211456299f * _S285 + 0.28885999321937561f * - s_diff_U_T_0 + 0.58700001239776611f * s_diff_Y_T_0, 0.10001000016927719f * _S285 + 0.43599998950958252f * s_diff_U_T_0 + 0.11400000005960464f * s_diff_Y_T_0);
    float3  _S287 = - _S281 + make_float3 (0.61500000953674316f * _S282 + -0.14712999761104584f * _S283 + 0.29899999499320984f * _S284, 0.51498997211456299f * s_diff_V_T_0 + 0.28885999321937561f * s_diff_U_T_0 + 0.58700001239776611f * _S284, 0.10001000016927719f * s_diff_V_T_0 + 0.43599998950958252f * _S283 + 0.11400000005960464f * _S284);
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S226;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S167;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S227;
    dprender_Ts_0->primal_0 = (*dprender_Ts_0).primal_0;
    dprender_Ts_0->differential_0 = _S202;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S236;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S243;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S250;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S256.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S260.differential_0;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S286;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S287;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S288, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S289, DiffPair_float_0 * _S290, DiffPair_float_0 * _S291, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S292, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S293, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S294, DiffPair_float_0 * _S295, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S296, DiffPair_float_0 * _S297, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S298, bool _S299, bool _S300, FixedArray<float, 15>  * _S301, FixedArray<float, 23>  * _S302)
{
    s_bwd_prop_per_pixel_losses_0(_S288, _S289, _S290, _S291, _S292, _S293, _S294, _S295, _S296, _S297, _S298, _S299, _S300, _S301, _S302);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_Ts_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, bool ref_alpha_2, bool has_mask_2, FixedArray<float, 15>  weights_2, FixedArray<float, 23>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_Ts_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0)
{
    float3  _S303 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S303;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S303;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S303;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S303;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S303;
    DiffPair_float_0 dp_render_Ts_0;
    (&dp_render_Ts_0)->primal_0 = render_Ts_1;
    (&dp_render_Ts_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S303;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S303;
    FixedArray<float, 15>  _S304 = weights_2;
    FixedArray<float, 23>  _S305 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_Ts_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, ref_alpha_2, has_mask_2, &_S304, &_S305);
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

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_9, float dOut_9)
{
    float _S306 = 1.0f / ((*dpx_9).primal_0 * 2.30258512496948242f) * dOut_9;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S306;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 23>  raw_losses_0, FixedArray<float, 15>  weights_3, FixedArray<float, 10>  * _S307)
{
    FixedArray<float, 10>  losses_1;
    float _S308 = (F32_max((raw_losses_0[int(17)]), (1.0f)));
    losses_1[int(0)] = raw_losses_0[int(0)] / _S308;
    losses_1[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S308)));
    bool _S309;
    if((raw_losses_0[int(18)]) > 0.0f)
    {
        _S309 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S309 = false;
    }
    float _S310;
    if(_S309)
    {
        _S310 = weights_3[int(6)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(18)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(18)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S310 = 0.0f;
    }
    losses_1[int(2)] = _S310;
    losses_1[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(19)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(20)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(19)]) > 0.5f) + int((raw_losses_0[int(20)]) > 0.5f)), (int(1)))));
    losses_1[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(22)]), (1.0f)));
    losses_1[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(21)]), (1.0f)));
    float _S311 = (F32_max((raw_losses_0[int(16)]), (1.0f)));
    losses_1[int(6)] = raw_losses_0[int(12)] / _S311;
    losses_1[int(7)] = raw_losses_0[int(13)] / _S311;
    losses_1[int(8)] = raw_losses_0[int(14)] / _S311;
    losses_1[int(9)] = raw_losses_0[int(15)] / _S311;
    *_S307 = losses_1;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S312)
{
    return (F32_sqrt((_S312)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S313, float _S314)
{
    _d_log10_0(_S313, _S314);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_0, FixedArray<float, 15>  * weights_4, FixedArray<float, 10>  * _s_dOut_1)
{
    FixedArray<float, 23>  _S315 = dpraw_losses_0->primal_0;
    float _S316 = (F32_max((dpraw_losses_0->primal_0[int(17)]), (1.0f)));
    float _S317 = _S316 * _S316;
    float _S318 = dpraw_losses_0->primal_0[int(1)] / _S316;
    bool _S319 = (dpraw_losses_0->primal_0[int(18)]) > 0.0f;
    bool _S320;
    if(_S319)
    {
        _S320 = (_S315[int(3)]) != 0.0f;
    }
    else
    {
        _S320 = false;
    }
    float _S321;
    float _S322;
    float _S323;
    float _S324;
    float _S325;
    float _S326;
    float _S327;
    float _S328;
    float _S329;
    float _S330;
    float _S331;
    float _S332;
    float _S333;
    float _S334;
    float _S335;
    if(_S320)
    {
        float _S336 = _S315[int(2)] * _S315[int(3)];
        float _S337 = _S315[int(18)] * _S315[int(18)];
        float _S338 = _S315[int(6)] - _S336 / _S315[int(18)];
        float _S339 = _S315[int(2)] * _S315[int(2)];
        float _S340 = _S315[int(4)] - _S339 / _S315[int(18)];
        float _S341 = _S315[int(3)] * _S315[int(3)];
        float _S342 = _S315[int(5)] - _S341 / _S315[int(18)];
        float _S343 = _S340 * _S342 + 1.0f;
        float _S344 = (F32_max((9.999999960041972e-13f), (_S343)));
        float _S345 = s_primal_ctx_sqrt_0(_S344);
        float _S346 = _S345 * _S345;
        float _S347 = 1.0f - _S338 / _S345;
        _S321 = (*weights_4)[int(6)];
        _S322 = _S347;
        _S323 = _S346;
        _S324 = _S338;
        _S325 = _S345;
        _S326 = _S344;
        _S327 = _S343;
        _S328 = _S340;
        _S329 = _S342;
        _S330 = _S337;
        _S331 = _S341;
        _S332 = _S315[int(3)];
        _S333 = _S339;
        _S334 = _S315[int(2)];
        _S335 = _S336;
    }
    else
    {
        _S321 = 0.0f;
        _S322 = 0.0f;
        _S323 = 0.0f;
        _S324 = 0.0f;
        _S325 = 0.0f;
        _S326 = 0.0f;
        _S327 = 0.0f;
        _S328 = 0.0f;
        _S329 = 0.0f;
        _S330 = 0.0f;
        _S331 = 0.0f;
        _S332 = 0.0f;
        _S333 = 0.0f;
        _S334 = 0.0f;
        _S335 = 0.0f;
    }
    float _S348 = (F32_max((_S315[int(19)]), (1.0f)));
    float _S349 = _S348 * _S348;
    float _S350 = (F32_max((_S315[int(20)]), (1.0f)));
    float _S351 = _S350 * _S350;
    float _S352 = float((I32_max((int((_S315[int(19)]) > 0.5f) + int((_S315[int(20)]) > 0.5f)), (int(1)))));
    float _S353 = _S315[int(9)] + _S315[int(10)];
    float _S354 = (F32_max((_S315[int(22)]), (1.0f)));
    float _S355 = _S354 * _S354;
    float _S356 = (F32_max((_S315[int(21)]), (1.0f)));
    float _S357 = _S356 * _S356;
    float _S358 = (F32_max((_S315[int(16)]), (1.0f)));
    float _S359 = _S358 * _S358;
    float _S360 = (*_s_dOut_1)[int(0)];
    float _S361 = (*_s_dOut_1)[int(1)];
    float _S362 = (*_s_dOut_1)[int(2)];
    float _S363 = (*_s_dOut_1)[int(9)] / _S359;
    float _S364 = _S358 * _S363;
    float _S365 = (*_s_dOut_1)[int(8)] / _S359;
    float _S366 = _S358 * _S365;
    float _S367 = (*_s_dOut_1)[int(7)] / _S359;
    float _S368 = _S358 * _S367;
    float _S369 = (*_s_dOut_1)[int(6)] / _S359;
    float _S370 = _S358 * _S369;
    float _S371 = _S315[int(15)] * - _S363 + _S315[int(14)] * - _S365 + _S315[int(13)] * - _S367 + _S315[int(12)] * - _S369;
    DiffPair_float_0 _S372;
    (&_S372)->primal_0 = _S315[int(16)];
    (&_S372)->differential_0 = 0.0f;
    DiffPair_float_0 _S373;
    (&_S373)->primal_0 = 1.0f;
    (&_S373)->differential_0 = 0.0f;
    _d_max_0(&_S372, &_S373, _S371);
    float _S374 = (*_s_dOut_1)[int(5)] / _S357;
    float _S375 = _S315[int(11)] * - _S374;
    float _S376 = _S356 * _S374;
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = _S315[int(21)];
    (&_S377)->differential_0 = 0.0f;
    DiffPair_float_0 _S378;
    (&_S378)->primal_0 = 1.0f;
    (&_S378)->differential_0 = 0.0f;
    _d_max_0(&_S377, &_S378, _S375);
    float _S379 = (*_s_dOut_1)[int(4)] / _S355;
    float _S380 = _S353 * - _S379;
    float _S381 = _S354 * _S379;
    DiffPair_float_0 _S382;
    (&_S382)->primal_0 = _S315[int(22)];
    (&_S382)->differential_0 = 0.0f;
    DiffPair_float_0 _S383;
    (&_S383)->primal_0 = 1.0f;
    (&_S383)->differential_0 = 0.0f;
    _d_max_0(&_S382, &_S383, _S380);
    float _S384 = (*_s_dOut_1)[int(3)] / _S352;
    float _S385 = _S384 / _S351;
    float _S386 = _S315[int(8)] * - _S385;
    float _S387 = _S350 * _S385;
    DiffPair_float_0 _S388;
    (&_S388)->primal_0 = _S315[int(20)];
    (&_S388)->differential_0 = 0.0f;
    DiffPair_float_0 _S389;
    (&_S389)->primal_0 = 1.0f;
    (&_S389)->differential_0 = 0.0f;
    _d_max_0(&_S388, &_S389, _S386);
    float _S390 = _S384 / _S349;
    float _S391 = _S315[int(7)] * - _S390;
    float _S392 = _S348 * _S390;
    DiffPair_float_0 _S393;
    (&_S393)->primal_0 = _S315[int(19)];
    (&_S393)->differential_0 = 0.0f;
    DiffPair_float_0 _S394;
    (&_S394)->primal_0 = 1.0f;
    (&_S394)->differential_0 = 0.0f;
    _d_max_0(&_S393, &_S394, _S391);
    FixedArray<float, 23>  _S395;
    _S395[int(0)] = 0.0f;
    _S395[int(1)] = 0.0f;
    _S395[int(2)] = 0.0f;
    _S395[int(3)] = 0.0f;
    _S395[int(4)] = 0.0f;
    _S395[int(5)] = 0.0f;
    _S395[int(6)] = 0.0f;
    _S395[int(7)] = 0.0f;
    _S395[int(8)] = 0.0f;
    _S395[int(9)] = 0.0f;
    _S395[int(10)] = 0.0f;
    _S395[int(11)] = 0.0f;
    _S395[int(12)] = 0.0f;
    _S395[int(13)] = 0.0f;
    _S395[int(14)] = 0.0f;
    _S395[int(15)] = 0.0f;
    _S395[int(16)] = 0.0f;
    _S395[int(17)] = 0.0f;
    _S395[int(18)] = 0.0f;
    _S395[int(19)] = 0.0f;
    _S395[int(20)] = 0.0f;
    _S395[int(21)] = 0.0f;
    _S395[int(22)] = 0.0f;
    _S395[int(15)] = _S364;
    _S395[int(14)] = _S366;
    _S395[int(13)] = _S368;
    _S395[int(16)] = _S372.differential_0;
    _S395[int(12)] = _S370;
    _S395[int(21)] = _S377.differential_0;
    _S395[int(11)] = _S376;
    _S395[int(22)] = _S382.differential_0;
    _S395[int(10)] = _S381;
    _S395[int(9)] = _S381;
    _S395[int(20)] = _S388.differential_0;
    _S395[int(8)] = _S387;
    _S395[int(19)] = _S393.differential_0;
    _S395[int(7)] = _S392;
    float _S396 = _S395[int(0)];
    float _S397 = _S395[int(1)];
    float _S398 = _S395[int(2)];
    float _S399 = _S395[int(3)];
    float _S400 = _S395[int(4)];
    float _S401 = _S395[int(5)];
    float _S402 = _S395[int(6)];
    float _S403 = _S395[int(7)];
    float _S404 = _S395[int(8)];
    float _S405 = _S395[int(9)];
    float _S406 = _S395[int(10)];
    float _S407 = _S395[int(11)];
    float _S408 = _S395[int(12)];
    float _S409 = _S395[int(13)];
    float _S410 = _S395[int(14)];
    float _S411 = _S395[int(15)];
    float _S412 = _S395[int(16)];
    float _S413 = _S395[int(17)];
    float _S414 = _S395[int(18)];
    float _S415 = _S395[int(19)];
    float _S416 = _S395[int(20)];
    float _S417 = _S395[int(21)];
    float _S418 = _S395[int(22)];
    FixedArray<float, 23>  _S419;
    if(_S320)
    {
        float _S420 = _S321 * _S362;
        DiffPair_float_0 _S421;
        (&_S421)->primal_0 = _S322;
        (&_S421)->differential_0 = 0.0f;
        DiffPair_float_0 _S422;
        (&_S422)->primal_0 = 0.0f;
        (&_S422)->differential_0 = 0.0f;
        DiffPair_float_0 _S423;
        (&_S423)->primal_0 = 2.0f;
        (&_S423)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S421, &_S422, &_S423, _S420);
        float _S424 = - _S421.differential_0 / _S323;
        float _S425 = _S324 * - _S424;
        float _S426 = _S325 * _S424;
        DiffPair_float_0 _S427;
        (&_S427)->primal_0 = _S326;
        (&_S427)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S427, _S425);
        DiffPair_float_0 _S428;
        (&_S428)->primal_0 = 9.999999960041972e-13f;
        (&_S428)->differential_0 = 0.0f;
        DiffPair_float_0 _S429;
        (&_S429)->primal_0 = _S327;
        (&_S429)->differential_0 = 0.0f;
        _d_max_0(&_S428, &_S429, _S427.differential_0);
        float _S430 = _S328 * _S429.differential_0;
        float _S431 = _S329 * _S429.differential_0;
        float _S432 = - _S430 / _S330;
        float _S433 = _S332 * (_S315[int(18)] * _S432);
        float _S434 = - _S431 / _S330;
        float _S435 = _S334 * (_S315[int(18)] * _S434);
        float _S436 = - _S426 / _S330;
        float _S437 = _S315[int(18)] * _S436;
        float _S438 = _S433 + _S433 + _S334 * _S437;
        float _S439 = _S435 + _S435 + _S332 * _S437;
        float _S440 = _S331 * - _S432 + _S333 * - _S434 + _S335 * - _S436;
        FixedArray<float, 23>  _S441;
        _S441[int(0)] = 0.0f;
        _S441[int(1)] = 0.0f;
        _S441[int(2)] = 0.0f;
        _S441[int(3)] = 0.0f;
        _S441[int(4)] = 0.0f;
        _S441[int(5)] = 0.0f;
        _S441[int(6)] = 0.0f;
        _S441[int(7)] = 0.0f;
        _S441[int(8)] = 0.0f;
        _S441[int(9)] = 0.0f;
        _S441[int(10)] = 0.0f;
        _S441[int(11)] = 0.0f;
        _S441[int(12)] = 0.0f;
        _S441[int(13)] = 0.0f;
        _S441[int(14)] = 0.0f;
        _S441[int(15)] = 0.0f;
        _S441[int(16)] = 0.0f;
        _S441[int(17)] = 0.0f;
        _S441[int(18)] = 0.0f;
        _S441[int(19)] = 0.0f;
        _S441[int(20)] = 0.0f;
        _S441[int(21)] = 0.0f;
        _S441[int(22)] = 0.0f;
        _S441[int(5)] = _S430;
        _S441[int(4)] = _S431;
        _S441[int(3)] = _S438;
        _S441[int(2)] = _S439;
        _S441[int(6)] = _S426;
        float _S442 = _S397 + _S441[int(1)];
        float _S443 = _S398 + _S441[int(2)];
        float _S444 = _S399 + _S441[int(3)];
        float _S445 = _S400 + _S441[int(4)];
        float _S446 = _S401 + _S441[int(5)];
        float _S447 = _S402 + _S441[int(6)];
        float _S448 = _S403 + _S441[int(7)];
        float _S449 = _S404 + _S441[int(8)];
        float _S450 = _S405 + _S441[int(9)];
        float _S451 = _S406 + _S441[int(10)];
        float _S452 = _S407 + _S441[int(11)];
        float _S453 = _S408 + _S441[int(12)];
        float _S454 = _S409 + _S441[int(13)];
        float _S455 = _S410 + _S441[int(14)];
        float _S456 = _S411 + _S441[int(15)];
        float _S457 = _S412 + _S441[int(16)];
        float _S458 = _S413 + _S441[int(17)];
        float _S459 = _S414 + _S441[int(18)];
        float _S460 = _S415 + _S441[int(19)];
        float _S461 = _S416 + _S441[int(20)];
        float _S462 = _S417 + _S441[int(21)];
        float _S463 = _S418 + _S441[int(22)];
        _S419[int(0)] = _S396 + _S441[int(0)];
        _S419[int(1)] = _S442;
        _S419[int(2)] = _S443;
        _S419[int(3)] = _S444;
        _S419[int(4)] = _S445;
        _S419[int(5)] = _S446;
        _S419[int(6)] = _S447;
        _S419[int(7)] = _S448;
        _S419[int(8)] = _S449;
        _S419[int(9)] = _S450;
        _S419[int(10)] = _S451;
        _S419[int(11)] = _S452;
        _S419[int(12)] = _S453;
        _S419[int(13)] = _S454;
        _S419[int(14)] = _S455;
        _S419[int(15)] = _S456;
        _S419[int(16)] = _S457;
        _S419[int(17)] = _S458;
        _S419[int(18)] = _S459;
        _S419[int(19)] = _S460;
        _S419[int(20)] = _S461;
        _S419[int(21)] = _S462;
        _S419[int(22)] = _S463;
        _S321 = _S440;
    }
    else
    {
        _S419[int(0)] = _S396;
        _S419[int(1)] = _S397;
        _S419[int(2)] = _S398;
        _S419[int(3)] = _S399;
        _S419[int(4)] = _S400;
        _S419[int(5)] = _S401;
        _S419[int(6)] = _S402;
        _S419[int(7)] = _S403;
        _S419[int(8)] = _S404;
        _S419[int(9)] = _S405;
        _S419[int(10)] = _S406;
        _S419[int(11)] = _S407;
        _S419[int(12)] = _S408;
        _S419[int(13)] = _S409;
        _S419[int(14)] = _S410;
        _S419[int(15)] = _S411;
        _S419[int(16)] = _S412;
        _S419[int(17)] = _S413;
        _S419[int(18)] = _S414;
        _S419[int(19)] = _S415;
        _S419[int(20)] = _S416;
        _S419[int(21)] = _S417;
        _S419[int(22)] = _S418;
        _S321 = 0.0f;
    }
    if(_S319)
    {
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
        _S464[int(3)] = 0.0f;
        float _S465 = _S419[int(1)] + _S464[int(1)];
        float _S466 = _S419[int(2)] + _S464[int(2)];
        float _S467 = _S419[int(3)] + _S464[int(3)];
        float _S468 = _S419[int(4)] + _S464[int(4)];
        float _S469 = _S419[int(5)] + _S464[int(5)];
        float _S470 = _S419[int(6)] + _S464[int(6)];
        float _S471 = _S419[int(7)] + _S464[int(7)];
        float _S472 = _S419[int(8)] + _S464[int(8)];
        float _S473 = _S419[int(9)] + _S464[int(9)];
        float _S474 = _S419[int(10)] + _S464[int(10)];
        float _S475 = _S419[int(11)] + _S464[int(11)];
        float _S476 = _S419[int(12)] + _S464[int(12)];
        float _S477 = _S419[int(13)] + _S464[int(13)];
        float _S478 = _S419[int(14)] + _S464[int(14)];
        float _S479 = _S419[int(15)] + _S464[int(15)];
        float _S480 = _S419[int(16)] + _S464[int(16)];
        float _S481 = _S419[int(17)] + _S464[int(17)];
        float _S482 = _S419[int(18)] + _S464[int(18)];
        float _S483 = _S419[int(19)] + _S464[int(19)];
        float _S484 = _S419[int(20)] + _S464[int(20)];
        float _S485 = _S419[int(21)] + _S464[int(21)];
        float _S486 = _S419[int(22)] + _S464[int(22)];
        _S419[int(0)] = _S419[int(0)] + _S464[int(0)];
        _S419[int(1)] = _S465;
        _S419[int(2)] = _S466;
        _S419[int(3)] = _S467;
        _S419[int(4)] = _S468;
        _S419[int(5)] = _S469;
        _S419[int(6)] = _S470;
        _S419[int(7)] = _S471;
        _S419[int(8)] = _S472;
        _S419[int(9)] = _S473;
        _S419[int(10)] = _S474;
        _S419[int(11)] = _S475;
        _S419[int(12)] = _S476;
        _S419[int(13)] = _S477;
        _S419[int(14)] = _S478;
        _S419[int(15)] = _S479;
        _S419[int(16)] = _S480;
        _S419[int(17)] = _S481;
        _S419[int(18)] = _S482;
        _S419[int(19)] = _S483;
        _S419[int(20)] = _S484;
        _S419[int(21)] = _S485;
        _S419[int(22)] = _S486;
    }
    float _S487 = -10.0f * _S361;
    DiffPair_float_0 _S488;
    (&_S488)->primal_0 = _S318;
    (&_S488)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S488, _S487);
    float _S489 = _S488.differential_0 / _S317;
    float _S490 = _S316 * _S489;
    float _S491 = _S360 / _S317;
    float _S492 = _S316 * _S491;
    float _S493 = _S315[int(1)] * - _S489 + _S315[int(0)] * - _S491;
    DiffPair_float_0 _S494;
    (&_S494)->primal_0 = _S315[int(17)];
    (&_S494)->differential_0 = 0.0f;
    DiffPair_float_0 _S495;
    (&_S495)->primal_0 = 1.0f;
    (&_S495)->differential_0 = 0.0f;
    _d_max_0(&_S494, &_S495, _S493);
    FixedArray<float, 23>  _S496;
    _S496[int(0)] = 0.0f;
    _S496[int(1)] = 0.0f;
    _S496[int(2)] = 0.0f;
    _S496[int(3)] = 0.0f;
    _S496[int(4)] = 0.0f;
    _S496[int(5)] = 0.0f;
    _S496[int(6)] = 0.0f;
    _S496[int(7)] = 0.0f;
    _S496[int(8)] = 0.0f;
    _S496[int(9)] = 0.0f;
    _S496[int(10)] = 0.0f;
    _S496[int(11)] = 0.0f;
    _S496[int(12)] = 0.0f;
    _S496[int(13)] = 0.0f;
    _S496[int(14)] = 0.0f;
    _S496[int(15)] = 0.0f;
    _S496[int(16)] = 0.0f;
    _S496[int(17)] = 0.0f;
    _S496[int(18)] = 0.0f;
    _S496[int(19)] = 0.0f;
    _S496[int(20)] = 0.0f;
    _S496[int(21)] = 0.0f;
    _S496[int(22)] = 0.0f;
    _S496[int(18)] = _S321;
    _S496[int(1)] = _S490;
    _S496[int(17)] = _S494.differential_0;
    _S496[int(0)] = _S492;
    FixedArray<float, 23>  _S497 = {
        _S419[int(0)] + _S496[int(0)], _S419[int(1)] + _S496[int(1)], _S419[int(2)] + _S496[int(2)], _S419[int(3)] + _S496[int(3)], _S419[int(4)] + _S496[int(4)], _S419[int(5)] + _S496[int(5)], _S419[int(6)] + _S496[int(6)], _S419[int(7)] + _S496[int(7)], _S419[int(8)] + _S496[int(8)], _S419[int(9)] + _S496[int(9)], _S419[int(10)] + _S496[int(10)], _S419[int(11)] + _S496[int(11)], _S419[int(12)] + _S496[int(12)], _S419[int(13)] + _S496[int(13)], _S419[int(14)] + _S496[int(14)], _S419[int(15)] + _S496[int(15)], _S419[int(16)] + _S496[int(16)], _S419[int(17)] + _S496[int(17)], _S419[int(18)] + _S496[int(18)], _S419[int(19)] + _S496[int(19)], _S419[int(20)] + _S496[int(20)], _S419[int(21)] + _S496[int(21)], _S419[int(22)] + _S496[int(22)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S497;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S498, FixedArray<float, 15>  * _S499, FixedArray<float, 10>  * _S500)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S498, _S499, _S500);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 23>  raw_losses_1, FixedArray<float, 15>  weights_5, FixedArray<float, 10>  v_losses_1, FixedArray<float, 23>  * _S501)
{
    FixedArray<float, 23>  _S502 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S502;
    FixedArray<float, 15>  _S503 = weights_5;
    FixedArray<float, 10>  _S504 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S503, &_S504);
    *_S501 = (&dp_raw_losses_0)->differential_0;
    return;
}

