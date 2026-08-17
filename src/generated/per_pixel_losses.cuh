#pragma once

#include "generated/slang.cuh"

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

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_7, float dOut_7)
{
    float _S13 = 1.0f / (*dpx_7).primal_0 * dOut_7;
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

inline __device__ void per_pixel_losses(float3  render_rgb_0, float3  ref_rgb_0, float render_depth_0, float ref_depth_0, float3  render_normal_0, float3  depth_normal_0, float3  ref_normal_0, float render_Ts_0, float3  rgb_dist_0, float depth_dist_0, float3  normal_dist_0, float median_depth_0, float3  median_normal_0, bool ref_alpha_0, bool has_mask_0, FixedArray<float, 19>  weights_0, FixedArray<float, 31>  * _S18)
{
    bool _S19;
    bool _S20;
    bool _S21;
    FixedArray<float, 31>  losses_0;
    bool mask_0;
    if(has_mask_0)
    {
        mask_0 = ref_alpha_0;
    }
    else
    {
        mask_0 = true;
    }
    float3  _S22;
    bool depth_mask_0 = ref_depth_0 != 0.0f;
    bool normal_mask_0 = (ref_normal_0.x + ref_normal_0.y + ref_normal_0.z) > -2.36599993705749512f;
    float _S23 = render_rgb_0.x;
    float _S24 = render_rgb_0.y;
    float _S25 = render_rgb_0.z;
    float _S26 = ref_rgb_0.x;
    float _S27 = ref_rgb_0.y;
    float _S28 = ref_rgb_0.z;
    float dY_0 = 0.29899999499320984f * _S23 + 0.58700001239776611f * _S24 + 0.11400000005960464f * _S25 - (0.29899999499320984f * _S26 + 0.58700001239776611f * _S27 + 0.11400000005960464f * _S28);
    float dU_0 = -0.14712999761104584f * _S23 - 0.28885999321937561f * _S24 + 0.43599998950958252f * _S25 - (-0.14712999761104584f * _S26 - 0.28885999321937561f * _S27 + 0.43599998950958252f * _S28);
    float dV_0 = 0.61500000953674316f * _S23 - 0.51498997211456299f * _S24 - 0.10001000016927719f * _S25 - (0.61500000953674316f * _S26 - 0.51498997211456299f * _S27 - 0.10001000016927719f * _S28);
    float _S29 = float(mask_0);
    float3  _S30 = ref_rgb_0 - render_rgb_0;
    float3  _S31 = abs_0(_S30);
    float _S32 = dot_0(_S30, _S30) * 0.3333333432674408f;
    losses_0[int(0)] = _S29 * (weights_0[int(0)] * ((_S31.x + _S31.y + _S31.z) * 0.3333333432674408f) + weights_0[int(1)] * _S32 + weights_0[int(2)] * (F32_abs((dY_0))) + weights_0[int(3)] * dY_0 * dY_0 + weights_0[int(4)] * dU_0 * dU_0 + weights_0[int(5)] * dV_0 * dV_0);
    losses_0[int(1)] = _S29 * clamp_0(_S32, 0.0f, 1.0f);
    float _S33 = float(depth_mask_0 & mask_0);
    float _S34 = _S33 * (F32_max((render_depth_0), (0.00009999999747379f)));
    float _S35 = _S33 * (F32_max((ref_depth_0), (0.00009999999747379f)));
    losses_0[int(2)] = _S34;
    losses_0[int(3)] = _S35;
    losses_0[int(4)] = _S34 * _S34;
    losses_0[int(5)] = _S35 * _S35;
    losses_0[int(6)] = _S34 * _S35;
    bool _S36 = normal_mask_0 & mask_0;
    for(;;)
    {
        float norm2_0 = dot_0(render_normal_0, render_normal_0);
        bool _S37 = norm2_0 == 0.0f;
        _S19 = _S37;
        if(_S37)
        {
            _S22 = make_float3 (0.0f);
            break;
        }
        _S22 = render_normal_0 * make_float3 ((F32_rsqrt((norm2_0))));
        break;
    }
    float3  _S38;
    bool _S39 = !_S19;
    for(;;)
    {
        float norm2_1 = dot_0(depth_normal_0, depth_normal_0);
        bool _S40 = norm2_1 == 0.0f;
        _S20 = _S40;
        if(_S40)
        {
            _S38 = make_float3 (0.0f);
            break;
        }
        _S38 = depth_normal_0 * make_float3 ((F32_rsqrt((norm2_1))));
        break;
    }
    float3  _S41;
    bool normal_mask_1;
    bool _S42 = !_S20;
    for(;;)
    {
        float norm2_2 = dot_0(ref_normal_0, ref_normal_0);
        if(norm2_2 == 0.0f)
        {
            _S41 = make_float3 (0.0f);
            normal_mask_1 = false;
            break;
        }
        _S41 = ref_normal_0 * make_float3 ((F32_rsqrt((norm2_2))));
        normal_mask_1 = _S36;
        break;
    }
    float3  _S43;
    float _S44 = float(_S39 & normal_mask_1);
    float cos_sim_loss_0 = 0.5f - 0.5f * dot_0(_S22, _S41);
    losses_0[int(7)] = weights_0[int(7)] * _S44 * (cos_sim_loss_0 + (F32_sqrt(((F32_max((cos_sim_loss_0), (9.999999960041972e-13f)))))));
    float _S45 = float(_S42 & normal_mask_1);
    float cos_sim_loss_1 = 0.5f - 0.5f * dot_0(_S38, _S41);
    losses_0[int(8)] = weights_0[int(7)] * _S45 * (cos_sim_loss_1 + (F32_sqrt(((F32_max((cos_sim_loss_1), (9.999999960041972e-13f)))))));
    float _S46 = float(_S39 & _S42);
    float cos_sim_loss_2 = 0.5f - 0.5f * dot_0(_S22, _S38);
    losses_0[int(11)] = weights_0[int(10)] * _S46 * (cos_sim_loss_2 + (F32_sqrt(((F32_max((cos_sim_loss_2), (9.999999960041972e-13f)))))));
    for(;;)
    {
        float norm2_3 = dot_0(median_normal_0, median_normal_0);
        bool _S47 = norm2_3 == 0.0f;
        _S21 = _S47;
        if(_S47)
        {
            _S43 = make_float3 (0.0f);
            break;
        }
        _S43 = median_normal_0 * make_float3 ((F32_rsqrt((norm2_3))));
        break;
    }
    bool _S48 = !_S21;
    if(mask_0)
    {
        mask_0 = render_depth_0 > 1.00000001335143196e-10f;
    }
    else
    {
        mask_0 = false;
    }
    bool mean_median_mask_0;
    if(mask_0)
    {
        mean_median_mask_0 = median_depth_0 > 1.00000001335143196e-10f;
    }
    else
    {
        mean_median_mask_0 = false;
    }
    float _S49 = float(mean_median_mask_0);
    losses_0[int(16)] = weights_0[int(15)] * _S49 * (F32_abs(((F32_log((render_depth_0))) - (F32_log((median_depth_0))))));
    float _S50 = float(_S48 & _S42);
    float cos_sim_loss_3 = 0.5f - 0.5f * dot_0(_S43, _S38);
    losses_0[int(17)] = weights_0[int(16)] * _S50 * (cos_sim_loss_3 + (F32_sqrt(((F32_max((cos_sim_loss_3), (9.999999960041972e-13f)))))));
    float _S51 = float(_S48 & normal_mask_1);
    float cos_sim_loss_4 = 0.5f - 0.5f * dot_0(_S43, _S41);
    losses_0[int(18)] = weights_0[int(17)] * _S51 * (cos_sim_loss_4 + (F32_sqrt(((F32_max((cos_sim_loss_4), (9.999999960041972e-13f)))))));
    float _S52 = float(_S48 & _S39);
    float cos_sim_loss_5 = 0.5f - 0.5f * dot_0(_S43, _S22);
    losses_0[int(19)] = weights_0[int(18)] * _S52 * (cos_sim_loss_5 + (F32_sqrt(((F32_max((cos_sim_loss_5), (9.999999960041972e-13f)))))));
    float render_alpha_0 = clamp_0(1.0f - render_Ts_0, 0.0f, 1.0f);
    float _S53 = float(depth_mask_0);
    float _S54 = float(ref_alpha_0);
    float _S55 = (F32_max((render_alpha_0), (_S54)));
    losses_0[int(9)] = weights_0[int(8)] * _S53 * - lerp_0((F32_log(((F32_max((1.0f - _S55), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S55), (9.99999997475242708e-07f)))))), _S54);
    float _S56 = 1.0f - render_alpha_0;
    float _S57 = 1.0f - _S54;
    float _S58 = (F32_max((_S56), (_S57)));
    losses_0[int(10)] = weights_0[int(9)] * _S53 * - lerp_0((F32_log(((F32_max((1.0f - _S58), (9.99999997475242708e-07f)))))), (F32_log(((F32_max((_S58), (9.99999997475242708e-07f)))))), _S57);
    losses_0[int(12)] = weights_0[int(11)] * 4.0f * render_alpha_0 * _S56;
    losses_0[int(13)] = weights_0[int(12)] * ((rgb_dist_0.x + rgb_dist_0.y + rgb_dist_0.z) * 0.3333333432674408f);
    losses_0[int(14)] = weights_0[int(13)] * depth_dist_0;
    losses_0[int(15)] = weights_0[int(14)] * ((normal_dist_0.x + normal_dist_0.y + normal_dist_0.z) * 0.3333333432674408f);
    losses_0[int(20)] = 1.0f;
    losses_0[int(21)] = _S29;
    losses_0[int(22)] = _S33;
    losses_0[int(23)] = _S44;
    losses_0[int(24)] = _S45;
    losses_0[int(25)] = _S46;
    losses_0[int(26)] = _S53;
    losses_0[int(27)] = _S49;
    losses_0[int(28)] = _S50;
    losses_0[int(29)] = _S51;
    losses_0[int(30)] = _S52;
    *_S18 = losses_0;
    return;
}

inline __device__ float s_primal_ctx_dot_0(float3  _S59, float3  _S60)
{
    return dot_0(_S59, _S60);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S61)
{
    return (F32_rsqrt((_S61)));
}

inline __device__ float s_primal_ctx_log_0(float _S62)
{
    return (F32_log((_S62)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S63, float _S64, float _S65)
{
    return clamp_0(_S63, _S64, _S65);
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S66, DiffPair_float_0 * _S67, DiffPair_float_0 * _S68, float _S69)
{
    _d_lerp_0(_S66, _S67, _S68, _S69);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S70, float _S71)
{
    _d_log_0(_S70, _S71);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S72, DiffPair_float_0 * _S73, DiffPair_float_0 * _S74, float _S75)
{
    _d_clamp_0(_S72, _S73, _S74, _S75);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S76, float _S77)
{
    _d_sqrt_0(_S76, _S77);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S78, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S79, float _S80)
{
    _d_dot_0(_S78, _S79, _S80);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S81, float _S82)
{
    _d_abs_0(_S81, _S82);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S83, float _S84)
{
    _d_rsqrt_0(_S83, _S84);
    return;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S85, float3  _S86)
{
    _d_abs_vector_0(_S85, _S86);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_rgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_rgb_0, DiffPair_float_0 * dprender_depth_0, DiffPair_float_0 * dpref_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprender_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepth_normal_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpref_normal_0, DiffPair_float_0 * dprender_Ts_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_dist_0, DiffPair_float_0 * dpdepth_dist_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpnormal_dist_0, DiffPair_float_0 * dpmedian_depth_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmedian_normal_0, bool ref_alpha_1, bool has_mask_1, FixedArray<float, 19>  * weights_1, FixedArray<float, 31>  * _s_dOut_0)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S87 = *dprender_rgb_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S88 = *dpref_rgb_0;
    DiffPair_float_0 _S89 = *dprender_depth_0;
    DiffPair_float_0 _S90 = *dpref_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S91 = *dprender_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S92 = *dpdepth_normal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S93 = *dpref_normal_0;
    DiffPair_float_0 _S94 = *dprender_Ts_0;
    DiffPair_float_0 _S95 = *dpmedian_depth_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S96 = *dpmedian_normal_0;
    float3  _S97 = make_float3 (0.0f);
    bool mask_1;
    if(has_mask_1)
    {
        mask_1 = ref_alpha_1;
    }
    else
    {
        mask_1 = true;
    }
    bool depth_mask_1 = (_S90.primal_0) != 0.0f;
    float _S98 = _S87.primal_0.x;
    float _S99 = _S87.primal_0.y;
    float _S100 = _S87.primal_0.z;
    float _S101 = _S88.primal_0.x;
    float _S102 = _S88.primal_0.y;
    float _S103 = _S88.primal_0.z;
    float dY_1 = 0.29899999499320984f * _S98 + 0.58700001239776611f * _S99 + 0.11400000005960464f * _S100 - (0.29899999499320984f * _S101 + 0.58700001239776611f * _S102 + 0.11400000005960464f * _S103);
    float dU_1 = -0.14712999761104584f * _S98 - 0.28885999321937561f * _S99 + 0.43599998950958252f * _S100 - (-0.14712999761104584f * _S101 - 0.28885999321937561f * _S102 + 0.43599998950958252f * _S103);
    float dV_1 = 0.61500000953674316f * _S98 - 0.51498997211456299f * _S99 - 0.10001000016927719f * _S100 - (0.61500000953674316f * _S101 - 0.51498997211456299f * _S102 - 0.10001000016927719f * _S103);
    float _S104 = float(mask_1);
    float _S105 = (*weights_1)[int(0)];
    float3  _S106 = _S88.primal_0 - _S87.primal_0;
    float _S107 = (*weights_1)[int(1)];
    float _S108 = s_primal_ctx_dot_0(_S106, _S106) * 0.3333333432674408f;
    float _S109 = (*weights_1)[int(2)];
    float _S110 = (*weights_1)[int(3)];
    float _S111 = (*weights_1)[int(3)] * dY_1;
    float _S112 = (*weights_1)[int(4)];
    float _S113 = (*weights_1)[int(4)] * dU_1;
    float _S114 = (*weights_1)[int(5)];
    float _S115 = (*weights_1)[int(5)] * dV_1;
    float _S116 = float(depth_mask_1 & mask_1);
    float _S117 = _S116 * (F32_max((_S89.primal_0), (0.00009999999747379f)));
    float _S118 = _S116 * (F32_max((_S90.primal_0), (0.00009999999747379f)));
    bool _S119 = ((_S93.primal_0.x + _S93.primal_0.y + _S93.primal_0.z) > -2.36599993705749512f) & mask_1;
    float _S120 = s_primal_ctx_dot_0(_S91.primal_0, _S91.primal_0);
    bool _S121 = _S120 == 0.0f;
    float3  _S122;
    if(_S121)
    {
        _S122 = make_float3 (0.0f);
    }
    bool _S123 = !_S121;
    float3  _S124;
    if(_S123)
    {
        float _S125 = s_primal_ctx_rsqrt_0(_S120);
        float3  _S126 = make_float3 (_S125);
        _S122 = _S91.primal_0 * make_float3 (_S125);
        _S124 = _S126;
    }
    else
    {
        _S124 = _S97;
    }
    float _S127 = s_primal_ctx_dot_0(_S92.primal_0, _S92.primal_0);
    bool _S128 = _S127 == 0.0f;
    float3  _S129;
    if(_S128)
    {
        _S129 = make_float3 (0.0f);
    }
    bool _S130 = !_S128;
    float3  _S131;
    if(_S130)
    {
        float _S132 = s_primal_ctx_rsqrt_0(_S127);
        float3  _S133 = make_float3 (_S132);
        _S129 = _S92.primal_0 * make_float3 (_S132);
        _S131 = _S133;
    }
    else
    {
        _S131 = _S97;
    }
    float _S134 = s_primal_ctx_dot_0(_S93.primal_0, _S93.primal_0);
    bool _S135 = _S134 == 0.0f;
    bool normal_mask_2;
    float3  _S136;
    if(_S135)
    {
        float3  _S137 = make_float3 (0.0f);
        normal_mask_2 = false;
        _S136 = _S137;
    }
    else
    {
        normal_mask_2 = _S119;
    }
    bool _S138 = !_S135;
    float3  _S139;
    if(_S138)
    {
        float _S140 = s_primal_ctx_rsqrt_0(_S134);
        float3  _S141 = make_float3 (_S140);
        _S136 = _S93.primal_0 * make_float3 (_S140);
        _S139 = _S141;
    }
    else
    {
        _S139 = _S97;
    }
    float _S142 = (*weights_1)[int(7)] * float(_S123 & normal_mask_2);
    float cos_sim_loss_6 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S122, _S136);
    float _S143 = (F32_max((cos_sim_loss_6), (9.999999960041972e-13f)));
    float _S144 = (*weights_1)[int(7)] * float(_S130 & normal_mask_2);
    float cos_sim_loss_7 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S129, _S136);
    float _S145 = (F32_max((cos_sim_loss_7), (9.999999960041972e-13f)));
    float _S146 = (*weights_1)[int(10)] * float(_S123 & _S130);
    float cos_sim_loss_8 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S122, _S129);
    float _S147 = (F32_max((cos_sim_loss_8), (9.999999960041972e-13f)));
    float _S148 = s_primal_ctx_dot_0(_S96.primal_0, _S96.primal_0);
    bool _S149 = _S148 == 0.0f;
    float3  _S150;
    if(_S149)
    {
        _S150 = make_float3 (0.0f);
    }
    bool _S151 = !_S149;
    float3  _S152;
    if(_S151)
    {
        float _S153 = s_primal_ctx_rsqrt_0(_S148);
        float3  _S154 = make_float3 (_S153);
        _S150 = _S96.primal_0 * make_float3 (_S153);
        _S152 = _S154;
    }
    else
    {
        _S152 = _S97;
    }
    if(mask_1)
    {
        mask_1 = (_S89.primal_0) > 1.00000001335143196e-10f;
    }
    else
    {
        mask_1 = false;
    }
    bool mean_median_mask_1;
    if(mask_1)
    {
        mean_median_mask_1 = (_S95.primal_0) > 1.00000001335143196e-10f;
    }
    else
    {
        mean_median_mask_1 = false;
    }
    float _S155 = (*weights_1)[int(15)] * float(mean_median_mask_1);
    float _S156 = s_primal_ctx_log_0(_S89.primal_0) - s_primal_ctx_log_0(_S95.primal_0);
    float _S157 = (*weights_1)[int(16)] * float(_S151 & _S130);
    float cos_sim_loss_9 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S150, _S129);
    float _S158 = (F32_max((cos_sim_loss_9), (9.999999960041972e-13f)));
    float _S159 = (*weights_1)[int(17)] * float(_S151 & normal_mask_2);
    float cos_sim_loss_10 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S150, _S136);
    float _S160 = (F32_max((cos_sim_loss_10), (9.999999960041972e-13f)));
    float _S161 = (*weights_1)[int(18)] * float(_S151 & _S123);
    float cos_sim_loss_11 = 0.5f - 0.5f * s_primal_ctx_dot_0(_S150, _S122);
    float _S162 = (F32_max((cos_sim_loss_11), (9.999999960041972e-13f)));
    float _S163 = 1.0f - _S94.primal_0;
    float _S164 = s_primal_ctx_clamp_0(_S163, 0.0f, 1.0f);
    float _S165 = float(depth_mask_1);
    float _S166 = (*weights_1)[int(8)] * _S165;
    float _S167 = float(ref_alpha_1);
    float _S168 = (F32_max((_S164), (_S167)));
    float _S169 = 1.0f - _S168;
    float _S170 = (F32_max((_S169), (9.99999997475242708e-07f)));
    float _S171 = s_primal_ctx_log_0(_S170);
    float _S172 = (F32_max((_S168), (9.99999997475242708e-07f)));
    float _S173 = s_primal_ctx_log_0(_S172);
    float _S174 = (*weights_1)[int(9)] * _S165;
    float _S175 = 1.0f - _S164;
    float _S176 = 1.0f - _S167;
    float _S177 = (F32_max((_S175), (_S176)));
    float _S178 = 1.0f - _S177;
    float _S179 = (F32_max((_S178), (9.99999997475242708e-07f)));
    float _S180 = s_primal_ctx_log_0(_S179);
    float _S181 = (F32_max((_S177), (9.99999997475242708e-07f)));
    float _S182 = s_primal_ctx_log_0(_S181);
    float _S183 = (*weights_1)[int(11)] * 4.0f;
    float _S184 = _S183 * _S164;
    float _S185 = (*weights_1)[int(12)];
    float _S186 = (*weights_1)[int(13)];
    float _S187 = (*weights_1)[int(14)];
    float _S188 = (*_s_dOut_0)[int(0)];
    float _S189 = (*_s_dOut_0)[int(1)];
    float _S190 = (*_s_dOut_0)[int(2)];
    float _S191 = (*_s_dOut_0)[int(3)];
    float _S192 = (*_s_dOut_0)[int(4)];
    float _S193 = (*_s_dOut_0)[int(5)];
    float _S194 = (*_s_dOut_0)[int(6)];
    float _S195 = (*_s_dOut_0)[int(7)];
    float _S196 = (*_s_dOut_0)[int(8)];
    float _S197 = (*_s_dOut_0)[int(11)];
    float _S198 = 0.3333333432674408f * (_S187 * (*_s_dOut_0)[int(15)]);
    float _S199 = _S186 * (*_s_dOut_0)[int(14)];
    float _S200 = 0.3333333432674408f * (_S185 * (*_s_dOut_0)[int(13)]);
    float _S201 = _S184 * (*_s_dOut_0)[int(12)];
    float _S202 = _S183 * (_S175 * (*_s_dOut_0)[int(12)]);
    float _S203 = - (_S174 * (*_s_dOut_0)[int(10)]);
    DiffPair_float_0 _S204;
    (&_S204)->primal_0 = _S180;
    (&_S204)->differential_0 = 0.0f;
    DiffPair_float_0 _S205;
    (&_S205)->primal_0 = _S182;
    (&_S205)->differential_0 = 0.0f;
    DiffPair_float_0 _S206;
    (&_S206)->primal_0 = _S176;
    (&_S206)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S204, &_S205, &_S206, _S203);
    DiffPair_float_0 _S207;
    (&_S207)->primal_0 = _S181;
    (&_S207)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S207, _S205.differential_0);
    DiffPair_float_0 _S208;
    (&_S208)->primal_0 = _S177;
    (&_S208)->differential_0 = 0.0f;
    DiffPair_float_0 _S209;
    (&_S209)->primal_0 = 9.99999997475242708e-07f;
    (&_S209)->differential_0 = 0.0f;
    _d_max_0(&_S208, &_S209, _S207.differential_0);
    DiffPair_float_0 _S210;
    (&_S210)->primal_0 = _S179;
    (&_S210)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S210, _S204.differential_0);
    DiffPair_float_0 _S211;
    (&_S211)->primal_0 = _S178;
    (&_S211)->differential_0 = 0.0f;
    DiffPair_float_0 _S212;
    (&_S212)->primal_0 = 9.99999997475242708e-07f;
    (&_S212)->differential_0 = 0.0f;
    _d_max_0(&_S211, &_S212, _S210.differential_0);
    float _S213 = _S208.differential_0 + - _S211.differential_0;
    DiffPair_float_0 _S214;
    (&_S214)->primal_0 = _S175;
    (&_S214)->differential_0 = 0.0f;
    DiffPair_float_0 _S215;
    (&_S215)->primal_0 = _S176;
    (&_S215)->differential_0 = 0.0f;
    _d_max_0(&_S214, &_S215, _S213);
    float _S216 = - (_S201 + _S214.differential_0);
    float _S217 = - (_S166 * (*_s_dOut_0)[int(9)]);
    DiffPair_float_0 _S218;
    (&_S218)->primal_0 = _S171;
    (&_S218)->differential_0 = 0.0f;
    DiffPair_float_0 _S219;
    (&_S219)->primal_0 = _S173;
    (&_S219)->differential_0 = 0.0f;
    DiffPair_float_0 _S220;
    (&_S220)->primal_0 = _S167;
    (&_S220)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S218, &_S219, &_S220, _S217);
    DiffPair_float_0 _S221;
    (&_S221)->primal_0 = _S172;
    (&_S221)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S221, _S219.differential_0);
    DiffPair_float_0 _S222;
    (&_S222)->primal_0 = _S168;
    (&_S222)->differential_0 = 0.0f;
    DiffPair_float_0 _S223;
    (&_S223)->primal_0 = 9.99999997475242708e-07f;
    (&_S223)->differential_0 = 0.0f;
    _d_max_0(&_S222, &_S223, _S221.differential_0);
    DiffPair_float_0 _S224;
    (&_S224)->primal_0 = _S170;
    (&_S224)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S224, _S218.differential_0);
    DiffPair_float_0 _S225;
    (&_S225)->primal_0 = _S169;
    (&_S225)->differential_0 = 0.0f;
    DiffPair_float_0 _S226;
    (&_S226)->primal_0 = 9.99999997475242708e-07f;
    (&_S226)->differential_0 = 0.0f;
    _d_max_0(&_S225, &_S226, _S224.differential_0);
    float _S227 = _S222.differential_0 + - _S225.differential_0;
    DiffPair_float_0 _S228;
    (&_S228)->primal_0 = _S164;
    (&_S228)->differential_0 = 0.0f;
    DiffPair_float_0 _S229;
    (&_S229)->primal_0 = _S167;
    (&_S229)->differential_0 = 0.0f;
    _d_max_0(&_S228, &_S229, _S227);
    float _S230 = _S202 + _S216 + _S228.differential_0;
    DiffPair_float_0 _S231;
    (&_S231)->primal_0 = _S163;
    (&_S231)->differential_0 = 0.0f;
    DiffPair_float_0 _S232;
    (&_S232)->primal_0 = 0.0f;
    (&_S232)->differential_0 = 0.0f;
    DiffPair_float_0 _S233;
    (&_S233)->primal_0 = 1.0f;
    (&_S233)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S231, &_S232, &_S233, _S230);
    float _S234 = - _S231.differential_0;
    float _S235 = _S161 * (*_s_dOut_0)[int(19)];
    DiffPair_float_0 _S236;
    (&_S236)->primal_0 = _S162;
    (&_S236)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S236, _S235);
    DiffPair_float_0 _S237;
    (&_S237)->primal_0 = cos_sim_loss_11;
    (&_S237)->differential_0 = 0.0f;
    DiffPair_float_0 _S238;
    (&_S238)->primal_0 = 9.999999960041972e-13f;
    (&_S238)->differential_0 = 0.0f;
    _d_max_0(&_S237, &_S238, _S236.differential_0);
    float _S239 = 0.5f * - (_S235 + _S237.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S240;
    (&_S240)->primal_0 = _S150;
    (&_S240)->differential_0 = _S97;
    float3  _S241 = _S122;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S122;
    (&_S242)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S240, &_S242, _S239);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S243 = _S242;
    float _S244 = _S159 * (*_s_dOut_0)[int(18)];
    DiffPair_float_0 _S245;
    (&_S245)->primal_0 = _S160;
    (&_S245)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S245, _S244);
    DiffPair_float_0 _S246;
    (&_S246)->primal_0 = cos_sim_loss_10;
    (&_S246)->differential_0 = 0.0f;
    DiffPair_float_0 _S247;
    (&_S247)->primal_0 = 9.999999960041972e-13f;
    (&_S247)->differential_0 = 0.0f;
    _d_max_0(&_S246, &_S247, _S245.differential_0);
    float _S248 = 0.5f * - (_S244 + _S246.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S249;
    (&_S249)->primal_0 = _S150;
    (&_S249)->differential_0 = _S97;
    float3  _S250 = _S136;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S251;
    (&_S251)->primal_0 = _S136;
    (&_S251)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S249, &_S251, _S248);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S252 = _S251;
    float _S253 = _S157 * (*_s_dOut_0)[int(17)];
    DiffPair_float_0 _S254;
    (&_S254)->primal_0 = _S158;
    (&_S254)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S254, _S253);
    DiffPair_float_0 _S255;
    (&_S255)->primal_0 = cos_sim_loss_9;
    (&_S255)->differential_0 = 0.0f;
    DiffPair_float_0 _S256;
    (&_S256)->primal_0 = 9.999999960041972e-13f;
    (&_S256)->differential_0 = 0.0f;
    _d_max_0(&_S255, &_S256, _S254.differential_0);
    float _S257 = 0.5f * - (_S253 + _S255.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S258;
    (&_S258)->primal_0 = _S150;
    (&_S258)->differential_0 = _S97;
    float3  _S259 = _S129;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S260;
    (&_S260)->primal_0 = _S129;
    (&_S260)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S258, &_S260, _S257);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S261 = _S260;
    float _S262 = _S155 * (*_s_dOut_0)[int(16)];
    DiffPair_float_0 _S263;
    (&_S263)->primal_0 = _S156;
    (&_S263)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S263, _S262);
    float _S264 = - _S263.differential_0;
    DiffPair_float_0 _S265;
    (&_S265)->primal_0 = _S95.primal_0;
    (&_S265)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S265, _S264);
    DiffPair_float_0 _S266 = _S265;
    DiffPair_float_0 _S267;
    (&_S267)->primal_0 = _S89.primal_0;
    (&_S267)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S267, _S263.differential_0);
    DiffPair_float_0 _S268 = _S267;
    float3  _S269 = _S240.differential_0 + _S249.differential_0 + _S258.differential_0;
    float3  _S270 = make_float3 (_S200, _S200, _S200);
    float3  _S271 = make_float3 (_S198, _S198, _S198);
    float _S272;
    if(_S151)
    {
        float3  _S273 = _S96.primal_0 * _S269;
        float3  _S274 = _S152 * _S269;
        float _S275 = _S273.x + _S273.y + _S273.z;
        DiffPair_float_0 _S276;
        (&_S276)->primal_0 = _S148;
        (&_S276)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S276, _S275);
        _S272 = _S276.differential_0;
        _S122 = _S274;
    }
    else
    {
        _S272 = 0.0f;
        _S122 = _S97;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S277;
    (&_S277)->primal_0 = _S96.primal_0;
    (&_S277)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S278;
    (&_S278)->primal_0 = _S96.primal_0;
    (&_S278)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S277, &_S278, _S272);
    float _S279 = _S146 * _S197;
    DiffPair_float_0 _S280;
    (&_S280)->primal_0 = _S147;
    (&_S280)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S280, _S279);
    DiffPair_float_0 _S281;
    (&_S281)->primal_0 = cos_sim_loss_8;
    (&_S281)->differential_0 = 0.0f;
    DiffPair_float_0 _S282;
    (&_S282)->primal_0 = 9.999999960041972e-13f;
    (&_S282)->differential_0 = 0.0f;
    _d_max_0(&_S281, &_S282, _S280.differential_0);
    float _S283 = 0.5f * - (_S279 + _S281.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S284;
    (&_S284)->primal_0 = _S241;
    (&_S284)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S285;
    (&_S285)->primal_0 = _S259;
    (&_S285)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S284, &_S285, _S283);
    float _S286 = _S144 * _S196;
    DiffPair_float_0 _S287;
    (&_S287)->primal_0 = _S145;
    (&_S287)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S287, _S286);
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = cos_sim_loss_7;
    (&_S288)->differential_0 = 0.0f;
    DiffPair_float_0 _S289;
    (&_S289)->primal_0 = 9.999999960041972e-13f;
    (&_S289)->differential_0 = 0.0f;
    _d_max_0(&_S288, &_S289, _S287.differential_0);
    float _S290 = 0.5f * - (_S286 + _S288.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S291;
    (&_S291)->primal_0 = _S259;
    (&_S291)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S292;
    (&_S292)->primal_0 = _S250;
    (&_S292)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S291, &_S292, _S290);
    float _S293 = _S142 * _S195;
    DiffPair_float_0 _S294;
    (&_S294)->primal_0 = _S143;
    (&_S294)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S294, _S293);
    DiffPair_float_0 _S295;
    (&_S295)->primal_0 = cos_sim_loss_6;
    (&_S295)->differential_0 = 0.0f;
    DiffPair_float_0 _S296;
    (&_S296)->primal_0 = 9.999999960041972e-13f;
    (&_S296)->differential_0 = 0.0f;
    _d_max_0(&_S295, &_S296, _S294.differential_0);
    float _S297 = 0.5f * - (_S293 + _S295.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S298;
    (&_S298)->primal_0 = _S241;
    (&_S298)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S299;
    (&_S299)->primal_0 = _S250;
    (&_S299)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S298, &_S299, _S297);
    float3  _S300 = _S278.differential_0 + _S277.differential_0 + _S122;
    float3  _S301 = _S285.differential_0 + _S291.differential_0 + _S261.differential_0;
    float3  _S302 = _S284.differential_0 + _S298.differential_0 + _S243.differential_0;
    float3  _S303 = _S292.differential_0 + _S299.differential_0 + _S252.differential_0;
    if(_S138)
    {
        float3  _S304 = _S93.primal_0 * _S303;
        float3  _S305 = _S139 * _S303;
        float _S306 = _S304.x + _S304.y + _S304.z;
        DiffPair_float_0 _S307;
        (&_S307)->primal_0 = _S134;
        (&_S307)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S307, _S306);
        _S272 = _S307.differential_0;
        _S122 = _S305;
    }
    else
    {
        _S272 = 0.0f;
        _S122 = _S97;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S308;
    (&_S308)->primal_0 = _S93.primal_0;
    (&_S308)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S309;
    (&_S309)->primal_0 = _S93.primal_0;
    (&_S309)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S308, &_S309, _S272);
    float3  _S310 = _S309.differential_0 + _S308.differential_0 + _S122;
    if(_S130)
    {
        float3  _S311 = _S92.primal_0 * _S301;
        float3  _S312 = _S131 * _S301;
        float _S313 = _S311.x + _S311.y + _S311.z;
        DiffPair_float_0 _S314;
        (&_S314)->primal_0 = _S127;
        (&_S314)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S314, _S313);
        _S272 = _S314.differential_0;
        _S122 = _S312;
    }
    else
    {
        _S272 = 0.0f;
        _S122 = _S97;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S315;
    (&_S315)->primal_0 = _S92.primal_0;
    (&_S315)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S316;
    (&_S316)->primal_0 = _S92.primal_0;
    (&_S316)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S315, &_S316, _S272);
    float3  _S317 = _S316.differential_0 + _S315.differential_0 + _S122;
    if(_S123)
    {
        float3  _S318 = _S91.primal_0 * _S302;
        float3  _S319 = _S124 * _S302;
        float _S320 = _S318.x + _S318.y + _S318.z;
        DiffPair_float_0 _S321;
        (&_S321)->primal_0 = _S120;
        (&_S321)->differential_0 = 0.0f;
        s_bwd_prop_rsqrt_0(&_S321, _S320);
        _S272 = _S321.differential_0;
        _S122 = _S319;
    }
    else
    {
        _S272 = 0.0f;
        _S122 = _S97;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S322;
    (&_S322)->primal_0 = _S91.primal_0;
    (&_S322)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S323;
    (&_S323)->primal_0 = _S91.primal_0;
    (&_S323)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S322, &_S323, _S272);
    float3  _S324 = _S323.differential_0 + _S322.differential_0 + _S122;
    float _S325 = _S118 * _S194;
    float _S326 = _S118 * _S193;
    float _S327 = _S117 * _S192;
    float _S328 = _S116 * (_S117 * _S194 + _S326 + _S326 + _S191);
    DiffPair_float_0 _S329;
    (&_S329)->primal_0 = _S90.primal_0;
    (&_S329)->differential_0 = 0.0f;
    DiffPair_float_0 _S330;
    (&_S330)->primal_0 = 0.00009999999747379f;
    (&_S330)->differential_0 = 0.0f;
    _d_max_0(&_S329, &_S330, _S328);
    float _S331 = _S116 * (_S325 + _S327 + _S327 + _S190);
    DiffPair_float_0 _S332;
    (&_S332)->primal_0 = _S89.primal_0;
    (&_S332)->differential_0 = 0.0f;
    DiffPair_float_0 _S333;
    (&_S333)->primal_0 = 0.00009999999747379f;
    (&_S333)->differential_0 = 0.0f;
    _d_max_0(&_S332, &_S333, _S331);
    float _S334 = _S104 * _S189;
    DiffPair_float_0 _S335;
    (&_S335)->primal_0 = _S108;
    (&_S335)->differential_0 = 0.0f;
    DiffPair_float_0 _S336;
    (&_S336)->primal_0 = 0.0f;
    (&_S336)->differential_0 = 0.0f;
    DiffPair_float_0 _S337;
    (&_S337)->primal_0 = 1.0f;
    (&_S337)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S335, &_S336, &_S337, _S334);
    float _S338 = _S104 * _S188;
    float _S339 = _S115 * _S338;
    float _S340 = _S114 * (dV_1 * _S338);
    float _S341 = _S113 * _S338;
    float _S342 = _S112 * (dU_1 * _S338);
    float _S343 = _S111 * _S338;
    float _S344 = _S110 * (dY_1 * _S338);
    float _S345 = _S109 * _S338;
    DiffPair_float_0 _S346;
    (&_S346)->primal_0 = dY_1;
    (&_S346)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S346, _S345);
    float _S347 = 0.3333333432674408f * (_S335.differential_0 + _S107 * _S338);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S348;
    (&_S348)->primal_0 = _S106;
    (&_S348)->differential_0 = _S97;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S349;
    (&_S349)->primal_0 = _S106;
    (&_S349)->differential_0 = _S97;
    s_bwd_prop_dot_0(&_S348, &_S349, _S347);
    float _S350 = 0.3333333432674408f * (_S105 * _S338);
    float3  _S351 = make_float3 (_S350, _S350, _S350);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S352;
    (&_S352)->primal_0 = _S106;
    (&_S352)->differential_0 = _S97;
    s_bwd_prop_abs_1(&_S352, _S351);
    float3  _S353 = _S349.differential_0 + _S348.differential_0 + _S352.differential_0;
    float _S354 = _S339 + _S340;
    float s_diff_V_T_0 = - _S354;
    float _S355 = _S341 + _S342;
    float s_diff_U_T_0 = - _S355;
    float _S356 = _S343 + _S344 + _S346.differential_0;
    float s_diff_Y_T_0 = - _S356;
    float _S357 = - s_diff_V_T_0;
    float _S358 = _S332.differential_0 + _S268.differential_0;
    float3  _S359 = _S353 + make_float3 (0.61500000953674316f * s_diff_V_T_0 + -0.14712999761104584f * s_diff_U_T_0 + 0.29899999499320984f * s_diff_Y_T_0, 0.51498997211456299f * _S357 + 0.28885999321937561f * - s_diff_U_T_0 + 0.58700001239776611f * s_diff_Y_T_0, 0.10001000016927719f * _S357 + 0.43599998950958252f * s_diff_U_T_0 + 0.11400000005960464f * s_diff_Y_T_0);
    float3  _S360 = - _S353 + make_float3 (0.61500000953674316f * _S354 + -0.14712999761104584f * _S355 + 0.29899999499320984f * _S356, 0.51498997211456299f * s_diff_V_T_0 + 0.28885999321937561f * s_diff_U_T_0 + 0.58700001239776611f * _S356, 0.10001000016927719f * s_diff_V_T_0 + 0.43599998950958252f * _S355 + 0.11400000005960464f * _S356);
    dpmedian_normal_0->primal_0 = (*dpmedian_normal_0).primal_0;
    dpmedian_normal_0->differential_0 = _S300;
    dpmedian_depth_0->primal_0 = (*dpmedian_depth_0).primal_0;
    dpmedian_depth_0->differential_0 = _S266.differential_0;
    dpnormal_dist_0->primal_0 = (*dpnormal_dist_0).primal_0;
    dpnormal_dist_0->differential_0 = _S271;
    dpdepth_dist_0->primal_0 = (*dpdepth_dist_0).primal_0;
    dpdepth_dist_0->differential_0 = _S199;
    dprgb_dist_0->primal_0 = (*dprgb_dist_0).primal_0;
    dprgb_dist_0->differential_0 = _S270;
    dprender_Ts_0->primal_0 = (*dprender_Ts_0).primal_0;
    dprender_Ts_0->differential_0 = _S234;
    dpref_normal_0->primal_0 = (*dpref_normal_0).primal_0;
    dpref_normal_0->differential_0 = _S310;
    dpdepth_normal_0->primal_0 = (*dpdepth_normal_0).primal_0;
    dpdepth_normal_0->differential_0 = _S317;
    dprender_normal_0->primal_0 = (*dprender_normal_0).primal_0;
    dprender_normal_0->differential_0 = _S324;
    dpref_depth_0->primal_0 = (*dpref_depth_0).primal_0;
    dpref_depth_0->differential_0 = _S329.differential_0;
    dprender_depth_0->primal_0 = (*dprender_depth_0).primal_0;
    dprender_depth_0->differential_0 = _S358;
    dpref_rgb_0->primal_0 = (*dpref_rgb_0).primal_0;
    dpref_rgb_0->differential_0 = _S359;
    dprender_rgb_0->primal_0 = (*dprender_rgb_0).primal_0;
    dprender_rgb_0->differential_0 = _S360;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S361, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S362, DiffPair_float_0 * _S363, DiffPair_float_0 * _S364, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S365, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S366, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S367, DiffPair_float_0 * _S368, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S369, DiffPair_float_0 * _S370, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S371, DiffPair_float_0 * _S372, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S373, bool _S374, bool _S375, FixedArray<float, 19>  * _S376, FixedArray<float, 31>  * _S377)
{
    s_bwd_prop_per_pixel_losses_0(_S361, _S362, _S363, _S364, _S365, _S366, _S367, _S368, _S369, _S370, _S371, _S372, _S373, _S374, _S375, _S376, _S377);
    return;
}

inline __device__ void per_pixel_losses_bwd(float3  render_rgb_1, float3  ref_rgb_1, float render_depth_1, float ref_depth_1, float3  render_normal_1, float3  depth_normal_1, float3  ref_normal_1, float render_Ts_1, float3  rgb_dist_1, float depth_dist_1, float3  normal_dist_1, float median_depth_1, float3  median_normal_1, bool ref_alpha_2, bool has_mask_2, FixedArray<float, 19>  weights_2, FixedArray<float, 31>  v_losses_0, float3  * v_render_rgb_0, float3  * v_ref_rgb_0, float * v_render_depth_0, float * v_ref_depth_0, float3  * v_render_normal_0, float3  * v_depth_normal_0, float3  * v_ref_normal_0, float * v_render_Ts_0, float3  * v_rgb_dist_0, float * v_depth_dist_0, float3  * v_normal_dist_0, float * v_median_depth_0, float3  * v_median_normal_0)
{
    float3  _S378 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_rgb_0;
    (&dp_render_rgb_0)->primal_0 = render_rgb_1;
    (&dp_render_rgb_0)->differential_0 = _S378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_rgb_0;
    (&dp_ref_rgb_0)->primal_0 = ref_rgb_1;
    (&dp_ref_rgb_0)->differential_0 = _S378;
    DiffPair_float_0 dp_render_depth_0;
    (&dp_render_depth_0)->primal_0 = render_depth_1;
    (&dp_render_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 dp_ref_depth_0;
    (&dp_ref_depth_0)->primal_0 = ref_depth_1;
    (&dp_ref_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_render_normal_0;
    (&dp_render_normal_0)->primal_0 = render_normal_1;
    (&dp_render_normal_0)->differential_0 = _S378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depth_normal_0;
    (&dp_depth_normal_0)->primal_0 = depth_normal_1;
    (&dp_depth_normal_0)->differential_0 = _S378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ref_normal_0;
    (&dp_ref_normal_0)->primal_0 = ref_normal_1;
    (&dp_ref_normal_0)->differential_0 = _S378;
    DiffPair_float_0 dp_render_Ts_0;
    (&dp_render_Ts_0)->primal_0 = render_Ts_1;
    (&dp_render_Ts_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_dist_0;
    (&dp_rgb_dist_0)->primal_0 = rgb_dist_1;
    (&dp_rgb_dist_0)->differential_0 = _S378;
    DiffPair_float_0 dp_depth_dist_0;
    (&dp_depth_dist_0)->primal_0 = depth_dist_1;
    (&dp_depth_dist_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_normal_dist_0;
    (&dp_normal_dist_0)->primal_0 = normal_dist_1;
    (&dp_normal_dist_0)->differential_0 = _S378;
    DiffPair_float_0 dp_median_depth_0;
    (&dp_median_depth_0)->primal_0 = median_depth_1;
    (&dp_median_depth_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_median_normal_0;
    (&dp_median_normal_0)->primal_0 = median_normal_1;
    (&dp_median_normal_0)->differential_0 = _S378;
    FixedArray<float, 19>  _S379 = weights_2;
    FixedArray<float, 31>  _S380 = v_losses_0;
    s_bwd_per_pixel_losses_0(&dp_render_rgb_0, &dp_ref_rgb_0, &dp_render_depth_0, &dp_ref_depth_0, &dp_render_normal_0, &dp_depth_normal_0, &dp_ref_normal_0, &dp_render_Ts_0, &dp_rgb_dist_0, &dp_depth_dist_0, &dp_normal_dist_0, &dp_median_depth_0, &dp_median_normal_0, ref_alpha_2, has_mask_2, &_S379, &_S380);
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
    *v_median_depth_0 = dp_median_depth_0.differential_0;
    *v_median_normal_0 = dp_median_normal_0.differential_0;
    return;
}

inline __device__ void _d_log10_0(DiffPair_float_0 * dpx_9, float dOut_9)
{
    float _S381 = 1.0f / ((*dpx_9).primal_0 * 2.30258512496948242f) * dOut_9;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S381;
    return;
}

inline __device__ void per_pixel_losses_reduce(FixedArray<float, 31>  raw_losses_0, FixedArray<float, 19>  weights_3, FixedArray<float, 14>  * _S382)
{
    FixedArray<float, 14>  losses_1;
    float _S383 = (F32_max((raw_losses_0[int(21)]), (1.0f)));
    losses_1[int(0)] = raw_losses_0[int(0)] / _S383;
    losses_1[int(1)] = -10.0f * (F32_log10((raw_losses_0[int(1)] / _S383)));
    bool _S384;
    if((raw_losses_0[int(22)]) > 0.0f)
    {
        _S384 = (raw_losses_0[int(3)]) != 0.0f;
    }
    else
    {
        _S384 = false;
    }
    float _S385;
    if(_S384)
    {
        _S385 = weights_3[int(6)] * clamp_0(1.0f - (raw_losses_0[int(6)] - raw_losses_0[int(2)] * raw_losses_0[int(3)] / raw_losses_0[int(22)]) / (F32_sqrt(((F32_max((9.999999960041972e-13f), ((raw_losses_0[int(4)] - raw_losses_0[int(2)] * raw_losses_0[int(2)] / raw_losses_0[int(22)]) * (raw_losses_0[int(5)] - raw_losses_0[int(3)] * raw_losses_0[int(3)] / raw_losses_0[int(22)]) + 1.0f)))))), 0.0f, 2.0f);
    }
    else
    {
        _S385 = 0.0f;
    }
    losses_1[int(2)] = _S385;
    losses_1[int(3)] = (raw_losses_0[int(7)] / (F32_max((raw_losses_0[int(23)]), (1.0f))) + raw_losses_0[int(8)] / (F32_max((raw_losses_0[int(24)]), (1.0f)))) / float((I32_max((int((raw_losses_0[int(23)]) > 0.5f) + int((raw_losses_0[int(24)]) > 0.5f)), (int(1)))));
    losses_1[int(4)] = (raw_losses_0[int(9)] + raw_losses_0[int(10)]) / (F32_max((raw_losses_0[int(26)]), (1.0f)));
    losses_1[int(5)] = raw_losses_0[int(11)] / (F32_max((raw_losses_0[int(25)]), (1.0f)));
    float _S386 = (F32_max((raw_losses_0[int(20)]), (1.0f)));
    losses_1[int(6)] = raw_losses_0[int(12)] / _S386;
    losses_1[int(7)] = raw_losses_0[int(13)] / _S386;
    losses_1[int(8)] = raw_losses_0[int(14)] / _S386;
    losses_1[int(9)] = raw_losses_0[int(15)] / _S386;
    losses_1[int(10)] = raw_losses_0[int(16)] / (F32_max((raw_losses_0[int(27)]), (1.0f)));
    losses_1[int(11)] = raw_losses_0[int(17)] / (F32_max((raw_losses_0[int(28)]), (1.0f)));
    losses_1[int(12)] = raw_losses_0[int(18)] / (F32_max((raw_losses_0[int(29)]), (1.0f)));
    losses_1[int(13)] = raw_losses_0[int(19)] / (F32_max((raw_losses_0[int(30)]), (1.0f)));
    *_S382 = losses_1;
    return;
}

struct DiffPair_arrayx3Cfloatx2C31x3E_0
{
    FixedArray<float, 31>  primal_0;
    FixedArray<float, 31>  differential_0;
};

inline __device__ float s_primal_ctx_sqrt_0(float _S387)
{
    return (F32_sqrt((_S387)));
}

inline __device__ void s_bwd_prop_log10_0(DiffPair_float_0 * _S388, float _S389)
{
    _d_log10_0(_S388, _S389);
    return;
}

inline __device__ void s_bwd_prop_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C31x3E_0 * dpraw_losses_0, FixedArray<float, 19>  * weights_4, FixedArray<float, 14>  * _s_dOut_1)
{
    FixedArray<float, 31>  _S390 = dpraw_losses_0->primal_0;
    float _S391 = (F32_max((dpraw_losses_0->primal_0[int(21)]), (1.0f)));
    float _S392 = _S391 * _S391;
    float _S393 = dpraw_losses_0->primal_0[int(1)] / _S391;
    bool _S394 = (dpraw_losses_0->primal_0[int(22)]) > 0.0f;
    bool _S395;
    if(_S394)
    {
        _S395 = (_S390[int(3)]) != 0.0f;
    }
    else
    {
        _S395 = false;
    }
    float _S396;
    float _S397;
    float _S398;
    float _S399;
    float _S400;
    float _S401;
    float _S402;
    float _S403;
    float _S404;
    float _S405;
    float _S406;
    float _S407;
    float _S408;
    float _S409;
    float _S410;
    if(_S395)
    {
        float _S411 = _S390[int(2)] * _S390[int(3)];
        float _S412 = _S390[int(22)] * _S390[int(22)];
        float _S413 = _S390[int(6)] - _S411 / _S390[int(22)];
        float _S414 = _S390[int(2)] * _S390[int(2)];
        float _S415 = _S390[int(4)] - _S414 / _S390[int(22)];
        float _S416 = _S390[int(3)] * _S390[int(3)];
        float _S417 = _S390[int(5)] - _S416 / _S390[int(22)];
        float _S418 = _S415 * _S417 + 1.0f;
        float _S419 = (F32_max((9.999999960041972e-13f), (_S418)));
        float _S420 = s_primal_ctx_sqrt_0(_S419);
        float _S421 = _S420 * _S420;
        float _S422 = 1.0f - _S413 / _S420;
        _S396 = (*weights_4)[int(6)];
        _S397 = _S422;
        _S398 = _S421;
        _S399 = _S413;
        _S400 = _S420;
        _S401 = _S419;
        _S402 = _S418;
        _S403 = _S415;
        _S404 = _S417;
        _S405 = _S412;
        _S406 = _S416;
        _S407 = _S390[int(3)];
        _S408 = _S414;
        _S409 = _S390[int(2)];
        _S410 = _S411;
    }
    else
    {
        _S396 = 0.0f;
        _S397 = 0.0f;
        _S398 = 0.0f;
        _S399 = 0.0f;
        _S400 = 0.0f;
        _S401 = 0.0f;
        _S402 = 0.0f;
        _S403 = 0.0f;
        _S404 = 0.0f;
        _S405 = 0.0f;
        _S406 = 0.0f;
        _S407 = 0.0f;
        _S408 = 0.0f;
        _S409 = 0.0f;
        _S410 = 0.0f;
    }
    float _S423 = (F32_max((_S390[int(23)]), (1.0f)));
    float _S424 = _S423 * _S423;
    float _S425 = (F32_max((_S390[int(24)]), (1.0f)));
    float _S426 = _S425 * _S425;
    float _S427 = float((I32_max((int((_S390[int(23)]) > 0.5f) + int((_S390[int(24)]) > 0.5f)), (int(1)))));
    float _S428 = _S390[int(9)] + _S390[int(10)];
    float _S429 = (F32_max((_S390[int(26)]), (1.0f)));
    float _S430 = _S429 * _S429;
    float _S431 = (F32_max((_S390[int(25)]), (1.0f)));
    float _S432 = _S431 * _S431;
    float _S433 = (F32_max((_S390[int(20)]), (1.0f)));
    float _S434 = _S433 * _S433;
    float _S435 = (F32_max((_S390[int(27)]), (1.0f)));
    float _S436 = _S435 * _S435;
    float _S437 = (F32_max((_S390[int(28)]), (1.0f)));
    float _S438 = _S437 * _S437;
    float _S439 = (F32_max((_S390[int(29)]), (1.0f)));
    float _S440 = _S439 * _S439;
    float _S441 = (F32_max((_S390[int(30)]), (1.0f)));
    float _S442 = _S441 * _S441;
    float _S443 = (*_s_dOut_1)[int(0)];
    float _S444 = (*_s_dOut_1)[int(1)];
    float _S445 = (*_s_dOut_1)[int(2)];
    float _S446 = (*_s_dOut_1)[int(13)] / _S442;
    float _S447 = _S390[int(19)] * - _S446;
    float _S448 = _S441 * _S446;
    DiffPair_float_0 _S449;
    (&_S449)->primal_0 = _S390[int(30)];
    (&_S449)->differential_0 = 0.0f;
    DiffPair_float_0 _S450;
    (&_S450)->primal_0 = 1.0f;
    (&_S450)->differential_0 = 0.0f;
    _d_max_0(&_S449, &_S450, _S447);
    float _S451 = (*_s_dOut_1)[int(12)] / _S440;
    float _S452 = _S390[int(18)] * - _S451;
    float _S453 = _S439 * _S451;
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S390[int(29)];
    (&_S454)->differential_0 = 0.0f;
    DiffPair_float_0 _S455;
    (&_S455)->primal_0 = 1.0f;
    (&_S455)->differential_0 = 0.0f;
    _d_max_0(&_S454, &_S455, _S452);
    float _S456 = (*_s_dOut_1)[int(11)] / _S438;
    float _S457 = _S390[int(17)] * - _S456;
    float _S458 = _S437 * _S456;
    DiffPair_float_0 _S459;
    (&_S459)->primal_0 = _S390[int(28)];
    (&_S459)->differential_0 = 0.0f;
    DiffPair_float_0 _S460;
    (&_S460)->primal_0 = 1.0f;
    (&_S460)->differential_0 = 0.0f;
    _d_max_0(&_S459, &_S460, _S457);
    float _S461 = (*_s_dOut_1)[int(10)] / _S436;
    float _S462 = _S390[int(16)] * - _S461;
    float _S463 = _S435 * _S461;
    DiffPair_float_0 _S464;
    (&_S464)->primal_0 = _S390[int(27)];
    (&_S464)->differential_0 = 0.0f;
    DiffPair_float_0 _S465;
    (&_S465)->primal_0 = 1.0f;
    (&_S465)->differential_0 = 0.0f;
    _d_max_0(&_S464, &_S465, _S462);
    float _S466 = (*_s_dOut_1)[int(9)] / _S434;
    float _S467 = _S433 * _S466;
    float _S468 = (*_s_dOut_1)[int(8)] / _S434;
    float _S469 = _S433 * _S468;
    float _S470 = (*_s_dOut_1)[int(7)] / _S434;
    float _S471 = _S433 * _S470;
    float _S472 = (*_s_dOut_1)[int(6)] / _S434;
    float _S473 = _S433 * _S472;
    float _S474 = _S390[int(15)] * - _S466 + _S390[int(14)] * - _S468 + _S390[int(13)] * - _S470 + _S390[int(12)] * - _S472;
    DiffPair_float_0 _S475;
    (&_S475)->primal_0 = _S390[int(20)];
    (&_S475)->differential_0 = 0.0f;
    DiffPair_float_0 _S476;
    (&_S476)->primal_0 = 1.0f;
    (&_S476)->differential_0 = 0.0f;
    _d_max_0(&_S475, &_S476, _S474);
    float _S477 = (*_s_dOut_1)[int(5)] / _S432;
    float _S478 = _S390[int(11)] * - _S477;
    float _S479 = _S431 * _S477;
    DiffPair_float_0 _S480;
    (&_S480)->primal_0 = _S390[int(25)];
    (&_S480)->differential_0 = 0.0f;
    DiffPair_float_0 _S481;
    (&_S481)->primal_0 = 1.0f;
    (&_S481)->differential_0 = 0.0f;
    _d_max_0(&_S480, &_S481, _S478);
    float _S482 = (*_s_dOut_1)[int(4)] / _S430;
    float _S483 = _S428 * - _S482;
    float _S484 = _S429 * _S482;
    DiffPair_float_0 _S485;
    (&_S485)->primal_0 = _S390[int(26)];
    (&_S485)->differential_0 = 0.0f;
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = 1.0f;
    (&_S486)->differential_0 = 0.0f;
    _d_max_0(&_S485, &_S486, _S483);
    float _S487 = (*_s_dOut_1)[int(3)] / _S427;
    float _S488 = _S487 / _S426;
    float _S489 = _S390[int(8)] * - _S488;
    float _S490 = _S425 * _S488;
    DiffPair_float_0 _S491;
    (&_S491)->primal_0 = _S390[int(24)];
    (&_S491)->differential_0 = 0.0f;
    DiffPair_float_0 _S492;
    (&_S492)->primal_0 = 1.0f;
    (&_S492)->differential_0 = 0.0f;
    _d_max_0(&_S491, &_S492, _S489);
    float _S493 = _S487 / _S424;
    float _S494 = _S390[int(7)] * - _S493;
    float _S495 = _S423 * _S493;
    DiffPair_float_0 _S496;
    (&_S496)->primal_0 = _S390[int(23)];
    (&_S496)->differential_0 = 0.0f;
    DiffPair_float_0 _S497;
    (&_S497)->primal_0 = 1.0f;
    (&_S497)->differential_0 = 0.0f;
    _d_max_0(&_S496, &_S497, _S494);
    FixedArray<float, 31>  _S498;
    _S498[int(0)] = 0.0f;
    _S498[int(1)] = 0.0f;
    _S498[int(2)] = 0.0f;
    _S498[int(3)] = 0.0f;
    _S498[int(4)] = 0.0f;
    _S498[int(5)] = 0.0f;
    _S498[int(6)] = 0.0f;
    _S498[int(7)] = 0.0f;
    _S498[int(8)] = 0.0f;
    _S498[int(9)] = 0.0f;
    _S498[int(10)] = 0.0f;
    _S498[int(11)] = 0.0f;
    _S498[int(12)] = 0.0f;
    _S498[int(13)] = 0.0f;
    _S498[int(14)] = 0.0f;
    _S498[int(15)] = 0.0f;
    _S498[int(16)] = 0.0f;
    _S498[int(17)] = 0.0f;
    _S498[int(18)] = 0.0f;
    _S498[int(19)] = 0.0f;
    _S498[int(20)] = 0.0f;
    _S498[int(21)] = 0.0f;
    _S498[int(22)] = 0.0f;
    _S498[int(23)] = 0.0f;
    _S498[int(24)] = 0.0f;
    _S498[int(25)] = 0.0f;
    _S498[int(26)] = 0.0f;
    _S498[int(27)] = 0.0f;
    _S498[int(28)] = 0.0f;
    _S498[int(29)] = 0.0f;
    _S498[int(30)] = 0.0f;
    _S498[int(20)] = _S475.differential_0;
    _S498[int(7)] = _S495;
    _S498[int(23)] = _S496.differential_0;
    _S498[int(8)] = _S490;
    _S498[int(24)] = _S491.differential_0;
    _S498[int(9)] = _S484;
    _S498[int(10)] = _S484;
    _S498[int(26)] = _S485.differential_0;
    _S498[int(11)] = _S479;
    _S498[int(25)] = _S480.differential_0;
    _S498[int(12)] = _S473;
    _S498[int(30)] = _S449.differential_0;
    _S498[int(13)] = _S471;
    _S498[int(14)] = _S469;
    _S498[int(15)] = _S467;
    _S498[int(16)] = _S463;
    _S498[int(27)] = _S464.differential_0;
    _S498[int(17)] = _S458;
    _S498[int(28)] = _S459.differential_0;
    _S498[int(18)] = _S453;
    _S498[int(29)] = _S454.differential_0;
    _S498[int(19)] = _S448;
    float _S499 = _S498[int(0)];
    float _S500 = _S498[int(1)];
    float _S501 = _S498[int(2)];
    float _S502 = _S498[int(3)];
    float _S503 = _S498[int(4)];
    float _S504 = _S498[int(5)];
    float _S505 = _S498[int(6)];
    float _S506 = _S498[int(7)];
    float _S507 = _S498[int(8)];
    float _S508 = _S498[int(9)];
    float _S509 = _S498[int(10)];
    float _S510 = _S498[int(11)];
    float _S511 = _S498[int(12)];
    float _S512 = _S498[int(13)];
    float _S513 = _S498[int(14)];
    float _S514 = _S498[int(15)];
    float _S515 = _S498[int(16)];
    float _S516 = _S498[int(17)];
    float _S517 = _S498[int(18)];
    float _S518 = _S498[int(19)];
    float _S519 = _S498[int(20)];
    float _S520 = _S498[int(21)];
    float _S521 = _S498[int(22)];
    float _S522 = _S498[int(23)];
    float _S523 = _S498[int(24)];
    float _S524 = _S498[int(25)];
    float _S525 = _S498[int(26)];
    float _S526 = _S498[int(27)];
    float _S527 = _S498[int(28)];
    float _S528 = _S498[int(29)];
    float _S529 = _S498[int(30)];
    FixedArray<float, 31>  _S530;
    if(_S395)
    {
        float _S531 = _S396 * _S445;
        DiffPair_float_0 _S532;
        (&_S532)->primal_0 = _S397;
        (&_S532)->differential_0 = 0.0f;
        DiffPair_float_0 _S533;
        (&_S533)->primal_0 = 0.0f;
        (&_S533)->differential_0 = 0.0f;
        DiffPair_float_0 _S534;
        (&_S534)->primal_0 = 2.0f;
        (&_S534)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S532, &_S533, &_S534, _S531);
        float _S535 = - _S532.differential_0 / _S398;
        float _S536 = _S399 * - _S535;
        float _S537 = _S400 * _S535;
        DiffPair_float_0 _S538;
        (&_S538)->primal_0 = _S401;
        (&_S538)->differential_0 = 0.0f;
        s_bwd_prop_sqrt_0(&_S538, _S536);
        DiffPair_float_0 _S539;
        (&_S539)->primal_0 = 9.999999960041972e-13f;
        (&_S539)->differential_0 = 0.0f;
        DiffPair_float_0 _S540;
        (&_S540)->primal_0 = _S402;
        (&_S540)->differential_0 = 0.0f;
        _d_max_0(&_S539, &_S540, _S538.differential_0);
        float _S541 = _S403 * _S540.differential_0;
        float _S542 = _S404 * _S540.differential_0;
        float _S543 = - _S541 / _S405;
        float _S544 = _S407 * (_S390[int(22)] * _S543);
        float _S545 = - _S542 / _S405;
        float _S546 = _S409 * (_S390[int(22)] * _S545);
        float _S547 = - _S537 / _S405;
        float _S548 = _S390[int(22)] * _S547;
        float _S549 = _S544 + _S544 + _S409 * _S548;
        float _S550 = _S546 + _S546 + _S407 * _S548;
        float _S551 = _S406 * - _S543 + _S408 * - _S545 + _S410 * - _S547;
        FixedArray<float, 31>  _S552;
        _S552[int(0)] = 0.0f;
        _S552[int(1)] = 0.0f;
        _S552[int(2)] = 0.0f;
        _S552[int(3)] = 0.0f;
        _S552[int(4)] = 0.0f;
        _S552[int(5)] = 0.0f;
        _S552[int(6)] = 0.0f;
        _S552[int(7)] = 0.0f;
        _S552[int(8)] = 0.0f;
        _S552[int(9)] = 0.0f;
        _S552[int(10)] = 0.0f;
        _S552[int(11)] = 0.0f;
        _S552[int(12)] = 0.0f;
        _S552[int(13)] = 0.0f;
        _S552[int(14)] = 0.0f;
        _S552[int(15)] = 0.0f;
        _S552[int(16)] = 0.0f;
        _S552[int(17)] = 0.0f;
        _S552[int(18)] = 0.0f;
        _S552[int(19)] = 0.0f;
        _S552[int(20)] = 0.0f;
        _S552[int(21)] = 0.0f;
        _S552[int(22)] = 0.0f;
        _S552[int(23)] = 0.0f;
        _S552[int(24)] = 0.0f;
        _S552[int(25)] = 0.0f;
        _S552[int(26)] = 0.0f;
        _S552[int(27)] = 0.0f;
        _S552[int(28)] = 0.0f;
        _S552[int(29)] = 0.0f;
        _S552[int(30)] = 0.0f;
        _S552[int(5)] = _S541;
        _S552[int(4)] = _S542;
        _S552[int(3)] = _S549;
        _S552[int(2)] = _S550;
        _S552[int(6)] = _S537;
        float _S553 = _S500 + _S552[int(1)];
        float _S554 = _S501 + _S552[int(2)];
        float _S555 = _S502 + _S552[int(3)];
        float _S556 = _S503 + _S552[int(4)];
        float _S557 = _S504 + _S552[int(5)];
        float _S558 = _S505 + _S552[int(6)];
        float _S559 = _S506 + _S552[int(7)];
        float _S560 = _S507 + _S552[int(8)];
        float _S561 = _S508 + _S552[int(9)];
        float _S562 = _S509 + _S552[int(10)];
        float _S563 = _S510 + _S552[int(11)];
        float _S564 = _S511 + _S552[int(12)];
        float _S565 = _S512 + _S552[int(13)];
        float _S566 = _S513 + _S552[int(14)];
        float _S567 = _S514 + _S552[int(15)];
        float _S568 = _S515 + _S552[int(16)];
        float _S569 = _S516 + _S552[int(17)];
        float _S570 = _S517 + _S552[int(18)];
        float _S571 = _S518 + _S552[int(19)];
        float _S572 = _S519 + _S552[int(20)];
        float _S573 = _S520 + _S552[int(21)];
        float _S574 = _S521 + _S552[int(22)];
        float _S575 = _S522 + _S552[int(23)];
        float _S576 = _S523 + _S552[int(24)];
        float _S577 = _S524 + _S552[int(25)];
        float _S578 = _S525 + _S552[int(26)];
        float _S579 = _S526 + _S552[int(27)];
        float _S580 = _S527 + _S552[int(28)];
        float _S581 = _S528 + _S552[int(29)];
        float _S582 = _S529 + _S552[int(30)];
        _S530[int(0)] = _S499 + _S552[int(0)];
        _S530[int(1)] = _S553;
        _S530[int(2)] = _S554;
        _S530[int(3)] = _S555;
        _S530[int(4)] = _S556;
        _S530[int(5)] = _S557;
        _S530[int(6)] = _S558;
        _S530[int(7)] = _S559;
        _S530[int(8)] = _S560;
        _S530[int(9)] = _S561;
        _S530[int(10)] = _S562;
        _S530[int(11)] = _S563;
        _S530[int(12)] = _S564;
        _S530[int(13)] = _S565;
        _S530[int(14)] = _S566;
        _S530[int(15)] = _S567;
        _S530[int(16)] = _S568;
        _S530[int(17)] = _S569;
        _S530[int(18)] = _S570;
        _S530[int(19)] = _S571;
        _S530[int(20)] = _S572;
        _S530[int(21)] = _S573;
        _S530[int(22)] = _S574;
        _S530[int(23)] = _S575;
        _S530[int(24)] = _S576;
        _S530[int(25)] = _S577;
        _S530[int(26)] = _S578;
        _S530[int(27)] = _S579;
        _S530[int(28)] = _S580;
        _S530[int(29)] = _S581;
        _S530[int(30)] = _S582;
        _S396 = _S551;
    }
    else
    {
        _S530[int(0)] = _S499;
        _S530[int(1)] = _S500;
        _S530[int(2)] = _S501;
        _S530[int(3)] = _S502;
        _S530[int(4)] = _S503;
        _S530[int(5)] = _S504;
        _S530[int(6)] = _S505;
        _S530[int(7)] = _S506;
        _S530[int(8)] = _S507;
        _S530[int(9)] = _S508;
        _S530[int(10)] = _S509;
        _S530[int(11)] = _S510;
        _S530[int(12)] = _S511;
        _S530[int(13)] = _S512;
        _S530[int(14)] = _S513;
        _S530[int(15)] = _S514;
        _S530[int(16)] = _S515;
        _S530[int(17)] = _S516;
        _S530[int(18)] = _S517;
        _S530[int(19)] = _S518;
        _S530[int(20)] = _S519;
        _S530[int(21)] = _S520;
        _S530[int(22)] = _S521;
        _S530[int(23)] = _S522;
        _S530[int(24)] = _S523;
        _S530[int(25)] = _S524;
        _S530[int(26)] = _S525;
        _S530[int(27)] = _S526;
        _S530[int(28)] = _S527;
        _S530[int(29)] = _S528;
        _S530[int(30)] = _S529;
        _S396 = 0.0f;
    }
    if(_S394)
    {
        FixedArray<float, 31>  _S583;
        _S583[int(0)] = 0.0f;
        _S583[int(1)] = 0.0f;
        _S583[int(2)] = 0.0f;
        _S583[int(3)] = 0.0f;
        _S583[int(4)] = 0.0f;
        _S583[int(5)] = 0.0f;
        _S583[int(6)] = 0.0f;
        _S583[int(7)] = 0.0f;
        _S583[int(8)] = 0.0f;
        _S583[int(9)] = 0.0f;
        _S583[int(10)] = 0.0f;
        _S583[int(11)] = 0.0f;
        _S583[int(12)] = 0.0f;
        _S583[int(13)] = 0.0f;
        _S583[int(14)] = 0.0f;
        _S583[int(15)] = 0.0f;
        _S583[int(16)] = 0.0f;
        _S583[int(17)] = 0.0f;
        _S583[int(18)] = 0.0f;
        _S583[int(19)] = 0.0f;
        _S583[int(20)] = 0.0f;
        _S583[int(21)] = 0.0f;
        _S583[int(22)] = 0.0f;
        _S583[int(23)] = 0.0f;
        _S583[int(24)] = 0.0f;
        _S583[int(25)] = 0.0f;
        _S583[int(26)] = 0.0f;
        _S583[int(27)] = 0.0f;
        _S583[int(28)] = 0.0f;
        _S583[int(29)] = 0.0f;
        _S583[int(30)] = 0.0f;
        _S583[int(3)] = 0.0f;
        float _S584 = _S530[int(1)] + _S583[int(1)];
        float _S585 = _S530[int(2)] + _S583[int(2)];
        float _S586 = _S530[int(3)] + _S583[int(3)];
        float _S587 = _S530[int(4)] + _S583[int(4)];
        float _S588 = _S530[int(5)] + _S583[int(5)];
        float _S589 = _S530[int(6)] + _S583[int(6)];
        float _S590 = _S530[int(7)] + _S583[int(7)];
        float _S591 = _S530[int(8)] + _S583[int(8)];
        float _S592 = _S530[int(9)] + _S583[int(9)];
        float _S593 = _S530[int(10)] + _S583[int(10)];
        float _S594 = _S530[int(11)] + _S583[int(11)];
        float _S595 = _S530[int(12)] + _S583[int(12)];
        float _S596 = _S530[int(13)] + _S583[int(13)];
        float _S597 = _S530[int(14)] + _S583[int(14)];
        float _S598 = _S530[int(15)] + _S583[int(15)];
        float _S599 = _S530[int(16)] + _S583[int(16)];
        float _S600 = _S530[int(17)] + _S583[int(17)];
        float _S601 = _S530[int(18)] + _S583[int(18)];
        float _S602 = _S530[int(19)] + _S583[int(19)];
        float _S603 = _S530[int(20)] + _S583[int(20)];
        float _S604 = _S530[int(21)] + _S583[int(21)];
        float _S605 = _S530[int(22)] + _S583[int(22)];
        float _S606 = _S530[int(23)] + _S583[int(23)];
        float _S607 = _S530[int(24)] + _S583[int(24)];
        float _S608 = _S530[int(25)] + _S583[int(25)];
        float _S609 = _S530[int(26)] + _S583[int(26)];
        float _S610 = _S530[int(27)] + _S583[int(27)];
        float _S611 = _S530[int(28)] + _S583[int(28)];
        float _S612 = _S530[int(29)] + _S583[int(29)];
        float _S613 = _S530[int(30)] + _S583[int(30)];
        _S530[int(0)] = _S530[int(0)] + _S583[int(0)];
        _S530[int(1)] = _S584;
        _S530[int(2)] = _S585;
        _S530[int(3)] = _S586;
        _S530[int(4)] = _S587;
        _S530[int(5)] = _S588;
        _S530[int(6)] = _S589;
        _S530[int(7)] = _S590;
        _S530[int(8)] = _S591;
        _S530[int(9)] = _S592;
        _S530[int(10)] = _S593;
        _S530[int(11)] = _S594;
        _S530[int(12)] = _S595;
        _S530[int(13)] = _S596;
        _S530[int(14)] = _S597;
        _S530[int(15)] = _S598;
        _S530[int(16)] = _S599;
        _S530[int(17)] = _S600;
        _S530[int(18)] = _S601;
        _S530[int(19)] = _S602;
        _S530[int(20)] = _S603;
        _S530[int(21)] = _S604;
        _S530[int(22)] = _S605;
        _S530[int(23)] = _S606;
        _S530[int(24)] = _S607;
        _S530[int(25)] = _S608;
        _S530[int(26)] = _S609;
        _S530[int(27)] = _S610;
        _S530[int(28)] = _S611;
        _S530[int(29)] = _S612;
        _S530[int(30)] = _S613;
    }
    float _S614 = -10.0f * _S444;
    DiffPair_float_0 _S615;
    (&_S615)->primal_0 = _S393;
    (&_S615)->differential_0 = 0.0f;
    s_bwd_prop_log10_0(&_S615, _S614);
    float _S616 = _S615.differential_0 / _S392;
    float _S617 = _S391 * _S616;
    float _S618 = _S443 / _S392;
    float _S619 = _S391 * _S618;
    float _S620 = _S390[int(1)] * - _S616 + _S390[int(0)] * - _S618;
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S390[int(21)];
    (&_S621)->differential_0 = 0.0f;
    DiffPair_float_0 _S622;
    (&_S622)->primal_0 = 1.0f;
    (&_S622)->differential_0 = 0.0f;
    _d_max_0(&_S621, &_S622, _S620);
    FixedArray<float, 31>  _S623;
    _S623[int(0)] = 0.0f;
    _S623[int(1)] = 0.0f;
    _S623[int(2)] = 0.0f;
    _S623[int(3)] = 0.0f;
    _S623[int(4)] = 0.0f;
    _S623[int(5)] = 0.0f;
    _S623[int(6)] = 0.0f;
    _S623[int(7)] = 0.0f;
    _S623[int(8)] = 0.0f;
    _S623[int(9)] = 0.0f;
    _S623[int(10)] = 0.0f;
    _S623[int(11)] = 0.0f;
    _S623[int(12)] = 0.0f;
    _S623[int(13)] = 0.0f;
    _S623[int(14)] = 0.0f;
    _S623[int(15)] = 0.0f;
    _S623[int(16)] = 0.0f;
    _S623[int(17)] = 0.0f;
    _S623[int(18)] = 0.0f;
    _S623[int(19)] = 0.0f;
    _S623[int(20)] = 0.0f;
    _S623[int(21)] = 0.0f;
    _S623[int(22)] = 0.0f;
    _S623[int(23)] = 0.0f;
    _S623[int(24)] = 0.0f;
    _S623[int(25)] = 0.0f;
    _S623[int(26)] = 0.0f;
    _S623[int(27)] = 0.0f;
    _S623[int(28)] = 0.0f;
    _S623[int(29)] = 0.0f;
    _S623[int(30)] = 0.0f;
    _S623[int(22)] = _S396;
    _S623[int(1)] = _S617;
    _S623[int(21)] = _S621.differential_0;
    _S623[int(0)] = _S619;
    FixedArray<float, 31>  _S624 = {
        _S530[int(0)] + _S623[int(0)], _S530[int(1)] + _S623[int(1)], _S530[int(2)] + _S623[int(2)], _S530[int(3)] + _S623[int(3)], _S530[int(4)] + _S623[int(4)], _S530[int(5)] + _S623[int(5)], _S530[int(6)] + _S623[int(6)], _S530[int(7)] + _S623[int(7)], _S530[int(8)] + _S623[int(8)], _S530[int(9)] + _S623[int(9)], _S530[int(10)] + _S623[int(10)], _S530[int(11)] + _S623[int(11)], _S530[int(12)] + _S623[int(12)], _S530[int(13)] + _S623[int(13)], _S530[int(14)] + _S623[int(14)], _S530[int(15)] + _S623[int(15)], _S530[int(16)] + _S623[int(16)], _S530[int(17)] + _S623[int(17)], _S530[int(18)] + _S623[int(18)], _S530[int(19)] + _S623[int(19)], _S530[int(20)] + _S623[int(20)], _S530[int(21)] + _S623[int(21)], _S530[int(22)] + _S623[int(22)], _S530[int(23)] + _S623[int(23)], _S530[int(24)] + _S623[int(24)], _S530[int(25)] + _S623[int(25)], _S530[int(26)] + _S623[int(26)], _S530[int(27)] + _S623[int(27)], _S530[int(28)] + _S623[int(28)], _S530[int(29)] + _S623[int(29)], _S530[int(30)] + _S623[int(30)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S624;
    return;
}

inline __device__ void s_bwd_per_pixel_losses_reduce_0(DiffPair_arrayx3Cfloatx2C31x3E_0 * _S625, FixedArray<float, 19>  * _S626, FixedArray<float, 14>  * _S627)
{
    s_bwd_prop_per_pixel_losses_reduce_0(_S625, _S626, _S627);
    return;
}

inline __device__ void per_pixel_losses_reduce_bwd(FixedArray<float, 31>  raw_losses_1, FixedArray<float, 19>  weights_5, FixedArray<float, 14>  v_losses_1, FixedArray<float, 31>  * _S628)
{
    FixedArray<float, 31>  _S629 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C31x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_1;
    (&dp_raw_losses_0)->differential_0 = _S629;
    FixedArray<float, 19>  _S630 = weights_5;
    FixedArray<float, 14>  _S631 = v_losses_1;
    s_bwd_per_pixel_losses_reduce_0(&dp_raw_losses_0, &_S630, &_S631);
    *_S628 = (&dp_raw_losses_0)->differential_0;
    return;
}

