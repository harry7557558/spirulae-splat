#pragma once

#include "slang.cuh"

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_0)
{
    DiffPair_float_0 _S1 = *dpx_0;
    float _S2;
    if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
        {
            _S2 = 0.0f;
        }
        else
        {
            _S2 = 0.5f * dOut_0;
        }
    }
    dpx_0->primal_0 = _S1.primal_0;
    dpx_0->differential_0 = _S2;
    DiffPair_float_0 _S3 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S1.primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_0).primal_0))
        {
            _S2 = 0.0f;
        }
        else
        {
            _S2 = 0.5f * dOut_0;
        }
    }
    dpy_0->primal_0 = _S3.primal_0;
    dpy_0->differential_0 = _S2;
    return;
}

inline __device__ float rendered_depth_to_expected_depth(float depth_0, float transmittance_0)
{
    return depth_0 / (F32_max((1.0f - transmittance_0), (1.00000001335143196e-10f)));
}

inline __device__ void s_bwd_prop_rendered_depth_to_expected_depth_0(DiffPair_float_0 * dpdepth_0, DiffPair_float_0 * dptransmittance_0, float _s_dOut_0)
{
    float _S4 = 1.0f - (*dptransmittance_0).primal_0;
    float _S5 = (F32_max((_S4), (1.00000001335143196e-10f)));
    float _S6 = _s_dOut_0 / (_S5 * _S5);
    float _S7 = (*dpdepth_0).primal_0 * - _S6;
    float _S8 = _S5 * _S6;
    DiffPair_float_0 _S9;
    (&_S9)->primal_0 = _S4;
    (&_S9)->differential_0 = 0.0f;
    DiffPair_float_0 _S10;
    (&_S10)->primal_0 = 1.00000001335143196e-10f;
    (&_S10)->differential_0 = 0.0f;
    _d_max_0(&_S9, &_S10, _S7);
    float _S11 = - _S9.differential_0;
    dptransmittance_0->primal_0 = (*dptransmittance_0).primal_0;
    dptransmittance_0->differential_0 = _S11;
    dpdepth_0->primal_0 = (*dpdepth_0).primal_0;
    dpdepth_0->differential_0 = _S8;
    return;
}

inline __device__ void s_bwd_rendered_depth_to_expected_depth_0(DiffPair_float_0 * _S12, DiffPair_float_0 * _S13, float _S14)
{
    s_bwd_prop_rendered_depth_to_expected_depth_0(_S12, _S13, _S14);
    return;
}

inline __device__ void rendered_depth_to_expected_depth_bwd(float depth_1, float transmittance_1, float v_out_depth_0, float * v_depth_0, float * v_transmittance_0)
{
    DiffPair_float_0 p_depth_0;
    (&p_depth_0)->primal_0 = depth_1;
    (&p_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 p_transmittance_0;
    (&p_transmittance_0)->primal_0 = transmittance_1;
    (&p_transmittance_0)->differential_0 = 0.0f;
    s_bwd_rendered_depth_to_expected_depth_0(&p_depth_0, &p_transmittance_0, v_out_depth_0);
    *v_depth_0 = p_depth_0.differential_0;
    *v_transmittance_0 = p_transmittance_0.differential_0;
    return;
}

inline __device__ float3  min_0(float3  x_0, float3  y_0)
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
        *_slang_vector_get_element_ptr(&result_0, i_0) = (F32_min((_slang_vector_get_element(x_0, i_0)), (_slang_vector_get_element(y_0, i_0))));
        i_0 = i_0 + int(1);
    }
    return result_0;
}

inline __device__ float3  max_0(float3  x_1, float3  y_1)
{
    float3  result_1;
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
        *_slang_vector_get_element_ptr(&result_1, i_1) = (F32_max((_slang_vector_get_element(x_1, i_1)), (_slang_vector_get_element(y_1, i_1))));
        i_1 = i_1 + int(1);
    }
    return result_1;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_1)
{
    DiffPair_float_0 _S15 = *dpx_1;
    bool _S16;
    if(((*dpx_1).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S16 = ((*dpx_1).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S16 = false;
    }
    float _S17;
    if(_S16)
    {
        _S17 = dOut_1;
    }
    else
    {
        _S17 = 0.0f;
    }
    dpx_1->primal_0 = _S15.primal_0;
    dpx_1->differential_0 = _S17;
    DiffPair_float_0 _S18 = *dpMin_0;
    if((_S15.primal_0) < ((*dpMin_0).primal_0))
    {
        _S17 = dOut_1;
    }
    else
    {
        _S17 = 0.0f;
    }
    dpMin_0->primal_0 = _S18.primal_0;
    dpMin_0->differential_0 = _S17;
    DiffPair_float_0 _S19 = *dpMax_0;
    if(((*dpx_1).primal_0) > ((*dpMax_0).primal_0))
    {
        _S17 = dOut_1;
    }
    else
    {
        _S17 = 0.0f;
    }
    dpMax_0->primal_0 = _S19.primal_0;
    dpMax_0->differential_0 = _S17;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_2)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_2).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_1).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_2.x);
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_2).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_1).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_2.y);
    *&((&left_d_result_0)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_0)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_2).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_1).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_2.z);
    *&((&left_d_result_0)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_0)->z) = right_dp_2.differential_0;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = left_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  clamp_0(float3  x_2, float3  minBound_0, float3  maxBound_0)
{
    return min_0(max_0(x_2, minBound_0), maxBound_0);
}

inline __device__ float3  blend_background(float3  rgb_0, float transmittance_2, float3  background_0)
{
    return clamp_0(rgb_0 + make_float3 (transmittance_2) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S20, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S21, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S22, float3  _S23)
{
    _d_clamp_vector_0(_S20, _S21, _S22, _S23);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dptransmittance_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_1)
{
    float3  _S24 = make_float3 ((*dptransmittance_1).primal_0);
    float3  _S25 = make_float3 (0.0f);
    float3  _S26 = make_float3 (1.0f);
    float3  _S27 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S28;
    (&_S28)->primal_0 = (*dprgb_0).primal_0 + make_float3 ((*dptransmittance_1).primal_0) * (*dpbackground_0).primal_0;
    (&_S28)->differential_0 = _S27;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S29;
    (&_S29)->primal_0 = _S25;
    (&_S29)->differential_0 = _S27;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S30;
    (&_S30)->primal_0 = _S26;
    (&_S30)->differential_0 = _S27;
    s_bwd_prop_clamp_0(&_S28, &_S29, &_S30, _s_dOut_1);
    float3  _S31 = _S24 * _S28.differential_0;
    float3  _S32 = (*dpbackground_0).primal_0 * _S28.differential_0;
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S31;
    float _S33 = _S32.x + _S32.y + _S32.z;
    dptransmittance_1->primal_0 = (*dptransmittance_1).primal_0;
    dptransmittance_1->differential_0 = _S33;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S28.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S34, DiffPair_float_0 * _S35, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S36, float3  _S37)
{
    s_bwd_prop_blend_background_0(_S34, _S35, _S36, _S37);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float transmittance_3, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_transmittance_1, float3  * v_background_0)
{
    float3  _S38 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S38;
    DiffPair_float_0 p_transmittance_1;
    (&p_transmittance_1)->primal_0 = transmittance_3;
    (&p_transmittance_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S38;
    s_bwd_blend_background_0(&p_rgb_0, &p_transmittance_1, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_transmittance_1 = p_transmittance_1.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_2, float dOut_3)
{
    if(((*dpx_3).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_3->primal_0 = (*dpx_3).primal_0;
        dpx_3->differential_0 = 0.0f;
        dpy_2->primal_0 = (*dpy_2).primal_0;
        dpy_2->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_3).primal_0), ((*dpy_2).primal_0)));
        DiffPair_float_0 _S39 = *dpx_3;
        float _S40 = val_0 * (*dpy_2).primal_0 / (*dpx_3).primal_0 * dOut_3;
        dpx_3->primal_0 = (*dpx_3).primal_0;
        dpx_3->differential_0 = _S40;
        float _S41 = val_0 * (F32_log((_S39.primal_0))) * dOut_3;
        dpy_2->primal_0 = (*dpy_2).primal_0;
        dpy_2->differential_0 = _S41;
    }
    return;
}

inline __device__ DiffPair_float_0 _d_pow_1(DiffPair_float_0 * dpx_4, DiffPair_float_0 * dpy_3)
{
    float _S42 = dpx_4->primal_0;
    if((dpx_4->primal_0) < 9.99999997475242708e-07f)
    {
        DiffPair_float_0 _S43 = { 0.0f, 0.0f };
        return _S43;
    }
    float val_1 = (F32_pow((_S42), (dpy_3->primal_0)));
    DiffPair_float_0 _S44 = { val_1, val_1 * (F32_log((_S42))) * dpy_3->differential_0 + val_1 * dpy_3->primal_0 / _S42 * dpx_4->differential_0 };
    return _S44;
}

inline __device__ float linear_rgb_to_srgb(float x_3)
{
    float _S45;
    if(x_3 < 0.00313080009073019f)
    {
        _S45 = x_3 * 12.92000007629394531f;
    }
    else
    {
        _S45 = 1.0549999475479126f * (F32_pow((x_3), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return _S45;
}

inline __device__ float linear_rgb_to_srgb_grad(float x_4)
{
    float _S46;
    if(x_4 < 0.00313080009073019f)
    {
        _S46 = 12.92000007629394531f;
    }
    else
    {
        DiffPair_float_0 _S47;
        (&_S47)->primal_0 = x_4;
        (&_S47)->differential_0 = 1.0f;
        DiffPair_float_0 _S48;
        (&_S48)->primal_0 = 0.4166666567325592f;
        (&_S48)->differential_0 = 0.0f;
        DiffPair_float_0 _S49 = _d_pow_1(&_S47, &_S48);
        _S46 = _S49.differential_0 * 1.0549999475479126f;
    }
    return _S46;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float _S50 = rgb_2.x;
    float _S51;
    if(_S50 < 0.00313080009073019f)
    {
        _S51 = _S50 * 12.92000007629394531f;
    }
    else
    {
        _S51 = 1.0549999475479126f * (F32_pow((_S50), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S52 = rgb_2.y;
    float _S53;
    if(_S52 < 0.00313080009073019f)
    {
        _S53 = _S52 * 12.92000007629394531f;
    }
    else
    {
        _S53 = 1.0549999475479126f * (F32_pow((_S52), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S54 = rgb_2.z;
    float _S55;
    if(_S54 < 0.00313080009073019f)
    {
        _S55 = _S54 * 12.92000007629394531f;
    }
    else
    {
        _S55 = 1.0549999475479126f * (F32_pow((_S54), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S51, _S53, _S55);
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S56, DiffPair_float_0 * _S57, float _S58)
{
    _d_pow_0(_S56, _S57, _S58);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_2)
{
    float _S59 = (*dprgb_1).primal_0.x;
    float _S60 = (*dprgb_1).primal_0.y;
    float _S61 = (*dprgb_1).primal_0.z;
    float _S62;
    if(_S61 < 0.00313080009073019f)
    {
        _S62 = 12.92000007629394531f * _s_dOut_2.z;
    }
    else
    {
        float _S63 = 1.0549999475479126f * _s_dOut_2.z;
        DiffPair_float_0 _S64;
        (&_S64)->primal_0 = _S61;
        (&_S64)->differential_0 = 0.0f;
        DiffPair_float_0 _S65;
        (&_S65)->primal_0 = 0.4166666567325592f;
        (&_S65)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S64, &_S65, _S63);
        _S62 = _S64.differential_0;
    }
    float _S66;
    if(_S60 < 0.00313080009073019f)
    {
        _S66 = 12.92000007629394531f * _s_dOut_2.y;
    }
    else
    {
        float _S67 = 1.0549999475479126f * _s_dOut_2.y;
        DiffPair_float_0 _S68;
        (&_S68)->primal_0 = _S60;
        (&_S68)->differential_0 = 0.0f;
        DiffPair_float_0 _S69;
        (&_S69)->primal_0 = 0.4166666567325592f;
        (&_S69)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S68, &_S69, _S67);
        _S66 = _S68.differential_0;
    }
    float _S70;
    if(_S59 < 0.00313080009073019f)
    {
        _S70 = 12.92000007629394531f * _s_dOut_2.x;
    }
    else
    {
        float _S71 = 1.0549999475479126f * _s_dOut_2.x;
        DiffPair_float_0 _S72;
        (&_S72)->primal_0 = _S59;
        (&_S72)->differential_0 = 0.0f;
        DiffPair_float_0 _S73;
        (&_S73)->primal_0 = 0.4166666567325592f;
        (&_S73)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S72, &_S73, _S71);
        _S70 = _S72.differential_0;
    }
    float3  _S74 = make_float3 (_S70, _S66, _S62);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S74;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S75, float3  _S76)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S75, _S76);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S77 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S77;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_5)
{
    Matrix<float, 2, 2>  result_2;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(2))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_0)), c_0) = _slang_vector_get_element(x_5.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_2;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_5, float dOut_4)
{
    float _S78 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_5).primal_0)))))) * dOut_4;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S78;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, float dOut_5)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_4).primal_0.x * dOut_5;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_6).primal_0.x * dOut_5;
    *&((&x_d_result_0)->y) = (*dpy_4).primal_0.y * dOut_5;
    *&((&y_d_result_0)->y) = (*dpx_6).primal_0.y * dOut_5;
    *&((&x_d_result_0)->z) = (*dpy_4).primal_0.z * dOut_5;
    *&((&y_d_result_0)->z) = (*dpx_6).primal_0.z * dOut_5;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = x_d_result_0;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_6, float3  y_2)
{
    int i_2 = int(0);
    float result_3 = 0.0f;
    for(;;)
    {
        if(i_2 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_4 = result_3 + _slang_vector_get_element(x_6, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_3 = result_4;
    }
    return result_3;
}

inline __device__ float dot_1(float2  x_7, float2  y_3)
{
    int i_3 = int(0);
    float result_5 = 0.0f;
    for(;;)
    {
        if(i_3 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_6 = result_5 + _slang_vector_get_element(x_7, i_3) * _slang_vector_get_element(y_3, i_3);
        i_3 = i_3 + int(1);
        result_5 = result_6;
    }
    return result_5;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_1(x_8, x_8))));
}

inline __device__ float length_1(float3  x_9)
{
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ bool undistort_point_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0, int maxiter_0, float2  * uv_undist_0)
{
    int i_4 = int(0);
    float2  q_0 = uv_0;
    for(;;)
    {
        if(i_4 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S79 = (*dist_coeffs_0)[int(3)];
        float _S80 = (*dist_coeffs_0)[int(4)];
        float _S81 = (*dist_coeffs_0)[int(5)];
        float _S82 = (*dist_coeffs_0)[int(6)];
        float _S83 = (*dist_coeffs_0)[int(7)];
        float _S84 = (*dist_coeffs_0)[int(8)];
        float _S85 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S86 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S87 = (*dist_coeffs_0)[int(1)] + r2_0 * _S86;
        float _S88 = (*dist_coeffs_0)[int(0)] + r2_0 * _S87;
        float radial_0 = 1.0f + r2_0 * _S88;
        float _S89 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S90 = _S89 * u_0;
        float _S91 = 2.0f * u_0;
        float _S92 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S93 = _S92 * u_0;
        float _S94 = 2.0f * v_0;
        float2  _S95 = q_0 * make_float2 (radial_0) + make_float2 (_S90 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S91 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S93 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S94 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S95 + make_float2 ((*dist_coeffs_0)[int(8)] * _S95.x + (*dist_coeffs_0)[int(9)] * _S95.y, 0.0f) - uv_0;
        float _S96 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S96 + _S96);
        float2  _S97 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S88 + (s_diff_r2_0 * _S87 + (s_diff_r2_0 * _S86 + s_diff_r2_0 * _S79 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S89 * v_0 + 0.0f * _S90 + (s_diff_r2_0 + (_S91 + _S91)) * _S81 + s_diff_r2_0 * _S82, _S92 * v_0 + 0.0f * _S93 + (s_diff_r2_0 + (_S96 + 0.0f * _S94)) * _S80 + s_diff_r2_0 * _S83);
        float _S98 = 0.0f * u_0;
        float s_diff_r2_1 = _S98 + _S98 + (v_0 + v_0);
        float2  _S99 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S88 + (s_diff_r2_1 * _S87 + (s_diff_r2_1 * _S86 + s_diff_r2_1 * _S79 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S89 * v_0 + _S90 + (s_diff_r2_1 + (_S98 + 0.0f * _S91)) * _S81 + s_diff_r2_1 * _S82, 0.0f * _S92 * v_0 + _S93 + (s_diff_r2_1 + (_S94 + _S94)) * _S80 + s_diff_r2_1 * _S83);
        Matrix<float, 2, 2>  _S100 = transpose_0(makeMatrix<float, 2, 2> (_S97 + make_float2 (_S97.x * _S84 + _S97.y * _S85, 0.0f), _S99 + make_float2 (_S99.x * _S84 + _S99.y * _S85, 0.0f)));
        float inv_det_0 = 1.0f / (_S100.rows[int(0)].x * _S100.rows[int(1)].y - _S100.rows[int(0)].y * _S100.rows[int(1)].x);
        float _S101 = r_1.x;
        float _S102 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S101 * _S100.rows[int(1)].y - _S102 * _S100.rows[int(0)].y) * inv_det_0, (- _S101 * _S100.rows[int(1)].x + _S102 * _S100.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S103 = (*dist_coeffs_0)[int(0)];
    float _S104 = (*dist_coeffs_0)[int(1)];
    float _S105 = (*dist_coeffs_0)[int(2)];
    float _S106 = (*dist_coeffs_0)[int(3)];
    float _S107 = (*dist_coeffs_0)[int(4)];
    float _S108 = (*dist_coeffs_0)[int(5)];
    float _S109 = (*dist_coeffs_0)[int(6)];
    float _S110 = (*dist_coeffs_0)[int(7)];
    float _S111 = (*dist_coeffs_0)[int(8)];
    float _S112 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S113 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S113 + _S113);
    float _S114 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S115 = (*dist_coeffs_0)[int(1)] + r2_1 * _S114;
    float _S116 = (*dist_coeffs_0)[int(0)] + r2_1 * _S115;
    float radial_1 = 1.0f + r2_1 * _S116;
    float _S117 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S118 = _S117 * u_1;
    float _S119 = 2.0f * u_1;
    float _S120 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S121 = _S120 * u_1;
    float _S122 = 2.0f * v_1;
    float2  _S123 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S116 + (s_diff_r2_2 * _S115 + (s_diff_r2_2 * _S114 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S117 * v_1 + 0.0f * _S118 + (s_diff_r2_2 + (_S119 + _S119)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S120 * v_1 + 0.0f * _S121 + (s_diff_r2_2 + (_S113 + 0.0f * _S122)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S124 = 0.0f * u_1;
    float s_diff_r2_3 = _S124 + _S124 + (v_1 + v_1);
    float2  _S125 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S116 + (s_diff_r2_3 * _S115 + (s_diff_r2_3 * _S114 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S117 * v_1 + _S118 + (s_diff_r2_3 + (_S124 + 0.0f * _S119)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S120 * v_1 + _S121 + (s_diff_r2_3 + (_S122 + _S122)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S126 = transpose_0(makeMatrix<float, 2, 2> (_S123 + make_float2 (_S123.x * (*dist_coeffs_0)[int(8)] + _S123.y * (*dist_coeffs_0)[int(9)], 0.0f), _S125 + make_float2 (_S125.x * (*dist_coeffs_0)[int(8)] + _S125.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S127;
    if((F32_min((determinant_0(_S126)), ((F32_min((_S126.rows[int(0)].x), (_S126.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S128 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S103 + r2_2 * (_S104 + r2_2 * (_S105 + r2_2 * _S106)))) + make_float2 (_S117 * u_2 * v_2 + _S108 * (r2_2 + 2.0f * u_2 * u_2) + _S109 * r2_2, _S120 * u_2 * v_2 + _S107 * (r2_2 + 2.0f * v_2 * v_2) + _S110 * r2_2);
        _S127 = (length_0(_S128 + make_float2 (_S111 * _S128.x + _S112 * _S128.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S127 = false;
    }
    return _S127;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S129 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S129;
    FixedArray<float, 10>  _S130 = dist_coeffs_1;
    bool _S131 = undistort_point_0(_S129, &_S130, int(12), &uv_1);
    if(!_S131)
    {
        int3  _S132 = make_int3 (int(0));
        float3  _S133 = make_float3 ((float)_S132.x, (float)_S132.y, (float)_S132.z);
        return _S133;
    }
    float3  raydir_0;
    if(is_fisheye_0)
    {
        float theta_0 = length_0(uv_1);
        float3  raydir_1 = make_float3 ((uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).x, (uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).y, (F32_cos((theta_0))));
        if(!is_ray_depth_0)
        {
            raydir_0 = raydir_1 / make_float3 (raydir_1.z);
        }
        else
        {
            raydir_0 = raydir_1;
        }
    }
    else
    {
        float3  raydir_2 = make_float3 (uv_1.x, uv_1.y, 1.0f);
        if(is_ray_depth_0)
        {
            raydir_0 = normalize_0(raydir_2);
        }
        else
        {
            raydir_0 = raydir_2;
        }
    }
    return raydir_0;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_6)
{
    float _S134 = dOut_6.y;
    float _S135 = dOut_6.z;
    float _S136 = dOut_6.x;
    float _S137 = (*a_0).primal_0.z * _S134 + - (*a_0).primal_0.y * _S135;
    float _S138 = - (*a_0).primal_0.z * _S136 + (*a_0).primal_0.x * _S135;
    float _S139 = (*a_0).primal_0.y * _S136 + - (*a_0).primal_0.x * _S134;
    float3  _S140 = make_float3 (- (*b_0).primal_0.z * _S134 + (*b_0).primal_0.y * _S135, (*b_0).primal_0.z * _S136 + - (*b_0).primal_0.x * _S135, - (*b_0).primal_0.y * _S136 + (*b_0).primal_0.x * _S134);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S140;
    float3  _S141 = make_float3 (_S137, _S138, _S139);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S141;
    return;
}

inline __device__ float3  cross_0(float3  left_0, float3  right_0)
{
    float _S142 = left_0.y;
    float _S143 = right_0.z;
    float _S144 = left_0.z;
    float _S145 = right_0.y;
    float _S146 = right_0.x;
    float _S147 = left_0.x;
    return make_float3 (_S142 * _S143 - _S144 * _S145, _S144 * _S146 - _S147 * _S143, _S147 * _S145 - _S142 * _S146);
}

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  points_0)
{
    float3  normal_0 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    float3  normal_1;
    if((dot_0(normal_0, normal_0)) != 0.0f)
    {
        normal_1 = normal_0 / make_float3 (length_1(normal_0));
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S148, float3  _S149)
{
    return cross_0(_S148, _S149);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S150, float3  _S151)
{
    return dot_0(_S150, _S151);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S152, float _S153)
{
    _d_sqrt_0(_S152, _S153);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_3)
{
    float _S154 = (*dpx_7).primal_0.x;
    float _S155 = (*dpx_7).primal_0.y;
    float _S156 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S157;
    (&_S157)->primal_0 = _S154 * _S154 + _S155 * _S155 + _S156 * _S156;
    (&_S157)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S157, _s_dOut_3);
    float _S158 = (*dpx_7).primal_0.z * _S157.differential_0;
    float _S159 = _S158 + _S158;
    float _S160 = (*dpx_7).primal_0.y * _S157.differential_0;
    float _S161 = _S160 + _S160;
    float _S162 = (*dpx_7).primal_0.x * _S157.differential_0;
    float _S163 = _S162 + _S162;
    float3  _S164 = make_float3 (0.0f);
    *&((&_S164)->z) = _S159;
    *&((&_S164)->y) = _S161;
    *&((&_S164)->x) = _S163;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S164;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S165, float _S166)
{
    s_bwd_prop_length_impl_0(_S165, _S166);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S167, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S168, float _S169)
{
    _d_dot_0(_S167, _S168, _S169);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S170, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S171, float3  _S172)
{
    _d_cross_0(_S170, _S171, _S172);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_4)
{
    float3  _S173 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S174 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S175 = s_primal_ctx_cross_0(dx_0, _S174);
    bool _S176 = (s_primal_ctx_dot_0(_S175, _S175)) != 0.0f;
    float3  _S177;
    float3  _S178;
    if(_S176)
    {
        float _S179 = length_1(_S175);
        float3  _S180 = make_float3 (_S179);
        _S177 = make_float3 (_S179 * _S179);
        _S178 = _S180;
    }
    else
    {
        _S177 = _S173;
        _S178 = _S173;
    }
    if(_S176)
    {
        float3  _S181 = _s_dOut_4 / _S177;
        float3  _S182 = _S175 * - _S181;
        float3  _S183 = _S178 * _S181;
        float _S184 = _S182.x + _S182.y + _S182.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S185;
        (&_S185)->primal_0 = _S175;
        (&_S185)->differential_0 = _S173;
        s_bwd_length_impl_0(&_S185, _S184);
        _S177 = _S183 + _S185.differential_0;
    }
    else
    {
        _S177 = _s_dOut_4;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S186;
    (&_S186)->primal_0 = _S175;
    (&_S186)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S187;
    (&_S187)->primal_0 = _S175;
    (&_S187)->differential_0 = _S173;
    s_bwd_prop_dot_0(&_S186, &_S187, 0.0f);
    float3  _S188 = _S187.differential_0 + _S186.differential_0 + _S177;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S189;
    (&_S189)->primal_0 = dx_0;
    (&_S189)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S190;
    (&_S190)->primal_0 = _S174;
    (&_S190)->differential_0 = _S173;
    s_bwd_prop_cross_0(&_S189, &_S190, _S188);
    float3  s_diff_dy_T_0 = - _S190.differential_0;
    float3  _S191 = - s_diff_dy_T_0;
    float3  _S192 = - _S189.differential_0;
    FixedArray<float3 , 4>  _S193;
    _S193[int(0)] = _S173;
    _S193[int(1)] = _S173;
    _S193[int(2)] = _S173;
    _S193[int(3)] = _S173;
    _S193[int(2)] = _S191;
    _S193[int(3)] = s_diff_dy_T_0;
    _S193[int(0)] = _S192;
    _S193[int(1)] = _S189.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S193;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S194, float3  _S195)
{
    s_bwd_prop_points_to_normal_0(_S194, _S195);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S196 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S196;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S197 = float2 {intrins_1.z, intrins_1.w};
    float2  _S198 = float2 {intrins_1.x, intrins_1.y};
    float2  _S199 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S197) / _S198;
    float2  uv_2 = _S199;
    FixedArray<float, 10>  _S200 = dist_coeffs_2;
    bool _S201 = undistort_point_0(_S199, &_S200, int(12), &uv_2);
    if(!_S201)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_3;
    if(is_fisheye_1)
    {
        float theta_1 = length_0(uv_2);
        float3  raydir_4 = make_float3 ((uv_2 / make_float2 ((F32_max((theta_1), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_1))))).x, (uv_2 / make_float2 ((F32_max((theta_1), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_1))))).y, (F32_cos((theta_1))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_4 / make_float3 (raydir_4.z);
        }
        else
        {
            raydir_3 = raydir_4;
        }
    }
    else
    {
        float3  raydir_5 = make_float3 (uv_2.x, uv_2.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_0(raydir_5);
        }
        else
        {
            raydir_3 = raydir_5;
        }
    }
    points_2[int(0)] = make_float3 (depths_0.x) * raydir_3;
    float2  _S202 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S197) / _S198;
    float2  uv_3 = _S202;
    FixedArray<float, 10>  _S203 = dist_coeffs_2;
    bool _S204 = undistort_point_0(_S202, &_S203, int(12), &uv_3);
    if(!_S204)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_2 = length_0(uv_3);
        float3  raydir_6 = make_float3 ((uv_3 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).x, (uv_3 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).y, (F32_cos((theta_2))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_6 / make_float3 (raydir_6.z);
        }
        else
        {
            raydir_3 = raydir_6;
        }
    }
    else
    {
        float3  raydir_7 = make_float3 (uv_3.x, uv_3.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_0(raydir_7);
        }
        else
        {
            raydir_3 = raydir_7;
        }
    }
    points_2[int(1)] = make_float3 (depths_0.y) * raydir_3;
    float2  _S205 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S197) / _S198;
    float2  uv_4 = _S205;
    FixedArray<float, 10>  _S206 = dist_coeffs_2;
    bool _S207 = undistort_point_0(_S205, &_S206, int(12), &uv_4);
    if(!_S207)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_3 = length_0(uv_4);
        float3  raydir_8 = make_float3 ((uv_4 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).x, (uv_4 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).y, (F32_cos((theta_3))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_8 / make_float3 (raydir_8.z);
        }
        else
        {
            raydir_3 = raydir_8;
        }
    }
    else
    {
        float3  raydir_9 = make_float3 (uv_4.x, uv_4.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_0(raydir_9);
        }
        else
        {
            raydir_3 = raydir_9;
        }
    }
    points_2[int(2)] = make_float3 (depths_0.z) * raydir_3;
    float2  _S208 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S197) / _S198;
    float2  uv_5 = _S208;
    FixedArray<float, 10>  _S209 = dist_coeffs_2;
    bool _S210 = undistort_point_0(_S208, &_S209, int(12), &uv_5);
    if(!_S210)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_1)
    {
        float theta_4 = length_0(uv_5);
        float3  raydir_10 = make_float3 ((uv_5 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_5 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_1)
        {
            raydir_3 = raydir_10 / make_float3 (raydir_10.z);
        }
        else
        {
            raydir_3 = raydir_10;
        }
    }
    else
    {
        float3  raydir_11 = make_float3 (uv_5.x, uv_5.y, 1.0f);
        if(is_ray_depth_1)
        {
            raydir_3 = normalize_0(raydir_11);
        }
        else
        {
            raydir_3 = raydir_11;
        }
    }
    points_2[int(3)] = make_float3 (depths_0.w) * raydir_3;
    float3  normal_2 = cross_0(points_2[int(1)] - points_2[int(0)], - (points_2[int(3)] - points_2[int(2)]));
    float3  normal_3;
    if((dot_0(normal_2, normal_2)) != 0.0f)
    {
        normal_3 = normal_2 / make_float3 (length_1(normal_2));
    }
    else
    {
        normal_3 = normal_2;
    }
    return normal_3;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S211;
    bool _S212;
    float2  _S213;
    bool _S214;
    float2  _S215;
    bool _S216;
    float2  _S217;
    bool _S218;
};

inline __device__ float s_primal_ctx_sin_0(float _S219)
{
    return (F32_sin((_S219)));
}

inline __device__ float s_primal_ctx_cos_0(float _S220)
{
    return (F32_cos((_S220)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S221 = make_float2 (0.0f);
    _s_diff_ctx_0->_S211 = _S221;
    _s_diff_ctx_0->_S212 = false;
    _s_diff_ctx_0->_S213 = _S221;
    _s_diff_ctx_0->_S214 = false;
    _s_diff_ctx_0->_S215 = _S221;
    _s_diff_ctx_0->_S216 = false;
    _s_diff_ctx_0->_S217 = _S221;
    _s_diff_ctx_0->_S218 = false;
    _s_diff_ctx_0->_S213 = _S221;
    _s_diff_ctx_0->_S214 = false;
    _s_diff_ctx_0->_S215 = _S221;
    _s_diff_ctx_0->_S216 = false;
    _s_diff_ctx_0->_S217 = _S221;
    _s_diff_ctx_0->_S218 = false;
    float3  _S222 = make_float3 (0.0f);
    float2  _S223 = float2 {intrins_2.z, intrins_2.w};
    float2  _S224 = float2 {intrins_2.x, intrins_2.y};
    float2  _S225 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S223) / _S224;
    float2  _S226 = _S225;
    bool _S227 = undistort_point_0(_S225, dist_coeffs_3, int(12), &_S226);
    _s_diff_ctx_0->_S211 = _S226;
    _s_diff_ctx_0->_S212 = _S227;
    float2  uv_6 = _S226;
    bool _S228 = !_S227;
    float3  normal_4;
    if(_S228)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S229 = !_S228;
    int _S230;
    FixedArray<float3 , 4>  points_3;
    if(_S229)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S231 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S231), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S231))).x, (uv_6 / make_float2 ((F32_max((_S231), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S231))).y, s_primal_ctx_cos_0(_S231));
            if(!is_ray_depth_2)
            {
                raydir_12 = raydir_13 / make_float3 (raydir_13.z);
            }
            else
            {
                raydir_12 = raydir_13;
            }
        }
        else
        {
            float3  raydir_14 = make_float3 (uv_6.x, uv_6.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_12 = normalize_0(raydir_14);
            }
            else
            {
                raydir_12 = raydir_14;
            }
        }
        float3  _S232 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S233 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S223) / _S224;
        float2  _S234 = _S233;
        bool _S235 = undistort_point_0(_S233, dist_coeffs_3, int(12), &_S234);
        _s_diff_ctx_0->_S213 = _S234;
        _s_diff_ctx_0->_S214 = _S235;
        float2  uv_7 = _S234;
        bool _S236 = !_S235;
        if(_S236)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S237 = !_S236;
        if(_S237)
        {
            if(is_fisheye_2)
            {
                float _S238 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S238), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S238))).x, (uv_7 / make_float2 ((F32_max((_S238), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S238))).y, s_primal_ctx_cos_0(_S238));
                if(!is_ray_depth_2)
                {
                    raydir_12 = raydir_15 / make_float3 (raydir_15.z);
                }
                else
                {
                    raydir_12 = raydir_15;
                }
            }
            else
            {
                float3  raydir_16 = make_float3 (uv_7.x, uv_7.y, 1.0f);
                if(is_ray_depth_2)
                {
                    raydir_12 = normalize_0(raydir_16);
                }
                else
                {
                    raydir_12 = raydir_16;
                }
            }
            float3  _S239 = make_float3 (dpdepths_0.y) * raydir_12;
            _S230 = int(2);
            points_3[int(0)] = _S232;
            points_3[int(1)] = _S239;
            points_3[int(2)] = _S222;
            points_3[int(3)] = _S222;
        }
        else
        {
            _S230 = int(0);
            points_3[int(0)] = _S232;
            points_3[int(1)] = _S222;
            points_3[int(2)] = _S222;
            points_3[int(3)] = _S222;
        }
        bool _runFlag_0;
        if(_S230 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S229;
            _S230 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S240 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S223) / _S224;
            float2  _S241 = _S240;
            bool _S242 = undistort_point_0(_S240, dist_coeffs_3, int(12), &_S241);
            _s_diff_ctx_0->_S215 = _S241;
            _s_diff_ctx_0->_S216 = _S242;
            float2  uv_8 = _S241;
            if(!_S242)
            {
                float3  _S243 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S230 = int(0);
                normal_4 = _S243;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S244 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S244), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S244))).x, (uv_8 / make_float2 ((F32_max((_S244), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S244))).y, s_primal_ctx_cos_0(_S244));
                    if(!is_ray_depth_2)
                    {
                        raydir_12 = raydir_17 / make_float3 (raydir_17.z);
                    }
                    else
                    {
                        raydir_12 = raydir_17;
                    }
                }
                else
                {
                    float3  raydir_18 = make_float3 (uv_8.x, uv_8.y, 1.0f);
                    if(is_ray_depth_2)
                    {
                        raydir_12 = normalize_0(raydir_18);
                    }
                    else
                    {
                        raydir_12 = raydir_18;
                    }
                }
                points_3[int(2)] = make_float3 (dpdepths_0.z) * raydir_12;
                float2  _S245 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S223) / _S224;
                float2  _S246 = _S245;
                bool _S247 = undistort_point_0(_S245, dist_coeffs_3, int(12), &_S246);
                _s_diff_ctx_0->_S217 = _S246;
                _s_diff_ctx_0->_S218 = _S247;
                float2  uv_9 = _S246;
                bool _S248 = !_S247;
                if(_S248)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S249 = !_S248;
                int _S250;
                if(_S249)
                {
                    if(is_fisheye_2)
                    {
                        float _S251 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S251), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S251))).x, (uv_9 / make_float2 ((F32_max((_S251), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S251))).y, s_primal_ctx_cos_0(_S251));
                        if(!is_ray_depth_2)
                        {
                            raydir_12 = raydir_19 / make_float3 (raydir_19.z);
                        }
                        else
                        {
                            raydir_12 = raydir_19;
                        }
                    }
                    else
                    {
                        float3  raydir_20 = make_float3 (uv_9.x, uv_9.y, 1.0f);
                        if(is_ray_depth_2)
                        {
                            raydir_12 = normalize_0(raydir_20);
                        }
                        else
                        {
                            raydir_12 = raydir_20;
                        }
                    }
                    points_3[int(3)] = make_float3 (dpdepths_0.w) * raydir_12;
                    _S250 = int(2);
                }
                else
                {
                    _S250 = int(0);
                }
                if(_S250 != int(2))
                {
                    _runFlag_0 = false;
                    _S230 = _S250;
                }
                if(_runFlag_0)
                {
                    _S230 = int(1);
                }
            }
        }
    }
    else
    {
        _S230 = int(0);
        points_3[int(0)] = _S222;
        points_3[int(1)] = _S222;
        points_3[int(2)] = _S222;
        points_3[int(3)] = _S222;
    }
    if(!(_S230 != int(1)))
    {
        float3  _S252 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S252, _S252)) != 0.0f)
        {
            normal_4 = _S252 / make_float3 (length_1(_S252));
        }
        else
        {
            normal_4 = _S252;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_5, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S253 = *dpdepths_1;
    float3  _S254 = make_float3 (0.0f);
    float2  _S255 = _s_diff_ctx_1->_S211;
    bool _S256 = !!_s_diff_ctx_1->_S212;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S257;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S258;
    if(_S256)
    {
        if(is_fisheye_3)
        {
            float _S259 = length_0(_S255);
            float3  raydir_25 = make_float3 ((_S255 / make_float2 ((F32_max((_S259), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S259))).x, (_S255 / make_float2 ((F32_max((_S259), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S259))).y, s_primal_ctx_cos_0(_S259));
            if(!is_ray_depth_3)
            {
                raydir_21 = raydir_25 / make_float3 (raydir_25.z);
            }
            else
            {
                raydir_21 = raydir_25;
            }
        }
        else
        {
            float3  raydir_26 = make_float3 (_S255.x, _S255.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S260 = make_float3 (_S253.primal_0.x) * raydir_21;
        float2  _S261 = _s_diff_ctx_1->_S213;
        bool _S262 = !!_s_diff_ctx_1->_S214;
        if(_S262)
        {
            if(is_fisheye_3)
            {
                float _S263 = length_0(_S261);
                float3  raydir_27 = make_float3 ((_S261 / make_float2 ((F32_max((_S263), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S263))).x, (_S261 / make_float2 ((F32_max((_S263), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S263))).y, s_primal_ctx_cos_0(_S263));
                if(!is_ray_depth_3)
                {
                    raydir_22 = raydir_27 / make_float3 (raydir_27.z);
                }
                else
                {
                    raydir_22 = raydir_27;
                }
            }
            else
            {
                float3  raydir_28 = make_float3 (_S261.x, _S261.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S264 = make_float3 (_S253.primal_0.y) * raydir_22;
            _S257 = int(2);
            points_4[int(0)] = _S260;
            points_4[int(1)] = _S264;
            points_4[int(2)] = _S254;
            points_4[int(3)] = _S254;
        }
        else
        {
            _S257 = int(0);
            points_4[int(0)] = _S260;
            points_4[int(1)] = _S254;
            points_4[int(2)] = _S254;
            points_4[int(3)] = _S254;
            raydir_22 = _S254;
        }
        if(_S257 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S256;
            _S257 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S265 = _s_diff_ctx_1->_S215;
            if(!_s_diff_ctx_1->_S216)
            {
                _runFlag_2 = false;
                _S257 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S266 = length_0(_S265);
                    float3  raydir_29 = make_float3 ((_S265 / make_float2 ((F32_max((_S266), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S266))).x, (_S265 / make_float2 ((F32_max((_S266), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S266))).y, s_primal_ctx_cos_0(_S266));
                    if(!is_ray_depth_3)
                    {
                        raydir_23 = raydir_29 / make_float3 (raydir_29.z);
                    }
                    else
                    {
                        raydir_23 = raydir_29;
                    }
                }
                else
                {
                    float3  raydir_30 = make_float3 (_S265.x, _S265.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S253.primal_0.z) * raydir_23;
                float2  _S267 = _s_diff_ctx_1->_S217;
                bool _S268 = !!_s_diff_ctx_1->_S218;
                int _S269;
                if(_S268)
                {
                    if(is_fisheye_3)
                    {
                        float _S270 = length_0(_S267);
                        float3  raydir_31 = make_float3 ((_S267 / make_float2 ((F32_max((_S270), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S270))).x, (_S267 / make_float2 ((F32_max((_S270), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S270))).y, s_primal_ctx_cos_0(_S270));
                        if(!is_ray_depth_3)
                        {
                            raydir_24 = raydir_31 / make_float3 (raydir_31.z);
                        }
                        else
                        {
                            raydir_24 = raydir_31;
                        }
                    }
                    else
                    {
                        float3  raydir_32 = make_float3 (_S267.x, _S267.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S253.primal_0.w) * raydir_24;
                    _S269 = int(2);
                }
                else
                {
                    _S269 = int(0);
                    raydir_24 = _S254;
                }
                if(_S269 != int(2))
                {
                    _runFlag_3 = false;
                    _S257 = _S269;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S257 = int(1);
                }
                float3  _S271 = raydir_23;
                _runFlag_3 = _S268;
                raydir_23 = raydir_24;
                raydir_24 = _S271;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S254;
                raydir_24 = _S254;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S254;
            raydir_24 = _S254;
        }
        float3  _S272 = raydir_21;
        float3  _S273 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S258 = _S262;
        raydir_23 = _S273;
        raydir_24 = _S272;
    }
    else
    {
        _S257 = int(0);
        points_4[int(0)] = _S254;
        points_4[int(1)] = _S254;
        points_4[int(2)] = _S254;
        points_4[int(3)] = _S254;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S254;
        raydir_22 = _S254;
        _S258 = false;
        raydir_23 = _S254;
        raydir_24 = _S254;
    }
    bool _S274 = !(_S257 != int(1));
    float3  _S275;
    float3  _S276;
    float3  _S277;
    float3  _S278;
    float3  _S279;
    bool _S280;
    if(_S274)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S281 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S282 = s_primal_ctx_cross_0(dx_1, _S281);
        bool _S283 = (s_primal_ctx_dot_0(_S282, _S282)) != 0.0f;
        if(_S283)
        {
            float _S284 = length_1(_S282);
            float3  _S285 = make_float3 (_S284);
            _S275 = make_float3 (_S284 * _S284);
            _S276 = _S285;
        }
        else
        {
            _S275 = _S254;
            _S276 = _S254;
        }
        float3  _S286 = _S276;
        _S280 = _S283;
        _S276 = _S282;
        _S277 = _S286;
        _S278 = dx_1;
        _S279 = _S281;
    }
    else
    {
        _S280 = false;
        _S275 = _S254;
        _S276 = _S254;
        _S277 = _S254;
        _S278 = _S254;
        _S279 = _S254;
    }
    float4  _S287 = make_float4 (0.0f);
    if(_S274)
    {
        if(_S280)
        {
            float3  _S288 = _s_dOut_5 / _S275;
            float3  _S289 = _S276 * - _S288;
            float3  _S290 = _S277 * _S288;
            float _S291 = _S289.x + _S289.y + _S289.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S292;
            (&_S292)->primal_0 = _S276;
            (&_S292)->differential_0 = _S254;
            s_bwd_length_impl_0(&_S292, _S291);
            _S275 = _S290 + _S292.differential_0;
        }
        else
        {
            _S275 = _s_dOut_5;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S293;
        (&_S293)->primal_0 = _S276;
        (&_S293)->differential_0 = _S254;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S294;
        (&_S294)->primal_0 = _S276;
        (&_S294)->differential_0 = _S254;
        s_bwd_prop_dot_0(&_S293, &_S294, 0.0f);
        float3  _S295 = _S294.differential_0 + _S293.differential_0 + _S275;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S296;
        (&_S296)->primal_0 = _S278;
        (&_S296)->differential_0 = _S254;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S297;
        (&_S297)->primal_0 = _S279;
        (&_S297)->differential_0 = _S254;
        s_bwd_prop_cross_0(&_S296, &_S297, _S295);
        float3  s_diff_dy_T_1 = - _S297.differential_0;
        float3  _S298 = - s_diff_dy_T_1;
        float3  _S299 = - _S296.differential_0;
        FixedArray<float3 , 4>  _S300;
        _S300[int(0)] = _S254;
        _S300[int(1)] = _S254;
        _S300[int(2)] = _S254;
        _S300[int(3)] = _S254;
        _S300[int(2)] = _S298;
        _S300[int(3)] = s_diff_dy_T_1;
        _S300[int(0)] = _S299;
        _S300[int(1)] = _S296.differential_0;
        points_4[int(0)] = _S300[int(0)];
        points_4[int(1)] = _S300[int(1)];
        points_4[int(2)] = _S300[int(2)];
        points_4[int(3)] = _S300[int(3)];
    }
    else
    {
        points_4[int(0)] = _S254;
        points_4[int(1)] = _S254;
        points_4[int(2)] = _S254;
        points_4[int(3)] = _S254;
    }
    float4  _S301;
    if(_S256)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S302 = points_4;
                FixedArray<float3 , 4>  _S303 = points_4;
                FixedArray<float3 , 4>  _S304 = points_4;
                FixedArray<float3 , 4>  _S305 = points_4;
                if(_runFlag_3)
                {
                    float3  _S306 = raydir_21 * _S305[int(3)];
                    float _S307 = _S306.x + _S306.y + _S306.z;
                    float4  _S308 = _S287;
                    *&((&_S308)->w) = _S307;
                    points_4[int(0)] = _S302[int(0)];
                    points_4[int(1)] = _S303[int(1)];
                    points_4[int(2)] = _S304[int(2)];
                    points_4[int(3)] = _S254;
                    _S301 = _S308;
                }
                else
                {
                    points_4[int(0)] = _S302[int(0)];
                    points_4[int(1)] = _S303[int(1)];
                    points_4[int(2)] = _S304[int(2)];
                    points_4[int(3)] = _S305[int(3)];
                    _S301 = _S287;
                }
                float3  _S309 = raydir_22 * points_4[int(2)];
                float _S310 = _S309.x + _S309.y + _S309.z;
                FixedArray<float3 , 4>  _S311 = points_4;
                FixedArray<float3 , 4>  _S312 = points_4;
                float4  _S313 = _S287;
                *&((&_S313)->z) = _S310;
                float4  _S314 = _S301 + _S313;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S311[int(1)];
                points_4[int(2)] = _S254;
                points_4[int(3)] = _S312[int(3)];
                _S301 = _S314;
            }
            else
            {
                FixedArray<float3 , 4>  _S315 = points_4;
                FixedArray<float3 , 4>  _S316 = points_4;
                FixedArray<float3 , 4>  _S317 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S315[int(1)];
                points_4[int(2)] = _S316[int(2)];
                points_4[int(3)] = _S317[int(3)];
                _S301 = _S287;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S318 = points_4;
            FixedArray<float3 , 4>  _S319 = points_4;
            FixedArray<float3 , 4>  _S320 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S318[int(1)];
            points_4[int(2)] = _S319[int(2)];
            points_4[int(3)] = _S320[int(3)];
            _S301 = _S287;
        }
        if(_S258)
        {
            FixedArray<float3 , 4>  _S321 = points_4;
            float3  _S322 = raydir_23 * points_4[int(1)];
            float _S323 = _S322.x + _S322.y + _S322.z;
            float4  _S324 = _S287;
            *&((&_S324)->y) = _S323;
            float4  _S325 = _S301 + _S324;
            points_4[int(0)] = _S254;
            points_4[int(1)] = _S254;
            points_4[int(2)] = _S254;
            points_4[int(3)] = _S254;
            raydir_21 = _S321[int(0)];
            _S301 = _S325;
        }
        else
        {
            FixedArray<float3 , 4>  _S326 = points_4;
            FixedArray<float3 , 4>  _S327 = points_4;
            FixedArray<float3 , 4>  _S328 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S326[int(1)];
            points_4[int(2)] = _S327[int(2)];
            points_4[int(3)] = _S328[int(3)];
            raydir_21 = _S254;
        }
        float3  _S329 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S330 = _S329.x + _S329.y + _S329.z;
        float4  _S331 = _S287;
        *&((&_S331)->x) = _S330;
        _S301 = _S301 + _S331;
    }
    else
    {
        _S301 = _S287;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S301;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S332, float4  _S333, FixedArray<float, 10>  * _S334, bool _S335, bool _S336, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S337, float3  _S338)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S339;
    float3  _S340 = s_primal_ctx_depth_to_normal_0(_S332, _S333, _S334, _S335, _S336, (*_S337).primal_0, &_S339);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S341 = _S339;
    s_bwd_prop_depth_to_normal_0(_S332, _S333, _S334, _S335, _S336, _S337, _S338, &_S341);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S342 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S342;
    FixedArray<float, 10>  _S343 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S343, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S344 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S344;
    FixedArray<float, 10>  _S345 = dist_coeffs_6;
    bool _S346 = undistort_point_0(_S344, &_S345, int(12), &uv_10);
    if(!_S346)
    {
        return 0.0f;
    }
    float3  raydir_33;
    if(is_fisheye_5)
    {
        float theta_5 = length_0(uv_10);
        float3  raydir_34 = make_float3 ((uv_10 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_10 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        raydir_33 = raydir_34 / make_float3 (raydir_34.z);
    }
    else
    {
        raydir_33 = make_float3 (uv_10.x, uv_10.y, 1.0f);
    }
    return float((F32_sign((raydir_33.z)))) / length_1(raydir_33);
}

