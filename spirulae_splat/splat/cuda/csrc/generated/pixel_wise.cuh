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

inline __device__ float rendered_depth_to_expected_depth(float depth_0, float alpha_0)
{
    return depth_0 / (F32_max((alpha_0), (1.00000001335143196e-10f)));
}

inline __device__ void s_bwd_prop_rendered_depth_to_expected_depth_0(DiffPair_float_0 * dpdepth_0, DiffPair_float_0 * dpalpha_0, float _s_dOut_0)
{
    float _S4 = (F32_max(((*dpalpha_0).primal_0), (1.00000001335143196e-10f)));
    float _S5 = _s_dOut_0 / (_S4 * _S4);
    float _S6 = (*dpdepth_0).primal_0 * - _S5;
    float _S7 = _S4 * _S5;
    DiffPair_float_0 _S8;
    (&_S8)->primal_0 = (*dpalpha_0).primal_0;
    (&_S8)->differential_0 = 0.0f;
    DiffPair_float_0 _S9;
    (&_S9)->primal_0 = 1.00000001335143196e-10f;
    (&_S9)->differential_0 = 0.0f;
    _d_max_0(&_S8, &_S9, _S6);
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S8.differential_0;
    dpdepth_0->primal_0 = (*dpdepth_0).primal_0;
    dpdepth_0->differential_0 = _S7;
    return;
}

inline __device__ void s_bwd_rendered_depth_to_expected_depth_0(DiffPair_float_0 * _S10, DiffPair_float_0 * _S11, float _S12)
{
    s_bwd_prop_rendered_depth_to_expected_depth_0(_S10, _S11, _S12);
    return;
}

inline __device__ void rendered_depth_to_expected_depth_bwd(float depth_1, float alpha_1, float v_out_depth_0, float * v_depth_0, float * v_alpha_0)
{
    DiffPair_float_0 p_depth_0;
    (&p_depth_0)->primal_0 = depth_1;
    (&p_depth_0)->differential_0 = 0.0f;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    s_bwd_rendered_depth_to_expected_depth_0(&p_depth_0, &p_alpha_0, v_out_depth_0);
    *v_depth_0 = p_depth_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
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
    DiffPair_float_0 _S13 = *dpx_1;
    bool _S14;
    if(((*dpx_1).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S14 = ((*dpx_1).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S14 = false;
    }
    float _S15;
    if(_S14)
    {
        _S15 = dOut_1;
    }
    else
    {
        _S15 = 0.0f;
    }
    dpx_1->primal_0 = _S13.primal_0;
    dpx_1->differential_0 = _S15;
    DiffPair_float_0 _S16 = *dpMin_0;
    if((_S13.primal_0) < ((*dpMin_0).primal_0))
    {
        _S15 = dOut_1;
    }
    else
    {
        _S15 = 0.0f;
    }
    dpMin_0->primal_0 = _S16.primal_0;
    dpMin_0->differential_0 = _S15;
    DiffPair_float_0 _S17 = *dpMax_0;
    if(((*dpx_1).primal_0) > ((*dpMax_0).primal_0))
    {
        _S15 = dOut_1;
    }
    else
    {
        _S15 = 0.0f;
    }
    dpMax_0->primal_0 = _S17.primal_0;
    dpMax_0->differential_0 = _S15;
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

inline __device__ float3  blend_background(float3  rgb_0, float alpha_2, float3  background_0)
{
    return clamp_0(rgb_0 + make_float3 (1.0f - alpha_2) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S18, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S19, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S20, float3  _S21)
{
    _d_clamp_vector_0(_S18, _S19, _S20, _S21);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_1)
{
    float _S22 = 1.0f - (*dpalpha_1).primal_0;
    float3  _S23 = make_float3 (_S22);
    float3  _S24 = make_float3 (0.0f);
    float3  _S25 = make_float3 (1.0f);
    float3  _S26 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S27;
    (&_S27)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S22) * (*dpbackground_0).primal_0;
    (&_S27)->differential_0 = _S26;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S28;
    (&_S28)->primal_0 = _S24;
    (&_S28)->differential_0 = _S26;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S29;
    (&_S29)->primal_0 = _S25;
    (&_S29)->differential_0 = _S26;
    s_bwd_prop_clamp_0(&_S27, &_S28, &_S29, _s_dOut_1);
    float3  _S30 = _S23 * _S27.differential_0;
    float3  _S31 = (*dpbackground_0).primal_0 * _S27.differential_0;
    float _S32 = - (_S31.x + _S31.y + _S31.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S30;
    dpalpha_1->primal_0 = (*dpalpha_1).primal_0;
    dpalpha_1->differential_0 = _S32;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S27.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S33, DiffPair_float_0 * _S34, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S35, float3  _S36)
{
    s_bwd_prop_blend_background_0(_S33, _S34, _S35, _S36);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_3, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_1, float3  * v_background_0)
{
    float3  _S37 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S37;
    DiffPair_float_0 p_alpha_1;
    (&p_alpha_1)->primal_0 = alpha_3;
    (&p_alpha_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S37;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_1, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_1 = p_alpha_1.differential_0;
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
        DiffPair_float_0 _S38 = *dpx_3;
        float _S39 = val_0 * (*dpy_2).primal_0 / (*dpx_3).primal_0 * dOut_3;
        dpx_3->primal_0 = (*dpx_3).primal_0;
        dpx_3->differential_0 = _S39;
        float _S40 = val_0 * (F32_log((_S38.primal_0))) * dOut_3;
        dpy_2->primal_0 = (*dpy_2).primal_0;
        dpy_2->differential_0 = _S40;
    }
    return;
}

inline __device__ DiffPair_float_0 _d_pow_1(DiffPair_float_0 * dpx_4, DiffPair_float_0 * dpy_3)
{
    float _S41 = dpx_4->primal_0;
    if((dpx_4->primal_0) < 9.99999997475242708e-07f)
    {
        DiffPair_float_0 _S42 = { 0.0f, 0.0f };
        return _S42;
    }
    float val_1 = (F32_pow((_S41), (dpy_3->primal_0)));
    DiffPair_float_0 _S43 = { val_1, val_1 * (F32_log((_S41))) * dpy_3->differential_0 + val_1 * dpy_3->primal_0 / _S41 * dpx_4->differential_0 };
    return _S43;
}

inline __device__ float linear_rgb_to_srgb(float x_3)
{
    float _S44;
    if(x_3 < 0.00313080009073019f)
    {
        _S44 = x_3 * 12.92000007629394531f;
    }
    else
    {
        _S44 = 1.0549999475479126f * (F32_pow((x_3), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return _S44;
}

inline __device__ float linear_rgb_to_srgb_grad(float x_4)
{
    float _S45;
    if(x_4 < 0.00313080009073019f)
    {
        _S45 = 12.92000007629394531f;
    }
    else
    {
        DiffPair_float_0 _S46;
        (&_S46)->primal_0 = x_4;
        (&_S46)->differential_0 = 1.0f;
        DiffPair_float_0 _S47;
        (&_S47)->primal_0 = 0.4166666567325592f;
        (&_S47)->differential_0 = 0.0f;
        DiffPair_float_0 _S48 = _d_pow_1(&_S46, &_S47);
        _S45 = _S48.differential_0 * 1.0549999475479126f;
    }
    return _S45;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float _S49 = rgb_2.x;
    float _S50;
    if(_S49 < 0.00313080009073019f)
    {
        _S50 = _S49 * 12.92000007629394531f;
    }
    else
    {
        _S50 = 1.0549999475479126f * (F32_pow((_S49), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S51 = rgb_2.y;
    float _S52;
    if(_S51 < 0.00313080009073019f)
    {
        _S52 = _S51 * 12.92000007629394531f;
    }
    else
    {
        _S52 = 1.0549999475479126f * (F32_pow((_S51), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S53 = rgb_2.z;
    float _S54;
    if(_S53 < 0.00313080009073019f)
    {
        _S54 = _S53 * 12.92000007629394531f;
    }
    else
    {
        _S54 = 1.0549999475479126f * (F32_pow((_S53), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S50, _S52, _S54);
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S55, DiffPair_float_0 * _S56, float _S57)
{
    _d_pow_0(_S55, _S56, _S57);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_2)
{
    float _S58 = (*dprgb_1).primal_0.x;
    float _S59 = (*dprgb_1).primal_0.y;
    float _S60 = (*dprgb_1).primal_0.z;
    float _S61;
    if(_S60 < 0.00313080009073019f)
    {
        _S61 = 12.92000007629394531f * _s_dOut_2.z;
    }
    else
    {
        float _S62 = 1.0549999475479126f * _s_dOut_2.z;
        DiffPair_float_0 _S63;
        (&_S63)->primal_0 = _S60;
        (&_S63)->differential_0 = 0.0f;
        DiffPair_float_0 _S64;
        (&_S64)->primal_0 = 0.4166666567325592f;
        (&_S64)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S63, &_S64, _S62);
        _S61 = _S63.differential_0;
    }
    float _S65;
    if(_S59 < 0.00313080009073019f)
    {
        _S65 = 12.92000007629394531f * _s_dOut_2.y;
    }
    else
    {
        float _S66 = 1.0549999475479126f * _s_dOut_2.y;
        DiffPair_float_0 _S67;
        (&_S67)->primal_0 = _S59;
        (&_S67)->differential_0 = 0.0f;
        DiffPair_float_0 _S68;
        (&_S68)->primal_0 = 0.4166666567325592f;
        (&_S68)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S67, &_S68, _S66);
        _S65 = _S67.differential_0;
    }
    float _S69;
    if(_S58 < 0.00313080009073019f)
    {
        _S69 = 12.92000007629394531f * _s_dOut_2.x;
    }
    else
    {
        float _S70 = 1.0549999475479126f * _s_dOut_2.x;
        DiffPair_float_0 _S71;
        (&_S71)->primal_0 = _S58;
        (&_S71)->differential_0 = 0.0f;
        DiffPair_float_0 _S72;
        (&_S72)->primal_0 = 0.4166666567325592f;
        (&_S72)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S71, &_S72, _S70);
        _S69 = _S71.differential_0;
    }
    float3  _S73 = make_float3 (_S69, _S65, _S61);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S73;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S74, float3  _S75)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S74, _S75);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S76 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S76;
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
    float _S77 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_5).primal_0)))))) * dOut_4;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S77;
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
        float _S78 = (*dist_coeffs_0)[int(3)];
        float _S79 = (*dist_coeffs_0)[int(4)];
        float _S80 = (*dist_coeffs_0)[int(5)];
        float _S81 = (*dist_coeffs_0)[int(6)];
        float _S82 = (*dist_coeffs_0)[int(7)];
        float _S83 = (*dist_coeffs_0)[int(8)];
        float _S84 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S85 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S86 = (*dist_coeffs_0)[int(1)] + r2_0 * _S85;
        float _S87 = (*dist_coeffs_0)[int(0)] + r2_0 * _S86;
        float radial_0 = 1.0f + r2_0 * _S87;
        float _S88 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S89 = _S88 * u_0;
        float _S90 = 2.0f * u_0;
        float _S91 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S92 = _S91 * u_0;
        float _S93 = 2.0f * v_0;
        float2  _S94 = q_0 * make_float2 (radial_0) + make_float2 (_S89 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S90 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S92 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S93 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S94 + make_float2 ((*dist_coeffs_0)[int(8)] * _S94.x + (*dist_coeffs_0)[int(9)] * _S94.y, 0.0f) - uv_0;
        float _S95 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S95 + _S95);
        float2  _S96 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S87 + (s_diff_r2_0 * _S86 + (s_diff_r2_0 * _S85 + s_diff_r2_0 * _S78 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S88 * v_0 + 0.0f * _S89 + (s_diff_r2_0 + (_S90 + _S90)) * _S80 + s_diff_r2_0 * _S81, _S91 * v_0 + 0.0f * _S92 + (s_diff_r2_0 + (_S95 + 0.0f * _S93)) * _S79 + s_diff_r2_0 * _S82);
        float _S97 = 0.0f * u_0;
        float s_diff_r2_1 = _S97 + _S97 + (v_0 + v_0);
        float2  _S98 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S87 + (s_diff_r2_1 * _S86 + (s_diff_r2_1 * _S85 + s_diff_r2_1 * _S78 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S88 * v_0 + _S89 + (s_diff_r2_1 + (_S97 + 0.0f * _S90)) * _S80 + s_diff_r2_1 * _S81, 0.0f * _S91 * v_0 + _S92 + (s_diff_r2_1 + (_S93 + _S93)) * _S79 + s_diff_r2_1 * _S82);
        Matrix<float, 2, 2>  _S99 = transpose_0(makeMatrix<float, 2, 2> (_S96 + make_float2 (_S96.x * _S83 + _S96.y * _S84, 0.0f), _S98 + make_float2 (_S98.x * _S83 + _S98.y * _S84, 0.0f)));
        float inv_det_0 = 1.0f / (_S99.rows[int(0)].x * _S99.rows[int(1)].y - _S99.rows[int(0)].y * _S99.rows[int(1)].x);
        float _S100 = r_1.x;
        float _S101 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S100 * _S99.rows[int(1)].y - _S101 * _S99.rows[int(0)].y) * inv_det_0, (- _S100 * _S99.rows[int(1)].x + _S101 * _S99.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S102 = (*dist_coeffs_0)[int(0)];
    float _S103 = (*dist_coeffs_0)[int(1)];
    float _S104 = (*dist_coeffs_0)[int(2)];
    float _S105 = (*dist_coeffs_0)[int(3)];
    float _S106 = (*dist_coeffs_0)[int(4)];
    float _S107 = (*dist_coeffs_0)[int(5)];
    float _S108 = (*dist_coeffs_0)[int(6)];
    float _S109 = (*dist_coeffs_0)[int(7)];
    float _S110 = (*dist_coeffs_0)[int(8)];
    float _S111 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S112 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S112 + _S112);
    float _S113 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S114 = (*dist_coeffs_0)[int(1)] + r2_1 * _S113;
    float _S115 = (*dist_coeffs_0)[int(0)] + r2_1 * _S114;
    float radial_1 = 1.0f + r2_1 * _S115;
    float _S116 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S117 = _S116 * u_1;
    float _S118 = 2.0f * u_1;
    float _S119 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S120 = _S119 * u_1;
    float _S121 = 2.0f * v_1;
    float2  _S122 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S115 + (s_diff_r2_2 * _S114 + (s_diff_r2_2 * _S113 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S116 * v_1 + 0.0f * _S117 + (s_diff_r2_2 + (_S118 + _S118)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S119 * v_1 + 0.0f * _S120 + (s_diff_r2_2 + (_S112 + 0.0f * _S121)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S123 = 0.0f * u_1;
    float s_diff_r2_3 = _S123 + _S123 + (v_1 + v_1);
    float2  _S124 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S115 + (s_diff_r2_3 * _S114 + (s_diff_r2_3 * _S113 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S116 * v_1 + _S117 + (s_diff_r2_3 + (_S123 + 0.0f * _S118)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S119 * v_1 + _S120 + (s_diff_r2_3 + (_S121 + _S121)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S125 = transpose_0(makeMatrix<float, 2, 2> (_S122 + make_float2 (_S122.x * (*dist_coeffs_0)[int(8)] + _S122.y * (*dist_coeffs_0)[int(9)], 0.0f), _S124 + make_float2 (_S124.x * (*dist_coeffs_0)[int(8)] + _S124.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S126;
    if((F32_min((determinant_0(_S125)), ((F32_min((_S125.rows[int(0)].x), (_S125.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S127 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S102 + r2_2 * (_S103 + r2_2 * (_S104 + r2_2 * _S105)))) + make_float2 (_S116 * u_2 * v_2 + _S107 * (r2_2 + 2.0f * u_2 * u_2) + _S108 * r2_2, _S119 * u_2 * v_2 + _S106 * (r2_2 + 2.0f * v_2 * v_2) + _S109 * r2_2);
        _S126 = (length_0(_S127 + make_float2 (_S110 * _S127.x + _S111 * _S127.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S126 = false;
    }
    return _S126;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S128 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S128;
    FixedArray<float, 10>  _S129 = dist_coeffs_1;
    bool _S130 = undistort_point_0(_S128, &_S129, int(12), &uv_1);
    if(!_S130)
    {
        int3  _S131 = make_int3 (int(0));
        float3  _S132 = make_float3 ((float)_S131.x, (float)_S131.y, (float)_S131.z);
        return _S132;
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
    float _S133 = dOut_6.y;
    float _S134 = dOut_6.z;
    float _S135 = dOut_6.x;
    float _S136 = (*a_0).primal_0.z * _S133 + - (*a_0).primal_0.y * _S134;
    float _S137 = - (*a_0).primal_0.z * _S135 + (*a_0).primal_0.x * _S134;
    float _S138 = (*a_0).primal_0.y * _S135 + - (*a_0).primal_0.x * _S133;
    float3  _S139 = make_float3 (- (*b_0).primal_0.z * _S133 + (*b_0).primal_0.y * _S134, (*b_0).primal_0.z * _S135 + - (*b_0).primal_0.x * _S134, - (*b_0).primal_0.y * _S135 + (*b_0).primal_0.x * _S133);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S139;
    float3  _S140 = make_float3 (_S136, _S137, _S138);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S140;
    return;
}

inline __device__ float3  cross_0(float3  left_0, float3  right_0)
{
    float _S141 = left_0.y;
    float _S142 = right_0.z;
    float _S143 = left_0.z;
    float _S144 = right_0.y;
    float _S145 = right_0.x;
    float _S146 = left_0.x;
    return make_float3 (_S141 * _S142 - _S143 * _S144, _S143 * _S145 - _S146 * _S142, _S146 * _S144 - _S141 * _S145);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S147, float3  _S148)
{
    return cross_0(_S147, _S148);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S149, float3  _S150)
{
    return dot_0(_S149, _S150);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S151, float _S152)
{
    _d_sqrt_0(_S151, _S152);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_3)
{
    float _S153 = (*dpx_7).primal_0.x;
    float _S154 = (*dpx_7).primal_0.y;
    float _S155 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S156;
    (&_S156)->primal_0 = _S153 * _S153 + _S154 * _S154 + _S155 * _S155;
    (&_S156)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S156, _s_dOut_3);
    float _S157 = (*dpx_7).primal_0.z * _S156.differential_0;
    float _S158 = _S157 + _S157;
    float _S159 = (*dpx_7).primal_0.y * _S156.differential_0;
    float _S160 = _S159 + _S159;
    float _S161 = (*dpx_7).primal_0.x * _S156.differential_0;
    float _S162 = _S161 + _S161;
    float3  _S163 = make_float3 (0.0f);
    *&((&_S163)->z) = _S158;
    *&((&_S163)->y) = _S160;
    *&((&_S163)->x) = _S162;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S163;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S164, float _S165)
{
    s_bwd_prop_length_impl_0(_S164, _S165);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S166, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S167, float _S168)
{
    _d_dot_0(_S166, _S167, _S168);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S169, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S170, float3  _S171)
{
    _d_cross_0(_S169, _S170, _S171);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_4)
{
    float3  _S172 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S173 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S174 = s_primal_ctx_cross_0(dx_0, _S173);
    bool _S175 = (s_primal_ctx_dot_0(_S174, _S174)) != 0.0f;
    float3  _S176;
    float3  _S177;
    if(_S175)
    {
        float _S178 = length_1(_S174);
        float3  _S179 = make_float3 (_S178);
        _S176 = make_float3 (_S178 * _S178);
        _S177 = _S179;
    }
    else
    {
        _S176 = _S172;
        _S177 = _S172;
    }
    if(_S175)
    {
        float3  _S180 = _s_dOut_4 / _S176;
        float3  _S181 = _S174 * - _S180;
        float3  _S182 = _S177 * _S180;
        float _S183 = _S181.x + _S181.y + _S181.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S184;
        (&_S184)->primal_0 = _S174;
        (&_S184)->differential_0 = _S172;
        s_bwd_length_impl_0(&_S184, _S183);
        _S176 = _S182 + _S184.differential_0;
    }
    else
    {
        _S176 = _s_dOut_4;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S185;
    (&_S185)->primal_0 = _S174;
    (&_S185)->differential_0 = _S172;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S186;
    (&_S186)->primal_0 = _S174;
    (&_S186)->differential_0 = _S172;
    s_bwd_prop_dot_0(&_S185, &_S186, 0.0f);
    float3  _S187 = _S186.differential_0 + _S185.differential_0 + _S176;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S188;
    (&_S188)->primal_0 = dx_0;
    (&_S188)->differential_0 = _S172;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S189;
    (&_S189)->primal_0 = _S173;
    (&_S189)->differential_0 = _S172;
    s_bwd_prop_cross_0(&_S188, &_S189, _S187);
    float3  s_diff_dy_T_0 = - _S189.differential_0;
    float3  _S190 = - s_diff_dy_T_0;
    float3  _S191 = - _S188.differential_0;
    FixedArray<float3 , 4>  _S192;
    _S192[int(0)] = _S172;
    _S192[int(1)] = _S172;
    _S192[int(2)] = _S172;
    _S192[int(3)] = _S172;
    _S192[int(2)] = _S190;
    _S192[int(3)] = s_diff_dy_T_0;
    _S192[int(0)] = _S191;
    _S192[int(1)] = _S188.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S192;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S193, float3  _S194)
{
    s_bwd_prop_points_to_normal_0(_S193, _S194);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S195 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S195;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S196 = float2 {intrins_1.z, intrins_1.w};
    float2  _S197 = float2 {intrins_1.x, intrins_1.y};
    float2  _S198 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S196) / _S197;
    float2  uv_2 = _S198;
    FixedArray<float, 10>  _S199 = dist_coeffs_2;
    bool _S200 = undistort_point_0(_S198, &_S199, int(12), &uv_2);
    if(!_S200)
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
    float2  _S201 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S196) / _S197;
    float2  uv_3 = _S201;
    FixedArray<float, 10>  _S202 = dist_coeffs_2;
    bool _S203 = undistort_point_0(_S201, &_S202, int(12), &uv_3);
    if(!_S203)
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
    float2  _S204 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S196) / _S197;
    float2  uv_4 = _S204;
    FixedArray<float, 10>  _S205 = dist_coeffs_2;
    bool _S206 = undistort_point_0(_S204, &_S205, int(12), &uv_4);
    if(!_S206)
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
    float2  _S207 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S196) / _S197;
    float2  uv_5 = _S207;
    FixedArray<float, 10>  _S208 = dist_coeffs_2;
    bool _S209 = undistort_point_0(_S207, &_S208, int(12), &uv_5);
    if(!_S209)
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
    float2  _S210;
    bool _S211;
    float2  _S212;
    bool _S213;
    float2  _S214;
    bool _S215;
    float2  _S216;
    bool _S217;
};

inline __device__ float s_primal_ctx_sin_0(float _S218)
{
    return (F32_sin((_S218)));
}

inline __device__ float s_primal_ctx_cos_0(float _S219)
{
    return (F32_cos((_S219)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S220 = make_float2 (0.0f);
    _s_diff_ctx_0->_S210 = _S220;
    _s_diff_ctx_0->_S211 = false;
    _s_diff_ctx_0->_S212 = _S220;
    _s_diff_ctx_0->_S213 = false;
    _s_diff_ctx_0->_S214 = _S220;
    _s_diff_ctx_0->_S215 = false;
    _s_diff_ctx_0->_S216 = _S220;
    _s_diff_ctx_0->_S217 = false;
    _s_diff_ctx_0->_S212 = _S220;
    _s_diff_ctx_0->_S213 = false;
    _s_diff_ctx_0->_S214 = _S220;
    _s_diff_ctx_0->_S215 = false;
    _s_diff_ctx_0->_S216 = _S220;
    _s_diff_ctx_0->_S217 = false;
    float3  _S221 = make_float3 (0.0f);
    float2  _S222 = float2 {intrins_2.z, intrins_2.w};
    float2  _S223 = float2 {intrins_2.x, intrins_2.y};
    float2  _S224 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S222) / _S223;
    float2  _S225 = _S224;
    bool _S226 = undistort_point_0(_S224, dist_coeffs_3, int(12), &_S225);
    _s_diff_ctx_0->_S210 = _S225;
    _s_diff_ctx_0->_S211 = _S226;
    float2  uv_6 = _S225;
    bool _S227 = !_S226;
    float3  normal_4;
    if(_S227)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S228 = !_S227;
    int _S229;
    FixedArray<float3 , 4>  points_3;
    if(_S228)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S230 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S230), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S230))).x, (uv_6 / make_float2 ((F32_max((_S230), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S230))).y, s_primal_ctx_cos_0(_S230));
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
        float3  _S231 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S232 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S222) / _S223;
        float2  _S233 = _S232;
        bool _S234 = undistort_point_0(_S232, dist_coeffs_3, int(12), &_S233);
        _s_diff_ctx_0->_S212 = _S233;
        _s_diff_ctx_0->_S213 = _S234;
        float2  uv_7 = _S233;
        bool _S235 = !_S234;
        if(_S235)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S236 = !_S235;
        if(_S236)
        {
            if(is_fisheye_2)
            {
                float _S237 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S237), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S237))).x, (uv_7 / make_float2 ((F32_max((_S237), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S237))).y, s_primal_ctx_cos_0(_S237));
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
            float3  _S238 = make_float3 (dpdepths_0.y) * raydir_12;
            _S229 = int(2);
            points_3[int(0)] = _S231;
            points_3[int(1)] = _S238;
            points_3[int(2)] = _S221;
            points_3[int(3)] = _S221;
        }
        else
        {
            _S229 = int(0);
            points_3[int(0)] = _S231;
            points_3[int(1)] = _S221;
            points_3[int(2)] = _S221;
            points_3[int(3)] = _S221;
        }
        bool _runFlag_0;
        if(_S229 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S228;
            _S229 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S239 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S222) / _S223;
            float2  _S240 = _S239;
            bool _S241 = undistort_point_0(_S239, dist_coeffs_3, int(12), &_S240);
            _s_diff_ctx_0->_S214 = _S240;
            _s_diff_ctx_0->_S215 = _S241;
            float2  uv_8 = _S240;
            if(!_S241)
            {
                float3  _S242 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S229 = int(0);
                normal_4 = _S242;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S243 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S243), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S243))).x, (uv_8 / make_float2 ((F32_max((_S243), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S243))).y, s_primal_ctx_cos_0(_S243));
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
                float2  _S244 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S222) / _S223;
                float2  _S245 = _S244;
                bool _S246 = undistort_point_0(_S244, dist_coeffs_3, int(12), &_S245);
                _s_diff_ctx_0->_S216 = _S245;
                _s_diff_ctx_0->_S217 = _S246;
                float2  uv_9 = _S245;
                bool _S247 = !_S246;
                if(_S247)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S248 = !_S247;
                int _S249;
                if(_S248)
                {
                    if(is_fisheye_2)
                    {
                        float _S250 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S250), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S250))).x, (uv_9 / make_float2 ((F32_max((_S250), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S250))).y, s_primal_ctx_cos_0(_S250));
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
                    _S249 = int(2);
                }
                else
                {
                    _S249 = int(0);
                }
                if(_S249 != int(2))
                {
                    _runFlag_0 = false;
                    _S229 = _S249;
                }
                if(_runFlag_0)
                {
                    _S229 = int(1);
                }
            }
        }
    }
    else
    {
        _S229 = int(0);
        points_3[int(0)] = _S221;
        points_3[int(1)] = _S221;
        points_3[int(2)] = _S221;
        points_3[int(3)] = _S221;
    }
    if(!(_S229 != int(1)))
    {
        float3  _S251 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S251, _S251)) != 0.0f)
        {
            normal_4 = _S251 / make_float3 (length_1(_S251));
        }
        else
        {
            normal_4 = _S251;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_5, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S252 = *dpdepths_1;
    float3  _S253 = make_float3 (0.0f);
    float2  _S254 = _s_diff_ctx_1->_S210;
    bool _S255 = !!_s_diff_ctx_1->_S211;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S256;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S257;
    if(_S255)
    {
        if(is_fisheye_3)
        {
            float _S258 = length_0(_S254);
            float3  raydir_25 = make_float3 ((_S254 / make_float2 ((F32_max((_S258), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S258))).x, (_S254 / make_float2 ((F32_max((_S258), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S258))).y, s_primal_ctx_cos_0(_S258));
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
            float3  raydir_26 = make_float3 (_S254.x, _S254.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S259 = make_float3 (_S252.primal_0.x) * raydir_21;
        float2  _S260 = _s_diff_ctx_1->_S212;
        bool _S261 = !!_s_diff_ctx_1->_S213;
        if(_S261)
        {
            if(is_fisheye_3)
            {
                float _S262 = length_0(_S260);
                float3  raydir_27 = make_float3 ((_S260 / make_float2 ((F32_max((_S262), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S262))).x, (_S260 / make_float2 ((F32_max((_S262), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S262))).y, s_primal_ctx_cos_0(_S262));
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
                float3  raydir_28 = make_float3 (_S260.x, _S260.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S263 = make_float3 (_S252.primal_0.y) * raydir_22;
            _S256 = int(2);
            points_4[int(0)] = _S259;
            points_4[int(1)] = _S263;
            points_4[int(2)] = _S253;
            points_4[int(3)] = _S253;
        }
        else
        {
            _S256 = int(0);
            points_4[int(0)] = _S259;
            points_4[int(1)] = _S253;
            points_4[int(2)] = _S253;
            points_4[int(3)] = _S253;
            raydir_22 = _S253;
        }
        if(_S256 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S255;
            _S256 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S264 = _s_diff_ctx_1->_S214;
            if(!_s_diff_ctx_1->_S215)
            {
                _runFlag_2 = false;
                _S256 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S265 = length_0(_S264);
                    float3  raydir_29 = make_float3 ((_S264 / make_float2 ((F32_max((_S265), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S265))).x, (_S264 / make_float2 ((F32_max((_S265), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S265))).y, s_primal_ctx_cos_0(_S265));
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
                    float3  raydir_30 = make_float3 (_S264.x, _S264.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S252.primal_0.z) * raydir_23;
                float2  _S266 = _s_diff_ctx_1->_S216;
                bool _S267 = !!_s_diff_ctx_1->_S217;
                int _S268;
                if(_S267)
                {
                    if(is_fisheye_3)
                    {
                        float _S269 = length_0(_S266);
                        float3  raydir_31 = make_float3 ((_S266 / make_float2 ((F32_max((_S269), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S269))).x, (_S266 / make_float2 ((F32_max((_S269), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S269))).y, s_primal_ctx_cos_0(_S269));
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
                        float3  raydir_32 = make_float3 (_S266.x, _S266.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S252.primal_0.w) * raydir_24;
                    _S268 = int(2);
                }
                else
                {
                    _S268 = int(0);
                    raydir_24 = _S253;
                }
                if(_S268 != int(2))
                {
                    _runFlag_3 = false;
                    _S256 = _S268;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S256 = int(1);
                }
                float3  _S270 = raydir_23;
                _runFlag_3 = _S267;
                raydir_23 = raydir_24;
                raydir_24 = _S270;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S253;
                raydir_24 = _S253;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S253;
            raydir_24 = _S253;
        }
        float3  _S271 = raydir_21;
        float3  _S272 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S257 = _S261;
        raydir_23 = _S272;
        raydir_24 = _S271;
    }
    else
    {
        _S256 = int(0);
        points_4[int(0)] = _S253;
        points_4[int(1)] = _S253;
        points_4[int(2)] = _S253;
        points_4[int(3)] = _S253;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S253;
        raydir_22 = _S253;
        _S257 = false;
        raydir_23 = _S253;
        raydir_24 = _S253;
    }
    bool _S273 = !(_S256 != int(1));
    float3  _S274;
    float3  _S275;
    float3  _S276;
    float3  _S277;
    float3  _S278;
    bool _S279;
    if(_S273)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S280 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S281 = s_primal_ctx_cross_0(dx_1, _S280);
        bool _S282 = (s_primal_ctx_dot_0(_S281, _S281)) != 0.0f;
        if(_S282)
        {
            float _S283 = length_1(_S281);
            float3  _S284 = make_float3 (_S283);
            _S274 = make_float3 (_S283 * _S283);
            _S275 = _S284;
        }
        else
        {
            _S274 = _S253;
            _S275 = _S253;
        }
        float3  _S285 = _S275;
        _S279 = _S282;
        _S275 = _S281;
        _S276 = _S285;
        _S277 = dx_1;
        _S278 = _S280;
    }
    else
    {
        _S279 = false;
        _S274 = _S253;
        _S275 = _S253;
        _S276 = _S253;
        _S277 = _S253;
        _S278 = _S253;
    }
    float4  _S286 = make_float4 (0.0f);
    if(_S273)
    {
        if(_S279)
        {
            float3  _S287 = _s_dOut_5 / _S274;
            float3  _S288 = _S275 * - _S287;
            float3  _S289 = _S276 * _S287;
            float _S290 = _S288.x + _S288.y + _S288.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S291;
            (&_S291)->primal_0 = _S275;
            (&_S291)->differential_0 = _S253;
            s_bwd_length_impl_0(&_S291, _S290);
            _S274 = _S289 + _S291.differential_0;
        }
        else
        {
            _S274 = _s_dOut_5;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S292;
        (&_S292)->primal_0 = _S275;
        (&_S292)->differential_0 = _S253;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S293;
        (&_S293)->primal_0 = _S275;
        (&_S293)->differential_0 = _S253;
        s_bwd_prop_dot_0(&_S292, &_S293, 0.0f);
        float3  _S294 = _S293.differential_0 + _S292.differential_0 + _S274;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S295;
        (&_S295)->primal_0 = _S277;
        (&_S295)->differential_0 = _S253;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S296;
        (&_S296)->primal_0 = _S278;
        (&_S296)->differential_0 = _S253;
        s_bwd_prop_cross_0(&_S295, &_S296, _S294);
        float3  s_diff_dy_T_1 = - _S296.differential_0;
        float3  _S297 = - s_diff_dy_T_1;
        float3  _S298 = - _S295.differential_0;
        FixedArray<float3 , 4>  _S299;
        _S299[int(0)] = _S253;
        _S299[int(1)] = _S253;
        _S299[int(2)] = _S253;
        _S299[int(3)] = _S253;
        _S299[int(2)] = _S297;
        _S299[int(3)] = s_diff_dy_T_1;
        _S299[int(0)] = _S298;
        _S299[int(1)] = _S295.differential_0;
        points_4[int(0)] = _S299[int(0)];
        points_4[int(1)] = _S299[int(1)];
        points_4[int(2)] = _S299[int(2)];
        points_4[int(3)] = _S299[int(3)];
    }
    else
    {
        points_4[int(0)] = _S253;
        points_4[int(1)] = _S253;
        points_4[int(2)] = _S253;
        points_4[int(3)] = _S253;
    }
    float4  _S300;
    if(_S255)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S301 = points_4;
                FixedArray<float3 , 4>  _S302 = points_4;
                FixedArray<float3 , 4>  _S303 = points_4;
                FixedArray<float3 , 4>  _S304 = points_4;
                if(_runFlag_3)
                {
                    float3  _S305 = raydir_21 * _S304[int(3)];
                    float _S306 = _S305.x + _S305.y + _S305.z;
                    float4  _S307 = _S286;
                    *&((&_S307)->w) = _S306;
                    points_4[int(0)] = _S301[int(0)];
                    points_4[int(1)] = _S302[int(1)];
                    points_4[int(2)] = _S303[int(2)];
                    points_4[int(3)] = _S253;
                    _S300 = _S307;
                }
                else
                {
                    points_4[int(0)] = _S301[int(0)];
                    points_4[int(1)] = _S302[int(1)];
                    points_4[int(2)] = _S303[int(2)];
                    points_4[int(3)] = _S304[int(3)];
                    _S300 = _S286;
                }
                float3  _S308 = raydir_22 * points_4[int(2)];
                float _S309 = _S308.x + _S308.y + _S308.z;
                FixedArray<float3 , 4>  _S310 = points_4;
                FixedArray<float3 , 4>  _S311 = points_4;
                float4  _S312 = _S286;
                *&((&_S312)->z) = _S309;
                float4  _S313 = _S300 + _S312;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S310[int(1)];
                points_4[int(2)] = _S253;
                points_4[int(3)] = _S311[int(3)];
                _S300 = _S313;
            }
            else
            {
                FixedArray<float3 , 4>  _S314 = points_4;
                FixedArray<float3 , 4>  _S315 = points_4;
                FixedArray<float3 , 4>  _S316 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S314[int(1)];
                points_4[int(2)] = _S315[int(2)];
                points_4[int(3)] = _S316[int(3)];
                _S300 = _S286;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S317 = points_4;
            FixedArray<float3 , 4>  _S318 = points_4;
            FixedArray<float3 , 4>  _S319 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S317[int(1)];
            points_4[int(2)] = _S318[int(2)];
            points_4[int(3)] = _S319[int(3)];
            _S300 = _S286;
        }
        if(_S257)
        {
            FixedArray<float3 , 4>  _S320 = points_4;
            float3  _S321 = raydir_23 * points_4[int(1)];
            float _S322 = _S321.x + _S321.y + _S321.z;
            float4  _S323 = _S286;
            *&((&_S323)->y) = _S322;
            float4  _S324 = _S300 + _S323;
            points_4[int(0)] = _S253;
            points_4[int(1)] = _S253;
            points_4[int(2)] = _S253;
            points_4[int(3)] = _S253;
            raydir_21 = _S320[int(0)];
            _S300 = _S324;
        }
        else
        {
            FixedArray<float3 , 4>  _S325 = points_4;
            FixedArray<float3 , 4>  _S326 = points_4;
            FixedArray<float3 , 4>  _S327 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S325[int(1)];
            points_4[int(2)] = _S326[int(2)];
            points_4[int(3)] = _S327[int(3)];
            raydir_21 = _S253;
        }
        float3  _S328 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S329 = _S328.x + _S328.y + _S328.z;
        float4  _S330 = _S286;
        *&((&_S330)->x) = _S329;
        _S300 = _S300 + _S330;
    }
    else
    {
        _S300 = _S286;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S300;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S331, float4  _S332, FixedArray<float, 10>  * _S333, bool _S334, bool _S335, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S336, float3  _S337)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S338;
    float3  _S339 = s_primal_ctx_depth_to_normal_0(_S331, _S332, _S333, _S334, _S335, (*_S336).primal_0, &_S338);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S340 = _S338;
    s_bwd_prop_depth_to_normal_0(_S331, _S332, _S333, _S334, _S335, _S336, _S337, &_S340);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S341 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S341;
    FixedArray<float, 10>  _S342 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S342, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S343 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S343;
    FixedArray<float, 10>  _S344 = dist_coeffs_6;
    bool _S345 = undistort_point_0(_S343, &_S344, int(12), &uv_10);
    if(!_S345)
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

