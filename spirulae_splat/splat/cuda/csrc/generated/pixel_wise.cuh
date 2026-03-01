#pragma once

#include "slang.cuh"

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

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_0)
{
    DiffPair_float_0 _S1 = *dpx_0;
    bool _S2;
    if(((*dpx_0).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S2 = ((*dpx_0).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S2 = false;
    }
    float _S3;
    if(_S2)
    {
        _S3 = dOut_0;
    }
    else
    {
        _S3 = 0.0f;
    }
    dpx_0->primal_0 = _S1.primal_0;
    dpx_0->differential_0 = _S3;
    DiffPair_float_0 _S4 = *dpMin_0;
    if((_S1.primal_0) < ((*dpMin_0).primal_0))
    {
        _S3 = dOut_0;
    }
    else
    {
        _S3 = 0.0f;
    }
    dpMin_0->primal_0 = _S4.primal_0;
    dpMin_0->differential_0 = _S3;
    DiffPair_float_0 _S5 = *dpMax_0;
    if(((*dpx_0).primal_0) > ((*dpMax_0).primal_0))
    {
        _S3 = dOut_0;
    }
    else
    {
        _S3 = 0.0f;
    }
    dpMax_0->primal_0 = _S5.primal_0;
    dpMax_0->differential_0 = _S3;
    return;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_1)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_1).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_0).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_1.x);
    float3  left_d_result_0;
    *&((&left_d_result_0)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_1).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_0).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_1.y);
    *&((&left_d_result_0)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_0)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_1).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_0).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_1.z);
    *&((&left_d_result_0)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_0)->z) = right_dp_2.differential_0;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = left_d_result_0;
    dpy_0->primal_0 = (*dpy_0).primal_0;
    dpy_0->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  clamp_0(float3  x_2, float3  minBound_0, float3  maxBound_0)
{
    return min_0(max_0(x_2, minBound_0), maxBound_0);
}

inline __device__ float3  blend_background(float3  rgb_0, float alpha_0, float3  background_0)
{
    return clamp_0(rgb_0 + make_float3 (1.0f - alpha_0) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S6, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S8, float3  _S9)
{
    _d_clamp_vector_0(_S6, _S7, _S8, _S9);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dpalpha_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_0)
{
    float _S10 = 1.0f - (*dpalpha_0).primal_0;
    float3  _S11 = make_float3 (_S10);
    float3  _S12 = make_float3 (0.0f);
    float3  _S13 = make_float3 (1.0f);
    float3  _S14 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S15;
    (&_S15)->primal_0 = (*dprgb_0).primal_0 + make_float3 (_S10) * (*dpbackground_0).primal_0;
    (&_S15)->differential_0 = _S14;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S16;
    (&_S16)->primal_0 = _S12;
    (&_S16)->differential_0 = _S14;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S17;
    (&_S17)->primal_0 = _S13;
    (&_S17)->differential_0 = _S14;
    s_bwd_prop_clamp_0(&_S15, &_S16, &_S17, _s_dOut_0);
    float3  _S18 = _S11 * _S15.differential_0;
    float3  _S19 = (*dpbackground_0).primal_0 * _S15.differential_0;
    float _S20 = - (_S19.x + _S19.y + _S19.z);
    dpbackground_0->primal_0 = (*dpbackground_0).primal_0;
    dpbackground_0->differential_0 = _S18;
    dpalpha_0->primal_0 = (*dpalpha_0).primal_0;
    dpalpha_0->differential_0 = _S20;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = _S15.differential_0;
    return;
}

inline __device__ void s_bwd_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S21, DiffPair_float_0 * _S22, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S23, float3  _S24)
{
    s_bwd_prop_blend_background_0(_S21, _S22, _S23, _S24);
    return;
}

inline __device__ void blend_background_bwd(float3  rgb_1, float alpha_1, float3  background_1, float3  v_out_rgb_0, float3  * v_rgb_0, float * v_alpha_0, float3  * v_background_0)
{
    float3  _S25 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_0;
    (&p_rgb_0)->primal_0 = rgb_1;
    (&p_rgb_0)->differential_0 = _S25;
    DiffPair_float_0 p_alpha_0;
    (&p_alpha_0)->primal_0 = alpha_1;
    (&p_alpha_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_background_0;
    (&p_background_0)->primal_0 = background_1;
    (&p_background_0)->differential_0 = _S25;
    s_bwd_blend_background_0(&p_rgb_0, &p_alpha_0, &p_background_0, v_out_rgb_0);
    *v_rgb_0 = p_rgb_0.differential_0;
    *v_alpha_0 = p_alpha_0.differential_0;
    *v_background_0 = p_background_0.differential_0;
    return;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpy_1, float dOut_2)
{
    if(((*dpx_2).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_2->primal_0 = (*dpx_2).primal_0;
        dpx_2->differential_0 = 0.0f;
        dpy_1->primal_0 = (*dpy_1).primal_0;
        dpy_1->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_2).primal_0), ((*dpy_1).primal_0)));
        DiffPair_float_0 _S26 = *dpx_2;
        float _S27 = val_0 * (*dpy_1).primal_0 / (*dpx_2).primal_0 * dOut_2;
        dpx_2->primal_0 = (*dpx_2).primal_0;
        dpx_2->differential_0 = _S27;
        float _S28 = val_0 * (F32_log((_S26.primal_0))) * dOut_2;
        dpy_1->primal_0 = (*dpy_1).primal_0;
        dpy_1->differential_0 = _S28;
    }
    return;
}

inline __device__ DiffPair_float_0 _d_pow_1(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_2)
{
    float _S29 = dpx_3->primal_0;
    if((dpx_3->primal_0) < 9.99999997475242708e-07f)
    {
        DiffPair_float_0 _S30 = { 0.0f, 0.0f };
        return _S30;
    }
    float val_1 = (F32_pow((_S29), (dpy_2->primal_0)));
    DiffPair_float_0 _S31 = { val_1, val_1 * (F32_log((_S29))) * dpy_2->differential_0 + val_1 * dpy_2->primal_0 / _S29 * dpx_3->differential_0 };
    return _S31;
}

inline __device__ float linear_rgb_to_srgb(float x_3)
{
    float _S32;
    if(x_3 < 0.00313080009073019f)
    {
        _S32 = x_3 * 12.92000007629394531f;
    }
    else
    {
        _S32 = 1.0549999475479126f * (F32_pow((x_3), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return _S32;
}

inline __device__ float linear_rgb_to_srgb_grad(float x_4)
{
    float _S33;
    if(x_4 < 0.00313080009073019f)
    {
        _S33 = 12.92000007629394531f;
    }
    else
    {
        DiffPair_float_0 _S34;
        (&_S34)->primal_0 = x_4;
        (&_S34)->differential_0 = 1.0f;
        DiffPair_float_0 _S35;
        (&_S35)->primal_0 = 0.4166666567325592f;
        (&_S35)->differential_0 = 0.0f;
        DiffPair_float_0 _S36 = _d_pow_1(&_S34, &_S35);
        _S33 = _S36.differential_0 * 1.0549999475479126f;
    }
    return _S33;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float _S37 = rgb_2.x;
    float _S38;
    if(_S37 < 0.00313080009073019f)
    {
        _S38 = _S37 * 12.92000007629394531f;
    }
    else
    {
        _S38 = 1.0549999475479126f * (F32_pow((_S37), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S39 = rgb_2.y;
    float _S40;
    if(_S39 < 0.00313080009073019f)
    {
        _S40 = _S39 * 12.92000007629394531f;
    }
    else
    {
        _S40 = 1.0549999475479126f * (F32_pow((_S39), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S41 = rgb_2.z;
    float _S42;
    if(_S41 < 0.00313080009073019f)
    {
        _S42 = _S41 * 12.92000007629394531f;
    }
    else
    {
        _S42 = 1.0549999475479126f * (F32_pow((_S41), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S38, _S40, _S42);
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S43, DiffPair_float_0 * _S44, float _S45)
{
    _d_pow_0(_S43, _S44, _S45);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_1)
{
    float _S46 = (*dprgb_1).primal_0.x;
    float _S47 = (*dprgb_1).primal_0.y;
    float _S48 = (*dprgb_1).primal_0.z;
    float _S49;
    if(_S48 < 0.00313080009073019f)
    {
        _S49 = 12.92000007629394531f * _s_dOut_1.z;
    }
    else
    {
        float _S50 = 1.0549999475479126f * _s_dOut_1.z;
        DiffPair_float_0 _S51;
        (&_S51)->primal_0 = _S48;
        (&_S51)->differential_0 = 0.0f;
        DiffPair_float_0 _S52;
        (&_S52)->primal_0 = 0.4166666567325592f;
        (&_S52)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S51, &_S52, _S50);
        _S49 = _S51.differential_0;
    }
    float _S53;
    if(_S47 < 0.00313080009073019f)
    {
        _S53 = 12.92000007629394531f * _s_dOut_1.y;
    }
    else
    {
        float _S54 = 1.0549999475479126f * _s_dOut_1.y;
        DiffPair_float_0 _S55;
        (&_S55)->primal_0 = _S47;
        (&_S55)->differential_0 = 0.0f;
        DiffPair_float_0 _S56;
        (&_S56)->primal_0 = 0.4166666567325592f;
        (&_S56)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S55, &_S56, _S54);
        _S53 = _S55.differential_0;
    }
    float _S57;
    if(_S46 < 0.00313080009073019f)
    {
        _S57 = 12.92000007629394531f * _s_dOut_1.x;
    }
    else
    {
        float _S58 = 1.0549999475479126f * _s_dOut_1.x;
        DiffPair_float_0 _S59;
        (&_S59)->primal_0 = _S46;
        (&_S59)->differential_0 = 0.0f;
        DiffPair_float_0 _S60;
        (&_S60)->primal_0 = 0.4166666567325592f;
        (&_S60)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S59, &_S60, _S58);
        _S57 = _S59.differential_0;
    }
    float3  _S61 = make_float3 (_S57, _S53, _S49);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S61;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S62, float3  _S63)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S62, _S63);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S64 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S64;
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

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_4, float dOut_3)
{
    float _S65 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_4).primal_0)))))) * dOut_3;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S65;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_3, float dOut_4)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_3).primal_0.x * dOut_4;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_5).primal_0.x * dOut_4;
    *&((&x_d_result_0)->y) = (*dpy_3).primal_0.y * dOut_4;
    *&((&y_d_result_0)->y) = (*dpx_5).primal_0.y * dOut_4;
    *&((&x_d_result_0)->z) = (*dpy_3).primal_0.z * dOut_4;
    *&((&y_d_result_0)->z) = (*dpx_5).primal_0.z * dOut_4;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = x_d_result_0;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = y_d_result_0;
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
        float _S66 = (*dist_coeffs_0)[int(3)];
        float _S67 = (*dist_coeffs_0)[int(4)];
        float _S68 = (*dist_coeffs_0)[int(5)];
        float _S69 = (*dist_coeffs_0)[int(6)];
        float _S70 = (*dist_coeffs_0)[int(7)];
        float _S71 = (*dist_coeffs_0)[int(8)];
        float _S72 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S73 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S74 = (*dist_coeffs_0)[int(1)] + r2_0 * _S73;
        float _S75 = (*dist_coeffs_0)[int(0)] + r2_0 * _S74;
        float radial_0 = 1.0f + r2_0 * _S75;
        float _S76 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S77 = _S76 * u_0;
        float _S78 = 2.0f * u_0;
        float _S79 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S80 = _S79 * u_0;
        float _S81 = 2.0f * v_0;
        float2  _S82 = q_0 * make_float2 (radial_0) + make_float2 (_S77 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S78 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S80 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S81 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S82 + make_float2 ((*dist_coeffs_0)[int(8)] * _S82.x + (*dist_coeffs_0)[int(9)] * _S82.y, 0.0f) - uv_0;
        float _S83 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S83 + _S83);
        float2  _S84 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S75 + (s_diff_r2_0 * _S74 + (s_diff_r2_0 * _S73 + s_diff_r2_0 * _S66 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S76 * v_0 + 0.0f * _S77 + (s_diff_r2_0 + (_S78 + _S78)) * _S68 + s_diff_r2_0 * _S69, _S79 * v_0 + 0.0f * _S80 + (s_diff_r2_0 + (_S83 + 0.0f * _S81)) * _S67 + s_diff_r2_0 * _S70);
        float _S85 = 0.0f * u_0;
        float s_diff_r2_1 = _S85 + _S85 + (v_0 + v_0);
        float2  _S86 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S75 + (s_diff_r2_1 * _S74 + (s_diff_r2_1 * _S73 + s_diff_r2_1 * _S66 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S76 * v_0 + _S77 + (s_diff_r2_1 + (_S85 + 0.0f * _S78)) * _S68 + s_diff_r2_1 * _S69, 0.0f * _S79 * v_0 + _S80 + (s_diff_r2_1 + (_S81 + _S81)) * _S67 + s_diff_r2_1 * _S70);
        Matrix<float, 2, 2>  _S87 = transpose_0(makeMatrix<float, 2, 2> (_S84 + make_float2 (_S84.x * _S71 + _S84.y * _S72, 0.0f), _S86 + make_float2 (_S86.x * _S71 + _S86.y * _S72, 0.0f)));
        float inv_det_0 = 1.0f / (_S87.rows[int(0)].x * _S87.rows[int(1)].y - _S87.rows[int(0)].y * _S87.rows[int(1)].x);
        float _S88 = r_1.x;
        float _S89 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S88 * _S87.rows[int(1)].y - _S89 * _S87.rows[int(0)].y) * inv_det_0, (- _S88 * _S87.rows[int(1)].x + _S89 * _S87.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S90 = (*dist_coeffs_0)[int(0)];
    float _S91 = (*dist_coeffs_0)[int(1)];
    float _S92 = (*dist_coeffs_0)[int(2)];
    float _S93 = (*dist_coeffs_0)[int(3)];
    float _S94 = (*dist_coeffs_0)[int(4)];
    float _S95 = (*dist_coeffs_0)[int(5)];
    float _S96 = (*dist_coeffs_0)[int(6)];
    float _S97 = (*dist_coeffs_0)[int(7)];
    float _S98 = (*dist_coeffs_0)[int(8)];
    float _S99 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S100 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S100 + _S100);
    float _S101 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S102 = (*dist_coeffs_0)[int(1)] + r2_1 * _S101;
    float _S103 = (*dist_coeffs_0)[int(0)] + r2_1 * _S102;
    float radial_1 = 1.0f + r2_1 * _S103;
    float _S104 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S105 = _S104 * u_1;
    float _S106 = 2.0f * u_1;
    float _S107 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S108 = _S107 * u_1;
    float _S109 = 2.0f * v_1;
    float2  _S110 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S103 + (s_diff_r2_2 * _S102 + (s_diff_r2_2 * _S101 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S104 * v_1 + 0.0f * _S105 + (s_diff_r2_2 + (_S106 + _S106)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S107 * v_1 + 0.0f * _S108 + (s_diff_r2_2 + (_S100 + 0.0f * _S109)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S111 = 0.0f * u_1;
    float s_diff_r2_3 = _S111 + _S111 + (v_1 + v_1);
    float2  _S112 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S103 + (s_diff_r2_3 * _S102 + (s_diff_r2_3 * _S101 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S104 * v_1 + _S105 + (s_diff_r2_3 + (_S111 + 0.0f * _S106)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S107 * v_1 + _S108 + (s_diff_r2_3 + (_S109 + _S109)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S113 = transpose_0(makeMatrix<float, 2, 2> (_S110 + make_float2 (_S110.x * (*dist_coeffs_0)[int(8)] + _S110.y * (*dist_coeffs_0)[int(9)], 0.0f), _S112 + make_float2 (_S112.x * (*dist_coeffs_0)[int(8)] + _S112.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S114;
    if((F32_min((determinant_0(_S113)), ((F32_min((_S113.rows[int(0)].x), (_S113.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S115 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S90 + r2_2 * (_S91 + r2_2 * (_S92 + r2_2 * _S93)))) + make_float2 (_S104 * u_2 * v_2 + _S95 * (r2_2 + 2.0f * u_2 * u_2) + _S96 * r2_2, _S107 * u_2 * v_2 + _S94 * (r2_2 + 2.0f * v_2 * v_2) + _S97 * r2_2);
        _S114 = (length_0(_S115 + make_float2 (_S98 * _S115.x + _S99 * _S115.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S114 = false;
    }
    return _S114;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S116 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S116;
    FixedArray<float, 10>  _S117 = dist_coeffs_1;
    bool _S118 = undistort_point_0(_S116, &_S117, int(12), &uv_1);
    if(!_S118)
    {
        int3  _S119 = make_int3 (int(0));
        float3  _S120 = make_float3 ((float)_S119.x, (float)_S119.y, (float)_S119.z);
        return _S120;
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

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_5)
{
    float _S121 = dOut_5.y;
    float _S122 = dOut_5.z;
    float _S123 = dOut_5.x;
    float _S124 = (*a_0).primal_0.z * _S121 + - (*a_0).primal_0.y * _S122;
    float _S125 = - (*a_0).primal_0.z * _S123 + (*a_0).primal_0.x * _S122;
    float _S126 = (*a_0).primal_0.y * _S123 + - (*a_0).primal_0.x * _S121;
    float3  _S127 = make_float3 (- (*b_0).primal_0.z * _S121 + (*b_0).primal_0.y * _S122, (*b_0).primal_0.z * _S123 + - (*b_0).primal_0.x * _S122, - (*b_0).primal_0.y * _S123 + (*b_0).primal_0.x * _S121);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S127;
    float3  _S128 = make_float3 (_S124, _S125, _S126);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S128;
    return;
}

inline __device__ float3  cross_0(float3  left_0, float3  right_0)
{
    float _S129 = left_0.y;
    float _S130 = right_0.z;
    float _S131 = left_0.z;
    float _S132 = right_0.y;
    float _S133 = right_0.x;
    float _S134 = left_0.x;
    return make_float3 (_S129 * _S130 - _S131 * _S132, _S131 * _S133 - _S134 * _S130, _S134 * _S132 - _S129 * _S133);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S135, float3  _S136)
{
    return cross_0(_S135, _S136);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S137, float3  _S138)
{
    return dot_0(_S137, _S138);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S139, float _S140)
{
    _d_sqrt_0(_S139, _S140);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float _s_dOut_2)
{
    float _S141 = (*dpx_6).primal_0.x;
    float _S142 = (*dpx_6).primal_0.y;
    float _S143 = (*dpx_6).primal_0.z;
    DiffPair_float_0 _S144;
    (&_S144)->primal_0 = _S141 * _S141 + _S142 * _S142 + _S143 * _S143;
    (&_S144)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S144, _s_dOut_2);
    float _S145 = (*dpx_6).primal_0.z * _S144.differential_0;
    float _S146 = _S145 + _S145;
    float _S147 = (*dpx_6).primal_0.y * _S144.differential_0;
    float _S148 = _S147 + _S147;
    float _S149 = (*dpx_6).primal_0.x * _S144.differential_0;
    float _S150 = _S149 + _S149;
    float3  _S151 = make_float3 (0.0f);
    *&((&_S151)->z) = _S146;
    *&((&_S151)->y) = _S148;
    *&((&_S151)->x) = _S150;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S151;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S152, float _S153)
{
    s_bwd_prop_length_impl_0(_S152, _S153);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S154, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S155, float _S156)
{
    _d_dot_0(_S154, _S155, _S156);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S157, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S158, float3  _S159)
{
    _d_cross_0(_S157, _S158, _S159);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_3)
{
    float3  _S160 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S161 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S162 = s_primal_ctx_cross_0(dx_0, _S161);
    bool _S163 = (s_primal_ctx_dot_0(_S162, _S162)) != 0.0f;
    float3  _S164;
    float3  _S165;
    if(_S163)
    {
        float _S166 = length_1(_S162);
        float3  _S167 = make_float3 (_S166);
        _S164 = make_float3 (_S166 * _S166);
        _S165 = _S167;
    }
    else
    {
        _S164 = _S160;
        _S165 = _S160;
    }
    if(_S163)
    {
        float3  _S168 = _s_dOut_3 / _S164;
        float3  _S169 = _S162 * - _S168;
        float3  _S170 = _S165 * _S168;
        float _S171 = _S169.x + _S169.y + _S169.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S172;
        (&_S172)->primal_0 = _S162;
        (&_S172)->differential_0 = _S160;
        s_bwd_length_impl_0(&_S172, _S171);
        _S164 = _S170 + _S172.differential_0;
    }
    else
    {
        _S164 = _s_dOut_3;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S173;
    (&_S173)->primal_0 = _S162;
    (&_S173)->differential_0 = _S160;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S174;
    (&_S174)->primal_0 = _S162;
    (&_S174)->differential_0 = _S160;
    s_bwd_prop_dot_0(&_S173, &_S174, 0.0f);
    float3  _S175 = _S174.differential_0 + _S173.differential_0 + _S164;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S176;
    (&_S176)->primal_0 = dx_0;
    (&_S176)->differential_0 = _S160;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S177;
    (&_S177)->primal_0 = _S161;
    (&_S177)->differential_0 = _S160;
    s_bwd_prop_cross_0(&_S176, &_S177, _S175);
    float3  s_diff_dy_T_0 = - _S177.differential_0;
    float3  _S178 = - s_diff_dy_T_0;
    float3  _S179 = - _S176.differential_0;
    FixedArray<float3 , 4>  _S180;
    _S180[int(0)] = _S160;
    _S180[int(1)] = _S160;
    _S180[int(2)] = _S160;
    _S180[int(3)] = _S160;
    _S180[int(2)] = _S178;
    _S180[int(3)] = s_diff_dy_T_0;
    _S180[int(0)] = _S179;
    _S180[int(1)] = _S176.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S180;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S181, float3  _S182)
{
    s_bwd_prop_points_to_normal_0(_S181, _S182);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S183 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S183;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S184 = float2 {intrins_1.z, intrins_1.w};
    float2  _S185 = float2 {intrins_1.x, intrins_1.y};
    float2  _S186 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S184) / _S185;
    float2  uv_2 = _S186;
    FixedArray<float, 10>  _S187 = dist_coeffs_2;
    bool _S188 = undistort_point_0(_S186, &_S187, int(12), &uv_2);
    if(!_S188)
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
    float2  _S189 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S184) / _S185;
    float2  uv_3 = _S189;
    FixedArray<float, 10>  _S190 = dist_coeffs_2;
    bool _S191 = undistort_point_0(_S189, &_S190, int(12), &uv_3);
    if(!_S191)
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
    float2  _S192 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S184) / _S185;
    float2  uv_4 = _S192;
    FixedArray<float, 10>  _S193 = dist_coeffs_2;
    bool _S194 = undistort_point_0(_S192, &_S193, int(12), &uv_4);
    if(!_S194)
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
    float2  _S195 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S184) / _S185;
    float2  uv_5 = _S195;
    FixedArray<float, 10>  _S196 = dist_coeffs_2;
    bool _S197 = undistort_point_0(_S195, &_S196, int(12), &uv_5);
    if(!_S197)
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
    float2  _S198;
    bool _S199;
    float2  _S200;
    bool _S201;
    float2  _S202;
    bool _S203;
    float2  _S204;
    bool _S205;
};

inline __device__ float s_primal_ctx_sin_0(float _S206)
{
    return (F32_sin((_S206)));
}

inline __device__ float s_primal_ctx_cos_0(float _S207)
{
    return (F32_cos((_S207)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S208 = make_float2 (0.0f);
    _s_diff_ctx_0->_S198 = _S208;
    _s_diff_ctx_0->_S199 = false;
    _s_diff_ctx_0->_S200 = _S208;
    _s_diff_ctx_0->_S201 = false;
    _s_diff_ctx_0->_S202 = _S208;
    _s_diff_ctx_0->_S203 = false;
    _s_diff_ctx_0->_S204 = _S208;
    _s_diff_ctx_0->_S205 = false;
    _s_diff_ctx_0->_S200 = _S208;
    _s_diff_ctx_0->_S201 = false;
    _s_diff_ctx_0->_S202 = _S208;
    _s_diff_ctx_0->_S203 = false;
    _s_diff_ctx_0->_S204 = _S208;
    _s_diff_ctx_0->_S205 = false;
    float3  _S209 = make_float3 (0.0f);
    float2  _S210 = float2 {intrins_2.z, intrins_2.w};
    float2  _S211 = float2 {intrins_2.x, intrins_2.y};
    float2  _S212 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S210) / _S211;
    float2  _S213 = _S212;
    bool _S214 = undistort_point_0(_S212, dist_coeffs_3, int(12), &_S213);
    _s_diff_ctx_0->_S198 = _S213;
    _s_diff_ctx_0->_S199 = _S214;
    float2  uv_6 = _S213;
    bool _S215 = !_S214;
    float3  normal_4;
    if(_S215)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S216 = !_S215;
    int _S217;
    FixedArray<float3 , 4>  points_3;
    if(_S216)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S218 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S218), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S218))).x, (uv_6 / make_float2 ((F32_max((_S218), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S218))).y, s_primal_ctx_cos_0(_S218));
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
        float3  _S219 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S220 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S210) / _S211;
        float2  _S221 = _S220;
        bool _S222 = undistort_point_0(_S220, dist_coeffs_3, int(12), &_S221);
        _s_diff_ctx_0->_S200 = _S221;
        _s_diff_ctx_0->_S201 = _S222;
        float2  uv_7 = _S221;
        bool _S223 = !_S222;
        if(_S223)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S224 = !_S223;
        if(_S224)
        {
            if(is_fisheye_2)
            {
                float _S225 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S225), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S225))).x, (uv_7 / make_float2 ((F32_max((_S225), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S225))).y, s_primal_ctx_cos_0(_S225));
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
            float3  _S226 = make_float3 (dpdepths_0.y) * raydir_12;
            _S217 = int(2);
            points_3[int(0)] = _S219;
            points_3[int(1)] = _S226;
            points_3[int(2)] = _S209;
            points_3[int(3)] = _S209;
        }
        else
        {
            _S217 = int(0);
            points_3[int(0)] = _S219;
            points_3[int(1)] = _S209;
            points_3[int(2)] = _S209;
            points_3[int(3)] = _S209;
        }
        bool _runFlag_0;
        if(_S217 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S216;
            _S217 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S227 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S210) / _S211;
            float2  _S228 = _S227;
            bool _S229 = undistort_point_0(_S227, dist_coeffs_3, int(12), &_S228);
            _s_diff_ctx_0->_S202 = _S228;
            _s_diff_ctx_0->_S203 = _S229;
            float2  uv_8 = _S228;
            if(!_S229)
            {
                float3  _S230 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S217 = int(0);
                normal_4 = _S230;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S231 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S231), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S231))).x, (uv_8 / make_float2 ((F32_max((_S231), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S231))).y, s_primal_ctx_cos_0(_S231));
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
                float2  _S232 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S210) / _S211;
                float2  _S233 = _S232;
                bool _S234 = undistort_point_0(_S232, dist_coeffs_3, int(12), &_S233);
                _s_diff_ctx_0->_S204 = _S233;
                _s_diff_ctx_0->_S205 = _S234;
                float2  uv_9 = _S233;
                bool _S235 = !_S234;
                if(_S235)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S236 = !_S235;
                int _S237;
                if(_S236)
                {
                    if(is_fisheye_2)
                    {
                        float _S238 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S238), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S238))).x, (uv_9 / make_float2 ((F32_max((_S238), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S238))).y, s_primal_ctx_cos_0(_S238));
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
                    _S237 = int(2);
                }
                else
                {
                    _S237 = int(0);
                }
                if(_S237 != int(2))
                {
                    _runFlag_0 = false;
                    _S217 = _S237;
                }
                if(_runFlag_0)
                {
                    _S217 = int(1);
                }
            }
        }
    }
    else
    {
        _S217 = int(0);
        points_3[int(0)] = _S209;
        points_3[int(1)] = _S209;
        points_3[int(2)] = _S209;
        points_3[int(3)] = _S209;
    }
    if(!(_S217 != int(1)))
    {
        float3  _S239 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S239, _S239)) != 0.0f)
        {
            normal_4 = _S239 / make_float3 (length_1(_S239));
        }
        else
        {
            normal_4 = _S239;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_4, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S240 = *dpdepths_1;
    float3  _S241 = make_float3 (0.0f);
    float2  _S242 = _s_diff_ctx_1->_S198;
    bool _S243 = !!_s_diff_ctx_1->_S199;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S244;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S245;
    if(_S243)
    {
        if(is_fisheye_3)
        {
            float _S246 = length_0(_S242);
            float3  raydir_25 = make_float3 ((_S242 / make_float2 ((F32_max((_S246), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S246))).x, (_S242 / make_float2 ((F32_max((_S246), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S246))).y, s_primal_ctx_cos_0(_S246));
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
            float3  raydir_26 = make_float3 (_S242.x, _S242.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S247 = make_float3 (_S240.primal_0.x) * raydir_21;
        float2  _S248 = _s_diff_ctx_1->_S200;
        bool _S249 = !!_s_diff_ctx_1->_S201;
        if(_S249)
        {
            if(is_fisheye_3)
            {
                float _S250 = length_0(_S248);
                float3  raydir_27 = make_float3 ((_S248 / make_float2 ((F32_max((_S250), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S250))).x, (_S248 / make_float2 ((F32_max((_S250), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S250))).y, s_primal_ctx_cos_0(_S250));
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
                float3  raydir_28 = make_float3 (_S248.x, _S248.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S251 = make_float3 (_S240.primal_0.y) * raydir_22;
            _S244 = int(2);
            points_4[int(0)] = _S247;
            points_4[int(1)] = _S251;
            points_4[int(2)] = _S241;
            points_4[int(3)] = _S241;
        }
        else
        {
            _S244 = int(0);
            points_4[int(0)] = _S247;
            points_4[int(1)] = _S241;
            points_4[int(2)] = _S241;
            points_4[int(3)] = _S241;
            raydir_22 = _S241;
        }
        if(_S244 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S243;
            _S244 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S252 = _s_diff_ctx_1->_S202;
            if(!_s_diff_ctx_1->_S203)
            {
                _runFlag_2 = false;
                _S244 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S253 = length_0(_S252);
                    float3  raydir_29 = make_float3 ((_S252 / make_float2 ((F32_max((_S253), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S253))).x, (_S252 / make_float2 ((F32_max((_S253), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S253))).y, s_primal_ctx_cos_0(_S253));
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
                    float3  raydir_30 = make_float3 (_S252.x, _S252.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S240.primal_0.z) * raydir_23;
                float2  _S254 = _s_diff_ctx_1->_S204;
                bool _S255 = !!_s_diff_ctx_1->_S205;
                int _S256;
                if(_S255)
                {
                    if(is_fisheye_3)
                    {
                        float _S257 = length_0(_S254);
                        float3  raydir_31 = make_float3 ((_S254 / make_float2 ((F32_max((_S257), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S257))).x, (_S254 / make_float2 ((F32_max((_S257), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S257))).y, s_primal_ctx_cos_0(_S257));
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
                        float3  raydir_32 = make_float3 (_S254.x, _S254.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S240.primal_0.w) * raydir_24;
                    _S256 = int(2);
                }
                else
                {
                    _S256 = int(0);
                    raydir_24 = _S241;
                }
                if(_S256 != int(2))
                {
                    _runFlag_3 = false;
                    _S244 = _S256;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S244 = int(1);
                }
                float3  _S258 = raydir_23;
                _runFlag_3 = _S255;
                raydir_23 = raydir_24;
                raydir_24 = _S258;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S241;
                raydir_24 = _S241;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S241;
            raydir_24 = _S241;
        }
        float3  _S259 = raydir_21;
        float3  _S260 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S245 = _S249;
        raydir_23 = _S260;
        raydir_24 = _S259;
    }
    else
    {
        _S244 = int(0);
        points_4[int(0)] = _S241;
        points_4[int(1)] = _S241;
        points_4[int(2)] = _S241;
        points_4[int(3)] = _S241;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S241;
        raydir_22 = _S241;
        _S245 = false;
        raydir_23 = _S241;
        raydir_24 = _S241;
    }
    bool _S261 = !(_S244 != int(1));
    float3  _S262;
    float3  _S263;
    float3  _S264;
    float3  _S265;
    float3  _S266;
    bool _S267;
    if(_S261)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S268 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S269 = s_primal_ctx_cross_0(dx_1, _S268);
        bool _S270 = (s_primal_ctx_dot_0(_S269, _S269)) != 0.0f;
        if(_S270)
        {
            float _S271 = length_1(_S269);
            float3  _S272 = make_float3 (_S271);
            _S262 = make_float3 (_S271 * _S271);
            _S263 = _S272;
        }
        else
        {
            _S262 = _S241;
            _S263 = _S241;
        }
        float3  _S273 = _S263;
        _S267 = _S270;
        _S263 = _S269;
        _S264 = _S273;
        _S265 = dx_1;
        _S266 = _S268;
    }
    else
    {
        _S267 = false;
        _S262 = _S241;
        _S263 = _S241;
        _S264 = _S241;
        _S265 = _S241;
        _S266 = _S241;
    }
    float4  _S274 = make_float4 (0.0f);
    if(_S261)
    {
        if(_S267)
        {
            float3  _S275 = _s_dOut_4 / _S262;
            float3  _S276 = _S263 * - _S275;
            float3  _S277 = _S264 * _S275;
            float _S278 = _S276.x + _S276.y + _S276.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S279;
            (&_S279)->primal_0 = _S263;
            (&_S279)->differential_0 = _S241;
            s_bwd_length_impl_0(&_S279, _S278);
            _S262 = _S277 + _S279.differential_0;
        }
        else
        {
            _S262 = _s_dOut_4;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S280;
        (&_S280)->primal_0 = _S263;
        (&_S280)->differential_0 = _S241;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S281;
        (&_S281)->primal_0 = _S263;
        (&_S281)->differential_0 = _S241;
        s_bwd_prop_dot_0(&_S280, &_S281, 0.0f);
        float3  _S282 = _S281.differential_0 + _S280.differential_0 + _S262;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S283;
        (&_S283)->primal_0 = _S265;
        (&_S283)->differential_0 = _S241;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S284;
        (&_S284)->primal_0 = _S266;
        (&_S284)->differential_0 = _S241;
        s_bwd_prop_cross_0(&_S283, &_S284, _S282);
        float3  s_diff_dy_T_1 = - _S284.differential_0;
        float3  _S285 = - s_diff_dy_T_1;
        float3  _S286 = - _S283.differential_0;
        FixedArray<float3 , 4>  _S287;
        _S287[int(0)] = _S241;
        _S287[int(1)] = _S241;
        _S287[int(2)] = _S241;
        _S287[int(3)] = _S241;
        _S287[int(2)] = _S285;
        _S287[int(3)] = s_diff_dy_T_1;
        _S287[int(0)] = _S286;
        _S287[int(1)] = _S283.differential_0;
        points_4[int(0)] = _S287[int(0)];
        points_4[int(1)] = _S287[int(1)];
        points_4[int(2)] = _S287[int(2)];
        points_4[int(3)] = _S287[int(3)];
    }
    else
    {
        points_4[int(0)] = _S241;
        points_4[int(1)] = _S241;
        points_4[int(2)] = _S241;
        points_4[int(3)] = _S241;
    }
    float4  _S288;
    if(_S243)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S289 = points_4;
                FixedArray<float3 , 4>  _S290 = points_4;
                FixedArray<float3 , 4>  _S291 = points_4;
                FixedArray<float3 , 4>  _S292 = points_4;
                if(_runFlag_3)
                {
                    float3  _S293 = raydir_21 * _S292[int(3)];
                    float _S294 = _S293.x + _S293.y + _S293.z;
                    float4  _S295 = _S274;
                    *&((&_S295)->w) = _S294;
                    points_4[int(0)] = _S289[int(0)];
                    points_4[int(1)] = _S290[int(1)];
                    points_4[int(2)] = _S291[int(2)];
                    points_4[int(3)] = _S241;
                    _S288 = _S295;
                }
                else
                {
                    points_4[int(0)] = _S289[int(0)];
                    points_4[int(1)] = _S290[int(1)];
                    points_4[int(2)] = _S291[int(2)];
                    points_4[int(3)] = _S292[int(3)];
                    _S288 = _S274;
                }
                float3  _S296 = raydir_22 * points_4[int(2)];
                float _S297 = _S296.x + _S296.y + _S296.z;
                FixedArray<float3 , 4>  _S298 = points_4;
                FixedArray<float3 , 4>  _S299 = points_4;
                float4  _S300 = _S274;
                *&((&_S300)->z) = _S297;
                float4  _S301 = _S288 + _S300;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S298[int(1)];
                points_4[int(2)] = _S241;
                points_4[int(3)] = _S299[int(3)];
                _S288 = _S301;
            }
            else
            {
                FixedArray<float3 , 4>  _S302 = points_4;
                FixedArray<float3 , 4>  _S303 = points_4;
                FixedArray<float3 , 4>  _S304 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S302[int(1)];
                points_4[int(2)] = _S303[int(2)];
                points_4[int(3)] = _S304[int(3)];
                _S288 = _S274;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S305 = points_4;
            FixedArray<float3 , 4>  _S306 = points_4;
            FixedArray<float3 , 4>  _S307 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S305[int(1)];
            points_4[int(2)] = _S306[int(2)];
            points_4[int(3)] = _S307[int(3)];
            _S288 = _S274;
        }
        if(_S245)
        {
            FixedArray<float3 , 4>  _S308 = points_4;
            float3  _S309 = raydir_23 * points_4[int(1)];
            float _S310 = _S309.x + _S309.y + _S309.z;
            float4  _S311 = _S274;
            *&((&_S311)->y) = _S310;
            float4  _S312 = _S288 + _S311;
            points_4[int(0)] = _S241;
            points_4[int(1)] = _S241;
            points_4[int(2)] = _S241;
            points_4[int(3)] = _S241;
            raydir_21 = _S308[int(0)];
            _S288 = _S312;
        }
        else
        {
            FixedArray<float3 , 4>  _S313 = points_4;
            FixedArray<float3 , 4>  _S314 = points_4;
            FixedArray<float3 , 4>  _S315 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S313[int(1)];
            points_4[int(2)] = _S314[int(2)];
            points_4[int(3)] = _S315[int(3)];
            raydir_21 = _S241;
        }
        float3  _S316 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S317 = _S316.x + _S316.y + _S316.z;
        float4  _S318 = _S274;
        *&((&_S318)->x) = _S317;
        _S288 = _S288 + _S318;
    }
    else
    {
        _S288 = _S274;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S288;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S319, float4  _S320, FixedArray<float, 10>  * _S321, bool _S322, bool _S323, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S324, float3  _S325)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S326;
    float3  _S327 = s_primal_ctx_depth_to_normal_0(_S319, _S320, _S321, _S322, _S323, (*_S324).primal_0, &_S326);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S328 = _S326;
    s_bwd_prop_depth_to_normal_0(_S319, _S320, _S321, _S322, _S323, _S324, _S325, &_S328);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S329 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S329;
    FixedArray<float, 10>  _S330 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S330, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S331 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S331;
    FixedArray<float, 10>  _S332 = dist_coeffs_6;
    bool _S333 = undistort_point_0(_S331, &_S332, int(12), &uv_10);
    if(!_S333)
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

