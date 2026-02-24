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

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2)
{
    float3  _S29 = rgb_2;
    float _S30;
    if((rgb_2.x) < 0.00313080009073019f)
    {
        _S30 = _S29.x * 12.92000007629394531f;
    }
    else
    {
        _S30 = 1.0549999475479126f * (F32_pow((_S29.x), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S29)->x) = _S30;
    if((_S29.y) < 0.00313080009073019f)
    {
        _S30 = _S29.y * 12.92000007629394531f;
    }
    else
    {
        _S30 = 1.0549999475479126f * (F32_pow((_S29.y), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S29)->y) = _S30;
    if((_S29.z) < 0.00313080009073019f)
    {
        _S30 = _S29.z * 12.92000007629394531f;
    }
    else
    {
        _S30 = 1.0549999475479126f * (F32_pow((_S29.z), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    *&((&_S29)->z) = _S30;
    return _S29;
}

inline __device__ float s_primal_ctx_pow_0(float _S31, float _S32)
{
    return (F32_pow((_S31), (_S32)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S33, DiffPair_float_0 * _S34, float _S35)
{
    _d_pow_0(_S33, _S34, _S35);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, float3  _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S36 = *dprgb_1;
    float _S37 = (*dprgb_1).primal_0.x;
    bool _S38 = _S37 < 0.00313080009073019f;
    float _S39;
    if(_S38)
    {
        _S39 = _S37 * 12.92000007629394531f;
    }
    else
    {
        _S39 = 1.0549999475479126f * s_primal_ctx_pow_0(_S37, 0.4166666567325592f) - 0.05499999970197678f;
    }
    float3  _S40 = _S36.primal_0;
    *&((&_S40)->x) = _S39;
    float _S41 = _S40.y;
    bool _S42 = _S41 < 0.00313080009073019f;
    if(_S42)
    {
        _S39 = _S41 * 12.92000007629394531f;
    }
    else
    {
        _S39 = 1.0549999475479126f * s_primal_ctx_pow_0(_S41, 0.4166666567325592f) - 0.05499999970197678f;
    }
    *&((&_S40)->y) = _S39;
    float _S43 = _S40.z;
    bool _S44 = _S43 < 0.00313080009073019f;
    _S40 = _s_dOut_1;
    *&((&_S40)->z) = 0.0f;
    if(_S44)
    {
        _S39 = 12.92000007629394531f * _s_dOut_1.z;
    }
    else
    {
        float _S45 = 1.0549999475479126f * _s_dOut_1.z;
        DiffPair_float_0 _S46;
        (&_S46)->primal_0 = _S43;
        (&_S46)->differential_0 = 0.0f;
        DiffPair_float_0 _S47;
        (&_S47)->primal_0 = 0.4166666567325592f;
        (&_S47)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S46, &_S47, _S45);
        _S39 = _S46.differential_0;
    }
    float3  _S48 = _S40 + make_float3 (0.0f, 0.0f, _S39);
    _S40 = _S48;
    *&((&_S40)->y) = 0.0f;
    if(_S42)
    {
        _S39 = 12.92000007629394531f * _S48.y;
    }
    else
    {
        float _S49 = 1.0549999475479126f * _S48.y;
        DiffPair_float_0 _S50;
        (&_S50)->primal_0 = _S41;
        (&_S50)->differential_0 = 0.0f;
        DiffPair_float_0 _S51;
        (&_S51)->primal_0 = 0.4166666567325592f;
        (&_S51)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S50, &_S51, _S49);
        _S39 = _S50.differential_0;
    }
    float3  _S52 = _S40 + make_float3 (0.0f, _S39, 0.0f);
    _S40 = _S52;
    *&((&_S40)->x) = 0.0f;
    if(_S38)
    {
        _S39 = 12.92000007629394531f * _S52.x;
    }
    else
    {
        float _S53 = 1.0549999475479126f * _S52.x;
        DiffPair_float_0 _S54;
        (&_S54)->primal_0 = _S37;
        (&_S54)->differential_0 = 0.0f;
        DiffPair_float_0 _S55;
        (&_S55)->primal_0 = 0.4166666567325592f;
        (&_S55)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S54, &_S55, _S53);
        _S39 = _S54.differential_0;
    }
    float3  _S56 = _S40 + make_float3 (_S39, 0.0f, 0.0f);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S56;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S57, float3  _S58)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S57, _S58);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, float3  v_out_rgb_1)
{
    float3  _S59 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S59;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_3)
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
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_0)), c_0) = _slang_vector_get_element(x_3.rows[c_0], r_0);
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

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_3, float dOut_3)
{
    float _S60 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_3).primal_0)))))) * dOut_3;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S60;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, float dOut_4)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_4;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_4).primal_0.x * dOut_4;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_4;
    *&((&y_d_result_0)->y) = (*dpx_4).primal_0.y * dOut_4;
    *&((&x_d_result_0)->z) = (*dpy_2).primal_0.z * dOut_4;
    *&((&y_d_result_0)->z) = (*dpx_4).primal_0.z * dOut_4;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_4, float3  y_2)
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
        float result_4 = result_3 + _slang_vector_get_element(x_4, i_2) * _slang_vector_get_element(y_2, i_2);
        i_2 = i_2 + int(1);
        result_3 = result_4;
    }
    return result_3;
}

inline __device__ float dot_1(float2  x_5, float2  y_3)
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
        float result_6 = result_5 + _slang_vector_get_element(x_5, i_3) * _slang_vector_get_element(y_3, i_3);
        i_3 = i_3 + int(1);
        result_5 = result_6;
    }
    return result_5;
}

inline __device__ float length_0(float2  x_6)
{
    return (F32_sqrt((dot_1(x_6, x_6))));
}

inline __device__ float length_1(float3  x_7)
{
    return (F32_sqrt((dot_0(x_7, x_7))));
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
        float _S61 = (*dist_coeffs_0)[int(3)];
        float _S62 = (*dist_coeffs_0)[int(4)];
        float _S63 = (*dist_coeffs_0)[int(5)];
        float _S64 = (*dist_coeffs_0)[int(6)];
        float _S65 = (*dist_coeffs_0)[int(7)];
        float _S66 = (*dist_coeffs_0)[int(8)];
        float _S67 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S68 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S69 = (*dist_coeffs_0)[int(1)] + r2_0 * _S68;
        float _S70 = (*dist_coeffs_0)[int(0)] + r2_0 * _S69;
        float radial_0 = 1.0f + r2_0 * _S70;
        float _S71 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S72 = _S71 * u_0;
        float _S73 = 2.0f * u_0;
        float _S74 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S75 = _S74 * u_0;
        float _S76 = 2.0f * v_0;
        float2  _S77 = q_0 * make_float2 (radial_0) + make_float2 (_S72 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S73 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S75 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S76 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S77 + make_float2 ((*dist_coeffs_0)[int(8)] * _S77.x + (*dist_coeffs_0)[int(9)] * _S77.y, 0.0f) - uv_0;
        float _S78 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S78 + _S78);
        float2  _S79 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S70 + (s_diff_r2_0 * _S69 + (s_diff_r2_0 * _S68 + s_diff_r2_0 * _S61 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S71 * v_0 + 0.0f * _S72 + (s_diff_r2_0 + (_S73 + _S73)) * _S63 + s_diff_r2_0 * _S64, _S74 * v_0 + 0.0f * _S75 + (s_diff_r2_0 + (_S78 + 0.0f * _S76)) * _S62 + s_diff_r2_0 * _S65);
        float _S80 = 0.0f * u_0;
        float s_diff_r2_1 = _S80 + _S80 + (v_0 + v_0);
        float2  _S81 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S70 + (s_diff_r2_1 * _S69 + (s_diff_r2_1 * _S68 + s_diff_r2_1 * _S61 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S71 * v_0 + _S72 + (s_diff_r2_1 + (_S80 + 0.0f * _S73)) * _S63 + s_diff_r2_1 * _S64, 0.0f * _S74 * v_0 + _S75 + (s_diff_r2_1 + (_S76 + _S76)) * _S62 + s_diff_r2_1 * _S65);
        Matrix<float, 2, 2>  _S82 = transpose_0(makeMatrix<float, 2, 2> (_S79 + make_float2 (_S79.x * _S66 + _S79.y * _S67, 0.0f), _S81 + make_float2 (_S81.x * _S66 + _S81.y * _S67, 0.0f)));
        float inv_det_0 = 1.0f / (_S82.rows[int(0)].x * _S82.rows[int(1)].y - _S82.rows[int(0)].y * _S82.rows[int(1)].x);
        float _S83 = r_1.x;
        float _S84 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S83 * _S82.rows[int(1)].y - _S84 * _S82.rows[int(0)].y) * inv_det_0, (- _S83 * _S82.rows[int(1)].x + _S84 * _S82.rows[int(0)].x) * inv_det_0);
        i_4 = i_4 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S85 = (*dist_coeffs_0)[int(0)];
    float _S86 = (*dist_coeffs_0)[int(1)];
    float _S87 = (*dist_coeffs_0)[int(2)];
    float _S88 = (*dist_coeffs_0)[int(3)];
    float _S89 = (*dist_coeffs_0)[int(4)];
    float _S90 = (*dist_coeffs_0)[int(5)];
    float _S91 = (*dist_coeffs_0)[int(6)];
    float _S92 = (*dist_coeffs_0)[int(7)];
    float _S93 = (*dist_coeffs_0)[int(8)];
    float _S94 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S95 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S95 + _S95);
    float _S96 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S97 = (*dist_coeffs_0)[int(1)] + r2_1 * _S96;
    float _S98 = (*dist_coeffs_0)[int(0)] + r2_1 * _S97;
    float radial_1 = 1.0f + r2_1 * _S98;
    float _S99 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S100 = _S99 * u_1;
    float _S101 = 2.0f * u_1;
    float _S102 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S103 = _S102 * u_1;
    float _S104 = 2.0f * v_1;
    float2  _S105 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S98 + (s_diff_r2_2 * _S97 + (s_diff_r2_2 * _S96 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S99 * v_1 + 0.0f * _S100 + (s_diff_r2_2 + (_S101 + _S101)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S102 * v_1 + 0.0f * _S103 + (s_diff_r2_2 + (_S95 + 0.0f * _S104)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S106 = 0.0f * u_1;
    float s_diff_r2_3 = _S106 + _S106 + (v_1 + v_1);
    float2  _S107 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S98 + (s_diff_r2_3 * _S97 + (s_diff_r2_3 * _S96 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S99 * v_1 + _S100 + (s_diff_r2_3 + (_S106 + 0.0f * _S101)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S102 * v_1 + _S103 + (s_diff_r2_3 + (_S104 + _S104)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S108 = transpose_0(makeMatrix<float, 2, 2> (_S105 + make_float2 (_S105.x * (*dist_coeffs_0)[int(8)] + _S105.y * (*dist_coeffs_0)[int(9)], 0.0f), _S107 + make_float2 (_S107.x * (*dist_coeffs_0)[int(8)] + _S107.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S109;
    if((F32_min((determinant_0(_S108)), ((F32_min((_S108.rows[int(0)].x), (_S108.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S110 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S85 + r2_2 * (_S86 + r2_2 * (_S87 + r2_2 * _S88)))) + make_float2 (_S99 * u_2 * v_2 + _S90 * (r2_2 + 2.0f * u_2 * u_2) + _S91 * r2_2, _S102 * u_2 * v_2 + _S89 * (r2_2 + 2.0f * v_2 * v_2) + _S92 * r2_2);
        _S109 = (length_0(_S110 + make_float2 (_S93 * _S110.x + _S94 * _S110.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S109 = false;
    }
    return _S109;
}

inline __device__ float3  normalize_0(float3  x_8)
{
    return x_8 / make_float3 (length_1(x_8));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S111 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S111;
    FixedArray<float, 10>  _S112 = dist_coeffs_1;
    bool _S113 = undistort_point_0(_S111, &_S112, int(12), &uv_1);
    if(!_S113)
    {
        int3  _S114 = make_int3 (int(0));
        float3  _S115 = make_float3 ((float)_S114.x, (float)_S114.y, (float)_S114.z);
        return _S115;
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
    float _S116 = dOut_5.y;
    float _S117 = dOut_5.z;
    float _S118 = dOut_5.x;
    float _S119 = (*a_0).primal_0.z * _S116 + - (*a_0).primal_0.y * _S117;
    float _S120 = - (*a_0).primal_0.z * _S118 + (*a_0).primal_0.x * _S117;
    float _S121 = (*a_0).primal_0.y * _S118 + - (*a_0).primal_0.x * _S116;
    float3  _S122 = make_float3 (- (*b_0).primal_0.z * _S116 + (*b_0).primal_0.y * _S117, (*b_0).primal_0.z * _S118 + - (*b_0).primal_0.x * _S117, - (*b_0).primal_0.y * _S118 + (*b_0).primal_0.x * _S116);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S122;
    float3  _S123 = make_float3 (_S119, _S120, _S121);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S123;
    return;
}

inline __device__ float3  cross_0(float3  left_0, float3  right_0)
{
    float _S124 = left_0.y;
    float _S125 = right_0.z;
    float _S126 = left_0.z;
    float _S127 = right_0.y;
    float _S128 = right_0.x;
    float _S129 = left_0.x;
    return make_float3 (_S124 * _S125 - _S126 * _S127, _S126 * _S128 - _S129 * _S125, _S129 * _S127 - _S124 * _S128);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S130, float3  _S131)
{
    return cross_0(_S130, _S131);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S132, float3  _S133)
{
    return dot_0(_S132, _S133);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S134, float _S135)
{
    _d_sqrt_0(_S134, _S135);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float _s_dOut_2)
{
    float _S136 = (*dpx_5).primal_0.x;
    float _S137 = (*dpx_5).primal_0.y;
    float _S138 = (*dpx_5).primal_0.z;
    DiffPair_float_0 _S139;
    (&_S139)->primal_0 = _S136 * _S136 + _S137 * _S137 + _S138 * _S138;
    (&_S139)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S139, _s_dOut_2);
    float _S140 = (*dpx_5).primal_0.z * _S139.differential_0;
    float _S141 = _S140 + _S140;
    float _S142 = (*dpx_5).primal_0.y * _S139.differential_0;
    float _S143 = _S142 + _S142;
    float _S144 = (*dpx_5).primal_0.x * _S139.differential_0;
    float _S145 = _S144 + _S144;
    float3  _S146 = make_float3 (0.0f);
    *&((&_S146)->z) = _S141;
    *&((&_S146)->y) = _S143;
    *&((&_S146)->x) = _S145;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S146;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S147, float _S148)
{
    s_bwd_prop_length_impl_0(_S147, _S148);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S149, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S150, float _S151)
{
    _d_dot_0(_S149, _S150, _S151);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S152, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S153, float3  _S154)
{
    _d_cross_0(_S152, _S153, _S154);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_3)
{
    float3  _S155 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S156 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S157 = s_primal_ctx_cross_0(dx_0, _S156);
    bool _S158 = (s_primal_ctx_dot_0(_S157, _S157)) != 0.0f;
    float3  _S159;
    float3  _S160;
    if(_S158)
    {
        float _S161 = length_1(_S157);
        float3  _S162 = make_float3 (_S161);
        _S159 = make_float3 (_S161 * _S161);
        _S160 = _S162;
    }
    else
    {
        _S159 = _S155;
        _S160 = _S155;
    }
    if(_S158)
    {
        float3  _S163 = _s_dOut_3 / _S159;
        float3  _S164 = _S157 * - _S163;
        float3  _S165 = _S160 * _S163;
        float _S166 = _S164.x + _S164.y + _S164.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S167;
        (&_S167)->primal_0 = _S157;
        (&_S167)->differential_0 = _S155;
        s_bwd_length_impl_0(&_S167, _S166);
        _S159 = _S165 + _S167.differential_0;
    }
    else
    {
        _S159 = _s_dOut_3;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S168;
    (&_S168)->primal_0 = _S157;
    (&_S168)->differential_0 = _S155;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S169;
    (&_S169)->primal_0 = _S157;
    (&_S169)->differential_0 = _S155;
    s_bwd_prop_dot_0(&_S168, &_S169, 0.0f);
    float3  _S170 = _S169.differential_0 + _S168.differential_0 + _S159;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S171;
    (&_S171)->primal_0 = dx_0;
    (&_S171)->differential_0 = _S155;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S172;
    (&_S172)->primal_0 = _S156;
    (&_S172)->differential_0 = _S155;
    s_bwd_prop_cross_0(&_S171, &_S172, _S170);
    float3  s_diff_dy_T_0 = - _S172.differential_0;
    float3  _S173 = - s_diff_dy_T_0;
    float3  _S174 = - _S171.differential_0;
    FixedArray<float3 , 4>  _S175;
    _S175[int(0)] = _S155;
    _S175[int(1)] = _S155;
    _S175[int(2)] = _S155;
    _S175[int(3)] = _S155;
    _S175[int(2)] = _S173;
    _S175[int(3)] = s_diff_dy_T_0;
    _S175[int(0)] = _S174;
    _S175[int(1)] = _S171.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S175;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S176, float3  _S177)
{
    s_bwd_prop_points_to_normal_0(_S176, _S177);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S178 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S178;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S179 = float2 {intrins_1.z, intrins_1.w};
    float2  _S180 = float2 {intrins_1.x, intrins_1.y};
    float2  _S181 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S179) / _S180;
    float2  uv_2 = _S181;
    FixedArray<float, 10>  _S182 = dist_coeffs_2;
    bool _S183 = undistort_point_0(_S181, &_S182, int(12), &uv_2);
    if(!_S183)
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
    float2  _S184 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S179) / _S180;
    float2  uv_3 = _S184;
    FixedArray<float, 10>  _S185 = dist_coeffs_2;
    bool _S186 = undistort_point_0(_S184, &_S185, int(12), &uv_3);
    if(!_S186)
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
    float2  _S187 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S179) / _S180;
    float2  uv_4 = _S187;
    FixedArray<float, 10>  _S188 = dist_coeffs_2;
    bool _S189 = undistort_point_0(_S187, &_S188, int(12), &uv_4);
    if(!_S189)
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
    float2  _S190 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S179) / _S180;
    float2  uv_5 = _S190;
    FixedArray<float, 10>  _S191 = dist_coeffs_2;
    bool _S192 = undistort_point_0(_S190, &_S191, int(12), &uv_5);
    if(!_S192)
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
    float2  _S193;
    bool _S194;
    float2  _S195;
    bool _S196;
    float2  _S197;
    bool _S198;
    float2  _S199;
    bool _S200;
};

inline __device__ float s_primal_ctx_sin_0(float _S201)
{
    return (F32_sin((_S201)));
}

inline __device__ float s_primal_ctx_cos_0(float _S202)
{
    return (F32_cos((_S202)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S203 = make_float2 (0.0f);
    _s_diff_ctx_0->_S193 = _S203;
    _s_diff_ctx_0->_S194 = false;
    _s_diff_ctx_0->_S195 = _S203;
    _s_diff_ctx_0->_S196 = false;
    _s_diff_ctx_0->_S197 = _S203;
    _s_diff_ctx_0->_S198 = false;
    _s_diff_ctx_0->_S199 = _S203;
    _s_diff_ctx_0->_S200 = false;
    _s_diff_ctx_0->_S195 = _S203;
    _s_diff_ctx_0->_S196 = false;
    _s_diff_ctx_0->_S197 = _S203;
    _s_diff_ctx_0->_S198 = false;
    _s_diff_ctx_0->_S199 = _S203;
    _s_diff_ctx_0->_S200 = false;
    float3  _S204 = make_float3 (0.0f);
    float2  _S205 = float2 {intrins_2.z, intrins_2.w};
    float2  _S206 = float2 {intrins_2.x, intrins_2.y};
    float2  _S207 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S205) / _S206;
    float2  _S208 = _S207;
    bool _S209 = undistort_point_0(_S207, dist_coeffs_3, int(12), &_S208);
    _s_diff_ctx_0->_S193 = _S208;
    _s_diff_ctx_0->_S194 = _S209;
    float2  uv_6 = _S208;
    bool _S210 = !_S209;
    float3  normal_4;
    if(_S210)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S211 = !_S210;
    int _S212;
    FixedArray<float3 , 4>  points_3;
    if(_S211)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S213 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S213), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S213))).x, (uv_6 / make_float2 ((F32_max((_S213), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S213))).y, s_primal_ctx_cos_0(_S213));
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
        float3  _S214 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S215 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S205) / _S206;
        float2  _S216 = _S215;
        bool _S217 = undistort_point_0(_S215, dist_coeffs_3, int(12), &_S216);
        _s_diff_ctx_0->_S195 = _S216;
        _s_diff_ctx_0->_S196 = _S217;
        float2  uv_7 = _S216;
        bool _S218 = !_S217;
        if(_S218)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S219 = !_S218;
        if(_S219)
        {
            if(is_fisheye_2)
            {
                float _S220 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S220), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S220))).x, (uv_7 / make_float2 ((F32_max((_S220), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S220))).y, s_primal_ctx_cos_0(_S220));
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
            float3  _S221 = make_float3 (dpdepths_0.y) * raydir_12;
            _S212 = int(2);
            points_3[int(0)] = _S214;
            points_3[int(1)] = _S221;
            points_3[int(2)] = _S204;
            points_3[int(3)] = _S204;
        }
        else
        {
            _S212 = int(0);
            points_3[int(0)] = _S214;
            points_3[int(1)] = _S204;
            points_3[int(2)] = _S204;
            points_3[int(3)] = _S204;
        }
        bool _runFlag_0;
        if(_S212 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S211;
            _S212 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S222 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S205) / _S206;
            float2  _S223 = _S222;
            bool _S224 = undistort_point_0(_S222, dist_coeffs_3, int(12), &_S223);
            _s_diff_ctx_0->_S197 = _S223;
            _s_diff_ctx_0->_S198 = _S224;
            float2  uv_8 = _S223;
            if(!_S224)
            {
                float3  _S225 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S212 = int(0);
                normal_4 = _S225;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S226 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S226), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S226))).x, (uv_8 / make_float2 ((F32_max((_S226), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S226))).y, s_primal_ctx_cos_0(_S226));
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
                float2  _S227 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S205) / _S206;
                float2  _S228 = _S227;
                bool _S229 = undistort_point_0(_S227, dist_coeffs_3, int(12), &_S228);
                _s_diff_ctx_0->_S199 = _S228;
                _s_diff_ctx_0->_S200 = _S229;
                float2  uv_9 = _S228;
                bool _S230 = !_S229;
                if(_S230)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S231 = !_S230;
                int _S232;
                if(_S231)
                {
                    if(is_fisheye_2)
                    {
                        float _S233 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S233), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S233))).x, (uv_9 / make_float2 ((F32_max((_S233), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S233))).y, s_primal_ctx_cos_0(_S233));
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
                    _S232 = int(2);
                }
                else
                {
                    _S232 = int(0);
                }
                if(_S232 != int(2))
                {
                    _runFlag_0 = false;
                    _S212 = _S232;
                }
                if(_runFlag_0)
                {
                    _S212 = int(1);
                }
            }
        }
    }
    else
    {
        _S212 = int(0);
        points_3[int(0)] = _S204;
        points_3[int(1)] = _S204;
        points_3[int(2)] = _S204;
        points_3[int(3)] = _S204;
    }
    if(!(_S212 != int(1)))
    {
        float3  _S234 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S234, _S234)) != 0.0f)
        {
            normal_4 = _S234 / make_float3 (length_1(_S234));
        }
        else
        {
            normal_4 = _S234;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_4, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S235 = *dpdepths_1;
    float3  _S236 = make_float3 (0.0f);
    float2  _S237 = _s_diff_ctx_1->_S193;
    bool _S238 = !!_s_diff_ctx_1->_S194;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S239;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S240;
    if(_S238)
    {
        if(is_fisheye_3)
        {
            float _S241 = length_0(_S237);
            float3  raydir_25 = make_float3 ((_S237 / make_float2 ((F32_max((_S241), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S241))).x, (_S237 / make_float2 ((F32_max((_S241), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S241))).y, s_primal_ctx_cos_0(_S241));
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
            float3  raydir_26 = make_float3 (_S237.x, _S237.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S242 = make_float3 (_S235.primal_0.x) * raydir_21;
        float2  _S243 = _s_diff_ctx_1->_S195;
        bool _S244 = !!_s_diff_ctx_1->_S196;
        if(_S244)
        {
            if(is_fisheye_3)
            {
                float _S245 = length_0(_S243);
                float3  raydir_27 = make_float3 ((_S243 / make_float2 ((F32_max((_S245), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S245))).x, (_S243 / make_float2 ((F32_max((_S245), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S245))).y, s_primal_ctx_cos_0(_S245));
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
                float3  raydir_28 = make_float3 (_S243.x, _S243.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S246 = make_float3 (_S235.primal_0.y) * raydir_22;
            _S239 = int(2);
            points_4[int(0)] = _S242;
            points_4[int(1)] = _S246;
            points_4[int(2)] = _S236;
            points_4[int(3)] = _S236;
        }
        else
        {
            _S239 = int(0);
            points_4[int(0)] = _S242;
            points_4[int(1)] = _S236;
            points_4[int(2)] = _S236;
            points_4[int(3)] = _S236;
            raydir_22 = _S236;
        }
        if(_S239 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S238;
            _S239 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S247 = _s_diff_ctx_1->_S197;
            if(!_s_diff_ctx_1->_S198)
            {
                _runFlag_2 = false;
                _S239 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S248 = length_0(_S247);
                    float3  raydir_29 = make_float3 ((_S247 / make_float2 ((F32_max((_S248), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S248))).x, (_S247 / make_float2 ((F32_max((_S248), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S248))).y, s_primal_ctx_cos_0(_S248));
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
                    float3  raydir_30 = make_float3 (_S247.x, _S247.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S235.primal_0.z) * raydir_23;
                float2  _S249 = _s_diff_ctx_1->_S199;
                bool _S250 = !!_s_diff_ctx_1->_S200;
                int _S251;
                if(_S250)
                {
                    if(is_fisheye_3)
                    {
                        float _S252 = length_0(_S249);
                        float3  raydir_31 = make_float3 ((_S249 / make_float2 ((F32_max((_S252), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S252))).x, (_S249 / make_float2 ((F32_max((_S252), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S252))).y, s_primal_ctx_cos_0(_S252));
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
                        float3  raydir_32 = make_float3 (_S249.x, _S249.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S235.primal_0.w) * raydir_24;
                    _S251 = int(2);
                }
                else
                {
                    _S251 = int(0);
                    raydir_24 = _S236;
                }
                if(_S251 != int(2))
                {
                    _runFlag_3 = false;
                    _S239 = _S251;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S239 = int(1);
                }
                float3  _S253 = raydir_23;
                _runFlag_3 = _S250;
                raydir_23 = raydir_24;
                raydir_24 = _S253;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S236;
                raydir_24 = _S236;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S236;
            raydir_24 = _S236;
        }
        float3  _S254 = raydir_21;
        float3  _S255 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S240 = _S244;
        raydir_23 = _S255;
        raydir_24 = _S254;
    }
    else
    {
        _S239 = int(0);
        points_4[int(0)] = _S236;
        points_4[int(1)] = _S236;
        points_4[int(2)] = _S236;
        points_4[int(3)] = _S236;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S236;
        raydir_22 = _S236;
        _S240 = false;
        raydir_23 = _S236;
        raydir_24 = _S236;
    }
    bool _S256 = !(_S239 != int(1));
    float3  _S257;
    float3  _S258;
    float3  _S259;
    float3  _S260;
    float3  _S261;
    bool _S262;
    if(_S256)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S263 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S264 = s_primal_ctx_cross_0(dx_1, _S263);
        bool _S265 = (s_primal_ctx_dot_0(_S264, _S264)) != 0.0f;
        if(_S265)
        {
            float _S266 = length_1(_S264);
            float3  _S267 = make_float3 (_S266);
            _S257 = make_float3 (_S266 * _S266);
            _S258 = _S267;
        }
        else
        {
            _S257 = _S236;
            _S258 = _S236;
        }
        float3  _S268 = _S258;
        _S262 = _S265;
        _S258 = _S264;
        _S259 = _S268;
        _S260 = dx_1;
        _S261 = _S263;
    }
    else
    {
        _S262 = false;
        _S257 = _S236;
        _S258 = _S236;
        _S259 = _S236;
        _S260 = _S236;
        _S261 = _S236;
    }
    float4  _S269 = make_float4 (0.0f);
    if(_S256)
    {
        if(_S262)
        {
            float3  _S270 = _s_dOut_4 / _S257;
            float3  _S271 = _S258 * - _S270;
            float3  _S272 = _S259 * _S270;
            float _S273 = _S271.x + _S271.y + _S271.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S274;
            (&_S274)->primal_0 = _S258;
            (&_S274)->differential_0 = _S236;
            s_bwd_length_impl_0(&_S274, _S273);
            _S257 = _S272 + _S274.differential_0;
        }
        else
        {
            _S257 = _s_dOut_4;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S275;
        (&_S275)->primal_0 = _S258;
        (&_S275)->differential_0 = _S236;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S276;
        (&_S276)->primal_0 = _S258;
        (&_S276)->differential_0 = _S236;
        s_bwd_prop_dot_0(&_S275, &_S276, 0.0f);
        float3  _S277 = _S276.differential_0 + _S275.differential_0 + _S257;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S278;
        (&_S278)->primal_0 = _S260;
        (&_S278)->differential_0 = _S236;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S279;
        (&_S279)->primal_0 = _S261;
        (&_S279)->differential_0 = _S236;
        s_bwd_prop_cross_0(&_S278, &_S279, _S277);
        float3  s_diff_dy_T_1 = - _S279.differential_0;
        float3  _S280 = - s_diff_dy_T_1;
        float3  _S281 = - _S278.differential_0;
        FixedArray<float3 , 4>  _S282;
        _S282[int(0)] = _S236;
        _S282[int(1)] = _S236;
        _S282[int(2)] = _S236;
        _S282[int(3)] = _S236;
        _S282[int(2)] = _S280;
        _S282[int(3)] = s_diff_dy_T_1;
        _S282[int(0)] = _S281;
        _S282[int(1)] = _S278.differential_0;
        points_4[int(0)] = _S282[int(0)];
        points_4[int(1)] = _S282[int(1)];
        points_4[int(2)] = _S282[int(2)];
        points_4[int(3)] = _S282[int(3)];
    }
    else
    {
        points_4[int(0)] = _S236;
        points_4[int(1)] = _S236;
        points_4[int(2)] = _S236;
        points_4[int(3)] = _S236;
    }
    float4  _S283;
    if(_S238)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S284 = points_4;
                FixedArray<float3 , 4>  _S285 = points_4;
                FixedArray<float3 , 4>  _S286 = points_4;
                FixedArray<float3 , 4>  _S287 = points_4;
                if(_runFlag_3)
                {
                    float3  _S288 = raydir_21 * _S287[int(3)];
                    float _S289 = _S288.x + _S288.y + _S288.z;
                    float4  _S290 = _S269;
                    *&((&_S290)->w) = _S289;
                    points_4[int(0)] = _S284[int(0)];
                    points_4[int(1)] = _S285[int(1)];
                    points_4[int(2)] = _S286[int(2)];
                    points_4[int(3)] = _S236;
                    _S283 = _S290;
                }
                else
                {
                    points_4[int(0)] = _S284[int(0)];
                    points_4[int(1)] = _S285[int(1)];
                    points_4[int(2)] = _S286[int(2)];
                    points_4[int(3)] = _S287[int(3)];
                    _S283 = _S269;
                }
                float3  _S291 = raydir_22 * points_4[int(2)];
                float _S292 = _S291.x + _S291.y + _S291.z;
                FixedArray<float3 , 4>  _S293 = points_4;
                FixedArray<float3 , 4>  _S294 = points_4;
                float4  _S295 = _S269;
                *&((&_S295)->z) = _S292;
                float4  _S296 = _S283 + _S295;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S293[int(1)];
                points_4[int(2)] = _S236;
                points_4[int(3)] = _S294[int(3)];
                _S283 = _S296;
            }
            else
            {
                FixedArray<float3 , 4>  _S297 = points_4;
                FixedArray<float3 , 4>  _S298 = points_4;
                FixedArray<float3 , 4>  _S299 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S297[int(1)];
                points_4[int(2)] = _S298[int(2)];
                points_4[int(3)] = _S299[int(3)];
                _S283 = _S269;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S300 = points_4;
            FixedArray<float3 , 4>  _S301 = points_4;
            FixedArray<float3 , 4>  _S302 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S300[int(1)];
            points_4[int(2)] = _S301[int(2)];
            points_4[int(3)] = _S302[int(3)];
            _S283 = _S269;
        }
        if(_S240)
        {
            FixedArray<float3 , 4>  _S303 = points_4;
            float3  _S304 = raydir_23 * points_4[int(1)];
            float _S305 = _S304.x + _S304.y + _S304.z;
            float4  _S306 = _S269;
            *&((&_S306)->y) = _S305;
            float4  _S307 = _S283 + _S306;
            points_4[int(0)] = _S236;
            points_4[int(1)] = _S236;
            points_4[int(2)] = _S236;
            points_4[int(3)] = _S236;
            raydir_21 = _S303[int(0)];
            _S283 = _S307;
        }
        else
        {
            FixedArray<float3 , 4>  _S308 = points_4;
            FixedArray<float3 , 4>  _S309 = points_4;
            FixedArray<float3 , 4>  _S310 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S308[int(1)];
            points_4[int(2)] = _S309[int(2)];
            points_4[int(3)] = _S310[int(3)];
            raydir_21 = _S236;
        }
        float3  _S311 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S312 = _S311.x + _S311.y + _S311.z;
        float4  _S313 = _S269;
        *&((&_S313)->x) = _S312;
        _S283 = _S283 + _S313;
    }
    else
    {
        _S283 = _S269;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S283;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S314, float4  _S315, FixedArray<float, 10>  * _S316, bool _S317, bool _S318, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S319, float3  _S320)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S321;
    float3  _S322 = s_primal_ctx_depth_to_normal_0(_S314, _S315, _S316, _S317, _S318, (*_S319).primal_0, &_S321);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S323 = _S321;
    s_bwd_prop_depth_to_normal_0(_S314, _S315, _S316, _S317, _S318, _S319, _S320, &_S323);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S324 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S324;
    FixedArray<float, 10>  _S325 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S325, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S326 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S326;
    FixedArray<float, 10>  _S327 = dist_coeffs_6;
    bool _S328 = undistort_point_0(_S326, &_S327, int(12), &uv_10);
    if(!_S328)
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

