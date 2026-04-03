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

inline __device__ float srgb_to_linear_rgb(float x_5)
{
    float _S50;
    if(x_5 < 0.05499999970197678f)
    {
        _S50 = x_5 * 0.07739938050508499f;
    }
    else
    {
        _S50 = (F32_pow((0.94786733388900757f * (x_5 + 0.05499999970197678f)), (2.40000009536743164f)));
    }
    return _S50;
}

inline __device__ float srgb_to_linear_rgb_grad(float x_6)
{
    float _S51;
    if(x_6 < 0.05499999970197678f)
    {
        _S51 = 0.07739938050508499f;
    }
    else
    {
        DiffPair_float_0 _S52;
        (&_S52)->primal_0 = 0.94786733388900757f * (x_6 + 0.05499999970197678f);
        (&_S52)->differential_0 = 0.94786733388900757f;
        DiffPair_float_0 _S53;
        (&_S53)->primal_0 = 2.40000009536743164f;
        (&_S53)->differential_0 = 0.0f;
        DiffPair_float_0 _S54 = _d_pow_1(&_S52, &_S53);
        _S51 = _S54.differential_0;
    }
    return _S51;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_4)
{
    float _S55 = (*left_0).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_4.x;
    float sum_0 = _S55 + (*left_0).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_4.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_4.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_1;
    float _S56 = (*left_0).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_4.x;
    float sum_2 = _S56 + (*left_0).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_4.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_4.z;
    *&((&right_d_result_1)->y) = sum_3;
    float _S57 = (*left_0).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_4.x;
    float sum_4 = _S57 + (*left_0).primal_0.rows[int(1)].z * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_4.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_4.z;
    *&((&right_d_result_1)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_1;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_1;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_2;
    int i_2 = int(0);
    for(;;)
    {
        if(i_2 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_6 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_2], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_2, i_2) = sum_6;
        i_2 = i_2 + int(1);
    }
    return result_2;
}

inline __device__ float3  linear_rgb_to_srgb(float3  rgb_2, Matrix<float, 3, 3>  color_matrix_0)
{
    float3  _S58 = mul_0(color_matrix_0, rgb_2);
    float _S59 = _S58.x;
    float _S60;
    if(_S59 < 0.00313080009073019f)
    {
        _S60 = _S59 * 12.92000007629394531f;
    }
    else
    {
        _S60 = 1.0549999475479126f * (F32_pow((_S59), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S61 = _S58.y;
    float _S62;
    if(_S61 < 0.00313080009073019f)
    {
        _S62 = _S61 * 12.92000007629394531f;
    }
    else
    {
        _S62 = 1.0549999475479126f * (F32_pow((_S61), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S63 = _S58.z;
    float _S64;
    if(_S63 < 0.00313080009073019f)
    {
        _S64 = _S63 * 12.92000007629394531f;
    }
    else
    {
        _S64 = 1.0549999475479126f * (F32_pow((_S63), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S60, _S62, _S64);
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S65, float3  _S66)
{
    return mul_0(_S65, _S66);
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S67, DiffPair_float_0 * _S68, float _S69)
{
    _d_pow_0(_S67, _S68, _S69);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S70, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S71, float3  _S72)
{
    _d_mul_0(_S70, _S71, _S72);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, Matrix<float, 3, 3>  color_matrix_1, float3  _s_dOut_2)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S73 = *dprgb_1;
    float3  _S74 = s_primal_ctx_mul_0(color_matrix_1, (*dprgb_1).primal_0);
    float _S75 = _S74.x;
    float _S76 = _S74.y;
    float _S77 = _S74.z;
    float _S78;
    if(_S77 < 0.00313080009073019f)
    {
        _S78 = 12.92000007629394531f * _s_dOut_2.z;
    }
    else
    {
        float _S79 = 1.0549999475479126f * _s_dOut_2.z;
        DiffPair_float_0 _S80;
        (&_S80)->primal_0 = _S77;
        (&_S80)->differential_0 = 0.0f;
        DiffPair_float_0 _S81;
        (&_S81)->primal_0 = 0.4166666567325592f;
        (&_S81)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S80, &_S81, _S79);
        _S78 = _S80.differential_0;
    }
    float _S82;
    if(_S76 < 0.00313080009073019f)
    {
        _S82 = 12.92000007629394531f * _s_dOut_2.y;
    }
    else
    {
        float _S83 = 1.0549999475479126f * _s_dOut_2.y;
        DiffPair_float_0 _S84;
        (&_S84)->primal_0 = _S76;
        (&_S84)->differential_0 = 0.0f;
        DiffPair_float_0 _S85;
        (&_S85)->primal_0 = 0.4166666567325592f;
        (&_S85)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S84, &_S85, _S83);
        _S82 = _S84.differential_0;
    }
    float _S86;
    if(_S75 < 0.00313080009073019f)
    {
        _S86 = 12.92000007629394531f * _s_dOut_2.x;
    }
    else
    {
        float _S87 = 1.0549999475479126f * _s_dOut_2.x;
        DiffPair_float_0 _S88;
        (&_S88)->primal_0 = _S75;
        (&_S88)->differential_0 = 0.0f;
        DiffPair_float_0 _S89;
        (&_S89)->primal_0 = 0.4166666567325592f;
        (&_S89)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S88, &_S89, _S87);
        _S86 = _S88.differential_0;
    }
    float3  _S90 = make_float3 (_S86, _S82, _S78);
    Matrix<float, 3, 3>  _S91 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S92;
    (&_S92)->primal_0 = color_matrix_1;
    (&_S92)->differential_0 = _S91;
    float3  _S93 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S94;
    (&_S94)->primal_0 = _S73.primal_0;
    (&_S94)->differential_0 = _S93;
    s_bwd_prop_mul_0(&_S92, &_S94, _S90);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S94.differential_0;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S95, Matrix<float, 3, 3>  _S96, float3  _S97)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S95, _S96, _S97);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, Matrix<float, 3, 3>  color_matrix_2, float3  v_out_rgb_1)
{
    float3  _S98 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S98;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, color_matrix_2, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ float3  rgb_to_srgb(float3  rgb_4, Matrix<float, 3, 3>  color_matrix_3)
{
    float _S99 = rgb_4.x;
    float _S100;
    if(_S99 < 0.05499999970197678f)
    {
        _S100 = _S99 * 0.07739938050508499f;
    }
    else
    {
        _S100 = (F32_pow((0.94786733388900757f * (_S99 + 0.05499999970197678f)), (2.40000009536743164f)));
    }
    float _S101 = rgb_4.y;
    float _S102;
    if(_S101 < 0.05499999970197678f)
    {
        _S102 = _S101 * 0.07739938050508499f;
    }
    else
    {
        _S102 = (F32_pow((0.94786733388900757f * (_S101 + 0.05499999970197678f)), (2.40000009536743164f)));
    }
    float _S103 = rgb_4.z;
    float _S104;
    if(_S103 < 0.05499999970197678f)
    {
        _S104 = _S103 * 0.07739938050508499f;
    }
    else
    {
        _S104 = (F32_pow((0.94786733388900757f * (_S103 + 0.05499999970197678f)), (2.40000009536743164f)));
    }
    float3  _S105 = mul_0(color_matrix_3, make_float3 (_S100, _S102, _S104));
    float _S106 = _S105.x;
    if(_S106 < 0.00313080009073019f)
    {
        _S100 = _S106 * 12.92000007629394531f;
    }
    else
    {
        _S100 = 1.0549999475479126f * (F32_pow((_S106), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S107 = _S105.y;
    if(_S107 < 0.00313080009073019f)
    {
        _S102 = _S107 * 12.92000007629394531f;
    }
    else
    {
        _S102 = 1.0549999475479126f * (F32_pow((_S107), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S108 = _S105.z;
    if(_S108 < 0.00313080009073019f)
    {
        _S104 = _S108 * 12.92000007629394531f;
    }
    else
    {
        _S104 = 1.0549999475479126f * (F32_pow((_S108), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S100, _S102, _S104);
}

inline __device__ float s_primal_ctx_pow_0(float _S109, float _S110)
{
    return (F32_pow((_S109), (_S110)));
}

inline __device__ void s_bwd_prop_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_2, Matrix<float, 3, 3>  color_matrix_4, float3  _s_dOut_3)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S111 = *dprgb_2;
    float _S112 = (*dprgb_2).primal_0.x;
    bool _S113 = _S112 < 0.05499999970197678f;
    float _S114;
    if(_S113)
    {
        _S114 = _S112 * 0.07739938050508499f;
    }
    else
    {
        _S114 = s_primal_ctx_pow_0(0.94786733388900757f * (_S112 + 0.05499999970197678f), 2.40000009536743164f);
    }
    float _S115 = _S111.primal_0.y;
    bool _S116 = _S115 < 0.05499999970197678f;
    float _S117;
    if(_S116)
    {
        _S117 = _S115 * 0.07739938050508499f;
    }
    else
    {
        _S117 = s_primal_ctx_pow_0(0.94786733388900757f * (_S115 + 0.05499999970197678f), 2.40000009536743164f);
    }
    float _S118 = _S111.primal_0.z;
    bool _S119 = _S118 < 0.05499999970197678f;
    float _S120;
    if(_S119)
    {
        _S120 = _S118 * 0.07739938050508499f;
    }
    else
    {
        _S120 = s_primal_ctx_pow_0(0.94786733388900757f * (_S118 + 0.05499999970197678f), 2.40000009536743164f);
    }
    float3  _S121 = make_float3 (_S114, _S117, _S120);
    float3  _S122 = s_primal_ctx_mul_0(color_matrix_4, _S121);
    float _S123 = _S122.x;
    float _S124 = _S122.y;
    float _S125 = _S122.z;
    if(_S125 < 0.00313080009073019f)
    {
        _S114 = 12.92000007629394531f * _s_dOut_3.z;
    }
    else
    {
        float _S126 = 1.0549999475479126f * _s_dOut_3.z;
        DiffPair_float_0 _S127;
        (&_S127)->primal_0 = _S125;
        (&_S127)->differential_0 = 0.0f;
        DiffPair_float_0 _S128;
        (&_S128)->primal_0 = 0.4166666567325592f;
        (&_S128)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S127, &_S128, _S126);
        _S114 = _S127.differential_0;
    }
    if(_S124 < 0.00313080009073019f)
    {
        _S117 = 12.92000007629394531f * _s_dOut_3.y;
    }
    else
    {
        float _S129 = 1.0549999475479126f * _s_dOut_3.y;
        DiffPair_float_0 _S130;
        (&_S130)->primal_0 = _S124;
        (&_S130)->differential_0 = 0.0f;
        DiffPair_float_0 _S131;
        (&_S131)->primal_0 = 0.4166666567325592f;
        (&_S131)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S130, &_S131, _S129);
        _S117 = _S130.differential_0;
    }
    if(_S123 < 0.00313080009073019f)
    {
        _S120 = 12.92000007629394531f * _s_dOut_3.x;
    }
    else
    {
        float _S132 = 1.0549999475479126f * _s_dOut_3.x;
        DiffPair_float_0 _S133;
        (&_S133)->primal_0 = _S123;
        (&_S133)->differential_0 = 0.0f;
        DiffPair_float_0 _S134;
        (&_S134)->primal_0 = 0.4166666567325592f;
        (&_S134)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S133, &_S134, _S132);
        _S120 = _S133.differential_0;
    }
    float3  _S135 = make_float3 (_S120, _S117, _S114);
    Matrix<float, 3, 3>  _S136 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S137;
    (&_S137)->primal_0 = color_matrix_4;
    (&_S137)->differential_0 = _S136;
    float3  _S138 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S139;
    (&_S139)->primal_0 = _S121;
    (&_S139)->differential_0 = _S138;
    s_bwd_prop_mul_0(&_S137, &_S139, _S135);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S140 = _S139;
    if(_S119)
    {
        _S114 = 0.0f;
    }
    else
    {
        _S114 = 0.94786733388900757f * (_S118 + 0.05499999970197678f);
    }
    if(_S119)
    {
        _S114 = 0.07739938050508499f * _S140.differential_0.z;
    }
    else
    {
        DiffPair_float_0 _S141;
        (&_S141)->primal_0 = _S114;
        (&_S141)->differential_0 = 0.0f;
        DiffPair_float_0 _S142;
        (&_S142)->primal_0 = 2.40000009536743164f;
        (&_S142)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S141, &_S142, _S140.differential_0.z);
        _S114 = 0.94786733388900757f * _S141.differential_0;
    }
    if(_S116)
    {
        _S117 = 0.0f;
    }
    else
    {
        _S117 = 0.94786733388900757f * (_S115 + 0.05499999970197678f);
    }
    if(_S116)
    {
        _S117 = 0.07739938050508499f * _S140.differential_0.y;
    }
    else
    {
        DiffPair_float_0 _S143;
        (&_S143)->primal_0 = _S117;
        (&_S143)->differential_0 = 0.0f;
        DiffPair_float_0 _S144;
        (&_S144)->primal_0 = 2.40000009536743164f;
        (&_S144)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S143, &_S144, _S140.differential_0.y);
        _S117 = 0.94786733388900757f * _S143.differential_0;
    }
    if(_S113)
    {
        _S120 = 0.0f;
    }
    else
    {
        _S120 = 0.94786733388900757f * (_S112 + 0.05499999970197678f);
    }
    if(_S113)
    {
        _S120 = 0.07739938050508499f * _S140.differential_0.x;
    }
    else
    {
        DiffPair_float_0 _S145;
        (&_S145)->primal_0 = _S120;
        (&_S145)->differential_0 = 0.0f;
        DiffPair_float_0 _S146;
        (&_S146)->primal_0 = 2.40000009536743164f;
        (&_S146)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S145, &_S146, _S140.differential_0.x);
        _S120 = 0.94786733388900757f * _S145.differential_0;
    }
    float3  _S147 = make_float3 (_S120, _S117, _S114);
    dprgb_2->primal_0 = (*dprgb_2).primal_0;
    dprgb_2->differential_0 = _S147;
    return;
}

inline __device__ void s_bwd_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S148, Matrix<float, 3, 3>  _S149, float3  _S150)
{
    s_bwd_prop_rgb_to_srgb_0(_S148, _S149, _S150);
    return;
}

inline __device__ float3  rgb_to_srgb_bwd(float3  rgb_5, Matrix<float, 3, 3>  color_matrix_5, float3  v_out_rgb_2)
{
    float3  _S151 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_2;
    (&p_rgb_2)->primal_0 = rgb_5;
    (&p_rgb_2)->differential_0 = _S151;
    s_bwd_rgb_to_srgb_0(&p_rgb_2, color_matrix_5, v_out_rgb_2);
    return p_rgb_2.differential_0;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_7)
{
    Matrix<float, 2, 2>  result_3;
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
            *_slang_vector_get_element_ptr(((&result_3)->rows + (r_0)), c_0) = _slang_vector_get_element(x_7.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_3;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_5, float dOut_5)
{
    float _S152 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_5).primal_0)))))) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S152;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, float dOut_6)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_4).primal_0.x * dOut_6;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_6).primal_0.x * dOut_6;
    *&((&x_d_result_0)->y) = (*dpy_4).primal_0.y * dOut_6;
    *&((&y_d_result_0)->y) = (*dpx_6).primal_0.y * dOut_6;
    *&((&x_d_result_0)->z) = (*dpy_4).primal_0.z * dOut_6;
    *&((&y_d_result_0)->z) = (*dpx_6).primal_0.z * dOut_6;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = x_d_result_0;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_8, float3  y_2)
{
    int i_3 = int(0);
    float result_4 = 0.0f;
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_5 = result_4 + _slang_vector_get_element(x_8, i_3) * _slang_vector_get_element(y_2, i_3);
        i_3 = i_3 + int(1);
        result_4 = result_5;
    }
    return result_4;
}

inline __device__ float dot_1(float2  x_9, float2  y_3)
{
    int i_4 = int(0);
    float result_6 = 0.0f;
    for(;;)
    {
        if(i_4 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_7 = result_6 + _slang_vector_get_element(x_9, i_4) * _slang_vector_get_element(y_3, i_4);
        i_4 = i_4 + int(1);
        result_6 = result_7;
    }
    return result_6;
}

inline __device__ float length_0(float2  x_10)
{
    return (F32_sqrt((dot_1(x_10, x_10))));
}

inline __device__ float length_1(float3  x_11)
{
    return (F32_sqrt((dot_0(x_11, x_11))));
}

inline __device__ bool undistort_point_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0, int maxiter_0, float2  * uv_undist_0)
{
    int i_5 = int(0);
    float2  q_0 = uv_0;
    for(;;)
    {
        if(i_5 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float _S153 = (*dist_coeffs_0)[int(3)];
        float _S154 = (*dist_coeffs_0)[int(4)];
        float _S155 = (*dist_coeffs_0)[int(5)];
        float _S156 = (*dist_coeffs_0)[int(6)];
        float _S157 = (*dist_coeffs_0)[int(7)];
        float _S158 = (*dist_coeffs_0)[int(8)];
        float _S159 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S160 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S161 = (*dist_coeffs_0)[int(1)] + r2_0 * _S160;
        float _S162 = (*dist_coeffs_0)[int(0)] + r2_0 * _S161;
        float radial_0 = 1.0f + r2_0 * _S162;
        float _S163 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S164 = _S163 * u_0;
        float _S165 = 2.0f * u_0;
        float _S166 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S167 = _S166 * u_0;
        float _S168 = 2.0f * v_0;
        float2  _S169 = q_0 * make_float2 (radial_0) + make_float2 (_S164 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S165 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S167 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S168 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S169 + make_float2 ((*dist_coeffs_0)[int(8)] * _S169.x + (*dist_coeffs_0)[int(9)] * _S169.y, 0.0f) - uv_0;
        float _S170 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S170 + _S170);
        float2  _S171 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S162 + (s_diff_r2_0 * _S161 + (s_diff_r2_0 * _S160 + s_diff_r2_0 * _S153 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S163 * v_0 + 0.0f * _S164 + (s_diff_r2_0 + (_S165 + _S165)) * _S155 + s_diff_r2_0 * _S156, _S166 * v_0 + 0.0f * _S167 + (s_diff_r2_0 + (_S170 + 0.0f * _S168)) * _S154 + s_diff_r2_0 * _S157);
        float _S172 = 0.0f * u_0;
        float s_diff_r2_1 = _S172 + _S172 + (v_0 + v_0);
        float2  _S173 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S162 + (s_diff_r2_1 * _S161 + (s_diff_r2_1 * _S160 + s_diff_r2_1 * _S153 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S163 * v_0 + _S164 + (s_diff_r2_1 + (_S172 + 0.0f * _S165)) * _S155 + s_diff_r2_1 * _S156, 0.0f * _S166 * v_0 + _S167 + (s_diff_r2_1 + (_S168 + _S168)) * _S154 + s_diff_r2_1 * _S157);
        Matrix<float, 2, 2>  _S174 = transpose_0(makeMatrix<float, 2, 2> (_S171 + make_float2 (_S171.x * _S158 + _S171.y * _S159, 0.0f), _S173 + make_float2 (_S173.x * _S158 + _S173.y * _S159, 0.0f)));
        float inv_det_0 = 1.0f / (_S174.rows[int(0)].x * _S174.rows[int(1)].y - _S174.rows[int(0)].y * _S174.rows[int(1)].x);
        float _S175 = r_1.x;
        float _S176 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S175 * _S174.rows[int(1)].y - _S176 * _S174.rows[int(0)].y) * inv_det_0, (- _S175 * _S174.rows[int(1)].x + _S176 * _S174.rows[int(0)].x) * inv_det_0);
        i_5 = i_5 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S177 = (*dist_coeffs_0)[int(0)];
    float _S178 = (*dist_coeffs_0)[int(1)];
    float _S179 = (*dist_coeffs_0)[int(2)];
    float _S180 = (*dist_coeffs_0)[int(3)];
    float _S181 = (*dist_coeffs_0)[int(4)];
    float _S182 = (*dist_coeffs_0)[int(5)];
    float _S183 = (*dist_coeffs_0)[int(6)];
    float _S184 = (*dist_coeffs_0)[int(7)];
    float _S185 = (*dist_coeffs_0)[int(8)];
    float _S186 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S187 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S187 + _S187);
    float _S188 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S189 = (*dist_coeffs_0)[int(1)] + r2_1 * _S188;
    float _S190 = (*dist_coeffs_0)[int(0)] + r2_1 * _S189;
    float radial_1 = 1.0f + r2_1 * _S190;
    float _S191 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S192 = _S191 * u_1;
    float _S193 = 2.0f * u_1;
    float _S194 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S195 = _S194 * u_1;
    float _S196 = 2.0f * v_1;
    float2  _S197 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S190 + (s_diff_r2_2 * _S189 + (s_diff_r2_2 * _S188 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S191 * v_1 + 0.0f * _S192 + (s_diff_r2_2 + (_S193 + _S193)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S194 * v_1 + 0.0f * _S195 + (s_diff_r2_2 + (_S187 + 0.0f * _S196)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S198 = 0.0f * u_1;
    float s_diff_r2_3 = _S198 + _S198 + (v_1 + v_1);
    float2  _S199 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S190 + (s_diff_r2_3 * _S189 + (s_diff_r2_3 * _S188 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S191 * v_1 + _S192 + (s_diff_r2_3 + (_S198 + 0.0f * _S193)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S194 * v_1 + _S195 + (s_diff_r2_3 + (_S196 + _S196)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S200 = transpose_0(makeMatrix<float, 2, 2> (_S197 + make_float2 (_S197.x * (*dist_coeffs_0)[int(8)] + _S197.y * (*dist_coeffs_0)[int(9)], 0.0f), _S199 + make_float2 (_S199.x * (*dist_coeffs_0)[int(8)] + _S199.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S201;
    if((F32_min((determinant_0(_S200)), ((F32_min((_S200.rows[int(0)].x), (_S200.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S202 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S177 + r2_2 * (_S178 + r2_2 * (_S179 + r2_2 * _S180)))) + make_float2 (_S191 * u_2 * v_2 + _S182 * (r2_2 + 2.0f * u_2 * u_2) + _S183 * r2_2, _S194 * u_2 * v_2 + _S181 * (r2_2 + 2.0f * v_2 * v_2) + _S184 * r2_2);
        _S201 = (length_0(_S202 + make_float2 (_S185 * _S202.x + _S186 * _S202.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S201 = false;
    }
    return _S201;
}

inline __device__ float3  normalize_0(float3  x_12)
{
    return x_12 / make_float3 (length_1(x_12));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S203 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S203;
    FixedArray<float, 10>  _S204 = dist_coeffs_1;
    bool _S205 = undistort_point_0(_S203, &_S204, int(12), &uv_1);
    if(!_S205)
    {
        int3  _S206 = make_int3 (int(0));
        float3  _S207 = make_float3 ((float)_S206.x, (float)_S206.y, (float)_S206.z);
        return _S207;
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

inline __device__ float3  depth_to_point(float2  pix_pos_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float depth_2)
{
    float2  _S208 = (pix_pos_1 - float2 {intrins_1.z, intrins_1.w}) / float2 {intrins_1.x, intrins_1.y};
    float2  uv_2 = _S208;
    FixedArray<float, 10>  _S209 = dist_coeffs_2;
    bool _S210 = undistort_point_0(_S208, &_S209, int(12), &uv_2);
    if(!_S210)
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
    return make_float3 (depth_2) * raydir_3;
}

struct s_bwd_prop_depth_to_point_Intermediates_0
{
    float2  _S211;
    bool _S212;
};

inline __device__ float s_primal_ctx_sin_0(float _S213)
{
    return (F32_sin((_S213)));
}

inline __device__ float s_primal_ctx_cos_0(float _S214)
{
    return (F32_cos((_S214)));
}

inline __device__ float3  s_primal_ctx_depth_to_point_0(float2  pix_pos_2, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float dpdepth_1, s_bwd_prop_depth_to_point_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S211 = make_float2 (0.0f);
    _s_diff_ctx_0->_S212 = false;
    float2  _S215 = (pix_pos_2 - float2 {intrins_2.z, intrins_2.w}) / float2 {intrins_2.x, intrins_2.y};
    float2  _S216 = _S215;
    bool _S217 = undistort_point_0(_S215, dist_coeffs_3, int(12), &_S216);
    _s_diff_ctx_0->_S211 = _S216;
    _s_diff_ctx_0->_S212 = _S217;
    float2  uv_3 = _S216;
    bool _S218 = !_S217;
    float3  raydir_6;
    if(_S218)
    {
        raydir_6 = make_float3 (0.0f);
    }
    bool _S219 = !_S218;
    if(_S219)
    {
        if(is_fisheye_2)
        {
            float _S220 = length_0(uv_3);
            float3  raydir_7 = make_float3 ((uv_3 / make_float2 ((F32_max((_S220), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S220))).x, (uv_3 / make_float2 ((F32_max((_S220), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S220))).y, s_primal_ctx_cos_0(_S220));
            if(!is_ray_depth_2)
            {
                raydir_6 = raydir_7 / make_float3 (raydir_7.z);
            }
            else
            {
                raydir_6 = raydir_7;
            }
        }
        else
        {
            float3  raydir_8 = make_float3 (uv_3.x, uv_3.y, 1.0f);
            if(is_ray_depth_2)
            {
                raydir_6 = normalize_0(raydir_8);
            }
            else
            {
                raydir_6 = raydir_8;
            }
        }
        raydir_6 = make_float3 (dpdepth_1) * raydir_6;
    }
    return raydir_6;
}

inline __device__ void s_bwd_prop_depth_to_point_0(float2  pix_pos_3, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_float_0 * dpdepth_2, float3  _s_dOut_4, s_bwd_prop_depth_to_point_Intermediates_0 * _s_diff_ctx_1)
{
    float3  _S221 = make_float3 (0.0f);
    float2  _S222 = _s_diff_ctx_1->_S211;
    bool _S223 = !!_s_diff_ctx_1->_S212;
    float3  raydir_9;
    if(_S223)
    {
        if(is_fisheye_3)
        {
            float _S224 = length_0(_S222);
            float3  raydir_10 = make_float3 ((_S222 / make_float2 ((F32_max((_S224), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S224))).x, (_S222 / make_float2 ((F32_max((_S224), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S224))).y, s_primal_ctx_cos_0(_S224));
            if(!is_ray_depth_3)
            {
                raydir_9 = raydir_10 / make_float3 (raydir_10.z);
            }
            else
            {
                raydir_9 = raydir_10;
            }
        }
        else
        {
            float3  raydir_11 = make_float3 (_S222.x, _S222.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_9 = normalize_0(raydir_11);
            }
            else
            {
                raydir_9 = raydir_11;
            }
        }
    }
    else
    {
        raydir_9 = _S221;
    }
    if(_S223)
    {
        raydir_9 = raydir_9 * _s_dOut_4;
    }
    else
    {
        raydir_9 = _S221;
    }
    float _S225 = raydir_9.x + raydir_9.y + raydir_9.z;
    dpdepth_2->primal_0 = (*dpdepth_2).primal_0;
    dpdepth_2->differential_0 = _S225;
    return;
}

inline __device__ void s_bwd_depth_to_point_0(float2  _S226, float4  _S227, FixedArray<float, 10>  * _S228, bool _S229, bool _S230, DiffPair_float_0 * _S231, float3  _S232)
{
    s_bwd_prop_depth_to_point_Intermediates_0 _S233;
    float3  _S234 = s_primal_ctx_depth_to_point_0(_S226, _S227, _S228, _S229, _S230, (*_S231).primal_0, &_S233);
    s_bwd_prop_depth_to_point_Intermediates_0 _S235 = _S233;
    s_bwd_prop_depth_to_point_0(_S226, _S227, _S228, _S229, _S230, _S231, _S232, &_S235);
    return;
}

inline __device__ float depth_to_point_vjp(float2  pix_pos_4, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float depth_3, float3  v_point_0)
{
    DiffPair_float_0 dp_depth_0;
    (&dp_depth_0)->primal_0 = depth_3;
    (&dp_depth_0)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S236 = dist_coeffs_5;
    s_bwd_depth_to_point_0(pix_pos_4, intrins_4, &_S236, is_fisheye_4, is_ray_depth_4, &dp_depth_0, v_point_0);
    return dp_depth_0.differential_0;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_7)
{
    float _S237 = dOut_7.y;
    float _S238 = dOut_7.z;
    float _S239 = dOut_7.x;
    float _S240 = (*a_0).primal_0.z * _S237 + - (*a_0).primal_0.y * _S238;
    float _S241 = - (*a_0).primal_0.z * _S239 + (*a_0).primal_0.x * _S238;
    float _S242 = (*a_0).primal_0.y * _S239 + - (*a_0).primal_0.x * _S237;
    float3  _S243 = make_float3 (- (*b_0).primal_0.z * _S237 + (*b_0).primal_0.y * _S238, (*b_0).primal_0.z * _S239 + - (*b_0).primal_0.x * _S238, - (*b_0).primal_0.y * _S239 + (*b_0).primal_0.x * _S237);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S243;
    float3  _S244 = make_float3 (_S240, _S241, _S242);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S244;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S245 = left_2.y;
    float _S246 = right_2.z;
    float _S247 = left_2.z;
    float _S248 = right_2.y;
    float _S249 = right_2.x;
    float _S250 = left_2.x;
    return make_float3 (_S245 * _S246 - _S247 * _S248, _S247 * _S249 - _S250 * _S246, _S250 * _S248 - _S245 * _S249);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S251, float3  _S252)
{
    return cross_0(_S251, _S252);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S253, float3  _S254)
{
    return dot_0(_S253, _S254);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S255, float _S256)
{
    _d_sqrt_0(_S255, _S256);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_5)
{
    float _S257 = (*dpx_7).primal_0.x;
    float _S258 = (*dpx_7).primal_0.y;
    float _S259 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S260;
    (&_S260)->primal_0 = _S257 * _S257 + _S258 * _S258 + _S259 * _S259;
    (&_S260)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S260, _s_dOut_5);
    float _S261 = (*dpx_7).primal_0.z * _S260.differential_0;
    float _S262 = _S261 + _S261;
    float _S263 = (*dpx_7).primal_0.y * _S260.differential_0;
    float _S264 = _S263 + _S263;
    float _S265 = (*dpx_7).primal_0.x * _S260.differential_0;
    float _S266 = _S265 + _S265;
    float3  _S267 = make_float3 (0.0f);
    *&((&_S267)->z) = _S262;
    *&((&_S267)->y) = _S264;
    *&((&_S267)->x) = _S266;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S267;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S268, float _S269)
{
    s_bwd_prop_length_impl_0(_S268, _S269);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S270, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S271, float _S272)
{
    _d_dot_0(_S270, _S271, _S272);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S273, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S274, float3  _S275)
{
    _d_cross_0(_S273, _S274, _S275);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_6)
{
    float3  _S276 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S277 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S278 = s_primal_ctx_cross_0(dx_0, _S277);
    bool _S279 = (s_primal_ctx_dot_0(_S278, _S278)) != 0.0f;
    float3  _S280;
    float3  _S281;
    if(_S279)
    {
        float _S282 = length_1(_S278);
        float3  _S283 = make_float3 (_S282);
        _S280 = make_float3 (_S282 * _S282);
        _S281 = _S283;
    }
    else
    {
        _S280 = _S276;
        _S281 = _S276;
    }
    if(_S279)
    {
        float3  _S284 = _s_dOut_6 / _S280;
        float3  _S285 = _S278 * - _S284;
        float3  _S286 = _S281 * _S284;
        float _S287 = _S285.x + _S285.y + _S285.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S288;
        (&_S288)->primal_0 = _S278;
        (&_S288)->differential_0 = _S276;
        s_bwd_length_impl_0(&_S288, _S287);
        _S280 = _S286 + _S288.differential_0;
    }
    else
    {
        _S280 = _s_dOut_6;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S289;
    (&_S289)->primal_0 = _S278;
    (&_S289)->differential_0 = _S276;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S290;
    (&_S290)->primal_0 = _S278;
    (&_S290)->differential_0 = _S276;
    s_bwd_prop_dot_0(&_S289, &_S290, 0.0f);
    float3  _S291 = _S290.differential_0 + _S289.differential_0 + _S280;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S292;
    (&_S292)->primal_0 = dx_0;
    (&_S292)->differential_0 = _S276;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S293;
    (&_S293)->primal_0 = _S277;
    (&_S293)->differential_0 = _S276;
    s_bwd_prop_cross_0(&_S292, &_S293, _S291);
    float3  s_diff_dy_T_0 = - _S293.differential_0;
    float3  _S294 = - s_diff_dy_T_0;
    float3  _S295 = - _S292.differential_0;
    FixedArray<float3 , 4>  _S296;
    _S296[int(0)] = _S276;
    _S296[int(1)] = _S276;
    _S296[int(2)] = _S276;
    _S296[int(3)] = _S276;
    _S296[int(2)] = _S294;
    _S296[int(3)] = s_diff_dy_T_0;
    _S296[int(0)] = _S295;
    _S296[int(1)] = _S292.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S296;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S297, float3  _S298)
{
    s_bwd_prop_points_to_normal_0(_S297, _S298);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S299 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S299;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5, bool is_ray_depth_5, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S300 = float2 {intrins_5.z, intrins_5.w};
    float2  _S301 = float2 {intrins_5.x, intrins_5.y};
    float2  _S302 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S300) / _S301;
    float2  uv_4 = _S302;
    FixedArray<float, 10>  _S303 = dist_coeffs_6;
    bool _S304 = undistort_point_0(_S302, &_S303, int(12), &uv_4);
    if(!_S304)
    {
        return make_float3 (0.0f);
    }
    float3  raydir_12;
    if(is_fisheye_5)
    {
        float theta_2 = length_0(uv_4);
        float3  raydir_13 = make_float3 ((uv_4 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).x, (uv_4 / make_float2 ((F32_max((theta_2), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_2))))).y, (F32_cos((theta_2))));
        if(!is_ray_depth_5)
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
        float3  raydir_14 = make_float3 (uv_4.x, uv_4.y, 1.0f);
        if(is_ray_depth_5)
        {
            raydir_12 = normalize_0(raydir_14);
        }
        else
        {
            raydir_12 = raydir_14;
        }
    }
    points_2[int(0)] = make_float3 (depths_0.x) * raydir_12;
    float2  _S305 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S300) / _S301;
    float2  uv_5 = _S305;
    FixedArray<float, 10>  _S306 = dist_coeffs_6;
    bool _S307 = undistort_point_0(_S305, &_S306, int(12), &uv_5);
    if(!_S307)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_3 = length_0(uv_5);
        float3  raydir_15 = make_float3 ((uv_5 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).x, (uv_5 / make_float2 ((F32_max((theta_3), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_3))))).y, (F32_cos((theta_3))));
        if(!is_ray_depth_5)
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
        float3  raydir_16 = make_float3 (uv_5.x, uv_5.y, 1.0f);
        if(is_ray_depth_5)
        {
            raydir_12 = normalize_0(raydir_16);
        }
        else
        {
            raydir_12 = raydir_16;
        }
    }
    points_2[int(1)] = make_float3 (depths_0.y) * raydir_12;
    float2  _S308 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S300) / _S301;
    float2  uv_6 = _S308;
    FixedArray<float, 10>  _S309 = dist_coeffs_6;
    bool _S310 = undistort_point_0(_S308, &_S309, int(12), &uv_6);
    if(!_S310)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_4 = length_0(uv_6);
        float3  raydir_17 = make_float3 ((uv_6 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).x, (uv_6 / make_float2 ((F32_max((theta_4), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_4))))).y, (F32_cos((theta_4))));
        if(!is_ray_depth_5)
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
        float3  raydir_18 = make_float3 (uv_6.x, uv_6.y, 1.0f);
        if(is_ray_depth_5)
        {
            raydir_12 = normalize_0(raydir_18);
        }
        else
        {
            raydir_12 = raydir_18;
        }
    }
    points_2[int(2)] = make_float3 (depths_0.z) * raydir_12;
    float2  _S311 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S300) / _S301;
    float2  uv_7 = _S311;
    FixedArray<float, 10>  _S312 = dist_coeffs_6;
    bool _S313 = undistort_point_0(_S311, &_S312, int(12), &uv_7);
    if(!_S313)
    {
        return make_float3 (0.0f);
    }
    if(is_fisheye_5)
    {
        float theta_5 = length_0(uv_7);
        float3  raydir_19 = make_float3 ((uv_7 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).x, (uv_7 / make_float2 ((F32_max((theta_5), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_5))))).y, (F32_cos((theta_5))));
        if(!is_ray_depth_5)
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
        float3  raydir_20 = make_float3 (uv_7.x, uv_7.y, 1.0f);
        if(is_ray_depth_5)
        {
            raydir_12 = normalize_0(raydir_20);
        }
        else
        {
            raydir_12 = raydir_20;
        }
    }
    points_2[int(3)] = make_float3 (depths_0.w) * raydir_12;
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
    float2  _S314;
    bool _S315;
    float2  _S316;
    bool _S317;
    float2  _S318;
    bool _S319;
    float2  _S320;
    bool _S321;
};

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_6, FixedArray<float, 10>  * dist_coeffs_7, bool is_fisheye_6, bool is_ray_depth_6, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S322 = make_float2 (0.0f);
    _s_diff_ctx_2->_S314 = _S322;
    _s_diff_ctx_2->_S315 = false;
    _s_diff_ctx_2->_S316 = _S322;
    _s_diff_ctx_2->_S317 = false;
    _s_diff_ctx_2->_S318 = _S322;
    _s_diff_ctx_2->_S319 = false;
    _s_diff_ctx_2->_S320 = _S322;
    _s_diff_ctx_2->_S321 = false;
    _s_diff_ctx_2->_S316 = _S322;
    _s_diff_ctx_2->_S317 = false;
    _s_diff_ctx_2->_S318 = _S322;
    _s_diff_ctx_2->_S319 = false;
    _s_diff_ctx_2->_S320 = _S322;
    _s_diff_ctx_2->_S321 = false;
    float3  _S323 = make_float3 (0.0f);
    float2  _S324 = float2 {intrins_6.z, intrins_6.w};
    float2  _S325 = float2 {intrins_6.x, intrins_6.y};
    float2  _S326 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S324) / _S325;
    float2  _S327 = _S326;
    bool _S328 = undistort_point_0(_S326, dist_coeffs_7, int(12), &_S327);
    _s_diff_ctx_2->_S314 = _S327;
    _s_diff_ctx_2->_S315 = _S328;
    float2  uv_8 = _S327;
    bool _S329 = !_S328;
    float3  normal_4;
    if(_S329)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S330 = !_S329;
    int _S331;
    FixedArray<float3 , 4>  points_3;
    if(_S330)
    {
        float3  raydir_21;
        if(is_fisheye_6)
        {
            float _S332 = length_0(uv_8);
            float3  raydir_22 = make_float3 ((uv_8 / make_float2 ((F32_max((_S332), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S332))).x, (uv_8 / make_float2 ((F32_max((_S332), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S332))).y, s_primal_ctx_cos_0(_S332));
            if(!is_ray_depth_6)
            {
                raydir_21 = raydir_22 / make_float3 (raydir_22.z);
            }
            else
            {
                raydir_21 = raydir_22;
            }
        }
        else
        {
            float3  raydir_23 = make_float3 (uv_8.x, uv_8.y, 1.0f);
            if(is_ray_depth_6)
            {
                raydir_21 = normalize_0(raydir_23);
            }
            else
            {
                raydir_21 = raydir_23;
            }
        }
        float3  _S333 = make_float3 (dpdepths_0.x) * raydir_21;
        float2  _S334 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S324) / _S325;
        float2  _S335 = _S334;
        bool _S336 = undistort_point_0(_S334, dist_coeffs_7, int(12), &_S335);
        _s_diff_ctx_2->_S316 = _S335;
        _s_diff_ctx_2->_S317 = _S336;
        float2  uv_9 = _S335;
        bool _S337 = !_S336;
        if(_S337)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S338 = !_S337;
        if(_S338)
        {
            if(is_fisheye_6)
            {
                float _S339 = length_0(uv_9);
                float3  raydir_24 = make_float3 ((uv_9 / make_float2 ((F32_max((_S339), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S339))).x, (uv_9 / make_float2 ((F32_max((_S339), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S339))).y, s_primal_ctx_cos_0(_S339));
                if(!is_ray_depth_6)
                {
                    raydir_21 = raydir_24 / make_float3 (raydir_24.z);
                }
                else
                {
                    raydir_21 = raydir_24;
                }
            }
            else
            {
                float3  raydir_25 = make_float3 (uv_9.x, uv_9.y, 1.0f);
                if(is_ray_depth_6)
                {
                    raydir_21 = normalize_0(raydir_25);
                }
                else
                {
                    raydir_21 = raydir_25;
                }
            }
            float3  _S340 = make_float3 (dpdepths_0.y) * raydir_21;
            _S331 = int(2);
            points_3[int(0)] = _S333;
            points_3[int(1)] = _S340;
            points_3[int(2)] = _S323;
            points_3[int(3)] = _S323;
        }
        else
        {
            _S331 = int(0);
            points_3[int(0)] = _S333;
            points_3[int(1)] = _S323;
            points_3[int(2)] = _S323;
            points_3[int(3)] = _S323;
        }
        bool _runFlag_0;
        if(_S331 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S330;
            _S331 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S341 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S324) / _S325;
            float2  _S342 = _S341;
            bool _S343 = undistort_point_0(_S341, dist_coeffs_7, int(12), &_S342);
            _s_diff_ctx_2->_S318 = _S342;
            _s_diff_ctx_2->_S319 = _S343;
            float2  uv_10 = _S342;
            if(!_S343)
            {
                float3  _S344 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S331 = int(0);
                normal_4 = _S344;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_6)
                {
                    float _S345 = length_0(uv_10);
                    float3  raydir_26 = make_float3 ((uv_10 / make_float2 ((F32_max((_S345), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S345))).x, (uv_10 / make_float2 ((F32_max((_S345), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S345))).y, s_primal_ctx_cos_0(_S345));
                    if(!is_ray_depth_6)
                    {
                        raydir_21 = raydir_26 / make_float3 (raydir_26.z);
                    }
                    else
                    {
                        raydir_21 = raydir_26;
                    }
                }
                else
                {
                    float3  raydir_27 = make_float3 (uv_10.x, uv_10.y, 1.0f);
                    if(is_ray_depth_6)
                    {
                        raydir_21 = normalize_0(raydir_27);
                    }
                    else
                    {
                        raydir_21 = raydir_27;
                    }
                }
                points_3[int(2)] = make_float3 (dpdepths_0.z) * raydir_21;
                float2  _S346 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S324) / _S325;
                float2  _S347 = _S346;
                bool _S348 = undistort_point_0(_S346, dist_coeffs_7, int(12), &_S347);
                _s_diff_ctx_2->_S320 = _S347;
                _s_diff_ctx_2->_S321 = _S348;
                float2  uv_11 = _S347;
                bool _S349 = !_S348;
                if(_S349)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S350 = !_S349;
                int _S351;
                if(_S350)
                {
                    if(is_fisheye_6)
                    {
                        float _S352 = length_0(uv_11);
                        float3  raydir_28 = make_float3 ((uv_11 / make_float2 ((F32_max((_S352), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S352))).x, (uv_11 / make_float2 ((F32_max((_S352), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S352))).y, s_primal_ctx_cos_0(_S352));
                        if(!is_ray_depth_6)
                        {
                            raydir_21 = raydir_28 / make_float3 (raydir_28.z);
                        }
                        else
                        {
                            raydir_21 = raydir_28;
                        }
                    }
                    else
                    {
                        float3  raydir_29 = make_float3 (uv_11.x, uv_11.y, 1.0f);
                        if(is_ray_depth_6)
                        {
                            raydir_21 = normalize_0(raydir_29);
                        }
                        else
                        {
                            raydir_21 = raydir_29;
                        }
                    }
                    points_3[int(3)] = make_float3 (dpdepths_0.w) * raydir_21;
                    _S351 = int(2);
                }
                else
                {
                    _S351 = int(0);
                }
                if(_S351 != int(2))
                {
                    _runFlag_0 = false;
                    _S331 = _S351;
                }
                if(_runFlag_0)
                {
                    _S331 = int(1);
                }
            }
        }
    }
    else
    {
        _S331 = int(0);
        points_3[int(0)] = _S323;
        points_3[int(1)] = _S323;
        points_3[int(2)] = _S323;
        points_3[int(3)] = _S323;
    }
    if(!(_S331 != int(1)))
    {
        float3  _S353 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S353, _S353)) != 0.0f)
        {
            normal_4 = _S353 / make_float3 (length_1(_S353));
        }
        else
        {
            normal_4 = _S353;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_7, FixedArray<float, 10>  * dist_coeffs_8, bool is_fisheye_7, bool is_ray_depth_7, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_7, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S354 = *dpdepths_1;
    float3  _S355 = make_float3 (0.0f);
    float2  _S356 = _s_diff_ctx_3->_S314;
    bool _S357 = !!_s_diff_ctx_3->_S315;
    float3  raydir_30;
    float3  raydir_31;
    float3  raydir_32;
    float3  raydir_33;
    int _S358;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S359;
    if(_S357)
    {
        if(is_fisheye_7)
        {
            float _S360 = length_0(_S356);
            float3  raydir_34 = make_float3 ((_S356 / make_float2 ((F32_max((_S360), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S360))).x, (_S356 / make_float2 ((F32_max((_S360), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S360))).y, s_primal_ctx_cos_0(_S360));
            if(!is_ray_depth_7)
            {
                raydir_30 = raydir_34 / make_float3 (raydir_34.z);
            }
            else
            {
                raydir_30 = raydir_34;
            }
        }
        else
        {
            float3  raydir_35 = make_float3 (_S356.x, _S356.y, 1.0f);
            if(is_ray_depth_7)
            {
                raydir_30 = normalize_0(raydir_35);
            }
            else
            {
                raydir_30 = raydir_35;
            }
        }
        float3  _S361 = make_float3 (_S354.primal_0.x) * raydir_30;
        float2  _S362 = _s_diff_ctx_3->_S316;
        bool _S363 = !!_s_diff_ctx_3->_S317;
        if(_S363)
        {
            if(is_fisheye_7)
            {
                float _S364 = length_0(_S362);
                float3  raydir_36 = make_float3 ((_S362 / make_float2 ((F32_max((_S364), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S364))).x, (_S362 / make_float2 ((F32_max((_S364), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S364))).y, s_primal_ctx_cos_0(_S364));
                if(!is_ray_depth_7)
                {
                    raydir_31 = raydir_36 / make_float3 (raydir_36.z);
                }
                else
                {
                    raydir_31 = raydir_36;
                }
            }
            else
            {
                float3  raydir_37 = make_float3 (_S362.x, _S362.y, 1.0f);
                if(is_ray_depth_7)
                {
                    raydir_31 = normalize_0(raydir_37);
                }
                else
                {
                    raydir_31 = raydir_37;
                }
            }
            float3  _S365 = make_float3 (_S354.primal_0.y) * raydir_31;
            _S358 = int(2);
            points_4[int(0)] = _S361;
            points_4[int(1)] = _S365;
            points_4[int(2)] = _S355;
            points_4[int(3)] = _S355;
        }
        else
        {
            _S358 = int(0);
            points_4[int(0)] = _S361;
            points_4[int(1)] = _S355;
            points_4[int(2)] = _S355;
            points_4[int(3)] = _S355;
            raydir_31 = _S355;
        }
        if(_S358 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S357;
            _S358 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S366 = _s_diff_ctx_3->_S318;
            if(!_s_diff_ctx_3->_S319)
            {
                _runFlag_2 = false;
                _S358 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_7)
                {
                    float _S367 = length_0(_S366);
                    float3  raydir_38 = make_float3 ((_S366 / make_float2 ((F32_max((_S367), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S367))).x, (_S366 / make_float2 ((F32_max((_S367), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S367))).y, s_primal_ctx_cos_0(_S367));
                    if(!is_ray_depth_7)
                    {
                        raydir_32 = raydir_38 / make_float3 (raydir_38.z);
                    }
                    else
                    {
                        raydir_32 = raydir_38;
                    }
                }
                else
                {
                    float3  raydir_39 = make_float3 (_S366.x, _S366.y, 1.0f);
                    if(is_ray_depth_7)
                    {
                        raydir_32 = normalize_0(raydir_39);
                    }
                    else
                    {
                        raydir_32 = raydir_39;
                    }
                }
                points_4[int(2)] = make_float3 (_S354.primal_0.z) * raydir_32;
                float2  _S368 = _s_diff_ctx_3->_S320;
                bool _S369 = !!_s_diff_ctx_3->_S321;
                int _S370;
                if(_S369)
                {
                    if(is_fisheye_7)
                    {
                        float _S371 = length_0(_S368);
                        float3  raydir_40 = make_float3 ((_S368 / make_float2 ((F32_max((_S371), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S371))).x, (_S368 / make_float2 ((F32_max((_S371), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S371))).y, s_primal_ctx_cos_0(_S371));
                        if(!is_ray_depth_7)
                        {
                            raydir_33 = raydir_40 / make_float3 (raydir_40.z);
                        }
                        else
                        {
                            raydir_33 = raydir_40;
                        }
                    }
                    else
                    {
                        float3  raydir_41 = make_float3 (_S368.x, _S368.y, 1.0f);
                        if(is_ray_depth_7)
                        {
                            raydir_33 = normalize_0(raydir_41);
                        }
                        else
                        {
                            raydir_33 = raydir_41;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S354.primal_0.w) * raydir_33;
                    _S370 = int(2);
                }
                else
                {
                    _S370 = int(0);
                    raydir_33 = _S355;
                }
                if(_S370 != int(2))
                {
                    _runFlag_3 = false;
                    _S358 = _S370;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S358 = int(1);
                }
                float3  _S372 = raydir_32;
                _runFlag_3 = _S369;
                raydir_32 = raydir_33;
                raydir_33 = _S372;
            }
            else
            {
                _runFlag_3 = false;
                raydir_32 = _S355;
                raydir_33 = _S355;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_32 = _S355;
            raydir_33 = _S355;
        }
        float3  _S373 = raydir_30;
        float3  _S374 = raydir_31;
        raydir_30 = raydir_32;
        raydir_31 = raydir_33;
        _S359 = _S363;
        raydir_32 = _S374;
        raydir_33 = _S373;
    }
    else
    {
        _S358 = int(0);
        points_4[int(0)] = _S355;
        points_4[int(1)] = _S355;
        points_4[int(2)] = _S355;
        points_4[int(3)] = _S355;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_30 = _S355;
        raydir_31 = _S355;
        _S359 = false;
        raydir_32 = _S355;
        raydir_33 = _S355;
    }
    bool _S375 = !(_S358 != int(1));
    float3  _S376;
    float3  _S377;
    float3  _S378;
    float3  _S379;
    float3  _S380;
    bool _S381;
    if(_S375)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S382 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S383 = s_primal_ctx_cross_0(dx_1, _S382);
        bool _S384 = (s_primal_ctx_dot_0(_S383, _S383)) != 0.0f;
        if(_S384)
        {
            float _S385 = length_1(_S383);
            float3  _S386 = make_float3 (_S385);
            _S376 = make_float3 (_S385 * _S385);
            _S377 = _S386;
        }
        else
        {
            _S376 = _S355;
            _S377 = _S355;
        }
        float3  _S387 = _S377;
        _S381 = _S384;
        _S377 = _S383;
        _S378 = _S387;
        _S379 = dx_1;
        _S380 = _S382;
    }
    else
    {
        _S381 = false;
        _S376 = _S355;
        _S377 = _S355;
        _S378 = _S355;
        _S379 = _S355;
        _S380 = _S355;
    }
    float4  _S388 = make_float4 (0.0f);
    if(_S375)
    {
        if(_S381)
        {
            float3  _S389 = _s_dOut_7 / _S376;
            float3  _S390 = _S377 * - _S389;
            float3  _S391 = _S378 * _S389;
            float _S392 = _S390.x + _S390.y + _S390.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S393;
            (&_S393)->primal_0 = _S377;
            (&_S393)->differential_0 = _S355;
            s_bwd_length_impl_0(&_S393, _S392);
            _S376 = _S391 + _S393.differential_0;
        }
        else
        {
            _S376 = _s_dOut_7;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S394;
        (&_S394)->primal_0 = _S377;
        (&_S394)->differential_0 = _S355;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S395;
        (&_S395)->primal_0 = _S377;
        (&_S395)->differential_0 = _S355;
        s_bwd_prop_dot_0(&_S394, &_S395, 0.0f);
        float3  _S396 = _S395.differential_0 + _S394.differential_0 + _S376;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S397;
        (&_S397)->primal_0 = _S379;
        (&_S397)->differential_0 = _S355;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S398;
        (&_S398)->primal_0 = _S380;
        (&_S398)->differential_0 = _S355;
        s_bwd_prop_cross_0(&_S397, &_S398, _S396);
        float3  s_diff_dy_T_1 = - _S398.differential_0;
        float3  _S399 = - s_diff_dy_T_1;
        float3  _S400 = - _S397.differential_0;
        FixedArray<float3 , 4>  _S401;
        _S401[int(0)] = _S355;
        _S401[int(1)] = _S355;
        _S401[int(2)] = _S355;
        _S401[int(3)] = _S355;
        _S401[int(2)] = _S399;
        _S401[int(3)] = s_diff_dy_T_1;
        _S401[int(0)] = _S400;
        _S401[int(1)] = _S397.differential_0;
        points_4[int(0)] = _S401[int(0)];
        points_4[int(1)] = _S401[int(1)];
        points_4[int(2)] = _S401[int(2)];
        points_4[int(3)] = _S401[int(3)];
    }
    else
    {
        points_4[int(0)] = _S355;
        points_4[int(1)] = _S355;
        points_4[int(2)] = _S355;
        points_4[int(3)] = _S355;
    }
    float4  _S402;
    if(_S357)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S403 = points_4;
                FixedArray<float3 , 4>  _S404 = points_4;
                FixedArray<float3 , 4>  _S405 = points_4;
                FixedArray<float3 , 4>  _S406 = points_4;
                if(_runFlag_3)
                {
                    float3  _S407 = raydir_30 * _S406[int(3)];
                    float _S408 = _S407.x + _S407.y + _S407.z;
                    float4  _S409 = _S388;
                    *&((&_S409)->w) = _S408;
                    points_4[int(0)] = _S403[int(0)];
                    points_4[int(1)] = _S404[int(1)];
                    points_4[int(2)] = _S405[int(2)];
                    points_4[int(3)] = _S355;
                    _S402 = _S409;
                }
                else
                {
                    points_4[int(0)] = _S403[int(0)];
                    points_4[int(1)] = _S404[int(1)];
                    points_4[int(2)] = _S405[int(2)];
                    points_4[int(3)] = _S406[int(3)];
                    _S402 = _S388;
                }
                float3  _S410 = raydir_31 * points_4[int(2)];
                float _S411 = _S410.x + _S410.y + _S410.z;
                FixedArray<float3 , 4>  _S412 = points_4;
                FixedArray<float3 , 4>  _S413 = points_4;
                float4  _S414 = _S388;
                *&((&_S414)->z) = _S411;
                float4  _S415 = _S402 + _S414;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S412[int(1)];
                points_4[int(2)] = _S355;
                points_4[int(3)] = _S413[int(3)];
                _S402 = _S415;
            }
            else
            {
                FixedArray<float3 , 4>  _S416 = points_4;
                FixedArray<float3 , 4>  _S417 = points_4;
                FixedArray<float3 , 4>  _S418 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S416[int(1)];
                points_4[int(2)] = _S417[int(2)];
                points_4[int(3)] = _S418[int(3)];
                _S402 = _S388;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S419 = points_4;
            FixedArray<float3 , 4>  _S420 = points_4;
            FixedArray<float3 , 4>  _S421 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S419[int(1)];
            points_4[int(2)] = _S420[int(2)];
            points_4[int(3)] = _S421[int(3)];
            _S402 = _S388;
        }
        if(_S359)
        {
            FixedArray<float3 , 4>  _S422 = points_4;
            float3  _S423 = raydir_32 * points_4[int(1)];
            float _S424 = _S423.x + _S423.y + _S423.z;
            float4  _S425 = _S388;
            *&((&_S425)->y) = _S424;
            float4  _S426 = _S402 + _S425;
            points_4[int(0)] = _S355;
            points_4[int(1)] = _S355;
            points_4[int(2)] = _S355;
            points_4[int(3)] = _S355;
            raydir_30 = _S422[int(0)];
            _S402 = _S426;
        }
        else
        {
            FixedArray<float3 , 4>  _S427 = points_4;
            FixedArray<float3 , 4>  _S428 = points_4;
            FixedArray<float3 , 4>  _S429 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S427[int(1)];
            points_4[int(2)] = _S428[int(2)];
            points_4[int(3)] = _S429[int(3)];
            raydir_30 = _S355;
        }
        float3  _S430 = raydir_33 * (points_4[int(0)] + raydir_30);
        float _S431 = _S430.x + _S430.y + _S430.z;
        float4  _S432 = _S388;
        *&((&_S432)->x) = _S431;
        _S402 = _S402 + _S432;
    }
    else
    {
        _S402 = _S388;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S402;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S433, float4  _S434, FixedArray<float, 10>  * _S435, bool _S436, bool _S437, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S438, float3  _S439)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S440;
    float3  _S441 = s_primal_ctx_depth_to_normal_0(_S433, _S434, _S435, _S436, _S437, (*_S438).primal_0, &_S440);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S442 = _S440;
    s_bwd_prop_depth_to_normal_0(_S433, _S434, _S435, _S436, _S437, _S438, _S439, &_S442);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_8, FixedArray<float, 10>  dist_coeffs_9, bool is_fisheye_8, bool is_ray_depth_8, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S443 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S443;
    FixedArray<float, 10>  _S444 = dist_coeffs_9;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_8, &_S444, is_fisheye_8, is_ray_depth_8, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_9, FixedArray<float, 10>  dist_coeffs_10, bool is_fisheye_9)
{
    float2  _S445 = (pix_center_4 - float2 {intrins_9.z, intrins_9.w}) / float2 {intrins_9.x, intrins_9.y};
    float2  uv_12 = _S445;
    FixedArray<float, 10>  _S446 = dist_coeffs_10;
    bool _S447 = undistort_point_0(_S445, &_S446, int(12), &uv_12);
    if(!_S447)
    {
        return 0.0f;
    }
    float3  raydir_42;
    if(is_fisheye_9)
    {
        float theta_6 = length_0(uv_12);
        float3  raydir_43 = make_float3 ((uv_12 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).x, (uv_12 / make_float2 ((F32_max((theta_6), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_6))))).y, (F32_cos((theta_6))));
        raydir_42 = raydir_43 / make_float3 (raydir_43.z);
    }
    else
    {
        raydir_42 = make_float3 (uv_12.x, uv_12.y, 1.0f);
    }
    return float((F32_sign((raydir_42.z)))) / length_1(raydir_42);
}

inline __device__ float depth_normal_loss(float2  pix_center_5, float4  intrins_10, FixedArray<float, 10>  dist_coeffs_11, bool is_fisheye_10, bool is_ray_depth_9, float4  depths_2, float3  gt_normal_0)
{
    FixedArray<float3 , 5>  points_5;
    float2  _S448 = float2 {intrins_10.z, intrins_10.w};
    float2  _S449 = float2 {intrins_10.x, intrins_10.y};
    float2  _S450 = (pix_center_5 + make_float2 (-1.0f, -0.0f) - _S448) / _S449;
    float2  uv_13 = _S450;
    FixedArray<float, 10>  _S451 = dist_coeffs_11;
    bool _S452 = undistort_point_0(_S450, &_S451, int(12), &uv_13);
    if(!_S452)
    {
        return 0.0f;
    }
    float3  raydir_44;
    if(is_fisheye_10)
    {
        float theta_7 = length_0(uv_13);
        float3  raydir_45 = make_float3 ((uv_13 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).x, (uv_13 / make_float2 ((F32_max((theta_7), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_7))))).y, (F32_cos((theta_7))));
        if(!is_ray_depth_9)
        {
            raydir_44 = raydir_45 / make_float3 (raydir_45.z);
        }
        else
        {
            raydir_44 = raydir_45;
        }
    }
    else
    {
        float3  raydir_46 = make_float3 (uv_13.x, uv_13.y, 1.0f);
        if(is_ray_depth_9)
        {
            raydir_44 = normalize_0(raydir_46);
        }
        else
        {
            raydir_44 = raydir_46;
        }
    }
    points_5[int(0)] = make_float3 (depths_2.x) * raydir_44;
    float2  _S453 = (pix_center_5 + make_float2 (1.0f, -0.0f) - _S448) / _S449;
    float2  uv_14 = _S453;
    FixedArray<float, 10>  _S454 = dist_coeffs_11;
    bool _S455 = undistort_point_0(_S453, &_S454, int(12), &uv_14);
    if(!_S455)
    {
        return 0.0f;
    }
    if(is_fisheye_10)
    {
        float theta_8 = length_0(uv_14);
        float3  raydir_47 = make_float3 ((uv_14 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).x, (uv_14 / make_float2 ((F32_max((theta_8), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_8))))).y, (F32_cos((theta_8))));
        if(!is_ray_depth_9)
        {
            raydir_44 = raydir_47 / make_float3 (raydir_47.z);
        }
        else
        {
            raydir_44 = raydir_47;
        }
    }
    else
    {
        float3  raydir_48 = make_float3 (uv_14.x, uv_14.y, 1.0f);
        if(is_ray_depth_9)
        {
            raydir_44 = normalize_0(raydir_48);
        }
        else
        {
            raydir_44 = raydir_48;
        }
    }
    points_5[int(1)] = make_float3 (depths_2.y) * raydir_44;
    float2  _S456 = (pix_center_5 + make_float2 (0.0f, -1.0f) - _S448) / _S449;
    float2  uv_15 = _S456;
    FixedArray<float, 10>  _S457 = dist_coeffs_11;
    bool _S458 = undistort_point_0(_S456, &_S457, int(12), &uv_15);
    if(!_S458)
    {
        return 0.0f;
    }
    if(is_fisheye_10)
    {
        float theta_9 = length_0(uv_15);
        float3  raydir_49 = make_float3 ((uv_15 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).x, (uv_15 / make_float2 ((F32_max((theta_9), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_9))))).y, (F32_cos((theta_9))));
        if(!is_ray_depth_9)
        {
            raydir_44 = raydir_49 / make_float3 (raydir_49.z);
        }
        else
        {
            raydir_44 = raydir_49;
        }
    }
    else
    {
        float3  raydir_50 = make_float3 (uv_15.x, uv_15.y, 1.0f);
        if(is_ray_depth_9)
        {
            raydir_44 = normalize_0(raydir_50);
        }
        else
        {
            raydir_44 = raydir_50;
        }
    }
    points_5[int(2)] = make_float3 (depths_2.z) * raydir_44;
    float2  _S459 = (pix_center_5 + make_float2 (0.0f, 1.0f) - _S448) / _S449;
    float2  uv_16 = _S459;
    FixedArray<float, 10>  _S460 = dist_coeffs_11;
    bool _S461 = undistort_point_0(_S459, &_S460, int(12), &uv_16);
    if(!_S461)
    {
        return 0.0f;
    }
    if(is_fisheye_10)
    {
        float theta_10 = length_0(uv_16);
        float3  raydir_51 = make_float3 ((uv_16 / make_float2 ((F32_max((theta_10), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_10))))).x, (uv_16 / make_float2 ((F32_max((theta_10), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_10))))).y, (F32_cos((theta_10))));
        if(!is_ray_depth_9)
        {
            raydir_44 = raydir_51 / make_float3 (raydir_51.z);
        }
        else
        {
            raydir_44 = raydir_51;
        }
    }
    else
    {
        float3  raydir_52 = make_float3 (uv_16.x, uv_16.y, 1.0f);
        if(is_ray_depth_9)
        {
            raydir_44 = normalize_0(raydir_52);
        }
        else
        {
            raydir_44 = raydir_52;
        }
    }
    points_5[int(3)] = make_float3 (depths_2.w) * raydir_44;
    float2  _S462 = (pix_center_5 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S448) / _S449;
    float2  uv_17 = _S462;
    FixedArray<float, 10>  _S463 = dist_coeffs_11;
    bool _S464 = undistort_point_0(_S462, &_S463, int(12), &uv_17);
    if(!_S464)
    {
        return 0.0f;
    }
    if(is_fisheye_10)
    {
        float theta_11 = length_0(uv_17);
        float3  raydir_53 = make_float3 ((uv_17 / make_float2 ((F32_max((theta_11), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_11))))).x, (uv_17 / make_float2 ((F32_max((theta_11), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_11))))).y, (F32_cos((theta_11))));
        if(!is_ray_depth_9)
        {
            raydir_44 = raydir_53 / make_float3 (raydir_53.z);
        }
        else
        {
            raydir_44 = raydir_53;
        }
    }
    else
    {
        float3  raydir_54 = make_float3 (uv_17.x, uv_17.y, 1.0f);
        if(is_ray_depth_9)
        {
            raydir_44 = normalize_0(raydir_54);
        }
        else
        {
            raydir_44 = raydir_54;
        }
    }
    float3  normal_5 = cross_0(points_5[int(1)] - points_5[int(0)], - (points_5[int(3)] - points_5[int(2)]));
    float3  normal_6;
    if((dot_0(normal_5, normal_5)) != 0.0f)
    {
        normal_6 = normalize_0(normal_5);
    }
    else
    {
        normal_6 = normal_5;
    }
    float3  _S465;
    if((dot_0(gt_normal_0, gt_normal_0)) != 0.0f)
    {
        _S465 = normalize_0(gt_normal_0);
    }
    else
    {
        _S465 = gt_normal_0;
    }
    return (1.0f - dot_0(normal_6, _S465) + 0.00100000004749745f) / ((F32_max((dot_0(normal_6, - normalize_0(raydir_44))), (0.0f))) + 0.00100000004749745f);
}

struct s_bwd_prop_depth_normal_loss_Intermediates_0
{
    float2  _S466;
    bool _S467;
    float2  _S468;
    bool _S469;
    float2  _S470;
    bool _S471;
    float2  _S472;
    bool _S473;
    float2  _S474;
    bool _S475;
};

inline __device__ float s_primal_ctx_depth_normal_loss_0(float2  pix_center_6, float4  intrins_11, FixedArray<float, 10>  * dist_coeffs_12, bool is_fisheye_11, bool is_ray_depth_10, float4  dpdepths_2, float3  dpgt_normal_0, s_bwd_prop_depth_normal_loss_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S476 = make_float2 (0.0f);
    _s_diff_ctx_4->_S466 = _S476;
    _s_diff_ctx_4->_S467 = false;
    _s_diff_ctx_4->_S468 = _S476;
    _s_diff_ctx_4->_S469 = false;
    _s_diff_ctx_4->_S470 = _S476;
    _s_diff_ctx_4->_S471 = false;
    _s_diff_ctx_4->_S472 = _S476;
    _s_diff_ctx_4->_S473 = false;
    _s_diff_ctx_4->_S474 = _S476;
    _s_diff_ctx_4->_S475 = false;
    _s_diff_ctx_4->_S468 = _S476;
    _s_diff_ctx_4->_S469 = false;
    _s_diff_ctx_4->_S470 = _S476;
    _s_diff_ctx_4->_S471 = false;
    _s_diff_ctx_4->_S472 = _S476;
    _s_diff_ctx_4->_S473 = false;
    _s_diff_ctx_4->_S474 = _S476;
    _s_diff_ctx_4->_S475 = false;
    float3  _S477 = make_float3 (0.0f);
    float2  _S478 = float2 {intrins_11.z, intrins_11.w};
    float2  _S479 = float2 {intrins_11.x, intrins_11.y};
    float2  _S480 = (pix_center_6 + make_float2 (-1.0f, -0.0f) - _S478) / _S479;
    float2  _S481 = _S480;
    bool _S482 = undistort_point_0(_S480, dist_coeffs_12, int(12), &_S481);
    _s_diff_ctx_4->_S466 = _S481;
    _s_diff_ctx_4->_S467 = _S482;
    float2  uv_18 = _S481;
    bool _S483 = !!_S482;
    float3  raydir_55;
    int _S484;
    FixedArray<float3 , 5>  points_6;
    if(_S483)
    {
        if(is_fisheye_11)
        {
            float _S485 = length_0(uv_18);
            float3  raydir_56 = make_float3 ((uv_18 / make_float2 ((F32_max((_S485), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S485))).x, (uv_18 / make_float2 ((F32_max((_S485), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S485))).y, s_primal_ctx_cos_0(_S485));
            if(!is_ray_depth_10)
            {
                raydir_55 = raydir_56 / make_float3 (raydir_56.z);
            }
            else
            {
                raydir_55 = raydir_56;
            }
        }
        else
        {
            float3  raydir_57 = make_float3 (uv_18.x, uv_18.y, 1.0f);
            if(is_ray_depth_10)
            {
                raydir_55 = normalize_0(raydir_57);
            }
            else
            {
                raydir_55 = raydir_57;
            }
        }
        float3  _S486 = make_float3 (dpdepths_2.x) * raydir_55;
        float2  _S487 = (pix_center_6 + make_float2 (1.0f, -0.0f) - _S478) / _S479;
        float2  _S488 = _S487;
        bool _S489 = undistort_point_0(_S487, dist_coeffs_12, int(12), &_S488);
        _s_diff_ctx_4->_S468 = _S488;
        _s_diff_ctx_4->_S469 = _S489;
        float2  uv_19 = _S488;
        bool _runFlag_4;
        if(!_S489)
        {
            _runFlag_4 = false;
        }
        else
        {
            _runFlag_4 = _S483;
        }
        if(_runFlag_4)
        {
            if(is_fisheye_11)
            {
                float _S490 = length_0(uv_19);
                float3  raydir_58 = make_float3 ((uv_19 / make_float2 ((F32_max((_S490), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S490))).x, (uv_19 / make_float2 ((F32_max((_S490), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S490))).y, s_primal_ctx_cos_0(_S490));
                if(!is_ray_depth_10)
                {
                    raydir_55 = raydir_58 / make_float3 (raydir_58.z);
                }
                else
                {
                    raydir_55 = raydir_58;
                }
            }
            else
            {
                float3  raydir_59 = make_float3 (uv_19.x, uv_19.y, 1.0f);
                if(is_ray_depth_10)
                {
                    raydir_55 = normalize_0(raydir_59);
                }
                else
                {
                    raydir_55 = raydir_59;
                }
            }
            float3  _S491 = make_float3 (dpdepths_2.y) * raydir_55;
            float2  _S492 = (pix_center_6 + make_float2 (0.0f, -1.0f) - _S478) / _S479;
            float2  _S493 = _S492;
            bool _S494 = undistort_point_0(_S492, dist_coeffs_12, int(12), &_S493);
            _s_diff_ctx_4->_S470 = _S493;
            _s_diff_ctx_4->_S471 = _S494;
            float2  uv_20 = _S493;
            if(!_S494)
            {
                _runFlag_4 = false;
            }
            if(_runFlag_4)
            {
                if(is_fisheye_11)
                {
                    float _S495 = length_0(uv_20);
                    float3  raydir_60 = make_float3 ((uv_20 / make_float2 ((F32_max((_S495), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S495))).x, (uv_20 / make_float2 ((F32_max((_S495), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S495))).y, s_primal_ctx_cos_0(_S495));
                    if(!is_ray_depth_10)
                    {
                        raydir_55 = raydir_60 / make_float3 (raydir_60.z);
                    }
                    else
                    {
                        raydir_55 = raydir_60;
                    }
                }
                else
                {
                    float3  raydir_61 = make_float3 (uv_20.x, uv_20.y, 1.0f);
                    if(is_ray_depth_10)
                    {
                        raydir_55 = normalize_0(raydir_61);
                    }
                    else
                    {
                        raydir_55 = raydir_61;
                    }
                }
                float3  _S496 = make_float3 (dpdepths_2.z) * raydir_55;
                float2  _S497 = (pix_center_6 + make_float2 (0.0f, 1.0f) - _S478) / _S479;
                float2  _S498 = _S497;
                bool _S499 = undistort_point_0(_S497, dist_coeffs_12, int(12), &_S498);
                _s_diff_ctx_4->_S472 = _S498;
                _s_diff_ctx_4->_S473 = _S499;
                float2  uv_21 = _S498;
                if(!_S499)
                {
                    _runFlag_4 = false;
                }
                if(_runFlag_4)
                {
                    if(is_fisheye_11)
                    {
                        float _S500 = length_0(uv_21);
                        float3  raydir_62 = make_float3 ((uv_21 / make_float2 ((F32_max((_S500), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S500))).x, (uv_21 / make_float2 ((F32_max((_S500), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S500))).y, s_primal_ctx_cos_0(_S500));
                        if(!is_ray_depth_10)
                        {
                            raydir_55 = raydir_62 / make_float3 (raydir_62.z);
                        }
                        else
                        {
                            raydir_55 = raydir_62;
                        }
                    }
                    else
                    {
                        float3  raydir_63 = make_float3 (uv_21.x, uv_21.y, 1.0f);
                        if(is_ray_depth_10)
                        {
                            raydir_55 = normalize_0(raydir_63);
                        }
                        else
                        {
                            raydir_55 = raydir_63;
                        }
                    }
                    float3  _S501 = make_float3 (dpdepths_2.w) * raydir_55;
                    float2  _S502 = (pix_center_6 - _S478) / _S479;
                    float2  _S503 = _S502;
                    bool _S504 = undistort_point_0(_S502, dist_coeffs_12, int(12), &_S503);
                    _s_diff_ctx_4->_S474 = _S503;
                    _s_diff_ctx_4->_S475 = _S504;
                    float2  uv_22 = _S503;
                    if(!_S504)
                    {
                        _runFlag_4 = false;
                    }
                    if(_runFlag_4)
                    {
                        if(is_fisheye_11)
                        {
                            float _S505 = length_0(uv_22);
                            float3  raydir_64 = make_float3 ((uv_22 / make_float2 ((F32_max((_S505), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S505))).x, (uv_22 / make_float2 ((F32_max((_S505), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S505))).y, s_primal_ctx_cos_0(_S505));
                            if(!is_ray_depth_10)
                            {
                                raydir_55 = raydir_64 / make_float3 (raydir_64.z);
                            }
                            else
                            {
                                raydir_55 = raydir_64;
                            }
                        }
                        else
                        {
                            float3  raydir_65 = make_float3 (uv_22.x, uv_22.y, 1.0f);
                            if(is_ray_depth_10)
                            {
                                raydir_55 = normalize_0(raydir_65);
                            }
                            else
                            {
                                raydir_55 = raydir_65;
                            }
                        }
                        _S484 = int(1);
                    }
                    else
                    {
                        _S484 = int(0);
                    }
                    points_6[int(0)] = _S486;
                    points_6[int(1)] = _S491;
                    points_6[int(2)] = _S496;
                    points_6[int(3)] = _S501;
                    points_6[int(4)] = _S477;
                }
                else
                {
                    _S484 = int(0);
                    points_6[int(0)] = _S486;
                    points_6[int(1)] = _S491;
                    points_6[int(2)] = _S496;
                    points_6[int(3)] = _S477;
                    points_6[int(4)] = _S477;
                }
            }
            else
            {
                _S484 = int(0);
                points_6[int(0)] = _S486;
                points_6[int(1)] = _S491;
                points_6[int(2)] = _S477;
                points_6[int(3)] = _S477;
                points_6[int(4)] = _S477;
            }
        }
        else
        {
            _S484 = int(0);
            points_6[int(0)] = _S486;
            points_6[int(1)] = _S477;
            points_6[int(2)] = _S477;
            points_6[int(3)] = _S477;
            points_6[int(4)] = _S477;
        }
    }
    else
    {
        _S484 = int(0);
        points_6[int(0)] = _S477;
        points_6[int(1)] = _S477;
        points_6[int(2)] = _S477;
        points_6[int(3)] = _S477;
        points_6[int(4)] = _S477;
    }
    float _S506;
    if(!(_S484 != int(1)))
    {
        float3  _S507 = s_primal_ctx_cross_0(points_6[int(1)] - points_6[int(0)], - (points_6[int(3)] - points_6[int(2)]));
        float3  normal_7;
        if((s_primal_ctx_dot_0(_S507, _S507)) != 0.0f)
        {
            normal_7 = normalize_0(_S507);
        }
        else
        {
            normal_7 = _S507;
        }
        float3  _S508;
        if((s_primal_ctx_dot_0(dpgt_normal_0, dpgt_normal_0)) != 0.0f)
        {
            _S508 = normalize_0(dpgt_normal_0);
        }
        else
        {
            _S508 = dpgt_normal_0;
        }
        _S506 = (1.0f - s_primal_ctx_dot_0(normal_7, _S508) + 0.00100000004749745f) / ((F32_max((s_primal_ctx_dot_0(normal_7, - normalize_0(raydir_55))), (0.0f))) + 0.00100000004749745f);
    }
    else
    {
        _S506 = 0.0f;
    }
    return _S506;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  _s_dOut_8)
{
    float _S509 = length_1((*dpx_8).primal_0);
    float3  _S510 = (*dpx_8).primal_0 * _s_dOut_8;
    float3  _S511 = make_float3 (1.0f / _S509) * _s_dOut_8;
    float _S512 = - ((_S510.x + _S510.y + _S510.z) / (_S509 * _S509));
    float3  _S513 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S514;
    (&_S514)->primal_0 = (*dpx_8).primal_0;
    (&_S514)->differential_0 = _S513;
    s_bwd_length_impl_0(&_S514, _S512);
    float3  _S515 = _S511 + _S514.differential_0;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S515;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S516, float3  _S517)
{
    s_bwd_prop_normalize_impl_0(_S516, _S517);
    return;
}

inline __device__ void s_bwd_prop_depth_normal_loss_0(float2  pix_center_7, float4  intrins_12, FixedArray<float, 10>  * dist_coeffs_13, bool is_fisheye_12, bool is_ray_depth_11, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpgt_normal_1, float _s_dOut_9, s_bwd_prop_depth_normal_loss_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S518 = *dpdepths_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S519 = *dpgt_normal_1;
    float3  _S520 = make_float3 (0.0f);
    float2  _S521 = _s_diff_ctx_5->_S466;
    bool _S522 = !!_s_diff_ctx_5->_S467;
    float3  raydir_66;
    float3  raydir_67;
    float3  raydir_68;
    float3  raydir_69;
    float3  raydir_70;
    bool _runFlag_5;
    bool _runFlag_6;
    bool _runFlag_7;
    bool _runFlag_8;
    int _S523;
    FixedArray<float3 , 5>  points_7;
    if(_S522)
    {
        if(is_fisheye_12)
        {
            float _S524 = length_0(_S521);
            float3  raydir_71 = make_float3 ((_S521 / make_float2 ((F32_max((_S524), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S524))).x, (_S521 / make_float2 ((F32_max((_S524), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S524))).y, s_primal_ctx_cos_0(_S524));
            if(!is_ray_depth_11)
            {
                raydir_66 = raydir_71 / make_float3 (raydir_71.z);
            }
            else
            {
                raydir_66 = raydir_71;
            }
        }
        else
        {
            float3  raydir_72 = make_float3 (_S521.x, _S521.y, 1.0f);
            if(is_ray_depth_11)
            {
                raydir_66 = normalize_0(raydir_72);
            }
            else
            {
                raydir_66 = raydir_72;
            }
        }
        float3  _S525 = make_float3 (_S518.primal_0.x) * raydir_66;
        float2  _S526 = _s_diff_ctx_5->_S468;
        if(!_s_diff_ctx_5->_S469)
        {
            _runFlag_5 = false;
        }
        else
        {
            _runFlag_5 = _S522;
        }
        if(_runFlag_5)
        {
            if(is_fisheye_12)
            {
                float _S527 = length_0(_S526);
                float3  raydir_73 = make_float3 ((_S526 / make_float2 ((F32_max((_S527), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S527))).x, (_S526 / make_float2 ((F32_max((_S527), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S527))).y, s_primal_ctx_cos_0(_S527));
                if(!is_ray_depth_11)
                {
                    raydir_67 = raydir_73 / make_float3 (raydir_73.z);
                }
                else
                {
                    raydir_67 = raydir_73;
                }
            }
            else
            {
                float3  raydir_74 = make_float3 (_S526.x, _S526.y, 1.0f);
                if(is_ray_depth_11)
                {
                    raydir_67 = normalize_0(raydir_74);
                }
                else
                {
                    raydir_67 = raydir_74;
                }
            }
            float3  _S528 = make_float3 (_S518.primal_0.y) * raydir_67;
            float2  _S529 = _s_diff_ctx_5->_S470;
            if(!_s_diff_ctx_5->_S471)
            {
                _runFlag_6 = false;
            }
            else
            {
                _runFlag_6 = _runFlag_5;
            }
            if(_runFlag_6)
            {
                if(is_fisheye_12)
                {
                    float _S530 = length_0(_S529);
                    float3  raydir_75 = make_float3 ((_S529 / make_float2 ((F32_max((_S530), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S530))).x, (_S529 / make_float2 ((F32_max((_S530), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S530))).y, s_primal_ctx_cos_0(_S530));
                    if(!is_ray_depth_11)
                    {
                        raydir_68 = raydir_75 / make_float3 (raydir_75.z);
                    }
                    else
                    {
                        raydir_68 = raydir_75;
                    }
                }
                else
                {
                    float3  raydir_76 = make_float3 (_S529.x, _S529.y, 1.0f);
                    if(is_ray_depth_11)
                    {
                        raydir_68 = normalize_0(raydir_76);
                    }
                    else
                    {
                        raydir_68 = raydir_76;
                    }
                }
                float3  _S531 = make_float3 (_S518.primal_0.z) * raydir_68;
                float2  _S532 = _s_diff_ctx_5->_S472;
                if(!_s_diff_ctx_5->_S473)
                {
                    _runFlag_7 = false;
                }
                else
                {
                    _runFlag_7 = _runFlag_6;
                }
                if(_runFlag_7)
                {
                    if(is_fisheye_12)
                    {
                        float _S533 = length_0(_S532);
                        float3  raydir_77 = make_float3 ((_S532 / make_float2 ((F32_max((_S533), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S533))).x, (_S532 / make_float2 ((F32_max((_S533), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S533))).y, s_primal_ctx_cos_0(_S533));
                        if(!is_ray_depth_11)
                        {
                            raydir_69 = raydir_77 / make_float3 (raydir_77.z);
                        }
                        else
                        {
                            raydir_69 = raydir_77;
                        }
                    }
                    else
                    {
                        float3  raydir_78 = make_float3 (_S532.x, _S532.y, 1.0f);
                        if(is_ray_depth_11)
                        {
                            raydir_69 = normalize_0(raydir_78);
                        }
                        else
                        {
                            raydir_69 = raydir_78;
                        }
                    }
                    float3  _S534 = make_float3 (_S518.primal_0.w) * raydir_69;
                    float2  _S535 = _s_diff_ctx_5->_S474;
                    if(!_s_diff_ctx_5->_S475)
                    {
                        _runFlag_8 = false;
                    }
                    else
                    {
                        _runFlag_8 = _runFlag_7;
                    }
                    if(_runFlag_8)
                    {
                        if(is_fisheye_12)
                        {
                            float _S536 = length_0(_S535);
                            float3  raydir_79 = make_float3 ((_S535 / make_float2 ((F32_max((_S536), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S536))).x, (_S535 / make_float2 ((F32_max((_S536), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S536))).y, s_primal_ctx_cos_0(_S536));
                            if(!is_ray_depth_11)
                            {
                                raydir_70 = raydir_79 / make_float3 (raydir_79.z);
                            }
                            else
                            {
                                raydir_70 = raydir_79;
                            }
                        }
                        else
                        {
                            float3  raydir_80 = make_float3 (_S535.x, _S535.y, 1.0f);
                            if(is_ray_depth_11)
                            {
                                raydir_70 = normalize_0(raydir_80);
                            }
                            else
                            {
                                raydir_70 = raydir_80;
                            }
                        }
                        _S523 = int(1);
                    }
                    else
                    {
                        _S523 = int(0);
                        raydir_70 = raydir_69;
                    }
                    points_7[int(0)] = _S525;
                    points_7[int(1)] = _S528;
                    points_7[int(2)] = _S531;
                    points_7[int(3)] = _S534;
                    points_7[int(4)] = _S520;
                }
                else
                {
                    _S523 = int(0);
                    raydir_70 = raydir_68;
                    points_7[int(0)] = _S525;
                    points_7[int(1)] = _S528;
                    points_7[int(2)] = _S531;
                    points_7[int(3)] = _S520;
                    points_7[int(4)] = _S520;
                    raydir_69 = _S520;
                }
                float3  _S537 = raydir_68;
                raydir_68 = raydir_69;
                raydir_69 = _S537;
            }
            else
            {
                _S523 = int(0);
                raydir_70 = raydir_67;
                points_7[int(0)] = _S525;
                points_7[int(1)] = _S528;
                points_7[int(2)] = _S520;
                points_7[int(3)] = _S520;
                points_7[int(4)] = _S520;
                _runFlag_7 = false;
                raydir_68 = _S520;
                raydir_69 = _S520;
            }
            float3  _S538 = raydir_67;
            raydir_67 = raydir_68;
            raydir_68 = raydir_69;
            raydir_69 = _S538;
        }
        else
        {
            _S523 = int(0);
            raydir_70 = raydir_66;
            points_7[int(0)] = _S525;
            points_7[int(1)] = _S520;
            points_7[int(2)] = _S520;
            points_7[int(3)] = _S520;
            points_7[int(4)] = _S520;
            _runFlag_6 = false;
            _runFlag_7 = false;
            raydir_67 = _S520;
            raydir_68 = _S520;
            raydir_69 = _S520;
        }
        float3  _S539 = raydir_66;
        raydir_66 = raydir_67;
        raydir_67 = raydir_68;
        raydir_68 = raydir_69;
        raydir_69 = _S539;
    }
    else
    {
        _S523 = int(0);
        points_7[int(0)] = _S520;
        points_7[int(1)] = _S520;
        points_7[int(2)] = _S520;
        points_7[int(3)] = _S520;
        points_7[int(4)] = _S520;
        _runFlag_5 = false;
        _runFlag_6 = false;
        _runFlag_7 = false;
        raydir_66 = _S520;
        raydir_67 = _S520;
        raydir_68 = _S520;
        raydir_69 = _S520;
    }
    bool _S540 = !(_S523 != int(1));
    float3  normal_8;
    float3  _S541;
    float3  _S542;
    float3  _S543;
    float3  _S544;
    bool _S545;
    float _S546;
    float _S547;
    float _S548;
    float _S549;
    if(_S540)
    {
        float3  dx_2 = points_7[int(1)] - points_7[int(0)];
        float3  _S550 = - (points_7[int(3)] - points_7[int(2)]);
        float3  _S551 = s_primal_ctx_cross_0(dx_2, _S550);
        bool _S552 = (s_primal_ctx_dot_0(_S551, _S551)) != 0.0f;
        if(_S552)
        {
            normal_8 = normalize_0(_S551);
        }
        else
        {
            normal_8 = _S551;
        }
        bool _S553 = (s_primal_ctx_dot_0(_S519.primal_0, _S519.primal_0)) != 0.0f;
        if(_S553)
        {
            _S541 = normalize_0(_S519.primal_0);
        }
        else
        {
            _S541 = _S519.primal_0;
        }
        float3  _S554 = - normalize_0(raydir_70);
        float _S555 = s_primal_ctx_dot_0(normal_8, _S554);
        float _S556 = 1.0f - s_primal_ctx_dot_0(normal_8, _S541) + 0.00100000004749745f;
        float _S557 = (F32_max((_S555), (0.0f))) + 0.00100000004749745f;
        _S546 = _S557 * _S557;
        _S547 = _S556;
        _S548 = _S557;
        _S549 = _S555;
        raydir_70 = normal_8;
        normal_8 = _S554;
        _runFlag_8 = _S553;
        _S545 = _S552;
        _S542 = _S551;
        _S543 = dx_2;
        _S544 = _S550;
    }
    else
    {
        _S546 = 0.0f;
        _S547 = 0.0f;
        _S548 = 0.0f;
        _S549 = 0.0f;
        raydir_70 = _S520;
        normal_8 = _S520;
        _S541 = _S520;
        _runFlag_8 = false;
        _S545 = false;
        _S542 = _S520;
        _S543 = _S520;
        _S544 = _S520;
    }
    float4  _S558 = make_float4 (0.0f);
    if(_S540)
    {
        float _S559 = _s_dOut_9 / _S546;
        float _S560 = _S547 * - _S559;
        float s_diff_num_T_0 = _S548 * _S559;
        DiffPair_float_0 _S561;
        (&_S561)->primal_0 = _S549;
        (&_S561)->differential_0 = 0.0f;
        DiffPair_float_0 _S562;
        (&_S562)->primal_0 = 0.0f;
        (&_S562)->differential_0 = 0.0f;
        _d_max_0(&_S561, &_S562, _S560);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S563;
        (&_S563)->primal_0 = raydir_70;
        (&_S563)->differential_0 = _S520;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S564;
        (&_S564)->primal_0 = normal_8;
        (&_S564)->differential_0 = _S520;
        s_bwd_prop_dot_0(&_S563, &_S564, _S561.differential_0);
        float _S565 = - s_diff_num_T_0;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S566;
        (&_S566)->primal_0 = raydir_70;
        (&_S566)->differential_0 = _S520;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S567;
        (&_S567)->primal_0 = _S541;
        (&_S567)->differential_0 = _S520;
        s_bwd_prop_dot_0(&_S566, &_S567, _S565);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S568 = _S567;
        float3  _S569 = _S563.differential_0 + _S566.differential_0;
        if(_runFlag_8)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S570;
            (&_S570)->primal_0 = _S519.primal_0;
            (&_S570)->differential_0 = _S520;
            s_bwd_normalize_impl_0(&_S570, _S568.differential_0);
            raydir_70 = _S570.differential_0;
        }
        else
        {
            raydir_70 = _S568.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S571;
        (&_S571)->primal_0 = _S519.primal_0;
        (&_S571)->differential_0 = _S520;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S572;
        (&_S572)->primal_0 = _S519.primal_0;
        (&_S572)->differential_0 = _S520;
        s_bwd_prop_dot_0(&_S571, &_S572, 0.0f);
        float3  _S573 = _S572.differential_0 + _S571.differential_0 + raydir_70;
        if(_S545)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S574;
            (&_S574)->primal_0 = _S542;
            (&_S574)->differential_0 = _S520;
            s_bwd_normalize_impl_0(&_S574, _S569);
            raydir_70 = _S574.differential_0;
        }
        else
        {
            raydir_70 = _S569;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S575;
        (&_S575)->primal_0 = _S542;
        (&_S575)->differential_0 = _S520;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S576;
        (&_S576)->primal_0 = _S542;
        (&_S576)->differential_0 = _S520;
        s_bwd_prop_dot_0(&_S575, &_S576, 0.0f);
        float3  _S577 = _S576.differential_0 + _S575.differential_0 + raydir_70;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S578;
        (&_S578)->primal_0 = _S543;
        (&_S578)->differential_0 = _S520;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S579;
        (&_S579)->primal_0 = _S544;
        (&_S579)->differential_0 = _S520;
        s_bwd_prop_cross_0(&_S578, &_S579, _S577);
        float3  s_diff_dy_T_2 = - _S579.differential_0;
        float3  _S580 = - s_diff_dy_T_2;
        float3  _S581 = - _S578.differential_0;
        FixedArray<float3 , 5>  _S582;
        _S582[int(0)] = _S520;
        _S582[int(1)] = _S520;
        _S582[int(2)] = _S520;
        _S582[int(3)] = _S520;
        _S582[int(4)] = _S520;
        _S582[int(2)] = _S580;
        _S582[int(3)] = s_diff_dy_T_2;
        _S582[int(0)] = _S581;
        _S582[int(1)] = _S578.differential_0;
        points_7[int(0)] = _S582[int(0)];
        points_7[int(1)] = _S582[int(1)];
        points_7[int(2)] = _S582[int(2)];
        points_7[int(3)] = _S582[int(3)];
        points_7[int(4)] = _S582[int(4)];
        raydir_70 = _S573;
    }
    else
    {
        points_7[int(0)] = _S520;
        points_7[int(1)] = _S520;
        points_7[int(2)] = _S520;
        points_7[int(3)] = _S520;
        points_7[int(4)] = _S520;
        raydir_70 = _S520;
    }
    float4  _S583;
    if(_S522)
    {
        if(_runFlag_5)
        {
            if(_runFlag_6)
            {
                if(_runFlag_7)
                {
                    FixedArray<float3 , 5>  _S584 = points_7;
                    FixedArray<float3 , 5>  _S585 = points_7;
                    FixedArray<float3 , 5>  _S586 = points_7;
                    float3  _S587 = raydir_66 * points_7[int(3)];
                    float _S588 = _S587.x + _S587.y + _S587.z;
                    float4  _S589 = _S558;
                    *&((&_S589)->w) = _S588;
                    points_7[int(0)] = _S520;
                    points_7[int(1)] = _S520;
                    points_7[int(2)] = _S520;
                    points_7[int(3)] = _S520;
                    points_7[int(4)] = _S520;
                    raydir_66 = _S586[int(2)];
                    normal_8 = _S584[int(0)];
                    _S541 = _S585[int(1)];
                    _S583 = _S589;
                }
                else
                {
                    FixedArray<float3 , 5>  _S590 = points_7;
                    FixedArray<float3 , 5>  _S591 = points_7;
                    FixedArray<float3 , 5>  _S592 = points_7;
                    FixedArray<float3 , 5>  _S593 = points_7;
                    points_7[int(0)] = points_7[int(0)];
                    points_7[int(1)] = _S590[int(1)];
                    points_7[int(2)] = _S591[int(2)];
                    points_7[int(3)] = _S592[int(3)];
                    points_7[int(4)] = _S593[int(4)];
                    raydir_66 = _S520;
                    normal_8 = _S520;
                    _S541 = _S520;
                    _S583 = _S558;
                }
                float3  _S594 = raydir_67 * (points_7[int(2)] + raydir_66);
                float _S595 = _S594.x + _S594.y + _S594.z;
                float3  _S596 = points_7[int(0)] + normal_8;
                float3  _S597 = points_7[int(1)] + _S541;
                float4  _S598 = _S558;
                *&((&_S598)->z) = _S595;
                float4  _S599 = _S583 + _S598;
                points_7[int(0)] = _S520;
                points_7[int(1)] = _S520;
                points_7[int(2)] = _S520;
                points_7[int(3)] = _S520;
                points_7[int(4)] = _S520;
                raydir_66 = _S597;
                raydir_67 = _S596;
                _S583 = _S599;
            }
            else
            {
                FixedArray<float3 , 5>  _S600 = points_7;
                FixedArray<float3 , 5>  _S601 = points_7;
                FixedArray<float3 , 5>  _S602 = points_7;
                FixedArray<float3 , 5>  _S603 = points_7;
                points_7[int(0)] = points_7[int(0)];
                points_7[int(1)] = _S600[int(1)];
                points_7[int(2)] = _S601[int(2)];
                points_7[int(3)] = _S602[int(3)];
                points_7[int(4)] = _S603[int(4)];
                raydir_66 = _S520;
                raydir_67 = _S520;
                _S583 = _S558;
            }
            float3  _S604 = raydir_68 * (points_7[int(1)] + raydir_66);
            float _S605 = _S604.x + _S604.y + _S604.z;
            float3  _S606 = points_7[int(0)] + raydir_67;
            float4  _S607 = _S558;
            *&((&_S607)->y) = _S605;
            float4  _S608 = _S583 + _S607;
            points_7[int(0)] = _S520;
            points_7[int(1)] = _S520;
            points_7[int(2)] = _S520;
            points_7[int(3)] = _S520;
            points_7[int(4)] = _S520;
            raydir_66 = _S606;
            _S583 = _S608;
        }
        else
        {
            FixedArray<float3 , 5>  _S609 = points_7;
            FixedArray<float3 , 5>  _S610 = points_7;
            FixedArray<float3 , 5>  _S611 = points_7;
            FixedArray<float3 , 5>  _S612 = points_7;
            points_7[int(0)] = points_7[int(0)];
            points_7[int(1)] = _S609[int(1)];
            points_7[int(2)] = _S610[int(2)];
            points_7[int(3)] = _S611[int(3)];
            points_7[int(4)] = _S612[int(4)];
            raydir_66 = _S520;
            _S583 = _S558;
        }
        float3  _S613 = raydir_69 * (points_7[int(0)] + raydir_66);
        float _S614 = _S613.x + _S613.y + _S613.z;
        float4  _S615 = _S558;
        *&((&_S615)->x) = _S614;
        _S583 = _S583 + _S615;
    }
    else
    {
        _S583 = _S558;
    }
    dpgt_normal_1->primal_0 = (*dpgt_normal_1).primal_0;
    dpgt_normal_1->differential_0 = raydir_70;
    dpdepths_3->primal_0 = (*dpdepths_3).primal_0;
    dpdepths_3->differential_0 = _S583;
    return;
}

inline __device__ void s_bwd_depth_normal_loss_0(float2  _S616, float4  _S617, FixedArray<float, 10>  * _S618, bool _S619, bool _S620, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S621, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S622, float _S623)
{
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S624;
    float _S625 = s_primal_ctx_depth_normal_loss_0(_S616, _S617, _S618, _S619, _S620, (*_S621).primal_0, (*_S622).primal_0, &_S624);
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S626 = _S624;
    s_bwd_prop_depth_normal_loss_0(_S616, _S617, _S618, _S619, _S620, _S621, _S622, _S623, &_S626);
    return;
}

inline __device__ void depth_normal_loss_vjp(float2  pix_center_8, float4  intrins_13, FixedArray<float, 10>  dist_coeffs_14, bool is_fisheye_13, bool is_ray_depth_12, float4  depths_3, float3  gt_normal_1, float v_loss_0, float4  * v_depths_1, float3  * v_gt_normal_0)
{
    float4  _S627 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_1;
    (&dp_depths_1)->primal_0 = depths_3;
    (&dp_depths_1)->differential_0 = _S627;
    float3  _S628 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_gt_normal_0;
    (&dp_gt_normal_0)->primal_0 = gt_normal_1;
    (&dp_gt_normal_0)->differential_0 = _S628;
    FixedArray<float, 10>  _S629 = dist_coeffs_14;
    s_bwd_depth_normal_loss_0(pix_center_8, intrins_13, &_S629, is_fisheye_13, is_ray_depth_12, &dp_depths_1, &dp_gt_normal_0, v_loss_0);
    *v_depths_1 = dp_depths_1.differential_0;
    *v_gt_normal_0 = dp_gt_normal_0.differential_0;
    return;
}

