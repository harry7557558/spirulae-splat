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

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_7)
{
    float _S208 = dOut_7.y;
    float _S209 = dOut_7.z;
    float _S210 = dOut_7.x;
    float _S211 = (*a_0).primal_0.z * _S208 + - (*a_0).primal_0.y * _S209;
    float _S212 = - (*a_0).primal_0.z * _S210 + (*a_0).primal_0.x * _S209;
    float _S213 = (*a_0).primal_0.y * _S210 + - (*a_0).primal_0.x * _S208;
    float3  _S214 = make_float3 (- (*b_0).primal_0.z * _S208 + (*b_0).primal_0.y * _S209, (*b_0).primal_0.z * _S210 + - (*b_0).primal_0.x * _S209, - (*b_0).primal_0.y * _S210 + (*b_0).primal_0.x * _S208);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S214;
    float3  _S215 = make_float3 (_S211, _S212, _S213);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S215;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S216 = left_2.y;
    float _S217 = right_2.z;
    float _S218 = left_2.z;
    float _S219 = right_2.y;
    float _S220 = right_2.x;
    float _S221 = left_2.x;
    return make_float3 (_S216 * _S217 - _S218 * _S219, _S218 * _S220 - _S221 * _S217, _S221 * _S219 - _S216 * _S220);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S222, float3  _S223)
{
    return cross_0(_S222, _S223);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S224, float3  _S225)
{
    return dot_0(_S224, _S225);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S226, float _S227)
{
    _d_sqrt_0(_S226, _S227);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_4)
{
    float _S228 = (*dpx_7).primal_0.x;
    float _S229 = (*dpx_7).primal_0.y;
    float _S230 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S231;
    (&_S231)->primal_0 = _S228 * _S228 + _S229 * _S229 + _S230 * _S230;
    (&_S231)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S231, _s_dOut_4);
    float _S232 = (*dpx_7).primal_0.z * _S231.differential_0;
    float _S233 = _S232 + _S232;
    float _S234 = (*dpx_7).primal_0.y * _S231.differential_0;
    float _S235 = _S234 + _S234;
    float _S236 = (*dpx_7).primal_0.x * _S231.differential_0;
    float _S237 = _S236 + _S236;
    float3  _S238 = make_float3 (0.0f);
    *&((&_S238)->z) = _S233;
    *&((&_S238)->y) = _S235;
    *&((&_S238)->x) = _S237;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S238;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S239, float _S240)
{
    s_bwd_prop_length_impl_0(_S239, _S240);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S241, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S242, float _S243)
{
    _d_dot_0(_S241, _S242, _S243);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S244, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S245, float3  _S246)
{
    _d_cross_0(_S244, _S245, _S246);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_5)
{
    float3  _S247 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S248 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S249 = s_primal_ctx_cross_0(dx_0, _S248);
    bool _S250 = (s_primal_ctx_dot_0(_S249, _S249)) != 0.0f;
    float3  _S251;
    float3  _S252;
    if(_S250)
    {
        float _S253 = length_1(_S249);
        float3  _S254 = make_float3 (_S253);
        _S251 = make_float3 (_S253 * _S253);
        _S252 = _S254;
    }
    else
    {
        _S251 = _S247;
        _S252 = _S247;
    }
    if(_S250)
    {
        float3  _S255 = _s_dOut_5 / _S251;
        float3  _S256 = _S249 * - _S255;
        float3  _S257 = _S252 * _S255;
        float _S258 = _S256.x + _S256.y + _S256.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S259;
        (&_S259)->primal_0 = _S249;
        (&_S259)->differential_0 = _S247;
        s_bwd_length_impl_0(&_S259, _S258);
        _S251 = _S257 + _S259.differential_0;
    }
    else
    {
        _S251 = _s_dOut_5;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S260;
    (&_S260)->primal_0 = _S249;
    (&_S260)->differential_0 = _S247;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S261;
    (&_S261)->primal_0 = _S249;
    (&_S261)->differential_0 = _S247;
    s_bwd_prop_dot_0(&_S260, &_S261, 0.0f);
    float3  _S262 = _S261.differential_0 + _S260.differential_0 + _S251;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S263;
    (&_S263)->primal_0 = dx_0;
    (&_S263)->differential_0 = _S247;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S264;
    (&_S264)->primal_0 = _S248;
    (&_S264)->differential_0 = _S247;
    s_bwd_prop_cross_0(&_S263, &_S264, _S262);
    float3  s_diff_dy_T_0 = - _S264.differential_0;
    float3  _S265 = - s_diff_dy_T_0;
    float3  _S266 = - _S263.differential_0;
    FixedArray<float3 , 4>  _S267;
    _S267[int(0)] = _S247;
    _S267[int(1)] = _S247;
    _S267[int(2)] = _S247;
    _S267[int(3)] = _S247;
    _S267[int(2)] = _S265;
    _S267[int(3)] = s_diff_dy_T_0;
    _S267[int(0)] = _S266;
    _S267[int(1)] = _S263.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S267;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S268, float3  _S269)
{
    s_bwd_prop_points_to_normal_0(_S268, _S269);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S270 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S270;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S271 = float2 {intrins_1.z, intrins_1.w};
    float2  _S272 = float2 {intrins_1.x, intrins_1.y};
    float2  _S273 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S271) / _S272;
    float2  uv_2 = _S273;
    FixedArray<float, 10>  _S274 = dist_coeffs_2;
    bool _S275 = undistort_point_0(_S273, &_S274, int(12), &uv_2);
    if(!_S275)
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
    float2  _S276 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S271) / _S272;
    float2  uv_3 = _S276;
    FixedArray<float, 10>  _S277 = dist_coeffs_2;
    bool _S278 = undistort_point_0(_S276, &_S277, int(12), &uv_3);
    if(!_S278)
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
    float2  _S279 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S271) / _S272;
    float2  uv_4 = _S279;
    FixedArray<float, 10>  _S280 = dist_coeffs_2;
    bool _S281 = undistort_point_0(_S279, &_S280, int(12), &uv_4);
    if(!_S281)
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
    float2  _S282 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S271) / _S272;
    float2  uv_5 = _S282;
    FixedArray<float, 10>  _S283 = dist_coeffs_2;
    bool _S284 = undistort_point_0(_S282, &_S283, int(12), &uv_5);
    if(!_S284)
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
    float2  _S285;
    bool _S286;
    float2  _S287;
    bool _S288;
    float2  _S289;
    bool _S290;
    float2  _S291;
    bool _S292;
};

inline __device__ float s_primal_ctx_sin_0(float _S293)
{
    return (F32_sin((_S293)));
}

inline __device__ float s_primal_ctx_cos_0(float _S294)
{
    return (F32_cos((_S294)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S295 = make_float2 (0.0f);
    _s_diff_ctx_0->_S285 = _S295;
    _s_diff_ctx_0->_S286 = false;
    _s_diff_ctx_0->_S287 = _S295;
    _s_diff_ctx_0->_S288 = false;
    _s_diff_ctx_0->_S289 = _S295;
    _s_diff_ctx_0->_S290 = false;
    _s_diff_ctx_0->_S291 = _S295;
    _s_diff_ctx_0->_S292 = false;
    _s_diff_ctx_0->_S287 = _S295;
    _s_diff_ctx_0->_S288 = false;
    _s_diff_ctx_0->_S289 = _S295;
    _s_diff_ctx_0->_S290 = false;
    _s_diff_ctx_0->_S291 = _S295;
    _s_diff_ctx_0->_S292 = false;
    float3  _S296 = make_float3 (0.0f);
    float2  _S297 = float2 {intrins_2.z, intrins_2.w};
    float2  _S298 = float2 {intrins_2.x, intrins_2.y};
    float2  _S299 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S297) / _S298;
    float2  _S300 = _S299;
    bool _S301 = undistort_point_0(_S299, dist_coeffs_3, int(12), &_S300);
    _s_diff_ctx_0->_S285 = _S300;
    _s_diff_ctx_0->_S286 = _S301;
    float2  uv_6 = _S300;
    bool _S302 = !_S301;
    float3  normal_4;
    if(_S302)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S303 = !_S302;
    int _S304;
    FixedArray<float3 , 4>  points_3;
    if(_S303)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S305 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S305), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S305))).x, (uv_6 / make_float2 ((F32_max((_S305), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S305))).y, s_primal_ctx_cos_0(_S305));
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
        float3  _S306 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S307 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S297) / _S298;
        float2  _S308 = _S307;
        bool _S309 = undistort_point_0(_S307, dist_coeffs_3, int(12), &_S308);
        _s_diff_ctx_0->_S287 = _S308;
        _s_diff_ctx_0->_S288 = _S309;
        float2  uv_7 = _S308;
        bool _S310 = !_S309;
        if(_S310)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S311 = !_S310;
        if(_S311)
        {
            if(is_fisheye_2)
            {
                float _S312 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S312), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S312))).x, (uv_7 / make_float2 ((F32_max((_S312), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S312))).y, s_primal_ctx_cos_0(_S312));
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
            float3  _S313 = make_float3 (dpdepths_0.y) * raydir_12;
            _S304 = int(2);
            points_3[int(0)] = _S306;
            points_3[int(1)] = _S313;
            points_3[int(2)] = _S296;
            points_3[int(3)] = _S296;
        }
        else
        {
            _S304 = int(0);
            points_3[int(0)] = _S306;
            points_3[int(1)] = _S296;
            points_3[int(2)] = _S296;
            points_3[int(3)] = _S296;
        }
        bool _runFlag_0;
        if(_S304 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S303;
            _S304 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S314 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S297) / _S298;
            float2  _S315 = _S314;
            bool _S316 = undistort_point_0(_S314, dist_coeffs_3, int(12), &_S315);
            _s_diff_ctx_0->_S289 = _S315;
            _s_diff_ctx_0->_S290 = _S316;
            float2  uv_8 = _S315;
            if(!_S316)
            {
                float3  _S317 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S304 = int(0);
                normal_4 = _S317;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S318 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S318), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S318))).x, (uv_8 / make_float2 ((F32_max((_S318), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S318))).y, s_primal_ctx_cos_0(_S318));
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
                float2  _S319 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S297) / _S298;
                float2  _S320 = _S319;
                bool _S321 = undistort_point_0(_S319, dist_coeffs_3, int(12), &_S320);
                _s_diff_ctx_0->_S291 = _S320;
                _s_diff_ctx_0->_S292 = _S321;
                float2  uv_9 = _S320;
                bool _S322 = !_S321;
                if(_S322)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S323 = !_S322;
                int _S324;
                if(_S323)
                {
                    if(is_fisheye_2)
                    {
                        float _S325 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S325), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S325))).x, (uv_9 / make_float2 ((F32_max((_S325), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S325))).y, s_primal_ctx_cos_0(_S325));
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
                    _S324 = int(2);
                }
                else
                {
                    _S324 = int(0);
                }
                if(_S324 != int(2))
                {
                    _runFlag_0 = false;
                    _S304 = _S324;
                }
                if(_runFlag_0)
                {
                    _S304 = int(1);
                }
            }
        }
    }
    else
    {
        _S304 = int(0);
        points_3[int(0)] = _S296;
        points_3[int(1)] = _S296;
        points_3[int(2)] = _S296;
        points_3[int(3)] = _S296;
    }
    if(!(_S304 != int(1)))
    {
        float3  _S326 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S326, _S326)) != 0.0f)
        {
            normal_4 = _S326 / make_float3 (length_1(_S326));
        }
        else
        {
            normal_4 = _S326;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_6, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S327 = *dpdepths_1;
    float3  _S328 = make_float3 (0.0f);
    float2  _S329 = _s_diff_ctx_1->_S285;
    bool _S330 = !!_s_diff_ctx_1->_S286;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S331;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S332;
    if(_S330)
    {
        if(is_fisheye_3)
        {
            float _S333 = length_0(_S329);
            float3  raydir_25 = make_float3 ((_S329 / make_float2 ((F32_max((_S333), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S333))).x, (_S329 / make_float2 ((F32_max((_S333), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S333))).y, s_primal_ctx_cos_0(_S333));
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
            float3  raydir_26 = make_float3 (_S329.x, _S329.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S334 = make_float3 (_S327.primal_0.x) * raydir_21;
        float2  _S335 = _s_diff_ctx_1->_S287;
        bool _S336 = !!_s_diff_ctx_1->_S288;
        if(_S336)
        {
            if(is_fisheye_3)
            {
                float _S337 = length_0(_S335);
                float3  raydir_27 = make_float3 ((_S335 / make_float2 ((F32_max((_S337), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S337))).x, (_S335 / make_float2 ((F32_max((_S337), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S337))).y, s_primal_ctx_cos_0(_S337));
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
                float3  raydir_28 = make_float3 (_S335.x, _S335.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S338 = make_float3 (_S327.primal_0.y) * raydir_22;
            _S331 = int(2);
            points_4[int(0)] = _S334;
            points_4[int(1)] = _S338;
            points_4[int(2)] = _S328;
            points_4[int(3)] = _S328;
        }
        else
        {
            _S331 = int(0);
            points_4[int(0)] = _S334;
            points_4[int(1)] = _S328;
            points_4[int(2)] = _S328;
            points_4[int(3)] = _S328;
            raydir_22 = _S328;
        }
        if(_S331 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S330;
            _S331 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S339 = _s_diff_ctx_1->_S289;
            if(!_s_diff_ctx_1->_S290)
            {
                _runFlag_2 = false;
                _S331 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S340 = length_0(_S339);
                    float3  raydir_29 = make_float3 ((_S339 / make_float2 ((F32_max((_S340), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S340))).x, (_S339 / make_float2 ((F32_max((_S340), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S340))).y, s_primal_ctx_cos_0(_S340));
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
                    float3  raydir_30 = make_float3 (_S339.x, _S339.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S327.primal_0.z) * raydir_23;
                float2  _S341 = _s_diff_ctx_1->_S291;
                bool _S342 = !!_s_diff_ctx_1->_S292;
                int _S343;
                if(_S342)
                {
                    if(is_fisheye_3)
                    {
                        float _S344 = length_0(_S341);
                        float3  raydir_31 = make_float3 ((_S341 / make_float2 ((F32_max((_S344), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S344))).x, (_S341 / make_float2 ((F32_max((_S344), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S344))).y, s_primal_ctx_cos_0(_S344));
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
                        float3  raydir_32 = make_float3 (_S341.x, _S341.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S327.primal_0.w) * raydir_24;
                    _S343 = int(2);
                }
                else
                {
                    _S343 = int(0);
                    raydir_24 = _S328;
                }
                if(_S343 != int(2))
                {
                    _runFlag_3 = false;
                    _S331 = _S343;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S331 = int(1);
                }
                float3  _S345 = raydir_23;
                _runFlag_3 = _S342;
                raydir_23 = raydir_24;
                raydir_24 = _S345;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S328;
                raydir_24 = _S328;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S328;
            raydir_24 = _S328;
        }
        float3  _S346 = raydir_21;
        float3  _S347 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S332 = _S336;
        raydir_23 = _S347;
        raydir_24 = _S346;
    }
    else
    {
        _S331 = int(0);
        points_4[int(0)] = _S328;
        points_4[int(1)] = _S328;
        points_4[int(2)] = _S328;
        points_4[int(3)] = _S328;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S328;
        raydir_22 = _S328;
        _S332 = false;
        raydir_23 = _S328;
        raydir_24 = _S328;
    }
    bool _S348 = !(_S331 != int(1));
    float3  _S349;
    float3  _S350;
    float3  _S351;
    float3  _S352;
    float3  _S353;
    bool _S354;
    if(_S348)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S355 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S356 = s_primal_ctx_cross_0(dx_1, _S355);
        bool _S357 = (s_primal_ctx_dot_0(_S356, _S356)) != 0.0f;
        if(_S357)
        {
            float _S358 = length_1(_S356);
            float3  _S359 = make_float3 (_S358);
            _S349 = make_float3 (_S358 * _S358);
            _S350 = _S359;
        }
        else
        {
            _S349 = _S328;
            _S350 = _S328;
        }
        float3  _S360 = _S350;
        _S354 = _S357;
        _S350 = _S356;
        _S351 = _S360;
        _S352 = dx_1;
        _S353 = _S355;
    }
    else
    {
        _S354 = false;
        _S349 = _S328;
        _S350 = _S328;
        _S351 = _S328;
        _S352 = _S328;
        _S353 = _S328;
    }
    float4  _S361 = make_float4 (0.0f);
    if(_S348)
    {
        if(_S354)
        {
            float3  _S362 = _s_dOut_6 / _S349;
            float3  _S363 = _S350 * - _S362;
            float3  _S364 = _S351 * _S362;
            float _S365 = _S363.x + _S363.y + _S363.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S366;
            (&_S366)->primal_0 = _S350;
            (&_S366)->differential_0 = _S328;
            s_bwd_length_impl_0(&_S366, _S365);
            _S349 = _S364 + _S366.differential_0;
        }
        else
        {
            _S349 = _s_dOut_6;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S367;
        (&_S367)->primal_0 = _S350;
        (&_S367)->differential_0 = _S328;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S368;
        (&_S368)->primal_0 = _S350;
        (&_S368)->differential_0 = _S328;
        s_bwd_prop_dot_0(&_S367, &_S368, 0.0f);
        float3  _S369 = _S368.differential_0 + _S367.differential_0 + _S349;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S370;
        (&_S370)->primal_0 = _S352;
        (&_S370)->differential_0 = _S328;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S371;
        (&_S371)->primal_0 = _S353;
        (&_S371)->differential_0 = _S328;
        s_bwd_prop_cross_0(&_S370, &_S371, _S369);
        float3  s_diff_dy_T_1 = - _S371.differential_0;
        float3  _S372 = - s_diff_dy_T_1;
        float3  _S373 = - _S370.differential_0;
        FixedArray<float3 , 4>  _S374;
        _S374[int(0)] = _S328;
        _S374[int(1)] = _S328;
        _S374[int(2)] = _S328;
        _S374[int(3)] = _S328;
        _S374[int(2)] = _S372;
        _S374[int(3)] = s_diff_dy_T_1;
        _S374[int(0)] = _S373;
        _S374[int(1)] = _S370.differential_0;
        points_4[int(0)] = _S374[int(0)];
        points_4[int(1)] = _S374[int(1)];
        points_4[int(2)] = _S374[int(2)];
        points_4[int(3)] = _S374[int(3)];
    }
    else
    {
        points_4[int(0)] = _S328;
        points_4[int(1)] = _S328;
        points_4[int(2)] = _S328;
        points_4[int(3)] = _S328;
    }
    float4  _S375;
    if(_S330)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S376 = points_4;
                FixedArray<float3 , 4>  _S377 = points_4;
                FixedArray<float3 , 4>  _S378 = points_4;
                FixedArray<float3 , 4>  _S379 = points_4;
                if(_runFlag_3)
                {
                    float3  _S380 = raydir_21 * _S379[int(3)];
                    float _S381 = _S380.x + _S380.y + _S380.z;
                    float4  _S382 = _S361;
                    *&((&_S382)->w) = _S381;
                    points_4[int(0)] = _S376[int(0)];
                    points_4[int(1)] = _S377[int(1)];
                    points_4[int(2)] = _S378[int(2)];
                    points_4[int(3)] = _S328;
                    _S375 = _S382;
                }
                else
                {
                    points_4[int(0)] = _S376[int(0)];
                    points_4[int(1)] = _S377[int(1)];
                    points_4[int(2)] = _S378[int(2)];
                    points_4[int(3)] = _S379[int(3)];
                    _S375 = _S361;
                }
                float3  _S383 = raydir_22 * points_4[int(2)];
                float _S384 = _S383.x + _S383.y + _S383.z;
                FixedArray<float3 , 4>  _S385 = points_4;
                FixedArray<float3 , 4>  _S386 = points_4;
                float4  _S387 = _S361;
                *&((&_S387)->z) = _S384;
                float4  _S388 = _S375 + _S387;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S385[int(1)];
                points_4[int(2)] = _S328;
                points_4[int(3)] = _S386[int(3)];
                _S375 = _S388;
            }
            else
            {
                FixedArray<float3 , 4>  _S389 = points_4;
                FixedArray<float3 , 4>  _S390 = points_4;
                FixedArray<float3 , 4>  _S391 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S389[int(1)];
                points_4[int(2)] = _S390[int(2)];
                points_4[int(3)] = _S391[int(3)];
                _S375 = _S361;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S392 = points_4;
            FixedArray<float3 , 4>  _S393 = points_4;
            FixedArray<float3 , 4>  _S394 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S392[int(1)];
            points_4[int(2)] = _S393[int(2)];
            points_4[int(3)] = _S394[int(3)];
            _S375 = _S361;
        }
        if(_S332)
        {
            FixedArray<float3 , 4>  _S395 = points_4;
            float3  _S396 = raydir_23 * points_4[int(1)];
            float _S397 = _S396.x + _S396.y + _S396.z;
            float4  _S398 = _S361;
            *&((&_S398)->y) = _S397;
            float4  _S399 = _S375 + _S398;
            points_4[int(0)] = _S328;
            points_4[int(1)] = _S328;
            points_4[int(2)] = _S328;
            points_4[int(3)] = _S328;
            raydir_21 = _S395[int(0)];
            _S375 = _S399;
        }
        else
        {
            FixedArray<float3 , 4>  _S400 = points_4;
            FixedArray<float3 , 4>  _S401 = points_4;
            FixedArray<float3 , 4>  _S402 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S400[int(1)];
            points_4[int(2)] = _S401[int(2)];
            points_4[int(3)] = _S402[int(3)];
            raydir_21 = _S328;
        }
        float3  _S403 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S404 = _S403.x + _S403.y + _S403.z;
        float4  _S405 = _S361;
        *&((&_S405)->x) = _S404;
        _S375 = _S375 + _S405;
    }
    else
    {
        _S375 = _S361;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S375;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S406, float4  _S407, FixedArray<float, 10>  * _S408, bool _S409, bool _S410, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S411, float3  _S412)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S413;
    float3  _S414 = s_primal_ctx_depth_to_normal_0(_S406, _S407, _S408, _S409, _S410, (*_S411).primal_0, &_S413);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S415 = _S413;
    s_bwd_prop_depth_to_normal_0(_S406, _S407, _S408, _S409, _S410, _S411, _S412, &_S415);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S416 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S416;
    FixedArray<float, 10>  _S417 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S417, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S418 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S418;
    FixedArray<float, 10>  _S419 = dist_coeffs_6;
    bool _S420 = undistort_point_0(_S418, &_S419, int(12), &uv_10);
    if(!_S420)
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

