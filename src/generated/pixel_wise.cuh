#pragma once

#include "generated/slang.cuh"

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
    return clamp_0(rgb_0 + make_float3 (transmittance_2, transmittance_2, transmittance_2) * background_0, make_float3 (0.0f), make_float3 (1.0f));
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S20, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S21, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S22, float3  _S23)
{
    _d_clamp_vector_0(_S20, _S21, _S22, _S23);
    return;
}

inline __device__ void s_bwd_prop_blend_background_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_float_0 * dptransmittance_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpbackground_0, float3  _s_dOut_1)
{
    float3  _S24 = make_float3 ((*dptransmittance_1).primal_0, (*dptransmittance_1).primal_0, (*dptransmittance_1).primal_0);
    float3  _S25 = make_float3 (0.0f);
    float3  _S26 = make_float3 (1.0f);
    float3  _S27 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S28;
    (&_S28)->primal_0 = (*dprgb_0).primal_0 + _S24 * (*dpbackground_0).primal_0;
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

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_4, float dOut_5)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_4).primal_0.x * dOut_5;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_5).primal_0.x * dOut_5;
    *&((&x_d_result_0)->y) = (*dpy_4).primal_0.y * dOut_5;
    *&((&y_d_result_0)->y) = (*dpx_5).primal_0.y * dOut_5;
    *&((&x_d_result_0)->z) = (*dpy_4).primal_0.z * dOut_5;
    *&((&y_d_result_0)->z) = (*dpx_5).primal_0.z * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = x_d_result_0;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_7, float3  y_2)
{
    int i_3 = int(0);
    float result_3 = 0.0f;
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_4 = result_3 + _slang_vector_get_element(x_7, i_3) * _slang_vector_get_element(y_2, i_3);
        i_3 = i_3 + int(1);
        result_3 = result_4;
    }
    return result_3;
}

inline __device__ float dot_1(float2  x_8, float2  y_3)
{
    int i_4 = int(0);
    float result_5 = 0.0f;
    for(;;)
    {
        if(i_4 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_6 = result_5 + _slang_vector_get_element(x_8, i_4) * _slang_vector_get_element(y_3, i_4);
        i_4 = i_4 + int(1);
        result_5 = result_6;
    }
    return result_5;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_6)
{
    float _S152 = dOut_6.y;
    float _S153 = dOut_6.z;
    float _S154 = dOut_6.x;
    float _S155 = (*a_0).primal_0.z * _S152 + - (*a_0).primal_0.y * _S153;
    float _S156 = - (*a_0).primal_0.z * _S154 + (*a_0).primal_0.x * _S153;
    float _S157 = (*a_0).primal_0.y * _S154 + - (*a_0).primal_0.x * _S152;
    float3  _S158 = make_float3 (- (*b_0).primal_0.z * _S152 + (*b_0).primal_0.y * _S153, (*b_0).primal_0.z * _S154 + - (*b_0).primal_0.x * _S153, - (*b_0).primal_0.y * _S154 + (*b_0).primal_0.x * _S152);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S158;
    float3  _S159 = make_float3 (_S155, _S156, _S157);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S159;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S160 = left_2.y;
    float _S161 = right_2.z;
    float _S162 = left_2.z;
    float _S163 = right_2.y;
    float _S164 = right_2.x;
    float _S165 = left_2.x;
    return make_float3 (_S160 * _S161 - _S162 * _S163, _S162 * _S164 - _S165 * _S161, _S165 * _S163 - _S160 * _S164);
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_6, float dOut_7)
{
    float _S166 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_6).primal_0)))))) * dOut_7;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S166;
    return;
}

inline __device__ float length_0(float3  x_9)
{
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ float length_1(float2  x_10)
{
    return (F32_sqrt((dot_1(x_10, x_10))));
}

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  points_0)
{
    float3  _S167 = points_0[int(0)];
    bool _S168;
    if((dot_0(_S167, _S167)) == 0.0f)
    {
        _S168 = true;
    }
    else
    {
        float3  _S169 = points_0[int(1)];
        _S168 = (dot_0(_S169, _S169)) == 0.0f;
    }
    if(_S168)
    {
        _S168 = true;
    }
    else
    {
        float3  _S170 = points_0[int(2)];
        _S168 = (dot_0(_S170, _S170)) == 0.0f;
    }
    if(_S168)
    {
        _S168 = true;
    }
    else
    {
        float3  _S171 = points_0[int(3)];
        _S168 = (dot_0(_S171, _S171)) == 0.0f;
    }
    if(_S168)
    {
        return make_float3 (0.0f);
    }
    float3  normal_0 = cross_0(points_0[int(1)] - points_0[int(0)], - (points_0[int(3)] - points_0[int(2)]));
    float3  normal_1;
    if((dot_0(normal_0, normal_0)) != 0.0f)
    {
        normal_1 = normal_0 / make_float3 (length_0(normal_0));
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

inline __device__ float s_primal_ctx_dot_0(float3  _S172, float3  _S173)
{
    return dot_0(_S172, _S173);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S174, float3  _S175)
{
    return cross_0(_S174, _S175);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S176, float _S177)
{
    _d_sqrt_0(_S176, _S177);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_4)
{
    float _S178 = (*dpx_7).primal_0.x;
    float _S179 = (*dpx_7).primal_0.y;
    float _S180 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S181;
    (&_S181)->primal_0 = _S178 * _S178 + _S179 * _S179 + _S180 * _S180;
    (&_S181)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S181, _s_dOut_4);
    float _S182 = (*dpx_7).primal_0.z * _S181.differential_0;
    float _S183 = _S182 + _S182;
    float _S184 = (*dpx_7).primal_0.y * _S181.differential_0;
    float _S185 = _S184 + _S184;
    float _S186 = (*dpx_7).primal_0.x * _S181.differential_0;
    float _S187 = _S186 + _S186;
    float3  _S188 = make_float3 (0.0f);
    *&((&_S188)->z) = _S183;
    *&((&_S188)->y) = _S185;
    *&((&_S188)->x) = _S187;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S188;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S189, float _S190)
{
    s_bwd_prop_length_impl_0(_S189, _S190);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S191, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S192, float _S193)
{
    _d_dot_0(_S191, _S192, _S193);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S194, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S195, float3  _S196)
{
    _d_cross_0(_S194, _S195, _S196);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_5)
{
    FixedArray<float3 , 4>  _S197 = dppoints_0->primal_0;
    float3  _S198 = make_float3 (0.0f);
    float3  _S199 = dppoints_0->primal_0[int(0)];
    bool _S200 = (s_primal_ctx_dot_0(_S199, _S199)) == 0.0f;
    bool _S201;
    float3  _S202;
    if(_S200)
    {
        _S201 = true;
        _S202 = _S198;
    }
    else
    {
        float3  _S203 = _S197[int(1)];
        _S201 = (s_primal_ctx_dot_0(_S203, _S203)) == 0.0f;
        _S202 = _S197[int(1)];
    }
    bool _S204;
    float3  _S205;
    if(_S201)
    {
        _S204 = true;
        _S205 = _S198;
    }
    else
    {
        float3  _S206 = _S197[int(2)];
        _S204 = (s_primal_ctx_dot_0(_S206, _S206)) == 0.0f;
        _S205 = _S197[int(2)];
    }
    bool _S207;
    float3  _S208;
    if(_S204)
    {
        _S207 = true;
        _S208 = _S198;
    }
    else
    {
        float3  _S209 = _S197[int(3)];
        _S207 = (s_primal_ctx_dot_0(_S209, _S209)) == 0.0f;
        _S208 = _S197[int(3)];
    }
    bool _S210 = !_S207;
    float3  _S211;
    float3  _S212;
    float3  _S213;
    float3  _S214;
    float3  _S215;
    if(_S210)
    {
        float3  dx_0 = _S197[int(1)] - _S197[int(0)];
        float3  _S216 = - (_S197[int(3)] - _S197[int(2)]);
        float3  _S217 = s_primal_ctx_cross_0(dx_0, _S216);
        bool _S218 = (s_primal_ctx_dot_0(_S217, _S217)) != 0.0f;
        if(_S218)
        {
            float _S219 = length_0(_S217);
            float3  _S220 = make_float3 (_S219);
            _S211 = make_float3 (_S219 * _S219);
            _S212 = _S220;
        }
        else
        {
            _S211 = _S198;
            _S212 = _S198;
        }
        float3  _S221 = _S212;
        _S207 = _S218;
        _S212 = _S217;
        _S213 = _S221;
        _S214 = dx_0;
        _S215 = _S216;
    }
    else
    {
        _S207 = false;
        _S211 = _S198;
        _S212 = _S198;
        _S213 = _S198;
        _S214 = _S198;
        _S215 = _S198;
    }
    FixedArray<float3 , 4>  _S222;
    if(_S210)
    {
        if(_S207)
        {
            float3  _S223 = _s_dOut_5 / _S211;
            float3  _S224 = _S212 * - _S223;
            float3  _S225 = _S213 * _S223;
            float _S226 = _S224.x + _S224.y + _S224.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S227;
            (&_S227)->primal_0 = _S212;
            (&_S227)->differential_0 = _S198;
            s_bwd_length_impl_0(&_S227, _S226);
            _S211 = _S225 + _S227.differential_0;
        }
        else
        {
            _S211 = _s_dOut_5;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S228;
        (&_S228)->primal_0 = _S212;
        (&_S228)->differential_0 = _S198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S229;
        (&_S229)->primal_0 = _S212;
        (&_S229)->differential_0 = _S198;
        s_bwd_prop_dot_0(&_S228, &_S229, 0.0f);
        float3  _S230 = _S229.differential_0 + _S228.differential_0 + _S211;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S231;
        (&_S231)->primal_0 = _S214;
        (&_S231)->differential_0 = _S198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S232;
        (&_S232)->primal_0 = _S215;
        (&_S232)->differential_0 = _S198;
        s_bwd_prop_cross_0(&_S231, &_S232, _S230);
        float3  s_diff_dy_T_0 = - _S232.differential_0;
        float3  _S233 = - s_diff_dy_T_0;
        float3  _S234 = - _S231.differential_0;
        FixedArray<float3 , 4>  _S235;
        _S235[int(0)] = _S198;
        _S235[int(1)] = _S198;
        _S235[int(2)] = _S198;
        _S235[int(3)] = _S198;
        _S235[int(2)] = _S233;
        _S235[int(3)] = s_diff_dy_T_0;
        _S235[int(1)] = _S231.differential_0;
        _S222[int(0)] = _S235[int(0)];
        _S222[int(1)] = _S235[int(1)];
        _S222[int(2)] = _S235[int(2)];
        _S222[int(3)] = _S235[int(3)];
        _S211 = _S234;
    }
    else
    {
        _S222[int(0)] = _S198;
        _S222[int(1)] = _S198;
        _S222[int(2)] = _S198;
        _S222[int(3)] = _S198;
        _S211 = _S198;
    }
    if(_S204)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S236;
        (&_S236)->primal_0 = _S208;
        (&_S236)->differential_0 = _S198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S237;
        (&_S237)->primal_0 = _S208;
        (&_S237)->differential_0 = _S198;
        s_bwd_prop_dot_0(&_S236, &_S237, 0.0f);
        float3  _S238 = _S237.differential_0 + _S236.differential_0;
        FixedArray<float3 , 4>  _S239;
        _S239[int(0)] = _S198;
        _S239[int(1)] = _S198;
        _S239[int(2)] = _S198;
        _S239[int(3)] = _S198;
        _S239[int(3)] = _S238;
        float3  _S240 = _S222[int(1)] + _S239[int(1)];
        float3  _S241 = _S222[int(2)] + _S239[int(2)];
        float3  _S242 = _S222[int(3)] + _S239[int(3)];
        _S222[int(0)] = _S222[int(0)] + _S239[int(0)];
        _S222[int(1)] = _S240;
        _S222[int(2)] = _S241;
        _S222[int(3)] = _S242;
    }
    if(_S201)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S243;
        (&_S243)->primal_0 = _S205;
        (&_S243)->differential_0 = _S198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S244;
        (&_S244)->primal_0 = _S205;
        (&_S244)->differential_0 = _S198;
        s_bwd_prop_dot_0(&_S243, &_S244, 0.0f);
        float3  _S245 = _S244.differential_0 + _S243.differential_0;
        FixedArray<float3 , 4>  _S246;
        _S246[int(0)] = _S198;
        _S246[int(1)] = _S198;
        _S246[int(2)] = _S198;
        _S246[int(3)] = _S198;
        _S246[int(2)] = _S245;
        float3  _S247 = _S222[int(1)] + _S246[int(1)];
        float3  _S248 = _S222[int(2)] + _S246[int(2)];
        float3  _S249 = _S222[int(3)] + _S246[int(3)];
        _S222[int(0)] = _S222[int(0)] + _S246[int(0)];
        _S222[int(1)] = _S247;
        _S222[int(2)] = _S248;
        _S222[int(3)] = _S249;
    }
    if(_S200)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S250;
        (&_S250)->primal_0 = _S202;
        (&_S250)->differential_0 = _S198;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S251;
        (&_S251)->primal_0 = _S202;
        (&_S251)->differential_0 = _S198;
        s_bwd_prop_dot_0(&_S250, &_S251, 0.0f);
        float3  _S252 = _S251.differential_0 + _S250.differential_0;
        FixedArray<float3 , 4>  _S253;
        _S253[int(0)] = _S198;
        _S253[int(1)] = _S198;
        _S253[int(2)] = _S198;
        _S253[int(3)] = _S198;
        _S253[int(1)] = _S252;
        float3  _S254 = _S222[int(1)] + _S253[int(1)];
        float3  _S255 = _S222[int(2)] + _S253[int(2)];
        float3  _S256 = _S222[int(3)] + _S253[int(3)];
        _S222[int(0)] = _S222[int(0)] + _S253[int(0)];
        _S222[int(1)] = _S254;
        _S222[int(2)] = _S255;
        _S222[int(3)] = _S256;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S257;
    (&_S257)->primal_0 = _S197[int(0)];
    (&_S257)->differential_0 = _S198;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S258;
    (&_S258)->primal_0 = _S197[int(0)];
    (&_S258)->differential_0 = _S198;
    s_bwd_prop_dot_0(&_S257, &_S258, 0.0f);
    float3  _S259 = _S258.differential_0 + _S257.differential_0 + _S211;
    FixedArray<float3 , 4>  _S260;
    _S260[int(0)] = _S198;
    _S260[int(1)] = _S198;
    _S260[int(2)] = _S198;
    _S260[int(3)] = _S198;
    _S260[int(0)] = _S259;
    FixedArray<float3 , 4>  _S261 = {
        _S222[int(0)] + _S260[int(0)], _S222[int(1)] + _S260[int(1)], _S222[int(2)] + _S260[int(2)], _S222[int(3)] + _S260[int(3)]
    };
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S261;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S262, float3  _S263)
{
    s_bwd_prop_points_to_normal_0(_S262, _S263);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S264 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S264;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_11)
{
    Matrix<float, 2, 2>  result_7;
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
            *_slang_vector_get_element_ptr(((&result_7)->rows + (r_0)), c_0) = _slang_vector_get_element(x_11.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_7;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

inline __device__ bool undistort_point_0(float2  uv_0, FixedArray<float, 1>  * dist_coeffs_0, int maxiter_0, float2  * uv_undist_0)
{
    *uv_undist_0 = uv_0;
    return true;
}

inline __device__ float2  DistOpenCV_distort_0(float2  uv_1, FixedArray<float, 4>  * coeffs_0)
{
    float u_0 = uv_1.x;
    float v_0 = uv_1.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    return uv_1 * make_float2 (1.0f + r2_0 * ((*coeffs_0)[int(0)] + r2_0 * (*coeffs_0)[int(1)])) + make_float2 (2.0f * (*coeffs_0)[int(2)] * u_0 * v_0 + (*coeffs_0)[int(3)] * (r2_0 + 2.0f * u_0 * u_0), 2.0f * (*coeffs_0)[int(3)] * u_0 * v_0 + (*coeffs_0)[int(2)] * (r2_0 + 2.0f * v_0 * v_0));
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistOpenCV_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_0, FixedArray<float, 4>  * coeffs_1)
{
    float u_1 = dpuv_0->primal_0.x;
    float s_diff_u_0 = dpuv_0->differential_0.x;
    float v_1 = dpuv_0->primal_0.y;
    float s_diff_v_0 = dpuv_0->differential_0.y;
    float _S265 = s_diff_u_0 * u_1;
    float _S266 = s_diff_v_0 * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_0 = _S265 + _S265 + (_S266 + _S266);
    float _S267 = (*coeffs_1)[int(0)] + r2_1 * (*coeffs_1)[int(1)];
    float radial_0 = 1.0f + r2_1 * _S267;
    float _S268 = 2.0f * (*coeffs_1)[int(2)];
    float _S269 = _S268 * u_1;
    float _S270 = 2.0f * u_1;
    float _S271 = 2.0f * (*coeffs_1)[int(3)];
    float _S272 = _S271 * u_1;
    float _S273 = 2.0f * v_1;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S274 = { dpuv_0->primal_0 * make_float2 (radial_0) + make_float2 (_S269 * v_1 + (*coeffs_1)[int(3)] * (r2_1 + _S270 * u_1), _S272 * v_1 + (*coeffs_1)[int(2)] * (r2_1 + _S273 * v_1)), dpuv_0->differential_0 * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S267 + s_diff_r2_0 * (*coeffs_1)[int(1)] * r2_1) * dpuv_0->primal_0 + make_float2 (s_diff_u_0 * _S268 * v_1 + s_diff_v_0 * _S269 + (s_diff_r2_0 + (s_diff_u_0 * 2.0f * u_1 + s_diff_u_0 * _S270)) * (*coeffs_1)[int(3)], s_diff_u_0 * _S271 * v_1 + s_diff_v_0 * _S272 + (s_diff_r2_0 + (s_diff_v_0 * 2.0f * v_1 + s_diff_v_0 * _S273)) * (*coeffs_1)[int(2)]) };
    return _S274;
}

inline __device__ bool undistort_point_1(float2  uv_2, FixedArray<float, 4>  * dist_coeffs_1, int maxiter_1, float2  * uv_undist_1)
{
    int i_5 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_5 < maxiter_1)
        {
        }
        else
        {
            break;
        }
        float2  _S275 = DistOpenCV_distort_0(q_0, dist_coeffs_1);
        float2  r_1 = _S275 - uv_2;
        float2  _S276 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S277;
        (&_S277)->primal_0 = q_0;
        (&_S277)->differential_0 = _S276;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S278 = s_fwd_DistOpenCV_distort_0(&_S277, dist_coeffs_1);
        float2  _S279 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S280;
        (&_S280)->primal_0 = q_0;
        (&_S280)->differential_0 = _S279;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S281 = s_fwd_DistOpenCV_distort_0(&_S280, dist_coeffs_1);
        Matrix<float, 2, 2>  _S282 = transpose_0(makeMatrix<float, 2, 2> (_S278.differential_0, _S281.differential_0));
        float inv_det_0 = 1.0f / (_S282.rows[int(0)].x * _S282.rows[int(1)].y - _S282.rows[int(0)].y * _S282.rows[int(1)].x);
        float _S283 = r_1.x;
        float _S284 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S283 * _S282.rows[int(1)].y - _S284 * _S282.rows[int(0)].y) * inv_det_0, (- _S283 * _S282.rows[int(1)].x + _S284 * _S282.rows[int(0)].x) * inv_det_0);
        i_5 = i_5 + int(1);
        q_0 = q_1;
    }
    *uv_undist_1 = q_0;
    float2  _S285 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S286;
    (&_S286)->primal_0 = q_0;
    (&_S286)->differential_0 = _S285;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S287 = s_fwd_DistOpenCV_distort_0(&_S286, dist_coeffs_1);
    float2  _S288 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S289;
    (&_S289)->primal_0 = q_0;
    (&_S289)->differential_0 = _S288;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S290 = s_fwd_DistOpenCV_distort_0(&_S289, dist_coeffs_1);
    Matrix<float, 2, 2>  _S291 = transpose_0(makeMatrix<float, 2, 2> (_S287.differential_0, _S290.differential_0));
    float _S292 = (F32_min((determinant_0(_S291)), ((F32_min((_S291.rows[int(0)].x), (_S291.rows[int(1)].y))))));
    bool _S293;
    if(_S292 > 0.25f)
    {
        _S293 = _S292 < 4.0f;
    }
    else
    {
        _S293 = false;
    }
    if(_S293)
    {
        float2  _S294 = DistOpenCV_distort_0(q_0, dist_coeffs_1);
        _S293 = (dot_1(q_0, _S294)) >= 0.0f;
    }
    else
    {
        _S293 = false;
    }
    if(_S293)
    {
        float2  _S295 = DistOpenCV_distort_0(*uv_undist_1, dist_coeffs_1);
        _S293 = (length_1(_S295 - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S293 = false;
    }
    return _S293;
}

inline __device__ float2  DistThinPrism_distort_0(float2  uv_3, FixedArray<float, 8>  * coeffs_2)
{
    float u_2 = uv_3.x;
    float v_2 = uv_3.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    return uv_3 * make_float2 (1.0f + r2_2 * ((*coeffs_2)[int(0)] + r2_2 * ((*coeffs_2)[int(1)] + r2_2 * ((*coeffs_2)[int(2)] + r2_2 * (*coeffs_2)[int(3)])))) + make_float2 (2.0f * (*coeffs_2)[int(4)] * u_2 * v_2 + (*coeffs_2)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*coeffs_2)[int(6)] * r2_2, 2.0f * (*coeffs_2)[int(5)] * u_2 * v_2 + (*coeffs_2)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*coeffs_2)[int(7)] * r2_2);
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistThinPrism_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_1, FixedArray<float, 8>  * coeffs_3)
{
    float u_3 = dpuv_1->primal_0.x;
    float s_diff_u_1 = dpuv_1->differential_0.x;
    float v_3 = dpuv_1->primal_0.y;
    float s_diff_v_1 = dpuv_1->differential_0.y;
    float _S296 = s_diff_u_1 * u_3;
    float _S297 = s_diff_v_1 * v_3;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float s_diff_r2_1 = _S296 + _S296 + (_S297 + _S297);
    float _S298 = (*coeffs_3)[int(2)] + r2_3 * (*coeffs_3)[int(3)];
    float _S299 = (*coeffs_3)[int(1)] + r2_3 * _S298;
    float _S300 = (*coeffs_3)[int(0)] + r2_3 * _S299;
    float radial_1 = 1.0f + r2_3 * _S300;
    float _S301 = 2.0f * (*coeffs_3)[int(4)];
    float _S302 = _S301 * u_3;
    float _S303 = 2.0f * u_3;
    float _S304 = 2.0f * (*coeffs_3)[int(5)];
    float _S305 = _S304 * u_3;
    float _S306 = 2.0f * v_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S307 = { dpuv_1->primal_0 * make_float2 (radial_1) + make_float2 (_S302 * v_3 + (*coeffs_3)[int(5)] * (r2_3 + _S303 * u_3) + (*coeffs_3)[int(6)] * r2_3, _S305 * v_3 + (*coeffs_3)[int(4)] * (r2_3 + _S306 * v_3) + (*coeffs_3)[int(7)] * r2_3), dpuv_1->differential_0 * make_float2 (radial_1) + make_float2 (s_diff_r2_1 * _S300 + (s_diff_r2_1 * _S299 + (s_diff_r2_1 * _S298 + s_diff_r2_1 * (*coeffs_3)[int(3)] * r2_3) * r2_3) * r2_3) * dpuv_1->primal_0 + make_float2 (s_diff_u_1 * _S301 * v_3 + s_diff_v_1 * _S302 + (s_diff_r2_1 + (s_diff_u_1 * 2.0f * u_3 + s_diff_u_1 * _S303)) * (*coeffs_3)[int(5)] + s_diff_r2_1 * (*coeffs_3)[int(6)], s_diff_u_1 * _S304 * v_3 + s_diff_v_1 * _S305 + (s_diff_r2_1 + (s_diff_v_1 * 2.0f * v_3 + s_diff_v_1 * _S306)) * (*coeffs_3)[int(4)] + s_diff_r2_1 * (*coeffs_3)[int(7)]) };
    return _S307;
}

inline __device__ bool undistort_point_2(float2  uv_4, FixedArray<float, 8>  * dist_coeffs_2, int maxiter_2, float2  * uv_undist_2)
{
    int i_6 = int(0);
    float2  q_2 = uv_4;
    for(;;)
    {
        if(i_6 < maxiter_2)
        {
        }
        else
        {
            break;
        }
        float2  _S308 = DistThinPrism_distort_0(q_2, dist_coeffs_2);
        float2  r_2 = _S308 - uv_4;
        float2  _S309 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S310;
        (&_S310)->primal_0 = q_2;
        (&_S310)->differential_0 = _S309;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S311 = s_fwd_DistThinPrism_distort_0(&_S310, dist_coeffs_2);
        float2  _S312 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S313;
        (&_S313)->primal_0 = q_2;
        (&_S313)->differential_0 = _S312;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S314 = s_fwd_DistThinPrism_distort_0(&_S313, dist_coeffs_2);
        Matrix<float, 2, 2>  _S315 = transpose_0(makeMatrix<float, 2, 2> (_S311.differential_0, _S314.differential_0));
        float inv_det_1 = 1.0f / (_S315.rows[int(0)].x * _S315.rows[int(1)].y - _S315.rows[int(0)].y * _S315.rows[int(1)].x);
        float _S316 = r_2.x;
        float _S317 = r_2.y;
        float2  q_3 = q_2 - make_float2 ((_S316 * _S315.rows[int(1)].y - _S317 * _S315.rows[int(0)].y) * inv_det_1, (- _S316 * _S315.rows[int(1)].x + _S317 * _S315.rows[int(0)].x) * inv_det_1);
        i_6 = i_6 + int(1);
        q_2 = q_3;
    }
    *uv_undist_2 = q_2;
    float2  _S318 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S319;
    (&_S319)->primal_0 = q_2;
    (&_S319)->differential_0 = _S318;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S320 = s_fwd_DistThinPrism_distort_0(&_S319, dist_coeffs_2);
    float2  _S321 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S322;
    (&_S322)->primal_0 = q_2;
    (&_S322)->differential_0 = _S321;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S323 = s_fwd_DistThinPrism_distort_0(&_S322, dist_coeffs_2);
    Matrix<float, 2, 2>  _S324 = transpose_0(makeMatrix<float, 2, 2> (_S320.differential_0, _S323.differential_0));
    float _S325 = (F32_min((determinant_0(_S324)), ((F32_min((_S324.rows[int(0)].x), (_S324.rows[int(1)].y))))));
    bool _S326;
    if(_S325 > 0.25f)
    {
        _S326 = _S325 < 4.0f;
    }
    else
    {
        _S326 = false;
    }
    if(_S326)
    {
        float2  _S327 = DistThinPrism_distort_0(q_2, dist_coeffs_2);
        _S326 = (dot_1(q_2, _S327)) >= 0.0f;
    }
    else
    {
        _S326 = false;
    }
    if(_S326)
    {
        float2  _S328 = DistThinPrism_distort_0(*uv_undist_2, dist_coeffs_2);
        _S326 = (length_1(_S328 - uv_4)) < 0.00999999977648258f;
    }
    else
    {
        _S326 = false;
    }
    return _S326;
}

inline __device__ float2  DistRational_distort_0(float2  uv_5, FixedArray<float, 8>  * coeffs_4)
{
    float u_4 = uv_5.x;
    float v_4 = uv_5.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    return uv_5 * make_float2 ((1.0f + r2_4 * ((*coeffs_4)[int(0)] + r2_4 * ((*coeffs_4)[int(1)] + r2_4 * (*coeffs_4)[int(2)]))) / (1.0f + r2_4 * ((*coeffs_4)[int(3)] + r2_4 * ((*coeffs_4)[int(4)] + r2_4 * (*coeffs_4)[int(5)])))) + make_float2 (2.0f * (*coeffs_4)[int(6)] * u_4 * v_4 + (*coeffs_4)[int(7)] * (r2_4 + 2.0f * u_4 * u_4), 2.0f * (*coeffs_4)[int(7)] * u_4 * v_4 + (*coeffs_4)[int(6)] * (r2_4 + 2.0f * v_4 * v_4));
}

inline __device__ DiffPair_vectorx3Cfloatx2C2x3E_0 s_fwd_DistRational_distort_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpuv_2, FixedArray<float, 8>  * coeffs_5)
{
    float u_5 = dpuv_2->primal_0.x;
    float s_diff_u_2 = dpuv_2->differential_0.x;
    float v_5 = dpuv_2->primal_0.y;
    float s_diff_v_2 = dpuv_2->differential_0.y;
    float _S329 = s_diff_u_2 * u_5;
    float _S330 = s_diff_v_2 * v_5;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float s_diff_r2_2 = _S329 + _S329 + (_S330 + _S330);
    float _S331 = (*coeffs_5)[int(1)] + r2_5 * (*coeffs_5)[int(2)];
    float _S332 = (*coeffs_5)[int(0)] + r2_5 * _S331;
    float _S333 = 1.0f + r2_5 * _S332;
    float _S334 = (*coeffs_5)[int(4)] + r2_5 * (*coeffs_5)[int(5)];
    float _S335 = (*coeffs_5)[int(3)] + r2_5 * _S334;
    float _S336 = 1.0f + r2_5 * _S335;
    float radial_2 = _S333 / _S336;
    float _S337 = 2.0f * (*coeffs_5)[int(6)];
    float _S338 = _S337 * u_5;
    float _S339 = 2.0f * u_5;
    float _S340 = 2.0f * (*coeffs_5)[int(7)];
    float _S341 = _S340 * u_5;
    float _S342 = 2.0f * v_5;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S343 = { dpuv_2->primal_0 * make_float2 (radial_2) + make_float2 (_S338 * v_5 + (*coeffs_5)[int(7)] * (r2_5 + _S339 * u_5), _S341 * v_5 + (*coeffs_5)[int(6)] * (r2_5 + _S342 * v_5)), dpuv_2->differential_0 * make_float2 (radial_2) + make_float2 (((s_diff_r2_2 * _S332 + (s_diff_r2_2 * _S331 + s_diff_r2_2 * (*coeffs_5)[int(2)] * r2_5) * r2_5) * _S336 - _S333 * (s_diff_r2_2 * _S335 + (s_diff_r2_2 * _S334 + s_diff_r2_2 * (*coeffs_5)[int(5)] * r2_5) * r2_5)) / (_S336 * _S336)) * dpuv_2->primal_0 + make_float2 (s_diff_u_2 * _S337 * v_5 + s_diff_v_2 * _S338 + (s_diff_r2_2 + (s_diff_u_2 * 2.0f * u_5 + s_diff_u_2 * _S339)) * (*coeffs_5)[int(7)], s_diff_u_2 * _S340 * v_5 + s_diff_v_2 * _S341 + (s_diff_r2_2 + (s_diff_v_2 * 2.0f * v_5 + s_diff_v_2 * _S342)) * (*coeffs_5)[int(6)]) };
    return _S343;
}

inline __device__ bool undistort_point_3(float2  uv_6, FixedArray<float, 8>  * dist_coeffs_3, int maxiter_3, float2  * uv_undist_3)
{
    int i_7 = int(0);
    float2  q_4 = uv_6;
    for(;;)
    {
        if(i_7 < maxiter_3)
        {
        }
        else
        {
            break;
        }
        float2  _S344 = DistRational_distort_0(q_4, dist_coeffs_3);
        float2  r_3 = _S344 - uv_6;
        float2  _S345 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S346;
        (&_S346)->primal_0 = q_4;
        (&_S346)->differential_0 = _S345;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S347 = s_fwd_DistRational_distort_0(&_S346, dist_coeffs_3);
        float2  _S348 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S349;
        (&_S349)->primal_0 = q_4;
        (&_S349)->differential_0 = _S348;
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S350 = s_fwd_DistRational_distort_0(&_S349, dist_coeffs_3);
        Matrix<float, 2, 2>  _S351 = transpose_0(makeMatrix<float, 2, 2> (_S347.differential_0, _S350.differential_0));
        float inv_det_2 = 1.0f / (_S351.rows[int(0)].x * _S351.rows[int(1)].y - _S351.rows[int(0)].y * _S351.rows[int(1)].x);
        float _S352 = r_3.x;
        float _S353 = r_3.y;
        float2  q_5 = q_4 - make_float2 ((_S352 * _S351.rows[int(1)].y - _S353 * _S351.rows[int(0)].y) * inv_det_2, (- _S352 * _S351.rows[int(1)].x + _S353 * _S351.rows[int(0)].x) * inv_det_2);
        i_7 = i_7 + int(1);
        q_4 = q_5;
    }
    *uv_undist_3 = q_4;
    float2  _S354 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S355;
    (&_S355)->primal_0 = q_4;
    (&_S355)->differential_0 = _S354;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S356 = s_fwd_DistRational_distort_0(&_S355, dist_coeffs_3);
    float2  _S357 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S358;
    (&_S358)->primal_0 = q_4;
    (&_S358)->differential_0 = _S357;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S359 = s_fwd_DistRational_distort_0(&_S358, dist_coeffs_3);
    Matrix<float, 2, 2>  _S360 = transpose_0(makeMatrix<float, 2, 2> (_S356.differential_0, _S359.differential_0));
    float _S361 = (F32_min((determinant_0(_S360)), ((F32_min((_S360.rows[int(0)].x), (_S360.rows[int(1)].y))))));
    bool _S362;
    if(_S361 > 0.25f)
    {
        _S362 = _S361 < 4.0f;
    }
    else
    {
        _S362 = false;
    }
    if(_S362)
    {
        float2  _S363 = DistRational_distort_0(q_4, dist_coeffs_3);
        _S362 = (dot_1(q_4, _S363)) >= 0.0f;
    }
    else
    {
        _S362 = false;
    }
    if(_S362)
    {
        float2  _S364 = DistRational_distort_0(*uv_undist_3, dist_coeffs_3);
        _S362 = (length_1(_S364 - uv_6)) < 0.00999999977648258f;
    }
    else
    {
        _S362 = false;
    }
    return _S362;
}

inline __device__ float3  normalize_0(float3  x_12)
{
    return x_12 / make_float3 (length_0(x_12));
}

inline __device__ float3  unproject_raydir_0(float2  uv_7, int camera_model_0, bool is_ray_depth_0)
{
    float3  raydir_0;
    bool is_unit_0;
    if(camera_model_0 == int(1))
    {
        float theta_0 = length_1(uv_7);
        float3  _S365 = make_float3 ((uv_7 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).x, (uv_7 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).y, (F32_cos((theta_0))));
        is_unit_0 = true;
        raydir_0 = _S365;
    }
    else
    {
        bool _S366 = camera_model_0 == int(2);
        if(_S366)
        {
            float r_4 = length_1(uv_7);
            raydir_0 = make_float3 ((uv_7 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_4 * r_4)))))))).x, (uv_7 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_4 * r_4)))))))).y, 1.0f - 0.5f * r_4 * r_4);
        }
        else
        {
            raydir_0 = make_float3 (uv_7.x, uv_7.y, 1.0f);
        }
        is_unit_0 = _S366;
    }
    if(is_ray_depth_0)
    {
        if(is_unit_0)
        {
        }
        else
        {
            raydir_0 = normalize_0(raydir_0);
        }
    }
    else
    {
        raydir_0 = raydir_0 / make_float3 (raydir_0.z);
    }
    return raydir_0;
}

inline __device__ float3  generate_ray_d2n_none(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 1>  dist_coeffs_4, int camera_model_1, bool is_ray_depth_1)
{
    float3  _S367;
    for(;;)
    {
        float2  uv_8 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
        FixedArray<float, 1>  _S368 = dist_coeffs_4;
        float2  uv_u_0;
        bool _S369 = undistort_point_0(uv_8, &_S368, int(12), &uv_u_0);
        if(!_S369)
        {
            int3  _S370 = make_int3 (int(0));
            float3  _S371 = make_float3 ((float)_S370.x, (float)_S370.y, (float)_S370.z);
            _S367 = _S371;
            break;
        }
        _S367 = unproject_raydir_0(uv_u_0, camera_model_1, is_ray_depth_1);
        break;
    }
    return _S367;
}

inline __device__ float3  depth_to_point_none(float2  pix_pos_1, float4  intrins_1, FixedArray<float, 1>  dist_coeffs_5, int camera_model_2, bool is_ray_depth_2, float depth_2)
{
    float3  _S372;
    for(;;)
    {
        float2  uv_9 = (pix_pos_1 - float2 {intrins_1.z, intrins_1.w}) / float2 {intrins_1.x, intrins_1.y};
        FixedArray<float, 1>  _S373 = dist_coeffs_5;
        float2  uv_u_1;
        bool _S374 = undistort_point_0(uv_9, &_S373, int(12), &uv_u_1);
        if(!_S374)
        {
            _S372 = make_float3 (0.0f);
            break;
        }
        _S372 = make_float3 (depth_2) * unproject_raydir_0(uv_u_1, camera_model_2, is_ray_depth_2);
        break;
    }
    return _S372;
}

struct s_bwd_prop_depth_to_point_Intermediates_0
{
    float2  _S375;
    bool _S376;
};

inline __device__ float s_primal_ctx_sin_0(float _S377)
{
    return (F32_sin((_S377)));
}

inline __device__ float s_primal_ctx_cos_0(float _S378)
{
    return (F32_cos((_S378)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S379)
{
    return (F32_sqrt((_S379)));
}

inline __device__ float3  s_primal_ctx_unproject_raydir_0(float2  dpuv_3, int camera_model_3, bool is_ray_depth_3)
{
    float3  raydir_1;
    bool is_unit_1;
    if(camera_model_3 == int(1))
    {
        float _S380 = length_1(dpuv_3);
        float3  _S381 = make_float3 ((dpuv_3 / make_float2 ((F32_max((_S380), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S380))).x, (dpuv_3 / make_float2 ((F32_max((_S380), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S380))).y, s_primal_ctx_cos_0(_S380));
        is_unit_1 = true;
        raydir_1 = _S381;
    }
    else
    {
        bool _S382 = camera_model_3 == int(2);
        if(_S382)
        {
            float _S383 = length_1(dpuv_3);
            raydir_1 = make_float3 ((dpuv_3 * make_float2 (s_primal_ctx_sqrt_0((F32_max((0.0f), (1.0f - 0.25f * _S383 * _S383)))))).x, (dpuv_3 * make_float2 (s_primal_ctx_sqrt_0((F32_max((0.0f), (1.0f - 0.25f * _S383 * _S383)))))).y, 1.0f - 0.5f * _S383 * _S383);
        }
        else
        {
            raydir_1 = make_float3 (dpuv_3.x, dpuv_3.y, 1.0f);
        }
        is_unit_1 = _S382;
    }
    if(is_ray_depth_3)
    {
        if(is_unit_1)
        {
        }
        else
        {
            raydir_1 = normalize_0(raydir_1);
        }
    }
    else
    {
        raydir_1 = raydir_1 / make_float3 (raydir_1.z);
    }
    return raydir_1;
}

inline __device__ float depth_to_point_vjp_none(float2  pix_pos_2, float4  intrins_2, FixedArray<float, 1>  dist_coeffs_6, int camera_model_4, bool is_ray_depth_4, float depth_3, float3  v_point_0)
{
    float2  _S384 = make_float2 (0.0f);
    s_bwd_prop_depth_to_point_Intermediates_0 _S385;
    (&_S385)->_S375 = _S384;
    (&_S385)->_S376 = false;
    float2  uv_10 = (pix_pos_2 - float2 {intrins_2.z, intrins_2.w}) / float2 {intrins_2.x, intrins_2.y};
    float2  _S386 = _S384;
    FixedArray<float, 1>  _S387 = dist_coeffs_6;
    bool _S388 = undistort_point_0(uv_10, &_S387, int(12), &_S386);
    (&_S385)->_S375 = _S386;
    (&_S385)->_S376 = _S388;
    s_bwd_prop_depth_to_point_Intermediates_0 _S389 = _S385;
    float3  _S390 = make_float3 (0.0f);
    bool _S391 = !!_S385._S376;
    float3  _S392;
    if(_S391)
    {
        _S392 = s_primal_ctx_unproject_raydir_0(_S389._S375, camera_model_4, is_ray_depth_4);
    }
    else
    {
        _S392 = _S390;
    }
    if(_S391)
    {
        _S392 = _S392 * v_point_0;
    }
    else
    {
        _S392 = _S390;
    }
    return _S392.x + _S392.y + _S392.z;
}

inline __device__ float3  depth_to_normal_none(float2  pix_center_0, float4  intrins_3, FixedArray<float, 1>  dist_coeffs_7, int camera_model_5, bool is_ray_depth_5, float4  depths_0)
{
    float3  normal_2;
    for(;;)
    {
        bool _S393;
        if((depths_0.x) == 0.0f)
        {
            _S393 = true;
        }
        else
        {
            _S393 = (depths_0.y) == 0.0f;
        }
        if(_S393)
        {
            _S393 = true;
        }
        else
        {
            _S393 = (depths_0.z) == 0.0f;
        }
        if(_S393)
        {
            _S393 = true;
        }
        else
        {
            _S393 = (depths_0.w) == 0.0f;
        }
        if(_S393)
        {
            normal_2 = make_float3 (0.0f);
            break;
        }
        float3  * _S394;
        float3  * _S395;
        float3  * _S396;
        float3  * _S397;
        int _S398;
        FixedArray<float3 , 4>  points_2;
        for(;;)
        {
            float2  _S399 = float2 {intrins_3.z, intrins_3.w};
            float2  _S400 = float2 {intrins_3.x, intrins_3.y};
            float2  uv_11 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S399) / _S400;
            FixedArray<float, 1>  _S401 = dist_coeffs_7;
            float2  uv_u_2;
            bool _S402 = undistort_point_0(uv_11, &_S401, int(12), &uv_u_2);
            if(!_S402)
            {
                float3  _S403 = make_float3 (0.0f);
                _S398 = int(0);
                _S397 = nullptr;
                _S396 = nullptr;
                _S395 = nullptr;
                _S394 = nullptr;
                normal_2 = _S403;
                break;
            }
            points_2[int(0)] = make_float3 (depths_0.x) * unproject_raydir_0(uv_u_2, camera_model_5, is_ray_depth_5);
            for(;;)
            {
                float2  uv_12 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S399) / _S400;
                FixedArray<float, 1>  _S404 = dist_coeffs_7;
                float2  uv_u_3;
                bool _S405 = undistort_point_0(uv_12, &_S404, int(12), &uv_u_3);
                if(!_S405)
                {
                    float3  _S406 = make_float3 (0.0f);
                    _S398 = int(0);
                    _S397 = nullptr;
                    normal_2 = _S406;
                    break;
                }
                points_2[int(1)] = make_float3 (depths_0.y) * unproject_raydir_0(uv_u_3, camera_model_5, is_ray_depth_5);
                _S398 = int(2);
                _S397 = &points_2[int(1)];
                break;
            }
            if(_S398 != int(2))
            {
                _S396 = &points_2[int(0)];
                _S395 = nullptr;
                _S394 = nullptr;
                break;
            }
            float2  uv_13 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S399) / _S400;
            FixedArray<float, 1>  _S407 = dist_coeffs_7;
            float2  uv_u_4;
            bool _S408 = undistort_point_0(uv_13, &_S407, int(12), &uv_u_4);
            if(!_S408)
            {
                float3  _S409 = make_float3 (0.0f);
                _S398 = int(0);
                _S396 = &points_2[int(0)];
                _S395 = nullptr;
                _S394 = nullptr;
                normal_2 = _S409;
                break;
            }
            points_2[int(2)] = make_float3 (depths_0.z) * unproject_raydir_0(uv_u_4, camera_model_5, is_ray_depth_5);
            for(;;)
            {
                float2  uv_14 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S399) / _S400;
                FixedArray<float, 1>  _S410 = dist_coeffs_7;
                float2  uv_u_5;
                bool _S411 = undistort_point_0(uv_14, &_S410, int(12), &uv_u_5);
                if(!_S411)
                {
                    float3  _S412 = make_float3 (0.0f);
                    _S398 = int(0);
                    _S396 = nullptr;
                    normal_2 = _S412;
                    break;
                }
                points_2[int(3)] = make_float3 (depths_0.w) * unproject_raydir_0(uv_u_5, camera_model_5, is_ray_depth_5);
                _S398 = int(2);
                _S396 = &points_2[int(3)];
                break;
            }
            if(_S398 != int(2))
            {
                float3  * _S413 = _S396;
                _S396 = &points_2[int(0)];
                _S395 = _S413;
                _S394 = &points_2[int(2)];
                break;
            }
            float3  * _S414 = _S396;
            _S398 = int(1);
            _S396 = &points_2[int(0)];
            _S395 = _S414;
            _S394 = &points_2[int(2)];
            break;
        }
        if(_S398 != int(1))
        {
            break;
        }
        float3  normal_3 = cross_0(*_S397 - *_S396, - (*_S395 - *_S394));
        if((dot_0(normal_3, normal_3)) != 0.0f)
        {
            normal_2 = normal_3 / make_float3 (length_0(normal_3));
        }
        else
        {
            normal_2 = normal_3;
        }
        break;
    }
    return normal_2;
}

struct s_bwd_prop_depth_to_normal_Intermediates_0
{
    float2  _S415;
    bool _S416;
    float2  _S417;
    bool _S418;
    float2  _S419;
    bool _S420;
    float2  _S421;
    bool _S422;
};

inline __device__ void depth_to_normal_vjp_none(float2  pix_center_1, float4  intrins_4, FixedArray<float, 1>  dist_coeffs_8, int camera_model_6, bool is_ray_depth_6, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float2  _S423 = make_float2 (0.0f);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S424;
    (&_S424)->_S415 = _S423;
    (&_S424)->_S416 = false;
    (&_S424)->_S417 = _S423;
    (&_S424)->_S418 = false;
    (&_S424)->_S419 = _S423;
    (&_S424)->_S420 = false;
    (&_S424)->_S421 = _S423;
    (&_S424)->_S422 = false;
    (&_S424)->_S415 = _S423;
    (&_S424)->_S416 = false;
    (&_S424)->_S417 = _S423;
    (&_S424)->_S418 = false;
    (&_S424)->_S419 = _S423;
    (&_S424)->_S420 = false;
    (&_S424)->_S421 = _S423;
    (&_S424)->_S422 = false;
    bool _S425 = (depths_1.x) == 0.0f;
    bool _runFlag_0;
    if(_S425)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.y) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.z) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.w) == 0.0f;
    }
    int _S426;
    if(!_runFlag_0)
    {
        float2  _S427 = float2 {intrins_4.z, intrins_4.w};
        float2  _S428 = float2 {intrins_4.x, intrins_4.y};
        float2  uv_15 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S427) / _S428;
        float2  _S429 = _S423;
        FixedArray<float, 1>  _S430 = dist_coeffs_8;
        bool _S431 = undistort_point_0(uv_15, &_S430, int(12), &_S429);
        (&_S424)->_S415 = _S429;
        (&_S424)->_S416 = _S431;
        bool _S432 = !!_S431;
        if(_S432)
        {
            float2  uv_16 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S427) / _S428;
            float2  _S433 = _S423;
            FixedArray<float, 1>  _S434 = dist_coeffs_8;
            bool _S435 = undistort_point_0(uv_16, &_S434, int(12), &_S433);
            (&_S424)->_S417 = _S433;
            (&_S424)->_S418 = _S435;
            if(!!_S435)
            {
                _S426 = int(2);
            }
            else
            {
                _S426 = int(0);
            }
            if(_S426 != int(2))
            {
                _runFlag_0 = false;
            }
            else
            {
                _runFlag_0 = _S432;
            }
            if(_runFlag_0)
            {
                float2  uv_17 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S427) / _S428;
                float2  _S436 = _S423;
                FixedArray<float, 1>  _S437 = dist_coeffs_8;
                bool _S438 = undistort_point_0(uv_17, &_S437, int(12), &_S436);
                (&_S424)->_S419 = _S436;
                (&_S424)->_S420 = _S438;
                if(!_S438)
                {
                    _runFlag_0 = false;
                }
                if(_runFlag_0)
                {
                    float2  uv_18 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S427) / _S428;
                    float2  _S439 = _S423;
                    FixedArray<float, 1>  _S440 = dist_coeffs_8;
                    bool _S441 = undistort_point_0(uv_18, &_S440, int(12), &_S439);
                    (&_S424)->_S421 = _S439;
                    (&_S424)->_S422 = _S441;
                }
            }
        }
    }
    s_bwd_prop_depth_to_normal_Intermediates_0 _S442 = _S424;
    float3  _S443 = make_float3 (0.0f);
    if(_S425)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.y) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.z) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (depths_1.w) == 0.0f;
    }
    bool _S444 = !_runFlag_0;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _S445;
    bool _runFlag_3;
    bool _S446;
    bool _S447;
    FixedArray<float3 , 4>  points_3;
    float3  _S448;
    float3  _S449;
    float3  _S450;
    float3  _S451;
    float3  _S452;
    float3  _S453;
    float3  _S454;
    float3  _S455;
    float3  _S456;
    if(_S444)
    {
        bool _S457 = !!_S442._S416;
        if(_S457)
        {
            float3  _S458 = s_primal_ctx_unproject_raydir_0(_S442._S415, camera_model_6, is_ray_depth_6);
            float3  _S459 = make_float3 (depths_1.x) * _S458;
            bool _S460 = !!_S442._S418;
            if(_S460)
            {
                float3  _S461 = s_primal_ctx_unproject_raydir_0(_S442._S417, camera_model_6, is_ray_depth_6);
                float3  _S462 = make_float3 (depths_1.y) * _S461;
                _S426 = int(2);
                points_3[int(0)] = _S459;
                points_3[int(1)] = _S462;
                points_3[int(2)] = _S443;
                points_3[int(3)] = _S443;
                _S448 = _S461;
            }
            else
            {
                _S426 = int(0);
                points_3[int(0)] = _S459;
                points_3[int(1)] = _S443;
                points_3[int(2)] = _S443;
                points_3[int(3)] = _S443;
                _S448 = _S443;
            }
            if(_S426 != int(2))
            {
                _runFlag_0 = false;
            }
            else
            {
                _runFlag_0 = _S457;
                _S426 = int(0);
            }
            if(_runFlag_0)
            {
                if(!_S442._S420)
                {
                    _runFlag_1 = false;
                    _S426 = int(0);
                }
                else
                {
                    _runFlag_1 = _runFlag_0;
                }
                if(_runFlag_1)
                {
                    float3  _S463 = s_primal_ctx_unproject_raydir_0(_S442._S419, camera_model_6, is_ray_depth_6);
                    points_3[int(2)] = make_float3 (depths_1.z) * _S463;
                    bool _S464 = !!_S442._S422;
                    int _S465;
                    if(_S464)
                    {
                        float3  _S466 = s_primal_ctx_unproject_raydir_0(_S442._S421, camera_model_6, is_ray_depth_6);
                        points_3[int(3)] = make_float3 (depths_1.w) * _S466;
                        _S465 = int(2);
                        _S449 = _S466;
                    }
                    else
                    {
                        _S465 = int(0);
                        _S449 = _S443;
                    }
                    if(_S465 != int(2))
                    {
                        _runFlag_2 = false;
                        _S426 = _S465;
                    }
                    else
                    {
                        _runFlag_2 = _runFlag_1;
                    }
                    if(_runFlag_2)
                    {
                        _S426 = int(1);
                    }
                    _runFlag_2 = _S464;
                    _S450 = _S463;
                }
                else
                {
                    _runFlag_2 = false;
                    _S449 = _S443;
                    _S450 = _S443;
                }
            }
            else
            {
                _runFlag_1 = false;
                _runFlag_2 = false;
                _S449 = _S443;
                _S450 = _S443;
            }
            float3  _S467 = _S448;
            _S448 = _S449;
            _S449 = _S450;
            _S445 = _S460;
            _S450 = _S467;
            _S451 = _S458;
        }
        else
        {
            _S426 = int(0);
            points_3[int(0)] = _S443;
            points_3[int(1)] = _S443;
            points_3[int(2)] = _S443;
            points_3[int(3)] = _S443;
            _runFlag_0 = false;
            _runFlag_1 = false;
            _runFlag_2 = false;
            _S448 = _S443;
            _S449 = _S443;
            _S445 = false;
            _S450 = _S443;
            _S451 = _S443;
        }
        if(_S426 != int(1))
        {
            _runFlag_3 = false;
        }
        else
        {
            _runFlag_3 = _S444;
        }
        if(_runFlag_3)
        {
            float3  dx_1 = points_3[int(1)] - points_3[int(0)];
            float3  _S468 = - (points_3[int(3)] - points_3[int(2)]);
            float3  _S469 = s_primal_ctx_cross_0(dx_1, _S468);
            bool _S470 = (s_primal_ctx_dot_0(_S469, _S469)) != 0.0f;
            if(_S470)
            {
                float _S471 = length_0(_S469);
                float3  _S472 = make_float3 (_S471);
                _S452 = make_float3 (_S471 * _S471);
                _S453 = _S472;
            }
            else
            {
                _S452 = _S443;
                _S453 = _S443;
            }
            float3  _S473 = _S453;
            _S446 = _S470;
            _S453 = _S469;
            _S454 = _S473;
            _S455 = dx_1;
            _S456 = _S468;
        }
        else
        {
            _S446 = false;
            _S452 = _S443;
            _S453 = _S443;
            _S454 = _S443;
            _S455 = _S443;
            _S456 = _S443;
        }
        bool _S474 = _runFlag_0;
        bool _S475 = _runFlag_1;
        bool _S476 = _runFlag_2;
        float3  _S477 = _S448;
        float3  _S478 = _S449;
        bool _S479 = _S445;
        float3  _S480 = _S450;
        float3  _S481 = _S451;
        _runFlag_0 = _runFlag_3;
        _runFlag_1 = _S446;
        _S448 = _S452;
        _S449 = _S453;
        _S450 = _S454;
        _S451 = _S455;
        _S452 = _S456;
        _runFlag_2 = _S457;
        _S445 = _S474;
        _runFlag_3 = _S475;
        _S446 = _S476;
        _S453 = _S477;
        _S454 = _S478;
        _S447 = _S479;
        _S455 = _S480;
        _S456 = _S481;
    }
    else
    {
        _runFlag_0 = false;
        _runFlag_1 = false;
        _S448 = _S443;
        _S449 = _S443;
        _S450 = _S443;
        _S451 = _S443;
        _S452 = _S443;
        _runFlag_2 = false;
        _S445 = false;
        _runFlag_3 = false;
        _S446 = false;
        _S453 = _S443;
        _S454 = _S443;
        _S447 = false;
        _S455 = _S443;
        _S456 = _S443;
    }
    float4  _S482 = make_float4 (0.0f);
    float4  _S483;
    if(_S444)
    {
        if(_runFlag_0)
        {
            if(_runFlag_1)
            {
                float3  _S484 = v_normal_1 / _S448;
                float3  _S485 = _S449 * - _S484;
                float3  _S486 = _S450 * _S484;
                float _S487 = _S485.x + _S485.y + _S485.z;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S488;
                (&_S488)->primal_0 = _S449;
                (&_S488)->differential_0 = _S443;
                s_bwd_length_impl_0(&_S488, _S487);
                _S448 = _S486 + _S488.differential_0;
            }
            else
            {
                _S448 = v_normal_1;
            }
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S489;
            (&_S489)->primal_0 = _S449;
            (&_S489)->differential_0 = _S443;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S490;
            (&_S490)->primal_0 = _S449;
            (&_S490)->differential_0 = _S443;
            s_bwd_prop_dot_0(&_S489, &_S490, 0.0f);
            float3  _S491 = _S490.differential_0 + _S489.differential_0 + _S448;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S492;
            (&_S492)->primal_0 = _S451;
            (&_S492)->differential_0 = _S443;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S493;
            (&_S493)->primal_0 = _S452;
            (&_S493)->differential_0 = _S443;
            s_bwd_prop_cross_0(&_S492, &_S493, _S491);
            float3  s_diff_dy_T_1 = - _S493.differential_0;
            float3  _S494 = - s_diff_dy_T_1;
            float3  _S495 = - _S492.differential_0;
            FixedArray<float3 , 4>  _S496;
            _S496[int(0)] = _S443;
            _S496[int(1)] = _S443;
            _S496[int(2)] = _S443;
            _S496[int(3)] = _S443;
            _S496[int(2)] = _S494;
            _S496[int(3)] = s_diff_dy_T_1;
            _S496[int(0)] = _S495;
            _S496[int(1)] = _S492.differential_0;
            points_3[int(0)] = _S496[int(0)];
            points_3[int(1)] = _S496[int(1)];
            points_3[int(2)] = _S496[int(2)];
            points_3[int(3)] = _S496[int(3)];
        }
        else
        {
            points_3[int(0)] = _S443;
            points_3[int(1)] = _S443;
            points_3[int(2)] = _S443;
            points_3[int(3)] = _S443;
        }
        if(_runFlag_2)
        {
            if(_S445)
            {
                if(_runFlag_3)
                {
                    FixedArray<float3 , 4>  _S497 = points_3;
                    FixedArray<float3 , 4>  _S498 = points_3;
                    FixedArray<float3 , 4>  _S499 = points_3;
                    FixedArray<float3 , 4>  _S500 = points_3;
                    if(_S446)
                    {
                        float3  _S501 = _S453 * _S500[int(3)];
                        float _S502 = _S501.x + _S501.y + _S501.z;
                        float4  _S503 = _S482;
                        *&((&_S503)->w) = _S502;
                        points_3[int(0)] = _S497[int(0)];
                        points_3[int(1)] = _S498[int(1)];
                        points_3[int(2)] = _S499[int(2)];
                        points_3[int(3)] = _S443;
                        _S483 = _S503;
                    }
                    else
                    {
                        points_3[int(0)] = _S497[int(0)];
                        points_3[int(1)] = _S498[int(1)];
                        points_3[int(2)] = _S499[int(2)];
                        points_3[int(3)] = _S500[int(3)];
                        _S483 = _S482;
                    }
                    float3  _S504 = _S454 * points_3[int(2)];
                    float _S505 = _S504.x + _S504.y + _S504.z;
                    FixedArray<float3 , 4>  _S506 = points_3;
                    FixedArray<float3 , 4>  _S507 = points_3;
                    float4  _S508 = _S482;
                    *&((&_S508)->z) = _S505;
                    float4  _S509 = _S483 + _S508;
                    points_3[int(0)] = points_3[int(0)];
                    points_3[int(1)] = _S506[int(1)];
                    points_3[int(2)] = _S443;
                    points_3[int(3)] = _S507[int(3)];
                    _S483 = _S509;
                }
                else
                {
                    FixedArray<float3 , 4>  _S510 = points_3;
                    FixedArray<float3 , 4>  _S511 = points_3;
                    FixedArray<float3 , 4>  _S512 = points_3;
                    points_3[int(0)] = points_3[int(0)];
                    points_3[int(1)] = _S510[int(1)];
                    points_3[int(2)] = _S511[int(2)];
                    points_3[int(3)] = _S512[int(3)];
                    _S483 = _S482;
                }
            }
            else
            {
                FixedArray<float3 , 4>  _S513 = points_3;
                FixedArray<float3 , 4>  _S514 = points_3;
                FixedArray<float3 , 4>  _S515 = points_3;
                points_3[int(0)] = points_3[int(0)];
                points_3[int(1)] = _S513[int(1)];
                points_3[int(2)] = _S514[int(2)];
                points_3[int(3)] = _S515[int(3)];
                _S483 = _S482;
            }
            if(_S447)
            {
                FixedArray<float3 , 4>  _S516 = points_3;
                float3  _S517 = _S455 * points_3[int(1)];
                float _S518 = _S517.x + _S517.y + _S517.z;
                float4  _S519 = _S482;
                *&((&_S519)->y) = _S518;
                float4  _S520 = _S483 + _S519;
                points_3[int(0)] = _S443;
                points_3[int(1)] = _S443;
                points_3[int(2)] = _S443;
                points_3[int(3)] = _S443;
                _S448 = _S516[int(0)];
                _S483 = _S520;
            }
            else
            {
                FixedArray<float3 , 4>  _S521 = points_3;
                FixedArray<float3 , 4>  _S522 = points_3;
                FixedArray<float3 , 4>  _S523 = points_3;
                points_3[int(0)] = points_3[int(0)];
                points_3[int(1)] = _S521[int(1)];
                points_3[int(2)] = _S522[int(2)];
                points_3[int(3)] = _S523[int(3)];
                _S448 = _S443;
            }
            float3  _S524 = _S456 * (points_3[int(0)] + _S448);
            float _S525 = _S524.x + _S524.y + _S524.z;
            float4  _S526 = _S482;
            *&((&_S526)->x) = _S525;
            _S483 = _S483 + _S526;
        }
        else
        {
            _S483 = _S482;
        }
    }
    else
    {
        _S483 = _S482;
    }
    *v_depths_0 = _S483;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor_none(float2  pix_center_2, float4  intrins_5, FixedArray<float, 1>  dist_coeffs_9, int camera_model_7)
{
    float _S527;
    for(;;)
    {
        float2  uv_19 = (pix_center_2 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
        FixedArray<float, 1>  _S528 = dist_coeffs_9;
        float2  uv_u_6;
        bool _S529 = undistort_point_0(uv_19, &_S528, int(12), &uv_u_6);
        if(!_S529)
        {
            _S527 = 0.0f;
            break;
        }
        float3  raydir_2 = unproject_raydir_0(uv_u_6, camera_model_7, false);
        _S527 = float((F32_sign((raydir_2.z)))) / length_0(raydir_2);
        break;
    }
    return _S527;
}

inline __device__ float depth_normal_loss_none(float2  pix_center_3, float4  intrins_6, FixedArray<float, 1>  dist_coeffs_10, int camera_model_8, bool is_ray_depth_7, float4  depths_2, float3  gt_normal_0)
{
    float _S530;
    for(;;)
    {
        float3  _S531;
        float3  * _S532;
        float3  * _S533;
        float3  * _S534;
        float3  * _S535;
        int _S536;
        FixedArray<float3 , 5>  points_4;
        for(;;)
        {
            float2  _S537 = float2 {intrins_6.z, intrins_6.w};
            float2  _S538 = float2 {intrins_6.x, intrins_6.y};
            float2  uv_20 = (pix_center_3 + make_float2 (-1.0f, -0.0f) - _S537) / _S538;
            FixedArray<float, 1>  _S539 = dist_coeffs_10;
            float2  uv_u_7;
            bool _S540 = undistort_point_0(uv_20, &_S539, int(12), &uv_u_7);
            float3  _S541 = make_float3 (0.0f);
            if(!_S540)
            {
                _S536 = int(0);
                _S535 = nullptr;
                _S534 = nullptr;
                _S533 = nullptr;
                _S532 = nullptr;
                _S531 = _S541;
                break;
            }
            float3  raydir_3 = unproject_raydir_0(uv_u_7, camera_model_8, is_ray_depth_7);
            points_4[int(0)] = make_float3 (depths_2.x) * raydir_3;
            float2  uv_21 = (pix_center_3 + make_float2 (1.0f, -0.0f) - _S537) / _S538;
            FixedArray<float, 1>  _S542 = dist_coeffs_10;
            float2  uv_u_8;
            bool _S543 = undistort_point_0(uv_21, &_S542, int(12), &uv_u_8);
            if(!_S543)
            {
                _S536 = int(0);
                _S535 = nullptr;
                _S534 = &points_4[int(0)];
                _S533 = nullptr;
                _S532 = nullptr;
                _S531 = _S541;
                break;
            }
            float3  raydir_4 = unproject_raydir_0(uv_u_8, camera_model_8, is_ray_depth_7);
            points_4[int(1)] = make_float3 (depths_2.y) * raydir_4;
            float2  uv_22 = (pix_center_3 + make_float2 (0.0f, -1.0f) - _S537) / _S538;
            FixedArray<float, 1>  _S544 = dist_coeffs_10;
            float2  uv_u_9;
            bool _S545 = undistort_point_0(uv_22, &_S544, int(12), &uv_u_9);
            if(!_S545)
            {
                _S536 = int(0);
                _S535 = &points_4[int(1)];
                _S534 = &points_4[int(0)];
                _S533 = nullptr;
                _S532 = nullptr;
                _S531 = _S541;
                break;
            }
            float3  raydir_5 = unproject_raydir_0(uv_u_9, camera_model_8, is_ray_depth_7);
            points_4[int(2)] = make_float3 (depths_2.z) * raydir_5;
            float2  uv_23 = (pix_center_3 + make_float2 (0.0f, 1.0f) - _S537) / _S538;
            FixedArray<float, 1>  _S546 = dist_coeffs_10;
            float2  uv_u_10;
            bool _S547 = undistort_point_0(uv_23, &_S546, int(12), &uv_u_10);
            if(!_S547)
            {
                _S536 = int(0);
                _S535 = &points_4[int(1)];
                _S534 = &points_4[int(0)];
                _S533 = nullptr;
                _S532 = &points_4[int(2)];
                _S531 = _S541;
                break;
            }
            float3  raydir_6 = unproject_raydir_0(uv_u_10, camera_model_8, is_ray_depth_7);
            points_4[int(3)] = make_float3 (depths_2.w) * raydir_6;
            float2  uv_24 = (pix_center_3 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S537) / _S538;
            FixedArray<float, 1>  _S548 = dist_coeffs_10;
            float2  uv_u_11;
            bool _S549 = undistort_point_0(uv_24, &_S548, int(12), &uv_u_11);
            if(!_S549)
            {
                _S536 = int(0);
                _S535 = &points_4[int(1)];
                _S534 = &points_4[int(0)];
                _S533 = &points_4[int(3)];
                _S532 = &points_4[int(2)];
                _S531 = _S541;
                break;
            }
            float3  raydir_7 = unproject_raydir_0(uv_u_11, camera_model_8, is_ray_depth_7);
            _S536 = int(1);
            _S535 = &points_4[int(1)];
            _S534 = &points_4[int(0)];
            _S533 = &points_4[int(3)];
            _S532 = &points_4[int(2)];
            _S531 = raydir_7;
            break;
        }
        if(_S536 != int(1))
        {
            _S530 = 0.0f;
            break;
        }
        float3  normal_4 = cross_0(*_S535 - *_S534, - (*_S533 - *_S532));
        float3  normal_5;
        if((dot_0(normal_4, normal_4)) != 0.0f)
        {
            normal_5 = normalize_0(normal_4);
        }
        else
        {
            normal_5 = normal_4;
        }
        float3  _S550;
        if((dot_0(gt_normal_0, gt_normal_0)) != 0.0f)
        {
            _S550 = normalize_0(gt_normal_0);
        }
        else
        {
            _S550 = gt_normal_0;
        }
        _S530 = (1.0f - dot_0(normal_5, _S550) + 0.00100000004749745f) / ((F32_max((dot_0(normal_5, - normalize_0(_S531))), (0.0f))) + 0.00100000004749745f);
        break;
    }
    return _S530;
}

struct s_bwd_prop_depth_normal_loss_Intermediates_0
{
    float2  _S551;
    bool _S552;
    float2  _S553;
    bool _S554;
    float2  _S555;
    bool _S556;
    float2  _S557;
    bool _S558;
    float2  _S559;
    bool _S560;
};

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  _s_dOut_6)
{
    float _S561 = length_0((*dpx_8).primal_0);
    float3  _S562 = (*dpx_8).primal_0 * _s_dOut_6;
    float3  _S563 = make_float3 (1.0f / _S561) * _s_dOut_6;
    float _S564 = - ((_S562.x + _S562.y + _S562.z) / (_S561 * _S561));
    float3  _S565 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S566;
    (&_S566)->primal_0 = (*dpx_8).primal_0;
    (&_S566)->differential_0 = _S565;
    s_bwd_length_impl_0(&_S566, _S564);
    float3  _S567 = _S563 + _S566.differential_0;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S567;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S568, float3  _S569)
{
    s_bwd_prop_normalize_impl_0(_S568, _S569);
    return;
}

inline __device__ void depth_normal_loss_vjp_none(float2  pix_center_4, float4  intrins_7, FixedArray<float, 1>  dist_coeffs_11, int camera_model_9, bool is_ray_depth_8, float4  depths_3, float3  gt_normal_1, float v_loss_0, float4  * v_depths_1, float3  * v_gt_normal_0)
{
    float2  _S570 = make_float2 (0.0f);
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S571;
    (&_S571)->_S551 = _S570;
    (&_S571)->_S552 = false;
    (&_S571)->_S553 = _S570;
    (&_S571)->_S554 = false;
    (&_S571)->_S555 = _S570;
    (&_S571)->_S556 = false;
    (&_S571)->_S557 = _S570;
    (&_S571)->_S558 = false;
    (&_S571)->_S559 = _S570;
    (&_S571)->_S560 = false;
    (&_S571)->_S553 = _S570;
    (&_S571)->_S554 = false;
    (&_S571)->_S555 = _S570;
    (&_S571)->_S556 = false;
    (&_S571)->_S557 = _S570;
    (&_S571)->_S558 = false;
    (&_S571)->_S559 = _S570;
    (&_S571)->_S560 = false;
    float2  _S572 = float2 {intrins_7.z, intrins_7.w};
    float2  _S573 = float2 {intrins_7.x, intrins_7.y};
    float2  uv_25 = (pix_center_4 + make_float2 (-1.0f, -0.0f) - _S572) / _S573;
    float2  _S574 = _S570;
    FixedArray<float, 1>  _S575 = dist_coeffs_11;
    bool _S576 = undistort_point_0(uv_25, &_S575, int(12), &_S574);
    (&_S571)->_S551 = _S574;
    (&_S571)->_S552 = _S576;
    bool _S577 = !!_S576;
    bool _runFlag_4;
    if(_S577)
    {
        float2  uv_26 = (pix_center_4 + make_float2 (1.0f, -0.0f) - _S572) / _S573;
        float2  _S578 = _S570;
        FixedArray<float, 1>  _S579 = dist_coeffs_11;
        bool _S580 = undistort_point_0(uv_26, &_S579, int(12), &_S578);
        (&_S571)->_S553 = _S578;
        (&_S571)->_S554 = _S580;
        if(!_S580)
        {
            _runFlag_4 = false;
        }
        else
        {
            _runFlag_4 = _S577;
        }
        if(_runFlag_4)
        {
            float2  uv_27 = (pix_center_4 + make_float2 (0.0f, -1.0f) - _S572) / _S573;
            float2  _S581 = _S570;
            FixedArray<float, 1>  _S582 = dist_coeffs_11;
            bool _S583 = undistort_point_0(uv_27, &_S582, int(12), &_S581);
            (&_S571)->_S555 = _S581;
            (&_S571)->_S556 = _S583;
            if(!_S583)
            {
                _runFlag_4 = false;
            }
            if(_runFlag_4)
            {
                float2  uv_28 = (pix_center_4 + make_float2 (0.0f, 1.0f) - _S572) / _S573;
                float2  _S584 = _S570;
                FixedArray<float, 1>  _S585 = dist_coeffs_11;
                bool _S586 = undistort_point_0(uv_28, &_S585, int(12), &_S584);
                (&_S571)->_S557 = _S584;
                (&_S571)->_S558 = _S586;
                if(!_S586)
                {
                    _runFlag_4 = false;
                }
                if(_runFlag_4)
                {
                    float2  uv_29 = (pix_center_4 - _S572) / _S573;
                    float2  _S587 = _S570;
                    FixedArray<float, 1>  _S588 = dist_coeffs_11;
                    bool _S589 = undistort_point_0(uv_29, &_S588, int(12), &_S587);
                    (&_S571)->_S559 = _S587;
                    (&_S571)->_S560 = _S589;
                }
            }
        }
    }
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S590 = _S571;
    float3  _S591 = make_float3 (0.0f);
    bool _S592 = !!_S571._S552;
    bool _runFlag_5;
    bool _runFlag_6;
    bool _runFlag_7;
    int _S593;
    float3  raydir_8;
    float3  _S594;
    float3  _S595;
    float3  _S596;
    float3  _S597;
    FixedArray<float3 , 5>  points_5;
    if(_S592)
    {
        float3  _S598 = s_primal_ctx_unproject_raydir_0(_S590._S551, camera_model_9, is_ray_depth_8);
        float3  _S599 = make_float3 (depths_3.x) * _S598;
        if(!_S590._S554)
        {
            _runFlag_4 = false;
        }
        else
        {
            _runFlag_4 = _S592;
        }
        if(_runFlag_4)
        {
            float3  _S600 = s_primal_ctx_unproject_raydir_0(_S590._S553, camera_model_9, is_ray_depth_8);
            float3  _S601 = make_float3 (depths_3.y) * _S600;
            if(!_S590._S556)
            {
                _runFlag_5 = false;
            }
            else
            {
                _runFlag_5 = _runFlag_4;
            }
            if(_runFlag_5)
            {
                float3  _S602 = s_primal_ctx_unproject_raydir_0(_S590._S555, camera_model_9, is_ray_depth_8);
                float3  _S603 = make_float3 (depths_3.z) * _S602;
                if(!_S590._S558)
                {
                    _runFlag_6 = false;
                }
                else
                {
                    _runFlag_6 = _runFlag_5;
                }
                if(_runFlag_6)
                {
                    float3  _S604 = s_primal_ctx_unproject_raydir_0(_S590._S557, camera_model_9, is_ray_depth_8);
                    float3  _S605 = make_float3 (depths_3.w) * _S604;
                    if(!_S590._S560)
                    {
                        _runFlag_7 = false;
                    }
                    else
                    {
                        _runFlag_7 = _runFlag_6;
                    }
                    if(_runFlag_7)
                    {
                        float3  _S606 = s_primal_ctx_unproject_raydir_0(_S590._S559, camera_model_9, is_ray_depth_8);
                        _S593 = int(1);
                        raydir_8 = _S606;
                    }
                    else
                    {
                        _S593 = int(0);
                        raydir_8 = _S604;
                    }
                    points_5[int(0)] = _S599;
                    points_5[int(1)] = _S601;
                    points_5[int(2)] = _S603;
                    points_5[int(3)] = _S605;
                    points_5[int(4)] = _S591;
                    _S594 = _S604;
                }
                else
                {
                    _S593 = int(0);
                    raydir_8 = _S602;
                    points_5[int(0)] = _S599;
                    points_5[int(1)] = _S601;
                    points_5[int(2)] = _S603;
                    points_5[int(3)] = _S591;
                    points_5[int(4)] = _S591;
                    _S594 = _S591;
                }
                _S595 = _S602;
            }
            else
            {
                _S593 = int(0);
                raydir_8 = _S600;
                points_5[int(0)] = _S599;
                points_5[int(1)] = _S601;
                points_5[int(2)] = _S591;
                points_5[int(3)] = _S591;
                points_5[int(4)] = _S591;
                _runFlag_6 = false;
                _S594 = _S591;
                _S595 = _S591;
            }
            _S596 = _S600;
        }
        else
        {
            _S593 = int(0);
            raydir_8 = _S598;
            points_5[int(0)] = _S599;
            points_5[int(1)] = _S591;
            points_5[int(2)] = _S591;
            points_5[int(3)] = _S591;
            points_5[int(4)] = _S591;
            _runFlag_5 = false;
            _runFlag_6 = false;
            _S594 = _S591;
            _S595 = _S591;
            _S596 = _S591;
        }
        _S597 = _S598;
    }
    else
    {
        _S593 = int(0);
        points_5[int(0)] = _S591;
        points_5[int(1)] = _S591;
        points_5[int(2)] = _S591;
        points_5[int(3)] = _S591;
        points_5[int(4)] = _S591;
        _runFlag_4 = false;
        _runFlag_5 = false;
        _runFlag_6 = false;
        _S594 = _S591;
        _S595 = _S591;
        _S596 = _S591;
        _S597 = _S591;
    }
    bool _S607 = !(_S593 != int(1));
    bool _S608;
    float3  normal_6;
    float3  _S609;
    float3  _S610;
    float3  _S611;
    float3  _S612;
    float _S613;
    float _S614;
    float _S615;
    float _S616;
    if(_S607)
    {
        float3  dx_2 = points_5[int(1)] - points_5[int(0)];
        float3  _S617 = - (points_5[int(3)] - points_5[int(2)]);
        float3  _S618 = s_primal_ctx_cross_0(dx_2, _S617);
        bool _S619 = (s_primal_ctx_dot_0(_S618, _S618)) != 0.0f;
        if(_S619)
        {
            normal_6 = normalize_0(_S618);
        }
        else
        {
            normal_6 = _S618;
        }
        bool _S620 = (s_primal_ctx_dot_0(gt_normal_1, gt_normal_1)) != 0.0f;
        if(_S620)
        {
            _S609 = normalize_0(gt_normal_1);
        }
        else
        {
            _S609 = gt_normal_1;
        }
        float3  _S621 = - normalize_0(raydir_8);
        float _S622 = s_primal_ctx_dot_0(normal_6, _S621);
        float _S623 = 1.0f - s_primal_ctx_dot_0(normal_6, _S609) + 0.00100000004749745f;
        float _S624 = (F32_max((_S622), (0.0f))) + 0.00100000004749745f;
        _S613 = _S624 * _S624;
        _S614 = _S623;
        _S615 = _S624;
        _S616 = _S622;
        raydir_8 = normal_6;
        normal_6 = _S621;
        _runFlag_7 = _S620;
        _S608 = _S619;
        _S610 = _S618;
        _S611 = dx_2;
        _S612 = _S617;
    }
    else
    {
        _S613 = 0.0f;
        _S614 = 0.0f;
        _S615 = 0.0f;
        _S616 = 0.0f;
        raydir_8 = _S591;
        normal_6 = _S591;
        _S609 = _S591;
        _runFlag_7 = false;
        _S608 = false;
        _S610 = _S591;
        _S611 = _S591;
        _S612 = _S591;
    }
    float4  _S625 = make_float4 (0.0f);
    if(_S607)
    {
        float _S626 = v_loss_0 / _S613;
        float _S627 = _S614 * - _S626;
        float s_diff_num_T_0 = _S615 * _S626;
        DiffPair_float_0 _S628;
        (&_S628)->primal_0 = _S616;
        (&_S628)->differential_0 = 0.0f;
        DiffPair_float_0 _S629;
        (&_S629)->primal_0 = 0.0f;
        (&_S629)->differential_0 = 0.0f;
        _d_max_0(&_S628, &_S629, _S627);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S630;
        (&_S630)->primal_0 = raydir_8;
        (&_S630)->differential_0 = _S591;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S631;
        (&_S631)->primal_0 = normal_6;
        (&_S631)->differential_0 = _S591;
        s_bwd_prop_dot_0(&_S630, &_S631, _S628.differential_0);
        float _S632 = - s_diff_num_T_0;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S633;
        (&_S633)->primal_0 = raydir_8;
        (&_S633)->differential_0 = _S591;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S634;
        (&_S634)->primal_0 = _S609;
        (&_S634)->differential_0 = _S591;
        s_bwd_prop_dot_0(&_S633, &_S634, _S632);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S635 = _S634;
        float3  _S636 = _S630.differential_0 + _S633.differential_0;
        if(_runFlag_7)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S637;
            (&_S637)->primal_0 = gt_normal_1;
            (&_S637)->differential_0 = _S591;
            s_bwd_normalize_impl_0(&_S637, _S635.differential_0);
            raydir_8 = _S637.differential_0;
        }
        else
        {
            raydir_8 = _S635.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
        (&_S638)->primal_0 = gt_normal_1;
        (&_S638)->differential_0 = _S591;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S639;
        (&_S639)->primal_0 = gt_normal_1;
        (&_S639)->differential_0 = _S591;
        s_bwd_prop_dot_0(&_S638, &_S639, 0.0f);
        float3  _S640 = _S639.differential_0 + _S638.differential_0 + raydir_8;
        if(_S608)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S641;
            (&_S641)->primal_0 = _S610;
            (&_S641)->differential_0 = _S591;
            s_bwd_normalize_impl_0(&_S641, _S636);
            raydir_8 = _S641.differential_0;
        }
        else
        {
            raydir_8 = _S636;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S642;
        (&_S642)->primal_0 = _S610;
        (&_S642)->differential_0 = _S591;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S643;
        (&_S643)->primal_0 = _S610;
        (&_S643)->differential_0 = _S591;
        s_bwd_prop_dot_0(&_S642, &_S643, 0.0f);
        float3  _S644 = _S643.differential_0 + _S642.differential_0 + raydir_8;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S645;
        (&_S645)->primal_0 = _S611;
        (&_S645)->differential_0 = _S591;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S646;
        (&_S646)->primal_0 = _S612;
        (&_S646)->differential_0 = _S591;
        s_bwd_prop_cross_0(&_S645, &_S646, _S644);
        float3  s_diff_dy_T_2 = - _S646.differential_0;
        float3  _S647 = - s_diff_dy_T_2;
        float3  _S648 = - _S645.differential_0;
        FixedArray<float3 , 5>  _S649;
        _S649[int(0)] = _S591;
        _S649[int(1)] = _S591;
        _S649[int(2)] = _S591;
        _S649[int(3)] = _S591;
        _S649[int(4)] = _S591;
        _S649[int(2)] = _S647;
        _S649[int(3)] = s_diff_dy_T_2;
        _S649[int(0)] = _S648;
        _S649[int(1)] = _S645.differential_0;
        points_5[int(0)] = _S649[int(0)];
        points_5[int(1)] = _S649[int(1)];
        points_5[int(2)] = _S649[int(2)];
        points_5[int(3)] = _S649[int(3)];
        points_5[int(4)] = _S649[int(4)];
        raydir_8 = _S640;
    }
    else
    {
        points_5[int(0)] = _S591;
        points_5[int(1)] = _S591;
        points_5[int(2)] = _S591;
        points_5[int(3)] = _S591;
        points_5[int(4)] = _S591;
        raydir_8 = _S591;
    }
    float4  _S650;
    if(_S592)
    {
        if(_runFlag_4)
        {
            if(_runFlag_5)
            {
                if(_runFlag_6)
                {
                    FixedArray<float3 , 5>  _S651 = points_5;
                    FixedArray<float3 , 5>  _S652 = points_5;
                    FixedArray<float3 , 5>  _S653 = points_5;
                    float3  _S654 = _S594 * points_5[int(3)];
                    float _S655 = _S654.x + _S654.y + _S654.z;
                    float4  _S656 = _S625;
                    *&((&_S656)->w) = _S655;
                    points_5[int(0)] = _S591;
                    points_5[int(1)] = _S591;
                    points_5[int(2)] = _S591;
                    points_5[int(3)] = _S591;
                    points_5[int(4)] = _S591;
                    _S594 = _S653[int(2)];
                    normal_6 = _S651[int(0)];
                    _S609 = _S652[int(1)];
                    _S650 = _S656;
                }
                else
                {
                    FixedArray<float3 , 5>  _S657 = points_5;
                    FixedArray<float3 , 5>  _S658 = points_5;
                    FixedArray<float3 , 5>  _S659 = points_5;
                    FixedArray<float3 , 5>  _S660 = points_5;
                    points_5[int(0)] = points_5[int(0)];
                    points_5[int(1)] = _S657[int(1)];
                    points_5[int(2)] = _S658[int(2)];
                    points_5[int(3)] = _S659[int(3)];
                    points_5[int(4)] = _S660[int(4)];
                    _S594 = _S591;
                    normal_6 = _S591;
                    _S609 = _S591;
                    _S650 = _S625;
                }
                float3  _S661 = _S595 * (points_5[int(2)] + _S594);
                float _S662 = _S661.x + _S661.y + _S661.z;
                float3  _S663 = points_5[int(0)] + normal_6;
                float3  _S664 = points_5[int(1)] + _S609;
                float4  _S665 = _S625;
                *&((&_S665)->z) = _S662;
                float4  _S666 = _S650 + _S665;
                points_5[int(0)] = _S591;
                points_5[int(1)] = _S591;
                points_5[int(2)] = _S591;
                points_5[int(3)] = _S591;
                points_5[int(4)] = _S591;
                _S594 = _S664;
                _S595 = _S663;
                _S650 = _S666;
            }
            else
            {
                FixedArray<float3 , 5>  _S667 = points_5;
                FixedArray<float3 , 5>  _S668 = points_5;
                FixedArray<float3 , 5>  _S669 = points_5;
                FixedArray<float3 , 5>  _S670 = points_5;
                points_5[int(0)] = points_5[int(0)];
                points_5[int(1)] = _S667[int(1)];
                points_5[int(2)] = _S668[int(2)];
                points_5[int(3)] = _S669[int(3)];
                points_5[int(4)] = _S670[int(4)];
                _S594 = _S591;
                _S595 = _S591;
                _S650 = _S625;
            }
            float3  _S671 = _S596 * (points_5[int(1)] + _S594);
            float _S672 = _S671.x + _S671.y + _S671.z;
            float3  _S673 = points_5[int(0)] + _S595;
            float4  _S674 = _S625;
            *&((&_S674)->y) = _S672;
            float4  _S675 = _S650 + _S674;
            points_5[int(0)] = _S591;
            points_5[int(1)] = _S591;
            points_5[int(2)] = _S591;
            points_5[int(3)] = _S591;
            points_5[int(4)] = _S591;
            _S594 = _S673;
            _S650 = _S675;
        }
        else
        {
            FixedArray<float3 , 5>  _S676 = points_5;
            FixedArray<float3 , 5>  _S677 = points_5;
            FixedArray<float3 , 5>  _S678 = points_5;
            FixedArray<float3 , 5>  _S679 = points_5;
            points_5[int(0)] = points_5[int(0)];
            points_5[int(1)] = _S676[int(1)];
            points_5[int(2)] = _S677[int(2)];
            points_5[int(3)] = _S678[int(3)];
            points_5[int(4)] = _S679[int(4)];
            _S594 = _S591;
            _S650 = _S625;
        }
        float3  _S680 = _S597 * (points_5[int(0)] + _S594);
        float _S681 = _S680.x + _S680.y + _S680.z;
        float4  _S682 = _S625;
        *&((&_S682)->x) = _S681;
        _S650 = _S650 + _S682;
    }
    else
    {
        _S650 = _S625;
    }
    *v_depths_1 = _S650;
    *v_gt_normal_0 = raydir_8;
    return;
}

inline __device__ float3  generate_ray_d2n_opencv(float2  pix_pos_3, float4  intrins_8, FixedArray<float, 4>  dist_coeffs_12, int camera_model_10, bool is_ray_depth_9)
{
    float3  _S683;
    for(;;)
    {
        float2  uv_30 = (pix_pos_3 - float2 {intrins_8.z, intrins_8.w}) / float2 {intrins_8.x, intrins_8.y};
        FixedArray<float, 4>  _S684 = dist_coeffs_12;
        float2  uv_u_12;
        bool _S685 = undistort_point_1(uv_30, &_S684, int(12), &uv_u_12);
        if(!_S685)
        {
            int3  _S686 = make_int3 (int(0));
            float3  _S687 = make_float3 ((float)_S686.x, (float)_S686.y, (float)_S686.z);
            _S683 = _S687;
            break;
        }
        _S683 = unproject_raydir_0(uv_u_12, camera_model_10, is_ray_depth_9);
        break;
    }
    return _S683;
}

inline __device__ float3  depth_to_point_opencv(float2  pix_pos_4, float4  intrins_9, FixedArray<float, 4>  dist_coeffs_13, int camera_model_11, bool is_ray_depth_10, float depth_4)
{
    float3  _S688;
    for(;;)
    {
        float2  uv_31 = (pix_pos_4 - float2 {intrins_9.z, intrins_9.w}) / float2 {intrins_9.x, intrins_9.y};
        FixedArray<float, 4>  _S689 = dist_coeffs_13;
        float2  uv_u_13;
        bool _S690 = undistort_point_1(uv_31, &_S689, int(12), &uv_u_13);
        if(!_S690)
        {
            _S688 = make_float3 (0.0f);
            break;
        }
        _S688 = make_float3 (depth_4) * unproject_raydir_0(uv_u_13, camera_model_11, is_ray_depth_10);
        break;
    }
    return _S688;
}

struct s_bwd_prop_depth_to_point_Intermediates_1
{
    float2  _S691;
    bool _S692;
};

inline __device__ float depth_to_point_vjp_opencv(float2  pix_pos_5, float4  intrins_10, FixedArray<float, 4>  dist_coeffs_14, int camera_model_12, bool is_ray_depth_11, float depth_5, float3  v_point_1)
{
    float2  _S693 = make_float2 (0.0f);
    s_bwd_prop_depth_to_point_Intermediates_1 _S694;
    (&_S694)->_S691 = _S693;
    (&_S694)->_S692 = false;
    float2  uv_32 = (pix_pos_5 - float2 {intrins_10.z, intrins_10.w}) / float2 {intrins_10.x, intrins_10.y};
    float2  _S695 = _S693;
    FixedArray<float, 4>  _S696 = dist_coeffs_14;
    bool _S697 = undistort_point_1(uv_32, &_S696, int(12), &_S695);
    (&_S694)->_S691 = _S695;
    (&_S694)->_S692 = _S697;
    s_bwd_prop_depth_to_point_Intermediates_1 _S698 = _S694;
    float3  _S699 = make_float3 (0.0f);
    bool _S700 = !!_S694._S692;
    float3  _S701;
    if(_S700)
    {
        _S701 = s_primal_ctx_unproject_raydir_0(_S698._S691, camera_model_12, is_ray_depth_11);
    }
    else
    {
        _S701 = _S699;
    }
    if(_S700)
    {
        _S701 = _S701 * v_point_1;
    }
    else
    {
        _S701 = _S699;
    }
    return _S701.x + _S701.y + _S701.z;
}

inline __device__ float3  depth_to_normal_opencv(float2  pix_center_5, float4  intrins_11, FixedArray<float, 4>  dist_coeffs_15, int camera_model_13, bool is_ray_depth_12, float4  depths_4)
{
    float3  normal_7;
    for(;;)
    {
        bool _S702;
        if((depths_4.x) == 0.0f)
        {
            _S702 = true;
        }
        else
        {
            _S702 = (depths_4.y) == 0.0f;
        }
        if(_S702)
        {
            _S702 = true;
        }
        else
        {
            _S702 = (depths_4.z) == 0.0f;
        }
        if(_S702)
        {
            _S702 = true;
        }
        else
        {
            _S702 = (depths_4.w) == 0.0f;
        }
        if(_S702)
        {
            normal_7 = make_float3 (0.0f);
            break;
        }
        float3  * _S703;
        float3  * _S704;
        float3  * _S705;
        float3  * _S706;
        int _S707;
        FixedArray<float3 , 4>  points_6;
        for(;;)
        {
            float2  _S708 = float2 {intrins_11.z, intrins_11.w};
            float2  _S709 = float2 {intrins_11.x, intrins_11.y};
            float2  uv_33 = (pix_center_5 + make_float2 (-1.0f, -0.0f) - _S708) / _S709;
            FixedArray<float, 4>  _S710 = dist_coeffs_15;
            float2  uv_u_14;
            bool _S711 = undistort_point_1(uv_33, &_S710, int(12), &uv_u_14);
            if(!_S711)
            {
                float3  _S712 = make_float3 (0.0f);
                _S707 = int(0);
                _S706 = nullptr;
                _S705 = nullptr;
                _S704 = nullptr;
                _S703 = nullptr;
                normal_7 = _S712;
                break;
            }
            points_6[int(0)] = make_float3 (depths_4.x) * unproject_raydir_0(uv_u_14, camera_model_13, is_ray_depth_12);
            for(;;)
            {
                float2  uv_34 = (pix_center_5 + make_float2 (1.0f, -0.0f) - _S708) / _S709;
                FixedArray<float, 4>  _S713 = dist_coeffs_15;
                float2  uv_u_15;
                bool _S714 = undistort_point_1(uv_34, &_S713, int(12), &uv_u_15);
                if(!_S714)
                {
                    float3  _S715 = make_float3 (0.0f);
                    _S707 = int(0);
                    _S706 = nullptr;
                    normal_7 = _S715;
                    break;
                }
                points_6[int(1)] = make_float3 (depths_4.y) * unproject_raydir_0(uv_u_15, camera_model_13, is_ray_depth_12);
                _S707 = int(2);
                _S706 = &points_6[int(1)];
                break;
            }
            if(_S707 != int(2))
            {
                _S705 = &points_6[int(0)];
                _S704 = nullptr;
                _S703 = nullptr;
                break;
            }
            float2  uv_35 = (pix_center_5 + make_float2 (0.0f, -1.0f) - _S708) / _S709;
            FixedArray<float, 4>  _S716 = dist_coeffs_15;
            float2  uv_u_16;
            bool _S717 = undistort_point_1(uv_35, &_S716, int(12), &uv_u_16);
            if(!_S717)
            {
                float3  _S718 = make_float3 (0.0f);
                _S707 = int(0);
                _S705 = &points_6[int(0)];
                _S704 = nullptr;
                _S703 = nullptr;
                normal_7 = _S718;
                break;
            }
            points_6[int(2)] = make_float3 (depths_4.z) * unproject_raydir_0(uv_u_16, camera_model_13, is_ray_depth_12);
            for(;;)
            {
                float2  uv_36 = (pix_center_5 + make_float2 (0.0f, 1.0f) - _S708) / _S709;
                FixedArray<float, 4>  _S719 = dist_coeffs_15;
                float2  uv_u_17;
                bool _S720 = undistort_point_1(uv_36, &_S719, int(12), &uv_u_17);
                if(!_S720)
                {
                    float3  _S721 = make_float3 (0.0f);
                    _S707 = int(0);
                    _S705 = nullptr;
                    normal_7 = _S721;
                    break;
                }
                points_6[int(3)] = make_float3 (depths_4.w) * unproject_raydir_0(uv_u_17, camera_model_13, is_ray_depth_12);
                _S707 = int(2);
                _S705 = &points_6[int(3)];
                break;
            }
            if(_S707 != int(2))
            {
                float3  * _S722 = _S705;
                _S705 = &points_6[int(0)];
                _S704 = _S722;
                _S703 = &points_6[int(2)];
                break;
            }
            float3  * _S723 = _S705;
            _S707 = int(1);
            _S705 = &points_6[int(0)];
            _S704 = _S723;
            _S703 = &points_6[int(2)];
            break;
        }
        if(_S707 != int(1))
        {
            break;
        }
        float3  normal_8 = cross_0(*_S706 - *_S705, - (*_S704 - *_S703));
        if((dot_0(normal_8, normal_8)) != 0.0f)
        {
            normal_7 = normal_8 / make_float3 (length_0(normal_8));
        }
        else
        {
            normal_7 = normal_8;
        }
        break;
    }
    return normal_7;
}

struct s_bwd_prop_depth_to_normal_Intermediates_1
{
    float2  _S724;
    bool _S725;
    float2  _S726;
    bool _S727;
    float2  _S728;
    bool _S729;
    float2  _S730;
    bool _S731;
};

inline __device__ void depth_to_normal_vjp_opencv(float2  pix_center_6, float4  intrins_12, FixedArray<float, 4>  dist_coeffs_16, int camera_model_14, bool is_ray_depth_13, float4  depths_5, float3  v_normal_2, float4  * v_depths_2)
{
    float2  _S732 = make_float2 (0.0f);
    s_bwd_prop_depth_to_normal_Intermediates_1 _S733;
    (&_S733)->_S724 = _S732;
    (&_S733)->_S725 = false;
    (&_S733)->_S726 = _S732;
    (&_S733)->_S727 = false;
    (&_S733)->_S728 = _S732;
    (&_S733)->_S729 = false;
    (&_S733)->_S730 = _S732;
    (&_S733)->_S731 = false;
    (&_S733)->_S724 = _S732;
    (&_S733)->_S725 = false;
    (&_S733)->_S726 = _S732;
    (&_S733)->_S727 = false;
    (&_S733)->_S728 = _S732;
    (&_S733)->_S729 = false;
    (&_S733)->_S730 = _S732;
    (&_S733)->_S731 = false;
    bool _S734 = (depths_5.x) == 0.0f;
    bool _runFlag_8;
    if(_S734)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.y) == 0.0f;
    }
    if(_runFlag_8)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.z) == 0.0f;
    }
    if(_runFlag_8)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.w) == 0.0f;
    }
    int _S735;
    if(!_runFlag_8)
    {
        float2  _S736 = float2 {intrins_12.z, intrins_12.w};
        float2  _S737 = float2 {intrins_12.x, intrins_12.y};
        float2  uv_37 = (pix_center_6 + make_float2 (-1.0f, -0.0f) - _S736) / _S737;
        float2  _S738 = _S732;
        FixedArray<float, 4>  _S739 = dist_coeffs_16;
        bool _S740 = undistort_point_1(uv_37, &_S739, int(12), &_S738);
        (&_S733)->_S724 = _S738;
        (&_S733)->_S725 = _S740;
        bool _S741 = !!_S740;
        if(_S741)
        {
            float2  uv_38 = (pix_center_6 + make_float2 (1.0f, -0.0f) - _S736) / _S737;
            float2  _S742 = _S732;
            FixedArray<float, 4>  _S743 = dist_coeffs_16;
            bool _S744 = undistort_point_1(uv_38, &_S743, int(12), &_S742);
            (&_S733)->_S726 = _S742;
            (&_S733)->_S727 = _S744;
            if(!!_S744)
            {
                _S735 = int(2);
            }
            else
            {
                _S735 = int(0);
            }
            if(_S735 != int(2))
            {
                _runFlag_8 = false;
            }
            else
            {
                _runFlag_8 = _S741;
            }
            if(_runFlag_8)
            {
                float2  uv_39 = (pix_center_6 + make_float2 (0.0f, -1.0f) - _S736) / _S737;
                float2  _S745 = _S732;
                FixedArray<float, 4>  _S746 = dist_coeffs_16;
                bool _S747 = undistort_point_1(uv_39, &_S746, int(12), &_S745);
                (&_S733)->_S728 = _S745;
                (&_S733)->_S729 = _S747;
                if(!_S747)
                {
                    _runFlag_8 = false;
                }
                if(_runFlag_8)
                {
                    float2  uv_40 = (pix_center_6 + make_float2 (0.0f, 1.0f) - _S736) / _S737;
                    float2  _S748 = _S732;
                    FixedArray<float, 4>  _S749 = dist_coeffs_16;
                    bool _S750 = undistort_point_1(uv_40, &_S749, int(12), &_S748);
                    (&_S733)->_S730 = _S748;
                    (&_S733)->_S731 = _S750;
                }
            }
        }
    }
    s_bwd_prop_depth_to_normal_Intermediates_1 _S751 = _S733;
    float3  _S752 = make_float3 (0.0f);
    if(_S734)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.y) == 0.0f;
    }
    if(_runFlag_8)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.z) == 0.0f;
    }
    if(_runFlag_8)
    {
        _runFlag_8 = true;
    }
    else
    {
        _runFlag_8 = (depths_5.w) == 0.0f;
    }
    bool _S753 = !_runFlag_8;
    bool _runFlag_9;
    bool _runFlag_10;
    bool _S754;
    bool _runFlag_11;
    bool _S755;
    bool _S756;
    FixedArray<float3 , 4>  points_7;
    float3  _S757;
    float3  _S758;
    float3  _S759;
    float3  _S760;
    float3  _S761;
    float3  _S762;
    float3  _S763;
    float3  _S764;
    float3  _S765;
    if(_S753)
    {
        bool _S766 = !!_S751._S725;
        if(_S766)
        {
            float3  _S767 = s_primal_ctx_unproject_raydir_0(_S751._S724, camera_model_14, is_ray_depth_13);
            float3  _S768 = make_float3 (depths_5.x) * _S767;
            bool _S769 = !!_S751._S727;
            if(_S769)
            {
                float3  _S770 = s_primal_ctx_unproject_raydir_0(_S751._S726, camera_model_14, is_ray_depth_13);
                float3  _S771 = make_float3 (depths_5.y) * _S770;
                _S735 = int(2);
                points_7[int(0)] = _S768;
                points_7[int(1)] = _S771;
                points_7[int(2)] = _S752;
                points_7[int(3)] = _S752;
                _S757 = _S770;
            }
            else
            {
                _S735 = int(0);
                points_7[int(0)] = _S768;
                points_7[int(1)] = _S752;
                points_7[int(2)] = _S752;
                points_7[int(3)] = _S752;
                _S757 = _S752;
            }
            if(_S735 != int(2))
            {
                _runFlag_8 = false;
            }
            else
            {
                _runFlag_8 = _S766;
                _S735 = int(0);
            }
            if(_runFlag_8)
            {
                if(!_S751._S729)
                {
                    _runFlag_9 = false;
                    _S735 = int(0);
                }
                else
                {
                    _runFlag_9 = _runFlag_8;
                }
                if(_runFlag_9)
                {
                    float3  _S772 = s_primal_ctx_unproject_raydir_0(_S751._S728, camera_model_14, is_ray_depth_13);
                    points_7[int(2)] = make_float3 (depths_5.z) * _S772;
                    bool _S773 = !!_S751._S731;
                    int _S774;
                    if(_S773)
                    {
                        float3  _S775 = s_primal_ctx_unproject_raydir_0(_S751._S730, camera_model_14, is_ray_depth_13);
                        points_7[int(3)] = make_float3 (depths_5.w) * _S775;
                        _S774 = int(2);
                        _S758 = _S775;
                    }
                    else
                    {
                        _S774 = int(0);
                        _S758 = _S752;
                    }
                    if(_S774 != int(2))
                    {
                        _runFlag_10 = false;
                        _S735 = _S774;
                    }
                    else
                    {
                        _runFlag_10 = _runFlag_9;
                    }
                    if(_runFlag_10)
                    {
                        _S735 = int(1);
                    }
                    _runFlag_10 = _S773;
                    _S759 = _S772;
                }
                else
                {
                    _runFlag_10 = false;
                    _S758 = _S752;
                    _S759 = _S752;
                }
            }
            else
            {
                _runFlag_9 = false;
                _runFlag_10 = false;
                _S758 = _S752;
                _S759 = _S752;
            }
            float3  _S776 = _S757;
            _S757 = _S758;
            _S758 = _S759;
            _S754 = _S769;
            _S759 = _S776;
            _S760 = _S767;
        }
        else
        {
            _S735 = int(0);
            points_7[int(0)] = _S752;
            points_7[int(1)] = _S752;
            points_7[int(2)] = _S752;
            points_7[int(3)] = _S752;
            _runFlag_8 = false;
            _runFlag_9 = false;
            _runFlag_10 = false;
            _S757 = _S752;
            _S758 = _S752;
            _S754 = false;
            _S759 = _S752;
            _S760 = _S752;
        }
        if(_S735 != int(1))
        {
            _runFlag_11 = false;
        }
        else
        {
            _runFlag_11 = _S753;
        }
        if(_runFlag_11)
        {
            float3  dx_3 = points_7[int(1)] - points_7[int(0)];
            float3  _S777 = - (points_7[int(3)] - points_7[int(2)]);
            float3  _S778 = s_primal_ctx_cross_0(dx_3, _S777);
            bool _S779 = (s_primal_ctx_dot_0(_S778, _S778)) != 0.0f;
            if(_S779)
            {
                float _S780 = length_0(_S778);
                float3  _S781 = make_float3 (_S780);
                _S761 = make_float3 (_S780 * _S780);
                _S762 = _S781;
            }
            else
            {
                _S761 = _S752;
                _S762 = _S752;
            }
            float3  _S782 = _S762;
            _S755 = _S779;
            _S762 = _S778;
            _S763 = _S782;
            _S764 = dx_3;
            _S765 = _S777;
        }
        else
        {
            _S755 = false;
            _S761 = _S752;
            _S762 = _S752;
            _S763 = _S752;
            _S764 = _S752;
            _S765 = _S752;
        }
        bool _S783 = _runFlag_8;
        bool _S784 = _runFlag_9;
        bool _S785 = _runFlag_10;
        float3  _S786 = _S757;
        float3  _S787 = _S758;
        bool _S788 = _S754;
        float3  _S789 = _S759;
        float3  _S790 = _S760;
        _runFlag_8 = _runFlag_11;
        _runFlag_9 = _S755;
        _S757 = _S761;
        _S758 = _S762;
        _S759 = _S763;
        _S760 = _S764;
        _S761 = _S765;
        _runFlag_10 = _S766;
        _S754 = _S783;
        _runFlag_11 = _S784;
        _S755 = _S785;
        _S762 = _S786;
        _S763 = _S787;
        _S756 = _S788;
        _S764 = _S789;
        _S765 = _S790;
    }
    else
    {
        _runFlag_8 = false;
        _runFlag_9 = false;
        _S757 = _S752;
        _S758 = _S752;
        _S759 = _S752;
        _S760 = _S752;
        _S761 = _S752;
        _runFlag_10 = false;
        _S754 = false;
        _runFlag_11 = false;
        _S755 = false;
        _S762 = _S752;
        _S763 = _S752;
        _S756 = false;
        _S764 = _S752;
        _S765 = _S752;
    }
    float4  _S791 = make_float4 (0.0f);
    float4  _S792;
    if(_S753)
    {
        if(_runFlag_8)
        {
            if(_runFlag_9)
            {
                float3  _S793 = v_normal_2 / _S757;
                float3  _S794 = _S758 * - _S793;
                float3  _S795 = _S759 * _S793;
                float _S796 = _S794.x + _S794.y + _S794.z;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S797;
                (&_S797)->primal_0 = _S758;
                (&_S797)->differential_0 = _S752;
                s_bwd_length_impl_0(&_S797, _S796);
                _S757 = _S795 + _S797.differential_0;
            }
            else
            {
                _S757 = v_normal_2;
            }
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S798;
            (&_S798)->primal_0 = _S758;
            (&_S798)->differential_0 = _S752;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S799;
            (&_S799)->primal_0 = _S758;
            (&_S799)->differential_0 = _S752;
            s_bwd_prop_dot_0(&_S798, &_S799, 0.0f);
            float3  _S800 = _S799.differential_0 + _S798.differential_0 + _S757;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S801;
            (&_S801)->primal_0 = _S760;
            (&_S801)->differential_0 = _S752;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S802;
            (&_S802)->primal_0 = _S761;
            (&_S802)->differential_0 = _S752;
            s_bwd_prop_cross_0(&_S801, &_S802, _S800);
            float3  s_diff_dy_T_3 = - _S802.differential_0;
            float3  _S803 = - s_diff_dy_T_3;
            float3  _S804 = - _S801.differential_0;
            FixedArray<float3 , 4>  _S805;
            _S805[int(0)] = _S752;
            _S805[int(1)] = _S752;
            _S805[int(2)] = _S752;
            _S805[int(3)] = _S752;
            _S805[int(2)] = _S803;
            _S805[int(3)] = s_diff_dy_T_3;
            _S805[int(0)] = _S804;
            _S805[int(1)] = _S801.differential_0;
            points_7[int(0)] = _S805[int(0)];
            points_7[int(1)] = _S805[int(1)];
            points_7[int(2)] = _S805[int(2)];
            points_7[int(3)] = _S805[int(3)];
        }
        else
        {
            points_7[int(0)] = _S752;
            points_7[int(1)] = _S752;
            points_7[int(2)] = _S752;
            points_7[int(3)] = _S752;
        }
        if(_runFlag_10)
        {
            if(_S754)
            {
                if(_runFlag_11)
                {
                    FixedArray<float3 , 4>  _S806 = points_7;
                    FixedArray<float3 , 4>  _S807 = points_7;
                    FixedArray<float3 , 4>  _S808 = points_7;
                    FixedArray<float3 , 4>  _S809 = points_7;
                    if(_S755)
                    {
                        float3  _S810 = _S762 * _S809[int(3)];
                        float _S811 = _S810.x + _S810.y + _S810.z;
                        float4  _S812 = _S791;
                        *&((&_S812)->w) = _S811;
                        points_7[int(0)] = _S806[int(0)];
                        points_7[int(1)] = _S807[int(1)];
                        points_7[int(2)] = _S808[int(2)];
                        points_7[int(3)] = _S752;
                        _S792 = _S812;
                    }
                    else
                    {
                        points_7[int(0)] = _S806[int(0)];
                        points_7[int(1)] = _S807[int(1)];
                        points_7[int(2)] = _S808[int(2)];
                        points_7[int(3)] = _S809[int(3)];
                        _S792 = _S791;
                    }
                    float3  _S813 = _S763 * points_7[int(2)];
                    float _S814 = _S813.x + _S813.y + _S813.z;
                    FixedArray<float3 , 4>  _S815 = points_7;
                    FixedArray<float3 , 4>  _S816 = points_7;
                    float4  _S817 = _S791;
                    *&((&_S817)->z) = _S814;
                    float4  _S818 = _S792 + _S817;
                    points_7[int(0)] = points_7[int(0)];
                    points_7[int(1)] = _S815[int(1)];
                    points_7[int(2)] = _S752;
                    points_7[int(3)] = _S816[int(3)];
                    _S792 = _S818;
                }
                else
                {
                    FixedArray<float3 , 4>  _S819 = points_7;
                    FixedArray<float3 , 4>  _S820 = points_7;
                    FixedArray<float3 , 4>  _S821 = points_7;
                    points_7[int(0)] = points_7[int(0)];
                    points_7[int(1)] = _S819[int(1)];
                    points_7[int(2)] = _S820[int(2)];
                    points_7[int(3)] = _S821[int(3)];
                    _S792 = _S791;
                }
            }
            else
            {
                FixedArray<float3 , 4>  _S822 = points_7;
                FixedArray<float3 , 4>  _S823 = points_7;
                FixedArray<float3 , 4>  _S824 = points_7;
                points_7[int(0)] = points_7[int(0)];
                points_7[int(1)] = _S822[int(1)];
                points_7[int(2)] = _S823[int(2)];
                points_7[int(3)] = _S824[int(3)];
                _S792 = _S791;
            }
            if(_S756)
            {
                FixedArray<float3 , 4>  _S825 = points_7;
                float3  _S826 = _S764 * points_7[int(1)];
                float _S827 = _S826.x + _S826.y + _S826.z;
                float4  _S828 = _S791;
                *&((&_S828)->y) = _S827;
                float4  _S829 = _S792 + _S828;
                points_7[int(0)] = _S752;
                points_7[int(1)] = _S752;
                points_7[int(2)] = _S752;
                points_7[int(3)] = _S752;
                _S757 = _S825[int(0)];
                _S792 = _S829;
            }
            else
            {
                FixedArray<float3 , 4>  _S830 = points_7;
                FixedArray<float3 , 4>  _S831 = points_7;
                FixedArray<float3 , 4>  _S832 = points_7;
                points_7[int(0)] = points_7[int(0)];
                points_7[int(1)] = _S830[int(1)];
                points_7[int(2)] = _S831[int(2)];
                points_7[int(3)] = _S832[int(3)];
                _S757 = _S752;
            }
            float3  _S833 = _S765 * (points_7[int(0)] + _S757);
            float _S834 = _S833.x + _S833.y + _S833.z;
            float4  _S835 = _S791;
            *&((&_S835)->x) = _S834;
            _S792 = _S792 + _S835;
        }
        else
        {
            _S792 = _S791;
        }
    }
    else
    {
        _S792 = _S791;
    }
    *v_depths_2 = _S792;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor_opencv(float2  pix_center_7, float4  intrins_13, FixedArray<float, 4>  dist_coeffs_17, int camera_model_15)
{
    float _S836;
    for(;;)
    {
        float2  uv_41 = (pix_center_7 - float2 {intrins_13.z, intrins_13.w}) / float2 {intrins_13.x, intrins_13.y};
        FixedArray<float, 4>  _S837 = dist_coeffs_17;
        float2  uv_u_18;
        bool _S838 = undistort_point_1(uv_41, &_S837, int(12), &uv_u_18);
        if(!_S838)
        {
            _S836 = 0.0f;
            break;
        }
        float3  raydir_9 = unproject_raydir_0(uv_u_18, camera_model_15, false);
        _S836 = float((F32_sign((raydir_9.z)))) / length_0(raydir_9);
        break;
    }
    return _S836;
}

inline __device__ float depth_normal_loss_opencv(float2  pix_center_8, float4  intrins_14, FixedArray<float, 4>  dist_coeffs_18, int camera_model_16, bool is_ray_depth_14, float4  depths_6, float3  gt_normal_2)
{
    float _S839;
    for(;;)
    {
        float3  _S840;
        float3  * _S841;
        float3  * _S842;
        float3  * _S843;
        float3  * _S844;
        int _S845;
        FixedArray<float3 , 5>  points_8;
        for(;;)
        {
            float2  _S846 = float2 {intrins_14.z, intrins_14.w};
            float2  _S847 = float2 {intrins_14.x, intrins_14.y};
            float2  uv_42 = (pix_center_8 + make_float2 (-1.0f, -0.0f) - _S846) / _S847;
            FixedArray<float, 4>  _S848 = dist_coeffs_18;
            float2  uv_u_19;
            bool _S849 = undistort_point_1(uv_42, &_S848, int(12), &uv_u_19);
            float3  _S850 = make_float3 (0.0f);
            if(!_S849)
            {
                _S845 = int(0);
                _S844 = nullptr;
                _S843 = nullptr;
                _S842 = nullptr;
                _S841 = nullptr;
                _S840 = _S850;
                break;
            }
            float3  raydir_10 = unproject_raydir_0(uv_u_19, camera_model_16, is_ray_depth_14);
            points_8[int(0)] = make_float3 (depths_6.x) * raydir_10;
            float2  uv_43 = (pix_center_8 + make_float2 (1.0f, -0.0f) - _S846) / _S847;
            FixedArray<float, 4>  _S851 = dist_coeffs_18;
            float2  uv_u_20;
            bool _S852 = undistort_point_1(uv_43, &_S851, int(12), &uv_u_20);
            if(!_S852)
            {
                _S845 = int(0);
                _S844 = nullptr;
                _S843 = &points_8[int(0)];
                _S842 = nullptr;
                _S841 = nullptr;
                _S840 = _S850;
                break;
            }
            float3  raydir_11 = unproject_raydir_0(uv_u_20, camera_model_16, is_ray_depth_14);
            points_8[int(1)] = make_float3 (depths_6.y) * raydir_11;
            float2  uv_44 = (pix_center_8 + make_float2 (0.0f, -1.0f) - _S846) / _S847;
            FixedArray<float, 4>  _S853 = dist_coeffs_18;
            float2  uv_u_21;
            bool _S854 = undistort_point_1(uv_44, &_S853, int(12), &uv_u_21);
            if(!_S854)
            {
                _S845 = int(0);
                _S844 = &points_8[int(1)];
                _S843 = &points_8[int(0)];
                _S842 = nullptr;
                _S841 = nullptr;
                _S840 = _S850;
                break;
            }
            float3  raydir_12 = unproject_raydir_0(uv_u_21, camera_model_16, is_ray_depth_14);
            points_8[int(2)] = make_float3 (depths_6.z) * raydir_12;
            float2  uv_45 = (pix_center_8 + make_float2 (0.0f, 1.0f) - _S846) / _S847;
            FixedArray<float, 4>  _S855 = dist_coeffs_18;
            float2  uv_u_22;
            bool _S856 = undistort_point_1(uv_45, &_S855, int(12), &uv_u_22);
            if(!_S856)
            {
                _S845 = int(0);
                _S844 = &points_8[int(1)];
                _S843 = &points_8[int(0)];
                _S842 = nullptr;
                _S841 = &points_8[int(2)];
                _S840 = _S850;
                break;
            }
            float3  raydir_13 = unproject_raydir_0(uv_u_22, camera_model_16, is_ray_depth_14);
            points_8[int(3)] = make_float3 (depths_6.w) * raydir_13;
            float2  uv_46 = (pix_center_8 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S846) / _S847;
            FixedArray<float, 4>  _S857 = dist_coeffs_18;
            float2  uv_u_23;
            bool _S858 = undistort_point_1(uv_46, &_S857, int(12), &uv_u_23);
            if(!_S858)
            {
                _S845 = int(0);
                _S844 = &points_8[int(1)];
                _S843 = &points_8[int(0)];
                _S842 = &points_8[int(3)];
                _S841 = &points_8[int(2)];
                _S840 = _S850;
                break;
            }
            float3  raydir_14 = unproject_raydir_0(uv_u_23, camera_model_16, is_ray_depth_14);
            _S845 = int(1);
            _S844 = &points_8[int(1)];
            _S843 = &points_8[int(0)];
            _S842 = &points_8[int(3)];
            _S841 = &points_8[int(2)];
            _S840 = raydir_14;
            break;
        }
        if(_S845 != int(1))
        {
            _S839 = 0.0f;
            break;
        }
        float3  normal_9 = cross_0(*_S844 - *_S843, - (*_S842 - *_S841));
        float3  normal_10;
        if((dot_0(normal_9, normal_9)) != 0.0f)
        {
            normal_10 = normalize_0(normal_9);
        }
        else
        {
            normal_10 = normal_9;
        }
        float3  _S859;
        if((dot_0(gt_normal_2, gt_normal_2)) != 0.0f)
        {
            _S859 = normalize_0(gt_normal_2);
        }
        else
        {
            _S859 = gt_normal_2;
        }
        _S839 = (1.0f - dot_0(normal_10, _S859) + 0.00100000004749745f) / ((F32_max((dot_0(normal_10, - normalize_0(_S840))), (0.0f))) + 0.00100000004749745f);
        break;
    }
    return _S839;
}

struct s_bwd_prop_depth_normal_loss_Intermediates_1
{
    float2  _S860;
    bool _S861;
    float2  _S862;
    bool _S863;
    float2  _S864;
    bool _S865;
    float2  _S866;
    bool _S867;
    float2  _S868;
    bool _S869;
};

inline __device__ void depth_normal_loss_vjp_opencv(float2  pix_center_9, float4  intrins_15, FixedArray<float, 4>  dist_coeffs_19, int camera_model_17, bool is_ray_depth_15, float4  depths_7, float3  gt_normal_3, float v_loss_1, float4  * v_depths_3, float3  * v_gt_normal_1)
{
    float2  _S870 = make_float2 (0.0f);
    s_bwd_prop_depth_normal_loss_Intermediates_1 _S871;
    (&_S871)->_S860 = _S870;
    (&_S871)->_S861 = false;
    (&_S871)->_S862 = _S870;
    (&_S871)->_S863 = false;
    (&_S871)->_S864 = _S870;
    (&_S871)->_S865 = false;
    (&_S871)->_S866 = _S870;
    (&_S871)->_S867 = false;
    (&_S871)->_S868 = _S870;
    (&_S871)->_S869 = false;
    (&_S871)->_S862 = _S870;
    (&_S871)->_S863 = false;
    (&_S871)->_S864 = _S870;
    (&_S871)->_S865 = false;
    (&_S871)->_S866 = _S870;
    (&_S871)->_S867 = false;
    (&_S871)->_S868 = _S870;
    (&_S871)->_S869 = false;
    float2  _S872 = float2 {intrins_15.z, intrins_15.w};
    float2  _S873 = float2 {intrins_15.x, intrins_15.y};
    float2  uv_47 = (pix_center_9 + make_float2 (-1.0f, -0.0f) - _S872) / _S873;
    float2  _S874 = _S870;
    FixedArray<float, 4>  _S875 = dist_coeffs_19;
    bool _S876 = undistort_point_1(uv_47, &_S875, int(12), &_S874);
    (&_S871)->_S860 = _S874;
    (&_S871)->_S861 = _S876;
    bool _S877 = !!_S876;
    bool _runFlag_12;
    if(_S877)
    {
        float2  uv_48 = (pix_center_9 + make_float2 (1.0f, -0.0f) - _S872) / _S873;
        float2  _S878 = _S870;
        FixedArray<float, 4>  _S879 = dist_coeffs_19;
        bool _S880 = undistort_point_1(uv_48, &_S879, int(12), &_S878);
        (&_S871)->_S862 = _S878;
        (&_S871)->_S863 = _S880;
        if(!_S880)
        {
            _runFlag_12 = false;
        }
        else
        {
            _runFlag_12 = _S877;
        }
        if(_runFlag_12)
        {
            float2  uv_49 = (pix_center_9 + make_float2 (0.0f, -1.0f) - _S872) / _S873;
            float2  _S881 = _S870;
            FixedArray<float, 4>  _S882 = dist_coeffs_19;
            bool _S883 = undistort_point_1(uv_49, &_S882, int(12), &_S881);
            (&_S871)->_S864 = _S881;
            (&_S871)->_S865 = _S883;
            if(!_S883)
            {
                _runFlag_12 = false;
            }
            if(_runFlag_12)
            {
                float2  uv_50 = (pix_center_9 + make_float2 (0.0f, 1.0f) - _S872) / _S873;
                float2  _S884 = _S870;
                FixedArray<float, 4>  _S885 = dist_coeffs_19;
                bool _S886 = undistort_point_1(uv_50, &_S885, int(12), &_S884);
                (&_S871)->_S866 = _S884;
                (&_S871)->_S867 = _S886;
                if(!_S886)
                {
                    _runFlag_12 = false;
                }
                if(_runFlag_12)
                {
                    float2  uv_51 = (pix_center_9 - _S872) / _S873;
                    float2  _S887 = _S870;
                    FixedArray<float, 4>  _S888 = dist_coeffs_19;
                    bool _S889 = undistort_point_1(uv_51, &_S888, int(12), &_S887);
                    (&_S871)->_S868 = _S887;
                    (&_S871)->_S869 = _S889;
                }
            }
        }
    }
    s_bwd_prop_depth_normal_loss_Intermediates_1 _S890 = _S871;
    float3  _S891 = make_float3 (0.0f);
    bool _S892 = !!_S871._S861;
    bool _runFlag_13;
    bool _runFlag_14;
    bool _runFlag_15;
    int _S893;
    float3  raydir_15;
    float3  _S894;
    float3  _S895;
    float3  _S896;
    float3  _S897;
    FixedArray<float3 , 5>  points_9;
    if(_S892)
    {
        float3  _S898 = s_primal_ctx_unproject_raydir_0(_S890._S860, camera_model_17, is_ray_depth_15);
        float3  _S899 = make_float3 (depths_7.x) * _S898;
        if(!_S890._S863)
        {
            _runFlag_12 = false;
        }
        else
        {
            _runFlag_12 = _S892;
        }
        if(_runFlag_12)
        {
            float3  _S900 = s_primal_ctx_unproject_raydir_0(_S890._S862, camera_model_17, is_ray_depth_15);
            float3  _S901 = make_float3 (depths_7.y) * _S900;
            if(!_S890._S865)
            {
                _runFlag_13 = false;
            }
            else
            {
                _runFlag_13 = _runFlag_12;
            }
            if(_runFlag_13)
            {
                float3  _S902 = s_primal_ctx_unproject_raydir_0(_S890._S864, camera_model_17, is_ray_depth_15);
                float3  _S903 = make_float3 (depths_7.z) * _S902;
                if(!_S890._S867)
                {
                    _runFlag_14 = false;
                }
                else
                {
                    _runFlag_14 = _runFlag_13;
                }
                if(_runFlag_14)
                {
                    float3  _S904 = s_primal_ctx_unproject_raydir_0(_S890._S866, camera_model_17, is_ray_depth_15);
                    float3  _S905 = make_float3 (depths_7.w) * _S904;
                    if(!_S890._S869)
                    {
                        _runFlag_15 = false;
                    }
                    else
                    {
                        _runFlag_15 = _runFlag_14;
                    }
                    if(_runFlag_15)
                    {
                        float3  _S906 = s_primal_ctx_unproject_raydir_0(_S890._S868, camera_model_17, is_ray_depth_15);
                        _S893 = int(1);
                        raydir_15 = _S906;
                    }
                    else
                    {
                        _S893 = int(0);
                        raydir_15 = _S904;
                    }
                    points_9[int(0)] = _S899;
                    points_9[int(1)] = _S901;
                    points_9[int(2)] = _S903;
                    points_9[int(3)] = _S905;
                    points_9[int(4)] = _S891;
                    _S894 = _S904;
                }
                else
                {
                    _S893 = int(0);
                    raydir_15 = _S902;
                    points_9[int(0)] = _S899;
                    points_9[int(1)] = _S901;
                    points_9[int(2)] = _S903;
                    points_9[int(3)] = _S891;
                    points_9[int(4)] = _S891;
                    _S894 = _S891;
                }
                _S895 = _S902;
            }
            else
            {
                _S893 = int(0);
                raydir_15 = _S900;
                points_9[int(0)] = _S899;
                points_9[int(1)] = _S901;
                points_9[int(2)] = _S891;
                points_9[int(3)] = _S891;
                points_9[int(4)] = _S891;
                _runFlag_14 = false;
                _S894 = _S891;
                _S895 = _S891;
            }
            _S896 = _S900;
        }
        else
        {
            _S893 = int(0);
            raydir_15 = _S898;
            points_9[int(0)] = _S899;
            points_9[int(1)] = _S891;
            points_9[int(2)] = _S891;
            points_9[int(3)] = _S891;
            points_9[int(4)] = _S891;
            _runFlag_13 = false;
            _runFlag_14 = false;
            _S894 = _S891;
            _S895 = _S891;
            _S896 = _S891;
        }
        _S897 = _S898;
    }
    else
    {
        _S893 = int(0);
        points_9[int(0)] = _S891;
        points_9[int(1)] = _S891;
        points_9[int(2)] = _S891;
        points_9[int(3)] = _S891;
        points_9[int(4)] = _S891;
        _runFlag_12 = false;
        _runFlag_13 = false;
        _runFlag_14 = false;
        _S894 = _S891;
        _S895 = _S891;
        _S896 = _S891;
        _S897 = _S891;
    }
    bool _S907 = !(_S893 != int(1));
    bool _S908;
    float3  normal_11;
    float3  _S909;
    float3  _S910;
    float3  _S911;
    float3  _S912;
    float _S913;
    float _S914;
    float _S915;
    float _S916;
    if(_S907)
    {
        float3  dx_4 = points_9[int(1)] - points_9[int(0)];
        float3  _S917 = - (points_9[int(3)] - points_9[int(2)]);
        float3  _S918 = s_primal_ctx_cross_0(dx_4, _S917);
        bool _S919 = (s_primal_ctx_dot_0(_S918, _S918)) != 0.0f;
        if(_S919)
        {
            normal_11 = normalize_0(_S918);
        }
        else
        {
            normal_11 = _S918;
        }
        bool _S920 = (s_primal_ctx_dot_0(gt_normal_3, gt_normal_3)) != 0.0f;
        if(_S920)
        {
            _S909 = normalize_0(gt_normal_3);
        }
        else
        {
            _S909 = gt_normal_3;
        }
        float3  _S921 = - normalize_0(raydir_15);
        float _S922 = s_primal_ctx_dot_0(normal_11, _S921);
        float _S923 = 1.0f - s_primal_ctx_dot_0(normal_11, _S909) + 0.00100000004749745f;
        float _S924 = (F32_max((_S922), (0.0f))) + 0.00100000004749745f;
        _S913 = _S924 * _S924;
        _S914 = _S923;
        _S915 = _S924;
        _S916 = _S922;
        raydir_15 = normal_11;
        normal_11 = _S921;
        _runFlag_15 = _S920;
        _S908 = _S919;
        _S910 = _S918;
        _S911 = dx_4;
        _S912 = _S917;
    }
    else
    {
        _S913 = 0.0f;
        _S914 = 0.0f;
        _S915 = 0.0f;
        _S916 = 0.0f;
        raydir_15 = _S891;
        normal_11 = _S891;
        _S909 = _S891;
        _runFlag_15 = false;
        _S908 = false;
        _S910 = _S891;
        _S911 = _S891;
        _S912 = _S891;
    }
    float4  _S925 = make_float4 (0.0f);
    if(_S907)
    {
        float _S926 = v_loss_1 / _S913;
        float _S927 = _S914 * - _S926;
        float s_diff_num_T_1 = _S915 * _S926;
        DiffPair_float_0 _S928;
        (&_S928)->primal_0 = _S916;
        (&_S928)->differential_0 = 0.0f;
        DiffPair_float_0 _S929;
        (&_S929)->primal_0 = 0.0f;
        (&_S929)->differential_0 = 0.0f;
        _d_max_0(&_S928, &_S929, _S927);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S930;
        (&_S930)->primal_0 = raydir_15;
        (&_S930)->differential_0 = _S891;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S931;
        (&_S931)->primal_0 = normal_11;
        (&_S931)->differential_0 = _S891;
        s_bwd_prop_dot_0(&_S930, &_S931, _S928.differential_0);
        float _S932 = - s_diff_num_T_1;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S933;
        (&_S933)->primal_0 = raydir_15;
        (&_S933)->differential_0 = _S891;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S934;
        (&_S934)->primal_0 = _S909;
        (&_S934)->differential_0 = _S891;
        s_bwd_prop_dot_0(&_S933, &_S934, _S932);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S935 = _S934;
        float3  _S936 = _S930.differential_0 + _S933.differential_0;
        if(_runFlag_15)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S937;
            (&_S937)->primal_0 = gt_normal_3;
            (&_S937)->differential_0 = _S891;
            s_bwd_normalize_impl_0(&_S937, _S935.differential_0);
            raydir_15 = _S937.differential_0;
        }
        else
        {
            raydir_15 = _S935.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S938;
        (&_S938)->primal_0 = gt_normal_3;
        (&_S938)->differential_0 = _S891;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S939;
        (&_S939)->primal_0 = gt_normal_3;
        (&_S939)->differential_0 = _S891;
        s_bwd_prop_dot_0(&_S938, &_S939, 0.0f);
        float3  _S940 = _S939.differential_0 + _S938.differential_0 + raydir_15;
        if(_S908)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S941;
            (&_S941)->primal_0 = _S910;
            (&_S941)->differential_0 = _S891;
            s_bwd_normalize_impl_0(&_S941, _S936);
            raydir_15 = _S941.differential_0;
        }
        else
        {
            raydir_15 = _S936;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S942;
        (&_S942)->primal_0 = _S910;
        (&_S942)->differential_0 = _S891;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S943;
        (&_S943)->primal_0 = _S910;
        (&_S943)->differential_0 = _S891;
        s_bwd_prop_dot_0(&_S942, &_S943, 0.0f);
        float3  _S944 = _S943.differential_0 + _S942.differential_0 + raydir_15;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S945;
        (&_S945)->primal_0 = _S911;
        (&_S945)->differential_0 = _S891;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S946;
        (&_S946)->primal_0 = _S912;
        (&_S946)->differential_0 = _S891;
        s_bwd_prop_cross_0(&_S945, &_S946, _S944);
        float3  s_diff_dy_T_4 = - _S946.differential_0;
        float3  _S947 = - s_diff_dy_T_4;
        float3  _S948 = - _S945.differential_0;
        FixedArray<float3 , 5>  _S949;
        _S949[int(0)] = _S891;
        _S949[int(1)] = _S891;
        _S949[int(2)] = _S891;
        _S949[int(3)] = _S891;
        _S949[int(4)] = _S891;
        _S949[int(2)] = _S947;
        _S949[int(3)] = s_diff_dy_T_4;
        _S949[int(0)] = _S948;
        _S949[int(1)] = _S945.differential_0;
        points_9[int(0)] = _S949[int(0)];
        points_9[int(1)] = _S949[int(1)];
        points_9[int(2)] = _S949[int(2)];
        points_9[int(3)] = _S949[int(3)];
        points_9[int(4)] = _S949[int(4)];
        raydir_15 = _S940;
    }
    else
    {
        points_9[int(0)] = _S891;
        points_9[int(1)] = _S891;
        points_9[int(2)] = _S891;
        points_9[int(3)] = _S891;
        points_9[int(4)] = _S891;
        raydir_15 = _S891;
    }
    float4  _S950;
    if(_S892)
    {
        if(_runFlag_12)
        {
            if(_runFlag_13)
            {
                if(_runFlag_14)
                {
                    FixedArray<float3 , 5>  _S951 = points_9;
                    FixedArray<float3 , 5>  _S952 = points_9;
                    FixedArray<float3 , 5>  _S953 = points_9;
                    float3  _S954 = _S894 * points_9[int(3)];
                    float _S955 = _S954.x + _S954.y + _S954.z;
                    float4  _S956 = _S925;
                    *&((&_S956)->w) = _S955;
                    points_9[int(0)] = _S891;
                    points_9[int(1)] = _S891;
                    points_9[int(2)] = _S891;
                    points_9[int(3)] = _S891;
                    points_9[int(4)] = _S891;
                    _S894 = _S953[int(2)];
                    normal_11 = _S951[int(0)];
                    _S909 = _S952[int(1)];
                    _S950 = _S956;
                }
                else
                {
                    FixedArray<float3 , 5>  _S957 = points_9;
                    FixedArray<float3 , 5>  _S958 = points_9;
                    FixedArray<float3 , 5>  _S959 = points_9;
                    FixedArray<float3 , 5>  _S960 = points_9;
                    points_9[int(0)] = points_9[int(0)];
                    points_9[int(1)] = _S957[int(1)];
                    points_9[int(2)] = _S958[int(2)];
                    points_9[int(3)] = _S959[int(3)];
                    points_9[int(4)] = _S960[int(4)];
                    _S894 = _S891;
                    normal_11 = _S891;
                    _S909 = _S891;
                    _S950 = _S925;
                }
                float3  _S961 = _S895 * (points_9[int(2)] + _S894);
                float _S962 = _S961.x + _S961.y + _S961.z;
                float3  _S963 = points_9[int(0)] + normal_11;
                float3  _S964 = points_9[int(1)] + _S909;
                float4  _S965 = _S925;
                *&((&_S965)->z) = _S962;
                float4  _S966 = _S950 + _S965;
                points_9[int(0)] = _S891;
                points_9[int(1)] = _S891;
                points_9[int(2)] = _S891;
                points_9[int(3)] = _S891;
                points_9[int(4)] = _S891;
                _S894 = _S964;
                _S895 = _S963;
                _S950 = _S966;
            }
            else
            {
                FixedArray<float3 , 5>  _S967 = points_9;
                FixedArray<float3 , 5>  _S968 = points_9;
                FixedArray<float3 , 5>  _S969 = points_9;
                FixedArray<float3 , 5>  _S970 = points_9;
                points_9[int(0)] = points_9[int(0)];
                points_9[int(1)] = _S967[int(1)];
                points_9[int(2)] = _S968[int(2)];
                points_9[int(3)] = _S969[int(3)];
                points_9[int(4)] = _S970[int(4)];
                _S894 = _S891;
                _S895 = _S891;
                _S950 = _S925;
            }
            float3  _S971 = _S896 * (points_9[int(1)] + _S894);
            float _S972 = _S971.x + _S971.y + _S971.z;
            float3  _S973 = points_9[int(0)] + _S895;
            float4  _S974 = _S925;
            *&((&_S974)->y) = _S972;
            float4  _S975 = _S950 + _S974;
            points_9[int(0)] = _S891;
            points_9[int(1)] = _S891;
            points_9[int(2)] = _S891;
            points_9[int(3)] = _S891;
            points_9[int(4)] = _S891;
            _S894 = _S973;
            _S950 = _S975;
        }
        else
        {
            FixedArray<float3 , 5>  _S976 = points_9;
            FixedArray<float3 , 5>  _S977 = points_9;
            FixedArray<float3 , 5>  _S978 = points_9;
            FixedArray<float3 , 5>  _S979 = points_9;
            points_9[int(0)] = points_9[int(0)];
            points_9[int(1)] = _S976[int(1)];
            points_9[int(2)] = _S977[int(2)];
            points_9[int(3)] = _S978[int(3)];
            points_9[int(4)] = _S979[int(4)];
            _S894 = _S891;
            _S950 = _S925;
        }
        float3  _S980 = _S897 * (points_9[int(0)] + _S894);
        float _S981 = _S980.x + _S980.y + _S980.z;
        float4  _S982 = _S925;
        *&((&_S982)->x) = _S981;
        _S950 = _S950 + _S982;
    }
    else
    {
        _S950 = _S925;
    }
    *v_depths_3 = _S950;
    *v_gt_normal_1 = raydir_15;
    return;
}

inline __device__ float3  generate_ray_d2n_prism(float2  pix_pos_6, float4  intrins_16, FixedArray<float, 8>  dist_coeffs_20, int camera_model_18, bool is_ray_depth_16)
{
    float3  _S983;
    for(;;)
    {
        float2  uv_52 = (pix_pos_6 - float2 {intrins_16.z, intrins_16.w}) / float2 {intrins_16.x, intrins_16.y};
        FixedArray<float, 8>  _S984 = dist_coeffs_20;
        float2  uv_u_24;
        bool _S985 = undistort_point_2(uv_52, &_S984, int(12), &uv_u_24);
        if(!_S985)
        {
            int3  _S986 = make_int3 (int(0));
            float3  _S987 = make_float3 ((float)_S986.x, (float)_S986.y, (float)_S986.z);
            _S983 = _S987;
            break;
        }
        _S983 = unproject_raydir_0(uv_u_24, camera_model_18, is_ray_depth_16);
        break;
    }
    return _S983;
}

inline __device__ float3  depth_to_point_prism(float2  pix_pos_7, float4  intrins_17, FixedArray<float, 8>  dist_coeffs_21, int camera_model_19, bool is_ray_depth_17, float depth_6)
{
    float3  _S988;
    for(;;)
    {
        float2  uv_53 = (pix_pos_7 - float2 {intrins_17.z, intrins_17.w}) / float2 {intrins_17.x, intrins_17.y};
        FixedArray<float, 8>  _S989 = dist_coeffs_21;
        float2  uv_u_25;
        bool _S990 = undistort_point_2(uv_53, &_S989, int(12), &uv_u_25);
        if(!_S990)
        {
            _S988 = make_float3 (0.0f);
            break;
        }
        _S988 = make_float3 (depth_6) * unproject_raydir_0(uv_u_25, camera_model_19, is_ray_depth_17);
        break;
    }
    return _S988;
}

struct s_bwd_prop_depth_to_point_Intermediates_2
{
    float2  _S991;
    bool _S992;
};

inline __device__ float depth_to_point_vjp_prism(float2  pix_pos_8, float4  intrins_18, FixedArray<float, 8>  dist_coeffs_22, int camera_model_20, bool is_ray_depth_18, float depth_7, float3  v_point_2)
{
    float2  _S993 = make_float2 (0.0f);
    s_bwd_prop_depth_to_point_Intermediates_2 _S994;
    (&_S994)->_S991 = _S993;
    (&_S994)->_S992 = false;
    float2  uv_54 = (pix_pos_8 - float2 {intrins_18.z, intrins_18.w}) / float2 {intrins_18.x, intrins_18.y};
    float2  _S995 = _S993;
    FixedArray<float, 8>  _S996 = dist_coeffs_22;
    bool _S997 = undistort_point_2(uv_54, &_S996, int(12), &_S995);
    (&_S994)->_S991 = _S995;
    (&_S994)->_S992 = _S997;
    s_bwd_prop_depth_to_point_Intermediates_2 _S998 = _S994;
    float3  _S999 = make_float3 (0.0f);
    bool _S1000 = !!_S994._S992;
    float3  _S1001;
    if(_S1000)
    {
        _S1001 = s_primal_ctx_unproject_raydir_0(_S998._S991, camera_model_20, is_ray_depth_18);
    }
    else
    {
        _S1001 = _S999;
    }
    if(_S1000)
    {
        _S1001 = _S1001 * v_point_2;
    }
    else
    {
        _S1001 = _S999;
    }
    return _S1001.x + _S1001.y + _S1001.z;
}

inline __device__ float3  depth_to_normal_prism(float2  pix_center_10, float4  intrins_19, FixedArray<float, 8>  dist_coeffs_23, int camera_model_21, bool is_ray_depth_19, float4  depths_8)
{
    float3  normal_12;
    for(;;)
    {
        bool _S1002;
        if((depths_8.x) == 0.0f)
        {
            _S1002 = true;
        }
        else
        {
            _S1002 = (depths_8.y) == 0.0f;
        }
        if(_S1002)
        {
            _S1002 = true;
        }
        else
        {
            _S1002 = (depths_8.z) == 0.0f;
        }
        if(_S1002)
        {
            _S1002 = true;
        }
        else
        {
            _S1002 = (depths_8.w) == 0.0f;
        }
        if(_S1002)
        {
            normal_12 = make_float3 (0.0f);
            break;
        }
        float3  * _S1003;
        float3  * _S1004;
        float3  * _S1005;
        float3  * _S1006;
        int _S1007;
        FixedArray<float3 , 4>  points_10;
        for(;;)
        {
            float2  _S1008 = float2 {intrins_19.z, intrins_19.w};
            float2  _S1009 = float2 {intrins_19.x, intrins_19.y};
            float2  uv_55 = (pix_center_10 + make_float2 (-1.0f, -0.0f) - _S1008) / _S1009;
            FixedArray<float, 8>  _S1010 = dist_coeffs_23;
            float2  uv_u_26;
            bool _S1011 = undistort_point_2(uv_55, &_S1010, int(12), &uv_u_26);
            if(!_S1011)
            {
                float3  _S1012 = make_float3 (0.0f);
                _S1007 = int(0);
                _S1006 = nullptr;
                _S1005 = nullptr;
                _S1004 = nullptr;
                _S1003 = nullptr;
                normal_12 = _S1012;
                break;
            }
            points_10[int(0)] = make_float3 (depths_8.x) * unproject_raydir_0(uv_u_26, camera_model_21, is_ray_depth_19);
            for(;;)
            {
                float2  uv_56 = (pix_center_10 + make_float2 (1.0f, -0.0f) - _S1008) / _S1009;
                FixedArray<float, 8>  _S1013 = dist_coeffs_23;
                float2  uv_u_27;
                bool _S1014 = undistort_point_2(uv_56, &_S1013, int(12), &uv_u_27);
                if(!_S1014)
                {
                    float3  _S1015 = make_float3 (0.0f);
                    _S1007 = int(0);
                    _S1006 = nullptr;
                    normal_12 = _S1015;
                    break;
                }
                points_10[int(1)] = make_float3 (depths_8.y) * unproject_raydir_0(uv_u_27, camera_model_21, is_ray_depth_19);
                _S1007 = int(2);
                _S1006 = &points_10[int(1)];
                break;
            }
            if(_S1007 != int(2))
            {
                _S1005 = &points_10[int(0)];
                _S1004 = nullptr;
                _S1003 = nullptr;
                break;
            }
            float2  uv_57 = (pix_center_10 + make_float2 (0.0f, -1.0f) - _S1008) / _S1009;
            FixedArray<float, 8>  _S1016 = dist_coeffs_23;
            float2  uv_u_28;
            bool _S1017 = undistort_point_2(uv_57, &_S1016, int(12), &uv_u_28);
            if(!_S1017)
            {
                float3  _S1018 = make_float3 (0.0f);
                _S1007 = int(0);
                _S1005 = &points_10[int(0)];
                _S1004 = nullptr;
                _S1003 = nullptr;
                normal_12 = _S1018;
                break;
            }
            points_10[int(2)] = make_float3 (depths_8.z) * unproject_raydir_0(uv_u_28, camera_model_21, is_ray_depth_19);
            for(;;)
            {
                float2  uv_58 = (pix_center_10 + make_float2 (0.0f, 1.0f) - _S1008) / _S1009;
                FixedArray<float, 8>  _S1019 = dist_coeffs_23;
                float2  uv_u_29;
                bool _S1020 = undistort_point_2(uv_58, &_S1019, int(12), &uv_u_29);
                if(!_S1020)
                {
                    float3  _S1021 = make_float3 (0.0f);
                    _S1007 = int(0);
                    _S1005 = nullptr;
                    normal_12 = _S1021;
                    break;
                }
                points_10[int(3)] = make_float3 (depths_8.w) * unproject_raydir_0(uv_u_29, camera_model_21, is_ray_depth_19);
                _S1007 = int(2);
                _S1005 = &points_10[int(3)];
                break;
            }
            if(_S1007 != int(2))
            {
                float3  * _S1022 = _S1005;
                _S1005 = &points_10[int(0)];
                _S1004 = _S1022;
                _S1003 = &points_10[int(2)];
                break;
            }
            float3  * _S1023 = _S1005;
            _S1007 = int(1);
            _S1005 = &points_10[int(0)];
            _S1004 = _S1023;
            _S1003 = &points_10[int(2)];
            break;
        }
        if(_S1007 != int(1))
        {
            break;
        }
        float3  normal_13 = cross_0(*_S1006 - *_S1005, - (*_S1004 - *_S1003));
        if((dot_0(normal_13, normal_13)) != 0.0f)
        {
            normal_12 = normal_13 / make_float3 (length_0(normal_13));
        }
        else
        {
            normal_12 = normal_13;
        }
        break;
    }
    return normal_12;
}

struct s_bwd_prop_depth_to_normal_Intermediates_2
{
    float2  _S1024;
    bool _S1025;
    float2  _S1026;
    bool _S1027;
    float2  _S1028;
    bool _S1029;
    float2  _S1030;
    bool _S1031;
};

inline __device__ void depth_to_normal_vjp_prism(float2  pix_center_11, float4  intrins_20, FixedArray<float, 8>  dist_coeffs_24, int camera_model_22, bool is_ray_depth_20, float4  depths_9, float3  v_normal_3, float4  * v_depths_4)
{
    float2  _S1032 = make_float2 (0.0f);
    s_bwd_prop_depth_to_normal_Intermediates_2 _S1033;
    (&_S1033)->_S1024 = _S1032;
    (&_S1033)->_S1025 = false;
    (&_S1033)->_S1026 = _S1032;
    (&_S1033)->_S1027 = false;
    (&_S1033)->_S1028 = _S1032;
    (&_S1033)->_S1029 = false;
    (&_S1033)->_S1030 = _S1032;
    (&_S1033)->_S1031 = false;
    (&_S1033)->_S1024 = _S1032;
    (&_S1033)->_S1025 = false;
    (&_S1033)->_S1026 = _S1032;
    (&_S1033)->_S1027 = false;
    (&_S1033)->_S1028 = _S1032;
    (&_S1033)->_S1029 = false;
    (&_S1033)->_S1030 = _S1032;
    (&_S1033)->_S1031 = false;
    bool _S1034 = (depths_9.x) == 0.0f;
    bool _runFlag_16;
    if(_S1034)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.y) == 0.0f;
    }
    if(_runFlag_16)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.z) == 0.0f;
    }
    if(_runFlag_16)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.w) == 0.0f;
    }
    int _S1035;
    if(!_runFlag_16)
    {
        float2  _S1036 = float2 {intrins_20.z, intrins_20.w};
        float2  _S1037 = float2 {intrins_20.x, intrins_20.y};
        float2  uv_59 = (pix_center_11 + make_float2 (-1.0f, -0.0f) - _S1036) / _S1037;
        float2  _S1038 = _S1032;
        FixedArray<float, 8>  _S1039 = dist_coeffs_24;
        bool _S1040 = undistort_point_2(uv_59, &_S1039, int(12), &_S1038);
        (&_S1033)->_S1024 = _S1038;
        (&_S1033)->_S1025 = _S1040;
        bool _S1041 = !!_S1040;
        if(_S1041)
        {
            float2  uv_60 = (pix_center_11 + make_float2 (1.0f, -0.0f) - _S1036) / _S1037;
            float2  _S1042 = _S1032;
            FixedArray<float, 8>  _S1043 = dist_coeffs_24;
            bool _S1044 = undistort_point_2(uv_60, &_S1043, int(12), &_S1042);
            (&_S1033)->_S1026 = _S1042;
            (&_S1033)->_S1027 = _S1044;
            if(!!_S1044)
            {
                _S1035 = int(2);
            }
            else
            {
                _S1035 = int(0);
            }
            if(_S1035 != int(2))
            {
                _runFlag_16 = false;
            }
            else
            {
                _runFlag_16 = _S1041;
            }
            if(_runFlag_16)
            {
                float2  uv_61 = (pix_center_11 + make_float2 (0.0f, -1.0f) - _S1036) / _S1037;
                float2  _S1045 = _S1032;
                FixedArray<float, 8>  _S1046 = dist_coeffs_24;
                bool _S1047 = undistort_point_2(uv_61, &_S1046, int(12), &_S1045);
                (&_S1033)->_S1028 = _S1045;
                (&_S1033)->_S1029 = _S1047;
                if(!_S1047)
                {
                    _runFlag_16 = false;
                }
                if(_runFlag_16)
                {
                    float2  uv_62 = (pix_center_11 + make_float2 (0.0f, 1.0f) - _S1036) / _S1037;
                    float2  _S1048 = _S1032;
                    FixedArray<float, 8>  _S1049 = dist_coeffs_24;
                    bool _S1050 = undistort_point_2(uv_62, &_S1049, int(12), &_S1048);
                    (&_S1033)->_S1030 = _S1048;
                    (&_S1033)->_S1031 = _S1050;
                }
            }
        }
    }
    s_bwd_prop_depth_to_normal_Intermediates_2 _S1051 = _S1033;
    float3  _S1052 = make_float3 (0.0f);
    if(_S1034)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.y) == 0.0f;
    }
    if(_runFlag_16)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.z) == 0.0f;
    }
    if(_runFlag_16)
    {
        _runFlag_16 = true;
    }
    else
    {
        _runFlag_16 = (depths_9.w) == 0.0f;
    }
    bool _S1053 = !_runFlag_16;
    bool _runFlag_17;
    bool _runFlag_18;
    bool _S1054;
    bool _runFlag_19;
    bool _S1055;
    bool _S1056;
    FixedArray<float3 , 4>  points_11;
    float3  _S1057;
    float3  _S1058;
    float3  _S1059;
    float3  _S1060;
    float3  _S1061;
    float3  _S1062;
    float3  _S1063;
    float3  _S1064;
    float3  _S1065;
    if(_S1053)
    {
        bool _S1066 = !!_S1051._S1025;
        if(_S1066)
        {
            float3  _S1067 = s_primal_ctx_unproject_raydir_0(_S1051._S1024, camera_model_22, is_ray_depth_20);
            float3  _S1068 = make_float3 (depths_9.x) * _S1067;
            bool _S1069 = !!_S1051._S1027;
            if(_S1069)
            {
                float3  _S1070 = s_primal_ctx_unproject_raydir_0(_S1051._S1026, camera_model_22, is_ray_depth_20);
                float3  _S1071 = make_float3 (depths_9.y) * _S1070;
                _S1035 = int(2);
                points_11[int(0)] = _S1068;
                points_11[int(1)] = _S1071;
                points_11[int(2)] = _S1052;
                points_11[int(3)] = _S1052;
                _S1057 = _S1070;
            }
            else
            {
                _S1035 = int(0);
                points_11[int(0)] = _S1068;
                points_11[int(1)] = _S1052;
                points_11[int(2)] = _S1052;
                points_11[int(3)] = _S1052;
                _S1057 = _S1052;
            }
            if(_S1035 != int(2))
            {
                _runFlag_16 = false;
            }
            else
            {
                _runFlag_16 = _S1066;
                _S1035 = int(0);
            }
            if(_runFlag_16)
            {
                if(!_S1051._S1029)
                {
                    _runFlag_17 = false;
                    _S1035 = int(0);
                }
                else
                {
                    _runFlag_17 = _runFlag_16;
                }
                if(_runFlag_17)
                {
                    float3  _S1072 = s_primal_ctx_unproject_raydir_0(_S1051._S1028, camera_model_22, is_ray_depth_20);
                    points_11[int(2)] = make_float3 (depths_9.z) * _S1072;
                    bool _S1073 = !!_S1051._S1031;
                    int _S1074;
                    if(_S1073)
                    {
                        float3  _S1075 = s_primal_ctx_unproject_raydir_0(_S1051._S1030, camera_model_22, is_ray_depth_20);
                        points_11[int(3)] = make_float3 (depths_9.w) * _S1075;
                        _S1074 = int(2);
                        _S1058 = _S1075;
                    }
                    else
                    {
                        _S1074 = int(0);
                        _S1058 = _S1052;
                    }
                    if(_S1074 != int(2))
                    {
                        _runFlag_18 = false;
                        _S1035 = _S1074;
                    }
                    else
                    {
                        _runFlag_18 = _runFlag_17;
                    }
                    if(_runFlag_18)
                    {
                        _S1035 = int(1);
                    }
                    _runFlag_18 = _S1073;
                    _S1059 = _S1072;
                }
                else
                {
                    _runFlag_18 = false;
                    _S1058 = _S1052;
                    _S1059 = _S1052;
                }
            }
            else
            {
                _runFlag_17 = false;
                _runFlag_18 = false;
                _S1058 = _S1052;
                _S1059 = _S1052;
            }
            float3  _S1076 = _S1057;
            _S1057 = _S1058;
            _S1058 = _S1059;
            _S1054 = _S1069;
            _S1059 = _S1076;
            _S1060 = _S1067;
        }
        else
        {
            _S1035 = int(0);
            points_11[int(0)] = _S1052;
            points_11[int(1)] = _S1052;
            points_11[int(2)] = _S1052;
            points_11[int(3)] = _S1052;
            _runFlag_16 = false;
            _runFlag_17 = false;
            _runFlag_18 = false;
            _S1057 = _S1052;
            _S1058 = _S1052;
            _S1054 = false;
            _S1059 = _S1052;
            _S1060 = _S1052;
        }
        if(_S1035 != int(1))
        {
            _runFlag_19 = false;
        }
        else
        {
            _runFlag_19 = _S1053;
        }
        if(_runFlag_19)
        {
            float3  dx_5 = points_11[int(1)] - points_11[int(0)];
            float3  _S1077 = - (points_11[int(3)] - points_11[int(2)]);
            float3  _S1078 = s_primal_ctx_cross_0(dx_5, _S1077);
            bool _S1079 = (s_primal_ctx_dot_0(_S1078, _S1078)) != 0.0f;
            if(_S1079)
            {
                float _S1080 = length_0(_S1078);
                float3  _S1081 = make_float3 (_S1080);
                _S1061 = make_float3 (_S1080 * _S1080);
                _S1062 = _S1081;
            }
            else
            {
                _S1061 = _S1052;
                _S1062 = _S1052;
            }
            float3  _S1082 = _S1062;
            _S1055 = _S1079;
            _S1062 = _S1078;
            _S1063 = _S1082;
            _S1064 = dx_5;
            _S1065 = _S1077;
        }
        else
        {
            _S1055 = false;
            _S1061 = _S1052;
            _S1062 = _S1052;
            _S1063 = _S1052;
            _S1064 = _S1052;
            _S1065 = _S1052;
        }
        bool _S1083 = _runFlag_16;
        bool _S1084 = _runFlag_17;
        bool _S1085 = _runFlag_18;
        float3  _S1086 = _S1057;
        float3  _S1087 = _S1058;
        bool _S1088 = _S1054;
        float3  _S1089 = _S1059;
        float3  _S1090 = _S1060;
        _runFlag_16 = _runFlag_19;
        _runFlag_17 = _S1055;
        _S1057 = _S1061;
        _S1058 = _S1062;
        _S1059 = _S1063;
        _S1060 = _S1064;
        _S1061 = _S1065;
        _runFlag_18 = _S1066;
        _S1054 = _S1083;
        _runFlag_19 = _S1084;
        _S1055 = _S1085;
        _S1062 = _S1086;
        _S1063 = _S1087;
        _S1056 = _S1088;
        _S1064 = _S1089;
        _S1065 = _S1090;
    }
    else
    {
        _runFlag_16 = false;
        _runFlag_17 = false;
        _S1057 = _S1052;
        _S1058 = _S1052;
        _S1059 = _S1052;
        _S1060 = _S1052;
        _S1061 = _S1052;
        _runFlag_18 = false;
        _S1054 = false;
        _runFlag_19 = false;
        _S1055 = false;
        _S1062 = _S1052;
        _S1063 = _S1052;
        _S1056 = false;
        _S1064 = _S1052;
        _S1065 = _S1052;
    }
    float4  _S1091 = make_float4 (0.0f);
    float4  _S1092;
    if(_S1053)
    {
        if(_runFlag_16)
        {
            if(_runFlag_17)
            {
                float3  _S1093 = v_normal_3 / _S1057;
                float3  _S1094 = _S1058 * - _S1093;
                float3  _S1095 = _S1059 * _S1093;
                float _S1096 = _S1094.x + _S1094.y + _S1094.z;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1097;
                (&_S1097)->primal_0 = _S1058;
                (&_S1097)->differential_0 = _S1052;
                s_bwd_length_impl_0(&_S1097, _S1096);
                _S1057 = _S1095 + _S1097.differential_0;
            }
            else
            {
                _S1057 = v_normal_3;
            }
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1098;
            (&_S1098)->primal_0 = _S1058;
            (&_S1098)->differential_0 = _S1052;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1099;
            (&_S1099)->primal_0 = _S1058;
            (&_S1099)->differential_0 = _S1052;
            s_bwd_prop_dot_0(&_S1098, &_S1099, 0.0f);
            float3  _S1100 = _S1099.differential_0 + _S1098.differential_0 + _S1057;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1101;
            (&_S1101)->primal_0 = _S1060;
            (&_S1101)->differential_0 = _S1052;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1102;
            (&_S1102)->primal_0 = _S1061;
            (&_S1102)->differential_0 = _S1052;
            s_bwd_prop_cross_0(&_S1101, &_S1102, _S1100);
            float3  s_diff_dy_T_5 = - _S1102.differential_0;
            float3  _S1103 = - s_diff_dy_T_5;
            float3  _S1104 = - _S1101.differential_0;
            FixedArray<float3 , 4>  _S1105;
            _S1105[int(0)] = _S1052;
            _S1105[int(1)] = _S1052;
            _S1105[int(2)] = _S1052;
            _S1105[int(3)] = _S1052;
            _S1105[int(2)] = _S1103;
            _S1105[int(3)] = s_diff_dy_T_5;
            _S1105[int(0)] = _S1104;
            _S1105[int(1)] = _S1101.differential_0;
            points_11[int(0)] = _S1105[int(0)];
            points_11[int(1)] = _S1105[int(1)];
            points_11[int(2)] = _S1105[int(2)];
            points_11[int(3)] = _S1105[int(3)];
        }
        else
        {
            points_11[int(0)] = _S1052;
            points_11[int(1)] = _S1052;
            points_11[int(2)] = _S1052;
            points_11[int(3)] = _S1052;
        }
        if(_runFlag_18)
        {
            if(_S1054)
            {
                if(_runFlag_19)
                {
                    FixedArray<float3 , 4>  _S1106 = points_11;
                    FixedArray<float3 , 4>  _S1107 = points_11;
                    FixedArray<float3 , 4>  _S1108 = points_11;
                    FixedArray<float3 , 4>  _S1109 = points_11;
                    if(_S1055)
                    {
                        float3  _S1110 = _S1062 * _S1109[int(3)];
                        float _S1111 = _S1110.x + _S1110.y + _S1110.z;
                        float4  _S1112 = _S1091;
                        *&((&_S1112)->w) = _S1111;
                        points_11[int(0)] = _S1106[int(0)];
                        points_11[int(1)] = _S1107[int(1)];
                        points_11[int(2)] = _S1108[int(2)];
                        points_11[int(3)] = _S1052;
                        _S1092 = _S1112;
                    }
                    else
                    {
                        points_11[int(0)] = _S1106[int(0)];
                        points_11[int(1)] = _S1107[int(1)];
                        points_11[int(2)] = _S1108[int(2)];
                        points_11[int(3)] = _S1109[int(3)];
                        _S1092 = _S1091;
                    }
                    float3  _S1113 = _S1063 * points_11[int(2)];
                    float _S1114 = _S1113.x + _S1113.y + _S1113.z;
                    FixedArray<float3 , 4>  _S1115 = points_11;
                    FixedArray<float3 , 4>  _S1116 = points_11;
                    float4  _S1117 = _S1091;
                    *&((&_S1117)->z) = _S1114;
                    float4  _S1118 = _S1092 + _S1117;
                    points_11[int(0)] = points_11[int(0)];
                    points_11[int(1)] = _S1115[int(1)];
                    points_11[int(2)] = _S1052;
                    points_11[int(3)] = _S1116[int(3)];
                    _S1092 = _S1118;
                }
                else
                {
                    FixedArray<float3 , 4>  _S1119 = points_11;
                    FixedArray<float3 , 4>  _S1120 = points_11;
                    FixedArray<float3 , 4>  _S1121 = points_11;
                    points_11[int(0)] = points_11[int(0)];
                    points_11[int(1)] = _S1119[int(1)];
                    points_11[int(2)] = _S1120[int(2)];
                    points_11[int(3)] = _S1121[int(3)];
                    _S1092 = _S1091;
                }
            }
            else
            {
                FixedArray<float3 , 4>  _S1122 = points_11;
                FixedArray<float3 , 4>  _S1123 = points_11;
                FixedArray<float3 , 4>  _S1124 = points_11;
                points_11[int(0)] = points_11[int(0)];
                points_11[int(1)] = _S1122[int(1)];
                points_11[int(2)] = _S1123[int(2)];
                points_11[int(3)] = _S1124[int(3)];
                _S1092 = _S1091;
            }
            if(_S1056)
            {
                FixedArray<float3 , 4>  _S1125 = points_11;
                float3  _S1126 = _S1064 * points_11[int(1)];
                float _S1127 = _S1126.x + _S1126.y + _S1126.z;
                float4  _S1128 = _S1091;
                *&((&_S1128)->y) = _S1127;
                float4  _S1129 = _S1092 + _S1128;
                points_11[int(0)] = _S1052;
                points_11[int(1)] = _S1052;
                points_11[int(2)] = _S1052;
                points_11[int(3)] = _S1052;
                _S1057 = _S1125[int(0)];
                _S1092 = _S1129;
            }
            else
            {
                FixedArray<float3 , 4>  _S1130 = points_11;
                FixedArray<float3 , 4>  _S1131 = points_11;
                FixedArray<float3 , 4>  _S1132 = points_11;
                points_11[int(0)] = points_11[int(0)];
                points_11[int(1)] = _S1130[int(1)];
                points_11[int(2)] = _S1131[int(2)];
                points_11[int(3)] = _S1132[int(3)];
                _S1057 = _S1052;
            }
            float3  _S1133 = _S1065 * (points_11[int(0)] + _S1057);
            float _S1134 = _S1133.x + _S1133.y + _S1133.z;
            float4  _S1135 = _S1091;
            *&((&_S1135)->x) = _S1134;
            _S1092 = _S1092 + _S1135;
        }
        else
        {
            _S1092 = _S1091;
        }
    }
    else
    {
        _S1092 = _S1091;
    }
    *v_depths_4 = _S1092;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor_prism(float2  pix_center_12, float4  intrins_21, FixedArray<float, 8>  dist_coeffs_25, int camera_model_23)
{
    float _S1136;
    for(;;)
    {
        float2  uv_63 = (pix_center_12 - float2 {intrins_21.z, intrins_21.w}) / float2 {intrins_21.x, intrins_21.y};
        FixedArray<float, 8>  _S1137 = dist_coeffs_25;
        float2  uv_u_30;
        bool _S1138 = undistort_point_2(uv_63, &_S1137, int(12), &uv_u_30);
        if(!_S1138)
        {
            _S1136 = 0.0f;
            break;
        }
        float3  raydir_16 = unproject_raydir_0(uv_u_30, camera_model_23, false);
        _S1136 = float((F32_sign((raydir_16.z)))) / length_0(raydir_16);
        break;
    }
    return _S1136;
}

inline __device__ float depth_normal_loss_prism(float2  pix_center_13, float4  intrins_22, FixedArray<float, 8>  dist_coeffs_26, int camera_model_24, bool is_ray_depth_21, float4  depths_10, float3  gt_normal_4)
{
    float _S1139;
    for(;;)
    {
        float3  _S1140;
        float3  * _S1141;
        float3  * _S1142;
        float3  * _S1143;
        float3  * _S1144;
        int _S1145;
        FixedArray<float3 , 5>  points_12;
        for(;;)
        {
            float2  _S1146 = float2 {intrins_22.z, intrins_22.w};
            float2  _S1147 = float2 {intrins_22.x, intrins_22.y};
            float2  uv_64 = (pix_center_13 + make_float2 (-1.0f, -0.0f) - _S1146) / _S1147;
            FixedArray<float, 8>  _S1148 = dist_coeffs_26;
            float2  uv_u_31;
            bool _S1149 = undistort_point_2(uv_64, &_S1148, int(12), &uv_u_31);
            float3  _S1150 = make_float3 (0.0f);
            if(!_S1149)
            {
                _S1145 = int(0);
                _S1144 = nullptr;
                _S1143 = nullptr;
                _S1142 = nullptr;
                _S1141 = nullptr;
                _S1140 = _S1150;
                break;
            }
            float3  raydir_17 = unproject_raydir_0(uv_u_31, camera_model_24, is_ray_depth_21);
            points_12[int(0)] = make_float3 (depths_10.x) * raydir_17;
            float2  uv_65 = (pix_center_13 + make_float2 (1.0f, -0.0f) - _S1146) / _S1147;
            FixedArray<float, 8>  _S1151 = dist_coeffs_26;
            float2  uv_u_32;
            bool _S1152 = undistort_point_2(uv_65, &_S1151, int(12), &uv_u_32);
            if(!_S1152)
            {
                _S1145 = int(0);
                _S1144 = nullptr;
                _S1143 = &points_12[int(0)];
                _S1142 = nullptr;
                _S1141 = nullptr;
                _S1140 = _S1150;
                break;
            }
            float3  raydir_18 = unproject_raydir_0(uv_u_32, camera_model_24, is_ray_depth_21);
            points_12[int(1)] = make_float3 (depths_10.y) * raydir_18;
            float2  uv_66 = (pix_center_13 + make_float2 (0.0f, -1.0f) - _S1146) / _S1147;
            FixedArray<float, 8>  _S1153 = dist_coeffs_26;
            float2  uv_u_33;
            bool _S1154 = undistort_point_2(uv_66, &_S1153, int(12), &uv_u_33);
            if(!_S1154)
            {
                _S1145 = int(0);
                _S1144 = &points_12[int(1)];
                _S1143 = &points_12[int(0)];
                _S1142 = nullptr;
                _S1141 = nullptr;
                _S1140 = _S1150;
                break;
            }
            float3  raydir_19 = unproject_raydir_0(uv_u_33, camera_model_24, is_ray_depth_21);
            points_12[int(2)] = make_float3 (depths_10.z) * raydir_19;
            float2  uv_67 = (pix_center_13 + make_float2 (0.0f, 1.0f) - _S1146) / _S1147;
            FixedArray<float, 8>  _S1155 = dist_coeffs_26;
            float2  uv_u_34;
            bool _S1156 = undistort_point_2(uv_67, &_S1155, int(12), &uv_u_34);
            if(!_S1156)
            {
                _S1145 = int(0);
                _S1144 = &points_12[int(1)];
                _S1143 = &points_12[int(0)];
                _S1142 = nullptr;
                _S1141 = &points_12[int(2)];
                _S1140 = _S1150;
                break;
            }
            float3  raydir_20 = unproject_raydir_0(uv_u_34, camera_model_24, is_ray_depth_21);
            points_12[int(3)] = make_float3 (depths_10.w) * raydir_20;
            float2  uv_68 = (pix_center_13 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S1146) / _S1147;
            FixedArray<float, 8>  _S1157 = dist_coeffs_26;
            float2  uv_u_35;
            bool _S1158 = undistort_point_2(uv_68, &_S1157, int(12), &uv_u_35);
            if(!_S1158)
            {
                _S1145 = int(0);
                _S1144 = &points_12[int(1)];
                _S1143 = &points_12[int(0)];
                _S1142 = &points_12[int(3)];
                _S1141 = &points_12[int(2)];
                _S1140 = _S1150;
                break;
            }
            float3  raydir_21 = unproject_raydir_0(uv_u_35, camera_model_24, is_ray_depth_21);
            _S1145 = int(1);
            _S1144 = &points_12[int(1)];
            _S1143 = &points_12[int(0)];
            _S1142 = &points_12[int(3)];
            _S1141 = &points_12[int(2)];
            _S1140 = raydir_21;
            break;
        }
        if(_S1145 != int(1))
        {
            _S1139 = 0.0f;
            break;
        }
        float3  normal_14 = cross_0(*_S1144 - *_S1143, - (*_S1142 - *_S1141));
        float3  normal_15;
        if((dot_0(normal_14, normal_14)) != 0.0f)
        {
            normal_15 = normalize_0(normal_14);
        }
        else
        {
            normal_15 = normal_14;
        }
        float3  _S1159;
        if((dot_0(gt_normal_4, gt_normal_4)) != 0.0f)
        {
            _S1159 = normalize_0(gt_normal_4);
        }
        else
        {
            _S1159 = gt_normal_4;
        }
        _S1139 = (1.0f - dot_0(normal_15, _S1159) + 0.00100000004749745f) / ((F32_max((dot_0(normal_15, - normalize_0(_S1140))), (0.0f))) + 0.00100000004749745f);
        break;
    }
    return _S1139;
}

struct s_bwd_prop_depth_normal_loss_Intermediates_2
{
    float2  _S1160;
    bool _S1161;
    float2  _S1162;
    bool _S1163;
    float2  _S1164;
    bool _S1165;
    float2  _S1166;
    bool _S1167;
    float2  _S1168;
    bool _S1169;
};

inline __device__ void depth_normal_loss_vjp_prism(float2  pix_center_14, float4  intrins_23, FixedArray<float, 8>  dist_coeffs_27, int camera_model_25, bool is_ray_depth_22, float4  depths_11, float3  gt_normal_5, float v_loss_2, float4  * v_depths_5, float3  * v_gt_normal_2)
{
    float2  _S1170 = make_float2 (0.0f);
    s_bwd_prop_depth_normal_loss_Intermediates_2 _S1171;
    (&_S1171)->_S1160 = _S1170;
    (&_S1171)->_S1161 = false;
    (&_S1171)->_S1162 = _S1170;
    (&_S1171)->_S1163 = false;
    (&_S1171)->_S1164 = _S1170;
    (&_S1171)->_S1165 = false;
    (&_S1171)->_S1166 = _S1170;
    (&_S1171)->_S1167 = false;
    (&_S1171)->_S1168 = _S1170;
    (&_S1171)->_S1169 = false;
    (&_S1171)->_S1162 = _S1170;
    (&_S1171)->_S1163 = false;
    (&_S1171)->_S1164 = _S1170;
    (&_S1171)->_S1165 = false;
    (&_S1171)->_S1166 = _S1170;
    (&_S1171)->_S1167 = false;
    (&_S1171)->_S1168 = _S1170;
    (&_S1171)->_S1169 = false;
    float2  _S1172 = float2 {intrins_23.z, intrins_23.w};
    float2  _S1173 = float2 {intrins_23.x, intrins_23.y};
    float2  uv_69 = (pix_center_14 + make_float2 (-1.0f, -0.0f) - _S1172) / _S1173;
    float2  _S1174 = _S1170;
    FixedArray<float, 8>  _S1175 = dist_coeffs_27;
    bool _S1176 = undistort_point_2(uv_69, &_S1175, int(12), &_S1174);
    (&_S1171)->_S1160 = _S1174;
    (&_S1171)->_S1161 = _S1176;
    bool _S1177 = !!_S1176;
    bool _runFlag_20;
    if(_S1177)
    {
        float2  uv_70 = (pix_center_14 + make_float2 (1.0f, -0.0f) - _S1172) / _S1173;
        float2  _S1178 = _S1170;
        FixedArray<float, 8>  _S1179 = dist_coeffs_27;
        bool _S1180 = undistort_point_2(uv_70, &_S1179, int(12), &_S1178);
        (&_S1171)->_S1162 = _S1178;
        (&_S1171)->_S1163 = _S1180;
        if(!_S1180)
        {
            _runFlag_20 = false;
        }
        else
        {
            _runFlag_20 = _S1177;
        }
        if(_runFlag_20)
        {
            float2  uv_71 = (pix_center_14 + make_float2 (0.0f, -1.0f) - _S1172) / _S1173;
            float2  _S1181 = _S1170;
            FixedArray<float, 8>  _S1182 = dist_coeffs_27;
            bool _S1183 = undistort_point_2(uv_71, &_S1182, int(12), &_S1181);
            (&_S1171)->_S1164 = _S1181;
            (&_S1171)->_S1165 = _S1183;
            if(!_S1183)
            {
                _runFlag_20 = false;
            }
            if(_runFlag_20)
            {
                float2  uv_72 = (pix_center_14 + make_float2 (0.0f, 1.0f) - _S1172) / _S1173;
                float2  _S1184 = _S1170;
                FixedArray<float, 8>  _S1185 = dist_coeffs_27;
                bool _S1186 = undistort_point_2(uv_72, &_S1185, int(12), &_S1184);
                (&_S1171)->_S1166 = _S1184;
                (&_S1171)->_S1167 = _S1186;
                if(!_S1186)
                {
                    _runFlag_20 = false;
                }
                if(_runFlag_20)
                {
                    float2  uv_73 = (pix_center_14 - _S1172) / _S1173;
                    float2  _S1187 = _S1170;
                    FixedArray<float, 8>  _S1188 = dist_coeffs_27;
                    bool _S1189 = undistort_point_2(uv_73, &_S1188, int(12), &_S1187);
                    (&_S1171)->_S1168 = _S1187;
                    (&_S1171)->_S1169 = _S1189;
                }
            }
        }
    }
    s_bwd_prop_depth_normal_loss_Intermediates_2 _S1190 = _S1171;
    float3  _S1191 = make_float3 (0.0f);
    bool _S1192 = !!_S1171._S1161;
    bool _runFlag_21;
    bool _runFlag_22;
    bool _runFlag_23;
    int _S1193;
    float3  raydir_22;
    float3  _S1194;
    float3  _S1195;
    float3  _S1196;
    float3  _S1197;
    FixedArray<float3 , 5>  points_13;
    if(_S1192)
    {
        float3  _S1198 = s_primal_ctx_unproject_raydir_0(_S1190._S1160, camera_model_25, is_ray_depth_22);
        float3  _S1199 = make_float3 (depths_11.x) * _S1198;
        if(!_S1190._S1163)
        {
            _runFlag_20 = false;
        }
        else
        {
            _runFlag_20 = _S1192;
        }
        if(_runFlag_20)
        {
            float3  _S1200 = s_primal_ctx_unproject_raydir_0(_S1190._S1162, camera_model_25, is_ray_depth_22);
            float3  _S1201 = make_float3 (depths_11.y) * _S1200;
            if(!_S1190._S1165)
            {
                _runFlag_21 = false;
            }
            else
            {
                _runFlag_21 = _runFlag_20;
            }
            if(_runFlag_21)
            {
                float3  _S1202 = s_primal_ctx_unproject_raydir_0(_S1190._S1164, camera_model_25, is_ray_depth_22);
                float3  _S1203 = make_float3 (depths_11.z) * _S1202;
                if(!_S1190._S1167)
                {
                    _runFlag_22 = false;
                }
                else
                {
                    _runFlag_22 = _runFlag_21;
                }
                if(_runFlag_22)
                {
                    float3  _S1204 = s_primal_ctx_unproject_raydir_0(_S1190._S1166, camera_model_25, is_ray_depth_22);
                    float3  _S1205 = make_float3 (depths_11.w) * _S1204;
                    if(!_S1190._S1169)
                    {
                        _runFlag_23 = false;
                    }
                    else
                    {
                        _runFlag_23 = _runFlag_22;
                    }
                    if(_runFlag_23)
                    {
                        float3  _S1206 = s_primal_ctx_unproject_raydir_0(_S1190._S1168, camera_model_25, is_ray_depth_22);
                        _S1193 = int(1);
                        raydir_22 = _S1206;
                    }
                    else
                    {
                        _S1193 = int(0);
                        raydir_22 = _S1204;
                    }
                    points_13[int(0)] = _S1199;
                    points_13[int(1)] = _S1201;
                    points_13[int(2)] = _S1203;
                    points_13[int(3)] = _S1205;
                    points_13[int(4)] = _S1191;
                    _S1194 = _S1204;
                }
                else
                {
                    _S1193 = int(0);
                    raydir_22 = _S1202;
                    points_13[int(0)] = _S1199;
                    points_13[int(1)] = _S1201;
                    points_13[int(2)] = _S1203;
                    points_13[int(3)] = _S1191;
                    points_13[int(4)] = _S1191;
                    _S1194 = _S1191;
                }
                _S1195 = _S1202;
            }
            else
            {
                _S1193 = int(0);
                raydir_22 = _S1200;
                points_13[int(0)] = _S1199;
                points_13[int(1)] = _S1201;
                points_13[int(2)] = _S1191;
                points_13[int(3)] = _S1191;
                points_13[int(4)] = _S1191;
                _runFlag_22 = false;
                _S1194 = _S1191;
                _S1195 = _S1191;
            }
            _S1196 = _S1200;
        }
        else
        {
            _S1193 = int(0);
            raydir_22 = _S1198;
            points_13[int(0)] = _S1199;
            points_13[int(1)] = _S1191;
            points_13[int(2)] = _S1191;
            points_13[int(3)] = _S1191;
            points_13[int(4)] = _S1191;
            _runFlag_21 = false;
            _runFlag_22 = false;
            _S1194 = _S1191;
            _S1195 = _S1191;
            _S1196 = _S1191;
        }
        _S1197 = _S1198;
    }
    else
    {
        _S1193 = int(0);
        points_13[int(0)] = _S1191;
        points_13[int(1)] = _S1191;
        points_13[int(2)] = _S1191;
        points_13[int(3)] = _S1191;
        points_13[int(4)] = _S1191;
        _runFlag_20 = false;
        _runFlag_21 = false;
        _runFlag_22 = false;
        _S1194 = _S1191;
        _S1195 = _S1191;
        _S1196 = _S1191;
        _S1197 = _S1191;
    }
    bool _S1207 = !(_S1193 != int(1));
    bool _S1208;
    float3  normal_16;
    float3  _S1209;
    float3  _S1210;
    float3  _S1211;
    float3  _S1212;
    float _S1213;
    float _S1214;
    float _S1215;
    float _S1216;
    if(_S1207)
    {
        float3  dx_6 = points_13[int(1)] - points_13[int(0)];
        float3  _S1217 = - (points_13[int(3)] - points_13[int(2)]);
        float3  _S1218 = s_primal_ctx_cross_0(dx_6, _S1217);
        bool _S1219 = (s_primal_ctx_dot_0(_S1218, _S1218)) != 0.0f;
        if(_S1219)
        {
            normal_16 = normalize_0(_S1218);
        }
        else
        {
            normal_16 = _S1218;
        }
        bool _S1220 = (s_primal_ctx_dot_0(gt_normal_5, gt_normal_5)) != 0.0f;
        if(_S1220)
        {
            _S1209 = normalize_0(gt_normal_5);
        }
        else
        {
            _S1209 = gt_normal_5;
        }
        float3  _S1221 = - normalize_0(raydir_22);
        float _S1222 = s_primal_ctx_dot_0(normal_16, _S1221);
        float _S1223 = 1.0f - s_primal_ctx_dot_0(normal_16, _S1209) + 0.00100000004749745f;
        float _S1224 = (F32_max((_S1222), (0.0f))) + 0.00100000004749745f;
        _S1213 = _S1224 * _S1224;
        _S1214 = _S1223;
        _S1215 = _S1224;
        _S1216 = _S1222;
        raydir_22 = normal_16;
        normal_16 = _S1221;
        _runFlag_23 = _S1220;
        _S1208 = _S1219;
        _S1210 = _S1218;
        _S1211 = dx_6;
        _S1212 = _S1217;
    }
    else
    {
        _S1213 = 0.0f;
        _S1214 = 0.0f;
        _S1215 = 0.0f;
        _S1216 = 0.0f;
        raydir_22 = _S1191;
        normal_16 = _S1191;
        _S1209 = _S1191;
        _runFlag_23 = false;
        _S1208 = false;
        _S1210 = _S1191;
        _S1211 = _S1191;
        _S1212 = _S1191;
    }
    float4  _S1225 = make_float4 (0.0f);
    if(_S1207)
    {
        float _S1226 = v_loss_2 / _S1213;
        float _S1227 = _S1214 * - _S1226;
        float s_diff_num_T_2 = _S1215 * _S1226;
        DiffPair_float_0 _S1228;
        (&_S1228)->primal_0 = _S1216;
        (&_S1228)->differential_0 = 0.0f;
        DiffPair_float_0 _S1229;
        (&_S1229)->primal_0 = 0.0f;
        (&_S1229)->differential_0 = 0.0f;
        _d_max_0(&_S1228, &_S1229, _S1227);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1230;
        (&_S1230)->primal_0 = raydir_22;
        (&_S1230)->differential_0 = _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1231;
        (&_S1231)->primal_0 = normal_16;
        (&_S1231)->differential_0 = _S1191;
        s_bwd_prop_dot_0(&_S1230, &_S1231, _S1228.differential_0);
        float _S1232 = - s_diff_num_T_2;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1233;
        (&_S1233)->primal_0 = raydir_22;
        (&_S1233)->differential_0 = _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1234;
        (&_S1234)->primal_0 = _S1209;
        (&_S1234)->differential_0 = _S1191;
        s_bwd_prop_dot_0(&_S1233, &_S1234, _S1232);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1235 = _S1234;
        float3  _S1236 = _S1230.differential_0 + _S1233.differential_0;
        if(_runFlag_23)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1237;
            (&_S1237)->primal_0 = gt_normal_5;
            (&_S1237)->differential_0 = _S1191;
            s_bwd_normalize_impl_0(&_S1237, _S1235.differential_0);
            raydir_22 = _S1237.differential_0;
        }
        else
        {
            raydir_22 = _S1235.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1238;
        (&_S1238)->primal_0 = gt_normal_5;
        (&_S1238)->differential_0 = _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239;
        (&_S1239)->primal_0 = gt_normal_5;
        (&_S1239)->differential_0 = _S1191;
        s_bwd_prop_dot_0(&_S1238, &_S1239, 0.0f);
        float3  _S1240 = _S1239.differential_0 + _S1238.differential_0 + raydir_22;
        if(_S1208)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1241;
            (&_S1241)->primal_0 = _S1210;
            (&_S1241)->differential_0 = _S1191;
            s_bwd_normalize_impl_0(&_S1241, _S1236);
            raydir_22 = _S1241.differential_0;
        }
        else
        {
            raydir_22 = _S1236;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1242;
        (&_S1242)->primal_0 = _S1210;
        (&_S1242)->differential_0 = _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1243;
        (&_S1243)->primal_0 = _S1210;
        (&_S1243)->differential_0 = _S1191;
        s_bwd_prop_dot_0(&_S1242, &_S1243, 0.0f);
        float3  _S1244 = _S1243.differential_0 + _S1242.differential_0 + raydir_22;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1245;
        (&_S1245)->primal_0 = _S1211;
        (&_S1245)->differential_0 = _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1246;
        (&_S1246)->primal_0 = _S1212;
        (&_S1246)->differential_0 = _S1191;
        s_bwd_prop_cross_0(&_S1245, &_S1246, _S1244);
        float3  s_diff_dy_T_6 = - _S1246.differential_0;
        float3  _S1247 = - s_diff_dy_T_6;
        float3  _S1248 = - _S1245.differential_0;
        FixedArray<float3 , 5>  _S1249;
        _S1249[int(0)] = _S1191;
        _S1249[int(1)] = _S1191;
        _S1249[int(2)] = _S1191;
        _S1249[int(3)] = _S1191;
        _S1249[int(4)] = _S1191;
        _S1249[int(2)] = _S1247;
        _S1249[int(3)] = s_diff_dy_T_6;
        _S1249[int(0)] = _S1248;
        _S1249[int(1)] = _S1245.differential_0;
        points_13[int(0)] = _S1249[int(0)];
        points_13[int(1)] = _S1249[int(1)];
        points_13[int(2)] = _S1249[int(2)];
        points_13[int(3)] = _S1249[int(3)];
        points_13[int(4)] = _S1249[int(4)];
        raydir_22 = _S1240;
    }
    else
    {
        points_13[int(0)] = _S1191;
        points_13[int(1)] = _S1191;
        points_13[int(2)] = _S1191;
        points_13[int(3)] = _S1191;
        points_13[int(4)] = _S1191;
        raydir_22 = _S1191;
    }
    float4  _S1250;
    if(_S1192)
    {
        if(_runFlag_20)
        {
            if(_runFlag_21)
            {
                if(_runFlag_22)
                {
                    FixedArray<float3 , 5>  _S1251 = points_13;
                    FixedArray<float3 , 5>  _S1252 = points_13;
                    FixedArray<float3 , 5>  _S1253 = points_13;
                    float3  _S1254 = _S1194 * points_13[int(3)];
                    float _S1255 = _S1254.x + _S1254.y + _S1254.z;
                    float4  _S1256 = _S1225;
                    *&((&_S1256)->w) = _S1255;
                    points_13[int(0)] = _S1191;
                    points_13[int(1)] = _S1191;
                    points_13[int(2)] = _S1191;
                    points_13[int(3)] = _S1191;
                    points_13[int(4)] = _S1191;
                    _S1194 = _S1253[int(2)];
                    normal_16 = _S1251[int(0)];
                    _S1209 = _S1252[int(1)];
                    _S1250 = _S1256;
                }
                else
                {
                    FixedArray<float3 , 5>  _S1257 = points_13;
                    FixedArray<float3 , 5>  _S1258 = points_13;
                    FixedArray<float3 , 5>  _S1259 = points_13;
                    FixedArray<float3 , 5>  _S1260 = points_13;
                    points_13[int(0)] = points_13[int(0)];
                    points_13[int(1)] = _S1257[int(1)];
                    points_13[int(2)] = _S1258[int(2)];
                    points_13[int(3)] = _S1259[int(3)];
                    points_13[int(4)] = _S1260[int(4)];
                    _S1194 = _S1191;
                    normal_16 = _S1191;
                    _S1209 = _S1191;
                    _S1250 = _S1225;
                }
                float3  _S1261 = _S1195 * (points_13[int(2)] + _S1194);
                float _S1262 = _S1261.x + _S1261.y + _S1261.z;
                float3  _S1263 = points_13[int(0)] + normal_16;
                float3  _S1264 = points_13[int(1)] + _S1209;
                float4  _S1265 = _S1225;
                *&((&_S1265)->z) = _S1262;
                float4  _S1266 = _S1250 + _S1265;
                points_13[int(0)] = _S1191;
                points_13[int(1)] = _S1191;
                points_13[int(2)] = _S1191;
                points_13[int(3)] = _S1191;
                points_13[int(4)] = _S1191;
                _S1194 = _S1264;
                _S1195 = _S1263;
                _S1250 = _S1266;
            }
            else
            {
                FixedArray<float3 , 5>  _S1267 = points_13;
                FixedArray<float3 , 5>  _S1268 = points_13;
                FixedArray<float3 , 5>  _S1269 = points_13;
                FixedArray<float3 , 5>  _S1270 = points_13;
                points_13[int(0)] = points_13[int(0)];
                points_13[int(1)] = _S1267[int(1)];
                points_13[int(2)] = _S1268[int(2)];
                points_13[int(3)] = _S1269[int(3)];
                points_13[int(4)] = _S1270[int(4)];
                _S1194 = _S1191;
                _S1195 = _S1191;
                _S1250 = _S1225;
            }
            float3  _S1271 = _S1196 * (points_13[int(1)] + _S1194);
            float _S1272 = _S1271.x + _S1271.y + _S1271.z;
            float3  _S1273 = points_13[int(0)] + _S1195;
            float4  _S1274 = _S1225;
            *&((&_S1274)->y) = _S1272;
            float4  _S1275 = _S1250 + _S1274;
            points_13[int(0)] = _S1191;
            points_13[int(1)] = _S1191;
            points_13[int(2)] = _S1191;
            points_13[int(3)] = _S1191;
            points_13[int(4)] = _S1191;
            _S1194 = _S1273;
            _S1250 = _S1275;
        }
        else
        {
            FixedArray<float3 , 5>  _S1276 = points_13;
            FixedArray<float3 , 5>  _S1277 = points_13;
            FixedArray<float3 , 5>  _S1278 = points_13;
            FixedArray<float3 , 5>  _S1279 = points_13;
            points_13[int(0)] = points_13[int(0)];
            points_13[int(1)] = _S1276[int(1)];
            points_13[int(2)] = _S1277[int(2)];
            points_13[int(3)] = _S1278[int(3)];
            points_13[int(4)] = _S1279[int(4)];
            _S1194 = _S1191;
            _S1250 = _S1225;
        }
        float3  _S1280 = _S1197 * (points_13[int(0)] + _S1194);
        float _S1281 = _S1280.x + _S1280.y + _S1280.z;
        float4  _S1282 = _S1225;
        *&((&_S1282)->x) = _S1281;
        _S1250 = _S1250 + _S1282;
    }
    else
    {
        _S1250 = _S1225;
    }
    *v_depths_5 = _S1250;
    *v_gt_normal_2 = raydir_22;
    return;
}

inline __device__ float3  generate_ray_d2n_rational(float2  pix_pos_9, float4  intrins_24, FixedArray<float, 8>  dist_coeffs_28, int camera_model_26, bool is_ray_depth_23)
{
    float3  _S1283;
    for(;;)
    {
        float2  uv_74 = (pix_pos_9 - float2 {intrins_24.z, intrins_24.w}) / float2 {intrins_24.x, intrins_24.y};
        FixedArray<float, 8>  _S1284 = dist_coeffs_28;
        float2  uv_u_36;
        bool _S1285 = undistort_point_3(uv_74, &_S1284, int(12), &uv_u_36);
        if(!_S1285)
        {
            int3  _S1286 = make_int3 (int(0));
            float3  _S1287 = make_float3 ((float)_S1286.x, (float)_S1286.y, (float)_S1286.z);
            _S1283 = _S1287;
            break;
        }
        _S1283 = unproject_raydir_0(uv_u_36, camera_model_26, is_ray_depth_23);
        break;
    }
    return _S1283;
}

inline __device__ float3  depth_to_point_rational(float2  pix_pos_10, float4  intrins_25, FixedArray<float, 8>  dist_coeffs_29, int camera_model_27, bool is_ray_depth_24, float depth_8)
{
    float3  _S1288;
    for(;;)
    {
        float2  uv_75 = (pix_pos_10 - float2 {intrins_25.z, intrins_25.w}) / float2 {intrins_25.x, intrins_25.y};
        FixedArray<float, 8>  _S1289 = dist_coeffs_29;
        float2  uv_u_37;
        bool _S1290 = undistort_point_3(uv_75, &_S1289, int(12), &uv_u_37);
        if(!_S1290)
        {
            _S1288 = make_float3 (0.0f);
            break;
        }
        _S1288 = make_float3 (depth_8) * unproject_raydir_0(uv_u_37, camera_model_27, is_ray_depth_24);
        break;
    }
    return _S1288;
}

struct s_bwd_prop_depth_to_point_Intermediates_3
{
    float2  _S1291;
    bool _S1292;
};

inline __device__ float depth_to_point_vjp_rational(float2  pix_pos_11, float4  intrins_26, FixedArray<float, 8>  dist_coeffs_30, int camera_model_28, bool is_ray_depth_25, float depth_9, float3  v_point_3)
{
    float2  _S1293 = make_float2 (0.0f);
    s_bwd_prop_depth_to_point_Intermediates_3 _S1294;
    (&_S1294)->_S1291 = _S1293;
    (&_S1294)->_S1292 = false;
    float2  uv_76 = (pix_pos_11 - float2 {intrins_26.z, intrins_26.w}) / float2 {intrins_26.x, intrins_26.y};
    float2  _S1295 = _S1293;
    FixedArray<float, 8>  _S1296 = dist_coeffs_30;
    bool _S1297 = undistort_point_3(uv_76, &_S1296, int(12), &_S1295);
    (&_S1294)->_S1291 = _S1295;
    (&_S1294)->_S1292 = _S1297;
    s_bwd_prop_depth_to_point_Intermediates_3 _S1298 = _S1294;
    float3  _S1299 = make_float3 (0.0f);
    bool _S1300 = !!_S1294._S1292;
    float3  _S1301;
    if(_S1300)
    {
        _S1301 = s_primal_ctx_unproject_raydir_0(_S1298._S1291, camera_model_28, is_ray_depth_25);
    }
    else
    {
        _S1301 = _S1299;
    }
    if(_S1300)
    {
        _S1301 = _S1301 * v_point_3;
    }
    else
    {
        _S1301 = _S1299;
    }
    return _S1301.x + _S1301.y + _S1301.z;
}

inline __device__ float3  depth_to_normal_rational(float2  pix_center_15, float4  intrins_27, FixedArray<float, 8>  dist_coeffs_31, int camera_model_29, bool is_ray_depth_26, float4  depths_12)
{
    float3  normal_17;
    for(;;)
    {
        bool _S1302;
        if((depths_12.x) == 0.0f)
        {
            _S1302 = true;
        }
        else
        {
            _S1302 = (depths_12.y) == 0.0f;
        }
        if(_S1302)
        {
            _S1302 = true;
        }
        else
        {
            _S1302 = (depths_12.z) == 0.0f;
        }
        if(_S1302)
        {
            _S1302 = true;
        }
        else
        {
            _S1302 = (depths_12.w) == 0.0f;
        }
        if(_S1302)
        {
            normal_17 = make_float3 (0.0f);
            break;
        }
        float3  * _S1303;
        float3  * _S1304;
        float3  * _S1305;
        float3  * _S1306;
        int _S1307;
        FixedArray<float3 , 4>  points_14;
        for(;;)
        {
            float2  _S1308 = float2 {intrins_27.z, intrins_27.w};
            float2  _S1309 = float2 {intrins_27.x, intrins_27.y};
            float2  uv_77 = (pix_center_15 + make_float2 (-1.0f, -0.0f) - _S1308) / _S1309;
            FixedArray<float, 8>  _S1310 = dist_coeffs_31;
            float2  uv_u_38;
            bool _S1311 = undistort_point_3(uv_77, &_S1310, int(12), &uv_u_38);
            if(!_S1311)
            {
                float3  _S1312 = make_float3 (0.0f);
                _S1307 = int(0);
                _S1306 = nullptr;
                _S1305 = nullptr;
                _S1304 = nullptr;
                _S1303 = nullptr;
                normal_17 = _S1312;
                break;
            }
            points_14[int(0)] = make_float3 (depths_12.x) * unproject_raydir_0(uv_u_38, camera_model_29, is_ray_depth_26);
            for(;;)
            {
                float2  uv_78 = (pix_center_15 + make_float2 (1.0f, -0.0f) - _S1308) / _S1309;
                FixedArray<float, 8>  _S1313 = dist_coeffs_31;
                float2  uv_u_39;
                bool _S1314 = undistort_point_3(uv_78, &_S1313, int(12), &uv_u_39);
                if(!_S1314)
                {
                    float3  _S1315 = make_float3 (0.0f);
                    _S1307 = int(0);
                    _S1306 = nullptr;
                    normal_17 = _S1315;
                    break;
                }
                points_14[int(1)] = make_float3 (depths_12.y) * unproject_raydir_0(uv_u_39, camera_model_29, is_ray_depth_26);
                _S1307 = int(2);
                _S1306 = &points_14[int(1)];
                break;
            }
            if(_S1307 != int(2))
            {
                _S1305 = &points_14[int(0)];
                _S1304 = nullptr;
                _S1303 = nullptr;
                break;
            }
            float2  uv_79 = (pix_center_15 + make_float2 (0.0f, -1.0f) - _S1308) / _S1309;
            FixedArray<float, 8>  _S1316 = dist_coeffs_31;
            float2  uv_u_40;
            bool _S1317 = undistort_point_3(uv_79, &_S1316, int(12), &uv_u_40);
            if(!_S1317)
            {
                float3  _S1318 = make_float3 (0.0f);
                _S1307 = int(0);
                _S1305 = &points_14[int(0)];
                _S1304 = nullptr;
                _S1303 = nullptr;
                normal_17 = _S1318;
                break;
            }
            points_14[int(2)] = make_float3 (depths_12.z) * unproject_raydir_0(uv_u_40, camera_model_29, is_ray_depth_26);
            for(;;)
            {
                float2  uv_80 = (pix_center_15 + make_float2 (0.0f, 1.0f) - _S1308) / _S1309;
                FixedArray<float, 8>  _S1319 = dist_coeffs_31;
                float2  uv_u_41;
                bool _S1320 = undistort_point_3(uv_80, &_S1319, int(12), &uv_u_41);
                if(!_S1320)
                {
                    float3  _S1321 = make_float3 (0.0f);
                    _S1307 = int(0);
                    _S1305 = nullptr;
                    normal_17 = _S1321;
                    break;
                }
                points_14[int(3)] = make_float3 (depths_12.w) * unproject_raydir_0(uv_u_41, camera_model_29, is_ray_depth_26);
                _S1307 = int(2);
                _S1305 = &points_14[int(3)];
                break;
            }
            if(_S1307 != int(2))
            {
                float3  * _S1322 = _S1305;
                _S1305 = &points_14[int(0)];
                _S1304 = _S1322;
                _S1303 = &points_14[int(2)];
                break;
            }
            float3  * _S1323 = _S1305;
            _S1307 = int(1);
            _S1305 = &points_14[int(0)];
            _S1304 = _S1323;
            _S1303 = &points_14[int(2)];
            break;
        }
        if(_S1307 != int(1))
        {
            break;
        }
        float3  normal_18 = cross_0(*_S1306 - *_S1305, - (*_S1304 - *_S1303));
        if((dot_0(normal_18, normal_18)) != 0.0f)
        {
            normal_17 = normal_18 / make_float3 (length_0(normal_18));
        }
        else
        {
            normal_17 = normal_18;
        }
        break;
    }
    return normal_17;
}

struct s_bwd_prop_depth_to_normal_Intermediates_3
{
    float2  _S1324;
    bool _S1325;
    float2  _S1326;
    bool _S1327;
    float2  _S1328;
    bool _S1329;
    float2  _S1330;
    bool _S1331;
};

inline __device__ void depth_to_normal_vjp_rational(float2  pix_center_16, float4  intrins_28, FixedArray<float, 8>  dist_coeffs_32, int camera_model_30, bool is_ray_depth_27, float4  depths_13, float3  v_normal_4, float4  * v_depths_6)
{
    float2  _S1332 = make_float2 (0.0f);
    s_bwd_prop_depth_to_normal_Intermediates_3 _S1333;
    (&_S1333)->_S1324 = _S1332;
    (&_S1333)->_S1325 = false;
    (&_S1333)->_S1326 = _S1332;
    (&_S1333)->_S1327 = false;
    (&_S1333)->_S1328 = _S1332;
    (&_S1333)->_S1329 = false;
    (&_S1333)->_S1330 = _S1332;
    (&_S1333)->_S1331 = false;
    (&_S1333)->_S1324 = _S1332;
    (&_S1333)->_S1325 = false;
    (&_S1333)->_S1326 = _S1332;
    (&_S1333)->_S1327 = false;
    (&_S1333)->_S1328 = _S1332;
    (&_S1333)->_S1329 = false;
    (&_S1333)->_S1330 = _S1332;
    (&_S1333)->_S1331 = false;
    bool _S1334 = (depths_13.x) == 0.0f;
    bool _runFlag_24;
    if(_S1334)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.y) == 0.0f;
    }
    if(_runFlag_24)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.z) == 0.0f;
    }
    if(_runFlag_24)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.w) == 0.0f;
    }
    int _S1335;
    if(!_runFlag_24)
    {
        float2  _S1336 = float2 {intrins_28.z, intrins_28.w};
        float2  _S1337 = float2 {intrins_28.x, intrins_28.y};
        float2  uv_81 = (pix_center_16 + make_float2 (-1.0f, -0.0f) - _S1336) / _S1337;
        float2  _S1338 = _S1332;
        FixedArray<float, 8>  _S1339 = dist_coeffs_32;
        bool _S1340 = undistort_point_3(uv_81, &_S1339, int(12), &_S1338);
        (&_S1333)->_S1324 = _S1338;
        (&_S1333)->_S1325 = _S1340;
        bool _S1341 = !!_S1340;
        if(_S1341)
        {
            float2  uv_82 = (pix_center_16 + make_float2 (1.0f, -0.0f) - _S1336) / _S1337;
            float2  _S1342 = _S1332;
            FixedArray<float, 8>  _S1343 = dist_coeffs_32;
            bool _S1344 = undistort_point_3(uv_82, &_S1343, int(12), &_S1342);
            (&_S1333)->_S1326 = _S1342;
            (&_S1333)->_S1327 = _S1344;
            if(!!_S1344)
            {
                _S1335 = int(2);
            }
            else
            {
                _S1335 = int(0);
            }
            if(_S1335 != int(2))
            {
                _runFlag_24 = false;
            }
            else
            {
                _runFlag_24 = _S1341;
            }
            if(_runFlag_24)
            {
                float2  uv_83 = (pix_center_16 + make_float2 (0.0f, -1.0f) - _S1336) / _S1337;
                float2  _S1345 = _S1332;
                FixedArray<float, 8>  _S1346 = dist_coeffs_32;
                bool _S1347 = undistort_point_3(uv_83, &_S1346, int(12), &_S1345);
                (&_S1333)->_S1328 = _S1345;
                (&_S1333)->_S1329 = _S1347;
                if(!_S1347)
                {
                    _runFlag_24 = false;
                }
                if(_runFlag_24)
                {
                    float2  uv_84 = (pix_center_16 + make_float2 (0.0f, 1.0f) - _S1336) / _S1337;
                    float2  _S1348 = _S1332;
                    FixedArray<float, 8>  _S1349 = dist_coeffs_32;
                    bool _S1350 = undistort_point_3(uv_84, &_S1349, int(12), &_S1348);
                    (&_S1333)->_S1330 = _S1348;
                    (&_S1333)->_S1331 = _S1350;
                }
            }
        }
    }
    s_bwd_prop_depth_to_normal_Intermediates_3 _S1351 = _S1333;
    float3  _S1352 = make_float3 (0.0f);
    if(_S1334)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.y) == 0.0f;
    }
    if(_runFlag_24)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.z) == 0.0f;
    }
    if(_runFlag_24)
    {
        _runFlag_24 = true;
    }
    else
    {
        _runFlag_24 = (depths_13.w) == 0.0f;
    }
    bool _S1353 = !_runFlag_24;
    bool _runFlag_25;
    bool _runFlag_26;
    bool _S1354;
    bool _runFlag_27;
    bool _S1355;
    bool _S1356;
    FixedArray<float3 , 4>  points_15;
    float3  _S1357;
    float3  _S1358;
    float3  _S1359;
    float3  _S1360;
    float3  _S1361;
    float3  _S1362;
    float3  _S1363;
    float3  _S1364;
    float3  _S1365;
    if(_S1353)
    {
        bool _S1366 = !!_S1351._S1325;
        if(_S1366)
        {
            float3  _S1367 = s_primal_ctx_unproject_raydir_0(_S1351._S1324, camera_model_30, is_ray_depth_27);
            float3  _S1368 = make_float3 (depths_13.x) * _S1367;
            bool _S1369 = !!_S1351._S1327;
            if(_S1369)
            {
                float3  _S1370 = s_primal_ctx_unproject_raydir_0(_S1351._S1326, camera_model_30, is_ray_depth_27);
                float3  _S1371 = make_float3 (depths_13.y) * _S1370;
                _S1335 = int(2);
                points_15[int(0)] = _S1368;
                points_15[int(1)] = _S1371;
                points_15[int(2)] = _S1352;
                points_15[int(3)] = _S1352;
                _S1357 = _S1370;
            }
            else
            {
                _S1335 = int(0);
                points_15[int(0)] = _S1368;
                points_15[int(1)] = _S1352;
                points_15[int(2)] = _S1352;
                points_15[int(3)] = _S1352;
                _S1357 = _S1352;
            }
            if(_S1335 != int(2))
            {
                _runFlag_24 = false;
            }
            else
            {
                _runFlag_24 = _S1366;
                _S1335 = int(0);
            }
            if(_runFlag_24)
            {
                if(!_S1351._S1329)
                {
                    _runFlag_25 = false;
                    _S1335 = int(0);
                }
                else
                {
                    _runFlag_25 = _runFlag_24;
                }
                if(_runFlag_25)
                {
                    float3  _S1372 = s_primal_ctx_unproject_raydir_0(_S1351._S1328, camera_model_30, is_ray_depth_27);
                    points_15[int(2)] = make_float3 (depths_13.z) * _S1372;
                    bool _S1373 = !!_S1351._S1331;
                    int _S1374;
                    if(_S1373)
                    {
                        float3  _S1375 = s_primal_ctx_unproject_raydir_0(_S1351._S1330, camera_model_30, is_ray_depth_27);
                        points_15[int(3)] = make_float3 (depths_13.w) * _S1375;
                        _S1374 = int(2);
                        _S1358 = _S1375;
                    }
                    else
                    {
                        _S1374 = int(0);
                        _S1358 = _S1352;
                    }
                    if(_S1374 != int(2))
                    {
                        _runFlag_26 = false;
                        _S1335 = _S1374;
                    }
                    else
                    {
                        _runFlag_26 = _runFlag_25;
                    }
                    if(_runFlag_26)
                    {
                        _S1335 = int(1);
                    }
                    _runFlag_26 = _S1373;
                    _S1359 = _S1372;
                }
                else
                {
                    _runFlag_26 = false;
                    _S1358 = _S1352;
                    _S1359 = _S1352;
                }
            }
            else
            {
                _runFlag_25 = false;
                _runFlag_26 = false;
                _S1358 = _S1352;
                _S1359 = _S1352;
            }
            float3  _S1376 = _S1357;
            _S1357 = _S1358;
            _S1358 = _S1359;
            _S1354 = _S1369;
            _S1359 = _S1376;
            _S1360 = _S1367;
        }
        else
        {
            _S1335 = int(0);
            points_15[int(0)] = _S1352;
            points_15[int(1)] = _S1352;
            points_15[int(2)] = _S1352;
            points_15[int(3)] = _S1352;
            _runFlag_24 = false;
            _runFlag_25 = false;
            _runFlag_26 = false;
            _S1357 = _S1352;
            _S1358 = _S1352;
            _S1354 = false;
            _S1359 = _S1352;
            _S1360 = _S1352;
        }
        if(_S1335 != int(1))
        {
            _runFlag_27 = false;
        }
        else
        {
            _runFlag_27 = _S1353;
        }
        if(_runFlag_27)
        {
            float3  dx_7 = points_15[int(1)] - points_15[int(0)];
            float3  _S1377 = - (points_15[int(3)] - points_15[int(2)]);
            float3  _S1378 = s_primal_ctx_cross_0(dx_7, _S1377);
            bool _S1379 = (s_primal_ctx_dot_0(_S1378, _S1378)) != 0.0f;
            if(_S1379)
            {
                float _S1380 = length_0(_S1378);
                float3  _S1381 = make_float3 (_S1380);
                _S1361 = make_float3 (_S1380 * _S1380);
                _S1362 = _S1381;
            }
            else
            {
                _S1361 = _S1352;
                _S1362 = _S1352;
            }
            float3  _S1382 = _S1362;
            _S1355 = _S1379;
            _S1362 = _S1378;
            _S1363 = _S1382;
            _S1364 = dx_7;
            _S1365 = _S1377;
        }
        else
        {
            _S1355 = false;
            _S1361 = _S1352;
            _S1362 = _S1352;
            _S1363 = _S1352;
            _S1364 = _S1352;
            _S1365 = _S1352;
        }
        bool _S1383 = _runFlag_24;
        bool _S1384 = _runFlag_25;
        bool _S1385 = _runFlag_26;
        float3  _S1386 = _S1357;
        float3  _S1387 = _S1358;
        bool _S1388 = _S1354;
        float3  _S1389 = _S1359;
        float3  _S1390 = _S1360;
        _runFlag_24 = _runFlag_27;
        _runFlag_25 = _S1355;
        _S1357 = _S1361;
        _S1358 = _S1362;
        _S1359 = _S1363;
        _S1360 = _S1364;
        _S1361 = _S1365;
        _runFlag_26 = _S1366;
        _S1354 = _S1383;
        _runFlag_27 = _S1384;
        _S1355 = _S1385;
        _S1362 = _S1386;
        _S1363 = _S1387;
        _S1356 = _S1388;
        _S1364 = _S1389;
        _S1365 = _S1390;
    }
    else
    {
        _runFlag_24 = false;
        _runFlag_25 = false;
        _S1357 = _S1352;
        _S1358 = _S1352;
        _S1359 = _S1352;
        _S1360 = _S1352;
        _S1361 = _S1352;
        _runFlag_26 = false;
        _S1354 = false;
        _runFlag_27 = false;
        _S1355 = false;
        _S1362 = _S1352;
        _S1363 = _S1352;
        _S1356 = false;
        _S1364 = _S1352;
        _S1365 = _S1352;
    }
    float4  _S1391 = make_float4 (0.0f);
    float4  _S1392;
    if(_S1353)
    {
        if(_runFlag_24)
        {
            if(_runFlag_25)
            {
                float3  _S1393 = v_normal_4 / _S1357;
                float3  _S1394 = _S1358 * - _S1393;
                float3  _S1395 = _S1359 * _S1393;
                float _S1396 = _S1394.x + _S1394.y + _S1394.z;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S1397;
                (&_S1397)->primal_0 = _S1358;
                (&_S1397)->differential_0 = _S1352;
                s_bwd_length_impl_0(&_S1397, _S1396);
                _S1357 = _S1395 + _S1397.differential_0;
            }
            else
            {
                _S1357 = v_normal_4;
            }
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1398;
            (&_S1398)->primal_0 = _S1358;
            (&_S1398)->differential_0 = _S1352;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1399;
            (&_S1399)->primal_0 = _S1358;
            (&_S1399)->differential_0 = _S1352;
            s_bwd_prop_dot_0(&_S1398, &_S1399, 0.0f);
            float3  _S1400 = _S1399.differential_0 + _S1398.differential_0 + _S1357;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1401;
            (&_S1401)->primal_0 = _S1360;
            (&_S1401)->differential_0 = _S1352;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1402;
            (&_S1402)->primal_0 = _S1361;
            (&_S1402)->differential_0 = _S1352;
            s_bwd_prop_cross_0(&_S1401, &_S1402, _S1400);
            float3  s_diff_dy_T_7 = - _S1402.differential_0;
            float3  _S1403 = - s_diff_dy_T_7;
            float3  _S1404 = - _S1401.differential_0;
            FixedArray<float3 , 4>  _S1405;
            _S1405[int(0)] = _S1352;
            _S1405[int(1)] = _S1352;
            _S1405[int(2)] = _S1352;
            _S1405[int(3)] = _S1352;
            _S1405[int(2)] = _S1403;
            _S1405[int(3)] = s_diff_dy_T_7;
            _S1405[int(0)] = _S1404;
            _S1405[int(1)] = _S1401.differential_0;
            points_15[int(0)] = _S1405[int(0)];
            points_15[int(1)] = _S1405[int(1)];
            points_15[int(2)] = _S1405[int(2)];
            points_15[int(3)] = _S1405[int(3)];
        }
        else
        {
            points_15[int(0)] = _S1352;
            points_15[int(1)] = _S1352;
            points_15[int(2)] = _S1352;
            points_15[int(3)] = _S1352;
        }
        if(_runFlag_26)
        {
            if(_S1354)
            {
                if(_runFlag_27)
                {
                    FixedArray<float3 , 4>  _S1406 = points_15;
                    FixedArray<float3 , 4>  _S1407 = points_15;
                    FixedArray<float3 , 4>  _S1408 = points_15;
                    FixedArray<float3 , 4>  _S1409 = points_15;
                    if(_S1355)
                    {
                        float3  _S1410 = _S1362 * _S1409[int(3)];
                        float _S1411 = _S1410.x + _S1410.y + _S1410.z;
                        float4  _S1412 = _S1391;
                        *&((&_S1412)->w) = _S1411;
                        points_15[int(0)] = _S1406[int(0)];
                        points_15[int(1)] = _S1407[int(1)];
                        points_15[int(2)] = _S1408[int(2)];
                        points_15[int(3)] = _S1352;
                        _S1392 = _S1412;
                    }
                    else
                    {
                        points_15[int(0)] = _S1406[int(0)];
                        points_15[int(1)] = _S1407[int(1)];
                        points_15[int(2)] = _S1408[int(2)];
                        points_15[int(3)] = _S1409[int(3)];
                        _S1392 = _S1391;
                    }
                    float3  _S1413 = _S1363 * points_15[int(2)];
                    float _S1414 = _S1413.x + _S1413.y + _S1413.z;
                    FixedArray<float3 , 4>  _S1415 = points_15;
                    FixedArray<float3 , 4>  _S1416 = points_15;
                    float4  _S1417 = _S1391;
                    *&((&_S1417)->z) = _S1414;
                    float4  _S1418 = _S1392 + _S1417;
                    points_15[int(0)] = points_15[int(0)];
                    points_15[int(1)] = _S1415[int(1)];
                    points_15[int(2)] = _S1352;
                    points_15[int(3)] = _S1416[int(3)];
                    _S1392 = _S1418;
                }
                else
                {
                    FixedArray<float3 , 4>  _S1419 = points_15;
                    FixedArray<float3 , 4>  _S1420 = points_15;
                    FixedArray<float3 , 4>  _S1421 = points_15;
                    points_15[int(0)] = points_15[int(0)];
                    points_15[int(1)] = _S1419[int(1)];
                    points_15[int(2)] = _S1420[int(2)];
                    points_15[int(3)] = _S1421[int(3)];
                    _S1392 = _S1391;
                }
            }
            else
            {
                FixedArray<float3 , 4>  _S1422 = points_15;
                FixedArray<float3 , 4>  _S1423 = points_15;
                FixedArray<float3 , 4>  _S1424 = points_15;
                points_15[int(0)] = points_15[int(0)];
                points_15[int(1)] = _S1422[int(1)];
                points_15[int(2)] = _S1423[int(2)];
                points_15[int(3)] = _S1424[int(3)];
                _S1392 = _S1391;
            }
            if(_S1356)
            {
                FixedArray<float3 , 4>  _S1425 = points_15;
                float3  _S1426 = _S1364 * points_15[int(1)];
                float _S1427 = _S1426.x + _S1426.y + _S1426.z;
                float4  _S1428 = _S1391;
                *&((&_S1428)->y) = _S1427;
                float4  _S1429 = _S1392 + _S1428;
                points_15[int(0)] = _S1352;
                points_15[int(1)] = _S1352;
                points_15[int(2)] = _S1352;
                points_15[int(3)] = _S1352;
                _S1357 = _S1425[int(0)];
                _S1392 = _S1429;
            }
            else
            {
                FixedArray<float3 , 4>  _S1430 = points_15;
                FixedArray<float3 , 4>  _S1431 = points_15;
                FixedArray<float3 , 4>  _S1432 = points_15;
                points_15[int(0)] = points_15[int(0)];
                points_15[int(1)] = _S1430[int(1)];
                points_15[int(2)] = _S1431[int(2)];
                points_15[int(3)] = _S1432[int(3)];
                _S1357 = _S1352;
            }
            float3  _S1433 = _S1365 * (points_15[int(0)] + _S1357);
            float _S1434 = _S1433.x + _S1433.y + _S1433.z;
            float4  _S1435 = _S1391;
            *&((&_S1435)->x) = _S1434;
            _S1392 = _S1392 + _S1435;
        }
        else
        {
            _S1392 = _S1391;
        }
    }
    else
    {
        _S1392 = _S1391;
    }
    *v_depths_6 = _S1392;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor_rational(float2  pix_center_17, float4  intrins_29, FixedArray<float, 8>  dist_coeffs_33, int camera_model_31)
{
    float _S1436;
    for(;;)
    {
        float2  uv_85 = (pix_center_17 - float2 {intrins_29.z, intrins_29.w}) / float2 {intrins_29.x, intrins_29.y};
        FixedArray<float, 8>  _S1437 = dist_coeffs_33;
        float2  uv_u_42;
        bool _S1438 = undistort_point_3(uv_85, &_S1437, int(12), &uv_u_42);
        if(!_S1438)
        {
            _S1436 = 0.0f;
            break;
        }
        float3  raydir_23 = unproject_raydir_0(uv_u_42, camera_model_31, false);
        _S1436 = float((F32_sign((raydir_23.z)))) / length_0(raydir_23);
        break;
    }
    return _S1436;
}

inline __device__ float depth_normal_loss_rational(float2  pix_center_18, float4  intrins_30, FixedArray<float, 8>  dist_coeffs_34, int camera_model_32, bool is_ray_depth_28, float4  depths_14, float3  gt_normal_6)
{
    float _S1439;
    for(;;)
    {
        float3  _S1440;
        float3  * _S1441;
        float3  * _S1442;
        float3  * _S1443;
        float3  * _S1444;
        int _S1445;
        FixedArray<float3 , 5>  points_16;
        for(;;)
        {
            float2  _S1446 = float2 {intrins_30.z, intrins_30.w};
            float2  _S1447 = float2 {intrins_30.x, intrins_30.y};
            float2  uv_86 = (pix_center_18 + make_float2 (-1.0f, -0.0f) - _S1446) / _S1447;
            FixedArray<float, 8>  _S1448 = dist_coeffs_34;
            float2  uv_u_43;
            bool _S1449 = undistort_point_3(uv_86, &_S1448, int(12), &uv_u_43);
            float3  _S1450 = make_float3 (0.0f);
            if(!_S1449)
            {
                _S1445 = int(0);
                _S1444 = nullptr;
                _S1443 = nullptr;
                _S1442 = nullptr;
                _S1441 = nullptr;
                _S1440 = _S1450;
                break;
            }
            float3  raydir_24 = unproject_raydir_0(uv_u_43, camera_model_32, is_ray_depth_28);
            points_16[int(0)] = make_float3 (depths_14.x) * raydir_24;
            float2  uv_87 = (pix_center_18 + make_float2 (1.0f, -0.0f) - _S1446) / _S1447;
            FixedArray<float, 8>  _S1451 = dist_coeffs_34;
            float2  uv_u_44;
            bool _S1452 = undistort_point_3(uv_87, &_S1451, int(12), &uv_u_44);
            if(!_S1452)
            {
                _S1445 = int(0);
                _S1444 = nullptr;
                _S1443 = &points_16[int(0)];
                _S1442 = nullptr;
                _S1441 = nullptr;
                _S1440 = _S1450;
                break;
            }
            float3  raydir_25 = unproject_raydir_0(uv_u_44, camera_model_32, is_ray_depth_28);
            points_16[int(1)] = make_float3 (depths_14.y) * raydir_25;
            float2  uv_88 = (pix_center_18 + make_float2 (0.0f, -1.0f) - _S1446) / _S1447;
            FixedArray<float, 8>  _S1453 = dist_coeffs_34;
            float2  uv_u_45;
            bool _S1454 = undistort_point_3(uv_88, &_S1453, int(12), &uv_u_45);
            if(!_S1454)
            {
                _S1445 = int(0);
                _S1444 = &points_16[int(1)];
                _S1443 = &points_16[int(0)];
                _S1442 = nullptr;
                _S1441 = nullptr;
                _S1440 = _S1450;
                break;
            }
            float3  raydir_26 = unproject_raydir_0(uv_u_45, camera_model_32, is_ray_depth_28);
            points_16[int(2)] = make_float3 (depths_14.z) * raydir_26;
            float2  uv_89 = (pix_center_18 + make_float2 (0.0f, 1.0f) - _S1446) / _S1447;
            FixedArray<float, 8>  _S1455 = dist_coeffs_34;
            float2  uv_u_46;
            bool _S1456 = undistort_point_3(uv_89, &_S1455, int(12), &uv_u_46);
            if(!_S1456)
            {
                _S1445 = int(0);
                _S1444 = &points_16[int(1)];
                _S1443 = &points_16[int(0)];
                _S1442 = nullptr;
                _S1441 = &points_16[int(2)];
                _S1440 = _S1450;
                break;
            }
            float3  raydir_27 = unproject_raydir_0(uv_u_46, camera_model_32, is_ray_depth_28);
            points_16[int(3)] = make_float3 (depths_14.w) * raydir_27;
            float2  uv_90 = (pix_center_18 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S1446) / _S1447;
            FixedArray<float, 8>  _S1457 = dist_coeffs_34;
            float2  uv_u_47;
            bool _S1458 = undistort_point_3(uv_90, &_S1457, int(12), &uv_u_47);
            if(!_S1458)
            {
                _S1445 = int(0);
                _S1444 = &points_16[int(1)];
                _S1443 = &points_16[int(0)];
                _S1442 = &points_16[int(3)];
                _S1441 = &points_16[int(2)];
                _S1440 = _S1450;
                break;
            }
            float3  raydir_28 = unproject_raydir_0(uv_u_47, camera_model_32, is_ray_depth_28);
            _S1445 = int(1);
            _S1444 = &points_16[int(1)];
            _S1443 = &points_16[int(0)];
            _S1442 = &points_16[int(3)];
            _S1441 = &points_16[int(2)];
            _S1440 = raydir_28;
            break;
        }
        if(_S1445 != int(1))
        {
            _S1439 = 0.0f;
            break;
        }
        float3  normal_19 = cross_0(*_S1444 - *_S1443, - (*_S1442 - *_S1441));
        float3  normal_20;
        if((dot_0(normal_19, normal_19)) != 0.0f)
        {
            normal_20 = normalize_0(normal_19);
        }
        else
        {
            normal_20 = normal_19;
        }
        float3  _S1459;
        if((dot_0(gt_normal_6, gt_normal_6)) != 0.0f)
        {
            _S1459 = normalize_0(gt_normal_6);
        }
        else
        {
            _S1459 = gt_normal_6;
        }
        _S1439 = (1.0f - dot_0(normal_20, _S1459) + 0.00100000004749745f) / ((F32_max((dot_0(normal_20, - normalize_0(_S1440))), (0.0f))) + 0.00100000004749745f);
        break;
    }
    return _S1439;
}

struct s_bwd_prop_depth_normal_loss_Intermediates_3
{
    float2  _S1460;
    bool _S1461;
    float2  _S1462;
    bool _S1463;
    float2  _S1464;
    bool _S1465;
    float2  _S1466;
    bool _S1467;
    float2  _S1468;
    bool _S1469;
};

inline __device__ void depth_normal_loss_vjp_rational(float2  pix_center_19, float4  intrins_31, FixedArray<float, 8>  dist_coeffs_35, int camera_model_33, bool is_ray_depth_29, float4  depths_15, float3  gt_normal_7, float v_loss_3, float4  * v_depths_7, float3  * v_gt_normal_3)
{
    float2  _S1470 = make_float2 (0.0f);
    s_bwd_prop_depth_normal_loss_Intermediates_3 _S1471;
    (&_S1471)->_S1460 = _S1470;
    (&_S1471)->_S1461 = false;
    (&_S1471)->_S1462 = _S1470;
    (&_S1471)->_S1463 = false;
    (&_S1471)->_S1464 = _S1470;
    (&_S1471)->_S1465 = false;
    (&_S1471)->_S1466 = _S1470;
    (&_S1471)->_S1467 = false;
    (&_S1471)->_S1468 = _S1470;
    (&_S1471)->_S1469 = false;
    (&_S1471)->_S1462 = _S1470;
    (&_S1471)->_S1463 = false;
    (&_S1471)->_S1464 = _S1470;
    (&_S1471)->_S1465 = false;
    (&_S1471)->_S1466 = _S1470;
    (&_S1471)->_S1467 = false;
    (&_S1471)->_S1468 = _S1470;
    (&_S1471)->_S1469 = false;
    float2  _S1472 = float2 {intrins_31.z, intrins_31.w};
    float2  _S1473 = float2 {intrins_31.x, intrins_31.y};
    float2  uv_91 = (pix_center_19 + make_float2 (-1.0f, -0.0f) - _S1472) / _S1473;
    float2  _S1474 = _S1470;
    FixedArray<float, 8>  _S1475 = dist_coeffs_35;
    bool _S1476 = undistort_point_3(uv_91, &_S1475, int(12), &_S1474);
    (&_S1471)->_S1460 = _S1474;
    (&_S1471)->_S1461 = _S1476;
    bool _S1477 = !!_S1476;
    bool _runFlag_28;
    if(_S1477)
    {
        float2  uv_92 = (pix_center_19 + make_float2 (1.0f, -0.0f) - _S1472) / _S1473;
        float2  _S1478 = _S1470;
        FixedArray<float, 8>  _S1479 = dist_coeffs_35;
        bool _S1480 = undistort_point_3(uv_92, &_S1479, int(12), &_S1478);
        (&_S1471)->_S1462 = _S1478;
        (&_S1471)->_S1463 = _S1480;
        if(!_S1480)
        {
            _runFlag_28 = false;
        }
        else
        {
            _runFlag_28 = _S1477;
        }
        if(_runFlag_28)
        {
            float2  uv_93 = (pix_center_19 + make_float2 (0.0f, -1.0f) - _S1472) / _S1473;
            float2  _S1481 = _S1470;
            FixedArray<float, 8>  _S1482 = dist_coeffs_35;
            bool _S1483 = undistort_point_3(uv_93, &_S1482, int(12), &_S1481);
            (&_S1471)->_S1464 = _S1481;
            (&_S1471)->_S1465 = _S1483;
            if(!_S1483)
            {
                _runFlag_28 = false;
            }
            if(_runFlag_28)
            {
                float2  uv_94 = (pix_center_19 + make_float2 (0.0f, 1.0f) - _S1472) / _S1473;
                float2  _S1484 = _S1470;
                FixedArray<float, 8>  _S1485 = dist_coeffs_35;
                bool _S1486 = undistort_point_3(uv_94, &_S1485, int(12), &_S1484);
                (&_S1471)->_S1466 = _S1484;
                (&_S1471)->_S1467 = _S1486;
                if(!_S1486)
                {
                    _runFlag_28 = false;
                }
                if(_runFlag_28)
                {
                    float2  uv_95 = (pix_center_19 - _S1472) / _S1473;
                    float2  _S1487 = _S1470;
                    FixedArray<float, 8>  _S1488 = dist_coeffs_35;
                    bool _S1489 = undistort_point_3(uv_95, &_S1488, int(12), &_S1487);
                    (&_S1471)->_S1468 = _S1487;
                    (&_S1471)->_S1469 = _S1489;
                }
            }
        }
    }
    s_bwd_prop_depth_normal_loss_Intermediates_3 _S1490 = _S1471;
    float3  _S1491 = make_float3 (0.0f);
    bool _S1492 = !!_S1471._S1461;
    bool _runFlag_29;
    bool _runFlag_30;
    bool _runFlag_31;
    int _S1493;
    float3  raydir_29;
    float3  _S1494;
    float3  _S1495;
    float3  _S1496;
    float3  _S1497;
    FixedArray<float3 , 5>  points_17;
    if(_S1492)
    {
        float3  _S1498 = s_primal_ctx_unproject_raydir_0(_S1490._S1460, camera_model_33, is_ray_depth_29);
        float3  _S1499 = make_float3 (depths_15.x) * _S1498;
        if(!_S1490._S1463)
        {
            _runFlag_28 = false;
        }
        else
        {
            _runFlag_28 = _S1492;
        }
        if(_runFlag_28)
        {
            float3  _S1500 = s_primal_ctx_unproject_raydir_0(_S1490._S1462, camera_model_33, is_ray_depth_29);
            float3  _S1501 = make_float3 (depths_15.y) * _S1500;
            if(!_S1490._S1465)
            {
                _runFlag_29 = false;
            }
            else
            {
                _runFlag_29 = _runFlag_28;
            }
            if(_runFlag_29)
            {
                float3  _S1502 = s_primal_ctx_unproject_raydir_0(_S1490._S1464, camera_model_33, is_ray_depth_29);
                float3  _S1503 = make_float3 (depths_15.z) * _S1502;
                if(!_S1490._S1467)
                {
                    _runFlag_30 = false;
                }
                else
                {
                    _runFlag_30 = _runFlag_29;
                }
                if(_runFlag_30)
                {
                    float3  _S1504 = s_primal_ctx_unproject_raydir_0(_S1490._S1466, camera_model_33, is_ray_depth_29);
                    float3  _S1505 = make_float3 (depths_15.w) * _S1504;
                    if(!_S1490._S1469)
                    {
                        _runFlag_31 = false;
                    }
                    else
                    {
                        _runFlag_31 = _runFlag_30;
                    }
                    if(_runFlag_31)
                    {
                        float3  _S1506 = s_primal_ctx_unproject_raydir_0(_S1490._S1468, camera_model_33, is_ray_depth_29);
                        _S1493 = int(1);
                        raydir_29 = _S1506;
                    }
                    else
                    {
                        _S1493 = int(0);
                        raydir_29 = _S1504;
                    }
                    points_17[int(0)] = _S1499;
                    points_17[int(1)] = _S1501;
                    points_17[int(2)] = _S1503;
                    points_17[int(3)] = _S1505;
                    points_17[int(4)] = _S1491;
                    _S1494 = _S1504;
                }
                else
                {
                    _S1493 = int(0);
                    raydir_29 = _S1502;
                    points_17[int(0)] = _S1499;
                    points_17[int(1)] = _S1501;
                    points_17[int(2)] = _S1503;
                    points_17[int(3)] = _S1491;
                    points_17[int(4)] = _S1491;
                    _S1494 = _S1491;
                }
                _S1495 = _S1502;
            }
            else
            {
                _S1493 = int(0);
                raydir_29 = _S1500;
                points_17[int(0)] = _S1499;
                points_17[int(1)] = _S1501;
                points_17[int(2)] = _S1491;
                points_17[int(3)] = _S1491;
                points_17[int(4)] = _S1491;
                _runFlag_30 = false;
                _S1494 = _S1491;
                _S1495 = _S1491;
            }
            _S1496 = _S1500;
        }
        else
        {
            _S1493 = int(0);
            raydir_29 = _S1498;
            points_17[int(0)] = _S1499;
            points_17[int(1)] = _S1491;
            points_17[int(2)] = _S1491;
            points_17[int(3)] = _S1491;
            points_17[int(4)] = _S1491;
            _runFlag_29 = false;
            _runFlag_30 = false;
            _S1494 = _S1491;
            _S1495 = _S1491;
            _S1496 = _S1491;
        }
        _S1497 = _S1498;
    }
    else
    {
        _S1493 = int(0);
        points_17[int(0)] = _S1491;
        points_17[int(1)] = _S1491;
        points_17[int(2)] = _S1491;
        points_17[int(3)] = _S1491;
        points_17[int(4)] = _S1491;
        _runFlag_28 = false;
        _runFlag_29 = false;
        _runFlag_30 = false;
        _S1494 = _S1491;
        _S1495 = _S1491;
        _S1496 = _S1491;
        _S1497 = _S1491;
    }
    bool _S1507 = !(_S1493 != int(1));
    bool _S1508;
    float3  normal_21;
    float3  _S1509;
    float3  _S1510;
    float3  _S1511;
    float3  _S1512;
    float _S1513;
    float _S1514;
    float _S1515;
    float _S1516;
    if(_S1507)
    {
        float3  dx_8 = points_17[int(1)] - points_17[int(0)];
        float3  _S1517 = - (points_17[int(3)] - points_17[int(2)]);
        float3  _S1518 = s_primal_ctx_cross_0(dx_8, _S1517);
        bool _S1519 = (s_primal_ctx_dot_0(_S1518, _S1518)) != 0.0f;
        if(_S1519)
        {
            normal_21 = normalize_0(_S1518);
        }
        else
        {
            normal_21 = _S1518;
        }
        bool _S1520 = (s_primal_ctx_dot_0(gt_normal_7, gt_normal_7)) != 0.0f;
        if(_S1520)
        {
            _S1509 = normalize_0(gt_normal_7);
        }
        else
        {
            _S1509 = gt_normal_7;
        }
        float3  _S1521 = - normalize_0(raydir_29);
        float _S1522 = s_primal_ctx_dot_0(normal_21, _S1521);
        float _S1523 = 1.0f - s_primal_ctx_dot_0(normal_21, _S1509) + 0.00100000004749745f;
        float _S1524 = (F32_max((_S1522), (0.0f))) + 0.00100000004749745f;
        _S1513 = _S1524 * _S1524;
        _S1514 = _S1523;
        _S1515 = _S1524;
        _S1516 = _S1522;
        raydir_29 = normal_21;
        normal_21 = _S1521;
        _runFlag_31 = _S1520;
        _S1508 = _S1519;
        _S1510 = _S1518;
        _S1511 = dx_8;
        _S1512 = _S1517;
    }
    else
    {
        _S1513 = 0.0f;
        _S1514 = 0.0f;
        _S1515 = 0.0f;
        _S1516 = 0.0f;
        raydir_29 = _S1491;
        normal_21 = _S1491;
        _S1509 = _S1491;
        _runFlag_31 = false;
        _S1508 = false;
        _S1510 = _S1491;
        _S1511 = _S1491;
        _S1512 = _S1491;
    }
    float4  _S1525 = make_float4 (0.0f);
    if(_S1507)
    {
        float _S1526 = v_loss_3 / _S1513;
        float _S1527 = _S1514 * - _S1526;
        float s_diff_num_T_3 = _S1515 * _S1526;
        DiffPair_float_0 _S1528;
        (&_S1528)->primal_0 = _S1516;
        (&_S1528)->differential_0 = 0.0f;
        DiffPair_float_0 _S1529;
        (&_S1529)->primal_0 = 0.0f;
        (&_S1529)->differential_0 = 0.0f;
        _d_max_0(&_S1528, &_S1529, _S1527);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1530;
        (&_S1530)->primal_0 = raydir_29;
        (&_S1530)->differential_0 = _S1491;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1531;
        (&_S1531)->primal_0 = normal_21;
        (&_S1531)->differential_0 = _S1491;
        s_bwd_prop_dot_0(&_S1530, &_S1531, _S1528.differential_0);
        float _S1532 = - s_diff_num_T_3;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1533;
        (&_S1533)->primal_0 = raydir_29;
        (&_S1533)->differential_0 = _S1491;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1534;
        (&_S1534)->primal_0 = _S1509;
        (&_S1534)->differential_0 = _S1491;
        s_bwd_prop_dot_0(&_S1533, &_S1534, _S1532);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1535 = _S1534;
        float3  _S1536 = _S1530.differential_0 + _S1533.differential_0;
        if(_runFlag_31)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1537;
            (&_S1537)->primal_0 = gt_normal_7;
            (&_S1537)->differential_0 = _S1491;
            s_bwd_normalize_impl_0(&_S1537, _S1535.differential_0);
            raydir_29 = _S1537.differential_0;
        }
        else
        {
            raydir_29 = _S1535.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1538;
        (&_S1538)->primal_0 = gt_normal_7;
        (&_S1538)->differential_0 = _S1491;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1539;
        (&_S1539)->primal_0 = gt_normal_7;
        (&_S1539)->differential_0 = _S1491;
        s_bwd_prop_dot_0(&_S1538, &_S1539, 0.0f);
        float3  _S1540 = _S1539.differential_0 + _S1538.differential_0 + raydir_29;
        if(_S1508)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1541;
            (&_S1541)->primal_0 = _S1510;
            (&_S1541)->differential_0 = _S1491;
            s_bwd_normalize_impl_0(&_S1541, _S1536);
            raydir_29 = _S1541.differential_0;
        }
        else
        {
            raydir_29 = _S1536;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1542;
        (&_S1542)->primal_0 = _S1510;
        (&_S1542)->differential_0 = _S1491;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1543;
        (&_S1543)->primal_0 = _S1510;
        (&_S1543)->differential_0 = _S1491;
        s_bwd_prop_dot_0(&_S1542, &_S1543, 0.0f);
        float3  _S1544 = _S1543.differential_0 + _S1542.differential_0 + raydir_29;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1545;
        (&_S1545)->primal_0 = _S1511;
        (&_S1545)->differential_0 = _S1491;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1546;
        (&_S1546)->primal_0 = _S1512;
        (&_S1546)->differential_0 = _S1491;
        s_bwd_prop_cross_0(&_S1545, &_S1546, _S1544);
        float3  s_diff_dy_T_8 = - _S1546.differential_0;
        float3  _S1547 = - s_diff_dy_T_8;
        float3  _S1548 = - _S1545.differential_0;
        FixedArray<float3 , 5>  _S1549;
        _S1549[int(0)] = _S1491;
        _S1549[int(1)] = _S1491;
        _S1549[int(2)] = _S1491;
        _S1549[int(3)] = _S1491;
        _S1549[int(4)] = _S1491;
        _S1549[int(2)] = _S1547;
        _S1549[int(3)] = s_diff_dy_T_8;
        _S1549[int(0)] = _S1548;
        _S1549[int(1)] = _S1545.differential_0;
        points_17[int(0)] = _S1549[int(0)];
        points_17[int(1)] = _S1549[int(1)];
        points_17[int(2)] = _S1549[int(2)];
        points_17[int(3)] = _S1549[int(3)];
        points_17[int(4)] = _S1549[int(4)];
        raydir_29 = _S1540;
    }
    else
    {
        points_17[int(0)] = _S1491;
        points_17[int(1)] = _S1491;
        points_17[int(2)] = _S1491;
        points_17[int(3)] = _S1491;
        points_17[int(4)] = _S1491;
        raydir_29 = _S1491;
    }
    float4  _S1550;
    if(_S1492)
    {
        if(_runFlag_28)
        {
            if(_runFlag_29)
            {
                if(_runFlag_30)
                {
                    FixedArray<float3 , 5>  _S1551 = points_17;
                    FixedArray<float3 , 5>  _S1552 = points_17;
                    FixedArray<float3 , 5>  _S1553 = points_17;
                    float3  _S1554 = _S1494 * points_17[int(3)];
                    float _S1555 = _S1554.x + _S1554.y + _S1554.z;
                    float4  _S1556 = _S1525;
                    *&((&_S1556)->w) = _S1555;
                    points_17[int(0)] = _S1491;
                    points_17[int(1)] = _S1491;
                    points_17[int(2)] = _S1491;
                    points_17[int(3)] = _S1491;
                    points_17[int(4)] = _S1491;
                    _S1494 = _S1553[int(2)];
                    normal_21 = _S1551[int(0)];
                    _S1509 = _S1552[int(1)];
                    _S1550 = _S1556;
                }
                else
                {
                    FixedArray<float3 , 5>  _S1557 = points_17;
                    FixedArray<float3 , 5>  _S1558 = points_17;
                    FixedArray<float3 , 5>  _S1559 = points_17;
                    FixedArray<float3 , 5>  _S1560 = points_17;
                    points_17[int(0)] = points_17[int(0)];
                    points_17[int(1)] = _S1557[int(1)];
                    points_17[int(2)] = _S1558[int(2)];
                    points_17[int(3)] = _S1559[int(3)];
                    points_17[int(4)] = _S1560[int(4)];
                    _S1494 = _S1491;
                    normal_21 = _S1491;
                    _S1509 = _S1491;
                    _S1550 = _S1525;
                }
                float3  _S1561 = _S1495 * (points_17[int(2)] + _S1494);
                float _S1562 = _S1561.x + _S1561.y + _S1561.z;
                float3  _S1563 = points_17[int(0)] + normal_21;
                float3  _S1564 = points_17[int(1)] + _S1509;
                float4  _S1565 = _S1525;
                *&((&_S1565)->z) = _S1562;
                float4  _S1566 = _S1550 + _S1565;
                points_17[int(0)] = _S1491;
                points_17[int(1)] = _S1491;
                points_17[int(2)] = _S1491;
                points_17[int(3)] = _S1491;
                points_17[int(4)] = _S1491;
                _S1494 = _S1564;
                _S1495 = _S1563;
                _S1550 = _S1566;
            }
            else
            {
                FixedArray<float3 , 5>  _S1567 = points_17;
                FixedArray<float3 , 5>  _S1568 = points_17;
                FixedArray<float3 , 5>  _S1569 = points_17;
                FixedArray<float3 , 5>  _S1570 = points_17;
                points_17[int(0)] = points_17[int(0)];
                points_17[int(1)] = _S1567[int(1)];
                points_17[int(2)] = _S1568[int(2)];
                points_17[int(3)] = _S1569[int(3)];
                points_17[int(4)] = _S1570[int(4)];
                _S1494 = _S1491;
                _S1495 = _S1491;
                _S1550 = _S1525;
            }
            float3  _S1571 = _S1496 * (points_17[int(1)] + _S1494);
            float _S1572 = _S1571.x + _S1571.y + _S1571.z;
            float3  _S1573 = points_17[int(0)] + _S1495;
            float4  _S1574 = _S1525;
            *&((&_S1574)->y) = _S1572;
            float4  _S1575 = _S1550 + _S1574;
            points_17[int(0)] = _S1491;
            points_17[int(1)] = _S1491;
            points_17[int(2)] = _S1491;
            points_17[int(3)] = _S1491;
            points_17[int(4)] = _S1491;
            _S1494 = _S1573;
            _S1550 = _S1575;
        }
        else
        {
            FixedArray<float3 , 5>  _S1576 = points_17;
            FixedArray<float3 , 5>  _S1577 = points_17;
            FixedArray<float3 , 5>  _S1578 = points_17;
            FixedArray<float3 , 5>  _S1579 = points_17;
            points_17[int(0)] = points_17[int(0)];
            points_17[int(1)] = _S1576[int(1)];
            points_17[int(2)] = _S1577[int(2)];
            points_17[int(3)] = _S1578[int(3)];
            points_17[int(4)] = _S1579[int(4)];
            _S1494 = _S1491;
            _S1550 = _S1525;
        }
        float3  _S1580 = _S1497 * (points_17[int(0)] + _S1494);
        float _S1581 = _S1580.x + _S1580.y + _S1580.z;
        float4  _S1582 = _S1525;
        *&((&_S1582)->x) = _S1581;
        _S1550 = _S1550 + _S1582;
    }
    else
    {
        _S1550 = _S1525;
    }
    *v_depths_7 = _S1550;
    *v_gt_normal_3 = raydir_29;
    return;
}

