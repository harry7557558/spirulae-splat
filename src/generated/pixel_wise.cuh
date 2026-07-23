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

inline __device__ float3  unproject_raydir_0(float2  uv_1, int camera_model_0, bool is_ray_depth_0)
{
    float3  raydir_0;
    bool is_unit_0;
    if(camera_model_0 == int(1))
    {
        float theta_0 = length_0(uv_1);
        float3  _S203 = make_float3 ((uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).x, (uv_1 / make_float2 ((F32_max((theta_0), (1.00000001168609742e-07f)))) * make_float2 ((F32_sin((theta_0))))).y, (F32_cos((theta_0))));
        is_unit_0 = true;
        raydir_0 = _S203;
    }
    else
    {
        bool _S204 = camera_model_0 == int(2);
        if(_S204)
        {
            float r_2 = length_0(uv_1);
            raydir_0 = make_float3 ((uv_1 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_2 * r_2)))))))).x, (uv_1 * make_float2 ((F32_sqrt(((F32_max((0.0f), (1.0f - 0.25f * r_2 * r_2)))))))).y, 1.0f - 0.5f * r_2 * r_2);
        }
        else
        {
            raydir_0 = make_float3 (uv_1.x, uv_1.y, 1.0f);
        }
        is_unit_0 = _S204;
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

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, int camera_model_1, bool is_ray_depth_1)
{
    float2  _S205 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_2 = _S205;
    FixedArray<float, 10>  _S206 = dist_coeffs_1;
    bool _S207 = undistort_point_0(_S205, &_S206, int(12), &uv_2);
    if(!_S207)
    {
        int3  _S208 = make_int3 (int(0));
        float3  _S209 = make_float3 ((float)_S208.x, (float)_S208.y, (float)_S208.z);
        return _S209;
    }
    return unproject_raydir_0(uv_2, camera_model_1, is_ray_depth_1);
}

inline __device__ float3  depth_to_point(float2  pix_pos_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, int camera_model_2, bool is_ray_depth_2, float depth_2)
{
    float2  _S210 = (pix_pos_1 - float2 {intrins_1.z, intrins_1.w}) / float2 {intrins_1.x, intrins_1.y};
    float2  uv_3 = _S210;
    FixedArray<float, 10>  _S211 = dist_coeffs_2;
    bool _S212 = undistort_point_0(_S210, &_S211, int(12), &uv_3);
    if(!_S212)
    {
        return make_float3 (0.0f);
    }
    return make_float3 (depth_2) * unproject_raydir_0(uv_3, camera_model_2, is_ray_depth_2);
}

struct s_bwd_prop_depth_to_point_Intermediates_0
{
    float2  _S213;
    bool _S214;
};

inline __device__ float s_primal_ctx_sin_0(float _S215)
{
    return (F32_sin((_S215)));
}

inline __device__ float s_primal_ctx_cos_0(float _S216)
{
    return (F32_cos((_S216)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S217)
{
    return (F32_sqrt((_S217)));
}

inline __device__ float3  s_primal_ctx_unproject_raydir_0(float2  dpuv_0, int camera_model_3, bool is_ray_depth_3)
{
    float3  raydir_1;
    bool is_unit_1;
    if(camera_model_3 == int(1))
    {
        float _S218 = length_0(dpuv_0);
        float3  _S219 = make_float3 ((dpuv_0 / make_float2 ((F32_max((_S218), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S218))).x, (dpuv_0 / make_float2 ((F32_max((_S218), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S218))).y, s_primal_ctx_cos_0(_S218));
        is_unit_1 = true;
        raydir_1 = _S219;
    }
    else
    {
        bool _S220 = camera_model_3 == int(2);
        if(_S220)
        {
            float _S221 = length_0(dpuv_0);
            raydir_1 = make_float3 ((dpuv_0 * make_float2 (s_primal_ctx_sqrt_0((F32_max((0.0f), (1.0f - 0.25f * _S221 * _S221)))))).x, (dpuv_0 * make_float2 (s_primal_ctx_sqrt_0((F32_max((0.0f), (1.0f - 0.25f * _S221 * _S221)))))).y, 1.0f - 0.5f * _S221 * _S221);
        }
        else
        {
            raydir_1 = make_float3 (dpuv_0.x, dpuv_0.y, 1.0f);
        }
        is_unit_1 = _S220;
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

inline __device__ float3  s_primal_ctx_depth_to_point_0(float2  pix_pos_2, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, int camera_model_4, bool is_ray_depth_4, float dpdepth_1, s_bwd_prop_depth_to_point_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S213 = make_float2 (0.0f);
    _s_diff_ctx_0->_S214 = false;
    float2  _S222 = (pix_pos_2 - float2 {intrins_2.z, intrins_2.w}) / float2 {intrins_2.x, intrins_2.y};
    float2  _S223 = _S222;
    bool _S224 = undistort_point_0(_S222, dist_coeffs_3, int(12), &_S223);
    _s_diff_ctx_0->_S213 = _S223;
    _s_diff_ctx_0->_S214 = _S224;
    float2  uv_4 = _S223;
    bool _S225 = !_S224;
    float3  _S226;
    if(_S225)
    {
        _S226 = make_float3 (0.0f);
    }
    bool _S227 = !_S225;
    if(_S227)
    {
        _S226 = make_float3 (dpdepth_1) * s_primal_ctx_unproject_raydir_0(uv_4, camera_model_4, is_ray_depth_4);
    }
    return _S226;
}

inline __device__ void s_bwd_prop_depth_to_point_0(float2  pix_pos_3, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, int camera_model_5, bool is_ray_depth_5, DiffPair_float_0 * dpdepth_2, float3  _s_dOut_4, s_bwd_prop_depth_to_point_Intermediates_0 * _s_diff_ctx_1)
{
    float3  _S228 = make_float3 (0.0f);
    float2  _S229 = _s_diff_ctx_1->_S213;
    bool _S230 = !!_s_diff_ctx_1->_S214;
    float3  _S231;
    if(_S230)
    {
        _S231 = s_primal_ctx_unproject_raydir_0(_S229, camera_model_5, is_ray_depth_5);
    }
    else
    {
        _S231 = _S228;
    }
    if(_S230)
    {
        _S231 = _S231 * _s_dOut_4;
    }
    else
    {
        _S231 = _S228;
    }
    float _S232 = _S231.x + _S231.y + _S231.z;
    dpdepth_2->primal_0 = (*dpdepth_2).primal_0;
    dpdepth_2->differential_0 = _S232;
    return;
}

inline __device__ void s_bwd_depth_to_point_0(float2  _S233, float4  _S234, FixedArray<float, 10>  * _S235, int _S236, bool _S237, DiffPair_float_0 * _S238, float3  _S239)
{
    s_bwd_prop_depth_to_point_Intermediates_0 _S240;
    float3  _S241 = s_primal_ctx_depth_to_point_0(_S233, _S234, _S235, _S236, _S237, (*_S238).primal_0, &_S240);
    s_bwd_prop_depth_to_point_Intermediates_0 _S242 = _S240;
    s_bwd_prop_depth_to_point_0(_S233, _S234, _S235, _S236, _S237, _S238, _S239, &_S242);
    return;
}

inline __device__ float depth_to_point_vjp(float2  pix_pos_4, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, int camera_model_6, bool is_ray_depth_6, float depth_3, float3  v_point_0)
{
    DiffPair_float_0 dp_depth_0;
    (&dp_depth_0)->primal_0 = depth_3;
    (&dp_depth_0)->differential_0 = 0.0f;
    FixedArray<float, 10>  _S243 = dist_coeffs_5;
    s_bwd_depth_to_point_0(pix_pos_4, intrins_4, &_S243, camera_model_6, is_ray_depth_6, &dp_depth_0, v_point_0);
    return dp_depth_0.differential_0;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_7)
{
    float _S244 = dOut_7.y;
    float _S245 = dOut_7.z;
    float _S246 = dOut_7.x;
    float _S247 = (*a_0).primal_0.z * _S244 + - (*a_0).primal_0.y * _S245;
    float _S248 = - (*a_0).primal_0.z * _S246 + (*a_0).primal_0.x * _S245;
    float _S249 = (*a_0).primal_0.y * _S246 + - (*a_0).primal_0.x * _S244;
    float3  _S250 = make_float3 (- (*b_0).primal_0.z * _S244 + (*b_0).primal_0.y * _S245, (*b_0).primal_0.z * _S246 + - (*b_0).primal_0.x * _S245, - (*b_0).primal_0.y * _S246 + (*b_0).primal_0.x * _S244);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S250;
    float3  _S251 = make_float3 (_S247, _S248, _S249);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S251;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S252 = left_2.y;
    float _S253 = right_2.z;
    float _S254 = left_2.z;
    float _S255 = right_2.y;
    float _S256 = right_2.x;
    float _S257 = left_2.x;
    return make_float3 (_S252 * _S253 - _S254 * _S255, _S254 * _S256 - _S257 * _S253, _S257 * _S255 - _S252 * _S256);
}

inline __device__ float3  points_to_normal(FixedArray<float3 , 4>  points_0)
{
    float3  _S258 = points_0[int(0)];
    bool _S259;
    if((dot_0(_S258, _S258)) == 0.0f)
    {
        _S259 = true;
    }
    else
    {
        float3  _S260 = points_0[int(1)];
        _S259 = (dot_0(_S260, _S260)) == 0.0f;
    }
    if(_S259)
    {
        _S259 = true;
    }
    else
    {
        float3  _S261 = points_0[int(2)];
        _S259 = (dot_0(_S261, _S261)) == 0.0f;
    }
    if(_S259)
    {
        _S259 = true;
    }
    else
    {
        float3  _S262 = points_0[int(3)];
        _S259 = (dot_0(_S262, _S262)) == 0.0f;
    }
    if(_S259)
    {
        return make_float3 (0.0f);
    }
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

inline __device__ float s_primal_ctx_dot_0(float3  _S263, float3  _S264)
{
    return dot_0(_S263, _S264);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S265, float3  _S266)
{
    return cross_0(_S265, _S266);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S267, float _S268)
{
    _d_sqrt_0(_S267, _S268);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_5)
{
    float _S269 = (*dpx_7).primal_0.x;
    float _S270 = (*dpx_7).primal_0.y;
    float _S271 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S272;
    (&_S272)->primal_0 = _S269 * _S269 + _S270 * _S270 + _S271 * _S271;
    (&_S272)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S272, _s_dOut_5);
    float _S273 = (*dpx_7).primal_0.z * _S272.differential_0;
    float _S274 = _S273 + _S273;
    float _S275 = (*dpx_7).primal_0.y * _S272.differential_0;
    float _S276 = _S275 + _S275;
    float _S277 = (*dpx_7).primal_0.x * _S272.differential_0;
    float _S278 = _S277 + _S277;
    float3  _S279 = make_float3 (0.0f);
    *&((&_S279)->z) = _S274;
    *&((&_S279)->y) = _S276;
    *&((&_S279)->x) = _S278;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S279;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S280, float _S281)
{
    s_bwd_prop_length_impl_0(_S280, _S281);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S282, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S283, float _S284)
{
    _d_dot_0(_S282, _S283, _S284);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S285, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S286, float3  _S287)
{
    _d_cross_0(_S285, _S286, _S287);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_6)
{
    FixedArray<float3 , 4>  _S288 = dppoints_0->primal_0;
    float3  _S289 = make_float3 (0.0f);
    float3  _S290 = dppoints_0->primal_0[int(0)];
    bool _S291 = (s_primal_ctx_dot_0(_S290, _S290)) == 0.0f;
    bool _S292;
    float3  _S293;
    if(_S291)
    {
        _S292 = true;
        _S293 = _S289;
    }
    else
    {
        float3  _S294 = _S288[int(1)];
        _S292 = (s_primal_ctx_dot_0(_S294, _S294)) == 0.0f;
        _S293 = _S288[int(1)];
    }
    bool _S295;
    float3  _S296;
    if(_S292)
    {
        _S295 = true;
        _S296 = _S289;
    }
    else
    {
        float3  _S297 = _S288[int(2)];
        _S295 = (s_primal_ctx_dot_0(_S297, _S297)) == 0.0f;
        _S296 = _S288[int(2)];
    }
    bool _S298;
    float3  _S299;
    if(_S295)
    {
        _S298 = true;
        _S299 = _S289;
    }
    else
    {
        float3  _S300 = _S288[int(3)];
        _S298 = (s_primal_ctx_dot_0(_S300, _S300)) == 0.0f;
        _S299 = _S288[int(3)];
    }
    bool _S301 = !_S298;
    float3  _S302;
    float3  _S303;
    float3  _S304;
    float3  _S305;
    float3  _S306;
    if(_S301)
    {
        float3  dx_0 = _S288[int(1)] - _S288[int(0)];
        float3  _S307 = - (_S288[int(3)] - _S288[int(2)]);
        float3  _S308 = s_primal_ctx_cross_0(dx_0, _S307);
        bool _S309 = (s_primal_ctx_dot_0(_S308, _S308)) != 0.0f;
        if(_S309)
        {
            float _S310 = length_1(_S308);
            float3  _S311 = make_float3 (_S310);
            _S302 = make_float3 (_S310 * _S310);
            _S303 = _S311;
        }
        else
        {
            _S302 = _S289;
            _S303 = _S289;
        }
        float3  _S312 = _S303;
        _S298 = _S309;
        _S303 = _S308;
        _S304 = _S312;
        _S305 = dx_0;
        _S306 = _S307;
    }
    else
    {
        _S298 = false;
        _S302 = _S289;
        _S303 = _S289;
        _S304 = _S289;
        _S305 = _S289;
        _S306 = _S289;
    }
    FixedArray<float3 , 4>  _S313;
    if(_S301)
    {
        if(_S298)
        {
            float3  _S314 = _s_dOut_6 / _S302;
            float3  _S315 = _S303 * - _S314;
            float3  _S316 = _S304 * _S314;
            float _S317 = _S315.x + _S315.y + _S315.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S318;
            (&_S318)->primal_0 = _S303;
            (&_S318)->differential_0 = _S289;
            s_bwd_length_impl_0(&_S318, _S317);
            _S302 = _S316 + _S318.differential_0;
        }
        else
        {
            _S302 = _s_dOut_6;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S319;
        (&_S319)->primal_0 = _S303;
        (&_S319)->differential_0 = _S289;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S320;
        (&_S320)->primal_0 = _S303;
        (&_S320)->differential_0 = _S289;
        s_bwd_prop_dot_0(&_S319, &_S320, 0.0f);
        float3  _S321 = _S320.differential_0 + _S319.differential_0 + _S302;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S322;
        (&_S322)->primal_0 = _S305;
        (&_S322)->differential_0 = _S289;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S323;
        (&_S323)->primal_0 = _S306;
        (&_S323)->differential_0 = _S289;
        s_bwd_prop_cross_0(&_S322, &_S323, _S321);
        float3  s_diff_dy_T_0 = - _S323.differential_0;
        float3  _S324 = - s_diff_dy_T_0;
        float3  _S325 = - _S322.differential_0;
        FixedArray<float3 , 4>  _S326;
        _S326[int(0)] = _S289;
        _S326[int(1)] = _S289;
        _S326[int(2)] = _S289;
        _S326[int(3)] = _S289;
        _S326[int(2)] = _S324;
        _S326[int(3)] = s_diff_dy_T_0;
        _S326[int(1)] = _S322.differential_0;
        _S313[int(0)] = _S326[int(0)];
        _S313[int(1)] = _S326[int(1)];
        _S313[int(2)] = _S326[int(2)];
        _S313[int(3)] = _S326[int(3)];
        _S302 = _S325;
    }
    else
    {
        _S313[int(0)] = _S289;
        _S313[int(1)] = _S289;
        _S313[int(2)] = _S289;
        _S313[int(3)] = _S289;
        _S302 = _S289;
    }
    if(_S295)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S327;
        (&_S327)->primal_0 = _S299;
        (&_S327)->differential_0 = _S289;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S328;
        (&_S328)->primal_0 = _S299;
        (&_S328)->differential_0 = _S289;
        s_bwd_prop_dot_0(&_S327, &_S328, 0.0f);
        float3  _S329 = _S328.differential_0 + _S327.differential_0;
        FixedArray<float3 , 4>  _S330;
        _S330[int(0)] = _S289;
        _S330[int(1)] = _S289;
        _S330[int(2)] = _S289;
        _S330[int(3)] = _S289;
        _S330[int(3)] = _S329;
        float3  _S331 = _S313[int(1)] + _S330[int(1)];
        float3  _S332 = _S313[int(2)] + _S330[int(2)];
        float3  _S333 = _S313[int(3)] + _S330[int(3)];
        _S313[int(0)] = _S313[int(0)] + _S330[int(0)];
        _S313[int(1)] = _S331;
        _S313[int(2)] = _S332;
        _S313[int(3)] = _S333;
    }
    if(_S292)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S334;
        (&_S334)->primal_0 = _S296;
        (&_S334)->differential_0 = _S289;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S335;
        (&_S335)->primal_0 = _S296;
        (&_S335)->differential_0 = _S289;
        s_bwd_prop_dot_0(&_S334, &_S335, 0.0f);
        float3  _S336 = _S335.differential_0 + _S334.differential_0;
        FixedArray<float3 , 4>  _S337;
        _S337[int(0)] = _S289;
        _S337[int(1)] = _S289;
        _S337[int(2)] = _S289;
        _S337[int(3)] = _S289;
        _S337[int(2)] = _S336;
        float3  _S338 = _S313[int(1)] + _S337[int(1)];
        float3  _S339 = _S313[int(2)] + _S337[int(2)];
        float3  _S340 = _S313[int(3)] + _S337[int(3)];
        _S313[int(0)] = _S313[int(0)] + _S337[int(0)];
        _S313[int(1)] = _S338;
        _S313[int(2)] = _S339;
        _S313[int(3)] = _S340;
    }
    if(_S291)
    {
    }
    else
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S341;
        (&_S341)->primal_0 = _S293;
        (&_S341)->differential_0 = _S289;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S342;
        (&_S342)->primal_0 = _S293;
        (&_S342)->differential_0 = _S289;
        s_bwd_prop_dot_0(&_S341, &_S342, 0.0f);
        float3  _S343 = _S342.differential_0 + _S341.differential_0;
        FixedArray<float3 , 4>  _S344;
        _S344[int(0)] = _S289;
        _S344[int(1)] = _S289;
        _S344[int(2)] = _S289;
        _S344[int(3)] = _S289;
        _S344[int(1)] = _S343;
        float3  _S345 = _S313[int(1)] + _S344[int(1)];
        float3  _S346 = _S313[int(2)] + _S344[int(2)];
        float3  _S347 = _S313[int(3)] + _S344[int(3)];
        _S313[int(0)] = _S313[int(0)] + _S344[int(0)];
        _S313[int(1)] = _S345;
        _S313[int(2)] = _S346;
        _S313[int(3)] = _S347;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S348;
    (&_S348)->primal_0 = _S288[int(0)];
    (&_S348)->differential_0 = _S289;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S349;
    (&_S349)->primal_0 = _S288[int(0)];
    (&_S349)->differential_0 = _S289;
    s_bwd_prop_dot_0(&_S348, &_S349, 0.0f);
    float3  _S350 = _S349.differential_0 + _S348.differential_0 + _S302;
    FixedArray<float3 , 4>  _S351;
    _S351[int(0)] = _S289;
    _S351[int(1)] = _S289;
    _S351[int(2)] = _S289;
    _S351[int(3)] = _S289;
    _S351[int(0)] = _S350;
    FixedArray<float3 , 4>  _S352 = {
        _S313[int(0)] + _S351[int(0)], _S313[int(1)] + _S351[int(1)], _S313[int(2)] + _S351[int(2)], _S313[int(3)] + _S351[int(3)]
    };
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S352;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S353, float3  _S354)
{
    s_bwd_prop_points_to_normal_0(_S353, _S354);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S355 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S355;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, int camera_model_7, bool is_ray_depth_7, float4  depths_0)
{
    bool _S356;
    if((depths_0.x) == 0.0f)
    {
        _S356 = true;
    }
    else
    {
        _S356 = (depths_0.y) == 0.0f;
    }
    if(_S356)
    {
        _S356 = true;
    }
    else
    {
        _S356 = (depths_0.z) == 0.0f;
    }
    if(_S356)
    {
        _S356 = true;
    }
    else
    {
        _S356 = (depths_0.w) == 0.0f;
    }
    if(_S356)
    {
        return make_float3 (0.0f);
    }
    FixedArray<float3 , 4>  points_2;
    float2  _S357 = float2 {intrins_5.z, intrins_5.w};
    float2  _S358 = float2 {intrins_5.x, intrins_5.y};
    float2  _S359 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S357) / _S358;
    float2  uv_5 = _S359;
    FixedArray<float, 10>  _S360 = dist_coeffs_6;
    bool _S361 = undistort_point_0(_S359, &_S360, int(12), &uv_5);
    if(!_S361)
    {
        return make_float3 (0.0f);
    }
    points_2[int(0)] = make_float3 (depths_0.x) * unproject_raydir_0(uv_5, camera_model_7, is_ray_depth_7);
    float2  _S362 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S357) / _S358;
    float2  uv_6 = _S362;
    FixedArray<float, 10>  _S363 = dist_coeffs_6;
    bool _S364 = undistort_point_0(_S362, &_S363, int(12), &uv_6);
    if(!_S364)
    {
        return make_float3 (0.0f);
    }
    points_2[int(1)] = make_float3 (depths_0.y) * unproject_raydir_0(uv_6, camera_model_7, is_ray_depth_7);
    float2  _S365 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S357) / _S358;
    float2  uv_7 = _S365;
    FixedArray<float, 10>  _S366 = dist_coeffs_6;
    bool _S367 = undistort_point_0(_S365, &_S366, int(12), &uv_7);
    if(!_S367)
    {
        return make_float3 (0.0f);
    }
    points_2[int(2)] = make_float3 (depths_0.z) * unproject_raydir_0(uv_7, camera_model_7, is_ray_depth_7);
    float2  _S368 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S357) / _S358;
    float2  uv_8 = _S368;
    FixedArray<float, 10>  _S369 = dist_coeffs_6;
    bool _S370 = undistort_point_0(_S368, &_S369, int(12), &uv_8);
    if(!_S370)
    {
        return make_float3 (0.0f);
    }
    points_2[int(3)] = make_float3 (depths_0.w) * unproject_raydir_0(uv_8, camera_model_7, is_ray_depth_7);
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
    float2  _S371;
    bool _S372;
    float2  _S373;
    bool _S374;
    float2  _S375;
    bool _S376;
    float2  _S377;
    bool _S378;
};

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_6, FixedArray<float, 10>  * dist_coeffs_7, int camera_model_8, bool is_ray_depth_8, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_2)
{
    float2  _S379 = make_float2 (0.0f);
    _s_diff_ctx_2->_S371 = _S379;
    _s_diff_ctx_2->_S372 = false;
    _s_diff_ctx_2->_S373 = _S379;
    _s_diff_ctx_2->_S374 = false;
    _s_diff_ctx_2->_S375 = _S379;
    _s_diff_ctx_2->_S376 = false;
    _s_diff_ctx_2->_S377 = _S379;
    _s_diff_ctx_2->_S378 = false;
    _s_diff_ctx_2->_S371 = _S379;
    _s_diff_ctx_2->_S372 = false;
    _s_diff_ctx_2->_S373 = _S379;
    _s_diff_ctx_2->_S374 = false;
    _s_diff_ctx_2->_S375 = _S379;
    _s_diff_ctx_2->_S376 = false;
    _s_diff_ctx_2->_S377 = _S379;
    _s_diff_ctx_2->_S378 = false;
    float3  _S380 = make_float3 (0.0f);
    bool _runFlag_0;
    if((dpdepths_0.x) == 0.0f)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (dpdepths_0.y) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (dpdepths_0.z) == 0.0f;
    }
    if(_runFlag_0)
    {
        _runFlag_0 = true;
    }
    else
    {
        _runFlag_0 = (dpdepths_0.w) == 0.0f;
    }
    float3  normal_4;
    if(_runFlag_0)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S381 = !_runFlag_0;
    if(_S381)
    {
        float2  _S382 = float2 {intrins_6.z, intrins_6.w};
        float2  _S383 = float2 {intrins_6.x, intrins_6.y};
        float2  _S384 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S382) / _S383;
        float2  _S385 = _S384;
        bool _S386 = undistort_point_0(_S384, dist_coeffs_7, int(12), &_S385);
        _s_diff_ctx_2->_S371 = _S385;
        _s_diff_ctx_2->_S372 = _S386;
        float2  uv_9 = _S385;
        bool _S387 = !_S386;
        if(_S387)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S388 = !_S387;
        int _S389;
        FixedArray<float3 , 4>  points_3;
        if(_S388)
        {
            float3  _S390 = make_float3 (dpdepths_0.x) * s_primal_ctx_unproject_raydir_0(uv_9, camera_model_8, is_ray_depth_8);
            float2  _S391 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S382) / _S383;
            float2  _S392 = _S391;
            bool _S393 = undistort_point_0(_S391, dist_coeffs_7, int(12), &_S392);
            _s_diff_ctx_2->_S373 = _S392;
            _s_diff_ctx_2->_S374 = _S393;
            float2  uv_10 = _S392;
            bool _S394 = !_S393;
            if(_S394)
            {
                normal_4 = make_float3 (0.0f);
            }
            bool _S395 = !_S394;
            if(_S395)
            {
                float3  _S396 = make_float3 (dpdepths_0.y) * s_primal_ctx_unproject_raydir_0(uv_10, camera_model_8, is_ray_depth_8);
                _S389 = int(2);
                points_3[int(0)] = _S390;
                points_3[int(1)] = _S396;
                points_3[int(2)] = _S380;
                points_3[int(3)] = _S380;
            }
            else
            {
                _S389 = int(0);
                points_3[int(0)] = _S390;
                points_3[int(1)] = _S380;
                points_3[int(2)] = _S380;
                points_3[int(3)] = _S380;
            }
            if(_S389 != int(2))
            {
                _runFlag_0 = false;
            }
            else
            {
                _runFlag_0 = _S388;
                _S389 = int(0);
            }
            if(_runFlag_0)
            {
                float2  _S397 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S382) / _S383;
                float2  _S398 = _S397;
                bool _S399 = undistort_point_0(_S397, dist_coeffs_7, int(12), &_S398);
                _s_diff_ctx_2->_S375 = _S398;
                _s_diff_ctx_2->_S376 = _S399;
                float2  uv_11 = _S398;
                if(!_S399)
                {
                    float3  _S400 = make_float3 (0.0f);
                    _runFlag_0 = false;
                    _S389 = int(0);
                    normal_4 = _S400;
                }
                if(_runFlag_0)
                {
                    points_3[int(2)] = make_float3 (dpdepths_0.z) * s_primal_ctx_unproject_raydir_0(uv_11, camera_model_8, is_ray_depth_8);
                    float2  _S401 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S382) / _S383;
                    float2  _S402 = _S401;
                    bool _S403 = undistort_point_0(_S401, dist_coeffs_7, int(12), &_S402);
                    _s_diff_ctx_2->_S377 = _S402;
                    _s_diff_ctx_2->_S378 = _S403;
                    float2  uv_12 = _S402;
                    bool _S404 = !_S403;
                    if(_S404)
                    {
                        normal_4 = make_float3 (0.0f);
                    }
                    bool _S405 = !_S404;
                    int _S406;
                    if(_S405)
                    {
                        points_3[int(3)] = make_float3 (dpdepths_0.w) * s_primal_ctx_unproject_raydir_0(uv_12, camera_model_8, is_ray_depth_8);
                        _S406 = int(2);
                    }
                    else
                    {
                        _S406 = int(0);
                    }
                    if(_S406 != int(2))
                    {
                        _runFlag_0 = false;
                        _S389 = _S406;
                    }
                    if(_runFlag_0)
                    {
                        _S389 = int(1);
                    }
                }
            }
        }
        else
        {
            _S389 = int(0);
            points_3[int(0)] = _S380;
            points_3[int(1)] = _S380;
            points_3[int(2)] = _S380;
            points_3[int(3)] = _S380;
        }
        if(_S389 != int(1))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S381;
        }
        if(_runFlag_0)
        {
            float3  _S407 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
            if((s_primal_ctx_dot_0(_S407, _S407)) != 0.0f)
            {
                normal_4 = _S407 / make_float3 (length_1(_S407));
            }
            else
            {
                normal_4 = _S407;
            }
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_7, FixedArray<float, 10>  * dist_coeffs_8, int camera_model_9, bool is_ray_depth_9, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_7, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S408 = *dpdepths_1;
    float3  _S409 = make_float3 (0.0f);
    bool _runFlag_1;
    if(((*dpdepths_1).primal_0.x) == 0.0f)
    {
        _runFlag_1 = true;
    }
    else
    {
        _runFlag_1 = (_S408.primal_0.y) == 0.0f;
    }
    if(_runFlag_1)
    {
        _runFlag_1 = true;
    }
    else
    {
        _runFlag_1 = (_S408.primal_0.z) == 0.0f;
    }
    if(_runFlag_1)
    {
        _runFlag_1 = true;
    }
    else
    {
        _runFlag_1 = (_S408.primal_0.w) == 0.0f;
    }
    bool _S410 = !_runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S411;
    bool _runFlag_4;
    bool _S412;
    bool _S413;
    FixedArray<float3 , 4>  points_4;
    float3  _S414;
    float3  _S415;
    float3  _S416;
    float3  _S417;
    float3  _S418;
    float3  _S419;
    float3  _S420;
    float3  _S421;
    float3  _S422;
    if(_S410)
    {
        float2  _S423 = _s_diff_ctx_3->_S371;
        bool _S424 = !!_s_diff_ctx_3->_S372;
        int _S425;
        if(_S424)
        {
            float3  _S426 = s_primal_ctx_unproject_raydir_0(_S423, camera_model_9, is_ray_depth_9);
            float3  _S427 = make_float3 (_S408.primal_0.x) * _S426;
            float2  _S428 = _s_diff_ctx_3->_S373;
            bool _S429 = !!_s_diff_ctx_3->_S374;
            if(_S429)
            {
                float3  _S430 = s_primal_ctx_unproject_raydir_0(_S428, camera_model_9, is_ray_depth_9);
                float3  _S431 = make_float3 (_S408.primal_0.y) * _S430;
                _S425 = int(2);
                points_4[int(0)] = _S427;
                points_4[int(1)] = _S431;
                points_4[int(2)] = _S409;
                points_4[int(3)] = _S409;
                _S414 = _S430;
            }
            else
            {
                _S425 = int(0);
                points_4[int(0)] = _S427;
                points_4[int(1)] = _S409;
                points_4[int(2)] = _S409;
                points_4[int(3)] = _S409;
                _S414 = _S409;
            }
            if(_S425 != int(2))
            {
                _runFlag_1 = false;
            }
            else
            {
                _runFlag_1 = _S424;
                _S425 = int(0);
            }
            if(_runFlag_1)
            {
                float2  _S432 = _s_diff_ctx_3->_S375;
                if(!_s_diff_ctx_3->_S376)
                {
                    _runFlag_2 = false;
                    _S425 = int(0);
                }
                else
                {
                    _runFlag_2 = _runFlag_1;
                }
                if(_runFlag_2)
                {
                    float3  _S433 = s_primal_ctx_unproject_raydir_0(_S432, camera_model_9, is_ray_depth_9);
                    points_4[int(2)] = make_float3 (_S408.primal_0.z) * _S433;
                    float2  _S434 = _s_diff_ctx_3->_S377;
                    bool _S435 = !!_s_diff_ctx_3->_S378;
                    int _S436;
                    if(_S435)
                    {
                        float3  _S437 = s_primal_ctx_unproject_raydir_0(_S434, camera_model_9, is_ray_depth_9);
                        points_4[int(3)] = make_float3 (_S408.primal_0.w) * _S437;
                        _S436 = int(2);
                        _S415 = _S437;
                    }
                    else
                    {
                        _S436 = int(0);
                        _S415 = _S409;
                    }
                    if(_S436 != int(2))
                    {
                        _runFlag_3 = false;
                        _S425 = _S436;
                    }
                    else
                    {
                        _runFlag_3 = _runFlag_2;
                    }
                    if(_runFlag_3)
                    {
                        _S425 = int(1);
                    }
                    _runFlag_3 = _S435;
                    _S416 = _S433;
                }
                else
                {
                    _runFlag_3 = false;
                    _S415 = _S409;
                    _S416 = _S409;
                }
            }
            else
            {
                _runFlag_2 = false;
                _runFlag_3 = false;
                _S415 = _S409;
                _S416 = _S409;
            }
            float3  _S438 = _S414;
            _S414 = _S415;
            _S415 = _S416;
            _S411 = _S429;
            _S416 = _S438;
            _S417 = _S426;
        }
        else
        {
            _S425 = int(0);
            points_4[int(0)] = _S409;
            points_4[int(1)] = _S409;
            points_4[int(2)] = _S409;
            points_4[int(3)] = _S409;
            _runFlag_1 = false;
            _runFlag_2 = false;
            _runFlag_3 = false;
            _S414 = _S409;
            _S415 = _S409;
            _S411 = false;
            _S416 = _S409;
            _S417 = _S409;
        }
        if(_S425 != int(1))
        {
            _runFlag_4 = false;
        }
        else
        {
            _runFlag_4 = _S410;
        }
        if(_runFlag_4)
        {
            float3  dx_1 = points_4[int(1)] - points_4[int(0)];
            float3  _S439 = - (points_4[int(3)] - points_4[int(2)]);
            float3  _S440 = s_primal_ctx_cross_0(dx_1, _S439);
            bool _S441 = (s_primal_ctx_dot_0(_S440, _S440)) != 0.0f;
            if(_S441)
            {
                float _S442 = length_1(_S440);
                float3  _S443 = make_float3 (_S442);
                _S418 = make_float3 (_S442 * _S442);
                _S419 = _S443;
            }
            else
            {
                _S418 = _S409;
                _S419 = _S409;
            }
            float3  _S444 = _S419;
            _S412 = _S441;
            _S419 = _S440;
            _S420 = _S444;
            _S421 = dx_1;
            _S422 = _S439;
        }
        else
        {
            _S412 = false;
            _S418 = _S409;
            _S419 = _S409;
            _S420 = _S409;
            _S421 = _S409;
            _S422 = _S409;
        }
        bool _S445 = _runFlag_1;
        bool _S446 = _runFlag_2;
        bool _S447 = _runFlag_3;
        float3  _S448 = _S414;
        float3  _S449 = _S415;
        bool _S450 = _S411;
        float3  _S451 = _S416;
        float3  _S452 = _S417;
        _runFlag_1 = _runFlag_4;
        _runFlag_2 = _S412;
        _S414 = _S418;
        _S415 = _S419;
        _S416 = _S420;
        _S417 = _S421;
        _S418 = _S422;
        _runFlag_3 = _S424;
        _S411 = _S445;
        _runFlag_4 = _S446;
        _S412 = _S447;
        _S419 = _S448;
        _S420 = _S449;
        _S413 = _S450;
        _S421 = _S451;
        _S422 = _S452;
    }
    else
    {
        _runFlag_1 = false;
        _runFlag_2 = false;
        _S414 = _S409;
        _S415 = _S409;
        _S416 = _S409;
        _S417 = _S409;
        _S418 = _S409;
        _runFlag_3 = false;
        _S411 = false;
        _runFlag_4 = false;
        _S412 = false;
        _S419 = _S409;
        _S420 = _S409;
        _S413 = false;
        _S421 = _S409;
        _S422 = _S409;
    }
    float4  _S453 = make_float4 (0.0f);
    float4  _S454;
    if(_S410)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                float3  _S455 = _s_dOut_7 / _S414;
                float3  _S456 = _S415 * - _S455;
                float3  _S457 = _S416 * _S455;
                float _S458 = _S456.x + _S456.y + _S456.z;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S459;
                (&_S459)->primal_0 = _S415;
                (&_S459)->differential_0 = _S409;
                s_bwd_length_impl_0(&_S459, _S458);
                _S414 = _S457 + _S459.differential_0;
            }
            else
            {
                _S414 = _s_dOut_7;
            }
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S460;
            (&_S460)->primal_0 = _S415;
            (&_S460)->differential_0 = _S409;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S461;
            (&_S461)->primal_0 = _S415;
            (&_S461)->differential_0 = _S409;
            s_bwd_prop_dot_0(&_S460, &_S461, 0.0f);
            float3  _S462 = _S461.differential_0 + _S460.differential_0 + _S414;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S463;
            (&_S463)->primal_0 = _S417;
            (&_S463)->differential_0 = _S409;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S464;
            (&_S464)->primal_0 = _S418;
            (&_S464)->differential_0 = _S409;
            s_bwd_prop_cross_0(&_S463, &_S464, _S462);
            float3  s_diff_dy_T_1 = - _S464.differential_0;
            float3  _S465 = - s_diff_dy_T_1;
            float3  _S466 = - _S463.differential_0;
            FixedArray<float3 , 4>  _S467;
            _S467[int(0)] = _S409;
            _S467[int(1)] = _S409;
            _S467[int(2)] = _S409;
            _S467[int(3)] = _S409;
            _S467[int(2)] = _S465;
            _S467[int(3)] = s_diff_dy_T_1;
            _S467[int(0)] = _S466;
            _S467[int(1)] = _S463.differential_0;
            points_4[int(0)] = _S467[int(0)];
            points_4[int(1)] = _S467[int(1)];
            points_4[int(2)] = _S467[int(2)];
            points_4[int(3)] = _S467[int(3)];
        }
        else
        {
            points_4[int(0)] = _S409;
            points_4[int(1)] = _S409;
            points_4[int(2)] = _S409;
            points_4[int(3)] = _S409;
        }
        if(_runFlag_3)
        {
            if(_S411)
            {
                if(_runFlag_4)
                {
                    FixedArray<float3 , 4>  _S468 = points_4;
                    FixedArray<float3 , 4>  _S469 = points_4;
                    FixedArray<float3 , 4>  _S470 = points_4;
                    FixedArray<float3 , 4>  _S471 = points_4;
                    if(_S412)
                    {
                        float3  _S472 = _S419 * _S471[int(3)];
                        float _S473 = _S472.x + _S472.y + _S472.z;
                        float4  _S474 = _S453;
                        *&((&_S474)->w) = _S473;
                        points_4[int(0)] = _S468[int(0)];
                        points_4[int(1)] = _S469[int(1)];
                        points_4[int(2)] = _S470[int(2)];
                        points_4[int(3)] = _S409;
                        _S454 = _S474;
                    }
                    else
                    {
                        points_4[int(0)] = _S468[int(0)];
                        points_4[int(1)] = _S469[int(1)];
                        points_4[int(2)] = _S470[int(2)];
                        points_4[int(3)] = _S471[int(3)];
                        _S454 = _S453;
                    }
                    float3  _S475 = _S420 * points_4[int(2)];
                    float _S476 = _S475.x + _S475.y + _S475.z;
                    FixedArray<float3 , 4>  _S477 = points_4;
                    FixedArray<float3 , 4>  _S478 = points_4;
                    float4  _S479 = _S453;
                    *&((&_S479)->z) = _S476;
                    float4  _S480 = _S454 + _S479;
                    points_4[int(0)] = points_4[int(0)];
                    points_4[int(1)] = _S477[int(1)];
                    points_4[int(2)] = _S409;
                    points_4[int(3)] = _S478[int(3)];
                    _S454 = _S480;
                }
                else
                {
                    FixedArray<float3 , 4>  _S481 = points_4;
                    FixedArray<float3 , 4>  _S482 = points_4;
                    FixedArray<float3 , 4>  _S483 = points_4;
                    points_4[int(0)] = points_4[int(0)];
                    points_4[int(1)] = _S481[int(1)];
                    points_4[int(2)] = _S482[int(2)];
                    points_4[int(3)] = _S483[int(3)];
                    _S454 = _S453;
                }
            }
            else
            {
                FixedArray<float3 , 4>  _S484 = points_4;
                FixedArray<float3 , 4>  _S485 = points_4;
                FixedArray<float3 , 4>  _S486 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S484[int(1)];
                points_4[int(2)] = _S485[int(2)];
                points_4[int(3)] = _S486[int(3)];
                _S454 = _S453;
            }
            if(_S413)
            {
                FixedArray<float3 , 4>  _S487 = points_4;
                float3  _S488 = _S421 * points_4[int(1)];
                float _S489 = _S488.x + _S488.y + _S488.z;
                float4  _S490 = _S453;
                *&((&_S490)->y) = _S489;
                float4  _S491 = _S454 + _S490;
                points_4[int(0)] = _S409;
                points_4[int(1)] = _S409;
                points_4[int(2)] = _S409;
                points_4[int(3)] = _S409;
                _S414 = _S487[int(0)];
                _S454 = _S491;
            }
            else
            {
                FixedArray<float3 , 4>  _S492 = points_4;
                FixedArray<float3 , 4>  _S493 = points_4;
                FixedArray<float3 , 4>  _S494 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S492[int(1)];
                points_4[int(2)] = _S493[int(2)];
                points_4[int(3)] = _S494[int(3)];
                _S414 = _S409;
            }
            float3  _S495 = _S422 * (points_4[int(0)] + _S414);
            float _S496 = _S495.x + _S495.y + _S495.z;
            float4  _S497 = _S453;
            *&((&_S497)->x) = _S496;
            _S454 = _S454 + _S497;
        }
        else
        {
            _S454 = _S453;
        }
    }
    else
    {
        _S454 = _S453;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S454;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S498, float4  _S499, FixedArray<float, 10>  * _S500, int _S501, bool _S502, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S503, float3  _S504)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S505;
    float3  _S506 = s_primal_ctx_depth_to_normal_0(_S498, _S499, _S500, _S501, _S502, (*_S503).primal_0, &_S505);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S507 = _S505;
    s_bwd_prop_depth_to_normal_0(_S498, _S499, _S500, _S501, _S502, _S503, _S504, &_S507);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_8, FixedArray<float, 10>  dist_coeffs_9, int camera_model_10, bool is_ray_depth_10, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S508 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S508;
    FixedArray<float, 10>  _S509 = dist_coeffs_9;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_8, &_S509, camera_model_10, is_ray_depth_10, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_9, FixedArray<float, 10>  dist_coeffs_10, int camera_model_11)
{
    float2  _S510 = (pix_center_4 - float2 {intrins_9.z, intrins_9.w}) / float2 {intrins_9.x, intrins_9.y};
    float2  uv_13 = _S510;
    FixedArray<float, 10>  _S511 = dist_coeffs_10;
    bool _S512 = undistort_point_0(_S510, &_S511, int(12), &uv_13);
    if(!_S512)
    {
        return 0.0f;
    }
    float3  raydir_2 = unproject_raydir_0(uv_13, camera_model_11, false);
    return float((F32_sign((raydir_2.z)))) / length_1(raydir_2);
}

inline __device__ float depth_normal_loss(float2  pix_center_5, float4  intrins_10, FixedArray<float, 10>  dist_coeffs_11, int camera_model_12, bool is_ray_depth_11, float4  depths_2, float3  gt_normal_0)
{
    FixedArray<float3 , 5>  points_5;
    float2  _S513 = float2 {intrins_10.z, intrins_10.w};
    float2  _S514 = float2 {intrins_10.x, intrins_10.y};
    float2  _S515 = (pix_center_5 + make_float2 (-1.0f, -0.0f) - _S513) / _S514;
    float2  uv_14 = _S515;
    FixedArray<float, 10>  _S516 = dist_coeffs_11;
    bool _S517 = undistort_point_0(_S515, &_S516, int(12), &uv_14);
    if(!_S517)
    {
        return 0.0f;
    }
    float3  raydir_3 = unproject_raydir_0(uv_14, camera_model_12, is_ray_depth_11);
    points_5[int(0)] = make_float3 (depths_2.x) * raydir_3;
    float2  _S518 = (pix_center_5 + make_float2 (1.0f, -0.0f) - _S513) / _S514;
    float2  uv_15 = _S518;
    FixedArray<float, 10>  _S519 = dist_coeffs_11;
    bool _S520 = undistort_point_0(_S518, &_S519, int(12), &uv_15);
    if(!_S520)
    {
        return 0.0f;
    }
    float3  raydir_4 = unproject_raydir_0(uv_15, camera_model_12, is_ray_depth_11);
    points_5[int(1)] = make_float3 (depths_2.y) * raydir_4;
    float2  _S521 = (pix_center_5 + make_float2 (0.0f, -1.0f) - _S513) / _S514;
    float2  uv_16 = _S521;
    FixedArray<float, 10>  _S522 = dist_coeffs_11;
    bool _S523 = undistort_point_0(_S521, &_S522, int(12), &uv_16);
    if(!_S523)
    {
        return 0.0f;
    }
    float3  raydir_5 = unproject_raydir_0(uv_16, camera_model_12, is_ray_depth_11);
    points_5[int(2)] = make_float3 (depths_2.z) * raydir_5;
    float2  _S524 = (pix_center_5 + make_float2 (0.0f, 1.0f) - _S513) / _S514;
    float2  uv_17 = _S524;
    FixedArray<float, 10>  _S525 = dist_coeffs_11;
    bool _S526 = undistort_point_0(_S524, &_S525, int(12), &uv_17);
    if(!_S526)
    {
        return 0.0f;
    }
    float3  raydir_6 = unproject_raydir_0(uv_17, camera_model_12, is_ray_depth_11);
    points_5[int(3)] = make_float3 (depths_2.w) * raydir_6;
    float2  _S527 = (pix_center_5 + make_float2 (0.0f) * make_float2 (0.0f, 3.0f) - _S513) / _S514;
    float2  uv_18 = _S527;
    FixedArray<float, 10>  _S528 = dist_coeffs_11;
    bool _S529 = undistort_point_0(_S527, &_S528, int(12), &uv_18);
    if(!_S529)
    {
        return 0.0f;
    }
    float3  raydir_7 = unproject_raydir_0(uv_18, camera_model_12, is_ray_depth_11);
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
    float3  _S530;
    if((dot_0(gt_normal_0, gt_normal_0)) != 0.0f)
    {
        _S530 = normalize_0(gt_normal_0);
    }
    else
    {
        _S530 = gt_normal_0;
    }
    return (1.0f - dot_0(normal_6, _S530) + 0.00100000004749745f) / ((F32_max((dot_0(normal_6, - normalize_0(raydir_7))), (0.0f))) + 0.00100000004749745f);
}

struct s_bwd_prop_depth_normal_loss_Intermediates_0
{
    float2  _S531;
    bool _S532;
    float2  _S533;
    bool _S534;
    float2  _S535;
    bool _S536;
    float2  _S537;
    bool _S538;
    float2  _S539;
    bool _S540;
};

inline __device__ float s_primal_ctx_depth_normal_loss_0(float2  pix_center_6, float4  intrins_11, FixedArray<float, 10>  * dist_coeffs_12, int camera_model_13, bool is_ray_depth_12, float4  dpdepths_2, float3  dpgt_normal_0, s_bwd_prop_depth_normal_loss_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S541 = make_float2 (0.0f);
    _s_diff_ctx_4->_S531 = _S541;
    _s_diff_ctx_4->_S532 = false;
    _s_diff_ctx_4->_S533 = _S541;
    _s_diff_ctx_4->_S534 = false;
    _s_diff_ctx_4->_S535 = _S541;
    _s_diff_ctx_4->_S536 = false;
    _s_diff_ctx_4->_S537 = _S541;
    _s_diff_ctx_4->_S538 = false;
    _s_diff_ctx_4->_S539 = _S541;
    _s_diff_ctx_4->_S540 = false;
    _s_diff_ctx_4->_S533 = _S541;
    _s_diff_ctx_4->_S534 = false;
    _s_diff_ctx_4->_S535 = _S541;
    _s_diff_ctx_4->_S536 = false;
    _s_diff_ctx_4->_S537 = _S541;
    _s_diff_ctx_4->_S538 = false;
    _s_diff_ctx_4->_S539 = _S541;
    _s_diff_ctx_4->_S540 = false;
    float3  _S542 = make_float3 (0.0f);
    float2  _S543 = float2 {intrins_11.z, intrins_11.w};
    float2  _S544 = float2 {intrins_11.x, intrins_11.y};
    float2  _S545 = (pix_center_6 + make_float2 (-1.0f, -0.0f) - _S543) / _S544;
    float2  _S546 = _S545;
    bool _S547 = undistort_point_0(_S545, dist_coeffs_12, int(12), &_S546);
    _s_diff_ctx_4->_S531 = _S546;
    _s_diff_ctx_4->_S532 = _S547;
    float2  uv_19 = _S546;
    bool _S548 = !!_S547;
    int _S549;
    float3  raydir_8;
    FixedArray<float3 , 5>  points_6;
    if(_S548)
    {
        float3  _S550 = s_primal_ctx_unproject_raydir_0(uv_19, camera_model_13, is_ray_depth_12);
        float3  _S551 = make_float3 (dpdepths_2.x) * _S550;
        float2  _S552 = (pix_center_6 + make_float2 (1.0f, -0.0f) - _S543) / _S544;
        float2  _S553 = _S552;
        bool _S554 = undistort_point_0(_S552, dist_coeffs_12, int(12), &_S553);
        _s_diff_ctx_4->_S533 = _S553;
        _s_diff_ctx_4->_S534 = _S554;
        float2  uv_20 = _S553;
        bool _runFlag_5;
        if(!_S554)
        {
            _runFlag_5 = false;
        }
        else
        {
            _runFlag_5 = _S548;
        }
        if(_runFlag_5)
        {
            float3  _S555 = s_primal_ctx_unproject_raydir_0(uv_20, camera_model_13, is_ray_depth_12);
            float3  _S556 = make_float3 (dpdepths_2.y) * _S555;
            float2  _S557 = (pix_center_6 + make_float2 (0.0f, -1.0f) - _S543) / _S544;
            float2  _S558 = _S557;
            bool _S559 = undistort_point_0(_S557, dist_coeffs_12, int(12), &_S558);
            _s_diff_ctx_4->_S535 = _S558;
            _s_diff_ctx_4->_S536 = _S559;
            float2  uv_21 = _S558;
            if(!_S559)
            {
                _runFlag_5 = false;
            }
            if(_runFlag_5)
            {
                float3  _S560 = s_primal_ctx_unproject_raydir_0(uv_21, camera_model_13, is_ray_depth_12);
                float3  _S561 = make_float3 (dpdepths_2.z) * _S560;
                float2  _S562 = (pix_center_6 + make_float2 (0.0f, 1.0f) - _S543) / _S544;
                float2  _S563 = _S562;
                bool _S564 = undistort_point_0(_S562, dist_coeffs_12, int(12), &_S563);
                _s_diff_ctx_4->_S537 = _S563;
                _s_diff_ctx_4->_S538 = _S564;
                float2  uv_22 = _S563;
                if(!_S564)
                {
                    _runFlag_5 = false;
                }
                if(_runFlag_5)
                {
                    float3  _S565 = s_primal_ctx_unproject_raydir_0(uv_22, camera_model_13, is_ray_depth_12);
                    float3  _S566 = make_float3 (dpdepths_2.w) * _S565;
                    float2  _S567 = (pix_center_6 - _S543) / _S544;
                    float2  _S568 = _S567;
                    bool _S569 = undistort_point_0(_S567, dist_coeffs_12, int(12), &_S568);
                    _s_diff_ctx_4->_S539 = _S568;
                    _s_diff_ctx_4->_S540 = _S569;
                    float2  uv_23 = _S568;
                    if(!_S569)
                    {
                        _runFlag_5 = false;
                    }
                    if(_runFlag_5)
                    {
                        float3  _S570 = s_primal_ctx_unproject_raydir_0(uv_23, camera_model_13, is_ray_depth_12);
                        _S549 = int(1);
                        raydir_8 = _S570;
                    }
                    else
                    {
                        _S549 = int(0);
                        raydir_8 = _S565;
                    }
                    points_6[int(0)] = _S551;
                    points_6[int(1)] = _S556;
                    points_6[int(2)] = _S561;
                    points_6[int(3)] = _S566;
                    points_6[int(4)] = _S542;
                }
                else
                {
                    _S549 = int(0);
                    raydir_8 = _S560;
                    points_6[int(0)] = _S551;
                    points_6[int(1)] = _S556;
                    points_6[int(2)] = _S561;
                    points_6[int(3)] = _S542;
                    points_6[int(4)] = _S542;
                }
            }
            else
            {
                _S549 = int(0);
                raydir_8 = _S555;
                points_6[int(0)] = _S551;
                points_6[int(1)] = _S556;
                points_6[int(2)] = _S542;
                points_6[int(3)] = _S542;
                points_6[int(4)] = _S542;
            }
        }
        else
        {
            _S549 = int(0);
            raydir_8 = _S550;
            points_6[int(0)] = _S551;
            points_6[int(1)] = _S542;
            points_6[int(2)] = _S542;
            points_6[int(3)] = _S542;
            points_6[int(4)] = _S542;
        }
    }
    else
    {
        _S549 = int(0);
        points_6[int(0)] = _S542;
        points_6[int(1)] = _S542;
        points_6[int(2)] = _S542;
        points_6[int(3)] = _S542;
        points_6[int(4)] = _S542;
    }
    float _S571;
    if(!(_S549 != int(1)))
    {
        float3  _S572 = s_primal_ctx_cross_0(points_6[int(1)] - points_6[int(0)], - (points_6[int(3)] - points_6[int(2)]));
        float3  normal_7;
        if((s_primal_ctx_dot_0(_S572, _S572)) != 0.0f)
        {
            normal_7 = normalize_0(_S572);
        }
        else
        {
            normal_7 = _S572;
        }
        float3  _S573;
        if((s_primal_ctx_dot_0(dpgt_normal_0, dpgt_normal_0)) != 0.0f)
        {
            _S573 = normalize_0(dpgt_normal_0);
        }
        else
        {
            _S573 = dpgt_normal_0;
        }
        _S571 = (1.0f - s_primal_ctx_dot_0(normal_7, _S573) + 0.00100000004749745f) / ((F32_max((s_primal_ctx_dot_0(normal_7, - normalize_0(raydir_8))), (0.0f))) + 0.00100000004749745f);
    }
    else
    {
        _S571 = 0.0f;
    }
    return _S571;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_8, float3  _s_dOut_8)
{
    float _S574 = length_1((*dpx_8).primal_0);
    float3  _S575 = (*dpx_8).primal_0 * _s_dOut_8;
    float3  _S576 = make_float3 (1.0f / _S574) * _s_dOut_8;
    float _S577 = - ((_S575.x + _S575.y + _S575.z) / (_S574 * _S574));
    float3  _S578 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S579;
    (&_S579)->primal_0 = (*dpx_8).primal_0;
    (&_S579)->differential_0 = _S578;
    s_bwd_length_impl_0(&_S579, _S577);
    float3  _S580 = _S576 + _S579.differential_0;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S580;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S581, float3  _S582)
{
    s_bwd_prop_normalize_impl_0(_S581, _S582);
    return;
}

inline __device__ void s_bwd_prop_depth_normal_loss_0(float2  pix_center_7, float4  intrins_12, FixedArray<float, 10>  * dist_coeffs_13, int camera_model_14, bool is_ray_depth_13, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpgt_normal_1, float _s_dOut_9, s_bwd_prop_depth_normal_loss_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S583 = *dpdepths_3;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S584 = *dpgt_normal_1;
    float3  _S585 = make_float3 (0.0f);
    float2  _S586 = _s_diff_ctx_5->_S531;
    bool _S587 = !!_s_diff_ctx_5->_S532;
    bool _runFlag_6;
    bool _runFlag_7;
    bool _runFlag_8;
    bool _runFlag_9;
    int _S588;
    float3  raydir_9;
    float3  _S589;
    float3  _S590;
    float3  _S591;
    float3  _S592;
    FixedArray<float3 , 5>  points_7;
    if(_S587)
    {
        float3  _S593 = s_primal_ctx_unproject_raydir_0(_S586, camera_model_14, is_ray_depth_13);
        float3  _S594 = make_float3 (_S583.primal_0.x) * _S593;
        float2  _S595 = _s_diff_ctx_5->_S533;
        if(!_s_diff_ctx_5->_S534)
        {
            _runFlag_6 = false;
        }
        else
        {
            _runFlag_6 = _S587;
        }
        if(_runFlag_6)
        {
            float3  _S596 = s_primal_ctx_unproject_raydir_0(_S595, camera_model_14, is_ray_depth_13);
            float3  _S597 = make_float3 (_S583.primal_0.y) * _S596;
            float2  _S598 = _s_diff_ctx_5->_S535;
            if(!_s_diff_ctx_5->_S536)
            {
                _runFlag_7 = false;
            }
            else
            {
                _runFlag_7 = _runFlag_6;
            }
            if(_runFlag_7)
            {
                float3  _S599 = s_primal_ctx_unproject_raydir_0(_S598, camera_model_14, is_ray_depth_13);
                float3  _S600 = make_float3 (_S583.primal_0.z) * _S599;
                float2  _S601 = _s_diff_ctx_5->_S537;
                if(!_s_diff_ctx_5->_S538)
                {
                    _runFlag_8 = false;
                }
                else
                {
                    _runFlag_8 = _runFlag_7;
                }
                if(_runFlag_8)
                {
                    float3  _S602 = s_primal_ctx_unproject_raydir_0(_S601, camera_model_14, is_ray_depth_13);
                    float3  _S603 = make_float3 (_S583.primal_0.w) * _S602;
                    float2  _S604 = _s_diff_ctx_5->_S539;
                    if(!_s_diff_ctx_5->_S540)
                    {
                        _runFlag_9 = false;
                    }
                    else
                    {
                        _runFlag_9 = _runFlag_8;
                    }
                    if(_runFlag_9)
                    {
                        float3  _S605 = s_primal_ctx_unproject_raydir_0(_S604, camera_model_14, is_ray_depth_13);
                        _S588 = int(1);
                        raydir_9 = _S605;
                    }
                    else
                    {
                        _S588 = int(0);
                        raydir_9 = _S602;
                    }
                    points_7[int(0)] = _S594;
                    points_7[int(1)] = _S597;
                    points_7[int(2)] = _S600;
                    points_7[int(3)] = _S603;
                    points_7[int(4)] = _S585;
                    _S589 = _S602;
                }
                else
                {
                    _S588 = int(0);
                    raydir_9 = _S599;
                    points_7[int(0)] = _S594;
                    points_7[int(1)] = _S597;
                    points_7[int(2)] = _S600;
                    points_7[int(3)] = _S585;
                    points_7[int(4)] = _S585;
                    _S589 = _S585;
                }
                _S590 = _S599;
            }
            else
            {
                _S588 = int(0);
                raydir_9 = _S596;
                points_7[int(0)] = _S594;
                points_7[int(1)] = _S597;
                points_7[int(2)] = _S585;
                points_7[int(3)] = _S585;
                points_7[int(4)] = _S585;
                _runFlag_8 = false;
                _S589 = _S585;
                _S590 = _S585;
            }
            _S591 = _S596;
        }
        else
        {
            _S588 = int(0);
            raydir_9 = _S593;
            points_7[int(0)] = _S594;
            points_7[int(1)] = _S585;
            points_7[int(2)] = _S585;
            points_7[int(3)] = _S585;
            points_7[int(4)] = _S585;
            _runFlag_7 = false;
            _runFlag_8 = false;
            _S589 = _S585;
            _S590 = _S585;
            _S591 = _S585;
        }
        _S592 = _S593;
    }
    else
    {
        _S588 = int(0);
        points_7[int(0)] = _S585;
        points_7[int(1)] = _S585;
        points_7[int(2)] = _S585;
        points_7[int(3)] = _S585;
        points_7[int(4)] = _S585;
        _runFlag_6 = false;
        _runFlag_7 = false;
        _runFlag_8 = false;
        _S589 = _S585;
        _S590 = _S585;
        _S591 = _S585;
        _S592 = _S585;
    }
    bool _S606 = !(_S588 != int(1));
    bool _S607;
    float3  normal_8;
    float3  _S608;
    float3  _S609;
    float3  _S610;
    float3  _S611;
    float _S612;
    float _S613;
    float _S614;
    float _S615;
    if(_S606)
    {
        float3  dx_2 = points_7[int(1)] - points_7[int(0)];
        float3  _S616 = - (points_7[int(3)] - points_7[int(2)]);
        float3  _S617 = s_primal_ctx_cross_0(dx_2, _S616);
        bool _S618 = (s_primal_ctx_dot_0(_S617, _S617)) != 0.0f;
        if(_S618)
        {
            normal_8 = normalize_0(_S617);
        }
        else
        {
            normal_8 = _S617;
        }
        bool _S619 = (s_primal_ctx_dot_0(_S584.primal_0, _S584.primal_0)) != 0.0f;
        if(_S619)
        {
            _S608 = normalize_0(_S584.primal_0);
        }
        else
        {
            _S608 = _S584.primal_0;
        }
        float3  _S620 = - normalize_0(raydir_9);
        float _S621 = s_primal_ctx_dot_0(normal_8, _S620);
        float _S622 = 1.0f - s_primal_ctx_dot_0(normal_8, _S608) + 0.00100000004749745f;
        float _S623 = (F32_max((_S621), (0.0f))) + 0.00100000004749745f;
        _S612 = _S623 * _S623;
        _S613 = _S622;
        _S614 = _S623;
        _S615 = _S621;
        raydir_9 = normal_8;
        normal_8 = _S620;
        _runFlag_9 = _S619;
        _S607 = _S618;
        _S609 = _S617;
        _S610 = dx_2;
        _S611 = _S616;
    }
    else
    {
        _S612 = 0.0f;
        _S613 = 0.0f;
        _S614 = 0.0f;
        _S615 = 0.0f;
        raydir_9 = _S585;
        normal_8 = _S585;
        _S608 = _S585;
        _runFlag_9 = false;
        _S607 = false;
        _S609 = _S585;
        _S610 = _S585;
        _S611 = _S585;
    }
    float4  _S624 = make_float4 (0.0f);
    if(_S606)
    {
        float _S625 = _s_dOut_9 / _S612;
        float _S626 = _S613 * - _S625;
        float s_diff_num_T_0 = _S614 * _S625;
        DiffPair_float_0 _S627;
        (&_S627)->primal_0 = _S615;
        (&_S627)->differential_0 = 0.0f;
        DiffPair_float_0 _S628;
        (&_S628)->primal_0 = 0.0f;
        (&_S628)->differential_0 = 0.0f;
        _d_max_0(&_S627, &_S628, _S626);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S629;
        (&_S629)->primal_0 = raydir_9;
        (&_S629)->differential_0 = _S585;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S630;
        (&_S630)->primal_0 = normal_8;
        (&_S630)->differential_0 = _S585;
        s_bwd_prop_dot_0(&_S629, &_S630, _S627.differential_0);
        float _S631 = - s_diff_num_T_0;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S632;
        (&_S632)->primal_0 = raydir_9;
        (&_S632)->differential_0 = _S585;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S633;
        (&_S633)->primal_0 = _S608;
        (&_S633)->differential_0 = _S585;
        s_bwd_prop_dot_0(&_S632, &_S633, _S631);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S634 = _S633;
        float3  _S635 = _S629.differential_0 + _S632.differential_0;
        if(_runFlag_9)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S636;
            (&_S636)->primal_0 = _S584.primal_0;
            (&_S636)->differential_0 = _S585;
            s_bwd_normalize_impl_0(&_S636, _S634.differential_0);
            raydir_9 = _S636.differential_0;
        }
        else
        {
            raydir_9 = _S634.differential_0;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S637;
        (&_S637)->primal_0 = _S584.primal_0;
        (&_S637)->differential_0 = _S585;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S638;
        (&_S638)->primal_0 = _S584.primal_0;
        (&_S638)->differential_0 = _S585;
        s_bwd_prop_dot_0(&_S637, &_S638, 0.0f);
        float3  _S639 = _S638.differential_0 + _S637.differential_0 + raydir_9;
        if(_S607)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S640;
            (&_S640)->primal_0 = _S609;
            (&_S640)->differential_0 = _S585;
            s_bwd_normalize_impl_0(&_S640, _S635);
            raydir_9 = _S640.differential_0;
        }
        else
        {
            raydir_9 = _S635;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S641;
        (&_S641)->primal_0 = _S609;
        (&_S641)->differential_0 = _S585;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S642;
        (&_S642)->primal_0 = _S609;
        (&_S642)->differential_0 = _S585;
        s_bwd_prop_dot_0(&_S641, &_S642, 0.0f);
        float3  _S643 = _S642.differential_0 + _S641.differential_0 + raydir_9;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S644;
        (&_S644)->primal_0 = _S610;
        (&_S644)->differential_0 = _S585;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S645;
        (&_S645)->primal_0 = _S611;
        (&_S645)->differential_0 = _S585;
        s_bwd_prop_cross_0(&_S644, &_S645, _S643);
        float3  s_diff_dy_T_2 = - _S645.differential_0;
        float3  _S646 = - s_diff_dy_T_2;
        float3  _S647 = - _S644.differential_0;
        FixedArray<float3 , 5>  _S648;
        _S648[int(0)] = _S585;
        _S648[int(1)] = _S585;
        _S648[int(2)] = _S585;
        _S648[int(3)] = _S585;
        _S648[int(4)] = _S585;
        _S648[int(2)] = _S646;
        _S648[int(3)] = s_diff_dy_T_2;
        _S648[int(0)] = _S647;
        _S648[int(1)] = _S644.differential_0;
        points_7[int(0)] = _S648[int(0)];
        points_7[int(1)] = _S648[int(1)];
        points_7[int(2)] = _S648[int(2)];
        points_7[int(3)] = _S648[int(3)];
        points_7[int(4)] = _S648[int(4)];
        raydir_9 = _S639;
    }
    else
    {
        points_7[int(0)] = _S585;
        points_7[int(1)] = _S585;
        points_7[int(2)] = _S585;
        points_7[int(3)] = _S585;
        points_7[int(4)] = _S585;
        raydir_9 = _S585;
    }
    float4  _S649;
    if(_S587)
    {
        if(_runFlag_6)
        {
            if(_runFlag_7)
            {
                if(_runFlag_8)
                {
                    FixedArray<float3 , 5>  _S650 = points_7;
                    FixedArray<float3 , 5>  _S651 = points_7;
                    FixedArray<float3 , 5>  _S652 = points_7;
                    float3  _S653 = _S589 * points_7[int(3)];
                    float _S654 = _S653.x + _S653.y + _S653.z;
                    float4  _S655 = _S624;
                    *&((&_S655)->w) = _S654;
                    points_7[int(0)] = _S585;
                    points_7[int(1)] = _S585;
                    points_7[int(2)] = _S585;
                    points_7[int(3)] = _S585;
                    points_7[int(4)] = _S585;
                    _S589 = _S652[int(2)];
                    normal_8 = _S650[int(0)];
                    _S608 = _S651[int(1)];
                    _S649 = _S655;
                }
                else
                {
                    FixedArray<float3 , 5>  _S656 = points_7;
                    FixedArray<float3 , 5>  _S657 = points_7;
                    FixedArray<float3 , 5>  _S658 = points_7;
                    FixedArray<float3 , 5>  _S659 = points_7;
                    points_7[int(0)] = points_7[int(0)];
                    points_7[int(1)] = _S656[int(1)];
                    points_7[int(2)] = _S657[int(2)];
                    points_7[int(3)] = _S658[int(3)];
                    points_7[int(4)] = _S659[int(4)];
                    _S589 = _S585;
                    normal_8 = _S585;
                    _S608 = _S585;
                    _S649 = _S624;
                }
                float3  _S660 = _S590 * (points_7[int(2)] + _S589);
                float _S661 = _S660.x + _S660.y + _S660.z;
                float3  _S662 = points_7[int(0)] + normal_8;
                float3  _S663 = points_7[int(1)] + _S608;
                float4  _S664 = _S624;
                *&((&_S664)->z) = _S661;
                float4  _S665 = _S649 + _S664;
                points_7[int(0)] = _S585;
                points_7[int(1)] = _S585;
                points_7[int(2)] = _S585;
                points_7[int(3)] = _S585;
                points_7[int(4)] = _S585;
                _S589 = _S663;
                _S590 = _S662;
                _S649 = _S665;
            }
            else
            {
                FixedArray<float3 , 5>  _S666 = points_7;
                FixedArray<float3 , 5>  _S667 = points_7;
                FixedArray<float3 , 5>  _S668 = points_7;
                FixedArray<float3 , 5>  _S669 = points_7;
                points_7[int(0)] = points_7[int(0)];
                points_7[int(1)] = _S666[int(1)];
                points_7[int(2)] = _S667[int(2)];
                points_7[int(3)] = _S668[int(3)];
                points_7[int(4)] = _S669[int(4)];
                _S589 = _S585;
                _S590 = _S585;
                _S649 = _S624;
            }
            float3  _S670 = _S591 * (points_7[int(1)] + _S589);
            float _S671 = _S670.x + _S670.y + _S670.z;
            float3  _S672 = points_7[int(0)] + _S590;
            float4  _S673 = _S624;
            *&((&_S673)->y) = _S671;
            float4  _S674 = _S649 + _S673;
            points_7[int(0)] = _S585;
            points_7[int(1)] = _S585;
            points_7[int(2)] = _S585;
            points_7[int(3)] = _S585;
            points_7[int(4)] = _S585;
            _S589 = _S672;
            _S649 = _S674;
        }
        else
        {
            FixedArray<float3 , 5>  _S675 = points_7;
            FixedArray<float3 , 5>  _S676 = points_7;
            FixedArray<float3 , 5>  _S677 = points_7;
            FixedArray<float3 , 5>  _S678 = points_7;
            points_7[int(0)] = points_7[int(0)];
            points_7[int(1)] = _S675[int(1)];
            points_7[int(2)] = _S676[int(2)];
            points_7[int(3)] = _S677[int(3)];
            points_7[int(4)] = _S678[int(4)];
            _S589 = _S585;
            _S649 = _S624;
        }
        float3  _S679 = _S592 * (points_7[int(0)] + _S589);
        float _S680 = _S679.x + _S679.y + _S679.z;
        float4  _S681 = _S624;
        *&((&_S681)->x) = _S680;
        _S649 = _S649 + _S681;
    }
    else
    {
        _S649 = _S624;
    }
    dpgt_normal_1->primal_0 = (*dpgt_normal_1).primal_0;
    dpgt_normal_1->differential_0 = raydir_9;
    dpdepths_3->primal_0 = (*dpdepths_3).primal_0;
    dpdepths_3->differential_0 = _S649;
    return;
}

inline __device__ void s_bwd_depth_normal_loss_0(float2  _S682, float4  _S683, FixedArray<float, 10>  * _S684, int _S685, bool _S686, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S687, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S688, float _S689)
{
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S690;
    float _S691 = s_primal_ctx_depth_normal_loss_0(_S682, _S683, _S684, _S685, _S686, (*_S687).primal_0, (*_S688).primal_0, &_S690);
    s_bwd_prop_depth_normal_loss_Intermediates_0 _S692 = _S690;
    s_bwd_prop_depth_normal_loss_0(_S682, _S683, _S684, _S685, _S686, _S687, _S688, _S689, &_S692);
    return;
}

inline __device__ void depth_normal_loss_vjp(float2  pix_center_8, float4  intrins_13, FixedArray<float, 10>  dist_coeffs_14, int camera_model_15, bool is_ray_depth_14, float4  depths_3, float3  gt_normal_1, float v_loss_0, float4  * v_depths_1, float3  * v_gt_normal_0)
{
    float4  _S693 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_1;
    (&dp_depths_1)->primal_0 = depths_3;
    (&dp_depths_1)->differential_0 = _S693;
    float3  _S694 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_gt_normal_0;
    (&dp_gt_normal_0)->primal_0 = gt_normal_1;
    (&dp_gt_normal_0)->differential_0 = _S694;
    FixedArray<float, 10>  _S695 = dist_coeffs_14;
    s_bwd_depth_normal_loss_0(pix_center_8, intrins_13, &_S695, camera_model_15, is_ray_depth_14, &dp_depths_1, &dp_gt_normal_0, v_loss_0);
    *v_depths_1 = dp_depths_1.differential_0;
    *v_gt_normal_0 = dp_gt_normal_0.differential_0;
    return;
}

