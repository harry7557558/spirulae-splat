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

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_4)
{
    float _S50 = (*left_0).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_4.x;
    float sum_0 = _S50 + (*left_0).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_4.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_4.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_1;
    float _S51 = (*left_0).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_4.x;
    float sum_2 = _S51 + (*left_0).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_4.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_4.z;
    *&((&right_d_result_1)->y) = sum_3;
    float _S52 = (*left_0).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_4.x;
    float sum_4 = _S52 + (*left_0).primal_0.rows[int(1)].z * dOut_4.y;
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
    float3  _S53 = mul_0(color_matrix_0, rgb_2);
    float _S54 = _S53.x;
    float _S55;
    if(_S54 < 0.00313080009073019f)
    {
        _S55 = _S54 * 12.92000007629394531f;
    }
    else
    {
        _S55 = 1.0549999475479126f * (F32_pow((_S54), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S56 = _S53.y;
    float _S57;
    if(_S56 < 0.00313080009073019f)
    {
        _S57 = _S56 * 12.92000007629394531f;
    }
    else
    {
        _S57 = 1.0549999475479126f * (F32_pow((_S56), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    float _S58 = _S53.z;
    float _S59;
    if(_S58 < 0.00313080009073019f)
    {
        _S59 = _S58 * 12.92000007629394531f;
    }
    else
    {
        _S59 = 1.0549999475479126f * (F32_pow((_S58), (0.4166666567325592f))) - 0.05499999970197678f;
    }
    return make_float3 (_S55, _S57, _S59);
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S60, float3  _S61)
{
    return mul_0(_S60, _S61);
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S62, DiffPair_float_0 * _S63, float _S64)
{
    _d_pow_0(_S62, _S63, _S64);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S65, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S66, float3  _S67)
{
    _d_mul_0(_S65, _S66, _S67);
    return;
}

inline __device__ void s_bwd_prop_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, Matrix<float, 3, 3>  color_matrix_1, float3  _s_dOut_2)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S68 = *dprgb_1;
    float3  _S69 = s_primal_ctx_mul_0(color_matrix_1, (*dprgb_1).primal_0);
    float _S70 = _S69.x;
    float _S71 = _S69.y;
    float _S72 = _S69.z;
    float _S73;
    if(_S72 < 0.00313080009073019f)
    {
        _S73 = 12.92000007629394531f * _s_dOut_2.z;
    }
    else
    {
        float _S74 = 1.0549999475479126f * _s_dOut_2.z;
        DiffPair_float_0 _S75;
        (&_S75)->primal_0 = _S72;
        (&_S75)->differential_0 = 0.0f;
        DiffPair_float_0 _S76;
        (&_S76)->primal_0 = 0.4166666567325592f;
        (&_S76)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S75, &_S76, _S74);
        _S73 = _S75.differential_0;
    }
    float _S77;
    if(_S71 < 0.00313080009073019f)
    {
        _S77 = 12.92000007629394531f * _s_dOut_2.y;
    }
    else
    {
        float _S78 = 1.0549999475479126f * _s_dOut_2.y;
        DiffPair_float_0 _S79;
        (&_S79)->primal_0 = _S71;
        (&_S79)->differential_0 = 0.0f;
        DiffPair_float_0 _S80;
        (&_S80)->primal_0 = 0.4166666567325592f;
        (&_S80)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S79, &_S80, _S78);
        _S77 = _S79.differential_0;
    }
    float _S81;
    if(_S70 < 0.00313080009073019f)
    {
        _S81 = 12.92000007629394531f * _s_dOut_2.x;
    }
    else
    {
        float _S82 = 1.0549999475479126f * _s_dOut_2.x;
        DiffPair_float_0 _S83;
        (&_S83)->primal_0 = _S70;
        (&_S83)->differential_0 = 0.0f;
        DiffPair_float_0 _S84;
        (&_S84)->primal_0 = 0.4166666567325592f;
        (&_S84)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S83, &_S84, _S82);
        _S81 = _S83.differential_0;
    }
    float3  _S85 = make_float3 (_S81, _S77, _S73);
    Matrix<float, 3, 3>  _S86 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S87;
    (&_S87)->primal_0 = color_matrix_1;
    (&_S87)->differential_0 = _S86;
    float3  _S88 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S89;
    (&_S89)->primal_0 = _S68.primal_0;
    (&_S89)->differential_0 = _S88;
    s_bwd_prop_mul_0(&_S87, &_S89, _S85);
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = _S89.differential_0;
    return;
}

inline __device__ void s_bwd_linear_rgb_to_srgb_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S90, Matrix<float, 3, 3>  _S91, float3  _S92)
{
    s_bwd_prop_linear_rgb_to_srgb_0(_S90, _S91, _S92);
    return;
}

inline __device__ float3  linear_rgb_to_srgb_bwd(float3  rgb_3, Matrix<float, 3, 3>  color_matrix_2, float3  v_out_rgb_1)
{
    float3  _S93 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 p_rgb_1;
    (&p_rgb_1)->primal_0 = rgb_3;
    (&p_rgb_1)->differential_0 = _S93;
    s_bwd_linear_rgb_to_srgb_0(&p_rgb_1, color_matrix_2, v_out_rgb_1);
    return p_rgb_1.differential_0;
}

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_5)
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
            *_slang_vector_get_element_ptr(((&result_3)->rows + (r_0)), c_0) = _slang_vector_get_element(x_5.rows[c_0], r_0);
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
    float _S94 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_5).primal_0)))))) * dOut_5;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S94;
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

inline __device__ float dot_0(float3  x_6, float3  y_2)
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
        float result_5 = result_4 + _slang_vector_get_element(x_6, i_3) * _slang_vector_get_element(y_2, i_3);
        i_3 = i_3 + int(1);
        result_4 = result_5;
    }
    return result_4;
}

inline __device__ float dot_1(float2  x_7, float2  y_3)
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
        float result_7 = result_6 + _slang_vector_get_element(x_7, i_4) * _slang_vector_get_element(y_3, i_4);
        i_4 = i_4 + int(1);
        result_6 = result_7;
    }
    return result_6;
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
        float _S95 = (*dist_coeffs_0)[int(3)];
        float _S96 = (*dist_coeffs_0)[int(4)];
        float _S97 = (*dist_coeffs_0)[int(5)];
        float _S98 = (*dist_coeffs_0)[int(6)];
        float _S99 = (*dist_coeffs_0)[int(7)];
        float _S100 = (*dist_coeffs_0)[int(8)];
        float _S101 = (*dist_coeffs_0)[int(9)];
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S102 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
        float _S103 = (*dist_coeffs_0)[int(1)] + r2_0 * _S102;
        float _S104 = (*dist_coeffs_0)[int(0)] + r2_0 * _S103;
        float radial_0 = 1.0f + r2_0 * _S104;
        float _S105 = 2.0f * (*dist_coeffs_0)[int(4)];
        float _S106 = _S105 * u_0;
        float _S107 = 2.0f * u_0;
        float _S108 = 2.0f * (*dist_coeffs_0)[int(5)];
        float _S109 = _S108 * u_0;
        float _S110 = 2.0f * v_0;
        float2  _S111 = q_0 * make_float2 (radial_0) + make_float2 (_S106 * v_0 + (*dist_coeffs_0)[int(5)] * (r2_0 + _S107 * u_0) + (*dist_coeffs_0)[int(6)] * r2_0, _S109 * v_0 + (*dist_coeffs_0)[int(4)] * (r2_0 + _S110 * v_0) + (*dist_coeffs_0)[int(7)] * r2_0);
        float2  r_1 = _S111 + make_float2 ((*dist_coeffs_0)[int(8)] * _S111.x + (*dist_coeffs_0)[int(9)] * _S111.y, 0.0f) - uv_0;
        float _S112 = 0.0f * v_0;
        float s_diff_r2_0 = u_0 + u_0 + (_S112 + _S112);
        float2  _S113 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S104 + (s_diff_r2_0 * _S103 + (s_diff_r2_0 * _S102 + s_diff_r2_0 * _S95 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (_S105 * v_0 + 0.0f * _S106 + (s_diff_r2_0 + (_S107 + _S107)) * _S97 + s_diff_r2_0 * _S98, _S108 * v_0 + 0.0f * _S109 + (s_diff_r2_0 + (_S112 + 0.0f * _S110)) * _S96 + s_diff_r2_0 * _S99);
        float _S114 = 0.0f * u_0;
        float s_diff_r2_1 = _S114 + _S114 + (v_0 + v_0);
        float2  _S115 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S104 + (s_diff_r2_1 * _S103 + (s_diff_r2_1 * _S102 + s_diff_r2_1 * _S95 * r2_0) * r2_0) * r2_0) * q_0 + make_float2 (0.0f * _S105 * v_0 + _S106 + (s_diff_r2_1 + (_S114 + 0.0f * _S107)) * _S97 + s_diff_r2_1 * _S98, 0.0f * _S108 * v_0 + _S109 + (s_diff_r2_1 + (_S110 + _S110)) * _S96 + s_diff_r2_1 * _S99);
        Matrix<float, 2, 2>  _S116 = transpose_0(makeMatrix<float, 2, 2> (_S113 + make_float2 (_S113.x * _S100 + _S113.y * _S101, 0.0f), _S115 + make_float2 (_S115.x * _S100 + _S115.y * _S101, 0.0f)));
        float inv_det_0 = 1.0f / (_S116.rows[int(0)].x * _S116.rows[int(1)].y - _S116.rows[int(0)].y * _S116.rows[int(1)].x);
        float _S117 = r_1.x;
        float _S118 = r_1.y;
        float2  q_1 = q_0 - make_float2 ((_S117 * _S116.rows[int(1)].y - _S118 * _S116.rows[int(0)].y) * inv_det_0, (- _S117 * _S116.rows[int(1)].x + _S118 * _S116.rows[int(0)].x) * inv_det_0);
        i_5 = i_5 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    float _S119 = (*dist_coeffs_0)[int(0)];
    float _S120 = (*dist_coeffs_0)[int(1)];
    float _S121 = (*dist_coeffs_0)[int(2)];
    float _S122 = (*dist_coeffs_0)[int(3)];
    float _S123 = (*dist_coeffs_0)[int(4)];
    float _S124 = (*dist_coeffs_0)[int(5)];
    float _S125 = (*dist_coeffs_0)[int(6)];
    float _S126 = (*dist_coeffs_0)[int(7)];
    float _S127 = (*dist_coeffs_0)[int(8)];
    float _S128 = (*dist_coeffs_0)[int(9)];
    float u_1 = q_0.x;
    float v_1 = q_0.y;
    float _S129 = 0.0f * v_1;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float s_diff_r2_2 = u_1 + u_1 + (_S129 + _S129);
    float _S130 = (*dist_coeffs_0)[int(2)] + r2_1 * (*dist_coeffs_0)[int(3)];
    float _S131 = (*dist_coeffs_0)[int(1)] + r2_1 * _S130;
    float _S132 = (*dist_coeffs_0)[int(0)] + r2_1 * _S131;
    float radial_1 = 1.0f + r2_1 * _S132;
    float _S133 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S134 = _S133 * u_1;
    float _S135 = 2.0f * u_1;
    float _S136 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S137 = _S136 * u_1;
    float _S138 = 2.0f * v_1;
    float2  _S139 = make_float2 (1.0f, 0.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_2 * _S132 + (s_diff_r2_2 * _S131 + (s_diff_r2_2 * _S130 + s_diff_r2_2 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (_S133 * v_1 + 0.0f * _S134 + (s_diff_r2_2 + (_S135 + _S135)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_2 * (*dist_coeffs_0)[int(6)], _S136 * v_1 + 0.0f * _S137 + (s_diff_r2_2 + (_S129 + 0.0f * _S138)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_2 * (*dist_coeffs_0)[int(7)]);
    float _S140 = 0.0f * u_1;
    float s_diff_r2_3 = _S140 + _S140 + (v_1 + v_1);
    float2  _S141 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_3 * _S132 + (s_diff_r2_3 * _S131 + (s_diff_r2_3 * _S130 + s_diff_r2_3 * (*dist_coeffs_0)[int(3)] * r2_1) * r2_1) * r2_1) * q_0 + make_float2 (0.0f * _S133 * v_1 + _S134 + (s_diff_r2_3 + (_S140 + 0.0f * _S135)) * (*dist_coeffs_0)[int(5)] + s_diff_r2_3 * (*dist_coeffs_0)[int(6)], 0.0f * _S136 * v_1 + _S137 + (s_diff_r2_3 + (_S138 + _S138)) * (*dist_coeffs_0)[int(4)] + s_diff_r2_3 * (*dist_coeffs_0)[int(7)]);
    Matrix<float, 2, 2>  _S142 = transpose_0(makeMatrix<float, 2, 2> (_S139 + make_float2 (_S139.x * (*dist_coeffs_0)[int(8)] + _S139.y * (*dist_coeffs_0)[int(9)], 0.0f), _S141 + make_float2 (_S141.x * (*dist_coeffs_0)[int(8)] + _S141.y * (*dist_coeffs_0)[int(9)], 0.0f)));
    bool _S143;
    if((F32_min((determinant_0(_S142)), ((F32_min((_S142.rows[int(0)].x), (_S142.rows[int(1)].y)))))) > 0.0f)
    {
        float u_2 = (*uv_undist_0).x;
        float v_2 = (*uv_undist_0).y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float2  _S144 = *uv_undist_0 * make_float2 (1.0f + r2_2 * (_S119 + r2_2 * (_S120 + r2_2 * (_S121 + r2_2 * _S122)))) + make_float2 (_S133 * u_2 * v_2 + _S124 * (r2_2 + 2.0f * u_2 * u_2) + _S125 * r2_2, _S136 * u_2 * v_2 + _S123 * (r2_2 + 2.0f * v_2 * v_2) + _S126 * r2_2);
        _S143 = (length_0(_S144 + make_float2 (_S127 * _S144.x + _S128 * _S144.y, 0.0f) - uv_0)) < 0.00999999977648258f;
    }
    else
    {
        _S143 = false;
    }
    return _S143;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float3  generate_ray_d2n(float2  pix_pos_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_1, bool is_fisheye_0, bool is_ray_depth_0)
{
    float2  _S145 = (pix_pos_0 - float2 {intrins_0.z, intrins_0.w}) / float2 {intrins_0.x, intrins_0.y};
    float2  uv_1 = _S145;
    FixedArray<float, 10>  _S146 = dist_coeffs_1;
    bool _S147 = undistort_point_0(_S145, &_S146, int(12), &uv_1);
    if(!_S147)
    {
        int3  _S148 = make_int3 (int(0));
        float3  _S149 = make_float3 ((float)_S148.x, (float)_S148.y, (float)_S148.z);
        return _S149;
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
    float _S150 = dOut_7.y;
    float _S151 = dOut_7.z;
    float _S152 = dOut_7.x;
    float _S153 = (*a_0).primal_0.z * _S150 + - (*a_0).primal_0.y * _S151;
    float _S154 = - (*a_0).primal_0.z * _S152 + (*a_0).primal_0.x * _S151;
    float _S155 = (*a_0).primal_0.y * _S152 + - (*a_0).primal_0.x * _S150;
    float3  _S156 = make_float3 (- (*b_0).primal_0.z * _S150 + (*b_0).primal_0.y * _S151, (*b_0).primal_0.z * _S152 + - (*b_0).primal_0.x * _S151, - (*b_0).primal_0.y * _S152 + (*b_0).primal_0.x * _S150);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S156;
    float3  _S157 = make_float3 (_S153, _S154, _S155);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S157;
    return;
}

inline __device__ float3  cross_0(float3  left_2, float3  right_2)
{
    float _S158 = left_2.y;
    float _S159 = right_2.z;
    float _S160 = left_2.z;
    float _S161 = right_2.y;
    float _S162 = right_2.x;
    float _S163 = left_2.x;
    return make_float3 (_S158 * _S159 - _S160 * _S161, _S160 * _S162 - _S163 * _S159, _S163 * _S161 - _S158 * _S162);
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

inline __device__ float3  s_primal_ctx_cross_0(float3  _S164, float3  _S165)
{
    return cross_0(_S164, _S165);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S166, float3  _S167)
{
    return dot_0(_S166, _S167);
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S168, float _S169)
{
    _d_sqrt_0(_S168, _S169);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float _s_dOut_3)
{
    float _S170 = (*dpx_7).primal_0.x;
    float _S171 = (*dpx_7).primal_0.y;
    float _S172 = (*dpx_7).primal_0.z;
    DiffPair_float_0 _S173;
    (&_S173)->primal_0 = _S170 * _S170 + _S171 * _S171 + _S172 * _S172;
    (&_S173)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S173, _s_dOut_3);
    float _S174 = (*dpx_7).primal_0.z * _S173.differential_0;
    float _S175 = _S174 + _S174;
    float _S176 = (*dpx_7).primal_0.y * _S173.differential_0;
    float _S177 = _S176 + _S176;
    float _S178 = (*dpx_7).primal_0.x * _S173.differential_0;
    float _S179 = _S178 + _S178;
    float3  _S180 = make_float3 (0.0f);
    *&((&_S180)->z) = _S175;
    *&((&_S180)->y) = _S177;
    *&((&_S180)->x) = _S179;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S180;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S181, float _S182)
{
    s_bwd_prop_length_impl_0(_S181, _S182);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S183, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S184, float _S185)
{
    _d_dot_0(_S183, _S184, _S185);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S186, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S187, float3  _S188)
{
    _d_cross_0(_S186, _S187, _S188);
    return;
}

inline __device__ void s_bwd_prop_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * dppoints_0, float3  _s_dOut_4)
{
    float3  _S189 = make_float3 (0.0f);
    float3  dx_0 = dppoints_0->primal_0[int(1)] - dppoints_0->primal_0[int(0)];
    float3  _S190 = - (dppoints_0->primal_0[int(3)] - dppoints_0->primal_0[int(2)]);
    float3  _S191 = s_primal_ctx_cross_0(dx_0, _S190);
    bool _S192 = (s_primal_ctx_dot_0(_S191, _S191)) != 0.0f;
    float3  _S193;
    float3  _S194;
    if(_S192)
    {
        float _S195 = length_1(_S191);
        float3  _S196 = make_float3 (_S195);
        _S193 = make_float3 (_S195 * _S195);
        _S194 = _S196;
    }
    else
    {
        _S193 = _S189;
        _S194 = _S189;
    }
    if(_S192)
    {
        float3  _S197 = _s_dOut_4 / _S193;
        float3  _S198 = _S191 * - _S197;
        float3  _S199 = _S194 * _S197;
        float _S200 = _S198.x + _S198.y + _S198.z;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S201;
        (&_S201)->primal_0 = _S191;
        (&_S201)->differential_0 = _S189;
        s_bwd_length_impl_0(&_S201, _S200);
        _S193 = _S199 + _S201.differential_0;
    }
    else
    {
        _S193 = _s_dOut_4;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S202;
    (&_S202)->primal_0 = _S191;
    (&_S202)->differential_0 = _S189;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S203;
    (&_S203)->primal_0 = _S191;
    (&_S203)->differential_0 = _S189;
    s_bwd_prop_dot_0(&_S202, &_S203, 0.0f);
    float3  _S204 = _S203.differential_0 + _S202.differential_0 + _S193;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S205;
    (&_S205)->primal_0 = dx_0;
    (&_S205)->differential_0 = _S189;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S206;
    (&_S206)->primal_0 = _S190;
    (&_S206)->differential_0 = _S189;
    s_bwd_prop_cross_0(&_S205, &_S206, _S204);
    float3  s_diff_dy_T_0 = - _S206.differential_0;
    float3  _S207 = - s_diff_dy_T_0;
    float3  _S208 = - _S205.differential_0;
    FixedArray<float3 , 4>  _S209;
    _S209[int(0)] = _S189;
    _S209[int(1)] = _S189;
    _S209[int(2)] = _S189;
    _S209[int(3)] = _S189;
    _S209[int(2)] = _S207;
    _S209[int(3)] = s_diff_dy_T_0;
    _S209[int(0)] = _S208;
    _S209[int(1)] = _S205.differential_0;
    dppoints_0->primal_0 = dppoints_0->primal_0;
    dppoints_0->differential_0 = _S209;
    return;
}

inline __device__ void s_bwd_points_to_normal_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 * _S210, float3  _S211)
{
    s_bwd_prop_points_to_normal_0(_S210, _S211);
    return;
}

inline __device__ void points_to_normal_vjp(FixedArray<float3 , 4>  points_1, float3  v_normal_0, FixedArray<float3 , 4>  * v_points_0)
{
    FixedArray<float3 , 4>  _S212 = { make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f), make_float3 (0.0f) };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C4x3E_0 dp_points_0;
    (&dp_points_0)->primal_0 = points_1;
    (&dp_points_0)->differential_0 = _S212;
    s_bwd_points_to_normal_0(&dp_points_0, v_normal_0);
    *v_points_0 = (&dp_points_0)->differential_0;
    return;
}

inline __device__ float3  depth_to_normal(float2  pix_center_0, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_2, bool is_fisheye_1, bool is_ray_depth_1, float4  depths_0)
{
    FixedArray<float3 , 4>  points_2;
    float2  _S213 = float2 {intrins_1.z, intrins_1.w};
    float2  _S214 = float2 {intrins_1.x, intrins_1.y};
    float2  _S215 = (pix_center_0 + make_float2 (-1.0f, -0.0f) - _S213) / _S214;
    float2  uv_2 = _S215;
    FixedArray<float, 10>  _S216 = dist_coeffs_2;
    bool _S217 = undistort_point_0(_S215, &_S216, int(12), &uv_2);
    if(!_S217)
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
    float2  _S218 = (pix_center_0 + make_float2 (1.0f, -0.0f) - _S213) / _S214;
    float2  uv_3 = _S218;
    FixedArray<float, 10>  _S219 = dist_coeffs_2;
    bool _S220 = undistort_point_0(_S218, &_S219, int(12), &uv_3);
    if(!_S220)
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
    float2  _S221 = (pix_center_0 + make_float2 (0.0f, -1.0f) - _S213) / _S214;
    float2  uv_4 = _S221;
    FixedArray<float, 10>  _S222 = dist_coeffs_2;
    bool _S223 = undistort_point_0(_S221, &_S222, int(12), &uv_4);
    if(!_S223)
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
    float2  _S224 = (pix_center_0 + make_float2 (0.0f, 1.0f) - _S213) / _S214;
    float2  uv_5 = _S224;
    FixedArray<float, 10>  _S225 = dist_coeffs_2;
    bool _S226 = undistort_point_0(_S224, &_S225, int(12), &uv_5);
    if(!_S226)
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
    float2  _S227;
    bool _S228;
    float2  _S229;
    bool _S230;
    float2  _S231;
    bool _S232;
    float2  _S233;
    bool _S234;
};

inline __device__ float s_primal_ctx_sin_0(float _S235)
{
    return (F32_sin((_S235)));
}

inline __device__ float s_primal_ctx_cos_0(float _S236)
{
    return (F32_cos((_S236)));
}

inline __device__ float3  s_primal_ctx_depth_to_normal_0(float2  pix_center_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_3, bool is_fisheye_2, bool is_ray_depth_2, float4  dpdepths_0, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_0)
{
    float2  _S237 = make_float2 (0.0f);
    _s_diff_ctx_0->_S227 = _S237;
    _s_diff_ctx_0->_S228 = false;
    _s_diff_ctx_0->_S229 = _S237;
    _s_diff_ctx_0->_S230 = false;
    _s_diff_ctx_0->_S231 = _S237;
    _s_diff_ctx_0->_S232 = false;
    _s_diff_ctx_0->_S233 = _S237;
    _s_diff_ctx_0->_S234 = false;
    _s_diff_ctx_0->_S229 = _S237;
    _s_diff_ctx_0->_S230 = false;
    _s_diff_ctx_0->_S231 = _S237;
    _s_diff_ctx_0->_S232 = false;
    _s_diff_ctx_0->_S233 = _S237;
    _s_diff_ctx_0->_S234 = false;
    float3  _S238 = make_float3 (0.0f);
    float2  _S239 = float2 {intrins_2.z, intrins_2.w};
    float2  _S240 = float2 {intrins_2.x, intrins_2.y};
    float2  _S241 = (pix_center_1 + make_float2 (-1.0f, -0.0f) - _S239) / _S240;
    float2  _S242 = _S241;
    bool _S243 = undistort_point_0(_S241, dist_coeffs_3, int(12), &_S242);
    _s_diff_ctx_0->_S227 = _S242;
    _s_diff_ctx_0->_S228 = _S243;
    float2  uv_6 = _S242;
    bool _S244 = !_S243;
    float3  normal_4;
    if(_S244)
    {
        normal_4 = make_float3 (0.0f);
    }
    bool _S245 = !_S244;
    int _S246;
    FixedArray<float3 , 4>  points_3;
    if(_S245)
    {
        float3  raydir_12;
        if(is_fisheye_2)
        {
            float _S247 = length_0(uv_6);
            float3  raydir_13 = make_float3 ((uv_6 / make_float2 ((F32_max((_S247), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S247))).x, (uv_6 / make_float2 ((F32_max((_S247), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S247))).y, s_primal_ctx_cos_0(_S247));
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
        float3  _S248 = make_float3 (dpdepths_0.x) * raydir_12;
        float2  _S249 = (pix_center_1 + make_float2 (1.0f, -0.0f) - _S239) / _S240;
        float2  _S250 = _S249;
        bool _S251 = undistort_point_0(_S249, dist_coeffs_3, int(12), &_S250);
        _s_diff_ctx_0->_S229 = _S250;
        _s_diff_ctx_0->_S230 = _S251;
        float2  uv_7 = _S250;
        bool _S252 = !_S251;
        if(_S252)
        {
            normal_4 = make_float3 (0.0f);
        }
        bool _S253 = !_S252;
        if(_S253)
        {
            if(is_fisheye_2)
            {
                float _S254 = length_0(uv_7);
                float3  raydir_15 = make_float3 ((uv_7 / make_float2 ((F32_max((_S254), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S254))).x, (uv_7 / make_float2 ((F32_max((_S254), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S254))).y, s_primal_ctx_cos_0(_S254));
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
            float3  _S255 = make_float3 (dpdepths_0.y) * raydir_12;
            _S246 = int(2);
            points_3[int(0)] = _S248;
            points_3[int(1)] = _S255;
            points_3[int(2)] = _S238;
            points_3[int(3)] = _S238;
        }
        else
        {
            _S246 = int(0);
            points_3[int(0)] = _S248;
            points_3[int(1)] = _S238;
            points_3[int(2)] = _S238;
            points_3[int(3)] = _S238;
        }
        bool _runFlag_0;
        if(_S246 != int(2))
        {
            _runFlag_0 = false;
        }
        else
        {
            _runFlag_0 = _S245;
            _S246 = int(0);
        }
        if(_runFlag_0)
        {
            float2  _S256 = (pix_center_1 + make_float2 (0.0f, -1.0f) - _S239) / _S240;
            float2  _S257 = _S256;
            bool _S258 = undistort_point_0(_S256, dist_coeffs_3, int(12), &_S257);
            _s_diff_ctx_0->_S231 = _S257;
            _s_diff_ctx_0->_S232 = _S258;
            float2  uv_8 = _S257;
            if(!_S258)
            {
                float3  _S259 = make_float3 (0.0f);
                _runFlag_0 = false;
                _S246 = int(0);
                normal_4 = _S259;
            }
            if(_runFlag_0)
            {
                if(is_fisheye_2)
                {
                    float _S260 = length_0(uv_8);
                    float3  raydir_17 = make_float3 ((uv_8 / make_float2 ((F32_max((_S260), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S260))).x, (uv_8 / make_float2 ((F32_max((_S260), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S260))).y, s_primal_ctx_cos_0(_S260));
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
                float2  _S261 = (pix_center_1 + make_float2 (0.0f, 1.0f) - _S239) / _S240;
                float2  _S262 = _S261;
                bool _S263 = undistort_point_0(_S261, dist_coeffs_3, int(12), &_S262);
                _s_diff_ctx_0->_S233 = _S262;
                _s_diff_ctx_0->_S234 = _S263;
                float2  uv_9 = _S262;
                bool _S264 = !_S263;
                if(_S264)
                {
                    normal_4 = make_float3 (0.0f);
                }
                bool _S265 = !_S264;
                int _S266;
                if(_S265)
                {
                    if(is_fisheye_2)
                    {
                        float _S267 = length_0(uv_9);
                        float3  raydir_19 = make_float3 ((uv_9 / make_float2 ((F32_max((_S267), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S267))).x, (uv_9 / make_float2 ((F32_max((_S267), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S267))).y, s_primal_ctx_cos_0(_S267));
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
                    _S266 = int(2);
                }
                else
                {
                    _S266 = int(0);
                }
                if(_S266 != int(2))
                {
                    _runFlag_0 = false;
                    _S246 = _S266;
                }
                if(_runFlag_0)
                {
                    _S246 = int(1);
                }
            }
        }
    }
    else
    {
        _S246 = int(0);
        points_3[int(0)] = _S238;
        points_3[int(1)] = _S238;
        points_3[int(2)] = _S238;
        points_3[int(3)] = _S238;
    }
    if(!(_S246 != int(1)))
    {
        float3  _S268 = s_primal_ctx_cross_0(points_3[int(1)] - points_3[int(0)], - (points_3[int(3)] - points_3[int(2)]));
        if((s_primal_ctx_dot_0(_S268, _S268)) != 0.0f)
        {
            normal_4 = _S268 / make_float3 (length_1(_S268));
        }
        else
        {
            normal_4 = _S268;
        }
    }
    return normal_4;
}

inline __device__ void s_bwd_prop_depth_to_normal_0(float2  pix_center_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_4, bool is_fisheye_3, bool is_ray_depth_3, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpdepths_1, float3  _s_dOut_5, s_bwd_prop_depth_to_normal_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S269 = *dpdepths_1;
    float3  _S270 = make_float3 (0.0f);
    float2  _S271 = _s_diff_ctx_1->_S227;
    bool _S272 = !!_s_diff_ctx_1->_S228;
    float3  raydir_21;
    float3  raydir_22;
    float3  raydir_23;
    float3  raydir_24;
    int _S273;
    FixedArray<float3 , 4>  points_4;
    bool _runFlag_1;
    bool _runFlag_2;
    bool _runFlag_3;
    bool _S274;
    if(_S272)
    {
        if(is_fisheye_3)
        {
            float _S275 = length_0(_S271);
            float3  raydir_25 = make_float3 ((_S271 / make_float2 ((F32_max((_S275), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S275))).x, (_S271 / make_float2 ((F32_max((_S275), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S275))).y, s_primal_ctx_cos_0(_S275));
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
            float3  raydir_26 = make_float3 (_S271.x, _S271.y, 1.0f);
            if(is_ray_depth_3)
            {
                raydir_21 = normalize_0(raydir_26);
            }
            else
            {
                raydir_21 = raydir_26;
            }
        }
        float3  _S276 = make_float3 (_S269.primal_0.x) * raydir_21;
        float2  _S277 = _s_diff_ctx_1->_S229;
        bool _S278 = !!_s_diff_ctx_1->_S230;
        if(_S278)
        {
            if(is_fisheye_3)
            {
                float _S279 = length_0(_S277);
                float3  raydir_27 = make_float3 ((_S277 / make_float2 ((F32_max((_S279), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S279))).x, (_S277 / make_float2 ((F32_max((_S279), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S279))).y, s_primal_ctx_cos_0(_S279));
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
                float3  raydir_28 = make_float3 (_S277.x, _S277.y, 1.0f);
                if(is_ray_depth_3)
                {
                    raydir_22 = normalize_0(raydir_28);
                }
                else
                {
                    raydir_22 = raydir_28;
                }
            }
            float3  _S280 = make_float3 (_S269.primal_0.y) * raydir_22;
            _S273 = int(2);
            points_4[int(0)] = _S276;
            points_4[int(1)] = _S280;
            points_4[int(2)] = _S270;
            points_4[int(3)] = _S270;
        }
        else
        {
            _S273 = int(0);
            points_4[int(0)] = _S276;
            points_4[int(1)] = _S270;
            points_4[int(2)] = _S270;
            points_4[int(3)] = _S270;
            raydir_22 = _S270;
        }
        if(_S273 != int(2))
        {
            _runFlag_1 = false;
        }
        else
        {
            _runFlag_1 = _S272;
            _S273 = int(0);
        }
        if(_runFlag_1)
        {
            float2  _S281 = _s_diff_ctx_1->_S231;
            if(!_s_diff_ctx_1->_S232)
            {
                _runFlag_2 = false;
                _S273 = int(0);
            }
            else
            {
                _runFlag_2 = _runFlag_1;
            }
            if(_runFlag_2)
            {
                if(is_fisheye_3)
                {
                    float _S282 = length_0(_S281);
                    float3  raydir_29 = make_float3 ((_S281 / make_float2 ((F32_max((_S282), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S282))).x, (_S281 / make_float2 ((F32_max((_S282), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S282))).y, s_primal_ctx_cos_0(_S282));
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
                    float3  raydir_30 = make_float3 (_S281.x, _S281.y, 1.0f);
                    if(is_ray_depth_3)
                    {
                        raydir_23 = normalize_0(raydir_30);
                    }
                    else
                    {
                        raydir_23 = raydir_30;
                    }
                }
                points_4[int(2)] = make_float3 (_S269.primal_0.z) * raydir_23;
                float2  _S283 = _s_diff_ctx_1->_S233;
                bool _S284 = !!_s_diff_ctx_1->_S234;
                int _S285;
                if(_S284)
                {
                    if(is_fisheye_3)
                    {
                        float _S286 = length_0(_S283);
                        float3  raydir_31 = make_float3 ((_S283 / make_float2 ((F32_max((_S286), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S286))).x, (_S283 / make_float2 ((F32_max((_S286), (1.00000001168609742e-07f)))) * make_float2 (s_primal_ctx_sin_0(_S286))).y, s_primal_ctx_cos_0(_S286));
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
                        float3  raydir_32 = make_float3 (_S283.x, _S283.y, 1.0f);
                        if(is_ray_depth_3)
                        {
                            raydir_24 = normalize_0(raydir_32);
                        }
                        else
                        {
                            raydir_24 = raydir_32;
                        }
                    }
                    points_4[int(3)] = make_float3 (_S269.primal_0.w) * raydir_24;
                    _S285 = int(2);
                }
                else
                {
                    _S285 = int(0);
                    raydir_24 = _S270;
                }
                if(_S285 != int(2))
                {
                    _runFlag_3 = false;
                    _S273 = _S285;
                }
                else
                {
                    _runFlag_3 = _runFlag_2;
                }
                if(_runFlag_3)
                {
                    _S273 = int(1);
                }
                float3  _S287 = raydir_23;
                _runFlag_3 = _S284;
                raydir_23 = raydir_24;
                raydir_24 = _S287;
            }
            else
            {
                _runFlag_3 = false;
                raydir_23 = _S270;
                raydir_24 = _S270;
            }
        }
        else
        {
            _runFlag_2 = false;
            _runFlag_3 = false;
            raydir_23 = _S270;
            raydir_24 = _S270;
        }
        float3  _S288 = raydir_21;
        float3  _S289 = raydir_22;
        raydir_21 = raydir_23;
        raydir_22 = raydir_24;
        _S274 = _S278;
        raydir_23 = _S289;
        raydir_24 = _S288;
    }
    else
    {
        _S273 = int(0);
        points_4[int(0)] = _S270;
        points_4[int(1)] = _S270;
        points_4[int(2)] = _S270;
        points_4[int(3)] = _S270;
        _runFlag_1 = false;
        _runFlag_2 = false;
        _runFlag_3 = false;
        raydir_21 = _S270;
        raydir_22 = _S270;
        _S274 = false;
        raydir_23 = _S270;
        raydir_24 = _S270;
    }
    bool _S290 = !(_S273 != int(1));
    float3  _S291;
    float3  _S292;
    float3  _S293;
    float3  _S294;
    float3  _S295;
    bool _S296;
    if(_S290)
    {
        float3  dx_1 = points_4[int(1)] - points_4[int(0)];
        float3  _S297 = - (points_4[int(3)] - points_4[int(2)]);
        float3  _S298 = s_primal_ctx_cross_0(dx_1, _S297);
        bool _S299 = (s_primal_ctx_dot_0(_S298, _S298)) != 0.0f;
        if(_S299)
        {
            float _S300 = length_1(_S298);
            float3  _S301 = make_float3 (_S300);
            _S291 = make_float3 (_S300 * _S300);
            _S292 = _S301;
        }
        else
        {
            _S291 = _S270;
            _S292 = _S270;
        }
        float3  _S302 = _S292;
        _S296 = _S299;
        _S292 = _S298;
        _S293 = _S302;
        _S294 = dx_1;
        _S295 = _S297;
    }
    else
    {
        _S296 = false;
        _S291 = _S270;
        _S292 = _S270;
        _S293 = _S270;
        _S294 = _S270;
        _S295 = _S270;
    }
    float4  _S303 = make_float4 (0.0f);
    if(_S290)
    {
        if(_S296)
        {
            float3  _S304 = _s_dOut_5 / _S291;
            float3  _S305 = _S292 * - _S304;
            float3  _S306 = _S293 * _S304;
            float _S307 = _S305.x + _S305.y + _S305.z;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S308;
            (&_S308)->primal_0 = _S292;
            (&_S308)->differential_0 = _S270;
            s_bwd_length_impl_0(&_S308, _S307);
            _S291 = _S306 + _S308.differential_0;
        }
        else
        {
            _S291 = _s_dOut_5;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S309;
        (&_S309)->primal_0 = _S292;
        (&_S309)->differential_0 = _S270;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S310;
        (&_S310)->primal_0 = _S292;
        (&_S310)->differential_0 = _S270;
        s_bwd_prop_dot_0(&_S309, &_S310, 0.0f);
        float3  _S311 = _S310.differential_0 + _S309.differential_0 + _S291;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S312;
        (&_S312)->primal_0 = _S294;
        (&_S312)->differential_0 = _S270;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S313;
        (&_S313)->primal_0 = _S295;
        (&_S313)->differential_0 = _S270;
        s_bwd_prop_cross_0(&_S312, &_S313, _S311);
        float3  s_diff_dy_T_1 = - _S313.differential_0;
        float3  _S314 = - s_diff_dy_T_1;
        float3  _S315 = - _S312.differential_0;
        FixedArray<float3 , 4>  _S316;
        _S316[int(0)] = _S270;
        _S316[int(1)] = _S270;
        _S316[int(2)] = _S270;
        _S316[int(3)] = _S270;
        _S316[int(2)] = _S314;
        _S316[int(3)] = s_diff_dy_T_1;
        _S316[int(0)] = _S315;
        _S316[int(1)] = _S312.differential_0;
        points_4[int(0)] = _S316[int(0)];
        points_4[int(1)] = _S316[int(1)];
        points_4[int(2)] = _S316[int(2)];
        points_4[int(3)] = _S316[int(3)];
    }
    else
    {
        points_4[int(0)] = _S270;
        points_4[int(1)] = _S270;
        points_4[int(2)] = _S270;
        points_4[int(3)] = _S270;
    }
    float4  _S317;
    if(_S272)
    {
        if(_runFlag_1)
        {
            if(_runFlag_2)
            {
                FixedArray<float3 , 4>  _S318 = points_4;
                FixedArray<float3 , 4>  _S319 = points_4;
                FixedArray<float3 , 4>  _S320 = points_4;
                FixedArray<float3 , 4>  _S321 = points_4;
                if(_runFlag_3)
                {
                    float3  _S322 = raydir_21 * _S321[int(3)];
                    float _S323 = _S322.x + _S322.y + _S322.z;
                    float4  _S324 = _S303;
                    *&((&_S324)->w) = _S323;
                    points_4[int(0)] = _S318[int(0)];
                    points_4[int(1)] = _S319[int(1)];
                    points_4[int(2)] = _S320[int(2)];
                    points_4[int(3)] = _S270;
                    _S317 = _S324;
                }
                else
                {
                    points_4[int(0)] = _S318[int(0)];
                    points_4[int(1)] = _S319[int(1)];
                    points_4[int(2)] = _S320[int(2)];
                    points_4[int(3)] = _S321[int(3)];
                    _S317 = _S303;
                }
                float3  _S325 = raydir_22 * points_4[int(2)];
                float _S326 = _S325.x + _S325.y + _S325.z;
                FixedArray<float3 , 4>  _S327 = points_4;
                FixedArray<float3 , 4>  _S328 = points_4;
                float4  _S329 = _S303;
                *&((&_S329)->z) = _S326;
                float4  _S330 = _S317 + _S329;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S327[int(1)];
                points_4[int(2)] = _S270;
                points_4[int(3)] = _S328[int(3)];
                _S317 = _S330;
            }
            else
            {
                FixedArray<float3 , 4>  _S331 = points_4;
                FixedArray<float3 , 4>  _S332 = points_4;
                FixedArray<float3 , 4>  _S333 = points_4;
                points_4[int(0)] = points_4[int(0)];
                points_4[int(1)] = _S331[int(1)];
                points_4[int(2)] = _S332[int(2)];
                points_4[int(3)] = _S333[int(3)];
                _S317 = _S303;
            }
        }
        else
        {
            FixedArray<float3 , 4>  _S334 = points_4;
            FixedArray<float3 , 4>  _S335 = points_4;
            FixedArray<float3 , 4>  _S336 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S334[int(1)];
            points_4[int(2)] = _S335[int(2)];
            points_4[int(3)] = _S336[int(3)];
            _S317 = _S303;
        }
        if(_S274)
        {
            FixedArray<float3 , 4>  _S337 = points_4;
            float3  _S338 = raydir_23 * points_4[int(1)];
            float _S339 = _S338.x + _S338.y + _S338.z;
            float4  _S340 = _S303;
            *&((&_S340)->y) = _S339;
            float4  _S341 = _S317 + _S340;
            points_4[int(0)] = _S270;
            points_4[int(1)] = _S270;
            points_4[int(2)] = _S270;
            points_4[int(3)] = _S270;
            raydir_21 = _S337[int(0)];
            _S317 = _S341;
        }
        else
        {
            FixedArray<float3 , 4>  _S342 = points_4;
            FixedArray<float3 , 4>  _S343 = points_4;
            FixedArray<float3 , 4>  _S344 = points_4;
            points_4[int(0)] = points_4[int(0)];
            points_4[int(1)] = _S342[int(1)];
            points_4[int(2)] = _S343[int(2)];
            points_4[int(3)] = _S344[int(3)];
            raydir_21 = _S270;
        }
        float3  _S345 = raydir_24 * (points_4[int(0)] + raydir_21);
        float _S346 = _S345.x + _S345.y + _S345.z;
        float4  _S347 = _S303;
        *&((&_S347)->x) = _S346;
        _S317 = _S317 + _S347;
    }
    else
    {
        _S317 = _S303;
    }
    dpdepths_1->primal_0 = (*dpdepths_1).primal_0;
    dpdepths_1->differential_0 = _S317;
    return;
}

inline __device__ void s_bwd_depth_to_normal_0(float2  _S348, float4  _S349, FixedArray<float, 10>  * _S350, bool _S351, bool _S352, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S353, float3  _S354)
{
    s_bwd_prop_depth_to_normal_Intermediates_0 _S355;
    float3  _S356 = s_primal_ctx_depth_to_normal_0(_S348, _S349, _S350, _S351, _S352, (*_S353).primal_0, &_S355);
    s_bwd_prop_depth_to_normal_Intermediates_0 _S357 = _S355;
    s_bwd_prop_depth_to_normal_0(_S348, _S349, _S350, _S351, _S352, _S353, _S354, &_S357);
    return;
}

inline __device__ void depth_to_normal_vjp(float2  pix_center_3, float4  intrins_4, FixedArray<float, 10>  dist_coeffs_5, bool is_fisheye_4, bool is_ray_depth_4, float4  depths_1, float3  v_normal_1, float4  * v_depths_0)
{
    float4  _S358 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S358;
    FixedArray<float, 10>  _S359 = dist_coeffs_5;
    s_bwd_depth_to_normal_0(pix_center_3, intrins_4, &_S359, is_fisheye_4, is_ray_depth_4, &dp_depths_0, v_normal_1);
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ float ray_depth_to_linear_depth_factor(float2  pix_center_4, float4  intrins_5, FixedArray<float, 10>  dist_coeffs_6, bool is_fisheye_5)
{
    float2  _S360 = (pix_center_4 - float2 {intrins_5.z, intrins_5.w}) / float2 {intrins_5.x, intrins_5.y};
    float2  uv_10 = _S360;
    FixedArray<float, 10>  _S361 = dist_coeffs_6;
    bool _S362 = undistort_point_0(_S360, &_S361, int(12), &uv_10);
    if(!_S362)
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

