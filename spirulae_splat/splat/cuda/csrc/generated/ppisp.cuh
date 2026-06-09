#pragma once

#include "slang.cuh"

struct VignettingChannelParams_0
{
    float cx_0;
    float cy_0;
    float alpha0_0;
    float alpha1_0;
    float alpha2_0;
};

inline __device__ VignettingChannelParams_0 VignettingChannelParams_x24_syn_dzero_0()
{
    VignettingChannelParams_0 result_0;
    (&result_0)->cx_0 = 0.0f;
    (&result_0)->cy_0 = 0.0f;
    (&result_0)->alpha0_0 = 0.0f;
    (&result_0)->alpha1_0 = 0.0f;
    (&result_0)->alpha2_0 = 0.0f;
    return result_0;
}

struct ColorPPISPParams_0
{
    float2  b_0;
    float2  r_0;
    float2  g_0;
    float2  n_0;
};

inline __device__ ColorPPISPParams_0 ColorPPISPParams_x24_syn_dzero_0()
{
    ColorPPISPParams_0 result_1;
    float2  _S1 = make_float2 (0.0f);
    (&result_1)->b_0 = _S1;
    (&result_1)->r_0 = _S1;
    (&result_1)->g_0 = _S1;
    (&result_1)->n_0 = _S1;
    return result_1;
}

struct PPISPParamsNoCRF_0
{
    float exposure_0;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_0;
    ColorPPISPParams_0 color_params_0;
};

inline __device__ PPISPParamsNoCRF_0 PPISPParamsNoCRF_x24_syn_dzero_0()
{
    PPISPParamsNoCRF_0 result_2;
    (&result_2)->exposure_0 = 0.0f;
    VignettingChannelParams_0 _S2 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_2)->vignette_params_0[int(0)] = _S2;
    (&result_2)->vignette_params_0[int(1)] = _S2;
    (&result_2)->vignette_params_0[int(2)] = _S2;
    (&result_2)->color_params_0 = ColorPPISPParams_x24_syn_dzero_0();
    return result_2;
}

struct RQSCRFPPISPChannelParams_0
{
    float g0_0;
    float g1_0;
    float x0_0;
    float y0_0;
    float gc_0;
};

inline __device__ RQSCRFPPISPChannelParams_0 RQSCRFPPISPChannelParams_x24_syn_dzero_0()
{
    RQSCRFPPISPChannelParams_0 result_3;
    (&result_3)->g0_0 = 0.0f;
    (&result_3)->g1_0 = 0.0f;
    (&result_3)->x0_0 = 0.0f;
    (&result_3)->y0_0 = 0.0f;
    (&result_3)->gc_0 = 0.0f;
    return result_3;
}

struct PPISPParamsRQS_0
{
    float exposure_1;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_1;
    ColorPPISPParams_0 color_params_1;
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  crf_params_0;
};

inline __device__ PPISPParamsRQS_0 PPISPParamsRQS_x24_syn_dzero_0()
{
    PPISPParamsRQS_0 result_4;
    (&result_4)->exposure_1 = 0.0f;
    VignettingChannelParams_0 _S3 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_4)->vignette_params_1[int(0)] = _S3;
    (&result_4)->vignette_params_1[int(1)] = _S3;
    (&result_4)->vignette_params_1[int(2)] = _S3;
    (&result_4)->color_params_1 = ColorPPISPParams_x24_syn_dzero_0();
    RQSCRFPPISPChannelParams_0 _S4 = RQSCRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_4)->crf_params_0[int(0)] = _S4;
    (&result_4)->crf_params_0[int(1)] = _S4;
    (&result_4)->crf_params_0[int(2)] = _S4;
    return result_4;
}

struct CRFPPISPChannelParams_0
{
    float toe_0;
    float shoulder_0;
    float gamma_0;
    float center_0;
};

inline __device__ CRFPPISPChannelParams_0 CRFPPISPChannelParams_x24_syn_dzero_0()
{
    CRFPPISPChannelParams_0 result_5;
    (&result_5)->toe_0 = 0.0f;
    (&result_5)->shoulder_0 = 0.0f;
    (&result_5)->gamma_0 = 0.0f;
    (&result_5)->center_0 = 0.0f;
    return result_5;
}

struct PPISPParams_0
{
    float exposure_2;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_2;
    ColorPPISPParams_0 color_params_2;
    FixedArray<CRFPPISPChannelParams_0, 3>  crf_params_1;
};

inline __device__ PPISPParams_0 PPISPParams_x24_syn_dzero_0()
{
    PPISPParams_0 result_6;
    (&result_6)->exposure_2 = 0.0f;
    VignettingChannelParams_0 _S5 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_6)->vignette_params_2[int(0)] = _S5;
    (&result_6)->vignette_params_2[int(1)] = _S5;
    (&result_6)->vignette_params_2[int(2)] = _S5;
    (&result_6)->color_params_2 = ColorPPISPParams_x24_syn_dzero_0();
    CRFPPISPChannelParams_0 _S6 = CRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_6)->crf_params_1[int(0)] = _S6;
    (&result_6)->crf_params_1[int(1)] = _S6;
    (&result_6)->crf_params_1[int(2)] = _S6;
    return result_6;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S7 = (F32_exp2(((*dpx_0).primal_0))) * 0.69314718246459961f * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S7;
    return;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S8 = *dpx_1;
    float _S9;
    if(((*dpx_1).primal_0) > ((*dpy_0).primal_0))
    {
        _S9 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_0).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S8.primal_0;
    dpx_1->differential_0 = _S9;
    DiffPair_float_0 _S10 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S8.primal_0))
    {
        _S9 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_1).primal_0))
        {
            _S9 = 0.0f;
        }
        else
        {
            _S9 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S10.primal_0;
    dpy_0->differential_0 = _S9;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_2)
{
    DiffPair_float_0 _S11 = *dpx_2;
    bool _S12;
    if(((*dpx_2).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S12 = ((*dpx_2).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S12 = false;
    }
    float _S13;
    if(_S12)
    {
        _S13 = dOut_2;
    }
    else
    {
        _S13 = 0.0f;
    }
    dpx_2->primal_0 = _S11.primal_0;
    dpx_2->differential_0 = _S13;
    DiffPair_float_0 _S14 = *dpMin_0;
    if((_S11.primal_0) < ((*dpMin_0).primal_0))
    {
        _S13 = dOut_2;
    }
    else
    {
        _S13 = 0.0f;
    }
    dpMin_0->primal_0 = _S14.primal_0;
    dpMin_0->differential_0 = _S13;
    DiffPair_float_0 _S15 = *dpMax_0;
    if(((*dpx_2).primal_0) > ((*dpMax_0).primal_0))
    {
        _S13 = dOut_2;
    }
    else
    {
        _S13 = 0.0f;
    }
    dpMax_0->primal_0 = _S15.primal_0;
    dpMax_0->differential_0 = _S13;
    return;
}

inline __device__ float clamp_0(float x_0, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_0), (minBound_0)))), (maxBound_0)));
}

struct DiffPair_matrixx3Cfloatx2C2x2C2x3E_0
{
    Matrix<float, 2, 2>  primal_0;
    Matrix<float, 2, 2>  differential_0;
};

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * right_0, float2  dOut_3)
{
    float _S16 = (*left_0).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 2, 2>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_3.x;
    float sum_0 = _S16 + (*left_0).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_3.y;
    float2  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_0;
    float _S17 = (*left_0).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_3.x;
    float sum_1 = _S17 + (*left_0).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_3.y;
    *&((&right_d_result_0)->y) = sum_1;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_1, float3  dOut_4)
{
    float _S18 = (*left_1).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_1).primal_0.x * dOut_4.x;
    float sum_2 = _S18 + (*left_1).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_1).primal_0.x * dOut_4.y;
    float sum_3 = sum_2 + (*left_1).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_1).primal_0.x * dOut_4.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_3;
    float _S19 = (*left_1).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_1).primal_0.y * dOut_4.x;
    float sum_4 = _S19 + (*left_1).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_1).primal_0.y * dOut_4.y;
    float sum_5 = sum_4 + (*left_1).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_1).primal_0.y * dOut_4.z;
    *&((&right_d_result_1)->y) = sum_5;
    float _S20 = (*left_1).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_1).primal_0.z * dOut_4.x;
    float sum_6 = _S20 + (*left_1).primal_0.rows[int(1)].z * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = (*right_1).primal_0.z * dOut_4.y;
    float sum_7 = sum_6 + (*left_1).primal_0.rows[int(2)].z * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = (*right_1).primal_0.z * dOut_4.z;
    *&((&right_d_result_1)->z) = sum_7;
    left_1->primal_0 = (*left_1).primal_0;
    left_1->differential_0 = left_d_result_1;
    right_1->primal_0 = (*right_1).primal_0;
    right_1->differential_0 = right_d_result_1;
    return;
}

inline __device__ float2  mul_0(Matrix<float, 2, 2>  left_2, float2  right_2)
{
    float2  result_7;
    int i_0 = int(0);
    for(;;)
    {
        if(i_0 < int(2))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_8 = 0.0f;
        for(;;)
        {
            if(j_0 < int(2))
            {
            }
            else
            {
                break;
            }
            float sum_9 = sum_8 + _slang_vector_get_element(left_2.rows[i_0], j_0) * _slang_vector_get_element(right_2, j_0);
            j_0 = j_0 + int(1);
            sum_8 = sum_9;
        }
        *_slang_vector_get_element_ptr(&result_7, i_0) = sum_8;
        i_0 = i_0 + int(1);
    }
    return result_7;
}

inline __device__ float3  mul_1(Matrix<float, 3, 3>  left_3, float3  right_3)
{
    float3  result_8;
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
        int j_1 = int(0);
        float sum_10 = 0.0f;
        for(;;)
        {
            if(j_1 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_11 = sum_10 + _slang_vector_get_element(left_3.rows[i_1], j_1) * _slang_vector_get_element(right_3, j_1);
            j_1 = j_1 + int(1);
            sum_10 = sum_11;
        }
        *_slang_vector_get_element_ptr(&result_8, i_1) = sum_10;
        i_1 = i_1 + int(1);
    }
    return result_8;
}

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_4, Matrix<float, 3, 3>  dOut_5)
{
    Matrix<float, 3, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_2;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_4).primal_0.rows[int(0)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_4).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_4).primal_0.rows[int(1)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_4).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_4).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_4).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].z;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_2;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_2;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_3(Matrix<float, 3, 3>  left_5, Matrix<float, 3, 3>  right_5)
{
    Matrix<float, 3, 3>  result_9;
    int r_1 = int(0);
    for(;;)
    {
        if(r_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_2 = int(0);
            float sum_12 = 0.0f;
            for(;;)
            {
                if(i_2 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_13 = sum_12 + _slang_vector_get_element(left_5.rows[r_1], i_2) * _slang_vector_get_element(right_5.rows[i_2], c_0);
                i_2 = i_2 + int(1);
                sum_12 = sum_13;
            }
            *_slang_vector_get_element_ptr(((&result_9)->rows + (r_1)), c_0) = sum_12;
            c_0 = c_0 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_9;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_6)
{
    float _S21 = dOut_6.y;
    float _S22 = dOut_6.z;
    float _S23 = dOut_6.x;
    float _S24 = (*a_0).primal_0.z * _S21 + - (*a_0).primal_0.y * _S22;
    float _S25 = - (*a_0).primal_0.z * _S23 + (*a_0).primal_0.x * _S22;
    float _S26 = (*a_0).primal_0.y * _S23 + - (*a_0).primal_0.x * _S21;
    float3  _S27 = make_float3 (- (*b_1).primal_0.z * _S21 + (*b_1).primal_0.y * _S22, (*b_1).primal_0.z * _S23 + - (*b_1).primal_0.x * _S22, - (*b_1).primal_0.y * _S23 + (*b_1).primal_0.x * _S21);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S27;
    float3  _S28 = make_float3 (_S24, _S25, _S26);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S28;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S29 = left_6.y;
    float _S30 = right_6.z;
    float _S31 = left_6.z;
    float _S32 = right_6.y;
    float _S33 = right_6.x;
    float _S34 = left_6.x;
    return make_float3 (_S29 * _S30 - _S31 * _S32, _S31 * _S33 - _S34 * _S30, _S34 * _S32 - _S29 * _S33);
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_7)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_7;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_3).primal_0.x * dOut_7;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_7;
    *&((&y_d_result_0)->y) = (*dpx_3).primal_0.y * dOut_7;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_7;
    *&((&y_d_result_0)->z) = (*dpx_3).primal_0.z * dOut_7;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float3  x_1, float3  y_0)
{
    int i_3 = int(0);
    float result_10 = 0.0f;
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_11 = result_10 + _slang_vector_get_element(x_1, i_3) * _slang_vector_get_element(y_0, i_3);
        i_3 = i_3 + int(1);
        result_10 = result_11;
    }
    return result_10;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_4, float dOut_8)
{
    float _S35 = _slang_select(((*dpx_4).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_4).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S35;
    return;
}

inline __device__ float3  min_0(float3  x_2, float3  y_1)
{
    float3  result_12;
    int i_4 = int(0);
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_12, i_4) = (F32_min((_slang_vector_get_element(x_2, i_4)), (_slang_vector_get_element(y_1, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_12;
}

inline __device__ float3  max_0(float3  x_3, float3  y_2)
{
    float3  result_13;
    int i_5 = int(0);
    for(;;)
    {
        if(i_5 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_13, i_5) = (F32_max((_slang_vector_get_element(x_3, i_5)), (_slang_vector_get_element(y_2, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_13;
}

inline __device__ void _d_clamp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpz_0, float3  dOut_9)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_5).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_0;
    (&middle_dp_0)->primal_0 = (*dpy_2).primal_0.x;
    (&middle_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpz_0).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_0, &middle_dp_0, &right_dp_0, dOut_9.x);
    float3  left_d_result_3;
    *&((&left_d_result_3)->x) = left_dp_0.differential_0;
    float3  middle_d_result_0;
    *&((&middle_d_result_0)->x) = middle_dp_0.differential_0;
    float3  right_d_result_3;
    *&((&right_d_result_3)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_5).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_1;
    (&middle_dp_1)->primal_0 = (*dpy_2).primal_0.y;
    (&middle_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpz_0).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_1, &middle_dp_1, &right_dp_1, dOut_9.y);
    *&((&left_d_result_3)->y) = left_dp_1.differential_0;
    *&((&middle_d_result_0)->y) = middle_dp_1.differential_0;
    *&((&right_d_result_3)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_5).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 middle_dp_2;
    (&middle_dp_2)->primal_0 = (*dpy_2).primal_0.z;
    (&middle_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpz_0).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_clamp_0(&left_dp_2, &middle_dp_2, &right_dp_2, dOut_9.z);
    *&((&left_d_result_3)->z) = left_dp_2.differential_0;
    *&((&middle_d_result_0)->z) = middle_dp_2.differential_0;
    *&((&right_d_result_3)->z) = right_dp_2.differential_0;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = left_d_result_3;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = middle_d_result_0;
    dpz_0->primal_0 = (*dpz_0).primal_0;
    dpz_0->differential_0 = right_d_result_3;
    return;
}

inline __device__ float3  clamp_1(float3  x_4, float3  minBound_1, float3  maxBound_1)
{
    return min_0(max_0(x_4, minBound_1), maxBound_1);
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_6, float dOut_10)
{
    float _S36 = (F32_exp(((*dpx_6).primal_0))) * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S36;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_7, float dOut_11)
{
    float _S37 = 1.0f / (*dpx_7).primal_0 * dOut_11;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S37;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_12)
{
    float _S38 = (1.0f - (*dps_0).primal_0) * dOut_12;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S38;
    DiffPair_float_0 _S39 = *dpy_3;
    float _S40 = (*dps_0).primal_0 * dOut_12;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S40;
    float _S41 = (_S39.primal_0 - (*dpx_8).primal_0) * dOut_12;
    dps_0->primal_0 = _S39.primal_0;
    dps_0->differential_0 = _S41;
    return;
}

inline __device__ float lerp_0(float x_5, float y_3, float s_0)
{
    return x_5 + (y_3 - x_5) * s_0;
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_9, DiffPair_float_0 * dpy_4, float dOut_13)
{
    if(((*dpx_9).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_9->primal_0 = (*dpx_9).primal_0;
        dpx_9->differential_0 = 0.0f;
        dpy_4->primal_0 = (*dpy_4).primal_0;
        dpy_4->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_9).primal_0), ((*dpy_4).primal_0)));
        DiffPair_float_0 _S42 = *dpx_9;
        float _S43 = val_0 * (*dpy_4).primal_0 / (*dpx_9).primal_0 * dOut_13;
        dpx_9->primal_0 = (*dpx_9).primal_0;
        dpx_9->differential_0 = _S43;
        float _S44 = val_0 * (F32_log((_S42.primal_0))) * dOut_13;
        dpy_4->primal_0 = (*dpy_4).primal_0;
        dpy_4->differential_0 = _S44;
    }
    return;
}

inline __device__ float3  apply_ppisp(float3  rgb_in_0, float2  pix_coord_0, float2  image_center_0, float2  img_size_0, FixedArray<float, 36>  params_0)
{
    PPISPParams_0 p_0;
    (&p_0)->exposure_2 = params_0[int(0)];
    (&(&p_0)->vignette_params_2[int(0)])->cx_0 = params_0[int(1)];
    (&(&p_0)->vignette_params_2[int(0)])->cy_0 = params_0[int(2)];
    (&(&p_0)->vignette_params_2[int(0)])->alpha0_0 = params_0[int(3)];
    (&(&p_0)->vignette_params_2[int(0)])->alpha1_0 = params_0[int(4)];
    (&(&p_0)->vignette_params_2[int(0)])->alpha2_0 = params_0[int(5)];
    (&(&p_0)->vignette_params_2[int(1)])->cx_0 = params_0[int(6)];
    (&(&p_0)->vignette_params_2[int(1)])->cy_0 = params_0[int(7)];
    (&(&p_0)->vignette_params_2[int(1)])->alpha0_0 = params_0[int(8)];
    (&(&p_0)->vignette_params_2[int(1)])->alpha1_0 = params_0[int(9)];
    (&(&p_0)->vignette_params_2[int(1)])->alpha2_0 = params_0[int(10)];
    (&(&p_0)->vignette_params_2[int(2)])->cx_0 = params_0[int(11)];
    (&(&p_0)->vignette_params_2[int(2)])->cy_0 = params_0[int(12)];
    (&(&p_0)->vignette_params_2[int(2)])->alpha0_0 = params_0[int(13)];
    (&(&p_0)->vignette_params_2[int(2)])->alpha1_0 = params_0[int(14)];
    (&(&p_0)->vignette_params_2[int(2)])->alpha2_0 = params_0[int(15)];
    *&((&(&(&p_0)->color_params_2)->b_0)->x) = params_0[int(16)];
    *&((&(&(&p_0)->color_params_2)->b_0)->y) = params_0[int(17)];
    *&((&(&(&p_0)->color_params_2)->r_0)->x) = params_0[int(18)];
    *&((&(&(&p_0)->color_params_2)->r_0)->y) = params_0[int(19)];
    *&((&(&(&p_0)->color_params_2)->g_0)->x) = params_0[int(20)];
    *&((&(&(&p_0)->color_params_2)->g_0)->y) = params_0[int(21)];
    *&((&(&(&p_0)->color_params_2)->n_0)->x) = params_0[int(22)];
    *&((&(&(&p_0)->color_params_2)->n_0)->y) = params_0[int(23)];
    (&(&p_0)->crf_params_1[int(0)])->toe_0 = params_0[int(24)];
    (&(&p_0)->crf_params_1[int(0)])->shoulder_0 = params_0[int(25)];
    (&(&p_0)->crf_params_1[int(0)])->gamma_0 = params_0[int(26)];
    (&(&p_0)->crf_params_1[int(0)])->center_0 = params_0[int(27)];
    (&(&p_0)->crf_params_1[int(1)])->toe_0 = params_0[int(28)];
    (&(&p_0)->crf_params_1[int(1)])->shoulder_0 = params_0[int(29)];
    (&(&p_0)->crf_params_1[int(1)])->gamma_0 = params_0[int(30)];
    (&(&p_0)->crf_params_1[int(1)])->center_0 = params_0[int(31)];
    (&(&p_0)->crf_params_1[int(2)])->toe_0 = params_0[int(32)];
    (&(&p_0)->crf_params_1[int(2)])->shoulder_0 = params_0[int(33)];
    (&(&p_0)->crf_params_1[int(2)])->gamma_0 = params_0[int(34)];
    (&(&p_0)->crf_params_1[int(2)])->center_0 = params_0[int(35)];
    PPISPParams_0 _S45 = p_0;
    float _S46 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S47 = (pix_coord_0.x - image_center_0.x) / _S46;
    float _S48 = (pix_coord_0.y - image_center_0.y) / _S46;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_2))));
    float dx_0 = _S47 - p_0.vignette_params_2[int(0)].cx_0;
    float dy_0 = _S48 - p_0.vignette_params_2[int(0)].cy_0;
    float r2_0 = dx_0 * dx_0 + dy_0 * dy_0;
    float r4_0 = r2_0 * r2_0;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_2[int(0)].alpha2_0 * (r4_0 * r2_0) + p_0.vignette_params_2[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_2[int(0)].alpha0_0 * r2_0 + 1.0f, 0.0f, 1.0f);
    float dx_1 = _S47 - p_0.vignette_params_2[int(1)].cx_0;
    float dy_1 = _S48 - p_0.vignette_params_2[int(1)].cy_0;
    float r2_1 = dx_1 * dx_1 + dy_1 * dy_1;
    float r4_1 = r2_1 * r2_1;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_2[int(1)].alpha2_0 * (r4_1 * r2_1) + p_0.vignette_params_2[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_2[int(1)].alpha0_0 * r2_1 + 1.0f, 0.0f, 1.0f);
    float dx_2 = _S47 - p_0.vignette_params_2[int(2)].cx_0;
    float dy_2 = _S48 - p_0.vignette_params_2[int(2)].cy_0;
    float r2_2 = dx_2 * dx_2 + dy_2 * dy_2;
    float r4_2 = r2_2 * r2_2;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_2[int(2)].alpha2_0 * (r4_2 * r2_2) + p_0.vignette_params_2[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_2[int(2)].alpha0_0 * r2_2 + 1.0f, 0.0f, 1.0f);
    float3  _S49 = rgb_out_0;
    float2  bd_0 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_2.b_0);
    float2  rd_0 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_2.r_0);
    float2  gd_0 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_2.g_0);
    float2  nd_0 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_2.n_0);
    float _S50 = 0.3333333432674408f + nd_0.x;
    float _S51 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_0 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S51, 1.0f, 0.0f, - _S50, - _S51, _S50, 0.0f), T_0);
    float3  r0_0 = make_float3 (M_0.rows[int(0)].x, M_0.rows[int(0)].y, M_0.rows[int(0)].z);
    float3  r1_0 = make_float3 (M_0.rows[int(1)].x, M_0.rows[int(1)].y, M_0.rows[int(1)].z);
    float3  r2_3 = make_float3 (M_0.rows[int(2)].x, M_0.rows[int(2)].y, M_0.rows[int(2)].z);
    float3  lambda_v_0 = cross_0(r0_0, r1_0);
    float3  lambda_v_1;
    if((dot_0(lambda_v_0, lambda_v_0)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_2 = cross_0(r0_0, r2_3);
        if((dot_0(lambda_v_2, lambda_v_2)) < 9.99999968265522539e-21f)
        {
            lambda_v_1 = cross_0(r1_0, r2_3);
        }
        else
        {
            lambda_v_1 = lambda_v_2;
        }
    }
    else
    {
        lambda_v_1 = lambda_v_0;
    }
    Matrix<float, 3, 3>  H_0 = mul_3(mul_3(T_0, makeMatrix<float, 3, 3> (lambda_v_1.x, 0.0f, 0.0f, 0.0f, lambda_v_1.y, 0.0f, 0.0f, 0.0f, lambda_v_1.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_1;
    if((F32_abs((H_0.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_1 = H_0 * makeMatrix<float, 3, 3> (1.0f / H_0.rows[int(2)].z);
    }
    else
    {
        H_1 = H_0;
    }
    float _S52 = _S49.x;
    float _S53 = _S49.y;
    float intensity_0 = _S52 + _S53 + _S49.z;
    float3  rgi_out_0 = mul_1(H_1, make_float3 (_S52, _S53, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / ((F32_max((rgi_out_0.z), (0.0f))) + 0.00000999999974738f));
    float _S54 = rgi_out_1.x;
    float _S55 = rgi_out_1.y;
    float3  _S56 = clamp_1(make_float3 (_S54, _S55, rgi_out_1.z - _S54 - _S55), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S57 = _S56.x;
    float _S58 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(0)].toe_0))))));
    float _S59 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(0)].shoulder_0))))));
    float _S60 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(0)].gamma_0))))));
    float _S61 = 1.0f / (1.0f + (F32_exp((- _S45.crf_params_1[int(0)].center_0))));
    float a_1 = _S59 * _S61 / lerp_0(_S58, _S59, _S61);
    float b_2 = 1.0f - a_1;
    float y_4;
    if(_S57 <= _S61)
    {
        y_4 = a_1 * (F32_pow((_S57 / _S61), (_S58)));
    }
    else
    {
        y_4 = 1.0f - b_2 * (F32_pow(((1.0f - _S57) / (1.0f - _S61)), (_S59)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S60)));
    float _S62 = _S56.y;
    float _S63 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(1)].toe_0))))));
    float _S64 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(1)].shoulder_0))))));
    float _S65 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(1)].gamma_0))))));
    float _S66 = 1.0f / (1.0f + (F32_exp((- _S45.crf_params_1[int(1)].center_0))));
    float a_2 = _S64 * _S66 / lerp_0(_S63, _S64, _S66);
    float b_3 = 1.0f - a_2;
    if(_S62 <= _S66)
    {
        y_4 = a_2 * (F32_pow((_S62 / _S66), (_S63)));
    }
    else
    {
        y_4 = 1.0f - b_3 * (F32_pow(((1.0f - _S62) / (1.0f - _S66)), (_S64)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S65)));
    float _S67 = _S56.z;
    float _S68 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(2)].toe_0))))));
    float _S69 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(2)].shoulder_0))))));
    float _S70 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S45.crf_params_1[int(2)].gamma_0))))));
    float _S71 = 1.0f / (1.0f + (F32_exp((- _S45.crf_params_1[int(2)].center_0))));
    float a_3 = _S69 * _S71 / lerp_0(_S68, _S69, _S71);
    float b_4 = 1.0f - a_3;
    if(_S67 <= _S71)
    {
        y_4 = a_3 * (F32_pow((_S67 / _S71), (_S68)));
    }
    else
    {
        y_4 = 1.0f - b_4 * (F32_pow(((1.0f - _S67) / (1.0f - _S71)), (_S69)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S70)));
    return rgb_out_1;
}

inline __device__ float3  apply_ppisp_rqs(float3  rgb_in_1, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, FixedArray<float, 39>  params_1)
{
    PPISPParamsRQS_0 p_1;
    (&p_1)->exposure_1 = params_1[int(0)];
    (&(&p_1)->vignette_params_1[int(0)])->cx_0 = params_1[int(1)];
    (&(&p_1)->vignette_params_1[int(0)])->cy_0 = params_1[int(2)];
    (&(&p_1)->vignette_params_1[int(0)])->alpha0_0 = params_1[int(3)];
    (&(&p_1)->vignette_params_1[int(0)])->alpha1_0 = params_1[int(4)];
    (&(&p_1)->vignette_params_1[int(0)])->alpha2_0 = params_1[int(5)];
    (&(&p_1)->vignette_params_1[int(1)])->cx_0 = params_1[int(6)];
    (&(&p_1)->vignette_params_1[int(1)])->cy_0 = params_1[int(7)];
    (&(&p_1)->vignette_params_1[int(1)])->alpha0_0 = params_1[int(8)];
    (&(&p_1)->vignette_params_1[int(1)])->alpha1_0 = params_1[int(9)];
    (&(&p_1)->vignette_params_1[int(1)])->alpha2_0 = params_1[int(10)];
    (&(&p_1)->vignette_params_1[int(2)])->cx_0 = params_1[int(11)];
    (&(&p_1)->vignette_params_1[int(2)])->cy_0 = params_1[int(12)];
    (&(&p_1)->vignette_params_1[int(2)])->alpha0_0 = params_1[int(13)];
    (&(&p_1)->vignette_params_1[int(2)])->alpha1_0 = params_1[int(14)];
    (&(&p_1)->vignette_params_1[int(2)])->alpha2_0 = params_1[int(15)];
    *&((&(&(&p_1)->color_params_1)->b_0)->x) = params_1[int(16)];
    *&((&(&(&p_1)->color_params_1)->b_0)->y) = params_1[int(17)];
    *&((&(&(&p_1)->color_params_1)->r_0)->x) = params_1[int(18)];
    *&((&(&(&p_1)->color_params_1)->r_0)->y) = params_1[int(19)];
    *&((&(&(&p_1)->color_params_1)->g_0)->x) = params_1[int(20)];
    *&((&(&(&p_1)->color_params_1)->g_0)->y) = params_1[int(21)];
    *&((&(&(&p_1)->color_params_1)->n_0)->x) = params_1[int(22)];
    *&((&(&(&p_1)->color_params_1)->n_0)->y) = params_1[int(23)];
    (&(&p_1)->crf_params_0[int(0)])->g0_0 = params_1[int(24)];
    (&(&p_1)->crf_params_0[int(0)])->g1_0 = params_1[int(25)];
    (&(&p_1)->crf_params_0[int(0)])->x0_0 = params_1[int(26)];
    (&(&p_1)->crf_params_0[int(0)])->y0_0 = params_1[int(27)];
    (&(&p_1)->crf_params_0[int(0)])->gc_0 = params_1[int(28)];
    (&(&p_1)->crf_params_0[int(1)])->g0_0 = params_1[int(29)];
    (&(&p_1)->crf_params_0[int(1)])->g1_0 = params_1[int(30)];
    (&(&p_1)->crf_params_0[int(1)])->x0_0 = params_1[int(31)];
    (&(&p_1)->crf_params_0[int(1)])->y0_0 = params_1[int(32)];
    (&(&p_1)->crf_params_0[int(1)])->gc_0 = params_1[int(33)];
    (&(&p_1)->crf_params_0[int(2)])->g0_0 = params_1[int(34)];
    (&(&p_1)->crf_params_0[int(2)])->g1_0 = params_1[int(35)];
    (&(&p_1)->crf_params_0[int(2)])->x0_0 = params_1[int(36)];
    (&(&p_1)->crf_params_0[int(2)])->y0_0 = params_1[int(37)];
    (&(&p_1)->crf_params_0[int(2)])->gc_0 = params_1[int(38)];
    PPISPParamsRQS_0 _S72 = p_1;
    float _S73 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S74 = (pix_coord_1.x - image_center_1.x) / _S73;
    float _S75 = (pix_coord_1.y - image_center_1.y) / _S73;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_1))));
    float dx_3 = _S74 - p_1.vignette_params_1[int(0)].cx_0;
    float dy_3 = _S75 - p_1.vignette_params_1[int(0)].cy_0;
    float r2_4 = dx_3 * dx_3 + dy_3 * dy_3;
    float r4_3 = r2_4 * r2_4;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_1[int(0)].alpha2_0 * (r4_3 * r2_4) + p_1.vignette_params_1[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_1[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S74 - p_1.vignette_params_1[int(1)].cx_0;
    float dy_4 = _S75 - p_1.vignette_params_1[int(1)].cy_0;
    float r2_5 = dx_4 * dx_4 + dy_4 * dy_4;
    float r4_4 = r2_5 * r2_5;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_1[int(1)].alpha2_0 * (r4_4 * r2_5) + p_1.vignette_params_1[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_1[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_5 = _S74 - p_1.vignette_params_1[int(2)].cx_0;
    float dy_5 = _S75 - p_1.vignette_params_1[int(2)].cy_0;
    float r2_6 = dx_5 * dx_5 + dy_5 * dy_5;
    float r4_5 = r2_6 * r2_6;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_1[int(2)].alpha2_0 * (r4_5 * r2_6) + p_1.vignette_params_1[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_1[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S76 = rgb_out_2;
    float2  bd_1 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_1.b_0);
    float2  rd_1 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_1.r_0);
    float2  gd_1 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_1.g_0);
    float2  nd_1 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_1.n_0);
    float _S77 = 0.3333333432674408f + nd_1.x;
    float _S78 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_1 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S78, 1.0f, 0.0f, - _S77, - _S78, _S77, 0.0f), T_1);
    float3  r0_1 = make_float3 (M_1.rows[int(0)].x, M_1.rows[int(0)].y, M_1.rows[int(0)].z);
    float3  r1_1 = make_float3 (M_1.rows[int(1)].x, M_1.rows[int(1)].y, M_1.rows[int(1)].z);
    float3  r2_7 = make_float3 (M_1.rows[int(2)].x, M_1.rows[int(2)].y, M_1.rows[int(2)].z);
    float3  lambda_v_3 = cross_0(r0_1, r1_1);
    float3  lambda_v_4;
    if((dot_0(lambda_v_3, lambda_v_3)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_5 = cross_0(r0_1, r2_7);
        if((dot_0(lambda_v_5, lambda_v_5)) < 9.99999968265522539e-21f)
        {
            lambda_v_4 = cross_0(r1_1, r2_7);
        }
        else
        {
            lambda_v_4 = lambda_v_5;
        }
    }
    else
    {
        lambda_v_4 = lambda_v_3;
    }
    Matrix<float, 3, 3>  H_2 = mul_3(mul_3(T_1, makeMatrix<float, 3, 3> (lambda_v_4.x, 0.0f, 0.0f, 0.0f, lambda_v_4.y, 0.0f, 0.0f, 0.0f, lambda_v_4.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_3;
    if((F32_abs((H_2.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_3 = H_2 * makeMatrix<float, 3, 3> (1.0f / H_2.rows[int(2)].z);
    }
    else
    {
        H_3 = H_2;
    }
    float _S79 = _S76.x;
    float _S80 = _S76.y;
    float intensity_1 = _S79 + _S80 + _S76.z;
    float3  rgi_out_2 = mul_1(H_3, make_float3 (_S79, _S80, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / ((F32_max((rgi_out_2.z), (0.0f))) + 0.00000999999974738f));
    float _S81 = rgi_out_3.x;
    float _S82 = rgi_out_3.y;
    float3  _S83 = clamp_1(make_float3 (_S81, _S82, rgi_out_3.z - _S81 - _S82), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S84 = _S83.x;
    float g0_1 = (F32_exp((_S72.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S72.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S72.crf_params_0[int(0)].gc_0)));
    float y_5;
    if(_S84 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S84 / x0_1;
        float _S85 = 1.0f - t0_0;
        y_5 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S85) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S85);
    }
    else
    {
        float _S86 = 1.0f - y0_1;
        float _S87 = 1.0f - x0_1;
        float s1_0 = _S86 / _S87;
        float t1_0 = (_S84 - x0_1) / _S87;
        float _S88 = 1.0f - t1_0;
        y_5 = y0_1 + _S86 * (s1_0 * t1_0 * t1_0 + gc_1 * t1_0 * _S88) / (s1_0 + (gc_1 + g1_1 - 2.0f * s1_0) * t1_0 * _S88);
    }
    *&((&rgb_out_3)->x) = y_5;
    float _S89 = _S83.y;
    float g0_2 = (F32_exp((_S72.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S72.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S72.crf_params_0[int(1)].gc_0)));
    if(_S89 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S89 / x0_2;
        float _S90 = 1.0f - t0_1;
        y_5 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S90) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S90);
    }
    else
    {
        float _S91 = 1.0f - y0_2;
        float _S92 = 1.0f - x0_2;
        float s1_1 = _S91 / _S92;
        float t1_1 = (_S89 - x0_2) / _S92;
        float _S93 = 1.0f - t1_1;
        y_5 = y0_2 + _S91 * (s1_1 * t1_1 * t1_1 + gc_2 * t1_1 * _S93) / (s1_1 + (gc_2 + g1_2 - 2.0f * s1_1) * t1_1 * _S93);
    }
    *&((&rgb_out_3)->y) = y_5;
    float _S94 = _S83.z;
    float g0_3 = (F32_exp((_S72.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S72.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S72.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S72.crf_params_0[int(2)].gc_0)));
    if(_S94 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S94 / x0_3;
        float _S95 = 1.0f - t0_2;
        y_5 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S95) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S95);
    }
    else
    {
        float _S96 = 1.0f - y0_3;
        float _S97 = 1.0f - x0_3;
        float s1_2 = _S96 / _S97;
        float t1_2 = (_S94 - x0_3) / _S97;
        float _S98 = 1.0f - t1_2;
        y_5 = y0_3 + _S96 * (s1_2 * t1_2 * t1_2 + gc_3 * t1_2 * _S98) / (s1_2 + (gc_3 + g1_3 - 2.0f * s1_2) * t1_2 * _S98);
    }
    *&((&rgb_out_3)->z) = y_5;
    return rgb_out_3;
}

inline __device__ float3  apply_ppisp_no_crf(float3  rgb_in_2, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, FixedArray<float, 24>  params_2)
{
    PPISPParamsNoCRF_0 p_2;
    (&p_2)->exposure_0 = params_2[int(0)];
    (&(&p_2)->vignette_params_0[int(0)])->cx_0 = params_2[int(1)];
    (&(&p_2)->vignette_params_0[int(0)])->cy_0 = params_2[int(2)];
    (&(&p_2)->vignette_params_0[int(0)])->alpha0_0 = params_2[int(3)];
    (&(&p_2)->vignette_params_0[int(0)])->alpha1_0 = params_2[int(4)];
    (&(&p_2)->vignette_params_0[int(0)])->alpha2_0 = params_2[int(5)];
    (&(&p_2)->vignette_params_0[int(1)])->cx_0 = params_2[int(6)];
    (&(&p_2)->vignette_params_0[int(1)])->cy_0 = params_2[int(7)];
    (&(&p_2)->vignette_params_0[int(1)])->alpha0_0 = params_2[int(8)];
    (&(&p_2)->vignette_params_0[int(1)])->alpha1_0 = params_2[int(9)];
    (&(&p_2)->vignette_params_0[int(1)])->alpha2_0 = params_2[int(10)];
    (&(&p_2)->vignette_params_0[int(2)])->cx_0 = params_2[int(11)];
    (&(&p_2)->vignette_params_0[int(2)])->cy_0 = params_2[int(12)];
    (&(&p_2)->vignette_params_0[int(2)])->alpha0_0 = params_2[int(13)];
    (&(&p_2)->vignette_params_0[int(2)])->alpha1_0 = params_2[int(14)];
    (&(&p_2)->vignette_params_0[int(2)])->alpha2_0 = params_2[int(15)];
    *&((&(&(&p_2)->color_params_0)->b_0)->x) = params_2[int(16)];
    *&((&(&(&p_2)->color_params_0)->b_0)->y) = params_2[int(17)];
    *&((&(&(&p_2)->color_params_0)->r_0)->x) = params_2[int(18)];
    *&((&(&(&p_2)->color_params_0)->r_0)->y) = params_2[int(19)];
    *&((&(&(&p_2)->color_params_0)->g_0)->x) = params_2[int(20)];
    *&((&(&(&p_2)->color_params_0)->g_0)->y) = params_2[int(21)];
    *&((&(&(&p_2)->color_params_0)->n_0)->x) = params_2[int(22)];
    *&((&(&(&p_2)->color_params_0)->n_0)->y) = params_2[int(23)];
    float _S99 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S100 = (pix_coord_2.x - image_center_2.x) / _S99;
    float _S101 = (pix_coord_2.y - image_center_2.y) / _S99;
    float3  rgb_out_4 = rgb_in_2 * make_float3 ((F32_exp2((p_2.exposure_0))));
    float dx_6 = _S100 - p_2.vignette_params_0[int(0)].cx_0;
    float dy_6 = _S101 - p_2.vignette_params_0[int(0)].cy_0;
    float r2_8 = dx_6 * dx_6 + dy_6 * dy_6;
    float r4_6 = r2_8 * r2_8;
    *&((&rgb_out_4)->x) = *&((&rgb_out_4)->x) * clamp_0(p_2.vignette_params_0[int(0)].alpha2_0 * (r4_6 * r2_8) + p_2.vignette_params_0[int(0)].alpha1_0 * r4_6 + p_2.vignette_params_0[int(0)].alpha0_0 * r2_8 + 1.0f, 0.0f, 1.0f);
    float dx_7 = _S100 - p_2.vignette_params_0[int(1)].cx_0;
    float dy_7 = _S101 - p_2.vignette_params_0[int(1)].cy_0;
    float r2_9 = dx_7 * dx_7 + dy_7 * dy_7;
    float r4_7 = r2_9 * r2_9;
    *&((&rgb_out_4)->y) = *&((&rgb_out_4)->y) * clamp_0(p_2.vignette_params_0[int(1)].alpha2_0 * (r4_7 * r2_9) + p_2.vignette_params_0[int(1)].alpha1_0 * r4_7 + p_2.vignette_params_0[int(1)].alpha0_0 * r2_9 + 1.0f, 0.0f, 1.0f);
    float dx_8 = _S100 - p_2.vignette_params_0[int(2)].cx_0;
    float dy_8 = _S101 - p_2.vignette_params_0[int(2)].cy_0;
    float r2_10 = dx_8 * dx_8 + dy_8 * dy_8;
    float r4_8 = r2_10 * r2_10;
    *&((&rgb_out_4)->z) = *&((&rgb_out_4)->z) * clamp_0(p_2.vignette_params_0[int(2)].alpha2_0 * (r4_8 * r2_10) + p_2.vignette_params_0[int(2)].alpha1_0 * r4_8 + p_2.vignette_params_0[int(2)].alpha0_0 * r2_10 + 1.0f, 0.0f, 1.0f);
    float3  _S102 = rgb_out_4;
    float2  bd_2 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_2.color_params_0.b_0);
    float2  rd_2 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_2.color_params_0.r_0);
    float2  gd_2 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_2.color_params_0.g_0);
    float2  nd_2 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_2.color_params_0.n_0);
    float _S103 = 0.3333333432674408f + nd_2.x;
    float _S104 = 0.3333333432674408f + nd_2.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (bd_2.x, 1.0f + rd_2.x, gd_2.x, bd_2.y, rd_2.y, 1.0f + gd_2.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_2 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S104, 1.0f, 0.0f, - _S103, - _S104, _S103, 0.0f), T_2);
    float3  r0_2 = make_float3 (M_2.rows[int(0)].x, M_2.rows[int(0)].y, M_2.rows[int(0)].z);
    float3  r1_2 = make_float3 (M_2.rows[int(1)].x, M_2.rows[int(1)].y, M_2.rows[int(1)].z);
    float3  r2_11 = make_float3 (M_2.rows[int(2)].x, M_2.rows[int(2)].y, M_2.rows[int(2)].z);
    float3  lambda_v_6 = cross_0(r0_2, r1_2);
    float3  lambda_v_7;
    if((dot_0(lambda_v_6, lambda_v_6)) < 9.99999968265522539e-21f)
    {
        float3  lambda_v_8 = cross_0(r0_2, r2_11);
        if((dot_0(lambda_v_8, lambda_v_8)) < 9.99999968265522539e-21f)
        {
            lambda_v_7 = cross_0(r1_2, r2_11);
        }
        else
        {
            lambda_v_7 = lambda_v_8;
        }
    }
    else
    {
        lambda_v_7 = lambda_v_6;
    }
    Matrix<float, 3, 3>  H_4 = mul_3(mul_3(T_2, makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z)), makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f));
    Matrix<float, 3, 3>  H_5;
    if((F32_abs((H_4.rows[int(2)].z))) > 9.99999968265522539e-21f)
    {
        H_5 = H_4 * makeMatrix<float, 3, 3> (1.0f / H_4.rows[int(2)].z);
    }
    else
    {
        H_5 = H_4;
    }
    float _S105 = _S102.x;
    float _S106 = _S102.y;
    float intensity_2 = _S105 + _S106 + _S102.z;
    float3  rgi_out_4 = mul_1(H_5, make_float3 (_S105, _S106, intensity_2));
    float3  rgi_out_5 = rgi_out_4 * make_float3 (intensity_2 / ((F32_max((rgi_out_4.z), (0.0f))) + 0.00000999999974738f));
    float _S107 = rgi_out_5.x;
    float _S108 = rgi_out_5.y;
    return clamp_1(make_float3 (_S107, _S108, rgi_out_5.z - _S107 - _S108), make_float3 (0.0f), make_float3 (1.0f));
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S109)
{
    return (F32_exp2((_S109)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S110, float _S111, float _S112)
{
    return clamp_0(_S110, _S111, _S112);
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S113, float2  _S114)
{
    return mul_0(_S113, _S114);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S115, Matrix<float, 3, 3>  _S116)
{
    return mul_3(_S115, _S116);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S117, float3  _S118)
{
    return cross_0(_S117, _S118);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S119, float3  _S120)
{
    return dot_0(_S119, _S120);
}

inline __device__ float s_primal_ctx_abs_0(float _S121)
{
    return (F32_abs((_S121)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S122, float3  _S123)
{
    return mul_1(_S122, _S123);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S124, float3  _S125, float3  _S126)
{
    return clamp_1(_S124, _S125, _S126);
}

inline __device__ float s_primal_ctx_exp_0(float _S127)
{
    return (F32_exp((_S127)));
}

inline __device__ float s_primal_ctx_log_0(float _S128)
{
    return (F32_log((_S128)));
}

inline __device__ float s_primal_ctx_lerp_0(float _S129, float _S130, float _S131)
{
    return lerp_0(_S129, _S130, _S131);
}

inline __device__ float s_primal_ctx_pow_0(float _S132, float _S133)
{
    return (F32_pow((_S132), (_S133)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S134, DiffPair_float_0 * _S135, float _S136)
{
    _d_pow_0(_S134, _S135, _S136);
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S137, DiffPair_float_0 * _S138, DiffPair_float_0 * _S139, float _S140)
{
    _d_lerp_0(_S137, _S138, _S139, _S140);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S141, float _S142)
{
    _d_exp_0(_S141, _S142);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S143, float _S144)
{
    _d_log_0(_S143, _S144);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S145, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S146, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S147, float3  _S148)
{
    _d_clamp_vector_0(_S145, _S146, _S147, _S148);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S149, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S150, float3  _S151)
{
    _d_mul_1(_S149, _S150, _S151);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S152, float _S153)
{
    _d_abs_0(_S152, _S153);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S154, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S155, Matrix<float, 3, 3>  _S156)
{
    mul_2(_S154, _S155, _S156);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S157, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S158, float3  _S159)
{
    _d_cross_0(_S157, _S158, _S159);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S160, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S161, float _S162)
{
    _d_dot_0(_S160, _S161, _S162);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S163, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S164, float2  _S165)
{
    _d_mul_0(_S163, _S164, _S165);
    return;
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_float_0 * _S166, DiffPair_float_0 * _S167, DiffPair_float_0 * _S168, float _S169)
{
    _d_clamp_0(_S166, _S167, _S168, _S169);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S170, float _S171)
{
    _d_exp2_0(_S170, _S171);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_0)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S172 = *dprgb_in_0;
    float3  _S173 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S174 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S175 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S176 = {
        _S175, _S175, _S175
    };
    float2  _S177 = make_float2 (0.0f);
    ColorPPISPParams_0 _S178 = { _S177, _S177, _S177, _S177 };
    CRFPPISPChannelParams_0 _S179 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S180 = {
        _S179, _S179, _S179
    };
    PPISPParams_0 _S181;
    (&_S181)->exposure_2 = dpparams_0->primal_0[int(0)];
    (&_S181)->vignette_params_2 = _S176;
    (&_S181)->color_params_2 = _S178;
    (&_S181)->crf_params_1 = _S180;
    (&(&_S181)->vignette_params_2[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S181)->vignette_params_2[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S182 = dpparams_0->primal_0[int(3)];
    (&(&_S181)->vignette_params_2[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S183 = dpparams_0->primal_0[int(4)];
    (&(&_S181)->vignette_params_2[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S184 = dpparams_0->primal_0[int(5)];
    (&(&_S181)->vignette_params_2[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S181)->vignette_params_2[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S181)->vignette_params_2[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S185 = dpparams_0->primal_0[int(8)];
    (&(&_S181)->vignette_params_2[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S186 = dpparams_0->primal_0[int(9)];
    (&(&_S181)->vignette_params_2[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S187 = dpparams_0->primal_0[int(10)];
    (&(&_S181)->vignette_params_2[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S181)->vignette_params_2[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S181)->vignette_params_2[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S188 = dpparams_0->primal_0[int(13)];
    (&(&_S181)->vignette_params_2[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S189 = dpparams_0->primal_0[int(14)];
    (&(&_S181)->vignette_params_2[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S190 = dpparams_0->primal_0[int(15)];
    (&(&_S181)->vignette_params_2[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S181)->color_params_2)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S181)->color_params_2)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S181)->color_params_2)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S181)->color_params_2)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S181)->color_params_2)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S181)->color_params_2)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S181)->color_params_2)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S181)->color_params_2)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S191 = dpparams_0->primal_0[int(24)];
    (&(&_S181)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S192 = dpparams_0->primal_0[int(25)];
    (&(&_S181)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S193 = dpparams_0->primal_0[int(26)];
    (&(&_S181)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S194 = dpparams_0->primal_0[int(27)];
    (&(&_S181)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S195 = dpparams_0->primal_0[int(28)];
    (&(&_S181)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S196 = dpparams_0->primal_0[int(29)];
    (&(&_S181)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S197 = dpparams_0->primal_0[int(30)];
    (&(&_S181)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S198 = dpparams_0->primal_0[int(31)];
    (&(&_S181)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S199 = dpparams_0->primal_0[int(32)];
    (&(&_S181)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S200 = dpparams_0->primal_0[int(33)];
    (&(&_S181)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S201 = dpparams_0->primal_0[int(34)];
    (&(&_S181)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S202 = dpparams_0->primal_0[int(35)];
    (&(&_S181)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S203 = _S181;
    float _S204 = s_primal_ctx_exp2_0(_S181.exposure_2);
    float3  _S205 = make_float3 (_S204);
    float3  rgb_out_5 = (*dprgb_in_0).primal_0 * make_float3 (_S204);
    float _S206 = (F32_max((img_size_3.x), (img_size_3.y)));
    float _S207 = (pix_coord_3.x - image_center_3.x) / _S206;
    float _S208 = (pix_coord_3.y - image_center_3.y) / _S206;
    float dx_9 = _S207 - dpparams_0->primal_0[int(1)];
    float dy_9 = _S208 - dpparams_0->primal_0[int(2)];
    float r2_12 = dx_9 * dx_9 + dy_9 * dy_9;
    float r4_9 = r2_12 * r2_12;
    float r6_0 = r4_9 * r2_12;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_9 + dpparams_0->primal_0[int(3)] * r2_12 + 1.0f;
    float _S209 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S210 = rgb_out_5.x * _S209;
    float3  _S211 = rgb_out_5;
    *&((&_S211)->x) = _S210;
    float dx_10 = _S207 - dpparams_0->primal_0[int(6)];
    float dy_10 = _S208 - dpparams_0->primal_0[int(7)];
    float r2_13 = dx_10 * dx_10 + dy_10 * dy_10;
    float r4_10 = r2_13 * r2_13;
    float r6_1 = r4_10 * r2_13;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_10 + dpparams_0->primal_0[int(8)] * r2_13 + 1.0f;
    float _S212 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S211)->y) = rgb_out_5.y * _S212;
    float dx_11 = _S207 - dpparams_0->primal_0[int(11)];
    float dy_11 = _S208 - dpparams_0->primal_0[int(12)];
    float r2_14 = dx_11 * dx_11 + dy_11 * dy_11;
    float r4_11 = r2_14 * r2_14;
    float r6_2 = r4_11 * r2_14;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_11 + dpparams_0->primal_0[int(13)] * r2_14 + 1.0f;
    float _S213 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S211)->z) = rgb_out_5.z * _S213;
    PPISPParams_0 _S214 = _S181;
    float2  _S215 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S181.color_params_2.b_0);
    float2  _S216 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S181.color_params_2.r_0);
    float2  _S217 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S181.color_params_2.g_0);
    float2  _S218 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S181.color_params_2.n_0);
    float _S219 = 0.3333333432674408f + _S218.x;
    float _S220 = 0.3333333432674408f + _S218.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S215.x, 1.0f + _S216.x, _S217.x, _S215.y, _S216.y, 1.0f + _S217.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S220, 1.0f, 0.0f, - _S219, - _S220, _S219, 0.0f);
    Matrix<float, 3, 3>  _S221 = s_primal_ctx_mul_1(skew_0, T_3);
    float3  r0_3 = make_float3 (_S221.rows[int(0)].x, _S221.rows[int(0)].y, _S221.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S221.rows[int(1)].x, _S221.rows[int(1)].y, _S221.rows[int(1)].z);
    float3  r2_15 = make_float3 (_S221.rows[int(2)].x, _S221.rows[int(2)].y, _S221.rows[int(2)].z);
    float3  _S222 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S223 = (s_primal_ctx_dot_0(_S222, _S222)) < 9.99999968265522539e-21f;
    float3  lambda_v_9;
    float3  _S224;
    bool _S225;
    if(_S223)
    {
        float3  _S226 = s_primal_ctx_cross_0(r0_3, r2_15);
        bool _S227 = (s_primal_ctx_dot_0(_S226, _S226)) < 9.99999968265522539e-21f;
        if(_S227)
        {
            lambda_v_9 = s_primal_ctx_cross_0(r1_3, r2_15);
        }
        else
        {
            lambda_v_9 = _S226;
        }
        _S225 = _S227;
        _S224 = _S226;
    }
    else
    {
        lambda_v_9 = _S222;
        _S225 = false;
        _S224 = _S173;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_9.x, 0.0f, 0.0f, 0.0f, lambda_v_9.y, 0.0f, 0.0f, 0.0f, lambda_v_9.z);
    Matrix<float, 3, 3>  _S228 = s_primal_ctx_mul_1(T_3, D_0);
    Matrix<float, 3, 3>  _S229 = s_primal_ctx_mul_1(_S228, S_inv_0);
    bool _S230 = (s_primal_ctx_abs_0(_S229.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_6;
    Matrix<float, 3, 3>  _S231;
    float _S232;
    if(_S230)
    {
        float inv_s_0 = 1.0f / _S229.rows[int(2)].z;
        Matrix<float, 3, 3>  _S233 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S234 = _S229.rows[int(2)].z * _S229.rows[int(2)].z;
        H_6 = _S229 * makeMatrix<float, 3, 3> (inv_s_0);
        _S231 = _S233;
        _S232 = _S234;
    }
    else
    {
        H_6 = _S229;
        _S231 = _S174;
        _S232 = 0.0f;
    }
    float _S235 = _S211.x;
    float _S236 = _S211.y;
    float intensity_3 = _S235 + _S236 + _S211.z;
    float3  rgi_in_0 = make_float3 (_S235, _S236, intensity_3);
    float3  _S237 = s_primal_ctx_mul_2(H_6, rgi_in_0);
    float _S238 = _S237.z;
    float _S239 = (F32_max((_S238), (0.0f))) + 0.00000999999974738f;
    float norm_factor_0 = intensity_3 / _S239;
    float3  _S240 = make_float3 (norm_factor_0);
    float _S241 = _S239 * _S239;
    float3  rgi_out_6 = _S237 * make_float3 (norm_factor_0);
    float _S242 = rgi_out_6.x;
    float _S243 = rgi_out_6.y;
    float3  _S244 = make_float3 (_S242, _S243, rgi_out_6.z - _S242 - _S243);
    float3  _S245 = make_float3 (0.0f);
    float3  _S246 = make_float3 (1.0f);
    float3  _S247 = s_primal_ctx_clamp_1(_S244, _S245, _S246);
    float _S248 = _S247.x;
    float _S249 = 1.0f + s_primal_ctx_exp_0(_S191);
    float _S250 = 0.30000001192092896f + s_primal_ctx_log_0(_S249);
    float _S251 = 1.0f + s_primal_ctx_exp_0(_S192);
    float _S252 = 0.30000001192092896f + s_primal_ctx_log_0(_S251);
    float _S253 = 1.0f + s_primal_ctx_exp_0(_S193);
    float _S254 = 0.10000000149011612f + s_primal_ctx_log_0(_S253);
    float _S255 = - _S194;
    float _S256 = 1.0f + s_primal_ctx_exp_0(_S255);
    float _S257 = 1.0f / _S256;
    float _S258 = _S256 * _S256;
    float _S259 = s_primal_ctx_lerp_0(_S250, _S252, _S257);
    float _S260 = _S252 * _S257;
    float a_4 = _S260 / _S259;
    float _S261 = _S259 * _S259;
    float b_5 = 1.0f - a_4;
    bool _S262 = _S248 <= _S257;
    float y_6;
    float _S263;
    float _S264;
    float _S265;
    float _S266;
    float _S267;
    float _S268;
    float _S269;
    float _S270;
    if(_S262)
    {
        float _S271 = _S248 / _S257;
        float _S272 = _S257 * _S257;
        float _S273 = s_primal_ctx_pow_0(_S271, _S250);
        y_6 = a_4 * _S273;
        _S263 = _S273;
        _S264 = _S271;
        _S265 = _S272;
        _S266 = 0.0f;
        _S267 = 0.0f;
        _S268 = 0.0f;
        _S269 = 0.0f;
        _S270 = 0.0f;
    }
    else
    {
        float _S274 = 1.0f - _S248;
        float _S275 = 1.0f - _S257;
        float _S276 = _S274 / _S275;
        float _S277 = _S275 * _S275;
        float _S278 = s_primal_ctx_pow_0(_S276, _S252);
        y_6 = 1.0f - b_5 * _S278;
        _S263 = 0.0f;
        _S264 = 0.0f;
        _S265 = 0.0f;
        _S266 = _S278;
        _S267 = _S276;
        _S268 = _S277;
        _S269 = _S274;
        _S270 = _S275;
    }
    float _S279 = (F32_max((0.0f), (y_6)));
    float _S280 = _S247.y;
    float _S281 = 1.0f + s_primal_ctx_exp_0(_S195);
    float _S282 = 0.30000001192092896f + s_primal_ctx_log_0(_S281);
    float _S283 = 1.0f + s_primal_ctx_exp_0(_S196);
    float _S284 = 0.30000001192092896f + s_primal_ctx_log_0(_S283);
    float _S285 = 1.0f + s_primal_ctx_exp_0(_S197);
    float _S286 = 0.10000000149011612f + s_primal_ctx_log_0(_S285);
    float _S287 = - _S198;
    float _S288 = 1.0f + s_primal_ctx_exp_0(_S287);
    float _S289 = 1.0f / _S288;
    float _S290 = _S288 * _S288;
    float _S291 = s_primal_ctx_lerp_0(_S282, _S284, _S289);
    float _S292 = _S284 * _S289;
    float a_5 = _S292 / _S291;
    float _S293 = _S291 * _S291;
    float b_6 = 1.0f - a_5;
    bool _S294 = _S280 <= _S289;
    float y_7;
    float _S295;
    float _S296;
    float _S297;
    float _S298;
    float _S299;
    float _S300;
    float _S301;
    float _S302;
    if(_S294)
    {
        float _S303 = _S280 / _S289;
        float _S304 = _S289 * _S289;
        float _S305 = s_primal_ctx_pow_0(_S303, _S282);
        y_7 = a_5 * _S305;
        _S295 = _S305;
        _S296 = _S303;
        _S297 = _S304;
        _S298 = 0.0f;
        _S299 = 0.0f;
        _S300 = 0.0f;
        _S301 = 0.0f;
        _S302 = 0.0f;
    }
    else
    {
        float _S306 = 1.0f - _S280;
        float _S307 = 1.0f - _S289;
        float _S308 = _S306 / _S307;
        float _S309 = _S307 * _S307;
        float _S310 = s_primal_ctx_pow_0(_S308, _S284);
        y_7 = 1.0f - b_6 * _S310;
        _S295 = 0.0f;
        _S296 = 0.0f;
        _S297 = 0.0f;
        _S298 = _S310;
        _S299 = _S308;
        _S300 = _S309;
        _S301 = _S306;
        _S302 = _S307;
    }
    float _S311 = (F32_max((0.0f), (y_7)));
    float _S312 = _S247.z;
    float _S313 = 1.0f + s_primal_ctx_exp_0(_S199);
    float _S314 = 0.30000001192092896f + s_primal_ctx_log_0(_S313);
    float _S315 = 1.0f + s_primal_ctx_exp_0(_S200);
    float _S316 = 0.30000001192092896f + s_primal_ctx_log_0(_S315);
    float _S317 = 1.0f + s_primal_ctx_exp_0(_S201);
    float _S318 = 0.10000000149011612f + s_primal_ctx_log_0(_S317);
    float _S319 = - _S202;
    float _S320 = 1.0f + s_primal_ctx_exp_0(_S319);
    float _S321 = 1.0f / _S320;
    float _S322 = _S320 * _S320;
    float _S323 = s_primal_ctx_lerp_0(_S314, _S316, _S321);
    float _S324 = _S316 * _S321;
    float a_6 = _S324 / _S323;
    float _S325 = _S323 * _S323;
    float b_7 = 1.0f - a_6;
    bool _S326 = _S312 <= _S321;
    float y_8;
    float _S327;
    float _S328;
    float _S329;
    float _S330;
    float _S331;
    float _S332;
    float _S333;
    float _S334;
    if(_S326)
    {
        float _S335 = _S312 / _S321;
        float _S336 = _S321 * _S321;
        float _S337 = s_primal_ctx_pow_0(_S335, _S314);
        y_8 = a_6 * _S337;
        _S327 = _S337;
        _S328 = _S335;
        _S329 = _S336;
        _S330 = 0.0f;
        _S331 = 0.0f;
        _S332 = 0.0f;
        _S333 = 0.0f;
        _S334 = 0.0f;
    }
    else
    {
        float _S338 = 1.0f - _S312;
        float _S339 = 1.0f - _S321;
        float _S340 = _S338 / _S339;
        float _S341 = _S339 * _S339;
        float _S342 = s_primal_ctx_pow_0(_S340, _S316);
        y_8 = 1.0f - b_7 * _S342;
        _S327 = 0.0f;
        _S328 = 0.0f;
        _S329 = 0.0f;
        _S330 = _S342;
        _S331 = _S340;
        _S332 = _S341;
        _S333 = _S338;
        _S334 = _S339;
    }
    float _S343 = (F32_max((0.0f), (y_8)));
    DiffPair_float_0 _S344;
    (&_S344)->primal_0 = _S343;
    (&_S344)->differential_0 = 0.0f;
    DiffPair_float_0 _S345;
    (&_S345)->primal_0 = _S318;
    (&_S345)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S344, &_S345, _s_dOut_0.z);
    DiffPair_float_0 _S346 = _S345;
    DiffPair_float_0 _S347;
    (&_S347)->primal_0 = 0.0f;
    (&_S347)->differential_0 = 0.0f;
    DiffPair_float_0 _S348;
    (&_S348)->primal_0 = y_8;
    (&_S348)->differential_0 = 0.0f;
    _d_max_0(&_S347, &_S348, _S344.differential_0);
    DiffPair_float_0 _S349 = _S348;
    if(_S326)
    {
        float _S350 = a_6 * _S349.differential_0;
        float _S351 = _S327 * _S349.differential_0;
        DiffPair_float_0 _S352;
        (&_S352)->primal_0 = _S328;
        (&_S352)->differential_0 = 0.0f;
        DiffPair_float_0 _S353;
        (&_S353)->primal_0 = _S314;
        (&_S353)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S352, &_S353, _S350);
        float _S354 = _S352.differential_0 / _S329;
        float _S355 = _S312 * - _S354;
        float _S356 = _S321 * _S354;
        y_8 = 0.0f;
        _S327 = _S351;
        _S328 = _S355;
        _S329 = 0.0f;
        _S330 = _S353.differential_0;
        _S331 = _S356;
    }
    else
    {
        float _S357 = - _S349.differential_0;
        float _S358 = b_7 * _S357;
        float _S359 = _S330 * _S357;
        DiffPair_float_0 _S360;
        (&_S360)->primal_0 = _S331;
        (&_S360)->differential_0 = 0.0f;
        DiffPair_float_0 _S361;
        (&_S361)->primal_0 = _S316;
        (&_S361)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S360, &_S361, _S358);
        float _S362 = _S360.differential_0 / _S332;
        float _S363 = - (_S333 * - _S362);
        float _S364 = - (_S334 * _S362);
        y_8 = _S359;
        _S327 = 0.0f;
        _S328 = _S363;
        _S329 = _S361.differential_0;
        _S330 = 0.0f;
        _S331 = _S364;
    }
    float _S365 = (- y_8 + _S327) / _S325;
    float _S366 = _S324 * - _S365;
    float _S367 = _S323 * _S365;
    float _S368 = _S316 * _S367;
    float _S369 = _S321 * _S367;
    DiffPair_float_0 _S370;
    (&_S370)->primal_0 = _S314;
    (&_S370)->differential_0 = 0.0f;
    DiffPair_float_0 _S371;
    (&_S371)->primal_0 = _S316;
    (&_S371)->differential_0 = 0.0f;
    DiffPair_float_0 _S372;
    (&_S372)->primal_0 = _S321;
    (&_S372)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S370, &_S371, &_S372, _S366);
    float _S373 = - ((_S368 + _S372.differential_0 + _S328) / _S322);
    DiffPair_float_0 _S374;
    (&_S374)->primal_0 = _S319;
    (&_S374)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S374, _S373);
    float _S375 = - _S374.differential_0;
    DiffPair_float_0 _S376;
    (&_S376)->primal_0 = _S317;
    (&_S376)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S376, _S346.differential_0);
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = _S201;
    (&_S377)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S377, _S376.differential_0);
    DiffPair_float_0 _S378 = _S377;
    float _S379 = _S369 + _S371.differential_0 + _S329;
    DiffPair_float_0 _S380;
    (&_S380)->primal_0 = _S315;
    (&_S380)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S380, _S379);
    DiffPair_float_0 _S381;
    (&_S381)->primal_0 = _S200;
    (&_S381)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S381, _S380.differential_0);
    DiffPair_float_0 _S382 = _S381;
    float _S383 = _S370.differential_0 + _S330;
    DiffPair_float_0 _S384;
    (&_S384)->primal_0 = _S313;
    (&_S384)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S384, _S383);
    DiffPair_float_0 _S385;
    (&_S385)->primal_0 = _S199;
    (&_S385)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S385, _S384.differential_0);
    DiffPair_float_0 _S386 = _S385;
    float3  _S387 = make_float3 (0.0f, 0.0f, _S331);
    DiffPair_float_0 _S388;
    (&_S388)->primal_0 = _S311;
    (&_S388)->differential_0 = 0.0f;
    DiffPair_float_0 _S389;
    (&_S389)->primal_0 = _S286;
    (&_S389)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S388, &_S389, _s_dOut_0.y);
    DiffPair_float_0 _S390 = _S389;
    DiffPair_float_0 _S391;
    (&_S391)->primal_0 = 0.0f;
    (&_S391)->differential_0 = 0.0f;
    DiffPair_float_0 _S392;
    (&_S392)->primal_0 = y_7;
    (&_S392)->differential_0 = 0.0f;
    _d_max_0(&_S391, &_S392, _S388.differential_0);
    DiffPair_float_0 _S393 = _S392;
    if(_S294)
    {
        float _S394 = a_5 * _S393.differential_0;
        float _S395 = _S295 * _S393.differential_0;
        DiffPair_float_0 _S396;
        (&_S396)->primal_0 = _S296;
        (&_S396)->differential_0 = 0.0f;
        DiffPair_float_0 _S397;
        (&_S397)->primal_0 = _S282;
        (&_S397)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S396, &_S397, _S394);
        float _S398 = _S396.differential_0 / _S297;
        float _S399 = _S280 * - _S398;
        float _S400 = _S289 * _S398;
        y_7 = 0.0f;
        _S295 = _S395;
        _S296 = _S399;
        _S297 = 0.0f;
        _S298 = _S397.differential_0;
        _S299 = _S400;
    }
    else
    {
        float _S401 = - _S393.differential_0;
        float _S402 = b_6 * _S401;
        float _S403 = _S298 * _S401;
        DiffPair_float_0 _S404;
        (&_S404)->primal_0 = _S299;
        (&_S404)->differential_0 = 0.0f;
        DiffPair_float_0 _S405;
        (&_S405)->primal_0 = _S284;
        (&_S405)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S404, &_S405, _S402);
        float _S406 = _S404.differential_0 / _S300;
        float _S407 = - (_S301 * - _S406);
        float _S408 = - (_S302 * _S406);
        y_7 = _S403;
        _S295 = 0.0f;
        _S296 = _S407;
        _S297 = _S405.differential_0;
        _S298 = 0.0f;
        _S299 = _S408;
    }
    float _S409 = (- y_7 + _S295) / _S293;
    float _S410 = _S292 * - _S409;
    float _S411 = _S291 * _S409;
    float _S412 = _S284 * _S411;
    float _S413 = _S289 * _S411;
    DiffPair_float_0 _S414;
    (&_S414)->primal_0 = _S282;
    (&_S414)->differential_0 = 0.0f;
    DiffPair_float_0 _S415;
    (&_S415)->primal_0 = _S284;
    (&_S415)->differential_0 = 0.0f;
    DiffPair_float_0 _S416;
    (&_S416)->primal_0 = _S289;
    (&_S416)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S414, &_S415, &_S416, _S410);
    float _S417 = - ((_S412 + _S416.differential_0 + _S296) / _S290);
    DiffPair_float_0 _S418;
    (&_S418)->primal_0 = _S287;
    (&_S418)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S418, _S417);
    float _S419 = - _S418.differential_0;
    DiffPair_float_0 _S420;
    (&_S420)->primal_0 = _S285;
    (&_S420)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S420, _S390.differential_0);
    DiffPair_float_0 _S421;
    (&_S421)->primal_0 = _S197;
    (&_S421)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S421, _S420.differential_0);
    DiffPair_float_0 _S422 = _S421;
    float _S423 = _S413 + _S415.differential_0 + _S297;
    DiffPair_float_0 _S424;
    (&_S424)->primal_0 = _S283;
    (&_S424)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S424, _S423);
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = _S196;
    (&_S425)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S425, _S424.differential_0);
    DiffPair_float_0 _S426 = _S425;
    float _S427 = _S414.differential_0 + _S298;
    DiffPair_float_0 _S428;
    (&_S428)->primal_0 = _S281;
    (&_S428)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S428, _S427);
    DiffPair_float_0 _S429;
    (&_S429)->primal_0 = _S195;
    (&_S429)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S429, _S428.differential_0);
    DiffPair_float_0 _S430 = _S429;
    float3  _S431 = _S387 + make_float3 (0.0f, _S299, 0.0f);
    DiffPair_float_0 _S432;
    (&_S432)->primal_0 = _S279;
    (&_S432)->differential_0 = 0.0f;
    DiffPair_float_0 _S433;
    (&_S433)->primal_0 = _S254;
    (&_S433)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S432, &_S433, _s_dOut_0.x);
    DiffPair_float_0 _S434 = _S433;
    DiffPair_float_0 _S435;
    (&_S435)->primal_0 = 0.0f;
    (&_S435)->differential_0 = 0.0f;
    DiffPair_float_0 _S436;
    (&_S436)->primal_0 = y_6;
    (&_S436)->differential_0 = 0.0f;
    _d_max_0(&_S435, &_S436, _S432.differential_0);
    DiffPair_float_0 _S437 = _S436;
    if(_S262)
    {
        float _S438 = a_4 * _S437.differential_0;
        float _S439 = _S263 * _S437.differential_0;
        DiffPair_float_0 _S440;
        (&_S440)->primal_0 = _S264;
        (&_S440)->differential_0 = 0.0f;
        DiffPair_float_0 _S441;
        (&_S441)->primal_0 = _S250;
        (&_S441)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S440, &_S441, _S438);
        float _S442 = _S440.differential_0 / _S265;
        float _S443 = _S248 * - _S442;
        float _S444 = _S257 * _S442;
        y_6 = 0.0f;
        _S263 = _S439;
        _S264 = _S443;
        _S265 = 0.0f;
        _S266 = _S441.differential_0;
        _S267 = _S444;
    }
    else
    {
        float _S445 = - _S437.differential_0;
        float _S446 = b_5 * _S445;
        float _S447 = _S266 * _S445;
        DiffPair_float_0 _S448;
        (&_S448)->primal_0 = _S267;
        (&_S448)->differential_0 = 0.0f;
        DiffPair_float_0 _S449;
        (&_S449)->primal_0 = _S252;
        (&_S449)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S448, &_S449, _S446);
        float _S450 = _S448.differential_0 / _S268;
        float _S451 = - (_S269 * - _S450);
        float _S452 = - (_S270 * _S450);
        y_6 = _S447;
        _S263 = 0.0f;
        _S264 = _S451;
        _S265 = _S449.differential_0;
        _S266 = 0.0f;
        _S267 = _S452;
    }
    float _S453 = (- y_6 + _S263) / _S261;
    float _S454 = _S260 * - _S453;
    float _S455 = _S259 * _S453;
    float _S456 = _S252 * _S455;
    float _S457 = _S257 * _S455;
    DiffPair_float_0 _S458;
    (&_S458)->primal_0 = _S250;
    (&_S458)->differential_0 = 0.0f;
    DiffPair_float_0 _S459;
    (&_S459)->primal_0 = _S252;
    (&_S459)->differential_0 = 0.0f;
    DiffPair_float_0 _S460;
    (&_S460)->primal_0 = _S257;
    (&_S460)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S458, &_S459, &_S460, _S454);
    float _S461 = - ((_S456 + _S460.differential_0 + _S264) / _S258);
    DiffPair_float_0 _S462;
    (&_S462)->primal_0 = _S255;
    (&_S462)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S462, _S461);
    float _S463 = - _S462.differential_0;
    DiffPair_float_0 _S464;
    (&_S464)->primal_0 = _S253;
    (&_S464)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S464, _S434.differential_0);
    DiffPair_float_0 _S465;
    (&_S465)->primal_0 = _S193;
    (&_S465)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S465, _S464.differential_0);
    DiffPair_float_0 _S466 = _S465;
    float _S467 = _S457 + _S459.differential_0 + _S265;
    DiffPair_float_0 _S468;
    (&_S468)->primal_0 = _S251;
    (&_S468)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S468, _S467);
    DiffPair_float_0 _S469;
    (&_S469)->primal_0 = _S192;
    (&_S469)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S469, _S468.differential_0);
    DiffPair_float_0 _S470 = _S469;
    float _S471 = _S458.differential_0 + _S266;
    DiffPair_float_0 _S472;
    (&_S472)->primal_0 = _S249;
    (&_S472)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S472, _S471);
    DiffPair_float_0 _S473;
    (&_S473)->primal_0 = _S191;
    (&_S473)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S473, _S472.differential_0);
    DiffPair_float_0 _S474 = _S473;
    float3  _S475 = _S431 + make_float3 (_S267, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S476;
    (&_S476)->primal_0 = _S244;
    (&_S476)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S477;
    (&_S477)->primal_0 = _S245;
    (&_S477)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S478;
    (&_S478)->primal_0 = _S246;
    (&_S478)->differential_0 = _S173;
    s_bwd_prop_clamp_0(&_S476, &_S477, &_S478, _S475);
    float _S479 = - _S476.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S476.differential_0.x + _S479, _S476.differential_0.y + _S479, _S476.differential_0.z);
    float3  _S480 = _S237 * s_diff_rgi_out_T_0;
    float3  _S481 = _S240 * s_diff_rgi_out_T_0;
    float _S482 = (_S480.x + _S480.y + _S480.z) / _S241;
    float _S483 = intensity_3 * - _S482;
    float _S484 = _S239 * _S482;
    DiffPair_float_0 _S485;
    (&_S485)->primal_0 = _S238;
    (&_S485)->differential_0 = 0.0f;
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = 0.0f;
    (&_S486)->differential_0 = 0.0f;
    _d_max_0(&_S485, &_S486, _S483);
    float3  _S487 = _S481 + make_float3 (0.0f, 0.0f, _S485.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S488;
    (&_S488)->primal_0 = H_6;
    (&_S488)->differential_0 = _S174;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S489;
    (&_S489)->primal_0 = rgi_in_0;
    (&_S489)->differential_0 = _S173;
    s_bwd_prop_mul_0(&_S488, &_S489, _S487);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S490 = _S488;
    float _S491 = _S484 + _S489.differential_0.z;
    float _S492 = _S489.differential_0.y + _S491;
    float _S493 = _S489.differential_0.x + _S491;
    float3  _S494 = make_float3 (_S493, _S492, _S491);
    if(_S230)
    {
        Matrix<float, 3, 3>  _S495 = _S229 * _S490.differential_0;
        Matrix<float, 3, 3>  _S496 = _S231 * _S490.differential_0;
        _S232 = - ((_S495.rows[int(0)].x + _S495.rows[int(0)].y + _S495.rows[int(0)].z + _S495.rows[int(1)].x + _S495.rows[int(1)].y + _S495.rows[int(1)].z + _S495.rows[int(2)].x + _S495.rows[int(2)].y + _S495.rows[int(2)].z) / _S232);
        H_6 = _S496;
    }
    else
    {
        _S232 = 0.0f;
        H_6 = _S490.differential_0;
    }
    DiffPair_float_0 _S497;
    (&_S497)->primal_0 = _S229.rows[int(2)].z;
    (&_S497)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S497, 0.0f);
    float _S498 = _S497.differential_0 + _S232;
    float3  _S499 = _S173;
    *&((&_S499)->z) = _S498;
    Matrix<float, 3, 3>  _S500 = _S174;
    _S500[int(2)] = _S499;
    Matrix<float, 3, 3>  _S501 = H_6 + _S500;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S502;
    (&_S502)->primal_0 = _S228;
    (&_S502)->differential_0 = _S174;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S503;
    (&_S503)->primal_0 = S_inv_0;
    (&_S503)->differential_0 = _S174;
    s_bwd_prop_mul_1(&_S502, &_S503, _S501);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S504;
    (&_S504)->primal_0 = T_3;
    (&_S504)->differential_0 = _S174;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S505;
    (&_S505)->primal_0 = D_0;
    (&_S505)->differential_0 = _S174;
    s_bwd_prop_mul_1(&_S504, &_S505, _S502.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S506 = _S504;
    float3  _S507 = make_float3 (_S505.differential_0.rows[int(0)].x, _S505.differential_0.rows[int(1)].y, _S505.differential_0.rows[int(2)].z);
    float3  _S508;
    if(_S223)
    {
        if(_S225)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S509;
            (&_S509)->primal_0 = r1_3;
            (&_S509)->differential_0 = _S173;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S510;
            (&_S510)->primal_0 = r2_15;
            (&_S510)->differential_0 = _S173;
            s_bwd_prop_cross_0(&_S509, &_S510, _S507);
            _S211 = _S173;
            lambda_v_9 = _S510.differential_0;
            _S508 = _S509.differential_0;
        }
        else
        {
            _S211 = _S507;
            lambda_v_9 = _S173;
            _S508 = _S173;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S511;
        (&_S511)->primal_0 = _S224;
        (&_S511)->differential_0 = _S173;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S512;
        (&_S512)->primal_0 = _S224;
        (&_S512)->differential_0 = _S173;
        s_bwd_prop_dot_0(&_S511, &_S512, 0.0f);
        float3  _S513 = _S512.differential_0 + _S511.differential_0 + _S211;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S514;
        (&_S514)->primal_0 = r0_3;
        (&_S514)->differential_0 = _S173;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S515;
        (&_S515)->primal_0 = r2_15;
        (&_S515)->differential_0 = _S173;
        s_bwd_prop_cross_0(&_S514, &_S515, _S513);
        float3  _S516 = _S515.differential_0 + lambda_v_9;
        _S211 = _S173;
        lambda_v_9 = _S516;
        _S224 = _S508;
        _S508 = _S514.differential_0;
    }
    else
    {
        _S211 = _S507;
        lambda_v_9 = _S173;
        _S224 = _S173;
        _S508 = _S173;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S517;
    (&_S517)->primal_0 = _S222;
    (&_S517)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S518;
    (&_S518)->primal_0 = _S222;
    (&_S518)->differential_0 = _S173;
    s_bwd_prop_dot_0(&_S517, &_S518, 0.0f);
    float3  _S519 = _S518.differential_0 + _S517.differential_0 + _S211;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S520;
    (&_S520)->primal_0 = r0_3;
    (&_S520)->differential_0 = _S173;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S521;
    (&_S521)->primal_0 = r1_3;
    (&_S521)->differential_0 = _S173;
    s_bwd_prop_cross_0(&_S520, &_S521, _S519);
    float3  _S522 = _S173;
    *&((&_S522)->z) = lambda_v_9.z;
    *&((&_S522)->y) = lambda_v_9.y;
    *&((&_S522)->x) = lambda_v_9.x;
    float3  _S523 = _S521.differential_0 + _S224;
    float3  _S524 = _S173;
    *&((&_S524)->z) = _S523.z;
    *&((&_S524)->y) = _S523.y;
    *&((&_S524)->x) = _S523.x;
    float3  _S525 = _S520.differential_0 + _S508;
    float3  _S526 = _S173;
    *&((&_S526)->z) = _S525.z;
    *&((&_S526)->y) = _S525.y;
    *&((&_S526)->x) = _S525.x;
    Matrix<float, 3, 3>  _S527 = _S174;
    _S527[int(2)] = _S522;
    _S527[int(1)] = _S524;
    _S527[int(0)] = _S526;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S528;
    (&_S528)->primal_0 = skew_0;
    (&_S528)->differential_0 = _S174;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S529;
    (&_S529)->primal_0 = T_3;
    (&_S529)->differential_0 = _S174;
    s_bwd_prop_mul_1(&_S528, &_S529, _S527);
    Matrix<float, 3, 3>  _S530 = _S529.differential_0 + _S506.differential_0;
    float2  _S531 = make_float2 (_S528.differential_0.rows[int(2)].y + - _S528.differential_0.rows[int(1)].z, _S528.differential_0.rows[int(0)].z + - _S528.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S532 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S533;
    (&_S533)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S533)->differential_0 = _S532;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S534;
    (&_S534)->primal_0 = _S214.color_params_2.n_0;
    (&_S534)->differential_0 = _S177;
    s_bwd_prop_mul_2(&_S533, &_S534, _S531);
    float2  _S535 = make_float2 (_S530.rows[int(0)].z, _S530.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S536;
    (&_S536)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S536)->differential_0 = _S532;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S537;
    (&_S537)->primal_0 = _S214.color_params_2.g_0;
    (&_S537)->differential_0 = _S177;
    s_bwd_prop_mul_2(&_S536, &_S537, _S535);
    float2  _S538 = make_float2 (_S530.rows[int(0)].y, _S530.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S539;
    (&_S539)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S539)->differential_0 = _S532;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S540;
    (&_S540)->primal_0 = _S214.color_params_2.r_0;
    (&_S540)->differential_0 = _S177;
    s_bwd_prop_mul_2(&_S539, &_S540, _S538);
    float2  _S541 = make_float2 (_S530.rows[int(0)].x, _S530.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S542;
    (&_S542)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S542)->differential_0 = _S532;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S543;
    (&_S543)->primal_0 = _S214.color_params_2.b_0;
    (&_S543)->differential_0 = _S177;
    s_bwd_prop_mul_2(&_S542, &_S543, _S541);
    ColorPPISPParams_0 _S544 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S544)->n_0 = _S534.differential_0;
    (&_S544)->g_0 = _S537.differential_0;
    (&_S544)->r_0 = _S540.differential_0;
    (&_S544)->b_0 = _S543.differential_0;
    _S211 = _S494;
    *&((&_S211)->z) = 0.0f;
    float _S545 = rgb_out_5.z * _S491;
    float _S546 = _S213 * _S491;
    DiffPair_float_0 _S547;
    (&_S547)->primal_0 = falloff_2;
    (&_S547)->differential_0 = 0.0f;
    DiffPair_float_0 _S548;
    (&_S548)->primal_0 = 0.0f;
    (&_S548)->differential_0 = 0.0f;
    DiffPair_float_0 _S549;
    (&_S549)->primal_0 = 1.0f;
    (&_S549)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S547, &_S548, &_S549, _S545);
    float _S550 = r2_14 * _S547.differential_0;
    float _S551 = r4_11 * _S547.differential_0;
    float s_diff_r6_T_0 = _S190 * _S547.differential_0;
    float _S552 = r6_2 * _S547.differential_0;
    float _S553 = r2_14 * (_S189 * _S547.differential_0 + r2_14 * s_diff_r6_T_0);
    float _S554 = _S188 * _S547.differential_0 + r4_11 * s_diff_r6_T_0 + _S553 + _S553;
    float _S555 = dy_11 * _S554;
    float _S556 = dx_11 * _S554;
    float _S557 = - (_S555 + _S555);
    float _S558 = - (_S556 + _S556);
    *&((&_S211)->y) = 0.0f;
    float _S559 = rgb_out_5.y * _S492;
    float _S560 = _S212 * _S492;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = falloff_1;
    (&_S561)->differential_0 = 0.0f;
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = 0.0f;
    (&_S562)->differential_0 = 0.0f;
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = 1.0f;
    (&_S563)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S561, &_S562, &_S563, _S559);
    float _S564 = r2_13 * _S561.differential_0;
    float _S565 = r4_10 * _S561.differential_0;
    float s_diff_r6_T_1 = _S187 * _S561.differential_0;
    float _S566 = r6_1 * _S561.differential_0;
    float _S567 = r2_13 * (_S186 * _S561.differential_0 + r2_13 * s_diff_r6_T_1);
    float _S568 = _S185 * _S561.differential_0 + r4_10 * s_diff_r6_T_1 + _S567 + _S567;
    float _S569 = dy_10 * _S568;
    float _S570 = dx_10 * _S568;
    float _S571 = - (_S569 + _S569);
    float _S572 = - (_S570 + _S570);
    *&((&_S211)->x) = 0.0f;
    float _S573 = rgb_out_5.x * _S493;
    float _S574 = _S209 * _S493;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = falloff_0;
    (&_S575)->differential_0 = 0.0f;
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = 0.0f;
    (&_S576)->differential_0 = 0.0f;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = 1.0f;
    (&_S577)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S575, &_S576, &_S577, _S573);
    float _S578 = r2_12 * _S575.differential_0;
    float _S579 = r4_9 * _S575.differential_0;
    float s_diff_r6_T_2 = _S184 * _S575.differential_0;
    float _S580 = r6_0 * _S575.differential_0;
    float _S581 = r2_12 * (_S183 * _S575.differential_0 + r2_12 * s_diff_r6_T_2);
    float _S582 = _S182 * _S575.differential_0 + r4_9 * s_diff_r6_T_2 + _S581 + _S581;
    float _S583 = dy_9 * _S582;
    float _S584 = dx_9 * _S582;
    float _S585 = - (_S583 + _S583);
    float _S586 = - (_S584 + _S584);
    float3  _S587 = _S173;
    *&((&_S587)->z) = _S546;
    *&((&_S587)->y) = _S560;
    *&((&_S587)->x) = _S574;
    float3  _S588 = _S211 + _S587;
    float3  _S589 = _S172.primal_0 * _S588;
    float3  _S590 = _S205 * _S588;
    float _S591 = _S589.x + _S589.y + _S589.z;
    DiffPair_float_0 _S592;
    (&_S592)->primal_0 = _S203.exposure_2;
    (&_S592)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S592, _S591);
    PPISPParams_0 _S593 = PPISPParams_x24_syn_dzero_0();
    (&_S593)->color_params_2 = _S544;
    (&_S593)->exposure_2 = _S592.differential_0;
    _S181 = _S593;
    (&(&_S181)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S594 = _S593.crf_params_1[int(2)].center_0 + _S375;
    (&(&_S181)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S595 = _S593.crf_params_1[int(2)].gamma_0 + _S378.differential_0;
    (&(&_S181)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S596 = _S593.crf_params_1[int(2)].shoulder_0 + _S382.differential_0;
    (&(&_S181)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S597 = _S593.crf_params_1[int(2)].toe_0 + _S386.differential_0;
    (&(&_S181)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S598 = _S593.crf_params_1[int(1)].center_0 + _S419;
    (&(&_S181)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S599 = _S593.crf_params_1[int(1)].gamma_0 + _S422.differential_0;
    (&(&_S181)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S600 = _S593.crf_params_1[int(1)].shoulder_0 + _S426.differential_0;
    (&(&_S181)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S601 = _S593.crf_params_1[int(1)].toe_0 + _S430.differential_0;
    (&(&_S181)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S602 = _S593.crf_params_1[int(0)].center_0 + _S463;
    (&(&_S181)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S603 = _S593.crf_params_1[int(0)].gamma_0 + _S466.differential_0;
    (&(&_S181)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S604 = _S593.crf_params_1[int(0)].shoulder_0 + _S470.differential_0;
    (&(&_S181)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S605 = _S593.crf_params_1[int(0)].toe_0 + _S474.differential_0;
    *&((&(&(&_S181)->color_params_2)->n_0)->y) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->n_0)->x) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->g_0)->y) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->g_0)->x) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->r_0)->y) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->r_0)->x) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->b_0)->y) = 0.0f;
    *&((&(&(&_S181)->color_params_2)->b_0)->x) = 0.0f;
    (&(&_S181)->vignette_params_2[int(2)])->alpha2_0 = 0.0f;
    float _S606 = _S552 + _S593.vignette_params_2[int(2)].alpha2_0;
    (&(&_S181)->vignette_params_2[int(2)])->alpha1_0 = 0.0f;
    float _S607 = _S551 + _S593.vignette_params_2[int(2)].alpha1_0;
    (&(&_S181)->vignette_params_2[int(2)])->alpha0_0 = 0.0f;
    float _S608 = _S550 + _S593.vignette_params_2[int(2)].alpha0_0;
    (&(&_S181)->vignette_params_2[int(2)])->cy_0 = 0.0f;
    float _S609 = _S557 + _S593.vignette_params_2[int(2)].cy_0;
    (&(&_S181)->vignette_params_2[int(2)])->cx_0 = 0.0f;
    float _S610 = _S558 + _S593.vignette_params_2[int(2)].cx_0;
    (&(&_S181)->vignette_params_2[int(1)])->alpha2_0 = 0.0f;
    float _S611 = _S566 + _S593.vignette_params_2[int(1)].alpha2_0;
    (&(&_S181)->vignette_params_2[int(1)])->alpha1_0 = 0.0f;
    float _S612 = _S565 + _S593.vignette_params_2[int(1)].alpha1_0;
    (&(&_S181)->vignette_params_2[int(1)])->alpha0_0 = 0.0f;
    float _S613 = _S564 + _S593.vignette_params_2[int(1)].alpha0_0;
    (&(&_S181)->vignette_params_2[int(1)])->cy_0 = 0.0f;
    float _S614 = _S571 + _S593.vignette_params_2[int(1)].cy_0;
    (&(&_S181)->vignette_params_2[int(1)])->cx_0 = 0.0f;
    float _S615 = _S572 + _S593.vignette_params_2[int(1)].cx_0;
    (&(&_S181)->vignette_params_2[int(0)])->alpha2_0 = 0.0f;
    float _S616 = _S580 + _S593.vignette_params_2[int(0)].alpha2_0;
    (&(&_S181)->vignette_params_2[int(0)])->alpha1_0 = 0.0f;
    float _S617 = _S579 + _S593.vignette_params_2[int(0)].alpha1_0;
    (&(&_S181)->vignette_params_2[int(0)])->alpha0_0 = 0.0f;
    float _S618 = _S578 + _S593.vignette_params_2[int(0)].alpha0_0;
    (&(&_S181)->vignette_params_2[int(0)])->cy_0 = 0.0f;
    float _S619 = _S585 + _S593.vignette_params_2[int(0)].cy_0;
    (&(&_S181)->vignette_params_2[int(0)])->cx_0 = 0.0f;
    float _S620 = _S586 + _S593.vignette_params_2[int(0)].cx_0;
    FixedArray<float, 36>  _S621;
    _S621[int(0)] = 0.0f;
    _S621[int(1)] = 0.0f;
    _S621[int(2)] = 0.0f;
    _S621[int(3)] = 0.0f;
    _S621[int(4)] = 0.0f;
    _S621[int(5)] = 0.0f;
    _S621[int(6)] = 0.0f;
    _S621[int(7)] = 0.0f;
    _S621[int(8)] = 0.0f;
    _S621[int(9)] = 0.0f;
    _S621[int(10)] = 0.0f;
    _S621[int(11)] = 0.0f;
    _S621[int(12)] = 0.0f;
    _S621[int(13)] = 0.0f;
    _S621[int(14)] = 0.0f;
    _S621[int(15)] = 0.0f;
    _S621[int(16)] = 0.0f;
    _S621[int(17)] = 0.0f;
    _S621[int(18)] = 0.0f;
    _S621[int(19)] = 0.0f;
    _S621[int(20)] = 0.0f;
    _S621[int(21)] = 0.0f;
    _S621[int(22)] = 0.0f;
    _S621[int(23)] = 0.0f;
    _S621[int(24)] = 0.0f;
    _S621[int(25)] = 0.0f;
    _S621[int(26)] = 0.0f;
    _S621[int(27)] = 0.0f;
    _S621[int(28)] = 0.0f;
    _S621[int(29)] = 0.0f;
    _S621[int(30)] = 0.0f;
    _S621[int(31)] = 0.0f;
    _S621[int(32)] = 0.0f;
    _S621[int(33)] = 0.0f;
    _S621[int(34)] = 0.0f;
    _S621[int(35)] = 0.0f;
    _S621[int(8)] = _S613;
    _S621[int(16)] = _S593.color_params_2.b_0.x;
    _S621[int(15)] = _S606;
    _S621[int(14)] = _S607;
    _S621[int(13)] = _S608;
    _S621[int(12)] = _S609;
    _S621[int(11)] = _S610;
    _S621[int(10)] = _S611;
    _S621[int(9)] = _S612;
    _S621[int(17)] = _S593.color_params_2.b_0.y;
    _S621[int(7)] = _S614;
    _S621[int(6)] = _S615;
    _S621[int(5)] = _S616;
    _S621[int(4)] = _S617;
    _S621[int(3)] = _S618;
    _S621[int(2)] = _S619;
    _S621[int(1)] = _S620;
    _S621[int(0)] = _S181.exposure_2;
    _S621[int(26)] = _S603;
    _S621[int(34)] = _S595;
    _S621[int(33)] = _S596;
    _S621[int(32)] = _S597;
    _S621[int(31)] = _S598;
    _S621[int(30)] = _S599;
    _S621[int(29)] = _S600;
    _S621[int(28)] = _S601;
    _S621[int(27)] = _S602;
    _S621[int(35)] = _S594;
    _S621[int(25)] = _S604;
    _S621[int(24)] = _S605;
    _S621[int(23)] = _S593.color_params_2.n_0.y;
    _S621[int(22)] = _S593.color_params_2.n_0.x;
    _S621[int(21)] = _S593.color_params_2.g_0.y;
    _S621[int(20)] = _S593.color_params_2.g_0.x;
    _S621[int(19)] = _S593.color_params_2.r_0.y;
    _S621[int(18)] = _S593.color_params_2.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S621;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S590;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S622, float2  _S623, float2  _S624, float2  _S625, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S626, float3  _S627)
{
    s_bwd_prop_apply_ppisp_0(_S622, _S623, _S624, _S625, _S626, _S627);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_3, float2  pix_coord_4, float2  image_center_4, float2  img_size_4, FixedArray<float, 36>  params_3, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S628 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_3;
    (&dp_rgb_in_0)->differential_0 = _S628;
    FixedArray<float, 36>  _S629 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_3;
    (&dp_params_0)->differential_0 = _S629;
    s_bwd_apply_ppisp_0(&dp_rgb_in_0, pix_coord_4, image_center_4, img_size_4, &dp_params_0, grad_out_0);
    *grad_rgb_in_0 = dp_rgb_in_0.differential_0;
    *grad_params_0 = (&dp_params_0)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C39x3E_0
{
    FixedArray<float, 39>  primal_0;
    FixedArray<float, 39>  differential_0;
};

inline __device__ void s_bwd_prop_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_1, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_1, float3  _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S630 = *dprgb_in_1;
    float3  _S631 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S632 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S633 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S634 = {
        _S633, _S633, _S633
    };
    float2  _S635 = make_float2 (0.0f);
    ColorPPISPParams_0 _S636 = { _S635, _S635, _S635, _S635 };
    RQSCRFPPISPChannelParams_0 _S637 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S638 = {
        _S637, _S637, _S637
    };
    PPISPParamsRQS_0 _S639;
    (&_S639)->exposure_1 = dpparams_1->primal_0[int(0)];
    (&_S639)->vignette_params_1 = _S634;
    (&_S639)->color_params_1 = _S636;
    (&_S639)->crf_params_0 = _S638;
    (&(&_S639)->vignette_params_1[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S639)->vignette_params_1[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S640 = dpparams_1->primal_0[int(3)];
    (&(&_S639)->vignette_params_1[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S641 = dpparams_1->primal_0[int(4)];
    (&(&_S639)->vignette_params_1[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S642 = dpparams_1->primal_0[int(5)];
    (&(&_S639)->vignette_params_1[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S639)->vignette_params_1[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S639)->vignette_params_1[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S643 = dpparams_1->primal_0[int(8)];
    (&(&_S639)->vignette_params_1[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S644 = dpparams_1->primal_0[int(9)];
    (&(&_S639)->vignette_params_1[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S645 = dpparams_1->primal_0[int(10)];
    (&(&_S639)->vignette_params_1[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S639)->vignette_params_1[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S639)->vignette_params_1[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S646 = dpparams_1->primal_0[int(13)];
    (&(&_S639)->vignette_params_1[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S647 = dpparams_1->primal_0[int(14)];
    (&(&_S639)->vignette_params_1[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S648 = dpparams_1->primal_0[int(15)];
    (&(&_S639)->vignette_params_1[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S639)->color_params_1)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S639)->color_params_1)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S639)->color_params_1)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S639)->color_params_1)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S639)->color_params_1)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S639)->color_params_1)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S639)->color_params_1)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S639)->color_params_1)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S649 = dpparams_1->primal_0[int(24)];
    (&(&_S639)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S650 = dpparams_1->primal_0[int(25)];
    (&(&_S639)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S651 = dpparams_1->primal_0[int(26)];
    (&(&_S639)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S652 = dpparams_1->primal_0[int(27)];
    (&(&_S639)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S653 = dpparams_1->primal_0[int(28)];
    (&(&_S639)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S654 = dpparams_1->primal_0[int(29)];
    (&(&_S639)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S655 = dpparams_1->primal_0[int(30)];
    (&(&_S639)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S656 = dpparams_1->primal_0[int(31)];
    (&(&_S639)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S657 = dpparams_1->primal_0[int(32)];
    (&(&_S639)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S658 = dpparams_1->primal_0[int(33)];
    (&(&_S639)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S659 = dpparams_1->primal_0[int(34)];
    (&(&_S639)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S660 = dpparams_1->primal_0[int(35)];
    (&(&_S639)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S661 = dpparams_1->primal_0[int(36)];
    (&(&_S639)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S662 = dpparams_1->primal_0[int(37)];
    (&(&_S639)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S663 = dpparams_1->primal_0[int(38)];
    (&(&_S639)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S664 = _S639;
    float _S665 = s_primal_ctx_exp2_0(_S639.exposure_1);
    float3  _S666 = make_float3 (_S665);
    float3  rgb_out_6 = (*dprgb_in_1).primal_0 * make_float3 (_S665);
    float _S667 = (F32_max((img_size_5.x), (img_size_5.y)));
    float _S668 = (pix_coord_5.x - image_center_5.x) / _S667;
    float _S669 = (pix_coord_5.y - image_center_5.y) / _S667;
    float dx_12 = _S668 - dpparams_1->primal_0[int(1)];
    float dy_12 = _S669 - dpparams_1->primal_0[int(2)];
    float r2_16 = dx_12 * dx_12 + dy_12 * dy_12;
    float r4_12 = r2_16 * r2_16;
    float r6_3 = r4_12 * r2_16;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_12 + dpparams_1->primal_0[int(3)] * r2_16 + 1.0f;
    float _S670 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S671 = rgb_out_6.x * _S670;
    float3  _S672 = rgb_out_6;
    *&((&_S672)->x) = _S671;
    float dx_13 = _S668 - dpparams_1->primal_0[int(6)];
    float dy_13 = _S669 - dpparams_1->primal_0[int(7)];
    float r2_17 = dx_13 * dx_13 + dy_13 * dy_13;
    float r4_13 = r2_17 * r2_17;
    float r6_4 = r4_13 * r2_17;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_13 + dpparams_1->primal_0[int(8)] * r2_17 + 1.0f;
    float _S673 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S672)->y) = rgb_out_6.y * _S673;
    float dx_14 = _S668 - dpparams_1->primal_0[int(11)];
    float dy_14 = _S669 - dpparams_1->primal_0[int(12)];
    float r2_18 = dx_14 * dx_14 + dy_14 * dy_14;
    float r4_14 = r2_18 * r2_18;
    float r6_5 = r4_14 * r2_18;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_14 + dpparams_1->primal_0[int(13)] * r2_18 + 1.0f;
    float _S674 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S672)->z) = rgb_out_6.z * _S674;
    PPISPParamsRQS_0 _S675 = _S639;
    float2  _S676 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S639.color_params_1.b_0);
    float2  _S677 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S639.color_params_1.r_0);
    float2  _S678 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S639.color_params_1.g_0);
    float2  _S679 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S639.color_params_1.n_0);
    float _S680 = 0.3333333432674408f + _S679.x;
    float _S681 = 0.3333333432674408f + _S679.y;
    Matrix<float, 3, 3>  T_4 = makeMatrix<float, 3, 3> (_S676.x, 1.0f + _S677.x, _S678.x, _S676.y, _S677.y, 1.0f + _S678.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S681, 1.0f, 0.0f, - _S680, - _S681, _S680, 0.0f);
    Matrix<float, 3, 3>  _S682 = s_primal_ctx_mul_1(skew_1, T_4);
    float3  r0_4 = make_float3 (_S682.rows[int(0)].x, _S682.rows[int(0)].y, _S682.rows[int(0)].z);
    float3  r1_4 = make_float3 (_S682.rows[int(1)].x, _S682.rows[int(1)].y, _S682.rows[int(1)].z);
    float3  r2_19 = make_float3 (_S682.rows[int(2)].x, _S682.rows[int(2)].y, _S682.rows[int(2)].z);
    float3  _S683 = s_primal_ctx_cross_0(r0_4, r1_4);
    bool _S684 = (s_primal_ctx_dot_0(_S683, _S683)) < 9.99999968265522539e-21f;
    float3  lambda_v_10;
    float3  _S685;
    bool _S686;
    if(_S684)
    {
        float3  _S687 = s_primal_ctx_cross_0(r0_4, r2_19);
        bool _S688 = (s_primal_ctx_dot_0(_S687, _S687)) < 9.99999968265522539e-21f;
        if(_S688)
        {
            lambda_v_10 = s_primal_ctx_cross_0(r1_4, r2_19);
        }
        else
        {
            lambda_v_10 = _S687;
        }
        _S686 = _S688;
        _S685 = _S687;
    }
    else
    {
        lambda_v_10 = _S683;
        _S686 = false;
        _S685 = _S631;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_10.x, 0.0f, 0.0f, 0.0f, lambda_v_10.y, 0.0f, 0.0f, 0.0f, lambda_v_10.z);
    Matrix<float, 3, 3>  _S689 = s_primal_ctx_mul_1(T_4, D_1);
    Matrix<float, 3, 3>  _S690 = s_primal_ctx_mul_1(_S689, S_inv_1);
    bool _S691 = (s_primal_ctx_abs_0(_S690.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_7;
    Matrix<float, 3, 3>  _S692;
    float _S693;
    if(_S691)
    {
        float inv_s_1 = 1.0f / _S690.rows[int(2)].z;
        Matrix<float, 3, 3>  _S694 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S695 = _S690.rows[int(2)].z * _S690.rows[int(2)].z;
        H_7 = _S690 * makeMatrix<float, 3, 3> (inv_s_1);
        _S692 = _S694;
        _S693 = _S695;
    }
    else
    {
        H_7 = _S690;
        _S692 = _S632;
        _S693 = 0.0f;
    }
    float _S696 = _S672.x;
    float _S697 = _S672.y;
    float intensity_4 = _S696 + _S697 + _S672.z;
    float3  rgi_in_1 = make_float3 (_S696, _S697, intensity_4);
    float3  _S698 = s_primal_ctx_mul_2(H_7, rgi_in_1);
    float _S699 = _S698.z;
    float _S700 = (F32_max((_S699), (0.0f))) + 0.00000999999974738f;
    float norm_factor_1 = intensity_4 / _S700;
    float3  _S701 = make_float3 (norm_factor_1);
    float _S702 = _S700 * _S700;
    float3  rgi_out_7 = _S698 * make_float3 (norm_factor_1);
    float _S703 = rgi_out_7.x;
    float _S704 = rgi_out_7.y;
    float3  _S705 = make_float3 (_S703, _S704, rgi_out_7.z - _S703 - _S704);
    float3  _S706 = make_float3 (0.0f);
    float3  _S707 = make_float3 (1.0f);
    float3  _S708 = s_primal_ctx_clamp_1(_S705, _S706, _S707);
    float _S709 = _S708.x;
    float _S710 = s_primal_ctx_exp_0(_S649);
    float _S711 = s_primal_ctx_exp_0(_S650);
    float _S712 = - _S651;
    float _S713 = 1.0f + s_primal_ctx_exp_0(_S712);
    float x0_4 = 1.0f / _S713;
    float _S714 = _S713 * _S713;
    float _S715 = - _S652;
    float _S716 = 1.0f + s_primal_ctx_exp_0(_S715);
    float y0_4 = 1.0f / _S716;
    float _S717 = _S716 * _S716;
    float _S718 = s_primal_ctx_exp_0(_S653);
    bool _S719 = _S709 < x0_4;
    float _S720;
    float _S721;
    float _S722;
    float _S723;
    float _S724;
    float _S725;
    float _S726;
    float _S727;
    float _S728;
    float _S729;
    float _S730;
    float _S731;
    float _S732;
    float _S733;
    float _S734;
    float _S735;
    float _S736;
    float _S737;
    float _S738;
    float _S739;
    float _S740;
    float _S741;
    float _S742;
    float _S743;
    float _S744;
    float _S745;
    float _S746;
    if(_S719)
    {
        float s0_3 = y0_4 / x0_4;
        float _S747 = x0_4 * x0_4;
        float t0_3 = _S709 / x0_4;
        float _S748 = s0_3 * t0_3;
        float _S749 = _S710 * t0_3;
        float _S750 = 1.0f - t0_3;
        float _S751 = _S748 * t0_3 + _S749 * _S750;
        float _S752 = y0_4 * _S751;
        float _S753 = _S710 + _S718 - 2.0f * s0_3;
        float _S754 = _S753 * t0_3;
        float _S755 = s0_3 + _S754 * _S750;
        _S720 = _S755 * _S755;
        _S721 = _S752;
        _S722 = _S755;
        _S723 = _S754;
        _S724 = _S750;
        _S725 = _S753;
        _S726 = t0_3;
        _S727 = _S751;
        _S728 = _S749;
        _S729 = _S748;
        _S730 = s0_3;
        _S731 = _S747;
        _S732 = 0.0f;
        _S733 = 0.0f;
        _S734 = 0.0f;
        _S735 = 0.0f;
        _S736 = 0.0f;
        _S737 = 0.0f;
        _S738 = 0.0f;
        _S739 = 0.0f;
        _S740 = 0.0f;
        _S741 = 0.0f;
        _S742 = 0.0f;
        _S743 = 0.0f;
        _S744 = 0.0f;
        _S745 = 0.0f;
        _S746 = 0.0f;
    }
    else
    {
        float _S756 = 1.0f - y0_4;
        float _S757 = 1.0f - x0_4;
        float s1_3 = _S756 / _S757;
        float _S758 = _S757 * _S757;
        float _S759 = _S709 - x0_4;
        float t1_3 = _S759 / _S757;
        float _S760 = s1_3 * t1_3;
        float _S761 = _S718 * t1_3;
        float _S762 = 1.0f - t1_3;
        float _S763 = _S760 * t1_3 + _S761 * _S762;
        float _S764 = _S756 * _S763;
        float _S765 = _S718 + _S711 - 2.0f * s1_3;
        float _S766 = _S765 * t1_3;
        float _S767 = s1_3 + _S766 * _S762;
        float _S768 = _S767 * _S767;
        _S720 = 0.0f;
        _S721 = 0.0f;
        _S722 = 0.0f;
        _S723 = 0.0f;
        _S724 = 0.0f;
        _S725 = 0.0f;
        _S726 = 0.0f;
        _S727 = 0.0f;
        _S728 = 0.0f;
        _S729 = 0.0f;
        _S730 = 0.0f;
        _S731 = 0.0f;
        _S732 = _S768;
        _S733 = _S764;
        _S734 = _S767;
        _S735 = _S766;
        _S736 = _S762;
        _S737 = _S765;
        _S738 = t1_3;
        _S739 = _S756;
        _S740 = _S763;
        _S741 = _S761;
        _S742 = _S760;
        _S743 = s1_3;
        _S744 = _S758;
        _S745 = _S759;
        _S746 = _S757;
    }
    float _S769 = _S708.y;
    float _S770 = s_primal_ctx_exp_0(_S654);
    float _S771 = s_primal_ctx_exp_0(_S655);
    float _S772 = - _S656;
    float _S773 = 1.0f + s_primal_ctx_exp_0(_S772);
    float x0_5 = 1.0f / _S773;
    float _S774 = _S773 * _S773;
    float _S775 = - _S657;
    float _S776 = 1.0f + s_primal_ctx_exp_0(_S775);
    float y0_5 = 1.0f / _S776;
    float _S777 = _S776 * _S776;
    float _S778 = s_primal_ctx_exp_0(_S658);
    bool _S779 = _S769 < x0_5;
    float _S780;
    float _S781;
    float _S782;
    float _S783;
    float _S784;
    float _S785;
    float _S786;
    float _S787;
    float _S788;
    float _S789;
    float _S790;
    float _S791;
    float _S792;
    float _S793;
    float _S794;
    float _S795;
    float _S796;
    float _S797;
    float _S798;
    float _S799;
    float _S800;
    float _S801;
    float _S802;
    float _S803;
    float _S804;
    float _S805;
    float _S806;
    if(_S779)
    {
        float s0_4 = y0_5 / x0_5;
        float _S807 = x0_5 * x0_5;
        float t0_4 = _S769 / x0_5;
        float _S808 = s0_4 * t0_4;
        float _S809 = _S770 * t0_4;
        float _S810 = 1.0f - t0_4;
        float _S811 = _S808 * t0_4 + _S809 * _S810;
        float _S812 = y0_5 * _S811;
        float _S813 = _S770 + _S778 - 2.0f * s0_4;
        float _S814 = _S813 * t0_4;
        float _S815 = s0_4 + _S814 * _S810;
        _S780 = _S815 * _S815;
        _S781 = _S812;
        _S782 = _S815;
        _S783 = _S814;
        _S784 = _S810;
        _S785 = _S813;
        _S786 = t0_4;
        _S787 = _S811;
        _S788 = _S809;
        _S789 = _S808;
        _S790 = s0_4;
        _S791 = _S807;
        _S792 = 0.0f;
        _S793 = 0.0f;
        _S794 = 0.0f;
        _S795 = 0.0f;
        _S796 = 0.0f;
        _S797 = 0.0f;
        _S798 = 0.0f;
        _S799 = 0.0f;
        _S800 = 0.0f;
        _S801 = 0.0f;
        _S802 = 0.0f;
        _S803 = 0.0f;
        _S804 = 0.0f;
        _S805 = 0.0f;
        _S806 = 0.0f;
    }
    else
    {
        float _S816 = 1.0f - y0_5;
        float _S817 = 1.0f - x0_5;
        float s1_4 = _S816 / _S817;
        float _S818 = _S817 * _S817;
        float _S819 = _S769 - x0_5;
        float t1_4 = _S819 / _S817;
        float _S820 = s1_4 * t1_4;
        float _S821 = _S778 * t1_4;
        float _S822 = 1.0f - t1_4;
        float _S823 = _S820 * t1_4 + _S821 * _S822;
        float _S824 = _S816 * _S823;
        float _S825 = _S778 + _S771 - 2.0f * s1_4;
        float _S826 = _S825 * t1_4;
        float _S827 = s1_4 + _S826 * _S822;
        float _S828 = _S827 * _S827;
        _S780 = 0.0f;
        _S781 = 0.0f;
        _S782 = 0.0f;
        _S783 = 0.0f;
        _S784 = 0.0f;
        _S785 = 0.0f;
        _S786 = 0.0f;
        _S787 = 0.0f;
        _S788 = 0.0f;
        _S789 = 0.0f;
        _S790 = 0.0f;
        _S791 = 0.0f;
        _S792 = _S828;
        _S793 = _S824;
        _S794 = _S827;
        _S795 = _S826;
        _S796 = _S822;
        _S797 = _S825;
        _S798 = t1_4;
        _S799 = _S816;
        _S800 = _S823;
        _S801 = _S821;
        _S802 = _S820;
        _S803 = s1_4;
        _S804 = _S818;
        _S805 = _S819;
        _S806 = _S817;
    }
    float _S829 = _S708.z;
    float _S830 = s_primal_ctx_exp_0(_S659);
    float _S831 = s_primal_ctx_exp_0(_S660);
    float _S832 = - _S661;
    float _S833 = 1.0f + s_primal_ctx_exp_0(_S832);
    float x0_6 = 1.0f / _S833;
    float _S834 = _S833 * _S833;
    float _S835 = - _S662;
    float _S836 = 1.0f + s_primal_ctx_exp_0(_S835);
    float y0_6 = 1.0f / _S836;
    float _S837 = _S836 * _S836;
    float _S838 = s_primal_ctx_exp_0(_S663);
    bool _S839 = _S829 < x0_6;
    float _S840;
    float _S841;
    float _S842;
    float _S843;
    float _S844;
    float _S845;
    float _S846;
    float _S847;
    float _S848;
    float _S849;
    float _S850;
    float _S851;
    float _S852;
    float _S853;
    float _S854;
    float _S855;
    float _S856;
    float _S857;
    float _S858;
    float _S859;
    float _S860;
    float _S861;
    float _S862;
    float _S863;
    float _S864;
    float _S865;
    float _S866;
    if(_S839)
    {
        float s0_5 = y0_6 / x0_6;
        float _S867 = x0_6 * x0_6;
        float t0_5 = _S829 / x0_6;
        float _S868 = s0_5 * t0_5;
        float _S869 = _S830 * t0_5;
        float _S870 = 1.0f - t0_5;
        float _S871 = _S868 * t0_5 + _S869 * _S870;
        float _S872 = y0_6 * _S871;
        float _S873 = _S830 + _S838 - 2.0f * s0_5;
        float _S874 = _S873 * t0_5;
        float _S875 = s0_5 + _S874 * _S870;
        _S840 = _S875 * _S875;
        _S841 = _S872;
        _S842 = _S875;
        _S843 = _S874;
        _S844 = _S870;
        _S845 = _S873;
        _S846 = t0_5;
        _S847 = _S871;
        _S848 = _S869;
        _S849 = _S868;
        _S850 = s0_5;
        _S851 = _S867;
        _S852 = 0.0f;
        _S853 = 0.0f;
        _S854 = 0.0f;
        _S855 = 0.0f;
        _S856 = 0.0f;
        _S857 = 0.0f;
        _S858 = 0.0f;
        _S859 = 0.0f;
        _S860 = 0.0f;
        _S861 = 0.0f;
        _S862 = 0.0f;
        _S863 = 0.0f;
        _S864 = 0.0f;
        _S865 = 0.0f;
        _S866 = 0.0f;
    }
    else
    {
        float _S876 = 1.0f - y0_6;
        float _S877 = 1.0f - x0_6;
        float s1_5 = _S876 / _S877;
        float _S878 = _S877 * _S877;
        float _S879 = _S829 - x0_6;
        float t1_5 = _S879 / _S877;
        float _S880 = s1_5 * t1_5;
        float _S881 = _S838 * t1_5;
        float _S882 = 1.0f - t1_5;
        float _S883 = _S880 * t1_5 + _S881 * _S882;
        float _S884 = _S876 * _S883;
        float _S885 = _S838 + _S831 - 2.0f * s1_5;
        float _S886 = _S885 * t1_5;
        float _S887 = s1_5 + _S886 * _S882;
        float _S888 = _S887 * _S887;
        _S840 = 0.0f;
        _S841 = 0.0f;
        _S842 = 0.0f;
        _S843 = 0.0f;
        _S844 = 0.0f;
        _S845 = 0.0f;
        _S846 = 0.0f;
        _S847 = 0.0f;
        _S848 = 0.0f;
        _S849 = 0.0f;
        _S850 = 0.0f;
        _S851 = 0.0f;
        _S852 = _S888;
        _S853 = _S884;
        _S854 = _S887;
        _S855 = _S886;
        _S856 = _S882;
        _S857 = _S885;
        _S858 = t1_5;
        _S859 = _S876;
        _S860 = _S883;
        _S861 = _S881;
        _S862 = _S880;
        _S863 = s1_5;
        _S864 = _S878;
        _S865 = _S879;
        _S866 = _S877;
    }
    if(_S839)
    {
        float _S889 = _s_dOut_1.z / _S840;
        float _S890 = _S841 * - _S889;
        float _S891 = _S842 * _S889;
        float _S892 = _S844 * _S890;
        float _S893 = _S846 * _S892;
        float _S894 = y0_6 * _S891;
        float _S895 = _S844 * _S894;
        float _S896 = _S846 * _S894;
        float _S897 = (_S845 * _S892 + - (_S843 * _S890 + _S848 * _S894) + _S830 * _S895 + _S849 * _S894 + _S850 * _S896) / _S851;
        float _S898 = x0_6 * _S897;
        float _S899 = (_S890 + 2.0f * - _S893 + _S846 * _S896) / _S851;
        float _S900 = _S847 * _S891 + x0_6 * _S899;
        float _S901 = _S893 + _S846 * _S895;
        float _S902 = _S829 * - _S897 + y0_6 * - _S899;
        _S840 = _S893;
        _S841 = _S900;
        _S842 = _S902;
        _S843 = 0.0f;
        _S844 = _S901;
        _S845 = _S898;
    }
    else
    {
        float _S903 = _s_dOut_1.z / _S852;
        float _S904 = _S853 * - _S903;
        float _S905 = _S854 * _S903;
        float _S906 = _S856 * _S904;
        float _S907 = _S858 * _S906;
        float _S908 = _S859 * _S905;
        float _S909 = _S856 * _S908;
        float _S910 = _S858 * _S908;
        float _S911 = (_S857 * _S906 + - (_S855 * _S904 + _S861 * _S908) + _S838 * _S909 + _S862 * _S908 + _S863 * _S910) / _S864;
        float _S912 = _S866 * _S911;
        float _S913 = (_S904 + 2.0f * - _S907 + _S858 * _S910) / _S864;
        float _S914 = _s_dOut_1.z + - (_S860 * _S905 + _S866 * _S913);
        float _S915 = - _S912 + - (_S865 * - _S911 + _S859 * - _S913);
        _S840 = _S907 + _S858 * _S909;
        _S841 = _S914;
        _S842 = _S915;
        _S843 = _S907;
        _S844 = 0.0f;
        _S845 = _S912;
    }
    DiffPair_float_0 _S916;
    (&_S916)->primal_0 = _S663;
    (&_S916)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S916, _S840);
    DiffPair_float_0 _S917 = _S916;
    float _S918 = - (_S841 / _S837);
    DiffPair_float_0 _S919;
    (&_S919)->primal_0 = _S835;
    (&_S919)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S919, _S918);
    float _S920 = - _S919.differential_0;
    float _S921 = - (_S842 / _S834);
    DiffPair_float_0 _S922;
    (&_S922)->primal_0 = _S832;
    (&_S922)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S922, _S921);
    float _S923 = - _S922.differential_0;
    DiffPair_float_0 _S924;
    (&_S924)->primal_0 = _S660;
    (&_S924)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S924, _S843);
    DiffPair_float_0 _S925 = _S924;
    DiffPair_float_0 _S926;
    (&_S926)->primal_0 = _S659;
    (&_S926)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S926, _S844);
    DiffPair_float_0 _S927 = _S926;
    float3  _S928 = make_float3 (0.0f, 0.0f, _S845);
    if(_S779)
    {
        float _S929 = _s_dOut_1.y / _S780;
        float _S930 = _S781 * - _S929;
        float _S931 = _S782 * _S929;
        float _S932 = _S784 * _S930;
        float _S933 = _S786 * _S932;
        float _S934 = y0_5 * _S931;
        float _S935 = _S784 * _S934;
        float _S936 = _S786 * _S934;
        float _S937 = (_S785 * _S932 + - (_S783 * _S930 + _S788 * _S934) + _S770 * _S935 + _S789 * _S934 + _S790 * _S936) / _S791;
        float _S938 = x0_5 * _S937;
        float _S939 = (_S930 + 2.0f * - _S933 + _S786 * _S936) / _S791;
        float _S940 = _S787 * _S931 + x0_5 * _S939;
        float _S941 = _S933 + _S786 * _S935;
        float _S942 = _S769 * - _S937 + y0_5 * - _S939;
        _S780 = _S933;
        _S781 = _S940;
        _S782 = _S942;
        _S783 = 0.0f;
        _S784 = _S941;
        _S785 = _S938;
    }
    else
    {
        float _S943 = _s_dOut_1.y / _S792;
        float _S944 = _S793 * - _S943;
        float _S945 = _S794 * _S943;
        float _S946 = _S796 * _S944;
        float _S947 = _S798 * _S946;
        float _S948 = _S799 * _S945;
        float _S949 = _S796 * _S948;
        float _S950 = _S798 * _S948;
        float _S951 = (_S797 * _S946 + - (_S795 * _S944 + _S801 * _S948) + _S778 * _S949 + _S802 * _S948 + _S803 * _S950) / _S804;
        float _S952 = _S806 * _S951;
        float _S953 = (_S944 + 2.0f * - _S947 + _S798 * _S950) / _S804;
        float _S954 = _s_dOut_1.y + - (_S800 * _S945 + _S806 * _S953);
        float _S955 = - _S952 + - (_S805 * - _S951 + _S799 * - _S953);
        _S780 = _S947 + _S798 * _S949;
        _S781 = _S954;
        _S782 = _S955;
        _S783 = _S947;
        _S784 = 0.0f;
        _S785 = _S952;
    }
    DiffPair_float_0 _S956;
    (&_S956)->primal_0 = _S658;
    (&_S956)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S956, _S780);
    DiffPair_float_0 _S957 = _S956;
    float _S958 = - (_S781 / _S777);
    DiffPair_float_0 _S959;
    (&_S959)->primal_0 = _S775;
    (&_S959)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S959, _S958);
    float _S960 = - _S959.differential_0;
    float _S961 = - (_S782 / _S774);
    DiffPair_float_0 _S962;
    (&_S962)->primal_0 = _S772;
    (&_S962)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S962, _S961);
    float _S963 = - _S962.differential_0;
    DiffPair_float_0 _S964;
    (&_S964)->primal_0 = _S655;
    (&_S964)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S964, _S783);
    DiffPair_float_0 _S965 = _S964;
    DiffPair_float_0 _S966;
    (&_S966)->primal_0 = _S654;
    (&_S966)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S966, _S784);
    DiffPair_float_0 _S967 = _S966;
    float3  _S968 = _S928 + make_float3 (0.0f, _S785, 0.0f);
    if(_S719)
    {
        float _S969 = _s_dOut_1.x / _S720;
        float _S970 = _S721 * - _S969;
        float _S971 = _S722 * _S969;
        float _S972 = _S724 * _S970;
        float _S973 = _S726 * _S972;
        float _S974 = y0_4 * _S971;
        float _S975 = _S724 * _S974;
        float _S976 = _S726 * _S974;
        float _S977 = (_S725 * _S972 + - (_S723 * _S970 + _S728 * _S974) + _S710 * _S975 + _S729 * _S974 + _S730 * _S976) / _S731;
        float _S978 = x0_4 * _S977;
        float _S979 = (_S970 + 2.0f * - _S973 + _S726 * _S976) / _S731;
        float _S980 = _S727 * _S971 + x0_4 * _S979;
        float _S981 = _S973 + _S726 * _S975;
        float _S982 = _S709 * - _S977 + y0_4 * - _S979;
        _S720 = _S973;
        _S721 = _S980;
        _S722 = _S982;
        _S723 = 0.0f;
        _S724 = _S981;
        _S725 = _S978;
    }
    else
    {
        float _S983 = _s_dOut_1.x / _S732;
        float _S984 = _S733 * - _S983;
        float _S985 = _S734 * _S983;
        float _S986 = _S736 * _S984;
        float _S987 = _S738 * _S986;
        float _S988 = _S739 * _S985;
        float _S989 = _S736 * _S988;
        float _S990 = _S738 * _S988;
        float _S991 = (_S737 * _S986 + - (_S735 * _S984 + _S741 * _S988) + _S718 * _S989 + _S742 * _S988 + _S743 * _S990) / _S744;
        float _S992 = _S746 * _S991;
        float _S993 = (_S984 + 2.0f * - _S987 + _S738 * _S990) / _S744;
        float _S994 = _s_dOut_1.x + - (_S740 * _S985 + _S746 * _S993);
        float _S995 = - _S992 + - (_S745 * - _S991 + _S739 * - _S993);
        _S720 = _S987 + _S738 * _S989;
        _S721 = _S994;
        _S722 = _S995;
        _S723 = _S987;
        _S724 = 0.0f;
        _S725 = _S992;
    }
    DiffPair_float_0 _S996;
    (&_S996)->primal_0 = _S653;
    (&_S996)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S996, _S720);
    DiffPair_float_0 _S997 = _S996;
    float _S998 = - (_S721 / _S717);
    DiffPair_float_0 _S999;
    (&_S999)->primal_0 = _S715;
    (&_S999)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S999, _S998);
    float _S1000 = - _S999.differential_0;
    float _S1001 = - (_S722 / _S714);
    DiffPair_float_0 _S1002;
    (&_S1002)->primal_0 = _S712;
    (&_S1002)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1002, _S1001);
    float _S1003 = - _S1002.differential_0;
    DiffPair_float_0 _S1004;
    (&_S1004)->primal_0 = _S650;
    (&_S1004)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1004, _S723);
    DiffPair_float_0 _S1005 = _S1004;
    DiffPair_float_0 _S1006;
    (&_S1006)->primal_0 = _S649;
    (&_S1006)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1006, _S724);
    DiffPair_float_0 _S1007 = _S1006;
    float3  _S1008 = _S968 + make_float3 (_S725, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1009;
    (&_S1009)->primal_0 = _S705;
    (&_S1009)->differential_0 = _S631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1010;
    (&_S1010)->primal_0 = _S706;
    (&_S1010)->differential_0 = _S631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1011;
    (&_S1011)->primal_0 = _S707;
    (&_S1011)->differential_0 = _S631;
    s_bwd_prop_clamp_0(&_S1009, &_S1010, &_S1011, _S1008);
    float _S1012 = - _S1009.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S1009.differential_0.x + _S1012, _S1009.differential_0.y + _S1012, _S1009.differential_0.z);
    float3  _S1013 = _S698 * s_diff_rgi_out_T_1;
    float3  _S1014 = _S701 * s_diff_rgi_out_T_1;
    float _S1015 = (_S1013.x + _S1013.y + _S1013.z) / _S702;
    float _S1016 = intensity_4 * - _S1015;
    float _S1017 = _S700 * _S1015;
    DiffPair_float_0 _S1018;
    (&_S1018)->primal_0 = _S699;
    (&_S1018)->differential_0 = 0.0f;
    DiffPair_float_0 _S1019;
    (&_S1019)->primal_0 = 0.0f;
    (&_S1019)->differential_0 = 0.0f;
    _d_max_0(&_S1018, &_S1019, _S1016);
    float3  _S1020 = _S1014 + make_float3 (0.0f, 0.0f, _S1018.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1021;
    (&_S1021)->primal_0 = H_7;
    (&_S1021)->differential_0 = _S632;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1022;
    (&_S1022)->primal_0 = rgi_in_1;
    (&_S1022)->differential_0 = _S631;
    s_bwd_prop_mul_0(&_S1021, &_S1022, _S1020);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1023 = _S1021;
    float _S1024 = _S1017 + _S1022.differential_0.z;
    float _S1025 = _S1022.differential_0.y + _S1024;
    float _S1026 = _S1022.differential_0.x + _S1024;
    float3  _S1027 = make_float3 (_S1026, _S1025, _S1024);
    if(_S691)
    {
        Matrix<float, 3, 3>  _S1028 = _S690 * _S1023.differential_0;
        Matrix<float, 3, 3>  _S1029 = _S692 * _S1023.differential_0;
        _S693 = - ((_S1028.rows[int(0)].x + _S1028.rows[int(0)].y + _S1028.rows[int(0)].z + _S1028.rows[int(1)].x + _S1028.rows[int(1)].y + _S1028.rows[int(1)].z + _S1028.rows[int(2)].x + _S1028.rows[int(2)].y + _S1028.rows[int(2)].z) / _S693);
        H_7 = _S1029;
    }
    else
    {
        _S693 = 0.0f;
        H_7 = _S1023.differential_0;
    }
    DiffPair_float_0 _S1030;
    (&_S1030)->primal_0 = _S690.rows[int(2)].z;
    (&_S1030)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1030, 0.0f);
    float _S1031 = _S1030.differential_0 + _S693;
    float3  _S1032 = _S631;
    *&((&_S1032)->z) = _S1031;
    Matrix<float, 3, 3>  _S1033 = _S632;
    _S1033[int(2)] = _S1032;
    Matrix<float, 3, 3>  _S1034 = H_7 + _S1033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1035;
    (&_S1035)->primal_0 = _S689;
    (&_S1035)->differential_0 = _S632;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1036;
    (&_S1036)->primal_0 = S_inv_1;
    (&_S1036)->differential_0 = _S632;
    s_bwd_prop_mul_1(&_S1035, &_S1036, _S1034);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1037;
    (&_S1037)->primal_0 = T_4;
    (&_S1037)->differential_0 = _S632;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1038;
    (&_S1038)->primal_0 = D_1;
    (&_S1038)->differential_0 = _S632;
    s_bwd_prop_mul_1(&_S1037, &_S1038, _S1035.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1039 = _S1037;
    float3  _S1040 = make_float3 (_S1038.differential_0.rows[int(0)].x, _S1038.differential_0.rows[int(1)].y, _S1038.differential_0.rows[int(2)].z);
    float3  _S1041;
    if(_S684)
    {
        if(_S686)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1042;
            (&_S1042)->primal_0 = r1_4;
            (&_S1042)->differential_0 = _S631;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1043;
            (&_S1043)->primal_0 = r2_19;
            (&_S1043)->differential_0 = _S631;
            s_bwd_prop_cross_0(&_S1042, &_S1043, _S1040);
            _S672 = _S631;
            lambda_v_10 = _S1043.differential_0;
            _S1041 = _S1042.differential_0;
        }
        else
        {
            _S672 = _S1040;
            lambda_v_10 = _S631;
            _S1041 = _S631;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1044;
        (&_S1044)->primal_0 = _S685;
        (&_S1044)->differential_0 = _S631;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1045;
        (&_S1045)->primal_0 = _S685;
        (&_S1045)->differential_0 = _S631;
        s_bwd_prop_dot_0(&_S1044, &_S1045, 0.0f);
        float3  _S1046 = _S1045.differential_0 + _S1044.differential_0 + _S672;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1047;
        (&_S1047)->primal_0 = r0_4;
        (&_S1047)->differential_0 = _S631;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1048;
        (&_S1048)->primal_0 = r2_19;
        (&_S1048)->differential_0 = _S631;
        s_bwd_prop_cross_0(&_S1047, &_S1048, _S1046);
        float3  _S1049 = _S1048.differential_0 + lambda_v_10;
        _S672 = _S631;
        lambda_v_10 = _S1049;
        _S685 = _S1041;
        _S1041 = _S1047.differential_0;
    }
    else
    {
        _S672 = _S1040;
        lambda_v_10 = _S631;
        _S685 = _S631;
        _S1041 = _S631;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1050;
    (&_S1050)->primal_0 = _S683;
    (&_S1050)->differential_0 = _S631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1051;
    (&_S1051)->primal_0 = _S683;
    (&_S1051)->differential_0 = _S631;
    s_bwd_prop_dot_0(&_S1050, &_S1051, 0.0f);
    float3  _S1052 = _S1051.differential_0 + _S1050.differential_0 + _S672;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1053;
    (&_S1053)->primal_0 = r0_4;
    (&_S1053)->differential_0 = _S631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1054;
    (&_S1054)->primal_0 = r1_4;
    (&_S1054)->differential_0 = _S631;
    s_bwd_prop_cross_0(&_S1053, &_S1054, _S1052);
    float3  _S1055 = _S631;
    *&((&_S1055)->z) = lambda_v_10.z;
    *&((&_S1055)->y) = lambda_v_10.y;
    *&((&_S1055)->x) = lambda_v_10.x;
    float3  _S1056 = _S1054.differential_0 + _S685;
    float3  _S1057 = _S631;
    *&((&_S1057)->z) = _S1056.z;
    *&((&_S1057)->y) = _S1056.y;
    *&((&_S1057)->x) = _S1056.x;
    float3  _S1058 = _S1053.differential_0 + _S1041;
    float3  _S1059 = _S631;
    *&((&_S1059)->z) = _S1058.z;
    *&((&_S1059)->y) = _S1058.y;
    *&((&_S1059)->x) = _S1058.x;
    Matrix<float, 3, 3>  _S1060 = _S632;
    _S1060[int(2)] = _S1055;
    _S1060[int(1)] = _S1057;
    _S1060[int(0)] = _S1059;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1061;
    (&_S1061)->primal_0 = skew_1;
    (&_S1061)->differential_0 = _S632;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1062;
    (&_S1062)->primal_0 = T_4;
    (&_S1062)->differential_0 = _S632;
    s_bwd_prop_mul_1(&_S1061, &_S1062, _S1060);
    Matrix<float, 3, 3>  _S1063 = _S1062.differential_0 + _S1039.differential_0;
    float2  _S1064 = make_float2 (_S1061.differential_0.rows[int(2)].y + - _S1061.differential_0.rows[int(1)].z, _S1061.differential_0.rows[int(0)].z + - _S1061.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1065 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1066;
    (&_S1066)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1066)->differential_0 = _S1065;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1067;
    (&_S1067)->primal_0 = _S675.color_params_1.n_0;
    (&_S1067)->differential_0 = _S635;
    s_bwd_prop_mul_2(&_S1066, &_S1067, _S1064);
    float2  _S1068 = make_float2 (_S1063.rows[int(0)].z, _S1063.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1069;
    (&_S1069)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1069)->differential_0 = _S1065;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1070;
    (&_S1070)->primal_0 = _S675.color_params_1.g_0;
    (&_S1070)->differential_0 = _S635;
    s_bwd_prop_mul_2(&_S1069, &_S1070, _S1068);
    float2  _S1071 = make_float2 (_S1063.rows[int(0)].y, _S1063.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1072;
    (&_S1072)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1072)->differential_0 = _S1065;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1073;
    (&_S1073)->primal_0 = _S675.color_params_1.r_0;
    (&_S1073)->differential_0 = _S635;
    s_bwd_prop_mul_2(&_S1072, &_S1073, _S1071);
    float2  _S1074 = make_float2 (_S1063.rows[int(0)].x, _S1063.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1075;
    (&_S1075)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1075)->differential_0 = _S1065;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1076;
    (&_S1076)->primal_0 = _S675.color_params_1.b_0;
    (&_S1076)->differential_0 = _S635;
    s_bwd_prop_mul_2(&_S1075, &_S1076, _S1074);
    ColorPPISPParams_0 _S1077 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1077)->n_0 = _S1067.differential_0;
    (&_S1077)->g_0 = _S1070.differential_0;
    (&_S1077)->r_0 = _S1073.differential_0;
    (&_S1077)->b_0 = _S1076.differential_0;
    _S672 = _S1027;
    *&((&_S672)->z) = 0.0f;
    float _S1078 = rgb_out_6.z * _S1024;
    float _S1079 = _S674 * _S1024;
    DiffPair_float_0 _S1080;
    (&_S1080)->primal_0 = falloff_5;
    (&_S1080)->differential_0 = 0.0f;
    DiffPair_float_0 _S1081;
    (&_S1081)->primal_0 = 0.0f;
    (&_S1081)->differential_0 = 0.0f;
    DiffPair_float_0 _S1082;
    (&_S1082)->primal_0 = 1.0f;
    (&_S1082)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1080, &_S1081, &_S1082, _S1078);
    float _S1083 = r2_18 * _S1080.differential_0;
    float _S1084 = r4_14 * _S1080.differential_0;
    float s_diff_r6_T_3 = _S648 * _S1080.differential_0;
    float _S1085 = r6_5 * _S1080.differential_0;
    float _S1086 = r2_18 * (_S647 * _S1080.differential_0 + r2_18 * s_diff_r6_T_3);
    float _S1087 = _S646 * _S1080.differential_0 + r4_14 * s_diff_r6_T_3 + _S1086 + _S1086;
    float _S1088 = dy_14 * _S1087;
    float _S1089 = dx_14 * _S1087;
    float _S1090 = - (_S1088 + _S1088);
    float _S1091 = - (_S1089 + _S1089);
    *&((&_S672)->y) = 0.0f;
    float _S1092 = rgb_out_6.y * _S1025;
    float _S1093 = _S673 * _S1025;
    DiffPair_float_0 _S1094;
    (&_S1094)->primal_0 = falloff_4;
    (&_S1094)->differential_0 = 0.0f;
    DiffPair_float_0 _S1095;
    (&_S1095)->primal_0 = 0.0f;
    (&_S1095)->differential_0 = 0.0f;
    DiffPair_float_0 _S1096;
    (&_S1096)->primal_0 = 1.0f;
    (&_S1096)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1094, &_S1095, &_S1096, _S1092);
    float _S1097 = r2_17 * _S1094.differential_0;
    float _S1098 = r4_13 * _S1094.differential_0;
    float s_diff_r6_T_4 = _S645 * _S1094.differential_0;
    float _S1099 = r6_4 * _S1094.differential_0;
    float _S1100 = r2_17 * (_S644 * _S1094.differential_0 + r2_17 * s_diff_r6_T_4);
    float _S1101 = _S643 * _S1094.differential_0 + r4_13 * s_diff_r6_T_4 + _S1100 + _S1100;
    float _S1102 = dy_13 * _S1101;
    float _S1103 = dx_13 * _S1101;
    float _S1104 = - (_S1102 + _S1102);
    float _S1105 = - (_S1103 + _S1103);
    *&((&_S672)->x) = 0.0f;
    float _S1106 = rgb_out_6.x * _S1026;
    float _S1107 = _S670 * _S1026;
    DiffPair_float_0 _S1108;
    (&_S1108)->primal_0 = falloff_3;
    (&_S1108)->differential_0 = 0.0f;
    DiffPair_float_0 _S1109;
    (&_S1109)->primal_0 = 0.0f;
    (&_S1109)->differential_0 = 0.0f;
    DiffPair_float_0 _S1110;
    (&_S1110)->primal_0 = 1.0f;
    (&_S1110)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1108, &_S1109, &_S1110, _S1106);
    float _S1111 = r2_16 * _S1108.differential_0;
    float _S1112 = r4_12 * _S1108.differential_0;
    float s_diff_r6_T_5 = _S642 * _S1108.differential_0;
    float _S1113 = r6_3 * _S1108.differential_0;
    float _S1114 = r2_16 * (_S641 * _S1108.differential_0 + r2_16 * s_diff_r6_T_5);
    float _S1115 = _S640 * _S1108.differential_0 + r4_12 * s_diff_r6_T_5 + _S1114 + _S1114;
    float _S1116 = dy_12 * _S1115;
    float _S1117 = dx_12 * _S1115;
    float _S1118 = - (_S1116 + _S1116);
    float _S1119 = - (_S1117 + _S1117);
    float3  _S1120 = _S631;
    *&((&_S1120)->z) = _S1079;
    *&((&_S1120)->y) = _S1093;
    *&((&_S1120)->x) = _S1107;
    float3  _S1121 = _S672 + _S1120;
    float3  _S1122 = _S630.primal_0 * _S1121;
    float3  _S1123 = _S666 * _S1121;
    float _S1124 = _S1122.x + _S1122.y + _S1122.z;
    DiffPair_float_0 _S1125;
    (&_S1125)->primal_0 = _S664.exposure_1;
    (&_S1125)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1125, _S1124);
    PPISPParamsRQS_0 _S1126 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S1126)->color_params_1 = _S1077;
    (&_S1126)->exposure_1 = _S1125.differential_0;
    _S639 = _S1126;
    (&(&_S639)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S1127 = _S1126.crf_params_0[int(2)].gc_0 + _S917.differential_0;
    (&(&_S639)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S1128 = _S1126.crf_params_0[int(2)].y0_0 + _S920;
    (&(&_S639)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S1129 = _S1126.crf_params_0[int(2)].x0_0 + _S923;
    (&(&_S639)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S1130 = _S1126.crf_params_0[int(2)].g1_0 + _S925.differential_0;
    (&(&_S639)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S1131 = _S1126.crf_params_0[int(2)].g0_0 + _S927.differential_0;
    (&(&_S639)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S1132 = _S1126.crf_params_0[int(1)].gc_0 + _S957.differential_0;
    (&(&_S639)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S1133 = _S1126.crf_params_0[int(1)].y0_0 + _S960;
    (&(&_S639)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S1134 = _S1126.crf_params_0[int(1)].x0_0 + _S963;
    (&(&_S639)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S1135 = _S1126.crf_params_0[int(1)].g1_0 + _S965.differential_0;
    (&(&_S639)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S1136 = _S1126.crf_params_0[int(1)].g0_0 + _S967.differential_0;
    (&(&_S639)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S1137 = _S1126.crf_params_0[int(0)].gc_0 + _S997.differential_0;
    (&(&_S639)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S1138 = _S1126.crf_params_0[int(0)].y0_0 + _S1000;
    (&(&_S639)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S1139 = _S1126.crf_params_0[int(0)].x0_0 + _S1003;
    (&(&_S639)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S1140 = _S1126.crf_params_0[int(0)].g1_0 + _S1005.differential_0;
    (&(&_S639)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S1141 = _S1126.crf_params_0[int(0)].g0_0 + _S1007.differential_0;
    *&((&(&(&_S639)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S639)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S639)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1142 = _S1085 + _S1126.vignette_params_1[int(2)].alpha2_0;
    (&(&_S639)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1143 = _S1084 + _S1126.vignette_params_1[int(2)].alpha1_0;
    (&(&_S639)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1144 = _S1083 + _S1126.vignette_params_1[int(2)].alpha0_0;
    (&(&_S639)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1145 = _S1090 + _S1126.vignette_params_1[int(2)].cy_0;
    (&(&_S639)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1146 = _S1091 + _S1126.vignette_params_1[int(2)].cx_0;
    (&(&_S639)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1147 = _S1099 + _S1126.vignette_params_1[int(1)].alpha2_0;
    (&(&_S639)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1148 = _S1098 + _S1126.vignette_params_1[int(1)].alpha1_0;
    (&(&_S639)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1149 = _S1097 + _S1126.vignette_params_1[int(1)].alpha0_0;
    (&(&_S639)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1150 = _S1104 + _S1126.vignette_params_1[int(1)].cy_0;
    (&(&_S639)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1151 = _S1105 + _S1126.vignette_params_1[int(1)].cx_0;
    (&(&_S639)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1152 = _S1113 + _S1126.vignette_params_1[int(0)].alpha2_0;
    (&(&_S639)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1153 = _S1112 + _S1126.vignette_params_1[int(0)].alpha1_0;
    (&(&_S639)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1154 = _S1111 + _S1126.vignette_params_1[int(0)].alpha0_0;
    (&(&_S639)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1155 = _S1118 + _S1126.vignette_params_1[int(0)].cy_0;
    (&(&_S639)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1156 = _S1119 + _S1126.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 39>  _S1157;
    _S1157[int(0)] = 0.0f;
    _S1157[int(1)] = 0.0f;
    _S1157[int(2)] = 0.0f;
    _S1157[int(3)] = 0.0f;
    _S1157[int(4)] = 0.0f;
    _S1157[int(5)] = 0.0f;
    _S1157[int(6)] = 0.0f;
    _S1157[int(7)] = 0.0f;
    _S1157[int(8)] = 0.0f;
    _S1157[int(9)] = 0.0f;
    _S1157[int(10)] = 0.0f;
    _S1157[int(11)] = 0.0f;
    _S1157[int(12)] = 0.0f;
    _S1157[int(13)] = 0.0f;
    _S1157[int(14)] = 0.0f;
    _S1157[int(15)] = 0.0f;
    _S1157[int(16)] = 0.0f;
    _S1157[int(17)] = 0.0f;
    _S1157[int(18)] = 0.0f;
    _S1157[int(19)] = 0.0f;
    _S1157[int(20)] = 0.0f;
    _S1157[int(21)] = 0.0f;
    _S1157[int(22)] = 0.0f;
    _S1157[int(23)] = 0.0f;
    _S1157[int(24)] = 0.0f;
    _S1157[int(25)] = 0.0f;
    _S1157[int(26)] = 0.0f;
    _S1157[int(27)] = 0.0f;
    _S1157[int(28)] = 0.0f;
    _S1157[int(29)] = 0.0f;
    _S1157[int(30)] = 0.0f;
    _S1157[int(31)] = 0.0f;
    _S1157[int(32)] = 0.0f;
    _S1157[int(33)] = 0.0f;
    _S1157[int(34)] = 0.0f;
    _S1157[int(35)] = 0.0f;
    _S1157[int(36)] = 0.0f;
    _S1157[int(37)] = 0.0f;
    _S1157[int(38)] = 0.0f;
    _S1157[int(9)] = _S1148;
    _S1157[int(18)] = _S1126.color_params_1.r_0.x;
    _S1157[int(17)] = _S1126.color_params_1.b_0.y;
    _S1157[int(16)] = _S1126.color_params_1.b_0.x;
    _S1157[int(15)] = _S1142;
    _S1157[int(14)] = _S1143;
    _S1157[int(13)] = _S1144;
    _S1157[int(12)] = _S1145;
    _S1157[int(11)] = _S1146;
    _S1157[int(10)] = _S1147;
    _S1157[int(19)] = _S1126.color_params_1.r_0.y;
    _S1157[int(8)] = _S1149;
    _S1157[int(7)] = _S1150;
    _S1157[int(6)] = _S1151;
    _S1157[int(5)] = _S1152;
    _S1157[int(4)] = _S1153;
    _S1157[int(3)] = _S1154;
    _S1157[int(2)] = _S1155;
    _S1157[int(1)] = _S1156;
    _S1157[int(0)] = _S639.exposure_1;
    _S1157[int(28)] = _S1137;
    _S1157[int(37)] = _S1128;
    _S1157[int(36)] = _S1129;
    _S1157[int(35)] = _S1130;
    _S1157[int(34)] = _S1131;
    _S1157[int(33)] = _S1132;
    _S1157[int(32)] = _S1133;
    _S1157[int(31)] = _S1134;
    _S1157[int(30)] = _S1135;
    _S1157[int(29)] = _S1136;
    _S1157[int(38)] = _S1127;
    _S1157[int(27)] = _S1138;
    _S1157[int(26)] = _S1139;
    _S1157[int(25)] = _S1140;
    _S1157[int(24)] = _S1141;
    _S1157[int(23)] = _S1126.color_params_1.n_0.y;
    _S1157[int(22)] = _S1126.color_params_1.n_0.x;
    _S1157[int(21)] = _S1126.color_params_1.g_0.y;
    _S1157[int(20)] = _S1126.color_params_1.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S1157;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S1123;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1158, float2  _S1159, float2  _S1160, float2  _S1161, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S1162, float3  _S1163)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S1158, _S1159, _S1160, _S1161, _S1162, _S1163);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_4, float2  pix_coord_6, float2  image_center_6, float2  img_size_6, FixedArray<float, 39>  params_4, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S1164 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_4;
    (&dp_rgb_in_1)->differential_0 = _S1164;
    FixedArray<float, 39>  _S1165 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_4;
    (&dp_params_1)->differential_0 = _S1165;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_6, image_center_6, img_size_6, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C24x3E_0
{
    FixedArray<float, 24>  primal_0;
    FixedArray<float, 24>  differential_0;
};

inline __device__ void s_bwd_prop_apply_ppisp_no_crf_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_2, float2  pix_coord_7, float2  image_center_7, float2  img_size_7, DiffPair_arrayx3Cfloatx2C24x3E_0 * dpparams_2, float3  _s_dOut_2)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1166 = *dprgb_in_2;
    float3  _S1167 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S1168 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S1169 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1170 = {
        _S1169, _S1169, _S1169
    };
    float2  _S1171 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1172 = { _S1171, _S1171, _S1171, _S1171 };
    PPISPParamsNoCRF_0 _S1173;
    (&_S1173)->exposure_0 = dpparams_2->primal_0[int(0)];
    (&_S1173)->vignette_params_0 = _S1170;
    (&_S1173)->color_params_0 = _S1172;
    (&(&_S1173)->vignette_params_0[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S1173)->vignette_params_0[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    float _S1174 = dpparams_2->primal_0[int(3)];
    (&(&_S1173)->vignette_params_0[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    float _S1175 = dpparams_2->primal_0[int(4)];
    (&(&_S1173)->vignette_params_0[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    float _S1176 = dpparams_2->primal_0[int(5)];
    (&(&_S1173)->vignette_params_0[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S1173)->vignette_params_0[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S1173)->vignette_params_0[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    float _S1177 = dpparams_2->primal_0[int(8)];
    (&(&_S1173)->vignette_params_0[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    float _S1178 = dpparams_2->primal_0[int(9)];
    (&(&_S1173)->vignette_params_0[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    float _S1179 = dpparams_2->primal_0[int(10)];
    (&(&_S1173)->vignette_params_0[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S1173)->vignette_params_0[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S1173)->vignette_params_0[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    float _S1180 = dpparams_2->primal_0[int(13)];
    (&(&_S1173)->vignette_params_0[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    float _S1181 = dpparams_2->primal_0[int(14)];
    (&(&_S1173)->vignette_params_0[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    float _S1182 = dpparams_2->primal_0[int(15)];
    (&(&_S1173)->vignette_params_0[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S1173)->color_params_0)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S1173)->color_params_0)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S1173)->color_params_0)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S1173)->color_params_0)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S1173)->color_params_0)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S1173)->color_params_0)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S1173)->color_params_0)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S1173)->color_params_0)->n_0)->y) = dpparams_2->primal_0[int(23)];
    PPISPParamsNoCRF_0 _S1183 = _S1173;
    float _S1184 = s_primal_ctx_exp2_0(_S1173.exposure_0);
    float3  _S1185 = make_float3 (_S1184);
    float3  rgb_out_7 = (*dprgb_in_2).primal_0 * make_float3 (_S1184);
    float _S1186 = (F32_max((img_size_7.x), (img_size_7.y)));
    float _S1187 = (pix_coord_7.x - image_center_7.x) / _S1186;
    float _S1188 = (pix_coord_7.y - image_center_7.y) / _S1186;
    float dx_15 = _S1187 - dpparams_2->primal_0[int(1)];
    float dy_15 = _S1188 - dpparams_2->primal_0[int(2)];
    float r2_20 = dx_15 * dx_15 + dy_15 * dy_15;
    float r4_15 = r2_20 * r2_20;
    float r6_6 = r4_15 * r2_20;
    float falloff_6 = dpparams_2->primal_0[int(5)] * r6_6 + dpparams_2->primal_0[int(4)] * r4_15 + dpparams_2->primal_0[int(3)] * r2_20 + 1.0f;
    float _S1189 = s_primal_ctx_clamp_0(falloff_6, 0.0f, 1.0f);
    float _S1190 = rgb_out_7.x * _S1189;
    float3  _S1191 = rgb_out_7;
    *&((&_S1191)->x) = _S1190;
    float dx_16 = _S1187 - dpparams_2->primal_0[int(6)];
    float dy_16 = _S1188 - dpparams_2->primal_0[int(7)];
    float r2_21 = dx_16 * dx_16 + dy_16 * dy_16;
    float r4_16 = r2_21 * r2_21;
    float r6_7 = r4_16 * r2_21;
    float falloff_7 = dpparams_2->primal_0[int(10)] * r6_7 + dpparams_2->primal_0[int(9)] * r4_16 + dpparams_2->primal_0[int(8)] * r2_21 + 1.0f;
    float _S1192 = s_primal_ctx_clamp_0(falloff_7, 0.0f, 1.0f);
    *&((&_S1191)->y) = rgb_out_7.y * _S1192;
    float dx_17 = _S1187 - dpparams_2->primal_0[int(11)];
    float dy_17 = _S1188 - dpparams_2->primal_0[int(12)];
    float r2_22 = dx_17 * dx_17 + dy_17 * dy_17;
    float r4_17 = r2_22 * r2_22;
    float r6_8 = r4_17 * r2_22;
    float falloff_8 = dpparams_2->primal_0[int(15)] * r6_8 + dpparams_2->primal_0[int(14)] * r4_17 + dpparams_2->primal_0[int(13)] * r2_22 + 1.0f;
    float _S1193 = s_primal_ctx_clamp_0(falloff_8, 0.0f, 1.0f);
    *&((&_S1191)->z) = rgb_out_7.z * _S1193;
    PPISPParamsNoCRF_0 _S1194 = _S1173;
    float2  _S1195 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S1173.color_params_0.b_0);
    float2  _S1196 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S1173.color_params_0.r_0);
    float2  _S1197 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S1173.color_params_0.g_0);
    float2  _S1198 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S1173.color_params_0.n_0);
    float _S1199 = 0.3333333432674408f + _S1198.x;
    float _S1200 = 0.3333333432674408f + _S1198.y;
    Matrix<float, 3, 3>  T_5 = makeMatrix<float, 3, 3> (_S1195.x, 1.0f + _S1196.x, _S1197.x, _S1195.y, _S1196.y, 1.0f + _S1197.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_2 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S1200, 1.0f, 0.0f, - _S1199, - _S1200, _S1199, 0.0f);
    Matrix<float, 3, 3>  _S1201 = s_primal_ctx_mul_1(skew_2, T_5);
    float3  r0_5 = make_float3 (_S1201.rows[int(0)].x, _S1201.rows[int(0)].y, _S1201.rows[int(0)].z);
    float3  r1_5 = make_float3 (_S1201.rows[int(1)].x, _S1201.rows[int(1)].y, _S1201.rows[int(1)].z);
    float3  r2_23 = make_float3 (_S1201.rows[int(2)].x, _S1201.rows[int(2)].y, _S1201.rows[int(2)].z);
    float3  _S1202 = s_primal_ctx_cross_0(r0_5, r1_5);
    bool _S1203 = (s_primal_ctx_dot_0(_S1202, _S1202)) < 9.99999968265522539e-21f;
    float3  lambda_v_11;
    float3  _S1204;
    bool _S1205;
    if(_S1203)
    {
        float3  _S1206 = s_primal_ctx_cross_0(r0_5, r2_23);
        bool _S1207 = (s_primal_ctx_dot_0(_S1206, _S1206)) < 9.99999968265522539e-21f;
        if(_S1207)
        {
            lambda_v_11 = s_primal_ctx_cross_0(r1_5, r2_23);
        }
        else
        {
            lambda_v_11 = _S1206;
        }
        _S1205 = _S1207;
        _S1204 = _S1206;
    }
    else
    {
        lambda_v_11 = _S1202;
        _S1205 = false;
        _S1204 = _S1167;
    }
    Matrix<float, 3, 3>  S_inv_2 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_2 = makeMatrix<float, 3, 3> (lambda_v_11.x, 0.0f, 0.0f, 0.0f, lambda_v_11.y, 0.0f, 0.0f, 0.0f, lambda_v_11.z);
    Matrix<float, 3, 3>  _S1208 = s_primal_ctx_mul_1(T_5, D_2);
    Matrix<float, 3, 3>  _S1209 = s_primal_ctx_mul_1(_S1208, S_inv_2);
    bool _S1210 = (s_primal_ctx_abs_0(_S1209.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_8;
    Matrix<float, 3, 3>  _S1211;
    float _S1212;
    if(_S1210)
    {
        float inv_s_2 = 1.0f / _S1209.rows[int(2)].z;
        Matrix<float, 3, 3>  _S1213 = makeMatrix<float, 3, 3> (inv_s_2);
        float _S1214 = _S1209.rows[int(2)].z * _S1209.rows[int(2)].z;
        H_8 = _S1209 * makeMatrix<float, 3, 3> (inv_s_2);
        _S1211 = _S1213;
        _S1212 = _S1214;
    }
    else
    {
        H_8 = _S1209;
        _S1211 = _S1168;
        _S1212 = 0.0f;
    }
    float _S1215 = _S1191.x;
    float _S1216 = _S1191.y;
    float intensity_5 = _S1215 + _S1216 + _S1191.z;
    float3  rgi_in_2 = make_float3 (_S1215, _S1216, intensity_5);
    float3  _S1217 = s_primal_ctx_mul_2(H_8, rgi_in_2);
    float _S1218 = _S1217.z;
    float _S1219 = (F32_max((_S1218), (0.0f))) + 0.00000999999974738f;
    float norm_factor_2 = intensity_5 / _S1219;
    float3  _S1220 = make_float3 (norm_factor_2);
    float _S1221 = _S1219 * _S1219;
    float3  rgi_out_8 = _S1217 * make_float3 (norm_factor_2);
    float _S1222 = rgi_out_8.x;
    float _S1223 = rgi_out_8.y;
    float3  _S1224 = make_float3 (0.0f);
    float3  _S1225 = make_float3 (1.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1226;
    (&_S1226)->primal_0 = make_float3 (_S1222, _S1223, rgi_out_8.z - _S1222 - _S1223);
    (&_S1226)->differential_0 = _S1167;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1227;
    (&_S1227)->primal_0 = _S1224;
    (&_S1227)->differential_0 = _S1167;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1228;
    (&_S1228)->primal_0 = _S1225;
    (&_S1228)->differential_0 = _S1167;
    s_bwd_prop_clamp_0(&_S1226, &_S1227, &_S1228, _s_dOut_2);
    float _S1229 = - _S1226.differential_0.z;
    float3  s_diff_rgi_out_T_2 = make_float3 (_S1226.differential_0.x + _S1229, _S1226.differential_0.y + _S1229, _S1226.differential_0.z);
    float3  _S1230 = _S1217 * s_diff_rgi_out_T_2;
    float3  _S1231 = _S1220 * s_diff_rgi_out_T_2;
    float _S1232 = (_S1230.x + _S1230.y + _S1230.z) / _S1221;
    float _S1233 = intensity_5 * - _S1232;
    float _S1234 = _S1219 * _S1232;
    DiffPair_float_0 _S1235;
    (&_S1235)->primal_0 = _S1218;
    (&_S1235)->differential_0 = 0.0f;
    DiffPair_float_0 _S1236;
    (&_S1236)->primal_0 = 0.0f;
    (&_S1236)->differential_0 = 0.0f;
    _d_max_0(&_S1235, &_S1236, _S1233);
    float3  _S1237 = _S1231 + make_float3 (0.0f, 0.0f, _S1235.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1238;
    (&_S1238)->primal_0 = H_8;
    (&_S1238)->differential_0 = _S1168;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239;
    (&_S1239)->primal_0 = rgi_in_2;
    (&_S1239)->differential_0 = _S1167;
    s_bwd_prop_mul_0(&_S1238, &_S1239, _S1237);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1240 = _S1238;
    float _S1241 = _S1234 + _S1239.differential_0.z;
    float _S1242 = _S1239.differential_0.y + _S1241;
    float _S1243 = _S1239.differential_0.x + _S1241;
    float3  _S1244 = make_float3 (_S1243, _S1242, _S1241);
    if(_S1210)
    {
        Matrix<float, 3, 3>  _S1245 = _S1209 * _S1240.differential_0;
        Matrix<float, 3, 3>  _S1246 = _S1211 * _S1240.differential_0;
        _S1212 = - ((_S1245.rows[int(0)].x + _S1245.rows[int(0)].y + _S1245.rows[int(0)].z + _S1245.rows[int(1)].x + _S1245.rows[int(1)].y + _S1245.rows[int(1)].z + _S1245.rows[int(2)].x + _S1245.rows[int(2)].y + _S1245.rows[int(2)].z) / _S1212);
        H_8 = _S1246;
    }
    else
    {
        _S1212 = 0.0f;
        H_8 = _S1240.differential_0;
    }
    DiffPair_float_0 _S1247;
    (&_S1247)->primal_0 = _S1209.rows[int(2)].z;
    (&_S1247)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1247, 0.0f);
    float _S1248 = _S1247.differential_0 + _S1212;
    float3  _S1249 = _S1167;
    *&((&_S1249)->z) = _S1248;
    Matrix<float, 3, 3>  _S1250 = _S1168;
    _S1250[int(2)] = _S1249;
    Matrix<float, 3, 3>  _S1251 = H_8 + _S1250;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1252;
    (&_S1252)->primal_0 = _S1208;
    (&_S1252)->differential_0 = _S1168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1253;
    (&_S1253)->primal_0 = S_inv_2;
    (&_S1253)->differential_0 = _S1168;
    s_bwd_prop_mul_1(&_S1252, &_S1253, _S1251);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1254;
    (&_S1254)->primal_0 = T_5;
    (&_S1254)->differential_0 = _S1168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1255;
    (&_S1255)->primal_0 = D_2;
    (&_S1255)->differential_0 = _S1168;
    s_bwd_prop_mul_1(&_S1254, &_S1255, _S1252.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1256 = _S1254;
    float3  _S1257 = make_float3 (_S1255.differential_0.rows[int(0)].x, _S1255.differential_0.rows[int(1)].y, _S1255.differential_0.rows[int(2)].z);
    float3  _S1258;
    if(_S1203)
    {
        if(_S1205)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1259;
            (&_S1259)->primal_0 = r1_5;
            (&_S1259)->differential_0 = _S1167;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1260;
            (&_S1260)->primal_0 = r2_23;
            (&_S1260)->differential_0 = _S1167;
            s_bwd_prop_cross_0(&_S1259, &_S1260, _S1257);
            _S1191 = _S1167;
            lambda_v_11 = _S1260.differential_0;
            _S1258 = _S1259.differential_0;
        }
        else
        {
            _S1191 = _S1257;
            lambda_v_11 = _S1167;
            _S1258 = _S1167;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1261;
        (&_S1261)->primal_0 = _S1204;
        (&_S1261)->differential_0 = _S1167;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1262;
        (&_S1262)->primal_0 = _S1204;
        (&_S1262)->differential_0 = _S1167;
        s_bwd_prop_dot_0(&_S1261, &_S1262, 0.0f);
        float3  _S1263 = _S1262.differential_0 + _S1261.differential_0 + _S1191;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1264;
        (&_S1264)->primal_0 = r0_5;
        (&_S1264)->differential_0 = _S1167;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1265;
        (&_S1265)->primal_0 = r2_23;
        (&_S1265)->differential_0 = _S1167;
        s_bwd_prop_cross_0(&_S1264, &_S1265, _S1263);
        float3  _S1266 = _S1265.differential_0 + lambda_v_11;
        _S1191 = _S1167;
        lambda_v_11 = _S1266;
        _S1204 = _S1258;
        _S1258 = _S1264.differential_0;
    }
    else
    {
        _S1191 = _S1257;
        lambda_v_11 = _S1167;
        _S1204 = _S1167;
        _S1258 = _S1167;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1267;
    (&_S1267)->primal_0 = _S1202;
    (&_S1267)->differential_0 = _S1167;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1268;
    (&_S1268)->primal_0 = _S1202;
    (&_S1268)->differential_0 = _S1167;
    s_bwd_prop_dot_0(&_S1267, &_S1268, 0.0f);
    float3  _S1269 = _S1268.differential_0 + _S1267.differential_0 + _S1191;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1270;
    (&_S1270)->primal_0 = r0_5;
    (&_S1270)->differential_0 = _S1167;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1271;
    (&_S1271)->primal_0 = r1_5;
    (&_S1271)->differential_0 = _S1167;
    s_bwd_prop_cross_0(&_S1270, &_S1271, _S1269);
    float3  _S1272 = _S1167;
    *&((&_S1272)->z) = lambda_v_11.z;
    *&((&_S1272)->y) = lambda_v_11.y;
    *&((&_S1272)->x) = lambda_v_11.x;
    float3  _S1273 = _S1271.differential_0 + _S1204;
    float3  _S1274 = _S1167;
    *&((&_S1274)->z) = _S1273.z;
    *&((&_S1274)->y) = _S1273.y;
    *&((&_S1274)->x) = _S1273.x;
    float3  _S1275 = _S1270.differential_0 + _S1258;
    float3  _S1276 = _S1167;
    *&((&_S1276)->z) = _S1275.z;
    *&((&_S1276)->y) = _S1275.y;
    *&((&_S1276)->x) = _S1275.x;
    Matrix<float, 3, 3>  _S1277 = _S1168;
    _S1277[int(2)] = _S1272;
    _S1277[int(1)] = _S1274;
    _S1277[int(0)] = _S1276;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1278;
    (&_S1278)->primal_0 = skew_2;
    (&_S1278)->differential_0 = _S1168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1279;
    (&_S1279)->primal_0 = T_5;
    (&_S1279)->differential_0 = _S1168;
    s_bwd_prop_mul_1(&_S1278, &_S1279, _S1277);
    Matrix<float, 3, 3>  _S1280 = _S1279.differential_0 + _S1256.differential_0;
    float2  _S1281 = make_float2 (_S1278.differential_0.rows[int(2)].y + - _S1278.differential_0.rows[int(1)].z, _S1278.differential_0.rows[int(0)].z + - _S1278.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1282 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1283;
    (&_S1283)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1283)->differential_0 = _S1282;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1284;
    (&_S1284)->primal_0 = _S1194.color_params_0.n_0;
    (&_S1284)->differential_0 = _S1171;
    s_bwd_prop_mul_2(&_S1283, &_S1284, _S1281);
    float2  _S1285 = make_float2 (_S1280.rows[int(0)].z, _S1280.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1286;
    (&_S1286)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1286)->differential_0 = _S1282;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1287;
    (&_S1287)->primal_0 = _S1194.color_params_0.g_0;
    (&_S1287)->differential_0 = _S1171;
    s_bwd_prop_mul_2(&_S1286, &_S1287, _S1285);
    float2  _S1288 = make_float2 (_S1280.rows[int(0)].y, _S1280.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1289;
    (&_S1289)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1289)->differential_0 = _S1282;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1290;
    (&_S1290)->primal_0 = _S1194.color_params_0.r_0;
    (&_S1290)->differential_0 = _S1171;
    s_bwd_prop_mul_2(&_S1289, &_S1290, _S1288);
    float2  _S1291 = make_float2 (_S1280.rows[int(0)].x, _S1280.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1292;
    (&_S1292)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1292)->differential_0 = _S1282;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1293;
    (&_S1293)->primal_0 = _S1194.color_params_0.b_0;
    (&_S1293)->differential_0 = _S1171;
    s_bwd_prop_mul_2(&_S1292, &_S1293, _S1291);
    ColorPPISPParams_0 _S1294 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1294)->n_0 = _S1284.differential_0;
    (&_S1294)->g_0 = _S1287.differential_0;
    (&_S1294)->r_0 = _S1290.differential_0;
    (&_S1294)->b_0 = _S1293.differential_0;
    _S1191 = _S1244;
    *&((&_S1191)->z) = 0.0f;
    float _S1295 = rgb_out_7.z * _S1241;
    float _S1296 = _S1193 * _S1241;
    DiffPair_float_0 _S1297;
    (&_S1297)->primal_0 = falloff_8;
    (&_S1297)->differential_0 = 0.0f;
    DiffPair_float_0 _S1298;
    (&_S1298)->primal_0 = 0.0f;
    (&_S1298)->differential_0 = 0.0f;
    DiffPair_float_0 _S1299;
    (&_S1299)->primal_0 = 1.0f;
    (&_S1299)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1297, &_S1298, &_S1299, _S1295);
    float _S1300 = r2_22 * _S1297.differential_0;
    float _S1301 = r4_17 * _S1297.differential_0;
    float s_diff_r6_T_6 = _S1182 * _S1297.differential_0;
    float _S1302 = r6_8 * _S1297.differential_0;
    float _S1303 = r2_22 * (_S1181 * _S1297.differential_0 + r2_22 * s_diff_r6_T_6);
    float _S1304 = _S1180 * _S1297.differential_0 + r4_17 * s_diff_r6_T_6 + _S1303 + _S1303;
    float _S1305 = dy_17 * _S1304;
    float _S1306 = dx_17 * _S1304;
    float _S1307 = - (_S1305 + _S1305);
    float _S1308 = - (_S1306 + _S1306);
    *&((&_S1191)->y) = 0.0f;
    float _S1309 = rgb_out_7.y * _S1242;
    float _S1310 = _S1192 * _S1242;
    DiffPair_float_0 _S1311;
    (&_S1311)->primal_0 = falloff_7;
    (&_S1311)->differential_0 = 0.0f;
    DiffPair_float_0 _S1312;
    (&_S1312)->primal_0 = 0.0f;
    (&_S1312)->differential_0 = 0.0f;
    DiffPair_float_0 _S1313;
    (&_S1313)->primal_0 = 1.0f;
    (&_S1313)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1311, &_S1312, &_S1313, _S1309);
    float _S1314 = r2_21 * _S1311.differential_0;
    float _S1315 = r4_16 * _S1311.differential_0;
    float s_diff_r6_T_7 = _S1179 * _S1311.differential_0;
    float _S1316 = r6_7 * _S1311.differential_0;
    float _S1317 = r2_21 * (_S1178 * _S1311.differential_0 + r2_21 * s_diff_r6_T_7);
    float _S1318 = _S1177 * _S1311.differential_0 + r4_16 * s_diff_r6_T_7 + _S1317 + _S1317;
    float _S1319 = dy_16 * _S1318;
    float _S1320 = dx_16 * _S1318;
    float _S1321 = - (_S1319 + _S1319);
    float _S1322 = - (_S1320 + _S1320);
    *&((&_S1191)->x) = 0.0f;
    float _S1323 = rgb_out_7.x * _S1243;
    float _S1324 = _S1189 * _S1243;
    DiffPair_float_0 _S1325;
    (&_S1325)->primal_0 = falloff_6;
    (&_S1325)->differential_0 = 0.0f;
    DiffPair_float_0 _S1326;
    (&_S1326)->primal_0 = 0.0f;
    (&_S1326)->differential_0 = 0.0f;
    DiffPair_float_0 _S1327;
    (&_S1327)->primal_0 = 1.0f;
    (&_S1327)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1325, &_S1326, &_S1327, _S1323);
    float _S1328 = r2_20 * _S1325.differential_0;
    float _S1329 = r4_15 * _S1325.differential_0;
    float s_diff_r6_T_8 = _S1176 * _S1325.differential_0;
    float _S1330 = r6_6 * _S1325.differential_0;
    float _S1331 = r2_20 * (_S1175 * _S1325.differential_0 + r2_20 * s_diff_r6_T_8);
    float _S1332 = _S1174 * _S1325.differential_0 + r4_15 * s_diff_r6_T_8 + _S1331 + _S1331;
    float _S1333 = dy_15 * _S1332;
    float _S1334 = dx_15 * _S1332;
    float _S1335 = - (_S1333 + _S1333);
    float _S1336 = - (_S1334 + _S1334);
    float3  _S1337 = _S1167;
    *&((&_S1337)->z) = _S1296;
    *&((&_S1337)->y) = _S1310;
    *&((&_S1337)->x) = _S1324;
    float3  _S1338 = _S1191 + _S1337;
    float3  _S1339 = _S1166.primal_0 * _S1338;
    float3  _S1340 = _S1185 * _S1338;
    float _S1341 = _S1339.x + _S1339.y + _S1339.z;
    DiffPair_float_0 _S1342;
    (&_S1342)->primal_0 = _S1183.exposure_0;
    (&_S1342)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1342, _S1341);
    PPISPParamsNoCRF_0 _S1343 = PPISPParamsNoCRF_x24_syn_dzero_0();
    (&_S1343)->color_params_0 = _S1294;
    (&_S1343)->exposure_0 = _S1342.differential_0;
    _S1173 = _S1343;
    *&((&(&(&_S1173)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1173)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1173)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1344 = _S1302 + _S1343.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1173)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1345 = _S1301 + _S1343.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1173)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1346 = _S1300 + _S1343.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1173)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1347 = _S1307 + _S1343.vignette_params_0[int(2)].cy_0;
    (&(&_S1173)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1348 = _S1308 + _S1343.vignette_params_0[int(2)].cx_0;
    (&(&_S1173)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1349 = _S1316 + _S1343.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1173)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1350 = _S1315 + _S1343.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1173)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1351 = _S1314 + _S1343.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1173)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1352 = _S1321 + _S1343.vignette_params_0[int(1)].cy_0;
    (&(&_S1173)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1353 = _S1322 + _S1343.vignette_params_0[int(1)].cx_0;
    (&(&_S1173)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1354 = _S1330 + _S1343.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1173)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1355 = _S1329 + _S1343.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1173)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1356 = _S1328 + _S1343.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1173)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1357 = _S1335 + _S1343.vignette_params_0[int(0)].cy_0;
    (&(&_S1173)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1358 = _S1336 + _S1343.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 24>  _S1359;
    _S1359[int(0)] = 0.0f;
    _S1359[int(1)] = 0.0f;
    _S1359[int(2)] = 0.0f;
    _S1359[int(3)] = 0.0f;
    _S1359[int(4)] = 0.0f;
    _S1359[int(5)] = 0.0f;
    _S1359[int(6)] = 0.0f;
    _S1359[int(7)] = 0.0f;
    _S1359[int(8)] = 0.0f;
    _S1359[int(9)] = 0.0f;
    _S1359[int(10)] = 0.0f;
    _S1359[int(11)] = 0.0f;
    _S1359[int(12)] = 0.0f;
    _S1359[int(13)] = 0.0f;
    _S1359[int(14)] = 0.0f;
    _S1359[int(15)] = 0.0f;
    _S1359[int(16)] = 0.0f;
    _S1359[int(17)] = 0.0f;
    _S1359[int(18)] = 0.0f;
    _S1359[int(19)] = 0.0f;
    _S1359[int(20)] = 0.0f;
    _S1359[int(21)] = 0.0f;
    _S1359[int(22)] = 0.0f;
    _S1359[int(23)] = 0.0f;
    _S1359[int(11)] = _S1348;
    _S1359[int(0)] = _S1173.exposure_0;
    _S1359[int(1)] = _S1358;
    _S1359[int(2)] = _S1357;
    _S1359[int(3)] = _S1356;
    _S1359[int(4)] = _S1355;
    _S1359[int(5)] = _S1354;
    _S1359[int(6)] = _S1353;
    _S1359[int(7)] = _S1352;
    _S1359[int(8)] = _S1351;
    _S1359[int(9)] = _S1350;
    _S1359[int(10)] = _S1349;
    _S1359[int(23)] = _S1343.color_params_0.n_0.y;
    _S1359[int(12)] = _S1347;
    _S1359[int(13)] = _S1346;
    _S1359[int(14)] = _S1345;
    _S1359[int(15)] = _S1344;
    _S1359[int(16)] = _S1343.color_params_0.b_0.x;
    _S1359[int(17)] = _S1343.color_params_0.b_0.y;
    _S1359[int(18)] = _S1343.color_params_0.r_0.x;
    _S1359[int(19)] = _S1343.color_params_0.r_0.y;
    _S1359[int(20)] = _S1343.color_params_0.g_0.x;
    _S1359[int(21)] = _S1343.color_params_0.g_0.y;
    _S1359[int(22)] = _S1343.color_params_0.n_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S1359;
    dprgb_in_2->primal_0 = (*dprgb_in_2).primal_0;
    dprgb_in_2->differential_0 = _S1340;
    return;
}

inline __device__ void s_bwd_apply_ppisp_no_crf_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1360, float2  _S1361, float2  _S1362, float2  _S1363, DiffPair_arrayx3Cfloatx2C24x3E_0 * _S1364, float3  _S1365)
{
    s_bwd_prop_apply_ppisp_no_crf_0(_S1360, _S1361, _S1362, _S1363, _S1364, _S1365);
    return;
}

inline __device__ void apply_ppisp_no_crf_vjp(float3  rgb_in_5, float2  pix_coord_8, float2  image_center_8, float2  img_size_8, FixedArray<float, 24>  params_5, float3  grad_out_2, float3  * grad_rgb_in_2, FixedArray<float, 24>  * grad_params_2)
{
    float3  _S1366 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_2;
    (&dp_rgb_in_2)->primal_0 = rgb_in_5;
    (&dp_rgb_in_2)->differential_0 = _S1366;
    FixedArray<float, 24>  _S1367 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C24x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_5;
    (&dp_params_2)->differential_0 = _S1367;
    s_bwd_apply_ppisp_no_crf_0(&dp_rgb_in_2, pix_coord_8, image_center_8, img_size_8, &dp_params_2, grad_out_2);
    *grad_rgb_in_2 = dp_rgb_in_2.differential_0;
    *grad_params_2 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_6, FixedArray<float, 22>  * _S1368)
{
    PPISPParams_0 p_3;
    (&p_3)->exposure_2 = params_6[int(0)];
    (&(&p_3)->vignette_params_2[int(0)])->cx_0 = params_6[int(1)];
    (&(&p_3)->vignette_params_2[int(0)])->cy_0 = params_6[int(2)];
    (&(&p_3)->vignette_params_2[int(0)])->alpha0_0 = params_6[int(3)];
    (&(&p_3)->vignette_params_2[int(0)])->alpha1_0 = params_6[int(4)];
    (&(&p_3)->vignette_params_2[int(0)])->alpha2_0 = params_6[int(5)];
    (&(&p_3)->vignette_params_2[int(1)])->cx_0 = params_6[int(6)];
    (&(&p_3)->vignette_params_2[int(1)])->cy_0 = params_6[int(7)];
    (&(&p_3)->vignette_params_2[int(1)])->alpha0_0 = params_6[int(8)];
    (&(&p_3)->vignette_params_2[int(1)])->alpha1_0 = params_6[int(9)];
    (&(&p_3)->vignette_params_2[int(1)])->alpha2_0 = params_6[int(10)];
    (&(&p_3)->vignette_params_2[int(2)])->cx_0 = params_6[int(11)];
    (&(&p_3)->vignette_params_2[int(2)])->cy_0 = params_6[int(12)];
    (&(&p_3)->vignette_params_2[int(2)])->alpha0_0 = params_6[int(13)];
    (&(&p_3)->vignette_params_2[int(2)])->alpha1_0 = params_6[int(14)];
    (&(&p_3)->vignette_params_2[int(2)])->alpha2_0 = params_6[int(15)];
    *&((&(&(&p_3)->color_params_2)->b_0)->x) = params_6[int(16)];
    *&((&(&(&p_3)->color_params_2)->b_0)->y) = params_6[int(17)];
    *&((&(&(&p_3)->color_params_2)->r_0)->x) = params_6[int(18)];
    *&((&(&(&p_3)->color_params_2)->r_0)->y) = params_6[int(19)];
    *&((&(&(&p_3)->color_params_2)->g_0)->x) = params_6[int(20)];
    *&((&(&(&p_3)->color_params_2)->g_0)->y) = params_6[int(21)];
    *&((&(&(&p_3)->color_params_2)->n_0)->x) = params_6[int(22)];
    *&((&(&(&p_3)->color_params_2)->n_0)->y) = params_6[int(23)];
    (&(&p_3)->crf_params_1[int(0)])->toe_0 = params_6[int(24)];
    (&(&p_3)->crf_params_1[int(0)])->shoulder_0 = params_6[int(25)];
    (&(&p_3)->crf_params_1[int(0)])->gamma_0 = params_6[int(26)];
    (&(&p_3)->crf_params_1[int(0)])->center_0 = params_6[int(27)];
    (&(&p_3)->crf_params_1[int(1)])->toe_0 = params_6[int(28)];
    (&(&p_3)->crf_params_1[int(1)])->shoulder_0 = params_6[int(29)];
    (&(&p_3)->crf_params_1[int(1)])->gamma_0 = params_6[int(30)];
    (&(&p_3)->crf_params_1[int(1)])->center_0 = params_6[int(31)];
    (&(&p_3)->crf_params_1[int(2)])->toe_0 = params_6[int(32)];
    (&(&p_3)->crf_params_1[int(2)])->shoulder_0 = params_6[int(33)];
    (&(&p_3)->crf_params_1[int(2)])->gamma_0 = params_6[int(34)];
    (&(&p_3)->crf_params_1[int(2)])->center_0 = params_6[int(35)];
    FixedArray<float, 22>  losses_0;
    losses_0[int(0)] = 0.0f;
    losses_0[int(1)] = 0.0f;
    losses_0[int(2)] = 0.0f;
    losses_0[int(3)] = 0.0f;
    losses_0[int(4)] = 0.0f;
    losses_0[int(5)] = 0.0f;
    losses_0[int(6)] = 0.0f;
    losses_0[int(7)] = 0.0f;
    losses_0[int(8)] = 0.0f;
    losses_0[int(9)] = 0.0f;
    losses_0[int(10)] = 0.0f;
    losses_0[int(11)] = 0.0f;
    losses_0[int(12)] = 0.0f;
    losses_0[int(13)] = 0.0f;
    losses_0[int(14)] = 0.0f;
    losses_0[int(15)] = 0.0f;
    losses_0[int(16)] = 0.0f;
    losses_0[int(17)] = 0.0f;
    losses_0[int(18)] = 0.0f;
    losses_0[int(19)] = 0.0f;
    losses_0[int(20)] = 0.0f;
    losses_0[int(21)] = 0.0f;
    losses_0[int(0)] = p_3.exposure_2;
    float _S1369 = p_3.vignette_params_2[int(0)].cx_0;
    float _S1370 = p_3.vignette_params_2[int(0)].cy_0;
    float _S1371 = p_3.vignette_params_2[int(1)].cx_0;
    float _S1372 = p_3.vignette_params_2[int(1)].cy_0;
    float _S1373 = p_3.vignette_params_2[int(2)].cx_0;
    float _S1374 = p_3.vignette_params_2[int(2)].cy_0;
    losses_0[int(1)] = _S1369 * _S1369 + _S1370 * _S1370 + _S1371 * _S1371 + _S1372 * _S1372 + _S1373 * _S1373 + _S1374 * _S1374;
    losses_0[int(2)] = (F32_max((0.0f), (p_3.vignette_params_2[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(2)].alpha0_0)));
    losses_0[int(3)] = (F32_max((0.0f), (p_3.vignette_params_2[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(2)].alpha1_0)));
    losses_0[int(4)] = (F32_max((0.0f), (p_3.vignette_params_2[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_2[int(2)].alpha2_0)));
    float mean_0 = (p_3.vignette_params_2[int(0)].cx_0 + p_3.vignette_params_2[int(1)].cx_0 + p_3.vignette_params_2[int(2)].cx_0) / 3.0f;
    float _S1375 = p_3.vignette_params_2[int(0)].cx_0 - mean_0;
    float _S1376 = p_3.vignette_params_2[int(1)].cx_0 - mean_0;
    float _S1377 = p_3.vignette_params_2[int(2)].cx_0 - mean_0;
    losses_0[int(5)] = (_S1375 * _S1375 + _S1376 * _S1376 + _S1377 * _S1377) / 3.0f;
    float mean_1 = (p_3.vignette_params_2[int(0)].cy_0 + p_3.vignette_params_2[int(1)].cy_0 + p_3.vignette_params_2[int(2)].cy_0) / 3.0f;
    float _S1378 = p_3.vignette_params_2[int(0)].cy_0 - mean_1;
    float _S1379 = p_3.vignette_params_2[int(1)].cy_0 - mean_1;
    float _S1380 = p_3.vignette_params_2[int(2)].cy_0 - mean_1;
    losses_0[int(6)] = (_S1378 * _S1378 + _S1379 * _S1379 + _S1380 * _S1380) / 3.0f;
    float mean_2 = (p_3.vignette_params_2[int(0)].alpha0_0 + p_3.vignette_params_2[int(1)].alpha0_0 + p_3.vignette_params_2[int(2)].alpha0_0) / 3.0f;
    float _S1381 = p_3.vignette_params_2[int(0)].alpha0_0 - mean_2;
    float _S1382 = p_3.vignette_params_2[int(1)].alpha0_0 - mean_2;
    float _S1383 = p_3.vignette_params_2[int(2)].alpha0_0 - mean_2;
    losses_0[int(7)] = (_S1381 * _S1381 + _S1382 * _S1382 + _S1383 * _S1383) / 3.0f;
    float mean_3 = (p_3.vignette_params_2[int(0)].alpha1_0 + p_3.vignette_params_2[int(1)].alpha1_0 + p_3.vignette_params_2[int(2)].alpha1_0) / 3.0f;
    float _S1384 = p_3.vignette_params_2[int(0)].alpha1_0 - mean_3;
    float _S1385 = p_3.vignette_params_2[int(1)].alpha1_0 - mean_3;
    float _S1386 = p_3.vignette_params_2[int(2)].alpha1_0 - mean_3;
    losses_0[int(8)] = (_S1384 * _S1384 + _S1385 * _S1385 + _S1386 * _S1386) / 3.0f;
    float mean_4 = (p_3.vignette_params_2[int(0)].alpha2_0 + p_3.vignette_params_2[int(1)].alpha2_0 + p_3.vignette_params_2[int(2)].alpha2_0) / 3.0f;
    float _S1387 = p_3.vignette_params_2[int(0)].alpha2_0 - mean_4;
    float _S1388 = p_3.vignette_params_2[int(1)].alpha2_0 - mean_4;
    float _S1389 = p_3.vignette_params_2[int(2)].alpha2_0 - mean_4;
    losses_0[int(9)] = (_S1387 * _S1387 + _S1388 * _S1388 + _S1389 * _S1389) / 3.0f;
    float2  bd_3 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_3.color_params_2.b_0);
    float2  rd_3 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_3.color_params_2.r_0);
    float2  gd_3 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_3.color_params_2.g_0);
    float2  nd_3 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_3.color_params_2.n_0);
    losses_0[int(10)] = bd_3.x;
    losses_0[int(11)] = bd_3.y;
    losses_0[int(12)] = rd_3.x;
    losses_0[int(13)] = rd_3.y;
    losses_0[int(14)] = gd_3.x;
    losses_0[int(15)] = gd_3.y;
    losses_0[int(16)] = nd_3.x;
    losses_0[int(17)] = nd_3.y;
    float mean_5 = (p_3.crf_params_1[int(0)].toe_0 + p_3.crf_params_1[int(1)].toe_0 + p_3.crf_params_1[int(2)].toe_0) / 3.0f;
    float _S1390 = p_3.crf_params_1[int(0)].toe_0 - mean_5;
    float _S1391 = p_3.crf_params_1[int(1)].toe_0 - mean_5;
    float _S1392 = p_3.crf_params_1[int(2)].toe_0 - mean_5;
    losses_0[int(18)] = (_S1390 * _S1390 + _S1391 * _S1391 + _S1392 * _S1392) / 3.0f;
    float mean_6 = (p_3.crf_params_1[int(0)].shoulder_0 + p_3.crf_params_1[int(1)].shoulder_0 + p_3.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S1393 = p_3.crf_params_1[int(0)].shoulder_0 - mean_6;
    float _S1394 = p_3.crf_params_1[int(1)].shoulder_0 - mean_6;
    float _S1395 = p_3.crf_params_1[int(2)].shoulder_0 - mean_6;
    losses_0[int(19)] = (_S1393 * _S1393 + _S1394 * _S1394 + _S1395 * _S1395) / 3.0f;
    float mean_7 = (p_3.crf_params_1[int(0)].gamma_0 + p_3.crf_params_1[int(1)].gamma_0 + p_3.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S1396 = p_3.crf_params_1[int(0)].gamma_0 - mean_7;
    float _S1397 = p_3.crf_params_1[int(1)].gamma_0 - mean_7;
    float _S1398 = p_3.crf_params_1[int(2)].gamma_0 - mean_7;
    losses_0[int(20)] = (_S1396 * _S1396 + _S1397 * _S1397 + _S1398 * _S1398) / 3.0f;
    float mean_8 = (p_3.crf_params_1[int(0)].center_0 + p_3.crf_params_1[int(1)].center_0 + p_3.crf_params_1[int(2)].center_0) / 3.0f;
    float _S1399 = p_3.crf_params_1[int(0)].center_0 - mean_8;
    float _S1400 = p_3.crf_params_1[int(1)].center_0 - mean_8;
    float _S1401 = p_3.crf_params_1[int(2)].center_0 - mean_8;
    losses_0[int(21)] = (_S1399 * _S1399 + _S1400 * _S1400 + _S1401 * _S1401) / 3.0f;
    *_S1368 = losses_0;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_7, FixedArray<float, 23>  * _S1402)
{
    PPISPParamsRQS_0 p_4;
    (&p_4)->exposure_1 = params_7[int(0)];
    (&(&p_4)->vignette_params_1[int(0)])->cx_0 = params_7[int(1)];
    (&(&p_4)->vignette_params_1[int(0)])->cy_0 = params_7[int(2)];
    (&(&p_4)->vignette_params_1[int(0)])->alpha0_0 = params_7[int(3)];
    (&(&p_4)->vignette_params_1[int(0)])->alpha1_0 = params_7[int(4)];
    (&(&p_4)->vignette_params_1[int(0)])->alpha2_0 = params_7[int(5)];
    (&(&p_4)->vignette_params_1[int(1)])->cx_0 = params_7[int(6)];
    (&(&p_4)->vignette_params_1[int(1)])->cy_0 = params_7[int(7)];
    (&(&p_4)->vignette_params_1[int(1)])->alpha0_0 = params_7[int(8)];
    (&(&p_4)->vignette_params_1[int(1)])->alpha1_0 = params_7[int(9)];
    (&(&p_4)->vignette_params_1[int(1)])->alpha2_0 = params_7[int(10)];
    (&(&p_4)->vignette_params_1[int(2)])->cx_0 = params_7[int(11)];
    (&(&p_4)->vignette_params_1[int(2)])->cy_0 = params_7[int(12)];
    (&(&p_4)->vignette_params_1[int(2)])->alpha0_0 = params_7[int(13)];
    (&(&p_4)->vignette_params_1[int(2)])->alpha1_0 = params_7[int(14)];
    (&(&p_4)->vignette_params_1[int(2)])->alpha2_0 = params_7[int(15)];
    *&((&(&(&p_4)->color_params_1)->b_0)->x) = params_7[int(16)];
    *&((&(&(&p_4)->color_params_1)->b_0)->y) = params_7[int(17)];
    *&((&(&(&p_4)->color_params_1)->r_0)->x) = params_7[int(18)];
    *&((&(&(&p_4)->color_params_1)->r_0)->y) = params_7[int(19)];
    *&((&(&(&p_4)->color_params_1)->g_0)->x) = params_7[int(20)];
    *&((&(&(&p_4)->color_params_1)->g_0)->y) = params_7[int(21)];
    *&((&(&(&p_4)->color_params_1)->n_0)->x) = params_7[int(22)];
    *&((&(&(&p_4)->color_params_1)->n_0)->y) = params_7[int(23)];
    (&(&p_4)->crf_params_0[int(0)])->g0_0 = params_7[int(24)];
    (&(&p_4)->crf_params_0[int(0)])->g1_0 = params_7[int(25)];
    (&(&p_4)->crf_params_0[int(0)])->x0_0 = params_7[int(26)];
    (&(&p_4)->crf_params_0[int(0)])->y0_0 = params_7[int(27)];
    (&(&p_4)->crf_params_0[int(0)])->gc_0 = params_7[int(28)];
    (&(&p_4)->crf_params_0[int(1)])->g0_0 = params_7[int(29)];
    (&(&p_4)->crf_params_0[int(1)])->g1_0 = params_7[int(30)];
    (&(&p_4)->crf_params_0[int(1)])->x0_0 = params_7[int(31)];
    (&(&p_4)->crf_params_0[int(1)])->y0_0 = params_7[int(32)];
    (&(&p_4)->crf_params_0[int(1)])->gc_0 = params_7[int(33)];
    (&(&p_4)->crf_params_0[int(2)])->g0_0 = params_7[int(34)];
    (&(&p_4)->crf_params_0[int(2)])->g1_0 = params_7[int(35)];
    (&(&p_4)->crf_params_0[int(2)])->x0_0 = params_7[int(36)];
    (&(&p_4)->crf_params_0[int(2)])->y0_0 = params_7[int(37)];
    (&(&p_4)->crf_params_0[int(2)])->gc_0 = params_7[int(38)];
    FixedArray<float, 23>  losses_1;
    losses_1[int(0)] = 0.0f;
    losses_1[int(1)] = 0.0f;
    losses_1[int(2)] = 0.0f;
    losses_1[int(3)] = 0.0f;
    losses_1[int(4)] = 0.0f;
    losses_1[int(5)] = 0.0f;
    losses_1[int(6)] = 0.0f;
    losses_1[int(7)] = 0.0f;
    losses_1[int(8)] = 0.0f;
    losses_1[int(9)] = 0.0f;
    losses_1[int(10)] = 0.0f;
    losses_1[int(11)] = 0.0f;
    losses_1[int(12)] = 0.0f;
    losses_1[int(13)] = 0.0f;
    losses_1[int(14)] = 0.0f;
    losses_1[int(15)] = 0.0f;
    losses_1[int(16)] = 0.0f;
    losses_1[int(17)] = 0.0f;
    losses_1[int(18)] = 0.0f;
    losses_1[int(19)] = 0.0f;
    losses_1[int(20)] = 0.0f;
    losses_1[int(21)] = 0.0f;
    losses_1[int(22)] = 0.0f;
    losses_1[int(0)] = p_4.exposure_1;
    float _S1403 = p_4.vignette_params_1[int(0)].cx_0;
    float _S1404 = p_4.vignette_params_1[int(0)].cy_0;
    float _S1405 = p_4.vignette_params_1[int(1)].cx_0;
    float _S1406 = p_4.vignette_params_1[int(1)].cy_0;
    float _S1407 = p_4.vignette_params_1[int(2)].cx_0;
    float _S1408 = p_4.vignette_params_1[int(2)].cy_0;
    losses_1[int(1)] = _S1403 * _S1403 + _S1404 * _S1404 + _S1405 * _S1405 + _S1406 * _S1406 + _S1407 * _S1407 + _S1408 * _S1408;
    losses_1[int(2)] = (F32_max((0.0f), (p_4.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(2)].alpha0_0)));
    losses_1[int(3)] = (F32_max((0.0f), (p_4.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(2)].alpha1_0)));
    losses_1[int(4)] = (F32_max((0.0f), (p_4.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_4.vignette_params_1[int(2)].alpha2_0)));
    float mean_9 = (p_4.vignette_params_1[int(0)].cx_0 + p_4.vignette_params_1[int(1)].cx_0 + p_4.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S1409 = p_4.vignette_params_1[int(0)].cx_0 - mean_9;
    float _S1410 = p_4.vignette_params_1[int(1)].cx_0 - mean_9;
    float _S1411 = p_4.vignette_params_1[int(2)].cx_0 - mean_9;
    losses_1[int(5)] = (_S1409 * _S1409 + _S1410 * _S1410 + _S1411 * _S1411) / 3.0f;
    float mean_10 = (p_4.vignette_params_1[int(0)].cy_0 + p_4.vignette_params_1[int(1)].cy_0 + p_4.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S1412 = p_4.vignette_params_1[int(0)].cy_0 - mean_10;
    float _S1413 = p_4.vignette_params_1[int(1)].cy_0 - mean_10;
    float _S1414 = p_4.vignette_params_1[int(2)].cy_0 - mean_10;
    losses_1[int(6)] = (_S1412 * _S1412 + _S1413 * _S1413 + _S1414 * _S1414) / 3.0f;
    float mean_11 = (p_4.vignette_params_1[int(0)].alpha0_0 + p_4.vignette_params_1[int(1)].alpha0_0 + p_4.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S1415 = p_4.vignette_params_1[int(0)].alpha0_0 - mean_11;
    float _S1416 = p_4.vignette_params_1[int(1)].alpha0_0 - mean_11;
    float _S1417 = p_4.vignette_params_1[int(2)].alpha0_0 - mean_11;
    losses_1[int(7)] = (_S1415 * _S1415 + _S1416 * _S1416 + _S1417 * _S1417) / 3.0f;
    float mean_12 = (p_4.vignette_params_1[int(0)].alpha1_0 + p_4.vignette_params_1[int(1)].alpha1_0 + p_4.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S1418 = p_4.vignette_params_1[int(0)].alpha1_0 - mean_12;
    float _S1419 = p_4.vignette_params_1[int(1)].alpha1_0 - mean_12;
    float _S1420 = p_4.vignette_params_1[int(2)].alpha1_0 - mean_12;
    losses_1[int(8)] = (_S1418 * _S1418 + _S1419 * _S1419 + _S1420 * _S1420) / 3.0f;
    float mean_13 = (p_4.vignette_params_1[int(0)].alpha2_0 + p_4.vignette_params_1[int(1)].alpha2_0 + p_4.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S1421 = p_4.vignette_params_1[int(0)].alpha2_0 - mean_13;
    float _S1422 = p_4.vignette_params_1[int(1)].alpha2_0 - mean_13;
    float _S1423 = p_4.vignette_params_1[int(2)].alpha2_0 - mean_13;
    losses_1[int(9)] = (_S1421 * _S1421 + _S1422 * _S1422 + _S1423 * _S1423) / 3.0f;
    float2  bd_4 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_4.color_params_1.b_0);
    float2  rd_4 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_4.color_params_1.r_0);
    float2  gd_4 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_4.color_params_1.g_0);
    float2  nd_4 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_4.color_params_1.n_0);
    losses_1[int(10)] = bd_4.x;
    losses_1[int(11)] = bd_4.y;
    losses_1[int(12)] = rd_4.x;
    losses_1[int(13)] = rd_4.y;
    losses_1[int(14)] = gd_4.x;
    losses_1[int(15)] = gd_4.y;
    losses_1[int(16)] = nd_4.x;
    losses_1[int(17)] = nd_4.y;
    float mean_14 = (p_4.crf_params_0[int(0)].g0_0 + p_4.crf_params_0[int(1)].g0_0 + p_4.crf_params_0[int(2)].g0_0) / 3.0f;
    float _S1424 = p_4.crf_params_0[int(0)].g0_0 - mean_14;
    float _S1425 = p_4.crf_params_0[int(1)].g0_0 - mean_14;
    float _S1426 = p_4.crf_params_0[int(2)].g0_0 - mean_14;
    losses_1[int(18)] = (_S1424 * _S1424 + _S1425 * _S1425 + _S1426 * _S1426) / 3.0f;
    float mean_15 = (p_4.crf_params_0[int(0)].g1_0 + p_4.crf_params_0[int(1)].g1_0 + p_4.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S1427 = p_4.crf_params_0[int(0)].g1_0 - mean_15;
    float _S1428 = p_4.crf_params_0[int(1)].g1_0 - mean_15;
    float _S1429 = p_4.crf_params_0[int(2)].g1_0 - mean_15;
    losses_1[int(19)] = (_S1427 * _S1427 + _S1428 * _S1428 + _S1429 * _S1429) / 3.0f;
    float mean_16 = (p_4.crf_params_0[int(0)].x0_0 + p_4.crf_params_0[int(1)].x0_0 + p_4.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S1430 = p_4.crf_params_0[int(0)].x0_0 - mean_16;
    float _S1431 = p_4.crf_params_0[int(1)].x0_0 - mean_16;
    float _S1432 = p_4.crf_params_0[int(2)].x0_0 - mean_16;
    losses_1[int(20)] = (_S1430 * _S1430 + _S1431 * _S1431 + _S1432 * _S1432) / 3.0f;
    float mean_17 = (p_4.crf_params_0[int(0)].y0_0 + p_4.crf_params_0[int(1)].y0_0 + p_4.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S1433 = p_4.crf_params_0[int(0)].y0_0 - mean_17;
    float _S1434 = p_4.crf_params_0[int(1)].y0_0 - mean_17;
    float _S1435 = p_4.crf_params_0[int(2)].y0_0 - mean_17;
    losses_1[int(21)] = (_S1433 * _S1433 + _S1434 * _S1434 + _S1435 * _S1435) / 3.0f;
    float mean_18 = (p_4.crf_params_0[int(0)].gc_0 + p_4.crf_params_0[int(1)].gc_0 + p_4.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S1436 = p_4.crf_params_0[int(0)].gc_0 - mean_18;
    float _S1437 = p_4.crf_params_0[int(1)].gc_0 - mean_18;
    float _S1438 = p_4.crf_params_0[int(2)].gc_0 - mean_18;
    losses_1[int(22)] = (_S1436 * _S1436 + _S1437 * _S1437 + _S1438 * _S1438) / 3.0f;
    *_S1402 = losses_1;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_3, FixedArray<float, 22>  * _s_dOut_3)
{
    VignettingChannelParams_0 _S1439 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1440 = {
        _S1439, _S1439, _S1439
    };
    float2  _S1441 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1442 = { _S1441, _S1441, _S1441, _S1441 };
    CRFPPISPChannelParams_0 _S1443 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1444 = {
        _S1443, _S1443, _S1443
    };
    PPISPParams_0 _S1445;
    (&_S1445)->exposure_2 = dpparams_3->primal_0[int(0)];
    (&_S1445)->vignette_params_2 = _S1440;
    (&_S1445)->color_params_2 = _S1442;
    (&_S1445)->crf_params_1 = _S1444;
    (&(&_S1445)->vignette_params_2[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S1445)->vignette_params_2[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S1445)->vignette_params_2[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S1445)->vignette_params_2[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S1445)->vignette_params_2[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S1445)->vignette_params_2[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S1445)->vignette_params_2[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S1445)->vignette_params_2[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S1445)->vignette_params_2[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S1445)->vignette_params_2[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S1445)->vignette_params_2[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S1445)->vignette_params_2[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S1445)->vignette_params_2[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S1445)->vignette_params_2[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S1445)->vignette_params_2[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S1445)->color_params_2)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S1445)->color_params_2)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S1445)->color_params_2)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S1445)->color_params_2)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S1445)->color_params_2)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S1445)->color_params_2)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S1445)->color_params_2)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S1445)->color_params_2)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S1445)->crf_params_1[int(0)])->toe_0 = dpparams_3->primal_0[int(24)];
    (&(&_S1445)->crf_params_1[int(0)])->shoulder_0 = dpparams_3->primal_0[int(25)];
    (&(&_S1445)->crf_params_1[int(0)])->gamma_0 = dpparams_3->primal_0[int(26)];
    (&(&_S1445)->crf_params_1[int(0)])->center_0 = dpparams_3->primal_0[int(27)];
    (&(&_S1445)->crf_params_1[int(1)])->toe_0 = dpparams_3->primal_0[int(28)];
    (&(&_S1445)->crf_params_1[int(1)])->shoulder_0 = dpparams_3->primal_0[int(29)];
    (&(&_S1445)->crf_params_1[int(1)])->gamma_0 = dpparams_3->primal_0[int(30)];
    (&(&_S1445)->crf_params_1[int(1)])->center_0 = dpparams_3->primal_0[int(31)];
    (&(&_S1445)->crf_params_1[int(2)])->toe_0 = dpparams_3->primal_0[int(32)];
    (&(&_S1445)->crf_params_1[int(2)])->shoulder_0 = dpparams_3->primal_0[int(33)];
    (&(&_S1445)->crf_params_1[int(2)])->gamma_0 = dpparams_3->primal_0[int(34)];
    (&(&_S1445)->crf_params_1[int(2)])->center_0 = dpparams_3->primal_0[int(35)];
    float mean_19 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S1446 = dpparams_3->primal_0[int(1)] - mean_19;
    float _S1447 = dpparams_3->primal_0[int(6)] - mean_19;
    float _S1448 = dpparams_3->primal_0[int(11)] - mean_19;
    float mean_20 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S1449 = dpparams_3->primal_0[int(2)] - mean_20;
    float _S1450 = dpparams_3->primal_0[int(7)] - mean_20;
    float _S1451 = dpparams_3->primal_0[int(12)] - mean_20;
    float mean_21 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S1452 = dpparams_3->primal_0[int(3)] - mean_21;
    float _S1453 = dpparams_3->primal_0[int(8)] - mean_21;
    float _S1454 = dpparams_3->primal_0[int(13)] - mean_21;
    float mean_22 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S1455 = dpparams_3->primal_0[int(4)] - mean_22;
    float _S1456 = dpparams_3->primal_0[int(9)] - mean_22;
    float _S1457 = dpparams_3->primal_0[int(14)] - mean_22;
    float mean_23 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S1458 = dpparams_3->primal_0[int(5)] - mean_23;
    float _S1459 = dpparams_3->primal_0[int(10)] - mean_23;
    float _S1460 = dpparams_3->primal_0[int(15)] - mean_23;
    float mean_24 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(32)]) / 3.0f;
    float mean_25 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(33)]) / 3.0f;
    float mean_26 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_27 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float _S1461 = 0.3333333432674408f * (*_s_dOut_3)[int(21)];
    float _S1462 = (dpparams_3->primal_0[int(35)] - mean_27) * _S1461;
    float _S1463 = _S1462 + _S1462;
    float _S1464 = (dpparams_3->primal_0[int(31)] - mean_27) * _S1461;
    float _S1465 = _S1464 + _S1464;
    float _S1466 = (dpparams_3->primal_0[int(27)] - mean_27) * _S1461;
    float _S1467 = _S1466 + _S1466;
    float _S1468 = 0.3333333432674408f * (- _S1463 + - _S1465 + - _S1467);
    float _S1469 = 0.3333333432674408f * (*_s_dOut_3)[int(20)];
    float _S1470 = (dpparams_3->primal_0[int(34)] - mean_26) * _S1469;
    float _S1471 = _S1470 + _S1470;
    float _S1472 = (dpparams_3->primal_0[int(30)] - mean_26) * _S1469;
    float _S1473 = _S1472 + _S1472;
    float _S1474 = (dpparams_3->primal_0[int(26)] - mean_26) * _S1469;
    float _S1475 = _S1474 + _S1474;
    float _S1476 = 0.3333333432674408f * (- _S1471 + - _S1473 + - _S1475);
    float _S1477 = 0.3333333432674408f * (*_s_dOut_3)[int(19)];
    float _S1478 = (dpparams_3->primal_0[int(33)] - mean_25) * _S1477;
    float _S1479 = _S1478 + _S1478;
    float _S1480 = (dpparams_3->primal_0[int(29)] - mean_25) * _S1477;
    float _S1481 = _S1480 + _S1480;
    float _S1482 = (dpparams_3->primal_0[int(25)] - mean_25) * _S1477;
    float _S1483 = _S1482 + _S1482;
    float _S1484 = 0.3333333432674408f * (- _S1479 + - _S1481 + - _S1483);
    float _S1485 = 0.3333333432674408f * (*_s_dOut_3)[int(18)];
    float _S1486 = (dpparams_3->primal_0[int(32)] - mean_24) * _S1485;
    float _S1487 = _S1486 + _S1486;
    float _S1488 = (dpparams_3->primal_0[int(28)] - mean_24) * _S1485;
    float _S1489 = _S1488 + _S1488;
    float _S1490 = (dpparams_3->primal_0[int(24)] - mean_24) * _S1485;
    float _S1491 = _S1490 + _S1490;
    float _S1492 = 0.3333333432674408f * (- _S1487 + - _S1489 + - _S1491);
    float2  _S1493 = make_float2 ((*_s_dOut_3)[int(16)], (*_s_dOut_3)[int(17)]);
    Matrix<float, 2, 2>  _S1494 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1495;
    (&_S1495)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1495)->differential_0 = _S1494;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1496;
    (&_S1496)->primal_0 = _S1445.color_params_2.n_0;
    (&_S1496)->differential_0 = _S1441;
    s_bwd_prop_mul_2(&_S1495, &_S1496, _S1493);
    float2  _S1497 = make_float2 ((*_s_dOut_3)[int(14)], (*_s_dOut_3)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1498;
    (&_S1498)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1498)->differential_0 = _S1494;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1499;
    (&_S1499)->primal_0 = _S1445.color_params_2.g_0;
    (&_S1499)->differential_0 = _S1441;
    s_bwd_prop_mul_2(&_S1498, &_S1499, _S1497);
    float2  _S1500 = make_float2 ((*_s_dOut_3)[int(12)], (*_s_dOut_3)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1501;
    (&_S1501)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1501)->differential_0 = _S1494;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1502;
    (&_S1502)->primal_0 = _S1445.color_params_2.r_0;
    (&_S1502)->differential_0 = _S1441;
    s_bwd_prop_mul_2(&_S1501, &_S1502, _S1500);
    float2  _S1503 = make_float2 ((*_s_dOut_3)[int(10)], (*_s_dOut_3)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1504;
    (&_S1504)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1504)->differential_0 = _S1494;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1505;
    (&_S1505)->primal_0 = _S1445.color_params_2.b_0;
    (&_S1505)->differential_0 = _S1441;
    s_bwd_prop_mul_2(&_S1504, &_S1505, _S1503);
    ColorPPISPParams_0 _S1506 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1506)->n_0 = _S1496.differential_0;
    (&_S1506)->g_0 = _S1499.differential_0;
    (&_S1506)->r_0 = _S1502.differential_0;
    (&_S1506)->b_0 = _S1505.differential_0;
    float _S1507 = 0.3333333432674408f * (*_s_dOut_3)[int(9)];
    float _S1508 = _S1460 * _S1507;
    float _S1509 = _S1508 + _S1508;
    float _S1510 = _S1459 * _S1507;
    float _S1511 = _S1510 + _S1510;
    float _S1512 = _S1458 * _S1507;
    float _S1513 = _S1512 + _S1512;
    float _S1514 = 0.3333333432674408f * (- _S1509 + - _S1511 + - _S1513);
    float _S1515 = 0.3333333432674408f * (*_s_dOut_3)[int(8)];
    float _S1516 = _S1457 * _S1515;
    float _S1517 = _S1516 + _S1516;
    float _S1518 = _S1456 * _S1515;
    float _S1519 = _S1518 + _S1518;
    float _S1520 = _S1455 * _S1515;
    float _S1521 = _S1520 + _S1520;
    float _S1522 = 0.3333333432674408f * (- _S1517 + - _S1519 + - _S1521);
    float _S1523 = 0.3333333432674408f * (*_s_dOut_3)[int(7)];
    float _S1524 = _S1454 * _S1523;
    float _S1525 = _S1524 + _S1524;
    float _S1526 = _S1453 * _S1523;
    float _S1527 = _S1526 + _S1526;
    float _S1528 = _S1452 * _S1523;
    float _S1529 = _S1528 + _S1528;
    float _S1530 = 0.3333333432674408f * (- _S1525 + - _S1527 + - _S1529);
    float _S1531 = 0.3333333432674408f * (*_s_dOut_3)[int(6)];
    float _S1532 = _S1451 * _S1531;
    float _S1533 = _S1532 + _S1532;
    float _S1534 = _S1450 * _S1531;
    float _S1535 = _S1534 + _S1534;
    float _S1536 = _S1449 * _S1531;
    float _S1537 = _S1536 + _S1536;
    float _S1538 = 0.3333333432674408f * (- _S1533 + - _S1535 + - _S1537);
    float _S1539 = 0.3333333432674408f * (*_s_dOut_3)[int(5)];
    float _S1540 = _S1448 * _S1539;
    float _S1541 = _S1540 + _S1540;
    float _S1542 = _S1447 * _S1539;
    float _S1543 = _S1542 + _S1542;
    float _S1544 = _S1446 * _S1539;
    float _S1545 = _S1544 + _S1544;
    float _S1546 = 0.3333333432674408f * (- _S1541 + - _S1543 + - _S1545);
    DiffPair_float_0 _S1547;
    (&_S1547)->primal_0 = 0.0f;
    (&_S1547)->differential_0 = 0.0f;
    DiffPair_float_0 _S1548;
    (&_S1548)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S1548)->differential_0 = 0.0f;
    _d_max_0(&_S1547, &_S1548, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1549;
    (&_S1549)->primal_0 = 0.0f;
    (&_S1549)->differential_0 = 0.0f;
    DiffPair_float_0 _S1550;
    (&_S1550)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S1550)->differential_0 = 0.0f;
    _d_max_0(&_S1549, &_S1550, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1551;
    (&_S1551)->primal_0 = 0.0f;
    (&_S1551)->differential_0 = 0.0f;
    DiffPair_float_0 _S1552;
    (&_S1552)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S1552)->differential_0 = 0.0f;
    _d_max_0(&_S1551, &_S1552, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1553;
    (&_S1553)->primal_0 = 0.0f;
    (&_S1553)->differential_0 = 0.0f;
    DiffPair_float_0 _S1554;
    (&_S1554)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S1554)->differential_0 = 0.0f;
    _d_max_0(&_S1553, &_S1554, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1555;
    (&_S1555)->primal_0 = 0.0f;
    (&_S1555)->differential_0 = 0.0f;
    DiffPair_float_0 _S1556;
    (&_S1556)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S1556)->differential_0 = 0.0f;
    _d_max_0(&_S1555, &_S1556, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1557;
    (&_S1557)->primal_0 = 0.0f;
    (&_S1557)->differential_0 = 0.0f;
    DiffPair_float_0 _S1558;
    (&_S1558)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S1558)->differential_0 = 0.0f;
    _d_max_0(&_S1557, &_S1558, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1559;
    (&_S1559)->primal_0 = 0.0f;
    (&_S1559)->differential_0 = 0.0f;
    DiffPair_float_0 _S1560;
    (&_S1560)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S1560)->differential_0 = 0.0f;
    _d_max_0(&_S1559, &_S1560, (*_s_dOut_3)[int(2)]);
    DiffPair_float_0 _S1561;
    (&_S1561)->primal_0 = 0.0f;
    (&_S1561)->differential_0 = 0.0f;
    DiffPair_float_0 _S1562;
    (&_S1562)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S1562)->differential_0 = 0.0f;
    _d_max_0(&_S1561, &_S1562, (*_s_dOut_3)[int(2)]);
    DiffPair_float_0 _S1563;
    (&_S1563)->primal_0 = 0.0f;
    (&_S1563)->differential_0 = 0.0f;
    DiffPair_float_0 _S1564;
    (&_S1564)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S1564)->differential_0 = 0.0f;
    _d_max_0(&_S1563, &_S1564, (*_s_dOut_3)[int(2)]);
    float _S1565 = dpparams_3->primal_0[int(12)] * (*_s_dOut_3)[int(1)];
    float _S1566 = dpparams_3->primal_0[int(11)] * (*_s_dOut_3)[int(1)];
    float _S1567 = dpparams_3->primal_0[int(7)] * (*_s_dOut_3)[int(1)];
    float _S1568 = dpparams_3->primal_0[int(6)] * (*_s_dOut_3)[int(1)];
    float _S1569 = dpparams_3->primal_0[int(2)] * (*_s_dOut_3)[int(1)];
    float _S1570 = dpparams_3->primal_0[int(1)] * (*_s_dOut_3)[int(1)];
    PPISPParams_0 _S1571 = PPISPParams_x24_syn_dzero_0();
    (&_S1571)->color_params_2 = _S1506;
    (&_S1571)->exposure_2 = (*_s_dOut_3)[int(0)];
    _S1445 = _S1571;
    (&(&_S1445)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1572 = _S1463 + _S1468 + _S1571.crf_params_1[int(2)].center_0;
    (&(&_S1445)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1573 = _S1471 + _S1476 + _S1571.crf_params_1[int(2)].gamma_0;
    (&(&_S1445)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1574 = _S1479 + _S1484 + _S1571.crf_params_1[int(2)].shoulder_0;
    (&(&_S1445)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1575 = _S1487 + _S1492 + _S1571.crf_params_1[int(2)].toe_0;
    (&(&_S1445)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1576 = _S1465 + _S1468 + _S1571.crf_params_1[int(1)].center_0;
    (&(&_S1445)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1577 = _S1473 + _S1476 + _S1571.crf_params_1[int(1)].gamma_0;
    (&(&_S1445)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1578 = _S1481 + _S1484 + _S1571.crf_params_1[int(1)].shoulder_0;
    (&(&_S1445)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1579 = _S1489 + _S1492 + _S1571.crf_params_1[int(1)].toe_0;
    (&(&_S1445)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1580 = _S1467 + _S1468 + _S1571.crf_params_1[int(0)].center_0;
    (&(&_S1445)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1581 = _S1475 + _S1476 + _S1571.crf_params_1[int(0)].gamma_0;
    (&(&_S1445)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1582 = _S1483 + _S1484 + _S1571.crf_params_1[int(0)].shoulder_0;
    (&(&_S1445)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1583 = _S1491 + _S1492 + _S1571.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S1445)->color_params_2)->n_0)->y) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->n_0)->x) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->g_0)->y) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->g_0)->x) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->r_0)->y) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->r_0)->x) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->b_0)->y) = 0.0f;
    *&((&(&(&_S1445)->color_params_2)->b_0)->x) = 0.0f;
    (&(&_S1445)->vignette_params_2[int(2)])->alpha2_0 = 0.0f;
    float _S1584 = _S1509 + _S1514 + _S1548.differential_0 + _S1571.vignette_params_2[int(2)].alpha2_0;
    (&(&_S1445)->vignette_params_2[int(2)])->alpha1_0 = 0.0f;
    float _S1585 = _S1517 + _S1522 + _S1554.differential_0 + _S1571.vignette_params_2[int(2)].alpha1_0;
    (&(&_S1445)->vignette_params_2[int(2)])->alpha0_0 = 0.0f;
    float _S1586 = _S1525 + _S1530 + _S1560.differential_0 + _S1571.vignette_params_2[int(2)].alpha0_0;
    (&(&_S1445)->vignette_params_2[int(2)])->cy_0 = 0.0f;
    float _S1587 = _S1533 + _S1538 + _S1565 + _S1565 + _S1571.vignette_params_2[int(2)].cy_0;
    (&(&_S1445)->vignette_params_2[int(2)])->cx_0 = 0.0f;
    float _S1588 = _S1541 + _S1546 + _S1566 + _S1566 + _S1571.vignette_params_2[int(2)].cx_0;
    (&(&_S1445)->vignette_params_2[int(1)])->alpha2_0 = 0.0f;
    float _S1589 = _S1511 + _S1514 + _S1550.differential_0 + _S1571.vignette_params_2[int(1)].alpha2_0;
    (&(&_S1445)->vignette_params_2[int(1)])->alpha1_0 = 0.0f;
    float _S1590 = _S1519 + _S1522 + _S1556.differential_0 + _S1571.vignette_params_2[int(1)].alpha1_0;
    (&(&_S1445)->vignette_params_2[int(1)])->alpha0_0 = 0.0f;
    float _S1591 = _S1527 + _S1530 + _S1562.differential_0 + _S1571.vignette_params_2[int(1)].alpha0_0;
    (&(&_S1445)->vignette_params_2[int(1)])->cy_0 = 0.0f;
    float _S1592 = _S1535 + _S1538 + _S1567 + _S1567 + _S1571.vignette_params_2[int(1)].cy_0;
    (&(&_S1445)->vignette_params_2[int(1)])->cx_0 = 0.0f;
    float _S1593 = _S1543 + _S1546 + _S1568 + _S1568 + _S1571.vignette_params_2[int(1)].cx_0;
    (&(&_S1445)->vignette_params_2[int(0)])->alpha2_0 = 0.0f;
    float _S1594 = _S1513 + _S1514 + _S1552.differential_0 + _S1571.vignette_params_2[int(0)].alpha2_0;
    (&(&_S1445)->vignette_params_2[int(0)])->alpha1_0 = 0.0f;
    float _S1595 = _S1521 + _S1522 + _S1558.differential_0 + _S1571.vignette_params_2[int(0)].alpha1_0;
    (&(&_S1445)->vignette_params_2[int(0)])->alpha0_0 = 0.0f;
    float _S1596 = _S1529 + _S1530 + _S1564.differential_0 + _S1571.vignette_params_2[int(0)].alpha0_0;
    (&(&_S1445)->vignette_params_2[int(0)])->cy_0 = 0.0f;
    float _S1597 = _S1537 + _S1538 + _S1569 + _S1569 + _S1571.vignette_params_2[int(0)].cy_0;
    (&(&_S1445)->vignette_params_2[int(0)])->cx_0 = 0.0f;
    float _S1598 = _S1545 + _S1546 + _S1570 + _S1570 + _S1571.vignette_params_2[int(0)].cx_0;
    FixedArray<float, 36>  _S1599;
    _S1599[int(0)] = 0.0f;
    _S1599[int(1)] = 0.0f;
    _S1599[int(2)] = 0.0f;
    _S1599[int(3)] = 0.0f;
    _S1599[int(4)] = 0.0f;
    _S1599[int(5)] = 0.0f;
    _S1599[int(6)] = 0.0f;
    _S1599[int(7)] = 0.0f;
    _S1599[int(8)] = 0.0f;
    _S1599[int(9)] = 0.0f;
    _S1599[int(10)] = 0.0f;
    _S1599[int(11)] = 0.0f;
    _S1599[int(12)] = 0.0f;
    _S1599[int(13)] = 0.0f;
    _S1599[int(14)] = 0.0f;
    _S1599[int(15)] = 0.0f;
    _S1599[int(16)] = 0.0f;
    _S1599[int(17)] = 0.0f;
    _S1599[int(18)] = 0.0f;
    _S1599[int(19)] = 0.0f;
    _S1599[int(20)] = 0.0f;
    _S1599[int(21)] = 0.0f;
    _S1599[int(22)] = 0.0f;
    _S1599[int(23)] = 0.0f;
    _S1599[int(24)] = 0.0f;
    _S1599[int(25)] = 0.0f;
    _S1599[int(26)] = 0.0f;
    _S1599[int(27)] = 0.0f;
    _S1599[int(28)] = 0.0f;
    _S1599[int(29)] = 0.0f;
    _S1599[int(30)] = 0.0f;
    _S1599[int(31)] = 0.0f;
    _S1599[int(32)] = 0.0f;
    _S1599[int(33)] = 0.0f;
    _S1599[int(34)] = 0.0f;
    _S1599[int(35)] = 0.0f;
    _S1599[int(8)] = _S1591;
    _S1599[int(16)] = _S1571.color_params_2.b_0.x;
    _S1599[int(15)] = _S1584;
    _S1599[int(14)] = _S1585;
    _S1599[int(13)] = _S1586;
    _S1599[int(12)] = _S1587;
    _S1599[int(11)] = _S1588;
    _S1599[int(10)] = _S1589;
    _S1599[int(9)] = _S1590;
    _S1599[int(17)] = _S1571.color_params_2.b_0.y;
    _S1599[int(7)] = _S1592;
    _S1599[int(6)] = _S1593;
    _S1599[int(5)] = _S1594;
    _S1599[int(4)] = _S1595;
    _S1599[int(3)] = _S1596;
    _S1599[int(2)] = _S1597;
    _S1599[int(1)] = _S1598;
    _S1599[int(0)] = _S1445.exposure_2;
    _S1599[int(26)] = _S1581;
    _S1599[int(34)] = _S1573;
    _S1599[int(33)] = _S1574;
    _S1599[int(32)] = _S1575;
    _S1599[int(31)] = _S1576;
    _S1599[int(30)] = _S1577;
    _S1599[int(29)] = _S1578;
    _S1599[int(28)] = _S1579;
    _S1599[int(27)] = _S1580;
    _S1599[int(35)] = _S1572;
    _S1599[int(25)] = _S1582;
    _S1599[int(24)] = _S1583;
    _S1599[int(23)] = _S1571.color_params_2.n_0.y;
    _S1599[int(22)] = _S1571.color_params_2.n_0.x;
    _S1599[int(21)] = _S1571.color_params_2.g_0.y;
    _S1599[int(20)] = _S1571.color_params_2.g_0.x;
    _S1599[int(19)] = _S1571.color_params_2.r_0.y;
    _S1599[int(18)] = _S1571.color_params_2.r_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S1599;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1600, FixedArray<float, 22>  * _S1601)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S1600, _S1601);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_8, FixedArray<float, 22>  grad_out_3, FixedArray<float, 36>  * _S1602)
{
    FixedArray<float, 36>  _S1603 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_8;
    (&dp_params_3)->differential_0 = _S1603;
    FixedArray<float, 22>  _S1604 = grad_out_3;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_3, &_S1604);
    *_S1602 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_4, FixedArray<float, 23>  * _s_dOut_4)
{
    VignettingChannelParams_0 _S1605 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1606 = {
        _S1605, _S1605, _S1605
    };
    float2  _S1607 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1608 = { _S1607, _S1607, _S1607, _S1607 };
    RQSCRFPPISPChannelParams_0 _S1609 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1610 = {
        _S1609, _S1609, _S1609
    };
    PPISPParamsRQS_0 _S1611;
    (&_S1611)->exposure_1 = dpparams_4->primal_0[int(0)];
    (&_S1611)->vignette_params_1 = _S1606;
    (&_S1611)->color_params_1 = _S1608;
    (&_S1611)->crf_params_0 = _S1610;
    (&(&_S1611)->vignette_params_1[int(0)])->cx_0 = dpparams_4->primal_0[int(1)];
    (&(&_S1611)->vignette_params_1[int(0)])->cy_0 = dpparams_4->primal_0[int(2)];
    (&(&_S1611)->vignette_params_1[int(0)])->alpha0_0 = dpparams_4->primal_0[int(3)];
    (&(&_S1611)->vignette_params_1[int(0)])->alpha1_0 = dpparams_4->primal_0[int(4)];
    (&(&_S1611)->vignette_params_1[int(0)])->alpha2_0 = dpparams_4->primal_0[int(5)];
    (&(&_S1611)->vignette_params_1[int(1)])->cx_0 = dpparams_4->primal_0[int(6)];
    (&(&_S1611)->vignette_params_1[int(1)])->cy_0 = dpparams_4->primal_0[int(7)];
    (&(&_S1611)->vignette_params_1[int(1)])->alpha0_0 = dpparams_4->primal_0[int(8)];
    (&(&_S1611)->vignette_params_1[int(1)])->alpha1_0 = dpparams_4->primal_0[int(9)];
    (&(&_S1611)->vignette_params_1[int(1)])->alpha2_0 = dpparams_4->primal_0[int(10)];
    (&(&_S1611)->vignette_params_1[int(2)])->cx_0 = dpparams_4->primal_0[int(11)];
    (&(&_S1611)->vignette_params_1[int(2)])->cy_0 = dpparams_4->primal_0[int(12)];
    (&(&_S1611)->vignette_params_1[int(2)])->alpha0_0 = dpparams_4->primal_0[int(13)];
    (&(&_S1611)->vignette_params_1[int(2)])->alpha1_0 = dpparams_4->primal_0[int(14)];
    (&(&_S1611)->vignette_params_1[int(2)])->alpha2_0 = dpparams_4->primal_0[int(15)];
    *&((&(&(&_S1611)->color_params_1)->b_0)->x) = dpparams_4->primal_0[int(16)];
    *&((&(&(&_S1611)->color_params_1)->b_0)->y) = dpparams_4->primal_0[int(17)];
    *&((&(&(&_S1611)->color_params_1)->r_0)->x) = dpparams_4->primal_0[int(18)];
    *&((&(&(&_S1611)->color_params_1)->r_0)->y) = dpparams_4->primal_0[int(19)];
    *&((&(&(&_S1611)->color_params_1)->g_0)->x) = dpparams_4->primal_0[int(20)];
    *&((&(&(&_S1611)->color_params_1)->g_0)->y) = dpparams_4->primal_0[int(21)];
    *&((&(&(&_S1611)->color_params_1)->n_0)->x) = dpparams_4->primal_0[int(22)];
    *&((&(&(&_S1611)->color_params_1)->n_0)->y) = dpparams_4->primal_0[int(23)];
    (&(&_S1611)->crf_params_0[int(0)])->g0_0 = dpparams_4->primal_0[int(24)];
    (&(&_S1611)->crf_params_0[int(0)])->g1_0 = dpparams_4->primal_0[int(25)];
    (&(&_S1611)->crf_params_0[int(0)])->x0_0 = dpparams_4->primal_0[int(26)];
    (&(&_S1611)->crf_params_0[int(0)])->y0_0 = dpparams_4->primal_0[int(27)];
    (&(&_S1611)->crf_params_0[int(0)])->gc_0 = dpparams_4->primal_0[int(28)];
    (&(&_S1611)->crf_params_0[int(1)])->g0_0 = dpparams_4->primal_0[int(29)];
    (&(&_S1611)->crf_params_0[int(1)])->g1_0 = dpparams_4->primal_0[int(30)];
    (&(&_S1611)->crf_params_0[int(1)])->x0_0 = dpparams_4->primal_0[int(31)];
    (&(&_S1611)->crf_params_0[int(1)])->y0_0 = dpparams_4->primal_0[int(32)];
    (&(&_S1611)->crf_params_0[int(1)])->gc_0 = dpparams_4->primal_0[int(33)];
    (&(&_S1611)->crf_params_0[int(2)])->g0_0 = dpparams_4->primal_0[int(34)];
    (&(&_S1611)->crf_params_0[int(2)])->g1_0 = dpparams_4->primal_0[int(35)];
    (&(&_S1611)->crf_params_0[int(2)])->x0_0 = dpparams_4->primal_0[int(36)];
    (&(&_S1611)->crf_params_0[int(2)])->y0_0 = dpparams_4->primal_0[int(37)];
    (&(&_S1611)->crf_params_0[int(2)])->gc_0 = dpparams_4->primal_0[int(38)];
    float mean_28 = (dpparams_4->primal_0[int(1)] + dpparams_4->primal_0[int(6)] + dpparams_4->primal_0[int(11)]) / 3.0f;
    float _S1612 = dpparams_4->primal_0[int(1)] - mean_28;
    float _S1613 = dpparams_4->primal_0[int(6)] - mean_28;
    float _S1614 = dpparams_4->primal_0[int(11)] - mean_28;
    float mean_29 = (dpparams_4->primal_0[int(2)] + dpparams_4->primal_0[int(7)] + dpparams_4->primal_0[int(12)]) / 3.0f;
    float _S1615 = dpparams_4->primal_0[int(2)] - mean_29;
    float _S1616 = dpparams_4->primal_0[int(7)] - mean_29;
    float _S1617 = dpparams_4->primal_0[int(12)] - mean_29;
    float mean_30 = (dpparams_4->primal_0[int(3)] + dpparams_4->primal_0[int(8)] + dpparams_4->primal_0[int(13)]) / 3.0f;
    float _S1618 = dpparams_4->primal_0[int(3)] - mean_30;
    float _S1619 = dpparams_4->primal_0[int(8)] - mean_30;
    float _S1620 = dpparams_4->primal_0[int(13)] - mean_30;
    float mean_31 = (dpparams_4->primal_0[int(4)] + dpparams_4->primal_0[int(9)] + dpparams_4->primal_0[int(14)]) / 3.0f;
    float _S1621 = dpparams_4->primal_0[int(4)] - mean_31;
    float _S1622 = dpparams_4->primal_0[int(9)] - mean_31;
    float _S1623 = dpparams_4->primal_0[int(14)] - mean_31;
    float mean_32 = (dpparams_4->primal_0[int(5)] + dpparams_4->primal_0[int(10)] + dpparams_4->primal_0[int(15)]) / 3.0f;
    float _S1624 = dpparams_4->primal_0[int(5)] - mean_32;
    float _S1625 = dpparams_4->primal_0[int(10)] - mean_32;
    float _S1626 = dpparams_4->primal_0[int(15)] - mean_32;
    float mean_33 = (dpparams_4->primal_0[int(24)] + dpparams_4->primal_0[int(29)] + dpparams_4->primal_0[int(34)]) / 3.0f;
    float mean_34 = (dpparams_4->primal_0[int(25)] + dpparams_4->primal_0[int(30)] + dpparams_4->primal_0[int(35)]) / 3.0f;
    float mean_35 = (dpparams_4->primal_0[int(26)] + dpparams_4->primal_0[int(31)] + dpparams_4->primal_0[int(36)]) / 3.0f;
    float mean_36 = (dpparams_4->primal_0[int(27)] + dpparams_4->primal_0[int(32)] + dpparams_4->primal_0[int(37)]) / 3.0f;
    float mean_37 = (dpparams_4->primal_0[int(28)] + dpparams_4->primal_0[int(33)] + dpparams_4->primal_0[int(38)]) / 3.0f;
    float _S1627 = 0.3333333432674408f * (*_s_dOut_4)[int(22)];
    float _S1628 = (dpparams_4->primal_0[int(38)] - mean_37) * _S1627;
    float _S1629 = _S1628 + _S1628;
    float _S1630 = (dpparams_4->primal_0[int(33)] - mean_37) * _S1627;
    float _S1631 = _S1630 + _S1630;
    float _S1632 = (dpparams_4->primal_0[int(28)] - mean_37) * _S1627;
    float _S1633 = _S1632 + _S1632;
    float _S1634 = 0.3333333432674408f * (- _S1629 + - _S1631 + - _S1633);
    float _S1635 = 0.3333333432674408f * (*_s_dOut_4)[int(21)];
    float _S1636 = (dpparams_4->primal_0[int(37)] - mean_36) * _S1635;
    float _S1637 = _S1636 + _S1636;
    float _S1638 = (dpparams_4->primal_0[int(32)] - mean_36) * _S1635;
    float _S1639 = _S1638 + _S1638;
    float _S1640 = (dpparams_4->primal_0[int(27)] - mean_36) * _S1635;
    float _S1641 = _S1640 + _S1640;
    float _S1642 = 0.3333333432674408f * (- _S1637 + - _S1639 + - _S1641);
    float _S1643 = 0.3333333432674408f * (*_s_dOut_4)[int(20)];
    float _S1644 = (dpparams_4->primal_0[int(36)] - mean_35) * _S1643;
    float _S1645 = _S1644 + _S1644;
    float _S1646 = (dpparams_4->primal_0[int(31)] - mean_35) * _S1643;
    float _S1647 = _S1646 + _S1646;
    float _S1648 = (dpparams_4->primal_0[int(26)] - mean_35) * _S1643;
    float _S1649 = _S1648 + _S1648;
    float _S1650 = 0.3333333432674408f * (- _S1645 + - _S1647 + - _S1649);
    float _S1651 = 0.3333333432674408f * (*_s_dOut_4)[int(19)];
    float _S1652 = (dpparams_4->primal_0[int(35)] - mean_34) * _S1651;
    float _S1653 = _S1652 + _S1652;
    float _S1654 = (dpparams_4->primal_0[int(30)] - mean_34) * _S1651;
    float _S1655 = _S1654 + _S1654;
    float _S1656 = (dpparams_4->primal_0[int(25)] - mean_34) * _S1651;
    float _S1657 = _S1656 + _S1656;
    float _S1658 = 0.3333333432674408f * (- _S1653 + - _S1655 + - _S1657);
    float _S1659 = 0.3333333432674408f * (*_s_dOut_4)[int(18)];
    float _S1660 = (dpparams_4->primal_0[int(34)] - mean_33) * _S1659;
    float _S1661 = _S1660 + _S1660;
    float _S1662 = (dpparams_4->primal_0[int(29)] - mean_33) * _S1659;
    float _S1663 = _S1662 + _S1662;
    float _S1664 = (dpparams_4->primal_0[int(24)] - mean_33) * _S1659;
    float _S1665 = _S1664 + _S1664;
    float _S1666 = 0.3333333432674408f * (- _S1661 + - _S1663 + - _S1665);
    float2  _S1667 = make_float2 ((*_s_dOut_4)[int(16)], (*_s_dOut_4)[int(17)]);
    Matrix<float, 2, 2>  _S1668 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1669;
    (&_S1669)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1669)->differential_0 = _S1668;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1670;
    (&_S1670)->primal_0 = _S1611.color_params_1.n_0;
    (&_S1670)->differential_0 = _S1607;
    s_bwd_prop_mul_2(&_S1669, &_S1670, _S1667);
    float2  _S1671 = make_float2 ((*_s_dOut_4)[int(14)], (*_s_dOut_4)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1672;
    (&_S1672)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1672)->differential_0 = _S1668;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1673;
    (&_S1673)->primal_0 = _S1611.color_params_1.g_0;
    (&_S1673)->differential_0 = _S1607;
    s_bwd_prop_mul_2(&_S1672, &_S1673, _S1671);
    float2  _S1674 = make_float2 ((*_s_dOut_4)[int(12)], (*_s_dOut_4)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1675;
    (&_S1675)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1675)->differential_0 = _S1668;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1676;
    (&_S1676)->primal_0 = _S1611.color_params_1.r_0;
    (&_S1676)->differential_0 = _S1607;
    s_bwd_prop_mul_2(&_S1675, &_S1676, _S1674);
    float2  _S1677 = make_float2 ((*_s_dOut_4)[int(10)], (*_s_dOut_4)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1678;
    (&_S1678)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1678)->differential_0 = _S1668;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1679;
    (&_S1679)->primal_0 = _S1611.color_params_1.b_0;
    (&_S1679)->differential_0 = _S1607;
    s_bwd_prop_mul_2(&_S1678, &_S1679, _S1677);
    ColorPPISPParams_0 _S1680 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1680)->n_0 = _S1670.differential_0;
    (&_S1680)->g_0 = _S1673.differential_0;
    (&_S1680)->r_0 = _S1676.differential_0;
    (&_S1680)->b_0 = _S1679.differential_0;
    float _S1681 = 0.3333333432674408f * (*_s_dOut_4)[int(9)];
    float _S1682 = _S1626 * _S1681;
    float _S1683 = _S1682 + _S1682;
    float _S1684 = _S1625 * _S1681;
    float _S1685 = _S1684 + _S1684;
    float _S1686 = _S1624 * _S1681;
    float _S1687 = _S1686 + _S1686;
    float _S1688 = 0.3333333432674408f * (- _S1683 + - _S1685 + - _S1687);
    float _S1689 = 0.3333333432674408f * (*_s_dOut_4)[int(8)];
    float _S1690 = _S1623 * _S1689;
    float _S1691 = _S1690 + _S1690;
    float _S1692 = _S1622 * _S1689;
    float _S1693 = _S1692 + _S1692;
    float _S1694 = _S1621 * _S1689;
    float _S1695 = _S1694 + _S1694;
    float _S1696 = 0.3333333432674408f * (- _S1691 + - _S1693 + - _S1695);
    float _S1697 = 0.3333333432674408f * (*_s_dOut_4)[int(7)];
    float _S1698 = _S1620 * _S1697;
    float _S1699 = _S1698 + _S1698;
    float _S1700 = _S1619 * _S1697;
    float _S1701 = _S1700 + _S1700;
    float _S1702 = _S1618 * _S1697;
    float _S1703 = _S1702 + _S1702;
    float _S1704 = 0.3333333432674408f * (- _S1699 + - _S1701 + - _S1703);
    float _S1705 = 0.3333333432674408f * (*_s_dOut_4)[int(6)];
    float _S1706 = _S1617 * _S1705;
    float _S1707 = _S1706 + _S1706;
    float _S1708 = _S1616 * _S1705;
    float _S1709 = _S1708 + _S1708;
    float _S1710 = _S1615 * _S1705;
    float _S1711 = _S1710 + _S1710;
    float _S1712 = 0.3333333432674408f * (- _S1707 + - _S1709 + - _S1711);
    float _S1713 = 0.3333333432674408f * (*_s_dOut_4)[int(5)];
    float _S1714 = _S1614 * _S1713;
    float _S1715 = _S1714 + _S1714;
    float _S1716 = _S1613 * _S1713;
    float _S1717 = _S1716 + _S1716;
    float _S1718 = _S1612 * _S1713;
    float _S1719 = _S1718 + _S1718;
    float _S1720 = 0.3333333432674408f * (- _S1715 + - _S1717 + - _S1719);
    DiffPair_float_0 _S1721;
    (&_S1721)->primal_0 = 0.0f;
    (&_S1721)->differential_0 = 0.0f;
    DiffPair_float_0 _S1722;
    (&_S1722)->primal_0 = dpparams_4->primal_0[int(15)];
    (&_S1722)->differential_0 = 0.0f;
    _d_max_0(&_S1721, &_S1722, (*_s_dOut_4)[int(4)]);
    DiffPair_float_0 _S1723;
    (&_S1723)->primal_0 = 0.0f;
    (&_S1723)->differential_0 = 0.0f;
    DiffPair_float_0 _S1724;
    (&_S1724)->primal_0 = dpparams_4->primal_0[int(10)];
    (&_S1724)->differential_0 = 0.0f;
    _d_max_0(&_S1723, &_S1724, (*_s_dOut_4)[int(4)]);
    DiffPair_float_0 _S1725;
    (&_S1725)->primal_0 = 0.0f;
    (&_S1725)->differential_0 = 0.0f;
    DiffPair_float_0 _S1726;
    (&_S1726)->primal_0 = dpparams_4->primal_0[int(5)];
    (&_S1726)->differential_0 = 0.0f;
    _d_max_0(&_S1725, &_S1726, (*_s_dOut_4)[int(4)]);
    DiffPair_float_0 _S1727;
    (&_S1727)->primal_0 = 0.0f;
    (&_S1727)->differential_0 = 0.0f;
    DiffPair_float_0 _S1728;
    (&_S1728)->primal_0 = dpparams_4->primal_0[int(14)];
    (&_S1728)->differential_0 = 0.0f;
    _d_max_0(&_S1727, &_S1728, (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S1729;
    (&_S1729)->primal_0 = 0.0f;
    (&_S1729)->differential_0 = 0.0f;
    DiffPair_float_0 _S1730;
    (&_S1730)->primal_0 = dpparams_4->primal_0[int(9)];
    (&_S1730)->differential_0 = 0.0f;
    _d_max_0(&_S1729, &_S1730, (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S1731;
    (&_S1731)->primal_0 = 0.0f;
    (&_S1731)->differential_0 = 0.0f;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = dpparams_4->primal_0[int(4)];
    (&_S1732)->differential_0 = 0.0f;
    _d_max_0(&_S1731, &_S1732, (*_s_dOut_4)[int(3)]);
    DiffPair_float_0 _S1733;
    (&_S1733)->primal_0 = 0.0f;
    (&_S1733)->differential_0 = 0.0f;
    DiffPair_float_0 _S1734;
    (&_S1734)->primal_0 = dpparams_4->primal_0[int(13)];
    (&_S1734)->differential_0 = 0.0f;
    _d_max_0(&_S1733, &_S1734, (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S1735;
    (&_S1735)->primal_0 = 0.0f;
    (&_S1735)->differential_0 = 0.0f;
    DiffPair_float_0 _S1736;
    (&_S1736)->primal_0 = dpparams_4->primal_0[int(8)];
    (&_S1736)->differential_0 = 0.0f;
    _d_max_0(&_S1735, &_S1736, (*_s_dOut_4)[int(2)]);
    DiffPair_float_0 _S1737;
    (&_S1737)->primal_0 = 0.0f;
    (&_S1737)->differential_0 = 0.0f;
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = dpparams_4->primal_0[int(3)];
    (&_S1738)->differential_0 = 0.0f;
    _d_max_0(&_S1737, &_S1738, (*_s_dOut_4)[int(2)]);
    float _S1739 = dpparams_4->primal_0[int(12)] * (*_s_dOut_4)[int(1)];
    float _S1740 = dpparams_4->primal_0[int(11)] * (*_s_dOut_4)[int(1)];
    float _S1741 = dpparams_4->primal_0[int(7)] * (*_s_dOut_4)[int(1)];
    float _S1742 = dpparams_4->primal_0[int(6)] * (*_s_dOut_4)[int(1)];
    float _S1743 = dpparams_4->primal_0[int(2)] * (*_s_dOut_4)[int(1)];
    float _S1744 = dpparams_4->primal_0[int(1)] * (*_s_dOut_4)[int(1)];
    PPISPParamsRQS_0 _S1745 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S1745)->color_params_1 = _S1680;
    (&_S1745)->exposure_1 = (*_s_dOut_4)[int(0)];
    _S1611 = _S1745;
    (&(&_S1611)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S1746 = _S1629 + _S1634 + _S1745.crf_params_0[int(2)].gc_0;
    (&(&_S1611)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S1747 = _S1637 + _S1642 + _S1745.crf_params_0[int(2)].y0_0;
    (&(&_S1611)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S1748 = _S1645 + _S1650 + _S1745.crf_params_0[int(2)].x0_0;
    (&(&_S1611)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S1749 = _S1653 + _S1658 + _S1745.crf_params_0[int(2)].g1_0;
    (&(&_S1611)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S1750 = _S1661 + _S1666 + _S1745.crf_params_0[int(2)].g0_0;
    (&(&_S1611)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S1751 = _S1631 + _S1634 + _S1745.crf_params_0[int(1)].gc_0;
    (&(&_S1611)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S1752 = _S1639 + _S1642 + _S1745.crf_params_0[int(1)].y0_0;
    (&(&_S1611)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S1753 = _S1647 + _S1650 + _S1745.crf_params_0[int(1)].x0_0;
    (&(&_S1611)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S1754 = _S1655 + _S1658 + _S1745.crf_params_0[int(1)].g1_0;
    (&(&_S1611)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S1755 = _S1663 + _S1666 + _S1745.crf_params_0[int(1)].g0_0;
    (&(&_S1611)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S1756 = _S1633 + _S1634 + _S1745.crf_params_0[int(0)].gc_0;
    (&(&_S1611)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S1757 = _S1641 + _S1642 + _S1745.crf_params_0[int(0)].y0_0;
    (&(&_S1611)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S1758 = _S1649 + _S1650 + _S1745.crf_params_0[int(0)].x0_0;
    (&(&_S1611)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S1759 = _S1657 + _S1658 + _S1745.crf_params_0[int(0)].g1_0;
    (&(&_S1611)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S1760 = _S1665 + _S1666 + _S1745.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S1611)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1611)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1611)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1761 = _S1683 + _S1688 + _S1722.differential_0 + _S1745.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1611)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1762 = _S1691 + _S1696 + _S1728.differential_0 + _S1745.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1611)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1763 = _S1699 + _S1704 + _S1734.differential_0 + _S1745.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1611)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1764 = _S1707 + _S1712 + _S1739 + _S1739 + _S1745.vignette_params_1[int(2)].cy_0;
    (&(&_S1611)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1765 = _S1715 + _S1720 + _S1740 + _S1740 + _S1745.vignette_params_1[int(2)].cx_0;
    (&(&_S1611)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1766 = _S1685 + _S1688 + _S1724.differential_0 + _S1745.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1611)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1767 = _S1693 + _S1696 + _S1730.differential_0 + _S1745.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1611)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1768 = _S1701 + _S1704 + _S1736.differential_0 + _S1745.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1611)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1769 = _S1709 + _S1712 + _S1741 + _S1741 + _S1745.vignette_params_1[int(1)].cy_0;
    (&(&_S1611)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1770 = _S1717 + _S1720 + _S1742 + _S1742 + _S1745.vignette_params_1[int(1)].cx_0;
    (&(&_S1611)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1771 = _S1687 + _S1688 + _S1726.differential_0 + _S1745.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1611)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1772 = _S1695 + _S1696 + _S1732.differential_0 + _S1745.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1611)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1773 = _S1703 + _S1704 + _S1738.differential_0 + _S1745.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1611)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1774 = _S1711 + _S1712 + _S1743 + _S1743 + _S1745.vignette_params_1[int(0)].cy_0;
    (&(&_S1611)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1775 = _S1719 + _S1720 + _S1744 + _S1744 + _S1745.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 39>  _S1776;
    _S1776[int(0)] = 0.0f;
    _S1776[int(1)] = 0.0f;
    _S1776[int(2)] = 0.0f;
    _S1776[int(3)] = 0.0f;
    _S1776[int(4)] = 0.0f;
    _S1776[int(5)] = 0.0f;
    _S1776[int(6)] = 0.0f;
    _S1776[int(7)] = 0.0f;
    _S1776[int(8)] = 0.0f;
    _S1776[int(9)] = 0.0f;
    _S1776[int(10)] = 0.0f;
    _S1776[int(11)] = 0.0f;
    _S1776[int(12)] = 0.0f;
    _S1776[int(13)] = 0.0f;
    _S1776[int(14)] = 0.0f;
    _S1776[int(15)] = 0.0f;
    _S1776[int(16)] = 0.0f;
    _S1776[int(17)] = 0.0f;
    _S1776[int(18)] = 0.0f;
    _S1776[int(19)] = 0.0f;
    _S1776[int(20)] = 0.0f;
    _S1776[int(21)] = 0.0f;
    _S1776[int(22)] = 0.0f;
    _S1776[int(23)] = 0.0f;
    _S1776[int(24)] = 0.0f;
    _S1776[int(25)] = 0.0f;
    _S1776[int(26)] = 0.0f;
    _S1776[int(27)] = 0.0f;
    _S1776[int(28)] = 0.0f;
    _S1776[int(29)] = 0.0f;
    _S1776[int(30)] = 0.0f;
    _S1776[int(31)] = 0.0f;
    _S1776[int(32)] = 0.0f;
    _S1776[int(33)] = 0.0f;
    _S1776[int(34)] = 0.0f;
    _S1776[int(35)] = 0.0f;
    _S1776[int(36)] = 0.0f;
    _S1776[int(37)] = 0.0f;
    _S1776[int(38)] = 0.0f;
    _S1776[int(9)] = _S1767;
    _S1776[int(18)] = _S1745.color_params_1.r_0.x;
    _S1776[int(17)] = _S1745.color_params_1.b_0.y;
    _S1776[int(16)] = _S1745.color_params_1.b_0.x;
    _S1776[int(15)] = _S1761;
    _S1776[int(14)] = _S1762;
    _S1776[int(13)] = _S1763;
    _S1776[int(12)] = _S1764;
    _S1776[int(11)] = _S1765;
    _S1776[int(10)] = _S1766;
    _S1776[int(19)] = _S1745.color_params_1.r_0.y;
    _S1776[int(8)] = _S1768;
    _S1776[int(7)] = _S1769;
    _S1776[int(6)] = _S1770;
    _S1776[int(5)] = _S1771;
    _S1776[int(4)] = _S1772;
    _S1776[int(3)] = _S1773;
    _S1776[int(2)] = _S1774;
    _S1776[int(1)] = _S1775;
    _S1776[int(0)] = _S1611.exposure_1;
    _S1776[int(28)] = _S1756;
    _S1776[int(37)] = _S1747;
    _S1776[int(36)] = _S1748;
    _S1776[int(35)] = _S1749;
    _S1776[int(34)] = _S1750;
    _S1776[int(33)] = _S1751;
    _S1776[int(32)] = _S1752;
    _S1776[int(31)] = _S1753;
    _S1776[int(30)] = _S1754;
    _S1776[int(29)] = _S1755;
    _S1776[int(38)] = _S1746;
    _S1776[int(27)] = _S1757;
    _S1776[int(26)] = _S1758;
    _S1776[int(25)] = _S1759;
    _S1776[int(24)] = _S1760;
    _S1776[int(23)] = _S1745.color_params_1.n_0.y;
    _S1776[int(22)] = _S1745.color_params_1.n_0.x;
    _S1776[int(21)] = _S1745.color_params_1.g_0.y;
    _S1776[int(20)] = _S1745.color_params_1.g_0.x;
    dpparams_4->primal_0 = dpparams_4->primal_0;
    dpparams_4->differential_0 = _S1776;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S1777, FixedArray<float, 23>  * _S1778)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S1777, _S1778);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_9, FixedArray<float, 23>  grad_out_4, FixedArray<float, 39>  * _S1779)
{
    FixedArray<float, 39>  _S1780 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_4;
    (&dp_params_4)->primal_0 = params_9;
    (&dp_params_4)->differential_0 = _S1780;
    FixedArray<float, 23>  _S1781 = grad_out_4;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_4, &_S1781);
    *_S1779 = (&dp_params_4)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_no_crf_regularization_loss(FixedArray<float, 24>  params_10, FixedArray<float, 18>  * _S1782)
{
    PPISPParamsNoCRF_0 p_5;
    (&p_5)->exposure_0 = params_10[int(0)];
    (&(&p_5)->vignette_params_0[int(0)])->cx_0 = params_10[int(1)];
    (&(&p_5)->vignette_params_0[int(0)])->cy_0 = params_10[int(2)];
    (&(&p_5)->vignette_params_0[int(0)])->alpha0_0 = params_10[int(3)];
    (&(&p_5)->vignette_params_0[int(0)])->alpha1_0 = params_10[int(4)];
    (&(&p_5)->vignette_params_0[int(0)])->alpha2_0 = params_10[int(5)];
    (&(&p_5)->vignette_params_0[int(1)])->cx_0 = params_10[int(6)];
    (&(&p_5)->vignette_params_0[int(1)])->cy_0 = params_10[int(7)];
    (&(&p_5)->vignette_params_0[int(1)])->alpha0_0 = params_10[int(8)];
    (&(&p_5)->vignette_params_0[int(1)])->alpha1_0 = params_10[int(9)];
    (&(&p_5)->vignette_params_0[int(1)])->alpha2_0 = params_10[int(10)];
    (&(&p_5)->vignette_params_0[int(2)])->cx_0 = params_10[int(11)];
    (&(&p_5)->vignette_params_0[int(2)])->cy_0 = params_10[int(12)];
    (&(&p_5)->vignette_params_0[int(2)])->alpha0_0 = params_10[int(13)];
    (&(&p_5)->vignette_params_0[int(2)])->alpha1_0 = params_10[int(14)];
    (&(&p_5)->vignette_params_0[int(2)])->alpha2_0 = params_10[int(15)];
    *&((&(&(&p_5)->color_params_0)->b_0)->x) = params_10[int(16)];
    *&((&(&(&p_5)->color_params_0)->b_0)->y) = params_10[int(17)];
    *&((&(&(&p_5)->color_params_0)->r_0)->x) = params_10[int(18)];
    *&((&(&(&p_5)->color_params_0)->r_0)->y) = params_10[int(19)];
    *&((&(&(&p_5)->color_params_0)->g_0)->x) = params_10[int(20)];
    *&((&(&(&p_5)->color_params_0)->g_0)->y) = params_10[int(21)];
    *&((&(&(&p_5)->color_params_0)->n_0)->x) = params_10[int(22)];
    *&((&(&(&p_5)->color_params_0)->n_0)->y) = params_10[int(23)];
    FixedArray<float, 18>  losses_2;
    losses_2[int(0)] = 0.0f;
    losses_2[int(1)] = 0.0f;
    losses_2[int(2)] = 0.0f;
    losses_2[int(3)] = 0.0f;
    losses_2[int(4)] = 0.0f;
    losses_2[int(5)] = 0.0f;
    losses_2[int(6)] = 0.0f;
    losses_2[int(7)] = 0.0f;
    losses_2[int(8)] = 0.0f;
    losses_2[int(9)] = 0.0f;
    losses_2[int(10)] = 0.0f;
    losses_2[int(11)] = 0.0f;
    losses_2[int(12)] = 0.0f;
    losses_2[int(13)] = 0.0f;
    losses_2[int(14)] = 0.0f;
    losses_2[int(15)] = 0.0f;
    losses_2[int(16)] = 0.0f;
    losses_2[int(17)] = 0.0f;
    losses_2[int(0)] = p_5.exposure_0;
    float _S1783 = p_5.vignette_params_0[int(0)].cx_0;
    float _S1784 = p_5.vignette_params_0[int(0)].cy_0;
    float _S1785 = p_5.vignette_params_0[int(1)].cx_0;
    float _S1786 = p_5.vignette_params_0[int(1)].cy_0;
    float _S1787 = p_5.vignette_params_0[int(2)].cx_0;
    float _S1788 = p_5.vignette_params_0[int(2)].cy_0;
    losses_2[int(1)] = _S1783 * _S1783 + _S1784 * _S1784 + _S1785 * _S1785 + _S1786 * _S1786 + _S1787 * _S1787 + _S1788 * _S1788;
    losses_2[int(2)] = (F32_max((0.0f), (p_5.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(2)].alpha0_0)));
    losses_2[int(3)] = (F32_max((0.0f), (p_5.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(2)].alpha1_0)));
    losses_2[int(4)] = (F32_max((0.0f), (p_5.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_5.vignette_params_0[int(2)].alpha2_0)));
    float mean_38 = (p_5.vignette_params_0[int(0)].cx_0 + p_5.vignette_params_0[int(1)].cx_0 + p_5.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S1789 = p_5.vignette_params_0[int(0)].cx_0 - mean_38;
    float _S1790 = p_5.vignette_params_0[int(1)].cx_0 - mean_38;
    float _S1791 = p_5.vignette_params_0[int(2)].cx_0 - mean_38;
    losses_2[int(5)] = (_S1789 * _S1789 + _S1790 * _S1790 + _S1791 * _S1791) / 3.0f;
    float mean_39 = (p_5.vignette_params_0[int(0)].cy_0 + p_5.vignette_params_0[int(1)].cy_0 + p_5.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S1792 = p_5.vignette_params_0[int(0)].cy_0 - mean_39;
    float _S1793 = p_5.vignette_params_0[int(1)].cy_0 - mean_39;
    float _S1794 = p_5.vignette_params_0[int(2)].cy_0 - mean_39;
    losses_2[int(6)] = (_S1792 * _S1792 + _S1793 * _S1793 + _S1794 * _S1794) / 3.0f;
    float mean_40 = (p_5.vignette_params_0[int(0)].alpha0_0 + p_5.vignette_params_0[int(1)].alpha0_0 + p_5.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S1795 = p_5.vignette_params_0[int(0)].alpha0_0 - mean_40;
    float _S1796 = p_5.vignette_params_0[int(1)].alpha0_0 - mean_40;
    float _S1797 = p_5.vignette_params_0[int(2)].alpha0_0 - mean_40;
    losses_2[int(7)] = (_S1795 * _S1795 + _S1796 * _S1796 + _S1797 * _S1797) / 3.0f;
    float mean_41 = (p_5.vignette_params_0[int(0)].alpha1_0 + p_5.vignette_params_0[int(1)].alpha1_0 + p_5.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S1798 = p_5.vignette_params_0[int(0)].alpha1_0 - mean_41;
    float _S1799 = p_5.vignette_params_0[int(1)].alpha1_0 - mean_41;
    float _S1800 = p_5.vignette_params_0[int(2)].alpha1_0 - mean_41;
    losses_2[int(8)] = (_S1798 * _S1798 + _S1799 * _S1799 + _S1800 * _S1800) / 3.0f;
    float mean_42 = (p_5.vignette_params_0[int(0)].alpha2_0 + p_5.vignette_params_0[int(1)].alpha2_0 + p_5.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S1801 = p_5.vignette_params_0[int(0)].alpha2_0 - mean_42;
    float _S1802 = p_5.vignette_params_0[int(1)].alpha2_0 - mean_42;
    float _S1803 = p_5.vignette_params_0[int(2)].alpha2_0 - mean_42;
    losses_2[int(9)] = (_S1801 * _S1801 + _S1802 * _S1802 + _S1803 * _S1803) / 3.0f;
    float2  bd_5 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_5.color_params_0.b_0);
    float2  rd_5 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_5.color_params_0.r_0);
    float2  gd_5 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_5.color_params_0.g_0);
    float2  nd_5 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_5.color_params_0.n_0);
    losses_2[int(10)] = bd_5.x;
    losses_2[int(11)] = bd_5.y;
    losses_2[int(12)] = rd_5.x;
    losses_2[int(13)] = rd_5.y;
    losses_2[int(14)] = gd_5.x;
    losses_2[int(15)] = gd_5.y;
    losses_2[int(16)] = nd_5.x;
    losses_2[int(17)] = nd_5.y;
    *_S1782 = losses_2;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_no_crf_regularization_loss_0(DiffPair_arrayx3Cfloatx2C24x3E_0 * dpparams_5, FixedArray<float, 18>  * _s_dOut_5)
{
    VignettingChannelParams_0 _S1804 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1805 = {
        _S1804, _S1804, _S1804
    };
    float2  _S1806 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1807 = { _S1806, _S1806, _S1806, _S1806 };
    PPISPParamsNoCRF_0 _S1808;
    (&_S1808)->exposure_0 = dpparams_5->primal_0[int(0)];
    (&_S1808)->vignette_params_0 = _S1805;
    (&_S1808)->color_params_0 = _S1807;
    (&(&_S1808)->vignette_params_0[int(0)])->cx_0 = dpparams_5->primal_0[int(1)];
    (&(&_S1808)->vignette_params_0[int(0)])->cy_0 = dpparams_5->primal_0[int(2)];
    (&(&_S1808)->vignette_params_0[int(0)])->alpha0_0 = dpparams_5->primal_0[int(3)];
    (&(&_S1808)->vignette_params_0[int(0)])->alpha1_0 = dpparams_5->primal_0[int(4)];
    (&(&_S1808)->vignette_params_0[int(0)])->alpha2_0 = dpparams_5->primal_0[int(5)];
    (&(&_S1808)->vignette_params_0[int(1)])->cx_0 = dpparams_5->primal_0[int(6)];
    (&(&_S1808)->vignette_params_0[int(1)])->cy_0 = dpparams_5->primal_0[int(7)];
    (&(&_S1808)->vignette_params_0[int(1)])->alpha0_0 = dpparams_5->primal_0[int(8)];
    (&(&_S1808)->vignette_params_0[int(1)])->alpha1_0 = dpparams_5->primal_0[int(9)];
    (&(&_S1808)->vignette_params_0[int(1)])->alpha2_0 = dpparams_5->primal_0[int(10)];
    (&(&_S1808)->vignette_params_0[int(2)])->cx_0 = dpparams_5->primal_0[int(11)];
    (&(&_S1808)->vignette_params_0[int(2)])->cy_0 = dpparams_5->primal_0[int(12)];
    (&(&_S1808)->vignette_params_0[int(2)])->alpha0_0 = dpparams_5->primal_0[int(13)];
    (&(&_S1808)->vignette_params_0[int(2)])->alpha1_0 = dpparams_5->primal_0[int(14)];
    (&(&_S1808)->vignette_params_0[int(2)])->alpha2_0 = dpparams_5->primal_0[int(15)];
    *&((&(&(&_S1808)->color_params_0)->b_0)->x) = dpparams_5->primal_0[int(16)];
    *&((&(&(&_S1808)->color_params_0)->b_0)->y) = dpparams_5->primal_0[int(17)];
    *&((&(&(&_S1808)->color_params_0)->r_0)->x) = dpparams_5->primal_0[int(18)];
    *&((&(&(&_S1808)->color_params_0)->r_0)->y) = dpparams_5->primal_0[int(19)];
    *&((&(&(&_S1808)->color_params_0)->g_0)->x) = dpparams_5->primal_0[int(20)];
    *&((&(&(&_S1808)->color_params_0)->g_0)->y) = dpparams_5->primal_0[int(21)];
    *&((&(&(&_S1808)->color_params_0)->n_0)->x) = dpparams_5->primal_0[int(22)];
    *&((&(&(&_S1808)->color_params_0)->n_0)->y) = dpparams_5->primal_0[int(23)];
    float mean_43 = (dpparams_5->primal_0[int(1)] + dpparams_5->primal_0[int(6)] + dpparams_5->primal_0[int(11)]) / 3.0f;
    float _S1809 = dpparams_5->primal_0[int(1)] - mean_43;
    float _S1810 = dpparams_5->primal_0[int(6)] - mean_43;
    float _S1811 = dpparams_5->primal_0[int(11)] - mean_43;
    float mean_44 = (dpparams_5->primal_0[int(2)] + dpparams_5->primal_0[int(7)] + dpparams_5->primal_0[int(12)]) / 3.0f;
    float _S1812 = dpparams_5->primal_0[int(2)] - mean_44;
    float _S1813 = dpparams_5->primal_0[int(7)] - mean_44;
    float _S1814 = dpparams_5->primal_0[int(12)] - mean_44;
    float mean_45 = (dpparams_5->primal_0[int(3)] + dpparams_5->primal_0[int(8)] + dpparams_5->primal_0[int(13)]) / 3.0f;
    float _S1815 = dpparams_5->primal_0[int(3)] - mean_45;
    float _S1816 = dpparams_5->primal_0[int(8)] - mean_45;
    float _S1817 = dpparams_5->primal_0[int(13)] - mean_45;
    float mean_46 = (dpparams_5->primal_0[int(4)] + dpparams_5->primal_0[int(9)] + dpparams_5->primal_0[int(14)]) / 3.0f;
    float _S1818 = dpparams_5->primal_0[int(4)] - mean_46;
    float _S1819 = dpparams_5->primal_0[int(9)] - mean_46;
    float _S1820 = dpparams_5->primal_0[int(14)] - mean_46;
    float mean_47 = (dpparams_5->primal_0[int(5)] + dpparams_5->primal_0[int(10)] + dpparams_5->primal_0[int(15)]) / 3.0f;
    float _S1821 = dpparams_5->primal_0[int(5)] - mean_47;
    float _S1822 = dpparams_5->primal_0[int(10)] - mean_47;
    float _S1823 = dpparams_5->primal_0[int(15)] - mean_47;
    float2  _S1824 = make_float2 ((*_s_dOut_5)[int(16)], (*_s_dOut_5)[int(17)]);
    Matrix<float, 2, 2>  _S1825 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1826;
    (&_S1826)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1826)->differential_0 = _S1825;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1827;
    (&_S1827)->primal_0 = _S1808.color_params_0.n_0;
    (&_S1827)->differential_0 = _S1806;
    s_bwd_prop_mul_2(&_S1826, &_S1827, _S1824);
    float2  _S1828 = make_float2 ((*_s_dOut_5)[int(14)], (*_s_dOut_5)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1829;
    (&_S1829)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1829)->differential_0 = _S1825;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1830;
    (&_S1830)->primal_0 = _S1808.color_params_0.g_0;
    (&_S1830)->differential_0 = _S1806;
    s_bwd_prop_mul_2(&_S1829, &_S1830, _S1828);
    float2  _S1831 = make_float2 ((*_s_dOut_5)[int(12)], (*_s_dOut_5)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1832;
    (&_S1832)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1832)->differential_0 = _S1825;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1833;
    (&_S1833)->primal_0 = _S1808.color_params_0.r_0;
    (&_S1833)->differential_0 = _S1806;
    s_bwd_prop_mul_2(&_S1832, &_S1833, _S1831);
    float2  _S1834 = make_float2 ((*_s_dOut_5)[int(10)], (*_s_dOut_5)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1835;
    (&_S1835)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1835)->differential_0 = _S1825;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1836;
    (&_S1836)->primal_0 = _S1808.color_params_0.b_0;
    (&_S1836)->differential_0 = _S1806;
    s_bwd_prop_mul_2(&_S1835, &_S1836, _S1834);
    ColorPPISPParams_0 _S1837 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1837)->n_0 = _S1827.differential_0;
    (&_S1837)->g_0 = _S1830.differential_0;
    (&_S1837)->r_0 = _S1833.differential_0;
    (&_S1837)->b_0 = _S1836.differential_0;
    float _S1838 = 0.3333333432674408f * (*_s_dOut_5)[int(9)];
    float _S1839 = _S1823 * _S1838;
    float _S1840 = _S1839 + _S1839;
    float _S1841 = _S1822 * _S1838;
    float _S1842 = _S1841 + _S1841;
    float _S1843 = _S1821 * _S1838;
    float _S1844 = _S1843 + _S1843;
    float _S1845 = 0.3333333432674408f * (- _S1840 + - _S1842 + - _S1844);
    float _S1846 = 0.3333333432674408f * (*_s_dOut_5)[int(8)];
    float _S1847 = _S1820 * _S1846;
    float _S1848 = _S1847 + _S1847;
    float _S1849 = _S1819 * _S1846;
    float _S1850 = _S1849 + _S1849;
    float _S1851 = _S1818 * _S1846;
    float _S1852 = _S1851 + _S1851;
    float _S1853 = 0.3333333432674408f * (- _S1848 + - _S1850 + - _S1852);
    float _S1854 = 0.3333333432674408f * (*_s_dOut_5)[int(7)];
    float _S1855 = _S1817 * _S1854;
    float _S1856 = _S1855 + _S1855;
    float _S1857 = _S1816 * _S1854;
    float _S1858 = _S1857 + _S1857;
    float _S1859 = _S1815 * _S1854;
    float _S1860 = _S1859 + _S1859;
    float _S1861 = 0.3333333432674408f * (- _S1856 + - _S1858 + - _S1860);
    float _S1862 = 0.3333333432674408f * (*_s_dOut_5)[int(6)];
    float _S1863 = _S1814 * _S1862;
    float _S1864 = _S1863 + _S1863;
    float _S1865 = _S1813 * _S1862;
    float _S1866 = _S1865 + _S1865;
    float _S1867 = _S1812 * _S1862;
    float _S1868 = _S1867 + _S1867;
    float _S1869 = 0.3333333432674408f * (- _S1864 + - _S1866 + - _S1868);
    float _S1870 = 0.3333333432674408f * (*_s_dOut_5)[int(5)];
    float _S1871 = _S1811 * _S1870;
    float _S1872 = _S1871 + _S1871;
    float _S1873 = _S1810 * _S1870;
    float _S1874 = _S1873 + _S1873;
    float _S1875 = _S1809 * _S1870;
    float _S1876 = _S1875 + _S1875;
    float _S1877 = 0.3333333432674408f * (- _S1872 + - _S1874 + - _S1876);
    DiffPair_float_0 _S1878;
    (&_S1878)->primal_0 = 0.0f;
    (&_S1878)->differential_0 = 0.0f;
    DiffPair_float_0 _S1879;
    (&_S1879)->primal_0 = dpparams_5->primal_0[int(15)];
    (&_S1879)->differential_0 = 0.0f;
    _d_max_0(&_S1878, &_S1879, (*_s_dOut_5)[int(4)]);
    DiffPair_float_0 _S1880;
    (&_S1880)->primal_0 = 0.0f;
    (&_S1880)->differential_0 = 0.0f;
    DiffPair_float_0 _S1881;
    (&_S1881)->primal_0 = dpparams_5->primal_0[int(10)];
    (&_S1881)->differential_0 = 0.0f;
    _d_max_0(&_S1880, &_S1881, (*_s_dOut_5)[int(4)]);
    DiffPair_float_0 _S1882;
    (&_S1882)->primal_0 = 0.0f;
    (&_S1882)->differential_0 = 0.0f;
    DiffPair_float_0 _S1883;
    (&_S1883)->primal_0 = dpparams_5->primal_0[int(5)];
    (&_S1883)->differential_0 = 0.0f;
    _d_max_0(&_S1882, &_S1883, (*_s_dOut_5)[int(4)]);
    DiffPair_float_0 _S1884;
    (&_S1884)->primal_0 = 0.0f;
    (&_S1884)->differential_0 = 0.0f;
    DiffPair_float_0 _S1885;
    (&_S1885)->primal_0 = dpparams_5->primal_0[int(14)];
    (&_S1885)->differential_0 = 0.0f;
    _d_max_0(&_S1884, &_S1885, (*_s_dOut_5)[int(3)]);
    DiffPair_float_0 _S1886;
    (&_S1886)->primal_0 = 0.0f;
    (&_S1886)->differential_0 = 0.0f;
    DiffPair_float_0 _S1887;
    (&_S1887)->primal_0 = dpparams_5->primal_0[int(9)];
    (&_S1887)->differential_0 = 0.0f;
    _d_max_0(&_S1886, &_S1887, (*_s_dOut_5)[int(3)]);
    DiffPair_float_0 _S1888;
    (&_S1888)->primal_0 = 0.0f;
    (&_S1888)->differential_0 = 0.0f;
    DiffPair_float_0 _S1889;
    (&_S1889)->primal_0 = dpparams_5->primal_0[int(4)];
    (&_S1889)->differential_0 = 0.0f;
    _d_max_0(&_S1888, &_S1889, (*_s_dOut_5)[int(3)]);
    DiffPair_float_0 _S1890;
    (&_S1890)->primal_0 = 0.0f;
    (&_S1890)->differential_0 = 0.0f;
    DiffPair_float_0 _S1891;
    (&_S1891)->primal_0 = dpparams_5->primal_0[int(13)];
    (&_S1891)->differential_0 = 0.0f;
    _d_max_0(&_S1890, &_S1891, (*_s_dOut_5)[int(2)]);
    DiffPair_float_0 _S1892;
    (&_S1892)->primal_0 = 0.0f;
    (&_S1892)->differential_0 = 0.0f;
    DiffPair_float_0 _S1893;
    (&_S1893)->primal_0 = dpparams_5->primal_0[int(8)];
    (&_S1893)->differential_0 = 0.0f;
    _d_max_0(&_S1892, &_S1893, (*_s_dOut_5)[int(2)]);
    DiffPair_float_0 _S1894;
    (&_S1894)->primal_0 = 0.0f;
    (&_S1894)->differential_0 = 0.0f;
    DiffPair_float_0 _S1895;
    (&_S1895)->primal_0 = dpparams_5->primal_0[int(3)];
    (&_S1895)->differential_0 = 0.0f;
    _d_max_0(&_S1894, &_S1895, (*_s_dOut_5)[int(2)]);
    float _S1896 = dpparams_5->primal_0[int(12)] * (*_s_dOut_5)[int(1)];
    float _S1897 = dpparams_5->primal_0[int(11)] * (*_s_dOut_5)[int(1)];
    float _S1898 = dpparams_5->primal_0[int(7)] * (*_s_dOut_5)[int(1)];
    float _S1899 = dpparams_5->primal_0[int(6)] * (*_s_dOut_5)[int(1)];
    float _S1900 = dpparams_5->primal_0[int(2)] * (*_s_dOut_5)[int(1)];
    float _S1901 = dpparams_5->primal_0[int(1)] * (*_s_dOut_5)[int(1)];
    PPISPParamsNoCRF_0 _S1902 = PPISPParamsNoCRF_x24_syn_dzero_0();
    (&_S1902)->color_params_0 = _S1837;
    (&_S1902)->exposure_0 = (*_s_dOut_5)[int(0)];
    _S1808 = _S1902;
    *&((&(&(&_S1808)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1808)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1808)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1903 = _S1840 + _S1845 + _S1879.differential_0 + _S1902.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1808)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1904 = _S1848 + _S1853 + _S1885.differential_0 + _S1902.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1808)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1905 = _S1856 + _S1861 + _S1891.differential_0 + _S1902.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1808)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1906 = _S1864 + _S1869 + _S1896 + _S1896 + _S1902.vignette_params_0[int(2)].cy_0;
    (&(&_S1808)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1907 = _S1872 + _S1877 + _S1897 + _S1897 + _S1902.vignette_params_0[int(2)].cx_0;
    (&(&_S1808)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1908 = _S1842 + _S1845 + _S1881.differential_0 + _S1902.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1808)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1909 = _S1850 + _S1853 + _S1887.differential_0 + _S1902.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1808)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1910 = _S1858 + _S1861 + _S1893.differential_0 + _S1902.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1808)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1911 = _S1866 + _S1869 + _S1898 + _S1898 + _S1902.vignette_params_0[int(1)].cy_0;
    (&(&_S1808)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1912 = _S1874 + _S1877 + _S1899 + _S1899 + _S1902.vignette_params_0[int(1)].cx_0;
    (&(&_S1808)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1913 = _S1844 + _S1845 + _S1883.differential_0 + _S1902.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1808)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1914 = _S1852 + _S1853 + _S1889.differential_0 + _S1902.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1808)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1915 = _S1860 + _S1861 + _S1895.differential_0 + _S1902.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1808)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1916 = _S1868 + _S1869 + _S1900 + _S1900 + _S1902.vignette_params_0[int(0)].cy_0;
    (&(&_S1808)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1917 = _S1876 + _S1877 + _S1901 + _S1901 + _S1902.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 24>  _S1918;
    _S1918[int(0)] = 0.0f;
    _S1918[int(1)] = 0.0f;
    _S1918[int(2)] = 0.0f;
    _S1918[int(3)] = 0.0f;
    _S1918[int(4)] = 0.0f;
    _S1918[int(5)] = 0.0f;
    _S1918[int(6)] = 0.0f;
    _S1918[int(7)] = 0.0f;
    _S1918[int(8)] = 0.0f;
    _S1918[int(9)] = 0.0f;
    _S1918[int(10)] = 0.0f;
    _S1918[int(11)] = 0.0f;
    _S1918[int(12)] = 0.0f;
    _S1918[int(13)] = 0.0f;
    _S1918[int(14)] = 0.0f;
    _S1918[int(15)] = 0.0f;
    _S1918[int(16)] = 0.0f;
    _S1918[int(17)] = 0.0f;
    _S1918[int(18)] = 0.0f;
    _S1918[int(19)] = 0.0f;
    _S1918[int(20)] = 0.0f;
    _S1918[int(21)] = 0.0f;
    _S1918[int(22)] = 0.0f;
    _S1918[int(23)] = 0.0f;
    _S1918[int(11)] = _S1907;
    _S1918[int(0)] = _S1808.exposure_0;
    _S1918[int(1)] = _S1917;
    _S1918[int(2)] = _S1916;
    _S1918[int(3)] = _S1915;
    _S1918[int(4)] = _S1914;
    _S1918[int(5)] = _S1913;
    _S1918[int(6)] = _S1912;
    _S1918[int(7)] = _S1911;
    _S1918[int(8)] = _S1910;
    _S1918[int(9)] = _S1909;
    _S1918[int(10)] = _S1908;
    _S1918[int(23)] = _S1902.color_params_0.n_0.y;
    _S1918[int(12)] = _S1906;
    _S1918[int(13)] = _S1905;
    _S1918[int(14)] = _S1904;
    _S1918[int(15)] = _S1903;
    _S1918[int(16)] = _S1902.color_params_0.b_0.x;
    _S1918[int(17)] = _S1902.color_params_0.b_0.y;
    _S1918[int(18)] = _S1902.color_params_0.r_0.x;
    _S1918[int(19)] = _S1902.color_params_0.r_0.y;
    _S1918[int(20)] = _S1902.color_params_0.g_0.x;
    _S1918[int(21)] = _S1902.color_params_0.g_0.y;
    _S1918[int(22)] = _S1902.color_params_0.n_0.x;
    dpparams_5->primal_0 = dpparams_5->primal_0;
    dpparams_5->differential_0 = _S1918;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_no_crf_regularization_loss_0(DiffPair_arrayx3Cfloatx2C24x3E_0 * _S1919, FixedArray<float, 18>  * _S1920)
{
    s_bwd_prop_compute_raw_ppisp_no_crf_regularization_loss_0(_S1919, _S1920);
    return;
}

inline __device__ void compute_raw_ppisp_no_crf_regularization_loss_vjp(FixedArray<float, 24>  params_11, FixedArray<float, 18>  grad_out_5, FixedArray<float, 24>  * _S1921)
{
    FixedArray<float, 24>  _S1922 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C24x3E_0 dp_params_5;
    (&dp_params_5)->primal_0 = params_11;
    (&dp_params_5)->differential_0 = _S1922;
    FixedArray<float, 18>  _S1923 = grad_out_5;
    s_bwd_compute_raw_ppisp_no_crf_regularization_loss_0(&dp_params_5, &_S1923);
    *_S1921 = (&dp_params_5)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_0, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S1924)
{
    float _S1925;
    FixedArray<float, 6>  losses_3;
    float _S1926 = float(num_cameras_0);
    float _S1927 = raw_losses_0[int(0)] / _S1926;
    for(;;)
    {
        float _S1928 = (F32_abs((_S1927)));
        if(_S1928 < 0.10000000149011612f)
        {
            _S1925 = 0.5f * _S1927 * _S1927 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1925 = _S1928 - 0.05000000074505806f;
            break;
        }
    }
    losses_3[int(0)] = _S1925;
    losses_3[int(1)] = raw_losses_0[int(1)] / (3.0f * _S1926);
    losses_3[int(2)] = (raw_losses_0[int(2)] + raw_losses_0[int(3)] + raw_losses_0[int(4)]) / (9.0f * _S1926);
    losses_3[int(3)] = (raw_losses_0[int(5)] + raw_losses_0[int(6)] + raw_losses_0[int(7)] + raw_losses_0[int(8)] + raw_losses_0[int(9)]) / (5.0f * _S1926);
    float _S1929 = raw_losses_0[int(10)] / _S1926;
    for(;;)
    {
        float _S1930 = (F32_abs((_S1929)));
        if(_S1930 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1929 * _S1929 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1930 - 0.00249999994412065f;
            break;
        }
    }
    float _S1931;
    float _S1932 = raw_losses_0[int(11)] / _S1926;
    for(;;)
    {
        float _S1933 = (F32_abs((_S1932)));
        if(_S1933 < 0.00499999988824129f)
        {
            _S1931 = 0.5f * _S1932 * _S1932 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1931 = _S1933 - 0.00249999994412065f;
            break;
        }
    }
    float _S1934 = _S1925 + _S1931;
    float _S1935 = raw_losses_0[int(12)] / _S1926;
    for(;;)
    {
        float _S1936 = (F32_abs((_S1935)));
        if(_S1936 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1935 * _S1935 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1936 - 0.00249999994412065f;
            break;
        }
    }
    float _S1937 = _S1934 + _S1925;
    float _S1938 = raw_losses_0[int(13)] / _S1926;
    for(;;)
    {
        float _S1939 = (F32_abs((_S1938)));
        if(_S1939 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1938 * _S1938 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1939 - 0.00249999994412065f;
            break;
        }
    }
    float _S1940 = _S1937 + _S1925;
    float _S1941 = raw_losses_0[int(14)] / _S1926;
    for(;;)
    {
        float _S1942 = (F32_abs((_S1941)));
        if(_S1942 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1941 * _S1941 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1942 - 0.00249999994412065f;
            break;
        }
    }
    float _S1943 = _S1940 + _S1925;
    float _S1944 = raw_losses_0[int(15)] / _S1926;
    for(;;)
    {
        float _S1945 = (F32_abs((_S1944)));
        if(_S1945 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1944 * _S1944 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1945 - 0.00249999994412065f;
            break;
        }
    }
    float _S1946 = _S1943 + _S1925;
    float _S1947 = raw_losses_0[int(16)] / _S1926;
    for(;;)
    {
        float _S1948 = (F32_abs((_S1947)));
        if(_S1948 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1947 * _S1947 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1948 - 0.00249999994412065f;
            break;
        }
    }
    float _S1949 = _S1946 + _S1925;
    float _S1950 = raw_losses_0[int(17)] / _S1926;
    for(;;)
    {
        float _S1951 = (F32_abs((_S1950)));
        if(_S1951 < 0.00499999988824129f)
        {
            _S1925 = 0.5f * _S1950 * _S1950 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1925 = _S1951 - 0.00249999994412065f;
            break;
        }
    }
    float _S1952 = (_S1949 + _S1925) / 8.0f;
    float _S1953 = (raw_losses_0[int(18)] + raw_losses_0[int(19)] + raw_losses_0[int(20)] + raw_losses_0[int(21)]) / (4.0f * _S1926);
    losses_3[int(0)] = losses_3[int(0)] * loss_weights_0[int(0)];
    losses_3[int(1)] = losses_3[int(1)] * loss_weights_0[int(1)];
    losses_3[int(2)] = losses_3[int(2)] * loss_weights_0[int(2)];
    losses_3[int(3)] = losses_3[int(3)] * loss_weights_0[int(3)];
    losses_3[int(4)] = _S1952 * loss_weights_0[int(4)];
    losses_3[int(5)] = _S1953 * loss_weights_0[int(5)];
    *_S1924 = losses_3;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_1, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S1954)
{
    float _S1955;
    FixedArray<float, 6>  losses_4;
    float _S1956 = float(num_cameras_1);
    float _S1957 = raw_losses_1[int(0)] / _S1956;
    for(;;)
    {
        float _S1958 = (F32_abs((_S1957)));
        if(_S1958 < 0.10000000149011612f)
        {
            _S1955 = 0.5f * _S1957 * _S1957 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1955 = _S1958 - 0.05000000074505806f;
            break;
        }
    }
    losses_4[int(0)] = _S1955;
    losses_4[int(1)] = raw_losses_1[int(1)] / (3.0f * _S1956);
    losses_4[int(2)] = (raw_losses_1[int(2)] + raw_losses_1[int(3)] + raw_losses_1[int(4)]) / (9.0f * _S1956);
    float _S1959 = 5.0f * _S1956;
    losses_4[int(3)] = (raw_losses_1[int(5)] + raw_losses_1[int(6)] + raw_losses_1[int(7)] + raw_losses_1[int(8)] + raw_losses_1[int(9)]) / _S1959;
    float _S1960 = raw_losses_1[int(10)] / _S1956;
    for(;;)
    {
        float _S1961 = (F32_abs((_S1960)));
        if(_S1961 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1960 * _S1960 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1961 - 0.00249999994412065f;
            break;
        }
    }
    float _S1962;
    float _S1963 = raw_losses_1[int(11)] / _S1956;
    for(;;)
    {
        float _S1964 = (F32_abs((_S1963)));
        if(_S1964 < 0.00499999988824129f)
        {
            _S1962 = 0.5f * _S1963 * _S1963 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1962 = _S1964 - 0.00249999994412065f;
            break;
        }
    }
    float _S1965 = _S1955 + _S1962;
    float _S1966 = raw_losses_1[int(12)] / _S1956;
    for(;;)
    {
        float _S1967 = (F32_abs((_S1966)));
        if(_S1967 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1966 * _S1966 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1967 - 0.00249999994412065f;
            break;
        }
    }
    float _S1968 = _S1965 + _S1955;
    float _S1969 = raw_losses_1[int(13)] / _S1956;
    for(;;)
    {
        float _S1970 = (F32_abs((_S1969)));
        if(_S1970 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1969 * _S1969 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1970 - 0.00249999994412065f;
            break;
        }
    }
    float _S1971 = _S1968 + _S1955;
    float _S1972 = raw_losses_1[int(14)] / _S1956;
    for(;;)
    {
        float _S1973 = (F32_abs((_S1972)));
        if(_S1973 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1972 * _S1972 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1973 - 0.00249999994412065f;
            break;
        }
    }
    float _S1974 = _S1971 + _S1955;
    float _S1975 = raw_losses_1[int(15)] / _S1956;
    for(;;)
    {
        float _S1976 = (F32_abs((_S1975)));
        if(_S1976 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1975 * _S1975 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1976 - 0.00249999994412065f;
            break;
        }
    }
    float _S1977 = _S1974 + _S1955;
    float _S1978 = raw_losses_1[int(16)] / _S1956;
    for(;;)
    {
        float _S1979 = (F32_abs((_S1978)));
        if(_S1979 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1978 * _S1978 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1979 - 0.00249999994412065f;
            break;
        }
    }
    float _S1980 = _S1977 + _S1955;
    float _S1981 = raw_losses_1[int(17)] / _S1956;
    for(;;)
    {
        float _S1982 = (F32_abs((_S1981)));
        if(_S1982 < 0.00499999988824129f)
        {
            _S1955 = 0.5f * _S1981 * _S1981 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1955 = _S1982 - 0.00249999994412065f;
            break;
        }
    }
    float _S1983 = (_S1980 + _S1955) / 8.0f;
    float _S1984 = (raw_losses_1[int(18)] + raw_losses_1[int(19)] + raw_losses_1[int(20)] + raw_losses_1[int(21)] + raw_losses_1[int(22)]) / _S1959;
    losses_4[int(0)] = losses_4[int(0)] * loss_weights_1[int(0)];
    losses_4[int(1)] = losses_4[int(1)] * loss_weights_1[int(1)];
    losses_4[int(2)] = losses_4[int(2)] * loss_weights_1[int(2)];
    losses_4[int(3)] = losses_4[int(3)] * loss_weights_1[int(3)];
    losses_4[int(4)] = _S1983 * loss_weights_1[int(4)];
    losses_4[int(5)] = _S1984 * loss_weights_1[int(5)];
    *_S1954 = losses_4;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_0, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_6)
{
    FixedArray<float, 22>  _S1985 = dpraw_losses_0->primal_0;
    float _S1986 = float(num_cameras_2);
    float _S1987 = dpraw_losses_0->primal_0[int(0)] / _S1986;
    bool _S1988 = (s_primal_ctx_abs_0(_S1987)) < 0.10000000149011612f;
    float _S1989;
    if(_S1988)
    {
        _S1989 = 0.5f * _S1987;
    }
    else
    {
        _S1989 = 0.0f;
    }
    float _S1990 = 3.0f * _S1986;
    float _S1991 = 9.0f * _S1986;
    float _S1992 = 5.0f * _S1986;
    float _S1993 = _S1985[int(10)] / _S1986;
    bool _S1994 = (s_primal_ctx_abs_0(_S1993)) < 0.00499999988824129f;
    float _S1995;
    if(_S1994)
    {
        _S1995 = 0.5f * _S1993;
    }
    else
    {
        _S1995 = 0.0f;
    }
    float _S1996 = _S1985[int(11)] / _S1986;
    bool _S1997 = (s_primal_ctx_abs_0(_S1996)) < 0.00499999988824129f;
    float _S1998;
    if(_S1997)
    {
        _S1998 = 0.5f * _S1996;
    }
    else
    {
        _S1998 = 0.0f;
    }
    float _S1999 = _S1985[int(12)] / _S1986;
    bool _S2000 = (s_primal_ctx_abs_0(_S1999)) < 0.00499999988824129f;
    float _S2001;
    if(_S2000)
    {
        _S2001 = 0.5f * _S1999;
    }
    else
    {
        _S2001 = 0.0f;
    }
    float _S2002 = _S1985[int(13)] / _S1986;
    bool _S2003 = (s_primal_ctx_abs_0(_S2002)) < 0.00499999988824129f;
    float _S2004;
    if(_S2003)
    {
        _S2004 = 0.5f * _S2002;
    }
    else
    {
        _S2004 = 0.0f;
    }
    float _S2005 = _S1985[int(14)] / _S1986;
    bool _S2006 = (s_primal_ctx_abs_0(_S2005)) < 0.00499999988824129f;
    float _S2007;
    if(_S2006)
    {
        _S2007 = 0.5f * _S2005;
    }
    else
    {
        _S2007 = 0.0f;
    }
    float _S2008 = _S1985[int(15)] / _S1986;
    bool _S2009 = (s_primal_ctx_abs_0(_S2008)) < 0.00499999988824129f;
    float _S2010;
    if(_S2009)
    {
        _S2010 = 0.5f * _S2008;
    }
    else
    {
        _S2010 = 0.0f;
    }
    float _S2011 = _S1985[int(16)] / _S1986;
    bool _S2012 = (s_primal_ctx_abs_0(_S2011)) < 0.00499999988824129f;
    float _S2013;
    if(_S2012)
    {
        _S2013 = 0.5f * _S2011;
    }
    else
    {
        _S2013 = 0.0f;
    }
    float _S2014 = _S1985[int(17)] / _S1986;
    bool _S2015 = (s_primal_ctx_abs_0(_S2014)) < 0.00499999988824129f;
    float _S2016;
    if(_S2015)
    {
        _S2016 = 0.5f * _S2014;
    }
    else
    {
        _S2016 = 0.0f;
    }
    float _S2017 = (*loss_weights_2)[int(3)] * (*_s_dOut_6)[int(3)];
    float _S2018 = (*loss_weights_2)[int(2)] * (*_s_dOut_6)[int(2)];
    float _S2019 = (*loss_weights_2)[int(1)] * (*_s_dOut_6)[int(1)];
    float _S2020 = (*loss_weights_2)[int(0)] * (*_s_dOut_6)[int(0)];
    float _S2021 = (*loss_weights_2)[int(5)] * (*_s_dOut_6)[int(5)] / (4.0f * _S1986);
    float _S2022 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_6)[int(4)]);
    FixedArray<float, 22>  _S2023;
    _S2023[int(0)] = 0.0f;
    _S2023[int(1)] = 0.0f;
    _S2023[int(2)] = 0.0f;
    _S2023[int(3)] = 0.0f;
    _S2023[int(4)] = 0.0f;
    _S2023[int(5)] = 0.0f;
    _S2023[int(6)] = 0.0f;
    _S2023[int(7)] = 0.0f;
    _S2023[int(8)] = 0.0f;
    _S2023[int(9)] = 0.0f;
    _S2023[int(10)] = 0.0f;
    _S2023[int(11)] = 0.0f;
    _S2023[int(12)] = 0.0f;
    _S2023[int(13)] = 0.0f;
    _S2023[int(14)] = 0.0f;
    _S2023[int(15)] = 0.0f;
    _S2023[int(16)] = 0.0f;
    _S2023[int(17)] = 0.0f;
    _S2023[int(18)] = 0.0f;
    _S2023[int(19)] = 0.0f;
    _S2023[int(20)] = 0.0f;
    _S2023[int(21)] = 0.0f;
    _S2023[int(21)] = _S2021;
    _S2023[int(20)] = _S2021;
    _S2023[int(19)] = _S2021;
    _S2023[int(18)] = _S2021;
    float _S2024 = _S2023[int(0)];
    float _S2025 = _S2023[int(1)];
    float _S2026 = _S2023[int(2)];
    float _S2027 = _S2023[int(3)];
    float _S2028 = _S2023[int(4)];
    float _S2029 = _S2023[int(5)];
    float _S2030 = _S2023[int(6)];
    float _S2031 = _S2023[int(7)];
    float _S2032 = _S2023[int(8)];
    float _S2033 = _S2023[int(9)];
    float _S2034 = _S2023[int(10)];
    float _S2035 = _S2023[int(11)];
    float _S2036 = _S2023[int(12)];
    float _S2037 = _S2023[int(13)];
    float _S2038 = _S2023[int(14)];
    float _S2039 = _S2023[int(15)];
    float _S2040 = _S2023[int(16)];
    float _S2041 = _S2023[int(17)];
    float _S2042 = _S2023[int(18)];
    float _S2043 = _S2023[int(19)];
    float _S2044 = _S2023[int(20)];
    float _S2045 = _S2023[int(21)];
    float _S2046;
    if(_S2015)
    {
        float _S2047 = 200.0f * _S2022;
        float _S2048 = _S2016 * _S2047 + 0.5f * (_S2014 * _S2047);
        _S2016 = 0.0f;
        _S2046 = _S2048;
    }
    else
    {
        _S2016 = _S2022;
        _S2046 = 0.0f;
    }
    DiffPair_float_0 _S2049;
    (&_S2049)->primal_0 = _S2014;
    (&_S2049)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2049, _S2016);
    float _S2050 = (_S2049.differential_0 + _S2046) / _S1986;
    FixedArray<float, 22>  _S2051;
    _S2051[int(0)] = 0.0f;
    _S2051[int(1)] = 0.0f;
    _S2051[int(2)] = 0.0f;
    _S2051[int(3)] = 0.0f;
    _S2051[int(4)] = 0.0f;
    _S2051[int(5)] = 0.0f;
    _S2051[int(6)] = 0.0f;
    _S2051[int(7)] = 0.0f;
    _S2051[int(8)] = 0.0f;
    _S2051[int(9)] = 0.0f;
    _S2051[int(10)] = 0.0f;
    _S2051[int(11)] = 0.0f;
    _S2051[int(12)] = 0.0f;
    _S2051[int(13)] = 0.0f;
    _S2051[int(14)] = 0.0f;
    _S2051[int(15)] = 0.0f;
    _S2051[int(16)] = 0.0f;
    _S2051[int(17)] = 0.0f;
    _S2051[int(18)] = 0.0f;
    _S2051[int(19)] = 0.0f;
    _S2051[int(20)] = 0.0f;
    _S2051[int(21)] = 0.0f;
    _S2051[int(17)] = _S2050;
    float _S2052 = _S2024 + _S2051[int(0)];
    float _S2053 = _S2025 + _S2051[int(1)];
    float _S2054 = _S2026 + _S2051[int(2)];
    float _S2055 = _S2027 + _S2051[int(3)];
    float _S2056 = _S2028 + _S2051[int(4)];
    float _S2057 = _S2029 + _S2051[int(5)];
    float _S2058 = _S2030 + _S2051[int(6)];
    float _S2059 = _S2031 + _S2051[int(7)];
    float _S2060 = _S2032 + _S2051[int(8)];
    float _S2061 = _S2033 + _S2051[int(9)];
    float _S2062 = _S2034 + _S2051[int(10)];
    float _S2063 = _S2035 + _S2051[int(11)];
    float _S2064 = _S2036 + _S2051[int(12)];
    float _S2065 = _S2037 + _S2051[int(13)];
    float _S2066 = _S2038 + _S2051[int(14)];
    float _S2067 = _S2039 + _S2051[int(15)];
    float _S2068 = _S2040 + _S2051[int(16)];
    float _S2069 = _S2041 + _S2051[int(17)];
    float _S2070 = _S2042 + _S2051[int(18)];
    float _S2071 = _S2043 + _S2051[int(19)];
    float _S2072 = _S2044 + _S2051[int(20)];
    float _S2073 = _S2045 + _S2051[int(21)];
    if(_S2012)
    {
        float _S2074 = 200.0f * _S2022;
        float _S2075 = _S2013 * _S2074 + 0.5f * (_S2011 * _S2074);
        _S2013 = 0.0f;
        _S2016 = _S2075;
    }
    else
    {
        _S2013 = _S2022;
        _S2016 = 0.0f;
    }
    DiffPair_float_0 _S2076;
    (&_S2076)->primal_0 = _S2011;
    (&_S2076)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2076, _S2013);
    float _S2077 = (_S2076.differential_0 + _S2016) / _S1986;
    FixedArray<float, 22>  _S2078;
    _S2078[int(0)] = 0.0f;
    _S2078[int(1)] = 0.0f;
    _S2078[int(2)] = 0.0f;
    _S2078[int(3)] = 0.0f;
    _S2078[int(4)] = 0.0f;
    _S2078[int(5)] = 0.0f;
    _S2078[int(6)] = 0.0f;
    _S2078[int(7)] = 0.0f;
    _S2078[int(8)] = 0.0f;
    _S2078[int(9)] = 0.0f;
    _S2078[int(10)] = 0.0f;
    _S2078[int(11)] = 0.0f;
    _S2078[int(12)] = 0.0f;
    _S2078[int(13)] = 0.0f;
    _S2078[int(14)] = 0.0f;
    _S2078[int(15)] = 0.0f;
    _S2078[int(16)] = 0.0f;
    _S2078[int(17)] = 0.0f;
    _S2078[int(18)] = 0.0f;
    _S2078[int(19)] = 0.0f;
    _S2078[int(20)] = 0.0f;
    _S2078[int(21)] = 0.0f;
    _S2078[int(16)] = _S2077;
    float _S2079 = _S2052 + _S2078[int(0)];
    float _S2080 = _S2053 + _S2078[int(1)];
    float _S2081 = _S2054 + _S2078[int(2)];
    float _S2082 = _S2055 + _S2078[int(3)];
    float _S2083 = _S2056 + _S2078[int(4)];
    float _S2084 = _S2057 + _S2078[int(5)];
    float _S2085 = _S2058 + _S2078[int(6)];
    float _S2086 = _S2059 + _S2078[int(7)];
    float _S2087 = _S2060 + _S2078[int(8)];
    float _S2088 = _S2061 + _S2078[int(9)];
    float _S2089 = _S2062 + _S2078[int(10)];
    float _S2090 = _S2063 + _S2078[int(11)];
    float _S2091 = _S2064 + _S2078[int(12)];
    float _S2092 = _S2065 + _S2078[int(13)];
    float _S2093 = _S2066 + _S2078[int(14)];
    float _S2094 = _S2067 + _S2078[int(15)];
    float _S2095 = _S2068 + _S2078[int(16)];
    float _S2096 = _S2069 + _S2078[int(17)];
    float _S2097 = _S2070 + _S2078[int(18)];
    float _S2098 = _S2071 + _S2078[int(19)];
    float _S2099 = _S2072 + _S2078[int(20)];
    float _S2100 = _S2073 + _S2078[int(21)];
    if(_S2009)
    {
        float _S2101 = 200.0f * _S2022;
        float _S2102 = _S2010 * _S2101 + 0.5f * (_S2008 * _S2101);
        _S2010 = 0.0f;
        _S2013 = _S2102;
    }
    else
    {
        _S2010 = _S2022;
        _S2013 = 0.0f;
    }
    DiffPair_float_0 _S2103;
    (&_S2103)->primal_0 = _S2008;
    (&_S2103)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2103, _S2010);
    float _S2104 = (_S2103.differential_0 + _S2013) / _S1986;
    FixedArray<float, 22>  _S2105;
    _S2105[int(0)] = 0.0f;
    _S2105[int(1)] = 0.0f;
    _S2105[int(2)] = 0.0f;
    _S2105[int(3)] = 0.0f;
    _S2105[int(4)] = 0.0f;
    _S2105[int(5)] = 0.0f;
    _S2105[int(6)] = 0.0f;
    _S2105[int(7)] = 0.0f;
    _S2105[int(8)] = 0.0f;
    _S2105[int(9)] = 0.0f;
    _S2105[int(10)] = 0.0f;
    _S2105[int(11)] = 0.0f;
    _S2105[int(12)] = 0.0f;
    _S2105[int(13)] = 0.0f;
    _S2105[int(14)] = 0.0f;
    _S2105[int(15)] = 0.0f;
    _S2105[int(16)] = 0.0f;
    _S2105[int(17)] = 0.0f;
    _S2105[int(18)] = 0.0f;
    _S2105[int(19)] = 0.0f;
    _S2105[int(20)] = 0.0f;
    _S2105[int(21)] = 0.0f;
    _S2105[int(15)] = _S2104;
    float _S2106 = _S2079 + _S2105[int(0)];
    float _S2107 = _S2080 + _S2105[int(1)];
    float _S2108 = _S2081 + _S2105[int(2)];
    float _S2109 = _S2082 + _S2105[int(3)];
    float _S2110 = _S2083 + _S2105[int(4)];
    float _S2111 = _S2084 + _S2105[int(5)];
    float _S2112 = _S2085 + _S2105[int(6)];
    float _S2113 = _S2086 + _S2105[int(7)];
    float _S2114 = _S2087 + _S2105[int(8)];
    float _S2115 = _S2088 + _S2105[int(9)];
    float _S2116 = _S2089 + _S2105[int(10)];
    float _S2117 = _S2090 + _S2105[int(11)];
    float _S2118 = _S2091 + _S2105[int(12)];
    float _S2119 = _S2092 + _S2105[int(13)];
    float _S2120 = _S2093 + _S2105[int(14)];
    float _S2121 = _S2094 + _S2105[int(15)];
    float _S2122 = _S2095 + _S2105[int(16)];
    float _S2123 = _S2096 + _S2105[int(17)];
    float _S2124 = _S2097 + _S2105[int(18)];
    float _S2125 = _S2098 + _S2105[int(19)];
    float _S2126 = _S2099 + _S2105[int(20)];
    float _S2127 = _S2100 + _S2105[int(21)];
    if(_S2006)
    {
        float _S2128 = 200.0f * _S2022;
        float _S2129 = _S2007 * _S2128 + 0.5f * (_S2005 * _S2128);
        _S2007 = 0.0f;
        _S2010 = _S2129;
    }
    else
    {
        _S2007 = _S2022;
        _S2010 = 0.0f;
    }
    DiffPair_float_0 _S2130;
    (&_S2130)->primal_0 = _S2005;
    (&_S2130)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2130, _S2007);
    float _S2131 = (_S2130.differential_0 + _S2010) / _S1986;
    FixedArray<float, 22>  _S2132;
    _S2132[int(0)] = 0.0f;
    _S2132[int(1)] = 0.0f;
    _S2132[int(2)] = 0.0f;
    _S2132[int(3)] = 0.0f;
    _S2132[int(4)] = 0.0f;
    _S2132[int(5)] = 0.0f;
    _S2132[int(6)] = 0.0f;
    _S2132[int(7)] = 0.0f;
    _S2132[int(8)] = 0.0f;
    _S2132[int(9)] = 0.0f;
    _S2132[int(10)] = 0.0f;
    _S2132[int(11)] = 0.0f;
    _S2132[int(12)] = 0.0f;
    _S2132[int(13)] = 0.0f;
    _S2132[int(14)] = 0.0f;
    _S2132[int(15)] = 0.0f;
    _S2132[int(16)] = 0.0f;
    _S2132[int(17)] = 0.0f;
    _S2132[int(18)] = 0.0f;
    _S2132[int(19)] = 0.0f;
    _S2132[int(20)] = 0.0f;
    _S2132[int(21)] = 0.0f;
    _S2132[int(14)] = _S2131;
    float _S2133 = _S2106 + _S2132[int(0)];
    float _S2134 = _S2107 + _S2132[int(1)];
    float _S2135 = _S2108 + _S2132[int(2)];
    float _S2136 = _S2109 + _S2132[int(3)];
    float _S2137 = _S2110 + _S2132[int(4)];
    float _S2138 = _S2111 + _S2132[int(5)];
    float _S2139 = _S2112 + _S2132[int(6)];
    float _S2140 = _S2113 + _S2132[int(7)];
    float _S2141 = _S2114 + _S2132[int(8)];
    float _S2142 = _S2115 + _S2132[int(9)];
    float _S2143 = _S2116 + _S2132[int(10)];
    float _S2144 = _S2117 + _S2132[int(11)];
    float _S2145 = _S2118 + _S2132[int(12)];
    float _S2146 = _S2119 + _S2132[int(13)];
    float _S2147 = _S2120 + _S2132[int(14)];
    float _S2148 = _S2121 + _S2132[int(15)];
    float _S2149 = _S2122 + _S2132[int(16)];
    float _S2150 = _S2123 + _S2132[int(17)];
    float _S2151 = _S2124 + _S2132[int(18)];
    float _S2152 = _S2125 + _S2132[int(19)];
    float _S2153 = _S2126 + _S2132[int(20)];
    float _S2154 = _S2127 + _S2132[int(21)];
    if(_S2003)
    {
        float _S2155 = 200.0f * _S2022;
        float _S2156 = _S2004 * _S2155 + 0.5f * (_S2002 * _S2155);
        _S2004 = 0.0f;
        _S2007 = _S2156;
    }
    else
    {
        _S2004 = _S2022;
        _S2007 = 0.0f;
    }
    DiffPair_float_0 _S2157;
    (&_S2157)->primal_0 = _S2002;
    (&_S2157)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2157, _S2004);
    float _S2158 = (_S2157.differential_0 + _S2007) / _S1986;
    FixedArray<float, 22>  _S2159;
    _S2159[int(0)] = 0.0f;
    _S2159[int(1)] = 0.0f;
    _S2159[int(2)] = 0.0f;
    _S2159[int(3)] = 0.0f;
    _S2159[int(4)] = 0.0f;
    _S2159[int(5)] = 0.0f;
    _S2159[int(6)] = 0.0f;
    _S2159[int(7)] = 0.0f;
    _S2159[int(8)] = 0.0f;
    _S2159[int(9)] = 0.0f;
    _S2159[int(10)] = 0.0f;
    _S2159[int(11)] = 0.0f;
    _S2159[int(12)] = 0.0f;
    _S2159[int(13)] = 0.0f;
    _S2159[int(14)] = 0.0f;
    _S2159[int(15)] = 0.0f;
    _S2159[int(16)] = 0.0f;
    _S2159[int(17)] = 0.0f;
    _S2159[int(18)] = 0.0f;
    _S2159[int(19)] = 0.0f;
    _S2159[int(20)] = 0.0f;
    _S2159[int(21)] = 0.0f;
    _S2159[int(13)] = _S2158;
    float _S2160 = _S2133 + _S2159[int(0)];
    float _S2161 = _S2134 + _S2159[int(1)];
    float _S2162 = _S2135 + _S2159[int(2)];
    float _S2163 = _S2136 + _S2159[int(3)];
    float _S2164 = _S2137 + _S2159[int(4)];
    float _S2165 = _S2138 + _S2159[int(5)];
    float _S2166 = _S2139 + _S2159[int(6)];
    float _S2167 = _S2140 + _S2159[int(7)];
    float _S2168 = _S2141 + _S2159[int(8)];
    float _S2169 = _S2142 + _S2159[int(9)];
    float _S2170 = _S2143 + _S2159[int(10)];
    float _S2171 = _S2144 + _S2159[int(11)];
    float _S2172 = _S2145 + _S2159[int(12)];
    float _S2173 = _S2146 + _S2159[int(13)];
    float _S2174 = _S2147 + _S2159[int(14)];
    float _S2175 = _S2148 + _S2159[int(15)];
    float _S2176 = _S2149 + _S2159[int(16)];
    float _S2177 = _S2150 + _S2159[int(17)];
    float _S2178 = _S2151 + _S2159[int(18)];
    float _S2179 = _S2152 + _S2159[int(19)];
    float _S2180 = _S2153 + _S2159[int(20)];
    float _S2181 = _S2154 + _S2159[int(21)];
    if(_S2000)
    {
        float _S2182 = 200.0f * _S2022;
        float _S2183 = _S2001 * _S2182 + 0.5f * (_S1999 * _S2182);
        _S2001 = 0.0f;
        _S2004 = _S2183;
    }
    else
    {
        _S2001 = _S2022;
        _S2004 = 0.0f;
    }
    DiffPair_float_0 _S2184;
    (&_S2184)->primal_0 = _S1999;
    (&_S2184)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2184, _S2001);
    float _S2185 = (_S2184.differential_0 + _S2004) / _S1986;
    FixedArray<float, 22>  _S2186;
    _S2186[int(0)] = 0.0f;
    _S2186[int(1)] = 0.0f;
    _S2186[int(2)] = 0.0f;
    _S2186[int(3)] = 0.0f;
    _S2186[int(4)] = 0.0f;
    _S2186[int(5)] = 0.0f;
    _S2186[int(6)] = 0.0f;
    _S2186[int(7)] = 0.0f;
    _S2186[int(8)] = 0.0f;
    _S2186[int(9)] = 0.0f;
    _S2186[int(10)] = 0.0f;
    _S2186[int(11)] = 0.0f;
    _S2186[int(12)] = 0.0f;
    _S2186[int(13)] = 0.0f;
    _S2186[int(14)] = 0.0f;
    _S2186[int(15)] = 0.0f;
    _S2186[int(16)] = 0.0f;
    _S2186[int(17)] = 0.0f;
    _S2186[int(18)] = 0.0f;
    _S2186[int(19)] = 0.0f;
    _S2186[int(20)] = 0.0f;
    _S2186[int(21)] = 0.0f;
    _S2186[int(12)] = _S2185;
    float _S2187 = _S2160 + _S2186[int(0)];
    float _S2188 = _S2161 + _S2186[int(1)];
    float _S2189 = _S2162 + _S2186[int(2)];
    float _S2190 = _S2163 + _S2186[int(3)];
    float _S2191 = _S2164 + _S2186[int(4)];
    float _S2192 = _S2165 + _S2186[int(5)];
    float _S2193 = _S2166 + _S2186[int(6)];
    float _S2194 = _S2167 + _S2186[int(7)];
    float _S2195 = _S2168 + _S2186[int(8)];
    float _S2196 = _S2169 + _S2186[int(9)];
    float _S2197 = _S2170 + _S2186[int(10)];
    float _S2198 = _S2171 + _S2186[int(11)];
    float _S2199 = _S2172 + _S2186[int(12)];
    float _S2200 = _S2173 + _S2186[int(13)];
    float _S2201 = _S2174 + _S2186[int(14)];
    float _S2202 = _S2175 + _S2186[int(15)];
    float _S2203 = _S2176 + _S2186[int(16)];
    float _S2204 = _S2177 + _S2186[int(17)];
    float _S2205 = _S2178 + _S2186[int(18)];
    float _S2206 = _S2179 + _S2186[int(19)];
    float _S2207 = _S2180 + _S2186[int(20)];
    float _S2208 = _S2181 + _S2186[int(21)];
    if(_S1997)
    {
        float _S2209 = 200.0f * _S2022;
        float _S2210 = _S1998 * _S2209 + 0.5f * (_S1996 * _S2209);
        _S1998 = 0.0f;
        _S2001 = _S2210;
    }
    else
    {
        _S1998 = _S2022;
        _S2001 = 0.0f;
    }
    DiffPair_float_0 _S2211;
    (&_S2211)->primal_0 = _S1996;
    (&_S2211)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2211, _S1998);
    float _S2212 = (_S2211.differential_0 + _S2001) / _S1986;
    FixedArray<float, 22>  _S2213;
    _S2213[int(0)] = 0.0f;
    _S2213[int(1)] = 0.0f;
    _S2213[int(2)] = 0.0f;
    _S2213[int(3)] = 0.0f;
    _S2213[int(4)] = 0.0f;
    _S2213[int(5)] = 0.0f;
    _S2213[int(6)] = 0.0f;
    _S2213[int(7)] = 0.0f;
    _S2213[int(8)] = 0.0f;
    _S2213[int(9)] = 0.0f;
    _S2213[int(10)] = 0.0f;
    _S2213[int(11)] = 0.0f;
    _S2213[int(12)] = 0.0f;
    _S2213[int(13)] = 0.0f;
    _S2213[int(14)] = 0.0f;
    _S2213[int(15)] = 0.0f;
    _S2213[int(16)] = 0.0f;
    _S2213[int(17)] = 0.0f;
    _S2213[int(18)] = 0.0f;
    _S2213[int(19)] = 0.0f;
    _S2213[int(20)] = 0.0f;
    _S2213[int(21)] = 0.0f;
    _S2213[int(11)] = _S2212;
    float _S2214 = _S2187 + _S2213[int(0)];
    float _S2215 = _S2188 + _S2213[int(1)];
    float _S2216 = _S2189 + _S2213[int(2)];
    float _S2217 = _S2190 + _S2213[int(3)];
    float _S2218 = _S2191 + _S2213[int(4)];
    float _S2219 = _S2192 + _S2213[int(5)];
    float _S2220 = _S2193 + _S2213[int(6)];
    float _S2221 = _S2194 + _S2213[int(7)];
    float _S2222 = _S2195 + _S2213[int(8)];
    float _S2223 = _S2196 + _S2213[int(9)];
    float _S2224 = _S2197 + _S2213[int(10)];
    float _S2225 = _S2198 + _S2213[int(11)];
    float _S2226 = _S2199 + _S2213[int(12)];
    float _S2227 = _S2200 + _S2213[int(13)];
    float _S2228 = _S2201 + _S2213[int(14)];
    float _S2229 = _S2202 + _S2213[int(15)];
    float _S2230 = _S2203 + _S2213[int(16)];
    float _S2231 = _S2204 + _S2213[int(17)];
    float _S2232 = _S2205 + _S2213[int(18)];
    float _S2233 = _S2206 + _S2213[int(19)];
    float _S2234 = _S2207 + _S2213[int(20)];
    float _S2235 = _S2208 + _S2213[int(21)];
    if(_S1994)
    {
        float _S2236 = 200.0f * _S2022;
        float _S2237 = _S1995 * _S2236 + 0.5f * (_S1993 * _S2236);
        _S1995 = 0.0f;
        _S1998 = _S2237;
    }
    else
    {
        _S1995 = _S2022;
        _S1998 = 0.0f;
    }
    DiffPair_float_0 _S2238;
    (&_S2238)->primal_0 = _S1993;
    (&_S2238)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2238, _S1995);
    float _S2239 = (_S2238.differential_0 + _S1998) / _S1986;
    float _S2240 = _S2017 / _S1992;
    float _S2241 = _S2018 / _S1991;
    float _S2242 = _S2019 / _S1990;
    FixedArray<float, 22>  _S2243;
    _S2243[int(0)] = 0.0f;
    _S2243[int(1)] = 0.0f;
    _S2243[int(2)] = 0.0f;
    _S2243[int(3)] = 0.0f;
    _S2243[int(4)] = 0.0f;
    _S2243[int(5)] = 0.0f;
    _S2243[int(6)] = 0.0f;
    _S2243[int(7)] = 0.0f;
    _S2243[int(8)] = 0.0f;
    _S2243[int(9)] = 0.0f;
    _S2243[int(10)] = 0.0f;
    _S2243[int(11)] = 0.0f;
    _S2243[int(12)] = 0.0f;
    _S2243[int(13)] = 0.0f;
    _S2243[int(14)] = 0.0f;
    _S2243[int(15)] = 0.0f;
    _S2243[int(16)] = 0.0f;
    _S2243[int(17)] = 0.0f;
    _S2243[int(18)] = 0.0f;
    _S2243[int(19)] = 0.0f;
    _S2243[int(20)] = 0.0f;
    _S2243[int(21)] = 0.0f;
    _S2243[int(10)] = _S2239;
    _S2243[int(9)] = _S2240;
    _S2243[int(8)] = _S2240;
    _S2243[int(7)] = _S2240;
    _S2243[int(6)] = _S2240;
    _S2243[int(5)] = _S2240;
    _S2243[int(4)] = _S2241;
    _S2243[int(3)] = _S2241;
    _S2243[int(2)] = _S2241;
    _S2243[int(1)] = _S2242;
    float _S2244 = _S2214 + _S2243[int(0)];
    float _S2245 = _S2215 + _S2243[int(1)];
    float _S2246 = _S2216 + _S2243[int(2)];
    float _S2247 = _S2217 + _S2243[int(3)];
    float _S2248 = _S2218 + _S2243[int(4)];
    float _S2249 = _S2219 + _S2243[int(5)];
    float _S2250 = _S2220 + _S2243[int(6)];
    float _S2251 = _S2221 + _S2243[int(7)];
    float _S2252 = _S2222 + _S2243[int(8)];
    float _S2253 = _S2223 + _S2243[int(9)];
    float _S2254 = _S2224 + _S2243[int(10)];
    float _S2255 = _S2225 + _S2243[int(11)];
    float _S2256 = _S2226 + _S2243[int(12)];
    float _S2257 = _S2227 + _S2243[int(13)];
    float _S2258 = _S2228 + _S2243[int(14)];
    float _S2259 = _S2229 + _S2243[int(15)];
    float _S2260 = _S2230 + _S2243[int(16)];
    float _S2261 = _S2231 + _S2243[int(17)];
    float _S2262 = _S2232 + _S2243[int(18)];
    float _S2263 = _S2233 + _S2243[int(19)];
    float _S2264 = _S2234 + _S2243[int(20)];
    float _S2265 = _S2235 + _S2243[int(21)];
    if(_S1988)
    {
        float _S2266 = 10.0f * _S2020;
        float _S2267 = _S1989 * _S2266 + 0.5f * (_S1987 * _S2266);
        _S1989 = 0.0f;
        _S1995 = _S2267;
    }
    else
    {
        _S1989 = _S2020;
        _S1995 = 0.0f;
    }
    DiffPair_float_0 _S2268;
    (&_S2268)->primal_0 = _S1987;
    (&_S2268)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2268, _S1989);
    float _S2269 = (_S2268.differential_0 + _S1995) / _S1986;
    FixedArray<float, 22>  _S2270;
    _S2270[int(0)] = 0.0f;
    _S2270[int(1)] = 0.0f;
    _S2270[int(2)] = 0.0f;
    _S2270[int(3)] = 0.0f;
    _S2270[int(4)] = 0.0f;
    _S2270[int(5)] = 0.0f;
    _S2270[int(6)] = 0.0f;
    _S2270[int(7)] = 0.0f;
    _S2270[int(8)] = 0.0f;
    _S2270[int(9)] = 0.0f;
    _S2270[int(10)] = 0.0f;
    _S2270[int(11)] = 0.0f;
    _S2270[int(12)] = 0.0f;
    _S2270[int(13)] = 0.0f;
    _S2270[int(14)] = 0.0f;
    _S2270[int(15)] = 0.0f;
    _S2270[int(16)] = 0.0f;
    _S2270[int(17)] = 0.0f;
    _S2270[int(18)] = 0.0f;
    _S2270[int(19)] = 0.0f;
    _S2270[int(20)] = 0.0f;
    _S2270[int(21)] = 0.0f;
    _S2270[int(0)] = _S2269;
    FixedArray<float, 22>  _S2271 = {
        _S2244 + _S2270[int(0)], _S2245 + _S2270[int(1)], _S2246 + _S2270[int(2)], _S2247 + _S2270[int(3)], _S2248 + _S2270[int(4)], _S2249 + _S2270[int(5)], _S2250 + _S2270[int(6)], _S2251 + _S2270[int(7)], _S2252 + _S2270[int(8)], _S2253 + _S2270[int(9)], _S2254 + _S2270[int(10)], _S2255 + _S2270[int(11)], _S2256 + _S2270[int(12)], _S2257 + _S2270[int(13)], _S2258 + _S2270[int(14)], _S2259 + _S2270[int(15)], _S2260 + _S2270[int(16)], _S2261 + _S2270[int(17)], _S2262 + _S2270[int(18)], _S2263 + _S2270[int(19)], _S2264 + _S2270[int(20)], _S2265 + _S2270[int(21)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S2271;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S2272, int _S2273, FixedArray<float, 6>  * _S2274, FixedArray<float, 6>  * _S2275)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S2272, _S2273, _S2274, _S2275);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_2, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_6, FixedArray<float, 22>  * _S2276)
{
    FixedArray<float, 22>  _S2277 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_2;
    (&dp_raw_losses_0)->differential_0 = _S2277;
    FixedArray<float, 6>  _S2278 = loss_weights_3;
    FixedArray<float, 6>  _S2279 = grad_out_6;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_0, num_cameras_3, &_S2278, &_S2279);
    *_S2276 = (&dp_raw_losses_0)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_1, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_7)
{
    FixedArray<float, 23>  _S2280 = dpraw_losses_1->primal_0;
    float _S2281 = float(num_cameras_4);
    float _S2282 = dpraw_losses_1->primal_0[int(0)] / _S2281;
    bool _S2283 = (s_primal_ctx_abs_0(_S2282)) < 0.10000000149011612f;
    float _S2284;
    if(_S2283)
    {
        _S2284 = 0.5f * _S2282;
    }
    else
    {
        _S2284 = 0.0f;
    }
    float _S2285 = 3.0f * _S2281;
    float _S2286 = 9.0f * _S2281;
    float _S2287 = 5.0f * _S2281;
    float _S2288 = _S2280[int(10)] / _S2281;
    bool _S2289 = (s_primal_ctx_abs_0(_S2288)) < 0.00499999988824129f;
    float _S2290;
    if(_S2289)
    {
        _S2290 = 0.5f * _S2288;
    }
    else
    {
        _S2290 = 0.0f;
    }
    float _S2291 = _S2280[int(11)] / _S2281;
    bool _S2292 = (s_primal_ctx_abs_0(_S2291)) < 0.00499999988824129f;
    float _S2293;
    if(_S2292)
    {
        _S2293 = 0.5f * _S2291;
    }
    else
    {
        _S2293 = 0.0f;
    }
    float _S2294 = _S2280[int(12)] / _S2281;
    bool _S2295 = (s_primal_ctx_abs_0(_S2294)) < 0.00499999988824129f;
    float _S2296;
    if(_S2295)
    {
        _S2296 = 0.5f * _S2294;
    }
    else
    {
        _S2296 = 0.0f;
    }
    float _S2297 = _S2280[int(13)] / _S2281;
    bool _S2298 = (s_primal_ctx_abs_0(_S2297)) < 0.00499999988824129f;
    float _S2299;
    if(_S2298)
    {
        _S2299 = 0.5f * _S2297;
    }
    else
    {
        _S2299 = 0.0f;
    }
    float _S2300 = _S2280[int(14)] / _S2281;
    bool _S2301 = (s_primal_ctx_abs_0(_S2300)) < 0.00499999988824129f;
    float _S2302;
    if(_S2301)
    {
        _S2302 = 0.5f * _S2300;
    }
    else
    {
        _S2302 = 0.0f;
    }
    float _S2303 = _S2280[int(15)] / _S2281;
    bool _S2304 = (s_primal_ctx_abs_0(_S2303)) < 0.00499999988824129f;
    float _S2305;
    if(_S2304)
    {
        _S2305 = 0.5f * _S2303;
    }
    else
    {
        _S2305 = 0.0f;
    }
    float _S2306 = _S2280[int(16)] / _S2281;
    bool _S2307 = (s_primal_ctx_abs_0(_S2306)) < 0.00499999988824129f;
    float _S2308;
    if(_S2307)
    {
        _S2308 = 0.5f * _S2306;
    }
    else
    {
        _S2308 = 0.0f;
    }
    float _S2309 = _S2280[int(17)] / _S2281;
    bool _S2310 = (s_primal_ctx_abs_0(_S2309)) < 0.00499999988824129f;
    float _S2311;
    if(_S2310)
    {
        _S2311 = 0.5f * _S2309;
    }
    else
    {
        _S2311 = 0.0f;
    }
    float _S2312 = (*loss_weights_4)[int(3)] * (*_s_dOut_7)[int(3)];
    float _S2313 = (*loss_weights_4)[int(2)] * (*_s_dOut_7)[int(2)];
    float _S2314 = (*loss_weights_4)[int(1)] * (*_s_dOut_7)[int(1)];
    float _S2315 = (*loss_weights_4)[int(0)] * (*_s_dOut_7)[int(0)];
    float _S2316 = (*loss_weights_4)[int(5)] * (*_s_dOut_7)[int(5)] / _S2287;
    float _S2317 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_7)[int(4)]);
    FixedArray<float, 23>  _S2318;
    _S2318[int(0)] = 0.0f;
    _S2318[int(1)] = 0.0f;
    _S2318[int(2)] = 0.0f;
    _S2318[int(3)] = 0.0f;
    _S2318[int(4)] = 0.0f;
    _S2318[int(5)] = 0.0f;
    _S2318[int(6)] = 0.0f;
    _S2318[int(7)] = 0.0f;
    _S2318[int(8)] = 0.0f;
    _S2318[int(9)] = 0.0f;
    _S2318[int(10)] = 0.0f;
    _S2318[int(11)] = 0.0f;
    _S2318[int(12)] = 0.0f;
    _S2318[int(13)] = 0.0f;
    _S2318[int(14)] = 0.0f;
    _S2318[int(15)] = 0.0f;
    _S2318[int(16)] = 0.0f;
    _S2318[int(17)] = 0.0f;
    _S2318[int(18)] = 0.0f;
    _S2318[int(19)] = 0.0f;
    _S2318[int(20)] = 0.0f;
    _S2318[int(21)] = 0.0f;
    _S2318[int(22)] = 0.0f;
    _S2318[int(22)] = _S2316;
    _S2318[int(21)] = _S2316;
    _S2318[int(20)] = _S2316;
    _S2318[int(19)] = _S2316;
    _S2318[int(18)] = _S2316;
    float _S2319 = _S2318[int(0)];
    float _S2320 = _S2318[int(1)];
    float _S2321 = _S2318[int(2)];
    float _S2322 = _S2318[int(3)];
    float _S2323 = _S2318[int(4)];
    float _S2324 = _S2318[int(5)];
    float _S2325 = _S2318[int(6)];
    float _S2326 = _S2318[int(7)];
    float _S2327 = _S2318[int(8)];
    float _S2328 = _S2318[int(9)];
    float _S2329 = _S2318[int(10)];
    float _S2330 = _S2318[int(11)];
    float _S2331 = _S2318[int(12)];
    float _S2332 = _S2318[int(13)];
    float _S2333 = _S2318[int(14)];
    float _S2334 = _S2318[int(15)];
    float _S2335 = _S2318[int(16)];
    float _S2336 = _S2318[int(17)];
    float _S2337 = _S2318[int(18)];
    float _S2338 = _S2318[int(19)];
    float _S2339 = _S2318[int(20)];
    float _S2340 = _S2318[int(21)];
    float _S2341 = _S2318[int(22)];
    float _S2342;
    if(_S2310)
    {
        float _S2343 = 200.0f * _S2317;
        float _S2344 = _S2311 * _S2343 + 0.5f * (_S2309 * _S2343);
        _S2311 = 0.0f;
        _S2342 = _S2344;
    }
    else
    {
        _S2311 = _S2317;
        _S2342 = 0.0f;
    }
    DiffPair_float_0 _S2345;
    (&_S2345)->primal_0 = _S2309;
    (&_S2345)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2345, _S2311);
    float _S2346 = (_S2345.differential_0 + _S2342) / _S2281;
    FixedArray<float, 23>  _S2347;
    _S2347[int(0)] = 0.0f;
    _S2347[int(1)] = 0.0f;
    _S2347[int(2)] = 0.0f;
    _S2347[int(3)] = 0.0f;
    _S2347[int(4)] = 0.0f;
    _S2347[int(5)] = 0.0f;
    _S2347[int(6)] = 0.0f;
    _S2347[int(7)] = 0.0f;
    _S2347[int(8)] = 0.0f;
    _S2347[int(9)] = 0.0f;
    _S2347[int(10)] = 0.0f;
    _S2347[int(11)] = 0.0f;
    _S2347[int(12)] = 0.0f;
    _S2347[int(13)] = 0.0f;
    _S2347[int(14)] = 0.0f;
    _S2347[int(15)] = 0.0f;
    _S2347[int(16)] = 0.0f;
    _S2347[int(17)] = 0.0f;
    _S2347[int(18)] = 0.0f;
    _S2347[int(19)] = 0.0f;
    _S2347[int(20)] = 0.0f;
    _S2347[int(21)] = 0.0f;
    _S2347[int(22)] = 0.0f;
    _S2347[int(17)] = _S2346;
    float _S2348 = _S2319 + _S2347[int(0)];
    float _S2349 = _S2320 + _S2347[int(1)];
    float _S2350 = _S2321 + _S2347[int(2)];
    float _S2351 = _S2322 + _S2347[int(3)];
    float _S2352 = _S2323 + _S2347[int(4)];
    float _S2353 = _S2324 + _S2347[int(5)];
    float _S2354 = _S2325 + _S2347[int(6)];
    float _S2355 = _S2326 + _S2347[int(7)];
    float _S2356 = _S2327 + _S2347[int(8)];
    float _S2357 = _S2328 + _S2347[int(9)];
    float _S2358 = _S2329 + _S2347[int(10)];
    float _S2359 = _S2330 + _S2347[int(11)];
    float _S2360 = _S2331 + _S2347[int(12)];
    float _S2361 = _S2332 + _S2347[int(13)];
    float _S2362 = _S2333 + _S2347[int(14)];
    float _S2363 = _S2334 + _S2347[int(15)];
    float _S2364 = _S2335 + _S2347[int(16)];
    float _S2365 = _S2336 + _S2347[int(17)];
    float _S2366 = _S2337 + _S2347[int(18)];
    float _S2367 = _S2338 + _S2347[int(19)];
    float _S2368 = _S2339 + _S2347[int(20)];
    float _S2369 = _S2340 + _S2347[int(21)];
    float _S2370 = _S2341 + _S2347[int(22)];
    if(_S2307)
    {
        float _S2371 = 200.0f * _S2317;
        float _S2372 = _S2308 * _S2371 + 0.5f * (_S2306 * _S2371);
        _S2308 = 0.0f;
        _S2311 = _S2372;
    }
    else
    {
        _S2308 = _S2317;
        _S2311 = 0.0f;
    }
    DiffPair_float_0 _S2373;
    (&_S2373)->primal_0 = _S2306;
    (&_S2373)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2373, _S2308);
    float _S2374 = (_S2373.differential_0 + _S2311) / _S2281;
    FixedArray<float, 23>  _S2375;
    _S2375[int(0)] = 0.0f;
    _S2375[int(1)] = 0.0f;
    _S2375[int(2)] = 0.0f;
    _S2375[int(3)] = 0.0f;
    _S2375[int(4)] = 0.0f;
    _S2375[int(5)] = 0.0f;
    _S2375[int(6)] = 0.0f;
    _S2375[int(7)] = 0.0f;
    _S2375[int(8)] = 0.0f;
    _S2375[int(9)] = 0.0f;
    _S2375[int(10)] = 0.0f;
    _S2375[int(11)] = 0.0f;
    _S2375[int(12)] = 0.0f;
    _S2375[int(13)] = 0.0f;
    _S2375[int(14)] = 0.0f;
    _S2375[int(15)] = 0.0f;
    _S2375[int(16)] = 0.0f;
    _S2375[int(17)] = 0.0f;
    _S2375[int(18)] = 0.0f;
    _S2375[int(19)] = 0.0f;
    _S2375[int(20)] = 0.0f;
    _S2375[int(21)] = 0.0f;
    _S2375[int(22)] = 0.0f;
    _S2375[int(16)] = _S2374;
    float _S2376 = _S2348 + _S2375[int(0)];
    float _S2377 = _S2349 + _S2375[int(1)];
    float _S2378 = _S2350 + _S2375[int(2)];
    float _S2379 = _S2351 + _S2375[int(3)];
    float _S2380 = _S2352 + _S2375[int(4)];
    float _S2381 = _S2353 + _S2375[int(5)];
    float _S2382 = _S2354 + _S2375[int(6)];
    float _S2383 = _S2355 + _S2375[int(7)];
    float _S2384 = _S2356 + _S2375[int(8)];
    float _S2385 = _S2357 + _S2375[int(9)];
    float _S2386 = _S2358 + _S2375[int(10)];
    float _S2387 = _S2359 + _S2375[int(11)];
    float _S2388 = _S2360 + _S2375[int(12)];
    float _S2389 = _S2361 + _S2375[int(13)];
    float _S2390 = _S2362 + _S2375[int(14)];
    float _S2391 = _S2363 + _S2375[int(15)];
    float _S2392 = _S2364 + _S2375[int(16)];
    float _S2393 = _S2365 + _S2375[int(17)];
    float _S2394 = _S2366 + _S2375[int(18)];
    float _S2395 = _S2367 + _S2375[int(19)];
    float _S2396 = _S2368 + _S2375[int(20)];
    float _S2397 = _S2369 + _S2375[int(21)];
    float _S2398 = _S2370 + _S2375[int(22)];
    if(_S2304)
    {
        float _S2399 = 200.0f * _S2317;
        float _S2400 = _S2305 * _S2399 + 0.5f * (_S2303 * _S2399);
        _S2305 = 0.0f;
        _S2308 = _S2400;
    }
    else
    {
        _S2305 = _S2317;
        _S2308 = 0.0f;
    }
    DiffPair_float_0 _S2401;
    (&_S2401)->primal_0 = _S2303;
    (&_S2401)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2401, _S2305);
    float _S2402 = (_S2401.differential_0 + _S2308) / _S2281;
    FixedArray<float, 23>  _S2403;
    _S2403[int(0)] = 0.0f;
    _S2403[int(1)] = 0.0f;
    _S2403[int(2)] = 0.0f;
    _S2403[int(3)] = 0.0f;
    _S2403[int(4)] = 0.0f;
    _S2403[int(5)] = 0.0f;
    _S2403[int(6)] = 0.0f;
    _S2403[int(7)] = 0.0f;
    _S2403[int(8)] = 0.0f;
    _S2403[int(9)] = 0.0f;
    _S2403[int(10)] = 0.0f;
    _S2403[int(11)] = 0.0f;
    _S2403[int(12)] = 0.0f;
    _S2403[int(13)] = 0.0f;
    _S2403[int(14)] = 0.0f;
    _S2403[int(15)] = 0.0f;
    _S2403[int(16)] = 0.0f;
    _S2403[int(17)] = 0.0f;
    _S2403[int(18)] = 0.0f;
    _S2403[int(19)] = 0.0f;
    _S2403[int(20)] = 0.0f;
    _S2403[int(21)] = 0.0f;
    _S2403[int(22)] = 0.0f;
    _S2403[int(15)] = _S2402;
    float _S2404 = _S2376 + _S2403[int(0)];
    float _S2405 = _S2377 + _S2403[int(1)];
    float _S2406 = _S2378 + _S2403[int(2)];
    float _S2407 = _S2379 + _S2403[int(3)];
    float _S2408 = _S2380 + _S2403[int(4)];
    float _S2409 = _S2381 + _S2403[int(5)];
    float _S2410 = _S2382 + _S2403[int(6)];
    float _S2411 = _S2383 + _S2403[int(7)];
    float _S2412 = _S2384 + _S2403[int(8)];
    float _S2413 = _S2385 + _S2403[int(9)];
    float _S2414 = _S2386 + _S2403[int(10)];
    float _S2415 = _S2387 + _S2403[int(11)];
    float _S2416 = _S2388 + _S2403[int(12)];
    float _S2417 = _S2389 + _S2403[int(13)];
    float _S2418 = _S2390 + _S2403[int(14)];
    float _S2419 = _S2391 + _S2403[int(15)];
    float _S2420 = _S2392 + _S2403[int(16)];
    float _S2421 = _S2393 + _S2403[int(17)];
    float _S2422 = _S2394 + _S2403[int(18)];
    float _S2423 = _S2395 + _S2403[int(19)];
    float _S2424 = _S2396 + _S2403[int(20)];
    float _S2425 = _S2397 + _S2403[int(21)];
    float _S2426 = _S2398 + _S2403[int(22)];
    if(_S2301)
    {
        float _S2427 = 200.0f * _S2317;
        float _S2428 = _S2302 * _S2427 + 0.5f * (_S2300 * _S2427);
        _S2302 = 0.0f;
        _S2305 = _S2428;
    }
    else
    {
        _S2302 = _S2317;
        _S2305 = 0.0f;
    }
    DiffPair_float_0 _S2429;
    (&_S2429)->primal_0 = _S2300;
    (&_S2429)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2429, _S2302);
    float _S2430 = (_S2429.differential_0 + _S2305) / _S2281;
    FixedArray<float, 23>  _S2431;
    _S2431[int(0)] = 0.0f;
    _S2431[int(1)] = 0.0f;
    _S2431[int(2)] = 0.0f;
    _S2431[int(3)] = 0.0f;
    _S2431[int(4)] = 0.0f;
    _S2431[int(5)] = 0.0f;
    _S2431[int(6)] = 0.0f;
    _S2431[int(7)] = 0.0f;
    _S2431[int(8)] = 0.0f;
    _S2431[int(9)] = 0.0f;
    _S2431[int(10)] = 0.0f;
    _S2431[int(11)] = 0.0f;
    _S2431[int(12)] = 0.0f;
    _S2431[int(13)] = 0.0f;
    _S2431[int(14)] = 0.0f;
    _S2431[int(15)] = 0.0f;
    _S2431[int(16)] = 0.0f;
    _S2431[int(17)] = 0.0f;
    _S2431[int(18)] = 0.0f;
    _S2431[int(19)] = 0.0f;
    _S2431[int(20)] = 0.0f;
    _S2431[int(21)] = 0.0f;
    _S2431[int(22)] = 0.0f;
    _S2431[int(14)] = _S2430;
    float _S2432 = _S2404 + _S2431[int(0)];
    float _S2433 = _S2405 + _S2431[int(1)];
    float _S2434 = _S2406 + _S2431[int(2)];
    float _S2435 = _S2407 + _S2431[int(3)];
    float _S2436 = _S2408 + _S2431[int(4)];
    float _S2437 = _S2409 + _S2431[int(5)];
    float _S2438 = _S2410 + _S2431[int(6)];
    float _S2439 = _S2411 + _S2431[int(7)];
    float _S2440 = _S2412 + _S2431[int(8)];
    float _S2441 = _S2413 + _S2431[int(9)];
    float _S2442 = _S2414 + _S2431[int(10)];
    float _S2443 = _S2415 + _S2431[int(11)];
    float _S2444 = _S2416 + _S2431[int(12)];
    float _S2445 = _S2417 + _S2431[int(13)];
    float _S2446 = _S2418 + _S2431[int(14)];
    float _S2447 = _S2419 + _S2431[int(15)];
    float _S2448 = _S2420 + _S2431[int(16)];
    float _S2449 = _S2421 + _S2431[int(17)];
    float _S2450 = _S2422 + _S2431[int(18)];
    float _S2451 = _S2423 + _S2431[int(19)];
    float _S2452 = _S2424 + _S2431[int(20)];
    float _S2453 = _S2425 + _S2431[int(21)];
    float _S2454 = _S2426 + _S2431[int(22)];
    if(_S2298)
    {
        float _S2455 = 200.0f * _S2317;
        float _S2456 = _S2299 * _S2455 + 0.5f * (_S2297 * _S2455);
        _S2299 = 0.0f;
        _S2302 = _S2456;
    }
    else
    {
        _S2299 = _S2317;
        _S2302 = 0.0f;
    }
    DiffPair_float_0 _S2457;
    (&_S2457)->primal_0 = _S2297;
    (&_S2457)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2457, _S2299);
    float _S2458 = (_S2457.differential_0 + _S2302) / _S2281;
    FixedArray<float, 23>  _S2459;
    _S2459[int(0)] = 0.0f;
    _S2459[int(1)] = 0.0f;
    _S2459[int(2)] = 0.0f;
    _S2459[int(3)] = 0.0f;
    _S2459[int(4)] = 0.0f;
    _S2459[int(5)] = 0.0f;
    _S2459[int(6)] = 0.0f;
    _S2459[int(7)] = 0.0f;
    _S2459[int(8)] = 0.0f;
    _S2459[int(9)] = 0.0f;
    _S2459[int(10)] = 0.0f;
    _S2459[int(11)] = 0.0f;
    _S2459[int(12)] = 0.0f;
    _S2459[int(13)] = 0.0f;
    _S2459[int(14)] = 0.0f;
    _S2459[int(15)] = 0.0f;
    _S2459[int(16)] = 0.0f;
    _S2459[int(17)] = 0.0f;
    _S2459[int(18)] = 0.0f;
    _S2459[int(19)] = 0.0f;
    _S2459[int(20)] = 0.0f;
    _S2459[int(21)] = 0.0f;
    _S2459[int(22)] = 0.0f;
    _S2459[int(13)] = _S2458;
    float _S2460 = _S2432 + _S2459[int(0)];
    float _S2461 = _S2433 + _S2459[int(1)];
    float _S2462 = _S2434 + _S2459[int(2)];
    float _S2463 = _S2435 + _S2459[int(3)];
    float _S2464 = _S2436 + _S2459[int(4)];
    float _S2465 = _S2437 + _S2459[int(5)];
    float _S2466 = _S2438 + _S2459[int(6)];
    float _S2467 = _S2439 + _S2459[int(7)];
    float _S2468 = _S2440 + _S2459[int(8)];
    float _S2469 = _S2441 + _S2459[int(9)];
    float _S2470 = _S2442 + _S2459[int(10)];
    float _S2471 = _S2443 + _S2459[int(11)];
    float _S2472 = _S2444 + _S2459[int(12)];
    float _S2473 = _S2445 + _S2459[int(13)];
    float _S2474 = _S2446 + _S2459[int(14)];
    float _S2475 = _S2447 + _S2459[int(15)];
    float _S2476 = _S2448 + _S2459[int(16)];
    float _S2477 = _S2449 + _S2459[int(17)];
    float _S2478 = _S2450 + _S2459[int(18)];
    float _S2479 = _S2451 + _S2459[int(19)];
    float _S2480 = _S2452 + _S2459[int(20)];
    float _S2481 = _S2453 + _S2459[int(21)];
    float _S2482 = _S2454 + _S2459[int(22)];
    if(_S2295)
    {
        float _S2483 = 200.0f * _S2317;
        float _S2484 = _S2296 * _S2483 + 0.5f * (_S2294 * _S2483);
        _S2296 = 0.0f;
        _S2299 = _S2484;
    }
    else
    {
        _S2296 = _S2317;
        _S2299 = 0.0f;
    }
    DiffPair_float_0 _S2485;
    (&_S2485)->primal_0 = _S2294;
    (&_S2485)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2485, _S2296);
    float _S2486 = (_S2485.differential_0 + _S2299) / _S2281;
    FixedArray<float, 23>  _S2487;
    _S2487[int(0)] = 0.0f;
    _S2487[int(1)] = 0.0f;
    _S2487[int(2)] = 0.0f;
    _S2487[int(3)] = 0.0f;
    _S2487[int(4)] = 0.0f;
    _S2487[int(5)] = 0.0f;
    _S2487[int(6)] = 0.0f;
    _S2487[int(7)] = 0.0f;
    _S2487[int(8)] = 0.0f;
    _S2487[int(9)] = 0.0f;
    _S2487[int(10)] = 0.0f;
    _S2487[int(11)] = 0.0f;
    _S2487[int(12)] = 0.0f;
    _S2487[int(13)] = 0.0f;
    _S2487[int(14)] = 0.0f;
    _S2487[int(15)] = 0.0f;
    _S2487[int(16)] = 0.0f;
    _S2487[int(17)] = 0.0f;
    _S2487[int(18)] = 0.0f;
    _S2487[int(19)] = 0.0f;
    _S2487[int(20)] = 0.0f;
    _S2487[int(21)] = 0.0f;
    _S2487[int(22)] = 0.0f;
    _S2487[int(12)] = _S2486;
    float _S2488 = _S2460 + _S2487[int(0)];
    float _S2489 = _S2461 + _S2487[int(1)];
    float _S2490 = _S2462 + _S2487[int(2)];
    float _S2491 = _S2463 + _S2487[int(3)];
    float _S2492 = _S2464 + _S2487[int(4)];
    float _S2493 = _S2465 + _S2487[int(5)];
    float _S2494 = _S2466 + _S2487[int(6)];
    float _S2495 = _S2467 + _S2487[int(7)];
    float _S2496 = _S2468 + _S2487[int(8)];
    float _S2497 = _S2469 + _S2487[int(9)];
    float _S2498 = _S2470 + _S2487[int(10)];
    float _S2499 = _S2471 + _S2487[int(11)];
    float _S2500 = _S2472 + _S2487[int(12)];
    float _S2501 = _S2473 + _S2487[int(13)];
    float _S2502 = _S2474 + _S2487[int(14)];
    float _S2503 = _S2475 + _S2487[int(15)];
    float _S2504 = _S2476 + _S2487[int(16)];
    float _S2505 = _S2477 + _S2487[int(17)];
    float _S2506 = _S2478 + _S2487[int(18)];
    float _S2507 = _S2479 + _S2487[int(19)];
    float _S2508 = _S2480 + _S2487[int(20)];
    float _S2509 = _S2481 + _S2487[int(21)];
    float _S2510 = _S2482 + _S2487[int(22)];
    if(_S2292)
    {
        float _S2511 = 200.0f * _S2317;
        float _S2512 = _S2293 * _S2511 + 0.5f * (_S2291 * _S2511);
        _S2293 = 0.0f;
        _S2296 = _S2512;
    }
    else
    {
        _S2293 = _S2317;
        _S2296 = 0.0f;
    }
    DiffPair_float_0 _S2513;
    (&_S2513)->primal_0 = _S2291;
    (&_S2513)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2513, _S2293);
    float _S2514 = (_S2513.differential_0 + _S2296) / _S2281;
    FixedArray<float, 23>  _S2515;
    _S2515[int(0)] = 0.0f;
    _S2515[int(1)] = 0.0f;
    _S2515[int(2)] = 0.0f;
    _S2515[int(3)] = 0.0f;
    _S2515[int(4)] = 0.0f;
    _S2515[int(5)] = 0.0f;
    _S2515[int(6)] = 0.0f;
    _S2515[int(7)] = 0.0f;
    _S2515[int(8)] = 0.0f;
    _S2515[int(9)] = 0.0f;
    _S2515[int(10)] = 0.0f;
    _S2515[int(11)] = 0.0f;
    _S2515[int(12)] = 0.0f;
    _S2515[int(13)] = 0.0f;
    _S2515[int(14)] = 0.0f;
    _S2515[int(15)] = 0.0f;
    _S2515[int(16)] = 0.0f;
    _S2515[int(17)] = 0.0f;
    _S2515[int(18)] = 0.0f;
    _S2515[int(19)] = 0.0f;
    _S2515[int(20)] = 0.0f;
    _S2515[int(21)] = 0.0f;
    _S2515[int(22)] = 0.0f;
    _S2515[int(11)] = _S2514;
    float _S2516 = _S2488 + _S2515[int(0)];
    float _S2517 = _S2489 + _S2515[int(1)];
    float _S2518 = _S2490 + _S2515[int(2)];
    float _S2519 = _S2491 + _S2515[int(3)];
    float _S2520 = _S2492 + _S2515[int(4)];
    float _S2521 = _S2493 + _S2515[int(5)];
    float _S2522 = _S2494 + _S2515[int(6)];
    float _S2523 = _S2495 + _S2515[int(7)];
    float _S2524 = _S2496 + _S2515[int(8)];
    float _S2525 = _S2497 + _S2515[int(9)];
    float _S2526 = _S2498 + _S2515[int(10)];
    float _S2527 = _S2499 + _S2515[int(11)];
    float _S2528 = _S2500 + _S2515[int(12)];
    float _S2529 = _S2501 + _S2515[int(13)];
    float _S2530 = _S2502 + _S2515[int(14)];
    float _S2531 = _S2503 + _S2515[int(15)];
    float _S2532 = _S2504 + _S2515[int(16)];
    float _S2533 = _S2505 + _S2515[int(17)];
    float _S2534 = _S2506 + _S2515[int(18)];
    float _S2535 = _S2507 + _S2515[int(19)];
    float _S2536 = _S2508 + _S2515[int(20)];
    float _S2537 = _S2509 + _S2515[int(21)];
    float _S2538 = _S2510 + _S2515[int(22)];
    if(_S2289)
    {
        float _S2539 = 200.0f * _S2317;
        float _S2540 = _S2290 * _S2539 + 0.5f * (_S2288 * _S2539);
        _S2290 = 0.0f;
        _S2293 = _S2540;
    }
    else
    {
        _S2290 = _S2317;
        _S2293 = 0.0f;
    }
    DiffPair_float_0 _S2541;
    (&_S2541)->primal_0 = _S2288;
    (&_S2541)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2541, _S2290);
    float _S2542 = (_S2541.differential_0 + _S2293) / _S2281;
    float _S2543 = _S2312 / _S2287;
    float _S2544 = _S2313 / _S2286;
    float _S2545 = _S2314 / _S2285;
    FixedArray<float, 23>  _S2546;
    _S2546[int(0)] = 0.0f;
    _S2546[int(1)] = 0.0f;
    _S2546[int(2)] = 0.0f;
    _S2546[int(3)] = 0.0f;
    _S2546[int(4)] = 0.0f;
    _S2546[int(5)] = 0.0f;
    _S2546[int(6)] = 0.0f;
    _S2546[int(7)] = 0.0f;
    _S2546[int(8)] = 0.0f;
    _S2546[int(9)] = 0.0f;
    _S2546[int(10)] = 0.0f;
    _S2546[int(11)] = 0.0f;
    _S2546[int(12)] = 0.0f;
    _S2546[int(13)] = 0.0f;
    _S2546[int(14)] = 0.0f;
    _S2546[int(15)] = 0.0f;
    _S2546[int(16)] = 0.0f;
    _S2546[int(17)] = 0.0f;
    _S2546[int(18)] = 0.0f;
    _S2546[int(19)] = 0.0f;
    _S2546[int(20)] = 0.0f;
    _S2546[int(21)] = 0.0f;
    _S2546[int(22)] = 0.0f;
    _S2546[int(10)] = _S2542;
    _S2546[int(9)] = _S2543;
    _S2546[int(8)] = _S2543;
    _S2546[int(7)] = _S2543;
    _S2546[int(6)] = _S2543;
    _S2546[int(5)] = _S2543;
    _S2546[int(4)] = _S2544;
    _S2546[int(3)] = _S2544;
    _S2546[int(2)] = _S2544;
    _S2546[int(1)] = _S2545;
    float _S2547 = _S2516 + _S2546[int(0)];
    float _S2548 = _S2517 + _S2546[int(1)];
    float _S2549 = _S2518 + _S2546[int(2)];
    float _S2550 = _S2519 + _S2546[int(3)];
    float _S2551 = _S2520 + _S2546[int(4)];
    float _S2552 = _S2521 + _S2546[int(5)];
    float _S2553 = _S2522 + _S2546[int(6)];
    float _S2554 = _S2523 + _S2546[int(7)];
    float _S2555 = _S2524 + _S2546[int(8)];
    float _S2556 = _S2525 + _S2546[int(9)];
    float _S2557 = _S2526 + _S2546[int(10)];
    float _S2558 = _S2527 + _S2546[int(11)];
    float _S2559 = _S2528 + _S2546[int(12)];
    float _S2560 = _S2529 + _S2546[int(13)];
    float _S2561 = _S2530 + _S2546[int(14)];
    float _S2562 = _S2531 + _S2546[int(15)];
    float _S2563 = _S2532 + _S2546[int(16)];
    float _S2564 = _S2533 + _S2546[int(17)];
    float _S2565 = _S2534 + _S2546[int(18)];
    float _S2566 = _S2535 + _S2546[int(19)];
    float _S2567 = _S2536 + _S2546[int(20)];
    float _S2568 = _S2537 + _S2546[int(21)];
    float _S2569 = _S2538 + _S2546[int(22)];
    if(_S2283)
    {
        float _S2570 = 10.0f * _S2315;
        float _S2571 = _S2284 * _S2570 + 0.5f * (_S2282 * _S2570);
        _S2284 = 0.0f;
        _S2290 = _S2571;
    }
    else
    {
        _S2284 = _S2315;
        _S2290 = 0.0f;
    }
    DiffPair_float_0 _S2572;
    (&_S2572)->primal_0 = _S2282;
    (&_S2572)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2572, _S2284);
    float _S2573 = (_S2572.differential_0 + _S2290) / _S2281;
    FixedArray<float, 23>  _S2574;
    _S2574[int(0)] = 0.0f;
    _S2574[int(1)] = 0.0f;
    _S2574[int(2)] = 0.0f;
    _S2574[int(3)] = 0.0f;
    _S2574[int(4)] = 0.0f;
    _S2574[int(5)] = 0.0f;
    _S2574[int(6)] = 0.0f;
    _S2574[int(7)] = 0.0f;
    _S2574[int(8)] = 0.0f;
    _S2574[int(9)] = 0.0f;
    _S2574[int(10)] = 0.0f;
    _S2574[int(11)] = 0.0f;
    _S2574[int(12)] = 0.0f;
    _S2574[int(13)] = 0.0f;
    _S2574[int(14)] = 0.0f;
    _S2574[int(15)] = 0.0f;
    _S2574[int(16)] = 0.0f;
    _S2574[int(17)] = 0.0f;
    _S2574[int(18)] = 0.0f;
    _S2574[int(19)] = 0.0f;
    _S2574[int(20)] = 0.0f;
    _S2574[int(21)] = 0.0f;
    _S2574[int(22)] = 0.0f;
    _S2574[int(0)] = _S2573;
    FixedArray<float, 23>  _S2575 = {
        _S2547 + _S2574[int(0)], _S2548 + _S2574[int(1)], _S2549 + _S2574[int(2)], _S2550 + _S2574[int(3)], _S2551 + _S2574[int(4)], _S2552 + _S2574[int(5)], _S2553 + _S2574[int(6)], _S2554 + _S2574[int(7)], _S2555 + _S2574[int(8)], _S2556 + _S2574[int(9)], _S2557 + _S2574[int(10)], _S2558 + _S2574[int(11)], _S2559 + _S2574[int(12)], _S2560 + _S2574[int(13)], _S2561 + _S2574[int(14)], _S2562 + _S2574[int(15)], _S2563 + _S2574[int(16)], _S2564 + _S2574[int(17)], _S2565 + _S2574[int(18)], _S2566 + _S2574[int(19)], _S2567 + _S2574[int(20)], _S2568 + _S2574[int(21)], _S2569 + _S2574[int(22)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2575;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S2576, int _S2577, FixedArray<float, 6>  * _S2578, FixedArray<float, 6>  * _S2579)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S2576, _S2577, _S2578, _S2579);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_3, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_7, FixedArray<float, 23>  * _S2580)
{
    FixedArray<float, 23>  _S2581 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_3;
    (&dp_raw_losses_1)->differential_0 = _S2581;
    FixedArray<float, 6>  _S2582 = loss_weights_5;
    FixedArray<float, 6>  _S2583 = grad_out_7;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_1, num_cameras_5, &_S2582, &_S2583);
    *_S2580 = (&dp_raw_losses_1)->differential_0;
    return;
}

inline __device__ void compute_ppisp_no_crf_regularization_loss(FixedArray<float, 18>  raw_losses_4, int num_cameras_6, FixedArray<float, 6>  loss_weights_6, FixedArray<float, 6>  * _S2584)
{
    float _S2585;
    FixedArray<float, 6>  losses_5;
    float _S2586 = float(num_cameras_6);
    float _S2587 = raw_losses_4[int(0)] / _S2586;
    for(;;)
    {
        float _S2588 = (F32_abs((_S2587)));
        if(_S2588 < 0.10000000149011612f)
        {
            _S2585 = 0.5f * _S2587 * _S2587 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S2585 = _S2588 - 0.05000000074505806f;
            break;
        }
    }
    losses_5[int(0)] = _S2585;
    losses_5[int(1)] = raw_losses_4[int(1)] / (3.0f * _S2586);
    losses_5[int(2)] = (raw_losses_4[int(2)] + raw_losses_4[int(3)] + raw_losses_4[int(4)]) / (9.0f * _S2586);
    losses_5[int(3)] = (raw_losses_4[int(5)] + raw_losses_4[int(6)] + raw_losses_4[int(7)] + raw_losses_4[int(8)] + raw_losses_4[int(9)]) / (5.0f * _S2586);
    float _S2589 = raw_losses_4[int(10)] / _S2586;
    for(;;)
    {
        float _S2590 = (F32_abs((_S2589)));
        if(_S2590 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2589 * _S2589 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2590 - 0.00249999994412065f;
            break;
        }
    }
    float _S2591;
    float _S2592 = raw_losses_4[int(11)] / _S2586;
    for(;;)
    {
        float _S2593 = (F32_abs((_S2592)));
        if(_S2593 < 0.00499999988824129f)
        {
            _S2591 = 0.5f * _S2592 * _S2592 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2591 = _S2593 - 0.00249999994412065f;
            break;
        }
    }
    float _S2594 = _S2585 + _S2591;
    float _S2595 = raw_losses_4[int(12)] / _S2586;
    for(;;)
    {
        float _S2596 = (F32_abs((_S2595)));
        if(_S2596 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2595 * _S2595 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2596 - 0.00249999994412065f;
            break;
        }
    }
    float _S2597 = _S2594 + _S2585;
    float _S2598 = raw_losses_4[int(13)] / _S2586;
    for(;;)
    {
        float _S2599 = (F32_abs((_S2598)));
        if(_S2599 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2598 * _S2598 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2599 - 0.00249999994412065f;
            break;
        }
    }
    float _S2600 = _S2597 + _S2585;
    float _S2601 = raw_losses_4[int(14)] / _S2586;
    for(;;)
    {
        float _S2602 = (F32_abs((_S2601)));
        if(_S2602 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2601 * _S2601 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2602 - 0.00249999994412065f;
            break;
        }
    }
    float _S2603 = _S2600 + _S2585;
    float _S2604 = raw_losses_4[int(15)] / _S2586;
    for(;;)
    {
        float _S2605 = (F32_abs((_S2604)));
        if(_S2605 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2604 * _S2604 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2605 - 0.00249999994412065f;
            break;
        }
    }
    float _S2606 = _S2603 + _S2585;
    float _S2607 = raw_losses_4[int(16)] / _S2586;
    for(;;)
    {
        float _S2608 = (F32_abs((_S2607)));
        if(_S2608 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2607 * _S2607 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2608 - 0.00249999994412065f;
            break;
        }
    }
    float _S2609 = _S2606 + _S2585;
    float _S2610 = raw_losses_4[int(17)] / _S2586;
    for(;;)
    {
        float _S2611 = (F32_abs((_S2610)));
        if(_S2611 < 0.00499999988824129f)
        {
            _S2585 = 0.5f * _S2610 * _S2610 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S2585 = _S2611 - 0.00249999994412065f;
            break;
        }
    }
    float _S2612 = (_S2609 + _S2585) / 8.0f;
    losses_5[int(5)] = 0.0f;
    losses_5[int(0)] = losses_5[int(0)] * loss_weights_6[int(0)];
    losses_5[int(1)] = losses_5[int(1)] * loss_weights_6[int(1)];
    losses_5[int(2)] = losses_5[int(2)] * loss_weights_6[int(2)];
    losses_5[int(3)] = losses_5[int(3)] * loss_weights_6[int(3)];
    losses_5[int(4)] = _S2612 * loss_weights_6[int(4)];
    *_S2584 = losses_5;
    return;
}

struct DiffPair_arrayx3Cfloatx2C18x3E_0
{
    FixedArray<float, 18>  primal_0;
    FixedArray<float, 18>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_no_crf_regularization_loss_0(DiffPair_arrayx3Cfloatx2C18x3E_0 * dpraw_losses_2, int num_cameras_7, FixedArray<float, 6>  * loss_weights_7, FixedArray<float, 6>  * _s_dOut_8)
{
    FixedArray<float, 18>  _S2613 = dpraw_losses_2->primal_0;
    float _S2614 = float(num_cameras_7);
    float _S2615 = dpraw_losses_2->primal_0[int(0)] / _S2614;
    bool _S2616 = (s_primal_ctx_abs_0(_S2615)) < 0.10000000149011612f;
    float _S2617;
    if(_S2616)
    {
        _S2617 = 0.5f * _S2615;
    }
    else
    {
        _S2617 = 0.0f;
    }
    float _S2618 = 3.0f * _S2614;
    float _S2619 = 9.0f * _S2614;
    float _S2620 = 5.0f * _S2614;
    float _S2621 = _S2613[int(10)] / _S2614;
    bool _S2622 = (s_primal_ctx_abs_0(_S2621)) < 0.00499999988824129f;
    float _S2623;
    if(_S2622)
    {
        _S2623 = 0.5f * _S2621;
    }
    else
    {
        _S2623 = 0.0f;
    }
    float _S2624 = _S2613[int(11)] / _S2614;
    bool _S2625 = (s_primal_ctx_abs_0(_S2624)) < 0.00499999988824129f;
    float _S2626;
    if(_S2625)
    {
        _S2626 = 0.5f * _S2624;
    }
    else
    {
        _S2626 = 0.0f;
    }
    float _S2627 = _S2613[int(12)] / _S2614;
    bool _S2628 = (s_primal_ctx_abs_0(_S2627)) < 0.00499999988824129f;
    float _S2629;
    if(_S2628)
    {
        _S2629 = 0.5f * _S2627;
    }
    else
    {
        _S2629 = 0.0f;
    }
    float _S2630 = _S2613[int(13)] / _S2614;
    bool _S2631 = (s_primal_ctx_abs_0(_S2630)) < 0.00499999988824129f;
    float _S2632;
    if(_S2631)
    {
        _S2632 = 0.5f * _S2630;
    }
    else
    {
        _S2632 = 0.0f;
    }
    float _S2633 = _S2613[int(14)] / _S2614;
    bool _S2634 = (s_primal_ctx_abs_0(_S2633)) < 0.00499999988824129f;
    float _S2635;
    if(_S2634)
    {
        _S2635 = 0.5f * _S2633;
    }
    else
    {
        _S2635 = 0.0f;
    }
    float _S2636 = _S2613[int(15)] / _S2614;
    bool _S2637 = (s_primal_ctx_abs_0(_S2636)) < 0.00499999988824129f;
    float _S2638;
    if(_S2637)
    {
        _S2638 = 0.5f * _S2636;
    }
    else
    {
        _S2638 = 0.0f;
    }
    float _S2639 = _S2613[int(16)] / _S2614;
    bool _S2640 = (s_primal_ctx_abs_0(_S2639)) < 0.00499999988824129f;
    float _S2641;
    if(_S2640)
    {
        _S2641 = 0.5f * _S2639;
    }
    else
    {
        _S2641 = 0.0f;
    }
    float _S2642 = _S2613[int(17)] / _S2614;
    bool _S2643 = (s_primal_ctx_abs_0(_S2642)) < 0.00499999988824129f;
    float _S2644;
    if(_S2643)
    {
        _S2644 = 0.5f * _S2642;
    }
    else
    {
        _S2644 = 0.0f;
    }
    float _S2645 = (*loss_weights_7)[int(3)] * (*_s_dOut_8)[int(3)];
    float _S2646 = (*loss_weights_7)[int(2)] * (*_s_dOut_8)[int(2)];
    float _S2647 = (*loss_weights_7)[int(1)] * (*_s_dOut_8)[int(1)];
    float _S2648 = (*loss_weights_7)[int(0)] * (*_s_dOut_8)[int(0)];
    float _S2649 = 0.125f * ((*loss_weights_7)[int(4)] * (*_s_dOut_8)[int(4)]);
    float _S2650;
    if(_S2643)
    {
        float _S2651 = 200.0f * _S2649;
        float _S2652 = _S2644 * _S2651 + 0.5f * (_S2642 * _S2651);
        _S2644 = 0.0f;
        _S2650 = _S2652;
    }
    else
    {
        _S2644 = _S2649;
        _S2650 = 0.0f;
    }
    DiffPair_float_0 _S2653;
    (&_S2653)->primal_0 = _S2642;
    (&_S2653)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2653, _S2644);
    float _S2654 = (_S2653.differential_0 + _S2650) / _S2614;
    FixedArray<float, 18>  _S2655;
    _S2655[int(0)] = 0.0f;
    _S2655[int(1)] = 0.0f;
    _S2655[int(2)] = 0.0f;
    _S2655[int(3)] = 0.0f;
    _S2655[int(4)] = 0.0f;
    _S2655[int(5)] = 0.0f;
    _S2655[int(6)] = 0.0f;
    _S2655[int(7)] = 0.0f;
    _S2655[int(8)] = 0.0f;
    _S2655[int(9)] = 0.0f;
    _S2655[int(10)] = 0.0f;
    _S2655[int(11)] = 0.0f;
    _S2655[int(12)] = 0.0f;
    _S2655[int(13)] = 0.0f;
    _S2655[int(14)] = 0.0f;
    _S2655[int(15)] = 0.0f;
    _S2655[int(16)] = 0.0f;
    _S2655[int(17)] = 0.0f;
    _S2655[int(17)] = _S2654;
    float _S2656 = _S2655[int(0)];
    float _S2657 = _S2655[int(1)];
    float _S2658 = _S2655[int(2)];
    float _S2659 = _S2655[int(3)];
    float _S2660 = _S2655[int(4)];
    float _S2661 = _S2655[int(5)];
    float _S2662 = _S2655[int(6)];
    float _S2663 = _S2655[int(7)];
    float _S2664 = _S2655[int(8)];
    float _S2665 = _S2655[int(9)];
    float _S2666 = _S2655[int(10)];
    float _S2667 = _S2655[int(11)];
    float _S2668 = _S2655[int(12)];
    float _S2669 = _S2655[int(13)];
    float _S2670 = _S2655[int(14)];
    float _S2671 = _S2655[int(15)];
    float _S2672 = _S2655[int(16)];
    float _S2673 = _S2655[int(17)];
    if(_S2640)
    {
        float _S2674 = 200.0f * _S2649;
        float _S2675 = _S2641 * _S2674 + 0.5f * (_S2639 * _S2674);
        _S2641 = 0.0f;
        _S2644 = _S2675;
    }
    else
    {
        _S2641 = _S2649;
        _S2644 = 0.0f;
    }
    DiffPair_float_0 _S2676;
    (&_S2676)->primal_0 = _S2639;
    (&_S2676)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2676, _S2641);
    float _S2677 = (_S2676.differential_0 + _S2644) / _S2614;
    FixedArray<float, 18>  _S2678;
    _S2678[int(0)] = 0.0f;
    _S2678[int(1)] = 0.0f;
    _S2678[int(2)] = 0.0f;
    _S2678[int(3)] = 0.0f;
    _S2678[int(4)] = 0.0f;
    _S2678[int(5)] = 0.0f;
    _S2678[int(6)] = 0.0f;
    _S2678[int(7)] = 0.0f;
    _S2678[int(8)] = 0.0f;
    _S2678[int(9)] = 0.0f;
    _S2678[int(10)] = 0.0f;
    _S2678[int(11)] = 0.0f;
    _S2678[int(12)] = 0.0f;
    _S2678[int(13)] = 0.0f;
    _S2678[int(14)] = 0.0f;
    _S2678[int(15)] = 0.0f;
    _S2678[int(16)] = 0.0f;
    _S2678[int(17)] = 0.0f;
    _S2678[int(16)] = _S2677;
    float _S2679 = _S2656 + _S2678[int(0)];
    float _S2680 = _S2657 + _S2678[int(1)];
    float _S2681 = _S2658 + _S2678[int(2)];
    float _S2682 = _S2659 + _S2678[int(3)];
    float _S2683 = _S2660 + _S2678[int(4)];
    float _S2684 = _S2661 + _S2678[int(5)];
    float _S2685 = _S2662 + _S2678[int(6)];
    float _S2686 = _S2663 + _S2678[int(7)];
    float _S2687 = _S2664 + _S2678[int(8)];
    float _S2688 = _S2665 + _S2678[int(9)];
    float _S2689 = _S2666 + _S2678[int(10)];
    float _S2690 = _S2667 + _S2678[int(11)];
    float _S2691 = _S2668 + _S2678[int(12)];
    float _S2692 = _S2669 + _S2678[int(13)];
    float _S2693 = _S2670 + _S2678[int(14)];
    float _S2694 = _S2671 + _S2678[int(15)];
    float _S2695 = _S2672 + _S2678[int(16)];
    float _S2696 = _S2673 + _S2678[int(17)];
    if(_S2637)
    {
        float _S2697 = 200.0f * _S2649;
        float _S2698 = _S2638 * _S2697 + 0.5f * (_S2636 * _S2697);
        _S2638 = 0.0f;
        _S2641 = _S2698;
    }
    else
    {
        _S2638 = _S2649;
        _S2641 = 0.0f;
    }
    DiffPair_float_0 _S2699;
    (&_S2699)->primal_0 = _S2636;
    (&_S2699)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2699, _S2638);
    float _S2700 = (_S2699.differential_0 + _S2641) / _S2614;
    FixedArray<float, 18>  _S2701;
    _S2701[int(0)] = 0.0f;
    _S2701[int(1)] = 0.0f;
    _S2701[int(2)] = 0.0f;
    _S2701[int(3)] = 0.0f;
    _S2701[int(4)] = 0.0f;
    _S2701[int(5)] = 0.0f;
    _S2701[int(6)] = 0.0f;
    _S2701[int(7)] = 0.0f;
    _S2701[int(8)] = 0.0f;
    _S2701[int(9)] = 0.0f;
    _S2701[int(10)] = 0.0f;
    _S2701[int(11)] = 0.0f;
    _S2701[int(12)] = 0.0f;
    _S2701[int(13)] = 0.0f;
    _S2701[int(14)] = 0.0f;
    _S2701[int(15)] = 0.0f;
    _S2701[int(16)] = 0.0f;
    _S2701[int(17)] = 0.0f;
    _S2701[int(15)] = _S2700;
    float _S2702 = _S2679 + _S2701[int(0)];
    float _S2703 = _S2680 + _S2701[int(1)];
    float _S2704 = _S2681 + _S2701[int(2)];
    float _S2705 = _S2682 + _S2701[int(3)];
    float _S2706 = _S2683 + _S2701[int(4)];
    float _S2707 = _S2684 + _S2701[int(5)];
    float _S2708 = _S2685 + _S2701[int(6)];
    float _S2709 = _S2686 + _S2701[int(7)];
    float _S2710 = _S2687 + _S2701[int(8)];
    float _S2711 = _S2688 + _S2701[int(9)];
    float _S2712 = _S2689 + _S2701[int(10)];
    float _S2713 = _S2690 + _S2701[int(11)];
    float _S2714 = _S2691 + _S2701[int(12)];
    float _S2715 = _S2692 + _S2701[int(13)];
    float _S2716 = _S2693 + _S2701[int(14)];
    float _S2717 = _S2694 + _S2701[int(15)];
    float _S2718 = _S2695 + _S2701[int(16)];
    float _S2719 = _S2696 + _S2701[int(17)];
    if(_S2634)
    {
        float _S2720 = 200.0f * _S2649;
        float _S2721 = _S2635 * _S2720 + 0.5f * (_S2633 * _S2720);
        _S2635 = 0.0f;
        _S2638 = _S2721;
    }
    else
    {
        _S2635 = _S2649;
        _S2638 = 0.0f;
    }
    DiffPair_float_0 _S2722;
    (&_S2722)->primal_0 = _S2633;
    (&_S2722)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2722, _S2635);
    float _S2723 = (_S2722.differential_0 + _S2638) / _S2614;
    FixedArray<float, 18>  _S2724;
    _S2724[int(0)] = 0.0f;
    _S2724[int(1)] = 0.0f;
    _S2724[int(2)] = 0.0f;
    _S2724[int(3)] = 0.0f;
    _S2724[int(4)] = 0.0f;
    _S2724[int(5)] = 0.0f;
    _S2724[int(6)] = 0.0f;
    _S2724[int(7)] = 0.0f;
    _S2724[int(8)] = 0.0f;
    _S2724[int(9)] = 0.0f;
    _S2724[int(10)] = 0.0f;
    _S2724[int(11)] = 0.0f;
    _S2724[int(12)] = 0.0f;
    _S2724[int(13)] = 0.0f;
    _S2724[int(14)] = 0.0f;
    _S2724[int(15)] = 0.0f;
    _S2724[int(16)] = 0.0f;
    _S2724[int(17)] = 0.0f;
    _S2724[int(14)] = _S2723;
    float _S2725 = _S2702 + _S2724[int(0)];
    float _S2726 = _S2703 + _S2724[int(1)];
    float _S2727 = _S2704 + _S2724[int(2)];
    float _S2728 = _S2705 + _S2724[int(3)];
    float _S2729 = _S2706 + _S2724[int(4)];
    float _S2730 = _S2707 + _S2724[int(5)];
    float _S2731 = _S2708 + _S2724[int(6)];
    float _S2732 = _S2709 + _S2724[int(7)];
    float _S2733 = _S2710 + _S2724[int(8)];
    float _S2734 = _S2711 + _S2724[int(9)];
    float _S2735 = _S2712 + _S2724[int(10)];
    float _S2736 = _S2713 + _S2724[int(11)];
    float _S2737 = _S2714 + _S2724[int(12)];
    float _S2738 = _S2715 + _S2724[int(13)];
    float _S2739 = _S2716 + _S2724[int(14)];
    float _S2740 = _S2717 + _S2724[int(15)];
    float _S2741 = _S2718 + _S2724[int(16)];
    float _S2742 = _S2719 + _S2724[int(17)];
    if(_S2631)
    {
        float _S2743 = 200.0f * _S2649;
        float _S2744 = _S2632 * _S2743 + 0.5f * (_S2630 * _S2743);
        _S2632 = 0.0f;
        _S2635 = _S2744;
    }
    else
    {
        _S2632 = _S2649;
        _S2635 = 0.0f;
    }
    DiffPair_float_0 _S2745;
    (&_S2745)->primal_0 = _S2630;
    (&_S2745)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2745, _S2632);
    float _S2746 = (_S2745.differential_0 + _S2635) / _S2614;
    FixedArray<float, 18>  _S2747;
    _S2747[int(0)] = 0.0f;
    _S2747[int(1)] = 0.0f;
    _S2747[int(2)] = 0.0f;
    _S2747[int(3)] = 0.0f;
    _S2747[int(4)] = 0.0f;
    _S2747[int(5)] = 0.0f;
    _S2747[int(6)] = 0.0f;
    _S2747[int(7)] = 0.0f;
    _S2747[int(8)] = 0.0f;
    _S2747[int(9)] = 0.0f;
    _S2747[int(10)] = 0.0f;
    _S2747[int(11)] = 0.0f;
    _S2747[int(12)] = 0.0f;
    _S2747[int(13)] = 0.0f;
    _S2747[int(14)] = 0.0f;
    _S2747[int(15)] = 0.0f;
    _S2747[int(16)] = 0.0f;
    _S2747[int(17)] = 0.0f;
    _S2747[int(13)] = _S2746;
    float _S2748 = _S2725 + _S2747[int(0)];
    float _S2749 = _S2726 + _S2747[int(1)];
    float _S2750 = _S2727 + _S2747[int(2)];
    float _S2751 = _S2728 + _S2747[int(3)];
    float _S2752 = _S2729 + _S2747[int(4)];
    float _S2753 = _S2730 + _S2747[int(5)];
    float _S2754 = _S2731 + _S2747[int(6)];
    float _S2755 = _S2732 + _S2747[int(7)];
    float _S2756 = _S2733 + _S2747[int(8)];
    float _S2757 = _S2734 + _S2747[int(9)];
    float _S2758 = _S2735 + _S2747[int(10)];
    float _S2759 = _S2736 + _S2747[int(11)];
    float _S2760 = _S2737 + _S2747[int(12)];
    float _S2761 = _S2738 + _S2747[int(13)];
    float _S2762 = _S2739 + _S2747[int(14)];
    float _S2763 = _S2740 + _S2747[int(15)];
    float _S2764 = _S2741 + _S2747[int(16)];
    float _S2765 = _S2742 + _S2747[int(17)];
    if(_S2628)
    {
        float _S2766 = 200.0f * _S2649;
        float _S2767 = _S2629 * _S2766 + 0.5f * (_S2627 * _S2766);
        _S2629 = 0.0f;
        _S2632 = _S2767;
    }
    else
    {
        _S2629 = _S2649;
        _S2632 = 0.0f;
    }
    DiffPair_float_0 _S2768;
    (&_S2768)->primal_0 = _S2627;
    (&_S2768)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2768, _S2629);
    float _S2769 = (_S2768.differential_0 + _S2632) / _S2614;
    FixedArray<float, 18>  _S2770;
    _S2770[int(0)] = 0.0f;
    _S2770[int(1)] = 0.0f;
    _S2770[int(2)] = 0.0f;
    _S2770[int(3)] = 0.0f;
    _S2770[int(4)] = 0.0f;
    _S2770[int(5)] = 0.0f;
    _S2770[int(6)] = 0.0f;
    _S2770[int(7)] = 0.0f;
    _S2770[int(8)] = 0.0f;
    _S2770[int(9)] = 0.0f;
    _S2770[int(10)] = 0.0f;
    _S2770[int(11)] = 0.0f;
    _S2770[int(12)] = 0.0f;
    _S2770[int(13)] = 0.0f;
    _S2770[int(14)] = 0.0f;
    _S2770[int(15)] = 0.0f;
    _S2770[int(16)] = 0.0f;
    _S2770[int(17)] = 0.0f;
    _S2770[int(12)] = _S2769;
    float _S2771 = _S2748 + _S2770[int(0)];
    float _S2772 = _S2749 + _S2770[int(1)];
    float _S2773 = _S2750 + _S2770[int(2)];
    float _S2774 = _S2751 + _S2770[int(3)];
    float _S2775 = _S2752 + _S2770[int(4)];
    float _S2776 = _S2753 + _S2770[int(5)];
    float _S2777 = _S2754 + _S2770[int(6)];
    float _S2778 = _S2755 + _S2770[int(7)];
    float _S2779 = _S2756 + _S2770[int(8)];
    float _S2780 = _S2757 + _S2770[int(9)];
    float _S2781 = _S2758 + _S2770[int(10)];
    float _S2782 = _S2759 + _S2770[int(11)];
    float _S2783 = _S2760 + _S2770[int(12)];
    float _S2784 = _S2761 + _S2770[int(13)];
    float _S2785 = _S2762 + _S2770[int(14)];
    float _S2786 = _S2763 + _S2770[int(15)];
    float _S2787 = _S2764 + _S2770[int(16)];
    float _S2788 = _S2765 + _S2770[int(17)];
    if(_S2625)
    {
        float _S2789 = 200.0f * _S2649;
        float _S2790 = _S2626 * _S2789 + 0.5f * (_S2624 * _S2789);
        _S2626 = 0.0f;
        _S2629 = _S2790;
    }
    else
    {
        _S2626 = _S2649;
        _S2629 = 0.0f;
    }
    DiffPair_float_0 _S2791;
    (&_S2791)->primal_0 = _S2624;
    (&_S2791)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2791, _S2626);
    float _S2792 = (_S2791.differential_0 + _S2629) / _S2614;
    FixedArray<float, 18>  _S2793;
    _S2793[int(0)] = 0.0f;
    _S2793[int(1)] = 0.0f;
    _S2793[int(2)] = 0.0f;
    _S2793[int(3)] = 0.0f;
    _S2793[int(4)] = 0.0f;
    _S2793[int(5)] = 0.0f;
    _S2793[int(6)] = 0.0f;
    _S2793[int(7)] = 0.0f;
    _S2793[int(8)] = 0.0f;
    _S2793[int(9)] = 0.0f;
    _S2793[int(10)] = 0.0f;
    _S2793[int(11)] = 0.0f;
    _S2793[int(12)] = 0.0f;
    _S2793[int(13)] = 0.0f;
    _S2793[int(14)] = 0.0f;
    _S2793[int(15)] = 0.0f;
    _S2793[int(16)] = 0.0f;
    _S2793[int(17)] = 0.0f;
    _S2793[int(11)] = _S2792;
    float _S2794 = _S2771 + _S2793[int(0)];
    float _S2795 = _S2772 + _S2793[int(1)];
    float _S2796 = _S2773 + _S2793[int(2)];
    float _S2797 = _S2774 + _S2793[int(3)];
    float _S2798 = _S2775 + _S2793[int(4)];
    float _S2799 = _S2776 + _S2793[int(5)];
    float _S2800 = _S2777 + _S2793[int(6)];
    float _S2801 = _S2778 + _S2793[int(7)];
    float _S2802 = _S2779 + _S2793[int(8)];
    float _S2803 = _S2780 + _S2793[int(9)];
    float _S2804 = _S2781 + _S2793[int(10)];
    float _S2805 = _S2782 + _S2793[int(11)];
    float _S2806 = _S2783 + _S2793[int(12)];
    float _S2807 = _S2784 + _S2793[int(13)];
    float _S2808 = _S2785 + _S2793[int(14)];
    float _S2809 = _S2786 + _S2793[int(15)];
    float _S2810 = _S2787 + _S2793[int(16)];
    float _S2811 = _S2788 + _S2793[int(17)];
    if(_S2622)
    {
        float _S2812 = 200.0f * _S2649;
        float _S2813 = _S2623 * _S2812 + 0.5f * (_S2621 * _S2812);
        _S2623 = 0.0f;
        _S2626 = _S2813;
    }
    else
    {
        _S2623 = _S2649;
        _S2626 = 0.0f;
    }
    DiffPair_float_0 _S2814;
    (&_S2814)->primal_0 = _S2621;
    (&_S2814)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2814, _S2623);
    float _S2815 = (_S2814.differential_0 + _S2626) / _S2614;
    float _S2816 = _S2645 / _S2620;
    float _S2817 = _S2646 / _S2619;
    float _S2818 = _S2647 / _S2618;
    FixedArray<float, 18>  _S2819;
    _S2819[int(0)] = 0.0f;
    _S2819[int(1)] = 0.0f;
    _S2819[int(2)] = 0.0f;
    _S2819[int(3)] = 0.0f;
    _S2819[int(4)] = 0.0f;
    _S2819[int(5)] = 0.0f;
    _S2819[int(6)] = 0.0f;
    _S2819[int(7)] = 0.0f;
    _S2819[int(8)] = 0.0f;
    _S2819[int(9)] = 0.0f;
    _S2819[int(10)] = 0.0f;
    _S2819[int(11)] = 0.0f;
    _S2819[int(12)] = 0.0f;
    _S2819[int(13)] = 0.0f;
    _S2819[int(14)] = 0.0f;
    _S2819[int(15)] = 0.0f;
    _S2819[int(16)] = 0.0f;
    _S2819[int(17)] = 0.0f;
    _S2819[int(10)] = _S2815;
    _S2819[int(9)] = _S2816;
    _S2819[int(8)] = _S2816;
    _S2819[int(7)] = _S2816;
    _S2819[int(6)] = _S2816;
    _S2819[int(5)] = _S2816;
    _S2819[int(4)] = _S2817;
    _S2819[int(3)] = _S2817;
    _S2819[int(2)] = _S2817;
    _S2819[int(1)] = _S2818;
    float _S2820 = _S2794 + _S2819[int(0)];
    float _S2821 = _S2795 + _S2819[int(1)];
    float _S2822 = _S2796 + _S2819[int(2)];
    float _S2823 = _S2797 + _S2819[int(3)];
    float _S2824 = _S2798 + _S2819[int(4)];
    float _S2825 = _S2799 + _S2819[int(5)];
    float _S2826 = _S2800 + _S2819[int(6)];
    float _S2827 = _S2801 + _S2819[int(7)];
    float _S2828 = _S2802 + _S2819[int(8)];
    float _S2829 = _S2803 + _S2819[int(9)];
    float _S2830 = _S2804 + _S2819[int(10)];
    float _S2831 = _S2805 + _S2819[int(11)];
    float _S2832 = _S2806 + _S2819[int(12)];
    float _S2833 = _S2807 + _S2819[int(13)];
    float _S2834 = _S2808 + _S2819[int(14)];
    float _S2835 = _S2809 + _S2819[int(15)];
    float _S2836 = _S2810 + _S2819[int(16)];
    float _S2837 = _S2811 + _S2819[int(17)];
    if(_S2616)
    {
        float _S2838 = 10.0f * _S2648;
        float _S2839 = _S2617 * _S2838 + 0.5f * (_S2615 * _S2838);
        _S2617 = 0.0f;
        _S2623 = _S2839;
    }
    else
    {
        _S2617 = _S2648;
        _S2623 = 0.0f;
    }
    DiffPair_float_0 _S2840;
    (&_S2840)->primal_0 = _S2615;
    (&_S2840)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2840, _S2617);
    float _S2841 = (_S2840.differential_0 + _S2623) / _S2614;
    FixedArray<float, 18>  _S2842;
    _S2842[int(0)] = 0.0f;
    _S2842[int(1)] = 0.0f;
    _S2842[int(2)] = 0.0f;
    _S2842[int(3)] = 0.0f;
    _S2842[int(4)] = 0.0f;
    _S2842[int(5)] = 0.0f;
    _S2842[int(6)] = 0.0f;
    _S2842[int(7)] = 0.0f;
    _S2842[int(8)] = 0.0f;
    _S2842[int(9)] = 0.0f;
    _S2842[int(10)] = 0.0f;
    _S2842[int(11)] = 0.0f;
    _S2842[int(12)] = 0.0f;
    _S2842[int(13)] = 0.0f;
    _S2842[int(14)] = 0.0f;
    _S2842[int(15)] = 0.0f;
    _S2842[int(16)] = 0.0f;
    _S2842[int(17)] = 0.0f;
    _S2842[int(0)] = _S2841;
    FixedArray<float, 18>  _S2843 = {
        _S2820 + _S2842[int(0)], _S2821 + _S2842[int(1)], _S2822 + _S2842[int(2)], _S2823 + _S2842[int(3)], _S2824 + _S2842[int(4)], _S2825 + _S2842[int(5)], _S2826 + _S2842[int(6)], _S2827 + _S2842[int(7)], _S2828 + _S2842[int(8)], _S2829 + _S2842[int(9)], _S2830 + _S2842[int(10)], _S2831 + _S2842[int(11)], _S2832 + _S2842[int(12)], _S2833 + _S2842[int(13)], _S2834 + _S2842[int(14)], _S2835 + _S2842[int(15)], _S2836 + _S2842[int(16)], _S2837 + _S2842[int(17)]
    };
    dpraw_losses_2->primal_0 = dpraw_losses_2->primal_0;
    dpraw_losses_2->differential_0 = _S2843;
    return;
}

inline __device__ void s_bwd_compute_ppisp_no_crf_regularization_loss_0(DiffPair_arrayx3Cfloatx2C18x3E_0 * _S2844, int _S2845, FixedArray<float, 6>  * _S2846, FixedArray<float, 6>  * _S2847)
{
    s_bwd_prop_compute_ppisp_no_crf_regularization_loss_0(_S2844, _S2845, _S2846, _S2847);
    return;
}

inline __device__ void compute_ppisp_no_crf_regularization_loss_vjp(FixedArray<float, 18>  raw_losses_5, int num_cameras_8, FixedArray<float, 6>  loss_weights_8, FixedArray<float, 6>  grad_out_8, FixedArray<float, 18>  * _S2848)
{
    FixedArray<float, 18>  _S2849 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C18x3E_0 dp_raw_losses_2;
    (&dp_raw_losses_2)->primal_0 = raw_losses_5;
    (&dp_raw_losses_2)->differential_0 = _S2849;
    FixedArray<float, 6>  _S2850 = loss_weights_8;
    FixedArray<float, 6>  _S2851 = grad_out_8;
    s_bwd_compute_ppisp_no_crf_regularization_loss_0(&dp_raw_losses_2, num_cameras_8, &_S2850, &_S2851);
    *_S2848 = (&dp_raw_losses_2)->differential_0;
    return;
}

