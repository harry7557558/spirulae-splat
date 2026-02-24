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
    RQSCRFPPISPChannelParams_0 result_2;
    (&result_2)->g0_0 = 0.0f;
    (&result_2)->g1_0 = 0.0f;
    (&result_2)->x0_0 = 0.0f;
    (&result_2)->y0_0 = 0.0f;
    (&result_2)->gc_0 = 0.0f;
    return result_2;
}

struct PPISPParamsRQS_0
{
    float exposure_0;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_0;
    ColorPPISPParams_0 color_params_0;
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  crf_params_0;
};

inline __device__ PPISPParamsRQS_0 PPISPParamsRQS_x24_syn_dzero_0()
{
    PPISPParamsRQS_0 result_3;
    (&result_3)->exposure_0 = 0.0f;
    VignettingChannelParams_0 _S2 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_3)->vignette_params_0[int(0)] = _S2;
    (&result_3)->vignette_params_0[int(1)] = _S2;
    (&result_3)->vignette_params_0[int(2)] = _S2;
    (&result_3)->color_params_0 = ColorPPISPParams_x24_syn_dzero_0();
    RQSCRFPPISPChannelParams_0 _S3 = RQSCRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_3)->crf_params_0[int(0)] = _S3;
    (&result_3)->crf_params_0[int(1)] = _S3;
    (&result_3)->crf_params_0[int(2)] = _S3;
    return result_3;
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
    CRFPPISPChannelParams_0 result_4;
    (&result_4)->toe_0 = 0.0f;
    (&result_4)->shoulder_0 = 0.0f;
    (&result_4)->gamma_0 = 0.0f;
    (&result_4)->center_0 = 0.0f;
    return result_4;
}

struct PPISPParams_0
{
    float exposure_1;
    FixedArray<VignettingChannelParams_0, 3>  vignette_params_1;
    ColorPPISPParams_0 color_params_1;
    FixedArray<CRFPPISPChannelParams_0, 3>  crf_params_1;
};

inline __device__ PPISPParams_0 PPISPParams_x24_syn_dzero_0()
{
    PPISPParams_0 result_5;
    (&result_5)->exposure_1 = 0.0f;
    VignettingChannelParams_0 _S4 = VignettingChannelParams_x24_syn_dzero_0();
    (&result_5)->vignette_params_1[int(0)] = _S4;
    (&result_5)->vignette_params_1[int(1)] = _S4;
    (&result_5)->vignette_params_1[int(2)] = _S4;
    (&result_5)->color_params_1 = ColorPPISPParams_x24_syn_dzero_0();
    CRFPPISPChannelParams_0 _S5 = CRFPPISPChannelParams_x24_syn_dzero_0();
    (&result_5)->crf_params_1[int(0)] = _S5;
    (&result_5)->crf_params_1[int(1)] = _S5;
    (&result_5)->crf_params_1[int(2)] = _S5;
    return result_5;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_0, float dOut_0)
{
    float _S6 = (F32_exp2(((*dpx_0).primal_0))) * 0.69314718246459961f * dOut_0;
    dpx_0->primal_0 = (*dpx_0).primal_0;
    dpx_0->differential_0 = _S6;
    return;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_0, float dOut_1)
{
    DiffPair_float_0 _S7 = *dpx_1;
    float _S8;
    if(((*dpx_1).primal_0) > ((*dpy_0).primal_0))
    {
        _S8 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_0).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S7.primal_0;
    dpx_1->differential_0 = _S8;
    DiffPair_float_0 _S9 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S7.primal_0))
    {
        _S8 = dOut_1;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_1).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_1;
        }
    }
    dpy_0->primal_0 = _S9.primal_0;
    dpy_0->differential_0 = _S8;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_2)
{
    DiffPair_float_0 _S10 = *dpx_2;
    bool _S11;
    if(((*dpx_2).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S11 = ((*dpx_2).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S11 = false;
    }
    float _S12;
    if(_S11)
    {
        _S12 = dOut_2;
    }
    else
    {
        _S12 = 0.0f;
    }
    dpx_2->primal_0 = _S10.primal_0;
    dpx_2->differential_0 = _S12;
    DiffPair_float_0 _S13 = *dpMin_0;
    if((_S10.primal_0) < ((*dpMin_0).primal_0))
    {
        _S12 = dOut_2;
    }
    else
    {
        _S12 = 0.0f;
    }
    dpMin_0->primal_0 = _S13.primal_0;
    dpMin_0->differential_0 = _S12;
    DiffPair_float_0 _S14 = *dpMax_0;
    if(((*dpx_2).primal_0) > ((*dpMax_0).primal_0))
    {
        _S12 = dOut_2;
    }
    else
    {
        _S12 = 0.0f;
    }
    dpMax_0->primal_0 = _S14.primal_0;
    dpMax_0->differential_0 = _S12;
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
    float _S15 = (*left_0).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 2, 2>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_3.x;
    float sum_0 = _S15 + (*left_0).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_3.y;
    float2  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_0;
    float _S16 = (*left_0).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_3.x;
    float sum_1 = _S16 + (*left_0).primal_0.rows[int(1)].y * dOut_3.y;
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
    float _S17 = (*left_1).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = (*right_1).primal_0.x * dOut_4.x;
    float sum_2 = _S17 + (*left_1).primal_0.rows[int(1)].x * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = (*right_1).primal_0.x * dOut_4.y;
    float sum_3 = sum_2 + (*left_1).primal_0.rows[int(2)].x * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = (*right_1).primal_0.x * dOut_4.z;
    float3  right_d_result_1;
    *&((&right_d_result_1)->x) = sum_3;
    float _S18 = (*left_1).primal_0.rows[int(0)].y * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = (*right_1).primal_0.y * dOut_4.x;
    float sum_4 = _S18 + (*left_1).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = (*right_1).primal_0.y * dOut_4.y;
    float sum_5 = sum_4 + (*left_1).primal_0.rows[int(2)].y * dOut_4.z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = (*right_1).primal_0.y * dOut_4.z;
    *&((&right_d_result_1)->y) = sum_5;
    float _S19 = (*left_1).primal_0.rows[int(0)].z * dOut_4.x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = (*right_1).primal_0.z * dOut_4.x;
    float sum_6 = _S19 + (*left_1).primal_0.rows[int(1)].z * dOut_4.y;
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
    float2  result_6;
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
        *_slang_vector_get_element_ptr(&result_6, i_0) = sum_8;
        i_0 = i_0 + int(1);
    }
    return result_6;
}

inline __device__ float3  mul_1(Matrix<float, 3, 3>  left_3, float3  right_3)
{
    float3  result_7;
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
        *_slang_vector_get_element_ptr(&result_7, i_1) = sum_10;
        i_1 = i_1 + int(1);
    }
    return result_7;
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
    Matrix<float, 3, 3>  result_8;
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
            *_slang_vector_get_element_ptr(((&result_8)->rows + (r_1)), c_0) = sum_12;
            c_0 = c_0 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_8;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_1, float3  dOut_6)
{
    float _S20 = dOut_6.y;
    float _S21 = dOut_6.z;
    float _S22 = dOut_6.x;
    float _S23 = (*a_0).primal_0.z * _S20 + - (*a_0).primal_0.y * _S21;
    float _S24 = - (*a_0).primal_0.z * _S22 + (*a_0).primal_0.x * _S21;
    float _S25 = (*a_0).primal_0.y * _S22 + - (*a_0).primal_0.x * _S20;
    float3  _S26 = make_float3 (- (*b_1).primal_0.z * _S20 + (*b_1).primal_0.y * _S21, (*b_1).primal_0.z * _S22 + - (*b_1).primal_0.x * _S21, - (*b_1).primal_0.y * _S22 + (*b_1).primal_0.x * _S20);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S26;
    float3  _S27 = make_float3 (_S23, _S24, _S25);
    b_1->primal_0 = (*b_1).primal_0;
    b_1->differential_0 = _S27;
    return;
}

inline __device__ float3  cross_0(float3  left_6, float3  right_6)
{
    float _S28 = left_6.y;
    float _S29 = right_6.z;
    float _S30 = left_6.z;
    float _S31 = right_6.y;
    float _S32 = right_6.x;
    float _S33 = left_6.x;
    return make_float3 (_S28 * _S29 - _S30 * _S31, _S30 * _S32 - _S33 * _S29, _S33 * _S31 - _S28 * _S32);
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
    float result_9 = 0.0f;
    for(;;)
    {
        if(i_3 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_10 = result_9 + _slang_vector_get_element(x_1, i_3) * _slang_vector_get_element(y_0, i_3);
        i_3 = i_3 + int(1);
        result_9 = result_10;
    }
    return result_9;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_4, float dOut_8)
{
    float _S34 = _slang_select(((*dpx_4).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_4).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S34;
    return;
}

inline __device__ float3  min_0(float3  x_2, float3  y_1)
{
    float3  result_11;
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
        *_slang_vector_get_element_ptr(&result_11, i_4) = (F32_min((_slang_vector_get_element(x_2, i_4)), (_slang_vector_get_element(y_1, i_4))));
        i_4 = i_4 + int(1);
    }
    return result_11;
}

inline __device__ float3  max_0(float3  x_3, float3  y_2)
{
    float3  result_12;
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
        *_slang_vector_get_element_ptr(&result_12, i_5) = (F32_max((_slang_vector_get_element(x_3, i_5)), (_slang_vector_get_element(y_2, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_12;
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
    float _S35 = (F32_exp(((*dpx_6).primal_0))) * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S35;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_7, float dOut_11)
{
    float _S36 = 1.0f / (*dpx_7).primal_0 * dOut_11;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S36;
    return;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpy_3, DiffPair_float_0 * dps_0, float dOut_12)
{
    float _S37 = (1.0f - (*dps_0).primal_0) * dOut_12;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S37;
    DiffPair_float_0 _S38 = *dpy_3;
    float _S39 = (*dps_0).primal_0 * dOut_12;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S39;
    float _S40 = (_S38.primal_0 - (*dpx_8).primal_0) * dOut_12;
    dps_0->primal_0 = _S38.primal_0;
    dps_0->differential_0 = _S40;
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
        DiffPair_float_0 _S41 = *dpx_9;
        float _S42 = val_0 * (*dpy_4).primal_0 / (*dpx_9).primal_0 * dOut_13;
        dpx_9->primal_0 = (*dpx_9).primal_0;
        dpx_9->differential_0 = _S42;
        float _S43 = val_0 * (F32_log((_S41.primal_0))) * dOut_13;
        dpy_4->primal_0 = (*dpy_4).primal_0;
        dpy_4->differential_0 = _S43;
    }
    return;
}

inline __device__ float3  apply_ppisp(float3  rgb_in_0, float2  pix_coord_0, float2  image_center_0, float2  img_size_0, FixedArray<float, 36>  params_0)
{
    PPISPParams_0 p_0;
    (&p_0)->exposure_1 = params_0[int(0)];
    (&(&p_0)->vignette_params_1[int(0)])->cx_0 = params_0[int(1)];
    (&(&p_0)->vignette_params_1[int(0)])->cy_0 = params_0[int(2)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha0_0 = params_0[int(3)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha1_0 = params_0[int(4)];
    (&(&p_0)->vignette_params_1[int(0)])->alpha2_0 = params_0[int(5)];
    (&(&p_0)->vignette_params_1[int(1)])->cx_0 = params_0[int(6)];
    (&(&p_0)->vignette_params_1[int(1)])->cy_0 = params_0[int(7)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha0_0 = params_0[int(8)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha1_0 = params_0[int(9)];
    (&(&p_0)->vignette_params_1[int(1)])->alpha2_0 = params_0[int(10)];
    (&(&p_0)->vignette_params_1[int(2)])->cx_0 = params_0[int(11)];
    (&(&p_0)->vignette_params_1[int(2)])->cy_0 = params_0[int(12)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha0_0 = params_0[int(13)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha1_0 = params_0[int(14)];
    (&(&p_0)->vignette_params_1[int(2)])->alpha2_0 = params_0[int(15)];
    *&((&(&(&p_0)->color_params_1)->b_0)->x) = params_0[int(16)];
    *&((&(&(&p_0)->color_params_1)->b_0)->y) = params_0[int(17)];
    *&((&(&(&p_0)->color_params_1)->r_0)->x) = params_0[int(18)];
    *&((&(&(&p_0)->color_params_1)->r_0)->y) = params_0[int(19)];
    *&((&(&(&p_0)->color_params_1)->g_0)->x) = params_0[int(20)];
    *&((&(&(&p_0)->color_params_1)->g_0)->y) = params_0[int(21)];
    *&((&(&(&p_0)->color_params_1)->n_0)->x) = params_0[int(22)];
    *&((&(&(&p_0)->color_params_1)->n_0)->y) = params_0[int(23)];
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
    PPISPParams_0 _S44 = p_0;
    float _S45 = (F32_max((img_size_0.x), (img_size_0.y)));
    float _S46 = (pix_coord_0.x - image_center_0.x) / _S45;
    float _S47 = (pix_coord_0.y - image_center_0.y) / _S45;
    float3  rgb_out_0 = rgb_in_0 * make_float3 ((F32_exp2((p_0.exposure_1))));
    float dx_0 = _S46 - p_0.vignette_params_1[int(0)].cx_0;
    float dy_0 = _S47 - p_0.vignette_params_1[int(0)].cy_0;
    float r2_0 = dx_0 * dx_0 + dy_0 * dy_0;
    float r4_0 = r2_0 * r2_0;
    *&((&rgb_out_0)->x) = *&((&rgb_out_0)->x) * clamp_0(p_0.vignette_params_1[int(0)].alpha2_0 * (r4_0 * r2_0) + p_0.vignette_params_1[int(0)].alpha1_0 * r4_0 + p_0.vignette_params_1[int(0)].alpha0_0 * r2_0 + 1.0f, 0.0f, 1.0f);
    float dx_1 = _S46 - p_0.vignette_params_1[int(1)].cx_0;
    float dy_1 = _S47 - p_0.vignette_params_1[int(1)].cy_0;
    float r2_1 = dx_1 * dx_1 + dy_1 * dy_1;
    float r4_1 = r2_1 * r2_1;
    *&((&rgb_out_0)->y) = *&((&rgb_out_0)->y) * clamp_0(p_0.vignette_params_1[int(1)].alpha2_0 * (r4_1 * r2_1) + p_0.vignette_params_1[int(1)].alpha1_0 * r4_1 + p_0.vignette_params_1[int(1)].alpha0_0 * r2_1 + 1.0f, 0.0f, 1.0f);
    float dx_2 = _S46 - p_0.vignette_params_1[int(2)].cx_0;
    float dy_2 = _S47 - p_0.vignette_params_1[int(2)].cy_0;
    float r2_2 = dx_2 * dx_2 + dy_2 * dy_2;
    float r4_2 = r2_2 * r2_2;
    *&((&rgb_out_0)->z) = *&((&rgb_out_0)->z) * clamp_0(p_0.vignette_params_1[int(2)].alpha2_0 * (r4_2 * r2_2) + p_0.vignette_params_1[int(2)].alpha1_0 * r4_2 + p_0.vignette_params_1[int(2)].alpha0_0 * r2_2 + 1.0f, 0.0f, 1.0f);
    float3  _S48 = rgb_out_0;
    float2  bd_0 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_0.color_params_1.b_0);
    float2  rd_0 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_0.color_params_1.r_0);
    float2  gd_0 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_0.color_params_1.g_0);
    float2  nd_0 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_0.color_params_1.n_0);
    float _S49 = 0.3333333432674408f + nd_0.x;
    float _S50 = 0.3333333432674408f + nd_0.y;
    Matrix<float, 3, 3>  T_0 = makeMatrix<float, 3, 3> (bd_0.x, 1.0f + rd_0.x, gd_0.x, bd_0.y, rd_0.y, 1.0f + gd_0.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_0 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S50, 1.0f, 0.0f, - _S49, - _S50, _S49, 0.0f), T_0);
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
    float _S51 = _S48.x;
    float _S52 = _S48.y;
    float intensity_0 = _S51 + _S52 + _S48.z;
    float3  rgi_out_0 = mul_1(H_1, make_float3 (_S51, _S52, intensity_0));
    float3  rgi_out_1 = rgi_out_0 * make_float3 (intensity_0 / ((F32_max((rgi_out_0.z), (0.0f))) + 0.00000999999974738f));
    float _S53 = rgi_out_1.x;
    float _S54 = rgi_out_1.y;
    float3  _S55 = clamp_1(make_float3 (_S53, _S54, rgi_out_1.z - _S53 - _S54), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_1;
    float _S56 = _S55.x;
    float _S57 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(0)].toe_0))))));
    float _S58 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(0)].shoulder_0))))));
    float _S59 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(0)].gamma_0))))));
    float _S60 = 1.0f / (1.0f + (F32_exp((- _S44.crf_params_1[int(0)].center_0))));
    float a_1 = _S58 * _S60 / lerp_0(_S57, _S58, _S60);
    float b_2 = 1.0f - a_1;
    float y_4;
    if(_S56 <= _S60)
    {
        y_4 = a_1 * (F32_pow((_S56 / _S60), (_S57)));
    }
    else
    {
        y_4 = 1.0f - b_2 * (F32_pow(((1.0f - _S56) / (1.0f - _S60)), (_S58)));
    }
    *&((&rgb_out_1)->x) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S59)));
    float _S61 = _S55.y;
    float _S62 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(1)].toe_0))))));
    float _S63 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(1)].shoulder_0))))));
    float _S64 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(1)].gamma_0))))));
    float _S65 = 1.0f / (1.0f + (F32_exp((- _S44.crf_params_1[int(1)].center_0))));
    float a_2 = _S63 * _S65 / lerp_0(_S62, _S63, _S65);
    float b_3 = 1.0f - a_2;
    if(_S61 <= _S65)
    {
        y_4 = a_2 * (F32_pow((_S61 / _S65), (_S62)));
    }
    else
    {
        y_4 = 1.0f - b_3 * (F32_pow(((1.0f - _S61) / (1.0f - _S65)), (_S63)));
    }
    *&((&rgb_out_1)->y) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S64)));
    float _S66 = _S55.z;
    float _S67 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(2)].toe_0))))));
    float _S68 = 0.30000001192092896f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(2)].shoulder_0))))));
    float _S69 = 0.10000000149011612f + (F32_log((1.0f + (F32_exp((_S44.crf_params_1[int(2)].gamma_0))))));
    float _S70 = 1.0f / (1.0f + (F32_exp((- _S44.crf_params_1[int(2)].center_0))));
    float a_3 = _S68 * _S70 / lerp_0(_S67, _S68, _S70);
    float b_4 = 1.0f - a_3;
    if(_S66 <= _S70)
    {
        y_4 = a_3 * (F32_pow((_S66 / _S70), (_S67)));
    }
    else
    {
        y_4 = 1.0f - b_4 * (F32_pow(((1.0f - _S66) / (1.0f - _S70)), (_S68)));
    }
    *&((&rgb_out_1)->z) = (F32_pow(((F32_max((0.0f), (y_4)))), (_S69)));
    return rgb_out_1;
}

inline __device__ float3  apply_ppisp_rqs(float3  rgb_in_1, float2  pix_coord_1, float2  image_center_1, float2  img_size_1, FixedArray<float, 39>  params_1)
{
    PPISPParamsRQS_0 p_1;
    (&p_1)->exposure_0 = params_1[int(0)];
    (&(&p_1)->vignette_params_0[int(0)])->cx_0 = params_1[int(1)];
    (&(&p_1)->vignette_params_0[int(0)])->cy_0 = params_1[int(2)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha0_0 = params_1[int(3)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha1_0 = params_1[int(4)];
    (&(&p_1)->vignette_params_0[int(0)])->alpha2_0 = params_1[int(5)];
    (&(&p_1)->vignette_params_0[int(1)])->cx_0 = params_1[int(6)];
    (&(&p_1)->vignette_params_0[int(1)])->cy_0 = params_1[int(7)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha0_0 = params_1[int(8)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha1_0 = params_1[int(9)];
    (&(&p_1)->vignette_params_0[int(1)])->alpha2_0 = params_1[int(10)];
    (&(&p_1)->vignette_params_0[int(2)])->cx_0 = params_1[int(11)];
    (&(&p_1)->vignette_params_0[int(2)])->cy_0 = params_1[int(12)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha0_0 = params_1[int(13)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha1_0 = params_1[int(14)];
    (&(&p_1)->vignette_params_0[int(2)])->alpha2_0 = params_1[int(15)];
    *&((&(&(&p_1)->color_params_0)->b_0)->x) = params_1[int(16)];
    *&((&(&(&p_1)->color_params_0)->b_0)->y) = params_1[int(17)];
    *&((&(&(&p_1)->color_params_0)->r_0)->x) = params_1[int(18)];
    *&((&(&(&p_1)->color_params_0)->r_0)->y) = params_1[int(19)];
    *&((&(&(&p_1)->color_params_0)->g_0)->x) = params_1[int(20)];
    *&((&(&(&p_1)->color_params_0)->g_0)->y) = params_1[int(21)];
    *&((&(&(&p_1)->color_params_0)->n_0)->x) = params_1[int(22)];
    *&((&(&(&p_1)->color_params_0)->n_0)->y) = params_1[int(23)];
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
    PPISPParamsRQS_0 _S71 = p_1;
    float _S72 = (F32_max((img_size_1.x), (img_size_1.y)));
    float _S73 = (pix_coord_1.x - image_center_1.x) / _S72;
    float _S74 = (pix_coord_1.y - image_center_1.y) / _S72;
    float3  rgb_out_2 = rgb_in_1 * make_float3 ((F32_exp2((p_1.exposure_0))));
    float dx_3 = _S73 - p_1.vignette_params_0[int(0)].cx_0;
    float dy_3 = _S74 - p_1.vignette_params_0[int(0)].cy_0;
    float r2_4 = dx_3 * dx_3 + dy_3 * dy_3;
    float r4_3 = r2_4 * r2_4;
    *&((&rgb_out_2)->x) = *&((&rgb_out_2)->x) * clamp_0(p_1.vignette_params_0[int(0)].alpha2_0 * (r4_3 * r2_4) + p_1.vignette_params_0[int(0)].alpha1_0 * r4_3 + p_1.vignette_params_0[int(0)].alpha0_0 * r2_4 + 1.0f, 0.0f, 1.0f);
    float dx_4 = _S73 - p_1.vignette_params_0[int(1)].cx_0;
    float dy_4 = _S74 - p_1.vignette_params_0[int(1)].cy_0;
    float r2_5 = dx_4 * dx_4 + dy_4 * dy_4;
    float r4_4 = r2_5 * r2_5;
    *&((&rgb_out_2)->y) = *&((&rgb_out_2)->y) * clamp_0(p_1.vignette_params_0[int(1)].alpha2_0 * (r4_4 * r2_5) + p_1.vignette_params_0[int(1)].alpha1_0 * r4_4 + p_1.vignette_params_0[int(1)].alpha0_0 * r2_5 + 1.0f, 0.0f, 1.0f);
    float dx_5 = _S73 - p_1.vignette_params_0[int(2)].cx_0;
    float dy_5 = _S74 - p_1.vignette_params_0[int(2)].cy_0;
    float r2_6 = dx_5 * dx_5 + dy_5 * dy_5;
    float r4_5 = r2_6 * r2_6;
    *&((&rgb_out_2)->z) = *&((&rgb_out_2)->z) * clamp_0(p_1.vignette_params_0[int(2)].alpha2_0 * (r4_5 * r2_6) + p_1.vignette_params_0[int(2)].alpha1_0 * r4_5 + p_1.vignette_params_0[int(2)].alpha0_0 * r2_6 + 1.0f, 0.0f, 1.0f);
    float3  _S75 = rgb_out_2;
    float2  bd_1 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_1.color_params_0.b_0);
    float2  rd_1 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_1.color_params_0.r_0);
    float2  gd_1 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_1.color_params_0.g_0);
    float2  nd_1 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_1.color_params_0.n_0);
    float _S76 = 0.3333333432674408f + nd_1.x;
    float _S77 = 0.3333333432674408f + nd_1.y;
    Matrix<float, 3, 3>  T_1 = makeMatrix<float, 3, 3> (bd_1.x, 1.0f + rd_1.x, gd_1.x, bd_1.y, rd_1.y, 1.0f + gd_1.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  M_1 = mul_3(makeMatrix<float, 3, 3> (0.0f, -1.0f, _S77, 1.0f, 0.0f, - _S76, - _S77, _S76, 0.0f), T_1);
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
    float _S78 = _S75.x;
    float _S79 = _S75.y;
    float intensity_1 = _S78 + _S79 + _S75.z;
    float3  rgi_out_2 = mul_1(H_3, make_float3 (_S78, _S79, intensity_1));
    float3  rgi_out_3 = rgi_out_2 * make_float3 (intensity_1 / ((F32_max((rgi_out_2.z), (0.0f))) + 0.00000999999974738f));
    float _S80 = rgi_out_3.x;
    float _S81 = rgi_out_3.y;
    float3  _S82 = clamp_1(make_float3 (_S80, _S81, rgi_out_3.z - _S80 - _S81), make_float3 (0.0f), make_float3 (1.0f));
    float3  rgb_out_3;
    float _S83 = _S82.x;
    float g0_1 = (F32_exp((_S71.crf_params_0[int(0)].g0_0)));
    float g1_1 = (F32_exp((_S71.crf_params_0[int(0)].g1_0)));
    float x0_1 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(0)].x0_0))));
    float y0_1 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(0)].y0_0))));
    float gc_1 = (F32_exp((_S71.crf_params_0[int(0)].gc_0)));
    float y_5;
    if(_S83 < x0_1)
    {
        float s0_0 = y0_1 / x0_1;
        float t0_0 = _S83 / x0_1;
        float _S84 = 1.0f - t0_0;
        y_5 = y0_1 * (s0_0 * t0_0 * t0_0 + g0_1 * t0_0 * _S84) / (s0_0 + (g0_1 + gc_1 - 2.0f * s0_0) * t0_0 * _S84);
    }
    else
    {
        float _S85 = 1.0f - y0_1;
        float _S86 = 1.0f - x0_1;
        float s1_0 = _S85 / _S86;
        float t1_0 = (_S83 - x0_1) / _S86;
        float _S87 = 1.0f - t1_0;
        y_5 = y0_1 + _S85 * (s1_0 * t1_0 * t1_0 + gc_1 * t1_0 * _S87) / (s1_0 + (gc_1 + g1_1 - 2.0f * s1_0) * t1_0 * _S87);
    }
    *&((&rgb_out_3)->x) = y_5;
    float _S88 = _S82.y;
    float g0_2 = (F32_exp((_S71.crf_params_0[int(1)].g0_0)));
    float g1_2 = (F32_exp((_S71.crf_params_0[int(1)].g1_0)));
    float x0_2 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(1)].x0_0))));
    float y0_2 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(1)].y0_0))));
    float gc_2 = (F32_exp((_S71.crf_params_0[int(1)].gc_0)));
    if(_S88 < x0_2)
    {
        float s0_1 = y0_2 / x0_2;
        float t0_1 = _S88 / x0_2;
        float _S89 = 1.0f - t0_1;
        y_5 = y0_2 * (s0_1 * t0_1 * t0_1 + g0_2 * t0_1 * _S89) / (s0_1 + (g0_2 + gc_2 - 2.0f * s0_1) * t0_1 * _S89);
    }
    else
    {
        float _S90 = 1.0f - y0_2;
        float _S91 = 1.0f - x0_2;
        float s1_1 = _S90 / _S91;
        float t1_1 = (_S88 - x0_2) / _S91;
        float _S92 = 1.0f - t1_1;
        y_5 = y0_2 + _S90 * (s1_1 * t1_1 * t1_1 + gc_2 * t1_1 * _S92) / (s1_1 + (gc_2 + g1_2 - 2.0f * s1_1) * t1_1 * _S92);
    }
    *&((&rgb_out_3)->y) = y_5;
    float _S93 = _S82.z;
    float g0_3 = (F32_exp((_S71.crf_params_0[int(2)].g0_0)));
    float g1_3 = (F32_exp((_S71.crf_params_0[int(2)].g1_0)));
    float x0_3 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(2)].x0_0))));
    float y0_3 = 1.0f / (1.0f + (F32_exp((- _S71.crf_params_0[int(2)].y0_0))));
    float gc_3 = (F32_exp((_S71.crf_params_0[int(2)].gc_0)));
    if(_S93 < x0_3)
    {
        float s0_2 = y0_3 / x0_3;
        float t0_2 = _S93 / x0_3;
        float _S94 = 1.0f - t0_2;
        y_5 = y0_3 * (s0_2 * t0_2 * t0_2 + g0_3 * t0_2 * _S94) / (s0_2 + (g0_3 + gc_3 - 2.0f * s0_2) * t0_2 * _S94);
    }
    else
    {
        float _S95 = 1.0f - y0_3;
        float _S96 = 1.0f - x0_3;
        float s1_2 = _S95 / _S96;
        float t1_2 = (_S93 - x0_3) / _S96;
        float _S97 = 1.0f - t1_2;
        y_5 = y0_3 + _S95 * (s1_2 * t1_2 * t1_2 + gc_3 * t1_2 * _S97) / (s1_2 + (gc_3 + g1_3 - 2.0f * s1_2) * t1_2 * _S97);
    }
    *&((&rgb_out_3)->z) = y_5;
    return rgb_out_3;
}

struct DiffPair_arrayx3Cfloatx2C36x3E_0
{
    FixedArray<float, 36>  primal_0;
    FixedArray<float, 36>  differential_0;
};

inline __device__ float s_primal_ctx_exp2_0(float _S98)
{
    return (F32_exp2((_S98)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S99, float _S100, float _S101)
{
    return clamp_0(_S99, _S100, _S101);
}

inline __device__ float2  s_primal_ctx_mul_0(Matrix<float, 2, 2>  _S102, float2  _S103)
{
    return mul_0(_S102, _S103);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S104, Matrix<float, 3, 3>  _S105)
{
    return mul_3(_S104, _S105);
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S106, float3  _S107)
{
    return cross_0(_S106, _S107);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S108, float3  _S109)
{
    return dot_0(_S108, _S109);
}

inline __device__ float s_primal_ctx_abs_0(float _S110)
{
    return (F32_abs((_S110)));
}

inline __device__ float3  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S111, float3  _S112)
{
    return mul_1(_S111, _S112);
}

inline __device__ float3  s_primal_ctx_clamp_1(float3  _S113, float3  _S114, float3  _S115)
{
    return clamp_1(_S113, _S114, _S115);
}

inline __device__ float s_primal_ctx_exp_0(float _S116)
{
    return (F32_exp((_S116)));
}

inline __device__ float s_primal_ctx_log_0(float _S117)
{
    return (F32_log((_S117)));
}

inline __device__ float s_primal_ctx_lerp_0(float _S118, float _S119, float _S120)
{
    return lerp_0(_S118, _S119, _S120);
}

inline __device__ float s_primal_ctx_pow_0(float _S121, float _S122)
{
    return (F32_pow((_S121), (_S122)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S123, DiffPair_float_0 * _S124, float _S125)
{
    _d_pow_0(_S123, _S124, _S125);
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S126, DiffPair_float_0 * _S127, DiffPair_float_0 * _S128, float _S129)
{
    _d_lerp_0(_S126, _S127, _S128, _S129);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S130, float _S131)
{
    _d_exp_0(_S130, _S131);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S132, float _S133)
{
    _d_log_0(_S132, _S133);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S134, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S135, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S136, float3  _S137)
{
    _d_clamp_vector_0(_S134, _S135, _S136, _S137);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S138, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S139, float3  _S140)
{
    _d_mul_1(_S138, _S139, _S140);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S141, float _S142)
{
    _d_abs_0(_S141, _S142);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S143, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S144, Matrix<float, 3, 3>  _S145)
{
    mul_2(_S143, _S144, _S145);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S146, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S147, float3  _S148)
{
    _d_cross_0(_S146, _S147, _S148);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S149, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S150, float _S151)
{
    _d_dot_0(_S149, _S150, _S151);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 * _S152, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S153, float2  _S154)
{
    _d_mul_0(_S152, _S153, _S154);
    return;
}

inline __device__ void s_bwd_prop_clamp_1(DiffPair_float_0 * _S155, DiffPair_float_0 * _S156, DiffPair_float_0 * _S157, float _S158)
{
    _d_clamp_0(_S155, _S156, _S157, _S158);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S159, float _S160)
{
    _d_exp2_0(_S159, _S160);
    return;
}

inline __device__ void s_bwd_prop_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_0, float2  pix_coord_2, float2  image_center_2, float2  img_size_2, DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_0, float3  _s_dOut_0)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S161 = *dprgb_in_0;
    float3  _S162 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S163 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S164 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S165 = {
        _S164, _S164, _S164
    };
    float2  _S166 = make_float2 (0.0f);
    ColorPPISPParams_0 _S167 = { _S166, _S166, _S166, _S166 };
    CRFPPISPChannelParams_0 _S168 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S169 = {
        _S168, _S168, _S168
    };
    PPISPParams_0 _S170;
    (&_S170)->exposure_1 = dpparams_0->primal_0[int(0)];
    (&_S170)->vignette_params_1 = _S165;
    (&_S170)->color_params_1 = _S167;
    (&_S170)->crf_params_1 = _S169;
    (&(&_S170)->vignette_params_1[int(0)])->cx_0 = dpparams_0->primal_0[int(1)];
    (&(&_S170)->vignette_params_1[int(0)])->cy_0 = dpparams_0->primal_0[int(2)];
    float _S171 = dpparams_0->primal_0[int(3)];
    (&(&_S170)->vignette_params_1[int(0)])->alpha0_0 = dpparams_0->primal_0[int(3)];
    float _S172 = dpparams_0->primal_0[int(4)];
    (&(&_S170)->vignette_params_1[int(0)])->alpha1_0 = dpparams_0->primal_0[int(4)];
    float _S173 = dpparams_0->primal_0[int(5)];
    (&(&_S170)->vignette_params_1[int(0)])->alpha2_0 = dpparams_0->primal_0[int(5)];
    (&(&_S170)->vignette_params_1[int(1)])->cx_0 = dpparams_0->primal_0[int(6)];
    (&(&_S170)->vignette_params_1[int(1)])->cy_0 = dpparams_0->primal_0[int(7)];
    float _S174 = dpparams_0->primal_0[int(8)];
    (&(&_S170)->vignette_params_1[int(1)])->alpha0_0 = dpparams_0->primal_0[int(8)];
    float _S175 = dpparams_0->primal_0[int(9)];
    (&(&_S170)->vignette_params_1[int(1)])->alpha1_0 = dpparams_0->primal_0[int(9)];
    float _S176 = dpparams_0->primal_0[int(10)];
    (&(&_S170)->vignette_params_1[int(1)])->alpha2_0 = dpparams_0->primal_0[int(10)];
    (&(&_S170)->vignette_params_1[int(2)])->cx_0 = dpparams_0->primal_0[int(11)];
    (&(&_S170)->vignette_params_1[int(2)])->cy_0 = dpparams_0->primal_0[int(12)];
    float _S177 = dpparams_0->primal_0[int(13)];
    (&(&_S170)->vignette_params_1[int(2)])->alpha0_0 = dpparams_0->primal_0[int(13)];
    float _S178 = dpparams_0->primal_0[int(14)];
    (&(&_S170)->vignette_params_1[int(2)])->alpha1_0 = dpparams_0->primal_0[int(14)];
    float _S179 = dpparams_0->primal_0[int(15)];
    (&(&_S170)->vignette_params_1[int(2)])->alpha2_0 = dpparams_0->primal_0[int(15)];
    *&((&(&(&_S170)->color_params_1)->b_0)->x) = dpparams_0->primal_0[int(16)];
    *&((&(&(&_S170)->color_params_1)->b_0)->y) = dpparams_0->primal_0[int(17)];
    *&((&(&(&_S170)->color_params_1)->r_0)->x) = dpparams_0->primal_0[int(18)];
    *&((&(&(&_S170)->color_params_1)->r_0)->y) = dpparams_0->primal_0[int(19)];
    *&((&(&(&_S170)->color_params_1)->g_0)->x) = dpparams_0->primal_0[int(20)];
    *&((&(&(&_S170)->color_params_1)->g_0)->y) = dpparams_0->primal_0[int(21)];
    *&((&(&(&_S170)->color_params_1)->n_0)->x) = dpparams_0->primal_0[int(22)];
    *&((&(&(&_S170)->color_params_1)->n_0)->y) = dpparams_0->primal_0[int(23)];
    float _S180 = dpparams_0->primal_0[int(24)];
    (&(&_S170)->crf_params_1[int(0)])->toe_0 = dpparams_0->primal_0[int(24)];
    float _S181 = dpparams_0->primal_0[int(25)];
    (&(&_S170)->crf_params_1[int(0)])->shoulder_0 = dpparams_0->primal_0[int(25)];
    float _S182 = dpparams_0->primal_0[int(26)];
    (&(&_S170)->crf_params_1[int(0)])->gamma_0 = dpparams_0->primal_0[int(26)];
    float _S183 = dpparams_0->primal_0[int(27)];
    (&(&_S170)->crf_params_1[int(0)])->center_0 = dpparams_0->primal_0[int(27)];
    float _S184 = dpparams_0->primal_0[int(28)];
    (&(&_S170)->crf_params_1[int(1)])->toe_0 = dpparams_0->primal_0[int(28)];
    float _S185 = dpparams_0->primal_0[int(29)];
    (&(&_S170)->crf_params_1[int(1)])->shoulder_0 = dpparams_0->primal_0[int(29)];
    float _S186 = dpparams_0->primal_0[int(30)];
    (&(&_S170)->crf_params_1[int(1)])->gamma_0 = dpparams_0->primal_0[int(30)];
    float _S187 = dpparams_0->primal_0[int(31)];
    (&(&_S170)->crf_params_1[int(1)])->center_0 = dpparams_0->primal_0[int(31)];
    float _S188 = dpparams_0->primal_0[int(32)];
    (&(&_S170)->crf_params_1[int(2)])->toe_0 = dpparams_0->primal_0[int(32)];
    float _S189 = dpparams_0->primal_0[int(33)];
    (&(&_S170)->crf_params_1[int(2)])->shoulder_0 = dpparams_0->primal_0[int(33)];
    float _S190 = dpparams_0->primal_0[int(34)];
    (&(&_S170)->crf_params_1[int(2)])->gamma_0 = dpparams_0->primal_0[int(34)];
    float _S191 = dpparams_0->primal_0[int(35)];
    (&(&_S170)->crf_params_1[int(2)])->center_0 = dpparams_0->primal_0[int(35)];
    PPISPParams_0 _S192 = _S170;
    float _S193 = s_primal_ctx_exp2_0(_S170.exposure_1);
    float3  _S194 = make_float3 (_S193);
    float3  rgb_out_4 = (*dprgb_in_0).primal_0 * make_float3 (_S193);
    float _S195 = (F32_max((img_size_2.x), (img_size_2.y)));
    float _S196 = (pix_coord_2.x - image_center_2.x) / _S195;
    float _S197 = (pix_coord_2.y - image_center_2.y) / _S195;
    float dx_6 = _S196 - dpparams_0->primal_0[int(1)];
    float dy_6 = _S197 - dpparams_0->primal_0[int(2)];
    float r2_8 = dx_6 * dx_6 + dy_6 * dy_6;
    float r4_6 = r2_8 * r2_8;
    float r6_0 = r4_6 * r2_8;
    float falloff_0 = dpparams_0->primal_0[int(5)] * r6_0 + dpparams_0->primal_0[int(4)] * r4_6 + dpparams_0->primal_0[int(3)] * r2_8 + 1.0f;
    float _S198 = s_primal_ctx_clamp_0(falloff_0, 0.0f, 1.0f);
    float _S199 = rgb_out_4.x * _S198;
    float3  _S200 = rgb_out_4;
    *&((&_S200)->x) = _S199;
    float dx_7 = _S196 - dpparams_0->primal_0[int(6)];
    float dy_7 = _S197 - dpparams_0->primal_0[int(7)];
    float r2_9 = dx_7 * dx_7 + dy_7 * dy_7;
    float r4_7 = r2_9 * r2_9;
    float r6_1 = r4_7 * r2_9;
    float falloff_1 = dpparams_0->primal_0[int(10)] * r6_1 + dpparams_0->primal_0[int(9)] * r4_7 + dpparams_0->primal_0[int(8)] * r2_9 + 1.0f;
    float _S201 = s_primal_ctx_clamp_0(falloff_1, 0.0f, 1.0f);
    *&((&_S200)->y) = rgb_out_4.y * _S201;
    float dx_8 = _S196 - dpparams_0->primal_0[int(11)];
    float dy_8 = _S197 - dpparams_0->primal_0[int(12)];
    float r2_10 = dx_8 * dx_8 + dy_8 * dy_8;
    float r4_8 = r2_10 * r2_10;
    float r6_2 = r4_8 * r2_10;
    float falloff_2 = dpparams_0->primal_0[int(15)] * r6_2 + dpparams_0->primal_0[int(14)] * r4_8 + dpparams_0->primal_0[int(13)] * r2_10 + 1.0f;
    float _S202 = s_primal_ctx_clamp_0(falloff_2, 0.0f, 1.0f);
    *&((&_S200)->z) = rgb_out_4.z * _S202;
    PPISPParams_0 _S203 = _S170;
    float2  _S204 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S170.color_params_1.b_0);
    float2  _S205 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S170.color_params_1.r_0);
    float2  _S206 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S170.color_params_1.g_0);
    float2  _S207 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S170.color_params_1.n_0);
    float _S208 = 0.3333333432674408f + _S207.x;
    float _S209 = 0.3333333432674408f + _S207.y;
    Matrix<float, 3, 3>  T_2 = makeMatrix<float, 3, 3> (_S204.x, 1.0f + _S205.x, _S206.x, _S204.y, _S205.y, 1.0f + _S206.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_0 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S209, 1.0f, 0.0f, - _S208, - _S209, _S208, 0.0f);
    Matrix<float, 3, 3>  _S210 = s_primal_ctx_mul_1(skew_0, T_2);
    float3  r0_2 = make_float3 (_S210.rows[int(0)].x, _S210.rows[int(0)].y, _S210.rows[int(0)].z);
    float3  r1_2 = make_float3 (_S210.rows[int(1)].x, _S210.rows[int(1)].y, _S210.rows[int(1)].z);
    float3  r2_11 = make_float3 (_S210.rows[int(2)].x, _S210.rows[int(2)].y, _S210.rows[int(2)].z);
    float3  _S211 = s_primal_ctx_cross_0(r0_2, r1_2);
    bool _S212 = (s_primal_ctx_dot_0(_S211, _S211)) < 9.99999968265522539e-21f;
    float3  lambda_v_6;
    float3  _S213;
    bool _S214;
    if(_S212)
    {
        float3  _S215 = s_primal_ctx_cross_0(r0_2, r2_11);
        bool _S216 = (s_primal_ctx_dot_0(_S215, _S215)) < 9.99999968265522539e-21f;
        if(_S216)
        {
            lambda_v_6 = s_primal_ctx_cross_0(r1_2, r2_11);
        }
        else
        {
            lambda_v_6 = _S215;
        }
        _S214 = _S216;
        _S213 = _S215;
    }
    else
    {
        lambda_v_6 = _S211;
        _S214 = false;
        _S213 = _S162;
    }
    Matrix<float, 3, 3>  S_inv_0 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_0 = makeMatrix<float, 3, 3> (lambda_v_6.x, 0.0f, 0.0f, 0.0f, lambda_v_6.y, 0.0f, 0.0f, 0.0f, lambda_v_6.z);
    Matrix<float, 3, 3>  _S217 = s_primal_ctx_mul_1(T_2, D_0);
    Matrix<float, 3, 3>  _S218 = s_primal_ctx_mul_1(_S217, S_inv_0);
    bool _S219 = (s_primal_ctx_abs_0(_S218.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_4;
    Matrix<float, 3, 3>  _S220;
    float _S221;
    if(_S219)
    {
        float inv_s_0 = 1.0f / _S218.rows[int(2)].z;
        Matrix<float, 3, 3>  _S222 = makeMatrix<float, 3, 3> (inv_s_0);
        float _S223 = _S218.rows[int(2)].z * _S218.rows[int(2)].z;
        H_4 = _S218 * makeMatrix<float, 3, 3> (inv_s_0);
        _S220 = _S222;
        _S221 = _S223;
    }
    else
    {
        H_4 = _S218;
        _S220 = _S163;
        _S221 = 0.0f;
    }
    float _S224 = _S200.x;
    float _S225 = _S200.y;
    float intensity_2 = _S224 + _S225 + _S200.z;
    float3  rgi_in_0 = make_float3 (_S224, _S225, intensity_2);
    float3  _S226 = s_primal_ctx_mul_2(H_4, rgi_in_0);
    float _S227 = _S226.z;
    float _S228 = (F32_max((_S227), (0.0f))) + 0.00000999999974738f;
    float norm_factor_0 = intensity_2 / _S228;
    float3  _S229 = make_float3 (norm_factor_0);
    float _S230 = _S228 * _S228;
    float3  rgi_out_4 = _S226 * make_float3 (norm_factor_0);
    float _S231 = rgi_out_4.x;
    float _S232 = rgi_out_4.y;
    float3  _S233 = make_float3 (_S231, _S232, rgi_out_4.z - _S231 - _S232);
    float3  _S234 = make_float3 (0.0f);
    float3  _S235 = make_float3 (1.0f);
    float3  _S236 = s_primal_ctx_clamp_1(_S233, _S234, _S235);
    float _S237 = _S236.x;
    float _S238 = 1.0f + s_primal_ctx_exp_0(_S180);
    float _S239 = 0.30000001192092896f + s_primal_ctx_log_0(_S238);
    float _S240 = 1.0f + s_primal_ctx_exp_0(_S181);
    float _S241 = 0.30000001192092896f + s_primal_ctx_log_0(_S240);
    float _S242 = 1.0f + s_primal_ctx_exp_0(_S182);
    float _S243 = 0.10000000149011612f + s_primal_ctx_log_0(_S242);
    float _S244 = - _S183;
    float _S245 = 1.0f + s_primal_ctx_exp_0(_S244);
    float _S246 = 1.0f / _S245;
    float _S247 = _S245 * _S245;
    float _S248 = s_primal_ctx_lerp_0(_S239, _S241, _S246);
    float _S249 = _S241 * _S246;
    float a_4 = _S249 / _S248;
    float _S250 = _S248 * _S248;
    float b_5 = 1.0f - a_4;
    bool _S251 = _S237 <= _S246;
    float y_6;
    float _S252;
    float _S253;
    float _S254;
    float _S255;
    float _S256;
    float _S257;
    float _S258;
    float _S259;
    if(_S251)
    {
        float _S260 = _S237 / _S246;
        float _S261 = _S246 * _S246;
        float _S262 = s_primal_ctx_pow_0(_S260, _S239);
        y_6 = a_4 * _S262;
        _S252 = _S262;
        _S253 = _S260;
        _S254 = _S261;
        _S255 = 0.0f;
        _S256 = 0.0f;
        _S257 = 0.0f;
        _S258 = 0.0f;
        _S259 = 0.0f;
    }
    else
    {
        float _S263 = 1.0f - _S237;
        float _S264 = 1.0f - _S246;
        float _S265 = _S263 / _S264;
        float _S266 = _S264 * _S264;
        float _S267 = s_primal_ctx_pow_0(_S265, _S241);
        y_6 = 1.0f - b_5 * _S267;
        _S252 = 0.0f;
        _S253 = 0.0f;
        _S254 = 0.0f;
        _S255 = _S267;
        _S256 = _S265;
        _S257 = _S266;
        _S258 = _S263;
        _S259 = _S264;
    }
    float _S268 = (F32_max((0.0f), (y_6)));
    float _S269 = _S236.y;
    float _S270 = 1.0f + s_primal_ctx_exp_0(_S184);
    float _S271 = 0.30000001192092896f + s_primal_ctx_log_0(_S270);
    float _S272 = 1.0f + s_primal_ctx_exp_0(_S185);
    float _S273 = 0.30000001192092896f + s_primal_ctx_log_0(_S272);
    float _S274 = 1.0f + s_primal_ctx_exp_0(_S186);
    float _S275 = 0.10000000149011612f + s_primal_ctx_log_0(_S274);
    float _S276 = - _S187;
    float _S277 = 1.0f + s_primal_ctx_exp_0(_S276);
    float _S278 = 1.0f / _S277;
    float _S279 = _S277 * _S277;
    float _S280 = s_primal_ctx_lerp_0(_S271, _S273, _S278);
    float _S281 = _S273 * _S278;
    float a_5 = _S281 / _S280;
    float _S282 = _S280 * _S280;
    float b_6 = 1.0f - a_5;
    bool _S283 = _S269 <= _S278;
    float y_7;
    float _S284;
    float _S285;
    float _S286;
    float _S287;
    float _S288;
    float _S289;
    float _S290;
    float _S291;
    if(_S283)
    {
        float _S292 = _S269 / _S278;
        float _S293 = _S278 * _S278;
        float _S294 = s_primal_ctx_pow_0(_S292, _S271);
        y_7 = a_5 * _S294;
        _S284 = _S294;
        _S285 = _S292;
        _S286 = _S293;
        _S287 = 0.0f;
        _S288 = 0.0f;
        _S289 = 0.0f;
        _S290 = 0.0f;
        _S291 = 0.0f;
    }
    else
    {
        float _S295 = 1.0f - _S269;
        float _S296 = 1.0f - _S278;
        float _S297 = _S295 / _S296;
        float _S298 = _S296 * _S296;
        float _S299 = s_primal_ctx_pow_0(_S297, _S273);
        y_7 = 1.0f - b_6 * _S299;
        _S284 = 0.0f;
        _S285 = 0.0f;
        _S286 = 0.0f;
        _S287 = _S299;
        _S288 = _S297;
        _S289 = _S298;
        _S290 = _S295;
        _S291 = _S296;
    }
    float _S300 = (F32_max((0.0f), (y_7)));
    float _S301 = _S236.z;
    float _S302 = 1.0f + s_primal_ctx_exp_0(_S188);
    float _S303 = 0.30000001192092896f + s_primal_ctx_log_0(_S302);
    float _S304 = 1.0f + s_primal_ctx_exp_0(_S189);
    float _S305 = 0.30000001192092896f + s_primal_ctx_log_0(_S304);
    float _S306 = 1.0f + s_primal_ctx_exp_0(_S190);
    float _S307 = 0.10000000149011612f + s_primal_ctx_log_0(_S306);
    float _S308 = - _S191;
    float _S309 = 1.0f + s_primal_ctx_exp_0(_S308);
    float _S310 = 1.0f / _S309;
    float _S311 = _S309 * _S309;
    float _S312 = s_primal_ctx_lerp_0(_S303, _S305, _S310);
    float _S313 = _S305 * _S310;
    float a_6 = _S313 / _S312;
    float _S314 = _S312 * _S312;
    float b_7 = 1.0f - a_6;
    bool _S315 = _S301 <= _S310;
    float y_8;
    float _S316;
    float _S317;
    float _S318;
    float _S319;
    float _S320;
    float _S321;
    float _S322;
    float _S323;
    if(_S315)
    {
        float _S324 = _S301 / _S310;
        float _S325 = _S310 * _S310;
        float _S326 = s_primal_ctx_pow_0(_S324, _S303);
        y_8 = a_6 * _S326;
        _S316 = _S326;
        _S317 = _S324;
        _S318 = _S325;
        _S319 = 0.0f;
        _S320 = 0.0f;
        _S321 = 0.0f;
        _S322 = 0.0f;
        _S323 = 0.0f;
    }
    else
    {
        float _S327 = 1.0f - _S301;
        float _S328 = 1.0f - _S310;
        float _S329 = _S327 / _S328;
        float _S330 = _S328 * _S328;
        float _S331 = s_primal_ctx_pow_0(_S329, _S305);
        y_8 = 1.0f - b_7 * _S331;
        _S316 = 0.0f;
        _S317 = 0.0f;
        _S318 = 0.0f;
        _S319 = _S331;
        _S320 = _S329;
        _S321 = _S330;
        _S322 = _S327;
        _S323 = _S328;
    }
    float _S332 = (F32_max((0.0f), (y_8)));
    DiffPair_float_0 _S333;
    (&_S333)->primal_0 = _S332;
    (&_S333)->differential_0 = 0.0f;
    DiffPair_float_0 _S334;
    (&_S334)->primal_0 = _S307;
    (&_S334)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S333, &_S334, _s_dOut_0.z);
    DiffPair_float_0 _S335 = _S334;
    DiffPair_float_0 _S336;
    (&_S336)->primal_0 = 0.0f;
    (&_S336)->differential_0 = 0.0f;
    DiffPair_float_0 _S337;
    (&_S337)->primal_0 = y_8;
    (&_S337)->differential_0 = 0.0f;
    _d_max_0(&_S336, &_S337, _S333.differential_0);
    DiffPair_float_0 _S338 = _S337;
    if(_S315)
    {
        float _S339 = a_6 * _S338.differential_0;
        float _S340 = _S316 * _S338.differential_0;
        DiffPair_float_0 _S341;
        (&_S341)->primal_0 = _S317;
        (&_S341)->differential_0 = 0.0f;
        DiffPair_float_0 _S342;
        (&_S342)->primal_0 = _S303;
        (&_S342)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S341, &_S342, _S339);
        float _S343 = _S341.differential_0 / _S318;
        float _S344 = _S301 * - _S343;
        float _S345 = _S310 * _S343;
        y_8 = 0.0f;
        _S316 = _S340;
        _S317 = _S344;
        _S318 = 0.0f;
        _S319 = _S342.differential_0;
        _S320 = _S345;
    }
    else
    {
        float _S346 = - _S338.differential_0;
        float _S347 = b_7 * _S346;
        float _S348 = _S319 * _S346;
        DiffPair_float_0 _S349;
        (&_S349)->primal_0 = _S320;
        (&_S349)->differential_0 = 0.0f;
        DiffPair_float_0 _S350;
        (&_S350)->primal_0 = _S305;
        (&_S350)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S349, &_S350, _S347);
        float _S351 = _S349.differential_0 / _S321;
        float _S352 = - (_S322 * - _S351);
        float _S353 = - (_S323 * _S351);
        y_8 = _S348;
        _S316 = 0.0f;
        _S317 = _S352;
        _S318 = _S350.differential_0;
        _S319 = 0.0f;
        _S320 = _S353;
    }
    float _S354 = (- y_8 + _S316) / _S314;
    float _S355 = _S313 * - _S354;
    float _S356 = _S312 * _S354;
    float _S357 = _S305 * _S356;
    float _S358 = _S310 * _S356;
    DiffPair_float_0 _S359;
    (&_S359)->primal_0 = _S303;
    (&_S359)->differential_0 = 0.0f;
    DiffPair_float_0 _S360;
    (&_S360)->primal_0 = _S305;
    (&_S360)->differential_0 = 0.0f;
    DiffPair_float_0 _S361;
    (&_S361)->primal_0 = _S310;
    (&_S361)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S359, &_S360, &_S361, _S355);
    float _S362 = - ((_S357 + _S361.differential_0 + _S317) / _S311);
    DiffPair_float_0 _S363;
    (&_S363)->primal_0 = _S308;
    (&_S363)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S363, _S362);
    float _S364 = - _S363.differential_0;
    DiffPair_float_0 _S365;
    (&_S365)->primal_0 = _S306;
    (&_S365)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S365, _S335.differential_0);
    DiffPair_float_0 _S366;
    (&_S366)->primal_0 = _S190;
    (&_S366)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S366, _S365.differential_0);
    DiffPair_float_0 _S367 = _S366;
    float _S368 = _S358 + _S360.differential_0 + _S318;
    DiffPair_float_0 _S369;
    (&_S369)->primal_0 = _S304;
    (&_S369)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S369, _S368);
    DiffPair_float_0 _S370;
    (&_S370)->primal_0 = _S189;
    (&_S370)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S370, _S369.differential_0);
    DiffPair_float_0 _S371 = _S370;
    float _S372 = _S359.differential_0 + _S319;
    DiffPair_float_0 _S373;
    (&_S373)->primal_0 = _S302;
    (&_S373)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S373, _S372);
    DiffPair_float_0 _S374;
    (&_S374)->primal_0 = _S188;
    (&_S374)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S374, _S373.differential_0);
    DiffPair_float_0 _S375 = _S374;
    float3  _S376 = make_float3 (0.0f, 0.0f, _S320);
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = _S300;
    (&_S377)->differential_0 = 0.0f;
    DiffPair_float_0 _S378;
    (&_S378)->primal_0 = _S275;
    (&_S378)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S377, &_S378, _s_dOut_0.y);
    DiffPair_float_0 _S379 = _S378;
    DiffPair_float_0 _S380;
    (&_S380)->primal_0 = 0.0f;
    (&_S380)->differential_0 = 0.0f;
    DiffPair_float_0 _S381;
    (&_S381)->primal_0 = y_7;
    (&_S381)->differential_0 = 0.0f;
    _d_max_0(&_S380, &_S381, _S377.differential_0);
    DiffPair_float_0 _S382 = _S381;
    if(_S283)
    {
        float _S383 = a_5 * _S382.differential_0;
        float _S384 = _S284 * _S382.differential_0;
        DiffPair_float_0 _S385;
        (&_S385)->primal_0 = _S285;
        (&_S385)->differential_0 = 0.0f;
        DiffPair_float_0 _S386;
        (&_S386)->primal_0 = _S271;
        (&_S386)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S385, &_S386, _S383);
        float _S387 = _S385.differential_0 / _S286;
        float _S388 = _S269 * - _S387;
        float _S389 = _S278 * _S387;
        y_7 = 0.0f;
        _S284 = _S384;
        _S285 = _S388;
        _S286 = 0.0f;
        _S287 = _S386.differential_0;
        _S288 = _S389;
    }
    else
    {
        float _S390 = - _S382.differential_0;
        float _S391 = b_6 * _S390;
        float _S392 = _S287 * _S390;
        DiffPair_float_0 _S393;
        (&_S393)->primal_0 = _S288;
        (&_S393)->differential_0 = 0.0f;
        DiffPair_float_0 _S394;
        (&_S394)->primal_0 = _S273;
        (&_S394)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S393, &_S394, _S391);
        float _S395 = _S393.differential_0 / _S289;
        float _S396 = - (_S290 * - _S395);
        float _S397 = - (_S291 * _S395);
        y_7 = _S392;
        _S284 = 0.0f;
        _S285 = _S396;
        _S286 = _S394.differential_0;
        _S287 = 0.0f;
        _S288 = _S397;
    }
    float _S398 = (- y_7 + _S284) / _S282;
    float _S399 = _S281 * - _S398;
    float _S400 = _S280 * _S398;
    float _S401 = _S273 * _S400;
    float _S402 = _S278 * _S400;
    DiffPair_float_0 _S403;
    (&_S403)->primal_0 = _S271;
    (&_S403)->differential_0 = 0.0f;
    DiffPair_float_0 _S404;
    (&_S404)->primal_0 = _S273;
    (&_S404)->differential_0 = 0.0f;
    DiffPair_float_0 _S405;
    (&_S405)->primal_0 = _S278;
    (&_S405)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S403, &_S404, &_S405, _S399);
    float _S406 = - ((_S401 + _S405.differential_0 + _S285) / _S279);
    DiffPair_float_0 _S407;
    (&_S407)->primal_0 = _S276;
    (&_S407)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S407, _S406);
    float _S408 = - _S407.differential_0;
    DiffPair_float_0 _S409;
    (&_S409)->primal_0 = _S274;
    (&_S409)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S409, _S379.differential_0);
    DiffPair_float_0 _S410;
    (&_S410)->primal_0 = _S186;
    (&_S410)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S410, _S409.differential_0);
    DiffPair_float_0 _S411 = _S410;
    float _S412 = _S402 + _S404.differential_0 + _S286;
    DiffPair_float_0 _S413;
    (&_S413)->primal_0 = _S272;
    (&_S413)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S413, _S412);
    DiffPair_float_0 _S414;
    (&_S414)->primal_0 = _S185;
    (&_S414)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S414, _S413.differential_0);
    DiffPair_float_0 _S415 = _S414;
    float _S416 = _S403.differential_0 + _S287;
    DiffPair_float_0 _S417;
    (&_S417)->primal_0 = _S270;
    (&_S417)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S417, _S416);
    DiffPair_float_0 _S418;
    (&_S418)->primal_0 = _S184;
    (&_S418)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S418, _S417.differential_0);
    DiffPair_float_0 _S419 = _S418;
    float3  _S420 = _S376 + make_float3 (0.0f, _S288, 0.0f);
    DiffPair_float_0 _S421;
    (&_S421)->primal_0 = _S268;
    (&_S421)->differential_0 = 0.0f;
    DiffPair_float_0 _S422;
    (&_S422)->primal_0 = _S243;
    (&_S422)->differential_0 = 0.0f;
    s_bwd_prop_pow_0(&_S421, &_S422, _s_dOut_0.x);
    DiffPair_float_0 _S423 = _S422;
    DiffPair_float_0 _S424;
    (&_S424)->primal_0 = 0.0f;
    (&_S424)->differential_0 = 0.0f;
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = y_6;
    (&_S425)->differential_0 = 0.0f;
    _d_max_0(&_S424, &_S425, _S421.differential_0);
    DiffPair_float_0 _S426 = _S425;
    if(_S251)
    {
        float _S427 = a_4 * _S426.differential_0;
        float _S428 = _S252 * _S426.differential_0;
        DiffPair_float_0 _S429;
        (&_S429)->primal_0 = _S253;
        (&_S429)->differential_0 = 0.0f;
        DiffPair_float_0 _S430;
        (&_S430)->primal_0 = _S239;
        (&_S430)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S429, &_S430, _S427);
        float _S431 = _S429.differential_0 / _S254;
        float _S432 = _S237 * - _S431;
        float _S433 = _S246 * _S431;
        y_6 = 0.0f;
        _S252 = _S428;
        _S253 = _S432;
        _S254 = 0.0f;
        _S255 = _S430.differential_0;
        _S256 = _S433;
    }
    else
    {
        float _S434 = - _S426.differential_0;
        float _S435 = b_5 * _S434;
        float _S436 = _S255 * _S434;
        DiffPair_float_0 _S437;
        (&_S437)->primal_0 = _S256;
        (&_S437)->differential_0 = 0.0f;
        DiffPair_float_0 _S438;
        (&_S438)->primal_0 = _S241;
        (&_S438)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S437, &_S438, _S435);
        float _S439 = _S437.differential_0 / _S257;
        float _S440 = - (_S258 * - _S439);
        float _S441 = - (_S259 * _S439);
        y_6 = _S436;
        _S252 = 0.0f;
        _S253 = _S440;
        _S254 = _S438.differential_0;
        _S255 = 0.0f;
        _S256 = _S441;
    }
    float _S442 = (- y_6 + _S252) / _S250;
    float _S443 = _S249 * - _S442;
    float _S444 = _S248 * _S442;
    float _S445 = _S241 * _S444;
    float _S446 = _S246 * _S444;
    DiffPair_float_0 _S447;
    (&_S447)->primal_0 = _S239;
    (&_S447)->differential_0 = 0.0f;
    DiffPair_float_0 _S448;
    (&_S448)->primal_0 = _S241;
    (&_S448)->differential_0 = 0.0f;
    DiffPair_float_0 _S449;
    (&_S449)->primal_0 = _S246;
    (&_S449)->differential_0 = 0.0f;
    s_bwd_prop_lerp_0(&_S447, &_S448, &_S449, _S443);
    float _S450 = - ((_S445 + _S449.differential_0 + _S253) / _S247);
    DiffPair_float_0 _S451;
    (&_S451)->primal_0 = _S244;
    (&_S451)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S451, _S450);
    float _S452 = - _S451.differential_0;
    DiffPair_float_0 _S453;
    (&_S453)->primal_0 = _S242;
    (&_S453)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S453, _S423.differential_0);
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S182;
    (&_S454)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S454, _S453.differential_0);
    DiffPair_float_0 _S455 = _S454;
    float _S456 = _S446 + _S448.differential_0 + _S254;
    DiffPair_float_0 _S457;
    (&_S457)->primal_0 = _S240;
    (&_S457)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S457, _S456);
    DiffPair_float_0 _S458;
    (&_S458)->primal_0 = _S181;
    (&_S458)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S458, _S457.differential_0);
    DiffPair_float_0 _S459 = _S458;
    float _S460 = _S447.differential_0 + _S255;
    DiffPair_float_0 _S461;
    (&_S461)->primal_0 = _S238;
    (&_S461)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S461, _S460);
    DiffPair_float_0 _S462;
    (&_S462)->primal_0 = _S180;
    (&_S462)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S462, _S461.differential_0);
    DiffPair_float_0 _S463 = _S462;
    float3  _S464 = _S420 + make_float3 (_S256, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S465;
    (&_S465)->primal_0 = _S233;
    (&_S465)->differential_0 = _S162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S466;
    (&_S466)->primal_0 = _S234;
    (&_S466)->differential_0 = _S162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S467;
    (&_S467)->primal_0 = _S235;
    (&_S467)->differential_0 = _S162;
    s_bwd_prop_clamp_0(&_S465, &_S466, &_S467, _S464);
    float _S468 = - _S465.differential_0.z;
    float3  s_diff_rgi_out_T_0 = make_float3 (_S465.differential_0.x + _S468, _S465.differential_0.y + _S468, _S465.differential_0.z);
    float3  _S469 = _S226 * s_diff_rgi_out_T_0;
    float3  _S470 = _S229 * s_diff_rgi_out_T_0;
    float _S471 = (_S469.x + _S469.y + _S469.z) / _S230;
    float _S472 = intensity_2 * - _S471;
    float _S473 = _S228 * _S471;
    DiffPair_float_0 _S474;
    (&_S474)->primal_0 = _S227;
    (&_S474)->differential_0 = 0.0f;
    DiffPair_float_0 _S475;
    (&_S475)->primal_0 = 0.0f;
    (&_S475)->differential_0 = 0.0f;
    _d_max_0(&_S474, &_S475, _S472);
    float3  _S476 = _S470 + make_float3 (0.0f, 0.0f, _S474.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S477;
    (&_S477)->primal_0 = H_4;
    (&_S477)->differential_0 = _S163;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S478;
    (&_S478)->primal_0 = rgi_in_0;
    (&_S478)->differential_0 = _S162;
    s_bwd_prop_mul_0(&_S477, &_S478, _S476);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S479 = _S477;
    float _S480 = _S473 + _S478.differential_0.z;
    float _S481 = _S478.differential_0.y + _S480;
    float _S482 = _S478.differential_0.x + _S480;
    float3  _S483 = make_float3 (_S482, _S481, _S480);
    if(_S219)
    {
        Matrix<float, 3, 3>  _S484 = _S218 * _S479.differential_0;
        Matrix<float, 3, 3>  _S485 = _S220 * _S479.differential_0;
        _S221 = - ((_S484.rows[int(0)].x + _S484.rows[int(0)].y + _S484.rows[int(0)].z + _S484.rows[int(1)].x + _S484.rows[int(1)].y + _S484.rows[int(1)].z + _S484.rows[int(2)].x + _S484.rows[int(2)].y + _S484.rows[int(2)].z) / _S221);
        H_4 = _S485;
    }
    else
    {
        _S221 = 0.0f;
        H_4 = _S479.differential_0;
    }
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = _S218.rows[int(2)].z;
    (&_S486)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S486, 0.0f);
    float _S487 = _S486.differential_0 + _S221;
    float3  _S488 = _S162;
    *&((&_S488)->z) = _S487;
    Matrix<float, 3, 3>  _S489 = _S163;
    _S489[int(2)] = _S488;
    Matrix<float, 3, 3>  _S490 = H_4 + _S489;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S491;
    (&_S491)->primal_0 = _S217;
    (&_S491)->differential_0 = _S163;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S492;
    (&_S492)->primal_0 = S_inv_0;
    (&_S492)->differential_0 = _S163;
    s_bwd_prop_mul_1(&_S491, &_S492, _S490);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S493;
    (&_S493)->primal_0 = T_2;
    (&_S493)->differential_0 = _S163;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S494;
    (&_S494)->primal_0 = D_0;
    (&_S494)->differential_0 = _S163;
    s_bwd_prop_mul_1(&_S493, &_S494, _S491.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S495 = _S493;
    float3  _S496 = make_float3 (_S494.differential_0.rows[int(0)].x, _S494.differential_0.rows[int(1)].y, _S494.differential_0.rows[int(2)].z);
    float3  _S497;
    if(_S212)
    {
        if(_S214)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S498;
            (&_S498)->primal_0 = r1_2;
            (&_S498)->differential_0 = _S162;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S499;
            (&_S499)->primal_0 = r2_11;
            (&_S499)->differential_0 = _S162;
            s_bwd_prop_cross_0(&_S498, &_S499, _S496);
            _S200 = _S162;
            lambda_v_6 = _S499.differential_0;
            _S497 = _S498.differential_0;
        }
        else
        {
            _S200 = _S496;
            lambda_v_6 = _S162;
            _S497 = _S162;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S500;
        (&_S500)->primal_0 = _S213;
        (&_S500)->differential_0 = _S162;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S501;
        (&_S501)->primal_0 = _S213;
        (&_S501)->differential_0 = _S162;
        s_bwd_prop_dot_0(&_S500, &_S501, 0.0f);
        float3  _S502 = _S501.differential_0 + _S500.differential_0 + _S200;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S503;
        (&_S503)->primal_0 = r0_2;
        (&_S503)->differential_0 = _S162;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S504;
        (&_S504)->primal_0 = r2_11;
        (&_S504)->differential_0 = _S162;
        s_bwd_prop_cross_0(&_S503, &_S504, _S502);
        float3  _S505 = _S504.differential_0 + lambda_v_6;
        _S200 = _S162;
        lambda_v_6 = _S505;
        _S213 = _S497;
        _S497 = _S503.differential_0;
    }
    else
    {
        _S200 = _S496;
        lambda_v_6 = _S162;
        _S213 = _S162;
        _S497 = _S162;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S506;
    (&_S506)->primal_0 = _S211;
    (&_S506)->differential_0 = _S162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S507;
    (&_S507)->primal_0 = _S211;
    (&_S507)->differential_0 = _S162;
    s_bwd_prop_dot_0(&_S506, &_S507, 0.0f);
    float3  _S508 = _S507.differential_0 + _S506.differential_0 + _S200;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S509;
    (&_S509)->primal_0 = r0_2;
    (&_S509)->differential_0 = _S162;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S510;
    (&_S510)->primal_0 = r1_2;
    (&_S510)->differential_0 = _S162;
    s_bwd_prop_cross_0(&_S509, &_S510, _S508);
    float3  _S511 = _S162;
    *&((&_S511)->z) = lambda_v_6.z;
    *&((&_S511)->y) = lambda_v_6.y;
    *&((&_S511)->x) = lambda_v_6.x;
    float3  _S512 = _S510.differential_0 + _S213;
    float3  _S513 = _S162;
    *&((&_S513)->z) = _S512.z;
    *&((&_S513)->y) = _S512.y;
    *&((&_S513)->x) = _S512.x;
    float3  _S514 = _S509.differential_0 + _S497;
    float3  _S515 = _S162;
    *&((&_S515)->z) = _S514.z;
    *&((&_S515)->y) = _S514.y;
    *&((&_S515)->x) = _S514.x;
    Matrix<float, 3, 3>  _S516 = _S163;
    _S516[int(2)] = _S511;
    _S516[int(1)] = _S513;
    _S516[int(0)] = _S515;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S517;
    (&_S517)->primal_0 = skew_0;
    (&_S517)->differential_0 = _S163;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S518;
    (&_S518)->primal_0 = T_2;
    (&_S518)->differential_0 = _S163;
    s_bwd_prop_mul_1(&_S517, &_S518, _S516);
    Matrix<float, 3, 3>  _S519 = _S518.differential_0 + _S495.differential_0;
    float2  _S520 = make_float2 (_S517.differential_0.rows[int(2)].y + - _S517.differential_0.rows[int(1)].z, _S517.differential_0.rows[int(0)].z + - _S517.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S521 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S522;
    (&_S522)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S522)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S523;
    (&_S523)->primal_0 = _S203.color_params_1.n_0;
    (&_S523)->differential_0 = _S166;
    s_bwd_prop_mul_2(&_S522, &_S523, _S520);
    float2  _S524 = make_float2 (_S519.rows[int(0)].z, _S519.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S525;
    (&_S525)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S525)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S526;
    (&_S526)->primal_0 = _S203.color_params_1.g_0;
    (&_S526)->differential_0 = _S166;
    s_bwd_prop_mul_2(&_S525, &_S526, _S524);
    float2  _S527 = make_float2 (_S519.rows[int(0)].y, _S519.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S528;
    (&_S528)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S528)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S529;
    (&_S529)->primal_0 = _S203.color_params_1.r_0;
    (&_S529)->differential_0 = _S166;
    s_bwd_prop_mul_2(&_S528, &_S529, _S527);
    float2  _S530 = make_float2 (_S519.rows[int(0)].x, _S519.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S531;
    (&_S531)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S531)->differential_0 = _S521;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S532;
    (&_S532)->primal_0 = _S203.color_params_1.b_0;
    (&_S532)->differential_0 = _S166;
    s_bwd_prop_mul_2(&_S531, &_S532, _S530);
    ColorPPISPParams_0 _S533 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S533)->n_0 = _S523.differential_0;
    (&_S533)->g_0 = _S526.differential_0;
    (&_S533)->r_0 = _S529.differential_0;
    (&_S533)->b_0 = _S532.differential_0;
    _S200 = _S483;
    *&((&_S200)->z) = 0.0f;
    float _S534 = rgb_out_4.z * _S480;
    float _S535 = _S202 * _S480;
    DiffPair_float_0 _S536;
    (&_S536)->primal_0 = falloff_2;
    (&_S536)->differential_0 = 0.0f;
    DiffPair_float_0 _S537;
    (&_S537)->primal_0 = 0.0f;
    (&_S537)->differential_0 = 0.0f;
    DiffPair_float_0 _S538;
    (&_S538)->primal_0 = 1.0f;
    (&_S538)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S536, &_S537, &_S538, _S534);
    float _S539 = r2_10 * _S536.differential_0;
    float _S540 = r4_8 * _S536.differential_0;
    float s_diff_r6_T_0 = _S179 * _S536.differential_0;
    float _S541 = r6_2 * _S536.differential_0;
    float _S542 = r2_10 * (_S178 * _S536.differential_0 + r2_10 * s_diff_r6_T_0);
    float _S543 = _S177 * _S536.differential_0 + r4_8 * s_diff_r6_T_0 + _S542 + _S542;
    float _S544 = dy_8 * _S543;
    float _S545 = dx_8 * _S543;
    float _S546 = - (_S544 + _S544);
    float _S547 = - (_S545 + _S545);
    *&((&_S200)->y) = 0.0f;
    float _S548 = rgb_out_4.y * _S481;
    float _S549 = _S201 * _S481;
    DiffPair_float_0 _S550;
    (&_S550)->primal_0 = falloff_1;
    (&_S550)->differential_0 = 0.0f;
    DiffPair_float_0 _S551;
    (&_S551)->primal_0 = 0.0f;
    (&_S551)->differential_0 = 0.0f;
    DiffPair_float_0 _S552;
    (&_S552)->primal_0 = 1.0f;
    (&_S552)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S550, &_S551, &_S552, _S548);
    float _S553 = r2_9 * _S550.differential_0;
    float _S554 = r4_7 * _S550.differential_0;
    float s_diff_r6_T_1 = _S176 * _S550.differential_0;
    float _S555 = r6_1 * _S550.differential_0;
    float _S556 = r2_9 * (_S175 * _S550.differential_0 + r2_9 * s_diff_r6_T_1);
    float _S557 = _S174 * _S550.differential_0 + r4_7 * s_diff_r6_T_1 + _S556 + _S556;
    float _S558 = dy_7 * _S557;
    float _S559 = dx_7 * _S557;
    float _S560 = - (_S558 + _S558);
    float _S561 = - (_S559 + _S559);
    *&((&_S200)->x) = 0.0f;
    float _S562 = rgb_out_4.x * _S482;
    float _S563 = _S198 * _S482;
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = falloff_0;
    (&_S564)->differential_0 = 0.0f;
    DiffPair_float_0 _S565;
    (&_S565)->primal_0 = 0.0f;
    (&_S565)->differential_0 = 0.0f;
    DiffPair_float_0 _S566;
    (&_S566)->primal_0 = 1.0f;
    (&_S566)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S564, &_S565, &_S566, _S562);
    float _S567 = r2_8 * _S564.differential_0;
    float _S568 = r4_6 * _S564.differential_0;
    float s_diff_r6_T_2 = _S173 * _S564.differential_0;
    float _S569 = r6_0 * _S564.differential_0;
    float _S570 = r2_8 * (_S172 * _S564.differential_0 + r2_8 * s_diff_r6_T_2);
    float _S571 = _S171 * _S564.differential_0 + r4_6 * s_diff_r6_T_2 + _S570 + _S570;
    float _S572 = dy_6 * _S571;
    float _S573 = dx_6 * _S571;
    float _S574 = - (_S572 + _S572);
    float _S575 = - (_S573 + _S573);
    float3  _S576 = _S162;
    *&((&_S576)->z) = _S535;
    *&((&_S576)->y) = _S549;
    *&((&_S576)->x) = _S563;
    float3  _S577 = _S200 + _S576;
    float3  _S578 = _S161.primal_0 * _S577;
    float3  _S579 = _S194 * _S577;
    float _S580 = _S578.x + _S578.y + _S578.z;
    DiffPair_float_0 _S581;
    (&_S581)->primal_0 = _S192.exposure_1;
    (&_S581)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S581, _S580);
    PPISPParams_0 _S582 = PPISPParams_x24_syn_dzero_0();
    (&_S582)->color_params_1 = _S533;
    (&_S582)->exposure_1 = _S581.differential_0;
    _S170 = _S582;
    (&(&_S170)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S583 = _S582.crf_params_1[int(2)].center_0 + _S364;
    (&(&_S170)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S584 = _S582.crf_params_1[int(2)].gamma_0 + _S367.differential_0;
    (&(&_S170)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S585 = _S582.crf_params_1[int(2)].shoulder_0 + _S371.differential_0;
    (&(&_S170)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S586 = _S582.crf_params_1[int(2)].toe_0 + _S375.differential_0;
    (&(&_S170)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S587 = _S582.crf_params_1[int(1)].center_0 + _S408;
    (&(&_S170)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S588 = _S582.crf_params_1[int(1)].gamma_0 + _S411.differential_0;
    (&(&_S170)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S589 = _S582.crf_params_1[int(1)].shoulder_0 + _S415.differential_0;
    (&(&_S170)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S590 = _S582.crf_params_1[int(1)].toe_0 + _S419.differential_0;
    (&(&_S170)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S591 = _S582.crf_params_1[int(0)].center_0 + _S452;
    (&(&_S170)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S592 = _S582.crf_params_1[int(0)].gamma_0 + _S455.differential_0;
    (&(&_S170)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S593 = _S582.crf_params_1[int(0)].shoulder_0 + _S459.differential_0;
    (&(&_S170)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S594 = _S582.crf_params_1[int(0)].toe_0 + _S463.differential_0;
    *&((&(&(&_S170)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S170)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S170)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S595 = _S541 + _S582.vignette_params_1[int(2)].alpha2_0;
    (&(&_S170)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S596 = _S540 + _S582.vignette_params_1[int(2)].alpha1_0;
    (&(&_S170)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S597 = _S539 + _S582.vignette_params_1[int(2)].alpha0_0;
    (&(&_S170)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S598 = _S546 + _S582.vignette_params_1[int(2)].cy_0;
    (&(&_S170)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S599 = _S547 + _S582.vignette_params_1[int(2)].cx_0;
    (&(&_S170)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S600 = _S555 + _S582.vignette_params_1[int(1)].alpha2_0;
    (&(&_S170)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S601 = _S554 + _S582.vignette_params_1[int(1)].alpha1_0;
    (&(&_S170)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S602 = _S553 + _S582.vignette_params_1[int(1)].alpha0_0;
    (&(&_S170)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S603 = _S560 + _S582.vignette_params_1[int(1)].cy_0;
    (&(&_S170)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S604 = _S561 + _S582.vignette_params_1[int(1)].cx_0;
    (&(&_S170)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S605 = _S569 + _S582.vignette_params_1[int(0)].alpha2_0;
    (&(&_S170)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S606 = _S568 + _S582.vignette_params_1[int(0)].alpha1_0;
    (&(&_S170)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S607 = _S567 + _S582.vignette_params_1[int(0)].alpha0_0;
    (&(&_S170)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S608 = _S574 + _S582.vignette_params_1[int(0)].cy_0;
    (&(&_S170)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S609 = _S575 + _S582.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S610;
    _S610[int(0)] = 0.0f;
    _S610[int(1)] = 0.0f;
    _S610[int(2)] = 0.0f;
    _S610[int(3)] = 0.0f;
    _S610[int(4)] = 0.0f;
    _S610[int(5)] = 0.0f;
    _S610[int(6)] = 0.0f;
    _S610[int(7)] = 0.0f;
    _S610[int(8)] = 0.0f;
    _S610[int(9)] = 0.0f;
    _S610[int(10)] = 0.0f;
    _S610[int(11)] = 0.0f;
    _S610[int(12)] = 0.0f;
    _S610[int(13)] = 0.0f;
    _S610[int(14)] = 0.0f;
    _S610[int(15)] = 0.0f;
    _S610[int(16)] = 0.0f;
    _S610[int(17)] = 0.0f;
    _S610[int(18)] = 0.0f;
    _S610[int(19)] = 0.0f;
    _S610[int(20)] = 0.0f;
    _S610[int(21)] = 0.0f;
    _S610[int(22)] = 0.0f;
    _S610[int(23)] = 0.0f;
    _S610[int(24)] = 0.0f;
    _S610[int(25)] = 0.0f;
    _S610[int(26)] = 0.0f;
    _S610[int(27)] = 0.0f;
    _S610[int(28)] = 0.0f;
    _S610[int(29)] = 0.0f;
    _S610[int(30)] = 0.0f;
    _S610[int(31)] = 0.0f;
    _S610[int(32)] = 0.0f;
    _S610[int(33)] = 0.0f;
    _S610[int(34)] = 0.0f;
    _S610[int(35)] = 0.0f;
    _S610[int(8)] = _S602;
    _S610[int(16)] = _S582.color_params_1.b_0.x;
    _S610[int(15)] = _S595;
    _S610[int(14)] = _S596;
    _S610[int(13)] = _S597;
    _S610[int(12)] = _S598;
    _S610[int(11)] = _S599;
    _S610[int(10)] = _S600;
    _S610[int(9)] = _S601;
    _S610[int(17)] = _S582.color_params_1.b_0.y;
    _S610[int(7)] = _S603;
    _S610[int(6)] = _S604;
    _S610[int(5)] = _S605;
    _S610[int(4)] = _S606;
    _S610[int(3)] = _S607;
    _S610[int(2)] = _S608;
    _S610[int(1)] = _S609;
    _S610[int(0)] = _S170.exposure_1;
    _S610[int(26)] = _S592;
    _S610[int(34)] = _S584;
    _S610[int(33)] = _S585;
    _S610[int(32)] = _S586;
    _S610[int(31)] = _S587;
    _S610[int(30)] = _S588;
    _S610[int(29)] = _S589;
    _S610[int(28)] = _S590;
    _S610[int(27)] = _S591;
    _S610[int(35)] = _S583;
    _S610[int(25)] = _S593;
    _S610[int(24)] = _S594;
    _S610[int(23)] = _S582.color_params_1.n_0.y;
    _S610[int(22)] = _S582.color_params_1.n_0.x;
    _S610[int(21)] = _S582.color_params_1.g_0.y;
    _S610[int(20)] = _S582.color_params_1.g_0.x;
    _S610[int(19)] = _S582.color_params_1.r_0.y;
    _S610[int(18)] = _S582.color_params_1.r_0.x;
    dpparams_0->primal_0 = dpparams_0->primal_0;
    dpparams_0->differential_0 = _S610;
    dprgb_in_0->primal_0 = (*dprgb_in_0).primal_0;
    dprgb_in_0->differential_0 = _S579;
    return;
}

inline __device__ void s_bwd_apply_ppisp_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S611, float2  _S612, float2  _S613, float2  _S614, DiffPair_arrayx3Cfloatx2C36x3E_0 * _S615, float3  _S616)
{
    s_bwd_prop_apply_ppisp_0(_S611, _S612, _S613, _S614, _S615, _S616);
    return;
}

inline __device__ void apply_ppisp_vjp(float3  rgb_in_2, float2  pix_coord_3, float2  image_center_3, float2  img_size_3, FixedArray<float, 36>  params_2, float3  grad_out_0, float3  * grad_rgb_in_0, FixedArray<float, 36>  * grad_params_0)
{
    float3  _S617 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_0;
    (&dp_rgb_in_0)->primal_0 = rgb_in_2;
    (&dp_rgb_in_0)->differential_0 = _S617;
    FixedArray<float, 36>  _S618 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_0;
    (&dp_params_0)->primal_0 = params_2;
    (&dp_params_0)->differential_0 = _S618;
    s_bwd_apply_ppisp_0(&dp_rgb_in_0, pix_coord_3, image_center_3, img_size_3, &dp_params_0, grad_out_0);
    *grad_rgb_in_0 = dp_rgb_in_0.differential_0;
    *grad_params_0 = (&dp_params_0)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C39x3E_0
{
    FixedArray<float, 39>  primal_0;
    FixedArray<float, 39>  differential_0;
};

inline __device__ void s_bwd_prop_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_in_1, float2  pix_coord_4, float2  image_center_4, float2  img_size_4, DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_1, float3  _s_dOut_1)
{
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S619 = *dprgb_in_1;
    float3  _S620 = make_float3 (0.0f);
    Matrix<float, 3, 3>  _S621 = makeMatrix<float, 3, 3> (0.0f);
    VignettingChannelParams_0 _S622 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S623 = {
        _S622, _S622, _S622
    };
    float2  _S624 = make_float2 (0.0f);
    ColorPPISPParams_0 _S625 = { _S624, _S624, _S624, _S624 };
    RQSCRFPPISPChannelParams_0 _S626 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S627 = {
        _S626, _S626, _S626
    };
    PPISPParamsRQS_0 _S628;
    (&_S628)->exposure_0 = dpparams_1->primal_0[int(0)];
    (&_S628)->vignette_params_0 = _S623;
    (&_S628)->color_params_0 = _S625;
    (&_S628)->crf_params_0 = _S627;
    (&(&_S628)->vignette_params_0[int(0)])->cx_0 = dpparams_1->primal_0[int(1)];
    (&(&_S628)->vignette_params_0[int(0)])->cy_0 = dpparams_1->primal_0[int(2)];
    float _S629 = dpparams_1->primal_0[int(3)];
    (&(&_S628)->vignette_params_0[int(0)])->alpha0_0 = dpparams_1->primal_0[int(3)];
    float _S630 = dpparams_1->primal_0[int(4)];
    (&(&_S628)->vignette_params_0[int(0)])->alpha1_0 = dpparams_1->primal_0[int(4)];
    float _S631 = dpparams_1->primal_0[int(5)];
    (&(&_S628)->vignette_params_0[int(0)])->alpha2_0 = dpparams_1->primal_0[int(5)];
    (&(&_S628)->vignette_params_0[int(1)])->cx_0 = dpparams_1->primal_0[int(6)];
    (&(&_S628)->vignette_params_0[int(1)])->cy_0 = dpparams_1->primal_0[int(7)];
    float _S632 = dpparams_1->primal_0[int(8)];
    (&(&_S628)->vignette_params_0[int(1)])->alpha0_0 = dpparams_1->primal_0[int(8)];
    float _S633 = dpparams_1->primal_0[int(9)];
    (&(&_S628)->vignette_params_0[int(1)])->alpha1_0 = dpparams_1->primal_0[int(9)];
    float _S634 = dpparams_1->primal_0[int(10)];
    (&(&_S628)->vignette_params_0[int(1)])->alpha2_0 = dpparams_1->primal_0[int(10)];
    (&(&_S628)->vignette_params_0[int(2)])->cx_0 = dpparams_1->primal_0[int(11)];
    (&(&_S628)->vignette_params_0[int(2)])->cy_0 = dpparams_1->primal_0[int(12)];
    float _S635 = dpparams_1->primal_0[int(13)];
    (&(&_S628)->vignette_params_0[int(2)])->alpha0_0 = dpparams_1->primal_0[int(13)];
    float _S636 = dpparams_1->primal_0[int(14)];
    (&(&_S628)->vignette_params_0[int(2)])->alpha1_0 = dpparams_1->primal_0[int(14)];
    float _S637 = dpparams_1->primal_0[int(15)];
    (&(&_S628)->vignette_params_0[int(2)])->alpha2_0 = dpparams_1->primal_0[int(15)];
    *&((&(&(&_S628)->color_params_0)->b_0)->x) = dpparams_1->primal_0[int(16)];
    *&((&(&(&_S628)->color_params_0)->b_0)->y) = dpparams_1->primal_0[int(17)];
    *&((&(&(&_S628)->color_params_0)->r_0)->x) = dpparams_1->primal_0[int(18)];
    *&((&(&(&_S628)->color_params_0)->r_0)->y) = dpparams_1->primal_0[int(19)];
    *&((&(&(&_S628)->color_params_0)->g_0)->x) = dpparams_1->primal_0[int(20)];
    *&((&(&(&_S628)->color_params_0)->g_0)->y) = dpparams_1->primal_0[int(21)];
    *&((&(&(&_S628)->color_params_0)->n_0)->x) = dpparams_1->primal_0[int(22)];
    *&((&(&(&_S628)->color_params_0)->n_0)->y) = dpparams_1->primal_0[int(23)];
    float _S638 = dpparams_1->primal_0[int(24)];
    (&(&_S628)->crf_params_0[int(0)])->g0_0 = dpparams_1->primal_0[int(24)];
    float _S639 = dpparams_1->primal_0[int(25)];
    (&(&_S628)->crf_params_0[int(0)])->g1_0 = dpparams_1->primal_0[int(25)];
    float _S640 = dpparams_1->primal_0[int(26)];
    (&(&_S628)->crf_params_0[int(0)])->x0_0 = dpparams_1->primal_0[int(26)];
    float _S641 = dpparams_1->primal_0[int(27)];
    (&(&_S628)->crf_params_0[int(0)])->y0_0 = dpparams_1->primal_0[int(27)];
    float _S642 = dpparams_1->primal_0[int(28)];
    (&(&_S628)->crf_params_0[int(0)])->gc_0 = dpparams_1->primal_0[int(28)];
    float _S643 = dpparams_1->primal_0[int(29)];
    (&(&_S628)->crf_params_0[int(1)])->g0_0 = dpparams_1->primal_0[int(29)];
    float _S644 = dpparams_1->primal_0[int(30)];
    (&(&_S628)->crf_params_0[int(1)])->g1_0 = dpparams_1->primal_0[int(30)];
    float _S645 = dpparams_1->primal_0[int(31)];
    (&(&_S628)->crf_params_0[int(1)])->x0_0 = dpparams_1->primal_0[int(31)];
    float _S646 = dpparams_1->primal_0[int(32)];
    (&(&_S628)->crf_params_0[int(1)])->y0_0 = dpparams_1->primal_0[int(32)];
    float _S647 = dpparams_1->primal_0[int(33)];
    (&(&_S628)->crf_params_0[int(1)])->gc_0 = dpparams_1->primal_0[int(33)];
    float _S648 = dpparams_1->primal_0[int(34)];
    (&(&_S628)->crf_params_0[int(2)])->g0_0 = dpparams_1->primal_0[int(34)];
    float _S649 = dpparams_1->primal_0[int(35)];
    (&(&_S628)->crf_params_0[int(2)])->g1_0 = dpparams_1->primal_0[int(35)];
    float _S650 = dpparams_1->primal_0[int(36)];
    (&(&_S628)->crf_params_0[int(2)])->x0_0 = dpparams_1->primal_0[int(36)];
    float _S651 = dpparams_1->primal_0[int(37)];
    (&(&_S628)->crf_params_0[int(2)])->y0_0 = dpparams_1->primal_0[int(37)];
    float _S652 = dpparams_1->primal_0[int(38)];
    (&(&_S628)->crf_params_0[int(2)])->gc_0 = dpparams_1->primal_0[int(38)];
    PPISPParamsRQS_0 _S653 = _S628;
    float _S654 = s_primal_ctx_exp2_0(_S628.exposure_0);
    float3  _S655 = make_float3 (_S654);
    float3  rgb_out_5 = (*dprgb_in_1).primal_0 * make_float3 (_S654);
    float _S656 = (F32_max((img_size_4.x), (img_size_4.y)));
    float _S657 = (pix_coord_4.x - image_center_4.x) / _S656;
    float _S658 = (pix_coord_4.y - image_center_4.y) / _S656;
    float dx_9 = _S657 - dpparams_1->primal_0[int(1)];
    float dy_9 = _S658 - dpparams_1->primal_0[int(2)];
    float r2_12 = dx_9 * dx_9 + dy_9 * dy_9;
    float r4_9 = r2_12 * r2_12;
    float r6_3 = r4_9 * r2_12;
    float falloff_3 = dpparams_1->primal_0[int(5)] * r6_3 + dpparams_1->primal_0[int(4)] * r4_9 + dpparams_1->primal_0[int(3)] * r2_12 + 1.0f;
    float _S659 = s_primal_ctx_clamp_0(falloff_3, 0.0f, 1.0f);
    float _S660 = rgb_out_5.x * _S659;
    float3  _S661 = rgb_out_5;
    *&((&_S661)->x) = _S660;
    float dx_10 = _S657 - dpparams_1->primal_0[int(6)];
    float dy_10 = _S658 - dpparams_1->primal_0[int(7)];
    float r2_13 = dx_10 * dx_10 + dy_10 * dy_10;
    float r4_10 = r2_13 * r2_13;
    float r6_4 = r4_10 * r2_13;
    float falloff_4 = dpparams_1->primal_0[int(10)] * r6_4 + dpparams_1->primal_0[int(9)] * r4_10 + dpparams_1->primal_0[int(8)] * r2_13 + 1.0f;
    float _S662 = s_primal_ctx_clamp_0(falloff_4, 0.0f, 1.0f);
    *&((&_S661)->y) = rgb_out_5.y * _S662;
    float dx_11 = _S657 - dpparams_1->primal_0[int(11)];
    float dy_11 = _S658 - dpparams_1->primal_0[int(12)];
    float r2_14 = dx_11 * dx_11 + dy_11 * dy_11;
    float r4_11 = r2_14 * r2_14;
    float r6_5 = r4_11 * r2_14;
    float falloff_5 = dpparams_1->primal_0[int(15)] * r6_5 + dpparams_1->primal_0[int(14)] * r4_11 + dpparams_1->primal_0[int(13)] * r2_14 + 1.0f;
    float _S663 = s_primal_ctx_clamp_0(falloff_5, 0.0f, 1.0f);
    *&((&_S661)->z) = rgb_out_5.z * _S663;
    PPISPParamsRQS_0 _S664 = _S628;
    float2  _S665 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), _S628.color_params_0.b_0);
    float2  _S666 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), _S628.color_params_0.r_0);
    float2  _S667 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), _S628.color_params_0.g_0);
    float2  _S668 = s_primal_ctx_mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), _S628.color_params_0.n_0);
    float _S669 = 0.3333333432674408f + _S668.x;
    float _S670 = 0.3333333432674408f + _S668.y;
    Matrix<float, 3, 3>  T_3 = makeMatrix<float, 3, 3> (_S665.x, 1.0f + _S666.x, _S667.x, _S665.y, _S666.y, 1.0f + _S667.y, 1.0f, 1.0f, 1.0f);
    Matrix<float, 3, 3>  skew_1 = makeMatrix<float, 3, 3> (0.0f, -1.0f, _S670, 1.0f, 0.0f, - _S669, - _S670, _S669, 0.0f);
    Matrix<float, 3, 3>  _S671 = s_primal_ctx_mul_1(skew_1, T_3);
    float3  r0_3 = make_float3 (_S671.rows[int(0)].x, _S671.rows[int(0)].y, _S671.rows[int(0)].z);
    float3  r1_3 = make_float3 (_S671.rows[int(1)].x, _S671.rows[int(1)].y, _S671.rows[int(1)].z);
    float3  r2_15 = make_float3 (_S671.rows[int(2)].x, _S671.rows[int(2)].y, _S671.rows[int(2)].z);
    float3  _S672 = s_primal_ctx_cross_0(r0_3, r1_3);
    bool _S673 = (s_primal_ctx_dot_0(_S672, _S672)) < 9.99999968265522539e-21f;
    float3  lambda_v_7;
    float3  _S674;
    bool _S675;
    if(_S673)
    {
        float3  _S676 = s_primal_ctx_cross_0(r0_3, r2_15);
        bool _S677 = (s_primal_ctx_dot_0(_S676, _S676)) < 9.99999968265522539e-21f;
        if(_S677)
        {
            lambda_v_7 = s_primal_ctx_cross_0(r1_3, r2_15);
        }
        else
        {
            lambda_v_7 = _S676;
        }
        _S675 = _S677;
        _S674 = _S676;
    }
    else
    {
        lambda_v_7 = _S672;
        _S675 = false;
        _S674 = _S620;
    }
    Matrix<float, 3, 3>  S_inv_1 = makeMatrix<float, 3, 3> (-1.0f, -1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    Matrix<float, 3, 3>  D_1 = makeMatrix<float, 3, 3> (lambda_v_7.x, 0.0f, 0.0f, 0.0f, lambda_v_7.y, 0.0f, 0.0f, 0.0f, lambda_v_7.z);
    Matrix<float, 3, 3>  _S678 = s_primal_ctx_mul_1(T_3, D_1);
    Matrix<float, 3, 3>  _S679 = s_primal_ctx_mul_1(_S678, S_inv_1);
    bool _S680 = (s_primal_ctx_abs_0(_S679.rows[int(2)].z)) > 9.99999968265522539e-21f;
    Matrix<float, 3, 3>  H_5;
    Matrix<float, 3, 3>  _S681;
    float _S682;
    if(_S680)
    {
        float inv_s_1 = 1.0f / _S679.rows[int(2)].z;
        Matrix<float, 3, 3>  _S683 = makeMatrix<float, 3, 3> (inv_s_1);
        float _S684 = _S679.rows[int(2)].z * _S679.rows[int(2)].z;
        H_5 = _S679 * makeMatrix<float, 3, 3> (inv_s_1);
        _S681 = _S683;
        _S682 = _S684;
    }
    else
    {
        H_5 = _S679;
        _S681 = _S621;
        _S682 = 0.0f;
    }
    float _S685 = _S661.x;
    float _S686 = _S661.y;
    float intensity_3 = _S685 + _S686 + _S661.z;
    float3  rgi_in_1 = make_float3 (_S685, _S686, intensity_3);
    float3  _S687 = s_primal_ctx_mul_2(H_5, rgi_in_1);
    float _S688 = _S687.z;
    float _S689 = (F32_max((_S688), (0.0f))) + 0.00000999999974738f;
    float norm_factor_1 = intensity_3 / _S689;
    float3  _S690 = make_float3 (norm_factor_1);
    float _S691 = _S689 * _S689;
    float3  rgi_out_5 = _S687 * make_float3 (norm_factor_1);
    float _S692 = rgi_out_5.x;
    float _S693 = rgi_out_5.y;
    float3  _S694 = make_float3 (_S692, _S693, rgi_out_5.z - _S692 - _S693);
    float3  _S695 = make_float3 (0.0f);
    float3  _S696 = make_float3 (1.0f);
    float3  _S697 = s_primal_ctx_clamp_1(_S694, _S695, _S696);
    float _S698 = _S697.x;
    float _S699 = s_primal_ctx_exp_0(_S638);
    float _S700 = s_primal_ctx_exp_0(_S639);
    float _S701 = - _S640;
    float _S702 = 1.0f + s_primal_ctx_exp_0(_S701);
    float x0_4 = 1.0f / _S702;
    float _S703 = _S702 * _S702;
    float _S704 = - _S641;
    float _S705 = 1.0f + s_primal_ctx_exp_0(_S704);
    float y0_4 = 1.0f / _S705;
    float _S706 = _S705 * _S705;
    float _S707 = s_primal_ctx_exp_0(_S642);
    bool _S708 = _S698 < x0_4;
    float _S709;
    float _S710;
    float _S711;
    float _S712;
    float _S713;
    float _S714;
    float _S715;
    float _S716;
    float _S717;
    float _S718;
    float _S719;
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
    if(_S708)
    {
        float s0_3 = y0_4 / x0_4;
        float _S736 = x0_4 * x0_4;
        float t0_3 = _S698 / x0_4;
        float _S737 = s0_3 * t0_3;
        float _S738 = _S699 * t0_3;
        float _S739 = 1.0f - t0_3;
        float _S740 = _S737 * t0_3 + _S738 * _S739;
        float _S741 = y0_4 * _S740;
        float _S742 = _S699 + _S707 - 2.0f * s0_3;
        float _S743 = _S742 * t0_3;
        float _S744 = s0_3 + _S743 * _S739;
        _S709 = _S744 * _S744;
        _S710 = _S741;
        _S711 = _S744;
        _S712 = _S743;
        _S713 = _S739;
        _S714 = _S742;
        _S715 = t0_3;
        _S716 = _S740;
        _S717 = _S738;
        _S718 = _S737;
        _S719 = s0_3;
        _S720 = _S736;
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
        _S732 = 0.0f;
        _S733 = 0.0f;
        _S734 = 0.0f;
        _S735 = 0.0f;
    }
    else
    {
        float _S745 = 1.0f - y0_4;
        float _S746 = 1.0f - x0_4;
        float s1_3 = _S745 / _S746;
        float _S747 = _S746 * _S746;
        float _S748 = _S698 - x0_4;
        float t1_3 = _S748 / _S746;
        float _S749 = s1_3 * t1_3;
        float _S750 = _S707 * t1_3;
        float _S751 = 1.0f - t1_3;
        float _S752 = _S749 * t1_3 + _S750 * _S751;
        float _S753 = _S745 * _S752;
        float _S754 = _S707 + _S700 - 2.0f * s1_3;
        float _S755 = _S754 * t1_3;
        float _S756 = s1_3 + _S755 * _S751;
        float _S757 = _S756 * _S756;
        _S709 = 0.0f;
        _S710 = 0.0f;
        _S711 = 0.0f;
        _S712 = 0.0f;
        _S713 = 0.0f;
        _S714 = 0.0f;
        _S715 = 0.0f;
        _S716 = 0.0f;
        _S717 = 0.0f;
        _S718 = 0.0f;
        _S719 = 0.0f;
        _S720 = 0.0f;
        _S721 = _S757;
        _S722 = _S753;
        _S723 = _S756;
        _S724 = _S755;
        _S725 = _S751;
        _S726 = _S754;
        _S727 = t1_3;
        _S728 = _S745;
        _S729 = _S752;
        _S730 = _S750;
        _S731 = _S749;
        _S732 = s1_3;
        _S733 = _S747;
        _S734 = _S748;
        _S735 = _S746;
    }
    float _S758 = _S697.y;
    float _S759 = s_primal_ctx_exp_0(_S643);
    float _S760 = s_primal_ctx_exp_0(_S644);
    float _S761 = - _S645;
    float _S762 = 1.0f + s_primal_ctx_exp_0(_S761);
    float x0_5 = 1.0f / _S762;
    float _S763 = _S762 * _S762;
    float _S764 = - _S646;
    float _S765 = 1.0f + s_primal_ctx_exp_0(_S764);
    float y0_5 = 1.0f / _S765;
    float _S766 = _S765 * _S765;
    float _S767 = s_primal_ctx_exp_0(_S647);
    bool _S768 = _S758 < x0_5;
    float _S769;
    float _S770;
    float _S771;
    float _S772;
    float _S773;
    float _S774;
    float _S775;
    float _S776;
    float _S777;
    float _S778;
    float _S779;
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
    if(_S768)
    {
        float s0_4 = y0_5 / x0_5;
        float _S796 = x0_5 * x0_5;
        float t0_4 = _S758 / x0_5;
        float _S797 = s0_4 * t0_4;
        float _S798 = _S759 * t0_4;
        float _S799 = 1.0f - t0_4;
        float _S800 = _S797 * t0_4 + _S798 * _S799;
        float _S801 = y0_5 * _S800;
        float _S802 = _S759 + _S767 - 2.0f * s0_4;
        float _S803 = _S802 * t0_4;
        float _S804 = s0_4 + _S803 * _S799;
        _S769 = _S804 * _S804;
        _S770 = _S801;
        _S771 = _S804;
        _S772 = _S803;
        _S773 = _S799;
        _S774 = _S802;
        _S775 = t0_4;
        _S776 = _S800;
        _S777 = _S798;
        _S778 = _S797;
        _S779 = s0_4;
        _S780 = _S796;
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
        _S792 = 0.0f;
        _S793 = 0.0f;
        _S794 = 0.0f;
        _S795 = 0.0f;
    }
    else
    {
        float _S805 = 1.0f - y0_5;
        float _S806 = 1.0f - x0_5;
        float s1_4 = _S805 / _S806;
        float _S807 = _S806 * _S806;
        float _S808 = _S758 - x0_5;
        float t1_4 = _S808 / _S806;
        float _S809 = s1_4 * t1_4;
        float _S810 = _S767 * t1_4;
        float _S811 = 1.0f - t1_4;
        float _S812 = _S809 * t1_4 + _S810 * _S811;
        float _S813 = _S805 * _S812;
        float _S814 = _S767 + _S760 - 2.0f * s1_4;
        float _S815 = _S814 * t1_4;
        float _S816 = s1_4 + _S815 * _S811;
        float _S817 = _S816 * _S816;
        _S769 = 0.0f;
        _S770 = 0.0f;
        _S771 = 0.0f;
        _S772 = 0.0f;
        _S773 = 0.0f;
        _S774 = 0.0f;
        _S775 = 0.0f;
        _S776 = 0.0f;
        _S777 = 0.0f;
        _S778 = 0.0f;
        _S779 = 0.0f;
        _S780 = 0.0f;
        _S781 = _S817;
        _S782 = _S813;
        _S783 = _S816;
        _S784 = _S815;
        _S785 = _S811;
        _S786 = _S814;
        _S787 = t1_4;
        _S788 = _S805;
        _S789 = _S812;
        _S790 = _S810;
        _S791 = _S809;
        _S792 = s1_4;
        _S793 = _S807;
        _S794 = _S808;
        _S795 = _S806;
    }
    float _S818 = _S697.z;
    float _S819 = s_primal_ctx_exp_0(_S648);
    float _S820 = s_primal_ctx_exp_0(_S649);
    float _S821 = - _S650;
    float _S822 = 1.0f + s_primal_ctx_exp_0(_S821);
    float x0_6 = 1.0f / _S822;
    float _S823 = _S822 * _S822;
    float _S824 = - _S651;
    float _S825 = 1.0f + s_primal_ctx_exp_0(_S824);
    float y0_6 = 1.0f / _S825;
    float _S826 = _S825 * _S825;
    float _S827 = s_primal_ctx_exp_0(_S652);
    bool _S828 = _S818 < x0_6;
    float _S829;
    float _S830;
    float _S831;
    float _S832;
    float _S833;
    float _S834;
    float _S835;
    float _S836;
    float _S837;
    float _S838;
    float _S839;
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
    if(_S828)
    {
        float s0_5 = y0_6 / x0_6;
        float _S856 = x0_6 * x0_6;
        float t0_5 = _S818 / x0_6;
        float _S857 = s0_5 * t0_5;
        float _S858 = _S819 * t0_5;
        float _S859 = 1.0f - t0_5;
        float _S860 = _S857 * t0_5 + _S858 * _S859;
        float _S861 = y0_6 * _S860;
        float _S862 = _S819 + _S827 - 2.0f * s0_5;
        float _S863 = _S862 * t0_5;
        float _S864 = s0_5 + _S863 * _S859;
        _S829 = _S864 * _S864;
        _S830 = _S861;
        _S831 = _S864;
        _S832 = _S863;
        _S833 = _S859;
        _S834 = _S862;
        _S835 = t0_5;
        _S836 = _S860;
        _S837 = _S858;
        _S838 = _S857;
        _S839 = s0_5;
        _S840 = _S856;
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
        _S852 = 0.0f;
        _S853 = 0.0f;
        _S854 = 0.0f;
        _S855 = 0.0f;
    }
    else
    {
        float _S865 = 1.0f - y0_6;
        float _S866 = 1.0f - x0_6;
        float s1_5 = _S865 / _S866;
        float _S867 = _S866 * _S866;
        float _S868 = _S818 - x0_6;
        float t1_5 = _S868 / _S866;
        float _S869 = s1_5 * t1_5;
        float _S870 = _S827 * t1_5;
        float _S871 = 1.0f - t1_5;
        float _S872 = _S869 * t1_5 + _S870 * _S871;
        float _S873 = _S865 * _S872;
        float _S874 = _S827 + _S820 - 2.0f * s1_5;
        float _S875 = _S874 * t1_5;
        float _S876 = s1_5 + _S875 * _S871;
        float _S877 = _S876 * _S876;
        _S829 = 0.0f;
        _S830 = 0.0f;
        _S831 = 0.0f;
        _S832 = 0.0f;
        _S833 = 0.0f;
        _S834 = 0.0f;
        _S835 = 0.0f;
        _S836 = 0.0f;
        _S837 = 0.0f;
        _S838 = 0.0f;
        _S839 = 0.0f;
        _S840 = 0.0f;
        _S841 = _S877;
        _S842 = _S873;
        _S843 = _S876;
        _S844 = _S875;
        _S845 = _S871;
        _S846 = _S874;
        _S847 = t1_5;
        _S848 = _S865;
        _S849 = _S872;
        _S850 = _S870;
        _S851 = _S869;
        _S852 = s1_5;
        _S853 = _S867;
        _S854 = _S868;
        _S855 = _S866;
    }
    if(_S828)
    {
        float _S878 = _s_dOut_1.z / _S829;
        float _S879 = _S830 * - _S878;
        float _S880 = _S831 * _S878;
        float _S881 = _S833 * _S879;
        float _S882 = _S835 * _S881;
        float _S883 = y0_6 * _S880;
        float _S884 = _S833 * _S883;
        float _S885 = _S835 * _S883;
        float _S886 = (_S834 * _S881 + - (_S832 * _S879 + _S837 * _S883) + _S819 * _S884 + _S838 * _S883 + _S839 * _S885) / _S840;
        float _S887 = x0_6 * _S886;
        float _S888 = (_S879 + 2.0f * - _S882 + _S835 * _S885) / _S840;
        float _S889 = _S836 * _S880 + x0_6 * _S888;
        float _S890 = _S882 + _S835 * _S884;
        float _S891 = _S818 * - _S886 + y0_6 * - _S888;
        _S829 = _S882;
        _S830 = _S889;
        _S831 = _S891;
        _S832 = 0.0f;
        _S833 = _S890;
        _S834 = _S887;
    }
    else
    {
        float _S892 = _s_dOut_1.z / _S841;
        float _S893 = _S842 * - _S892;
        float _S894 = _S843 * _S892;
        float _S895 = _S845 * _S893;
        float _S896 = _S847 * _S895;
        float _S897 = _S848 * _S894;
        float _S898 = _S845 * _S897;
        float _S899 = _S847 * _S897;
        float _S900 = (_S846 * _S895 + - (_S844 * _S893 + _S850 * _S897) + _S827 * _S898 + _S851 * _S897 + _S852 * _S899) / _S853;
        float _S901 = _S855 * _S900;
        float _S902 = (_S893 + 2.0f * - _S896 + _S847 * _S899) / _S853;
        float _S903 = _s_dOut_1.z + - (_S849 * _S894 + _S855 * _S902);
        float _S904 = - _S901 + - (_S854 * - _S900 + _S848 * - _S902);
        _S829 = _S896 + _S847 * _S898;
        _S830 = _S903;
        _S831 = _S904;
        _S832 = _S896;
        _S833 = 0.0f;
        _S834 = _S901;
    }
    DiffPair_float_0 _S905;
    (&_S905)->primal_0 = _S652;
    (&_S905)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S905, _S829);
    DiffPair_float_0 _S906 = _S905;
    float _S907 = - (_S830 / _S826);
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S824;
    (&_S908)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S908, _S907);
    float _S909 = - _S908.differential_0;
    float _S910 = - (_S831 / _S823);
    DiffPair_float_0 _S911;
    (&_S911)->primal_0 = _S821;
    (&_S911)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S911, _S910);
    float _S912 = - _S911.differential_0;
    DiffPair_float_0 _S913;
    (&_S913)->primal_0 = _S649;
    (&_S913)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S913, _S832);
    DiffPair_float_0 _S914 = _S913;
    DiffPair_float_0 _S915;
    (&_S915)->primal_0 = _S648;
    (&_S915)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S915, _S833);
    DiffPair_float_0 _S916 = _S915;
    float3  _S917 = make_float3 (0.0f, 0.0f, _S834);
    if(_S768)
    {
        float _S918 = _s_dOut_1.y / _S769;
        float _S919 = _S770 * - _S918;
        float _S920 = _S771 * _S918;
        float _S921 = _S773 * _S919;
        float _S922 = _S775 * _S921;
        float _S923 = y0_5 * _S920;
        float _S924 = _S773 * _S923;
        float _S925 = _S775 * _S923;
        float _S926 = (_S774 * _S921 + - (_S772 * _S919 + _S777 * _S923) + _S759 * _S924 + _S778 * _S923 + _S779 * _S925) / _S780;
        float _S927 = x0_5 * _S926;
        float _S928 = (_S919 + 2.0f * - _S922 + _S775 * _S925) / _S780;
        float _S929 = _S776 * _S920 + x0_5 * _S928;
        float _S930 = _S922 + _S775 * _S924;
        float _S931 = _S758 * - _S926 + y0_5 * - _S928;
        _S769 = _S922;
        _S770 = _S929;
        _S771 = _S931;
        _S772 = 0.0f;
        _S773 = _S930;
        _S774 = _S927;
    }
    else
    {
        float _S932 = _s_dOut_1.y / _S781;
        float _S933 = _S782 * - _S932;
        float _S934 = _S783 * _S932;
        float _S935 = _S785 * _S933;
        float _S936 = _S787 * _S935;
        float _S937 = _S788 * _S934;
        float _S938 = _S785 * _S937;
        float _S939 = _S787 * _S937;
        float _S940 = (_S786 * _S935 + - (_S784 * _S933 + _S790 * _S937) + _S767 * _S938 + _S791 * _S937 + _S792 * _S939) / _S793;
        float _S941 = _S795 * _S940;
        float _S942 = (_S933 + 2.0f * - _S936 + _S787 * _S939) / _S793;
        float _S943 = _s_dOut_1.y + - (_S789 * _S934 + _S795 * _S942);
        float _S944 = - _S941 + - (_S794 * - _S940 + _S788 * - _S942);
        _S769 = _S936 + _S787 * _S938;
        _S770 = _S943;
        _S771 = _S944;
        _S772 = _S936;
        _S773 = 0.0f;
        _S774 = _S941;
    }
    DiffPair_float_0 _S945;
    (&_S945)->primal_0 = _S647;
    (&_S945)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S945, _S769);
    DiffPair_float_0 _S946 = _S945;
    float _S947 = - (_S770 / _S766);
    DiffPair_float_0 _S948;
    (&_S948)->primal_0 = _S764;
    (&_S948)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S948, _S947);
    float _S949 = - _S948.differential_0;
    float _S950 = - (_S771 / _S763);
    DiffPair_float_0 _S951;
    (&_S951)->primal_0 = _S761;
    (&_S951)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S951, _S950);
    float _S952 = - _S951.differential_0;
    DiffPair_float_0 _S953;
    (&_S953)->primal_0 = _S644;
    (&_S953)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S953, _S772);
    DiffPair_float_0 _S954 = _S953;
    DiffPair_float_0 _S955;
    (&_S955)->primal_0 = _S643;
    (&_S955)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S955, _S773);
    DiffPair_float_0 _S956 = _S955;
    float3  _S957 = _S917 + make_float3 (0.0f, _S774, 0.0f);
    if(_S708)
    {
        float _S958 = _s_dOut_1.x / _S709;
        float _S959 = _S710 * - _S958;
        float _S960 = _S711 * _S958;
        float _S961 = _S713 * _S959;
        float _S962 = _S715 * _S961;
        float _S963 = y0_4 * _S960;
        float _S964 = _S713 * _S963;
        float _S965 = _S715 * _S963;
        float _S966 = (_S714 * _S961 + - (_S712 * _S959 + _S717 * _S963) + _S699 * _S964 + _S718 * _S963 + _S719 * _S965) / _S720;
        float _S967 = x0_4 * _S966;
        float _S968 = (_S959 + 2.0f * - _S962 + _S715 * _S965) / _S720;
        float _S969 = _S716 * _S960 + x0_4 * _S968;
        float _S970 = _S962 + _S715 * _S964;
        float _S971 = _S698 * - _S966 + y0_4 * - _S968;
        _S709 = _S962;
        _S710 = _S969;
        _S711 = _S971;
        _S712 = 0.0f;
        _S713 = _S970;
        _S714 = _S967;
    }
    else
    {
        float _S972 = _s_dOut_1.x / _S721;
        float _S973 = _S722 * - _S972;
        float _S974 = _S723 * _S972;
        float _S975 = _S725 * _S973;
        float _S976 = _S727 * _S975;
        float _S977 = _S728 * _S974;
        float _S978 = _S725 * _S977;
        float _S979 = _S727 * _S977;
        float _S980 = (_S726 * _S975 + - (_S724 * _S973 + _S730 * _S977) + _S707 * _S978 + _S731 * _S977 + _S732 * _S979) / _S733;
        float _S981 = _S735 * _S980;
        float _S982 = (_S973 + 2.0f * - _S976 + _S727 * _S979) / _S733;
        float _S983 = _s_dOut_1.x + - (_S729 * _S974 + _S735 * _S982);
        float _S984 = - _S981 + - (_S734 * - _S980 + _S728 * - _S982);
        _S709 = _S976 + _S727 * _S978;
        _S710 = _S983;
        _S711 = _S984;
        _S712 = _S976;
        _S713 = 0.0f;
        _S714 = _S981;
    }
    DiffPair_float_0 _S985;
    (&_S985)->primal_0 = _S642;
    (&_S985)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S985, _S709);
    DiffPair_float_0 _S986 = _S985;
    float _S987 = - (_S710 / _S706);
    DiffPair_float_0 _S988;
    (&_S988)->primal_0 = _S704;
    (&_S988)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S988, _S987);
    float _S989 = - _S988.differential_0;
    float _S990 = - (_S711 / _S703);
    DiffPair_float_0 _S991;
    (&_S991)->primal_0 = _S701;
    (&_S991)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S991, _S990);
    float _S992 = - _S991.differential_0;
    DiffPair_float_0 _S993;
    (&_S993)->primal_0 = _S639;
    (&_S993)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S993, _S712);
    DiffPair_float_0 _S994 = _S993;
    DiffPair_float_0 _S995;
    (&_S995)->primal_0 = _S638;
    (&_S995)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S995, _S713);
    DiffPair_float_0 _S996 = _S995;
    float3  _S997 = _S957 + make_float3 (_S714, 0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S998;
    (&_S998)->primal_0 = _S694;
    (&_S998)->differential_0 = _S620;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S999;
    (&_S999)->primal_0 = _S695;
    (&_S999)->differential_0 = _S620;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1000;
    (&_S1000)->primal_0 = _S696;
    (&_S1000)->differential_0 = _S620;
    s_bwd_prop_clamp_0(&_S998, &_S999, &_S1000, _S997);
    float _S1001 = - _S998.differential_0.z;
    float3  s_diff_rgi_out_T_1 = make_float3 (_S998.differential_0.x + _S1001, _S998.differential_0.y + _S1001, _S998.differential_0.z);
    float3  _S1002 = _S687 * s_diff_rgi_out_T_1;
    float3  _S1003 = _S690 * s_diff_rgi_out_T_1;
    float _S1004 = (_S1002.x + _S1002.y + _S1002.z) / _S691;
    float _S1005 = intensity_3 * - _S1004;
    float _S1006 = _S689 * _S1004;
    DiffPair_float_0 _S1007;
    (&_S1007)->primal_0 = _S688;
    (&_S1007)->differential_0 = 0.0f;
    DiffPair_float_0 _S1008;
    (&_S1008)->primal_0 = 0.0f;
    (&_S1008)->differential_0 = 0.0f;
    _d_max_0(&_S1007, &_S1008, _S1005);
    float3  _S1009 = _S1003 + make_float3 (0.0f, 0.0f, _S1007.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1010;
    (&_S1010)->primal_0 = H_5;
    (&_S1010)->differential_0 = _S621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1011;
    (&_S1011)->primal_0 = rgi_in_1;
    (&_S1011)->differential_0 = _S620;
    s_bwd_prop_mul_0(&_S1010, &_S1011, _S1009);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1012 = _S1010;
    float _S1013 = _S1006 + _S1011.differential_0.z;
    float _S1014 = _S1011.differential_0.y + _S1013;
    float _S1015 = _S1011.differential_0.x + _S1013;
    float3  _S1016 = make_float3 (_S1015, _S1014, _S1013);
    if(_S680)
    {
        Matrix<float, 3, 3>  _S1017 = _S679 * _S1012.differential_0;
        Matrix<float, 3, 3>  _S1018 = _S681 * _S1012.differential_0;
        _S682 = - ((_S1017.rows[int(0)].x + _S1017.rows[int(0)].y + _S1017.rows[int(0)].z + _S1017.rows[int(1)].x + _S1017.rows[int(1)].y + _S1017.rows[int(1)].z + _S1017.rows[int(2)].x + _S1017.rows[int(2)].y + _S1017.rows[int(2)].z) / _S682);
        H_5 = _S1018;
    }
    else
    {
        _S682 = 0.0f;
        H_5 = _S1012.differential_0;
    }
    DiffPair_float_0 _S1019;
    (&_S1019)->primal_0 = _S679.rows[int(2)].z;
    (&_S1019)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1019, 0.0f);
    float _S1020 = _S1019.differential_0 + _S682;
    float3  _S1021 = _S620;
    *&((&_S1021)->z) = _S1020;
    Matrix<float, 3, 3>  _S1022 = _S621;
    _S1022[int(2)] = _S1021;
    Matrix<float, 3, 3>  _S1023 = H_5 + _S1022;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1024;
    (&_S1024)->primal_0 = _S678;
    (&_S1024)->differential_0 = _S621;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1025;
    (&_S1025)->primal_0 = S_inv_1;
    (&_S1025)->differential_0 = _S621;
    s_bwd_prop_mul_1(&_S1024, &_S1025, _S1023);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1026;
    (&_S1026)->primal_0 = T_3;
    (&_S1026)->differential_0 = _S621;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1027;
    (&_S1027)->primal_0 = D_1;
    (&_S1027)->differential_0 = _S621;
    s_bwd_prop_mul_1(&_S1026, &_S1027, _S1024.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1028 = _S1026;
    float3  _S1029 = make_float3 (_S1027.differential_0.rows[int(0)].x, _S1027.differential_0.rows[int(1)].y, _S1027.differential_0.rows[int(2)].z);
    float3  _S1030;
    if(_S673)
    {
        if(_S675)
        {
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1031;
            (&_S1031)->primal_0 = r1_3;
            (&_S1031)->differential_0 = _S620;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S1032;
            (&_S1032)->primal_0 = r2_15;
            (&_S1032)->differential_0 = _S620;
            s_bwd_prop_cross_0(&_S1031, &_S1032, _S1029);
            _S661 = _S620;
            lambda_v_7 = _S1032.differential_0;
            _S1030 = _S1031.differential_0;
        }
        else
        {
            _S661 = _S1029;
            lambda_v_7 = _S620;
            _S1030 = _S620;
        }
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1033;
        (&_S1033)->primal_0 = _S674;
        (&_S1033)->differential_0 = _S620;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1034;
        (&_S1034)->primal_0 = _S674;
        (&_S1034)->differential_0 = _S620;
        s_bwd_prop_dot_0(&_S1033, &_S1034, 0.0f);
        float3  _S1035 = _S1034.differential_0 + _S1033.differential_0 + _S661;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1036;
        (&_S1036)->primal_0 = r0_3;
        (&_S1036)->differential_0 = _S620;
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S1037;
        (&_S1037)->primal_0 = r2_15;
        (&_S1037)->differential_0 = _S620;
        s_bwd_prop_cross_0(&_S1036, &_S1037, _S1035);
        float3  _S1038 = _S1037.differential_0 + lambda_v_7;
        _S661 = _S620;
        lambda_v_7 = _S1038;
        _S674 = _S1030;
        _S1030 = _S1036.differential_0;
    }
    else
    {
        _S661 = _S1029;
        lambda_v_7 = _S620;
        _S674 = _S620;
        _S1030 = _S620;
    }
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1039;
    (&_S1039)->primal_0 = _S672;
    (&_S1039)->differential_0 = _S620;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1040;
    (&_S1040)->primal_0 = _S672;
    (&_S1040)->differential_0 = _S620;
    s_bwd_prop_dot_0(&_S1039, &_S1040, 0.0f);
    float3  _S1041 = _S1040.differential_0 + _S1039.differential_0 + _S661;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1042;
    (&_S1042)->primal_0 = r0_3;
    (&_S1042)->differential_0 = _S620;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1043;
    (&_S1043)->primal_0 = r1_3;
    (&_S1043)->differential_0 = _S620;
    s_bwd_prop_cross_0(&_S1042, &_S1043, _S1041);
    float3  _S1044 = _S620;
    *&((&_S1044)->z) = lambda_v_7.z;
    *&((&_S1044)->y) = lambda_v_7.y;
    *&((&_S1044)->x) = lambda_v_7.x;
    float3  _S1045 = _S1043.differential_0 + _S674;
    float3  _S1046 = _S620;
    *&((&_S1046)->z) = _S1045.z;
    *&((&_S1046)->y) = _S1045.y;
    *&((&_S1046)->x) = _S1045.x;
    float3  _S1047 = _S1042.differential_0 + _S1030;
    float3  _S1048 = _S620;
    *&((&_S1048)->z) = _S1047.z;
    *&((&_S1048)->y) = _S1047.y;
    *&((&_S1048)->x) = _S1047.x;
    Matrix<float, 3, 3>  _S1049 = _S621;
    _S1049[int(2)] = _S1044;
    _S1049[int(1)] = _S1046;
    _S1049[int(0)] = _S1048;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1050;
    (&_S1050)->primal_0 = skew_1;
    (&_S1050)->differential_0 = _S621;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1051;
    (&_S1051)->primal_0 = T_3;
    (&_S1051)->differential_0 = _S621;
    s_bwd_prop_mul_1(&_S1050, &_S1051, _S1049);
    Matrix<float, 3, 3>  _S1052 = _S1051.differential_0 + _S1028.differential_0;
    float2  _S1053 = make_float2 (_S1050.differential_0.rows[int(2)].y + - _S1050.differential_0.rows[int(1)].z, _S1050.differential_0.rows[int(0)].z + - _S1050.differential_0.rows[int(2)].x);
    Matrix<float, 2, 2>  _S1054 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1055;
    (&_S1055)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1055)->differential_0 = _S1054;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1056;
    (&_S1056)->primal_0 = _S664.color_params_0.n_0;
    (&_S1056)->differential_0 = _S624;
    s_bwd_prop_mul_2(&_S1055, &_S1056, _S1053);
    float2  _S1057 = make_float2 (_S1052.rows[int(0)].z, _S1052.rows[int(1)].z);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1058;
    (&_S1058)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1058)->differential_0 = _S1054;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1059;
    (&_S1059)->primal_0 = _S664.color_params_0.g_0;
    (&_S1059)->differential_0 = _S624;
    s_bwd_prop_mul_2(&_S1058, &_S1059, _S1057);
    float2  _S1060 = make_float2 (_S1052.rows[int(0)].y, _S1052.rows[int(1)].y);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1061;
    (&_S1061)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1061)->differential_0 = _S1054;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1062;
    (&_S1062)->primal_0 = _S664.color_params_0.r_0;
    (&_S1062)->differential_0 = _S624;
    s_bwd_prop_mul_2(&_S1061, &_S1062, _S1060);
    float2  _S1063 = make_float2 (_S1052.rows[int(0)].x, _S1052.rows[int(1)].x);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1064;
    (&_S1064)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1064)->differential_0 = _S1054;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1065;
    (&_S1065)->primal_0 = _S664.color_params_0.b_0;
    (&_S1065)->differential_0 = _S624;
    s_bwd_prop_mul_2(&_S1064, &_S1065, _S1063);
    ColorPPISPParams_0 _S1066 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1066)->n_0 = _S1056.differential_0;
    (&_S1066)->g_0 = _S1059.differential_0;
    (&_S1066)->r_0 = _S1062.differential_0;
    (&_S1066)->b_0 = _S1065.differential_0;
    _S661 = _S1016;
    *&((&_S661)->z) = 0.0f;
    float _S1067 = rgb_out_5.z * _S1013;
    float _S1068 = _S663 * _S1013;
    DiffPair_float_0 _S1069;
    (&_S1069)->primal_0 = falloff_5;
    (&_S1069)->differential_0 = 0.0f;
    DiffPair_float_0 _S1070;
    (&_S1070)->primal_0 = 0.0f;
    (&_S1070)->differential_0 = 0.0f;
    DiffPair_float_0 _S1071;
    (&_S1071)->primal_0 = 1.0f;
    (&_S1071)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1069, &_S1070, &_S1071, _S1067);
    float _S1072 = r2_14 * _S1069.differential_0;
    float _S1073 = r4_11 * _S1069.differential_0;
    float s_diff_r6_T_3 = _S637 * _S1069.differential_0;
    float _S1074 = r6_5 * _S1069.differential_0;
    float _S1075 = r2_14 * (_S636 * _S1069.differential_0 + r2_14 * s_diff_r6_T_3);
    float _S1076 = _S635 * _S1069.differential_0 + r4_11 * s_diff_r6_T_3 + _S1075 + _S1075;
    float _S1077 = dy_11 * _S1076;
    float _S1078 = dx_11 * _S1076;
    float _S1079 = - (_S1077 + _S1077);
    float _S1080 = - (_S1078 + _S1078);
    *&((&_S661)->y) = 0.0f;
    float _S1081 = rgb_out_5.y * _S1014;
    float _S1082 = _S662 * _S1014;
    DiffPair_float_0 _S1083;
    (&_S1083)->primal_0 = falloff_4;
    (&_S1083)->differential_0 = 0.0f;
    DiffPair_float_0 _S1084;
    (&_S1084)->primal_0 = 0.0f;
    (&_S1084)->differential_0 = 0.0f;
    DiffPair_float_0 _S1085;
    (&_S1085)->primal_0 = 1.0f;
    (&_S1085)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1083, &_S1084, &_S1085, _S1081);
    float _S1086 = r2_13 * _S1083.differential_0;
    float _S1087 = r4_10 * _S1083.differential_0;
    float s_diff_r6_T_4 = _S634 * _S1083.differential_0;
    float _S1088 = r6_4 * _S1083.differential_0;
    float _S1089 = r2_13 * (_S633 * _S1083.differential_0 + r2_13 * s_diff_r6_T_4);
    float _S1090 = _S632 * _S1083.differential_0 + r4_10 * s_diff_r6_T_4 + _S1089 + _S1089;
    float _S1091 = dy_10 * _S1090;
    float _S1092 = dx_10 * _S1090;
    float _S1093 = - (_S1091 + _S1091);
    float _S1094 = - (_S1092 + _S1092);
    *&((&_S661)->x) = 0.0f;
    float _S1095 = rgb_out_5.x * _S1015;
    float _S1096 = _S659 * _S1015;
    DiffPair_float_0 _S1097;
    (&_S1097)->primal_0 = falloff_3;
    (&_S1097)->differential_0 = 0.0f;
    DiffPair_float_0 _S1098;
    (&_S1098)->primal_0 = 0.0f;
    (&_S1098)->differential_0 = 0.0f;
    DiffPair_float_0 _S1099;
    (&_S1099)->primal_0 = 1.0f;
    (&_S1099)->differential_0 = 0.0f;
    s_bwd_prop_clamp_1(&_S1097, &_S1098, &_S1099, _S1095);
    float _S1100 = r2_12 * _S1097.differential_0;
    float _S1101 = r4_9 * _S1097.differential_0;
    float s_diff_r6_T_5 = _S631 * _S1097.differential_0;
    float _S1102 = r6_3 * _S1097.differential_0;
    float _S1103 = r2_12 * (_S630 * _S1097.differential_0 + r2_12 * s_diff_r6_T_5);
    float _S1104 = _S629 * _S1097.differential_0 + r4_9 * s_diff_r6_T_5 + _S1103 + _S1103;
    float _S1105 = dy_9 * _S1104;
    float _S1106 = dx_9 * _S1104;
    float _S1107 = - (_S1105 + _S1105);
    float _S1108 = - (_S1106 + _S1106);
    float3  _S1109 = _S620;
    *&((&_S1109)->z) = _S1068;
    *&((&_S1109)->y) = _S1082;
    *&((&_S1109)->x) = _S1096;
    float3  _S1110 = _S661 + _S1109;
    float3  _S1111 = _S619.primal_0 * _S1110;
    float3  _S1112 = _S655 * _S1110;
    float _S1113 = _S1111.x + _S1111.y + _S1111.z;
    DiffPair_float_0 _S1114;
    (&_S1114)->primal_0 = _S653.exposure_0;
    (&_S1114)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1114, _S1113);
    PPISPParamsRQS_0 _S1115 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S1115)->color_params_0 = _S1066;
    (&_S1115)->exposure_0 = _S1114.differential_0;
    _S628 = _S1115;
    (&(&_S628)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S1116 = _S1115.crf_params_0[int(2)].gc_0 + _S906.differential_0;
    (&(&_S628)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S1117 = _S1115.crf_params_0[int(2)].y0_0 + _S909;
    (&(&_S628)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S1118 = _S1115.crf_params_0[int(2)].x0_0 + _S912;
    (&(&_S628)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S1119 = _S1115.crf_params_0[int(2)].g1_0 + _S914.differential_0;
    (&(&_S628)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S1120 = _S1115.crf_params_0[int(2)].g0_0 + _S916.differential_0;
    (&(&_S628)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S1121 = _S1115.crf_params_0[int(1)].gc_0 + _S946.differential_0;
    (&(&_S628)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S1122 = _S1115.crf_params_0[int(1)].y0_0 + _S949;
    (&(&_S628)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S1123 = _S1115.crf_params_0[int(1)].x0_0 + _S952;
    (&(&_S628)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S1124 = _S1115.crf_params_0[int(1)].g1_0 + _S954.differential_0;
    (&(&_S628)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S1125 = _S1115.crf_params_0[int(1)].g0_0 + _S956.differential_0;
    (&(&_S628)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S1126 = _S1115.crf_params_0[int(0)].gc_0 + _S986.differential_0;
    (&(&_S628)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S1127 = _S1115.crf_params_0[int(0)].y0_0 + _S989;
    (&(&_S628)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S1128 = _S1115.crf_params_0[int(0)].x0_0 + _S992;
    (&(&_S628)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S1129 = _S1115.crf_params_0[int(0)].g1_0 + _S994.differential_0;
    (&(&_S628)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S1130 = _S1115.crf_params_0[int(0)].g0_0 + _S996.differential_0;
    *&((&(&(&_S628)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S628)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S628)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1131 = _S1074 + _S1115.vignette_params_0[int(2)].alpha2_0;
    (&(&_S628)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1132 = _S1073 + _S1115.vignette_params_0[int(2)].alpha1_0;
    (&(&_S628)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1133 = _S1072 + _S1115.vignette_params_0[int(2)].alpha0_0;
    (&(&_S628)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1134 = _S1079 + _S1115.vignette_params_0[int(2)].cy_0;
    (&(&_S628)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1135 = _S1080 + _S1115.vignette_params_0[int(2)].cx_0;
    (&(&_S628)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1136 = _S1088 + _S1115.vignette_params_0[int(1)].alpha2_0;
    (&(&_S628)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1137 = _S1087 + _S1115.vignette_params_0[int(1)].alpha1_0;
    (&(&_S628)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1138 = _S1086 + _S1115.vignette_params_0[int(1)].alpha0_0;
    (&(&_S628)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1139 = _S1093 + _S1115.vignette_params_0[int(1)].cy_0;
    (&(&_S628)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1140 = _S1094 + _S1115.vignette_params_0[int(1)].cx_0;
    (&(&_S628)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1141 = _S1102 + _S1115.vignette_params_0[int(0)].alpha2_0;
    (&(&_S628)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1142 = _S1101 + _S1115.vignette_params_0[int(0)].alpha1_0;
    (&(&_S628)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1143 = _S1100 + _S1115.vignette_params_0[int(0)].alpha0_0;
    (&(&_S628)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1144 = _S1107 + _S1115.vignette_params_0[int(0)].cy_0;
    (&(&_S628)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1145 = _S1108 + _S1115.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S1146;
    _S1146[int(0)] = 0.0f;
    _S1146[int(1)] = 0.0f;
    _S1146[int(2)] = 0.0f;
    _S1146[int(3)] = 0.0f;
    _S1146[int(4)] = 0.0f;
    _S1146[int(5)] = 0.0f;
    _S1146[int(6)] = 0.0f;
    _S1146[int(7)] = 0.0f;
    _S1146[int(8)] = 0.0f;
    _S1146[int(9)] = 0.0f;
    _S1146[int(10)] = 0.0f;
    _S1146[int(11)] = 0.0f;
    _S1146[int(12)] = 0.0f;
    _S1146[int(13)] = 0.0f;
    _S1146[int(14)] = 0.0f;
    _S1146[int(15)] = 0.0f;
    _S1146[int(16)] = 0.0f;
    _S1146[int(17)] = 0.0f;
    _S1146[int(18)] = 0.0f;
    _S1146[int(19)] = 0.0f;
    _S1146[int(20)] = 0.0f;
    _S1146[int(21)] = 0.0f;
    _S1146[int(22)] = 0.0f;
    _S1146[int(23)] = 0.0f;
    _S1146[int(24)] = 0.0f;
    _S1146[int(25)] = 0.0f;
    _S1146[int(26)] = 0.0f;
    _S1146[int(27)] = 0.0f;
    _S1146[int(28)] = 0.0f;
    _S1146[int(29)] = 0.0f;
    _S1146[int(30)] = 0.0f;
    _S1146[int(31)] = 0.0f;
    _S1146[int(32)] = 0.0f;
    _S1146[int(33)] = 0.0f;
    _S1146[int(34)] = 0.0f;
    _S1146[int(35)] = 0.0f;
    _S1146[int(36)] = 0.0f;
    _S1146[int(37)] = 0.0f;
    _S1146[int(38)] = 0.0f;
    _S1146[int(9)] = _S1137;
    _S1146[int(18)] = _S1115.color_params_0.r_0.x;
    _S1146[int(17)] = _S1115.color_params_0.b_0.y;
    _S1146[int(16)] = _S1115.color_params_0.b_0.x;
    _S1146[int(15)] = _S1131;
    _S1146[int(14)] = _S1132;
    _S1146[int(13)] = _S1133;
    _S1146[int(12)] = _S1134;
    _S1146[int(11)] = _S1135;
    _S1146[int(10)] = _S1136;
    _S1146[int(19)] = _S1115.color_params_0.r_0.y;
    _S1146[int(8)] = _S1138;
    _S1146[int(7)] = _S1139;
    _S1146[int(6)] = _S1140;
    _S1146[int(5)] = _S1141;
    _S1146[int(4)] = _S1142;
    _S1146[int(3)] = _S1143;
    _S1146[int(2)] = _S1144;
    _S1146[int(1)] = _S1145;
    _S1146[int(0)] = _S628.exposure_0;
    _S1146[int(28)] = _S1126;
    _S1146[int(37)] = _S1117;
    _S1146[int(36)] = _S1118;
    _S1146[int(35)] = _S1119;
    _S1146[int(34)] = _S1120;
    _S1146[int(33)] = _S1121;
    _S1146[int(32)] = _S1122;
    _S1146[int(31)] = _S1123;
    _S1146[int(30)] = _S1124;
    _S1146[int(29)] = _S1125;
    _S1146[int(38)] = _S1116;
    _S1146[int(27)] = _S1127;
    _S1146[int(26)] = _S1128;
    _S1146[int(25)] = _S1129;
    _S1146[int(24)] = _S1130;
    _S1146[int(23)] = _S1115.color_params_0.n_0.y;
    _S1146[int(22)] = _S1115.color_params_0.n_0.x;
    _S1146[int(21)] = _S1115.color_params_0.g_0.y;
    _S1146[int(20)] = _S1115.color_params_0.g_0.x;
    dpparams_1->primal_0 = dpparams_1->primal_0;
    dpparams_1->differential_0 = _S1146;
    dprgb_in_1->primal_0 = (*dprgb_in_1).primal_0;
    dprgb_in_1->differential_0 = _S1112;
    return;
}

inline __device__ void s_bwd_apply_ppisp_rqs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1147, float2  _S1148, float2  _S1149, float2  _S1150, DiffPair_arrayx3Cfloatx2C39x3E_0 * _S1151, float3  _S1152)
{
    s_bwd_prop_apply_ppisp_rqs_0(_S1147, _S1148, _S1149, _S1150, _S1151, _S1152);
    return;
}

inline __device__ void apply_ppisp_rqs_vjp(float3  rgb_in_3, float2  pix_coord_5, float2  image_center_5, float2  img_size_5, FixedArray<float, 39>  params_3, float3  grad_out_1, float3  * grad_rgb_in_1, FixedArray<float, 39>  * grad_params_1)
{
    float3  _S1153 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_in_1;
    (&dp_rgb_in_1)->primal_0 = rgb_in_3;
    (&dp_rgb_in_1)->differential_0 = _S1153;
    FixedArray<float, 39>  _S1154 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_1;
    (&dp_params_1)->primal_0 = params_3;
    (&dp_params_1)->differential_0 = _S1154;
    s_bwd_apply_ppisp_rqs_0(&dp_rgb_in_1, pix_coord_5, image_center_5, img_size_5, &dp_params_1, grad_out_1);
    *grad_rgb_in_1 = dp_rgb_in_1.differential_0;
    *grad_params_1 = (&dp_params_1)->differential_0;
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss(FixedArray<float, 36>  params_4, FixedArray<float, 22>  * _S1155)
{
    PPISPParams_0 p_2;
    (&p_2)->exposure_1 = params_4[int(0)];
    (&(&p_2)->vignette_params_1[int(0)])->cx_0 = params_4[int(1)];
    (&(&p_2)->vignette_params_1[int(0)])->cy_0 = params_4[int(2)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha0_0 = params_4[int(3)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha1_0 = params_4[int(4)];
    (&(&p_2)->vignette_params_1[int(0)])->alpha2_0 = params_4[int(5)];
    (&(&p_2)->vignette_params_1[int(1)])->cx_0 = params_4[int(6)];
    (&(&p_2)->vignette_params_1[int(1)])->cy_0 = params_4[int(7)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha0_0 = params_4[int(8)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha1_0 = params_4[int(9)];
    (&(&p_2)->vignette_params_1[int(1)])->alpha2_0 = params_4[int(10)];
    (&(&p_2)->vignette_params_1[int(2)])->cx_0 = params_4[int(11)];
    (&(&p_2)->vignette_params_1[int(2)])->cy_0 = params_4[int(12)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha0_0 = params_4[int(13)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha1_0 = params_4[int(14)];
    (&(&p_2)->vignette_params_1[int(2)])->alpha2_0 = params_4[int(15)];
    *&((&(&(&p_2)->color_params_1)->b_0)->x) = params_4[int(16)];
    *&((&(&(&p_2)->color_params_1)->b_0)->y) = params_4[int(17)];
    *&((&(&(&p_2)->color_params_1)->r_0)->x) = params_4[int(18)];
    *&((&(&(&p_2)->color_params_1)->r_0)->y) = params_4[int(19)];
    *&((&(&(&p_2)->color_params_1)->g_0)->x) = params_4[int(20)];
    *&((&(&(&p_2)->color_params_1)->g_0)->y) = params_4[int(21)];
    *&((&(&(&p_2)->color_params_1)->n_0)->x) = params_4[int(22)];
    *&((&(&(&p_2)->color_params_1)->n_0)->y) = params_4[int(23)];
    (&(&p_2)->crf_params_1[int(0)])->toe_0 = params_4[int(24)];
    (&(&p_2)->crf_params_1[int(0)])->shoulder_0 = params_4[int(25)];
    (&(&p_2)->crf_params_1[int(0)])->gamma_0 = params_4[int(26)];
    (&(&p_2)->crf_params_1[int(0)])->center_0 = params_4[int(27)];
    (&(&p_2)->crf_params_1[int(1)])->toe_0 = params_4[int(28)];
    (&(&p_2)->crf_params_1[int(1)])->shoulder_0 = params_4[int(29)];
    (&(&p_2)->crf_params_1[int(1)])->gamma_0 = params_4[int(30)];
    (&(&p_2)->crf_params_1[int(1)])->center_0 = params_4[int(31)];
    (&(&p_2)->crf_params_1[int(2)])->toe_0 = params_4[int(32)];
    (&(&p_2)->crf_params_1[int(2)])->shoulder_0 = params_4[int(33)];
    (&(&p_2)->crf_params_1[int(2)])->gamma_0 = params_4[int(34)];
    (&(&p_2)->crf_params_1[int(2)])->center_0 = params_4[int(35)];
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
    losses_0[int(0)] = p_2.exposure_1;
    float _S1156 = p_2.vignette_params_1[int(0)].cx_0;
    float _S1157 = p_2.vignette_params_1[int(0)].cy_0;
    float _S1158 = p_2.vignette_params_1[int(1)].cx_0;
    float _S1159 = p_2.vignette_params_1[int(1)].cy_0;
    float _S1160 = p_2.vignette_params_1[int(2)].cx_0;
    float _S1161 = p_2.vignette_params_1[int(2)].cy_0;
    losses_0[int(1)] = _S1156 * _S1156 + _S1157 * _S1157 + _S1158 * _S1158 + _S1159 * _S1159 + _S1160 * _S1160 + _S1161 * _S1161;
    losses_0[int(2)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha0_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha0_0)));
    losses_0[int(3)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha1_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha1_0)));
    losses_0[int(4)] = (F32_max((0.0f), (p_2.vignette_params_1[int(0)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(1)].alpha2_0))) + (F32_max((0.0f), (p_2.vignette_params_1[int(2)].alpha2_0)));
    float mean_0 = (p_2.vignette_params_1[int(0)].cx_0 + p_2.vignette_params_1[int(1)].cx_0 + p_2.vignette_params_1[int(2)].cx_0) / 3.0f;
    float _S1162 = p_2.vignette_params_1[int(0)].cx_0 - mean_0;
    float _S1163 = p_2.vignette_params_1[int(1)].cx_0 - mean_0;
    float _S1164 = p_2.vignette_params_1[int(2)].cx_0 - mean_0;
    losses_0[int(5)] = (_S1162 * _S1162 + _S1163 * _S1163 + _S1164 * _S1164) / 3.0f;
    float mean_1 = (p_2.vignette_params_1[int(0)].cy_0 + p_2.vignette_params_1[int(1)].cy_0 + p_2.vignette_params_1[int(2)].cy_0) / 3.0f;
    float _S1165 = p_2.vignette_params_1[int(0)].cy_0 - mean_1;
    float _S1166 = p_2.vignette_params_1[int(1)].cy_0 - mean_1;
    float _S1167 = p_2.vignette_params_1[int(2)].cy_0 - mean_1;
    losses_0[int(6)] = (_S1165 * _S1165 + _S1166 * _S1166 + _S1167 * _S1167) / 3.0f;
    float mean_2 = (p_2.vignette_params_1[int(0)].alpha0_0 + p_2.vignette_params_1[int(1)].alpha0_0 + p_2.vignette_params_1[int(2)].alpha0_0) / 3.0f;
    float _S1168 = p_2.vignette_params_1[int(0)].alpha0_0 - mean_2;
    float _S1169 = p_2.vignette_params_1[int(1)].alpha0_0 - mean_2;
    float _S1170 = p_2.vignette_params_1[int(2)].alpha0_0 - mean_2;
    losses_0[int(7)] = (_S1168 * _S1168 + _S1169 * _S1169 + _S1170 * _S1170) / 3.0f;
    float mean_3 = (p_2.vignette_params_1[int(0)].alpha1_0 + p_2.vignette_params_1[int(1)].alpha1_0 + p_2.vignette_params_1[int(2)].alpha1_0) / 3.0f;
    float _S1171 = p_2.vignette_params_1[int(0)].alpha1_0 - mean_3;
    float _S1172 = p_2.vignette_params_1[int(1)].alpha1_0 - mean_3;
    float _S1173 = p_2.vignette_params_1[int(2)].alpha1_0 - mean_3;
    losses_0[int(8)] = (_S1171 * _S1171 + _S1172 * _S1172 + _S1173 * _S1173) / 3.0f;
    float mean_4 = (p_2.vignette_params_1[int(0)].alpha2_0 + p_2.vignette_params_1[int(1)].alpha2_0 + p_2.vignette_params_1[int(2)].alpha2_0) / 3.0f;
    float _S1174 = p_2.vignette_params_1[int(0)].alpha2_0 - mean_4;
    float _S1175 = p_2.vignette_params_1[int(1)].alpha2_0 - mean_4;
    float _S1176 = p_2.vignette_params_1[int(2)].alpha2_0 - mean_4;
    losses_0[int(9)] = (_S1174 * _S1174 + _S1175 * _S1175 + _S1176 * _S1176) / 3.0f;
    float2  bd_2 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_2.color_params_1.b_0);
    float2  rd_2 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_2.color_params_1.r_0);
    float2  gd_2 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_2.color_params_1.g_0);
    float2  nd_2 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_2.color_params_1.n_0);
    losses_0[int(10)] = bd_2.x;
    losses_0[int(11)] = bd_2.y;
    losses_0[int(12)] = rd_2.x;
    losses_0[int(13)] = rd_2.y;
    losses_0[int(14)] = gd_2.x;
    losses_0[int(15)] = gd_2.y;
    losses_0[int(16)] = nd_2.x;
    losses_0[int(17)] = nd_2.y;
    float mean_5 = (p_2.crf_params_1[int(0)].toe_0 + p_2.crf_params_1[int(1)].toe_0 + p_2.crf_params_1[int(2)].toe_0) / 3.0f;
    float _S1177 = p_2.crf_params_1[int(0)].toe_0 - mean_5;
    float _S1178 = p_2.crf_params_1[int(1)].toe_0 - mean_5;
    float _S1179 = p_2.crf_params_1[int(2)].toe_0 - mean_5;
    losses_0[int(18)] = (_S1177 * _S1177 + _S1178 * _S1178 + _S1179 * _S1179) / 3.0f;
    float mean_6 = (p_2.crf_params_1[int(0)].shoulder_0 + p_2.crf_params_1[int(1)].shoulder_0 + p_2.crf_params_1[int(2)].shoulder_0) / 3.0f;
    float _S1180 = p_2.crf_params_1[int(0)].shoulder_0 - mean_6;
    float _S1181 = p_2.crf_params_1[int(1)].shoulder_0 - mean_6;
    float _S1182 = p_2.crf_params_1[int(2)].shoulder_0 - mean_6;
    losses_0[int(19)] = (_S1180 * _S1180 + _S1181 * _S1181 + _S1182 * _S1182) / 3.0f;
    float mean_7 = (p_2.crf_params_1[int(0)].gamma_0 + p_2.crf_params_1[int(1)].gamma_0 + p_2.crf_params_1[int(2)].gamma_0) / 3.0f;
    float _S1183 = p_2.crf_params_1[int(0)].gamma_0 - mean_7;
    float _S1184 = p_2.crf_params_1[int(1)].gamma_0 - mean_7;
    float _S1185 = p_2.crf_params_1[int(2)].gamma_0 - mean_7;
    losses_0[int(20)] = (_S1183 * _S1183 + _S1184 * _S1184 + _S1185 * _S1185) / 3.0f;
    float mean_8 = (p_2.crf_params_1[int(0)].center_0 + p_2.crf_params_1[int(1)].center_0 + p_2.crf_params_1[int(2)].center_0) / 3.0f;
    float _S1186 = p_2.crf_params_1[int(0)].center_0 - mean_8;
    float _S1187 = p_2.crf_params_1[int(1)].center_0 - mean_8;
    float _S1188 = p_2.crf_params_1[int(2)].center_0 - mean_8;
    losses_0[int(21)] = (_S1186 * _S1186 + _S1187 * _S1187 + _S1188 * _S1188) / 3.0f;
    *_S1155 = losses_0;
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss(FixedArray<float, 39>  params_5, FixedArray<float, 23>  * _S1189)
{
    PPISPParamsRQS_0 p_3;
    (&p_3)->exposure_0 = params_5[int(0)];
    (&(&p_3)->vignette_params_0[int(0)])->cx_0 = params_5[int(1)];
    (&(&p_3)->vignette_params_0[int(0)])->cy_0 = params_5[int(2)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha0_0 = params_5[int(3)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha1_0 = params_5[int(4)];
    (&(&p_3)->vignette_params_0[int(0)])->alpha2_0 = params_5[int(5)];
    (&(&p_3)->vignette_params_0[int(1)])->cx_0 = params_5[int(6)];
    (&(&p_3)->vignette_params_0[int(1)])->cy_0 = params_5[int(7)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha0_0 = params_5[int(8)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha1_0 = params_5[int(9)];
    (&(&p_3)->vignette_params_0[int(1)])->alpha2_0 = params_5[int(10)];
    (&(&p_3)->vignette_params_0[int(2)])->cx_0 = params_5[int(11)];
    (&(&p_3)->vignette_params_0[int(2)])->cy_0 = params_5[int(12)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha0_0 = params_5[int(13)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha1_0 = params_5[int(14)];
    (&(&p_3)->vignette_params_0[int(2)])->alpha2_0 = params_5[int(15)];
    *&((&(&(&p_3)->color_params_0)->b_0)->x) = params_5[int(16)];
    *&((&(&(&p_3)->color_params_0)->b_0)->y) = params_5[int(17)];
    *&((&(&(&p_3)->color_params_0)->r_0)->x) = params_5[int(18)];
    *&((&(&(&p_3)->color_params_0)->r_0)->y) = params_5[int(19)];
    *&((&(&(&p_3)->color_params_0)->g_0)->x) = params_5[int(20)];
    *&((&(&(&p_3)->color_params_0)->g_0)->y) = params_5[int(21)];
    *&((&(&(&p_3)->color_params_0)->n_0)->x) = params_5[int(22)];
    *&((&(&(&p_3)->color_params_0)->n_0)->y) = params_5[int(23)];
    (&(&p_3)->crf_params_0[int(0)])->g0_0 = params_5[int(24)];
    (&(&p_3)->crf_params_0[int(0)])->g1_0 = params_5[int(25)];
    (&(&p_3)->crf_params_0[int(0)])->x0_0 = params_5[int(26)];
    (&(&p_3)->crf_params_0[int(0)])->y0_0 = params_5[int(27)];
    (&(&p_3)->crf_params_0[int(0)])->gc_0 = params_5[int(28)];
    (&(&p_3)->crf_params_0[int(1)])->g0_0 = params_5[int(29)];
    (&(&p_3)->crf_params_0[int(1)])->g1_0 = params_5[int(30)];
    (&(&p_3)->crf_params_0[int(1)])->x0_0 = params_5[int(31)];
    (&(&p_3)->crf_params_0[int(1)])->y0_0 = params_5[int(32)];
    (&(&p_3)->crf_params_0[int(1)])->gc_0 = params_5[int(33)];
    (&(&p_3)->crf_params_0[int(2)])->g0_0 = params_5[int(34)];
    (&(&p_3)->crf_params_0[int(2)])->g1_0 = params_5[int(35)];
    (&(&p_3)->crf_params_0[int(2)])->x0_0 = params_5[int(36)];
    (&(&p_3)->crf_params_0[int(2)])->y0_0 = params_5[int(37)];
    (&(&p_3)->crf_params_0[int(2)])->gc_0 = params_5[int(38)];
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
    losses_1[int(0)] = p_3.exposure_0;
    float _S1190 = p_3.vignette_params_0[int(0)].cx_0;
    float _S1191 = p_3.vignette_params_0[int(0)].cy_0;
    float _S1192 = p_3.vignette_params_0[int(1)].cx_0;
    float _S1193 = p_3.vignette_params_0[int(1)].cy_0;
    float _S1194 = p_3.vignette_params_0[int(2)].cx_0;
    float _S1195 = p_3.vignette_params_0[int(2)].cy_0;
    losses_1[int(1)] = _S1190 * _S1190 + _S1191 * _S1191 + _S1192 * _S1192 + _S1193 * _S1193 + _S1194 * _S1194 + _S1195 * _S1195;
    losses_1[int(2)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha0_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha0_0)));
    losses_1[int(3)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha1_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha1_0)));
    losses_1[int(4)] = (F32_max((0.0f), (p_3.vignette_params_0[int(0)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(1)].alpha2_0))) + (F32_max((0.0f), (p_3.vignette_params_0[int(2)].alpha2_0)));
    float mean_9 = (p_3.vignette_params_0[int(0)].cx_0 + p_3.vignette_params_0[int(1)].cx_0 + p_3.vignette_params_0[int(2)].cx_0) / 3.0f;
    float _S1196 = p_3.vignette_params_0[int(0)].cx_0 - mean_9;
    float _S1197 = p_3.vignette_params_0[int(1)].cx_0 - mean_9;
    float _S1198 = p_3.vignette_params_0[int(2)].cx_0 - mean_9;
    losses_1[int(5)] = (_S1196 * _S1196 + _S1197 * _S1197 + _S1198 * _S1198) / 3.0f;
    float mean_10 = (p_3.vignette_params_0[int(0)].cy_0 + p_3.vignette_params_0[int(1)].cy_0 + p_3.vignette_params_0[int(2)].cy_0) / 3.0f;
    float _S1199 = p_3.vignette_params_0[int(0)].cy_0 - mean_10;
    float _S1200 = p_3.vignette_params_0[int(1)].cy_0 - mean_10;
    float _S1201 = p_3.vignette_params_0[int(2)].cy_0 - mean_10;
    losses_1[int(6)] = (_S1199 * _S1199 + _S1200 * _S1200 + _S1201 * _S1201) / 3.0f;
    float mean_11 = (p_3.vignette_params_0[int(0)].alpha0_0 + p_3.vignette_params_0[int(1)].alpha0_0 + p_3.vignette_params_0[int(2)].alpha0_0) / 3.0f;
    float _S1202 = p_3.vignette_params_0[int(0)].alpha0_0 - mean_11;
    float _S1203 = p_3.vignette_params_0[int(1)].alpha0_0 - mean_11;
    float _S1204 = p_3.vignette_params_0[int(2)].alpha0_0 - mean_11;
    losses_1[int(7)] = (_S1202 * _S1202 + _S1203 * _S1203 + _S1204 * _S1204) / 3.0f;
    float mean_12 = (p_3.vignette_params_0[int(0)].alpha1_0 + p_3.vignette_params_0[int(1)].alpha1_0 + p_3.vignette_params_0[int(2)].alpha1_0) / 3.0f;
    float _S1205 = p_3.vignette_params_0[int(0)].alpha1_0 - mean_12;
    float _S1206 = p_3.vignette_params_0[int(1)].alpha1_0 - mean_12;
    float _S1207 = p_3.vignette_params_0[int(2)].alpha1_0 - mean_12;
    losses_1[int(8)] = (_S1205 * _S1205 + _S1206 * _S1206 + _S1207 * _S1207) / 3.0f;
    float mean_13 = (p_3.vignette_params_0[int(0)].alpha2_0 + p_3.vignette_params_0[int(1)].alpha2_0 + p_3.vignette_params_0[int(2)].alpha2_0) / 3.0f;
    float _S1208 = p_3.vignette_params_0[int(0)].alpha2_0 - mean_13;
    float _S1209 = p_3.vignette_params_0[int(1)].alpha2_0 - mean_13;
    float _S1210 = p_3.vignette_params_0[int(2)].alpha2_0 - mean_13;
    losses_1[int(9)] = (_S1208 * _S1208 + _S1209 * _S1209 + _S1210 * _S1210) / 3.0f;
    float2  bd_3 = mul_0(makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f), p_3.color_params_0.b_0);
    float2  rd_3 = mul_0(makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f), p_3.color_params_0.r_0);
    float2  gd_3 = mul_0(makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f), p_3.color_params_0.g_0);
    float2  nd_3 = mul_0(makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f), p_3.color_params_0.n_0);
    losses_1[int(10)] = bd_3.x;
    losses_1[int(11)] = bd_3.y;
    losses_1[int(12)] = rd_3.x;
    losses_1[int(13)] = rd_3.y;
    losses_1[int(14)] = gd_3.x;
    losses_1[int(15)] = gd_3.y;
    losses_1[int(16)] = nd_3.x;
    losses_1[int(17)] = nd_3.y;
    float mean_14 = (p_3.crf_params_0[int(0)].g0_0 + p_3.crf_params_0[int(1)].g0_0 + p_3.crf_params_0[int(2)].g0_0) / 3.0f;
    float _S1211 = p_3.crf_params_0[int(0)].g0_0 - mean_14;
    float _S1212 = p_3.crf_params_0[int(1)].g0_0 - mean_14;
    float _S1213 = p_3.crf_params_0[int(2)].g0_0 - mean_14;
    losses_1[int(18)] = (_S1211 * _S1211 + _S1212 * _S1212 + _S1213 * _S1213) / 3.0f;
    float mean_15 = (p_3.crf_params_0[int(0)].g1_0 + p_3.crf_params_0[int(1)].g1_0 + p_3.crf_params_0[int(2)].g1_0) / 3.0f;
    float _S1214 = p_3.crf_params_0[int(0)].g1_0 - mean_15;
    float _S1215 = p_3.crf_params_0[int(1)].g1_0 - mean_15;
    float _S1216 = p_3.crf_params_0[int(2)].g1_0 - mean_15;
    losses_1[int(19)] = (_S1214 * _S1214 + _S1215 * _S1215 + _S1216 * _S1216) / 3.0f;
    float mean_16 = (p_3.crf_params_0[int(0)].x0_0 + p_3.crf_params_0[int(1)].x0_0 + p_3.crf_params_0[int(2)].x0_0) / 3.0f;
    float _S1217 = p_3.crf_params_0[int(0)].x0_0 - mean_16;
    float _S1218 = p_3.crf_params_0[int(1)].x0_0 - mean_16;
    float _S1219 = p_3.crf_params_0[int(2)].x0_0 - mean_16;
    losses_1[int(20)] = (_S1217 * _S1217 + _S1218 * _S1218 + _S1219 * _S1219) / 3.0f;
    float mean_17 = (p_3.crf_params_0[int(0)].y0_0 + p_3.crf_params_0[int(1)].y0_0 + p_3.crf_params_0[int(2)].y0_0) / 3.0f;
    float _S1220 = p_3.crf_params_0[int(0)].y0_0 - mean_17;
    float _S1221 = p_3.crf_params_0[int(1)].y0_0 - mean_17;
    float _S1222 = p_3.crf_params_0[int(2)].y0_0 - mean_17;
    losses_1[int(21)] = (_S1220 * _S1220 + _S1221 * _S1221 + _S1222 * _S1222) / 3.0f;
    float mean_18 = (p_3.crf_params_0[int(0)].gc_0 + p_3.crf_params_0[int(1)].gc_0 + p_3.crf_params_0[int(2)].gc_0) / 3.0f;
    float _S1223 = p_3.crf_params_0[int(0)].gc_0 - mean_18;
    float _S1224 = p_3.crf_params_0[int(1)].gc_0 - mean_18;
    float _S1225 = p_3.crf_params_0[int(2)].gc_0 - mean_18;
    losses_1[int(22)] = (_S1223 * _S1223 + _S1224 * _S1224 + _S1225 * _S1225) / 3.0f;
    *_S1189 = losses_1;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * dpparams_2, FixedArray<float, 22>  * _s_dOut_2)
{
    VignettingChannelParams_0 _S1226 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1227 = {
        _S1226, _S1226, _S1226
    };
    float2  _S1228 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1229 = { _S1228, _S1228, _S1228, _S1228 };
    CRFPPISPChannelParams_0 _S1230 = { 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<CRFPPISPChannelParams_0, 3>  _S1231 = {
        _S1230, _S1230, _S1230
    };
    PPISPParams_0 _S1232;
    (&_S1232)->exposure_1 = dpparams_2->primal_0[int(0)];
    (&_S1232)->vignette_params_1 = _S1227;
    (&_S1232)->color_params_1 = _S1229;
    (&_S1232)->crf_params_1 = _S1231;
    (&(&_S1232)->vignette_params_1[int(0)])->cx_0 = dpparams_2->primal_0[int(1)];
    (&(&_S1232)->vignette_params_1[int(0)])->cy_0 = dpparams_2->primal_0[int(2)];
    (&(&_S1232)->vignette_params_1[int(0)])->alpha0_0 = dpparams_2->primal_0[int(3)];
    (&(&_S1232)->vignette_params_1[int(0)])->alpha1_0 = dpparams_2->primal_0[int(4)];
    (&(&_S1232)->vignette_params_1[int(0)])->alpha2_0 = dpparams_2->primal_0[int(5)];
    (&(&_S1232)->vignette_params_1[int(1)])->cx_0 = dpparams_2->primal_0[int(6)];
    (&(&_S1232)->vignette_params_1[int(1)])->cy_0 = dpparams_2->primal_0[int(7)];
    (&(&_S1232)->vignette_params_1[int(1)])->alpha0_0 = dpparams_2->primal_0[int(8)];
    (&(&_S1232)->vignette_params_1[int(1)])->alpha1_0 = dpparams_2->primal_0[int(9)];
    (&(&_S1232)->vignette_params_1[int(1)])->alpha2_0 = dpparams_2->primal_0[int(10)];
    (&(&_S1232)->vignette_params_1[int(2)])->cx_0 = dpparams_2->primal_0[int(11)];
    (&(&_S1232)->vignette_params_1[int(2)])->cy_0 = dpparams_2->primal_0[int(12)];
    (&(&_S1232)->vignette_params_1[int(2)])->alpha0_0 = dpparams_2->primal_0[int(13)];
    (&(&_S1232)->vignette_params_1[int(2)])->alpha1_0 = dpparams_2->primal_0[int(14)];
    (&(&_S1232)->vignette_params_1[int(2)])->alpha2_0 = dpparams_2->primal_0[int(15)];
    *&((&(&(&_S1232)->color_params_1)->b_0)->x) = dpparams_2->primal_0[int(16)];
    *&((&(&(&_S1232)->color_params_1)->b_0)->y) = dpparams_2->primal_0[int(17)];
    *&((&(&(&_S1232)->color_params_1)->r_0)->x) = dpparams_2->primal_0[int(18)];
    *&((&(&(&_S1232)->color_params_1)->r_0)->y) = dpparams_2->primal_0[int(19)];
    *&((&(&(&_S1232)->color_params_1)->g_0)->x) = dpparams_2->primal_0[int(20)];
    *&((&(&(&_S1232)->color_params_1)->g_0)->y) = dpparams_2->primal_0[int(21)];
    *&((&(&(&_S1232)->color_params_1)->n_0)->x) = dpparams_2->primal_0[int(22)];
    *&((&(&(&_S1232)->color_params_1)->n_0)->y) = dpparams_2->primal_0[int(23)];
    (&(&_S1232)->crf_params_1[int(0)])->toe_0 = dpparams_2->primal_0[int(24)];
    (&(&_S1232)->crf_params_1[int(0)])->shoulder_0 = dpparams_2->primal_0[int(25)];
    (&(&_S1232)->crf_params_1[int(0)])->gamma_0 = dpparams_2->primal_0[int(26)];
    (&(&_S1232)->crf_params_1[int(0)])->center_0 = dpparams_2->primal_0[int(27)];
    (&(&_S1232)->crf_params_1[int(1)])->toe_0 = dpparams_2->primal_0[int(28)];
    (&(&_S1232)->crf_params_1[int(1)])->shoulder_0 = dpparams_2->primal_0[int(29)];
    (&(&_S1232)->crf_params_1[int(1)])->gamma_0 = dpparams_2->primal_0[int(30)];
    (&(&_S1232)->crf_params_1[int(1)])->center_0 = dpparams_2->primal_0[int(31)];
    (&(&_S1232)->crf_params_1[int(2)])->toe_0 = dpparams_2->primal_0[int(32)];
    (&(&_S1232)->crf_params_1[int(2)])->shoulder_0 = dpparams_2->primal_0[int(33)];
    (&(&_S1232)->crf_params_1[int(2)])->gamma_0 = dpparams_2->primal_0[int(34)];
    (&(&_S1232)->crf_params_1[int(2)])->center_0 = dpparams_2->primal_0[int(35)];
    float mean_19 = (dpparams_2->primal_0[int(1)] + dpparams_2->primal_0[int(6)] + dpparams_2->primal_0[int(11)]) / 3.0f;
    float _S1233 = dpparams_2->primal_0[int(1)] - mean_19;
    float _S1234 = dpparams_2->primal_0[int(6)] - mean_19;
    float _S1235 = dpparams_2->primal_0[int(11)] - mean_19;
    float mean_20 = (dpparams_2->primal_0[int(2)] + dpparams_2->primal_0[int(7)] + dpparams_2->primal_0[int(12)]) / 3.0f;
    float _S1236 = dpparams_2->primal_0[int(2)] - mean_20;
    float _S1237 = dpparams_2->primal_0[int(7)] - mean_20;
    float _S1238 = dpparams_2->primal_0[int(12)] - mean_20;
    float mean_21 = (dpparams_2->primal_0[int(3)] + dpparams_2->primal_0[int(8)] + dpparams_2->primal_0[int(13)]) / 3.0f;
    float _S1239 = dpparams_2->primal_0[int(3)] - mean_21;
    float _S1240 = dpparams_2->primal_0[int(8)] - mean_21;
    float _S1241 = dpparams_2->primal_0[int(13)] - mean_21;
    float mean_22 = (dpparams_2->primal_0[int(4)] + dpparams_2->primal_0[int(9)] + dpparams_2->primal_0[int(14)]) / 3.0f;
    float _S1242 = dpparams_2->primal_0[int(4)] - mean_22;
    float _S1243 = dpparams_2->primal_0[int(9)] - mean_22;
    float _S1244 = dpparams_2->primal_0[int(14)] - mean_22;
    float mean_23 = (dpparams_2->primal_0[int(5)] + dpparams_2->primal_0[int(10)] + dpparams_2->primal_0[int(15)]) / 3.0f;
    float _S1245 = dpparams_2->primal_0[int(5)] - mean_23;
    float _S1246 = dpparams_2->primal_0[int(10)] - mean_23;
    float _S1247 = dpparams_2->primal_0[int(15)] - mean_23;
    float mean_24 = (dpparams_2->primal_0[int(24)] + dpparams_2->primal_0[int(28)] + dpparams_2->primal_0[int(32)]) / 3.0f;
    float mean_25 = (dpparams_2->primal_0[int(25)] + dpparams_2->primal_0[int(29)] + dpparams_2->primal_0[int(33)]) / 3.0f;
    float mean_26 = (dpparams_2->primal_0[int(26)] + dpparams_2->primal_0[int(30)] + dpparams_2->primal_0[int(34)]) / 3.0f;
    float mean_27 = (dpparams_2->primal_0[int(27)] + dpparams_2->primal_0[int(31)] + dpparams_2->primal_0[int(35)]) / 3.0f;
    float _S1248 = 0.3333333432674408f * (*_s_dOut_2)[int(21)];
    float _S1249 = (dpparams_2->primal_0[int(35)] - mean_27) * _S1248;
    float _S1250 = _S1249 + _S1249;
    float _S1251 = (dpparams_2->primal_0[int(31)] - mean_27) * _S1248;
    float _S1252 = _S1251 + _S1251;
    float _S1253 = (dpparams_2->primal_0[int(27)] - mean_27) * _S1248;
    float _S1254 = _S1253 + _S1253;
    float _S1255 = 0.3333333432674408f * (- _S1250 + - _S1252 + - _S1254);
    float _S1256 = 0.3333333432674408f * (*_s_dOut_2)[int(20)];
    float _S1257 = (dpparams_2->primal_0[int(34)] - mean_26) * _S1256;
    float _S1258 = _S1257 + _S1257;
    float _S1259 = (dpparams_2->primal_0[int(30)] - mean_26) * _S1256;
    float _S1260 = _S1259 + _S1259;
    float _S1261 = (dpparams_2->primal_0[int(26)] - mean_26) * _S1256;
    float _S1262 = _S1261 + _S1261;
    float _S1263 = 0.3333333432674408f * (- _S1258 + - _S1260 + - _S1262);
    float _S1264 = 0.3333333432674408f * (*_s_dOut_2)[int(19)];
    float _S1265 = (dpparams_2->primal_0[int(33)] - mean_25) * _S1264;
    float _S1266 = _S1265 + _S1265;
    float _S1267 = (dpparams_2->primal_0[int(29)] - mean_25) * _S1264;
    float _S1268 = _S1267 + _S1267;
    float _S1269 = (dpparams_2->primal_0[int(25)] - mean_25) * _S1264;
    float _S1270 = _S1269 + _S1269;
    float _S1271 = 0.3333333432674408f * (- _S1266 + - _S1268 + - _S1270);
    float _S1272 = 0.3333333432674408f * (*_s_dOut_2)[int(18)];
    float _S1273 = (dpparams_2->primal_0[int(32)] - mean_24) * _S1272;
    float _S1274 = _S1273 + _S1273;
    float _S1275 = (dpparams_2->primal_0[int(28)] - mean_24) * _S1272;
    float _S1276 = _S1275 + _S1275;
    float _S1277 = (dpparams_2->primal_0[int(24)] - mean_24) * _S1272;
    float _S1278 = _S1277 + _S1277;
    float _S1279 = 0.3333333432674408f * (- _S1274 + - _S1276 + - _S1278);
    float2  _S1280 = make_float2 ((*_s_dOut_2)[int(16)], (*_s_dOut_2)[int(17)]);
    Matrix<float, 2, 2>  _S1281 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1282;
    (&_S1282)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1282)->differential_0 = _S1281;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1283;
    (&_S1283)->primal_0 = _S1232.color_params_1.n_0;
    (&_S1283)->differential_0 = _S1228;
    s_bwd_prop_mul_2(&_S1282, &_S1283, _S1280);
    float2  _S1284 = make_float2 ((*_s_dOut_2)[int(14)], (*_s_dOut_2)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1285;
    (&_S1285)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1285)->differential_0 = _S1281;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1286;
    (&_S1286)->primal_0 = _S1232.color_params_1.g_0;
    (&_S1286)->differential_0 = _S1228;
    s_bwd_prop_mul_2(&_S1285, &_S1286, _S1284);
    float2  _S1287 = make_float2 ((*_s_dOut_2)[int(12)], (*_s_dOut_2)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1288;
    (&_S1288)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1288)->differential_0 = _S1281;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1289;
    (&_S1289)->primal_0 = _S1232.color_params_1.r_0;
    (&_S1289)->differential_0 = _S1228;
    s_bwd_prop_mul_2(&_S1288, &_S1289, _S1287);
    float2  _S1290 = make_float2 ((*_s_dOut_2)[int(10)], (*_s_dOut_2)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1291;
    (&_S1291)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1291)->differential_0 = _S1281;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1292;
    (&_S1292)->primal_0 = _S1232.color_params_1.b_0;
    (&_S1292)->differential_0 = _S1228;
    s_bwd_prop_mul_2(&_S1291, &_S1292, _S1290);
    ColorPPISPParams_0 _S1293 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1293)->n_0 = _S1283.differential_0;
    (&_S1293)->g_0 = _S1286.differential_0;
    (&_S1293)->r_0 = _S1289.differential_0;
    (&_S1293)->b_0 = _S1292.differential_0;
    float _S1294 = 0.3333333432674408f * (*_s_dOut_2)[int(9)];
    float _S1295 = _S1247 * _S1294;
    float _S1296 = _S1295 + _S1295;
    float _S1297 = _S1246 * _S1294;
    float _S1298 = _S1297 + _S1297;
    float _S1299 = _S1245 * _S1294;
    float _S1300 = _S1299 + _S1299;
    float _S1301 = 0.3333333432674408f * (- _S1296 + - _S1298 + - _S1300);
    float _S1302 = 0.3333333432674408f * (*_s_dOut_2)[int(8)];
    float _S1303 = _S1244 * _S1302;
    float _S1304 = _S1303 + _S1303;
    float _S1305 = _S1243 * _S1302;
    float _S1306 = _S1305 + _S1305;
    float _S1307 = _S1242 * _S1302;
    float _S1308 = _S1307 + _S1307;
    float _S1309 = 0.3333333432674408f * (- _S1304 + - _S1306 + - _S1308);
    float _S1310 = 0.3333333432674408f * (*_s_dOut_2)[int(7)];
    float _S1311 = _S1241 * _S1310;
    float _S1312 = _S1311 + _S1311;
    float _S1313 = _S1240 * _S1310;
    float _S1314 = _S1313 + _S1313;
    float _S1315 = _S1239 * _S1310;
    float _S1316 = _S1315 + _S1315;
    float _S1317 = 0.3333333432674408f * (- _S1312 + - _S1314 + - _S1316);
    float _S1318 = 0.3333333432674408f * (*_s_dOut_2)[int(6)];
    float _S1319 = _S1238 * _S1318;
    float _S1320 = _S1319 + _S1319;
    float _S1321 = _S1237 * _S1318;
    float _S1322 = _S1321 + _S1321;
    float _S1323 = _S1236 * _S1318;
    float _S1324 = _S1323 + _S1323;
    float _S1325 = 0.3333333432674408f * (- _S1320 + - _S1322 + - _S1324);
    float _S1326 = 0.3333333432674408f * (*_s_dOut_2)[int(5)];
    float _S1327 = _S1235 * _S1326;
    float _S1328 = _S1327 + _S1327;
    float _S1329 = _S1234 * _S1326;
    float _S1330 = _S1329 + _S1329;
    float _S1331 = _S1233 * _S1326;
    float _S1332 = _S1331 + _S1331;
    float _S1333 = 0.3333333432674408f * (- _S1328 + - _S1330 + - _S1332);
    DiffPair_float_0 _S1334;
    (&_S1334)->primal_0 = 0.0f;
    (&_S1334)->differential_0 = 0.0f;
    DiffPair_float_0 _S1335;
    (&_S1335)->primal_0 = dpparams_2->primal_0[int(15)];
    (&_S1335)->differential_0 = 0.0f;
    _d_max_0(&_S1334, &_S1335, (*_s_dOut_2)[int(4)]);
    DiffPair_float_0 _S1336;
    (&_S1336)->primal_0 = 0.0f;
    (&_S1336)->differential_0 = 0.0f;
    DiffPair_float_0 _S1337;
    (&_S1337)->primal_0 = dpparams_2->primal_0[int(10)];
    (&_S1337)->differential_0 = 0.0f;
    _d_max_0(&_S1336, &_S1337, (*_s_dOut_2)[int(4)]);
    DiffPair_float_0 _S1338;
    (&_S1338)->primal_0 = 0.0f;
    (&_S1338)->differential_0 = 0.0f;
    DiffPair_float_0 _S1339;
    (&_S1339)->primal_0 = dpparams_2->primal_0[int(5)];
    (&_S1339)->differential_0 = 0.0f;
    _d_max_0(&_S1338, &_S1339, (*_s_dOut_2)[int(4)]);
    DiffPair_float_0 _S1340;
    (&_S1340)->primal_0 = 0.0f;
    (&_S1340)->differential_0 = 0.0f;
    DiffPair_float_0 _S1341;
    (&_S1341)->primal_0 = dpparams_2->primal_0[int(14)];
    (&_S1341)->differential_0 = 0.0f;
    _d_max_0(&_S1340, &_S1341, (*_s_dOut_2)[int(3)]);
    DiffPair_float_0 _S1342;
    (&_S1342)->primal_0 = 0.0f;
    (&_S1342)->differential_0 = 0.0f;
    DiffPair_float_0 _S1343;
    (&_S1343)->primal_0 = dpparams_2->primal_0[int(9)];
    (&_S1343)->differential_0 = 0.0f;
    _d_max_0(&_S1342, &_S1343, (*_s_dOut_2)[int(3)]);
    DiffPair_float_0 _S1344;
    (&_S1344)->primal_0 = 0.0f;
    (&_S1344)->differential_0 = 0.0f;
    DiffPair_float_0 _S1345;
    (&_S1345)->primal_0 = dpparams_2->primal_0[int(4)];
    (&_S1345)->differential_0 = 0.0f;
    _d_max_0(&_S1344, &_S1345, (*_s_dOut_2)[int(3)]);
    DiffPair_float_0 _S1346;
    (&_S1346)->primal_0 = 0.0f;
    (&_S1346)->differential_0 = 0.0f;
    DiffPair_float_0 _S1347;
    (&_S1347)->primal_0 = dpparams_2->primal_0[int(13)];
    (&_S1347)->differential_0 = 0.0f;
    _d_max_0(&_S1346, &_S1347, (*_s_dOut_2)[int(2)]);
    DiffPair_float_0 _S1348;
    (&_S1348)->primal_0 = 0.0f;
    (&_S1348)->differential_0 = 0.0f;
    DiffPair_float_0 _S1349;
    (&_S1349)->primal_0 = dpparams_2->primal_0[int(8)];
    (&_S1349)->differential_0 = 0.0f;
    _d_max_0(&_S1348, &_S1349, (*_s_dOut_2)[int(2)]);
    DiffPair_float_0 _S1350;
    (&_S1350)->primal_0 = 0.0f;
    (&_S1350)->differential_0 = 0.0f;
    DiffPair_float_0 _S1351;
    (&_S1351)->primal_0 = dpparams_2->primal_0[int(3)];
    (&_S1351)->differential_0 = 0.0f;
    _d_max_0(&_S1350, &_S1351, (*_s_dOut_2)[int(2)]);
    float _S1352 = dpparams_2->primal_0[int(12)] * (*_s_dOut_2)[int(1)];
    float _S1353 = dpparams_2->primal_0[int(11)] * (*_s_dOut_2)[int(1)];
    float _S1354 = dpparams_2->primal_0[int(7)] * (*_s_dOut_2)[int(1)];
    float _S1355 = dpparams_2->primal_0[int(6)] * (*_s_dOut_2)[int(1)];
    float _S1356 = dpparams_2->primal_0[int(2)] * (*_s_dOut_2)[int(1)];
    float _S1357 = dpparams_2->primal_0[int(1)] * (*_s_dOut_2)[int(1)];
    PPISPParams_0 _S1358 = PPISPParams_x24_syn_dzero_0();
    (&_S1358)->color_params_1 = _S1293;
    (&_S1358)->exposure_1 = (*_s_dOut_2)[int(0)];
    _S1232 = _S1358;
    (&(&_S1232)->crf_params_1[int(2)])->center_0 = 0.0f;
    float _S1359 = _S1250 + _S1255 + _S1358.crf_params_1[int(2)].center_0;
    (&(&_S1232)->crf_params_1[int(2)])->gamma_0 = 0.0f;
    float _S1360 = _S1258 + _S1263 + _S1358.crf_params_1[int(2)].gamma_0;
    (&(&_S1232)->crf_params_1[int(2)])->shoulder_0 = 0.0f;
    float _S1361 = _S1266 + _S1271 + _S1358.crf_params_1[int(2)].shoulder_0;
    (&(&_S1232)->crf_params_1[int(2)])->toe_0 = 0.0f;
    float _S1362 = _S1274 + _S1279 + _S1358.crf_params_1[int(2)].toe_0;
    (&(&_S1232)->crf_params_1[int(1)])->center_0 = 0.0f;
    float _S1363 = _S1252 + _S1255 + _S1358.crf_params_1[int(1)].center_0;
    (&(&_S1232)->crf_params_1[int(1)])->gamma_0 = 0.0f;
    float _S1364 = _S1260 + _S1263 + _S1358.crf_params_1[int(1)].gamma_0;
    (&(&_S1232)->crf_params_1[int(1)])->shoulder_0 = 0.0f;
    float _S1365 = _S1268 + _S1271 + _S1358.crf_params_1[int(1)].shoulder_0;
    (&(&_S1232)->crf_params_1[int(1)])->toe_0 = 0.0f;
    float _S1366 = _S1276 + _S1279 + _S1358.crf_params_1[int(1)].toe_0;
    (&(&_S1232)->crf_params_1[int(0)])->center_0 = 0.0f;
    float _S1367 = _S1254 + _S1255 + _S1358.crf_params_1[int(0)].center_0;
    (&(&_S1232)->crf_params_1[int(0)])->gamma_0 = 0.0f;
    float _S1368 = _S1262 + _S1263 + _S1358.crf_params_1[int(0)].gamma_0;
    (&(&_S1232)->crf_params_1[int(0)])->shoulder_0 = 0.0f;
    float _S1369 = _S1270 + _S1271 + _S1358.crf_params_1[int(0)].shoulder_0;
    (&(&_S1232)->crf_params_1[int(0)])->toe_0 = 0.0f;
    float _S1370 = _S1278 + _S1279 + _S1358.crf_params_1[int(0)].toe_0;
    *&((&(&(&_S1232)->color_params_1)->n_0)->y) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->n_0)->x) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->g_0)->y) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->g_0)->x) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->r_0)->y) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->r_0)->x) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->b_0)->y) = 0.0f;
    *&((&(&(&_S1232)->color_params_1)->b_0)->x) = 0.0f;
    (&(&_S1232)->vignette_params_1[int(2)])->alpha2_0 = 0.0f;
    float _S1371 = _S1296 + _S1301 + _S1335.differential_0 + _S1358.vignette_params_1[int(2)].alpha2_0;
    (&(&_S1232)->vignette_params_1[int(2)])->alpha1_0 = 0.0f;
    float _S1372 = _S1304 + _S1309 + _S1341.differential_0 + _S1358.vignette_params_1[int(2)].alpha1_0;
    (&(&_S1232)->vignette_params_1[int(2)])->alpha0_0 = 0.0f;
    float _S1373 = _S1312 + _S1317 + _S1347.differential_0 + _S1358.vignette_params_1[int(2)].alpha0_0;
    (&(&_S1232)->vignette_params_1[int(2)])->cy_0 = 0.0f;
    float _S1374 = _S1320 + _S1325 + _S1352 + _S1352 + _S1358.vignette_params_1[int(2)].cy_0;
    (&(&_S1232)->vignette_params_1[int(2)])->cx_0 = 0.0f;
    float _S1375 = _S1328 + _S1333 + _S1353 + _S1353 + _S1358.vignette_params_1[int(2)].cx_0;
    (&(&_S1232)->vignette_params_1[int(1)])->alpha2_0 = 0.0f;
    float _S1376 = _S1298 + _S1301 + _S1337.differential_0 + _S1358.vignette_params_1[int(1)].alpha2_0;
    (&(&_S1232)->vignette_params_1[int(1)])->alpha1_0 = 0.0f;
    float _S1377 = _S1306 + _S1309 + _S1343.differential_0 + _S1358.vignette_params_1[int(1)].alpha1_0;
    (&(&_S1232)->vignette_params_1[int(1)])->alpha0_0 = 0.0f;
    float _S1378 = _S1314 + _S1317 + _S1349.differential_0 + _S1358.vignette_params_1[int(1)].alpha0_0;
    (&(&_S1232)->vignette_params_1[int(1)])->cy_0 = 0.0f;
    float _S1379 = _S1322 + _S1325 + _S1354 + _S1354 + _S1358.vignette_params_1[int(1)].cy_0;
    (&(&_S1232)->vignette_params_1[int(1)])->cx_0 = 0.0f;
    float _S1380 = _S1330 + _S1333 + _S1355 + _S1355 + _S1358.vignette_params_1[int(1)].cx_0;
    (&(&_S1232)->vignette_params_1[int(0)])->alpha2_0 = 0.0f;
    float _S1381 = _S1300 + _S1301 + _S1339.differential_0 + _S1358.vignette_params_1[int(0)].alpha2_0;
    (&(&_S1232)->vignette_params_1[int(0)])->alpha1_0 = 0.0f;
    float _S1382 = _S1308 + _S1309 + _S1345.differential_0 + _S1358.vignette_params_1[int(0)].alpha1_0;
    (&(&_S1232)->vignette_params_1[int(0)])->alpha0_0 = 0.0f;
    float _S1383 = _S1316 + _S1317 + _S1351.differential_0 + _S1358.vignette_params_1[int(0)].alpha0_0;
    (&(&_S1232)->vignette_params_1[int(0)])->cy_0 = 0.0f;
    float _S1384 = _S1324 + _S1325 + _S1356 + _S1356 + _S1358.vignette_params_1[int(0)].cy_0;
    (&(&_S1232)->vignette_params_1[int(0)])->cx_0 = 0.0f;
    float _S1385 = _S1332 + _S1333 + _S1357 + _S1357 + _S1358.vignette_params_1[int(0)].cx_0;
    FixedArray<float, 36>  _S1386;
    _S1386[int(0)] = 0.0f;
    _S1386[int(1)] = 0.0f;
    _S1386[int(2)] = 0.0f;
    _S1386[int(3)] = 0.0f;
    _S1386[int(4)] = 0.0f;
    _S1386[int(5)] = 0.0f;
    _S1386[int(6)] = 0.0f;
    _S1386[int(7)] = 0.0f;
    _S1386[int(8)] = 0.0f;
    _S1386[int(9)] = 0.0f;
    _S1386[int(10)] = 0.0f;
    _S1386[int(11)] = 0.0f;
    _S1386[int(12)] = 0.0f;
    _S1386[int(13)] = 0.0f;
    _S1386[int(14)] = 0.0f;
    _S1386[int(15)] = 0.0f;
    _S1386[int(16)] = 0.0f;
    _S1386[int(17)] = 0.0f;
    _S1386[int(18)] = 0.0f;
    _S1386[int(19)] = 0.0f;
    _S1386[int(20)] = 0.0f;
    _S1386[int(21)] = 0.0f;
    _S1386[int(22)] = 0.0f;
    _S1386[int(23)] = 0.0f;
    _S1386[int(24)] = 0.0f;
    _S1386[int(25)] = 0.0f;
    _S1386[int(26)] = 0.0f;
    _S1386[int(27)] = 0.0f;
    _S1386[int(28)] = 0.0f;
    _S1386[int(29)] = 0.0f;
    _S1386[int(30)] = 0.0f;
    _S1386[int(31)] = 0.0f;
    _S1386[int(32)] = 0.0f;
    _S1386[int(33)] = 0.0f;
    _S1386[int(34)] = 0.0f;
    _S1386[int(35)] = 0.0f;
    _S1386[int(8)] = _S1378;
    _S1386[int(16)] = _S1358.color_params_1.b_0.x;
    _S1386[int(15)] = _S1371;
    _S1386[int(14)] = _S1372;
    _S1386[int(13)] = _S1373;
    _S1386[int(12)] = _S1374;
    _S1386[int(11)] = _S1375;
    _S1386[int(10)] = _S1376;
    _S1386[int(9)] = _S1377;
    _S1386[int(17)] = _S1358.color_params_1.b_0.y;
    _S1386[int(7)] = _S1379;
    _S1386[int(6)] = _S1380;
    _S1386[int(5)] = _S1381;
    _S1386[int(4)] = _S1382;
    _S1386[int(3)] = _S1383;
    _S1386[int(2)] = _S1384;
    _S1386[int(1)] = _S1385;
    _S1386[int(0)] = _S1232.exposure_1;
    _S1386[int(26)] = _S1368;
    _S1386[int(34)] = _S1360;
    _S1386[int(33)] = _S1361;
    _S1386[int(32)] = _S1362;
    _S1386[int(31)] = _S1363;
    _S1386[int(30)] = _S1364;
    _S1386[int(29)] = _S1365;
    _S1386[int(28)] = _S1366;
    _S1386[int(27)] = _S1367;
    _S1386[int(35)] = _S1359;
    _S1386[int(25)] = _S1369;
    _S1386[int(24)] = _S1370;
    _S1386[int(23)] = _S1358.color_params_1.n_0.y;
    _S1386[int(22)] = _S1358.color_params_1.n_0.x;
    _S1386[int(21)] = _S1358.color_params_1.g_0.y;
    _S1386[int(20)] = _S1358.color_params_1.g_0.x;
    _S1386[int(19)] = _S1358.color_params_1.r_0.y;
    _S1386[int(18)] = _S1358.color_params_1.r_0.x;
    dpparams_2->primal_0 = dpparams_2->primal_0;
    dpparams_2->differential_0 = _S1386;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C36x3E_0 * _S1387, FixedArray<float, 22>  * _S1388)
{
    s_bwd_prop_compute_raw_ppisp_regularization_loss_0(_S1387, _S1388);
    return;
}

inline __device__ void compute_raw_ppisp_regularization_loss_vjp(FixedArray<float, 36>  params_6, FixedArray<float, 22>  grad_out_2, FixedArray<float, 36>  * _S1389)
{
    FixedArray<float, 36>  _S1390 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C36x3E_0 dp_params_2;
    (&dp_params_2)->primal_0 = params_6;
    (&dp_params_2)->differential_0 = _S1390;
    FixedArray<float, 22>  _S1391 = grad_out_2;
    s_bwd_compute_raw_ppisp_regularization_loss_0(&dp_params_2, &_S1391);
    *_S1389 = (&dp_params_2)->differential_0;
    return;
}

inline __device__ void s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * dpparams_3, FixedArray<float, 23>  * _s_dOut_3)
{
    VignettingChannelParams_0 _S1392 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<VignettingChannelParams_0, 3>  _S1393 = {
        _S1392, _S1392, _S1392
    };
    float2  _S1394 = make_float2 (0.0f);
    ColorPPISPParams_0 _S1395 = { _S1394, _S1394, _S1394, _S1394 };
    RQSCRFPPISPChannelParams_0 _S1396 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    FixedArray<RQSCRFPPISPChannelParams_0, 3>  _S1397 = {
        _S1396, _S1396, _S1396
    };
    PPISPParamsRQS_0 _S1398;
    (&_S1398)->exposure_0 = dpparams_3->primal_0[int(0)];
    (&_S1398)->vignette_params_0 = _S1393;
    (&_S1398)->color_params_0 = _S1395;
    (&_S1398)->crf_params_0 = _S1397;
    (&(&_S1398)->vignette_params_0[int(0)])->cx_0 = dpparams_3->primal_0[int(1)];
    (&(&_S1398)->vignette_params_0[int(0)])->cy_0 = dpparams_3->primal_0[int(2)];
    (&(&_S1398)->vignette_params_0[int(0)])->alpha0_0 = dpparams_3->primal_0[int(3)];
    (&(&_S1398)->vignette_params_0[int(0)])->alpha1_0 = dpparams_3->primal_0[int(4)];
    (&(&_S1398)->vignette_params_0[int(0)])->alpha2_0 = dpparams_3->primal_0[int(5)];
    (&(&_S1398)->vignette_params_0[int(1)])->cx_0 = dpparams_3->primal_0[int(6)];
    (&(&_S1398)->vignette_params_0[int(1)])->cy_0 = dpparams_3->primal_0[int(7)];
    (&(&_S1398)->vignette_params_0[int(1)])->alpha0_0 = dpparams_3->primal_0[int(8)];
    (&(&_S1398)->vignette_params_0[int(1)])->alpha1_0 = dpparams_3->primal_0[int(9)];
    (&(&_S1398)->vignette_params_0[int(1)])->alpha2_0 = dpparams_3->primal_0[int(10)];
    (&(&_S1398)->vignette_params_0[int(2)])->cx_0 = dpparams_3->primal_0[int(11)];
    (&(&_S1398)->vignette_params_0[int(2)])->cy_0 = dpparams_3->primal_0[int(12)];
    (&(&_S1398)->vignette_params_0[int(2)])->alpha0_0 = dpparams_3->primal_0[int(13)];
    (&(&_S1398)->vignette_params_0[int(2)])->alpha1_0 = dpparams_3->primal_0[int(14)];
    (&(&_S1398)->vignette_params_0[int(2)])->alpha2_0 = dpparams_3->primal_0[int(15)];
    *&((&(&(&_S1398)->color_params_0)->b_0)->x) = dpparams_3->primal_0[int(16)];
    *&((&(&(&_S1398)->color_params_0)->b_0)->y) = dpparams_3->primal_0[int(17)];
    *&((&(&(&_S1398)->color_params_0)->r_0)->x) = dpparams_3->primal_0[int(18)];
    *&((&(&(&_S1398)->color_params_0)->r_0)->y) = dpparams_3->primal_0[int(19)];
    *&((&(&(&_S1398)->color_params_0)->g_0)->x) = dpparams_3->primal_0[int(20)];
    *&((&(&(&_S1398)->color_params_0)->g_0)->y) = dpparams_3->primal_0[int(21)];
    *&((&(&(&_S1398)->color_params_0)->n_0)->x) = dpparams_3->primal_0[int(22)];
    *&((&(&(&_S1398)->color_params_0)->n_0)->y) = dpparams_3->primal_0[int(23)];
    (&(&_S1398)->crf_params_0[int(0)])->g0_0 = dpparams_3->primal_0[int(24)];
    (&(&_S1398)->crf_params_0[int(0)])->g1_0 = dpparams_3->primal_0[int(25)];
    (&(&_S1398)->crf_params_0[int(0)])->x0_0 = dpparams_3->primal_0[int(26)];
    (&(&_S1398)->crf_params_0[int(0)])->y0_0 = dpparams_3->primal_0[int(27)];
    (&(&_S1398)->crf_params_0[int(0)])->gc_0 = dpparams_3->primal_0[int(28)];
    (&(&_S1398)->crf_params_0[int(1)])->g0_0 = dpparams_3->primal_0[int(29)];
    (&(&_S1398)->crf_params_0[int(1)])->g1_0 = dpparams_3->primal_0[int(30)];
    (&(&_S1398)->crf_params_0[int(1)])->x0_0 = dpparams_3->primal_0[int(31)];
    (&(&_S1398)->crf_params_0[int(1)])->y0_0 = dpparams_3->primal_0[int(32)];
    (&(&_S1398)->crf_params_0[int(1)])->gc_0 = dpparams_3->primal_0[int(33)];
    (&(&_S1398)->crf_params_0[int(2)])->g0_0 = dpparams_3->primal_0[int(34)];
    (&(&_S1398)->crf_params_0[int(2)])->g1_0 = dpparams_3->primal_0[int(35)];
    (&(&_S1398)->crf_params_0[int(2)])->x0_0 = dpparams_3->primal_0[int(36)];
    (&(&_S1398)->crf_params_0[int(2)])->y0_0 = dpparams_3->primal_0[int(37)];
    (&(&_S1398)->crf_params_0[int(2)])->gc_0 = dpparams_3->primal_0[int(38)];
    float mean_28 = (dpparams_3->primal_0[int(1)] + dpparams_3->primal_0[int(6)] + dpparams_3->primal_0[int(11)]) / 3.0f;
    float _S1399 = dpparams_3->primal_0[int(1)] - mean_28;
    float _S1400 = dpparams_3->primal_0[int(6)] - mean_28;
    float _S1401 = dpparams_3->primal_0[int(11)] - mean_28;
    float mean_29 = (dpparams_3->primal_0[int(2)] + dpparams_3->primal_0[int(7)] + dpparams_3->primal_0[int(12)]) / 3.0f;
    float _S1402 = dpparams_3->primal_0[int(2)] - mean_29;
    float _S1403 = dpparams_3->primal_0[int(7)] - mean_29;
    float _S1404 = dpparams_3->primal_0[int(12)] - mean_29;
    float mean_30 = (dpparams_3->primal_0[int(3)] + dpparams_3->primal_0[int(8)] + dpparams_3->primal_0[int(13)]) / 3.0f;
    float _S1405 = dpparams_3->primal_0[int(3)] - mean_30;
    float _S1406 = dpparams_3->primal_0[int(8)] - mean_30;
    float _S1407 = dpparams_3->primal_0[int(13)] - mean_30;
    float mean_31 = (dpparams_3->primal_0[int(4)] + dpparams_3->primal_0[int(9)] + dpparams_3->primal_0[int(14)]) / 3.0f;
    float _S1408 = dpparams_3->primal_0[int(4)] - mean_31;
    float _S1409 = dpparams_3->primal_0[int(9)] - mean_31;
    float _S1410 = dpparams_3->primal_0[int(14)] - mean_31;
    float mean_32 = (dpparams_3->primal_0[int(5)] + dpparams_3->primal_0[int(10)] + dpparams_3->primal_0[int(15)]) / 3.0f;
    float _S1411 = dpparams_3->primal_0[int(5)] - mean_32;
    float _S1412 = dpparams_3->primal_0[int(10)] - mean_32;
    float _S1413 = dpparams_3->primal_0[int(15)] - mean_32;
    float mean_33 = (dpparams_3->primal_0[int(24)] + dpparams_3->primal_0[int(29)] + dpparams_3->primal_0[int(34)]) / 3.0f;
    float mean_34 = (dpparams_3->primal_0[int(25)] + dpparams_3->primal_0[int(30)] + dpparams_3->primal_0[int(35)]) / 3.0f;
    float mean_35 = (dpparams_3->primal_0[int(26)] + dpparams_3->primal_0[int(31)] + dpparams_3->primal_0[int(36)]) / 3.0f;
    float mean_36 = (dpparams_3->primal_0[int(27)] + dpparams_3->primal_0[int(32)] + dpparams_3->primal_0[int(37)]) / 3.0f;
    float mean_37 = (dpparams_3->primal_0[int(28)] + dpparams_3->primal_0[int(33)] + dpparams_3->primal_0[int(38)]) / 3.0f;
    float _S1414 = 0.3333333432674408f * (*_s_dOut_3)[int(22)];
    float _S1415 = (dpparams_3->primal_0[int(38)] - mean_37) * _S1414;
    float _S1416 = _S1415 + _S1415;
    float _S1417 = (dpparams_3->primal_0[int(33)] - mean_37) * _S1414;
    float _S1418 = _S1417 + _S1417;
    float _S1419 = (dpparams_3->primal_0[int(28)] - mean_37) * _S1414;
    float _S1420 = _S1419 + _S1419;
    float _S1421 = 0.3333333432674408f * (- _S1416 + - _S1418 + - _S1420);
    float _S1422 = 0.3333333432674408f * (*_s_dOut_3)[int(21)];
    float _S1423 = (dpparams_3->primal_0[int(37)] - mean_36) * _S1422;
    float _S1424 = _S1423 + _S1423;
    float _S1425 = (dpparams_3->primal_0[int(32)] - mean_36) * _S1422;
    float _S1426 = _S1425 + _S1425;
    float _S1427 = (dpparams_3->primal_0[int(27)] - mean_36) * _S1422;
    float _S1428 = _S1427 + _S1427;
    float _S1429 = 0.3333333432674408f * (- _S1424 + - _S1426 + - _S1428);
    float _S1430 = 0.3333333432674408f * (*_s_dOut_3)[int(20)];
    float _S1431 = (dpparams_3->primal_0[int(36)] - mean_35) * _S1430;
    float _S1432 = _S1431 + _S1431;
    float _S1433 = (dpparams_3->primal_0[int(31)] - mean_35) * _S1430;
    float _S1434 = _S1433 + _S1433;
    float _S1435 = (dpparams_3->primal_0[int(26)] - mean_35) * _S1430;
    float _S1436 = _S1435 + _S1435;
    float _S1437 = 0.3333333432674408f * (- _S1432 + - _S1434 + - _S1436);
    float _S1438 = 0.3333333432674408f * (*_s_dOut_3)[int(19)];
    float _S1439 = (dpparams_3->primal_0[int(35)] - mean_34) * _S1438;
    float _S1440 = _S1439 + _S1439;
    float _S1441 = (dpparams_3->primal_0[int(30)] - mean_34) * _S1438;
    float _S1442 = _S1441 + _S1441;
    float _S1443 = (dpparams_3->primal_0[int(25)] - mean_34) * _S1438;
    float _S1444 = _S1443 + _S1443;
    float _S1445 = 0.3333333432674408f * (- _S1440 + - _S1442 + - _S1444);
    float _S1446 = 0.3333333432674408f * (*_s_dOut_3)[int(18)];
    float _S1447 = (dpparams_3->primal_0[int(34)] - mean_33) * _S1446;
    float _S1448 = _S1447 + _S1447;
    float _S1449 = (dpparams_3->primal_0[int(29)] - mean_33) * _S1446;
    float _S1450 = _S1449 + _S1449;
    float _S1451 = (dpparams_3->primal_0[int(24)] - mean_33) * _S1446;
    float _S1452 = _S1451 + _S1451;
    float _S1453 = 0.3333333432674408f * (- _S1448 + - _S1450 + - _S1452);
    float2  _S1454 = make_float2 ((*_s_dOut_3)[int(16)], (*_s_dOut_3)[int(17)]);
    Matrix<float, 2, 2>  _S1455 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1456;
    (&_S1456)->primal_0 = makeMatrix<float, 2, 2> (0.01283689960837364f, -0.00346540007740259f, -0.00346540007740259f, 0.01281579956412315f);
    (&_S1456)->differential_0 = _S1455;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1457;
    (&_S1457)->primal_0 = _S1398.color_params_0.n_0;
    (&_S1457)->differential_0 = _S1394;
    s_bwd_prop_mul_2(&_S1456, &_S1457, _S1454);
    float2  _S1458 = make_float2 ((*_s_dOut_3)[int(14)], (*_s_dOut_3)[int(15)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1459;
    (&_S1459)->primal_0 = makeMatrix<float, 2, 2> (0.04333360120654106f, -0.01805369928479195f, -0.01805369928479195f, 0.0580499991774559f);
    (&_S1459)->differential_0 = _S1455;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1460;
    (&_S1460)->primal_0 = _S1398.color_params_0.g_0;
    (&_S1460)->differential_0 = _S1394;
    s_bwd_prop_mul_2(&_S1459, &_S1460, _S1458);
    float2  _S1461 = make_float2 ((*_s_dOut_3)[int(12)], (*_s_dOut_3)[int(13)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1462;
    (&_S1462)->primal_0 = makeMatrix<float, 2, 2> (0.05805699899792671f, -0.0179871991276741f, -0.0179871991276741f, 0.04310610145330429f);
    (&_S1462)->differential_0 = _S1455;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1463;
    (&_S1463)->primal_0 = _S1398.color_params_0.r_0;
    (&_S1463)->differential_0 = _S1394;
    s_bwd_prop_mul_2(&_S1462, &_S1463, _S1461);
    float2  _S1464 = make_float2 ((*_s_dOut_3)[int(10)], (*_s_dOut_3)[int(11)]);
    DiffPair_matrixx3Cfloatx2C2x2C2x3E_0 _S1465;
    (&_S1465)->primal_0 = makeMatrix<float, 2, 2> (0.04805419966578484f, -0.0043631000444293f, -0.0043631000444293f, 0.04812829941511154f);
    (&_S1465)->differential_0 = _S1455;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1466;
    (&_S1466)->primal_0 = _S1398.color_params_0.b_0;
    (&_S1466)->differential_0 = _S1394;
    s_bwd_prop_mul_2(&_S1465, &_S1466, _S1464);
    ColorPPISPParams_0 _S1467 = ColorPPISPParams_x24_syn_dzero_0();
    (&_S1467)->n_0 = _S1457.differential_0;
    (&_S1467)->g_0 = _S1460.differential_0;
    (&_S1467)->r_0 = _S1463.differential_0;
    (&_S1467)->b_0 = _S1466.differential_0;
    float _S1468 = 0.3333333432674408f * (*_s_dOut_3)[int(9)];
    float _S1469 = _S1413 * _S1468;
    float _S1470 = _S1469 + _S1469;
    float _S1471 = _S1412 * _S1468;
    float _S1472 = _S1471 + _S1471;
    float _S1473 = _S1411 * _S1468;
    float _S1474 = _S1473 + _S1473;
    float _S1475 = 0.3333333432674408f * (- _S1470 + - _S1472 + - _S1474);
    float _S1476 = 0.3333333432674408f * (*_s_dOut_3)[int(8)];
    float _S1477 = _S1410 * _S1476;
    float _S1478 = _S1477 + _S1477;
    float _S1479 = _S1409 * _S1476;
    float _S1480 = _S1479 + _S1479;
    float _S1481 = _S1408 * _S1476;
    float _S1482 = _S1481 + _S1481;
    float _S1483 = 0.3333333432674408f * (- _S1478 + - _S1480 + - _S1482);
    float _S1484 = 0.3333333432674408f * (*_s_dOut_3)[int(7)];
    float _S1485 = _S1407 * _S1484;
    float _S1486 = _S1485 + _S1485;
    float _S1487 = _S1406 * _S1484;
    float _S1488 = _S1487 + _S1487;
    float _S1489 = _S1405 * _S1484;
    float _S1490 = _S1489 + _S1489;
    float _S1491 = 0.3333333432674408f * (- _S1486 + - _S1488 + - _S1490);
    float _S1492 = 0.3333333432674408f * (*_s_dOut_3)[int(6)];
    float _S1493 = _S1404 * _S1492;
    float _S1494 = _S1493 + _S1493;
    float _S1495 = _S1403 * _S1492;
    float _S1496 = _S1495 + _S1495;
    float _S1497 = _S1402 * _S1492;
    float _S1498 = _S1497 + _S1497;
    float _S1499 = 0.3333333432674408f * (- _S1494 + - _S1496 + - _S1498);
    float _S1500 = 0.3333333432674408f * (*_s_dOut_3)[int(5)];
    float _S1501 = _S1401 * _S1500;
    float _S1502 = _S1501 + _S1501;
    float _S1503 = _S1400 * _S1500;
    float _S1504 = _S1503 + _S1503;
    float _S1505 = _S1399 * _S1500;
    float _S1506 = _S1505 + _S1505;
    float _S1507 = 0.3333333432674408f * (- _S1502 + - _S1504 + - _S1506);
    DiffPair_float_0 _S1508;
    (&_S1508)->primal_0 = 0.0f;
    (&_S1508)->differential_0 = 0.0f;
    DiffPair_float_0 _S1509;
    (&_S1509)->primal_0 = dpparams_3->primal_0[int(15)];
    (&_S1509)->differential_0 = 0.0f;
    _d_max_0(&_S1508, &_S1509, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1510;
    (&_S1510)->primal_0 = 0.0f;
    (&_S1510)->differential_0 = 0.0f;
    DiffPair_float_0 _S1511;
    (&_S1511)->primal_0 = dpparams_3->primal_0[int(10)];
    (&_S1511)->differential_0 = 0.0f;
    _d_max_0(&_S1510, &_S1511, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1512;
    (&_S1512)->primal_0 = 0.0f;
    (&_S1512)->differential_0 = 0.0f;
    DiffPair_float_0 _S1513;
    (&_S1513)->primal_0 = dpparams_3->primal_0[int(5)];
    (&_S1513)->differential_0 = 0.0f;
    _d_max_0(&_S1512, &_S1513, (*_s_dOut_3)[int(4)]);
    DiffPair_float_0 _S1514;
    (&_S1514)->primal_0 = 0.0f;
    (&_S1514)->differential_0 = 0.0f;
    DiffPair_float_0 _S1515;
    (&_S1515)->primal_0 = dpparams_3->primal_0[int(14)];
    (&_S1515)->differential_0 = 0.0f;
    _d_max_0(&_S1514, &_S1515, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1516;
    (&_S1516)->primal_0 = 0.0f;
    (&_S1516)->differential_0 = 0.0f;
    DiffPair_float_0 _S1517;
    (&_S1517)->primal_0 = dpparams_3->primal_0[int(9)];
    (&_S1517)->differential_0 = 0.0f;
    _d_max_0(&_S1516, &_S1517, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1518;
    (&_S1518)->primal_0 = 0.0f;
    (&_S1518)->differential_0 = 0.0f;
    DiffPair_float_0 _S1519;
    (&_S1519)->primal_0 = dpparams_3->primal_0[int(4)];
    (&_S1519)->differential_0 = 0.0f;
    _d_max_0(&_S1518, &_S1519, (*_s_dOut_3)[int(3)]);
    DiffPair_float_0 _S1520;
    (&_S1520)->primal_0 = 0.0f;
    (&_S1520)->differential_0 = 0.0f;
    DiffPair_float_0 _S1521;
    (&_S1521)->primal_0 = dpparams_3->primal_0[int(13)];
    (&_S1521)->differential_0 = 0.0f;
    _d_max_0(&_S1520, &_S1521, (*_s_dOut_3)[int(2)]);
    DiffPair_float_0 _S1522;
    (&_S1522)->primal_0 = 0.0f;
    (&_S1522)->differential_0 = 0.0f;
    DiffPair_float_0 _S1523;
    (&_S1523)->primal_0 = dpparams_3->primal_0[int(8)];
    (&_S1523)->differential_0 = 0.0f;
    _d_max_0(&_S1522, &_S1523, (*_s_dOut_3)[int(2)]);
    DiffPair_float_0 _S1524;
    (&_S1524)->primal_0 = 0.0f;
    (&_S1524)->differential_0 = 0.0f;
    DiffPair_float_0 _S1525;
    (&_S1525)->primal_0 = dpparams_3->primal_0[int(3)];
    (&_S1525)->differential_0 = 0.0f;
    _d_max_0(&_S1524, &_S1525, (*_s_dOut_3)[int(2)]);
    float _S1526 = dpparams_3->primal_0[int(12)] * (*_s_dOut_3)[int(1)];
    float _S1527 = dpparams_3->primal_0[int(11)] * (*_s_dOut_3)[int(1)];
    float _S1528 = dpparams_3->primal_0[int(7)] * (*_s_dOut_3)[int(1)];
    float _S1529 = dpparams_3->primal_0[int(6)] * (*_s_dOut_3)[int(1)];
    float _S1530 = dpparams_3->primal_0[int(2)] * (*_s_dOut_3)[int(1)];
    float _S1531 = dpparams_3->primal_0[int(1)] * (*_s_dOut_3)[int(1)];
    PPISPParamsRQS_0 _S1532 = PPISPParamsRQS_x24_syn_dzero_0();
    (&_S1532)->color_params_0 = _S1467;
    (&_S1532)->exposure_0 = (*_s_dOut_3)[int(0)];
    _S1398 = _S1532;
    (&(&_S1398)->crf_params_0[int(2)])->gc_0 = 0.0f;
    float _S1533 = _S1416 + _S1421 + _S1532.crf_params_0[int(2)].gc_0;
    (&(&_S1398)->crf_params_0[int(2)])->y0_0 = 0.0f;
    float _S1534 = _S1424 + _S1429 + _S1532.crf_params_0[int(2)].y0_0;
    (&(&_S1398)->crf_params_0[int(2)])->x0_0 = 0.0f;
    float _S1535 = _S1432 + _S1437 + _S1532.crf_params_0[int(2)].x0_0;
    (&(&_S1398)->crf_params_0[int(2)])->g1_0 = 0.0f;
    float _S1536 = _S1440 + _S1445 + _S1532.crf_params_0[int(2)].g1_0;
    (&(&_S1398)->crf_params_0[int(2)])->g0_0 = 0.0f;
    float _S1537 = _S1448 + _S1453 + _S1532.crf_params_0[int(2)].g0_0;
    (&(&_S1398)->crf_params_0[int(1)])->gc_0 = 0.0f;
    float _S1538 = _S1418 + _S1421 + _S1532.crf_params_0[int(1)].gc_0;
    (&(&_S1398)->crf_params_0[int(1)])->y0_0 = 0.0f;
    float _S1539 = _S1426 + _S1429 + _S1532.crf_params_0[int(1)].y0_0;
    (&(&_S1398)->crf_params_0[int(1)])->x0_0 = 0.0f;
    float _S1540 = _S1434 + _S1437 + _S1532.crf_params_0[int(1)].x0_0;
    (&(&_S1398)->crf_params_0[int(1)])->g1_0 = 0.0f;
    float _S1541 = _S1442 + _S1445 + _S1532.crf_params_0[int(1)].g1_0;
    (&(&_S1398)->crf_params_0[int(1)])->g0_0 = 0.0f;
    float _S1542 = _S1450 + _S1453 + _S1532.crf_params_0[int(1)].g0_0;
    (&(&_S1398)->crf_params_0[int(0)])->gc_0 = 0.0f;
    float _S1543 = _S1420 + _S1421 + _S1532.crf_params_0[int(0)].gc_0;
    (&(&_S1398)->crf_params_0[int(0)])->y0_0 = 0.0f;
    float _S1544 = _S1428 + _S1429 + _S1532.crf_params_0[int(0)].y0_0;
    (&(&_S1398)->crf_params_0[int(0)])->x0_0 = 0.0f;
    float _S1545 = _S1436 + _S1437 + _S1532.crf_params_0[int(0)].x0_0;
    (&(&_S1398)->crf_params_0[int(0)])->g1_0 = 0.0f;
    float _S1546 = _S1444 + _S1445 + _S1532.crf_params_0[int(0)].g1_0;
    (&(&_S1398)->crf_params_0[int(0)])->g0_0 = 0.0f;
    float _S1547 = _S1452 + _S1453 + _S1532.crf_params_0[int(0)].g0_0;
    *&((&(&(&_S1398)->color_params_0)->n_0)->y) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->n_0)->x) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->g_0)->y) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->g_0)->x) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->r_0)->y) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->r_0)->x) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->b_0)->y) = 0.0f;
    *&((&(&(&_S1398)->color_params_0)->b_0)->x) = 0.0f;
    (&(&_S1398)->vignette_params_0[int(2)])->alpha2_0 = 0.0f;
    float _S1548 = _S1470 + _S1475 + _S1509.differential_0 + _S1532.vignette_params_0[int(2)].alpha2_0;
    (&(&_S1398)->vignette_params_0[int(2)])->alpha1_0 = 0.0f;
    float _S1549 = _S1478 + _S1483 + _S1515.differential_0 + _S1532.vignette_params_0[int(2)].alpha1_0;
    (&(&_S1398)->vignette_params_0[int(2)])->alpha0_0 = 0.0f;
    float _S1550 = _S1486 + _S1491 + _S1521.differential_0 + _S1532.vignette_params_0[int(2)].alpha0_0;
    (&(&_S1398)->vignette_params_0[int(2)])->cy_0 = 0.0f;
    float _S1551 = _S1494 + _S1499 + _S1526 + _S1526 + _S1532.vignette_params_0[int(2)].cy_0;
    (&(&_S1398)->vignette_params_0[int(2)])->cx_0 = 0.0f;
    float _S1552 = _S1502 + _S1507 + _S1527 + _S1527 + _S1532.vignette_params_0[int(2)].cx_0;
    (&(&_S1398)->vignette_params_0[int(1)])->alpha2_0 = 0.0f;
    float _S1553 = _S1472 + _S1475 + _S1511.differential_0 + _S1532.vignette_params_0[int(1)].alpha2_0;
    (&(&_S1398)->vignette_params_0[int(1)])->alpha1_0 = 0.0f;
    float _S1554 = _S1480 + _S1483 + _S1517.differential_0 + _S1532.vignette_params_0[int(1)].alpha1_0;
    (&(&_S1398)->vignette_params_0[int(1)])->alpha0_0 = 0.0f;
    float _S1555 = _S1488 + _S1491 + _S1523.differential_0 + _S1532.vignette_params_0[int(1)].alpha0_0;
    (&(&_S1398)->vignette_params_0[int(1)])->cy_0 = 0.0f;
    float _S1556 = _S1496 + _S1499 + _S1528 + _S1528 + _S1532.vignette_params_0[int(1)].cy_0;
    (&(&_S1398)->vignette_params_0[int(1)])->cx_0 = 0.0f;
    float _S1557 = _S1504 + _S1507 + _S1529 + _S1529 + _S1532.vignette_params_0[int(1)].cx_0;
    (&(&_S1398)->vignette_params_0[int(0)])->alpha2_0 = 0.0f;
    float _S1558 = _S1474 + _S1475 + _S1513.differential_0 + _S1532.vignette_params_0[int(0)].alpha2_0;
    (&(&_S1398)->vignette_params_0[int(0)])->alpha1_0 = 0.0f;
    float _S1559 = _S1482 + _S1483 + _S1519.differential_0 + _S1532.vignette_params_0[int(0)].alpha1_0;
    (&(&_S1398)->vignette_params_0[int(0)])->alpha0_0 = 0.0f;
    float _S1560 = _S1490 + _S1491 + _S1525.differential_0 + _S1532.vignette_params_0[int(0)].alpha0_0;
    (&(&_S1398)->vignette_params_0[int(0)])->cy_0 = 0.0f;
    float _S1561 = _S1498 + _S1499 + _S1530 + _S1530 + _S1532.vignette_params_0[int(0)].cy_0;
    (&(&_S1398)->vignette_params_0[int(0)])->cx_0 = 0.0f;
    float _S1562 = _S1506 + _S1507 + _S1531 + _S1531 + _S1532.vignette_params_0[int(0)].cx_0;
    FixedArray<float, 39>  _S1563;
    _S1563[int(0)] = 0.0f;
    _S1563[int(1)] = 0.0f;
    _S1563[int(2)] = 0.0f;
    _S1563[int(3)] = 0.0f;
    _S1563[int(4)] = 0.0f;
    _S1563[int(5)] = 0.0f;
    _S1563[int(6)] = 0.0f;
    _S1563[int(7)] = 0.0f;
    _S1563[int(8)] = 0.0f;
    _S1563[int(9)] = 0.0f;
    _S1563[int(10)] = 0.0f;
    _S1563[int(11)] = 0.0f;
    _S1563[int(12)] = 0.0f;
    _S1563[int(13)] = 0.0f;
    _S1563[int(14)] = 0.0f;
    _S1563[int(15)] = 0.0f;
    _S1563[int(16)] = 0.0f;
    _S1563[int(17)] = 0.0f;
    _S1563[int(18)] = 0.0f;
    _S1563[int(19)] = 0.0f;
    _S1563[int(20)] = 0.0f;
    _S1563[int(21)] = 0.0f;
    _S1563[int(22)] = 0.0f;
    _S1563[int(23)] = 0.0f;
    _S1563[int(24)] = 0.0f;
    _S1563[int(25)] = 0.0f;
    _S1563[int(26)] = 0.0f;
    _S1563[int(27)] = 0.0f;
    _S1563[int(28)] = 0.0f;
    _S1563[int(29)] = 0.0f;
    _S1563[int(30)] = 0.0f;
    _S1563[int(31)] = 0.0f;
    _S1563[int(32)] = 0.0f;
    _S1563[int(33)] = 0.0f;
    _S1563[int(34)] = 0.0f;
    _S1563[int(35)] = 0.0f;
    _S1563[int(36)] = 0.0f;
    _S1563[int(37)] = 0.0f;
    _S1563[int(38)] = 0.0f;
    _S1563[int(9)] = _S1554;
    _S1563[int(18)] = _S1532.color_params_0.r_0.x;
    _S1563[int(17)] = _S1532.color_params_0.b_0.y;
    _S1563[int(16)] = _S1532.color_params_0.b_0.x;
    _S1563[int(15)] = _S1548;
    _S1563[int(14)] = _S1549;
    _S1563[int(13)] = _S1550;
    _S1563[int(12)] = _S1551;
    _S1563[int(11)] = _S1552;
    _S1563[int(10)] = _S1553;
    _S1563[int(19)] = _S1532.color_params_0.r_0.y;
    _S1563[int(8)] = _S1555;
    _S1563[int(7)] = _S1556;
    _S1563[int(6)] = _S1557;
    _S1563[int(5)] = _S1558;
    _S1563[int(4)] = _S1559;
    _S1563[int(3)] = _S1560;
    _S1563[int(2)] = _S1561;
    _S1563[int(1)] = _S1562;
    _S1563[int(0)] = _S1398.exposure_0;
    _S1563[int(28)] = _S1543;
    _S1563[int(37)] = _S1534;
    _S1563[int(36)] = _S1535;
    _S1563[int(35)] = _S1536;
    _S1563[int(34)] = _S1537;
    _S1563[int(33)] = _S1538;
    _S1563[int(32)] = _S1539;
    _S1563[int(31)] = _S1540;
    _S1563[int(30)] = _S1541;
    _S1563[int(29)] = _S1542;
    _S1563[int(38)] = _S1533;
    _S1563[int(27)] = _S1544;
    _S1563[int(26)] = _S1545;
    _S1563[int(25)] = _S1546;
    _S1563[int(24)] = _S1547;
    _S1563[int(23)] = _S1532.color_params_0.n_0.y;
    _S1563[int(22)] = _S1532.color_params_0.n_0.x;
    _S1563[int(21)] = _S1532.color_params_0.g_0.y;
    _S1563[int(20)] = _S1532.color_params_0.g_0.x;
    dpparams_3->primal_0 = dpparams_3->primal_0;
    dpparams_3->differential_0 = _S1563;
    return;
}

inline __device__ void s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C39x3E_0 * _S1564, FixedArray<float, 23>  * _S1565)
{
    s_bwd_prop_compute_raw_ppisp_rqs_regularization_loss_0(_S1564, _S1565);
    return;
}

inline __device__ void compute_raw_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 39>  params_7, FixedArray<float, 23>  grad_out_3, FixedArray<float, 39>  * _S1566)
{
    FixedArray<float, 39>  _S1567 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C39x3E_0 dp_params_3;
    (&dp_params_3)->primal_0 = params_7;
    (&dp_params_3)->differential_0 = _S1567;
    FixedArray<float, 23>  _S1568 = grad_out_3;
    s_bwd_compute_raw_ppisp_rqs_regularization_loss_0(&dp_params_3, &_S1568);
    *_S1566 = (&dp_params_3)->differential_0;
    return;
}

inline __device__ void compute_ppisp_regularization_loss(FixedArray<float, 22>  raw_losses_0, int num_cameras_0, FixedArray<float, 6>  loss_weights_0, FixedArray<float, 6>  * _S1569)
{
    float _S1570;
    FixedArray<float, 6>  losses_2;
    float _S1571 = float(num_cameras_0);
    float _S1572 = raw_losses_0[int(0)] / _S1571;
    for(;;)
    {
        float _S1573 = (F32_abs((_S1572)));
        if(_S1573 < 0.10000000149011612f)
        {
            _S1570 = 0.5f * _S1572 * _S1572 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1570 = _S1573 - 0.05000000074505806f;
            break;
        }
    }
    losses_2[int(0)] = _S1570;
    losses_2[int(1)] = raw_losses_0[int(1)] / (3.0f * _S1571);
    losses_2[int(2)] = (raw_losses_0[int(2)] + raw_losses_0[int(3)] + raw_losses_0[int(4)]) / (9.0f * _S1571);
    losses_2[int(3)] = (raw_losses_0[int(5)] + raw_losses_0[int(6)] + raw_losses_0[int(7)] + raw_losses_0[int(8)] + raw_losses_0[int(9)]) / (5.0f * _S1571);
    float _S1574 = raw_losses_0[int(10)] / _S1571;
    for(;;)
    {
        float _S1575 = (F32_abs((_S1574)));
        if(_S1575 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1574 * _S1574 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1575 - 0.00249999994412065f;
            break;
        }
    }
    float _S1576;
    float _S1577 = raw_losses_0[int(11)] / _S1571;
    for(;;)
    {
        float _S1578 = (F32_abs((_S1577)));
        if(_S1578 < 0.00499999988824129f)
        {
            _S1576 = 0.5f * _S1577 * _S1577 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1576 = _S1578 - 0.00249999994412065f;
            break;
        }
    }
    float _S1579 = _S1570 + _S1576;
    float _S1580 = raw_losses_0[int(12)] / _S1571;
    for(;;)
    {
        float _S1581 = (F32_abs((_S1580)));
        if(_S1581 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1580 * _S1580 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1581 - 0.00249999994412065f;
            break;
        }
    }
    float _S1582 = _S1579 + _S1570;
    float _S1583 = raw_losses_0[int(13)] / _S1571;
    for(;;)
    {
        float _S1584 = (F32_abs((_S1583)));
        if(_S1584 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1583 * _S1583 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1584 - 0.00249999994412065f;
            break;
        }
    }
    float _S1585 = _S1582 + _S1570;
    float _S1586 = raw_losses_0[int(14)] / _S1571;
    for(;;)
    {
        float _S1587 = (F32_abs((_S1586)));
        if(_S1587 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1586 * _S1586 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1587 - 0.00249999994412065f;
            break;
        }
    }
    float _S1588 = _S1585 + _S1570;
    float _S1589 = raw_losses_0[int(15)] / _S1571;
    for(;;)
    {
        float _S1590 = (F32_abs((_S1589)));
        if(_S1590 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1589 * _S1589 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1590 - 0.00249999994412065f;
            break;
        }
    }
    float _S1591 = _S1588 + _S1570;
    float _S1592 = raw_losses_0[int(16)] / _S1571;
    for(;;)
    {
        float _S1593 = (F32_abs((_S1592)));
        if(_S1593 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1592 * _S1592 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1593 - 0.00249999994412065f;
            break;
        }
    }
    float _S1594 = _S1591 + _S1570;
    float _S1595 = raw_losses_0[int(17)] / _S1571;
    for(;;)
    {
        float _S1596 = (F32_abs((_S1595)));
        if(_S1596 < 0.00499999988824129f)
        {
            _S1570 = 0.5f * _S1595 * _S1595 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1570 = _S1596 - 0.00249999994412065f;
            break;
        }
    }
    float _S1597 = (_S1594 + _S1570) / 8.0f;
    float _S1598 = (raw_losses_0[int(18)] + raw_losses_0[int(19)] + raw_losses_0[int(20)] + raw_losses_0[int(21)]) / (4.0f * _S1571);
    losses_2[int(0)] = losses_2[int(0)] * loss_weights_0[int(0)];
    losses_2[int(1)] = losses_2[int(1)] * loss_weights_0[int(1)];
    losses_2[int(2)] = losses_2[int(2)] * loss_weights_0[int(2)];
    losses_2[int(3)] = losses_2[int(3)] * loss_weights_0[int(3)];
    losses_2[int(4)] = _S1597 * loss_weights_0[int(4)];
    losses_2[int(5)] = _S1598 * loss_weights_0[int(5)];
    *_S1569 = losses_2;
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss(FixedArray<float, 23>  raw_losses_1, int num_cameras_1, FixedArray<float, 6>  loss_weights_1, FixedArray<float, 6>  * _S1599)
{
    float _S1600;
    FixedArray<float, 6>  losses_3;
    float _S1601 = float(num_cameras_1);
    float _S1602 = raw_losses_1[int(0)] / _S1601;
    for(;;)
    {
        float _S1603 = (F32_abs((_S1602)));
        if(_S1603 < 0.10000000149011612f)
        {
            _S1600 = 0.5f * _S1602 * _S1602 / 0.10000000149011612f;
            break;
        }
        else
        {
            _S1600 = _S1603 - 0.05000000074505806f;
            break;
        }
    }
    losses_3[int(0)] = _S1600;
    losses_3[int(1)] = raw_losses_1[int(1)] / (3.0f * _S1601);
    losses_3[int(2)] = (raw_losses_1[int(2)] + raw_losses_1[int(3)] + raw_losses_1[int(4)]) / (9.0f * _S1601);
    float _S1604 = 5.0f * _S1601;
    losses_3[int(3)] = (raw_losses_1[int(5)] + raw_losses_1[int(6)] + raw_losses_1[int(7)] + raw_losses_1[int(8)] + raw_losses_1[int(9)]) / _S1604;
    float _S1605 = raw_losses_1[int(10)] / _S1601;
    for(;;)
    {
        float _S1606 = (F32_abs((_S1605)));
        if(_S1606 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1605 * _S1605 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1606 - 0.00249999994412065f;
            break;
        }
    }
    float _S1607;
    float _S1608 = raw_losses_1[int(11)] / _S1601;
    for(;;)
    {
        float _S1609 = (F32_abs((_S1608)));
        if(_S1609 < 0.00499999988824129f)
        {
            _S1607 = 0.5f * _S1608 * _S1608 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1607 = _S1609 - 0.00249999994412065f;
            break;
        }
    }
    float _S1610 = _S1600 + _S1607;
    float _S1611 = raw_losses_1[int(12)] / _S1601;
    for(;;)
    {
        float _S1612 = (F32_abs((_S1611)));
        if(_S1612 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1611 * _S1611 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1612 - 0.00249999994412065f;
            break;
        }
    }
    float _S1613 = _S1610 + _S1600;
    float _S1614 = raw_losses_1[int(13)] / _S1601;
    for(;;)
    {
        float _S1615 = (F32_abs((_S1614)));
        if(_S1615 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1614 * _S1614 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1615 - 0.00249999994412065f;
            break;
        }
    }
    float _S1616 = _S1613 + _S1600;
    float _S1617 = raw_losses_1[int(14)] / _S1601;
    for(;;)
    {
        float _S1618 = (F32_abs((_S1617)));
        if(_S1618 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1617 * _S1617 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1618 - 0.00249999994412065f;
            break;
        }
    }
    float _S1619 = _S1616 + _S1600;
    float _S1620 = raw_losses_1[int(15)] / _S1601;
    for(;;)
    {
        float _S1621 = (F32_abs((_S1620)));
        if(_S1621 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1620 * _S1620 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1621 - 0.00249999994412065f;
            break;
        }
    }
    float _S1622 = _S1619 + _S1600;
    float _S1623 = raw_losses_1[int(16)] / _S1601;
    for(;;)
    {
        float _S1624 = (F32_abs((_S1623)));
        if(_S1624 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1623 * _S1623 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1624 - 0.00249999994412065f;
            break;
        }
    }
    float _S1625 = _S1622 + _S1600;
    float _S1626 = raw_losses_1[int(17)] / _S1601;
    for(;;)
    {
        float _S1627 = (F32_abs((_S1626)));
        if(_S1627 < 0.00499999988824129f)
        {
            _S1600 = 0.5f * _S1626 * _S1626 / 0.00499999988824129f;
            break;
        }
        else
        {
            _S1600 = _S1627 - 0.00249999994412065f;
            break;
        }
    }
    float _S1628 = (_S1625 + _S1600) / 8.0f;
    float _S1629 = (raw_losses_1[int(18)] + raw_losses_1[int(19)] + raw_losses_1[int(20)] + raw_losses_1[int(21)] + raw_losses_1[int(22)]) / _S1604;
    losses_3[int(0)] = losses_3[int(0)] * loss_weights_1[int(0)];
    losses_3[int(1)] = losses_3[int(1)] * loss_weights_1[int(1)];
    losses_3[int(2)] = losses_3[int(2)] * loss_weights_1[int(2)];
    losses_3[int(3)] = losses_3[int(3)] * loss_weights_1[int(3)];
    losses_3[int(4)] = _S1628 * loss_weights_1[int(4)];
    losses_3[int(5)] = _S1629 * loss_weights_1[int(5)];
    *_S1599 = losses_3;
    return;
}

struct DiffPair_arrayx3Cfloatx2C22x3E_0
{
    FixedArray<float, 22>  primal_0;
    FixedArray<float, 22>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * dpraw_losses_0, int num_cameras_2, FixedArray<float, 6>  * loss_weights_2, FixedArray<float, 6>  * _s_dOut_4)
{
    FixedArray<float, 22>  _S1630 = dpraw_losses_0->primal_0;
    float _S1631 = float(num_cameras_2);
    float _S1632 = dpraw_losses_0->primal_0[int(0)] / _S1631;
    bool _S1633 = (s_primal_ctx_abs_0(_S1632)) < 0.10000000149011612f;
    float _S1634;
    if(_S1633)
    {
        _S1634 = 0.5f * _S1632;
    }
    else
    {
        _S1634 = 0.0f;
    }
    float _S1635 = 3.0f * _S1631;
    float _S1636 = 9.0f * _S1631;
    float _S1637 = 5.0f * _S1631;
    float _S1638 = _S1630[int(10)] / _S1631;
    bool _S1639 = (s_primal_ctx_abs_0(_S1638)) < 0.00499999988824129f;
    float _S1640;
    if(_S1639)
    {
        _S1640 = 0.5f * _S1638;
    }
    else
    {
        _S1640 = 0.0f;
    }
    float _S1641 = _S1630[int(11)] / _S1631;
    bool _S1642 = (s_primal_ctx_abs_0(_S1641)) < 0.00499999988824129f;
    float _S1643;
    if(_S1642)
    {
        _S1643 = 0.5f * _S1641;
    }
    else
    {
        _S1643 = 0.0f;
    }
    float _S1644 = _S1630[int(12)] / _S1631;
    bool _S1645 = (s_primal_ctx_abs_0(_S1644)) < 0.00499999988824129f;
    float _S1646;
    if(_S1645)
    {
        _S1646 = 0.5f * _S1644;
    }
    else
    {
        _S1646 = 0.0f;
    }
    float _S1647 = _S1630[int(13)] / _S1631;
    bool _S1648 = (s_primal_ctx_abs_0(_S1647)) < 0.00499999988824129f;
    float _S1649;
    if(_S1648)
    {
        _S1649 = 0.5f * _S1647;
    }
    else
    {
        _S1649 = 0.0f;
    }
    float _S1650 = _S1630[int(14)] / _S1631;
    bool _S1651 = (s_primal_ctx_abs_0(_S1650)) < 0.00499999988824129f;
    float _S1652;
    if(_S1651)
    {
        _S1652 = 0.5f * _S1650;
    }
    else
    {
        _S1652 = 0.0f;
    }
    float _S1653 = _S1630[int(15)] / _S1631;
    bool _S1654 = (s_primal_ctx_abs_0(_S1653)) < 0.00499999988824129f;
    float _S1655;
    if(_S1654)
    {
        _S1655 = 0.5f * _S1653;
    }
    else
    {
        _S1655 = 0.0f;
    }
    float _S1656 = _S1630[int(16)] / _S1631;
    bool _S1657 = (s_primal_ctx_abs_0(_S1656)) < 0.00499999988824129f;
    float _S1658;
    if(_S1657)
    {
        _S1658 = 0.5f * _S1656;
    }
    else
    {
        _S1658 = 0.0f;
    }
    float _S1659 = _S1630[int(17)] / _S1631;
    bool _S1660 = (s_primal_ctx_abs_0(_S1659)) < 0.00499999988824129f;
    float _S1661;
    if(_S1660)
    {
        _S1661 = 0.5f * _S1659;
    }
    else
    {
        _S1661 = 0.0f;
    }
    float _S1662 = (*loss_weights_2)[int(3)] * (*_s_dOut_4)[int(3)];
    float _S1663 = (*loss_weights_2)[int(2)] * (*_s_dOut_4)[int(2)];
    float _S1664 = (*loss_weights_2)[int(1)] * (*_s_dOut_4)[int(1)];
    float _S1665 = (*loss_weights_2)[int(0)] * (*_s_dOut_4)[int(0)];
    float _S1666 = (*loss_weights_2)[int(5)] * (*_s_dOut_4)[int(5)] / (4.0f * _S1631);
    float _S1667 = 0.125f * ((*loss_weights_2)[int(4)] * (*_s_dOut_4)[int(4)]);
    FixedArray<float, 22>  _S1668;
    _S1668[int(0)] = 0.0f;
    _S1668[int(1)] = 0.0f;
    _S1668[int(2)] = 0.0f;
    _S1668[int(3)] = 0.0f;
    _S1668[int(4)] = 0.0f;
    _S1668[int(5)] = 0.0f;
    _S1668[int(6)] = 0.0f;
    _S1668[int(7)] = 0.0f;
    _S1668[int(8)] = 0.0f;
    _S1668[int(9)] = 0.0f;
    _S1668[int(10)] = 0.0f;
    _S1668[int(11)] = 0.0f;
    _S1668[int(12)] = 0.0f;
    _S1668[int(13)] = 0.0f;
    _S1668[int(14)] = 0.0f;
    _S1668[int(15)] = 0.0f;
    _S1668[int(16)] = 0.0f;
    _S1668[int(17)] = 0.0f;
    _S1668[int(18)] = 0.0f;
    _S1668[int(19)] = 0.0f;
    _S1668[int(20)] = 0.0f;
    _S1668[int(21)] = 0.0f;
    _S1668[int(21)] = _S1666;
    _S1668[int(20)] = _S1666;
    _S1668[int(19)] = _S1666;
    _S1668[int(18)] = _S1666;
    float _S1669 = _S1668[int(0)];
    float _S1670 = _S1668[int(1)];
    float _S1671 = _S1668[int(2)];
    float _S1672 = _S1668[int(3)];
    float _S1673 = _S1668[int(4)];
    float _S1674 = _S1668[int(5)];
    float _S1675 = _S1668[int(6)];
    float _S1676 = _S1668[int(7)];
    float _S1677 = _S1668[int(8)];
    float _S1678 = _S1668[int(9)];
    float _S1679 = _S1668[int(10)];
    float _S1680 = _S1668[int(11)];
    float _S1681 = _S1668[int(12)];
    float _S1682 = _S1668[int(13)];
    float _S1683 = _S1668[int(14)];
    float _S1684 = _S1668[int(15)];
    float _S1685 = _S1668[int(16)];
    float _S1686 = _S1668[int(17)];
    float _S1687 = _S1668[int(18)];
    float _S1688 = _S1668[int(19)];
    float _S1689 = _S1668[int(20)];
    float _S1690 = _S1668[int(21)];
    float _S1691;
    if(_S1660)
    {
        float _S1692 = 200.0f * _S1667;
        float _S1693 = _S1661 * _S1692 + 0.5f * (_S1659 * _S1692);
        _S1661 = 0.0f;
        _S1691 = _S1693;
    }
    else
    {
        _S1661 = _S1667;
        _S1691 = 0.0f;
    }
    DiffPair_float_0 _S1694;
    (&_S1694)->primal_0 = _S1659;
    (&_S1694)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1694, _S1661);
    float _S1695 = (_S1694.differential_0 + _S1691) / _S1631;
    FixedArray<float, 22>  _S1696;
    _S1696[int(0)] = 0.0f;
    _S1696[int(1)] = 0.0f;
    _S1696[int(2)] = 0.0f;
    _S1696[int(3)] = 0.0f;
    _S1696[int(4)] = 0.0f;
    _S1696[int(5)] = 0.0f;
    _S1696[int(6)] = 0.0f;
    _S1696[int(7)] = 0.0f;
    _S1696[int(8)] = 0.0f;
    _S1696[int(9)] = 0.0f;
    _S1696[int(10)] = 0.0f;
    _S1696[int(11)] = 0.0f;
    _S1696[int(12)] = 0.0f;
    _S1696[int(13)] = 0.0f;
    _S1696[int(14)] = 0.0f;
    _S1696[int(15)] = 0.0f;
    _S1696[int(16)] = 0.0f;
    _S1696[int(17)] = 0.0f;
    _S1696[int(18)] = 0.0f;
    _S1696[int(19)] = 0.0f;
    _S1696[int(20)] = 0.0f;
    _S1696[int(21)] = 0.0f;
    _S1696[int(17)] = _S1695;
    float _S1697 = _S1669 + _S1696[int(0)];
    float _S1698 = _S1670 + _S1696[int(1)];
    float _S1699 = _S1671 + _S1696[int(2)];
    float _S1700 = _S1672 + _S1696[int(3)];
    float _S1701 = _S1673 + _S1696[int(4)];
    float _S1702 = _S1674 + _S1696[int(5)];
    float _S1703 = _S1675 + _S1696[int(6)];
    float _S1704 = _S1676 + _S1696[int(7)];
    float _S1705 = _S1677 + _S1696[int(8)];
    float _S1706 = _S1678 + _S1696[int(9)];
    float _S1707 = _S1679 + _S1696[int(10)];
    float _S1708 = _S1680 + _S1696[int(11)];
    float _S1709 = _S1681 + _S1696[int(12)];
    float _S1710 = _S1682 + _S1696[int(13)];
    float _S1711 = _S1683 + _S1696[int(14)];
    float _S1712 = _S1684 + _S1696[int(15)];
    float _S1713 = _S1685 + _S1696[int(16)];
    float _S1714 = _S1686 + _S1696[int(17)];
    float _S1715 = _S1687 + _S1696[int(18)];
    float _S1716 = _S1688 + _S1696[int(19)];
    float _S1717 = _S1689 + _S1696[int(20)];
    float _S1718 = _S1690 + _S1696[int(21)];
    if(_S1657)
    {
        float _S1719 = 200.0f * _S1667;
        float _S1720 = _S1658 * _S1719 + 0.5f * (_S1656 * _S1719);
        _S1658 = 0.0f;
        _S1661 = _S1720;
    }
    else
    {
        _S1658 = _S1667;
        _S1661 = 0.0f;
    }
    DiffPair_float_0 _S1721;
    (&_S1721)->primal_0 = _S1656;
    (&_S1721)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1721, _S1658);
    float _S1722 = (_S1721.differential_0 + _S1661) / _S1631;
    FixedArray<float, 22>  _S1723;
    _S1723[int(0)] = 0.0f;
    _S1723[int(1)] = 0.0f;
    _S1723[int(2)] = 0.0f;
    _S1723[int(3)] = 0.0f;
    _S1723[int(4)] = 0.0f;
    _S1723[int(5)] = 0.0f;
    _S1723[int(6)] = 0.0f;
    _S1723[int(7)] = 0.0f;
    _S1723[int(8)] = 0.0f;
    _S1723[int(9)] = 0.0f;
    _S1723[int(10)] = 0.0f;
    _S1723[int(11)] = 0.0f;
    _S1723[int(12)] = 0.0f;
    _S1723[int(13)] = 0.0f;
    _S1723[int(14)] = 0.0f;
    _S1723[int(15)] = 0.0f;
    _S1723[int(16)] = 0.0f;
    _S1723[int(17)] = 0.0f;
    _S1723[int(18)] = 0.0f;
    _S1723[int(19)] = 0.0f;
    _S1723[int(20)] = 0.0f;
    _S1723[int(21)] = 0.0f;
    _S1723[int(16)] = _S1722;
    float _S1724 = _S1697 + _S1723[int(0)];
    float _S1725 = _S1698 + _S1723[int(1)];
    float _S1726 = _S1699 + _S1723[int(2)];
    float _S1727 = _S1700 + _S1723[int(3)];
    float _S1728 = _S1701 + _S1723[int(4)];
    float _S1729 = _S1702 + _S1723[int(5)];
    float _S1730 = _S1703 + _S1723[int(6)];
    float _S1731 = _S1704 + _S1723[int(7)];
    float _S1732 = _S1705 + _S1723[int(8)];
    float _S1733 = _S1706 + _S1723[int(9)];
    float _S1734 = _S1707 + _S1723[int(10)];
    float _S1735 = _S1708 + _S1723[int(11)];
    float _S1736 = _S1709 + _S1723[int(12)];
    float _S1737 = _S1710 + _S1723[int(13)];
    float _S1738 = _S1711 + _S1723[int(14)];
    float _S1739 = _S1712 + _S1723[int(15)];
    float _S1740 = _S1713 + _S1723[int(16)];
    float _S1741 = _S1714 + _S1723[int(17)];
    float _S1742 = _S1715 + _S1723[int(18)];
    float _S1743 = _S1716 + _S1723[int(19)];
    float _S1744 = _S1717 + _S1723[int(20)];
    float _S1745 = _S1718 + _S1723[int(21)];
    if(_S1654)
    {
        float _S1746 = 200.0f * _S1667;
        float _S1747 = _S1655 * _S1746 + 0.5f * (_S1653 * _S1746);
        _S1655 = 0.0f;
        _S1658 = _S1747;
    }
    else
    {
        _S1655 = _S1667;
        _S1658 = 0.0f;
    }
    DiffPair_float_0 _S1748;
    (&_S1748)->primal_0 = _S1653;
    (&_S1748)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1748, _S1655);
    float _S1749 = (_S1748.differential_0 + _S1658) / _S1631;
    FixedArray<float, 22>  _S1750;
    _S1750[int(0)] = 0.0f;
    _S1750[int(1)] = 0.0f;
    _S1750[int(2)] = 0.0f;
    _S1750[int(3)] = 0.0f;
    _S1750[int(4)] = 0.0f;
    _S1750[int(5)] = 0.0f;
    _S1750[int(6)] = 0.0f;
    _S1750[int(7)] = 0.0f;
    _S1750[int(8)] = 0.0f;
    _S1750[int(9)] = 0.0f;
    _S1750[int(10)] = 0.0f;
    _S1750[int(11)] = 0.0f;
    _S1750[int(12)] = 0.0f;
    _S1750[int(13)] = 0.0f;
    _S1750[int(14)] = 0.0f;
    _S1750[int(15)] = 0.0f;
    _S1750[int(16)] = 0.0f;
    _S1750[int(17)] = 0.0f;
    _S1750[int(18)] = 0.0f;
    _S1750[int(19)] = 0.0f;
    _S1750[int(20)] = 0.0f;
    _S1750[int(21)] = 0.0f;
    _S1750[int(15)] = _S1749;
    float _S1751 = _S1724 + _S1750[int(0)];
    float _S1752 = _S1725 + _S1750[int(1)];
    float _S1753 = _S1726 + _S1750[int(2)];
    float _S1754 = _S1727 + _S1750[int(3)];
    float _S1755 = _S1728 + _S1750[int(4)];
    float _S1756 = _S1729 + _S1750[int(5)];
    float _S1757 = _S1730 + _S1750[int(6)];
    float _S1758 = _S1731 + _S1750[int(7)];
    float _S1759 = _S1732 + _S1750[int(8)];
    float _S1760 = _S1733 + _S1750[int(9)];
    float _S1761 = _S1734 + _S1750[int(10)];
    float _S1762 = _S1735 + _S1750[int(11)];
    float _S1763 = _S1736 + _S1750[int(12)];
    float _S1764 = _S1737 + _S1750[int(13)];
    float _S1765 = _S1738 + _S1750[int(14)];
    float _S1766 = _S1739 + _S1750[int(15)];
    float _S1767 = _S1740 + _S1750[int(16)];
    float _S1768 = _S1741 + _S1750[int(17)];
    float _S1769 = _S1742 + _S1750[int(18)];
    float _S1770 = _S1743 + _S1750[int(19)];
    float _S1771 = _S1744 + _S1750[int(20)];
    float _S1772 = _S1745 + _S1750[int(21)];
    if(_S1651)
    {
        float _S1773 = 200.0f * _S1667;
        float _S1774 = _S1652 * _S1773 + 0.5f * (_S1650 * _S1773);
        _S1652 = 0.0f;
        _S1655 = _S1774;
    }
    else
    {
        _S1652 = _S1667;
        _S1655 = 0.0f;
    }
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = _S1650;
    (&_S1775)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1775, _S1652);
    float _S1776 = (_S1775.differential_0 + _S1655) / _S1631;
    FixedArray<float, 22>  _S1777;
    _S1777[int(0)] = 0.0f;
    _S1777[int(1)] = 0.0f;
    _S1777[int(2)] = 0.0f;
    _S1777[int(3)] = 0.0f;
    _S1777[int(4)] = 0.0f;
    _S1777[int(5)] = 0.0f;
    _S1777[int(6)] = 0.0f;
    _S1777[int(7)] = 0.0f;
    _S1777[int(8)] = 0.0f;
    _S1777[int(9)] = 0.0f;
    _S1777[int(10)] = 0.0f;
    _S1777[int(11)] = 0.0f;
    _S1777[int(12)] = 0.0f;
    _S1777[int(13)] = 0.0f;
    _S1777[int(14)] = 0.0f;
    _S1777[int(15)] = 0.0f;
    _S1777[int(16)] = 0.0f;
    _S1777[int(17)] = 0.0f;
    _S1777[int(18)] = 0.0f;
    _S1777[int(19)] = 0.0f;
    _S1777[int(20)] = 0.0f;
    _S1777[int(21)] = 0.0f;
    _S1777[int(14)] = _S1776;
    float _S1778 = _S1751 + _S1777[int(0)];
    float _S1779 = _S1752 + _S1777[int(1)];
    float _S1780 = _S1753 + _S1777[int(2)];
    float _S1781 = _S1754 + _S1777[int(3)];
    float _S1782 = _S1755 + _S1777[int(4)];
    float _S1783 = _S1756 + _S1777[int(5)];
    float _S1784 = _S1757 + _S1777[int(6)];
    float _S1785 = _S1758 + _S1777[int(7)];
    float _S1786 = _S1759 + _S1777[int(8)];
    float _S1787 = _S1760 + _S1777[int(9)];
    float _S1788 = _S1761 + _S1777[int(10)];
    float _S1789 = _S1762 + _S1777[int(11)];
    float _S1790 = _S1763 + _S1777[int(12)];
    float _S1791 = _S1764 + _S1777[int(13)];
    float _S1792 = _S1765 + _S1777[int(14)];
    float _S1793 = _S1766 + _S1777[int(15)];
    float _S1794 = _S1767 + _S1777[int(16)];
    float _S1795 = _S1768 + _S1777[int(17)];
    float _S1796 = _S1769 + _S1777[int(18)];
    float _S1797 = _S1770 + _S1777[int(19)];
    float _S1798 = _S1771 + _S1777[int(20)];
    float _S1799 = _S1772 + _S1777[int(21)];
    if(_S1648)
    {
        float _S1800 = 200.0f * _S1667;
        float _S1801 = _S1649 * _S1800 + 0.5f * (_S1647 * _S1800);
        _S1649 = 0.0f;
        _S1652 = _S1801;
    }
    else
    {
        _S1649 = _S1667;
        _S1652 = 0.0f;
    }
    DiffPair_float_0 _S1802;
    (&_S1802)->primal_0 = _S1647;
    (&_S1802)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1802, _S1649);
    float _S1803 = (_S1802.differential_0 + _S1652) / _S1631;
    FixedArray<float, 22>  _S1804;
    _S1804[int(0)] = 0.0f;
    _S1804[int(1)] = 0.0f;
    _S1804[int(2)] = 0.0f;
    _S1804[int(3)] = 0.0f;
    _S1804[int(4)] = 0.0f;
    _S1804[int(5)] = 0.0f;
    _S1804[int(6)] = 0.0f;
    _S1804[int(7)] = 0.0f;
    _S1804[int(8)] = 0.0f;
    _S1804[int(9)] = 0.0f;
    _S1804[int(10)] = 0.0f;
    _S1804[int(11)] = 0.0f;
    _S1804[int(12)] = 0.0f;
    _S1804[int(13)] = 0.0f;
    _S1804[int(14)] = 0.0f;
    _S1804[int(15)] = 0.0f;
    _S1804[int(16)] = 0.0f;
    _S1804[int(17)] = 0.0f;
    _S1804[int(18)] = 0.0f;
    _S1804[int(19)] = 0.0f;
    _S1804[int(20)] = 0.0f;
    _S1804[int(21)] = 0.0f;
    _S1804[int(13)] = _S1803;
    float _S1805 = _S1778 + _S1804[int(0)];
    float _S1806 = _S1779 + _S1804[int(1)];
    float _S1807 = _S1780 + _S1804[int(2)];
    float _S1808 = _S1781 + _S1804[int(3)];
    float _S1809 = _S1782 + _S1804[int(4)];
    float _S1810 = _S1783 + _S1804[int(5)];
    float _S1811 = _S1784 + _S1804[int(6)];
    float _S1812 = _S1785 + _S1804[int(7)];
    float _S1813 = _S1786 + _S1804[int(8)];
    float _S1814 = _S1787 + _S1804[int(9)];
    float _S1815 = _S1788 + _S1804[int(10)];
    float _S1816 = _S1789 + _S1804[int(11)];
    float _S1817 = _S1790 + _S1804[int(12)];
    float _S1818 = _S1791 + _S1804[int(13)];
    float _S1819 = _S1792 + _S1804[int(14)];
    float _S1820 = _S1793 + _S1804[int(15)];
    float _S1821 = _S1794 + _S1804[int(16)];
    float _S1822 = _S1795 + _S1804[int(17)];
    float _S1823 = _S1796 + _S1804[int(18)];
    float _S1824 = _S1797 + _S1804[int(19)];
    float _S1825 = _S1798 + _S1804[int(20)];
    float _S1826 = _S1799 + _S1804[int(21)];
    if(_S1645)
    {
        float _S1827 = 200.0f * _S1667;
        float _S1828 = _S1646 * _S1827 + 0.5f * (_S1644 * _S1827);
        _S1646 = 0.0f;
        _S1649 = _S1828;
    }
    else
    {
        _S1646 = _S1667;
        _S1649 = 0.0f;
    }
    DiffPair_float_0 _S1829;
    (&_S1829)->primal_0 = _S1644;
    (&_S1829)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1829, _S1646);
    float _S1830 = (_S1829.differential_0 + _S1649) / _S1631;
    FixedArray<float, 22>  _S1831;
    _S1831[int(0)] = 0.0f;
    _S1831[int(1)] = 0.0f;
    _S1831[int(2)] = 0.0f;
    _S1831[int(3)] = 0.0f;
    _S1831[int(4)] = 0.0f;
    _S1831[int(5)] = 0.0f;
    _S1831[int(6)] = 0.0f;
    _S1831[int(7)] = 0.0f;
    _S1831[int(8)] = 0.0f;
    _S1831[int(9)] = 0.0f;
    _S1831[int(10)] = 0.0f;
    _S1831[int(11)] = 0.0f;
    _S1831[int(12)] = 0.0f;
    _S1831[int(13)] = 0.0f;
    _S1831[int(14)] = 0.0f;
    _S1831[int(15)] = 0.0f;
    _S1831[int(16)] = 0.0f;
    _S1831[int(17)] = 0.0f;
    _S1831[int(18)] = 0.0f;
    _S1831[int(19)] = 0.0f;
    _S1831[int(20)] = 0.0f;
    _S1831[int(21)] = 0.0f;
    _S1831[int(12)] = _S1830;
    float _S1832 = _S1805 + _S1831[int(0)];
    float _S1833 = _S1806 + _S1831[int(1)];
    float _S1834 = _S1807 + _S1831[int(2)];
    float _S1835 = _S1808 + _S1831[int(3)];
    float _S1836 = _S1809 + _S1831[int(4)];
    float _S1837 = _S1810 + _S1831[int(5)];
    float _S1838 = _S1811 + _S1831[int(6)];
    float _S1839 = _S1812 + _S1831[int(7)];
    float _S1840 = _S1813 + _S1831[int(8)];
    float _S1841 = _S1814 + _S1831[int(9)];
    float _S1842 = _S1815 + _S1831[int(10)];
    float _S1843 = _S1816 + _S1831[int(11)];
    float _S1844 = _S1817 + _S1831[int(12)];
    float _S1845 = _S1818 + _S1831[int(13)];
    float _S1846 = _S1819 + _S1831[int(14)];
    float _S1847 = _S1820 + _S1831[int(15)];
    float _S1848 = _S1821 + _S1831[int(16)];
    float _S1849 = _S1822 + _S1831[int(17)];
    float _S1850 = _S1823 + _S1831[int(18)];
    float _S1851 = _S1824 + _S1831[int(19)];
    float _S1852 = _S1825 + _S1831[int(20)];
    float _S1853 = _S1826 + _S1831[int(21)];
    if(_S1642)
    {
        float _S1854 = 200.0f * _S1667;
        float _S1855 = _S1643 * _S1854 + 0.5f * (_S1641 * _S1854);
        _S1643 = 0.0f;
        _S1646 = _S1855;
    }
    else
    {
        _S1643 = _S1667;
        _S1646 = 0.0f;
    }
    DiffPair_float_0 _S1856;
    (&_S1856)->primal_0 = _S1641;
    (&_S1856)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1856, _S1643);
    float _S1857 = (_S1856.differential_0 + _S1646) / _S1631;
    FixedArray<float, 22>  _S1858;
    _S1858[int(0)] = 0.0f;
    _S1858[int(1)] = 0.0f;
    _S1858[int(2)] = 0.0f;
    _S1858[int(3)] = 0.0f;
    _S1858[int(4)] = 0.0f;
    _S1858[int(5)] = 0.0f;
    _S1858[int(6)] = 0.0f;
    _S1858[int(7)] = 0.0f;
    _S1858[int(8)] = 0.0f;
    _S1858[int(9)] = 0.0f;
    _S1858[int(10)] = 0.0f;
    _S1858[int(11)] = 0.0f;
    _S1858[int(12)] = 0.0f;
    _S1858[int(13)] = 0.0f;
    _S1858[int(14)] = 0.0f;
    _S1858[int(15)] = 0.0f;
    _S1858[int(16)] = 0.0f;
    _S1858[int(17)] = 0.0f;
    _S1858[int(18)] = 0.0f;
    _S1858[int(19)] = 0.0f;
    _S1858[int(20)] = 0.0f;
    _S1858[int(21)] = 0.0f;
    _S1858[int(11)] = _S1857;
    float _S1859 = _S1832 + _S1858[int(0)];
    float _S1860 = _S1833 + _S1858[int(1)];
    float _S1861 = _S1834 + _S1858[int(2)];
    float _S1862 = _S1835 + _S1858[int(3)];
    float _S1863 = _S1836 + _S1858[int(4)];
    float _S1864 = _S1837 + _S1858[int(5)];
    float _S1865 = _S1838 + _S1858[int(6)];
    float _S1866 = _S1839 + _S1858[int(7)];
    float _S1867 = _S1840 + _S1858[int(8)];
    float _S1868 = _S1841 + _S1858[int(9)];
    float _S1869 = _S1842 + _S1858[int(10)];
    float _S1870 = _S1843 + _S1858[int(11)];
    float _S1871 = _S1844 + _S1858[int(12)];
    float _S1872 = _S1845 + _S1858[int(13)];
    float _S1873 = _S1846 + _S1858[int(14)];
    float _S1874 = _S1847 + _S1858[int(15)];
    float _S1875 = _S1848 + _S1858[int(16)];
    float _S1876 = _S1849 + _S1858[int(17)];
    float _S1877 = _S1850 + _S1858[int(18)];
    float _S1878 = _S1851 + _S1858[int(19)];
    float _S1879 = _S1852 + _S1858[int(20)];
    float _S1880 = _S1853 + _S1858[int(21)];
    if(_S1639)
    {
        float _S1881 = 200.0f * _S1667;
        float _S1882 = _S1640 * _S1881 + 0.5f * (_S1638 * _S1881);
        _S1640 = 0.0f;
        _S1643 = _S1882;
    }
    else
    {
        _S1640 = _S1667;
        _S1643 = 0.0f;
    }
    DiffPair_float_0 _S1883;
    (&_S1883)->primal_0 = _S1638;
    (&_S1883)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1883, _S1640);
    float _S1884 = (_S1883.differential_0 + _S1643) / _S1631;
    float _S1885 = _S1662 / _S1637;
    float _S1886 = _S1663 / _S1636;
    float _S1887 = _S1664 / _S1635;
    FixedArray<float, 22>  _S1888;
    _S1888[int(0)] = 0.0f;
    _S1888[int(1)] = 0.0f;
    _S1888[int(2)] = 0.0f;
    _S1888[int(3)] = 0.0f;
    _S1888[int(4)] = 0.0f;
    _S1888[int(5)] = 0.0f;
    _S1888[int(6)] = 0.0f;
    _S1888[int(7)] = 0.0f;
    _S1888[int(8)] = 0.0f;
    _S1888[int(9)] = 0.0f;
    _S1888[int(10)] = 0.0f;
    _S1888[int(11)] = 0.0f;
    _S1888[int(12)] = 0.0f;
    _S1888[int(13)] = 0.0f;
    _S1888[int(14)] = 0.0f;
    _S1888[int(15)] = 0.0f;
    _S1888[int(16)] = 0.0f;
    _S1888[int(17)] = 0.0f;
    _S1888[int(18)] = 0.0f;
    _S1888[int(19)] = 0.0f;
    _S1888[int(20)] = 0.0f;
    _S1888[int(21)] = 0.0f;
    _S1888[int(10)] = _S1884;
    _S1888[int(9)] = _S1885;
    _S1888[int(8)] = _S1885;
    _S1888[int(7)] = _S1885;
    _S1888[int(6)] = _S1885;
    _S1888[int(5)] = _S1885;
    _S1888[int(4)] = _S1886;
    _S1888[int(3)] = _S1886;
    _S1888[int(2)] = _S1886;
    _S1888[int(1)] = _S1887;
    float _S1889 = _S1859 + _S1888[int(0)];
    float _S1890 = _S1860 + _S1888[int(1)];
    float _S1891 = _S1861 + _S1888[int(2)];
    float _S1892 = _S1862 + _S1888[int(3)];
    float _S1893 = _S1863 + _S1888[int(4)];
    float _S1894 = _S1864 + _S1888[int(5)];
    float _S1895 = _S1865 + _S1888[int(6)];
    float _S1896 = _S1866 + _S1888[int(7)];
    float _S1897 = _S1867 + _S1888[int(8)];
    float _S1898 = _S1868 + _S1888[int(9)];
    float _S1899 = _S1869 + _S1888[int(10)];
    float _S1900 = _S1870 + _S1888[int(11)];
    float _S1901 = _S1871 + _S1888[int(12)];
    float _S1902 = _S1872 + _S1888[int(13)];
    float _S1903 = _S1873 + _S1888[int(14)];
    float _S1904 = _S1874 + _S1888[int(15)];
    float _S1905 = _S1875 + _S1888[int(16)];
    float _S1906 = _S1876 + _S1888[int(17)];
    float _S1907 = _S1877 + _S1888[int(18)];
    float _S1908 = _S1878 + _S1888[int(19)];
    float _S1909 = _S1879 + _S1888[int(20)];
    float _S1910 = _S1880 + _S1888[int(21)];
    if(_S1633)
    {
        float _S1911 = 10.0f * _S1665;
        float _S1912 = _S1634 * _S1911 + 0.5f * (_S1632 * _S1911);
        _S1634 = 0.0f;
        _S1640 = _S1912;
    }
    else
    {
        _S1634 = _S1665;
        _S1640 = 0.0f;
    }
    DiffPair_float_0 _S1913;
    (&_S1913)->primal_0 = _S1632;
    (&_S1913)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1913, _S1634);
    float _S1914 = (_S1913.differential_0 + _S1640) / _S1631;
    FixedArray<float, 22>  _S1915;
    _S1915[int(0)] = 0.0f;
    _S1915[int(1)] = 0.0f;
    _S1915[int(2)] = 0.0f;
    _S1915[int(3)] = 0.0f;
    _S1915[int(4)] = 0.0f;
    _S1915[int(5)] = 0.0f;
    _S1915[int(6)] = 0.0f;
    _S1915[int(7)] = 0.0f;
    _S1915[int(8)] = 0.0f;
    _S1915[int(9)] = 0.0f;
    _S1915[int(10)] = 0.0f;
    _S1915[int(11)] = 0.0f;
    _S1915[int(12)] = 0.0f;
    _S1915[int(13)] = 0.0f;
    _S1915[int(14)] = 0.0f;
    _S1915[int(15)] = 0.0f;
    _S1915[int(16)] = 0.0f;
    _S1915[int(17)] = 0.0f;
    _S1915[int(18)] = 0.0f;
    _S1915[int(19)] = 0.0f;
    _S1915[int(20)] = 0.0f;
    _S1915[int(21)] = 0.0f;
    _S1915[int(0)] = _S1914;
    FixedArray<float, 22>  _S1916 = {
        _S1889 + _S1915[int(0)], _S1890 + _S1915[int(1)], _S1891 + _S1915[int(2)], _S1892 + _S1915[int(3)], _S1893 + _S1915[int(4)], _S1894 + _S1915[int(5)], _S1895 + _S1915[int(6)], _S1896 + _S1915[int(7)], _S1897 + _S1915[int(8)], _S1898 + _S1915[int(9)], _S1899 + _S1915[int(10)], _S1900 + _S1915[int(11)], _S1901 + _S1915[int(12)], _S1902 + _S1915[int(13)], _S1903 + _S1915[int(14)], _S1904 + _S1915[int(15)], _S1905 + _S1915[int(16)], _S1906 + _S1915[int(17)], _S1907 + _S1915[int(18)], _S1908 + _S1915[int(19)], _S1909 + _S1915[int(20)], _S1910 + _S1915[int(21)]
    };
    dpraw_losses_0->primal_0 = dpraw_losses_0->primal_0;
    dpraw_losses_0->differential_0 = _S1916;
    return;
}

inline __device__ void s_bwd_compute_ppisp_regularization_loss_0(DiffPair_arrayx3Cfloatx2C22x3E_0 * _S1917, int _S1918, FixedArray<float, 6>  * _S1919, FixedArray<float, 6>  * _S1920)
{
    s_bwd_prop_compute_ppisp_regularization_loss_0(_S1917, _S1918, _S1919, _S1920);
    return;
}

inline __device__ void compute_ppisp_regularization_loss_vjp(FixedArray<float, 22>  raw_losses_2, int num_cameras_3, FixedArray<float, 6>  loss_weights_3, FixedArray<float, 6>  grad_out_4, FixedArray<float, 22>  * _S1921)
{
    FixedArray<float, 22>  _S1922 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C22x3E_0 dp_raw_losses_0;
    (&dp_raw_losses_0)->primal_0 = raw_losses_2;
    (&dp_raw_losses_0)->differential_0 = _S1922;
    FixedArray<float, 6>  _S1923 = loss_weights_3;
    FixedArray<float, 6>  _S1924 = grad_out_4;
    s_bwd_compute_ppisp_regularization_loss_0(&dp_raw_losses_0, num_cameras_3, &_S1923, &_S1924);
    *_S1921 = (&dp_raw_losses_0)->differential_0;
    return;
}

struct DiffPair_arrayx3Cfloatx2C23x3E_0
{
    FixedArray<float, 23>  primal_0;
    FixedArray<float, 23>  differential_0;
};

inline __device__ void s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * dpraw_losses_1, int num_cameras_4, FixedArray<float, 6>  * loss_weights_4, FixedArray<float, 6>  * _s_dOut_5)
{
    FixedArray<float, 23>  _S1925 = dpraw_losses_1->primal_0;
    float _S1926 = float(num_cameras_4);
    float _S1927 = dpraw_losses_1->primal_0[int(0)] / _S1926;
    bool _S1928 = (s_primal_ctx_abs_0(_S1927)) < 0.10000000149011612f;
    float _S1929;
    if(_S1928)
    {
        _S1929 = 0.5f * _S1927;
    }
    else
    {
        _S1929 = 0.0f;
    }
    float _S1930 = 3.0f * _S1926;
    float _S1931 = 9.0f * _S1926;
    float _S1932 = 5.0f * _S1926;
    float _S1933 = _S1925[int(10)] / _S1926;
    bool _S1934 = (s_primal_ctx_abs_0(_S1933)) < 0.00499999988824129f;
    float _S1935;
    if(_S1934)
    {
        _S1935 = 0.5f * _S1933;
    }
    else
    {
        _S1935 = 0.0f;
    }
    float _S1936 = _S1925[int(11)] / _S1926;
    bool _S1937 = (s_primal_ctx_abs_0(_S1936)) < 0.00499999988824129f;
    float _S1938;
    if(_S1937)
    {
        _S1938 = 0.5f * _S1936;
    }
    else
    {
        _S1938 = 0.0f;
    }
    float _S1939 = _S1925[int(12)] / _S1926;
    bool _S1940 = (s_primal_ctx_abs_0(_S1939)) < 0.00499999988824129f;
    float _S1941;
    if(_S1940)
    {
        _S1941 = 0.5f * _S1939;
    }
    else
    {
        _S1941 = 0.0f;
    }
    float _S1942 = _S1925[int(13)] / _S1926;
    bool _S1943 = (s_primal_ctx_abs_0(_S1942)) < 0.00499999988824129f;
    float _S1944;
    if(_S1943)
    {
        _S1944 = 0.5f * _S1942;
    }
    else
    {
        _S1944 = 0.0f;
    }
    float _S1945 = _S1925[int(14)] / _S1926;
    bool _S1946 = (s_primal_ctx_abs_0(_S1945)) < 0.00499999988824129f;
    float _S1947;
    if(_S1946)
    {
        _S1947 = 0.5f * _S1945;
    }
    else
    {
        _S1947 = 0.0f;
    }
    float _S1948 = _S1925[int(15)] / _S1926;
    bool _S1949 = (s_primal_ctx_abs_0(_S1948)) < 0.00499999988824129f;
    float _S1950;
    if(_S1949)
    {
        _S1950 = 0.5f * _S1948;
    }
    else
    {
        _S1950 = 0.0f;
    }
    float _S1951 = _S1925[int(16)] / _S1926;
    bool _S1952 = (s_primal_ctx_abs_0(_S1951)) < 0.00499999988824129f;
    float _S1953;
    if(_S1952)
    {
        _S1953 = 0.5f * _S1951;
    }
    else
    {
        _S1953 = 0.0f;
    }
    float _S1954 = _S1925[int(17)] / _S1926;
    bool _S1955 = (s_primal_ctx_abs_0(_S1954)) < 0.00499999988824129f;
    float _S1956;
    if(_S1955)
    {
        _S1956 = 0.5f * _S1954;
    }
    else
    {
        _S1956 = 0.0f;
    }
    float _S1957 = (*loss_weights_4)[int(3)] * (*_s_dOut_5)[int(3)];
    float _S1958 = (*loss_weights_4)[int(2)] * (*_s_dOut_5)[int(2)];
    float _S1959 = (*loss_weights_4)[int(1)] * (*_s_dOut_5)[int(1)];
    float _S1960 = (*loss_weights_4)[int(0)] * (*_s_dOut_5)[int(0)];
    float _S1961 = (*loss_weights_4)[int(5)] * (*_s_dOut_5)[int(5)] / _S1932;
    float _S1962 = 0.125f * ((*loss_weights_4)[int(4)] * (*_s_dOut_5)[int(4)]);
    FixedArray<float, 23>  _S1963;
    _S1963[int(0)] = 0.0f;
    _S1963[int(1)] = 0.0f;
    _S1963[int(2)] = 0.0f;
    _S1963[int(3)] = 0.0f;
    _S1963[int(4)] = 0.0f;
    _S1963[int(5)] = 0.0f;
    _S1963[int(6)] = 0.0f;
    _S1963[int(7)] = 0.0f;
    _S1963[int(8)] = 0.0f;
    _S1963[int(9)] = 0.0f;
    _S1963[int(10)] = 0.0f;
    _S1963[int(11)] = 0.0f;
    _S1963[int(12)] = 0.0f;
    _S1963[int(13)] = 0.0f;
    _S1963[int(14)] = 0.0f;
    _S1963[int(15)] = 0.0f;
    _S1963[int(16)] = 0.0f;
    _S1963[int(17)] = 0.0f;
    _S1963[int(18)] = 0.0f;
    _S1963[int(19)] = 0.0f;
    _S1963[int(20)] = 0.0f;
    _S1963[int(21)] = 0.0f;
    _S1963[int(22)] = 0.0f;
    _S1963[int(22)] = _S1961;
    _S1963[int(21)] = _S1961;
    _S1963[int(20)] = _S1961;
    _S1963[int(19)] = _S1961;
    _S1963[int(18)] = _S1961;
    float _S1964 = _S1963[int(0)];
    float _S1965 = _S1963[int(1)];
    float _S1966 = _S1963[int(2)];
    float _S1967 = _S1963[int(3)];
    float _S1968 = _S1963[int(4)];
    float _S1969 = _S1963[int(5)];
    float _S1970 = _S1963[int(6)];
    float _S1971 = _S1963[int(7)];
    float _S1972 = _S1963[int(8)];
    float _S1973 = _S1963[int(9)];
    float _S1974 = _S1963[int(10)];
    float _S1975 = _S1963[int(11)];
    float _S1976 = _S1963[int(12)];
    float _S1977 = _S1963[int(13)];
    float _S1978 = _S1963[int(14)];
    float _S1979 = _S1963[int(15)];
    float _S1980 = _S1963[int(16)];
    float _S1981 = _S1963[int(17)];
    float _S1982 = _S1963[int(18)];
    float _S1983 = _S1963[int(19)];
    float _S1984 = _S1963[int(20)];
    float _S1985 = _S1963[int(21)];
    float _S1986 = _S1963[int(22)];
    float _S1987;
    if(_S1955)
    {
        float _S1988 = 200.0f * _S1962;
        float _S1989 = _S1956 * _S1988 + 0.5f * (_S1954 * _S1988);
        _S1956 = 0.0f;
        _S1987 = _S1989;
    }
    else
    {
        _S1956 = _S1962;
        _S1987 = 0.0f;
    }
    DiffPair_float_0 _S1990;
    (&_S1990)->primal_0 = _S1954;
    (&_S1990)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1990, _S1956);
    float _S1991 = (_S1990.differential_0 + _S1987) / _S1926;
    FixedArray<float, 23>  _S1992;
    _S1992[int(0)] = 0.0f;
    _S1992[int(1)] = 0.0f;
    _S1992[int(2)] = 0.0f;
    _S1992[int(3)] = 0.0f;
    _S1992[int(4)] = 0.0f;
    _S1992[int(5)] = 0.0f;
    _S1992[int(6)] = 0.0f;
    _S1992[int(7)] = 0.0f;
    _S1992[int(8)] = 0.0f;
    _S1992[int(9)] = 0.0f;
    _S1992[int(10)] = 0.0f;
    _S1992[int(11)] = 0.0f;
    _S1992[int(12)] = 0.0f;
    _S1992[int(13)] = 0.0f;
    _S1992[int(14)] = 0.0f;
    _S1992[int(15)] = 0.0f;
    _S1992[int(16)] = 0.0f;
    _S1992[int(17)] = 0.0f;
    _S1992[int(18)] = 0.0f;
    _S1992[int(19)] = 0.0f;
    _S1992[int(20)] = 0.0f;
    _S1992[int(21)] = 0.0f;
    _S1992[int(22)] = 0.0f;
    _S1992[int(17)] = _S1991;
    float _S1993 = _S1964 + _S1992[int(0)];
    float _S1994 = _S1965 + _S1992[int(1)];
    float _S1995 = _S1966 + _S1992[int(2)];
    float _S1996 = _S1967 + _S1992[int(3)];
    float _S1997 = _S1968 + _S1992[int(4)];
    float _S1998 = _S1969 + _S1992[int(5)];
    float _S1999 = _S1970 + _S1992[int(6)];
    float _S2000 = _S1971 + _S1992[int(7)];
    float _S2001 = _S1972 + _S1992[int(8)];
    float _S2002 = _S1973 + _S1992[int(9)];
    float _S2003 = _S1974 + _S1992[int(10)];
    float _S2004 = _S1975 + _S1992[int(11)];
    float _S2005 = _S1976 + _S1992[int(12)];
    float _S2006 = _S1977 + _S1992[int(13)];
    float _S2007 = _S1978 + _S1992[int(14)];
    float _S2008 = _S1979 + _S1992[int(15)];
    float _S2009 = _S1980 + _S1992[int(16)];
    float _S2010 = _S1981 + _S1992[int(17)];
    float _S2011 = _S1982 + _S1992[int(18)];
    float _S2012 = _S1983 + _S1992[int(19)];
    float _S2013 = _S1984 + _S1992[int(20)];
    float _S2014 = _S1985 + _S1992[int(21)];
    float _S2015 = _S1986 + _S1992[int(22)];
    if(_S1952)
    {
        float _S2016 = 200.0f * _S1962;
        float _S2017 = _S1953 * _S2016 + 0.5f * (_S1951 * _S2016);
        _S1953 = 0.0f;
        _S1956 = _S2017;
    }
    else
    {
        _S1953 = _S1962;
        _S1956 = 0.0f;
    }
    DiffPair_float_0 _S2018;
    (&_S2018)->primal_0 = _S1951;
    (&_S2018)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2018, _S1953);
    float _S2019 = (_S2018.differential_0 + _S1956) / _S1926;
    FixedArray<float, 23>  _S2020;
    _S2020[int(0)] = 0.0f;
    _S2020[int(1)] = 0.0f;
    _S2020[int(2)] = 0.0f;
    _S2020[int(3)] = 0.0f;
    _S2020[int(4)] = 0.0f;
    _S2020[int(5)] = 0.0f;
    _S2020[int(6)] = 0.0f;
    _S2020[int(7)] = 0.0f;
    _S2020[int(8)] = 0.0f;
    _S2020[int(9)] = 0.0f;
    _S2020[int(10)] = 0.0f;
    _S2020[int(11)] = 0.0f;
    _S2020[int(12)] = 0.0f;
    _S2020[int(13)] = 0.0f;
    _S2020[int(14)] = 0.0f;
    _S2020[int(15)] = 0.0f;
    _S2020[int(16)] = 0.0f;
    _S2020[int(17)] = 0.0f;
    _S2020[int(18)] = 0.0f;
    _S2020[int(19)] = 0.0f;
    _S2020[int(20)] = 0.0f;
    _S2020[int(21)] = 0.0f;
    _S2020[int(22)] = 0.0f;
    _S2020[int(16)] = _S2019;
    float _S2021 = _S1993 + _S2020[int(0)];
    float _S2022 = _S1994 + _S2020[int(1)];
    float _S2023 = _S1995 + _S2020[int(2)];
    float _S2024 = _S1996 + _S2020[int(3)];
    float _S2025 = _S1997 + _S2020[int(4)];
    float _S2026 = _S1998 + _S2020[int(5)];
    float _S2027 = _S1999 + _S2020[int(6)];
    float _S2028 = _S2000 + _S2020[int(7)];
    float _S2029 = _S2001 + _S2020[int(8)];
    float _S2030 = _S2002 + _S2020[int(9)];
    float _S2031 = _S2003 + _S2020[int(10)];
    float _S2032 = _S2004 + _S2020[int(11)];
    float _S2033 = _S2005 + _S2020[int(12)];
    float _S2034 = _S2006 + _S2020[int(13)];
    float _S2035 = _S2007 + _S2020[int(14)];
    float _S2036 = _S2008 + _S2020[int(15)];
    float _S2037 = _S2009 + _S2020[int(16)];
    float _S2038 = _S2010 + _S2020[int(17)];
    float _S2039 = _S2011 + _S2020[int(18)];
    float _S2040 = _S2012 + _S2020[int(19)];
    float _S2041 = _S2013 + _S2020[int(20)];
    float _S2042 = _S2014 + _S2020[int(21)];
    float _S2043 = _S2015 + _S2020[int(22)];
    if(_S1949)
    {
        float _S2044 = 200.0f * _S1962;
        float _S2045 = _S1950 * _S2044 + 0.5f * (_S1948 * _S2044);
        _S1950 = 0.0f;
        _S1953 = _S2045;
    }
    else
    {
        _S1950 = _S1962;
        _S1953 = 0.0f;
    }
    DiffPair_float_0 _S2046;
    (&_S2046)->primal_0 = _S1948;
    (&_S2046)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2046, _S1950);
    float _S2047 = (_S2046.differential_0 + _S1953) / _S1926;
    FixedArray<float, 23>  _S2048;
    _S2048[int(0)] = 0.0f;
    _S2048[int(1)] = 0.0f;
    _S2048[int(2)] = 0.0f;
    _S2048[int(3)] = 0.0f;
    _S2048[int(4)] = 0.0f;
    _S2048[int(5)] = 0.0f;
    _S2048[int(6)] = 0.0f;
    _S2048[int(7)] = 0.0f;
    _S2048[int(8)] = 0.0f;
    _S2048[int(9)] = 0.0f;
    _S2048[int(10)] = 0.0f;
    _S2048[int(11)] = 0.0f;
    _S2048[int(12)] = 0.0f;
    _S2048[int(13)] = 0.0f;
    _S2048[int(14)] = 0.0f;
    _S2048[int(15)] = 0.0f;
    _S2048[int(16)] = 0.0f;
    _S2048[int(17)] = 0.0f;
    _S2048[int(18)] = 0.0f;
    _S2048[int(19)] = 0.0f;
    _S2048[int(20)] = 0.0f;
    _S2048[int(21)] = 0.0f;
    _S2048[int(22)] = 0.0f;
    _S2048[int(15)] = _S2047;
    float _S2049 = _S2021 + _S2048[int(0)];
    float _S2050 = _S2022 + _S2048[int(1)];
    float _S2051 = _S2023 + _S2048[int(2)];
    float _S2052 = _S2024 + _S2048[int(3)];
    float _S2053 = _S2025 + _S2048[int(4)];
    float _S2054 = _S2026 + _S2048[int(5)];
    float _S2055 = _S2027 + _S2048[int(6)];
    float _S2056 = _S2028 + _S2048[int(7)];
    float _S2057 = _S2029 + _S2048[int(8)];
    float _S2058 = _S2030 + _S2048[int(9)];
    float _S2059 = _S2031 + _S2048[int(10)];
    float _S2060 = _S2032 + _S2048[int(11)];
    float _S2061 = _S2033 + _S2048[int(12)];
    float _S2062 = _S2034 + _S2048[int(13)];
    float _S2063 = _S2035 + _S2048[int(14)];
    float _S2064 = _S2036 + _S2048[int(15)];
    float _S2065 = _S2037 + _S2048[int(16)];
    float _S2066 = _S2038 + _S2048[int(17)];
    float _S2067 = _S2039 + _S2048[int(18)];
    float _S2068 = _S2040 + _S2048[int(19)];
    float _S2069 = _S2041 + _S2048[int(20)];
    float _S2070 = _S2042 + _S2048[int(21)];
    float _S2071 = _S2043 + _S2048[int(22)];
    if(_S1946)
    {
        float _S2072 = 200.0f * _S1962;
        float _S2073 = _S1947 * _S2072 + 0.5f * (_S1945 * _S2072);
        _S1947 = 0.0f;
        _S1950 = _S2073;
    }
    else
    {
        _S1947 = _S1962;
        _S1950 = 0.0f;
    }
    DiffPair_float_0 _S2074;
    (&_S2074)->primal_0 = _S1945;
    (&_S2074)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2074, _S1947);
    float _S2075 = (_S2074.differential_0 + _S1950) / _S1926;
    FixedArray<float, 23>  _S2076;
    _S2076[int(0)] = 0.0f;
    _S2076[int(1)] = 0.0f;
    _S2076[int(2)] = 0.0f;
    _S2076[int(3)] = 0.0f;
    _S2076[int(4)] = 0.0f;
    _S2076[int(5)] = 0.0f;
    _S2076[int(6)] = 0.0f;
    _S2076[int(7)] = 0.0f;
    _S2076[int(8)] = 0.0f;
    _S2076[int(9)] = 0.0f;
    _S2076[int(10)] = 0.0f;
    _S2076[int(11)] = 0.0f;
    _S2076[int(12)] = 0.0f;
    _S2076[int(13)] = 0.0f;
    _S2076[int(14)] = 0.0f;
    _S2076[int(15)] = 0.0f;
    _S2076[int(16)] = 0.0f;
    _S2076[int(17)] = 0.0f;
    _S2076[int(18)] = 0.0f;
    _S2076[int(19)] = 0.0f;
    _S2076[int(20)] = 0.0f;
    _S2076[int(21)] = 0.0f;
    _S2076[int(22)] = 0.0f;
    _S2076[int(14)] = _S2075;
    float _S2077 = _S2049 + _S2076[int(0)];
    float _S2078 = _S2050 + _S2076[int(1)];
    float _S2079 = _S2051 + _S2076[int(2)];
    float _S2080 = _S2052 + _S2076[int(3)];
    float _S2081 = _S2053 + _S2076[int(4)];
    float _S2082 = _S2054 + _S2076[int(5)];
    float _S2083 = _S2055 + _S2076[int(6)];
    float _S2084 = _S2056 + _S2076[int(7)];
    float _S2085 = _S2057 + _S2076[int(8)];
    float _S2086 = _S2058 + _S2076[int(9)];
    float _S2087 = _S2059 + _S2076[int(10)];
    float _S2088 = _S2060 + _S2076[int(11)];
    float _S2089 = _S2061 + _S2076[int(12)];
    float _S2090 = _S2062 + _S2076[int(13)];
    float _S2091 = _S2063 + _S2076[int(14)];
    float _S2092 = _S2064 + _S2076[int(15)];
    float _S2093 = _S2065 + _S2076[int(16)];
    float _S2094 = _S2066 + _S2076[int(17)];
    float _S2095 = _S2067 + _S2076[int(18)];
    float _S2096 = _S2068 + _S2076[int(19)];
    float _S2097 = _S2069 + _S2076[int(20)];
    float _S2098 = _S2070 + _S2076[int(21)];
    float _S2099 = _S2071 + _S2076[int(22)];
    if(_S1943)
    {
        float _S2100 = 200.0f * _S1962;
        float _S2101 = _S1944 * _S2100 + 0.5f * (_S1942 * _S2100);
        _S1944 = 0.0f;
        _S1947 = _S2101;
    }
    else
    {
        _S1944 = _S1962;
        _S1947 = 0.0f;
    }
    DiffPair_float_0 _S2102;
    (&_S2102)->primal_0 = _S1942;
    (&_S2102)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2102, _S1944);
    float _S2103 = (_S2102.differential_0 + _S1947) / _S1926;
    FixedArray<float, 23>  _S2104;
    _S2104[int(0)] = 0.0f;
    _S2104[int(1)] = 0.0f;
    _S2104[int(2)] = 0.0f;
    _S2104[int(3)] = 0.0f;
    _S2104[int(4)] = 0.0f;
    _S2104[int(5)] = 0.0f;
    _S2104[int(6)] = 0.0f;
    _S2104[int(7)] = 0.0f;
    _S2104[int(8)] = 0.0f;
    _S2104[int(9)] = 0.0f;
    _S2104[int(10)] = 0.0f;
    _S2104[int(11)] = 0.0f;
    _S2104[int(12)] = 0.0f;
    _S2104[int(13)] = 0.0f;
    _S2104[int(14)] = 0.0f;
    _S2104[int(15)] = 0.0f;
    _S2104[int(16)] = 0.0f;
    _S2104[int(17)] = 0.0f;
    _S2104[int(18)] = 0.0f;
    _S2104[int(19)] = 0.0f;
    _S2104[int(20)] = 0.0f;
    _S2104[int(21)] = 0.0f;
    _S2104[int(22)] = 0.0f;
    _S2104[int(13)] = _S2103;
    float _S2105 = _S2077 + _S2104[int(0)];
    float _S2106 = _S2078 + _S2104[int(1)];
    float _S2107 = _S2079 + _S2104[int(2)];
    float _S2108 = _S2080 + _S2104[int(3)];
    float _S2109 = _S2081 + _S2104[int(4)];
    float _S2110 = _S2082 + _S2104[int(5)];
    float _S2111 = _S2083 + _S2104[int(6)];
    float _S2112 = _S2084 + _S2104[int(7)];
    float _S2113 = _S2085 + _S2104[int(8)];
    float _S2114 = _S2086 + _S2104[int(9)];
    float _S2115 = _S2087 + _S2104[int(10)];
    float _S2116 = _S2088 + _S2104[int(11)];
    float _S2117 = _S2089 + _S2104[int(12)];
    float _S2118 = _S2090 + _S2104[int(13)];
    float _S2119 = _S2091 + _S2104[int(14)];
    float _S2120 = _S2092 + _S2104[int(15)];
    float _S2121 = _S2093 + _S2104[int(16)];
    float _S2122 = _S2094 + _S2104[int(17)];
    float _S2123 = _S2095 + _S2104[int(18)];
    float _S2124 = _S2096 + _S2104[int(19)];
    float _S2125 = _S2097 + _S2104[int(20)];
    float _S2126 = _S2098 + _S2104[int(21)];
    float _S2127 = _S2099 + _S2104[int(22)];
    if(_S1940)
    {
        float _S2128 = 200.0f * _S1962;
        float _S2129 = _S1941 * _S2128 + 0.5f * (_S1939 * _S2128);
        _S1941 = 0.0f;
        _S1944 = _S2129;
    }
    else
    {
        _S1941 = _S1962;
        _S1944 = 0.0f;
    }
    DiffPair_float_0 _S2130;
    (&_S2130)->primal_0 = _S1939;
    (&_S2130)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2130, _S1941);
    float _S2131 = (_S2130.differential_0 + _S1944) / _S1926;
    FixedArray<float, 23>  _S2132;
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
    _S2132[int(22)] = 0.0f;
    _S2132[int(12)] = _S2131;
    float _S2133 = _S2105 + _S2132[int(0)];
    float _S2134 = _S2106 + _S2132[int(1)];
    float _S2135 = _S2107 + _S2132[int(2)];
    float _S2136 = _S2108 + _S2132[int(3)];
    float _S2137 = _S2109 + _S2132[int(4)];
    float _S2138 = _S2110 + _S2132[int(5)];
    float _S2139 = _S2111 + _S2132[int(6)];
    float _S2140 = _S2112 + _S2132[int(7)];
    float _S2141 = _S2113 + _S2132[int(8)];
    float _S2142 = _S2114 + _S2132[int(9)];
    float _S2143 = _S2115 + _S2132[int(10)];
    float _S2144 = _S2116 + _S2132[int(11)];
    float _S2145 = _S2117 + _S2132[int(12)];
    float _S2146 = _S2118 + _S2132[int(13)];
    float _S2147 = _S2119 + _S2132[int(14)];
    float _S2148 = _S2120 + _S2132[int(15)];
    float _S2149 = _S2121 + _S2132[int(16)];
    float _S2150 = _S2122 + _S2132[int(17)];
    float _S2151 = _S2123 + _S2132[int(18)];
    float _S2152 = _S2124 + _S2132[int(19)];
    float _S2153 = _S2125 + _S2132[int(20)];
    float _S2154 = _S2126 + _S2132[int(21)];
    float _S2155 = _S2127 + _S2132[int(22)];
    if(_S1937)
    {
        float _S2156 = 200.0f * _S1962;
        float _S2157 = _S1938 * _S2156 + 0.5f * (_S1936 * _S2156);
        _S1938 = 0.0f;
        _S1941 = _S2157;
    }
    else
    {
        _S1938 = _S1962;
        _S1941 = 0.0f;
    }
    DiffPair_float_0 _S2158;
    (&_S2158)->primal_0 = _S1936;
    (&_S2158)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2158, _S1938);
    float _S2159 = (_S2158.differential_0 + _S1941) / _S1926;
    FixedArray<float, 23>  _S2160;
    _S2160[int(0)] = 0.0f;
    _S2160[int(1)] = 0.0f;
    _S2160[int(2)] = 0.0f;
    _S2160[int(3)] = 0.0f;
    _S2160[int(4)] = 0.0f;
    _S2160[int(5)] = 0.0f;
    _S2160[int(6)] = 0.0f;
    _S2160[int(7)] = 0.0f;
    _S2160[int(8)] = 0.0f;
    _S2160[int(9)] = 0.0f;
    _S2160[int(10)] = 0.0f;
    _S2160[int(11)] = 0.0f;
    _S2160[int(12)] = 0.0f;
    _S2160[int(13)] = 0.0f;
    _S2160[int(14)] = 0.0f;
    _S2160[int(15)] = 0.0f;
    _S2160[int(16)] = 0.0f;
    _S2160[int(17)] = 0.0f;
    _S2160[int(18)] = 0.0f;
    _S2160[int(19)] = 0.0f;
    _S2160[int(20)] = 0.0f;
    _S2160[int(21)] = 0.0f;
    _S2160[int(22)] = 0.0f;
    _S2160[int(11)] = _S2159;
    float _S2161 = _S2133 + _S2160[int(0)];
    float _S2162 = _S2134 + _S2160[int(1)];
    float _S2163 = _S2135 + _S2160[int(2)];
    float _S2164 = _S2136 + _S2160[int(3)];
    float _S2165 = _S2137 + _S2160[int(4)];
    float _S2166 = _S2138 + _S2160[int(5)];
    float _S2167 = _S2139 + _S2160[int(6)];
    float _S2168 = _S2140 + _S2160[int(7)];
    float _S2169 = _S2141 + _S2160[int(8)];
    float _S2170 = _S2142 + _S2160[int(9)];
    float _S2171 = _S2143 + _S2160[int(10)];
    float _S2172 = _S2144 + _S2160[int(11)];
    float _S2173 = _S2145 + _S2160[int(12)];
    float _S2174 = _S2146 + _S2160[int(13)];
    float _S2175 = _S2147 + _S2160[int(14)];
    float _S2176 = _S2148 + _S2160[int(15)];
    float _S2177 = _S2149 + _S2160[int(16)];
    float _S2178 = _S2150 + _S2160[int(17)];
    float _S2179 = _S2151 + _S2160[int(18)];
    float _S2180 = _S2152 + _S2160[int(19)];
    float _S2181 = _S2153 + _S2160[int(20)];
    float _S2182 = _S2154 + _S2160[int(21)];
    float _S2183 = _S2155 + _S2160[int(22)];
    if(_S1934)
    {
        float _S2184 = 200.0f * _S1962;
        float _S2185 = _S1935 * _S2184 + 0.5f * (_S1933 * _S2184);
        _S1935 = 0.0f;
        _S1938 = _S2185;
    }
    else
    {
        _S1935 = _S1962;
        _S1938 = 0.0f;
    }
    DiffPair_float_0 _S2186;
    (&_S2186)->primal_0 = _S1933;
    (&_S2186)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2186, _S1935);
    float _S2187 = (_S2186.differential_0 + _S1938) / _S1926;
    float _S2188 = _S1957 / _S1932;
    float _S2189 = _S1958 / _S1931;
    float _S2190 = _S1959 / _S1930;
    FixedArray<float, 23>  _S2191;
    _S2191[int(0)] = 0.0f;
    _S2191[int(1)] = 0.0f;
    _S2191[int(2)] = 0.0f;
    _S2191[int(3)] = 0.0f;
    _S2191[int(4)] = 0.0f;
    _S2191[int(5)] = 0.0f;
    _S2191[int(6)] = 0.0f;
    _S2191[int(7)] = 0.0f;
    _S2191[int(8)] = 0.0f;
    _S2191[int(9)] = 0.0f;
    _S2191[int(10)] = 0.0f;
    _S2191[int(11)] = 0.0f;
    _S2191[int(12)] = 0.0f;
    _S2191[int(13)] = 0.0f;
    _S2191[int(14)] = 0.0f;
    _S2191[int(15)] = 0.0f;
    _S2191[int(16)] = 0.0f;
    _S2191[int(17)] = 0.0f;
    _S2191[int(18)] = 0.0f;
    _S2191[int(19)] = 0.0f;
    _S2191[int(20)] = 0.0f;
    _S2191[int(21)] = 0.0f;
    _S2191[int(22)] = 0.0f;
    _S2191[int(10)] = _S2187;
    _S2191[int(9)] = _S2188;
    _S2191[int(8)] = _S2188;
    _S2191[int(7)] = _S2188;
    _S2191[int(6)] = _S2188;
    _S2191[int(5)] = _S2188;
    _S2191[int(4)] = _S2189;
    _S2191[int(3)] = _S2189;
    _S2191[int(2)] = _S2189;
    _S2191[int(1)] = _S2190;
    float _S2192 = _S2161 + _S2191[int(0)];
    float _S2193 = _S2162 + _S2191[int(1)];
    float _S2194 = _S2163 + _S2191[int(2)];
    float _S2195 = _S2164 + _S2191[int(3)];
    float _S2196 = _S2165 + _S2191[int(4)];
    float _S2197 = _S2166 + _S2191[int(5)];
    float _S2198 = _S2167 + _S2191[int(6)];
    float _S2199 = _S2168 + _S2191[int(7)];
    float _S2200 = _S2169 + _S2191[int(8)];
    float _S2201 = _S2170 + _S2191[int(9)];
    float _S2202 = _S2171 + _S2191[int(10)];
    float _S2203 = _S2172 + _S2191[int(11)];
    float _S2204 = _S2173 + _S2191[int(12)];
    float _S2205 = _S2174 + _S2191[int(13)];
    float _S2206 = _S2175 + _S2191[int(14)];
    float _S2207 = _S2176 + _S2191[int(15)];
    float _S2208 = _S2177 + _S2191[int(16)];
    float _S2209 = _S2178 + _S2191[int(17)];
    float _S2210 = _S2179 + _S2191[int(18)];
    float _S2211 = _S2180 + _S2191[int(19)];
    float _S2212 = _S2181 + _S2191[int(20)];
    float _S2213 = _S2182 + _S2191[int(21)];
    float _S2214 = _S2183 + _S2191[int(22)];
    if(_S1928)
    {
        float _S2215 = 10.0f * _S1960;
        float _S2216 = _S1929 * _S2215 + 0.5f * (_S1927 * _S2215);
        _S1929 = 0.0f;
        _S1935 = _S2216;
    }
    else
    {
        _S1929 = _S1960;
        _S1935 = 0.0f;
    }
    DiffPair_float_0 _S2217;
    (&_S2217)->primal_0 = _S1927;
    (&_S2217)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2217, _S1929);
    float _S2218 = (_S2217.differential_0 + _S1935) / _S1926;
    FixedArray<float, 23>  _S2219;
    _S2219[int(0)] = 0.0f;
    _S2219[int(1)] = 0.0f;
    _S2219[int(2)] = 0.0f;
    _S2219[int(3)] = 0.0f;
    _S2219[int(4)] = 0.0f;
    _S2219[int(5)] = 0.0f;
    _S2219[int(6)] = 0.0f;
    _S2219[int(7)] = 0.0f;
    _S2219[int(8)] = 0.0f;
    _S2219[int(9)] = 0.0f;
    _S2219[int(10)] = 0.0f;
    _S2219[int(11)] = 0.0f;
    _S2219[int(12)] = 0.0f;
    _S2219[int(13)] = 0.0f;
    _S2219[int(14)] = 0.0f;
    _S2219[int(15)] = 0.0f;
    _S2219[int(16)] = 0.0f;
    _S2219[int(17)] = 0.0f;
    _S2219[int(18)] = 0.0f;
    _S2219[int(19)] = 0.0f;
    _S2219[int(20)] = 0.0f;
    _S2219[int(21)] = 0.0f;
    _S2219[int(22)] = 0.0f;
    _S2219[int(0)] = _S2218;
    FixedArray<float, 23>  _S2220 = {
        _S2192 + _S2219[int(0)], _S2193 + _S2219[int(1)], _S2194 + _S2219[int(2)], _S2195 + _S2219[int(3)], _S2196 + _S2219[int(4)], _S2197 + _S2219[int(5)], _S2198 + _S2219[int(6)], _S2199 + _S2219[int(7)], _S2200 + _S2219[int(8)], _S2201 + _S2219[int(9)], _S2202 + _S2219[int(10)], _S2203 + _S2219[int(11)], _S2204 + _S2219[int(12)], _S2205 + _S2219[int(13)], _S2206 + _S2219[int(14)], _S2207 + _S2219[int(15)], _S2208 + _S2219[int(16)], _S2209 + _S2219[int(17)], _S2210 + _S2219[int(18)], _S2211 + _S2219[int(19)], _S2212 + _S2219[int(20)], _S2213 + _S2219[int(21)], _S2214 + _S2219[int(22)]
    };
    dpraw_losses_1->primal_0 = dpraw_losses_1->primal_0;
    dpraw_losses_1->differential_0 = _S2220;
    return;
}

inline __device__ void s_bwd_compute_ppisp_rqs_regularization_loss_0(DiffPair_arrayx3Cfloatx2C23x3E_0 * _S2221, int _S2222, FixedArray<float, 6>  * _S2223, FixedArray<float, 6>  * _S2224)
{
    s_bwd_prop_compute_ppisp_rqs_regularization_loss_0(_S2221, _S2222, _S2223, _S2224);
    return;
}

inline __device__ void compute_ppisp_rqs_regularization_loss_vjp(FixedArray<float, 23>  raw_losses_3, int num_cameras_5, FixedArray<float, 6>  loss_weights_5, FixedArray<float, 6>  grad_out_5, FixedArray<float, 23>  * _S2225)
{
    FixedArray<float, 23>  _S2226 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C23x3E_0 dp_raw_losses_1;
    (&dp_raw_losses_1)->primal_0 = raw_losses_3;
    (&dp_raw_losses_1)->differential_0 = _S2226;
    FixedArray<float, 6>  _S2227 = loss_weights_5;
    FixedArray<float, 6>  _S2228 = grad_out_5;
    s_bwd_compute_ppisp_rqs_regularization_loss_0(&dp_raw_losses_1, num_cameras_5, &_S2227, &_S2228);
    *_S2225 = (&dp_raw_losses_1)->differential_0;
    return;
}

