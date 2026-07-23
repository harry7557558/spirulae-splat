#pragma once

#include "generated/slang.cuh"

inline __device__ Matrix<float, 2, 2>  transpose_0(Matrix<float, 2, 2>  x_0)
{
    Matrix<float, 2, 2>  result_0;
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
            *_slang_vector_get_element_ptr(((&result_0)->rows + (r_0)), c_0) = _slang_vector_get_element(x_0.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_0;
}

inline __device__ Matrix<float, 3, 2>  transpose_1(Matrix<float, 2, 3>  x_1)
{
    Matrix<float, 3, 2>  result_1;
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
        int c_1 = int(0);
        for(;;)
        {
            if(c_1 < int(2))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_1)->rows + (r_1)), c_1) = _slang_vector_get_element(x_1.rows[c_1], r_1);
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_1;
}

inline __device__ Matrix<float, 2, 3>  transpose_2(Matrix<float, 3, 2>  x_2)
{
    Matrix<float, 2, 3>  result_2;
    int r_2 = int(0);
    for(;;)
    {
        if(r_2 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_2 = int(0);
        for(;;)
        {
            if(c_2 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_2)), c_2) = _slang_vector_get_element(x_2.rows[c_2], r_2);
            c_2 = c_2 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_2;
}

inline __device__ Matrix<float, 3, 3>  transpose_3(Matrix<float, 3, 3>  x_3)
{
    Matrix<float, 3, 3>  result_3;
    int r_3 = int(0);
    for(;;)
    {
        if(r_3 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_3 = int(0);
        for(;;)
        {
            if(c_3 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_3)->rows + (r_3)), c_3) = _slang_vector_get_element(x_3.rows[c_3], r_3);
            c_3 = c_3 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_3;
}

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_0)
{
    return m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_0)
{
    DiffPair_float_0 _S1 = *dpx_0;
    float _S2;
    if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
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
    if(((*dpy_0).primal_0) < (_S1.primal_0))
    {
        _S2 = dOut_0;
    }
    else
    {
        if(((*dpy_0).primal_0) > ((*dpx_0).primal_0))
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

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_1, float dOut_1)
{
    DiffPair_float_0 _S4 = *dpx_1;
    float _S5;
    if(((*dpx_1).primal_0) > ((*dpy_1).primal_0))
    {
        _S5 = dOut_1;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_1).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_1;
        }
    }
    dpx_1->primal_0 = _S4.primal_0;
    dpx_1->differential_0 = _S5;
    DiffPair_float_0 _S6 = *dpy_1;
    if(((*dpy_1).primal_0) > (_S4.primal_0))
    {
        _S5 = dOut_1;
    }
    else
    {
        if(((*dpy_1).primal_0) < ((*dpx_1).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_1;
        }
    }
    dpy_1->primal_0 = _S6.primal_0;
    dpy_1->differential_0 = _S5;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_2, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_2)
{
    DiffPair_float_0 _S7 = *dpx_2;
    bool _S8;
    if(((*dpx_2).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S8 = ((*dpx_2).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S8 = false;
    }
    float _S9;
    if(_S8)
    {
        _S9 = dOut_2;
    }
    else
    {
        _S9 = 0.0f;
    }
    dpx_2->primal_0 = _S7.primal_0;
    dpx_2->differential_0 = _S9;
    DiffPair_float_0 _S10 = *dpMin_0;
    if((_S7.primal_0) < ((*dpMin_0).primal_0))
    {
        _S9 = dOut_2;
    }
    else
    {
        _S9 = 0.0f;
    }
    dpMin_0->primal_0 = _S10.primal_0;
    dpMin_0->differential_0 = _S9;
    DiffPair_float_0 _S11 = *dpMax_0;
    if(((*dpx_2).primal_0) > ((*dpMax_0).primal_0))
    {
        _S9 = dOut_2;
    }
    else
    {
        _S9 = 0.0f;
    }
    dpMax_0->primal_0 = _S11.primal_0;
    dpMax_0->differential_0 = _S9;
    return;
}

inline __device__ float clamp_0(float x_4, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_4), (minBound_0)))), (maxBound_0)));
}

struct DiffPair_matrixx3Cfloatx2C2x2C3x3E_0
{
    Matrix<float, 2, 3>  primal_0;
    Matrix<float, 2, 3>  differential_0;
};

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

inline __device__ void mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_0, Matrix<float, 2, 3>  dOut_3)
{
    Matrix<float, 2, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_0;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = *&(((&left_d_result_0)->rows + (int(0)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = *&(((&left_d_result_0)->rows + (int(0)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_3.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = *&(((&left_d_result_0)->rows + (int(0)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_3.rows[int(0)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(0)))->x) = *&(((&right_d_result_0)->rows + (int(0)))->x) + (*left_0).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(1)))->x) = *&(((&right_d_result_0)->rows + (int(1)))->x) + (*left_0).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_0)->rows + (int(2)))->x) = *&(((&right_d_result_0)->rows + (int(2)))->x) + (*left_0).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(0)))->y) = *&(((&right_d_result_0)->rows + (int(0)))->y) + (*left_0).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(1)))->y) = *&(((&right_d_result_0)->rows + (int(1)))->y) + (*left_0).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_0)->rows + (int(2)))->y) = *&(((&right_d_result_0)->rows + (int(2)))->y) + (*left_0).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = *&(((&left_d_result_0)->rows + (int(1)))->x) + (*right_0).primal_0.rows[int(0)].z * dOut_3.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(0)))->z) = *&(((&right_d_result_0)->rows + (int(0)))->z) + (*left_0).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = *&(((&left_d_result_0)->rows + (int(1)))->y) + (*right_0).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(1)))->z) = *&(((&right_d_result_0)->rows + (int(1)))->z) + (*left_0).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].z;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = *&(((&left_d_result_0)->rows + (int(1)))->z) + (*right_0).primal_0.rows[int(2)].z * dOut_3.rows[int(1)].z;
    *&(((&right_d_result_0)->rows + (int(2)))->z) = *&(((&right_d_result_0)->rows + (int(2)))->z) + (*left_0).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].z;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C2x3E_0
{
    Matrix<float, 3, 2>  primal_0;
    Matrix<float, 3, 2>  differential_0;
};

inline __device__ void mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_1, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * right_1, Matrix<float, 2, 2>  dOut_4)
{
    Matrix<float, 2, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = 0.0f;
    Matrix<float, 3, 2>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_1).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_1).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_1).primal_0.rows[int(1)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_1).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_1).primal_0.rows[int(2)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_1).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_1).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_1).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_1).primal_0.rows[int(1)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_1).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_1).primal_0.rows[int(2)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_1).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_1).primal_0.rows[int(0)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_1).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_1).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_1).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_1).primal_0.rows[int(2)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_1).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_1).primal_0.rows[int(0)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_1).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_1).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_1).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_1).primal_0.rows[int(2)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_1).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].y;
    left_1->primal_0 = (*left_1).primal_0;
    left_1->differential_0 = left_d_result_1;
    right_1->primal_0 = (*right_1).primal_0;
    right_1->differential_0 = right_d_result_1;
    return;
}

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, Matrix<float, 3, 3>  dOut_5)
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
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].x;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].y;
    *&(((&left_d_result_2)->rows + (int(2)))->x) = *&(((&left_d_result_2)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(2)].x * dOut_5.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->y) = *&(((&left_d_result_2)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(2)].y * dOut_5.rows[int(2)].z;
    *&(((&left_d_result_2)->rows + (int(2)))->z) = *&(((&left_d_result_2)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(2)].z * dOut_5.rows[int(2)].z;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_2;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_2;
    return;
}

inline __device__ Matrix<float, 2, 3>  mul_3(Matrix<float, 2, 3>  left_3, Matrix<float, 3, 3>  right_3)
{
    Matrix<float, 2, 3>  result_4;
    int r_4 = int(0);
    for(;;)
    {
        if(r_4 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_4 = int(0);
        for(;;)
        {
            if(c_4 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_0 = int(0);
            float sum_0 = 0.0f;
            for(;;)
            {
                if(i_0 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_1 = sum_0 + _slang_vector_get_element(left_3.rows[r_4], i_0) * _slang_vector_get_element(right_3.rows[i_0], c_4);
                i_0 = i_0 + int(1);
                sum_0 = sum_1;
            }
            *_slang_vector_get_element_ptr(((&result_4)->rows + (r_4)), c_4) = sum_0;
            c_4 = c_4 + int(1);
        }
        r_4 = r_4 + int(1);
    }
    return result_4;
}

inline __device__ Matrix<float, 2, 2>  mul_4(Matrix<float, 2, 3>  left_4, Matrix<float, 3, 2>  right_4)
{
    Matrix<float, 2, 2>  result_5;
    int r_5 = int(0);
    for(;;)
    {
        if(r_5 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_5 = int(0);
        for(;;)
        {
            if(c_5 < int(2))
            {
            }
            else
            {
                break;
            }
            int i_1 = int(0);
            float sum_2 = 0.0f;
            for(;;)
            {
                if(i_1 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_3 = sum_2 + _slang_vector_get_element(left_4.rows[r_5], i_1) * _slang_vector_get_element(right_4.rows[i_1], c_5);
                i_1 = i_1 + int(1);
                sum_2 = sum_3;
            }
            *_slang_vector_get_element_ptr(((&result_5)->rows + (r_5)), c_5) = sum_2;
            c_5 = c_5 + int(1);
        }
        r_5 = r_5 + int(1);
    }
    return result_5;
}

inline __device__ Matrix<float, 3, 3>  mul_5(Matrix<float, 3, 3>  left_5, Matrix<float, 3, 3>  right_5)
{
    Matrix<float, 3, 3>  result_6;
    int r_6 = int(0);
    for(;;)
    {
        if(r_6 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_6 = int(0);
        for(;;)
        {
            if(c_6 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_2 = int(0);
            float sum_4 = 0.0f;
            for(;;)
            {
                if(i_2 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_5 = sum_4 + _slang_vector_get_element(left_5.rows[r_6], i_2) * _slang_vector_get_element(right_5.rows[i_2], c_6);
                i_2 = i_2 + int(1);
                sum_4 = sum_5;
            }
            *_slang_vector_get_element_ptr(((&result_6)->rows + (r_6)), c_6) = sum_4;
            c_6 = c_6 + int(1);
        }
        r_6 = r_6 + int(1);
    }
    return result_6;
}

inline __device__ bool persp_proj_3dgs_nav(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float4  intrins_0, FixedArray<float, 10>  dist_coeffs_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    bool _S12;
    float2  _S13;
    float _S14;
    for(;;)
    {
        float cx_0 = intrins_0.z;
        float cy_0 = intrins_0.w;
        for(;;)
        {
            float2  _S15 = float2 {mean3d_0.x, mean3d_0.y};
            _S13 = _S15;
            float _S16 = mean3d_0.z;
            _S14 = _S16;
            *mean2d_0 = _S15 / make_float2 (_S16);
            if(_S16 < 0.0f)
            {
                _S12 = true;
            }
            else
            {
                float u_0 = (*mean2d_0).x;
                float v_0 = (*mean2d_0).y;
                float _S17 = u_0 + u_0;
                float r2_0 = u_0 * u_0 + v_0 * v_0;
                float _S18 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                float _S19 = dist_coeffs_0[int(1)] + r2_0 * _S18;
                float _S20 = dist_coeffs_0[int(0)] + r2_0 * _S19;
                float radial_0 = 1.0f + r2_0 * _S20;
                float _S21 = 2.0f * dist_coeffs_0[int(4)];
                float _S22 = 2.0f * u_0;
                float _S23 = 2.0f * dist_coeffs_0[int(5)];
                float _S24 = 2.0f * v_0;
                float2  _S25 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (_S17 * _S20 + (_S17 * _S19 + (_S17 * _S18 + _S17 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * *mean2d_0 + make_float2 (_S21 * v_0 + (_S17 + (_S22 + _S22)) * dist_coeffs_0[int(5)] + _S17 * dist_coeffs_0[int(6)], _S23 * v_0 + _S17 * dist_coeffs_0[int(4)] + _S17 * dist_coeffs_0[int(7)]);
                float _S26 = v_0 + v_0;
                float2  _S27 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (_S26 * _S20 + (_S26 * _S19 + (_S26 * _S18 + _S26 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * *mean2d_0 + make_float2 (_S21 * u_0 + _S26 * dist_coeffs_0[int(5)] + _S26 * dist_coeffs_0[int(6)], _S23 * u_0 + (_S26 + (_S24 + _S24)) * dist_coeffs_0[int(4)] + _S26 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S28 = transpose_0(makeMatrix<float, 2, 2> (_S25 + make_float2 (_S25.x * dist_coeffs_0[int(8)] + _S25.y * dist_coeffs_0[int(9)], 0.0f), _S27 + make_float2 (_S27.x * dist_coeffs_0[int(8)] + _S27.y * dist_coeffs_0[int(9)], 0.0f)));
                _S12 = !((F32_min((determinant_0(_S28)), ((F32_min((_S28.rows[int(0)].x), (_S28.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S12)
            {
                break;
            }
            float u_1 = (*mean2d_0).x;
            float v_1 = (*mean2d_0).y;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float2  _S29 = *mean2d_0 * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
            float2  _S30 = _S29 + make_float2 (dist_coeffs_0[int(8)] * _S29.x + dist_coeffs_0[int(9)] * _S29.y, 0.0f);
            *mean2d_0 = make_float2 (intrins_0.x * _S30.x + cx_0, intrins_0.y * _S30.y + cy_0);
            break;
        }
        if(!!_S12)
        {
            _S12 = false;
            break;
        }
        Matrix<float, 2, 3>  J_0;
        float2  _S31 = _S13 / make_float2 (_S14);
        float _S32 = _S14 * _S14;
        float2  _S33 = make_float2 (1.0f, 0.0f) * make_float2 (_S14) / make_float2 (_S32);
        float u_2 = _S31.x;
        float s_diff_u_0 = _S33.x;
        float v_2 = _S31.y;
        float s_diff_v_0 = _S33.y;
        float _S34 = s_diff_u_0 * u_2;
        float _S35 = s_diff_v_0 * v_2;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float s_diff_r2_0 = _S34 + _S34 + (_S35 + _S35);
        float _S36 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
        float _S37 = dist_coeffs_0[int(1)] + r2_2 * _S36;
        float _S38 = dist_coeffs_0[int(0)] + r2_2 * _S37;
        float _S39 = 2.0f * dist_coeffs_0[int(4)];
        float _S40 = 2.0f * dist_coeffs_0[int(5)];
        float2  _S41 = _S33 * make_float2 (1.0f + r2_2 * _S38) + make_float2 (s_diff_r2_0 * _S38 + (s_diff_r2_0 * _S37 + (s_diff_r2_0 * _S36 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * _S31 + make_float2 (s_diff_u_0 * _S39 * v_2 + s_diff_v_0 * (_S39 * u_2) + (s_diff_r2_0 + (s_diff_u_0 * 2.0f * u_2 + s_diff_u_0 * (2.0f * u_2))) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], s_diff_u_0 * _S40 * v_2 + s_diff_v_0 * (_S40 * u_2) + (s_diff_r2_0 + (s_diff_v_0 * 2.0f * v_2 + s_diff_v_0 * (2.0f * v_2))) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
        float2  _S42 = _S41 + make_float2 (_S41.x * dist_coeffs_0[int(8)] + _S41.y * dist_coeffs_0[int(9)], 0.0f);
        float fx_0 = intrins_0.x;
        float fy_0 = intrins_0.y;
        float _S43 = _S42.y * fy_0;
        Matrix<float, 2, 3>  J_1;
        *&(((&J_1)->rows + (int(0)))->x) = _S42.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->x) = _S43;
        float2  _S44 = _S13 / make_float2 (_S14);
        float2  _S45 = make_float2 (0.0f, 1.0f) * make_float2 (_S14) / make_float2 (_S32);
        float u_3 = _S44.x;
        float s_diff_u_1 = _S45.x;
        float v_3 = _S44.y;
        float s_diff_v_1 = _S45.y;
        float _S46 = s_diff_u_1 * u_3;
        float _S47 = s_diff_v_1 * v_3;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float s_diff_r2_1 = _S46 + _S46 + (_S47 + _S47);
        float _S48 = dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)];
        float _S49 = dist_coeffs_0[int(1)] + r2_3 * _S48;
        float _S50 = dist_coeffs_0[int(0)] + r2_3 * _S49;
        float2  _S51 = _S45 * make_float2 (1.0f + r2_3 * _S50) + make_float2 (s_diff_r2_1 * _S50 + (s_diff_r2_1 * _S49 + (s_diff_r2_1 * _S48 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_3) * r2_3) * r2_3) * _S44 + make_float2 (s_diff_u_1 * _S39 * v_3 + s_diff_v_1 * (_S39 * u_3) + (s_diff_r2_1 + (s_diff_u_1 * 2.0f * u_3 + s_diff_u_1 * (2.0f * u_3))) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], s_diff_u_1 * _S40 * v_3 + s_diff_v_1 * (_S40 * u_3) + (s_diff_r2_1 + (s_diff_v_1 * 2.0f * v_3 + s_diff_v_1 * (2.0f * v_3))) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
        float2  _S52 = _S51 + make_float2 (_S51.x * dist_coeffs_0[int(8)] + _S51.y * dist_coeffs_0[int(9)], 0.0f);
        float _S53 = _S52.y * fy_0;
        *&(((&J_1)->rows + (int(0)))->y) = _S52.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->y) = _S53;
        float2  _S54 = _S13 / make_float2 (_S14);
        float2  _S55 = (make_float2 (0.0f, 0.0f) - _S13) / make_float2 (_S32);
        float u_4 = _S54.x;
        float s_diff_u_2 = _S55.x;
        float v_4 = _S54.y;
        float s_diff_v_2 = _S55.y;
        float _S56 = s_diff_u_2 * u_4;
        float _S57 = s_diff_v_2 * v_4;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float s_diff_r2_2 = _S56 + _S56 + (_S57 + _S57);
        float _S58 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
        float _S59 = dist_coeffs_0[int(1)] + r2_4 * _S58;
        float _S60 = dist_coeffs_0[int(0)] + r2_4 * _S59;
        float2  _S61 = _S55 * make_float2 (1.0f + r2_4 * _S60) + make_float2 (s_diff_r2_2 * _S60 + (s_diff_r2_2 * _S59 + (s_diff_r2_2 * _S58 + s_diff_r2_2 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * _S54 + make_float2 (s_diff_u_2 * _S39 * v_4 + s_diff_v_2 * (_S39 * u_4) + (s_diff_r2_2 + (s_diff_u_2 * 2.0f * u_4 + s_diff_u_2 * (2.0f * u_4))) * dist_coeffs_0[int(5)] + s_diff_r2_2 * dist_coeffs_0[int(6)], s_diff_u_2 * _S40 * v_4 + s_diff_v_2 * (_S40 * u_4) + (s_diff_r2_2 + (s_diff_v_2 * 2.0f * v_4 + s_diff_v_2 * (2.0f * v_4))) * dist_coeffs_0[int(4)] + s_diff_r2_2 * dist_coeffs_0[int(7)]);
        float2  _S62 = _S61 + make_float2 (_S61.x * dist_coeffs_0[int(8)] + _S61.y * dist_coeffs_0[int(9)], 0.0f);
        float _S63 = _S62.y * fy_0;
        *&(((&J_1)->rows + (int(0)))->z) = _S62.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->z) = _S63;
        J_0 = J_1;
        float _S64 = float(width_0);
        float _S65 = 0.30000001192092896f * (0.5f * _S64);
        float lim_x_pos_0 = _S64 + _S65;
        float rz_0 = 1.0f / _S14;
        float _S66 = - _S65;
        float max_Jyz_0 = - (_S66 - cy_0) * rz_0;
        float min_Jyz_0 = - (lim_x_pos_0 - cy_0) * rz_0;
        *&(((&J_0)->rows + (int(0)))->z) = clamp_0(*&(((&J_0)->rows + (int(0)))->z), - (lim_x_pos_0 - cx_0) * rz_0, - (_S66 - cx_0) * rz_0);
        *&(((&J_0)->rows + (int(1)))->z) = clamp_0(*&(((&J_0)->rows + (int(1)))->z), min_Jyz_0, max_Jyz_0);
        *cov2d_0 = mul_4(mul_3(J_0, cov3d_0), transpose_1(J_0));
        _S12 = true;
        break;
    }
    return _S12;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_3, float dOut_6)
{
    float _S67 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_3).primal_0)))))) * dOut_6;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S67;
    return;
}

inline __device__ DiffPair_float_0 _d_sqrt_1(DiffPair_float_0 * dpx_4)
{
    DiffPair_float_0 _S68 = { (F32_sqrt((dpx_4->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_4->primal_0)))))) * dpx_4->differential_0 };
    return _S68;
}

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, float dOut_7)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_7;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_5).primal_0.x * dOut_7;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_7;
    *&((&y_d_result_0)->y) = (*dpx_5).primal_0.y * dOut_7;
    *&((&x_d_result_0)->z) = (*dpy_2).primal_0.z * dOut_7;
    *&((&y_d_result_0)->z) = (*dpx_5).primal_0.z * dOut_7;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float2  x_5, float2  y_0)
{
    int i_3 = int(0);
    float result_7 = 0.0f;
    for(;;)
    {
        if(i_3 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_8 = result_7 + _slang_vector_get_element(x_5, i_3) * _slang_vector_get_element(y_0, i_3);
        i_3 = i_3 + int(1);
        result_7 = result_8;
    }
    return result_7;
}

inline __device__ float dot_1(float3  x_6, float3  y_1)
{
    int i_4 = int(0);
    float result_9 = 0.0f;
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_10 = result_9 + _slang_vector_get_element(x_6, i_4) * _slang_vector_get_element(y_1, i_4);
        i_4 = i_4 + int(1);
        result_9 = result_10;
    }
    return result_9;
}

inline __device__ float dot_2(float4  x_7, float4  y_2)
{
    int i_5 = int(0);
    float result_11 = 0.0f;
    for(;;)
    {
        if(i_5 < int(4))
        {
        }
        else
        {
            break;
        }
        float result_12 = result_11 + _slang_vector_get_element(x_7, i_5) * _slang_vector_get_element(y_2, i_5);
        i_5 = i_5 + int(1);
        result_11 = result_12;
    }
    return result_11;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_0(x_8, x_8))));
}

inline __device__ float length_1(float3  x_9)
{
    return (F32_sqrt((dot_1(x_9, x_9))));
}

inline __device__ float length_2(float4  x_10)
{
    return (F32_sqrt((dot_2(x_10, x_10))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_3, DiffPair_float_0 * dpx_6, float dOut_8)
{
    DiffPair_float_0 _S69 = *dpx_6;
    float _S70 = - (*dpy_3).primal_0 / ((*dpx_6).primal_0 * (*dpx_6).primal_0 + (*dpy_3).primal_0 * (*dpy_3).primal_0) * dOut_8;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S70;
    float _S71 = _S69.primal_0 / (_S69.primal_0 * _S69.primal_0 + (*dpy_3).primal_0 * (*dpy_3).primal_0) * dOut_8;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S71;
    return;
}

inline __device__ DiffPair_float_0 _d_atan2_1(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_7)
{
    float _S72 = dpx_7->primal_0 * dpx_7->primal_0 + dpy_4->primal_0 * dpy_4->primal_0;
    DiffPair_float_0 _S73 = { (F32_atan2((dpy_4->primal_0), (dpx_7->primal_0))), - dpy_4->primal_0 / _S72 * dpx_7->differential_0 + dpx_7->primal_0 / _S72 * dpy_4->differential_0 };
    return _S73;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8)
{
    float _S74 = *&((&dpx_8->differential_0)->x) * *&((&dpx_8->primal_0)->x);
    float _S75 = *&((&dpx_8->differential_0)->y) * *&((&dpx_8->primal_0)->y);
    float s_diff_len_0 = _S74 + _S74 + (_S75 + _S75);
    DiffPair_float_0 _S76;
    (&_S76)->primal_0 = *&((&dpx_8->primal_0)->x) * *&((&dpx_8->primal_0)->x) + *&((&dpx_8->primal_0)->y) * *&((&dpx_8->primal_0)->y);
    (&_S76)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S77 = _d_sqrt_1(&_S76);
    DiffPair_float_0 _S78 = { _S77.primal_0, _S77.differential_0 };
    return _S78;
}

inline __device__ bool fisheye_proj_3dgs_nav(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    bool _S79;
    float2  _S80;
    float _S81;
    float _S82;
    float _S83;
    float _S84;
    float _S85;
    float _S86;
    float _S87;
    float _S88;
    float _S89;
    float _S90;
    float _S91;
    float _S92;
    float _S93;
    bool _S94;
    for(;;)
    {
        float k_0;
        for(;;)
        {
            float2  _S95 = float2 {mean3d_1.x, mean3d_1.y};
            _S80 = _S95;
            float r_7 = length_0(_S95);
            float _S96 = mean3d_1.z;
            _S81 = _S96;
            float theta_0 = (F32_atan2((r_7), (_S96)));
            if(theta_0 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S96;
            }
            else
            {
                k_0 = theta_0 / r_7;
            }
            float2  _S97 = _S95 * make_float2 (k_0);
            *mean2d_1 = _S97;
            float2  _S98 = make_float2 (1.0f, 0.0f);
            _S82 = dist_coeffs_1[int(0)];
            _S83 = dist_coeffs_1[int(1)];
            _S84 = dist_coeffs_1[int(2)];
            _S85 = dist_coeffs_1[int(3)];
            _S86 = dist_coeffs_1[int(4)];
            _S87 = dist_coeffs_1[int(5)];
            _S88 = dist_coeffs_1[int(6)];
            _S89 = dist_coeffs_1[int(7)];
            _S90 = dist_coeffs_1[int(8)];
            _S91 = dist_coeffs_1[int(9)];
            float u_5 = _S97.x;
            float v_5 = _S97.y;
            float _S99 = u_5 + u_5;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float _S100 = dist_coeffs_1[int(2)] + r2_5 * dist_coeffs_1[int(3)];
            float _S101 = dist_coeffs_1[int(1)] + r2_5 * _S100;
            float _S102 = dist_coeffs_1[int(0)] + r2_5 * _S101;
            float _S103 = _S99 * _S102 + (_S99 * _S101 + (_S99 * _S100 + _S99 * dist_coeffs_1[int(3)] * r2_5) * r2_5) * r2_5;
            float radial_1 = 1.0f + r2_5 * _S102;
            float _S104 = 2.0f * dist_coeffs_1[int(4)];
            _S92 = _S104;
            float _S105 = _S104 * u_5;
            float _S106 = 2.0f * u_5;
            float s_diff_du_0 = _S104 * v_5 + (_S99 + (_S106 + _S106)) * dist_coeffs_1[int(5)] + _S99 * dist_coeffs_1[int(6)];
            float _S107 = 2.0f * dist_coeffs_1[int(5)];
            _S93 = _S107;
            float _S108 = 2.0f * v_5;
            float2  _S109 = _S98 * make_float2 (radial_1) + make_float2 (_S103) * _S97 + make_float2 (s_diff_du_0, _S107 * v_5 + _S99 * dist_coeffs_1[int(4)] + _S99 * dist_coeffs_1[int(7)]);
            float _S110 = v_5 + v_5;
            float2  _S111 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (_S110 * _S102 + (_S110 * _S101 + (_S110 * _S100 + _S110 * dist_coeffs_1[int(3)] * r2_5) * r2_5) * r2_5) * _S97 + make_float2 (_S105 + _S110 * dist_coeffs_1[int(5)] + _S110 * dist_coeffs_1[int(6)], _S107 * u_5 + (_S110 + (_S108 + _S108)) * dist_coeffs_1[int(4)] + _S110 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S112 = transpose_0(makeMatrix<float, 2, 2> (_S109 + make_float2 (_S109.x * dist_coeffs_1[int(8)] + _S109.y * dist_coeffs_1[int(9)], 0.0f), _S111 + make_float2 (_S111.x * dist_coeffs_1[int(8)] + _S111.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S113 = !((F32_min((determinant_0(_S112)), ((F32_min((_S112.rows[int(0)].x), (_S112.rows[int(1)].y)))))) > 0.0f);
            _S94 = _S113;
            if(_S113)
            {
                break;
            }
            float u_6 = (*mean2d_1).x;
            float v_6 = (*mean2d_1).y;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float2  _S114 = *mean2d_1 * make_float2 (1.0f + r2_6 * (dist_coeffs_1[int(0)] + r2_6 * (dist_coeffs_1[int(1)] + r2_6 * (dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)])))) + make_float2 (_S104 * u_6 * v_6 + dist_coeffs_1[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + dist_coeffs_1[int(6)] * r2_6, _S107 * u_6 * v_6 + dist_coeffs_1[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + dist_coeffs_1[int(7)] * r2_6);
            float2  _S115 = _S114 + make_float2 (dist_coeffs_1[int(8)] * _S114.x + dist_coeffs_1[int(9)] * _S114.y, 0.0f);
            *mean2d_1 = make_float2 (intrins_1.x * _S115.x + intrins_1.z, intrins_1.y * _S115.y + intrins_1.w);
            break;
        }
        if(!!_S94)
        {
            _S79 = false;
            break;
        }
        Matrix<float, 2, 3>  J_2;
        float2  _S116 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S117;
        (&_S117)->primal_0 = _S80;
        (&_S117)->differential_0 = _S116;
        DiffPair_float_0 _S118 = s_fwd_length_impl_0(&_S117);
        float _S119 = _S81;
        DiffPair_float_0 _S120;
        (&_S120)->primal_0 = _S118.primal_0;
        (&_S120)->differential_0 = _S118.differential_0;
        DiffPair_float_0 _S121;
        (&_S121)->primal_0 = _S81;
        (&_S121)->differential_0 = 0.0f;
        DiffPair_float_0 _S122 = _d_atan2_1(&_S120, &_S121);
        float s_diff_k_0;
        if((_S122.primal_0) < 0.00100000004749745f)
        {
            float _S123 = _S122.differential_0 * _S122.primal_0;
            float _S124 = (0.0f - (_S123 + _S123) * 0.3333333432674408f) * _S81 / (_S81 * _S81);
            k_0 = (1.0f - _S122.primal_0 * _S122.primal_0 / 3.0f) / _S81;
            s_diff_k_0 = _S124;
        }
        else
        {
            float _S125 = (_S122.differential_0 * _S118.primal_0 - _S122.primal_0 * _S118.differential_0) / (_S118.primal_0 * _S118.primal_0);
            k_0 = _S122.primal_0 / _S118.primal_0;
            s_diff_k_0 = _S125;
        }
        float2  _S126 = _S80 * make_float2 (k_0);
        float2  _S127 = _S116 * make_float2 (k_0) + make_float2 (s_diff_k_0) * _S80;
        float u_7 = _S126.x;
        float s_diff_u_3 = _S127.x;
        float v_7 = _S126.y;
        float s_diff_v_3 = _S127.y;
        float _S128 = s_diff_u_3 * u_7;
        float _S129 = s_diff_v_3 * v_7;
        float r2_7 = u_7 * u_7 + v_7 * v_7;
        float s_diff_r2_3 = _S128 + _S128 + (_S129 + _S129);
        float _S130 = _S84 + r2_7 * _S85;
        float _S131 = _S83 + r2_7 * _S130;
        float _S132 = _S82 + r2_7 * _S131;
        float2  _S133 = _S127 * make_float2 (1.0f + r2_7 * _S132) + make_float2 (s_diff_r2_3 * _S132 + (s_diff_r2_3 * _S131 + (s_diff_r2_3 * _S130 + s_diff_r2_3 * _S85 * r2_7) * r2_7) * r2_7) * _S126 + make_float2 (s_diff_u_3 * _S92 * v_7 + s_diff_v_3 * (_S92 * u_7) + (s_diff_r2_3 + (s_diff_u_3 * 2.0f * u_7 + s_diff_u_3 * (2.0f * u_7))) * _S87 + s_diff_r2_3 * _S88, s_diff_u_3 * _S93 * v_7 + s_diff_v_3 * (_S93 * u_7) + (s_diff_r2_3 + (s_diff_v_3 * 2.0f * v_7 + s_diff_v_3 * (2.0f * v_7))) * _S86 + s_diff_r2_3 * _S89);
        float2  _S134 = _S133 + make_float2 (_S133.x * _S90 + _S133.y * _S91, 0.0f);
        float fx_1 = intrins_1.x;
        float fy_1 = intrins_1.y;
        float _S135 = _S134.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->x) = _S134.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->x) = _S135;
        float2  _S136 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S137;
        (&_S137)->primal_0 = _S80;
        (&_S137)->differential_0 = _S136;
        DiffPair_float_0 _S138 = s_fwd_length_impl_0(&_S137);
        DiffPair_float_0 _S139;
        (&_S139)->primal_0 = _S138.primal_0;
        (&_S139)->differential_0 = _S138.differential_0;
        DiffPair_float_0 _S140;
        (&_S140)->primal_0 = _S119;
        (&_S140)->differential_0 = 0.0f;
        DiffPair_float_0 _S141 = _d_atan2_1(&_S139, &_S140);
        if((_S141.primal_0) < 0.00100000004749745f)
        {
            float _S142 = _S141.differential_0 * _S141.primal_0;
            float _S143 = (0.0f - (_S142 + _S142) * 0.3333333432674408f) * _S81 / (_S81 * _S81);
            k_0 = (1.0f - _S141.primal_0 * _S141.primal_0 / 3.0f) / _S81;
            s_diff_k_0 = _S143;
        }
        else
        {
            float _S144 = (_S141.differential_0 * _S138.primal_0 - _S141.primal_0 * _S138.differential_0) / (_S138.primal_0 * _S138.primal_0);
            k_0 = _S141.primal_0 / _S138.primal_0;
            s_diff_k_0 = _S144;
        }
        float2  _S145 = _S80 * make_float2 (k_0);
        float2  _S146 = _S136 * make_float2 (k_0) + make_float2 (s_diff_k_0) * _S80;
        float u_8 = _S145.x;
        float s_diff_u_4 = _S146.x;
        float v_8 = _S145.y;
        float s_diff_v_4 = _S146.y;
        float _S147 = s_diff_u_4 * u_8;
        float _S148 = s_diff_v_4 * v_8;
        float r2_8 = u_8 * u_8 + v_8 * v_8;
        float s_diff_r2_4 = _S147 + _S147 + (_S148 + _S148);
        float _S149 = _S84 + r2_8 * _S85;
        float _S150 = _S83 + r2_8 * _S149;
        float _S151 = _S82 + r2_8 * _S150;
        float2  _S152 = _S146 * make_float2 (1.0f + r2_8 * _S151) + make_float2 (s_diff_r2_4 * _S151 + (s_diff_r2_4 * _S150 + (s_diff_r2_4 * _S149 + s_diff_r2_4 * _S85 * r2_8) * r2_8) * r2_8) * _S145 + make_float2 (s_diff_u_4 * _S92 * v_8 + s_diff_v_4 * (_S92 * u_8) + (s_diff_r2_4 + (s_diff_u_4 * 2.0f * u_8 + s_diff_u_4 * (2.0f * u_8))) * _S87 + s_diff_r2_4 * _S88, s_diff_u_4 * _S93 * v_8 + s_diff_v_4 * (_S93 * u_8) + (s_diff_r2_4 + (s_diff_v_4 * 2.0f * v_8 + s_diff_v_4 * (2.0f * v_8))) * _S86 + s_diff_r2_4 * _S89);
        float2  _S153 = _S152 + make_float2 (_S152.x * _S90 + _S152.y * _S91, 0.0f);
        float _S154 = _S153.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->y) = _S153.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->y) = _S154;
        float2  _S155 = make_float2 (0.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S156;
        (&_S156)->primal_0 = _S80;
        (&_S156)->differential_0 = _S155;
        DiffPair_float_0 _S157 = s_fwd_length_impl_0(&_S156);
        DiffPair_float_0 _S158;
        (&_S158)->primal_0 = _S157.primal_0;
        (&_S158)->differential_0 = _S157.differential_0;
        DiffPair_float_0 _S159;
        (&_S159)->primal_0 = _S81;
        (&_S159)->differential_0 = 1.0f;
        DiffPair_float_0 _S160 = _d_atan2_1(&_S158, &_S159);
        if((_S160.primal_0) < 0.00100000004749745f)
        {
            float _S161 = _S160.differential_0 * _S160.primal_0;
            float _S162 = 1.0f - _S160.primal_0 * _S160.primal_0 / 3.0f;
            float _S163 = ((0.0f - (_S161 + _S161) * 0.3333333432674408f) * _S81 - _S162) / (_S81 * _S81);
            k_0 = _S162 / _S81;
            s_diff_k_0 = _S163;
        }
        else
        {
            float _S164 = (_S160.differential_0 * _S157.primal_0 - _S160.primal_0 * _S157.differential_0) / (_S157.primal_0 * _S157.primal_0);
            k_0 = _S160.primal_0 / _S157.primal_0;
            s_diff_k_0 = _S164;
        }
        float2  _S165 = _S80 * make_float2 (k_0);
        float2  _S166 = make_float2 (s_diff_k_0) * _S80;
        float u_9 = _S165.x;
        float s_diff_u_5 = _S166.x;
        float v_9 = _S165.y;
        float s_diff_v_5 = _S166.y;
        float _S167 = s_diff_u_5 * u_9;
        float _S168 = s_diff_v_5 * v_9;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float s_diff_r2_5 = _S167 + _S167 + (_S168 + _S168);
        float _S169 = _S84 + r2_9 * _S85;
        float _S170 = _S83 + r2_9 * _S169;
        float _S171 = _S82 + r2_9 * _S170;
        float2  _S172 = _S166 * make_float2 (1.0f + r2_9 * _S171) + make_float2 (s_diff_r2_5 * _S171 + (s_diff_r2_5 * _S170 + (s_diff_r2_5 * _S169 + s_diff_r2_5 * _S85 * r2_9) * r2_9) * r2_9) * _S165 + make_float2 (s_diff_u_5 * _S92 * v_9 + s_diff_v_5 * (_S92 * u_9) + (s_diff_r2_5 + (s_diff_u_5 * 2.0f * u_9 + s_diff_u_5 * (2.0f * u_9))) * _S87 + s_diff_r2_5 * _S88, s_diff_u_5 * _S93 * v_9 + s_diff_v_5 * (_S93 * u_9) + (s_diff_r2_5 + (s_diff_v_5 * 2.0f * v_9 + s_diff_v_5 * (2.0f * v_9))) * _S86 + s_diff_r2_5 * _S89);
        float2  _S173 = _S172 + make_float2 (_S172.x * _S90 + _S172.y * _S91, 0.0f);
        float _S174 = _S173.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->z) = _S173.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->z) = _S174;
        *cov2d_1 = mul_4(mul_3(J_2, cov3d_1), transpose_1(J_2));
        _S79 = true;
        break;
    }
    return _S79;
}

inline __device__ void _d_cos_0(DiffPair_float_0 * dpx_9, float dOut_9)
{
    float _S175 = - (F32_sin(((*dpx_9).primal_0))) * dOut_9;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S175;
    return;
}

inline __device__ void _d_sin_0(DiffPair_float_0 * dpx_10, float dOut_10)
{
    float _S176 = (F32_cos(((*dpx_10).primal_0))) * dOut_10;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S176;
    return;
}

inline __device__ DiffPair_float_0 _d_sin_1(DiffPair_float_0 * dpx_11)
{
    DiffPair_float_0 _S177 = { (F32_sin((dpx_11->primal_0))), (F32_cos((dpx_11->primal_0))) * dpx_11->differential_0 };
    return _S177;
}

inline __device__ bool equisolid_proj_3dgs_nav(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_2, FixedArray<float, 10>  dist_coeffs_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    bool _S178;
    float2  _S179;
    float _S180;
    float _S181;
    float _S182;
    float _S183;
    float _S184;
    float _S185;
    float _S186;
    float _S187;
    float _S188;
    float _S189;
    float _S190;
    float _S191;
    float _S192;
    bool _S193;
    for(;;)
    {
        float k_1;
        for(;;)
        {
            float2  _S194 = float2 {mean3d_2.x, mean3d_2.y};
            _S179 = _S194;
            float r_8 = length_0(_S194);
            float _S195 = mean3d_2.z;
            _S180 = _S195;
            float theta_1 = (F32_atan2((r_8), (_S195)));
            if(r_8 < 9.99999997475242708e-07f)
            {
                k_1 = (1.0f - theta_1 * theta_1 / 24.0f) / _S195;
            }
            else
            {
                k_1 = 2.0f * (F32_sin((0.5f * theta_1))) / r_8;
            }
            float2  _S196 = _S194 * make_float2 (k_1);
            *mean2d_2 = _S196;
            float2  _S197 = make_float2 (1.0f, 0.0f);
            _S181 = dist_coeffs_2[int(0)];
            _S182 = dist_coeffs_2[int(1)];
            _S183 = dist_coeffs_2[int(2)];
            _S184 = dist_coeffs_2[int(3)];
            _S185 = dist_coeffs_2[int(4)];
            _S186 = dist_coeffs_2[int(5)];
            _S187 = dist_coeffs_2[int(6)];
            _S188 = dist_coeffs_2[int(7)];
            _S189 = dist_coeffs_2[int(8)];
            _S190 = dist_coeffs_2[int(9)];
            float u_10 = _S196.x;
            float v_10 = _S196.y;
            float _S198 = u_10 + u_10;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float _S199 = dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)];
            float _S200 = dist_coeffs_2[int(1)] + r2_10 * _S199;
            float _S201 = dist_coeffs_2[int(0)] + r2_10 * _S200;
            float _S202 = _S198 * _S201 + (_S198 * _S200 + (_S198 * _S199 + _S198 * dist_coeffs_2[int(3)] * r2_10) * r2_10) * r2_10;
            float radial_2 = 1.0f + r2_10 * _S201;
            float _S203 = 2.0f * dist_coeffs_2[int(4)];
            _S191 = _S203;
            float _S204 = _S203 * u_10;
            float _S205 = 2.0f * u_10;
            float s_diff_du_1 = _S203 * v_10 + (_S198 + (_S205 + _S205)) * dist_coeffs_2[int(5)] + _S198 * dist_coeffs_2[int(6)];
            float _S206 = 2.0f * dist_coeffs_2[int(5)];
            _S192 = _S206;
            float _S207 = 2.0f * v_10;
            float2  _S208 = _S197 * make_float2 (radial_2) + make_float2 (_S202) * _S196 + make_float2 (s_diff_du_1, _S206 * v_10 + _S198 * dist_coeffs_2[int(4)] + _S198 * dist_coeffs_2[int(7)]);
            float _S209 = v_10 + v_10;
            float2  _S210 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (_S209 * _S201 + (_S209 * _S200 + (_S209 * _S199 + _S209 * dist_coeffs_2[int(3)] * r2_10) * r2_10) * r2_10) * _S196 + make_float2 (_S204 + _S209 * dist_coeffs_2[int(5)] + _S209 * dist_coeffs_2[int(6)], _S206 * u_10 + (_S209 + (_S207 + _S207)) * dist_coeffs_2[int(4)] + _S209 * dist_coeffs_2[int(7)]);
            Matrix<float, 2, 2>  _S211 = transpose_0(makeMatrix<float, 2, 2> (_S208 + make_float2 (_S208.x * dist_coeffs_2[int(8)] + _S208.y * dist_coeffs_2[int(9)], 0.0f), _S210 + make_float2 (_S210.x * dist_coeffs_2[int(8)] + _S210.y * dist_coeffs_2[int(9)], 0.0f)));
            bool _S212 = !((F32_min((determinant_0(_S211)), ((F32_min((_S211.rows[int(0)].x), (_S211.rows[int(1)].y)))))) > 0.0f);
            _S193 = _S212;
            if(_S212)
            {
                break;
            }
            float u_11 = (*mean2d_2).x;
            float v_11 = (*mean2d_2).y;
            float r2_11 = u_11 * u_11 + v_11 * v_11;
            float2  _S213 = *mean2d_2 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S203 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S206 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
            float2  _S214 = _S213 + make_float2 (dist_coeffs_2[int(8)] * _S213.x + dist_coeffs_2[int(9)] * _S213.y, 0.0f);
            *mean2d_2 = make_float2 (intrins_2.x * _S214.x + intrins_2.z, intrins_2.y * _S214.y + intrins_2.w);
            break;
        }
        if(!!_S193)
        {
            _S178 = false;
            break;
        }
        Matrix<float, 2, 3>  J_3;
        float2  _S215 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S216;
        (&_S216)->primal_0 = _S179;
        (&_S216)->differential_0 = _S215;
        DiffPair_float_0 _S217 = s_fwd_length_impl_0(&_S216);
        float _S218 = _S180;
        DiffPair_float_0 _S219;
        (&_S219)->primal_0 = _S217.primal_0;
        (&_S219)->differential_0 = _S217.differential_0;
        DiffPair_float_0 _S220;
        (&_S220)->primal_0 = _S180;
        (&_S220)->differential_0 = 0.0f;
        DiffPair_float_0 _S221 = _d_atan2_1(&_S219, &_S220);
        float s_diff_k_1;
        if((_S217.primal_0) < 9.99999997475242708e-07f)
        {
            float _S222 = _S221.differential_0 * _S221.primal_0;
            float _S223 = (0.0f - (_S222 + _S222) * 0.0416666679084301f) * _S180 / (_S180 * _S180);
            k_1 = (1.0f - _S221.primal_0 * _S221.primal_0 / 24.0f) / _S180;
            s_diff_k_1 = _S223;
        }
        else
        {
            float _S224 = _S221.differential_0 * 0.5f;
            DiffPair_float_0 _S225;
            (&_S225)->primal_0 = 0.5f * _S221.primal_0;
            (&_S225)->differential_0 = _S224;
            DiffPair_float_0 _S226 = _d_sin_1(&_S225);
            float _S227 = 2.0f * _S226.primal_0;
            float _S228 = (_S226.differential_0 * 2.0f * _S217.primal_0 - _S227 * _S217.differential_0) / (_S217.primal_0 * _S217.primal_0);
            k_1 = _S227 / _S217.primal_0;
            s_diff_k_1 = _S228;
        }
        float2  _S229 = _S179 * make_float2 (k_1);
        float2  _S230 = _S215 * make_float2 (k_1) + make_float2 (s_diff_k_1) * _S179;
        float u_12 = _S229.x;
        float s_diff_u_6 = _S230.x;
        float v_12 = _S229.y;
        float s_diff_v_6 = _S230.y;
        float _S231 = s_diff_u_6 * u_12;
        float _S232 = s_diff_v_6 * v_12;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float s_diff_r2_6 = _S231 + _S231 + (_S232 + _S232);
        float _S233 = _S183 + r2_12 * _S184;
        float _S234 = _S182 + r2_12 * _S233;
        float _S235 = _S181 + r2_12 * _S234;
        float2  _S236 = _S230 * make_float2 (1.0f + r2_12 * _S235) + make_float2 (s_diff_r2_6 * _S235 + (s_diff_r2_6 * _S234 + (s_diff_r2_6 * _S233 + s_diff_r2_6 * _S184 * r2_12) * r2_12) * r2_12) * _S229 + make_float2 (s_diff_u_6 * _S191 * v_12 + s_diff_v_6 * (_S191 * u_12) + (s_diff_r2_6 + (s_diff_u_6 * 2.0f * u_12 + s_diff_u_6 * (2.0f * u_12))) * _S186 + s_diff_r2_6 * _S187, s_diff_u_6 * _S192 * v_12 + s_diff_v_6 * (_S192 * u_12) + (s_diff_r2_6 + (s_diff_v_6 * 2.0f * v_12 + s_diff_v_6 * (2.0f * v_12))) * _S185 + s_diff_r2_6 * _S188);
        float2  _S237 = _S236 + make_float2 (_S236.x * _S189 + _S236.y * _S190, 0.0f);
        float fx_2 = intrins_2.x;
        float fy_2 = intrins_2.y;
        float _S238 = _S237.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->x) = _S237.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->x) = _S238;
        float2  _S239 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S240;
        (&_S240)->primal_0 = _S179;
        (&_S240)->differential_0 = _S239;
        DiffPair_float_0 _S241 = s_fwd_length_impl_0(&_S240);
        DiffPair_float_0 _S242;
        (&_S242)->primal_0 = _S241.primal_0;
        (&_S242)->differential_0 = _S241.differential_0;
        DiffPair_float_0 _S243;
        (&_S243)->primal_0 = _S218;
        (&_S243)->differential_0 = 0.0f;
        DiffPair_float_0 _S244 = _d_atan2_1(&_S242, &_S243);
        if((_S241.primal_0) < 9.99999997475242708e-07f)
        {
            float _S245 = _S244.differential_0 * _S244.primal_0;
            float _S246 = (0.0f - (_S245 + _S245) * 0.0416666679084301f) * _S180 / (_S180 * _S180);
            k_1 = (1.0f - _S244.primal_0 * _S244.primal_0 / 24.0f) / _S180;
            s_diff_k_1 = _S246;
        }
        else
        {
            float _S247 = _S244.differential_0 * 0.5f;
            DiffPair_float_0 _S248;
            (&_S248)->primal_0 = 0.5f * _S244.primal_0;
            (&_S248)->differential_0 = _S247;
            DiffPair_float_0 _S249 = _d_sin_1(&_S248);
            float _S250 = 2.0f * _S249.primal_0;
            float _S251 = (_S249.differential_0 * 2.0f * _S241.primal_0 - _S250 * _S241.differential_0) / (_S241.primal_0 * _S241.primal_0);
            k_1 = _S250 / _S241.primal_0;
            s_diff_k_1 = _S251;
        }
        float2  _S252 = _S179 * make_float2 (k_1);
        float2  _S253 = _S239 * make_float2 (k_1) + make_float2 (s_diff_k_1) * _S179;
        float u_13 = _S252.x;
        float s_diff_u_7 = _S253.x;
        float v_13 = _S252.y;
        float s_diff_v_7 = _S253.y;
        float _S254 = s_diff_u_7 * u_13;
        float _S255 = s_diff_v_7 * v_13;
        float r2_13 = u_13 * u_13 + v_13 * v_13;
        float s_diff_r2_7 = _S254 + _S254 + (_S255 + _S255);
        float _S256 = _S183 + r2_13 * _S184;
        float _S257 = _S182 + r2_13 * _S256;
        float _S258 = _S181 + r2_13 * _S257;
        float2  _S259 = _S253 * make_float2 (1.0f + r2_13 * _S258) + make_float2 (s_diff_r2_7 * _S258 + (s_diff_r2_7 * _S257 + (s_diff_r2_7 * _S256 + s_diff_r2_7 * _S184 * r2_13) * r2_13) * r2_13) * _S252 + make_float2 (s_diff_u_7 * _S191 * v_13 + s_diff_v_7 * (_S191 * u_13) + (s_diff_r2_7 + (s_diff_u_7 * 2.0f * u_13 + s_diff_u_7 * (2.0f * u_13))) * _S186 + s_diff_r2_7 * _S187, s_diff_u_7 * _S192 * v_13 + s_diff_v_7 * (_S192 * u_13) + (s_diff_r2_7 + (s_diff_v_7 * 2.0f * v_13 + s_diff_v_7 * (2.0f * v_13))) * _S185 + s_diff_r2_7 * _S188);
        float2  _S260 = _S259 + make_float2 (_S259.x * _S189 + _S259.y * _S190, 0.0f);
        float _S261 = _S260.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->y) = _S260.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->y) = _S261;
        float2  _S262 = make_float2 (0.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S263;
        (&_S263)->primal_0 = _S179;
        (&_S263)->differential_0 = _S262;
        DiffPair_float_0 _S264 = s_fwd_length_impl_0(&_S263);
        DiffPair_float_0 _S265;
        (&_S265)->primal_0 = _S264.primal_0;
        (&_S265)->differential_0 = _S264.differential_0;
        DiffPair_float_0 _S266;
        (&_S266)->primal_0 = _S180;
        (&_S266)->differential_0 = 1.0f;
        DiffPair_float_0 _S267 = _d_atan2_1(&_S265, &_S266);
        if((_S264.primal_0) < 9.99999997475242708e-07f)
        {
            float _S268 = _S267.differential_0 * _S267.primal_0;
            float _S269 = 1.0f - _S267.primal_0 * _S267.primal_0 / 24.0f;
            float _S270 = ((0.0f - (_S268 + _S268) * 0.0416666679084301f) * _S180 - _S269) / (_S180 * _S180);
            k_1 = _S269 / _S180;
            s_diff_k_1 = _S270;
        }
        else
        {
            float _S271 = _S267.differential_0 * 0.5f;
            DiffPair_float_0 _S272;
            (&_S272)->primal_0 = 0.5f * _S267.primal_0;
            (&_S272)->differential_0 = _S271;
            DiffPair_float_0 _S273 = _d_sin_1(&_S272);
            float _S274 = 2.0f * _S273.primal_0;
            float _S275 = (_S273.differential_0 * 2.0f * _S264.primal_0 - _S274 * _S264.differential_0) / (_S264.primal_0 * _S264.primal_0);
            k_1 = _S274 / _S264.primal_0;
            s_diff_k_1 = _S275;
        }
        float2  _S276 = _S179 * make_float2 (k_1);
        float2  _S277 = make_float2 (s_diff_k_1) * _S179;
        float u_14 = _S276.x;
        float s_diff_u_8 = _S277.x;
        float v_14 = _S276.y;
        float s_diff_v_8 = _S277.y;
        float _S278 = s_diff_u_8 * u_14;
        float _S279 = s_diff_v_8 * v_14;
        float r2_14 = u_14 * u_14 + v_14 * v_14;
        float s_diff_r2_8 = _S278 + _S278 + (_S279 + _S279);
        float _S280 = _S183 + r2_14 * _S184;
        float _S281 = _S182 + r2_14 * _S280;
        float _S282 = _S181 + r2_14 * _S281;
        float2  _S283 = _S277 * make_float2 (1.0f + r2_14 * _S282) + make_float2 (s_diff_r2_8 * _S282 + (s_diff_r2_8 * _S281 + (s_diff_r2_8 * _S280 + s_diff_r2_8 * _S184 * r2_14) * r2_14) * r2_14) * _S276 + make_float2 (s_diff_u_8 * _S191 * v_14 + s_diff_v_8 * (_S191 * u_14) + (s_diff_r2_8 + (s_diff_u_8 * 2.0f * u_14 + s_diff_u_8 * (2.0f * u_14))) * _S186 + s_diff_r2_8 * _S187, s_diff_u_8 * _S192 * v_14 + s_diff_v_8 * (_S192 * u_14) + (s_diff_r2_8 + (s_diff_v_8 * 2.0f * v_14 + s_diff_v_8 * (2.0f * v_14))) * _S185 + s_diff_r2_8 * _S188);
        float2  _S284 = _S283 + make_float2 (_S283.x * _S189 + _S283.y * _S190, 0.0f);
        float _S285 = _S284.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->z) = _S284.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->z) = _S285;
        *cov2d_2 = mul_4(mul_3(J_3, cov3d_2), transpose_1(J_3));
        _S178 = true;
        break;
    }
    return _S178;
}

inline __device__ bool equirect_proj_3dgs_nav(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float4  intrins_3, FixedArray<float, 10>  dist_coeffs_3, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    float _S286 = mean3d_3.x;
    float _S287 = mean3d_3.z;
    float _S288 = mean3d_3.y;
    float2  _S289 = float2 {mean3d_3.x, mean3d_3.z};
    float fx_3 = intrins_3.x;
    float fy_3 = intrins_3.y;
    *mean2d_3 = make_float2 (fx_3 * (F32_atan2((_S286), (_S287))) + intrins_3.z, fy_3 * (F32_atan2((_S288), (length_0(_S289)))) + intrins_3.w);
    DiffPair_float_0 _S290;
    (&_S290)->primal_0 = _S286;
    (&_S290)->differential_0 = 1.0f;
    DiffPair_float_0 _S291;
    (&_S291)->primal_0 = _S287;
    (&_S291)->differential_0 = 0.0f;
    DiffPair_float_0 _S292 = _d_atan2_1(&_S290, &_S291);
    float2  _S293 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S294;
    (&_S294)->primal_0 = _S289;
    (&_S294)->differential_0 = _S293;
    DiffPair_float_0 _S295 = s_fwd_length_impl_0(&_S294);
    DiffPair_float_0 _S296;
    (&_S296)->primal_0 = _S288;
    (&_S296)->differential_0 = 0.0f;
    DiffPair_float_0 _S297;
    (&_S297)->primal_0 = _S295.primal_0;
    (&_S297)->differential_0 = _S295.differential_0;
    DiffPair_float_0 _S298 = _d_atan2_1(&_S296, &_S297);
    float _S299 = _S298.differential_0 * fy_3;
    Matrix<float, 2, 3>  J_4;
    *&(((&J_4)->rows + (int(0)))->x) = _S292.differential_0 * fx_3;
    *&(((&J_4)->rows + (int(1)))->x) = _S299;
    DiffPair_float_0 _S300;
    (&_S300)->primal_0 = _S286;
    (&_S300)->differential_0 = 0.0f;
    DiffPair_float_0 _S301;
    (&_S301)->primal_0 = _S287;
    (&_S301)->differential_0 = 0.0f;
    DiffPair_float_0 _S302 = _d_atan2_1(&_S300, &_S301);
    float2  _S303 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S304;
    (&_S304)->primal_0 = _S289;
    (&_S304)->differential_0 = _S303;
    DiffPair_float_0 _S305 = s_fwd_length_impl_0(&_S304);
    DiffPair_float_0 _S306;
    (&_S306)->primal_0 = _S288;
    (&_S306)->differential_0 = 1.0f;
    DiffPair_float_0 _S307;
    (&_S307)->primal_0 = _S305.primal_0;
    (&_S307)->differential_0 = _S305.differential_0;
    DiffPair_float_0 _S308 = _d_atan2_1(&_S306, &_S307);
    float _S309 = _S308.differential_0 * fy_3;
    *&(((&J_4)->rows + (int(0)))->y) = _S302.differential_0 * fx_3;
    *&(((&J_4)->rows + (int(1)))->y) = _S309;
    DiffPair_float_0 _S310;
    (&_S310)->primal_0 = _S286;
    (&_S310)->differential_0 = 0.0f;
    DiffPair_float_0 _S311;
    (&_S311)->primal_0 = _S287;
    (&_S311)->differential_0 = 1.0f;
    DiffPair_float_0 _S312 = _d_atan2_1(&_S310, &_S311);
    float2  _S313 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S314;
    (&_S314)->primal_0 = _S289;
    (&_S314)->differential_0 = _S313;
    DiffPair_float_0 _S315 = s_fwd_length_impl_0(&_S314);
    DiffPair_float_0 _S316;
    (&_S316)->primal_0 = _S288;
    (&_S316)->differential_0 = 0.0f;
    DiffPair_float_0 _S317;
    (&_S317)->primal_0 = _S315.primal_0;
    (&_S317)->differential_0 = _S315.differential_0;
    DiffPair_float_0 _S318 = _d_atan2_1(&_S316, &_S317);
    float _S319 = _S318.differential_0 * fy_3;
    *&(((&J_4)->rows + (int(0)))->z) = _S312.differential_0 * fx_3;
    *&(((&J_4)->rows + (int(1)))->z) = _S319;
    *cov2d_3 = mul_4(mul_3(J_4, cov3d_3), transpose_1(J_4));
    return true;
}

inline __device__ float2  min_0(float2  x_11, float2  y_3)
{
    float2  result_13;
    int i_6 = int(0);
    for(;;)
    {
        if(i_6 < int(2))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_13, i_6) = (F32_min((_slang_vector_get_element(x_11, i_6)), (_slang_vector_get_element(y_3, i_6))));
        i_6 = i_6 + int(1);
    }
    return result_13;
}

inline __device__ float2  max_0(float2  x_12, float2  y_4)
{
    float2  result_14;
    int i_7 = int(0);
    for(;;)
    {
        if(i_7 < int(2))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_14, i_7) = (F32_max((_slang_vector_get_element(x_12, i_7)), (_slang_vector_get_element(y_4, i_7))));
        i_7 = i_7 + int(1);
    }
    return result_14;
}

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_6, float3  dOut_11)
{
    float _S320 = (*left_6).primal_0.rows[int(0)].x * dOut_11.x;
    Matrix<float, 3, 3>  left_d_result_3;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = (*right_6).primal_0.x * dOut_11.x;
    float sum_6 = _S320 + (*left_6).primal_0.rows[int(1)].x * dOut_11.y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = (*right_6).primal_0.x * dOut_11.y;
    float sum_7 = sum_6 + (*left_6).primal_0.rows[int(2)].x * dOut_11.z;
    *&(((&left_d_result_3)->rows + (int(2)))->x) = (*right_6).primal_0.x * dOut_11.z;
    float3  right_d_result_3;
    *&((&right_d_result_3)->x) = sum_7;
    float _S321 = (*left_6).primal_0.rows[int(0)].y * dOut_11.x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = (*right_6).primal_0.y * dOut_11.x;
    float sum_8 = _S321 + (*left_6).primal_0.rows[int(1)].y * dOut_11.y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = (*right_6).primal_0.y * dOut_11.y;
    float sum_9 = sum_8 + (*left_6).primal_0.rows[int(2)].y * dOut_11.z;
    *&(((&left_d_result_3)->rows + (int(2)))->y) = (*right_6).primal_0.y * dOut_11.z;
    *&((&right_d_result_3)->y) = sum_9;
    float _S322 = (*left_6).primal_0.rows[int(0)].z * dOut_11.x;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = (*right_6).primal_0.z * dOut_11.x;
    float sum_10 = _S322 + (*left_6).primal_0.rows[int(1)].z * dOut_11.y;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = (*right_6).primal_0.z * dOut_11.y;
    float sum_11 = sum_10 + (*left_6).primal_0.rows[int(2)].z * dOut_11.z;
    *&(((&left_d_result_3)->rows + (int(2)))->z) = (*right_6).primal_0.z * dOut_11.z;
    *&((&right_d_result_3)->z) = sum_11;
    left_6->primal_0 = (*left_6).primal_0;
    left_6->differential_0 = left_d_result_3;
    right_6->primal_0 = (*right_6).primal_0;
    right_6->differential_0 = right_d_result_3;
    return;
}

inline __device__ float3  mul_6(Matrix<float, 3, 3>  left_7, float3  right_7)
{
    float3  result_15;
    int i_8 = int(0);
    for(;;)
    {
        if(i_8 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_12 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_13 = sum_12 + _slang_vector_get_element(left_7.rows[i_8], j_0) * _slang_vector_get_element(right_7, j_0);
            j_0 = j_0 + int(1);
            sum_12 = sum_13;
        }
        *_slang_vector_get_element_ptr(&result_15, i_8) = sum_12;
        i_8 = i_8 + int(1);
    }
    return result_15;
}

inline __device__ float2  mul_7(Matrix<float, 2, 2>  left_8, float2  right_8)
{
    float2  result_16;
    int i_9 = int(0);
    for(;;)
    {
        if(i_9 < int(2))
        {
        }
        else
        {
            break;
        }
        int j_1 = int(0);
        float sum_14 = 0.0f;
        for(;;)
        {
            if(j_1 < int(2))
            {
            }
            else
            {
                break;
            }
            float sum_15 = sum_14 + _slang_vector_get_element(left_8.rows[i_9], j_1) * _slang_vector_get_element(right_8, j_1);
            j_1 = j_1 + int(1);
            sum_14 = sum_15;
        }
        *_slang_vector_get_element_ptr(&result_16, i_9) = sum_14;
        i_9 = i_9 + int(1);
    }
    return result_16;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_12, float dOut_12)
{
    float _S323 = (F32_exp(((*dpx_12).primal_0))) * dOut_12;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S323;
    return;
}

inline __device__ float3  exp_0(float3  x_13)
{
    float3  result_17;
    int i_10 = int(0);
    for(;;)
    {
        if(i_10 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_17, i_10) = (F32_exp((_slang_vector_get_element(x_13, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_17;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_13, float3  dOut_13)
{
    float3  _S324 = exp_0((*dpx_13).primal_0) * dOut_13;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S324;
    return;
}

inline __device__ float4  normalize_0(float4  x_14)
{
    return x_14 / make_float4 (length_2(x_14));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_14, float dOut_14)
{
    float _S325 = 1.0f / (*dpx_14).primal_0 * dOut_14;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S325;
    return;
}

inline __device__ float view_radius_3dgs_0(float3  mean_0, float3  log_scale_0, float logit_opacity_0, float3  campos_0)
{
    float radius_0 = (F32_exp(((F32_max((log_scale_0.x), ((F32_max((log_scale_0.y), (log_scale_0.z))))))))) * (F32_sqrt((2.0f * (F32_log(((F32_max((255.0f / (1.0f + (F32_exp((- logit_opacity_0))))), (1.0f)))))))));
    float dist_0 = length_1(mean_0 - campos_0);
    return radius_0 / ((F32_max((dist_0), (radius_0))) + (F32_sqrt(((F32_max((dist_0 * dist_0 - radius_0 * radius_0), (0.0f)))))));
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_1, float4  quat_0, float3  scale_0, float in_opacity_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_4, float fy_4, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_4, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * sorting_depth_0, float * radius_1, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0)
{
    float2  _S326;
    for(;;)
    {
        float3  mean_c_0 = mul_6(R_0, mean_1) + t_0;
        float _S327 = mean_c_0.z;
        *depth_0 = length_1(mean_c_0);
        if(_S327 <= 0.0f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        *opacity_0 = 1.0f / (1.0f + (F32_exp((- in_opacity_0))));
        bool _S328;
        float4  _S329 = normalize_0(quat_0);
        float3  _S330 = exp_0(scale_0);
        float x_15 = _S329.y;
        float x2_0 = x_15 * x_15;
        float y2_0 = _S329.z * _S329.z;
        float z2_0 = _S329.w * _S329.w;
        float xy_0 = _S329.y * _S329.z;
        float xz_0 = _S329.y * _S329.w;
        float yz_0 = _S329.z * _S329.w;
        float wx_0 = _S329.x * _S329.y;
        float wy_0 = _S329.x * _S329.z;
        float wz_0 = _S329.x * _S329.w;
        Matrix<float, 3, 3>  M_0 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0))), makeMatrix<float, 3, 3> (_S330.x, 0.0f, 0.0f, 0.0f, _S330.y, 0.0f, 0.0f, 0.0f, _S330.z));
        Matrix<float, 3, 3>  _S331 = transpose_3(R_0);
        Matrix<float, 3, 3>  covar_c_0 = mul_5(mul_5(R_0, mul_5(M_0, transpose_3(M_0))), _S331);
        for(;;)
        {
            for(;;)
            {
                float2  _S332 = float2 {mean_c_0.x, mean_c_0.y};
                _S326 = _S332;
                *mean2d_4 = _S332 / make_float2 (_S327);
                if(_S327 < 0.0f)
                {
                    _S328 = true;
                }
                else
                {
                    float u_15 = (*mean2d_4).x;
                    float v_15 = (*mean2d_4).y;
                    float _S333 = u_15 + u_15;
                    float r2_15 = u_15 * u_15 + v_15 * v_15;
                    float _S334 = dist_coeffs_4[int(2)] + r2_15 * dist_coeffs_4[int(3)];
                    float _S335 = dist_coeffs_4[int(1)] + r2_15 * _S334;
                    float _S336 = dist_coeffs_4[int(0)] + r2_15 * _S335;
                    float radial_3 = 1.0f + r2_15 * _S336;
                    float _S337 = 2.0f * dist_coeffs_4[int(4)];
                    float _S338 = 2.0f * u_15;
                    float _S339 = 2.0f * dist_coeffs_4[int(5)];
                    float _S340 = 2.0f * v_15;
                    float2  _S341 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (_S333 * _S336 + (_S333 * _S335 + (_S333 * _S334 + _S333 * dist_coeffs_4[int(3)] * r2_15) * r2_15) * r2_15) * *mean2d_4 + make_float2 (_S337 * v_15 + (_S333 + (_S338 + _S338)) * dist_coeffs_4[int(5)] + _S333 * dist_coeffs_4[int(6)], _S339 * v_15 + _S333 * dist_coeffs_4[int(4)] + _S333 * dist_coeffs_4[int(7)]);
                    float _S342 = v_15 + v_15;
                    float2  _S343 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (_S342 * _S336 + (_S342 * _S335 + (_S342 * _S334 + _S342 * dist_coeffs_4[int(3)] * r2_15) * r2_15) * r2_15) * *mean2d_4 + make_float2 (_S337 * u_15 + _S342 * dist_coeffs_4[int(5)] + _S342 * dist_coeffs_4[int(6)], _S339 * u_15 + (_S342 + (_S340 + _S340)) * dist_coeffs_4[int(4)] + _S342 * dist_coeffs_4[int(7)]);
                    Matrix<float, 2, 2>  _S344 = transpose_0(makeMatrix<float, 2, 2> (_S341 + make_float2 (_S341.x * dist_coeffs_4[int(8)] + _S341.y * dist_coeffs_4[int(9)], 0.0f), _S343 + make_float2 (_S343.x * dist_coeffs_4[int(8)] + _S343.y * dist_coeffs_4[int(9)], 0.0f)));
                    _S328 = !((F32_min((determinant_0(_S344)), ((F32_min((_S344.rows[int(0)].x), (_S344.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S328)
                {
                    break;
                }
                float u_16 = (*mean2d_4).x;
                float v_16 = (*mean2d_4).y;
                float r2_16 = u_16 * u_16 + v_16 * v_16;
                float2  _S345 = *mean2d_4 * make_float2 (1.0f + r2_16 * (dist_coeffs_4[int(0)] + r2_16 * (dist_coeffs_4[int(1)] + r2_16 * (dist_coeffs_4[int(2)] + r2_16 * dist_coeffs_4[int(3)])))) + make_float2 (2.0f * dist_coeffs_4[int(4)] * u_16 * v_16 + dist_coeffs_4[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + dist_coeffs_4[int(6)] * r2_16, 2.0f * dist_coeffs_4[int(5)] * u_16 * v_16 + dist_coeffs_4[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + dist_coeffs_4[int(7)] * r2_16);
                float2  _S346 = _S345 + make_float2 (dist_coeffs_4[int(8)] * _S345.x + dist_coeffs_4[int(9)] * _S345.y, 0.0f);
                *mean2d_4 = make_float2 (fx_4 * _S346.x + cx_1, fy_4 * _S346.y + cy_1);
                break;
            }
            if(!!_S328)
            {
                _S328 = false;
                break;
            }
            Matrix<float, 2, 3>  J_5;
            float2  _S347 = _S326 / make_float2 (_S327);
            float _S348 = _S327 * _S327;
            float2  _S349 = make_float2 (1.0f, 0.0f) * make_float2 (_S327) / make_float2 (_S348);
            float u_17 = _S347.x;
            float s_diff_u_9 = _S349.x;
            float v_17 = _S347.y;
            float s_diff_v_9 = _S349.y;
            float _S350 = s_diff_u_9 * u_17;
            float _S351 = s_diff_v_9 * v_17;
            float r2_17 = u_17 * u_17 + v_17 * v_17;
            float s_diff_r2_9 = _S350 + _S350 + (_S351 + _S351);
            float _S352 = dist_coeffs_4[int(2)] + r2_17 * dist_coeffs_4[int(3)];
            float _S353 = dist_coeffs_4[int(1)] + r2_17 * _S352;
            float _S354 = dist_coeffs_4[int(0)] + r2_17 * _S353;
            float _S355 = 2.0f * dist_coeffs_4[int(4)];
            float _S356 = 2.0f * dist_coeffs_4[int(5)];
            float2  _S357 = _S349 * make_float2 (1.0f + r2_17 * _S354) + make_float2 (s_diff_r2_9 * _S354 + (s_diff_r2_9 * _S353 + (s_diff_r2_9 * _S352 + s_diff_r2_9 * dist_coeffs_4[int(3)] * r2_17) * r2_17) * r2_17) * _S347 + make_float2 (s_diff_u_9 * _S355 * v_17 + s_diff_v_9 * (_S355 * u_17) + (s_diff_r2_9 + (s_diff_u_9 * 2.0f * u_17 + s_diff_u_9 * (2.0f * u_17))) * dist_coeffs_4[int(5)] + s_diff_r2_9 * dist_coeffs_4[int(6)], s_diff_u_9 * _S356 * v_17 + s_diff_v_9 * (_S356 * u_17) + (s_diff_r2_9 + (s_diff_v_9 * 2.0f * v_17 + s_diff_v_9 * (2.0f * v_17))) * dist_coeffs_4[int(4)] + s_diff_r2_9 * dist_coeffs_4[int(7)]);
            float2  _S358 = _S357 + make_float2 (_S357.x * dist_coeffs_4[int(8)] + _S357.y * dist_coeffs_4[int(9)], 0.0f);
            float _S359 = _S358.y * fy_4;
            Matrix<float, 2, 3>  J_6;
            *&(((&J_6)->rows + (int(0)))->x) = _S358.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->x) = _S359;
            float2  _S360 = _S326 / make_float2 (_S327);
            float2  _S361 = make_float2 (0.0f, 1.0f) * make_float2 (_S327) / make_float2 (_S348);
            float u_18 = _S360.x;
            float s_diff_u_10 = _S361.x;
            float v_18 = _S360.y;
            float s_diff_v_10 = _S361.y;
            float _S362 = s_diff_u_10 * u_18;
            float _S363 = s_diff_v_10 * v_18;
            float r2_18 = u_18 * u_18 + v_18 * v_18;
            float s_diff_r2_10 = _S362 + _S362 + (_S363 + _S363);
            float _S364 = dist_coeffs_4[int(2)] + r2_18 * dist_coeffs_4[int(3)];
            float _S365 = dist_coeffs_4[int(1)] + r2_18 * _S364;
            float _S366 = dist_coeffs_4[int(0)] + r2_18 * _S365;
            float2  _S367 = _S361 * make_float2 (1.0f + r2_18 * _S366) + make_float2 (s_diff_r2_10 * _S366 + (s_diff_r2_10 * _S365 + (s_diff_r2_10 * _S364 + s_diff_r2_10 * dist_coeffs_4[int(3)] * r2_18) * r2_18) * r2_18) * _S360 + make_float2 (s_diff_u_10 * _S355 * v_18 + s_diff_v_10 * (_S355 * u_18) + (s_diff_r2_10 + (s_diff_u_10 * 2.0f * u_18 + s_diff_u_10 * (2.0f * u_18))) * dist_coeffs_4[int(5)] + s_diff_r2_10 * dist_coeffs_4[int(6)], s_diff_u_10 * _S356 * v_18 + s_diff_v_10 * (_S356 * u_18) + (s_diff_r2_10 + (s_diff_v_10 * 2.0f * v_18 + s_diff_v_10 * (2.0f * v_18))) * dist_coeffs_4[int(4)] + s_diff_r2_10 * dist_coeffs_4[int(7)]);
            float2  _S368 = _S367 + make_float2 (_S367.x * dist_coeffs_4[int(8)] + _S367.y * dist_coeffs_4[int(9)], 0.0f);
            float _S369 = _S368.y * fy_4;
            *&(((&J_6)->rows + (int(0)))->y) = _S368.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->y) = _S369;
            float2  _S370 = _S326 / make_float2 (_S327);
            float2  _S371 = (make_float2 (0.0f, 0.0f) - _S326) / make_float2 (_S348);
            float u_19 = _S370.x;
            float s_diff_u_11 = _S371.x;
            float v_19 = _S370.y;
            float s_diff_v_11 = _S371.y;
            float _S372 = s_diff_u_11 * u_19;
            float _S373 = s_diff_v_11 * v_19;
            float r2_19 = u_19 * u_19 + v_19 * v_19;
            float s_diff_r2_11 = _S372 + _S372 + (_S373 + _S373);
            float _S374 = dist_coeffs_4[int(2)] + r2_19 * dist_coeffs_4[int(3)];
            float _S375 = dist_coeffs_4[int(1)] + r2_19 * _S374;
            float _S376 = dist_coeffs_4[int(0)] + r2_19 * _S375;
            float2  _S377 = _S371 * make_float2 (1.0f + r2_19 * _S376) + make_float2 (s_diff_r2_11 * _S376 + (s_diff_r2_11 * _S375 + (s_diff_r2_11 * _S374 + s_diff_r2_11 * dist_coeffs_4[int(3)] * r2_19) * r2_19) * r2_19) * _S370 + make_float2 (s_diff_u_11 * _S355 * v_19 + s_diff_v_11 * (_S355 * u_19) + (s_diff_r2_11 + (s_diff_u_11 * 2.0f * u_19 + s_diff_u_11 * (2.0f * u_19))) * dist_coeffs_4[int(5)] + s_diff_r2_11 * dist_coeffs_4[int(6)], s_diff_u_11 * _S356 * v_19 + s_diff_v_11 * (_S356 * u_19) + (s_diff_r2_11 + (s_diff_v_11 * 2.0f * v_19 + s_diff_v_11 * (2.0f * v_19))) * dist_coeffs_4[int(4)] + s_diff_r2_11 * dist_coeffs_4[int(7)]);
            float2  _S378 = _S377 + make_float2 (_S377.x * dist_coeffs_4[int(8)] + _S377.y * dist_coeffs_4[int(9)], 0.0f);
            float _S379 = _S378.y * fy_4;
            *&(((&J_6)->rows + (int(0)))->z) = _S378.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->z) = _S379;
            J_5 = J_6;
            float _S380 = float(image_width_0);
            float _S381 = 0.30000001192092896f * (0.5f * _S380);
            float lim_x_pos_1 = _S380 + _S381;
            float rz_1 = 1.0f / _S327;
            float _S382 = - _S381;
            float max_Jyz_1 = - (_S382 - cy_1) * rz_1;
            float min_Jyz_1 = - (lim_x_pos_1 - cy_1) * rz_1;
            *&(((&J_5)->rows + (int(0)))->z) = clamp_0(*&(((&J_5)->rows + (int(0)))->z), - (lim_x_pos_1 - cx_1) * rz_1, - (_S382 - cx_1) * rz_1);
            *&(((&J_5)->rows + (int(1)))->z) = clamp_0(*&(((&J_5)->rows + (int(1)))->z), min_Jyz_1, max_Jyz_1);
            covar2d_0 = mul_4(mul_3(J_5, covar_c_0), transpose_1(J_5));
            _S328 = true;
            break;
        }
        if(!(true & _S328))
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float eps2d_0;
        if(antialiased_0)
        {
            eps2d_0 = 0.10000000149011612f;
        }
        else
        {
            eps2d_0 = 0.30000001192092896f;
        }
        float det_orig_0 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S383 = *&(((&covar2d_0)->rows + (int(0)))->x) + eps2d_0;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S383;
        float _S384 = *&(((&covar2d_0)->rows + (int(1)))->y) + eps2d_0;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S384;
        float det_blur_0 = _S383 * _S384 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
        if(det_blur_0 <= 0.0f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float invdet_0 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S385 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_0, - covar2d_0.rows[int(0)].y * invdet_0, - covar2d_0.rows[int(1)].x * invdet_0, covar2d_0.rows[int(0)].x * invdet_0);
        if(antialiased_0)
        {
            *opacity_0 = *opacity_0 * compensation_0;
        }
        if((*opacity_0) < 0.00392156885936856f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float _S386 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_0 / 0.00392156885936856f)))))))));
        float radius_x_0 = _S386 * (F32_sqrt((covar2d_0[int(0)].x)));
        float radius_y_0 = _S386 * (F32_sqrt((covar2d_0[int(1)].y)));
        float _S387 = (*mean2d_4).x - radius_x_0;
        float _S388 = (*mean2d_4).x + radius_x_0;
        float _S389 = (*mean2d_4).y - radius_y_0;
        float _S390 = (*mean2d_4).y + radius_y_0;
        if(_S388 <= 0.0f)
        {
            _S328 = true;
        }
        else
        {
            _S328 = _S387 >= float(image_width_0);
        }
        if(_S328)
        {
            _S328 = true;
        }
        else
        {
            _S328 = _S390 <= 0.0f;
        }
        if(_S328)
        {
            _S328 = true;
        }
        else
        {
            _S328 = _S389 >= float(image_height_0);
        }
        if(_S328)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (_S387, _S389, _S388, _S390);
        *sorting_depth_0 = _S327;
        *conic_0 = make_float3 (_S385.rows[int(0)].x, _S385.rows[int(0)].y, _S385.rows[int(1)].y);
        *radius_1 = view_radius_3dgs_0(mean_1, scale_0, in_opacity_0, - mul_6(_S331, t_0));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_2, float4  quat_1, float3  scale_1, float in_opacity_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_5, float fy_5, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_5, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * sorting_depth_1, float * radius_2, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1)
{
    float2  _S391;
    float _S392;
    float _S393;
    float _S394;
    float _S395;
    float _S396;
    float _S397;
    float _S398;
    float _S399;
    float _S400;
    float _S401;
    float _S402;
    float _S403;
    float _S404;
    bool _S405;
    for(;;)
    {
        float3  mean_c_1 = mul_6(R_1, mean_2) + t_1;
        float _S406 = length_1(mean_c_1);
        *depth_1 = _S406;
        if(_S406 <= 0.0f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_1;
        *opacity_1 = 1.0f / (1.0f + (F32_exp((- in_opacity_1))));
        bool is_valid_0;
        float eps2d_1;
        float4  _S407 = normalize_0(quat_1);
        float3  _S408 = exp_0(scale_1);
        float x_16 = _S407.y;
        float x2_1 = x_16 * x_16;
        float y2_1 = _S407.z * _S407.z;
        float z2_1 = _S407.w * _S407.w;
        float xy_1 = _S407.y * _S407.z;
        float xz_1 = _S407.y * _S407.w;
        float yz_1 = _S407.z * _S407.w;
        float wx_1 = _S407.x * _S407.y;
        float wy_1 = _S407.x * _S407.z;
        float wz_1 = _S407.x * _S407.w;
        Matrix<float, 3, 3>  M_1 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (_S408.x, 0.0f, 0.0f, 0.0f, _S408.y, 0.0f, 0.0f, 0.0f, _S408.z));
        Matrix<float, 3, 3>  _S409 = transpose_3(R_1);
        Matrix<float, 3, 3>  covar_c_1 = mul_5(mul_5(R_1, mul_5(M_1, transpose_3(M_1))), _S409);
        for(;;)
        {
            float k_2;
            for(;;)
            {
                float2  _S410 = float2 {mean_c_1.x, mean_c_1.y};
                _S391 = _S410;
                float r_9 = length_0(_S410);
                float _S411 = mean_c_1.z;
                _S392 = _S411;
                float theta_2 = (F32_atan2((r_9), (_S411)));
                if(theta_2 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_2 * theta_2 / 3.0f) / _S411;
                }
                else
                {
                    k_2 = theta_2 / r_9;
                }
                float2  _S412 = _S410 * make_float2 (k_2);
                *mean2d_5 = _S412;
                float2  _S413 = make_float2 (1.0f, 0.0f);
                _S393 = dist_coeffs_5[int(0)];
                _S394 = dist_coeffs_5[int(1)];
                _S395 = dist_coeffs_5[int(2)];
                _S396 = dist_coeffs_5[int(3)];
                _S397 = dist_coeffs_5[int(4)];
                _S398 = dist_coeffs_5[int(5)];
                _S399 = dist_coeffs_5[int(6)];
                _S400 = dist_coeffs_5[int(7)];
                _S401 = dist_coeffs_5[int(8)];
                _S402 = dist_coeffs_5[int(9)];
                float u_20 = _S412.x;
                float v_20 = _S412.y;
                float _S414 = u_20 + u_20;
                float r2_20 = u_20 * u_20 + v_20 * v_20;
                float _S415 = dist_coeffs_5[int(2)] + r2_20 * dist_coeffs_5[int(3)];
                float _S416 = dist_coeffs_5[int(1)] + r2_20 * _S415;
                float _S417 = dist_coeffs_5[int(0)] + r2_20 * _S416;
                float _S418 = _S414 * _S417 + (_S414 * _S416 + (_S414 * _S415 + _S414 * dist_coeffs_5[int(3)] * r2_20) * r2_20) * r2_20;
                float radial_4 = 1.0f + r2_20 * _S417;
                float _S419 = 2.0f * dist_coeffs_5[int(4)];
                _S403 = _S419;
                float _S420 = _S419 * u_20;
                float _S421 = 2.0f * u_20;
                float s_diff_du_2 = _S419 * v_20 + (_S414 + (_S421 + _S421)) * dist_coeffs_5[int(5)] + _S414 * dist_coeffs_5[int(6)];
                float _S422 = 2.0f * dist_coeffs_5[int(5)];
                _S404 = _S422;
                float _S423 = 2.0f * v_20;
                float2  _S424 = _S413 * make_float2 (radial_4) + make_float2 (_S418) * _S412 + make_float2 (s_diff_du_2, _S422 * v_20 + _S414 * dist_coeffs_5[int(4)] + _S414 * dist_coeffs_5[int(7)]);
                float _S425 = v_20 + v_20;
                float2  _S426 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (_S425 * _S417 + (_S425 * _S416 + (_S425 * _S415 + _S425 * dist_coeffs_5[int(3)] * r2_20) * r2_20) * r2_20) * _S412 + make_float2 (_S420 + _S425 * dist_coeffs_5[int(5)] + _S425 * dist_coeffs_5[int(6)], _S422 * u_20 + (_S425 + (_S423 + _S423)) * dist_coeffs_5[int(4)] + _S425 * dist_coeffs_5[int(7)]);
                Matrix<float, 2, 2>  _S427 = transpose_0(makeMatrix<float, 2, 2> (_S424 + make_float2 (_S424.x * dist_coeffs_5[int(8)] + _S424.y * dist_coeffs_5[int(9)], 0.0f), _S426 + make_float2 (_S426.x * dist_coeffs_5[int(8)] + _S426.y * dist_coeffs_5[int(9)], 0.0f)));
                bool _S428 = !((F32_min((determinant_0(_S427)), ((F32_min((_S427.rows[int(0)].x), (_S427.rows[int(1)].y)))))) > 0.0f);
                _S405 = _S428;
                if(_S428)
                {
                    break;
                }
                float u_21 = (*mean2d_5).x;
                float v_21 = (*mean2d_5).y;
                float r2_21 = u_21 * u_21 + v_21 * v_21;
                float2  _S429 = *mean2d_5 * make_float2 (1.0f + r2_21 * (dist_coeffs_5[int(0)] + r2_21 * (dist_coeffs_5[int(1)] + r2_21 * (dist_coeffs_5[int(2)] + r2_21 * dist_coeffs_5[int(3)])))) + make_float2 (_S419 * u_21 * v_21 + dist_coeffs_5[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + dist_coeffs_5[int(6)] * r2_21, _S422 * u_21 * v_21 + dist_coeffs_5[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + dist_coeffs_5[int(7)] * r2_21);
                float2  _S430 = _S429 + make_float2 (dist_coeffs_5[int(8)] * _S429.x + dist_coeffs_5[int(9)] * _S429.y, 0.0f);
                *mean2d_5 = make_float2 (fx_5 * _S430.x + cx_2, fy_5 * _S430.y + cy_2);
                break;
            }
            if(!!_S405)
            {
                is_valid_0 = false;
                break;
            }
            Matrix<float, 2, 3>  J_7;
            float2  _S431 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S432;
            (&_S432)->primal_0 = _S391;
            (&_S432)->differential_0 = _S431;
            DiffPair_float_0 _S433 = s_fwd_length_impl_0(&_S432);
            float _S434 = _S392;
            DiffPair_float_0 _S435;
            (&_S435)->primal_0 = _S433.primal_0;
            (&_S435)->differential_0 = _S433.differential_0;
            DiffPair_float_0 _S436;
            (&_S436)->primal_0 = _S392;
            (&_S436)->differential_0 = 0.0f;
            DiffPair_float_0 _S437 = _d_atan2_1(&_S435, &_S436);
            if((_S437.primal_0) < 0.00100000004749745f)
            {
                float _S438 = _S437.differential_0 * _S437.primal_0;
                float _S439 = (0.0f - (_S438 + _S438) * 0.3333333432674408f) * _S392 / (_S392 * _S392);
                k_2 = (1.0f - _S437.primal_0 * _S437.primal_0 / 3.0f) / _S392;
                eps2d_1 = _S439;
            }
            else
            {
                float _S440 = (_S437.differential_0 * _S433.primal_0 - _S437.primal_0 * _S433.differential_0) / (_S433.primal_0 * _S433.primal_0);
                k_2 = _S437.primal_0 / _S433.primal_0;
                eps2d_1 = _S440;
            }
            float2  _S441 = _S391 * make_float2 (k_2);
            float2  _S442 = _S431 * make_float2 (k_2) + make_float2 (eps2d_1) * _S391;
            float u_22 = _S441.x;
            float s_diff_u_12 = _S442.x;
            float v_22 = _S441.y;
            float s_diff_v_12 = _S442.y;
            float _S443 = s_diff_u_12 * u_22;
            float _S444 = s_diff_v_12 * v_22;
            float r2_22 = u_22 * u_22 + v_22 * v_22;
            float s_diff_r2_12 = _S443 + _S443 + (_S444 + _S444);
            float _S445 = _S395 + r2_22 * _S396;
            float _S446 = _S394 + r2_22 * _S445;
            float _S447 = _S393 + r2_22 * _S446;
            float2  _S448 = _S442 * make_float2 (1.0f + r2_22 * _S447) + make_float2 (s_diff_r2_12 * _S447 + (s_diff_r2_12 * _S446 + (s_diff_r2_12 * _S445 + s_diff_r2_12 * _S396 * r2_22) * r2_22) * r2_22) * _S441 + make_float2 (s_diff_u_12 * _S403 * v_22 + s_diff_v_12 * (_S403 * u_22) + (s_diff_r2_12 + (s_diff_u_12 * 2.0f * u_22 + s_diff_u_12 * (2.0f * u_22))) * _S398 + s_diff_r2_12 * _S399, s_diff_u_12 * _S404 * v_22 + s_diff_v_12 * (_S404 * u_22) + (s_diff_r2_12 + (s_diff_v_12 * 2.0f * v_22 + s_diff_v_12 * (2.0f * v_22))) * _S397 + s_diff_r2_12 * _S400);
            float2  _S449 = _S448 + make_float2 (_S448.x * _S401 + _S448.y * _S402, 0.0f);
            float _S450 = _S449.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->x) = _S449.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->x) = _S450;
            float2  _S451 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S452;
            (&_S452)->primal_0 = _S391;
            (&_S452)->differential_0 = _S451;
            DiffPair_float_0 _S453 = s_fwd_length_impl_0(&_S452);
            DiffPair_float_0 _S454;
            (&_S454)->primal_0 = _S453.primal_0;
            (&_S454)->differential_0 = _S453.differential_0;
            DiffPair_float_0 _S455;
            (&_S455)->primal_0 = _S434;
            (&_S455)->differential_0 = 0.0f;
            DiffPair_float_0 _S456 = _d_atan2_1(&_S454, &_S455);
            if((_S456.primal_0) < 0.00100000004749745f)
            {
                float _S457 = _S456.differential_0 * _S456.primal_0;
                float _S458 = (0.0f - (_S457 + _S457) * 0.3333333432674408f) * _S392 / (_S392 * _S392);
                k_2 = (1.0f - _S456.primal_0 * _S456.primal_0 / 3.0f) / _S392;
                eps2d_1 = _S458;
            }
            else
            {
                float _S459 = (_S456.differential_0 * _S453.primal_0 - _S456.primal_0 * _S453.differential_0) / (_S453.primal_0 * _S453.primal_0);
                k_2 = _S456.primal_0 / _S453.primal_0;
                eps2d_1 = _S459;
            }
            float2  _S460 = _S391 * make_float2 (k_2);
            float2  _S461 = _S451 * make_float2 (k_2) + make_float2 (eps2d_1) * _S391;
            float u_23 = _S460.x;
            float s_diff_u_13 = _S461.x;
            float v_23 = _S460.y;
            float s_diff_v_13 = _S461.y;
            float _S462 = s_diff_u_13 * u_23;
            float _S463 = s_diff_v_13 * v_23;
            float r2_23 = u_23 * u_23 + v_23 * v_23;
            float s_diff_r2_13 = _S462 + _S462 + (_S463 + _S463);
            float _S464 = _S395 + r2_23 * _S396;
            float _S465 = _S394 + r2_23 * _S464;
            float _S466 = _S393 + r2_23 * _S465;
            float2  _S467 = _S461 * make_float2 (1.0f + r2_23 * _S466) + make_float2 (s_diff_r2_13 * _S466 + (s_diff_r2_13 * _S465 + (s_diff_r2_13 * _S464 + s_diff_r2_13 * _S396 * r2_23) * r2_23) * r2_23) * _S460 + make_float2 (s_diff_u_13 * _S403 * v_23 + s_diff_v_13 * (_S403 * u_23) + (s_diff_r2_13 + (s_diff_u_13 * 2.0f * u_23 + s_diff_u_13 * (2.0f * u_23))) * _S398 + s_diff_r2_13 * _S399, s_diff_u_13 * _S404 * v_23 + s_diff_v_13 * (_S404 * u_23) + (s_diff_r2_13 + (s_diff_v_13 * 2.0f * v_23 + s_diff_v_13 * (2.0f * v_23))) * _S397 + s_diff_r2_13 * _S400);
            float2  _S468 = _S467 + make_float2 (_S467.x * _S401 + _S467.y * _S402, 0.0f);
            float _S469 = _S468.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->y) = _S468.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->y) = _S469;
            float2  _S470 = make_float2 (0.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S471;
            (&_S471)->primal_0 = _S391;
            (&_S471)->differential_0 = _S470;
            DiffPair_float_0 _S472 = s_fwd_length_impl_0(&_S471);
            DiffPair_float_0 _S473;
            (&_S473)->primal_0 = _S472.primal_0;
            (&_S473)->differential_0 = _S472.differential_0;
            DiffPair_float_0 _S474;
            (&_S474)->primal_0 = _S392;
            (&_S474)->differential_0 = 1.0f;
            DiffPair_float_0 _S475 = _d_atan2_1(&_S473, &_S474);
            if((_S475.primal_0) < 0.00100000004749745f)
            {
                float _S476 = _S475.differential_0 * _S475.primal_0;
                float _S477 = 1.0f - _S475.primal_0 * _S475.primal_0 / 3.0f;
                float _S478 = ((0.0f - (_S476 + _S476) * 0.3333333432674408f) * _S392 - _S477) / (_S392 * _S392);
                k_2 = _S477 / _S392;
                eps2d_1 = _S478;
            }
            else
            {
                float _S479 = (_S475.differential_0 * _S472.primal_0 - _S475.primal_0 * _S472.differential_0) / (_S472.primal_0 * _S472.primal_0);
                k_2 = _S475.primal_0 / _S472.primal_0;
                eps2d_1 = _S479;
            }
            float2  _S480 = _S391 * make_float2 (k_2);
            float2  _S481 = make_float2 (eps2d_1) * _S391;
            float u_24 = _S480.x;
            float s_diff_u_14 = _S481.x;
            float v_24 = _S480.y;
            float s_diff_v_14 = _S481.y;
            float _S482 = s_diff_u_14 * u_24;
            float _S483 = s_diff_v_14 * v_24;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float s_diff_r2_14 = _S482 + _S482 + (_S483 + _S483);
            float _S484 = _S395 + r2_24 * _S396;
            float _S485 = _S394 + r2_24 * _S484;
            float _S486 = _S393 + r2_24 * _S485;
            float2  _S487 = _S481 * make_float2 (1.0f + r2_24 * _S486) + make_float2 (s_diff_r2_14 * _S486 + (s_diff_r2_14 * _S485 + (s_diff_r2_14 * _S484 + s_diff_r2_14 * _S396 * r2_24) * r2_24) * r2_24) * _S480 + make_float2 (s_diff_u_14 * _S403 * v_24 + s_diff_v_14 * (_S403 * u_24) + (s_diff_r2_14 + (s_diff_u_14 * 2.0f * u_24 + s_diff_u_14 * (2.0f * u_24))) * _S398 + s_diff_r2_14 * _S399, s_diff_u_14 * _S404 * v_24 + s_diff_v_14 * (_S404 * u_24) + (s_diff_r2_14 + (s_diff_v_14 * 2.0f * v_24 + s_diff_v_14 * (2.0f * v_24))) * _S397 + s_diff_r2_14 * _S400);
            float2  _S488 = _S487 + make_float2 (_S487.x * _S401 + _S487.y * _S402, 0.0f);
            float _S489 = _S488.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->z) = _S488.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->z) = _S489;
            covar2d_1 = mul_4(mul_3(J_7, covar_c_1), transpose_1(J_7));
            is_valid_0 = true;
            break;
        }
        bool is_valid_1 = true & is_valid_0;
        float2  mean2d_c_0 = *mean2d_5 - make_float2 (cx_2, cy_2);
        float invdet_1 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        float opac_0 = *opacity_1 * (F32_exp((-0.5f * dot_0(mul_7(makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_1, - covar2d_1.rows[int(0)].y * invdet_1, - covar2d_1.rows[int(1)].x * invdet_1, covar2d_1.rows[int(0)].x * invdet_1), mean2d_c_0), mean2d_c_0))));
        if(_S392 < 0.0f)
        {
            is_valid_0 = opac_0 > 0.00392156885936856f;
        }
        else
        {
            is_valid_0 = false;
        }
        if(is_valid_0)
        {
            is_valid_0 = false;
        }
        else
        {
            is_valid_0 = is_valid_1;
        }
        if(!is_valid_0)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        if(antialiased_1)
        {
            eps2d_1 = 0.10000000149011612f;
        }
        else
        {
            eps2d_1 = 0.30000001192092896f;
        }
        float det_orig_1 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S490 = *&(((&covar2d_1)->rows + (int(0)))->x) + eps2d_1;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S490;
        float _S491 = *&(((&covar2d_1)->rows + (int(1)))->y) + eps2d_1;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S491;
        float det_blur_1 = _S490 * _S491 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S492 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
        if(antialiased_1)
        {
            *opacity_1 = *opacity_1 * compensation_1;
        }
        if((*opacity_1) < 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S493 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_1 / 0.00392156885936856f)))))))));
        float radius_x_1 = _S493 * (F32_sqrt((covar2d_1[int(0)].x)));
        float radius_y_1 = _S493 * (F32_sqrt((covar2d_1[int(1)].y)));
        float _S494 = (*mean2d_5).x - radius_x_1;
        float _S495 = (*mean2d_5).x + radius_x_1;
        float _S496 = (*mean2d_5).y - radius_y_1;
        float _S497 = (*mean2d_5).y + radius_y_1;
        if(_S495 <= 0.0f)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S494 >= float(image_width_1);
        }
        if(is_valid_0)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S497 <= 0.0f;
        }
        if(is_valid_0)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S496 >= float(image_height_1);
        }
        if(is_valid_0)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (_S494, _S496, _S495, _S497);
        float x_17 = mean_c_1.x;
        float y_5 = mean_c_1.y;
        float _S498 = x_17 * x_17 + y_5 * y_5;
        *sorting_depth_1 = _S392 * _S392 * _S392 * _S392 + 0.001953125f * _S498 * _S498;
        *conic_1 = make_float3 (_S492.rows[int(0)].x, _S492.rows[int(0)].y, _S492.rows[int(1)].y);
        *radius_2 = view_radius_3dgs_0(mean_2, scale_1, in_opacity_1, - mul_6(_S409, t_1));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_equisolid(bool antialiased_2, float3  mean_3, float4  quat_2, float3  scale_2, float in_opacity_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_6, float fy_6, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_6, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * sorting_depth_2, float * radius_3, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2)
{
    float2  _S499;
    float _S500;
    float _S501;
    float _S502;
    float _S503;
    float _S504;
    float _S505;
    float _S506;
    float _S507;
    float _S508;
    float _S509;
    float _S510;
    float _S511;
    float _S512;
    bool _S513;
    for(;;)
    {
        float3  mean_c_2 = mul_6(R_2, mean_3) + t_2;
        float _S514 = length_1(mean_c_2);
        *depth_2 = _S514;
        if(_S514 <= 0.0f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        *opacity_2 = 1.0f / (1.0f + (F32_exp((- in_opacity_2))));
        bool is_valid_2;
        float eps2d_2;
        float4  _S515 = normalize_0(quat_2);
        float3  _S516 = exp_0(scale_2);
        float x_18 = _S515.y;
        float x2_2 = x_18 * x_18;
        float y2_2 = _S515.z * _S515.z;
        float z2_2 = _S515.w * _S515.w;
        float xy_2 = _S515.y * _S515.z;
        float xz_2 = _S515.y * _S515.w;
        float yz_2 = _S515.z * _S515.w;
        float wx_2 = _S515.x * _S515.y;
        float wy_2 = _S515.x * _S515.z;
        float wz_2 = _S515.x * _S515.w;
        Matrix<float, 3, 3>  M_2 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (_S516.x, 0.0f, 0.0f, 0.0f, _S516.y, 0.0f, 0.0f, 0.0f, _S516.z));
        Matrix<float, 3, 3>  _S517 = transpose_3(R_2);
        Matrix<float, 3, 3>  covar_c_2 = mul_5(mul_5(R_2, mul_5(M_2, transpose_3(M_2))), _S517);
        for(;;)
        {
            float k_3;
            for(;;)
            {
                float2  _S518 = float2 {mean_c_2.x, mean_c_2.y};
                _S499 = _S518;
                float r_10 = length_0(_S518);
                float _S519 = mean_c_2.z;
                _S500 = _S519;
                float theta_3 = (F32_atan2((r_10), (_S519)));
                if(r_10 < 9.99999997475242708e-07f)
                {
                    k_3 = (1.0f - theta_3 * theta_3 / 24.0f) / _S519;
                }
                else
                {
                    k_3 = 2.0f * (F32_sin((0.5f * theta_3))) / r_10;
                }
                float2  _S520 = _S518 * make_float2 (k_3);
                *mean2d_6 = _S520;
                float2  _S521 = make_float2 (1.0f, 0.0f);
                _S501 = dist_coeffs_6[int(0)];
                _S502 = dist_coeffs_6[int(1)];
                _S503 = dist_coeffs_6[int(2)];
                _S504 = dist_coeffs_6[int(3)];
                _S505 = dist_coeffs_6[int(4)];
                _S506 = dist_coeffs_6[int(5)];
                _S507 = dist_coeffs_6[int(6)];
                _S508 = dist_coeffs_6[int(7)];
                _S509 = dist_coeffs_6[int(8)];
                _S510 = dist_coeffs_6[int(9)];
                float u_25 = _S520.x;
                float v_25 = _S520.y;
                float _S522 = u_25 + u_25;
                float r2_25 = u_25 * u_25 + v_25 * v_25;
                float _S523 = dist_coeffs_6[int(2)] + r2_25 * dist_coeffs_6[int(3)];
                float _S524 = dist_coeffs_6[int(1)] + r2_25 * _S523;
                float _S525 = dist_coeffs_6[int(0)] + r2_25 * _S524;
                float _S526 = _S522 * _S525 + (_S522 * _S524 + (_S522 * _S523 + _S522 * dist_coeffs_6[int(3)] * r2_25) * r2_25) * r2_25;
                float radial_5 = 1.0f + r2_25 * _S525;
                float _S527 = 2.0f * dist_coeffs_6[int(4)];
                _S511 = _S527;
                float _S528 = _S527 * u_25;
                float _S529 = 2.0f * u_25;
                float s_diff_du_3 = _S527 * v_25 + (_S522 + (_S529 + _S529)) * dist_coeffs_6[int(5)] + _S522 * dist_coeffs_6[int(6)];
                float _S530 = 2.0f * dist_coeffs_6[int(5)];
                _S512 = _S530;
                float _S531 = 2.0f * v_25;
                float2  _S532 = _S521 * make_float2 (radial_5) + make_float2 (_S526) * _S520 + make_float2 (s_diff_du_3, _S530 * v_25 + _S522 * dist_coeffs_6[int(4)] + _S522 * dist_coeffs_6[int(7)]);
                float _S533 = v_25 + v_25;
                float2  _S534 = make_float2 (0.0f, 1.0f) * make_float2 (radial_5) + make_float2 (_S533 * _S525 + (_S533 * _S524 + (_S533 * _S523 + _S533 * dist_coeffs_6[int(3)] * r2_25) * r2_25) * r2_25) * _S520 + make_float2 (_S528 + _S533 * dist_coeffs_6[int(5)] + _S533 * dist_coeffs_6[int(6)], _S530 * u_25 + (_S533 + (_S531 + _S531)) * dist_coeffs_6[int(4)] + _S533 * dist_coeffs_6[int(7)]);
                Matrix<float, 2, 2>  _S535 = transpose_0(makeMatrix<float, 2, 2> (_S532 + make_float2 (_S532.x * dist_coeffs_6[int(8)] + _S532.y * dist_coeffs_6[int(9)], 0.0f), _S534 + make_float2 (_S534.x * dist_coeffs_6[int(8)] + _S534.y * dist_coeffs_6[int(9)], 0.0f)));
                bool _S536 = !((F32_min((determinant_0(_S535)), ((F32_min((_S535.rows[int(0)].x), (_S535.rows[int(1)].y)))))) > 0.0f);
                _S513 = _S536;
                if(_S536)
                {
                    break;
                }
                float u_26 = (*mean2d_6).x;
                float v_26 = (*mean2d_6).y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S537 = *mean2d_6 * make_float2 (1.0f + r2_26 * (dist_coeffs_6[int(0)] + r2_26 * (dist_coeffs_6[int(1)] + r2_26 * (dist_coeffs_6[int(2)] + r2_26 * dist_coeffs_6[int(3)])))) + make_float2 (_S527 * u_26 * v_26 + dist_coeffs_6[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + dist_coeffs_6[int(6)] * r2_26, _S530 * u_26 * v_26 + dist_coeffs_6[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + dist_coeffs_6[int(7)] * r2_26);
                float2  _S538 = _S537 + make_float2 (dist_coeffs_6[int(8)] * _S537.x + dist_coeffs_6[int(9)] * _S537.y, 0.0f);
                *mean2d_6 = make_float2 (fx_6 * _S538.x + cx_3, fy_6 * _S538.y + cy_3);
                break;
            }
            if(!!_S513)
            {
                is_valid_2 = false;
                break;
            }
            Matrix<float, 2, 3>  J_8;
            float2  _S539 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S540;
            (&_S540)->primal_0 = _S499;
            (&_S540)->differential_0 = _S539;
            DiffPair_float_0 _S541 = s_fwd_length_impl_0(&_S540);
            float _S542 = _S500;
            DiffPair_float_0 _S543;
            (&_S543)->primal_0 = _S541.primal_0;
            (&_S543)->differential_0 = _S541.differential_0;
            DiffPair_float_0 _S544;
            (&_S544)->primal_0 = _S500;
            (&_S544)->differential_0 = 0.0f;
            DiffPair_float_0 _S545 = _d_atan2_1(&_S543, &_S544);
            if((_S541.primal_0) < 9.99999997475242708e-07f)
            {
                float _S546 = _S545.differential_0 * _S545.primal_0;
                float _S547 = (0.0f - (_S546 + _S546) * 0.0416666679084301f) * _S500 / (_S500 * _S500);
                k_3 = (1.0f - _S545.primal_0 * _S545.primal_0 / 24.0f) / _S500;
                eps2d_2 = _S547;
            }
            else
            {
                float _S548 = _S545.differential_0 * 0.5f;
                DiffPair_float_0 _S549;
                (&_S549)->primal_0 = 0.5f * _S545.primal_0;
                (&_S549)->differential_0 = _S548;
                DiffPair_float_0 _S550 = _d_sin_1(&_S549);
                float _S551 = 2.0f * _S550.primal_0;
                float _S552 = (_S550.differential_0 * 2.0f * _S541.primal_0 - _S551 * _S541.differential_0) / (_S541.primal_0 * _S541.primal_0);
                k_3 = _S551 / _S541.primal_0;
                eps2d_2 = _S552;
            }
            float2  _S553 = _S499 * make_float2 (k_3);
            float2  _S554 = _S539 * make_float2 (k_3) + make_float2 (eps2d_2) * _S499;
            float u_27 = _S553.x;
            float s_diff_u_15 = _S554.x;
            float v_27 = _S553.y;
            float s_diff_v_15 = _S554.y;
            float _S555 = s_diff_u_15 * u_27;
            float _S556 = s_diff_v_15 * v_27;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float s_diff_r2_15 = _S555 + _S555 + (_S556 + _S556);
            float _S557 = _S503 + r2_27 * _S504;
            float _S558 = _S502 + r2_27 * _S557;
            float _S559 = _S501 + r2_27 * _S558;
            float2  _S560 = _S554 * make_float2 (1.0f + r2_27 * _S559) + make_float2 (s_diff_r2_15 * _S559 + (s_diff_r2_15 * _S558 + (s_diff_r2_15 * _S557 + s_diff_r2_15 * _S504 * r2_27) * r2_27) * r2_27) * _S553 + make_float2 (s_diff_u_15 * _S511 * v_27 + s_diff_v_15 * (_S511 * u_27) + (s_diff_r2_15 + (s_diff_u_15 * 2.0f * u_27 + s_diff_u_15 * (2.0f * u_27))) * _S506 + s_diff_r2_15 * _S507, s_diff_u_15 * _S512 * v_27 + s_diff_v_15 * (_S512 * u_27) + (s_diff_r2_15 + (s_diff_v_15 * 2.0f * v_27 + s_diff_v_15 * (2.0f * v_27))) * _S505 + s_diff_r2_15 * _S508);
            float2  _S561 = _S560 + make_float2 (_S560.x * _S509 + _S560.y * _S510, 0.0f);
            float _S562 = _S561.y * fy_6;
            *&(((&J_8)->rows + (int(0)))->x) = _S561.x * fx_6;
            *&(((&J_8)->rows + (int(1)))->x) = _S562;
            float2  _S563 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S564;
            (&_S564)->primal_0 = _S499;
            (&_S564)->differential_0 = _S563;
            DiffPair_float_0 _S565 = s_fwd_length_impl_0(&_S564);
            DiffPair_float_0 _S566;
            (&_S566)->primal_0 = _S565.primal_0;
            (&_S566)->differential_0 = _S565.differential_0;
            DiffPair_float_0 _S567;
            (&_S567)->primal_0 = _S542;
            (&_S567)->differential_0 = 0.0f;
            DiffPair_float_0 _S568 = _d_atan2_1(&_S566, &_S567);
            if((_S565.primal_0) < 9.99999997475242708e-07f)
            {
                float _S569 = _S568.differential_0 * _S568.primal_0;
                float _S570 = (0.0f - (_S569 + _S569) * 0.0416666679084301f) * _S500 / (_S500 * _S500);
                k_3 = (1.0f - _S568.primal_0 * _S568.primal_0 / 24.0f) / _S500;
                eps2d_2 = _S570;
            }
            else
            {
                float _S571 = _S568.differential_0 * 0.5f;
                DiffPair_float_0 _S572;
                (&_S572)->primal_0 = 0.5f * _S568.primal_0;
                (&_S572)->differential_0 = _S571;
                DiffPair_float_0 _S573 = _d_sin_1(&_S572);
                float _S574 = 2.0f * _S573.primal_0;
                float _S575 = (_S573.differential_0 * 2.0f * _S565.primal_0 - _S574 * _S565.differential_0) / (_S565.primal_0 * _S565.primal_0);
                k_3 = _S574 / _S565.primal_0;
                eps2d_2 = _S575;
            }
            float2  _S576 = _S499 * make_float2 (k_3);
            float2  _S577 = _S563 * make_float2 (k_3) + make_float2 (eps2d_2) * _S499;
            float u_28 = _S576.x;
            float s_diff_u_16 = _S577.x;
            float v_28 = _S576.y;
            float s_diff_v_16 = _S577.y;
            float _S578 = s_diff_u_16 * u_28;
            float _S579 = s_diff_v_16 * v_28;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float s_diff_r2_16 = _S578 + _S578 + (_S579 + _S579);
            float _S580 = _S503 + r2_28 * _S504;
            float _S581 = _S502 + r2_28 * _S580;
            float _S582 = _S501 + r2_28 * _S581;
            float2  _S583 = _S577 * make_float2 (1.0f + r2_28 * _S582) + make_float2 (s_diff_r2_16 * _S582 + (s_diff_r2_16 * _S581 + (s_diff_r2_16 * _S580 + s_diff_r2_16 * _S504 * r2_28) * r2_28) * r2_28) * _S576 + make_float2 (s_diff_u_16 * _S511 * v_28 + s_diff_v_16 * (_S511 * u_28) + (s_diff_r2_16 + (s_diff_u_16 * 2.0f * u_28 + s_diff_u_16 * (2.0f * u_28))) * _S506 + s_diff_r2_16 * _S507, s_diff_u_16 * _S512 * v_28 + s_diff_v_16 * (_S512 * u_28) + (s_diff_r2_16 + (s_diff_v_16 * 2.0f * v_28 + s_diff_v_16 * (2.0f * v_28))) * _S505 + s_diff_r2_16 * _S508);
            float2  _S584 = _S583 + make_float2 (_S583.x * _S509 + _S583.y * _S510, 0.0f);
            float _S585 = _S584.y * fy_6;
            *&(((&J_8)->rows + (int(0)))->y) = _S584.x * fx_6;
            *&(((&J_8)->rows + (int(1)))->y) = _S585;
            float2  _S586 = make_float2 (0.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S587;
            (&_S587)->primal_0 = _S499;
            (&_S587)->differential_0 = _S586;
            DiffPair_float_0 _S588 = s_fwd_length_impl_0(&_S587);
            DiffPair_float_0 _S589;
            (&_S589)->primal_0 = _S588.primal_0;
            (&_S589)->differential_0 = _S588.differential_0;
            DiffPair_float_0 _S590;
            (&_S590)->primal_0 = _S500;
            (&_S590)->differential_0 = 1.0f;
            DiffPair_float_0 _S591 = _d_atan2_1(&_S589, &_S590);
            if((_S588.primal_0) < 9.99999997475242708e-07f)
            {
                float _S592 = _S591.differential_0 * _S591.primal_0;
                float _S593 = 1.0f - _S591.primal_0 * _S591.primal_0 / 24.0f;
                float _S594 = ((0.0f - (_S592 + _S592) * 0.0416666679084301f) * _S500 - _S593) / (_S500 * _S500);
                k_3 = _S593 / _S500;
                eps2d_2 = _S594;
            }
            else
            {
                float _S595 = _S591.differential_0 * 0.5f;
                DiffPair_float_0 _S596;
                (&_S596)->primal_0 = 0.5f * _S591.primal_0;
                (&_S596)->differential_0 = _S595;
                DiffPair_float_0 _S597 = _d_sin_1(&_S596);
                float _S598 = 2.0f * _S597.primal_0;
                float _S599 = (_S597.differential_0 * 2.0f * _S588.primal_0 - _S598 * _S588.differential_0) / (_S588.primal_0 * _S588.primal_0);
                k_3 = _S598 / _S588.primal_0;
                eps2d_2 = _S599;
            }
            float2  _S600 = _S499 * make_float2 (k_3);
            float2  _S601 = make_float2 (eps2d_2) * _S499;
            float u_29 = _S600.x;
            float s_diff_u_17 = _S601.x;
            float v_29 = _S600.y;
            float s_diff_v_17 = _S601.y;
            float _S602 = s_diff_u_17 * u_29;
            float _S603 = s_diff_v_17 * v_29;
            float r2_29 = u_29 * u_29 + v_29 * v_29;
            float s_diff_r2_17 = _S602 + _S602 + (_S603 + _S603);
            float _S604 = _S503 + r2_29 * _S504;
            float _S605 = _S502 + r2_29 * _S604;
            float _S606 = _S501 + r2_29 * _S605;
            float2  _S607 = _S601 * make_float2 (1.0f + r2_29 * _S606) + make_float2 (s_diff_r2_17 * _S606 + (s_diff_r2_17 * _S605 + (s_diff_r2_17 * _S604 + s_diff_r2_17 * _S504 * r2_29) * r2_29) * r2_29) * _S600 + make_float2 (s_diff_u_17 * _S511 * v_29 + s_diff_v_17 * (_S511 * u_29) + (s_diff_r2_17 + (s_diff_u_17 * 2.0f * u_29 + s_diff_u_17 * (2.0f * u_29))) * _S506 + s_diff_r2_17 * _S507, s_diff_u_17 * _S512 * v_29 + s_diff_v_17 * (_S512 * u_29) + (s_diff_r2_17 + (s_diff_v_17 * 2.0f * v_29 + s_diff_v_17 * (2.0f * v_29))) * _S505 + s_diff_r2_17 * _S508);
            float2  _S608 = _S607 + make_float2 (_S607.x * _S509 + _S607.y * _S510, 0.0f);
            float _S609 = _S608.y * fy_6;
            *&(((&J_8)->rows + (int(0)))->z) = _S608.x * fx_6;
            *&(((&J_8)->rows + (int(1)))->z) = _S609;
            covar2d_2 = mul_4(mul_3(J_8, covar_c_2), transpose_1(J_8));
            is_valid_2 = true;
            break;
        }
        bool is_valid_3 = true & is_valid_2;
        float2  mean2d_c_1 = *mean2d_6 - make_float2 (cx_3, cy_3);
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        float opac_1 = *opacity_2 * (F32_exp((-0.5f * dot_0(mul_7(makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3), mean2d_c_1), mean2d_c_1))));
        if(_S500 < 0.0f)
        {
            is_valid_2 = opac_1 > 0.00392156885936856f;
        }
        else
        {
            is_valid_2 = false;
        }
        if(is_valid_2)
        {
            is_valid_2 = false;
        }
        else
        {
            is_valid_2 = is_valid_3;
        }
        if(!is_valid_2)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        if(antialiased_2)
        {
            eps2d_2 = 0.10000000149011612f;
        }
        else
        {
            eps2d_2 = 0.30000001192092896f;
        }
        float det_orig_2 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S610 = *&(((&covar2d_2)->rows + (int(0)))->x) + eps2d_2;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S610;
        float _S611 = *&(((&covar2d_2)->rows + (int(1)))->y) + eps2d_2;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S611;
        float det_blur_2 = _S610 * _S611 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        float invdet_4 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S612 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_4, - covar2d_2.rows[int(0)].y * invdet_4, - covar2d_2.rows[int(1)].x * invdet_4, covar2d_2.rows[int(0)].x * invdet_4);
        if(antialiased_2)
        {
            *opacity_2 = *opacity_2 * compensation_2;
        }
        if((*opacity_2) < 0.00392156885936856f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        float _S613 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_2 / 0.00392156885936856f)))))))));
        float radius_x_2 = _S613 * (F32_sqrt((covar2d_2[int(0)].x)));
        float radius_y_2 = _S613 * (F32_sqrt((covar2d_2[int(1)].y)));
        float _S614 = (*mean2d_6).x - radius_x_2;
        float _S615 = (*mean2d_6).x + radius_x_2;
        float _S616 = (*mean2d_6).y - radius_y_2;
        float _S617 = (*mean2d_6).y + radius_y_2;
        if(_S615 <= 0.0f)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S614 >= float(image_width_2);
        }
        if(is_valid_2)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S617 <= 0.0f;
        }
        if(is_valid_2)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S616 >= float(image_height_2);
        }
        if(is_valid_2)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_2 = make_float4 (_S614, _S616, _S615, _S617);
        float x_19 = mean_c_2.x;
        float y_6 = mean_c_2.y;
        float _S618 = x_19 * x_19 + y_6 * y_6;
        *sorting_depth_2 = _S500 * _S500 * _S500 * _S500 + 0.001953125f * _S618 * _S618;
        *conic_2 = make_float3 (_S612.rows[int(0)].x, _S612.rows[int(0)].y, _S612.rows[int(1)].y);
        *radius_3 = view_radius_3dgs_0(mean_3, scale_2, in_opacity_2, - mul_6(_S517, t_2));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_equirect(bool antialiased_3, float3  mean_4, float4  quat_3, float3  scale_3, float in_opacity_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_7, float fy_7, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_7, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * sorting_depth_3, float * radius_4, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_6(R_3, mean_4) + t_3;
        float _S619 = length_1(mean_c_3);
        *depth_3 = _S619;
        if(_S619 <= 0.0f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_3;
        *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
        float4  _S620 = normalize_0(quat_3);
        float3  _S621 = exp_0(scale_3);
        float x_20 = _S620.y;
        float x2_3 = x_20 * x_20;
        float y2_3 = _S620.z * _S620.z;
        float z2_3 = _S620.w * _S620.w;
        float xy_3 = _S620.y * _S620.z;
        float xz_3 = _S620.y * _S620.w;
        float yz_3 = _S620.z * _S620.w;
        float wx_3 = _S620.x * _S620.y;
        float wy_3 = _S620.x * _S620.z;
        float wz_3 = _S620.x * _S620.w;
        Matrix<float, 3, 3>  M_3 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S621.x, 0.0f, 0.0f, 0.0f, _S621.y, 0.0f, 0.0f, 0.0f, _S621.z));
        Matrix<float, 3, 3>  _S622 = transpose_3(R_3);
        Matrix<float, 3, 3>  covar_c_3 = mul_5(mul_5(R_3, mul_5(M_3, transpose_3(M_3))), _S622);
        float _S623 = mean_c_3.x;
        float _S624 = mean_c_3.z;
        float _S625 = mean_c_3.y;
        float2  _S626 = float2 {mean_c_3.x, mean_c_3.z};
        *mean2d_7 = make_float2 (fx_7 * (F32_atan2((_S623), (_S624))) + cx_4, fy_7 * (F32_atan2((_S625), (length_0(_S626)))) + cy_4);
        DiffPair_float_0 _S627;
        (&_S627)->primal_0 = _S623;
        (&_S627)->differential_0 = 1.0f;
        DiffPair_float_0 _S628;
        (&_S628)->primal_0 = _S624;
        (&_S628)->differential_0 = 0.0f;
        DiffPair_float_0 _S629 = _d_atan2_1(&_S627, &_S628);
        float2  _S630 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S631;
        (&_S631)->primal_0 = _S626;
        (&_S631)->differential_0 = _S630;
        DiffPair_float_0 _S632 = s_fwd_length_impl_0(&_S631);
        DiffPair_float_0 _S633;
        (&_S633)->primal_0 = _S625;
        (&_S633)->differential_0 = 0.0f;
        DiffPair_float_0 _S634;
        (&_S634)->primal_0 = _S632.primal_0;
        (&_S634)->differential_0 = _S632.differential_0;
        DiffPair_float_0 _S635 = _d_atan2_1(&_S633, &_S634);
        float _S636 = _S635.differential_0 * fy_7;
        Matrix<float, 2, 3>  J_9;
        *&(((&J_9)->rows + (int(0)))->x) = _S629.differential_0 * fx_7;
        *&(((&J_9)->rows + (int(1)))->x) = _S636;
        DiffPair_float_0 _S637;
        (&_S637)->primal_0 = _S623;
        (&_S637)->differential_0 = 0.0f;
        DiffPair_float_0 _S638;
        (&_S638)->primal_0 = _S624;
        (&_S638)->differential_0 = 0.0f;
        DiffPair_float_0 _S639 = _d_atan2_1(&_S637, &_S638);
        float2  _S640 = make_float2 (0.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S641;
        (&_S641)->primal_0 = _S626;
        (&_S641)->differential_0 = _S640;
        DiffPair_float_0 _S642 = s_fwd_length_impl_0(&_S641);
        DiffPair_float_0 _S643;
        (&_S643)->primal_0 = _S625;
        (&_S643)->differential_0 = 1.0f;
        DiffPair_float_0 _S644;
        (&_S644)->primal_0 = _S642.primal_0;
        (&_S644)->differential_0 = _S642.differential_0;
        DiffPair_float_0 _S645 = _d_atan2_1(&_S643, &_S644);
        float _S646 = _S645.differential_0 * fy_7;
        *&(((&J_9)->rows + (int(0)))->y) = _S639.differential_0 * fx_7;
        *&(((&J_9)->rows + (int(1)))->y) = _S646;
        DiffPair_float_0 _S647;
        (&_S647)->primal_0 = _S623;
        (&_S647)->differential_0 = 0.0f;
        DiffPair_float_0 _S648;
        (&_S648)->primal_0 = _S624;
        (&_S648)->differential_0 = 1.0f;
        DiffPair_float_0 _S649 = _d_atan2_1(&_S647, &_S648);
        float2  _S650 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S651;
        (&_S651)->primal_0 = _S626;
        (&_S651)->differential_0 = _S650;
        DiffPair_float_0 _S652 = s_fwd_length_impl_0(&_S651);
        DiffPair_float_0 _S653;
        (&_S653)->primal_0 = _S625;
        (&_S653)->differential_0 = 0.0f;
        DiffPair_float_0 _S654;
        (&_S654)->primal_0 = _S652.primal_0;
        (&_S654)->differential_0 = _S652.differential_0;
        DiffPair_float_0 _S655 = _d_atan2_1(&_S653, &_S654);
        float _S656 = _S655.differential_0 * fy_7;
        *&(((&J_9)->rows + (int(0)))->z) = _S649.differential_0 * fx_7;
        *&(((&J_9)->rows + (int(1)))->z) = _S656;
        covar2d_3 = mul_4(mul_3(J_9, covar_c_3), transpose_1(J_9));
        float eps2d_3;
        if(antialiased_3)
        {
            eps2d_3 = 0.10000000149011612f;
        }
        else
        {
            eps2d_3 = 0.30000001192092896f;
        }
        float det_orig_3 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S657 = *&(((&covar2d_3)->rows + (int(0)))->x) + eps2d_3;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S657;
        float _S658 = *&(((&covar2d_3)->rows + (int(1)))->y) + eps2d_3;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S658;
        float det_blur_3 = _S657 * _S658 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        float invdet_5 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
        Matrix<float, 2, 2>  _S659 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_5, - covar2d_3.rows[int(0)].y * invdet_5, - covar2d_3.rows[int(1)].x * invdet_5, covar2d_3.rows[int(0)].x * invdet_5);
        if(antialiased_3)
        {
            *opacity_3 = *opacity_3 * compensation_3;
        }
        if((*opacity_3) < 0.00392156885936856f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        float _S660 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
        float radius_x_3 = _S660 * (F32_sqrt((covar2d_3[int(0)].x)));
        float radius_y_3 = _S660 * (F32_sqrt((covar2d_3[int(1)].y)));
        float _S661 = (*mean2d_7).x - radius_x_3;
        float _S662 = (*mean2d_7).x + radius_x_3;
        float _S663 = (*mean2d_7).y - radius_y_3;
        float _S664 = (*mean2d_7).y + radius_y_3;
        bool _S665;
        if(_S662 <= 0.0f)
        {
            _S665 = true;
        }
        else
        {
            _S665 = _S661 >= float(image_width_3);
        }
        if(_S665)
        {
            _S665 = true;
        }
        else
        {
            _S665 = _S664 <= 0.0f;
        }
        if(_S665)
        {
            _S665 = true;
        }
        else
        {
            _S665 = _S663 >= float(image_height_3);
        }
        if(_S665)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_3 = make_float4 (_S661, _S663, _S662, _S664);
        *sorting_depth_3 = dot_1(mean_c_3, mean_c_3);
        *conic_3 = make_float3 (_S659.rows[int(0)].x, _S659.rows[int(0)].y, _S659.rows[int(1)].y);
        *radius_4 = view_radius_3dgs_0(mean_4, scale_3, in_opacity_3, - mul_6(_S622, t_3));
        break;
    }
    return;
}

struct SigmaPoints_0
{
    FixedArray<float3 , 7>  p_0;
    FixedArray<float, 7>  w_mean_0;
    FixedArray<float, 7>  w_cov_0;
};

inline __device__ void projection_3dgut_persp(bool antialiased_4, float3  mean_5, float4  quat_4, float3  scale_4, float in_opacity_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_8, float fy_8, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_8, uint image_width_4, uint image_height_4, float4  * aabb_xyxy_4, float * sorting_depth_4, float * radius_5, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4)
{
    float _S666;
    float _S667;
    float2  * _S668;
    float2  * _S669;
    float2  * _S670;
    bool _S671;
    float2  * _S672;
    float2  * _S673;
    float2  * _S674;
    bool _S675;
    float2  * _S676;
    bool _S677;
    for(;;)
    {
        float3  mean_c_4 = mul_6(R_4, mean_5) + t_4;
        float _S678 = mean_c_4.z;
        *depth_4 = length_1(mean_c_4);
        if(_S678 <= 0.0f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_4;
        *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
        bool _S679;
        float3  _S680 = exp_0(scale_4);
        float4  _S681 = normalize_0(quat_4);
        float x_21 = _S681.y;
        float x2_4 = x_21 * x_21;
        float y2_4 = _S681.z * _S681.z;
        float z2_4 = _S681.w * _S681.w;
        float xy_4 = _S681.y * _S681.z;
        float xz_4 = _S681.y * _S681.w;
        float yz_4 = _S681.z * _S681.w;
        float wx_4 = _S681.x * _S681.y;
        float wy_4 = _S681.x * _S681.z;
        float wz_4 = _S681.x * _S681.w;
        Matrix<float, 3, 3>  _S682 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_5;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S683 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S683 * _S680.x) * _S682.rows[0U];
        float3  _S684 = mean_5 + delta_0;
        float3  _S685 = mean_5 - delta_0;
        float3  delta_1 = make_float3 (_S683 * _S680.y) * _S682.rows[1U];
        float3  _S686 = mean_5 + delta_1;
        float3  _S687 = mean_5 - delta_1;
        float3  delta_2 = make_float3 (_S683 * _S680.z) * _S682.rows[2U];
        float3  _S688 = mean_5 + delta_2;
        float3  _S689 = mean_5 - delta_2;
        (&ret_0)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_0)->p_0[0U] = mul_6(R_4, (&ret_0)->p_0[0U]) + t_4;
        (&ret_0)->p_0[1U] = mul_6(R_4, _S684) + t_4;
        (&ret_0)->p_0[2U] = mul_6(R_4, _S686) + t_4;
        (&ret_0)->p_0[3U] = mul_6(R_4, _S688) + t_4;
        (&ret_0)->p_0[4U] = mul_6(R_4, _S685) + t_4;
        (&ret_0)->p_0[5U] = mul_6(R_4, _S687) + t_4;
        (&ret_0)->p_0[6U] = mul_6(R_4, _S689) + t_4;
        SigmaPoints_0 _S690 = ret_0;
        for(;;)
        {
            int2  _S691 = make_int2 (int(0));
            float2  _S692 = make_float2 ((float)_S691.x, (float)_S691.y);
            *mean2d_8 = _S692;
            covar2d_4 = makeMatrix<float, 2, 2> (0.0f);
            float _S693 = float(image_width_4);
            _S666 = _S693;
            float tan_fovx_0 = 0.5f * _S693 / fx_8;
            float _S694 = float(image_height_4);
            _S667 = _S694;
            float _S695 = 0.30000001192092896f * tan_fovx_0 * fx_8;
            float lim_x_pos_2 = _S693 + _S695;
            float _S696 = 0.30000001192092896f * (0.5f * _S694 / fy_8) * fy_8;
            float lim_y_pos_0 = _S694 + _S696;
            FixedArray<float2 , 7>  proj_points_0;
            for(;;)
            {
                _S668 = &proj_points_0[int(0)];
                for(;;)
                {
                    float _S697 = _S690.p_0[int(0)].z;
                    proj_points_0[int(0)] = float2 {_S690.p_0[int(0)].x, _S690.p_0[int(0)].y} / make_float2 (_S697);
                    if(_S697 < 0.0f)
                    {
                        _S679 = true;
                    }
                    else
                    {
                        float u_30 = proj_points_0[int(0)].x;
                        float v_30 = proj_points_0[int(0)].y;
                        float _S698 = u_30 + u_30;
                        float r2_30 = u_30 * u_30 + v_30 * v_30;
                        float _S699 = dist_coeffs_8[int(2)] + r2_30 * dist_coeffs_8[int(3)];
                        float _S700 = dist_coeffs_8[int(1)] + r2_30 * _S699;
                        float _S701 = dist_coeffs_8[int(0)] + r2_30 * _S700;
                        float radial_6 = 1.0f + r2_30 * _S701;
                        float _S702 = 2.0f * dist_coeffs_8[int(4)];
                        float _S703 = 2.0f * u_30;
                        float _S704 = 2.0f * dist_coeffs_8[int(5)];
                        float _S705 = 2.0f * v_30;
                        float2  _S706 = make_float2 (1.0f, 0.0f) * make_float2 (radial_6) + make_float2 (_S698 * _S701 + (_S698 * _S700 + (_S698 * _S699 + _S698 * dist_coeffs_8[int(3)] * r2_30) * r2_30) * r2_30) * proj_points_0[int(0)] + make_float2 (_S702 * v_30 + (_S698 + (_S703 + _S703)) * dist_coeffs_8[int(5)] + _S698 * dist_coeffs_8[int(6)], _S704 * v_30 + _S698 * dist_coeffs_8[int(4)] + _S698 * dist_coeffs_8[int(7)]);
                        float _S707 = v_30 + v_30;
                        float2  _S708 = make_float2 (0.0f, 1.0f) * make_float2 (radial_6) + make_float2 (_S707 * _S701 + (_S707 * _S700 + (_S707 * _S699 + _S707 * dist_coeffs_8[int(3)] * r2_30) * r2_30) * r2_30) * proj_points_0[int(0)] + make_float2 (_S702 * u_30 + _S707 * dist_coeffs_8[int(5)] + _S707 * dist_coeffs_8[int(6)], _S704 * u_30 + (_S707 + (_S705 + _S705)) * dist_coeffs_8[int(4)] + _S707 * dist_coeffs_8[int(7)]);
                        Matrix<float, 2, 2>  _S709 = transpose_0(makeMatrix<float, 2, 2> (_S706 + make_float2 (_S706.x * dist_coeffs_8[int(8)] + _S706.y * dist_coeffs_8[int(9)], 0.0f), _S708 + make_float2 (_S708.x * dist_coeffs_8[int(8)] + _S708.y * dist_coeffs_8[int(9)], 0.0f)));
                        _S679 = !((F32_min((determinant_0(_S709)), ((F32_min((_S709.rows[int(0)].x), (_S709.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S679)
                    {
                        break;
                    }
                    float u_31 = proj_points_0[int(0)].x;
                    float v_31 = proj_points_0[int(0)].y;
                    float r2_31 = u_31 * u_31 + v_31 * v_31;
                    float2  _S710 = proj_points_0[int(0)] * make_float2 (1.0f + r2_31 * (dist_coeffs_8[int(0)] + r2_31 * (dist_coeffs_8[int(1)] + r2_31 * (dist_coeffs_8[int(2)] + r2_31 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_31 * v_31 + dist_coeffs_8[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + dist_coeffs_8[int(6)] * r2_31, 2.0f * dist_coeffs_8[int(5)] * u_31 * v_31 + dist_coeffs_8[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + dist_coeffs_8[int(7)] * r2_31);
                    float2  _S711 = _S710 + make_float2 (dist_coeffs_8[int(8)] * _S710.x + dist_coeffs_8[int(9)] * _S710.y, 0.0f);
                    proj_points_0[int(0)] = make_float2 (fx_8 * _S711.x + cx_5, fy_8 * _S711.y + cy_5);
                    break;
                }
                bool all_valid_0 = true & (!_S679);
                _S669 = &proj_points_0[int(1)];
                for(;;)
                {
                    float _S712 = _S690.p_0[int(1)].z;
                    proj_points_0[int(1)] = float2 {_S690.p_0[int(1)].x, _S690.p_0[int(1)].y} / make_float2 (_S712);
                    if(_S712 < 0.0f)
                    {
                        _S679 = true;
                    }
                    else
                    {
                        float u_32 = proj_points_0[int(1)].x;
                        float v_32 = proj_points_0[int(1)].y;
                        float _S713 = u_32 + u_32;
                        float r2_32 = u_32 * u_32 + v_32 * v_32;
                        float _S714 = dist_coeffs_8[int(2)] + r2_32 * dist_coeffs_8[int(3)];
                        float _S715 = dist_coeffs_8[int(1)] + r2_32 * _S714;
                        float _S716 = dist_coeffs_8[int(0)] + r2_32 * _S715;
                        float radial_7 = 1.0f + r2_32 * _S716;
                        float _S717 = 2.0f * dist_coeffs_8[int(4)];
                        float _S718 = 2.0f * u_32;
                        float _S719 = 2.0f * dist_coeffs_8[int(5)];
                        float _S720 = 2.0f * v_32;
                        float2  _S721 = make_float2 (1.0f, 0.0f) * make_float2 (radial_7) + make_float2 (_S713 * _S716 + (_S713 * _S715 + (_S713 * _S714 + _S713 * dist_coeffs_8[int(3)] * r2_32) * r2_32) * r2_32) * proj_points_0[int(1)] + make_float2 (_S717 * v_32 + (_S713 + (_S718 + _S718)) * dist_coeffs_8[int(5)] + _S713 * dist_coeffs_8[int(6)], _S719 * v_32 + _S713 * dist_coeffs_8[int(4)] + _S713 * dist_coeffs_8[int(7)]);
                        float _S722 = v_32 + v_32;
                        float2  _S723 = make_float2 (0.0f, 1.0f) * make_float2 (radial_7) + make_float2 (_S722 * _S716 + (_S722 * _S715 + (_S722 * _S714 + _S722 * dist_coeffs_8[int(3)] * r2_32) * r2_32) * r2_32) * proj_points_0[int(1)] + make_float2 (_S717 * u_32 + _S722 * dist_coeffs_8[int(5)] + _S722 * dist_coeffs_8[int(6)], _S719 * u_32 + (_S722 + (_S720 + _S720)) * dist_coeffs_8[int(4)] + _S722 * dist_coeffs_8[int(7)]);
                        Matrix<float, 2, 2>  _S724 = transpose_0(makeMatrix<float, 2, 2> (_S721 + make_float2 (_S721.x * dist_coeffs_8[int(8)] + _S721.y * dist_coeffs_8[int(9)], 0.0f), _S723 + make_float2 (_S723.x * dist_coeffs_8[int(8)] + _S723.y * dist_coeffs_8[int(9)], 0.0f)));
                        _S679 = !((F32_min((determinant_0(_S724)), ((F32_min((_S724.rows[int(0)].x), (_S724.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S679)
                    {
                        break;
                    }
                    float u_33 = proj_points_0[int(1)].x;
                    float v_33 = proj_points_0[int(1)].y;
                    float r2_33 = u_33 * u_33 + v_33 * v_33;
                    float2  _S725 = proj_points_0[int(1)] * make_float2 (1.0f + r2_33 * (dist_coeffs_8[int(0)] + r2_33 * (dist_coeffs_8[int(1)] + r2_33 * (dist_coeffs_8[int(2)] + r2_33 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_33 * v_33 + dist_coeffs_8[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + dist_coeffs_8[int(6)] * r2_33, 2.0f * dist_coeffs_8[int(5)] * u_33 * v_33 + dist_coeffs_8[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + dist_coeffs_8[int(7)] * r2_33);
                    float2  _S726 = _S725 + make_float2 (dist_coeffs_8[int(8)] * _S725.x + dist_coeffs_8[int(9)] * _S725.y, 0.0f);
                    proj_points_0[int(1)] = make_float2 (fx_8 * _S726.x + cx_5, fy_8 * _S726.y + cy_5);
                    break;
                }
                bool all_valid_1 = all_valid_0 & (!_S679);
                for(;;)
                {
                    _S670 = &proj_points_0[int(2)];
                    for(;;)
                    {
                        float _S727 = _S690.p_0[int(2)].z;
                        proj_points_0[int(2)] = float2 {_S690.p_0[int(2)].x, _S690.p_0[int(2)].y} / make_float2 (_S727);
                        if(_S727 < 0.0f)
                        {
                            _S679 = true;
                        }
                        else
                        {
                            float u_34 = proj_points_0[int(2)].x;
                            float v_34 = proj_points_0[int(2)].y;
                            float _S728 = u_34 + u_34;
                            float r2_34 = u_34 * u_34 + v_34 * v_34;
                            float _S729 = dist_coeffs_8[int(2)] + r2_34 * dist_coeffs_8[int(3)];
                            float _S730 = dist_coeffs_8[int(1)] + r2_34 * _S729;
                            float _S731 = dist_coeffs_8[int(0)] + r2_34 * _S730;
                            float radial_8 = 1.0f + r2_34 * _S731;
                            float _S732 = 2.0f * dist_coeffs_8[int(4)];
                            float _S733 = 2.0f * u_34;
                            float _S734 = 2.0f * dist_coeffs_8[int(5)];
                            float _S735 = 2.0f * v_34;
                            float2  _S736 = make_float2 (1.0f, 0.0f) * make_float2 (radial_8) + make_float2 (_S728 * _S731 + (_S728 * _S730 + (_S728 * _S729 + _S728 * dist_coeffs_8[int(3)] * r2_34) * r2_34) * r2_34) * proj_points_0[int(2)] + make_float2 (_S732 * v_34 + (_S728 + (_S733 + _S733)) * dist_coeffs_8[int(5)] + _S728 * dist_coeffs_8[int(6)], _S734 * v_34 + _S728 * dist_coeffs_8[int(4)] + _S728 * dist_coeffs_8[int(7)]);
                            float _S737 = v_34 + v_34;
                            float2  _S738 = make_float2 (0.0f, 1.0f) * make_float2 (radial_8) + make_float2 (_S737 * _S731 + (_S737 * _S730 + (_S737 * _S729 + _S737 * dist_coeffs_8[int(3)] * r2_34) * r2_34) * r2_34) * proj_points_0[int(2)] + make_float2 (_S732 * u_34 + _S737 * dist_coeffs_8[int(5)] + _S737 * dist_coeffs_8[int(6)], _S734 * u_34 + (_S737 + (_S735 + _S735)) * dist_coeffs_8[int(4)] + _S737 * dist_coeffs_8[int(7)]);
                            Matrix<float, 2, 2>  _S739 = transpose_0(makeMatrix<float, 2, 2> (_S736 + make_float2 (_S736.x * dist_coeffs_8[int(8)] + _S736.y * dist_coeffs_8[int(9)], 0.0f), _S738 + make_float2 (_S738.x * dist_coeffs_8[int(8)] + _S738.y * dist_coeffs_8[int(9)], 0.0f)));
                            _S679 = !((F32_min((determinant_0(_S739)), ((F32_min((_S739.rows[int(0)].x), (_S739.rows[int(1)].y)))))) > 0.0f);
                        }
                        if(_S679)
                        {
                            break;
                        }
                        float u_35 = proj_points_0[int(2)].x;
                        float v_35 = proj_points_0[int(2)].y;
                        float r2_35 = u_35 * u_35 + v_35 * v_35;
                        float2  _S740 = proj_points_0[int(2)] * make_float2 (1.0f + r2_35 * (dist_coeffs_8[int(0)] + r2_35 * (dist_coeffs_8[int(1)] + r2_35 * (dist_coeffs_8[int(2)] + r2_35 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_35 * v_35 + dist_coeffs_8[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + dist_coeffs_8[int(6)] * r2_35, 2.0f * dist_coeffs_8[int(5)] * u_35 * v_35 + dist_coeffs_8[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + dist_coeffs_8[int(7)] * r2_35);
                        float2  _S741 = _S740 + make_float2 (dist_coeffs_8[int(8)] * _S740.x + dist_coeffs_8[int(9)] * _S740.y, 0.0f);
                        proj_points_0[int(2)] = make_float2 (fx_8 * _S741.x + cx_5, fy_8 * _S741.y + cy_5);
                        break;
                    }
                    _S671 = all_valid_1 & (!_S679);
                    break;
                }
                _S672 = &proj_points_0[int(3)];
                for(;;)
                {
                    float _S742 = _S690.p_0[int(3)].z;
                    proj_points_0[int(3)] = float2 {_S690.p_0[int(3)].x, _S690.p_0[int(3)].y} / make_float2 (_S742);
                    if(_S742 < 0.0f)
                    {
                        _S679 = true;
                    }
                    else
                    {
                        float u_36 = proj_points_0[int(3)].x;
                        float v_36 = proj_points_0[int(3)].y;
                        float _S743 = u_36 + u_36;
                        float r2_36 = u_36 * u_36 + v_36 * v_36;
                        float _S744 = dist_coeffs_8[int(2)] + r2_36 * dist_coeffs_8[int(3)];
                        float _S745 = dist_coeffs_8[int(1)] + r2_36 * _S744;
                        float _S746 = dist_coeffs_8[int(0)] + r2_36 * _S745;
                        float radial_9 = 1.0f + r2_36 * _S746;
                        float _S747 = 2.0f * dist_coeffs_8[int(4)];
                        float _S748 = 2.0f * u_36;
                        float _S749 = 2.0f * dist_coeffs_8[int(5)];
                        float _S750 = 2.0f * v_36;
                        float2  _S751 = make_float2 (1.0f, 0.0f) * make_float2 (radial_9) + make_float2 (_S743 * _S746 + (_S743 * _S745 + (_S743 * _S744 + _S743 * dist_coeffs_8[int(3)] * r2_36) * r2_36) * r2_36) * proj_points_0[int(3)] + make_float2 (_S747 * v_36 + (_S743 + (_S748 + _S748)) * dist_coeffs_8[int(5)] + _S743 * dist_coeffs_8[int(6)], _S749 * v_36 + _S743 * dist_coeffs_8[int(4)] + _S743 * dist_coeffs_8[int(7)]);
                        float _S752 = v_36 + v_36;
                        float2  _S753 = make_float2 (0.0f, 1.0f) * make_float2 (radial_9) + make_float2 (_S752 * _S746 + (_S752 * _S745 + (_S752 * _S744 + _S752 * dist_coeffs_8[int(3)] * r2_36) * r2_36) * r2_36) * proj_points_0[int(3)] + make_float2 (_S747 * u_36 + _S752 * dist_coeffs_8[int(5)] + _S752 * dist_coeffs_8[int(6)], _S749 * u_36 + (_S752 + (_S750 + _S750)) * dist_coeffs_8[int(4)] + _S752 * dist_coeffs_8[int(7)]);
                        Matrix<float, 2, 2>  _S754 = transpose_0(makeMatrix<float, 2, 2> (_S751 + make_float2 (_S751.x * dist_coeffs_8[int(8)] + _S751.y * dist_coeffs_8[int(9)], 0.0f), _S753 + make_float2 (_S753.x * dist_coeffs_8[int(8)] + _S753.y * dist_coeffs_8[int(9)], 0.0f)));
                        _S679 = !((F32_min((determinant_0(_S754)), ((F32_min((_S754.rows[int(0)].x), (_S754.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S679)
                    {
                        break;
                    }
                    float u_37 = proj_points_0[int(3)].x;
                    float v_37 = proj_points_0[int(3)].y;
                    float r2_37 = u_37 * u_37 + v_37 * v_37;
                    float2  _S755 = proj_points_0[int(3)] * make_float2 (1.0f + r2_37 * (dist_coeffs_8[int(0)] + r2_37 * (dist_coeffs_8[int(1)] + r2_37 * (dist_coeffs_8[int(2)] + r2_37 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_37 * v_37 + dist_coeffs_8[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + dist_coeffs_8[int(6)] * r2_37, 2.0f * dist_coeffs_8[int(5)] * u_37 * v_37 + dist_coeffs_8[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + dist_coeffs_8[int(7)] * r2_37);
                    float2  _S756 = _S755 + make_float2 (dist_coeffs_8[int(8)] * _S755.x + dist_coeffs_8[int(9)] * _S755.y, 0.0f);
                    proj_points_0[int(3)] = make_float2 (fx_8 * _S756.x + cx_5, fy_8 * _S756.y + cy_5);
                    break;
                }
                bool all_valid_2 = _S671 & (!_S679);
                _S673 = &proj_points_0[int(4)];
                for(;;)
                {
                    float _S757 = _S690.p_0[int(4)].z;
                    proj_points_0[int(4)] = float2 {_S690.p_0[int(4)].x, _S690.p_0[int(4)].y} / make_float2 (_S757);
                    if(_S757 < 0.0f)
                    {
                        _S679 = true;
                    }
                    else
                    {
                        float u_38 = proj_points_0[int(4)].x;
                        float v_38 = proj_points_0[int(4)].y;
                        float _S758 = u_38 + u_38;
                        float r2_38 = u_38 * u_38 + v_38 * v_38;
                        float _S759 = dist_coeffs_8[int(2)] + r2_38 * dist_coeffs_8[int(3)];
                        float _S760 = dist_coeffs_8[int(1)] + r2_38 * _S759;
                        float _S761 = dist_coeffs_8[int(0)] + r2_38 * _S760;
                        float radial_10 = 1.0f + r2_38 * _S761;
                        float _S762 = 2.0f * dist_coeffs_8[int(4)];
                        float _S763 = 2.0f * u_38;
                        float _S764 = 2.0f * dist_coeffs_8[int(5)];
                        float _S765 = 2.0f * v_38;
                        float2  _S766 = make_float2 (1.0f, 0.0f) * make_float2 (radial_10) + make_float2 (_S758 * _S761 + (_S758 * _S760 + (_S758 * _S759 + _S758 * dist_coeffs_8[int(3)] * r2_38) * r2_38) * r2_38) * proj_points_0[int(4)] + make_float2 (_S762 * v_38 + (_S758 + (_S763 + _S763)) * dist_coeffs_8[int(5)] + _S758 * dist_coeffs_8[int(6)], _S764 * v_38 + _S758 * dist_coeffs_8[int(4)] + _S758 * dist_coeffs_8[int(7)]);
                        float _S767 = v_38 + v_38;
                        float2  _S768 = make_float2 (0.0f, 1.0f) * make_float2 (radial_10) + make_float2 (_S767 * _S761 + (_S767 * _S760 + (_S767 * _S759 + _S767 * dist_coeffs_8[int(3)] * r2_38) * r2_38) * r2_38) * proj_points_0[int(4)] + make_float2 (_S762 * u_38 + _S767 * dist_coeffs_8[int(5)] + _S767 * dist_coeffs_8[int(6)], _S764 * u_38 + (_S767 + (_S765 + _S765)) * dist_coeffs_8[int(4)] + _S767 * dist_coeffs_8[int(7)]);
                        Matrix<float, 2, 2>  _S769 = transpose_0(makeMatrix<float, 2, 2> (_S766 + make_float2 (_S766.x * dist_coeffs_8[int(8)] + _S766.y * dist_coeffs_8[int(9)], 0.0f), _S768 + make_float2 (_S768.x * dist_coeffs_8[int(8)] + _S768.y * dist_coeffs_8[int(9)], 0.0f)));
                        _S679 = !((F32_min((determinant_0(_S769)), ((F32_min((_S769.rows[int(0)].x), (_S769.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S679)
                    {
                        break;
                    }
                    float u_39 = proj_points_0[int(4)].x;
                    float v_39 = proj_points_0[int(4)].y;
                    float r2_39 = u_39 * u_39 + v_39 * v_39;
                    float2  _S770 = proj_points_0[int(4)] * make_float2 (1.0f + r2_39 * (dist_coeffs_8[int(0)] + r2_39 * (dist_coeffs_8[int(1)] + r2_39 * (dist_coeffs_8[int(2)] + r2_39 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_39 * v_39 + dist_coeffs_8[int(5)] * (r2_39 + 2.0f * u_39 * u_39) + dist_coeffs_8[int(6)] * r2_39, 2.0f * dist_coeffs_8[int(5)] * u_39 * v_39 + dist_coeffs_8[int(4)] * (r2_39 + 2.0f * v_39 * v_39) + dist_coeffs_8[int(7)] * r2_39);
                    float2  _S771 = _S770 + make_float2 (dist_coeffs_8[int(8)] * _S770.x + dist_coeffs_8[int(9)] * _S770.y, 0.0f);
                    proj_points_0[int(4)] = make_float2 (fx_8 * _S771.x + cx_5, fy_8 * _S771.y + cy_5);
                    break;
                }
                bool all_valid_3 = all_valid_2 & (!_S679);
                for(;;)
                {
                    _S674 = &proj_points_0[int(5)];
                    for(;;)
                    {
                        float _S772 = _S690.p_0[int(5)].z;
                        proj_points_0[int(5)] = float2 {_S690.p_0[int(5)].x, _S690.p_0[int(5)].y} / make_float2 (_S772);
                        if(_S772 < 0.0f)
                        {
                            _S679 = true;
                        }
                        else
                        {
                            float u_40 = proj_points_0[int(5)].x;
                            float v_40 = proj_points_0[int(5)].y;
                            float _S773 = u_40 + u_40;
                            float r2_40 = u_40 * u_40 + v_40 * v_40;
                            float _S774 = dist_coeffs_8[int(2)] + r2_40 * dist_coeffs_8[int(3)];
                            float _S775 = dist_coeffs_8[int(1)] + r2_40 * _S774;
                            float _S776 = dist_coeffs_8[int(0)] + r2_40 * _S775;
                            float radial_11 = 1.0f + r2_40 * _S776;
                            float _S777 = 2.0f * dist_coeffs_8[int(4)];
                            float _S778 = 2.0f * u_40;
                            float _S779 = 2.0f * dist_coeffs_8[int(5)];
                            float _S780 = 2.0f * v_40;
                            float2  _S781 = make_float2 (1.0f, 0.0f) * make_float2 (radial_11) + make_float2 (_S773 * _S776 + (_S773 * _S775 + (_S773 * _S774 + _S773 * dist_coeffs_8[int(3)] * r2_40) * r2_40) * r2_40) * proj_points_0[int(5)] + make_float2 (_S777 * v_40 + (_S773 + (_S778 + _S778)) * dist_coeffs_8[int(5)] + _S773 * dist_coeffs_8[int(6)], _S779 * v_40 + _S773 * dist_coeffs_8[int(4)] + _S773 * dist_coeffs_8[int(7)]);
                            float _S782 = v_40 + v_40;
                            float2  _S783 = make_float2 (0.0f, 1.0f) * make_float2 (radial_11) + make_float2 (_S782 * _S776 + (_S782 * _S775 + (_S782 * _S774 + _S782 * dist_coeffs_8[int(3)] * r2_40) * r2_40) * r2_40) * proj_points_0[int(5)] + make_float2 (_S777 * u_40 + _S782 * dist_coeffs_8[int(5)] + _S782 * dist_coeffs_8[int(6)], _S779 * u_40 + (_S782 + (_S780 + _S780)) * dist_coeffs_8[int(4)] + _S782 * dist_coeffs_8[int(7)]);
                            Matrix<float, 2, 2>  _S784 = transpose_0(makeMatrix<float, 2, 2> (_S781 + make_float2 (_S781.x * dist_coeffs_8[int(8)] + _S781.y * dist_coeffs_8[int(9)], 0.0f), _S783 + make_float2 (_S783.x * dist_coeffs_8[int(8)] + _S783.y * dist_coeffs_8[int(9)], 0.0f)));
                            _S679 = !((F32_min((determinant_0(_S784)), ((F32_min((_S784.rows[int(0)].x), (_S784.rows[int(1)].y)))))) > 0.0f);
                        }
                        if(_S679)
                        {
                            break;
                        }
                        float u_41 = proj_points_0[int(5)].x;
                        float v_41 = proj_points_0[int(5)].y;
                        float r2_41 = u_41 * u_41 + v_41 * v_41;
                        float2  _S785 = proj_points_0[int(5)] * make_float2 (1.0f + r2_41 * (dist_coeffs_8[int(0)] + r2_41 * (dist_coeffs_8[int(1)] + r2_41 * (dist_coeffs_8[int(2)] + r2_41 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_41 * v_41 + dist_coeffs_8[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + dist_coeffs_8[int(6)] * r2_41, 2.0f * dist_coeffs_8[int(5)] * u_41 * v_41 + dist_coeffs_8[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + dist_coeffs_8[int(7)] * r2_41);
                        float2  _S786 = _S785 + make_float2 (dist_coeffs_8[int(8)] * _S785.x + dist_coeffs_8[int(9)] * _S785.y, 0.0f);
                        proj_points_0[int(5)] = make_float2 (fx_8 * _S786.x + cx_5, fy_8 * _S786.y + cy_5);
                        break;
                    }
                    _S675 = all_valid_3 & (!_S679);
                    break;
                }
                _S676 = &proj_points_0[int(6)];
                for(;;)
                {
                    float _S787 = _S690.p_0[int(6)].z;
                    proj_points_0[int(6)] = float2 {_S690.p_0[int(6)].x, _S690.p_0[int(6)].y} / make_float2 (_S787);
                    if(_S787 < 0.0f)
                    {
                        _S679 = true;
                    }
                    else
                    {
                        float u_42 = proj_points_0[int(6)].x;
                        float v_42 = proj_points_0[int(6)].y;
                        float _S788 = u_42 + u_42;
                        float r2_42 = u_42 * u_42 + v_42 * v_42;
                        float _S789 = dist_coeffs_8[int(2)] + r2_42 * dist_coeffs_8[int(3)];
                        float _S790 = dist_coeffs_8[int(1)] + r2_42 * _S789;
                        float _S791 = dist_coeffs_8[int(0)] + r2_42 * _S790;
                        float radial_12 = 1.0f + r2_42 * _S791;
                        float _S792 = 2.0f * dist_coeffs_8[int(4)];
                        float _S793 = 2.0f * u_42;
                        float _S794 = 2.0f * dist_coeffs_8[int(5)];
                        float _S795 = 2.0f * v_42;
                        float2  _S796 = make_float2 (1.0f, 0.0f) * make_float2 (radial_12) + make_float2 (_S788 * _S791 + (_S788 * _S790 + (_S788 * _S789 + _S788 * dist_coeffs_8[int(3)] * r2_42) * r2_42) * r2_42) * proj_points_0[int(6)] + make_float2 (_S792 * v_42 + (_S788 + (_S793 + _S793)) * dist_coeffs_8[int(5)] + _S788 * dist_coeffs_8[int(6)], _S794 * v_42 + _S788 * dist_coeffs_8[int(4)] + _S788 * dist_coeffs_8[int(7)]);
                        float _S797 = v_42 + v_42;
                        float2  _S798 = make_float2 (0.0f, 1.0f) * make_float2 (radial_12) + make_float2 (_S797 * _S791 + (_S797 * _S790 + (_S797 * _S789 + _S797 * dist_coeffs_8[int(3)] * r2_42) * r2_42) * r2_42) * proj_points_0[int(6)] + make_float2 (_S792 * u_42 + _S797 * dist_coeffs_8[int(5)] + _S797 * dist_coeffs_8[int(6)], _S794 * u_42 + (_S797 + (_S795 + _S795)) * dist_coeffs_8[int(4)] + _S797 * dist_coeffs_8[int(7)]);
                        Matrix<float, 2, 2>  _S799 = transpose_0(makeMatrix<float, 2, 2> (_S796 + make_float2 (_S796.x * dist_coeffs_8[int(8)] + _S796.y * dist_coeffs_8[int(9)], 0.0f), _S798 + make_float2 (_S798.x * dist_coeffs_8[int(8)] + _S798.y * dist_coeffs_8[int(9)], 0.0f)));
                        _S679 = !((F32_min((determinant_0(_S799)), ((F32_min((_S799.rows[int(0)].x), (_S799.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S679)
                    {
                        break;
                    }
                    float u_43 = proj_points_0[int(6)].x;
                    float v_43 = proj_points_0[int(6)].y;
                    float r2_43 = u_43 * u_43 + v_43 * v_43;
                    float2  _S800 = proj_points_0[int(6)] * make_float2 (1.0f + r2_43 * (dist_coeffs_8[int(0)] + r2_43 * (dist_coeffs_8[int(1)] + r2_43 * (dist_coeffs_8[int(2)] + r2_43 * dist_coeffs_8[int(3)])))) + make_float2 (2.0f * dist_coeffs_8[int(4)] * u_43 * v_43 + dist_coeffs_8[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + dist_coeffs_8[int(6)] * r2_43, 2.0f * dist_coeffs_8[int(5)] * u_43 * v_43 + dist_coeffs_8[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + dist_coeffs_8[int(7)] * r2_43);
                    float2  _S801 = _S800 + make_float2 (dist_coeffs_8[int(8)] * _S800.x + dist_coeffs_8[int(9)] * _S800.y, 0.0f);
                    proj_points_0[int(6)] = make_float2 (fx_8 * _S801.x + cx_5, fy_8 * _S801.y + cy_5);
                    break;
                }
                _S677 = _S675 & (!_S679);
                break;
            }
            if(!_S677)
            {
                _S679 = false;
                break;
            }
            float2  _S802 = *mean2d_8 + make_float2 (_S690.w_mean_0[int(0)]) * *_S668 + make_float2 (_S690.w_mean_0[int(1)]) * *_S669 + make_float2 (_S690.w_mean_0[int(2)]) * *_S670 + make_float2 (_S690.w_mean_0[int(3)]) * *_S672 + make_float2 (_S690.w_mean_0[int(4)]) * *_S673 + make_float2 (_S690.w_mean_0[int(5)]) * *_S674 + make_float2 (_S690.w_mean_0[int(6)]) * *_S676;
            *mean2d_8 = _S802;
            float _S803 = - _S695;
            float _S804 = - _S696;
            float2  _S805 = make_float2 (clamp_0(_S802.x, _S803, lim_x_pos_2), clamp_0(_S802.y, _S804, lim_y_pos_0));
            float2  d_0 = make_float2 (clamp_0((*_S668).x, _S803, lim_x_pos_2), clamp_0((*_S668).y, _S804, lim_y_pos_0)) - _S805;
            float _S806 = d_0.x;
            float _S807 = d_0.y;
            float _S808 = _S806 * _S807;
            float2  d_1 = make_float2 (clamp_0((*_S669).x, _S803, lim_x_pos_2), clamp_0((*_S669).y, _S804, lim_y_pos_0)) - _S805;
            float _S809 = d_1.x;
            float _S810 = d_1.y;
            float _S811 = _S809 * _S810;
            float2  d_2 = make_float2 (clamp_0((*_S670).x, _S803, lim_x_pos_2), clamp_0((*_S670).y, _S804, lim_y_pos_0)) - _S805;
            float _S812 = d_2.x;
            float _S813 = d_2.y;
            float _S814 = _S812 * _S813;
            float2  d_3 = make_float2 (clamp_0((*_S672).x, _S803, lim_x_pos_2), clamp_0((*_S672).y, _S804, lim_y_pos_0)) - _S805;
            float _S815 = d_3.x;
            float _S816 = d_3.y;
            float _S817 = _S815 * _S816;
            float2  d_4 = make_float2 (clamp_0((*_S673).x, _S803, lim_x_pos_2), clamp_0((*_S673).y, _S804, lim_y_pos_0)) - _S805;
            float _S818 = d_4.x;
            float _S819 = d_4.y;
            float _S820 = _S818 * _S819;
            float2  d_5 = make_float2 (clamp_0((*_S674).x, _S803, lim_x_pos_2), clamp_0((*_S674).y, _S804, lim_y_pos_0)) - _S805;
            float _S821 = d_5.x;
            float _S822 = d_5.y;
            float _S823 = _S821 * _S822;
            float2  d_6 = make_float2 (clamp_0((*_S676).x, _S803, lim_x_pos_2), clamp_0((*_S676).y, _S804, lim_y_pos_0)) - _S805;
            float _S824 = d_6.x;
            float _S825 = d_6.y;
            float _S826 = _S824 * _S825;
            covar2d_4 = covar2d_4 + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S806 * _S806, _S808, _S808, _S807 * _S807) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S809 * _S809, _S811, _S811, _S810 * _S810) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S812 * _S812, _S814, _S814, _S813 * _S813) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S815 * _S815, _S817, _S817, _S816 * _S816) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S818 * _S818, _S820, _S820, _S819 * _S819) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S821 * _S821, _S823, _S823, _S822 * _S822) + makeMatrix<float, 2, 2> (_S690.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S824 * _S824, _S826, _S826, _S825 * _S825);
            _S679 = true;
            break;
        }
        if(!(true & _S679))
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        float eps2d_4;
        if(antialiased_4)
        {
            eps2d_4 = 0.10000000149011612f;
        }
        else
        {
            eps2d_4 = 0.30000001192092896f;
        }
        float det_orig_4 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S827 = *&(((&covar2d_4)->rows + (int(0)))->x) + eps2d_4;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S827;
        float _S828 = *&(((&covar2d_4)->rows + (int(1)))->y) + eps2d_4;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S828;
        float det_blur_4 = _S827 * _S828 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / det_blur_4))))));
        if(det_blur_4 <= 0.0f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        float invdet_6 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
        Matrix<float, 2, 2>  _S829 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_6, - covar2d_4.rows[int(0)].y * invdet_6, - covar2d_4.rows[int(1)].x * invdet_6, covar2d_4.rows[int(0)].x * invdet_6);
        if(antialiased_4)
        {
            *opacity_4 = *opacity_4 * compensation_4;
        }
        if((*opacity_4) < 0.00392156885936856f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        float _S830 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
        float radius_x_4 = _S830 * (F32_sqrt((covar2d_4[int(0)].x)));
        float radius_y_4 = _S830 * (F32_sqrt((covar2d_4[int(1)].y)));
        float _S831 = (*mean2d_8).x - radius_x_4;
        float _S832 = (*mean2d_8).x + radius_x_4;
        float _S833 = (*mean2d_8).y - radius_y_4;
        float _S834 = (*mean2d_8).y + radius_y_4;
        if(_S832 <= 0.0f)
        {
            _S679 = true;
        }
        else
        {
            _S679 = _S831 >= _S666;
        }
        if(_S679)
        {
            _S679 = true;
        }
        else
        {
            _S679 = _S834 <= 0.0f;
        }
        if(_S679)
        {
            _S679 = true;
        }
        else
        {
            _S679 = _S833 >= _S667;
        }
        if(_S679)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_4 = make_float4 (_S831, _S833, _S832, _S834);
        *sorting_depth_4 = _S678;
        *conic_4 = make_float3 (_S829.rows[int(0)].x, _S829.rows[int(0)].y, _S829.rows[int(1)].y);
        *radius_5 = view_radius_3dgs_0(mean_5, scale_4, in_opacity_4, - mul_6(transpose_3(R_4), t_4));
        break;
    }
    return;
}

inline __device__ void projection_3dgut_fisheye(bool antialiased_5, float3  mean_6, float4  quat_5, float3  scale_5, float in_opacity_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_9, float fy_9, float cx_6, float cy_6, FixedArray<float, 10>  dist_coeffs_9, uint image_width_5, uint image_height_5, float4  * aabb_xyxy_5, float * sorting_depth_5, float * radius_6, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5)
{
    float2  * _S835;
    float _S836;
    float2  _S837;
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
    float2  _S850;
    bool _S851;
    float2  * _S852;
    float _S853;
    bool _S854;
    float2  * _S855;
    float _S856;
    bool _S857;
    bool _S858;
    float2  * _S859;
    float _S860;
    bool _S861;
    float2  * _S862;
    float _S863;
    bool _S864;
    float2  * _S865;
    float _S866;
    bool _S867;
    bool _S868;
    float2  * _S869;
    float _S870;
    bool _S871;
    bool _S872;
    for(;;)
    {
        float3  mean_c_5 = mul_6(R_5, mean_6) + t_5;
        float _S873 = length_1(mean_c_5);
        *depth_5 = _S873;
        if(_S873 <= 0.0f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_5;
        *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
        bool _S874;
        float3  _S875 = exp_0(scale_5);
        float4  _S876 = normalize_0(quat_5);
        float x_22 = _S876.y;
        float x2_5 = x_22 * x_22;
        float y2_5 = _S876.z * _S876.z;
        float z2_5 = _S876.w * _S876.w;
        float xy_5 = _S876.y * _S876.z;
        float xz_5 = _S876.y * _S876.w;
        float yz_5 = _S876.z * _S876.w;
        float wx_5 = _S876.x * _S876.y;
        float wy_5 = _S876.x * _S876.z;
        float wz_5 = _S876.x * _S876.w;
        Matrix<float, 3, 3>  _S877 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_6;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S878 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S878 * _S875.x) * _S877.rows[0U];
        float3  _S879 = mean_6 + delta_3;
        float3  _S880 = mean_6 - delta_3;
        float3  delta_4 = make_float3 (_S878 * _S875.y) * _S877.rows[1U];
        float3  _S881 = mean_6 + delta_4;
        float3  _S882 = mean_6 - delta_4;
        float3  delta_5 = make_float3 (_S878 * _S875.z) * _S877.rows[2U];
        float3  _S883 = mean_6 + delta_5;
        float3  _S884 = mean_6 - delta_5;
        (&ret_1)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_1)->p_0[0U] = mul_6(R_5, (&ret_1)->p_0[0U]) + t_5;
        (&ret_1)->p_0[1U] = mul_6(R_5, _S879) + t_5;
        (&ret_1)->p_0[2U] = mul_6(R_5, _S881) + t_5;
        (&ret_1)->p_0[3U] = mul_6(R_5, _S883) + t_5;
        (&ret_1)->p_0[4U] = mul_6(R_5, _S880) + t_5;
        (&ret_1)->p_0[5U] = mul_6(R_5, _S882) + t_5;
        (&ret_1)->p_0[6U] = mul_6(R_5, _S884) + t_5;
        SigmaPoints_0 _S885 = ret_1;
        for(;;)
        {
            int2  _S886 = make_int2 (int(0));
            float2  _S887 = make_float2 ((float)_S886.x, (float)_S886.y);
            *mean2d_9 = _S887;
            covar2d_5 = makeMatrix<float, 2, 2> (0.0f);
            FixedArray<float2 , 7>  proj_points_1;
            for(;;)
            {
                float k_4;
                _S835 = &proj_points_1[int(0)];
                for(;;)
                {
                    float2  _S888 = float2 {_S885.p_0[int(0)].x, _S885.p_0[int(0)].y};
                    float r_11 = length_0(_S888);
                    float _S889 = _S885.p_0[int(0)].z;
                    _S836 = _S889;
                    float theta_4 = (F32_atan2((r_11), (_S889)));
                    if(theta_4 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S889;
                    }
                    else
                    {
                        k_4 = theta_4 / r_11;
                    }
                    float2  _S890 = _S888 * make_float2 (k_4);
                    proj_points_1[int(0)] = _S890;
                    float2  _S891 = make_float2 (1.0f, 0.0f);
                    _S837 = _S891;
                    _S838 = dist_coeffs_9[int(0)];
                    _S839 = dist_coeffs_9[int(1)];
                    _S840 = dist_coeffs_9[int(2)];
                    _S841 = dist_coeffs_9[int(3)];
                    _S842 = dist_coeffs_9[int(4)];
                    _S843 = dist_coeffs_9[int(5)];
                    _S844 = dist_coeffs_9[int(6)];
                    _S845 = dist_coeffs_9[int(7)];
                    _S846 = dist_coeffs_9[int(8)];
                    _S847 = dist_coeffs_9[int(9)];
                    float u_44 = _S890.x;
                    float v_44 = _S890.y;
                    float _S892 = u_44 + u_44;
                    float r2_44 = u_44 * u_44 + v_44 * v_44;
                    float _S893 = dist_coeffs_9[int(2)] + r2_44 * dist_coeffs_9[int(3)];
                    float _S894 = dist_coeffs_9[int(1)] + r2_44 * _S893;
                    float _S895 = dist_coeffs_9[int(0)] + r2_44 * _S894;
                    float _S896 = _S892 * _S895 + (_S892 * _S894 + (_S892 * _S893 + _S892 * dist_coeffs_9[int(3)] * r2_44) * r2_44) * r2_44;
                    float radial_13 = 1.0f + r2_44 * _S895;
                    float _S897 = 2.0f * dist_coeffs_9[int(4)];
                    _S848 = _S897;
                    float _S898 = _S897 * u_44;
                    float _S899 = 2.0f * u_44;
                    float s_diff_du_4 = _S897 * v_44 + (_S892 + (_S899 + _S899)) * dist_coeffs_9[int(5)] + _S892 * dist_coeffs_9[int(6)];
                    float _S900 = 2.0f * dist_coeffs_9[int(5)];
                    _S849 = _S900;
                    float _S901 = _S900 * u_44;
                    float _S902 = 2.0f * v_44;
                    float2  _S903 = _S891 * make_float2 (radial_13) + make_float2 (_S896) * _S890 + make_float2 (s_diff_du_4, _S900 * v_44 + _S892 * dist_coeffs_9[int(4)] + _S892 * dist_coeffs_9[int(7)]);
                    float2  _S904 = _S903 + make_float2 (_S903.x * dist_coeffs_9[int(8)] + _S903.y * dist_coeffs_9[int(9)], 0.0f);
                    float2  _S905 = make_float2 (0.0f, 1.0f);
                    _S850 = _S905;
                    float _S906 = v_44 + v_44;
                    float2  _S907 = _S905 * make_float2 (radial_13) + make_float2 (_S906 * _S895 + (_S906 * _S894 + (_S906 * _S893 + _S906 * dist_coeffs_9[int(3)] * r2_44) * r2_44) * r2_44) * _S890 + make_float2 (_S898 + _S906 * dist_coeffs_9[int(5)] + _S906 * dist_coeffs_9[int(6)], _S901 + (_S906 + (_S902 + _S902)) * dist_coeffs_9[int(4)] + _S906 * dist_coeffs_9[int(7)]);
                    Matrix<float, 2, 2>  _S908 = transpose_0(makeMatrix<float, 2, 2> (_S904, _S907 + make_float2 (_S907.x * dist_coeffs_9[int(8)] + _S907.y * dist_coeffs_9[int(9)], 0.0f)));
                    bool _S909 = !((F32_min((determinant_0(_S908)), ((F32_min((_S908.rows[int(0)].x), (_S908.rows[int(1)].y)))))) > 0.0f);
                    _S851 = _S909;
                    if(_S909)
                    {
                        break;
                    }
                    float u_45 = proj_points_1[int(0)].x;
                    float v_45 = proj_points_1[int(0)].y;
                    float r2_45 = u_45 * u_45 + v_45 * v_45;
                    float2  _S910 = proj_points_1[int(0)] * make_float2 (1.0f + r2_45 * (dist_coeffs_9[int(0)] + r2_45 * (dist_coeffs_9[int(1)] + r2_45 * (dist_coeffs_9[int(2)] + r2_45 * dist_coeffs_9[int(3)])))) + make_float2 (_S897 * u_45 * v_45 + dist_coeffs_9[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + dist_coeffs_9[int(6)] * r2_45, _S900 * u_45 * v_45 + dist_coeffs_9[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + dist_coeffs_9[int(7)] * r2_45);
                    float2  _S911 = _S910 + make_float2 (dist_coeffs_9[int(8)] * _S910.x + dist_coeffs_9[int(9)] * _S910.y, 0.0f);
                    proj_points_1[int(0)] = make_float2 (fx_9 * _S911.x + cx_6, fy_9 * _S911.y + cy_6);
                    break;
                }
                bool all_valid_4 = true & (!_S851);
                _S852 = &proj_points_1[int(1)];
                for(;;)
                {
                    float2  _S912 = float2 {_S885.p_0[int(1)].x, _S885.p_0[int(1)].y};
                    float r_12 = length_0(_S912);
                    float _S913 = _S885.p_0[int(1)].z;
                    _S853 = _S913;
                    float theta_5 = (F32_atan2((r_12), (_S913)));
                    if(theta_5 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_5 * theta_5 / 3.0f) / _S913;
                    }
                    else
                    {
                        k_4 = theta_5 / r_12;
                    }
                    float2  _S914 = _S912 * make_float2 (k_4);
                    proj_points_1[int(1)] = _S914;
                    float u_46 = _S914.x;
                    float v_46 = _S914.y;
                    float _S915 = u_46 + u_46;
                    float r2_46 = u_46 * u_46 + v_46 * v_46;
                    float _S916 = _S840 + r2_46 * _S841;
                    float _S917 = _S839 + r2_46 * _S916;
                    float _S918 = _S838 + r2_46 * _S917;
                    float radial_14 = 1.0f + r2_46 * _S918;
                    float _S919 = 2.0f * u_46;
                    float _S920 = 2.0f * v_46;
                    float2  _S921 = _S837 * make_float2 (radial_14) + make_float2 (_S915 * _S918 + (_S915 * _S917 + (_S915 * _S916 + _S915 * _S841 * r2_46) * r2_46) * r2_46) * _S914 + make_float2 (_S848 * v_46 + (_S915 + (_S919 + _S919)) * _S843 + _S915 * _S844, _S849 * v_46 + _S915 * _S842 + _S915 * _S845);
                    float _S922 = v_46 + v_46;
                    float2  _S923 = _S850 * make_float2 (radial_14) + make_float2 (_S922 * _S918 + (_S922 * _S917 + (_S922 * _S916 + _S922 * _S841 * r2_46) * r2_46) * r2_46) * _S914 + make_float2 (_S848 * u_46 + _S922 * _S843 + _S922 * _S844, _S849 * u_46 + (_S922 + (_S920 + _S920)) * _S842 + _S922 * _S845);
                    Matrix<float, 2, 2>  _S924 = transpose_0(makeMatrix<float, 2, 2> (_S921 + make_float2 (_S921.x * _S846 + _S921.y * _S847, 0.0f), _S923 + make_float2 (_S923.x * _S846 + _S923.y * _S847, 0.0f)));
                    bool _S925 = !((F32_min((determinant_0(_S924)), ((F32_min((_S924.rows[int(0)].x), (_S924.rows[int(1)].y)))))) > 0.0f);
                    _S854 = _S925;
                    if(_S925)
                    {
                        break;
                    }
                    float u_47 = proj_points_1[int(1)].x;
                    float v_47 = proj_points_1[int(1)].y;
                    float r2_47 = u_47 * u_47 + v_47 * v_47;
                    float2  _S926 = proj_points_1[int(1)] * make_float2 (1.0f + r2_47 * (_S838 + r2_47 * (_S839 + r2_47 * (_S840 + r2_47 * _S841)))) + make_float2 (_S848 * u_47 * v_47 + _S843 * (r2_47 + 2.0f * u_47 * u_47) + _S844 * r2_47, _S849 * u_47 * v_47 + _S842 * (r2_47 + 2.0f * v_47 * v_47) + _S845 * r2_47);
                    float2  _S927 = _S926 + make_float2 (_S846 * _S926.x + _S847 * _S926.y, 0.0f);
                    proj_points_1[int(1)] = make_float2 (fx_9 * _S927.x + cx_6, fy_9 * _S927.y + cy_6);
                    break;
                }
                bool all_valid_5 = all_valid_4 & (!_S854);
                for(;;)
                {
                    _S855 = &proj_points_1[int(2)];
                    for(;;)
                    {
                        float2  _S928 = float2 {_S885.p_0[int(2)].x, _S885.p_0[int(2)].y};
                        float r_13 = length_0(_S928);
                        float _S929 = _S885.p_0[int(2)].z;
                        _S856 = _S929;
                        float theta_6 = (F32_atan2((r_13), (_S929)));
                        if(theta_6 < 0.00100000004749745f)
                        {
                            k_4 = (1.0f - theta_6 * theta_6 / 3.0f) / _S929;
                        }
                        else
                        {
                            k_4 = theta_6 / r_13;
                        }
                        float2  _S930 = _S928 * make_float2 (k_4);
                        proj_points_1[int(2)] = _S930;
                        float u_48 = _S930.x;
                        float v_48 = _S930.y;
                        float _S931 = u_48 + u_48;
                        float r2_48 = u_48 * u_48 + v_48 * v_48;
                        float _S932 = _S840 + r2_48 * _S841;
                        float _S933 = _S839 + r2_48 * _S932;
                        float _S934 = _S838 + r2_48 * _S933;
                        float radial_15 = 1.0f + r2_48 * _S934;
                        float _S935 = 2.0f * u_48;
                        float _S936 = 2.0f * v_48;
                        float2  _S937 = _S837 * make_float2 (radial_15) + make_float2 (_S931 * _S934 + (_S931 * _S933 + (_S931 * _S932 + _S931 * _S841 * r2_48) * r2_48) * r2_48) * _S930 + make_float2 (_S848 * v_48 + (_S931 + (_S935 + _S935)) * _S843 + _S931 * _S844, _S849 * v_48 + _S931 * _S842 + _S931 * _S845);
                        float _S938 = v_48 + v_48;
                        float2  _S939 = _S850 * make_float2 (radial_15) + make_float2 (_S938 * _S934 + (_S938 * _S933 + (_S938 * _S932 + _S938 * _S841 * r2_48) * r2_48) * r2_48) * _S930 + make_float2 (_S848 * u_48 + _S938 * _S843 + _S938 * _S844, _S849 * u_48 + (_S938 + (_S936 + _S936)) * _S842 + _S938 * _S845);
                        Matrix<float, 2, 2>  _S940 = transpose_0(makeMatrix<float, 2, 2> (_S937 + make_float2 (_S937.x * _S846 + _S937.y * _S847, 0.0f), _S939 + make_float2 (_S939.x * _S846 + _S939.y * _S847, 0.0f)));
                        bool _S941 = !((F32_min((determinant_0(_S940)), ((F32_min((_S940.rows[int(0)].x), (_S940.rows[int(1)].y)))))) > 0.0f);
                        _S857 = _S941;
                        if(_S941)
                        {
                            break;
                        }
                        float u_49 = proj_points_1[int(2)].x;
                        float v_49 = proj_points_1[int(2)].y;
                        float r2_49 = u_49 * u_49 + v_49 * v_49;
                        float2  _S942 = proj_points_1[int(2)] * make_float2 (1.0f + r2_49 * (_S838 + r2_49 * (_S839 + r2_49 * (_S840 + r2_49 * _S841)))) + make_float2 (_S848 * u_49 * v_49 + _S843 * (r2_49 + 2.0f * u_49 * u_49) + _S844 * r2_49, _S849 * u_49 * v_49 + _S842 * (r2_49 + 2.0f * v_49 * v_49) + _S845 * r2_49);
                        float2  _S943 = _S942 + make_float2 (_S846 * _S942.x + _S847 * _S942.y, 0.0f);
                        proj_points_1[int(2)] = make_float2 (fx_9 * _S943.x + cx_6, fy_9 * _S943.y + cy_6);
                        break;
                    }
                    _S858 = all_valid_5 & (!_S857);
                    break;
                }
                _S859 = &proj_points_1[int(3)];
                for(;;)
                {
                    float2  _S944 = float2 {_S885.p_0[int(3)].x, _S885.p_0[int(3)].y};
                    float r_14 = length_0(_S944);
                    float _S945 = _S885.p_0[int(3)].z;
                    _S860 = _S945;
                    float theta_7 = (F32_atan2((r_14), (_S945)));
                    if(theta_7 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_7 * theta_7 / 3.0f) / _S945;
                    }
                    else
                    {
                        k_4 = theta_7 / r_14;
                    }
                    float2  _S946 = _S944 * make_float2 (k_4);
                    proj_points_1[int(3)] = _S946;
                    float u_50 = _S946.x;
                    float v_50 = _S946.y;
                    float _S947 = u_50 + u_50;
                    float r2_50 = u_50 * u_50 + v_50 * v_50;
                    float _S948 = _S840 + r2_50 * _S841;
                    float _S949 = _S839 + r2_50 * _S948;
                    float _S950 = _S838 + r2_50 * _S949;
                    float radial_16 = 1.0f + r2_50 * _S950;
                    float _S951 = 2.0f * u_50;
                    float _S952 = 2.0f * v_50;
                    float2  _S953 = _S837 * make_float2 (radial_16) + make_float2 (_S947 * _S950 + (_S947 * _S949 + (_S947 * _S948 + _S947 * _S841 * r2_50) * r2_50) * r2_50) * _S946 + make_float2 (_S848 * v_50 + (_S947 + (_S951 + _S951)) * _S843 + _S947 * _S844, _S849 * v_50 + _S947 * _S842 + _S947 * _S845);
                    float _S954 = v_50 + v_50;
                    float2  _S955 = _S850 * make_float2 (radial_16) + make_float2 (_S954 * _S950 + (_S954 * _S949 + (_S954 * _S948 + _S954 * _S841 * r2_50) * r2_50) * r2_50) * _S946 + make_float2 (_S848 * u_50 + _S954 * _S843 + _S954 * _S844, _S849 * u_50 + (_S954 + (_S952 + _S952)) * _S842 + _S954 * _S845);
                    Matrix<float, 2, 2>  _S956 = transpose_0(makeMatrix<float, 2, 2> (_S953 + make_float2 (_S953.x * _S846 + _S953.y * _S847, 0.0f), _S955 + make_float2 (_S955.x * _S846 + _S955.y * _S847, 0.0f)));
                    bool _S957 = !((F32_min((determinant_0(_S956)), ((F32_min((_S956.rows[int(0)].x), (_S956.rows[int(1)].y)))))) > 0.0f);
                    _S861 = _S957;
                    if(_S957)
                    {
                        break;
                    }
                    float u_51 = proj_points_1[int(3)].x;
                    float v_51 = proj_points_1[int(3)].y;
                    float r2_51 = u_51 * u_51 + v_51 * v_51;
                    float2  _S958 = proj_points_1[int(3)] * make_float2 (1.0f + r2_51 * (_S838 + r2_51 * (_S839 + r2_51 * (_S840 + r2_51 * _S841)))) + make_float2 (_S848 * u_51 * v_51 + _S843 * (r2_51 + 2.0f * u_51 * u_51) + _S844 * r2_51, _S849 * u_51 * v_51 + _S842 * (r2_51 + 2.0f * v_51 * v_51) + _S845 * r2_51);
                    float2  _S959 = _S958 + make_float2 (_S846 * _S958.x + _S847 * _S958.y, 0.0f);
                    proj_points_1[int(3)] = make_float2 (fx_9 * _S959.x + cx_6, fy_9 * _S959.y + cy_6);
                    break;
                }
                bool all_valid_6 = _S858 & (!_S861);
                _S862 = &proj_points_1[int(4)];
                for(;;)
                {
                    float2  _S960 = float2 {_S885.p_0[int(4)].x, _S885.p_0[int(4)].y};
                    float r_15 = length_0(_S960);
                    float _S961 = _S885.p_0[int(4)].z;
                    _S863 = _S961;
                    float theta_8 = (F32_atan2((r_15), (_S961)));
                    if(theta_8 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_8 * theta_8 / 3.0f) / _S961;
                    }
                    else
                    {
                        k_4 = theta_8 / r_15;
                    }
                    float2  _S962 = _S960 * make_float2 (k_4);
                    proj_points_1[int(4)] = _S962;
                    float u_52 = _S962.x;
                    float v_52 = _S962.y;
                    float _S963 = u_52 + u_52;
                    float r2_52 = u_52 * u_52 + v_52 * v_52;
                    float _S964 = _S840 + r2_52 * _S841;
                    float _S965 = _S839 + r2_52 * _S964;
                    float _S966 = _S838 + r2_52 * _S965;
                    float radial_17 = 1.0f + r2_52 * _S966;
                    float _S967 = 2.0f * u_52;
                    float _S968 = 2.0f * v_52;
                    float2  _S969 = _S837 * make_float2 (radial_17) + make_float2 (_S963 * _S966 + (_S963 * _S965 + (_S963 * _S964 + _S963 * _S841 * r2_52) * r2_52) * r2_52) * _S962 + make_float2 (_S848 * v_52 + (_S963 + (_S967 + _S967)) * _S843 + _S963 * _S844, _S849 * v_52 + _S963 * _S842 + _S963 * _S845);
                    float _S970 = v_52 + v_52;
                    float2  _S971 = _S850 * make_float2 (radial_17) + make_float2 (_S970 * _S966 + (_S970 * _S965 + (_S970 * _S964 + _S970 * _S841 * r2_52) * r2_52) * r2_52) * _S962 + make_float2 (_S848 * u_52 + _S970 * _S843 + _S970 * _S844, _S849 * u_52 + (_S970 + (_S968 + _S968)) * _S842 + _S970 * _S845);
                    Matrix<float, 2, 2>  _S972 = transpose_0(makeMatrix<float, 2, 2> (_S969 + make_float2 (_S969.x * _S846 + _S969.y * _S847, 0.0f), _S971 + make_float2 (_S971.x * _S846 + _S971.y * _S847, 0.0f)));
                    bool _S973 = !((F32_min((determinant_0(_S972)), ((F32_min((_S972.rows[int(0)].x), (_S972.rows[int(1)].y)))))) > 0.0f);
                    _S864 = _S973;
                    if(_S973)
                    {
                        break;
                    }
                    float u_53 = proj_points_1[int(4)].x;
                    float v_53 = proj_points_1[int(4)].y;
                    float r2_53 = u_53 * u_53 + v_53 * v_53;
                    float2  _S974 = proj_points_1[int(4)] * make_float2 (1.0f + r2_53 * (_S838 + r2_53 * (_S839 + r2_53 * (_S840 + r2_53 * _S841)))) + make_float2 (_S848 * u_53 * v_53 + _S843 * (r2_53 + 2.0f * u_53 * u_53) + _S844 * r2_53, _S849 * u_53 * v_53 + _S842 * (r2_53 + 2.0f * v_53 * v_53) + _S845 * r2_53);
                    float2  _S975 = _S974 + make_float2 (_S846 * _S974.x + _S847 * _S974.y, 0.0f);
                    proj_points_1[int(4)] = make_float2 (fx_9 * _S975.x + cx_6, fy_9 * _S975.y + cy_6);
                    break;
                }
                bool all_valid_7 = all_valid_6 & (!_S864);
                for(;;)
                {
                    _S865 = &proj_points_1[int(5)];
                    for(;;)
                    {
                        float2  _S976 = float2 {_S885.p_0[int(5)].x, _S885.p_0[int(5)].y};
                        float r_16 = length_0(_S976);
                        float _S977 = _S885.p_0[int(5)].z;
                        _S866 = _S977;
                        float theta_9 = (F32_atan2((r_16), (_S977)));
                        if(theta_9 < 0.00100000004749745f)
                        {
                            k_4 = (1.0f - theta_9 * theta_9 / 3.0f) / _S977;
                        }
                        else
                        {
                            k_4 = theta_9 / r_16;
                        }
                        float2  _S978 = _S976 * make_float2 (k_4);
                        proj_points_1[int(5)] = _S978;
                        float u_54 = _S978.x;
                        float v_54 = _S978.y;
                        float _S979 = u_54 + u_54;
                        float r2_54 = u_54 * u_54 + v_54 * v_54;
                        float _S980 = _S840 + r2_54 * _S841;
                        float _S981 = _S839 + r2_54 * _S980;
                        float _S982 = _S838 + r2_54 * _S981;
                        float radial_18 = 1.0f + r2_54 * _S982;
                        float _S983 = 2.0f * u_54;
                        float _S984 = 2.0f * v_54;
                        float2  _S985 = _S837 * make_float2 (radial_18) + make_float2 (_S979 * _S982 + (_S979 * _S981 + (_S979 * _S980 + _S979 * _S841 * r2_54) * r2_54) * r2_54) * _S978 + make_float2 (_S848 * v_54 + (_S979 + (_S983 + _S983)) * _S843 + _S979 * _S844, _S849 * v_54 + _S979 * _S842 + _S979 * _S845);
                        float _S986 = v_54 + v_54;
                        float2  _S987 = _S850 * make_float2 (radial_18) + make_float2 (_S986 * _S982 + (_S986 * _S981 + (_S986 * _S980 + _S986 * _S841 * r2_54) * r2_54) * r2_54) * _S978 + make_float2 (_S848 * u_54 + _S986 * _S843 + _S986 * _S844, _S849 * u_54 + (_S986 + (_S984 + _S984)) * _S842 + _S986 * _S845);
                        Matrix<float, 2, 2>  _S988 = transpose_0(makeMatrix<float, 2, 2> (_S985 + make_float2 (_S985.x * _S846 + _S985.y * _S847, 0.0f), _S987 + make_float2 (_S987.x * _S846 + _S987.y * _S847, 0.0f)));
                        bool _S989 = !((F32_min((determinant_0(_S988)), ((F32_min((_S988.rows[int(0)].x), (_S988.rows[int(1)].y)))))) > 0.0f);
                        _S867 = _S989;
                        if(_S989)
                        {
                            break;
                        }
                        float u_55 = proj_points_1[int(5)].x;
                        float v_55 = proj_points_1[int(5)].y;
                        float r2_55 = u_55 * u_55 + v_55 * v_55;
                        float2  _S990 = proj_points_1[int(5)] * make_float2 (1.0f + r2_55 * (_S838 + r2_55 * (_S839 + r2_55 * (_S840 + r2_55 * _S841)))) + make_float2 (_S848 * u_55 * v_55 + _S843 * (r2_55 + 2.0f * u_55 * u_55) + _S844 * r2_55, _S849 * u_55 * v_55 + _S842 * (r2_55 + 2.0f * v_55 * v_55) + _S845 * r2_55);
                        float2  _S991 = _S990 + make_float2 (_S846 * _S990.x + _S847 * _S990.y, 0.0f);
                        proj_points_1[int(5)] = make_float2 (fx_9 * _S991.x + cx_6, fy_9 * _S991.y + cy_6);
                        break;
                    }
                    _S868 = all_valid_7 & (!_S867);
                    break;
                }
                _S869 = &proj_points_1[int(6)];
                for(;;)
                {
                    float2  _S992 = float2 {_S885.p_0[int(6)].x, _S885.p_0[int(6)].y};
                    float r_17 = length_0(_S992);
                    float _S993 = _S885.p_0[int(6)].z;
                    _S870 = _S993;
                    float theta_10 = (F32_atan2((r_17), (_S993)));
                    if(theta_10 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_10 * theta_10 / 3.0f) / _S993;
                    }
                    else
                    {
                        k_4 = theta_10 / r_17;
                    }
                    float2  _S994 = _S992 * make_float2 (k_4);
                    proj_points_1[int(6)] = _S994;
                    float u_56 = _S994.x;
                    float v_56 = _S994.y;
                    float _S995 = u_56 + u_56;
                    float r2_56 = u_56 * u_56 + v_56 * v_56;
                    float _S996 = _S840 + r2_56 * _S841;
                    float _S997 = _S839 + r2_56 * _S996;
                    float _S998 = _S838 + r2_56 * _S997;
                    float radial_19 = 1.0f + r2_56 * _S998;
                    float _S999 = 2.0f * u_56;
                    float _S1000 = 2.0f * v_56;
                    float2  _S1001 = _S837 * make_float2 (radial_19) + make_float2 (_S995 * _S998 + (_S995 * _S997 + (_S995 * _S996 + _S995 * _S841 * r2_56) * r2_56) * r2_56) * _S994 + make_float2 (_S848 * v_56 + (_S995 + (_S999 + _S999)) * _S843 + _S995 * _S844, _S849 * v_56 + _S995 * _S842 + _S995 * _S845);
                    float _S1002 = v_56 + v_56;
                    float2  _S1003 = _S850 * make_float2 (radial_19) + make_float2 (_S1002 * _S998 + (_S1002 * _S997 + (_S1002 * _S996 + _S1002 * _S841 * r2_56) * r2_56) * r2_56) * _S994 + make_float2 (_S848 * u_56 + _S1002 * _S843 + _S1002 * _S844, _S849 * u_56 + (_S1002 + (_S1000 + _S1000)) * _S842 + _S1002 * _S845);
                    Matrix<float, 2, 2>  _S1004 = transpose_0(makeMatrix<float, 2, 2> (_S1001 + make_float2 (_S1001.x * _S846 + _S1001.y * _S847, 0.0f), _S1003 + make_float2 (_S1003.x * _S846 + _S1003.y * _S847, 0.0f)));
                    bool _S1005 = !((F32_min((determinant_0(_S1004)), ((F32_min((_S1004.rows[int(0)].x), (_S1004.rows[int(1)].y)))))) > 0.0f);
                    _S871 = _S1005;
                    if(_S1005)
                    {
                        break;
                    }
                    float u_57 = proj_points_1[int(6)].x;
                    float v_57 = proj_points_1[int(6)].y;
                    float r2_57 = u_57 * u_57 + v_57 * v_57;
                    float2  _S1006 = proj_points_1[int(6)] * make_float2 (1.0f + r2_57 * (_S838 + r2_57 * (_S839 + r2_57 * (_S840 + r2_57 * _S841)))) + make_float2 (_S848 * u_57 * v_57 + _S843 * (r2_57 + 2.0f * u_57 * u_57) + _S844 * r2_57, _S849 * u_57 * v_57 + _S842 * (r2_57 + 2.0f * v_57 * v_57) + _S845 * r2_57);
                    float2  _S1007 = _S1006 + make_float2 (_S846 * _S1006.x + _S847 * _S1006.y, 0.0f);
                    proj_points_1[int(6)] = make_float2 (fx_9 * _S1007.x + cx_6, fy_9 * _S1007.y + cy_6);
                    break;
                }
                _S872 = _S868 & (!_S871);
                break;
            }
            if(!_S872)
            {
                _S874 = false;
                break;
            }
            float2  p_1 = *_S835 + (*_S852 - *_S835) * make_float2 (3.32899999618530273f);
            float2  p_2 = *_S835 + (*_S855 - *_S835) * make_float2 (3.32899999618530273f);
            float2  p_3 = *_S835 + (*_S859 - *_S835) * make_float2 (3.32899999618530273f);
            float2  p_4 = *_S835 + (*_S862 - *_S835) * make_float2 (3.32899999618530273f);
            float2  p_5 = *_S835 + (*_S865 - *_S835) * make_float2 (3.32899999618530273f);
            float2  p_6 = *_S835 + (*_S869 - *_S835) * make_float2 (3.32899999618530273f);
            float2  _S1008 = make_float2 (cx_6, cy_6);
            float2  min_p_0 = min_0(min_0(min_0(min_0(min_0(min_0(*_S835, p_1), p_2), p_3), p_4), p_5), p_6) - _S1008;
            float2  max_p_0 = max_0(max_0(max_0(max_0(max_0(max_0(*_S835, p_1), p_2), p_3), p_4), p_5), p_6) - _S1008;
            if((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S836), (_S853)))), (_S856)))), (_S860)))), (_S863)))), (_S866)))), (_S870))) <= 0.0f)
            {
                _S874 = (min_p_0.x * max_p_0.x) < 0.0f;
            }
            else
            {
                _S874 = false;
            }
            if(_S874)
            {
                _S874 = (min_p_0.y * max_p_0.y) < 0.0f;
            }
            else
            {
                _S874 = false;
            }
            if(_S874)
            {
                _S874 = false;
                break;
            }
            float2  _S1009 = *mean2d_9 + make_float2 (_S885.w_mean_0[int(0)]) * *_S835 + make_float2 (_S885.w_mean_0[int(1)]) * *_S852 + make_float2 (_S885.w_mean_0[int(2)]) * *_S855 + make_float2 (_S885.w_mean_0[int(3)]) * *_S859 + make_float2 (_S885.w_mean_0[int(4)]) * *_S862 + make_float2 (_S885.w_mean_0[int(5)]) * *_S865 + make_float2 (_S885.w_mean_0[int(6)]) * *_S869;
            *mean2d_9 = _S1009;
            float2  d_7 = *_S835 - _S1009;
            float _S1010 = d_7.x;
            float _S1011 = d_7.y;
            float _S1012 = _S1010 * _S1011;
            float2  d_8 = *_S852 - _S1009;
            float _S1013 = d_8.x;
            float _S1014 = d_8.y;
            float _S1015 = _S1013 * _S1014;
            float2  d_9 = *_S855 - _S1009;
            float _S1016 = d_9.x;
            float _S1017 = d_9.y;
            float _S1018 = _S1016 * _S1017;
            float2  d_10 = *_S859 - _S1009;
            float _S1019 = d_10.x;
            float _S1020 = d_10.y;
            float _S1021 = _S1019 * _S1020;
            float2  d_11 = *_S862 - _S1009;
            float _S1022 = d_11.x;
            float _S1023 = d_11.y;
            float _S1024 = _S1022 * _S1023;
            float2  d_12 = *_S865 - _S1009;
            float _S1025 = d_12.x;
            float _S1026 = d_12.y;
            float _S1027 = _S1025 * _S1026;
            float2  d_13 = *_S869 - _S1009;
            float _S1028 = d_13.x;
            float _S1029 = d_13.y;
            float _S1030 = _S1028 * _S1029;
            covar2d_5 = covar2d_5 + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S1010 * _S1010, _S1012, _S1012, _S1011 * _S1011) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S1013 * _S1013, _S1015, _S1015, _S1014 * _S1014) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S1016 * _S1016, _S1018, _S1018, _S1017 * _S1017) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S1019 * _S1019, _S1021, _S1021, _S1020 * _S1020) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S1022 * _S1022, _S1024, _S1024, _S1023 * _S1023) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S1025 * _S1025, _S1027, _S1027, _S1026 * _S1026) + makeMatrix<float, 2, 2> (_S885.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S1028 * _S1028, _S1030, _S1030, _S1029 * _S1029);
            _S874 = true;
            break;
        }
        if(!(true & _S874))
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        float eps2d_5;
        if(antialiased_5)
        {
            eps2d_5 = 0.10000000149011612f;
        }
        else
        {
            eps2d_5 = 0.30000001192092896f;
        }
        float det_orig_5 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
        float _S1031 = *&(((&covar2d_5)->rows + (int(0)))->x) + eps2d_5;
        *&(((&covar2d_5)->rows + (int(0)))->x) = _S1031;
        float _S1032 = *&(((&covar2d_5)->rows + (int(1)))->y) + eps2d_5;
        *&(((&covar2d_5)->rows + (int(1)))->y) = _S1032;
        float det_blur_5 = _S1031 * _S1032 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
        float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / det_blur_5))))));
        if(det_blur_5 <= 0.0f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        float invdet_7 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
        Matrix<float, 2, 2>  _S1033 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_7, - covar2d_5.rows[int(0)].y * invdet_7, - covar2d_5.rows[int(1)].x * invdet_7, covar2d_5.rows[int(0)].x * invdet_7);
        if(antialiased_5)
        {
            *opacity_5 = *opacity_5 * compensation_5;
        }
        if((*opacity_5) < 0.00392156885936856f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        float _S1034 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
        float radius_x_5 = _S1034 * (F32_sqrt((covar2d_5[int(0)].x)));
        float radius_y_5 = _S1034 * (F32_sqrt((covar2d_5[int(1)].y)));
        float _S1035 = (*mean2d_9).x - radius_x_5;
        float _S1036 = (*mean2d_9).x + radius_x_5;
        float _S1037 = (*mean2d_9).y - radius_y_5;
        float _S1038 = (*mean2d_9).y + radius_y_5;
        if(_S1036 <= 0.0f)
        {
            _S874 = true;
        }
        else
        {
            _S874 = _S1035 >= float(image_width_5);
        }
        if(_S874)
        {
            _S874 = true;
        }
        else
        {
            _S874 = _S1038 <= 0.0f;
        }
        if(_S874)
        {
            _S874 = true;
        }
        else
        {
            _S874 = _S1037 >= float(image_height_5);
        }
        if(_S874)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_5 = make_float4 (_S1035, _S1037, _S1036, _S1038);
        float x_23 = mean_c_5.x;
        float y_7 = mean_c_5.y;
        float z_0 = mean_c_5.z;
        float _S1039 = x_23 * x_23 + y_7 * y_7;
        *sorting_depth_5 = z_0 * z_0 * z_0 * z_0 + 0.001953125f * _S1039 * _S1039;
        *conic_5 = make_float3 (_S1033.rows[int(0)].x, _S1033.rows[int(0)].y, _S1033.rows[int(1)].y);
        *radius_6 = view_radius_3dgs_0(mean_6, scale_5, in_opacity_5, - mul_6(transpose_3(R_5), t_5));
        break;
    }
    return;
}

inline __device__ void projection_3dgut_equisolid(bool antialiased_6, float3  mean_7, float4  quat_6, float3  scale_6, float in_opacity_6, Matrix<float, 3, 3>  R_6, float3  t_6, float fx_10, float fy_10, float cx_7, float cy_7, FixedArray<float, 10>  dist_coeffs_10, uint image_width_6, uint image_height_6, float4  * aabb_xyxy_6, float * sorting_depth_6, float * radius_7, float2  * mean2d_10, float * depth_6, float3  * conic_6, float * opacity_6)
{
    float2  * _S1040;
    float _S1041;
    float2  _S1042;
    float _S1043;
    float _S1044;
    float _S1045;
    float _S1046;
    float _S1047;
    float _S1048;
    float _S1049;
    float _S1050;
    float _S1051;
    float _S1052;
    float _S1053;
    float _S1054;
    float2  _S1055;
    bool _S1056;
    float2  * _S1057;
    float _S1058;
    bool _S1059;
    float2  * _S1060;
    float _S1061;
    bool _S1062;
    bool _S1063;
    float2  * _S1064;
    float _S1065;
    bool _S1066;
    float2  * _S1067;
    float _S1068;
    bool _S1069;
    float2  * _S1070;
    float _S1071;
    bool _S1072;
    bool _S1073;
    float2  * _S1074;
    float _S1075;
    bool _S1076;
    bool _S1077;
    for(;;)
    {
        float3  mean_c_6 = mul_6(R_6, mean_7) + t_6;
        float _S1078 = length_1(mean_c_6);
        *depth_6 = _S1078;
        if(_S1078 <= 0.0f)
        {
            *aabb_xyxy_6 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_6;
        *opacity_6 = 1.0f / (1.0f + (F32_exp((- in_opacity_6))));
        bool _S1079;
        float3  _S1080 = exp_0(scale_6);
        float4  _S1081 = normalize_0(quat_6);
        float x_24 = _S1081.y;
        float x2_6 = x_24 * x_24;
        float y2_6 = _S1081.z * _S1081.z;
        float z2_6 = _S1081.w * _S1081.w;
        float xy_6 = _S1081.y * _S1081.z;
        float xz_6 = _S1081.y * _S1081.w;
        float yz_6 = _S1081.z * _S1081.w;
        float wx_6 = _S1081.x * _S1081.y;
        float wy_6 = _S1081.x * _S1081.z;
        float wz_6 = _S1081.x * _S1081.w;
        Matrix<float, 3, 3>  _S1082 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))));
        SigmaPoints_0 ret_2;
        (&ret_2)->p_0[int(0)] = mean_7;
        (&ret_2)->w_mean_0[int(0)] = 0.0f;
        (&ret_2)->w_cov_0[int(0)] = 2.0f;
        float _S1083 = (F32_sqrt((3.0f)));
        float3  delta_6 = make_float3 (_S1083 * _S1080.x) * _S1082.rows[0U];
        float3  _S1084 = mean_7 + delta_6;
        float3  _S1085 = mean_7 - delta_6;
        float3  delta_7 = make_float3 (_S1083 * _S1080.y) * _S1082.rows[1U];
        float3  _S1086 = mean_7 + delta_7;
        float3  _S1087 = mean_7 - delta_7;
        float3  delta_8 = make_float3 (_S1083 * _S1080.z) * _S1082.rows[2U];
        float3  _S1088 = mean_7 + delta_8;
        float3  _S1089 = mean_7 - delta_8;
        (&ret_2)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_2)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_2)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_2)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_2)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_2)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_2)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_2)->p_0[0U] = mul_6(R_6, (&ret_2)->p_0[0U]) + t_6;
        (&ret_2)->p_0[1U] = mul_6(R_6, _S1084) + t_6;
        (&ret_2)->p_0[2U] = mul_6(R_6, _S1086) + t_6;
        (&ret_2)->p_0[3U] = mul_6(R_6, _S1088) + t_6;
        (&ret_2)->p_0[4U] = mul_6(R_6, _S1085) + t_6;
        (&ret_2)->p_0[5U] = mul_6(R_6, _S1087) + t_6;
        (&ret_2)->p_0[6U] = mul_6(R_6, _S1089) + t_6;
        SigmaPoints_0 _S1090 = ret_2;
        for(;;)
        {
            int2  _S1091 = make_int2 (int(0));
            float2  _S1092 = make_float2 ((float)_S1091.x, (float)_S1091.y);
            *mean2d_10 = _S1092;
            covar2d_6 = makeMatrix<float, 2, 2> (0.0f);
            FixedArray<float2 , 7>  proj_points_2;
            for(;;)
            {
                float k_5;
                _S1040 = &proj_points_2[int(0)];
                for(;;)
                {
                    float2  _S1093 = float2 {_S1090.p_0[int(0)].x, _S1090.p_0[int(0)].y};
                    float r_18 = length_0(_S1093);
                    float _S1094 = _S1090.p_0[int(0)].z;
                    _S1041 = _S1094;
                    float theta_11 = (F32_atan2((r_18), (_S1094)));
                    if(r_18 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_11 * theta_11 / 24.0f) / _S1094;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_11))) / r_18;
                    }
                    float2  _S1095 = _S1093 * make_float2 (k_5);
                    proj_points_2[int(0)] = _S1095;
                    float2  _S1096 = make_float2 (1.0f, 0.0f);
                    _S1042 = _S1096;
                    _S1043 = dist_coeffs_10[int(0)];
                    _S1044 = dist_coeffs_10[int(1)];
                    _S1045 = dist_coeffs_10[int(2)];
                    _S1046 = dist_coeffs_10[int(3)];
                    _S1047 = dist_coeffs_10[int(4)];
                    _S1048 = dist_coeffs_10[int(5)];
                    _S1049 = dist_coeffs_10[int(6)];
                    _S1050 = dist_coeffs_10[int(7)];
                    _S1051 = dist_coeffs_10[int(8)];
                    _S1052 = dist_coeffs_10[int(9)];
                    float u_58 = _S1095.x;
                    float v_58 = _S1095.y;
                    float _S1097 = u_58 + u_58;
                    float r2_58 = u_58 * u_58 + v_58 * v_58;
                    float _S1098 = dist_coeffs_10[int(2)] + r2_58 * dist_coeffs_10[int(3)];
                    float _S1099 = dist_coeffs_10[int(1)] + r2_58 * _S1098;
                    float _S1100 = dist_coeffs_10[int(0)] + r2_58 * _S1099;
                    float _S1101 = _S1097 * _S1100 + (_S1097 * _S1099 + (_S1097 * _S1098 + _S1097 * dist_coeffs_10[int(3)] * r2_58) * r2_58) * r2_58;
                    float radial_20 = 1.0f + r2_58 * _S1100;
                    float _S1102 = 2.0f * dist_coeffs_10[int(4)];
                    _S1053 = _S1102;
                    float _S1103 = _S1102 * u_58;
                    float _S1104 = 2.0f * u_58;
                    float s_diff_du_5 = _S1102 * v_58 + (_S1097 + (_S1104 + _S1104)) * dist_coeffs_10[int(5)] + _S1097 * dist_coeffs_10[int(6)];
                    float _S1105 = 2.0f * dist_coeffs_10[int(5)];
                    _S1054 = _S1105;
                    float _S1106 = _S1105 * u_58;
                    float _S1107 = 2.0f * v_58;
                    float2  _S1108 = _S1096 * make_float2 (radial_20) + make_float2 (_S1101) * _S1095 + make_float2 (s_diff_du_5, _S1105 * v_58 + _S1097 * dist_coeffs_10[int(4)] + _S1097 * dist_coeffs_10[int(7)]);
                    float2  _S1109 = _S1108 + make_float2 (_S1108.x * dist_coeffs_10[int(8)] + _S1108.y * dist_coeffs_10[int(9)], 0.0f);
                    float2  _S1110 = make_float2 (0.0f, 1.0f);
                    _S1055 = _S1110;
                    float _S1111 = v_58 + v_58;
                    float2  _S1112 = _S1110 * make_float2 (radial_20) + make_float2 (_S1111 * _S1100 + (_S1111 * _S1099 + (_S1111 * _S1098 + _S1111 * dist_coeffs_10[int(3)] * r2_58) * r2_58) * r2_58) * _S1095 + make_float2 (_S1103 + _S1111 * dist_coeffs_10[int(5)] + _S1111 * dist_coeffs_10[int(6)], _S1106 + (_S1111 + (_S1107 + _S1107)) * dist_coeffs_10[int(4)] + _S1111 * dist_coeffs_10[int(7)]);
                    Matrix<float, 2, 2>  _S1113 = transpose_0(makeMatrix<float, 2, 2> (_S1109, _S1112 + make_float2 (_S1112.x * dist_coeffs_10[int(8)] + _S1112.y * dist_coeffs_10[int(9)], 0.0f)));
                    bool _S1114 = !((F32_min((determinant_0(_S1113)), ((F32_min((_S1113.rows[int(0)].x), (_S1113.rows[int(1)].y)))))) > 0.0f);
                    _S1056 = _S1114;
                    if(_S1114)
                    {
                        break;
                    }
                    float u_59 = proj_points_2[int(0)].x;
                    float v_59 = proj_points_2[int(0)].y;
                    float r2_59 = u_59 * u_59 + v_59 * v_59;
                    float2  _S1115 = proj_points_2[int(0)] * make_float2 (1.0f + r2_59 * (dist_coeffs_10[int(0)] + r2_59 * (dist_coeffs_10[int(1)] + r2_59 * (dist_coeffs_10[int(2)] + r2_59 * dist_coeffs_10[int(3)])))) + make_float2 (_S1102 * u_59 * v_59 + dist_coeffs_10[int(5)] * (r2_59 + 2.0f * u_59 * u_59) + dist_coeffs_10[int(6)] * r2_59, _S1105 * u_59 * v_59 + dist_coeffs_10[int(4)] * (r2_59 + 2.0f * v_59 * v_59) + dist_coeffs_10[int(7)] * r2_59);
                    float2  _S1116 = _S1115 + make_float2 (dist_coeffs_10[int(8)] * _S1115.x + dist_coeffs_10[int(9)] * _S1115.y, 0.0f);
                    proj_points_2[int(0)] = make_float2 (fx_10 * _S1116.x + cx_7, fy_10 * _S1116.y + cy_7);
                    break;
                }
                bool all_valid_8 = true & (!_S1056);
                _S1057 = &proj_points_2[int(1)];
                for(;;)
                {
                    float2  _S1117 = float2 {_S1090.p_0[int(1)].x, _S1090.p_0[int(1)].y};
                    float r_19 = length_0(_S1117);
                    float _S1118 = _S1090.p_0[int(1)].z;
                    _S1058 = _S1118;
                    float theta_12 = (F32_atan2((r_19), (_S1118)));
                    if(r_19 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_12 * theta_12 / 24.0f) / _S1118;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_12))) / r_19;
                    }
                    float2  _S1119 = _S1117 * make_float2 (k_5);
                    proj_points_2[int(1)] = _S1119;
                    float u_60 = _S1119.x;
                    float v_60 = _S1119.y;
                    float _S1120 = u_60 + u_60;
                    float r2_60 = u_60 * u_60 + v_60 * v_60;
                    float _S1121 = _S1045 + r2_60 * _S1046;
                    float _S1122 = _S1044 + r2_60 * _S1121;
                    float _S1123 = _S1043 + r2_60 * _S1122;
                    float radial_21 = 1.0f + r2_60 * _S1123;
                    float _S1124 = 2.0f * u_60;
                    float _S1125 = 2.0f * v_60;
                    float2  _S1126 = _S1042 * make_float2 (radial_21) + make_float2 (_S1120 * _S1123 + (_S1120 * _S1122 + (_S1120 * _S1121 + _S1120 * _S1046 * r2_60) * r2_60) * r2_60) * _S1119 + make_float2 (_S1053 * v_60 + (_S1120 + (_S1124 + _S1124)) * _S1048 + _S1120 * _S1049, _S1054 * v_60 + _S1120 * _S1047 + _S1120 * _S1050);
                    float _S1127 = v_60 + v_60;
                    float2  _S1128 = _S1055 * make_float2 (radial_21) + make_float2 (_S1127 * _S1123 + (_S1127 * _S1122 + (_S1127 * _S1121 + _S1127 * _S1046 * r2_60) * r2_60) * r2_60) * _S1119 + make_float2 (_S1053 * u_60 + _S1127 * _S1048 + _S1127 * _S1049, _S1054 * u_60 + (_S1127 + (_S1125 + _S1125)) * _S1047 + _S1127 * _S1050);
                    Matrix<float, 2, 2>  _S1129 = transpose_0(makeMatrix<float, 2, 2> (_S1126 + make_float2 (_S1126.x * _S1051 + _S1126.y * _S1052, 0.0f), _S1128 + make_float2 (_S1128.x * _S1051 + _S1128.y * _S1052, 0.0f)));
                    bool _S1130 = !((F32_min((determinant_0(_S1129)), ((F32_min((_S1129.rows[int(0)].x), (_S1129.rows[int(1)].y)))))) > 0.0f);
                    _S1059 = _S1130;
                    if(_S1130)
                    {
                        break;
                    }
                    float u_61 = proj_points_2[int(1)].x;
                    float v_61 = proj_points_2[int(1)].y;
                    float r2_61 = u_61 * u_61 + v_61 * v_61;
                    float2  _S1131 = proj_points_2[int(1)] * make_float2 (1.0f + r2_61 * (_S1043 + r2_61 * (_S1044 + r2_61 * (_S1045 + r2_61 * _S1046)))) + make_float2 (_S1053 * u_61 * v_61 + _S1048 * (r2_61 + 2.0f * u_61 * u_61) + _S1049 * r2_61, _S1054 * u_61 * v_61 + _S1047 * (r2_61 + 2.0f * v_61 * v_61) + _S1050 * r2_61);
                    float2  _S1132 = _S1131 + make_float2 (_S1051 * _S1131.x + _S1052 * _S1131.y, 0.0f);
                    proj_points_2[int(1)] = make_float2 (fx_10 * _S1132.x + cx_7, fy_10 * _S1132.y + cy_7);
                    break;
                }
                bool all_valid_9 = all_valid_8 & (!_S1059);
                for(;;)
                {
                    _S1060 = &proj_points_2[int(2)];
                    for(;;)
                    {
                        float2  _S1133 = float2 {_S1090.p_0[int(2)].x, _S1090.p_0[int(2)].y};
                        float r_20 = length_0(_S1133);
                        float _S1134 = _S1090.p_0[int(2)].z;
                        _S1061 = _S1134;
                        float theta_13 = (F32_atan2((r_20), (_S1134)));
                        if(r_20 < 9.99999997475242708e-07f)
                        {
                            k_5 = (1.0f - theta_13 * theta_13 / 24.0f) / _S1134;
                        }
                        else
                        {
                            k_5 = 2.0f * (F32_sin((0.5f * theta_13))) / r_20;
                        }
                        float2  _S1135 = _S1133 * make_float2 (k_5);
                        proj_points_2[int(2)] = _S1135;
                        float u_62 = _S1135.x;
                        float v_62 = _S1135.y;
                        float _S1136 = u_62 + u_62;
                        float r2_62 = u_62 * u_62 + v_62 * v_62;
                        float _S1137 = _S1045 + r2_62 * _S1046;
                        float _S1138 = _S1044 + r2_62 * _S1137;
                        float _S1139 = _S1043 + r2_62 * _S1138;
                        float radial_22 = 1.0f + r2_62 * _S1139;
                        float _S1140 = 2.0f * u_62;
                        float _S1141 = 2.0f * v_62;
                        float2  _S1142 = _S1042 * make_float2 (radial_22) + make_float2 (_S1136 * _S1139 + (_S1136 * _S1138 + (_S1136 * _S1137 + _S1136 * _S1046 * r2_62) * r2_62) * r2_62) * _S1135 + make_float2 (_S1053 * v_62 + (_S1136 + (_S1140 + _S1140)) * _S1048 + _S1136 * _S1049, _S1054 * v_62 + _S1136 * _S1047 + _S1136 * _S1050);
                        float _S1143 = v_62 + v_62;
                        float2  _S1144 = _S1055 * make_float2 (radial_22) + make_float2 (_S1143 * _S1139 + (_S1143 * _S1138 + (_S1143 * _S1137 + _S1143 * _S1046 * r2_62) * r2_62) * r2_62) * _S1135 + make_float2 (_S1053 * u_62 + _S1143 * _S1048 + _S1143 * _S1049, _S1054 * u_62 + (_S1143 + (_S1141 + _S1141)) * _S1047 + _S1143 * _S1050);
                        Matrix<float, 2, 2>  _S1145 = transpose_0(makeMatrix<float, 2, 2> (_S1142 + make_float2 (_S1142.x * _S1051 + _S1142.y * _S1052, 0.0f), _S1144 + make_float2 (_S1144.x * _S1051 + _S1144.y * _S1052, 0.0f)));
                        bool _S1146 = !((F32_min((determinant_0(_S1145)), ((F32_min((_S1145.rows[int(0)].x), (_S1145.rows[int(1)].y)))))) > 0.0f);
                        _S1062 = _S1146;
                        if(_S1146)
                        {
                            break;
                        }
                        float u_63 = proj_points_2[int(2)].x;
                        float v_63 = proj_points_2[int(2)].y;
                        float r2_63 = u_63 * u_63 + v_63 * v_63;
                        float2  _S1147 = proj_points_2[int(2)] * make_float2 (1.0f + r2_63 * (_S1043 + r2_63 * (_S1044 + r2_63 * (_S1045 + r2_63 * _S1046)))) + make_float2 (_S1053 * u_63 * v_63 + _S1048 * (r2_63 + 2.0f * u_63 * u_63) + _S1049 * r2_63, _S1054 * u_63 * v_63 + _S1047 * (r2_63 + 2.0f * v_63 * v_63) + _S1050 * r2_63);
                        float2  _S1148 = _S1147 + make_float2 (_S1051 * _S1147.x + _S1052 * _S1147.y, 0.0f);
                        proj_points_2[int(2)] = make_float2 (fx_10 * _S1148.x + cx_7, fy_10 * _S1148.y + cy_7);
                        break;
                    }
                    _S1063 = all_valid_9 & (!_S1062);
                    break;
                }
                _S1064 = &proj_points_2[int(3)];
                for(;;)
                {
                    float2  _S1149 = float2 {_S1090.p_0[int(3)].x, _S1090.p_0[int(3)].y};
                    float r_21 = length_0(_S1149);
                    float _S1150 = _S1090.p_0[int(3)].z;
                    _S1065 = _S1150;
                    float theta_14 = (F32_atan2((r_21), (_S1150)));
                    if(r_21 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_14 * theta_14 / 24.0f) / _S1150;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_14))) / r_21;
                    }
                    float2  _S1151 = _S1149 * make_float2 (k_5);
                    proj_points_2[int(3)] = _S1151;
                    float u_64 = _S1151.x;
                    float v_64 = _S1151.y;
                    float _S1152 = u_64 + u_64;
                    float r2_64 = u_64 * u_64 + v_64 * v_64;
                    float _S1153 = _S1045 + r2_64 * _S1046;
                    float _S1154 = _S1044 + r2_64 * _S1153;
                    float _S1155 = _S1043 + r2_64 * _S1154;
                    float radial_23 = 1.0f + r2_64 * _S1155;
                    float _S1156 = 2.0f * u_64;
                    float _S1157 = 2.0f * v_64;
                    float2  _S1158 = _S1042 * make_float2 (radial_23) + make_float2 (_S1152 * _S1155 + (_S1152 * _S1154 + (_S1152 * _S1153 + _S1152 * _S1046 * r2_64) * r2_64) * r2_64) * _S1151 + make_float2 (_S1053 * v_64 + (_S1152 + (_S1156 + _S1156)) * _S1048 + _S1152 * _S1049, _S1054 * v_64 + _S1152 * _S1047 + _S1152 * _S1050);
                    float _S1159 = v_64 + v_64;
                    float2  _S1160 = _S1055 * make_float2 (radial_23) + make_float2 (_S1159 * _S1155 + (_S1159 * _S1154 + (_S1159 * _S1153 + _S1159 * _S1046 * r2_64) * r2_64) * r2_64) * _S1151 + make_float2 (_S1053 * u_64 + _S1159 * _S1048 + _S1159 * _S1049, _S1054 * u_64 + (_S1159 + (_S1157 + _S1157)) * _S1047 + _S1159 * _S1050);
                    Matrix<float, 2, 2>  _S1161 = transpose_0(makeMatrix<float, 2, 2> (_S1158 + make_float2 (_S1158.x * _S1051 + _S1158.y * _S1052, 0.0f), _S1160 + make_float2 (_S1160.x * _S1051 + _S1160.y * _S1052, 0.0f)));
                    bool _S1162 = !((F32_min((determinant_0(_S1161)), ((F32_min((_S1161.rows[int(0)].x), (_S1161.rows[int(1)].y)))))) > 0.0f);
                    _S1066 = _S1162;
                    if(_S1162)
                    {
                        break;
                    }
                    float u_65 = proj_points_2[int(3)].x;
                    float v_65 = proj_points_2[int(3)].y;
                    float r2_65 = u_65 * u_65 + v_65 * v_65;
                    float2  _S1163 = proj_points_2[int(3)] * make_float2 (1.0f + r2_65 * (_S1043 + r2_65 * (_S1044 + r2_65 * (_S1045 + r2_65 * _S1046)))) + make_float2 (_S1053 * u_65 * v_65 + _S1048 * (r2_65 + 2.0f * u_65 * u_65) + _S1049 * r2_65, _S1054 * u_65 * v_65 + _S1047 * (r2_65 + 2.0f * v_65 * v_65) + _S1050 * r2_65);
                    float2  _S1164 = _S1163 + make_float2 (_S1051 * _S1163.x + _S1052 * _S1163.y, 0.0f);
                    proj_points_2[int(3)] = make_float2 (fx_10 * _S1164.x + cx_7, fy_10 * _S1164.y + cy_7);
                    break;
                }
                bool all_valid_10 = _S1063 & (!_S1066);
                _S1067 = &proj_points_2[int(4)];
                for(;;)
                {
                    float2  _S1165 = float2 {_S1090.p_0[int(4)].x, _S1090.p_0[int(4)].y};
                    float r_22 = length_0(_S1165);
                    float _S1166 = _S1090.p_0[int(4)].z;
                    _S1068 = _S1166;
                    float theta_15 = (F32_atan2((r_22), (_S1166)));
                    if(r_22 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_15 * theta_15 / 24.0f) / _S1166;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_15))) / r_22;
                    }
                    float2  _S1167 = _S1165 * make_float2 (k_5);
                    proj_points_2[int(4)] = _S1167;
                    float u_66 = _S1167.x;
                    float v_66 = _S1167.y;
                    float _S1168 = u_66 + u_66;
                    float r2_66 = u_66 * u_66 + v_66 * v_66;
                    float _S1169 = _S1045 + r2_66 * _S1046;
                    float _S1170 = _S1044 + r2_66 * _S1169;
                    float _S1171 = _S1043 + r2_66 * _S1170;
                    float radial_24 = 1.0f + r2_66 * _S1171;
                    float _S1172 = 2.0f * u_66;
                    float _S1173 = 2.0f * v_66;
                    float2  _S1174 = _S1042 * make_float2 (radial_24) + make_float2 (_S1168 * _S1171 + (_S1168 * _S1170 + (_S1168 * _S1169 + _S1168 * _S1046 * r2_66) * r2_66) * r2_66) * _S1167 + make_float2 (_S1053 * v_66 + (_S1168 + (_S1172 + _S1172)) * _S1048 + _S1168 * _S1049, _S1054 * v_66 + _S1168 * _S1047 + _S1168 * _S1050);
                    float _S1175 = v_66 + v_66;
                    float2  _S1176 = _S1055 * make_float2 (radial_24) + make_float2 (_S1175 * _S1171 + (_S1175 * _S1170 + (_S1175 * _S1169 + _S1175 * _S1046 * r2_66) * r2_66) * r2_66) * _S1167 + make_float2 (_S1053 * u_66 + _S1175 * _S1048 + _S1175 * _S1049, _S1054 * u_66 + (_S1175 + (_S1173 + _S1173)) * _S1047 + _S1175 * _S1050);
                    Matrix<float, 2, 2>  _S1177 = transpose_0(makeMatrix<float, 2, 2> (_S1174 + make_float2 (_S1174.x * _S1051 + _S1174.y * _S1052, 0.0f), _S1176 + make_float2 (_S1176.x * _S1051 + _S1176.y * _S1052, 0.0f)));
                    bool _S1178 = !((F32_min((determinant_0(_S1177)), ((F32_min((_S1177.rows[int(0)].x), (_S1177.rows[int(1)].y)))))) > 0.0f);
                    _S1069 = _S1178;
                    if(_S1178)
                    {
                        break;
                    }
                    float u_67 = proj_points_2[int(4)].x;
                    float v_67 = proj_points_2[int(4)].y;
                    float r2_67 = u_67 * u_67 + v_67 * v_67;
                    float2  _S1179 = proj_points_2[int(4)] * make_float2 (1.0f + r2_67 * (_S1043 + r2_67 * (_S1044 + r2_67 * (_S1045 + r2_67 * _S1046)))) + make_float2 (_S1053 * u_67 * v_67 + _S1048 * (r2_67 + 2.0f * u_67 * u_67) + _S1049 * r2_67, _S1054 * u_67 * v_67 + _S1047 * (r2_67 + 2.0f * v_67 * v_67) + _S1050 * r2_67);
                    float2  _S1180 = _S1179 + make_float2 (_S1051 * _S1179.x + _S1052 * _S1179.y, 0.0f);
                    proj_points_2[int(4)] = make_float2 (fx_10 * _S1180.x + cx_7, fy_10 * _S1180.y + cy_7);
                    break;
                }
                bool all_valid_11 = all_valid_10 & (!_S1069);
                for(;;)
                {
                    _S1070 = &proj_points_2[int(5)];
                    for(;;)
                    {
                        float2  _S1181 = float2 {_S1090.p_0[int(5)].x, _S1090.p_0[int(5)].y};
                        float r_23 = length_0(_S1181);
                        float _S1182 = _S1090.p_0[int(5)].z;
                        _S1071 = _S1182;
                        float theta_16 = (F32_atan2((r_23), (_S1182)));
                        if(r_23 < 9.99999997475242708e-07f)
                        {
                            k_5 = (1.0f - theta_16 * theta_16 / 24.0f) / _S1182;
                        }
                        else
                        {
                            k_5 = 2.0f * (F32_sin((0.5f * theta_16))) / r_23;
                        }
                        float2  _S1183 = _S1181 * make_float2 (k_5);
                        proj_points_2[int(5)] = _S1183;
                        float u_68 = _S1183.x;
                        float v_68 = _S1183.y;
                        float _S1184 = u_68 + u_68;
                        float r2_68 = u_68 * u_68 + v_68 * v_68;
                        float _S1185 = _S1045 + r2_68 * _S1046;
                        float _S1186 = _S1044 + r2_68 * _S1185;
                        float _S1187 = _S1043 + r2_68 * _S1186;
                        float radial_25 = 1.0f + r2_68 * _S1187;
                        float _S1188 = 2.0f * u_68;
                        float _S1189 = 2.0f * v_68;
                        float2  _S1190 = _S1042 * make_float2 (radial_25) + make_float2 (_S1184 * _S1187 + (_S1184 * _S1186 + (_S1184 * _S1185 + _S1184 * _S1046 * r2_68) * r2_68) * r2_68) * _S1183 + make_float2 (_S1053 * v_68 + (_S1184 + (_S1188 + _S1188)) * _S1048 + _S1184 * _S1049, _S1054 * v_68 + _S1184 * _S1047 + _S1184 * _S1050);
                        float _S1191 = v_68 + v_68;
                        float2  _S1192 = _S1055 * make_float2 (radial_25) + make_float2 (_S1191 * _S1187 + (_S1191 * _S1186 + (_S1191 * _S1185 + _S1191 * _S1046 * r2_68) * r2_68) * r2_68) * _S1183 + make_float2 (_S1053 * u_68 + _S1191 * _S1048 + _S1191 * _S1049, _S1054 * u_68 + (_S1191 + (_S1189 + _S1189)) * _S1047 + _S1191 * _S1050);
                        Matrix<float, 2, 2>  _S1193 = transpose_0(makeMatrix<float, 2, 2> (_S1190 + make_float2 (_S1190.x * _S1051 + _S1190.y * _S1052, 0.0f), _S1192 + make_float2 (_S1192.x * _S1051 + _S1192.y * _S1052, 0.0f)));
                        bool _S1194 = !((F32_min((determinant_0(_S1193)), ((F32_min((_S1193.rows[int(0)].x), (_S1193.rows[int(1)].y)))))) > 0.0f);
                        _S1072 = _S1194;
                        if(_S1194)
                        {
                            break;
                        }
                        float u_69 = proj_points_2[int(5)].x;
                        float v_69 = proj_points_2[int(5)].y;
                        float r2_69 = u_69 * u_69 + v_69 * v_69;
                        float2  _S1195 = proj_points_2[int(5)] * make_float2 (1.0f + r2_69 * (_S1043 + r2_69 * (_S1044 + r2_69 * (_S1045 + r2_69 * _S1046)))) + make_float2 (_S1053 * u_69 * v_69 + _S1048 * (r2_69 + 2.0f * u_69 * u_69) + _S1049 * r2_69, _S1054 * u_69 * v_69 + _S1047 * (r2_69 + 2.0f * v_69 * v_69) + _S1050 * r2_69);
                        float2  _S1196 = _S1195 + make_float2 (_S1051 * _S1195.x + _S1052 * _S1195.y, 0.0f);
                        proj_points_2[int(5)] = make_float2 (fx_10 * _S1196.x + cx_7, fy_10 * _S1196.y + cy_7);
                        break;
                    }
                    _S1073 = all_valid_11 & (!_S1072);
                    break;
                }
                _S1074 = &proj_points_2[int(6)];
                for(;;)
                {
                    float2  _S1197 = float2 {_S1090.p_0[int(6)].x, _S1090.p_0[int(6)].y};
                    float r_24 = length_0(_S1197);
                    float _S1198 = _S1090.p_0[int(6)].z;
                    _S1075 = _S1198;
                    float theta_17 = (F32_atan2((r_24), (_S1198)));
                    if(r_24 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_17 * theta_17 / 24.0f) / _S1198;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_17))) / r_24;
                    }
                    float2  _S1199 = _S1197 * make_float2 (k_5);
                    proj_points_2[int(6)] = _S1199;
                    float u_70 = _S1199.x;
                    float v_70 = _S1199.y;
                    float _S1200 = u_70 + u_70;
                    float r2_70 = u_70 * u_70 + v_70 * v_70;
                    float _S1201 = _S1045 + r2_70 * _S1046;
                    float _S1202 = _S1044 + r2_70 * _S1201;
                    float _S1203 = _S1043 + r2_70 * _S1202;
                    float radial_26 = 1.0f + r2_70 * _S1203;
                    float _S1204 = 2.0f * u_70;
                    float _S1205 = 2.0f * v_70;
                    float2  _S1206 = _S1042 * make_float2 (radial_26) + make_float2 (_S1200 * _S1203 + (_S1200 * _S1202 + (_S1200 * _S1201 + _S1200 * _S1046 * r2_70) * r2_70) * r2_70) * _S1199 + make_float2 (_S1053 * v_70 + (_S1200 + (_S1204 + _S1204)) * _S1048 + _S1200 * _S1049, _S1054 * v_70 + _S1200 * _S1047 + _S1200 * _S1050);
                    float _S1207 = v_70 + v_70;
                    float2  _S1208 = _S1055 * make_float2 (radial_26) + make_float2 (_S1207 * _S1203 + (_S1207 * _S1202 + (_S1207 * _S1201 + _S1207 * _S1046 * r2_70) * r2_70) * r2_70) * _S1199 + make_float2 (_S1053 * u_70 + _S1207 * _S1048 + _S1207 * _S1049, _S1054 * u_70 + (_S1207 + (_S1205 + _S1205)) * _S1047 + _S1207 * _S1050);
                    Matrix<float, 2, 2>  _S1209 = transpose_0(makeMatrix<float, 2, 2> (_S1206 + make_float2 (_S1206.x * _S1051 + _S1206.y * _S1052, 0.0f), _S1208 + make_float2 (_S1208.x * _S1051 + _S1208.y * _S1052, 0.0f)));
                    bool _S1210 = !((F32_min((determinant_0(_S1209)), ((F32_min((_S1209.rows[int(0)].x), (_S1209.rows[int(1)].y)))))) > 0.0f);
                    _S1076 = _S1210;
                    if(_S1210)
                    {
                        break;
                    }
                    float u_71 = proj_points_2[int(6)].x;
                    float v_71 = proj_points_2[int(6)].y;
                    float r2_71 = u_71 * u_71 + v_71 * v_71;
                    float2  _S1211 = proj_points_2[int(6)] * make_float2 (1.0f + r2_71 * (_S1043 + r2_71 * (_S1044 + r2_71 * (_S1045 + r2_71 * _S1046)))) + make_float2 (_S1053 * u_71 * v_71 + _S1048 * (r2_71 + 2.0f * u_71 * u_71) + _S1049 * r2_71, _S1054 * u_71 * v_71 + _S1047 * (r2_71 + 2.0f * v_71 * v_71) + _S1050 * r2_71);
                    float2  _S1212 = _S1211 + make_float2 (_S1051 * _S1211.x + _S1052 * _S1211.y, 0.0f);
                    proj_points_2[int(6)] = make_float2 (fx_10 * _S1212.x + cx_7, fy_10 * _S1212.y + cy_7);
                    break;
                }
                _S1077 = _S1073 & (!_S1076);
                break;
            }
            if(!_S1077)
            {
                _S1079 = false;
                break;
            }
            float2  p_7 = *_S1040 + (*_S1057 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  p_8 = *_S1040 + (*_S1060 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  p_9 = *_S1040 + (*_S1064 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  p_10 = *_S1040 + (*_S1067 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  p_11 = *_S1040 + (*_S1070 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  p_12 = *_S1040 + (*_S1074 - *_S1040) * make_float2 (3.32899999618530273f);
            float2  _S1213 = make_float2 (cx_7, cy_7);
            float2  min_p_1 = min_0(min_0(min_0(min_0(min_0(min_0(*_S1040, p_7), p_8), p_9), p_10), p_11), p_12) - _S1213;
            float2  max_p_1 = max_0(max_0(max_0(max_0(max_0(max_0(*_S1040, p_7), p_8), p_9), p_10), p_11), p_12) - _S1213;
            if((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S1041), (_S1058)))), (_S1061)))), (_S1065)))), (_S1068)))), (_S1071)))), (_S1075))) <= 0.0f)
            {
                _S1079 = (min_p_1.x * max_p_1.x) < 0.0f;
            }
            else
            {
                _S1079 = false;
            }
            if(_S1079)
            {
                _S1079 = (min_p_1.y * max_p_1.y) < 0.0f;
            }
            else
            {
                _S1079 = false;
            }
            if(_S1079)
            {
                _S1079 = false;
                break;
            }
            float2  _S1214 = *mean2d_10 + make_float2 (_S1090.w_mean_0[int(0)]) * *_S1040 + make_float2 (_S1090.w_mean_0[int(1)]) * *_S1057 + make_float2 (_S1090.w_mean_0[int(2)]) * *_S1060 + make_float2 (_S1090.w_mean_0[int(3)]) * *_S1064 + make_float2 (_S1090.w_mean_0[int(4)]) * *_S1067 + make_float2 (_S1090.w_mean_0[int(5)]) * *_S1070 + make_float2 (_S1090.w_mean_0[int(6)]) * *_S1074;
            *mean2d_10 = _S1214;
            float2  d_14 = *_S1040 - _S1214;
            float _S1215 = d_14.x;
            float _S1216 = d_14.y;
            float _S1217 = _S1215 * _S1216;
            float2  d_15 = *_S1057 - _S1214;
            float _S1218 = d_15.x;
            float _S1219 = d_15.y;
            float _S1220 = _S1218 * _S1219;
            float2  d_16 = *_S1060 - _S1214;
            float _S1221 = d_16.x;
            float _S1222 = d_16.y;
            float _S1223 = _S1221 * _S1222;
            float2  d_17 = *_S1064 - _S1214;
            float _S1224 = d_17.x;
            float _S1225 = d_17.y;
            float _S1226 = _S1224 * _S1225;
            float2  d_18 = *_S1067 - _S1214;
            float _S1227 = d_18.x;
            float _S1228 = d_18.y;
            float _S1229 = _S1227 * _S1228;
            float2  d_19 = *_S1070 - _S1214;
            float _S1230 = d_19.x;
            float _S1231 = d_19.y;
            float _S1232 = _S1230 * _S1231;
            float2  d_20 = *_S1074 - _S1214;
            float _S1233 = d_20.x;
            float _S1234 = d_20.y;
            float _S1235 = _S1233 * _S1234;
            covar2d_6 = covar2d_6 + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S1215 * _S1215, _S1217, _S1217, _S1216 * _S1216) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S1218 * _S1218, _S1220, _S1220, _S1219 * _S1219) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S1221 * _S1221, _S1223, _S1223, _S1222 * _S1222) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S1224 * _S1224, _S1226, _S1226, _S1225 * _S1225) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S1227 * _S1227, _S1229, _S1229, _S1228 * _S1228) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S1230 * _S1230, _S1232, _S1232, _S1231 * _S1231) + makeMatrix<float, 2, 2> (_S1090.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S1233 * _S1233, _S1235, _S1235, _S1234 * _S1234);
            _S1079 = true;
            break;
        }
        if(!(true & _S1079))
        {
            *aabb_xyxy_6 = make_float4 (0.0f);
            break;
        }
        float eps2d_6;
        if(antialiased_6)
        {
            eps2d_6 = 0.10000000149011612f;
        }
        else
        {
            eps2d_6 = 0.30000001192092896f;
        }
        float det_orig_6 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
        float _S1236 = *&(((&covar2d_6)->rows + (int(0)))->x) + eps2d_6;
        *&(((&covar2d_6)->rows + (int(0)))->x) = _S1236;
        float _S1237 = *&(((&covar2d_6)->rows + (int(1)))->y) + eps2d_6;
        *&(((&covar2d_6)->rows + (int(1)))->y) = _S1237;
        float det_blur_6 = _S1236 * _S1237 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
        float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / det_blur_6))))));
        if(det_blur_6 <= 0.0f)
        {
            *aabb_xyxy_6 = make_float4 (0.0f);
            break;
        }
        float invdet_8 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
        Matrix<float, 2, 2>  _S1238 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_8, - covar2d_6.rows[int(0)].y * invdet_8, - covar2d_6.rows[int(1)].x * invdet_8, covar2d_6.rows[int(0)].x * invdet_8);
        if(antialiased_6)
        {
            *opacity_6 = *opacity_6 * compensation_6;
        }
        if((*opacity_6) < 0.00392156885936856f)
        {
            *aabb_xyxy_6 = make_float4 (0.0f);
            break;
        }
        float _S1239 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_6 / 0.00392156885936856f)))))))));
        float radius_x_6 = _S1239 * (F32_sqrt((covar2d_6[int(0)].x)));
        float radius_y_6 = _S1239 * (F32_sqrt((covar2d_6[int(1)].y)));
        float _S1240 = (*mean2d_10).x - radius_x_6;
        float _S1241 = (*mean2d_10).x + radius_x_6;
        float _S1242 = (*mean2d_10).y - radius_y_6;
        float _S1243 = (*mean2d_10).y + radius_y_6;
        if(_S1241 <= 0.0f)
        {
            _S1079 = true;
        }
        else
        {
            _S1079 = _S1240 >= float(image_width_6);
        }
        if(_S1079)
        {
            _S1079 = true;
        }
        else
        {
            _S1079 = _S1243 <= 0.0f;
        }
        if(_S1079)
        {
            _S1079 = true;
        }
        else
        {
            _S1079 = _S1242 >= float(image_height_6);
        }
        if(_S1079)
        {
            *aabb_xyxy_6 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_6 = make_float4 (_S1240, _S1242, _S1241, _S1243);
        float x_25 = mean_c_6.x;
        float y_8 = mean_c_6.y;
        float z_1 = mean_c_6.z;
        float _S1244 = x_25 * x_25 + y_8 * y_8;
        *sorting_depth_6 = z_1 * z_1 * z_1 * z_1 + 0.001953125f * _S1244 * _S1244;
        *conic_6 = make_float3 (_S1238.rows[int(0)].x, _S1238.rows[int(0)].y, _S1238.rows[int(1)].y);
        *radius_7 = view_radius_3dgs_0(mean_7, scale_6, in_opacity_6, - mul_6(transpose_3(R_6), t_6));
        break;
    }
    return;
}

inline __device__ void projection_3dgut_equirect(bool antialiased_7, float3  mean_8, float4  quat_7, float3  scale_7, float in_opacity_7, Matrix<float, 3, 3>  R_7, float3  t_7, float fx_11, float fy_11, float cx_8, float cy_8, FixedArray<float, 10>  dist_coeffs_11, uint image_width_7, uint image_height_7, float4  * aabb_xyxy_7, float * sorting_depth_7, float * radius_8, float2  * mean2d_11, float * depth_7, float3  * conic_7, float * opacity_7)
{
    for(;;)
    {
        float3  mean_c_7 = mul_6(R_7, mean_8) + t_7;
        float _S1245 = length_1(mean_c_7);
        *depth_7 = _S1245;
        if(_S1245 <= 0.0f)
        {
            *aabb_xyxy_7 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_7;
        *opacity_7 = 1.0f / (1.0f + (F32_exp((- in_opacity_7))));
        float3  _S1246 = exp_0(scale_7);
        float4  _S1247 = normalize_0(quat_7);
        float x_26 = _S1247.y;
        float x2_7 = x_26 * x_26;
        float y2_7 = _S1247.z * _S1247.z;
        float z2_7 = _S1247.w * _S1247.w;
        float xy_7 = _S1247.y * _S1247.z;
        float xz_7 = _S1247.y * _S1247.w;
        float yz_7 = _S1247.z * _S1247.w;
        float wx_7 = _S1247.x * _S1247.y;
        float wy_7 = _S1247.x * _S1247.z;
        float wz_7 = _S1247.x * _S1247.w;
        Matrix<float, 3, 3>  _S1248 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_7), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))));
        SigmaPoints_0 ret_3;
        (&ret_3)->p_0[int(0)] = mean_8;
        (&ret_3)->w_mean_0[int(0)] = 0.0f;
        (&ret_3)->w_cov_0[int(0)] = 2.0f;
        float _S1249 = (F32_sqrt((3.0f)));
        float3  delta_9 = make_float3 (_S1249 * _S1246.x) * _S1248.rows[0U];
        float3  _S1250 = mean_8 + delta_9;
        float3  _S1251 = mean_8 - delta_9;
        float3  delta_10 = make_float3 (_S1249 * _S1246.y) * _S1248.rows[1U];
        float3  _S1252 = mean_8 + delta_10;
        float3  _S1253 = mean_8 - delta_10;
        float3  delta_11 = make_float3 (_S1249 * _S1246.z) * _S1248.rows[2U];
        float3  _S1254 = mean_8 + delta_11;
        float3  _S1255 = mean_8 - delta_11;
        (&ret_3)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_3)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_3)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_3)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_3)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_3)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_3)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_3)->p_0[0U] = mul_6(R_7, (&ret_3)->p_0[0U]) + t_7;
        (&ret_3)->p_0[1U] = mul_6(R_7, _S1250) + t_7;
        (&ret_3)->p_0[2U] = mul_6(R_7, _S1252) + t_7;
        (&ret_3)->p_0[3U] = mul_6(R_7, _S1254) + t_7;
        (&ret_3)->p_0[4U] = mul_6(R_7, _S1251) + t_7;
        (&ret_3)->p_0[5U] = mul_6(R_7, _S1253) + t_7;
        (&ret_3)->p_0[6U] = mul_6(R_7, _S1255) + t_7;
        FixedArray<float2 , 7>  proj_points_3;
        float _S1256 = fx_11 * (F32_atan2((ret_3.p_0[int(0)].x), (ret_3.p_0[int(0)].z))) + cx_8;
        float2  _S1257 = make_float2 (_S1256, fy_11 * (F32_atan2((ret_3.p_0[int(0)].y), (length_0(float2 {ret_3.p_0[int(0)].x, ret_3.p_0[int(0)].z})))) + cy_8);
        proj_points_3[int(0)] = _S1257;
        float _S1258 = fx_11 * (F32_atan2((ret_3.p_0[int(1)].x), (ret_3.p_0[int(1)].z))) + cx_8;
        proj_points_3[int(1)] = make_float2 (_S1258, fy_11 * (F32_atan2((ret_3.p_0[int(1)].y), (length_0(float2 {ret_3.p_0[int(1)].x, ret_3.p_0[int(1)].z})))) + cy_8);
        float _S1259 = fx_11 * (F32_atan2((ret_3.p_0[int(2)].x), (ret_3.p_0[int(2)].z))) + cx_8;
        proj_points_3[int(2)] = make_float2 (_S1259, fy_11 * (F32_atan2((ret_3.p_0[int(2)].y), (length_0(float2 {ret_3.p_0[int(2)].x, ret_3.p_0[int(2)].z})))) + cy_8);
        float _S1260 = fx_11 * (F32_atan2((ret_3.p_0[int(3)].x), (ret_3.p_0[int(3)].z))) + cx_8;
        proj_points_3[int(3)] = make_float2 (_S1260, fy_11 * (F32_atan2((ret_3.p_0[int(3)].y), (length_0(float2 {ret_3.p_0[int(3)].x, ret_3.p_0[int(3)].z})))) + cy_8);
        float _S1261 = fx_11 * (F32_atan2((ret_3.p_0[int(4)].x), (ret_3.p_0[int(4)].z))) + cx_8;
        proj_points_3[int(4)] = make_float2 (_S1261, fy_11 * (F32_atan2((ret_3.p_0[int(4)].y), (length_0(float2 {ret_3.p_0[int(4)].x, ret_3.p_0[int(4)].z})))) + cy_8);
        float _S1262 = fx_11 * (F32_atan2((ret_3.p_0[int(5)].x), (ret_3.p_0[int(5)].z))) + cx_8;
        proj_points_3[int(5)] = make_float2 (_S1262, fy_11 * (F32_atan2((ret_3.p_0[int(5)].y), (length_0(float2 {ret_3.p_0[int(5)].x, ret_3.p_0[int(5)].z})))) + cy_8);
        float _S1263 = fx_11 * (F32_atan2((ret_3.p_0[int(6)].x), (ret_3.p_0[int(6)].z))) + cx_8;
        proj_points_3[int(6)] = make_float2 (_S1263, fy_11 * (F32_atan2((ret_3.p_0[int(6)].y), (length_0(float2 {ret_3.p_0[int(6)].x, ret_3.p_0[int(6)].z})))) + cy_8);
        float _S1264 = fx_11 * 6.28318548202514648f;
        float du_0 = _S1258 - _S1256;
        *&((&proj_points_3[int(1)])->x) = _S1256 + (du_0 - _S1264 * (F32_round((du_0 / _S1264))));
        float du_1 = _S1259 - _S1256;
        *&((&proj_points_3[int(2)])->x) = _S1256 + (du_1 - _S1264 * (F32_round((du_1 / _S1264))));
        float du_2 = _S1260 - _S1256;
        *&((&proj_points_3[int(3)])->x) = _S1256 + (du_2 - _S1264 * (F32_round((du_2 / _S1264))));
        float du_3 = _S1261 - _S1256;
        *&((&proj_points_3[int(4)])->x) = _S1256 + (du_3 - _S1264 * (F32_round((du_3 / _S1264))));
        float du_4 = _S1262 - _S1256;
        *&((&proj_points_3[int(5)])->x) = _S1256 + (du_4 - _S1264 * (F32_round((du_4 / _S1264))));
        float du_5 = _S1263 - _S1256;
        *&((&proj_points_3[int(6)])->x) = _S1256 + (du_5 - _S1264 * (F32_round((du_5 / _S1264))));
        float2  _S1265 = make_float2 (ret_3.w_mean_0[int(0)]) * _S1257 + make_float2 (ret_3.w_mean_0[int(1)]) * proj_points_3[int(1)] + make_float2 (ret_3.w_mean_0[int(2)]) * proj_points_3[int(2)] + make_float2 (ret_3.w_mean_0[int(3)]) * proj_points_3[int(3)] + make_float2 (ret_3.w_mean_0[int(4)]) * proj_points_3[int(4)] + make_float2 (ret_3.w_mean_0[int(5)]) * proj_points_3[int(5)] + make_float2 (ret_3.w_mean_0[int(6)]) * proj_points_3[int(6)];
        *mean2d_11 = _S1265;
        float2  d_21 = _S1257 - _S1265;
        float _S1266 = d_21.x;
        float _S1267 = d_21.y;
        float _S1268 = _S1266 * _S1267;
        float2  d_22 = proj_points_3[int(1)] - _S1265;
        float _S1269 = d_22.x;
        float _S1270 = d_22.y;
        float _S1271 = _S1269 * _S1270;
        float2  d_23 = proj_points_3[int(2)] - _S1265;
        float _S1272 = d_23.x;
        float _S1273 = d_23.y;
        float _S1274 = _S1272 * _S1273;
        float2  d_24 = proj_points_3[int(3)] - _S1265;
        float _S1275 = d_24.x;
        float _S1276 = d_24.y;
        float _S1277 = _S1275 * _S1276;
        float2  d_25 = proj_points_3[int(4)] - _S1265;
        float _S1278 = d_25.x;
        float _S1279 = d_25.y;
        float _S1280 = _S1278 * _S1279;
        float2  d_26 = proj_points_3[int(5)] - _S1265;
        float _S1281 = d_26.x;
        float _S1282 = d_26.y;
        float _S1283 = _S1281 * _S1282;
        float2  d_27 = proj_points_3[int(6)] - _S1265;
        float _S1284 = d_27.x;
        float _S1285 = d_27.y;
        float _S1286 = _S1284 * _S1285;
        covar2d_7 = makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S1266 * _S1266, _S1268, _S1268, _S1267 * _S1267) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S1269 * _S1269, _S1271, _S1271, _S1270 * _S1270) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S1272 * _S1272, _S1274, _S1274, _S1273 * _S1273) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S1275 * _S1275, _S1277, _S1277, _S1276 * _S1276) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S1278 * _S1278, _S1280, _S1280, _S1279 * _S1279) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S1281 * _S1281, _S1283, _S1283, _S1282 * _S1282) + makeMatrix<float, 2, 2> (ret_3.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S1284 * _S1284, _S1286, _S1286, _S1285 * _S1285);
        float eps2d_7;
        if(antialiased_7)
        {
            eps2d_7 = 0.10000000149011612f;
        }
        else
        {
            eps2d_7 = 0.30000001192092896f;
        }
        float det_orig_7 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
        float _S1287 = *&(((&covar2d_7)->rows + (int(0)))->x) + eps2d_7;
        *&(((&covar2d_7)->rows + (int(0)))->x) = _S1287;
        float _S1288 = *&(((&covar2d_7)->rows + (int(1)))->y) + eps2d_7;
        *&(((&covar2d_7)->rows + (int(1)))->y) = _S1288;
        float det_blur_7 = _S1287 * _S1288 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
        float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / det_blur_7))))));
        if(det_blur_7 <= 0.0f)
        {
            *aabb_xyxy_7 = make_float4 (0.0f);
            break;
        }
        float invdet_9 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
        Matrix<float, 2, 2>  _S1289 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_9, - covar2d_7.rows[int(0)].y * invdet_9, - covar2d_7.rows[int(1)].x * invdet_9, covar2d_7.rows[int(0)].x * invdet_9);
        if(antialiased_7)
        {
            *opacity_7 = *opacity_7 * compensation_7;
        }
        if((*opacity_7) < 0.00392156885936856f)
        {
            *aabb_xyxy_7 = make_float4 (0.0f);
            break;
        }
        float _S1290 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_7 / 0.00392156885936856f)))))))));
        float radius_x_7 = _S1290 * (F32_sqrt((covar2d_7[int(0)].x)));
        float radius_y_7 = _S1290 * (F32_sqrt((covar2d_7[int(1)].y)));
        float _S1291 = (*mean2d_11).x - radius_x_7;
        float _S1292 = (*mean2d_11).x + radius_x_7;
        float _S1293 = (*mean2d_11).y - radius_y_7;
        float _S1294 = (*mean2d_11).y + radius_y_7;
        bool _S1295;
        if(_S1292 <= 0.0f)
        {
            _S1295 = true;
        }
        else
        {
            _S1295 = _S1291 >= float(image_width_7);
        }
        if(_S1295)
        {
            _S1295 = true;
        }
        else
        {
            _S1295 = _S1294 <= 0.0f;
        }
        if(_S1295)
        {
            _S1295 = true;
        }
        else
        {
            _S1295 = _S1293 >= float(image_height_7);
        }
        if(_S1295)
        {
            *aabb_xyxy_7 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_7 = make_float4 (_S1291, _S1293, _S1292, _S1294);
        *sorting_depth_7 = dot_1(mean_c_7, mean_c_7);
        *conic_7 = make_float3 (_S1289.rows[int(0)].x, _S1289.rows[int(0)].y, _S1289.rows[int(1)].y);
        *radius_8 = view_radius_3dgs_0(mean_8, scale_7, in_opacity_7, - mul_6(transpose_3(R_7), t_7));
        break;
    }
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S1296, float3  _S1297)
{
    return mul_6(_S1296, _S1297);
}

inline __device__ float s_primal_ctx_exp_0(float _S1298)
{
    return (F32_exp((_S1298)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S1299)
{
    return exp_0(_S1299);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1300, Matrix<float, 3, 3>  _S1301)
{
    return mul_5(_S1300, _S1301);
}

inline __device__ float s_primal_ctx_clamp_0(float _S1302, float _S1303, float _S1304)
{
    return clamp_0(_S1302, _S1303, _S1304);
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S1305, Matrix<float, 3, 3>  _S1306)
{
    return mul_3(_S1305, _S1306);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S1307, Matrix<float, 3, 2>  _S1308)
{
    return mul_4(_S1307, _S1308);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S1309)
{
    return (F32_sqrt((_S1309)));
}

inline __device__ float s_primal_ctx_log_0(float _S1310)
{
    return (F32_log((_S1310)));
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1311, float _S1312)
{
    _d_sqrt_0(_S1311, _S1312);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float _s_dOut_0)
{
    float _S1313 = (*dpx_15).primal_0.x;
    float _S1314 = (*dpx_15).primal_0.y;
    float _S1315 = (*dpx_15).primal_0.z;
    DiffPair_float_0 _S1316;
    (&_S1316)->primal_0 = _S1313 * _S1313 + _S1314 * _S1314 + _S1315 * _S1315;
    (&_S1316)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1316, _s_dOut_0);
    float _S1317 = (*dpx_15).primal_0.z * _S1316.differential_0;
    float _S1318 = _S1317 + _S1317;
    float _S1319 = (*dpx_15).primal_0.y * _S1316.differential_0;
    float _S1320 = _S1319 + _S1319;
    float _S1321 = (*dpx_15).primal_0.x * _S1316.differential_0;
    float _S1322 = _S1321 + _S1321;
    float3  _S1323 = make_float3 (0.0f);
    *&((&_S1323)->z) = _S1318;
    *&((&_S1323)->y) = _S1320;
    *&((&_S1323)->x) = _S1322;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S1323;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1324, float _S1325)
{
    s_bwd_prop_length_impl_0(_S1324, _S1325);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1326, float _S1327)
{
    _d_exp_0(_S1326, _S1327);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S1328, float _S1329)
{
    _d_log_0(_S1328, _S1329);
    return;
}

inline __device__ void s_bwd_prop_view_radius_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dplog_scale_0, DiffPair_float_0 * dplogit_opacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpcampos_0, float _s_dOut_1)
{
    float _S1330 = - (*dplogit_opacity_0).primal_0;
    float _S1331 = 1.0f + s_primal_ctx_exp_0(_S1330);
    float _S1332 = 255.0f / _S1331;
    float _S1333 = _S1331 * _S1331;
    float _S1334 = (F32_max((_S1332), (1.0f)));
    float _S1335 = 2.0f * s_primal_ctx_log_0(_S1334);
    float _S1336 = s_primal_ctx_sqrt_0(_S1335);
    float _S1337 = (*dplog_scale_0).primal_0.x;
    float _S1338 = (*dplog_scale_0).primal_0.y;
    float _S1339 = (*dplog_scale_0).primal_0.z;
    float _S1340 = (F32_max((_S1338), (_S1339)));
    float _S1341 = (F32_max((_S1337), (_S1340)));
    float _S1342 = s_primal_ctx_exp_0(_S1341);
    float radius_9 = _S1342 * _S1336;
    float3  _S1343 = (*dpmean_0).primal_0 - (*dpcampos_0).primal_0;
    float _S1344 = length_1(_S1343);
    float _S1345 = _S1344 * _S1344 - radius_9 * radius_9;
    float _S1346 = (F32_max((_S1345), (0.0f)));
    float _S1347 = (F32_max((_S1344), (radius_9))) + s_primal_ctx_sqrt_0(_S1346);
    float _S1348 = _s_dOut_1 / (_S1347 * _S1347);
    float _S1349 = radius_9 * - _S1348;
    float _S1350 = _S1347 * _S1348;
    DiffPair_float_0 _S1351;
    (&_S1351)->primal_0 = _S1346;
    (&_S1351)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1351, _S1349);
    DiffPair_float_0 _S1352;
    (&_S1352)->primal_0 = _S1345;
    (&_S1352)->differential_0 = 0.0f;
    DiffPair_float_0 _S1353;
    (&_S1353)->primal_0 = 0.0f;
    (&_S1353)->differential_0 = 0.0f;
    _d_max_0(&_S1352, &_S1353, _S1351.differential_0);
    float _S1354 = radius_9 * - _S1352.differential_0;
    float _S1355 = _S1344 * _S1352.differential_0;
    DiffPair_float_0 _S1356;
    (&_S1356)->primal_0 = _S1344;
    (&_S1356)->differential_0 = 0.0f;
    DiffPair_float_0 _S1357;
    (&_S1357)->primal_0 = radius_9;
    (&_S1357)->differential_0 = 0.0f;
    _d_max_0(&_S1356, &_S1357, _S1349);
    float _S1358 = _S1355 + _S1355 + _S1356.differential_0;
    float3  _S1359 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1360;
    (&_S1360)->primal_0 = _S1343;
    (&_S1360)->differential_0 = _S1359;
    s_bwd_length_impl_0(&_S1360, _S1358);
    float3  _S1361 = - _S1360.differential_0;
    float _S1362 = _S1350 + _S1354 + _S1354 + _S1357.differential_0;
    float _S1363 = _S1342 * _S1362;
    float _S1364 = _S1336 * _S1362;
    DiffPair_float_0 _S1365;
    (&_S1365)->primal_0 = _S1341;
    (&_S1365)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1365, _S1364);
    DiffPair_float_0 _S1366;
    (&_S1366)->primal_0 = _S1337;
    (&_S1366)->differential_0 = 0.0f;
    DiffPair_float_0 _S1367;
    (&_S1367)->primal_0 = _S1340;
    (&_S1367)->differential_0 = 0.0f;
    _d_max_0(&_S1366, &_S1367, _S1365.differential_0);
    DiffPair_float_0 _S1368;
    (&_S1368)->primal_0 = _S1338;
    (&_S1368)->differential_0 = 0.0f;
    DiffPair_float_0 _S1369;
    (&_S1369)->primal_0 = _S1339;
    (&_S1369)->differential_0 = 0.0f;
    _d_max_0(&_S1368, &_S1369, _S1367.differential_0);
    DiffPair_float_0 _S1370;
    (&_S1370)->primal_0 = _S1335;
    (&_S1370)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1370, _S1363);
    float _S1371 = 2.0f * _S1370.differential_0;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1334;
    (&_S1372)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1372, _S1371);
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S1332;
    (&_S1373)->differential_0 = 0.0f;
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = 1.0f;
    (&_S1374)->differential_0 = 0.0f;
    _d_max_0(&_S1373, &_S1374, _S1372.differential_0);
    float _S1375 = 255.0f * - (_S1373.differential_0 / _S1333);
    DiffPair_float_0 _S1376;
    (&_S1376)->primal_0 = _S1330;
    (&_S1376)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1376, _S1375);
    float _S1377 = - _S1376.differential_0;
    dpcampos_0->primal_0 = (*dpcampos_0).primal_0;
    dpcampos_0->differential_0 = _S1361;
    dplogit_opacity_0->primal_0 = (*dplogit_opacity_0).primal_0;
    dplogit_opacity_0->differential_0 = _S1377;
    float3  _S1378 = make_float3 (_S1366.differential_0, _S1368.differential_0, _S1369.differential_0);
    dplog_scale_0->primal_0 = (*dplog_scale_0).primal_0;
    dplog_scale_0->differential_0 = _S1378;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S1360.differential_0;
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S1379, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S1380, Matrix<float, 2, 2>  _S1381)
{
    mul_1(_S1379, _S1380, _S1381);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S1382, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1383, Matrix<float, 2, 3>  _S1384)
{
    mul_0(_S1382, _S1383, _S1384);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1385, DiffPair_float_0 * _S1386, DiffPair_float_0 * _S1387, float _S1388)
{
    _d_clamp_0(_S1385, _S1386, _S1387, _S1388);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1389, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1390, Matrix<float, 3, 3>  _S1391)
{
    mul_2(_S1389, _S1390, _S1391);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1392, float3  _S1393)
{
    _d_exp_vector_0(_S1392, _S1393);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_16, float _s_dOut_2)
{
    float _S1394 = (*dpx_16).primal_0.x;
    float _S1395 = (*dpx_16).primal_0.y;
    float _S1396 = (*dpx_16).primal_0.z;
    float _S1397 = (*dpx_16).primal_0.w;
    DiffPair_float_0 _S1398;
    (&_S1398)->primal_0 = _S1394 * _S1394 + _S1395 * _S1395 + _S1396 * _S1396 + _S1397 * _S1397;
    (&_S1398)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1398, _s_dOut_2);
    float _S1399 = (*dpx_16).primal_0.w * _S1398.differential_0;
    float _S1400 = _S1399 + _S1399;
    float _S1401 = (*dpx_16).primal_0.z * _S1398.differential_0;
    float _S1402 = _S1401 + _S1401;
    float _S1403 = (*dpx_16).primal_0.y * _S1398.differential_0;
    float _S1404 = _S1403 + _S1403;
    float _S1405 = (*dpx_16).primal_0.x * _S1398.differential_0;
    float _S1406 = _S1405 + _S1405;
    float4  _S1407 = make_float4 (0.0f);
    *&((&_S1407)->w) = _S1400;
    *&((&_S1407)->z) = _S1402;
    *&((&_S1407)->y) = _S1404;
    *&((&_S1407)->x) = _S1406;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S1407;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1408, float _S1409)
{
    s_bwd_prop_length_impl_1(_S1408, _S1409);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_17, float4  _s_dOut_3)
{
    float _S1410 = length_2((*dpx_17).primal_0);
    float4  _S1411 = (*dpx_17).primal_0 * _s_dOut_3;
    float4  _S1412 = make_float4 (1.0f / _S1410) * _s_dOut_3;
    float _S1413 = - ((_S1411.x + _S1411.y + _S1411.z + _S1411.w) / (_S1410 * _S1410));
    float4  _S1414 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1415;
    (&_S1415)->primal_0 = (*dpx_17).primal_0;
    (&_S1415)->differential_0 = _S1414;
    s_bwd_length_impl_1(&_S1415, _S1413);
    float4  _S1416 = _S1412 + _S1415.differential_0;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S1416;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1417, float4  _S1418)
{
    s_bwd_prop_normalize_impl_0(_S1417, _S1418);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1419, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1420, float3  _S1421)
{
    _d_mul_0(_S1419, _S1420, _S1421);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_8, float3  mean_9, float4  quat_8, float3  scale_8, float in_opacity_8, Matrix<float, 3, 3>  R_8, float3  t_8, float fx_12, float fy_12, float cx_9, float cy_9, FixedArray<float, 10>  dist_coeffs_12, uint image_width_8, uint image_height_8, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_8 = s_primal_ctx_mul_0(R_8, mean_9) + t_8;
    float _S1422 = - in_opacity_8;
    float _S1423 = 1.0f + s_primal_ctx_exp_0(_S1422);
    float _S1424 = 1.0f / _S1423;
    float _S1425 = _S1423 * _S1423;
    float4  _S1426 = normalize_0(quat_8);
    float3  _S1427 = s_primal_ctx_exp_1(scale_8);
    float _S1428 = _S1426.y;
    float x2_8 = _S1428 * _S1428;
    float y2_8 = _S1426.z * _S1426.z;
    float z2_8 = _S1426.w * _S1426.w;
    float xy_8 = _S1426.y * _S1426.z;
    float xz_8 = _S1426.y * _S1426.w;
    float yz_8 = _S1426.z * _S1426.w;
    float wx_8 = _S1426.x * _S1426.y;
    float wy_8 = _S1426.x * _S1426.z;
    float wz_8 = _S1426.x * _S1426.w;
    Matrix<float, 3, 3>  _S1429 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_8), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_8), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S1427.x, 0.0f, 0.0f, 0.0f, _S1427.y, 0.0f, 0.0f, 0.0f, _S1427.z);
    Matrix<float, 3, 3>  _S1430 = s_primal_ctx_mul_1(_S1429, S_0);
    Matrix<float, 3, 3>  _S1431 = transpose_3(_S1430);
    Matrix<float, 3, 3>  _S1432 = s_primal_ctx_mul_1(_S1430, _S1431);
    Matrix<float, 3, 3>  _S1433 = s_primal_ctx_mul_1(R_8, _S1432);
    Matrix<float, 3, 3>  _S1434 = transpose_3(R_8);
    Matrix<float, 3, 3>  _S1435 = s_primal_ctx_mul_1(_S1433, _S1434);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S1436 = float2 {mean_c_8.x, mean_c_8.y};
    float2  _S1437 = make_float2 (1.0f, 0.0f);
    float _S1438 = mean_c_8.z;
    float2  _S1439 = make_float2 (_S1438);
    float2  _S1440 = _S1436 / make_float2 (_S1438);
    float _S1441 = _S1438 * _S1438;
    float2  _S1442 = make_float2 (_S1441);
    float2  _S1443 = _S1437 * make_float2 (_S1438);
    float2  _S1444 = _S1443 / make_float2 (_S1441);
    float2  _S1445 = make_float2 (_S1441 * _S1441);
    float u_72 = _S1440.x;
    float s_diff_u_18 = _S1444.x;
    float v_72 = _S1440.y;
    float s_diff_v_18 = _S1444.y;
    float _S1446 = s_diff_u_18 * u_72;
    float _S1447 = s_diff_v_18 * v_72;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float s_diff_r2_18 = _S1446 + _S1446 + (_S1447 + _S1447);
    float _S1448 = s_diff_r2_18 * dist_coeffs_12[int(3)];
    float _S1449 = dist_coeffs_12[int(2)] + r2_72 * dist_coeffs_12[int(3)];
    float _S1450 = s_diff_r2_18 * _S1449 + _S1448 * r2_72;
    float _S1451 = dist_coeffs_12[int(1)] + r2_72 * _S1449;
    float _S1452 = s_diff_r2_18 * _S1451 + _S1450 * r2_72;
    float _S1453 = dist_coeffs_12[int(0)] + r2_72 * _S1451;
    float _S1454 = s_diff_r2_18 * _S1453 + _S1452 * r2_72;
    float2  _S1455 = make_float2 (_S1454);
    float radial_27 = 1.0f + r2_72 * _S1453;
    float2  _S1456 = make_float2 (radial_27);
    float _S1457 = 2.0f * dist_coeffs_12[int(4)];
    float _S1458 = _S1457 * u_72;
    float _S1459 = s_diff_u_18 * _S1457;
    float _S1460 = 2.0f * u_72;
    float _S1461 = s_diff_u_18 * 2.0f;
    float _S1462 = 2.0f * dist_coeffs_12[int(5)];
    float _S1463 = _S1462 * u_72;
    float _S1464 = s_diff_u_18 * _S1462;
    float _S1465 = 2.0f * v_72;
    float _S1466 = s_diff_v_18 * 2.0f;
    float2  _S1467 = _S1444 * make_float2 (radial_27) + make_float2 (_S1454) * _S1440 + make_float2 (_S1459 * v_72 + s_diff_v_18 * _S1458 + (s_diff_r2_18 + (_S1461 * u_72 + s_diff_u_18 * _S1460)) * dist_coeffs_12[int(5)] + s_diff_r2_18 * dist_coeffs_12[int(6)], _S1464 * v_72 + s_diff_v_18 * _S1463 + (s_diff_r2_18 + (_S1466 * v_72 + s_diff_v_18 * _S1465)) * dist_coeffs_12[int(4)] + s_diff_r2_18 * dist_coeffs_12[int(7)]);
    float2  _S1468 = _S1467 + make_float2 (_S1467.x * dist_coeffs_12[int(8)] + _S1467.y * dist_coeffs_12[int(9)], 0.0f);
    float _S1469 = _S1468.x * fx_12;
    float _S1470 = _S1468.y * fy_12;
    Matrix<float, 2, 3>  _S1471 = J_10;
    *&(((&_S1471)->rows + (int(0)))->x) = _S1469;
    *&(((&_S1471)->rows + (int(1)))->x) = _S1470;
    float2  _S1472 = make_float2 (0.0f, 1.0f);
    float2  _S1473 = _S1436 / make_float2 (_S1438);
    float2  _S1474 = _S1472 * make_float2 (_S1438);
    float2  _S1475 = _S1474 / make_float2 (_S1441);
    float u_73 = _S1473.x;
    float s_diff_u_19 = _S1475.x;
    float v_73 = _S1473.y;
    float s_diff_v_19 = _S1475.y;
    float _S1476 = s_diff_u_19 * u_73;
    float _S1477 = s_diff_v_19 * v_73;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float s_diff_r2_19 = _S1476 + _S1476 + (_S1477 + _S1477);
    float _S1478 = s_diff_r2_19 * dist_coeffs_12[int(3)];
    float _S1479 = dist_coeffs_12[int(2)] + r2_73 * dist_coeffs_12[int(3)];
    float _S1480 = s_diff_r2_19 * _S1479 + _S1478 * r2_73;
    float _S1481 = dist_coeffs_12[int(1)] + r2_73 * _S1479;
    float _S1482 = s_diff_r2_19 * _S1481 + _S1480 * r2_73;
    float _S1483 = dist_coeffs_12[int(0)] + r2_73 * _S1481;
    float _S1484 = s_diff_r2_19 * _S1483 + _S1482 * r2_73;
    float2  _S1485 = make_float2 (_S1484);
    float radial_28 = 1.0f + r2_73 * _S1483;
    float2  _S1486 = make_float2 (radial_28);
    float _S1487 = _S1457 * u_73;
    float _S1488 = s_diff_u_19 * _S1457;
    float _S1489 = 2.0f * u_73;
    float _S1490 = s_diff_u_19 * 2.0f;
    float _S1491 = _S1462 * u_73;
    float _S1492 = s_diff_u_19 * _S1462;
    float _S1493 = 2.0f * v_73;
    float _S1494 = s_diff_v_19 * 2.0f;
    float2  _S1495 = _S1475 * make_float2 (radial_28) + make_float2 (_S1484) * _S1473 + make_float2 (_S1488 * v_73 + s_diff_v_19 * _S1487 + (s_diff_r2_19 + (_S1490 * u_73 + s_diff_u_19 * _S1489)) * dist_coeffs_12[int(5)] + s_diff_r2_19 * dist_coeffs_12[int(6)], _S1492 * v_73 + s_diff_v_19 * _S1491 + (s_diff_r2_19 + (_S1494 * v_73 + s_diff_v_19 * _S1493)) * dist_coeffs_12[int(4)] + s_diff_r2_19 * dist_coeffs_12[int(7)]);
    float2  _S1496 = _S1495 + make_float2 (_S1495.x * dist_coeffs_12[int(8)] + _S1495.y * dist_coeffs_12[int(9)], 0.0f);
    float _S1497 = _S1496.y * fy_12;
    *&(((&_S1471)->rows + (int(0)))->y) = _S1496.x * fx_12;
    *&(((&_S1471)->rows + (int(1)))->y) = _S1497;
    float2  _S1498 = _S1436 / make_float2 (_S1438);
    float2  _S1499 = make_float2 (0.0f, 0.0f) - _S1436;
    float2  _S1500 = _S1499 / make_float2 (_S1441);
    float u_74 = _S1498.x;
    float s_diff_u_20 = _S1500.x;
    float v_74 = _S1498.y;
    float s_diff_v_20 = _S1500.y;
    float _S1501 = s_diff_u_20 * u_74;
    float _S1502 = s_diff_v_20 * v_74;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float s_diff_r2_20 = _S1501 + _S1501 + (_S1502 + _S1502);
    float _S1503 = s_diff_r2_20 * dist_coeffs_12[int(3)];
    float _S1504 = dist_coeffs_12[int(2)] + r2_74 * dist_coeffs_12[int(3)];
    float _S1505 = s_diff_r2_20 * _S1504 + _S1503 * r2_74;
    float _S1506 = dist_coeffs_12[int(1)] + r2_74 * _S1504;
    float _S1507 = s_diff_r2_20 * _S1506 + _S1505 * r2_74;
    float _S1508 = dist_coeffs_12[int(0)] + r2_74 * _S1506;
    float _S1509 = s_diff_r2_20 * _S1508 + _S1507 * r2_74;
    float2  _S1510 = make_float2 (_S1509);
    float radial_29 = 1.0f + r2_74 * _S1508;
    float2  _S1511 = make_float2 (radial_29);
    float _S1512 = _S1457 * u_74;
    float _S1513 = s_diff_u_20 * _S1457;
    float _S1514 = 2.0f * u_74;
    float _S1515 = s_diff_u_20 * 2.0f;
    float _S1516 = _S1462 * u_74;
    float _S1517 = s_diff_u_20 * _S1462;
    float _S1518 = 2.0f * v_74;
    float _S1519 = s_diff_v_20 * 2.0f;
    float2  _S1520 = _S1500 * make_float2 (radial_29) + make_float2 (_S1509) * _S1498 + make_float2 (_S1513 * v_74 + s_diff_v_20 * _S1512 + (s_diff_r2_20 + (_S1515 * u_74 + s_diff_u_20 * _S1514)) * dist_coeffs_12[int(5)] + s_diff_r2_20 * dist_coeffs_12[int(6)], _S1517 * v_74 + s_diff_v_20 * _S1516 + (s_diff_r2_20 + (_S1519 * v_74 + s_diff_v_20 * _S1518)) * dist_coeffs_12[int(4)] + s_diff_r2_20 * dist_coeffs_12[int(7)]);
    float2  _S1521 = _S1520 + make_float2 (_S1520.x * dist_coeffs_12[int(8)] + _S1520.y * dist_coeffs_12[int(9)], 0.0f);
    float _S1522 = _S1521.x * fx_12;
    float _S1523 = _S1521.y * fy_12;
    float _S1524 = float(image_width_8);
    float _S1525 = 0.30000001192092896f * (0.5f * _S1524);
    float lim_x_pos_3 = _S1524 + _S1525;
    float rz_2 = 1.0f / _S1438;
    float _S1526 = - _S1525;
    float _S1527 = - (_S1526 - cx_9);
    float max_Jxz_0 = _S1527 * rz_2;
    float _S1528 = - (lim_x_pos_3 - cx_9);
    float min_Jxz_0 = _S1528 * rz_2;
    float _S1529 = - (_S1526 - cy_9);
    float max_Jyz_2 = _S1529 * rz_2;
    float _S1530 = - (lim_x_pos_3 - cy_9);
    float min_Jyz_2 = _S1530 * rz_2;
    *&(((&_S1471)->rows + (int(0)))->z) = s_primal_ctx_clamp_0(_S1522, min_Jxz_0, max_Jxz_0);
    *&(((&_S1471)->rows + (int(1)))->z) = s_primal_ctx_clamp_0(_S1523, min_Jyz_2, max_Jyz_2);
    Matrix<float, 2, 3>  _S1531 = s_primal_ctx_mul_2(_S1471, _S1435);
    Matrix<float, 3, 2>  _S1532 = transpose_1(_S1471);
    Matrix<float, 2, 2>  _S1533 = s_primal_ctx_mul_3(_S1531, _S1532);
    float eps2d_8;
    if(antialiased_8)
    {
        eps2d_8 = 0.10000000149011612f;
    }
    else
    {
        eps2d_8 = 0.30000001192092896f;
    }
    float _S1534 = _S1533.rows[int(0)].y * _S1533.rows[int(1)].x;
    float det_orig_8 = _S1533.rows[int(0)].x * _S1533.rows[int(1)].y - _S1534;
    float _S1535 = _S1533.rows[int(0)].x + eps2d_8;
    Matrix<float, 2, 2>  _S1536 = _S1533;
    *&(((&_S1536)->rows + (int(0)))->x) = _S1535;
    float _S1537 = _S1533.rows[int(1)].y + eps2d_8;
    *&(((&_S1536)->rows + (int(1)))->y) = _S1537;
    Matrix<float, 2, 2>  _S1538 = _S1536;
    Matrix<float, 2, 2>  _S1539 = _S1536;
    float det_blur_8 = _S1535 * _S1537 - _S1534;
    float _S1540 = det_orig_8 / det_blur_8;
    float _S1541 = det_blur_8 * det_blur_8;
    float _S1542 = (F32_max((0.0f), (_S1540)));
    float _S1543 = s_primal_ctx_sqrt_0(_S1542);
    float invdet_10 = 1.0f / det_blur_8;
    float _S1544 = - _S1533.rows[int(0)].y;
    float _S1545 = - _S1533.rows[int(1)].x;
    if(antialiased_8)
    {
        eps2d_8 = _S1424 * _S1543;
    }
    else
    {
        eps2d_8 = _S1424;
    }
    float _S1546 = eps2d_8 / 0.00392156885936856f;
    float _S1547 = 2.0f * s_primal_ctx_log_0(_S1546);
    float _S1548 = s_primal_ctx_sqrt_0(_S1547);
    float _S1549 = _S1538.rows[int(0)].x;
    float _S1550 = _S1539.rows[int(1)].y;
    float3  campos_1 = - s_primal_ctx_mul_0(_S1434, t_8);
    float3  _S1551 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1552;
    (&_S1552)->primal_0 = mean_9;
    (&_S1552)->differential_0 = _S1551;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1553;
    (&_S1553)->primal_0 = scale_8;
    (&_S1553)->differential_0 = _S1551;
    DiffPair_float_0 _S1554;
    (&_S1554)->primal_0 = in_opacity_8;
    (&_S1554)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1555;
    (&_S1555)->primal_0 = campos_1;
    (&_S1555)->differential_0 = _S1551;
    s_bwd_prop_view_radius_3dgs_0(&_S1552, &_S1553, &_S1554, &_S1555, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1556 = _S1552;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1557 = _S1553;
    DiffPair_float_0 _S1558 = _S1554;
    float2  _S1559 = make_float2 (0.0f);
    float2  _S1560 = _S1559;
    *&((&_S1560)->y) = v_conic_0.z;
    float2  _S1561 = _S1559;
    *&((&_S1561)->y) = v_conic_0.y;
    *&((&_S1561)->x) = v_conic_0.x;
    DiffPair_float_0 _S1562;
    (&_S1562)->primal_0 = _S1550;
    (&_S1562)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1562, 0.0f);
    DiffPair_float_0 _S1563;
    (&_S1563)->primal_0 = _S1549;
    (&_S1563)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1563, 0.0f);
    DiffPair_float_0 _S1564;
    (&_S1564)->primal_0 = 3.32999992370605469f;
    (&_S1564)->differential_0 = 0.0f;
    DiffPair_float_0 _S1565;
    (&_S1565)->primal_0 = _S1548;
    (&_S1565)->differential_0 = 0.0f;
    _d_min_0(&_S1564, &_S1565, 0.0f);
    DiffPair_float_0 _S1566;
    (&_S1566)->primal_0 = _S1547;
    (&_S1566)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1566, _S1565.differential_0);
    float _S1567 = 2.0f * _S1566.differential_0;
    DiffPair_float_0 _S1568;
    (&_S1568)->primal_0 = _S1546;
    (&_S1568)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1568, _S1567);
    float _S1569 = v_opacity_0 + 254.9999847412109375f * _S1568.differential_0;
    float2  _S1570 = make_float2 (_S1563.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S1571 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1572 = _S1571;
    _S1572[int(1)] = _S1560;
    _S1572[int(0)] = _S1561;
    Matrix<float, 2, 2>  _S1573 = _S1572;
    float2  _S1574 = make_float2 (0.0f, _S1562.differential_0);
    float _S1575;
    if(antialiased_8)
    {
        float _S1576 = _S1543 * _S1569;
        eps2d_8 = _S1424 * _S1569;
        _S1575 = _S1576;
    }
    else
    {
        eps2d_8 = 0.0f;
        _S1575 = _S1569;
    }
    float _S1577 = invdet_10 * _S1573.rows[int(1)].y;
    float _S1578 = - (invdet_10 * _S1573.rows[int(1)].x);
    float _S1579 = - (invdet_10 * _S1573.rows[int(0)].y);
    float _S1580 = invdet_10 * _S1573.rows[int(0)].x;
    float _S1581 = - ((_S1535 * _S1573.rows[int(1)].y + _S1545 * _S1573.rows[int(1)].x + _S1544 * _S1573.rows[int(0)].y + _S1537 * _S1573.rows[int(0)].x) / _S1541);
    DiffPair_float_0 _S1582;
    (&_S1582)->primal_0 = _S1542;
    (&_S1582)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1582, eps2d_8);
    DiffPair_float_0 _S1583;
    (&_S1583)->primal_0 = 0.0f;
    (&_S1583)->differential_0 = 0.0f;
    DiffPair_float_0 _S1584;
    (&_S1584)->primal_0 = _S1540;
    (&_S1584)->differential_0 = 0.0f;
    _d_max_0(&_S1583, &_S1584, _S1582.differential_0);
    float _S1585 = _S1584.differential_0 / _S1541;
    float s_diff_det_orig_T_0 = det_blur_8 * _S1585;
    float _S1586 = det_orig_8 * - _S1585 + _S1581;
    float _S1587 = - _S1586;
    float _S1588 = _S1535 * _S1586;
    float _S1589 = _S1537 * _S1586;
    Matrix<float, 2, 2>  _S1590 = _S1571;
    _S1590[int(1)] = _S1574;
    _S1590[int(0)] = _S1570;
    _S1536 = _S1590;
    *&(((&_S1536)->rows + (int(1)))->y) = 0.0f;
    float _S1591 = _S1588 + _S1590.rows[int(1)].y + _S1580;
    *&(((&_S1536)->rows + (int(0)))->x) = 0.0f;
    float _S1592 = _S1589 + _S1590.rows[int(0)].x + _S1577;
    float _S1593 = _S1587 + - s_diff_det_orig_T_0;
    float _S1594 = _S1533.rows[int(0)].y * _S1593 + _S1578;
    float _S1595 = _S1533.rows[int(1)].x * _S1593 + _S1579;
    float _S1596 = _S1533.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1597 = _S1591 + _S1533.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1598 = _S1559;
    *&((&_S1598)->x) = _S1594;
    *&((&_S1598)->y) = _S1597;
    float _S1599 = _S1592 + _S1596;
    float2  _S1600 = _S1559;
    *&((&_S1600)->y) = _S1595;
    *&((&_S1600)->x) = _S1599;
    Matrix<float, 2, 2>  _S1601 = _S1571;
    _S1601[int(1)] = _S1598;
    _S1601[int(0)] = _S1600;
    Matrix<float, 2, 2>  _S1602 = _S1536 + _S1601;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1603;
    (&_S1603)->primal_0 = _S1531;
    (&_S1603)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S1604 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1605;
    (&_S1605)->primal_0 = _S1532;
    (&_S1605)->differential_0 = _S1604;
    s_bwd_prop_mul_0(&_S1603, &_S1605, _S1602);
    Matrix<float, 2, 3>  _S1606 = transpose_2(_S1605.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1607;
    (&_S1607)->primal_0 = _S1471;
    (&_S1607)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S1608 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1609;
    (&_S1609)->primal_0 = _S1435;
    (&_S1609)->differential_0 = _S1608;
    s_bwd_prop_mul_1(&_S1607, &_S1609, _S1603.differential_0);
    Matrix<float, 2, 3>  _S1610 = _S1606 + _S1607.differential_0;
    DiffPair_float_0 _S1611;
    (&_S1611)->primal_0 = _S1523;
    (&_S1611)->differential_0 = 0.0f;
    DiffPair_float_0 _S1612;
    (&_S1612)->primal_0 = min_Jyz_2;
    (&_S1612)->differential_0 = 0.0f;
    DiffPair_float_0 _S1613;
    (&_S1613)->primal_0 = max_Jyz_2;
    (&_S1613)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1611, &_S1612, &_S1613, _S1610.rows[int(1)].z);
    DiffPair_float_0 _S1614;
    (&_S1614)->primal_0 = _S1522;
    (&_S1614)->differential_0 = 0.0f;
    DiffPair_float_0 _S1615;
    (&_S1615)->primal_0 = min_Jxz_0;
    (&_S1615)->differential_0 = 0.0f;
    DiffPair_float_0 _S1616;
    (&_S1616)->primal_0 = max_Jxz_0;
    (&_S1616)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1614, &_S1615, &_S1616, _S1610.rows[int(0)].z);
    float _S1617 = - ((_S1530 * _S1612.differential_0 + _S1529 * _S1613.differential_0 + _S1528 * _S1615.differential_0 + _S1527 * _S1616.differential_0) / _S1441);
    float2  _S1618 = make_float2 (0.0f, _S1611.differential_0) + make_float2 (_S1614.differential_0, 0.0f);
    float _S1619 = fx_12 * _S1618.x;
    float2  _S1620 = make_float2 (_S1619, fy_12 * _S1618.y) + make_float2 (dist_coeffs_12[int(8)] * _S1619, dist_coeffs_12[int(9)] * _S1619);
    float2  _S1621 = _S1498 * _S1620;
    float2  _S1622 = _S1500 * _S1620;
    float _S1623 = dist_coeffs_12[int(4)] * _S1620.y;
    float _S1624 = dist_coeffs_12[int(5)] * _S1620.x;
    float _S1625 = _S1622.x + _S1622.y;
    float _S1626 = _S1621.x + _S1621.y;
    float _S1627 = r2_74 * _S1626;
    float _S1628 = s_diff_r2_20 * _S1626 + r2_74 * _S1625;
    float _S1629 = r2_74 * _S1627;
    float _S1630 = s_diff_r2_20 * _S1627 + r2_74 * _S1628;
    float _S1631 = dist_coeffs_12[int(7)] * _S1620.y + _S1623 + dist_coeffs_12[int(6)] * _S1620.x + _S1624 + _S1508 * _S1626 + _S1506 * _S1627 + _S1504 * _S1629 + dist_coeffs_12[int(3)] * (r2_74 * _S1629);
    float _S1632 = _S1507 * _S1626 + _S1508 * _S1625 + _S1505 * _S1627 + _S1506 * _S1628 + _S1503 * _S1629 + _S1504 * _S1630 + dist_coeffs_12[int(3)] * (s_diff_r2_20 * _S1629 + r2_74 * _S1630);
    float _S1633 = _S1631 + _S1631;
    float _S1634 = v_74 * _S1632;
    float _S1635 = u_74 * _S1632;
    float2  _S1636 = (_S1511 * _S1620 + make_float2 (_S1462 * (v_74 * _S1620.y) + _S1514 * _S1624 + 2.0f * (u_74 * _S1624) + _S1457 * (v_74 * _S1620.x) + u_74 * _S1633, _S1518 * _S1623 + 2.0f * (v_74 * _S1623) + _S1516 * _S1620.y + _S1512 * _S1620.x + v_74 * _S1633)) / _S1445;
    float2  _S1637 = _S1499 * - _S1636;
    float _S1638 = _S1438 * (_S1637.x + _S1637.y);
    float2  _S1639 = (_S1510 * _S1620 + make_float2 (_S1462 * (s_diff_v_20 * _S1620.y) + _S1515 * _S1624 + 2.0f * (s_diff_u_20 * _S1624) + _S1457 * (s_diff_v_20 * _S1620.x) + s_diff_u_20 * _S1633 + _S1635 + _S1635, _S1519 * _S1623 + 2.0f * (s_diff_v_20 * _S1623) + _S1517 * _S1620.y + _S1513 * _S1620.x + s_diff_v_20 * _S1633 + _S1634 + _S1634)) / _S1442;
    float2  _S1640 = _S1436 * - _S1639;
    float2  _S1641 = - (_S1442 * _S1636) + _S1439 * _S1639;
    float3  _S1642 = make_float3 (_S1641.x, _S1641.y, _S1638 + _S1638 + _S1640.x + _S1640.y);
    float2  _S1643 = make_float2 (0.0f, _S1610.rows[int(1)].y) + make_float2 (_S1610.rows[int(0)].y, 0.0f);
    float _S1644 = fx_12 * _S1643.x;
    float2  _S1645 = make_float2 (_S1644, fy_12 * _S1643.y) + make_float2 (dist_coeffs_12[int(8)] * _S1644, dist_coeffs_12[int(9)] * _S1644);
    float2  _S1646 = _S1473 * _S1645;
    float2  _S1647 = _S1475 * _S1645;
    float _S1648 = dist_coeffs_12[int(4)] * _S1645.y;
    float _S1649 = dist_coeffs_12[int(5)] * _S1645.x;
    float _S1650 = _S1647.x + _S1647.y;
    float _S1651 = _S1646.x + _S1646.y;
    float _S1652 = r2_73 * _S1651;
    float _S1653 = s_diff_r2_19 * _S1651 + r2_73 * _S1650;
    float _S1654 = r2_73 * _S1652;
    float _S1655 = s_diff_r2_19 * _S1652 + r2_73 * _S1653;
    float _S1656 = dist_coeffs_12[int(7)] * _S1645.y + _S1648 + dist_coeffs_12[int(6)] * _S1645.x + _S1649 + _S1483 * _S1651 + _S1481 * _S1652 + _S1479 * _S1654 + dist_coeffs_12[int(3)] * (r2_73 * _S1654);
    float _S1657 = _S1482 * _S1651 + _S1483 * _S1650 + _S1480 * _S1652 + _S1481 * _S1653 + _S1478 * _S1654 + _S1479 * _S1655 + dist_coeffs_12[int(3)] * (s_diff_r2_19 * _S1654 + r2_73 * _S1655);
    float _S1658 = _S1656 + _S1656;
    float _S1659 = v_73 * _S1657;
    float _S1660 = u_73 * _S1657;
    float2  _S1661 = (_S1486 * _S1645 + make_float2 (_S1462 * (v_73 * _S1645.y) + _S1489 * _S1649 + 2.0f * (u_73 * _S1649) + _S1457 * (v_73 * _S1645.x) + u_73 * _S1658, _S1493 * _S1648 + 2.0f * (v_73 * _S1648) + _S1491 * _S1645.y + _S1487 * _S1645.x + v_73 * _S1658)) / _S1445;
    float2  _S1662 = _S1474 * - _S1661;
    float _S1663 = _S1438 * (_S1662.x + _S1662.y);
    float2  _S1664 = _S1472 * (_S1442 * _S1661);
    float2  _S1665 = (_S1485 * _S1645 + make_float2 (_S1462 * (s_diff_v_19 * _S1645.y) + _S1490 * _S1649 + 2.0f * (s_diff_u_19 * _S1649) + _S1457 * (s_diff_v_19 * _S1645.x) + s_diff_u_19 * _S1658 + _S1660 + _S1660, _S1494 * _S1648 + 2.0f * (s_diff_v_19 * _S1648) + _S1492 * _S1645.y + _S1488 * _S1645.x + s_diff_v_19 * _S1658 + _S1659 + _S1659)) / _S1442;
    float2  _S1666 = _S1436 * - _S1665;
    float2  _S1667 = _S1439 * _S1665;
    float3  _S1668 = make_float3 (_S1667.x, _S1667.y, _S1663 + _S1663 + _S1664.x + _S1664.y + _S1666.x + _S1666.y);
    float2  _S1669 = make_float2 (0.0f, _S1610.rows[int(1)].x) + make_float2 (_S1610.rows[int(0)].x, 0.0f);
    float _S1670 = fx_12 * _S1669.x;
    float2  _S1671 = make_float2 (_S1670, fy_12 * _S1669.y) + make_float2 (dist_coeffs_12[int(8)] * _S1670, dist_coeffs_12[int(9)] * _S1670);
    float2  _S1672 = _S1440 * _S1671;
    float2  _S1673 = _S1444 * _S1671;
    float _S1674 = dist_coeffs_12[int(4)] * _S1671.y;
    float _S1675 = dist_coeffs_12[int(5)] * _S1671.x;
    float _S1676 = _S1673.x + _S1673.y;
    float _S1677 = _S1672.x + _S1672.y;
    float _S1678 = r2_72 * _S1677;
    float _S1679 = s_diff_r2_18 * _S1677 + r2_72 * _S1676;
    float _S1680 = r2_72 * _S1678;
    float _S1681 = s_diff_r2_18 * _S1678 + r2_72 * _S1679;
    float _S1682 = dist_coeffs_12[int(7)] * _S1671.y + _S1674 + dist_coeffs_12[int(6)] * _S1671.x + _S1675 + _S1453 * _S1677 + _S1451 * _S1678 + _S1449 * _S1680 + dist_coeffs_12[int(3)] * (r2_72 * _S1680);
    float _S1683 = _S1452 * _S1677 + _S1453 * _S1676 + _S1450 * _S1678 + _S1451 * _S1679 + _S1448 * _S1680 + _S1449 * _S1681 + dist_coeffs_12[int(3)] * (s_diff_r2_18 * _S1680 + r2_72 * _S1681);
    float _S1684 = _S1682 + _S1682;
    float _S1685 = v_72 * _S1683;
    float _S1686 = u_72 * _S1683;
    float2  _S1687 = (_S1456 * _S1671 + make_float2 (_S1462 * (v_72 * _S1671.y) + _S1460 * _S1675 + 2.0f * (u_72 * _S1675) + _S1457 * (v_72 * _S1671.x) + u_72 * _S1684, _S1465 * _S1674 + 2.0f * (v_72 * _S1674) + _S1463 * _S1671.y + _S1458 * _S1671.x + v_72 * _S1684)) / _S1445;
    float2  _S1688 = _S1443 * - _S1687;
    float _S1689 = _S1438 * (_S1688.x + _S1688.y);
    float2  _S1690 = _S1437 * (_S1442 * _S1687);
    float2  _S1691 = (_S1455 * _S1671 + make_float2 (_S1462 * (s_diff_v_18 * _S1671.y) + _S1461 * _S1675 + 2.0f * (s_diff_u_18 * _S1675) + _S1457 * (s_diff_v_18 * _S1671.x) + s_diff_u_18 * _S1684 + _S1686 + _S1686, _S1466 * _S1674 + 2.0f * (s_diff_v_18 * _S1674) + _S1464 * _S1671.y + _S1459 * _S1671.x + s_diff_v_18 * _S1684 + _S1685 + _S1685)) / _S1442;
    float2  _S1692 = _S1436 * - _S1691;
    float2  _S1693 = _S1439 * _S1691;
    float3  _S1694 = make_float3 (_S1693.x, _S1693.y, _S1689 + _S1689 + _S1690.x + _S1690.y + _S1692.x + _S1692.y);
    float2  _S1695 = _S1436 / make_float2 (_S1438);
    float _S1696 = fx_12 * v_mean2d_0.x;
    float u_75 = _S1695.x;
    float v_75 = _S1695.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S1697 = dist_coeffs_12[int(2)] + r2_75 * dist_coeffs_12[int(3)];
    float _S1698 = dist_coeffs_12[int(1)] + r2_75 * _S1697;
    float _S1699 = dist_coeffs_12[int(0)] + r2_75 * _S1698;
    float2  _S1700 = make_float2 (_S1696, fy_12 * v_mean2d_0.y) + make_float2 (dist_coeffs_12[int(8)] * _S1696, dist_coeffs_12[int(9)] * _S1696);
    float2  _S1701 = _S1695 * _S1700;
    float _S1702 = dist_coeffs_12[int(4)] * _S1700.y;
    float _S1703 = dist_coeffs_12[int(5)] * _S1700.x;
    float _S1704 = _S1701.x + _S1701.y;
    float _S1705 = r2_75 * _S1704;
    float _S1706 = r2_75 * _S1705;
    float _S1707 = dist_coeffs_12[int(7)] * _S1700.y + _S1702 + dist_coeffs_12[int(6)] * _S1700.x + _S1703 + _S1699 * _S1704 + _S1698 * _S1705 + _S1697 * _S1706 + dist_coeffs_12[int(3)] * (r2_75 * _S1706);
    float _S1708 = v_75 * _S1707;
    float _S1709 = u_75 * _S1707;
    float2  _S1710 = (make_float2 (1.0f + r2_75 * _S1699) * _S1700 + make_float2 (_S1462 * (v_75 * _S1700.y) + 2.0f * u_75 * _S1703 + 2.0f * (u_75 * _S1703) + _S1457 * (v_75 * _S1700.x) + _S1709 + _S1709, 2.0f * v_75 * _S1702 + 2.0f * (v_75 * _S1702) + _S1462 * u_75 * _S1700.y + _S1457 * u_75 * _S1700.x + _S1708 + _S1708)) / _S1442;
    float2  _S1711 = _S1436 * - _S1710;
    float2  _S1712 = _S1439 * _S1710;
    float3  _S1713 = make_float3 (_S1712.x, _S1712.y, _S1711.x + _S1711.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1714;
    (&_S1714)->primal_0 = _S1433;
    (&_S1714)->differential_0 = _S1608;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1715;
    (&_S1715)->primal_0 = _S1434;
    (&_S1715)->differential_0 = _S1608;
    s_bwd_prop_mul_2(&_S1714, &_S1715, _S1609.differential_0);
    Matrix<float, 3, 3>  _S1716 = transpose_3(_S1715.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1717;
    (&_S1717)->primal_0 = R_8;
    (&_S1717)->differential_0 = _S1608;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1718;
    (&_S1718)->primal_0 = _S1432;
    (&_S1718)->differential_0 = _S1608;
    s_bwd_prop_mul_2(&_S1717, &_S1718, _S1714.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1719;
    (&_S1719)->primal_0 = _S1430;
    (&_S1719)->differential_0 = _S1608;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1720;
    (&_S1720)->primal_0 = _S1431;
    (&_S1720)->differential_0 = _S1608;
    s_bwd_prop_mul_2(&_S1719, &_S1720, _S1718.differential_0);
    Matrix<float, 3, 3>  _S1721 = _S1719.differential_0 + transpose_3(_S1720.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1722;
    (&_S1722)->primal_0 = _S1429;
    (&_S1722)->differential_0 = _S1608;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1723;
    (&_S1723)->primal_0 = S_0;
    (&_S1723)->differential_0 = _S1608;
    s_bwd_prop_mul_2(&_S1722, &_S1723, _S1721);
    Matrix<float, 3, 3>  _S1724 = transpose_3(_S1722.differential_0);
    float _S1725 = 2.0f * - _S1724.rows[int(2)].z;
    float _S1726 = 2.0f * _S1724.rows[int(2)].y;
    float _S1727 = 2.0f * _S1724.rows[int(2)].x;
    float _S1728 = 2.0f * _S1724.rows[int(1)].z;
    float _S1729 = 2.0f * - _S1724.rows[int(1)].y;
    float _S1730 = 2.0f * _S1724.rows[int(1)].x;
    float _S1731 = 2.0f * _S1724.rows[int(0)].z;
    float _S1732 = 2.0f * _S1724.rows[int(0)].y;
    float _S1733 = 2.0f * - _S1724.rows[int(0)].x;
    float _S1734 = - _S1730 + _S1732;
    float _S1735 = _S1727 + - _S1731;
    float _S1736 = - _S1726 + _S1728;
    float _S1737 = _S1726 + _S1728;
    float _S1738 = _S1727 + _S1731;
    float _S1739 = _S1730 + _S1732;
    float _S1740 = _S1426.w * (_S1729 + _S1733);
    float _S1741 = _S1426.z * (_S1725 + _S1733);
    float _S1742 = _S1426.y * (_S1725 + _S1729);
    float _S1743 = _S1426.x * _S1734 + _S1426.z * _S1737 + _S1426.y * _S1738 + _S1740 + _S1740;
    float _S1744 = _S1426.x * _S1735 + _S1426.w * _S1737 + _S1426.y * _S1739 + _S1741 + _S1741;
    float _S1745 = _S1426.x * _S1736 + _S1426.w * _S1738 + _S1426.z * _S1739 + _S1742 + _S1742;
    float _S1746 = _S1426.w * _S1734 + _S1426.z * _S1735 + _S1426.y * _S1736;
    float3  _S1747 = _S1551;
    *&((&_S1747)->z) = _S1723.differential_0.rows[int(2)].z;
    *&((&_S1747)->y) = _S1723.differential_0.rows[int(1)].y;
    *&((&_S1747)->x) = _S1723.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1748;
    (&_S1748)->primal_0 = scale_8;
    (&_S1748)->differential_0 = _S1551;
    s_bwd_prop_exp_1(&_S1748, _S1747);
    float4  _S1749 = make_float4 (0.0f);
    float4  _S1750 = _S1749;
    *&((&_S1750)->w) = _S1743;
    *&((&_S1750)->z) = _S1744;
    *&((&_S1750)->y) = _S1745;
    *&((&_S1750)->x) = _S1746;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1751;
    (&_S1751)->primal_0 = quat_8;
    (&_S1751)->differential_0 = _S1749;
    s_bwd_normalize_impl_0(&_S1751, _S1750);
    float _S1752 = - (_S1575 / _S1425);
    DiffPair_float_0 _S1753;
    (&_S1753)->primal_0 = _S1422;
    (&_S1753)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1753, _S1752);
    float _S1754 = - _S1753.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1755;
    (&_S1755)->primal_0 = mean_c_8;
    (&_S1755)->differential_0 = _S1551;
    s_bwd_length_impl_0(&_S1755, v_depth_0);
    float3  _S1756 = _S1642 + _S1668 + _S1694 + _S1713 + _S1755.differential_0 + make_float3 (0.0f, 0.0f, _S1617);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1757;
    (&_S1757)->primal_0 = R_8;
    (&_S1757)->differential_0 = _S1608;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1758;
    (&_S1758)->primal_0 = mean_9;
    (&_S1758)->differential_0 = _S1551;
    s_bwd_prop_mul_3(&_S1757, &_S1758, _S1756);
    Matrix<float, 3, 3>  _S1759 = _S1716 + _S1717.differential_0 + _S1757.differential_0;
    float _S1760 = _S1754 + _S1558.differential_0;
    float3  _S1761 = _S1748.differential_0 + _S1557.differential_0;
    *v_mean_0 = *v_mean_0 + (_S1758.differential_0 + _S1556.differential_0);
    *v_quat_0 = *v_quat_0 + _S1751.differential_0;
    *v_scale_0 = *v_scale_0 + _S1761;
    *v_in_opacity_0 = *v_in_opacity_0 + _S1760;
    *v_R_0 = *v_R_0 + _S1759;
    *v_t_0 = *v_t_0 + _S1756;
    return;
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_0)
{
    DiffPair_float_0 _S1762 = { s_primal_ctx_sqrt_0(dpdpx_0->primal_0), 0.5f / s_primal_ctx_sqrt_0((F32_max((1.00000001168609742e-07f), (dpdpx_0->primal_0)))) * dpdpx_0->differential_0 };
    return _S1762;
}

inline __device__ DiffPair_float_0 s_primal_ctx_s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_1)
{
    float _S1763 = *&((&dpdpx_1->differential_0)->x) * *&((&dpdpx_1->primal_0)->x);
    float _S1764 = *&((&dpdpx_1->differential_0)->y) * *&((&dpdpx_1->primal_0)->y);
    float s_diff_len_1 = _S1763 + _S1763 + (_S1764 + _S1764);
    DiffPair_float_0 _S1765;
    (&_S1765)->primal_0 = *&((&dpdpx_1->primal_0)->x) * *&((&dpdpx_1->primal_0)->x) + *&((&dpdpx_1->primal_0)->y) * *&((&dpdpx_1->primal_0)->y);
    (&_S1765)->differential_0 = s_diff_len_1;
    DiffPair_float_0 _S1766 = s_primal_ctx_d_sqrt_0(&_S1765);
    DiffPair_float_0 _S1767 = { _S1766.primal_0, _S1766.differential_0 };
    return _S1767;
}

inline __device__ float s_primal_ctx_atan2_0(float _S1768, float _S1769)
{
    return (F32_atan2((_S1768), (_S1769)));
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_2)
{
    float _S1770 = dpdpx_2->primal_0 * dpdpx_2->primal_0 + dpdpy_0->primal_0 * dpdpy_0->primal_0;
    DiffPair_float_0 _S1771 = { s_primal_ctx_atan2_0(dpdpy_0->primal_0, dpdpx_2->primal_0), - dpdpy_0->primal_0 / _S1770 * dpdpx_2->differential_0 + dpdpx_2->primal_0 / _S1770 * dpdpy_0->differential_0 };
    return _S1771;
}

struct DiffPair_0
{
    DiffPair_float_0 primal_0;
    DiffPair_float_0 differential_0;
};

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S1772, DiffPair_float_0 * _S1773, float _S1774)
{
    _d_atan2_0(_S1772, _S1773, _S1774);
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_0 * dpdpy_1, DiffPair_0 * dpdpx_3, DiffPair_float_0 * _s_dOut_4)
{
    float _S1775 = - (*dpdpy_1).primal_0.primal_0;
    float _S1776 = (*dpdpx_3).primal_0.primal_0;
    float _S1777 = _S1776 * _S1776 + (*dpdpy_1).primal_0.primal_0 * (*dpdpy_1).primal_0.primal_0;
    float _S1778 = _S1775 / _S1777;
    float _S1779 = _S1777 * _S1777;
    float _S1780 = (*dpdpx_3).primal_0.primal_0 / _S1777;
    DiffPair_float_0 _S1781;
    (&_S1781)->primal_0 = (*dpdpy_1).primal_0.primal_0;
    (&_S1781)->differential_0 = 0.0f;
    DiffPair_float_0 _S1782;
    (&_S1782)->primal_0 = (*dpdpx_3).primal_0.primal_0;
    (&_S1782)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1781, &_S1782, _s_dOut_4->primal_0);
    float _S1783 = _S1780 * _s_dOut_4->differential_0;
    float _S1784 = (*dpdpy_1).primal_0.differential_0 * _s_dOut_4->differential_0 / _S1779;
    float _S1785 = (*dpdpx_3).primal_0.primal_0 * - _S1784;
    float _S1786 = (*dpdpy_1).primal_0.primal_0 * _S1785;
    float _S1787 = (*dpdpx_3).primal_0.primal_0 * _S1785;
    float _S1788 = (*dpdpx_3).primal_0.differential_0 * _s_dOut_4->differential_0 / _S1779;
    float _S1789 = _S1775 * - _S1788;
    float _S1790 = (*dpdpy_1).primal_0.primal_0 * _S1789;
    float _S1791 = (*dpdpx_3).primal_0.primal_0 * _S1789;
    float _S1792 = - (_S1777 * _S1788);
    DiffPair_float_0 _S1793 = { _S1782.differential_0 + _S1787 + _S1787 + _S1777 * _S1784 + _S1791 + _S1791, _S1778 * _s_dOut_4->differential_0 };
    dpdpx_3->primal_0 = (*dpdpx_3).primal_0;
    dpdpx_3->differential_0 = _S1793;
    DiffPair_float_0 _S1794 = { _S1781.differential_0 + _S1786 + _S1786 + _S1790 + _S1790 + _S1792, _S1783 };
    dpdpy_1->primal_0 = (*dpdpy_1).primal_0;
    dpdpy_1->differential_0 = _S1794;
    return;
}

struct DiffPair_1
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 differential_0;
};

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_0 * dpdpx_4, DiffPair_float_0 * _s_dOut_5)
{
    float _S1795 = (F32_max((1.00000001168609742e-07f), ((*dpdpx_4).primal_0.primal_0)));
    float _S1796 = s_primal_ctx_sqrt_0(_S1795);
    float _S1797 = 0.5f / _S1796 * _s_dOut_5->differential_0;
    float _S1798 = 0.5f * - ((*dpdpx_4).primal_0.differential_0 * _s_dOut_5->differential_0 / (_S1796 * _S1796));
    DiffPair_float_0 _S1799;
    (&_S1799)->primal_0 = _S1795;
    (&_S1799)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1799, _S1798);
    DiffPair_float_0 _S1800;
    (&_S1800)->primal_0 = 1.00000001168609742e-07f;
    (&_S1800)->differential_0 = 0.0f;
    DiffPair_float_0 _S1801;
    (&_S1801)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S1801)->differential_0 = 0.0f;
    _d_max_0(&_S1800, &_S1801, _S1799.differential_0);
    DiffPair_float_0 _S1802;
    (&_S1802)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S1802)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1802, _s_dOut_5->primal_0);
    DiffPair_float_0 _S1803 = { _S1801.differential_0 + _S1802.differential_0, _S1797 };
    dpdpx_4->primal_0 = (*dpdpx_4).primal_0;
    dpdpx_4->differential_0 = _S1803;
    return;
}

inline __device__ void s_bwd_prop_s_fwd_length_impl_0(DiffPair_1 * dpdpx_5, DiffPair_float_0 * _s_dOut_6)
{
    float _S1804 = (*dpdpx_5).primal_0.primal_0.x;
    float _S1805 = (*dpdpx_5).primal_0.differential_0.x * (*dpdpx_5).primal_0.primal_0.x;
    float _S1806 = (*dpdpx_5).primal_0.primal_0.y;
    float _S1807 = (*dpdpx_5).primal_0.differential_0.y * (*dpdpx_5).primal_0.primal_0.y;
    DiffPair_float_0 _S1808 = { _S1804 * _S1804 + _S1806 * _S1806, _S1805 + _S1805 + (_S1807 + _S1807) };
    DiffPair_float_0 _S1809 = { 0.0f, 0.0f };
    DiffPair_0 _S1810;
    (&_S1810)->primal_0 = _S1808;
    (&_S1810)->differential_0 = _S1809;
    DiffPair_float_0 _S1811;
    (&_S1811)->primal_0 = _s_dOut_6->primal_0;
    (&_S1811)->differential_0 = _s_dOut_6->differential_0;
    s_bwd_prop_d_sqrt_0(&_S1810, &_S1811);
    float _S1812 = _S1810.differential_0.differential_0;
    float _S1813 = _S1812 + _S1812;
    float _S1814 = (*dpdpx_5).primal_0.primal_0.y * _S1813;
    float _S1815 = (*dpdpx_5).primal_0.primal_0.y * _S1810.differential_0.primal_0;
    float _S1816 = (*dpdpx_5).primal_0.differential_0.y * _S1813 + _S1815 + _S1815;
    float _S1817 = (*dpdpx_5).primal_0.primal_0.x * _S1813;
    float _S1818 = (*dpdpx_5).primal_0.primal_0.x * _S1810.differential_0.primal_0;
    float _S1819 = (*dpdpx_5).primal_0.differential_0.x * _S1813 + _S1818 + _S1818;
    float2  _S1820 = make_float2 (0.0f);
    float2  _S1821 = _S1820;
    *&((&_S1821)->y) = _S1816;
    *&((&_S1821)->x) = _S1819;
    float2  _S1822 = _S1820;
    *&((&_S1822)->y) = _S1814;
    *&((&_S1822)->x) = _S1817;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1823 = { _S1821, _S1822 };
    dpdpx_5->primal_0 = (*dpdpx_5).primal_0;
    dpdpx_5->differential_0 = _S1823;
    return;
}

inline __device__ void s_bwd_prop_length_impl_2(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float _s_dOut_7)
{
    float _S1824 = (*dpx_18).primal_0.x;
    float _S1825 = (*dpx_18).primal_0.y;
    DiffPair_float_0 _S1826;
    (&_S1826)->primal_0 = _S1824 * _S1824 + _S1825 * _S1825;
    (&_S1826)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1826, _s_dOut_7);
    float _S1827 = (*dpx_18).primal_0.y * _S1826.differential_0;
    float _S1828 = _S1827 + _S1827;
    float _S1829 = (*dpx_18).primal_0.x * _S1826.differential_0;
    float _S1830 = _S1829 + _S1829;
    float2  _S1831 = make_float2 (0.0f);
    *&((&_S1831)->y) = _S1828;
    *&((&_S1831)->x) = _S1830;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S1831;
    return;
}

inline __device__ void s_bwd_length_impl_2(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1832, float _S1833)
{
    s_bwd_prop_length_impl_2(_S1832, _S1833);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_9, float3  mean_10, float4  quat_9, float3  scale_9, float in_opacity_9, Matrix<float, 3, 3>  R_9, float3  t_9, float fx_13, float fy_13, float cx_10, float cy_10, FixedArray<float, 10>  dist_coeffs_13, uint image_width_9, uint image_height_9, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_9 = s_primal_ctx_mul_0(R_9, mean_10) + t_9;
    float _S1834 = - in_opacity_9;
    float _S1835 = 1.0f + s_primal_ctx_exp_0(_S1834);
    float _S1836 = 1.0f / _S1835;
    float _S1837 = _S1835 * _S1835;
    float4  _S1838 = normalize_0(quat_9);
    float3  _S1839 = s_primal_ctx_exp_1(scale_9);
    float _S1840 = _S1838.y;
    float x2_9 = _S1840 * _S1840;
    float y2_9 = _S1838.z * _S1838.z;
    float z2_9 = _S1838.w * _S1838.w;
    float xy_9 = _S1838.y * _S1838.z;
    float xz_9 = _S1838.y * _S1838.w;
    float yz_9 = _S1838.z * _S1838.w;
    float wx_9 = _S1838.x * _S1838.y;
    float wy_9 = _S1838.x * _S1838.z;
    float wz_9 = _S1838.x * _S1838.w;
    Matrix<float, 3, 3>  _S1841 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_9), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_9), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1839.x, 0.0f, 0.0f, 0.0f, _S1839.y, 0.0f, 0.0f, 0.0f, _S1839.z);
    Matrix<float, 3, 3>  _S1842 = s_primal_ctx_mul_1(_S1841, S_1);
    Matrix<float, 3, 3>  _S1843 = transpose_3(_S1842);
    Matrix<float, 3, 3>  _S1844 = s_primal_ctx_mul_1(_S1842, _S1843);
    Matrix<float, 3, 3>  _S1845 = s_primal_ctx_mul_1(R_9, _S1844);
    Matrix<float, 3, 3>  _S1846 = transpose_3(R_9);
    Matrix<float, 3, 3>  _S1847 = s_primal_ctx_mul_1(_S1845, _S1846);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S1848 = float2 {mean_c_9.x, mean_c_9.y};
    float2  _S1849 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1850 = { _S1848, _S1849 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1851;
    (&_S1851)->primal_0 = _S1848;
    (&_S1851)->differential_0 = _S1849;
    DiffPair_float_0 _S1852 = s_primal_ctx_s_fwd_length_impl_0(&_S1851);
    float _S1853 = mean_c_9.z;
    DiffPair_float_0 _S1854 = { _S1852.primal_0, _S1852.differential_0 };
    DiffPair_float_0 _S1855 = { _S1853, 0.0f };
    DiffPair_float_0 _S1856;
    (&_S1856)->primal_0 = _S1852.primal_0;
    (&_S1856)->differential_0 = _S1852.differential_0;
    DiffPair_float_0 _S1857;
    (&_S1857)->primal_0 = _S1853;
    (&_S1857)->differential_0 = 0.0f;
    DiffPair_float_0 _S1858 = s_primal_ctx_d_atan2_0(&_S1856, &_S1857);
    bool _S1859 = (_S1858.primal_0) < 0.00100000004749745f;
    float k_6;
    float s_diff_k_2;
    float _S1860;
    float _S1861;
    float _S1862;
    float _S1863;
    float _S1864;
    float _S1865;
    float _S1866;
    float _S1867;
    if(_S1859)
    {
        float _S1868 = _S1858.differential_0 * _S1858.primal_0;
        float _S1869 = 1.0f - _S1858.primal_0 * _S1858.primal_0 / 3.0f;
        float _S1870 = 0.0f - (_S1868 + _S1868) * 0.3333333432674408f;
        float _S1871 = _S1853 * _S1853;
        float _S1872 = _S1870 * _S1853;
        float _S1873 = _S1872 / _S1871;
        float _S1874 = _S1871 * _S1871;
        k_6 = _S1869 / _S1853;
        s_diff_k_2 = _S1873;
        _S1860 = _S1874;
        _S1861 = _S1872;
        _S1862 = _S1871;
        _S1863 = _S1869;
        _S1864 = _S1870;
        _S1865 = 0.0f;
        _S1866 = 0.0f;
        _S1867 = 0.0f;
    }
    else
    {
        float _S1875 = _S1852.primal_0 * _S1852.primal_0;
        float _S1876 = _S1858.differential_0 * _S1852.primal_0 - _S1858.primal_0 * _S1852.differential_0;
        float _S1877 = _S1876 / _S1875;
        float _S1878 = _S1875 * _S1875;
        k_6 = _S1858.primal_0 / _S1852.primal_0;
        s_diff_k_2 = _S1877;
        _S1860 = 0.0f;
        _S1861 = 0.0f;
        _S1862 = 0.0f;
        _S1863 = 0.0f;
        _S1864 = 0.0f;
        _S1865 = _S1878;
        _S1866 = _S1876;
        _S1867 = _S1875;
    }
    float2  _S1879 = make_float2 (k_6);
    float2  _S1880 = make_float2 (s_diff_k_2);
    float2  _S1881 = _S1848 * make_float2 (k_6);
    float2  _S1882 = _S1849 * make_float2 (k_6) + make_float2 (s_diff_k_2) * _S1848;
    float u_76 = _S1881.x;
    float s_diff_u_21 = _S1882.x;
    float v_76 = _S1881.y;
    float s_diff_v_21 = _S1882.y;
    float _S1883 = s_diff_u_21 * u_76;
    float _S1884 = s_diff_v_21 * v_76;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float s_diff_r2_21 = _S1883 + _S1883 + (_S1884 + _S1884);
    float _S1885 = s_diff_r2_21 * dist_coeffs_13[int(3)];
    float _S1886 = dist_coeffs_13[int(2)] + r2_76 * dist_coeffs_13[int(3)];
    float _S1887 = s_diff_r2_21 * _S1886 + _S1885 * r2_76;
    float _S1888 = dist_coeffs_13[int(1)] + r2_76 * _S1886;
    float _S1889 = s_diff_r2_21 * _S1888 + _S1887 * r2_76;
    float _S1890 = dist_coeffs_13[int(0)] + r2_76 * _S1888;
    float _S1891 = s_diff_r2_21 * _S1890 + _S1889 * r2_76;
    float2  _S1892 = make_float2 (_S1891);
    float radial_30 = 1.0f + r2_76 * _S1890;
    float2  _S1893 = make_float2 (radial_30);
    float _S1894 = 2.0f * dist_coeffs_13[int(4)];
    float _S1895 = _S1894 * u_76;
    float _S1896 = s_diff_u_21 * _S1894;
    float _S1897 = 2.0f * u_76;
    float _S1898 = s_diff_u_21 * 2.0f;
    float _S1899 = 2.0f * dist_coeffs_13[int(5)];
    float _S1900 = _S1899 * u_76;
    float _S1901 = s_diff_u_21 * _S1899;
    float _S1902 = 2.0f * v_76;
    float _S1903 = s_diff_v_21 * 2.0f;
    float2  _S1904 = _S1882 * make_float2 (radial_30) + make_float2 (_S1891) * _S1881 + make_float2 (_S1896 * v_76 + s_diff_v_21 * _S1895 + (s_diff_r2_21 + (_S1898 * u_76 + s_diff_u_21 * _S1897)) * dist_coeffs_13[int(5)] + s_diff_r2_21 * dist_coeffs_13[int(6)], _S1901 * v_76 + s_diff_v_21 * _S1900 + (s_diff_r2_21 + (_S1903 * v_76 + s_diff_v_21 * _S1902)) * dist_coeffs_13[int(4)] + s_diff_r2_21 * dist_coeffs_13[int(7)]);
    float2  _S1905 = _S1904 + make_float2 (_S1904.x * dist_coeffs_13[int(8)] + _S1904.y * dist_coeffs_13[int(9)], 0.0f);
    float _S1906 = _S1905.x * fx_13;
    float _S1907 = _S1905.y * fy_13;
    Matrix<float, 2, 3>  _S1908 = J_11;
    *&(((&_S1908)->rows + (int(0)))->x) = _S1906;
    *&(((&_S1908)->rows + (int(1)))->x) = _S1907;
    float2  _S1909 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1910 = { _S1848, _S1909 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1911;
    (&_S1911)->primal_0 = _S1848;
    (&_S1911)->differential_0 = _S1909;
    DiffPair_float_0 _S1912 = s_primal_ctx_s_fwd_length_impl_0(&_S1911);
    DiffPair_float_0 _S1913 = { _S1912.primal_0, _S1912.differential_0 };
    DiffPair_float_0 _S1914;
    (&_S1914)->primal_0 = _S1912.primal_0;
    (&_S1914)->differential_0 = _S1912.differential_0;
    DiffPair_float_0 _S1915;
    (&_S1915)->primal_0 = _S1853;
    (&_S1915)->differential_0 = 0.0f;
    DiffPair_float_0 _S1916 = s_primal_ctx_d_atan2_0(&_S1914, &_S1915);
    bool _S1917 = (_S1916.primal_0) < 0.00100000004749745f;
    float _S1918;
    float _S1919;
    float _S1920;
    float _S1921;
    float _S1922;
    float _S1923;
    float _S1924;
    float _S1925;
    if(_S1917)
    {
        float _S1926 = _S1916.differential_0 * _S1916.primal_0;
        float _S1927 = 1.0f - _S1916.primal_0 * _S1916.primal_0 / 3.0f;
        float _S1928 = 0.0f - (_S1926 + _S1926) * 0.3333333432674408f;
        float _S1929 = _S1853 * _S1853;
        float _S1930 = _S1928 * _S1853;
        float _S1931 = _S1930 / _S1929;
        float _S1932 = _S1929 * _S1929;
        k_6 = _S1927 / _S1853;
        s_diff_k_2 = _S1931;
        _S1918 = _S1932;
        _S1919 = _S1930;
        _S1920 = _S1929;
        _S1921 = _S1927;
        _S1922 = _S1928;
        _S1923 = 0.0f;
        _S1924 = 0.0f;
        _S1925 = 0.0f;
    }
    else
    {
        float _S1933 = _S1912.primal_0 * _S1912.primal_0;
        float _S1934 = _S1916.differential_0 * _S1912.primal_0 - _S1916.primal_0 * _S1912.differential_0;
        float _S1935 = _S1934 / _S1933;
        float _S1936 = _S1933 * _S1933;
        k_6 = _S1916.primal_0 / _S1912.primal_0;
        s_diff_k_2 = _S1935;
        _S1918 = 0.0f;
        _S1919 = 0.0f;
        _S1920 = 0.0f;
        _S1921 = 0.0f;
        _S1922 = 0.0f;
        _S1923 = _S1936;
        _S1924 = _S1934;
        _S1925 = _S1933;
    }
    float2  _S1937 = make_float2 (k_6);
    float2  _S1938 = make_float2 (s_diff_k_2);
    float2  _S1939 = _S1848 * make_float2 (k_6);
    float2  _S1940 = _S1909 * make_float2 (k_6) + make_float2 (s_diff_k_2) * _S1848;
    float u_77 = _S1939.x;
    float s_diff_u_22 = _S1940.x;
    float v_77 = _S1939.y;
    float s_diff_v_22 = _S1940.y;
    float _S1941 = s_diff_u_22 * u_77;
    float _S1942 = s_diff_v_22 * v_77;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float s_diff_r2_22 = _S1941 + _S1941 + (_S1942 + _S1942);
    float _S1943 = s_diff_r2_22 * dist_coeffs_13[int(3)];
    float _S1944 = dist_coeffs_13[int(2)] + r2_77 * dist_coeffs_13[int(3)];
    float _S1945 = s_diff_r2_22 * _S1944 + _S1943 * r2_77;
    float _S1946 = dist_coeffs_13[int(1)] + r2_77 * _S1944;
    float _S1947 = s_diff_r2_22 * _S1946 + _S1945 * r2_77;
    float _S1948 = dist_coeffs_13[int(0)] + r2_77 * _S1946;
    float _S1949 = s_diff_r2_22 * _S1948 + _S1947 * r2_77;
    float2  _S1950 = make_float2 (_S1949);
    float radial_31 = 1.0f + r2_77 * _S1948;
    float2  _S1951 = make_float2 (radial_31);
    float _S1952 = _S1894 * u_77;
    float _S1953 = s_diff_u_22 * _S1894;
    float _S1954 = 2.0f * u_77;
    float _S1955 = s_diff_u_22 * 2.0f;
    float _S1956 = _S1899 * u_77;
    float _S1957 = s_diff_u_22 * _S1899;
    float _S1958 = 2.0f * v_77;
    float _S1959 = s_diff_v_22 * 2.0f;
    float2  _S1960 = _S1940 * make_float2 (radial_31) + make_float2 (_S1949) * _S1939 + make_float2 (_S1953 * v_77 + s_diff_v_22 * _S1952 + (s_diff_r2_22 + (_S1955 * u_77 + s_diff_u_22 * _S1954)) * dist_coeffs_13[int(5)] + s_diff_r2_22 * dist_coeffs_13[int(6)], _S1957 * v_77 + s_diff_v_22 * _S1956 + (s_diff_r2_22 + (_S1959 * v_77 + s_diff_v_22 * _S1958)) * dist_coeffs_13[int(4)] + s_diff_r2_22 * dist_coeffs_13[int(7)]);
    float2  _S1961 = _S1960 + make_float2 (_S1960.x * dist_coeffs_13[int(8)] + _S1960.y * dist_coeffs_13[int(9)], 0.0f);
    float _S1962 = _S1961.y * fy_13;
    *&(((&_S1908)->rows + (int(0)))->y) = _S1961.x * fx_13;
    *&(((&_S1908)->rows + (int(1)))->y) = _S1962;
    float2  _S1963 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1964 = { _S1848, _S1963 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1965;
    (&_S1965)->primal_0 = _S1848;
    (&_S1965)->differential_0 = _S1963;
    DiffPair_float_0 _S1966 = s_primal_ctx_s_fwd_length_impl_0(&_S1965);
    DiffPair_float_0 _S1967 = { _S1966.primal_0, _S1966.differential_0 };
    DiffPair_float_0 _S1968 = { _S1853, 1.0f };
    DiffPair_float_0 _S1969;
    (&_S1969)->primal_0 = _S1966.primal_0;
    (&_S1969)->differential_0 = _S1966.differential_0;
    DiffPair_float_0 _S1970;
    (&_S1970)->primal_0 = _S1853;
    (&_S1970)->differential_0 = 1.0f;
    DiffPair_float_0 _S1971 = s_primal_ctx_d_atan2_0(&_S1969, &_S1970);
    bool _S1972 = (_S1971.primal_0) < 0.00100000004749745f;
    float _S1973;
    float _S1974;
    float _S1975;
    float _S1976;
    float _S1977;
    float _S1978;
    float _S1979;
    float _S1980;
    if(_S1972)
    {
        float _S1981 = _S1971.differential_0 * _S1971.primal_0;
        float _S1982 = 1.0f - _S1971.primal_0 * _S1971.primal_0 / 3.0f;
        float _S1983 = 0.0f - (_S1981 + _S1981) * 0.3333333432674408f;
        float _S1984 = _S1853 * _S1853;
        float _S1985 = _S1983 * _S1853 - _S1982;
        float _S1986 = _S1985 / _S1984;
        float _S1987 = _S1984 * _S1984;
        k_6 = _S1982 / _S1853;
        s_diff_k_2 = _S1986;
        _S1973 = _S1987;
        _S1974 = _S1985;
        _S1975 = _S1984;
        _S1976 = _S1982;
        _S1977 = _S1983;
        _S1978 = 0.0f;
        _S1979 = 0.0f;
        _S1980 = 0.0f;
    }
    else
    {
        float _S1988 = _S1966.primal_0 * _S1966.primal_0;
        float _S1989 = _S1971.differential_0 * _S1966.primal_0 - _S1971.primal_0 * _S1966.differential_0;
        float _S1990 = _S1989 / _S1988;
        float _S1991 = _S1988 * _S1988;
        k_6 = _S1971.primal_0 / _S1966.primal_0;
        s_diff_k_2 = _S1990;
        _S1973 = 0.0f;
        _S1974 = 0.0f;
        _S1975 = 0.0f;
        _S1976 = 0.0f;
        _S1977 = 0.0f;
        _S1978 = _S1991;
        _S1979 = _S1989;
        _S1980 = _S1988;
    }
    float2  _S1992 = make_float2 (k_6);
    float2  _S1993 = make_float2 (s_diff_k_2);
    float2  _S1994 = _S1848 * make_float2 (k_6);
    float2  _S1995 = make_float2 (s_diff_k_2) * _S1848;
    float u_78 = _S1994.x;
    float s_diff_u_23 = _S1995.x;
    float v_78 = _S1994.y;
    float s_diff_v_23 = _S1995.y;
    float _S1996 = s_diff_u_23 * u_78;
    float _S1997 = s_diff_v_23 * v_78;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float s_diff_r2_23 = _S1996 + _S1996 + (_S1997 + _S1997);
    float _S1998 = s_diff_r2_23 * dist_coeffs_13[int(3)];
    float _S1999 = dist_coeffs_13[int(2)] + r2_78 * dist_coeffs_13[int(3)];
    float _S2000 = s_diff_r2_23 * _S1999 + _S1998 * r2_78;
    float _S2001 = dist_coeffs_13[int(1)] + r2_78 * _S1999;
    float _S2002 = s_diff_r2_23 * _S2001 + _S2000 * r2_78;
    float _S2003 = dist_coeffs_13[int(0)] + r2_78 * _S2001;
    float _S2004 = s_diff_r2_23 * _S2003 + _S2002 * r2_78;
    float2  _S2005 = make_float2 (_S2004);
    float radial_32 = 1.0f + r2_78 * _S2003;
    float2  _S2006 = make_float2 (radial_32);
    float _S2007 = _S1894 * u_78;
    float _S2008 = s_diff_u_23 * _S1894;
    float _S2009 = 2.0f * u_78;
    float _S2010 = s_diff_u_23 * 2.0f;
    float _S2011 = _S1899 * u_78;
    float _S2012 = s_diff_u_23 * _S1899;
    float _S2013 = 2.0f * v_78;
    float _S2014 = s_diff_v_23 * 2.0f;
    float2  _S2015 = _S1995 * make_float2 (radial_32) + make_float2 (_S2004) * _S1994 + make_float2 (_S2008 * v_78 + s_diff_v_23 * _S2007 + (s_diff_r2_23 + (_S2010 * u_78 + s_diff_u_23 * _S2009)) * dist_coeffs_13[int(5)] + s_diff_r2_23 * dist_coeffs_13[int(6)], _S2012 * v_78 + s_diff_v_23 * _S2011 + (s_diff_r2_23 + (_S2014 * v_78 + s_diff_v_23 * _S2013)) * dist_coeffs_13[int(4)] + s_diff_r2_23 * dist_coeffs_13[int(7)]);
    float2  _S2016 = _S2015 + make_float2 (_S2015.x * dist_coeffs_13[int(8)] + _S2015.y * dist_coeffs_13[int(9)], 0.0f);
    float _S2017 = _S2016.y * fy_13;
    *&(((&_S1908)->rows + (int(0)))->z) = _S2016.x * fx_13;
    *&(((&_S1908)->rows + (int(1)))->z) = _S2017;
    Matrix<float, 2, 3>  _S2018 = s_primal_ctx_mul_2(_S1908, _S1847);
    Matrix<float, 3, 2>  _S2019 = transpose_1(_S1908);
    Matrix<float, 2, 2>  _S2020 = s_primal_ctx_mul_3(_S2018, _S2019);
    float eps2d_9;
    if(antialiased_9)
    {
        eps2d_9 = 0.10000000149011612f;
    }
    else
    {
        eps2d_9 = 0.30000001192092896f;
    }
    float _S2021 = _S2020.rows[int(0)].y * _S2020.rows[int(1)].x;
    float det_orig_9 = _S2020.rows[int(0)].x * _S2020.rows[int(1)].y - _S2021;
    float _S2022 = _S2020.rows[int(0)].x + eps2d_9;
    Matrix<float, 2, 2>  _S2023 = _S2020;
    *&(((&_S2023)->rows + (int(0)))->x) = _S2022;
    float _S2024 = _S2020.rows[int(1)].y + eps2d_9;
    *&(((&_S2023)->rows + (int(1)))->y) = _S2024;
    Matrix<float, 2, 2>  _S2025 = _S2023;
    Matrix<float, 2, 2>  _S2026 = _S2023;
    float det_blur_9 = _S2022 * _S2024 - _S2021;
    float _S2027 = det_orig_9 / det_blur_9;
    float _S2028 = det_blur_9 * det_blur_9;
    float _S2029 = (F32_max((0.0f), (_S2027)));
    float _S2030 = s_primal_ctx_sqrt_0(_S2029);
    float invdet_11 = 1.0f / det_blur_9;
    float _S2031 = - _S2020.rows[int(0)].y;
    float _S2032 = - _S2020.rows[int(1)].x;
    if(antialiased_9)
    {
        k_6 = _S1836 * _S2030;
    }
    else
    {
        k_6 = _S1836;
    }
    float _S2033 = k_6 / 0.00392156885936856f;
    float _S2034 = 2.0f * s_primal_ctx_log_0(_S2033);
    float _S2035 = s_primal_ctx_sqrt_0(_S2034);
    float _S2036 = _S2025.rows[int(0)].x;
    float _S2037 = _S2026.rows[int(1)].y;
    float3  campos_2 = - s_primal_ctx_mul_0(_S1846, t_9);
    float3  _S2038 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2039;
    (&_S2039)->primal_0 = mean_10;
    (&_S2039)->differential_0 = _S2038;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2040;
    (&_S2040)->primal_0 = scale_9;
    (&_S2040)->differential_0 = _S2038;
    DiffPair_float_0 _S2041;
    (&_S2041)->primal_0 = in_opacity_9;
    (&_S2041)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2042;
    (&_S2042)->primal_0 = campos_2;
    (&_S2042)->differential_0 = _S2038;
    s_bwd_prop_view_radius_3dgs_0(&_S2039, &_S2040, &_S2041, &_S2042, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2043 = _S2039;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2044 = _S2040;
    DiffPair_float_0 _S2045 = _S2041;
    float2  _S2046 = make_float2 (0.0f);
    float2  _S2047 = _S2046;
    *&((&_S2047)->y) = v_conic_1.z;
    float2  _S2048 = _S2046;
    *&((&_S2048)->y) = v_conic_1.y;
    *&((&_S2048)->x) = v_conic_1.x;
    DiffPair_float_0 _S2049;
    (&_S2049)->primal_0 = _S2037;
    (&_S2049)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2049, 0.0f);
    DiffPair_float_0 _S2050;
    (&_S2050)->primal_0 = _S2036;
    (&_S2050)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2050, 0.0f);
    DiffPair_float_0 _S2051;
    (&_S2051)->primal_0 = 3.32999992370605469f;
    (&_S2051)->differential_0 = 0.0f;
    DiffPair_float_0 _S2052;
    (&_S2052)->primal_0 = _S2035;
    (&_S2052)->differential_0 = 0.0f;
    _d_min_0(&_S2051, &_S2052, 0.0f);
    DiffPair_float_0 _S2053;
    (&_S2053)->primal_0 = _S2034;
    (&_S2053)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2053, _S2052.differential_0);
    float _S2054 = 2.0f * _S2053.differential_0;
    DiffPair_float_0 _S2055;
    (&_S2055)->primal_0 = _S2033;
    (&_S2055)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2055, _S2054);
    float _S2056 = v_opacity_1 + 254.9999847412109375f * _S2055.differential_0;
    float2  _S2057 = make_float2 (_S2050.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S2058 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2059 = _S2058;
    _S2059[int(1)] = _S2047;
    _S2059[int(0)] = _S2048;
    Matrix<float, 2, 2>  _S2060 = _S2059;
    float2  _S2061 = make_float2 (0.0f, _S2049.differential_0);
    if(antialiased_9)
    {
        float _S2062 = _S2030 * _S2056;
        k_6 = _S1836 * _S2056;
        s_diff_k_2 = _S2062;
    }
    else
    {
        k_6 = 0.0f;
        s_diff_k_2 = _S2056;
    }
    float _S2063 = invdet_11 * _S2060.rows[int(1)].y;
    float _S2064 = - (invdet_11 * _S2060.rows[int(1)].x);
    float _S2065 = - (invdet_11 * _S2060.rows[int(0)].y);
    float _S2066 = invdet_11 * _S2060.rows[int(0)].x;
    float _S2067 = - ((_S2022 * _S2060.rows[int(1)].y + _S2032 * _S2060.rows[int(1)].x + _S2031 * _S2060.rows[int(0)].y + _S2024 * _S2060.rows[int(0)].x) / _S2028);
    DiffPair_float_0 _S2068;
    (&_S2068)->primal_0 = _S2029;
    (&_S2068)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2068, k_6);
    DiffPair_float_0 _S2069 = { 0.0f, 0.0f };
    DiffPair_float_0 _S2070;
    (&_S2070)->primal_0 = 0.0f;
    (&_S2070)->differential_0 = 0.0f;
    DiffPair_float_0 _S2071;
    (&_S2071)->primal_0 = _S2027;
    (&_S2071)->differential_0 = 0.0f;
    _d_max_0(&_S2070, &_S2071, _S2068.differential_0);
    float _S2072 = _S2071.differential_0 / _S2028;
    float s_diff_det_orig_T_1 = det_blur_9 * _S2072;
    float _S2073 = det_orig_9 * - _S2072 + _S2067;
    float _S2074 = - _S2073;
    float _S2075 = _S2022 * _S2073;
    float _S2076 = _S2024 * _S2073;
    Matrix<float, 2, 2>  _S2077 = _S2058;
    _S2077[int(1)] = _S2061;
    _S2077[int(0)] = _S2057;
    _S2023 = _S2077;
    *&(((&_S2023)->rows + (int(1)))->y) = 0.0f;
    float _S2078 = _S2075 + _S2077.rows[int(1)].y + _S2066;
    *&(((&_S2023)->rows + (int(0)))->x) = 0.0f;
    float _S2079 = _S2076 + _S2077.rows[int(0)].x + _S2063;
    float _S2080 = _S2074 + - s_diff_det_orig_T_1;
    float _S2081 = _S2020.rows[int(0)].y * _S2080 + _S2064;
    float _S2082 = _S2020.rows[int(1)].x * _S2080 + _S2065;
    float _S2083 = _S2020.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2084 = _S2078 + _S2020.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2085 = _S2046;
    *&((&_S2085)->x) = _S2081;
    *&((&_S2085)->y) = _S2084;
    float _S2086 = _S2079 + _S2083;
    float2  _S2087 = _S2046;
    *&((&_S2087)->y) = _S2082;
    *&((&_S2087)->x) = _S2086;
    Matrix<float, 2, 2>  _S2088 = _S2058;
    _S2088[int(1)] = _S2085;
    _S2088[int(0)] = _S2087;
    Matrix<float, 2, 2>  _S2089 = _S2023 + _S2088;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2090;
    (&_S2090)->primal_0 = _S2018;
    (&_S2090)->differential_0 = J_11;
    Matrix<float, 3, 2>  _S2091 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2092;
    (&_S2092)->primal_0 = _S2019;
    (&_S2092)->differential_0 = _S2091;
    s_bwd_prop_mul_0(&_S2090, &_S2092, _S2089);
    Matrix<float, 2, 3>  _S2093 = transpose_2(_S2092.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2094;
    (&_S2094)->primal_0 = _S1908;
    (&_S2094)->differential_0 = J_11;
    Matrix<float, 3, 3>  _S2095 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2096;
    (&_S2096)->primal_0 = _S1847;
    (&_S2096)->differential_0 = _S2095;
    s_bwd_prop_mul_1(&_S2094, &_S2096, _S2090.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2097 = _S2096;
    Matrix<float, 2, 3>  _S2098 = _S2093 + _S2094.differential_0;
    float2  _S2099 = make_float2 (0.0f, _S2098.rows[int(1)].z) + make_float2 (_S2098.rows[int(0)].z, 0.0f);
    float _S2100 = fx_13 * _S2099.x;
    float2  _S2101 = make_float2 (_S2100, fy_13 * _S2099.y) + make_float2 (dist_coeffs_13[int(8)] * _S2100, dist_coeffs_13[int(9)] * _S2100);
    float2  _S2102 = _S1994 * _S2101;
    float2  _S2103 = _S1995 * _S2101;
    float _S2104 = dist_coeffs_13[int(4)] * _S2101.y;
    float _S2105 = dist_coeffs_13[int(5)] * _S2101.x;
    float _S2106 = _S2103.x + _S2103.y;
    float _S2107 = _S2102.x + _S2102.y;
    float _S2108 = r2_78 * _S2107;
    float _S2109 = s_diff_r2_23 * _S2107 + r2_78 * _S2106;
    float _S2110 = r2_78 * _S2108;
    float _S2111 = s_diff_r2_23 * _S2108 + r2_78 * _S2109;
    float _S2112 = dist_coeffs_13[int(7)] * _S2101.y + _S2104 + dist_coeffs_13[int(6)] * _S2101.x + _S2105 + _S2003 * _S2107 + _S2001 * _S2108 + _S1999 * _S2110 + dist_coeffs_13[int(3)] * (r2_78 * _S2110);
    float _S2113 = _S2002 * _S2107 + _S2003 * _S2106 + _S2000 * _S2108 + _S2001 * _S2109 + _S1998 * _S2110 + _S1999 * _S2111 + dist_coeffs_13[int(3)] * (s_diff_r2_23 * _S2110 + r2_78 * _S2111);
    float _S2114 = _S2112 + _S2112;
    float _S2115 = v_78 * _S2113;
    float _S2116 = u_78 * _S2113;
    float2  _S2117 = _S2005 * _S2101 + make_float2 (_S1899 * (s_diff_v_23 * _S2101.y) + _S2010 * _S2105 + 2.0f * (s_diff_u_23 * _S2105) + _S1894 * (s_diff_v_23 * _S2101.x) + s_diff_u_23 * _S2114 + _S2116 + _S2116, _S2014 * _S2104 + 2.0f * (s_diff_v_23 * _S2104) + _S2012 * _S2101.y + _S2008 * _S2101.x + s_diff_v_23 * _S2114 + _S2115 + _S2115);
    float2  _S2118 = _S2006 * _S2101 + make_float2 (_S1899 * (v_78 * _S2101.y) + _S2009 * _S2105 + 2.0f * (u_78 * _S2105) + _S1894 * (v_78 * _S2101.x) + u_78 * _S2114, _S2013 * _S2104 + 2.0f * (v_78 * _S2104) + _S2011 * _S2101.y + _S2007 * _S2101.x + v_78 * _S2114);
    float2  _S2119 = _S1848 * _S2118;
    float2  _S2120 = _S1848 * _S2117;
    float _S2121 = _S2120.x + _S2120.y;
    float _S2122 = _S2119.x + _S2119.y;
    float2  _S2123 = _S1993 * _S2118 + _S1992 * _S2117;
    if(_S1972)
    {
        float _S2124 = _S2122 / _S1973;
        float _S2125 = _S1975 * _S2124;
        float _S2126 = _S1853 * (_S1974 * - _S2124);
        float _S2127 = _S2121 / _S1975;
        float _S2128 = 0.3333333432674408f * - (_S1853 * _S2125);
        float _S2129 = _S2128 + _S2128;
        float _S2130 = _S1971.primal_0 * (0.3333333432674408f * - (- _S2125 + _S1853 * _S2127));
        float _S2131 = _S2126 + _S2126 + _S1977 * _S2125 + _S1976 * - _S2127;
        float _S2132 = _S1971.differential_0 * _S2129 + _S2130 + _S2130;
        k_6 = _S1971.primal_0 * _S2129;
        _S1973 = _S2132;
        _S1974 = _S2131;
        _S1975 = 0.0f;
        _S1976 = 0.0f;
    }
    else
    {
        float _S2133 = _S2122 / _S1978;
        float _S2134 = _S1980 * _S2133;
        float _S2135 = _S1966.primal_0 * (_S1979 * - _S2133);
        float _S2136 = - _S2134;
        float _S2137 = _S1971.primal_0 * _S2136;
        float _S2138 = _S2121 / _S1980;
        float _S2139 = _S2135 + _S2135 + _S1971.differential_0 * _S2134 + _S1971.primal_0 * - _S2138;
        float _S2140 = _S1966.differential_0 * _S2136 + _S1966.primal_0 * _S2138;
        k_6 = _S1966.primal_0 * _S2134;
        _S1973 = _S2140;
        _S1974 = 0.0f;
        _S1975 = _S2137;
        _S1976 = _S2139;
    }
    DiffPair_0 _S2141;
    (&_S2141)->primal_0 = _S1967;
    (&_S2141)->differential_0 = _S2069;
    DiffPair_0 _S2142;
    (&_S2142)->primal_0 = _S1968;
    (&_S2142)->differential_0 = _S2069;
    DiffPair_float_0 _S2143;
    (&_S2143)->primal_0 = _S1973;
    (&_S2143)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2141, &_S2142, &_S2143);
    float _S2144 = _S2142.differential_0.primal_0 + _S1974;
    float _S2145 = _S2141.differential_0.differential_0 + _S1975;
    float _S2146 = _S2141.differential_0.primal_0 + _S1976;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2147 = { _S2046, _S2046 };
    DiffPair_1 _S2148;
    (&_S2148)->primal_0 = _S1964;
    (&_S2148)->differential_0 = _S2147;
    DiffPair_float_0 _S2149;
    (&_S2149)->primal_0 = _S2146;
    (&_S2149)->differential_0 = _S2145;
    s_bwd_prop_s_fwd_length_impl_0(&_S2148, &_S2149);
    float2  _S2150 = _S2148.differential_0.primal_0 + _S2123;
    float3  _S2151 = make_float3 (_S2150.x, _S2150.y, _S2144);
    float2  _S2152 = make_float2 (0.0f, _S2098.rows[int(1)].y) + make_float2 (_S2098.rows[int(0)].y, 0.0f);
    float _S2153 = fx_13 * _S2152.x;
    float2  _S2154 = make_float2 (_S2153, fy_13 * _S2152.y) + make_float2 (dist_coeffs_13[int(8)] * _S2153, dist_coeffs_13[int(9)] * _S2153);
    float2  _S2155 = _S1939 * _S2154;
    float2  _S2156 = _S1940 * _S2154;
    float _S2157 = dist_coeffs_13[int(4)] * _S2154.y;
    float _S2158 = dist_coeffs_13[int(5)] * _S2154.x;
    float _S2159 = _S2156.x + _S2156.y;
    float _S2160 = _S2155.x + _S2155.y;
    float _S2161 = r2_77 * _S2160;
    float _S2162 = s_diff_r2_22 * _S2160 + r2_77 * _S2159;
    float _S2163 = r2_77 * _S2161;
    float _S2164 = s_diff_r2_22 * _S2161 + r2_77 * _S2162;
    float _S2165 = dist_coeffs_13[int(7)] * _S2154.y + _S2157 + dist_coeffs_13[int(6)] * _S2154.x + _S2158 + _S1948 * _S2160 + _S1946 * _S2161 + _S1944 * _S2163 + dist_coeffs_13[int(3)] * (r2_77 * _S2163);
    float _S2166 = _S1947 * _S2160 + _S1948 * _S2159 + _S1945 * _S2161 + _S1946 * _S2162 + _S1943 * _S2163 + _S1944 * _S2164 + dist_coeffs_13[int(3)] * (s_diff_r2_22 * _S2163 + r2_77 * _S2164);
    float _S2167 = _S2165 + _S2165;
    float _S2168 = v_77 * _S2166;
    float _S2169 = u_77 * _S2166;
    float2  _S2170 = _S1950 * _S2154 + make_float2 (_S1899 * (s_diff_v_22 * _S2154.y) + _S1955 * _S2158 + 2.0f * (s_diff_u_22 * _S2158) + _S1894 * (s_diff_v_22 * _S2154.x) + s_diff_u_22 * _S2167 + _S2169 + _S2169, _S1959 * _S2157 + 2.0f * (s_diff_v_22 * _S2157) + _S1957 * _S2154.y + _S1953 * _S2154.x + s_diff_v_22 * _S2167 + _S2168 + _S2168);
    float2  _S2171 = _S1951 * _S2154 + make_float2 (_S1899 * (v_77 * _S2154.y) + _S1954 * _S2158 + 2.0f * (u_77 * _S2158) + _S1894 * (v_77 * _S2154.x) + u_77 * _S2167, _S1958 * _S2157 + 2.0f * (v_77 * _S2157) + _S1956 * _S2154.y + _S1952 * _S2154.x + v_77 * _S2167);
    float2  _S2172 = _S1848 * _S2171;
    float2  _S2173 = _S1909 * _S2171;
    float2  _S2174 = _S1848 * _S2170;
    float _S2175 = _S2173.x + _S2173.y + _S2174.x + _S2174.y;
    float _S2176 = _S2172.x + _S2172.y;
    float2  _S2177 = _S1938 * _S2171 + _S1937 * _S2170;
    if(_S1917)
    {
        float _S2178 = _S2176 / _S1918;
        float _S2179 = _S1920 * _S2178;
        float _S2180 = _S1853 * (_S1919 * - _S2178);
        float _S2181 = _S2175 / _S1920;
        float _S2182 = 0.3333333432674408f * - (_S1853 * _S2179);
        float _S2183 = _S2182 + _S2182;
        float _S2184 = _S1916.primal_0 * (0.3333333432674408f * - (_S1853 * _S2181));
        float _S2185 = _S2180 + _S2180 + _S1922 * _S2179 + _S1921 * - _S2181;
        float _S2186 = _S1916.differential_0 * _S2183 + _S2184 + _S2184;
        k_6 = _S1916.primal_0 * _S2183;
        _S1918 = _S2186;
        _S1919 = _S2185;
        _S1920 = 0.0f;
        _S1921 = 0.0f;
    }
    else
    {
        float _S2187 = _S2176 / _S1923;
        float _S2188 = _S1925 * _S2187;
        float _S2189 = _S1912.primal_0 * (_S1924 * - _S2187);
        float _S2190 = - _S2188;
        float _S2191 = _S1916.primal_0 * _S2190;
        float _S2192 = _S2175 / _S1925;
        float _S2193 = _S2189 + _S2189 + _S1916.differential_0 * _S2188 + _S1916.primal_0 * - _S2192;
        float _S2194 = _S1912.differential_0 * _S2190 + _S1912.primal_0 * _S2192;
        k_6 = _S1912.primal_0 * _S2188;
        _S1918 = _S2194;
        _S1919 = 0.0f;
        _S1920 = _S2191;
        _S1921 = _S2193;
    }
    DiffPair_0 _S2195;
    (&_S2195)->primal_0 = _S1913;
    (&_S2195)->differential_0 = _S2069;
    DiffPair_0 _S2196;
    (&_S2196)->primal_0 = _S1855;
    (&_S2196)->differential_0 = _S2069;
    DiffPair_float_0 _S2197;
    (&_S2197)->primal_0 = _S1918;
    (&_S2197)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2195, &_S2196, &_S2197);
    float _S2198 = _S2196.differential_0.primal_0 + _S1919;
    float _S2199 = _S2195.differential_0.differential_0 + _S1920;
    float _S2200 = _S2195.differential_0.primal_0 + _S1921;
    DiffPair_1 _S2201;
    (&_S2201)->primal_0 = _S1910;
    (&_S2201)->differential_0 = _S2147;
    DiffPair_float_0 _S2202;
    (&_S2202)->primal_0 = _S2200;
    (&_S2202)->differential_0 = _S2199;
    s_bwd_prop_s_fwd_length_impl_0(&_S2201, &_S2202);
    float2  _S2203 = _S2201.differential_0.primal_0 + _S2177;
    float2  _S2204 = make_float2 (0.0f, _S2098.rows[int(1)].x) + make_float2 (_S2098.rows[int(0)].x, 0.0f);
    float _S2205 = fx_13 * _S2204.x;
    float2  _S2206 = make_float2 (_S2205, fy_13 * _S2204.y) + make_float2 (dist_coeffs_13[int(8)] * _S2205, dist_coeffs_13[int(9)] * _S2205);
    float2  _S2207 = _S1881 * _S2206;
    float2  _S2208 = _S1882 * _S2206;
    float _S2209 = dist_coeffs_13[int(4)] * _S2206.y;
    float _S2210 = dist_coeffs_13[int(5)] * _S2206.x;
    float _S2211 = _S2208.x + _S2208.y;
    float _S2212 = _S2207.x + _S2207.y;
    float _S2213 = r2_76 * _S2212;
    float _S2214 = s_diff_r2_21 * _S2212 + r2_76 * _S2211;
    float _S2215 = r2_76 * _S2213;
    float _S2216 = s_diff_r2_21 * _S2213 + r2_76 * _S2214;
    float _S2217 = dist_coeffs_13[int(7)] * _S2206.y + _S2209 + dist_coeffs_13[int(6)] * _S2206.x + _S2210 + _S1890 * _S2212 + _S1888 * _S2213 + _S1886 * _S2215 + dist_coeffs_13[int(3)] * (r2_76 * _S2215);
    float _S2218 = _S1889 * _S2212 + _S1890 * _S2211 + _S1887 * _S2213 + _S1888 * _S2214 + _S1885 * _S2215 + _S1886 * _S2216 + dist_coeffs_13[int(3)] * (s_diff_r2_21 * _S2215 + r2_76 * _S2216);
    float _S2219 = _S2217 + _S2217;
    float _S2220 = v_76 * _S2218;
    float _S2221 = u_76 * _S2218;
    float2  _S2222 = _S1892 * _S2206 + make_float2 (_S1899 * (s_diff_v_21 * _S2206.y) + _S1898 * _S2210 + 2.0f * (s_diff_u_21 * _S2210) + _S1894 * (s_diff_v_21 * _S2206.x) + s_diff_u_21 * _S2219 + _S2221 + _S2221, _S1903 * _S2209 + 2.0f * (s_diff_v_21 * _S2209) + _S1901 * _S2206.y + _S1896 * _S2206.x + s_diff_v_21 * _S2219 + _S2220 + _S2220);
    float2  _S2223 = _S1893 * _S2206 + make_float2 (_S1899 * (v_76 * _S2206.y) + _S1897 * _S2210 + 2.0f * (u_76 * _S2210) + _S1894 * (v_76 * _S2206.x) + u_76 * _S2219, _S1902 * _S2209 + 2.0f * (v_76 * _S2209) + _S1900 * _S2206.y + _S1895 * _S2206.x + v_76 * _S2219);
    float3  _S2224 = make_float3 (_S2203.x, _S2203.y, _S2198) + _S2151;
    float2  _S2225 = _S1848 * _S2223;
    float2  _S2226 = _S1849 * _S2223;
    float2  _S2227 = _S1848 * _S2222;
    float _S2228 = _S2226.x + _S2226.y + _S2227.x + _S2227.y;
    float _S2229 = _S2225.x + _S2225.y;
    float2  _S2230 = _S1880 * _S2223 + _S1879 * _S2222;
    if(_S1859)
    {
        float _S2231 = _S2229 / _S1860;
        float _S2232 = _S1862 * _S2231;
        float _S2233 = _S1853 * (_S1861 * - _S2231);
        float _S2234 = _S2228 / _S1862;
        float _S2235 = 0.3333333432674408f * - (_S1853 * _S2232);
        float _S2236 = _S2235 + _S2235;
        float _S2237 = _S1858.primal_0 * (0.3333333432674408f * - (_S1853 * _S2234));
        float _S2238 = _S2233 + _S2233 + _S1864 * _S2232 + _S1863 * - _S2234;
        float _S2239 = _S1858.differential_0 * _S2236 + _S2237 + _S2237;
        k_6 = _S1858.primal_0 * _S2236;
        _S1860 = _S2239;
        _S1861 = _S2238;
        _S1862 = 0.0f;
        _S1863 = 0.0f;
    }
    else
    {
        float _S2240 = _S2229 / _S1865;
        float _S2241 = _S1867 * _S2240;
        float _S2242 = _S1852.primal_0 * (_S1866 * - _S2240);
        float _S2243 = - _S2241;
        float _S2244 = _S1858.primal_0 * _S2243;
        float _S2245 = _S2228 / _S1867;
        float _S2246 = _S2242 + _S2242 + _S1858.differential_0 * _S2241 + _S1858.primal_0 * - _S2245;
        float _S2247 = _S1852.differential_0 * _S2243 + _S1852.primal_0 * _S2245;
        k_6 = _S1852.primal_0 * _S2241;
        _S1860 = _S2247;
        _S1861 = 0.0f;
        _S1862 = _S2244;
        _S1863 = _S2246;
    }
    DiffPair_0 _S2248;
    (&_S2248)->primal_0 = _S1854;
    (&_S2248)->differential_0 = _S2069;
    DiffPair_0 _S2249;
    (&_S2249)->primal_0 = _S1855;
    (&_S2249)->differential_0 = _S2069;
    DiffPair_float_0 _S2250;
    (&_S2250)->primal_0 = _S1860;
    (&_S2250)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2248, &_S2249, &_S2250);
    float _S2251 = _S2249.differential_0.primal_0 + _S1861;
    float _S2252 = _S2248.differential_0.differential_0 + _S1862;
    float _S2253 = _S2248.differential_0.primal_0 + _S1863;
    DiffPair_1 _S2254;
    (&_S2254)->primal_0 = _S1850;
    (&_S2254)->differential_0 = _S2147;
    DiffPair_float_0 _S2255;
    (&_S2255)->primal_0 = _S2253;
    (&_S2255)->differential_0 = _S2252;
    s_bwd_prop_s_fwd_length_impl_0(&_S2254, &_S2255);
    float2  _S2256 = _S2254.differential_0.primal_0 + _S2230;
    float3  _S2257 = make_float3 (_S2256.x, _S2256.y, _S2251);
    float _S2258 = length_0(_S1848);
    float _S2259 = s_primal_ctx_atan2_0(_S2258, _S1853);
    bool _S2260 = _S2259 < 0.00100000004749745f;
    if(_S2260)
    {
        float _S2261 = 1.0f - _S2259 * _S2259 / 3.0f;
        float _S2262 = _S1853 * _S1853;
        k_6 = _S2261 / _S1853;
        _S1860 = _S2262;
        _S1861 = _S2261;
        _S1862 = 0.0f;
    }
    else
    {
        float _S2263 = _S2258 * _S2258;
        k_6 = _S2259 / _S2258;
        _S1860 = 0.0f;
        _S1861 = 0.0f;
        _S1862 = _S2263;
    }
    float2  _S2264 = make_float2 (k_6);
    float2  _S2265 = _S1848 * make_float2 (k_6);
    float _S2266 = fx_13 * v_mean2d_1.x;
    float u_79 = _S2265.x;
    float v_79 = _S2265.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S2267 = dist_coeffs_13[int(2)] + r2_79 * dist_coeffs_13[int(3)];
    float _S2268 = dist_coeffs_13[int(1)] + r2_79 * _S2267;
    float _S2269 = dist_coeffs_13[int(0)] + r2_79 * _S2268;
    float2  _S2270 = make_float2 (_S2266, fy_13 * v_mean2d_1.y) + make_float2 (dist_coeffs_13[int(8)] * _S2266, dist_coeffs_13[int(9)] * _S2266);
    float2  _S2271 = _S2265 * _S2270;
    float _S2272 = dist_coeffs_13[int(4)] * _S2270.y;
    float _S2273 = dist_coeffs_13[int(5)] * _S2270.x;
    float _S2274 = _S2271.x + _S2271.y;
    float _S2275 = r2_79 * _S2274;
    float _S2276 = r2_79 * _S2275;
    float _S2277 = dist_coeffs_13[int(7)] * _S2270.y + _S2272 + dist_coeffs_13[int(6)] * _S2270.x + _S2273 + _S2269 * _S2274 + _S2268 * _S2275 + _S2267 * _S2276 + dist_coeffs_13[int(3)] * (r2_79 * _S2276);
    float _S2278 = v_79 * _S2277;
    float _S2279 = u_79 * _S2277;
    float2  _S2280 = make_float2 (1.0f + r2_79 * _S2269) * _S2270 + make_float2 (_S1899 * (v_79 * _S2270.y) + 2.0f * u_79 * _S2273 + 2.0f * (u_79 * _S2273) + _S1894 * (v_79 * _S2270.x) + _S2279 + _S2279, 2.0f * v_79 * _S2272 + 2.0f * (v_79 * _S2272) + _S1899 * u_79 * _S2270.y + _S1894 * u_79 * _S2270.x + _S2278 + _S2278);
    float2  _S2281 = _S1848 * _S2280;
    float2  _S2282 = _S2264 * _S2280;
    float _S2283 = _S2281.x + _S2281.y;
    if(_S2260)
    {
        float _S2284 = _S2283 / _S1860;
        float _S2285 = _S1861 * - _S2284;
        float _S2286 = _S2259 * (0.3333333432674408f * - (_S1853 * _S2284));
        k_6 = _S2286 + _S2286;
        _S1860 = _S2285;
        _S1861 = 0.0f;
    }
    else
    {
        float _S2287 = _S2283 / _S1862;
        float _S2288 = _S2259 * - _S2287;
        k_6 = _S2258 * _S2287;
        _S1860 = 0.0f;
        _S1861 = _S2288;
    }
    DiffPair_float_0 _S2289;
    (&_S2289)->primal_0 = _S2258;
    (&_S2289)->differential_0 = 0.0f;
    DiffPair_float_0 _S2290;
    (&_S2290)->primal_0 = _S1853;
    (&_S2290)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2289, &_S2290, k_6);
    float _S2291 = _S2290.differential_0 + _S1860;
    float _S2292 = _S2289.differential_0 + _S1861;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2293;
    (&_S2293)->primal_0 = _S1848;
    (&_S2293)->differential_0 = _S2046;
    s_bwd_length_impl_2(&_S2293, _S2292);
    float2  _S2294 = _S2293.differential_0 + _S2282;
    float3  _S2295 = make_float3 (_S2294.x, _S2294.y, _S2291);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2296;
    (&_S2296)->primal_0 = _S1845;
    (&_S2296)->differential_0 = _S2095;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2297;
    (&_S2297)->primal_0 = _S1846;
    (&_S2297)->differential_0 = _S2095;
    s_bwd_prop_mul_2(&_S2296, &_S2297, _S2097.differential_0);
    Matrix<float, 3, 3>  _S2298 = transpose_3(_S2297.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2299;
    (&_S2299)->primal_0 = R_9;
    (&_S2299)->differential_0 = _S2095;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2300;
    (&_S2300)->primal_0 = _S1844;
    (&_S2300)->differential_0 = _S2095;
    s_bwd_prop_mul_2(&_S2299, &_S2300, _S2296.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2301;
    (&_S2301)->primal_0 = _S1842;
    (&_S2301)->differential_0 = _S2095;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2302;
    (&_S2302)->primal_0 = _S1843;
    (&_S2302)->differential_0 = _S2095;
    s_bwd_prop_mul_2(&_S2301, &_S2302, _S2300.differential_0);
    Matrix<float, 3, 3>  _S2303 = _S2301.differential_0 + transpose_3(_S2302.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2304;
    (&_S2304)->primal_0 = _S1841;
    (&_S2304)->differential_0 = _S2095;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2305;
    (&_S2305)->primal_0 = S_1;
    (&_S2305)->differential_0 = _S2095;
    s_bwd_prop_mul_2(&_S2304, &_S2305, _S2303);
    Matrix<float, 3, 3>  _S2306 = transpose_3(_S2304.differential_0);
    float _S2307 = 2.0f * - _S2306.rows[int(2)].z;
    float _S2308 = 2.0f * _S2306.rows[int(2)].y;
    float _S2309 = 2.0f * _S2306.rows[int(2)].x;
    float _S2310 = 2.0f * _S2306.rows[int(1)].z;
    float _S2311 = 2.0f * - _S2306.rows[int(1)].y;
    float _S2312 = 2.0f * _S2306.rows[int(1)].x;
    float _S2313 = 2.0f * _S2306.rows[int(0)].z;
    float _S2314 = 2.0f * _S2306.rows[int(0)].y;
    float _S2315 = 2.0f * - _S2306.rows[int(0)].x;
    float _S2316 = - _S2312 + _S2314;
    float _S2317 = _S2309 + - _S2313;
    float _S2318 = - _S2308 + _S2310;
    float _S2319 = _S2308 + _S2310;
    float _S2320 = _S2309 + _S2313;
    float _S2321 = _S2312 + _S2314;
    float _S2322 = _S1838.w * (_S2311 + _S2315);
    float _S2323 = _S1838.z * (_S2307 + _S2315);
    float _S2324 = _S1838.y * (_S2307 + _S2311);
    float _S2325 = _S1838.x * _S2316 + _S1838.z * _S2319 + _S1838.y * _S2320 + _S2322 + _S2322;
    float _S2326 = _S1838.x * _S2317 + _S1838.w * _S2319 + _S1838.y * _S2321 + _S2323 + _S2323;
    float _S2327 = _S1838.x * _S2318 + _S1838.w * _S2320 + _S1838.z * _S2321 + _S2324 + _S2324;
    float _S2328 = _S1838.w * _S2316 + _S1838.z * _S2317 + _S1838.y * _S2318;
    float3  _S2329 = _S2038;
    *&((&_S2329)->z) = _S2305.differential_0.rows[int(2)].z;
    *&((&_S2329)->y) = _S2305.differential_0.rows[int(1)].y;
    *&((&_S2329)->x) = _S2305.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2330;
    (&_S2330)->primal_0 = scale_9;
    (&_S2330)->differential_0 = _S2038;
    s_bwd_prop_exp_1(&_S2330, _S2329);
    float4  _S2331 = make_float4 (0.0f);
    float4  _S2332 = _S2331;
    *&((&_S2332)->w) = _S2325;
    *&((&_S2332)->z) = _S2326;
    *&((&_S2332)->y) = _S2327;
    *&((&_S2332)->x) = _S2328;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2333;
    (&_S2333)->primal_0 = quat_9;
    (&_S2333)->differential_0 = _S2331;
    s_bwd_normalize_impl_0(&_S2333, _S2332);
    float _S2334 = - (s_diff_k_2 / _S1837);
    DiffPair_float_0 _S2335;
    (&_S2335)->primal_0 = _S1834;
    (&_S2335)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2335, _S2334);
    float _S2336 = - _S2335.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2337;
    (&_S2337)->primal_0 = mean_c_9;
    (&_S2337)->differential_0 = _S2038;
    s_bwd_length_impl_0(&_S2337, v_depth_1);
    float3  _S2338 = _S2257 + _S2295 + _S2337.differential_0 + _S2224;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2339;
    (&_S2339)->primal_0 = R_9;
    (&_S2339)->differential_0 = _S2095;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2340;
    (&_S2340)->primal_0 = mean_10;
    (&_S2340)->differential_0 = _S2038;
    s_bwd_prop_mul_3(&_S2339, &_S2340, _S2338);
    Matrix<float, 3, 3>  _S2341 = _S2298 + _S2299.differential_0 + _S2339.differential_0;
    float _S2342 = _S2336 + _S2045.differential_0;
    float3  _S2343 = _S2330.differential_0 + _S2044.differential_0;
    *v_mean_1 = *v_mean_1 + (_S2340.differential_0 + _S2043.differential_0);
    *v_quat_1 = *v_quat_1 + _S2333.differential_0;
    *v_scale_1 = *v_scale_1 + _S2343;
    *v_in_opacity_1 = *v_in_opacity_1 + _S2342;
    *v_R_1 = *v_R_1 + _S2341;
    *v_t_1 = *v_t_1 + _S2338;
    return;
}

inline __device__ float s_primal_ctx_sin_0(float _S2344)
{
    return (F32_sin((_S2344)));
}

inline __device__ float s_primal_ctx_cos_0(float _S2345)
{
    return (F32_cos((_S2345)));
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_sin_0(DiffPair_float_0 * dpdpx_6)
{
    DiffPair_float_0 _S2346 = { s_primal_ctx_sin_0(dpdpx_6->primal_0), s_primal_ctx_cos_0(dpdpx_6->primal_0) * dpdpx_6->differential_0 };
    return _S2346;
}

inline __device__ void s_bwd_prop_cos_0(DiffPair_float_0 * _S2347, float _S2348)
{
    _d_cos_0(_S2347, _S2348);
    return;
}

inline __device__ void s_bwd_prop_sin_0(DiffPair_float_0 * _S2349, float _S2350)
{
    _d_sin_0(_S2349, _S2350);
    return;
}

inline __device__ void s_bwd_prop_d_sin_0(DiffPair_0 * dpdpx_7, DiffPair_float_0 * _s_dOut_8)
{
    float _S2351 = s_primal_ctx_cos_0((*dpdpx_7).primal_0.primal_0) * _s_dOut_8->differential_0;
    float _S2352 = (*dpdpx_7).primal_0.differential_0 * _s_dOut_8->differential_0;
    DiffPair_float_0 _S2353;
    (&_S2353)->primal_0 = (*dpdpx_7).primal_0.primal_0;
    (&_S2353)->differential_0 = 0.0f;
    s_bwd_prop_cos_0(&_S2353, _S2352);
    DiffPair_float_0 _S2354;
    (&_S2354)->primal_0 = (*dpdpx_7).primal_0.primal_0;
    (&_S2354)->differential_0 = 0.0f;
    s_bwd_prop_sin_0(&_S2354, _s_dOut_8->primal_0);
    DiffPair_float_0 _S2355 = { _S2353.differential_0 + _S2354.differential_0, _S2351 };
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S2355;
    return;
}

inline __device__ void projection_3dgs_equisolid_vjp(bool antialiased_10, float3  mean_11, float4  quat_10, float3  scale_10, float in_opacity_10, Matrix<float, 3, 3>  R_10, float3  t_10, float fx_14, float fy_14, float cx_11, float cy_11, FixedArray<float, 10>  dist_coeffs_14, uint image_width_10, uint image_height_10, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    DiffPair_float_0 _S2356 = { 0.0f, 0.0f };
    float3  mean_c_10 = s_primal_ctx_mul_0(R_10, mean_11) + t_10;
    float _S2357 = - in_opacity_10;
    float _S2358 = 1.0f + s_primal_ctx_exp_0(_S2357);
    float _S2359 = 1.0f / _S2358;
    float _S2360 = _S2358 * _S2358;
    float4  _S2361 = normalize_0(quat_10);
    float3  _S2362 = s_primal_ctx_exp_1(scale_10);
    float _S2363 = _S2361.y;
    float x2_10 = _S2363 * _S2363;
    float y2_10 = _S2361.z * _S2361.z;
    float z2_10 = _S2361.w * _S2361.w;
    float xy_10 = _S2361.y * _S2361.z;
    float xz_10 = _S2361.y * _S2361.w;
    float yz_10 = _S2361.z * _S2361.w;
    float wx_10 = _S2361.x * _S2361.y;
    float wy_10 = _S2361.x * _S2361.z;
    float wz_10 = _S2361.x * _S2361.w;
    Matrix<float, 3, 3>  _S2364 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_10), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_10), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2362.x, 0.0f, 0.0f, 0.0f, _S2362.y, 0.0f, 0.0f, 0.0f, _S2362.z);
    Matrix<float, 3, 3>  _S2365 = s_primal_ctx_mul_1(_S2364, S_2);
    Matrix<float, 3, 3>  _S2366 = transpose_3(_S2365);
    Matrix<float, 3, 3>  _S2367 = s_primal_ctx_mul_1(_S2365, _S2366);
    Matrix<float, 3, 3>  _S2368 = s_primal_ctx_mul_1(R_10, _S2367);
    Matrix<float, 3, 3>  _S2369 = transpose_3(R_10);
    Matrix<float, 3, 3>  _S2370 = s_primal_ctx_mul_1(_S2368, _S2369);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S2371 = float2 {mean_c_10.x, mean_c_10.y};
    float2  _S2372 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2373 = { _S2371, _S2372 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2374;
    (&_S2374)->primal_0 = _S2371;
    (&_S2374)->differential_0 = _S2372;
    DiffPair_float_0 _S2375 = s_primal_ctx_s_fwd_length_impl_0(&_S2374);
    float _S2376 = mean_c_10.z;
    DiffPair_float_0 _S2377 = { _S2375.primal_0, _S2375.differential_0 };
    DiffPair_float_0 _S2378 = { _S2376, 0.0f };
    DiffPair_float_0 _S2379;
    (&_S2379)->primal_0 = _S2375.primal_0;
    (&_S2379)->differential_0 = _S2375.differential_0;
    DiffPair_float_0 _S2380;
    (&_S2380)->primal_0 = _S2376;
    (&_S2380)->differential_0 = 0.0f;
    DiffPair_float_0 _S2381 = s_primal_ctx_d_atan2_0(&_S2379, &_S2380);
    bool _S2382 = (_S2375.primal_0) < 9.99999997475242708e-07f;
    float k_7;
    float s_diff_k_3;
    float _S2383;
    float _S2384;
    float _S2385;
    float _S2386;
    float _S2387;
    float _S2388;
    float _S2389;
    float _S2390;
    float _S2391;
    float _S2392;
    DiffPair_float_0 _S2393;
    if(_S2382)
    {
        float _S2394 = _S2381.differential_0 * _S2381.primal_0;
        float _S2395 = 1.0f - _S2381.primal_0 * _S2381.primal_0 / 24.0f;
        float _S2396 = 0.0f - (_S2394 + _S2394) * 0.0416666679084301f;
        float _S2397 = _S2376 * _S2376;
        float _S2398 = _S2396 * _S2376;
        float _S2399 = _S2398 / _S2397;
        float _S2400 = _S2397 * _S2397;
        k_7 = _S2395 / _S2376;
        s_diff_k_3 = _S2399;
        _S2383 = _S2400;
        _S2384 = _S2398;
        _S2385 = _S2397;
        _S2386 = _S2395;
        _S2387 = _S2396;
        _S2388 = 0.0f;
        _S2389 = 0.0f;
        _S2390 = 0.0f;
        _S2391 = 0.0f;
        _S2392 = 0.0f;
        (&_S2393)->primal_0 = 0.0f;
        (&_S2393)->differential_0 = 0.0f;
    }
    else
    {
        float _S2401 = 0.5f * _S2381.primal_0;
        float _S2402 = _S2381.differential_0 * 0.5f;
        DiffPair_float_0 _S2403;
        (&_S2403)->primal_0 = _S2401;
        (&_S2403)->differential_0 = _S2402;
        DiffPair_float_0 _S2404 = s_primal_ctx_d_sin_0(&_S2403);
        float _S2405 = 2.0f * _S2404.primal_0;
        float _S2406 = _S2404.differential_0 * 2.0f;
        float _S2407 = _S2375.primal_0 * _S2375.primal_0;
        float _S2408 = _S2406 * _S2375.primal_0 - _S2405 * _S2375.differential_0;
        float _S2409 = _S2408 / _S2407;
        float _S2410 = _S2407 * _S2407;
        k_7 = _S2405 / _S2375.primal_0;
        s_diff_k_3 = _S2409;
        _S2383 = 0.0f;
        _S2384 = 0.0f;
        _S2385 = 0.0f;
        _S2386 = 0.0f;
        _S2387 = 0.0f;
        _S2388 = _S2410;
        _S2389 = _S2408;
        _S2390 = _S2407;
        _S2391 = _S2405;
        _S2392 = _S2406;
        (&_S2393)->primal_0 = _S2401;
        (&_S2393)->differential_0 = _S2402;
    }
    float2  _S2411 = make_float2 (k_7);
    float2  _S2412 = make_float2 (s_diff_k_3);
    float2  _S2413 = _S2371 * make_float2 (k_7);
    float2  _S2414 = _S2372 * make_float2 (k_7) + make_float2 (s_diff_k_3) * _S2371;
    float u_80 = _S2413.x;
    float s_diff_u_24 = _S2414.x;
    float v_80 = _S2413.y;
    float s_diff_v_24 = _S2414.y;
    float _S2415 = s_diff_u_24 * u_80;
    float _S2416 = s_diff_v_24 * v_80;
    float r2_80 = u_80 * u_80 + v_80 * v_80;
    float s_diff_r2_24 = _S2415 + _S2415 + (_S2416 + _S2416);
    float _S2417 = s_diff_r2_24 * dist_coeffs_14[int(3)];
    float _S2418 = dist_coeffs_14[int(2)] + r2_80 * dist_coeffs_14[int(3)];
    float _S2419 = s_diff_r2_24 * _S2418 + _S2417 * r2_80;
    float _S2420 = dist_coeffs_14[int(1)] + r2_80 * _S2418;
    float _S2421 = s_diff_r2_24 * _S2420 + _S2419 * r2_80;
    float _S2422 = dist_coeffs_14[int(0)] + r2_80 * _S2420;
    float _S2423 = s_diff_r2_24 * _S2422 + _S2421 * r2_80;
    float2  _S2424 = make_float2 (_S2423);
    float radial_33 = 1.0f + r2_80 * _S2422;
    float2  _S2425 = make_float2 (radial_33);
    float _S2426 = 2.0f * dist_coeffs_14[int(4)];
    float _S2427 = _S2426 * u_80;
    float _S2428 = s_diff_u_24 * _S2426;
    float _S2429 = 2.0f * u_80;
    float _S2430 = s_diff_u_24 * 2.0f;
    float _S2431 = 2.0f * dist_coeffs_14[int(5)];
    float _S2432 = _S2431 * u_80;
    float _S2433 = s_diff_u_24 * _S2431;
    float _S2434 = 2.0f * v_80;
    float _S2435 = s_diff_v_24 * 2.0f;
    float2  _S2436 = _S2414 * make_float2 (radial_33) + make_float2 (_S2423) * _S2413 + make_float2 (_S2428 * v_80 + s_diff_v_24 * _S2427 + (s_diff_r2_24 + (_S2430 * u_80 + s_diff_u_24 * _S2429)) * dist_coeffs_14[int(5)] + s_diff_r2_24 * dist_coeffs_14[int(6)], _S2433 * v_80 + s_diff_v_24 * _S2432 + (s_diff_r2_24 + (_S2435 * v_80 + s_diff_v_24 * _S2434)) * dist_coeffs_14[int(4)] + s_diff_r2_24 * dist_coeffs_14[int(7)]);
    float2  _S2437 = _S2436 + make_float2 (_S2436.x * dist_coeffs_14[int(8)] + _S2436.y * dist_coeffs_14[int(9)], 0.0f);
    float _S2438 = _S2437.x * fx_14;
    float _S2439 = _S2437.y * fy_14;
    Matrix<float, 2, 3>  _S2440 = J_12;
    *&(((&_S2440)->rows + (int(0)))->x) = _S2438;
    *&(((&_S2440)->rows + (int(1)))->x) = _S2439;
    float2  _S2441 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2442 = { _S2371, _S2441 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2443;
    (&_S2443)->primal_0 = _S2371;
    (&_S2443)->differential_0 = _S2441;
    DiffPair_float_0 _S2444 = s_primal_ctx_s_fwd_length_impl_0(&_S2443);
    DiffPair_float_0 _S2445 = { _S2444.primal_0, _S2444.differential_0 };
    DiffPair_float_0 _S2446;
    (&_S2446)->primal_0 = _S2444.primal_0;
    (&_S2446)->differential_0 = _S2444.differential_0;
    DiffPair_float_0 _S2447;
    (&_S2447)->primal_0 = _S2376;
    (&_S2447)->differential_0 = 0.0f;
    DiffPair_float_0 _S2448 = s_primal_ctx_d_atan2_0(&_S2446, &_S2447);
    bool _S2449 = (_S2444.primal_0) < 9.99999997475242708e-07f;
    float _S2450;
    float _S2451;
    float _S2452;
    float _S2453;
    float _S2454;
    float _S2455;
    float _S2456;
    float _S2457;
    float _S2458;
    float _S2459;
    DiffPair_float_0 _S2460;
    if(_S2449)
    {
        float _S2461 = _S2448.differential_0 * _S2448.primal_0;
        float _S2462 = 1.0f - _S2448.primal_0 * _S2448.primal_0 / 24.0f;
        float _S2463 = 0.0f - (_S2461 + _S2461) * 0.0416666679084301f;
        float _S2464 = _S2376 * _S2376;
        float _S2465 = _S2463 * _S2376;
        float _S2466 = _S2465 / _S2464;
        float _S2467 = _S2464 * _S2464;
        k_7 = _S2462 / _S2376;
        s_diff_k_3 = _S2466;
        _S2450 = _S2467;
        _S2451 = _S2465;
        _S2452 = _S2464;
        _S2453 = _S2462;
        _S2454 = _S2463;
        _S2455 = 0.0f;
        _S2456 = 0.0f;
        _S2457 = 0.0f;
        _S2458 = 0.0f;
        _S2459 = 0.0f;
        (&_S2460)->primal_0 = 0.0f;
        (&_S2460)->differential_0 = 0.0f;
    }
    else
    {
        float _S2468 = 0.5f * _S2448.primal_0;
        float _S2469 = _S2448.differential_0 * 0.5f;
        DiffPair_float_0 _S2470;
        (&_S2470)->primal_0 = _S2468;
        (&_S2470)->differential_0 = _S2469;
        DiffPair_float_0 _S2471 = s_primal_ctx_d_sin_0(&_S2470);
        float _S2472 = 2.0f * _S2471.primal_0;
        float _S2473 = _S2471.differential_0 * 2.0f;
        float _S2474 = _S2444.primal_0 * _S2444.primal_0;
        float _S2475 = _S2473 * _S2444.primal_0 - _S2472 * _S2444.differential_0;
        float _S2476 = _S2475 / _S2474;
        float _S2477 = _S2474 * _S2474;
        k_7 = _S2472 / _S2444.primal_0;
        s_diff_k_3 = _S2476;
        _S2450 = 0.0f;
        _S2451 = 0.0f;
        _S2452 = 0.0f;
        _S2453 = 0.0f;
        _S2454 = 0.0f;
        _S2455 = _S2477;
        _S2456 = _S2475;
        _S2457 = _S2474;
        _S2458 = _S2472;
        _S2459 = _S2473;
        (&_S2460)->primal_0 = _S2468;
        (&_S2460)->differential_0 = _S2469;
    }
    float2  _S2478 = make_float2 (k_7);
    float2  _S2479 = make_float2 (s_diff_k_3);
    float2  _S2480 = _S2371 * make_float2 (k_7);
    float2  _S2481 = _S2441 * make_float2 (k_7) + make_float2 (s_diff_k_3) * _S2371;
    float u_81 = _S2480.x;
    float s_diff_u_25 = _S2481.x;
    float v_81 = _S2480.y;
    float s_diff_v_25 = _S2481.y;
    float _S2482 = s_diff_u_25 * u_81;
    float _S2483 = s_diff_v_25 * v_81;
    float r2_81 = u_81 * u_81 + v_81 * v_81;
    float s_diff_r2_25 = _S2482 + _S2482 + (_S2483 + _S2483);
    float _S2484 = s_diff_r2_25 * dist_coeffs_14[int(3)];
    float _S2485 = dist_coeffs_14[int(2)] + r2_81 * dist_coeffs_14[int(3)];
    float _S2486 = s_diff_r2_25 * _S2485 + _S2484 * r2_81;
    float _S2487 = dist_coeffs_14[int(1)] + r2_81 * _S2485;
    float _S2488 = s_diff_r2_25 * _S2487 + _S2486 * r2_81;
    float _S2489 = dist_coeffs_14[int(0)] + r2_81 * _S2487;
    float _S2490 = s_diff_r2_25 * _S2489 + _S2488 * r2_81;
    float2  _S2491 = make_float2 (_S2490);
    float radial_34 = 1.0f + r2_81 * _S2489;
    float2  _S2492 = make_float2 (radial_34);
    float _S2493 = _S2426 * u_81;
    float _S2494 = s_diff_u_25 * _S2426;
    float _S2495 = 2.0f * u_81;
    float _S2496 = s_diff_u_25 * 2.0f;
    float _S2497 = _S2431 * u_81;
    float _S2498 = s_diff_u_25 * _S2431;
    float _S2499 = 2.0f * v_81;
    float _S2500 = s_diff_v_25 * 2.0f;
    float2  _S2501 = _S2481 * make_float2 (radial_34) + make_float2 (_S2490) * _S2480 + make_float2 (_S2494 * v_81 + s_diff_v_25 * _S2493 + (s_diff_r2_25 + (_S2496 * u_81 + s_diff_u_25 * _S2495)) * dist_coeffs_14[int(5)] + s_diff_r2_25 * dist_coeffs_14[int(6)], _S2498 * v_81 + s_diff_v_25 * _S2497 + (s_diff_r2_25 + (_S2500 * v_81 + s_diff_v_25 * _S2499)) * dist_coeffs_14[int(4)] + s_diff_r2_25 * dist_coeffs_14[int(7)]);
    float2  _S2502 = _S2501 + make_float2 (_S2501.x * dist_coeffs_14[int(8)] + _S2501.y * dist_coeffs_14[int(9)], 0.0f);
    float _S2503 = _S2502.y * fy_14;
    *&(((&_S2440)->rows + (int(0)))->y) = _S2502.x * fx_14;
    *&(((&_S2440)->rows + (int(1)))->y) = _S2503;
    float2  _S2504 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2505 = { _S2371, _S2504 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2506;
    (&_S2506)->primal_0 = _S2371;
    (&_S2506)->differential_0 = _S2504;
    DiffPair_float_0 _S2507 = s_primal_ctx_s_fwd_length_impl_0(&_S2506);
    DiffPair_float_0 _S2508 = { _S2507.primal_0, _S2507.differential_0 };
    DiffPair_float_0 _S2509 = { _S2376, 1.0f };
    DiffPair_float_0 _S2510;
    (&_S2510)->primal_0 = _S2507.primal_0;
    (&_S2510)->differential_0 = _S2507.differential_0;
    DiffPair_float_0 _S2511;
    (&_S2511)->primal_0 = _S2376;
    (&_S2511)->differential_0 = 1.0f;
    DiffPair_float_0 _S2512 = s_primal_ctx_d_atan2_0(&_S2510, &_S2511);
    bool _S2513 = (_S2507.primal_0) < 9.99999997475242708e-07f;
    float _S2514;
    float _S2515;
    float _S2516;
    float _S2517;
    float _S2518;
    float _S2519;
    float _S2520;
    float _S2521;
    float _S2522;
    float _S2523;
    DiffPair_float_0 _S2524;
    if(_S2513)
    {
        float _S2525 = _S2512.differential_0 * _S2512.primal_0;
        float _S2526 = 1.0f - _S2512.primal_0 * _S2512.primal_0 / 24.0f;
        float _S2527 = 0.0f - (_S2525 + _S2525) * 0.0416666679084301f;
        float _S2528 = _S2376 * _S2376;
        float _S2529 = _S2527 * _S2376 - _S2526;
        float _S2530 = _S2529 / _S2528;
        float _S2531 = _S2528 * _S2528;
        k_7 = _S2526 / _S2376;
        s_diff_k_3 = _S2530;
        _S2514 = _S2531;
        _S2515 = _S2529;
        _S2516 = _S2528;
        _S2517 = _S2526;
        _S2518 = _S2527;
        _S2519 = 0.0f;
        _S2520 = 0.0f;
        _S2521 = 0.0f;
        _S2522 = 0.0f;
        _S2523 = 0.0f;
        (&_S2524)->primal_0 = 0.0f;
        (&_S2524)->differential_0 = 0.0f;
    }
    else
    {
        float _S2532 = 0.5f * _S2512.primal_0;
        float _S2533 = _S2512.differential_0 * 0.5f;
        DiffPair_float_0 _S2534;
        (&_S2534)->primal_0 = _S2532;
        (&_S2534)->differential_0 = _S2533;
        DiffPair_float_0 _S2535 = s_primal_ctx_d_sin_0(&_S2534);
        float _S2536 = 2.0f * _S2535.primal_0;
        float _S2537 = _S2535.differential_0 * 2.0f;
        float _S2538 = _S2507.primal_0 * _S2507.primal_0;
        float _S2539 = _S2537 * _S2507.primal_0 - _S2536 * _S2507.differential_0;
        float _S2540 = _S2539 / _S2538;
        float _S2541 = _S2538 * _S2538;
        k_7 = _S2536 / _S2507.primal_0;
        s_diff_k_3 = _S2540;
        _S2514 = 0.0f;
        _S2515 = 0.0f;
        _S2516 = 0.0f;
        _S2517 = 0.0f;
        _S2518 = 0.0f;
        _S2519 = _S2541;
        _S2520 = _S2539;
        _S2521 = _S2538;
        _S2522 = _S2536;
        _S2523 = _S2537;
        (&_S2524)->primal_0 = _S2532;
        (&_S2524)->differential_0 = _S2533;
    }
    float2  _S2542 = make_float2 (k_7);
    float2  _S2543 = make_float2 (s_diff_k_3);
    float2  _S2544 = _S2371 * make_float2 (k_7);
    float2  _S2545 = make_float2 (s_diff_k_3) * _S2371;
    float u_82 = _S2544.x;
    float s_diff_u_26 = _S2545.x;
    float v_82 = _S2544.y;
    float s_diff_v_26 = _S2545.y;
    float _S2546 = s_diff_u_26 * u_82;
    float _S2547 = s_diff_v_26 * v_82;
    float r2_82 = u_82 * u_82 + v_82 * v_82;
    float s_diff_r2_26 = _S2546 + _S2546 + (_S2547 + _S2547);
    float _S2548 = s_diff_r2_26 * dist_coeffs_14[int(3)];
    float _S2549 = dist_coeffs_14[int(2)] + r2_82 * dist_coeffs_14[int(3)];
    float _S2550 = s_diff_r2_26 * _S2549 + _S2548 * r2_82;
    float _S2551 = dist_coeffs_14[int(1)] + r2_82 * _S2549;
    float _S2552 = s_diff_r2_26 * _S2551 + _S2550 * r2_82;
    float _S2553 = dist_coeffs_14[int(0)] + r2_82 * _S2551;
    float _S2554 = s_diff_r2_26 * _S2553 + _S2552 * r2_82;
    float2  _S2555 = make_float2 (_S2554);
    float radial_35 = 1.0f + r2_82 * _S2553;
    float2  _S2556 = make_float2 (radial_35);
    float _S2557 = _S2426 * u_82;
    float _S2558 = s_diff_u_26 * _S2426;
    float _S2559 = 2.0f * u_82;
    float _S2560 = s_diff_u_26 * 2.0f;
    float _S2561 = _S2431 * u_82;
    float _S2562 = s_diff_u_26 * _S2431;
    float _S2563 = 2.0f * v_82;
    float _S2564 = s_diff_v_26 * 2.0f;
    float2  _S2565 = _S2545 * make_float2 (radial_35) + make_float2 (_S2554) * _S2544 + make_float2 (_S2558 * v_82 + s_diff_v_26 * _S2557 + (s_diff_r2_26 + (_S2560 * u_82 + s_diff_u_26 * _S2559)) * dist_coeffs_14[int(5)] + s_diff_r2_26 * dist_coeffs_14[int(6)], _S2562 * v_82 + s_diff_v_26 * _S2561 + (s_diff_r2_26 + (_S2564 * v_82 + s_diff_v_26 * _S2563)) * dist_coeffs_14[int(4)] + s_diff_r2_26 * dist_coeffs_14[int(7)]);
    float2  _S2566 = _S2565 + make_float2 (_S2565.x * dist_coeffs_14[int(8)] + _S2565.y * dist_coeffs_14[int(9)], 0.0f);
    float _S2567 = _S2566.y * fy_14;
    *&(((&_S2440)->rows + (int(0)))->z) = _S2566.x * fx_14;
    *&(((&_S2440)->rows + (int(1)))->z) = _S2567;
    Matrix<float, 2, 3>  _S2568 = s_primal_ctx_mul_2(_S2440, _S2370);
    Matrix<float, 3, 2>  _S2569 = transpose_1(_S2440);
    Matrix<float, 2, 2>  _S2570 = s_primal_ctx_mul_3(_S2568, _S2569);
    float eps2d_10;
    if(antialiased_10)
    {
        eps2d_10 = 0.10000000149011612f;
    }
    else
    {
        eps2d_10 = 0.30000001192092896f;
    }
    float _S2571 = _S2570.rows[int(0)].y * _S2570.rows[int(1)].x;
    float det_orig_10 = _S2570.rows[int(0)].x * _S2570.rows[int(1)].y - _S2571;
    float _S2572 = _S2570.rows[int(0)].x + eps2d_10;
    Matrix<float, 2, 2>  _S2573 = _S2570;
    *&(((&_S2573)->rows + (int(0)))->x) = _S2572;
    float _S2574 = _S2570.rows[int(1)].y + eps2d_10;
    *&(((&_S2573)->rows + (int(1)))->y) = _S2574;
    Matrix<float, 2, 2>  _S2575 = _S2573;
    Matrix<float, 2, 2>  _S2576 = _S2573;
    float det_blur_10 = _S2572 * _S2574 - _S2571;
    float _S2577 = det_orig_10 / det_blur_10;
    float _S2578 = det_blur_10 * det_blur_10;
    float _S2579 = (F32_max((0.0f), (_S2577)));
    float _S2580 = s_primal_ctx_sqrt_0(_S2579);
    float invdet_12 = 1.0f / det_blur_10;
    float _S2581 = - _S2570.rows[int(0)].y;
    float _S2582 = - _S2570.rows[int(1)].x;
    if(antialiased_10)
    {
        k_7 = _S2359 * _S2580;
    }
    else
    {
        k_7 = _S2359;
    }
    float _S2583 = k_7 / 0.00392156885936856f;
    float _S2584 = 2.0f * s_primal_ctx_log_0(_S2583);
    float _S2585 = s_primal_ctx_sqrt_0(_S2584);
    float _S2586 = _S2575.rows[int(0)].x;
    float _S2587 = _S2576.rows[int(1)].y;
    float3  campos_3 = - s_primal_ctx_mul_0(_S2369, t_10);
    float3  _S2588 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2589;
    (&_S2589)->primal_0 = mean_11;
    (&_S2589)->differential_0 = _S2588;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2590;
    (&_S2590)->primal_0 = scale_10;
    (&_S2590)->differential_0 = _S2588;
    DiffPair_float_0 _S2591;
    (&_S2591)->primal_0 = in_opacity_10;
    (&_S2591)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2592;
    (&_S2592)->primal_0 = campos_3;
    (&_S2592)->differential_0 = _S2588;
    s_bwd_prop_view_radius_3dgs_0(&_S2589, &_S2590, &_S2591, &_S2592, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2593 = _S2589;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2594 = _S2590;
    DiffPair_float_0 _S2595 = _S2591;
    float2  _S2596 = make_float2 (0.0f);
    float2  _S2597 = _S2596;
    *&((&_S2597)->y) = v_conic_2.z;
    float2  _S2598 = _S2596;
    *&((&_S2598)->y) = v_conic_2.y;
    *&((&_S2598)->x) = v_conic_2.x;
    DiffPair_float_0 _S2599;
    (&_S2599)->primal_0 = _S2587;
    (&_S2599)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2599, 0.0f);
    DiffPair_float_0 _S2600;
    (&_S2600)->primal_0 = _S2586;
    (&_S2600)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2600, 0.0f);
    DiffPair_float_0 _S2601;
    (&_S2601)->primal_0 = 3.32999992370605469f;
    (&_S2601)->differential_0 = 0.0f;
    DiffPair_float_0 _S2602;
    (&_S2602)->primal_0 = _S2585;
    (&_S2602)->differential_0 = 0.0f;
    _d_min_0(&_S2601, &_S2602, 0.0f);
    DiffPair_float_0 _S2603;
    (&_S2603)->primal_0 = _S2584;
    (&_S2603)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2603, _S2602.differential_0);
    float _S2604 = 2.0f * _S2603.differential_0;
    DiffPair_float_0 _S2605;
    (&_S2605)->primal_0 = _S2583;
    (&_S2605)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2605, _S2604);
    float _S2606 = v_opacity_2 + 254.9999847412109375f * _S2605.differential_0;
    float2  _S2607 = make_float2 (_S2600.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S2608 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2609 = _S2608;
    _S2609[int(1)] = _S2597;
    _S2609[int(0)] = _S2598;
    Matrix<float, 2, 2>  _S2610 = _S2609;
    float2  _S2611 = make_float2 (0.0f, _S2599.differential_0);
    if(antialiased_10)
    {
        float _S2612 = _S2580 * _S2606;
        k_7 = _S2359 * _S2606;
        s_diff_k_3 = _S2612;
    }
    else
    {
        k_7 = 0.0f;
        s_diff_k_3 = _S2606;
    }
    float _S2613 = invdet_12 * _S2610.rows[int(1)].y;
    float _S2614 = - (invdet_12 * _S2610.rows[int(1)].x);
    float _S2615 = - (invdet_12 * _S2610.rows[int(0)].y);
    float _S2616 = invdet_12 * _S2610.rows[int(0)].x;
    float _S2617 = - ((_S2572 * _S2610.rows[int(1)].y + _S2582 * _S2610.rows[int(1)].x + _S2581 * _S2610.rows[int(0)].y + _S2574 * _S2610.rows[int(0)].x) / _S2578);
    DiffPair_float_0 _S2618;
    (&_S2618)->primal_0 = _S2579;
    (&_S2618)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2618, k_7);
    DiffPair_float_0 _S2619;
    (&_S2619)->primal_0 = 0.0f;
    (&_S2619)->differential_0 = 0.0f;
    DiffPair_float_0 _S2620;
    (&_S2620)->primal_0 = _S2577;
    (&_S2620)->differential_0 = 0.0f;
    _d_max_0(&_S2619, &_S2620, _S2618.differential_0);
    float _S2621 = _S2620.differential_0 / _S2578;
    float s_diff_det_orig_T_2 = det_blur_10 * _S2621;
    float _S2622 = det_orig_10 * - _S2621 + _S2617;
    float _S2623 = - _S2622;
    float _S2624 = _S2572 * _S2622;
    float _S2625 = _S2574 * _S2622;
    Matrix<float, 2, 2>  _S2626 = _S2608;
    _S2626[int(1)] = _S2611;
    _S2626[int(0)] = _S2607;
    _S2573 = _S2626;
    *&(((&_S2573)->rows + (int(1)))->y) = 0.0f;
    float _S2627 = _S2624 + _S2626.rows[int(1)].y + _S2616;
    *&(((&_S2573)->rows + (int(0)))->x) = 0.0f;
    float _S2628 = _S2625 + _S2626.rows[int(0)].x + _S2613;
    float _S2629 = _S2623 + - s_diff_det_orig_T_2;
    float _S2630 = _S2570.rows[int(0)].y * _S2629 + _S2614;
    float _S2631 = _S2570.rows[int(1)].x * _S2629 + _S2615;
    float _S2632 = _S2570.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2633 = _S2627 + _S2570.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2634 = _S2596;
    *&((&_S2634)->x) = _S2630;
    *&((&_S2634)->y) = _S2633;
    float _S2635 = _S2628 + _S2632;
    float2  _S2636 = _S2596;
    *&((&_S2636)->y) = _S2631;
    *&((&_S2636)->x) = _S2635;
    Matrix<float, 2, 2>  _S2637 = _S2608;
    _S2637[int(1)] = _S2634;
    _S2637[int(0)] = _S2636;
    Matrix<float, 2, 2>  _S2638 = _S2573 + _S2637;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2639;
    (&_S2639)->primal_0 = _S2568;
    (&_S2639)->differential_0 = J_12;
    Matrix<float, 3, 2>  _S2640 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2641;
    (&_S2641)->primal_0 = _S2569;
    (&_S2641)->differential_0 = _S2640;
    s_bwd_prop_mul_0(&_S2639, &_S2641, _S2638);
    Matrix<float, 2, 3>  _S2642 = transpose_2(_S2641.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2643;
    (&_S2643)->primal_0 = _S2440;
    (&_S2643)->differential_0 = J_12;
    Matrix<float, 3, 3>  _S2644 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2645;
    (&_S2645)->primal_0 = _S2370;
    (&_S2645)->differential_0 = _S2644;
    s_bwd_prop_mul_1(&_S2643, &_S2645, _S2639.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2646 = _S2645;
    Matrix<float, 2, 3>  _S2647 = _S2642 + _S2643.differential_0;
    float2  _S2648 = make_float2 (0.0f, _S2647.rows[int(1)].z) + make_float2 (_S2647.rows[int(0)].z, 0.0f);
    float _S2649 = fx_14 * _S2648.x;
    float2  _S2650 = make_float2 (_S2649, fy_14 * _S2648.y) + make_float2 (dist_coeffs_14[int(8)] * _S2649, dist_coeffs_14[int(9)] * _S2649);
    float2  _S2651 = _S2544 * _S2650;
    float2  _S2652 = _S2545 * _S2650;
    float _S2653 = dist_coeffs_14[int(4)] * _S2650.y;
    float _S2654 = dist_coeffs_14[int(5)] * _S2650.x;
    float _S2655 = _S2652.x + _S2652.y;
    float _S2656 = _S2651.x + _S2651.y;
    float _S2657 = r2_82 * _S2656;
    float _S2658 = s_diff_r2_26 * _S2656 + r2_82 * _S2655;
    float _S2659 = r2_82 * _S2657;
    float _S2660 = s_diff_r2_26 * _S2657 + r2_82 * _S2658;
    float _S2661 = dist_coeffs_14[int(7)] * _S2650.y + _S2653 + dist_coeffs_14[int(6)] * _S2650.x + _S2654 + _S2553 * _S2656 + _S2551 * _S2657 + _S2549 * _S2659 + dist_coeffs_14[int(3)] * (r2_82 * _S2659);
    float _S2662 = _S2552 * _S2656 + _S2553 * _S2655 + _S2550 * _S2657 + _S2551 * _S2658 + _S2548 * _S2659 + _S2549 * _S2660 + dist_coeffs_14[int(3)] * (s_diff_r2_26 * _S2659 + r2_82 * _S2660);
    float _S2663 = _S2661 + _S2661;
    float _S2664 = v_82 * _S2662;
    float _S2665 = u_82 * _S2662;
    float2  _S2666 = _S2555 * _S2650 + make_float2 (_S2431 * (s_diff_v_26 * _S2650.y) + _S2560 * _S2654 + 2.0f * (s_diff_u_26 * _S2654) + _S2426 * (s_diff_v_26 * _S2650.x) + s_diff_u_26 * _S2663 + _S2665 + _S2665, _S2564 * _S2653 + 2.0f * (s_diff_v_26 * _S2653) + _S2562 * _S2650.y + _S2558 * _S2650.x + s_diff_v_26 * _S2663 + _S2664 + _S2664);
    float2  _S2667 = _S2556 * _S2650 + make_float2 (_S2431 * (v_82 * _S2650.y) + _S2559 * _S2654 + 2.0f * (u_82 * _S2654) + _S2426 * (v_82 * _S2650.x) + u_82 * _S2663, _S2563 * _S2653 + 2.0f * (v_82 * _S2653) + _S2561 * _S2650.y + _S2557 * _S2650.x + v_82 * _S2663);
    float2  _S2668 = _S2371 * _S2667;
    float2  _S2669 = _S2371 * _S2666;
    float _S2670 = _S2669.x + _S2669.y;
    float _S2671 = _S2668.x + _S2668.y;
    float2  _S2672 = _S2543 * _S2667 + _S2542 * _S2666;
    if(_S2513)
    {
        float _S2673 = _S2671 / _S2514;
        float _S2674 = _S2516 * _S2673;
        float _S2675 = _S2376 * (_S2515 * - _S2673);
        float _S2676 = _S2670 / _S2516;
        float _S2677 = 0.0416666679084301f * - (_S2376 * _S2674);
        float _S2678 = _S2677 + _S2677;
        float _S2679 = _S2512.primal_0 * (0.0416666679084301f * - (- _S2674 + _S2376 * _S2676));
        float _S2680 = _S2675 + _S2675 + _S2518 * _S2674 + _S2517 * - _S2676;
        float _S2681 = _S2512.differential_0 * _S2678 + _S2679 + _S2679;
        k_7 = _S2512.primal_0 * _S2678;
        _S2514 = _S2681;
        _S2515 = _S2680;
        _S2516 = 0.0f;
        _S2517 = 0.0f;
    }
    else
    {
        float _S2682 = _S2671 / _S2519;
        float _S2683 = _S2521 * _S2682;
        float _S2684 = _S2507.primal_0 * (_S2520 * - _S2682);
        float _S2685 = - _S2683;
        float _S2686 = _S2522 * _S2685;
        float _S2687 = _S2523 * _S2683;
        float _S2688 = _S2670 / _S2521;
        float _S2689 = _S2522 * - _S2688;
        float _S2690 = 2.0f * (_S2507.primal_0 * _S2683);
        float _S2691 = 2.0f * (_S2507.differential_0 * _S2685 + _S2507.primal_0 * _S2688);
        DiffPair_0 _S2692;
        (&_S2692)->primal_0 = _S2524;
        (&_S2692)->differential_0 = _S2356;
        DiffPair_float_0 _S2693;
        (&_S2693)->primal_0 = _S2691;
        (&_S2693)->differential_0 = _S2690;
        s_bwd_prop_d_sin_0(&_S2692, &_S2693);
        float _S2694 = 0.5f * _S2692.differential_0.primal_0;
        float _S2695 = _S2684 + _S2684 + _S2687 + _S2689;
        k_7 = 0.5f * _S2692.differential_0.differential_0;
        _S2514 = _S2694;
        _S2515 = 0.0f;
        _S2516 = _S2686;
        _S2517 = _S2695;
    }
    DiffPair_0 _S2696;
    (&_S2696)->primal_0 = _S2508;
    (&_S2696)->differential_0 = _S2356;
    DiffPair_0 _S2697;
    (&_S2697)->primal_0 = _S2509;
    (&_S2697)->differential_0 = _S2356;
    DiffPair_float_0 _S2698;
    (&_S2698)->primal_0 = _S2514;
    (&_S2698)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2696, &_S2697, &_S2698);
    float _S2699 = _S2697.differential_0.primal_0 + _S2515;
    float _S2700 = _S2696.differential_0.differential_0 + _S2516;
    float _S2701 = _S2696.differential_0.primal_0 + _S2517;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2702 = { _S2596, _S2596 };
    DiffPair_1 _S2703;
    (&_S2703)->primal_0 = _S2505;
    (&_S2703)->differential_0 = _S2702;
    DiffPair_float_0 _S2704;
    (&_S2704)->primal_0 = _S2701;
    (&_S2704)->differential_0 = _S2700;
    s_bwd_prop_s_fwd_length_impl_0(&_S2703, &_S2704);
    float2  _S2705 = _S2703.differential_0.primal_0 + _S2672;
    float3  _S2706 = make_float3 (_S2705.x, _S2705.y, _S2699);
    float2  _S2707 = make_float2 (0.0f, _S2647.rows[int(1)].y) + make_float2 (_S2647.rows[int(0)].y, 0.0f);
    float _S2708 = fx_14 * _S2707.x;
    float2  _S2709 = make_float2 (_S2708, fy_14 * _S2707.y) + make_float2 (dist_coeffs_14[int(8)] * _S2708, dist_coeffs_14[int(9)] * _S2708);
    float2  _S2710 = _S2480 * _S2709;
    float2  _S2711 = _S2481 * _S2709;
    float _S2712 = dist_coeffs_14[int(4)] * _S2709.y;
    float _S2713 = dist_coeffs_14[int(5)] * _S2709.x;
    float _S2714 = _S2711.x + _S2711.y;
    float _S2715 = _S2710.x + _S2710.y;
    float _S2716 = r2_81 * _S2715;
    float _S2717 = s_diff_r2_25 * _S2715 + r2_81 * _S2714;
    float _S2718 = r2_81 * _S2716;
    float _S2719 = s_diff_r2_25 * _S2716 + r2_81 * _S2717;
    float _S2720 = dist_coeffs_14[int(7)] * _S2709.y + _S2712 + dist_coeffs_14[int(6)] * _S2709.x + _S2713 + _S2489 * _S2715 + _S2487 * _S2716 + _S2485 * _S2718 + dist_coeffs_14[int(3)] * (r2_81 * _S2718);
    float _S2721 = _S2488 * _S2715 + _S2489 * _S2714 + _S2486 * _S2716 + _S2487 * _S2717 + _S2484 * _S2718 + _S2485 * _S2719 + dist_coeffs_14[int(3)] * (s_diff_r2_25 * _S2718 + r2_81 * _S2719);
    float _S2722 = _S2720 + _S2720;
    float _S2723 = v_81 * _S2721;
    float _S2724 = u_81 * _S2721;
    float2  _S2725 = _S2491 * _S2709 + make_float2 (_S2431 * (s_diff_v_25 * _S2709.y) + _S2496 * _S2713 + 2.0f * (s_diff_u_25 * _S2713) + _S2426 * (s_diff_v_25 * _S2709.x) + s_diff_u_25 * _S2722 + _S2724 + _S2724, _S2500 * _S2712 + 2.0f * (s_diff_v_25 * _S2712) + _S2498 * _S2709.y + _S2494 * _S2709.x + s_diff_v_25 * _S2722 + _S2723 + _S2723);
    float2  _S2726 = _S2492 * _S2709 + make_float2 (_S2431 * (v_81 * _S2709.y) + _S2495 * _S2713 + 2.0f * (u_81 * _S2713) + _S2426 * (v_81 * _S2709.x) + u_81 * _S2722, _S2499 * _S2712 + 2.0f * (v_81 * _S2712) + _S2497 * _S2709.y + _S2493 * _S2709.x + v_81 * _S2722);
    float2  _S2727 = _S2371 * _S2726;
    float2  _S2728 = _S2441 * _S2726;
    float2  _S2729 = _S2371 * _S2725;
    float _S2730 = _S2728.x + _S2728.y + _S2729.x + _S2729.y;
    float _S2731 = _S2727.x + _S2727.y;
    float2  _S2732 = _S2479 * _S2726 + _S2478 * _S2725;
    if(_S2449)
    {
        float _S2733 = _S2731 / _S2450;
        float _S2734 = _S2452 * _S2733;
        float _S2735 = _S2376 * (_S2451 * - _S2733);
        float _S2736 = _S2730 / _S2452;
        float _S2737 = 0.0416666679084301f * - (_S2376 * _S2734);
        float _S2738 = _S2737 + _S2737;
        float _S2739 = _S2448.primal_0 * (0.0416666679084301f * - (_S2376 * _S2736));
        float _S2740 = _S2735 + _S2735 + _S2454 * _S2734 + _S2453 * - _S2736;
        float _S2741 = _S2448.differential_0 * _S2738 + _S2739 + _S2739;
        k_7 = _S2448.primal_0 * _S2738;
        _S2450 = _S2741;
        _S2451 = _S2740;
        _S2452 = 0.0f;
        _S2453 = 0.0f;
    }
    else
    {
        float _S2742 = _S2731 / _S2455;
        float _S2743 = _S2457 * _S2742;
        float _S2744 = _S2444.primal_0 * (_S2456 * - _S2742);
        float _S2745 = - _S2743;
        float _S2746 = _S2458 * _S2745;
        float _S2747 = _S2459 * _S2743;
        float _S2748 = _S2730 / _S2457;
        float _S2749 = _S2458 * - _S2748;
        float _S2750 = 2.0f * (_S2444.primal_0 * _S2743);
        float _S2751 = 2.0f * (_S2444.differential_0 * _S2745 + _S2444.primal_0 * _S2748);
        DiffPair_0 _S2752;
        (&_S2752)->primal_0 = _S2460;
        (&_S2752)->differential_0 = _S2356;
        DiffPair_float_0 _S2753;
        (&_S2753)->primal_0 = _S2751;
        (&_S2753)->differential_0 = _S2750;
        s_bwd_prop_d_sin_0(&_S2752, &_S2753);
        float _S2754 = 0.5f * _S2752.differential_0.primal_0;
        float _S2755 = _S2744 + _S2744 + _S2747 + _S2749;
        k_7 = 0.5f * _S2752.differential_0.differential_0;
        _S2450 = _S2754;
        _S2451 = 0.0f;
        _S2452 = _S2746;
        _S2453 = _S2755;
    }
    DiffPair_0 _S2756;
    (&_S2756)->primal_0 = _S2445;
    (&_S2756)->differential_0 = _S2356;
    DiffPair_0 _S2757;
    (&_S2757)->primal_0 = _S2378;
    (&_S2757)->differential_0 = _S2356;
    DiffPair_float_0 _S2758;
    (&_S2758)->primal_0 = _S2450;
    (&_S2758)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2756, &_S2757, &_S2758);
    float _S2759 = _S2757.differential_0.primal_0 + _S2451;
    float _S2760 = _S2756.differential_0.differential_0 + _S2452;
    float _S2761 = _S2756.differential_0.primal_0 + _S2453;
    DiffPair_1 _S2762;
    (&_S2762)->primal_0 = _S2442;
    (&_S2762)->differential_0 = _S2702;
    DiffPair_float_0 _S2763;
    (&_S2763)->primal_0 = _S2761;
    (&_S2763)->differential_0 = _S2760;
    s_bwd_prop_s_fwd_length_impl_0(&_S2762, &_S2763);
    float2  _S2764 = _S2762.differential_0.primal_0 + _S2732;
    float2  _S2765 = make_float2 (0.0f, _S2647.rows[int(1)].x) + make_float2 (_S2647.rows[int(0)].x, 0.0f);
    float _S2766 = fx_14 * _S2765.x;
    float2  _S2767 = make_float2 (_S2766, fy_14 * _S2765.y) + make_float2 (dist_coeffs_14[int(8)] * _S2766, dist_coeffs_14[int(9)] * _S2766);
    float2  _S2768 = _S2413 * _S2767;
    float2  _S2769 = _S2414 * _S2767;
    float _S2770 = dist_coeffs_14[int(4)] * _S2767.y;
    float _S2771 = dist_coeffs_14[int(5)] * _S2767.x;
    float _S2772 = _S2769.x + _S2769.y;
    float _S2773 = _S2768.x + _S2768.y;
    float _S2774 = r2_80 * _S2773;
    float _S2775 = s_diff_r2_24 * _S2773 + r2_80 * _S2772;
    float _S2776 = r2_80 * _S2774;
    float _S2777 = s_diff_r2_24 * _S2774 + r2_80 * _S2775;
    float _S2778 = dist_coeffs_14[int(7)] * _S2767.y + _S2770 + dist_coeffs_14[int(6)] * _S2767.x + _S2771 + _S2422 * _S2773 + _S2420 * _S2774 + _S2418 * _S2776 + dist_coeffs_14[int(3)] * (r2_80 * _S2776);
    float _S2779 = _S2421 * _S2773 + _S2422 * _S2772 + _S2419 * _S2774 + _S2420 * _S2775 + _S2417 * _S2776 + _S2418 * _S2777 + dist_coeffs_14[int(3)] * (s_diff_r2_24 * _S2776 + r2_80 * _S2777);
    float _S2780 = _S2778 + _S2778;
    float _S2781 = v_80 * _S2779;
    float _S2782 = u_80 * _S2779;
    float2  _S2783 = _S2424 * _S2767 + make_float2 (_S2431 * (s_diff_v_24 * _S2767.y) + _S2430 * _S2771 + 2.0f * (s_diff_u_24 * _S2771) + _S2426 * (s_diff_v_24 * _S2767.x) + s_diff_u_24 * _S2780 + _S2782 + _S2782, _S2435 * _S2770 + 2.0f * (s_diff_v_24 * _S2770) + _S2433 * _S2767.y + _S2428 * _S2767.x + s_diff_v_24 * _S2780 + _S2781 + _S2781);
    float2  _S2784 = _S2425 * _S2767 + make_float2 (_S2431 * (v_80 * _S2767.y) + _S2429 * _S2771 + 2.0f * (u_80 * _S2771) + _S2426 * (v_80 * _S2767.x) + u_80 * _S2780, _S2434 * _S2770 + 2.0f * (v_80 * _S2770) + _S2432 * _S2767.y + _S2427 * _S2767.x + v_80 * _S2780);
    float3  _S2785 = make_float3 (_S2764.x, _S2764.y, _S2759) + _S2706;
    float2  _S2786 = _S2371 * _S2784;
    float2  _S2787 = _S2372 * _S2784;
    float2  _S2788 = _S2371 * _S2783;
    float _S2789 = _S2787.x + _S2787.y + _S2788.x + _S2788.y;
    float _S2790 = _S2786.x + _S2786.y;
    float2  _S2791 = _S2412 * _S2784 + _S2411 * _S2783;
    if(_S2382)
    {
        float _S2792 = _S2790 / _S2383;
        float _S2793 = _S2385 * _S2792;
        float _S2794 = _S2376 * (_S2384 * - _S2792);
        float _S2795 = _S2789 / _S2385;
        float _S2796 = 0.0416666679084301f * - (_S2376 * _S2793);
        float _S2797 = _S2796 + _S2796;
        float _S2798 = _S2381.primal_0 * (0.0416666679084301f * - (_S2376 * _S2795));
        float _S2799 = _S2794 + _S2794 + _S2387 * _S2793 + _S2386 * - _S2795;
        float _S2800 = _S2381.differential_0 * _S2797 + _S2798 + _S2798;
        k_7 = _S2381.primal_0 * _S2797;
        _S2383 = _S2800;
        _S2384 = _S2799;
        _S2385 = 0.0f;
        _S2386 = 0.0f;
    }
    else
    {
        float _S2801 = _S2790 / _S2388;
        float _S2802 = _S2390 * _S2801;
        float _S2803 = _S2375.primal_0 * (_S2389 * - _S2801);
        float _S2804 = - _S2802;
        float _S2805 = _S2391 * _S2804;
        float _S2806 = _S2392 * _S2802;
        float _S2807 = _S2789 / _S2390;
        float _S2808 = _S2391 * - _S2807;
        float _S2809 = 2.0f * (_S2375.primal_0 * _S2802);
        float _S2810 = 2.0f * (_S2375.differential_0 * _S2804 + _S2375.primal_0 * _S2807);
        DiffPair_0 _S2811;
        (&_S2811)->primal_0 = _S2393;
        (&_S2811)->differential_0 = _S2356;
        DiffPair_float_0 _S2812;
        (&_S2812)->primal_0 = _S2810;
        (&_S2812)->differential_0 = _S2809;
        s_bwd_prop_d_sin_0(&_S2811, &_S2812);
        float _S2813 = 0.5f * _S2811.differential_0.primal_0;
        float _S2814 = _S2803 + _S2803 + _S2806 + _S2808;
        k_7 = 0.5f * _S2811.differential_0.differential_0;
        _S2383 = _S2813;
        _S2384 = 0.0f;
        _S2385 = _S2805;
        _S2386 = _S2814;
    }
    DiffPair_0 _S2815;
    (&_S2815)->primal_0 = _S2377;
    (&_S2815)->differential_0 = _S2356;
    DiffPair_0 _S2816;
    (&_S2816)->primal_0 = _S2378;
    (&_S2816)->differential_0 = _S2356;
    DiffPair_float_0 _S2817;
    (&_S2817)->primal_0 = _S2383;
    (&_S2817)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2815, &_S2816, &_S2817);
    float _S2818 = _S2816.differential_0.primal_0 + _S2384;
    float _S2819 = _S2815.differential_0.differential_0 + _S2385;
    float _S2820 = _S2815.differential_0.primal_0 + _S2386;
    DiffPair_1 _S2821;
    (&_S2821)->primal_0 = _S2373;
    (&_S2821)->differential_0 = _S2702;
    DiffPair_float_0 _S2822;
    (&_S2822)->primal_0 = _S2820;
    (&_S2822)->differential_0 = _S2819;
    s_bwd_prop_s_fwd_length_impl_0(&_S2821, &_S2822);
    float2  _S2823 = _S2821.differential_0.primal_0 + _S2791;
    float3  _S2824 = make_float3 (_S2823.x, _S2823.y, _S2818);
    float _S2825 = length_0(_S2371);
    float _S2826 = s_primal_ctx_atan2_0(_S2825, _S2376);
    bool _S2827 = _S2825 < 9.99999997475242708e-07f;
    if(_S2827)
    {
        float _S2828 = 1.0f - _S2826 * _S2826 / 24.0f;
        float _S2829 = _S2376 * _S2376;
        k_7 = _S2828 / _S2376;
        _S2383 = _S2829;
        _S2384 = _S2828;
        _S2385 = 0.0f;
        _S2386 = 0.0f;
        _S2387 = 0.0f;
    }
    else
    {
        float _S2830 = 0.5f * _S2826;
        float _S2831 = 2.0f * s_primal_ctx_sin_0(_S2830);
        float _S2832 = _S2825 * _S2825;
        k_7 = _S2831 / _S2825;
        _S2383 = 0.0f;
        _S2384 = 0.0f;
        _S2385 = _S2832;
        _S2386 = _S2831;
        _S2387 = _S2830;
    }
    float2  _S2833 = make_float2 (k_7);
    float2  _S2834 = _S2371 * make_float2 (k_7);
    float _S2835 = fx_14 * v_mean2d_2.x;
    float u_83 = _S2834.x;
    float v_83 = _S2834.y;
    float r2_83 = u_83 * u_83 + v_83 * v_83;
    float _S2836 = dist_coeffs_14[int(2)] + r2_83 * dist_coeffs_14[int(3)];
    float _S2837 = dist_coeffs_14[int(1)] + r2_83 * _S2836;
    float _S2838 = dist_coeffs_14[int(0)] + r2_83 * _S2837;
    float2  _S2839 = make_float2 (_S2835, fy_14 * v_mean2d_2.y) + make_float2 (dist_coeffs_14[int(8)] * _S2835, dist_coeffs_14[int(9)] * _S2835);
    float2  _S2840 = _S2834 * _S2839;
    float _S2841 = dist_coeffs_14[int(4)] * _S2839.y;
    float _S2842 = dist_coeffs_14[int(5)] * _S2839.x;
    float _S2843 = _S2840.x + _S2840.y;
    float _S2844 = r2_83 * _S2843;
    float _S2845 = r2_83 * _S2844;
    float _S2846 = dist_coeffs_14[int(7)] * _S2839.y + _S2841 + dist_coeffs_14[int(6)] * _S2839.x + _S2842 + _S2838 * _S2843 + _S2837 * _S2844 + _S2836 * _S2845 + dist_coeffs_14[int(3)] * (r2_83 * _S2845);
    float _S2847 = v_83 * _S2846;
    float _S2848 = u_83 * _S2846;
    float2  _S2849 = make_float2 (1.0f + r2_83 * _S2838) * _S2839 + make_float2 (_S2431 * (v_83 * _S2839.y) + 2.0f * u_83 * _S2842 + 2.0f * (u_83 * _S2842) + _S2426 * (v_83 * _S2839.x) + _S2848 + _S2848, 2.0f * v_83 * _S2841 + 2.0f * (v_83 * _S2841) + _S2431 * u_83 * _S2839.y + _S2426 * u_83 * _S2839.x + _S2847 + _S2847);
    float2  _S2850 = _S2371 * _S2849;
    float2  _S2851 = _S2833 * _S2849;
    float _S2852 = _S2850.x + _S2850.y;
    if(_S2827)
    {
        float _S2853 = _S2852 / _S2383;
        float _S2854 = _S2384 * - _S2853;
        float _S2855 = _S2826 * (0.0416666679084301f * - (_S2376 * _S2853));
        k_7 = _S2855 + _S2855;
        _S2383 = _S2854;
        _S2384 = 0.0f;
    }
    else
    {
        float _S2856 = _S2852 / _S2385;
        float _S2857 = _S2386 * - _S2856;
        float _S2858 = 2.0f * (_S2825 * _S2856);
        DiffPair_float_0 _S2859;
        (&_S2859)->primal_0 = _S2387;
        (&_S2859)->differential_0 = 0.0f;
        s_bwd_prop_sin_0(&_S2859, _S2858);
        k_7 = 0.5f * _S2859.differential_0;
        _S2383 = 0.0f;
        _S2384 = _S2857;
    }
    DiffPair_float_0 _S2860;
    (&_S2860)->primal_0 = _S2825;
    (&_S2860)->differential_0 = 0.0f;
    DiffPair_float_0 _S2861;
    (&_S2861)->primal_0 = _S2376;
    (&_S2861)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2860, &_S2861, k_7);
    float _S2862 = _S2861.differential_0 + _S2383;
    float _S2863 = _S2860.differential_0 + _S2384;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2864;
    (&_S2864)->primal_0 = _S2371;
    (&_S2864)->differential_0 = _S2596;
    s_bwd_length_impl_2(&_S2864, _S2863);
    float2  _S2865 = _S2864.differential_0 + _S2851;
    float3  _S2866 = make_float3 (_S2865.x, _S2865.y, _S2862);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2867;
    (&_S2867)->primal_0 = _S2368;
    (&_S2867)->differential_0 = _S2644;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2868;
    (&_S2868)->primal_0 = _S2369;
    (&_S2868)->differential_0 = _S2644;
    s_bwd_prop_mul_2(&_S2867, &_S2868, _S2646.differential_0);
    Matrix<float, 3, 3>  _S2869 = transpose_3(_S2868.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2870;
    (&_S2870)->primal_0 = R_10;
    (&_S2870)->differential_0 = _S2644;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2871;
    (&_S2871)->primal_0 = _S2367;
    (&_S2871)->differential_0 = _S2644;
    s_bwd_prop_mul_2(&_S2870, &_S2871, _S2867.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2872;
    (&_S2872)->primal_0 = _S2365;
    (&_S2872)->differential_0 = _S2644;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2873;
    (&_S2873)->primal_0 = _S2366;
    (&_S2873)->differential_0 = _S2644;
    s_bwd_prop_mul_2(&_S2872, &_S2873, _S2871.differential_0);
    Matrix<float, 3, 3>  _S2874 = _S2872.differential_0 + transpose_3(_S2873.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2875;
    (&_S2875)->primal_0 = _S2364;
    (&_S2875)->differential_0 = _S2644;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2876;
    (&_S2876)->primal_0 = S_2;
    (&_S2876)->differential_0 = _S2644;
    s_bwd_prop_mul_2(&_S2875, &_S2876, _S2874);
    Matrix<float, 3, 3>  _S2877 = transpose_3(_S2875.differential_0);
    float _S2878 = 2.0f * - _S2877.rows[int(2)].z;
    float _S2879 = 2.0f * _S2877.rows[int(2)].y;
    float _S2880 = 2.0f * _S2877.rows[int(2)].x;
    float _S2881 = 2.0f * _S2877.rows[int(1)].z;
    float _S2882 = 2.0f * - _S2877.rows[int(1)].y;
    float _S2883 = 2.0f * _S2877.rows[int(1)].x;
    float _S2884 = 2.0f * _S2877.rows[int(0)].z;
    float _S2885 = 2.0f * _S2877.rows[int(0)].y;
    float _S2886 = 2.0f * - _S2877.rows[int(0)].x;
    float _S2887 = - _S2883 + _S2885;
    float _S2888 = _S2880 + - _S2884;
    float _S2889 = - _S2879 + _S2881;
    float _S2890 = _S2879 + _S2881;
    float _S2891 = _S2880 + _S2884;
    float _S2892 = _S2883 + _S2885;
    float _S2893 = _S2361.w * (_S2882 + _S2886);
    float _S2894 = _S2361.z * (_S2878 + _S2886);
    float _S2895 = _S2361.y * (_S2878 + _S2882);
    float _S2896 = _S2361.x * _S2887 + _S2361.z * _S2890 + _S2361.y * _S2891 + _S2893 + _S2893;
    float _S2897 = _S2361.x * _S2888 + _S2361.w * _S2890 + _S2361.y * _S2892 + _S2894 + _S2894;
    float _S2898 = _S2361.x * _S2889 + _S2361.w * _S2891 + _S2361.z * _S2892 + _S2895 + _S2895;
    float _S2899 = _S2361.w * _S2887 + _S2361.z * _S2888 + _S2361.y * _S2889;
    float3  _S2900 = _S2588;
    *&((&_S2900)->z) = _S2876.differential_0.rows[int(2)].z;
    *&((&_S2900)->y) = _S2876.differential_0.rows[int(1)].y;
    *&((&_S2900)->x) = _S2876.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2901;
    (&_S2901)->primal_0 = scale_10;
    (&_S2901)->differential_0 = _S2588;
    s_bwd_prop_exp_1(&_S2901, _S2900);
    float4  _S2902 = make_float4 (0.0f);
    float4  _S2903 = _S2902;
    *&((&_S2903)->w) = _S2896;
    *&((&_S2903)->z) = _S2897;
    *&((&_S2903)->y) = _S2898;
    *&((&_S2903)->x) = _S2899;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2904;
    (&_S2904)->primal_0 = quat_10;
    (&_S2904)->differential_0 = _S2902;
    s_bwd_normalize_impl_0(&_S2904, _S2903);
    float _S2905 = - (s_diff_k_3 / _S2360);
    DiffPair_float_0 _S2906;
    (&_S2906)->primal_0 = _S2357;
    (&_S2906)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2906, _S2905);
    float _S2907 = - _S2906.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2908;
    (&_S2908)->primal_0 = mean_c_10;
    (&_S2908)->differential_0 = _S2588;
    s_bwd_length_impl_0(&_S2908, v_depth_2);
    float3  _S2909 = _S2824 + _S2866 + _S2908.differential_0 + _S2785;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2910;
    (&_S2910)->primal_0 = R_10;
    (&_S2910)->differential_0 = _S2644;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2911;
    (&_S2911)->primal_0 = mean_11;
    (&_S2911)->differential_0 = _S2588;
    s_bwd_prop_mul_3(&_S2910, &_S2911, _S2909);
    Matrix<float, 3, 3>  _S2912 = _S2869 + _S2870.differential_0 + _S2910.differential_0;
    float _S2913 = _S2907 + _S2595.differential_0;
    float3  _S2914 = _S2901.differential_0 + _S2594.differential_0;
    *v_mean_2 = *v_mean_2 + (_S2911.differential_0 + _S2593.differential_0);
    *v_quat_2 = *v_quat_2 + _S2904.differential_0;
    *v_scale_2 = *v_scale_2 + _S2914;
    *v_in_opacity_2 = *v_in_opacity_2 + _S2913;
    *v_R_2 = *v_R_2 + _S2912;
    *v_t_2 = *v_t_2 + _S2909;
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2915, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2916, float _S2917)
{
    _d_dot_0(_S2915, _S2916, _S2917);
    return;
}

inline __device__ void projection_3dgs_equirect_vjp(bool antialiased_11, float3  mean_12, float4  quat_11, float3  scale_11, float in_opacity_11, Matrix<float, 3, 3>  R_11, float3  t_11, float fx_15, float fy_15, float cx_12, float cy_12, FixedArray<float, 10>  dist_coeffs_15, uint image_width_11, uint image_height_11, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_11 = s_primal_ctx_mul_0(R_11, mean_12) + t_11;
    float _S2918 = - in_opacity_11;
    float _S2919 = 1.0f + s_primal_ctx_exp_0(_S2918);
    float _S2920 = 1.0f / _S2919;
    float _S2921 = _S2919 * _S2919;
    float4  _S2922 = normalize_0(quat_11);
    float3  _S2923 = s_primal_ctx_exp_1(scale_11);
    float _S2924 = _S2922.y;
    float x2_11 = _S2924 * _S2924;
    float y2_11 = _S2922.z * _S2922.z;
    float z2_11 = _S2922.w * _S2922.w;
    float xy_11 = _S2922.y * _S2922.z;
    float xz_11 = _S2922.y * _S2922.w;
    float yz_11 = _S2922.z * _S2922.w;
    float wx_11 = _S2922.x * _S2922.y;
    float wy_11 = _S2922.x * _S2922.z;
    float wz_11 = _S2922.x * _S2922.w;
    Matrix<float, 3, 3>  _S2925 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_11), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_11), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S2923.x, 0.0f, 0.0f, 0.0f, _S2923.y, 0.0f, 0.0f, 0.0f, _S2923.z);
    Matrix<float, 3, 3>  _S2926 = s_primal_ctx_mul_1(_S2925, S_3);
    Matrix<float, 3, 3>  _S2927 = transpose_3(_S2926);
    Matrix<float, 3, 3>  _S2928 = s_primal_ctx_mul_1(_S2926, _S2927);
    Matrix<float, 3, 3>  _S2929 = s_primal_ctx_mul_1(R_11, _S2928);
    Matrix<float, 3, 3>  _S2930 = transpose_3(R_11);
    Matrix<float, 3, 3>  _S2931 = s_primal_ctx_mul_1(_S2929, _S2930);
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (0.0f);
    float _S2932 = mean_c_11.x;
    float _S2933 = mean_c_11.z;
    DiffPair_float_0 _S2934 = { _S2932, 1.0f };
    DiffPair_float_0 _S2935 = { _S2933, 0.0f };
    DiffPair_float_0 _S2936;
    (&_S2936)->primal_0 = _S2932;
    (&_S2936)->differential_0 = 1.0f;
    DiffPair_float_0 _S2937;
    (&_S2937)->primal_0 = _S2933;
    (&_S2937)->differential_0 = 0.0f;
    DiffPair_float_0 _S2938 = s_primal_ctx_d_atan2_0(&_S2936, &_S2937);
    float _S2939 = mean_c_11.y;
    float2  _S2940 = float2 {mean_c_11.x, mean_c_11.z};
    float2  _S2941 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2942 = { _S2940, _S2941 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2943;
    (&_S2943)->primal_0 = _S2940;
    (&_S2943)->differential_0 = _S2941;
    DiffPair_float_0 _S2944 = s_primal_ctx_s_fwd_length_impl_0(&_S2943);
    DiffPair_float_0 _S2945 = { _S2939, 0.0f };
    DiffPair_float_0 _S2946 = { _S2944.primal_0, _S2944.differential_0 };
    DiffPair_float_0 _S2947;
    (&_S2947)->primal_0 = _S2939;
    (&_S2947)->differential_0 = 0.0f;
    DiffPair_float_0 _S2948;
    (&_S2948)->primal_0 = _S2944.primal_0;
    (&_S2948)->differential_0 = _S2944.differential_0;
    DiffPair_float_0 _S2949 = s_primal_ctx_d_atan2_0(&_S2947, &_S2948);
    float _S2950 = _S2938.differential_0 * fx_15;
    float _S2951 = _S2949.differential_0 * fy_15;
    Matrix<float, 2, 3>  _S2952 = J_13;
    *&(((&_S2952)->rows + (int(0)))->x) = _S2950;
    *&(((&_S2952)->rows + (int(1)))->x) = _S2951;
    DiffPair_float_0 _S2953 = { _S2932, 0.0f };
    DiffPair_float_0 _S2954;
    (&_S2954)->primal_0 = _S2932;
    (&_S2954)->differential_0 = 0.0f;
    DiffPair_float_0 _S2955;
    (&_S2955)->primal_0 = _S2933;
    (&_S2955)->differential_0 = 0.0f;
    DiffPair_float_0 _S2956 = s_primal_ctx_d_atan2_0(&_S2954, &_S2955);
    float2  _S2957 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2958 = { _S2940, _S2957 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2959;
    (&_S2959)->primal_0 = _S2940;
    (&_S2959)->differential_0 = _S2957;
    DiffPair_float_0 _S2960 = s_primal_ctx_s_fwd_length_impl_0(&_S2959);
    DiffPair_float_0 _S2961 = { _S2939, 1.0f };
    DiffPair_float_0 _S2962 = { _S2960.primal_0, _S2960.differential_0 };
    DiffPair_float_0 _S2963;
    (&_S2963)->primal_0 = _S2939;
    (&_S2963)->differential_0 = 1.0f;
    DiffPair_float_0 _S2964;
    (&_S2964)->primal_0 = _S2960.primal_0;
    (&_S2964)->differential_0 = _S2960.differential_0;
    DiffPair_float_0 _S2965 = s_primal_ctx_d_atan2_0(&_S2963, &_S2964);
    float _S2966 = _S2965.differential_0 * fy_15;
    *&(((&_S2952)->rows + (int(0)))->y) = _S2956.differential_0 * fx_15;
    *&(((&_S2952)->rows + (int(1)))->y) = _S2966;
    DiffPair_float_0 _S2967 = { _S2933, 1.0f };
    DiffPair_float_0 _S2968;
    (&_S2968)->primal_0 = _S2932;
    (&_S2968)->differential_0 = 0.0f;
    DiffPair_float_0 _S2969;
    (&_S2969)->primal_0 = _S2933;
    (&_S2969)->differential_0 = 1.0f;
    DiffPair_float_0 _S2970 = s_primal_ctx_d_atan2_0(&_S2968, &_S2969);
    float2  _S2971 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2972 = { _S2940, _S2971 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2973;
    (&_S2973)->primal_0 = _S2940;
    (&_S2973)->differential_0 = _S2971;
    DiffPair_float_0 _S2974 = s_primal_ctx_s_fwd_length_impl_0(&_S2973);
    DiffPair_float_0 _S2975 = { _S2974.primal_0, _S2974.differential_0 };
    DiffPair_float_0 _S2976;
    (&_S2976)->primal_0 = _S2939;
    (&_S2976)->differential_0 = 0.0f;
    DiffPair_float_0 _S2977;
    (&_S2977)->primal_0 = _S2974.primal_0;
    (&_S2977)->differential_0 = _S2974.differential_0;
    DiffPair_float_0 _S2978 = s_primal_ctx_d_atan2_0(&_S2976, &_S2977);
    float _S2979 = _S2978.differential_0 * fy_15;
    *&(((&_S2952)->rows + (int(0)))->z) = _S2970.differential_0 * fx_15;
    *&(((&_S2952)->rows + (int(1)))->z) = _S2979;
    Matrix<float, 2, 3>  _S2980 = s_primal_ctx_mul_2(_S2952, _S2931);
    Matrix<float, 3, 2>  _S2981 = transpose_1(_S2952);
    Matrix<float, 2, 2>  _S2982 = s_primal_ctx_mul_3(_S2980, _S2981);
    float eps2d_11;
    if(antialiased_11)
    {
        eps2d_11 = 0.10000000149011612f;
    }
    else
    {
        eps2d_11 = 0.30000001192092896f;
    }
    float _S2983 = _S2982.rows[int(0)].y * _S2982.rows[int(1)].x;
    float det_orig_11 = _S2982.rows[int(0)].x * _S2982.rows[int(1)].y - _S2983;
    float _S2984 = _S2982.rows[int(0)].x + eps2d_11;
    Matrix<float, 2, 2>  _S2985 = _S2982;
    *&(((&_S2985)->rows + (int(0)))->x) = _S2984;
    float _S2986 = _S2982.rows[int(1)].y + eps2d_11;
    *&(((&_S2985)->rows + (int(1)))->y) = _S2986;
    Matrix<float, 2, 2>  _S2987 = _S2985;
    Matrix<float, 2, 2>  _S2988 = _S2985;
    float det_blur_11 = _S2984 * _S2986 - _S2983;
    float _S2989 = det_orig_11 / det_blur_11;
    float _S2990 = det_blur_11 * det_blur_11;
    float _S2991 = (F32_max((0.0f), (_S2989)));
    float _S2992 = s_primal_ctx_sqrt_0(_S2991);
    float invdet_13 = 1.0f / det_blur_11;
    float _S2993 = - _S2982.rows[int(0)].y;
    float _S2994 = - _S2982.rows[int(1)].x;
    if(antialiased_11)
    {
        eps2d_11 = _S2920 * _S2992;
    }
    else
    {
        eps2d_11 = _S2920;
    }
    float _S2995 = eps2d_11 / 0.00392156885936856f;
    float _S2996 = 2.0f * s_primal_ctx_log_0(_S2995);
    float _S2997 = s_primal_ctx_sqrt_0(_S2996);
    float _S2998 = _S2987.rows[int(0)].x;
    float _S2999 = _S2988.rows[int(1)].y;
    float3  campos_4 = - s_primal_ctx_mul_0(_S2930, t_11);
    float3  _S3000 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3001;
    (&_S3001)->primal_0 = mean_12;
    (&_S3001)->differential_0 = _S3000;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3002;
    (&_S3002)->primal_0 = scale_11;
    (&_S3002)->differential_0 = _S3000;
    DiffPair_float_0 _S3003;
    (&_S3003)->primal_0 = in_opacity_11;
    (&_S3003)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3004;
    (&_S3004)->primal_0 = campos_4;
    (&_S3004)->differential_0 = _S3000;
    s_bwd_prop_view_radius_3dgs_0(&_S3001, &_S3002, &_S3003, &_S3004, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3005 = _S3001;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3006 = _S3002;
    DiffPair_float_0 _S3007 = _S3003;
    float2  _S3008 = make_float2 (0.0f);
    float2  _S3009 = _S3008;
    *&((&_S3009)->y) = v_conic_3.z;
    float2  _S3010 = _S3008;
    *&((&_S3010)->y) = v_conic_3.y;
    *&((&_S3010)->x) = v_conic_3.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3011;
    (&_S3011)->primal_0 = mean_c_11;
    (&_S3011)->differential_0 = _S3000;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3012;
    (&_S3012)->primal_0 = mean_c_11;
    (&_S3012)->differential_0 = _S3000;
    s_bwd_prop_dot_0(&_S3011, &_S3012, 0.0f);
    DiffPair_float_0 _S3013;
    (&_S3013)->primal_0 = _S2999;
    (&_S3013)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3013, 0.0f);
    DiffPair_float_0 _S3014;
    (&_S3014)->primal_0 = _S2998;
    (&_S3014)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3014, 0.0f);
    DiffPair_float_0 _S3015;
    (&_S3015)->primal_0 = 3.32999992370605469f;
    (&_S3015)->differential_0 = 0.0f;
    DiffPair_float_0 _S3016;
    (&_S3016)->primal_0 = _S2997;
    (&_S3016)->differential_0 = 0.0f;
    _d_min_0(&_S3015, &_S3016, 0.0f);
    DiffPair_float_0 _S3017;
    (&_S3017)->primal_0 = _S2996;
    (&_S3017)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3017, _S3016.differential_0);
    float _S3018 = 2.0f * _S3017.differential_0;
    DiffPair_float_0 _S3019;
    (&_S3019)->primal_0 = _S2995;
    (&_S3019)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3019, _S3018);
    float _S3020 = v_opacity_3 + 254.9999847412109375f * _S3019.differential_0;
    float2  _S3021 = make_float2 (_S3014.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S3022 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3023 = _S3022;
    _S3023[int(1)] = _S3009;
    _S3023[int(0)] = _S3010;
    Matrix<float, 2, 2>  _S3024 = _S3023;
    float3  _S3025 = _S3012.differential_0 + _S3011.differential_0;
    float2  _S3026 = make_float2 (0.0f, _S3013.differential_0);
    float _S3027;
    if(antialiased_11)
    {
        float _S3028 = _S2992 * _S3020;
        eps2d_11 = _S2920 * _S3020;
        _S3027 = _S3028;
    }
    else
    {
        eps2d_11 = 0.0f;
        _S3027 = _S3020;
    }
    float _S3029 = invdet_13 * _S3024.rows[int(1)].y;
    float _S3030 = - (invdet_13 * _S3024.rows[int(1)].x);
    float _S3031 = - (invdet_13 * _S3024.rows[int(0)].y);
    float _S3032 = invdet_13 * _S3024.rows[int(0)].x;
    float _S3033 = - ((_S2984 * _S3024.rows[int(1)].y + _S2994 * _S3024.rows[int(1)].x + _S2993 * _S3024.rows[int(0)].y + _S2986 * _S3024.rows[int(0)].x) / _S2990);
    DiffPair_float_0 _S3034;
    (&_S3034)->primal_0 = _S2991;
    (&_S3034)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3034, eps2d_11);
    DiffPair_float_0 _S3035 = { 0.0f, 0.0f };
    DiffPair_float_0 _S3036;
    (&_S3036)->primal_0 = 0.0f;
    (&_S3036)->differential_0 = 0.0f;
    DiffPair_float_0 _S3037;
    (&_S3037)->primal_0 = _S2989;
    (&_S3037)->differential_0 = 0.0f;
    _d_max_0(&_S3036, &_S3037, _S3034.differential_0);
    float _S3038 = _S3037.differential_0 / _S2990;
    float s_diff_det_orig_T_3 = det_blur_11 * _S3038;
    float _S3039 = det_orig_11 * - _S3038 + _S3033;
    float _S3040 = - _S3039;
    float _S3041 = _S2984 * _S3039;
    float _S3042 = _S2986 * _S3039;
    Matrix<float, 2, 2>  _S3043 = _S3022;
    _S3043[int(1)] = _S3026;
    _S3043[int(0)] = _S3021;
    _S2985 = _S3043;
    *&(((&_S2985)->rows + (int(1)))->y) = 0.0f;
    float _S3044 = _S3041 + _S3043.rows[int(1)].y + _S3032;
    *&(((&_S2985)->rows + (int(0)))->x) = 0.0f;
    float _S3045 = _S3042 + _S3043.rows[int(0)].x + _S3029;
    float _S3046 = _S3040 + - s_diff_det_orig_T_3;
    float _S3047 = _S2982.rows[int(0)].y * _S3046 + _S3030;
    float _S3048 = _S2982.rows[int(1)].x * _S3046 + _S3031;
    float _S3049 = _S2982.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S3050 = _S3044 + _S2982.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S3051 = _S3008;
    *&((&_S3051)->x) = _S3047;
    *&((&_S3051)->y) = _S3050;
    float _S3052 = _S3045 + _S3049;
    float2  _S3053 = _S3008;
    *&((&_S3053)->y) = _S3048;
    *&((&_S3053)->x) = _S3052;
    Matrix<float, 2, 2>  _S3054 = _S3022;
    _S3054[int(1)] = _S3051;
    _S3054[int(0)] = _S3053;
    Matrix<float, 2, 2>  _S3055 = _S2985 + _S3054;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S3056;
    (&_S3056)->primal_0 = _S2980;
    (&_S3056)->differential_0 = J_13;
    Matrix<float, 3, 2>  _S3057 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S3058;
    (&_S3058)->primal_0 = _S2981;
    (&_S3058)->differential_0 = _S3057;
    s_bwd_prop_mul_0(&_S3056, &_S3058, _S3055);
    Matrix<float, 2, 3>  _S3059 = transpose_2(_S3058.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S3060;
    (&_S3060)->primal_0 = _S2952;
    (&_S3060)->differential_0 = J_13;
    Matrix<float, 3, 3>  _S3061 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3062;
    (&_S3062)->primal_0 = _S2931;
    (&_S3062)->differential_0 = _S3061;
    s_bwd_prop_mul_1(&_S3060, &_S3062, _S3056.differential_0);
    Matrix<float, 2, 3>  _S3063 = _S3059 + _S3060.differential_0;
    float2  _S3064 = make_float2 (0.0f, _S3063.rows[int(1)].z) + make_float2 (_S3063.rows[int(0)].z, 0.0f);
    float _S3065 = fy_15 * _S3064.y;
    float _S3066 = fx_15 * _S3064.x;
    DiffPair_0 _S3067;
    (&_S3067)->primal_0 = _S2945;
    (&_S3067)->differential_0 = _S3035;
    DiffPair_0 _S3068;
    (&_S3068)->primal_0 = _S2975;
    (&_S3068)->differential_0 = _S3035;
    DiffPair_float_0 _S3069;
    (&_S3069)->primal_0 = 0.0f;
    (&_S3069)->differential_0 = _S3065;
    s_bwd_prop_d_atan2_0(&_S3067, &_S3068, &_S3069);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3070 = { _S3008, _S3008 };
    DiffPair_1 _S3071;
    (&_S3071)->primal_0 = _S2972;
    (&_S3071)->differential_0 = _S3070;
    DiffPair_float_0 _S3072;
    (&_S3072)->primal_0 = _S3068.differential_0.primal_0;
    (&_S3072)->differential_0 = _S3068.differential_0.differential_0;
    s_bwd_prop_s_fwd_length_impl_0(&_S3071, &_S3072);
    DiffPair_0 _S3073;
    (&_S3073)->primal_0 = _S2953;
    (&_S3073)->differential_0 = _S3035;
    DiffPair_0 _S3074;
    (&_S3074)->primal_0 = _S2967;
    (&_S3074)->differential_0 = _S3035;
    DiffPair_float_0 _S3075;
    (&_S3075)->primal_0 = 0.0f;
    (&_S3075)->differential_0 = _S3066;
    s_bwd_prop_d_atan2_0(&_S3073, &_S3074, &_S3075);
    float3  _S3076 = make_float3 (_S3071.differential_0.primal_0.x + _S3073.differential_0.primal_0, _S3067.differential_0.primal_0, _S3071.differential_0.primal_0.y + _S3074.differential_0.primal_0);
    float2  _S3077 = make_float2 (0.0f, _S3063.rows[int(1)].y) + make_float2 (_S3063.rows[int(0)].y, 0.0f);
    float _S3078 = fy_15 * _S3077.y;
    float _S3079 = fx_15 * _S3077.x;
    DiffPair_0 _S3080;
    (&_S3080)->primal_0 = _S2961;
    (&_S3080)->differential_0 = _S3035;
    DiffPair_0 _S3081;
    (&_S3081)->primal_0 = _S2962;
    (&_S3081)->differential_0 = _S3035;
    DiffPair_float_0 _S3082;
    (&_S3082)->primal_0 = 0.0f;
    (&_S3082)->differential_0 = _S3078;
    s_bwd_prop_d_atan2_0(&_S3080, &_S3081, &_S3082);
    DiffPair_1 _S3083;
    (&_S3083)->primal_0 = _S2958;
    (&_S3083)->differential_0 = _S3070;
    DiffPair_float_0 _S3084;
    (&_S3084)->primal_0 = _S3081.differential_0.primal_0;
    (&_S3084)->differential_0 = _S3081.differential_0.differential_0;
    s_bwd_prop_s_fwd_length_impl_0(&_S3083, &_S3084);
    DiffPair_0 _S3085;
    (&_S3085)->primal_0 = _S2953;
    (&_S3085)->differential_0 = _S3035;
    DiffPair_0 _S3086;
    (&_S3086)->primal_0 = _S2935;
    (&_S3086)->differential_0 = _S3035;
    DiffPair_float_0 _S3087;
    (&_S3087)->primal_0 = 0.0f;
    (&_S3087)->differential_0 = _S3079;
    s_bwd_prop_d_atan2_0(&_S3085, &_S3086, &_S3087);
    float3  _S3088 = make_float3 (_S3083.differential_0.primal_0.x + _S3085.differential_0.primal_0, _S3080.differential_0.primal_0, _S3083.differential_0.primal_0.y + _S3086.differential_0.primal_0);
    float2  _S3089 = make_float2 (0.0f, _S3063.rows[int(1)].x) + make_float2 (_S3063.rows[int(0)].x, 0.0f);
    float _S3090 = fy_15 * _S3089.y;
    float _S3091 = fx_15 * _S3089.x;
    DiffPair_0 _S3092;
    (&_S3092)->primal_0 = _S2945;
    (&_S3092)->differential_0 = _S3035;
    DiffPair_0 _S3093;
    (&_S3093)->primal_0 = _S2946;
    (&_S3093)->differential_0 = _S3035;
    DiffPair_float_0 _S3094;
    (&_S3094)->primal_0 = 0.0f;
    (&_S3094)->differential_0 = _S3090;
    s_bwd_prop_d_atan2_0(&_S3092, &_S3093, &_S3094);
    DiffPair_1 _S3095;
    (&_S3095)->primal_0 = _S2942;
    (&_S3095)->differential_0 = _S3070;
    DiffPair_float_0 _S3096;
    (&_S3096)->primal_0 = _S3093.differential_0.primal_0;
    (&_S3096)->differential_0 = _S3093.differential_0.differential_0;
    s_bwd_prop_s_fwd_length_impl_0(&_S3095, &_S3096);
    DiffPair_0 _S3097;
    (&_S3097)->primal_0 = _S2934;
    (&_S3097)->differential_0 = _S3035;
    DiffPair_0 _S3098;
    (&_S3098)->primal_0 = _S2935;
    (&_S3098)->differential_0 = _S3035;
    DiffPair_float_0 _S3099;
    (&_S3099)->primal_0 = 0.0f;
    (&_S3099)->differential_0 = _S3091;
    s_bwd_prop_d_atan2_0(&_S3097, &_S3098, &_S3099);
    float3  _S3100 = make_float3 (_S3095.differential_0.primal_0.x + _S3097.differential_0.primal_0, _S3092.differential_0.primal_0, _S3095.differential_0.primal_0.y + _S3098.differential_0.primal_0);
    float _S3101 = length_0(_S2940);
    float _S3102 = fy_15 * v_mean2d_3.y;
    float _S3103 = fx_15 * v_mean2d_3.x;
    DiffPair_float_0 _S3104;
    (&_S3104)->primal_0 = _S2939;
    (&_S3104)->differential_0 = 0.0f;
    DiffPair_float_0 _S3105;
    (&_S3105)->primal_0 = _S3101;
    (&_S3105)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3104, &_S3105, _S3102);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3106;
    (&_S3106)->primal_0 = _S2940;
    (&_S3106)->differential_0 = _S3008;
    s_bwd_length_impl_2(&_S3106, _S3105.differential_0);
    DiffPair_float_0 _S3107;
    (&_S3107)->primal_0 = _S2932;
    (&_S3107)->differential_0 = 0.0f;
    DiffPair_float_0 _S3108;
    (&_S3108)->primal_0 = _S2933;
    (&_S3108)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3107, &_S3108, _S3103);
    float3  _S3109 = make_float3 (_S3106.differential_0.x + _S3107.differential_0, _S3104.differential_0, _S3106.differential_0.y + _S3108.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3110;
    (&_S3110)->primal_0 = _S2929;
    (&_S3110)->differential_0 = _S3061;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3111;
    (&_S3111)->primal_0 = _S2930;
    (&_S3111)->differential_0 = _S3061;
    s_bwd_prop_mul_2(&_S3110, &_S3111, _S3062.differential_0);
    Matrix<float, 3, 3>  _S3112 = transpose_3(_S3111.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3113;
    (&_S3113)->primal_0 = R_11;
    (&_S3113)->differential_0 = _S3061;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3114;
    (&_S3114)->primal_0 = _S2928;
    (&_S3114)->differential_0 = _S3061;
    s_bwd_prop_mul_2(&_S3113, &_S3114, _S3110.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3115;
    (&_S3115)->primal_0 = _S2926;
    (&_S3115)->differential_0 = _S3061;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3116;
    (&_S3116)->primal_0 = _S2927;
    (&_S3116)->differential_0 = _S3061;
    s_bwd_prop_mul_2(&_S3115, &_S3116, _S3114.differential_0);
    Matrix<float, 3, 3>  _S3117 = _S3115.differential_0 + transpose_3(_S3116.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3118;
    (&_S3118)->primal_0 = _S2925;
    (&_S3118)->differential_0 = _S3061;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3119;
    (&_S3119)->primal_0 = S_3;
    (&_S3119)->differential_0 = _S3061;
    s_bwd_prop_mul_2(&_S3118, &_S3119, _S3117);
    Matrix<float, 3, 3>  _S3120 = transpose_3(_S3118.differential_0);
    float _S3121 = 2.0f * - _S3120.rows[int(2)].z;
    float _S3122 = 2.0f * _S3120.rows[int(2)].y;
    float _S3123 = 2.0f * _S3120.rows[int(2)].x;
    float _S3124 = 2.0f * _S3120.rows[int(1)].z;
    float _S3125 = 2.0f * - _S3120.rows[int(1)].y;
    float _S3126 = 2.0f * _S3120.rows[int(1)].x;
    float _S3127 = 2.0f * _S3120.rows[int(0)].z;
    float _S3128 = 2.0f * _S3120.rows[int(0)].y;
    float _S3129 = 2.0f * - _S3120.rows[int(0)].x;
    float _S3130 = - _S3126 + _S3128;
    float _S3131 = _S3123 + - _S3127;
    float _S3132 = - _S3122 + _S3124;
    float _S3133 = _S3122 + _S3124;
    float _S3134 = _S3123 + _S3127;
    float _S3135 = _S3126 + _S3128;
    float _S3136 = _S2922.w * (_S3125 + _S3129);
    float _S3137 = _S2922.z * (_S3121 + _S3129);
    float _S3138 = _S2922.y * (_S3121 + _S3125);
    float _S3139 = _S2922.x * _S3130 + _S2922.z * _S3133 + _S2922.y * _S3134 + _S3136 + _S3136;
    float _S3140 = _S2922.x * _S3131 + _S2922.w * _S3133 + _S2922.y * _S3135 + _S3137 + _S3137;
    float _S3141 = _S2922.x * _S3132 + _S2922.w * _S3134 + _S2922.z * _S3135 + _S3138 + _S3138;
    float _S3142 = _S2922.w * _S3130 + _S2922.z * _S3131 + _S2922.y * _S3132;
    float3  _S3143 = _S3000;
    *&((&_S3143)->z) = _S3119.differential_0.rows[int(2)].z;
    *&((&_S3143)->y) = _S3119.differential_0.rows[int(1)].y;
    *&((&_S3143)->x) = _S3119.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3144;
    (&_S3144)->primal_0 = scale_11;
    (&_S3144)->differential_0 = _S3000;
    s_bwd_prop_exp_1(&_S3144, _S3143);
    float4  _S3145 = make_float4 (0.0f);
    float4  _S3146 = _S3145;
    *&((&_S3146)->w) = _S3139;
    *&((&_S3146)->z) = _S3140;
    *&((&_S3146)->y) = _S3141;
    *&((&_S3146)->x) = _S3142;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3147;
    (&_S3147)->primal_0 = quat_11;
    (&_S3147)->differential_0 = _S3145;
    s_bwd_normalize_impl_0(&_S3147, _S3146);
    float _S3148 = - (_S3027 / _S2921);
    DiffPair_float_0 _S3149;
    (&_S3149)->primal_0 = _S2918;
    (&_S3149)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3149, _S3148);
    float _S3150 = - _S3149.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3151;
    (&_S3151)->primal_0 = mean_c_11;
    (&_S3151)->differential_0 = _S3000;
    s_bwd_length_impl_0(&_S3151, v_depth_3);
    float3  _S3152 = _S3076 + _S3088 + _S3100 + _S3109 + _S3151.differential_0 + _S3025;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3153;
    (&_S3153)->primal_0 = R_11;
    (&_S3153)->differential_0 = _S3061;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3154;
    (&_S3154)->primal_0 = mean_12;
    (&_S3154)->differential_0 = _S3000;
    s_bwd_prop_mul_3(&_S3153, &_S3154, _S3152);
    Matrix<float, 3, 3>  _S3155 = _S3112 + _S3113.differential_0 + _S3153.differential_0;
    float _S3156 = _S3150 + _S3007.differential_0;
    float3  _S3157 = _S3144.differential_0 + _S3006.differential_0;
    *v_mean_3 = *v_mean_3 + (_S3154.differential_0 + _S3005.differential_0);
    *v_quat_3 = *v_quat_3 + _S3147.differential_0;
    *v_scale_3 = *v_scale_3 + _S3157;
    *v_in_opacity_3 = *v_in_opacity_3 + _S3156;
    *v_R_3 = *v_R_3 + _S3155;
    *v_t_3 = *v_t_3 + _S3152;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_persp_projection_Intermediates_0
{
    float2  _S3158;
    float2  _S3159;
    float2  _S3160;
    float2  _S3161;
    float2  _S3162;
    float2  _S3163;
    float2  _S3164;
};

inline __device__ void projection_3dgut_persp_vjp(bool antialiased_12, float3  mean_13, float4  quat_12, float3  scale_12, float in_opacity_12, Matrix<float, 3, 3>  R_12, float3  t_12, float fx_16, float fy_16, float cx_13, float cy_13, FixedArray<float, 10>  dist_coeffs_16, uint image_width_12, uint image_height_12, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float2  _S3165 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_persp_projection_Intermediates_0 _S3166;
    (&_S3166)->_S3158 = _S3165;
    (&_S3166)->_S3159 = _S3165;
    (&_S3166)->_S3160 = _S3165;
    (&_S3166)->_S3161 = _S3165;
    (&_S3166)->_S3162 = _S3165;
    (&_S3166)->_S3163 = _S3165;
    (&_S3166)->_S3164 = _S3165;
    float3  _S3167 = make_float3 (0.0f);
    float3  _S3168 = s_primal_ctx_exp_1(scale_12);
    float4  _S3169 = normalize_0(quat_12);
    float _S3170 = _S3169.y;
    float x2_12 = _S3170 * _S3170;
    float y2_12 = _S3169.z * _S3169.z;
    float z2_12 = _S3169.w * _S3169.w;
    float xy_12 = _S3169.y * _S3169.z;
    float xz_12 = _S3169.y * _S3169.w;
    float yz_12 = _S3169.z * _S3169.w;
    float wx_12 = _S3169.x * _S3169.y;
    float wy_12 = _S3169.x * _S3169.z;
    float wz_12 = _S3169.x * _S3169.w;
    Matrix<float, 3, 3>  _S3171 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_12), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_12), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))));
    FixedArray<float3 , 7>  _S3172 = {
        _S3167, _S3167, _S3167, _S3167, _S3167, _S3167, _S3167
    };
    FixedArray<float, 7>  _S3173 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3174;
    (&_S3174)->p_0 = _S3172;
    (&_S3174)->w_mean_0 = _S3173;
    (&_S3174)->w_cov_0 = _S3173;
    (&_S3174)->p_0[int(0)] = mean_13;
    SigmaPoints_0 _S3175 = _S3174;
    (&_S3175)->w_mean_0[int(0)] = 0.0f;
    (&_S3175)->w_cov_0[int(0)] = 2.0f;
    float _S3176 = s_primal_ctx_sqrt_0(3.0f);
    float _S3177 = _S3176 * _S3168.x;
    float3  delta_12 = make_float3 (_S3177) * _S3171.rows[0U];
    float3  _S3178 = mean_13 + delta_12;
    (&_S3175)->p_0[1U] = _S3178;
    float3  _S3179 = mean_13 - delta_12;
    (&_S3175)->p_0[4U] = _S3179;
    float _S3180 = _S3176 * _S3168.y;
    float3  delta_13 = make_float3 (_S3180) * _S3171.rows[1U];
    float3  _S3181 = mean_13 + delta_13;
    (&_S3175)->p_0[2U] = _S3181;
    float3  _S3182 = mean_13 - delta_13;
    (&_S3175)->p_0[5U] = _S3182;
    float _S3183 = _S3176 * _S3168.z;
    float3  delta_14 = make_float3 (_S3183) * _S3171.rows[2U];
    float3  _S3184 = mean_13 + delta_14;
    (&_S3175)->p_0[3U] = _S3184;
    float3  _S3185 = mean_13 - delta_14;
    (&_S3175)->p_0[6U] = _S3185;
    (&_S3175)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3186 = _S3175;
    (&_S3186)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3187 = _S3186;
    (&_S3187)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3188 = _S3187;
    (&_S3188)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3189 = _S3188;
    (&_S3189)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3190 = _S3189;
    (&_S3190)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3191 = _S3190;
    (&_S3191)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3192 = _S3191;
    (&_S3192)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3193 = _S3192;
    (&_S3193)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3194 = _S3193;
    (&_S3194)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3195 = _S3194;
    (&_S3195)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3196 = _S3195;
    (&_S3196)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3197 = _S3174;
    float3  _S3198 = s_primal_ctx_mul_0(R_12, _S3174.p_0[0U]) + t_12;
    _S3174 = _S3196;
    (&_S3174)->p_0[0U] = _S3198;
    SigmaPoints_0 _S3199 = _S3174;
    (&_S3174)->p_0[1U] = s_primal_ctx_mul_0(R_12, _S3178) + t_12;
    SigmaPoints_0 _S3200 = _S3174;
    (&_S3174)->p_0[2U] = s_primal_ctx_mul_0(R_12, _S3181) + t_12;
    SigmaPoints_0 _S3201 = _S3174;
    (&_S3174)->p_0[3U] = s_primal_ctx_mul_0(R_12, _S3184) + t_12;
    SigmaPoints_0 _S3202 = _S3174;
    (&_S3174)->p_0[4U] = s_primal_ctx_mul_0(R_12, _S3179) + t_12;
    SigmaPoints_0 _S3203 = _S3174;
    (&_S3174)->p_0[5U] = s_primal_ctx_mul_0(R_12, _S3182) + t_12;
    SigmaPoints_0 _S3204 = _S3174;
    (&_S3174)->p_0[6U] = s_primal_ctx_mul_0(R_12, _S3185) + t_12;
    float2  _S3205 = float2 {_S3199.p_0[int(0)].x, _S3199.p_0[int(0)].y} / make_float2 (_S3199.p_0[int(0)].z);
    float u_84 = _S3205.x;
    float v_84 = _S3205.y;
    float r2_84 = u_84 * u_84 + v_84 * v_84;
    float _S3206 = 2.0f * dist_coeffs_16[int(4)];
    float _S3207 = 2.0f * dist_coeffs_16[int(5)];
    float2  _S3208 = _S3205 * make_float2 (1.0f + r2_84 * (dist_coeffs_16[int(0)] + r2_84 * (dist_coeffs_16[int(1)] + r2_84 * (dist_coeffs_16[int(2)] + r2_84 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_84 * v_84 + dist_coeffs_16[int(5)] * (r2_84 + 2.0f * u_84 * u_84) + dist_coeffs_16[int(6)] * r2_84, _S3207 * u_84 * v_84 + dist_coeffs_16[int(4)] * (r2_84 + 2.0f * v_84 * v_84) + dist_coeffs_16[int(7)] * r2_84);
    float2  _S3209 = _S3208 + make_float2 (dist_coeffs_16[int(8)] * _S3208.x + dist_coeffs_16[int(9)] * _S3208.y, 0.0f);
    (&_S3166)->_S3158 = make_float2 (fx_16 * _S3209.x + cx_13, fy_16 * _S3209.y + cy_13);
    float2  _S3210 = float2 {_S3200.p_0[int(1)].x, _S3200.p_0[int(1)].y} / make_float2 (_S3200.p_0[int(1)].z);
    float u_85 = _S3210.x;
    float v_85 = _S3210.y;
    float r2_85 = u_85 * u_85 + v_85 * v_85;
    float2  _S3211 = _S3210 * make_float2 (1.0f + r2_85 * (dist_coeffs_16[int(0)] + r2_85 * (dist_coeffs_16[int(1)] + r2_85 * (dist_coeffs_16[int(2)] + r2_85 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_85 * v_85 + dist_coeffs_16[int(5)] * (r2_85 + 2.0f * u_85 * u_85) + dist_coeffs_16[int(6)] * r2_85, _S3207 * u_85 * v_85 + dist_coeffs_16[int(4)] * (r2_85 + 2.0f * v_85 * v_85) + dist_coeffs_16[int(7)] * r2_85);
    float2  _S3212 = _S3211 + make_float2 (dist_coeffs_16[int(8)] * _S3211.x + dist_coeffs_16[int(9)] * _S3211.y, 0.0f);
    (&_S3166)->_S3159 = make_float2 (fx_16 * _S3212.x + cx_13, fy_16 * _S3212.y + cy_13);
    float2  _S3213 = float2 {_S3201.p_0[int(2)].x, _S3201.p_0[int(2)].y} / make_float2 (_S3201.p_0[int(2)].z);
    float u_86 = _S3213.x;
    float v_86 = _S3213.y;
    float r2_86 = u_86 * u_86 + v_86 * v_86;
    float2  _S3214 = _S3213 * make_float2 (1.0f + r2_86 * (dist_coeffs_16[int(0)] + r2_86 * (dist_coeffs_16[int(1)] + r2_86 * (dist_coeffs_16[int(2)] + r2_86 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_86 * v_86 + dist_coeffs_16[int(5)] * (r2_86 + 2.0f * u_86 * u_86) + dist_coeffs_16[int(6)] * r2_86, _S3207 * u_86 * v_86 + dist_coeffs_16[int(4)] * (r2_86 + 2.0f * v_86 * v_86) + dist_coeffs_16[int(7)] * r2_86);
    float2  _S3215 = _S3214 + make_float2 (dist_coeffs_16[int(8)] * _S3214.x + dist_coeffs_16[int(9)] * _S3214.y, 0.0f);
    (&_S3166)->_S3160 = make_float2 (fx_16 * _S3215.x + cx_13, fy_16 * _S3215.y + cy_13);
    float2  _S3216 = float2 {_S3202.p_0[int(3)].x, _S3202.p_0[int(3)].y} / make_float2 (_S3202.p_0[int(3)].z);
    float u_87 = _S3216.x;
    float v_87 = _S3216.y;
    float r2_87 = u_87 * u_87 + v_87 * v_87;
    float2  _S3217 = _S3216 * make_float2 (1.0f + r2_87 * (dist_coeffs_16[int(0)] + r2_87 * (dist_coeffs_16[int(1)] + r2_87 * (dist_coeffs_16[int(2)] + r2_87 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_87 * v_87 + dist_coeffs_16[int(5)] * (r2_87 + 2.0f * u_87 * u_87) + dist_coeffs_16[int(6)] * r2_87, _S3207 * u_87 * v_87 + dist_coeffs_16[int(4)] * (r2_87 + 2.0f * v_87 * v_87) + dist_coeffs_16[int(7)] * r2_87);
    float2  _S3218 = _S3217 + make_float2 (dist_coeffs_16[int(8)] * _S3217.x + dist_coeffs_16[int(9)] * _S3217.y, 0.0f);
    (&_S3166)->_S3161 = make_float2 (fx_16 * _S3218.x + cx_13, fy_16 * _S3218.y + cy_13);
    float2  _S3219 = float2 {_S3203.p_0[int(4)].x, _S3203.p_0[int(4)].y} / make_float2 (_S3203.p_0[int(4)].z);
    float u_88 = _S3219.x;
    float v_88 = _S3219.y;
    float r2_88 = u_88 * u_88 + v_88 * v_88;
    float2  _S3220 = _S3219 * make_float2 (1.0f + r2_88 * (dist_coeffs_16[int(0)] + r2_88 * (dist_coeffs_16[int(1)] + r2_88 * (dist_coeffs_16[int(2)] + r2_88 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_88 * v_88 + dist_coeffs_16[int(5)] * (r2_88 + 2.0f * u_88 * u_88) + dist_coeffs_16[int(6)] * r2_88, _S3207 * u_88 * v_88 + dist_coeffs_16[int(4)] * (r2_88 + 2.0f * v_88 * v_88) + dist_coeffs_16[int(7)] * r2_88);
    float2  _S3221 = _S3220 + make_float2 (dist_coeffs_16[int(8)] * _S3220.x + dist_coeffs_16[int(9)] * _S3220.y, 0.0f);
    (&_S3166)->_S3162 = make_float2 (fx_16 * _S3221.x + cx_13, fy_16 * _S3221.y + cy_13);
    float2  _S3222 = float2 {_S3204.p_0[int(5)].x, _S3204.p_0[int(5)].y} / make_float2 (_S3204.p_0[int(5)].z);
    float u_89 = _S3222.x;
    float v_89 = _S3222.y;
    float r2_89 = u_89 * u_89 + v_89 * v_89;
    float2  _S3223 = _S3222 * make_float2 (1.0f + r2_89 * (dist_coeffs_16[int(0)] + r2_89 * (dist_coeffs_16[int(1)] + r2_89 * (dist_coeffs_16[int(2)] + r2_89 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_89 * v_89 + dist_coeffs_16[int(5)] * (r2_89 + 2.0f * u_89 * u_89) + dist_coeffs_16[int(6)] * r2_89, _S3207 * u_89 * v_89 + dist_coeffs_16[int(4)] * (r2_89 + 2.0f * v_89 * v_89) + dist_coeffs_16[int(7)] * r2_89);
    float2  _S3224 = _S3223 + make_float2 (dist_coeffs_16[int(8)] * _S3223.x + dist_coeffs_16[int(9)] * _S3223.y, 0.0f);
    (&_S3166)->_S3163 = make_float2 (fx_16 * _S3224.x + cx_13, fy_16 * _S3224.y + cy_13);
    float2  _S3225 = float2 {_S3174.p_0[int(6)].x, _S3174.p_0[int(6)].y} / make_float2 (_S3174.p_0[int(6)].z);
    float u_90 = _S3225.x;
    float v_90 = _S3225.y;
    float r2_90 = u_90 * u_90 + v_90 * v_90;
    float2  _S3226 = _S3225 * make_float2 (1.0f + r2_90 * (dist_coeffs_16[int(0)] + r2_90 * (dist_coeffs_16[int(1)] + r2_90 * (dist_coeffs_16[int(2)] + r2_90 * dist_coeffs_16[int(3)])))) + make_float2 (_S3206 * u_90 * v_90 + dist_coeffs_16[int(5)] * (r2_90 + 2.0f * u_90 * u_90) + dist_coeffs_16[int(6)] * r2_90, _S3207 * u_90 * v_90 + dist_coeffs_16[int(4)] * (r2_90 + 2.0f * v_90 * v_90) + dist_coeffs_16[int(7)] * r2_90);
    float2  _S3227 = _S3226 + make_float2 (dist_coeffs_16[int(8)] * _S3226.x + dist_coeffs_16[int(9)] * _S3226.y, 0.0f);
    (&_S3166)->_S3164 = make_float2 (fx_16 * _S3227.x + cx_13, fy_16 * _S3227.y + cy_13);
    float3  mean_c_12 = s_primal_ctx_mul_0(R_12, mean_13) + t_12;
    float _S3228 = - in_opacity_12;
    float _S3229 = 1.0f + s_primal_ctx_exp_0(_S3228);
    float _S3230 = 1.0f / _S3229;
    float _S3231 = _S3229 * _S3229;
    float3  _S3232 = make_float3 (_S3177);
    float3  _S3233 = make_float3 (_S3180);
    float3  _S3234 = make_float3 (_S3183);
    float _S3235 = float(image_width_12);
    float _S3236 = float(image_height_12);
    float _S3237 = 0.30000001192092896f * (0.5f * _S3235 / fx_16) * fx_16;
    float lim_x_pos_4 = _S3235 + _S3237;
    float _S3238 = 0.30000001192092896f * (0.5f * _S3236 / fy_16) * fy_16;
    float lim_y_pos_1 = _S3236 + _S3238;
    float2  _S3239 = make_float2 (_S3175.w_mean_0[int(1)]) * _S3166._S3159 + make_float2 (_S3187.w_mean_0[int(2)]) * _S3166._S3160 + make_float2 (_S3189.w_mean_0[int(3)]) * _S3166._S3161 + make_float2 (_S3191.w_mean_0[int(4)]) * _S3166._S3162 + make_float2 (_S3193.w_mean_0[int(5)]) * _S3166._S3163 + make_float2 (_S3195.w_mean_0[int(6)]) * _S3166._S3164;
    float _S3240 = - _S3237;
    float _S3241 = - _S3238;
    float2  _S3242 = make_float2 (s_primal_ctx_clamp_0(_S3239.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3239.y, _S3241, lim_y_pos_1));
    float2  d_28 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3158.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3158.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3243 = d_28.x;
    float _S3244 = d_28.y;
    float _S3245 = _S3243 * _S3244;
    float2  d_29 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3159.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3159.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3246 = d_29.x;
    float _S3247 = d_29.y;
    float _S3248 = _S3246 * _S3247;
    float2  d_30 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3160.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3160.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3249 = d_30.x;
    float _S3250 = d_30.y;
    float _S3251 = _S3249 * _S3250;
    float2  d_31 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3161.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3161.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3252 = d_31.x;
    float _S3253 = d_31.y;
    float _S3254 = _S3252 * _S3253;
    float2  d_32 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3162.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3162.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3255 = d_32.x;
    float _S3256 = d_32.y;
    float _S3257 = _S3255 * _S3256;
    float2  d_33 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3163.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3163.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3258 = d_33.x;
    float _S3259 = d_33.y;
    float _S3260 = _S3258 * _S3259;
    float2  d_34 = make_float2 (s_primal_ctx_clamp_0(_S3166._S3164.x, _S3240, lim_x_pos_4), s_primal_ctx_clamp_0(_S3166._S3164.y, _S3241, lim_y_pos_1)) - _S3242;
    float _S3261 = d_34.x;
    float _S3262 = d_34.y;
    float _S3263 = _S3261 * _S3262;
    Matrix<float, 2, 2>  covar2d_8 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3243 * _S3243, _S3245, _S3245, _S3244 * _S3244) + makeMatrix<float, 2, 2> (_S3186.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3246 * _S3246, _S3248, _S3248, _S3247 * _S3247) + makeMatrix<float, 2, 2> (_S3188.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3249 * _S3249, _S3251, _S3251, _S3250 * _S3250) + makeMatrix<float, 2, 2> (_S3190.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3252 * _S3252, _S3254, _S3254, _S3253 * _S3253) + makeMatrix<float, 2, 2> (_S3192.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3255 * _S3255, _S3257, _S3257, _S3256 * _S3256) + makeMatrix<float, 2, 2> (_S3194.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3258 * _S3258, _S3260, _S3260, _S3259 * _S3259) + makeMatrix<float, 2, 2> (_S3196.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3261 * _S3261, _S3263, _S3263, _S3262 * _S3262);
    float eps2d_12;
    if(antialiased_12)
    {
        eps2d_12 = 0.10000000149011612f;
    }
    else
    {
        eps2d_12 = 0.30000001192092896f;
    }
    float _S3264 = covar2d_8.rows[int(0)].y * covar2d_8.rows[int(1)].x;
    float det_orig_12 = covar2d_8.rows[int(0)].x * covar2d_8.rows[int(1)].y - _S3264;
    float _S3265 = covar2d_8.rows[int(0)].x + eps2d_12;
    Matrix<float, 2, 2>  _S3266 = covar2d_8;
    *&(((&_S3266)->rows + (int(0)))->x) = _S3265;
    float _S3267 = covar2d_8.rows[int(1)].y + eps2d_12;
    *&(((&_S3266)->rows + (int(1)))->y) = _S3267;
    Matrix<float, 2, 2>  _S3268 = _S3266;
    Matrix<float, 2, 2>  _S3269 = _S3266;
    float det_blur_12 = _S3265 * _S3267 - _S3264;
    float _S3270 = det_orig_12 / det_blur_12;
    float _S3271 = det_blur_12 * det_blur_12;
    float _S3272 = (F32_max((0.0f), (_S3270)));
    float _S3273 = s_primal_ctx_sqrt_0(_S3272);
    float invdet_14 = 1.0f / det_blur_12;
    float _S3274 = - covar2d_8.rows[int(0)].y;
    float _S3275 = - covar2d_8.rows[int(1)].x;
    if(antialiased_12)
    {
        eps2d_12 = _S3230 * _S3273;
    }
    else
    {
        eps2d_12 = _S3230;
    }
    float _S3276 = eps2d_12 / 0.00392156885936856f;
    float _S3277 = 2.0f * s_primal_ctx_log_0(_S3276);
    float _S3278 = s_primal_ctx_sqrt_0(_S3277);
    float _S3279 = _S3268.rows[int(0)].x;
    float _S3280 = _S3269.rows[int(1)].y;
    float3  campos_5 = - s_primal_ctx_mul_0(transpose_3(R_12), t_12);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3281;
    (&_S3281)->primal_0 = mean_13;
    (&_S3281)->differential_0 = _S3167;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3282;
    (&_S3282)->primal_0 = scale_12;
    (&_S3282)->differential_0 = _S3167;
    DiffPair_float_0 _S3283;
    (&_S3283)->primal_0 = in_opacity_12;
    (&_S3283)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3284;
    (&_S3284)->primal_0 = campos_5;
    (&_S3284)->differential_0 = _S3167;
    s_bwd_prop_view_radius_3dgs_0(&_S3281, &_S3282, &_S3283, &_S3284, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3285 = _S3281;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3286 = _S3282;
    DiffPair_float_0 _S3287 = _S3283;
    float2  _S3288 = _S3165;
    *&((&_S3288)->y) = v_conic_4.z;
    float2  _S3289 = _S3165;
    *&((&_S3289)->y) = v_conic_4.y;
    *&((&_S3289)->x) = v_conic_4.x;
    DiffPair_float_0 _S3290;
    (&_S3290)->primal_0 = _S3280;
    (&_S3290)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3290, 0.0f);
    DiffPair_float_0 _S3291;
    (&_S3291)->primal_0 = _S3279;
    (&_S3291)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3291, 0.0f);
    DiffPair_float_0 _S3292;
    (&_S3292)->primal_0 = 3.32999992370605469f;
    (&_S3292)->differential_0 = 0.0f;
    DiffPair_float_0 _S3293;
    (&_S3293)->primal_0 = _S3278;
    (&_S3293)->differential_0 = 0.0f;
    _d_min_0(&_S3292, &_S3293, 0.0f);
    DiffPair_float_0 _S3294;
    (&_S3294)->primal_0 = _S3277;
    (&_S3294)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3294, _S3293.differential_0);
    float _S3295 = 2.0f * _S3294.differential_0;
    DiffPair_float_0 _S3296;
    (&_S3296)->primal_0 = _S3276;
    (&_S3296)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3296, _S3295);
    float _S3297 = v_opacity_4 + 254.9999847412109375f * _S3296.differential_0;
    Matrix<float, 2, 2>  _S3298 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3299 = _S3298;
    _S3299[int(1)] = _S3288;
    _S3299[int(0)] = _S3289;
    Matrix<float, 2, 2>  _S3300 = _S3299;
    float2  _S3301 = make_float2 (0.0f, _S3290.differential_0);
    float2  _S3302 = make_float2 (_S3291.differential_0, 0.0f);
    float _S3303;
    if(antialiased_12)
    {
        float _S3304 = _S3273 * _S3297;
        eps2d_12 = _S3230 * _S3297;
        _S3303 = _S3304;
    }
    else
    {
        eps2d_12 = 0.0f;
        _S3303 = _S3297;
    }
    float _S3305 = invdet_14 * _S3300.rows[int(1)].y;
    float _S3306 = - (invdet_14 * _S3300.rows[int(1)].x);
    float _S3307 = - (invdet_14 * _S3300.rows[int(0)].y);
    float _S3308 = invdet_14 * _S3300.rows[int(0)].x;
    float _S3309 = - ((_S3265 * _S3300.rows[int(1)].y + _S3275 * _S3300.rows[int(1)].x + _S3274 * _S3300.rows[int(0)].y + _S3267 * _S3300.rows[int(0)].x) / _S3271);
    DiffPair_float_0 _S3310;
    (&_S3310)->primal_0 = _S3272;
    (&_S3310)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3310, eps2d_12);
    DiffPair_float_0 _S3311;
    (&_S3311)->primal_0 = 0.0f;
    (&_S3311)->differential_0 = 0.0f;
    DiffPair_float_0 _S3312;
    (&_S3312)->primal_0 = _S3270;
    (&_S3312)->differential_0 = 0.0f;
    _d_max_0(&_S3311, &_S3312, _S3310.differential_0);
    float _S3313 = _S3312.differential_0 / _S3271;
    float s_diff_det_orig_T_4 = det_blur_12 * _S3313;
    float _S3314 = det_orig_12 * - _S3313 + _S3309;
    float _S3315 = - _S3314;
    float _S3316 = _S3265 * _S3314;
    float _S3317 = _S3267 * _S3314;
    Matrix<float, 2, 2>  _S3318 = _S3298;
    _S3318[int(1)] = _S3301;
    _S3318[int(0)] = _S3302;
    float _S3319 = _S3317 + _S3318.rows[int(0)].x + _S3305;
    float _S3320 = _S3315 + - s_diff_det_orig_T_4;
    float _S3321 = covar2d_8.rows[int(0)].y * _S3320 + _S3306;
    float _S3322 = covar2d_8.rows[int(1)].x * _S3320 + _S3307;
    float _S3323 = covar2d_8.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S3324 = _S3316 + _S3318.rows[int(1)].y + _S3308 + covar2d_8.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S3325 = _S3165;
    *&((&_S3325)->x) = _S3321;
    *&((&_S3325)->y) = _S3324;
    float _S3326 = _S3319 + _S3323;
    float2  _S3327 = _S3165;
    *&((&_S3327)->y) = _S3322;
    *&((&_S3327)->x) = _S3326;
    Matrix<float, 2, 2>  _S3328 = _S3298;
    _S3328[int(1)] = _S3325;
    _S3328[int(0)] = _S3327;
    Matrix<float, 3, 3>  _S3329 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3330;
    (&_S3330)->primal_0 = R_12;
    (&_S3330)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3331;
    (&_S3331)->primal_0 = _S3185;
    (&_S3331)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3330, &_S3331, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3332;
    (&_S3332)->primal_0 = R_12;
    (&_S3332)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3333;
    (&_S3333)->primal_0 = _S3182;
    (&_S3333)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3332, &_S3333, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3334;
    (&_S3334)->primal_0 = R_12;
    (&_S3334)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3335;
    (&_S3335)->primal_0 = _S3179;
    (&_S3335)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3334, &_S3335, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3336;
    (&_S3336)->primal_0 = R_12;
    (&_S3336)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3337;
    (&_S3337)->primal_0 = _S3184;
    (&_S3337)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3336, &_S3337, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3338;
    (&_S3338)->primal_0 = R_12;
    (&_S3338)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3339;
    (&_S3339)->primal_0 = _S3181;
    (&_S3339)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3338, &_S3339, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3340;
    (&_S3340)->primal_0 = R_12;
    (&_S3340)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3341;
    (&_S3341)->primal_0 = _S3178;
    (&_S3341)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3340, &_S3341, _S3167);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3342;
    (&_S3342)->primal_0 = R_12;
    (&_S3342)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3343;
    (&_S3343)->primal_0 = _S3197.p_0[0U];
    (&_S3343)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3342, &_S3343, _S3167);
    float3  _S3344 = - _S3331.differential_0 + _S3337.differential_0;
    float3  _S3345 = _S3234 * _S3344;
    float3  _S3346 = _S3171.rows[2U] * _S3344;
    float _S3347 = _S3176 * (_S3346.x + _S3346.y + _S3346.z);
    float3  _S3348 = - _S3333.differential_0 + _S3339.differential_0;
    float3  _S3349 = _S3233 * _S3348;
    float3  _S3350 = _S3171.rows[1U] * _S3348;
    float _S3351 = _S3176 * (_S3350.x + _S3350.y + _S3350.z);
    float3  _S3352 = - _S3335.differential_0 + _S3341.differential_0;
    float3  _S3353 = _S3232 * _S3352;
    float3  _S3354 = _S3171.rows[0U] * _S3352;
    float _S3355 = _S3176 * (_S3354.x + _S3354.y + _S3354.z);
    Matrix<float, 3, 3>  _S3356 = _S3329;
    _S3356[2U] = _S3345;
    _S3356[1U] = _S3349;
    _S3356[0U] = _S3353;
    Matrix<float, 3, 3>  _S3357 = transpose_3(transpose_3(_S3356));
    float _S3358 = 2.0f * - _S3357.rows[int(2)].z;
    float _S3359 = 2.0f * _S3357.rows[int(2)].y;
    float _S3360 = 2.0f * _S3357.rows[int(2)].x;
    float _S3361 = 2.0f * _S3357.rows[int(1)].z;
    float _S3362 = 2.0f * - _S3357.rows[int(1)].y;
    float _S3363 = 2.0f * _S3357.rows[int(1)].x;
    float _S3364 = 2.0f * _S3357.rows[int(0)].z;
    float _S3365 = 2.0f * _S3357.rows[int(0)].y;
    float _S3366 = 2.0f * - _S3357.rows[int(0)].x;
    float _S3367 = - _S3363 + _S3365;
    float _S3368 = _S3360 + - _S3364;
    float _S3369 = - _S3359 + _S3361;
    float _S3370 = _S3359 + _S3361;
    float _S3371 = _S3360 + _S3364;
    float _S3372 = _S3363 + _S3365;
    float _S3373 = _S3169.w * (_S3362 + _S3366);
    float _S3374 = _S3169.z * (_S3358 + _S3366);
    float _S3375 = _S3169.y * (_S3358 + _S3362);
    float _S3376 = _S3169.x * _S3367 + _S3169.z * _S3370 + _S3169.y * _S3371 + _S3373 + _S3373;
    float _S3377 = _S3169.x * _S3368 + _S3169.w * _S3370 + _S3169.y * _S3372 + _S3374 + _S3374;
    float _S3378 = _S3169.x * _S3369 + _S3169.w * _S3371 + _S3169.z * _S3372 + _S3375 + _S3375;
    float _S3379 = _S3169.w * _S3367 + _S3169.z * _S3368 + _S3169.y * _S3369;
    float4  _S3380 = make_float4 (0.0f);
    float4  _S3381 = _S3380;
    *&((&_S3381)->w) = _S3376;
    *&((&_S3381)->z) = _S3377;
    *&((&_S3381)->y) = _S3378;
    *&((&_S3381)->x) = _S3379;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3382;
    (&_S3382)->primal_0 = quat_12;
    (&_S3382)->differential_0 = _S3380;
    s_bwd_normalize_impl_0(&_S3382, _S3381);
    float3  _S3383 = _S3167;
    *&((&_S3383)->z) = _S3347;
    *&((&_S3383)->y) = _S3351;
    *&((&_S3383)->x) = _S3355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3384;
    (&_S3384)->primal_0 = scale_12;
    (&_S3384)->differential_0 = _S3167;
    s_bwd_prop_exp_1(&_S3384, _S3383);
    float _S3385 = - (_S3303 / _S3231);
    DiffPair_float_0 _S3386;
    (&_S3386)->primal_0 = _S3228;
    (&_S3386)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3386, _S3385);
    float _S3387 = - _S3386.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3388;
    (&_S3388)->primal_0 = mean_c_12;
    (&_S3388)->differential_0 = _S3167;
    s_bwd_length_impl_0(&_S3388, v_depth_4);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3389;
    (&_S3389)->primal_0 = R_12;
    (&_S3389)->differential_0 = _S3329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3390;
    (&_S3390)->primal_0 = mean_13;
    (&_S3390)->differential_0 = _S3167;
    s_bwd_prop_mul_3(&_S3389, &_S3390, _S3388.differential_0);
    Matrix<float, 3, 3>  _S3391 = _S3330.differential_0 + _S3332.differential_0 + _S3334.differential_0 + _S3336.differential_0 + _S3338.differential_0 + _S3340.differential_0 + _S3342.differential_0 + _S3389.differential_0;
    float _S3392 = _S3387 + _S3287.differential_0;
    float3  _S3393 = _S3384.differential_0 + _S3286.differential_0;
    *v_mean_4 = *v_mean_4 + (_S3331.differential_0 + _S3337.differential_0 + _S3333.differential_0 + _S3339.differential_0 + _S3335.differential_0 + _S3341.differential_0 + _S3390.differential_0 + _S3285.differential_0);
    *v_quat_4 = *v_quat_4 + _S3382.differential_0;
    *v_scale_4 = *v_scale_4 + _S3393;
    *v_in_opacity_4 = *v_in_opacity_4 + _S3392;
    *v_R_4 = *v_R_4 + _S3391;
    *v_t_4 = *v_t_4 + _S3388.differential_0;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_fisheye_projection_Intermediates_0
{
    float2  _S3394;
    float2  _S3395;
    float2  _S3396;
    float2  _S3397;
    float2  _S3398;
    float2  _S3399;
    float2  _S3400;
};

inline __device__ void projection_3dgut_fisheye_vjp(bool antialiased_13, float3  mean_14, float4  quat_13, float3  scale_13, float in_opacity_13, Matrix<float, 3, 3>  R_13, float3  t_13, float fx_17, float fy_17, float cx_14, float cy_14, FixedArray<float, 10>  dist_coeffs_17, uint image_width_13, uint image_height_13, float2  v_mean2d_5, float v_depth_5, float3  v_conic_5, float v_opacity_5, float3  * v_mean_5, float4  * v_quat_5, float3  * v_scale_5, float * v_in_opacity_5, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    float2  _S3401 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_fisheye_projection_Intermediates_0 _S3402;
    (&_S3402)->_S3394 = _S3401;
    (&_S3402)->_S3395 = _S3401;
    (&_S3402)->_S3396 = _S3401;
    (&_S3402)->_S3397 = _S3401;
    (&_S3402)->_S3398 = _S3401;
    (&_S3402)->_S3399 = _S3401;
    (&_S3402)->_S3400 = _S3401;
    (&_S3402)->_S3394 = _S3401;
    (&_S3402)->_S3395 = _S3401;
    (&_S3402)->_S3396 = _S3401;
    (&_S3402)->_S3397 = _S3401;
    (&_S3402)->_S3398 = _S3401;
    (&_S3402)->_S3399 = _S3401;
    (&_S3402)->_S3400 = _S3401;
    float3  _S3403 = make_float3 (0.0f);
    float3  _S3404 = s_primal_ctx_exp_1(scale_13);
    float4  _S3405 = normalize_0(quat_13);
    float _S3406 = _S3405.y;
    float x2_13 = _S3406 * _S3406;
    float y2_13 = _S3405.z * _S3405.z;
    float z2_13 = _S3405.w * _S3405.w;
    float xy_13 = _S3405.y * _S3405.z;
    float xz_13 = _S3405.y * _S3405.w;
    float yz_13 = _S3405.z * _S3405.w;
    float wx_13 = _S3405.x * _S3405.y;
    float wy_13 = _S3405.x * _S3405.z;
    float wz_13 = _S3405.x * _S3405.w;
    Matrix<float, 3, 3>  _S3407 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_13), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_13), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13))));
    FixedArray<float3 , 7>  _S3408 = {
        _S3403, _S3403, _S3403, _S3403, _S3403, _S3403, _S3403
    };
    FixedArray<float, 7>  _S3409 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3410;
    (&_S3410)->p_0 = _S3408;
    (&_S3410)->w_mean_0 = _S3409;
    (&_S3410)->w_cov_0 = _S3409;
    (&_S3410)->p_0[int(0)] = mean_14;
    SigmaPoints_0 _S3411 = _S3410;
    (&_S3411)->w_mean_0[int(0)] = 0.0f;
    (&_S3411)->w_cov_0[int(0)] = 2.0f;
    float _S3412 = s_primal_ctx_sqrt_0(3.0f);
    float _S3413 = _S3412 * _S3404.x;
    float3  delta_15 = make_float3 (_S3413) * _S3407.rows[0U];
    float3  _S3414 = mean_14 + delta_15;
    (&_S3411)->p_0[1U] = _S3414;
    float3  _S3415 = mean_14 - delta_15;
    (&_S3411)->p_0[4U] = _S3415;
    float _S3416 = _S3412 * _S3404.y;
    float3  delta_16 = make_float3 (_S3416) * _S3407.rows[1U];
    float3  _S3417 = mean_14 + delta_16;
    (&_S3411)->p_0[2U] = _S3417;
    float3  _S3418 = mean_14 - delta_16;
    (&_S3411)->p_0[5U] = _S3418;
    float _S3419 = _S3412 * _S3404.z;
    float3  delta_17 = make_float3 (_S3419) * _S3407.rows[2U];
    float3  _S3420 = mean_14 + delta_17;
    (&_S3411)->p_0[3U] = _S3420;
    float3  _S3421 = mean_14 - delta_17;
    (&_S3411)->p_0[6U] = _S3421;
    (&_S3411)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3422 = _S3411;
    (&_S3422)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3423 = _S3422;
    (&_S3423)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3424 = _S3423;
    (&_S3424)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3425 = _S3424;
    (&_S3425)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3426 = _S3425;
    (&_S3426)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3427 = _S3426;
    (&_S3427)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3428 = _S3427;
    (&_S3428)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3429 = _S3428;
    (&_S3429)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3430 = _S3429;
    (&_S3430)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3431 = _S3430;
    (&_S3431)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3432 = _S3431;
    (&_S3432)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3433 = _S3410;
    float3  _S3434 = s_primal_ctx_mul_0(R_13, _S3410.p_0[0U]) + t_13;
    _S3410 = _S3432;
    (&_S3410)->p_0[0U] = _S3434;
    SigmaPoints_0 _S3435 = _S3410;
    (&_S3410)->p_0[1U] = s_primal_ctx_mul_0(R_13, _S3414) + t_13;
    SigmaPoints_0 _S3436 = _S3410;
    (&_S3410)->p_0[2U] = s_primal_ctx_mul_0(R_13, _S3417) + t_13;
    SigmaPoints_0 _S3437 = _S3410;
    (&_S3410)->p_0[3U] = s_primal_ctx_mul_0(R_13, _S3420) + t_13;
    SigmaPoints_0 _S3438 = _S3410;
    (&_S3410)->p_0[4U] = s_primal_ctx_mul_0(R_13, _S3415) + t_13;
    SigmaPoints_0 _S3439 = _S3410;
    (&_S3410)->p_0[5U] = s_primal_ctx_mul_0(R_13, _S3418) + t_13;
    SigmaPoints_0 _S3440 = _S3410;
    (&_S3410)->p_0[6U] = s_primal_ctx_mul_0(R_13, _S3421) + t_13;
    SigmaPoints_0 _S3441 = _S3410;
    float2  _S3442 = float2 {_S3435.p_0[int(0)].x, _S3435.p_0[int(0)].y};
    float _S3443 = length_0(_S3442);
    float _S3444 = _S3435.p_0[int(0)].z;
    float _S3445 = s_primal_ctx_atan2_0(_S3443, _S3444);
    float k_8;
    if(_S3445 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3445 * _S3445 / 3.0f) / _S3444;
    }
    else
    {
        k_8 = _S3445 / _S3443;
    }
    float2  _S3446 = _S3442 * make_float2 (k_8);
    float u_91 = _S3446.x;
    float v_91 = _S3446.y;
    float r2_91 = u_91 * u_91 + v_91 * v_91;
    float _S3447 = 2.0f * dist_coeffs_17[int(4)];
    float _S3448 = 2.0f * dist_coeffs_17[int(5)];
    float2  _S3449 = _S3446 * make_float2 (1.0f + r2_91 * (dist_coeffs_17[int(0)] + r2_91 * (dist_coeffs_17[int(1)] + r2_91 * (dist_coeffs_17[int(2)] + r2_91 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_91 * v_91 + dist_coeffs_17[int(5)] * (r2_91 + 2.0f * u_91 * u_91) + dist_coeffs_17[int(6)] * r2_91, _S3448 * u_91 * v_91 + dist_coeffs_17[int(4)] * (r2_91 + 2.0f * v_91 * v_91) + dist_coeffs_17[int(7)] * r2_91);
    float2  _S3450 = _S3449 + make_float2 (dist_coeffs_17[int(8)] * _S3449.x + dist_coeffs_17[int(9)] * _S3449.y, 0.0f);
    (&_S3402)->_S3394 = make_float2 (fx_17 * _S3450.x + cx_14, fy_17 * _S3450.y + cy_14);
    float2  _S3451 = float2 {_S3436.p_0[int(1)].x, _S3436.p_0[int(1)].y};
    float _S3452 = length_0(_S3451);
    float _S3453 = _S3436.p_0[int(1)].z;
    float _S3454 = s_primal_ctx_atan2_0(_S3452, _S3453);
    if(_S3454 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3454 * _S3454 / 3.0f) / _S3453;
    }
    else
    {
        k_8 = _S3454 / _S3452;
    }
    float2  _S3455 = _S3451 * make_float2 (k_8);
    float u_92 = _S3455.x;
    float v_92 = _S3455.y;
    float r2_92 = u_92 * u_92 + v_92 * v_92;
    float2  _S3456 = _S3455 * make_float2 (1.0f + r2_92 * (dist_coeffs_17[int(0)] + r2_92 * (dist_coeffs_17[int(1)] + r2_92 * (dist_coeffs_17[int(2)] + r2_92 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_92 * v_92 + dist_coeffs_17[int(5)] * (r2_92 + 2.0f * u_92 * u_92) + dist_coeffs_17[int(6)] * r2_92, _S3448 * u_92 * v_92 + dist_coeffs_17[int(4)] * (r2_92 + 2.0f * v_92 * v_92) + dist_coeffs_17[int(7)] * r2_92);
    float2  _S3457 = _S3456 + make_float2 (dist_coeffs_17[int(8)] * _S3456.x + dist_coeffs_17[int(9)] * _S3456.y, 0.0f);
    (&_S3402)->_S3395 = make_float2 (fx_17 * _S3457.x + cx_14, fy_17 * _S3457.y + cy_14);
    float2  _S3458 = float2 {_S3437.p_0[int(2)].x, _S3437.p_0[int(2)].y};
    float _S3459 = length_0(_S3458);
    float _S3460 = _S3437.p_0[int(2)].z;
    float _S3461 = s_primal_ctx_atan2_0(_S3459, _S3460);
    if(_S3461 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3461 * _S3461 / 3.0f) / _S3460;
    }
    else
    {
        k_8 = _S3461 / _S3459;
    }
    float2  _S3462 = _S3458 * make_float2 (k_8);
    float u_93 = _S3462.x;
    float v_93 = _S3462.y;
    float r2_93 = u_93 * u_93 + v_93 * v_93;
    float2  _S3463 = _S3462 * make_float2 (1.0f + r2_93 * (dist_coeffs_17[int(0)] + r2_93 * (dist_coeffs_17[int(1)] + r2_93 * (dist_coeffs_17[int(2)] + r2_93 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_93 * v_93 + dist_coeffs_17[int(5)] * (r2_93 + 2.0f * u_93 * u_93) + dist_coeffs_17[int(6)] * r2_93, _S3448 * u_93 * v_93 + dist_coeffs_17[int(4)] * (r2_93 + 2.0f * v_93 * v_93) + dist_coeffs_17[int(7)] * r2_93);
    float2  _S3464 = _S3463 + make_float2 (dist_coeffs_17[int(8)] * _S3463.x + dist_coeffs_17[int(9)] * _S3463.y, 0.0f);
    (&_S3402)->_S3396 = make_float2 (fx_17 * _S3464.x + cx_14, fy_17 * _S3464.y + cy_14);
    float2  _S3465 = float2 {_S3438.p_0[int(3)].x, _S3438.p_0[int(3)].y};
    float _S3466 = length_0(_S3465);
    float _S3467 = _S3438.p_0[int(3)].z;
    float _S3468 = s_primal_ctx_atan2_0(_S3466, _S3467);
    if(_S3468 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3468 * _S3468 / 3.0f) / _S3467;
    }
    else
    {
        k_8 = _S3468 / _S3466;
    }
    float2  _S3469 = _S3465 * make_float2 (k_8);
    float u_94 = _S3469.x;
    float v_94 = _S3469.y;
    float r2_94 = u_94 * u_94 + v_94 * v_94;
    float2  _S3470 = _S3469 * make_float2 (1.0f + r2_94 * (dist_coeffs_17[int(0)] + r2_94 * (dist_coeffs_17[int(1)] + r2_94 * (dist_coeffs_17[int(2)] + r2_94 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_94 * v_94 + dist_coeffs_17[int(5)] * (r2_94 + 2.0f * u_94 * u_94) + dist_coeffs_17[int(6)] * r2_94, _S3448 * u_94 * v_94 + dist_coeffs_17[int(4)] * (r2_94 + 2.0f * v_94 * v_94) + dist_coeffs_17[int(7)] * r2_94);
    float2  _S3471 = _S3470 + make_float2 (dist_coeffs_17[int(8)] * _S3470.x + dist_coeffs_17[int(9)] * _S3470.y, 0.0f);
    (&_S3402)->_S3397 = make_float2 (fx_17 * _S3471.x + cx_14, fy_17 * _S3471.y + cy_14);
    float2  _S3472 = float2 {_S3439.p_0[int(4)].x, _S3439.p_0[int(4)].y};
    float _S3473 = length_0(_S3472);
    float _S3474 = _S3439.p_0[int(4)].z;
    float _S3475 = s_primal_ctx_atan2_0(_S3473, _S3474);
    if(_S3475 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3475 * _S3475 / 3.0f) / _S3474;
    }
    else
    {
        k_8 = _S3475 / _S3473;
    }
    float2  _S3476 = _S3472 * make_float2 (k_8);
    float u_95 = _S3476.x;
    float v_95 = _S3476.y;
    float r2_95 = u_95 * u_95 + v_95 * v_95;
    float2  _S3477 = _S3476 * make_float2 (1.0f + r2_95 * (dist_coeffs_17[int(0)] + r2_95 * (dist_coeffs_17[int(1)] + r2_95 * (dist_coeffs_17[int(2)] + r2_95 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_95 * v_95 + dist_coeffs_17[int(5)] * (r2_95 + 2.0f * u_95 * u_95) + dist_coeffs_17[int(6)] * r2_95, _S3448 * u_95 * v_95 + dist_coeffs_17[int(4)] * (r2_95 + 2.0f * v_95 * v_95) + dist_coeffs_17[int(7)] * r2_95);
    float2  _S3478 = _S3477 + make_float2 (dist_coeffs_17[int(8)] * _S3477.x + dist_coeffs_17[int(9)] * _S3477.y, 0.0f);
    (&_S3402)->_S3398 = make_float2 (fx_17 * _S3478.x + cx_14, fy_17 * _S3478.y + cy_14);
    float2  _S3479 = float2 {_S3440.p_0[int(5)].x, _S3440.p_0[int(5)].y};
    float _S3480 = length_0(_S3479);
    float _S3481 = _S3440.p_0[int(5)].z;
    float _S3482 = s_primal_ctx_atan2_0(_S3480, _S3481);
    if(_S3482 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3482 * _S3482 / 3.0f) / _S3481;
    }
    else
    {
        k_8 = _S3482 / _S3480;
    }
    float2  _S3483 = _S3479 * make_float2 (k_8);
    float u_96 = _S3483.x;
    float v_96 = _S3483.y;
    float r2_96 = u_96 * u_96 + v_96 * v_96;
    float2  _S3484 = _S3483 * make_float2 (1.0f + r2_96 * (dist_coeffs_17[int(0)] + r2_96 * (dist_coeffs_17[int(1)] + r2_96 * (dist_coeffs_17[int(2)] + r2_96 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_96 * v_96 + dist_coeffs_17[int(5)] * (r2_96 + 2.0f * u_96 * u_96) + dist_coeffs_17[int(6)] * r2_96, _S3448 * u_96 * v_96 + dist_coeffs_17[int(4)] * (r2_96 + 2.0f * v_96 * v_96) + dist_coeffs_17[int(7)] * r2_96);
    float2  _S3485 = _S3484 + make_float2 (dist_coeffs_17[int(8)] * _S3484.x + dist_coeffs_17[int(9)] * _S3484.y, 0.0f);
    (&_S3402)->_S3399 = make_float2 (fx_17 * _S3485.x + cx_14, fy_17 * _S3485.y + cy_14);
    float2  _S3486 = float2 {_S3441.p_0[int(6)].x, _S3441.p_0[int(6)].y};
    float _S3487 = length_0(_S3486);
    float _S3488 = _S3441.p_0[int(6)].z;
    float _S3489 = s_primal_ctx_atan2_0(_S3487, _S3488);
    if(_S3489 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3489 * _S3489 / 3.0f) / _S3488;
    }
    else
    {
        k_8 = _S3489 / _S3487;
    }
    float2  _S3490 = _S3486 * make_float2 (k_8);
    float u_97 = _S3490.x;
    float v_97 = _S3490.y;
    float r2_97 = u_97 * u_97 + v_97 * v_97;
    float2  _S3491 = _S3490 * make_float2 (1.0f + r2_97 * (dist_coeffs_17[int(0)] + r2_97 * (dist_coeffs_17[int(1)] + r2_97 * (dist_coeffs_17[int(2)] + r2_97 * dist_coeffs_17[int(3)])))) + make_float2 (_S3447 * u_97 * v_97 + dist_coeffs_17[int(5)] * (r2_97 + 2.0f * u_97 * u_97) + dist_coeffs_17[int(6)] * r2_97, _S3448 * u_97 * v_97 + dist_coeffs_17[int(4)] * (r2_97 + 2.0f * v_97 * v_97) + dist_coeffs_17[int(7)] * r2_97);
    float2  _S3492 = _S3491 + make_float2 (dist_coeffs_17[int(8)] * _S3491.x + dist_coeffs_17[int(9)] * _S3491.y, 0.0f);
    (&_S3402)->_S3400 = make_float2 (fx_17 * _S3492.x + cx_14, fy_17 * _S3492.y + cy_14);
    float3  mean_c_13 = s_primal_ctx_mul_0(R_13, mean_14) + t_13;
    float _S3493 = - in_opacity_13;
    float _S3494 = 1.0f + s_primal_ctx_exp_0(_S3493);
    float _S3495 = 1.0f / _S3494;
    float _S3496 = _S3494 * _S3494;
    float3  _S3497 = make_float3 (_S3413);
    float3  _S3498 = make_float3 (_S3416);
    float3  _S3499 = make_float3 (_S3419);
    float2  _S3500 = make_float2 (_S3411.w_mean_0[int(1)]) * _S3402._S3395 + make_float2 (_S3423.w_mean_0[int(2)]) * _S3402._S3396 + make_float2 (_S3425.w_mean_0[int(3)]) * _S3402._S3397 + make_float2 (_S3427.w_mean_0[int(4)]) * _S3402._S3398 + make_float2 (_S3429.w_mean_0[int(5)]) * _S3402._S3399 + make_float2 (_S3431.w_mean_0[int(6)]) * _S3402._S3400;
    float2  d_35 = _S3402._S3394 - _S3500;
    float _S3501 = d_35.x;
    float _S3502 = d_35.y;
    float _S3503 = _S3501 * _S3502;
    float2  d_36 = _S3402._S3395 - _S3500;
    float _S3504 = d_36.x;
    float _S3505 = d_36.y;
    float _S3506 = _S3504 * _S3505;
    float2  d_37 = _S3402._S3396 - _S3500;
    float _S3507 = d_37.x;
    float _S3508 = d_37.y;
    float _S3509 = _S3507 * _S3508;
    float2  d_38 = _S3402._S3397 - _S3500;
    float _S3510 = d_38.x;
    float _S3511 = d_38.y;
    float _S3512 = _S3510 * _S3511;
    float2  d_39 = _S3402._S3398 - _S3500;
    float _S3513 = d_39.x;
    float _S3514 = d_39.y;
    float _S3515 = _S3513 * _S3514;
    float2  d_40 = _S3402._S3399 - _S3500;
    float _S3516 = d_40.x;
    float _S3517 = d_40.y;
    float _S3518 = _S3516 * _S3517;
    float2  d_41 = _S3402._S3400 - _S3500;
    float _S3519 = d_41.x;
    float _S3520 = d_41.y;
    float _S3521 = _S3519 * _S3520;
    Matrix<float, 2, 2>  covar2d_9 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3501 * _S3501, _S3503, _S3503, _S3502 * _S3502) + makeMatrix<float, 2, 2> (_S3422.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3504 * _S3504, _S3506, _S3506, _S3505 * _S3505) + makeMatrix<float, 2, 2> (_S3424.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3507 * _S3507, _S3509, _S3509, _S3508 * _S3508) + makeMatrix<float, 2, 2> (_S3426.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3510 * _S3510, _S3512, _S3512, _S3511 * _S3511) + makeMatrix<float, 2, 2> (_S3428.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3513 * _S3513, _S3515, _S3515, _S3514 * _S3514) + makeMatrix<float, 2, 2> (_S3430.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3516 * _S3516, _S3518, _S3518, _S3517 * _S3517) + makeMatrix<float, 2, 2> (_S3432.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3519 * _S3519, _S3521, _S3521, _S3520 * _S3520);
    float eps2d_13;
    if(antialiased_13)
    {
        eps2d_13 = 0.10000000149011612f;
    }
    else
    {
        eps2d_13 = 0.30000001192092896f;
    }
    float _S3522 = covar2d_9.rows[int(0)].y * covar2d_9.rows[int(1)].x;
    float det_orig_13 = covar2d_9.rows[int(0)].x * covar2d_9.rows[int(1)].y - _S3522;
    float _S3523 = covar2d_9.rows[int(0)].x + eps2d_13;
    Matrix<float, 2, 2>  _S3524 = covar2d_9;
    *&(((&_S3524)->rows + (int(0)))->x) = _S3523;
    float _S3525 = covar2d_9.rows[int(1)].y + eps2d_13;
    *&(((&_S3524)->rows + (int(1)))->y) = _S3525;
    Matrix<float, 2, 2>  _S3526 = _S3524;
    Matrix<float, 2, 2>  _S3527 = _S3524;
    float det_blur_13 = _S3523 * _S3525 - _S3522;
    float _S3528 = det_orig_13 / det_blur_13;
    float _S3529 = det_blur_13 * det_blur_13;
    float _S3530 = (F32_max((0.0f), (_S3528)));
    float _S3531 = s_primal_ctx_sqrt_0(_S3530);
    float invdet_15 = 1.0f / det_blur_13;
    float _S3532 = - covar2d_9.rows[int(0)].y;
    float _S3533 = - covar2d_9.rows[int(1)].x;
    if(antialiased_13)
    {
        k_8 = _S3495 * _S3531;
    }
    else
    {
        k_8 = _S3495;
    }
    float _S3534 = k_8 / 0.00392156885936856f;
    float _S3535 = 2.0f * s_primal_ctx_log_0(_S3534);
    float _S3536 = s_primal_ctx_sqrt_0(_S3535);
    float _S3537 = _S3526.rows[int(0)].x;
    float _S3538 = _S3527.rows[int(1)].y;
    float3  campos_6 = - s_primal_ctx_mul_0(transpose_3(R_13), t_13);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3539;
    (&_S3539)->primal_0 = mean_14;
    (&_S3539)->differential_0 = _S3403;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3540;
    (&_S3540)->primal_0 = scale_13;
    (&_S3540)->differential_0 = _S3403;
    DiffPair_float_0 _S3541;
    (&_S3541)->primal_0 = in_opacity_13;
    (&_S3541)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3542;
    (&_S3542)->primal_0 = campos_6;
    (&_S3542)->differential_0 = _S3403;
    s_bwd_prop_view_radius_3dgs_0(&_S3539, &_S3540, &_S3541, &_S3542, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3543 = _S3539;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3544 = _S3540;
    DiffPair_float_0 _S3545 = _S3541;
    float2  _S3546 = _S3401;
    *&((&_S3546)->y) = v_conic_5.z;
    float2  _S3547 = _S3401;
    *&((&_S3547)->y) = v_conic_5.y;
    *&((&_S3547)->x) = v_conic_5.x;
    DiffPair_float_0 _S3548;
    (&_S3548)->primal_0 = _S3538;
    (&_S3548)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3548, 0.0f);
    DiffPair_float_0 _S3549;
    (&_S3549)->primal_0 = _S3537;
    (&_S3549)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3549, 0.0f);
    DiffPair_float_0 _S3550;
    (&_S3550)->primal_0 = 3.32999992370605469f;
    (&_S3550)->differential_0 = 0.0f;
    DiffPair_float_0 _S3551;
    (&_S3551)->primal_0 = _S3536;
    (&_S3551)->differential_0 = 0.0f;
    _d_min_0(&_S3550, &_S3551, 0.0f);
    DiffPair_float_0 _S3552;
    (&_S3552)->primal_0 = _S3535;
    (&_S3552)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3552, _S3551.differential_0);
    float _S3553 = 2.0f * _S3552.differential_0;
    DiffPair_float_0 _S3554;
    (&_S3554)->primal_0 = _S3534;
    (&_S3554)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3554, _S3553);
    float _S3555 = v_opacity_5 + 254.9999847412109375f * _S3554.differential_0;
    Matrix<float, 2, 2>  _S3556 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3557 = _S3556;
    _S3557[int(1)] = _S3546;
    _S3557[int(0)] = _S3547;
    Matrix<float, 2, 2>  _S3558 = _S3557;
    float2  _S3559 = make_float2 (0.0f, _S3548.differential_0);
    float2  _S3560 = make_float2 (_S3549.differential_0, 0.0f);
    if(antialiased_13)
    {
        float _S3561 = _S3531 * _S3555;
        k_8 = _S3495 * _S3555;
        eps2d_13 = _S3561;
    }
    else
    {
        k_8 = 0.0f;
        eps2d_13 = _S3555;
    }
    float _S3562 = invdet_15 * _S3558.rows[int(1)].y;
    float _S3563 = - (invdet_15 * _S3558.rows[int(1)].x);
    float _S3564 = - (invdet_15 * _S3558.rows[int(0)].y);
    float _S3565 = invdet_15 * _S3558.rows[int(0)].x;
    float _S3566 = - ((_S3523 * _S3558.rows[int(1)].y + _S3533 * _S3558.rows[int(1)].x + _S3532 * _S3558.rows[int(0)].y + _S3525 * _S3558.rows[int(0)].x) / _S3529);
    DiffPair_float_0 _S3567;
    (&_S3567)->primal_0 = _S3530;
    (&_S3567)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3567, k_8);
    DiffPair_float_0 _S3568;
    (&_S3568)->primal_0 = 0.0f;
    (&_S3568)->differential_0 = 0.0f;
    DiffPair_float_0 _S3569;
    (&_S3569)->primal_0 = _S3528;
    (&_S3569)->differential_0 = 0.0f;
    _d_max_0(&_S3568, &_S3569, _S3567.differential_0);
    float _S3570 = _S3569.differential_0 / _S3529;
    float s_diff_det_orig_T_5 = det_blur_13 * _S3570;
    float _S3571 = det_orig_13 * - _S3570 + _S3566;
    float _S3572 = - _S3571;
    float _S3573 = _S3523 * _S3571;
    float _S3574 = _S3525 * _S3571;
    Matrix<float, 2, 2>  _S3575 = _S3556;
    _S3575[int(1)] = _S3559;
    _S3575[int(0)] = _S3560;
    float _S3576 = _S3574 + _S3575.rows[int(0)].x + _S3562;
    float _S3577 = _S3572 + - s_diff_det_orig_T_5;
    float _S3578 = covar2d_9.rows[int(0)].y * _S3577 + _S3563;
    float _S3579 = covar2d_9.rows[int(1)].x * _S3577 + _S3564;
    float _S3580 = covar2d_9.rows[int(1)].y * s_diff_det_orig_T_5;
    float _S3581 = _S3573 + _S3575.rows[int(1)].y + _S3565 + covar2d_9.rows[int(0)].x * s_diff_det_orig_T_5;
    float2  _S3582 = _S3401;
    *&((&_S3582)->x) = _S3578;
    *&((&_S3582)->y) = _S3581;
    float _S3583 = _S3576 + _S3580;
    float2  _S3584 = _S3401;
    *&((&_S3584)->y) = _S3579;
    *&((&_S3584)->x) = _S3583;
    Matrix<float, 2, 2>  _S3585 = _S3556;
    _S3585[int(1)] = _S3582;
    _S3585[int(0)] = _S3584;
    Matrix<float, 3, 3>  _S3586 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3587;
    (&_S3587)->primal_0 = R_13;
    (&_S3587)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3588;
    (&_S3588)->primal_0 = _S3421;
    (&_S3588)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3587, &_S3588, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3589;
    (&_S3589)->primal_0 = R_13;
    (&_S3589)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3590;
    (&_S3590)->primal_0 = _S3418;
    (&_S3590)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3589, &_S3590, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3591;
    (&_S3591)->primal_0 = R_13;
    (&_S3591)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3592;
    (&_S3592)->primal_0 = _S3415;
    (&_S3592)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3591, &_S3592, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3593;
    (&_S3593)->primal_0 = R_13;
    (&_S3593)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3594;
    (&_S3594)->primal_0 = _S3420;
    (&_S3594)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3593, &_S3594, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3595;
    (&_S3595)->primal_0 = R_13;
    (&_S3595)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3596;
    (&_S3596)->primal_0 = _S3417;
    (&_S3596)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3595, &_S3596, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3597;
    (&_S3597)->primal_0 = R_13;
    (&_S3597)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3598;
    (&_S3598)->primal_0 = _S3414;
    (&_S3598)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3597, &_S3598, _S3403);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3599;
    (&_S3599)->primal_0 = R_13;
    (&_S3599)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3600;
    (&_S3600)->primal_0 = _S3433.p_0[0U];
    (&_S3600)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3599, &_S3600, _S3403);
    float3  _S3601 = - _S3588.differential_0 + _S3594.differential_0;
    float3  _S3602 = _S3499 * _S3601;
    float3  _S3603 = _S3407.rows[2U] * _S3601;
    float _S3604 = _S3412 * (_S3603.x + _S3603.y + _S3603.z);
    float3  _S3605 = - _S3590.differential_0 + _S3596.differential_0;
    float3  _S3606 = _S3498 * _S3605;
    float3  _S3607 = _S3407.rows[1U] * _S3605;
    float _S3608 = _S3412 * (_S3607.x + _S3607.y + _S3607.z);
    float3  _S3609 = - _S3592.differential_0 + _S3598.differential_0;
    float3  _S3610 = _S3497 * _S3609;
    float3  _S3611 = _S3407.rows[0U] * _S3609;
    float _S3612 = _S3412 * (_S3611.x + _S3611.y + _S3611.z);
    Matrix<float, 3, 3>  _S3613 = _S3586;
    _S3613[2U] = _S3602;
    _S3613[1U] = _S3606;
    _S3613[0U] = _S3610;
    Matrix<float, 3, 3>  _S3614 = transpose_3(transpose_3(_S3613));
    float _S3615 = 2.0f * - _S3614.rows[int(2)].z;
    float _S3616 = 2.0f * _S3614.rows[int(2)].y;
    float _S3617 = 2.0f * _S3614.rows[int(2)].x;
    float _S3618 = 2.0f * _S3614.rows[int(1)].z;
    float _S3619 = 2.0f * - _S3614.rows[int(1)].y;
    float _S3620 = 2.0f * _S3614.rows[int(1)].x;
    float _S3621 = 2.0f * _S3614.rows[int(0)].z;
    float _S3622 = 2.0f * _S3614.rows[int(0)].y;
    float _S3623 = 2.0f * - _S3614.rows[int(0)].x;
    float _S3624 = - _S3620 + _S3622;
    float _S3625 = _S3617 + - _S3621;
    float _S3626 = - _S3616 + _S3618;
    float _S3627 = _S3616 + _S3618;
    float _S3628 = _S3617 + _S3621;
    float _S3629 = _S3620 + _S3622;
    float _S3630 = _S3405.w * (_S3619 + _S3623);
    float _S3631 = _S3405.z * (_S3615 + _S3623);
    float _S3632 = _S3405.y * (_S3615 + _S3619);
    float _S3633 = _S3405.x * _S3624 + _S3405.z * _S3627 + _S3405.y * _S3628 + _S3630 + _S3630;
    float _S3634 = _S3405.x * _S3625 + _S3405.w * _S3627 + _S3405.y * _S3629 + _S3631 + _S3631;
    float _S3635 = _S3405.x * _S3626 + _S3405.w * _S3628 + _S3405.z * _S3629 + _S3632 + _S3632;
    float _S3636 = _S3405.w * _S3624 + _S3405.z * _S3625 + _S3405.y * _S3626;
    float4  _S3637 = make_float4 (0.0f);
    float4  _S3638 = _S3637;
    *&((&_S3638)->w) = _S3633;
    *&((&_S3638)->z) = _S3634;
    *&((&_S3638)->y) = _S3635;
    *&((&_S3638)->x) = _S3636;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3639;
    (&_S3639)->primal_0 = quat_13;
    (&_S3639)->differential_0 = _S3637;
    s_bwd_normalize_impl_0(&_S3639, _S3638);
    float3  _S3640 = _S3403;
    *&((&_S3640)->z) = _S3604;
    *&((&_S3640)->y) = _S3608;
    *&((&_S3640)->x) = _S3612;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3641;
    (&_S3641)->primal_0 = scale_13;
    (&_S3641)->differential_0 = _S3403;
    s_bwd_prop_exp_1(&_S3641, _S3640);
    float _S3642 = - (eps2d_13 / _S3496);
    DiffPair_float_0 _S3643;
    (&_S3643)->primal_0 = _S3493;
    (&_S3643)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3643, _S3642);
    float _S3644 = - _S3643.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3645;
    (&_S3645)->primal_0 = mean_c_13;
    (&_S3645)->differential_0 = _S3403;
    s_bwd_length_impl_0(&_S3645, v_depth_5);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3646;
    (&_S3646)->primal_0 = R_13;
    (&_S3646)->differential_0 = _S3586;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3647;
    (&_S3647)->primal_0 = mean_14;
    (&_S3647)->differential_0 = _S3403;
    s_bwd_prop_mul_3(&_S3646, &_S3647, _S3645.differential_0);
    Matrix<float, 3, 3>  _S3648 = _S3587.differential_0 + _S3589.differential_0 + _S3591.differential_0 + _S3593.differential_0 + _S3595.differential_0 + _S3597.differential_0 + _S3599.differential_0 + _S3646.differential_0;
    float _S3649 = _S3644 + _S3545.differential_0;
    float3  _S3650 = _S3641.differential_0 + _S3544.differential_0;
    *v_mean_5 = *v_mean_5 + (_S3588.differential_0 + _S3594.differential_0 + _S3590.differential_0 + _S3596.differential_0 + _S3592.differential_0 + _S3598.differential_0 + _S3647.differential_0 + _S3543.differential_0);
    *v_quat_5 = *v_quat_5 + _S3639.differential_0;
    *v_scale_5 = *v_scale_5 + _S3650;
    *v_in_opacity_5 = *v_in_opacity_5 + _S3649;
    *v_R_5 = *v_R_5 + _S3648;
    *v_t_5 = *v_t_5 + _S3645.differential_0;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_equisolid_projection_Intermediates_0
{
    float2  _S3651;
    float2  _S3652;
    float2  _S3653;
    float2  _S3654;
    float2  _S3655;
    float2  _S3656;
    float2  _S3657;
};

inline __device__ void projection_3dgut_equisolid_vjp(bool antialiased_14, float3  mean_15, float4  quat_14, float3  scale_14, float in_opacity_14, Matrix<float, 3, 3>  R_14, float3  t_14, float fx_18, float fy_18, float cx_15, float cy_15, FixedArray<float, 10>  dist_coeffs_18, uint image_width_14, uint image_height_14, float2  v_mean2d_6, float v_depth_6, float3  v_conic_6, float v_opacity_6, float3  * v_mean_6, float4  * v_quat_6, float3  * v_scale_6, float * v_in_opacity_6, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float2  _S3658 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_equisolid_projection_Intermediates_0 _S3659;
    (&_S3659)->_S3651 = _S3658;
    (&_S3659)->_S3652 = _S3658;
    (&_S3659)->_S3653 = _S3658;
    (&_S3659)->_S3654 = _S3658;
    (&_S3659)->_S3655 = _S3658;
    (&_S3659)->_S3656 = _S3658;
    (&_S3659)->_S3657 = _S3658;
    (&_S3659)->_S3651 = _S3658;
    (&_S3659)->_S3652 = _S3658;
    (&_S3659)->_S3653 = _S3658;
    (&_S3659)->_S3654 = _S3658;
    (&_S3659)->_S3655 = _S3658;
    (&_S3659)->_S3656 = _S3658;
    (&_S3659)->_S3657 = _S3658;
    float3  _S3660 = make_float3 (0.0f);
    float3  _S3661 = s_primal_ctx_exp_1(scale_14);
    float4  _S3662 = normalize_0(quat_14);
    float _S3663 = _S3662.y;
    float x2_14 = _S3663 * _S3663;
    float y2_14 = _S3662.z * _S3662.z;
    float z2_14 = _S3662.w * _S3662.w;
    float xy_14 = _S3662.y * _S3662.z;
    float xz_14 = _S3662.y * _S3662.w;
    float yz_14 = _S3662.z * _S3662.w;
    float wx_14 = _S3662.x * _S3662.y;
    float wy_14 = _S3662.x * _S3662.z;
    float wz_14 = _S3662.x * _S3662.w;
    Matrix<float, 3, 3>  _S3664 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_14), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_14), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14))));
    FixedArray<float3 , 7>  _S3665 = {
        _S3660, _S3660, _S3660, _S3660, _S3660, _S3660, _S3660
    };
    FixedArray<float, 7>  _S3666 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3667;
    (&_S3667)->p_0 = _S3665;
    (&_S3667)->w_mean_0 = _S3666;
    (&_S3667)->w_cov_0 = _S3666;
    (&_S3667)->p_0[int(0)] = mean_15;
    SigmaPoints_0 _S3668 = _S3667;
    (&_S3668)->w_mean_0[int(0)] = 0.0f;
    (&_S3668)->w_cov_0[int(0)] = 2.0f;
    float _S3669 = s_primal_ctx_sqrt_0(3.0f);
    float _S3670 = _S3669 * _S3661.x;
    float3  delta_18 = make_float3 (_S3670) * _S3664.rows[0U];
    float3  _S3671 = mean_15 + delta_18;
    (&_S3668)->p_0[1U] = _S3671;
    float3  _S3672 = mean_15 - delta_18;
    (&_S3668)->p_0[4U] = _S3672;
    float _S3673 = _S3669 * _S3661.y;
    float3  delta_19 = make_float3 (_S3673) * _S3664.rows[1U];
    float3  _S3674 = mean_15 + delta_19;
    (&_S3668)->p_0[2U] = _S3674;
    float3  _S3675 = mean_15 - delta_19;
    (&_S3668)->p_0[5U] = _S3675;
    float _S3676 = _S3669 * _S3661.z;
    float3  delta_20 = make_float3 (_S3676) * _S3664.rows[2U];
    float3  _S3677 = mean_15 + delta_20;
    (&_S3668)->p_0[3U] = _S3677;
    float3  _S3678 = mean_15 - delta_20;
    (&_S3668)->p_0[6U] = _S3678;
    (&_S3668)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3679 = _S3668;
    (&_S3679)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3680 = _S3679;
    (&_S3680)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3681 = _S3680;
    (&_S3681)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3682 = _S3681;
    (&_S3682)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3683 = _S3682;
    (&_S3683)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3684 = _S3683;
    (&_S3684)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3685 = _S3684;
    (&_S3685)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3686 = _S3685;
    (&_S3686)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3687 = _S3686;
    (&_S3687)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3688 = _S3687;
    (&_S3688)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3689 = _S3688;
    (&_S3689)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3690 = _S3667;
    float3  _S3691 = s_primal_ctx_mul_0(R_14, _S3667.p_0[0U]) + t_14;
    _S3667 = _S3689;
    (&_S3667)->p_0[0U] = _S3691;
    SigmaPoints_0 _S3692 = _S3667;
    (&_S3667)->p_0[1U] = s_primal_ctx_mul_0(R_14, _S3671) + t_14;
    SigmaPoints_0 _S3693 = _S3667;
    (&_S3667)->p_0[2U] = s_primal_ctx_mul_0(R_14, _S3674) + t_14;
    SigmaPoints_0 _S3694 = _S3667;
    (&_S3667)->p_0[3U] = s_primal_ctx_mul_0(R_14, _S3677) + t_14;
    SigmaPoints_0 _S3695 = _S3667;
    (&_S3667)->p_0[4U] = s_primal_ctx_mul_0(R_14, _S3672) + t_14;
    SigmaPoints_0 _S3696 = _S3667;
    (&_S3667)->p_0[5U] = s_primal_ctx_mul_0(R_14, _S3675) + t_14;
    SigmaPoints_0 _S3697 = _S3667;
    (&_S3667)->p_0[6U] = s_primal_ctx_mul_0(R_14, _S3678) + t_14;
    SigmaPoints_0 _S3698 = _S3667;
    float2  _S3699 = float2 {_S3692.p_0[int(0)].x, _S3692.p_0[int(0)].y};
    float _S3700 = length_0(_S3699);
    float _S3701 = _S3692.p_0[int(0)].z;
    float _S3702 = s_primal_ctx_atan2_0(_S3700, _S3701);
    float k_9;
    if(_S3700 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3702 * _S3702 / 24.0f) / _S3701;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3702) / _S3700;
    }
    float2  _S3703 = _S3699 * make_float2 (k_9);
    float u_98 = _S3703.x;
    float v_98 = _S3703.y;
    float r2_98 = u_98 * u_98 + v_98 * v_98;
    float _S3704 = 2.0f * dist_coeffs_18[int(4)];
    float _S3705 = 2.0f * dist_coeffs_18[int(5)];
    float2  _S3706 = _S3703 * make_float2 (1.0f + r2_98 * (dist_coeffs_18[int(0)] + r2_98 * (dist_coeffs_18[int(1)] + r2_98 * (dist_coeffs_18[int(2)] + r2_98 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_98 * v_98 + dist_coeffs_18[int(5)] * (r2_98 + 2.0f * u_98 * u_98) + dist_coeffs_18[int(6)] * r2_98, _S3705 * u_98 * v_98 + dist_coeffs_18[int(4)] * (r2_98 + 2.0f * v_98 * v_98) + dist_coeffs_18[int(7)] * r2_98);
    float2  _S3707 = _S3706 + make_float2 (dist_coeffs_18[int(8)] * _S3706.x + dist_coeffs_18[int(9)] * _S3706.y, 0.0f);
    (&_S3659)->_S3651 = make_float2 (fx_18 * _S3707.x + cx_15, fy_18 * _S3707.y + cy_15);
    float2  _S3708 = float2 {_S3693.p_0[int(1)].x, _S3693.p_0[int(1)].y};
    float _S3709 = length_0(_S3708);
    float _S3710 = _S3693.p_0[int(1)].z;
    float _S3711 = s_primal_ctx_atan2_0(_S3709, _S3710);
    if(_S3709 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3711 * _S3711 / 24.0f) / _S3710;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3711) / _S3709;
    }
    float2  _S3712 = _S3708 * make_float2 (k_9);
    float u_99 = _S3712.x;
    float v_99 = _S3712.y;
    float r2_99 = u_99 * u_99 + v_99 * v_99;
    float2  _S3713 = _S3712 * make_float2 (1.0f + r2_99 * (dist_coeffs_18[int(0)] + r2_99 * (dist_coeffs_18[int(1)] + r2_99 * (dist_coeffs_18[int(2)] + r2_99 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_99 * v_99 + dist_coeffs_18[int(5)] * (r2_99 + 2.0f * u_99 * u_99) + dist_coeffs_18[int(6)] * r2_99, _S3705 * u_99 * v_99 + dist_coeffs_18[int(4)] * (r2_99 + 2.0f * v_99 * v_99) + dist_coeffs_18[int(7)] * r2_99);
    float2  _S3714 = _S3713 + make_float2 (dist_coeffs_18[int(8)] * _S3713.x + dist_coeffs_18[int(9)] * _S3713.y, 0.0f);
    (&_S3659)->_S3652 = make_float2 (fx_18 * _S3714.x + cx_15, fy_18 * _S3714.y + cy_15);
    float2  _S3715 = float2 {_S3694.p_0[int(2)].x, _S3694.p_0[int(2)].y};
    float _S3716 = length_0(_S3715);
    float _S3717 = _S3694.p_0[int(2)].z;
    float _S3718 = s_primal_ctx_atan2_0(_S3716, _S3717);
    if(_S3716 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3718 * _S3718 / 24.0f) / _S3717;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3718) / _S3716;
    }
    float2  _S3719 = _S3715 * make_float2 (k_9);
    float u_100 = _S3719.x;
    float v_100 = _S3719.y;
    float r2_100 = u_100 * u_100 + v_100 * v_100;
    float2  _S3720 = _S3719 * make_float2 (1.0f + r2_100 * (dist_coeffs_18[int(0)] + r2_100 * (dist_coeffs_18[int(1)] + r2_100 * (dist_coeffs_18[int(2)] + r2_100 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_100 * v_100 + dist_coeffs_18[int(5)] * (r2_100 + 2.0f * u_100 * u_100) + dist_coeffs_18[int(6)] * r2_100, _S3705 * u_100 * v_100 + dist_coeffs_18[int(4)] * (r2_100 + 2.0f * v_100 * v_100) + dist_coeffs_18[int(7)] * r2_100);
    float2  _S3721 = _S3720 + make_float2 (dist_coeffs_18[int(8)] * _S3720.x + dist_coeffs_18[int(9)] * _S3720.y, 0.0f);
    (&_S3659)->_S3653 = make_float2 (fx_18 * _S3721.x + cx_15, fy_18 * _S3721.y + cy_15);
    float2  _S3722 = float2 {_S3695.p_0[int(3)].x, _S3695.p_0[int(3)].y};
    float _S3723 = length_0(_S3722);
    float _S3724 = _S3695.p_0[int(3)].z;
    float _S3725 = s_primal_ctx_atan2_0(_S3723, _S3724);
    if(_S3723 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3725 * _S3725 / 24.0f) / _S3724;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3725) / _S3723;
    }
    float2  _S3726 = _S3722 * make_float2 (k_9);
    float u_101 = _S3726.x;
    float v_101 = _S3726.y;
    float r2_101 = u_101 * u_101 + v_101 * v_101;
    float2  _S3727 = _S3726 * make_float2 (1.0f + r2_101 * (dist_coeffs_18[int(0)] + r2_101 * (dist_coeffs_18[int(1)] + r2_101 * (dist_coeffs_18[int(2)] + r2_101 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_101 * v_101 + dist_coeffs_18[int(5)] * (r2_101 + 2.0f * u_101 * u_101) + dist_coeffs_18[int(6)] * r2_101, _S3705 * u_101 * v_101 + dist_coeffs_18[int(4)] * (r2_101 + 2.0f * v_101 * v_101) + dist_coeffs_18[int(7)] * r2_101);
    float2  _S3728 = _S3727 + make_float2 (dist_coeffs_18[int(8)] * _S3727.x + dist_coeffs_18[int(9)] * _S3727.y, 0.0f);
    (&_S3659)->_S3654 = make_float2 (fx_18 * _S3728.x + cx_15, fy_18 * _S3728.y + cy_15);
    float2  _S3729 = float2 {_S3696.p_0[int(4)].x, _S3696.p_0[int(4)].y};
    float _S3730 = length_0(_S3729);
    float _S3731 = _S3696.p_0[int(4)].z;
    float _S3732 = s_primal_ctx_atan2_0(_S3730, _S3731);
    if(_S3730 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3732 * _S3732 / 24.0f) / _S3731;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3732) / _S3730;
    }
    float2  _S3733 = _S3729 * make_float2 (k_9);
    float u_102 = _S3733.x;
    float v_102 = _S3733.y;
    float r2_102 = u_102 * u_102 + v_102 * v_102;
    float2  _S3734 = _S3733 * make_float2 (1.0f + r2_102 * (dist_coeffs_18[int(0)] + r2_102 * (dist_coeffs_18[int(1)] + r2_102 * (dist_coeffs_18[int(2)] + r2_102 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_102 * v_102 + dist_coeffs_18[int(5)] * (r2_102 + 2.0f * u_102 * u_102) + dist_coeffs_18[int(6)] * r2_102, _S3705 * u_102 * v_102 + dist_coeffs_18[int(4)] * (r2_102 + 2.0f * v_102 * v_102) + dist_coeffs_18[int(7)] * r2_102);
    float2  _S3735 = _S3734 + make_float2 (dist_coeffs_18[int(8)] * _S3734.x + dist_coeffs_18[int(9)] * _S3734.y, 0.0f);
    (&_S3659)->_S3655 = make_float2 (fx_18 * _S3735.x + cx_15, fy_18 * _S3735.y + cy_15);
    float2  _S3736 = float2 {_S3697.p_0[int(5)].x, _S3697.p_0[int(5)].y};
    float _S3737 = length_0(_S3736);
    float _S3738 = _S3697.p_0[int(5)].z;
    float _S3739 = s_primal_ctx_atan2_0(_S3737, _S3738);
    if(_S3737 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3739 * _S3739 / 24.0f) / _S3738;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3739) / _S3737;
    }
    float2  _S3740 = _S3736 * make_float2 (k_9);
    float u_103 = _S3740.x;
    float v_103 = _S3740.y;
    float r2_103 = u_103 * u_103 + v_103 * v_103;
    float2  _S3741 = _S3740 * make_float2 (1.0f + r2_103 * (dist_coeffs_18[int(0)] + r2_103 * (dist_coeffs_18[int(1)] + r2_103 * (dist_coeffs_18[int(2)] + r2_103 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_103 * v_103 + dist_coeffs_18[int(5)] * (r2_103 + 2.0f * u_103 * u_103) + dist_coeffs_18[int(6)] * r2_103, _S3705 * u_103 * v_103 + dist_coeffs_18[int(4)] * (r2_103 + 2.0f * v_103 * v_103) + dist_coeffs_18[int(7)] * r2_103);
    float2  _S3742 = _S3741 + make_float2 (dist_coeffs_18[int(8)] * _S3741.x + dist_coeffs_18[int(9)] * _S3741.y, 0.0f);
    (&_S3659)->_S3656 = make_float2 (fx_18 * _S3742.x + cx_15, fy_18 * _S3742.y + cy_15);
    float2  _S3743 = float2 {_S3698.p_0[int(6)].x, _S3698.p_0[int(6)].y};
    float _S3744 = length_0(_S3743);
    float _S3745 = _S3698.p_0[int(6)].z;
    float _S3746 = s_primal_ctx_atan2_0(_S3744, _S3745);
    if(_S3744 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3746 * _S3746 / 24.0f) / _S3745;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3746) / _S3744;
    }
    float2  _S3747 = _S3743 * make_float2 (k_9);
    float u_104 = _S3747.x;
    float v_104 = _S3747.y;
    float r2_104 = u_104 * u_104 + v_104 * v_104;
    float2  _S3748 = _S3747 * make_float2 (1.0f + r2_104 * (dist_coeffs_18[int(0)] + r2_104 * (dist_coeffs_18[int(1)] + r2_104 * (dist_coeffs_18[int(2)] + r2_104 * dist_coeffs_18[int(3)])))) + make_float2 (_S3704 * u_104 * v_104 + dist_coeffs_18[int(5)] * (r2_104 + 2.0f * u_104 * u_104) + dist_coeffs_18[int(6)] * r2_104, _S3705 * u_104 * v_104 + dist_coeffs_18[int(4)] * (r2_104 + 2.0f * v_104 * v_104) + dist_coeffs_18[int(7)] * r2_104);
    float2  _S3749 = _S3748 + make_float2 (dist_coeffs_18[int(8)] * _S3748.x + dist_coeffs_18[int(9)] * _S3748.y, 0.0f);
    (&_S3659)->_S3657 = make_float2 (fx_18 * _S3749.x + cx_15, fy_18 * _S3749.y + cy_15);
    float3  mean_c_14 = s_primal_ctx_mul_0(R_14, mean_15) + t_14;
    float _S3750 = - in_opacity_14;
    float _S3751 = 1.0f + s_primal_ctx_exp_0(_S3750);
    float _S3752 = 1.0f / _S3751;
    float _S3753 = _S3751 * _S3751;
    float3  _S3754 = make_float3 (_S3670);
    float3  _S3755 = make_float3 (_S3673);
    float3  _S3756 = make_float3 (_S3676);
    float2  _S3757 = make_float2 (_S3668.w_mean_0[int(1)]) * _S3659._S3652 + make_float2 (_S3680.w_mean_0[int(2)]) * _S3659._S3653 + make_float2 (_S3682.w_mean_0[int(3)]) * _S3659._S3654 + make_float2 (_S3684.w_mean_0[int(4)]) * _S3659._S3655 + make_float2 (_S3686.w_mean_0[int(5)]) * _S3659._S3656 + make_float2 (_S3688.w_mean_0[int(6)]) * _S3659._S3657;
    float2  d_42 = _S3659._S3651 - _S3757;
    float _S3758 = d_42.x;
    float _S3759 = d_42.y;
    float _S3760 = _S3758 * _S3759;
    float2  d_43 = _S3659._S3652 - _S3757;
    float _S3761 = d_43.x;
    float _S3762 = d_43.y;
    float _S3763 = _S3761 * _S3762;
    float2  d_44 = _S3659._S3653 - _S3757;
    float _S3764 = d_44.x;
    float _S3765 = d_44.y;
    float _S3766 = _S3764 * _S3765;
    float2  d_45 = _S3659._S3654 - _S3757;
    float _S3767 = d_45.x;
    float _S3768 = d_45.y;
    float _S3769 = _S3767 * _S3768;
    float2  d_46 = _S3659._S3655 - _S3757;
    float _S3770 = d_46.x;
    float _S3771 = d_46.y;
    float _S3772 = _S3770 * _S3771;
    float2  d_47 = _S3659._S3656 - _S3757;
    float _S3773 = d_47.x;
    float _S3774 = d_47.y;
    float _S3775 = _S3773 * _S3774;
    float2  d_48 = _S3659._S3657 - _S3757;
    float _S3776 = d_48.x;
    float _S3777 = d_48.y;
    float _S3778 = _S3776 * _S3777;
    Matrix<float, 2, 2>  covar2d_10 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3758 * _S3758, _S3760, _S3760, _S3759 * _S3759) + makeMatrix<float, 2, 2> (_S3679.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3761 * _S3761, _S3763, _S3763, _S3762 * _S3762) + makeMatrix<float, 2, 2> (_S3681.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3764 * _S3764, _S3766, _S3766, _S3765 * _S3765) + makeMatrix<float, 2, 2> (_S3683.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3767 * _S3767, _S3769, _S3769, _S3768 * _S3768) + makeMatrix<float, 2, 2> (_S3685.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3770 * _S3770, _S3772, _S3772, _S3771 * _S3771) + makeMatrix<float, 2, 2> (_S3687.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3773 * _S3773, _S3775, _S3775, _S3774 * _S3774) + makeMatrix<float, 2, 2> (_S3689.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3776 * _S3776, _S3778, _S3778, _S3777 * _S3777);
    float eps2d_14;
    if(antialiased_14)
    {
        eps2d_14 = 0.10000000149011612f;
    }
    else
    {
        eps2d_14 = 0.30000001192092896f;
    }
    float _S3779 = covar2d_10.rows[int(0)].y * covar2d_10.rows[int(1)].x;
    float det_orig_14 = covar2d_10.rows[int(0)].x * covar2d_10.rows[int(1)].y - _S3779;
    float _S3780 = covar2d_10.rows[int(0)].x + eps2d_14;
    Matrix<float, 2, 2>  _S3781 = covar2d_10;
    *&(((&_S3781)->rows + (int(0)))->x) = _S3780;
    float _S3782 = covar2d_10.rows[int(1)].y + eps2d_14;
    *&(((&_S3781)->rows + (int(1)))->y) = _S3782;
    Matrix<float, 2, 2>  _S3783 = _S3781;
    Matrix<float, 2, 2>  _S3784 = _S3781;
    float det_blur_14 = _S3780 * _S3782 - _S3779;
    float _S3785 = det_orig_14 / det_blur_14;
    float _S3786 = det_blur_14 * det_blur_14;
    float _S3787 = (F32_max((0.0f), (_S3785)));
    float _S3788 = s_primal_ctx_sqrt_0(_S3787);
    float invdet_16 = 1.0f / det_blur_14;
    float _S3789 = - covar2d_10.rows[int(0)].y;
    float _S3790 = - covar2d_10.rows[int(1)].x;
    if(antialiased_14)
    {
        k_9 = _S3752 * _S3788;
    }
    else
    {
        k_9 = _S3752;
    }
    float _S3791 = k_9 / 0.00392156885936856f;
    float _S3792 = 2.0f * s_primal_ctx_log_0(_S3791);
    float _S3793 = s_primal_ctx_sqrt_0(_S3792);
    float _S3794 = _S3783.rows[int(0)].x;
    float _S3795 = _S3784.rows[int(1)].y;
    float3  campos_7 = - s_primal_ctx_mul_0(transpose_3(R_14), t_14);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3796;
    (&_S3796)->primal_0 = mean_15;
    (&_S3796)->differential_0 = _S3660;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3797;
    (&_S3797)->primal_0 = scale_14;
    (&_S3797)->differential_0 = _S3660;
    DiffPair_float_0 _S3798;
    (&_S3798)->primal_0 = in_opacity_14;
    (&_S3798)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3799;
    (&_S3799)->primal_0 = campos_7;
    (&_S3799)->differential_0 = _S3660;
    s_bwd_prop_view_radius_3dgs_0(&_S3796, &_S3797, &_S3798, &_S3799, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3800 = _S3796;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3801 = _S3797;
    DiffPair_float_0 _S3802 = _S3798;
    float2  _S3803 = _S3658;
    *&((&_S3803)->y) = v_conic_6.z;
    float2  _S3804 = _S3658;
    *&((&_S3804)->y) = v_conic_6.y;
    *&((&_S3804)->x) = v_conic_6.x;
    DiffPair_float_0 _S3805;
    (&_S3805)->primal_0 = _S3795;
    (&_S3805)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3805, 0.0f);
    DiffPair_float_0 _S3806;
    (&_S3806)->primal_0 = _S3794;
    (&_S3806)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3806, 0.0f);
    DiffPair_float_0 _S3807;
    (&_S3807)->primal_0 = 3.32999992370605469f;
    (&_S3807)->differential_0 = 0.0f;
    DiffPair_float_0 _S3808;
    (&_S3808)->primal_0 = _S3793;
    (&_S3808)->differential_0 = 0.0f;
    _d_min_0(&_S3807, &_S3808, 0.0f);
    DiffPair_float_0 _S3809;
    (&_S3809)->primal_0 = _S3792;
    (&_S3809)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3809, _S3808.differential_0);
    float _S3810 = 2.0f * _S3809.differential_0;
    DiffPair_float_0 _S3811;
    (&_S3811)->primal_0 = _S3791;
    (&_S3811)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3811, _S3810);
    float _S3812 = v_opacity_6 + 254.9999847412109375f * _S3811.differential_0;
    Matrix<float, 2, 2>  _S3813 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3814 = _S3813;
    _S3814[int(1)] = _S3803;
    _S3814[int(0)] = _S3804;
    Matrix<float, 2, 2>  _S3815 = _S3814;
    float2  _S3816 = make_float2 (0.0f, _S3805.differential_0);
    float2  _S3817 = make_float2 (_S3806.differential_0, 0.0f);
    if(antialiased_14)
    {
        float _S3818 = _S3788 * _S3812;
        k_9 = _S3752 * _S3812;
        eps2d_14 = _S3818;
    }
    else
    {
        k_9 = 0.0f;
        eps2d_14 = _S3812;
    }
    float _S3819 = invdet_16 * _S3815.rows[int(1)].y;
    float _S3820 = - (invdet_16 * _S3815.rows[int(1)].x);
    float _S3821 = - (invdet_16 * _S3815.rows[int(0)].y);
    float _S3822 = invdet_16 * _S3815.rows[int(0)].x;
    float _S3823 = - ((_S3780 * _S3815.rows[int(1)].y + _S3790 * _S3815.rows[int(1)].x + _S3789 * _S3815.rows[int(0)].y + _S3782 * _S3815.rows[int(0)].x) / _S3786);
    DiffPair_float_0 _S3824;
    (&_S3824)->primal_0 = _S3787;
    (&_S3824)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3824, k_9);
    DiffPair_float_0 _S3825;
    (&_S3825)->primal_0 = 0.0f;
    (&_S3825)->differential_0 = 0.0f;
    DiffPair_float_0 _S3826;
    (&_S3826)->primal_0 = _S3785;
    (&_S3826)->differential_0 = 0.0f;
    _d_max_0(&_S3825, &_S3826, _S3824.differential_0);
    float _S3827 = _S3826.differential_0 / _S3786;
    float s_diff_det_orig_T_6 = det_blur_14 * _S3827;
    float _S3828 = det_orig_14 * - _S3827 + _S3823;
    float _S3829 = - _S3828;
    float _S3830 = _S3780 * _S3828;
    float _S3831 = _S3782 * _S3828;
    Matrix<float, 2, 2>  _S3832 = _S3813;
    _S3832[int(1)] = _S3816;
    _S3832[int(0)] = _S3817;
    float _S3833 = _S3831 + _S3832.rows[int(0)].x + _S3819;
    float _S3834 = _S3829 + - s_diff_det_orig_T_6;
    float _S3835 = covar2d_10.rows[int(0)].y * _S3834 + _S3820;
    float _S3836 = covar2d_10.rows[int(1)].x * _S3834 + _S3821;
    float _S3837 = covar2d_10.rows[int(1)].y * s_diff_det_orig_T_6;
    float _S3838 = _S3830 + _S3832.rows[int(1)].y + _S3822 + covar2d_10.rows[int(0)].x * s_diff_det_orig_T_6;
    float2  _S3839 = _S3658;
    *&((&_S3839)->x) = _S3835;
    *&((&_S3839)->y) = _S3838;
    float _S3840 = _S3833 + _S3837;
    float2  _S3841 = _S3658;
    *&((&_S3841)->y) = _S3836;
    *&((&_S3841)->x) = _S3840;
    Matrix<float, 2, 2>  _S3842 = _S3813;
    _S3842[int(1)] = _S3839;
    _S3842[int(0)] = _S3841;
    Matrix<float, 3, 3>  _S3843 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3844;
    (&_S3844)->primal_0 = R_14;
    (&_S3844)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3845;
    (&_S3845)->primal_0 = _S3678;
    (&_S3845)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3844, &_S3845, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3846;
    (&_S3846)->primal_0 = R_14;
    (&_S3846)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3847;
    (&_S3847)->primal_0 = _S3675;
    (&_S3847)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3846, &_S3847, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3848;
    (&_S3848)->primal_0 = R_14;
    (&_S3848)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3849;
    (&_S3849)->primal_0 = _S3672;
    (&_S3849)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3848, &_S3849, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3850;
    (&_S3850)->primal_0 = R_14;
    (&_S3850)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3851;
    (&_S3851)->primal_0 = _S3677;
    (&_S3851)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3850, &_S3851, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3852;
    (&_S3852)->primal_0 = R_14;
    (&_S3852)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3853;
    (&_S3853)->primal_0 = _S3674;
    (&_S3853)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3852, &_S3853, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3854;
    (&_S3854)->primal_0 = R_14;
    (&_S3854)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3855;
    (&_S3855)->primal_0 = _S3671;
    (&_S3855)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3854, &_S3855, _S3660);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3856;
    (&_S3856)->primal_0 = R_14;
    (&_S3856)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3857;
    (&_S3857)->primal_0 = _S3690.p_0[0U];
    (&_S3857)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3856, &_S3857, _S3660);
    float3  _S3858 = - _S3845.differential_0 + _S3851.differential_0;
    float3  _S3859 = _S3756 * _S3858;
    float3  _S3860 = _S3664.rows[2U] * _S3858;
    float _S3861 = _S3669 * (_S3860.x + _S3860.y + _S3860.z);
    float3  _S3862 = - _S3847.differential_0 + _S3853.differential_0;
    float3  _S3863 = _S3755 * _S3862;
    float3  _S3864 = _S3664.rows[1U] * _S3862;
    float _S3865 = _S3669 * (_S3864.x + _S3864.y + _S3864.z);
    float3  _S3866 = - _S3849.differential_0 + _S3855.differential_0;
    float3  _S3867 = _S3754 * _S3866;
    float3  _S3868 = _S3664.rows[0U] * _S3866;
    float _S3869 = _S3669 * (_S3868.x + _S3868.y + _S3868.z);
    Matrix<float, 3, 3>  _S3870 = _S3843;
    _S3870[2U] = _S3859;
    _S3870[1U] = _S3863;
    _S3870[0U] = _S3867;
    Matrix<float, 3, 3>  _S3871 = transpose_3(transpose_3(_S3870));
    float _S3872 = 2.0f * - _S3871.rows[int(2)].z;
    float _S3873 = 2.0f * _S3871.rows[int(2)].y;
    float _S3874 = 2.0f * _S3871.rows[int(2)].x;
    float _S3875 = 2.0f * _S3871.rows[int(1)].z;
    float _S3876 = 2.0f * - _S3871.rows[int(1)].y;
    float _S3877 = 2.0f * _S3871.rows[int(1)].x;
    float _S3878 = 2.0f * _S3871.rows[int(0)].z;
    float _S3879 = 2.0f * _S3871.rows[int(0)].y;
    float _S3880 = 2.0f * - _S3871.rows[int(0)].x;
    float _S3881 = - _S3877 + _S3879;
    float _S3882 = _S3874 + - _S3878;
    float _S3883 = - _S3873 + _S3875;
    float _S3884 = _S3873 + _S3875;
    float _S3885 = _S3874 + _S3878;
    float _S3886 = _S3877 + _S3879;
    float _S3887 = _S3662.w * (_S3876 + _S3880);
    float _S3888 = _S3662.z * (_S3872 + _S3880);
    float _S3889 = _S3662.y * (_S3872 + _S3876);
    float _S3890 = _S3662.x * _S3881 + _S3662.z * _S3884 + _S3662.y * _S3885 + _S3887 + _S3887;
    float _S3891 = _S3662.x * _S3882 + _S3662.w * _S3884 + _S3662.y * _S3886 + _S3888 + _S3888;
    float _S3892 = _S3662.x * _S3883 + _S3662.w * _S3885 + _S3662.z * _S3886 + _S3889 + _S3889;
    float _S3893 = _S3662.w * _S3881 + _S3662.z * _S3882 + _S3662.y * _S3883;
    float4  _S3894 = make_float4 (0.0f);
    float4  _S3895 = _S3894;
    *&((&_S3895)->w) = _S3890;
    *&((&_S3895)->z) = _S3891;
    *&((&_S3895)->y) = _S3892;
    *&((&_S3895)->x) = _S3893;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3896;
    (&_S3896)->primal_0 = quat_14;
    (&_S3896)->differential_0 = _S3894;
    s_bwd_normalize_impl_0(&_S3896, _S3895);
    float3  _S3897 = _S3660;
    *&((&_S3897)->z) = _S3861;
    *&((&_S3897)->y) = _S3865;
    *&((&_S3897)->x) = _S3869;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3898;
    (&_S3898)->primal_0 = scale_14;
    (&_S3898)->differential_0 = _S3660;
    s_bwd_prop_exp_1(&_S3898, _S3897);
    float _S3899 = - (eps2d_14 / _S3753);
    DiffPair_float_0 _S3900;
    (&_S3900)->primal_0 = _S3750;
    (&_S3900)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3900, _S3899);
    float _S3901 = - _S3900.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3902;
    (&_S3902)->primal_0 = mean_c_14;
    (&_S3902)->differential_0 = _S3660;
    s_bwd_length_impl_0(&_S3902, v_depth_6);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3903;
    (&_S3903)->primal_0 = R_14;
    (&_S3903)->differential_0 = _S3843;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3904;
    (&_S3904)->primal_0 = mean_15;
    (&_S3904)->differential_0 = _S3660;
    s_bwd_prop_mul_3(&_S3903, &_S3904, _S3902.differential_0);
    Matrix<float, 3, 3>  _S3905 = _S3844.differential_0 + _S3846.differential_0 + _S3848.differential_0 + _S3850.differential_0 + _S3852.differential_0 + _S3854.differential_0 + _S3856.differential_0 + _S3903.differential_0;
    float _S3906 = _S3901 + _S3802.differential_0;
    float3  _S3907 = _S3898.differential_0 + _S3801.differential_0;
    *v_mean_6 = *v_mean_6 + (_S3845.differential_0 + _S3851.differential_0 + _S3847.differential_0 + _S3853.differential_0 + _S3849.differential_0 + _S3855.differential_0 + _S3904.differential_0 + _S3800.differential_0);
    *v_quat_6 = *v_quat_6 + _S3896.differential_0;
    *v_scale_6 = *v_scale_6 + _S3907;
    *v_in_opacity_6 = *v_in_opacity_6 + _S3906;
    *v_R_6 = *v_R_6 + _S3905;
    *v_t_6 = *v_t_6 + _S3902.differential_0;
    return;
}

inline __device__ void projection_3dgut_equirect_vjp(bool antialiased_15, float3  mean_16, float4  quat_15, float3  scale_15, float in_opacity_15, Matrix<float, 3, 3>  R_15, float3  t_15, float fx_19, float fy_19, float cx_16, float cy_16, FixedArray<float, 10>  dist_coeffs_19, uint image_width_15, uint image_height_15, float2  v_mean2d_7, float v_depth_7, float3  v_conic_7, float v_opacity_7, float3  * v_mean_7, float4  * v_quat_7, float3  * v_scale_7, float * v_in_opacity_7, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  _S3908 = make_float3 (0.0f);
    float3  mean_c_15 = s_primal_ctx_mul_0(R_15, mean_16) + t_15;
    float _S3909 = - in_opacity_15;
    float _S3910 = 1.0f + s_primal_ctx_exp_0(_S3909);
    float _S3911 = 1.0f / _S3910;
    float _S3912 = _S3910 * _S3910;
    float3  _S3913 = s_primal_ctx_exp_1(scale_15);
    float4  _S3914 = normalize_0(quat_15);
    float _S3915 = _S3914.y;
    float x2_15 = _S3915 * _S3915;
    float y2_15 = _S3914.z * _S3914.z;
    float z2_15 = _S3914.w * _S3914.w;
    float xy_15 = _S3914.y * _S3914.z;
    float xz_15 = _S3914.y * _S3914.w;
    float yz_15 = _S3914.z * _S3914.w;
    float wx_15 = _S3914.x * _S3914.y;
    float wy_15 = _S3914.x * _S3914.z;
    float wz_15 = _S3914.x * _S3914.w;
    Matrix<float, 3, 3>  _S3916 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_15), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_15), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15))));
    FixedArray<float3 , 7>  _S3917 = {
        _S3908, _S3908, _S3908, _S3908, _S3908, _S3908, _S3908
    };
    FixedArray<float, 7>  _S3918 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3919;
    (&_S3919)->p_0 = _S3917;
    (&_S3919)->w_mean_0 = _S3918;
    (&_S3919)->w_cov_0 = _S3918;
    (&_S3919)->p_0[int(0)] = mean_16;
    SigmaPoints_0 _S3920 = _S3919;
    (&_S3920)->w_mean_0[int(0)] = 0.0f;
    (&_S3920)->w_cov_0[int(0)] = 2.0f;
    float _S3921 = s_primal_ctx_sqrt_0(3.0f);
    float _S3922 = _S3921 * _S3913.x;
    float3  _S3923 = make_float3 (_S3922);
    float3  delta_21 = make_float3 (_S3922) * _S3916.rows[0U];
    float3  _S3924 = mean_16 + delta_21;
    (&_S3920)->p_0[1U] = _S3924;
    float3  _S3925 = mean_16 - delta_21;
    (&_S3920)->p_0[4U] = _S3925;
    float _S3926 = _S3921 * _S3913.y;
    float3  _S3927 = make_float3 (_S3926);
    float3  delta_22 = make_float3 (_S3926) * _S3916.rows[1U];
    float3  _S3928 = mean_16 + delta_22;
    (&_S3920)->p_0[2U] = _S3928;
    float3  _S3929 = mean_16 - delta_22;
    (&_S3920)->p_0[5U] = _S3929;
    float _S3930 = _S3921 * _S3913.z;
    float3  _S3931 = make_float3 (_S3930);
    float3  delta_23 = make_float3 (_S3930) * _S3916.rows[2U];
    float3  _S3932 = mean_16 + delta_23;
    (&_S3920)->p_0[3U] = _S3932;
    float3  _S3933 = mean_16 - delta_23;
    (&_S3920)->p_0[6U] = _S3933;
    (&_S3920)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3934 = _S3920;
    (&_S3934)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3935 = _S3934;
    (&_S3935)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3936 = _S3935;
    (&_S3936)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3937 = _S3936;
    (&_S3937)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3938 = _S3937;
    (&_S3938)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3939 = _S3938;
    (&_S3939)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3940 = _S3939;
    (&_S3940)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3941 = _S3940;
    (&_S3941)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3942 = _S3941;
    (&_S3942)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3943 = _S3942;
    (&_S3943)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3944 = _S3943;
    (&_S3944)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3945 = _S3919;
    float3  _S3946 = s_primal_ctx_mul_0(R_15, _S3919.p_0[0U]) + t_15;
    _S3919 = _S3944;
    (&_S3919)->p_0[0U] = _S3946;
    SigmaPoints_0 _S3947 = _S3919;
    (&_S3919)->p_0[1U] = s_primal_ctx_mul_0(R_15, _S3924) + t_15;
    SigmaPoints_0 _S3948 = _S3919;
    (&_S3919)->p_0[2U] = s_primal_ctx_mul_0(R_15, _S3928) + t_15;
    SigmaPoints_0 _S3949 = _S3919;
    (&_S3919)->p_0[3U] = s_primal_ctx_mul_0(R_15, _S3932) + t_15;
    SigmaPoints_0 _S3950 = _S3919;
    (&_S3919)->p_0[4U] = s_primal_ctx_mul_0(R_15, _S3925) + t_15;
    SigmaPoints_0 _S3951 = _S3919;
    (&_S3919)->p_0[5U] = s_primal_ctx_mul_0(R_15, _S3929) + t_15;
    SigmaPoints_0 _S3952 = _S3919;
    (&_S3919)->p_0[6U] = s_primal_ctx_mul_0(R_15, _S3933) + t_15;
    float _S3953 = fx_19 * s_primal_ctx_atan2_0(_S3947.p_0[int(0)].x, _S3947.p_0[int(0)].z) + cx_16;
    float _S3954 = fx_19 * s_primal_ctx_atan2_0(_S3948.p_0[int(1)].x, _S3948.p_0[int(1)].z) + cx_16;
    float2  _S3955 = make_float2 (_S3954, fy_19 * s_primal_ctx_atan2_0(_S3948.p_0[int(1)].y, length_0(float2 {_S3948.p_0[int(1)].x, _S3948.p_0[int(1)].z})) + cy_16);
    float2  _S3956 = make_float2 (fx_19 * s_primal_ctx_atan2_0(_S3949.p_0[int(2)].x, _S3949.p_0[int(2)].z) + cx_16, fy_19 * s_primal_ctx_atan2_0(_S3949.p_0[int(2)].y, length_0(float2 {_S3949.p_0[int(2)].x, _S3949.p_0[int(2)].z})) + cy_16);
    float2  _S3957 = make_float2 (fx_19 * s_primal_ctx_atan2_0(_S3950.p_0[int(3)].x, _S3950.p_0[int(3)].z) + cx_16, fy_19 * s_primal_ctx_atan2_0(_S3950.p_0[int(3)].y, length_0(float2 {_S3950.p_0[int(3)].x, _S3950.p_0[int(3)].z})) + cy_16);
    float2  _S3958 = make_float2 (fx_19 * s_primal_ctx_atan2_0(_S3951.p_0[int(4)].x, _S3951.p_0[int(4)].z) + cx_16, fy_19 * s_primal_ctx_atan2_0(_S3951.p_0[int(4)].y, length_0(float2 {_S3951.p_0[int(4)].x, _S3951.p_0[int(4)].z})) + cy_16);
    float2  _S3959 = make_float2 (fx_19 * s_primal_ctx_atan2_0(_S3952.p_0[int(5)].x, _S3952.p_0[int(5)].z) + cx_16, fy_19 * s_primal_ctx_atan2_0(_S3952.p_0[int(5)].y, length_0(float2 {_S3952.p_0[int(5)].x, _S3952.p_0[int(5)].z})) + cy_16);
    float2  _S3960 = make_float2 (fx_19 * s_primal_ctx_atan2_0(_S3919.p_0[int(6)].x, _S3919.p_0[int(6)].z) + cx_16, fy_19 * s_primal_ctx_atan2_0(_S3919.p_0[int(6)].y, length_0(float2 {_S3919.p_0[int(6)].x, _S3919.p_0[int(6)].z})) + cy_16);
    float _S3961 = fx_19 * 6.28318548202514648f;
    float du_6 = _S3954 - _S3953;
    float _S3962 = _S3953 + (du_6 - _S3961 * (F32_round((du_6 / _S3961))));
    FixedArray<float2 , 7>  _S3963;
    _S3963[int(0)] = make_float2 (_S3953, fy_19 * s_primal_ctx_atan2_0(_S3947.p_0[int(0)].y, length_0(float2 {_S3947.p_0[int(0)].x, _S3947.p_0[int(0)].z})) + cy_16);
    _S3963[int(1)] = _S3955;
    _S3963[int(2)] = _S3956;
    _S3963[int(3)] = _S3957;
    _S3963[int(4)] = _S3958;
    _S3963[int(5)] = _S3959;
    _S3963[int(6)] = _S3960;
    *&((&_S3963[int(1)])->x) = _S3962;
    float du_7 = _S3963[int(2)].x - _S3953;
    *&((&_S3963[int(2)])->x) = _S3953 + (du_7 - _S3961 * (F32_round((du_7 / _S3961))));
    float du_8 = _S3963[int(3)].x - _S3953;
    *&((&_S3963[int(3)])->x) = _S3953 + (du_8 - _S3961 * (F32_round((du_8 / _S3961))));
    float du_9 = _S3963[int(4)].x - _S3953;
    *&((&_S3963[int(4)])->x) = _S3953 + (du_9 - _S3961 * (F32_round((du_9 / _S3961))));
    float du_10 = _S3963[int(5)].x - _S3953;
    *&((&_S3963[int(5)])->x) = _S3953 + (du_10 - _S3961 * (F32_round((du_10 / _S3961))));
    float du_11 = _S3963[int(6)].x - _S3953;
    *&((&_S3963[int(6)])->x) = _S3953 + (du_11 - _S3961 * (F32_round((du_11 / _S3961))));
    float2  _S3964 = make_float2 (_S3920.w_mean_0[int(1)]);
    float2  _S3965 = make_float2 (_S3935.w_mean_0[int(2)]);
    float2  _S3966 = make_float2 (_S3937.w_mean_0[int(3)]);
    float2  _S3967 = make_float2 (_S3939.w_mean_0[int(4)]);
    float2  _S3968 = make_float2 (_S3941.w_mean_0[int(5)]);
    float2  _S3969 = make_float2 (_S3943.w_mean_0[int(6)]);
    float2  _S3970 = make_float2 (_S3920.w_mean_0[int(1)]) * _S3963[int(1)] + make_float2 (_S3935.w_mean_0[int(2)]) * _S3963[int(2)] + make_float2 (_S3937.w_mean_0[int(3)]) * _S3963[int(3)] + make_float2 (_S3939.w_mean_0[int(4)]) * _S3963[int(4)] + make_float2 (_S3941.w_mean_0[int(5)]) * _S3963[int(5)] + make_float2 (_S3943.w_mean_0[int(6)]) * _S3963[int(6)];
    float2  d_49 = _S3963[int(0)] - _S3970;
    float _S3971 = d_49.x;
    float _S3972 = d_49.y;
    float _S3973 = _S3971 * _S3972;
    float2  d_50 = _S3963[int(1)] - _S3970;
    Matrix<float, 2, 2>  _S3974 = makeMatrix<float, 2, 2> (_S3934.w_cov_0[int(1)]);
    float _S3975 = d_50.x;
    float _S3976 = d_50.y;
    float _S3977 = _S3975 * _S3976;
    float2  d_51 = _S3963[int(2)] - _S3970;
    Matrix<float, 2, 2>  _S3978 = makeMatrix<float, 2, 2> (_S3936.w_cov_0[int(2)]);
    float _S3979 = d_51.x;
    float _S3980 = d_51.y;
    float _S3981 = _S3979 * _S3980;
    float2  d_52 = _S3963[int(3)] - _S3970;
    Matrix<float, 2, 2>  _S3982 = makeMatrix<float, 2, 2> (_S3938.w_cov_0[int(3)]);
    float _S3983 = d_52.x;
    float _S3984 = d_52.y;
    float _S3985 = _S3983 * _S3984;
    float2  d_53 = _S3963[int(4)] - _S3970;
    Matrix<float, 2, 2>  _S3986 = makeMatrix<float, 2, 2> (_S3940.w_cov_0[int(4)]);
    float _S3987 = d_53.x;
    float _S3988 = d_53.y;
    float _S3989 = _S3987 * _S3988;
    float2  d_54 = _S3963[int(5)] - _S3970;
    Matrix<float, 2, 2>  _S3990 = makeMatrix<float, 2, 2> (_S3942.w_cov_0[int(5)]);
    float _S3991 = d_54.x;
    float _S3992 = d_54.y;
    float _S3993 = _S3991 * _S3992;
    float2  d_55 = _S3963[int(6)] - _S3970;
    Matrix<float, 2, 2>  _S3994 = makeMatrix<float, 2, 2> (_S3944.w_cov_0[int(6)]);
    float _S3995 = d_55.x;
    float _S3996 = d_55.y;
    float _S3997 = _S3995 * _S3996;
    Matrix<float, 2, 2>  covar2d_11 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3971 * _S3971, _S3973, _S3973, _S3972 * _S3972) + makeMatrix<float, 2, 2> (_S3934.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3975 * _S3975, _S3977, _S3977, _S3976 * _S3976) + makeMatrix<float, 2, 2> (_S3936.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3979 * _S3979, _S3981, _S3981, _S3980 * _S3980) + makeMatrix<float, 2, 2> (_S3938.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3983 * _S3983, _S3985, _S3985, _S3984 * _S3984) + makeMatrix<float, 2, 2> (_S3940.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3987 * _S3987, _S3989, _S3989, _S3988 * _S3988) + makeMatrix<float, 2, 2> (_S3942.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3991 * _S3991, _S3993, _S3993, _S3992 * _S3992) + makeMatrix<float, 2, 2> (_S3944.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3995 * _S3995, _S3997, _S3997, _S3996 * _S3996);
    float eps2d_15;
    if(antialiased_15)
    {
        eps2d_15 = 0.10000000149011612f;
    }
    else
    {
        eps2d_15 = 0.30000001192092896f;
    }
    float _S3998 = covar2d_11.rows[int(0)].y * covar2d_11.rows[int(1)].x;
    float det_orig_15 = covar2d_11.rows[int(0)].x * covar2d_11.rows[int(1)].y - _S3998;
    float _S3999 = covar2d_11.rows[int(0)].x + eps2d_15;
    Matrix<float, 2, 2>  _S4000 = covar2d_11;
    *&(((&_S4000)->rows + (int(0)))->x) = _S3999;
    float _S4001 = covar2d_11.rows[int(1)].y + eps2d_15;
    *&(((&_S4000)->rows + (int(1)))->y) = _S4001;
    Matrix<float, 2, 2>  _S4002 = _S4000;
    Matrix<float, 2, 2>  _S4003 = _S4000;
    float det_blur_15 = _S3999 * _S4001 - _S3998;
    float _S4004 = det_orig_15 / det_blur_15;
    float _S4005 = det_blur_15 * det_blur_15;
    float _S4006 = (F32_max((0.0f), (_S4004)));
    float _S4007 = s_primal_ctx_sqrt_0(_S4006);
    float invdet_17 = 1.0f / det_blur_15;
    float _S4008 = - covar2d_11.rows[int(0)].y;
    float _S4009 = - covar2d_11.rows[int(1)].x;
    if(antialiased_15)
    {
        eps2d_15 = _S3911 * _S4007;
    }
    else
    {
        eps2d_15 = _S3911;
    }
    float _S4010 = eps2d_15 / 0.00392156885936856f;
    float _S4011 = 2.0f * s_primal_ctx_log_0(_S4010);
    float _S4012 = s_primal_ctx_sqrt_0(_S4011);
    float _S4013 = _S4002.rows[int(0)].x;
    float _S4014 = _S4003.rows[int(1)].y;
    float3  campos_8 = - s_primal_ctx_mul_0(transpose_3(R_15), t_15);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4015;
    (&_S4015)->primal_0 = mean_16;
    (&_S4015)->differential_0 = _S3908;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4016;
    (&_S4016)->primal_0 = scale_15;
    (&_S4016)->differential_0 = _S3908;
    DiffPair_float_0 _S4017;
    (&_S4017)->primal_0 = in_opacity_15;
    (&_S4017)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4018;
    (&_S4018)->primal_0 = campos_8;
    (&_S4018)->differential_0 = _S3908;
    s_bwd_prop_view_radius_3dgs_0(&_S4015, &_S4016, &_S4017, &_S4018, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4019 = _S4015;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4020 = _S4016;
    DiffPair_float_0 _S4021 = _S4017;
    float2  _S4022 = make_float2 (0.0f);
    float2  _S4023 = _S4022;
    *&((&_S4023)->y) = v_conic_7.z;
    float2  _S4024 = _S4022;
    *&((&_S4024)->y) = v_conic_7.y;
    *&((&_S4024)->x) = v_conic_7.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4025;
    (&_S4025)->primal_0 = mean_c_15;
    (&_S4025)->differential_0 = _S3908;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4026;
    (&_S4026)->primal_0 = mean_c_15;
    (&_S4026)->differential_0 = _S3908;
    s_bwd_prop_dot_0(&_S4025, &_S4026, 0.0f);
    DiffPair_float_0 _S4027;
    (&_S4027)->primal_0 = _S4014;
    (&_S4027)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4027, 0.0f);
    DiffPair_float_0 _S4028;
    (&_S4028)->primal_0 = _S4013;
    (&_S4028)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4028, 0.0f);
    DiffPair_float_0 _S4029;
    (&_S4029)->primal_0 = 3.32999992370605469f;
    (&_S4029)->differential_0 = 0.0f;
    DiffPair_float_0 _S4030;
    (&_S4030)->primal_0 = _S4012;
    (&_S4030)->differential_0 = 0.0f;
    _d_min_0(&_S4029, &_S4030, 0.0f);
    DiffPair_float_0 _S4031;
    (&_S4031)->primal_0 = _S4011;
    (&_S4031)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4031, _S4030.differential_0);
    float _S4032 = 2.0f * _S4031.differential_0;
    DiffPair_float_0 _S4033;
    (&_S4033)->primal_0 = _S4010;
    (&_S4033)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4033, _S4032);
    float _S4034 = v_opacity_7 + 254.9999847412109375f * _S4033.differential_0;
    Matrix<float, 2, 2>  _S4035 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S4036 = _S4035;
    _S4036[int(1)] = _S4023;
    _S4036[int(0)] = _S4024;
    Matrix<float, 2, 2>  _S4037 = _S4036;
    float3  _S4038 = _S4026.differential_0 + _S4025.differential_0;
    float2  _S4039 = make_float2 (0.0f, _S4027.differential_0);
    float2  _S4040 = make_float2 (_S4028.differential_0, 0.0f);
    float _S4041;
    if(antialiased_15)
    {
        float _S4042 = _S4007 * _S4034;
        eps2d_15 = _S3911 * _S4034;
        _S4041 = _S4042;
    }
    else
    {
        eps2d_15 = 0.0f;
        _S4041 = _S4034;
    }
    float _S4043 = invdet_17 * _S4037.rows[int(1)].y;
    float _S4044 = - (invdet_17 * _S4037.rows[int(1)].x);
    float _S4045 = - (invdet_17 * _S4037.rows[int(0)].y);
    float _S4046 = invdet_17 * _S4037.rows[int(0)].x;
    float _S4047 = - ((_S3999 * _S4037.rows[int(1)].y + _S4009 * _S4037.rows[int(1)].x + _S4008 * _S4037.rows[int(0)].y + _S4001 * _S4037.rows[int(0)].x) / _S4005);
    DiffPair_float_0 _S4048;
    (&_S4048)->primal_0 = _S4006;
    (&_S4048)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4048, eps2d_15);
    DiffPair_float_0 _S4049;
    (&_S4049)->primal_0 = 0.0f;
    (&_S4049)->differential_0 = 0.0f;
    DiffPair_float_0 _S4050;
    (&_S4050)->primal_0 = _S4004;
    (&_S4050)->differential_0 = 0.0f;
    _d_max_0(&_S4049, &_S4050, _S4048.differential_0);
    float _S4051 = _S4050.differential_0 / _S4005;
    float s_diff_det_orig_T_7 = det_blur_15 * _S4051;
    float _S4052 = det_orig_15 * - _S4051 + _S4047;
    float _S4053 = - _S4052;
    float _S4054 = _S3999 * _S4052;
    float _S4055 = _S4001 * _S4052;
    Matrix<float, 2, 2>  _S4056 = _S4035;
    _S4056[int(1)] = _S4039;
    _S4056[int(0)] = _S4040;
    _S4000 = _S4056;
    *&(((&_S4000)->rows + (int(1)))->y) = 0.0f;
    float _S4057 = _S4054 + _S4056.rows[int(1)].y + _S4046;
    *&(((&_S4000)->rows + (int(0)))->x) = 0.0f;
    float _S4058 = _S4055 + _S4056.rows[int(0)].x + _S4043;
    float _S4059 = _S4053 + - s_diff_det_orig_T_7;
    float _S4060 = covar2d_11.rows[int(0)].y * _S4059 + _S4044;
    float _S4061 = covar2d_11.rows[int(1)].x * _S4059 + _S4045;
    float _S4062 = covar2d_11.rows[int(1)].y * s_diff_det_orig_T_7;
    float _S4063 = _S4057 + covar2d_11.rows[int(0)].x * s_diff_det_orig_T_7;
    float2  _S4064 = _S4022;
    *&((&_S4064)->x) = _S4060;
    *&((&_S4064)->y) = _S4063;
    float _S4065 = _S4058 + _S4062;
    float2  _S4066 = _S4022;
    *&((&_S4066)->y) = _S4061;
    *&((&_S4066)->x) = _S4065;
    Matrix<float, 2, 2>  _S4067 = _S4035;
    _S4067[int(1)] = _S4064;
    _S4067[int(0)] = _S4066;
    Matrix<float, 2, 2>  _S4068 = _S4000 + _S4067;
    Matrix<float, 2, 2>  _S4069 = _S3994 * _S4068;
    float _S4070 = _S3996 * _S4069.rows[int(1)].y;
    float _S4071 = _S4069.rows[int(0)].y + _S4069.rows[int(1)].x;
    float _S4072 = _S3995 * _S4069.rows[int(0)].x;
    float2  s_diff_d_T_0 = make_float2 (_S3996 * _S4071 + _S4072 + _S4072, _S4070 + _S4070 + _S3995 * _S4071);
    Matrix<float, 2, 2>  _S4073 = _S3990 * _S4068;
    float _S4074 = _S3992 * _S4073.rows[int(1)].y;
    float _S4075 = _S4073.rows[int(0)].y + _S4073.rows[int(1)].x;
    float _S4076 = _S3991 * _S4073.rows[int(0)].x;
    float2  s_diff_d_T_1 = make_float2 (_S3992 * _S4075 + _S4076 + _S4076, _S4074 + _S4074 + _S3991 * _S4075);
    Matrix<float, 2, 2>  _S4077 = _S3986 * _S4068;
    float _S4078 = _S3988 * _S4077.rows[int(1)].y;
    float _S4079 = _S4077.rows[int(0)].y + _S4077.rows[int(1)].x;
    float _S4080 = _S3987 * _S4077.rows[int(0)].x;
    float2  s_diff_d_T_2 = make_float2 (_S3988 * _S4079 + _S4080 + _S4080, _S4078 + _S4078 + _S3987 * _S4079);
    Matrix<float, 2, 2>  _S4081 = _S3982 * _S4068;
    float _S4082 = _S3984 * _S4081.rows[int(1)].y;
    float _S4083 = _S4081.rows[int(0)].y + _S4081.rows[int(1)].x;
    float _S4084 = _S3983 * _S4081.rows[int(0)].x;
    float2  s_diff_d_T_3 = make_float2 (_S3984 * _S4083 + _S4084 + _S4084, _S4082 + _S4082 + _S3983 * _S4083);
    Matrix<float, 2, 2>  _S4085 = _S3978 * _S4068;
    float _S4086 = _S3980 * _S4085.rows[int(1)].y;
    float _S4087 = _S4085.rows[int(0)].y + _S4085.rows[int(1)].x;
    float _S4088 = _S3979 * _S4085.rows[int(0)].x;
    float2  s_diff_d_T_4 = make_float2 (_S3980 * _S4087 + _S4088 + _S4088, _S4086 + _S4086 + _S3979 * _S4087);
    Matrix<float, 2, 2>  _S4089 = _S3974 * _S4068;
    float _S4090 = _S3976 * _S4089.rows[int(1)].y;
    float _S4091 = _S4089.rows[int(0)].y + _S4089.rows[int(1)].x;
    float _S4092 = _S3975 * _S4089.rows[int(0)].x;
    float2  s_diff_d_T_5 = make_float2 (_S3976 * _S4091 + _S4092 + _S4092, _S4090 + _S4090 + _S3975 * _S4091);
    Matrix<float, 2, 2>  _S4093 = makeMatrix<float, 2, 2> (2.0f) * _S4068;
    float _S4094 = _S3972 * _S4093.rows[int(1)].y;
    float _S4095 = _S4093.rows[int(0)].y + _S4093.rows[int(1)].x;
    float _S4096 = _S3971 * _S4093.rows[int(0)].x;
    float2  s_diff_d_T_6 = make_float2 (_S3972 * _S4095 + _S4096 + _S4096, _S4094 + _S4094 + _S3971 * _S4095);
    float2  _S4097 = - s_diff_d_T_0 + - s_diff_d_T_1 + - s_diff_d_T_2 + - s_diff_d_T_3 + - s_diff_d_T_4 + - s_diff_d_T_5 + - s_diff_d_T_6 + v_mean2d_7;
    float2  _S4098 = s_diff_d_T_0 + _S3969 * _S4097;
    float2  _S4099 = s_diff_d_T_1 + _S3968 * _S4097;
    float2  _S4100 = s_diff_d_T_2 + _S3967 * _S4097;
    float2  _S4101 = s_diff_d_T_3 + _S3966 * _S4097;
    float2  _S4102 = s_diff_d_T_4 + _S3965 * _S4097;
    float2  _S4103 = s_diff_d_T_5 + _S3964 * _S4097;
    FixedArray<float2 , 7>  _S4104;
    _S4104[int(0)] = _S4022;
    _S4104[int(1)] = _S4022;
    _S4104[int(2)] = _S4022;
    _S4104[int(3)] = _S4022;
    _S4104[int(4)] = _S4022;
    _S4104[int(5)] = _S4022;
    _S4104[int(6)] = _S4022;
    _S4104[int(6)] = _S4098;
    _S4104[int(5)] = _S4099;
    _S4104[int(4)] = _S4100;
    _S4104[int(3)] = _S4101;
    _S4104[int(2)] = _S4102;
    _S4104[int(1)] = _S4103;
    _S4104[int(0)] = s_diff_d_T_6;
    _S3963 = _S4104;
    *&((&_S3963[int(6)])->x) = 0.0f;
    float2  _S4105 = make_float2 (_S4104[int(6)].x, 0.0f);
    FixedArray<float2 , 7>  _S4106;
    _S4106[int(0)] = _S4022;
    _S4106[int(1)] = _S4022;
    _S4106[int(2)] = _S4022;
    _S4106[int(3)] = _S4022;
    _S4106[int(4)] = _S4022;
    _S4106[int(5)] = _S4022;
    _S4106[int(6)] = _S4022;
    _S4106[int(6)] = _S4105;
    float2  _S4107 = _S3963[int(1)] + _S4106[int(1)];
    float2  _S4108 = _S3963[int(2)] + _S4106[int(2)];
    float2  _S4109 = _S3963[int(3)] + _S4106[int(3)];
    float2  _S4110 = _S3963[int(4)] + _S4106[int(4)];
    float2  _S4111 = _S3963[int(5)] + _S4106[int(5)];
    float2  _S4112 = _S3963[int(6)] + _S4106[int(6)];
    _S3963[int(0)] = _S3963[int(0)] + _S4106[int(0)];
    _S3963[int(1)] = _S4107;
    _S3963[int(2)] = _S4108;
    _S3963[int(3)] = _S4109;
    _S3963[int(4)] = _S4110;
    _S3963[int(5)] = _S4111;
    _S3963[int(6)] = _S4112;
    *&((&_S3963[int(5)])->x) = 0.0f;
    float2  _S4113 = make_float2 (_S4111.x, 0.0f);
    FixedArray<float2 , 7>  _S4114;
    _S4114[int(0)] = _S4022;
    _S4114[int(1)] = _S4022;
    _S4114[int(2)] = _S4022;
    _S4114[int(3)] = _S4022;
    _S4114[int(4)] = _S4022;
    _S4114[int(5)] = _S4022;
    _S4114[int(6)] = _S4022;
    _S4114[int(5)] = _S4113;
    float2  _S4115 = _S3963[int(1)] + _S4114[int(1)];
    float2  _S4116 = _S3963[int(2)] + _S4114[int(2)];
    float2  _S4117 = _S3963[int(3)] + _S4114[int(3)];
    float2  _S4118 = _S3963[int(4)] + _S4114[int(4)];
    float2  _S4119 = _S3963[int(5)] + _S4114[int(5)];
    float2  _S4120 = _S3963[int(6)] + _S4114[int(6)];
    _S3963[int(0)] = _S3963[int(0)] + _S4114[int(0)];
    _S3963[int(1)] = _S4115;
    _S3963[int(2)] = _S4116;
    _S3963[int(3)] = _S4117;
    _S3963[int(4)] = _S4118;
    _S3963[int(5)] = _S4119;
    _S3963[int(6)] = _S4120;
    *&((&_S3963[int(4)])->x) = 0.0f;
    Matrix<float, 3, 3>  _S4121 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4122;
    (&_S4122)->primal_0 = R_15;
    (&_S4122)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4123;
    (&_S4123)->primal_0 = _S3933;
    (&_S4123)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4122, &_S4123, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4124;
    (&_S4124)->primal_0 = R_15;
    (&_S4124)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4125;
    (&_S4125)->primal_0 = _S3929;
    (&_S4125)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4124, &_S4125, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4126;
    (&_S4126)->primal_0 = R_15;
    (&_S4126)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4127;
    (&_S4127)->primal_0 = _S3925;
    (&_S4127)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4126, &_S4127, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4128;
    (&_S4128)->primal_0 = R_15;
    (&_S4128)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4129;
    (&_S4129)->primal_0 = _S3932;
    (&_S4129)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4128, &_S4129, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4130;
    (&_S4130)->primal_0 = R_15;
    (&_S4130)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4131;
    (&_S4131)->primal_0 = _S3928;
    (&_S4131)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4130, &_S4131, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4132;
    (&_S4132)->primal_0 = R_15;
    (&_S4132)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4133;
    (&_S4133)->primal_0 = _S3924;
    (&_S4133)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4132, &_S4133, _S3908);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4134;
    (&_S4134)->primal_0 = R_15;
    (&_S4134)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4135;
    (&_S4135)->primal_0 = _S3945.p_0[0U];
    (&_S4135)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4134, &_S4135, _S3908);
    float3  _S4136 = - _S4123.differential_0 + _S4129.differential_0;
    float3  _S4137 = _S3931 * _S4136;
    float3  _S4138 = _S3916.rows[2U] * _S4136;
    float _S4139 = _S3921 * (_S4138.x + _S4138.y + _S4138.z);
    float3  _S4140 = - _S4125.differential_0 + _S4131.differential_0;
    float3  _S4141 = _S3927 * _S4140;
    float3  _S4142 = _S3916.rows[1U] * _S4140;
    float _S4143 = _S3921 * (_S4142.x + _S4142.y + _S4142.z);
    float3  _S4144 = - _S4127.differential_0 + _S4133.differential_0;
    float3  _S4145 = _S3923 * _S4144;
    float3  _S4146 = _S3916.rows[0U] * _S4144;
    float _S4147 = _S3921 * (_S4146.x + _S4146.y + _S4146.z);
    Matrix<float, 3, 3>  _S4148 = _S4121;
    _S4148[2U] = _S4137;
    _S4148[1U] = _S4141;
    _S4148[0U] = _S4145;
    Matrix<float, 3, 3>  _S4149 = transpose_3(transpose_3(_S4148));
    float _S4150 = 2.0f * - _S4149.rows[int(2)].z;
    float _S4151 = 2.0f * _S4149.rows[int(2)].y;
    float _S4152 = 2.0f * _S4149.rows[int(2)].x;
    float _S4153 = 2.0f * _S4149.rows[int(1)].z;
    float _S4154 = 2.0f * - _S4149.rows[int(1)].y;
    float _S4155 = 2.0f * _S4149.rows[int(1)].x;
    float _S4156 = 2.0f * _S4149.rows[int(0)].z;
    float _S4157 = 2.0f * _S4149.rows[int(0)].y;
    float _S4158 = 2.0f * - _S4149.rows[int(0)].x;
    float _S4159 = - _S4155 + _S4157;
    float _S4160 = _S4152 + - _S4156;
    float _S4161 = - _S4151 + _S4153;
    float _S4162 = _S4151 + _S4153;
    float _S4163 = _S4152 + _S4156;
    float _S4164 = _S4155 + _S4157;
    float _S4165 = _S3914.w * (_S4154 + _S4158);
    float _S4166 = _S3914.z * (_S4150 + _S4158);
    float _S4167 = _S3914.y * (_S4150 + _S4154);
    float _S4168 = _S3914.x * _S4159 + _S3914.z * _S4162 + _S3914.y * _S4163 + _S4165 + _S4165;
    float _S4169 = _S3914.x * _S4160 + _S3914.w * _S4162 + _S3914.y * _S4164 + _S4166 + _S4166;
    float _S4170 = _S3914.x * _S4161 + _S3914.w * _S4163 + _S3914.z * _S4164 + _S4167 + _S4167;
    float _S4171 = _S3914.w * _S4159 + _S3914.z * _S4160 + _S3914.y * _S4161;
    float4  _S4172 = make_float4 (0.0f);
    float4  _S4173 = _S4172;
    *&((&_S4173)->w) = _S4168;
    *&((&_S4173)->z) = _S4169;
    *&((&_S4173)->y) = _S4170;
    *&((&_S4173)->x) = _S4171;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S4174;
    (&_S4174)->primal_0 = quat_15;
    (&_S4174)->differential_0 = _S4172;
    s_bwd_normalize_impl_0(&_S4174, _S4173);
    float3  _S4175 = _S3908;
    *&((&_S4175)->z) = _S4139;
    *&((&_S4175)->y) = _S4143;
    *&((&_S4175)->x) = _S4147;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4176;
    (&_S4176)->primal_0 = scale_15;
    (&_S4176)->differential_0 = _S3908;
    s_bwd_prop_exp_1(&_S4176, _S4175);
    float _S4177 = - (_S4041 / _S3912);
    DiffPair_float_0 _S4178;
    (&_S4178)->primal_0 = _S3909;
    (&_S4178)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4178, _S4177);
    float _S4179 = - _S4178.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4180;
    (&_S4180)->primal_0 = mean_c_15;
    (&_S4180)->differential_0 = _S3908;
    s_bwd_length_impl_0(&_S4180, v_depth_7);
    float3  _S4181 = _S4180.differential_0 + _S4038;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4182;
    (&_S4182)->primal_0 = R_15;
    (&_S4182)->differential_0 = _S4121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4183;
    (&_S4183)->primal_0 = mean_16;
    (&_S4183)->differential_0 = _S3908;
    s_bwd_prop_mul_3(&_S4182, &_S4183, _S4181);
    Matrix<float, 3, 3>  _S4184 = _S4122.differential_0 + _S4124.differential_0 + _S4126.differential_0 + _S4128.differential_0 + _S4130.differential_0 + _S4132.differential_0 + _S4134.differential_0 + _S4182.differential_0;
    float _S4185 = _S4179 + _S4021.differential_0;
    float3  _S4186 = _S4176.differential_0 + _S4020.differential_0;
    *v_mean_7 = *v_mean_7 + (_S4123.differential_0 + _S4129.differential_0 + _S4125.differential_0 + _S4131.differential_0 + _S4127.differential_0 + _S4133.differential_0 + _S4183.differential_0 + _S4019.differential_0);
    *v_quat_7 = *v_quat_7 + _S4174.differential_0;
    *v_scale_7 = *v_scale_7 + _S4186;
    *v_in_opacity_7 = *v_in_opacity_7 + _S4185;
    *v_R_7 = *v_R_7 + _S4184;
    *v_t_7 = *v_t_7 + _S4181;
    return;
}

