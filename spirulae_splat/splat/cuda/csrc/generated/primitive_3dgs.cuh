#pragma once

#include "slang.cuh"

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
                float _S17 = 0.0f * v_0;
                float r2_0 = u_0 * u_0 + v_0 * v_0;
                float s_diff_r2_0 = u_0 + u_0 + (_S17 + _S17);
                float _S18 = dist_coeffs_0[int(2)] + r2_0 * dist_coeffs_0[int(3)];
                float _S19 = dist_coeffs_0[int(1)] + r2_0 * _S18;
                float _S20 = dist_coeffs_0[int(0)] + r2_0 * _S19;
                float radial_0 = 1.0f + r2_0 * _S20;
                float _S21 = 2.0f * dist_coeffs_0[int(4)];
                float _S22 = _S21 * u_0;
                float _S23 = 2.0f * u_0;
                float _S24 = 2.0f * dist_coeffs_0[int(5)];
                float _S25 = _S24 * u_0;
                float _S26 = 2.0f * v_0;
                float2  _S27 = make_float2 (1.0f, 0.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_0 * _S20 + (s_diff_r2_0 * _S19 + (s_diff_r2_0 * _S18 + s_diff_r2_0 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * *mean2d_0 + make_float2 (_S21 * v_0 + 0.0f * _S22 + (s_diff_r2_0 + (_S23 + _S23)) * dist_coeffs_0[int(5)] + s_diff_r2_0 * dist_coeffs_0[int(6)], _S24 * v_0 + 0.0f * _S25 + (s_diff_r2_0 + (_S17 + 0.0f * _S26)) * dist_coeffs_0[int(4)] + s_diff_r2_0 * dist_coeffs_0[int(7)]);
                float _S28 = 0.0f * u_0;
                float s_diff_r2_1 = _S28 + _S28 + (v_0 + v_0);
                float2  _S29 = make_float2 (0.0f, 1.0f) * make_float2 (radial_0) + make_float2 (s_diff_r2_1 * _S20 + (s_diff_r2_1 * _S19 + (s_diff_r2_1 * _S18 + s_diff_r2_1 * dist_coeffs_0[int(3)] * r2_0) * r2_0) * r2_0) * *mean2d_0 + make_float2 (0.0f * _S21 * v_0 + _S22 + (s_diff_r2_1 + (_S28 + 0.0f * _S23)) * dist_coeffs_0[int(5)] + s_diff_r2_1 * dist_coeffs_0[int(6)], 0.0f * _S24 * v_0 + _S25 + (s_diff_r2_1 + (_S26 + _S26)) * dist_coeffs_0[int(4)] + s_diff_r2_1 * dist_coeffs_0[int(7)]);
                Matrix<float, 2, 2>  _S30 = transpose_0(makeMatrix<float, 2, 2> (_S27 + make_float2 (_S27.x * dist_coeffs_0[int(8)] + _S27.y * dist_coeffs_0[int(9)], 0.0f), _S29 + make_float2 (_S29.x * dist_coeffs_0[int(8)] + _S29.y * dist_coeffs_0[int(9)], 0.0f)));
                _S12 = !((F32_min((determinant_0(_S30)), ((F32_min((_S30.rows[int(0)].x), (_S30.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S12)
            {
                break;
            }
            float u_1 = (*mean2d_0).x;
            float v_1 = (*mean2d_0).y;
            float r2_1 = u_1 * u_1 + v_1 * v_1;
            float2  _S31 = *mean2d_0 * make_float2 (1.0f + r2_1 * (dist_coeffs_0[int(0)] + r2_1 * (dist_coeffs_0[int(1)] + r2_1 * (dist_coeffs_0[int(2)] + r2_1 * dist_coeffs_0[int(3)])))) + make_float2 (2.0f * dist_coeffs_0[int(4)] * u_1 * v_1 + dist_coeffs_0[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + dist_coeffs_0[int(6)] * r2_1, 2.0f * dist_coeffs_0[int(5)] * u_1 * v_1 + dist_coeffs_0[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + dist_coeffs_0[int(7)] * r2_1);
            float2  _S32 = _S31 + make_float2 (dist_coeffs_0[int(8)] * _S31.x + dist_coeffs_0[int(9)] * _S31.y, 0.0f);
            *mean2d_0 = make_float2 (intrins_0.x * _S32.x + cx_0, intrins_0.y * _S32.y + cy_0);
            break;
        }
        if(!!_S12)
        {
            _S12 = false;
            break;
        }
        Matrix<float, 2, 3>  J_0;
        float2  _S33 = _S13 / make_float2 (_S14);
        float2  _S34 = _S13 * make_float2 (0.0f);
        float _S35 = _S14 * _S14;
        float2  _S36 = (make_float2 (1.0f, 0.0f) * make_float2 (_S14) - _S34) / make_float2 (_S35);
        float u_2 = _S33.x;
        float s_diff_u_0 = _S36.x;
        float v_2 = _S33.y;
        float s_diff_v_0 = _S36.y;
        float _S37 = s_diff_u_0 * u_2;
        float _S38 = s_diff_v_0 * v_2;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float s_diff_r2_2 = _S37 + _S37 + (_S38 + _S38);
        float _S39 = dist_coeffs_0[int(2)] + r2_2 * dist_coeffs_0[int(3)];
        float _S40 = dist_coeffs_0[int(1)] + r2_2 * _S39;
        float _S41 = dist_coeffs_0[int(0)] + r2_2 * _S40;
        float _S42 = 2.0f * dist_coeffs_0[int(4)];
        float _S43 = 2.0f * dist_coeffs_0[int(5)];
        float2  _S44 = _S36 * make_float2 (1.0f + r2_2 * _S41) + make_float2 (s_diff_r2_2 * _S41 + (s_diff_r2_2 * _S40 + (s_diff_r2_2 * _S39 + s_diff_r2_2 * dist_coeffs_0[int(3)] * r2_2) * r2_2) * r2_2) * _S33 + make_float2 (s_diff_u_0 * _S42 * v_2 + s_diff_v_0 * (_S42 * u_2) + (s_diff_r2_2 + (s_diff_u_0 * 2.0f * u_2 + s_diff_u_0 * (2.0f * u_2))) * dist_coeffs_0[int(5)] + s_diff_r2_2 * dist_coeffs_0[int(6)], s_diff_u_0 * _S43 * v_2 + s_diff_v_0 * (_S43 * u_2) + (s_diff_r2_2 + (s_diff_v_0 * 2.0f * v_2 + s_diff_v_0 * (2.0f * v_2))) * dist_coeffs_0[int(4)] + s_diff_r2_2 * dist_coeffs_0[int(7)]);
        float2  _S45 = _S44 + make_float2 (_S44.x * dist_coeffs_0[int(8)] + _S44.y * dist_coeffs_0[int(9)], 0.0f);
        float fx_0 = intrins_0.x;
        float fy_0 = intrins_0.y;
        float _S46 = _S45.y * fy_0;
        Matrix<float, 2, 3>  J_1;
        *&(((&J_1)->rows + (int(0)))->x) = _S45.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->x) = _S46;
        float2  _S47 = _S13 / make_float2 (_S14);
        float2  _S48 = (make_float2 (0.0f, 1.0f) * make_float2 (_S14) - _S34) / make_float2 (_S35);
        float u_3 = _S47.x;
        float s_diff_u_1 = _S48.x;
        float v_3 = _S47.y;
        float s_diff_v_1 = _S48.y;
        float _S49 = s_diff_u_1 * u_3;
        float _S50 = s_diff_v_1 * v_3;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float s_diff_r2_3 = _S49 + _S49 + (_S50 + _S50);
        float _S51 = dist_coeffs_0[int(2)] + r2_3 * dist_coeffs_0[int(3)];
        float _S52 = dist_coeffs_0[int(1)] + r2_3 * _S51;
        float _S53 = dist_coeffs_0[int(0)] + r2_3 * _S52;
        float2  _S54 = _S48 * make_float2 (1.0f + r2_3 * _S53) + make_float2 (s_diff_r2_3 * _S53 + (s_diff_r2_3 * _S52 + (s_diff_r2_3 * _S51 + s_diff_r2_3 * dist_coeffs_0[int(3)] * r2_3) * r2_3) * r2_3) * _S47 + make_float2 (s_diff_u_1 * _S42 * v_3 + s_diff_v_1 * (_S42 * u_3) + (s_diff_r2_3 + (s_diff_u_1 * 2.0f * u_3 + s_diff_u_1 * (2.0f * u_3))) * dist_coeffs_0[int(5)] + s_diff_r2_3 * dist_coeffs_0[int(6)], s_diff_u_1 * _S43 * v_3 + s_diff_v_1 * (_S43 * u_3) + (s_diff_r2_3 + (s_diff_v_1 * 2.0f * v_3 + s_diff_v_1 * (2.0f * v_3))) * dist_coeffs_0[int(4)] + s_diff_r2_3 * dist_coeffs_0[int(7)]);
        float2  _S55 = _S54 + make_float2 (_S54.x * dist_coeffs_0[int(8)] + _S54.y * dist_coeffs_0[int(9)], 0.0f);
        float _S56 = _S55.y * fy_0;
        *&(((&J_1)->rows + (int(0)))->y) = _S55.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->y) = _S56;
        float2  _S57 = _S13 / make_float2 (_S14);
        float2  _S58 = (make_float2 (0.0f, 0.0f) * make_float2 (_S14) - _S13) / make_float2 (_S35);
        float u_4 = _S57.x;
        float s_diff_u_2 = _S58.x;
        float v_4 = _S57.y;
        float s_diff_v_2 = _S58.y;
        float _S59 = s_diff_u_2 * u_4;
        float _S60 = s_diff_v_2 * v_4;
        float r2_4 = u_4 * u_4 + v_4 * v_4;
        float s_diff_r2_4 = _S59 + _S59 + (_S60 + _S60);
        float _S61 = dist_coeffs_0[int(2)] + r2_4 * dist_coeffs_0[int(3)];
        float _S62 = dist_coeffs_0[int(1)] + r2_4 * _S61;
        float _S63 = dist_coeffs_0[int(0)] + r2_4 * _S62;
        float2  _S64 = _S58 * make_float2 (1.0f + r2_4 * _S63) + make_float2 (s_diff_r2_4 * _S63 + (s_diff_r2_4 * _S62 + (s_diff_r2_4 * _S61 + s_diff_r2_4 * dist_coeffs_0[int(3)] * r2_4) * r2_4) * r2_4) * _S57 + make_float2 (s_diff_u_2 * _S42 * v_4 + s_diff_v_2 * (_S42 * u_4) + (s_diff_r2_4 + (s_diff_u_2 * 2.0f * u_4 + s_diff_u_2 * (2.0f * u_4))) * dist_coeffs_0[int(5)] + s_diff_r2_4 * dist_coeffs_0[int(6)], s_diff_u_2 * _S43 * v_4 + s_diff_v_2 * (_S43 * u_4) + (s_diff_r2_4 + (s_diff_v_2 * 2.0f * v_4 + s_diff_v_2 * (2.0f * v_4))) * dist_coeffs_0[int(4)] + s_diff_r2_4 * dist_coeffs_0[int(7)]);
        float2  _S65 = _S64 + make_float2 (_S64.x * dist_coeffs_0[int(8)] + _S64.y * dist_coeffs_0[int(9)], 0.0f);
        float _S66 = _S65.y * fy_0;
        *&(((&J_1)->rows + (int(0)))->z) = _S65.x * fx_0;
        *&(((&J_1)->rows + (int(1)))->z) = _S66;
        J_0 = J_1;
        float _S67 = float(width_0);
        float _S68 = 0.30000001192092896f * (0.5f * _S67);
        float lim_x_pos_0 = _S67 + _S68;
        float rz_0 = 1.0f / _S14;
        float _S69 = - _S68;
        float max_Jyz_0 = - (_S69 - cy_0) * rz_0;
        float min_Jyz_0 = - (lim_x_pos_0 - cy_0) * rz_0;
        *&(((&J_0)->rows + (int(0)))->z) = clamp_0(*&(((&J_0)->rows + (int(0)))->z), - (lim_x_pos_0 - cx_0) * rz_0, - (_S69 - cx_0) * rz_0);
        *&(((&J_0)->rows + (int(1)))->z) = clamp_0(*&(((&J_0)->rows + (int(1)))->z), min_Jyz_0, max_Jyz_0);
        *cov2d_0 = mul_4(mul_3(J_0, cov3d_0), transpose_1(J_0));
        _S12 = true;
        break;
    }
    return _S12;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_3, float dOut_6)
{
    float _S70 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_3).primal_0)))))) * dOut_6;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = _S70;
    return;
}

inline __device__ DiffPair_float_0 _d_sqrt_1(DiffPair_float_0 * dpx_4)
{
    DiffPair_float_0 _S71 = { (F32_sqrt((dpx_4->primal_0))), 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), (dpx_4->primal_0)))))) * dpx_4->differential_0 };
    return _S71;
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

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_2, DiffPair_float_0 * dpx_5, float dOut_7)
{
    DiffPair_float_0 _S72 = *dpx_5;
    float _S73 = - (*dpy_2).primal_0 / ((*dpx_5).primal_0 * (*dpx_5).primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_7;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S73;
    float _S74 = _S72.primal_0 / (_S72.primal_0 * _S72.primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_7;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = _S74;
    return;
}

inline __device__ DiffPair_float_0 _d_atan2_1(DiffPair_float_0 * dpy_3, DiffPair_float_0 * dpx_6)
{
    float _S75 = dpx_6->primal_0 * dpx_6->primal_0 + dpy_3->primal_0 * dpy_3->primal_0;
    DiffPair_float_0 _S76 = { (F32_atan2((dpy_3->primal_0), (dpx_6->primal_0))), - dpy_3->primal_0 / _S75 * dpx_6->differential_0 + dpx_6->primal_0 / _S75 * dpy_3->differential_0 };
    return _S76;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ DiffPair_float_0 s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_7)
{
    float _S77 = *&((&dpx_7->differential_0)->x) * *&((&dpx_7->primal_0)->x);
    float _S78 = *&((&dpx_7->differential_0)->y) * *&((&dpx_7->primal_0)->y);
    float s_diff_len_0 = _S77 + _S77 + (_S78 + _S78);
    DiffPair_float_0 _S79;
    (&_S79)->primal_0 = *&((&dpx_7->primal_0)->x) * *&((&dpx_7->primal_0)->x) + *&((&dpx_7->primal_0)->y) * *&((&dpx_7->primal_0)->y);
    (&_S79)->differential_0 = s_diff_len_0;
    DiffPair_float_0 _S80 = _d_sqrt_1(&_S79);
    DiffPair_float_0 _S81 = { _S80.primal_0, _S80.differential_0 };
    return _S81;
}

inline __device__ bool fisheye_proj_3dgs_nav(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_1, FixedArray<float, 10>  dist_coeffs_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    bool _S82;
    float2  _S83;
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
    float _S94;
    float _S95;
    float _S96;
    bool _S97;
    for(;;)
    {
        float k_0;
        for(;;)
        {
            float2  _S98 = float2 {mean3d_1.x, mean3d_1.y};
            _S83 = _S98;
            float r_7 = length_0(_S98);
            float _S99 = mean3d_1.z;
            _S84 = _S99;
            float theta_0 = (F32_atan2((r_7), (_S99)));
            if(theta_0 < 0.00100000004749745f)
            {
                k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S99;
            }
            else
            {
                k_0 = theta_0 / r_7;
            }
            float2  _S100 = _S98 * make_float2 (k_0);
            *mean2d_1 = _S100;
            float2  _S101 = make_float2 (1.0f, 0.0f);
            _S85 = dist_coeffs_1[int(0)];
            _S86 = dist_coeffs_1[int(1)];
            _S87 = dist_coeffs_1[int(2)];
            _S88 = dist_coeffs_1[int(3)];
            _S89 = dist_coeffs_1[int(4)];
            _S90 = dist_coeffs_1[int(5)];
            _S91 = dist_coeffs_1[int(6)];
            _S92 = dist_coeffs_1[int(7)];
            _S93 = dist_coeffs_1[int(8)];
            _S94 = dist_coeffs_1[int(9)];
            float u_5 = _S100.x;
            float v_5 = _S100.y;
            float _S102 = 0.0f * v_5;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float s_diff_r2_5 = u_5 + u_5 + (_S102 + _S102);
            float _S103 = dist_coeffs_1[int(2)] + r2_5 * dist_coeffs_1[int(3)];
            float _S104 = dist_coeffs_1[int(1)] + r2_5 * _S103;
            float _S105 = dist_coeffs_1[int(0)] + r2_5 * _S104;
            float _S106 = s_diff_r2_5 * _S105 + (s_diff_r2_5 * _S104 + (s_diff_r2_5 * _S103 + s_diff_r2_5 * dist_coeffs_1[int(3)] * r2_5) * r2_5) * r2_5;
            float radial_1 = 1.0f + r2_5 * _S105;
            float _S107 = 2.0f * dist_coeffs_1[int(4)];
            _S95 = _S107;
            float _S108 = _S107 * u_5;
            float _S109 = 2.0f * u_5;
            float s_diff_du_0 = _S107 * v_5 + 0.0f * _S108 + (s_diff_r2_5 + (_S109 + _S109)) * dist_coeffs_1[int(5)] + s_diff_r2_5 * dist_coeffs_1[int(6)];
            float _S110 = 2.0f * dist_coeffs_1[int(5)];
            _S96 = _S110;
            float _S111 = _S110 * u_5;
            float _S112 = 2.0f * v_5;
            float2  _S113 = _S101 * make_float2 (radial_1) + make_float2 (_S106) * _S100 + make_float2 (s_diff_du_0, _S110 * v_5 + 0.0f * _S111 + (s_diff_r2_5 + (_S102 + 0.0f * _S112)) * dist_coeffs_1[int(4)] + s_diff_r2_5 * dist_coeffs_1[int(7)]);
            float _S114 = 0.0f * u_5;
            float s_diff_r2_6 = _S114 + _S114 + (v_5 + v_5);
            float2  _S115 = make_float2 (0.0f, 1.0f) * make_float2 (radial_1) + make_float2 (s_diff_r2_6 * _S105 + (s_diff_r2_6 * _S104 + (s_diff_r2_6 * _S103 + s_diff_r2_6 * dist_coeffs_1[int(3)] * r2_5) * r2_5) * r2_5) * _S100 + make_float2 (0.0f * _S107 * v_5 + _S108 + (s_diff_r2_6 + (_S114 + 0.0f * _S109)) * dist_coeffs_1[int(5)] + s_diff_r2_6 * dist_coeffs_1[int(6)], 0.0f * _S110 * v_5 + _S111 + (s_diff_r2_6 + (_S112 + _S112)) * dist_coeffs_1[int(4)] + s_diff_r2_6 * dist_coeffs_1[int(7)]);
            Matrix<float, 2, 2>  _S116 = transpose_0(makeMatrix<float, 2, 2> (_S113 + make_float2 (_S113.x * dist_coeffs_1[int(8)] + _S113.y * dist_coeffs_1[int(9)], 0.0f), _S115 + make_float2 (_S115.x * dist_coeffs_1[int(8)] + _S115.y * dist_coeffs_1[int(9)], 0.0f)));
            bool _S117 = !((F32_min((determinant_0(_S116)), ((F32_min((_S116.rows[int(0)].x), (_S116.rows[int(1)].y)))))) > 0.0f);
            _S97 = _S117;
            if(_S117)
            {
                break;
            }
            float u_6 = (*mean2d_1).x;
            float v_6 = (*mean2d_1).y;
            float r2_6 = u_6 * u_6 + v_6 * v_6;
            float2  _S118 = *mean2d_1 * make_float2 (1.0f + r2_6 * (dist_coeffs_1[int(0)] + r2_6 * (dist_coeffs_1[int(1)] + r2_6 * (dist_coeffs_1[int(2)] + r2_6 * dist_coeffs_1[int(3)])))) + make_float2 (_S107 * u_6 * v_6 + dist_coeffs_1[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + dist_coeffs_1[int(6)] * r2_6, _S110 * u_6 * v_6 + dist_coeffs_1[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + dist_coeffs_1[int(7)] * r2_6);
            float2  _S119 = _S118 + make_float2 (dist_coeffs_1[int(8)] * _S118.x + dist_coeffs_1[int(9)] * _S118.y, 0.0f);
            *mean2d_1 = make_float2 (intrins_1.x * _S119.x + intrins_1.z, intrins_1.y * _S119.y + intrins_1.w);
            break;
        }
        if(!!_S97)
        {
            _S82 = false;
            break;
        }
        Matrix<float, 2, 3>  J_2;
        float2  _S120 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S121;
        (&_S121)->primal_0 = _S83;
        (&_S121)->differential_0 = _S120;
        DiffPair_float_0 _S122 = s_fwd_length_impl_0(&_S121);
        float _S123 = _S84;
        DiffPair_float_0 _S124;
        (&_S124)->primal_0 = _S122.primal_0;
        (&_S124)->differential_0 = _S122.differential_0;
        DiffPair_float_0 _S125;
        (&_S125)->primal_0 = _S84;
        (&_S125)->differential_0 = 0.0f;
        DiffPair_float_0 _S126 = _d_atan2_1(&_S124, &_S125);
        float s_diff_k_0;
        if((_S126.primal_0) < 0.00100000004749745f)
        {
            float _S127 = _S126.differential_0 * _S126.primal_0;
            float _S128 = 1.0f - _S126.primal_0 * _S126.primal_0 / 3.0f;
            float _S129 = ((0.0f - (_S127 + _S127) * 0.3333333432674408f) * _S84 - _S128 * 0.0f) / (_S84 * _S84);
            k_0 = _S128 / _S84;
            s_diff_k_0 = _S129;
        }
        else
        {
            float _S130 = (_S126.differential_0 * _S122.primal_0 - _S126.primal_0 * _S122.differential_0) / (_S122.primal_0 * _S122.primal_0);
            k_0 = _S126.primal_0 / _S122.primal_0;
            s_diff_k_0 = _S130;
        }
        float2  _S131 = _S83 * make_float2 (k_0);
        float2  _S132 = _S120 * make_float2 (k_0) + make_float2 (s_diff_k_0) * _S83;
        float u_7 = _S131.x;
        float s_diff_u_3 = _S132.x;
        float v_7 = _S131.y;
        float s_diff_v_3 = _S132.y;
        float _S133 = s_diff_u_3 * u_7;
        float _S134 = s_diff_v_3 * v_7;
        float r2_7 = u_7 * u_7 + v_7 * v_7;
        float s_diff_r2_7 = _S133 + _S133 + (_S134 + _S134);
        float _S135 = _S87 + r2_7 * _S88;
        float _S136 = _S86 + r2_7 * _S135;
        float _S137 = _S85 + r2_7 * _S136;
        float2  _S138 = _S132 * make_float2 (1.0f + r2_7 * _S137) + make_float2 (s_diff_r2_7 * _S137 + (s_diff_r2_7 * _S136 + (s_diff_r2_7 * _S135 + s_diff_r2_7 * _S88 * r2_7) * r2_7) * r2_7) * _S131 + make_float2 (s_diff_u_3 * _S95 * v_7 + s_diff_v_3 * (_S95 * u_7) + (s_diff_r2_7 + (s_diff_u_3 * 2.0f * u_7 + s_diff_u_3 * (2.0f * u_7))) * _S90 + s_diff_r2_7 * _S91, s_diff_u_3 * _S96 * v_7 + s_diff_v_3 * (_S96 * u_7) + (s_diff_r2_7 + (s_diff_v_3 * 2.0f * v_7 + s_diff_v_3 * (2.0f * v_7))) * _S89 + s_diff_r2_7 * _S92);
        float2  _S139 = _S138 + make_float2 (_S138.x * _S93 + _S138.y * _S94, 0.0f);
        float fx_1 = intrins_1.x;
        float fy_1 = intrins_1.y;
        float _S140 = _S139.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->x) = _S139.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->x) = _S140;
        float2  _S141 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S142;
        (&_S142)->primal_0 = _S83;
        (&_S142)->differential_0 = _S141;
        DiffPair_float_0 _S143 = s_fwd_length_impl_0(&_S142);
        DiffPair_float_0 _S144;
        (&_S144)->primal_0 = _S143.primal_0;
        (&_S144)->differential_0 = _S143.differential_0;
        DiffPair_float_0 _S145;
        (&_S145)->primal_0 = _S123;
        (&_S145)->differential_0 = 0.0f;
        DiffPair_float_0 _S146 = _d_atan2_1(&_S144, &_S145);
        if((_S146.primal_0) < 0.00100000004749745f)
        {
            float _S147 = _S146.differential_0 * _S146.primal_0;
            float _S148 = 1.0f - _S146.primal_0 * _S146.primal_0 / 3.0f;
            float _S149 = ((0.0f - (_S147 + _S147) * 0.3333333432674408f) * _S84 - _S148 * 0.0f) / (_S84 * _S84);
            k_0 = _S148 / _S84;
            s_diff_k_0 = _S149;
        }
        else
        {
            float _S150 = (_S146.differential_0 * _S143.primal_0 - _S146.primal_0 * _S143.differential_0) / (_S143.primal_0 * _S143.primal_0);
            k_0 = _S146.primal_0 / _S143.primal_0;
            s_diff_k_0 = _S150;
        }
        float2  _S151 = _S83 * make_float2 (k_0);
        float2  _S152 = _S141 * make_float2 (k_0) + make_float2 (s_diff_k_0) * _S83;
        float u_8 = _S151.x;
        float s_diff_u_4 = _S152.x;
        float v_8 = _S151.y;
        float s_diff_v_4 = _S152.y;
        float _S153 = s_diff_u_4 * u_8;
        float _S154 = s_diff_v_4 * v_8;
        float r2_8 = u_8 * u_8 + v_8 * v_8;
        float s_diff_r2_8 = _S153 + _S153 + (_S154 + _S154);
        float _S155 = _S87 + r2_8 * _S88;
        float _S156 = _S86 + r2_8 * _S155;
        float _S157 = _S85 + r2_8 * _S156;
        float2  _S158 = _S152 * make_float2 (1.0f + r2_8 * _S157) + make_float2 (s_diff_r2_8 * _S157 + (s_diff_r2_8 * _S156 + (s_diff_r2_8 * _S155 + s_diff_r2_8 * _S88 * r2_8) * r2_8) * r2_8) * _S151 + make_float2 (s_diff_u_4 * _S95 * v_8 + s_diff_v_4 * (_S95 * u_8) + (s_diff_r2_8 + (s_diff_u_4 * 2.0f * u_8 + s_diff_u_4 * (2.0f * u_8))) * _S90 + s_diff_r2_8 * _S91, s_diff_u_4 * _S96 * v_8 + s_diff_v_4 * (_S96 * u_8) + (s_diff_r2_8 + (s_diff_v_4 * 2.0f * v_8 + s_diff_v_4 * (2.0f * v_8))) * _S89 + s_diff_r2_8 * _S92);
        float2  _S159 = _S158 + make_float2 (_S158.x * _S93 + _S158.y * _S94, 0.0f);
        float _S160 = _S159.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->y) = _S159.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->y) = _S160;
        float2  _S161 = make_float2 (0.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S162;
        (&_S162)->primal_0 = _S83;
        (&_S162)->differential_0 = _S161;
        DiffPair_float_0 _S163 = s_fwd_length_impl_0(&_S162);
        DiffPair_float_0 _S164;
        (&_S164)->primal_0 = _S163.primal_0;
        (&_S164)->differential_0 = _S163.differential_0;
        DiffPair_float_0 _S165;
        (&_S165)->primal_0 = _S84;
        (&_S165)->differential_0 = 1.0f;
        DiffPair_float_0 _S166 = _d_atan2_1(&_S164, &_S165);
        if((_S166.primal_0) < 0.00100000004749745f)
        {
            float _S167 = _S166.differential_0 * _S166.primal_0;
            float _S168 = 1.0f - _S166.primal_0 * _S166.primal_0 / 3.0f;
            float _S169 = ((0.0f - (_S167 + _S167) * 0.3333333432674408f) * _S84 - _S168) / (_S84 * _S84);
            k_0 = _S168 / _S84;
            s_diff_k_0 = _S169;
        }
        else
        {
            float _S170 = (_S166.differential_0 * _S163.primal_0 - _S166.primal_0 * _S163.differential_0) / (_S163.primal_0 * _S163.primal_0);
            k_0 = _S166.primal_0 / _S163.primal_0;
            s_diff_k_0 = _S170;
        }
        float2  _S171 = _S83 * make_float2 (k_0);
        float2  _S172 = _S161 * make_float2 (k_0) + make_float2 (s_diff_k_0) * _S83;
        float u_9 = _S171.x;
        float s_diff_u_5 = _S172.x;
        float v_9 = _S171.y;
        float s_diff_v_5 = _S172.y;
        float _S173 = s_diff_u_5 * u_9;
        float _S174 = s_diff_v_5 * v_9;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float s_diff_r2_9 = _S173 + _S173 + (_S174 + _S174);
        float _S175 = _S87 + r2_9 * _S88;
        float _S176 = _S86 + r2_9 * _S175;
        float _S177 = _S85 + r2_9 * _S176;
        float2  _S178 = _S172 * make_float2 (1.0f + r2_9 * _S177) + make_float2 (s_diff_r2_9 * _S177 + (s_diff_r2_9 * _S176 + (s_diff_r2_9 * _S175 + s_diff_r2_9 * _S88 * r2_9) * r2_9) * r2_9) * _S171 + make_float2 (s_diff_u_5 * _S95 * v_9 + s_diff_v_5 * (_S95 * u_9) + (s_diff_r2_9 + (s_diff_u_5 * 2.0f * u_9 + s_diff_u_5 * (2.0f * u_9))) * _S90 + s_diff_r2_9 * _S91, s_diff_u_5 * _S96 * v_9 + s_diff_v_5 * (_S96 * u_9) + (s_diff_r2_9 + (s_diff_v_5 * 2.0f * v_9 + s_diff_v_5 * (2.0f * v_9))) * _S89 + s_diff_r2_9 * _S92);
        float2  _S179 = _S178 + make_float2 (_S178.x * _S93 + _S178.y * _S94, 0.0f);
        float _S180 = _S179.y * fy_1;
        *&(((&J_2)->rows + (int(0)))->z) = _S179.x * fx_1;
        *&(((&J_2)->rows + (int(1)))->z) = _S180;
        *cov2d_1 = mul_4(mul_3(J_2, cov3d_1), transpose_1(J_2));
        _S82 = true;
        break;
    }
    return _S82;
}

inline __device__ void _d_cos_0(DiffPair_float_0 * dpx_8, float dOut_8)
{
    float _S181 = - (F32_sin(((*dpx_8).primal_0))) * dOut_8;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S181;
    return;
}

inline __device__ void _d_sin_0(DiffPair_float_0 * dpx_9, float dOut_9)
{
    float _S182 = (F32_cos(((*dpx_9).primal_0))) * dOut_9;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S182;
    return;
}

inline __device__ DiffPair_float_0 _d_sin_1(DiffPair_float_0 * dpx_10)
{
    DiffPair_float_0 _S183 = { (F32_sin((dpx_10->primal_0))), (F32_cos((dpx_10->primal_0))) * dpx_10->differential_0 };
    return _S183;
}

inline __device__ bool equisolid_proj_3dgs_nav(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_2, FixedArray<float, 10>  dist_coeffs_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    bool _S184;
    float2  _S185;
    float _S186;
    float _S187;
    float _S188;
    float _S189;
    float _S190;
    float _S191;
    float _S192;
    float _S193;
    float _S194;
    float _S195;
    float _S196;
    float _S197;
    float _S198;
    bool _S199;
    for(;;)
    {
        float k_1;
        for(;;)
        {
            float2  _S200 = float2 {mean3d_2.x, mean3d_2.y};
            _S185 = _S200;
            float r_8 = length_0(_S200);
            float _S201 = mean3d_2.z;
            _S186 = _S201;
            float theta_1 = (F32_atan2((r_8), (_S201)));
            if(r_8 < 9.99999997475242708e-07f)
            {
                k_1 = (1.0f - theta_1 * theta_1 / 24.0f) / _S201;
            }
            else
            {
                k_1 = 2.0f * (F32_sin((0.5f * theta_1))) / r_8;
            }
            float2  _S202 = _S200 * make_float2 (k_1);
            *mean2d_2 = _S202;
            float2  _S203 = make_float2 (1.0f, 0.0f);
            _S187 = dist_coeffs_2[int(0)];
            _S188 = dist_coeffs_2[int(1)];
            _S189 = dist_coeffs_2[int(2)];
            _S190 = dist_coeffs_2[int(3)];
            _S191 = dist_coeffs_2[int(4)];
            _S192 = dist_coeffs_2[int(5)];
            _S193 = dist_coeffs_2[int(6)];
            _S194 = dist_coeffs_2[int(7)];
            _S195 = dist_coeffs_2[int(8)];
            _S196 = dist_coeffs_2[int(9)];
            float u_10 = _S202.x;
            float v_10 = _S202.y;
            float _S204 = 0.0f * v_10;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float s_diff_r2_10 = u_10 + u_10 + (_S204 + _S204);
            float _S205 = dist_coeffs_2[int(2)] + r2_10 * dist_coeffs_2[int(3)];
            float _S206 = dist_coeffs_2[int(1)] + r2_10 * _S205;
            float _S207 = dist_coeffs_2[int(0)] + r2_10 * _S206;
            float _S208 = s_diff_r2_10 * _S207 + (s_diff_r2_10 * _S206 + (s_diff_r2_10 * _S205 + s_diff_r2_10 * dist_coeffs_2[int(3)] * r2_10) * r2_10) * r2_10;
            float radial_2 = 1.0f + r2_10 * _S207;
            float _S209 = 2.0f * dist_coeffs_2[int(4)];
            _S197 = _S209;
            float _S210 = _S209 * u_10;
            float _S211 = 2.0f * u_10;
            float s_diff_du_1 = _S209 * v_10 + 0.0f * _S210 + (s_diff_r2_10 + (_S211 + _S211)) * dist_coeffs_2[int(5)] + s_diff_r2_10 * dist_coeffs_2[int(6)];
            float _S212 = 2.0f * dist_coeffs_2[int(5)];
            _S198 = _S212;
            float _S213 = _S212 * u_10;
            float _S214 = 2.0f * v_10;
            float2  _S215 = _S203 * make_float2 (radial_2) + make_float2 (_S208) * _S202 + make_float2 (s_diff_du_1, _S212 * v_10 + 0.0f * _S213 + (s_diff_r2_10 + (_S204 + 0.0f * _S214)) * dist_coeffs_2[int(4)] + s_diff_r2_10 * dist_coeffs_2[int(7)]);
            float _S216 = 0.0f * u_10;
            float s_diff_r2_11 = _S216 + _S216 + (v_10 + v_10);
            float2  _S217 = make_float2 (0.0f, 1.0f) * make_float2 (radial_2) + make_float2 (s_diff_r2_11 * _S207 + (s_diff_r2_11 * _S206 + (s_diff_r2_11 * _S205 + s_diff_r2_11 * dist_coeffs_2[int(3)] * r2_10) * r2_10) * r2_10) * _S202 + make_float2 (0.0f * _S209 * v_10 + _S210 + (s_diff_r2_11 + (_S216 + 0.0f * _S211)) * dist_coeffs_2[int(5)] + s_diff_r2_11 * dist_coeffs_2[int(6)], 0.0f * _S212 * v_10 + _S213 + (s_diff_r2_11 + (_S214 + _S214)) * dist_coeffs_2[int(4)] + s_diff_r2_11 * dist_coeffs_2[int(7)]);
            Matrix<float, 2, 2>  _S218 = transpose_0(makeMatrix<float, 2, 2> (_S215 + make_float2 (_S215.x * dist_coeffs_2[int(8)] + _S215.y * dist_coeffs_2[int(9)], 0.0f), _S217 + make_float2 (_S217.x * dist_coeffs_2[int(8)] + _S217.y * dist_coeffs_2[int(9)], 0.0f)));
            bool _S219 = !((F32_min((determinant_0(_S218)), ((F32_min((_S218.rows[int(0)].x), (_S218.rows[int(1)].y)))))) > 0.0f);
            _S199 = _S219;
            if(_S219)
            {
                break;
            }
            float u_11 = (*mean2d_2).x;
            float v_11 = (*mean2d_2).y;
            float r2_11 = u_11 * u_11 + v_11 * v_11;
            float2  _S220 = *mean2d_2 * make_float2 (1.0f + r2_11 * (dist_coeffs_2[int(0)] + r2_11 * (dist_coeffs_2[int(1)] + r2_11 * (dist_coeffs_2[int(2)] + r2_11 * dist_coeffs_2[int(3)])))) + make_float2 (_S209 * u_11 * v_11 + dist_coeffs_2[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + dist_coeffs_2[int(6)] * r2_11, _S212 * u_11 * v_11 + dist_coeffs_2[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + dist_coeffs_2[int(7)] * r2_11);
            float2  _S221 = _S220 + make_float2 (dist_coeffs_2[int(8)] * _S220.x + dist_coeffs_2[int(9)] * _S220.y, 0.0f);
            *mean2d_2 = make_float2 (intrins_2.x * _S221.x + intrins_2.z, intrins_2.y * _S221.y + intrins_2.w);
            break;
        }
        if(!!_S199)
        {
            _S184 = false;
            break;
        }
        Matrix<float, 2, 3>  J_3;
        float2  _S222 = make_float2 (1.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S223;
        (&_S223)->primal_0 = _S185;
        (&_S223)->differential_0 = _S222;
        DiffPair_float_0 _S224 = s_fwd_length_impl_0(&_S223);
        float _S225 = _S186;
        DiffPair_float_0 _S226;
        (&_S226)->primal_0 = _S224.primal_0;
        (&_S226)->differential_0 = _S224.differential_0;
        DiffPair_float_0 _S227;
        (&_S227)->primal_0 = _S186;
        (&_S227)->differential_0 = 0.0f;
        DiffPair_float_0 _S228 = _d_atan2_1(&_S226, &_S227);
        float s_diff_k_1;
        if((_S224.primal_0) < 9.99999997475242708e-07f)
        {
            float _S229 = _S228.differential_0 * _S228.primal_0;
            float _S230 = 1.0f - _S228.primal_0 * _S228.primal_0 / 24.0f;
            float _S231 = ((0.0f - (_S229 + _S229) * 0.0416666679084301f) * _S186 - _S230 * 0.0f) / (_S186 * _S186);
            k_1 = _S230 / _S186;
            s_diff_k_1 = _S231;
        }
        else
        {
            float _S232 = _S228.differential_0 * 0.5f;
            DiffPair_float_0 _S233;
            (&_S233)->primal_0 = 0.5f * _S228.primal_0;
            (&_S233)->differential_0 = _S232;
            DiffPair_float_0 _S234 = _d_sin_1(&_S233);
            float _S235 = 2.0f * _S234.primal_0;
            float _S236 = (_S234.differential_0 * 2.0f * _S224.primal_0 - _S235 * _S224.differential_0) / (_S224.primal_0 * _S224.primal_0);
            k_1 = _S235 / _S224.primal_0;
            s_diff_k_1 = _S236;
        }
        float2  _S237 = _S185 * make_float2 (k_1);
        float2  _S238 = _S222 * make_float2 (k_1) + make_float2 (s_diff_k_1) * _S185;
        float u_12 = _S237.x;
        float s_diff_u_6 = _S238.x;
        float v_12 = _S237.y;
        float s_diff_v_6 = _S238.y;
        float _S239 = s_diff_u_6 * u_12;
        float _S240 = s_diff_v_6 * v_12;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float s_diff_r2_12 = _S239 + _S239 + (_S240 + _S240);
        float _S241 = _S189 + r2_12 * _S190;
        float _S242 = _S188 + r2_12 * _S241;
        float _S243 = _S187 + r2_12 * _S242;
        float2  _S244 = _S238 * make_float2 (1.0f + r2_12 * _S243) + make_float2 (s_diff_r2_12 * _S243 + (s_diff_r2_12 * _S242 + (s_diff_r2_12 * _S241 + s_diff_r2_12 * _S190 * r2_12) * r2_12) * r2_12) * _S237 + make_float2 (s_diff_u_6 * _S197 * v_12 + s_diff_v_6 * (_S197 * u_12) + (s_diff_r2_12 + (s_diff_u_6 * 2.0f * u_12 + s_diff_u_6 * (2.0f * u_12))) * _S192 + s_diff_r2_12 * _S193, s_diff_u_6 * _S198 * v_12 + s_diff_v_6 * (_S198 * u_12) + (s_diff_r2_12 + (s_diff_v_6 * 2.0f * v_12 + s_diff_v_6 * (2.0f * v_12))) * _S191 + s_diff_r2_12 * _S194);
        float2  _S245 = _S244 + make_float2 (_S244.x * _S195 + _S244.y * _S196, 0.0f);
        float fx_2 = intrins_2.x;
        float fy_2 = intrins_2.y;
        float _S246 = _S245.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->x) = _S245.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->x) = _S246;
        float2  _S247 = make_float2 (0.0f, 1.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S248;
        (&_S248)->primal_0 = _S185;
        (&_S248)->differential_0 = _S247;
        DiffPair_float_0 _S249 = s_fwd_length_impl_0(&_S248);
        DiffPair_float_0 _S250;
        (&_S250)->primal_0 = _S249.primal_0;
        (&_S250)->differential_0 = _S249.differential_0;
        DiffPair_float_0 _S251;
        (&_S251)->primal_0 = _S225;
        (&_S251)->differential_0 = 0.0f;
        DiffPair_float_0 _S252 = _d_atan2_1(&_S250, &_S251);
        if((_S249.primal_0) < 9.99999997475242708e-07f)
        {
            float _S253 = _S252.differential_0 * _S252.primal_0;
            float _S254 = 1.0f - _S252.primal_0 * _S252.primal_0 / 24.0f;
            float _S255 = ((0.0f - (_S253 + _S253) * 0.0416666679084301f) * _S186 - _S254 * 0.0f) / (_S186 * _S186);
            k_1 = _S254 / _S186;
            s_diff_k_1 = _S255;
        }
        else
        {
            float _S256 = _S252.differential_0 * 0.5f;
            DiffPair_float_0 _S257;
            (&_S257)->primal_0 = 0.5f * _S252.primal_0;
            (&_S257)->differential_0 = _S256;
            DiffPair_float_0 _S258 = _d_sin_1(&_S257);
            float _S259 = 2.0f * _S258.primal_0;
            float _S260 = (_S258.differential_0 * 2.0f * _S249.primal_0 - _S259 * _S249.differential_0) / (_S249.primal_0 * _S249.primal_0);
            k_1 = _S259 / _S249.primal_0;
            s_diff_k_1 = _S260;
        }
        float2  _S261 = _S185 * make_float2 (k_1);
        float2  _S262 = _S247 * make_float2 (k_1) + make_float2 (s_diff_k_1) * _S185;
        float u_13 = _S261.x;
        float s_diff_u_7 = _S262.x;
        float v_13 = _S261.y;
        float s_diff_v_7 = _S262.y;
        float _S263 = s_diff_u_7 * u_13;
        float _S264 = s_diff_v_7 * v_13;
        float r2_13 = u_13 * u_13 + v_13 * v_13;
        float s_diff_r2_13 = _S263 + _S263 + (_S264 + _S264);
        float _S265 = _S189 + r2_13 * _S190;
        float _S266 = _S188 + r2_13 * _S265;
        float _S267 = _S187 + r2_13 * _S266;
        float2  _S268 = _S262 * make_float2 (1.0f + r2_13 * _S267) + make_float2 (s_diff_r2_13 * _S267 + (s_diff_r2_13 * _S266 + (s_diff_r2_13 * _S265 + s_diff_r2_13 * _S190 * r2_13) * r2_13) * r2_13) * _S261 + make_float2 (s_diff_u_7 * _S197 * v_13 + s_diff_v_7 * (_S197 * u_13) + (s_diff_r2_13 + (s_diff_u_7 * 2.0f * u_13 + s_diff_u_7 * (2.0f * u_13))) * _S192 + s_diff_r2_13 * _S193, s_diff_u_7 * _S198 * v_13 + s_diff_v_7 * (_S198 * u_13) + (s_diff_r2_13 + (s_diff_v_7 * 2.0f * v_13 + s_diff_v_7 * (2.0f * v_13))) * _S191 + s_diff_r2_13 * _S194);
        float2  _S269 = _S268 + make_float2 (_S268.x * _S195 + _S268.y * _S196, 0.0f);
        float _S270 = _S269.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->y) = _S269.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->y) = _S270;
        float2  _S271 = make_float2 (0.0f, 0.0f);
        DiffPair_vectorx3Cfloatx2C2x3E_0 _S272;
        (&_S272)->primal_0 = _S185;
        (&_S272)->differential_0 = _S271;
        DiffPair_float_0 _S273 = s_fwd_length_impl_0(&_S272);
        DiffPair_float_0 _S274;
        (&_S274)->primal_0 = _S273.primal_0;
        (&_S274)->differential_0 = _S273.differential_0;
        DiffPair_float_0 _S275;
        (&_S275)->primal_0 = _S186;
        (&_S275)->differential_0 = 1.0f;
        DiffPair_float_0 _S276 = _d_atan2_1(&_S274, &_S275);
        if((_S273.primal_0) < 9.99999997475242708e-07f)
        {
            float _S277 = _S276.differential_0 * _S276.primal_0;
            float _S278 = 1.0f - _S276.primal_0 * _S276.primal_0 / 24.0f;
            float _S279 = ((0.0f - (_S277 + _S277) * 0.0416666679084301f) * _S186 - _S278) / (_S186 * _S186);
            k_1 = _S278 / _S186;
            s_diff_k_1 = _S279;
        }
        else
        {
            float _S280 = _S276.differential_0 * 0.5f;
            DiffPair_float_0 _S281;
            (&_S281)->primal_0 = 0.5f * _S276.primal_0;
            (&_S281)->differential_0 = _S280;
            DiffPair_float_0 _S282 = _d_sin_1(&_S281);
            float _S283 = 2.0f * _S282.primal_0;
            float _S284 = (_S282.differential_0 * 2.0f * _S273.primal_0 - _S283 * _S273.differential_0) / (_S273.primal_0 * _S273.primal_0);
            k_1 = _S283 / _S273.primal_0;
            s_diff_k_1 = _S284;
        }
        float2  _S285 = _S185 * make_float2 (k_1);
        float2  _S286 = _S271 * make_float2 (k_1) + make_float2 (s_diff_k_1) * _S185;
        float u_14 = _S285.x;
        float s_diff_u_8 = _S286.x;
        float v_14 = _S285.y;
        float s_diff_v_8 = _S286.y;
        float _S287 = s_diff_u_8 * u_14;
        float _S288 = s_diff_v_8 * v_14;
        float r2_14 = u_14 * u_14 + v_14 * v_14;
        float s_diff_r2_14 = _S287 + _S287 + (_S288 + _S288);
        float _S289 = _S189 + r2_14 * _S190;
        float _S290 = _S188 + r2_14 * _S289;
        float _S291 = _S187 + r2_14 * _S290;
        float2  _S292 = _S286 * make_float2 (1.0f + r2_14 * _S291) + make_float2 (s_diff_r2_14 * _S291 + (s_diff_r2_14 * _S290 + (s_diff_r2_14 * _S289 + s_diff_r2_14 * _S190 * r2_14) * r2_14) * r2_14) * _S285 + make_float2 (s_diff_u_8 * _S197 * v_14 + s_diff_v_8 * (_S197 * u_14) + (s_diff_r2_14 + (s_diff_u_8 * 2.0f * u_14 + s_diff_u_8 * (2.0f * u_14))) * _S192 + s_diff_r2_14 * _S193, s_diff_u_8 * _S198 * v_14 + s_diff_v_8 * (_S198 * u_14) + (s_diff_r2_14 + (s_diff_v_8 * 2.0f * v_14 + s_diff_v_8 * (2.0f * v_14))) * _S191 + s_diff_r2_14 * _S194);
        float2  _S293 = _S292 + make_float2 (_S292.x * _S195 + _S292.y * _S196, 0.0f);
        float _S294 = _S293.y * fy_2;
        *&(((&J_3)->rows + (int(0)))->z) = _S293.x * fx_2;
        *&(((&J_3)->rows + (int(1)))->z) = _S294;
        *cov2d_2 = mul_4(mul_3(J_3, cov3d_2), transpose_1(J_3));
        _S184 = true;
        break;
    }
    return _S184;
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

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_6, float3  dOut_10)
{
    float _S295 = (*left_6).primal_0.rows[int(0)].x * dOut_10.x;
    Matrix<float, 3, 3>  left_d_result_3;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = (*right_6).primal_0.x * dOut_10.x;
    float sum_6 = _S295 + (*left_6).primal_0.rows[int(1)].x * dOut_10.y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = (*right_6).primal_0.x * dOut_10.y;
    float sum_7 = sum_6 + (*left_6).primal_0.rows[int(2)].x * dOut_10.z;
    *&(((&left_d_result_3)->rows + (int(2)))->x) = (*right_6).primal_0.x * dOut_10.z;
    float3  right_d_result_3;
    *&((&right_d_result_3)->x) = sum_7;
    float _S296 = (*left_6).primal_0.rows[int(0)].y * dOut_10.x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = (*right_6).primal_0.y * dOut_10.x;
    float sum_8 = _S296 + (*left_6).primal_0.rows[int(1)].y * dOut_10.y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = (*right_6).primal_0.y * dOut_10.y;
    float sum_9 = sum_8 + (*left_6).primal_0.rows[int(2)].y * dOut_10.z;
    *&(((&left_d_result_3)->rows + (int(2)))->y) = (*right_6).primal_0.y * dOut_10.z;
    *&((&right_d_result_3)->y) = sum_9;
    float _S297 = (*left_6).primal_0.rows[int(0)].z * dOut_10.x;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = (*right_6).primal_0.z * dOut_10.x;
    float sum_10 = _S297 + (*left_6).primal_0.rows[int(1)].z * dOut_10.y;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = (*right_6).primal_0.z * dOut_10.y;
    float sum_11 = sum_10 + (*left_6).primal_0.rows[int(2)].z * dOut_10.z;
    *&(((&left_d_result_3)->rows + (int(2)))->z) = (*right_6).primal_0.z * dOut_10.z;
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

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_11, float dOut_11)
{
    float _S298 = (F32_exp(((*dpx_11).primal_0))) * dOut_11;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S298;
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

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, float3  dOut_12)
{
    float3  _S299 = exp_0((*dpx_12).primal_0) * dOut_12;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S299;
    return;
}

inline __device__ float4  normalize_0(float4  x_14)
{
    return x_14 / make_float4 (length_2(x_14));
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_13, float dOut_13)
{
    float _S300 = 1.0f / (*dpx_13).primal_0 * dOut_13;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S300;
    return;
}

inline __device__ float view_radius_3dgs_0(float3  mean_0, float3  log_scale_0, float logit_opacity_0, float3  campos_0)
{
    float radius_0 = (F32_exp(((F32_max((log_scale_0.x), ((F32_max((log_scale_0.y), (log_scale_0.z))))))))) * (F32_sqrt((2.0f * (F32_log(((F32_max((255.0f / (1.0f + (F32_exp((- logit_opacity_0))))), (1.0f)))))))));
    float dist_0 = length_1(mean_0 - campos_0);
    return radius_0 / ((F32_max((dist_0), (radius_0))) + (F32_sqrt(((F32_max((dist_0 * dist_0 - radius_0 * radius_0), (0.0f)))))));
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_1, float4  quat_0, float3  scale_0, float in_opacity_0, Matrix<float, 3, 3>  R_0, float3  t_0, float fx_3, float fy_3, float cx_1, float cy_1, FixedArray<float, 10>  dist_coeffs_3, uint image_width_0, uint image_height_0, float4  * aabb_xyxy_0, float * sorting_depth_0, float * radius_1, float2  * mean2d_3, float * depth_0, float3  * conic_0, float * opacity_0)
{
    float2  _S301;
    for(;;)
    {
        float3  mean_c_0 = mul_6(R_0, mean_1) + t_0;
        float _S302 = mean_c_0.z;
        *depth_0 = _S302;
        if(_S302 <= 0.0f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        *opacity_0 = 1.0f / (1.0f + (F32_exp((- in_opacity_0))));
        bool _S303;
        float4  _S304 = normalize_0(quat_0);
        float3  _S305 = exp_0(scale_0);
        float x_15 = _S304.y;
        float x2_0 = x_15 * x_15;
        float y2_0 = _S304.z * _S304.z;
        float z2_0 = _S304.w * _S304.w;
        float xy_0 = _S304.y * _S304.z;
        float xz_0 = _S304.y * _S304.w;
        float yz_0 = _S304.z * _S304.w;
        float wx_0 = _S304.x * _S304.y;
        float wy_0 = _S304.x * _S304.z;
        float wz_0 = _S304.x * _S304.w;
        Matrix<float, 3, 3>  M_0 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0))), makeMatrix<float, 3, 3> (_S305.x, 0.0f, 0.0f, 0.0f, _S305.y, 0.0f, 0.0f, 0.0f, _S305.z));
        Matrix<float, 3, 3>  _S306 = transpose_3(R_0);
        Matrix<float, 3, 3>  covar_c_0 = mul_5(mul_5(R_0, mul_5(M_0, transpose_3(M_0))), _S306);
        for(;;)
        {
            for(;;)
            {
                float2  _S307 = float2 {mean_c_0.x, mean_c_0.y};
                _S301 = _S307;
                *mean2d_3 = _S307 / make_float2 (_S302);
                if(_S302 < 0.0f)
                {
                    _S303 = true;
                }
                else
                {
                    float u_15 = (*mean2d_3).x;
                    float v_15 = (*mean2d_3).y;
                    float _S308 = 0.0f * v_15;
                    float r2_15 = u_15 * u_15 + v_15 * v_15;
                    float s_diff_r2_15 = u_15 + u_15 + (_S308 + _S308);
                    float _S309 = dist_coeffs_3[int(2)] + r2_15 * dist_coeffs_3[int(3)];
                    float _S310 = dist_coeffs_3[int(1)] + r2_15 * _S309;
                    float _S311 = dist_coeffs_3[int(0)] + r2_15 * _S310;
                    float radial_3 = 1.0f + r2_15 * _S311;
                    float _S312 = 2.0f * dist_coeffs_3[int(4)];
                    float _S313 = _S312 * u_15;
                    float _S314 = 2.0f * u_15;
                    float _S315 = 2.0f * dist_coeffs_3[int(5)];
                    float _S316 = _S315 * u_15;
                    float _S317 = 2.0f * v_15;
                    float2  _S318 = make_float2 (1.0f, 0.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_15 * _S311 + (s_diff_r2_15 * _S310 + (s_diff_r2_15 * _S309 + s_diff_r2_15 * dist_coeffs_3[int(3)] * r2_15) * r2_15) * r2_15) * *mean2d_3 + make_float2 (_S312 * v_15 + 0.0f * _S313 + (s_diff_r2_15 + (_S314 + _S314)) * dist_coeffs_3[int(5)] + s_diff_r2_15 * dist_coeffs_3[int(6)], _S315 * v_15 + 0.0f * _S316 + (s_diff_r2_15 + (_S308 + 0.0f * _S317)) * dist_coeffs_3[int(4)] + s_diff_r2_15 * dist_coeffs_3[int(7)]);
                    float _S319 = 0.0f * u_15;
                    float s_diff_r2_16 = _S319 + _S319 + (v_15 + v_15);
                    float2  _S320 = make_float2 (0.0f, 1.0f) * make_float2 (radial_3) + make_float2 (s_diff_r2_16 * _S311 + (s_diff_r2_16 * _S310 + (s_diff_r2_16 * _S309 + s_diff_r2_16 * dist_coeffs_3[int(3)] * r2_15) * r2_15) * r2_15) * *mean2d_3 + make_float2 (0.0f * _S312 * v_15 + _S313 + (s_diff_r2_16 + (_S319 + 0.0f * _S314)) * dist_coeffs_3[int(5)] + s_diff_r2_16 * dist_coeffs_3[int(6)], 0.0f * _S315 * v_15 + _S316 + (s_diff_r2_16 + (_S317 + _S317)) * dist_coeffs_3[int(4)] + s_diff_r2_16 * dist_coeffs_3[int(7)]);
                    Matrix<float, 2, 2>  _S321 = transpose_0(makeMatrix<float, 2, 2> (_S318 + make_float2 (_S318.x * dist_coeffs_3[int(8)] + _S318.y * dist_coeffs_3[int(9)], 0.0f), _S320 + make_float2 (_S320.x * dist_coeffs_3[int(8)] + _S320.y * dist_coeffs_3[int(9)], 0.0f)));
                    _S303 = !((F32_min((determinant_0(_S321)), ((F32_min((_S321.rows[int(0)].x), (_S321.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S303)
                {
                    break;
                }
                float u_16 = (*mean2d_3).x;
                float v_16 = (*mean2d_3).y;
                float r2_16 = u_16 * u_16 + v_16 * v_16;
                float2  _S322 = *mean2d_3 * make_float2 (1.0f + r2_16 * (dist_coeffs_3[int(0)] + r2_16 * (dist_coeffs_3[int(1)] + r2_16 * (dist_coeffs_3[int(2)] + r2_16 * dist_coeffs_3[int(3)])))) + make_float2 (2.0f * dist_coeffs_3[int(4)] * u_16 * v_16 + dist_coeffs_3[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + dist_coeffs_3[int(6)] * r2_16, 2.0f * dist_coeffs_3[int(5)] * u_16 * v_16 + dist_coeffs_3[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + dist_coeffs_3[int(7)] * r2_16);
                float2  _S323 = _S322 + make_float2 (dist_coeffs_3[int(8)] * _S322.x + dist_coeffs_3[int(9)] * _S322.y, 0.0f);
                *mean2d_3 = make_float2 (fx_3 * _S323.x + cx_1, fy_3 * _S323.y + cy_1);
                break;
            }
            if(!!_S303)
            {
                _S303 = false;
                break;
            }
            Matrix<float, 2, 3>  J_4;
            float2  _S324 = _S301 / make_float2 (_S302);
            float2  _S325 = _S301 * make_float2 (0.0f);
            float _S326 = _S302 * _S302;
            float2  _S327 = (make_float2 (1.0f, 0.0f) * make_float2 (_S302) - _S325) / make_float2 (_S326);
            float u_17 = _S324.x;
            float s_diff_u_9 = _S327.x;
            float v_17 = _S324.y;
            float s_diff_v_9 = _S327.y;
            float _S328 = s_diff_u_9 * u_17;
            float _S329 = s_diff_v_9 * v_17;
            float r2_17 = u_17 * u_17 + v_17 * v_17;
            float s_diff_r2_17 = _S328 + _S328 + (_S329 + _S329);
            float _S330 = dist_coeffs_3[int(2)] + r2_17 * dist_coeffs_3[int(3)];
            float _S331 = dist_coeffs_3[int(1)] + r2_17 * _S330;
            float _S332 = dist_coeffs_3[int(0)] + r2_17 * _S331;
            float _S333 = 2.0f * dist_coeffs_3[int(4)];
            float _S334 = 2.0f * dist_coeffs_3[int(5)];
            float2  _S335 = _S327 * make_float2 (1.0f + r2_17 * _S332) + make_float2 (s_diff_r2_17 * _S332 + (s_diff_r2_17 * _S331 + (s_diff_r2_17 * _S330 + s_diff_r2_17 * dist_coeffs_3[int(3)] * r2_17) * r2_17) * r2_17) * _S324 + make_float2 (s_diff_u_9 * _S333 * v_17 + s_diff_v_9 * (_S333 * u_17) + (s_diff_r2_17 + (s_diff_u_9 * 2.0f * u_17 + s_diff_u_9 * (2.0f * u_17))) * dist_coeffs_3[int(5)] + s_diff_r2_17 * dist_coeffs_3[int(6)], s_diff_u_9 * _S334 * v_17 + s_diff_v_9 * (_S334 * u_17) + (s_diff_r2_17 + (s_diff_v_9 * 2.0f * v_17 + s_diff_v_9 * (2.0f * v_17))) * dist_coeffs_3[int(4)] + s_diff_r2_17 * dist_coeffs_3[int(7)]);
            float2  _S336 = _S335 + make_float2 (_S335.x * dist_coeffs_3[int(8)] + _S335.y * dist_coeffs_3[int(9)], 0.0f);
            float _S337 = _S336.y * fy_3;
            Matrix<float, 2, 3>  J_5;
            *&(((&J_5)->rows + (int(0)))->x) = _S336.x * fx_3;
            *&(((&J_5)->rows + (int(1)))->x) = _S337;
            float2  _S338 = _S301 / make_float2 (_S302);
            float2  _S339 = (make_float2 (0.0f, 1.0f) * make_float2 (_S302) - _S325) / make_float2 (_S326);
            float u_18 = _S338.x;
            float s_diff_u_10 = _S339.x;
            float v_18 = _S338.y;
            float s_diff_v_10 = _S339.y;
            float _S340 = s_diff_u_10 * u_18;
            float _S341 = s_diff_v_10 * v_18;
            float r2_18 = u_18 * u_18 + v_18 * v_18;
            float s_diff_r2_18 = _S340 + _S340 + (_S341 + _S341);
            float _S342 = dist_coeffs_3[int(2)] + r2_18 * dist_coeffs_3[int(3)];
            float _S343 = dist_coeffs_3[int(1)] + r2_18 * _S342;
            float _S344 = dist_coeffs_3[int(0)] + r2_18 * _S343;
            float2  _S345 = _S339 * make_float2 (1.0f + r2_18 * _S344) + make_float2 (s_diff_r2_18 * _S344 + (s_diff_r2_18 * _S343 + (s_diff_r2_18 * _S342 + s_diff_r2_18 * dist_coeffs_3[int(3)] * r2_18) * r2_18) * r2_18) * _S338 + make_float2 (s_diff_u_10 * _S333 * v_18 + s_diff_v_10 * (_S333 * u_18) + (s_diff_r2_18 + (s_diff_u_10 * 2.0f * u_18 + s_diff_u_10 * (2.0f * u_18))) * dist_coeffs_3[int(5)] + s_diff_r2_18 * dist_coeffs_3[int(6)], s_diff_u_10 * _S334 * v_18 + s_diff_v_10 * (_S334 * u_18) + (s_diff_r2_18 + (s_diff_v_10 * 2.0f * v_18 + s_diff_v_10 * (2.0f * v_18))) * dist_coeffs_3[int(4)] + s_diff_r2_18 * dist_coeffs_3[int(7)]);
            float2  _S346 = _S345 + make_float2 (_S345.x * dist_coeffs_3[int(8)] + _S345.y * dist_coeffs_3[int(9)], 0.0f);
            float _S347 = _S346.y * fy_3;
            *&(((&J_5)->rows + (int(0)))->y) = _S346.x * fx_3;
            *&(((&J_5)->rows + (int(1)))->y) = _S347;
            float2  _S348 = _S301 / make_float2 (_S302);
            float2  _S349 = (make_float2 (0.0f, 0.0f) * make_float2 (_S302) - _S301) / make_float2 (_S326);
            float u_19 = _S348.x;
            float s_diff_u_11 = _S349.x;
            float v_19 = _S348.y;
            float s_diff_v_11 = _S349.y;
            float _S350 = s_diff_u_11 * u_19;
            float _S351 = s_diff_v_11 * v_19;
            float r2_19 = u_19 * u_19 + v_19 * v_19;
            float s_diff_r2_19 = _S350 + _S350 + (_S351 + _S351);
            float _S352 = dist_coeffs_3[int(2)] + r2_19 * dist_coeffs_3[int(3)];
            float _S353 = dist_coeffs_3[int(1)] + r2_19 * _S352;
            float _S354 = dist_coeffs_3[int(0)] + r2_19 * _S353;
            float2  _S355 = _S349 * make_float2 (1.0f + r2_19 * _S354) + make_float2 (s_diff_r2_19 * _S354 + (s_diff_r2_19 * _S353 + (s_diff_r2_19 * _S352 + s_diff_r2_19 * dist_coeffs_3[int(3)] * r2_19) * r2_19) * r2_19) * _S348 + make_float2 (s_diff_u_11 * _S333 * v_19 + s_diff_v_11 * (_S333 * u_19) + (s_diff_r2_19 + (s_diff_u_11 * 2.0f * u_19 + s_diff_u_11 * (2.0f * u_19))) * dist_coeffs_3[int(5)] + s_diff_r2_19 * dist_coeffs_3[int(6)], s_diff_u_11 * _S334 * v_19 + s_diff_v_11 * (_S334 * u_19) + (s_diff_r2_19 + (s_diff_v_11 * 2.0f * v_19 + s_diff_v_11 * (2.0f * v_19))) * dist_coeffs_3[int(4)] + s_diff_r2_19 * dist_coeffs_3[int(7)]);
            float2  _S356 = _S355 + make_float2 (_S355.x * dist_coeffs_3[int(8)] + _S355.y * dist_coeffs_3[int(9)], 0.0f);
            float _S357 = _S356.y * fy_3;
            *&(((&J_5)->rows + (int(0)))->z) = _S356.x * fx_3;
            *&(((&J_5)->rows + (int(1)))->z) = _S357;
            J_4 = J_5;
            float _S358 = float(image_width_0);
            float _S359 = 0.30000001192092896f * (0.5f * _S358);
            float lim_x_pos_1 = _S358 + _S359;
            float rz_1 = 1.0f / _S302;
            float _S360 = - _S359;
            float max_Jyz_1 = - (_S360 - cy_1) * rz_1;
            float min_Jyz_1 = - (lim_x_pos_1 - cy_1) * rz_1;
            *&(((&J_4)->rows + (int(0)))->z) = clamp_0(*&(((&J_4)->rows + (int(0)))->z), - (lim_x_pos_1 - cx_1) * rz_1, - (_S360 - cx_1) * rz_1);
            *&(((&J_4)->rows + (int(1)))->z) = clamp_0(*&(((&J_4)->rows + (int(1)))->z), min_Jyz_1, max_Jyz_1);
            covar2d_0 = mul_4(mul_3(J_4, covar_c_0), transpose_1(J_4));
            _S303 = true;
            break;
        }
        if(!(true & _S303))
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
        float _S361 = *&(((&covar2d_0)->rows + (int(0)))->x) + eps2d_0;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S361;
        float _S362 = *&(((&covar2d_0)->rows + (int(1)))->y) + eps2d_0;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S362;
        float det_blur_0 = _S361 * _S362 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
        if(det_blur_0 <= 0.0f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float invdet_0 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S363 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_0, - covar2d_0.rows[int(0)].y * invdet_0, - covar2d_0.rows[int(1)].x * invdet_0, covar2d_0.rows[int(0)].x * invdet_0);
        if(antialiased_0)
        {
            *opacity_0 = *opacity_0 * compensation_0;
        }
        if((*opacity_0) < 0.00392156885936856f)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        float _S364 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_0 / 0.00392156885936856f)))))))));
        float radius_x_0 = _S364 * (F32_sqrt((covar2d_0[int(0)].x)));
        float radius_y_0 = _S364 * (F32_sqrt((covar2d_0[int(1)].y)));
        float _S365 = (*mean2d_3).x - radius_x_0;
        float _S366 = (*mean2d_3).x + radius_x_0;
        float _S367 = (*mean2d_3).y - radius_y_0;
        float _S368 = (*mean2d_3).y + radius_y_0;
        if(_S366 <= 0.0f)
        {
            _S303 = true;
        }
        else
        {
            _S303 = _S365 >= float(image_width_0);
        }
        if(_S303)
        {
            _S303 = true;
        }
        else
        {
            _S303 = _S368 <= 0.0f;
        }
        if(_S303)
        {
            _S303 = true;
        }
        else
        {
            _S303 = _S367 >= float(image_height_0);
        }
        if(_S303)
        {
            *aabb_xyxy_0 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_0 = make_float4 (_S365, _S367, _S366, _S368);
        *sorting_depth_0 = _S302;
        *conic_0 = make_float3 (_S363.rows[int(0)].x, _S363.rows[int(0)].y, _S363.rows[int(1)].y);
        *radius_1 = view_radius_3dgs_0(mean_1, scale_0, in_opacity_0, - mul_6(_S306, t_0));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_2, float4  quat_1, float3  scale_1, float in_opacity_1, Matrix<float, 3, 3>  R_1, float3  t_1, float fx_4, float fy_4, float cx_2, float cy_2, FixedArray<float, 10>  dist_coeffs_4, uint image_width_1, uint image_height_1, float4  * aabb_xyxy_1, float * sorting_depth_1, float * radius_2, float2  * mean2d_4, float * depth_1, float3  * conic_1, float * opacity_1)
{
    float2  _S369;
    float _S370;
    float _S371;
    float _S372;
    float _S373;
    float _S374;
    float _S375;
    float _S376;
    float _S377;
    float _S378;
    float _S379;
    float _S380;
    float _S381;
    bool _S382;
    for(;;)
    {
        float3  mean_c_1 = mul_6(R_1, mean_2) + t_1;
        float _S383 = length_1(mean_c_1);
        float _S384 = mean_c_1.z;
        *depth_1 = _S384;
        if(_S383 <= 0.0f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_1;
        *opacity_1 = 1.0f / (1.0f + (F32_exp((- in_opacity_1))));
        bool is_valid_0;
        float eps2d_1;
        float4  _S385 = normalize_0(quat_1);
        float3  _S386 = exp_0(scale_1);
        float x_16 = _S385.y;
        float x2_1 = x_16 * x_16;
        float y2_1 = _S385.z * _S385.z;
        float z2_1 = _S385.w * _S385.w;
        float xy_1 = _S385.y * _S385.z;
        float xz_1 = _S385.y * _S385.w;
        float yz_1 = _S385.z * _S385.w;
        float wx_1 = _S385.x * _S385.y;
        float wy_1 = _S385.x * _S385.z;
        float wz_1 = _S385.x * _S385.w;
        Matrix<float, 3, 3>  M_1 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (_S386.x, 0.0f, 0.0f, 0.0f, _S386.y, 0.0f, 0.0f, 0.0f, _S386.z));
        Matrix<float, 3, 3>  _S387 = transpose_3(R_1);
        Matrix<float, 3, 3>  covar_c_1 = mul_5(mul_5(R_1, mul_5(M_1, transpose_3(M_1))), _S387);
        for(;;)
        {
            float k_2;
            for(;;)
            {
                float2  _S388 = float2 {mean_c_1.x, mean_c_1.y};
                _S369 = _S388;
                float r_9 = length_0(_S388);
                float theta_2 = (F32_atan2((r_9), (_S384)));
                if(theta_2 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_2 * theta_2 / 3.0f) / _S384;
                }
                else
                {
                    k_2 = theta_2 / r_9;
                }
                float2  _S389 = _S388 * make_float2 (k_2);
                *mean2d_4 = _S389;
                float2  _S390 = make_float2 (1.0f, 0.0f);
                _S370 = dist_coeffs_4[int(0)];
                _S371 = dist_coeffs_4[int(1)];
                _S372 = dist_coeffs_4[int(2)];
                _S373 = dist_coeffs_4[int(3)];
                _S374 = dist_coeffs_4[int(4)];
                _S375 = dist_coeffs_4[int(5)];
                _S376 = dist_coeffs_4[int(6)];
                _S377 = dist_coeffs_4[int(7)];
                _S378 = dist_coeffs_4[int(8)];
                _S379 = dist_coeffs_4[int(9)];
                float u_20 = _S389.x;
                float v_20 = _S389.y;
                float _S391 = 0.0f * v_20;
                float r2_20 = u_20 * u_20 + v_20 * v_20;
                float s_diff_r2_20 = u_20 + u_20 + (_S391 + _S391);
                float _S392 = dist_coeffs_4[int(2)] + r2_20 * dist_coeffs_4[int(3)];
                float _S393 = dist_coeffs_4[int(1)] + r2_20 * _S392;
                float _S394 = dist_coeffs_4[int(0)] + r2_20 * _S393;
                float _S395 = s_diff_r2_20 * _S394 + (s_diff_r2_20 * _S393 + (s_diff_r2_20 * _S392 + s_diff_r2_20 * dist_coeffs_4[int(3)] * r2_20) * r2_20) * r2_20;
                float radial_4 = 1.0f + r2_20 * _S394;
                float _S396 = 2.0f * dist_coeffs_4[int(4)];
                _S380 = _S396;
                float _S397 = _S396 * u_20;
                float _S398 = 2.0f * u_20;
                float s_diff_du_2 = _S396 * v_20 + 0.0f * _S397 + (s_diff_r2_20 + (_S398 + _S398)) * dist_coeffs_4[int(5)] + s_diff_r2_20 * dist_coeffs_4[int(6)];
                float _S399 = 2.0f * dist_coeffs_4[int(5)];
                _S381 = _S399;
                float _S400 = _S399 * u_20;
                float _S401 = 2.0f * v_20;
                float2  _S402 = _S390 * make_float2 (radial_4) + make_float2 (_S395) * _S389 + make_float2 (s_diff_du_2, _S399 * v_20 + 0.0f * _S400 + (s_diff_r2_20 + (_S391 + 0.0f * _S401)) * dist_coeffs_4[int(4)] + s_diff_r2_20 * dist_coeffs_4[int(7)]);
                float _S403 = 0.0f * u_20;
                float s_diff_r2_21 = _S403 + _S403 + (v_20 + v_20);
                float2  _S404 = make_float2 (0.0f, 1.0f) * make_float2 (radial_4) + make_float2 (s_diff_r2_21 * _S394 + (s_diff_r2_21 * _S393 + (s_diff_r2_21 * _S392 + s_diff_r2_21 * dist_coeffs_4[int(3)] * r2_20) * r2_20) * r2_20) * _S389 + make_float2 (0.0f * _S396 * v_20 + _S397 + (s_diff_r2_21 + (_S403 + 0.0f * _S398)) * dist_coeffs_4[int(5)] + s_diff_r2_21 * dist_coeffs_4[int(6)], 0.0f * _S399 * v_20 + _S400 + (s_diff_r2_21 + (_S401 + _S401)) * dist_coeffs_4[int(4)] + s_diff_r2_21 * dist_coeffs_4[int(7)]);
                Matrix<float, 2, 2>  _S405 = transpose_0(makeMatrix<float, 2, 2> (_S402 + make_float2 (_S402.x * dist_coeffs_4[int(8)] + _S402.y * dist_coeffs_4[int(9)], 0.0f), _S404 + make_float2 (_S404.x * dist_coeffs_4[int(8)] + _S404.y * dist_coeffs_4[int(9)], 0.0f)));
                bool _S406 = !((F32_min((determinant_0(_S405)), ((F32_min((_S405.rows[int(0)].x), (_S405.rows[int(1)].y)))))) > 0.0f);
                _S382 = _S406;
                if(_S406)
                {
                    break;
                }
                float u_21 = (*mean2d_4).x;
                float v_21 = (*mean2d_4).y;
                float r2_21 = u_21 * u_21 + v_21 * v_21;
                float2  _S407 = *mean2d_4 * make_float2 (1.0f + r2_21 * (dist_coeffs_4[int(0)] + r2_21 * (dist_coeffs_4[int(1)] + r2_21 * (dist_coeffs_4[int(2)] + r2_21 * dist_coeffs_4[int(3)])))) + make_float2 (_S396 * u_21 * v_21 + dist_coeffs_4[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + dist_coeffs_4[int(6)] * r2_21, _S399 * u_21 * v_21 + dist_coeffs_4[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + dist_coeffs_4[int(7)] * r2_21);
                float2  _S408 = _S407 + make_float2 (dist_coeffs_4[int(8)] * _S407.x + dist_coeffs_4[int(9)] * _S407.y, 0.0f);
                *mean2d_4 = make_float2 (fx_4 * _S408.x + cx_2, fy_4 * _S408.y + cy_2);
                break;
            }
            if(!!_S382)
            {
                is_valid_0 = false;
                break;
            }
            Matrix<float, 2, 3>  J_6;
            float2  _S409 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S410;
            (&_S410)->primal_0 = _S369;
            (&_S410)->differential_0 = _S409;
            DiffPair_float_0 _S411 = s_fwd_length_impl_0(&_S410);
            DiffPair_float_0 _S412;
            (&_S412)->primal_0 = _S411.primal_0;
            (&_S412)->differential_0 = _S411.differential_0;
            DiffPair_float_0 _S413;
            (&_S413)->primal_0 = _S384;
            (&_S413)->differential_0 = 0.0f;
            DiffPair_float_0 _S414 = _d_atan2_1(&_S412, &_S413);
            if((_S414.primal_0) < 0.00100000004749745f)
            {
                float _S415 = _S414.differential_0 * _S414.primal_0;
                float _S416 = 1.0f - _S414.primal_0 * _S414.primal_0 / 3.0f;
                float _S417 = ((0.0f - (_S415 + _S415) * 0.3333333432674408f) * _S384 - _S416 * 0.0f) / (_S384 * _S384);
                k_2 = _S416 / _S384;
                eps2d_1 = _S417;
            }
            else
            {
                float _S418 = (_S414.differential_0 * _S411.primal_0 - _S414.primal_0 * _S411.differential_0) / (_S411.primal_0 * _S411.primal_0);
                k_2 = _S414.primal_0 / _S411.primal_0;
                eps2d_1 = _S418;
            }
            float2  _S419 = _S369 * make_float2 (k_2);
            float2  _S420 = _S409 * make_float2 (k_2) + make_float2 (eps2d_1) * _S369;
            float u_22 = _S419.x;
            float s_diff_u_12 = _S420.x;
            float v_22 = _S419.y;
            float s_diff_v_12 = _S420.y;
            float _S421 = s_diff_u_12 * u_22;
            float _S422 = s_diff_v_12 * v_22;
            float r2_22 = u_22 * u_22 + v_22 * v_22;
            float s_diff_r2_22 = _S421 + _S421 + (_S422 + _S422);
            float _S423 = _S372 + r2_22 * _S373;
            float _S424 = _S371 + r2_22 * _S423;
            float _S425 = _S370 + r2_22 * _S424;
            float2  _S426 = _S420 * make_float2 (1.0f + r2_22 * _S425) + make_float2 (s_diff_r2_22 * _S425 + (s_diff_r2_22 * _S424 + (s_diff_r2_22 * _S423 + s_diff_r2_22 * _S373 * r2_22) * r2_22) * r2_22) * _S419 + make_float2 (s_diff_u_12 * _S380 * v_22 + s_diff_v_12 * (_S380 * u_22) + (s_diff_r2_22 + (s_diff_u_12 * 2.0f * u_22 + s_diff_u_12 * (2.0f * u_22))) * _S375 + s_diff_r2_22 * _S376, s_diff_u_12 * _S381 * v_22 + s_diff_v_12 * (_S381 * u_22) + (s_diff_r2_22 + (s_diff_v_12 * 2.0f * v_22 + s_diff_v_12 * (2.0f * v_22))) * _S374 + s_diff_r2_22 * _S377);
            float2  _S427 = _S426 + make_float2 (_S426.x * _S378 + _S426.y * _S379, 0.0f);
            float _S428 = _S427.y * fy_4;
            *&(((&J_6)->rows + (int(0)))->x) = _S427.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->x) = _S428;
            float2  _S429 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S430;
            (&_S430)->primal_0 = _S369;
            (&_S430)->differential_0 = _S429;
            DiffPair_float_0 _S431 = s_fwd_length_impl_0(&_S430);
            DiffPair_float_0 _S432;
            (&_S432)->primal_0 = _S431.primal_0;
            (&_S432)->differential_0 = _S431.differential_0;
            DiffPair_float_0 _S433;
            (&_S433)->primal_0 = _S384;
            (&_S433)->differential_0 = 0.0f;
            DiffPair_float_0 _S434 = _d_atan2_1(&_S432, &_S433);
            if((_S434.primal_0) < 0.00100000004749745f)
            {
                float _S435 = _S434.differential_0 * _S434.primal_0;
                float _S436 = 1.0f - _S434.primal_0 * _S434.primal_0 / 3.0f;
                float _S437 = ((0.0f - (_S435 + _S435) * 0.3333333432674408f) * _S384 - _S436 * 0.0f) / (_S384 * _S384);
                k_2 = _S436 / _S384;
                eps2d_1 = _S437;
            }
            else
            {
                float _S438 = (_S434.differential_0 * _S431.primal_0 - _S434.primal_0 * _S431.differential_0) / (_S431.primal_0 * _S431.primal_0);
                k_2 = _S434.primal_0 / _S431.primal_0;
                eps2d_1 = _S438;
            }
            float2  _S439 = _S369 * make_float2 (k_2);
            float2  _S440 = _S429 * make_float2 (k_2) + make_float2 (eps2d_1) * _S369;
            float u_23 = _S439.x;
            float s_diff_u_13 = _S440.x;
            float v_23 = _S439.y;
            float s_diff_v_13 = _S440.y;
            float _S441 = s_diff_u_13 * u_23;
            float _S442 = s_diff_v_13 * v_23;
            float r2_23 = u_23 * u_23 + v_23 * v_23;
            float s_diff_r2_23 = _S441 + _S441 + (_S442 + _S442);
            float _S443 = _S372 + r2_23 * _S373;
            float _S444 = _S371 + r2_23 * _S443;
            float _S445 = _S370 + r2_23 * _S444;
            float2  _S446 = _S440 * make_float2 (1.0f + r2_23 * _S445) + make_float2 (s_diff_r2_23 * _S445 + (s_diff_r2_23 * _S444 + (s_diff_r2_23 * _S443 + s_diff_r2_23 * _S373 * r2_23) * r2_23) * r2_23) * _S439 + make_float2 (s_diff_u_13 * _S380 * v_23 + s_diff_v_13 * (_S380 * u_23) + (s_diff_r2_23 + (s_diff_u_13 * 2.0f * u_23 + s_diff_u_13 * (2.0f * u_23))) * _S375 + s_diff_r2_23 * _S376, s_diff_u_13 * _S381 * v_23 + s_diff_v_13 * (_S381 * u_23) + (s_diff_r2_23 + (s_diff_v_13 * 2.0f * v_23 + s_diff_v_13 * (2.0f * v_23))) * _S374 + s_diff_r2_23 * _S377);
            float2  _S447 = _S446 + make_float2 (_S446.x * _S378 + _S446.y * _S379, 0.0f);
            float _S448 = _S447.y * fy_4;
            *&(((&J_6)->rows + (int(0)))->y) = _S447.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->y) = _S448;
            float2  _S449 = make_float2 (0.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S450;
            (&_S450)->primal_0 = _S369;
            (&_S450)->differential_0 = _S449;
            DiffPair_float_0 _S451 = s_fwd_length_impl_0(&_S450);
            DiffPair_float_0 _S452;
            (&_S452)->primal_0 = _S451.primal_0;
            (&_S452)->differential_0 = _S451.differential_0;
            DiffPair_float_0 _S453;
            (&_S453)->primal_0 = _S384;
            (&_S453)->differential_0 = 1.0f;
            DiffPair_float_0 _S454 = _d_atan2_1(&_S452, &_S453);
            if((_S454.primal_0) < 0.00100000004749745f)
            {
                float _S455 = _S454.differential_0 * _S454.primal_0;
                float _S456 = 1.0f - _S454.primal_0 * _S454.primal_0 / 3.0f;
                float _S457 = ((0.0f - (_S455 + _S455) * 0.3333333432674408f) * _S384 - _S456) / (_S384 * _S384);
                k_2 = _S456 / _S384;
                eps2d_1 = _S457;
            }
            else
            {
                float _S458 = (_S454.differential_0 * _S451.primal_0 - _S454.primal_0 * _S451.differential_0) / (_S451.primal_0 * _S451.primal_0);
                k_2 = _S454.primal_0 / _S451.primal_0;
                eps2d_1 = _S458;
            }
            float2  _S459 = _S369 * make_float2 (k_2);
            float2  _S460 = _S449 * make_float2 (k_2) + make_float2 (eps2d_1) * _S369;
            float u_24 = _S459.x;
            float s_diff_u_14 = _S460.x;
            float v_24 = _S459.y;
            float s_diff_v_14 = _S460.y;
            float _S461 = s_diff_u_14 * u_24;
            float _S462 = s_diff_v_14 * v_24;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float s_diff_r2_24 = _S461 + _S461 + (_S462 + _S462);
            float _S463 = _S372 + r2_24 * _S373;
            float _S464 = _S371 + r2_24 * _S463;
            float _S465 = _S370 + r2_24 * _S464;
            float2  _S466 = _S460 * make_float2 (1.0f + r2_24 * _S465) + make_float2 (s_diff_r2_24 * _S465 + (s_diff_r2_24 * _S464 + (s_diff_r2_24 * _S463 + s_diff_r2_24 * _S373 * r2_24) * r2_24) * r2_24) * _S459 + make_float2 (s_diff_u_14 * _S380 * v_24 + s_diff_v_14 * (_S380 * u_24) + (s_diff_r2_24 + (s_diff_u_14 * 2.0f * u_24 + s_diff_u_14 * (2.0f * u_24))) * _S375 + s_diff_r2_24 * _S376, s_diff_u_14 * _S381 * v_24 + s_diff_v_14 * (_S381 * u_24) + (s_diff_r2_24 + (s_diff_v_14 * 2.0f * v_24 + s_diff_v_14 * (2.0f * v_24))) * _S374 + s_diff_r2_24 * _S377);
            float2  _S467 = _S466 + make_float2 (_S466.x * _S378 + _S466.y * _S379, 0.0f);
            float _S468 = _S467.y * fy_4;
            *&(((&J_6)->rows + (int(0)))->z) = _S467.x * fx_4;
            *&(((&J_6)->rows + (int(1)))->z) = _S468;
            covar2d_1 = mul_4(mul_3(J_6, covar_c_1), transpose_1(J_6));
            is_valid_0 = true;
            break;
        }
        bool is_valid_1 = true & is_valid_0;
        float2  mean2d_c_0 = *mean2d_4 - make_float2 (cx_2, cy_2);
        float invdet_1 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        float opac_0 = *opacity_1 * (F32_exp((-0.5f * dot_0(mul_7(makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_1, - covar2d_1.rows[int(0)].y * invdet_1, - covar2d_1.rows[int(1)].x * invdet_1, covar2d_1.rows[int(0)].x * invdet_1), mean2d_c_0), mean2d_c_0))));
        if(_S384 < 0.0f)
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
        float _S469 = *&(((&covar2d_1)->rows + (int(0)))->x) + eps2d_1;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S469;
        float _S470 = *&(((&covar2d_1)->rows + (int(1)))->y) + eps2d_1;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S470;
        float det_blur_1 = _S469 * _S470 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S471 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
        if(antialiased_1)
        {
            *opacity_1 = *opacity_1 * compensation_1;
        }
        if((*opacity_1) < 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        float _S472 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_1 / 0.00392156885936856f)))))))));
        float radius_x_1 = _S472 * (F32_sqrt((covar2d_1[int(0)].x)));
        float radius_y_1 = _S472 * (F32_sqrt((covar2d_1[int(1)].y)));
        float _S473 = (*mean2d_4).x - radius_x_1;
        float _S474 = (*mean2d_4).x + radius_x_1;
        float _S475 = (*mean2d_4).y - radius_y_1;
        float _S476 = (*mean2d_4).y + radius_y_1;
        if(_S474 <= 0.0f)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S473 >= float(image_width_1);
        }
        if(is_valid_0)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S476 <= 0.0f;
        }
        if(is_valid_0)
        {
            is_valid_0 = true;
        }
        else
        {
            is_valid_0 = _S475 >= float(image_height_1);
        }
        if(is_valid_0)
        {
            *aabb_xyxy_1 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_1 = make_float4 (_S473, _S475, _S474, _S476);
        float x_17 = mean_c_1.x;
        float y_5 = mean_c_1.y;
        float _S477 = x_17 * x_17 + y_5 * y_5;
        *sorting_depth_1 = _S384 * _S384 * _S384 * _S384 + 0.001953125f * _S477 * _S477;
        *conic_1 = make_float3 (_S471.rows[int(0)].x, _S471.rows[int(0)].y, _S471.rows[int(1)].y);
        *radius_2 = view_radius_3dgs_0(mean_2, scale_1, in_opacity_1, - mul_6(_S387, t_1));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_equisolid(bool antialiased_2, float3  mean_3, float4  quat_2, float3  scale_2, float in_opacity_2, Matrix<float, 3, 3>  R_2, float3  t_2, float fx_5, float fy_5, float cx_3, float cy_3, FixedArray<float, 10>  dist_coeffs_5, uint image_width_2, uint image_height_2, float4  * aabb_xyxy_2, float * sorting_depth_2, float * radius_3, float2  * mean2d_5, float * depth_2, float3  * conic_2, float * opacity_2)
{
    float2  _S478;
    float _S479;
    float _S480;
    float _S481;
    float _S482;
    float _S483;
    float _S484;
    float _S485;
    float _S486;
    float _S487;
    float _S488;
    float _S489;
    float _S490;
    bool _S491;
    for(;;)
    {
        float3  mean_c_2 = mul_6(R_2, mean_3) + t_2;
        float _S492 = length_1(mean_c_2);
        float _S493 = mean_c_2.z;
        *depth_2 = _S493;
        if(_S492 <= 0.0f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        *opacity_2 = 1.0f / (1.0f + (F32_exp((- in_opacity_2))));
        bool is_valid_2;
        float eps2d_2;
        float4  _S494 = normalize_0(quat_2);
        float3  _S495 = exp_0(scale_2);
        float x_18 = _S494.y;
        float x2_2 = x_18 * x_18;
        float y2_2 = _S494.z * _S494.z;
        float z2_2 = _S494.w * _S494.w;
        float xy_2 = _S494.y * _S494.z;
        float xz_2 = _S494.y * _S494.w;
        float yz_2 = _S494.z * _S494.w;
        float wx_2 = _S494.x * _S494.y;
        float wy_2 = _S494.x * _S494.z;
        float wz_2 = _S494.x * _S494.w;
        Matrix<float, 3, 3>  M_2 = mul_5(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (_S495.x, 0.0f, 0.0f, 0.0f, _S495.y, 0.0f, 0.0f, 0.0f, _S495.z));
        Matrix<float, 3, 3>  _S496 = transpose_3(R_2);
        Matrix<float, 3, 3>  covar_c_2 = mul_5(mul_5(R_2, mul_5(M_2, transpose_3(M_2))), _S496);
        for(;;)
        {
            float k_3;
            for(;;)
            {
                float2  _S497 = float2 {mean_c_2.x, mean_c_2.y};
                _S478 = _S497;
                float r_10 = length_0(_S497);
                float theta_3 = (F32_atan2((r_10), (_S493)));
                if(r_10 < 9.99999997475242708e-07f)
                {
                    k_3 = (1.0f - theta_3 * theta_3 / 24.0f) / _S493;
                }
                else
                {
                    k_3 = 2.0f * (F32_sin((0.5f * theta_3))) / r_10;
                }
                float2  _S498 = _S497 * make_float2 (k_3);
                *mean2d_5 = _S498;
                float2  _S499 = make_float2 (1.0f, 0.0f);
                _S479 = dist_coeffs_5[int(0)];
                _S480 = dist_coeffs_5[int(1)];
                _S481 = dist_coeffs_5[int(2)];
                _S482 = dist_coeffs_5[int(3)];
                _S483 = dist_coeffs_5[int(4)];
                _S484 = dist_coeffs_5[int(5)];
                _S485 = dist_coeffs_5[int(6)];
                _S486 = dist_coeffs_5[int(7)];
                _S487 = dist_coeffs_5[int(8)];
                _S488 = dist_coeffs_5[int(9)];
                float u_25 = _S498.x;
                float v_25 = _S498.y;
                float _S500 = 0.0f * v_25;
                float r2_25 = u_25 * u_25 + v_25 * v_25;
                float s_diff_r2_25 = u_25 + u_25 + (_S500 + _S500);
                float _S501 = dist_coeffs_5[int(2)] + r2_25 * dist_coeffs_5[int(3)];
                float _S502 = dist_coeffs_5[int(1)] + r2_25 * _S501;
                float _S503 = dist_coeffs_5[int(0)] + r2_25 * _S502;
                float _S504 = s_diff_r2_25 * _S503 + (s_diff_r2_25 * _S502 + (s_diff_r2_25 * _S501 + s_diff_r2_25 * dist_coeffs_5[int(3)] * r2_25) * r2_25) * r2_25;
                float radial_5 = 1.0f + r2_25 * _S503;
                float _S505 = 2.0f * dist_coeffs_5[int(4)];
                _S489 = _S505;
                float _S506 = _S505 * u_25;
                float _S507 = 2.0f * u_25;
                float s_diff_du_3 = _S505 * v_25 + 0.0f * _S506 + (s_diff_r2_25 + (_S507 + _S507)) * dist_coeffs_5[int(5)] + s_diff_r2_25 * dist_coeffs_5[int(6)];
                float _S508 = 2.0f * dist_coeffs_5[int(5)];
                _S490 = _S508;
                float _S509 = _S508 * u_25;
                float _S510 = 2.0f * v_25;
                float2  _S511 = _S499 * make_float2 (radial_5) + make_float2 (_S504) * _S498 + make_float2 (s_diff_du_3, _S508 * v_25 + 0.0f * _S509 + (s_diff_r2_25 + (_S500 + 0.0f * _S510)) * dist_coeffs_5[int(4)] + s_diff_r2_25 * dist_coeffs_5[int(7)]);
                float _S512 = 0.0f * u_25;
                float s_diff_r2_26 = _S512 + _S512 + (v_25 + v_25);
                float2  _S513 = make_float2 (0.0f, 1.0f) * make_float2 (radial_5) + make_float2 (s_diff_r2_26 * _S503 + (s_diff_r2_26 * _S502 + (s_diff_r2_26 * _S501 + s_diff_r2_26 * dist_coeffs_5[int(3)] * r2_25) * r2_25) * r2_25) * _S498 + make_float2 (0.0f * _S505 * v_25 + _S506 + (s_diff_r2_26 + (_S512 + 0.0f * _S507)) * dist_coeffs_5[int(5)] + s_diff_r2_26 * dist_coeffs_5[int(6)], 0.0f * _S508 * v_25 + _S509 + (s_diff_r2_26 + (_S510 + _S510)) * dist_coeffs_5[int(4)] + s_diff_r2_26 * dist_coeffs_5[int(7)]);
                Matrix<float, 2, 2>  _S514 = transpose_0(makeMatrix<float, 2, 2> (_S511 + make_float2 (_S511.x * dist_coeffs_5[int(8)] + _S511.y * dist_coeffs_5[int(9)], 0.0f), _S513 + make_float2 (_S513.x * dist_coeffs_5[int(8)] + _S513.y * dist_coeffs_5[int(9)], 0.0f)));
                bool _S515 = !((F32_min((determinant_0(_S514)), ((F32_min((_S514.rows[int(0)].x), (_S514.rows[int(1)].y)))))) > 0.0f);
                _S491 = _S515;
                if(_S515)
                {
                    break;
                }
                float u_26 = (*mean2d_5).x;
                float v_26 = (*mean2d_5).y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S516 = *mean2d_5 * make_float2 (1.0f + r2_26 * (dist_coeffs_5[int(0)] + r2_26 * (dist_coeffs_5[int(1)] + r2_26 * (dist_coeffs_5[int(2)] + r2_26 * dist_coeffs_5[int(3)])))) + make_float2 (_S505 * u_26 * v_26 + dist_coeffs_5[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + dist_coeffs_5[int(6)] * r2_26, _S508 * u_26 * v_26 + dist_coeffs_5[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + dist_coeffs_5[int(7)] * r2_26);
                float2  _S517 = _S516 + make_float2 (dist_coeffs_5[int(8)] * _S516.x + dist_coeffs_5[int(9)] * _S516.y, 0.0f);
                *mean2d_5 = make_float2 (fx_5 * _S517.x + cx_3, fy_5 * _S517.y + cy_3);
                break;
            }
            if(!!_S491)
            {
                is_valid_2 = false;
                break;
            }
            Matrix<float, 2, 3>  J_7;
            float2  _S518 = make_float2 (1.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S519;
            (&_S519)->primal_0 = _S478;
            (&_S519)->differential_0 = _S518;
            DiffPair_float_0 _S520 = s_fwd_length_impl_0(&_S519);
            DiffPair_float_0 _S521;
            (&_S521)->primal_0 = _S520.primal_0;
            (&_S521)->differential_0 = _S520.differential_0;
            DiffPair_float_0 _S522;
            (&_S522)->primal_0 = _S493;
            (&_S522)->differential_0 = 0.0f;
            DiffPair_float_0 _S523 = _d_atan2_1(&_S521, &_S522);
            if((_S520.primal_0) < 9.99999997475242708e-07f)
            {
                float _S524 = _S523.differential_0 * _S523.primal_0;
                float _S525 = 1.0f - _S523.primal_0 * _S523.primal_0 / 24.0f;
                float _S526 = ((0.0f - (_S524 + _S524) * 0.0416666679084301f) * _S493 - _S525 * 0.0f) / (_S493 * _S493);
                k_3 = _S525 / _S493;
                eps2d_2 = _S526;
            }
            else
            {
                float _S527 = _S523.differential_0 * 0.5f;
                DiffPair_float_0 _S528;
                (&_S528)->primal_0 = 0.5f * _S523.primal_0;
                (&_S528)->differential_0 = _S527;
                DiffPair_float_0 _S529 = _d_sin_1(&_S528);
                float _S530 = 2.0f * _S529.primal_0;
                float _S531 = (_S529.differential_0 * 2.0f * _S520.primal_0 - _S530 * _S520.differential_0) / (_S520.primal_0 * _S520.primal_0);
                k_3 = _S530 / _S520.primal_0;
                eps2d_2 = _S531;
            }
            float2  _S532 = _S478 * make_float2 (k_3);
            float2  _S533 = _S518 * make_float2 (k_3) + make_float2 (eps2d_2) * _S478;
            float u_27 = _S532.x;
            float s_diff_u_15 = _S533.x;
            float v_27 = _S532.y;
            float s_diff_v_15 = _S533.y;
            float _S534 = s_diff_u_15 * u_27;
            float _S535 = s_diff_v_15 * v_27;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float s_diff_r2_27 = _S534 + _S534 + (_S535 + _S535);
            float _S536 = _S481 + r2_27 * _S482;
            float _S537 = _S480 + r2_27 * _S536;
            float _S538 = _S479 + r2_27 * _S537;
            float2  _S539 = _S533 * make_float2 (1.0f + r2_27 * _S538) + make_float2 (s_diff_r2_27 * _S538 + (s_diff_r2_27 * _S537 + (s_diff_r2_27 * _S536 + s_diff_r2_27 * _S482 * r2_27) * r2_27) * r2_27) * _S532 + make_float2 (s_diff_u_15 * _S489 * v_27 + s_diff_v_15 * (_S489 * u_27) + (s_diff_r2_27 + (s_diff_u_15 * 2.0f * u_27 + s_diff_u_15 * (2.0f * u_27))) * _S484 + s_diff_r2_27 * _S485, s_diff_u_15 * _S490 * v_27 + s_diff_v_15 * (_S490 * u_27) + (s_diff_r2_27 + (s_diff_v_15 * 2.0f * v_27 + s_diff_v_15 * (2.0f * v_27))) * _S483 + s_diff_r2_27 * _S486);
            float2  _S540 = _S539 + make_float2 (_S539.x * _S487 + _S539.y * _S488, 0.0f);
            float _S541 = _S540.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->x) = _S540.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->x) = _S541;
            float2  _S542 = make_float2 (0.0f, 1.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S543;
            (&_S543)->primal_0 = _S478;
            (&_S543)->differential_0 = _S542;
            DiffPair_float_0 _S544 = s_fwd_length_impl_0(&_S543);
            DiffPair_float_0 _S545;
            (&_S545)->primal_0 = _S544.primal_0;
            (&_S545)->differential_0 = _S544.differential_0;
            DiffPair_float_0 _S546;
            (&_S546)->primal_0 = _S493;
            (&_S546)->differential_0 = 0.0f;
            DiffPair_float_0 _S547 = _d_atan2_1(&_S545, &_S546);
            if((_S544.primal_0) < 9.99999997475242708e-07f)
            {
                float _S548 = _S547.differential_0 * _S547.primal_0;
                float _S549 = 1.0f - _S547.primal_0 * _S547.primal_0 / 24.0f;
                float _S550 = ((0.0f - (_S548 + _S548) * 0.0416666679084301f) * _S493 - _S549 * 0.0f) / (_S493 * _S493);
                k_3 = _S549 / _S493;
                eps2d_2 = _S550;
            }
            else
            {
                float _S551 = _S547.differential_0 * 0.5f;
                DiffPair_float_0 _S552;
                (&_S552)->primal_0 = 0.5f * _S547.primal_0;
                (&_S552)->differential_0 = _S551;
                DiffPair_float_0 _S553 = _d_sin_1(&_S552);
                float _S554 = 2.0f * _S553.primal_0;
                float _S555 = (_S553.differential_0 * 2.0f * _S544.primal_0 - _S554 * _S544.differential_0) / (_S544.primal_0 * _S544.primal_0);
                k_3 = _S554 / _S544.primal_0;
                eps2d_2 = _S555;
            }
            float2  _S556 = _S478 * make_float2 (k_3);
            float2  _S557 = _S542 * make_float2 (k_3) + make_float2 (eps2d_2) * _S478;
            float u_28 = _S556.x;
            float s_diff_u_16 = _S557.x;
            float v_28 = _S556.y;
            float s_diff_v_16 = _S557.y;
            float _S558 = s_diff_u_16 * u_28;
            float _S559 = s_diff_v_16 * v_28;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float s_diff_r2_28 = _S558 + _S558 + (_S559 + _S559);
            float _S560 = _S481 + r2_28 * _S482;
            float _S561 = _S480 + r2_28 * _S560;
            float _S562 = _S479 + r2_28 * _S561;
            float2  _S563 = _S557 * make_float2 (1.0f + r2_28 * _S562) + make_float2 (s_diff_r2_28 * _S562 + (s_diff_r2_28 * _S561 + (s_diff_r2_28 * _S560 + s_diff_r2_28 * _S482 * r2_28) * r2_28) * r2_28) * _S556 + make_float2 (s_diff_u_16 * _S489 * v_28 + s_diff_v_16 * (_S489 * u_28) + (s_diff_r2_28 + (s_diff_u_16 * 2.0f * u_28 + s_diff_u_16 * (2.0f * u_28))) * _S484 + s_diff_r2_28 * _S485, s_diff_u_16 * _S490 * v_28 + s_diff_v_16 * (_S490 * u_28) + (s_diff_r2_28 + (s_diff_v_16 * 2.0f * v_28 + s_diff_v_16 * (2.0f * v_28))) * _S483 + s_diff_r2_28 * _S486);
            float2  _S564 = _S563 + make_float2 (_S563.x * _S487 + _S563.y * _S488, 0.0f);
            float _S565 = _S564.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->y) = _S564.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->y) = _S565;
            float2  _S566 = make_float2 (0.0f, 0.0f);
            DiffPair_vectorx3Cfloatx2C2x3E_0 _S567;
            (&_S567)->primal_0 = _S478;
            (&_S567)->differential_0 = _S566;
            DiffPair_float_0 _S568 = s_fwd_length_impl_0(&_S567);
            DiffPair_float_0 _S569;
            (&_S569)->primal_0 = _S568.primal_0;
            (&_S569)->differential_0 = _S568.differential_0;
            DiffPair_float_0 _S570;
            (&_S570)->primal_0 = _S493;
            (&_S570)->differential_0 = 1.0f;
            DiffPair_float_0 _S571 = _d_atan2_1(&_S569, &_S570);
            if((_S568.primal_0) < 9.99999997475242708e-07f)
            {
                float _S572 = _S571.differential_0 * _S571.primal_0;
                float _S573 = 1.0f - _S571.primal_0 * _S571.primal_0 / 24.0f;
                float _S574 = ((0.0f - (_S572 + _S572) * 0.0416666679084301f) * _S493 - _S573) / (_S493 * _S493);
                k_3 = _S573 / _S493;
                eps2d_2 = _S574;
            }
            else
            {
                float _S575 = _S571.differential_0 * 0.5f;
                DiffPair_float_0 _S576;
                (&_S576)->primal_0 = 0.5f * _S571.primal_0;
                (&_S576)->differential_0 = _S575;
                DiffPair_float_0 _S577 = _d_sin_1(&_S576);
                float _S578 = 2.0f * _S577.primal_0;
                float _S579 = (_S577.differential_0 * 2.0f * _S568.primal_0 - _S578 * _S568.differential_0) / (_S568.primal_0 * _S568.primal_0);
                k_3 = _S578 / _S568.primal_0;
                eps2d_2 = _S579;
            }
            float2  _S580 = _S478 * make_float2 (k_3);
            float2  _S581 = _S566 * make_float2 (k_3) + make_float2 (eps2d_2) * _S478;
            float u_29 = _S580.x;
            float s_diff_u_17 = _S581.x;
            float v_29 = _S580.y;
            float s_diff_v_17 = _S581.y;
            float _S582 = s_diff_u_17 * u_29;
            float _S583 = s_diff_v_17 * v_29;
            float r2_29 = u_29 * u_29 + v_29 * v_29;
            float s_diff_r2_29 = _S582 + _S582 + (_S583 + _S583);
            float _S584 = _S481 + r2_29 * _S482;
            float _S585 = _S480 + r2_29 * _S584;
            float _S586 = _S479 + r2_29 * _S585;
            float2  _S587 = _S581 * make_float2 (1.0f + r2_29 * _S586) + make_float2 (s_diff_r2_29 * _S586 + (s_diff_r2_29 * _S585 + (s_diff_r2_29 * _S584 + s_diff_r2_29 * _S482 * r2_29) * r2_29) * r2_29) * _S580 + make_float2 (s_diff_u_17 * _S489 * v_29 + s_diff_v_17 * (_S489 * u_29) + (s_diff_r2_29 + (s_diff_u_17 * 2.0f * u_29 + s_diff_u_17 * (2.0f * u_29))) * _S484 + s_diff_r2_29 * _S485, s_diff_u_17 * _S490 * v_29 + s_diff_v_17 * (_S490 * u_29) + (s_diff_r2_29 + (s_diff_v_17 * 2.0f * v_29 + s_diff_v_17 * (2.0f * v_29))) * _S483 + s_diff_r2_29 * _S486);
            float2  _S588 = _S587 + make_float2 (_S587.x * _S487 + _S587.y * _S488, 0.0f);
            float _S589 = _S588.y * fy_5;
            *&(((&J_7)->rows + (int(0)))->z) = _S588.x * fx_5;
            *&(((&J_7)->rows + (int(1)))->z) = _S589;
            covar2d_2 = mul_4(mul_3(J_7, covar_c_2), transpose_1(J_7));
            is_valid_2 = true;
            break;
        }
        bool is_valid_3 = true & is_valid_2;
        float2  mean2d_c_1 = *mean2d_5 - make_float2 (cx_3, cy_3);
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        float opac_1 = *opacity_2 * (F32_exp((-0.5f * dot_0(mul_7(makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3), mean2d_c_1), mean2d_c_1))));
        if(_S493 < 0.0f)
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
        float _S590 = *&(((&covar2d_2)->rows + (int(0)))->x) + eps2d_2;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S590;
        float _S591 = *&(((&covar2d_2)->rows + (int(1)))->y) + eps2d_2;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S591;
        float det_blur_2 = _S590 * _S591 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        float invdet_4 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S592 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_4, - covar2d_2.rows[int(0)].y * invdet_4, - covar2d_2.rows[int(1)].x * invdet_4, covar2d_2.rows[int(0)].x * invdet_4);
        if(antialiased_2)
        {
            *opacity_2 = *opacity_2 * compensation_2;
        }
        if((*opacity_2) < 0.00392156885936856f)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        float _S593 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_2 / 0.00392156885936856f)))))))));
        float radius_x_2 = _S593 * (F32_sqrt((covar2d_2[int(0)].x)));
        float radius_y_2 = _S593 * (F32_sqrt((covar2d_2[int(1)].y)));
        float _S594 = (*mean2d_5).x - radius_x_2;
        float _S595 = (*mean2d_5).x + radius_x_2;
        float _S596 = (*mean2d_5).y - radius_y_2;
        float _S597 = (*mean2d_5).y + radius_y_2;
        if(_S595 <= 0.0f)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S594 >= float(image_width_2);
        }
        if(is_valid_2)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S597 <= 0.0f;
        }
        if(is_valid_2)
        {
            is_valid_2 = true;
        }
        else
        {
            is_valid_2 = _S596 >= float(image_height_2);
        }
        if(is_valid_2)
        {
            *aabb_xyxy_2 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_2 = make_float4 (_S594, _S596, _S595, _S597);
        float x_19 = mean_c_2.x;
        float y_6 = mean_c_2.y;
        float _S598 = x_19 * x_19 + y_6 * y_6;
        *sorting_depth_2 = _S493 * _S493 * _S493 * _S493 + 0.001953125f * _S598 * _S598;
        *conic_2 = make_float3 (_S592.rows[int(0)].x, _S592.rows[int(0)].y, _S592.rows[int(1)].y);
        *radius_3 = view_radius_3dgs_0(mean_3, scale_2, in_opacity_2, - mul_6(_S496, t_2));
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

inline __device__ void projection_3dgut_persp(bool antialiased_3, float3  mean_4, float4  quat_3, float3  scale_3, float in_opacity_3, Matrix<float, 3, 3>  R_3, float3  t_3, float fx_6, float fy_6, float cx_4, float cy_4, FixedArray<float, 10>  dist_coeffs_6, uint image_width_3, uint image_height_3, float4  * aabb_xyxy_3, float * sorting_depth_3, float * radius_4, float2  * mean2d_6, float * depth_3, float3  * conic_3, float * opacity_3)
{
    float _S599;
    float _S600;
    float2  * _S601;
    float2  * _S602;
    float2  * _S603;
    bool _S604;
    float2  * _S605;
    float2  * _S606;
    float2  * _S607;
    bool _S608;
    float2  * _S609;
    bool _S610;
    for(;;)
    {
        float _S611 = (mul_6(R_3, mean_4) + t_3).z;
        *depth_3 = _S611;
        if(_S611 <= 0.0f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_3;
        *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
        bool _S612;
        float3  _S613 = exp_0(scale_3);
        float4  _S614 = normalize_0(quat_3);
        float x_20 = _S614.y;
        float x2_3 = x_20 * x_20;
        float y2_3 = _S614.z * _S614.z;
        float z2_3 = _S614.w * _S614.w;
        float xy_3 = _S614.y * _S614.z;
        float xz_3 = _S614.y * _S614.w;
        float yz_3 = _S614.z * _S614.w;
        float wx_3 = _S614.x * _S614.y;
        float wy_3 = _S614.x * _S614.z;
        float wz_3 = _S614.x * _S614.w;
        Matrix<float, 3, 3>  _S615 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_4;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S616 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S616 * _S613.x) * _S615.rows[0U];
        float3  _S617 = mean_4 + delta_0;
        float3  _S618 = mean_4 - delta_0;
        float3  delta_1 = make_float3 (_S616 * _S613.y) * _S615.rows[1U];
        float3  _S619 = mean_4 + delta_1;
        float3  _S620 = mean_4 - delta_1;
        float3  delta_2 = make_float3 (_S616 * _S613.z) * _S615.rows[2U];
        float3  _S621 = mean_4 + delta_2;
        float3  _S622 = mean_4 - delta_2;
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
        (&ret_0)->p_0[0U] = mul_6(R_3, (&ret_0)->p_0[0U]) + t_3;
        (&ret_0)->p_0[1U] = mul_6(R_3, _S617) + t_3;
        (&ret_0)->p_0[2U] = mul_6(R_3, _S619) + t_3;
        (&ret_0)->p_0[3U] = mul_6(R_3, _S621) + t_3;
        (&ret_0)->p_0[4U] = mul_6(R_3, _S618) + t_3;
        (&ret_0)->p_0[5U] = mul_6(R_3, _S620) + t_3;
        (&ret_0)->p_0[6U] = mul_6(R_3, _S622) + t_3;
        SigmaPoints_0 _S623 = ret_0;
        for(;;)
        {
            int2  _S624 = make_int2 (int(0));
            float2  _S625 = make_float2 ((float)_S624.x, (float)_S624.y);
            *mean2d_6 = _S625;
            covar2d_3 = makeMatrix<float, 2, 2> (0.0f);
            float _S626 = float(image_width_3);
            _S599 = _S626;
            float tan_fovx_0 = 0.5f * _S626 / fx_6;
            float _S627 = float(image_height_3);
            _S600 = _S627;
            float _S628 = 0.30000001192092896f * tan_fovx_0 * fx_6;
            float lim_x_pos_2 = _S626 + _S628;
            float _S629 = 0.30000001192092896f * (0.5f * _S627 / fy_6) * fy_6;
            float lim_y_pos_0 = _S627 + _S629;
            FixedArray<float2 , 7>  proj_points_0;
            for(;;)
            {
                _S601 = &proj_points_0[int(0)];
                for(;;)
                {
                    float _S630 = _S623.p_0[int(0)].z;
                    proj_points_0[int(0)] = float2 {_S623.p_0[int(0)].x, _S623.p_0[int(0)].y} / make_float2 (_S630);
                    if(_S630 < 0.0f)
                    {
                        _S612 = true;
                    }
                    else
                    {
                        float u_30 = proj_points_0[int(0)].x;
                        float v_30 = proj_points_0[int(0)].y;
                        float _S631 = 0.0f * v_30;
                        float r2_30 = u_30 * u_30 + v_30 * v_30;
                        float s_diff_r2_30 = u_30 + u_30 + (_S631 + _S631);
                        float _S632 = dist_coeffs_6[int(2)] + r2_30 * dist_coeffs_6[int(3)];
                        float _S633 = dist_coeffs_6[int(1)] + r2_30 * _S632;
                        float _S634 = dist_coeffs_6[int(0)] + r2_30 * _S633;
                        float radial_6 = 1.0f + r2_30 * _S634;
                        float _S635 = 2.0f * dist_coeffs_6[int(4)];
                        float _S636 = _S635 * u_30;
                        float _S637 = 2.0f * u_30;
                        float _S638 = 2.0f * dist_coeffs_6[int(5)];
                        float _S639 = _S638 * u_30;
                        float _S640 = 2.0f * v_30;
                        float2  _S641 = make_float2 (1.0f, 0.0f) * make_float2 (radial_6) + make_float2 (s_diff_r2_30 * _S634 + (s_diff_r2_30 * _S633 + (s_diff_r2_30 * _S632 + s_diff_r2_30 * dist_coeffs_6[int(3)] * r2_30) * r2_30) * r2_30) * proj_points_0[int(0)] + make_float2 (_S635 * v_30 + 0.0f * _S636 + (s_diff_r2_30 + (_S637 + _S637)) * dist_coeffs_6[int(5)] + s_diff_r2_30 * dist_coeffs_6[int(6)], _S638 * v_30 + 0.0f * _S639 + (s_diff_r2_30 + (_S631 + 0.0f * _S640)) * dist_coeffs_6[int(4)] + s_diff_r2_30 * dist_coeffs_6[int(7)]);
                        float _S642 = 0.0f * u_30;
                        float s_diff_r2_31 = _S642 + _S642 + (v_30 + v_30);
                        float2  _S643 = make_float2 (0.0f, 1.0f) * make_float2 (radial_6) + make_float2 (s_diff_r2_31 * _S634 + (s_diff_r2_31 * _S633 + (s_diff_r2_31 * _S632 + s_diff_r2_31 * dist_coeffs_6[int(3)] * r2_30) * r2_30) * r2_30) * proj_points_0[int(0)] + make_float2 (0.0f * _S635 * v_30 + _S636 + (s_diff_r2_31 + (_S642 + 0.0f * _S637)) * dist_coeffs_6[int(5)] + s_diff_r2_31 * dist_coeffs_6[int(6)], 0.0f * _S638 * v_30 + _S639 + (s_diff_r2_31 + (_S640 + _S640)) * dist_coeffs_6[int(4)] + s_diff_r2_31 * dist_coeffs_6[int(7)]);
                        Matrix<float, 2, 2>  _S644 = transpose_0(makeMatrix<float, 2, 2> (_S641 + make_float2 (_S641.x * dist_coeffs_6[int(8)] + _S641.y * dist_coeffs_6[int(9)], 0.0f), _S643 + make_float2 (_S643.x * dist_coeffs_6[int(8)] + _S643.y * dist_coeffs_6[int(9)], 0.0f)));
                        _S612 = !((F32_min((determinant_0(_S644)), ((F32_min((_S644.rows[int(0)].x), (_S644.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S612)
                    {
                        break;
                    }
                    float u_31 = proj_points_0[int(0)].x;
                    float v_31 = proj_points_0[int(0)].y;
                    float r2_31 = u_31 * u_31 + v_31 * v_31;
                    float2  _S645 = proj_points_0[int(0)] * make_float2 (1.0f + r2_31 * (dist_coeffs_6[int(0)] + r2_31 * (dist_coeffs_6[int(1)] + r2_31 * (dist_coeffs_6[int(2)] + r2_31 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_31 * v_31 + dist_coeffs_6[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + dist_coeffs_6[int(6)] * r2_31, 2.0f * dist_coeffs_6[int(5)] * u_31 * v_31 + dist_coeffs_6[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + dist_coeffs_6[int(7)] * r2_31);
                    float2  _S646 = _S645 + make_float2 (dist_coeffs_6[int(8)] * _S645.x + dist_coeffs_6[int(9)] * _S645.y, 0.0f);
                    proj_points_0[int(0)] = make_float2 (fx_6 * _S646.x + cx_4, fy_6 * _S646.y + cy_4);
                    break;
                }
                bool all_valid_0 = true & (!_S612);
                _S602 = &proj_points_0[int(1)];
                for(;;)
                {
                    float _S647 = _S623.p_0[int(1)].z;
                    proj_points_0[int(1)] = float2 {_S623.p_0[int(1)].x, _S623.p_0[int(1)].y} / make_float2 (_S647);
                    if(_S647 < 0.0f)
                    {
                        _S612 = true;
                    }
                    else
                    {
                        float u_32 = proj_points_0[int(1)].x;
                        float v_32 = proj_points_0[int(1)].y;
                        float _S648 = 0.0f * v_32;
                        float r2_32 = u_32 * u_32 + v_32 * v_32;
                        float s_diff_r2_32 = u_32 + u_32 + (_S648 + _S648);
                        float _S649 = dist_coeffs_6[int(2)] + r2_32 * dist_coeffs_6[int(3)];
                        float _S650 = dist_coeffs_6[int(1)] + r2_32 * _S649;
                        float _S651 = dist_coeffs_6[int(0)] + r2_32 * _S650;
                        float radial_7 = 1.0f + r2_32 * _S651;
                        float _S652 = 2.0f * dist_coeffs_6[int(4)];
                        float _S653 = _S652 * u_32;
                        float _S654 = 2.0f * u_32;
                        float _S655 = 2.0f * dist_coeffs_6[int(5)];
                        float _S656 = _S655 * u_32;
                        float _S657 = 2.0f * v_32;
                        float2  _S658 = make_float2 (1.0f, 0.0f) * make_float2 (radial_7) + make_float2 (s_diff_r2_32 * _S651 + (s_diff_r2_32 * _S650 + (s_diff_r2_32 * _S649 + s_diff_r2_32 * dist_coeffs_6[int(3)] * r2_32) * r2_32) * r2_32) * proj_points_0[int(1)] + make_float2 (_S652 * v_32 + 0.0f * _S653 + (s_diff_r2_32 + (_S654 + _S654)) * dist_coeffs_6[int(5)] + s_diff_r2_32 * dist_coeffs_6[int(6)], _S655 * v_32 + 0.0f * _S656 + (s_diff_r2_32 + (_S648 + 0.0f * _S657)) * dist_coeffs_6[int(4)] + s_diff_r2_32 * dist_coeffs_6[int(7)]);
                        float _S659 = 0.0f * u_32;
                        float s_diff_r2_33 = _S659 + _S659 + (v_32 + v_32);
                        float2  _S660 = make_float2 (0.0f, 1.0f) * make_float2 (radial_7) + make_float2 (s_diff_r2_33 * _S651 + (s_diff_r2_33 * _S650 + (s_diff_r2_33 * _S649 + s_diff_r2_33 * dist_coeffs_6[int(3)] * r2_32) * r2_32) * r2_32) * proj_points_0[int(1)] + make_float2 (0.0f * _S652 * v_32 + _S653 + (s_diff_r2_33 + (_S659 + 0.0f * _S654)) * dist_coeffs_6[int(5)] + s_diff_r2_33 * dist_coeffs_6[int(6)], 0.0f * _S655 * v_32 + _S656 + (s_diff_r2_33 + (_S657 + _S657)) * dist_coeffs_6[int(4)] + s_diff_r2_33 * dist_coeffs_6[int(7)]);
                        Matrix<float, 2, 2>  _S661 = transpose_0(makeMatrix<float, 2, 2> (_S658 + make_float2 (_S658.x * dist_coeffs_6[int(8)] + _S658.y * dist_coeffs_6[int(9)], 0.0f), _S660 + make_float2 (_S660.x * dist_coeffs_6[int(8)] + _S660.y * dist_coeffs_6[int(9)], 0.0f)));
                        _S612 = !((F32_min((determinant_0(_S661)), ((F32_min((_S661.rows[int(0)].x), (_S661.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S612)
                    {
                        break;
                    }
                    float u_33 = proj_points_0[int(1)].x;
                    float v_33 = proj_points_0[int(1)].y;
                    float r2_33 = u_33 * u_33 + v_33 * v_33;
                    float2  _S662 = proj_points_0[int(1)] * make_float2 (1.0f + r2_33 * (dist_coeffs_6[int(0)] + r2_33 * (dist_coeffs_6[int(1)] + r2_33 * (dist_coeffs_6[int(2)] + r2_33 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_33 * v_33 + dist_coeffs_6[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + dist_coeffs_6[int(6)] * r2_33, 2.0f * dist_coeffs_6[int(5)] * u_33 * v_33 + dist_coeffs_6[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + dist_coeffs_6[int(7)] * r2_33);
                    float2  _S663 = _S662 + make_float2 (dist_coeffs_6[int(8)] * _S662.x + dist_coeffs_6[int(9)] * _S662.y, 0.0f);
                    proj_points_0[int(1)] = make_float2 (fx_6 * _S663.x + cx_4, fy_6 * _S663.y + cy_4);
                    break;
                }
                bool all_valid_1 = all_valid_0 & (!_S612);
                for(;;)
                {
                    _S603 = &proj_points_0[int(2)];
                    for(;;)
                    {
                        float _S664 = _S623.p_0[int(2)].z;
                        proj_points_0[int(2)] = float2 {_S623.p_0[int(2)].x, _S623.p_0[int(2)].y} / make_float2 (_S664);
                        if(_S664 < 0.0f)
                        {
                            _S612 = true;
                        }
                        else
                        {
                            float u_34 = proj_points_0[int(2)].x;
                            float v_34 = proj_points_0[int(2)].y;
                            float _S665 = 0.0f * v_34;
                            float r2_34 = u_34 * u_34 + v_34 * v_34;
                            float s_diff_r2_34 = u_34 + u_34 + (_S665 + _S665);
                            float _S666 = dist_coeffs_6[int(2)] + r2_34 * dist_coeffs_6[int(3)];
                            float _S667 = dist_coeffs_6[int(1)] + r2_34 * _S666;
                            float _S668 = dist_coeffs_6[int(0)] + r2_34 * _S667;
                            float radial_8 = 1.0f + r2_34 * _S668;
                            float _S669 = 2.0f * dist_coeffs_6[int(4)];
                            float _S670 = _S669 * u_34;
                            float _S671 = 2.0f * u_34;
                            float _S672 = 2.0f * dist_coeffs_6[int(5)];
                            float _S673 = _S672 * u_34;
                            float _S674 = 2.0f * v_34;
                            float2  _S675 = make_float2 (1.0f, 0.0f) * make_float2 (radial_8) + make_float2 (s_diff_r2_34 * _S668 + (s_diff_r2_34 * _S667 + (s_diff_r2_34 * _S666 + s_diff_r2_34 * dist_coeffs_6[int(3)] * r2_34) * r2_34) * r2_34) * proj_points_0[int(2)] + make_float2 (_S669 * v_34 + 0.0f * _S670 + (s_diff_r2_34 + (_S671 + _S671)) * dist_coeffs_6[int(5)] + s_diff_r2_34 * dist_coeffs_6[int(6)], _S672 * v_34 + 0.0f * _S673 + (s_diff_r2_34 + (_S665 + 0.0f * _S674)) * dist_coeffs_6[int(4)] + s_diff_r2_34 * dist_coeffs_6[int(7)]);
                            float _S676 = 0.0f * u_34;
                            float s_diff_r2_35 = _S676 + _S676 + (v_34 + v_34);
                            float2  _S677 = make_float2 (0.0f, 1.0f) * make_float2 (radial_8) + make_float2 (s_diff_r2_35 * _S668 + (s_diff_r2_35 * _S667 + (s_diff_r2_35 * _S666 + s_diff_r2_35 * dist_coeffs_6[int(3)] * r2_34) * r2_34) * r2_34) * proj_points_0[int(2)] + make_float2 (0.0f * _S669 * v_34 + _S670 + (s_diff_r2_35 + (_S676 + 0.0f * _S671)) * dist_coeffs_6[int(5)] + s_diff_r2_35 * dist_coeffs_6[int(6)], 0.0f * _S672 * v_34 + _S673 + (s_diff_r2_35 + (_S674 + _S674)) * dist_coeffs_6[int(4)] + s_diff_r2_35 * dist_coeffs_6[int(7)]);
                            Matrix<float, 2, 2>  _S678 = transpose_0(makeMatrix<float, 2, 2> (_S675 + make_float2 (_S675.x * dist_coeffs_6[int(8)] + _S675.y * dist_coeffs_6[int(9)], 0.0f), _S677 + make_float2 (_S677.x * dist_coeffs_6[int(8)] + _S677.y * dist_coeffs_6[int(9)], 0.0f)));
                            _S612 = !((F32_min((determinant_0(_S678)), ((F32_min((_S678.rows[int(0)].x), (_S678.rows[int(1)].y)))))) > 0.0f);
                        }
                        if(_S612)
                        {
                            break;
                        }
                        float u_35 = proj_points_0[int(2)].x;
                        float v_35 = proj_points_0[int(2)].y;
                        float r2_35 = u_35 * u_35 + v_35 * v_35;
                        float2  _S679 = proj_points_0[int(2)] * make_float2 (1.0f + r2_35 * (dist_coeffs_6[int(0)] + r2_35 * (dist_coeffs_6[int(1)] + r2_35 * (dist_coeffs_6[int(2)] + r2_35 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_35 * v_35 + dist_coeffs_6[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + dist_coeffs_6[int(6)] * r2_35, 2.0f * dist_coeffs_6[int(5)] * u_35 * v_35 + dist_coeffs_6[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + dist_coeffs_6[int(7)] * r2_35);
                        float2  _S680 = _S679 + make_float2 (dist_coeffs_6[int(8)] * _S679.x + dist_coeffs_6[int(9)] * _S679.y, 0.0f);
                        proj_points_0[int(2)] = make_float2 (fx_6 * _S680.x + cx_4, fy_6 * _S680.y + cy_4);
                        break;
                    }
                    _S604 = all_valid_1 & (!_S612);
                    break;
                }
                _S605 = &proj_points_0[int(3)];
                for(;;)
                {
                    float _S681 = _S623.p_0[int(3)].z;
                    proj_points_0[int(3)] = float2 {_S623.p_0[int(3)].x, _S623.p_0[int(3)].y} / make_float2 (_S681);
                    if(_S681 < 0.0f)
                    {
                        _S612 = true;
                    }
                    else
                    {
                        float u_36 = proj_points_0[int(3)].x;
                        float v_36 = proj_points_0[int(3)].y;
                        float _S682 = 0.0f * v_36;
                        float r2_36 = u_36 * u_36 + v_36 * v_36;
                        float s_diff_r2_36 = u_36 + u_36 + (_S682 + _S682);
                        float _S683 = dist_coeffs_6[int(2)] + r2_36 * dist_coeffs_6[int(3)];
                        float _S684 = dist_coeffs_6[int(1)] + r2_36 * _S683;
                        float _S685 = dist_coeffs_6[int(0)] + r2_36 * _S684;
                        float radial_9 = 1.0f + r2_36 * _S685;
                        float _S686 = 2.0f * dist_coeffs_6[int(4)];
                        float _S687 = _S686 * u_36;
                        float _S688 = 2.0f * u_36;
                        float _S689 = 2.0f * dist_coeffs_6[int(5)];
                        float _S690 = _S689 * u_36;
                        float _S691 = 2.0f * v_36;
                        float2  _S692 = make_float2 (1.0f, 0.0f) * make_float2 (radial_9) + make_float2 (s_diff_r2_36 * _S685 + (s_diff_r2_36 * _S684 + (s_diff_r2_36 * _S683 + s_diff_r2_36 * dist_coeffs_6[int(3)] * r2_36) * r2_36) * r2_36) * proj_points_0[int(3)] + make_float2 (_S686 * v_36 + 0.0f * _S687 + (s_diff_r2_36 + (_S688 + _S688)) * dist_coeffs_6[int(5)] + s_diff_r2_36 * dist_coeffs_6[int(6)], _S689 * v_36 + 0.0f * _S690 + (s_diff_r2_36 + (_S682 + 0.0f * _S691)) * dist_coeffs_6[int(4)] + s_diff_r2_36 * dist_coeffs_6[int(7)]);
                        float _S693 = 0.0f * u_36;
                        float s_diff_r2_37 = _S693 + _S693 + (v_36 + v_36);
                        float2  _S694 = make_float2 (0.0f, 1.0f) * make_float2 (radial_9) + make_float2 (s_diff_r2_37 * _S685 + (s_diff_r2_37 * _S684 + (s_diff_r2_37 * _S683 + s_diff_r2_37 * dist_coeffs_6[int(3)] * r2_36) * r2_36) * r2_36) * proj_points_0[int(3)] + make_float2 (0.0f * _S686 * v_36 + _S687 + (s_diff_r2_37 + (_S693 + 0.0f * _S688)) * dist_coeffs_6[int(5)] + s_diff_r2_37 * dist_coeffs_6[int(6)], 0.0f * _S689 * v_36 + _S690 + (s_diff_r2_37 + (_S691 + _S691)) * dist_coeffs_6[int(4)] + s_diff_r2_37 * dist_coeffs_6[int(7)]);
                        Matrix<float, 2, 2>  _S695 = transpose_0(makeMatrix<float, 2, 2> (_S692 + make_float2 (_S692.x * dist_coeffs_6[int(8)] + _S692.y * dist_coeffs_6[int(9)], 0.0f), _S694 + make_float2 (_S694.x * dist_coeffs_6[int(8)] + _S694.y * dist_coeffs_6[int(9)], 0.0f)));
                        _S612 = !((F32_min((determinant_0(_S695)), ((F32_min((_S695.rows[int(0)].x), (_S695.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S612)
                    {
                        break;
                    }
                    float u_37 = proj_points_0[int(3)].x;
                    float v_37 = proj_points_0[int(3)].y;
                    float r2_37 = u_37 * u_37 + v_37 * v_37;
                    float2  _S696 = proj_points_0[int(3)] * make_float2 (1.0f + r2_37 * (dist_coeffs_6[int(0)] + r2_37 * (dist_coeffs_6[int(1)] + r2_37 * (dist_coeffs_6[int(2)] + r2_37 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_37 * v_37 + dist_coeffs_6[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + dist_coeffs_6[int(6)] * r2_37, 2.0f * dist_coeffs_6[int(5)] * u_37 * v_37 + dist_coeffs_6[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + dist_coeffs_6[int(7)] * r2_37);
                    float2  _S697 = _S696 + make_float2 (dist_coeffs_6[int(8)] * _S696.x + dist_coeffs_6[int(9)] * _S696.y, 0.0f);
                    proj_points_0[int(3)] = make_float2 (fx_6 * _S697.x + cx_4, fy_6 * _S697.y + cy_4);
                    break;
                }
                bool all_valid_2 = _S604 & (!_S612);
                _S606 = &proj_points_0[int(4)];
                for(;;)
                {
                    float _S698 = _S623.p_0[int(4)].z;
                    proj_points_0[int(4)] = float2 {_S623.p_0[int(4)].x, _S623.p_0[int(4)].y} / make_float2 (_S698);
                    if(_S698 < 0.0f)
                    {
                        _S612 = true;
                    }
                    else
                    {
                        float u_38 = proj_points_0[int(4)].x;
                        float v_38 = proj_points_0[int(4)].y;
                        float _S699 = 0.0f * v_38;
                        float r2_38 = u_38 * u_38 + v_38 * v_38;
                        float s_diff_r2_38 = u_38 + u_38 + (_S699 + _S699);
                        float _S700 = dist_coeffs_6[int(2)] + r2_38 * dist_coeffs_6[int(3)];
                        float _S701 = dist_coeffs_6[int(1)] + r2_38 * _S700;
                        float _S702 = dist_coeffs_6[int(0)] + r2_38 * _S701;
                        float radial_10 = 1.0f + r2_38 * _S702;
                        float _S703 = 2.0f * dist_coeffs_6[int(4)];
                        float _S704 = _S703 * u_38;
                        float _S705 = 2.0f * u_38;
                        float _S706 = 2.0f * dist_coeffs_6[int(5)];
                        float _S707 = _S706 * u_38;
                        float _S708 = 2.0f * v_38;
                        float2  _S709 = make_float2 (1.0f, 0.0f) * make_float2 (radial_10) + make_float2 (s_diff_r2_38 * _S702 + (s_diff_r2_38 * _S701 + (s_diff_r2_38 * _S700 + s_diff_r2_38 * dist_coeffs_6[int(3)] * r2_38) * r2_38) * r2_38) * proj_points_0[int(4)] + make_float2 (_S703 * v_38 + 0.0f * _S704 + (s_diff_r2_38 + (_S705 + _S705)) * dist_coeffs_6[int(5)] + s_diff_r2_38 * dist_coeffs_6[int(6)], _S706 * v_38 + 0.0f * _S707 + (s_diff_r2_38 + (_S699 + 0.0f * _S708)) * dist_coeffs_6[int(4)] + s_diff_r2_38 * dist_coeffs_6[int(7)]);
                        float _S710 = 0.0f * u_38;
                        float s_diff_r2_39 = _S710 + _S710 + (v_38 + v_38);
                        float2  _S711 = make_float2 (0.0f, 1.0f) * make_float2 (radial_10) + make_float2 (s_diff_r2_39 * _S702 + (s_diff_r2_39 * _S701 + (s_diff_r2_39 * _S700 + s_diff_r2_39 * dist_coeffs_6[int(3)] * r2_38) * r2_38) * r2_38) * proj_points_0[int(4)] + make_float2 (0.0f * _S703 * v_38 + _S704 + (s_diff_r2_39 + (_S710 + 0.0f * _S705)) * dist_coeffs_6[int(5)] + s_diff_r2_39 * dist_coeffs_6[int(6)], 0.0f * _S706 * v_38 + _S707 + (s_diff_r2_39 + (_S708 + _S708)) * dist_coeffs_6[int(4)] + s_diff_r2_39 * dist_coeffs_6[int(7)]);
                        Matrix<float, 2, 2>  _S712 = transpose_0(makeMatrix<float, 2, 2> (_S709 + make_float2 (_S709.x * dist_coeffs_6[int(8)] + _S709.y * dist_coeffs_6[int(9)], 0.0f), _S711 + make_float2 (_S711.x * dist_coeffs_6[int(8)] + _S711.y * dist_coeffs_6[int(9)], 0.0f)));
                        _S612 = !((F32_min((determinant_0(_S712)), ((F32_min((_S712.rows[int(0)].x), (_S712.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S612)
                    {
                        break;
                    }
                    float u_39 = proj_points_0[int(4)].x;
                    float v_39 = proj_points_0[int(4)].y;
                    float r2_39 = u_39 * u_39 + v_39 * v_39;
                    float2  _S713 = proj_points_0[int(4)] * make_float2 (1.0f + r2_39 * (dist_coeffs_6[int(0)] + r2_39 * (dist_coeffs_6[int(1)] + r2_39 * (dist_coeffs_6[int(2)] + r2_39 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_39 * v_39 + dist_coeffs_6[int(5)] * (r2_39 + 2.0f * u_39 * u_39) + dist_coeffs_6[int(6)] * r2_39, 2.0f * dist_coeffs_6[int(5)] * u_39 * v_39 + dist_coeffs_6[int(4)] * (r2_39 + 2.0f * v_39 * v_39) + dist_coeffs_6[int(7)] * r2_39);
                    float2  _S714 = _S713 + make_float2 (dist_coeffs_6[int(8)] * _S713.x + dist_coeffs_6[int(9)] * _S713.y, 0.0f);
                    proj_points_0[int(4)] = make_float2 (fx_6 * _S714.x + cx_4, fy_6 * _S714.y + cy_4);
                    break;
                }
                bool all_valid_3 = all_valid_2 & (!_S612);
                for(;;)
                {
                    _S607 = &proj_points_0[int(5)];
                    for(;;)
                    {
                        float _S715 = _S623.p_0[int(5)].z;
                        proj_points_0[int(5)] = float2 {_S623.p_0[int(5)].x, _S623.p_0[int(5)].y} / make_float2 (_S715);
                        if(_S715 < 0.0f)
                        {
                            _S612 = true;
                        }
                        else
                        {
                            float u_40 = proj_points_0[int(5)].x;
                            float v_40 = proj_points_0[int(5)].y;
                            float _S716 = 0.0f * v_40;
                            float r2_40 = u_40 * u_40 + v_40 * v_40;
                            float s_diff_r2_40 = u_40 + u_40 + (_S716 + _S716);
                            float _S717 = dist_coeffs_6[int(2)] + r2_40 * dist_coeffs_6[int(3)];
                            float _S718 = dist_coeffs_6[int(1)] + r2_40 * _S717;
                            float _S719 = dist_coeffs_6[int(0)] + r2_40 * _S718;
                            float radial_11 = 1.0f + r2_40 * _S719;
                            float _S720 = 2.0f * dist_coeffs_6[int(4)];
                            float _S721 = _S720 * u_40;
                            float _S722 = 2.0f * u_40;
                            float _S723 = 2.0f * dist_coeffs_6[int(5)];
                            float _S724 = _S723 * u_40;
                            float _S725 = 2.0f * v_40;
                            float2  _S726 = make_float2 (1.0f, 0.0f) * make_float2 (radial_11) + make_float2 (s_diff_r2_40 * _S719 + (s_diff_r2_40 * _S718 + (s_diff_r2_40 * _S717 + s_diff_r2_40 * dist_coeffs_6[int(3)] * r2_40) * r2_40) * r2_40) * proj_points_0[int(5)] + make_float2 (_S720 * v_40 + 0.0f * _S721 + (s_diff_r2_40 + (_S722 + _S722)) * dist_coeffs_6[int(5)] + s_diff_r2_40 * dist_coeffs_6[int(6)], _S723 * v_40 + 0.0f * _S724 + (s_diff_r2_40 + (_S716 + 0.0f * _S725)) * dist_coeffs_6[int(4)] + s_diff_r2_40 * dist_coeffs_6[int(7)]);
                            float _S727 = 0.0f * u_40;
                            float s_diff_r2_41 = _S727 + _S727 + (v_40 + v_40);
                            float2  _S728 = make_float2 (0.0f, 1.0f) * make_float2 (radial_11) + make_float2 (s_diff_r2_41 * _S719 + (s_diff_r2_41 * _S718 + (s_diff_r2_41 * _S717 + s_diff_r2_41 * dist_coeffs_6[int(3)] * r2_40) * r2_40) * r2_40) * proj_points_0[int(5)] + make_float2 (0.0f * _S720 * v_40 + _S721 + (s_diff_r2_41 + (_S727 + 0.0f * _S722)) * dist_coeffs_6[int(5)] + s_diff_r2_41 * dist_coeffs_6[int(6)], 0.0f * _S723 * v_40 + _S724 + (s_diff_r2_41 + (_S725 + _S725)) * dist_coeffs_6[int(4)] + s_diff_r2_41 * dist_coeffs_6[int(7)]);
                            Matrix<float, 2, 2>  _S729 = transpose_0(makeMatrix<float, 2, 2> (_S726 + make_float2 (_S726.x * dist_coeffs_6[int(8)] + _S726.y * dist_coeffs_6[int(9)], 0.0f), _S728 + make_float2 (_S728.x * dist_coeffs_6[int(8)] + _S728.y * dist_coeffs_6[int(9)], 0.0f)));
                            _S612 = !((F32_min((determinant_0(_S729)), ((F32_min((_S729.rows[int(0)].x), (_S729.rows[int(1)].y)))))) > 0.0f);
                        }
                        if(_S612)
                        {
                            break;
                        }
                        float u_41 = proj_points_0[int(5)].x;
                        float v_41 = proj_points_0[int(5)].y;
                        float r2_41 = u_41 * u_41 + v_41 * v_41;
                        float2  _S730 = proj_points_0[int(5)] * make_float2 (1.0f + r2_41 * (dist_coeffs_6[int(0)] + r2_41 * (dist_coeffs_6[int(1)] + r2_41 * (dist_coeffs_6[int(2)] + r2_41 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_41 * v_41 + dist_coeffs_6[int(5)] * (r2_41 + 2.0f * u_41 * u_41) + dist_coeffs_6[int(6)] * r2_41, 2.0f * dist_coeffs_6[int(5)] * u_41 * v_41 + dist_coeffs_6[int(4)] * (r2_41 + 2.0f * v_41 * v_41) + dist_coeffs_6[int(7)] * r2_41);
                        float2  _S731 = _S730 + make_float2 (dist_coeffs_6[int(8)] * _S730.x + dist_coeffs_6[int(9)] * _S730.y, 0.0f);
                        proj_points_0[int(5)] = make_float2 (fx_6 * _S731.x + cx_4, fy_6 * _S731.y + cy_4);
                        break;
                    }
                    _S608 = all_valid_3 & (!_S612);
                    break;
                }
                _S609 = &proj_points_0[int(6)];
                for(;;)
                {
                    float _S732 = _S623.p_0[int(6)].z;
                    proj_points_0[int(6)] = float2 {_S623.p_0[int(6)].x, _S623.p_0[int(6)].y} / make_float2 (_S732);
                    if(_S732 < 0.0f)
                    {
                        _S612 = true;
                    }
                    else
                    {
                        float u_42 = proj_points_0[int(6)].x;
                        float v_42 = proj_points_0[int(6)].y;
                        float _S733 = 0.0f * v_42;
                        float r2_42 = u_42 * u_42 + v_42 * v_42;
                        float s_diff_r2_42 = u_42 + u_42 + (_S733 + _S733);
                        float _S734 = dist_coeffs_6[int(2)] + r2_42 * dist_coeffs_6[int(3)];
                        float _S735 = dist_coeffs_6[int(1)] + r2_42 * _S734;
                        float _S736 = dist_coeffs_6[int(0)] + r2_42 * _S735;
                        float radial_12 = 1.0f + r2_42 * _S736;
                        float _S737 = 2.0f * dist_coeffs_6[int(4)];
                        float _S738 = _S737 * u_42;
                        float _S739 = 2.0f * u_42;
                        float _S740 = 2.0f * dist_coeffs_6[int(5)];
                        float _S741 = _S740 * u_42;
                        float _S742 = 2.0f * v_42;
                        float2  _S743 = make_float2 (1.0f, 0.0f) * make_float2 (radial_12) + make_float2 (s_diff_r2_42 * _S736 + (s_diff_r2_42 * _S735 + (s_diff_r2_42 * _S734 + s_diff_r2_42 * dist_coeffs_6[int(3)] * r2_42) * r2_42) * r2_42) * proj_points_0[int(6)] + make_float2 (_S737 * v_42 + 0.0f * _S738 + (s_diff_r2_42 + (_S739 + _S739)) * dist_coeffs_6[int(5)] + s_diff_r2_42 * dist_coeffs_6[int(6)], _S740 * v_42 + 0.0f * _S741 + (s_diff_r2_42 + (_S733 + 0.0f * _S742)) * dist_coeffs_6[int(4)] + s_diff_r2_42 * dist_coeffs_6[int(7)]);
                        float _S744 = 0.0f * u_42;
                        float s_diff_r2_43 = _S744 + _S744 + (v_42 + v_42);
                        float2  _S745 = make_float2 (0.0f, 1.0f) * make_float2 (radial_12) + make_float2 (s_diff_r2_43 * _S736 + (s_diff_r2_43 * _S735 + (s_diff_r2_43 * _S734 + s_diff_r2_43 * dist_coeffs_6[int(3)] * r2_42) * r2_42) * r2_42) * proj_points_0[int(6)] + make_float2 (0.0f * _S737 * v_42 + _S738 + (s_diff_r2_43 + (_S744 + 0.0f * _S739)) * dist_coeffs_6[int(5)] + s_diff_r2_43 * dist_coeffs_6[int(6)], 0.0f * _S740 * v_42 + _S741 + (s_diff_r2_43 + (_S742 + _S742)) * dist_coeffs_6[int(4)] + s_diff_r2_43 * dist_coeffs_6[int(7)]);
                        Matrix<float, 2, 2>  _S746 = transpose_0(makeMatrix<float, 2, 2> (_S743 + make_float2 (_S743.x * dist_coeffs_6[int(8)] + _S743.y * dist_coeffs_6[int(9)], 0.0f), _S745 + make_float2 (_S745.x * dist_coeffs_6[int(8)] + _S745.y * dist_coeffs_6[int(9)], 0.0f)));
                        _S612 = !((F32_min((determinant_0(_S746)), ((F32_min((_S746.rows[int(0)].x), (_S746.rows[int(1)].y)))))) > 0.0f);
                    }
                    if(_S612)
                    {
                        break;
                    }
                    float u_43 = proj_points_0[int(6)].x;
                    float v_43 = proj_points_0[int(6)].y;
                    float r2_43 = u_43 * u_43 + v_43 * v_43;
                    float2  _S747 = proj_points_0[int(6)] * make_float2 (1.0f + r2_43 * (dist_coeffs_6[int(0)] + r2_43 * (dist_coeffs_6[int(1)] + r2_43 * (dist_coeffs_6[int(2)] + r2_43 * dist_coeffs_6[int(3)])))) + make_float2 (2.0f * dist_coeffs_6[int(4)] * u_43 * v_43 + dist_coeffs_6[int(5)] * (r2_43 + 2.0f * u_43 * u_43) + dist_coeffs_6[int(6)] * r2_43, 2.0f * dist_coeffs_6[int(5)] * u_43 * v_43 + dist_coeffs_6[int(4)] * (r2_43 + 2.0f * v_43 * v_43) + dist_coeffs_6[int(7)] * r2_43);
                    float2  _S748 = _S747 + make_float2 (dist_coeffs_6[int(8)] * _S747.x + dist_coeffs_6[int(9)] * _S747.y, 0.0f);
                    proj_points_0[int(6)] = make_float2 (fx_6 * _S748.x + cx_4, fy_6 * _S748.y + cy_4);
                    break;
                }
                _S610 = _S608 & (!_S612);
                break;
            }
            if(!_S610)
            {
                _S612 = false;
                break;
            }
            float2  _S749 = *mean2d_6 + make_float2 (_S623.w_mean_0[int(0)]) * *_S601 + make_float2 (_S623.w_mean_0[int(1)]) * *_S602 + make_float2 (_S623.w_mean_0[int(2)]) * *_S603 + make_float2 (_S623.w_mean_0[int(3)]) * *_S605 + make_float2 (_S623.w_mean_0[int(4)]) * *_S606 + make_float2 (_S623.w_mean_0[int(5)]) * *_S607 + make_float2 (_S623.w_mean_0[int(6)]) * *_S609;
            *mean2d_6 = _S749;
            float _S750 = - _S628;
            float _S751 = - _S629;
            float2  _S752 = make_float2 (clamp_0(_S749.x, _S750, lim_x_pos_2), clamp_0(_S749.y, _S751, lim_y_pos_0));
            float2  d_0 = make_float2 (clamp_0((*_S601).x, _S750, lim_x_pos_2), clamp_0((*_S601).y, _S751, lim_y_pos_0)) - _S752;
            float _S753 = d_0.x;
            float _S754 = d_0.y;
            float _S755 = _S753 * _S754;
            float2  d_1 = make_float2 (clamp_0((*_S602).x, _S750, lim_x_pos_2), clamp_0((*_S602).y, _S751, lim_y_pos_0)) - _S752;
            float _S756 = d_1.x;
            float _S757 = d_1.y;
            float _S758 = _S756 * _S757;
            float2  d_2 = make_float2 (clamp_0((*_S603).x, _S750, lim_x_pos_2), clamp_0((*_S603).y, _S751, lim_y_pos_0)) - _S752;
            float _S759 = d_2.x;
            float _S760 = d_2.y;
            float _S761 = _S759 * _S760;
            float2  d_3 = make_float2 (clamp_0((*_S605).x, _S750, lim_x_pos_2), clamp_0((*_S605).y, _S751, lim_y_pos_0)) - _S752;
            float _S762 = d_3.x;
            float _S763 = d_3.y;
            float _S764 = _S762 * _S763;
            float2  d_4 = make_float2 (clamp_0((*_S606).x, _S750, lim_x_pos_2), clamp_0((*_S606).y, _S751, lim_y_pos_0)) - _S752;
            float _S765 = d_4.x;
            float _S766 = d_4.y;
            float _S767 = _S765 * _S766;
            float2  d_5 = make_float2 (clamp_0((*_S607).x, _S750, lim_x_pos_2), clamp_0((*_S607).y, _S751, lim_y_pos_0)) - _S752;
            float _S768 = d_5.x;
            float _S769 = d_5.y;
            float _S770 = _S768 * _S769;
            float2  d_6 = make_float2 (clamp_0((*_S609).x, _S750, lim_x_pos_2), clamp_0((*_S609).y, _S751, lim_y_pos_0)) - _S752;
            float _S771 = d_6.x;
            float _S772 = d_6.y;
            float _S773 = _S771 * _S772;
            covar2d_3 = covar2d_3 + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S753 * _S753, _S755, _S755, _S754 * _S754) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S756 * _S756, _S758, _S758, _S757 * _S757) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S759 * _S759, _S761, _S761, _S760 * _S760) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S762 * _S762, _S764, _S764, _S763 * _S763) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S765 * _S765, _S767, _S767, _S766 * _S766) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S768 * _S768, _S770, _S770, _S769 * _S769) + makeMatrix<float, 2, 2> (_S623.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S771 * _S771, _S773, _S773, _S772 * _S772);
            _S612 = true;
            break;
        }
        if(!(true & _S612))
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
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
        float _S774 = *&(((&covar2d_3)->rows + (int(0)))->x) + eps2d_3;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S774;
        float _S775 = *&(((&covar2d_3)->rows + (int(1)))->y) + eps2d_3;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S775;
        float det_blur_3 = _S774 * _S775 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        float invdet_5 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
        Matrix<float, 2, 2>  _S776 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_5, - covar2d_3.rows[int(0)].y * invdet_5, - covar2d_3.rows[int(1)].x * invdet_5, covar2d_3.rows[int(0)].x * invdet_5);
        if(antialiased_3)
        {
            *opacity_3 = *opacity_3 * compensation_3;
        }
        if((*opacity_3) < 0.00392156885936856f)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        float _S777 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
        float radius_x_3 = _S777 * (F32_sqrt((covar2d_3[int(0)].x)));
        float radius_y_3 = _S777 * (F32_sqrt((covar2d_3[int(1)].y)));
        float _S778 = (*mean2d_6).x - radius_x_3;
        float _S779 = (*mean2d_6).x + radius_x_3;
        float _S780 = (*mean2d_6).y - radius_y_3;
        float _S781 = (*mean2d_6).y + radius_y_3;
        if(_S779 <= 0.0f)
        {
            _S612 = true;
        }
        else
        {
            _S612 = _S778 >= _S599;
        }
        if(_S612)
        {
            _S612 = true;
        }
        else
        {
            _S612 = _S781 <= 0.0f;
        }
        if(_S612)
        {
            _S612 = true;
        }
        else
        {
            _S612 = _S780 >= _S600;
        }
        if(_S612)
        {
            *aabb_xyxy_3 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_3 = make_float4 (_S778, _S780, _S779, _S781);
        *sorting_depth_3 = _S611;
        *conic_3 = make_float3 (_S776.rows[int(0)].x, _S776.rows[int(0)].y, _S776.rows[int(1)].y);
        *radius_4 = view_radius_3dgs_0(mean_4, scale_3, in_opacity_3, - mul_6(transpose_3(R_3), t_3));
        break;
    }
    return;
}

inline __device__ void projection_3dgut_fisheye(bool antialiased_4, float3  mean_5, float4  quat_4, float3  scale_4, float in_opacity_4, Matrix<float, 3, 3>  R_4, float3  t_4, float fx_7, float fy_7, float cx_5, float cy_5, FixedArray<float, 10>  dist_coeffs_7, uint image_width_4, uint image_height_4, float4  * aabb_xyxy_4, float * sorting_depth_4, float * radius_5, float2  * mean2d_7, float * depth_4, float3  * conic_4, float * opacity_4)
{
    float2  * _S782;
    float _S783;
    float2  _S784;
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
    float2  _S797;
    float _S798;
    float _S799;
    bool _S800;
    float2  * _S801;
    float _S802;
    bool _S803;
    float2  * _S804;
    float _S805;
    bool _S806;
    bool _S807;
    float2  * _S808;
    float _S809;
    bool _S810;
    float2  * _S811;
    float _S812;
    bool _S813;
    float2  * _S814;
    float _S815;
    bool _S816;
    bool _S817;
    float2  * _S818;
    float _S819;
    bool _S820;
    bool _S821;
    for(;;)
    {
        float3  mean_c_3 = mul_6(R_4, mean_5) + t_4;
        float _S822 = length_1(mean_c_3);
        float _S823 = mean_c_3.z;
        *depth_4 = _S823;
        if(_S822 <= 0.0f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_4;
        *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
        bool _S824;
        float3  _S825 = exp_0(scale_4);
        float4  _S826 = normalize_0(quat_4);
        float x_21 = _S826.y;
        float x2_4 = x_21 * x_21;
        float y2_4 = _S826.z * _S826.z;
        float z2_4 = _S826.w * _S826.w;
        float xy_4 = _S826.y * _S826.z;
        float xz_4 = _S826.y * _S826.w;
        float yz_4 = _S826.z * _S826.w;
        float wx_4 = _S826.x * _S826.y;
        float wy_4 = _S826.x * _S826.z;
        float wz_4 = _S826.x * _S826.w;
        Matrix<float, 3, 3>  _S827 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_5;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S828 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S828 * _S825.x) * _S827.rows[0U];
        float3  _S829 = mean_5 + delta_3;
        float3  _S830 = mean_5 - delta_3;
        float3  delta_4 = make_float3 (_S828 * _S825.y) * _S827.rows[1U];
        float3  _S831 = mean_5 + delta_4;
        float3  _S832 = mean_5 - delta_4;
        float3  delta_5 = make_float3 (_S828 * _S825.z) * _S827.rows[2U];
        float3  _S833 = mean_5 + delta_5;
        float3  _S834 = mean_5 - delta_5;
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
        (&ret_1)->p_0[0U] = mul_6(R_4, (&ret_1)->p_0[0U]) + t_4;
        (&ret_1)->p_0[1U] = mul_6(R_4, _S829) + t_4;
        (&ret_1)->p_0[2U] = mul_6(R_4, _S831) + t_4;
        (&ret_1)->p_0[3U] = mul_6(R_4, _S833) + t_4;
        (&ret_1)->p_0[4U] = mul_6(R_4, _S830) + t_4;
        (&ret_1)->p_0[5U] = mul_6(R_4, _S832) + t_4;
        (&ret_1)->p_0[6U] = mul_6(R_4, _S834) + t_4;
        SigmaPoints_0 _S835 = ret_1;
        for(;;)
        {
            int2  _S836 = make_int2 (int(0));
            float2  _S837 = make_float2 ((float)_S836.x, (float)_S836.y);
            *mean2d_7 = _S837;
            covar2d_4 = makeMatrix<float, 2, 2> (0.0f);
            FixedArray<float2 , 7>  proj_points_1;
            for(;;)
            {
                float k_4;
                _S782 = &proj_points_1[int(0)];
                for(;;)
                {
                    float2  _S838 = float2 {_S835.p_0[int(0)].x, _S835.p_0[int(0)].y};
                    float r_11 = length_0(_S838);
                    float _S839 = _S835.p_0[int(0)].z;
                    _S783 = _S839;
                    float theta_4 = (F32_atan2((r_11), (_S839)));
                    if(theta_4 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S839;
                    }
                    else
                    {
                        k_4 = theta_4 / r_11;
                    }
                    float2  _S840 = _S838 * make_float2 (k_4);
                    proj_points_1[int(0)] = _S840;
                    float2  _S841 = make_float2 (1.0f, 0.0f);
                    _S784 = _S841;
                    _S785 = dist_coeffs_7[int(0)];
                    _S786 = dist_coeffs_7[int(1)];
                    _S787 = dist_coeffs_7[int(2)];
                    _S788 = dist_coeffs_7[int(3)];
                    _S789 = dist_coeffs_7[int(4)];
                    _S790 = dist_coeffs_7[int(5)];
                    _S791 = dist_coeffs_7[int(6)];
                    _S792 = dist_coeffs_7[int(7)];
                    _S793 = dist_coeffs_7[int(8)];
                    _S794 = dist_coeffs_7[int(9)];
                    float u_44 = _S840.x;
                    float v_44 = _S840.y;
                    float _S842 = 0.0f * v_44;
                    float r2_44 = u_44 * u_44 + v_44 * v_44;
                    float s_diff_r2_44 = u_44 + u_44 + (_S842 + _S842);
                    float _S843 = dist_coeffs_7[int(2)] + r2_44 * dist_coeffs_7[int(3)];
                    float _S844 = dist_coeffs_7[int(1)] + r2_44 * _S843;
                    float _S845 = dist_coeffs_7[int(0)] + r2_44 * _S844;
                    float _S846 = s_diff_r2_44 * _S845 + (s_diff_r2_44 * _S844 + (s_diff_r2_44 * _S843 + s_diff_r2_44 * dist_coeffs_7[int(3)] * r2_44) * r2_44) * r2_44;
                    float radial_13 = 1.0f + r2_44 * _S845;
                    float _S847 = 2.0f * dist_coeffs_7[int(4)];
                    _S795 = _S847;
                    float _S848 = _S847 * u_44;
                    float _S849 = 2.0f * u_44;
                    float s_diff_du_4 = _S847 * v_44 + 0.0f * _S848 + (s_diff_r2_44 + (_S849 + _S849)) * dist_coeffs_7[int(5)] + s_diff_r2_44 * dist_coeffs_7[int(6)];
                    float _S850 = 2.0f * dist_coeffs_7[int(5)];
                    _S796 = _S850;
                    float _S851 = _S850 * u_44;
                    float _S852 = 2.0f * v_44;
                    float2  _S853 = _S841 * make_float2 (radial_13) + make_float2 (_S846) * _S840 + make_float2 (s_diff_du_4, _S850 * v_44 + 0.0f * _S851 + (s_diff_r2_44 + (_S842 + 0.0f * _S852)) * dist_coeffs_7[int(4)] + s_diff_r2_44 * dist_coeffs_7[int(7)]);
                    float2  _S854 = _S853 + make_float2 (_S853.x * dist_coeffs_7[int(8)] + _S853.y * dist_coeffs_7[int(9)], 0.0f);
                    float2  _S855 = make_float2 (0.0f, 1.0f);
                    _S797 = _S855;
                    float _S856 = 0.0f * u_44;
                    float s_diff_r2_45 = _S856 + _S856 + (v_44 + v_44);
                    float _S857 = s_diff_r2_45 * _S845 + (s_diff_r2_45 * _S844 + (s_diff_r2_45 * _S843 + s_diff_r2_45 * dist_coeffs_7[int(3)] * r2_44) * r2_44) * r2_44;
                    float _S858 = 0.0f * _S847;
                    _S798 = _S858;
                    float s_diff_du_5 = _S858 * v_44 + _S848 + (s_diff_r2_45 + (_S856 + 0.0f * _S849)) * dist_coeffs_7[int(5)] + s_diff_r2_45 * dist_coeffs_7[int(6)];
                    float _S859 = 0.0f * _S850;
                    _S799 = _S859;
                    float2  _S860 = _S855 * make_float2 (radial_13) + make_float2 (_S857) * _S840 + make_float2 (s_diff_du_5, _S859 * v_44 + _S851 + (s_diff_r2_45 + (_S852 + _S852)) * dist_coeffs_7[int(4)] + s_diff_r2_45 * dist_coeffs_7[int(7)]);
                    Matrix<float, 2, 2>  _S861 = transpose_0(makeMatrix<float, 2, 2> (_S854, _S860 + make_float2 (_S860.x * dist_coeffs_7[int(8)] + _S860.y * dist_coeffs_7[int(9)], 0.0f)));
                    bool _S862 = !((F32_min((determinant_0(_S861)), ((F32_min((_S861.rows[int(0)].x), (_S861.rows[int(1)].y)))))) > 0.0f);
                    _S800 = _S862;
                    if(_S862)
                    {
                        break;
                    }
                    float u_45 = proj_points_1[int(0)].x;
                    float v_45 = proj_points_1[int(0)].y;
                    float r2_45 = u_45 * u_45 + v_45 * v_45;
                    float2  _S863 = proj_points_1[int(0)] * make_float2 (1.0f + r2_45 * (dist_coeffs_7[int(0)] + r2_45 * (dist_coeffs_7[int(1)] + r2_45 * (dist_coeffs_7[int(2)] + r2_45 * dist_coeffs_7[int(3)])))) + make_float2 (_S847 * u_45 * v_45 + dist_coeffs_7[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + dist_coeffs_7[int(6)] * r2_45, _S850 * u_45 * v_45 + dist_coeffs_7[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + dist_coeffs_7[int(7)] * r2_45);
                    float2  _S864 = _S863 + make_float2 (dist_coeffs_7[int(8)] * _S863.x + dist_coeffs_7[int(9)] * _S863.y, 0.0f);
                    proj_points_1[int(0)] = make_float2 (fx_7 * _S864.x + cx_5, fy_7 * _S864.y + cy_5);
                    break;
                }
                bool all_valid_4 = true & (!_S800);
                _S801 = &proj_points_1[int(1)];
                for(;;)
                {
                    float2  _S865 = float2 {_S835.p_0[int(1)].x, _S835.p_0[int(1)].y};
                    float r_12 = length_0(_S865);
                    float _S866 = _S835.p_0[int(1)].z;
                    _S802 = _S866;
                    float theta_5 = (F32_atan2((r_12), (_S866)));
                    if(theta_5 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_5 * theta_5 / 3.0f) / _S866;
                    }
                    else
                    {
                        k_4 = theta_5 / r_12;
                    }
                    float2  _S867 = _S865 * make_float2 (k_4);
                    proj_points_1[int(1)] = _S867;
                    float u_46 = _S867.x;
                    float v_46 = _S867.y;
                    float _S868 = 0.0f * v_46;
                    float r2_46 = u_46 * u_46 + v_46 * v_46;
                    float s_diff_r2_46 = u_46 + u_46 + (_S868 + _S868);
                    float _S869 = _S787 + r2_46 * _S788;
                    float _S870 = _S786 + r2_46 * _S869;
                    float _S871 = _S785 + r2_46 * _S870;
                    float radial_14 = 1.0f + r2_46 * _S871;
                    float _S872 = _S795 * u_46;
                    float _S873 = 2.0f * u_46;
                    float _S874 = _S796 * u_46;
                    float _S875 = 2.0f * v_46;
                    float2  _S876 = _S784 * make_float2 (radial_14) + make_float2 (s_diff_r2_46 * _S871 + (s_diff_r2_46 * _S870 + (s_diff_r2_46 * _S869 + s_diff_r2_46 * _S788 * r2_46) * r2_46) * r2_46) * _S867 + make_float2 (_S795 * v_46 + 0.0f * _S872 + (s_diff_r2_46 + (_S873 + _S873)) * _S790 + s_diff_r2_46 * _S791, _S796 * v_46 + 0.0f * _S874 + (s_diff_r2_46 + (_S868 + 0.0f * _S875)) * _S789 + s_diff_r2_46 * _S792);
                    float _S877 = 0.0f * u_46;
                    float s_diff_r2_47 = _S877 + _S877 + (v_46 + v_46);
                    float2  _S878 = _S797 * make_float2 (radial_14) + make_float2 (s_diff_r2_47 * _S871 + (s_diff_r2_47 * _S870 + (s_diff_r2_47 * _S869 + s_diff_r2_47 * _S788 * r2_46) * r2_46) * r2_46) * _S867 + make_float2 (_S798 * v_46 + _S872 + (s_diff_r2_47 + (_S877 + 0.0f * _S873)) * _S790 + s_diff_r2_47 * _S791, _S799 * v_46 + _S874 + (s_diff_r2_47 + (_S875 + _S875)) * _S789 + s_diff_r2_47 * _S792);
                    Matrix<float, 2, 2>  _S879 = transpose_0(makeMatrix<float, 2, 2> (_S876 + make_float2 (_S876.x * _S793 + _S876.y * _S794, 0.0f), _S878 + make_float2 (_S878.x * _S793 + _S878.y * _S794, 0.0f)));
                    bool _S880 = !((F32_min((determinant_0(_S879)), ((F32_min((_S879.rows[int(0)].x), (_S879.rows[int(1)].y)))))) > 0.0f);
                    _S803 = _S880;
                    if(_S880)
                    {
                        break;
                    }
                    float u_47 = proj_points_1[int(1)].x;
                    float v_47 = proj_points_1[int(1)].y;
                    float r2_47 = u_47 * u_47 + v_47 * v_47;
                    float2  _S881 = proj_points_1[int(1)] * make_float2 (1.0f + r2_47 * (_S785 + r2_47 * (_S786 + r2_47 * (_S787 + r2_47 * _S788)))) + make_float2 (_S795 * u_47 * v_47 + _S790 * (r2_47 + 2.0f * u_47 * u_47) + _S791 * r2_47, _S796 * u_47 * v_47 + _S789 * (r2_47 + 2.0f * v_47 * v_47) + _S792 * r2_47);
                    float2  _S882 = _S881 + make_float2 (_S793 * _S881.x + _S794 * _S881.y, 0.0f);
                    proj_points_1[int(1)] = make_float2 (fx_7 * _S882.x + cx_5, fy_7 * _S882.y + cy_5);
                    break;
                }
                bool all_valid_5 = all_valid_4 & (!_S803);
                for(;;)
                {
                    _S804 = &proj_points_1[int(2)];
                    for(;;)
                    {
                        float2  _S883 = float2 {_S835.p_0[int(2)].x, _S835.p_0[int(2)].y};
                        float r_13 = length_0(_S883);
                        float _S884 = _S835.p_0[int(2)].z;
                        _S805 = _S884;
                        float theta_6 = (F32_atan2((r_13), (_S884)));
                        if(theta_6 < 0.00100000004749745f)
                        {
                            k_4 = (1.0f - theta_6 * theta_6 / 3.0f) / _S884;
                        }
                        else
                        {
                            k_4 = theta_6 / r_13;
                        }
                        float2  _S885 = _S883 * make_float2 (k_4);
                        proj_points_1[int(2)] = _S885;
                        float u_48 = _S885.x;
                        float v_48 = _S885.y;
                        float _S886 = 0.0f * v_48;
                        float r2_48 = u_48 * u_48 + v_48 * v_48;
                        float s_diff_r2_48 = u_48 + u_48 + (_S886 + _S886);
                        float _S887 = _S787 + r2_48 * _S788;
                        float _S888 = _S786 + r2_48 * _S887;
                        float _S889 = _S785 + r2_48 * _S888;
                        float radial_15 = 1.0f + r2_48 * _S889;
                        float _S890 = _S795 * u_48;
                        float _S891 = 2.0f * u_48;
                        float _S892 = _S796 * u_48;
                        float _S893 = 2.0f * v_48;
                        float2  _S894 = _S784 * make_float2 (radial_15) + make_float2 (s_diff_r2_48 * _S889 + (s_diff_r2_48 * _S888 + (s_diff_r2_48 * _S887 + s_diff_r2_48 * _S788 * r2_48) * r2_48) * r2_48) * _S885 + make_float2 (_S795 * v_48 + 0.0f * _S890 + (s_diff_r2_48 + (_S891 + _S891)) * _S790 + s_diff_r2_48 * _S791, _S796 * v_48 + 0.0f * _S892 + (s_diff_r2_48 + (_S886 + 0.0f * _S893)) * _S789 + s_diff_r2_48 * _S792);
                        float _S895 = 0.0f * u_48;
                        float s_diff_r2_49 = _S895 + _S895 + (v_48 + v_48);
                        float2  _S896 = _S797 * make_float2 (radial_15) + make_float2 (s_diff_r2_49 * _S889 + (s_diff_r2_49 * _S888 + (s_diff_r2_49 * _S887 + s_diff_r2_49 * _S788 * r2_48) * r2_48) * r2_48) * _S885 + make_float2 (_S798 * v_48 + _S890 + (s_diff_r2_49 + (_S895 + 0.0f * _S891)) * _S790 + s_diff_r2_49 * _S791, _S799 * v_48 + _S892 + (s_diff_r2_49 + (_S893 + _S893)) * _S789 + s_diff_r2_49 * _S792);
                        Matrix<float, 2, 2>  _S897 = transpose_0(makeMatrix<float, 2, 2> (_S894 + make_float2 (_S894.x * _S793 + _S894.y * _S794, 0.0f), _S896 + make_float2 (_S896.x * _S793 + _S896.y * _S794, 0.0f)));
                        bool _S898 = !((F32_min((determinant_0(_S897)), ((F32_min((_S897.rows[int(0)].x), (_S897.rows[int(1)].y)))))) > 0.0f);
                        _S806 = _S898;
                        if(_S898)
                        {
                            break;
                        }
                        float u_49 = proj_points_1[int(2)].x;
                        float v_49 = proj_points_1[int(2)].y;
                        float r2_49 = u_49 * u_49 + v_49 * v_49;
                        float2  _S899 = proj_points_1[int(2)] * make_float2 (1.0f + r2_49 * (_S785 + r2_49 * (_S786 + r2_49 * (_S787 + r2_49 * _S788)))) + make_float2 (_S795 * u_49 * v_49 + _S790 * (r2_49 + 2.0f * u_49 * u_49) + _S791 * r2_49, _S796 * u_49 * v_49 + _S789 * (r2_49 + 2.0f * v_49 * v_49) + _S792 * r2_49);
                        float2  _S900 = _S899 + make_float2 (_S793 * _S899.x + _S794 * _S899.y, 0.0f);
                        proj_points_1[int(2)] = make_float2 (fx_7 * _S900.x + cx_5, fy_7 * _S900.y + cy_5);
                        break;
                    }
                    _S807 = all_valid_5 & (!_S806);
                    break;
                }
                _S808 = &proj_points_1[int(3)];
                for(;;)
                {
                    float2  _S901 = float2 {_S835.p_0[int(3)].x, _S835.p_0[int(3)].y};
                    float r_14 = length_0(_S901);
                    float _S902 = _S835.p_0[int(3)].z;
                    _S809 = _S902;
                    float theta_7 = (F32_atan2((r_14), (_S902)));
                    if(theta_7 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_7 * theta_7 / 3.0f) / _S902;
                    }
                    else
                    {
                        k_4 = theta_7 / r_14;
                    }
                    float2  _S903 = _S901 * make_float2 (k_4);
                    proj_points_1[int(3)] = _S903;
                    float u_50 = _S903.x;
                    float v_50 = _S903.y;
                    float _S904 = 0.0f * v_50;
                    float r2_50 = u_50 * u_50 + v_50 * v_50;
                    float s_diff_r2_50 = u_50 + u_50 + (_S904 + _S904);
                    float _S905 = _S787 + r2_50 * _S788;
                    float _S906 = _S786 + r2_50 * _S905;
                    float _S907 = _S785 + r2_50 * _S906;
                    float radial_16 = 1.0f + r2_50 * _S907;
                    float _S908 = _S795 * u_50;
                    float _S909 = 2.0f * u_50;
                    float _S910 = _S796 * u_50;
                    float _S911 = 2.0f * v_50;
                    float2  _S912 = _S784 * make_float2 (radial_16) + make_float2 (s_diff_r2_50 * _S907 + (s_diff_r2_50 * _S906 + (s_diff_r2_50 * _S905 + s_diff_r2_50 * _S788 * r2_50) * r2_50) * r2_50) * _S903 + make_float2 (_S795 * v_50 + 0.0f * _S908 + (s_diff_r2_50 + (_S909 + _S909)) * _S790 + s_diff_r2_50 * _S791, _S796 * v_50 + 0.0f * _S910 + (s_diff_r2_50 + (_S904 + 0.0f * _S911)) * _S789 + s_diff_r2_50 * _S792);
                    float _S913 = 0.0f * u_50;
                    float s_diff_r2_51 = _S913 + _S913 + (v_50 + v_50);
                    float2  _S914 = _S797 * make_float2 (radial_16) + make_float2 (s_diff_r2_51 * _S907 + (s_diff_r2_51 * _S906 + (s_diff_r2_51 * _S905 + s_diff_r2_51 * _S788 * r2_50) * r2_50) * r2_50) * _S903 + make_float2 (_S798 * v_50 + _S908 + (s_diff_r2_51 + (_S913 + 0.0f * _S909)) * _S790 + s_diff_r2_51 * _S791, _S799 * v_50 + _S910 + (s_diff_r2_51 + (_S911 + _S911)) * _S789 + s_diff_r2_51 * _S792);
                    Matrix<float, 2, 2>  _S915 = transpose_0(makeMatrix<float, 2, 2> (_S912 + make_float2 (_S912.x * _S793 + _S912.y * _S794, 0.0f), _S914 + make_float2 (_S914.x * _S793 + _S914.y * _S794, 0.0f)));
                    bool _S916 = !((F32_min((determinant_0(_S915)), ((F32_min((_S915.rows[int(0)].x), (_S915.rows[int(1)].y)))))) > 0.0f);
                    _S810 = _S916;
                    if(_S916)
                    {
                        break;
                    }
                    float u_51 = proj_points_1[int(3)].x;
                    float v_51 = proj_points_1[int(3)].y;
                    float r2_51 = u_51 * u_51 + v_51 * v_51;
                    float2  _S917 = proj_points_1[int(3)] * make_float2 (1.0f + r2_51 * (_S785 + r2_51 * (_S786 + r2_51 * (_S787 + r2_51 * _S788)))) + make_float2 (_S795 * u_51 * v_51 + _S790 * (r2_51 + 2.0f * u_51 * u_51) + _S791 * r2_51, _S796 * u_51 * v_51 + _S789 * (r2_51 + 2.0f * v_51 * v_51) + _S792 * r2_51);
                    float2  _S918 = _S917 + make_float2 (_S793 * _S917.x + _S794 * _S917.y, 0.0f);
                    proj_points_1[int(3)] = make_float2 (fx_7 * _S918.x + cx_5, fy_7 * _S918.y + cy_5);
                    break;
                }
                bool all_valid_6 = _S807 & (!_S810);
                _S811 = &proj_points_1[int(4)];
                for(;;)
                {
                    float2  _S919 = float2 {_S835.p_0[int(4)].x, _S835.p_0[int(4)].y};
                    float r_15 = length_0(_S919);
                    float _S920 = _S835.p_0[int(4)].z;
                    _S812 = _S920;
                    float theta_8 = (F32_atan2((r_15), (_S920)));
                    if(theta_8 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_8 * theta_8 / 3.0f) / _S920;
                    }
                    else
                    {
                        k_4 = theta_8 / r_15;
                    }
                    float2  _S921 = _S919 * make_float2 (k_4);
                    proj_points_1[int(4)] = _S921;
                    float u_52 = _S921.x;
                    float v_52 = _S921.y;
                    float _S922 = 0.0f * v_52;
                    float r2_52 = u_52 * u_52 + v_52 * v_52;
                    float s_diff_r2_52 = u_52 + u_52 + (_S922 + _S922);
                    float _S923 = _S787 + r2_52 * _S788;
                    float _S924 = _S786 + r2_52 * _S923;
                    float _S925 = _S785 + r2_52 * _S924;
                    float radial_17 = 1.0f + r2_52 * _S925;
                    float _S926 = _S795 * u_52;
                    float _S927 = 2.0f * u_52;
                    float _S928 = _S796 * u_52;
                    float _S929 = 2.0f * v_52;
                    float2  _S930 = _S784 * make_float2 (radial_17) + make_float2 (s_diff_r2_52 * _S925 + (s_diff_r2_52 * _S924 + (s_diff_r2_52 * _S923 + s_diff_r2_52 * _S788 * r2_52) * r2_52) * r2_52) * _S921 + make_float2 (_S795 * v_52 + 0.0f * _S926 + (s_diff_r2_52 + (_S927 + _S927)) * _S790 + s_diff_r2_52 * _S791, _S796 * v_52 + 0.0f * _S928 + (s_diff_r2_52 + (_S922 + 0.0f * _S929)) * _S789 + s_diff_r2_52 * _S792);
                    float _S931 = 0.0f * u_52;
                    float s_diff_r2_53 = _S931 + _S931 + (v_52 + v_52);
                    float2  _S932 = _S797 * make_float2 (radial_17) + make_float2 (s_diff_r2_53 * _S925 + (s_diff_r2_53 * _S924 + (s_diff_r2_53 * _S923 + s_diff_r2_53 * _S788 * r2_52) * r2_52) * r2_52) * _S921 + make_float2 (_S798 * v_52 + _S926 + (s_diff_r2_53 + (_S931 + 0.0f * _S927)) * _S790 + s_diff_r2_53 * _S791, _S799 * v_52 + _S928 + (s_diff_r2_53 + (_S929 + _S929)) * _S789 + s_diff_r2_53 * _S792);
                    Matrix<float, 2, 2>  _S933 = transpose_0(makeMatrix<float, 2, 2> (_S930 + make_float2 (_S930.x * _S793 + _S930.y * _S794, 0.0f), _S932 + make_float2 (_S932.x * _S793 + _S932.y * _S794, 0.0f)));
                    bool _S934 = !((F32_min((determinant_0(_S933)), ((F32_min((_S933.rows[int(0)].x), (_S933.rows[int(1)].y)))))) > 0.0f);
                    _S813 = _S934;
                    if(_S934)
                    {
                        break;
                    }
                    float u_53 = proj_points_1[int(4)].x;
                    float v_53 = proj_points_1[int(4)].y;
                    float r2_53 = u_53 * u_53 + v_53 * v_53;
                    float2  _S935 = proj_points_1[int(4)] * make_float2 (1.0f + r2_53 * (_S785 + r2_53 * (_S786 + r2_53 * (_S787 + r2_53 * _S788)))) + make_float2 (_S795 * u_53 * v_53 + _S790 * (r2_53 + 2.0f * u_53 * u_53) + _S791 * r2_53, _S796 * u_53 * v_53 + _S789 * (r2_53 + 2.0f * v_53 * v_53) + _S792 * r2_53);
                    float2  _S936 = _S935 + make_float2 (_S793 * _S935.x + _S794 * _S935.y, 0.0f);
                    proj_points_1[int(4)] = make_float2 (fx_7 * _S936.x + cx_5, fy_7 * _S936.y + cy_5);
                    break;
                }
                bool all_valid_7 = all_valid_6 & (!_S813);
                for(;;)
                {
                    _S814 = &proj_points_1[int(5)];
                    for(;;)
                    {
                        float2  _S937 = float2 {_S835.p_0[int(5)].x, _S835.p_0[int(5)].y};
                        float r_16 = length_0(_S937);
                        float _S938 = _S835.p_0[int(5)].z;
                        _S815 = _S938;
                        float theta_9 = (F32_atan2((r_16), (_S938)));
                        if(theta_9 < 0.00100000004749745f)
                        {
                            k_4 = (1.0f - theta_9 * theta_9 / 3.0f) / _S938;
                        }
                        else
                        {
                            k_4 = theta_9 / r_16;
                        }
                        float2  _S939 = _S937 * make_float2 (k_4);
                        proj_points_1[int(5)] = _S939;
                        float u_54 = _S939.x;
                        float v_54 = _S939.y;
                        float _S940 = 0.0f * v_54;
                        float r2_54 = u_54 * u_54 + v_54 * v_54;
                        float s_diff_r2_54 = u_54 + u_54 + (_S940 + _S940);
                        float _S941 = _S787 + r2_54 * _S788;
                        float _S942 = _S786 + r2_54 * _S941;
                        float _S943 = _S785 + r2_54 * _S942;
                        float radial_18 = 1.0f + r2_54 * _S943;
                        float _S944 = _S795 * u_54;
                        float _S945 = 2.0f * u_54;
                        float _S946 = _S796 * u_54;
                        float _S947 = 2.0f * v_54;
                        float2  _S948 = _S784 * make_float2 (radial_18) + make_float2 (s_diff_r2_54 * _S943 + (s_diff_r2_54 * _S942 + (s_diff_r2_54 * _S941 + s_diff_r2_54 * _S788 * r2_54) * r2_54) * r2_54) * _S939 + make_float2 (_S795 * v_54 + 0.0f * _S944 + (s_diff_r2_54 + (_S945 + _S945)) * _S790 + s_diff_r2_54 * _S791, _S796 * v_54 + 0.0f * _S946 + (s_diff_r2_54 + (_S940 + 0.0f * _S947)) * _S789 + s_diff_r2_54 * _S792);
                        float _S949 = 0.0f * u_54;
                        float s_diff_r2_55 = _S949 + _S949 + (v_54 + v_54);
                        float2  _S950 = _S797 * make_float2 (radial_18) + make_float2 (s_diff_r2_55 * _S943 + (s_diff_r2_55 * _S942 + (s_diff_r2_55 * _S941 + s_diff_r2_55 * _S788 * r2_54) * r2_54) * r2_54) * _S939 + make_float2 (_S798 * v_54 + _S944 + (s_diff_r2_55 + (_S949 + 0.0f * _S945)) * _S790 + s_diff_r2_55 * _S791, _S799 * v_54 + _S946 + (s_diff_r2_55 + (_S947 + _S947)) * _S789 + s_diff_r2_55 * _S792);
                        Matrix<float, 2, 2>  _S951 = transpose_0(makeMatrix<float, 2, 2> (_S948 + make_float2 (_S948.x * _S793 + _S948.y * _S794, 0.0f), _S950 + make_float2 (_S950.x * _S793 + _S950.y * _S794, 0.0f)));
                        bool _S952 = !((F32_min((determinant_0(_S951)), ((F32_min((_S951.rows[int(0)].x), (_S951.rows[int(1)].y)))))) > 0.0f);
                        _S816 = _S952;
                        if(_S952)
                        {
                            break;
                        }
                        float u_55 = proj_points_1[int(5)].x;
                        float v_55 = proj_points_1[int(5)].y;
                        float r2_55 = u_55 * u_55 + v_55 * v_55;
                        float2  _S953 = proj_points_1[int(5)] * make_float2 (1.0f + r2_55 * (_S785 + r2_55 * (_S786 + r2_55 * (_S787 + r2_55 * _S788)))) + make_float2 (_S795 * u_55 * v_55 + _S790 * (r2_55 + 2.0f * u_55 * u_55) + _S791 * r2_55, _S796 * u_55 * v_55 + _S789 * (r2_55 + 2.0f * v_55 * v_55) + _S792 * r2_55);
                        float2  _S954 = _S953 + make_float2 (_S793 * _S953.x + _S794 * _S953.y, 0.0f);
                        proj_points_1[int(5)] = make_float2 (fx_7 * _S954.x + cx_5, fy_7 * _S954.y + cy_5);
                        break;
                    }
                    _S817 = all_valid_7 & (!_S816);
                    break;
                }
                _S818 = &proj_points_1[int(6)];
                for(;;)
                {
                    float2  _S955 = float2 {_S835.p_0[int(6)].x, _S835.p_0[int(6)].y};
                    float r_17 = length_0(_S955);
                    float _S956 = _S835.p_0[int(6)].z;
                    _S819 = _S956;
                    float theta_10 = (F32_atan2((r_17), (_S956)));
                    if(theta_10 < 0.00100000004749745f)
                    {
                        k_4 = (1.0f - theta_10 * theta_10 / 3.0f) / _S956;
                    }
                    else
                    {
                        k_4 = theta_10 / r_17;
                    }
                    float2  _S957 = _S955 * make_float2 (k_4);
                    proj_points_1[int(6)] = _S957;
                    float u_56 = _S957.x;
                    float v_56 = _S957.y;
                    float _S958 = 0.0f * v_56;
                    float r2_56 = u_56 * u_56 + v_56 * v_56;
                    float s_diff_r2_56 = u_56 + u_56 + (_S958 + _S958);
                    float _S959 = _S787 + r2_56 * _S788;
                    float _S960 = _S786 + r2_56 * _S959;
                    float _S961 = _S785 + r2_56 * _S960;
                    float radial_19 = 1.0f + r2_56 * _S961;
                    float _S962 = _S795 * u_56;
                    float _S963 = 2.0f * u_56;
                    float _S964 = _S796 * u_56;
                    float _S965 = 2.0f * v_56;
                    float2  _S966 = _S784 * make_float2 (radial_19) + make_float2 (s_diff_r2_56 * _S961 + (s_diff_r2_56 * _S960 + (s_diff_r2_56 * _S959 + s_diff_r2_56 * _S788 * r2_56) * r2_56) * r2_56) * _S957 + make_float2 (_S795 * v_56 + 0.0f * _S962 + (s_diff_r2_56 + (_S963 + _S963)) * _S790 + s_diff_r2_56 * _S791, _S796 * v_56 + 0.0f * _S964 + (s_diff_r2_56 + (_S958 + 0.0f * _S965)) * _S789 + s_diff_r2_56 * _S792);
                    float _S967 = 0.0f * u_56;
                    float s_diff_r2_57 = _S967 + _S967 + (v_56 + v_56);
                    float2  _S968 = _S797 * make_float2 (radial_19) + make_float2 (s_diff_r2_57 * _S961 + (s_diff_r2_57 * _S960 + (s_diff_r2_57 * _S959 + s_diff_r2_57 * _S788 * r2_56) * r2_56) * r2_56) * _S957 + make_float2 (_S798 * v_56 + _S962 + (s_diff_r2_57 + (_S967 + 0.0f * _S963)) * _S790 + s_diff_r2_57 * _S791, _S799 * v_56 + _S964 + (s_diff_r2_57 + (_S965 + _S965)) * _S789 + s_diff_r2_57 * _S792);
                    Matrix<float, 2, 2>  _S969 = transpose_0(makeMatrix<float, 2, 2> (_S966 + make_float2 (_S966.x * _S793 + _S966.y * _S794, 0.0f), _S968 + make_float2 (_S968.x * _S793 + _S968.y * _S794, 0.0f)));
                    bool _S970 = !((F32_min((determinant_0(_S969)), ((F32_min((_S969.rows[int(0)].x), (_S969.rows[int(1)].y)))))) > 0.0f);
                    _S820 = _S970;
                    if(_S970)
                    {
                        break;
                    }
                    float u_57 = proj_points_1[int(6)].x;
                    float v_57 = proj_points_1[int(6)].y;
                    float r2_57 = u_57 * u_57 + v_57 * v_57;
                    float2  _S971 = proj_points_1[int(6)] * make_float2 (1.0f + r2_57 * (_S785 + r2_57 * (_S786 + r2_57 * (_S787 + r2_57 * _S788)))) + make_float2 (_S795 * u_57 * v_57 + _S790 * (r2_57 + 2.0f * u_57 * u_57) + _S791 * r2_57, _S796 * u_57 * v_57 + _S789 * (r2_57 + 2.0f * v_57 * v_57) + _S792 * r2_57);
                    float2  _S972 = _S971 + make_float2 (_S793 * _S971.x + _S794 * _S971.y, 0.0f);
                    proj_points_1[int(6)] = make_float2 (fx_7 * _S972.x + cx_5, fy_7 * _S972.y + cy_5);
                    break;
                }
                _S821 = _S817 & (!_S820);
                break;
            }
            if(!_S821)
            {
                _S824 = false;
                break;
            }
            float2  p_1 = *_S782 + (*_S801 - *_S782) * make_float2 (3.32899999618530273f);
            float2  p_2 = *_S782 + (*_S804 - *_S782) * make_float2 (3.32899999618530273f);
            float2  p_3 = *_S782 + (*_S808 - *_S782) * make_float2 (3.32899999618530273f);
            float2  p_4 = *_S782 + (*_S811 - *_S782) * make_float2 (3.32899999618530273f);
            float2  p_5 = *_S782 + (*_S814 - *_S782) * make_float2 (3.32899999618530273f);
            float2  p_6 = *_S782 + (*_S818 - *_S782) * make_float2 (3.32899999618530273f);
            float2  _S973 = make_float2 (cx_5, cy_5);
            float2  min_p_0 = min_0(min_0(min_0(min_0(min_0(min_0(*_S782, p_1), p_2), p_3), p_4), p_5), p_6) - _S973;
            float2  max_p_0 = max_0(max_0(max_0(max_0(max_0(max_0(*_S782, p_1), p_2), p_3), p_4), p_5), p_6) - _S973;
            if((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S783), (_S802)))), (_S805)))), (_S809)))), (_S812)))), (_S815)))), (_S819))) <= 0.0f)
            {
                _S824 = (min_p_0.x * max_p_0.x) < 0.0f;
            }
            else
            {
                _S824 = false;
            }
            if(_S824)
            {
                _S824 = (min_p_0.y * max_p_0.y) < 0.0f;
            }
            else
            {
                _S824 = false;
            }
            if(_S824)
            {
                _S824 = false;
                break;
            }
            float2  _S974 = *mean2d_7 + make_float2 (_S835.w_mean_0[int(0)]) * *_S782 + make_float2 (_S835.w_mean_0[int(1)]) * *_S801 + make_float2 (_S835.w_mean_0[int(2)]) * *_S804 + make_float2 (_S835.w_mean_0[int(3)]) * *_S808 + make_float2 (_S835.w_mean_0[int(4)]) * *_S811 + make_float2 (_S835.w_mean_0[int(5)]) * *_S814 + make_float2 (_S835.w_mean_0[int(6)]) * *_S818;
            *mean2d_7 = _S974;
            float2  d_7 = *_S782 - _S974;
            float _S975 = d_7.x;
            float _S976 = d_7.y;
            float _S977 = _S975 * _S976;
            float2  d_8 = *_S801 - _S974;
            float _S978 = d_8.x;
            float _S979 = d_8.y;
            float _S980 = _S978 * _S979;
            float2  d_9 = *_S804 - _S974;
            float _S981 = d_9.x;
            float _S982 = d_9.y;
            float _S983 = _S981 * _S982;
            float2  d_10 = *_S808 - _S974;
            float _S984 = d_10.x;
            float _S985 = d_10.y;
            float _S986 = _S984 * _S985;
            float2  d_11 = *_S811 - _S974;
            float _S987 = d_11.x;
            float _S988 = d_11.y;
            float _S989 = _S987 * _S988;
            float2  d_12 = *_S814 - _S974;
            float _S990 = d_12.x;
            float _S991 = d_12.y;
            float _S992 = _S990 * _S991;
            float2  d_13 = *_S818 - _S974;
            float _S993 = d_13.x;
            float _S994 = d_13.y;
            float _S995 = _S993 * _S994;
            covar2d_4 = covar2d_4 + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S975 * _S975, _S977, _S977, _S976 * _S976) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S978 * _S978, _S980, _S980, _S979 * _S979) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S981 * _S981, _S983, _S983, _S982 * _S982) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S984 * _S984, _S986, _S986, _S985 * _S985) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S987 * _S987, _S989, _S989, _S988 * _S988) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S990 * _S990, _S992, _S992, _S991 * _S991) + makeMatrix<float, 2, 2> (_S835.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S993 * _S993, _S995, _S995, _S994 * _S994);
            _S824 = true;
            break;
        }
        if(!(true & _S824))
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
        float _S996 = *&(((&covar2d_4)->rows + (int(0)))->x) + eps2d_4;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S996;
        float _S997 = *&(((&covar2d_4)->rows + (int(1)))->y) + eps2d_4;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S997;
        float det_blur_4 = _S996 * _S997 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / det_blur_4))))));
        if(det_blur_4 <= 0.0f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        float invdet_6 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
        Matrix<float, 2, 2>  _S998 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_6, - covar2d_4.rows[int(0)].y * invdet_6, - covar2d_4.rows[int(1)].x * invdet_6, covar2d_4.rows[int(0)].x * invdet_6);
        if(antialiased_4)
        {
            *opacity_4 = *opacity_4 * compensation_4;
        }
        if((*opacity_4) < 0.00392156885936856f)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        float _S999 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
        float radius_x_4 = _S999 * (F32_sqrt((covar2d_4[int(0)].x)));
        float radius_y_4 = _S999 * (F32_sqrt((covar2d_4[int(1)].y)));
        float _S1000 = (*mean2d_7).x - radius_x_4;
        float _S1001 = (*mean2d_7).x + radius_x_4;
        float _S1002 = (*mean2d_7).y - radius_y_4;
        float _S1003 = (*mean2d_7).y + radius_y_4;
        if(_S1001 <= 0.0f)
        {
            _S824 = true;
        }
        else
        {
            _S824 = _S1000 >= float(image_width_4);
        }
        if(_S824)
        {
            _S824 = true;
        }
        else
        {
            _S824 = _S1003 <= 0.0f;
        }
        if(_S824)
        {
            _S824 = true;
        }
        else
        {
            _S824 = _S1002 >= float(image_height_4);
        }
        if(_S824)
        {
            *aabb_xyxy_4 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_4 = make_float4 (_S1000, _S1002, _S1001, _S1003);
        float x_22 = mean_c_3.x;
        float y_7 = mean_c_3.y;
        float _S1004 = x_22 * x_22 + y_7 * y_7;
        *sorting_depth_4 = _S823 * _S823 * _S823 * _S823 + 0.001953125f * _S1004 * _S1004;
        *conic_4 = make_float3 (_S998.rows[int(0)].x, _S998.rows[int(0)].y, _S998.rows[int(1)].y);
        *radius_5 = view_radius_3dgs_0(mean_5, scale_4, in_opacity_4, - mul_6(transpose_3(R_4), t_4));
        break;
    }
    return;
}

inline __device__ void projection_3dgut_equisolid(bool antialiased_5, float3  mean_6, float4  quat_5, float3  scale_5, float in_opacity_5, Matrix<float, 3, 3>  R_5, float3  t_5, float fx_8, float fy_8, float cx_6, float cy_6, FixedArray<float, 10>  dist_coeffs_8, uint image_width_5, uint image_height_5, float4  * aabb_xyxy_5, float * sorting_depth_5, float * radius_6, float2  * mean2d_8, float * depth_5, float3  * conic_5, float * opacity_5)
{
    float2  * _S1005;
    float _S1006;
    float2  _S1007;
    float _S1008;
    float _S1009;
    float _S1010;
    float _S1011;
    float _S1012;
    float _S1013;
    float _S1014;
    float _S1015;
    float _S1016;
    float _S1017;
    float _S1018;
    float _S1019;
    float2  _S1020;
    float _S1021;
    float _S1022;
    bool _S1023;
    float2  * _S1024;
    float _S1025;
    bool _S1026;
    float2  * _S1027;
    float _S1028;
    bool _S1029;
    bool _S1030;
    float2  * _S1031;
    float _S1032;
    bool _S1033;
    float2  * _S1034;
    float _S1035;
    bool _S1036;
    float2  * _S1037;
    float _S1038;
    bool _S1039;
    bool _S1040;
    float2  * _S1041;
    float _S1042;
    bool _S1043;
    bool _S1044;
    for(;;)
    {
        float3  mean_c_4 = mul_6(R_5, mean_6) + t_5;
        float _S1045 = length_1(mean_c_4);
        float _S1046 = mean_c_4.z;
        *depth_5 = _S1046;
        if(_S1045 <= 0.0f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        Matrix<float, 2, 2>  covar2d_5;
        *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
        bool _S1047;
        float3  _S1048 = exp_0(scale_5);
        float4  _S1049 = normalize_0(quat_5);
        float x_23 = _S1049.y;
        float x2_5 = x_23 * x_23;
        float y2_5 = _S1049.z * _S1049.z;
        float z2_5 = _S1049.w * _S1049.w;
        float xy_5 = _S1049.y * _S1049.z;
        float xz_5 = _S1049.y * _S1049.w;
        float yz_5 = _S1049.z * _S1049.w;
        float wx_5 = _S1049.x * _S1049.y;
        float wy_5 = _S1049.x * _S1049.z;
        float wz_5 = _S1049.x * _S1049.w;
        Matrix<float, 3, 3>  _S1050 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))));
        SigmaPoints_0 ret_2;
        (&ret_2)->p_0[int(0)] = mean_6;
        (&ret_2)->w_mean_0[int(0)] = 0.0f;
        (&ret_2)->w_cov_0[int(0)] = 2.0f;
        float _S1051 = (F32_sqrt((3.0f)));
        float3  delta_6 = make_float3 (_S1051 * _S1048.x) * _S1050.rows[0U];
        float3  _S1052 = mean_6 + delta_6;
        float3  _S1053 = mean_6 - delta_6;
        float3  delta_7 = make_float3 (_S1051 * _S1048.y) * _S1050.rows[1U];
        float3  _S1054 = mean_6 + delta_7;
        float3  _S1055 = mean_6 - delta_7;
        float3  delta_8 = make_float3 (_S1051 * _S1048.z) * _S1050.rows[2U];
        float3  _S1056 = mean_6 + delta_8;
        float3  _S1057 = mean_6 - delta_8;
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
        (&ret_2)->p_0[0U] = mul_6(R_5, (&ret_2)->p_0[0U]) + t_5;
        (&ret_2)->p_0[1U] = mul_6(R_5, _S1052) + t_5;
        (&ret_2)->p_0[2U] = mul_6(R_5, _S1054) + t_5;
        (&ret_2)->p_0[3U] = mul_6(R_5, _S1056) + t_5;
        (&ret_2)->p_0[4U] = mul_6(R_5, _S1053) + t_5;
        (&ret_2)->p_0[5U] = mul_6(R_5, _S1055) + t_5;
        (&ret_2)->p_0[6U] = mul_6(R_5, _S1057) + t_5;
        SigmaPoints_0 _S1058 = ret_2;
        for(;;)
        {
            int2  _S1059 = make_int2 (int(0));
            float2  _S1060 = make_float2 ((float)_S1059.x, (float)_S1059.y);
            *mean2d_8 = _S1060;
            covar2d_5 = makeMatrix<float, 2, 2> (0.0f);
            FixedArray<float2 , 7>  proj_points_2;
            for(;;)
            {
                float k_5;
                _S1005 = &proj_points_2[int(0)];
                for(;;)
                {
                    float2  _S1061 = float2 {_S1058.p_0[int(0)].x, _S1058.p_0[int(0)].y};
                    float r_18 = length_0(_S1061);
                    float _S1062 = _S1058.p_0[int(0)].z;
                    _S1006 = _S1062;
                    float theta_11 = (F32_atan2((r_18), (_S1062)));
                    if(r_18 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_11 * theta_11 / 24.0f) / _S1062;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_11))) / r_18;
                    }
                    float2  _S1063 = _S1061 * make_float2 (k_5);
                    proj_points_2[int(0)] = _S1063;
                    float2  _S1064 = make_float2 (1.0f, 0.0f);
                    _S1007 = _S1064;
                    _S1008 = dist_coeffs_8[int(0)];
                    _S1009 = dist_coeffs_8[int(1)];
                    _S1010 = dist_coeffs_8[int(2)];
                    _S1011 = dist_coeffs_8[int(3)];
                    _S1012 = dist_coeffs_8[int(4)];
                    _S1013 = dist_coeffs_8[int(5)];
                    _S1014 = dist_coeffs_8[int(6)];
                    _S1015 = dist_coeffs_8[int(7)];
                    _S1016 = dist_coeffs_8[int(8)];
                    _S1017 = dist_coeffs_8[int(9)];
                    float u_58 = _S1063.x;
                    float v_58 = _S1063.y;
                    float _S1065 = 0.0f * v_58;
                    float r2_58 = u_58 * u_58 + v_58 * v_58;
                    float s_diff_r2_58 = u_58 + u_58 + (_S1065 + _S1065);
                    float _S1066 = dist_coeffs_8[int(2)] + r2_58 * dist_coeffs_8[int(3)];
                    float _S1067 = dist_coeffs_8[int(1)] + r2_58 * _S1066;
                    float _S1068 = dist_coeffs_8[int(0)] + r2_58 * _S1067;
                    float _S1069 = s_diff_r2_58 * _S1068 + (s_diff_r2_58 * _S1067 + (s_diff_r2_58 * _S1066 + s_diff_r2_58 * dist_coeffs_8[int(3)] * r2_58) * r2_58) * r2_58;
                    float radial_20 = 1.0f + r2_58 * _S1068;
                    float _S1070 = 2.0f * dist_coeffs_8[int(4)];
                    _S1018 = _S1070;
                    float _S1071 = _S1070 * u_58;
                    float _S1072 = 2.0f * u_58;
                    float s_diff_du_6 = _S1070 * v_58 + 0.0f * _S1071 + (s_diff_r2_58 + (_S1072 + _S1072)) * dist_coeffs_8[int(5)] + s_diff_r2_58 * dist_coeffs_8[int(6)];
                    float _S1073 = 2.0f * dist_coeffs_8[int(5)];
                    _S1019 = _S1073;
                    float _S1074 = _S1073 * u_58;
                    float _S1075 = 2.0f * v_58;
                    float2  _S1076 = _S1064 * make_float2 (radial_20) + make_float2 (_S1069) * _S1063 + make_float2 (s_diff_du_6, _S1073 * v_58 + 0.0f * _S1074 + (s_diff_r2_58 + (_S1065 + 0.0f * _S1075)) * dist_coeffs_8[int(4)] + s_diff_r2_58 * dist_coeffs_8[int(7)]);
                    float2  _S1077 = _S1076 + make_float2 (_S1076.x * dist_coeffs_8[int(8)] + _S1076.y * dist_coeffs_8[int(9)], 0.0f);
                    float2  _S1078 = make_float2 (0.0f, 1.0f);
                    _S1020 = _S1078;
                    float _S1079 = 0.0f * u_58;
                    float s_diff_r2_59 = _S1079 + _S1079 + (v_58 + v_58);
                    float _S1080 = s_diff_r2_59 * _S1068 + (s_diff_r2_59 * _S1067 + (s_diff_r2_59 * _S1066 + s_diff_r2_59 * dist_coeffs_8[int(3)] * r2_58) * r2_58) * r2_58;
                    float _S1081 = 0.0f * _S1070;
                    _S1021 = _S1081;
                    float s_diff_du_7 = _S1081 * v_58 + _S1071 + (s_diff_r2_59 + (_S1079 + 0.0f * _S1072)) * dist_coeffs_8[int(5)] + s_diff_r2_59 * dist_coeffs_8[int(6)];
                    float _S1082 = 0.0f * _S1073;
                    _S1022 = _S1082;
                    float2  _S1083 = _S1078 * make_float2 (radial_20) + make_float2 (_S1080) * _S1063 + make_float2 (s_diff_du_7, _S1082 * v_58 + _S1074 + (s_diff_r2_59 + (_S1075 + _S1075)) * dist_coeffs_8[int(4)] + s_diff_r2_59 * dist_coeffs_8[int(7)]);
                    Matrix<float, 2, 2>  _S1084 = transpose_0(makeMatrix<float, 2, 2> (_S1077, _S1083 + make_float2 (_S1083.x * dist_coeffs_8[int(8)] + _S1083.y * dist_coeffs_8[int(9)], 0.0f)));
                    bool _S1085 = !((F32_min((determinant_0(_S1084)), ((F32_min((_S1084.rows[int(0)].x), (_S1084.rows[int(1)].y)))))) > 0.0f);
                    _S1023 = _S1085;
                    if(_S1085)
                    {
                        break;
                    }
                    float u_59 = proj_points_2[int(0)].x;
                    float v_59 = proj_points_2[int(0)].y;
                    float r2_59 = u_59 * u_59 + v_59 * v_59;
                    float2  _S1086 = proj_points_2[int(0)] * make_float2 (1.0f + r2_59 * (dist_coeffs_8[int(0)] + r2_59 * (dist_coeffs_8[int(1)] + r2_59 * (dist_coeffs_8[int(2)] + r2_59 * dist_coeffs_8[int(3)])))) + make_float2 (_S1070 * u_59 * v_59 + dist_coeffs_8[int(5)] * (r2_59 + 2.0f * u_59 * u_59) + dist_coeffs_8[int(6)] * r2_59, _S1073 * u_59 * v_59 + dist_coeffs_8[int(4)] * (r2_59 + 2.0f * v_59 * v_59) + dist_coeffs_8[int(7)] * r2_59);
                    float2  _S1087 = _S1086 + make_float2 (dist_coeffs_8[int(8)] * _S1086.x + dist_coeffs_8[int(9)] * _S1086.y, 0.0f);
                    proj_points_2[int(0)] = make_float2 (fx_8 * _S1087.x + cx_6, fy_8 * _S1087.y + cy_6);
                    break;
                }
                bool all_valid_8 = true & (!_S1023);
                _S1024 = &proj_points_2[int(1)];
                for(;;)
                {
                    float2  _S1088 = float2 {_S1058.p_0[int(1)].x, _S1058.p_0[int(1)].y};
                    float r_19 = length_0(_S1088);
                    float _S1089 = _S1058.p_0[int(1)].z;
                    _S1025 = _S1089;
                    float theta_12 = (F32_atan2((r_19), (_S1089)));
                    if(r_19 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_12 * theta_12 / 24.0f) / _S1089;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_12))) / r_19;
                    }
                    float2  _S1090 = _S1088 * make_float2 (k_5);
                    proj_points_2[int(1)] = _S1090;
                    float u_60 = _S1090.x;
                    float v_60 = _S1090.y;
                    float _S1091 = 0.0f * v_60;
                    float r2_60 = u_60 * u_60 + v_60 * v_60;
                    float s_diff_r2_60 = u_60 + u_60 + (_S1091 + _S1091);
                    float _S1092 = _S1010 + r2_60 * _S1011;
                    float _S1093 = _S1009 + r2_60 * _S1092;
                    float _S1094 = _S1008 + r2_60 * _S1093;
                    float radial_21 = 1.0f + r2_60 * _S1094;
                    float _S1095 = _S1018 * u_60;
                    float _S1096 = 2.0f * u_60;
                    float _S1097 = _S1019 * u_60;
                    float _S1098 = 2.0f * v_60;
                    float2  _S1099 = _S1007 * make_float2 (radial_21) + make_float2 (s_diff_r2_60 * _S1094 + (s_diff_r2_60 * _S1093 + (s_diff_r2_60 * _S1092 + s_diff_r2_60 * _S1011 * r2_60) * r2_60) * r2_60) * _S1090 + make_float2 (_S1018 * v_60 + 0.0f * _S1095 + (s_diff_r2_60 + (_S1096 + _S1096)) * _S1013 + s_diff_r2_60 * _S1014, _S1019 * v_60 + 0.0f * _S1097 + (s_diff_r2_60 + (_S1091 + 0.0f * _S1098)) * _S1012 + s_diff_r2_60 * _S1015);
                    float _S1100 = 0.0f * u_60;
                    float s_diff_r2_61 = _S1100 + _S1100 + (v_60 + v_60);
                    float2  _S1101 = _S1020 * make_float2 (radial_21) + make_float2 (s_diff_r2_61 * _S1094 + (s_diff_r2_61 * _S1093 + (s_diff_r2_61 * _S1092 + s_diff_r2_61 * _S1011 * r2_60) * r2_60) * r2_60) * _S1090 + make_float2 (_S1021 * v_60 + _S1095 + (s_diff_r2_61 + (_S1100 + 0.0f * _S1096)) * _S1013 + s_diff_r2_61 * _S1014, _S1022 * v_60 + _S1097 + (s_diff_r2_61 + (_S1098 + _S1098)) * _S1012 + s_diff_r2_61 * _S1015);
                    Matrix<float, 2, 2>  _S1102 = transpose_0(makeMatrix<float, 2, 2> (_S1099 + make_float2 (_S1099.x * _S1016 + _S1099.y * _S1017, 0.0f), _S1101 + make_float2 (_S1101.x * _S1016 + _S1101.y * _S1017, 0.0f)));
                    bool _S1103 = !((F32_min((determinant_0(_S1102)), ((F32_min((_S1102.rows[int(0)].x), (_S1102.rows[int(1)].y)))))) > 0.0f);
                    _S1026 = _S1103;
                    if(_S1103)
                    {
                        break;
                    }
                    float u_61 = proj_points_2[int(1)].x;
                    float v_61 = proj_points_2[int(1)].y;
                    float r2_61 = u_61 * u_61 + v_61 * v_61;
                    float2  _S1104 = proj_points_2[int(1)] * make_float2 (1.0f + r2_61 * (_S1008 + r2_61 * (_S1009 + r2_61 * (_S1010 + r2_61 * _S1011)))) + make_float2 (_S1018 * u_61 * v_61 + _S1013 * (r2_61 + 2.0f * u_61 * u_61) + _S1014 * r2_61, _S1019 * u_61 * v_61 + _S1012 * (r2_61 + 2.0f * v_61 * v_61) + _S1015 * r2_61);
                    float2  _S1105 = _S1104 + make_float2 (_S1016 * _S1104.x + _S1017 * _S1104.y, 0.0f);
                    proj_points_2[int(1)] = make_float2 (fx_8 * _S1105.x + cx_6, fy_8 * _S1105.y + cy_6);
                    break;
                }
                bool all_valid_9 = all_valid_8 & (!_S1026);
                for(;;)
                {
                    _S1027 = &proj_points_2[int(2)];
                    for(;;)
                    {
                        float2  _S1106 = float2 {_S1058.p_0[int(2)].x, _S1058.p_0[int(2)].y};
                        float r_20 = length_0(_S1106);
                        float _S1107 = _S1058.p_0[int(2)].z;
                        _S1028 = _S1107;
                        float theta_13 = (F32_atan2((r_20), (_S1107)));
                        if(r_20 < 9.99999997475242708e-07f)
                        {
                            k_5 = (1.0f - theta_13 * theta_13 / 24.0f) / _S1107;
                        }
                        else
                        {
                            k_5 = 2.0f * (F32_sin((0.5f * theta_13))) / r_20;
                        }
                        float2  _S1108 = _S1106 * make_float2 (k_5);
                        proj_points_2[int(2)] = _S1108;
                        float u_62 = _S1108.x;
                        float v_62 = _S1108.y;
                        float _S1109 = 0.0f * v_62;
                        float r2_62 = u_62 * u_62 + v_62 * v_62;
                        float s_diff_r2_62 = u_62 + u_62 + (_S1109 + _S1109);
                        float _S1110 = _S1010 + r2_62 * _S1011;
                        float _S1111 = _S1009 + r2_62 * _S1110;
                        float _S1112 = _S1008 + r2_62 * _S1111;
                        float radial_22 = 1.0f + r2_62 * _S1112;
                        float _S1113 = _S1018 * u_62;
                        float _S1114 = 2.0f * u_62;
                        float _S1115 = _S1019 * u_62;
                        float _S1116 = 2.0f * v_62;
                        float2  _S1117 = _S1007 * make_float2 (radial_22) + make_float2 (s_diff_r2_62 * _S1112 + (s_diff_r2_62 * _S1111 + (s_diff_r2_62 * _S1110 + s_diff_r2_62 * _S1011 * r2_62) * r2_62) * r2_62) * _S1108 + make_float2 (_S1018 * v_62 + 0.0f * _S1113 + (s_diff_r2_62 + (_S1114 + _S1114)) * _S1013 + s_diff_r2_62 * _S1014, _S1019 * v_62 + 0.0f * _S1115 + (s_diff_r2_62 + (_S1109 + 0.0f * _S1116)) * _S1012 + s_diff_r2_62 * _S1015);
                        float _S1118 = 0.0f * u_62;
                        float s_diff_r2_63 = _S1118 + _S1118 + (v_62 + v_62);
                        float2  _S1119 = _S1020 * make_float2 (radial_22) + make_float2 (s_diff_r2_63 * _S1112 + (s_diff_r2_63 * _S1111 + (s_diff_r2_63 * _S1110 + s_diff_r2_63 * _S1011 * r2_62) * r2_62) * r2_62) * _S1108 + make_float2 (_S1021 * v_62 + _S1113 + (s_diff_r2_63 + (_S1118 + 0.0f * _S1114)) * _S1013 + s_diff_r2_63 * _S1014, _S1022 * v_62 + _S1115 + (s_diff_r2_63 + (_S1116 + _S1116)) * _S1012 + s_diff_r2_63 * _S1015);
                        Matrix<float, 2, 2>  _S1120 = transpose_0(makeMatrix<float, 2, 2> (_S1117 + make_float2 (_S1117.x * _S1016 + _S1117.y * _S1017, 0.0f), _S1119 + make_float2 (_S1119.x * _S1016 + _S1119.y * _S1017, 0.0f)));
                        bool _S1121 = !((F32_min((determinant_0(_S1120)), ((F32_min((_S1120.rows[int(0)].x), (_S1120.rows[int(1)].y)))))) > 0.0f);
                        _S1029 = _S1121;
                        if(_S1121)
                        {
                            break;
                        }
                        float u_63 = proj_points_2[int(2)].x;
                        float v_63 = proj_points_2[int(2)].y;
                        float r2_63 = u_63 * u_63 + v_63 * v_63;
                        float2  _S1122 = proj_points_2[int(2)] * make_float2 (1.0f + r2_63 * (_S1008 + r2_63 * (_S1009 + r2_63 * (_S1010 + r2_63 * _S1011)))) + make_float2 (_S1018 * u_63 * v_63 + _S1013 * (r2_63 + 2.0f * u_63 * u_63) + _S1014 * r2_63, _S1019 * u_63 * v_63 + _S1012 * (r2_63 + 2.0f * v_63 * v_63) + _S1015 * r2_63);
                        float2  _S1123 = _S1122 + make_float2 (_S1016 * _S1122.x + _S1017 * _S1122.y, 0.0f);
                        proj_points_2[int(2)] = make_float2 (fx_8 * _S1123.x + cx_6, fy_8 * _S1123.y + cy_6);
                        break;
                    }
                    _S1030 = all_valid_9 & (!_S1029);
                    break;
                }
                _S1031 = &proj_points_2[int(3)];
                for(;;)
                {
                    float2  _S1124 = float2 {_S1058.p_0[int(3)].x, _S1058.p_0[int(3)].y};
                    float r_21 = length_0(_S1124);
                    float _S1125 = _S1058.p_0[int(3)].z;
                    _S1032 = _S1125;
                    float theta_14 = (F32_atan2((r_21), (_S1125)));
                    if(r_21 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_14 * theta_14 / 24.0f) / _S1125;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_14))) / r_21;
                    }
                    float2  _S1126 = _S1124 * make_float2 (k_5);
                    proj_points_2[int(3)] = _S1126;
                    float u_64 = _S1126.x;
                    float v_64 = _S1126.y;
                    float _S1127 = 0.0f * v_64;
                    float r2_64 = u_64 * u_64 + v_64 * v_64;
                    float s_diff_r2_64 = u_64 + u_64 + (_S1127 + _S1127);
                    float _S1128 = _S1010 + r2_64 * _S1011;
                    float _S1129 = _S1009 + r2_64 * _S1128;
                    float _S1130 = _S1008 + r2_64 * _S1129;
                    float radial_23 = 1.0f + r2_64 * _S1130;
                    float _S1131 = _S1018 * u_64;
                    float _S1132 = 2.0f * u_64;
                    float _S1133 = _S1019 * u_64;
                    float _S1134 = 2.0f * v_64;
                    float2  _S1135 = _S1007 * make_float2 (radial_23) + make_float2 (s_diff_r2_64 * _S1130 + (s_diff_r2_64 * _S1129 + (s_diff_r2_64 * _S1128 + s_diff_r2_64 * _S1011 * r2_64) * r2_64) * r2_64) * _S1126 + make_float2 (_S1018 * v_64 + 0.0f * _S1131 + (s_diff_r2_64 + (_S1132 + _S1132)) * _S1013 + s_diff_r2_64 * _S1014, _S1019 * v_64 + 0.0f * _S1133 + (s_diff_r2_64 + (_S1127 + 0.0f * _S1134)) * _S1012 + s_diff_r2_64 * _S1015);
                    float _S1136 = 0.0f * u_64;
                    float s_diff_r2_65 = _S1136 + _S1136 + (v_64 + v_64);
                    float2  _S1137 = _S1020 * make_float2 (radial_23) + make_float2 (s_diff_r2_65 * _S1130 + (s_diff_r2_65 * _S1129 + (s_diff_r2_65 * _S1128 + s_diff_r2_65 * _S1011 * r2_64) * r2_64) * r2_64) * _S1126 + make_float2 (_S1021 * v_64 + _S1131 + (s_diff_r2_65 + (_S1136 + 0.0f * _S1132)) * _S1013 + s_diff_r2_65 * _S1014, _S1022 * v_64 + _S1133 + (s_diff_r2_65 + (_S1134 + _S1134)) * _S1012 + s_diff_r2_65 * _S1015);
                    Matrix<float, 2, 2>  _S1138 = transpose_0(makeMatrix<float, 2, 2> (_S1135 + make_float2 (_S1135.x * _S1016 + _S1135.y * _S1017, 0.0f), _S1137 + make_float2 (_S1137.x * _S1016 + _S1137.y * _S1017, 0.0f)));
                    bool _S1139 = !((F32_min((determinant_0(_S1138)), ((F32_min((_S1138.rows[int(0)].x), (_S1138.rows[int(1)].y)))))) > 0.0f);
                    _S1033 = _S1139;
                    if(_S1139)
                    {
                        break;
                    }
                    float u_65 = proj_points_2[int(3)].x;
                    float v_65 = proj_points_2[int(3)].y;
                    float r2_65 = u_65 * u_65 + v_65 * v_65;
                    float2  _S1140 = proj_points_2[int(3)] * make_float2 (1.0f + r2_65 * (_S1008 + r2_65 * (_S1009 + r2_65 * (_S1010 + r2_65 * _S1011)))) + make_float2 (_S1018 * u_65 * v_65 + _S1013 * (r2_65 + 2.0f * u_65 * u_65) + _S1014 * r2_65, _S1019 * u_65 * v_65 + _S1012 * (r2_65 + 2.0f * v_65 * v_65) + _S1015 * r2_65);
                    float2  _S1141 = _S1140 + make_float2 (_S1016 * _S1140.x + _S1017 * _S1140.y, 0.0f);
                    proj_points_2[int(3)] = make_float2 (fx_8 * _S1141.x + cx_6, fy_8 * _S1141.y + cy_6);
                    break;
                }
                bool all_valid_10 = _S1030 & (!_S1033);
                _S1034 = &proj_points_2[int(4)];
                for(;;)
                {
                    float2  _S1142 = float2 {_S1058.p_0[int(4)].x, _S1058.p_0[int(4)].y};
                    float r_22 = length_0(_S1142);
                    float _S1143 = _S1058.p_0[int(4)].z;
                    _S1035 = _S1143;
                    float theta_15 = (F32_atan2((r_22), (_S1143)));
                    if(r_22 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_15 * theta_15 / 24.0f) / _S1143;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_15))) / r_22;
                    }
                    float2  _S1144 = _S1142 * make_float2 (k_5);
                    proj_points_2[int(4)] = _S1144;
                    float u_66 = _S1144.x;
                    float v_66 = _S1144.y;
                    float _S1145 = 0.0f * v_66;
                    float r2_66 = u_66 * u_66 + v_66 * v_66;
                    float s_diff_r2_66 = u_66 + u_66 + (_S1145 + _S1145);
                    float _S1146 = _S1010 + r2_66 * _S1011;
                    float _S1147 = _S1009 + r2_66 * _S1146;
                    float _S1148 = _S1008 + r2_66 * _S1147;
                    float radial_24 = 1.0f + r2_66 * _S1148;
                    float _S1149 = _S1018 * u_66;
                    float _S1150 = 2.0f * u_66;
                    float _S1151 = _S1019 * u_66;
                    float _S1152 = 2.0f * v_66;
                    float2  _S1153 = _S1007 * make_float2 (radial_24) + make_float2 (s_diff_r2_66 * _S1148 + (s_diff_r2_66 * _S1147 + (s_diff_r2_66 * _S1146 + s_diff_r2_66 * _S1011 * r2_66) * r2_66) * r2_66) * _S1144 + make_float2 (_S1018 * v_66 + 0.0f * _S1149 + (s_diff_r2_66 + (_S1150 + _S1150)) * _S1013 + s_diff_r2_66 * _S1014, _S1019 * v_66 + 0.0f * _S1151 + (s_diff_r2_66 + (_S1145 + 0.0f * _S1152)) * _S1012 + s_diff_r2_66 * _S1015);
                    float _S1154 = 0.0f * u_66;
                    float s_diff_r2_67 = _S1154 + _S1154 + (v_66 + v_66);
                    float2  _S1155 = _S1020 * make_float2 (radial_24) + make_float2 (s_diff_r2_67 * _S1148 + (s_diff_r2_67 * _S1147 + (s_diff_r2_67 * _S1146 + s_diff_r2_67 * _S1011 * r2_66) * r2_66) * r2_66) * _S1144 + make_float2 (_S1021 * v_66 + _S1149 + (s_diff_r2_67 + (_S1154 + 0.0f * _S1150)) * _S1013 + s_diff_r2_67 * _S1014, _S1022 * v_66 + _S1151 + (s_diff_r2_67 + (_S1152 + _S1152)) * _S1012 + s_diff_r2_67 * _S1015);
                    Matrix<float, 2, 2>  _S1156 = transpose_0(makeMatrix<float, 2, 2> (_S1153 + make_float2 (_S1153.x * _S1016 + _S1153.y * _S1017, 0.0f), _S1155 + make_float2 (_S1155.x * _S1016 + _S1155.y * _S1017, 0.0f)));
                    bool _S1157 = !((F32_min((determinant_0(_S1156)), ((F32_min((_S1156.rows[int(0)].x), (_S1156.rows[int(1)].y)))))) > 0.0f);
                    _S1036 = _S1157;
                    if(_S1157)
                    {
                        break;
                    }
                    float u_67 = proj_points_2[int(4)].x;
                    float v_67 = proj_points_2[int(4)].y;
                    float r2_67 = u_67 * u_67 + v_67 * v_67;
                    float2  _S1158 = proj_points_2[int(4)] * make_float2 (1.0f + r2_67 * (_S1008 + r2_67 * (_S1009 + r2_67 * (_S1010 + r2_67 * _S1011)))) + make_float2 (_S1018 * u_67 * v_67 + _S1013 * (r2_67 + 2.0f * u_67 * u_67) + _S1014 * r2_67, _S1019 * u_67 * v_67 + _S1012 * (r2_67 + 2.0f * v_67 * v_67) + _S1015 * r2_67);
                    float2  _S1159 = _S1158 + make_float2 (_S1016 * _S1158.x + _S1017 * _S1158.y, 0.0f);
                    proj_points_2[int(4)] = make_float2 (fx_8 * _S1159.x + cx_6, fy_8 * _S1159.y + cy_6);
                    break;
                }
                bool all_valid_11 = all_valid_10 & (!_S1036);
                for(;;)
                {
                    _S1037 = &proj_points_2[int(5)];
                    for(;;)
                    {
                        float2  _S1160 = float2 {_S1058.p_0[int(5)].x, _S1058.p_0[int(5)].y};
                        float r_23 = length_0(_S1160);
                        float _S1161 = _S1058.p_0[int(5)].z;
                        _S1038 = _S1161;
                        float theta_16 = (F32_atan2((r_23), (_S1161)));
                        if(r_23 < 9.99999997475242708e-07f)
                        {
                            k_5 = (1.0f - theta_16 * theta_16 / 24.0f) / _S1161;
                        }
                        else
                        {
                            k_5 = 2.0f * (F32_sin((0.5f * theta_16))) / r_23;
                        }
                        float2  _S1162 = _S1160 * make_float2 (k_5);
                        proj_points_2[int(5)] = _S1162;
                        float u_68 = _S1162.x;
                        float v_68 = _S1162.y;
                        float _S1163 = 0.0f * v_68;
                        float r2_68 = u_68 * u_68 + v_68 * v_68;
                        float s_diff_r2_68 = u_68 + u_68 + (_S1163 + _S1163);
                        float _S1164 = _S1010 + r2_68 * _S1011;
                        float _S1165 = _S1009 + r2_68 * _S1164;
                        float _S1166 = _S1008 + r2_68 * _S1165;
                        float radial_25 = 1.0f + r2_68 * _S1166;
                        float _S1167 = _S1018 * u_68;
                        float _S1168 = 2.0f * u_68;
                        float _S1169 = _S1019 * u_68;
                        float _S1170 = 2.0f * v_68;
                        float2  _S1171 = _S1007 * make_float2 (radial_25) + make_float2 (s_diff_r2_68 * _S1166 + (s_diff_r2_68 * _S1165 + (s_diff_r2_68 * _S1164 + s_diff_r2_68 * _S1011 * r2_68) * r2_68) * r2_68) * _S1162 + make_float2 (_S1018 * v_68 + 0.0f * _S1167 + (s_diff_r2_68 + (_S1168 + _S1168)) * _S1013 + s_diff_r2_68 * _S1014, _S1019 * v_68 + 0.0f * _S1169 + (s_diff_r2_68 + (_S1163 + 0.0f * _S1170)) * _S1012 + s_diff_r2_68 * _S1015);
                        float _S1172 = 0.0f * u_68;
                        float s_diff_r2_69 = _S1172 + _S1172 + (v_68 + v_68);
                        float2  _S1173 = _S1020 * make_float2 (radial_25) + make_float2 (s_diff_r2_69 * _S1166 + (s_diff_r2_69 * _S1165 + (s_diff_r2_69 * _S1164 + s_diff_r2_69 * _S1011 * r2_68) * r2_68) * r2_68) * _S1162 + make_float2 (_S1021 * v_68 + _S1167 + (s_diff_r2_69 + (_S1172 + 0.0f * _S1168)) * _S1013 + s_diff_r2_69 * _S1014, _S1022 * v_68 + _S1169 + (s_diff_r2_69 + (_S1170 + _S1170)) * _S1012 + s_diff_r2_69 * _S1015);
                        Matrix<float, 2, 2>  _S1174 = transpose_0(makeMatrix<float, 2, 2> (_S1171 + make_float2 (_S1171.x * _S1016 + _S1171.y * _S1017, 0.0f), _S1173 + make_float2 (_S1173.x * _S1016 + _S1173.y * _S1017, 0.0f)));
                        bool _S1175 = !((F32_min((determinant_0(_S1174)), ((F32_min((_S1174.rows[int(0)].x), (_S1174.rows[int(1)].y)))))) > 0.0f);
                        _S1039 = _S1175;
                        if(_S1175)
                        {
                            break;
                        }
                        float u_69 = proj_points_2[int(5)].x;
                        float v_69 = proj_points_2[int(5)].y;
                        float r2_69 = u_69 * u_69 + v_69 * v_69;
                        float2  _S1176 = proj_points_2[int(5)] * make_float2 (1.0f + r2_69 * (_S1008 + r2_69 * (_S1009 + r2_69 * (_S1010 + r2_69 * _S1011)))) + make_float2 (_S1018 * u_69 * v_69 + _S1013 * (r2_69 + 2.0f * u_69 * u_69) + _S1014 * r2_69, _S1019 * u_69 * v_69 + _S1012 * (r2_69 + 2.0f * v_69 * v_69) + _S1015 * r2_69);
                        float2  _S1177 = _S1176 + make_float2 (_S1016 * _S1176.x + _S1017 * _S1176.y, 0.0f);
                        proj_points_2[int(5)] = make_float2 (fx_8 * _S1177.x + cx_6, fy_8 * _S1177.y + cy_6);
                        break;
                    }
                    _S1040 = all_valid_11 & (!_S1039);
                    break;
                }
                _S1041 = &proj_points_2[int(6)];
                for(;;)
                {
                    float2  _S1178 = float2 {_S1058.p_0[int(6)].x, _S1058.p_0[int(6)].y};
                    float r_24 = length_0(_S1178);
                    float _S1179 = _S1058.p_0[int(6)].z;
                    _S1042 = _S1179;
                    float theta_17 = (F32_atan2((r_24), (_S1179)));
                    if(r_24 < 9.99999997475242708e-07f)
                    {
                        k_5 = (1.0f - theta_17 * theta_17 / 24.0f) / _S1179;
                    }
                    else
                    {
                        k_5 = 2.0f * (F32_sin((0.5f * theta_17))) / r_24;
                    }
                    float2  _S1180 = _S1178 * make_float2 (k_5);
                    proj_points_2[int(6)] = _S1180;
                    float u_70 = _S1180.x;
                    float v_70 = _S1180.y;
                    float _S1181 = 0.0f * v_70;
                    float r2_70 = u_70 * u_70 + v_70 * v_70;
                    float s_diff_r2_70 = u_70 + u_70 + (_S1181 + _S1181);
                    float _S1182 = _S1010 + r2_70 * _S1011;
                    float _S1183 = _S1009 + r2_70 * _S1182;
                    float _S1184 = _S1008 + r2_70 * _S1183;
                    float radial_26 = 1.0f + r2_70 * _S1184;
                    float _S1185 = _S1018 * u_70;
                    float _S1186 = 2.0f * u_70;
                    float _S1187 = _S1019 * u_70;
                    float _S1188 = 2.0f * v_70;
                    float2  _S1189 = _S1007 * make_float2 (radial_26) + make_float2 (s_diff_r2_70 * _S1184 + (s_diff_r2_70 * _S1183 + (s_diff_r2_70 * _S1182 + s_diff_r2_70 * _S1011 * r2_70) * r2_70) * r2_70) * _S1180 + make_float2 (_S1018 * v_70 + 0.0f * _S1185 + (s_diff_r2_70 + (_S1186 + _S1186)) * _S1013 + s_diff_r2_70 * _S1014, _S1019 * v_70 + 0.0f * _S1187 + (s_diff_r2_70 + (_S1181 + 0.0f * _S1188)) * _S1012 + s_diff_r2_70 * _S1015);
                    float _S1190 = 0.0f * u_70;
                    float s_diff_r2_71 = _S1190 + _S1190 + (v_70 + v_70);
                    float2  _S1191 = _S1020 * make_float2 (radial_26) + make_float2 (s_diff_r2_71 * _S1184 + (s_diff_r2_71 * _S1183 + (s_diff_r2_71 * _S1182 + s_diff_r2_71 * _S1011 * r2_70) * r2_70) * r2_70) * _S1180 + make_float2 (_S1021 * v_70 + _S1185 + (s_diff_r2_71 + (_S1190 + 0.0f * _S1186)) * _S1013 + s_diff_r2_71 * _S1014, _S1022 * v_70 + _S1187 + (s_diff_r2_71 + (_S1188 + _S1188)) * _S1012 + s_diff_r2_71 * _S1015);
                    Matrix<float, 2, 2>  _S1192 = transpose_0(makeMatrix<float, 2, 2> (_S1189 + make_float2 (_S1189.x * _S1016 + _S1189.y * _S1017, 0.0f), _S1191 + make_float2 (_S1191.x * _S1016 + _S1191.y * _S1017, 0.0f)));
                    bool _S1193 = !((F32_min((determinant_0(_S1192)), ((F32_min((_S1192.rows[int(0)].x), (_S1192.rows[int(1)].y)))))) > 0.0f);
                    _S1043 = _S1193;
                    if(_S1193)
                    {
                        break;
                    }
                    float u_71 = proj_points_2[int(6)].x;
                    float v_71 = proj_points_2[int(6)].y;
                    float r2_71 = u_71 * u_71 + v_71 * v_71;
                    float2  _S1194 = proj_points_2[int(6)] * make_float2 (1.0f + r2_71 * (_S1008 + r2_71 * (_S1009 + r2_71 * (_S1010 + r2_71 * _S1011)))) + make_float2 (_S1018 * u_71 * v_71 + _S1013 * (r2_71 + 2.0f * u_71 * u_71) + _S1014 * r2_71, _S1019 * u_71 * v_71 + _S1012 * (r2_71 + 2.0f * v_71 * v_71) + _S1015 * r2_71);
                    float2  _S1195 = _S1194 + make_float2 (_S1016 * _S1194.x + _S1017 * _S1194.y, 0.0f);
                    proj_points_2[int(6)] = make_float2 (fx_8 * _S1195.x + cx_6, fy_8 * _S1195.y + cy_6);
                    break;
                }
                _S1044 = _S1040 & (!_S1043);
                break;
            }
            if(!_S1044)
            {
                _S1047 = false;
                break;
            }
            float2  p_7 = *_S1005 + (*_S1024 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  p_8 = *_S1005 + (*_S1027 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  p_9 = *_S1005 + (*_S1031 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  p_10 = *_S1005 + (*_S1034 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  p_11 = *_S1005 + (*_S1037 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  p_12 = *_S1005 + (*_S1041 - *_S1005) * make_float2 (3.32899999618530273f);
            float2  _S1196 = make_float2 (cx_6, cy_6);
            float2  min_p_1 = min_0(min_0(min_0(min_0(min_0(min_0(*_S1005, p_7), p_8), p_9), p_10), p_11), p_12) - _S1196;
            float2  max_p_1 = max_0(max_0(max_0(max_0(max_0(max_0(*_S1005, p_7), p_8), p_9), p_10), p_11), p_12) - _S1196;
            if((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S1006), (_S1025)))), (_S1028)))), (_S1032)))), (_S1035)))), (_S1038)))), (_S1042))) <= 0.0f)
            {
                _S1047 = (min_p_1.x * max_p_1.x) < 0.0f;
            }
            else
            {
                _S1047 = false;
            }
            if(_S1047)
            {
                _S1047 = (min_p_1.y * max_p_1.y) < 0.0f;
            }
            else
            {
                _S1047 = false;
            }
            if(_S1047)
            {
                _S1047 = false;
                break;
            }
            float2  _S1197 = *mean2d_8 + make_float2 (_S1058.w_mean_0[int(0)]) * *_S1005 + make_float2 (_S1058.w_mean_0[int(1)]) * *_S1024 + make_float2 (_S1058.w_mean_0[int(2)]) * *_S1027 + make_float2 (_S1058.w_mean_0[int(3)]) * *_S1031 + make_float2 (_S1058.w_mean_0[int(4)]) * *_S1034 + make_float2 (_S1058.w_mean_0[int(5)]) * *_S1037 + make_float2 (_S1058.w_mean_0[int(6)]) * *_S1041;
            *mean2d_8 = _S1197;
            float2  d_14 = *_S1005 - _S1197;
            float _S1198 = d_14.x;
            float _S1199 = d_14.y;
            float _S1200 = _S1198 * _S1199;
            float2  d_15 = *_S1024 - _S1197;
            float _S1201 = d_15.x;
            float _S1202 = d_15.y;
            float _S1203 = _S1201 * _S1202;
            float2  d_16 = *_S1027 - _S1197;
            float _S1204 = d_16.x;
            float _S1205 = d_16.y;
            float _S1206 = _S1204 * _S1205;
            float2  d_17 = *_S1031 - _S1197;
            float _S1207 = d_17.x;
            float _S1208 = d_17.y;
            float _S1209 = _S1207 * _S1208;
            float2  d_18 = *_S1034 - _S1197;
            float _S1210 = d_18.x;
            float _S1211 = d_18.y;
            float _S1212 = _S1210 * _S1211;
            float2  d_19 = *_S1037 - _S1197;
            float _S1213 = d_19.x;
            float _S1214 = d_19.y;
            float _S1215 = _S1213 * _S1214;
            float2  d_20 = *_S1041 - _S1197;
            float _S1216 = d_20.x;
            float _S1217 = d_20.y;
            float _S1218 = _S1216 * _S1217;
            covar2d_5 = covar2d_5 + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S1198 * _S1198, _S1200, _S1200, _S1199 * _S1199) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S1201 * _S1201, _S1203, _S1203, _S1202 * _S1202) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S1204 * _S1204, _S1206, _S1206, _S1205 * _S1205) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S1207 * _S1207, _S1209, _S1209, _S1208 * _S1208) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S1210 * _S1210, _S1212, _S1212, _S1211 * _S1211) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S1213 * _S1213, _S1215, _S1215, _S1214 * _S1214) + makeMatrix<float, 2, 2> (_S1058.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S1216 * _S1216, _S1218, _S1218, _S1217 * _S1217);
            _S1047 = true;
            break;
        }
        if(!(true & _S1047))
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
        float _S1219 = *&(((&covar2d_5)->rows + (int(0)))->x) + eps2d_5;
        *&(((&covar2d_5)->rows + (int(0)))->x) = _S1219;
        float _S1220 = *&(((&covar2d_5)->rows + (int(1)))->y) + eps2d_5;
        *&(((&covar2d_5)->rows + (int(1)))->y) = _S1220;
        float det_blur_5 = _S1219 * _S1220 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
        float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / det_blur_5))))));
        if(det_blur_5 <= 0.0f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        float invdet_7 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
        Matrix<float, 2, 2>  _S1221 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_7, - covar2d_5.rows[int(0)].y * invdet_7, - covar2d_5.rows[int(1)].x * invdet_7, covar2d_5.rows[int(0)].x * invdet_7);
        if(antialiased_5)
        {
            *opacity_5 = *opacity_5 * compensation_5;
        }
        if((*opacity_5) < 0.00392156885936856f)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        float _S1222 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
        float radius_x_5 = _S1222 * (F32_sqrt((covar2d_5[int(0)].x)));
        float radius_y_5 = _S1222 * (F32_sqrt((covar2d_5[int(1)].y)));
        float _S1223 = (*mean2d_8).x - radius_x_5;
        float _S1224 = (*mean2d_8).x + radius_x_5;
        float _S1225 = (*mean2d_8).y - radius_y_5;
        float _S1226 = (*mean2d_8).y + radius_y_5;
        if(_S1224 <= 0.0f)
        {
            _S1047 = true;
        }
        else
        {
            _S1047 = _S1223 >= float(image_width_5);
        }
        if(_S1047)
        {
            _S1047 = true;
        }
        else
        {
            _S1047 = _S1226 <= 0.0f;
        }
        if(_S1047)
        {
            _S1047 = true;
        }
        else
        {
            _S1047 = _S1225 >= float(image_height_5);
        }
        if(_S1047)
        {
            *aabb_xyxy_5 = make_float4 (0.0f);
            break;
        }
        *aabb_xyxy_5 = make_float4 (_S1223, _S1225, _S1224, _S1226);
        float x_24 = mean_c_4.x;
        float y_8 = mean_c_4.y;
        float _S1227 = x_24 * x_24 + y_8 * y_8;
        *sorting_depth_5 = _S1046 * _S1046 * _S1046 * _S1046 + 0.001953125f * _S1227 * _S1227;
        *conic_5 = make_float3 (_S1221.rows[int(0)].x, _S1221.rows[int(0)].y, _S1221.rows[int(1)].y);
        *radius_6 = view_radius_3dgs_0(mean_6, scale_5, in_opacity_5, - mul_6(transpose_3(R_5), t_5));
        break;
    }
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S1228, float3  _S1229)
{
    return mul_6(_S1228, _S1229);
}

inline __device__ float s_primal_ctx_exp_0(float _S1230)
{
    return (F32_exp((_S1230)));
}

inline __device__ float3  s_primal_ctx_exp_1(float3  _S1231)
{
    return exp_0(_S1231);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S1232, Matrix<float, 3, 3>  _S1233)
{
    return mul_5(_S1232, _S1233);
}

inline __device__ float s_primal_ctx_clamp_0(float _S1234, float _S1235, float _S1236)
{
    return clamp_0(_S1234, _S1235, _S1236);
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S1237, Matrix<float, 3, 3>  _S1238)
{
    return mul_3(_S1237, _S1238);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S1239, Matrix<float, 3, 2>  _S1240)
{
    return mul_4(_S1239, _S1240);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S1241)
{
    return (F32_sqrt((_S1241)));
}

inline __device__ float s_primal_ctx_log_0(float _S1242)
{
    return (F32_log((_S1242)));
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1243, float _S1244)
{
    _d_sqrt_0(_S1243, _S1244);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_14, float _s_dOut_0)
{
    float _S1245 = (*dpx_14).primal_0.x;
    float _S1246 = (*dpx_14).primal_0.y;
    float _S1247 = (*dpx_14).primal_0.z;
    DiffPair_float_0 _S1248;
    (&_S1248)->primal_0 = _S1245 * _S1245 + _S1246 * _S1246 + _S1247 * _S1247;
    (&_S1248)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1248, _s_dOut_0);
    float _S1249 = (*dpx_14).primal_0.z * _S1248.differential_0;
    float _S1250 = _S1249 + _S1249;
    float _S1251 = (*dpx_14).primal_0.y * _S1248.differential_0;
    float _S1252 = _S1251 + _S1251;
    float _S1253 = (*dpx_14).primal_0.x * _S1248.differential_0;
    float _S1254 = _S1253 + _S1253;
    float3  _S1255 = make_float3 (0.0f);
    *&((&_S1255)->z) = _S1250;
    *&((&_S1255)->y) = _S1252;
    *&((&_S1255)->x) = _S1254;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S1255;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1256, float _S1257)
{
    s_bwd_prop_length_impl_0(_S1256, _S1257);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S1258, float _S1259)
{
    _d_exp_0(_S1258, _S1259);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S1260, float _S1261)
{
    _d_log_0(_S1260, _S1261);
    return;
}

inline __device__ void s_bwd_prop_view_radius_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dplog_scale_0, DiffPair_float_0 * dplogit_opacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpcampos_0, float _s_dOut_1)
{
    float _S1262 = - (*dplogit_opacity_0).primal_0;
    float _S1263 = 1.0f + s_primal_ctx_exp_0(_S1262);
    float _S1264 = 255.0f / _S1263;
    float _S1265 = _S1263 * _S1263;
    float _S1266 = (F32_max((_S1264), (1.0f)));
    float _S1267 = 2.0f * s_primal_ctx_log_0(_S1266);
    float _S1268 = s_primal_ctx_sqrt_0(_S1267);
    float _S1269 = (*dplog_scale_0).primal_0.x;
    float _S1270 = (*dplog_scale_0).primal_0.y;
    float _S1271 = (*dplog_scale_0).primal_0.z;
    float _S1272 = (F32_max((_S1270), (_S1271)));
    float _S1273 = (F32_max((_S1269), (_S1272)));
    float _S1274 = s_primal_ctx_exp_0(_S1273);
    float radius_7 = _S1274 * _S1268;
    float3  _S1275 = (*dpmean_0).primal_0 - (*dpcampos_0).primal_0;
    float _S1276 = length_1(_S1275);
    float _S1277 = _S1276 * _S1276 - radius_7 * radius_7;
    float _S1278 = (F32_max((_S1277), (0.0f)));
    float _S1279 = (F32_max((_S1276), (radius_7))) + s_primal_ctx_sqrt_0(_S1278);
    float _S1280 = _s_dOut_1 / (_S1279 * _S1279);
    float _S1281 = radius_7 * - _S1280;
    float _S1282 = _S1279 * _S1280;
    DiffPair_float_0 _S1283;
    (&_S1283)->primal_0 = _S1278;
    (&_S1283)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1283, _S1281);
    DiffPair_float_0 _S1284;
    (&_S1284)->primal_0 = _S1277;
    (&_S1284)->differential_0 = 0.0f;
    DiffPair_float_0 _S1285;
    (&_S1285)->primal_0 = 0.0f;
    (&_S1285)->differential_0 = 0.0f;
    _d_max_0(&_S1284, &_S1285, _S1283.differential_0);
    float _S1286 = radius_7 * - _S1284.differential_0;
    float _S1287 = _S1276 * _S1284.differential_0;
    DiffPair_float_0 _S1288;
    (&_S1288)->primal_0 = _S1276;
    (&_S1288)->differential_0 = 0.0f;
    DiffPair_float_0 _S1289;
    (&_S1289)->primal_0 = radius_7;
    (&_S1289)->differential_0 = 0.0f;
    _d_max_0(&_S1288, &_S1289, _S1281);
    float _S1290 = _S1287 + _S1287 + _S1288.differential_0;
    float3  _S1291 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1292;
    (&_S1292)->primal_0 = _S1275;
    (&_S1292)->differential_0 = _S1291;
    s_bwd_length_impl_0(&_S1292, _S1290);
    float3  _S1293 = - _S1292.differential_0;
    float _S1294 = _S1282 + _S1286 + _S1286 + _S1289.differential_0;
    float _S1295 = _S1274 * _S1294;
    float _S1296 = _S1268 * _S1294;
    DiffPair_float_0 _S1297;
    (&_S1297)->primal_0 = _S1273;
    (&_S1297)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1297, _S1296);
    DiffPair_float_0 _S1298;
    (&_S1298)->primal_0 = _S1269;
    (&_S1298)->differential_0 = 0.0f;
    DiffPair_float_0 _S1299;
    (&_S1299)->primal_0 = _S1272;
    (&_S1299)->differential_0 = 0.0f;
    _d_max_0(&_S1298, &_S1299, _S1297.differential_0);
    DiffPair_float_0 _S1300;
    (&_S1300)->primal_0 = _S1270;
    (&_S1300)->differential_0 = 0.0f;
    DiffPair_float_0 _S1301;
    (&_S1301)->primal_0 = _S1271;
    (&_S1301)->differential_0 = 0.0f;
    _d_max_0(&_S1300, &_S1301, _S1299.differential_0);
    DiffPair_float_0 _S1302;
    (&_S1302)->primal_0 = _S1267;
    (&_S1302)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1302, _S1295);
    float _S1303 = 2.0f * _S1302.differential_0;
    DiffPair_float_0 _S1304;
    (&_S1304)->primal_0 = _S1266;
    (&_S1304)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1304, _S1303);
    DiffPair_float_0 _S1305;
    (&_S1305)->primal_0 = _S1264;
    (&_S1305)->differential_0 = 0.0f;
    DiffPair_float_0 _S1306;
    (&_S1306)->primal_0 = 1.0f;
    (&_S1306)->differential_0 = 0.0f;
    _d_max_0(&_S1305, &_S1306, _S1304.differential_0);
    float _S1307 = 255.0f * - (_S1305.differential_0 / _S1265);
    DiffPair_float_0 _S1308;
    (&_S1308)->primal_0 = _S1262;
    (&_S1308)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1308, _S1307);
    float _S1309 = - _S1308.differential_0;
    dpcampos_0->primal_0 = (*dpcampos_0).primal_0;
    dpcampos_0->differential_0 = _S1293;
    dplogit_opacity_0->primal_0 = (*dplogit_opacity_0).primal_0;
    dplogit_opacity_0->differential_0 = _S1309;
    float3  _S1310 = make_float3 (_S1298.differential_0, _S1300.differential_0, _S1301.differential_0);
    dplog_scale_0->primal_0 = (*dplog_scale_0).primal_0;
    dplog_scale_0->differential_0 = _S1310;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S1292.differential_0;
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S1311, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S1312, Matrix<float, 2, 2>  _S1313)
{
    mul_1(_S1311, _S1312, _S1313);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S1314, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1315, Matrix<float, 2, 3>  _S1316)
{
    mul_0(_S1314, _S1315, _S1316);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S1317, DiffPair_float_0 * _S1318, DiffPair_float_0 * _S1319, float _S1320)
{
    _d_clamp_0(_S1317, _S1318, _S1319, _S1320);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1321, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1322, Matrix<float, 3, 3>  _S1323)
{
    mul_2(_S1321, _S1322, _S1323);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1324, float3  _S1325)
{
    _d_exp_vector_0(_S1324, _S1325);
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_15, float _s_dOut_2)
{
    float _S1326 = (*dpx_15).primal_0.x;
    float _S1327 = (*dpx_15).primal_0.y;
    float _S1328 = (*dpx_15).primal_0.z;
    float _S1329 = (*dpx_15).primal_0.w;
    DiffPair_float_0 _S1330;
    (&_S1330)->primal_0 = _S1326 * _S1326 + _S1327 * _S1327 + _S1328 * _S1328 + _S1329 * _S1329;
    (&_S1330)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1330, _s_dOut_2);
    float _S1331 = (*dpx_15).primal_0.w * _S1330.differential_0;
    float _S1332 = _S1331 + _S1331;
    float _S1333 = (*dpx_15).primal_0.z * _S1330.differential_0;
    float _S1334 = _S1333 + _S1333;
    float _S1335 = (*dpx_15).primal_0.y * _S1330.differential_0;
    float _S1336 = _S1335 + _S1335;
    float _S1337 = (*dpx_15).primal_0.x * _S1330.differential_0;
    float _S1338 = _S1337 + _S1337;
    float4  _S1339 = make_float4 (0.0f);
    *&((&_S1339)->w) = _S1332;
    *&((&_S1339)->z) = _S1334;
    *&((&_S1339)->y) = _S1336;
    *&((&_S1339)->x) = _S1338;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S1339;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1340, float _S1341)
{
    s_bwd_prop_length_impl_1(_S1340, _S1341);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpx_16, float4  _s_dOut_3)
{
    float _S1342 = length_2((*dpx_16).primal_0);
    float4  _S1343 = (*dpx_16).primal_0 * _s_dOut_3;
    float4  _S1344 = make_float4 (1.0f / _S1342) * _s_dOut_3;
    float _S1345 = - ((_S1343.x + _S1343.y + _S1343.z + _S1343.w) / (_S1342 * _S1342));
    float4  _S1346 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1347;
    (&_S1347)->primal_0 = (*dpx_16).primal_0;
    (&_S1347)->differential_0 = _S1346;
    s_bwd_length_impl_1(&_S1347, _S1345);
    float4  _S1348 = _S1344 + _S1347.differential_0;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S1348;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1349, float4  _S1350)
{
    s_bwd_prop_normalize_impl_0(_S1349, _S1350);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S1351, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1352, float3  _S1353)
{
    _d_mul_0(_S1351, _S1352, _S1353);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_7, float4  quat_6, float3  scale_6, float in_opacity_6, Matrix<float, 3, 3>  R_6, float3  t_6, float fx_9, float fy_9, float cx_7, float cy_7, FixedArray<float, 10>  dist_coeffs_9, uint image_width_6, uint image_height_6, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_5 = s_primal_ctx_mul_0(R_6, mean_7) + t_6;
    float _S1354 = mean_c_5.z;
    float _S1355 = - in_opacity_6;
    float _S1356 = 1.0f + s_primal_ctx_exp_0(_S1355);
    float _S1357 = 1.0f / _S1356;
    float _S1358 = _S1356 * _S1356;
    float4  _S1359 = normalize_0(quat_6);
    float3  _S1360 = s_primal_ctx_exp_1(scale_6);
    float _S1361 = _S1359.y;
    float x2_6 = _S1361 * _S1361;
    float y2_6 = _S1359.z * _S1359.z;
    float z2_6 = _S1359.w * _S1359.w;
    float xy_6 = _S1359.y * _S1359.z;
    float xz_6 = _S1359.y * _S1359.w;
    float yz_6 = _S1359.z * _S1359.w;
    float wx_6 = _S1359.x * _S1359.y;
    float wy_6 = _S1359.x * _S1359.z;
    float wz_6 = _S1359.x * _S1359.w;
    Matrix<float, 3, 3>  _S1362 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S1360.x, 0.0f, 0.0f, 0.0f, _S1360.y, 0.0f, 0.0f, 0.0f, _S1360.z);
    Matrix<float, 3, 3>  _S1363 = s_primal_ctx_mul_1(_S1362, S_0);
    Matrix<float, 3, 3>  _S1364 = transpose_3(_S1363);
    Matrix<float, 3, 3>  _S1365 = s_primal_ctx_mul_1(_S1363, _S1364);
    Matrix<float, 3, 3>  _S1366 = s_primal_ctx_mul_1(R_6, _S1365);
    Matrix<float, 3, 3>  _S1367 = transpose_3(R_6);
    Matrix<float, 3, 3>  _S1368 = s_primal_ctx_mul_1(_S1366, _S1367);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S1369 = float2 {mean_c_5.x, mean_c_5.y};
    float2  _S1370 = make_float2 (1.0f, 0.0f);
    float2  _S1371 = make_float2 (_S1354);
    float2  _S1372 = _S1369 / make_float2 (_S1354);
    float _S1373 = _S1354 * _S1354;
    float2  _S1374 = make_float2 (_S1373);
    float2  _S1375 = _S1370 * make_float2 (_S1354);
    float2  _S1376 = _S1375 / make_float2 (_S1373);
    float2  _S1377 = make_float2 (_S1373 * _S1373);
    float u_72 = _S1372.x;
    float s_diff_u_18 = _S1376.x;
    float v_72 = _S1372.y;
    float s_diff_v_18 = _S1376.y;
    float _S1378 = s_diff_u_18 * u_72;
    float _S1379 = s_diff_v_18 * v_72;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float s_diff_r2_72 = _S1378 + _S1378 + (_S1379 + _S1379);
    float _S1380 = s_diff_r2_72 * dist_coeffs_9[int(3)];
    float _S1381 = dist_coeffs_9[int(2)] + r2_72 * dist_coeffs_9[int(3)];
    float _S1382 = s_diff_r2_72 * _S1381 + _S1380 * r2_72;
    float _S1383 = dist_coeffs_9[int(1)] + r2_72 * _S1381;
    float _S1384 = s_diff_r2_72 * _S1383 + _S1382 * r2_72;
    float _S1385 = dist_coeffs_9[int(0)] + r2_72 * _S1383;
    float _S1386 = s_diff_r2_72 * _S1385 + _S1384 * r2_72;
    float2  _S1387 = make_float2 (_S1386);
    float radial_27 = 1.0f + r2_72 * _S1385;
    float2  _S1388 = make_float2 (radial_27);
    float _S1389 = 2.0f * dist_coeffs_9[int(4)];
    float _S1390 = _S1389 * u_72;
    float _S1391 = s_diff_u_18 * _S1389;
    float _S1392 = 2.0f * u_72;
    float _S1393 = s_diff_u_18 * 2.0f;
    float _S1394 = 2.0f * dist_coeffs_9[int(5)];
    float _S1395 = _S1394 * u_72;
    float _S1396 = s_diff_u_18 * _S1394;
    float _S1397 = 2.0f * v_72;
    float _S1398 = s_diff_v_18 * 2.0f;
    float2  _S1399 = _S1376 * make_float2 (radial_27) + make_float2 (_S1386) * _S1372 + make_float2 (_S1391 * v_72 + s_diff_v_18 * _S1390 + (s_diff_r2_72 + (_S1393 * u_72 + s_diff_u_18 * _S1392)) * dist_coeffs_9[int(5)] + s_diff_r2_72 * dist_coeffs_9[int(6)], _S1396 * v_72 + s_diff_v_18 * _S1395 + (s_diff_r2_72 + (_S1398 * v_72 + s_diff_v_18 * _S1397)) * dist_coeffs_9[int(4)] + s_diff_r2_72 * dist_coeffs_9[int(7)]);
    float2  _S1400 = _S1399 + make_float2 (_S1399.x * dist_coeffs_9[int(8)] + _S1399.y * dist_coeffs_9[int(9)], 0.0f);
    float _S1401 = _S1400.x * fx_9;
    float _S1402 = _S1400.y * fy_9;
    Matrix<float, 2, 3>  _S1403 = J_8;
    *&(((&_S1403)->rows + (int(0)))->x) = _S1401;
    *&(((&_S1403)->rows + (int(1)))->x) = _S1402;
    float2  _S1404 = make_float2 (0.0f, 1.0f);
    float2  _S1405 = _S1369 / make_float2 (_S1354);
    float2  _S1406 = _S1404 * make_float2 (_S1354);
    float2  _S1407 = _S1406 / make_float2 (_S1373);
    float u_73 = _S1405.x;
    float s_diff_u_19 = _S1407.x;
    float v_73 = _S1405.y;
    float s_diff_v_19 = _S1407.y;
    float _S1408 = s_diff_u_19 * u_73;
    float _S1409 = s_diff_v_19 * v_73;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float s_diff_r2_73 = _S1408 + _S1408 + (_S1409 + _S1409);
    float _S1410 = s_diff_r2_73 * dist_coeffs_9[int(3)];
    float _S1411 = dist_coeffs_9[int(2)] + r2_73 * dist_coeffs_9[int(3)];
    float _S1412 = s_diff_r2_73 * _S1411 + _S1410 * r2_73;
    float _S1413 = dist_coeffs_9[int(1)] + r2_73 * _S1411;
    float _S1414 = s_diff_r2_73 * _S1413 + _S1412 * r2_73;
    float _S1415 = dist_coeffs_9[int(0)] + r2_73 * _S1413;
    float _S1416 = s_diff_r2_73 * _S1415 + _S1414 * r2_73;
    float2  _S1417 = make_float2 (_S1416);
    float radial_28 = 1.0f + r2_73 * _S1415;
    float2  _S1418 = make_float2 (radial_28);
    float _S1419 = _S1389 * u_73;
    float _S1420 = s_diff_u_19 * _S1389;
    float _S1421 = 2.0f * u_73;
    float _S1422 = s_diff_u_19 * 2.0f;
    float _S1423 = _S1394 * u_73;
    float _S1424 = s_diff_u_19 * _S1394;
    float _S1425 = 2.0f * v_73;
    float _S1426 = s_diff_v_19 * 2.0f;
    float2  _S1427 = _S1407 * make_float2 (radial_28) + make_float2 (_S1416) * _S1405 + make_float2 (_S1420 * v_73 + s_diff_v_19 * _S1419 + (s_diff_r2_73 + (_S1422 * u_73 + s_diff_u_19 * _S1421)) * dist_coeffs_9[int(5)] + s_diff_r2_73 * dist_coeffs_9[int(6)], _S1424 * v_73 + s_diff_v_19 * _S1423 + (s_diff_r2_73 + (_S1426 * v_73 + s_diff_v_19 * _S1425)) * dist_coeffs_9[int(4)] + s_diff_r2_73 * dist_coeffs_9[int(7)]);
    float2  _S1428 = _S1427 + make_float2 (_S1427.x * dist_coeffs_9[int(8)] + _S1427.y * dist_coeffs_9[int(9)], 0.0f);
    float _S1429 = _S1428.y * fy_9;
    *&(((&_S1403)->rows + (int(0)))->y) = _S1428.x * fx_9;
    *&(((&_S1403)->rows + (int(1)))->y) = _S1429;
    float2  _S1430 = _S1369 / make_float2 (_S1354);
    float2  _S1431 = make_float2 (0.0f, 0.0f) - _S1369;
    float2  _S1432 = _S1431 / make_float2 (_S1373);
    float u_74 = _S1430.x;
    float s_diff_u_20 = _S1432.x;
    float v_74 = _S1430.y;
    float s_diff_v_20 = _S1432.y;
    float _S1433 = s_diff_u_20 * u_74;
    float _S1434 = s_diff_v_20 * v_74;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float s_diff_r2_74 = _S1433 + _S1433 + (_S1434 + _S1434);
    float _S1435 = s_diff_r2_74 * dist_coeffs_9[int(3)];
    float _S1436 = dist_coeffs_9[int(2)] + r2_74 * dist_coeffs_9[int(3)];
    float _S1437 = s_diff_r2_74 * _S1436 + _S1435 * r2_74;
    float _S1438 = dist_coeffs_9[int(1)] + r2_74 * _S1436;
    float _S1439 = s_diff_r2_74 * _S1438 + _S1437 * r2_74;
    float _S1440 = dist_coeffs_9[int(0)] + r2_74 * _S1438;
    float _S1441 = s_diff_r2_74 * _S1440 + _S1439 * r2_74;
    float2  _S1442 = make_float2 (_S1441);
    float radial_29 = 1.0f + r2_74 * _S1440;
    float2  _S1443 = make_float2 (radial_29);
    float _S1444 = _S1389 * u_74;
    float _S1445 = s_diff_u_20 * _S1389;
    float _S1446 = 2.0f * u_74;
    float _S1447 = s_diff_u_20 * 2.0f;
    float _S1448 = _S1394 * u_74;
    float _S1449 = s_diff_u_20 * _S1394;
    float _S1450 = 2.0f * v_74;
    float _S1451 = s_diff_v_20 * 2.0f;
    float2  _S1452 = _S1432 * make_float2 (radial_29) + make_float2 (_S1441) * _S1430 + make_float2 (_S1445 * v_74 + s_diff_v_20 * _S1444 + (s_diff_r2_74 + (_S1447 * u_74 + s_diff_u_20 * _S1446)) * dist_coeffs_9[int(5)] + s_diff_r2_74 * dist_coeffs_9[int(6)], _S1449 * v_74 + s_diff_v_20 * _S1448 + (s_diff_r2_74 + (_S1451 * v_74 + s_diff_v_20 * _S1450)) * dist_coeffs_9[int(4)] + s_diff_r2_74 * dist_coeffs_9[int(7)]);
    float2  _S1453 = _S1452 + make_float2 (_S1452.x * dist_coeffs_9[int(8)] + _S1452.y * dist_coeffs_9[int(9)], 0.0f);
    float _S1454 = _S1453.x * fx_9;
    float _S1455 = _S1453.y * fy_9;
    float _S1456 = float(image_width_6);
    float _S1457 = 0.30000001192092896f * (0.5f * _S1456);
    float lim_x_pos_3 = _S1456 + _S1457;
    float rz_2 = 1.0f / _S1354;
    float _S1458 = - _S1457;
    float _S1459 = - (_S1458 - cx_7);
    float max_Jxz_0 = _S1459 * rz_2;
    float _S1460 = - (lim_x_pos_3 - cx_7);
    float min_Jxz_0 = _S1460 * rz_2;
    float _S1461 = - (_S1458 - cy_7);
    float max_Jyz_2 = _S1461 * rz_2;
    float _S1462 = - (lim_x_pos_3 - cy_7);
    float min_Jyz_2 = _S1462 * rz_2;
    *&(((&_S1403)->rows + (int(0)))->z) = s_primal_ctx_clamp_0(_S1454, min_Jxz_0, max_Jxz_0);
    *&(((&_S1403)->rows + (int(1)))->z) = s_primal_ctx_clamp_0(_S1455, min_Jyz_2, max_Jyz_2);
    Matrix<float, 2, 3>  _S1463 = s_primal_ctx_mul_2(_S1403, _S1368);
    Matrix<float, 3, 2>  _S1464 = transpose_1(_S1403);
    Matrix<float, 2, 2>  _S1465 = s_primal_ctx_mul_3(_S1463, _S1464);
    float eps2d_6;
    if(antialiased_6)
    {
        eps2d_6 = 0.10000000149011612f;
    }
    else
    {
        eps2d_6 = 0.30000001192092896f;
    }
    float _S1466 = _S1465.rows[int(0)].y * _S1465.rows[int(1)].x;
    float det_orig_6 = _S1465.rows[int(0)].x * _S1465.rows[int(1)].y - _S1466;
    float _S1467 = _S1465.rows[int(0)].x + eps2d_6;
    Matrix<float, 2, 2>  _S1468 = _S1465;
    *&(((&_S1468)->rows + (int(0)))->x) = _S1467;
    float _S1469 = _S1465.rows[int(1)].y + eps2d_6;
    *&(((&_S1468)->rows + (int(1)))->y) = _S1469;
    Matrix<float, 2, 2>  _S1470 = _S1468;
    Matrix<float, 2, 2>  _S1471 = _S1468;
    float det_blur_6 = _S1467 * _S1469 - _S1466;
    float _S1472 = det_orig_6 / det_blur_6;
    float _S1473 = det_blur_6 * det_blur_6;
    float _S1474 = (F32_max((0.0f), (_S1472)));
    float _S1475 = s_primal_ctx_sqrt_0(_S1474);
    float invdet_8 = 1.0f / det_blur_6;
    float _S1476 = - _S1465.rows[int(0)].y;
    float _S1477 = - _S1465.rows[int(1)].x;
    if(antialiased_6)
    {
        eps2d_6 = _S1357 * _S1475;
    }
    else
    {
        eps2d_6 = _S1357;
    }
    float _S1478 = eps2d_6 / 0.00392156885936856f;
    float _S1479 = 2.0f * s_primal_ctx_log_0(_S1478);
    float _S1480 = s_primal_ctx_sqrt_0(_S1479);
    float _S1481 = _S1470.rows[int(0)].x;
    float _S1482 = _S1471.rows[int(1)].y;
    float3  campos_1 = - s_primal_ctx_mul_0(_S1367, t_6);
    float3  _S1483 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1484;
    (&_S1484)->primal_0 = mean_7;
    (&_S1484)->differential_0 = _S1483;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1485;
    (&_S1485)->primal_0 = scale_6;
    (&_S1485)->differential_0 = _S1483;
    DiffPair_float_0 _S1486;
    (&_S1486)->primal_0 = in_opacity_6;
    (&_S1486)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1487;
    (&_S1487)->primal_0 = campos_1;
    (&_S1487)->differential_0 = _S1483;
    s_bwd_prop_view_radius_3dgs_0(&_S1484, &_S1485, &_S1486, &_S1487, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1488 = _S1484;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1489 = _S1485;
    DiffPair_float_0 _S1490 = _S1486;
    float2  _S1491 = make_float2 (0.0f);
    float2  _S1492 = _S1491;
    *&((&_S1492)->y) = v_conic_0.z;
    float2  _S1493 = _S1491;
    *&((&_S1493)->y) = v_conic_0.y;
    *&((&_S1493)->x) = v_conic_0.x;
    DiffPair_float_0 _S1494;
    (&_S1494)->primal_0 = _S1482;
    (&_S1494)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1494, 0.0f);
    DiffPair_float_0 _S1495;
    (&_S1495)->primal_0 = _S1481;
    (&_S1495)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1495, 0.0f);
    DiffPair_float_0 _S1496;
    (&_S1496)->primal_0 = 3.32999992370605469f;
    (&_S1496)->differential_0 = 0.0f;
    DiffPair_float_0 _S1497;
    (&_S1497)->primal_0 = _S1480;
    (&_S1497)->differential_0 = 0.0f;
    _d_min_0(&_S1496, &_S1497, 0.0f);
    DiffPair_float_0 _S1498;
    (&_S1498)->primal_0 = _S1479;
    (&_S1498)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1498, _S1497.differential_0);
    float _S1499 = 2.0f * _S1498.differential_0;
    DiffPair_float_0 _S1500;
    (&_S1500)->primal_0 = _S1478;
    (&_S1500)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1500, _S1499);
    float _S1501 = v_opacity_0 + 254.9999847412109375f * _S1500.differential_0;
    float2  _S1502 = make_float2 (_S1495.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S1503 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1504 = _S1503;
    _S1504[int(1)] = _S1492;
    _S1504[int(0)] = _S1493;
    Matrix<float, 2, 2>  _S1505 = _S1504;
    float2  _S1506 = make_float2 (0.0f, _S1494.differential_0);
    float _S1507;
    if(antialiased_6)
    {
        float _S1508 = _S1475 * _S1501;
        eps2d_6 = _S1357 * _S1501;
        _S1507 = _S1508;
    }
    else
    {
        eps2d_6 = 0.0f;
        _S1507 = _S1501;
    }
    float _S1509 = invdet_8 * _S1505.rows[int(1)].y;
    float _S1510 = - (invdet_8 * _S1505.rows[int(1)].x);
    float _S1511 = - (invdet_8 * _S1505.rows[int(0)].y);
    float _S1512 = invdet_8 * _S1505.rows[int(0)].x;
    float _S1513 = - ((_S1467 * _S1505.rows[int(1)].y + _S1477 * _S1505.rows[int(1)].x + _S1476 * _S1505.rows[int(0)].y + _S1469 * _S1505.rows[int(0)].x) / _S1473);
    DiffPair_float_0 _S1514;
    (&_S1514)->primal_0 = _S1474;
    (&_S1514)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1514, eps2d_6);
    DiffPair_float_0 _S1515;
    (&_S1515)->primal_0 = 0.0f;
    (&_S1515)->differential_0 = 0.0f;
    DiffPair_float_0 _S1516;
    (&_S1516)->primal_0 = _S1472;
    (&_S1516)->differential_0 = 0.0f;
    _d_max_0(&_S1515, &_S1516, _S1514.differential_0);
    float _S1517 = _S1516.differential_0 / _S1473;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1517;
    float _S1518 = det_orig_6 * - _S1517 + _S1513;
    float _S1519 = - _S1518;
    float _S1520 = _S1467 * _S1518;
    float _S1521 = _S1469 * _S1518;
    Matrix<float, 2, 2>  _S1522 = _S1503;
    _S1522[int(1)] = _S1506;
    _S1522[int(0)] = _S1502;
    _S1468 = _S1522;
    *&(((&_S1468)->rows + (int(1)))->y) = 0.0f;
    float _S1523 = _S1520 + _S1522.rows[int(1)].y + _S1512;
    *&(((&_S1468)->rows + (int(0)))->x) = 0.0f;
    float _S1524 = _S1521 + _S1522.rows[int(0)].x + _S1509;
    float _S1525 = _S1519 + - s_diff_det_orig_T_0;
    float _S1526 = _S1465.rows[int(0)].y * _S1525 + _S1510;
    float _S1527 = _S1465.rows[int(1)].x * _S1525 + _S1511;
    float _S1528 = _S1465.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1529 = _S1523 + _S1465.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1530 = _S1491;
    *&((&_S1530)->x) = _S1526;
    *&((&_S1530)->y) = _S1529;
    float _S1531 = _S1524 + _S1528;
    float2  _S1532 = _S1491;
    *&((&_S1532)->y) = _S1527;
    *&((&_S1532)->x) = _S1531;
    Matrix<float, 2, 2>  _S1533 = _S1503;
    _S1533[int(1)] = _S1530;
    _S1533[int(0)] = _S1532;
    Matrix<float, 2, 2>  _S1534 = _S1468 + _S1533;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1535;
    (&_S1535)->primal_0 = _S1463;
    (&_S1535)->differential_0 = J_8;
    Matrix<float, 3, 2>  _S1536 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1537;
    (&_S1537)->primal_0 = _S1464;
    (&_S1537)->differential_0 = _S1536;
    s_bwd_prop_mul_0(&_S1535, &_S1537, _S1534);
    Matrix<float, 2, 3>  _S1538 = transpose_2(_S1537.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1539;
    (&_S1539)->primal_0 = _S1403;
    (&_S1539)->differential_0 = J_8;
    Matrix<float, 3, 3>  _S1540 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1541;
    (&_S1541)->primal_0 = _S1368;
    (&_S1541)->differential_0 = _S1540;
    s_bwd_prop_mul_1(&_S1539, &_S1541, _S1535.differential_0);
    Matrix<float, 2, 3>  _S1542 = _S1538 + _S1539.differential_0;
    DiffPair_float_0 _S1543;
    (&_S1543)->primal_0 = _S1455;
    (&_S1543)->differential_0 = 0.0f;
    DiffPair_float_0 _S1544;
    (&_S1544)->primal_0 = min_Jyz_2;
    (&_S1544)->differential_0 = 0.0f;
    DiffPair_float_0 _S1545;
    (&_S1545)->primal_0 = max_Jyz_2;
    (&_S1545)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1543, &_S1544, &_S1545, _S1542.rows[int(1)].z);
    DiffPair_float_0 _S1546;
    (&_S1546)->primal_0 = _S1454;
    (&_S1546)->differential_0 = 0.0f;
    DiffPair_float_0 _S1547;
    (&_S1547)->primal_0 = min_Jxz_0;
    (&_S1547)->differential_0 = 0.0f;
    DiffPair_float_0 _S1548;
    (&_S1548)->primal_0 = max_Jxz_0;
    (&_S1548)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1546, &_S1547, &_S1548, _S1542.rows[int(0)].z);
    float _S1549 = - ((_S1462 * _S1544.differential_0 + _S1461 * _S1545.differential_0 + _S1460 * _S1547.differential_0 + _S1459 * _S1548.differential_0) / _S1373);
    float2  _S1550 = make_float2 (0.0f, _S1543.differential_0) + make_float2 (_S1546.differential_0, 0.0f);
    float _S1551 = fx_9 * _S1550.x;
    float2  _S1552 = make_float2 (_S1551, fy_9 * _S1550.y) + make_float2 (dist_coeffs_9[int(8)] * _S1551, dist_coeffs_9[int(9)] * _S1551);
    float2  _S1553 = _S1430 * _S1552;
    float2  _S1554 = _S1432 * _S1552;
    float _S1555 = dist_coeffs_9[int(4)] * _S1552.y;
    float _S1556 = dist_coeffs_9[int(5)] * _S1552.x;
    float _S1557 = _S1554.x + _S1554.y;
    float _S1558 = _S1553.x + _S1553.y;
    float _S1559 = r2_74 * _S1558;
    float _S1560 = s_diff_r2_74 * _S1558 + r2_74 * _S1557;
    float _S1561 = r2_74 * _S1559;
    float _S1562 = s_diff_r2_74 * _S1559 + r2_74 * _S1560;
    float _S1563 = dist_coeffs_9[int(7)] * _S1552.y + _S1555 + dist_coeffs_9[int(6)] * _S1552.x + _S1556 + _S1440 * _S1558 + _S1438 * _S1559 + _S1436 * _S1561 + dist_coeffs_9[int(3)] * (r2_74 * _S1561);
    float _S1564 = _S1439 * _S1558 + _S1440 * _S1557 + _S1437 * _S1559 + _S1438 * _S1560 + _S1435 * _S1561 + _S1436 * _S1562 + dist_coeffs_9[int(3)] * (s_diff_r2_74 * _S1561 + r2_74 * _S1562);
    float _S1565 = _S1563 + _S1563;
    float _S1566 = v_74 * _S1564;
    float _S1567 = u_74 * _S1564;
    float2  _S1568 = (_S1443 * _S1552 + make_float2 (_S1394 * (v_74 * _S1552.y) + _S1446 * _S1556 + 2.0f * (u_74 * _S1556) + _S1389 * (v_74 * _S1552.x) + u_74 * _S1565, _S1450 * _S1555 + 2.0f * (v_74 * _S1555) + _S1448 * _S1552.y + _S1444 * _S1552.x + v_74 * _S1565)) / _S1377;
    float2  _S1569 = _S1431 * - _S1568;
    float _S1570 = _S1354 * (_S1569.x + _S1569.y);
    float2  _S1571 = (_S1442 * _S1552 + make_float2 (_S1394 * (s_diff_v_20 * _S1552.y) + _S1447 * _S1556 + 2.0f * (s_diff_u_20 * _S1556) + _S1389 * (s_diff_v_20 * _S1552.x) + s_diff_u_20 * _S1565 + _S1567 + _S1567, _S1451 * _S1555 + 2.0f * (s_diff_v_20 * _S1555) + _S1449 * _S1552.y + _S1445 * _S1552.x + s_diff_v_20 * _S1565 + _S1566 + _S1566)) / _S1374;
    float2  _S1572 = _S1369 * - _S1571;
    float2  _S1573 = - (_S1374 * _S1568) + _S1371 * _S1571;
    float3  _S1574 = make_float3 (_S1573.x, _S1573.y, _S1570 + _S1570 + _S1572.x + _S1572.y);
    float2  _S1575 = make_float2 (0.0f, _S1542.rows[int(1)].y) + make_float2 (_S1542.rows[int(0)].y, 0.0f);
    float _S1576 = fx_9 * _S1575.x;
    float2  _S1577 = make_float2 (_S1576, fy_9 * _S1575.y) + make_float2 (dist_coeffs_9[int(8)] * _S1576, dist_coeffs_9[int(9)] * _S1576);
    float2  _S1578 = _S1405 * _S1577;
    float2  _S1579 = _S1407 * _S1577;
    float _S1580 = dist_coeffs_9[int(4)] * _S1577.y;
    float _S1581 = dist_coeffs_9[int(5)] * _S1577.x;
    float _S1582 = _S1579.x + _S1579.y;
    float _S1583 = _S1578.x + _S1578.y;
    float _S1584 = r2_73 * _S1583;
    float _S1585 = s_diff_r2_73 * _S1583 + r2_73 * _S1582;
    float _S1586 = r2_73 * _S1584;
    float _S1587 = s_diff_r2_73 * _S1584 + r2_73 * _S1585;
    float _S1588 = dist_coeffs_9[int(7)] * _S1577.y + _S1580 + dist_coeffs_9[int(6)] * _S1577.x + _S1581 + _S1415 * _S1583 + _S1413 * _S1584 + _S1411 * _S1586 + dist_coeffs_9[int(3)] * (r2_73 * _S1586);
    float _S1589 = _S1414 * _S1583 + _S1415 * _S1582 + _S1412 * _S1584 + _S1413 * _S1585 + _S1410 * _S1586 + _S1411 * _S1587 + dist_coeffs_9[int(3)] * (s_diff_r2_73 * _S1586 + r2_73 * _S1587);
    float _S1590 = _S1588 + _S1588;
    float _S1591 = v_73 * _S1589;
    float _S1592 = u_73 * _S1589;
    float2  _S1593 = (_S1418 * _S1577 + make_float2 (_S1394 * (v_73 * _S1577.y) + _S1421 * _S1581 + 2.0f * (u_73 * _S1581) + _S1389 * (v_73 * _S1577.x) + u_73 * _S1590, _S1425 * _S1580 + 2.0f * (v_73 * _S1580) + _S1423 * _S1577.y + _S1419 * _S1577.x + v_73 * _S1590)) / _S1377;
    float2  _S1594 = _S1406 * - _S1593;
    float _S1595 = _S1354 * (_S1594.x + _S1594.y);
    float2  _S1596 = _S1404 * (_S1374 * _S1593);
    float2  _S1597 = (_S1417 * _S1577 + make_float2 (_S1394 * (s_diff_v_19 * _S1577.y) + _S1422 * _S1581 + 2.0f * (s_diff_u_19 * _S1581) + _S1389 * (s_diff_v_19 * _S1577.x) + s_diff_u_19 * _S1590 + _S1592 + _S1592, _S1426 * _S1580 + 2.0f * (s_diff_v_19 * _S1580) + _S1424 * _S1577.y + _S1420 * _S1577.x + s_diff_v_19 * _S1590 + _S1591 + _S1591)) / _S1374;
    float2  _S1598 = _S1369 * - _S1597;
    float2  _S1599 = _S1371 * _S1597;
    float3  _S1600 = make_float3 (_S1599.x, _S1599.y, _S1595 + _S1595 + _S1596.x + _S1596.y + _S1598.x + _S1598.y);
    float2  _S1601 = make_float2 (0.0f, _S1542.rows[int(1)].x) + make_float2 (_S1542.rows[int(0)].x, 0.0f);
    float _S1602 = fx_9 * _S1601.x;
    float2  _S1603 = make_float2 (_S1602, fy_9 * _S1601.y) + make_float2 (dist_coeffs_9[int(8)] * _S1602, dist_coeffs_9[int(9)] * _S1602);
    float2  _S1604 = _S1372 * _S1603;
    float2  _S1605 = _S1376 * _S1603;
    float _S1606 = dist_coeffs_9[int(4)] * _S1603.y;
    float _S1607 = dist_coeffs_9[int(5)] * _S1603.x;
    float _S1608 = _S1605.x + _S1605.y;
    float _S1609 = _S1604.x + _S1604.y;
    float _S1610 = r2_72 * _S1609;
    float _S1611 = s_diff_r2_72 * _S1609 + r2_72 * _S1608;
    float _S1612 = r2_72 * _S1610;
    float _S1613 = s_diff_r2_72 * _S1610 + r2_72 * _S1611;
    float _S1614 = dist_coeffs_9[int(7)] * _S1603.y + _S1606 + dist_coeffs_9[int(6)] * _S1603.x + _S1607 + _S1385 * _S1609 + _S1383 * _S1610 + _S1381 * _S1612 + dist_coeffs_9[int(3)] * (r2_72 * _S1612);
    float _S1615 = _S1384 * _S1609 + _S1385 * _S1608 + _S1382 * _S1610 + _S1383 * _S1611 + _S1380 * _S1612 + _S1381 * _S1613 + dist_coeffs_9[int(3)] * (s_diff_r2_72 * _S1612 + r2_72 * _S1613);
    float _S1616 = _S1614 + _S1614;
    float _S1617 = v_72 * _S1615;
    float _S1618 = u_72 * _S1615;
    float2  _S1619 = (_S1388 * _S1603 + make_float2 (_S1394 * (v_72 * _S1603.y) + _S1392 * _S1607 + 2.0f * (u_72 * _S1607) + _S1389 * (v_72 * _S1603.x) + u_72 * _S1616, _S1397 * _S1606 + 2.0f * (v_72 * _S1606) + _S1395 * _S1603.y + _S1390 * _S1603.x + v_72 * _S1616)) / _S1377;
    float2  _S1620 = _S1375 * - _S1619;
    float _S1621 = _S1354 * (_S1620.x + _S1620.y);
    float2  _S1622 = _S1370 * (_S1374 * _S1619);
    float2  _S1623 = (_S1387 * _S1603 + make_float2 (_S1394 * (s_diff_v_18 * _S1603.y) + _S1393 * _S1607 + 2.0f * (s_diff_u_18 * _S1607) + _S1389 * (s_diff_v_18 * _S1603.x) + s_diff_u_18 * _S1616 + _S1618 + _S1618, _S1398 * _S1606 + 2.0f * (s_diff_v_18 * _S1606) + _S1396 * _S1603.y + _S1391 * _S1603.x + s_diff_v_18 * _S1616 + _S1617 + _S1617)) / _S1374;
    float2  _S1624 = _S1369 * - _S1623;
    float2  _S1625 = _S1371 * _S1623;
    float3  _S1626 = make_float3 (_S1625.x, _S1625.y, _S1621 + _S1621 + _S1622.x + _S1622.y + _S1624.x + _S1624.y);
    float2  _S1627 = _S1369 / make_float2 (_S1354);
    float _S1628 = fx_9 * v_mean2d_0.x;
    float u_75 = _S1627.x;
    float v_75 = _S1627.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S1629 = dist_coeffs_9[int(2)] + r2_75 * dist_coeffs_9[int(3)];
    float _S1630 = dist_coeffs_9[int(1)] + r2_75 * _S1629;
    float _S1631 = dist_coeffs_9[int(0)] + r2_75 * _S1630;
    float2  _S1632 = make_float2 (_S1628, fy_9 * v_mean2d_0.y) + make_float2 (dist_coeffs_9[int(8)] * _S1628, dist_coeffs_9[int(9)] * _S1628);
    float2  _S1633 = _S1627 * _S1632;
    float _S1634 = dist_coeffs_9[int(4)] * _S1632.y;
    float _S1635 = dist_coeffs_9[int(5)] * _S1632.x;
    float _S1636 = _S1633.x + _S1633.y;
    float _S1637 = r2_75 * _S1636;
    float _S1638 = r2_75 * _S1637;
    float _S1639 = dist_coeffs_9[int(7)] * _S1632.y + _S1634 + dist_coeffs_9[int(6)] * _S1632.x + _S1635 + _S1631 * _S1636 + _S1630 * _S1637 + _S1629 * _S1638 + dist_coeffs_9[int(3)] * (r2_75 * _S1638);
    float _S1640 = v_75 * _S1639;
    float _S1641 = u_75 * _S1639;
    float2  _S1642 = (make_float2 (1.0f + r2_75 * _S1631) * _S1632 + make_float2 (_S1394 * (v_75 * _S1632.y) + 2.0f * u_75 * _S1635 + 2.0f * (u_75 * _S1635) + _S1389 * (v_75 * _S1632.x) + _S1641 + _S1641, 2.0f * v_75 * _S1634 + 2.0f * (v_75 * _S1634) + _S1394 * u_75 * _S1632.y + _S1389 * u_75 * _S1632.x + _S1640 + _S1640)) / _S1374;
    float2  _S1643 = _S1369 * - _S1642;
    float2  _S1644 = _S1371 * _S1642;
    float3  _S1645 = make_float3 (_S1644.x, _S1644.y, _S1643.x + _S1643.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1646;
    (&_S1646)->primal_0 = _S1366;
    (&_S1646)->differential_0 = _S1540;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1647;
    (&_S1647)->primal_0 = _S1367;
    (&_S1647)->differential_0 = _S1540;
    s_bwd_prop_mul_2(&_S1646, &_S1647, _S1541.differential_0);
    Matrix<float, 3, 3>  _S1648 = transpose_3(_S1647.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1649;
    (&_S1649)->primal_0 = R_6;
    (&_S1649)->differential_0 = _S1540;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1650;
    (&_S1650)->primal_0 = _S1365;
    (&_S1650)->differential_0 = _S1540;
    s_bwd_prop_mul_2(&_S1649, &_S1650, _S1646.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1651;
    (&_S1651)->primal_0 = _S1363;
    (&_S1651)->differential_0 = _S1540;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1652;
    (&_S1652)->primal_0 = _S1364;
    (&_S1652)->differential_0 = _S1540;
    s_bwd_prop_mul_2(&_S1651, &_S1652, _S1650.differential_0);
    Matrix<float, 3, 3>  _S1653 = _S1651.differential_0 + transpose_3(_S1652.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1654;
    (&_S1654)->primal_0 = _S1362;
    (&_S1654)->differential_0 = _S1540;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1655;
    (&_S1655)->primal_0 = S_0;
    (&_S1655)->differential_0 = _S1540;
    s_bwd_prop_mul_2(&_S1654, &_S1655, _S1653);
    Matrix<float, 3, 3>  _S1656 = transpose_3(_S1654.differential_0);
    float _S1657 = 2.0f * - _S1656.rows[int(2)].z;
    float _S1658 = 2.0f * _S1656.rows[int(2)].y;
    float _S1659 = 2.0f * _S1656.rows[int(2)].x;
    float _S1660 = 2.0f * _S1656.rows[int(1)].z;
    float _S1661 = 2.0f * - _S1656.rows[int(1)].y;
    float _S1662 = 2.0f * _S1656.rows[int(1)].x;
    float _S1663 = 2.0f * _S1656.rows[int(0)].z;
    float _S1664 = 2.0f * _S1656.rows[int(0)].y;
    float _S1665 = 2.0f * - _S1656.rows[int(0)].x;
    float _S1666 = - _S1662 + _S1664;
    float _S1667 = _S1659 + - _S1663;
    float _S1668 = - _S1658 + _S1660;
    float _S1669 = _S1658 + _S1660;
    float _S1670 = _S1659 + _S1663;
    float _S1671 = _S1662 + _S1664;
    float _S1672 = _S1359.w * (_S1661 + _S1665);
    float _S1673 = _S1359.z * (_S1657 + _S1665);
    float _S1674 = _S1359.y * (_S1657 + _S1661);
    float _S1675 = _S1359.x * _S1666 + _S1359.z * _S1669 + _S1359.y * _S1670 + _S1672 + _S1672;
    float _S1676 = _S1359.x * _S1667 + _S1359.w * _S1669 + _S1359.y * _S1671 + _S1673 + _S1673;
    float _S1677 = _S1359.x * _S1668 + _S1359.w * _S1670 + _S1359.z * _S1671 + _S1674 + _S1674;
    float _S1678 = _S1359.w * _S1666 + _S1359.z * _S1667 + _S1359.y * _S1668;
    float3  _S1679 = _S1483;
    *&((&_S1679)->z) = _S1655.differential_0.rows[int(2)].z;
    *&((&_S1679)->y) = _S1655.differential_0.rows[int(1)].y;
    *&((&_S1679)->x) = _S1655.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1680;
    (&_S1680)->primal_0 = scale_6;
    (&_S1680)->differential_0 = _S1483;
    s_bwd_prop_exp_1(&_S1680, _S1679);
    float4  _S1681 = make_float4 (0.0f);
    float4  _S1682 = _S1681;
    *&((&_S1682)->w) = _S1675;
    *&((&_S1682)->z) = _S1676;
    *&((&_S1682)->y) = _S1677;
    *&((&_S1682)->x) = _S1678;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1683;
    (&_S1683)->primal_0 = quat_6;
    (&_S1683)->differential_0 = _S1681;
    s_bwd_normalize_impl_0(&_S1683, _S1682);
    float _S1684 = - (_S1507 / _S1358);
    DiffPair_float_0 _S1685;
    (&_S1685)->primal_0 = _S1355;
    (&_S1685)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1685, _S1684);
    float _S1686 = - _S1685.differential_0;
    float3  _S1687 = _S1574 + _S1600 + _S1626 + _S1645 + make_float3 (0.0f, 0.0f, _S1549 + v_depth_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1688;
    (&_S1688)->primal_0 = R_6;
    (&_S1688)->differential_0 = _S1540;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1689;
    (&_S1689)->primal_0 = mean_7;
    (&_S1689)->differential_0 = _S1483;
    s_bwd_prop_mul_3(&_S1688, &_S1689, _S1687);
    Matrix<float, 3, 3>  _S1690 = _S1648 + _S1649.differential_0 + _S1688.differential_0;
    float _S1691 = _S1686 + _S1490.differential_0;
    float3  _S1692 = _S1680.differential_0 + _S1489.differential_0;
    *v_mean_0 = *v_mean_0 + (_S1689.differential_0 + _S1488.differential_0);
    *v_quat_0 = *v_quat_0 + _S1683.differential_0;
    *v_scale_0 = *v_scale_0 + _S1692;
    *v_in_opacity_0 = *v_in_opacity_0 + _S1691;
    *v_R_0 = *v_R_0 + _S1690;
    *v_t_0 = *v_t_0 + _S1687;
    return;
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_0)
{
    DiffPair_float_0 _S1693 = { s_primal_ctx_sqrt_0(dpdpx_0->primal_0), 0.5f / s_primal_ctx_sqrt_0((F32_max((1.00000001168609742e-07f), (dpdpx_0->primal_0)))) * dpdpx_0->differential_0 };
    return _S1693;
}

inline __device__ DiffPair_float_0 s_primal_ctx_s_fwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_1)
{
    float _S1694 = *&((&dpdpx_1->differential_0)->x) * *&((&dpdpx_1->primal_0)->x);
    float _S1695 = *&((&dpdpx_1->differential_0)->y) * *&((&dpdpx_1->primal_0)->y);
    float s_diff_len_1 = _S1694 + _S1694 + (_S1695 + _S1695);
    DiffPair_float_0 _S1696;
    (&_S1696)->primal_0 = *&((&dpdpx_1->primal_0)->x) * *&((&dpdpx_1->primal_0)->x) + *&((&dpdpx_1->primal_0)->y) * *&((&dpdpx_1->primal_0)->y);
    (&_S1696)->differential_0 = s_diff_len_1;
    DiffPair_float_0 _S1697 = s_primal_ctx_d_sqrt_0(&_S1696);
    DiffPair_float_0 _S1698 = { _S1697.primal_0, _S1697.differential_0 };
    return _S1698;
}

inline __device__ float s_primal_ctx_atan2_0(float _S1699, float _S1700)
{
    return (F32_atan2((_S1699), (_S1700)));
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_2)
{
    float _S1701 = dpdpx_2->primal_0 * dpdpx_2->primal_0 + dpdpy_0->primal_0 * dpdpy_0->primal_0;
    DiffPair_float_0 _S1702 = { s_primal_ctx_atan2_0(dpdpy_0->primal_0, dpdpx_2->primal_0), - dpdpy_0->primal_0 / _S1701 * dpdpx_2->differential_0 + dpdpx_2->primal_0 / _S1701 * dpdpy_0->differential_0 };
    return _S1702;
}

struct DiffPair_0
{
    DiffPair_float_0 primal_0;
    DiffPair_float_0 differential_0;
};

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S1703, DiffPair_float_0 * _S1704, float _S1705)
{
    _d_atan2_0(_S1703, _S1704, _S1705);
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_0 * dpdpy_1, DiffPair_0 * dpdpx_3, DiffPair_float_0 * _s_dOut_4)
{
    float _S1706 = - (*dpdpy_1).primal_0.primal_0;
    float _S1707 = (*dpdpx_3).primal_0.primal_0;
    float _S1708 = _S1707 * _S1707 + (*dpdpy_1).primal_0.primal_0 * (*dpdpy_1).primal_0.primal_0;
    float _S1709 = _S1706 / _S1708;
    float _S1710 = _S1708 * _S1708;
    float _S1711 = (*dpdpx_3).primal_0.primal_0 / _S1708;
    DiffPair_float_0 _S1712;
    (&_S1712)->primal_0 = (*dpdpy_1).primal_0.primal_0;
    (&_S1712)->differential_0 = 0.0f;
    DiffPair_float_0 _S1713;
    (&_S1713)->primal_0 = (*dpdpx_3).primal_0.primal_0;
    (&_S1713)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1712, &_S1713, _s_dOut_4->primal_0);
    float _S1714 = _S1711 * _s_dOut_4->differential_0;
    float _S1715 = (*dpdpy_1).primal_0.differential_0 * _s_dOut_4->differential_0 / _S1710;
    float _S1716 = (*dpdpx_3).primal_0.primal_0 * - _S1715;
    float _S1717 = (*dpdpy_1).primal_0.primal_0 * _S1716;
    float _S1718 = (*dpdpx_3).primal_0.primal_0 * _S1716;
    float _S1719 = (*dpdpx_3).primal_0.differential_0 * _s_dOut_4->differential_0 / _S1710;
    float _S1720 = _S1706 * - _S1719;
    float _S1721 = (*dpdpy_1).primal_0.primal_0 * _S1720;
    float _S1722 = (*dpdpx_3).primal_0.primal_0 * _S1720;
    float _S1723 = - (_S1708 * _S1719);
    DiffPair_float_0 _S1724 = { _S1713.differential_0 + _S1718 + _S1718 + _S1708 * _S1715 + _S1722 + _S1722, _S1709 * _s_dOut_4->differential_0 };
    dpdpx_3->primal_0 = (*dpdpx_3).primal_0;
    dpdpx_3->differential_0 = _S1724;
    DiffPair_float_0 _S1725 = { _S1712.differential_0 + _S1717 + _S1717 + _S1721 + _S1721 + _S1723, _S1714 };
    dpdpy_1->primal_0 = (*dpdpy_1).primal_0;
    dpdpy_1->differential_0 = _S1725;
    return;
}

struct DiffPair_1
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 differential_0;
};

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_0 * dpdpx_4, DiffPair_float_0 * _s_dOut_5)
{
    float _S1726 = (F32_max((1.00000001168609742e-07f), ((*dpdpx_4).primal_0.primal_0)));
    float _S1727 = s_primal_ctx_sqrt_0(_S1726);
    float _S1728 = 0.5f / _S1727 * _s_dOut_5->differential_0;
    float _S1729 = 0.5f * - ((*dpdpx_4).primal_0.differential_0 * _s_dOut_5->differential_0 / (_S1727 * _S1727));
    DiffPair_float_0 _S1730;
    (&_S1730)->primal_0 = _S1726;
    (&_S1730)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1730, _S1729);
    DiffPair_float_0 _S1731;
    (&_S1731)->primal_0 = 1.00000001168609742e-07f;
    (&_S1731)->differential_0 = 0.0f;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S1732)->differential_0 = 0.0f;
    _d_max_0(&_S1731, &_S1732, _S1730.differential_0);
    DiffPair_float_0 _S1733;
    (&_S1733)->primal_0 = (*dpdpx_4).primal_0.primal_0;
    (&_S1733)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1733, _s_dOut_5->primal_0);
    DiffPair_float_0 _S1734 = { _S1732.differential_0 + _S1733.differential_0, _S1728 };
    dpdpx_4->primal_0 = (*dpdpx_4).primal_0;
    dpdpx_4->differential_0 = _S1734;
    return;
}

inline __device__ void s_bwd_prop_s_fwd_length_impl_0(DiffPair_1 * dpdpx_5, DiffPair_float_0 * _s_dOut_6)
{
    float _S1735 = (*dpdpx_5).primal_0.primal_0.x;
    float _S1736 = (*dpdpx_5).primal_0.differential_0.x * (*dpdpx_5).primal_0.primal_0.x;
    float _S1737 = (*dpdpx_5).primal_0.primal_0.y;
    float _S1738 = (*dpdpx_5).primal_0.differential_0.y * (*dpdpx_5).primal_0.primal_0.y;
    DiffPair_float_0 _S1739 = { _S1735 * _S1735 + _S1737 * _S1737, _S1736 + _S1736 + (_S1738 + _S1738) };
    DiffPair_float_0 _S1740 = { 0.0f, 0.0f };
    DiffPair_0 _S1741;
    (&_S1741)->primal_0 = _S1739;
    (&_S1741)->differential_0 = _S1740;
    DiffPair_float_0 _S1742;
    (&_S1742)->primal_0 = _s_dOut_6->primal_0;
    (&_S1742)->differential_0 = _s_dOut_6->differential_0;
    s_bwd_prop_d_sqrt_0(&_S1741, &_S1742);
    float _S1743 = _S1741.differential_0.differential_0;
    float _S1744 = _S1743 + _S1743;
    float _S1745 = (*dpdpx_5).primal_0.primal_0.y * _S1744;
    float _S1746 = (*dpdpx_5).primal_0.primal_0.y * _S1741.differential_0.primal_0;
    float _S1747 = (*dpdpx_5).primal_0.differential_0.y * _S1744 + _S1746 + _S1746;
    float _S1748 = (*dpdpx_5).primal_0.primal_0.x * _S1744;
    float _S1749 = (*dpdpx_5).primal_0.primal_0.x * _S1741.differential_0.primal_0;
    float _S1750 = (*dpdpx_5).primal_0.differential_0.x * _S1744 + _S1749 + _S1749;
    float2  _S1751 = make_float2 (0.0f);
    float2  _S1752 = _S1751;
    *&((&_S1752)->y) = _S1747;
    *&((&_S1752)->x) = _S1750;
    float2  _S1753 = _S1751;
    *&((&_S1753)->y) = _S1745;
    *&((&_S1753)->x) = _S1748;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1754 = { _S1752, _S1753 };
    dpdpx_5->primal_0 = (*dpdpx_5).primal_0;
    dpdpx_5->differential_0 = _S1754;
    return;
}

inline __device__ void s_bwd_prop_length_impl_2(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_17, float _s_dOut_7)
{
    float _S1755 = (*dpx_17).primal_0.x;
    float _S1756 = (*dpx_17).primal_0.y;
    DiffPair_float_0 _S1757;
    (&_S1757)->primal_0 = _S1755 * _S1755 + _S1756 * _S1756;
    (&_S1757)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1757, _s_dOut_7);
    float _S1758 = (*dpx_17).primal_0.y * _S1757.differential_0;
    float _S1759 = _S1758 + _S1758;
    float _S1760 = (*dpx_17).primal_0.x * _S1757.differential_0;
    float _S1761 = _S1760 + _S1760;
    float2  _S1762 = make_float2 (0.0f);
    *&((&_S1762)->y) = _S1759;
    *&((&_S1762)->x) = _S1761;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S1762;
    return;
}

inline __device__ void s_bwd_length_impl_2(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1763, float _S1764)
{
    s_bwd_prop_length_impl_2(_S1763, _S1764);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_8, float4  quat_7, float3  scale_7, float in_opacity_7, Matrix<float, 3, 3>  R_7, float3  t_7, float fx_10, float fy_10, float cx_8, float cy_8, FixedArray<float, 10>  dist_coeffs_10, uint image_width_7, uint image_height_7, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_7, mean_8) + t_7;
    float _S1765 = - in_opacity_7;
    float _S1766 = 1.0f + s_primal_ctx_exp_0(_S1765);
    float _S1767 = 1.0f / _S1766;
    float _S1768 = _S1766 * _S1766;
    float4  _S1769 = normalize_0(quat_7);
    float3  _S1770 = s_primal_ctx_exp_1(scale_7);
    float _S1771 = _S1769.y;
    float x2_7 = _S1771 * _S1771;
    float y2_7 = _S1769.z * _S1769.z;
    float z2_7 = _S1769.w * _S1769.w;
    float xy_7 = _S1769.y * _S1769.z;
    float xz_7 = _S1769.y * _S1769.w;
    float yz_7 = _S1769.z * _S1769.w;
    float wx_7 = _S1769.x * _S1769.y;
    float wy_7 = _S1769.x * _S1769.z;
    float wz_7 = _S1769.x * _S1769.w;
    Matrix<float, 3, 3>  _S1772 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_7), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1770.x, 0.0f, 0.0f, 0.0f, _S1770.y, 0.0f, 0.0f, 0.0f, _S1770.z);
    Matrix<float, 3, 3>  _S1773 = s_primal_ctx_mul_1(_S1772, S_1);
    Matrix<float, 3, 3>  _S1774 = transpose_3(_S1773);
    Matrix<float, 3, 3>  _S1775 = s_primal_ctx_mul_1(_S1773, _S1774);
    Matrix<float, 3, 3>  _S1776 = s_primal_ctx_mul_1(R_7, _S1775);
    Matrix<float, 3, 3>  _S1777 = transpose_3(R_7);
    Matrix<float, 3, 3>  _S1778 = s_primal_ctx_mul_1(_S1776, _S1777);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S1779 = float2 {mean_c_6.x, mean_c_6.y};
    float2  _S1780 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1781 = { _S1779, _S1780 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1782;
    (&_S1782)->primal_0 = _S1779;
    (&_S1782)->differential_0 = _S1780;
    DiffPair_float_0 _S1783 = s_primal_ctx_s_fwd_length_impl_0(&_S1782);
    float _S1784 = mean_c_6.z;
    DiffPair_float_0 _S1785 = { _S1783.primal_0, _S1783.differential_0 };
    DiffPair_float_0 _S1786 = { _S1784, 0.0f };
    DiffPair_float_0 _S1787;
    (&_S1787)->primal_0 = _S1783.primal_0;
    (&_S1787)->differential_0 = _S1783.differential_0;
    DiffPair_float_0 _S1788;
    (&_S1788)->primal_0 = _S1784;
    (&_S1788)->differential_0 = 0.0f;
    DiffPair_float_0 _S1789 = s_primal_ctx_d_atan2_0(&_S1787, &_S1788);
    bool _S1790 = (_S1789.primal_0) < 0.00100000004749745f;
    float k_6;
    float s_diff_k_2;
    float _S1791;
    float _S1792;
    float _S1793;
    float _S1794;
    float _S1795;
    float _S1796;
    float _S1797;
    float _S1798;
    if(_S1790)
    {
        float _S1799 = _S1789.differential_0 * _S1789.primal_0;
        float _S1800 = 1.0f - _S1789.primal_0 * _S1789.primal_0 / 3.0f;
        float _S1801 = 0.0f - (_S1799 + _S1799) * 0.3333333432674408f;
        float _S1802 = _S1784 * _S1784;
        float _S1803 = _S1801 * _S1784;
        float _S1804 = _S1803 / _S1802;
        float _S1805 = _S1802 * _S1802;
        k_6 = _S1800 / _S1784;
        s_diff_k_2 = _S1804;
        _S1791 = _S1805;
        _S1792 = _S1803;
        _S1793 = _S1802;
        _S1794 = _S1800;
        _S1795 = _S1801;
        _S1796 = 0.0f;
        _S1797 = 0.0f;
        _S1798 = 0.0f;
    }
    else
    {
        float _S1806 = _S1783.primal_0 * _S1783.primal_0;
        float _S1807 = _S1789.differential_0 * _S1783.primal_0 - _S1789.primal_0 * _S1783.differential_0;
        float _S1808 = _S1807 / _S1806;
        float _S1809 = _S1806 * _S1806;
        k_6 = _S1789.primal_0 / _S1783.primal_0;
        s_diff_k_2 = _S1808;
        _S1791 = 0.0f;
        _S1792 = 0.0f;
        _S1793 = 0.0f;
        _S1794 = 0.0f;
        _S1795 = 0.0f;
        _S1796 = _S1809;
        _S1797 = _S1807;
        _S1798 = _S1806;
    }
    float2  _S1810 = make_float2 (k_6);
    float2  _S1811 = make_float2 (s_diff_k_2);
    float2  _S1812 = _S1779 * make_float2 (k_6);
    float2  _S1813 = _S1780 * make_float2 (k_6) + make_float2 (s_diff_k_2) * _S1779;
    float u_76 = _S1812.x;
    float s_diff_u_21 = _S1813.x;
    float v_76 = _S1812.y;
    float s_diff_v_21 = _S1813.y;
    float _S1814 = s_diff_u_21 * u_76;
    float _S1815 = s_diff_v_21 * v_76;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float s_diff_r2_75 = _S1814 + _S1814 + (_S1815 + _S1815);
    float _S1816 = s_diff_r2_75 * dist_coeffs_10[int(3)];
    float _S1817 = dist_coeffs_10[int(2)] + r2_76 * dist_coeffs_10[int(3)];
    float _S1818 = s_diff_r2_75 * _S1817 + _S1816 * r2_76;
    float _S1819 = dist_coeffs_10[int(1)] + r2_76 * _S1817;
    float _S1820 = s_diff_r2_75 * _S1819 + _S1818 * r2_76;
    float _S1821 = dist_coeffs_10[int(0)] + r2_76 * _S1819;
    float _S1822 = s_diff_r2_75 * _S1821 + _S1820 * r2_76;
    float2  _S1823 = make_float2 (_S1822);
    float radial_30 = 1.0f + r2_76 * _S1821;
    float2  _S1824 = make_float2 (radial_30);
    float _S1825 = 2.0f * dist_coeffs_10[int(4)];
    float _S1826 = _S1825 * u_76;
    float _S1827 = s_diff_u_21 * _S1825;
    float _S1828 = 2.0f * u_76;
    float _S1829 = s_diff_u_21 * 2.0f;
    float _S1830 = 2.0f * dist_coeffs_10[int(5)];
    float _S1831 = _S1830 * u_76;
    float _S1832 = s_diff_u_21 * _S1830;
    float _S1833 = 2.0f * v_76;
    float _S1834 = s_diff_v_21 * 2.0f;
    float2  _S1835 = _S1813 * make_float2 (radial_30) + make_float2 (_S1822) * _S1812 + make_float2 (_S1827 * v_76 + s_diff_v_21 * _S1826 + (s_diff_r2_75 + (_S1829 * u_76 + s_diff_u_21 * _S1828)) * dist_coeffs_10[int(5)] + s_diff_r2_75 * dist_coeffs_10[int(6)], _S1832 * v_76 + s_diff_v_21 * _S1831 + (s_diff_r2_75 + (_S1834 * v_76 + s_diff_v_21 * _S1833)) * dist_coeffs_10[int(4)] + s_diff_r2_75 * dist_coeffs_10[int(7)]);
    float2  _S1836 = _S1835 + make_float2 (_S1835.x * dist_coeffs_10[int(8)] + _S1835.y * dist_coeffs_10[int(9)], 0.0f);
    float _S1837 = _S1836.x * fx_10;
    float _S1838 = _S1836.y * fy_10;
    Matrix<float, 2, 3>  _S1839 = J_9;
    *&(((&_S1839)->rows + (int(0)))->x) = _S1837;
    *&(((&_S1839)->rows + (int(1)))->x) = _S1838;
    float2  _S1840 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1841 = { _S1779, _S1840 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1842;
    (&_S1842)->primal_0 = _S1779;
    (&_S1842)->differential_0 = _S1840;
    DiffPair_float_0 _S1843 = s_primal_ctx_s_fwd_length_impl_0(&_S1842);
    DiffPair_float_0 _S1844 = { _S1843.primal_0, _S1843.differential_0 };
    DiffPair_float_0 _S1845;
    (&_S1845)->primal_0 = _S1843.primal_0;
    (&_S1845)->differential_0 = _S1843.differential_0;
    DiffPair_float_0 _S1846;
    (&_S1846)->primal_0 = _S1784;
    (&_S1846)->differential_0 = 0.0f;
    DiffPair_float_0 _S1847 = s_primal_ctx_d_atan2_0(&_S1845, &_S1846);
    bool _S1848 = (_S1847.primal_0) < 0.00100000004749745f;
    float _S1849;
    float _S1850;
    float _S1851;
    float _S1852;
    float _S1853;
    float _S1854;
    float _S1855;
    float _S1856;
    if(_S1848)
    {
        float _S1857 = _S1847.differential_0 * _S1847.primal_0;
        float _S1858 = 1.0f - _S1847.primal_0 * _S1847.primal_0 / 3.0f;
        float _S1859 = 0.0f - (_S1857 + _S1857) * 0.3333333432674408f;
        float _S1860 = _S1784 * _S1784;
        float _S1861 = _S1859 * _S1784;
        float _S1862 = _S1861 / _S1860;
        float _S1863 = _S1860 * _S1860;
        k_6 = _S1858 / _S1784;
        s_diff_k_2 = _S1862;
        _S1849 = _S1863;
        _S1850 = _S1861;
        _S1851 = _S1860;
        _S1852 = _S1858;
        _S1853 = _S1859;
        _S1854 = 0.0f;
        _S1855 = 0.0f;
        _S1856 = 0.0f;
    }
    else
    {
        float _S1864 = _S1843.primal_0 * _S1843.primal_0;
        float _S1865 = _S1847.differential_0 * _S1843.primal_0 - _S1847.primal_0 * _S1843.differential_0;
        float _S1866 = _S1865 / _S1864;
        float _S1867 = _S1864 * _S1864;
        k_6 = _S1847.primal_0 / _S1843.primal_0;
        s_diff_k_2 = _S1866;
        _S1849 = 0.0f;
        _S1850 = 0.0f;
        _S1851 = 0.0f;
        _S1852 = 0.0f;
        _S1853 = 0.0f;
        _S1854 = _S1867;
        _S1855 = _S1865;
        _S1856 = _S1864;
    }
    float2  _S1868 = make_float2 (k_6);
    float2  _S1869 = make_float2 (s_diff_k_2);
    float2  _S1870 = _S1779 * make_float2 (k_6);
    float2  _S1871 = _S1840 * make_float2 (k_6) + make_float2 (s_diff_k_2) * _S1779;
    float u_77 = _S1870.x;
    float s_diff_u_22 = _S1871.x;
    float v_77 = _S1870.y;
    float s_diff_v_22 = _S1871.y;
    float _S1872 = s_diff_u_22 * u_77;
    float _S1873 = s_diff_v_22 * v_77;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float s_diff_r2_76 = _S1872 + _S1872 + (_S1873 + _S1873);
    float _S1874 = s_diff_r2_76 * dist_coeffs_10[int(3)];
    float _S1875 = dist_coeffs_10[int(2)] + r2_77 * dist_coeffs_10[int(3)];
    float _S1876 = s_diff_r2_76 * _S1875 + _S1874 * r2_77;
    float _S1877 = dist_coeffs_10[int(1)] + r2_77 * _S1875;
    float _S1878 = s_diff_r2_76 * _S1877 + _S1876 * r2_77;
    float _S1879 = dist_coeffs_10[int(0)] + r2_77 * _S1877;
    float _S1880 = s_diff_r2_76 * _S1879 + _S1878 * r2_77;
    float2  _S1881 = make_float2 (_S1880);
    float radial_31 = 1.0f + r2_77 * _S1879;
    float2  _S1882 = make_float2 (radial_31);
    float _S1883 = _S1825 * u_77;
    float _S1884 = s_diff_u_22 * _S1825;
    float _S1885 = 2.0f * u_77;
    float _S1886 = s_diff_u_22 * 2.0f;
    float _S1887 = _S1830 * u_77;
    float _S1888 = s_diff_u_22 * _S1830;
    float _S1889 = 2.0f * v_77;
    float _S1890 = s_diff_v_22 * 2.0f;
    float2  _S1891 = _S1871 * make_float2 (radial_31) + make_float2 (_S1880) * _S1870 + make_float2 (_S1884 * v_77 + s_diff_v_22 * _S1883 + (s_diff_r2_76 + (_S1886 * u_77 + s_diff_u_22 * _S1885)) * dist_coeffs_10[int(5)] + s_diff_r2_76 * dist_coeffs_10[int(6)], _S1888 * v_77 + s_diff_v_22 * _S1887 + (s_diff_r2_76 + (_S1890 * v_77 + s_diff_v_22 * _S1889)) * dist_coeffs_10[int(4)] + s_diff_r2_76 * dist_coeffs_10[int(7)]);
    float2  _S1892 = _S1891 + make_float2 (_S1891.x * dist_coeffs_10[int(8)] + _S1891.y * dist_coeffs_10[int(9)], 0.0f);
    float _S1893 = _S1892.y * fy_10;
    *&(((&_S1839)->rows + (int(0)))->y) = _S1892.x * fx_10;
    *&(((&_S1839)->rows + (int(1)))->y) = _S1893;
    float2  _S1894 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1895 = { _S1779, _S1894 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1896;
    (&_S1896)->primal_0 = _S1779;
    (&_S1896)->differential_0 = _S1894;
    DiffPair_float_0 _S1897 = s_primal_ctx_s_fwd_length_impl_0(&_S1896);
    DiffPair_float_0 _S1898 = { _S1897.primal_0, _S1897.differential_0 };
    DiffPair_float_0 _S1899 = { _S1784, 1.0f };
    DiffPair_float_0 _S1900;
    (&_S1900)->primal_0 = _S1897.primal_0;
    (&_S1900)->differential_0 = _S1897.differential_0;
    DiffPair_float_0 _S1901;
    (&_S1901)->primal_0 = _S1784;
    (&_S1901)->differential_0 = 1.0f;
    DiffPair_float_0 _S1902 = s_primal_ctx_d_atan2_0(&_S1900, &_S1901);
    bool _S1903 = (_S1902.primal_0) < 0.00100000004749745f;
    float _S1904;
    float _S1905;
    float _S1906;
    float _S1907;
    float _S1908;
    float _S1909;
    float _S1910;
    float _S1911;
    if(_S1903)
    {
        float _S1912 = _S1902.differential_0 * _S1902.primal_0;
        float _S1913 = 1.0f - _S1902.primal_0 * _S1902.primal_0 / 3.0f;
        float _S1914 = 0.0f - (_S1912 + _S1912) * 0.3333333432674408f;
        float _S1915 = _S1784 * _S1784;
        float _S1916 = _S1914 * _S1784 - _S1913;
        float _S1917 = _S1916 / _S1915;
        float _S1918 = _S1915 * _S1915;
        k_6 = _S1913 / _S1784;
        s_diff_k_2 = _S1917;
        _S1904 = _S1918;
        _S1905 = _S1916;
        _S1906 = _S1915;
        _S1907 = _S1913;
        _S1908 = _S1914;
        _S1909 = 0.0f;
        _S1910 = 0.0f;
        _S1911 = 0.0f;
    }
    else
    {
        float _S1919 = _S1897.primal_0 * _S1897.primal_0;
        float _S1920 = _S1902.differential_0 * _S1897.primal_0 - _S1902.primal_0 * _S1897.differential_0;
        float _S1921 = _S1920 / _S1919;
        float _S1922 = _S1919 * _S1919;
        k_6 = _S1902.primal_0 / _S1897.primal_0;
        s_diff_k_2 = _S1921;
        _S1904 = 0.0f;
        _S1905 = 0.0f;
        _S1906 = 0.0f;
        _S1907 = 0.0f;
        _S1908 = 0.0f;
        _S1909 = _S1922;
        _S1910 = _S1920;
        _S1911 = _S1919;
    }
    float2  _S1923 = make_float2 (k_6);
    float2  _S1924 = make_float2 (s_diff_k_2);
    float2  _S1925 = _S1779 * make_float2 (k_6);
    float2  _S1926 = make_float2 (s_diff_k_2) * _S1779;
    float u_78 = _S1925.x;
    float s_diff_u_23 = _S1926.x;
    float v_78 = _S1925.y;
    float s_diff_v_23 = _S1926.y;
    float _S1927 = s_diff_u_23 * u_78;
    float _S1928 = s_diff_v_23 * v_78;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float s_diff_r2_77 = _S1927 + _S1927 + (_S1928 + _S1928);
    float _S1929 = s_diff_r2_77 * dist_coeffs_10[int(3)];
    float _S1930 = dist_coeffs_10[int(2)] + r2_78 * dist_coeffs_10[int(3)];
    float _S1931 = s_diff_r2_77 * _S1930 + _S1929 * r2_78;
    float _S1932 = dist_coeffs_10[int(1)] + r2_78 * _S1930;
    float _S1933 = s_diff_r2_77 * _S1932 + _S1931 * r2_78;
    float _S1934 = dist_coeffs_10[int(0)] + r2_78 * _S1932;
    float _S1935 = s_diff_r2_77 * _S1934 + _S1933 * r2_78;
    float2  _S1936 = make_float2 (_S1935);
    float radial_32 = 1.0f + r2_78 * _S1934;
    float2  _S1937 = make_float2 (radial_32);
    float _S1938 = _S1825 * u_78;
    float _S1939 = s_diff_u_23 * _S1825;
    float _S1940 = 2.0f * u_78;
    float _S1941 = s_diff_u_23 * 2.0f;
    float _S1942 = _S1830 * u_78;
    float _S1943 = s_diff_u_23 * _S1830;
    float _S1944 = 2.0f * v_78;
    float _S1945 = s_diff_v_23 * 2.0f;
    float2  _S1946 = _S1926 * make_float2 (radial_32) + make_float2 (_S1935) * _S1925 + make_float2 (_S1939 * v_78 + s_diff_v_23 * _S1938 + (s_diff_r2_77 + (_S1941 * u_78 + s_diff_u_23 * _S1940)) * dist_coeffs_10[int(5)] + s_diff_r2_77 * dist_coeffs_10[int(6)], _S1943 * v_78 + s_diff_v_23 * _S1942 + (s_diff_r2_77 + (_S1945 * v_78 + s_diff_v_23 * _S1944)) * dist_coeffs_10[int(4)] + s_diff_r2_77 * dist_coeffs_10[int(7)]);
    float2  _S1947 = _S1946 + make_float2 (_S1946.x * dist_coeffs_10[int(8)] + _S1946.y * dist_coeffs_10[int(9)], 0.0f);
    float _S1948 = _S1947.y * fy_10;
    *&(((&_S1839)->rows + (int(0)))->z) = _S1947.x * fx_10;
    *&(((&_S1839)->rows + (int(1)))->z) = _S1948;
    Matrix<float, 2, 3>  _S1949 = s_primal_ctx_mul_2(_S1839, _S1778);
    Matrix<float, 3, 2>  _S1950 = transpose_1(_S1839);
    Matrix<float, 2, 2>  _S1951 = s_primal_ctx_mul_3(_S1949, _S1950);
    float eps2d_7;
    if(antialiased_7)
    {
        eps2d_7 = 0.10000000149011612f;
    }
    else
    {
        eps2d_7 = 0.30000001192092896f;
    }
    float _S1952 = _S1951.rows[int(0)].y * _S1951.rows[int(1)].x;
    float det_orig_7 = _S1951.rows[int(0)].x * _S1951.rows[int(1)].y - _S1952;
    float _S1953 = _S1951.rows[int(0)].x + eps2d_7;
    Matrix<float, 2, 2>  _S1954 = _S1951;
    *&(((&_S1954)->rows + (int(0)))->x) = _S1953;
    float _S1955 = _S1951.rows[int(1)].y + eps2d_7;
    *&(((&_S1954)->rows + (int(1)))->y) = _S1955;
    Matrix<float, 2, 2>  _S1956 = _S1954;
    Matrix<float, 2, 2>  _S1957 = _S1954;
    float det_blur_7 = _S1953 * _S1955 - _S1952;
    float _S1958 = det_orig_7 / det_blur_7;
    float _S1959 = det_blur_7 * det_blur_7;
    float _S1960 = (F32_max((0.0f), (_S1958)));
    float _S1961 = s_primal_ctx_sqrt_0(_S1960);
    float invdet_9 = 1.0f / det_blur_7;
    float _S1962 = - _S1951.rows[int(0)].y;
    float _S1963 = - _S1951.rows[int(1)].x;
    if(antialiased_7)
    {
        k_6 = _S1767 * _S1961;
    }
    else
    {
        k_6 = _S1767;
    }
    float _S1964 = k_6 / 0.00392156885936856f;
    float _S1965 = 2.0f * s_primal_ctx_log_0(_S1964);
    float _S1966 = s_primal_ctx_sqrt_0(_S1965);
    float _S1967 = _S1956.rows[int(0)].x;
    float _S1968 = _S1957.rows[int(1)].y;
    float3  campos_2 = - s_primal_ctx_mul_0(_S1777, t_7);
    float3  _S1969 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1970;
    (&_S1970)->primal_0 = mean_8;
    (&_S1970)->differential_0 = _S1969;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1971;
    (&_S1971)->primal_0 = scale_7;
    (&_S1971)->differential_0 = _S1969;
    DiffPair_float_0 _S1972;
    (&_S1972)->primal_0 = in_opacity_7;
    (&_S1972)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1973;
    (&_S1973)->primal_0 = campos_2;
    (&_S1973)->differential_0 = _S1969;
    s_bwd_prop_view_radius_3dgs_0(&_S1970, &_S1971, &_S1972, &_S1973, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1974 = _S1970;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1975 = _S1971;
    DiffPair_float_0 _S1976 = _S1972;
    float2  _S1977 = make_float2 (0.0f);
    float2  _S1978 = _S1977;
    *&((&_S1978)->y) = v_conic_1.z;
    float2  _S1979 = _S1977;
    *&((&_S1979)->y) = v_conic_1.y;
    *&((&_S1979)->x) = v_conic_1.x;
    DiffPair_float_0 _S1980;
    (&_S1980)->primal_0 = _S1968;
    (&_S1980)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1980, 0.0f);
    DiffPair_float_0 _S1981;
    (&_S1981)->primal_0 = _S1967;
    (&_S1981)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1981, 0.0f);
    DiffPair_float_0 _S1982;
    (&_S1982)->primal_0 = 3.32999992370605469f;
    (&_S1982)->differential_0 = 0.0f;
    DiffPair_float_0 _S1983;
    (&_S1983)->primal_0 = _S1966;
    (&_S1983)->differential_0 = 0.0f;
    _d_min_0(&_S1982, &_S1983, 0.0f);
    DiffPair_float_0 _S1984;
    (&_S1984)->primal_0 = _S1965;
    (&_S1984)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1984, _S1983.differential_0);
    float _S1985 = 2.0f * _S1984.differential_0;
    DiffPair_float_0 _S1986;
    (&_S1986)->primal_0 = _S1964;
    (&_S1986)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1986, _S1985);
    float _S1987 = v_opacity_1 + 254.9999847412109375f * _S1986.differential_0;
    float2  _S1988 = make_float2 (_S1981.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S1989 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1990 = _S1989;
    _S1990[int(1)] = _S1978;
    _S1990[int(0)] = _S1979;
    Matrix<float, 2, 2>  _S1991 = _S1990;
    float2  _S1992 = make_float2 (0.0f, _S1980.differential_0);
    if(antialiased_7)
    {
        float _S1993 = _S1961 * _S1987;
        k_6 = _S1767 * _S1987;
        s_diff_k_2 = _S1993;
    }
    else
    {
        k_6 = 0.0f;
        s_diff_k_2 = _S1987;
    }
    float _S1994 = invdet_9 * _S1991.rows[int(1)].y;
    float _S1995 = - (invdet_9 * _S1991.rows[int(1)].x);
    float _S1996 = - (invdet_9 * _S1991.rows[int(0)].y);
    float _S1997 = invdet_9 * _S1991.rows[int(0)].x;
    float _S1998 = - ((_S1953 * _S1991.rows[int(1)].y + _S1963 * _S1991.rows[int(1)].x + _S1962 * _S1991.rows[int(0)].y + _S1955 * _S1991.rows[int(0)].x) / _S1959);
    DiffPair_float_0 _S1999;
    (&_S1999)->primal_0 = _S1960;
    (&_S1999)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1999, k_6);
    DiffPair_float_0 _S2000 = { 0.0f, 0.0f };
    DiffPair_float_0 _S2001;
    (&_S2001)->primal_0 = 0.0f;
    (&_S2001)->differential_0 = 0.0f;
    DiffPair_float_0 _S2002;
    (&_S2002)->primal_0 = _S1958;
    (&_S2002)->differential_0 = 0.0f;
    _d_max_0(&_S2001, &_S2002, _S1999.differential_0);
    float _S2003 = _S2002.differential_0 / _S1959;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2003;
    float _S2004 = det_orig_7 * - _S2003 + _S1998;
    float _S2005 = - _S2004;
    float _S2006 = _S1953 * _S2004;
    float _S2007 = _S1955 * _S2004;
    Matrix<float, 2, 2>  _S2008 = _S1989;
    _S2008[int(1)] = _S1992;
    _S2008[int(0)] = _S1988;
    _S1954 = _S2008;
    *&(((&_S1954)->rows + (int(1)))->y) = 0.0f;
    float _S2009 = _S2006 + _S2008.rows[int(1)].y + _S1997;
    *&(((&_S1954)->rows + (int(0)))->x) = 0.0f;
    float _S2010 = _S2007 + _S2008.rows[int(0)].x + _S1994;
    float _S2011 = _S2005 + - s_diff_det_orig_T_1;
    float _S2012 = _S1951.rows[int(0)].y * _S2011 + _S1995;
    float _S2013 = _S1951.rows[int(1)].x * _S2011 + _S1996;
    float _S2014 = _S1951.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2015 = _S2009 + _S1951.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2016 = _S1977;
    *&((&_S2016)->x) = _S2012;
    *&((&_S2016)->y) = _S2015;
    float _S2017 = _S2010 + _S2014;
    float2  _S2018 = _S1977;
    *&((&_S2018)->y) = _S2013;
    *&((&_S2018)->x) = _S2017;
    Matrix<float, 2, 2>  _S2019 = _S1989;
    _S2019[int(1)] = _S2016;
    _S2019[int(0)] = _S2018;
    Matrix<float, 2, 2>  _S2020 = _S1954 + _S2019;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2021;
    (&_S2021)->primal_0 = _S1949;
    (&_S2021)->differential_0 = J_9;
    Matrix<float, 3, 2>  _S2022 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2023;
    (&_S2023)->primal_0 = _S1950;
    (&_S2023)->differential_0 = _S2022;
    s_bwd_prop_mul_0(&_S2021, &_S2023, _S2020);
    Matrix<float, 2, 3>  _S2024 = transpose_2(_S2023.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2025;
    (&_S2025)->primal_0 = _S1839;
    (&_S2025)->differential_0 = J_9;
    Matrix<float, 3, 3>  _S2026 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2027;
    (&_S2027)->primal_0 = _S1778;
    (&_S2027)->differential_0 = _S2026;
    s_bwd_prop_mul_1(&_S2025, &_S2027, _S2021.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2028 = _S2027;
    Matrix<float, 2, 3>  _S2029 = _S2024 + _S2025.differential_0;
    float2  _S2030 = make_float2 (0.0f, _S2029.rows[int(1)].z) + make_float2 (_S2029.rows[int(0)].z, 0.0f);
    float _S2031 = fx_10 * _S2030.x;
    float2  _S2032 = make_float2 (_S2031, fy_10 * _S2030.y) + make_float2 (dist_coeffs_10[int(8)] * _S2031, dist_coeffs_10[int(9)] * _S2031);
    float2  _S2033 = _S1925 * _S2032;
    float2  _S2034 = _S1926 * _S2032;
    float _S2035 = dist_coeffs_10[int(4)] * _S2032.y;
    float _S2036 = dist_coeffs_10[int(5)] * _S2032.x;
    float _S2037 = _S2034.x + _S2034.y;
    float _S2038 = _S2033.x + _S2033.y;
    float _S2039 = r2_78 * _S2038;
    float _S2040 = s_diff_r2_77 * _S2038 + r2_78 * _S2037;
    float _S2041 = r2_78 * _S2039;
    float _S2042 = s_diff_r2_77 * _S2039 + r2_78 * _S2040;
    float _S2043 = dist_coeffs_10[int(7)] * _S2032.y + _S2035 + dist_coeffs_10[int(6)] * _S2032.x + _S2036 + _S1934 * _S2038 + _S1932 * _S2039 + _S1930 * _S2041 + dist_coeffs_10[int(3)] * (r2_78 * _S2041);
    float _S2044 = _S1933 * _S2038 + _S1934 * _S2037 + _S1931 * _S2039 + _S1932 * _S2040 + _S1929 * _S2041 + _S1930 * _S2042 + dist_coeffs_10[int(3)] * (s_diff_r2_77 * _S2041 + r2_78 * _S2042);
    float _S2045 = _S2043 + _S2043;
    float _S2046 = v_78 * _S2044;
    float _S2047 = u_78 * _S2044;
    float2  _S2048 = _S1936 * _S2032 + make_float2 (_S1830 * (s_diff_v_23 * _S2032.y) + _S1941 * _S2036 + 2.0f * (s_diff_u_23 * _S2036) + _S1825 * (s_diff_v_23 * _S2032.x) + s_diff_u_23 * _S2045 + _S2047 + _S2047, _S1945 * _S2035 + 2.0f * (s_diff_v_23 * _S2035) + _S1943 * _S2032.y + _S1939 * _S2032.x + s_diff_v_23 * _S2045 + _S2046 + _S2046);
    float2  _S2049 = _S1937 * _S2032 + make_float2 (_S1830 * (v_78 * _S2032.y) + _S1940 * _S2036 + 2.0f * (u_78 * _S2036) + _S1825 * (v_78 * _S2032.x) + u_78 * _S2045, _S1944 * _S2035 + 2.0f * (v_78 * _S2035) + _S1942 * _S2032.y + _S1938 * _S2032.x + v_78 * _S2045);
    float2  _S2050 = _S1779 * _S2049;
    float2  _S2051 = _S1779 * _S2048;
    float _S2052 = _S2051.x + _S2051.y;
    float _S2053 = _S2050.x + _S2050.y;
    float2  _S2054 = _S1924 * _S2049 + _S1923 * _S2048;
    if(_S1903)
    {
        float _S2055 = _S2053 / _S1904;
        float _S2056 = _S1906 * _S2055;
        float _S2057 = _S1784 * (_S1905 * - _S2055);
        float _S2058 = _S2052 / _S1906;
        float _S2059 = 0.3333333432674408f * - (_S1784 * _S2056);
        float _S2060 = _S2059 + _S2059;
        float _S2061 = _S1902.primal_0 * (0.3333333432674408f * - (- _S2056 + _S1784 * _S2058));
        float _S2062 = _S2057 + _S2057 + _S1908 * _S2056 + _S1907 * - _S2058;
        float _S2063 = _S1902.differential_0 * _S2060 + _S2061 + _S2061;
        k_6 = _S1902.primal_0 * _S2060;
        _S1904 = _S2063;
        _S1905 = _S2062;
        _S1906 = 0.0f;
        _S1907 = 0.0f;
    }
    else
    {
        float _S2064 = _S2053 / _S1909;
        float _S2065 = _S1911 * _S2064;
        float _S2066 = _S1897.primal_0 * (_S1910 * - _S2064);
        float _S2067 = - _S2065;
        float _S2068 = _S1902.primal_0 * _S2067;
        float _S2069 = _S2052 / _S1911;
        float _S2070 = _S2066 + _S2066 + _S1902.differential_0 * _S2065 + _S1902.primal_0 * - _S2069;
        float _S2071 = _S1897.differential_0 * _S2067 + _S1897.primal_0 * _S2069;
        k_6 = _S1897.primal_0 * _S2065;
        _S1904 = _S2071;
        _S1905 = 0.0f;
        _S1906 = _S2068;
        _S1907 = _S2070;
    }
    DiffPair_0 _S2072;
    (&_S2072)->primal_0 = _S1898;
    (&_S2072)->differential_0 = _S2000;
    DiffPair_0 _S2073;
    (&_S2073)->primal_0 = _S1899;
    (&_S2073)->differential_0 = _S2000;
    DiffPair_float_0 _S2074;
    (&_S2074)->primal_0 = _S1904;
    (&_S2074)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2072, &_S2073, &_S2074);
    float _S2075 = _S2073.differential_0.primal_0 + _S1905;
    float _S2076 = _S2072.differential_0.differential_0 + _S1906;
    float _S2077 = _S2072.differential_0.primal_0 + _S1907;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2078 = { _S1977, _S1977 };
    DiffPair_1 _S2079;
    (&_S2079)->primal_0 = _S1895;
    (&_S2079)->differential_0 = _S2078;
    DiffPair_float_0 _S2080;
    (&_S2080)->primal_0 = _S2077;
    (&_S2080)->differential_0 = _S2076;
    s_bwd_prop_s_fwd_length_impl_0(&_S2079, &_S2080);
    float2  _S2081 = _S2079.differential_0.primal_0 + _S2054;
    float3  _S2082 = make_float3 (_S2081.x, _S2081.y, _S2075);
    float2  _S2083 = make_float2 (0.0f, _S2029.rows[int(1)].y) + make_float2 (_S2029.rows[int(0)].y, 0.0f);
    float _S2084 = fx_10 * _S2083.x;
    float2  _S2085 = make_float2 (_S2084, fy_10 * _S2083.y) + make_float2 (dist_coeffs_10[int(8)] * _S2084, dist_coeffs_10[int(9)] * _S2084);
    float2  _S2086 = _S1870 * _S2085;
    float2  _S2087 = _S1871 * _S2085;
    float _S2088 = dist_coeffs_10[int(4)] * _S2085.y;
    float _S2089 = dist_coeffs_10[int(5)] * _S2085.x;
    float _S2090 = _S2087.x + _S2087.y;
    float _S2091 = _S2086.x + _S2086.y;
    float _S2092 = r2_77 * _S2091;
    float _S2093 = s_diff_r2_76 * _S2091 + r2_77 * _S2090;
    float _S2094 = r2_77 * _S2092;
    float _S2095 = s_diff_r2_76 * _S2092 + r2_77 * _S2093;
    float _S2096 = dist_coeffs_10[int(7)] * _S2085.y + _S2088 + dist_coeffs_10[int(6)] * _S2085.x + _S2089 + _S1879 * _S2091 + _S1877 * _S2092 + _S1875 * _S2094 + dist_coeffs_10[int(3)] * (r2_77 * _S2094);
    float _S2097 = _S1878 * _S2091 + _S1879 * _S2090 + _S1876 * _S2092 + _S1877 * _S2093 + _S1874 * _S2094 + _S1875 * _S2095 + dist_coeffs_10[int(3)] * (s_diff_r2_76 * _S2094 + r2_77 * _S2095);
    float _S2098 = _S2096 + _S2096;
    float _S2099 = v_77 * _S2097;
    float _S2100 = u_77 * _S2097;
    float2  _S2101 = _S1881 * _S2085 + make_float2 (_S1830 * (s_diff_v_22 * _S2085.y) + _S1886 * _S2089 + 2.0f * (s_diff_u_22 * _S2089) + _S1825 * (s_diff_v_22 * _S2085.x) + s_diff_u_22 * _S2098 + _S2100 + _S2100, _S1890 * _S2088 + 2.0f * (s_diff_v_22 * _S2088) + _S1888 * _S2085.y + _S1884 * _S2085.x + s_diff_v_22 * _S2098 + _S2099 + _S2099);
    float2  _S2102 = _S1882 * _S2085 + make_float2 (_S1830 * (v_77 * _S2085.y) + _S1885 * _S2089 + 2.0f * (u_77 * _S2089) + _S1825 * (v_77 * _S2085.x) + u_77 * _S2098, _S1889 * _S2088 + 2.0f * (v_77 * _S2088) + _S1887 * _S2085.y + _S1883 * _S2085.x + v_77 * _S2098);
    float2  _S2103 = _S1779 * _S2102;
    float2  _S2104 = _S1840 * _S2102;
    float2  _S2105 = _S1779 * _S2101;
    float _S2106 = _S2104.x + _S2104.y + _S2105.x + _S2105.y;
    float _S2107 = _S2103.x + _S2103.y;
    float2  _S2108 = _S1869 * _S2102 + _S1868 * _S2101;
    if(_S1848)
    {
        float _S2109 = _S2107 / _S1849;
        float _S2110 = _S1851 * _S2109;
        float _S2111 = _S1784 * (_S1850 * - _S2109);
        float _S2112 = _S2106 / _S1851;
        float _S2113 = 0.3333333432674408f * - (_S1784 * _S2110);
        float _S2114 = _S2113 + _S2113;
        float _S2115 = _S1847.primal_0 * (0.3333333432674408f * - (_S1784 * _S2112));
        float _S2116 = _S2111 + _S2111 + _S1853 * _S2110 + _S1852 * - _S2112;
        float _S2117 = _S1847.differential_0 * _S2114 + _S2115 + _S2115;
        k_6 = _S1847.primal_0 * _S2114;
        _S1849 = _S2117;
        _S1850 = _S2116;
        _S1851 = 0.0f;
        _S1852 = 0.0f;
    }
    else
    {
        float _S2118 = _S2107 / _S1854;
        float _S2119 = _S1856 * _S2118;
        float _S2120 = _S1843.primal_0 * (_S1855 * - _S2118);
        float _S2121 = - _S2119;
        float _S2122 = _S1847.primal_0 * _S2121;
        float _S2123 = _S2106 / _S1856;
        float _S2124 = _S2120 + _S2120 + _S1847.differential_0 * _S2119 + _S1847.primal_0 * - _S2123;
        float _S2125 = _S1843.differential_0 * _S2121 + _S1843.primal_0 * _S2123;
        k_6 = _S1843.primal_0 * _S2119;
        _S1849 = _S2125;
        _S1850 = 0.0f;
        _S1851 = _S2122;
        _S1852 = _S2124;
    }
    DiffPair_0 _S2126;
    (&_S2126)->primal_0 = _S1844;
    (&_S2126)->differential_0 = _S2000;
    DiffPair_0 _S2127;
    (&_S2127)->primal_0 = _S1786;
    (&_S2127)->differential_0 = _S2000;
    DiffPair_float_0 _S2128;
    (&_S2128)->primal_0 = _S1849;
    (&_S2128)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2126, &_S2127, &_S2128);
    float _S2129 = _S2127.differential_0.primal_0 + _S1850;
    float _S2130 = _S2126.differential_0.differential_0 + _S1851;
    float _S2131 = _S2126.differential_0.primal_0 + _S1852;
    DiffPair_1 _S2132;
    (&_S2132)->primal_0 = _S1841;
    (&_S2132)->differential_0 = _S2078;
    DiffPair_float_0 _S2133;
    (&_S2133)->primal_0 = _S2131;
    (&_S2133)->differential_0 = _S2130;
    s_bwd_prop_s_fwd_length_impl_0(&_S2132, &_S2133);
    float2  _S2134 = _S2132.differential_0.primal_0 + _S2108;
    float2  _S2135 = make_float2 (0.0f, _S2029.rows[int(1)].x) + make_float2 (_S2029.rows[int(0)].x, 0.0f);
    float _S2136 = fx_10 * _S2135.x;
    float2  _S2137 = make_float2 (_S2136, fy_10 * _S2135.y) + make_float2 (dist_coeffs_10[int(8)] * _S2136, dist_coeffs_10[int(9)] * _S2136);
    float2  _S2138 = _S1812 * _S2137;
    float2  _S2139 = _S1813 * _S2137;
    float _S2140 = dist_coeffs_10[int(4)] * _S2137.y;
    float _S2141 = dist_coeffs_10[int(5)] * _S2137.x;
    float _S2142 = _S2139.x + _S2139.y;
    float _S2143 = _S2138.x + _S2138.y;
    float _S2144 = r2_76 * _S2143;
    float _S2145 = s_diff_r2_75 * _S2143 + r2_76 * _S2142;
    float _S2146 = r2_76 * _S2144;
    float _S2147 = s_diff_r2_75 * _S2144 + r2_76 * _S2145;
    float _S2148 = dist_coeffs_10[int(7)] * _S2137.y + _S2140 + dist_coeffs_10[int(6)] * _S2137.x + _S2141 + _S1821 * _S2143 + _S1819 * _S2144 + _S1817 * _S2146 + dist_coeffs_10[int(3)] * (r2_76 * _S2146);
    float _S2149 = _S1820 * _S2143 + _S1821 * _S2142 + _S1818 * _S2144 + _S1819 * _S2145 + _S1816 * _S2146 + _S1817 * _S2147 + dist_coeffs_10[int(3)] * (s_diff_r2_75 * _S2146 + r2_76 * _S2147);
    float _S2150 = _S2148 + _S2148;
    float _S2151 = v_76 * _S2149;
    float _S2152 = u_76 * _S2149;
    float2  _S2153 = _S1823 * _S2137 + make_float2 (_S1830 * (s_diff_v_21 * _S2137.y) + _S1829 * _S2141 + 2.0f * (s_diff_u_21 * _S2141) + _S1825 * (s_diff_v_21 * _S2137.x) + s_diff_u_21 * _S2150 + _S2152 + _S2152, _S1834 * _S2140 + 2.0f * (s_diff_v_21 * _S2140) + _S1832 * _S2137.y + _S1827 * _S2137.x + s_diff_v_21 * _S2150 + _S2151 + _S2151);
    float2  _S2154 = _S1824 * _S2137 + make_float2 (_S1830 * (v_76 * _S2137.y) + _S1828 * _S2141 + 2.0f * (u_76 * _S2141) + _S1825 * (v_76 * _S2137.x) + u_76 * _S2150, _S1833 * _S2140 + 2.0f * (v_76 * _S2140) + _S1831 * _S2137.y + _S1826 * _S2137.x + v_76 * _S2150);
    float3  _S2155 = make_float3 (_S2134.x, _S2134.y, _S2129) + _S2082;
    float2  _S2156 = _S1779 * _S2154;
    float2  _S2157 = _S1780 * _S2154;
    float2  _S2158 = _S1779 * _S2153;
    float _S2159 = _S2157.x + _S2157.y + _S2158.x + _S2158.y;
    float _S2160 = _S2156.x + _S2156.y;
    float2  _S2161 = _S1811 * _S2154 + _S1810 * _S2153;
    if(_S1790)
    {
        float _S2162 = _S2160 / _S1791;
        float _S2163 = _S1793 * _S2162;
        float _S2164 = _S1784 * (_S1792 * - _S2162);
        float _S2165 = _S2159 / _S1793;
        float _S2166 = 0.3333333432674408f * - (_S1784 * _S2163);
        float _S2167 = _S2166 + _S2166;
        float _S2168 = _S1789.primal_0 * (0.3333333432674408f * - (_S1784 * _S2165));
        float _S2169 = _S2164 + _S2164 + _S1795 * _S2163 + _S1794 * - _S2165;
        float _S2170 = _S1789.differential_0 * _S2167 + _S2168 + _S2168;
        k_6 = _S1789.primal_0 * _S2167;
        _S1791 = _S2170;
        _S1792 = _S2169;
        _S1793 = 0.0f;
        _S1794 = 0.0f;
    }
    else
    {
        float _S2171 = _S2160 / _S1796;
        float _S2172 = _S1798 * _S2171;
        float _S2173 = _S1783.primal_0 * (_S1797 * - _S2171);
        float _S2174 = - _S2172;
        float _S2175 = _S1789.primal_0 * _S2174;
        float _S2176 = _S2159 / _S1798;
        float _S2177 = _S2173 + _S2173 + _S1789.differential_0 * _S2172 + _S1789.primal_0 * - _S2176;
        float _S2178 = _S1783.differential_0 * _S2174 + _S1783.primal_0 * _S2176;
        k_6 = _S1783.primal_0 * _S2172;
        _S1791 = _S2178;
        _S1792 = 0.0f;
        _S1793 = _S2175;
        _S1794 = _S2177;
    }
    DiffPair_0 _S2179;
    (&_S2179)->primal_0 = _S1785;
    (&_S2179)->differential_0 = _S2000;
    DiffPair_0 _S2180;
    (&_S2180)->primal_0 = _S1786;
    (&_S2180)->differential_0 = _S2000;
    DiffPair_float_0 _S2181;
    (&_S2181)->primal_0 = _S1791;
    (&_S2181)->differential_0 = k_6;
    s_bwd_prop_d_atan2_0(&_S2179, &_S2180, &_S2181);
    float _S2182 = _S2180.differential_0.primal_0 + _S1792;
    float _S2183 = _S2179.differential_0.differential_0 + _S1793;
    float _S2184 = _S2179.differential_0.primal_0 + _S1794;
    DiffPair_1 _S2185;
    (&_S2185)->primal_0 = _S1781;
    (&_S2185)->differential_0 = _S2078;
    DiffPair_float_0 _S2186;
    (&_S2186)->primal_0 = _S2184;
    (&_S2186)->differential_0 = _S2183;
    s_bwd_prop_s_fwd_length_impl_0(&_S2185, &_S2186);
    float2  _S2187 = _S2185.differential_0.primal_0 + _S2161;
    float3  _S2188 = make_float3 (_S2187.x, _S2187.y, _S2182);
    float _S2189 = length_0(_S1779);
    float _S2190 = s_primal_ctx_atan2_0(_S2189, _S1784);
    bool _S2191 = _S2190 < 0.00100000004749745f;
    if(_S2191)
    {
        float _S2192 = 1.0f - _S2190 * _S2190 / 3.0f;
        float _S2193 = _S1784 * _S1784;
        k_6 = _S2192 / _S1784;
        _S1791 = _S2193;
        _S1792 = _S2192;
        _S1793 = 0.0f;
    }
    else
    {
        float _S2194 = _S2189 * _S2189;
        k_6 = _S2190 / _S2189;
        _S1791 = 0.0f;
        _S1792 = 0.0f;
        _S1793 = _S2194;
    }
    float2  _S2195 = make_float2 (k_6);
    float2  _S2196 = _S1779 * make_float2 (k_6);
    float _S2197 = fx_10 * v_mean2d_1.x;
    float u_79 = _S2196.x;
    float v_79 = _S2196.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S2198 = dist_coeffs_10[int(2)] + r2_79 * dist_coeffs_10[int(3)];
    float _S2199 = dist_coeffs_10[int(1)] + r2_79 * _S2198;
    float _S2200 = dist_coeffs_10[int(0)] + r2_79 * _S2199;
    float2  _S2201 = make_float2 (_S2197, fy_10 * v_mean2d_1.y) + make_float2 (dist_coeffs_10[int(8)] * _S2197, dist_coeffs_10[int(9)] * _S2197);
    float2  _S2202 = _S2196 * _S2201;
    float _S2203 = dist_coeffs_10[int(4)] * _S2201.y;
    float _S2204 = dist_coeffs_10[int(5)] * _S2201.x;
    float _S2205 = _S2202.x + _S2202.y;
    float _S2206 = r2_79 * _S2205;
    float _S2207 = r2_79 * _S2206;
    float _S2208 = dist_coeffs_10[int(7)] * _S2201.y + _S2203 + dist_coeffs_10[int(6)] * _S2201.x + _S2204 + _S2200 * _S2205 + _S2199 * _S2206 + _S2198 * _S2207 + dist_coeffs_10[int(3)] * (r2_79 * _S2207);
    float _S2209 = v_79 * _S2208;
    float _S2210 = u_79 * _S2208;
    float2  _S2211 = make_float2 (1.0f + r2_79 * _S2200) * _S2201 + make_float2 (_S1830 * (v_79 * _S2201.y) + 2.0f * u_79 * _S2204 + 2.0f * (u_79 * _S2204) + _S1825 * (v_79 * _S2201.x) + _S2210 + _S2210, 2.0f * v_79 * _S2203 + 2.0f * (v_79 * _S2203) + _S1830 * u_79 * _S2201.y + _S1825 * u_79 * _S2201.x + _S2209 + _S2209);
    float2  _S2212 = _S1779 * _S2211;
    float2  _S2213 = _S2195 * _S2211;
    float _S2214 = _S2212.x + _S2212.y;
    if(_S2191)
    {
        float _S2215 = _S2214 / _S1791;
        float _S2216 = _S1792 * - _S2215;
        float _S2217 = _S2190 * (0.3333333432674408f * - (_S1784 * _S2215));
        k_6 = _S2217 + _S2217;
        _S1791 = _S2216;
        _S1792 = 0.0f;
    }
    else
    {
        float _S2218 = _S2214 / _S1793;
        float _S2219 = _S2190 * - _S2218;
        k_6 = _S2189 * _S2218;
        _S1791 = 0.0f;
        _S1792 = _S2219;
    }
    DiffPair_float_0 _S2220;
    (&_S2220)->primal_0 = _S2189;
    (&_S2220)->differential_0 = 0.0f;
    DiffPair_float_0 _S2221;
    (&_S2221)->primal_0 = _S1784;
    (&_S2221)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2220, &_S2221, k_6);
    float _S2222 = _S2221.differential_0 + _S1791;
    float _S2223 = _S2220.differential_0 + _S1792;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2224;
    (&_S2224)->primal_0 = _S1779;
    (&_S2224)->differential_0 = _S1977;
    s_bwd_length_impl_2(&_S2224, _S2223);
    float2  _S2225 = _S2224.differential_0 + _S2213;
    float3  _S2226 = make_float3 (_S2225.x, _S2225.y, _S2222);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2227;
    (&_S2227)->primal_0 = _S1776;
    (&_S2227)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2228;
    (&_S2228)->primal_0 = _S1777;
    (&_S2228)->differential_0 = _S2026;
    s_bwd_prop_mul_2(&_S2227, &_S2228, _S2028.differential_0);
    Matrix<float, 3, 3>  _S2229 = transpose_3(_S2228.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2230;
    (&_S2230)->primal_0 = R_7;
    (&_S2230)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2231;
    (&_S2231)->primal_0 = _S1775;
    (&_S2231)->differential_0 = _S2026;
    s_bwd_prop_mul_2(&_S2230, &_S2231, _S2227.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2232;
    (&_S2232)->primal_0 = _S1773;
    (&_S2232)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2233;
    (&_S2233)->primal_0 = _S1774;
    (&_S2233)->differential_0 = _S2026;
    s_bwd_prop_mul_2(&_S2232, &_S2233, _S2231.differential_0);
    Matrix<float, 3, 3>  _S2234 = _S2232.differential_0 + transpose_3(_S2233.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2235;
    (&_S2235)->primal_0 = _S1772;
    (&_S2235)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2236;
    (&_S2236)->primal_0 = S_1;
    (&_S2236)->differential_0 = _S2026;
    s_bwd_prop_mul_2(&_S2235, &_S2236, _S2234);
    Matrix<float, 3, 3>  _S2237 = transpose_3(_S2235.differential_0);
    float _S2238 = 2.0f * - _S2237.rows[int(2)].z;
    float _S2239 = 2.0f * _S2237.rows[int(2)].y;
    float _S2240 = 2.0f * _S2237.rows[int(2)].x;
    float _S2241 = 2.0f * _S2237.rows[int(1)].z;
    float _S2242 = 2.0f * - _S2237.rows[int(1)].y;
    float _S2243 = 2.0f * _S2237.rows[int(1)].x;
    float _S2244 = 2.0f * _S2237.rows[int(0)].z;
    float _S2245 = 2.0f * _S2237.rows[int(0)].y;
    float _S2246 = 2.0f * - _S2237.rows[int(0)].x;
    float _S2247 = - _S2243 + _S2245;
    float _S2248 = _S2240 + - _S2244;
    float _S2249 = - _S2239 + _S2241;
    float _S2250 = _S2239 + _S2241;
    float _S2251 = _S2240 + _S2244;
    float _S2252 = _S2243 + _S2245;
    float _S2253 = _S1769.w * (_S2242 + _S2246);
    float _S2254 = _S1769.z * (_S2238 + _S2246);
    float _S2255 = _S1769.y * (_S2238 + _S2242);
    float _S2256 = _S1769.x * _S2247 + _S1769.z * _S2250 + _S1769.y * _S2251 + _S2253 + _S2253;
    float _S2257 = _S1769.x * _S2248 + _S1769.w * _S2250 + _S1769.y * _S2252 + _S2254 + _S2254;
    float _S2258 = _S1769.x * _S2249 + _S1769.w * _S2251 + _S1769.z * _S2252 + _S2255 + _S2255;
    float _S2259 = _S1769.w * _S2247 + _S1769.z * _S2248 + _S1769.y * _S2249;
    float3  _S2260 = _S1969;
    *&((&_S2260)->z) = _S2236.differential_0.rows[int(2)].z;
    *&((&_S2260)->y) = _S2236.differential_0.rows[int(1)].y;
    *&((&_S2260)->x) = _S2236.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2261;
    (&_S2261)->primal_0 = scale_7;
    (&_S2261)->differential_0 = _S1969;
    s_bwd_prop_exp_1(&_S2261, _S2260);
    float4  _S2262 = make_float4 (0.0f);
    float4  _S2263 = _S2262;
    *&((&_S2263)->w) = _S2256;
    *&((&_S2263)->z) = _S2257;
    *&((&_S2263)->y) = _S2258;
    *&((&_S2263)->x) = _S2259;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2264;
    (&_S2264)->primal_0 = quat_7;
    (&_S2264)->differential_0 = _S2262;
    s_bwd_normalize_impl_0(&_S2264, _S2263);
    float _S2265 = - (s_diff_k_2 / _S1768);
    DiffPair_float_0 _S2266;
    (&_S2266)->primal_0 = _S1765;
    (&_S2266)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2266, _S2265);
    float _S2267 = - _S2266.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2268;
    (&_S2268)->primal_0 = mean_c_6;
    (&_S2268)->differential_0 = _S1969;
    s_bwd_length_impl_0(&_S2268, 0.0f);
    float3  _S2269 = _S2188 + _S2226 + _S2268.differential_0 + _S2155 + make_float3 (0.0f, 0.0f, v_depth_1);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2270;
    (&_S2270)->primal_0 = R_7;
    (&_S2270)->differential_0 = _S2026;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2271;
    (&_S2271)->primal_0 = mean_8;
    (&_S2271)->differential_0 = _S1969;
    s_bwd_prop_mul_3(&_S2270, &_S2271, _S2269);
    Matrix<float, 3, 3>  _S2272 = _S2229 + _S2230.differential_0 + _S2270.differential_0;
    float _S2273 = _S2267 + _S1976.differential_0;
    float3  _S2274 = _S2261.differential_0 + _S1975.differential_0;
    *v_mean_1 = *v_mean_1 + (_S2271.differential_0 + _S1974.differential_0);
    *v_quat_1 = *v_quat_1 + _S2264.differential_0;
    *v_scale_1 = *v_scale_1 + _S2274;
    *v_in_opacity_1 = *v_in_opacity_1 + _S2273;
    *v_R_1 = *v_R_1 + _S2272;
    *v_t_1 = *v_t_1 + _S2269;
    return;
}

inline __device__ float s_primal_ctx_sin_0(float _S2275)
{
    return (F32_sin((_S2275)));
}

inline __device__ float s_primal_ctx_cos_0(float _S2276)
{
    return (F32_cos((_S2276)));
}

inline __device__ DiffPair_float_0 s_primal_ctx_d_sin_0(DiffPair_float_0 * dpdpx_6)
{
    DiffPair_float_0 _S2277 = { s_primal_ctx_sin_0(dpdpx_6->primal_0), s_primal_ctx_cos_0(dpdpx_6->primal_0) * dpdpx_6->differential_0 };
    return _S2277;
}

inline __device__ void s_bwd_prop_cos_0(DiffPair_float_0 * _S2278, float _S2279)
{
    _d_cos_0(_S2278, _S2279);
    return;
}

inline __device__ void s_bwd_prop_sin_0(DiffPair_float_0 * _S2280, float _S2281)
{
    _d_sin_0(_S2280, _S2281);
    return;
}

inline __device__ void s_bwd_prop_d_sin_0(DiffPair_0 * dpdpx_7, DiffPair_float_0 * _s_dOut_8)
{
    float _S2282 = s_primal_ctx_cos_0((*dpdpx_7).primal_0.primal_0) * _s_dOut_8->differential_0;
    float _S2283 = (*dpdpx_7).primal_0.differential_0 * _s_dOut_8->differential_0;
    DiffPair_float_0 _S2284;
    (&_S2284)->primal_0 = (*dpdpx_7).primal_0.primal_0;
    (&_S2284)->differential_0 = 0.0f;
    s_bwd_prop_cos_0(&_S2284, _S2283);
    DiffPair_float_0 _S2285;
    (&_S2285)->primal_0 = (*dpdpx_7).primal_0.primal_0;
    (&_S2285)->differential_0 = 0.0f;
    s_bwd_prop_sin_0(&_S2285, _s_dOut_8->primal_0);
    DiffPair_float_0 _S2286 = { _S2284.differential_0 + _S2285.differential_0, _S2282 };
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S2286;
    return;
}

inline __device__ void projection_3dgs_equisolid_vjp(bool antialiased_8, float3  mean_9, float4  quat_8, float3  scale_8, float in_opacity_8, Matrix<float, 3, 3>  R_8, float3  t_8, float fx_11, float fy_11, float cx_9, float cy_9, FixedArray<float, 10>  dist_coeffs_11, uint image_width_8, uint image_height_8, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    DiffPair_float_0 _S2287 = { 0.0f, 0.0f };
    float3  mean_c_7 = s_primal_ctx_mul_0(R_8, mean_9) + t_8;
    float _S2288 = - in_opacity_8;
    float _S2289 = 1.0f + s_primal_ctx_exp_0(_S2288);
    float _S2290 = 1.0f / _S2289;
    float _S2291 = _S2289 * _S2289;
    float4  _S2292 = normalize_0(quat_8);
    float3  _S2293 = s_primal_ctx_exp_1(scale_8);
    float _S2294 = _S2292.y;
    float x2_8 = _S2294 * _S2294;
    float y2_8 = _S2292.z * _S2292.z;
    float z2_8 = _S2292.w * _S2292.w;
    float xy_8 = _S2292.y * _S2292.z;
    float xz_8 = _S2292.y * _S2292.w;
    float yz_8 = _S2292.z * _S2292.w;
    float wx_8 = _S2292.x * _S2292.y;
    float wy_8 = _S2292.x * _S2292.z;
    float wz_8 = _S2292.x * _S2292.w;
    Matrix<float, 3, 3>  _S2295 = transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_8), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_8), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2293.x, 0.0f, 0.0f, 0.0f, _S2293.y, 0.0f, 0.0f, 0.0f, _S2293.z);
    Matrix<float, 3, 3>  _S2296 = s_primal_ctx_mul_1(_S2295, S_2);
    Matrix<float, 3, 3>  _S2297 = transpose_3(_S2296);
    Matrix<float, 3, 3>  _S2298 = s_primal_ctx_mul_1(_S2296, _S2297);
    Matrix<float, 3, 3>  _S2299 = s_primal_ctx_mul_1(R_8, _S2298);
    Matrix<float, 3, 3>  _S2300 = transpose_3(R_8);
    Matrix<float, 3, 3>  _S2301 = s_primal_ctx_mul_1(_S2299, _S2300);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float2  _S2302 = float2 {mean_c_7.x, mean_c_7.y};
    float2  _S2303 = make_float2 (1.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2304 = { _S2302, _S2303 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2305;
    (&_S2305)->primal_0 = _S2302;
    (&_S2305)->differential_0 = _S2303;
    DiffPair_float_0 _S2306 = s_primal_ctx_s_fwd_length_impl_0(&_S2305);
    float _S2307 = mean_c_7.z;
    DiffPair_float_0 _S2308 = { _S2306.primal_0, _S2306.differential_0 };
    DiffPair_float_0 _S2309 = { _S2307, 0.0f };
    DiffPair_float_0 _S2310;
    (&_S2310)->primal_0 = _S2306.primal_0;
    (&_S2310)->differential_0 = _S2306.differential_0;
    DiffPair_float_0 _S2311;
    (&_S2311)->primal_0 = _S2307;
    (&_S2311)->differential_0 = 0.0f;
    DiffPair_float_0 _S2312 = s_primal_ctx_d_atan2_0(&_S2310, &_S2311);
    bool _S2313 = (_S2306.primal_0) < 9.99999997475242708e-07f;
    float k_7;
    float s_diff_k_3;
    float _S2314;
    float _S2315;
    float _S2316;
    float _S2317;
    float _S2318;
    float _S2319;
    float _S2320;
    float _S2321;
    float _S2322;
    float _S2323;
    DiffPair_float_0 _S2324;
    if(_S2313)
    {
        float _S2325 = _S2312.differential_0 * _S2312.primal_0;
        float _S2326 = 1.0f - _S2312.primal_0 * _S2312.primal_0 / 24.0f;
        float _S2327 = 0.0f - (_S2325 + _S2325) * 0.0416666679084301f;
        float _S2328 = _S2307 * _S2307;
        float _S2329 = _S2327 * _S2307;
        float _S2330 = _S2329 / _S2328;
        float _S2331 = _S2328 * _S2328;
        k_7 = _S2326 / _S2307;
        s_diff_k_3 = _S2330;
        _S2314 = _S2331;
        _S2315 = _S2329;
        _S2316 = _S2328;
        _S2317 = _S2326;
        _S2318 = _S2327;
        _S2319 = 0.0f;
        _S2320 = 0.0f;
        _S2321 = 0.0f;
        _S2322 = 0.0f;
        _S2323 = 0.0f;
        (&_S2324)->primal_0 = 0.0f;
        (&_S2324)->differential_0 = 0.0f;
    }
    else
    {
        float _S2332 = 0.5f * _S2312.primal_0;
        float _S2333 = _S2312.differential_0 * 0.5f;
        DiffPair_float_0 _S2334;
        (&_S2334)->primal_0 = _S2332;
        (&_S2334)->differential_0 = _S2333;
        DiffPair_float_0 _S2335 = s_primal_ctx_d_sin_0(&_S2334);
        float _S2336 = 2.0f * _S2335.primal_0;
        float _S2337 = _S2335.differential_0 * 2.0f;
        float _S2338 = _S2306.primal_0 * _S2306.primal_0;
        float _S2339 = _S2337 * _S2306.primal_0 - _S2336 * _S2306.differential_0;
        float _S2340 = _S2339 / _S2338;
        float _S2341 = _S2338 * _S2338;
        k_7 = _S2336 / _S2306.primal_0;
        s_diff_k_3 = _S2340;
        _S2314 = 0.0f;
        _S2315 = 0.0f;
        _S2316 = 0.0f;
        _S2317 = 0.0f;
        _S2318 = 0.0f;
        _S2319 = _S2341;
        _S2320 = _S2339;
        _S2321 = _S2338;
        _S2322 = _S2336;
        _S2323 = _S2337;
        (&_S2324)->primal_0 = _S2332;
        (&_S2324)->differential_0 = _S2333;
    }
    float2  _S2342 = make_float2 (k_7);
    float2  _S2343 = make_float2 (s_diff_k_3);
    float2  _S2344 = _S2302 * make_float2 (k_7);
    float2  _S2345 = _S2303 * make_float2 (k_7) + make_float2 (s_diff_k_3) * _S2302;
    float u_80 = _S2344.x;
    float s_diff_u_24 = _S2345.x;
    float v_80 = _S2344.y;
    float s_diff_v_24 = _S2345.y;
    float _S2346 = s_diff_u_24 * u_80;
    float _S2347 = s_diff_v_24 * v_80;
    float r2_80 = u_80 * u_80 + v_80 * v_80;
    float s_diff_r2_78 = _S2346 + _S2346 + (_S2347 + _S2347);
    float _S2348 = s_diff_r2_78 * dist_coeffs_11[int(3)];
    float _S2349 = dist_coeffs_11[int(2)] + r2_80 * dist_coeffs_11[int(3)];
    float _S2350 = s_diff_r2_78 * _S2349 + _S2348 * r2_80;
    float _S2351 = dist_coeffs_11[int(1)] + r2_80 * _S2349;
    float _S2352 = s_diff_r2_78 * _S2351 + _S2350 * r2_80;
    float _S2353 = dist_coeffs_11[int(0)] + r2_80 * _S2351;
    float _S2354 = s_diff_r2_78 * _S2353 + _S2352 * r2_80;
    float2  _S2355 = make_float2 (_S2354);
    float radial_33 = 1.0f + r2_80 * _S2353;
    float2  _S2356 = make_float2 (radial_33);
    float _S2357 = 2.0f * dist_coeffs_11[int(4)];
    float _S2358 = _S2357 * u_80;
    float _S2359 = s_diff_u_24 * _S2357;
    float _S2360 = 2.0f * u_80;
    float _S2361 = s_diff_u_24 * 2.0f;
    float _S2362 = 2.0f * dist_coeffs_11[int(5)];
    float _S2363 = _S2362 * u_80;
    float _S2364 = s_diff_u_24 * _S2362;
    float _S2365 = 2.0f * v_80;
    float _S2366 = s_diff_v_24 * 2.0f;
    float2  _S2367 = _S2345 * make_float2 (radial_33) + make_float2 (_S2354) * _S2344 + make_float2 (_S2359 * v_80 + s_diff_v_24 * _S2358 + (s_diff_r2_78 + (_S2361 * u_80 + s_diff_u_24 * _S2360)) * dist_coeffs_11[int(5)] + s_diff_r2_78 * dist_coeffs_11[int(6)], _S2364 * v_80 + s_diff_v_24 * _S2363 + (s_diff_r2_78 + (_S2366 * v_80 + s_diff_v_24 * _S2365)) * dist_coeffs_11[int(4)] + s_diff_r2_78 * dist_coeffs_11[int(7)]);
    float2  _S2368 = _S2367 + make_float2 (_S2367.x * dist_coeffs_11[int(8)] + _S2367.y * dist_coeffs_11[int(9)], 0.0f);
    float _S2369 = _S2368.x * fx_11;
    float _S2370 = _S2368.y * fy_11;
    Matrix<float, 2, 3>  _S2371 = J_10;
    *&(((&_S2371)->rows + (int(0)))->x) = _S2369;
    *&(((&_S2371)->rows + (int(1)))->x) = _S2370;
    float2  _S2372 = make_float2 (0.0f, 1.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2373 = { _S2302, _S2372 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2374;
    (&_S2374)->primal_0 = _S2302;
    (&_S2374)->differential_0 = _S2372;
    DiffPair_float_0 _S2375 = s_primal_ctx_s_fwd_length_impl_0(&_S2374);
    DiffPair_float_0 _S2376 = { _S2375.primal_0, _S2375.differential_0 };
    DiffPair_float_0 _S2377;
    (&_S2377)->primal_0 = _S2375.primal_0;
    (&_S2377)->differential_0 = _S2375.differential_0;
    DiffPair_float_0 _S2378;
    (&_S2378)->primal_0 = _S2307;
    (&_S2378)->differential_0 = 0.0f;
    DiffPair_float_0 _S2379 = s_primal_ctx_d_atan2_0(&_S2377, &_S2378);
    bool _S2380 = (_S2375.primal_0) < 9.99999997475242708e-07f;
    float _S2381;
    float _S2382;
    float _S2383;
    float _S2384;
    float _S2385;
    float _S2386;
    float _S2387;
    float _S2388;
    float _S2389;
    float _S2390;
    DiffPair_float_0 _S2391;
    if(_S2380)
    {
        float _S2392 = _S2379.differential_0 * _S2379.primal_0;
        float _S2393 = 1.0f - _S2379.primal_0 * _S2379.primal_0 / 24.0f;
        float _S2394 = 0.0f - (_S2392 + _S2392) * 0.0416666679084301f;
        float _S2395 = _S2307 * _S2307;
        float _S2396 = _S2394 * _S2307;
        float _S2397 = _S2396 / _S2395;
        float _S2398 = _S2395 * _S2395;
        k_7 = _S2393 / _S2307;
        s_diff_k_3 = _S2397;
        _S2381 = _S2398;
        _S2382 = _S2396;
        _S2383 = _S2395;
        _S2384 = _S2393;
        _S2385 = _S2394;
        _S2386 = 0.0f;
        _S2387 = 0.0f;
        _S2388 = 0.0f;
        _S2389 = 0.0f;
        _S2390 = 0.0f;
        (&_S2391)->primal_0 = 0.0f;
        (&_S2391)->differential_0 = 0.0f;
    }
    else
    {
        float _S2399 = 0.5f * _S2379.primal_0;
        float _S2400 = _S2379.differential_0 * 0.5f;
        DiffPair_float_0 _S2401;
        (&_S2401)->primal_0 = _S2399;
        (&_S2401)->differential_0 = _S2400;
        DiffPair_float_0 _S2402 = s_primal_ctx_d_sin_0(&_S2401);
        float _S2403 = 2.0f * _S2402.primal_0;
        float _S2404 = _S2402.differential_0 * 2.0f;
        float _S2405 = _S2375.primal_0 * _S2375.primal_0;
        float _S2406 = _S2404 * _S2375.primal_0 - _S2403 * _S2375.differential_0;
        float _S2407 = _S2406 / _S2405;
        float _S2408 = _S2405 * _S2405;
        k_7 = _S2403 / _S2375.primal_0;
        s_diff_k_3 = _S2407;
        _S2381 = 0.0f;
        _S2382 = 0.0f;
        _S2383 = 0.0f;
        _S2384 = 0.0f;
        _S2385 = 0.0f;
        _S2386 = _S2408;
        _S2387 = _S2406;
        _S2388 = _S2405;
        _S2389 = _S2403;
        _S2390 = _S2404;
        (&_S2391)->primal_0 = _S2399;
        (&_S2391)->differential_0 = _S2400;
    }
    float2  _S2409 = make_float2 (k_7);
    float2  _S2410 = make_float2 (s_diff_k_3);
    float2  _S2411 = _S2302 * make_float2 (k_7);
    float2  _S2412 = _S2372 * make_float2 (k_7) + make_float2 (s_diff_k_3) * _S2302;
    float u_81 = _S2411.x;
    float s_diff_u_25 = _S2412.x;
    float v_81 = _S2411.y;
    float s_diff_v_25 = _S2412.y;
    float _S2413 = s_diff_u_25 * u_81;
    float _S2414 = s_diff_v_25 * v_81;
    float r2_81 = u_81 * u_81 + v_81 * v_81;
    float s_diff_r2_79 = _S2413 + _S2413 + (_S2414 + _S2414);
    float _S2415 = s_diff_r2_79 * dist_coeffs_11[int(3)];
    float _S2416 = dist_coeffs_11[int(2)] + r2_81 * dist_coeffs_11[int(3)];
    float _S2417 = s_diff_r2_79 * _S2416 + _S2415 * r2_81;
    float _S2418 = dist_coeffs_11[int(1)] + r2_81 * _S2416;
    float _S2419 = s_diff_r2_79 * _S2418 + _S2417 * r2_81;
    float _S2420 = dist_coeffs_11[int(0)] + r2_81 * _S2418;
    float _S2421 = s_diff_r2_79 * _S2420 + _S2419 * r2_81;
    float2  _S2422 = make_float2 (_S2421);
    float radial_34 = 1.0f + r2_81 * _S2420;
    float2  _S2423 = make_float2 (radial_34);
    float _S2424 = _S2357 * u_81;
    float _S2425 = s_diff_u_25 * _S2357;
    float _S2426 = 2.0f * u_81;
    float _S2427 = s_diff_u_25 * 2.0f;
    float _S2428 = _S2362 * u_81;
    float _S2429 = s_diff_u_25 * _S2362;
    float _S2430 = 2.0f * v_81;
    float _S2431 = s_diff_v_25 * 2.0f;
    float2  _S2432 = _S2412 * make_float2 (radial_34) + make_float2 (_S2421) * _S2411 + make_float2 (_S2425 * v_81 + s_diff_v_25 * _S2424 + (s_diff_r2_79 + (_S2427 * u_81 + s_diff_u_25 * _S2426)) * dist_coeffs_11[int(5)] + s_diff_r2_79 * dist_coeffs_11[int(6)], _S2429 * v_81 + s_diff_v_25 * _S2428 + (s_diff_r2_79 + (_S2431 * v_81 + s_diff_v_25 * _S2430)) * dist_coeffs_11[int(4)] + s_diff_r2_79 * dist_coeffs_11[int(7)]);
    float2  _S2433 = _S2432 + make_float2 (_S2432.x * dist_coeffs_11[int(8)] + _S2432.y * dist_coeffs_11[int(9)], 0.0f);
    float _S2434 = _S2433.y * fy_11;
    *&(((&_S2371)->rows + (int(0)))->y) = _S2433.x * fx_11;
    *&(((&_S2371)->rows + (int(1)))->y) = _S2434;
    float2  _S2435 = make_float2 (0.0f, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2436 = { _S2302, _S2435 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2437;
    (&_S2437)->primal_0 = _S2302;
    (&_S2437)->differential_0 = _S2435;
    DiffPair_float_0 _S2438 = s_primal_ctx_s_fwd_length_impl_0(&_S2437);
    DiffPair_float_0 _S2439 = { _S2438.primal_0, _S2438.differential_0 };
    DiffPair_float_0 _S2440 = { _S2307, 1.0f };
    DiffPair_float_0 _S2441;
    (&_S2441)->primal_0 = _S2438.primal_0;
    (&_S2441)->differential_0 = _S2438.differential_0;
    DiffPair_float_0 _S2442;
    (&_S2442)->primal_0 = _S2307;
    (&_S2442)->differential_0 = 1.0f;
    DiffPair_float_0 _S2443 = s_primal_ctx_d_atan2_0(&_S2441, &_S2442);
    bool _S2444 = (_S2438.primal_0) < 9.99999997475242708e-07f;
    float _S2445;
    float _S2446;
    float _S2447;
    float _S2448;
    float _S2449;
    float _S2450;
    float _S2451;
    float _S2452;
    float _S2453;
    float _S2454;
    DiffPair_float_0 _S2455;
    if(_S2444)
    {
        float _S2456 = _S2443.differential_0 * _S2443.primal_0;
        float _S2457 = 1.0f - _S2443.primal_0 * _S2443.primal_0 / 24.0f;
        float _S2458 = 0.0f - (_S2456 + _S2456) * 0.0416666679084301f;
        float _S2459 = _S2307 * _S2307;
        float _S2460 = _S2458 * _S2307 - _S2457;
        float _S2461 = _S2460 / _S2459;
        float _S2462 = _S2459 * _S2459;
        k_7 = _S2457 / _S2307;
        s_diff_k_3 = _S2461;
        _S2445 = _S2462;
        _S2446 = _S2460;
        _S2447 = _S2459;
        _S2448 = _S2457;
        _S2449 = _S2458;
        _S2450 = 0.0f;
        _S2451 = 0.0f;
        _S2452 = 0.0f;
        _S2453 = 0.0f;
        _S2454 = 0.0f;
        (&_S2455)->primal_0 = 0.0f;
        (&_S2455)->differential_0 = 0.0f;
    }
    else
    {
        float _S2463 = 0.5f * _S2443.primal_0;
        float _S2464 = _S2443.differential_0 * 0.5f;
        DiffPair_float_0 _S2465;
        (&_S2465)->primal_0 = _S2463;
        (&_S2465)->differential_0 = _S2464;
        DiffPair_float_0 _S2466 = s_primal_ctx_d_sin_0(&_S2465);
        float _S2467 = 2.0f * _S2466.primal_0;
        float _S2468 = _S2466.differential_0 * 2.0f;
        float _S2469 = _S2438.primal_0 * _S2438.primal_0;
        float _S2470 = _S2468 * _S2438.primal_0 - _S2467 * _S2438.differential_0;
        float _S2471 = _S2470 / _S2469;
        float _S2472 = _S2469 * _S2469;
        k_7 = _S2467 / _S2438.primal_0;
        s_diff_k_3 = _S2471;
        _S2445 = 0.0f;
        _S2446 = 0.0f;
        _S2447 = 0.0f;
        _S2448 = 0.0f;
        _S2449 = 0.0f;
        _S2450 = _S2472;
        _S2451 = _S2470;
        _S2452 = _S2469;
        _S2453 = _S2467;
        _S2454 = _S2468;
        (&_S2455)->primal_0 = _S2463;
        (&_S2455)->differential_0 = _S2464;
    }
    float2  _S2473 = make_float2 (k_7);
    float2  _S2474 = make_float2 (s_diff_k_3);
    float2  _S2475 = _S2302 * make_float2 (k_7);
    float2  _S2476 = make_float2 (s_diff_k_3) * _S2302;
    float u_82 = _S2475.x;
    float s_diff_u_26 = _S2476.x;
    float v_82 = _S2475.y;
    float s_diff_v_26 = _S2476.y;
    float _S2477 = s_diff_u_26 * u_82;
    float _S2478 = s_diff_v_26 * v_82;
    float r2_82 = u_82 * u_82 + v_82 * v_82;
    float s_diff_r2_80 = _S2477 + _S2477 + (_S2478 + _S2478);
    float _S2479 = s_diff_r2_80 * dist_coeffs_11[int(3)];
    float _S2480 = dist_coeffs_11[int(2)] + r2_82 * dist_coeffs_11[int(3)];
    float _S2481 = s_diff_r2_80 * _S2480 + _S2479 * r2_82;
    float _S2482 = dist_coeffs_11[int(1)] + r2_82 * _S2480;
    float _S2483 = s_diff_r2_80 * _S2482 + _S2481 * r2_82;
    float _S2484 = dist_coeffs_11[int(0)] + r2_82 * _S2482;
    float _S2485 = s_diff_r2_80 * _S2484 + _S2483 * r2_82;
    float2  _S2486 = make_float2 (_S2485);
    float radial_35 = 1.0f + r2_82 * _S2484;
    float2  _S2487 = make_float2 (radial_35);
    float _S2488 = _S2357 * u_82;
    float _S2489 = s_diff_u_26 * _S2357;
    float _S2490 = 2.0f * u_82;
    float _S2491 = s_diff_u_26 * 2.0f;
    float _S2492 = _S2362 * u_82;
    float _S2493 = s_diff_u_26 * _S2362;
    float _S2494 = 2.0f * v_82;
    float _S2495 = s_diff_v_26 * 2.0f;
    float2  _S2496 = _S2476 * make_float2 (radial_35) + make_float2 (_S2485) * _S2475 + make_float2 (_S2489 * v_82 + s_diff_v_26 * _S2488 + (s_diff_r2_80 + (_S2491 * u_82 + s_diff_u_26 * _S2490)) * dist_coeffs_11[int(5)] + s_diff_r2_80 * dist_coeffs_11[int(6)], _S2493 * v_82 + s_diff_v_26 * _S2492 + (s_diff_r2_80 + (_S2495 * v_82 + s_diff_v_26 * _S2494)) * dist_coeffs_11[int(4)] + s_diff_r2_80 * dist_coeffs_11[int(7)]);
    float2  _S2497 = _S2496 + make_float2 (_S2496.x * dist_coeffs_11[int(8)] + _S2496.y * dist_coeffs_11[int(9)], 0.0f);
    float _S2498 = _S2497.y * fy_11;
    *&(((&_S2371)->rows + (int(0)))->z) = _S2497.x * fx_11;
    *&(((&_S2371)->rows + (int(1)))->z) = _S2498;
    Matrix<float, 2, 3>  _S2499 = s_primal_ctx_mul_2(_S2371, _S2301);
    Matrix<float, 3, 2>  _S2500 = transpose_1(_S2371);
    Matrix<float, 2, 2>  _S2501 = s_primal_ctx_mul_3(_S2499, _S2500);
    float eps2d_8;
    if(antialiased_8)
    {
        eps2d_8 = 0.10000000149011612f;
    }
    else
    {
        eps2d_8 = 0.30000001192092896f;
    }
    float _S2502 = _S2501.rows[int(0)].y * _S2501.rows[int(1)].x;
    float det_orig_8 = _S2501.rows[int(0)].x * _S2501.rows[int(1)].y - _S2502;
    float _S2503 = _S2501.rows[int(0)].x + eps2d_8;
    Matrix<float, 2, 2>  _S2504 = _S2501;
    *&(((&_S2504)->rows + (int(0)))->x) = _S2503;
    float _S2505 = _S2501.rows[int(1)].y + eps2d_8;
    *&(((&_S2504)->rows + (int(1)))->y) = _S2505;
    Matrix<float, 2, 2>  _S2506 = _S2504;
    Matrix<float, 2, 2>  _S2507 = _S2504;
    float det_blur_8 = _S2503 * _S2505 - _S2502;
    float _S2508 = det_orig_8 / det_blur_8;
    float _S2509 = det_blur_8 * det_blur_8;
    float _S2510 = (F32_max((0.0f), (_S2508)));
    float _S2511 = s_primal_ctx_sqrt_0(_S2510);
    float invdet_10 = 1.0f / det_blur_8;
    float _S2512 = - _S2501.rows[int(0)].y;
    float _S2513 = - _S2501.rows[int(1)].x;
    if(antialiased_8)
    {
        k_7 = _S2290 * _S2511;
    }
    else
    {
        k_7 = _S2290;
    }
    float _S2514 = k_7 / 0.00392156885936856f;
    float _S2515 = 2.0f * s_primal_ctx_log_0(_S2514);
    float _S2516 = s_primal_ctx_sqrt_0(_S2515);
    float _S2517 = _S2506.rows[int(0)].x;
    float _S2518 = _S2507.rows[int(1)].y;
    float3  campos_3 = - s_primal_ctx_mul_0(_S2300, t_8);
    float3  _S2519 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2520;
    (&_S2520)->primal_0 = mean_9;
    (&_S2520)->differential_0 = _S2519;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2521;
    (&_S2521)->primal_0 = scale_8;
    (&_S2521)->differential_0 = _S2519;
    DiffPair_float_0 _S2522;
    (&_S2522)->primal_0 = in_opacity_8;
    (&_S2522)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2523;
    (&_S2523)->primal_0 = campos_3;
    (&_S2523)->differential_0 = _S2519;
    s_bwd_prop_view_radius_3dgs_0(&_S2520, &_S2521, &_S2522, &_S2523, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2524 = _S2520;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2525 = _S2521;
    DiffPair_float_0 _S2526 = _S2522;
    float2  _S2527 = make_float2 (0.0f);
    float2  _S2528 = _S2527;
    *&((&_S2528)->y) = v_conic_2.z;
    float2  _S2529 = _S2527;
    *&((&_S2529)->y) = v_conic_2.y;
    *&((&_S2529)->x) = v_conic_2.x;
    DiffPair_float_0 _S2530;
    (&_S2530)->primal_0 = _S2518;
    (&_S2530)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2530, 0.0f);
    DiffPair_float_0 _S2531;
    (&_S2531)->primal_0 = _S2517;
    (&_S2531)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2531, 0.0f);
    DiffPair_float_0 _S2532;
    (&_S2532)->primal_0 = 3.32999992370605469f;
    (&_S2532)->differential_0 = 0.0f;
    DiffPair_float_0 _S2533;
    (&_S2533)->primal_0 = _S2516;
    (&_S2533)->differential_0 = 0.0f;
    _d_min_0(&_S2532, &_S2533, 0.0f);
    DiffPair_float_0 _S2534;
    (&_S2534)->primal_0 = _S2515;
    (&_S2534)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2534, _S2533.differential_0);
    float _S2535 = 2.0f * _S2534.differential_0;
    DiffPair_float_0 _S2536;
    (&_S2536)->primal_0 = _S2514;
    (&_S2536)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2536, _S2535);
    float _S2537 = v_opacity_2 + 254.9999847412109375f * _S2536.differential_0;
    float2  _S2538 = make_float2 (_S2531.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S2539 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2540 = _S2539;
    _S2540[int(1)] = _S2528;
    _S2540[int(0)] = _S2529;
    Matrix<float, 2, 2>  _S2541 = _S2540;
    float2  _S2542 = make_float2 (0.0f, _S2530.differential_0);
    if(antialiased_8)
    {
        float _S2543 = _S2511 * _S2537;
        k_7 = _S2290 * _S2537;
        s_diff_k_3 = _S2543;
    }
    else
    {
        k_7 = 0.0f;
        s_diff_k_3 = _S2537;
    }
    float _S2544 = invdet_10 * _S2541.rows[int(1)].y;
    float _S2545 = - (invdet_10 * _S2541.rows[int(1)].x);
    float _S2546 = - (invdet_10 * _S2541.rows[int(0)].y);
    float _S2547 = invdet_10 * _S2541.rows[int(0)].x;
    float _S2548 = - ((_S2503 * _S2541.rows[int(1)].y + _S2513 * _S2541.rows[int(1)].x + _S2512 * _S2541.rows[int(0)].y + _S2505 * _S2541.rows[int(0)].x) / _S2509);
    DiffPair_float_0 _S2549;
    (&_S2549)->primal_0 = _S2510;
    (&_S2549)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2549, k_7);
    DiffPair_float_0 _S2550;
    (&_S2550)->primal_0 = 0.0f;
    (&_S2550)->differential_0 = 0.0f;
    DiffPair_float_0 _S2551;
    (&_S2551)->primal_0 = _S2508;
    (&_S2551)->differential_0 = 0.0f;
    _d_max_0(&_S2550, &_S2551, _S2549.differential_0);
    float _S2552 = _S2551.differential_0 / _S2509;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2552;
    float _S2553 = det_orig_8 * - _S2552 + _S2548;
    float _S2554 = - _S2553;
    float _S2555 = _S2503 * _S2553;
    float _S2556 = _S2505 * _S2553;
    Matrix<float, 2, 2>  _S2557 = _S2539;
    _S2557[int(1)] = _S2542;
    _S2557[int(0)] = _S2538;
    _S2504 = _S2557;
    *&(((&_S2504)->rows + (int(1)))->y) = 0.0f;
    float _S2558 = _S2555 + _S2557.rows[int(1)].y + _S2547;
    *&(((&_S2504)->rows + (int(0)))->x) = 0.0f;
    float _S2559 = _S2556 + _S2557.rows[int(0)].x + _S2544;
    float _S2560 = _S2554 + - s_diff_det_orig_T_2;
    float _S2561 = _S2501.rows[int(0)].y * _S2560 + _S2545;
    float _S2562 = _S2501.rows[int(1)].x * _S2560 + _S2546;
    float _S2563 = _S2501.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2564 = _S2558 + _S2501.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2565 = _S2527;
    *&((&_S2565)->x) = _S2561;
    *&((&_S2565)->y) = _S2564;
    float _S2566 = _S2559 + _S2563;
    float2  _S2567 = _S2527;
    *&((&_S2567)->y) = _S2562;
    *&((&_S2567)->x) = _S2566;
    Matrix<float, 2, 2>  _S2568 = _S2539;
    _S2568[int(1)] = _S2565;
    _S2568[int(0)] = _S2567;
    Matrix<float, 2, 2>  _S2569 = _S2504 + _S2568;
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2570;
    (&_S2570)->primal_0 = _S2499;
    (&_S2570)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S2571 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2572;
    (&_S2572)->primal_0 = _S2500;
    (&_S2572)->differential_0 = _S2571;
    s_bwd_prop_mul_0(&_S2570, &_S2572, _S2569);
    Matrix<float, 2, 3>  _S2573 = transpose_2(_S2572.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2574;
    (&_S2574)->primal_0 = _S2371;
    (&_S2574)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S2575 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2576;
    (&_S2576)->primal_0 = _S2301;
    (&_S2576)->differential_0 = _S2575;
    s_bwd_prop_mul_1(&_S2574, &_S2576, _S2570.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2577 = _S2576;
    Matrix<float, 2, 3>  _S2578 = _S2573 + _S2574.differential_0;
    float2  _S2579 = make_float2 (0.0f, _S2578.rows[int(1)].z) + make_float2 (_S2578.rows[int(0)].z, 0.0f);
    float _S2580 = fx_11 * _S2579.x;
    float2  _S2581 = make_float2 (_S2580, fy_11 * _S2579.y) + make_float2 (dist_coeffs_11[int(8)] * _S2580, dist_coeffs_11[int(9)] * _S2580);
    float2  _S2582 = _S2475 * _S2581;
    float2  _S2583 = _S2476 * _S2581;
    float _S2584 = dist_coeffs_11[int(4)] * _S2581.y;
    float _S2585 = dist_coeffs_11[int(5)] * _S2581.x;
    float _S2586 = _S2583.x + _S2583.y;
    float _S2587 = _S2582.x + _S2582.y;
    float _S2588 = r2_82 * _S2587;
    float _S2589 = s_diff_r2_80 * _S2587 + r2_82 * _S2586;
    float _S2590 = r2_82 * _S2588;
    float _S2591 = s_diff_r2_80 * _S2588 + r2_82 * _S2589;
    float _S2592 = dist_coeffs_11[int(7)] * _S2581.y + _S2584 + dist_coeffs_11[int(6)] * _S2581.x + _S2585 + _S2484 * _S2587 + _S2482 * _S2588 + _S2480 * _S2590 + dist_coeffs_11[int(3)] * (r2_82 * _S2590);
    float _S2593 = _S2483 * _S2587 + _S2484 * _S2586 + _S2481 * _S2588 + _S2482 * _S2589 + _S2479 * _S2590 + _S2480 * _S2591 + dist_coeffs_11[int(3)] * (s_diff_r2_80 * _S2590 + r2_82 * _S2591);
    float _S2594 = _S2592 + _S2592;
    float _S2595 = v_82 * _S2593;
    float _S2596 = u_82 * _S2593;
    float2  _S2597 = _S2486 * _S2581 + make_float2 (_S2362 * (s_diff_v_26 * _S2581.y) + _S2491 * _S2585 + 2.0f * (s_diff_u_26 * _S2585) + _S2357 * (s_diff_v_26 * _S2581.x) + s_diff_u_26 * _S2594 + _S2596 + _S2596, _S2495 * _S2584 + 2.0f * (s_diff_v_26 * _S2584) + _S2493 * _S2581.y + _S2489 * _S2581.x + s_diff_v_26 * _S2594 + _S2595 + _S2595);
    float2  _S2598 = _S2487 * _S2581 + make_float2 (_S2362 * (v_82 * _S2581.y) + _S2490 * _S2585 + 2.0f * (u_82 * _S2585) + _S2357 * (v_82 * _S2581.x) + u_82 * _S2594, _S2494 * _S2584 + 2.0f * (v_82 * _S2584) + _S2492 * _S2581.y + _S2488 * _S2581.x + v_82 * _S2594);
    float2  _S2599 = _S2302 * _S2598;
    float2  _S2600 = _S2302 * _S2597;
    float _S2601 = _S2600.x + _S2600.y;
    float _S2602 = _S2599.x + _S2599.y;
    float2  _S2603 = _S2474 * _S2598 + _S2473 * _S2597;
    if(_S2444)
    {
        float _S2604 = _S2602 / _S2445;
        float _S2605 = _S2447 * _S2604;
        float _S2606 = _S2307 * (_S2446 * - _S2604);
        float _S2607 = _S2601 / _S2447;
        float _S2608 = 0.0416666679084301f * - (_S2307 * _S2605);
        float _S2609 = _S2608 + _S2608;
        float _S2610 = _S2443.primal_0 * (0.0416666679084301f * - (- _S2605 + _S2307 * _S2607));
        float _S2611 = _S2606 + _S2606 + _S2449 * _S2605 + _S2448 * - _S2607;
        float _S2612 = _S2443.differential_0 * _S2609 + _S2610 + _S2610;
        k_7 = _S2443.primal_0 * _S2609;
        _S2445 = _S2612;
        _S2446 = _S2611;
        _S2447 = 0.0f;
        _S2448 = 0.0f;
    }
    else
    {
        float _S2613 = _S2602 / _S2450;
        float _S2614 = _S2452 * _S2613;
        float _S2615 = _S2438.primal_0 * (_S2451 * - _S2613);
        float _S2616 = - _S2614;
        float _S2617 = _S2453 * _S2616;
        float _S2618 = _S2454 * _S2614;
        float _S2619 = _S2601 / _S2452;
        float _S2620 = _S2453 * - _S2619;
        float _S2621 = 2.0f * (_S2438.primal_0 * _S2614);
        float _S2622 = 2.0f * (_S2438.differential_0 * _S2616 + _S2438.primal_0 * _S2619);
        DiffPair_0 _S2623;
        (&_S2623)->primal_0 = _S2455;
        (&_S2623)->differential_0 = _S2287;
        DiffPair_float_0 _S2624;
        (&_S2624)->primal_0 = _S2622;
        (&_S2624)->differential_0 = _S2621;
        s_bwd_prop_d_sin_0(&_S2623, &_S2624);
        float _S2625 = 0.5f * _S2623.differential_0.primal_0;
        float _S2626 = _S2615 + _S2615 + _S2618 + _S2620;
        k_7 = 0.5f * _S2623.differential_0.differential_0;
        _S2445 = _S2625;
        _S2446 = 0.0f;
        _S2447 = _S2617;
        _S2448 = _S2626;
    }
    DiffPair_0 _S2627;
    (&_S2627)->primal_0 = _S2439;
    (&_S2627)->differential_0 = _S2287;
    DiffPair_0 _S2628;
    (&_S2628)->primal_0 = _S2440;
    (&_S2628)->differential_0 = _S2287;
    DiffPair_float_0 _S2629;
    (&_S2629)->primal_0 = _S2445;
    (&_S2629)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2627, &_S2628, &_S2629);
    float _S2630 = _S2628.differential_0.primal_0 + _S2446;
    float _S2631 = _S2627.differential_0.differential_0 + _S2447;
    float _S2632 = _S2627.differential_0.primal_0 + _S2448;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2633 = { _S2527, _S2527 };
    DiffPair_1 _S2634;
    (&_S2634)->primal_0 = _S2436;
    (&_S2634)->differential_0 = _S2633;
    DiffPair_float_0 _S2635;
    (&_S2635)->primal_0 = _S2632;
    (&_S2635)->differential_0 = _S2631;
    s_bwd_prop_s_fwd_length_impl_0(&_S2634, &_S2635);
    float2  _S2636 = _S2634.differential_0.primal_0 + _S2603;
    float3  _S2637 = make_float3 (_S2636.x, _S2636.y, _S2630);
    float2  _S2638 = make_float2 (0.0f, _S2578.rows[int(1)].y) + make_float2 (_S2578.rows[int(0)].y, 0.0f);
    float _S2639 = fx_11 * _S2638.x;
    float2  _S2640 = make_float2 (_S2639, fy_11 * _S2638.y) + make_float2 (dist_coeffs_11[int(8)] * _S2639, dist_coeffs_11[int(9)] * _S2639);
    float2  _S2641 = _S2411 * _S2640;
    float2  _S2642 = _S2412 * _S2640;
    float _S2643 = dist_coeffs_11[int(4)] * _S2640.y;
    float _S2644 = dist_coeffs_11[int(5)] * _S2640.x;
    float _S2645 = _S2642.x + _S2642.y;
    float _S2646 = _S2641.x + _S2641.y;
    float _S2647 = r2_81 * _S2646;
    float _S2648 = s_diff_r2_79 * _S2646 + r2_81 * _S2645;
    float _S2649 = r2_81 * _S2647;
    float _S2650 = s_diff_r2_79 * _S2647 + r2_81 * _S2648;
    float _S2651 = dist_coeffs_11[int(7)] * _S2640.y + _S2643 + dist_coeffs_11[int(6)] * _S2640.x + _S2644 + _S2420 * _S2646 + _S2418 * _S2647 + _S2416 * _S2649 + dist_coeffs_11[int(3)] * (r2_81 * _S2649);
    float _S2652 = _S2419 * _S2646 + _S2420 * _S2645 + _S2417 * _S2647 + _S2418 * _S2648 + _S2415 * _S2649 + _S2416 * _S2650 + dist_coeffs_11[int(3)] * (s_diff_r2_79 * _S2649 + r2_81 * _S2650);
    float _S2653 = _S2651 + _S2651;
    float _S2654 = v_81 * _S2652;
    float _S2655 = u_81 * _S2652;
    float2  _S2656 = _S2422 * _S2640 + make_float2 (_S2362 * (s_diff_v_25 * _S2640.y) + _S2427 * _S2644 + 2.0f * (s_diff_u_25 * _S2644) + _S2357 * (s_diff_v_25 * _S2640.x) + s_diff_u_25 * _S2653 + _S2655 + _S2655, _S2431 * _S2643 + 2.0f * (s_diff_v_25 * _S2643) + _S2429 * _S2640.y + _S2425 * _S2640.x + s_diff_v_25 * _S2653 + _S2654 + _S2654);
    float2  _S2657 = _S2423 * _S2640 + make_float2 (_S2362 * (v_81 * _S2640.y) + _S2426 * _S2644 + 2.0f * (u_81 * _S2644) + _S2357 * (v_81 * _S2640.x) + u_81 * _S2653, _S2430 * _S2643 + 2.0f * (v_81 * _S2643) + _S2428 * _S2640.y + _S2424 * _S2640.x + v_81 * _S2653);
    float2  _S2658 = _S2302 * _S2657;
    float2  _S2659 = _S2372 * _S2657;
    float2  _S2660 = _S2302 * _S2656;
    float _S2661 = _S2659.x + _S2659.y + _S2660.x + _S2660.y;
    float _S2662 = _S2658.x + _S2658.y;
    float2  _S2663 = _S2410 * _S2657 + _S2409 * _S2656;
    if(_S2380)
    {
        float _S2664 = _S2662 / _S2381;
        float _S2665 = _S2383 * _S2664;
        float _S2666 = _S2307 * (_S2382 * - _S2664);
        float _S2667 = _S2661 / _S2383;
        float _S2668 = 0.0416666679084301f * - (_S2307 * _S2665);
        float _S2669 = _S2668 + _S2668;
        float _S2670 = _S2379.primal_0 * (0.0416666679084301f * - (_S2307 * _S2667));
        float _S2671 = _S2666 + _S2666 + _S2385 * _S2665 + _S2384 * - _S2667;
        float _S2672 = _S2379.differential_0 * _S2669 + _S2670 + _S2670;
        k_7 = _S2379.primal_0 * _S2669;
        _S2381 = _S2672;
        _S2382 = _S2671;
        _S2383 = 0.0f;
        _S2384 = 0.0f;
    }
    else
    {
        float _S2673 = _S2662 / _S2386;
        float _S2674 = _S2388 * _S2673;
        float _S2675 = _S2375.primal_0 * (_S2387 * - _S2673);
        float _S2676 = - _S2674;
        float _S2677 = _S2389 * _S2676;
        float _S2678 = _S2390 * _S2674;
        float _S2679 = _S2661 / _S2388;
        float _S2680 = _S2389 * - _S2679;
        float _S2681 = 2.0f * (_S2375.primal_0 * _S2674);
        float _S2682 = 2.0f * (_S2375.differential_0 * _S2676 + _S2375.primal_0 * _S2679);
        DiffPair_0 _S2683;
        (&_S2683)->primal_0 = _S2391;
        (&_S2683)->differential_0 = _S2287;
        DiffPair_float_0 _S2684;
        (&_S2684)->primal_0 = _S2682;
        (&_S2684)->differential_0 = _S2681;
        s_bwd_prop_d_sin_0(&_S2683, &_S2684);
        float _S2685 = 0.5f * _S2683.differential_0.primal_0;
        float _S2686 = _S2675 + _S2675 + _S2678 + _S2680;
        k_7 = 0.5f * _S2683.differential_0.differential_0;
        _S2381 = _S2685;
        _S2382 = 0.0f;
        _S2383 = _S2677;
        _S2384 = _S2686;
    }
    DiffPair_0 _S2687;
    (&_S2687)->primal_0 = _S2376;
    (&_S2687)->differential_0 = _S2287;
    DiffPair_0 _S2688;
    (&_S2688)->primal_0 = _S2309;
    (&_S2688)->differential_0 = _S2287;
    DiffPair_float_0 _S2689;
    (&_S2689)->primal_0 = _S2381;
    (&_S2689)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2687, &_S2688, &_S2689);
    float _S2690 = _S2688.differential_0.primal_0 + _S2382;
    float _S2691 = _S2687.differential_0.differential_0 + _S2383;
    float _S2692 = _S2687.differential_0.primal_0 + _S2384;
    DiffPair_1 _S2693;
    (&_S2693)->primal_0 = _S2373;
    (&_S2693)->differential_0 = _S2633;
    DiffPair_float_0 _S2694;
    (&_S2694)->primal_0 = _S2692;
    (&_S2694)->differential_0 = _S2691;
    s_bwd_prop_s_fwd_length_impl_0(&_S2693, &_S2694);
    float2  _S2695 = _S2693.differential_0.primal_0 + _S2663;
    float2  _S2696 = make_float2 (0.0f, _S2578.rows[int(1)].x) + make_float2 (_S2578.rows[int(0)].x, 0.0f);
    float _S2697 = fx_11 * _S2696.x;
    float2  _S2698 = make_float2 (_S2697, fy_11 * _S2696.y) + make_float2 (dist_coeffs_11[int(8)] * _S2697, dist_coeffs_11[int(9)] * _S2697);
    float2  _S2699 = _S2344 * _S2698;
    float2  _S2700 = _S2345 * _S2698;
    float _S2701 = dist_coeffs_11[int(4)] * _S2698.y;
    float _S2702 = dist_coeffs_11[int(5)] * _S2698.x;
    float _S2703 = _S2700.x + _S2700.y;
    float _S2704 = _S2699.x + _S2699.y;
    float _S2705 = r2_80 * _S2704;
    float _S2706 = s_diff_r2_78 * _S2704 + r2_80 * _S2703;
    float _S2707 = r2_80 * _S2705;
    float _S2708 = s_diff_r2_78 * _S2705 + r2_80 * _S2706;
    float _S2709 = dist_coeffs_11[int(7)] * _S2698.y + _S2701 + dist_coeffs_11[int(6)] * _S2698.x + _S2702 + _S2353 * _S2704 + _S2351 * _S2705 + _S2349 * _S2707 + dist_coeffs_11[int(3)] * (r2_80 * _S2707);
    float _S2710 = _S2352 * _S2704 + _S2353 * _S2703 + _S2350 * _S2705 + _S2351 * _S2706 + _S2348 * _S2707 + _S2349 * _S2708 + dist_coeffs_11[int(3)] * (s_diff_r2_78 * _S2707 + r2_80 * _S2708);
    float _S2711 = _S2709 + _S2709;
    float _S2712 = v_80 * _S2710;
    float _S2713 = u_80 * _S2710;
    float2  _S2714 = _S2355 * _S2698 + make_float2 (_S2362 * (s_diff_v_24 * _S2698.y) + _S2361 * _S2702 + 2.0f * (s_diff_u_24 * _S2702) + _S2357 * (s_diff_v_24 * _S2698.x) + s_diff_u_24 * _S2711 + _S2713 + _S2713, _S2366 * _S2701 + 2.0f * (s_diff_v_24 * _S2701) + _S2364 * _S2698.y + _S2359 * _S2698.x + s_diff_v_24 * _S2711 + _S2712 + _S2712);
    float2  _S2715 = _S2356 * _S2698 + make_float2 (_S2362 * (v_80 * _S2698.y) + _S2360 * _S2702 + 2.0f * (u_80 * _S2702) + _S2357 * (v_80 * _S2698.x) + u_80 * _S2711, _S2365 * _S2701 + 2.0f * (v_80 * _S2701) + _S2363 * _S2698.y + _S2358 * _S2698.x + v_80 * _S2711);
    float3  _S2716 = make_float3 (_S2695.x, _S2695.y, _S2690) + _S2637;
    float2  _S2717 = _S2302 * _S2715;
    float2  _S2718 = _S2303 * _S2715;
    float2  _S2719 = _S2302 * _S2714;
    float _S2720 = _S2718.x + _S2718.y + _S2719.x + _S2719.y;
    float _S2721 = _S2717.x + _S2717.y;
    float2  _S2722 = _S2343 * _S2715 + _S2342 * _S2714;
    if(_S2313)
    {
        float _S2723 = _S2721 / _S2314;
        float _S2724 = _S2316 * _S2723;
        float _S2725 = _S2307 * (_S2315 * - _S2723);
        float _S2726 = _S2720 / _S2316;
        float _S2727 = 0.0416666679084301f * - (_S2307 * _S2724);
        float _S2728 = _S2727 + _S2727;
        float _S2729 = _S2312.primal_0 * (0.0416666679084301f * - (_S2307 * _S2726));
        float _S2730 = _S2725 + _S2725 + _S2318 * _S2724 + _S2317 * - _S2726;
        float _S2731 = _S2312.differential_0 * _S2728 + _S2729 + _S2729;
        k_7 = _S2312.primal_0 * _S2728;
        _S2314 = _S2731;
        _S2315 = _S2730;
        _S2316 = 0.0f;
        _S2317 = 0.0f;
    }
    else
    {
        float _S2732 = _S2721 / _S2319;
        float _S2733 = _S2321 * _S2732;
        float _S2734 = _S2306.primal_0 * (_S2320 * - _S2732);
        float _S2735 = - _S2733;
        float _S2736 = _S2322 * _S2735;
        float _S2737 = _S2323 * _S2733;
        float _S2738 = _S2720 / _S2321;
        float _S2739 = _S2322 * - _S2738;
        float _S2740 = 2.0f * (_S2306.primal_0 * _S2733);
        float _S2741 = 2.0f * (_S2306.differential_0 * _S2735 + _S2306.primal_0 * _S2738);
        DiffPair_0 _S2742;
        (&_S2742)->primal_0 = _S2324;
        (&_S2742)->differential_0 = _S2287;
        DiffPair_float_0 _S2743;
        (&_S2743)->primal_0 = _S2741;
        (&_S2743)->differential_0 = _S2740;
        s_bwd_prop_d_sin_0(&_S2742, &_S2743);
        float _S2744 = 0.5f * _S2742.differential_0.primal_0;
        float _S2745 = _S2734 + _S2734 + _S2737 + _S2739;
        k_7 = 0.5f * _S2742.differential_0.differential_0;
        _S2314 = _S2744;
        _S2315 = 0.0f;
        _S2316 = _S2736;
        _S2317 = _S2745;
    }
    DiffPair_0 _S2746;
    (&_S2746)->primal_0 = _S2308;
    (&_S2746)->differential_0 = _S2287;
    DiffPair_0 _S2747;
    (&_S2747)->primal_0 = _S2309;
    (&_S2747)->differential_0 = _S2287;
    DiffPair_float_0 _S2748;
    (&_S2748)->primal_0 = _S2314;
    (&_S2748)->differential_0 = k_7;
    s_bwd_prop_d_atan2_0(&_S2746, &_S2747, &_S2748);
    float _S2749 = _S2747.differential_0.primal_0 + _S2315;
    float _S2750 = _S2746.differential_0.differential_0 + _S2316;
    float _S2751 = _S2746.differential_0.primal_0 + _S2317;
    DiffPair_1 _S2752;
    (&_S2752)->primal_0 = _S2304;
    (&_S2752)->differential_0 = _S2633;
    DiffPair_float_0 _S2753;
    (&_S2753)->primal_0 = _S2751;
    (&_S2753)->differential_0 = _S2750;
    s_bwd_prop_s_fwd_length_impl_0(&_S2752, &_S2753);
    float2  _S2754 = _S2752.differential_0.primal_0 + _S2722;
    float3  _S2755 = make_float3 (_S2754.x, _S2754.y, _S2749);
    float _S2756 = length_0(_S2302);
    float _S2757 = s_primal_ctx_atan2_0(_S2756, _S2307);
    bool _S2758 = _S2756 < 9.99999997475242708e-07f;
    if(_S2758)
    {
        float _S2759 = 1.0f - _S2757 * _S2757 / 24.0f;
        float _S2760 = _S2307 * _S2307;
        k_7 = _S2759 / _S2307;
        _S2314 = _S2760;
        _S2315 = _S2759;
        _S2316 = 0.0f;
        _S2317 = 0.0f;
        _S2318 = 0.0f;
    }
    else
    {
        float _S2761 = 0.5f * _S2757;
        float _S2762 = 2.0f * s_primal_ctx_sin_0(_S2761);
        float _S2763 = _S2756 * _S2756;
        k_7 = _S2762 / _S2756;
        _S2314 = 0.0f;
        _S2315 = 0.0f;
        _S2316 = _S2763;
        _S2317 = _S2762;
        _S2318 = _S2761;
    }
    float2  _S2764 = make_float2 (k_7);
    float2  _S2765 = _S2302 * make_float2 (k_7);
    float _S2766 = fx_11 * v_mean2d_2.x;
    float u_83 = _S2765.x;
    float v_83 = _S2765.y;
    float r2_83 = u_83 * u_83 + v_83 * v_83;
    float _S2767 = dist_coeffs_11[int(2)] + r2_83 * dist_coeffs_11[int(3)];
    float _S2768 = dist_coeffs_11[int(1)] + r2_83 * _S2767;
    float _S2769 = dist_coeffs_11[int(0)] + r2_83 * _S2768;
    float2  _S2770 = make_float2 (_S2766, fy_11 * v_mean2d_2.y) + make_float2 (dist_coeffs_11[int(8)] * _S2766, dist_coeffs_11[int(9)] * _S2766);
    float2  _S2771 = _S2765 * _S2770;
    float _S2772 = dist_coeffs_11[int(4)] * _S2770.y;
    float _S2773 = dist_coeffs_11[int(5)] * _S2770.x;
    float _S2774 = _S2771.x + _S2771.y;
    float _S2775 = r2_83 * _S2774;
    float _S2776 = r2_83 * _S2775;
    float _S2777 = dist_coeffs_11[int(7)] * _S2770.y + _S2772 + dist_coeffs_11[int(6)] * _S2770.x + _S2773 + _S2769 * _S2774 + _S2768 * _S2775 + _S2767 * _S2776 + dist_coeffs_11[int(3)] * (r2_83 * _S2776);
    float _S2778 = v_83 * _S2777;
    float _S2779 = u_83 * _S2777;
    float2  _S2780 = make_float2 (1.0f + r2_83 * _S2769) * _S2770 + make_float2 (_S2362 * (v_83 * _S2770.y) + 2.0f * u_83 * _S2773 + 2.0f * (u_83 * _S2773) + _S2357 * (v_83 * _S2770.x) + _S2779 + _S2779, 2.0f * v_83 * _S2772 + 2.0f * (v_83 * _S2772) + _S2362 * u_83 * _S2770.y + _S2357 * u_83 * _S2770.x + _S2778 + _S2778);
    float2  _S2781 = _S2302 * _S2780;
    float2  _S2782 = _S2764 * _S2780;
    float _S2783 = _S2781.x + _S2781.y;
    if(_S2758)
    {
        float _S2784 = _S2783 / _S2314;
        float _S2785 = _S2315 * - _S2784;
        float _S2786 = _S2757 * (0.0416666679084301f * - (_S2307 * _S2784));
        k_7 = _S2786 + _S2786;
        _S2314 = _S2785;
        _S2315 = 0.0f;
    }
    else
    {
        float _S2787 = _S2783 / _S2316;
        float _S2788 = _S2317 * - _S2787;
        float _S2789 = 2.0f * (_S2756 * _S2787);
        DiffPair_float_0 _S2790;
        (&_S2790)->primal_0 = _S2318;
        (&_S2790)->differential_0 = 0.0f;
        s_bwd_prop_sin_0(&_S2790, _S2789);
        k_7 = 0.5f * _S2790.differential_0;
        _S2314 = 0.0f;
        _S2315 = _S2788;
    }
    DiffPair_float_0 _S2791;
    (&_S2791)->primal_0 = _S2756;
    (&_S2791)->differential_0 = 0.0f;
    DiffPair_float_0 _S2792;
    (&_S2792)->primal_0 = _S2307;
    (&_S2792)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2791, &_S2792, k_7);
    float _S2793 = _S2792.differential_0 + _S2314;
    float _S2794 = _S2791.differential_0 + _S2315;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2795;
    (&_S2795)->primal_0 = _S2302;
    (&_S2795)->differential_0 = _S2527;
    s_bwd_length_impl_2(&_S2795, _S2794);
    float2  _S2796 = _S2795.differential_0 + _S2782;
    float3  _S2797 = make_float3 (_S2796.x, _S2796.y, _S2793);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2798;
    (&_S2798)->primal_0 = _S2299;
    (&_S2798)->differential_0 = _S2575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2799;
    (&_S2799)->primal_0 = _S2300;
    (&_S2799)->differential_0 = _S2575;
    s_bwd_prop_mul_2(&_S2798, &_S2799, _S2577.differential_0);
    Matrix<float, 3, 3>  _S2800 = transpose_3(_S2799.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2801;
    (&_S2801)->primal_0 = R_8;
    (&_S2801)->differential_0 = _S2575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2802;
    (&_S2802)->primal_0 = _S2298;
    (&_S2802)->differential_0 = _S2575;
    s_bwd_prop_mul_2(&_S2801, &_S2802, _S2798.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2803;
    (&_S2803)->primal_0 = _S2296;
    (&_S2803)->differential_0 = _S2575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2804;
    (&_S2804)->primal_0 = _S2297;
    (&_S2804)->differential_0 = _S2575;
    s_bwd_prop_mul_2(&_S2803, &_S2804, _S2802.differential_0);
    Matrix<float, 3, 3>  _S2805 = _S2803.differential_0 + transpose_3(_S2804.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2806;
    (&_S2806)->primal_0 = _S2295;
    (&_S2806)->differential_0 = _S2575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2807;
    (&_S2807)->primal_0 = S_2;
    (&_S2807)->differential_0 = _S2575;
    s_bwd_prop_mul_2(&_S2806, &_S2807, _S2805);
    Matrix<float, 3, 3>  _S2808 = transpose_3(_S2806.differential_0);
    float _S2809 = 2.0f * - _S2808.rows[int(2)].z;
    float _S2810 = 2.0f * _S2808.rows[int(2)].y;
    float _S2811 = 2.0f * _S2808.rows[int(2)].x;
    float _S2812 = 2.0f * _S2808.rows[int(1)].z;
    float _S2813 = 2.0f * - _S2808.rows[int(1)].y;
    float _S2814 = 2.0f * _S2808.rows[int(1)].x;
    float _S2815 = 2.0f * _S2808.rows[int(0)].z;
    float _S2816 = 2.0f * _S2808.rows[int(0)].y;
    float _S2817 = 2.0f * - _S2808.rows[int(0)].x;
    float _S2818 = - _S2814 + _S2816;
    float _S2819 = _S2811 + - _S2815;
    float _S2820 = - _S2810 + _S2812;
    float _S2821 = _S2810 + _S2812;
    float _S2822 = _S2811 + _S2815;
    float _S2823 = _S2814 + _S2816;
    float _S2824 = _S2292.w * (_S2813 + _S2817);
    float _S2825 = _S2292.z * (_S2809 + _S2817);
    float _S2826 = _S2292.y * (_S2809 + _S2813);
    float _S2827 = _S2292.x * _S2818 + _S2292.z * _S2821 + _S2292.y * _S2822 + _S2824 + _S2824;
    float _S2828 = _S2292.x * _S2819 + _S2292.w * _S2821 + _S2292.y * _S2823 + _S2825 + _S2825;
    float _S2829 = _S2292.x * _S2820 + _S2292.w * _S2822 + _S2292.z * _S2823 + _S2826 + _S2826;
    float _S2830 = _S2292.w * _S2818 + _S2292.z * _S2819 + _S2292.y * _S2820;
    float3  _S2831 = _S2519;
    *&((&_S2831)->z) = _S2807.differential_0.rows[int(2)].z;
    *&((&_S2831)->y) = _S2807.differential_0.rows[int(1)].y;
    *&((&_S2831)->x) = _S2807.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2832;
    (&_S2832)->primal_0 = scale_8;
    (&_S2832)->differential_0 = _S2519;
    s_bwd_prop_exp_1(&_S2832, _S2831);
    float4  _S2833 = make_float4 (0.0f);
    float4  _S2834 = _S2833;
    *&((&_S2834)->w) = _S2827;
    *&((&_S2834)->z) = _S2828;
    *&((&_S2834)->y) = _S2829;
    *&((&_S2834)->x) = _S2830;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2835;
    (&_S2835)->primal_0 = quat_8;
    (&_S2835)->differential_0 = _S2833;
    s_bwd_normalize_impl_0(&_S2835, _S2834);
    float _S2836 = - (s_diff_k_3 / _S2291);
    DiffPair_float_0 _S2837;
    (&_S2837)->primal_0 = _S2288;
    (&_S2837)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2837, _S2836);
    float _S2838 = - _S2837.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2839;
    (&_S2839)->primal_0 = mean_c_7;
    (&_S2839)->differential_0 = _S2519;
    s_bwd_length_impl_0(&_S2839, 0.0f);
    float3  _S2840 = _S2755 + _S2797 + _S2839.differential_0 + _S2716 + make_float3 (0.0f, 0.0f, v_depth_2);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2841;
    (&_S2841)->primal_0 = R_8;
    (&_S2841)->differential_0 = _S2575;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2842;
    (&_S2842)->primal_0 = mean_9;
    (&_S2842)->differential_0 = _S2519;
    s_bwd_prop_mul_3(&_S2841, &_S2842, _S2840);
    Matrix<float, 3, 3>  _S2843 = _S2800 + _S2801.differential_0 + _S2841.differential_0;
    float _S2844 = _S2838 + _S2526.differential_0;
    float3  _S2845 = _S2832.differential_0 + _S2525.differential_0;
    *v_mean_2 = *v_mean_2 + (_S2842.differential_0 + _S2524.differential_0);
    *v_quat_2 = *v_quat_2 + _S2835.differential_0;
    *v_scale_2 = *v_scale_2 + _S2845;
    *v_in_opacity_2 = *v_in_opacity_2 + _S2844;
    *v_R_2 = *v_R_2 + _S2843;
    *v_t_2 = *v_t_2 + _S2840;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_persp_projection_Intermediates_0
{
    float2  _S2846;
    float2  _S2847;
    float2  _S2848;
    float2  _S2849;
    float2  _S2850;
    float2  _S2851;
    float2  _S2852;
};

inline __device__ void projection_3dgut_persp_vjp(bool antialiased_9, float3  mean_10, float4  quat_9, float3  scale_9, float in_opacity_9, Matrix<float, 3, 3>  R_9, float3  t_9, float fx_12, float fy_12, float cx_10, float cy_10, FixedArray<float, 10>  dist_coeffs_12, uint image_width_9, uint image_height_9, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float2  _S2853 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_persp_projection_Intermediates_0 _S2854;
    (&_S2854)->_S2846 = _S2853;
    (&_S2854)->_S2847 = _S2853;
    (&_S2854)->_S2848 = _S2853;
    (&_S2854)->_S2849 = _S2853;
    (&_S2854)->_S2850 = _S2853;
    (&_S2854)->_S2851 = _S2853;
    (&_S2854)->_S2852 = _S2853;
    float3  _S2855 = make_float3 (0.0f);
    float3  _S2856 = s_primal_ctx_exp_1(scale_9);
    float4  _S2857 = normalize_0(quat_9);
    float _S2858 = _S2857.y;
    float x2_9 = _S2858 * _S2858;
    float y2_9 = _S2857.z * _S2857.z;
    float z2_9 = _S2857.w * _S2857.w;
    float xy_9 = _S2857.y * _S2857.z;
    float xz_9 = _S2857.y * _S2857.w;
    float yz_9 = _S2857.z * _S2857.w;
    float wx_9 = _S2857.x * _S2857.y;
    float wy_9 = _S2857.x * _S2857.z;
    float wz_9 = _S2857.x * _S2857.w;
    Matrix<float, 3, 3>  _S2859 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_9), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_9), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))));
    FixedArray<float3 , 7>  _S2860 = {
        _S2855, _S2855, _S2855, _S2855, _S2855, _S2855, _S2855
    };
    FixedArray<float, 7>  _S2861 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2862;
    (&_S2862)->p_0 = _S2860;
    (&_S2862)->w_mean_0 = _S2861;
    (&_S2862)->w_cov_0 = _S2861;
    (&_S2862)->p_0[int(0)] = mean_10;
    SigmaPoints_0 _S2863 = _S2862;
    (&_S2863)->w_mean_0[int(0)] = 0.0f;
    (&_S2863)->w_cov_0[int(0)] = 2.0f;
    float _S2864 = s_primal_ctx_sqrt_0(3.0f);
    float _S2865 = _S2864 * _S2856.x;
    float3  delta_9 = make_float3 (_S2865) * _S2859.rows[0U];
    float3  _S2866 = mean_10 + delta_9;
    (&_S2863)->p_0[1U] = _S2866;
    float3  _S2867 = mean_10 - delta_9;
    (&_S2863)->p_0[4U] = _S2867;
    float _S2868 = _S2864 * _S2856.y;
    float3  delta_10 = make_float3 (_S2868) * _S2859.rows[1U];
    float3  _S2869 = mean_10 + delta_10;
    (&_S2863)->p_0[2U] = _S2869;
    float3  _S2870 = mean_10 - delta_10;
    (&_S2863)->p_0[5U] = _S2870;
    float _S2871 = _S2864 * _S2856.z;
    float3  delta_11 = make_float3 (_S2871) * _S2859.rows[2U];
    float3  _S2872 = mean_10 + delta_11;
    (&_S2863)->p_0[3U] = _S2872;
    float3  _S2873 = mean_10 - delta_11;
    (&_S2863)->p_0[6U] = _S2873;
    (&_S2863)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S2874 = _S2863;
    (&_S2874)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S2875 = _S2874;
    (&_S2875)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S2876 = _S2875;
    (&_S2876)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S2877 = _S2876;
    (&_S2877)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S2878 = _S2877;
    (&_S2878)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S2879 = _S2878;
    (&_S2879)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S2880 = _S2879;
    (&_S2880)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S2881 = _S2880;
    (&_S2881)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S2882 = _S2881;
    (&_S2882)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S2883 = _S2882;
    (&_S2883)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2884 = _S2883;
    (&_S2884)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2885 = _S2862;
    float3  _S2886 = s_primal_ctx_mul_0(R_9, _S2862.p_0[0U]) + t_9;
    _S2862 = _S2884;
    (&_S2862)->p_0[0U] = _S2886;
    SigmaPoints_0 _S2887 = _S2862;
    (&_S2862)->p_0[1U] = s_primal_ctx_mul_0(R_9, _S2866) + t_9;
    SigmaPoints_0 _S2888 = _S2862;
    (&_S2862)->p_0[2U] = s_primal_ctx_mul_0(R_9, _S2869) + t_9;
    SigmaPoints_0 _S2889 = _S2862;
    (&_S2862)->p_0[3U] = s_primal_ctx_mul_0(R_9, _S2872) + t_9;
    SigmaPoints_0 _S2890 = _S2862;
    (&_S2862)->p_0[4U] = s_primal_ctx_mul_0(R_9, _S2867) + t_9;
    SigmaPoints_0 _S2891 = _S2862;
    (&_S2862)->p_0[5U] = s_primal_ctx_mul_0(R_9, _S2870) + t_9;
    SigmaPoints_0 _S2892 = _S2862;
    (&_S2862)->p_0[6U] = s_primal_ctx_mul_0(R_9, _S2873) + t_9;
    float2  _S2893 = float2 {_S2887.p_0[int(0)].x, _S2887.p_0[int(0)].y} / make_float2 (_S2887.p_0[int(0)].z);
    float u_84 = _S2893.x;
    float v_84 = _S2893.y;
    float r2_84 = u_84 * u_84 + v_84 * v_84;
    float _S2894 = 2.0f * dist_coeffs_12[int(4)];
    float _S2895 = 2.0f * dist_coeffs_12[int(5)];
    float2  _S2896 = _S2893 * make_float2 (1.0f + r2_84 * (dist_coeffs_12[int(0)] + r2_84 * (dist_coeffs_12[int(1)] + r2_84 * (dist_coeffs_12[int(2)] + r2_84 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_84 * v_84 + dist_coeffs_12[int(5)] * (r2_84 + 2.0f * u_84 * u_84) + dist_coeffs_12[int(6)] * r2_84, _S2895 * u_84 * v_84 + dist_coeffs_12[int(4)] * (r2_84 + 2.0f * v_84 * v_84) + dist_coeffs_12[int(7)] * r2_84);
    float2  _S2897 = _S2896 + make_float2 (dist_coeffs_12[int(8)] * _S2896.x + dist_coeffs_12[int(9)] * _S2896.y, 0.0f);
    (&_S2854)->_S2846 = make_float2 (fx_12 * _S2897.x + cx_10, fy_12 * _S2897.y + cy_10);
    float2  _S2898 = float2 {_S2888.p_0[int(1)].x, _S2888.p_0[int(1)].y} / make_float2 (_S2888.p_0[int(1)].z);
    float u_85 = _S2898.x;
    float v_85 = _S2898.y;
    float r2_85 = u_85 * u_85 + v_85 * v_85;
    float2  _S2899 = _S2898 * make_float2 (1.0f + r2_85 * (dist_coeffs_12[int(0)] + r2_85 * (dist_coeffs_12[int(1)] + r2_85 * (dist_coeffs_12[int(2)] + r2_85 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_85 * v_85 + dist_coeffs_12[int(5)] * (r2_85 + 2.0f * u_85 * u_85) + dist_coeffs_12[int(6)] * r2_85, _S2895 * u_85 * v_85 + dist_coeffs_12[int(4)] * (r2_85 + 2.0f * v_85 * v_85) + dist_coeffs_12[int(7)] * r2_85);
    float2  _S2900 = _S2899 + make_float2 (dist_coeffs_12[int(8)] * _S2899.x + dist_coeffs_12[int(9)] * _S2899.y, 0.0f);
    (&_S2854)->_S2847 = make_float2 (fx_12 * _S2900.x + cx_10, fy_12 * _S2900.y + cy_10);
    float2  _S2901 = float2 {_S2889.p_0[int(2)].x, _S2889.p_0[int(2)].y} / make_float2 (_S2889.p_0[int(2)].z);
    float u_86 = _S2901.x;
    float v_86 = _S2901.y;
    float r2_86 = u_86 * u_86 + v_86 * v_86;
    float2  _S2902 = _S2901 * make_float2 (1.0f + r2_86 * (dist_coeffs_12[int(0)] + r2_86 * (dist_coeffs_12[int(1)] + r2_86 * (dist_coeffs_12[int(2)] + r2_86 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_86 * v_86 + dist_coeffs_12[int(5)] * (r2_86 + 2.0f * u_86 * u_86) + dist_coeffs_12[int(6)] * r2_86, _S2895 * u_86 * v_86 + dist_coeffs_12[int(4)] * (r2_86 + 2.0f * v_86 * v_86) + dist_coeffs_12[int(7)] * r2_86);
    float2  _S2903 = _S2902 + make_float2 (dist_coeffs_12[int(8)] * _S2902.x + dist_coeffs_12[int(9)] * _S2902.y, 0.0f);
    (&_S2854)->_S2848 = make_float2 (fx_12 * _S2903.x + cx_10, fy_12 * _S2903.y + cy_10);
    float2  _S2904 = float2 {_S2890.p_0[int(3)].x, _S2890.p_0[int(3)].y} / make_float2 (_S2890.p_0[int(3)].z);
    float u_87 = _S2904.x;
    float v_87 = _S2904.y;
    float r2_87 = u_87 * u_87 + v_87 * v_87;
    float2  _S2905 = _S2904 * make_float2 (1.0f + r2_87 * (dist_coeffs_12[int(0)] + r2_87 * (dist_coeffs_12[int(1)] + r2_87 * (dist_coeffs_12[int(2)] + r2_87 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_87 * v_87 + dist_coeffs_12[int(5)] * (r2_87 + 2.0f * u_87 * u_87) + dist_coeffs_12[int(6)] * r2_87, _S2895 * u_87 * v_87 + dist_coeffs_12[int(4)] * (r2_87 + 2.0f * v_87 * v_87) + dist_coeffs_12[int(7)] * r2_87);
    float2  _S2906 = _S2905 + make_float2 (dist_coeffs_12[int(8)] * _S2905.x + dist_coeffs_12[int(9)] * _S2905.y, 0.0f);
    (&_S2854)->_S2849 = make_float2 (fx_12 * _S2906.x + cx_10, fy_12 * _S2906.y + cy_10);
    float2  _S2907 = float2 {_S2891.p_0[int(4)].x, _S2891.p_0[int(4)].y} / make_float2 (_S2891.p_0[int(4)].z);
    float u_88 = _S2907.x;
    float v_88 = _S2907.y;
    float r2_88 = u_88 * u_88 + v_88 * v_88;
    float2  _S2908 = _S2907 * make_float2 (1.0f + r2_88 * (dist_coeffs_12[int(0)] + r2_88 * (dist_coeffs_12[int(1)] + r2_88 * (dist_coeffs_12[int(2)] + r2_88 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_88 * v_88 + dist_coeffs_12[int(5)] * (r2_88 + 2.0f * u_88 * u_88) + dist_coeffs_12[int(6)] * r2_88, _S2895 * u_88 * v_88 + dist_coeffs_12[int(4)] * (r2_88 + 2.0f * v_88 * v_88) + dist_coeffs_12[int(7)] * r2_88);
    float2  _S2909 = _S2908 + make_float2 (dist_coeffs_12[int(8)] * _S2908.x + dist_coeffs_12[int(9)] * _S2908.y, 0.0f);
    (&_S2854)->_S2850 = make_float2 (fx_12 * _S2909.x + cx_10, fy_12 * _S2909.y + cy_10);
    float2  _S2910 = float2 {_S2892.p_0[int(5)].x, _S2892.p_0[int(5)].y} / make_float2 (_S2892.p_0[int(5)].z);
    float u_89 = _S2910.x;
    float v_89 = _S2910.y;
    float r2_89 = u_89 * u_89 + v_89 * v_89;
    float2  _S2911 = _S2910 * make_float2 (1.0f + r2_89 * (dist_coeffs_12[int(0)] + r2_89 * (dist_coeffs_12[int(1)] + r2_89 * (dist_coeffs_12[int(2)] + r2_89 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_89 * v_89 + dist_coeffs_12[int(5)] * (r2_89 + 2.0f * u_89 * u_89) + dist_coeffs_12[int(6)] * r2_89, _S2895 * u_89 * v_89 + dist_coeffs_12[int(4)] * (r2_89 + 2.0f * v_89 * v_89) + dist_coeffs_12[int(7)] * r2_89);
    float2  _S2912 = _S2911 + make_float2 (dist_coeffs_12[int(8)] * _S2911.x + dist_coeffs_12[int(9)] * _S2911.y, 0.0f);
    (&_S2854)->_S2851 = make_float2 (fx_12 * _S2912.x + cx_10, fy_12 * _S2912.y + cy_10);
    float2  _S2913 = float2 {_S2862.p_0[int(6)].x, _S2862.p_0[int(6)].y} / make_float2 (_S2862.p_0[int(6)].z);
    float u_90 = _S2913.x;
    float v_90 = _S2913.y;
    float r2_90 = u_90 * u_90 + v_90 * v_90;
    float2  _S2914 = _S2913 * make_float2 (1.0f + r2_90 * (dist_coeffs_12[int(0)] + r2_90 * (dist_coeffs_12[int(1)] + r2_90 * (dist_coeffs_12[int(2)] + r2_90 * dist_coeffs_12[int(3)])))) + make_float2 (_S2894 * u_90 * v_90 + dist_coeffs_12[int(5)] * (r2_90 + 2.0f * u_90 * u_90) + dist_coeffs_12[int(6)] * r2_90, _S2895 * u_90 * v_90 + dist_coeffs_12[int(4)] * (r2_90 + 2.0f * v_90 * v_90) + dist_coeffs_12[int(7)] * r2_90);
    float2  _S2915 = _S2914 + make_float2 (dist_coeffs_12[int(8)] * _S2914.x + dist_coeffs_12[int(9)] * _S2914.y, 0.0f);
    (&_S2854)->_S2852 = make_float2 (fx_12 * _S2915.x + cx_10, fy_12 * _S2915.y + cy_10);
    float _S2916 = - in_opacity_9;
    float _S2917 = 1.0f + s_primal_ctx_exp_0(_S2916);
    float _S2918 = 1.0f / _S2917;
    float _S2919 = _S2917 * _S2917;
    float3  _S2920 = make_float3 (_S2865);
    float3  _S2921 = make_float3 (_S2868);
    float3  _S2922 = make_float3 (_S2871);
    float _S2923 = float(image_width_9);
    float _S2924 = float(image_height_9);
    float _S2925 = 0.30000001192092896f * (0.5f * _S2923 / fx_12) * fx_12;
    float lim_x_pos_4 = _S2923 + _S2925;
    float _S2926 = 0.30000001192092896f * (0.5f * _S2924 / fy_12) * fy_12;
    float lim_y_pos_1 = _S2924 + _S2926;
    float2  _S2927 = make_float2 (_S2863.w_mean_0[int(1)]) * _S2854._S2847 + make_float2 (_S2875.w_mean_0[int(2)]) * _S2854._S2848 + make_float2 (_S2877.w_mean_0[int(3)]) * _S2854._S2849 + make_float2 (_S2879.w_mean_0[int(4)]) * _S2854._S2850 + make_float2 (_S2881.w_mean_0[int(5)]) * _S2854._S2851 + make_float2 (_S2883.w_mean_0[int(6)]) * _S2854._S2852;
    float _S2928 = - _S2925;
    float _S2929 = - _S2926;
    float2  _S2930 = make_float2 (s_primal_ctx_clamp_0(_S2927.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2927.y, _S2929, lim_y_pos_1));
    float2  d_21 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2846.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2846.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2931 = d_21.x;
    float _S2932 = d_21.y;
    float _S2933 = _S2931 * _S2932;
    float2  d_22 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2847.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2847.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2934 = d_22.x;
    float _S2935 = d_22.y;
    float _S2936 = _S2934 * _S2935;
    float2  d_23 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2848.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2848.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2937 = d_23.x;
    float _S2938 = d_23.y;
    float _S2939 = _S2937 * _S2938;
    float2  d_24 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2849.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2849.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2940 = d_24.x;
    float _S2941 = d_24.y;
    float _S2942 = _S2940 * _S2941;
    float2  d_25 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2850.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2850.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2943 = d_25.x;
    float _S2944 = d_25.y;
    float _S2945 = _S2943 * _S2944;
    float2  d_26 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2851.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2851.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2946 = d_26.x;
    float _S2947 = d_26.y;
    float _S2948 = _S2946 * _S2947;
    float2  d_27 = make_float2 (s_primal_ctx_clamp_0(_S2854._S2852.x, _S2928, lim_x_pos_4), s_primal_ctx_clamp_0(_S2854._S2852.y, _S2929, lim_y_pos_1)) - _S2930;
    float _S2949 = d_27.x;
    float _S2950 = d_27.y;
    float _S2951 = _S2949 * _S2950;
    Matrix<float, 2, 2>  covar2d_6 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S2931 * _S2931, _S2933, _S2933, _S2932 * _S2932) + makeMatrix<float, 2, 2> (_S2874.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S2934 * _S2934, _S2936, _S2936, _S2935 * _S2935) + makeMatrix<float, 2, 2> (_S2876.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S2937 * _S2937, _S2939, _S2939, _S2938 * _S2938) + makeMatrix<float, 2, 2> (_S2878.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S2940 * _S2940, _S2942, _S2942, _S2941 * _S2941) + makeMatrix<float, 2, 2> (_S2880.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S2943 * _S2943, _S2945, _S2945, _S2944 * _S2944) + makeMatrix<float, 2, 2> (_S2882.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S2946 * _S2946, _S2948, _S2948, _S2947 * _S2947) + makeMatrix<float, 2, 2> (_S2884.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S2949 * _S2949, _S2951, _S2951, _S2950 * _S2950);
    float eps2d_9;
    if(antialiased_9)
    {
        eps2d_9 = 0.10000000149011612f;
    }
    else
    {
        eps2d_9 = 0.30000001192092896f;
    }
    float _S2952 = covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x;
    float det_orig_9 = covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - _S2952;
    float _S2953 = covar2d_6.rows[int(0)].x + eps2d_9;
    Matrix<float, 2, 2>  _S2954 = covar2d_6;
    *&(((&_S2954)->rows + (int(0)))->x) = _S2953;
    float _S2955 = covar2d_6.rows[int(1)].y + eps2d_9;
    *&(((&_S2954)->rows + (int(1)))->y) = _S2955;
    Matrix<float, 2, 2>  _S2956 = _S2954;
    Matrix<float, 2, 2>  _S2957 = _S2954;
    float det_blur_9 = _S2953 * _S2955 - _S2952;
    float _S2958 = det_orig_9 / det_blur_9;
    float _S2959 = det_blur_9 * det_blur_9;
    float _S2960 = (F32_max((0.0f), (_S2958)));
    float _S2961 = s_primal_ctx_sqrt_0(_S2960);
    float invdet_11 = 1.0f / det_blur_9;
    float _S2962 = - covar2d_6.rows[int(0)].y;
    float _S2963 = - covar2d_6.rows[int(1)].x;
    if(antialiased_9)
    {
        eps2d_9 = _S2918 * _S2961;
    }
    else
    {
        eps2d_9 = _S2918;
    }
    float _S2964 = eps2d_9 / 0.00392156885936856f;
    float _S2965 = 2.0f * s_primal_ctx_log_0(_S2964);
    float _S2966 = s_primal_ctx_sqrt_0(_S2965);
    float _S2967 = _S2956.rows[int(0)].x;
    float _S2968 = _S2957.rows[int(1)].y;
    float3  campos_4 = - s_primal_ctx_mul_0(transpose_3(R_9), t_9);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2969;
    (&_S2969)->primal_0 = mean_10;
    (&_S2969)->differential_0 = _S2855;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2970;
    (&_S2970)->primal_0 = scale_9;
    (&_S2970)->differential_0 = _S2855;
    DiffPair_float_0 _S2971;
    (&_S2971)->primal_0 = in_opacity_9;
    (&_S2971)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2972;
    (&_S2972)->primal_0 = campos_4;
    (&_S2972)->differential_0 = _S2855;
    s_bwd_prop_view_radius_3dgs_0(&_S2969, &_S2970, &_S2971, &_S2972, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2973 = _S2969;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2974 = _S2970;
    DiffPair_float_0 _S2975 = _S2971;
    float2  _S2976 = _S2853;
    *&((&_S2976)->y) = v_conic_3.z;
    float2  _S2977 = _S2853;
    *&((&_S2977)->y) = v_conic_3.y;
    *&((&_S2977)->x) = v_conic_3.x;
    DiffPair_float_0 _S2978;
    (&_S2978)->primal_0 = _S2968;
    (&_S2978)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2978, 0.0f);
    DiffPair_float_0 _S2979;
    (&_S2979)->primal_0 = _S2967;
    (&_S2979)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2979, 0.0f);
    DiffPair_float_0 _S2980;
    (&_S2980)->primal_0 = 3.32999992370605469f;
    (&_S2980)->differential_0 = 0.0f;
    DiffPair_float_0 _S2981;
    (&_S2981)->primal_0 = _S2966;
    (&_S2981)->differential_0 = 0.0f;
    _d_min_0(&_S2980, &_S2981, 0.0f);
    DiffPair_float_0 _S2982;
    (&_S2982)->primal_0 = _S2965;
    (&_S2982)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2982, _S2981.differential_0);
    float _S2983 = 2.0f * _S2982.differential_0;
    DiffPair_float_0 _S2984;
    (&_S2984)->primal_0 = _S2964;
    (&_S2984)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2984, _S2983);
    float _S2985 = v_opacity_3 + 254.9999847412109375f * _S2984.differential_0;
    Matrix<float, 2, 2>  _S2986 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2987 = _S2986;
    _S2987[int(1)] = _S2976;
    _S2987[int(0)] = _S2977;
    Matrix<float, 2, 2>  _S2988 = _S2987;
    float2  _S2989 = make_float2 (0.0f, _S2978.differential_0);
    float2  _S2990 = make_float2 (_S2979.differential_0, 0.0f);
    float _S2991;
    if(antialiased_9)
    {
        float _S2992 = _S2961 * _S2985;
        eps2d_9 = _S2918 * _S2985;
        _S2991 = _S2992;
    }
    else
    {
        eps2d_9 = 0.0f;
        _S2991 = _S2985;
    }
    float _S2993 = invdet_11 * _S2988.rows[int(1)].y;
    float _S2994 = - (invdet_11 * _S2988.rows[int(1)].x);
    float _S2995 = - (invdet_11 * _S2988.rows[int(0)].y);
    float _S2996 = invdet_11 * _S2988.rows[int(0)].x;
    float _S2997 = - ((_S2953 * _S2988.rows[int(1)].y + _S2963 * _S2988.rows[int(1)].x + _S2962 * _S2988.rows[int(0)].y + _S2955 * _S2988.rows[int(0)].x) / _S2959);
    DiffPair_float_0 _S2998;
    (&_S2998)->primal_0 = _S2960;
    (&_S2998)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2998, eps2d_9);
    DiffPair_float_0 _S2999;
    (&_S2999)->primal_0 = 0.0f;
    (&_S2999)->differential_0 = 0.0f;
    DiffPair_float_0 _S3000;
    (&_S3000)->primal_0 = _S2958;
    (&_S3000)->differential_0 = 0.0f;
    _d_max_0(&_S2999, &_S3000, _S2998.differential_0);
    float _S3001 = _S3000.differential_0 / _S2959;
    float s_diff_det_orig_T_3 = det_blur_9 * _S3001;
    float _S3002 = det_orig_9 * - _S3001 + _S2997;
    float _S3003 = - _S3002;
    float _S3004 = _S2953 * _S3002;
    float _S3005 = _S2955 * _S3002;
    Matrix<float, 2, 2>  _S3006 = _S2986;
    _S3006[int(1)] = _S2989;
    _S3006[int(0)] = _S2990;
    float _S3007 = _S3005 + _S3006.rows[int(0)].x + _S2993;
    float _S3008 = _S3003 + - s_diff_det_orig_T_3;
    float _S3009 = covar2d_6.rows[int(0)].y * _S3008 + _S2994;
    float _S3010 = covar2d_6.rows[int(1)].x * _S3008 + _S2995;
    float _S3011 = covar2d_6.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S3012 = _S3004 + _S3006.rows[int(1)].y + _S2996 + covar2d_6.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S3013 = _S2853;
    *&((&_S3013)->x) = _S3009;
    *&((&_S3013)->y) = _S3012;
    float _S3014 = _S3007 + _S3011;
    float2  _S3015 = _S2853;
    *&((&_S3015)->y) = _S3010;
    *&((&_S3015)->x) = _S3014;
    Matrix<float, 2, 2>  _S3016 = _S2986;
    _S3016[int(1)] = _S3013;
    _S3016[int(0)] = _S3015;
    Matrix<float, 3, 3>  _S3017 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3018;
    (&_S3018)->primal_0 = R_9;
    (&_S3018)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3019;
    (&_S3019)->primal_0 = _S2873;
    (&_S3019)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3018, &_S3019, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3020;
    (&_S3020)->primal_0 = R_9;
    (&_S3020)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3021;
    (&_S3021)->primal_0 = _S2870;
    (&_S3021)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3020, &_S3021, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3022;
    (&_S3022)->primal_0 = R_9;
    (&_S3022)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3023;
    (&_S3023)->primal_0 = _S2867;
    (&_S3023)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3022, &_S3023, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3024;
    (&_S3024)->primal_0 = R_9;
    (&_S3024)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3025;
    (&_S3025)->primal_0 = _S2872;
    (&_S3025)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3024, &_S3025, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3026;
    (&_S3026)->primal_0 = R_9;
    (&_S3026)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3027;
    (&_S3027)->primal_0 = _S2869;
    (&_S3027)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3026, &_S3027, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3028;
    (&_S3028)->primal_0 = R_9;
    (&_S3028)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3029;
    (&_S3029)->primal_0 = _S2866;
    (&_S3029)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3028, &_S3029, _S2855);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3030;
    (&_S3030)->primal_0 = R_9;
    (&_S3030)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3031;
    (&_S3031)->primal_0 = _S2885.p_0[0U];
    (&_S3031)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3030, &_S3031, _S2855);
    float3  _S3032 = - _S3019.differential_0 + _S3025.differential_0;
    float3  _S3033 = _S2922 * _S3032;
    float3  _S3034 = _S2859.rows[2U] * _S3032;
    float _S3035 = _S2864 * (_S3034.x + _S3034.y + _S3034.z);
    float3  _S3036 = - _S3021.differential_0 + _S3027.differential_0;
    float3  _S3037 = _S2921 * _S3036;
    float3  _S3038 = _S2859.rows[1U] * _S3036;
    float _S3039 = _S2864 * (_S3038.x + _S3038.y + _S3038.z);
    float3  _S3040 = - _S3023.differential_0 + _S3029.differential_0;
    float3  _S3041 = _S2920 * _S3040;
    float3  _S3042 = _S2859.rows[0U] * _S3040;
    float _S3043 = _S2864 * (_S3042.x + _S3042.y + _S3042.z);
    Matrix<float, 3, 3>  _S3044 = _S3017;
    _S3044[2U] = _S3033;
    _S3044[1U] = _S3037;
    _S3044[0U] = _S3041;
    Matrix<float, 3, 3>  _S3045 = transpose_3(transpose_3(_S3044));
    float _S3046 = 2.0f * - _S3045.rows[int(2)].z;
    float _S3047 = 2.0f * _S3045.rows[int(2)].y;
    float _S3048 = 2.0f * _S3045.rows[int(2)].x;
    float _S3049 = 2.0f * _S3045.rows[int(1)].z;
    float _S3050 = 2.0f * - _S3045.rows[int(1)].y;
    float _S3051 = 2.0f * _S3045.rows[int(1)].x;
    float _S3052 = 2.0f * _S3045.rows[int(0)].z;
    float _S3053 = 2.0f * _S3045.rows[int(0)].y;
    float _S3054 = 2.0f * - _S3045.rows[int(0)].x;
    float _S3055 = - _S3051 + _S3053;
    float _S3056 = _S3048 + - _S3052;
    float _S3057 = - _S3047 + _S3049;
    float _S3058 = _S3047 + _S3049;
    float _S3059 = _S3048 + _S3052;
    float _S3060 = _S3051 + _S3053;
    float _S3061 = _S2857.w * (_S3050 + _S3054);
    float _S3062 = _S2857.z * (_S3046 + _S3054);
    float _S3063 = _S2857.y * (_S3046 + _S3050);
    float _S3064 = _S2857.x * _S3055 + _S2857.z * _S3058 + _S2857.y * _S3059 + _S3061 + _S3061;
    float _S3065 = _S2857.x * _S3056 + _S2857.w * _S3058 + _S2857.y * _S3060 + _S3062 + _S3062;
    float _S3066 = _S2857.x * _S3057 + _S2857.w * _S3059 + _S2857.z * _S3060 + _S3063 + _S3063;
    float _S3067 = _S2857.w * _S3055 + _S2857.z * _S3056 + _S2857.y * _S3057;
    float4  _S3068 = make_float4 (0.0f);
    float4  _S3069 = _S3068;
    *&((&_S3069)->w) = _S3064;
    *&((&_S3069)->z) = _S3065;
    *&((&_S3069)->y) = _S3066;
    *&((&_S3069)->x) = _S3067;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3070;
    (&_S3070)->primal_0 = quat_9;
    (&_S3070)->differential_0 = _S3068;
    s_bwd_normalize_impl_0(&_S3070, _S3069);
    float3  _S3071 = _S2855;
    *&((&_S3071)->z) = _S3035;
    *&((&_S3071)->y) = _S3039;
    *&((&_S3071)->x) = _S3043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3072;
    (&_S3072)->primal_0 = scale_9;
    (&_S3072)->differential_0 = _S2855;
    s_bwd_prop_exp_1(&_S3072, _S3071);
    float _S3073 = - (_S2991 / _S2919);
    DiffPair_float_0 _S3074;
    (&_S3074)->primal_0 = _S2916;
    (&_S3074)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3074, _S3073);
    float _S3075 = - _S3074.differential_0;
    float3  s_diff_mean_c_T_0 = make_float3 (0.0f, 0.0f, v_depth_3);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3076;
    (&_S3076)->primal_0 = R_9;
    (&_S3076)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3077;
    (&_S3077)->primal_0 = mean_10;
    (&_S3077)->differential_0 = _S2855;
    s_bwd_prop_mul_3(&_S3076, &_S3077, s_diff_mean_c_T_0);
    Matrix<float, 3, 3>  _S3078 = _S3018.differential_0 + _S3020.differential_0 + _S3022.differential_0 + _S3024.differential_0 + _S3026.differential_0 + _S3028.differential_0 + _S3030.differential_0 + _S3076.differential_0;
    float _S3079 = _S3075 + _S2975.differential_0;
    float3  _S3080 = _S3072.differential_0 + _S2974.differential_0;
    *v_mean_3 = *v_mean_3 + (_S3019.differential_0 + _S3025.differential_0 + _S3021.differential_0 + _S3027.differential_0 + _S3023.differential_0 + _S3029.differential_0 + _S3077.differential_0 + _S2973.differential_0);
    *v_quat_3 = *v_quat_3 + _S3070.differential_0;
    *v_scale_3 = *v_scale_3 + _S3080;
    *v_in_opacity_3 = *v_in_opacity_3 + _S3079;
    *v_R_3 = *v_R_3 + _S3078;
    *v_t_3 = *v_t_3 + s_diff_mean_c_T_0;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_fisheye_projection_Intermediates_0
{
    float2  _S3081;
    float2  _S3082;
    float2  _S3083;
    float2  _S3084;
    float2  _S3085;
    float2  _S3086;
    float2  _S3087;
};

inline __device__ void projection_3dgut_fisheye_vjp(bool antialiased_10, float3  mean_11, float4  quat_10, float3  scale_10, float in_opacity_10, Matrix<float, 3, 3>  R_10, float3  t_10, float fx_13, float fy_13, float cx_11, float cy_11, FixedArray<float, 10>  dist_coeffs_13, uint image_width_10, uint image_height_10, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float2  _S3088 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_fisheye_projection_Intermediates_0 _S3089;
    (&_S3089)->_S3081 = _S3088;
    (&_S3089)->_S3082 = _S3088;
    (&_S3089)->_S3083 = _S3088;
    (&_S3089)->_S3084 = _S3088;
    (&_S3089)->_S3085 = _S3088;
    (&_S3089)->_S3086 = _S3088;
    (&_S3089)->_S3087 = _S3088;
    (&_S3089)->_S3081 = _S3088;
    (&_S3089)->_S3082 = _S3088;
    (&_S3089)->_S3083 = _S3088;
    (&_S3089)->_S3084 = _S3088;
    (&_S3089)->_S3085 = _S3088;
    (&_S3089)->_S3086 = _S3088;
    (&_S3089)->_S3087 = _S3088;
    float3  _S3090 = make_float3 (0.0f);
    float3  _S3091 = s_primal_ctx_exp_1(scale_10);
    float4  _S3092 = normalize_0(quat_10);
    float _S3093 = _S3092.y;
    float x2_10 = _S3093 * _S3093;
    float y2_10 = _S3092.z * _S3092.z;
    float z2_10 = _S3092.w * _S3092.w;
    float xy_10 = _S3092.y * _S3092.z;
    float xz_10 = _S3092.y * _S3092.w;
    float yz_10 = _S3092.z * _S3092.w;
    float wx_10 = _S3092.x * _S3092.y;
    float wy_10 = _S3092.x * _S3092.z;
    float wz_10 = _S3092.x * _S3092.w;
    Matrix<float, 3, 3>  _S3094 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_10), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_10), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))));
    FixedArray<float3 , 7>  _S3095 = {
        _S3090, _S3090, _S3090, _S3090, _S3090, _S3090, _S3090
    };
    FixedArray<float, 7>  _S3096 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3097;
    (&_S3097)->p_0 = _S3095;
    (&_S3097)->w_mean_0 = _S3096;
    (&_S3097)->w_cov_0 = _S3096;
    (&_S3097)->p_0[int(0)] = mean_11;
    SigmaPoints_0 _S3098 = _S3097;
    (&_S3098)->w_mean_0[int(0)] = 0.0f;
    (&_S3098)->w_cov_0[int(0)] = 2.0f;
    float _S3099 = s_primal_ctx_sqrt_0(3.0f);
    float _S3100 = _S3099 * _S3091.x;
    float3  delta_12 = make_float3 (_S3100) * _S3094.rows[0U];
    float3  _S3101 = mean_11 + delta_12;
    (&_S3098)->p_0[1U] = _S3101;
    float3  _S3102 = mean_11 - delta_12;
    (&_S3098)->p_0[4U] = _S3102;
    float _S3103 = _S3099 * _S3091.y;
    float3  delta_13 = make_float3 (_S3103) * _S3094.rows[1U];
    float3  _S3104 = mean_11 + delta_13;
    (&_S3098)->p_0[2U] = _S3104;
    float3  _S3105 = mean_11 - delta_13;
    (&_S3098)->p_0[5U] = _S3105;
    float _S3106 = _S3099 * _S3091.z;
    float3  delta_14 = make_float3 (_S3106) * _S3094.rows[2U];
    float3  _S3107 = mean_11 + delta_14;
    (&_S3098)->p_0[3U] = _S3107;
    float3  _S3108 = mean_11 - delta_14;
    (&_S3098)->p_0[6U] = _S3108;
    (&_S3098)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3109 = _S3098;
    (&_S3109)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3110 = _S3109;
    (&_S3110)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3111 = _S3110;
    (&_S3111)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3112 = _S3111;
    (&_S3112)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3113 = _S3112;
    (&_S3113)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3114 = _S3113;
    (&_S3114)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3115 = _S3114;
    (&_S3115)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3116 = _S3115;
    (&_S3116)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3117 = _S3116;
    (&_S3117)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3118 = _S3117;
    (&_S3118)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3119 = _S3118;
    (&_S3119)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3120 = _S3097;
    float3  _S3121 = s_primal_ctx_mul_0(R_10, _S3097.p_0[0U]) + t_10;
    _S3097 = _S3119;
    (&_S3097)->p_0[0U] = _S3121;
    SigmaPoints_0 _S3122 = _S3097;
    (&_S3097)->p_0[1U] = s_primal_ctx_mul_0(R_10, _S3101) + t_10;
    SigmaPoints_0 _S3123 = _S3097;
    (&_S3097)->p_0[2U] = s_primal_ctx_mul_0(R_10, _S3104) + t_10;
    SigmaPoints_0 _S3124 = _S3097;
    (&_S3097)->p_0[3U] = s_primal_ctx_mul_0(R_10, _S3107) + t_10;
    SigmaPoints_0 _S3125 = _S3097;
    (&_S3097)->p_0[4U] = s_primal_ctx_mul_0(R_10, _S3102) + t_10;
    SigmaPoints_0 _S3126 = _S3097;
    (&_S3097)->p_0[5U] = s_primal_ctx_mul_0(R_10, _S3105) + t_10;
    SigmaPoints_0 _S3127 = _S3097;
    (&_S3097)->p_0[6U] = s_primal_ctx_mul_0(R_10, _S3108) + t_10;
    SigmaPoints_0 _S3128 = _S3097;
    float2  _S3129 = float2 {_S3122.p_0[int(0)].x, _S3122.p_0[int(0)].y};
    float _S3130 = length_0(_S3129);
    float _S3131 = _S3122.p_0[int(0)].z;
    float _S3132 = s_primal_ctx_atan2_0(_S3130, _S3131);
    float k_8;
    if(_S3132 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3132 * _S3132 / 3.0f) / _S3131;
    }
    else
    {
        k_8 = _S3132 / _S3130;
    }
    float2  _S3133 = _S3129 * make_float2 (k_8);
    float u_91 = _S3133.x;
    float v_91 = _S3133.y;
    float r2_91 = u_91 * u_91 + v_91 * v_91;
    float _S3134 = 2.0f * dist_coeffs_13[int(4)];
    float _S3135 = 2.0f * dist_coeffs_13[int(5)];
    float2  _S3136 = _S3133 * make_float2 (1.0f + r2_91 * (dist_coeffs_13[int(0)] + r2_91 * (dist_coeffs_13[int(1)] + r2_91 * (dist_coeffs_13[int(2)] + r2_91 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_91 * v_91 + dist_coeffs_13[int(5)] * (r2_91 + 2.0f * u_91 * u_91) + dist_coeffs_13[int(6)] * r2_91, _S3135 * u_91 * v_91 + dist_coeffs_13[int(4)] * (r2_91 + 2.0f * v_91 * v_91) + dist_coeffs_13[int(7)] * r2_91);
    float2  _S3137 = _S3136 + make_float2 (dist_coeffs_13[int(8)] * _S3136.x + dist_coeffs_13[int(9)] * _S3136.y, 0.0f);
    (&_S3089)->_S3081 = make_float2 (fx_13 * _S3137.x + cx_11, fy_13 * _S3137.y + cy_11);
    float2  _S3138 = float2 {_S3123.p_0[int(1)].x, _S3123.p_0[int(1)].y};
    float _S3139 = length_0(_S3138);
    float _S3140 = _S3123.p_0[int(1)].z;
    float _S3141 = s_primal_ctx_atan2_0(_S3139, _S3140);
    if(_S3141 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3141 * _S3141 / 3.0f) / _S3140;
    }
    else
    {
        k_8 = _S3141 / _S3139;
    }
    float2  _S3142 = _S3138 * make_float2 (k_8);
    float u_92 = _S3142.x;
    float v_92 = _S3142.y;
    float r2_92 = u_92 * u_92 + v_92 * v_92;
    float2  _S3143 = _S3142 * make_float2 (1.0f + r2_92 * (dist_coeffs_13[int(0)] + r2_92 * (dist_coeffs_13[int(1)] + r2_92 * (dist_coeffs_13[int(2)] + r2_92 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_92 * v_92 + dist_coeffs_13[int(5)] * (r2_92 + 2.0f * u_92 * u_92) + dist_coeffs_13[int(6)] * r2_92, _S3135 * u_92 * v_92 + dist_coeffs_13[int(4)] * (r2_92 + 2.0f * v_92 * v_92) + dist_coeffs_13[int(7)] * r2_92);
    float2  _S3144 = _S3143 + make_float2 (dist_coeffs_13[int(8)] * _S3143.x + dist_coeffs_13[int(9)] * _S3143.y, 0.0f);
    (&_S3089)->_S3082 = make_float2 (fx_13 * _S3144.x + cx_11, fy_13 * _S3144.y + cy_11);
    float2  _S3145 = float2 {_S3124.p_0[int(2)].x, _S3124.p_0[int(2)].y};
    float _S3146 = length_0(_S3145);
    float _S3147 = _S3124.p_0[int(2)].z;
    float _S3148 = s_primal_ctx_atan2_0(_S3146, _S3147);
    if(_S3148 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3148 * _S3148 / 3.0f) / _S3147;
    }
    else
    {
        k_8 = _S3148 / _S3146;
    }
    float2  _S3149 = _S3145 * make_float2 (k_8);
    float u_93 = _S3149.x;
    float v_93 = _S3149.y;
    float r2_93 = u_93 * u_93 + v_93 * v_93;
    float2  _S3150 = _S3149 * make_float2 (1.0f + r2_93 * (dist_coeffs_13[int(0)] + r2_93 * (dist_coeffs_13[int(1)] + r2_93 * (dist_coeffs_13[int(2)] + r2_93 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_93 * v_93 + dist_coeffs_13[int(5)] * (r2_93 + 2.0f * u_93 * u_93) + dist_coeffs_13[int(6)] * r2_93, _S3135 * u_93 * v_93 + dist_coeffs_13[int(4)] * (r2_93 + 2.0f * v_93 * v_93) + dist_coeffs_13[int(7)] * r2_93);
    float2  _S3151 = _S3150 + make_float2 (dist_coeffs_13[int(8)] * _S3150.x + dist_coeffs_13[int(9)] * _S3150.y, 0.0f);
    (&_S3089)->_S3083 = make_float2 (fx_13 * _S3151.x + cx_11, fy_13 * _S3151.y + cy_11);
    float2  _S3152 = float2 {_S3125.p_0[int(3)].x, _S3125.p_0[int(3)].y};
    float _S3153 = length_0(_S3152);
    float _S3154 = _S3125.p_0[int(3)].z;
    float _S3155 = s_primal_ctx_atan2_0(_S3153, _S3154);
    if(_S3155 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3155 * _S3155 / 3.0f) / _S3154;
    }
    else
    {
        k_8 = _S3155 / _S3153;
    }
    float2  _S3156 = _S3152 * make_float2 (k_8);
    float u_94 = _S3156.x;
    float v_94 = _S3156.y;
    float r2_94 = u_94 * u_94 + v_94 * v_94;
    float2  _S3157 = _S3156 * make_float2 (1.0f + r2_94 * (dist_coeffs_13[int(0)] + r2_94 * (dist_coeffs_13[int(1)] + r2_94 * (dist_coeffs_13[int(2)] + r2_94 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_94 * v_94 + dist_coeffs_13[int(5)] * (r2_94 + 2.0f * u_94 * u_94) + dist_coeffs_13[int(6)] * r2_94, _S3135 * u_94 * v_94 + dist_coeffs_13[int(4)] * (r2_94 + 2.0f * v_94 * v_94) + dist_coeffs_13[int(7)] * r2_94);
    float2  _S3158 = _S3157 + make_float2 (dist_coeffs_13[int(8)] * _S3157.x + dist_coeffs_13[int(9)] * _S3157.y, 0.0f);
    (&_S3089)->_S3084 = make_float2 (fx_13 * _S3158.x + cx_11, fy_13 * _S3158.y + cy_11);
    float2  _S3159 = float2 {_S3126.p_0[int(4)].x, _S3126.p_0[int(4)].y};
    float _S3160 = length_0(_S3159);
    float _S3161 = _S3126.p_0[int(4)].z;
    float _S3162 = s_primal_ctx_atan2_0(_S3160, _S3161);
    if(_S3162 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3162 * _S3162 / 3.0f) / _S3161;
    }
    else
    {
        k_8 = _S3162 / _S3160;
    }
    float2  _S3163 = _S3159 * make_float2 (k_8);
    float u_95 = _S3163.x;
    float v_95 = _S3163.y;
    float r2_95 = u_95 * u_95 + v_95 * v_95;
    float2  _S3164 = _S3163 * make_float2 (1.0f + r2_95 * (dist_coeffs_13[int(0)] + r2_95 * (dist_coeffs_13[int(1)] + r2_95 * (dist_coeffs_13[int(2)] + r2_95 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_95 * v_95 + dist_coeffs_13[int(5)] * (r2_95 + 2.0f * u_95 * u_95) + dist_coeffs_13[int(6)] * r2_95, _S3135 * u_95 * v_95 + dist_coeffs_13[int(4)] * (r2_95 + 2.0f * v_95 * v_95) + dist_coeffs_13[int(7)] * r2_95);
    float2  _S3165 = _S3164 + make_float2 (dist_coeffs_13[int(8)] * _S3164.x + dist_coeffs_13[int(9)] * _S3164.y, 0.0f);
    (&_S3089)->_S3085 = make_float2 (fx_13 * _S3165.x + cx_11, fy_13 * _S3165.y + cy_11);
    float2  _S3166 = float2 {_S3127.p_0[int(5)].x, _S3127.p_0[int(5)].y};
    float _S3167 = length_0(_S3166);
    float _S3168 = _S3127.p_0[int(5)].z;
    float _S3169 = s_primal_ctx_atan2_0(_S3167, _S3168);
    if(_S3169 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3169 * _S3169 / 3.0f) / _S3168;
    }
    else
    {
        k_8 = _S3169 / _S3167;
    }
    float2  _S3170 = _S3166 * make_float2 (k_8);
    float u_96 = _S3170.x;
    float v_96 = _S3170.y;
    float r2_96 = u_96 * u_96 + v_96 * v_96;
    float2  _S3171 = _S3170 * make_float2 (1.0f + r2_96 * (dist_coeffs_13[int(0)] + r2_96 * (dist_coeffs_13[int(1)] + r2_96 * (dist_coeffs_13[int(2)] + r2_96 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_96 * v_96 + dist_coeffs_13[int(5)] * (r2_96 + 2.0f * u_96 * u_96) + dist_coeffs_13[int(6)] * r2_96, _S3135 * u_96 * v_96 + dist_coeffs_13[int(4)] * (r2_96 + 2.0f * v_96 * v_96) + dist_coeffs_13[int(7)] * r2_96);
    float2  _S3172 = _S3171 + make_float2 (dist_coeffs_13[int(8)] * _S3171.x + dist_coeffs_13[int(9)] * _S3171.y, 0.0f);
    (&_S3089)->_S3086 = make_float2 (fx_13 * _S3172.x + cx_11, fy_13 * _S3172.y + cy_11);
    float2  _S3173 = float2 {_S3128.p_0[int(6)].x, _S3128.p_0[int(6)].y};
    float _S3174 = length_0(_S3173);
    float _S3175 = _S3128.p_0[int(6)].z;
    float _S3176 = s_primal_ctx_atan2_0(_S3174, _S3175);
    if(_S3176 < 0.00100000004749745f)
    {
        k_8 = (1.0f - _S3176 * _S3176 / 3.0f) / _S3175;
    }
    else
    {
        k_8 = _S3176 / _S3174;
    }
    float2  _S3177 = _S3173 * make_float2 (k_8);
    float u_97 = _S3177.x;
    float v_97 = _S3177.y;
    float r2_97 = u_97 * u_97 + v_97 * v_97;
    float2  _S3178 = _S3177 * make_float2 (1.0f + r2_97 * (dist_coeffs_13[int(0)] + r2_97 * (dist_coeffs_13[int(1)] + r2_97 * (dist_coeffs_13[int(2)] + r2_97 * dist_coeffs_13[int(3)])))) + make_float2 (_S3134 * u_97 * v_97 + dist_coeffs_13[int(5)] * (r2_97 + 2.0f * u_97 * u_97) + dist_coeffs_13[int(6)] * r2_97, _S3135 * u_97 * v_97 + dist_coeffs_13[int(4)] * (r2_97 + 2.0f * v_97 * v_97) + dist_coeffs_13[int(7)] * r2_97);
    float2  _S3179 = _S3178 + make_float2 (dist_coeffs_13[int(8)] * _S3178.x + dist_coeffs_13[int(9)] * _S3178.y, 0.0f);
    (&_S3089)->_S3087 = make_float2 (fx_13 * _S3179.x + cx_11, fy_13 * _S3179.y + cy_11);
    float3  mean_c_8 = s_primal_ctx_mul_0(R_10, mean_11) + t_10;
    float _S3180 = - in_opacity_10;
    float _S3181 = 1.0f + s_primal_ctx_exp_0(_S3180);
    float _S3182 = 1.0f / _S3181;
    float _S3183 = _S3181 * _S3181;
    float3  _S3184 = make_float3 (_S3100);
    float3  _S3185 = make_float3 (_S3103);
    float3  _S3186 = make_float3 (_S3106);
    float2  _S3187 = make_float2 (_S3098.w_mean_0[int(1)]) * _S3089._S3082 + make_float2 (_S3110.w_mean_0[int(2)]) * _S3089._S3083 + make_float2 (_S3112.w_mean_0[int(3)]) * _S3089._S3084 + make_float2 (_S3114.w_mean_0[int(4)]) * _S3089._S3085 + make_float2 (_S3116.w_mean_0[int(5)]) * _S3089._S3086 + make_float2 (_S3118.w_mean_0[int(6)]) * _S3089._S3087;
    float2  d_28 = _S3089._S3081 - _S3187;
    float _S3188 = d_28.x;
    float _S3189 = d_28.y;
    float _S3190 = _S3188 * _S3189;
    float2  d_29 = _S3089._S3082 - _S3187;
    float _S3191 = d_29.x;
    float _S3192 = d_29.y;
    float _S3193 = _S3191 * _S3192;
    float2  d_30 = _S3089._S3083 - _S3187;
    float _S3194 = d_30.x;
    float _S3195 = d_30.y;
    float _S3196 = _S3194 * _S3195;
    float2  d_31 = _S3089._S3084 - _S3187;
    float _S3197 = d_31.x;
    float _S3198 = d_31.y;
    float _S3199 = _S3197 * _S3198;
    float2  d_32 = _S3089._S3085 - _S3187;
    float _S3200 = d_32.x;
    float _S3201 = d_32.y;
    float _S3202 = _S3200 * _S3201;
    float2  d_33 = _S3089._S3086 - _S3187;
    float _S3203 = d_33.x;
    float _S3204 = d_33.y;
    float _S3205 = _S3203 * _S3204;
    float2  d_34 = _S3089._S3087 - _S3187;
    float _S3206 = d_34.x;
    float _S3207 = d_34.y;
    float _S3208 = _S3206 * _S3207;
    Matrix<float, 2, 2>  covar2d_7 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3188 * _S3188, _S3190, _S3190, _S3189 * _S3189) + makeMatrix<float, 2, 2> (_S3109.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3191 * _S3191, _S3193, _S3193, _S3192 * _S3192) + makeMatrix<float, 2, 2> (_S3111.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3194 * _S3194, _S3196, _S3196, _S3195 * _S3195) + makeMatrix<float, 2, 2> (_S3113.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3197 * _S3197, _S3199, _S3199, _S3198 * _S3198) + makeMatrix<float, 2, 2> (_S3115.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3200 * _S3200, _S3202, _S3202, _S3201 * _S3201) + makeMatrix<float, 2, 2> (_S3117.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3203 * _S3203, _S3205, _S3205, _S3204 * _S3204) + makeMatrix<float, 2, 2> (_S3119.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3206 * _S3206, _S3208, _S3208, _S3207 * _S3207);
    float eps2d_10;
    if(antialiased_10)
    {
        eps2d_10 = 0.10000000149011612f;
    }
    else
    {
        eps2d_10 = 0.30000001192092896f;
    }
    float _S3209 = covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x;
    float det_orig_10 = covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - _S3209;
    float _S3210 = covar2d_7.rows[int(0)].x + eps2d_10;
    Matrix<float, 2, 2>  _S3211 = covar2d_7;
    *&(((&_S3211)->rows + (int(0)))->x) = _S3210;
    float _S3212 = covar2d_7.rows[int(1)].y + eps2d_10;
    *&(((&_S3211)->rows + (int(1)))->y) = _S3212;
    Matrix<float, 2, 2>  _S3213 = _S3211;
    Matrix<float, 2, 2>  _S3214 = _S3211;
    float det_blur_10 = _S3210 * _S3212 - _S3209;
    float _S3215 = det_orig_10 / det_blur_10;
    float _S3216 = det_blur_10 * det_blur_10;
    float _S3217 = (F32_max((0.0f), (_S3215)));
    float _S3218 = s_primal_ctx_sqrt_0(_S3217);
    float invdet_12 = 1.0f / det_blur_10;
    float _S3219 = - covar2d_7.rows[int(0)].y;
    float _S3220 = - covar2d_7.rows[int(1)].x;
    if(antialiased_10)
    {
        k_8 = _S3182 * _S3218;
    }
    else
    {
        k_8 = _S3182;
    }
    float _S3221 = k_8 / 0.00392156885936856f;
    float _S3222 = 2.0f * s_primal_ctx_log_0(_S3221);
    float _S3223 = s_primal_ctx_sqrt_0(_S3222);
    float _S3224 = _S3213.rows[int(0)].x;
    float _S3225 = _S3214.rows[int(1)].y;
    float3  campos_5 = - s_primal_ctx_mul_0(transpose_3(R_10), t_10);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3226;
    (&_S3226)->primal_0 = mean_11;
    (&_S3226)->differential_0 = _S3090;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3227;
    (&_S3227)->primal_0 = scale_10;
    (&_S3227)->differential_0 = _S3090;
    DiffPair_float_0 _S3228;
    (&_S3228)->primal_0 = in_opacity_10;
    (&_S3228)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3229;
    (&_S3229)->primal_0 = campos_5;
    (&_S3229)->differential_0 = _S3090;
    s_bwd_prop_view_radius_3dgs_0(&_S3226, &_S3227, &_S3228, &_S3229, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3230 = _S3226;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3231 = _S3227;
    DiffPair_float_0 _S3232 = _S3228;
    float2  _S3233 = _S3088;
    *&((&_S3233)->y) = v_conic_4.z;
    float2  _S3234 = _S3088;
    *&((&_S3234)->y) = v_conic_4.y;
    *&((&_S3234)->x) = v_conic_4.x;
    DiffPair_float_0 _S3235;
    (&_S3235)->primal_0 = _S3225;
    (&_S3235)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3235, 0.0f);
    DiffPair_float_0 _S3236;
    (&_S3236)->primal_0 = _S3224;
    (&_S3236)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3236, 0.0f);
    DiffPair_float_0 _S3237;
    (&_S3237)->primal_0 = 3.32999992370605469f;
    (&_S3237)->differential_0 = 0.0f;
    DiffPair_float_0 _S3238;
    (&_S3238)->primal_0 = _S3223;
    (&_S3238)->differential_0 = 0.0f;
    _d_min_0(&_S3237, &_S3238, 0.0f);
    DiffPair_float_0 _S3239;
    (&_S3239)->primal_0 = _S3222;
    (&_S3239)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3239, _S3238.differential_0);
    float _S3240 = 2.0f * _S3239.differential_0;
    DiffPair_float_0 _S3241;
    (&_S3241)->primal_0 = _S3221;
    (&_S3241)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3241, _S3240);
    float _S3242 = v_opacity_4 + 254.9999847412109375f * _S3241.differential_0;
    Matrix<float, 2, 2>  _S3243 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3244 = _S3243;
    _S3244[int(1)] = _S3233;
    _S3244[int(0)] = _S3234;
    Matrix<float, 2, 2>  _S3245 = _S3244;
    float2  _S3246 = make_float2 (0.0f, _S3235.differential_0);
    float2  _S3247 = make_float2 (_S3236.differential_0, 0.0f);
    if(antialiased_10)
    {
        float _S3248 = _S3218 * _S3242;
        k_8 = _S3182 * _S3242;
        eps2d_10 = _S3248;
    }
    else
    {
        k_8 = 0.0f;
        eps2d_10 = _S3242;
    }
    float _S3249 = invdet_12 * _S3245.rows[int(1)].y;
    float _S3250 = - (invdet_12 * _S3245.rows[int(1)].x);
    float _S3251 = - (invdet_12 * _S3245.rows[int(0)].y);
    float _S3252 = invdet_12 * _S3245.rows[int(0)].x;
    float _S3253 = - ((_S3210 * _S3245.rows[int(1)].y + _S3220 * _S3245.rows[int(1)].x + _S3219 * _S3245.rows[int(0)].y + _S3212 * _S3245.rows[int(0)].x) / _S3216);
    DiffPair_float_0 _S3254;
    (&_S3254)->primal_0 = _S3217;
    (&_S3254)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3254, k_8);
    DiffPair_float_0 _S3255;
    (&_S3255)->primal_0 = 0.0f;
    (&_S3255)->differential_0 = 0.0f;
    DiffPair_float_0 _S3256;
    (&_S3256)->primal_0 = _S3215;
    (&_S3256)->differential_0 = 0.0f;
    _d_max_0(&_S3255, &_S3256, _S3254.differential_0);
    float _S3257 = _S3256.differential_0 / _S3216;
    float s_diff_det_orig_T_4 = det_blur_10 * _S3257;
    float _S3258 = det_orig_10 * - _S3257 + _S3253;
    float _S3259 = - _S3258;
    float _S3260 = _S3210 * _S3258;
    float _S3261 = _S3212 * _S3258;
    Matrix<float, 2, 2>  _S3262 = _S3243;
    _S3262[int(1)] = _S3246;
    _S3262[int(0)] = _S3247;
    float _S3263 = _S3261 + _S3262.rows[int(0)].x + _S3249;
    float _S3264 = _S3259 + - s_diff_det_orig_T_4;
    float _S3265 = covar2d_7.rows[int(0)].y * _S3264 + _S3250;
    float _S3266 = covar2d_7.rows[int(1)].x * _S3264 + _S3251;
    float _S3267 = covar2d_7.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S3268 = _S3260 + _S3262.rows[int(1)].y + _S3252 + covar2d_7.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S3269 = _S3088;
    *&((&_S3269)->x) = _S3265;
    *&((&_S3269)->y) = _S3268;
    float _S3270 = _S3263 + _S3267;
    float2  _S3271 = _S3088;
    *&((&_S3271)->y) = _S3266;
    *&((&_S3271)->x) = _S3270;
    Matrix<float, 2, 2>  _S3272 = _S3243;
    _S3272[int(1)] = _S3269;
    _S3272[int(0)] = _S3271;
    Matrix<float, 3, 3>  _S3273 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3274;
    (&_S3274)->primal_0 = R_10;
    (&_S3274)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3275;
    (&_S3275)->primal_0 = _S3108;
    (&_S3275)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3274, &_S3275, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3276;
    (&_S3276)->primal_0 = R_10;
    (&_S3276)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3277;
    (&_S3277)->primal_0 = _S3105;
    (&_S3277)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3276, &_S3277, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3278;
    (&_S3278)->primal_0 = R_10;
    (&_S3278)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3279;
    (&_S3279)->primal_0 = _S3102;
    (&_S3279)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3278, &_S3279, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3280;
    (&_S3280)->primal_0 = R_10;
    (&_S3280)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3281;
    (&_S3281)->primal_0 = _S3107;
    (&_S3281)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3280, &_S3281, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3282;
    (&_S3282)->primal_0 = R_10;
    (&_S3282)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3283;
    (&_S3283)->primal_0 = _S3104;
    (&_S3283)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3282, &_S3283, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3284;
    (&_S3284)->primal_0 = R_10;
    (&_S3284)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3285;
    (&_S3285)->primal_0 = _S3101;
    (&_S3285)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3284, &_S3285, _S3090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3286;
    (&_S3286)->primal_0 = R_10;
    (&_S3286)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3287;
    (&_S3287)->primal_0 = _S3120.p_0[0U];
    (&_S3287)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3286, &_S3287, _S3090);
    float3  _S3288 = - _S3275.differential_0 + _S3281.differential_0;
    float3  _S3289 = _S3186 * _S3288;
    float3  _S3290 = _S3094.rows[2U] * _S3288;
    float _S3291 = _S3099 * (_S3290.x + _S3290.y + _S3290.z);
    float3  _S3292 = - _S3277.differential_0 + _S3283.differential_0;
    float3  _S3293 = _S3185 * _S3292;
    float3  _S3294 = _S3094.rows[1U] * _S3292;
    float _S3295 = _S3099 * (_S3294.x + _S3294.y + _S3294.z);
    float3  _S3296 = - _S3279.differential_0 + _S3285.differential_0;
    float3  _S3297 = _S3184 * _S3296;
    float3  _S3298 = _S3094.rows[0U] * _S3296;
    float _S3299 = _S3099 * (_S3298.x + _S3298.y + _S3298.z);
    Matrix<float, 3, 3>  _S3300 = _S3273;
    _S3300[2U] = _S3289;
    _S3300[1U] = _S3293;
    _S3300[0U] = _S3297;
    Matrix<float, 3, 3>  _S3301 = transpose_3(transpose_3(_S3300));
    float _S3302 = 2.0f * - _S3301.rows[int(2)].z;
    float _S3303 = 2.0f * _S3301.rows[int(2)].y;
    float _S3304 = 2.0f * _S3301.rows[int(2)].x;
    float _S3305 = 2.0f * _S3301.rows[int(1)].z;
    float _S3306 = 2.0f * - _S3301.rows[int(1)].y;
    float _S3307 = 2.0f * _S3301.rows[int(1)].x;
    float _S3308 = 2.0f * _S3301.rows[int(0)].z;
    float _S3309 = 2.0f * _S3301.rows[int(0)].y;
    float _S3310 = 2.0f * - _S3301.rows[int(0)].x;
    float _S3311 = - _S3307 + _S3309;
    float _S3312 = _S3304 + - _S3308;
    float _S3313 = - _S3303 + _S3305;
    float _S3314 = _S3303 + _S3305;
    float _S3315 = _S3304 + _S3308;
    float _S3316 = _S3307 + _S3309;
    float _S3317 = _S3092.w * (_S3306 + _S3310);
    float _S3318 = _S3092.z * (_S3302 + _S3310);
    float _S3319 = _S3092.y * (_S3302 + _S3306);
    float _S3320 = _S3092.x * _S3311 + _S3092.z * _S3314 + _S3092.y * _S3315 + _S3317 + _S3317;
    float _S3321 = _S3092.x * _S3312 + _S3092.w * _S3314 + _S3092.y * _S3316 + _S3318 + _S3318;
    float _S3322 = _S3092.x * _S3313 + _S3092.w * _S3315 + _S3092.z * _S3316 + _S3319 + _S3319;
    float _S3323 = _S3092.w * _S3311 + _S3092.z * _S3312 + _S3092.y * _S3313;
    float4  _S3324 = make_float4 (0.0f);
    float4  _S3325 = _S3324;
    *&((&_S3325)->w) = _S3320;
    *&((&_S3325)->z) = _S3321;
    *&((&_S3325)->y) = _S3322;
    *&((&_S3325)->x) = _S3323;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3326;
    (&_S3326)->primal_0 = quat_10;
    (&_S3326)->differential_0 = _S3324;
    s_bwd_normalize_impl_0(&_S3326, _S3325);
    float3  _S3327 = _S3090;
    *&((&_S3327)->z) = _S3291;
    *&((&_S3327)->y) = _S3295;
    *&((&_S3327)->x) = _S3299;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3328;
    (&_S3328)->primal_0 = scale_10;
    (&_S3328)->differential_0 = _S3090;
    s_bwd_prop_exp_1(&_S3328, _S3327);
    float _S3329 = - (eps2d_10 / _S3183);
    DiffPair_float_0 _S3330;
    (&_S3330)->primal_0 = _S3180;
    (&_S3330)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3330, _S3329);
    float _S3331 = - _S3330.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3332;
    (&_S3332)->primal_0 = mean_c_8;
    (&_S3332)->differential_0 = _S3090;
    s_bwd_length_impl_0(&_S3332, 0.0f);
    float3  _S3333 = _S3332.differential_0 + make_float3 (0.0f, 0.0f, v_depth_4);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3334;
    (&_S3334)->primal_0 = R_10;
    (&_S3334)->differential_0 = _S3273;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3335;
    (&_S3335)->primal_0 = mean_11;
    (&_S3335)->differential_0 = _S3090;
    s_bwd_prop_mul_3(&_S3334, &_S3335, _S3333);
    Matrix<float, 3, 3>  _S3336 = _S3274.differential_0 + _S3276.differential_0 + _S3278.differential_0 + _S3280.differential_0 + _S3282.differential_0 + _S3284.differential_0 + _S3286.differential_0 + _S3334.differential_0;
    float _S3337 = _S3331 + _S3232.differential_0;
    float3  _S3338 = _S3328.differential_0 + _S3231.differential_0;
    *v_mean_4 = *v_mean_4 + (_S3275.differential_0 + _S3281.differential_0 + _S3277.differential_0 + _S3283.differential_0 + _S3279.differential_0 + _S3285.differential_0 + _S3335.differential_0 + _S3230.differential_0);
    *v_quat_4 = *v_quat_4 + _S3326.differential_0;
    *v_scale_4 = *v_scale_4 + _S3338;
    *v_in_opacity_4 = *v_in_opacity_4 + _S3337;
    *v_R_4 = *v_R_4 + _S3336;
    *v_t_4 = *v_t_4 + _S3333;
    return;
}

struct s_bwd_prop_DiffProjection3DGS_3dgut_equisolid_projection_Intermediates_0
{
    float2  _S3339;
    float2  _S3340;
    float2  _S3341;
    float2  _S3342;
    float2  _S3343;
    float2  _S3344;
    float2  _S3345;
};

inline __device__ void projection_3dgut_equisolid_vjp(bool antialiased_11, float3  mean_12, float4  quat_11, float3  scale_11, float in_opacity_11, Matrix<float, 3, 3>  R_11, float3  t_11, float fx_14, float fy_14, float cx_12, float cy_12, FixedArray<float, 10>  dist_coeffs_14, uint image_width_11, uint image_height_11, float2  v_mean2d_5, float v_depth_5, float3  v_conic_5, float v_opacity_5, float3  * v_mean_5, float4  * v_quat_5, float3  * v_scale_5, float * v_in_opacity_5, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    float2  _S3346 = make_float2 (0.0f);
    s_bwd_prop_DiffProjection3DGS_3dgut_equisolid_projection_Intermediates_0 _S3347;
    (&_S3347)->_S3339 = _S3346;
    (&_S3347)->_S3340 = _S3346;
    (&_S3347)->_S3341 = _S3346;
    (&_S3347)->_S3342 = _S3346;
    (&_S3347)->_S3343 = _S3346;
    (&_S3347)->_S3344 = _S3346;
    (&_S3347)->_S3345 = _S3346;
    (&_S3347)->_S3339 = _S3346;
    (&_S3347)->_S3340 = _S3346;
    (&_S3347)->_S3341 = _S3346;
    (&_S3347)->_S3342 = _S3346;
    (&_S3347)->_S3343 = _S3346;
    (&_S3347)->_S3344 = _S3346;
    (&_S3347)->_S3345 = _S3346;
    float3  _S3348 = make_float3 (0.0f);
    float3  _S3349 = s_primal_ctx_exp_1(scale_11);
    float4  _S3350 = normalize_0(quat_11);
    float _S3351 = _S3350.y;
    float x2_11 = _S3351 * _S3351;
    float y2_11 = _S3350.z * _S3350.z;
    float z2_11 = _S3350.w * _S3350.w;
    float xy_11 = _S3350.y * _S3350.z;
    float xz_11 = _S3350.y * _S3350.w;
    float yz_11 = _S3350.z * _S3350.w;
    float wx_11 = _S3350.x * _S3350.y;
    float wy_11 = _S3350.x * _S3350.z;
    float wz_11 = _S3350.x * _S3350.w;
    Matrix<float, 3, 3>  _S3352 = transpose_3(transpose_3(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_11), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_11), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))));
    FixedArray<float3 , 7>  _S3353 = {
        _S3348, _S3348, _S3348, _S3348, _S3348, _S3348, _S3348
    };
    FixedArray<float, 7>  _S3354 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S3355;
    (&_S3355)->p_0 = _S3353;
    (&_S3355)->w_mean_0 = _S3354;
    (&_S3355)->w_cov_0 = _S3354;
    (&_S3355)->p_0[int(0)] = mean_12;
    SigmaPoints_0 _S3356 = _S3355;
    (&_S3356)->w_mean_0[int(0)] = 0.0f;
    (&_S3356)->w_cov_0[int(0)] = 2.0f;
    float _S3357 = s_primal_ctx_sqrt_0(3.0f);
    float _S3358 = _S3357 * _S3349.x;
    float3  delta_15 = make_float3 (_S3358) * _S3352.rows[0U];
    float3  _S3359 = mean_12 + delta_15;
    (&_S3356)->p_0[1U] = _S3359;
    float3  _S3360 = mean_12 - delta_15;
    (&_S3356)->p_0[4U] = _S3360;
    float _S3361 = _S3357 * _S3349.y;
    float3  delta_16 = make_float3 (_S3361) * _S3352.rows[1U];
    float3  _S3362 = mean_12 + delta_16;
    (&_S3356)->p_0[2U] = _S3362;
    float3  _S3363 = mean_12 - delta_16;
    (&_S3356)->p_0[5U] = _S3363;
    float _S3364 = _S3357 * _S3349.z;
    float3  delta_17 = make_float3 (_S3364) * _S3352.rows[2U];
    float3  _S3365 = mean_12 + delta_17;
    (&_S3356)->p_0[3U] = _S3365;
    float3  _S3366 = mean_12 - delta_17;
    (&_S3356)->p_0[6U] = _S3366;
    (&_S3356)->w_mean_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3367 = _S3356;
    (&_S3367)->w_cov_0[1U] = 0.1666666716337204f;
    SigmaPoints_0 _S3368 = _S3367;
    (&_S3368)->w_mean_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3369 = _S3368;
    (&_S3369)->w_cov_0[2U] = 0.1666666716337204f;
    SigmaPoints_0 _S3370 = _S3369;
    (&_S3370)->w_mean_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3371 = _S3370;
    (&_S3371)->w_cov_0[3U] = 0.1666666716337204f;
    SigmaPoints_0 _S3372 = _S3371;
    (&_S3372)->w_mean_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3373 = _S3372;
    (&_S3373)->w_cov_0[4U] = 0.1666666716337204f;
    SigmaPoints_0 _S3374 = _S3373;
    (&_S3374)->w_mean_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3375 = _S3374;
    (&_S3375)->w_cov_0[5U] = 0.1666666716337204f;
    SigmaPoints_0 _S3376 = _S3375;
    (&_S3376)->w_mean_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3377 = _S3376;
    (&_S3377)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S3378 = _S3355;
    float3  _S3379 = s_primal_ctx_mul_0(R_11, _S3355.p_0[0U]) + t_11;
    _S3355 = _S3377;
    (&_S3355)->p_0[0U] = _S3379;
    SigmaPoints_0 _S3380 = _S3355;
    (&_S3355)->p_0[1U] = s_primal_ctx_mul_0(R_11, _S3359) + t_11;
    SigmaPoints_0 _S3381 = _S3355;
    (&_S3355)->p_0[2U] = s_primal_ctx_mul_0(R_11, _S3362) + t_11;
    SigmaPoints_0 _S3382 = _S3355;
    (&_S3355)->p_0[3U] = s_primal_ctx_mul_0(R_11, _S3365) + t_11;
    SigmaPoints_0 _S3383 = _S3355;
    (&_S3355)->p_0[4U] = s_primal_ctx_mul_0(R_11, _S3360) + t_11;
    SigmaPoints_0 _S3384 = _S3355;
    (&_S3355)->p_0[5U] = s_primal_ctx_mul_0(R_11, _S3363) + t_11;
    SigmaPoints_0 _S3385 = _S3355;
    (&_S3355)->p_0[6U] = s_primal_ctx_mul_0(R_11, _S3366) + t_11;
    SigmaPoints_0 _S3386 = _S3355;
    float2  _S3387 = float2 {_S3380.p_0[int(0)].x, _S3380.p_0[int(0)].y};
    float _S3388 = length_0(_S3387);
    float _S3389 = _S3380.p_0[int(0)].z;
    float _S3390 = s_primal_ctx_atan2_0(_S3388, _S3389);
    float k_9;
    if(_S3388 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3390 * _S3390 / 24.0f) / _S3389;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3390) / _S3388;
    }
    float2  _S3391 = _S3387 * make_float2 (k_9);
    float u_98 = _S3391.x;
    float v_98 = _S3391.y;
    float r2_98 = u_98 * u_98 + v_98 * v_98;
    float _S3392 = 2.0f * dist_coeffs_14[int(4)];
    float _S3393 = 2.0f * dist_coeffs_14[int(5)];
    float2  _S3394 = _S3391 * make_float2 (1.0f + r2_98 * (dist_coeffs_14[int(0)] + r2_98 * (dist_coeffs_14[int(1)] + r2_98 * (dist_coeffs_14[int(2)] + r2_98 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_98 * v_98 + dist_coeffs_14[int(5)] * (r2_98 + 2.0f * u_98 * u_98) + dist_coeffs_14[int(6)] * r2_98, _S3393 * u_98 * v_98 + dist_coeffs_14[int(4)] * (r2_98 + 2.0f * v_98 * v_98) + dist_coeffs_14[int(7)] * r2_98);
    float2  _S3395 = _S3394 + make_float2 (dist_coeffs_14[int(8)] * _S3394.x + dist_coeffs_14[int(9)] * _S3394.y, 0.0f);
    (&_S3347)->_S3339 = make_float2 (fx_14 * _S3395.x + cx_12, fy_14 * _S3395.y + cy_12);
    float2  _S3396 = float2 {_S3381.p_0[int(1)].x, _S3381.p_0[int(1)].y};
    float _S3397 = length_0(_S3396);
    float _S3398 = _S3381.p_0[int(1)].z;
    float _S3399 = s_primal_ctx_atan2_0(_S3397, _S3398);
    if(_S3397 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3399 * _S3399 / 24.0f) / _S3398;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3399) / _S3397;
    }
    float2  _S3400 = _S3396 * make_float2 (k_9);
    float u_99 = _S3400.x;
    float v_99 = _S3400.y;
    float r2_99 = u_99 * u_99 + v_99 * v_99;
    float2  _S3401 = _S3400 * make_float2 (1.0f + r2_99 * (dist_coeffs_14[int(0)] + r2_99 * (dist_coeffs_14[int(1)] + r2_99 * (dist_coeffs_14[int(2)] + r2_99 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_99 * v_99 + dist_coeffs_14[int(5)] * (r2_99 + 2.0f * u_99 * u_99) + dist_coeffs_14[int(6)] * r2_99, _S3393 * u_99 * v_99 + dist_coeffs_14[int(4)] * (r2_99 + 2.0f * v_99 * v_99) + dist_coeffs_14[int(7)] * r2_99);
    float2  _S3402 = _S3401 + make_float2 (dist_coeffs_14[int(8)] * _S3401.x + dist_coeffs_14[int(9)] * _S3401.y, 0.0f);
    (&_S3347)->_S3340 = make_float2 (fx_14 * _S3402.x + cx_12, fy_14 * _S3402.y + cy_12);
    float2  _S3403 = float2 {_S3382.p_0[int(2)].x, _S3382.p_0[int(2)].y};
    float _S3404 = length_0(_S3403);
    float _S3405 = _S3382.p_0[int(2)].z;
    float _S3406 = s_primal_ctx_atan2_0(_S3404, _S3405);
    if(_S3404 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3406 * _S3406 / 24.0f) / _S3405;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3406) / _S3404;
    }
    float2  _S3407 = _S3403 * make_float2 (k_9);
    float u_100 = _S3407.x;
    float v_100 = _S3407.y;
    float r2_100 = u_100 * u_100 + v_100 * v_100;
    float2  _S3408 = _S3407 * make_float2 (1.0f + r2_100 * (dist_coeffs_14[int(0)] + r2_100 * (dist_coeffs_14[int(1)] + r2_100 * (dist_coeffs_14[int(2)] + r2_100 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_100 * v_100 + dist_coeffs_14[int(5)] * (r2_100 + 2.0f * u_100 * u_100) + dist_coeffs_14[int(6)] * r2_100, _S3393 * u_100 * v_100 + dist_coeffs_14[int(4)] * (r2_100 + 2.0f * v_100 * v_100) + dist_coeffs_14[int(7)] * r2_100);
    float2  _S3409 = _S3408 + make_float2 (dist_coeffs_14[int(8)] * _S3408.x + dist_coeffs_14[int(9)] * _S3408.y, 0.0f);
    (&_S3347)->_S3341 = make_float2 (fx_14 * _S3409.x + cx_12, fy_14 * _S3409.y + cy_12);
    float2  _S3410 = float2 {_S3383.p_0[int(3)].x, _S3383.p_0[int(3)].y};
    float _S3411 = length_0(_S3410);
    float _S3412 = _S3383.p_0[int(3)].z;
    float _S3413 = s_primal_ctx_atan2_0(_S3411, _S3412);
    if(_S3411 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3413 * _S3413 / 24.0f) / _S3412;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3413) / _S3411;
    }
    float2  _S3414 = _S3410 * make_float2 (k_9);
    float u_101 = _S3414.x;
    float v_101 = _S3414.y;
    float r2_101 = u_101 * u_101 + v_101 * v_101;
    float2  _S3415 = _S3414 * make_float2 (1.0f + r2_101 * (dist_coeffs_14[int(0)] + r2_101 * (dist_coeffs_14[int(1)] + r2_101 * (dist_coeffs_14[int(2)] + r2_101 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_101 * v_101 + dist_coeffs_14[int(5)] * (r2_101 + 2.0f * u_101 * u_101) + dist_coeffs_14[int(6)] * r2_101, _S3393 * u_101 * v_101 + dist_coeffs_14[int(4)] * (r2_101 + 2.0f * v_101 * v_101) + dist_coeffs_14[int(7)] * r2_101);
    float2  _S3416 = _S3415 + make_float2 (dist_coeffs_14[int(8)] * _S3415.x + dist_coeffs_14[int(9)] * _S3415.y, 0.0f);
    (&_S3347)->_S3342 = make_float2 (fx_14 * _S3416.x + cx_12, fy_14 * _S3416.y + cy_12);
    float2  _S3417 = float2 {_S3384.p_0[int(4)].x, _S3384.p_0[int(4)].y};
    float _S3418 = length_0(_S3417);
    float _S3419 = _S3384.p_0[int(4)].z;
    float _S3420 = s_primal_ctx_atan2_0(_S3418, _S3419);
    if(_S3418 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3420 * _S3420 / 24.0f) / _S3419;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3420) / _S3418;
    }
    float2  _S3421 = _S3417 * make_float2 (k_9);
    float u_102 = _S3421.x;
    float v_102 = _S3421.y;
    float r2_102 = u_102 * u_102 + v_102 * v_102;
    float2  _S3422 = _S3421 * make_float2 (1.0f + r2_102 * (dist_coeffs_14[int(0)] + r2_102 * (dist_coeffs_14[int(1)] + r2_102 * (dist_coeffs_14[int(2)] + r2_102 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_102 * v_102 + dist_coeffs_14[int(5)] * (r2_102 + 2.0f * u_102 * u_102) + dist_coeffs_14[int(6)] * r2_102, _S3393 * u_102 * v_102 + dist_coeffs_14[int(4)] * (r2_102 + 2.0f * v_102 * v_102) + dist_coeffs_14[int(7)] * r2_102);
    float2  _S3423 = _S3422 + make_float2 (dist_coeffs_14[int(8)] * _S3422.x + dist_coeffs_14[int(9)] * _S3422.y, 0.0f);
    (&_S3347)->_S3343 = make_float2 (fx_14 * _S3423.x + cx_12, fy_14 * _S3423.y + cy_12);
    float2  _S3424 = float2 {_S3385.p_0[int(5)].x, _S3385.p_0[int(5)].y};
    float _S3425 = length_0(_S3424);
    float _S3426 = _S3385.p_0[int(5)].z;
    float _S3427 = s_primal_ctx_atan2_0(_S3425, _S3426);
    if(_S3425 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3427 * _S3427 / 24.0f) / _S3426;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3427) / _S3425;
    }
    float2  _S3428 = _S3424 * make_float2 (k_9);
    float u_103 = _S3428.x;
    float v_103 = _S3428.y;
    float r2_103 = u_103 * u_103 + v_103 * v_103;
    float2  _S3429 = _S3428 * make_float2 (1.0f + r2_103 * (dist_coeffs_14[int(0)] + r2_103 * (dist_coeffs_14[int(1)] + r2_103 * (dist_coeffs_14[int(2)] + r2_103 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_103 * v_103 + dist_coeffs_14[int(5)] * (r2_103 + 2.0f * u_103 * u_103) + dist_coeffs_14[int(6)] * r2_103, _S3393 * u_103 * v_103 + dist_coeffs_14[int(4)] * (r2_103 + 2.0f * v_103 * v_103) + dist_coeffs_14[int(7)] * r2_103);
    float2  _S3430 = _S3429 + make_float2 (dist_coeffs_14[int(8)] * _S3429.x + dist_coeffs_14[int(9)] * _S3429.y, 0.0f);
    (&_S3347)->_S3344 = make_float2 (fx_14 * _S3430.x + cx_12, fy_14 * _S3430.y + cy_12);
    float2  _S3431 = float2 {_S3386.p_0[int(6)].x, _S3386.p_0[int(6)].y};
    float _S3432 = length_0(_S3431);
    float _S3433 = _S3386.p_0[int(6)].z;
    float _S3434 = s_primal_ctx_atan2_0(_S3432, _S3433);
    if(_S3432 < 9.99999997475242708e-07f)
    {
        k_9 = (1.0f - _S3434 * _S3434 / 24.0f) / _S3433;
    }
    else
    {
        k_9 = 2.0f * s_primal_ctx_sin_0(0.5f * _S3434) / _S3432;
    }
    float2  _S3435 = _S3431 * make_float2 (k_9);
    float u_104 = _S3435.x;
    float v_104 = _S3435.y;
    float r2_104 = u_104 * u_104 + v_104 * v_104;
    float2  _S3436 = _S3435 * make_float2 (1.0f + r2_104 * (dist_coeffs_14[int(0)] + r2_104 * (dist_coeffs_14[int(1)] + r2_104 * (dist_coeffs_14[int(2)] + r2_104 * dist_coeffs_14[int(3)])))) + make_float2 (_S3392 * u_104 * v_104 + dist_coeffs_14[int(5)] * (r2_104 + 2.0f * u_104 * u_104) + dist_coeffs_14[int(6)] * r2_104, _S3393 * u_104 * v_104 + dist_coeffs_14[int(4)] * (r2_104 + 2.0f * v_104 * v_104) + dist_coeffs_14[int(7)] * r2_104);
    float2  _S3437 = _S3436 + make_float2 (dist_coeffs_14[int(8)] * _S3436.x + dist_coeffs_14[int(9)] * _S3436.y, 0.0f);
    (&_S3347)->_S3345 = make_float2 (fx_14 * _S3437.x + cx_12, fy_14 * _S3437.y + cy_12);
    float3  mean_c_9 = s_primal_ctx_mul_0(R_11, mean_12) + t_11;
    float _S3438 = - in_opacity_11;
    float _S3439 = 1.0f + s_primal_ctx_exp_0(_S3438);
    float _S3440 = 1.0f / _S3439;
    float _S3441 = _S3439 * _S3439;
    float3  _S3442 = make_float3 (_S3358);
    float3  _S3443 = make_float3 (_S3361);
    float3  _S3444 = make_float3 (_S3364);
    float2  _S3445 = make_float2 (_S3356.w_mean_0[int(1)]) * _S3347._S3340 + make_float2 (_S3368.w_mean_0[int(2)]) * _S3347._S3341 + make_float2 (_S3370.w_mean_0[int(3)]) * _S3347._S3342 + make_float2 (_S3372.w_mean_0[int(4)]) * _S3347._S3343 + make_float2 (_S3374.w_mean_0[int(5)]) * _S3347._S3344 + make_float2 (_S3376.w_mean_0[int(6)]) * _S3347._S3345;
    float2  d_35 = _S3347._S3339 - _S3445;
    float _S3446 = d_35.x;
    float _S3447 = d_35.y;
    float _S3448 = _S3446 * _S3447;
    float2  d_36 = _S3347._S3340 - _S3445;
    float _S3449 = d_36.x;
    float _S3450 = d_36.y;
    float _S3451 = _S3449 * _S3450;
    float2  d_37 = _S3347._S3341 - _S3445;
    float _S3452 = d_37.x;
    float _S3453 = d_37.y;
    float _S3454 = _S3452 * _S3453;
    float2  d_38 = _S3347._S3342 - _S3445;
    float _S3455 = d_38.x;
    float _S3456 = d_38.y;
    float _S3457 = _S3455 * _S3456;
    float2  d_39 = _S3347._S3343 - _S3445;
    float _S3458 = d_39.x;
    float _S3459 = d_39.y;
    float _S3460 = _S3458 * _S3459;
    float2  d_40 = _S3347._S3344 - _S3445;
    float _S3461 = d_40.x;
    float _S3462 = d_40.y;
    float _S3463 = _S3461 * _S3462;
    float2  d_41 = _S3347._S3345 - _S3445;
    float _S3464 = d_41.x;
    float _S3465 = d_41.y;
    float _S3466 = _S3464 * _S3465;
    Matrix<float, 2, 2>  covar2d_8 = makeMatrix<float, 2, 2> (2.0f) * makeMatrix<float, 2, 2> (_S3446 * _S3446, _S3448, _S3448, _S3447 * _S3447) + makeMatrix<float, 2, 2> (_S3367.w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S3449 * _S3449, _S3451, _S3451, _S3450 * _S3450) + makeMatrix<float, 2, 2> (_S3369.w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S3452 * _S3452, _S3454, _S3454, _S3453 * _S3453) + makeMatrix<float, 2, 2> (_S3371.w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S3455 * _S3455, _S3457, _S3457, _S3456 * _S3456) + makeMatrix<float, 2, 2> (_S3373.w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S3458 * _S3458, _S3460, _S3460, _S3459 * _S3459) + makeMatrix<float, 2, 2> (_S3375.w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S3461 * _S3461, _S3463, _S3463, _S3462 * _S3462) + makeMatrix<float, 2, 2> (_S3377.w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S3464 * _S3464, _S3466, _S3466, _S3465 * _S3465);
    float eps2d_11;
    if(antialiased_11)
    {
        eps2d_11 = 0.10000000149011612f;
    }
    else
    {
        eps2d_11 = 0.30000001192092896f;
    }
    float _S3467 = covar2d_8.rows[int(0)].y * covar2d_8.rows[int(1)].x;
    float det_orig_11 = covar2d_8.rows[int(0)].x * covar2d_8.rows[int(1)].y - _S3467;
    float _S3468 = covar2d_8.rows[int(0)].x + eps2d_11;
    Matrix<float, 2, 2>  _S3469 = covar2d_8;
    *&(((&_S3469)->rows + (int(0)))->x) = _S3468;
    float _S3470 = covar2d_8.rows[int(1)].y + eps2d_11;
    *&(((&_S3469)->rows + (int(1)))->y) = _S3470;
    Matrix<float, 2, 2>  _S3471 = _S3469;
    Matrix<float, 2, 2>  _S3472 = _S3469;
    float det_blur_11 = _S3468 * _S3470 - _S3467;
    float _S3473 = det_orig_11 / det_blur_11;
    float _S3474 = det_blur_11 * det_blur_11;
    float _S3475 = (F32_max((0.0f), (_S3473)));
    float _S3476 = s_primal_ctx_sqrt_0(_S3475);
    float invdet_13 = 1.0f / det_blur_11;
    float _S3477 = - covar2d_8.rows[int(0)].y;
    float _S3478 = - covar2d_8.rows[int(1)].x;
    if(antialiased_11)
    {
        k_9 = _S3440 * _S3476;
    }
    else
    {
        k_9 = _S3440;
    }
    float _S3479 = k_9 / 0.00392156885936856f;
    float _S3480 = 2.0f * s_primal_ctx_log_0(_S3479);
    float _S3481 = s_primal_ctx_sqrt_0(_S3480);
    float _S3482 = _S3471.rows[int(0)].x;
    float _S3483 = _S3472.rows[int(1)].y;
    float3  campos_6 = - s_primal_ctx_mul_0(transpose_3(R_11), t_11);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3484;
    (&_S3484)->primal_0 = mean_12;
    (&_S3484)->differential_0 = _S3348;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3485;
    (&_S3485)->primal_0 = scale_11;
    (&_S3485)->differential_0 = _S3348;
    DiffPair_float_0 _S3486;
    (&_S3486)->primal_0 = in_opacity_11;
    (&_S3486)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3487;
    (&_S3487)->primal_0 = campos_6;
    (&_S3487)->differential_0 = _S3348;
    s_bwd_prop_view_radius_3dgs_0(&_S3484, &_S3485, &_S3486, &_S3487, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3488 = _S3484;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3489 = _S3485;
    DiffPair_float_0 _S3490 = _S3486;
    float2  _S3491 = _S3346;
    *&((&_S3491)->y) = v_conic_5.z;
    float2  _S3492 = _S3346;
    *&((&_S3492)->y) = v_conic_5.y;
    *&((&_S3492)->x) = v_conic_5.x;
    DiffPair_float_0 _S3493;
    (&_S3493)->primal_0 = _S3483;
    (&_S3493)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3493, 0.0f);
    DiffPair_float_0 _S3494;
    (&_S3494)->primal_0 = _S3482;
    (&_S3494)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3494, 0.0f);
    DiffPair_float_0 _S3495;
    (&_S3495)->primal_0 = 3.32999992370605469f;
    (&_S3495)->differential_0 = 0.0f;
    DiffPair_float_0 _S3496;
    (&_S3496)->primal_0 = _S3481;
    (&_S3496)->differential_0 = 0.0f;
    _d_min_0(&_S3495, &_S3496, 0.0f);
    DiffPair_float_0 _S3497;
    (&_S3497)->primal_0 = _S3480;
    (&_S3497)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3497, _S3496.differential_0);
    float _S3498 = 2.0f * _S3497.differential_0;
    DiffPair_float_0 _S3499;
    (&_S3499)->primal_0 = _S3479;
    (&_S3499)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3499, _S3498);
    float _S3500 = v_opacity_5 + 254.9999847412109375f * _S3499.differential_0;
    Matrix<float, 2, 2>  _S3501 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S3502 = _S3501;
    _S3502[int(1)] = _S3491;
    _S3502[int(0)] = _S3492;
    Matrix<float, 2, 2>  _S3503 = _S3502;
    float2  _S3504 = make_float2 (0.0f, _S3493.differential_0);
    float2  _S3505 = make_float2 (_S3494.differential_0, 0.0f);
    if(antialiased_11)
    {
        float _S3506 = _S3476 * _S3500;
        k_9 = _S3440 * _S3500;
        eps2d_11 = _S3506;
    }
    else
    {
        k_9 = 0.0f;
        eps2d_11 = _S3500;
    }
    float _S3507 = invdet_13 * _S3503.rows[int(1)].y;
    float _S3508 = - (invdet_13 * _S3503.rows[int(1)].x);
    float _S3509 = - (invdet_13 * _S3503.rows[int(0)].y);
    float _S3510 = invdet_13 * _S3503.rows[int(0)].x;
    float _S3511 = - ((_S3468 * _S3503.rows[int(1)].y + _S3478 * _S3503.rows[int(1)].x + _S3477 * _S3503.rows[int(0)].y + _S3470 * _S3503.rows[int(0)].x) / _S3474);
    DiffPair_float_0 _S3512;
    (&_S3512)->primal_0 = _S3475;
    (&_S3512)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3512, k_9);
    DiffPair_float_0 _S3513;
    (&_S3513)->primal_0 = 0.0f;
    (&_S3513)->differential_0 = 0.0f;
    DiffPair_float_0 _S3514;
    (&_S3514)->primal_0 = _S3473;
    (&_S3514)->differential_0 = 0.0f;
    _d_max_0(&_S3513, &_S3514, _S3512.differential_0);
    float _S3515 = _S3514.differential_0 / _S3474;
    float s_diff_det_orig_T_5 = det_blur_11 * _S3515;
    float _S3516 = det_orig_11 * - _S3515 + _S3511;
    float _S3517 = - _S3516;
    float _S3518 = _S3468 * _S3516;
    float _S3519 = _S3470 * _S3516;
    Matrix<float, 2, 2>  _S3520 = _S3501;
    _S3520[int(1)] = _S3504;
    _S3520[int(0)] = _S3505;
    float _S3521 = _S3519 + _S3520.rows[int(0)].x + _S3507;
    float _S3522 = _S3517 + - s_diff_det_orig_T_5;
    float _S3523 = covar2d_8.rows[int(0)].y * _S3522 + _S3508;
    float _S3524 = covar2d_8.rows[int(1)].x * _S3522 + _S3509;
    float _S3525 = covar2d_8.rows[int(1)].y * s_diff_det_orig_T_5;
    float _S3526 = _S3518 + _S3520.rows[int(1)].y + _S3510 + covar2d_8.rows[int(0)].x * s_diff_det_orig_T_5;
    float2  _S3527 = _S3346;
    *&((&_S3527)->x) = _S3523;
    *&((&_S3527)->y) = _S3526;
    float _S3528 = _S3521 + _S3525;
    float2  _S3529 = _S3346;
    *&((&_S3529)->y) = _S3524;
    *&((&_S3529)->x) = _S3528;
    Matrix<float, 2, 2>  _S3530 = _S3501;
    _S3530[int(1)] = _S3527;
    _S3530[int(0)] = _S3529;
    Matrix<float, 3, 3>  _S3531 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3532;
    (&_S3532)->primal_0 = R_11;
    (&_S3532)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3533;
    (&_S3533)->primal_0 = _S3366;
    (&_S3533)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3532, &_S3533, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3534;
    (&_S3534)->primal_0 = R_11;
    (&_S3534)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3535;
    (&_S3535)->primal_0 = _S3363;
    (&_S3535)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3534, &_S3535, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3536;
    (&_S3536)->primal_0 = R_11;
    (&_S3536)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3537;
    (&_S3537)->primal_0 = _S3360;
    (&_S3537)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3536, &_S3537, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3538;
    (&_S3538)->primal_0 = R_11;
    (&_S3538)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3539;
    (&_S3539)->primal_0 = _S3365;
    (&_S3539)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3538, &_S3539, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3540;
    (&_S3540)->primal_0 = R_11;
    (&_S3540)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3541;
    (&_S3541)->primal_0 = _S3362;
    (&_S3541)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3540, &_S3541, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3542;
    (&_S3542)->primal_0 = R_11;
    (&_S3542)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3543;
    (&_S3543)->primal_0 = _S3359;
    (&_S3543)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3542, &_S3543, _S3348);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3544;
    (&_S3544)->primal_0 = R_11;
    (&_S3544)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3545;
    (&_S3545)->primal_0 = _S3378.p_0[0U];
    (&_S3545)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3544, &_S3545, _S3348);
    float3  _S3546 = - _S3533.differential_0 + _S3539.differential_0;
    float3  _S3547 = _S3444 * _S3546;
    float3  _S3548 = _S3352.rows[2U] * _S3546;
    float _S3549 = _S3357 * (_S3548.x + _S3548.y + _S3548.z);
    float3  _S3550 = - _S3535.differential_0 + _S3541.differential_0;
    float3  _S3551 = _S3443 * _S3550;
    float3  _S3552 = _S3352.rows[1U] * _S3550;
    float _S3553 = _S3357 * (_S3552.x + _S3552.y + _S3552.z);
    float3  _S3554 = - _S3537.differential_0 + _S3543.differential_0;
    float3  _S3555 = _S3442 * _S3554;
    float3  _S3556 = _S3352.rows[0U] * _S3554;
    float _S3557 = _S3357 * (_S3556.x + _S3556.y + _S3556.z);
    Matrix<float, 3, 3>  _S3558 = _S3531;
    _S3558[2U] = _S3547;
    _S3558[1U] = _S3551;
    _S3558[0U] = _S3555;
    Matrix<float, 3, 3>  _S3559 = transpose_3(transpose_3(_S3558));
    float _S3560 = 2.0f * - _S3559.rows[int(2)].z;
    float _S3561 = 2.0f * _S3559.rows[int(2)].y;
    float _S3562 = 2.0f * _S3559.rows[int(2)].x;
    float _S3563 = 2.0f * _S3559.rows[int(1)].z;
    float _S3564 = 2.0f * - _S3559.rows[int(1)].y;
    float _S3565 = 2.0f * _S3559.rows[int(1)].x;
    float _S3566 = 2.0f * _S3559.rows[int(0)].z;
    float _S3567 = 2.0f * _S3559.rows[int(0)].y;
    float _S3568 = 2.0f * - _S3559.rows[int(0)].x;
    float _S3569 = - _S3565 + _S3567;
    float _S3570 = _S3562 + - _S3566;
    float _S3571 = - _S3561 + _S3563;
    float _S3572 = _S3561 + _S3563;
    float _S3573 = _S3562 + _S3566;
    float _S3574 = _S3565 + _S3567;
    float _S3575 = _S3350.w * (_S3564 + _S3568);
    float _S3576 = _S3350.z * (_S3560 + _S3568);
    float _S3577 = _S3350.y * (_S3560 + _S3564);
    float _S3578 = _S3350.x * _S3569 + _S3350.z * _S3572 + _S3350.y * _S3573 + _S3575 + _S3575;
    float _S3579 = _S3350.x * _S3570 + _S3350.w * _S3572 + _S3350.y * _S3574 + _S3576 + _S3576;
    float _S3580 = _S3350.x * _S3571 + _S3350.w * _S3573 + _S3350.z * _S3574 + _S3577 + _S3577;
    float _S3581 = _S3350.w * _S3569 + _S3350.z * _S3570 + _S3350.y * _S3571;
    float4  _S3582 = make_float4 (0.0f);
    float4  _S3583 = _S3582;
    *&((&_S3583)->w) = _S3578;
    *&((&_S3583)->z) = _S3579;
    *&((&_S3583)->y) = _S3580;
    *&((&_S3583)->x) = _S3581;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S3584;
    (&_S3584)->primal_0 = quat_11;
    (&_S3584)->differential_0 = _S3582;
    s_bwd_normalize_impl_0(&_S3584, _S3583);
    float3  _S3585 = _S3348;
    *&((&_S3585)->z) = _S3549;
    *&((&_S3585)->y) = _S3553;
    *&((&_S3585)->x) = _S3557;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3586;
    (&_S3586)->primal_0 = scale_11;
    (&_S3586)->differential_0 = _S3348;
    s_bwd_prop_exp_1(&_S3586, _S3585);
    float _S3587 = - (eps2d_11 / _S3441);
    DiffPair_float_0 _S3588;
    (&_S3588)->primal_0 = _S3438;
    (&_S3588)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3588, _S3587);
    float _S3589 = - _S3588.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3590;
    (&_S3590)->primal_0 = mean_c_9;
    (&_S3590)->differential_0 = _S3348;
    s_bwd_length_impl_0(&_S3590, 0.0f);
    float3  _S3591 = _S3590.differential_0 + make_float3 (0.0f, 0.0f, v_depth_5);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3592;
    (&_S3592)->primal_0 = R_11;
    (&_S3592)->differential_0 = _S3531;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3593;
    (&_S3593)->primal_0 = mean_12;
    (&_S3593)->differential_0 = _S3348;
    s_bwd_prop_mul_3(&_S3592, &_S3593, _S3591);
    Matrix<float, 3, 3>  _S3594 = _S3532.differential_0 + _S3534.differential_0 + _S3536.differential_0 + _S3538.differential_0 + _S3540.differential_0 + _S3542.differential_0 + _S3544.differential_0 + _S3592.differential_0;
    float _S3595 = _S3589 + _S3490.differential_0;
    float3  _S3596 = _S3586.differential_0 + _S3489.differential_0;
    *v_mean_5 = *v_mean_5 + (_S3533.differential_0 + _S3539.differential_0 + _S3535.differential_0 + _S3541.differential_0 + _S3537.differential_0 + _S3543.differential_0 + _S3593.differential_0 + _S3488.differential_0);
    *v_quat_5 = *v_quat_5 + _S3584.differential_0;
    *v_scale_5 = *v_scale_5 + _S3596;
    *v_in_opacity_5 = *v_in_opacity_5 + _S3595;
    *v_R_5 = *v_R_5 + _S3594;
    *v_t_5 = *v_t_5 + _S3591;
    return;
}

