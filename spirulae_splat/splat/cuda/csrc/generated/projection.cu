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

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_1, float dOut_1)
{
    float _S4 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_1).primal_0)))))) * dOut_1;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S4;
    return;
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S5 = -0.5f / ((*dpx_2).primal_0 * (F32_sqrt(((*dpx_2).primal_0)))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S5;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_0)
{
    Matrix<float, 3, 3>  result_0;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(3))
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

inline __device__ Matrix<float, 3, 3>  quat_to_rotmat(float4  quat_0)
{
    float x_3 = quat_0.y;
    float inv_norm_0 = (F32_rsqrt((x_3 * x_3 + quat_0.z * quat_0.z + quat_0.w * quat_0.w + quat_0.x * quat_0.x)));
    float x_4 = quat_0.y * inv_norm_0;
    float y_0 = quat_0.z * inv_norm_0;
    float z_0 = quat_0.w * inv_norm_0;
    float w_0 = quat_0.x * inv_norm_0;
    float x2_0 = x_4 * x_4;
    float y2_0 = y_0 * y_0;
    float z2_0 = z_0 * z_0;
    float xy_0 = x_4 * y_0;
    float xz_0 = x_4 * z_0;
    float yz_0 = y_0 * z_0;
    float wx_0 = w_0 * x_4;
    float wy_0 = w_0 * y_0;
    float wz_0 = w_0 * z_0;
    return transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0)));
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_3)
{
    float _S6 = (*left_0).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_3.x;
    float sum_0 = _S6 + (*left_0).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_3.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_3.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S7 = (*left_0).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_3.x;
    float sum_2 = _S7 + (*left_0).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_3.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_3.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S8 = (*left_0).primal_0.rows[int(0)].z * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_3.x;
    float sum_4 = _S8 + (*left_0).primal_0.rows[int(1)].z * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_3.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_3.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_3;
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
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_0], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_3, i_0) = sum_6;
        i_0 = i_0 + int(1);
    }
    return result_3;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ void mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, Matrix<float, 3, 3>  dOut_4)
{
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].z;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

struct DiffPair_matrixx3Cfloatx2C2x2C3x3E_0
{
    Matrix<float, 2, 3>  primal_0;
    Matrix<float, 2, 3>  differential_0;
};

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_3, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_3, Matrix<float, 2, 3>  dOut_5)
{
    Matrix<float, 2, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = 0.0f;
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
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    left_3->primal_0 = (*left_3).primal_0;
    left_3->differential_0 = left_d_result_2;
    right_3->primal_0 = (*right_3).primal_0;
    right_3->differential_0 = right_d_result_2;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C2x3E_0
{
    Matrix<float, 3, 2>  primal_0;
    Matrix<float, 3, 2>  differential_0;
};

inline __device__ void mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * right_4, Matrix<float, 2, 2>  dOut_6)
{
    Matrix<float, 2, 3>  left_d_result_3;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = 0.0f;
    Matrix<float, 3, 2>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(0)].z * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(0)].z * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(1)].z * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(1)].z * dOut_6.rows[int(1)].y;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_3;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_3;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_4(Matrix<float, 3, 3>  left_5, Matrix<float, 3, 3>  right_5)
{
    Matrix<float, 3, 3>  result_4;
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
            int i_1 = int(0);
            float sum_8 = 0.0f;
            for(;;)
            {
                if(i_1 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_9 = sum_8 + _slang_vector_get_element(left_5.rows[r_3], i_1) * _slang_vector_get_element(right_5.rows[i_1], c_3);
                i_1 = i_1 + int(1);
                sum_8 = sum_9;
            }
            *_slang_vector_get_element_ptr(((&result_4)->rows + (r_3)), c_3) = sum_8;
            c_3 = c_3 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_4;
}

inline __device__ Matrix<float, 2, 3>  mul_5(Matrix<float, 2, 3>  left_6, Matrix<float, 3, 3>  right_6)
{
    Matrix<float, 2, 3>  result_5;
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
            int i_2 = int(0);
            float sum_10 = 0.0f;
            for(;;)
            {
                if(i_2 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_11 = sum_10 + _slang_vector_get_element(left_6.rows[r_4], i_2) * _slang_vector_get_element(right_6.rows[i_2], c_4);
                i_2 = i_2 + int(1);
                sum_10 = sum_11;
            }
            *_slang_vector_get_element_ptr(((&result_5)->rows + (r_4)), c_4) = sum_10;
            c_4 = c_4 + int(1);
        }
        r_4 = r_4 + int(1);
    }
    return result_5;
}

inline __device__ Matrix<float, 2, 2>  mul_6(Matrix<float, 2, 3>  left_7, Matrix<float, 3, 2>  right_7)
{
    Matrix<float, 2, 2>  result_6;
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
            int i_3 = int(0);
            float sum_12 = 0.0f;
            for(;;)
            {
                if(i_3 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_13 = sum_12 + _slang_vector_get_element(left_7.rows[r_5], i_3) * _slang_vector_get_element(right_7.rows[i_3], c_5);
                i_3 = i_3 + int(1);
                sum_12 = sum_13;
            }
            *_slang_vector_get_element_ptr(((&result_6)->rows + (r_5)), c_5) = sum_12;
            c_5 = c_5 + int(1);
        }
        r_5 = r_5 + int(1);
    }
    return result_6;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_4(mul_4(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_1, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_5 = quat_1.y;
    float inv_norm_1 = (F32_rsqrt((x_5 * x_5 + quat_1.z * quat_1.z + quat_1.w * quat_1.w + quat_1.x * quat_1.x)));
    float x_6 = quat_1.y * inv_norm_1;
    float y_1 = quat_1.z * inv_norm_1;
    float z_1 = quat_1.w * inv_norm_1;
    float w_1 = quat_1.x * inv_norm_1;
    float x2_1 = x_6 * x_6;
    float y2_1 = y_1 * y_1;
    float z2_1 = z_1 * z_1;
    float xy_1 = x_6 * y_1;
    float xz_1 = x_6 * z_1;
    float yz_1 = y_1 * z_1;
    float wx_1 = w_1 * x_6;
    float wy_1 = w_1 * y_1;
    float wz_1 = w_1 * z_1;
    Matrix<float, 3, 3>  M_0 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_4(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_2, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_7 = quat_2.y;
    float inv_norm_2 = (F32_rsqrt((x_7 * x_7 + quat_2.z * quat_2.z + quat_2.w * quat_2.w + quat_2.x * quat_2.x)));
    float x_8 = quat_2.y * inv_norm_2;
    float y_2 = quat_2.z * inv_norm_2;
    float z_2 = quat_2.w * inv_norm_2;
    float w_2 = quat_2.x * inv_norm_2;
    float x2_2 = x_8 * x_8;
    float y2_2 = y_2 * y_2;
    float z2_2 = z_2 * z_2;
    float xy_2 = x_8 * y_2;
    float xz_2 = x_8 * z_2;
    float yz_2 = y_2 * z_2;
    float wx_2 = w_2 * x_8;
    float wy_2 = w_2 * y_2;
    float wz_2 = w_2 * z_2;
    *M_1 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_1, float dOut_7)
{
    DiffPair_float_0 _S9 = *dpx_3;
    float _S10;
    if(((*dpx_3).primal_0) < ((*dpy_1).primal_0))
    {
        _S10 = dOut_7;
    }
    else
    {
        if(((*dpx_3).primal_0) > ((*dpy_1).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_7;
        }
    }
    dpx_3->primal_0 = _S9.primal_0;
    dpx_3->differential_0 = _S10;
    DiffPair_float_0 _S11 = *dpy_1;
    if(((*dpy_1).primal_0) < (_S9.primal_0))
    {
        _S10 = dOut_7;
    }
    else
    {
        if(((*dpy_1).primal_0) > ((*dpx_3).primal_0))
        {
            _S10 = 0.0f;
        }
        else
        {
            _S10 = 0.5f * dOut_7;
        }
    }
    dpy_1->primal_0 = _S11.primal_0;
    dpy_1->differential_0 = _S10;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S12 = float(width_0);
    float _S13 = float(height_0);
    float _S14 = 0.30000001192092896f * (0.5f * _S12 / fx_0);
    float _S15 = 0.30000001192092896f * (0.5f * _S13 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_0 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S12 - cx_0) / fx_0 + _S14), ((F32_max((- (cx_0 / fx_0 + _S14)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S13 - cy_0) / fy_0 + _S15), ((F32_max((- (cy_0 / fy_0 + _S15)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_0, cov3d_0), transpose_1(J_0));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float fx_1, float fy_1, float cx_1, float cy_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float rz_1 = 1.0f / mean3d_1.z;
    float rz2_1 = rz_1 * rz_1;
    Matrix<float, 2, 3>  J_1 = makeMatrix<float, 2, 3> (fx_1 * rz_1, 0.0f, - fx_1 * mean3d_1.x * rz2_1, 0.0f, fy_1 * rz_1, - fy_1 * mean3d_1.y * rz2_1);
    *cov2d_1 = mul_6(mul_5(J_1, cov3d_1), transpose_1(J_1));
    *mean2d_1 = make_float2 (fx_1 * mean3d_1.x * rz_1 + cx_1, fy_1 * mean3d_1.y * rz_1 + cy_1);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_4, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_2, float dOut_8)
{
    float2  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_8;
    float2  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_4).primal_0.x * dOut_8;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_8;
    *&((&y_d_result_0)->y) = (*dpx_4).primal_0.y * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
    return;
}

inline __device__ float dot_0(float2  x_9, float2  y_3)
{
    int i_4 = int(0);
    float result_7 = 0.0f;
    for(;;)
    {
        if(i_4 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_8 = result_7 + _slang_vector_get_element(x_9, i_4) * _slang_vector_get_element(y_3, i_4);
        i_4 = i_4 + int(1);
        result_7 = result_8;
    }
    return result_7;
}

inline __device__ float length_0(float2  x_10)
{
    return (F32_sqrt((dot_0(x_10, x_10))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_3, DiffPair_float_0 * dpx_5, float dOut_9)
{
    DiffPair_float_0 _S16 = *dpx_5;
    float _S17 = - (*dpy_3).primal_0 / ((*dpx_5).primal_0 * (*dpx_5).primal_0 + (*dpy_3).primal_0 * (*dpy_3).primal_0) * dOut_9;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S17;
    float _S18 = _S16.primal_0 / (_S16.primal_0 * _S16.primal_0 + (*dpy_3).primal_0 * (*dpy_3).primal_0) * dOut_9;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = _S18;
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    float xy_len_0 = length_0(make_float2 (mean3d_2.x, mean3d_2.y)) + 1.00000001168609742e-07f;
    float theta_0 = (F32_atan2((xy_len_0), (mean3d_2.z + 1.00000001168609742e-07f)));
    *mean2d_2 = make_float2 (mean3d_2.x * fx_2 * theta_0 / xy_len_0 + cx_2, mean3d_2.y * fy_2 * theta_0 / xy_len_0 + cy_2);
    float x2_3 = mean3d_2.x * mean3d_2.x + 1.00000001168609742e-07f;
    float y2_3 = mean3d_2.y * mean3d_2.y;
    float xy_3 = mean3d_2.x * mean3d_2.y;
    float x2y2_0 = x2_3 + y2_3;
    float x2y2z2_inv_0 = 1.0f / (x2y2_0 + mean3d_2.z * mean3d_2.z);
    float b_0 = (F32_atan2((xy_len_0), (mean3d_2.z))) / xy_len_0 / x2y2_0;
    float a_0 = mean3d_2.z * x2y2z2_inv_0 / x2y2_0;
    float _S19 = a_0 - b_0;
    Matrix<float, 2, 3>  J_2 = makeMatrix<float, 2, 3> (fx_2 * (x2_3 * a_0 + y2_3 * b_0), fx_2 * xy_3 * _S19, - fx_2 * mean3d_2.x * x2y2z2_inv_0, fy_2 * xy_3 * _S19, fy_2 * (y2_3 * a_0 + x2_3 * b_0), - fy_2 * mean3d_2.y * x2y2z2_inv_0);
    *cov2d_2 = mul_6(mul_5(J_2, cov3d_2), transpose_1(J_2));
    return;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_3, float fy_3, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    Matrix<float, 2, 3>  J_3 = makeMatrix<float, 2, 3> (fx_3, 0.0f, 0.0f, 0.0f, fy_3, 0.0f);
    *cov2d_3 = mul_6(mul_5(J_3, cov3d_3), transpose_1(J_3));
    *mean2d_3 = make_float2 (fx_3 * mean3d_3.x + cx_3, fy_3 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *&((covar_1->rows + (int(0)))->x) = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    float _S20 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S20;
    float det_blur_0 = *&((covar_1->rows + (int(0)))->x) * _S20 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_6, float dOut_10)
{
    float _S21 = (F32_exp(((*dpx_6).primal_0))) * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S21;
    return;
}

inline __device__ float3  exp_0(float3  x_11)
{
    float3  result_9;
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
        *_slang_vector_get_element_ptr(&result_9, i_5) = (F32_exp((_slang_vector_get_element(x_11, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_9;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_7, float3  dOut_11)
{
    float3  _S22 = exp_0((*dpx_7).primal_0) * dOut_11;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S22;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_8, float dOut_12)
{
    float _S23 = 1.0f / (*dpx_8).primal_0 * dOut_12;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S23;
    return;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_4, float fy_4, float cx_4, float cy_4, uint image_width_0, uint image_height_0, float eps2d_1, float near_plane_0, float far_plane_0, float radius_clip_0, int2  * radii_0, float * depth_0, float2  * mean2d_4, float3  * conic_0, float * opacity_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_2, mean_0) + t_1;
        float _S24 = mean_c_0.z;
        bool _S25;
        if(_S24 < near_plane_0)
        {
            _S25 = true;
        }
        else
        {
            _S25 = _S24 > far_plane_0;
        }
        if(_S25)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        float3  _S26 = exp_0(scale_2);
        float x_12 = quat_3.y;
        float inv_norm_3 = (F32_rsqrt((x_12 * x_12 + quat_3.z * quat_3.z + quat_3.w * quat_3.w + quat_3.x * quat_3.x)));
        float x_13 = quat_3.y * inv_norm_3;
        float y_4 = quat_3.z * inv_norm_3;
        float z_3 = quat_3.w * inv_norm_3;
        float w_3 = quat_3.x * inv_norm_3;
        float x2_4 = x_13 * x_13;
        float y2_4 = y_4 * y_4;
        float z2_3 = z_3 * z_3;
        float xy_4 = x_13 * y_4;
        float xz_3 = x_13 * z_3;
        float yz_3 = y_4 * z_3;
        float wx_3 = w_3 * x_13;
        float wy_3 = w_3 * y_4;
        float wz_3 = w_3 * z_3;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_3), 2.0f * (xy_4 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_4 - wz_3), 1.0f - 2.0f * (x2_4 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S26.x, 0.0f, 0.0f, 0.0f, _S26.y, 0.0f, 0.0f, 0.0f, _S26.z));
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_2, mul_4(M_2, transpose_0(M_2))), transpose_0(R_2));
        Matrix<float, 2, 2>  covar2d_0;
        float _S27 = float(image_width_0);
        float _S28 = float(image_height_0);
        float _S29 = 0.30000001192092896f * (0.5f * _S27 / fx_4);
        float _S30 = 0.30000001192092896f * (0.5f * _S28 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S27 - cx_4) / fx_4 + _S29), ((F32_max((- (cx_4 / fx_4 + _S29)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S28 - cy_4) / fy_4 + _S30), ((F32_max((- (cy_4 / fy_4 + _S30)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_4, covar_c_0), transpose_1(J_4));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        *&(((&covar2d_0)->rows + (int(0)))->x) = *&(((&covar2d_0)->rows + (int(0)))->x) + eps2d_1;
        float _S31 = *&(((&covar2d_0)->rows + (int(1)))->y) + eps2d_1;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S31;
        float det_blur_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * _S31 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S32 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
        *opacity_0 = 1.0f / (1.0f + (F32_exp((- in_opacity_0))));
        if(antialiased_0)
        {
            *opacity_0 = *opacity_0 * compensation_1;
        }
        if((*opacity_0) < 0.00392156885936856f)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        float extend_0 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_0 / 0.00392156885936856f)))))))));
        float radius_x_0 = (F32_ceil((extend_0 * (F32_sqrt((covar2d_0[int(0)].x))))));
        float radius_y_0 = (F32_ceil((extend_0 * (F32_sqrt((covar2d_0[int(1)].y))))));
        if(radius_x_0 <= radius_clip_0)
        {
            _S25 = radius_y_0 <= radius_clip_0;
        }
        else
        {
            _S25 = false;
        }
        if(_S25)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_4).x + radius_x_0) <= 0.0f)
        {
            _S25 = true;
        }
        else
        {
            _S25 = ((*mean2d_4).x - radius_x_0) >= _S27;
        }
        if(_S25)
        {
            _S25 = true;
        }
        else
        {
            _S25 = ((*mean2d_4).y + radius_y_0) <= 0.0f;
        }
        if(_S25)
        {
            _S25 = true;
        }
        else
        {
            _S25 = ((*mean2d_4).y - radius_y_0) >= _S28;
        }
        if(_S25)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        *radii_0 = make_int2 (int(radius_x_0), int(radius_y_0));
        *depth_0 = _S24;
        *conic_0 = make_float3 (_S32.rows[int(0)].x, _S32.rows[int(0)].y, _S32.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_5, float fy_5, float cx_5, float cy_5, uint image_width_1, uint image_height_1, float eps2d_2, float near_plane_1, float far_plane_1, float radius_clip_1, int2  * radii_1, float * depth_1, float2  * mean2d_5, float3  * conic_1, float * opacity_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_3, mean_1) + t_2;
        float _S33 = mean_c_1.z;
        bool _S34;
        if(_S33 < near_plane_1)
        {
            _S34 = true;
        }
        else
        {
            _S34 = _S33 > far_plane_1;
        }
        if(_S34)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        float3  _S35 = exp_0(scale_3);
        float x_14 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_14 * x_14 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_15 = quat_4.y * inv_norm_4;
        float y_5 = quat_4.z * inv_norm_4;
        float z_4 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_5 = x_15 * x_15;
        float y2_5 = y_5 * y_5;
        float z2_4 = z_4 * z_4;
        float xy_5 = x_15 * y_5;
        float xz_4 = x_15 * z_4;
        float yz_4 = y_5 * z_4;
        float wx_4 = w_4 * x_15;
        float wy_4 = w_4 * y_5;
        float wz_4 = w_4 * z_4;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_4), 2.0f * (xy_5 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_5 - wz_4), 1.0f - 2.0f * (x2_5 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S35.x, 0.0f, 0.0f, 0.0f, _S35.y, 0.0f, 0.0f, 0.0f, _S35.z));
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_3, mul_4(M_3, transpose_0(M_3))), transpose_0(R_3));
        Matrix<float, 2, 2>  covar2d_1;
        float xy_len_1 = length_0(make_float2 (mean_c_1.x, mean_c_1.y)) + 1.00000001168609742e-07f;
        float theta_1 = (F32_atan2((xy_len_1), (mean_c_1.z + 1.00000001168609742e-07f)));
        *mean2d_5 = make_float2 (mean_c_1.x * fx_5 * theta_1 / xy_len_1 + cx_5, mean_c_1.y * fy_5 * theta_1 / xy_len_1 + cy_5);
        float x2_6 = mean_c_1.x * mean_c_1.x + 1.00000001168609742e-07f;
        float y2_6 = mean_c_1.y * mean_c_1.y;
        float xy_6 = mean_c_1.x * mean_c_1.y;
        float x2y2_1 = x2_6 + y2_6;
        float x2y2z2_inv_1 = 1.0f / (x2y2_1 + mean_c_1.z * mean_c_1.z);
        float b_1 = (F32_atan2((xy_len_1), (mean_c_1.z))) / xy_len_1 / x2y2_1;
        float a_1 = mean_c_1.z * x2y2z2_inv_1 / x2y2_1;
        float _S36 = a_1 - b_1;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_5 * (x2_6 * a_1 + y2_6 * b_1), fx_5 * xy_6 * _S36, - fx_5 * mean_c_1.x * x2y2z2_inv_1, fy_5 * xy_6 * _S36, fy_5 * (y2_6 * a_1 + x2_6 * b_1), - fy_5 * mean_c_1.y * x2y2z2_inv_1);
        covar2d_1 = mul_6(mul_5(J_5, covar_c_1), transpose_1(J_5));
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        *&(((&covar2d_1)->rows + (int(0)))->x) = *&(((&covar2d_1)->rows + (int(0)))->x) + eps2d_2;
        float _S37 = *&(((&covar2d_1)->rows + (int(1)))->y) + eps2d_2;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S37;
        float det_blur_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * _S37 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S38 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
        *opacity_1 = 1.0f / (1.0f + (F32_exp((- in_opacity_1))));
        if(antialiased_1)
        {
            *opacity_1 = *opacity_1 * compensation_2;
        }
        if((*opacity_1) < 0.00392156885936856f)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        float extend_1 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_1 / 0.00392156885936856f)))))))));
        float radius_x_1 = (F32_ceil((extend_1 * (F32_sqrt((covar2d_1[int(0)].x))))));
        float radius_y_1 = (F32_ceil((extend_1 * (F32_sqrt((covar2d_1[int(1)].y))))));
        if(radius_x_1 <= radius_clip_1)
        {
            _S34 = radius_y_1 <= radius_clip_1;
        }
        else
        {
            _S34 = false;
        }
        if(_S34)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_5).x + radius_x_1) <= 0.0f)
        {
            _S34 = true;
        }
        else
        {
            _S34 = ((*mean2d_5).x - radius_x_1) >= float(image_width_1);
        }
        if(_S34)
        {
            _S34 = true;
        }
        else
        {
            _S34 = ((*mean2d_5).y + radius_y_1) <= 0.0f;
        }
        if(_S34)
        {
            _S34 = true;
        }
        else
        {
            _S34 = ((*mean2d_5).y - radius_y_1) >= float(image_height_1);
        }
        if(_S34)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        *radii_1 = make_int2 (int(radius_x_1), int(radius_y_1));
        *depth_1 = _S33;
        *conic_1 = make_float3 (_S38.rows[int(0)].x, _S38.rows[int(0)].y, _S38.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_6, float fy_6, float cx_6, float cy_6, uint image_width_2, uint image_height_2, float eps2d_3, float near_plane_2, float far_plane_2, float radius_clip_2, int2  * radii_2, float * depth_2, float2  * mean2d_6, float3  * conic_2, float * opacity_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_4, mean_2) + t_3;
        float _S39 = mean_c_2.z;
        bool _S40;
        if(_S39 < near_plane_2)
        {
            _S40 = true;
        }
        else
        {
            _S40 = _S39 > far_plane_2;
        }
        if(_S40)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        float3  _S41 = exp_0(scale_4);
        float x_16 = quat_5.y;
        float inv_norm_5 = (F32_rsqrt((x_16 * x_16 + quat_5.z * quat_5.z + quat_5.w * quat_5.w + quat_5.x * quat_5.x)));
        float x_17 = quat_5.y * inv_norm_5;
        float y_6 = quat_5.z * inv_norm_5;
        float z_5 = quat_5.w * inv_norm_5;
        float w_5 = quat_5.x * inv_norm_5;
        float x2_7 = x_17 * x_17;
        float y2_7 = y_6 * y_6;
        float z2_5 = z_5 * z_5;
        float xy_7 = x_17 * y_6;
        float xz_5 = x_17 * z_5;
        float yz_5 = y_6 * z_5;
        float wx_5 = w_5 * x_17;
        float wy_5 = w_5 * y_6;
        float wz_5 = w_5 * z_5;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_5), 2.0f * (xy_7 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_7 - wz_5), 1.0f - 2.0f * (x2_7 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S41.x, 0.0f, 0.0f, 0.0f, _S41.y, 0.0f, 0.0f, 0.0f, _S41.z));
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_4, mul_4(M_4, transpose_0(M_4))), transpose_0(R_4));
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        *&(((&covar2d_2)->rows + (int(0)))->x) = *&(((&covar2d_2)->rows + (int(0)))->x) + eps2d_3;
        float _S42 = *&(((&covar2d_2)->rows + (int(1)))->y) + eps2d_3;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S42;
        float det_blur_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * _S42 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S43 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
        *opacity_2 = 1.0f / (1.0f + (F32_exp((- in_opacity_2))));
        if(antialiased_2)
        {
            *opacity_2 = *opacity_2 * compensation_3;
        }
        if((*opacity_2) < 0.00392156885936856f)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        float extend_2 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_2 / 0.00392156885936856f)))))))));
        float radius_x_2 = (F32_ceil((extend_2 * (F32_sqrt((covar2d_2[int(0)].x))))));
        float radius_y_2 = (F32_ceil((extend_2 * (F32_sqrt((covar2d_2[int(1)].y))))));
        if(radius_x_2 <= radius_clip_2)
        {
            _S40 = radius_y_2 <= radius_clip_2;
        }
        else
        {
            _S40 = false;
        }
        if(_S40)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_6).x + radius_x_2) <= 0.0f)
        {
            _S40 = true;
        }
        else
        {
            _S40 = ((*mean2d_6).x - radius_x_2) >= float(image_width_2);
        }
        if(_S40)
        {
            _S40 = true;
        }
        else
        {
            _S40 = ((*mean2d_6).y + radius_y_2) <= 0.0f;
        }
        if(_S40)
        {
            _S40 = true;
        }
        else
        {
            _S40 = ((*mean2d_6).y - radius_y_2) >= float(image_height_2);
        }
        if(_S40)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        *radii_2 = make_int2 (int(radius_x_2), int(radius_y_2));
        *depth_2 = _S39;
        *conic_2 = make_float3 (_S43.rows[int(0)].x, _S43.rows[int(0)].y, _S43.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_7, float fy_7, float cx_7, float cy_7, uint image_width_3, uint image_height_3, float eps2d_4, float near_plane_3, float far_plane_3, float radius_clip_3, int2  * radii_3, float * depth_3, float2  * mean2d_7, float3  * conic_3, float * opacity_3)
{
    float3  mean_c_3 = mul_0(R_5, mean_3) + t_4;
    float3  _S44 = exp_0(scale_5);
    float x_18 = quat_6.y;
    float inv_norm_6 = (F32_rsqrt((x_18 * x_18 + quat_6.z * quat_6.z + quat_6.w * quat_6.w + quat_6.x * quat_6.x)));
    float x_19 = quat_6.y * inv_norm_6;
    float y_7 = quat_6.z * inv_norm_6;
    float z_6 = quat_6.w * inv_norm_6;
    float w_6 = quat_6.x * inv_norm_6;
    float x2_8 = x_19 * x_19;
    float y2_8 = y_7 * y_7;
    float z2_6 = z_6 * z_6;
    float xy_8 = x_19 * y_7;
    float xz_6 = x_19 * z_6;
    float yz_6 = y_7 * z_6;
    float wx_6 = w_6 * x_19;
    float wy_6 = w_6 * y_7;
    float wz_6 = w_6 * z_6;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_6), 2.0f * (xy_8 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_8 - wz_6), 1.0f - 2.0f * (x2_8 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S44.x, 0.0f, 0.0f, 0.0f, _S44.y, 0.0f, 0.0f, 0.0f, _S44.z));
    float _S45 = float(image_width_3);
    float _S46 = float(image_height_3);
    float _S47 = 0.30000001192092896f * (0.5f * _S45 / fx_7);
    float _S48 = 0.30000001192092896f * (0.5f * _S46 / fy_7);
    float rz_3 = 1.0f / mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S45 - cx_7) / fx_7 + _S47), ((F32_max((- (cx_7 / fx_7 + _S47)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S46 - cy_7) / fy_7 + _S48), ((F32_max((- (cy_7 / fy_7 + _S48)), (mean_c_3.y * rz_3))))))) * rz2_3);
    Matrix<float, 2, 2>  covar2d_3 = mul_6(mul_5(J_7, mul_4(mul_4(R_5, mul_4(M_5, transpose_0(M_5))), transpose_0(R_5))), transpose_1(J_7));
    *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
    float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
    *&(((&covar2d_3)->rows + (int(0)))->x) = *&(((&covar2d_3)->rows + (int(0)))->x) + eps2d_4;
    float _S49 = *&(((&covar2d_3)->rows + (int(1)))->y) + eps2d_4;
    *&(((&covar2d_3)->rows + (int(1)))->y) = _S49;
    float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / (*&(((&covar2d_3)->rows + (int(0)))->x) * _S49 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
    Matrix<float, 2, 2>  _S50 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_4, - covar2d_3.rows[int(0)].y * invdet_4, - covar2d_3.rows[int(1)].x * invdet_4, covar2d_3.rows[int(0)].x * invdet_4);
    *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
    if(antialiased_3)
    {
        *opacity_3 = *opacity_3 * compensation_4;
    }
    float extend_3 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
    *radii_3 = make_int2 (int((F32_ceil((extend_3 * (F32_sqrt((covar2d_3[int(0)].x))))))), int((F32_ceil((extend_3 * (F32_sqrt((covar2d_3[int(1)].y))))))));
    *depth_3 = mean_c_3.z;
    *conic_3 = make_float3 (_S50.rows[int(0)].x, _S50.rows[int(0)].y, _S50.rows[int(1)].y);
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_8, float fy_8, float cx_8, float cy_8, uint image_width_4, uint image_height_4, float eps2d_5, float near_plane_4, float far_plane_4, float radius_clip_4, int2  * radii_4, float * depth_4, float2  * mean2d_8, float3  * conic_4, float * opacity_4)
{
    float3  mean_c_4 = mul_0(R_6, mean_4) + t_5;
    float3  _S51 = exp_0(scale_6);
    float x_20 = quat_7.y;
    float inv_norm_7 = (F32_rsqrt((x_20 * x_20 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x)));
    float x_21 = quat_7.y * inv_norm_7;
    float y_8 = quat_7.z * inv_norm_7;
    float z_7 = quat_7.w * inv_norm_7;
    float w_7 = quat_7.x * inv_norm_7;
    float x2_9 = x_21 * x_21;
    float y2_9 = y_8 * y_8;
    float z2_7 = z_7 * z_7;
    float xy_9 = x_21 * y_8;
    float xz_7 = x_21 * z_7;
    float yz_7 = y_8 * z_7;
    float wx_7 = w_7 * x_21;
    float wy_7 = w_7 * y_8;
    float wz_7 = w_7 * z_7;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_7), 2.0f * (xy_9 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_9 - wz_7), 1.0f - 2.0f * (x2_9 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S51.x, 0.0f, 0.0f, 0.0f, _S51.y, 0.0f, 0.0f, 0.0f, _S51.z));
    Matrix<float, 3, 3>  covar_c_3 = mul_4(mul_4(R_6, mul_4(M_6, transpose_0(M_6))), transpose_0(R_6));
    float xy_len_2 = length_0(make_float2 (mean_c_4.x, mean_c_4.y)) + 1.00000001168609742e-07f;
    float theta_2 = (F32_atan2((xy_len_2), (mean_c_4.z + 1.00000001168609742e-07f)));
    *mean2d_8 = make_float2 (mean_c_4.x * fx_8 * theta_2 / xy_len_2 + cx_8, mean_c_4.y * fy_8 * theta_2 / xy_len_2 + cy_8);
    float x2_10 = mean_c_4.x * mean_c_4.x + 1.00000001168609742e-07f;
    float y2_10 = mean_c_4.y * mean_c_4.y;
    float xy_10 = mean_c_4.x * mean_c_4.y;
    float x2y2_2 = x2_10 + y2_10;
    float x2y2z2_inv_2 = 1.0f / (x2y2_2 + mean_c_4.z * mean_c_4.z);
    float b_2 = (F32_atan2((xy_len_2), (mean_c_4.z))) / xy_len_2 / x2y2_2;
    float a_2 = mean_c_4.z * x2y2z2_inv_2 / x2y2_2;
    float _S52 = a_2 - b_2;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_8 * (x2_10 * a_2 + y2_10 * b_2), fx_8 * xy_10 * _S52, - fx_8 * mean_c_4.x * x2y2z2_inv_2, fy_8 * xy_10 * _S52, fy_8 * (y2_10 * a_2 + x2_10 * b_2), - fy_8 * mean_c_4.y * x2y2z2_inv_2);
    Matrix<float, 2, 2>  covar2d_4 = mul_6(mul_5(J_8, covar_c_3), transpose_1(J_8));
    float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
    *&(((&covar2d_4)->rows + (int(0)))->x) = *&(((&covar2d_4)->rows + (int(0)))->x) + eps2d_5;
    float _S53 = *&(((&covar2d_4)->rows + (int(1)))->y) + eps2d_5;
    *&(((&covar2d_4)->rows + (int(1)))->y) = _S53;
    float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / (*&(((&covar2d_4)->rows + (int(0)))->x) * _S53 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
    Matrix<float, 2, 2>  _S54 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_5, - covar2d_4.rows[int(0)].y * invdet_5, - covar2d_4.rows[int(1)].x * invdet_5, covar2d_4.rows[int(0)].x * invdet_5);
    *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
    if(antialiased_4)
    {
        *opacity_4 = *opacity_4 * compensation_5;
    }
    float extend_4 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
    *radii_4 = make_int2 (int((F32_ceil((extend_4 * (F32_sqrt((covar2d_4[int(0)].x))))))), int((F32_ceil((extend_4 * (F32_sqrt((covar2d_4[int(1)].y))))))));
    *depth_4 = mean_c_4.z;
    *conic_4 = make_float3 (_S54.rows[int(0)].x, _S54.rows[int(0)].y, _S54.rows[int(1)].y);
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_9, float fy_9, float cx_9, float cy_9, uint image_width_5, uint image_height_5, float eps2d_6, float near_plane_5, float far_plane_5, float radius_clip_5, int2  * radii_5, float * depth_5, float2  * mean2d_9, float3  * conic_5, float * opacity_5)
{
    float3  mean_c_5 = mul_0(R_7, mean_5) + t_6;
    float3  _S55 = exp_0(scale_7);
    float x_22 = quat_8.y;
    float inv_norm_8 = (F32_rsqrt((x_22 * x_22 + quat_8.z * quat_8.z + quat_8.w * quat_8.w + quat_8.x * quat_8.x)));
    float x_23 = quat_8.y * inv_norm_8;
    float y_9 = quat_8.z * inv_norm_8;
    float z_8 = quat_8.w * inv_norm_8;
    float w_8 = quat_8.x * inv_norm_8;
    float x2_11 = x_23 * x_23;
    float y2_11 = y_9 * y_9;
    float z2_8 = z_8 * z_8;
    float xy_11 = x_23 * y_9;
    float xz_8 = x_23 * z_8;
    float yz_8 = y_9 * z_8;
    float wx_8 = w_8 * x_23;
    float wy_8 = w_8 * y_9;
    float wz_8 = w_8 * z_8;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_8), 2.0f * (xy_11 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_11 - wz_8), 1.0f - 2.0f * (x2_11 + z2_8), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S55.x, 0.0f, 0.0f, 0.0f, _S55.y, 0.0f, 0.0f, 0.0f, _S55.z));
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_9, mul_4(mul_4(R_7, mul_4(M_7, transpose_0(M_7))), transpose_0(R_7))), transpose_1(J_9));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x + cx_9, fy_9 * mean_c_5.y + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    *&(((&covar2d_5)->rows + (int(0)))->x) = *&(((&covar2d_5)->rows + (int(0)))->x) + eps2d_6;
    float _S56 = *&(((&covar2d_5)->rows + (int(1)))->y) + eps2d_6;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S56;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (*&(((&covar2d_5)->rows + (int(0)))->x) * _S56 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S57 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_6, - covar2d_5.rows[int(0)].y * invdet_6, - covar2d_5.rows[int(1)].x * invdet_6, covar2d_5.rows[int(0)].x * invdet_6);
    *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
    if(antialiased_5)
    {
        *opacity_5 = *opacity_5 * compensation_6;
    }
    float extend_5 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
    *radii_5 = make_int2 (int((F32_ceil((extend_5 * (F32_sqrt((covar2d_5[int(0)].x))))))), int((F32_ceil((extend_5 * (F32_sqrt((covar2d_5[int(1)].y))))))));
    *depth_5 = mean_c_5.z;
    *conic_5 = make_float3 (_S57.rows[int(0)].x, _S57.rows[int(0)].y, _S57.rows[int(1)].y);
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S58, float3  _S59)
{
    return mul_0(_S58, _S59);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S60)
{
    return exp_0(_S60);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S61)
{
    return (F32_rsqrt((_S61)));
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S62, Matrix<float, 3, 3>  _S63)
{
    return mul_4(_S62, _S63);
}

inline __device__ float s_primal_ctx_max_0(float _S64, float _S65)
{
    return (F32_max((_S64), (_S65)));
}

inline __device__ float s_primal_ctx_min_0(float _S66, float _S67)
{
    return (F32_min((_S66), (_S67)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S68, Matrix<float, 3, 3>  _S69)
{
    return mul_5(_S68, _S69);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S70, Matrix<float, 3, 2>  _S71)
{
    return mul_6(_S70, _S71);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S72)
{
    return (F32_sqrt((_S72)));
}

inline __device__ float s_primal_ctx_exp_1(float _S73)
{
    return (F32_exp((_S73)));
}

inline __device__ float s_primal_ctx_log_0(float _S74)
{
    return (F32_log((_S74)));
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S75, float _S76)
{
    _d_sqrt_0(_S75, _S76);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S77, DiffPair_float_0 * _S78, float _S79)
{
    _d_min_0(_S77, _S78, _S79);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S80, float _S81)
{
    _d_log_0(_S80, _S81);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S82, float _S83)
{
    _d_exp_0(_S82, _S83);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S84, DiffPair_float_0 * _S85, float _S86)
{
    _d_max_0(_S84, _S85, _S86);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S87, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S88, Matrix<float, 2, 2>  _S89)
{
    mul_3(_S87, _S88, _S89);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S90, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S91, Matrix<float, 2, 3>  _S92)
{
    mul_2(_S90, _S91, _S92);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S93, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S94, Matrix<float, 3, 3>  _S95)
{
    mul_1(_S93, _S94, _S95);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S96, float _S97)
{
    _d_rsqrt_0(_S96, _S97);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S98, float3  _S99)
{
    _d_exp_vector_0(_S98, _S99);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S100, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S101, float3  _S102)
{
    _d_mul_0(_S100, _S101, _S102);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_10, float fy_10, float cx_10, float cy_10, uint image_width_6, uint image_height_6, float eps2d_7, float v_depth_0, float2  v_mean2d_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_8, mean_6) + t_7;
    float3  _S103 = s_primal_ctx_exp_0(scale_8);
    float _S104 = quat_9.y;
    float _S105 = _S104 * _S104 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S106 = s_primal_ctx_rsqrt_0(_S105);
    float x_24 = quat_9.y * _S106;
    float y_10 = quat_9.z * _S106;
    float z_9 = quat_9.w * _S106;
    float w_9 = quat_9.x * _S106;
    float x2_12 = x_24 * x_24;
    float y2_12 = y_10 * y_10;
    float z2_9 = z_9 * z_9;
    float xy_12 = x_24 * y_10;
    float xz_9 = x_24 * z_9;
    float yz_9 = y_10 * z_9;
    float wx_9 = w_9 * x_24;
    float wy_9 = w_9 * y_10;
    float wz_9 = w_9 * z_9;
    Matrix<float, 3, 3>  _S107 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_9), 2.0f * (xy_12 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_12 - wz_9), 1.0f - 2.0f * (x2_12 + z2_9), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_12 + y2_12)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S103.x, 0.0f, 0.0f, 0.0f, _S103.y, 0.0f, 0.0f, 0.0f, _S103.z);
    Matrix<float, 3, 3>  _S108 = s_primal_ctx_mul_1(_S107, S_0);
    Matrix<float, 3, 3>  _S109 = transpose_0(_S108);
    Matrix<float, 3, 3>  _S110 = s_primal_ctx_mul_1(_S108, _S109);
    Matrix<float, 3, 3>  _S111 = s_primal_ctx_mul_1(R_8, _S110);
    Matrix<float, 3, 3>  _S112 = transpose_0(R_8);
    Matrix<float, 3, 3>  _S113 = s_primal_ctx_mul_1(_S111, _S112);
    float _S114 = float(image_width_6);
    float _S115 = float(image_height_6);
    float _S116 = 0.30000001192092896f * (0.5f * _S114 / fx_10);
    float lim_x_pos_0 = (_S114 - cx_10) / fx_10 + _S116;
    float _S117 = 0.30000001192092896f * (0.5f * _S115 / fy_10);
    float lim_y_pos_0 = (_S115 - cy_10) / fy_10 + _S117;
    float rz_4 = 1.0f / mean_c_6.z;
    float _S118 = mean_c_6.z * mean_c_6.z;
    float rz2_4 = rz_4 * rz_4;
    float _S119 = - (cx_10 / fx_10 + _S116);
    float _S120 = mean_c_6.x * rz_4;
    float _S121 = s_primal_ctx_max_0(_S119, _S120);
    float _S122 = s_primal_ctx_min_0(lim_x_pos_0, _S121);
    float _S123 = - (cy_10 / fy_10 + _S117);
    float _S124 = mean_c_6.y * rz_4;
    float _S125 = s_primal_ctx_max_0(_S123, _S124);
    float _S126 = s_primal_ctx_min_0(lim_y_pos_0, _S125);
    float _S127 = - fx_10;
    float _S128 = _S127 * (mean_c_6.z * _S122);
    float _S129 = - fy_10;
    float _S130 = _S129 * (mean_c_6.z * _S126);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (fx_10 * rz_4, 0.0f, _S128 * rz2_4, 0.0f, fy_10 * rz_4, _S130 * rz2_4);
    Matrix<float, 2, 3>  _S131 = s_primal_ctx_mul_2(J_10, _S113);
    Matrix<float, 3, 2>  _S132 = transpose_1(J_10);
    Matrix<float, 2, 2>  _S133 = s_primal_ctx_mul_3(_S131, _S132);
    float _S134 = fx_10 * mean_c_6.x;
    float _S135 = fy_10 * mean_c_6.y;
    float _S136 = _S133.rows[int(0)].y * _S133.rows[int(1)].x;
    float det_orig_7 = _S133.rows[int(0)].x * _S133.rows[int(1)].y - _S136;
    float _S137 = _S133.rows[int(0)].x + eps2d_7;
    Matrix<float, 2, 2>  _S138 = _S133;
    *&(((&_S138)->rows + (int(0)))->x) = _S137;
    float _S139 = _S133.rows[int(1)].y + eps2d_7;
    *&(((&_S138)->rows + (int(1)))->y) = _S139;
    Matrix<float, 2, 2>  _S140 = _S138;
    Matrix<float, 2, 2>  _S141 = _S138;
    float det_blur_4 = _S137 * _S139 - _S136;
    float _S142 = det_orig_7 / det_blur_4;
    float _S143 = det_blur_4 * det_blur_4;
    float _S144 = s_primal_ctx_max_0(0.0f, _S142);
    float _S145 = s_primal_ctx_sqrt_0(_S144);
    float invdet_7 = 1.0f / det_blur_4;
    float _S146 = - _S133.rows[int(0)].y;
    float _S147 = - _S133.rows[int(1)].x;
    float _S148 = - in_opacity_6;
    float _S149 = 1.0f + s_primal_ctx_exp_1(_S148);
    float _S150 = 1.0f / _S149;
    float _S151 = _S149 * _S149;
    float _S152;
    if(antialiased_6)
    {
        _S152 = _S150 * _S145;
    }
    else
    {
        _S152 = _S150;
    }
    float _S153 = _S152 / 0.00392156885936856f;
    float _S154 = 2.0f * s_primal_ctx_log_0(_S153);
    float _S155 = s_primal_ctx_sqrt_0(_S154);
    float _S156 = _S140.rows[int(0)].x;
    float _S157 = _S141.rows[int(1)].y;
    float2  _S158 = make_float2 (0.0f);
    float2  _S159 = _S158;
    *&((&_S159)->y) = v_conic_0.z;
    float2  _S160 = _S158;
    *&((&_S160)->y) = v_conic_0.y;
    *&((&_S160)->x) = v_conic_0.x;
    DiffPair_float_0 _S161;
    (&_S161)->primal_0 = _S157;
    (&_S161)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S161, 0.0f);
    DiffPair_float_0 _S162;
    (&_S162)->primal_0 = _S156;
    (&_S162)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S162, 0.0f);
    DiffPair_float_0 _S163;
    (&_S163)->primal_0 = 3.32999992370605469f;
    (&_S163)->differential_0 = 0.0f;
    DiffPair_float_0 _S164;
    (&_S164)->primal_0 = _S155;
    (&_S164)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S163, &_S164, 0.0f);
    DiffPair_float_0 _S165;
    (&_S165)->primal_0 = _S154;
    (&_S165)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S165, _S164.differential_0);
    float _S166 = 2.0f * _S165.differential_0;
    DiffPair_float_0 _S167;
    (&_S167)->primal_0 = _S153;
    (&_S167)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S167, _S166);
    float _S168 = v_opacity_0 + 254.9999847412109375f * _S167.differential_0;
    Matrix<float, 2, 2>  _S169 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S170 = _S169;
    _S170[int(1)] = _S159;
    _S170[int(0)] = _S160;
    Matrix<float, 2, 2>  _S171 = _S170;
    float3  _S172 = make_float3 (0.0f, 0.0f, v_depth_0);
    float2  _S173 = make_float2 (_S162.differential_0, 0.0f);
    float2  _S174 = make_float2 (0.0f, _S161.differential_0);
    float _S175;
    if(antialiased_6)
    {
        float _S176 = _S150 * _S168;
        _S152 = _S145 * _S168;
        _S175 = _S176;
    }
    else
    {
        _S152 = _S168;
        _S175 = 0.0f;
    }
    float _S177 = - (_S152 / _S151);
    DiffPair_float_0 _S178;
    (&_S178)->primal_0 = _S148;
    (&_S178)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S178, _S177);
    float _S179 = - _S178.differential_0;
    float _S180 = invdet_7 * _S171.rows[int(1)].y;
    float _S181 = - (invdet_7 * _S171.rows[int(1)].x);
    float _S182 = - (invdet_7 * _S171.rows[int(0)].y);
    float _S183 = invdet_7 * _S171.rows[int(0)].x;
    float _S184 = - ((_S137 * _S171.rows[int(1)].y + _S147 * _S171.rows[int(1)].x + _S146 * _S171.rows[int(0)].y + _S139 * _S171.rows[int(0)].x) / _S143);
    DiffPair_float_0 _S185;
    (&_S185)->primal_0 = _S144;
    (&_S185)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S185, _S175);
    DiffPair_float_0 _S186;
    (&_S186)->primal_0 = 0.0f;
    (&_S186)->differential_0 = 0.0f;
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = _S142;
    (&_S187)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S186, &_S187, _S185.differential_0);
    float _S188 = _S187.differential_0 / _S143;
    float s_diff_det_orig_T_0 = det_blur_4 * _S188;
    float _S189 = _S184 + det_orig_7 * - _S188;
    float _S190 = - _S189;
    float _S191 = _S137 * _S189;
    float _S192 = _S139 * _S189;
    Matrix<float, 2, 2>  _S193 = _S169;
    _S193[int(1)] = _S174;
    _S193[int(0)] = _S173;
    _S138 = _S193;
    *&(((&_S138)->rows + (int(1)))->y) = 0.0f;
    float _S194 = _S183 + _S191 + _S193.rows[int(1)].y;
    *&(((&_S138)->rows + (int(0)))->x) = 0.0f;
    float _S195 = _S180 + _S192 + _S193.rows[int(0)].x;
    float _S196 = _S190 + - s_diff_det_orig_T_0;
    float _S197 = _S181 + _S133.rows[int(0)].y * _S196;
    float _S198 = _S182 + _S133.rows[int(1)].x * _S196;
    float _S199 = _S133.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S200 = _S194 + _S133.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S201 = _S158;
    *&((&_S201)->x) = _S197;
    *&((&_S201)->y) = _S200;
    float _S202 = _S195 + _S199;
    float2  _S203 = _S158;
    *&((&_S203)->y) = _S198;
    *&((&_S203)->x) = _S202;
    float _S204 = _S135 * v_mean2d_0.y;
    float _S205 = fy_10 * (rz_4 * v_mean2d_0.y);
    float _S206 = _S134 * v_mean2d_0.x;
    float _S207 = fx_10 * (rz_4 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S208 = _S169;
    _S208[int(1)] = _S201;
    _S208[int(0)] = _S203;
    Matrix<float, 2, 2>  _S209 = _S138 + _S208;
    Matrix<float, 2, 3>  _S210 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S211;
    (&_S211)->primal_0 = _S131;
    (&_S211)->differential_0 = _S210;
    Matrix<float, 3, 2>  _S212 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S213;
    (&_S213)->primal_0 = _S132;
    (&_S213)->differential_0 = _S212;
    s_bwd_prop_mul_0(&_S211, &_S213, _S209);
    Matrix<float, 2, 3>  _S214 = transpose_2(_S213.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S215;
    (&_S215)->primal_0 = J_10;
    (&_S215)->differential_0 = _S210;
    Matrix<float, 3, 3>  _S216 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S217;
    (&_S217)->primal_0 = _S113;
    (&_S217)->differential_0 = _S216;
    s_bwd_prop_mul_1(&_S215, &_S217, _S211.differential_0);
    Matrix<float, 2, 3>  _S218 = _S214 + _S215.differential_0;
    float _S219 = _S130 * _S218.rows[int(1)].z;
    float s_diff_ty_T_0 = _S129 * (rz2_4 * _S218.rows[int(1)].z);
    float _S220 = fy_10 * _S218.rows[int(1)].y;
    float _S221 = _S128 * _S218.rows[int(0)].z;
    float s_diff_tx_T_0 = _S127 * (rz2_4 * _S218.rows[int(0)].z);
    float _S222 = fx_10 * _S218.rows[int(0)].x;
    float _S223 = mean_c_6.z * s_diff_ty_T_0;
    float _S224 = _S126 * s_diff_ty_T_0;
    DiffPair_float_0 _S225;
    (&_S225)->primal_0 = lim_y_pos_0;
    (&_S225)->differential_0 = 0.0f;
    DiffPair_float_0 _S226;
    (&_S226)->primal_0 = _S125;
    (&_S226)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S225, &_S226, _S223);
    DiffPair_float_0 _S227;
    (&_S227)->primal_0 = _S123;
    (&_S227)->differential_0 = 0.0f;
    DiffPair_float_0 _S228;
    (&_S228)->primal_0 = _S124;
    (&_S228)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S227, &_S228, _S226.differential_0);
    float _S229 = mean_c_6.y * _S228.differential_0;
    float _S230 = rz_4 * _S228.differential_0;
    float _S231 = mean_c_6.z * s_diff_tx_T_0;
    float _S232 = _S122 * s_diff_tx_T_0;
    DiffPair_float_0 _S233;
    (&_S233)->primal_0 = lim_x_pos_0;
    (&_S233)->differential_0 = 0.0f;
    DiffPair_float_0 _S234;
    (&_S234)->primal_0 = _S121;
    (&_S234)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S233, &_S234, _S231);
    DiffPair_float_0 _S235;
    (&_S235)->primal_0 = _S119;
    (&_S235)->differential_0 = 0.0f;
    DiffPair_float_0 _S236;
    (&_S236)->primal_0 = _S120;
    (&_S236)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S235, &_S236, _S234.differential_0);
    float _S237 = rz_4 * (_S219 + _S221);
    float _S238 = _S224 + _S232 + - ((_S204 + _S206 + _S220 + _S222 + _S229 + mean_c_6.x * _S236.differential_0 + _S237 + _S237) / _S118);
    float _S239 = _S205 + _S230;
    float _S240 = _S207 + rz_4 * _S236.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S241;
    (&_S241)->primal_0 = _S111;
    (&_S241)->differential_0 = _S216;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S242;
    (&_S242)->primal_0 = _S112;
    (&_S242)->differential_0 = _S216;
    s_bwd_prop_mul_2(&_S241, &_S242, _S217.differential_0);
    Matrix<float, 3, 3>  _S243 = transpose_0(_S242.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S244;
    (&_S244)->primal_0 = R_8;
    (&_S244)->differential_0 = _S216;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S245;
    (&_S245)->primal_0 = _S110;
    (&_S245)->differential_0 = _S216;
    s_bwd_prop_mul_2(&_S244, &_S245, _S241.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S246;
    (&_S246)->primal_0 = _S108;
    (&_S246)->differential_0 = _S216;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S247;
    (&_S247)->primal_0 = _S109;
    (&_S247)->differential_0 = _S216;
    s_bwd_prop_mul_2(&_S246, &_S247, _S245.differential_0);
    Matrix<float, 3, 3>  _S248 = _S246.differential_0 + transpose_0(_S247.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S249;
    (&_S249)->primal_0 = _S107;
    (&_S249)->differential_0 = _S216;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S250;
    (&_S250)->primal_0 = S_0;
    (&_S250)->differential_0 = _S216;
    s_bwd_prop_mul_2(&_S249, &_S250, _S248);
    Matrix<float, 3, 3>  _S251 = transpose_0(_S249.differential_0);
    float _S252 = 2.0f * - _S251.rows[int(2)].z;
    float _S253 = 2.0f * _S251.rows[int(2)].y;
    float _S254 = 2.0f * _S251.rows[int(2)].x;
    float _S255 = 2.0f * _S251.rows[int(1)].z;
    float _S256 = 2.0f * - _S251.rows[int(1)].y;
    float _S257 = 2.0f * _S251.rows[int(1)].x;
    float _S258 = 2.0f * _S251.rows[int(0)].z;
    float _S259 = 2.0f * _S251.rows[int(0)].y;
    float _S260 = 2.0f * - _S251.rows[int(0)].x;
    float _S261 = - _S257 + _S259;
    float _S262 = _S254 + - _S258;
    float _S263 = - _S253 + _S255;
    float _S264 = _S253 + _S255;
    float _S265 = _S254 + _S258;
    float _S266 = _S257 + _S259;
    float _S267 = z_9 * (_S256 + _S260);
    float _S268 = y_10 * (_S252 + _S260);
    float _S269 = x_24 * (_S252 + _S256);
    float _S270 = z_9 * _S261 + y_10 * _S262 + x_24 * _S263;
    float _S271 = _S106 * _S270;
    float _S272 = w_9 * _S261 + y_10 * _S264 + x_24 * _S265 + _S267 + _S267;
    float _S273 = _S106 * _S272;
    float _S274 = w_9 * _S262 + z_9 * _S264 + x_24 * _S266 + _S268 + _S268;
    float _S275 = _S106 * _S274;
    float _S276 = w_9 * _S263 + z_9 * _S265 + y_10 * _S266 + _S269 + _S269;
    float _S277 = _S106 * _S276;
    float _S278 = quat_9.x * _S270 + quat_9.w * _S272 + quat_9.z * _S274 + quat_9.y * _S276;
    DiffPair_float_0 _S279;
    (&_S279)->primal_0 = _S105;
    (&_S279)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S279, _S278);
    float _S280 = quat_9.x * _S279.differential_0;
    float _S281 = quat_9.w * _S279.differential_0;
    float _S282 = quat_9.z * _S279.differential_0;
    float _S283 = quat_9.y * _S279.differential_0;
    float _S284 = _S273 + _S281 + _S281;
    float _S285 = _S275 + _S282 + _S282;
    float _S286 = _S277 + _S283 + _S283;
    float _S287 = _S271 + _S280 + _S280;
    float3  _S288 = make_float3 (0.0f);
    float3  _S289 = _S288;
    *&((&_S289)->z) = _S250.differential_0.rows[int(2)].z;
    *&((&_S289)->y) = _S250.differential_0.rows[int(1)].y;
    *&((&_S289)->x) = _S250.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S290;
    (&_S290)->primal_0 = scale_8;
    (&_S290)->differential_0 = _S288;
    s_bwd_prop_exp_1(&_S290, _S289);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S291 = _S290;
    float3  _S292 = _S288;
    *&((&_S292)->z) = _S238;
    *&((&_S292)->y) = _S239;
    *&((&_S292)->x) = _S240;
    float3  _S293 = _S172 + _S292;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S294;
    (&_S294)->primal_0 = R_8;
    (&_S294)->differential_0 = _S216;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S295;
    (&_S295)->primal_0 = mean_6;
    (&_S295)->differential_0 = _S288;
    s_bwd_prop_mul_3(&_S294, &_S295, _S293);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S296 = _S295;
    Matrix<float, 3, 3>  _S297 = _S243 + _S244.differential_0 + _S294.differential_0;
    float4  _S298 = make_float4 (0.0f);
    *&((&_S298)->w) = _S284;
    *&((&_S298)->z) = _S285;
    *&((&_S298)->y) = _S286;
    *&((&_S298)->x) = _S287;
    float4  _S299 = _S298;
    *v_mean_0 = _S296.differential_0;
    *v_quat_0 = _S299;
    *v_scale_0 = _S291.differential_0;
    *v_in_opacity_0 = _S179;
    *v_R_0 = _S297;
    *v_t_0 = _S293;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S300, float _S301)
{
    return (F32_atan2((_S300), (_S301)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S302, DiffPair_float_0 * _S303, float _S304)
{
    _d_atan2_0(_S302, _S303, _S304);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_9, float _s_dOut_0)
{
    float _S305 = (*dpx_9).primal_0.x;
    float _S306 = (*dpx_9).primal_0.y;
    DiffPair_float_0 _S307;
    (&_S307)->primal_0 = _S305 * _S305 + _S306 * _S306;
    (&_S307)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S307, _s_dOut_0);
    float _S308 = (*dpx_9).primal_0.y * _S307.differential_0;
    float _S309 = _S308 + _S308;
    float _S310 = (*dpx_9).primal_0.x * _S307.differential_0;
    float _S311 = _S310 + _S310;
    float2  _S312 = make_float2 (0.0f);
    *&((&_S312)->y) = _S309;
    *&((&_S312)->x) = _S311;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S312;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S313, float _S314)
{
    s_bwd_prop_length_impl_0(_S313, _S314);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_11, float fy_11, float cx_11, float cy_11, uint image_width_7, uint image_height_7, float eps2d_8, float v_depth_1, float2  v_mean2d_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_7 = s_primal_ctx_mul_0(R_9, mean_7) + t_8;
    float3  _S315 = s_primal_ctx_exp_0(scale_9);
    float _S316 = quat_10.y;
    float _S317 = _S316 * _S316 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S318 = s_primal_ctx_rsqrt_0(_S317);
    float x_25 = quat_10.y * _S318;
    float y_11 = quat_10.z * _S318;
    float z_10 = quat_10.w * _S318;
    float w_10 = quat_10.x * _S318;
    float x2_13 = x_25 * x_25;
    float y2_13 = y_11 * y_11;
    float z2_10 = z_10 * z_10;
    float xy_13 = x_25 * y_11;
    float xz_10 = x_25 * z_10;
    float yz_10 = y_11 * z_10;
    float wx_10 = w_10 * x_25;
    float wy_10 = w_10 * y_11;
    float wz_10 = w_10 * z_10;
    Matrix<float, 3, 3>  _S319 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_10), 2.0f * (xy_13 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_13 - wz_10), 1.0f - 2.0f * (x2_13 + z2_10), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S315.x, 0.0f, 0.0f, 0.0f, _S315.y, 0.0f, 0.0f, 0.0f, _S315.z);
    Matrix<float, 3, 3>  _S320 = s_primal_ctx_mul_1(_S319, S_1);
    Matrix<float, 3, 3>  _S321 = transpose_0(_S320);
    Matrix<float, 3, 3>  _S322 = s_primal_ctx_mul_1(_S320, _S321);
    Matrix<float, 3, 3>  _S323 = s_primal_ctx_mul_1(R_9, _S322);
    Matrix<float, 3, 3>  _S324 = transpose_0(R_9);
    Matrix<float, 3, 3>  _S325 = s_primal_ctx_mul_1(_S323, _S324);
    float2  _S326 = make_float2 (mean_c_7.x, mean_c_7.y);
    float xy_len_3 = length_0(_S326) + 1.00000001168609742e-07f;
    float _S327 = mean_c_7.z + 1.00000001168609742e-07f;
    float _S328 = s_primal_ctx_atan2_0(xy_len_3, _S327);
    float _S329 = mean_c_7.x * fx_11;
    float _S330 = _S329 * _S328;
    float _S331 = xy_len_3 * xy_len_3;
    float _S332 = mean_c_7.y * fy_11;
    float _S333 = _S332 * _S328;
    float x2_14 = mean_c_7.x * mean_c_7.x + 1.00000001168609742e-07f;
    float y2_14 = mean_c_7.y * mean_c_7.y;
    float xy_14 = mean_c_7.x * mean_c_7.y;
    float x2y2_3 = x2_14 + y2_14;
    float _S334 = x2y2_3 + mean_c_7.z * mean_c_7.z;
    float x2y2z2_inv_3 = 1.0f / _S334;
    float _S335 = _S334 * _S334;
    float _S336 = s_primal_ctx_atan2_0(xy_len_3, mean_c_7.z);
    float _S337 = _S336 / xy_len_3;
    float b_3 = _S337 / x2y2_3;
    float _S338 = x2y2_3 * x2y2_3;
    float _S339 = mean_c_7.z * x2y2z2_inv_3;
    float a_3 = _S339 / x2y2_3;
    float _S340 = fx_11 * xy_14;
    float _S341 = a_3 - b_3;
    float _S342 = - fx_11;
    float _S343 = _S342 * mean_c_7.x;
    float _S344 = fy_11 * xy_14;
    float _S345 = - fy_11;
    float _S346 = _S345 * mean_c_7.y;
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_11 * (x2_14 * a_3 + y2_14 * b_3), _S340 * _S341, _S343 * x2y2z2_inv_3, _S344 * _S341, fy_11 * (y2_14 * a_3 + x2_14 * b_3), _S346 * x2y2z2_inv_3);
    Matrix<float, 2, 3>  _S347 = s_primal_ctx_mul_2(J_11, _S325);
    Matrix<float, 3, 2>  _S348 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S349 = s_primal_ctx_mul_3(_S347, _S348);
    float _S350 = _S349.rows[int(0)].y * _S349.rows[int(1)].x;
    float det_orig_8 = _S349.rows[int(0)].x * _S349.rows[int(1)].y - _S350;
    float _S351 = _S349.rows[int(0)].x + eps2d_8;
    Matrix<float, 2, 2>  _S352 = _S349;
    *&(((&_S352)->rows + (int(0)))->x) = _S351;
    float _S353 = _S349.rows[int(1)].y + eps2d_8;
    *&(((&_S352)->rows + (int(1)))->y) = _S353;
    Matrix<float, 2, 2>  _S354 = _S352;
    Matrix<float, 2, 2>  _S355 = _S352;
    float det_blur_5 = _S351 * _S353 - _S350;
    float _S356 = det_orig_8 / det_blur_5;
    float _S357 = det_blur_5 * det_blur_5;
    float _S358 = s_primal_ctx_max_0(0.0f, _S356);
    float _S359 = s_primal_ctx_sqrt_0(_S358);
    float invdet_8 = 1.0f / det_blur_5;
    float _S360 = - _S349.rows[int(0)].y;
    float _S361 = - _S349.rows[int(1)].x;
    float _S362 = - in_opacity_7;
    float _S363 = 1.0f + s_primal_ctx_exp_1(_S362);
    float _S364 = 1.0f / _S363;
    float _S365 = _S363 * _S363;
    float _S366;
    if(antialiased_7)
    {
        _S366 = _S364 * _S359;
    }
    else
    {
        _S366 = _S364;
    }
    float _S367 = _S366 / 0.00392156885936856f;
    float _S368 = 2.0f * s_primal_ctx_log_0(_S367);
    float _S369 = s_primal_ctx_sqrt_0(_S368);
    float _S370 = _S354.rows[int(0)].x;
    float _S371 = _S355.rows[int(1)].y;
    float2  _S372 = make_float2 (0.0f);
    float2  _S373 = _S372;
    *&((&_S373)->y) = v_conic_1.z;
    float2  _S374 = _S372;
    *&((&_S374)->y) = v_conic_1.y;
    *&((&_S374)->x) = v_conic_1.x;
    DiffPair_float_0 _S375;
    (&_S375)->primal_0 = _S371;
    (&_S375)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S375, 0.0f);
    DiffPair_float_0 _S376;
    (&_S376)->primal_0 = _S370;
    (&_S376)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S376, 0.0f);
    DiffPair_float_0 _S377;
    (&_S377)->primal_0 = 3.32999992370605469f;
    (&_S377)->differential_0 = 0.0f;
    DiffPair_float_0 _S378;
    (&_S378)->primal_0 = _S369;
    (&_S378)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S377, &_S378, 0.0f);
    DiffPair_float_0 _S379;
    (&_S379)->primal_0 = _S368;
    (&_S379)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S379, _S378.differential_0);
    float _S380 = 2.0f * _S379.differential_0;
    DiffPair_float_0 _S381;
    (&_S381)->primal_0 = _S367;
    (&_S381)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S381, _S380);
    float _S382 = v_opacity_1 + 254.9999847412109375f * _S381.differential_0;
    Matrix<float, 2, 2>  _S383 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S384 = _S383;
    _S384[int(1)] = _S373;
    _S384[int(0)] = _S374;
    Matrix<float, 2, 2>  _S385 = _S384;
    float3  _S386 = make_float3 (0.0f, 0.0f, v_depth_1);
    float2  _S387 = make_float2 (_S376.differential_0, 0.0f);
    float2  _S388 = make_float2 (0.0f, _S375.differential_0);
    float _S389;
    if(antialiased_7)
    {
        float _S390 = _S364 * _S382;
        _S366 = _S359 * _S382;
        _S389 = _S390;
    }
    else
    {
        _S366 = _S382;
        _S389 = 0.0f;
    }
    float _S391 = - (_S366 / _S365);
    DiffPair_float_0 _S392;
    (&_S392)->primal_0 = _S362;
    (&_S392)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S392, _S391);
    float _S393 = - _S392.differential_0;
    float _S394 = invdet_8 * _S385.rows[int(1)].y;
    float _S395 = - (invdet_8 * _S385.rows[int(1)].x);
    float _S396 = - (invdet_8 * _S385.rows[int(0)].y);
    float _S397 = invdet_8 * _S385.rows[int(0)].x;
    float _S398 = - ((_S351 * _S385.rows[int(1)].y + _S361 * _S385.rows[int(1)].x + _S360 * _S385.rows[int(0)].y + _S353 * _S385.rows[int(0)].x) / _S357);
    DiffPair_float_0 _S399;
    (&_S399)->primal_0 = _S358;
    (&_S399)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S399, _S389);
    DiffPair_float_0 _S400;
    (&_S400)->primal_0 = 0.0f;
    (&_S400)->differential_0 = 0.0f;
    DiffPair_float_0 _S401;
    (&_S401)->primal_0 = _S356;
    (&_S401)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S400, &_S401, _S399.differential_0);
    float _S402 = _S401.differential_0 / _S357;
    float s_diff_det_orig_T_1 = det_blur_5 * _S402;
    float _S403 = _S398 + det_orig_8 * - _S402;
    float _S404 = - _S403;
    float _S405 = _S351 * _S403;
    float _S406 = _S353 * _S403;
    Matrix<float, 2, 2>  _S407 = _S383;
    _S407[int(1)] = _S388;
    _S407[int(0)] = _S387;
    _S352 = _S407;
    *&(((&_S352)->rows + (int(1)))->y) = 0.0f;
    float _S408 = _S397 + _S405 + _S407.rows[int(1)].y;
    *&(((&_S352)->rows + (int(0)))->x) = 0.0f;
    float _S409 = _S394 + _S406 + _S407.rows[int(0)].x;
    float _S410 = _S404 + - s_diff_det_orig_T_1;
    float _S411 = _S395 + _S349.rows[int(0)].y * _S410;
    float _S412 = _S396 + _S349.rows[int(1)].x * _S410;
    float _S413 = _S349.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S414 = _S408 + _S349.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S415 = _S372;
    *&((&_S415)->x) = _S411;
    *&((&_S415)->y) = _S414;
    float _S416 = _S409 + _S413;
    float2  _S417 = _S372;
    *&((&_S417)->y) = _S412;
    *&((&_S417)->x) = _S416;
    Matrix<float, 2, 2>  _S418 = _S383;
    _S418[int(1)] = _S415;
    _S418[int(0)] = _S417;
    Matrix<float, 2, 2>  _S419 = _S352 + _S418;
    Matrix<float, 2, 3>  _S420 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S421;
    (&_S421)->primal_0 = _S347;
    (&_S421)->differential_0 = _S420;
    Matrix<float, 3, 2>  _S422 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S423;
    (&_S423)->primal_0 = _S348;
    (&_S423)->differential_0 = _S422;
    s_bwd_prop_mul_0(&_S421, &_S423, _S419);
    Matrix<float, 2, 3>  _S424 = transpose_2(_S423.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S425;
    (&_S425)->primal_0 = J_11;
    (&_S425)->differential_0 = _S420;
    Matrix<float, 3, 3>  _S426 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S427;
    (&_S427)->primal_0 = _S325;
    (&_S427)->differential_0 = _S426;
    s_bwd_prop_mul_1(&_S425, &_S427, _S421.differential_0);
    Matrix<float, 2, 3>  _S428 = _S424 + _S425.differential_0;
    float _S429 = _S346 * _S428.rows[int(1)].z;
    float _S430 = _S345 * (x2y2z2_inv_3 * _S428.rows[int(1)].z);
    float _S431 = fy_11 * _S428.rows[int(1)].y;
    float _S432 = b_3 * _S431;
    float _S433 = a_3 * _S431;
    float _S434 = fy_11 * (_S341 * _S428.rows[int(1)].x);
    float _S435 = _S343 * _S428.rows[int(0)].z;
    float _S436 = _S342 * (x2y2z2_inv_3 * _S428.rows[int(0)].z);
    float _S437 = _S344 * _S428.rows[int(1)].x + _S340 * _S428.rows[int(0)].y;
    float _S438 = fx_11 * (_S341 * _S428.rows[int(0)].y);
    float _S439 = fx_11 * _S428.rows[int(0)].x;
    float _S440 = b_3 * _S439;
    float _S441 = a_3 * _S439;
    float _S442 = (y2_14 * _S431 + _S437 + x2_14 * _S439) / _S338;
    float _S443 = _S339 * - _S442;
    float _S444 = x2y2_3 * _S442;
    float _S445 = mean_c_7.z * _S444;
    float _S446 = x2y2z2_inv_3 * _S444;
    float _S447 = (x2_14 * _S431 + - _S437 + y2_14 * _S439) / _S338;
    float _S448 = _S337 * - _S447;
    float _S449 = x2y2_3 * _S447 / _S331;
    float _S450 = _S336 * - _S449;
    float _S451 = xy_len_3 * _S449;
    DiffPair_float_0 _S452;
    (&_S452)->primal_0 = xy_len_3;
    (&_S452)->differential_0 = 0.0f;
    DiffPair_float_0 _S453;
    (&_S453)->primal_0 = mean_c_7.z;
    (&_S453)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S452, &_S453, _S451);
    float _S454 = - ((_S429 + _S435 + _S445) / _S335);
    float _S455 = mean_c_7.z * _S454;
    float _S456 = _S443 + _S448 + _S454;
    float _S457 = _S434 + _S438;
    float _S458 = mean_c_7.x * _S457;
    float _S459 = mean_c_7.y * _S457;
    float _S460 = mean_c_7.y * (_S433 + _S440 + _S456);
    float _S461 = mean_c_7.x * (_S432 + _S441 + _S456);
    float _S462 = v_mean2d_1.y / _S331;
    float _S463 = _S333 * - _S462;
    float _S464 = xy_len_3 * _S462;
    float _S465 = fy_11 * (_S328 * _S464);
    float _S466 = v_mean2d_1.x / _S331;
    float _S467 = _S330 * - _S466;
    float _S468 = xy_len_3 * _S466;
    float _S469 = fx_11 * (_S328 * _S468);
    float _S470 = _S332 * _S464 + _S329 * _S468;
    DiffPair_float_0 _S471;
    (&_S471)->primal_0 = xy_len_3;
    (&_S471)->differential_0 = 0.0f;
    DiffPair_float_0 _S472;
    (&_S472)->primal_0 = _S327;
    (&_S472)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S471, &_S472, _S470);
    float _S473 = _S450 + _S452.differential_0 + _S463 + _S467 + _S471.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S474;
    (&_S474)->primal_0 = _S326;
    (&_S474)->differential_0 = _S372;
    s_bwd_length_impl_0(&_S474, _S473);
    float _S475 = _S446 + _S453.differential_0 + _S455 + _S455 + _S472.differential_0;
    float _S476 = _S430 + _S458 + _S460 + _S460 + _S465 + _S474.differential_0.y;
    float _S477 = _S436 + _S459 + _S461 + _S461 + _S469 + _S474.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S478;
    (&_S478)->primal_0 = _S323;
    (&_S478)->differential_0 = _S426;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S479;
    (&_S479)->primal_0 = _S324;
    (&_S479)->differential_0 = _S426;
    s_bwd_prop_mul_2(&_S478, &_S479, _S427.differential_0);
    Matrix<float, 3, 3>  _S480 = transpose_0(_S479.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S481;
    (&_S481)->primal_0 = R_9;
    (&_S481)->differential_0 = _S426;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S482;
    (&_S482)->primal_0 = _S322;
    (&_S482)->differential_0 = _S426;
    s_bwd_prop_mul_2(&_S481, &_S482, _S478.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S483;
    (&_S483)->primal_0 = _S320;
    (&_S483)->differential_0 = _S426;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S484;
    (&_S484)->primal_0 = _S321;
    (&_S484)->differential_0 = _S426;
    s_bwd_prop_mul_2(&_S483, &_S484, _S482.differential_0);
    Matrix<float, 3, 3>  _S485 = _S483.differential_0 + transpose_0(_S484.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S486;
    (&_S486)->primal_0 = _S319;
    (&_S486)->differential_0 = _S426;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S487;
    (&_S487)->primal_0 = S_1;
    (&_S487)->differential_0 = _S426;
    s_bwd_prop_mul_2(&_S486, &_S487, _S485);
    Matrix<float, 3, 3>  _S488 = transpose_0(_S486.differential_0);
    float _S489 = 2.0f * - _S488.rows[int(2)].z;
    float _S490 = 2.0f * _S488.rows[int(2)].y;
    float _S491 = 2.0f * _S488.rows[int(2)].x;
    float _S492 = 2.0f * _S488.rows[int(1)].z;
    float _S493 = 2.0f * - _S488.rows[int(1)].y;
    float _S494 = 2.0f * _S488.rows[int(1)].x;
    float _S495 = 2.0f * _S488.rows[int(0)].z;
    float _S496 = 2.0f * _S488.rows[int(0)].y;
    float _S497 = 2.0f * - _S488.rows[int(0)].x;
    float _S498 = - _S494 + _S496;
    float _S499 = _S491 + - _S495;
    float _S500 = - _S490 + _S492;
    float _S501 = _S490 + _S492;
    float _S502 = _S491 + _S495;
    float _S503 = _S494 + _S496;
    float _S504 = z_10 * (_S493 + _S497);
    float _S505 = y_11 * (_S489 + _S497);
    float _S506 = x_25 * (_S489 + _S493);
    float _S507 = z_10 * _S498 + y_11 * _S499 + x_25 * _S500;
    float _S508 = _S318 * _S507;
    float _S509 = w_10 * _S498 + y_11 * _S501 + x_25 * _S502 + _S504 + _S504;
    float _S510 = _S318 * _S509;
    float _S511 = w_10 * _S499 + z_10 * _S501 + x_25 * _S503 + _S505 + _S505;
    float _S512 = _S318 * _S511;
    float _S513 = w_10 * _S500 + z_10 * _S502 + y_11 * _S503 + _S506 + _S506;
    float _S514 = _S318 * _S513;
    float _S515 = quat_10.x * _S507 + quat_10.w * _S509 + quat_10.z * _S511 + quat_10.y * _S513;
    DiffPair_float_0 _S516;
    (&_S516)->primal_0 = _S317;
    (&_S516)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S516, _S515);
    float _S517 = quat_10.x * _S516.differential_0;
    float _S518 = quat_10.w * _S516.differential_0;
    float _S519 = quat_10.z * _S516.differential_0;
    float _S520 = quat_10.y * _S516.differential_0;
    float _S521 = _S510 + _S518 + _S518;
    float _S522 = _S512 + _S519 + _S519;
    float _S523 = _S514 + _S520 + _S520;
    float _S524 = _S508 + _S517 + _S517;
    float3  _S525 = make_float3 (0.0f);
    float3  _S526 = _S525;
    *&((&_S526)->z) = _S487.differential_0.rows[int(2)].z;
    *&((&_S526)->y) = _S487.differential_0.rows[int(1)].y;
    *&((&_S526)->x) = _S487.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S527;
    (&_S527)->primal_0 = scale_9;
    (&_S527)->differential_0 = _S525;
    s_bwd_prop_exp_1(&_S527, _S526);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S528 = _S527;
    float3  _S529 = _S525;
    *&((&_S529)->z) = _S475;
    *&((&_S529)->y) = _S476;
    *&((&_S529)->x) = _S477;
    float3  _S530 = _S386 + _S529;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S531;
    (&_S531)->primal_0 = R_9;
    (&_S531)->differential_0 = _S426;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S532;
    (&_S532)->primal_0 = mean_7;
    (&_S532)->differential_0 = _S525;
    s_bwd_prop_mul_3(&_S531, &_S532, _S530);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S533 = _S532;
    Matrix<float, 3, 3>  _S534 = _S480 + _S481.differential_0 + _S531.differential_0;
    float4  _S535 = make_float4 (0.0f);
    *&((&_S535)->w) = _S521;
    *&((&_S535)->z) = _S522;
    *&((&_S535)->y) = _S523;
    *&((&_S535)->x) = _S524;
    float4  _S536 = _S535;
    *v_mean_1 = _S533.differential_0;
    *v_quat_1 = _S536;
    *v_scale_1 = _S528.differential_0;
    *v_in_opacity_1 = _S393;
    *v_R_1 = _S534;
    *v_t_1 = _S530;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_12, float fy_12, float cx_12, float cy_12, uint image_width_8, uint image_height_8, float eps2d_9, float v_depth_2, float2  v_mean2d_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float3  _S537 = s_primal_ctx_exp_0(scale_10);
    float _S538 = quat_11.y;
    float _S539 = _S538 * _S538 + quat_11.z * quat_11.z + quat_11.w * quat_11.w + quat_11.x * quat_11.x;
    float _S540 = s_primal_ctx_rsqrt_0(_S539);
    float x_26 = quat_11.y * _S540;
    float y_12 = quat_11.z * _S540;
    float z_11 = quat_11.w * _S540;
    float w_11 = quat_11.x * _S540;
    float x2_15 = x_26 * x_26;
    float y2_15 = y_12 * y_12;
    float z2_11 = z_11 * z_11;
    float xy_15 = x_26 * y_12;
    float xz_11 = x_26 * z_11;
    float yz_11 = y_12 * z_11;
    float wx_11 = w_11 * x_26;
    float wy_11 = w_11 * y_12;
    float wz_11 = w_11 * z_11;
    Matrix<float, 3, 3>  _S541 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_11), 2.0f * (xy_15 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_15 - wz_11), 1.0f - 2.0f * (x2_15 + z2_11), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S537.x, 0.0f, 0.0f, 0.0f, _S537.y, 0.0f, 0.0f, 0.0f, _S537.z);
    Matrix<float, 3, 3>  _S542 = s_primal_ctx_mul_1(_S541, S_2);
    Matrix<float, 3, 3>  _S543 = transpose_0(_S542);
    Matrix<float, 3, 3>  _S544 = s_primal_ctx_mul_1(_S542, _S543);
    Matrix<float, 3, 3>  _S545 = s_primal_ctx_mul_1(R_10, _S544);
    Matrix<float, 3, 3>  _S546 = transpose_0(R_10);
    Matrix<float, 3, 3>  _S547 = s_primal_ctx_mul_1(_S545, _S546);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (fx_12, 0.0f, 0.0f, 0.0f, fy_12, 0.0f);
    Matrix<float, 2, 3>  _S548 = s_primal_ctx_mul_2(J_12, _S547);
    Matrix<float, 3, 2>  _S549 = transpose_1(J_12);
    Matrix<float, 2, 2>  _S550 = s_primal_ctx_mul_3(_S548, _S549);
    float _S551 = _S550.rows[int(0)].y * _S550.rows[int(1)].x;
    float det_orig_9 = _S550.rows[int(0)].x * _S550.rows[int(1)].y - _S551;
    float _S552 = _S550.rows[int(0)].x + eps2d_9;
    Matrix<float, 2, 2>  _S553 = _S550;
    *&(((&_S553)->rows + (int(0)))->x) = _S552;
    float _S554 = _S550.rows[int(1)].y + eps2d_9;
    *&(((&_S553)->rows + (int(1)))->y) = _S554;
    Matrix<float, 2, 2>  _S555 = _S553;
    Matrix<float, 2, 2>  _S556 = _S553;
    float det_blur_6 = _S552 * _S554 - _S551;
    float _S557 = det_orig_9 / det_blur_6;
    float _S558 = det_blur_6 * det_blur_6;
    float _S559 = s_primal_ctx_max_0(0.0f, _S557);
    float _S560 = s_primal_ctx_sqrt_0(_S559);
    float invdet_9 = 1.0f / det_blur_6;
    float _S561 = - _S550.rows[int(0)].y;
    float _S562 = - _S550.rows[int(1)].x;
    float _S563 = - in_opacity_8;
    float _S564 = 1.0f + s_primal_ctx_exp_1(_S563);
    float _S565 = 1.0f / _S564;
    float _S566 = _S564 * _S564;
    float _S567;
    if(antialiased_8)
    {
        _S567 = _S565 * _S560;
    }
    else
    {
        _S567 = _S565;
    }
    float _S568 = _S567 / 0.00392156885936856f;
    float _S569 = 2.0f * s_primal_ctx_log_0(_S568);
    float _S570 = s_primal_ctx_sqrt_0(_S569);
    float _S571 = _S555.rows[int(0)].x;
    float _S572 = _S556.rows[int(1)].y;
    float2  _S573 = make_float2 (0.0f);
    float2  _S574 = _S573;
    *&((&_S574)->y) = v_conic_2.z;
    float2  _S575 = _S573;
    *&((&_S575)->y) = v_conic_2.y;
    *&((&_S575)->x) = v_conic_2.x;
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S572;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S576, 0.0f);
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = _S571;
    (&_S577)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S577, 0.0f);
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = 3.32999992370605469f;
    (&_S578)->differential_0 = 0.0f;
    DiffPair_float_0 _S579;
    (&_S579)->primal_0 = _S570;
    (&_S579)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S578, &_S579, 0.0f);
    DiffPair_float_0 _S580;
    (&_S580)->primal_0 = _S569;
    (&_S580)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S580, _S579.differential_0);
    float _S581 = 2.0f * _S580.differential_0;
    DiffPair_float_0 _S582;
    (&_S582)->primal_0 = _S568;
    (&_S582)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S582, _S581);
    float _S583 = v_opacity_2 + 254.9999847412109375f * _S582.differential_0;
    Matrix<float, 2, 2>  _S584 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S585 = _S584;
    _S585[int(1)] = _S574;
    _S585[int(0)] = _S575;
    Matrix<float, 2, 2>  _S586 = _S585;
    float3  _S587 = make_float3 (0.0f, 0.0f, v_depth_2);
    float2  _S588 = make_float2 (_S577.differential_0, 0.0f);
    float2  _S589 = make_float2 (0.0f, _S576.differential_0);
    float _S590;
    if(antialiased_8)
    {
        float _S591 = _S565 * _S583;
        _S567 = _S560 * _S583;
        _S590 = _S591;
    }
    else
    {
        _S567 = _S583;
        _S590 = 0.0f;
    }
    float _S592 = - (_S567 / _S566);
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = _S563;
    (&_S593)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S593, _S592);
    float _S594 = - _S593.differential_0;
    float _S595 = invdet_9 * _S586.rows[int(1)].y;
    float _S596 = - (invdet_9 * _S586.rows[int(1)].x);
    float _S597 = - (invdet_9 * _S586.rows[int(0)].y);
    float _S598 = invdet_9 * _S586.rows[int(0)].x;
    float _S599 = - ((_S552 * _S586.rows[int(1)].y + _S562 * _S586.rows[int(1)].x + _S561 * _S586.rows[int(0)].y + _S554 * _S586.rows[int(0)].x) / _S558);
    DiffPair_float_0 _S600;
    (&_S600)->primal_0 = _S559;
    (&_S600)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S600, _S590);
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = 0.0f;
    (&_S601)->differential_0 = 0.0f;
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = _S557;
    (&_S602)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S601, &_S602, _S600.differential_0);
    float _S603 = _S602.differential_0 / _S558;
    float s_diff_det_orig_T_2 = det_blur_6 * _S603;
    float _S604 = _S599 + det_orig_9 * - _S603;
    float _S605 = - _S604;
    float _S606 = _S552 * _S604;
    float _S607 = _S554 * _S604;
    Matrix<float, 2, 2>  _S608 = _S584;
    _S608[int(1)] = _S589;
    _S608[int(0)] = _S588;
    _S553 = _S608;
    *&(((&_S553)->rows + (int(1)))->y) = 0.0f;
    float _S609 = _S598 + _S606 + _S608.rows[int(1)].y;
    *&(((&_S553)->rows + (int(0)))->x) = 0.0f;
    float _S610 = _S595 + _S607 + _S608.rows[int(0)].x;
    float _S611 = _S605 + - s_diff_det_orig_T_2;
    float _S612 = _S596 + _S550.rows[int(0)].y * _S611;
    float _S613 = _S597 + _S550.rows[int(1)].x * _S611;
    float _S614 = _S550.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S615 = _S609 + _S550.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S616 = _S573;
    *&((&_S616)->x) = _S612;
    *&((&_S616)->y) = _S615;
    float _S617 = _S610 + _S614;
    float2  _S618 = _S573;
    *&((&_S618)->y) = _S613;
    *&((&_S618)->x) = _S617;
    float _S619 = fy_12 * v_mean2d_2.y;
    float _S620 = fx_12 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S621 = _S584;
    _S621[int(1)] = _S616;
    _S621[int(0)] = _S618;
    Matrix<float, 2, 2>  _S622 = _S553 + _S621;
    Matrix<float, 2, 3>  _S623 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S624;
    (&_S624)->primal_0 = _S548;
    (&_S624)->differential_0 = _S623;
    Matrix<float, 3, 2>  _S625 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S626;
    (&_S626)->primal_0 = _S549;
    (&_S626)->differential_0 = _S625;
    s_bwd_prop_mul_0(&_S624, &_S626, _S622);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S627;
    (&_S627)->primal_0 = J_12;
    (&_S627)->differential_0 = _S623;
    Matrix<float, 3, 3>  _S628 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S629;
    (&_S629)->primal_0 = _S547;
    (&_S629)->differential_0 = _S628;
    s_bwd_prop_mul_1(&_S627, &_S629, _S624.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S630;
    (&_S630)->primal_0 = _S545;
    (&_S630)->differential_0 = _S628;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S631;
    (&_S631)->primal_0 = _S546;
    (&_S631)->differential_0 = _S628;
    s_bwd_prop_mul_2(&_S630, &_S631, _S629.differential_0);
    Matrix<float, 3, 3>  _S632 = transpose_0(_S631.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S633;
    (&_S633)->primal_0 = R_10;
    (&_S633)->differential_0 = _S628;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S634;
    (&_S634)->primal_0 = _S544;
    (&_S634)->differential_0 = _S628;
    s_bwd_prop_mul_2(&_S633, &_S634, _S630.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S635;
    (&_S635)->primal_0 = _S542;
    (&_S635)->differential_0 = _S628;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S636;
    (&_S636)->primal_0 = _S543;
    (&_S636)->differential_0 = _S628;
    s_bwd_prop_mul_2(&_S635, &_S636, _S634.differential_0);
    Matrix<float, 3, 3>  _S637 = _S635.differential_0 + transpose_0(_S636.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S541;
    (&_S638)->differential_0 = _S628;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S639;
    (&_S639)->primal_0 = S_2;
    (&_S639)->differential_0 = _S628;
    s_bwd_prop_mul_2(&_S638, &_S639, _S637);
    Matrix<float, 3, 3>  _S640 = transpose_0(_S638.differential_0);
    float _S641 = 2.0f * - _S640.rows[int(2)].z;
    float _S642 = 2.0f * _S640.rows[int(2)].y;
    float _S643 = 2.0f * _S640.rows[int(2)].x;
    float _S644 = 2.0f * _S640.rows[int(1)].z;
    float _S645 = 2.0f * - _S640.rows[int(1)].y;
    float _S646 = 2.0f * _S640.rows[int(1)].x;
    float _S647 = 2.0f * _S640.rows[int(0)].z;
    float _S648 = 2.0f * _S640.rows[int(0)].y;
    float _S649 = 2.0f * - _S640.rows[int(0)].x;
    float _S650 = - _S646 + _S648;
    float _S651 = _S643 + - _S647;
    float _S652 = - _S642 + _S644;
    float _S653 = _S642 + _S644;
    float _S654 = _S643 + _S647;
    float _S655 = _S646 + _S648;
    float _S656 = z_11 * (_S645 + _S649);
    float _S657 = y_12 * (_S641 + _S649);
    float _S658 = x_26 * (_S641 + _S645);
    float _S659 = z_11 * _S650 + y_12 * _S651 + x_26 * _S652;
    float _S660 = _S540 * _S659;
    float _S661 = w_11 * _S650 + y_12 * _S653 + x_26 * _S654 + _S656 + _S656;
    float _S662 = _S540 * _S661;
    float _S663 = w_11 * _S651 + z_11 * _S653 + x_26 * _S655 + _S657 + _S657;
    float _S664 = _S540 * _S663;
    float _S665 = w_11 * _S652 + z_11 * _S654 + y_12 * _S655 + _S658 + _S658;
    float _S666 = _S540 * _S665;
    float _S667 = quat_11.x * _S659 + quat_11.w * _S661 + quat_11.z * _S663 + quat_11.y * _S665;
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = _S539;
    (&_S668)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S668, _S667);
    float _S669 = quat_11.x * _S668.differential_0;
    float _S670 = quat_11.w * _S668.differential_0;
    float _S671 = quat_11.z * _S668.differential_0;
    float _S672 = quat_11.y * _S668.differential_0;
    float _S673 = _S662 + _S670 + _S670;
    float _S674 = _S664 + _S671 + _S671;
    float _S675 = _S666 + _S672 + _S672;
    float _S676 = _S660 + _S669 + _S669;
    float3  _S677 = make_float3 (0.0f);
    float3  _S678 = _S677;
    *&((&_S678)->z) = _S639.differential_0.rows[int(2)].z;
    *&((&_S678)->y) = _S639.differential_0.rows[int(1)].y;
    *&((&_S678)->x) = _S639.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S679;
    (&_S679)->primal_0 = scale_10;
    (&_S679)->differential_0 = _S677;
    s_bwd_prop_exp_1(&_S679, _S678);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S680 = _S679;
    float3  _S681 = _S677;
    *&((&_S681)->y) = _S619;
    *&((&_S681)->x) = _S620;
    float3  _S682 = _S587 + _S681;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S683;
    (&_S683)->primal_0 = R_10;
    (&_S683)->differential_0 = _S628;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S684;
    (&_S684)->primal_0 = mean_8;
    (&_S684)->differential_0 = _S677;
    s_bwd_prop_mul_3(&_S683, &_S684, _S682);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S685 = _S684;
    Matrix<float, 3, 3>  _S686 = _S632 + _S633.differential_0 + _S683.differential_0;
    float4  _S687 = make_float4 (0.0f);
    *&((&_S687)->w) = _S673;
    *&((&_S687)->z) = _S674;
    *&((&_S687)->y) = _S675;
    *&((&_S687)->x) = _S676;
    float4  _S688 = _S687;
    *v_mean_2 = _S685.differential_0;
    *v_quat_2 = _S688;
    *v_scale_2 = _S680.differential_0;
    *v_in_opacity_2 = _S594;
    *v_R_2 = _S686;
    *v_t_2 = _S682;
    return;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_9, float4  quat_12, float3  scale_11, float hardness_0, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_13, float fy_13, float cx_13, float cy_13, uint image_width_9, uint image_height_9, float near_plane_6, float far_plane_6, float radius_clip_6, int2  * radii_6, float * depth_6, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float * out_hardness_0)
{
    for(;;)
    {
        float _S689 = (mul_0(R_11, mean_9) + t_10).z;
        bool _S690;
        if(_S689 < near_plane_6)
        {
            _S690 = true;
        }
        else
        {
            _S690 = _S689 > far_plane_6;
        }
        if(_S690)
        {
            *radii_6 = make_int2 (int(0), int(0));
            break;
        }
        float3  _S691 = exp_0(scale_11);
        float x_27 = quat_12.y;
        float inv_norm_9 = (F32_rsqrt((x_27 * x_27 + quat_12.z * quat_12.z + quat_12.w * quat_12.w + quat_12.x * quat_12.x)));
        float x_28 = quat_12.y * inv_norm_9;
        float y_13 = quat_12.z * inv_norm_9;
        float z_12 = quat_12.w * inv_norm_9;
        float w_12 = quat_12.x * inv_norm_9;
        float x2_16 = x_28 * x_28;
        float y2_16 = y_13 * y_13;
        float z2_12 = z_12 * z_12;
        float xy_16 = x_28 * y_13;
        float xz_12 = x_28 * z_12;
        float yz_12 = y_13 * z_12;
        float wx_12 = w_12 * x_28;
        float wy_12 = w_12 * y_13;
        float wz_12 = w_12 * z_12;
        Matrix<float, 3, 3>  M_8 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_12), 2.0f * (xy_16 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_16 - wz_12), 1.0f - 2.0f * (x2_16 + z2_12), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_16 + y2_16))), makeMatrix<float, 3, 3> (_S691.x, 0.0f, 0.0f, 0.0f, _S691.y, 0.0f, 0.0f, 0.0f, _S691.z));
        float _S692 = (F32_sqrt((0.75f)));
        float3  vert0_c_0 = mul_0(R_11, mul_0(M_8, make_float3 (1.0f, 0.0f, 0.0f)) + mean_9) + t_10;
        float3  vert1_c_0 = mul_0(R_11, mul_0(M_8, make_float3 (-0.5f, _S692, 0.0f)) + mean_9) + t_10;
        float3  vert2_c_0 = mul_0(R_11, mul_0(M_8, make_float3 (-0.5f, - _S692, 0.0f)) + mean_9) + t_10;
        float _S693 = vert0_c_0.z;
        if(_S693 < near_plane_6)
        {
            _S690 = true;
        }
        else
        {
            _S690 = _S693 > far_plane_6;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = (vert1_c_0.z) < near_plane_6;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = (vert1_c_0.z) > far_plane_6;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = (vert2_c_0.z) < near_plane_6;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = (vert2_c_0.z) > far_plane_6;
        }
        if(_S690)
        {
            *radii_6 = make_int2 (int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S693);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (vert1_c_0.z);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (vert2_c_0.z);
        float2  _S694 = make_float2 (fx_13, fy_13);
        float2  _S695 = make_float2 (cx_13, cy_13);
        *uv0_0 = _S694 * *uv0_0 + _S695;
        *uv1_0 = _S694 * *uv1_0 + _S695;
        float2  _S696 = _S694 * *uv2_0 + _S695;
        *uv2_0 = _S696;
        float _S697 = _S696.x;
        float _S698 = float(image_width_9);
        float x_max_0 = (F32_min(((F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S697))) + 0.5f), (_S698)));
        float x_min_0 = (F32_max(((F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S697))) - 0.5f), (0.0f)));
        float _S699 = _S696.y;
        float _S700 = float(image_height_9);
        float y_max_0 = (F32_min(((F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S699))) + 0.5f), (_S700)));
        float y_min_0 = (F32_max(((F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S699))) - 0.5f), (0.0f)));
        float radius_x_3 = (F32_ceil((x_max_0 - x_min_0)));
        float radius_y_3 = (F32_ceil((y_max_0 - y_min_0)));
        if(radius_x_3 <= radius_clip_6)
        {
            _S690 = radius_y_3 <= radius_clip_6;
        }
        else
        {
            _S690 = false;
        }
        if(_S690)
        {
            *radii_6 = make_int2 (int(0), int(0));
            break;
        }
        if(x_max_0 <= 0.0f)
        {
            _S690 = true;
        }
        else
        {
            _S690 = x_min_0 >= _S698;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = y_max_0 <= 0.0f;
        }
        if(_S690)
        {
            _S690 = true;
        }
        else
        {
            _S690 = y_min_0 >= _S700;
        }
        if(_S690)
        {
            *radii_6 = make_int2 (int(0), int(0));
            break;
        }
        *radii_6 = make_int2 (int(radius_x_3), int(radius_y_3));
        *depth_6 = _S689;
        *out_hardness_0 = hardness_0;
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_10, float4  quat_13, float3  scale_12, float hardness_1, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_14, float fy_14, float cx_14, float cy_14, uint image_width_10, uint image_height_10, float near_plane_7, float far_plane_7, float radius_clip_7, int2  * radii_7, float * depth_7, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float * out_hardness_1)
{
    float3  mean_c_8 = mul_0(R_12, mean_10) + t_11;
    float3  _S701 = exp_0(scale_12);
    float x_29 = quat_13.y;
    float inv_norm_10 = (F32_rsqrt((x_29 * x_29 + quat_13.z * quat_13.z + quat_13.w * quat_13.w + quat_13.x * quat_13.x)));
    float x_30 = quat_13.y * inv_norm_10;
    float y_14 = quat_13.z * inv_norm_10;
    float z_13 = quat_13.w * inv_norm_10;
    float w_13 = quat_13.x * inv_norm_10;
    float x2_17 = x_30 * x_30;
    float y2_17 = y_14 * y_14;
    float z2_13 = z_13 * z_13;
    float xy_17 = x_30 * y_14;
    float xz_13 = x_30 * z_13;
    float yz_13 = y_14 * z_13;
    float wx_13 = w_13 * x_30;
    float wy_13 = w_13 * y_14;
    float wz_13 = w_13 * z_13;
    Matrix<float, 3, 3>  M_9 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_13), 2.0f * (xy_17 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_17 - wz_13), 1.0f - 2.0f * (x2_17 + z2_13), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_17 + y2_17))), makeMatrix<float, 3, 3> (_S701.x, 0.0f, 0.0f, 0.0f, _S701.y, 0.0f, 0.0f, 0.0f, _S701.z));
    float _S702 = (F32_sqrt((0.75f)));
    float3  vert0_c_1 = mul_0(R_12, mul_0(M_9, make_float3 (1.0f, 0.0f, 0.0f)) + mean_10) + t_11;
    float3  vert1_c_1 = mul_0(R_12, mul_0(M_9, make_float3 (-0.5f, _S702, 0.0f)) + mean_10) + t_11;
    float3  vert2_c_1 = mul_0(R_12, mul_0(M_9, make_float3 (-0.5f, - _S702, 0.0f)) + mean_10) + t_11;
    *uv0_1 = float2 {vert0_c_1.x, vert0_c_1.y} / make_float2 (vert0_c_1.z);
    *uv1_1 = float2 {vert1_c_1.x, vert1_c_1.y} / make_float2 (vert1_c_1.z);
    *uv2_1 = float2 {vert2_c_1.x, vert2_c_1.y} / make_float2 (vert2_c_1.z);
    float2  _S703 = make_float2 (fx_14, fy_14);
    float2  _S704 = make_float2 (cx_14, cy_14);
    *uv0_1 = _S703 * *uv0_1 + _S704;
    *uv1_1 = _S703 * *uv1_1 + _S704;
    float2  _S705 = _S703 * *uv2_1 + _S704;
    *uv2_1 = _S705;
    float _S706 = _S705.x;
    float _S707 = _S705.y;
    *radii_7 = make_int2 (int((F32_ceil(((F32_min(((F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S706))) + 0.5f), (float(image_width_10)))) - (F32_max(((F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S706))) - 0.5f), (0.0f))))))), int((F32_ceil(((F32_min(((F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S707))) + 0.5f), (float(image_height_10)))) - (F32_max(((F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S707))) - 0.5f), (0.0f))))))));
    *depth_7 = mean_c_8.z;
    *out_hardness_1 = hardness_1;
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_11, float4  quat_14, float3  scale_13, float hardness_2, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_15, float fy_15, float cx_15, float cy_15, uint image_width_11, uint image_height_11, float v_depth_3, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float v_out_hardness_0, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_hardness_0, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  _S708 = s_primal_ctx_exp_0(scale_13);
    float _S709 = quat_14.y;
    float _S710 = _S709 * _S709 + quat_14.z * quat_14.z + quat_14.w * quat_14.w + quat_14.x * quat_14.x;
    float _S711 = s_primal_ctx_rsqrt_0(_S710);
    float x_31 = quat_14.y * _S711;
    float y_15 = quat_14.z * _S711;
    float z_14 = quat_14.w * _S711;
    float w_14 = quat_14.x * _S711;
    float x2_18 = x_31 * x_31;
    float y2_18 = y_15 * y_15;
    float z2_14 = z_14 * z_14;
    float xy_18 = x_31 * y_15;
    float xz_14 = x_31 * z_14;
    float yz_14 = y_15 * z_14;
    float wx_14 = w_14 * x_31;
    float wy_14 = w_14 * y_15;
    float wz_14 = w_14 * z_14;
    Matrix<float, 3, 3>  _S712 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_18 + z2_14), 2.0f * (xy_18 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_18 - wz_14), 1.0f - 2.0f * (x2_18 + z2_14), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_18 + y2_18)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S708.x, 0.0f, 0.0f, 0.0f, _S708.y, 0.0f, 0.0f, 0.0f, _S708.z);
    Matrix<float, 3, 3>  _S713 = s_primal_ctx_mul_1(_S712, S_3);
    float3  _S714 = make_float3 (1.0f, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_0(_S713, _S714) + mean_11;
    float _S715 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S716 = make_float3 (-0.5f, _S715, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_0(_S713, _S716) + mean_11;
    float3  _S717 = make_float3 (-0.5f, - _S715, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_0(_S713, _S717) + mean_11;
    float3  vert0_c_2 = s_primal_ctx_mul_0(R_13, vert0_0) + t_12;
    float3  vert1_c_2 = s_primal_ctx_mul_0(R_13, vert1_0) + t_12;
    float3  vert2_c_2 = s_primal_ctx_mul_0(R_13, vert2_0) + t_12;
    float2  _S718 = float2 {vert0_c_2.x, vert0_c_2.y};
    float _S719 = vert0_c_2.z;
    float2  _S720 = make_float2 (_S719);
    float2  _S721 = make_float2 (_S719 * _S719);
    float2  _S722 = float2 {vert1_c_2.x, vert1_c_2.y};
    float _S723 = vert1_c_2.z;
    float2  _S724 = make_float2 (_S723);
    float2  _S725 = make_float2 (_S723 * _S723);
    float2  _S726 = float2 {vert2_c_2.x, vert2_c_2.y};
    float _S727 = vert2_c_2.z;
    float2  _S728 = make_float2 (_S727);
    float2  _S729 = make_float2 (_S727 * _S727);
    float2  _S730 = make_float2 (fx_15, fy_15);
    float2  _S731 = make_float2 (cx_15, cy_15);
    float2  _S732 = _S730 * (_S718 / make_float2 (_S719)) + _S731;
    float2  _S733 = _S730 * (_S722 / make_float2 (_S723)) + _S731;
    float2  _S734 = _S730 * (_S726 / make_float2 (_S727)) + _S731;
    float _S735 = _S732.x;
    float _S736 = _S733.x;
    float _S737 = s_primal_ctx_max_0(_S735, _S736);
    float _S738 = _S734.x;
    float _S739 = s_primal_ctx_max_0(_S737, _S738) + 0.5f;
    float _S740 = float(image_width_11);
    float _S741 = s_primal_ctx_min_0(_S735, _S736);
    float _S742 = s_primal_ctx_min_0(_S741, _S738) - 0.5f;
    float _S743 = _S732.y;
    float _S744 = _S733.y;
    float _S745 = s_primal_ctx_max_0(_S743, _S744);
    float _S746 = _S734.y;
    float _S747 = s_primal_ctx_max_0(_S745, _S746) + 0.5f;
    float _S748 = float(image_height_11);
    float _S749 = s_primal_ctx_min_0(_S743, _S744);
    DiffPair_float_0 _S750;
    (&_S750)->primal_0 = s_primal_ctx_min_0(_S749, _S746) - 0.5f;
    (&_S750)->differential_0 = 0.0f;
    DiffPair_float_0 _S751;
    (&_S751)->primal_0 = 0.0f;
    (&_S751)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S750, &_S751, -0.0f);
    DiffPair_float_0 _S752;
    (&_S752)->primal_0 = _S749;
    (&_S752)->differential_0 = 0.0f;
    DiffPair_float_0 _S753;
    (&_S753)->primal_0 = _S746;
    (&_S753)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S752, &_S753, _S750.differential_0);
    DiffPair_float_0 _S754;
    (&_S754)->primal_0 = _S743;
    (&_S754)->differential_0 = 0.0f;
    DiffPair_float_0 _S755;
    (&_S755)->primal_0 = _S744;
    (&_S755)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S754, &_S755, _S752.differential_0);
    DiffPair_float_0 _S756;
    (&_S756)->primal_0 = _S747;
    (&_S756)->differential_0 = 0.0f;
    DiffPair_float_0 _S757;
    (&_S757)->primal_0 = _S748;
    (&_S757)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S756, &_S757, 0.0f);
    DiffPair_float_0 _S758;
    (&_S758)->primal_0 = _S745;
    (&_S758)->differential_0 = 0.0f;
    DiffPair_float_0 _S759;
    (&_S759)->primal_0 = _S746;
    (&_S759)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S758, &_S759, _S756.differential_0);
    float _S760 = _S753.differential_0 + _S759.differential_0;
    DiffPair_float_0 _S761;
    (&_S761)->primal_0 = _S743;
    (&_S761)->differential_0 = 0.0f;
    DiffPair_float_0 _S762;
    (&_S762)->primal_0 = _S744;
    (&_S762)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S761, &_S762, _S758.differential_0);
    float _S763 = _S755.differential_0 + _S762.differential_0;
    float _S764 = _S754.differential_0 + _S761.differential_0;
    DiffPair_float_0 _S765;
    (&_S765)->primal_0 = _S742;
    (&_S765)->differential_0 = 0.0f;
    DiffPair_float_0 _S766;
    (&_S766)->primal_0 = 0.0f;
    (&_S766)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S765, &_S766, -0.0f);
    DiffPair_float_0 _S767;
    (&_S767)->primal_0 = _S741;
    (&_S767)->differential_0 = 0.0f;
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = _S738;
    (&_S768)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S767, &_S768, _S765.differential_0);
    DiffPair_float_0 _S769;
    (&_S769)->primal_0 = _S735;
    (&_S769)->differential_0 = 0.0f;
    DiffPair_float_0 _S770;
    (&_S770)->primal_0 = _S736;
    (&_S770)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S769, &_S770, _S767.differential_0);
    DiffPair_float_0 _S771;
    (&_S771)->primal_0 = _S739;
    (&_S771)->differential_0 = 0.0f;
    DiffPair_float_0 _S772;
    (&_S772)->primal_0 = _S740;
    (&_S772)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S771, &_S772, 0.0f);
    DiffPair_float_0 _S773;
    (&_S773)->primal_0 = _S737;
    (&_S773)->differential_0 = 0.0f;
    DiffPair_float_0 _S774;
    (&_S774)->primal_0 = _S738;
    (&_S774)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S773, &_S774, _S771.differential_0);
    float _S775 = _S768.differential_0 + _S774.differential_0;
    DiffPair_float_0 _S776;
    (&_S776)->primal_0 = _S735;
    (&_S776)->differential_0 = 0.0f;
    DiffPair_float_0 _S777;
    (&_S777)->primal_0 = _S736;
    (&_S777)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S776, &_S777, _S773.differential_0);
    float2  _S778 = _S730 * (v_uv2_0 + make_float2 (_S775, _S760)) / _S729;
    float2  _S779 = _S726 * - _S778;
    float2  _S780 = _S728 * _S778;
    float2  _S781 = _S730 * (v_uv1_0 + make_float2 (_S770.differential_0 + _S777.differential_0, _S763)) / _S725;
    float2  _S782 = _S722 * - _S781;
    float2  _S783 = _S724 * _S781;
    float _S784 = _S782.x + _S782.y;
    float2  _S785 = _S730 * (v_uv0_0 + make_float2 (_S769.differential_0 + _S776.differential_0, _S764)) / _S721;
    float2  _S786 = _S718 * - _S785;
    float2  _S787 = _S720 * _S785;
    float _S788 = _S786.x + _S786.y;
    float3  s_diff_vert2_c_T_0 = make_float3 (_S780.x, _S780.y, _S779.x + _S779.y);
    Matrix<float, 3, 3>  _S789 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S790;
    (&_S790)->primal_0 = R_13;
    (&_S790)->differential_0 = _S789;
    float3  _S791 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S792;
    (&_S792)->primal_0 = vert2_0;
    (&_S792)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S790, &_S792, s_diff_vert2_c_T_0);
    float3  s_diff_vert1_c_T_0 = make_float3 (_S783.x, _S783.y, _S784);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S793;
    (&_S793)->primal_0 = R_13;
    (&_S793)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S794;
    (&_S794)->primal_0 = vert1_0;
    (&_S794)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S793, &_S794, s_diff_vert1_c_T_0);
    float3  s_diff_vert0_c_T_0 = make_float3 (_S787.x, _S787.y, _S788);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S795;
    (&_S795)->primal_0 = R_13;
    (&_S795)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S796;
    (&_S796)->primal_0 = vert0_0;
    (&_S796)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S795, &_S796, s_diff_vert0_c_T_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S797;
    (&_S797)->primal_0 = _S713;
    (&_S797)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S798;
    (&_S798)->primal_0 = _S717;
    (&_S798)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S797, &_S798, _S792.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S799;
    (&_S799)->primal_0 = _S713;
    (&_S799)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S800;
    (&_S800)->primal_0 = _S716;
    (&_S800)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S799, &_S800, _S794.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S801;
    (&_S801)->primal_0 = _S713;
    (&_S801)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S802;
    (&_S802)->primal_0 = _S714;
    (&_S802)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S801, &_S802, _S796.differential_0);
    Matrix<float, 3, 3>  _S803 = _S797.differential_0 + _S799.differential_0 + _S801.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S804;
    (&_S804)->primal_0 = _S712;
    (&_S804)->differential_0 = _S789;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S805;
    (&_S805)->primal_0 = S_3;
    (&_S805)->differential_0 = _S789;
    s_bwd_prop_mul_2(&_S804, &_S805, _S803);
    Matrix<float, 3, 3>  _S806 = transpose_0(_S804.differential_0);
    float _S807 = 2.0f * - _S806.rows[int(2)].z;
    float _S808 = 2.0f * _S806.rows[int(2)].y;
    float _S809 = 2.0f * _S806.rows[int(2)].x;
    float _S810 = 2.0f * _S806.rows[int(1)].z;
    float _S811 = 2.0f * - _S806.rows[int(1)].y;
    float _S812 = 2.0f * _S806.rows[int(1)].x;
    float _S813 = 2.0f * _S806.rows[int(0)].z;
    float _S814 = 2.0f * _S806.rows[int(0)].y;
    float _S815 = 2.0f * - _S806.rows[int(0)].x;
    float _S816 = - _S812 + _S814;
    float _S817 = _S809 + - _S813;
    float _S818 = - _S808 + _S810;
    float _S819 = _S808 + _S810;
    float _S820 = _S809 + _S813;
    float _S821 = _S812 + _S814;
    float _S822 = z_14 * (_S811 + _S815);
    float _S823 = y_15 * (_S807 + _S815);
    float _S824 = x_31 * (_S807 + _S811);
    float _S825 = z_14 * _S816 + y_15 * _S817 + x_31 * _S818;
    float _S826 = _S711 * _S825;
    float _S827 = w_14 * _S816 + y_15 * _S819 + x_31 * _S820 + _S822 + _S822;
    float _S828 = _S711 * _S827;
    float _S829 = w_14 * _S817 + z_14 * _S819 + x_31 * _S821 + _S823 + _S823;
    float _S830 = _S711 * _S829;
    float _S831 = w_14 * _S818 + z_14 * _S820 + y_15 * _S821 + _S824 + _S824;
    float _S832 = _S711 * _S831;
    float _S833 = quat_14.x * _S825 + quat_14.w * _S827 + quat_14.z * _S829 + quat_14.y * _S831;
    DiffPair_float_0 _S834;
    (&_S834)->primal_0 = _S710;
    (&_S834)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S834, _S833);
    float _S835 = quat_14.x * _S834.differential_0;
    float _S836 = quat_14.w * _S834.differential_0;
    float _S837 = quat_14.z * _S834.differential_0;
    float _S838 = quat_14.y * _S834.differential_0;
    float _S839 = _S828 + _S836 + _S836;
    float _S840 = _S830 + _S837 + _S837;
    float _S841 = _S832 + _S838 + _S838;
    float _S842 = _S826 + _S835 + _S835;
    float3  _S843 = _S791;
    *&((&_S843)->z) = _S805.differential_0.rows[int(2)].z;
    *&((&_S843)->y) = _S805.differential_0.rows[int(1)].y;
    *&((&_S843)->x) = _S805.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S844;
    (&_S844)->primal_0 = scale_13;
    (&_S844)->differential_0 = _S791;
    s_bwd_prop_exp_1(&_S844, _S843);
    float3  s_diff_mean_c_T_0 = make_float3 (0.0f, 0.0f, v_depth_3);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S845;
    (&_S845)->primal_0 = R_13;
    (&_S845)->differential_0 = _S789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S846;
    (&_S846)->primal_0 = mean_11;
    (&_S846)->differential_0 = _S791;
    s_bwd_prop_mul_3(&_S845, &_S846, s_diff_mean_c_T_0);
    float3  _S847 = s_diff_vert2_c_T_0 + s_diff_vert1_c_T_0 + s_diff_vert0_c_T_0 + s_diff_mean_c_T_0;
    Matrix<float, 3, 3>  _S848 = _S790.differential_0 + _S793.differential_0 + _S795.differential_0 + _S845.differential_0;
    float4  _S849 = make_float4 (0.0f);
    *&((&_S849)->w) = _S839;
    *&((&_S849)->z) = _S840;
    *&((&_S849)->y) = _S841;
    *&((&_S849)->x) = _S842;
    *v_mean_3 = _S792.differential_0 + _S794.differential_0 + _S796.differential_0 + _S846.differential_0;
    *v_quat_3 = _S849;
    *v_scale_3 = _S844.differential_0;
    *v_hardness_0 = v_out_hardness_0;
    *v_R_3 = _S848;
    *v_t_3 = _S847;
    return;
}

inline __device__ float2  normalize_0(float2  x_32)
{
    return x_32 / make_float2 (length_0(x_32));
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_10, float dOut_13)
{
    float _S850 = _slang_select(((*dpx_10).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_10).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_13;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S850;
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float hardness_3, float2  p_0)
{
    float2  e0_0 = v1_0 - v0_0;
    float2  e1_0 = v2_0 - v1_0;
    float2  e2_0 = v0_0 - v2_0;
    float _S851 = e0_0.x * e1_0.y - e0_0.y * e1_0.x;
    float se_0 = float((F32_sign((_S851))));
    float2  _S852 = p_0 - v0_0;
    float2  _S853 = normalize_0(e0_0);
    float2  _S854 = p_0 - v1_0;
    float2  _S855 = normalize_0(e1_0);
    float2  _S856 = p_0 - v2_0;
    float2  _S857 = normalize_0(e2_0);
    return 0.99900001287460327f / (1.0f + (F32_exp((4.0f * ((F32_max(((F32_max((se_0 * (_S852.x * _S853.y - _S852.y * _S853.x)), (se_0 * (_S854.x * _S855.y - _S854.y * _S855.x))))), (se_0 * (_S856.x * _S857.y - _S856.y * _S857.x)))) / (2.0f * ((F32_abs((_S851))) / (length_0(e0_0) + length_0(e1_0) + length_0(e2_0))) * (1.0f - hardness_3) + hardness_3))))));
}

inline __device__ float s_primal_ctx_abs_0(float _S858)
{
    return (F32_abs((_S858)));
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S859, float _S860)
{
    _d_abs_0(_S859, _S860);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_11, float2  _s_dOut_1)
{
    float _S861 = length_0((*dpx_11).primal_0);
    float2  _S862 = (*dpx_11).primal_0 * _s_dOut_1;
    float2  _S863 = make_float2 (1.0f / _S861) * _s_dOut_1;
    float _S864 = - ((_S862.x + _S862.y) / (_S861 * _S861));
    float2  _S865 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S866;
    (&_S866)->primal_0 = (*dpx_11).primal_0;
    (&_S866)->differential_0 = _S865;
    s_bwd_length_impl_0(&_S866, _S864);
    float2  _S867 = _S863 + _S866.differential_0;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S867;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S868, float2  _S869)
{
    s_bwd_prop_normalize_impl_0(_S868, _S869);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_float_0 * dphardness_0, float2  p_1, float _s_dOut_2)
{
    float2  e0_1 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_1 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_1 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S870 = e0_1.x;
    float _S871 = e1_1.y;
    float _S872 = e0_1.y;
    float _S873 = e1_1.x;
    float _S874 = _S870 * _S871 - _S872 * _S873;
    float se_1 = float((F32_sign((_S874))));
    float2  _S875 = p_1 - (*dpv0_0).primal_0;
    float2  _S876 = normalize_0(e0_1);
    float _S877 = _S875.x;
    float _S878 = _S876.y;
    float _S879 = _S875.y;
    float _S880 = _S876.x;
    float de0_0 = se_1 * (_S877 * _S878 - _S879 * _S880);
    float2  _S881 = p_1 - (*dpv1_0).primal_0;
    float2  _S882 = normalize_0(e1_1);
    float _S883 = _S881.x;
    float _S884 = _S882.y;
    float _S885 = _S881.y;
    float _S886 = _S882.x;
    float de1_0 = se_1 * (_S883 * _S884 - _S885 * _S886);
    float2  _S887 = p_1 - (*dpv2_0).primal_0;
    float2  _S888 = normalize_0(e2_1);
    float _S889 = _S887.x;
    float _S890 = _S888.y;
    float _S891 = _S887.y;
    float _S892 = _S888.x;
    float de2_0 = se_1 * (_S889 * _S890 - _S891 * _S892);
    float _S893 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S894 = s_primal_ctx_max_0(_S893, de2_0);
    float _S895 = s_primal_ctx_abs_0(_S874);
    float _S896 = length_0(e0_1) + length_0(e1_1) + length_0(e2_1);
    float _S897 = _S896 * _S896;
    float _S898 = 2.0f * (_S895 / _S896);
    float _S899 = 1.0f - (*dphardness_0).primal_0;
    float _S900 = _S898 * _S899 + (*dphardness_0).primal_0;
    float _S901 = _S900 * _S900;
    float _S902 = 4.0f * (_S894 / _S900);
    float _S903 = 1.0f + s_primal_ctx_exp_1(_S902);
    float _S904 = 0.99900001287460327f * - (_s_dOut_2 / (_S903 * _S903));
    DiffPair_float_0 _S905;
    (&_S905)->primal_0 = _S902;
    (&_S905)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S905, _S904);
    float _S906 = 4.0f * _S905.differential_0 / _S901;
    float _S907 = _S894 * - _S906;
    float _S908 = _S900 * _S906;
    float _S909 = - (_S898 * _S907);
    float _S910 = 2.0f * (_S899 * _S907) / _S897;
    float _S911 = _S895 * - _S910;
    float _S912 = _S896 * _S910;
    float2  _S913 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S914;
    (&_S914)->primal_0 = e2_1;
    (&_S914)->differential_0 = _S913;
    s_bwd_length_impl_0(&_S914, _S911);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S915;
    (&_S915)->primal_0 = e1_1;
    (&_S915)->differential_0 = _S913;
    s_bwd_length_impl_0(&_S915, _S911);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S916;
    (&_S916)->primal_0 = e0_1;
    (&_S916)->differential_0 = _S913;
    s_bwd_length_impl_0(&_S916, _S911);
    DiffPair_float_0 _S917;
    (&_S917)->primal_0 = _S874;
    (&_S917)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S917, _S912);
    DiffPair_float_0 _S918;
    (&_S918)->primal_0 = _S893;
    (&_S918)->differential_0 = 0.0f;
    DiffPair_float_0 _S919;
    (&_S919)->primal_0 = de2_0;
    (&_S919)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S918, &_S919, _S908);
    DiffPair_float_0 _S920;
    (&_S920)->primal_0 = de0_0;
    (&_S920)->differential_0 = 0.0f;
    DiffPair_float_0 _S921;
    (&_S921)->primal_0 = de1_0;
    (&_S921)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S920, &_S921, _S918.differential_0);
    float _S922 = se_1 * _S919.differential_0;
    float _S923 = - _S922;
    float _S924 = _S892 * _S923;
    float _S925 = _S890 * _S922;
    float2  _S926 = make_float2 (_S891 * _S923, _S889 * _S922);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S927;
    (&_S927)->primal_0 = e2_1;
    (&_S927)->differential_0 = _S913;
    s_bwd_normalize_impl_0(&_S927, _S926);
    float2  _S928 = - make_float2 (_S925, _S924);
    float _S929 = se_1 * _S921.differential_0;
    float _S930 = - _S929;
    float _S931 = _S886 * _S930;
    float _S932 = _S884 * _S929;
    float2  _S933 = make_float2 (_S885 * _S930, _S883 * _S929);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S934;
    (&_S934)->primal_0 = e1_1;
    (&_S934)->differential_0 = _S913;
    s_bwd_normalize_impl_0(&_S934, _S933);
    float2  _S935 = - make_float2 (_S932, _S931);
    float _S936 = se_1 * _S920.differential_0;
    float _S937 = - _S936;
    float _S938 = _S880 * _S937;
    float _S939 = _S878 * _S936;
    float2  _S940 = make_float2 (_S879 * _S937, _S877 * _S936);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S941;
    (&_S941)->primal_0 = e0_1;
    (&_S941)->differential_0 = _S913;
    s_bwd_normalize_impl_0(&_S941, _S940);
    float2  _S942 = - make_float2 (_S939, _S938);
    float _S943 = - _S917.differential_0;
    float2  _S944 = _S914.differential_0 + _S927.differential_0;
    float2  _S945 = - _S944;
    float2  _S946 = _S915.differential_0 + _S934.differential_0 + make_float2 (_S872 * _S943, _S870 * _S917.differential_0);
    float2  _S947 = - _S946;
    float2  _S948 = _S916.differential_0 + _S941.differential_0 + make_float2 (_S871 * _S917.differential_0, _S873 * _S943);
    float2  _S949 = - _S948;
    float _S950 = _S907 + _S909;
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S950;
    float2  _S951 = _S928 + _S945 + _S946;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S951;
    float2  _S952 = _S935 + _S947 + _S948;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S952;
    float2  _S953 = _S942 + _S944 + _S949;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S953;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S954, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S955, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S956, DiffPair_float_0 * _S957, float2  _S958, float _S959)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S954, _S955, _S956, _S957, _S958, _S959);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float hardness_4, float2  p_2, float v_alpha_0, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float * v_hardness_1)
{
    float2  _S960 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S960;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S960;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S960;
    DiffPair_float_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_4;
    (&dp_hardness_0)->differential_0 = 0.0f;
    s_bwd_evaluate_alpha_opaque_triangle_fast_0(&dp_v0_0, &dp_v1_0, &dp_v2_0, &dp_hardness_0, p_2, v_alpha_0);
    *v_v0_0 = dp_v0_0.differential_0;
    *v_v1_0 = dp_v2_0.differential_0;
    *v_v2_0 = dp_v1_0.differential_0;
    *v_hardness_1 = dp_hardness_0.differential_0;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_12, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_14)
{
    DiffPair_float_0 _S961 = *dpx_12;
    bool _S962;
    if(((*dpx_12).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S962 = ((*dpx_12).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S962 = false;
    }
    float _S963;
    if(_S962)
    {
        _S963 = dOut_14;
    }
    else
    {
        _S963 = 0.0f;
    }
    dpx_12->primal_0 = _S961.primal_0;
    dpx_12->differential_0 = _S963;
    DiffPair_float_0 _S964 = *dpMin_0;
    if((_S961.primal_0) < ((*dpMin_0).primal_0))
    {
        _S963 = dOut_14;
    }
    else
    {
        _S963 = 0.0f;
    }
    dpMin_0->primal_0 = _S964.primal_0;
    dpMin_0->differential_0 = _S963;
    DiffPair_float_0 _S965 = *dpMax_0;
    if(((*dpx_12).primal_0) > ((*dpMax_0).primal_0))
    {
        _S963 = dOut_14;
    }
    else
    {
        _S963 = 0.0f;
    }
    dpMax_0->primal_0 = _S965.primal_0;
    dpMax_0->differential_0 = _S963;
    return;
}

inline __device__ float clamp_0(float x_33, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_33), (minBound_0)))), (maxBound_0)));
}

inline __device__ float evaluate_alpha_opaque_triangle_precise(float2  v0_2, float2  v1_2, float2  v2_2, float hardness_5, float2  p_3)
{
    float2  e0_2 = v1_2 - v0_2;
    float2  e1_2 = v2_2 - v1_2;
    float2  e2_2 = v0_2 - v2_2;
    float2  _S966 = p_3 - v0_2;
    float2  _S967 = p_3 - v1_2;
    float2  _S968 = p_3 - v2_2;
    float _S969 = e0_2.x;
    float _S970 = e1_2.y;
    float _S971 = e0_2.y;
    float _S972 = e1_2.x;
    float _S973 = _S969 * _S970 - _S971 * _S972;
    float se_2 = float((F32_sign((_S973))));
    return 0.99900001287460327f / (1.0f + (F32_exp((4.0f * (float((F32_sign(((F32_max(((F32_max((se_2 * (_S966.x * _S971 - _S966.y * _S969)), (se_2 * (_S967.x * _S970 - _S967.y * _S972))))), (se_2 * (_S968.x * e2_2.y - _S968.y * e2_2.x)))))))) * (F32_min(((F32_min((length_0(_S966 - e0_2 * make_float2 (clamp_0(dot_0(_S966, e0_2) / dot_0(e0_2, e0_2), 0.0f, 1.0f)))), (length_0(_S967 - e1_2 * make_float2 (clamp_0(dot_0(_S967, e1_2) / dot_0(e1_2, e1_2), 0.0f, 1.0f))))))), (length_0(_S968 - e2_2 * make_float2 (clamp_0(dot_0(_S968, e2_2) / dot_0(e2_2, e2_2), 0.0f, 1.0f)))))) / (2.0f * ((F32_abs((_S973))) / (length_0(e0_2) + length_0(e1_2) + length_0(e2_2))) * (1.0f - hardness_5) + hardness_5))))));
}

inline __device__ float s_primal_ctx_dot_0(float2  _S974, float2  _S975)
{
    return dot_0(_S974, _S975);
}

inline __device__ float s_primal_ctx_clamp_0(float _S976, float _S977, float _S978)
{
    return clamp_0(_S976, _S977, _S978);
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S979, DiffPair_float_0 * _S980, DiffPair_float_0 * _S981, float _S982)
{
    _d_clamp_0(_S979, _S980, _S981, _S982);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S983, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S984, float _S985)
{
    _d_dot_0(_S983, _S984, _S985);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_float_0 * dphardness_1, float2  p_4, float _s_dOut_3)
{
    float2  e0_3 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_3 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_3 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S986 = p_4 - (*dpv0_1).primal_0;
    float _S987 = s_primal_ctx_dot_0(_S986, e0_3);
    float _S988 = s_primal_ctx_dot_0(e0_3, e0_3);
    float _S989 = _S987 / _S988;
    float _S990 = _S988 * _S988;
    float _S991 = s_primal_ctx_clamp_0(_S989, 0.0f, 1.0f);
    float2  _S992 = make_float2 (_S991);
    float2  _S993 = _S986 - e0_3 * make_float2 (_S991);
    float _S994 = length_0(_S993);
    float2  _S995 = p_4 - (*dpv1_1).primal_0;
    float _S996 = s_primal_ctx_dot_0(_S995, e1_3);
    float _S997 = s_primal_ctx_dot_0(e1_3, e1_3);
    float _S998 = _S996 / _S997;
    float _S999 = _S997 * _S997;
    float _S1000 = s_primal_ctx_clamp_0(_S998, 0.0f, 1.0f);
    float2  _S1001 = make_float2 (_S1000);
    float2  _S1002 = _S995 - e1_3 * make_float2 (_S1000);
    float _S1003 = length_0(_S1002);
    float2  _S1004 = p_4 - (*dpv2_1).primal_0;
    float _S1005 = s_primal_ctx_dot_0(_S1004, e2_3);
    float _S1006 = s_primal_ctx_dot_0(e2_3, e2_3);
    float _S1007 = _S1005 / _S1006;
    float _S1008 = _S1006 * _S1006;
    float _S1009 = s_primal_ctx_clamp_0(_S1007, 0.0f, 1.0f);
    float2  _S1010 = make_float2 (_S1009);
    float2  _S1011 = _S1004 - e2_3 * make_float2 (_S1009);
    float _S1012 = length_0(_S1011);
    float _S1013 = e0_3.x;
    float _S1014 = e1_3.y;
    float _S1015 = e0_3.y;
    float _S1016 = e1_3.x;
    float _S1017 = _S1013 * _S1014 - _S1015 * _S1016;
    float se_3 = float((F32_sign((_S1017))));
    float _S1018 = _S986.x;
    float _S1019 = _S986.y;
    float s0_0 = se_3 * (_S1018 * _S1015 - _S1019 * _S1013);
    float _S1020 = _S995.x;
    float _S1021 = _S995.y;
    float s1_0 = se_3 * (_S1020 * _S1014 - _S1021 * _S1016);
    float _S1022 = _S1004.x;
    float _S1023 = e2_3.y;
    float _S1024 = _S1004.y;
    float _S1025 = e2_3.x;
    float s2_0 = se_3 * (_S1022 * _S1023 - _S1024 * _S1025);
    float _S1026 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S1026, s2_0)))));
    float _S1027 = s_primal_ctx_min_0(_S994, _S1003);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S1027, _S1012);
    float _S1028 = s_primal_ctx_abs_0(_S1017);
    float _S1029 = length_0(e0_3) + length_0(e1_3) + length_0(e2_3);
    float _S1030 = _S1029 * _S1029;
    float _S1031 = 2.0f * (_S1028 / _S1029);
    float _S1032 = 1.0f - (*dphardness_1).primal_0;
    float _S1033 = _S1031 * _S1032 + (*dphardness_1).primal_0;
    float _S1034 = _S1033 * _S1033;
    float _S1035 = 4.0f * (dv_0 / _S1033);
    float _S1036 = 1.0f + s_primal_ctx_exp_1(_S1035);
    float _S1037 = 0.99900001287460327f * - (_s_dOut_3 / (_S1036 * _S1036));
    DiffPair_float_0 _S1038;
    (&_S1038)->primal_0 = _S1035;
    (&_S1038)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1038, _S1037);
    float _S1039 = 4.0f * _S1038.differential_0 / _S1034;
    float _S1040 = dv_0 * - _S1039;
    float s_diff_dv_T_0 = _S1033 * _S1039;
    float _S1041 = - (_S1031 * _S1040);
    float _S1042 = 2.0f * (_S1032 * _S1040) / _S1030;
    float _S1043 = _S1028 * - _S1042;
    float _S1044 = _S1029 * _S1042;
    float2  _S1045 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1046;
    (&_S1046)->primal_0 = e2_3;
    (&_S1046)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1046, _S1043);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1047;
    (&_S1047)->primal_0 = e1_3;
    (&_S1047)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1047, _S1043);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1048;
    (&_S1048)->primal_0 = e0_3;
    (&_S1048)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1048, _S1043);
    DiffPair_float_0 _S1049;
    (&_S1049)->primal_0 = _S1017;
    (&_S1049)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1049, _S1044);
    float _S1050 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S1051;
    (&_S1051)->primal_0 = _S1027;
    (&_S1051)->differential_0 = 0.0f;
    DiffPair_float_0 _S1052;
    (&_S1052)->primal_0 = _S1012;
    (&_S1052)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1051, &_S1052, _S1050);
    DiffPair_float_0 _S1053;
    (&_S1053)->primal_0 = _S994;
    (&_S1053)->differential_0 = 0.0f;
    DiffPair_float_0 _S1054;
    (&_S1054)->primal_0 = _S1003;
    (&_S1054)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1053, &_S1054, _S1051.differential_0);
    DiffPair_float_0 _S1055;
    (&_S1055)->primal_0 = _S1026;
    (&_S1055)->differential_0 = 0.0f;
    DiffPair_float_0 _S1056;
    (&_S1056)->primal_0 = s2_0;
    (&_S1056)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1055, &_S1056, 0.0f);
    DiffPair_float_0 _S1057;
    (&_S1057)->primal_0 = s0_0;
    (&_S1057)->differential_0 = 0.0f;
    DiffPair_float_0 _S1058;
    (&_S1058)->primal_0 = s1_0;
    (&_S1058)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1057, &_S1058, _S1055.differential_0);
    float _S1059 = se_3 * _S1056.differential_0;
    float _S1060 = - _S1059;
    float _S1061 = _S1024 * _S1060;
    float _S1062 = _S1025 * _S1060;
    float _S1063 = _S1022 * _S1059;
    float _S1064 = _S1023 * _S1059;
    float _S1065 = se_3 * _S1058.differential_0;
    float _S1066 = - _S1065;
    float _S1067 = _S1016 * _S1066;
    float _S1068 = _S1014 * _S1065;
    float _S1069 = se_3 * _S1057.differential_0;
    float _S1070 = - _S1069;
    float _S1071 = _S1013 * _S1070;
    float _S1072 = _S1015 * _S1069;
    float _S1073 = - _S1049.differential_0;
    float _S1074 = _S1021 * _S1066 + _S1015 * _S1073;
    float _S1075 = _S1018 * _S1069 + _S1016 * _S1073;
    float _S1076 = _S1020 * _S1065 + _S1013 * _S1049.differential_0;
    float _S1077 = _S1019 * _S1070 + _S1014 * _S1049.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1078;
    (&_S1078)->primal_0 = _S1011;
    (&_S1078)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1078, _S1052.differential_0);
    float2  _S1079 = - _S1078.differential_0;
    float2  _S1080 = e2_3 * _S1079;
    float2  _S1081 = _S1010 * _S1079;
    float _S1082 = _S1080.x + _S1080.y;
    DiffPair_float_0 _S1083;
    (&_S1083)->primal_0 = _S1007;
    (&_S1083)->differential_0 = 0.0f;
    DiffPair_float_0 _S1084;
    (&_S1084)->primal_0 = 0.0f;
    (&_S1084)->differential_0 = 0.0f;
    DiffPair_float_0 _S1085;
    (&_S1085)->primal_0 = 1.0f;
    (&_S1085)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1083, &_S1084, &_S1085, _S1082);
    float _S1086 = _S1083.differential_0 / _S1008;
    float _S1087 = _S1005 * - _S1086;
    float _S1088 = _S1006 * _S1086;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1089;
    (&_S1089)->primal_0 = e2_3;
    (&_S1089)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1090;
    (&_S1090)->primal_0 = e2_3;
    (&_S1090)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1089, &_S1090, _S1087);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1091;
    (&_S1091)->primal_0 = _S1004;
    (&_S1091)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1092;
    (&_S1092)->primal_0 = e2_3;
    (&_S1092)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1091, &_S1092, _S1088);
    float2  _S1093 = - (_S1078.differential_0 + _S1091.differential_0 + make_float2 (_S1064, _S1062));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1094;
    (&_S1094)->primal_0 = _S1002;
    (&_S1094)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1094, _S1054.differential_0);
    float2  _S1095 = - _S1094.differential_0;
    float2  _S1096 = e1_3 * _S1095;
    float2  _S1097 = _S1001 * _S1095;
    float _S1098 = _S1096.x + _S1096.y;
    DiffPair_float_0 _S1099;
    (&_S1099)->primal_0 = _S998;
    (&_S1099)->differential_0 = 0.0f;
    DiffPair_float_0 _S1100;
    (&_S1100)->primal_0 = 0.0f;
    (&_S1100)->differential_0 = 0.0f;
    DiffPair_float_0 _S1101;
    (&_S1101)->primal_0 = 1.0f;
    (&_S1101)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1099, &_S1100, &_S1101, _S1098);
    float _S1102 = _S1099.differential_0 / _S999;
    float _S1103 = _S996 * - _S1102;
    float _S1104 = _S997 * _S1102;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1105;
    (&_S1105)->primal_0 = e1_3;
    (&_S1105)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1106;
    (&_S1106)->primal_0 = e1_3;
    (&_S1106)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1105, &_S1106, _S1103);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1107;
    (&_S1107)->primal_0 = _S995;
    (&_S1107)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1108;
    (&_S1108)->primal_0 = e1_3;
    (&_S1108)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1107, &_S1108, _S1104);
    float2  _S1109 = - (_S1094.differential_0 + _S1107.differential_0 + make_float2 (_S1068, _S1067));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1110;
    (&_S1110)->primal_0 = _S993;
    (&_S1110)->differential_0 = _S1045;
    s_bwd_length_impl_0(&_S1110, _S1053.differential_0);
    float2  _S1111 = - _S1110.differential_0;
    float2  _S1112 = e0_3 * _S1111;
    float2  _S1113 = _S992 * _S1111;
    float _S1114 = _S1112.x + _S1112.y;
    DiffPair_float_0 _S1115;
    (&_S1115)->primal_0 = _S989;
    (&_S1115)->differential_0 = 0.0f;
    DiffPair_float_0 _S1116;
    (&_S1116)->primal_0 = 0.0f;
    (&_S1116)->differential_0 = 0.0f;
    DiffPair_float_0 _S1117;
    (&_S1117)->primal_0 = 1.0f;
    (&_S1117)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S1115, &_S1116, &_S1117, _S1114);
    float _S1118 = _S1115.differential_0 / _S990;
    float _S1119 = _S987 * - _S1118;
    float _S1120 = _S988 * _S1118;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1121;
    (&_S1121)->primal_0 = e0_3;
    (&_S1121)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1122;
    (&_S1122)->primal_0 = e0_3;
    (&_S1122)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1121, &_S1122, _S1119);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1123;
    (&_S1123)->primal_0 = _S986;
    (&_S1123)->differential_0 = _S1045;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1124;
    (&_S1124)->primal_0 = e0_3;
    (&_S1124)->differential_0 = _S1045;
    s_bwd_prop_dot_0(&_S1123, &_S1124, _S1120);
    float2  _S1125 = - (_S1110.differential_0 + _S1123.differential_0 + make_float2 (_S1072, _S1071));
    float2  _S1126 = _S1046.differential_0 + _S1081 + _S1090.differential_0 + _S1089.differential_0 + _S1092.differential_0 + make_float2 (_S1061, _S1063);
    float2  _S1127 = - _S1126;
    float2  _S1128 = _S1047.differential_0 + _S1097 + _S1106.differential_0 + _S1105.differential_0 + _S1108.differential_0 + make_float2 (_S1074, _S1076);
    float2  _S1129 = - _S1128;
    float2  _S1130 = _S1048.differential_0 + _S1113 + _S1122.differential_0 + _S1121.differential_0 + _S1124.differential_0 + make_float2 (_S1077, _S1075);
    float2  _S1131 = - _S1130;
    float _S1132 = _S1040 + _S1041;
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S1132;
    float2  _S1133 = _S1093 + _S1127 + _S1128;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S1133;
    float2  _S1134 = _S1109 + _S1129 + _S1130;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S1134;
    float2  _S1135 = _S1125 + _S1126 + _S1131;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S1135;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1136, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1137, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1138, DiffPair_float_0 * _S1139, float2  _S1140, float _S1141)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S1136, _S1137, _S1138, _S1139, _S1140, _S1141);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float hardness_6, float2  p_5, float v_alpha_1, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float * v_hardness_2)
{
    float2  _S1142 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S1142;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S1142;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S1142;
    DiffPair_float_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_6;
    (&dp_hardness_1)->differential_0 = 0.0f;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_1);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_2 = dp_hardness_1.differential_0;
    return;
}

