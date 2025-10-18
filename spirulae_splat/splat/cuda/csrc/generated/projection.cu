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

inline __device__ float dot_0(float2  x_7, float2  y_2)
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
        float result_8 = result_7 + _slang_vector_get_element(x_7, i_4) * _slang_vector_get_element(y_2, i_4);
        i_4 = i_4 + int(1);
        result_7 = result_8;
    }
    return result_7;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_0(x_8, x_8))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_2, DiffPair_float_0 * dpx_4, float dOut_8)
{
    DiffPair_float_0 _S16 = *dpx_4;
    float _S17 = - (*dpy_2).primal_0 / ((*dpx_4).primal_0 * (*dpx_4).primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S17;
    float _S18 = _S16.primal_0 / (_S16.primal_0 * _S16.primal_0 + (*dpy_2).primal_0 * (*dpy_2).primal_0) * dOut_8;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = _S18;
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    float xy_len_0 = length_0(make_float2 (mean3d_2.x, mean3d_2.y)) + 1.00000001168609742e-07f;
    float theta_0 = (F32_atan2((xy_len_0), (mean3d_2.z + 1.00000001168609742e-07f)));
    *mean2d_2 = make_float2 (mean3d_2.x * fx_2 * theta_0 / xy_len_0 + cx_2, mean3d_2.y * fy_2 * theta_0 / xy_len_0 + cy_2);
    float x2_2 = mean3d_2.x * mean3d_2.x + 1.00000001168609742e-07f;
    float y2_2 = mean3d_2.y * mean3d_2.y;
    float xy_2 = mean3d_2.x * mean3d_2.y;
    float x2y2_0 = x2_2 + y2_2;
    float x2y2z2_inv_0 = 1.0f / (x2y2_0 + mean3d_2.z * mean3d_2.z);
    float b_0 = (F32_atan2((xy_len_0), (mean3d_2.z))) / xy_len_0 / x2y2_0;
    float a_0 = mean3d_2.z * x2y2z2_inv_0 / x2y2_0;
    float _S19 = a_0 - b_0;
    Matrix<float, 2, 3>  J_2 = makeMatrix<float, 2, 3> (fx_2 * (x2_2 * a_0 + y2_2 * b_0), fx_2 * xy_2 * _S19, - fx_2 * mean3d_2.x * x2y2z2_inv_0, fy_2 * xy_2 * _S19, fy_2 * (y2_2 * a_0 + x2_2 * b_0), - fy_2 * mean3d_2.y * x2y2z2_inv_0);
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

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_9)
{
    float _S21 = (F32_exp(((*dpx_5).primal_0))) * dOut_9;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S21;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_6, float dOut_10)
{
    float _S22 = 1.0f / (*dpx_6).primal_0 * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S22;
    return;
}

inline __device__ void projection_3dgs_persp(float3  mean_0, float4  quat_2, float3  scale_1, float in_opacity_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_4, float fy_4, float cx_4, float cy_4, uint image_width_0, uint image_height_0, float eps2d_1, float near_plane_0, float far_plane_0, float radius_clip_0, int2  * radii_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_2, mean_0) + t_1;
        float _S23 = mean_c_0.z;
        bool _S24;
        if(_S23 < near_plane_0)
        {
            _S24 = true;
        }
        else
        {
            _S24 = _S23 > far_plane_0;
        }
        if(_S24)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        float x_9 = quat_2.y;
        float inv_norm_2 = (F32_rsqrt((x_9 * x_9 + quat_2.z * quat_2.z + quat_2.w * quat_2.w + quat_2.x * quat_2.x)));
        float x_10 = quat_2.y * inv_norm_2;
        float y_3 = quat_2.z * inv_norm_2;
        float z_2 = quat_2.w * inv_norm_2;
        float w_2 = quat_2.x * inv_norm_2;
        float x2_3 = x_10 * x_10;
        float y2_3 = y_3 * y_3;
        float z2_2 = z_2 * z_2;
        float xy_3 = x_10 * y_3;
        float xz_2 = x_10 * z_2;
        float yz_2 = y_3 * z_2;
        float wx_2 = w_2 * x_10;
        float wy_2 = w_2 * y_3;
        float wz_2 = w_2 * z_2;
        Matrix<float, 3, 3>  M_1 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_2), 2.0f * (xy_3 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_3 - wz_2), 1.0f - 2.0f * (x2_3 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_2, mul_4(M_1, transpose_0(M_1))), transpose_0(R_2));
        Matrix<float, 2, 2>  covar2d_0;
        float _S25 = float(image_width_0);
        float _S26 = float(image_height_0);
        float _S27 = 0.30000001192092896f * (0.5f * _S25 / fx_4);
        float _S28 = 0.30000001192092896f * (0.5f * _S26 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S25 - cx_4) / fx_4 + _S27), ((F32_max((- (cx_4 / fx_4 + _S27)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S26 - cy_4) / fy_4 + _S28), ((F32_max((- (cy_4 / fy_4 + _S28)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_4, covar_c_0), transpose_1(J_4));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        *&(((&covar2d_0)->rows + (int(0)))->x) = *&(((&covar2d_0)->rows + (int(0)))->x) + eps2d_1;
        float _S29 = *&(((&covar2d_0)->rows + (int(1)))->y) + eps2d_1;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S29;
        if((*&(((&covar2d_0)->rows + (int(0)))->x) * _S29 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x)) <= 0.0f)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S30 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
        *opacity_0 = 1.0f / (1.0f + (F32_exp((- in_opacity_0))));
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
            _S24 = radius_y_0 <= radius_clip_0;
        }
        else
        {
            _S24 = false;
        }
        if(_S24)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_4).x + radius_x_0) <= 0.0f)
        {
            _S24 = true;
        }
        else
        {
            _S24 = ((*mean2d_4).x - radius_x_0) >= _S25;
        }
        if(_S24)
        {
            _S24 = true;
        }
        else
        {
            _S24 = ((*mean2d_4).y + radius_y_0) <= 0.0f;
        }
        if(_S24)
        {
            _S24 = true;
        }
        else
        {
            _S24 = ((*mean2d_4).y - radius_y_0) >= _S26;
        }
        if(_S24)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        *radii_0 = make_int2 (int(radius_x_0), int(radius_y_0));
        *depth_0 = _S23;
        *conic_0 = make_float3 (_S30.rows[int(0)].x, _S30.rows[int(0)].y, _S30.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(float3  mean_1, float4  quat_3, float3  scale_2, float in_opacity_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_5, float fy_5, float cx_5, float cy_5, uint image_width_1, uint image_height_1, float eps2d_2, float near_plane_1, float far_plane_1, float radius_clip_1, int2  * radii_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_3, mean_1) + t_2;
        float _S31 = mean_c_1.z;
        bool _S32;
        if(_S31 < near_plane_1)
        {
            _S32 = true;
        }
        else
        {
            _S32 = _S31 > far_plane_1;
        }
        if(_S32)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        float x_11 = quat_3.y;
        float inv_norm_3 = (F32_rsqrt((x_11 * x_11 + quat_3.z * quat_3.z + quat_3.w * quat_3.w + quat_3.x * quat_3.x)));
        float x_12 = quat_3.y * inv_norm_3;
        float y_4 = quat_3.z * inv_norm_3;
        float z_3 = quat_3.w * inv_norm_3;
        float w_3 = quat_3.x * inv_norm_3;
        float x2_4 = x_12 * x_12;
        float y2_4 = y_4 * y_4;
        float z2_3 = z_3 * z_3;
        float xy_4 = x_12 * y_4;
        float xz_3 = x_12 * z_3;
        float yz_3 = y_4 * z_3;
        float wx_3 = w_3 * x_12;
        float wy_3 = w_3 * y_4;
        float wz_3 = w_3 * z_3;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_3), 2.0f * (xy_4 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_4 - wz_3), 1.0f - 2.0f * (x2_4 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (scale_2.x, 0.0f, 0.0f, 0.0f, scale_2.y, 0.0f, 0.0f, 0.0f, scale_2.z));
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_3, mul_4(M_2, transpose_0(M_2))), transpose_0(R_3));
        Matrix<float, 2, 2>  covar2d_1;
        float xy_len_1 = length_0(make_float2 (mean_c_1.x, mean_c_1.y)) + 1.00000001168609742e-07f;
        float theta_1 = (F32_atan2((xy_len_1), (mean_c_1.z + 1.00000001168609742e-07f)));
        *mean2d_5 = make_float2 (mean_c_1.x * fx_5 * theta_1 / xy_len_1 + cx_5, mean_c_1.y * fy_5 * theta_1 / xy_len_1 + cy_5);
        float x2_5 = mean_c_1.x * mean_c_1.x + 1.00000001168609742e-07f;
        float y2_5 = mean_c_1.y * mean_c_1.y;
        float xy_5 = mean_c_1.x * mean_c_1.y;
        float x2y2_1 = x2_5 + y2_5;
        float x2y2z2_inv_1 = 1.0f / (x2y2_1 + mean_c_1.z * mean_c_1.z);
        float b_1 = (F32_atan2((xy_len_1), (mean_c_1.z))) / xy_len_1 / x2y2_1;
        float a_1 = mean_c_1.z * x2y2z2_inv_1 / x2y2_1;
        float _S33 = a_1 - b_1;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_5 * (x2_5 * a_1 + y2_5 * b_1), fx_5 * xy_5 * _S33, - fx_5 * mean_c_1.x * x2y2z2_inv_1, fy_5 * xy_5 * _S33, fy_5 * (y2_5 * a_1 + x2_5 * b_1), - fy_5 * mean_c_1.y * x2y2z2_inv_1);
        covar2d_1 = mul_6(mul_5(J_5, covar_c_1), transpose_1(J_5));
        *&(((&covar2d_1)->rows + (int(0)))->x) = *&(((&covar2d_1)->rows + (int(0)))->x) + eps2d_2;
        float _S34 = *&(((&covar2d_1)->rows + (int(1)))->y) + eps2d_2;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S34;
        if((*&(((&covar2d_1)->rows + (int(0)))->x) * _S34 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x)) <= 0.0f)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S35 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
        *opacity_1 = 1.0f / (1.0f + (F32_exp((- in_opacity_1))));
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
            _S32 = radius_y_1 <= radius_clip_1;
        }
        else
        {
            _S32 = false;
        }
        if(_S32)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_5).x + radius_x_1) <= 0.0f)
        {
            _S32 = true;
        }
        else
        {
            _S32 = ((*mean2d_5).x - radius_x_1) >= float(image_width_1);
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = ((*mean2d_5).y + radius_y_1) <= 0.0f;
        }
        if(_S32)
        {
            _S32 = true;
        }
        else
        {
            _S32 = ((*mean2d_5).y - radius_y_1) >= float(image_height_1);
        }
        if(_S32)
        {
            *radii_1 = make_int2 (int(0), int(0));
            break;
        }
        *radii_1 = make_int2 (int(radius_x_1), int(radius_y_1));
        *depth_1 = _S31;
        *conic_1 = make_float3 (_S35.rows[int(0)].x, _S35.rows[int(0)].y, _S35.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(float3  mean_2, float4  quat_4, float3  scale_3, float in_opacity_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_6, float fy_6, float cx_6, float cy_6, uint image_width_2, uint image_height_2, float eps2d_3, float near_plane_2, float far_plane_2, float radius_clip_2, int2  * radii_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_4, mean_2) + t_3;
        float _S36 = mean_c_2.z;
        bool _S37;
        if(_S36 < near_plane_2)
        {
            _S37 = true;
        }
        else
        {
            _S37 = _S36 > far_plane_2;
        }
        if(_S37)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        float x_13 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_13 * x_13 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_14 = quat_4.y * inv_norm_4;
        float y_5 = quat_4.z * inv_norm_4;
        float z_4 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_6 = x_14 * x_14;
        float y2_6 = y_5 * y_5;
        float z2_4 = z_4 * z_4;
        float xy_6 = x_14 * y_5;
        float xz_4 = x_14 * z_4;
        float yz_4 = y_5 * z_4;
        float wx_4 = w_4 * x_14;
        float wy_4 = w_4 * y_5;
        float wz_4 = w_4 * z_4;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_4), 2.0f * (xy_6 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_6 - wz_4), 1.0f - 2.0f * (x2_6 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (scale_3.x, 0.0f, 0.0f, 0.0f, scale_3.y, 0.0f, 0.0f, 0.0f, scale_3.z));
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_4, mul_4(M_3, transpose_0(M_3))), transpose_0(R_4));
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        *&(((&covar2d_2)->rows + (int(0)))->x) = *&(((&covar2d_2)->rows + (int(0)))->x) + eps2d_3;
        float _S38 = *&(((&covar2d_2)->rows + (int(1)))->y) + eps2d_3;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S38;
        if((*&(((&covar2d_2)->rows + (int(0)))->x) * _S38 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x)) <= 0.0f)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S39 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
        *opacity_2 = 1.0f / (1.0f + (F32_exp((- in_opacity_2))));
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
            _S37 = radius_y_2 <= radius_clip_2;
        }
        else
        {
            _S37 = false;
        }
        if(_S37)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        if(((*mean2d_6).x + radius_x_2) <= 0.0f)
        {
            _S37 = true;
        }
        else
        {
            _S37 = ((*mean2d_6).x - radius_x_2) >= float(image_width_2);
        }
        if(_S37)
        {
            _S37 = true;
        }
        else
        {
            _S37 = ((*mean2d_6).y + radius_y_2) <= 0.0f;
        }
        if(_S37)
        {
            _S37 = true;
        }
        else
        {
            _S37 = ((*mean2d_6).y - radius_y_2) >= float(image_height_2);
        }
        if(_S37)
        {
            *radii_2 = make_int2 (int(0), int(0));
            break;
        }
        *radii_2 = make_int2 (int(radius_x_2), int(radius_y_2));
        *depth_2 = _S36;
        *conic_2 = make_float3 (_S39.rows[int(0)].x, _S39.rows[int(0)].y, _S39.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S40, float3  _S41)
{
    return mul_0(_S40, _S41);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S42)
{
    return (F32_rsqrt((_S42)));
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S43, Matrix<float, 3, 3>  _S44)
{
    return mul_4(_S43, _S44);
}

inline __device__ float s_primal_ctx_max_0(float _S45, float _S46)
{
    return (F32_max((_S45), (_S46)));
}

inline __device__ float s_primal_ctx_min_0(float _S47, float _S48)
{
    return (F32_min((_S47), (_S48)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S49, Matrix<float, 3, 3>  _S50)
{
    return mul_5(_S49, _S50);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S51, Matrix<float, 3, 2>  _S52)
{
    return mul_6(_S51, _S52);
}

inline __device__ float s_primal_ctx_exp_0(float _S53)
{
    return (F32_exp((_S53)));
}

inline __device__ float s_primal_ctx_log_0(float _S54)
{
    return (F32_log((_S54)));
}

inline __device__ float s_primal_ctx_sqrt_0(float _S55)
{
    return (F32_sqrt((_S55)));
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S56, float _S57)
{
    _d_sqrt_0(_S56, _S57);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S58, DiffPair_float_0 * _S59, float _S60)
{
    _d_min_0(_S58, _S59, _S60);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S61, float _S62)
{
    _d_log_0(_S61, _S62);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S63, float _S64)
{
    _d_exp_0(_S63, _S64);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S65, DiffPair_float_0 * _S66, float _S67)
{
    _d_max_0(_S65, _S66, _S67);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S68, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S69, Matrix<float, 2, 2>  _S70)
{
    mul_3(_S68, _S69, _S70);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S71, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S72, Matrix<float, 2, 3>  _S73)
{
    mul_2(_S71, _S72, _S73);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S74, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S75, Matrix<float, 3, 3>  _S76)
{
    mul_1(_S74, _S75, _S76);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S77, float _S78)
{
    _d_rsqrt_0(_S77, _S78);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S79, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S80, float3  _S81)
{
    _d_mul_0(_S79, _S80, _S81);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(float3  mean_3, float4  quat_5, float3  scale_4, float in_opacity_3, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_7, float fy_7, float cx_7, float cy_7, uint image_width_3, uint image_height_3, float eps2d_4, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_3 = s_primal_ctx_mul_0(R_5, mean_3) + t_4;
    float _S82 = quat_5.y;
    float _S83 = _S82 * _S82 + quat_5.z * quat_5.z + quat_5.w * quat_5.w + quat_5.x * quat_5.x;
    float _S84 = s_primal_ctx_rsqrt_0(_S83);
    float x_15 = quat_5.y * _S84;
    float y_6 = quat_5.z * _S84;
    float z_5 = quat_5.w * _S84;
    float w_5 = quat_5.x * _S84;
    float x2_7 = x_15 * x_15;
    float y2_7 = y_6 * y_6;
    float z2_5 = z_5 * z_5;
    float xy_7 = x_15 * y_6;
    float xz_5 = x_15 * z_5;
    float yz_5 = y_6 * z_5;
    float wx_5 = w_5 * x_15;
    float wy_5 = w_5 * y_6;
    float wz_5 = w_5 * z_5;
    Matrix<float, 3, 3>  _S85 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_5), 2.0f * (xy_7 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_7 - wz_5), 1.0f - 2.0f * (x2_7 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_7 + y2_7)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (scale_4.x, 0.0f, 0.0f, 0.0f, scale_4.y, 0.0f, 0.0f, 0.0f, scale_4.z);
    Matrix<float, 3, 3>  _S86 = s_primal_ctx_mul_1(_S85, S_0);
    Matrix<float, 3, 3>  _S87 = transpose_0(_S86);
    Matrix<float, 3, 3>  _S88 = s_primal_ctx_mul_1(_S86, _S87);
    Matrix<float, 3, 3>  _S89 = s_primal_ctx_mul_1(R_5, _S88);
    Matrix<float, 3, 3>  _S90 = transpose_0(R_5);
    Matrix<float, 3, 3>  _S91 = s_primal_ctx_mul_1(_S89, _S90);
    float _S92 = float(image_width_3);
    float _S93 = float(image_height_3);
    float _S94 = 0.30000001192092896f * (0.5f * _S92 / fx_7);
    float lim_x_pos_0 = (_S92 - cx_7) / fx_7 + _S94;
    float _S95 = 0.30000001192092896f * (0.5f * _S93 / fy_7);
    float lim_y_pos_0 = (_S93 - cy_7) / fy_7 + _S95;
    float rz_3 = 1.0f / mean_c_3.z;
    float _S96 = mean_c_3.z * mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    float _S97 = - (cx_7 / fx_7 + _S94);
    float _S98 = mean_c_3.x * rz_3;
    float _S99 = s_primal_ctx_max_0(_S97, _S98);
    float _S100 = s_primal_ctx_min_0(lim_x_pos_0, _S99);
    float _S101 = - (cy_7 / fy_7 + _S95);
    float _S102 = mean_c_3.y * rz_3;
    float _S103 = s_primal_ctx_max_0(_S101, _S102);
    float _S104 = s_primal_ctx_min_0(lim_y_pos_0, _S103);
    float _S105 = - fx_7;
    float _S106 = _S105 * (mean_c_3.z * _S100);
    float _S107 = - fy_7;
    float _S108 = _S107 * (mean_c_3.z * _S104);
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, _S106 * rz2_3, 0.0f, fy_7 * rz_3, _S108 * rz2_3);
    Matrix<float, 2, 3>  _S109 = s_primal_ctx_mul_2(J_7, _S91);
    Matrix<float, 3, 2>  _S110 = transpose_1(J_7);
    Matrix<float, 2, 2>  _S111 = s_primal_ctx_mul_3(_S109, _S110);
    float _S112 = fx_7 * mean_c_3.x;
    float _S113 = fy_7 * mean_c_3.y;
    float _S114 = _S111.rows[int(0)].y * _S111.rows[int(1)].x;
    float det_orig_1 = _S111.rows[int(0)].x * _S111.rows[int(1)].y - _S114;
    float _S115 = _S111.rows[int(0)].x + eps2d_4;
    Matrix<float, 2, 2>  _S116 = _S111;
    *&(((&_S116)->rows + (int(0)))->x) = _S115;
    float _S117 = _S111.rows[int(1)].y + eps2d_4;
    *&(((&_S116)->rows + (int(1)))->y) = _S117;
    float det_blur_1 = _S115 * _S117 - _S114;
    float _S118 = det_orig_1 / det_blur_1;
    float _S119 = det_blur_1 * det_blur_1;
    float _S120 = s_primal_ctx_max_0(0.0f, _S118);
    float invdet_4 = 1.0f / det_blur_1;
    float _S121 = - _S111.rows[int(0)].y;
    float _S122 = - _S111.rows[int(1)].x;
    float _S123 = - in_opacity_3;
    float _S124 = 1.0f + s_primal_ctx_exp_0(_S123);
    float _S125 = _S124 * _S124;
    float _S126 = 1.0f / _S124 / 0.00392156885936856f;
    float _S127 = 2.0f * s_primal_ctx_log_0(_S126);
    float _S128 = s_primal_ctx_sqrt_0(_S127);
    float _S129 = _S116.rows[int(0)].x;
    float _S130 = _S116.rows[int(1)].y;
    float2  _S131 = make_float2 (0.0f);
    float2  _S132 = _S131;
    *&((&_S132)->y) = v_conic_0.z;
    float2  _S133 = _S131;
    *&((&_S133)->y) = v_conic_0.y;
    *&((&_S133)->x) = v_conic_0.x;
    DiffPair_float_0 _S134;
    (&_S134)->primal_0 = _S130;
    (&_S134)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S134, 0.0f);
    DiffPair_float_0 _S135;
    (&_S135)->primal_0 = _S129;
    (&_S135)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S135, 0.0f);
    DiffPair_float_0 _S136;
    (&_S136)->primal_0 = 3.32999992370605469f;
    (&_S136)->differential_0 = 0.0f;
    DiffPair_float_0 _S137;
    (&_S137)->primal_0 = _S128;
    (&_S137)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S136, &_S137, 0.0f);
    DiffPair_float_0 _S138;
    (&_S138)->primal_0 = _S127;
    (&_S138)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S138, _S137.differential_0);
    float _S139 = 2.0f * _S138.differential_0;
    DiffPair_float_0 _S140;
    (&_S140)->primal_0 = _S126;
    (&_S140)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S140, _S139);
    float _S141 = - ((v_opacity_0 + 254.9999847412109375f * _S140.differential_0) / _S125);
    DiffPair_float_0 _S142;
    (&_S142)->primal_0 = _S123;
    (&_S142)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S142, _S141);
    float _S143 = - _S142.differential_0;
    Matrix<float, 2, 2>  _S144 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S145 = _S144;
    _S145[int(1)] = _S132;
    _S145[int(0)] = _S133;
    float _S146 = invdet_4 * _S145.rows[int(1)].y;
    float _S147 = - (invdet_4 * _S145.rows[int(1)].x);
    float _S148 = - (invdet_4 * _S145.rows[int(0)].y);
    float _S149 = invdet_4 * _S145.rows[int(0)].x;
    float _S150 = - ((_S115 * _S145.rows[int(1)].y + _S122 * _S145.rows[int(1)].x + _S121 * _S145.rows[int(0)].y + _S117 * _S145.rows[int(0)].x) / _S119);
    DiffPair_float_0 _S151;
    (&_S151)->primal_0 = _S120;
    (&_S151)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S151, 0.0f);
    DiffPair_float_0 _S152;
    (&_S152)->primal_0 = 0.0f;
    (&_S152)->differential_0 = 0.0f;
    DiffPair_float_0 _S153;
    (&_S153)->primal_0 = _S118;
    (&_S153)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S152, &_S153, _S151.differential_0);
    float _S154 = _S153.differential_0 / _S119;
    float s_diff_det_orig_T_0 = det_blur_1 * _S154;
    float _S155 = _S150 + det_orig_1 * - _S154;
    float _S156 = - _S155;
    float _S157 = _S115 * _S155;
    float _S158 = _S117 * _S155;
    float2  _S159 = make_float2 (0.0f, _S134.differential_0);
    float2  _S160 = make_float2 (_S135.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S161 = _S144;
    _S161[int(1)] = _S159;
    _S161[int(0)] = _S160;
    _S116 = _S161;
    *&(((&_S116)->rows + (int(1)))->y) = 0.0f;
    float _S162 = _S149 + _S157 + _S161.rows[int(1)].y;
    *&(((&_S116)->rows + (int(0)))->x) = 0.0f;
    float _S163 = _S146 + _S158 + _S161.rows[int(0)].x;
    float _S164 = _S156 + - s_diff_det_orig_T_0;
    float _S165 = _S147 + _S111.rows[int(0)].y * _S164;
    float _S166 = _S148 + _S111.rows[int(1)].x * _S164;
    float _S167 = _S111.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S168 = _S162 + _S111.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S169 = _S131;
    *&((&_S169)->x) = _S165;
    *&((&_S169)->y) = _S168;
    float _S170 = _S163 + _S167;
    float2  _S171 = _S131;
    *&((&_S171)->y) = _S166;
    *&((&_S171)->x) = _S170;
    float _S172 = _S113 * v_mean2d_0.y;
    float _S173 = fy_7 * (rz_3 * v_mean2d_0.y);
    float _S174 = _S112 * v_mean2d_0.x;
    float _S175 = fx_7 * (rz_3 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S176 = _S144;
    _S176[int(1)] = _S169;
    _S176[int(0)] = _S171;
    Matrix<float, 2, 2>  _S177 = _S116 + _S176;
    Matrix<float, 2, 3>  _S178 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S179;
    (&_S179)->primal_0 = _S109;
    (&_S179)->differential_0 = _S178;
    Matrix<float, 3, 2>  _S180 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S181;
    (&_S181)->primal_0 = _S110;
    (&_S181)->differential_0 = _S180;
    s_bwd_prop_mul_0(&_S179, &_S181, _S177);
    Matrix<float, 2, 3>  _S182 = transpose_2(_S181.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S183;
    (&_S183)->primal_0 = J_7;
    (&_S183)->differential_0 = _S178;
    Matrix<float, 3, 3>  _S184 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S185;
    (&_S185)->primal_0 = _S91;
    (&_S185)->differential_0 = _S184;
    s_bwd_prop_mul_1(&_S183, &_S185, _S179.differential_0);
    Matrix<float, 2, 3>  _S186 = _S182 + _S183.differential_0;
    float _S187 = _S108 * _S186.rows[int(1)].z;
    float s_diff_ty_T_0 = _S107 * (rz2_3 * _S186.rows[int(1)].z);
    float _S188 = fy_7 * _S186.rows[int(1)].y;
    float _S189 = _S106 * _S186.rows[int(0)].z;
    float s_diff_tx_T_0 = _S105 * (rz2_3 * _S186.rows[int(0)].z);
    float _S190 = fx_7 * _S186.rows[int(0)].x;
    float _S191 = mean_c_3.z * s_diff_ty_T_0;
    float _S192 = _S104 * s_diff_ty_T_0;
    DiffPair_float_0 _S193;
    (&_S193)->primal_0 = lim_y_pos_0;
    (&_S193)->differential_0 = 0.0f;
    DiffPair_float_0 _S194;
    (&_S194)->primal_0 = _S103;
    (&_S194)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S193, &_S194, _S191);
    DiffPair_float_0 _S195;
    (&_S195)->primal_0 = _S101;
    (&_S195)->differential_0 = 0.0f;
    DiffPair_float_0 _S196;
    (&_S196)->primal_0 = _S102;
    (&_S196)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S195, &_S196, _S194.differential_0);
    float _S197 = mean_c_3.y * _S196.differential_0;
    float _S198 = rz_3 * _S196.differential_0;
    float _S199 = mean_c_3.z * s_diff_tx_T_0;
    float _S200 = _S100 * s_diff_tx_T_0;
    DiffPair_float_0 _S201;
    (&_S201)->primal_0 = lim_x_pos_0;
    (&_S201)->differential_0 = 0.0f;
    DiffPair_float_0 _S202;
    (&_S202)->primal_0 = _S99;
    (&_S202)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S201, &_S202, _S199);
    DiffPair_float_0 _S203;
    (&_S203)->primal_0 = _S97;
    (&_S203)->differential_0 = 0.0f;
    DiffPair_float_0 _S204;
    (&_S204)->primal_0 = _S98;
    (&_S204)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S203, &_S204, _S202.differential_0);
    float _S205 = rz_3 * (_S187 + _S189);
    float _S206 = _S192 + _S200 + - ((_S172 + _S174 + _S188 + _S190 + _S197 + mean_c_3.x * _S204.differential_0 + _S205 + _S205) / _S96);
    float _S207 = _S173 + _S198;
    float _S208 = _S175 + rz_3 * _S204.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S209;
    (&_S209)->primal_0 = _S89;
    (&_S209)->differential_0 = _S184;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S210;
    (&_S210)->primal_0 = _S90;
    (&_S210)->differential_0 = _S184;
    s_bwd_prop_mul_2(&_S209, &_S210, _S185.differential_0);
    Matrix<float, 3, 3>  _S211 = transpose_0(_S210.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S212;
    (&_S212)->primal_0 = R_5;
    (&_S212)->differential_0 = _S184;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S213;
    (&_S213)->primal_0 = _S88;
    (&_S213)->differential_0 = _S184;
    s_bwd_prop_mul_2(&_S212, &_S213, _S209.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S214;
    (&_S214)->primal_0 = _S86;
    (&_S214)->differential_0 = _S184;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S215;
    (&_S215)->primal_0 = _S87;
    (&_S215)->differential_0 = _S184;
    s_bwd_prop_mul_2(&_S214, &_S215, _S213.differential_0);
    Matrix<float, 3, 3>  _S216 = _S214.differential_0 + transpose_0(_S215.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S217;
    (&_S217)->primal_0 = _S85;
    (&_S217)->differential_0 = _S184;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S218;
    (&_S218)->primal_0 = S_0;
    (&_S218)->differential_0 = _S184;
    s_bwd_prop_mul_2(&_S217, &_S218, _S216);
    Matrix<float, 3, 3>  _S219 = transpose_0(_S217.differential_0);
    float _S220 = 2.0f * - _S219.rows[int(2)].z;
    float _S221 = 2.0f * _S219.rows[int(2)].y;
    float _S222 = 2.0f * _S219.rows[int(2)].x;
    float _S223 = 2.0f * _S219.rows[int(1)].z;
    float _S224 = 2.0f * - _S219.rows[int(1)].y;
    float _S225 = 2.0f * _S219.rows[int(1)].x;
    float _S226 = 2.0f * _S219.rows[int(0)].z;
    float _S227 = 2.0f * _S219.rows[int(0)].y;
    float _S228 = 2.0f * - _S219.rows[int(0)].x;
    float _S229 = - _S225 + _S227;
    float _S230 = _S222 + - _S226;
    float _S231 = - _S221 + _S223;
    float _S232 = _S221 + _S223;
    float _S233 = _S222 + _S226;
    float _S234 = _S225 + _S227;
    float _S235 = z_5 * (_S224 + _S228);
    float _S236 = y_6 * (_S220 + _S228);
    float _S237 = x_15 * (_S220 + _S224);
    float _S238 = z_5 * _S229 + y_6 * _S230 + x_15 * _S231;
    float _S239 = _S84 * _S238;
    float _S240 = w_5 * _S229 + y_6 * _S232 + x_15 * _S233 + _S235 + _S235;
    float _S241 = _S84 * _S240;
    float _S242 = w_5 * _S230 + z_5 * _S232 + x_15 * _S234 + _S236 + _S236;
    float _S243 = _S84 * _S242;
    float _S244 = w_5 * _S231 + z_5 * _S233 + y_6 * _S234 + _S237 + _S237;
    float _S245 = _S84 * _S244;
    float _S246 = quat_5.x * _S238 + quat_5.w * _S240 + quat_5.z * _S242 + quat_5.y * _S244;
    DiffPair_float_0 _S247;
    (&_S247)->primal_0 = _S83;
    (&_S247)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S247, _S246);
    float _S248 = quat_5.x * _S247.differential_0;
    float _S249 = quat_5.w * _S247.differential_0;
    float _S250 = quat_5.z * _S247.differential_0;
    float _S251 = quat_5.y * _S247.differential_0;
    float _S252 = _S241 + _S249 + _S249;
    float _S253 = _S243 + _S250 + _S250;
    float _S254 = _S245 + _S251 + _S251;
    float _S255 = _S239 + _S248 + _S248;
    float3  _S256 = make_float3 (0.0f, 0.0f, v_depth_0);
    float3  _S257 = make_float3 (0.0f);
    float3  _S258 = _S257;
    *&((&_S258)->z) = _S206;
    *&((&_S258)->y) = _S207;
    *&((&_S258)->x) = _S208;
    float3  _S259 = _S256 + _S258;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S260;
    (&_S260)->primal_0 = R_5;
    (&_S260)->differential_0 = _S184;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S261;
    (&_S261)->primal_0 = mean_3;
    (&_S261)->differential_0 = _S257;
    s_bwd_prop_mul_3(&_S260, &_S261, _S259);
    Matrix<float, 3, 3>  _S262 = _S211 + _S212.differential_0 + _S260.differential_0;
    float3  _S263 = _S257;
    *&((&_S263)->z) = _S218.differential_0.rows[int(2)].z;
    *&((&_S263)->y) = _S218.differential_0.rows[int(1)].y;
    *&((&_S263)->x) = _S218.differential_0.rows[int(0)].x;
    float4  _S264 = make_float4 (0.0f);
    *&((&_S264)->w) = _S252;
    *&((&_S264)->z) = _S253;
    *&((&_S264)->y) = _S254;
    *&((&_S264)->x) = _S255;
    *v_mean_0 = _S261.differential_0;
    *v_quat_0 = _S264;
    *v_scale_0 = _S263;
    *v_in_opacity_0 = _S143;
    *v_R_0 = _S262;
    *v_t_0 = _S259;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S265, float _S266)
{
    return (F32_atan2((_S265), (_S266)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S267, DiffPair_float_0 * _S268, float _S269)
{
    _d_atan2_0(_S267, _S268, _S269);
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_7, float _s_dOut_0)
{
    float _S270 = (*dpx_7).primal_0.x;
    float _S271 = (*dpx_7).primal_0.y;
    DiffPair_float_0 _S272;
    (&_S272)->primal_0 = _S270 * _S270 + _S271 * _S271;
    (&_S272)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S272, _s_dOut_0);
    float _S273 = (*dpx_7).primal_0.y * _S272.differential_0;
    float _S274 = _S273 + _S273;
    float _S275 = (*dpx_7).primal_0.x * _S272.differential_0;
    float _S276 = _S275 + _S275;
    float2  _S277 = make_float2 (0.0f);
    *&((&_S277)->y) = _S274;
    *&((&_S277)->x) = _S276;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S277;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S278, float _S279)
{
    s_bwd_prop_length_impl_0(_S278, _S279);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(float3  mean_4, float4  quat_6, float3  scale_5, float in_opacity_4, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_8, float fy_8, float cx_8, float cy_8, uint image_width_4, uint image_height_4, float eps2d_5, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_4 = s_primal_ctx_mul_0(R_6, mean_4) + t_5;
    float _S280 = quat_6.y;
    float _S281 = _S280 * _S280 + quat_6.z * quat_6.z + quat_6.w * quat_6.w + quat_6.x * quat_6.x;
    float _S282 = s_primal_ctx_rsqrt_0(_S281);
    float x_16 = quat_6.y * _S282;
    float y_7 = quat_6.z * _S282;
    float z_6 = quat_6.w * _S282;
    float w_6 = quat_6.x * _S282;
    float x2_8 = x_16 * x_16;
    float y2_8 = y_7 * y_7;
    float z2_6 = z_6 * z_6;
    float xy_8 = x_16 * y_7;
    float xz_6 = x_16 * z_6;
    float yz_6 = y_7 * z_6;
    float wx_6 = w_6 * x_16;
    float wy_6 = w_6 * y_7;
    float wz_6 = w_6 * z_6;
    Matrix<float, 3, 3>  _S283 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_6), 2.0f * (xy_8 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_8 - wz_6), 1.0f - 2.0f * (x2_8 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_8 + y2_8)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (scale_5.x, 0.0f, 0.0f, 0.0f, scale_5.y, 0.0f, 0.0f, 0.0f, scale_5.z);
    Matrix<float, 3, 3>  _S284 = s_primal_ctx_mul_1(_S283, S_1);
    Matrix<float, 3, 3>  _S285 = transpose_0(_S284);
    Matrix<float, 3, 3>  _S286 = s_primal_ctx_mul_1(_S284, _S285);
    Matrix<float, 3, 3>  _S287 = s_primal_ctx_mul_1(R_6, _S286);
    Matrix<float, 3, 3>  _S288 = transpose_0(R_6);
    Matrix<float, 3, 3>  _S289 = s_primal_ctx_mul_1(_S287, _S288);
    float2  _S290 = make_float2 (mean_c_4.x, mean_c_4.y);
    float xy_len_2 = length_0(_S290) + 1.00000001168609742e-07f;
    float _S291 = mean_c_4.z + 1.00000001168609742e-07f;
    float _S292 = s_primal_ctx_atan2_0(xy_len_2, _S291);
    float _S293 = mean_c_4.x * fx_8;
    float _S294 = _S293 * _S292;
    float _S295 = xy_len_2 * xy_len_2;
    float _S296 = mean_c_4.y * fy_8;
    float _S297 = _S296 * _S292;
    float x2_9 = mean_c_4.x * mean_c_4.x + 1.00000001168609742e-07f;
    float y2_9 = mean_c_4.y * mean_c_4.y;
    float xy_9 = mean_c_4.x * mean_c_4.y;
    float x2y2_2 = x2_9 + y2_9;
    float _S298 = x2y2_2 + mean_c_4.z * mean_c_4.z;
    float x2y2z2_inv_2 = 1.0f / _S298;
    float _S299 = _S298 * _S298;
    float _S300 = s_primal_ctx_atan2_0(xy_len_2, mean_c_4.z);
    float _S301 = _S300 / xy_len_2;
    float b_2 = _S301 / x2y2_2;
    float _S302 = x2y2_2 * x2y2_2;
    float _S303 = mean_c_4.z * x2y2z2_inv_2;
    float a_2 = _S303 / x2y2_2;
    float _S304 = fx_8 * xy_9;
    float _S305 = a_2 - b_2;
    float _S306 = - fx_8;
    float _S307 = _S306 * mean_c_4.x;
    float _S308 = fy_8 * xy_9;
    float _S309 = - fy_8;
    float _S310 = _S309 * mean_c_4.y;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_8 * (x2_9 * a_2 + y2_9 * b_2), _S304 * _S305, _S307 * x2y2z2_inv_2, _S308 * _S305, fy_8 * (y2_9 * a_2 + x2_9 * b_2), _S310 * x2y2z2_inv_2);
    Matrix<float, 2, 3>  _S311 = s_primal_ctx_mul_2(J_8, _S289);
    Matrix<float, 3, 2>  _S312 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S313 = s_primal_ctx_mul_3(_S311, _S312);
    float _S314 = _S313.rows[int(0)].y * _S313.rows[int(1)].x;
    float det_orig_2 = _S313.rows[int(0)].x * _S313.rows[int(1)].y - _S314;
    float _S315 = _S313.rows[int(0)].x + eps2d_5;
    Matrix<float, 2, 2>  _S316 = _S313;
    *&(((&_S316)->rows + (int(0)))->x) = _S315;
    float _S317 = _S313.rows[int(1)].y + eps2d_5;
    *&(((&_S316)->rows + (int(1)))->y) = _S317;
    float det_blur_2 = _S315 * _S317 - _S314;
    float _S318 = det_orig_2 / det_blur_2;
    float _S319 = det_blur_2 * det_blur_2;
    float _S320 = s_primal_ctx_max_0(0.0f, _S318);
    float invdet_5 = 1.0f / det_blur_2;
    float _S321 = - _S313.rows[int(0)].y;
    float _S322 = - _S313.rows[int(1)].x;
    float _S323 = - in_opacity_4;
    float _S324 = 1.0f + s_primal_ctx_exp_0(_S323);
    float _S325 = _S324 * _S324;
    float _S326 = 1.0f / _S324 / 0.00392156885936856f;
    float _S327 = 2.0f * s_primal_ctx_log_0(_S326);
    float _S328 = s_primal_ctx_sqrt_0(_S327);
    float _S329 = _S316.rows[int(0)].x;
    float _S330 = _S316.rows[int(1)].y;
    float2  _S331 = make_float2 (0.0f);
    float2  _S332 = _S331;
    *&((&_S332)->y) = v_conic_1.z;
    float2  _S333 = _S331;
    *&((&_S333)->y) = v_conic_1.y;
    *&((&_S333)->x) = v_conic_1.x;
    DiffPair_float_0 _S334;
    (&_S334)->primal_0 = _S330;
    (&_S334)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S334, 0.0f);
    DiffPair_float_0 _S335;
    (&_S335)->primal_0 = _S329;
    (&_S335)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S335, 0.0f);
    DiffPair_float_0 _S336;
    (&_S336)->primal_0 = 3.32999992370605469f;
    (&_S336)->differential_0 = 0.0f;
    DiffPair_float_0 _S337;
    (&_S337)->primal_0 = _S328;
    (&_S337)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S336, &_S337, 0.0f);
    DiffPair_float_0 _S338;
    (&_S338)->primal_0 = _S327;
    (&_S338)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S338, _S337.differential_0);
    float _S339 = 2.0f * _S338.differential_0;
    DiffPair_float_0 _S340;
    (&_S340)->primal_0 = _S326;
    (&_S340)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S340, _S339);
    float _S341 = - ((v_opacity_1 + 254.9999847412109375f * _S340.differential_0) / _S325);
    DiffPair_float_0 _S342;
    (&_S342)->primal_0 = _S323;
    (&_S342)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S342, _S341);
    float _S343 = - _S342.differential_0;
    Matrix<float, 2, 2>  _S344 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S345 = _S344;
    _S345[int(1)] = _S332;
    _S345[int(0)] = _S333;
    float _S346 = invdet_5 * _S345.rows[int(1)].y;
    float _S347 = - (invdet_5 * _S345.rows[int(1)].x);
    float _S348 = - (invdet_5 * _S345.rows[int(0)].y);
    float _S349 = invdet_5 * _S345.rows[int(0)].x;
    float _S350 = - ((_S315 * _S345.rows[int(1)].y + _S322 * _S345.rows[int(1)].x + _S321 * _S345.rows[int(0)].y + _S317 * _S345.rows[int(0)].x) / _S319);
    DiffPair_float_0 _S351;
    (&_S351)->primal_0 = _S320;
    (&_S351)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S351, 0.0f);
    DiffPair_float_0 _S352;
    (&_S352)->primal_0 = 0.0f;
    (&_S352)->differential_0 = 0.0f;
    DiffPair_float_0 _S353;
    (&_S353)->primal_0 = _S318;
    (&_S353)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S352, &_S353, _S351.differential_0);
    float _S354 = _S353.differential_0 / _S319;
    float s_diff_det_orig_T_1 = det_blur_2 * _S354;
    float _S355 = _S350 + det_orig_2 * - _S354;
    float _S356 = - _S355;
    float _S357 = _S315 * _S355;
    float _S358 = _S317 * _S355;
    float2  _S359 = make_float2 (0.0f, _S334.differential_0);
    float2  _S360 = make_float2 (_S335.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S361 = _S344;
    _S361[int(1)] = _S359;
    _S361[int(0)] = _S360;
    _S316 = _S361;
    *&(((&_S316)->rows + (int(1)))->y) = 0.0f;
    float _S362 = _S349 + _S357 + _S361.rows[int(1)].y;
    *&(((&_S316)->rows + (int(0)))->x) = 0.0f;
    float _S363 = _S346 + _S358 + _S361.rows[int(0)].x;
    float _S364 = _S356 + - s_diff_det_orig_T_1;
    float _S365 = _S347 + _S313.rows[int(0)].y * _S364;
    float _S366 = _S348 + _S313.rows[int(1)].x * _S364;
    float _S367 = _S313.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S368 = _S362 + _S313.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S369 = _S331;
    *&((&_S369)->x) = _S365;
    *&((&_S369)->y) = _S368;
    float _S370 = _S363 + _S367;
    float2  _S371 = _S331;
    *&((&_S371)->y) = _S366;
    *&((&_S371)->x) = _S370;
    Matrix<float, 2, 2>  _S372 = _S344;
    _S372[int(1)] = _S369;
    _S372[int(0)] = _S371;
    Matrix<float, 2, 2>  _S373 = _S316 + _S372;
    Matrix<float, 2, 3>  _S374 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S375;
    (&_S375)->primal_0 = _S311;
    (&_S375)->differential_0 = _S374;
    Matrix<float, 3, 2>  _S376 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S377;
    (&_S377)->primal_0 = _S312;
    (&_S377)->differential_0 = _S376;
    s_bwd_prop_mul_0(&_S375, &_S377, _S373);
    Matrix<float, 2, 3>  _S378 = transpose_2(_S377.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S379;
    (&_S379)->primal_0 = J_8;
    (&_S379)->differential_0 = _S374;
    Matrix<float, 3, 3>  _S380 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S381;
    (&_S381)->primal_0 = _S289;
    (&_S381)->differential_0 = _S380;
    s_bwd_prop_mul_1(&_S379, &_S381, _S375.differential_0);
    Matrix<float, 2, 3>  _S382 = _S378 + _S379.differential_0;
    float _S383 = _S310 * _S382.rows[int(1)].z;
    float _S384 = _S309 * (x2y2z2_inv_2 * _S382.rows[int(1)].z);
    float _S385 = fy_8 * _S382.rows[int(1)].y;
    float _S386 = b_2 * _S385;
    float _S387 = a_2 * _S385;
    float _S388 = fy_8 * (_S305 * _S382.rows[int(1)].x);
    float _S389 = _S307 * _S382.rows[int(0)].z;
    float _S390 = _S306 * (x2y2z2_inv_2 * _S382.rows[int(0)].z);
    float _S391 = _S308 * _S382.rows[int(1)].x + _S304 * _S382.rows[int(0)].y;
    float _S392 = fx_8 * (_S305 * _S382.rows[int(0)].y);
    float _S393 = fx_8 * _S382.rows[int(0)].x;
    float _S394 = b_2 * _S393;
    float _S395 = a_2 * _S393;
    float _S396 = (y2_9 * _S385 + _S391 + x2_9 * _S393) / _S302;
    float _S397 = _S303 * - _S396;
    float _S398 = x2y2_2 * _S396;
    float _S399 = mean_c_4.z * _S398;
    float _S400 = x2y2z2_inv_2 * _S398;
    float _S401 = (x2_9 * _S385 + - _S391 + y2_9 * _S393) / _S302;
    float _S402 = _S301 * - _S401;
    float _S403 = x2y2_2 * _S401 / _S295;
    float _S404 = _S300 * - _S403;
    float _S405 = xy_len_2 * _S403;
    DiffPair_float_0 _S406;
    (&_S406)->primal_0 = xy_len_2;
    (&_S406)->differential_0 = 0.0f;
    DiffPair_float_0 _S407;
    (&_S407)->primal_0 = mean_c_4.z;
    (&_S407)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S406, &_S407, _S405);
    float _S408 = - ((_S383 + _S389 + _S399) / _S299);
    float _S409 = mean_c_4.z * _S408;
    float _S410 = _S397 + _S402 + _S408;
    float _S411 = _S388 + _S392;
    float _S412 = mean_c_4.x * _S411;
    float _S413 = mean_c_4.y * _S411;
    float _S414 = mean_c_4.y * (_S387 + _S394 + _S410);
    float _S415 = mean_c_4.x * (_S386 + _S395 + _S410);
    float _S416 = v_mean2d_1.y / _S295;
    float _S417 = _S297 * - _S416;
    float _S418 = xy_len_2 * _S416;
    float _S419 = fy_8 * (_S292 * _S418);
    float _S420 = v_mean2d_1.x / _S295;
    float _S421 = _S294 * - _S420;
    float _S422 = xy_len_2 * _S420;
    float _S423 = fx_8 * (_S292 * _S422);
    float _S424 = _S296 * _S418 + _S293 * _S422;
    DiffPair_float_0 _S425;
    (&_S425)->primal_0 = xy_len_2;
    (&_S425)->differential_0 = 0.0f;
    DiffPair_float_0 _S426;
    (&_S426)->primal_0 = _S291;
    (&_S426)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S425, &_S426, _S424);
    float _S427 = _S404 + _S406.differential_0 + _S417 + _S421 + _S425.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S428;
    (&_S428)->primal_0 = _S290;
    (&_S428)->differential_0 = _S331;
    s_bwd_length_impl_0(&_S428, _S427);
    float _S429 = _S400 + _S407.differential_0 + _S409 + _S409 + _S426.differential_0;
    float _S430 = _S384 + _S412 + _S414 + _S414 + _S419 + _S428.differential_0.y;
    float _S431 = _S390 + _S413 + _S415 + _S415 + _S423 + _S428.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S432;
    (&_S432)->primal_0 = _S287;
    (&_S432)->differential_0 = _S380;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S433;
    (&_S433)->primal_0 = _S288;
    (&_S433)->differential_0 = _S380;
    s_bwd_prop_mul_2(&_S432, &_S433, _S381.differential_0);
    Matrix<float, 3, 3>  _S434 = transpose_0(_S433.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S435;
    (&_S435)->primal_0 = R_6;
    (&_S435)->differential_0 = _S380;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S436;
    (&_S436)->primal_0 = _S286;
    (&_S436)->differential_0 = _S380;
    s_bwd_prop_mul_2(&_S435, &_S436, _S432.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S437;
    (&_S437)->primal_0 = _S284;
    (&_S437)->differential_0 = _S380;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S438;
    (&_S438)->primal_0 = _S285;
    (&_S438)->differential_0 = _S380;
    s_bwd_prop_mul_2(&_S437, &_S438, _S436.differential_0);
    Matrix<float, 3, 3>  _S439 = _S437.differential_0 + transpose_0(_S438.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S440;
    (&_S440)->primal_0 = _S283;
    (&_S440)->differential_0 = _S380;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S441;
    (&_S441)->primal_0 = S_1;
    (&_S441)->differential_0 = _S380;
    s_bwd_prop_mul_2(&_S440, &_S441, _S439);
    Matrix<float, 3, 3>  _S442 = transpose_0(_S440.differential_0);
    float _S443 = 2.0f * - _S442.rows[int(2)].z;
    float _S444 = 2.0f * _S442.rows[int(2)].y;
    float _S445 = 2.0f * _S442.rows[int(2)].x;
    float _S446 = 2.0f * _S442.rows[int(1)].z;
    float _S447 = 2.0f * - _S442.rows[int(1)].y;
    float _S448 = 2.0f * _S442.rows[int(1)].x;
    float _S449 = 2.0f * _S442.rows[int(0)].z;
    float _S450 = 2.0f * _S442.rows[int(0)].y;
    float _S451 = 2.0f * - _S442.rows[int(0)].x;
    float _S452 = - _S448 + _S450;
    float _S453 = _S445 + - _S449;
    float _S454 = - _S444 + _S446;
    float _S455 = _S444 + _S446;
    float _S456 = _S445 + _S449;
    float _S457 = _S448 + _S450;
    float _S458 = z_6 * (_S447 + _S451);
    float _S459 = y_7 * (_S443 + _S451);
    float _S460 = x_16 * (_S443 + _S447);
    float _S461 = z_6 * _S452 + y_7 * _S453 + x_16 * _S454;
    float _S462 = _S282 * _S461;
    float _S463 = w_6 * _S452 + y_7 * _S455 + x_16 * _S456 + _S458 + _S458;
    float _S464 = _S282 * _S463;
    float _S465 = w_6 * _S453 + z_6 * _S455 + x_16 * _S457 + _S459 + _S459;
    float _S466 = _S282 * _S465;
    float _S467 = w_6 * _S454 + z_6 * _S456 + y_7 * _S457 + _S460 + _S460;
    float _S468 = _S282 * _S467;
    float _S469 = quat_6.x * _S461 + quat_6.w * _S463 + quat_6.z * _S465 + quat_6.y * _S467;
    DiffPair_float_0 _S470;
    (&_S470)->primal_0 = _S281;
    (&_S470)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S470, _S469);
    float _S471 = quat_6.x * _S470.differential_0;
    float _S472 = quat_6.w * _S470.differential_0;
    float _S473 = quat_6.z * _S470.differential_0;
    float _S474 = quat_6.y * _S470.differential_0;
    float _S475 = _S464 + _S472 + _S472;
    float _S476 = _S466 + _S473 + _S473;
    float _S477 = _S468 + _S474 + _S474;
    float _S478 = _S462 + _S471 + _S471;
    float3  _S479 = make_float3 (0.0f, 0.0f, v_depth_1);
    float3  _S480 = make_float3 (0.0f);
    float3  _S481 = _S480;
    *&((&_S481)->z) = _S429;
    *&((&_S481)->y) = _S430;
    *&((&_S481)->x) = _S431;
    float3  _S482 = _S479 + _S481;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S483;
    (&_S483)->primal_0 = R_6;
    (&_S483)->differential_0 = _S380;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S484;
    (&_S484)->primal_0 = mean_4;
    (&_S484)->differential_0 = _S480;
    s_bwd_prop_mul_3(&_S483, &_S484, _S482);
    Matrix<float, 3, 3>  _S485 = _S434 + _S435.differential_0 + _S483.differential_0;
    float3  _S486 = _S480;
    *&((&_S486)->z) = _S441.differential_0.rows[int(2)].z;
    *&((&_S486)->y) = _S441.differential_0.rows[int(1)].y;
    *&((&_S486)->x) = _S441.differential_0.rows[int(0)].x;
    float4  _S487 = make_float4 (0.0f);
    *&((&_S487)->w) = _S475;
    *&((&_S487)->z) = _S476;
    *&((&_S487)->y) = _S477;
    *&((&_S487)->x) = _S478;
    *v_mean_1 = _S484.differential_0;
    *v_quat_1 = _S487;
    *v_scale_1 = _S486;
    *v_in_opacity_1 = _S343;
    *v_R_1 = _S485;
    *v_t_1 = _S482;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(float3  mean_5, float4  quat_7, float3  scale_6, float in_opacity_5, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_9, float fy_9, float cx_9, float cy_9, uint image_width_5, uint image_height_5, float eps2d_6, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float _S488 = quat_7.y;
    float _S489 = _S488 * _S488 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x;
    float _S490 = s_primal_ctx_rsqrt_0(_S489);
    float x_17 = quat_7.y * _S490;
    float y_8 = quat_7.z * _S490;
    float z_7 = quat_7.w * _S490;
    float w_7 = quat_7.x * _S490;
    float x2_10 = x_17 * x_17;
    float y2_10 = y_8 * y_8;
    float z2_7 = z_7 * z_7;
    float xy_10 = x_17 * y_8;
    float xz_7 = x_17 * z_7;
    float yz_7 = y_8 * z_7;
    float wx_7 = w_7 * x_17;
    float wy_7 = w_7 * y_8;
    float wz_7 = w_7 * z_7;
    Matrix<float, 3, 3>  _S491 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_7), 2.0f * (xy_10 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_10 - wz_7), 1.0f - 2.0f * (x2_10 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (scale_6.x, 0.0f, 0.0f, 0.0f, scale_6.y, 0.0f, 0.0f, 0.0f, scale_6.z);
    Matrix<float, 3, 3>  _S492 = s_primal_ctx_mul_1(_S491, S_2);
    Matrix<float, 3, 3>  _S493 = transpose_0(_S492);
    Matrix<float, 3, 3>  _S494 = s_primal_ctx_mul_1(_S492, _S493);
    Matrix<float, 3, 3>  _S495 = s_primal_ctx_mul_1(R_7, _S494);
    Matrix<float, 3, 3>  _S496 = transpose_0(R_7);
    Matrix<float, 3, 3>  _S497 = s_primal_ctx_mul_1(_S495, _S496);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 3>  _S498 = s_primal_ctx_mul_2(J_9, _S497);
    Matrix<float, 3, 2>  _S499 = transpose_1(J_9);
    Matrix<float, 2, 2>  _S500 = s_primal_ctx_mul_3(_S498, _S499);
    float _S501 = _S500.rows[int(0)].y * _S500.rows[int(1)].x;
    float det_orig_3 = _S500.rows[int(0)].x * _S500.rows[int(1)].y - _S501;
    float _S502 = _S500.rows[int(0)].x + eps2d_6;
    Matrix<float, 2, 2>  _S503 = _S500;
    *&(((&_S503)->rows + (int(0)))->x) = _S502;
    float _S504 = _S500.rows[int(1)].y + eps2d_6;
    *&(((&_S503)->rows + (int(1)))->y) = _S504;
    float det_blur_3 = _S502 * _S504 - _S501;
    float _S505 = det_orig_3 / det_blur_3;
    float _S506 = det_blur_3 * det_blur_3;
    float _S507 = s_primal_ctx_max_0(0.0f, _S505);
    float invdet_6 = 1.0f / det_blur_3;
    float _S508 = - _S500.rows[int(0)].y;
    float _S509 = - _S500.rows[int(1)].x;
    float _S510 = - in_opacity_5;
    float _S511 = 1.0f + s_primal_ctx_exp_0(_S510);
    float _S512 = _S511 * _S511;
    float _S513 = 1.0f / _S511 / 0.00392156885936856f;
    float _S514 = 2.0f * s_primal_ctx_log_0(_S513);
    float _S515 = s_primal_ctx_sqrt_0(_S514);
    float _S516 = _S503.rows[int(0)].x;
    float _S517 = _S503.rows[int(1)].y;
    float2  _S518 = make_float2 (0.0f);
    float2  _S519 = _S518;
    *&((&_S519)->y) = v_conic_2.z;
    float2  _S520 = _S518;
    *&((&_S520)->y) = v_conic_2.y;
    *&((&_S520)->x) = v_conic_2.x;
    DiffPair_float_0 _S521;
    (&_S521)->primal_0 = _S517;
    (&_S521)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S521, 0.0f);
    DiffPair_float_0 _S522;
    (&_S522)->primal_0 = _S516;
    (&_S522)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S522, 0.0f);
    DiffPair_float_0 _S523;
    (&_S523)->primal_0 = 3.32999992370605469f;
    (&_S523)->differential_0 = 0.0f;
    DiffPair_float_0 _S524;
    (&_S524)->primal_0 = _S515;
    (&_S524)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S523, &_S524, 0.0f);
    DiffPair_float_0 _S525;
    (&_S525)->primal_0 = _S514;
    (&_S525)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S525, _S524.differential_0);
    float _S526 = 2.0f * _S525.differential_0;
    DiffPair_float_0 _S527;
    (&_S527)->primal_0 = _S513;
    (&_S527)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S527, _S526);
    float _S528 = - ((v_opacity_2 + 254.9999847412109375f * _S527.differential_0) / _S512);
    DiffPair_float_0 _S529;
    (&_S529)->primal_0 = _S510;
    (&_S529)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S529, _S528);
    float _S530 = - _S529.differential_0;
    Matrix<float, 2, 2>  _S531 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S532 = _S531;
    _S532[int(1)] = _S519;
    _S532[int(0)] = _S520;
    float _S533 = invdet_6 * _S532.rows[int(1)].y;
    float _S534 = - (invdet_6 * _S532.rows[int(1)].x);
    float _S535 = - (invdet_6 * _S532.rows[int(0)].y);
    float _S536 = invdet_6 * _S532.rows[int(0)].x;
    float _S537 = - ((_S502 * _S532.rows[int(1)].y + _S509 * _S532.rows[int(1)].x + _S508 * _S532.rows[int(0)].y + _S504 * _S532.rows[int(0)].x) / _S506);
    DiffPair_float_0 _S538;
    (&_S538)->primal_0 = _S507;
    (&_S538)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S538, 0.0f);
    DiffPair_float_0 _S539;
    (&_S539)->primal_0 = 0.0f;
    (&_S539)->differential_0 = 0.0f;
    DiffPair_float_0 _S540;
    (&_S540)->primal_0 = _S505;
    (&_S540)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S539, &_S540, _S538.differential_0);
    float _S541 = _S540.differential_0 / _S506;
    float s_diff_det_orig_T_2 = det_blur_3 * _S541;
    float _S542 = _S537 + det_orig_3 * - _S541;
    float _S543 = - _S542;
    float _S544 = _S502 * _S542;
    float _S545 = _S504 * _S542;
    float2  _S546 = make_float2 (0.0f, _S521.differential_0);
    float2  _S547 = make_float2 (_S522.differential_0, 0.0f);
    Matrix<float, 2, 2>  _S548 = _S531;
    _S548[int(1)] = _S546;
    _S548[int(0)] = _S547;
    _S503 = _S548;
    *&(((&_S503)->rows + (int(1)))->y) = 0.0f;
    float _S549 = _S536 + _S544 + _S548.rows[int(1)].y;
    *&(((&_S503)->rows + (int(0)))->x) = 0.0f;
    float _S550 = _S533 + _S545 + _S548.rows[int(0)].x;
    float _S551 = _S543 + - s_diff_det_orig_T_2;
    float _S552 = _S534 + _S500.rows[int(0)].y * _S551;
    float _S553 = _S535 + _S500.rows[int(1)].x * _S551;
    float _S554 = _S500.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S555 = _S549 + _S500.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S556 = _S518;
    *&((&_S556)->x) = _S552;
    *&((&_S556)->y) = _S555;
    float _S557 = _S550 + _S554;
    float2  _S558 = _S518;
    *&((&_S558)->y) = _S553;
    *&((&_S558)->x) = _S557;
    float _S559 = fy_9 * v_mean2d_2.y;
    float _S560 = fx_9 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S561 = _S531;
    _S561[int(1)] = _S556;
    _S561[int(0)] = _S558;
    Matrix<float, 2, 2>  _S562 = _S503 + _S561;
    Matrix<float, 2, 3>  _S563 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S564;
    (&_S564)->primal_0 = _S498;
    (&_S564)->differential_0 = _S563;
    Matrix<float, 3, 2>  _S565 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S566;
    (&_S566)->primal_0 = _S499;
    (&_S566)->differential_0 = _S565;
    s_bwd_prop_mul_0(&_S564, &_S566, _S562);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S567;
    (&_S567)->primal_0 = J_9;
    (&_S567)->differential_0 = _S563;
    Matrix<float, 3, 3>  _S568 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S569;
    (&_S569)->primal_0 = _S497;
    (&_S569)->differential_0 = _S568;
    s_bwd_prop_mul_1(&_S567, &_S569, _S564.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S570;
    (&_S570)->primal_0 = _S495;
    (&_S570)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S571;
    (&_S571)->primal_0 = _S496;
    (&_S571)->differential_0 = _S568;
    s_bwd_prop_mul_2(&_S570, &_S571, _S569.differential_0);
    Matrix<float, 3, 3>  _S572 = transpose_0(_S571.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S573;
    (&_S573)->primal_0 = R_7;
    (&_S573)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S574;
    (&_S574)->primal_0 = _S494;
    (&_S574)->differential_0 = _S568;
    s_bwd_prop_mul_2(&_S573, &_S574, _S570.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S575;
    (&_S575)->primal_0 = _S492;
    (&_S575)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S576;
    (&_S576)->primal_0 = _S493;
    (&_S576)->differential_0 = _S568;
    s_bwd_prop_mul_2(&_S575, &_S576, _S574.differential_0);
    Matrix<float, 3, 3>  _S577 = _S575.differential_0 + transpose_0(_S576.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S578;
    (&_S578)->primal_0 = _S491;
    (&_S578)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S579;
    (&_S579)->primal_0 = S_2;
    (&_S579)->differential_0 = _S568;
    s_bwd_prop_mul_2(&_S578, &_S579, _S577);
    Matrix<float, 3, 3>  _S580 = transpose_0(_S578.differential_0);
    float _S581 = 2.0f * - _S580.rows[int(2)].z;
    float _S582 = 2.0f * _S580.rows[int(2)].y;
    float _S583 = 2.0f * _S580.rows[int(2)].x;
    float _S584 = 2.0f * _S580.rows[int(1)].z;
    float _S585 = 2.0f * - _S580.rows[int(1)].y;
    float _S586 = 2.0f * _S580.rows[int(1)].x;
    float _S587 = 2.0f * _S580.rows[int(0)].z;
    float _S588 = 2.0f * _S580.rows[int(0)].y;
    float _S589 = 2.0f * - _S580.rows[int(0)].x;
    float _S590 = - _S586 + _S588;
    float _S591 = _S583 + - _S587;
    float _S592 = - _S582 + _S584;
    float _S593 = _S582 + _S584;
    float _S594 = _S583 + _S587;
    float _S595 = _S586 + _S588;
    float _S596 = z_7 * (_S585 + _S589);
    float _S597 = y_8 * (_S581 + _S589);
    float _S598 = x_17 * (_S581 + _S585);
    float _S599 = z_7 * _S590 + y_8 * _S591 + x_17 * _S592;
    float _S600 = _S490 * _S599;
    float _S601 = w_7 * _S590 + y_8 * _S593 + x_17 * _S594 + _S596 + _S596;
    float _S602 = _S490 * _S601;
    float _S603 = w_7 * _S591 + z_7 * _S593 + x_17 * _S595 + _S597 + _S597;
    float _S604 = _S490 * _S603;
    float _S605 = w_7 * _S592 + z_7 * _S594 + y_8 * _S595 + _S598 + _S598;
    float _S606 = _S490 * _S605;
    float _S607 = quat_7.x * _S599 + quat_7.w * _S601 + quat_7.z * _S603 + quat_7.y * _S605;
    DiffPair_float_0 _S608;
    (&_S608)->primal_0 = _S489;
    (&_S608)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S608, _S607);
    float _S609 = quat_7.x * _S608.differential_0;
    float _S610 = quat_7.w * _S608.differential_0;
    float _S611 = quat_7.z * _S608.differential_0;
    float _S612 = quat_7.y * _S608.differential_0;
    float _S613 = _S602 + _S610 + _S610;
    float _S614 = _S604 + _S611 + _S611;
    float _S615 = _S606 + _S612 + _S612;
    float _S616 = _S600 + _S609 + _S609;
    float3  _S617 = make_float3 (0.0f, 0.0f, v_depth_2);
    float3  _S618 = make_float3 (0.0f);
    float3  _S619 = _S618;
    *&((&_S619)->y) = _S559;
    *&((&_S619)->x) = _S560;
    float3  _S620 = _S617 + _S619;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S621;
    (&_S621)->primal_0 = R_7;
    (&_S621)->differential_0 = _S568;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S622;
    (&_S622)->primal_0 = mean_5;
    (&_S622)->differential_0 = _S618;
    s_bwd_prop_mul_3(&_S621, &_S622, _S620);
    Matrix<float, 3, 3>  _S623 = _S572 + _S573.differential_0 + _S621.differential_0;
    float3  _S624 = _S618;
    *&((&_S624)->z) = _S579.differential_0.rows[int(2)].z;
    *&((&_S624)->y) = _S579.differential_0.rows[int(1)].y;
    *&((&_S624)->x) = _S579.differential_0.rows[int(0)].x;
    float4  _S625 = make_float4 (0.0f);
    *&((&_S625)->w) = _S613;
    *&((&_S625)->z) = _S614;
    *&((&_S625)->y) = _S615;
    *&((&_S625)->x) = _S616;
    *v_mean_2 = _S622.differential_0;
    *v_quat_2 = _S625;
    *v_scale_2 = _S624;
    *v_in_opacity_2 = _S530;
    *v_R_2 = _S623;
    *v_t_2 = _S620;
    return;
}

