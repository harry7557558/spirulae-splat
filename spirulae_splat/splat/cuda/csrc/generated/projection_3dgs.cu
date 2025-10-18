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

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_9)
{
    float _S21 = (F32_exp(((*dpx_5).primal_0))) * dOut_9;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S21;
    return;
}

inline __device__ float3  exp_0(float3  x_9)
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
        *_slang_vector_get_element_ptr(&result_9, i_5) = (F32_exp((_slang_vector_get_element(x_9, i_5))));
        i_5 = i_5 + int(1);
    }
    return result_9;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  dOut_10)
{
    float3  _S22 = exp_0((*dpx_6).primal_0) * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S22;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_7, float dOut_11)
{
    float _S23 = 1.0f / (*dpx_7).primal_0 * dOut_11;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S23;
    return;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_2, float3  scale_1, float in_opacity_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_4, float fy_4, float cx_4, float cy_4, uint image_width_0, uint image_height_0, float eps2d_1, float near_plane_0, float far_plane_0, float radius_clip_0, int2  * radii_0, float * depth_0, float2  * mean2d_4, float3  * conic_0, float * opacity_0)
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
        float3  _S26 = exp_0(scale_1);
        float x_10 = quat_2.y;
        float inv_norm_2 = (F32_rsqrt((x_10 * x_10 + quat_2.z * quat_2.z + quat_2.w * quat_2.w + quat_2.x * quat_2.x)));
        float x_11 = quat_2.y * inv_norm_2;
        float y_3 = quat_2.z * inv_norm_2;
        float z_2 = quat_2.w * inv_norm_2;
        float w_2 = quat_2.x * inv_norm_2;
        float x2_3 = x_11 * x_11;
        float y2_3 = y_3 * y_3;
        float z2_2 = z_2 * z_2;
        float xy_3 = x_11 * y_3;
        float xz_2 = x_11 * z_2;
        float yz_2 = y_3 * z_2;
        float wx_2 = w_2 * x_11;
        float wy_2 = w_2 * y_3;
        float wz_2 = w_2 * z_2;
        Matrix<float, 3, 3>  M_1 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_2), 2.0f * (xy_3 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_3 - wz_2), 1.0f - 2.0f * (x2_3 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S26.x, 0.0f, 0.0f, 0.0f, _S26.y, 0.0f, 0.0f, 0.0f, _S26.z));
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_2, mul_4(M_1, transpose_0(M_1))), transpose_0(R_2));
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

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_3, float3  scale_2, float in_opacity_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_5, float fy_5, float cx_5, float cy_5, uint image_width_1, uint image_height_1, float eps2d_2, float near_plane_1, float far_plane_1, float radius_clip_1, int2  * radii_1, float * depth_1, float2  * mean2d_5, float3  * conic_1, float * opacity_1)
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
        float3  _S35 = exp_0(scale_2);
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
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_3), 2.0f * (xy_4 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_4 - wz_3), 1.0f - 2.0f * (x2_4 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S35.x, 0.0f, 0.0f, 0.0f, _S35.y, 0.0f, 0.0f, 0.0f, _S35.z));
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
        float _S36 = a_1 - b_1;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_5 * (x2_5 * a_1 + y2_5 * b_1), fx_5 * xy_5 * _S36, - fx_5 * mean_c_1.x * x2y2z2_inv_1, fy_5 * xy_5 * _S36, fy_5 * (y2_5 * a_1 + x2_5 * b_1), - fy_5 * mean_c_1.y * x2y2z2_inv_1);
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

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_4, float3  scale_3, float in_opacity_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_6, float fy_6, float cx_6, float cy_6, uint image_width_2, uint image_height_2, float eps2d_3, float near_plane_2, float far_plane_2, float radius_clip_2, int2  * radii_2, float * depth_2, float2  * mean2d_6, float3  * conic_2, float * opacity_2)
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
        float3  _S41 = exp_0(scale_3);
        float x_14 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_14 * x_14 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_15 = quat_4.y * inv_norm_4;
        float y_5 = quat_4.z * inv_norm_4;
        float z_4 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_6 = x_15 * x_15;
        float y2_6 = y_5 * y_5;
        float z2_4 = z_4 * z_4;
        float xy_6 = x_15 * y_5;
        float xz_4 = x_15 * z_4;
        float yz_4 = y_5 * z_4;
        float wx_4 = w_4 * x_15;
        float wy_4 = w_4 * y_5;
        float wz_4 = w_4 * z_4;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_4), 2.0f * (xy_6 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_6 - wz_4), 1.0f - 2.0f * (x2_6 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S41.x, 0.0f, 0.0f, 0.0f, _S41.y, 0.0f, 0.0f, 0.0f, _S41.z));
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_4, mul_4(M_3, transpose_0(M_3))), transpose_0(R_4));
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

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_3, float3  mean_3, float4  quat_5, float3  scale_4, float in_opacity_3, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_7, float fy_7, float cx_7, float cy_7, uint image_width_3, uint image_height_3, float eps2d_4, float near_plane_3, float far_plane_3, float radius_clip_3, int2  * radii_3, float * depth_3, float2  * mean2d_7, float3  * conic_3, float * opacity_3)
{
    float3  mean_c_3 = mul_0(R_5, mean_3) + t_4;
    float3  _S44 = exp_0(scale_4);
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
    Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_5), 2.0f * (xy_7 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_7 - wz_5), 1.0f - 2.0f * (x2_7 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S44.x, 0.0f, 0.0f, 0.0f, _S44.y, 0.0f, 0.0f, 0.0f, _S44.z));
    float _S45 = float(image_width_3);
    float _S46 = float(image_height_3);
    float _S47 = 0.30000001192092896f * (0.5f * _S45 / fx_7);
    float _S48 = 0.30000001192092896f * (0.5f * _S46 / fy_7);
    float rz_3 = 1.0f / mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S45 - cx_7) / fx_7 + _S47), ((F32_max((- (cx_7 / fx_7 + _S47)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S46 - cy_7) / fy_7 + _S48), ((F32_max((- (cy_7 / fy_7 + _S48)), (mean_c_3.y * rz_3))))))) * rz2_3);
    Matrix<float, 2, 2>  covar2d_3 = mul_6(mul_5(J_7, mul_4(mul_4(R_5, mul_4(M_4, transpose_0(M_4))), transpose_0(R_5))), transpose_1(J_7));
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

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_4, float3  mean_4, float4  quat_6, float3  scale_5, float in_opacity_4, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_8, float fy_8, float cx_8, float cy_8, uint image_width_4, uint image_height_4, float eps2d_5, float near_plane_4, float far_plane_4, float radius_clip_4, int2  * radii_4, float * depth_4, float2  * mean2d_8, float3  * conic_4, float * opacity_4)
{
    float3  mean_c_4 = mul_0(R_6, mean_4) + t_5;
    float3  _S51 = exp_0(scale_5);
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
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_6), 2.0f * (xy_8 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_8 - wz_6), 1.0f - 2.0f * (x2_8 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S51.x, 0.0f, 0.0f, 0.0f, _S51.y, 0.0f, 0.0f, 0.0f, _S51.z));
    Matrix<float, 3, 3>  covar_c_3 = mul_4(mul_4(R_6, mul_4(M_5, transpose_0(M_5))), transpose_0(R_6));
    float xy_len_2 = length_0(make_float2 (mean_c_4.x, mean_c_4.y)) + 1.00000001168609742e-07f;
    float theta_2 = (F32_atan2((xy_len_2), (mean_c_4.z + 1.00000001168609742e-07f)));
    *mean2d_8 = make_float2 (mean_c_4.x * fx_8 * theta_2 / xy_len_2 + cx_8, mean_c_4.y * fy_8 * theta_2 / xy_len_2 + cy_8);
    float x2_9 = mean_c_4.x * mean_c_4.x + 1.00000001168609742e-07f;
    float y2_9 = mean_c_4.y * mean_c_4.y;
    float xy_9 = mean_c_4.x * mean_c_4.y;
    float x2y2_2 = x2_9 + y2_9;
    float x2y2z2_inv_2 = 1.0f / (x2y2_2 + mean_c_4.z * mean_c_4.z);
    float b_2 = (F32_atan2((xy_len_2), (mean_c_4.z))) / xy_len_2 / x2y2_2;
    float a_2 = mean_c_4.z * x2y2z2_inv_2 / x2y2_2;
    float _S52 = a_2 - b_2;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_8 * (x2_9 * a_2 + y2_9 * b_2), fx_8 * xy_9 * _S52, - fx_8 * mean_c_4.x * x2y2z2_inv_2, fy_8 * xy_9 * _S52, fy_8 * (y2_9 * a_2 + x2_9 * b_2), - fy_8 * mean_c_4.y * x2y2z2_inv_2);
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

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_5, float3  mean_5, float4  quat_7, float3  scale_6, float in_opacity_5, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_9, float fy_9, float cx_9, float cy_9, uint image_width_5, uint image_height_5, float eps2d_6, float near_plane_5, float far_plane_5, float radius_clip_5, int2  * radii_5, float * depth_5, float2  * mean2d_9, float3  * conic_5, float * opacity_5)
{
    float3  mean_c_5 = mul_0(R_7, mean_5) + t_6;
    float3  _S55 = exp_0(scale_6);
    float x_20 = quat_7.y;
    float inv_norm_7 = (F32_rsqrt((x_20 * x_20 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x)));
    float x_21 = quat_7.y * inv_norm_7;
    float y_8 = quat_7.z * inv_norm_7;
    float z_7 = quat_7.w * inv_norm_7;
    float w_7 = quat_7.x * inv_norm_7;
    float x2_10 = x_21 * x_21;
    float y2_10 = y_8 * y_8;
    float z2_7 = z_7 * z_7;
    float xy_10 = x_21 * y_8;
    float xz_7 = x_21 * z_7;
    float yz_7 = y_8 * z_7;
    float wx_7 = w_7 * x_21;
    float wy_7 = w_7 * y_8;
    float wz_7 = w_7 * z_7;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_7), 2.0f * (xy_10 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_10 - wz_7), 1.0f - 2.0f * (x2_10 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S55.x, 0.0f, 0.0f, 0.0f, _S55.y, 0.0f, 0.0f, 0.0f, _S55.z));
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_9, mul_4(mul_4(R_7, mul_4(M_6, transpose_0(M_6))), transpose_0(R_7))), transpose_1(J_9));
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

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_8, float3  scale_7, float in_opacity_6, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_10, float fy_10, float cx_10, float cy_10, uint image_width_6, uint image_height_6, float eps2d_7, float v_depth_0, float2  v_mean2d_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_8, mean_6) + t_7;
    float3  _S103 = s_primal_ctx_exp_0(scale_7);
    float _S104 = quat_8.y;
    float _S105 = _S104 * _S104 + quat_8.z * quat_8.z + quat_8.w * quat_8.w + quat_8.x * quat_8.x;
    float _S106 = s_primal_ctx_rsqrt_0(_S105);
    float x_22 = quat_8.y * _S106;
    float y_9 = quat_8.z * _S106;
    float z_8 = quat_8.w * _S106;
    float w_8 = quat_8.x * _S106;
    float x2_11 = x_22 * x_22;
    float y2_11 = y_9 * y_9;
    float z2_8 = z_8 * z_8;
    float xy_11 = x_22 * y_9;
    float xz_8 = x_22 * z_8;
    float yz_8 = y_9 * z_8;
    float wx_8 = w_8 * x_22;
    float wy_8 = w_8 * y_9;
    float wz_8 = w_8 * z_8;
    Matrix<float, 3, 3>  _S107 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_8), 2.0f * (xy_11 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_11 - wz_8), 1.0f - 2.0f * (x2_11 + z2_8), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_11 + y2_11)));
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
    float _S267 = z_8 * (_S256 + _S260);
    float _S268 = y_9 * (_S252 + _S260);
    float _S269 = x_22 * (_S252 + _S256);
    float _S270 = z_8 * _S261 + y_9 * _S262 + x_22 * _S263;
    float _S271 = _S106 * _S270;
    float _S272 = w_8 * _S261 + y_9 * _S264 + x_22 * _S265 + _S267 + _S267;
    float _S273 = _S106 * _S272;
    float _S274 = w_8 * _S262 + z_8 * _S264 + x_22 * _S266 + _S268 + _S268;
    float _S275 = _S106 * _S274;
    float _S276 = w_8 * _S263 + z_8 * _S265 + y_9 * _S266 + _S269 + _S269;
    float _S277 = _S106 * _S276;
    float _S278 = quat_8.x * _S270 + quat_8.w * _S272 + quat_8.z * _S274 + quat_8.y * _S276;
    DiffPair_float_0 _S279;
    (&_S279)->primal_0 = _S105;
    (&_S279)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S279, _S278);
    float _S280 = quat_8.x * _S279.differential_0;
    float _S281 = quat_8.w * _S279.differential_0;
    float _S282 = quat_8.z * _S279.differential_0;
    float _S283 = quat_8.y * _S279.differential_0;
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
    (&_S290)->primal_0 = scale_7;
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

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8, float _s_dOut_0)
{
    float _S305 = (*dpx_8).primal_0.x;
    float _S306 = (*dpx_8).primal_0.y;
    DiffPair_float_0 _S307;
    (&_S307)->primal_0 = _S305 * _S305 + _S306 * _S306;
    (&_S307)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S307, _s_dOut_0);
    float _S308 = (*dpx_8).primal_0.y * _S307.differential_0;
    float _S309 = _S308 + _S308;
    float _S310 = (*dpx_8).primal_0.x * _S307.differential_0;
    float _S311 = _S310 + _S310;
    float2  _S312 = make_float2 (0.0f);
    *&((&_S312)->y) = _S309;
    *&((&_S312)->x) = _S311;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S312;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S313, float _S314)
{
    s_bwd_prop_length_impl_0(_S313, _S314);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_9, float3  scale_8, float in_opacity_7, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_11, float fy_11, float cx_11, float cy_11, uint image_width_7, uint image_height_7, float eps2d_8, float v_depth_1, float2  v_mean2d_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_7 = s_primal_ctx_mul_0(R_9, mean_7) + t_8;
    float3  _S315 = s_primal_ctx_exp_0(scale_8);
    float _S316 = quat_9.y;
    float _S317 = _S316 * _S316 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S318 = s_primal_ctx_rsqrt_0(_S317);
    float x_23 = quat_9.y * _S318;
    float y_10 = quat_9.z * _S318;
    float z_9 = quat_9.w * _S318;
    float w_9 = quat_9.x * _S318;
    float x2_12 = x_23 * x_23;
    float y2_12 = y_10 * y_10;
    float z2_9 = z_9 * z_9;
    float xy_12 = x_23 * y_10;
    float xz_9 = x_23 * z_9;
    float yz_9 = y_10 * z_9;
    float wx_9 = w_9 * x_23;
    float wy_9 = w_9 * y_10;
    float wz_9 = w_9 * z_9;
    Matrix<float, 3, 3>  _S319 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_9), 2.0f * (xy_12 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_12 - wz_9), 1.0f - 2.0f * (x2_12 + z2_9), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_12 + y2_12)));
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
    float x2_13 = mean_c_7.x * mean_c_7.x + 1.00000001168609742e-07f;
    float y2_13 = mean_c_7.y * mean_c_7.y;
    float xy_13 = mean_c_7.x * mean_c_7.y;
    float x2y2_3 = x2_13 + y2_13;
    float _S334 = x2y2_3 + mean_c_7.z * mean_c_7.z;
    float x2y2z2_inv_3 = 1.0f / _S334;
    float _S335 = _S334 * _S334;
    float _S336 = s_primal_ctx_atan2_0(xy_len_3, mean_c_7.z);
    float _S337 = _S336 / xy_len_3;
    float b_3 = _S337 / x2y2_3;
    float _S338 = x2y2_3 * x2y2_3;
    float _S339 = mean_c_7.z * x2y2z2_inv_3;
    float a_3 = _S339 / x2y2_3;
    float _S340 = fx_11 * xy_13;
    float _S341 = a_3 - b_3;
    float _S342 = - fx_11;
    float _S343 = _S342 * mean_c_7.x;
    float _S344 = fy_11 * xy_13;
    float _S345 = - fy_11;
    float _S346 = _S345 * mean_c_7.y;
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_11 * (x2_13 * a_3 + y2_13 * b_3), _S340 * _S341, _S343 * x2y2z2_inv_3, _S344 * _S341, fy_11 * (y2_13 * a_3 + x2_13 * b_3), _S346 * x2y2z2_inv_3);
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
    float _S442 = (y2_13 * _S431 + _S437 + x2_13 * _S439) / _S338;
    float _S443 = _S339 * - _S442;
    float _S444 = x2y2_3 * _S442;
    float _S445 = mean_c_7.z * _S444;
    float _S446 = x2y2z2_inv_3 * _S444;
    float _S447 = (x2_13 * _S431 + - _S437 + y2_13 * _S439) / _S338;
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
    float _S504 = z_9 * (_S493 + _S497);
    float _S505 = y_10 * (_S489 + _S497);
    float _S506 = x_23 * (_S489 + _S493);
    float _S507 = z_9 * _S498 + y_10 * _S499 + x_23 * _S500;
    float _S508 = _S318 * _S507;
    float _S509 = w_9 * _S498 + y_10 * _S501 + x_23 * _S502 + _S504 + _S504;
    float _S510 = _S318 * _S509;
    float _S511 = w_9 * _S499 + z_9 * _S501 + x_23 * _S503 + _S505 + _S505;
    float _S512 = _S318 * _S511;
    float _S513 = w_9 * _S500 + z_9 * _S502 + y_10 * _S503 + _S506 + _S506;
    float _S514 = _S318 * _S513;
    float _S515 = quat_9.x * _S507 + quat_9.w * _S509 + quat_9.z * _S511 + quat_9.y * _S513;
    DiffPair_float_0 _S516;
    (&_S516)->primal_0 = _S317;
    (&_S516)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S516, _S515);
    float _S517 = quat_9.x * _S516.differential_0;
    float _S518 = quat_9.w * _S516.differential_0;
    float _S519 = quat_9.z * _S516.differential_0;
    float _S520 = quat_9.y * _S516.differential_0;
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
    (&_S527)->primal_0 = scale_8;
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

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_10, float3  scale_9, float in_opacity_8, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_12, float fy_12, float cx_12, float cy_12, uint image_width_8, uint image_height_8, float eps2d_9, float v_depth_2, float2  v_mean2d_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float3  _S537 = s_primal_ctx_exp_0(scale_9);
    float _S538 = quat_10.y;
    float _S539 = _S538 * _S538 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S540 = s_primal_ctx_rsqrt_0(_S539);
    float x_24 = quat_10.y * _S540;
    float y_11 = quat_10.z * _S540;
    float z_10 = quat_10.w * _S540;
    float w_10 = quat_10.x * _S540;
    float x2_14 = x_24 * x_24;
    float y2_14 = y_11 * y_11;
    float z2_10 = z_10 * z_10;
    float xy_14 = x_24 * y_11;
    float xz_10 = x_24 * z_10;
    float yz_10 = y_11 * z_10;
    float wx_10 = w_10 * x_24;
    float wy_10 = w_10 * y_11;
    float wz_10 = w_10 * z_10;
    Matrix<float, 3, 3>  _S541 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_10), 2.0f * (xy_14 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_14 - wz_10), 1.0f - 2.0f * (x2_14 + z2_10), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_14 + y2_14)));
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
    float _S656 = z_10 * (_S645 + _S649);
    float _S657 = y_11 * (_S641 + _S649);
    float _S658 = x_24 * (_S641 + _S645);
    float _S659 = z_10 * _S650 + y_11 * _S651 + x_24 * _S652;
    float _S660 = _S540 * _S659;
    float _S661 = w_10 * _S650 + y_11 * _S653 + x_24 * _S654 + _S656 + _S656;
    float _S662 = _S540 * _S661;
    float _S663 = w_10 * _S651 + z_10 * _S653 + x_24 * _S655 + _S657 + _S657;
    float _S664 = _S540 * _S663;
    float _S665 = w_10 * _S652 + z_10 * _S654 + y_11 * _S655 + _S658 + _S658;
    float _S666 = _S540 * _S665;
    float _S667 = quat_10.x * _S659 + quat_10.w * _S661 + quat_10.z * _S663 + quat_10.y * _S665;
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = _S539;
    (&_S668)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S668, _S667);
    float _S669 = quat_10.x * _S668.differential_0;
    float _S670 = quat_10.w * _S668.differential_0;
    float _S671 = quat_10.z * _S668.differential_0;
    float _S672 = quat_10.y * _S668.differential_0;
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
    (&_S679)->primal_0 = scale_9;
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

