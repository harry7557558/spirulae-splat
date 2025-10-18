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

inline __device__ Matrix<float, 3, 3>  quat_to_rotmat(float4  quat_0)
{
    float x_1 = quat_0.y;
    float inv_norm_0 = (F32_rsqrt((x_1 * x_1 + quat_0.z * quat_0.z + quat_0.w * quat_0.w + quat_0.x * quat_0.x)));
    float x_2 = quat_0.y * inv_norm_0;
    float y_0 = quat_0.z * inv_norm_0;
    float z_0 = quat_0.w * inv_norm_0;
    float w_0 = quat_0.x * inv_norm_0;
    float x2_0 = x_2 * x_2;
    float y2_0 = y_0 * y_0;
    float z2_0 = z_0 * z_0;
    float xy_0 = x_2 * y_0;
    float xz_0 = x_2 * z_0;
    float yz_0 = y_0 * z_0;
    float wx_0 = w_0 * x_2;
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_1)
{
    float _S4 = (*left_0).primal_0.rows[int(0)].x * dOut_1.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_1.x;
    float sum_0 = _S4 + (*left_0).primal_0.rows[int(1)].x * dOut_1.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_1.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_1.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_1.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S5 = (*left_0).primal_0.rows[int(0)].y * dOut_1.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_1.x;
    float sum_2 = _S5 + (*left_0).primal_0.rows[int(1)].y * dOut_1.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_1.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_1.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_1.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S6 = (*left_0).primal_0.rows[int(0)].z * dOut_1.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_1.x;
    float sum_4 = _S6 + (*left_0).primal_0.rows[int(1)].z * dOut_1.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_1.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_1.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_1.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_1;
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
        *_slang_vector_get_element_ptr(&result_1, i_0) = sum_6;
        i_0 = i_0 + int(1);
    }
    return result_1;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_1(Matrix<float, 3, 3>  left_2, Matrix<float, 3, 3>  right_2)
{
    Matrix<float, 3, 3>  result_2;
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
            if(c_1 < int(3))
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
                float sum_9 = sum_8 + _slang_vector_get_element(left_2.rows[r_1], i_1) * _slang_vector_get_element(right_2.rows[i_1], c_1);
                i_1 = i_1 + int(1);
                sum_8 = sum_9;
            }
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_1)), c_1) = sum_8;
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_2;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_1(mul_1(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_1, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_3 = quat_1.y;
    float inv_norm_1 = (F32_rsqrt((x_3 * x_3 + quat_1.z * quat_1.z + quat_1.w * quat_1.w + quat_1.x * quat_1.x)));
    float x_4 = quat_1.y * inv_norm_1;
    float y_1 = quat_1.z * inv_norm_1;
    float z_1 = quat_1.w * inv_norm_1;
    float w_1 = quat_1.x * inv_norm_1;
    float x2_1 = x_4 * x_4;
    float y2_1 = y_1 * y_1;
    float z2_1 = z_1 * z_1;
    float xy_1 = x_4 * y_1;
    float xz_1 = x_4 * z_1;
    float yz_1 = y_1 * z_1;
    float wx_1 = w_1 * x_4;
    float wy_1 = w_1 * y_1;
    float wz_1 = w_1 * z_1;
    Matrix<float, 3, 3>  M_0 = mul_1(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_1(M_0, transpose_0(M_0));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_1, float dOut_2)
{
    DiffPair_float_0 _S7 = *dpx_1;
    float _S8;
    if(((*dpx_1).primal_0) < ((*dpy_1).primal_0))
    {
        _S8 = dOut_2;
    }
    else
    {
        if(((*dpx_1).primal_0) > ((*dpy_1).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_2;
        }
    }
    dpx_1->primal_0 = _S7.primal_0;
    dpx_1->differential_0 = _S8;
    DiffPair_float_0 _S9 = *dpy_1;
    if(((*dpy_1).primal_0) < (_S7.primal_0))
    {
        _S8 = dOut_2;
    }
    else
    {
        if(((*dpy_1).primal_0) > ((*dpx_1).primal_0))
        {
            _S8 = 0.0f;
        }
        else
        {
            _S8 = 0.5f * dOut_2;
        }
    }
    dpy_1->primal_0 = _S9.primal_0;
    dpy_1->differential_0 = _S8;
    return;
}

inline __device__ void projection_opaque_triangle_persp(float3  vert0_0, float3  vert1_0, float3  vert2_0, float hardness_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_0, float fy_0, float cx_0, float cy_0, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, float radius_clip_0, int2  * radii_0, float * depth_0, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float * out_hardness_0)
{
    for(;;)
    {
        float3  vert0_c_0 = mul_0(R_2, vert0_0) + t_1;
        float3  vert1_c_0 = mul_0(R_2, vert1_0) + t_1;
        float3  vert2_c_0 = mul_0(R_2, vert2_0) + t_1;
        float3  vert_mean_c_0 = mul_0(R_2, (vert0_0 + vert1_0 + vert2_0) / make_float3 (3.0f)) + t_1;
        float _S10 = vert0_c_0.z;
        bool _S11;
        if(_S10 < near_plane_0)
        {
            _S11 = true;
        }
        else
        {
            _S11 = _S10 > far_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert1_c_0.z) < near_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert1_c_0.z) > far_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert2_c_0.z) < near_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert2_c_0.z) > far_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert_mean_c_0.z) < near_plane_0;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = (vert_mean_c_0.z) > far_plane_0;
        }
        if(_S11)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S10);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (vert1_c_0.z);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (vert2_c_0.z);
        float2  _S12 = make_float2 (fx_0, fy_0);
        float2  _S13 = make_float2 (cx_0, cy_0);
        *uv0_0 = _S12 * *uv0_0 + _S13;
        *uv1_0 = _S12 * *uv1_0 + _S13;
        float2  _S14 = _S12 * *uv2_0 + _S13;
        *uv2_0 = _S14;
        float _S15 = _S14.x;
        float _S16 = float(image_width_0);
        float x_max_0 = (F32_min(((F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S15))) + 0.5f), (_S16)));
        float x_min_0 = (F32_max(((F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S15))) - 0.5f), (0.0f)));
        float _S17 = _S14.y;
        float _S18 = float(image_height_0);
        float y_max_0 = (F32_min(((F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S17))) + 0.5f), (_S18)));
        float y_min_0 = (F32_max(((F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S17))) - 0.5f), (0.0f)));
        float radius_x_0 = (F32_ceil((x_max_0 - x_min_0)));
        float radius_y_0 = (F32_ceil((y_max_0 - y_min_0)));
        if(radius_x_0 <= radius_clip_0)
        {
            _S11 = radius_y_0 <= radius_clip_0;
        }
        else
        {
            _S11 = false;
        }
        if(_S11)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        if(x_max_0 <= 0.0f)
        {
            _S11 = true;
        }
        else
        {
            _S11 = x_min_0 >= _S16;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = y_max_0 <= 0.0f;
        }
        if(_S11)
        {
            _S11 = true;
        }
        else
        {
            _S11 = y_min_0 >= _S18;
        }
        if(_S11)
        {
            *radii_0 = make_int2 (int(0), int(0));
            break;
        }
        *radii_0 = make_int2 (int(radius_x_0), int(radius_y_0));
        *depth_0 = vert_mean_c_0.z;
        *out_hardness_0 = hardness_0;
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  vert0_1, float3  vert1_1, float3  vert2_1, float hardness_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_1, float fy_1, float cx_1, float cy_1, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, float radius_clip_1, int2  * radii_1, float * depth_1, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float * out_hardness_1)
{
    float3  vert0_c_1 = mul_0(R_3, vert0_1) + t_2;
    float3  vert1_c_1 = mul_0(R_3, vert1_1) + t_2;
    float3  vert2_c_1 = mul_0(R_3, vert2_1) + t_2;
    float3  vert_mean_c_1 = mul_0(R_3, (vert0_1 + vert1_1 + vert2_1) / make_float3 (3.0f)) + t_2;
    *uv0_1 = float2 {vert0_c_1.x, vert0_c_1.y} / make_float2 (vert0_c_1.z);
    *uv1_1 = float2 {vert1_c_1.x, vert1_c_1.y} / make_float2 (vert1_c_1.z);
    *uv2_1 = float2 {vert2_c_1.x, vert2_c_1.y} / make_float2 (vert2_c_1.z);
    float2  _S19 = make_float2 (fx_1, fy_1);
    float2  _S20 = make_float2 (cx_1, cy_1);
    *uv0_1 = _S19 * *uv0_1 + _S20;
    *uv1_1 = _S19 * *uv1_1 + _S20;
    float2  _S21 = _S19 * *uv2_1 + _S20;
    *uv2_1 = _S21;
    float _S22 = _S21.x;
    float _S23 = _S21.y;
    *radii_1 = make_int2 (int((F32_ceil(((F32_min(((F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S22))) + 0.5f), (float(image_width_1)))) - (F32_max(((F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S22))) - 0.5f), (0.0f))))))), int((F32_ceil(((F32_min(((F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S23))) + 0.5f), (float(image_height_1)))) - (F32_max(((F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S23))) - 0.5f), (0.0f))))))));
    *depth_1 = vert_mean_c_1.z;
    *out_hardness_1 = hardness_1;
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S24, float3  _S25)
{
    return mul_0(_S24, _S25);
}

inline __device__ float s_primal_ctx_max_0(float _S26, float _S27)
{
    return (F32_max((_S26), (_S27)));
}

inline __device__ float s_primal_ctx_min_0(float _S28, float _S29)
{
    return (F32_min((_S28), (_S29)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S30, DiffPair_float_0 * _S31, float _S32)
{
    _d_max_0(_S30, _S31, _S32);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S33, DiffPair_float_0 * _S34, float _S35)
{
    _d_min_0(_S33, _S34, _S35);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S36, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S37, float3  _S38)
{
    _d_mul_0(_S36, _S37, _S38);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  vert0_2, float3  vert1_2, float3  vert2_2, float hardness_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_2, float fy_2, float cx_2, float cy_2, uint image_width_2, uint image_height_2, float v_depth_0, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float v_out_hardness_0, float3  * v_vert0_0, float3  * v_vert1_0, float3  * v_vert2_0, float * v_hardness_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  vert0_c_2 = s_primal_ctx_mul_0(R_4, vert0_2) + t_3;
    float3  vert_mean_0 = (vert0_2 + vert0_2 + vert0_2) / make_float3 (3.0f);
    float2  _S39 = float2 {vert0_c_2.x, vert0_c_2.y};
    float _S40 = vert0_c_2.z;
    float2  _S41 = make_float2 (_S40);
    float2  _S42 = make_float2 (_S40 * _S40);
    float2  _S43 = make_float2 (fx_2, fy_2);
    float2  _S44 = make_float2 (cx_2, cy_2);
    float2  _S45 = _S43 * (_S39 / make_float2 (_S40)) + _S44;
    float2  _S46 = _S43 * (_S39 / make_float2 (_S40)) + _S44;
    float2  _S47 = _S43 * (_S39 / make_float2 (_S40)) + _S44;
    float _S48 = _S45.x;
    float _S49 = _S46.x;
    float _S50 = s_primal_ctx_max_0(_S48, _S49);
    float _S51 = _S47.x;
    float _S52 = s_primal_ctx_max_0(_S50, _S51) + 0.5f;
    float _S53 = float(image_width_2);
    float _S54 = s_primal_ctx_min_0(_S48, _S49);
    float _S55 = s_primal_ctx_min_0(_S54, _S51) - 0.5f;
    float _S56 = _S45.y;
    float _S57 = _S46.y;
    float _S58 = s_primal_ctx_max_0(_S56, _S57);
    float _S59 = _S47.y;
    float _S60 = s_primal_ctx_max_0(_S58, _S59) + 0.5f;
    float _S61 = float(image_height_2);
    float _S62 = s_primal_ctx_min_0(_S56, _S57);
    DiffPair_float_0 _S63;
    (&_S63)->primal_0 = s_primal_ctx_min_0(_S62, _S59) - 0.5f;
    (&_S63)->differential_0 = 0.0f;
    DiffPair_float_0 _S64;
    (&_S64)->primal_0 = 0.0f;
    (&_S64)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S63, &_S64, -0.0f);
    DiffPair_float_0 _S65;
    (&_S65)->primal_0 = _S62;
    (&_S65)->differential_0 = 0.0f;
    DiffPair_float_0 _S66;
    (&_S66)->primal_0 = _S59;
    (&_S66)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S65, &_S66, _S63.differential_0);
    DiffPair_float_0 _S67;
    (&_S67)->primal_0 = _S56;
    (&_S67)->differential_0 = 0.0f;
    DiffPair_float_0 _S68;
    (&_S68)->primal_0 = _S57;
    (&_S68)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S67, &_S68, _S65.differential_0);
    DiffPair_float_0 _S69;
    (&_S69)->primal_0 = _S60;
    (&_S69)->differential_0 = 0.0f;
    DiffPair_float_0 _S70;
    (&_S70)->primal_0 = _S61;
    (&_S70)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S69, &_S70, 0.0f);
    DiffPair_float_0 _S71;
    (&_S71)->primal_0 = _S58;
    (&_S71)->differential_0 = 0.0f;
    DiffPair_float_0 _S72;
    (&_S72)->primal_0 = _S59;
    (&_S72)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S71, &_S72, _S69.differential_0);
    float _S73 = _S66.differential_0 + _S72.differential_0;
    DiffPair_float_0 _S74;
    (&_S74)->primal_0 = _S56;
    (&_S74)->differential_0 = 0.0f;
    DiffPair_float_0 _S75;
    (&_S75)->primal_0 = _S57;
    (&_S75)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S74, &_S75, _S71.differential_0);
    float _S76 = _S68.differential_0 + _S75.differential_0;
    float _S77 = _S67.differential_0 + _S74.differential_0;
    DiffPair_float_0 _S78;
    (&_S78)->primal_0 = _S55;
    (&_S78)->differential_0 = 0.0f;
    DiffPair_float_0 _S79;
    (&_S79)->primal_0 = 0.0f;
    (&_S79)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S78, &_S79, -0.0f);
    DiffPair_float_0 _S80;
    (&_S80)->primal_0 = _S54;
    (&_S80)->differential_0 = 0.0f;
    DiffPair_float_0 _S81;
    (&_S81)->primal_0 = _S51;
    (&_S81)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S80, &_S81, _S78.differential_0);
    DiffPair_float_0 _S82;
    (&_S82)->primal_0 = _S48;
    (&_S82)->differential_0 = 0.0f;
    DiffPair_float_0 _S83;
    (&_S83)->primal_0 = _S49;
    (&_S83)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S82, &_S83, _S80.differential_0);
    DiffPair_float_0 _S84;
    (&_S84)->primal_0 = _S52;
    (&_S84)->differential_0 = 0.0f;
    DiffPair_float_0 _S85;
    (&_S85)->primal_0 = _S53;
    (&_S85)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S84, &_S85, 0.0f);
    DiffPair_float_0 _S86;
    (&_S86)->primal_0 = _S50;
    (&_S86)->differential_0 = 0.0f;
    DiffPair_float_0 _S87;
    (&_S87)->primal_0 = _S51;
    (&_S87)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S86, &_S87, _S84.differential_0);
    float _S88 = _S81.differential_0 + _S87.differential_0;
    DiffPair_float_0 _S89;
    (&_S89)->primal_0 = _S48;
    (&_S89)->differential_0 = 0.0f;
    DiffPair_float_0 _S90;
    (&_S90)->primal_0 = _S49;
    (&_S90)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S89, &_S90, _S86.differential_0);
    float2  _S91 = _S43 * (v_uv2_0 + make_float2 (_S88, _S73)) / _S42;
    float2  _S92 = _S39 * - _S91;
    float2  _S93 = _S41 * _S91;
    float _S94 = _S92.x + _S92.y;
    float2  _S95 = _S43 * (v_uv1_0 + make_float2 (_S83.differential_0 + _S90.differential_0, _S76)) / _S42;
    float2  _S96 = _S39 * - _S95;
    float2  _S97 = _S41 * _S95;
    float _S98 = _S96.x + _S96.y;
    float2  _S99 = _S43 * (v_uv0_0 + make_float2 (_S82.differential_0 + _S89.differential_0, _S77)) / _S42;
    float2  _S100 = _S39 * - _S99;
    float2  _S101 = _S41 * _S99;
    float _S102 = _S100.x + _S100.y;
    float3  s_diff_vert_mean_c_T_0 = make_float3 (0.0f, 0.0f, v_depth_0);
    Matrix<float, 3, 3>  _S103 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S104;
    (&_S104)->primal_0 = R_4;
    (&_S104)->differential_0 = _S103;
    float3  _S105 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S106;
    (&_S106)->primal_0 = vert_mean_0;
    (&_S106)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S104, &_S106, s_diff_vert_mean_c_T_0);
    float3  _S107 = make_float3 (0.3333333432674408f) * _S106.differential_0;
    float3  s_diff_vert2_c_T_0 = make_float3 (_S93.x, _S93.y, _S94);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S108;
    (&_S108)->primal_0 = R_4;
    (&_S108)->differential_0 = _S103;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S109;
    (&_S109)->primal_0 = vert0_2;
    (&_S109)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S108, &_S109, s_diff_vert2_c_T_0);
    float3  s_diff_vert1_c_T_0 = make_float3 (_S97.x, _S97.y, _S98);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S110;
    (&_S110)->primal_0 = R_4;
    (&_S110)->differential_0 = _S103;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S111;
    (&_S111)->primal_0 = vert0_2;
    (&_S111)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S110, &_S111, s_diff_vert1_c_T_0);
    float3  s_diff_vert0_c_T_0 = make_float3 (_S101.x, _S101.y, _S102);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S112;
    (&_S112)->primal_0 = R_4;
    (&_S112)->differential_0 = _S103;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S113;
    (&_S113)->primal_0 = vert0_2;
    (&_S113)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S112, &_S113, s_diff_vert0_c_T_0);
    float3  _S114 = s_diff_vert_mean_c_T_0 + s_diff_vert2_c_T_0 + s_diff_vert1_c_T_0 + s_diff_vert0_c_T_0;
    Matrix<float, 3, 3>  _S115 = _S104.differential_0 + _S108.differential_0 + _S110.differential_0 + _S112.differential_0;
    float3  _S116 = _S107 + _S109.differential_0;
    float3  _S117 = _S107 + _S111.differential_0;
    *v_vert0_0 = _S107 + _S113.differential_0;
    *v_vert1_0 = _S117;
    *v_vert2_0 = _S116;
    *v_hardness_0 = v_out_hardness_0;
    *v_R_0 = _S115;
    *v_t_0 = _S114;
    return;
}

