#pragma once

#include "slang.cuh"

inline __device__ float3  sh_coeffs_to_color(int degree_0, float3  viewdir_0, float3  coeff_dc_0, float3  * coeffs_0)
{
    float3  colors_0 = make_float3 (0.282094806432724f) * coeff_dc_0;
    if(degree_0 <= int(0))
    {
        return colors_0;
    }
    float _S1 = viewdir_0.x;
    float _S2 = viewdir_0.y;
    float _S3 = viewdir_0.z;
    float inv_norm_0 = (F32_rsqrt((_S1 * _S1 + _S2 * _S2 + _S3 * _S3)));
    float x_0 = _S1 * inv_norm_0;
    float y_0 = _S2 * inv_norm_0;
    float z_0 = _S3 * inv_norm_0;
    float3  colors_1 = colors_0 + make_float3 (0.48860251903533936f) * (make_float3 (- y_0) * *(coeffs_0 + int(0)) + make_float3 (z_0) * *(coeffs_0 + int(1)) - make_float3 (x_0) * *(coeffs_0 + int(2)));
    if(degree_0 <= int(1))
    {
        return colors_1;
    }
    float z2_0 = z_0 * z_0;
    float fTmp0B_0 = -1.09254848957061768f * z_0;
    float fC1_0 = x_0 * x_0 - y_0 * y_0;
    float fS1_0 = 2.0f * x_0 * y_0;
    float pSH6_0 = 0.94617468118667603f * z2_0 - 0.31539157032966614f;
    float3  colors_2 = colors_1 + (make_float3 (0.54627424478530884f * fS1_0) * *(coeffs_0 + int(3)) + make_float3 (fTmp0B_0 * y_0) * *(coeffs_0 + int(4)) + make_float3 (pSH6_0) * *(coeffs_0 + int(5)) + make_float3 (fTmp0B_0 * x_0) * *(coeffs_0 + int(6)) + make_float3 (0.54627424478530884f * fC1_0) * *(coeffs_0 + int(7)));
    if(degree_0 <= int(2))
    {
        return colors_2;
    }
    float fTmp0C_0 = -2.28522896766662598f * z2_0 + 0.4570457935333252f;
    float fTmp1B_0 = 1.44530570507049561f * z_0;
    float fC2_0 = x_0 * fC1_0 - y_0 * fS1_0;
    float fS2_0 = x_0 * fS1_0 + y_0 * fC1_0;
    float pSH12_0 = z_0 * (1.86588168144226074f * z2_0 - 1.11952900886535645f);
    float3  colors_3 = colors_2 + (make_float3 (-0.59004360437393188f * fS2_0) * *(coeffs_0 + int(8)) + make_float3 (fTmp1B_0 * fS1_0) * *(coeffs_0 + int(9)) + make_float3 (fTmp0C_0 * y_0) * *(coeffs_0 + int(10)) + make_float3 (pSH12_0) * *(coeffs_0 + int(11)) + make_float3 (fTmp0C_0 * x_0) * *(coeffs_0 + int(12)) + make_float3 (fTmp1B_0 * fC1_0) * *(coeffs_0 + int(13)) + make_float3 (-0.59004360437393188f * fC2_0) * *(coeffs_0 + int(14)));
    if(degree_0 <= int(3))
    {
        return colors_3;
    }
    float fTmp0D_0 = z_0 * (-4.68332576751708984f * z2_0 + 2.00713968276977539f);
    float fTmp1C_0 = 3.31161141395568848f * z2_0 - 0.47308734059333801f;
    float fTmp2B_0 = -1.77013075351715088f * z_0;
    return colors_3 + (make_float3 (0.62583571672439575f * (x_0 * fS2_0 + y_0 * fC2_0)) * *(coeffs_0 + int(15)) + make_float3 (fTmp2B_0 * fS2_0) * *(coeffs_0 + int(16)) + make_float3 (fTmp1C_0 * fS1_0) * *(coeffs_0 + int(17)) + make_float3 (fTmp0D_0 * y_0) * *(coeffs_0 + int(18)) + make_float3 (1.9843134880065918f * z_0 * pSH12_0 - 1.00623059272766113f * pSH6_0) * *(coeffs_0 + int(19)) + make_float3 (fTmp0D_0 * x_0) * *(coeffs_0 + int(20)) + make_float3 (fTmp1C_0 * fC1_0) * *(coeffs_0 + int(21)) + make_float3 (fTmp2B_0 * fC2_0) * *(coeffs_0 + int(22)) + make_float3 (0.62583571672439575f * (x_0 * fC2_0 - y_0 * fS2_0)) * *(coeffs_0 + int(23)));
}

inline __device__ float dot_0(float3  x_1, float3  y_1)
{
    int i_0 = int(0);
    float result_0 = 0.0f;
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_1 = result_0 + _slang_vector_get_element(x_1, i_0) * _slang_vector_get_element(y_1, i_0);
        i_0 = i_0 + int(1);
        result_0 = result_1;
    }
    return result_0;
}

inline __device__ void sh_coeffs_to_color_vjp(int degree_1, float3  dir_0, float3  coeff_dc_1, float3  * coeffs_1, float3  v_colors_0, float3  * v_coeff_dc_0, float3  * v_coeffs_0, float3  * v_dir_0)
{
    for(;;)
    {
        *v_coeff_dc_0 = make_float3 (0.282094806432724f) * v_colors_0;
        if(degree_1 <= int(0))
        {
            int3  _S4 = make_int3 (int(0));
            float3  _S5 = make_float3 ((float)_S4.x, (float)_S4.y, (float)_S4.z);
            *v_dir_0 = _S5;
            break;
        }
        float _S6 = dir_0.x;
        float _S7 = dir_0.y;
        float _S8 = dir_0.z;
        float inorm_0 = (F32_rsqrt((_S6 * _S6 + _S7 * _S7 + _S8 * _S8)));
        float x_2 = _S6 * inorm_0;
        float y_2 = _S7 * inorm_0;
        float z_1 = _S8 * inorm_0;
        *(v_coeffs_0 + int(0)) = make_float3 (-0.48860251903533936f * y_2) * v_colors_0;
        *(v_coeffs_0 + int(1)) = make_float3 (0.48860251903533936f * z_1) * v_colors_0;
        *(v_coeffs_0 + int(2)) = make_float3 (-0.48860251903533936f * x_2) * v_colors_0;
        float _S9 = -0.48860251903533936f * dot_0(*(coeffs_1 + int(2)), v_colors_0);
        float _S10 = -0.48860251903533936f * dot_0(*(coeffs_1 + int(0)), v_colors_0);
        float _S11 = 0.48860251903533936f * dot_0(*(coeffs_1 + int(1)), v_colors_0);
        if(degree_1 <= int(1))
        {
            float3  dir_n_0 = make_float3 (x_2, y_2, z_1);
            float3  v_dir_n_0 = make_float3 (_S9, _S10, _S11);
            *v_dir_0 = (v_dir_n_0 - make_float3 (dot_0(v_dir_n_0, dir_n_0)) * dir_n_0) * make_float3 (inorm_0);
            break;
        }
        float z2_1 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_2 * x_2 - y_2 * y_2;
        float _S12 = 2.0f * x_2;
        float fS1_1 = _S12 * y_2;
        float pSH6_1 = 0.94617468118667603f * z2_1 - 0.31539157032966614f;
        float pSH7_0 = fTmp0B_1 * x_2;
        float pSH5_0 = fTmp0B_1 * y_2;
        float pSH8_0 = 0.54627424478530884f * fC1_1;
        float pSH4_0 = 0.54627424478530884f * fS1_1;
        *(v_coeffs_0 + int(3)) = make_float3 (pSH4_0) * v_colors_0;
        *(v_coeffs_0 + int(4)) = make_float3 (pSH5_0) * v_colors_0;
        *(v_coeffs_0 + int(5)) = make_float3 (pSH6_1) * v_colors_0;
        *(v_coeffs_0 + int(6)) = make_float3 (pSH7_0) * v_colors_0;
        *(v_coeffs_0 + int(7)) = make_float3 (pSH8_0) * v_colors_0;
        float fC1_y_0 = -2.0f * y_2;
        float fS1_x_0 = 2.0f * y_2;
        float pSH6_z_0 = 1.89234936237335205f * z_1;
        float pSH8_x_0 = 0.54627424478530884f * _S12;
        float3  * _S13 = coeffs_1 + int(3);
        float3  * _S14 = coeffs_1 + int(7);
        float3  * _S15 = coeffs_1 + int(6);
        float v_x_0 = _S9 + dot_0(v_colors_0, make_float3 (0.54627424478530884f * fS1_x_0) * *_S13 + make_float3 (pSH8_x_0) * *_S14 + make_float3 (fTmp0B_1) * *_S15);
        float3  * _S16 = coeffs_1 + int(4);
        float v_y_0 = _S10 + dot_0(v_colors_0, make_float3 (pSH8_x_0) * *_S13 + make_float3 (0.54627424478530884f * fC1_y_0) * *_S14 + make_float3 (fTmp0B_1) * *_S16);
        float v_z_0 = _S11 + dot_0(v_colors_0, make_float3 (pSH6_z_0) * *(coeffs_1 + int(5)) + make_float3 (-1.09254848957061768f * x_2) * *_S15 + make_float3 (-1.09254848957061768f * y_2) * *_S16);
        if(degree_1 <= int(2))
        {
            float3  dir_n_1 = make_float3 (x_2, y_2, z_1);
            float3  v_dir_n_1 = make_float3 (v_x_0, v_y_0, v_z_0);
            *v_dir_0 = (v_dir_n_1 - make_float3 (dot_0(v_dir_n_1, dir_n_1)) * dir_n_1) * make_float3 (inorm_0);
            break;
        }
        float fTmp0C_1 = -2.28522896766662598f * z2_1 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        float fC2_1 = x_2 * fC1_1 - y_2 * fS1_1;
        float fS2_1 = x_2 * fS1_1 + y_2 * fC1_1;
        float pSH12_1 = z_1 * (1.86588168144226074f * z2_1 - 1.11952900886535645f);
        float pSH13_0 = fTmp0C_1 * x_2;
        float pSH11_0 = fTmp0C_1 * y_2;
        float pSH14_0 = fTmp1B_1 * fC1_1;
        float pSH10_0 = fTmp1B_1 * fS1_1;
        float pSH15_0 = -0.59004360437393188f * fC2_1;
        float pSH9_0 = -0.59004360437393188f * fS2_1;
        *(v_coeffs_0 + int(8)) = make_float3 (pSH9_0) * v_colors_0;
        *(v_coeffs_0 + int(9)) = make_float3 (pSH10_0) * v_colors_0;
        *(v_coeffs_0 + int(10)) = make_float3 (pSH11_0) * v_colors_0;
        *(v_coeffs_0 + int(11)) = make_float3 (pSH12_1) * v_colors_0;
        *(v_coeffs_0 + int(12)) = make_float3 (pSH13_0) * v_colors_0;
        *(v_coeffs_0 + int(13)) = make_float3 (pSH14_0) * v_colors_0;
        *(v_coeffs_0 + int(14)) = make_float3 (pSH15_0) * v_colors_0;
        float fTmp0C_z_0 = -4.57045793533325195f * z_1;
        float _S17 = x_2 * _S12;
        float fC2_x_0 = fC1_1 + _S17 - y_2 * fS1_x_0;
        float _S18 = y_2 * _S12;
        float fC2_y_0 = x_2 * fC1_y_0 - fS1_1 - _S18;
        float fS2_x_0 = fS1_1 + x_2 * fS1_x_0 + _S18;
        float fS2_y_0 = _S17 + fC1_1 + y_2 * fC1_y_0;
        float pSH12_z_0 = 5.59764480590820312f * z2_1 - 1.11952900886535645f;
        float pSH14_x_0 = fTmp1B_1 * _S12;
        float3  * _S19 = coeffs_1 + int(8);
        float3  * _S20 = coeffs_1 + int(14);
        float3  * _S21 = coeffs_1 + int(9);
        float3  * _S22 = coeffs_1 + int(13);
        float3  * _S23 = coeffs_1 + int(12);
        float v_x_1 = v_x_0 + dot_0(v_colors_0, make_float3 (-0.59004360437393188f * fS2_x_0) * *_S19 + make_float3 (-0.59004360437393188f * fC2_x_0) * *_S20 + make_float3 (fTmp1B_1 * fS1_x_0) * *_S21 + make_float3 (pSH14_x_0) * *_S22 + make_float3 (fTmp0C_1) * *_S23);
        float3  * _S24 = coeffs_1 + int(10);
        float v_y_1 = v_y_0 + dot_0(v_colors_0, make_float3 (-0.59004360437393188f * fS2_y_0) * *_S19 + make_float3 (-0.59004360437393188f * fC2_y_0) * *_S20 + make_float3 (pSH14_x_0) * *_S21 + make_float3 (fTmp1B_1 * fC1_y_0) * *_S22 + make_float3 (fTmp0C_1) * *_S24);
        float v_z_1 = v_z_0 + dot_0(v_colors_0, make_float3 (pSH12_z_0) * *(coeffs_1 + int(11)) + make_float3 (fTmp0C_z_0 * x_2) * *_S23 + make_float3 (fTmp0C_z_0 * y_2) * *_S24 + make_float3 (1.44530570507049561f * fC1_1) * *_S22 + make_float3 (1.44530570507049561f * fS1_1) * *_S21);
        if(degree_1 <= int(3))
        {
            float3  dir_n_2 = make_float3 (x_2, y_2, z_1);
            float3  v_dir_n_2 = make_float3 (v_x_1, v_y_1, v_z_1);
            *v_dir_0 = (v_dir_n_2 - make_float3 (dot_0(v_dir_n_2, dir_n_2)) * dir_n_2) * make_float3 (inorm_0);
            break;
        }
        float fTmp0D_1 = z_1 * (-4.68332576751708984f * z2_1 + 2.00713968276977539f);
        float fTmp1C_1 = 3.31161141395568848f * z2_1 - 0.47308734059333801f;
        float fTmp2B_1 = -1.77013075351715088f * z_1;
        float pSH20_0 = 1.9843134880065918f * z_1 * pSH12_1 + -1.00623059272766113f * pSH6_1;
        float pSH21_0 = fTmp0D_1 * x_2;
        float pSH19_0 = fTmp0D_1 * y_2;
        float pSH22_0 = fTmp1C_1 * fC1_1;
        float pSH18_0 = fTmp1C_1 * fS1_1;
        float pSH23_0 = fTmp2B_1 * fC2_1;
        float pSH17_0 = fTmp2B_1 * fS2_1;
        float pSH24_0 = 0.62583571672439575f * (x_2 * fC2_1 - y_2 * fS2_1);
        float pSH16_0 = 0.62583571672439575f * (x_2 * fS2_1 + y_2 * fC2_1);
        *(v_coeffs_0 + int(15)) = make_float3 (pSH16_0) * v_colors_0;
        *(v_coeffs_0 + int(16)) = make_float3 (pSH17_0) * v_colors_0;
        *(v_coeffs_0 + int(17)) = make_float3 (pSH18_0) * v_colors_0;
        *(v_coeffs_0 + int(18)) = make_float3 (pSH19_0) * v_colors_0;
        *(v_coeffs_0 + int(19)) = make_float3 (pSH20_0) * v_colors_0;
        *(v_coeffs_0 + int(20)) = make_float3 (pSH21_0) * v_colors_0;
        *(v_coeffs_0 + int(21)) = make_float3 (pSH22_0) * v_colors_0;
        *(v_coeffs_0 + int(22)) = make_float3 (pSH23_0) * v_colors_0;
        *(v_coeffs_0 + int(23)) = make_float3 (pSH24_0) * v_colors_0;
        float fTmp0D_z_0 = -14.04997730255126953f * z2_1 + 2.00713968276977539f;
        float fTmp1C_z_0 = 6.62322282791137695f * z_1;
        float pSH22_x_0 = fTmp1C_1 * _S12;
        float3  * _S25 = coeffs_1 + int(15);
        float3  * _S26 = coeffs_1 + int(23);
        float3  * _S27 = coeffs_1 + int(16);
        float3  * _S28 = coeffs_1 + int(22);
        float3  * _S29 = coeffs_1 + int(17);
        float3  * _S30 = coeffs_1 + int(21);
        float3  * _S31 = coeffs_1 + int(20);
        float3  * _S32 = coeffs_1 + int(18);
        float3  dir_n_3 = make_float3 (x_2, y_2, z_1);
        float3  v_dir_n_3 = make_float3 (v_x_1 + dot_0(v_colors_0, make_float3 (0.62583571672439575f * (fS2_1 + y_2 * fC2_x_0 + x_2 * fS2_x_0)) * *_S25 + make_float3 (0.62583571672439575f * (fC2_1 + x_2 * fC2_x_0 - y_2 * fS2_x_0)) * *_S26 + make_float3 (fTmp2B_1 * fS2_x_0) * *_S27 + make_float3 (fTmp2B_1 * fC2_x_0) * *_S28 + make_float3 (fTmp1C_1 * fS1_x_0) * *_S29 + make_float3 (pSH22_x_0) * *_S30 + make_float3 (fTmp0D_1) * *_S31), v_y_1 + dot_0(v_colors_0, make_float3 (0.62583571672439575f * (x_2 * fS2_y_0 + fC2_1 + y_2 * fC2_y_0)) * *_S25 + make_float3 (0.62583571672439575f * (x_2 * fC2_y_0 - fS2_1 - y_2 * fS2_y_0)) * *_S26 + make_float3 (fTmp2B_1 * fS2_y_0) * *_S27 + make_float3 (fTmp2B_1 * fC2_y_0) * *_S28 + make_float3 (pSH22_x_0) * *_S29 + make_float3 (fTmp1C_1 * fC1_y_0) * *_S30 + make_float3 (fTmp0D_1) * *_S32), v_z_1 + dot_0(v_colors_0, make_float3 (1.9843134880065918f * (pSH12_1 + z_1 * pSH12_z_0) + -1.00623059272766113f * pSH6_z_0) * *(coeffs_1 + int(19)) + make_float3 (fTmp0D_z_0 * x_2) * *_S31 + make_float3 (fTmp0D_z_0 * y_2) * *_S32 + make_float3 (fTmp1C_z_0 * fC1_1) * *_S30 + make_float3 (fTmp1C_z_0 * fS1_1) * *_S29 + make_float3 (-1.77013075351715088f * fC2_1) * *_S28 + make_float3 (-1.77013075351715088f * fS2_1) * *_S27));
        *v_dir_0 = (v_dir_n_3 - make_float3 (dot_0(v_dir_n_3, dir_n_3)) * dir_n_3) * make_float3 (inorm_0);
        break;
    }
    return;
}

inline __device__ void sh_coeffs_to_color_vjp_inplace(int degree_2, float3  dir_1, float3  coeff_dc_2, float3  * coeffs_2, float3  v_colors_1, float3  * v_coeff_dc_1, float3  * v_coeffs_1, float3  * v_dir_1)
{
    for(;;)
    {
        *v_coeff_dc_1 = *v_coeff_dc_1 + make_float3 (0.282094806432724f) * v_colors_1;
        if(degree_2 <= int(0))
        {
            break;
        }
        float _S33 = dir_1.x;
        float _S34 = dir_1.y;
        float _S35 = dir_1.z;
        float inorm_1 = (F32_rsqrt((_S33 * _S33 + _S34 * _S34 + _S35 * _S35)));
        float x_3 = _S33 * inorm_1;
        float y_3 = _S34 * inorm_1;
        float z_2 = _S35 * inorm_1;
        float3  * _S36 = v_coeffs_1 + int(0);
        *_S36 = *_S36 + make_float3 (-0.48860251903533936f * y_3) * v_colors_1;
        float3  * _S37 = v_coeffs_1 + int(1);
        *_S37 = *_S37 + make_float3 (0.48860251903533936f * z_2) * v_colors_1;
        float3  * _S38 = v_coeffs_1 + int(2);
        *_S38 = *_S38 + make_float3 (-0.48860251903533936f * x_3) * v_colors_1;
        float _S39 = -0.48860251903533936f * dot_0(*(coeffs_2 + int(2)), v_colors_1);
        float _S40 = -0.48860251903533936f * dot_0(*(coeffs_2 + int(0)), v_colors_1);
        float _S41 = 0.48860251903533936f * dot_0(*(coeffs_2 + int(1)), v_colors_1);
        if(degree_2 <= int(1))
        {
            float3  dir_n_4 = make_float3 (x_3, y_3, z_2);
            float3  v_dir_n_4 = make_float3 (_S39, _S40, _S41);
            *v_dir_1 = *v_dir_1 + (v_dir_n_4 - make_float3 (dot_0(v_dir_n_4, dir_n_4)) * dir_n_4) * make_float3 (inorm_1);
            break;
        }
        float z2_2 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_3 * x_3 - y_3 * y_3;
        float _S42 = 2.0f * x_3;
        float fS1_2 = _S42 * y_3;
        float pSH6_2 = 0.94617468118667603f * z2_2 - 0.31539157032966614f;
        float pSH7_1 = fTmp0B_2 * x_3;
        float pSH5_1 = fTmp0B_2 * y_3;
        float pSH8_1 = 0.54627424478530884f * fC1_2;
        float pSH4_1 = 0.54627424478530884f * fS1_2;
        float3  * _S43 = v_coeffs_1 + int(3);
        *_S43 = *_S43 + make_float3 (pSH4_1) * v_colors_1;
        float3  * _S44 = v_coeffs_1 + int(4);
        *_S44 = *_S44 + make_float3 (pSH5_1) * v_colors_1;
        float3  * _S45 = v_coeffs_1 + int(5);
        *_S45 = *_S45 + make_float3 (pSH6_2) * v_colors_1;
        float3  * _S46 = v_coeffs_1 + int(6);
        *_S46 = *_S46 + make_float3 (pSH7_1) * v_colors_1;
        float3  * _S47 = v_coeffs_1 + int(7);
        *_S47 = *_S47 + make_float3 (pSH8_1) * v_colors_1;
        float fC1_y_1 = -2.0f * y_3;
        float fS1_x_1 = 2.0f * y_3;
        float pSH6_z_1 = 1.89234936237335205f * z_2;
        float pSH8_x_1 = 0.54627424478530884f * _S42;
        float3  * _S48 = coeffs_2 + int(3);
        float3  * _S49 = coeffs_2 + int(7);
        float3  * _S50 = coeffs_2 + int(6);
        float v_x_2 = _S39 + dot_0(v_colors_1, make_float3 (0.54627424478530884f * fS1_x_1) * *_S48 + make_float3 (pSH8_x_1) * *_S49 + make_float3 (fTmp0B_2) * *_S50);
        float3  * _S51 = coeffs_2 + int(4);
        float v_y_2 = _S40 + dot_0(v_colors_1, make_float3 (pSH8_x_1) * *_S48 + make_float3 (0.54627424478530884f * fC1_y_1) * *_S49 + make_float3 (fTmp0B_2) * *_S51);
        float v_z_2 = _S41 + dot_0(v_colors_1, make_float3 (pSH6_z_1) * *(coeffs_2 + int(5)) + make_float3 (-1.09254848957061768f * x_3) * *_S50 + make_float3 (-1.09254848957061768f * y_3) * *_S51);
        if(degree_2 <= int(2))
        {
            float3  dir_n_5 = make_float3 (x_3, y_3, z_2);
            float3  v_dir_n_5 = make_float3 (v_x_2, v_y_2, v_z_2);
            *v_dir_1 = *v_dir_1 + (v_dir_n_5 - make_float3 (dot_0(v_dir_n_5, dir_n_5)) * dir_n_5) * make_float3 (inorm_1);
            break;
        }
        float fTmp0C_2 = -2.28522896766662598f * z2_2 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        float fC2_2 = x_3 * fC1_2 - y_3 * fS1_2;
        float fS2_2 = x_3 * fS1_2 + y_3 * fC1_2;
        float pSH12_2 = z_2 * (1.86588168144226074f * z2_2 - 1.11952900886535645f);
        float pSH13_1 = fTmp0C_2 * x_3;
        float pSH11_1 = fTmp0C_2 * y_3;
        float pSH14_1 = fTmp1B_2 * fC1_2;
        float pSH10_1 = fTmp1B_2 * fS1_2;
        float pSH15_1 = -0.59004360437393188f * fC2_2;
        float pSH9_1 = -0.59004360437393188f * fS2_2;
        float3  * _S52 = v_coeffs_1 + int(8);
        *_S52 = *_S52 + make_float3 (pSH9_1) * v_colors_1;
        float3  * _S53 = v_coeffs_1 + int(9);
        *_S53 = *_S53 + make_float3 (pSH10_1) * v_colors_1;
        float3  * _S54 = v_coeffs_1 + int(10);
        *_S54 = *_S54 + make_float3 (pSH11_1) * v_colors_1;
        float3  * _S55 = v_coeffs_1 + int(11);
        *_S55 = *_S55 + make_float3 (pSH12_2) * v_colors_1;
        float3  * _S56 = v_coeffs_1 + int(12);
        *_S56 = *_S56 + make_float3 (pSH13_1) * v_colors_1;
        float3  * _S57 = v_coeffs_1 + int(13);
        *_S57 = *_S57 + make_float3 (pSH14_1) * v_colors_1;
        float3  * _S58 = v_coeffs_1 + int(14);
        *_S58 = *_S58 + make_float3 (pSH15_1) * v_colors_1;
        float fTmp0C_z_1 = -4.57045793533325195f * z_2;
        float _S59 = x_3 * _S42;
        float fC2_x_1 = fC1_2 + _S59 - y_3 * fS1_x_1;
        float _S60 = y_3 * _S42;
        float fC2_y_1 = x_3 * fC1_y_1 - fS1_2 - _S60;
        float fS2_x_1 = fS1_2 + x_3 * fS1_x_1 + _S60;
        float fS2_y_1 = _S59 + fC1_2 + y_3 * fC1_y_1;
        float pSH12_z_1 = 5.59764480590820312f * z2_2 - 1.11952900886535645f;
        float pSH14_x_1 = fTmp1B_2 * _S42;
        float3  * _S61 = coeffs_2 + int(8);
        float3  * _S62 = coeffs_2 + int(14);
        float3  * _S63 = coeffs_2 + int(9);
        float3  * _S64 = coeffs_2 + int(13);
        float3  * _S65 = coeffs_2 + int(12);
        float v_x_3 = v_x_2 + dot_0(v_colors_1, make_float3 (-0.59004360437393188f * fS2_x_1) * *_S61 + make_float3 (-0.59004360437393188f * fC2_x_1) * *_S62 + make_float3 (fTmp1B_2 * fS1_x_1) * *_S63 + make_float3 (pSH14_x_1) * *_S64 + make_float3 (fTmp0C_2) * *_S65);
        float3  * _S66 = coeffs_2 + int(10);
        float v_y_3 = v_y_2 + dot_0(v_colors_1, make_float3 (-0.59004360437393188f * fS2_y_1) * *_S61 + make_float3 (-0.59004360437393188f * fC2_y_1) * *_S62 + make_float3 (pSH14_x_1) * *_S63 + make_float3 (fTmp1B_2 * fC1_y_1) * *_S64 + make_float3 (fTmp0C_2) * *_S66);
        float v_z_3 = v_z_2 + dot_0(v_colors_1, make_float3 (pSH12_z_1) * *(coeffs_2 + int(11)) + make_float3 (fTmp0C_z_1 * x_3) * *_S65 + make_float3 (fTmp0C_z_1 * y_3) * *_S66 + make_float3 (1.44530570507049561f * fC1_2) * *_S64 + make_float3 (1.44530570507049561f * fS1_2) * *_S63);
        if(degree_2 <= int(3))
        {
            float3  dir_n_6 = make_float3 (x_3, y_3, z_2);
            float3  v_dir_n_6 = make_float3 (v_x_3, v_y_3, v_z_3);
            *v_dir_1 = *v_dir_1 + (v_dir_n_6 - make_float3 (dot_0(v_dir_n_6, dir_n_6)) * dir_n_6) * make_float3 (inorm_1);
            break;
        }
        float fTmp0D_2 = z_2 * (-4.68332576751708984f * z2_2 + 2.00713968276977539f);
        float fTmp1C_2 = 3.31161141395568848f * z2_2 - 0.47308734059333801f;
        float fTmp2B_2 = -1.77013075351715088f * z_2;
        float pSH20_1 = 1.9843134880065918f * z_2 * pSH12_2 + -1.00623059272766113f * pSH6_2;
        float pSH21_1 = fTmp0D_2 * x_3;
        float pSH19_1 = fTmp0D_2 * y_3;
        float pSH22_1 = fTmp1C_2 * fC1_2;
        float pSH18_1 = fTmp1C_2 * fS1_2;
        float pSH23_1 = fTmp2B_2 * fC2_2;
        float pSH17_1 = fTmp2B_2 * fS2_2;
        float pSH24_1 = 0.62583571672439575f * (x_3 * fC2_2 - y_3 * fS2_2);
        float pSH16_1 = 0.62583571672439575f * (x_3 * fS2_2 + y_3 * fC2_2);
        float3  * _S67 = v_coeffs_1 + int(15);
        *_S67 = *_S67 + make_float3 (pSH16_1) * v_colors_1;
        float3  * _S68 = v_coeffs_1 + int(16);
        *_S68 = *_S68 + make_float3 (pSH17_1) * v_colors_1;
        float3  * _S69 = v_coeffs_1 + int(17);
        *_S69 = *_S69 + make_float3 (pSH18_1) * v_colors_1;
        float3  * _S70 = v_coeffs_1 + int(18);
        *_S70 = *_S70 + make_float3 (pSH19_1) * v_colors_1;
        float3  * _S71 = v_coeffs_1 + int(19);
        *_S71 = *_S71 + make_float3 (pSH20_1) * v_colors_1;
        float3  * _S72 = v_coeffs_1 + int(20);
        *_S72 = *_S72 + make_float3 (pSH21_1) * v_colors_1;
        float3  * _S73 = v_coeffs_1 + int(21);
        *_S73 = *_S73 + make_float3 (pSH22_1) * v_colors_1;
        float3  * _S74 = v_coeffs_1 + int(22);
        *_S74 = *_S74 + make_float3 (pSH23_1) * v_colors_1;
        float3  * _S75 = v_coeffs_1 + int(23);
        *_S75 = *_S75 + make_float3 (pSH24_1) * v_colors_1;
        float fTmp0D_z_1 = -14.04997730255126953f * z2_2 + 2.00713968276977539f;
        float fTmp1C_z_1 = 6.62322282791137695f * z_2;
        float pSH22_x_1 = fTmp1C_2 * _S42;
        float3  * _S76 = coeffs_2 + int(15);
        float3  * _S77 = coeffs_2 + int(23);
        float3  * _S78 = coeffs_2 + int(16);
        float3  * _S79 = coeffs_2 + int(22);
        float3  * _S80 = coeffs_2 + int(17);
        float3  * _S81 = coeffs_2 + int(21);
        float3  * _S82 = coeffs_2 + int(20);
        float3  * _S83 = coeffs_2 + int(18);
        float3  dir_n_7 = make_float3 (x_3, y_3, z_2);
        float3  v_dir_n_7 = make_float3 (v_x_3 + dot_0(v_colors_1, make_float3 (0.62583571672439575f * (fS2_2 + y_3 * fC2_x_1 + x_3 * fS2_x_1)) * *_S76 + make_float3 (0.62583571672439575f * (fC2_2 + x_3 * fC2_x_1 - y_3 * fS2_x_1)) * *_S77 + make_float3 (fTmp2B_2 * fS2_x_1) * *_S78 + make_float3 (fTmp2B_2 * fC2_x_1) * *_S79 + make_float3 (fTmp1C_2 * fS1_x_1) * *_S80 + make_float3 (pSH22_x_1) * *_S81 + make_float3 (fTmp0D_2) * *_S82), v_y_3 + dot_0(v_colors_1, make_float3 (0.62583571672439575f * (x_3 * fS2_y_1 + fC2_2 + y_3 * fC2_y_1)) * *_S76 + make_float3 (0.62583571672439575f * (x_3 * fC2_y_1 - fS2_2 - y_3 * fS2_y_1)) * *_S77 + make_float3 (fTmp2B_2 * fS2_y_1) * *_S78 + make_float3 (fTmp2B_2 * fC2_y_1) * *_S79 + make_float3 (pSH22_x_1) * *_S80 + make_float3 (fTmp1C_2 * fC1_y_1) * *_S81 + make_float3 (fTmp0D_2) * *_S83), v_z_3 + dot_0(v_colors_1, make_float3 (1.9843134880065918f * (pSH12_2 + z_2 * pSH12_z_1) + -1.00623059272766113f * pSH6_z_1) * *(coeffs_2 + int(19)) + make_float3 (fTmp0D_z_1 * x_3) * *_S82 + make_float3 (fTmp0D_z_1 * y_3) * *_S83 + make_float3 (fTmp1C_z_1 * fC1_2) * *_S81 + make_float3 (fTmp1C_z_1 * fS1_2) * *_S80 + make_float3 (-1.77013075351715088f * fC2_2) * *_S79 + make_float3 (-1.77013075351715088f * fS2_2) * *_S78));
        *v_dir_1 = *v_dir_1 + (v_dir_n_7 - make_float3 (dot_0(v_dir_n_7, dir_n_7)) * dir_n_7) * make_float3 (inorm_1);
        break;
    }
    return;
}

inline __device__ void sh_coeffs_to_color_vjp_atomic(int degree_3, float3  dir_2, float3  coeff_dc_3, float3  * coeffs_3, float3  v_colors_2, float3  * v_coeff_dc_2, float3  * v_coeffs_2, float3  * v_dir_2)
{
    for(;;)
    {
        *v_coeff_dc_2 = *v_coeff_dc_2 + make_float3 (0.282094806432724f) * v_colors_2;
        if(degree_3 <= int(0))
        {
            break;
        }
        float _S84 = dir_2.x;
        float _S85 = dir_2.y;
        float _S86 = dir_2.z;
        float inorm_2 = (F32_rsqrt((_S84 * _S84 + _S85 * _S85 + _S86 * _S86)));
        float x_4 = _S84 * inorm_2;
        float y_4 = _S85 * inorm_2;
        float z_3 = _S86 * inorm_2;
        float3  temp_0 = make_float3 (-0.48860251903533936f * y_4) * v_colors_2;
        float _S87 = dot_0(temp_0, temp_0);
        bool _S88;
        if((F32_isfinite((_S87))))
        {
            _S88 = _S87 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S89 = v_coeffs_2 + int(0);
            float _S90 = atomicAdd(&(_S89->x), temp_0.x);
            float _S91 = atomicAdd(&(_S89->y), temp_0.y);
            float _S92 = atomicAdd(&(_S89->z), temp_0.z);
        }
        float3  temp_1 = make_float3 (0.48860251903533936f * z_3) * v_colors_2;
        float _S93 = dot_0(temp_1, temp_1);
        if((F32_isfinite((_S93))))
        {
            _S88 = _S93 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S94 = v_coeffs_2 + int(1);
            float _S95 = atomicAdd(&(_S94->x), temp_1.x);
            float _S96 = atomicAdd(&(_S94->y), temp_1.y);
            float _S97 = atomicAdd(&(_S94->z), temp_1.z);
        }
        float3  temp_2 = make_float3 (-0.48860251903533936f * x_4) * v_colors_2;
        float _S98 = dot_0(temp_2, temp_2);
        if((F32_isfinite((_S98))))
        {
            _S88 = _S98 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S99 = v_coeffs_2 + int(2);
            float _S100 = atomicAdd(&(_S99->x), temp_2.x);
            float _S101 = atomicAdd(&(_S99->y), temp_2.y);
            float _S102 = atomicAdd(&(_S99->z), temp_2.z);
        }
        float _S103 = -0.48860251903533936f * dot_0(*(coeffs_3 + int(2)), v_colors_2);
        float _S104 = -0.48860251903533936f * dot_0(*(coeffs_3 + int(0)), v_colors_2);
        float _S105 = 0.48860251903533936f * dot_0(*(coeffs_3 + int(1)), v_colors_2);
        if(degree_3 <= int(1))
        {
            float3  dir_n_8 = make_float3 (x_4, y_4, z_3);
            float3  v_dir_n_8 = make_float3 (_S103, _S104, _S105);
            *v_dir_2 = *v_dir_2 + (v_dir_n_8 - make_float3 (dot_0(v_dir_n_8, dir_n_8)) * dir_n_8) * make_float3 (inorm_2);
            break;
        }
        float z2_3 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_4 * x_4 - y_4 * y_4;
        float _S106 = 2.0f * x_4;
        float fS1_3 = _S106 * y_4;
        float pSH6_3 = 0.94617468118667603f * z2_3 - 0.31539157032966614f;
        float pSH7_2 = fTmp0B_3 * x_4;
        float pSH5_2 = fTmp0B_3 * y_4;
        float pSH8_2 = 0.54627424478530884f * fC1_3;
        float pSH4_2 = 0.54627424478530884f * fS1_3;
        float3  temp_3 = make_float3 (pSH4_2) * v_colors_2;
        float _S107 = dot_0(temp_3, temp_3);
        if((F32_isfinite((_S107))))
        {
            _S88 = _S107 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S108 = v_coeffs_2 + int(3);
            float _S109 = atomicAdd(&(_S108->x), temp_3.x);
            float _S110 = atomicAdd(&(_S108->y), temp_3.y);
            float _S111 = atomicAdd(&(_S108->z), temp_3.z);
        }
        float3  temp_4 = make_float3 (pSH5_2) * v_colors_2;
        float _S112 = dot_0(temp_4, temp_4);
        if((F32_isfinite((_S112))))
        {
            _S88 = _S112 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S113 = v_coeffs_2 + int(4);
            float _S114 = atomicAdd(&(_S113->x), temp_4.x);
            float _S115 = atomicAdd(&(_S113->y), temp_4.y);
            float _S116 = atomicAdd(&(_S113->z), temp_4.z);
        }
        float3  temp_5 = make_float3 (pSH6_3) * v_colors_2;
        float _S117 = dot_0(temp_5, temp_5);
        if((F32_isfinite((_S117))))
        {
            _S88 = _S117 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S118 = v_coeffs_2 + int(5);
            float _S119 = atomicAdd(&(_S118->x), temp_5.x);
            float _S120 = atomicAdd(&(_S118->y), temp_5.y);
            float _S121 = atomicAdd(&(_S118->z), temp_5.z);
        }
        float3  temp_6 = make_float3 (pSH7_2) * v_colors_2;
        float _S122 = dot_0(temp_6, temp_6);
        if((F32_isfinite((_S122))))
        {
            _S88 = _S122 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S123 = v_coeffs_2 + int(6);
            float _S124 = atomicAdd(&(_S123->x), temp_6.x);
            float _S125 = atomicAdd(&(_S123->y), temp_6.y);
            float _S126 = atomicAdd(&(_S123->z), temp_6.z);
        }
        float3  temp_7 = make_float3 (pSH8_2) * v_colors_2;
        float _S127 = dot_0(temp_7, temp_7);
        if((F32_isfinite((_S127))))
        {
            _S88 = _S127 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S128 = v_coeffs_2 + int(7);
            float _S129 = atomicAdd(&(_S128->x), temp_7.x);
            float _S130 = atomicAdd(&(_S128->y), temp_7.y);
            float _S131 = atomicAdd(&(_S128->z), temp_7.z);
        }
        float fC1_y_2 = -2.0f * y_4;
        float fS1_x_2 = 2.0f * y_4;
        float pSH6_z_2 = 1.89234936237335205f * z_3;
        float pSH8_x_2 = 0.54627424478530884f * _S106;
        float3  * _S132 = coeffs_3 + int(3);
        float3  * _S133 = coeffs_3 + int(7);
        float3  * _S134 = coeffs_3 + int(6);
        float v_x_4 = _S103 + dot_0(v_colors_2, make_float3 (0.54627424478530884f * fS1_x_2) * *_S132 + make_float3 (pSH8_x_2) * *_S133 + make_float3 (fTmp0B_3) * *_S134);
        float3  * _S135 = coeffs_3 + int(4);
        float v_y_4 = _S104 + dot_0(v_colors_2, make_float3 (pSH8_x_2) * *_S132 + make_float3 (0.54627424478530884f * fC1_y_2) * *_S133 + make_float3 (fTmp0B_3) * *_S135);
        float v_z_4 = _S105 + dot_0(v_colors_2, make_float3 (pSH6_z_2) * *(coeffs_3 + int(5)) + make_float3 (-1.09254848957061768f * x_4) * *_S134 + make_float3 (-1.09254848957061768f * y_4) * *_S135);
        if(degree_3 <= int(2))
        {
            float3  dir_n_9 = make_float3 (x_4, y_4, z_3);
            float3  v_dir_n_9 = make_float3 (v_x_4, v_y_4, v_z_4);
            *v_dir_2 = *v_dir_2 + (v_dir_n_9 - make_float3 (dot_0(v_dir_n_9, dir_n_9)) * dir_n_9) * make_float3 (inorm_2);
            break;
        }
        float fTmp0C_3 = -2.28522896766662598f * z2_3 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        float fC2_3 = x_4 * fC1_3 - y_4 * fS1_3;
        float fS2_3 = x_4 * fS1_3 + y_4 * fC1_3;
        float pSH12_3 = z_3 * (1.86588168144226074f * z2_3 - 1.11952900886535645f);
        float pSH13_2 = fTmp0C_3 * x_4;
        float pSH11_2 = fTmp0C_3 * y_4;
        float pSH14_2 = fTmp1B_3 * fC1_3;
        float pSH10_2 = fTmp1B_3 * fS1_3;
        float pSH15_2 = -0.59004360437393188f * fC2_3;
        float pSH9_2 = -0.59004360437393188f * fS2_3;
        float3  temp_8 = make_float3 (pSH9_2) * v_colors_2;
        float _S136 = dot_0(temp_8, temp_8);
        if((F32_isfinite((_S136))))
        {
            _S88 = _S136 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S137 = v_coeffs_2 + int(8);
            float _S138 = atomicAdd(&(_S137->x), temp_8.x);
            float _S139 = atomicAdd(&(_S137->y), temp_8.y);
            float _S140 = atomicAdd(&(_S137->z), temp_8.z);
        }
        float3  temp_9 = make_float3 (pSH10_2) * v_colors_2;
        float _S141 = dot_0(temp_9, temp_9);
        if((F32_isfinite((_S141))))
        {
            _S88 = _S141 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S142 = v_coeffs_2 + int(9);
            float _S143 = atomicAdd(&(_S142->x), temp_9.x);
            float _S144 = atomicAdd(&(_S142->y), temp_9.y);
            float _S145 = atomicAdd(&(_S142->z), temp_9.z);
        }
        float3  temp_10 = make_float3 (pSH11_2) * v_colors_2;
        float _S146 = dot_0(temp_10, temp_10);
        if((F32_isfinite((_S146))))
        {
            _S88 = _S146 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S147 = v_coeffs_2 + int(10);
            float _S148 = atomicAdd(&(_S147->x), temp_10.x);
            float _S149 = atomicAdd(&(_S147->y), temp_10.y);
            float _S150 = atomicAdd(&(_S147->z), temp_10.z);
        }
        float3  temp_11 = make_float3 (pSH12_3) * v_colors_2;
        float _S151 = dot_0(temp_11, temp_11);
        if((F32_isfinite((_S151))))
        {
            _S88 = _S151 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S152 = v_coeffs_2 + int(11);
            float _S153 = atomicAdd(&(_S152->x), temp_11.x);
            float _S154 = atomicAdd(&(_S152->y), temp_11.y);
            float _S155 = atomicAdd(&(_S152->z), temp_11.z);
        }
        float3  temp_12 = make_float3 (pSH13_2) * v_colors_2;
        float _S156 = dot_0(temp_12, temp_12);
        if((F32_isfinite((_S156))))
        {
            _S88 = _S156 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S157 = v_coeffs_2 + int(12);
            float _S158 = atomicAdd(&(_S157->x), temp_12.x);
            float _S159 = atomicAdd(&(_S157->y), temp_12.y);
            float _S160 = atomicAdd(&(_S157->z), temp_12.z);
        }
        float3  temp_13 = make_float3 (pSH14_2) * v_colors_2;
        float _S161 = dot_0(temp_13, temp_13);
        if((F32_isfinite((_S161))))
        {
            _S88 = _S161 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S162 = v_coeffs_2 + int(13);
            float _S163 = atomicAdd(&(_S162->x), temp_13.x);
            float _S164 = atomicAdd(&(_S162->y), temp_13.y);
            float _S165 = atomicAdd(&(_S162->z), temp_13.z);
        }
        float3  temp_14 = make_float3 (pSH15_2) * v_colors_2;
        float _S166 = dot_0(temp_14, temp_14);
        if((F32_isfinite((_S166))))
        {
            _S88 = _S166 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S167 = v_coeffs_2 + int(14);
            float _S168 = atomicAdd(&(_S167->x), temp_14.x);
            float _S169 = atomicAdd(&(_S167->y), temp_14.y);
            float _S170 = atomicAdd(&(_S167->z), temp_14.z);
        }
        float fTmp0C_z_2 = -4.57045793533325195f * z_3;
        float _S171 = x_4 * _S106;
        float fC2_x_2 = fC1_3 + _S171 - y_4 * fS1_x_2;
        float _S172 = y_4 * _S106;
        float fC2_y_2 = x_4 * fC1_y_2 - fS1_3 - _S172;
        float fS2_x_2 = fS1_3 + x_4 * fS1_x_2 + _S172;
        float fS2_y_2 = _S171 + fC1_3 + y_4 * fC1_y_2;
        float pSH12_z_2 = 5.59764480590820312f * z2_3 - 1.11952900886535645f;
        float pSH14_x_2 = fTmp1B_3 * _S106;
        float3  * _S173 = coeffs_3 + int(8);
        float3  * _S174 = coeffs_3 + int(14);
        float3  * _S175 = coeffs_3 + int(9);
        float3  * _S176 = coeffs_3 + int(13);
        float3  * _S177 = coeffs_3 + int(12);
        float v_x_5 = v_x_4 + dot_0(v_colors_2, make_float3 (-0.59004360437393188f * fS2_x_2) * *_S173 + make_float3 (-0.59004360437393188f * fC2_x_2) * *_S174 + make_float3 (fTmp1B_3 * fS1_x_2) * *_S175 + make_float3 (pSH14_x_2) * *_S176 + make_float3 (fTmp0C_3) * *_S177);
        float3  * _S178 = coeffs_3 + int(10);
        float v_y_5 = v_y_4 + dot_0(v_colors_2, make_float3 (-0.59004360437393188f * fS2_y_2) * *_S173 + make_float3 (-0.59004360437393188f * fC2_y_2) * *_S174 + make_float3 (pSH14_x_2) * *_S175 + make_float3 (fTmp1B_3 * fC1_y_2) * *_S176 + make_float3 (fTmp0C_3) * *_S178);
        float v_z_5 = v_z_4 + dot_0(v_colors_2, make_float3 (pSH12_z_2) * *(coeffs_3 + int(11)) + make_float3 (fTmp0C_z_2 * x_4) * *_S177 + make_float3 (fTmp0C_z_2 * y_4) * *_S178 + make_float3 (1.44530570507049561f * fC1_3) * *_S176 + make_float3 (1.44530570507049561f * fS1_3) * *_S175);
        if(degree_3 <= int(3))
        {
            float3  dir_n_10 = make_float3 (x_4, y_4, z_3);
            float3  v_dir_n_10 = make_float3 (v_x_5, v_y_5, v_z_5);
            *v_dir_2 = *v_dir_2 + (v_dir_n_10 - make_float3 (dot_0(v_dir_n_10, dir_n_10)) * dir_n_10) * make_float3 (inorm_2);
            break;
        }
        float fTmp0D_3 = z_3 * (-4.68332576751708984f * z2_3 + 2.00713968276977539f);
        float fTmp1C_3 = 3.31161141395568848f * z2_3 - 0.47308734059333801f;
        float fTmp2B_3 = -1.77013075351715088f * z_3;
        float pSH20_2 = 1.9843134880065918f * z_3 * pSH12_3 + -1.00623059272766113f * pSH6_3;
        float pSH21_2 = fTmp0D_3 * x_4;
        float pSH19_2 = fTmp0D_3 * y_4;
        float pSH22_2 = fTmp1C_3 * fC1_3;
        float pSH18_2 = fTmp1C_3 * fS1_3;
        float pSH23_2 = fTmp2B_3 * fC2_3;
        float pSH17_2 = fTmp2B_3 * fS2_3;
        float pSH24_2 = 0.62583571672439575f * (x_4 * fC2_3 - y_4 * fS2_3);
        float pSH16_2 = 0.62583571672439575f * (x_4 * fS2_3 + y_4 * fC2_3);
        float3  temp_15 = make_float3 (pSH16_2) * v_colors_2;
        float _S179 = dot_0(temp_15, temp_15);
        if((F32_isfinite((_S179))))
        {
            _S88 = _S179 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S180 = v_coeffs_2 + int(15);
            float _S181 = atomicAdd(&(_S180->x), temp_15.x);
            float _S182 = atomicAdd(&(_S180->y), temp_15.y);
            float _S183 = atomicAdd(&(_S180->z), temp_15.z);
        }
        float3  temp_16 = make_float3 (pSH17_2) * v_colors_2;
        float _S184 = dot_0(temp_16, temp_16);
        if((F32_isfinite((_S184))))
        {
            _S88 = _S184 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S185 = v_coeffs_2 + int(16);
            float _S186 = atomicAdd(&(_S185->x), temp_16.x);
            float _S187 = atomicAdd(&(_S185->y), temp_16.y);
            float _S188 = atomicAdd(&(_S185->z), temp_16.z);
        }
        float3  temp_17 = make_float3 (pSH18_2) * v_colors_2;
        float _S189 = dot_0(temp_17, temp_17);
        if((F32_isfinite((_S189))))
        {
            _S88 = _S189 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S190 = v_coeffs_2 + int(17);
            float _S191 = atomicAdd(&(_S190->x), temp_17.x);
            float _S192 = atomicAdd(&(_S190->y), temp_17.y);
            float _S193 = atomicAdd(&(_S190->z), temp_17.z);
        }
        float3  temp_18 = make_float3 (pSH19_2) * v_colors_2;
        float _S194 = dot_0(temp_18, temp_18);
        if((F32_isfinite((_S194))))
        {
            _S88 = _S194 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S195 = v_coeffs_2 + int(18);
            float _S196 = atomicAdd(&(_S195->x), temp_18.x);
            float _S197 = atomicAdd(&(_S195->y), temp_18.y);
            float _S198 = atomicAdd(&(_S195->z), temp_18.z);
        }
        float3  temp_19 = make_float3 (pSH20_2) * v_colors_2;
        float _S199 = dot_0(temp_19, temp_19);
        if((F32_isfinite((_S199))))
        {
            _S88 = _S199 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S200 = v_coeffs_2 + int(19);
            float _S201 = atomicAdd(&(_S200->x), temp_19.x);
            float _S202 = atomicAdd(&(_S200->y), temp_19.y);
            float _S203 = atomicAdd(&(_S200->z), temp_19.z);
        }
        float3  temp_20 = make_float3 (pSH21_2) * v_colors_2;
        float _S204 = dot_0(temp_20, temp_20);
        if((F32_isfinite((_S204))))
        {
            _S88 = _S204 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S205 = v_coeffs_2 + int(20);
            float _S206 = atomicAdd(&(_S205->x), temp_20.x);
            float _S207 = atomicAdd(&(_S205->y), temp_20.y);
            float _S208 = atomicAdd(&(_S205->z), temp_20.z);
        }
        float3  temp_21 = make_float3 (pSH22_2) * v_colors_2;
        float _S209 = dot_0(temp_21, temp_21);
        if((F32_isfinite((_S209))))
        {
            _S88 = _S209 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S210 = v_coeffs_2 + int(21);
            float _S211 = atomicAdd(&(_S210->x), temp_21.x);
            float _S212 = atomicAdd(&(_S210->y), temp_21.y);
            float _S213 = atomicAdd(&(_S210->z), temp_21.z);
        }
        float3  temp_22 = make_float3 (pSH23_2) * v_colors_2;
        float _S214 = dot_0(temp_22, temp_22);
        if((F32_isfinite((_S214))))
        {
            _S88 = _S214 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S215 = v_coeffs_2 + int(22);
            float _S216 = atomicAdd(&(_S215->x), temp_22.x);
            float _S217 = atomicAdd(&(_S215->y), temp_22.y);
            float _S218 = atomicAdd(&(_S215->z), temp_22.z);
        }
        float3  temp_23 = make_float3 (pSH24_2) * v_colors_2;
        float _S219 = dot_0(temp_23, temp_23);
        if((F32_isfinite((_S219))))
        {
            _S88 = _S219 != 0.0f;
        }
        else
        {
            _S88 = false;
        }
        if(_S88)
        {
            float3  * _S220 = v_coeffs_2 + int(23);
            float _S221 = atomicAdd(&(_S220->x), temp_23.x);
            float _S222 = atomicAdd(&(_S220->y), temp_23.y);
            float _S223 = atomicAdd(&(_S220->z), temp_23.z);
        }
        float fTmp0D_z_2 = -14.04997730255126953f * z2_3 + 2.00713968276977539f;
        float fTmp1C_z_2 = 6.62322282791137695f * z_3;
        float pSH22_x_2 = fTmp1C_3 * _S106;
        float3  * _S224 = coeffs_3 + int(15);
        float3  * _S225 = coeffs_3 + int(23);
        float3  * _S226 = coeffs_3 + int(16);
        float3  * _S227 = coeffs_3 + int(22);
        float3  * _S228 = coeffs_3 + int(17);
        float3  * _S229 = coeffs_3 + int(21);
        float3  * _S230 = coeffs_3 + int(20);
        float3  * _S231 = coeffs_3 + int(18);
        float3  dir_n_11 = make_float3 (x_4, y_4, z_3);
        float3  v_dir_n_11 = make_float3 (v_x_5 + dot_0(v_colors_2, make_float3 (0.62583571672439575f * (fS2_3 + y_4 * fC2_x_2 + x_4 * fS2_x_2)) * *_S224 + make_float3 (0.62583571672439575f * (fC2_3 + x_4 * fC2_x_2 - y_4 * fS2_x_2)) * *_S225 + make_float3 (fTmp2B_3 * fS2_x_2) * *_S226 + make_float3 (fTmp2B_3 * fC2_x_2) * *_S227 + make_float3 (fTmp1C_3 * fS1_x_2) * *_S228 + make_float3 (pSH22_x_2) * *_S229 + make_float3 (fTmp0D_3) * *_S230), v_y_5 + dot_0(v_colors_2, make_float3 (0.62583571672439575f * (x_4 * fS2_y_2 + fC2_3 + y_4 * fC2_y_2)) * *_S224 + make_float3 (0.62583571672439575f * (x_4 * fC2_y_2 - fS2_3 - y_4 * fS2_y_2)) * *_S225 + make_float3 (fTmp2B_3 * fS2_y_2) * *_S226 + make_float3 (fTmp2B_3 * fC2_y_2) * *_S227 + make_float3 (pSH22_x_2) * *_S228 + make_float3 (fTmp1C_3 * fC1_y_2) * *_S229 + make_float3 (fTmp0D_3) * *_S231), v_z_5 + dot_0(v_colors_2, make_float3 (1.9843134880065918f * (pSH12_3 + z_3 * pSH12_z_2) + -1.00623059272766113f * pSH6_z_2) * *(coeffs_3 + int(19)) + make_float3 (fTmp0D_z_2 * x_4) * *_S230 + make_float3 (fTmp0D_z_2 * y_4) * *_S231 + make_float3 (fTmp1C_z_2 * fC1_3) * *_S229 + make_float3 (fTmp1C_z_2 * fS1_3) * *_S228 + make_float3 (-1.77013075351715088f * fC2_3) * *_S227 + make_float3 (-1.77013075351715088f * fS2_3) * *_S226));
        *v_dir_2 = *v_dir_2 + (v_dir_n_11 - make_float3 (dot_0(v_dir_n_11, dir_n_11)) * dir_n_11) * make_float3 (inorm_2);
        break;
    }
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_5)
{
    Matrix<float, 3, 3>  result_2;
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
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_0)), c_0) = _slang_vector_get_element(x_5.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_2;
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_0)
{
    float _S232 = (*left_0).primal_0.rows[int(0)].x * dOut_0.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_0.x;
    float sum_0 = _S232 + (*left_0).primal_0.rows[int(1)].x * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_0.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_0.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S233 = (*left_0).primal_0.rows[int(0)].y * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_0.x;
    float sum_2 = _S233 + (*left_0).primal_0.rows[int(1)].y * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_0.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_0.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S234 = (*left_0).primal_0.rows[int(0)].z * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_0.x;
    float sum_4 = _S234 + (*left_0).primal_0.rows[int(1)].z * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_0.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_0.z;
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
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_1], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_3, i_1) = sum_6;
        i_1 = i_1 + int(1);
    }
    return result_3;
}

inline __device__ float3  max_0(float3  x_6, float3  y_5)
{
    float3  result_4;
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
        *_slang_vector_get_element_ptr(&result_4, i_2) = (F32_max((_slang_vector_get_element(x_6, i_2)), (_slang_vector_get_element(y_5, i_2))));
        i_2 = i_2 + int(1);
    }
    return result_4;
}

inline __device__ float3  sh_coeffs_to_color(int degree_4, float3  mean_0, Matrix<float, 3, 3>  R_0, float3  t_0, float3  coeff_dc_4, float3  * coeffs_4)
{
    float3  _S235;
    float3  _S236 = mean_0 + mul_0(transpose_0(R_0), t_0);
    for(;;)
    {
        float3  colors_4 = make_float3 (0.282094806432724f) * coeff_dc_4;
        if(degree_4 <= int(0))
        {
            _S235 = colors_4;
            break;
        }
        float _S237 = _S236.x;
        float _S238 = _S236.y;
        float _S239 = _S236.z;
        float inv_norm_1 = (F32_rsqrt((_S237 * _S237 + _S238 * _S238 + _S239 * _S239)));
        float x_7 = _S237 * inv_norm_1;
        float y_6 = _S238 * inv_norm_1;
        float z_4 = _S239 * inv_norm_1;
        float3  colors_5 = colors_4 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * *(coeffs_4 + int(0)) + make_float3 (z_4) * *(coeffs_4 + int(1)) - make_float3 (x_7) * *(coeffs_4 + int(2)));
        if(degree_4 <= int(1))
        {
            _S235 = colors_5;
            break;
        }
        float z2_4 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_7 * x_7 - y_6 * y_6;
        float fS1_4 = 2.0f * x_7 * y_6;
        float pSH6_4 = 0.94617468118667603f * z2_4 - 0.31539157032966614f;
        float3  colors_6 = colors_5 + (make_float3 (0.54627424478530884f * fS1_4) * *(coeffs_4 + int(3)) + make_float3 (fTmp0B_4 * y_6) * *(coeffs_4 + int(4)) + make_float3 (pSH6_4) * *(coeffs_4 + int(5)) + make_float3 (fTmp0B_4 * x_7) * *(coeffs_4 + int(6)) + make_float3 (0.54627424478530884f * fC1_4) * *(coeffs_4 + int(7)));
        if(degree_4 <= int(2))
        {
            _S235 = colors_6;
            break;
        }
        float fTmp0C_4 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        float fC2_4 = x_7 * fC1_4 - y_6 * fS1_4;
        float fS2_4 = x_7 * fS1_4 + y_6 * fC1_4;
        float pSH12_4 = z_4 * (1.86588168144226074f * z2_4 - 1.11952900886535645f);
        float3  colors_7 = colors_6 + (make_float3 (-0.59004360437393188f * fS2_4) * *(coeffs_4 + int(8)) + make_float3 (fTmp1B_4 * fS1_4) * *(coeffs_4 + int(9)) + make_float3 (fTmp0C_4 * y_6) * *(coeffs_4 + int(10)) + make_float3 (pSH12_4) * *(coeffs_4 + int(11)) + make_float3 (fTmp0C_4 * x_7) * *(coeffs_4 + int(12)) + make_float3 (fTmp1B_4 * fC1_4) * *(coeffs_4 + int(13)) + make_float3 (-0.59004360437393188f * fC2_4) * *(coeffs_4 + int(14)));
        if(degree_4 <= int(3))
        {
            _S235 = colors_7;
            break;
        }
        float fTmp0D_4 = z_4 * (-4.68332576751708984f * z2_4 + 2.00713968276977539f);
        float fTmp1C_4 = 3.31161141395568848f * z2_4 - 0.47308734059333801f;
        float fTmp2B_4 = -1.77013075351715088f * z_4;
        _S235 = colors_7 + (make_float3 (0.62583571672439575f * (x_7 * fS2_4 + y_6 * fC2_4)) * *(coeffs_4 + int(15)) + make_float3 (fTmp2B_4 * fS2_4) * *(coeffs_4 + int(16)) + make_float3 (fTmp1C_4 * fS1_4) * *(coeffs_4 + int(17)) + make_float3 (fTmp0D_4 * y_6) * *(coeffs_4 + int(18)) + make_float3 (1.9843134880065918f * z_4 * pSH12_4 - 1.00623059272766113f * pSH6_4) * *(coeffs_4 + int(19)) + make_float3 (fTmp0D_4 * x_7) * *(coeffs_4 + int(20)) + make_float3 (fTmp1C_4 * fC1_4) * *(coeffs_4 + int(21)) + make_float3 (fTmp2B_4 * fC2_4) * *(coeffs_4 + int(22)) + make_float3 (0.62583571672439575f * (x_7 * fC2_4 - y_6 * fS2_4)) * *(coeffs_4 + int(23)));
        break;
    }
    return max_0(_S235 + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S240, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S241, float3  _S242)
{
    _d_mul_0(_S240, _S241, _S242);
    return;
}

inline __device__ void sh_coeffs_to_color_vjp(int degree_5, float3  mean_1, Matrix<float, 3, 3>  R_1, float3  t_1, float3  coeff_dc_5, float3  * coeffs_5, float3  v_colors_3, float3  * v_coeff_dc_3, float3  * v_coeffs_3, float3  * v_mean_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    bool _S243;
    for(;;)
    {
        float3  v_viewdir_0;
        Matrix<float, 3, 3>  _S244 = transpose_0(R_1);
        float3  _S245 = mean_1 + mul_0(_S244, t_1);
        for(;;)
        {
            float3  colors_8 = make_float3 (0.282094806432724f) * coeff_dc_5;
            bool _S246 = degree_5 <= int(0);
            _S243 = _S246;
            if(_S246)
            {
                v_viewdir_0 = colors_8;
                break;
            }
            float _S247 = _S245.x;
            float _S248 = _S245.y;
            float _S249 = _S245.z;
            float inv_norm_2 = (F32_rsqrt((_S247 * _S247 + _S248 * _S248 + _S249 * _S249)));
            float x_8 = _S247 * inv_norm_2;
            float y_7 = _S248 * inv_norm_2;
            float z_5 = _S249 * inv_norm_2;
            float3  colors_9 = colors_8 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * *(coeffs_5 + int(0)) + make_float3 (z_5) * *(coeffs_5 + int(1)) - make_float3 (x_8) * *(coeffs_5 + int(2)));
            if(degree_5 <= int(1))
            {
                v_viewdir_0 = colors_9;
                break;
            }
            float z2_5 = z_5 * z_5;
            float fTmp0B_5 = -1.09254848957061768f * z_5;
            float fC1_5 = x_8 * x_8 - y_7 * y_7;
            float fS1_5 = 2.0f * x_8 * y_7;
            float pSH6_5 = 0.94617468118667603f * z2_5 - 0.31539157032966614f;
            float3  colors_10 = colors_9 + (make_float3 (0.54627424478530884f * fS1_5) * *(coeffs_5 + int(3)) + make_float3 (fTmp0B_5 * y_7) * *(coeffs_5 + int(4)) + make_float3 (pSH6_5) * *(coeffs_5 + int(5)) + make_float3 (fTmp0B_5 * x_8) * *(coeffs_5 + int(6)) + make_float3 (0.54627424478530884f * fC1_5) * *(coeffs_5 + int(7)));
            if(degree_5 <= int(2))
            {
                v_viewdir_0 = colors_10;
                break;
            }
            float fTmp0C_5 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
            float fTmp1B_5 = 1.44530570507049561f * z_5;
            float fC2_5 = x_8 * fC1_5 - y_7 * fS1_5;
            float fS2_5 = x_8 * fS1_5 + y_7 * fC1_5;
            float pSH12_5 = z_5 * (1.86588168144226074f * z2_5 - 1.11952900886535645f);
            float3  colors_11 = colors_10 + (make_float3 (-0.59004360437393188f * fS2_5) * *(coeffs_5 + int(8)) + make_float3 (fTmp1B_5 * fS1_5) * *(coeffs_5 + int(9)) + make_float3 (fTmp0C_5 * y_7) * *(coeffs_5 + int(10)) + make_float3 (pSH12_5) * *(coeffs_5 + int(11)) + make_float3 (fTmp0C_5 * x_8) * *(coeffs_5 + int(12)) + make_float3 (fTmp1B_5 * fC1_5) * *(coeffs_5 + int(13)) + make_float3 (-0.59004360437393188f * fC2_5) * *(coeffs_5 + int(14)));
            if(degree_5 <= int(3))
            {
                v_viewdir_0 = colors_11;
                break;
            }
            float fTmp0D_5 = z_5 * (-4.68332576751708984f * z2_5 + 2.00713968276977539f);
            float fTmp1C_5 = 3.31161141395568848f * z2_5 - 0.47308734059333801f;
            float fTmp2B_5 = -1.77013075351715088f * z_5;
            v_viewdir_0 = colors_11 + (make_float3 (0.62583571672439575f * (x_8 * fS2_5 + y_7 * fC2_5)) * *(coeffs_5 + int(15)) + make_float3 (fTmp2B_5 * fS2_5) * *(coeffs_5 + int(16)) + make_float3 (fTmp1C_5 * fS1_5) * *(coeffs_5 + int(17)) + make_float3 (fTmp0D_5 * y_7) * *(coeffs_5 + int(18)) + make_float3 (1.9843134880065918f * z_5 * pSH12_5 - 1.00623059272766113f * pSH6_5) * *(coeffs_5 + int(19)) + make_float3 (fTmp0D_5 * x_8) * *(coeffs_5 + int(20)) + make_float3 (fTmp1C_5 * fC1_5) * *(coeffs_5 + int(21)) + make_float3 (fTmp2B_5 * fC2_5) * *(coeffs_5 + int(22)) + make_float3 (0.62583571672439575f * (x_8 * fC2_5 - y_7 * fS2_5)) * *(coeffs_5 + int(23)));
            break;
        }
        float3  _S250 = v_colors_3 * make_float3 (float((v_viewdir_0.x) >= -0.5f), float((v_viewdir_0.y) >= -0.5f), float((v_viewdir_0.z) >= -0.5f));
        for(;;)
        {
            *v_coeff_dc_3 = make_float3 (0.282094806432724f) * _S250;
            if(_S243)
            {
                int3  _S251 = make_int3 (int(0));
                float3  _S252 = make_float3 ((float)_S251.x, (float)_S251.y, (float)_S251.z);
                v_viewdir_0 = _S252;
                break;
            }
            float _S253 = _S245.x;
            float _S254 = _S245.y;
            float _S255 = _S245.z;
            float inorm_3 = (F32_rsqrt((_S253 * _S253 + _S254 * _S254 + _S255 * _S255)));
            float x_9 = _S253 * inorm_3;
            float y_8 = _S254 * inorm_3;
            float z_6 = _S255 * inorm_3;
            *(v_coeffs_3 + int(0)) = make_float3 (-0.48860251903533936f * y_8) * _S250;
            *(v_coeffs_3 + int(1)) = make_float3 (0.48860251903533936f * z_6) * _S250;
            *(v_coeffs_3 + int(2)) = make_float3 (-0.48860251903533936f * x_9) * _S250;
            float _S256 = -0.48860251903533936f * dot_0(*(coeffs_5 + int(2)), _S250);
            float _S257 = -0.48860251903533936f * dot_0(*(coeffs_5 + int(0)), _S250);
            float _S258 = 0.48860251903533936f * dot_0(*(coeffs_5 + int(1)), _S250);
            if(degree_5 <= int(1))
            {
                float3  dir_n_12 = make_float3 (x_9, y_8, z_6);
                float3  v_dir_n_12 = make_float3 (_S256, _S257, _S258);
                v_viewdir_0 = (v_dir_n_12 - make_float3 (dot_0(v_dir_n_12, dir_n_12)) * dir_n_12) * make_float3 (inorm_3);
                break;
            }
            float z2_6 = z_6 * z_6;
            float fTmp0B_6 = -1.09254848957061768f * z_6;
            float fC1_6 = x_9 * x_9 - y_8 * y_8;
            float _S259 = 2.0f * x_9;
            float fS1_6 = _S259 * y_8;
            float pSH6_6 = 0.94617468118667603f * z2_6 - 0.31539157032966614f;
            float pSH7_3 = fTmp0B_6 * x_9;
            float pSH5_3 = fTmp0B_6 * y_8;
            float pSH8_3 = 0.54627424478530884f * fC1_6;
            float pSH4_3 = 0.54627424478530884f * fS1_6;
            *(v_coeffs_3 + int(3)) = make_float3 (pSH4_3) * _S250;
            *(v_coeffs_3 + int(4)) = make_float3 (pSH5_3) * _S250;
            *(v_coeffs_3 + int(5)) = make_float3 (pSH6_6) * _S250;
            *(v_coeffs_3 + int(6)) = make_float3 (pSH7_3) * _S250;
            *(v_coeffs_3 + int(7)) = make_float3 (pSH8_3) * _S250;
            float fC1_y_3 = -2.0f * y_8;
            float fS1_x_3 = 2.0f * y_8;
            float pSH6_z_3 = 1.89234936237335205f * z_6;
            float pSH8_x_3 = 0.54627424478530884f * _S259;
            float3  * _S260 = coeffs_5 + int(3);
            float3  * _S261 = coeffs_5 + int(7);
            float3  * _S262 = coeffs_5 + int(6);
            float v_x_6 = _S256 + dot_0(_S250, make_float3 (0.54627424478530884f * fS1_x_3) * *_S260 + make_float3 (pSH8_x_3) * *_S261 + make_float3 (fTmp0B_6) * *_S262);
            float3  * _S263 = coeffs_5 + int(4);
            float v_y_6 = _S257 + dot_0(_S250, make_float3 (pSH8_x_3) * *_S260 + make_float3 (0.54627424478530884f * fC1_y_3) * *_S261 + make_float3 (fTmp0B_6) * *_S263);
            float v_z_6 = _S258 + dot_0(_S250, make_float3 (pSH6_z_3) * *(coeffs_5 + int(5)) + make_float3 (-1.09254848957061768f * x_9) * *_S262 + make_float3 (-1.09254848957061768f * y_8) * *_S263);
            if(degree_5 <= int(2))
            {
                float3  dir_n_13 = make_float3 (x_9, y_8, z_6);
                float3  v_dir_n_13 = make_float3 (v_x_6, v_y_6, v_z_6);
                v_viewdir_0 = (v_dir_n_13 - make_float3 (dot_0(v_dir_n_13, dir_n_13)) * dir_n_13) * make_float3 (inorm_3);
                break;
            }
            float fTmp0C_6 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
            float fTmp1B_6 = 1.44530570507049561f * z_6;
            float fC2_6 = x_9 * fC1_6 - y_8 * fS1_6;
            float fS2_6 = x_9 * fS1_6 + y_8 * fC1_6;
            float pSH12_6 = z_6 * (1.86588168144226074f * z2_6 - 1.11952900886535645f);
            float pSH13_3 = fTmp0C_6 * x_9;
            float pSH11_3 = fTmp0C_6 * y_8;
            float pSH14_3 = fTmp1B_6 * fC1_6;
            float pSH10_3 = fTmp1B_6 * fS1_6;
            float pSH15_3 = -0.59004360437393188f * fC2_6;
            float pSH9_3 = -0.59004360437393188f * fS2_6;
            *(v_coeffs_3 + int(8)) = make_float3 (pSH9_3) * _S250;
            *(v_coeffs_3 + int(9)) = make_float3 (pSH10_3) * _S250;
            *(v_coeffs_3 + int(10)) = make_float3 (pSH11_3) * _S250;
            *(v_coeffs_3 + int(11)) = make_float3 (pSH12_6) * _S250;
            *(v_coeffs_3 + int(12)) = make_float3 (pSH13_3) * _S250;
            *(v_coeffs_3 + int(13)) = make_float3 (pSH14_3) * _S250;
            *(v_coeffs_3 + int(14)) = make_float3 (pSH15_3) * _S250;
            float fTmp0C_z_3 = -4.57045793533325195f * z_6;
            float _S264 = x_9 * _S259;
            float fC2_x_3 = fC1_6 + _S264 - y_8 * fS1_x_3;
            float _S265 = y_8 * _S259;
            float fC2_y_3 = x_9 * fC1_y_3 - fS1_6 - _S265;
            float fS2_x_3 = fS1_6 + x_9 * fS1_x_3 + _S265;
            float fS2_y_3 = _S264 + fC1_6 + y_8 * fC1_y_3;
            float pSH12_z_3 = 5.59764480590820312f * z2_6 - 1.11952900886535645f;
            float pSH14_x_3 = fTmp1B_6 * _S259;
            float3  * _S266 = coeffs_5 + int(8);
            float3  * _S267 = coeffs_5 + int(14);
            float3  * _S268 = coeffs_5 + int(9);
            float3  * _S269 = coeffs_5 + int(13);
            float3  * _S270 = coeffs_5 + int(12);
            float v_x_7 = v_x_6 + dot_0(_S250, make_float3 (-0.59004360437393188f * fS2_x_3) * *_S266 + make_float3 (-0.59004360437393188f * fC2_x_3) * *_S267 + make_float3 (fTmp1B_6 * fS1_x_3) * *_S268 + make_float3 (pSH14_x_3) * *_S269 + make_float3 (fTmp0C_6) * *_S270);
            float3  * _S271 = coeffs_5 + int(10);
            float v_y_7 = v_y_6 + dot_0(_S250, make_float3 (-0.59004360437393188f * fS2_y_3) * *_S266 + make_float3 (-0.59004360437393188f * fC2_y_3) * *_S267 + make_float3 (pSH14_x_3) * *_S268 + make_float3 (fTmp1B_6 * fC1_y_3) * *_S269 + make_float3 (fTmp0C_6) * *_S271);
            float v_z_7 = v_z_6 + dot_0(_S250, make_float3 (pSH12_z_3) * *(coeffs_5 + int(11)) + make_float3 (fTmp0C_z_3 * x_9) * *_S270 + make_float3 (fTmp0C_z_3 * y_8) * *_S271 + make_float3 (1.44530570507049561f * fC1_6) * *_S269 + make_float3 (1.44530570507049561f * fS1_6) * *_S268);
            if(degree_5 <= int(3))
            {
                float3  dir_n_14 = make_float3 (x_9, y_8, z_6);
                float3  v_dir_n_14 = make_float3 (v_x_7, v_y_7, v_z_7);
                v_viewdir_0 = (v_dir_n_14 - make_float3 (dot_0(v_dir_n_14, dir_n_14)) * dir_n_14) * make_float3 (inorm_3);
                break;
            }
            float fTmp0D_6 = z_6 * (-4.68332576751708984f * z2_6 + 2.00713968276977539f);
            float fTmp1C_6 = 3.31161141395568848f * z2_6 - 0.47308734059333801f;
            float fTmp2B_6 = -1.77013075351715088f * z_6;
            float pSH20_3 = 1.9843134880065918f * z_6 * pSH12_6 + -1.00623059272766113f * pSH6_6;
            float pSH21_3 = fTmp0D_6 * x_9;
            float pSH19_3 = fTmp0D_6 * y_8;
            float pSH22_3 = fTmp1C_6 * fC1_6;
            float pSH18_3 = fTmp1C_6 * fS1_6;
            float pSH23_3 = fTmp2B_6 * fC2_6;
            float pSH17_3 = fTmp2B_6 * fS2_6;
            float pSH24_3 = 0.62583571672439575f * (x_9 * fC2_6 - y_8 * fS2_6);
            float pSH16_3 = 0.62583571672439575f * (x_9 * fS2_6 + y_8 * fC2_6);
            *(v_coeffs_3 + int(15)) = make_float3 (pSH16_3) * _S250;
            *(v_coeffs_3 + int(16)) = make_float3 (pSH17_3) * _S250;
            *(v_coeffs_3 + int(17)) = make_float3 (pSH18_3) * _S250;
            *(v_coeffs_3 + int(18)) = make_float3 (pSH19_3) * _S250;
            *(v_coeffs_3 + int(19)) = make_float3 (pSH20_3) * _S250;
            *(v_coeffs_3 + int(20)) = make_float3 (pSH21_3) * _S250;
            *(v_coeffs_3 + int(21)) = make_float3 (pSH22_3) * _S250;
            *(v_coeffs_3 + int(22)) = make_float3 (pSH23_3) * _S250;
            *(v_coeffs_3 + int(23)) = make_float3 (pSH24_3) * _S250;
            float fTmp0D_z_3 = -14.04997730255126953f * z2_6 + 2.00713968276977539f;
            float fTmp1C_z_3 = 6.62322282791137695f * z_6;
            float pSH22_x_3 = fTmp1C_6 * _S259;
            float3  * _S272 = coeffs_5 + int(15);
            float3  * _S273 = coeffs_5 + int(23);
            float3  * _S274 = coeffs_5 + int(16);
            float3  * _S275 = coeffs_5 + int(22);
            float3  * _S276 = coeffs_5 + int(17);
            float3  * _S277 = coeffs_5 + int(21);
            float3  * _S278 = coeffs_5 + int(20);
            float3  * _S279 = coeffs_5 + int(18);
            float3  dir_n_15 = make_float3 (x_9, y_8, z_6);
            float3  v_dir_n_15 = make_float3 (v_x_7 + dot_0(_S250, make_float3 (0.62583571672439575f * (fS2_6 + y_8 * fC2_x_3 + x_9 * fS2_x_3)) * *_S272 + make_float3 (0.62583571672439575f * (fC2_6 + x_9 * fC2_x_3 - y_8 * fS2_x_3)) * *_S273 + make_float3 (fTmp2B_6 * fS2_x_3) * *_S274 + make_float3 (fTmp2B_6 * fC2_x_3) * *_S275 + make_float3 (fTmp1C_6 * fS1_x_3) * *_S276 + make_float3 (pSH22_x_3) * *_S277 + make_float3 (fTmp0D_6) * *_S278), v_y_7 + dot_0(_S250, make_float3 (0.62583571672439575f * (x_9 * fS2_y_3 + fC2_6 + y_8 * fC2_y_3)) * *_S272 + make_float3 (0.62583571672439575f * (x_9 * fC2_y_3 - fS2_6 - y_8 * fS2_y_3)) * *_S273 + make_float3 (fTmp2B_6 * fS2_y_3) * *_S274 + make_float3 (fTmp2B_6 * fC2_y_3) * *_S275 + make_float3 (pSH22_x_3) * *_S276 + make_float3 (fTmp1C_6 * fC1_y_3) * *_S277 + make_float3 (fTmp0D_6) * *_S279), v_z_7 + dot_0(_S250, make_float3 (1.9843134880065918f * (pSH12_6 + z_6 * pSH12_z_3) + -1.00623059272766113f * pSH6_z_3) * *(coeffs_5 + int(19)) + make_float3 (fTmp0D_z_3 * x_9) * *_S278 + make_float3 (fTmp0D_z_3 * y_8) * *_S279 + make_float3 (fTmp1C_z_3 * fC1_6) * *_S277 + make_float3 (fTmp1C_z_3 * fS1_6) * *_S276 + make_float3 (-1.77013075351715088f * fC2_6) * *_S275 + make_float3 (-1.77013075351715088f * fS2_6) * *_S274));
            v_viewdir_0 = (v_dir_n_15 - make_float3 (dot_0(v_dir_n_15, dir_n_15)) * dir_n_15) * make_float3 (inorm_3);
            break;
        }
        Matrix<float, 3, 3>  _S280 = makeMatrix<float, 3, 3> (0.0f);
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S281;
        (&_S281)->primal_0 = _S244;
        (&_S281)->differential_0 = _S280;
        float3  _S282 = make_float3 (0.0f);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S283;
        (&_S283)->primal_0 = t_1;
        (&_S283)->differential_0 = _S282;
        s_bwd_prop_mul_0(&_S281, &_S283, v_viewdir_0);
        Matrix<float, 3, 3>  _S284 = transpose_0(_S281.differential_0);
        *v_mean_0 = v_viewdir_0;
        *v_R_0 = _S284;
        *v_t_0 = _S283.differential_0;
        break;
    }
    return;
}

inline __device__ void sh_coeffs_to_color_vjp_inplace(int degree_6, float3  mean_2, Matrix<float, 3, 3>  R_2, float3  t_2, float3  coeff_dc_6, float3  * coeffs_6, float3  v_colors_4, float3  * v_coeff_dc_4, float3  * v_coeffs_4, float3  * v_mean_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    bool _S285;
    for(;;)
    {
        float3  v_viewdir_1;
        Matrix<float, 3, 3>  _S286 = transpose_0(R_2);
        float3  _S287 = mean_2 + mul_0(_S286, t_2);
        for(;;)
        {
            float3  colors_12 = make_float3 (0.282094806432724f) * coeff_dc_6;
            bool _S288 = degree_6 <= int(0);
            _S285 = _S288;
            if(_S288)
            {
                v_viewdir_1 = colors_12;
                break;
            }
            float _S289 = _S287.x;
            float _S290 = _S287.y;
            float _S291 = _S287.z;
            float inv_norm_3 = (F32_rsqrt((_S289 * _S289 + _S290 * _S290 + _S291 * _S291)));
            float x_10 = _S289 * inv_norm_3;
            float y_9 = _S290 * inv_norm_3;
            float z_7 = _S291 * inv_norm_3;
            float3  colors_13 = colors_12 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * *(coeffs_6 + int(0)) + make_float3 (z_7) * *(coeffs_6 + int(1)) - make_float3 (x_10) * *(coeffs_6 + int(2)));
            if(degree_6 <= int(1))
            {
                v_viewdir_1 = colors_13;
                break;
            }
            float z2_7 = z_7 * z_7;
            float fTmp0B_7 = -1.09254848957061768f * z_7;
            float fC1_7 = x_10 * x_10 - y_9 * y_9;
            float fS1_7 = 2.0f * x_10 * y_9;
            float pSH6_7 = 0.94617468118667603f * z2_7 - 0.31539157032966614f;
            float3  colors_14 = colors_13 + (make_float3 (0.54627424478530884f * fS1_7) * *(coeffs_6 + int(3)) + make_float3 (fTmp0B_7 * y_9) * *(coeffs_6 + int(4)) + make_float3 (pSH6_7) * *(coeffs_6 + int(5)) + make_float3 (fTmp0B_7 * x_10) * *(coeffs_6 + int(6)) + make_float3 (0.54627424478530884f * fC1_7) * *(coeffs_6 + int(7)));
            if(degree_6 <= int(2))
            {
                v_viewdir_1 = colors_14;
                break;
            }
            float fTmp0C_7 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
            float fTmp1B_7 = 1.44530570507049561f * z_7;
            float fC2_7 = x_10 * fC1_7 - y_9 * fS1_7;
            float fS2_7 = x_10 * fS1_7 + y_9 * fC1_7;
            float pSH12_7 = z_7 * (1.86588168144226074f * z2_7 - 1.11952900886535645f);
            float3  colors_15 = colors_14 + (make_float3 (-0.59004360437393188f * fS2_7) * *(coeffs_6 + int(8)) + make_float3 (fTmp1B_7 * fS1_7) * *(coeffs_6 + int(9)) + make_float3 (fTmp0C_7 * y_9) * *(coeffs_6 + int(10)) + make_float3 (pSH12_7) * *(coeffs_6 + int(11)) + make_float3 (fTmp0C_7 * x_10) * *(coeffs_6 + int(12)) + make_float3 (fTmp1B_7 * fC1_7) * *(coeffs_6 + int(13)) + make_float3 (-0.59004360437393188f * fC2_7) * *(coeffs_6 + int(14)));
            if(degree_6 <= int(3))
            {
                v_viewdir_1 = colors_15;
                break;
            }
            float fTmp0D_7 = z_7 * (-4.68332576751708984f * z2_7 + 2.00713968276977539f);
            float fTmp1C_7 = 3.31161141395568848f * z2_7 - 0.47308734059333801f;
            float fTmp2B_7 = -1.77013075351715088f * z_7;
            v_viewdir_1 = colors_15 + (make_float3 (0.62583571672439575f * (x_10 * fS2_7 + y_9 * fC2_7)) * *(coeffs_6 + int(15)) + make_float3 (fTmp2B_7 * fS2_7) * *(coeffs_6 + int(16)) + make_float3 (fTmp1C_7 * fS1_7) * *(coeffs_6 + int(17)) + make_float3 (fTmp0D_7 * y_9) * *(coeffs_6 + int(18)) + make_float3 (1.9843134880065918f * z_7 * pSH12_7 - 1.00623059272766113f * pSH6_7) * *(coeffs_6 + int(19)) + make_float3 (fTmp0D_7 * x_10) * *(coeffs_6 + int(20)) + make_float3 (fTmp1C_7 * fC1_7) * *(coeffs_6 + int(21)) + make_float3 (fTmp2B_7 * fC2_7) * *(coeffs_6 + int(22)) + make_float3 (0.62583571672439575f * (x_10 * fC2_7 - y_9 * fS2_7)) * *(coeffs_6 + int(23)));
            break;
        }
        float3  _S292 = v_colors_4 * make_float3 (float((v_viewdir_1.x) >= -0.5f), float((v_viewdir_1.y) >= -0.5f), float((v_viewdir_1.z) >= -0.5f));
        float3  v_viewdir_2 = {};
        for(;;)
        {
            *v_coeff_dc_4 = *v_coeff_dc_4 + make_float3 (0.282094806432724f) * _S292;
            if(_S285)
            {
                break;
            }
            float _S293 = _S287.x;
            float _S294 = _S287.y;
            float _S295 = _S287.z;
            float inorm_4 = (F32_rsqrt((_S293 * _S293 + _S294 * _S294 + _S295 * _S295)));
            float x_11 = _S293 * inorm_4;
            float y_10 = _S294 * inorm_4;
            float z_8 = _S295 * inorm_4;
            float3  * _S296 = v_coeffs_4 + int(0);
            *_S296 = *_S296 + make_float3 (-0.48860251903533936f * y_10) * _S292;
            float3  * _S297 = v_coeffs_4 + int(1);
            *_S297 = *_S297 + make_float3 (0.48860251903533936f * z_8) * _S292;
            float3  * _S298 = v_coeffs_4 + int(2);
            *_S298 = *_S298 + make_float3 (-0.48860251903533936f * x_11) * _S292;
            float _S299 = -0.48860251903533936f * dot_0(*(coeffs_6 + int(2)), _S292);
            float _S300 = -0.48860251903533936f * dot_0(*(coeffs_6 + int(0)), _S292);
            float _S301 = 0.48860251903533936f * dot_0(*(coeffs_6 + int(1)), _S292);
            if(degree_6 <= int(1))
            {
                float3  dir_n_16 = make_float3 (x_11, y_10, z_8);
                float3  v_dir_n_16 = make_float3 (_S299, _S300, _S301);
                v_viewdir_1 = v_viewdir_2 + (v_dir_n_16 - make_float3 (dot_0(v_dir_n_16, dir_n_16)) * dir_n_16) * make_float3 (inorm_4);
                break;
            }
            float z2_8 = z_8 * z_8;
            float fTmp0B_8 = -1.09254848957061768f * z_8;
            float fC1_8 = x_11 * x_11 - y_10 * y_10;
            float _S302 = 2.0f * x_11;
            float fS1_8 = _S302 * y_10;
            float pSH6_8 = 0.94617468118667603f * z2_8 - 0.31539157032966614f;
            float pSH7_4 = fTmp0B_8 * x_11;
            float pSH5_4 = fTmp0B_8 * y_10;
            float pSH8_4 = 0.54627424478530884f * fC1_8;
            float pSH4_4 = 0.54627424478530884f * fS1_8;
            float3  * _S303 = v_coeffs_4 + int(3);
            *_S303 = *_S303 + make_float3 (pSH4_4) * _S292;
            float3  * _S304 = v_coeffs_4 + int(4);
            *_S304 = *_S304 + make_float3 (pSH5_4) * _S292;
            float3  * _S305 = v_coeffs_4 + int(5);
            *_S305 = *_S305 + make_float3 (pSH6_8) * _S292;
            float3  * _S306 = v_coeffs_4 + int(6);
            *_S306 = *_S306 + make_float3 (pSH7_4) * _S292;
            float3  * _S307 = v_coeffs_4 + int(7);
            *_S307 = *_S307 + make_float3 (pSH8_4) * _S292;
            float fC1_y_4 = -2.0f * y_10;
            float fS1_x_4 = 2.0f * y_10;
            float pSH6_z_4 = 1.89234936237335205f * z_8;
            float pSH8_x_4 = 0.54627424478530884f * _S302;
            float3  * _S308 = coeffs_6 + int(3);
            float3  * _S309 = coeffs_6 + int(7);
            float3  * _S310 = coeffs_6 + int(6);
            float v_x_8 = _S299 + dot_0(_S292, make_float3 (0.54627424478530884f * fS1_x_4) * *_S308 + make_float3 (pSH8_x_4) * *_S309 + make_float3 (fTmp0B_8) * *_S310);
            float3  * _S311 = coeffs_6 + int(4);
            float v_y_8 = _S300 + dot_0(_S292, make_float3 (pSH8_x_4) * *_S308 + make_float3 (0.54627424478530884f * fC1_y_4) * *_S309 + make_float3 (fTmp0B_8) * *_S311);
            float v_z_8 = _S301 + dot_0(_S292, make_float3 (pSH6_z_4) * *(coeffs_6 + int(5)) + make_float3 (-1.09254848957061768f * x_11) * *_S310 + make_float3 (-1.09254848957061768f * y_10) * *_S311);
            if(degree_6 <= int(2))
            {
                float3  dir_n_17 = make_float3 (x_11, y_10, z_8);
                float3  v_dir_n_17 = make_float3 (v_x_8, v_y_8, v_z_8);
                v_viewdir_1 = v_viewdir_2 + (v_dir_n_17 - make_float3 (dot_0(v_dir_n_17, dir_n_17)) * dir_n_17) * make_float3 (inorm_4);
                break;
            }
            float fTmp0C_8 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
            float fTmp1B_8 = 1.44530570507049561f * z_8;
            float fC2_8 = x_11 * fC1_8 - y_10 * fS1_8;
            float fS2_8 = x_11 * fS1_8 + y_10 * fC1_8;
            float pSH12_8 = z_8 * (1.86588168144226074f * z2_8 - 1.11952900886535645f);
            float pSH13_4 = fTmp0C_8 * x_11;
            float pSH11_4 = fTmp0C_8 * y_10;
            float pSH14_4 = fTmp1B_8 * fC1_8;
            float pSH10_4 = fTmp1B_8 * fS1_8;
            float pSH15_4 = -0.59004360437393188f * fC2_8;
            float pSH9_4 = -0.59004360437393188f * fS2_8;
            float3  * _S312 = v_coeffs_4 + int(8);
            *_S312 = *_S312 + make_float3 (pSH9_4) * _S292;
            float3  * _S313 = v_coeffs_4 + int(9);
            *_S313 = *_S313 + make_float3 (pSH10_4) * _S292;
            float3  * _S314 = v_coeffs_4 + int(10);
            *_S314 = *_S314 + make_float3 (pSH11_4) * _S292;
            float3  * _S315 = v_coeffs_4 + int(11);
            *_S315 = *_S315 + make_float3 (pSH12_8) * _S292;
            float3  * _S316 = v_coeffs_4 + int(12);
            *_S316 = *_S316 + make_float3 (pSH13_4) * _S292;
            float3  * _S317 = v_coeffs_4 + int(13);
            *_S317 = *_S317 + make_float3 (pSH14_4) * _S292;
            float3  * _S318 = v_coeffs_4 + int(14);
            *_S318 = *_S318 + make_float3 (pSH15_4) * _S292;
            float fTmp0C_z_4 = -4.57045793533325195f * z_8;
            float _S319 = x_11 * _S302;
            float fC2_x_4 = fC1_8 + _S319 - y_10 * fS1_x_4;
            float _S320 = y_10 * _S302;
            float fC2_y_4 = x_11 * fC1_y_4 - fS1_8 - _S320;
            float fS2_x_4 = fS1_8 + x_11 * fS1_x_4 + _S320;
            float fS2_y_4 = _S319 + fC1_8 + y_10 * fC1_y_4;
            float pSH12_z_4 = 5.59764480590820312f * z2_8 - 1.11952900886535645f;
            float pSH14_x_4 = fTmp1B_8 * _S302;
            float3  * _S321 = coeffs_6 + int(8);
            float3  * _S322 = coeffs_6 + int(14);
            float3  * _S323 = coeffs_6 + int(9);
            float3  * _S324 = coeffs_6 + int(13);
            float3  * _S325 = coeffs_6 + int(12);
            float v_x_9 = v_x_8 + dot_0(_S292, make_float3 (-0.59004360437393188f * fS2_x_4) * *_S321 + make_float3 (-0.59004360437393188f * fC2_x_4) * *_S322 + make_float3 (fTmp1B_8 * fS1_x_4) * *_S323 + make_float3 (pSH14_x_4) * *_S324 + make_float3 (fTmp0C_8) * *_S325);
            float3  * _S326 = coeffs_6 + int(10);
            float v_y_9 = v_y_8 + dot_0(_S292, make_float3 (-0.59004360437393188f * fS2_y_4) * *_S321 + make_float3 (-0.59004360437393188f * fC2_y_4) * *_S322 + make_float3 (pSH14_x_4) * *_S323 + make_float3 (fTmp1B_8 * fC1_y_4) * *_S324 + make_float3 (fTmp0C_8) * *_S326);
            float v_z_9 = v_z_8 + dot_0(_S292, make_float3 (pSH12_z_4) * *(coeffs_6 + int(11)) + make_float3 (fTmp0C_z_4 * x_11) * *_S325 + make_float3 (fTmp0C_z_4 * y_10) * *_S326 + make_float3 (1.44530570507049561f * fC1_8) * *_S324 + make_float3 (1.44530570507049561f * fS1_8) * *_S323);
            if(degree_6 <= int(3))
            {
                float3  dir_n_18 = make_float3 (x_11, y_10, z_8);
                float3  v_dir_n_18 = make_float3 (v_x_9, v_y_9, v_z_9);
                v_viewdir_1 = v_viewdir_2 + (v_dir_n_18 - make_float3 (dot_0(v_dir_n_18, dir_n_18)) * dir_n_18) * make_float3 (inorm_4);
                break;
            }
            float fTmp0D_8 = z_8 * (-4.68332576751708984f * z2_8 + 2.00713968276977539f);
            float fTmp1C_8 = 3.31161141395568848f * z2_8 - 0.47308734059333801f;
            float fTmp2B_8 = -1.77013075351715088f * z_8;
            float pSH20_4 = 1.9843134880065918f * z_8 * pSH12_8 + -1.00623059272766113f * pSH6_8;
            float pSH21_4 = fTmp0D_8 * x_11;
            float pSH19_4 = fTmp0D_8 * y_10;
            float pSH22_4 = fTmp1C_8 * fC1_8;
            float pSH18_4 = fTmp1C_8 * fS1_8;
            float pSH23_4 = fTmp2B_8 * fC2_8;
            float pSH17_4 = fTmp2B_8 * fS2_8;
            float pSH24_4 = 0.62583571672439575f * (x_11 * fC2_8 - y_10 * fS2_8);
            float pSH16_4 = 0.62583571672439575f * (x_11 * fS2_8 + y_10 * fC2_8);
            float3  * _S327 = v_coeffs_4 + int(15);
            *_S327 = *_S327 + make_float3 (pSH16_4) * _S292;
            float3  * _S328 = v_coeffs_4 + int(16);
            *_S328 = *_S328 + make_float3 (pSH17_4) * _S292;
            float3  * _S329 = v_coeffs_4 + int(17);
            *_S329 = *_S329 + make_float3 (pSH18_4) * _S292;
            float3  * _S330 = v_coeffs_4 + int(18);
            *_S330 = *_S330 + make_float3 (pSH19_4) * _S292;
            float3  * _S331 = v_coeffs_4 + int(19);
            *_S331 = *_S331 + make_float3 (pSH20_4) * _S292;
            float3  * _S332 = v_coeffs_4 + int(20);
            *_S332 = *_S332 + make_float3 (pSH21_4) * _S292;
            float3  * _S333 = v_coeffs_4 + int(21);
            *_S333 = *_S333 + make_float3 (pSH22_4) * _S292;
            float3  * _S334 = v_coeffs_4 + int(22);
            *_S334 = *_S334 + make_float3 (pSH23_4) * _S292;
            float3  * _S335 = v_coeffs_4 + int(23);
            *_S335 = *_S335 + make_float3 (pSH24_4) * _S292;
            float fTmp0D_z_4 = -14.04997730255126953f * z2_8 + 2.00713968276977539f;
            float fTmp1C_z_4 = 6.62322282791137695f * z_8;
            float pSH22_x_4 = fTmp1C_8 * _S302;
            float3  * _S336 = coeffs_6 + int(15);
            float3  * _S337 = coeffs_6 + int(23);
            float3  * _S338 = coeffs_6 + int(16);
            float3  * _S339 = coeffs_6 + int(22);
            float3  * _S340 = coeffs_6 + int(17);
            float3  * _S341 = coeffs_6 + int(21);
            float3  * _S342 = coeffs_6 + int(20);
            float3  * _S343 = coeffs_6 + int(18);
            float3  dir_n_19 = make_float3 (x_11, y_10, z_8);
            float3  v_dir_n_19 = make_float3 (v_x_9 + dot_0(_S292, make_float3 (0.62583571672439575f * (fS2_8 + y_10 * fC2_x_4 + x_11 * fS2_x_4)) * *_S336 + make_float3 (0.62583571672439575f * (fC2_8 + x_11 * fC2_x_4 - y_10 * fS2_x_4)) * *_S337 + make_float3 (fTmp2B_8 * fS2_x_4) * *_S338 + make_float3 (fTmp2B_8 * fC2_x_4) * *_S339 + make_float3 (fTmp1C_8 * fS1_x_4) * *_S340 + make_float3 (pSH22_x_4) * *_S341 + make_float3 (fTmp0D_8) * *_S342), v_y_9 + dot_0(_S292, make_float3 (0.62583571672439575f * (x_11 * fS2_y_4 + fC2_8 + y_10 * fC2_y_4)) * *_S336 + make_float3 (0.62583571672439575f * (x_11 * fC2_y_4 - fS2_8 - y_10 * fS2_y_4)) * *_S337 + make_float3 (fTmp2B_8 * fS2_y_4) * *_S338 + make_float3 (fTmp2B_8 * fC2_y_4) * *_S339 + make_float3 (pSH22_x_4) * *_S340 + make_float3 (fTmp1C_8 * fC1_y_4) * *_S341 + make_float3 (fTmp0D_8) * *_S343), v_z_9 + dot_0(_S292, make_float3 (1.9843134880065918f * (pSH12_8 + z_8 * pSH12_z_4) + -1.00623059272766113f * pSH6_z_4) * *(coeffs_6 + int(19)) + make_float3 (fTmp0D_z_4 * x_11) * *_S342 + make_float3 (fTmp0D_z_4 * y_10) * *_S343 + make_float3 (fTmp1C_z_4 * fC1_8) * *_S341 + make_float3 (fTmp1C_z_4 * fS1_8) * *_S340 + make_float3 (-1.77013075351715088f * fC2_8) * *_S339 + make_float3 (-1.77013075351715088f * fS2_8) * *_S338));
            v_viewdir_1 = v_viewdir_2 + (v_dir_n_19 - make_float3 (dot_0(v_dir_n_19, dir_n_19)) * dir_n_19) * make_float3 (inorm_4);
            break;
        }
        Matrix<float, 3, 3>  _S344 = makeMatrix<float, 3, 3> (0.0f);
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S345;
        (&_S345)->primal_0 = _S286;
        (&_S345)->differential_0 = _S344;
        float3  _S346 = make_float3 (0.0f);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S347;
        (&_S347)->primal_0 = t_2;
        (&_S347)->differential_0 = _S346;
        s_bwd_prop_mul_0(&_S345, &_S347, v_viewdir_1);
        Matrix<float, 3, 3>  _S348 = transpose_0(_S345.differential_0);
        *v_mean_1 = *v_mean_1 + v_viewdir_1;
        *v_R_1 = *v_R_1 + _S348;
        *v_t_1 = *v_t_1 + _S347.differential_0;
        break;
    }
    return;
}

inline __device__ void sh_coeffs_to_color_vjp_atomic(int degree_7, float3  mean_3, Matrix<float, 3, 3>  R_3, float3  t_3, float3  coeff_dc_7, float3  * coeffs_7, float3  v_colors_5, float3  * v_coeff_dc_5, float3  * v_coeffs_5, float3  * v_mean_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    bool _S349;
    for(;;)
    {
        float3  v_viewdir_3;
        Matrix<float, 3, 3>  _S350 = transpose_0(R_3);
        float3  _S351 = mean_3 + mul_0(_S350, t_3);
        for(;;)
        {
            float3  colors_16 = make_float3 (0.282094806432724f) * coeff_dc_7;
            bool _S352 = degree_7 <= int(0);
            _S349 = _S352;
            if(_S352)
            {
                v_viewdir_3 = colors_16;
                break;
            }
            float _S353 = _S351.x;
            float _S354 = _S351.y;
            float _S355 = _S351.z;
            float inv_norm_4 = (F32_rsqrt((_S353 * _S353 + _S354 * _S354 + _S355 * _S355)));
            float x_12 = _S353 * inv_norm_4;
            float y_11 = _S354 * inv_norm_4;
            float z_9 = _S355 * inv_norm_4;
            float3  colors_17 = colors_16 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * *(coeffs_7 + int(0)) + make_float3 (z_9) * *(coeffs_7 + int(1)) - make_float3 (x_12) * *(coeffs_7 + int(2)));
            if(degree_7 <= int(1))
            {
                v_viewdir_3 = colors_17;
                break;
            }
            float z2_9 = z_9 * z_9;
            float fTmp0B_9 = -1.09254848957061768f * z_9;
            float fC1_9 = x_12 * x_12 - y_11 * y_11;
            float fS1_9 = 2.0f * x_12 * y_11;
            float pSH6_9 = 0.94617468118667603f * z2_9 - 0.31539157032966614f;
            float3  colors_18 = colors_17 + (make_float3 (0.54627424478530884f * fS1_9) * *(coeffs_7 + int(3)) + make_float3 (fTmp0B_9 * y_11) * *(coeffs_7 + int(4)) + make_float3 (pSH6_9) * *(coeffs_7 + int(5)) + make_float3 (fTmp0B_9 * x_12) * *(coeffs_7 + int(6)) + make_float3 (0.54627424478530884f * fC1_9) * *(coeffs_7 + int(7)));
            if(degree_7 <= int(2))
            {
                v_viewdir_3 = colors_18;
                break;
            }
            float fTmp0C_9 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
            float fTmp1B_9 = 1.44530570507049561f * z_9;
            float fC2_9 = x_12 * fC1_9 - y_11 * fS1_9;
            float fS2_9 = x_12 * fS1_9 + y_11 * fC1_9;
            float pSH12_9 = z_9 * (1.86588168144226074f * z2_9 - 1.11952900886535645f);
            float3  colors_19 = colors_18 + (make_float3 (-0.59004360437393188f * fS2_9) * *(coeffs_7 + int(8)) + make_float3 (fTmp1B_9 * fS1_9) * *(coeffs_7 + int(9)) + make_float3 (fTmp0C_9 * y_11) * *(coeffs_7 + int(10)) + make_float3 (pSH12_9) * *(coeffs_7 + int(11)) + make_float3 (fTmp0C_9 * x_12) * *(coeffs_7 + int(12)) + make_float3 (fTmp1B_9 * fC1_9) * *(coeffs_7 + int(13)) + make_float3 (-0.59004360437393188f * fC2_9) * *(coeffs_7 + int(14)));
            if(degree_7 <= int(3))
            {
                v_viewdir_3 = colors_19;
                break;
            }
            float fTmp0D_9 = z_9 * (-4.68332576751708984f * z2_9 + 2.00713968276977539f);
            float fTmp1C_9 = 3.31161141395568848f * z2_9 - 0.47308734059333801f;
            float fTmp2B_9 = -1.77013075351715088f * z_9;
            v_viewdir_3 = colors_19 + (make_float3 (0.62583571672439575f * (x_12 * fS2_9 + y_11 * fC2_9)) * *(coeffs_7 + int(15)) + make_float3 (fTmp2B_9 * fS2_9) * *(coeffs_7 + int(16)) + make_float3 (fTmp1C_9 * fS1_9) * *(coeffs_7 + int(17)) + make_float3 (fTmp0D_9 * y_11) * *(coeffs_7 + int(18)) + make_float3 (1.9843134880065918f * z_9 * pSH12_9 - 1.00623059272766113f * pSH6_9) * *(coeffs_7 + int(19)) + make_float3 (fTmp0D_9 * x_12) * *(coeffs_7 + int(20)) + make_float3 (fTmp1C_9 * fC1_9) * *(coeffs_7 + int(21)) + make_float3 (fTmp2B_9 * fC2_9) * *(coeffs_7 + int(22)) + make_float3 (0.62583571672439575f * (x_12 * fC2_9 - y_11 * fS2_9)) * *(coeffs_7 + int(23)));
            break;
        }
        float3  _S356 = v_colors_5 * make_float3 (float((v_viewdir_3.x) >= -0.5f), float((v_viewdir_3.y) >= -0.5f), float((v_viewdir_3.z) >= -0.5f));
        float3  v_viewdir_4 = {};
        for(;;)
        {
            *v_coeff_dc_5 = *v_coeff_dc_5 + make_float3 (0.282094806432724f) * _S356;
            if(_S349)
            {
                break;
            }
            float _S357 = _S351.x;
            float _S358 = _S351.y;
            float _S359 = _S351.z;
            float inorm_5 = (F32_rsqrt((_S357 * _S357 + _S358 * _S358 + _S359 * _S359)));
            float x_13 = _S357 * inorm_5;
            float y_12 = _S358 * inorm_5;
            float z_10 = _S359 * inorm_5;
            float3  temp_24 = make_float3 (-0.48860251903533936f * y_12) * _S356;
            float _S360 = dot_0(temp_24, temp_24);
            bool _S361;
            if((F32_isfinite((_S360))))
            {
                _S361 = _S360 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S362 = v_coeffs_5 + int(0);
                float _S363 = atomicAdd(&(_S362->x), temp_24.x);
                float _S364 = atomicAdd(&(_S362->y), temp_24.y);
                float _S365 = atomicAdd(&(_S362->z), temp_24.z);
            }
            float3  temp_25 = make_float3 (0.48860251903533936f * z_10) * _S356;
            float _S366 = dot_0(temp_25, temp_25);
            if((F32_isfinite((_S366))))
            {
                _S361 = _S366 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S367 = v_coeffs_5 + int(1);
                float _S368 = atomicAdd(&(_S367->x), temp_25.x);
                float _S369 = atomicAdd(&(_S367->y), temp_25.y);
                float _S370 = atomicAdd(&(_S367->z), temp_25.z);
            }
            float3  temp_26 = make_float3 (-0.48860251903533936f * x_13) * _S356;
            float _S371 = dot_0(temp_26, temp_26);
            if((F32_isfinite((_S371))))
            {
                _S361 = _S371 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S372 = v_coeffs_5 + int(2);
                float _S373 = atomicAdd(&(_S372->x), temp_26.x);
                float _S374 = atomicAdd(&(_S372->y), temp_26.y);
                float _S375 = atomicAdd(&(_S372->z), temp_26.z);
            }
            float _S376 = -0.48860251903533936f * dot_0(*(coeffs_7 + int(2)), _S356);
            float _S377 = -0.48860251903533936f * dot_0(*(coeffs_7 + int(0)), _S356);
            float _S378 = 0.48860251903533936f * dot_0(*(coeffs_7 + int(1)), _S356);
            if(degree_7 <= int(1))
            {
                float3  dir_n_20 = make_float3 (x_13, y_12, z_10);
                float3  v_dir_n_20 = make_float3 (_S376, _S377, _S378);
                v_viewdir_3 = v_viewdir_4 + (v_dir_n_20 - make_float3 (dot_0(v_dir_n_20, dir_n_20)) * dir_n_20) * make_float3 (inorm_5);
                break;
            }
            float z2_10 = z_10 * z_10;
            float fTmp0B_10 = -1.09254848957061768f * z_10;
            float fC1_10 = x_13 * x_13 - y_12 * y_12;
            float _S379 = 2.0f * x_13;
            float fS1_10 = _S379 * y_12;
            float pSH6_10 = 0.94617468118667603f * z2_10 - 0.31539157032966614f;
            float pSH7_5 = fTmp0B_10 * x_13;
            float pSH5_5 = fTmp0B_10 * y_12;
            float pSH8_5 = 0.54627424478530884f * fC1_10;
            float pSH4_5 = 0.54627424478530884f * fS1_10;
            float3  temp_27 = make_float3 (pSH4_5) * _S356;
            float _S380 = dot_0(temp_27, temp_27);
            if((F32_isfinite((_S380))))
            {
                _S361 = _S380 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S381 = v_coeffs_5 + int(3);
                float _S382 = atomicAdd(&(_S381->x), temp_27.x);
                float _S383 = atomicAdd(&(_S381->y), temp_27.y);
                float _S384 = atomicAdd(&(_S381->z), temp_27.z);
            }
            float3  temp_28 = make_float3 (pSH5_5) * _S356;
            float _S385 = dot_0(temp_28, temp_28);
            if((F32_isfinite((_S385))))
            {
                _S361 = _S385 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S386 = v_coeffs_5 + int(4);
                float _S387 = atomicAdd(&(_S386->x), temp_28.x);
                float _S388 = atomicAdd(&(_S386->y), temp_28.y);
                float _S389 = atomicAdd(&(_S386->z), temp_28.z);
            }
            float3  temp_29 = make_float3 (pSH6_10) * _S356;
            float _S390 = dot_0(temp_29, temp_29);
            if((F32_isfinite((_S390))))
            {
                _S361 = _S390 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S391 = v_coeffs_5 + int(5);
                float _S392 = atomicAdd(&(_S391->x), temp_29.x);
                float _S393 = atomicAdd(&(_S391->y), temp_29.y);
                float _S394 = atomicAdd(&(_S391->z), temp_29.z);
            }
            float3  temp_30 = make_float3 (pSH7_5) * _S356;
            float _S395 = dot_0(temp_30, temp_30);
            if((F32_isfinite((_S395))))
            {
                _S361 = _S395 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S396 = v_coeffs_5 + int(6);
                float _S397 = atomicAdd(&(_S396->x), temp_30.x);
                float _S398 = atomicAdd(&(_S396->y), temp_30.y);
                float _S399 = atomicAdd(&(_S396->z), temp_30.z);
            }
            float3  temp_31 = make_float3 (pSH8_5) * _S356;
            float _S400 = dot_0(temp_31, temp_31);
            if((F32_isfinite((_S400))))
            {
                _S361 = _S400 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S401 = v_coeffs_5 + int(7);
                float _S402 = atomicAdd(&(_S401->x), temp_31.x);
                float _S403 = atomicAdd(&(_S401->y), temp_31.y);
                float _S404 = atomicAdd(&(_S401->z), temp_31.z);
            }
            float fC1_y_5 = -2.0f * y_12;
            float fS1_x_5 = 2.0f * y_12;
            float pSH6_z_5 = 1.89234936237335205f * z_10;
            float pSH8_x_5 = 0.54627424478530884f * _S379;
            float3  * _S405 = coeffs_7 + int(3);
            float3  * _S406 = coeffs_7 + int(7);
            float3  * _S407 = coeffs_7 + int(6);
            float v_x_10 = _S376 + dot_0(_S356, make_float3 (0.54627424478530884f * fS1_x_5) * *_S405 + make_float3 (pSH8_x_5) * *_S406 + make_float3 (fTmp0B_10) * *_S407);
            float3  * _S408 = coeffs_7 + int(4);
            float v_y_10 = _S377 + dot_0(_S356, make_float3 (pSH8_x_5) * *_S405 + make_float3 (0.54627424478530884f * fC1_y_5) * *_S406 + make_float3 (fTmp0B_10) * *_S408);
            float v_z_10 = _S378 + dot_0(_S356, make_float3 (pSH6_z_5) * *(coeffs_7 + int(5)) + make_float3 (-1.09254848957061768f * x_13) * *_S407 + make_float3 (-1.09254848957061768f * y_12) * *_S408);
            if(degree_7 <= int(2))
            {
                float3  dir_n_21 = make_float3 (x_13, y_12, z_10);
                float3  v_dir_n_21 = make_float3 (v_x_10, v_y_10, v_z_10);
                v_viewdir_3 = v_viewdir_4 + (v_dir_n_21 - make_float3 (dot_0(v_dir_n_21, dir_n_21)) * dir_n_21) * make_float3 (inorm_5);
                break;
            }
            float fTmp0C_10 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
            float fTmp1B_10 = 1.44530570507049561f * z_10;
            float fC2_10 = x_13 * fC1_10 - y_12 * fS1_10;
            float fS2_10 = x_13 * fS1_10 + y_12 * fC1_10;
            float pSH12_10 = z_10 * (1.86588168144226074f * z2_10 - 1.11952900886535645f);
            float pSH13_5 = fTmp0C_10 * x_13;
            float pSH11_5 = fTmp0C_10 * y_12;
            float pSH14_5 = fTmp1B_10 * fC1_10;
            float pSH10_5 = fTmp1B_10 * fS1_10;
            float pSH15_5 = -0.59004360437393188f * fC2_10;
            float pSH9_5 = -0.59004360437393188f * fS2_10;
            float3  temp_32 = make_float3 (pSH9_5) * _S356;
            float _S409 = dot_0(temp_32, temp_32);
            if((F32_isfinite((_S409))))
            {
                _S361 = _S409 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S410 = v_coeffs_5 + int(8);
                float _S411 = atomicAdd(&(_S410->x), temp_32.x);
                float _S412 = atomicAdd(&(_S410->y), temp_32.y);
                float _S413 = atomicAdd(&(_S410->z), temp_32.z);
            }
            float3  temp_33 = make_float3 (pSH10_5) * _S356;
            float _S414 = dot_0(temp_33, temp_33);
            if((F32_isfinite((_S414))))
            {
                _S361 = _S414 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S415 = v_coeffs_5 + int(9);
                float _S416 = atomicAdd(&(_S415->x), temp_33.x);
                float _S417 = atomicAdd(&(_S415->y), temp_33.y);
                float _S418 = atomicAdd(&(_S415->z), temp_33.z);
            }
            float3  temp_34 = make_float3 (pSH11_5) * _S356;
            float _S419 = dot_0(temp_34, temp_34);
            if((F32_isfinite((_S419))))
            {
                _S361 = _S419 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S420 = v_coeffs_5 + int(10);
                float _S421 = atomicAdd(&(_S420->x), temp_34.x);
                float _S422 = atomicAdd(&(_S420->y), temp_34.y);
                float _S423 = atomicAdd(&(_S420->z), temp_34.z);
            }
            float3  temp_35 = make_float3 (pSH12_10) * _S356;
            float _S424 = dot_0(temp_35, temp_35);
            if((F32_isfinite((_S424))))
            {
                _S361 = _S424 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S425 = v_coeffs_5 + int(11);
                float _S426 = atomicAdd(&(_S425->x), temp_35.x);
                float _S427 = atomicAdd(&(_S425->y), temp_35.y);
                float _S428 = atomicAdd(&(_S425->z), temp_35.z);
            }
            float3  temp_36 = make_float3 (pSH13_5) * _S356;
            float _S429 = dot_0(temp_36, temp_36);
            if((F32_isfinite((_S429))))
            {
                _S361 = _S429 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S430 = v_coeffs_5 + int(12);
                float _S431 = atomicAdd(&(_S430->x), temp_36.x);
                float _S432 = atomicAdd(&(_S430->y), temp_36.y);
                float _S433 = atomicAdd(&(_S430->z), temp_36.z);
            }
            float3  temp_37 = make_float3 (pSH14_5) * _S356;
            float _S434 = dot_0(temp_37, temp_37);
            if((F32_isfinite((_S434))))
            {
                _S361 = _S434 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S435 = v_coeffs_5 + int(13);
                float _S436 = atomicAdd(&(_S435->x), temp_37.x);
                float _S437 = atomicAdd(&(_S435->y), temp_37.y);
                float _S438 = atomicAdd(&(_S435->z), temp_37.z);
            }
            float3  temp_38 = make_float3 (pSH15_5) * _S356;
            float _S439 = dot_0(temp_38, temp_38);
            if((F32_isfinite((_S439))))
            {
                _S361 = _S439 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S440 = v_coeffs_5 + int(14);
                float _S441 = atomicAdd(&(_S440->x), temp_38.x);
                float _S442 = atomicAdd(&(_S440->y), temp_38.y);
                float _S443 = atomicAdd(&(_S440->z), temp_38.z);
            }
            float fTmp0C_z_5 = -4.57045793533325195f * z_10;
            float _S444 = x_13 * _S379;
            float fC2_x_5 = fC1_10 + _S444 - y_12 * fS1_x_5;
            float _S445 = y_12 * _S379;
            float fC2_y_5 = x_13 * fC1_y_5 - fS1_10 - _S445;
            float fS2_x_5 = fS1_10 + x_13 * fS1_x_5 + _S445;
            float fS2_y_5 = _S444 + fC1_10 + y_12 * fC1_y_5;
            float pSH12_z_5 = 5.59764480590820312f * z2_10 - 1.11952900886535645f;
            float pSH14_x_5 = fTmp1B_10 * _S379;
            float3  * _S446 = coeffs_7 + int(8);
            float3  * _S447 = coeffs_7 + int(14);
            float3  * _S448 = coeffs_7 + int(9);
            float3  * _S449 = coeffs_7 + int(13);
            float3  * _S450 = coeffs_7 + int(12);
            float v_x_11 = v_x_10 + dot_0(_S356, make_float3 (-0.59004360437393188f * fS2_x_5) * *_S446 + make_float3 (-0.59004360437393188f * fC2_x_5) * *_S447 + make_float3 (fTmp1B_10 * fS1_x_5) * *_S448 + make_float3 (pSH14_x_5) * *_S449 + make_float3 (fTmp0C_10) * *_S450);
            float3  * _S451 = coeffs_7 + int(10);
            float v_y_11 = v_y_10 + dot_0(_S356, make_float3 (-0.59004360437393188f * fS2_y_5) * *_S446 + make_float3 (-0.59004360437393188f * fC2_y_5) * *_S447 + make_float3 (pSH14_x_5) * *_S448 + make_float3 (fTmp1B_10 * fC1_y_5) * *_S449 + make_float3 (fTmp0C_10) * *_S451);
            float v_z_11 = v_z_10 + dot_0(_S356, make_float3 (pSH12_z_5) * *(coeffs_7 + int(11)) + make_float3 (fTmp0C_z_5 * x_13) * *_S450 + make_float3 (fTmp0C_z_5 * y_12) * *_S451 + make_float3 (1.44530570507049561f * fC1_10) * *_S449 + make_float3 (1.44530570507049561f * fS1_10) * *_S448);
            if(degree_7 <= int(3))
            {
                float3  dir_n_22 = make_float3 (x_13, y_12, z_10);
                float3  v_dir_n_22 = make_float3 (v_x_11, v_y_11, v_z_11);
                v_viewdir_3 = v_viewdir_4 + (v_dir_n_22 - make_float3 (dot_0(v_dir_n_22, dir_n_22)) * dir_n_22) * make_float3 (inorm_5);
                break;
            }
            float fTmp0D_10 = z_10 * (-4.68332576751708984f * z2_10 + 2.00713968276977539f);
            float fTmp1C_10 = 3.31161141395568848f * z2_10 - 0.47308734059333801f;
            float fTmp2B_10 = -1.77013075351715088f * z_10;
            float pSH20_5 = 1.9843134880065918f * z_10 * pSH12_10 + -1.00623059272766113f * pSH6_10;
            float pSH21_5 = fTmp0D_10 * x_13;
            float pSH19_5 = fTmp0D_10 * y_12;
            float pSH22_5 = fTmp1C_10 * fC1_10;
            float pSH18_5 = fTmp1C_10 * fS1_10;
            float pSH23_5 = fTmp2B_10 * fC2_10;
            float pSH17_5 = fTmp2B_10 * fS2_10;
            float pSH24_5 = 0.62583571672439575f * (x_13 * fC2_10 - y_12 * fS2_10);
            float pSH16_5 = 0.62583571672439575f * (x_13 * fS2_10 + y_12 * fC2_10);
            float3  temp_39 = make_float3 (pSH16_5) * _S356;
            float _S452 = dot_0(temp_39, temp_39);
            if((F32_isfinite((_S452))))
            {
                _S361 = _S452 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S453 = v_coeffs_5 + int(15);
                float _S454 = atomicAdd(&(_S453->x), temp_39.x);
                float _S455 = atomicAdd(&(_S453->y), temp_39.y);
                float _S456 = atomicAdd(&(_S453->z), temp_39.z);
            }
            float3  temp_40 = make_float3 (pSH17_5) * _S356;
            float _S457 = dot_0(temp_40, temp_40);
            if((F32_isfinite((_S457))))
            {
                _S361 = _S457 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S458 = v_coeffs_5 + int(16);
                float _S459 = atomicAdd(&(_S458->x), temp_40.x);
                float _S460 = atomicAdd(&(_S458->y), temp_40.y);
                float _S461 = atomicAdd(&(_S458->z), temp_40.z);
            }
            float3  temp_41 = make_float3 (pSH18_5) * _S356;
            float _S462 = dot_0(temp_41, temp_41);
            if((F32_isfinite((_S462))))
            {
                _S361 = _S462 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S463 = v_coeffs_5 + int(17);
                float _S464 = atomicAdd(&(_S463->x), temp_41.x);
                float _S465 = atomicAdd(&(_S463->y), temp_41.y);
                float _S466 = atomicAdd(&(_S463->z), temp_41.z);
            }
            float3  temp_42 = make_float3 (pSH19_5) * _S356;
            float _S467 = dot_0(temp_42, temp_42);
            if((F32_isfinite((_S467))))
            {
                _S361 = _S467 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S468 = v_coeffs_5 + int(18);
                float _S469 = atomicAdd(&(_S468->x), temp_42.x);
                float _S470 = atomicAdd(&(_S468->y), temp_42.y);
                float _S471 = atomicAdd(&(_S468->z), temp_42.z);
            }
            float3  temp_43 = make_float3 (pSH20_5) * _S356;
            float _S472 = dot_0(temp_43, temp_43);
            if((F32_isfinite((_S472))))
            {
                _S361 = _S472 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S473 = v_coeffs_5 + int(19);
                float _S474 = atomicAdd(&(_S473->x), temp_43.x);
                float _S475 = atomicAdd(&(_S473->y), temp_43.y);
                float _S476 = atomicAdd(&(_S473->z), temp_43.z);
            }
            float3  temp_44 = make_float3 (pSH21_5) * _S356;
            float _S477 = dot_0(temp_44, temp_44);
            if((F32_isfinite((_S477))))
            {
                _S361 = _S477 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S478 = v_coeffs_5 + int(20);
                float _S479 = atomicAdd(&(_S478->x), temp_44.x);
                float _S480 = atomicAdd(&(_S478->y), temp_44.y);
                float _S481 = atomicAdd(&(_S478->z), temp_44.z);
            }
            float3  temp_45 = make_float3 (pSH22_5) * _S356;
            float _S482 = dot_0(temp_45, temp_45);
            if((F32_isfinite((_S482))))
            {
                _S361 = _S482 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S483 = v_coeffs_5 + int(21);
                float _S484 = atomicAdd(&(_S483->x), temp_45.x);
                float _S485 = atomicAdd(&(_S483->y), temp_45.y);
                float _S486 = atomicAdd(&(_S483->z), temp_45.z);
            }
            float3  temp_46 = make_float3 (pSH23_5) * _S356;
            float _S487 = dot_0(temp_46, temp_46);
            if((F32_isfinite((_S487))))
            {
                _S361 = _S487 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S488 = v_coeffs_5 + int(22);
                float _S489 = atomicAdd(&(_S488->x), temp_46.x);
                float _S490 = atomicAdd(&(_S488->y), temp_46.y);
                float _S491 = atomicAdd(&(_S488->z), temp_46.z);
            }
            float3  temp_47 = make_float3 (pSH24_5) * _S356;
            float _S492 = dot_0(temp_47, temp_47);
            if((F32_isfinite((_S492))))
            {
                _S361 = _S492 != 0.0f;
            }
            else
            {
                _S361 = false;
            }
            if(_S361)
            {
                float3  * _S493 = v_coeffs_5 + int(23);
                float _S494 = atomicAdd(&(_S493->x), temp_47.x);
                float _S495 = atomicAdd(&(_S493->y), temp_47.y);
                float _S496 = atomicAdd(&(_S493->z), temp_47.z);
            }
            float fTmp0D_z_5 = -14.04997730255126953f * z2_10 + 2.00713968276977539f;
            float fTmp1C_z_5 = 6.62322282791137695f * z_10;
            float pSH22_x_5 = fTmp1C_10 * _S379;
            float3  * _S497 = coeffs_7 + int(15);
            float3  * _S498 = coeffs_7 + int(23);
            float3  * _S499 = coeffs_7 + int(16);
            float3  * _S500 = coeffs_7 + int(22);
            float3  * _S501 = coeffs_7 + int(17);
            float3  * _S502 = coeffs_7 + int(21);
            float3  * _S503 = coeffs_7 + int(20);
            float3  * _S504 = coeffs_7 + int(18);
            float3  dir_n_23 = make_float3 (x_13, y_12, z_10);
            float3  v_dir_n_23 = make_float3 (v_x_11 + dot_0(_S356, make_float3 (0.62583571672439575f * (fS2_10 + y_12 * fC2_x_5 + x_13 * fS2_x_5)) * *_S497 + make_float3 (0.62583571672439575f * (fC2_10 + x_13 * fC2_x_5 - y_12 * fS2_x_5)) * *_S498 + make_float3 (fTmp2B_10 * fS2_x_5) * *_S499 + make_float3 (fTmp2B_10 * fC2_x_5) * *_S500 + make_float3 (fTmp1C_10 * fS1_x_5) * *_S501 + make_float3 (pSH22_x_5) * *_S502 + make_float3 (fTmp0D_10) * *_S503), v_y_11 + dot_0(_S356, make_float3 (0.62583571672439575f * (x_13 * fS2_y_5 + fC2_10 + y_12 * fC2_y_5)) * *_S497 + make_float3 (0.62583571672439575f * (x_13 * fC2_y_5 - fS2_10 - y_12 * fS2_y_5)) * *_S498 + make_float3 (fTmp2B_10 * fS2_y_5) * *_S499 + make_float3 (fTmp2B_10 * fC2_y_5) * *_S500 + make_float3 (pSH22_x_5) * *_S501 + make_float3 (fTmp1C_10 * fC1_y_5) * *_S502 + make_float3 (fTmp0D_10) * *_S504), v_z_11 + dot_0(_S356, make_float3 (1.9843134880065918f * (pSH12_10 + z_10 * pSH12_z_5) + -1.00623059272766113f * pSH6_z_5) * *(coeffs_7 + int(19)) + make_float3 (fTmp0D_z_5 * x_13) * *_S503 + make_float3 (fTmp0D_z_5 * y_12) * *_S504 + make_float3 (fTmp1C_z_5 * fC1_10) * *_S502 + make_float3 (fTmp1C_z_5 * fS1_10) * *_S501 + make_float3 (-1.77013075351715088f * fC2_10) * *_S500 + make_float3 (-1.77013075351715088f * fS2_10) * *_S499));
            v_viewdir_3 = v_viewdir_4 + (v_dir_n_23 - make_float3 (dot_0(v_dir_n_23, dir_n_23)) * dir_n_23) * make_float3 (inorm_5);
            break;
        }
        Matrix<float, 3, 3>  _S505 = makeMatrix<float, 3, 3> (0.0f);
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S506;
        (&_S506)->primal_0 = _S350;
        (&_S506)->differential_0 = _S505;
        float3  _S507 = make_float3 (0.0f);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S508;
        (&_S508)->primal_0 = t_3;
        (&_S508)->differential_0 = _S507;
        s_bwd_prop_mul_0(&_S506, &_S508, v_viewdir_3);
        Matrix<float, 3, 3>  _S509 = transpose_0(_S506.differential_0);
        *v_mean_2 = *v_mean_2 + v_viewdir_3;
        *v_R_2 = *v_R_2 + _S509;
        *v_t_2 = *v_t_2 + _S508.differential_0;
        break;
    }
    return;
}

inline __device__ float3  sh0_to_color(float3  viewdir_1, float3  coeff_dc_8, float3  * coeffs_8)
{
    return make_float3 (0.282094806432724f) * coeff_dc_8;
}

inline __device__ float3  sh1_to_color(float3  viewdir_2, float3  coeff_dc_9, float3  * coeffs_9)
{
    float _S510 = viewdir_2.x;
    float _S511 = viewdir_2.y;
    float _S512 = viewdir_2.z;
    float inv_norm_5 = (F32_rsqrt((_S510 * _S510 + _S511 * _S511 + _S512 * _S512)));
    return make_float3 (0.282094806432724f) * coeff_dc_9 + make_float3 (0.48860251903533936f) * (make_float3 (- (_S511 * inv_norm_5)) * *(coeffs_9 + int(0)) + make_float3 (_S512 * inv_norm_5) * *(coeffs_9 + int(1)) - make_float3 (_S510 * inv_norm_5) * *(coeffs_9 + int(2)));
}

inline __device__ float3  sh2_to_color(float3  viewdir_3, float3  coeff_dc_10, float3  * coeffs_10)
{
    float _S513 = viewdir_3.x;
    float _S514 = viewdir_3.y;
    float _S515 = viewdir_3.z;
    float inv_norm_6 = (F32_rsqrt((_S513 * _S513 + _S514 * _S514 + _S515 * _S515)));
    float x_14 = _S513 * inv_norm_6;
    float y_13 = _S514 * inv_norm_6;
    float z_11 = _S515 * inv_norm_6;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    return make_float3 (0.282094806432724f) * coeff_dc_10 + make_float3 (0.48860251903533936f) * (make_float3 (- y_13) * *(coeffs_10 + int(0)) + make_float3 (z_11) * *(coeffs_10 + int(1)) - make_float3 (x_14) * *(coeffs_10 + int(2))) + (make_float3 (0.54627424478530884f * (2.0f * x_14 * y_13)) * *(coeffs_10 + int(3)) + make_float3 (fTmp0B_11 * y_13) * *(coeffs_10 + int(4)) + make_float3 (0.94617468118667603f * (z_11 * z_11) - 0.31539157032966614f) * *(coeffs_10 + int(5)) + make_float3 (fTmp0B_11 * x_14) * *(coeffs_10 + int(6)) + make_float3 (0.54627424478530884f * (x_14 * x_14 - y_13 * y_13)) * *(coeffs_10 + int(7)));
}

inline __device__ float3  sh3_to_color(float3  viewdir_4, float3  coeff_dc_11, float3  * coeffs_11)
{
    float _S516 = viewdir_4.x;
    float _S517 = viewdir_4.y;
    float _S518 = viewdir_4.z;
    float inv_norm_7 = (F32_rsqrt((_S516 * _S516 + _S517 * _S517 + _S518 * _S518)));
    float x_15 = _S516 * inv_norm_7;
    float y_14 = _S517 * inv_norm_7;
    float z_12 = _S518 * inv_norm_7;
    float z2_11 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_11 = x_15 * x_15 - y_14 * y_14;
    float fS1_11 = 2.0f * x_15 * y_14;
    float fTmp0C_11 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_12;
    return make_float3 (0.282094806432724f) * coeff_dc_11 + make_float3 (0.48860251903533936f) * (make_float3 (- y_14) * *(coeffs_11 + int(0)) + make_float3 (z_12) * *(coeffs_11 + int(1)) - make_float3 (x_15) * *(coeffs_11 + int(2))) + (make_float3 (0.54627424478530884f * fS1_11) * *(coeffs_11 + int(3)) + make_float3 (fTmp0B_12 * y_14) * *(coeffs_11 + int(4)) + make_float3 (0.94617468118667603f * z2_11 - 0.31539157032966614f) * *(coeffs_11 + int(5)) + make_float3 (fTmp0B_12 * x_15) * *(coeffs_11 + int(6)) + make_float3 (0.54627424478530884f * fC1_11) * *(coeffs_11 + int(7))) + (make_float3 (-0.59004360437393188f * (x_15 * fS1_11 + y_14 * fC1_11)) * *(coeffs_11 + int(8)) + make_float3 (fTmp1B_11 * fS1_11) * *(coeffs_11 + int(9)) + make_float3 (fTmp0C_11 * y_14) * *(coeffs_11 + int(10)) + make_float3 (z_12 * (1.86588168144226074f * z2_11 - 1.11952900886535645f)) * *(coeffs_11 + int(11)) + make_float3 (fTmp0C_11 * x_15) * *(coeffs_11 + int(12)) + make_float3 (fTmp1B_11 * fC1_11) * *(coeffs_11 + int(13)) + make_float3 (-0.59004360437393188f * (x_15 * fC1_11 - y_14 * fS1_11)) * *(coeffs_11 + int(14)));
}

inline __device__ float3  sh4_to_color(float3  viewdir_5, float3  coeff_dc_12, float3  * coeffs_12)
{
    float _S519 = viewdir_5.x;
    float _S520 = viewdir_5.y;
    float _S521 = viewdir_5.z;
    float inv_norm_8 = (F32_rsqrt((_S519 * _S519 + _S520 * _S520 + _S521 * _S521)));
    float x_16 = _S519 * inv_norm_8;
    float y_15 = _S520 * inv_norm_8;
    float z_13 = _S521 * inv_norm_8;
    float z2_12 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_12 = x_16 * x_16 - y_15 * y_15;
    float fS1_12 = 2.0f * x_16 * y_15;
    float pSH6_11 = 0.94617468118667603f * z2_12 - 0.31539157032966614f;
    float fTmp0C_12 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_13;
    float fC2_11 = x_16 * fC1_12 - y_15 * fS1_12;
    float fS2_11 = x_16 * fS1_12 + y_15 * fC1_12;
    float pSH12_11 = z_13 * (1.86588168144226074f * z2_12 - 1.11952900886535645f);
    float fTmp0D_11 = z_13 * (-4.68332576751708984f * z2_12 + 2.00713968276977539f);
    float fTmp1C_11 = 3.31161141395568848f * z2_12 - 0.47308734059333801f;
    float fTmp2B_11 = -1.77013075351715088f * z_13;
    return make_float3 (0.282094806432724f) * coeff_dc_12 + make_float3 (0.48860251903533936f) * (make_float3 (- y_15) * *(coeffs_12 + int(0)) + make_float3 (z_13) * *(coeffs_12 + int(1)) - make_float3 (x_16) * *(coeffs_12 + int(2))) + (make_float3 (0.54627424478530884f * fS1_12) * *(coeffs_12 + int(3)) + make_float3 (fTmp0B_13 * y_15) * *(coeffs_12 + int(4)) + make_float3 (pSH6_11) * *(coeffs_12 + int(5)) + make_float3 (fTmp0B_13 * x_16) * *(coeffs_12 + int(6)) + make_float3 (0.54627424478530884f * fC1_12) * *(coeffs_12 + int(7))) + (make_float3 (-0.59004360437393188f * fS2_11) * *(coeffs_12 + int(8)) + make_float3 (fTmp1B_12 * fS1_12) * *(coeffs_12 + int(9)) + make_float3 (fTmp0C_12 * y_15) * *(coeffs_12 + int(10)) + make_float3 (pSH12_11) * *(coeffs_12 + int(11)) + make_float3 (fTmp0C_12 * x_16) * *(coeffs_12 + int(12)) + make_float3 (fTmp1B_12 * fC1_12) * *(coeffs_12 + int(13)) + make_float3 (-0.59004360437393188f * fC2_11) * *(coeffs_12 + int(14))) + (make_float3 (0.62583571672439575f * (x_16 * fS2_11 + y_15 * fC2_11)) * *(coeffs_12 + int(15)) + make_float3 (fTmp2B_11 * fS2_11) * *(coeffs_12 + int(16)) + make_float3 (fTmp1C_11 * fS1_12) * *(coeffs_12 + int(17)) + make_float3 (fTmp0D_11 * y_15) * *(coeffs_12 + int(18)) + make_float3 (1.9843134880065918f * z_13 * pSH12_11 - 1.00623059272766113f * pSH6_11) * *(coeffs_12 + int(19)) + make_float3 (fTmp0D_11 * x_16) * *(coeffs_12 + int(20)) + make_float3 (fTmp1C_11 * fC1_12) * *(coeffs_12 + int(21)) + make_float3 (fTmp2B_11 * fC2_11) * *(coeffs_12 + int(22)) + make_float3 (0.62583571672439575f * (x_16 * fC2_11 - y_15 * fS2_11)) * *(coeffs_12 + int(23)));
}

inline __device__ float3  sh0_to_color(float3  mean_4, Matrix<float, 3, 3>  R_4, float3  t_4, float3  coeff_dc_13, float3  * coeffs_13)
{
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_13 + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ float3  sh1_to_color(float3  mean_5, Matrix<float, 3, 3>  R_5, float3  t_5, float3  coeff_dc_14, float3  * coeffs_14)
{
    float3  _S522 = mean_5 + mul_0(transpose_0(R_5), t_5);
    float _S523 = _S522.x;
    float _S524 = _S522.y;
    float _S525 = _S522.z;
    float inv_norm_9 = (F32_rsqrt((_S523 * _S523 + _S524 * _S524 + _S525 * _S525)));
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_14 + make_float3 (0.48860251903533936f) * (make_float3 (- (_S524 * inv_norm_9)) * *(coeffs_14 + int(0)) + make_float3 (_S525 * inv_norm_9) * *(coeffs_14 + int(1)) - make_float3 (_S523 * inv_norm_9) * *(coeffs_14 + int(2))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ float3  sh2_to_color(float3  mean_6, Matrix<float, 3, 3>  R_6, float3  t_6, float3  coeff_dc_15, float3  * coeffs_15)
{
    float3  _S526 = mean_6 + mul_0(transpose_0(R_6), t_6);
    float _S527 = _S526.x;
    float _S528 = _S526.y;
    float _S529 = _S526.z;
    float inv_norm_10 = (F32_rsqrt((_S527 * _S527 + _S528 * _S528 + _S529 * _S529)));
    float x_17 = _S527 * inv_norm_10;
    float y_16 = _S528 * inv_norm_10;
    float z_14 = _S529 * inv_norm_10;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_15 + make_float3 (0.48860251903533936f) * (make_float3 (- y_16) * *(coeffs_15 + int(0)) + make_float3 (z_14) * *(coeffs_15 + int(1)) - make_float3 (x_17) * *(coeffs_15 + int(2))) + (make_float3 (0.54627424478530884f * (2.0f * x_17 * y_16)) * *(coeffs_15 + int(3)) + make_float3 (fTmp0B_14 * y_16) * *(coeffs_15 + int(4)) + make_float3 (0.94617468118667603f * (z_14 * z_14) - 0.31539157032966614f) * *(coeffs_15 + int(5)) + make_float3 (fTmp0B_14 * x_17) * *(coeffs_15 + int(6)) + make_float3 (0.54627424478530884f * (x_17 * x_17 - y_16 * y_16)) * *(coeffs_15 + int(7))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ float3  sh3_to_color(float3  mean_7, Matrix<float, 3, 3>  R_7, float3  t_7, float3  coeff_dc_16, float3  * coeffs_16)
{
    float3  _S530 = mean_7 + mul_0(transpose_0(R_7), t_7);
    float _S531 = _S530.x;
    float _S532 = _S530.y;
    float _S533 = _S530.z;
    float inv_norm_11 = (F32_rsqrt((_S531 * _S531 + _S532 * _S532 + _S533 * _S533)));
    float x_18 = _S531 * inv_norm_11;
    float y_17 = _S532 * inv_norm_11;
    float z_15 = _S533 * inv_norm_11;
    float z2_13 = z_15 * z_15;
    float fTmp0B_15 = -1.09254848957061768f * z_15;
    float fC1_13 = x_18 * x_18 - y_17 * y_17;
    float fS1_13 = 2.0f * x_18 * y_17;
    float fTmp0C_13 = -2.28522896766662598f * z2_13 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_15;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_16 + make_float3 (0.48860251903533936f) * (make_float3 (- y_17) * *(coeffs_16 + int(0)) + make_float3 (z_15) * *(coeffs_16 + int(1)) - make_float3 (x_18) * *(coeffs_16 + int(2))) + (make_float3 (0.54627424478530884f * fS1_13) * *(coeffs_16 + int(3)) + make_float3 (fTmp0B_15 * y_17) * *(coeffs_16 + int(4)) + make_float3 (0.94617468118667603f * z2_13 - 0.31539157032966614f) * *(coeffs_16 + int(5)) + make_float3 (fTmp0B_15 * x_18) * *(coeffs_16 + int(6)) + make_float3 (0.54627424478530884f * fC1_13) * *(coeffs_16 + int(7))) + (make_float3 (-0.59004360437393188f * (x_18 * fS1_13 + y_17 * fC1_13)) * *(coeffs_16 + int(8)) + make_float3 (fTmp1B_13 * fS1_13) * *(coeffs_16 + int(9)) + make_float3 (fTmp0C_13 * y_17) * *(coeffs_16 + int(10)) + make_float3 (z_15 * (1.86588168144226074f * z2_13 - 1.11952900886535645f)) * *(coeffs_16 + int(11)) + make_float3 (fTmp0C_13 * x_18) * *(coeffs_16 + int(12)) + make_float3 (fTmp1B_13 * fC1_13) * *(coeffs_16 + int(13)) + make_float3 (-0.59004360437393188f * (x_18 * fC1_13 - y_17 * fS1_13)) * *(coeffs_16 + int(14))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ float3  sh4_to_color(float3  mean_8, Matrix<float, 3, 3>  R_8, float3  t_8, float3  coeff_dc_17, float3  * coeffs_17)
{
    float3  _S534 = mean_8 + mul_0(transpose_0(R_8), t_8);
    float _S535 = _S534.x;
    float _S536 = _S534.y;
    float _S537 = _S534.z;
    float inv_norm_12 = (F32_rsqrt((_S535 * _S535 + _S536 * _S536 + _S537 * _S537)));
    float x_19 = _S535 * inv_norm_12;
    float y_18 = _S536 * inv_norm_12;
    float z_16 = _S537 * inv_norm_12;
    float z2_14 = z_16 * z_16;
    float fTmp0B_16 = -1.09254848957061768f * z_16;
    float fC1_14 = x_19 * x_19 - y_18 * y_18;
    float fS1_14 = 2.0f * x_19 * y_18;
    float pSH6_12 = 0.94617468118667603f * z2_14 - 0.31539157032966614f;
    float fTmp0C_14 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_16;
    float fC2_12 = x_19 * fC1_14 - y_18 * fS1_14;
    float fS2_12 = x_19 * fS1_14 + y_18 * fC1_14;
    float pSH12_12 = z_16 * (1.86588168144226074f * z2_14 - 1.11952900886535645f);
    float fTmp0D_12 = z_16 * (-4.68332576751708984f * z2_14 + 2.00713968276977539f);
    float fTmp1C_12 = 3.31161141395568848f * z2_14 - 0.47308734059333801f;
    float fTmp2B_12 = -1.77013075351715088f * z_16;
    return max_0(make_float3 (0.282094806432724f) * coeff_dc_17 + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * *(coeffs_17 + int(0)) + make_float3 (z_16) * *(coeffs_17 + int(1)) - make_float3 (x_19) * *(coeffs_17 + int(2))) + (make_float3 (0.54627424478530884f * fS1_14) * *(coeffs_17 + int(3)) + make_float3 (fTmp0B_16 * y_18) * *(coeffs_17 + int(4)) + make_float3 (pSH6_12) * *(coeffs_17 + int(5)) + make_float3 (fTmp0B_16 * x_19) * *(coeffs_17 + int(6)) + make_float3 (0.54627424478530884f * fC1_14) * *(coeffs_17 + int(7))) + (make_float3 (-0.59004360437393188f * fS2_12) * *(coeffs_17 + int(8)) + make_float3 (fTmp1B_14 * fS1_14) * *(coeffs_17 + int(9)) + make_float3 (fTmp0C_14 * y_18) * *(coeffs_17 + int(10)) + make_float3 (pSH12_12) * *(coeffs_17 + int(11)) + make_float3 (fTmp0C_14 * x_19) * *(coeffs_17 + int(12)) + make_float3 (fTmp1B_14 * fC1_14) * *(coeffs_17 + int(13)) + make_float3 (-0.59004360437393188f * fC2_12) * *(coeffs_17 + int(14))) + (make_float3 (0.62583571672439575f * (x_19 * fS2_12 + y_18 * fC2_12)) * *(coeffs_17 + int(15)) + make_float3 (fTmp2B_12 * fS2_12) * *(coeffs_17 + int(16)) + make_float3 (fTmp1C_12 * fS1_14) * *(coeffs_17 + int(17)) + make_float3 (fTmp0D_12 * y_18) * *(coeffs_17 + int(18)) + make_float3 (1.9843134880065918f * z_16 * pSH12_12 - 1.00623059272766113f * pSH6_12) * *(coeffs_17 + int(19)) + make_float3 (fTmp0D_12 * x_19) * *(coeffs_17 + int(20)) + make_float3 (fTmp1C_12 * fC1_14) * *(coeffs_17 + int(21)) + make_float3 (fTmp2B_12 * fC2_12) * *(coeffs_17 + int(22)) + make_float3 (0.62583571672439575f * (x_19 * fC2_12 - y_18 * fS2_12)) * *(coeffs_17 + int(23))) + make_float3 (0.5f), make_float3 (0.0f));
}

inline __device__ void sh0_to_color_vjp(float3  dir_3, float3  coeff_dc_18, float3  * coeffs_18, float3  v_colors_6, float3  * v_coeff_dc_6, float3  * v_coeffs_6, float3  * v_dir_3)
{
    *v_coeff_dc_6 = make_float3 (0.282094806432724f) * v_colors_6;
    int3  _S538 = make_int3 (int(0));
    float3  _S539 = make_float3 ((float)_S538.x, (float)_S538.y, (float)_S538.z);
    *v_dir_3 = _S539;
    return;
}

inline __device__ void sh1_to_color_vjp(float3  dir_4, float3  coeff_dc_19, float3  * coeffs_19, float3  v_colors_7, float3  * v_coeff_dc_7, float3  * v_coeffs_7, float3  * v_dir_4)
{
    *v_coeff_dc_7 = make_float3 (0.282094806432724f) * v_colors_7;
    float _S540 = dir_4.x;
    float _S541 = dir_4.y;
    float _S542 = dir_4.z;
    float inorm_6 = (F32_rsqrt((_S540 * _S540 + _S541 * _S541 + _S542 * _S542)));
    float x_20 = _S540 * inorm_6;
    float y_19 = _S541 * inorm_6;
    float z_17 = _S542 * inorm_6;
    *(v_coeffs_7 + int(0)) = make_float3 (-0.48860251903533936f * y_19) * v_colors_7;
    *(v_coeffs_7 + int(1)) = make_float3 (0.48860251903533936f * z_17) * v_colors_7;
    *(v_coeffs_7 + int(2)) = make_float3 (-0.48860251903533936f * x_20) * v_colors_7;
    float3  dir_n_24 = make_float3 (x_20, y_19, z_17);
    float3  v_dir_n_24 = make_float3 (-0.48860251903533936f * dot_0(*(coeffs_19 + int(2)), v_colors_7), -0.48860251903533936f * dot_0(*(coeffs_19 + int(0)), v_colors_7), 0.48860251903533936f * dot_0(*(coeffs_19 + int(1)), v_colors_7));
    *v_dir_4 = (v_dir_n_24 - make_float3 (dot_0(v_dir_n_24, dir_n_24)) * dir_n_24) * make_float3 (inorm_6);
    return;
}

inline __device__ void sh2_to_color_vjp(float3  dir_5, float3  coeff_dc_20, float3  * coeffs_20, float3  v_colors_8, float3  * v_coeff_dc_8, float3  * v_coeffs_8, float3  * v_dir_5)
{
    *v_coeff_dc_8 = make_float3 (0.282094806432724f) * v_colors_8;
    float _S543 = dir_5.x;
    float _S544 = dir_5.y;
    float _S545 = dir_5.z;
    float inorm_7 = (F32_rsqrt((_S543 * _S543 + _S544 * _S544 + _S545 * _S545)));
    float x_21 = _S543 * inorm_7;
    float y_20 = _S544 * inorm_7;
    float z_18 = _S545 * inorm_7;
    *(v_coeffs_8 + int(0)) = make_float3 (-0.48860251903533936f * y_20) * v_colors_8;
    *(v_coeffs_8 + int(1)) = make_float3 (0.48860251903533936f * z_18) * v_colors_8;
    *(v_coeffs_8 + int(2)) = make_float3 (-0.48860251903533936f * x_21) * v_colors_8;
    float _S546 = -0.48860251903533936f * dot_0(*(coeffs_20 + int(2)), v_colors_8);
    float _S547 = -0.48860251903533936f * dot_0(*(coeffs_20 + int(0)), v_colors_8);
    float _S548 = 0.48860251903533936f * dot_0(*(coeffs_20 + int(1)), v_colors_8);
    float fTmp0B_17 = -1.09254848957061768f * z_18;
    float _S549 = 2.0f * x_21;
    float pSH6_13 = 0.94617468118667603f * (z_18 * z_18) - 0.31539157032966614f;
    float pSH7_6 = fTmp0B_17 * x_21;
    float pSH5_6 = fTmp0B_17 * y_20;
    float pSH8_6 = 0.54627424478530884f * (x_21 * x_21 - y_20 * y_20);
    *(v_coeffs_8 + int(3)) = make_float3 (0.54627424478530884f * (_S549 * y_20)) * v_colors_8;
    *(v_coeffs_8 + int(4)) = make_float3 (pSH5_6) * v_colors_8;
    *(v_coeffs_8 + int(5)) = make_float3 (pSH6_13) * v_colors_8;
    *(v_coeffs_8 + int(6)) = make_float3 (pSH7_6) * v_colors_8;
    *(v_coeffs_8 + int(7)) = make_float3 (pSH8_6) * v_colors_8;
    float pSH8_x_6 = 0.54627424478530884f * _S549;
    float3  * _S550 = coeffs_20 + int(3);
    float3  * _S551 = coeffs_20 + int(7);
    float3  * _S552 = coeffs_20 + int(6);
    float3  * _S553 = coeffs_20 + int(4);
    float3  dir_n_25 = make_float3 (x_21, y_20, z_18);
    float3  v_dir_n_25 = make_float3 (_S546 + dot_0(v_colors_8, make_float3 (0.54627424478530884f * (2.0f * y_20)) * *_S550 + make_float3 (pSH8_x_6) * *_S551 + make_float3 (fTmp0B_17) * *_S552), _S547 + dot_0(v_colors_8, make_float3 (pSH8_x_6) * *_S550 + make_float3 (0.54627424478530884f * (-2.0f * y_20)) * *_S551 + make_float3 (fTmp0B_17) * *_S553), _S548 + dot_0(v_colors_8, make_float3 (1.89234936237335205f * z_18) * *(coeffs_20 + int(5)) + make_float3 (-1.09254848957061768f * x_21) * *_S552 + make_float3 (-1.09254848957061768f * y_20) * *_S553));
    *v_dir_5 = (v_dir_n_25 - make_float3 (dot_0(v_dir_n_25, dir_n_25)) * dir_n_25) * make_float3 (inorm_7);
    return;
}

inline __device__ void sh3_to_color_vjp(float3  dir_6, float3  coeff_dc_21, float3  * coeffs_21, float3  v_colors_9, float3  * v_coeff_dc_9, float3  * v_coeffs_9, float3  * v_dir_6)
{
    *v_coeff_dc_9 = make_float3 (0.282094806432724f) * v_colors_9;
    float _S554 = dir_6.x;
    float _S555 = dir_6.y;
    float _S556 = dir_6.z;
    float inorm_8 = (F32_rsqrt((_S554 * _S554 + _S555 * _S555 + _S556 * _S556)));
    float x_22 = _S554 * inorm_8;
    float y_21 = _S555 * inorm_8;
    float z_19 = _S556 * inorm_8;
    *(v_coeffs_9 + int(0)) = make_float3 (-0.48860251903533936f * y_21) * v_colors_9;
    *(v_coeffs_9 + int(1)) = make_float3 (0.48860251903533936f * z_19) * v_colors_9;
    *(v_coeffs_9 + int(2)) = make_float3 (-0.48860251903533936f * x_22) * v_colors_9;
    float _S557 = -0.48860251903533936f * dot_0(*(coeffs_21 + int(2)), v_colors_9);
    float _S558 = -0.48860251903533936f * dot_0(*(coeffs_21 + int(0)), v_colors_9);
    float _S559 = 0.48860251903533936f * dot_0(*(coeffs_21 + int(1)), v_colors_9);
    float z2_15 = z_19 * z_19;
    float fTmp0B_18 = -1.09254848957061768f * z_19;
    float fC1_15 = x_22 * x_22 - y_21 * y_21;
    float _S560 = 2.0f * x_22;
    float fS1_15 = _S560 * y_21;
    float pSH6_14 = 0.94617468118667603f * z2_15 - 0.31539157032966614f;
    float pSH7_7 = fTmp0B_18 * x_22;
    float pSH5_7 = fTmp0B_18 * y_21;
    float pSH8_7 = 0.54627424478530884f * fC1_15;
    *(v_coeffs_9 + int(3)) = make_float3 (0.54627424478530884f * fS1_15) * v_colors_9;
    *(v_coeffs_9 + int(4)) = make_float3 (pSH5_7) * v_colors_9;
    *(v_coeffs_9 + int(5)) = make_float3 (pSH6_14) * v_colors_9;
    *(v_coeffs_9 + int(6)) = make_float3 (pSH7_7) * v_colors_9;
    *(v_coeffs_9 + int(7)) = make_float3 (pSH8_7) * v_colors_9;
    float fC1_y_6 = -2.0f * y_21;
    float fS1_x_6 = 2.0f * y_21;
    float pSH8_x_7 = 0.54627424478530884f * _S560;
    float3  * _S561 = coeffs_21 + int(3);
    float3  * _S562 = coeffs_21 + int(7);
    float3  * _S563 = coeffs_21 + int(6);
    float v_x_12 = _S557 + dot_0(v_colors_9, make_float3 (0.54627424478530884f * fS1_x_6) * *_S561 + make_float3 (pSH8_x_7) * *_S562 + make_float3 (fTmp0B_18) * *_S563);
    float3  * _S564 = coeffs_21 + int(4);
    float v_y_12 = _S558 + dot_0(v_colors_9, make_float3 (pSH8_x_7) * *_S561 + make_float3 (0.54627424478530884f * fC1_y_6) * *_S562 + make_float3 (fTmp0B_18) * *_S564);
    float v_z_12 = _S559 + dot_0(v_colors_9, make_float3 (1.89234936237335205f * z_19) * *(coeffs_21 + int(5)) + make_float3 (-1.09254848957061768f * x_22) * *_S563 + make_float3 (-1.09254848957061768f * y_21) * *_S564);
    float fTmp0C_15 = -2.28522896766662598f * z2_15 + 0.4570457935333252f;
    float fTmp1B_15 = 1.44530570507049561f * z_19;
    float pSH12_13 = z_19 * (1.86588168144226074f * z2_15 - 1.11952900886535645f);
    float pSH13_6 = fTmp0C_15 * x_22;
    float pSH11_6 = fTmp0C_15 * y_21;
    float pSH14_6 = fTmp1B_15 * fC1_15;
    float pSH10_6 = fTmp1B_15 * fS1_15;
    float pSH15_6 = -0.59004360437393188f * (x_22 * fC1_15 - y_21 * fS1_15);
    *(v_coeffs_9 + int(8)) = make_float3 (-0.59004360437393188f * (x_22 * fS1_15 + y_21 * fC1_15)) * v_colors_9;
    *(v_coeffs_9 + int(9)) = make_float3 (pSH10_6) * v_colors_9;
    *(v_coeffs_9 + int(10)) = make_float3 (pSH11_6) * v_colors_9;
    *(v_coeffs_9 + int(11)) = make_float3 (pSH12_13) * v_colors_9;
    *(v_coeffs_9 + int(12)) = make_float3 (pSH13_6) * v_colors_9;
    *(v_coeffs_9 + int(13)) = make_float3 (pSH14_6) * v_colors_9;
    *(v_coeffs_9 + int(14)) = make_float3 (pSH15_6) * v_colors_9;
    float fTmp0C_z_6 = -4.57045793533325195f * z_19;
    float _S565 = x_22 * _S560;
    float _S566 = y_21 * _S560;
    float pSH14_x_6 = fTmp1B_15 * _S560;
    float3  * _S567 = coeffs_21 + int(8);
    float3  * _S568 = coeffs_21 + int(14);
    float3  * _S569 = coeffs_21 + int(9);
    float3  * _S570 = coeffs_21 + int(13);
    float3  * _S571 = coeffs_21 + int(12);
    float3  * _S572 = coeffs_21 + int(10);
    float3  dir_n_26 = make_float3 (x_22, y_21, z_19);
    float3  v_dir_n_26 = make_float3 (v_x_12 + dot_0(v_colors_9, make_float3 (-0.59004360437393188f * (fS1_15 + x_22 * fS1_x_6 + _S566)) * *_S567 + make_float3 (-0.59004360437393188f * (fC1_15 + _S565 - y_21 * fS1_x_6)) * *_S568 + make_float3 (fTmp1B_15 * fS1_x_6) * *_S569 + make_float3 (pSH14_x_6) * *_S570 + make_float3 (fTmp0C_15) * *_S571), v_y_12 + dot_0(v_colors_9, make_float3 (-0.59004360437393188f * (_S565 + fC1_15 + y_21 * fC1_y_6)) * *_S567 + make_float3 (-0.59004360437393188f * (x_22 * fC1_y_6 - fS1_15 - _S566)) * *_S568 + make_float3 (pSH14_x_6) * *_S569 + make_float3 (fTmp1B_15 * fC1_y_6) * *_S570 + make_float3 (fTmp0C_15) * *_S572), v_z_12 + dot_0(v_colors_9, make_float3 (5.59764480590820312f * z2_15 - 1.11952900886535645f) * *(coeffs_21 + int(11)) + make_float3 (fTmp0C_z_6 * x_22) * *_S571 + make_float3 (fTmp0C_z_6 * y_21) * *_S572 + make_float3 (1.44530570507049561f * fC1_15) * *_S570 + make_float3 (1.44530570507049561f * fS1_15) * *_S569));
    *v_dir_6 = (v_dir_n_26 - make_float3 (dot_0(v_dir_n_26, dir_n_26)) * dir_n_26) * make_float3 (inorm_8);
    return;
}

inline __device__ void sh4_to_color_vjp(float3  dir_7, float3  coeff_dc_22, float3  * coeffs_22, float3  v_colors_10, float3  * v_coeff_dc_10, float3  * v_coeffs_10, float3  * v_dir_7)
{
    *v_coeff_dc_10 = make_float3 (0.282094806432724f) * v_colors_10;
    float _S573 = dir_7.x;
    float _S574 = dir_7.y;
    float _S575 = dir_7.z;
    float inorm_9 = (F32_rsqrt((_S573 * _S573 + _S574 * _S574 + _S575 * _S575)));
    float x_23 = _S573 * inorm_9;
    float y_22 = _S574 * inorm_9;
    float z_20 = _S575 * inorm_9;
    *(v_coeffs_10 + int(0)) = make_float3 (-0.48860251903533936f * y_22) * v_colors_10;
    *(v_coeffs_10 + int(1)) = make_float3 (0.48860251903533936f * z_20) * v_colors_10;
    *(v_coeffs_10 + int(2)) = make_float3 (-0.48860251903533936f * x_23) * v_colors_10;
    float _S576 = -0.48860251903533936f * dot_0(*(coeffs_22 + int(2)), v_colors_10);
    float _S577 = -0.48860251903533936f * dot_0(*(coeffs_22 + int(0)), v_colors_10);
    float _S578 = 0.48860251903533936f * dot_0(*(coeffs_22 + int(1)), v_colors_10);
    float z2_16 = z_20 * z_20;
    float fTmp0B_19 = -1.09254848957061768f * z_20;
    float fC1_16 = x_23 * x_23 - y_22 * y_22;
    float _S579 = 2.0f * x_23;
    float fS1_16 = _S579 * y_22;
    float pSH6_15 = 0.94617468118667603f * z2_16 - 0.31539157032966614f;
    float pSH7_8 = fTmp0B_19 * x_23;
    float pSH5_8 = fTmp0B_19 * y_22;
    float pSH8_8 = 0.54627424478530884f * fC1_16;
    *(v_coeffs_10 + int(3)) = make_float3 (0.54627424478530884f * fS1_16) * v_colors_10;
    *(v_coeffs_10 + int(4)) = make_float3 (pSH5_8) * v_colors_10;
    *(v_coeffs_10 + int(5)) = make_float3 (pSH6_15) * v_colors_10;
    *(v_coeffs_10 + int(6)) = make_float3 (pSH7_8) * v_colors_10;
    *(v_coeffs_10 + int(7)) = make_float3 (pSH8_8) * v_colors_10;
    float fC1_y_7 = -2.0f * y_22;
    float fS1_x_7 = 2.0f * y_22;
    float pSH6_z_6 = 1.89234936237335205f * z_20;
    float pSH8_x_8 = 0.54627424478530884f * _S579;
    float3  * _S580 = coeffs_22 + int(3);
    float3  * _S581 = coeffs_22 + int(7);
    float3  * _S582 = coeffs_22 + int(6);
    float v_x_13 = _S576 + dot_0(v_colors_10, make_float3 (0.54627424478530884f * fS1_x_7) * *_S580 + make_float3 (pSH8_x_8) * *_S581 + make_float3 (fTmp0B_19) * *_S582);
    float3  * _S583 = coeffs_22 + int(4);
    float v_y_13 = _S577 + dot_0(v_colors_10, make_float3 (pSH8_x_8) * *_S580 + make_float3 (0.54627424478530884f * fC1_y_7) * *_S581 + make_float3 (fTmp0B_19) * *_S583);
    float v_z_13 = _S578 + dot_0(v_colors_10, make_float3 (pSH6_z_6) * *(coeffs_22 + int(5)) + make_float3 (-1.09254848957061768f * x_23) * *_S582 + make_float3 (-1.09254848957061768f * y_22) * *_S583);
    float fTmp0C_16 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_16 = 1.44530570507049561f * z_20;
    float fC2_13 = x_23 * fC1_16 - y_22 * fS1_16;
    float fS2_13 = x_23 * fS1_16 + y_22 * fC1_16;
    float pSH12_14 = z_20 * (1.86588168144226074f * z2_16 - 1.11952900886535645f);
    float pSH13_7 = fTmp0C_16 * x_23;
    float pSH11_7 = fTmp0C_16 * y_22;
    float pSH14_7 = fTmp1B_16 * fC1_16;
    float pSH10_7 = fTmp1B_16 * fS1_16;
    float pSH15_7 = -0.59004360437393188f * fC2_13;
    *(v_coeffs_10 + int(8)) = make_float3 (-0.59004360437393188f * fS2_13) * v_colors_10;
    *(v_coeffs_10 + int(9)) = make_float3 (pSH10_7) * v_colors_10;
    *(v_coeffs_10 + int(10)) = make_float3 (pSH11_7) * v_colors_10;
    *(v_coeffs_10 + int(11)) = make_float3 (pSH12_14) * v_colors_10;
    *(v_coeffs_10 + int(12)) = make_float3 (pSH13_7) * v_colors_10;
    *(v_coeffs_10 + int(13)) = make_float3 (pSH14_7) * v_colors_10;
    *(v_coeffs_10 + int(14)) = make_float3 (pSH15_7) * v_colors_10;
    float fTmp0C_z_7 = -4.57045793533325195f * z_20;
    float _S584 = x_23 * _S579;
    float fC2_x_6 = fC1_16 + _S584 - y_22 * fS1_x_7;
    float _S585 = y_22 * _S579;
    float fC2_y_6 = x_23 * fC1_y_7 - fS1_16 - _S585;
    float fS2_x_6 = fS1_16 + x_23 * fS1_x_7 + _S585;
    float fS2_y_6 = _S584 + fC1_16 + y_22 * fC1_y_7;
    float pSH12_z_6 = 5.59764480590820312f * z2_16 - 1.11952900886535645f;
    float pSH14_x_7 = fTmp1B_16 * _S579;
    float3  * _S586 = coeffs_22 + int(8);
    float3  * _S587 = coeffs_22 + int(14);
    float3  * _S588 = coeffs_22 + int(9);
    float3  * _S589 = coeffs_22 + int(13);
    float3  * _S590 = coeffs_22 + int(12);
    float v_x_14 = v_x_13 + dot_0(v_colors_10, make_float3 (-0.59004360437393188f * fS2_x_6) * *_S586 + make_float3 (-0.59004360437393188f * fC2_x_6) * *_S587 + make_float3 (fTmp1B_16 * fS1_x_7) * *_S588 + make_float3 (pSH14_x_7) * *_S589 + make_float3 (fTmp0C_16) * *_S590);
    float3  * _S591 = coeffs_22 + int(10);
    float v_y_14 = v_y_13 + dot_0(v_colors_10, make_float3 (-0.59004360437393188f * fS2_y_6) * *_S586 + make_float3 (-0.59004360437393188f * fC2_y_6) * *_S587 + make_float3 (pSH14_x_7) * *_S588 + make_float3 (fTmp1B_16 * fC1_y_7) * *_S589 + make_float3 (fTmp0C_16) * *_S591);
    float v_z_14 = v_z_13 + dot_0(v_colors_10, make_float3 (pSH12_z_6) * *(coeffs_22 + int(11)) + make_float3 (fTmp0C_z_7 * x_23) * *_S590 + make_float3 (fTmp0C_z_7 * y_22) * *_S591 + make_float3 (1.44530570507049561f * fC1_16) * *_S589 + make_float3 (1.44530570507049561f * fS1_16) * *_S588);
    float fTmp0D_13 = z_20 * (-4.68332576751708984f * z2_16 + 2.00713968276977539f);
    float fTmp1C_13 = 3.31161141395568848f * z2_16 - 0.47308734059333801f;
    float fTmp2B_13 = -1.77013075351715088f * z_20;
    float pSH20_6 = 1.9843134880065918f * z_20 * pSH12_14 + -1.00623059272766113f * pSH6_15;
    float pSH21_6 = fTmp0D_13 * x_23;
    float pSH19_6 = fTmp0D_13 * y_22;
    float pSH22_6 = fTmp1C_13 * fC1_16;
    float pSH18_6 = fTmp1C_13 * fS1_16;
    float pSH23_6 = fTmp2B_13 * fC2_13;
    float pSH17_6 = fTmp2B_13 * fS2_13;
    float pSH24_6 = 0.62583571672439575f * (x_23 * fC2_13 - y_22 * fS2_13);
    *(v_coeffs_10 + int(15)) = make_float3 (0.62583571672439575f * (x_23 * fS2_13 + y_22 * fC2_13)) * v_colors_10;
    *(v_coeffs_10 + int(16)) = make_float3 (pSH17_6) * v_colors_10;
    *(v_coeffs_10 + int(17)) = make_float3 (pSH18_6) * v_colors_10;
    *(v_coeffs_10 + int(18)) = make_float3 (pSH19_6) * v_colors_10;
    *(v_coeffs_10 + int(19)) = make_float3 (pSH20_6) * v_colors_10;
    *(v_coeffs_10 + int(20)) = make_float3 (pSH21_6) * v_colors_10;
    *(v_coeffs_10 + int(21)) = make_float3 (pSH22_6) * v_colors_10;
    *(v_coeffs_10 + int(22)) = make_float3 (pSH23_6) * v_colors_10;
    *(v_coeffs_10 + int(23)) = make_float3 (pSH24_6) * v_colors_10;
    float fTmp0D_z_6 = -14.04997730255126953f * z2_16 + 2.00713968276977539f;
    float fTmp1C_z_6 = 6.62322282791137695f * z_20;
    float pSH22_x_6 = fTmp1C_13 * _S579;
    float3  * _S592 = coeffs_22 + int(15);
    float3  * _S593 = coeffs_22 + int(23);
    float3  * _S594 = coeffs_22 + int(16);
    float3  * _S595 = coeffs_22 + int(22);
    float3  * _S596 = coeffs_22 + int(17);
    float3  * _S597 = coeffs_22 + int(21);
    float3  * _S598 = coeffs_22 + int(20);
    float3  * _S599 = coeffs_22 + int(18);
    float3  dir_n_27 = make_float3 (x_23, y_22, z_20);
    float3  v_dir_n_27 = make_float3 (v_x_14 + dot_0(v_colors_10, make_float3 (0.62583571672439575f * (fS2_13 + y_22 * fC2_x_6 + x_23 * fS2_x_6)) * *_S592 + make_float3 (0.62583571672439575f * (fC2_13 + x_23 * fC2_x_6 - y_22 * fS2_x_6)) * *_S593 + make_float3 (fTmp2B_13 * fS2_x_6) * *_S594 + make_float3 (fTmp2B_13 * fC2_x_6) * *_S595 + make_float3 (fTmp1C_13 * fS1_x_7) * *_S596 + make_float3 (pSH22_x_6) * *_S597 + make_float3 (fTmp0D_13) * *_S598), v_y_14 + dot_0(v_colors_10, make_float3 (0.62583571672439575f * (x_23 * fS2_y_6 + fC2_13 + y_22 * fC2_y_6)) * *_S592 + make_float3 (0.62583571672439575f * (x_23 * fC2_y_6 - fS2_13 - y_22 * fS2_y_6)) * *_S593 + make_float3 (fTmp2B_13 * fS2_y_6) * *_S594 + make_float3 (fTmp2B_13 * fC2_y_6) * *_S595 + make_float3 (pSH22_x_6) * *_S596 + make_float3 (fTmp1C_13 * fC1_y_7) * *_S597 + make_float3 (fTmp0D_13) * *_S599), v_z_14 + dot_0(v_colors_10, make_float3 (1.9843134880065918f * (pSH12_14 + z_20 * pSH12_z_6) + -1.00623059272766113f * pSH6_z_6) * *(coeffs_22 + int(19)) + make_float3 (fTmp0D_z_6 * x_23) * *_S598 + make_float3 (fTmp0D_z_6 * y_22) * *_S599 + make_float3 (fTmp1C_z_6 * fC1_16) * *_S597 + make_float3 (fTmp1C_z_6 * fS1_16) * *_S596 + make_float3 (-1.77013075351715088f * fC2_13) * *_S595 + make_float3 (-1.77013075351715088f * fS2_13) * *_S594));
    *v_dir_7 = (v_dir_n_27 - make_float3 (dot_0(v_dir_n_27, dir_n_27)) * dir_n_27) * make_float3 (inorm_9);
    return;
}

inline __device__ void sh0_to_color_vjp(float3  mean_9, Matrix<float, 3, 3>  R_9, float3  t_9, float3  coeff_dc_23, float3  * coeffs_23, float3  v_colors_11, float3  * v_coeff_dc_11, float3  * v_coeffs_11, float3  * v_mean_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  colors_20 = make_float3 (0.282094806432724f) * coeff_dc_23;
    *v_coeff_dc_11 = make_float3 (0.282094806432724f) * (v_colors_11 * make_float3 (float((colors_20.x) >= -0.5f), float((colors_20.y) >= -0.5f), float((colors_20.z) >= -0.5f)));
    int3  _S600 = make_int3 (int(0));
    float3  v_viewdir_5 = make_float3 ((float)_S600.x, (float)_S600.y, (float)_S600.z);
    Matrix<float, 3, 3>  _S601 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S602;
    (&_S602)->primal_0 = transpose_0(R_9);
    (&_S602)->differential_0 = _S601;
    float3  _S603 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S604;
    (&_S604)->primal_0 = t_9;
    (&_S604)->differential_0 = _S603;
    s_bwd_prop_mul_0(&_S602, &_S604, v_viewdir_5);
    Matrix<float, 3, 3>  _S605 = transpose_0(_S602.differential_0);
    *v_mean_3 = v_viewdir_5;
    *v_R_3 = _S605;
    *v_t_3 = _S604.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp(float3  mean_10, Matrix<float, 3, 3>  R_10, float3  t_10, float3  coeff_dc_24, float3  * coeffs_24, float3  v_colors_12, float3  * v_coeff_dc_12, float3  * v_coeffs_12, float3  * v_mean_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 3, 3>  _S606 = transpose_0(R_10);
    float3  _S607 = mean_10 + mul_0(_S606, t_10);
    float _S608 = _S607.x;
    float _S609 = _S607.y;
    float _S610 = _S607.z;
    float inv_norm_13 = (F32_rsqrt((_S608 * _S608 + _S609 * _S609 + _S610 * _S610)));
    float x_24 = _S608 * inv_norm_13;
    float y_23 = _S609 * inv_norm_13;
    float z_21 = _S610 * inv_norm_13;
    float3  * _S611 = coeffs_24 + int(0);
    float3  * _S612 = coeffs_24 + int(1);
    float3  * _S613 = coeffs_24 + int(2);
    float3  colors_21 = make_float3 (0.282094806432724f) * coeff_dc_24 + make_float3 (0.48860251903533936f) * (make_float3 (- y_23) * *_S611 + make_float3 (z_21) * *_S612 - make_float3 (x_24) * *_S613);
    float3  _S614 = v_colors_12 * make_float3 (float((colors_21.x) >= -0.5f), float((colors_21.y) >= -0.5f), float((colors_21.z) >= -0.5f));
    *v_coeff_dc_12 = make_float3 (0.282094806432724f) * _S614;
    *(v_coeffs_12 + int(0)) = make_float3 (-0.48860251903533936f * y_23) * _S614;
    *(v_coeffs_12 + int(1)) = make_float3 (0.48860251903533936f * z_21) * _S614;
    *(v_coeffs_12 + int(2)) = make_float3 (-0.48860251903533936f * x_24) * _S614;
    float3  dir_n_28 = make_float3 (x_24, y_23, z_21);
    float3  v_dir_n_28 = make_float3 (-0.48860251903533936f * dot_0(*_S613, _S614), -0.48860251903533936f * dot_0(*_S611, _S614), 0.48860251903533936f * dot_0(*_S612, _S614));
    float3  v_viewdir_6 = (v_dir_n_28 - make_float3 (dot_0(v_dir_n_28, dir_n_28)) * dir_n_28) * make_float3 (inv_norm_13);
    Matrix<float, 3, 3>  _S615 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S616;
    (&_S616)->primal_0 = _S606;
    (&_S616)->differential_0 = _S615;
    float3  _S617 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S618;
    (&_S618)->primal_0 = t_10;
    (&_S618)->differential_0 = _S617;
    s_bwd_prop_mul_0(&_S616, &_S618, v_viewdir_6);
    Matrix<float, 3, 3>  _S619 = transpose_0(_S616.differential_0);
    *v_mean_4 = v_viewdir_6;
    *v_R_4 = _S619;
    *v_t_4 = _S618.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp(float3  mean_11, Matrix<float, 3, 3>  R_11, float3  t_11, float3  coeff_dc_25, float3  * coeffs_25, float3  v_colors_13, float3  * v_coeff_dc_13, float3  * v_coeffs_13, float3  * v_mean_5, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 3, 3>  _S620 = transpose_0(R_11);
    float3  _S621 = mean_11 + mul_0(_S620, t_11);
    float _S622 = _S621.x;
    float _S623 = _S621.y;
    float _S624 = _S621.z;
    float inv_norm_14 = (F32_rsqrt((_S622 * _S622 + _S623 * _S623 + _S624 * _S624)));
    float x_25 = _S622 * inv_norm_14;
    float y_24 = _S623 * inv_norm_14;
    float z_22 = _S624 * inv_norm_14;
    float3  * _S625 = coeffs_25 + int(0);
    float3  * _S626 = coeffs_25 + int(1);
    float3  * _S627 = coeffs_25 + int(2);
    float fTmp0B_20 = -1.09254848957061768f * z_22;
    float _S628 = 2.0f * x_25;
    float pSH6_16 = 0.94617468118667603f * (z_22 * z_22) - 0.31539157032966614f;
    float pSH7_9 = fTmp0B_20 * x_25;
    float pSH5_9 = fTmp0B_20 * y_24;
    float pSH8_9 = 0.54627424478530884f * (x_25 * x_25 - y_24 * y_24);
    float pSH4_6 = 0.54627424478530884f * (_S628 * y_24);
    float3  * _S629 = coeffs_25 + int(3);
    float3  * _S630 = coeffs_25 + int(4);
    float3  * _S631 = coeffs_25 + int(5);
    float3  * _S632 = coeffs_25 + int(6);
    float3  * _S633 = coeffs_25 + int(7);
    float3  colors_22 = make_float3 (0.282094806432724f) * coeff_dc_25 + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * *_S625 + make_float3 (z_22) * *_S626 - make_float3 (x_25) * *_S627) + (make_float3 (pSH4_6) * *_S629 + make_float3 (pSH5_9) * *_S630 + make_float3 (pSH6_16) * *_S631 + make_float3 (pSH7_9) * *_S632 + make_float3 (pSH8_9) * *_S633);
    float3  _S634 = v_colors_13 * make_float3 (float((colors_22.x) >= -0.5f), float((colors_22.y) >= -0.5f), float((colors_22.z) >= -0.5f));
    *v_coeff_dc_13 = make_float3 (0.282094806432724f) * _S634;
    *(v_coeffs_13 + int(0)) = make_float3 (-0.48860251903533936f * y_24) * _S634;
    *(v_coeffs_13 + int(1)) = make_float3 (0.48860251903533936f * z_22) * _S634;
    *(v_coeffs_13 + int(2)) = make_float3 (-0.48860251903533936f * x_25) * _S634;
    float _S635 = -0.48860251903533936f * dot_0(*_S627, _S634);
    float _S636 = -0.48860251903533936f * dot_0(*_S625, _S634);
    float _S637 = 0.48860251903533936f * dot_0(*_S626, _S634);
    *(v_coeffs_13 + int(3)) = make_float3 (pSH4_6) * _S634;
    *(v_coeffs_13 + int(4)) = make_float3 (pSH5_9) * _S634;
    *(v_coeffs_13 + int(5)) = make_float3 (pSH6_16) * _S634;
    *(v_coeffs_13 + int(6)) = make_float3 (pSH7_9) * _S634;
    *(v_coeffs_13 + int(7)) = make_float3 (pSH8_9) * _S634;
    float pSH8_x_9 = 0.54627424478530884f * _S628;
    float3  dir_n_29 = make_float3 (x_25, y_24, z_22);
    float3  v_dir_n_29 = make_float3 (_S635 + dot_0(_S634, make_float3 (0.54627424478530884f * (2.0f * y_24)) * *_S629 + make_float3 (pSH8_x_9) * *_S633 + make_float3 (fTmp0B_20) * *_S632), _S636 + dot_0(_S634, make_float3 (pSH8_x_9) * *_S629 + make_float3 (0.54627424478530884f * (-2.0f * y_24)) * *_S633 + make_float3 (fTmp0B_20) * *_S630), _S637 + dot_0(_S634, make_float3 (1.89234936237335205f * z_22) * *_S631 + make_float3 (-1.09254848957061768f * x_25) * *_S632 + make_float3 (-1.09254848957061768f * y_24) * *_S630));
    float3  v_viewdir_7 = (v_dir_n_29 - make_float3 (dot_0(v_dir_n_29, dir_n_29)) * dir_n_29) * make_float3 (inv_norm_14);
    Matrix<float, 3, 3>  _S638 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S639;
    (&_S639)->primal_0 = _S620;
    (&_S639)->differential_0 = _S638;
    float3  _S640 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S641;
    (&_S641)->primal_0 = t_11;
    (&_S641)->differential_0 = _S640;
    s_bwd_prop_mul_0(&_S639, &_S641, v_viewdir_7);
    Matrix<float, 3, 3>  _S642 = transpose_0(_S639.differential_0);
    *v_mean_5 = v_viewdir_7;
    *v_R_5 = _S642;
    *v_t_5 = _S641.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp(float3  mean_12, Matrix<float, 3, 3>  R_12, float3  t_12, float3  coeff_dc_26, float3  * coeffs_26, float3  v_colors_14, float3  * v_coeff_dc_14, float3  * v_coeffs_14, float3  * v_mean_6, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    Matrix<float, 3, 3>  _S643 = transpose_0(R_12);
    float3  _S644 = mean_12 + mul_0(_S643, t_12);
    float _S645 = _S644.x;
    float _S646 = _S644.y;
    float _S647 = _S644.z;
    float inv_norm_15 = (F32_rsqrt((_S645 * _S645 + _S646 * _S646 + _S647 * _S647)));
    float x_26 = _S645 * inv_norm_15;
    float y_25 = _S646 * inv_norm_15;
    float z_23 = _S647 * inv_norm_15;
    float3  * _S648 = coeffs_26 + int(0);
    float3  * _S649 = coeffs_26 + int(1);
    float3  * _S650 = coeffs_26 + int(2);
    float z2_17 = z_23 * z_23;
    float fTmp0B_21 = -1.09254848957061768f * z_23;
    float fC1_17 = x_26 * x_26 - y_25 * y_25;
    float _S651 = 2.0f * x_26;
    float fS1_17 = _S651 * y_25;
    float pSH6_17 = 0.94617468118667603f * z2_17 - 0.31539157032966614f;
    float pSH7_10 = fTmp0B_21 * x_26;
    float pSH5_10 = fTmp0B_21 * y_25;
    float pSH8_10 = 0.54627424478530884f * fC1_17;
    float pSH4_7 = 0.54627424478530884f * fS1_17;
    float3  * _S652 = coeffs_26 + int(3);
    float3  * _S653 = coeffs_26 + int(4);
    float3  * _S654 = coeffs_26 + int(5);
    float3  * _S655 = coeffs_26 + int(6);
    float3  * _S656 = coeffs_26 + int(7);
    float fTmp0C_17 = -2.28522896766662598f * z2_17 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_23;
    float pSH12_15 = z_23 * (1.86588168144226074f * z2_17 - 1.11952900886535645f);
    float pSH13_8 = fTmp0C_17 * x_26;
    float pSH11_8 = fTmp0C_17 * y_25;
    float pSH14_8 = fTmp1B_17 * fC1_17;
    float pSH10_8 = fTmp1B_17 * fS1_17;
    float pSH15_8 = -0.59004360437393188f * (x_26 * fC1_17 - y_25 * fS1_17);
    float pSH9_6 = -0.59004360437393188f * (x_26 * fS1_17 + y_25 * fC1_17);
    float3  * _S657 = coeffs_26 + int(8);
    float3  * _S658 = coeffs_26 + int(9);
    float3  * _S659 = coeffs_26 + int(10);
    float3  * _S660 = coeffs_26 + int(11);
    float3  * _S661 = coeffs_26 + int(12);
    float3  * _S662 = coeffs_26 + int(13);
    float3  * _S663 = coeffs_26 + int(14);
    float3  colors_23 = make_float3 (0.282094806432724f) * coeff_dc_26 + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * *_S648 + make_float3 (z_23) * *_S649 - make_float3 (x_26) * *_S650) + (make_float3 (pSH4_7) * *_S652 + make_float3 (pSH5_10) * *_S653 + make_float3 (pSH6_17) * *_S654 + make_float3 (pSH7_10) * *_S655 + make_float3 (pSH8_10) * *_S656) + (make_float3 (pSH9_6) * *_S657 + make_float3 (pSH10_8) * *_S658 + make_float3 (pSH11_8) * *_S659 + make_float3 (pSH12_15) * *_S660 + make_float3 (pSH13_8) * *_S661 + make_float3 (pSH14_8) * *_S662 + make_float3 (pSH15_8) * *_S663);
    float3  _S664 = v_colors_14 * make_float3 (float((colors_23.x) >= -0.5f), float((colors_23.y) >= -0.5f), float((colors_23.z) >= -0.5f));
    *v_coeff_dc_14 = make_float3 (0.282094806432724f) * _S664;
    *(v_coeffs_14 + int(0)) = make_float3 (-0.48860251903533936f * y_25) * _S664;
    *(v_coeffs_14 + int(1)) = make_float3 (0.48860251903533936f * z_23) * _S664;
    *(v_coeffs_14 + int(2)) = make_float3 (-0.48860251903533936f * x_26) * _S664;
    float _S665 = -0.48860251903533936f * dot_0(*_S650, _S664);
    float _S666 = -0.48860251903533936f * dot_0(*_S648, _S664);
    float _S667 = 0.48860251903533936f * dot_0(*_S649, _S664);
    *(v_coeffs_14 + int(3)) = make_float3 (pSH4_7) * _S664;
    *(v_coeffs_14 + int(4)) = make_float3 (pSH5_10) * _S664;
    *(v_coeffs_14 + int(5)) = make_float3 (pSH6_17) * _S664;
    *(v_coeffs_14 + int(6)) = make_float3 (pSH7_10) * _S664;
    *(v_coeffs_14 + int(7)) = make_float3 (pSH8_10) * _S664;
    float fC1_y_8 = -2.0f * y_25;
    float fS1_x_8 = 2.0f * y_25;
    float pSH8_x_10 = 0.54627424478530884f * _S651;
    float v_x_15 = _S665 + dot_0(_S664, make_float3 (0.54627424478530884f * fS1_x_8) * *_S652 + make_float3 (pSH8_x_10) * *_S656 + make_float3 (fTmp0B_21) * *_S655);
    float v_y_15 = _S666 + dot_0(_S664, make_float3 (pSH8_x_10) * *_S652 + make_float3 (0.54627424478530884f * fC1_y_8) * *_S656 + make_float3 (fTmp0B_21) * *_S653);
    float v_z_15 = _S667 + dot_0(_S664, make_float3 (1.89234936237335205f * z_23) * *_S654 + make_float3 (-1.09254848957061768f * x_26) * *_S655 + make_float3 (-1.09254848957061768f * y_25) * *_S653);
    *(v_coeffs_14 + int(8)) = make_float3 (pSH9_6) * _S664;
    *(v_coeffs_14 + int(9)) = make_float3 (pSH10_8) * _S664;
    *(v_coeffs_14 + int(10)) = make_float3 (pSH11_8) * _S664;
    *(v_coeffs_14 + int(11)) = make_float3 (pSH12_15) * _S664;
    *(v_coeffs_14 + int(12)) = make_float3 (pSH13_8) * _S664;
    *(v_coeffs_14 + int(13)) = make_float3 (pSH14_8) * _S664;
    *(v_coeffs_14 + int(14)) = make_float3 (pSH15_8) * _S664;
    float fTmp0C_z_8 = -4.57045793533325195f * z_23;
    float _S668 = x_26 * _S651;
    float _S669 = y_25 * _S651;
    float pSH14_x_8 = fTmp1B_17 * _S651;
    float3  dir_n_30 = make_float3 (x_26, y_25, z_23);
    float3  v_dir_n_30 = make_float3 (v_x_15 + dot_0(_S664, make_float3 (-0.59004360437393188f * (fS1_17 + x_26 * fS1_x_8 + _S669)) * *_S657 + make_float3 (-0.59004360437393188f * (fC1_17 + _S668 - y_25 * fS1_x_8)) * *_S663 + make_float3 (fTmp1B_17 * fS1_x_8) * *_S658 + make_float3 (pSH14_x_8) * *_S662 + make_float3 (fTmp0C_17) * *_S661), v_y_15 + dot_0(_S664, make_float3 (-0.59004360437393188f * (_S668 + fC1_17 + y_25 * fC1_y_8)) * *_S657 + make_float3 (-0.59004360437393188f * (x_26 * fC1_y_8 - fS1_17 - _S669)) * *_S663 + make_float3 (pSH14_x_8) * *_S658 + make_float3 (fTmp1B_17 * fC1_y_8) * *_S662 + make_float3 (fTmp0C_17) * *_S659), v_z_15 + dot_0(_S664, make_float3 (5.59764480590820312f * z2_17 - 1.11952900886535645f) * *_S660 + make_float3 (fTmp0C_z_8 * x_26) * *_S661 + make_float3 (fTmp0C_z_8 * y_25) * *_S659 + make_float3 (1.44530570507049561f * fC1_17) * *_S662 + make_float3 (1.44530570507049561f * fS1_17) * *_S658));
    float3  v_viewdir_8 = (v_dir_n_30 - make_float3 (dot_0(v_dir_n_30, dir_n_30)) * dir_n_30) * make_float3 (inv_norm_15);
    Matrix<float, 3, 3>  _S670 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S671;
    (&_S671)->primal_0 = _S643;
    (&_S671)->differential_0 = _S670;
    float3  _S672 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S673;
    (&_S673)->primal_0 = t_12;
    (&_S673)->differential_0 = _S672;
    s_bwd_prop_mul_0(&_S671, &_S673, v_viewdir_8);
    Matrix<float, 3, 3>  _S674 = transpose_0(_S671.differential_0);
    *v_mean_6 = v_viewdir_8;
    *v_R_6 = _S674;
    *v_t_6 = _S673.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp(float3  mean_13, Matrix<float, 3, 3>  R_13, float3  t_13, float3  coeff_dc_27, float3  * coeffs_27, float3  v_colors_15, float3  * v_coeff_dc_15, float3  * v_coeffs_15, float3  * v_mean_7, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    Matrix<float, 3, 3>  _S675 = transpose_0(R_13);
    float3  _S676 = mean_13 + mul_0(_S675, t_13);
    float _S677 = _S676.x;
    float _S678 = _S676.y;
    float _S679 = _S676.z;
    float inv_norm_16 = (F32_rsqrt((_S677 * _S677 + _S678 * _S678 + _S679 * _S679)));
    float x_27 = _S677 * inv_norm_16;
    float y_26 = _S678 * inv_norm_16;
    float z_24 = _S679 * inv_norm_16;
    float3  * _S680 = coeffs_27 + int(0);
    float3  * _S681 = coeffs_27 + int(1);
    float3  * _S682 = coeffs_27 + int(2);
    float z2_18 = z_24 * z_24;
    float fTmp0B_22 = -1.09254848957061768f * z_24;
    float fC1_18 = x_27 * x_27 - y_26 * y_26;
    float _S683 = 2.0f * x_27;
    float fS1_18 = _S683 * y_26;
    float pSH6_18 = 0.94617468118667603f * z2_18 - 0.31539157032966614f;
    float pSH7_11 = fTmp0B_22 * x_27;
    float pSH5_11 = fTmp0B_22 * y_26;
    float pSH8_11 = 0.54627424478530884f * fC1_18;
    float pSH4_8 = 0.54627424478530884f * fS1_18;
    float3  * _S684 = coeffs_27 + int(3);
    float3  * _S685 = coeffs_27 + int(4);
    float3  * _S686 = coeffs_27 + int(5);
    float3  * _S687 = coeffs_27 + int(6);
    float3  * _S688 = coeffs_27 + int(7);
    float fTmp0C_18 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_24;
    float fC2_14 = x_27 * fC1_18 - y_26 * fS1_18;
    float fS2_14 = x_27 * fS1_18 + y_26 * fC1_18;
    float pSH12_16 = z_24 * (1.86588168144226074f * z2_18 - 1.11952900886535645f);
    float pSH13_9 = fTmp0C_18 * x_27;
    float pSH11_9 = fTmp0C_18 * y_26;
    float pSH14_9 = fTmp1B_18 * fC1_18;
    float pSH10_9 = fTmp1B_18 * fS1_18;
    float pSH15_9 = -0.59004360437393188f * fC2_14;
    float pSH9_7 = -0.59004360437393188f * fS2_14;
    float3  * _S689 = coeffs_27 + int(8);
    float3  * _S690 = coeffs_27 + int(9);
    float3  * _S691 = coeffs_27 + int(10);
    float3  * _S692 = coeffs_27 + int(11);
    float3  * _S693 = coeffs_27 + int(12);
    float3  * _S694 = coeffs_27 + int(13);
    float3  * _S695 = coeffs_27 + int(14);
    float fTmp0D_14 = z_24 * (-4.68332576751708984f * z2_18 + 2.00713968276977539f);
    float fTmp1C_14 = 3.31161141395568848f * z2_18 - 0.47308734059333801f;
    float fTmp2B_14 = -1.77013075351715088f * z_24;
    float _S696 = 1.9843134880065918f * z_24 * pSH12_16;
    float pSH21_7 = fTmp0D_14 * x_27;
    float pSH19_7 = fTmp0D_14 * y_26;
    float pSH22_7 = fTmp1C_14 * fC1_18;
    float pSH18_7 = fTmp1C_14 * fS1_18;
    float pSH23_7 = fTmp2B_14 * fC2_14;
    float pSH17_7 = fTmp2B_14 * fS2_14;
    float pSH24_7 = 0.62583571672439575f * (x_27 * fC2_14 - y_26 * fS2_14);
    float pSH16_6 = 0.62583571672439575f * (x_27 * fS2_14 + y_26 * fC2_14);
    float3  * _S697 = coeffs_27 + int(15);
    float3  * _S698 = coeffs_27 + int(16);
    float3  * _S699 = coeffs_27 + int(17);
    float3  * _S700 = coeffs_27 + int(18);
    float3  * _S701 = coeffs_27 + int(19);
    float3  * _S702 = coeffs_27 + int(20);
    float3  * _S703 = coeffs_27 + int(21);
    float3  * _S704 = coeffs_27 + int(22);
    float3  * _S705 = coeffs_27 + int(23);
    float3  colors_24 = make_float3 (0.282094806432724f) * coeff_dc_27 + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * *_S680 + make_float3 (z_24) * *_S681 - make_float3 (x_27) * *_S682) + (make_float3 (pSH4_8) * *_S684 + make_float3 (pSH5_11) * *_S685 + make_float3 (pSH6_18) * *_S686 + make_float3 (pSH7_11) * *_S687 + make_float3 (pSH8_11) * *_S688) + (make_float3 (pSH9_7) * *_S689 + make_float3 (pSH10_9) * *_S690 + make_float3 (pSH11_9) * *_S691 + make_float3 (pSH12_16) * *_S692 + make_float3 (pSH13_9) * *_S693 + make_float3 (pSH14_9) * *_S694 + make_float3 (pSH15_9) * *_S695) + (make_float3 (pSH16_6) * *_S697 + make_float3 (pSH17_7) * *_S698 + make_float3 (pSH18_7) * *_S699 + make_float3 (pSH19_7) * *_S700 + make_float3 (_S696 - 1.00623059272766113f * pSH6_18) * *_S701 + make_float3 (pSH21_7) * *_S702 + make_float3 (pSH22_7) * *_S703 + make_float3 (pSH23_7) * *_S704 + make_float3 (pSH24_7) * *_S705);
    float3  _S706 = v_colors_15 * make_float3 (float((colors_24.x) >= -0.5f), float((colors_24.y) >= -0.5f), float((colors_24.z) >= -0.5f));
    *v_coeff_dc_15 = make_float3 (0.282094806432724f) * _S706;
    *(v_coeffs_15 + int(0)) = make_float3 (-0.48860251903533936f * y_26) * _S706;
    *(v_coeffs_15 + int(1)) = make_float3 (0.48860251903533936f * z_24) * _S706;
    *(v_coeffs_15 + int(2)) = make_float3 (-0.48860251903533936f * x_27) * _S706;
    float _S707 = -0.48860251903533936f * dot_0(*_S682, _S706);
    float _S708 = -0.48860251903533936f * dot_0(*_S680, _S706);
    float _S709 = 0.48860251903533936f * dot_0(*_S681, _S706);
    *(v_coeffs_15 + int(3)) = make_float3 (pSH4_8) * _S706;
    *(v_coeffs_15 + int(4)) = make_float3 (pSH5_11) * _S706;
    *(v_coeffs_15 + int(5)) = make_float3 (pSH6_18) * _S706;
    *(v_coeffs_15 + int(6)) = make_float3 (pSH7_11) * _S706;
    *(v_coeffs_15 + int(7)) = make_float3 (pSH8_11) * _S706;
    float fC1_y_9 = -2.0f * y_26;
    float fS1_x_9 = 2.0f * y_26;
    float pSH6_z_7 = 1.89234936237335205f * z_24;
    float pSH8_x_11 = 0.54627424478530884f * _S683;
    float v_x_16 = _S707 + dot_0(_S706, make_float3 (0.54627424478530884f * fS1_x_9) * *_S684 + make_float3 (pSH8_x_11) * *_S688 + make_float3 (fTmp0B_22) * *_S687);
    float v_y_16 = _S708 + dot_0(_S706, make_float3 (pSH8_x_11) * *_S684 + make_float3 (0.54627424478530884f * fC1_y_9) * *_S688 + make_float3 (fTmp0B_22) * *_S685);
    float v_z_16 = _S709 + dot_0(_S706, make_float3 (pSH6_z_7) * *_S686 + make_float3 (-1.09254848957061768f * x_27) * *_S687 + make_float3 (-1.09254848957061768f * y_26) * *_S685);
    *(v_coeffs_15 + int(8)) = make_float3 (pSH9_7) * _S706;
    *(v_coeffs_15 + int(9)) = make_float3 (pSH10_9) * _S706;
    *(v_coeffs_15 + int(10)) = make_float3 (pSH11_9) * _S706;
    *(v_coeffs_15 + int(11)) = make_float3 (pSH12_16) * _S706;
    *(v_coeffs_15 + int(12)) = make_float3 (pSH13_9) * _S706;
    *(v_coeffs_15 + int(13)) = make_float3 (pSH14_9) * _S706;
    *(v_coeffs_15 + int(14)) = make_float3 (pSH15_9) * _S706;
    float fTmp0C_z_9 = -4.57045793533325195f * z_24;
    float _S710 = x_27 * _S683;
    float fC2_x_7 = fC1_18 + _S710 - y_26 * fS1_x_9;
    float _S711 = y_26 * _S683;
    float fC2_y_7 = x_27 * fC1_y_9 - fS1_18 - _S711;
    float fS2_x_7 = fS1_18 + x_27 * fS1_x_9 + _S711;
    float fS2_y_7 = _S710 + fC1_18 + y_26 * fC1_y_9;
    float pSH12_z_7 = 5.59764480590820312f * z2_18 - 1.11952900886535645f;
    float pSH14_x_9 = fTmp1B_18 * _S683;
    float v_x_17 = v_x_16 + dot_0(_S706, make_float3 (-0.59004360437393188f * fS2_x_7) * *_S689 + make_float3 (-0.59004360437393188f * fC2_x_7) * *_S695 + make_float3 (fTmp1B_18 * fS1_x_9) * *_S690 + make_float3 (pSH14_x_9) * *_S694 + make_float3 (fTmp0C_18) * *_S693);
    float v_y_17 = v_y_16 + dot_0(_S706, make_float3 (-0.59004360437393188f * fS2_y_7) * *_S689 + make_float3 (-0.59004360437393188f * fC2_y_7) * *_S695 + make_float3 (pSH14_x_9) * *_S690 + make_float3 (fTmp1B_18 * fC1_y_9) * *_S694 + make_float3 (fTmp0C_18) * *_S691);
    float v_z_17 = v_z_16 + dot_0(_S706, make_float3 (pSH12_z_7) * *_S692 + make_float3 (fTmp0C_z_9 * x_27) * *_S693 + make_float3 (fTmp0C_z_9 * y_26) * *_S691 + make_float3 (1.44530570507049561f * fC1_18) * *_S694 + make_float3 (1.44530570507049561f * fS1_18) * *_S690);
    float pSH20_7 = _S696 + -1.00623059272766113f * pSH6_18;
    *(v_coeffs_15 + int(15)) = make_float3 (pSH16_6) * _S706;
    *(v_coeffs_15 + int(16)) = make_float3 (pSH17_7) * _S706;
    *(v_coeffs_15 + int(17)) = make_float3 (pSH18_7) * _S706;
    *(v_coeffs_15 + int(18)) = make_float3 (pSH19_7) * _S706;
    *(v_coeffs_15 + int(19)) = make_float3 (pSH20_7) * _S706;
    *(v_coeffs_15 + int(20)) = make_float3 (pSH21_7) * _S706;
    *(v_coeffs_15 + int(21)) = make_float3 (pSH22_7) * _S706;
    *(v_coeffs_15 + int(22)) = make_float3 (pSH23_7) * _S706;
    *(v_coeffs_15 + int(23)) = make_float3 (pSH24_7) * _S706;
    float fTmp0D_z_7 = -14.04997730255126953f * z2_18 + 2.00713968276977539f;
    float fTmp1C_z_7 = 6.62322282791137695f * z_24;
    float pSH22_x_7 = fTmp1C_14 * _S683;
    float3  dir_n_31 = make_float3 (x_27, y_26, z_24);
    float3  v_dir_n_31 = make_float3 (v_x_17 + dot_0(_S706, make_float3 (0.62583571672439575f * (fS2_14 + y_26 * fC2_x_7 + x_27 * fS2_x_7)) * *_S697 + make_float3 (0.62583571672439575f * (fC2_14 + x_27 * fC2_x_7 - y_26 * fS2_x_7)) * *_S705 + make_float3 (fTmp2B_14 * fS2_x_7) * *_S698 + make_float3 (fTmp2B_14 * fC2_x_7) * *_S704 + make_float3 (fTmp1C_14 * fS1_x_9) * *_S699 + make_float3 (pSH22_x_7) * *_S703 + make_float3 (fTmp0D_14) * *_S702), v_y_17 + dot_0(_S706, make_float3 (0.62583571672439575f * (x_27 * fS2_y_7 + fC2_14 + y_26 * fC2_y_7)) * *_S697 + make_float3 (0.62583571672439575f * (x_27 * fC2_y_7 - fS2_14 - y_26 * fS2_y_7)) * *_S705 + make_float3 (fTmp2B_14 * fS2_y_7) * *_S698 + make_float3 (fTmp2B_14 * fC2_y_7) * *_S704 + make_float3 (pSH22_x_7) * *_S699 + make_float3 (fTmp1C_14 * fC1_y_9) * *_S703 + make_float3 (fTmp0D_14) * *_S700), v_z_17 + dot_0(_S706, make_float3 (1.9843134880065918f * (pSH12_16 + z_24 * pSH12_z_7) + -1.00623059272766113f * pSH6_z_7) * *_S701 + make_float3 (fTmp0D_z_7 * x_27) * *_S702 + make_float3 (fTmp0D_z_7 * y_26) * *_S700 + make_float3 (fTmp1C_z_7 * fC1_18) * *_S703 + make_float3 (fTmp1C_z_7 * fS1_18) * *_S699 + make_float3 (-1.77013075351715088f * fC2_14) * *_S704 + make_float3 (-1.77013075351715088f * fS2_14) * *_S698));
    float3  v_viewdir_9 = (v_dir_n_31 - make_float3 (dot_0(v_dir_n_31, dir_n_31)) * dir_n_31) * make_float3 (inv_norm_16);
    Matrix<float, 3, 3>  _S712 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S713;
    (&_S713)->primal_0 = _S675;
    (&_S713)->differential_0 = _S712;
    float3  _S714 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S715;
    (&_S715)->primal_0 = t_13;
    (&_S715)->differential_0 = _S714;
    s_bwd_prop_mul_0(&_S713, &_S715, v_viewdir_9);
    Matrix<float, 3, 3>  _S716 = transpose_0(_S713.differential_0);
    *v_mean_7 = v_viewdir_9;
    *v_R_7 = _S716;
    *v_t_7 = _S715.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_inplace(float3  mean_14, Matrix<float, 3, 3>  R_14, float3  t_14, float3  coeff_dc_28, float3  * coeffs_28, float3  v_colors_16, float3  * v_coeff_dc_16, float3  * v_coeffs_16, float3  * v_mean_8, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  colors_25 = make_float3 (0.282094806432724f) * coeff_dc_28;
    *v_coeff_dc_16 = *v_coeff_dc_16 + make_float3 (0.282094806432724f) * (v_colors_16 * make_float3 (float((colors_25.x) >= -0.5f), float((colors_25.y) >= -0.5f), float((colors_25.z) >= -0.5f)));
    float3  v_viewdir_10 = {};
    Matrix<float, 3, 3>  _S717 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S718;
    (&_S718)->primal_0 = transpose_0(R_14);
    (&_S718)->differential_0 = _S717;
    float3  _S719 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S720;
    (&_S720)->primal_0 = t_14;
    (&_S720)->differential_0 = _S719;
    s_bwd_prop_mul_0(&_S718, &_S720, v_viewdir_10);
    Matrix<float, 3, 3>  _S721 = transpose_0(_S718.differential_0);
    *v_mean_8 = *v_mean_8 + v_viewdir_10;
    *v_R_8 = *v_R_8 + _S721;
    *v_t_8 = *v_t_8 + _S720.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_inplace(float3  mean_15, Matrix<float, 3, 3>  R_15, float3  t_15, float3  coeff_dc_29, float3  * coeffs_29, float3  v_colors_17, float3  * v_coeff_dc_17, float3  * v_coeffs_17, float3  * v_mean_9, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    Matrix<float, 3, 3>  _S722 = transpose_0(R_15);
    float3  _S723 = mean_15 + mul_0(_S722, t_15);
    float _S724 = _S723.x;
    float _S725 = _S723.y;
    float _S726 = _S723.z;
    float inv_norm_17 = (F32_rsqrt((_S724 * _S724 + _S725 * _S725 + _S726 * _S726)));
    float x_28 = _S724 * inv_norm_17;
    float y_27 = _S725 * inv_norm_17;
    float z_25 = _S726 * inv_norm_17;
    float3  * _S727 = coeffs_29 + int(0);
    float3  * _S728 = coeffs_29 + int(1);
    float3  * _S729 = coeffs_29 + int(2);
    float3  colors_26 = make_float3 (0.282094806432724f) * coeff_dc_29 + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * *_S727 + make_float3 (z_25) * *_S728 - make_float3 (x_28) * *_S729);
    float3  _S730 = v_colors_17 * make_float3 (float((colors_26.x) >= -0.5f), float((colors_26.y) >= -0.5f), float((colors_26.z) >= -0.5f));
    float3  v_viewdir_11 = {};
    *v_coeff_dc_17 = *v_coeff_dc_17 + make_float3 (0.282094806432724f) * _S730;
    float3  * _S731 = v_coeffs_17 + int(0);
    *_S731 = *_S731 + make_float3 (-0.48860251903533936f * y_27) * _S730;
    float3  * _S732 = v_coeffs_17 + int(1);
    *_S732 = *_S732 + make_float3 (0.48860251903533936f * z_25) * _S730;
    float3  * _S733 = v_coeffs_17 + int(2);
    *_S733 = *_S733 + make_float3 (-0.48860251903533936f * x_28) * _S730;
    float3  dir_n_32 = make_float3 (x_28, y_27, z_25);
    float3  v_dir_n_32 = make_float3 (-0.48860251903533936f * dot_0(*_S729, _S730), -0.48860251903533936f * dot_0(*_S727, _S730), 0.48860251903533936f * dot_0(*_S728, _S730));
    float3  v_viewdir_12 = v_viewdir_11 + (v_dir_n_32 - make_float3 (dot_0(v_dir_n_32, dir_n_32)) * dir_n_32) * make_float3 (inv_norm_17);
    Matrix<float, 3, 3>  _S734 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S735;
    (&_S735)->primal_0 = _S722;
    (&_S735)->differential_0 = _S734;
    float3  _S736 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S737;
    (&_S737)->primal_0 = t_15;
    (&_S737)->differential_0 = _S736;
    s_bwd_prop_mul_0(&_S735, &_S737, v_viewdir_12);
    Matrix<float, 3, 3>  _S738 = transpose_0(_S735.differential_0);
    *v_mean_9 = *v_mean_9 + v_viewdir_12;
    *v_R_9 = *v_R_9 + _S738;
    *v_t_9 = *v_t_9 + _S737.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_inplace(float3  mean_16, Matrix<float, 3, 3>  R_16, float3  t_16, float3  coeff_dc_30, float3  * coeffs_30, float3  v_colors_18, float3  * v_coeff_dc_18, float3  * v_coeffs_18, float3  * v_mean_10, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_10)
{
    Matrix<float, 3, 3>  _S739 = transpose_0(R_16);
    float3  _S740 = mean_16 + mul_0(_S739, t_16);
    float _S741 = _S740.x;
    float _S742 = _S740.y;
    float _S743 = _S740.z;
    float inv_norm_18 = (F32_rsqrt((_S741 * _S741 + _S742 * _S742 + _S743 * _S743)));
    float x_29 = _S741 * inv_norm_18;
    float y_28 = _S742 * inv_norm_18;
    float z_26 = _S743 * inv_norm_18;
    float3  * _S744 = coeffs_30 + int(0);
    float3  * _S745 = coeffs_30 + int(1);
    float3  * _S746 = coeffs_30 + int(2);
    float fTmp0B_23 = -1.09254848957061768f * z_26;
    float _S747 = 2.0f * x_29;
    float pSH6_19 = 0.94617468118667603f * (z_26 * z_26) - 0.31539157032966614f;
    float pSH7_12 = fTmp0B_23 * x_29;
    float pSH5_12 = fTmp0B_23 * y_28;
    float pSH8_12 = 0.54627424478530884f * (x_29 * x_29 - y_28 * y_28);
    float pSH4_9 = 0.54627424478530884f * (_S747 * y_28);
    float3  * _S748 = coeffs_30 + int(3);
    float3  * _S749 = coeffs_30 + int(4);
    float3  * _S750 = coeffs_30 + int(5);
    float3  * _S751 = coeffs_30 + int(6);
    float3  * _S752 = coeffs_30 + int(7);
    float3  colors_27 = make_float3 (0.282094806432724f) * coeff_dc_30 + make_float3 (0.48860251903533936f) * (make_float3 (- y_28) * *_S744 + make_float3 (z_26) * *_S745 - make_float3 (x_29) * *_S746) + (make_float3 (pSH4_9) * *_S748 + make_float3 (pSH5_12) * *_S749 + make_float3 (pSH6_19) * *_S750 + make_float3 (pSH7_12) * *_S751 + make_float3 (pSH8_12) * *_S752);
    float3  _S753 = v_colors_18 * make_float3 (float((colors_27.x) >= -0.5f), float((colors_27.y) >= -0.5f), float((colors_27.z) >= -0.5f));
    *v_coeff_dc_18 = *v_coeff_dc_18 + make_float3 (0.282094806432724f) * _S753;
    float3  v_viewdir_13 = {};
    float3  * _S754 = v_coeffs_18 + int(0);
    *_S754 = *_S754 + make_float3 (-0.48860251903533936f * y_28) * _S753;
    float3  * _S755 = v_coeffs_18 + int(1);
    *_S755 = *_S755 + make_float3 (0.48860251903533936f * z_26) * _S753;
    float3  * _S756 = v_coeffs_18 + int(2);
    *_S756 = *_S756 + make_float3 (-0.48860251903533936f * x_29) * _S753;
    float _S757 = -0.48860251903533936f * dot_0(*_S746, _S753);
    float _S758 = -0.48860251903533936f * dot_0(*_S744, _S753);
    float _S759 = 0.48860251903533936f * dot_0(*_S745, _S753);
    float3  * _S760 = v_coeffs_18 + int(3);
    *_S760 = *_S760 + make_float3 (pSH4_9) * _S753;
    float3  * _S761 = v_coeffs_18 + int(4);
    *_S761 = *_S761 + make_float3 (pSH5_12) * _S753;
    float3  * _S762 = v_coeffs_18 + int(5);
    *_S762 = *_S762 + make_float3 (pSH6_19) * _S753;
    float3  * _S763 = v_coeffs_18 + int(6);
    *_S763 = *_S763 + make_float3 (pSH7_12) * _S753;
    float3  * _S764 = v_coeffs_18 + int(7);
    *_S764 = *_S764 + make_float3 (pSH8_12) * _S753;
    float pSH8_x_12 = 0.54627424478530884f * _S747;
    float3  dir_n_33 = make_float3 (x_29, y_28, z_26);
    float3  v_dir_n_33 = make_float3 (_S757 + dot_0(_S753, make_float3 (0.54627424478530884f * (2.0f * y_28)) * *_S748 + make_float3 (pSH8_x_12) * *_S752 + make_float3 (fTmp0B_23) * *_S751), _S758 + dot_0(_S753, make_float3 (pSH8_x_12) * *_S748 + make_float3 (0.54627424478530884f * (-2.0f * y_28)) * *_S752 + make_float3 (fTmp0B_23) * *_S749), _S759 + dot_0(_S753, make_float3 (1.89234936237335205f * z_26) * *_S750 + make_float3 (-1.09254848957061768f * x_29) * *_S751 + make_float3 (-1.09254848957061768f * y_28) * *_S749));
    float3  v_viewdir_14 = v_viewdir_13 + (v_dir_n_33 - make_float3 (dot_0(v_dir_n_33, dir_n_33)) * dir_n_33) * make_float3 (inv_norm_18);
    Matrix<float, 3, 3>  _S765 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S766;
    (&_S766)->primal_0 = _S739;
    (&_S766)->differential_0 = _S765;
    float3  _S767 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S768;
    (&_S768)->primal_0 = t_16;
    (&_S768)->differential_0 = _S767;
    s_bwd_prop_mul_0(&_S766, &_S768, v_viewdir_14);
    Matrix<float, 3, 3>  _S769 = transpose_0(_S766.differential_0);
    *v_mean_10 = *v_mean_10 + v_viewdir_14;
    *v_R_10 = *v_R_10 + _S769;
    *v_t_10 = *v_t_10 + _S768.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_inplace(float3  mean_17, Matrix<float, 3, 3>  R_17, float3  t_17, float3  coeff_dc_31, float3  * coeffs_31, float3  v_colors_19, float3  * v_coeff_dc_19, float3  * v_coeffs_19, float3  * v_mean_11, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_11)
{
    Matrix<float, 3, 3>  _S770 = transpose_0(R_17);
    float3  _S771 = mean_17 + mul_0(_S770, t_17);
    float _S772 = _S771.x;
    float _S773 = _S771.y;
    float _S774 = _S771.z;
    float inv_norm_19 = (F32_rsqrt((_S772 * _S772 + _S773 * _S773 + _S774 * _S774)));
    float x_30 = _S772 * inv_norm_19;
    float y_29 = _S773 * inv_norm_19;
    float z_27 = _S774 * inv_norm_19;
    float3  * _S775 = coeffs_31 + int(0);
    float3  * _S776 = coeffs_31 + int(1);
    float3  * _S777 = coeffs_31 + int(2);
    float z2_19 = z_27 * z_27;
    float fTmp0B_24 = -1.09254848957061768f * z_27;
    float fC1_19 = x_30 * x_30 - y_29 * y_29;
    float _S778 = 2.0f * x_30;
    float fS1_19 = _S778 * y_29;
    float pSH6_20 = 0.94617468118667603f * z2_19 - 0.31539157032966614f;
    float pSH7_13 = fTmp0B_24 * x_30;
    float pSH5_13 = fTmp0B_24 * y_29;
    float pSH8_13 = 0.54627424478530884f * fC1_19;
    float pSH4_10 = 0.54627424478530884f * fS1_19;
    float3  * _S779 = coeffs_31 + int(3);
    float3  * _S780 = coeffs_31 + int(4);
    float3  * _S781 = coeffs_31 + int(5);
    float3  * _S782 = coeffs_31 + int(6);
    float3  * _S783 = coeffs_31 + int(7);
    float fTmp0C_19 = -2.28522896766662598f * z2_19 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_27;
    float pSH12_17 = z_27 * (1.86588168144226074f * z2_19 - 1.11952900886535645f);
    float pSH13_10 = fTmp0C_19 * x_30;
    float pSH11_10 = fTmp0C_19 * y_29;
    float pSH14_10 = fTmp1B_19 * fC1_19;
    float pSH10_10 = fTmp1B_19 * fS1_19;
    float pSH15_10 = -0.59004360437393188f * (x_30 * fC1_19 - y_29 * fS1_19);
    float pSH9_8 = -0.59004360437393188f * (x_30 * fS1_19 + y_29 * fC1_19);
    float3  * _S784 = coeffs_31 + int(8);
    float3  * _S785 = coeffs_31 + int(9);
    float3  * _S786 = coeffs_31 + int(10);
    float3  * _S787 = coeffs_31 + int(11);
    float3  * _S788 = coeffs_31 + int(12);
    float3  * _S789 = coeffs_31 + int(13);
    float3  * _S790 = coeffs_31 + int(14);
    float3  colors_28 = make_float3 (0.282094806432724f) * coeff_dc_31 + make_float3 (0.48860251903533936f) * (make_float3 (- y_29) * *_S775 + make_float3 (z_27) * *_S776 - make_float3 (x_30) * *_S777) + (make_float3 (pSH4_10) * *_S779 + make_float3 (pSH5_13) * *_S780 + make_float3 (pSH6_20) * *_S781 + make_float3 (pSH7_13) * *_S782 + make_float3 (pSH8_13) * *_S783) + (make_float3 (pSH9_8) * *_S784 + make_float3 (pSH10_10) * *_S785 + make_float3 (pSH11_10) * *_S786 + make_float3 (pSH12_17) * *_S787 + make_float3 (pSH13_10) * *_S788 + make_float3 (pSH14_10) * *_S789 + make_float3 (pSH15_10) * *_S790);
    float3  _S791 = v_colors_19 * make_float3 (float((colors_28.x) >= -0.5f), float((colors_28.y) >= -0.5f), float((colors_28.z) >= -0.5f));
    float3  v_viewdir_15 = {};
    *v_coeff_dc_19 = *v_coeff_dc_19 + make_float3 (0.282094806432724f) * _S791;
    float3  * _S792 = v_coeffs_19 + int(0);
    *_S792 = *_S792 + make_float3 (-0.48860251903533936f * y_29) * _S791;
    float3  * _S793 = v_coeffs_19 + int(1);
    *_S793 = *_S793 + make_float3 (0.48860251903533936f * z_27) * _S791;
    float3  * _S794 = v_coeffs_19 + int(2);
    *_S794 = *_S794 + make_float3 (-0.48860251903533936f * x_30) * _S791;
    float _S795 = -0.48860251903533936f * dot_0(*_S777, _S791);
    float _S796 = -0.48860251903533936f * dot_0(*_S775, _S791);
    float _S797 = 0.48860251903533936f * dot_0(*_S776, _S791);
    float3  * _S798 = v_coeffs_19 + int(3);
    *_S798 = *_S798 + make_float3 (pSH4_10) * _S791;
    float3  * _S799 = v_coeffs_19 + int(4);
    *_S799 = *_S799 + make_float3 (pSH5_13) * _S791;
    float3  * _S800 = v_coeffs_19 + int(5);
    *_S800 = *_S800 + make_float3 (pSH6_20) * _S791;
    float3  * _S801 = v_coeffs_19 + int(6);
    *_S801 = *_S801 + make_float3 (pSH7_13) * _S791;
    float3  * _S802 = v_coeffs_19 + int(7);
    *_S802 = *_S802 + make_float3 (pSH8_13) * _S791;
    float fC1_y_10 = -2.0f * y_29;
    float fS1_x_10 = 2.0f * y_29;
    float pSH8_x_13 = 0.54627424478530884f * _S778;
    float v_x_18 = _S795 + dot_0(_S791, make_float3 (0.54627424478530884f * fS1_x_10) * *_S779 + make_float3 (pSH8_x_13) * *_S783 + make_float3 (fTmp0B_24) * *_S782);
    float v_y_18 = _S796 + dot_0(_S791, make_float3 (pSH8_x_13) * *_S779 + make_float3 (0.54627424478530884f * fC1_y_10) * *_S783 + make_float3 (fTmp0B_24) * *_S780);
    float v_z_18 = _S797 + dot_0(_S791, make_float3 (1.89234936237335205f * z_27) * *_S781 + make_float3 (-1.09254848957061768f * x_30) * *_S782 + make_float3 (-1.09254848957061768f * y_29) * *_S780);
    float3  * _S803 = v_coeffs_19 + int(8);
    *_S803 = *_S803 + make_float3 (pSH9_8) * _S791;
    float3  * _S804 = v_coeffs_19 + int(9);
    *_S804 = *_S804 + make_float3 (pSH10_10) * _S791;
    float3  * _S805 = v_coeffs_19 + int(10);
    *_S805 = *_S805 + make_float3 (pSH11_10) * _S791;
    float3  * _S806 = v_coeffs_19 + int(11);
    *_S806 = *_S806 + make_float3 (pSH12_17) * _S791;
    float3  * _S807 = v_coeffs_19 + int(12);
    *_S807 = *_S807 + make_float3 (pSH13_10) * _S791;
    float3  * _S808 = v_coeffs_19 + int(13);
    *_S808 = *_S808 + make_float3 (pSH14_10) * _S791;
    float3  * _S809 = v_coeffs_19 + int(14);
    *_S809 = *_S809 + make_float3 (pSH15_10) * _S791;
    float fTmp0C_z_10 = -4.57045793533325195f * z_27;
    float _S810 = x_30 * _S778;
    float _S811 = y_29 * _S778;
    float pSH14_x_10 = fTmp1B_19 * _S778;
    float3  dir_n_34 = make_float3 (x_30, y_29, z_27);
    float3  v_dir_n_34 = make_float3 (v_x_18 + dot_0(_S791, make_float3 (-0.59004360437393188f * (fS1_19 + x_30 * fS1_x_10 + _S811)) * *_S784 + make_float3 (-0.59004360437393188f * (fC1_19 + _S810 - y_29 * fS1_x_10)) * *_S790 + make_float3 (fTmp1B_19 * fS1_x_10) * *_S785 + make_float3 (pSH14_x_10) * *_S789 + make_float3 (fTmp0C_19) * *_S788), v_y_18 + dot_0(_S791, make_float3 (-0.59004360437393188f * (_S810 + fC1_19 + y_29 * fC1_y_10)) * *_S784 + make_float3 (-0.59004360437393188f * (x_30 * fC1_y_10 - fS1_19 - _S811)) * *_S790 + make_float3 (pSH14_x_10) * *_S785 + make_float3 (fTmp1B_19 * fC1_y_10) * *_S789 + make_float3 (fTmp0C_19) * *_S786), v_z_18 + dot_0(_S791, make_float3 (5.59764480590820312f * z2_19 - 1.11952900886535645f) * *_S787 + make_float3 (fTmp0C_z_10 * x_30) * *_S788 + make_float3 (fTmp0C_z_10 * y_29) * *_S786 + make_float3 (1.44530570507049561f * fC1_19) * *_S789 + make_float3 (1.44530570507049561f * fS1_19) * *_S785));
    float3  v_viewdir_16 = v_viewdir_15 + (v_dir_n_34 - make_float3 (dot_0(v_dir_n_34, dir_n_34)) * dir_n_34) * make_float3 (inv_norm_19);
    Matrix<float, 3, 3>  _S812 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S813;
    (&_S813)->primal_0 = _S770;
    (&_S813)->differential_0 = _S812;
    float3  _S814 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S815;
    (&_S815)->primal_0 = t_17;
    (&_S815)->differential_0 = _S814;
    s_bwd_prop_mul_0(&_S813, &_S815, v_viewdir_16);
    Matrix<float, 3, 3>  _S816 = transpose_0(_S813.differential_0);
    *v_mean_11 = *v_mean_11 + v_viewdir_16;
    *v_R_11 = *v_R_11 + _S816;
    *v_t_11 = *v_t_11 + _S815.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_inplace(float3  mean_18, Matrix<float, 3, 3>  R_18, float3  t_18, float3  coeff_dc_32, float3  * coeffs_32, float3  v_colors_20, float3  * v_coeff_dc_20, float3  * v_coeffs_20, float3  * v_mean_12, Matrix<float, 3, 3>  * v_R_12, float3  * v_t_12)
{
    Matrix<float, 3, 3>  _S817 = transpose_0(R_18);
    float3  _S818 = mean_18 + mul_0(_S817, t_18);
    float _S819 = _S818.x;
    float _S820 = _S818.y;
    float _S821 = _S818.z;
    float inv_norm_20 = (F32_rsqrt((_S819 * _S819 + _S820 * _S820 + _S821 * _S821)));
    float x_31 = _S819 * inv_norm_20;
    float y_30 = _S820 * inv_norm_20;
    float z_28 = _S821 * inv_norm_20;
    float3  * _S822 = coeffs_32 + int(0);
    float3  * _S823 = coeffs_32 + int(1);
    float3  * _S824 = coeffs_32 + int(2);
    float z2_20 = z_28 * z_28;
    float fTmp0B_25 = -1.09254848957061768f * z_28;
    float fC1_20 = x_31 * x_31 - y_30 * y_30;
    float _S825 = 2.0f * x_31;
    float fS1_20 = _S825 * y_30;
    float pSH6_21 = 0.94617468118667603f * z2_20 - 0.31539157032966614f;
    float pSH7_14 = fTmp0B_25 * x_31;
    float pSH5_14 = fTmp0B_25 * y_30;
    float pSH8_14 = 0.54627424478530884f * fC1_20;
    float pSH4_11 = 0.54627424478530884f * fS1_20;
    float3  * _S826 = coeffs_32 + int(3);
    float3  * _S827 = coeffs_32 + int(4);
    float3  * _S828 = coeffs_32 + int(5);
    float3  * _S829 = coeffs_32 + int(6);
    float3  * _S830 = coeffs_32 + int(7);
    float fTmp0C_20 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_28;
    float fC2_15 = x_31 * fC1_20 - y_30 * fS1_20;
    float fS2_15 = x_31 * fS1_20 + y_30 * fC1_20;
    float pSH12_18 = z_28 * (1.86588168144226074f * z2_20 - 1.11952900886535645f);
    float pSH13_11 = fTmp0C_20 * x_31;
    float pSH11_11 = fTmp0C_20 * y_30;
    float pSH14_11 = fTmp1B_20 * fC1_20;
    float pSH10_11 = fTmp1B_20 * fS1_20;
    float pSH15_11 = -0.59004360437393188f * fC2_15;
    float pSH9_9 = -0.59004360437393188f * fS2_15;
    float3  * _S831 = coeffs_32 + int(8);
    float3  * _S832 = coeffs_32 + int(9);
    float3  * _S833 = coeffs_32 + int(10);
    float3  * _S834 = coeffs_32 + int(11);
    float3  * _S835 = coeffs_32 + int(12);
    float3  * _S836 = coeffs_32 + int(13);
    float3  * _S837 = coeffs_32 + int(14);
    float fTmp0D_15 = z_28 * (-4.68332576751708984f * z2_20 + 2.00713968276977539f);
    float fTmp1C_15 = 3.31161141395568848f * z2_20 - 0.47308734059333801f;
    float fTmp2B_15 = -1.77013075351715088f * z_28;
    float _S838 = 1.9843134880065918f * z_28 * pSH12_18;
    float pSH21_8 = fTmp0D_15 * x_31;
    float pSH19_8 = fTmp0D_15 * y_30;
    float pSH22_8 = fTmp1C_15 * fC1_20;
    float pSH18_8 = fTmp1C_15 * fS1_20;
    float pSH23_8 = fTmp2B_15 * fC2_15;
    float pSH17_8 = fTmp2B_15 * fS2_15;
    float pSH24_8 = 0.62583571672439575f * (x_31 * fC2_15 - y_30 * fS2_15);
    float pSH16_7 = 0.62583571672439575f * (x_31 * fS2_15 + y_30 * fC2_15);
    float3  * _S839 = coeffs_32 + int(15);
    float3  * _S840 = coeffs_32 + int(16);
    float3  * _S841 = coeffs_32 + int(17);
    float3  * _S842 = coeffs_32 + int(18);
    float3  * _S843 = coeffs_32 + int(19);
    float3  * _S844 = coeffs_32 + int(20);
    float3  * _S845 = coeffs_32 + int(21);
    float3  * _S846 = coeffs_32 + int(22);
    float3  * _S847 = coeffs_32 + int(23);
    float3  colors_29 = make_float3 (0.282094806432724f) * coeff_dc_32 + make_float3 (0.48860251903533936f) * (make_float3 (- y_30) * *_S822 + make_float3 (z_28) * *_S823 - make_float3 (x_31) * *_S824) + (make_float3 (pSH4_11) * *_S826 + make_float3 (pSH5_14) * *_S827 + make_float3 (pSH6_21) * *_S828 + make_float3 (pSH7_14) * *_S829 + make_float3 (pSH8_14) * *_S830) + (make_float3 (pSH9_9) * *_S831 + make_float3 (pSH10_11) * *_S832 + make_float3 (pSH11_11) * *_S833 + make_float3 (pSH12_18) * *_S834 + make_float3 (pSH13_11) * *_S835 + make_float3 (pSH14_11) * *_S836 + make_float3 (pSH15_11) * *_S837) + (make_float3 (pSH16_7) * *_S839 + make_float3 (pSH17_8) * *_S840 + make_float3 (pSH18_8) * *_S841 + make_float3 (pSH19_8) * *_S842 + make_float3 (_S838 - 1.00623059272766113f * pSH6_21) * *_S843 + make_float3 (pSH21_8) * *_S844 + make_float3 (pSH22_8) * *_S845 + make_float3 (pSH23_8) * *_S846 + make_float3 (pSH24_8) * *_S847);
    float3  _S848 = v_colors_20 * make_float3 (float((colors_29.x) >= -0.5f), float((colors_29.y) >= -0.5f), float((colors_29.z) >= -0.5f));
    float3  v_viewdir_17 = {};
    *v_coeff_dc_20 = *v_coeff_dc_20 + make_float3 (0.282094806432724f) * _S848;
    float3  * _S849 = v_coeffs_20 + int(0);
    *_S849 = *_S849 + make_float3 (-0.48860251903533936f * y_30) * _S848;
    float3  * _S850 = v_coeffs_20 + int(1);
    *_S850 = *_S850 + make_float3 (0.48860251903533936f * z_28) * _S848;
    float3  * _S851 = v_coeffs_20 + int(2);
    *_S851 = *_S851 + make_float3 (-0.48860251903533936f * x_31) * _S848;
    float _S852 = -0.48860251903533936f * dot_0(*_S824, _S848);
    float _S853 = -0.48860251903533936f * dot_0(*_S822, _S848);
    float _S854 = 0.48860251903533936f * dot_0(*_S823, _S848);
    float3  * _S855 = v_coeffs_20 + int(3);
    *_S855 = *_S855 + make_float3 (pSH4_11) * _S848;
    float3  * _S856 = v_coeffs_20 + int(4);
    *_S856 = *_S856 + make_float3 (pSH5_14) * _S848;
    float3  * _S857 = v_coeffs_20 + int(5);
    *_S857 = *_S857 + make_float3 (pSH6_21) * _S848;
    float3  * _S858 = v_coeffs_20 + int(6);
    *_S858 = *_S858 + make_float3 (pSH7_14) * _S848;
    float3  * _S859 = v_coeffs_20 + int(7);
    *_S859 = *_S859 + make_float3 (pSH8_14) * _S848;
    float fC1_y_11 = -2.0f * y_30;
    float fS1_x_11 = 2.0f * y_30;
    float pSH6_z_8 = 1.89234936237335205f * z_28;
    float pSH8_x_14 = 0.54627424478530884f * _S825;
    float v_x_19 = _S852 + dot_0(_S848, make_float3 (0.54627424478530884f * fS1_x_11) * *_S826 + make_float3 (pSH8_x_14) * *_S830 + make_float3 (fTmp0B_25) * *_S829);
    float v_y_19 = _S853 + dot_0(_S848, make_float3 (pSH8_x_14) * *_S826 + make_float3 (0.54627424478530884f * fC1_y_11) * *_S830 + make_float3 (fTmp0B_25) * *_S827);
    float v_z_19 = _S854 + dot_0(_S848, make_float3 (pSH6_z_8) * *_S828 + make_float3 (-1.09254848957061768f * x_31) * *_S829 + make_float3 (-1.09254848957061768f * y_30) * *_S827);
    float3  * _S860 = v_coeffs_20 + int(8);
    *_S860 = *_S860 + make_float3 (pSH9_9) * _S848;
    float3  * _S861 = v_coeffs_20 + int(9);
    *_S861 = *_S861 + make_float3 (pSH10_11) * _S848;
    float3  * _S862 = v_coeffs_20 + int(10);
    *_S862 = *_S862 + make_float3 (pSH11_11) * _S848;
    float3  * _S863 = v_coeffs_20 + int(11);
    *_S863 = *_S863 + make_float3 (pSH12_18) * _S848;
    float3  * _S864 = v_coeffs_20 + int(12);
    *_S864 = *_S864 + make_float3 (pSH13_11) * _S848;
    float3  * _S865 = v_coeffs_20 + int(13);
    *_S865 = *_S865 + make_float3 (pSH14_11) * _S848;
    float3  * _S866 = v_coeffs_20 + int(14);
    *_S866 = *_S866 + make_float3 (pSH15_11) * _S848;
    float fTmp0C_z_11 = -4.57045793533325195f * z_28;
    float _S867 = x_31 * _S825;
    float fC2_x_8 = fC1_20 + _S867 - y_30 * fS1_x_11;
    float _S868 = y_30 * _S825;
    float fC2_y_8 = x_31 * fC1_y_11 - fS1_20 - _S868;
    float fS2_x_8 = fS1_20 + x_31 * fS1_x_11 + _S868;
    float fS2_y_8 = _S867 + fC1_20 + y_30 * fC1_y_11;
    float pSH12_z_8 = 5.59764480590820312f * z2_20 - 1.11952900886535645f;
    float pSH14_x_11 = fTmp1B_20 * _S825;
    float v_x_20 = v_x_19 + dot_0(_S848, make_float3 (-0.59004360437393188f * fS2_x_8) * *_S831 + make_float3 (-0.59004360437393188f * fC2_x_8) * *_S837 + make_float3 (fTmp1B_20 * fS1_x_11) * *_S832 + make_float3 (pSH14_x_11) * *_S836 + make_float3 (fTmp0C_20) * *_S835);
    float v_y_20 = v_y_19 + dot_0(_S848, make_float3 (-0.59004360437393188f * fS2_y_8) * *_S831 + make_float3 (-0.59004360437393188f * fC2_y_8) * *_S837 + make_float3 (pSH14_x_11) * *_S832 + make_float3 (fTmp1B_20 * fC1_y_11) * *_S836 + make_float3 (fTmp0C_20) * *_S833);
    float v_z_20 = v_z_19 + dot_0(_S848, make_float3 (pSH12_z_8) * *_S834 + make_float3 (fTmp0C_z_11 * x_31) * *_S835 + make_float3 (fTmp0C_z_11 * y_30) * *_S833 + make_float3 (1.44530570507049561f * fC1_20) * *_S836 + make_float3 (1.44530570507049561f * fS1_20) * *_S832);
    float pSH20_8 = _S838 + -1.00623059272766113f * pSH6_21;
    float3  * _S869 = v_coeffs_20 + int(15);
    *_S869 = *_S869 + make_float3 (pSH16_7) * _S848;
    float3  * _S870 = v_coeffs_20 + int(16);
    *_S870 = *_S870 + make_float3 (pSH17_8) * _S848;
    float3  * _S871 = v_coeffs_20 + int(17);
    *_S871 = *_S871 + make_float3 (pSH18_8) * _S848;
    float3  * _S872 = v_coeffs_20 + int(18);
    *_S872 = *_S872 + make_float3 (pSH19_8) * _S848;
    float3  * _S873 = v_coeffs_20 + int(19);
    *_S873 = *_S873 + make_float3 (pSH20_8) * _S848;
    float3  * _S874 = v_coeffs_20 + int(20);
    *_S874 = *_S874 + make_float3 (pSH21_8) * _S848;
    float3  * _S875 = v_coeffs_20 + int(21);
    *_S875 = *_S875 + make_float3 (pSH22_8) * _S848;
    float3  * _S876 = v_coeffs_20 + int(22);
    *_S876 = *_S876 + make_float3 (pSH23_8) * _S848;
    float3  * _S877 = v_coeffs_20 + int(23);
    *_S877 = *_S877 + make_float3 (pSH24_8) * _S848;
    float fTmp0D_z_8 = -14.04997730255126953f * z2_20 + 2.00713968276977539f;
    float fTmp1C_z_8 = 6.62322282791137695f * z_28;
    float pSH22_x_8 = fTmp1C_15 * _S825;
    float3  dir_n_35 = make_float3 (x_31, y_30, z_28);
    float3  v_dir_n_35 = make_float3 (v_x_20 + dot_0(_S848, make_float3 (0.62583571672439575f * (fS2_15 + y_30 * fC2_x_8 + x_31 * fS2_x_8)) * *_S839 + make_float3 (0.62583571672439575f * (fC2_15 + x_31 * fC2_x_8 - y_30 * fS2_x_8)) * *_S847 + make_float3 (fTmp2B_15 * fS2_x_8) * *_S840 + make_float3 (fTmp2B_15 * fC2_x_8) * *_S846 + make_float3 (fTmp1C_15 * fS1_x_11) * *_S841 + make_float3 (pSH22_x_8) * *_S845 + make_float3 (fTmp0D_15) * *_S844), v_y_20 + dot_0(_S848, make_float3 (0.62583571672439575f * (x_31 * fS2_y_8 + fC2_15 + y_30 * fC2_y_8)) * *_S839 + make_float3 (0.62583571672439575f * (x_31 * fC2_y_8 - fS2_15 - y_30 * fS2_y_8)) * *_S847 + make_float3 (fTmp2B_15 * fS2_y_8) * *_S840 + make_float3 (fTmp2B_15 * fC2_y_8) * *_S846 + make_float3 (pSH22_x_8) * *_S841 + make_float3 (fTmp1C_15 * fC1_y_11) * *_S845 + make_float3 (fTmp0D_15) * *_S842), v_z_20 + dot_0(_S848, make_float3 (1.9843134880065918f * (pSH12_18 + z_28 * pSH12_z_8) + -1.00623059272766113f * pSH6_z_8) * *_S843 + make_float3 (fTmp0D_z_8 * x_31) * *_S844 + make_float3 (fTmp0D_z_8 * y_30) * *_S842 + make_float3 (fTmp1C_z_8 * fC1_20) * *_S845 + make_float3 (fTmp1C_z_8 * fS1_20) * *_S841 + make_float3 (-1.77013075351715088f * fC2_15) * *_S846 + make_float3 (-1.77013075351715088f * fS2_15) * *_S840));
    float3  v_viewdir_18 = v_viewdir_17 + (v_dir_n_35 - make_float3 (dot_0(v_dir_n_35, dir_n_35)) * dir_n_35) * make_float3 (inv_norm_20);
    Matrix<float, 3, 3>  _S878 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S879;
    (&_S879)->primal_0 = _S817;
    (&_S879)->differential_0 = _S878;
    float3  _S880 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S881;
    (&_S881)->primal_0 = t_18;
    (&_S881)->differential_0 = _S880;
    s_bwd_prop_mul_0(&_S879, &_S881, v_viewdir_18);
    Matrix<float, 3, 3>  _S882 = transpose_0(_S879.differential_0);
    *v_mean_12 = *v_mean_12 + v_viewdir_18;
    *v_R_12 = *v_R_12 + _S882;
    *v_t_12 = *v_t_12 + _S881.differential_0;
    return;
}

inline __device__ void sh0_to_color_vjp_atomic(float3  mean_19, Matrix<float, 3, 3>  R_19, float3  t_19, float3  coeff_dc_33, float3  * coeffs_33, float3  v_colors_21, float3  * v_coeff_dc_21, float3  * v_coeffs_21, float3  * v_mean_13, Matrix<float, 3, 3>  * v_R_13, float3  * v_t_13)
{
    float3  colors_30 = make_float3 (0.282094806432724f) * coeff_dc_33;
    *v_coeff_dc_21 = *v_coeff_dc_21 + make_float3 (0.282094806432724f) * (v_colors_21 * make_float3 (float((colors_30.x) >= -0.5f), float((colors_30.y) >= -0.5f), float((colors_30.z) >= -0.5f)));
    float3  v_viewdir_19 = {};
    Matrix<float, 3, 3>  _S883 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S884;
    (&_S884)->primal_0 = transpose_0(R_19);
    (&_S884)->differential_0 = _S883;
    float3  _S885 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S886;
    (&_S886)->primal_0 = t_19;
    (&_S886)->differential_0 = _S885;
    s_bwd_prop_mul_0(&_S884, &_S886, v_viewdir_19);
    Matrix<float, 3, 3>  _S887 = transpose_0(_S884.differential_0);
    *v_mean_13 = *v_mean_13 + v_viewdir_19;
    *v_R_13 = *v_R_13 + _S887;
    *v_t_13 = *v_t_13 + _S886.differential_0;
    return;
}

inline __device__ void sh1_to_color_vjp_atomic(float3  mean_20, Matrix<float, 3, 3>  R_20, float3  t_20, float3  coeff_dc_34, float3  * coeffs_34, float3  v_colors_22, float3  * v_coeff_dc_22, float3  * v_coeffs_22, float3  * v_mean_14, Matrix<float, 3, 3>  * v_R_14, float3  * v_t_14)
{
    Matrix<float, 3, 3>  _S888 = transpose_0(R_20);
    float3  _S889 = mean_20 + mul_0(_S888, t_20);
    float _S890 = _S889.x;
    float _S891 = _S889.y;
    float _S892 = _S889.z;
    float inv_norm_21 = (F32_rsqrt((_S890 * _S890 + _S891 * _S891 + _S892 * _S892)));
    float x_32 = _S890 * inv_norm_21;
    float y_31 = _S891 * inv_norm_21;
    float z_29 = _S892 * inv_norm_21;
    float3  * _S893 = coeffs_34 + int(0);
    float3  * _S894 = coeffs_34 + int(1);
    float3  * _S895 = coeffs_34 + int(2);
    float3  colors_31 = make_float3 (0.282094806432724f) * coeff_dc_34 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * *_S893 + make_float3 (z_29) * *_S894 - make_float3 (x_32) * *_S895);
    float3  _S896 = v_colors_22 * make_float3 (float((colors_31.x) >= -0.5f), float((colors_31.y) >= -0.5f), float((colors_31.z) >= -0.5f));
    float3  v_viewdir_20 = {};
    *v_coeff_dc_22 = *v_coeff_dc_22 + make_float3 (0.282094806432724f) * _S896;
    float3  temp_48 = make_float3 (-0.48860251903533936f * y_31) * _S896;
    float _S897 = dot_0(temp_48, temp_48);
    bool _S898;
    if((F32_isfinite((_S897))))
    {
        _S898 = _S897 != 0.0f;
    }
    else
    {
        _S898 = false;
    }
    if(_S898)
    {
        float3  * _S899 = v_coeffs_22 + int(0);
        float _S900 = atomicAdd(&(_S899->x), temp_48.x);
        float _S901 = atomicAdd(&(_S899->y), temp_48.y);
        float _S902 = atomicAdd(&(_S899->z), temp_48.z);
    }
    float3  temp_49 = make_float3 (0.48860251903533936f * z_29) * _S896;
    float _S903 = dot_0(temp_49, temp_49);
    if((F32_isfinite((_S903))))
    {
        _S898 = _S903 != 0.0f;
    }
    else
    {
        _S898 = false;
    }
    if(_S898)
    {
        float3  * _S904 = v_coeffs_22 + int(1);
        float _S905 = atomicAdd(&(_S904->x), temp_49.x);
        float _S906 = atomicAdd(&(_S904->y), temp_49.y);
        float _S907 = atomicAdd(&(_S904->z), temp_49.z);
    }
    float3  temp_50 = make_float3 (-0.48860251903533936f * x_32) * _S896;
    float _S908 = dot_0(temp_50, temp_50);
    if((F32_isfinite((_S908))))
    {
        _S898 = _S908 != 0.0f;
    }
    else
    {
        _S898 = false;
    }
    if(_S898)
    {
        float3  * _S909 = v_coeffs_22 + int(2);
        float _S910 = atomicAdd(&(_S909->x), temp_50.x);
        float _S911 = atomicAdd(&(_S909->y), temp_50.y);
        float _S912 = atomicAdd(&(_S909->z), temp_50.z);
    }
    float3  dir_n_36 = make_float3 (x_32, y_31, z_29);
    float3  v_dir_n_36 = make_float3 (-0.48860251903533936f * dot_0(*_S895, _S896), -0.48860251903533936f * dot_0(*_S893, _S896), 0.48860251903533936f * dot_0(*_S894, _S896));
    float3  v_viewdir_21 = v_viewdir_20 + (v_dir_n_36 - make_float3 (dot_0(v_dir_n_36, dir_n_36)) * dir_n_36) * make_float3 (inv_norm_21);
    Matrix<float, 3, 3>  _S913 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S914;
    (&_S914)->primal_0 = _S888;
    (&_S914)->differential_0 = _S913;
    float3  _S915 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S916;
    (&_S916)->primal_0 = t_20;
    (&_S916)->differential_0 = _S915;
    s_bwd_prop_mul_0(&_S914, &_S916, v_viewdir_21);
    Matrix<float, 3, 3>  _S917 = transpose_0(_S914.differential_0);
    *v_mean_14 = *v_mean_14 + v_viewdir_21;
    *v_R_14 = *v_R_14 + _S917;
    *v_t_14 = *v_t_14 + _S916.differential_0;
    return;
}

inline __device__ void sh2_to_color_vjp_atomic(float3  mean_21, Matrix<float, 3, 3>  R_21, float3  t_21, float3  coeff_dc_35, float3  * coeffs_35, float3  v_colors_23, float3  * v_coeff_dc_23, float3  * v_coeffs_23, float3  * v_mean_15, Matrix<float, 3, 3>  * v_R_15, float3  * v_t_15)
{
    Matrix<float, 3, 3>  _S918 = transpose_0(R_21);
    float3  _S919 = mean_21 + mul_0(_S918, t_21);
    float _S920 = _S919.x;
    float _S921 = _S919.y;
    float _S922 = _S919.z;
    float inv_norm_22 = (F32_rsqrt((_S920 * _S920 + _S921 * _S921 + _S922 * _S922)));
    float x_33 = _S920 * inv_norm_22;
    float y_32 = _S921 * inv_norm_22;
    float z_30 = _S922 * inv_norm_22;
    float3  * _S923 = coeffs_35 + int(0);
    float3  * _S924 = coeffs_35 + int(1);
    float3  * _S925 = coeffs_35 + int(2);
    float fTmp0B_26 = -1.09254848957061768f * z_30;
    float _S926 = 2.0f * x_33;
    float pSH6_22 = 0.94617468118667603f * (z_30 * z_30) - 0.31539157032966614f;
    float pSH7_15 = fTmp0B_26 * x_33;
    float pSH5_15 = fTmp0B_26 * y_32;
    float pSH8_15 = 0.54627424478530884f * (x_33 * x_33 - y_32 * y_32);
    float pSH4_12 = 0.54627424478530884f * (_S926 * y_32);
    float3  * _S927 = coeffs_35 + int(3);
    float3  * _S928 = coeffs_35 + int(4);
    float3  * _S929 = coeffs_35 + int(5);
    float3  * _S930 = coeffs_35 + int(6);
    float3  * _S931 = coeffs_35 + int(7);
    float3  colors_32 = make_float3 (0.282094806432724f) * coeff_dc_35 + make_float3 (0.48860251903533936f) * (make_float3 (- y_32) * *_S923 + make_float3 (z_30) * *_S924 - make_float3 (x_33) * *_S925) + (make_float3 (pSH4_12) * *_S927 + make_float3 (pSH5_15) * *_S928 + make_float3 (pSH6_22) * *_S929 + make_float3 (pSH7_15) * *_S930 + make_float3 (pSH8_15) * *_S931);
    float3  _S932 = v_colors_23 * make_float3 (float((colors_32.x) >= -0.5f), float((colors_32.y) >= -0.5f), float((colors_32.z) >= -0.5f));
    *v_coeff_dc_23 = *v_coeff_dc_23 + make_float3 (0.282094806432724f) * _S932;
    float3  v_viewdir_22 = {};
    float3  temp_51 = make_float3 (-0.48860251903533936f * y_32) * _S932;
    float _S933 = dot_0(temp_51, temp_51);
    bool _S934;
    if((F32_isfinite((_S933))))
    {
        _S934 = _S933 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S935 = v_coeffs_23 + int(0);
        float _S936 = atomicAdd(&(_S935->x), temp_51.x);
        float _S937 = atomicAdd(&(_S935->y), temp_51.y);
        float _S938 = atomicAdd(&(_S935->z), temp_51.z);
    }
    float3  temp_52 = make_float3 (0.48860251903533936f * z_30) * _S932;
    float _S939 = dot_0(temp_52, temp_52);
    if((F32_isfinite((_S939))))
    {
        _S934 = _S939 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S940 = v_coeffs_23 + int(1);
        float _S941 = atomicAdd(&(_S940->x), temp_52.x);
        float _S942 = atomicAdd(&(_S940->y), temp_52.y);
        float _S943 = atomicAdd(&(_S940->z), temp_52.z);
    }
    float3  temp_53 = make_float3 (-0.48860251903533936f * x_33) * _S932;
    float _S944 = dot_0(temp_53, temp_53);
    if((F32_isfinite((_S944))))
    {
        _S934 = _S944 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S945 = v_coeffs_23 + int(2);
        float _S946 = atomicAdd(&(_S945->x), temp_53.x);
        float _S947 = atomicAdd(&(_S945->y), temp_53.y);
        float _S948 = atomicAdd(&(_S945->z), temp_53.z);
    }
    float _S949 = -0.48860251903533936f * dot_0(*_S925, _S932);
    float _S950 = -0.48860251903533936f * dot_0(*_S923, _S932);
    float _S951 = 0.48860251903533936f * dot_0(*_S924, _S932);
    float3  temp_54 = make_float3 (pSH4_12) * _S932;
    float _S952 = dot_0(temp_54, temp_54);
    if((F32_isfinite((_S952))))
    {
        _S934 = _S952 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S953 = v_coeffs_23 + int(3);
        float _S954 = atomicAdd(&(_S953->x), temp_54.x);
        float _S955 = atomicAdd(&(_S953->y), temp_54.y);
        float _S956 = atomicAdd(&(_S953->z), temp_54.z);
    }
    float3  temp_55 = make_float3 (pSH5_15) * _S932;
    float _S957 = dot_0(temp_55, temp_55);
    if((F32_isfinite((_S957))))
    {
        _S934 = _S957 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S958 = v_coeffs_23 + int(4);
        float _S959 = atomicAdd(&(_S958->x), temp_55.x);
        float _S960 = atomicAdd(&(_S958->y), temp_55.y);
        float _S961 = atomicAdd(&(_S958->z), temp_55.z);
    }
    float3  temp_56 = make_float3 (pSH6_22) * _S932;
    float _S962 = dot_0(temp_56, temp_56);
    if((F32_isfinite((_S962))))
    {
        _S934 = _S962 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S963 = v_coeffs_23 + int(5);
        float _S964 = atomicAdd(&(_S963->x), temp_56.x);
        float _S965 = atomicAdd(&(_S963->y), temp_56.y);
        float _S966 = atomicAdd(&(_S963->z), temp_56.z);
    }
    float3  temp_57 = make_float3 (pSH7_15) * _S932;
    float _S967 = dot_0(temp_57, temp_57);
    if((F32_isfinite((_S967))))
    {
        _S934 = _S967 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S968 = v_coeffs_23 + int(6);
        float _S969 = atomicAdd(&(_S968->x), temp_57.x);
        float _S970 = atomicAdd(&(_S968->y), temp_57.y);
        float _S971 = atomicAdd(&(_S968->z), temp_57.z);
    }
    float3  temp_58 = make_float3 (pSH8_15) * _S932;
    float _S972 = dot_0(temp_58, temp_58);
    if((F32_isfinite((_S972))))
    {
        _S934 = _S972 != 0.0f;
    }
    else
    {
        _S934 = false;
    }
    if(_S934)
    {
        float3  * _S973 = v_coeffs_23 + int(7);
        float _S974 = atomicAdd(&(_S973->x), temp_58.x);
        float _S975 = atomicAdd(&(_S973->y), temp_58.y);
        float _S976 = atomicAdd(&(_S973->z), temp_58.z);
    }
    float pSH8_x_15 = 0.54627424478530884f * _S926;
    float3  dir_n_37 = make_float3 (x_33, y_32, z_30);
    float3  v_dir_n_37 = make_float3 (_S949 + dot_0(_S932, make_float3 (0.54627424478530884f * (2.0f * y_32)) * *_S927 + make_float3 (pSH8_x_15) * *_S931 + make_float3 (fTmp0B_26) * *_S930), _S950 + dot_0(_S932, make_float3 (pSH8_x_15) * *_S927 + make_float3 (0.54627424478530884f * (-2.0f * y_32)) * *_S931 + make_float3 (fTmp0B_26) * *_S928), _S951 + dot_0(_S932, make_float3 (1.89234936237335205f * z_30) * *_S929 + make_float3 (-1.09254848957061768f * x_33) * *_S930 + make_float3 (-1.09254848957061768f * y_32) * *_S928));
    float3  v_viewdir_23 = v_viewdir_22 + (v_dir_n_37 - make_float3 (dot_0(v_dir_n_37, dir_n_37)) * dir_n_37) * make_float3 (inv_norm_22);
    Matrix<float, 3, 3>  _S977 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S978;
    (&_S978)->primal_0 = _S918;
    (&_S978)->differential_0 = _S977;
    float3  _S979 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S980;
    (&_S980)->primal_0 = t_21;
    (&_S980)->differential_0 = _S979;
    s_bwd_prop_mul_0(&_S978, &_S980, v_viewdir_23);
    Matrix<float, 3, 3>  _S981 = transpose_0(_S978.differential_0);
    *v_mean_15 = *v_mean_15 + v_viewdir_23;
    *v_R_15 = *v_R_15 + _S981;
    *v_t_15 = *v_t_15 + _S980.differential_0;
    return;
}

inline __device__ void sh3_to_color_vjp_atomic(float3  mean_22, Matrix<float, 3, 3>  R_22, float3  t_22, float3  coeff_dc_36, float3  * coeffs_36, float3  v_colors_24, float3  * v_coeff_dc_24, float3  * v_coeffs_24, float3  * v_mean_16, Matrix<float, 3, 3>  * v_R_16, float3  * v_t_16)
{
    Matrix<float, 3, 3>  _S982 = transpose_0(R_22);
    float3  _S983 = mean_22 + mul_0(_S982, t_22);
    float _S984 = _S983.x;
    float _S985 = _S983.y;
    float _S986 = _S983.z;
    float inv_norm_23 = (F32_rsqrt((_S984 * _S984 + _S985 * _S985 + _S986 * _S986)));
    float x_34 = _S984 * inv_norm_23;
    float y_33 = _S985 * inv_norm_23;
    float z_31 = _S986 * inv_norm_23;
    float3  * _S987 = coeffs_36 + int(0);
    float3  * _S988 = coeffs_36 + int(1);
    float3  * _S989 = coeffs_36 + int(2);
    float z2_21 = z_31 * z_31;
    float fTmp0B_27 = -1.09254848957061768f * z_31;
    float fC1_21 = x_34 * x_34 - y_33 * y_33;
    float _S990 = 2.0f * x_34;
    float fS1_21 = _S990 * y_33;
    float pSH6_23 = 0.94617468118667603f * z2_21 - 0.31539157032966614f;
    float pSH7_16 = fTmp0B_27 * x_34;
    float pSH5_16 = fTmp0B_27 * y_33;
    float pSH8_16 = 0.54627424478530884f * fC1_21;
    float pSH4_13 = 0.54627424478530884f * fS1_21;
    float3  * _S991 = coeffs_36 + int(3);
    float3  * _S992 = coeffs_36 + int(4);
    float3  * _S993 = coeffs_36 + int(5);
    float3  * _S994 = coeffs_36 + int(6);
    float3  * _S995 = coeffs_36 + int(7);
    float fTmp0C_21 = -2.28522896766662598f * z2_21 + 0.4570457935333252f;
    float fTmp1B_21 = 1.44530570507049561f * z_31;
    float pSH12_19 = z_31 * (1.86588168144226074f * z2_21 - 1.11952900886535645f);
    float pSH13_12 = fTmp0C_21 * x_34;
    float pSH11_12 = fTmp0C_21 * y_33;
    float pSH14_12 = fTmp1B_21 * fC1_21;
    float pSH10_12 = fTmp1B_21 * fS1_21;
    float pSH15_12 = -0.59004360437393188f * (x_34 * fC1_21 - y_33 * fS1_21);
    float pSH9_10 = -0.59004360437393188f * (x_34 * fS1_21 + y_33 * fC1_21);
    float3  * _S996 = coeffs_36 + int(8);
    float3  * _S997 = coeffs_36 + int(9);
    float3  * _S998 = coeffs_36 + int(10);
    float3  * _S999 = coeffs_36 + int(11);
    float3  * _S1000 = coeffs_36 + int(12);
    float3  * _S1001 = coeffs_36 + int(13);
    float3  * _S1002 = coeffs_36 + int(14);
    float3  colors_33 = make_float3 (0.282094806432724f) * coeff_dc_36 + make_float3 (0.48860251903533936f) * (make_float3 (- y_33) * *_S987 + make_float3 (z_31) * *_S988 - make_float3 (x_34) * *_S989) + (make_float3 (pSH4_13) * *_S991 + make_float3 (pSH5_16) * *_S992 + make_float3 (pSH6_23) * *_S993 + make_float3 (pSH7_16) * *_S994 + make_float3 (pSH8_16) * *_S995) + (make_float3 (pSH9_10) * *_S996 + make_float3 (pSH10_12) * *_S997 + make_float3 (pSH11_12) * *_S998 + make_float3 (pSH12_19) * *_S999 + make_float3 (pSH13_12) * *_S1000 + make_float3 (pSH14_12) * *_S1001 + make_float3 (pSH15_12) * *_S1002);
    float3  _S1003 = v_colors_24 * make_float3 (float((colors_33.x) >= -0.5f), float((colors_33.y) >= -0.5f), float((colors_33.z) >= -0.5f));
    float3  v_viewdir_24 = {};
    *v_coeff_dc_24 = *v_coeff_dc_24 + make_float3 (0.282094806432724f) * _S1003;
    float3  temp_59 = make_float3 (-0.48860251903533936f * y_33) * _S1003;
    float _S1004 = dot_0(temp_59, temp_59);
    bool _S1005;
    if((F32_isfinite((_S1004))))
    {
        _S1005 = _S1004 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1006 = v_coeffs_24 + int(0);
        float _S1007 = atomicAdd(&(_S1006->x), temp_59.x);
        float _S1008 = atomicAdd(&(_S1006->y), temp_59.y);
        float _S1009 = atomicAdd(&(_S1006->z), temp_59.z);
    }
    float3  temp_60 = make_float3 (0.48860251903533936f * z_31) * _S1003;
    float _S1010 = dot_0(temp_60, temp_60);
    if((F32_isfinite((_S1010))))
    {
        _S1005 = _S1010 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1011 = v_coeffs_24 + int(1);
        float _S1012 = atomicAdd(&(_S1011->x), temp_60.x);
        float _S1013 = atomicAdd(&(_S1011->y), temp_60.y);
        float _S1014 = atomicAdd(&(_S1011->z), temp_60.z);
    }
    float3  temp_61 = make_float3 (-0.48860251903533936f * x_34) * _S1003;
    float _S1015 = dot_0(temp_61, temp_61);
    if((F32_isfinite((_S1015))))
    {
        _S1005 = _S1015 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1016 = v_coeffs_24 + int(2);
        float _S1017 = atomicAdd(&(_S1016->x), temp_61.x);
        float _S1018 = atomicAdd(&(_S1016->y), temp_61.y);
        float _S1019 = atomicAdd(&(_S1016->z), temp_61.z);
    }
    float _S1020 = -0.48860251903533936f * dot_0(*_S989, _S1003);
    float _S1021 = -0.48860251903533936f * dot_0(*_S987, _S1003);
    float _S1022 = 0.48860251903533936f * dot_0(*_S988, _S1003);
    float3  temp_62 = make_float3 (pSH4_13) * _S1003;
    float _S1023 = dot_0(temp_62, temp_62);
    if((F32_isfinite((_S1023))))
    {
        _S1005 = _S1023 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1024 = v_coeffs_24 + int(3);
        float _S1025 = atomicAdd(&(_S1024->x), temp_62.x);
        float _S1026 = atomicAdd(&(_S1024->y), temp_62.y);
        float _S1027 = atomicAdd(&(_S1024->z), temp_62.z);
    }
    float3  temp_63 = make_float3 (pSH5_16) * _S1003;
    float _S1028 = dot_0(temp_63, temp_63);
    if((F32_isfinite((_S1028))))
    {
        _S1005 = _S1028 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1029 = v_coeffs_24 + int(4);
        float _S1030 = atomicAdd(&(_S1029->x), temp_63.x);
        float _S1031 = atomicAdd(&(_S1029->y), temp_63.y);
        float _S1032 = atomicAdd(&(_S1029->z), temp_63.z);
    }
    float3  temp_64 = make_float3 (pSH6_23) * _S1003;
    float _S1033 = dot_0(temp_64, temp_64);
    if((F32_isfinite((_S1033))))
    {
        _S1005 = _S1033 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1034 = v_coeffs_24 + int(5);
        float _S1035 = atomicAdd(&(_S1034->x), temp_64.x);
        float _S1036 = atomicAdd(&(_S1034->y), temp_64.y);
        float _S1037 = atomicAdd(&(_S1034->z), temp_64.z);
    }
    float3  temp_65 = make_float3 (pSH7_16) * _S1003;
    float _S1038 = dot_0(temp_65, temp_65);
    if((F32_isfinite((_S1038))))
    {
        _S1005 = _S1038 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1039 = v_coeffs_24 + int(6);
        float _S1040 = atomicAdd(&(_S1039->x), temp_65.x);
        float _S1041 = atomicAdd(&(_S1039->y), temp_65.y);
        float _S1042 = atomicAdd(&(_S1039->z), temp_65.z);
    }
    float3  temp_66 = make_float3 (pSH8_16) * _S1003;
    float _S1043 = dot_0(temp_66, temp_66);
    if((F32_isfinite((_S1043))))
    {
        _S1005 = _S1043 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1044 = v_coeffs_24 + int(7);
        float _S1045 = atomicAdd(&(_S1044->x), temp_66.x);
        float _S1046 = atomicAdd(&(_S1044->y), temp_66.y);
        float _S1047 = atomicAdd(&(_S1044->z), temp_66.z);
    }
    float fC1_y_12 = -2.0f * y_33;
    float fS1_x_12 = 2.0f * y_33;
    float pSH8_x_16 = 0.54627424478530884f * _S990;
    float v_x_21 = _S1020 + dot_0(_S1003, make_float3 (0.54627424478530884f * fS1_x_12) * *_S991 + make_float3 (pSH8_x_16) * *_S995 + make_float3 (fTmp0B_27) * *_S994);
    float v_y_21 = _S1021 + dot_0(_S1003, make_float3 (pSH8_x_16) * *_S991 + make_float3 (0.54627424478530884f * fC1_y_12) * *_S995 + make_float3 (fTmp0B_27) * *_S992);
    float v_z_21 = _S1022 + dot_0(_S1003, make_float3 (1.89234936237335205f * z_31) * *_S993 + make_float3 (-1.09254848957061768f * x_34) * *_S994 + make_float3 (-1.09254848957061768f * y_33) * *_S992);
    float3  temp_67 = make_float3 (pSH9_10) * _S1003;
    float _S1048 = dot_0(temp_67, temp_67);
    if((F32_isfinite((_S1048))))
    {
        _S1005 = _S1048 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1049 = v_coeffs_24 + int(8);
        float _S1050 = atomicAdd(&(_S1049->x), temp_67.x);
        float _S1051 = atomicAdd(&(_S1049->y), temp_67.y);
        float _S1052 = atomicAdd(&(_S1049->z), temp_67.z);
    }
    float3  temp_68 = make_float3 (pSH10_12) * _S1003;
    float _S1053 = dot_0(temp_68, temp_68);
    if((F32_isfinite((_S1053))))
    {
        _S1005 = _S1053 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1054 = v_coeffs_24 + int(9);
        float _S1055 = atomicAdd(&(_S1054->x), temp_68.x);
        float _S1056 = atomicAdd(&(_S1054->y), temp_68.y);
        float _S1057 = atomicAdd(&(_S1054->z), temp_68.z);
    }
    float3  temp_69 = make_float3 (pSH11_12) * _S1003;
    float _S1058 = dot_0(temp_69, temp_69);
    if((F32_isfinite((_S1058))))
    {
        _S1005 = _S1058 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1059 = v_coeffs_24 + int(10);
        float _S1060 = atomicAdd(&(_S1059->x), temp_69.x);
        float _S1061 = atomicAdd(&(_S1059->y), temp_69.y);
        float _S1062 = atomicAdd(&(_S1059->z), temp_69.z);
    }
    float3  temp_70 = make_float3 (pSH12_19) * _S1003;
    float _S1063 = dot_0(temp_70, temp_70);
    if((F32_isfinite((_S1063))))
    {
        _S1005 = _S1063 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1064 = v_coeffs_24 + int(11);
        float _S1065 = atomicAdd(&(_S1064->x), temp_70.x);
        float _S1066 = atomicAdd(&(_S1064->y), temp_70.y);
        float _S1067 = atomicAdd(&(_S1064->z), temp_70.z);
    }
    float3  temp_71 = make_float3 (pSH13_12) * _S1003;
    float _S1068 = dot_0(temp_71, temp_71);
    if((F32_isfinite((_S1068))))
    {
        _S1005 = _S1068 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1069 = v_coeffs_24 + int(12);
        float _S1070 = atomicAdd(&(_S1069->x), temp_71.x);
        float _S1071 = atomicAdd(&(_S1069->y), temp_71.y);
        float _S1072 = atomicAdd(&(_S1069->z), temp_71.z);
    }
    float3  temp_72 = make_float3 (pSH14_12) * _S1003;
    float _S1073 = dot_0(temp_72, temp_72);
    if((F32_isfinite((_S1073))))
    {
        _S1005 = _S1073 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1074 = v_coeffs_24 + int(13);
        float _S1075 = atomicAdd(&(_S1074->x), temp_72.x);
        float _S1076 = atomicAdd(&(_S1074->y), temp_72.y);
        float _S1077 = atomicAdd(&(_S1074->z), temp_72.z);
    }
    float3  temp_73 = make_float3 (pSH15_12) * _S1003;
    float _S1078 = dot_0(temp_73, temp_73);
    if((F32_isfinite((_S1078))))
    {
        _S1005 = _S1078 != 0.0f;
    }
    else
    {
        _S1005 = false;
    }
    if(_S1005)
    {
        float3  * _S1079 = v_coeffs_24 + int(14);
        float _S1080 = atomicAdd(&(_S1079->x), temp_73.x);
        float _S1081 = atomicAdd(&(_S1079->y), temp_73.y);
        float _S1082 = atomicAdd(&(_S1079->z), temp_73.z);
    }
    float fTmp0C_z_12 = -4.57045793533325195f * z_31;
    float _S1083 = x_34 * _S990;
    float _S1084 = y_33 * _S990;
    float pSH14_x_12 = fTmp1B_21 * _S990;
    float3  dir_n_38 = make_float3 (x_34, y_33, z_31);
    float3  v_dir_n_38 = make_float3 (v_x_21 + dot_0(_S1003, make_float3 (-0.59004360437393188f * (fS1_21 + x_34 * fS1_x_12 + _S1084)) * *_S996 + make_float3 (-0.59004360437393188f * (fC1_21 + _S1083 - y_33 * fS1_x_12)) * *_S1002 + make_float3 (fTmp1B_21 * fS1_x_12) * *_S997 + make_float3 (pSH14_x_12) * *_S1001 + make_float3 (fTmp0C_21) * *_S1000), v_y_21 + dot_0(_S1003, make_float3 (-0.59004360437393188f * (_S1083 + fC1_21 + y_33 * fC1_y_12)) * *_S996 + make_float3 (-0.59004360437393188f * (x_34 * fC1_y_12 - fS1_21 - _S1084)) * *_S1002 + make_float3 (pSH14_x_12) * *_S997 + make_float3 (fTmp1B_21 * fC1_y_12) * *_S1001 + make_float3 (fTmp0C_21) * *_S998), v_z_21 + dot_0(_S1003, make_float3 (5.59764480590820312f * z2_21 - 1.11952900886535645f) * *_S999 + make_float3 (fTmp0C_z_12 * x_34) * *_S1000 + make_float3 (fTmp0C_z_12 * y_33) * *_S998 + make_float3 (1.44530570507049561f * fC1_21) * *_S1001 + make_float3 (1.44530570507049561f * fS1_21) * *_S997));
    float3  v_viewdir_25 = v_viewdir_24 + (v_dir_n_38 - make_float3 (dot_0(v_dir_n_38, dir_n_38)) * dir_n_38) * make_float3 (inv_norm_23);
    Matrix<float, 3, 3>  _S1085 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1086;
    (&_S1086)->primal_0 = _S982;
    (&_S1086)->differential_0 = _S1085;
    float3  _S1087 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1088;
    (&_S1088)->primal_0 = t_22;
    (&_S1088)->differential_0 = _S1087;
    s_bwd_prop_mul_0(&_S1086, &_S1088, v_viewdir_25);
    Matrix<float, 3, 3>  _S1089 = transpose_0(_S1086.differential_0);
    *v_mean_16 = *v_mean_16 + v_viewdir_25;
    *v_R_16 = *v_R_16 + _S1089;
    *v_t_16 = *v_t_16 + _S1088.differential_0;
    return;
}

inline __device__ void sh4_to_color_vjp_atomic(float3  mean_23, Matrix<float, 3, 3>  R_23, float3  t_23, float3  coeff_dc_37, float3  * coeffs_37, float3  v_colors_25, float3  * v_coeff_dc_25, float3  * v_coeffs_25, float3  * v_mean_17, Matrix<float, 3, 3>  * v_R_17, float3  * v_t_17)
{
    Matrix<float, 3, 3>  _S1090 = transpose_0(R_23);
    float3  _S1091 = mean_23 + mul_0(_S1090, t_23);
    float _S1092 = _S1091.x;
    float _S1093 = _S1091.y;
    float _S1094 = _S1091.z;
    float inv_norm_24 = (F32_rsqrt((_S1092 * _S1092 + _S1093 * _S1093 + _S1094 * _S1094)));
    float x_35 = _S1092 * inv_norm_24;
    float y_34 = _S1093 * inv_norm_24;
    float z_32 = _S1094 * inv_norm_24;
    float3  * _S1095 = coeffs_37 + int(0);
    float3  * _S1096 = coeffs_37 + int(1);
    float3  * _S1097 = coeffs_37 + int(2);
    float z2_22 = z_32 * z_32;
    float fTmp0B_28 = -1.09254848957061768f * z_32;
    float fC1_22 = x_35 * x_35 - y_34 * y_34;
    float _S1098 = 2.0f * x_35;
    float fS1_22 = _S1098 * y_34;
    float pSH6_24 = 0.94617468118667603f * z2_22 - 0.31539157032966614f;
    float pSH7_17 = fTmp0B_28 * x_35;
    float pSH5_17 = fTmp0B_28 * y_34;
    float pSH8_17 = 0.54627424478530884f * fC1_22;
    float pSH4_14 = 0.54627424478530884f * fS1_22;
    float3  * _S1099 = coeffs_37 + int(3);
    float3  * _S1100 = coeffs_37 + int(4);
    float3  * _S1101 = coeffs_37 + int(5);
    float3  * _S1102 = coeffs_37 + int(6);
    float3  * _S1103 = coeffs_37 + int(7);
    float fTmp0C_22 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_22 = 1.44530570507049561f * z_32;
    float fC2_16 = x_35 * fC1_22 - y_34 * fS1_22;
    float fS2_16 = x_35 * fS1_22 + y_34 * fC1_22;
    float pSH12_20 = z_32 * (1.86588168144226074f * z2_22 - 1.11952900886535645f);
    float pSH13_13 = fTmp0C_22 * x_35;
    float pSH11_13 = fTmp0C_22 * y_34;
    float pSH14_13 = fTmp1B_22 * fC1_22;
    float pSH10_13 = fTmp1B_22 * fS1_22;
    float pSH15_13 = -0.59004360437393188f * fC2_16;
    float pSH9_11 = -0.59004360437393188f * fS2_16;
    float3  * _S1104 = coeffs_37 + int(8);
    float3  * _S1105 = coeffs_37 + int(9);
    float3  * _S1106 = coeffs_37 + int(10);
    float3  * _S1107 = coeffs_37 + int(11);
    float3  * _S1108 = coeffs_37 + int(12);
    float3  * _S1109 = coeffs_37 + int(13);
    float3  * _S1110 = coeffs_37 + int(14);
    float fTmp0D_16 = z_32 * (-4.68332576751708984f * z2_22 + 2.00713968276977539f);
    float fTmp1C_16 = 3.31161141395568848f * z2_22 - 0.47308734059333801f;
    float fTmp2B_16 = -1.77013075351715088f * z_32;
    float _S1111 = 1.9843134880065918f * z_32 * pSH12_20;
    float pSH21_9 = fTmp0D_16 * x_35;
    float pSH19_9 = fTmp0D_16 * y_34;
    float pSH22_9 = fTmp1C_16 * fC1_22;
    float pSH18_9 = fTmp1C_16 * fS1_22;
    float pSH23_9 = fTmp2B_16 * fC2_16;
    float pSH17_9 = fTmp2B_16 * fS2_16;
    float pSH24_9 = 0.62583571672439575f * (x_35 * fC2_16 - y_34 * fS2_16);
    float pSH16_8 = 0.62583571672439575f * (x_35 * fS2_16 + y_34 * fC2_16);
    float3  * _S1112 = coeffs_37 + int(15);
    float3  * _S1113 = coeffs_37 + int(16);
    float3  * _S1114 = coeffs_37 + int(17);
    float3  * _S1115 = coeffs_37 + int(18);
    float3  * _S1116 = coeffs_37 + int(19);
    float3  * _S1117 = coeffs_37 + int(20);
    float3  * _S1118 = coeffs_37 + int(21);
    float3  * _S1119 = coeffs_37 + int(22);
    float3  * _S1120 = coeffs_37 + int(23);
    float3  colors_34 = make_float3 (0.282094806432724f) * coeff_dc_37 + make_float3 (0.48860251903533936f) * (make_float3 (- y_34) * *_S1095 + make_float3 (z_32) * *_S1096 - make_float3 (x_35) * *_S1097) + (make_float3 (pSH4_14) * *_S1099 + make_float3 (pSH5_17) * *_S1100 + make_float3 (pSH6_24) * *_S1101 + make_float3 (pSH7_17) * *_S1102 + make_float3 (pSH8_17) * *_S1103) + (make_float3 (pSH9_11) * *_S1104 + make_float3 (pSH10_13) * *_S1105 + make_float3 (pSH11_13) * *_S1106 + make_float3 (pSH12_20) * *_S1107 + make_float3 (pSH13_13) * *_S1108 + make_float3 (pSH14_13) * *_S1109 + make_float3 (pSH15_13) * *_S1110) + (make_float3 (pSH16_8) * *_S1112 + make_float3 (pSH17_9) * *_S1113 + make_float3 (pSH18_9) * *_S1114 + make_float3 (pSH19_9) * *_S1115 + make_float3 (_S1111 - 1.00623059272766113f * pSH6_24) * *_S1116 + make_float3 (pSH21_9) * *_S1117 + make_float3 (pSH22_9) * *_S1118 + make_float3 (pSH23_9) * *_S1119 + make_float3 (pSH24_9) * *_S1120);
    float3  _S1121 = v_colors_25 * make_float3 (float((colors_34.x) >= -0.5f), float((colors_34.y) >= -0.5f), float((colors_34.z) >= -0.5f));
    float3  v_viewdir_26 = {};
    *v_coeff_dc_25 = *v_coeff_dc_25 + make_float3 (0.282094806432724f) * _S1121;
    float3  temp_74 = make_float3 (-0.48860251903533936f * y_34) * _S1121;
    float _S1122 = dot_0(temp_74, temp_74);
    bool _S1123;
    if((F32_isfinite((_S1122))))
    {
        _S1123 = _S1122 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1124 = v_coeffs_25 + int(0);
        float _S1125 = atomicAdd(&(_S1124->x), temp_74.x);
        float _S1126 = atomicAdd(&(_S1124->y), temp_74.y);
        float _S1127 = atomicAdd(&(_S1124->z), temp_74.z);
    }
    float3  temp_75 = make_float3 (0.48860251903533936f * z_32) * _S1121;
    float _S1128 = dot_0(temp_75, temp_75);
    if((F32_isfinite((_S1128))))
    {
        _S1123 = _S1128 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1129 = v_coeffs_25 + int(1);
        float _S1130 = atomicAdd(&(_S1129->x), temp_75.x);
        float _S1131 = atomicAdd(&(_S1129->y), temp_75.y);
        float _S1132 = atomicAdd(&(_S1129->z), temp_75.z);
    }
    float3  temp_76 = make_float3 (-0.48860251903533936f * x_35) * _S1121;
    float _S1133 = dot_0(temp_76, temp_76);
    if((F32_isfinite((_S1133))))
    {
        _S1123 = _S1133 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1134 = v_coeffs_25 + int(2);
        float _S1135 = atomicAdd(&(_S1134->x), temp_76.x);
        float _S1136 = atomicAdd(&(_S1134->y), temp_76.y);
        float _S1137 = atomicAdd(&(_S1134->z), temp_76.z);
    }
    float _S1138 = -0.48860251903533936f * dot_0(*_S1097, _S1121);
    float _S1139 = -0.48860251903533936f * dot_0(*_S1095, _S1121);
    float _S1140 = 0.48860251903533936f * dot_0(*_S1096, _S1121);
    float3  temp_77 = make_float3 (pSH4_14) * _S1121;
    float _S1141 = dot_0(temp_77, temp_77);
    if((F32_isfinite((_S1141))))
    {
        _S1123 = _S1141 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1142 = v_coeffs_25 + int(3);
        float _S1143 = atomicAdd(&(_S1142->x), temp_77.x);
        float _S1144 = atomicAdd(&(_S1142->y), temp_77.y);
        float _S1145 = atomicAdd(&(_S1142->z), temp_77.z);
    }
    float3  temp_78 = make_float3 (pSH5_17) * _S1121;
    float _S1146 = dot_0(temp_78, temp_78);
    if((F32_isfinite((_S1146))))
    {
        _S1123 = _S1146 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1147 = v_coeffs_25 + int(4);
        float _S1148 = atomicAdd(&(_S1147->x), temp_78.x);
        float _S1149 = atomicAdd(&(_S1147->y), temp_78.y);
        float _S1150 = atomicAdd(&(_S1147->z), temp_78.z);
    }
    float3  temp_79 = make_float3 (pSH6_24) * _S1121;
    float _S1151 = dot_0(temp_79, temp_79);
    if((F32_isfinite((_S1151))))
    {
        _S1123 = _S1151 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1152 = v_coeffs_25 + int(5);
        float _S1153 = atomicAdd(&(_S1152->x), temp_79.x);
        float _S1154 = atomicAdd(&(_S1152->y), temp_79.y);
        float _S1155 = atomicAdd(&(_S1152->z), temp_79.z);
    }
    float3  temp_80 = make_float3 (pSH7_17) * _S1121;
    float _S1156 = dot_0(temp_80, temp_80);
    if((F32_isfinite((_S1156))))
    {
        _S1123 = _S1156 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1157 = v_coeffs_25 + int(6);
        float _S1158 = atomicAdd(&(_S1157->x), temp_80.x);
        float _S1159 = atomicAdd(&(_S1157->y), temp_80.y);
        float _S1160 = atomicAdd(&(_S1157->z), temp_80.z);
    }
    float3  temp_81 = make_float3 (pSH8_17) * _S1121;
    float _S1161 = dot_0(temp_81, temp_81);
    if((F32_isfinite((_S1161))))
    {
        _S1123 = _S1161 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1162 = v_coeffs_25 + int(7);
        float _S1163 = atomicAdd(&(_S1162->x), temp_81.x);
        float _S1164 = atomicAdd(&(_S1162->y), temp_81.y);
        float _S1165 = atomicAdd(&(_S1162->z), temp_81.z);
    }
    float fC1_y_13 = -2.0f * y_34;
    float fS1_x_13 = 2.0f * y_34;
    float pSH6_z_9 = 1.89234936237335205f * z_32;
    float pSH8_x_17 = 0.54627424478530884f * _S1098;
    float v_x_22 = _S1138 + dot_0(_S1121, make_float3 (0.54627424478530884f * fS1_x_13) * *_S1099 + make_float3 (pSH8_x_17) * *_S1103 + make_float3 (fTmp0B_28) * *_S1102);
    float v_y_22 = _S1139 + dot_0(_S1121, make_float3 (pSH8_x_17) * *_S1099 + make_float3 (0.54627424478530884f * fC1_y_13) * *_S1103 + make_float3 (fTmp0B_28) * *_S1100);
    float v_z_22 = _S1140 + dot_0(_S1121, make_float3 (pSH6_z_9) * *_S1101 + make_float3 (-1.09254848957061768f * x_35) * *_S1102 + make_float3 (-1.09254848957061768f * y_34) * *_S1100);
    float3  temp_82 = make_float3 (pSH9_11) * _S1121;
    float _S1166 = dot_0(temp_82, temp_82);
    if((F32_isfinite((_S1166))))
    {
        _S1123 = _S1166 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1167 = v_coeffs_25 + int(8);
        float _S1168 = atomicAdd(&(_S1167->x), temp_82.x);
        float _S1169 = atomicAdd(&(_S1167->y), temp_82.y);
        float _S1170 = atomicAdd(&(_S1167->z), temp_82.z);
    }
    float3  temp_83 = make_float3 (pSH10_13) * _S1121;
    float _S1171 = dot_0(temp_83, temp_83);
    if((F32_isfinite((_S1171))))
    {
        _S1123 = _S1171 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1172 = v_coeffs_25 + int(9);
        float _S1173 = atomicAdd(&(_S1172->x), temp_83.x);
        float _S1174 = atomicAdd(&(_S1172->y), temp_83.y);
        float _S1175 = atomicAdd(&(_S1172->z), temp_83.z);
    }
    float3  temp_84 = make_float3 (pSH11_13) * _S1121;
    float _S1176 = dot_0(temp_84, temp_84);
    if((F32_isfinite((_S1176))))
    {
        _S1123 = _S1176 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1177 = v_coeffs_25 + int(10);
        float _S1178 = atomicAdd(&(_S1177->x), temp_84.x);
        float _S1179 = atomicAdd(&(_S1177->y), temp_84.y);
        float _S1180 = atomicAdd(&(_S1177->z), temp_84.z);
    }
    float3  temp_85 = make_float3 (pSH12_20) * _S1121;
    float _S1181 = dot_0(temp_85, temp_85);
    if((F32_isfinite((_S1181))))
    {
        _S1123 = _S1181 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1182 = v_coeffs_25 + int(11);
        float _S1183 = atomicAdd(&(_S1182->x), temp_85.x);
        float _S1184 = atomicAdd(&(_S1182->y), temp_85.y);
        float _S1185 = atomicAdd(&(_S1182->z), temp_85.z);
    }
    float3  temp_86 = make_float3 (pSH13_13) * _S1121;
    float _S1186 = dot_0(temp_86, temp_86);
    if((F32_isfinite((_S1186))))
    {
        _S1123 = _S1186 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1187 = v_coeffs_25 + int(12);
        float _S1188 = atomicAdd(&(_S1187->x), temp_86.x);
        float _S1189 = atomicAdd(&(_S1187->y), temp_86.y);
        float _S1190 = atomicAdd(&(_S1187->z), temp_86.z);
    }
    float3  temp_87 = make_float3 (pSH14_13) * _S1121;
    float _S1191 = dot_0(temp_87, temp_87);
    if((F32_isfinite((_S1191))))
    {
        _S1123 = _S1191 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1192 = v_coeffs_25 + int(13);
        float _S1193 = atomicAdd(&(_S1192->x), temp_87.x);
        float _S1194 = atomicAdd(&(_S1192->y), temp_87.y);
        float _S1195 = atomicAdd(&(_S1192->z), temp_87.z);
    }
    float3  temp_88 = make_float3 (pSH15_13) * _S1121;
    float _S1196 = dot_0(temp_88, temp_88);
    if((F32_isfinite((_S1196))))
    {
        _S1123 = _S1196 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1197 = v_coeffs_25 + int(14);
        float _S1198 = atomicAdd(&(_S1197->x), temp_88.x);
        float _S1199 = atomicAdd(&(_S1197->y), temp_88.y);
        float _S1200 = atomicAdd(&(_S1197->z), temp_88.z);
    }
    float fTmp0C_z_13 = -4.57045793533325195f * z_32;
    float _S1201 = x_35 * _S1098;
    float fC2_x_9 = fC1_22 + _S1201 - y_34 * fS1_x_13;
    float _S1202 = y_34 * _S1098;
    float fC2_y_9 = x_35 * fC1_y_13 - fS1_22 - _S1202;
    float fS2_x_9 = fS1_22 + x_35 * fS1_x_13 + _S1202;
    float fS2_y_9 = _S1201 + fC1_22 + y_34 * fC1_y_13;
    float pSH12_z_9 = 5.59764480590820312f * z2_22 - 1.11952900886535645f;
    float pSH14_x_13 = fTmp1B_22 * _S1098;
    float v_x_23 = v_x_22 + dot_0(_S1121, make_float3 (-0.59004360437393188f * fS2_x_9) * *_S1104 + make_float3 (-0.59004360437393188f * fC2_x_9) * *_S1110 + make_float3 (fTmp1B_22 * fS1_x_13) * *_S1105 + make_float3 (pSH14_x_13) * *_S1109 + make_float3 (fTmp0C_22) * *_S1108);
    float v_y_23 = v_y_22 + dot_0(_S1121, make_float3 (-0.59004360437393188f * fS2_y_9) * *_S1104 + make_float3 (-0.59004360437393188f * fC2_y_9) * *_S1110 + make_float3 (pSH14_x_13) * *_S1105 + make_float3 (fTmp1B_22 * fC1_y_13) * *_S1109 + make_float3 (fTmp0C_22) * *_S1106);
    float v_z_23 = v_z_22 + dot_0(_S1121, make_float3 (pSH12_z_9) * *_S1107 + make_float3 (fTmp0C_z_13 * x_35) * *_S1108 + make_float3 (fTmp0C_z_13 * y_34) * *_S1106 + make_float3 (1.44530570507049561f * fC1_22) * *_S1109 + make_float3 (1.44530570507049561f * fS1_22) * *_S1105);
    float pSH20_9 = _S1111 + -1.00623059272766113f * pSH6_24;
    float3  temp_89 = make_float3 (pSH16_8) * _S1121;
    float _S1203 = dot_0(temp_89, temp_89);
    if((F32_isfinite((_S1203))))
    {
        _S1123 = _S1203 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1204 = v_coeffs_25 + int(15);
        float _S1205 = atomicAdd(&(_S1204->x), temp_89.x);
        float _S1206 = atomicAdd(&(_S1204->y), temp_89.y);
        float _S1207 = atomicAdd(&(_S1204->z), temp_89.z);
    }
    float3  temp_90 = make_float3 (pSH17_9) * _S1121;
    float _S1208 = dot_0(temp_90, temp_90);
    if((F32_isfinite((_S1208))))
    {
        _S1123 = _S1208 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1209 = v_coeffs_25 + int(16);
        float _S1210 = atomicAdd(&(_S1209->x), temp_90.x);
        float _S1211 = atomicAdd(&(_S1209->y), temp_90.y);
        float _S1212 = atomicAdd(&(_S1209->z), temp_90.z);
    }
    float3  temp_91 = make_float3 (pSH18_9) * _S1121;
    float _S1213 = dot_0(temp_91, temp_91);
    if((F32_isfinite((_S1213))))
    {
        _S1123 = _S1213 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1214 = v_coeffs_25 + int(17);
        float _S1215 = atomicAdd(&(_S1214->x), temp_91.x);
        float _S1216 = atomicAdd(&(_S1214->y), temp_91.y);
        float _S1217 = atomicAdd(&(_S1214->z), temp_91.z);
    }
    float3  temp_92 = make_float3 (pSH19_9) * _S1121;
    float _S1218 = dot_0(temp_92, temp_92);
    if((F32_isfinite((_S1218))))
    {
        _S1123 = _S1218 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1219 = v_coeffs_25 + int(18);
        float _S1220 = atomicAdd(&(_S1219->x), temp_92.x);
        float _S1221 = atomicAdd(&(_S1219->y), temp_92.y);
        float _S1222 = atomicAdd(&(_S1219->z), temp_92.z);
    }
    float3  temp_93 = make_float3 (pSH20_9) * _S1121;
    float _S1223 = dot_0(temp_93, temp_93);
    if((F32_isfinite((_S1223))))
    {
        _S1123 = _S1223 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1224 = v_coeffs_25 + int(19);
        float _S1225 = atomicAdd(&(_S1224->x), temp_93.x);
        float _S1226 = atomicAdd(&(_S1224->y), temp_93.y);
        float _S1227 = atomicAdd(&(_S1224->z), temp_93.z);
    }
    float3  temp_94 = make_float3 (pSH21_9) * _S1121;
    float _S1228 = dot_0(temp_94, temp_94);
    if((F32_isfinite((_S1228))))
    {
        _S1123 = _S1228 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1229 = v_coeffs_25 + int(20);
        float _S1230 = atomicAdd(&(_S1229->x), temp_94.x);
        float _S1231 = atomicAdd(&(_S1229->y), temp_94.y);
        float _S1232 = atomicAdd(&(_S1229->z), temp_94.z);
    }
    float3  temp_95 = make_float3 (pSH22_9) * _S1121;
    float _S1233 = dot_0(temp_95, temp_95);
    if((F32_isfinite((_S1233))))
    {
        _S1123 = _S1233 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1234 = v_coeffs_25 + int(21);
        float _S1235 = atomicAdd(&(_S1234->x), temp_95.x);
        float _S1236 = atomicAdd(&(_S1234->y), temp_95.y);
        float _S1237 = atomicAdd(&(_S1234->z), temp_95.z);
    }
    float3  temp_96 = make_float3 (pSH23_9) * _S1121;
    float _S1238 = dot_0(temp_96, temp_96);
    if((F32_isfinite((_S1238))))
    {
        _S1123 = _S1238 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1239 = v_coeffs_25 + int(22);
        float _S1240 = atomicAdd(&(_S1239->x), temp_96.x);
        float _S1241 = atomicAdd(&(_S1239->y), temp_96.y);
        float _S1242 = atomicAdd(&(_S1239->z), temp_96.z);
    }
    float3  temp_97 = make_float3 (pSH24_9) * _S1121;
    float _S1243 = dot_0(temp_97, temp_97);
    if((F32_isfinite((_S1243))))
    {
        _S1123 = _S1243 != 0.0f;
    }
    else
    {
        _S1123 = false;
    }
    if(_S1123)
    {
        float3  * _S1244 = v_coeffs_25 + int(23);
        float _S1245 = atomicAdd(&(_S1244->x), temp_97.x);
        float _S1246 = atomicAdd(&(_S1244->y), temp_97.y);
        float _S1247 = atomicAdd(&(_S1244->z), temp_97.z);
    }
    float fTmp0D_z_9 = -14.04997730255126953f * z2_22 + 2.00713968276977539f;
    float fTmp1C_z_9 = 6.62322282791137695f * z_32;
    float pSH22_x_9 = fTmp1C_16 * _S1098;
    float3  dir_n_39 = make_float3 (x_35, y_34, z_32);
    float3  v_dir_n_39 = make_float3 (v_x_23 + dot_0(_S1121, make_float3 (0.62583571672439575f * (fS2_16 + y_34 * fC2_x_9 + x_35 * fS2_x_9)) * *_S1112 + make_float3 (0.62583571672439575f * (fC2_16 + x_35 * fC2_x_9 - y_34 * fS2_x_9)) * *_S1120 + make_float3 (fTmp2B_16 * fS2_x_9) * *_S1113 + make_float3 (fTmp2B_16 * fC2_x_9) * *_S1119 + make_float3 (fTmp1C_16 * fS1_x_13) * *_S1114 + make_float3 (pSH22_x_9) * *_S1118 + make_float3 (fTmp0D_16) * *_S1117), v_y_23 + dot_0(_S1121, make_float3 (0.62583571672439575f * (x_35 * fS2_y_9 + fC2_16 + y_34 * fC2_y_9)) * *_S1112 + make_float3 (0.62583571672439575f * (x_35 * fC2_y_9 - fS2_16 - y_34 * fS2_y_9)) * *_S1120 + make_float3 (fTmp2B_16 * fS2_y_9) * *_S1113 + make_float3 (fTmp2B_16 * fC2_y_9) * *_S1119 + make_float3 (pSH22_x_9) * *_S1114 + make_float3 (fTmp1C_16 * fC1_y_13) * *_S1118 + make_float3 (fTmp0D_16) * *_S1115), v_z_23 + dot_0(_S1121, make_float3 (1.9843134880065918f * (pSH12_20 + z_32 * pSH12_z_9) + -1.00623059272766113f * pSH6_z_9) * *_S1116 + make_float3 (fTmp0D_z_9 * x_35) * *_S1117 + make_float3 (fTmp0D_z_9 * y_34) * *_S1115 + make_float3 (fTmp1C_z_9 * fC1_22) * *_S1118 + make_float3 (fTmp1C_z_9 * fS1_22) * *_S1114 + make_float3 (-1.77013075351715088f * fC2_16) * *_S1119 + make_float3 (-1.77013075351715088f * fS2_16) * *_S1113));
    float3  v_viewdir_27 = v_viewdir_26 + (v_dir_n_39 - make_float3 (dot_0(v_dir_n_39, dir_n_39)) * dir_n_39) * make_float3 (inv_norm_24);
    Matrix<float, 3, 3>  _S1248 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1249;
    (&_S1249)->primal_0 = _S1090;
    (&_S1249)->differential_0 = _S1248;
    float3  _S1250 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1251;
    (&_S1251)->primal_0 = t_23;
    (&_S1251)->differential_0 = _S1250;
    s_bwd_prop_mul_0(&_S1249, &_S1251, v_viewdir_27);
    Matrix<float, 3, 3>  _S1252 = transpose_0(_S1249.differential_0);
    *v_mean_17 = *v_mean_17 + v_viewdir_27;
    *v_R_17 = *v_R_17 + _S1252;
    *v_t_17 = *v_t_17 + _S1251.differential_0;
    return;
}

