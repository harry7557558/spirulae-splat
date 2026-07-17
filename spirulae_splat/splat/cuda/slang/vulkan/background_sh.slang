// Background spherical-harmonics forward (Vulkan backend), mirroring the
// forward kernel in csrc/BackgroundSphericalHarmonics.cu: one thread per
// pixel, evaluate SH along the world-space ray direction. Ray generation and
// the camera-to-world rotation use the same projection_utils.slang exports
// the CUDA kernel wraps (generate_ray + transform_ray_d).
//
// kShDegree is a specialization constant (the CUDA kernel's template
// parameter); the camera model is a plain runtime value, exactly like the
// CUDA kernel's `int camera_model` argument.
//
// The backward (block-atomic SH gradient reduction) lands with the training
// phase.

#define CudaDeviceExport ForceInline
#include "projection_utils.slang"
#include "harmonics.slang"
#undef CudaDeviceExport

[SpecializationConstant]
const int kShDegree = 0;

struct BackgroundShFwdParams {  // 56 bytes: pushed directly
    float* viewmats;     // [B,16] row-major world->camera
    float4* intrins;     // [B] fx, fy, cx, cy
    float* dist_coeffs;  // [B,10] or null
    float* sh_coeffs;    // [(kShDegree+1)^2, 3]
    float* out_img;      // [B,H,W,3]
    uint32_t width;
    uint32_t height;
    uint32_t B;
    int camera_model;
};

// Workgroups are a flat 128x1 pixel-row segment (not the CUDA kernel's
// 16x16): REQUIRE_FULL_SUBGROUPS demands the X dimension be a multiple of
// the pinned subgroup size (README subgroup rules). Grid (ceil(W/128), H, B).
[shader("compute")]
[numthreads(128, 1, 1)]
void background_sh_fwd(uint3 gtid: SV_DispatchThreadID,
                       uniform BackgroundShFwdParams p) {
    uint px = gtid.x, py = gtid.y, bi = gtid.z;
    if (px >= p.width || py >= p.height || bi >= p.B) return;
    uint pix_id = (bi * p.height + py) * p.width + px;

    // Per-batch camera (indexed directly by bi; see the CUDA kernel's note
    // on why this is NOT a global camera id).
    float* vm = p.viewmats + bi * 16;
    float3x3 R;
    [ForceUnroll]
    for (int r = 0; r < 3; r++)
        R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
    float4 intrin = p.intrins[bi];
    float dist[10];
    [ForceUnroll]
    for (int i = 0; i < 10; i++)
        dist[i] = p.dist_coeffs == nullptr ? 0.0f
                                           : p.dist_coeffs[10 * bi + i];

    float2 uv = float2((float(px) + 0.5f - intrin.z) / intrin.x,
                       (float(py) + 0.5f - intrin.w) / intrin.y);

    float3 color = float3(0.0f);
    float3 raydir;
    if (generate_ray(uv, p.camera_model, dist, raydir)) {
        float3 dir = transform_ray_d(R, raydir);
        float3 dc = float3(p.sh_coeffs[0], p.sh_coeffs[1], p.sh_coeffs[2]);
        float3* coeffs_l1 = (float3*)p.sh_coeffs + 1;
        switch (kShDegree) {
        case 0: color = sh0_to_color_dir(dir, dc, coeffs_l1); break;
        case 1: color = sh1_to_color_dir(dir, dc, coeffs_l1); break;
        case 2: color = sh2_to_color_dir(dir, dc, coeffs_l1); break;
        case 3: color = sh3_to_color_dir(dir, dc, coeffs_l1); break;
        case 4: color = sh4_to_color_dir(dir, dc, coeffs_l1); break;
        }
        color = max(color + 0.5f, float3(0.0f));
    }
    p.out_img[3 * pix_id + 0] = color.x;
    p.out_img[3 * pix_id + 1] = color.y;
    p.out_img[3 * pix_id + 2] = color.z;
}
