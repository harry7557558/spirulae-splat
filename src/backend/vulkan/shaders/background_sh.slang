// Background spherical-harmonics forward (Vulkan backend), mirroring the
// forward kernel in kernels/background/BackgroundSphericalHarmonics.cu: one thread per
// pixel, evaluate SH along the world-space ray direction. Ray generation and
// the camera-to-world rotation use the same projection_utils.slang exports
// the CUDA kernel wraps (generate_ray + transform_ray_d).
//
// kShDegree and kDistortion (declared last) are specialization constants (the
// CUDA kernel's template parameters); the camera model is a plain runtime
// value, exactly like the CUDA kernel's `int camera_model` argument.
//
// Backward: CUDA's 512-thread block-atomic SH reduction
// (sh_block_atomic_add_f3_512) hardcodes warp==32 and needs float atomics,
// so the port replaces it with a fixed grid of kBgBwdSlices grid-striding
// threads, each accumulating its own scratch slice via the plain
// sh_coeffs_to_color_vjp_inplace (no atomics at all), followed by a small
// deterministic reduce pass that folds the slices into v_sh_coeffs.
//
// Coefficients are addressed as flat float*, not float3*, and the sh{N}_*_dir
// wrappers -- which only exist for float3* -- are bypassed for the underlying
// sh_coeffs_to_color(_vjp_inplace) they forward to. MoltenVK's SPIR-V -> MSL
// translation cannot do arithmetic on a packed_float3 physical-storage
// pointer: it emits a reinterpret_cast of the loaded packed_float3 to the
// pointer's ulong, which the Metal compiler rejects, and every background
// pipeline fails to create.

#define CudaDeviceExport ForceInline
#include "int64_compat.slang"
#include "shaders/projection_utils.slang"
#include "shaders/harmonics.slang"
#undef CudaDeviceExport
#include "dist_spec.slang"

[SpecializationConstant]
const int kShDegree = 0;

struct BackgroundShFwdParams {  // 56 bytes: pushed directly
    float* viewmats;     // [B,16] row-major world->camera
    float4* intrins;     // [B] fx, fy, cx, cy
    float* dist_coeffs;  // [B,8] or null
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
#define _bg_ray(D) \
    ray_ok = pixel_ray<D>(p.dist_coeffs, bi, p.camera_model, uv, raydir)

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

    float2 uv = float2((float(px) + 0.5f - intrin.z) / intrin.x,
                       (float(py) + 0.5f - intrin.w) / intrin.y);

    float3 color = float3(0.0f);
    float3 raydir = float3(0.0f);
    bool ray_ok = false;
    SS_DISPATCH_DIST(_bg_ray)
    if (ray_ok) {
        float3 dir = transform_ray_d(R, raydir);
        float3 dc = float3(p.sh_coeffs[0], p.sh_coeffs[1], p.sh_coeffs[2]);
        float* coeffs_l1 = p.sh_coeffs + 3;
        switch (kShDegree) {
        case 0: color = sh_coeffs_to_color(0, dir, dc, coeffs_l1); break;
        case 1: color = sh_coeffs_to_color(1, dir, dc, coeffs_l1); break;
        case 2: color = sh_coeffs_to_color(2, dir, dc, coeffs_l1); break;
        case 3: color = sh_coeffs_to_color(3, dir, dc, coeffs_l1); break;
        case 4: color = sh_coeffs_to_color(4, dir, dc, coeffs_l1); break;
        }
        color = max(color + 0.5f, float3(0.0f));
    }
    p.out_img[3 * pix_id + 0] = color.x;
    p.out_img[3 * pix_id + 1] = color.y;
    p.out_img[3 * pix_id + 2] = color.z;
}

// ---------------------------------------------------------------------------
// Backward
// ---------------------------------------------------------------------------

static const uint kBgBwdWgs = 64;
static const uint kBgBwdSlices = kBgBwdWgs * 256;

struct BackgroundShBwdParams {  // 80 bytes: pushed directly
    float* viewmats;     // [B,16]
    float4* intrins;     // [B]
    float* dist_coeffs;  // [B,8] (fallback-backed)
    float* sh_coeffs;    // [K,3], K = (kShDegree+1)^2
    float* out_color;    // [B,H,W,3] forward output (clamp mask)
    float* v_out_color;  // [B,H,W,3]
    float* scratch;      // [kBgBwdSlices, K, 3] zero-initialized partials
    uint32_t width;
    uint32_t height;
    uint32_t B;
    int camera_model;
};

[ForceInline]
bool _bg_finitef(float v) {
    return (asuint(v) & 0x7f800000u) != 0x7f800000u;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void background_sh_bwd(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                       uniform BackgroundShBwdParams p) {
    uint slice = wg.x * 256 + tid;
    const int K = (kShDegree + 1) * (kShDegree + 1);
    float* my = p.scratch + slice * K * 3;   // slot 0 = dc, 1..K-1 = bands
    float* my_coeffs = my + 3;
    uint total = p.B * p.height * p.width;

    float3 dc = float3(p.sh_coeffs[0], p.sh_coeffs[1], p.sh_coeffs[2]);
    float* coeffs_l1 = p.sh_coeffs + 3;

    float3 v_dc = float3(0.0f);
    for (uint pix = slice; pix < total; pix += kBgBwdSlices) {
        uint px = pix % p.width;
        uint rem = pix / p.width;
        uint py = rem % p.height;
        uint bi = rem / p.height;

        // dL/d(post-clamp) -> dL/d(pre-clamp): zero where the clamp kicked in.
        float3 color = float3(p.out_color[3 * pix], p.out_color[3 * pix + 1],
                              p.out_color[3 * pix + 2]);
        float3 v_color =
            float3(p.v_out_color[3 * pix], p.v_out_color[3 * pix + 1],
                   p.v_out_color[3 * pix + 2]);
        if (color.x == 0.0f || !_bg_finitef(v_color.x)) v_color.x = 0.0f;
        if (color.y == 0.0f || !_bg_finitef(v_color.y)) v_color.y = 0.0f;
        if (color.z == 0.0f || !_bg_finitef(v_color.z)) v_color.z = 0.0f;

        float* vm = p.viewmats + bi * 16;
        float3x3 R;
        [ForceUnroll]
        for (int r = 0; r < 3; r++)
            R[r] = float3(vm[4 * r + 0], vm[4 * r + 1], vm[4 * r + 2]);
        float4 intrin = p.intrins[bi];
        float2 uv = float2((float(px) + 0.5f - intrin.z) / intrin.x,
                           (float(py) + 0.5f - intrin.w) / intrin.y);
        float3 raydir = float3(0.0f);
        bool ray_ok = false;
        SS_DISPATCH_DIST(_bg_ray)
        if (!ray_ok) continue;
        float3 dir = transform_ray_d(R, raydir);

        // v_dir is discarded (engine holds cameras fixed, as in CUDA).
        float3 v_dir = float3(0.0f);
        switch (kShDegree) {
        case 0: sh_coeffs_to_color_vjp_inplace(0, dir, dc, coeffs_l1, v_color,
                                               v_dc, my_coeffs, v_dir); break;
        case 1: sh_coeffs_to_color_vjp_inplace(1, dir, dc, coeffs_l1, v_color,
                                               v_dc, my_coeffs, v_dir); break;
        case 2: sh_coeffs_to_color_vjp_inplace(2, dir, dc, coeffs_l1, v_color,
                                               v_dc, my_coeffs, v_dir); break;
        case 3: sh_coeffs_to_color_vjp_inplace(3, dir, dc, coeffs_l1, v_color,
                                               v_dc, my_coeffs, v_dir); break;
        case 4: sh_coeffs_to_color_vjp_inplace(4, dir, dc, coeffs_l1, v_color,
                                               v_dc, my_coeffs, v_dir); break;
        }
    }
    my[0] = v_dc.x;
    my[1] = v_dc.y;
    my[2] = v_dc.z;
}

struct BackgroundShBwdReduceParams {  // 20 bytes: pushed directly
    float* scratch;       // [kBgBwdSlices, K, 3], read flat
    float* v_sh_coeffs;   // [K,3] accumulated (caller zeroes)
    uint32_t k3;          // K * 3
};

groupshared float _bgsh_red[256];

// One workgroup PER output scalar (grid.x = k3); the block's 256 threads
// grid-stride the kBgBwdSlices partials for that column, then a groupshared
// tree folds them to one value. No atomics: each output column is owned by
// exactly one workgroup, so the write is a plain accumulate and the result
// is subgroup-size independent. Replaces the old single-workgroup kernel
// that folded all 16384 slices in one 48-wide serial loop (~1/60 of the GPU,
// 2.4 ms on RADV); this parallelizes across columns AND slices.
[shader("compute")]
[numthreads(256, 1, 1)]
void background_sh_bwd_reduce(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                              uniform BackgroundShBwdReduceParams p) {
    uint col = wg.x;
    if (col >= p.k3) return;  // whole workgroup shares col: uniform, no barrier split
    float s = 0.0f;
    for (uint i = tid; i < kBgBwdSlices; i += 256)
        s += p.scratch[i * p.k3 + col];

    GroupMemoryBarrierWithGroupSync();
    _bgsh_red[tid] = s;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll]
    for (uint st = 128; st > 0; st >>= 1) {
        if (tid < st) _bgsh_red[tid] += _bgsh_red[tid + st];
        GroupMemoryBarrierWithGroupSync();
    }
    if (tid == 0) p.v_sh_coeffs[col] += _bgsh_red[0];
}

#undef _bg_ray

// Declared last: spec IDs follow declaration order (see dist_spec.slang).
[SpecializationConstant]
const int kDistortion = 0;
