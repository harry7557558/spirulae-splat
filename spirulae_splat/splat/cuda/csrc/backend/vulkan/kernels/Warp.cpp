// Vulkan implementation of the dataset GT warp + byte->float conversion
// launch APIs (the launch_warp_* section of csrc/PixelWise.cuh and the
// uint8/16_*_to_float_raw converters of EngineInternal.h). Mirrors the
// CUDA launchers in csrc/ImageWarp.cu; the device work runs
// slang/vulkan/warp.slang.
//
// All launchers take raw device pointers (the engine's DataManager path
// drives them without tensor wrappers). Nullable pointers (dist_coeffs)
// go through vkk::or_fallback; input type tags map to the shader's
// elem_kind field instead of the CUDA template dispatch.

#include <PixelWise.cuh>
#include <EngineInternal.h>
#include <Common.cuh>

#include "KernelCommon.h"

namespace {

// Mirrors the elem-kind constants in slang/vulkan/warp.slang.
constexpr int kElemU8 = 0;
constexpr int kElemU16 = 1;
constexpr int kElemF32 = 2;

// Mirrors WarpParams in slang/vulkan/warp.slang.
struct WarpParams {
    uint64_t intrins, dist_coeffs, src, dst, axes;
    int32_t B, Hin, Win, C;
    int32_t K, Hout, Wout;
    int32_t elem_kind, camera_model;
    int32_t in_H, in_W;
    int32_t ray_depth;
    float norm_inv, decode_off;
};
static_assert(sizeof(WarpParams) == 5 * 8 + 14 * 4,
              "params layout must match the slang struct");

// Mirrors WarpMaskParams in slang/vulkan/warp.slang.
struct WarpMaskParams {
    uint64_t intrins, dist_coeffs, src, dst, axes;
    int32_t B, Hin, Win;
    int32_t K, Hout, Wout;
    int32_t camera_model;
    uint32_t wgs_per_row;
};
static_assert(sizeof(WarpMaskParams) == 5 * 8 + 8 * 4,
              "params layout must match the slang struct");

// Mirrors BytesToFloatParams in slang/vulkan/warp.slang.
struct BytesToFloatParams {
    uint64_t src, dst;
    uint32_t total;
    int32_t elem_kind;
    float scale, offset;
    uint32_t wgs_per_row;
    int32_t _pad0;
};
static_assert(sizeof(BytesToFloatParams) == 2 * 8 + 6 * 4,
              "params layout must match the slang struct");

// Common WarpParams assembly for the image / depth / normal warps.
WarpParams make_warp_params(const float* d_intrins, const float* d_dist,
                            const void* d_src, float* d_dst,
                            const float* d_axes, int B, int Hin, int Win,
                            int C, int K, int Hout, int Wout, int elem_kind,
                            CameraModelType cm, int in_H, int in_W,
                            bool ray_depth, float norm_inv,
                            float decode_off) {
    WarpParams p{};
    p.intrins = vkk::or_fallback(d_intrins);
    p.dist_coeffs = vkk::or_fallback(d_dist);
    p.src = vkk::or_fallback(d_src);
    p.dst = (uint64_t)d_dst;
    p.axes = vkk::or_fallback(d_axes);
    p.B = B; p.Hin = Hin; p.Win = Win; p.C = C;
    p.K = K; p.Hout = Hout; p.Wout = Wout;
    p.elem_kind = elem_kind;
    p.camera_model = (int)cm;
    p.in_H = in_H; p.in_W = in_W;
    p.ray_depth = ray_depth ? 1 : 0;
    p.norm_inv = norm_inv;
    p.decode_off = decode_off;
    return p;
}

// (ceil(Wout/16), ceil(Hout/16), B) grid of flat-256 workgroups (16x16
// logical tiles), mirroring the CUDA _LAUNCH_ARGS_3D(Wout, Hout, B, 16, 16, 1).
void dispatch_warp(const char* entry, const WarpParams& p) {
    if ((int64_t)p.B * p.K * p.Hout * p.Wout <= 0) return;
    if (p.B > 65535)
        throw std::runtime_error("warp: batch exceeds grid limit");
    vkk::dispatch(entry, {}, (uint32_t)((p.Wout + 15) / 16),
                  (uint32_t)((p.Hout + 15) / 16), (uint32_t)p.B, &p,
                  sizeof(p));
}

void run_bytes_to_float(const void* d_in, float* d_out, int64_t total,
                        int elem_kind, float scale, float offset) {
    BytesToFloatParams p{};
    p.src = vkk::or_fallback(d_in);
    p.dst = (uint64_t)d_out;
    p.total = (uint32_t)total;
    p.elem_kind = elem_kind;
    p.scale = scale;
    p.offset = offset;
    vkk::dispatch_flat("warp.bytes_to_float", {}, total, 256, &p, sizeof(p),
                       &p.wgs_per_row);
}

int resolve_depth_kind(uint32_t elem_size, const char* who) {
    if (elem_size == 2) return kElemU16;
    if (elem_size == 4) return kElemF32;
    throw std::runtime_error(std::string(who) +
                             ": depth must be uint16 or float32");
}

}  // namespace


// ---- Raw byte -> float converters (EngineInternal.h) -----------------------

void uint8_image_to_float_raw(const uint8_t* d_in, float* d_out, int B, int H,
                              int W, int C) {
    run_bytes_to_float(d_in, d_out, (int64_t)B * H * W * C, kElemU8,
                       1.0f / 255.0f, 0.0f);
}

void uint16_image_to_float_raw(const uint16_t* d_in, float* d_out, int B,
                               int H, int W, int C) {
    run_bytes_to_float(d_in, d_out, (int64_t)B * H * W * C, kElemU16,
                       1.0f / 65535.0f, 0.0f);
}

void uint8_normal_to_float_raw(const uint8_t* d_in, float* d_out, int B,
                               int H, int W, int C) {
    run_bytes_to_float(d_in, d_out, (int64_t)B * H * W * C, kElemU8,
                       1.0f / 127.5f, -1.0f);
}

void uint16_depth_to_float_raw(const uint16_t* d_in, float* d_out, int B,
                               int H, int W, int C) {
    run_bytes_to_float(d_in, d_out, (int64_t)B * H * W * C, kElemU16, 1.0f,
                       0.0f);
}

void uint8_image_to_float_tensor(DeviceTensor3D<uint8_t> img_in,
                                 DeviceTensor3D<float> img_out) {
    run_bytes_to_float(img_in.data_ptr(), img_out.data_ptr(),
                       img_in.size<0>() * img_in.size<1>() * img_in.size<2>(),
                       kElemU8, 1.0f / 255.0f, 0.0f);
}

void uint16_image_to_float_tensor(DeviceTensor3D<uint16_t> img_in,
                                  DeviceTensor3D<float> img_out) {
    run_bytes_to_float(img_in.data_ptr(), img_out.data_ptr(),
                       img_in.size<0>() * img_in.size<1>() * img_in.size<2>(),
                       kElemU16, 1.0f / 65535.0f, 0.0f);
}


// ---- Fused byte -> float GT warps (PixelWise.cuh) ---------------------------

void launch_warp_byte_to_float_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const void* d_byte, bool input_is_u16,
    int B, int Hin, int Win, int C,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    dispatch_warp("warp.warp_img_wide",
                  make_warp_params(
                      d_intrins, d_dist_coeffs, d_byte, d_float_out, d_axes,
                      B, Hin, Win, C, K, Hout, Wout,
                      input_is_u16 ? kElemU16 : kElemU8, cmt(camera_model),
                      Hin, Win, false,
                      input_is_u16 ? 1.0f / 65535.0f : 1.0f / 255.0f, 0.0f));
}

void launch_warp_byte_to_float_equi(
    const void* d_byte, bool input_is_u16,
    int B, int Hin, int Win, int C,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    dispatch_warp("warp.warp_img_equi",
                  make_warp_params(
                      nullptr, nullptr, d_byte, d_float_out, d_axes,
                      B, Hin, Win, C, K, Hout, Wout,
                      input_is_u16 ? kElemU16 : kElemU8,
                      CameraModelType::EQUIRECTANGULAR, Hin, Win, false,
                      input_is_u16 ? 1.0f / 65535.0f : 1.0f / 255.0f, 0.0f));
}

void launch_warp_mask_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const uint8_t* d_byte_mask,
    int B, int Hin, int Win,
    uint8_t* d_byte_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    WarpMaskParams p{};
    p.intrins = vkk::or_fallback(d_intrins);
    p.dist_coeffs = vkk::or_fallback(d_dist_coeffs);
    p.src = vkk::or_fallback(d_byte_mask);
    p.dst = (uint64_t)d_byte_out;
    p.axes = vkk::or_fallback(d_axes);
    p.B = B; p.Hin = Hin; p.Win = Win;
    p.K = K; p.Hout = Hout; p.Wout = Wout;
    p.camera_model = (int)cmt(camera_model);
    int64_t words = ((int64_t)B * K * Hout * Wout + 3) / 4;
    vkk::dispatch_flat("warp.warp_mask_wide", {}, words, 256, &p, sizeof(p),
                       &p.wgs_per_row);
}

void launch_warp_mask_equi(
    const uint8_t* d_byte_mask,
    int B, int Hin, int Win,
    uint8_t* d_byte_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    WarpMaskParams p{};
    p.intrins = vkk::or_fallback(nullptr);
    p.dist_coeffs = vkk::or_fallback(nullptr);
    p.src = vkk::or_fallback(d_byte_mask);
    p.dst = (uint64_t)d_byte_out;
    p.axes = vkk::or_fallback(d_axes);
    p.B = B; p.Hin = Hin; p.Win = Win;
    p.K = K; p.Hout = Hout; p.Wout = Wout;
    p.camera_model = (int)CameraModelType::EQUIRECTANGULAR;
    int64_t words = ((int64_t)B * K * Hout * Wout + 3) / 4;
    vkk::dispatch_flat("warp.warp_mask_equi", {}, words, 256, &p, sizeof(p),
                       &p.wgs_per_row);
}

void launch_warp_depth_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const void* d_depth, uint32_t elem_size,
    int B, int Hin, int Win,
    int in_H, int in_W,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes, bool input_is_ray_depth)
{
    dispatch_warp("warp.warp_depth_wide",
                  make_warp_params(
                      d_intrins, d_dist_coeffs, d_depth, d_float_out, d_axes,
                      B, Hin, Win, 1, K, Hout, Wout,
                      resolve_depth_kind(elem_size, "launch_warp_depth_wide"),
                      cmt(camera_model), in_H, in_W, input_is_ray_depth,
                      1.0f, 0.0f));
}

void launch_warp_depth_equi(
    const void* d_depth, uint32_t elem_size,
    int B, int Hin, int Win,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes, bool input_is_ray_depth)
{
    dispatch_warp("warp.warp_depth_equi",
                  make_warp_params(
                      nullptr, nullptr, d_depth, d_float_out, d_axes,
                      B, Hin, Win, 1, K, Hout, Wout,
                      resolve_depth_kind(elem_size, "launch_warp_depth_equi"),
                      CameraModelType::EQUIRECTANGULAR, Hin, Win,
                      input_is_ray_depth, 1.0f, 0.0f));
}

void launch_warp_normal_wide(
    std::string camera_model,
    const float* d_intrins,
    const float* d_dist_coeffs,
    const void* d_normal, uint32_t elem_size,
    int B, int Hin, int Win,
    int in_H, int in_W,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    if (elem_size != 1 && elem_size != 4)
        throw std::runtime_error(
            "launch_warp_normal_wide: normal must be uint8 or float32");
    bool u8 = elem_size == 1;
    dispatch_warp("warp.warp_normal_wide",
                  make_warp_params(
                      d_intrins, d_dist_coeffs, d_normal, d_float_out, d_axes,
                      B, Hin, Win, 3, K, Hout, Wout,
                      u8 ? kElemU8 : kElemF32, cmt(camera_model), in_H, in_W,
                      false, u8 ? 1.0f / 127.5f : 1.0f, u8 ? -1.0f : 0.0f));
}

void launch_warp_normal_equi(
    const void* d_normal, uint32_t elem_size,
    int B, int Hin, int Win,
    float* d_float_out, int K, int Hout, int Wout,
    const float* d_axes)
{
    if (elem_size != 1 && elem_size != 4)
        throw std::runtime_error(
            "launch_warp_normal_equi: normal must be uint8 or float32");
    bool u8 = elem_size == 1;
    dispatch_warp("warp.warp_normal_equi",
                  make_warp_params(
                      nullptr, nullptr, d_normal, d_float_out, d_axes,
                      B, Hin, Win, 3, K, Hout, Wout,
                      u8 ? kElemU8 : kElemF32,
                      CameraModelType::EQUIRECTANGULAR, Hin, Win, false,
                      u8 ? 1.0f / 127.5f : 1.0f, u8 ? -1.0f : 0.0f));
}
