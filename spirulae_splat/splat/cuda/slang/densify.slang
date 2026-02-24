// suppress CudaDeviceExport to reduce code size
#define CudaDeviceExport Differentiable
#include "projection_utils.slang"
#undef CudaDeviceExport

[CudaDeviceExport]
void mcmc_add_noise_3dgs(
    float scaler, float min_opacity,
    inout float3 mean, float3 scale, float4 quat, float opac
) {
    float k = 0.5 / min_opacity;
    float x0 = 1.0 - min_opacity;
    float op_sigmoid = 1.0 / (1.0 + exp(-k * ((1.0 - opac) - x0)));

    float3x3 covar;
    quat_scale_to_covar(normalize(quat), exp(scale), covar);

    float3 noise = randn3(float4(dot(mean, mean), dot(mean, scale), dot(scale, scale), dot(quat, float4(opac))));
    noise = noise * scaler * op_sigmoid;
    mean += mul(covar, noise);
}

[CudaDeviceExport]
void mcmc_add_noise_triangle(
    float scaler, float min_opacity,
    inout float3 mean, float3 scale, float4 quat, float opac
) {
    float k = 0.5 / min_opacity;
    float x0 = 1.0 - min_opacity;
    float op_sigmoid = 1.0 / (1.0 + exp(-k * ((1.0 - opac) - x0)));

    float3 vert0, vert1, vert2;
    map_opaque_triangle(mean, normalize(quat), scale, vert0, vert1, vert2);

    float3 vertc = (vert0 + vert1 + vert2) / 3.0f;
    float3 d0 = vert0 - vertc, d1 = vert1 - vertc, d2 = vert2 - vertc;
    float3 dn = 0.5f * min(min(length(d0), length(d1)), length(d2)) * normalize(cross(d0, d1));
    float3x3 covar = 0.5f * (outerProduct(d0) + outerProduct(d1) + outerProduct(d2) + outerProduct(dn)) / 3.5f;

    float3 noise = randn3(float4(dot(mean, mean), dot(mean, scale), dot(scale, scale), dot(quat, float4(opac))));
    noise = noise * scaler * op_sigmoid;
    mean += mul(covar, noise);
}

[CudaDeviceExport]
void long_axis_split_3dgs(
    float3 scale, float4 quat,
    out float3 new_scale, out float3 mean_delta
) {
    // https://arxiv.org/abs/2508.12313
    float l0 = max(max(scale.x, scale.y), scale.z);
    float d = 0.45 * l0;
    new_scale = scale * sqrt(1.0f - d * d / (l0 * l0));
    mean_delta = float3(0, 0, 0);
    if (scale.x > scale.y && scale.x > scale.z)
    { new_scale.x = l0 - d; mean_delta.x = d; }
    else if (scale.y > scale.z)
    { new_scale.y = l0 - d; mean_delta.y = d; }
    else
    { new_scale.z = l0 - d; mean_delta.z = d; }
    mean_delta = mul(normalized_quat_to_rotmat(normalize(quat)), mean_delta);
}
