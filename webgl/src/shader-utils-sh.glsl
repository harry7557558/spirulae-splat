#line 2

uniform ivec2 u_sh_config;

uniform highp usampler2D u_sh_texture;

#define use_sh bool(u_sh_config.x)
#define sh_degree u_sh_config.y
#define sh_dim (sh_degree * (sh_degree + 2))

ivec3 pix_sh;
uvec4 coeff_sh;

void init_coeffs(int splat_idx) {
    int width = textureSize(u_sh_texture, 0).x;
    int pixidx = splat_idx * ((sh_dim+1)/2);
    pix_sh = ivec3(pixidx % width, pixidx / width, 0);
}

vec3 next_coeff() {
    uvec2 packed;
    if (pix_sh.z == 0) {
        coeff_sh = texelFetch(u_sh_texture, pix_sh.xy, 0);
        packed = coeff_sh.xy;
        pix_sh.z += 1;
    }
    else {
        packed = coeff_sh.zw;
        pix_sh.x += 1;
        pix_sh.z = 0;
    }
    vec2 c1 = unpackHalf2x16(packed.x);
    vec2 c2 = unpackHalf2x16(packed.y);
    return vec3(c1, c2.x);
}

vec3 sh_coeffs_to_color_fast(
    vec3 base_color,
    int degree,
    vec3 viewdir
) {
    vec3 color = 0.2820947917738781 * base_color;
    if (degree < 1 || !use_sh) {
        return color;
    }

    float norm = length(viewdir);
    float x = viewdir.x / norm;
    float y = viewdir.y / norm;
    float z = viewdir.z / norm;

    float fTmp0A = 0.48860251190292;
    color += fTmp0A * (-y * next_coeff());
    color += fTmp0A * (z * next_coeff());
    color += fTmp0A * (-x * next_coeff());
    if (degree < 2) {
        return color;
    }
    float z2 = z * z;

    float fTmp0B = -1.092548430592079 * z;
    float fTmp1A = 0.5462742152960395;
    float fC1 = x * x - y * y;
    float fS1 = 2. * x * y;
    float pSH6 = (0.9461746957575601 * z2 - 0.3153915652525201);
    float pSH7 = fTmp0B * x;
    float pSH5 = fTmp0B * y;
    float pSH8 = fTmp1A * fC1;
    float pSH4 = fTmp1A * fS1;
    color += pSH4 * next_coeff();
    color += pSH5 * next_coeff();
    color += pSH6 * next_coeff();
    color += pSH7 * next_coeff();
    color += pSH8 * next_coeff();
    if (degree < 3) {
        return color;
    }

    float fTmp0C = -2.285228997322329 * z2 + 0.4570457994644658;
    float fTmp1B = 1.445305721320277 * z;
    float fTmp2A = -0.5900435899266435;
    float fC2 = x * fC1 - y * fS1;
    float fS2 = x * fS1 + y * fC1;
    float pSH12 = z * (1.865881662950577 * z2 - 1.119528997770346);
    float pSH13 = fTmp0C * x;
    float pSH11 = fTmp0C * y;
    float pSH14 = fTmp1B * fC1;
    float pSH10 = fTmp1B * fS1;
    float pSH15 = fTmp2A * fC2;
    float pSH9  = fTmp2A * fS2;
    color += pSH9 * next_coeff();
    color += pSH10 * next_coeff();
    color += pSH11 * next_coeff();
    color += pSH12 * next_coeff();
    color += pSH13 * next_coeff();
    color += pSH14 * next_coeff();
    color += pSH15 * next_coeff();
    if (degree < 4) {
        return color;
    }

    float fTmp0D = z * (-4.683325804901025 * z2 + 2.007139630671868);
    float fTmp1C = 3.31161143515146 * z2 - 0.47308734787878;
    float fTmp2B = -1.770130769779931 * z;
    float fTmp3A = 0.6258357354491763;
    float fC3 = x * fC2 - y * fS2;
    float fS3 = x * fS2 + y * fC2;
    float pSH20 = (1.984313483298443 * z * pSH12 - 1.006230589874905 * pSH6);
    float pSH21 = fTmp0D * x;
    float pSH19 = fTmp0D * y;
    float pSH22 = fTmp1C * fC1;
    float pSH18 = fTmp1C * fS1;
    float pSH23 = fTmp2B * fC2;
    float pSH17 = fTmp2B * fS2;
    float pSH24 = fTmp3A * fC3;
    float pSH16 = fTmp3A * fS3;
    color += pSH16 * next_coeff();
    color += pSH17 * next_coeff();
    color += pSH18 * next_coeff();
    color += pSH19 * next_coeff();
    color += pSH20 * next_coeff();
    color += pSH21 * next_coeff();
    color += pSH22 * next_coeff();
    color += pSH23 * next_coeff();
    color += pSH24 * next_coeff();

    return color;
}
