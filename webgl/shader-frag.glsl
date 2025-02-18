#version 300 es
precision highp float;
precision highp int;

uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;
uniform vec4 undistortion;

uniform ivec3 u_ch_config;

flat in vec4 vColor;
flat in vec3 vPosition;
flat in vec3 vAxesU;
flat in vec3 vAxesV;
flat in vec2 vAnisotropy;
flat in int vIndex;

uniform highp usampler2D u_ch_texture;
// uniform highp usampler2D u_sh_texture;

out vec4 fragColor;

#define PI 3.14159265


float bessel_j0(float x) {
    x = abs(x);
    const float A01 = -0.22334141550835468;
    const float A02 = -0.18054514613169334;
    const float A03 = 0.047655819492146555;
    const float A04 = -0.0024383224605644127;
    const float B01 = -0.2234104067240744;
    const float B02 = 0.06997488447829185;
    const float B03 = -0.009545108494228093;
    const float B04 = 0.0011020589286710896;
    const float C01 = -0.12645960753014673;
    const float C02 = -0.046982403187758204;
    const float D01 = 0.12546821167569544;
    const float D02 = -0.08179017949118085;
    if (x < 4.2) {
        return (1.0 + x * (A01 + x * (A02 + x * (A03 + x * A04)))) / 
               (1.0 + x * (B01 + x * (B02 + x * (B03 + x * B04))));
    } else {
        float inv_x = 1.0 / x;
        return sqrt(inv_x / PI) * 
               ((1.0 + inv_x * (C01 + inv_x * C02)) * cos(x) + 
                (1.0 + inv_x * (D01 + inv_x * D02)) * sin(x));
    }
}

float bessel_j1(float x) {
    float s = sign(x); x = abs(x);
    const float A11 = 0.5000532505934573;
    const float A12 = -0.07119189008625842;
    const float A13 = -0.03516544965310418;
    const float A14 = 0.005137712441402014;
    const float B11 = -0.14169754004287138;
    const float B12 = 0.05321374041943971;
    const float B13 = -0.006074191793869077;
    const float B14 = 0.0008890431150018836;
    const float C11 = -0.37578090667550257;
    const float C12 = 0.13415846822338878;
    const float D11 = 0.3771741874195154;
    const float D12 = 0.08328593955487182;
    if (x < 4.7) {
        return s * (x * (A11 + x * (A12 + x * (A13 + x * A14)))) / 
               (1.0 + x * (B11 + x * (B12 + x * (B13 + x * B14))));
    } else {
        float inv_x = 1.0 / x;
        return -s * sqrt(inv_x / PI) * 
               ((1.0 + inv_x * (C11 + inv_x * C12)) * cos(x) - 
                (1.0 + inv_x * (D11 + inv_x * D12)) * sin(x));
    }
}

float bessel_j2(float x) {
    x = abs(x);
    const float A22 = 0.12515704623184004;
    const float A23 = -0.029939064720165425;
    const float A24 = 0.0010672441356451344;
    const float B21 = -0.23414231622686957;
    const float B22 = 0.08321793931659045;
    const float B23 = -0.012670220025970099;
    const float B24 = 0.0015767563111494629;
    const float C21 = 1.874412399273724;
    const float C22 = -0.8101992991186221;
    const float C23 = 0.3015954731134034;
    const float D21 = -1.871274305839778;
    const float D22 = -0.8861009908575821;
    if (x < 4.0) {
        return (x * x * (A22 + x * (A23 + x * A24))) / 
               (1.0 + x * (B21 + x * (B22 + x * (B23 + x * B24))));
    } else {
        float inv_x = 1.0 / x;
        return -sqrt(inv_x / PI) * 
               ((1.0 + inv_x * (C21 + inv_x * (C22 + inv_x * C23))) * cos(x) + 
                (1.0 + inv_x * (D21 + inv_x * D22)) * sin(x));
    }
}

float bessel_j3(float x) {
    float s = sign(x); x = abs(x);
    const float A33 = 0.020910819133348472;
    const float A34 = -0.005072094376038419;
    const float A35 = 0.0002802765938927514;
    const float B31 = -0.2330616229268751;
    const float B32 = 0.06455328873550212;
    const float B33 = -0.008312028977714348;
    const float B34 = 0.0007466861514973682;
    const float C31 = -4.376749965939438;
    const float C32 = -7.327544311795212;
    const float C33 = 2.8595505732173425;
    const float D31 = 4.374149309521666;
    const float D32 = -7.3507982430716545;
    const float D33 = -3.7324735035522663;
    if (x < 5.0) {
        return s * (x * x * x * (A33 + x * (A34 + x * A35))) / 
               (1.0 + x * (B31 + x * (B32 + x * (B33 + x * B34))));
    } else {
        float inv_x = 1.0 / x;
        return s * sqrt(inv_x / PI) * 
               ((1.0 + inv_x * (C31 + inv_x * (C32 + inv_x * C33))) * cos(x) - 
                (1.0 + inv_x * (D31 + inv_x * (D32 + inv_x * D33))) * sin(x));
    }
}

float bessel_j4(float x) {
    x = abs(x);
    const float A44 = 0.002644492060608329;
    const float A45 = -0.0006004507528516955;
    const float A46 = 0.00003320308950860871;
    const float B41 = -0.20296731043978247;
    const float B42 = 0.04338600070178919;
    const float B43 = -0.0035265908540099847;
    const float B44 = 0.00013712907221840123;
    const float B45 = 0.000012746991211123013;
    const float C41 = 7.881701792737443;
    const float C42 = -27.37611266073206;
    const float C43 = -39.62023054032;
    const float D41 = -7.865445288974054;
    const float D42 = -27.479039704046176;
    const float D43 = 49.286435632834696;
    if (x < 8.0) {
        return (x * x * x * x * (A44 + x * (A45 + x * A46))) / 
               (1.0 + x * (B41 + x * (B42 + x * (B43 + x * (B44 + x * B45)))));
    } else {
        float inv_x = 1.0 / x;
        return sqrt(inv_x / PI) * 
               ((1.0 + inv_x * (C41 + inv_x * (C42 + inv_x * C43))) * cos(x) + 
                (1.0 + inv_x * (D41 + inv_x * (D42 + inv_x * D43))) * sin(x));
    }
}

float bessel_j(float m, float x) {
    if (m == 0.0)
        return bessel_j0(x);
    if (m == 1.0)
        return bessel_j1(x);
    if (m == 2.0)
        return bessel_j2(x);
    if (m == 3.0)
        return bessel_j3(x);
    if (m == 4.0)
        return bessel_j4(x);
    return 0.0;
}


bool get_intersection(
    vec2 pos_2d,
    out vec3 poi, out vec2 uv
) {
    const float radius = 1.0;
    mat3 A = mat3(
        vAxesU, vAxesV,
        vec3(pos_2d, 1.0)
    );
    if (determinant(A) == 0.0) {
        uv = vec2(-radius);
        return false;
    }
    vec3 uvt = -inverse(A) * vPosition;
    uv = uvt.xy;
    if (length(uv) > radius)
        return false;
    float t = -uvt.z;
    poi = vec3(pos_2d*t, t);
    return true;
}


float visibility_kernel(float r2) {
    return r2 > 1.0 ? 0.0 : 1.0-r2;
}

bool get_alpha(
    const vec2 uv,
    const float opac,
    const vec2 aniso,
    out float alpha
) {
    float r2 = dot(uv, uv);
    float vis = visibility_kernel(r2);
    float t = dot(uv, aniso);
    float m = t<0. ? 1. : t>1. ? 0. :
        t*t*(2.0*t-3.0) + 1.0;
    alpha = opac * vis * m;
    return r2 >= 0. && alpha > 1e-4;
}


#define use_ch bool(u_ch_config.x)
#define ch_degree_r u_ch_config.y
#define ch_degree_phi u_ch_config.z
#define ch_dim (ch_degree_r*(2*ch_degree_phi+1))

ivec3 pix_ch;
uvec4 coeff_ch;

void init_coeffs() {
    int width = textureSize(u_ch_texture, 0).x;
    int pixidx = vIndex * ((ch_dim+1)/2);
    pix_ch = ivec3(pixidx % width, pixidx / width, 0);
}

vec3 next_coeff() {
    uvec2 packed;
    if (pix_ch.z == 0) {
        coeff_ch = texelFetch(u_ch_texture, pix_ch.xy, 0);
        packed = coeff_ch.xy;
        pix_ch.z += 1;
    }
    else {
        packed = coeff_ch.zw;
        pix_ch.x += 1;
        pix_ch.z = 0;
    }
    vec2 c1 = unpackHalf2x16(packed.x);
    vec2 c2 = unpackHalf2x16(packed.y);
    return vec3(c1, c2.x);
}

vec3 ch_coeffs_to_color(vec2 uv) {
    if (!use_ch)
        return vec3(0);

    float degree_r = float(ch_degree_r);
    float degree_phi = float(ch_degree_phi);

    float r = length(uv);
    float phi = atan(uv.y, uv.x);
    float pi = 3.14159265358979;

    float idx = 0.0;
    vec3 color = vec3(0.0);

    for (float k = 1.0; k <= degree_r; k++) {
        // m == 0
        float w = bessel_j(0.0, k*pi*r);
        color += w * next_coeff();
        idx += 1.0;
        // m > 0
        for (float m = 1.0; m <= degree_phi; m++) {
            float wr = bessel_j(m, k*pi*r);
            float wc = wr * cos(m*phi);
            float ws = wr * sin(m*phi);
            color += wc * next_coeff();
            idx += 1.0;
            color += ws * next_coeff();
            idx += 1.0;
        }
    }
    return color;
}


// camera distortion models
// TODO: refactor into a separate file

float opencv_radius(
    in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;

    // calculate the extrema of dist(r)
    if (k2 == 0.0) {
        float b = -1.0/(3.0*k1);
        return b > 0.0 ? sqrt(b) : -1.0;
    }
    float disc = 9.0*k1*k1-20.0*k2;
    if (disc <= 0.0) return -1.0;
    disc = sqrt(disc);
    float u1 = (-3.0*k1 + disc) / (10.0*k2);
    float u2 = (-3.0*k1 - disc) / (10.0*k2);
    if (u1 <= 0.0) return u2>0.0 ? sqrt(u2) : -1.0;
    return u2>0.0 ? sqrt(min(u1,u2)) : sqrt(u1);
}

float fisheye_radius(
    in float theta, in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float theta2 = theta*theta;
    return theta*(1.0 + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
}

vec2 undistort_fisheye(in vec2 p, in vec4 undist_coeffs) {
    float l1 = undist_coeffs.x;
    float l2 = undist_coeffs.y;
    float l3 = undist_coeffs.z;
    float l4 = undist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrt(r2);
    float theta = r*(1.0 + r2*(l1 + r2*(l2 + r2*(l3 + r2*l4))));

    return p * (tan(theta)/r);
}


void main () {

    // uvec4 tex = texelFetch(u_ch_texture, ivec2(gl_FragCoord.xy), 0);
    // if (tex.x != 0u) {
    //     vec2 c1 = unpackHalf2x16(tex.x);
    //     vec2 c2 = unpackHalf2x16(tex.y);
    //     fragColor = vec4(c1, c2.x, 1.0);
    //     return;
    // }

    vec2 pos2d = vec2(1,-1)*(gl_FragCoord.xy-0.5*viewport)/focal;
    vec2 pos2d_undist = pos2d;

    float fr = -1.0;
    if (camera_model == 0) {
        fr = opencv_radius(distortion);
    }
    else if (camera_model == 1) {
        fr = fisheye_radius(0.5*PI, distortion);
        // pos2d_undist = undistort_fisheye(pos2d, undistortion);
    }
    const vec3 vignetting_background = vec3(0.1);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 0.0);
        return;
    }

    vec3 poi;
    vec2 uv;
    if (!get_intersection(pos2d_undist, poi, uv))
        discard;
    float alpha;
    if (!get_alpha(uv, vColor.a, vAnisotropy, alpha))
        discard;
    vec3 color = vColor.rgb;
    if (ch_dim > 0) {
        init_coeffs();
        color /= 1.0 + exp(-ch_coeffs_to_color(uv));
    }
    
    if (fr > 0.0) {
        float vignetting = min(length(pos2d)/fr, 1.0);
        vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
        color = mix(vignetting_background, color, vignetting);
    }

    fragColor = vec4(color, alpha);
}
