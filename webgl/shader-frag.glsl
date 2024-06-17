#version 300 es
precision highp float;

uniform vec2 focal;
uniform vec2 viewport;

in vec4 vColor;
in vec3 vPosition;
in vec3 vAxesU;
in vec3 vAxesV;

out vec4 fragColor;


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


void main () {
    vec2 pos_2d = vec2(1,-1)*(gl_FragCoord.xy-0.5*viewport)/focal;
    vec3 poi;
    vec2 uv;
    if (!get_intersection(pos_2d, poi, uv))
        discard;
    float alpha = 1.0-dot(uv,uv);
    vec3 color = vColor.rgb;
    float B = alpha * vColor.a;
    fragColor = vec4(B * color, B);
}
