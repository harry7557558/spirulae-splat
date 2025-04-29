#version 300 es
precision highp float;
precision highp int;


flat in uvec4 vColor0;
flat in uvec4 vColor1;
flat in uvec4 vColor2;

out uvec4 fragColor;


void main () {

    int texel = int(gl_FragCoord.x) % 3;
    fragColor = texel == 0 ? vColor0 : texel == 1 ? vColor1 : vColor2;

    return;

}
