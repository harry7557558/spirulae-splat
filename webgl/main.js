function getProjectionMatrix(fx, fy, width, height) {
    const znear = 0.2;
    const zfar = 200;
    return [
        [(2 * fx) / width, 0, 0, 0],
        [0, -(2 * fy) / height, 0, 0],
        [0, 0, zfar / (zfar - znear), 1],
        [0, 0, -(zfar * znear) / (zfar - znear), 0],
    ].flat();
}

function multiply4(a, b) {
    return [
        b[0] * a[0] + b[1] * a[4] + b[2] * a[8] + b[3] * a[12],
        b[0] * a[1] + b[1] * a[5] + b[2] * a[9] + b[3] * a[13],
        b[0] * a[2] + b[1] * a[6] + b[2] * a[10] + b[3] * a[14],
        b[0] * a[3] + b[1] * a[7] + b[2] * a[11] + b[3] * a[15],
        b[4] * a[0] + b[5] * a[4] + b[6] * a[8] + b[7] * a[12],
        b[4] * a[1] + b[5] * a[5] + b[6] * a[9] + b[7] * a[13],
        b[4] * a[2] + b[5] * a[6] + b[6] * a[10] + b[7] * a[14],
        b[4] * a[3] + b[5] * a[7] + b[6] * a[11] + b[7] * a[15],
        b[8] * a[0] + b[9] * a[4] + b[10] * a[8] + b[11] * a[12],
        b[8] * a[1] + b[9] * a[5] + b[10] * a[9] + b[11] * a[13],
        b[8] * a[2] + b[9] * a[6] + b[10] * a[10] + b[11] * a[14],
        b[8] * a[3] + b[9] * a[7] + b[10] * a[11] + b[11] * a[15],
        b[12] * a[0] + b[13] * a[4] + b[14] * a[8] + b[15] * a[12],
        b[12] * a[1] + b[13] * a[5] + b[14] * a[9] + b[15] * a[13],
        b[12] * a[2] + b[13] * a[6] + b[14] * a[10] + b[15] * a[14],
        b[12] * a[3] + b[13] * a[7] + b[14] * a[11] + b[15] * a[15],
    ];
}

function invert4(a) {
    let b00 = a[0] * a[5] - a[1] * a[4];
    let b01 = a[0] * a[6] - a[2] * a[4];
    let b02 = a[0] * a[7] - a[3] * a[4];
    let b03 = a[1] * a[6] - a[2] * a[5];
    let b04 = a[1] * a[7] - a[3] * a[5];
    let b05 = a[2] * a[7] - a[3] * a[6];
    let b06 = a[8] * a[13] - a[9] * a[12];
    let b07 = a[8] * a[14] - a[10] * a[12];
    let b08 = a[8] * a[15] - a[11] * a[12];
    let b09 = a[9] * a[14] - a[10] * a[13];
    let b10 = a[9] * a[15] - a[11] * a[13];
    let b11 = a[10] * a[15] - a[11] * a[14];
    let det =
        b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
    if (!det) return null;
    return [
        (a[5] * b11 - a[6] * b10 + a[7] * b09) / det,
        (a[2] * b10 - a[1] * b11 - a[3] * b09) / det,
        (a[13] * b05 - a[14] * b04 + a[15] * b03) / det,
        (a[10] * b04 - a[9] * b05 - a[11] * b03) / det,
        (a[6] * b08 - a[4] * b11 - a[7] * b07) / det,
        (a[0] * b11 - a[2] * b08 + a[3] * b07) / det,
        (a[14] * b02 - a[12] * b05 - a[15] * b01) / det,
        (a[8] * b05 - a[10] * b02 + a[11] * b01) / det,
        (a[4] * b10 - a[5] * b08 + a[7] * b06) / det,
        (a[1] * b08 - a[0] * b10 - a[3] * b06) / det,
        (a[12] * b04 - a[13] * b02 + a[15] * b00) / det,
        (a[9] * b02 - a[8] * b04 - a[11] * b00) / det,
        (a[5] * b07 - a[4] * b09 - a[6] * b06) / det,
        (a[0] * b09 - a[1] * b07 + a[2] * b06) / det,
        (a[13] * b01 - a[12] * b03 - a[14] * b00) / det,
        (a[8] * b03 - a[9] * b01 + a[10] * b00) / det,
    ];
}

function rotate4(a, rad, x, y, z, relative=false) {
    // let dist = 0.1;
    let inv = invert4(a);
    let dist = -1.0 * Math.sqrt(0.5 * (
        a[12]*a[12] + a[13]*a[13] + a[14]*a[14] + a[15]*a[15]
    ));
    let tr = [-inv[12]+dist*inv[8], -inv[13]+dist*inv[9], -inv[14]+dist*inv[10]];
    if (relative) {
        a = translate4(a, -tr[0], -tr[1], -tr[2], false);
    }
    let len = Math.hypot(x, y, z);
    x /= len;
    y /= len;
    z /= len;
    let s = Math.sin(rad);
    let c = Math.cos(rad);
    let t = 1 - c;
    let b00 = x * x * t + c;
    let b01 = y * x * t + z * s;
    let b02 = z * x * t - y * s;
    let b10 = x * y * t - z * s;
    let b11 = y * y * t + c;
    let b12 = z * y * t + x * s;
    let b20 = x * z * t + y * s;
    let b21 = y * z * t - x * s;
    let b22 = z * z * t + c;
    a = [
        a[0] * b00 + a[4] * b01 + a[8] * b02,
        a[1] * b00 + a[5] * b01 + a[9] * b02,
        a[2] * b00 + a[6] * b01 + a[10] * b02,
        a[3] * b00 + a[7] * b01 + a[11] * b02,
        a[0] * b10 + a[4] * b11 + a[8] * b12,
        a[1] * b10 + a[5] * b11 + a[9] * b12,
        a[2] * b10 + a[6] * b11 + a[10] * b12,
        a[3] * b10 + a[7] * b11 + a[11] * b12,
        a[0] * b20 + a[4] * b21 + a[8] * b22,
        a[1] * b20 + a[5] * b21 + a[9] * b22,
        a[2] * b20 + a[6] * b21 + a[10] * b22,
        a[3] * b20 + a[7] * b21 + a[11] * b22,
        ...a.slice(12, 16),
    ];
    if (relative) {
        a = translate4(a, tr[0], tr[1], tr[2], false);
    }
    return a;
}

function translate4(a, x, y, z, relative=true) {
    if (relative) {
        let sc = Math.sqrt(0.5 * (
            a[12]*a[12] + a[13]*a[13] + a[14]*a[14] + a[15]*a[15]
        ));
        x *= sc, y *= sc, z *= sc;
    }
    return [
        ...a.slice(0, 12),
        a[0] * x + a[4] * y + a[8] * z + a[12],
        a[1] * x + a[5] * y + a[9] * z + a[13],
        a[2] * x + a[6] * y + a[10] * z + a[14],
        a[3] * x + a[7] * y + a[11] * z + a[15],
    ];
}


let Module;

// config: JavaScript object
// packedData: ArrayBuffer
// buffers: Array of ArrayBuffers
function unpackComponents(config, packedData, buffers) {
    const componentLength = config.componentLength;
    const componentViews = config.componentViews;

    if (!packedData)
        return;
    const length = Math.min(config.length, packedData.byteLength/componentLength);
    
    const components = componentViews.map(view => {
        const dataType = view.type;
        if (dataType.startsWith('quat')) {
            const numElements = Math.max(Number(dataType.slice(4)), 1);
            return new Uint32Array(length * numElements);
        } else {
            // TODO: assert only quat is supported
        }
    });

    let packedDataByte = new Uint8Array(packedData);

    for (let i = 0; i < length; i++) {
        componentViews.forEach((view, j) => {
            const bitLength = view.bitLength;
            const bitOffset = view.bitOffset;
            const dataType = view.type;
            
            if (dataType.startsWith('quat')) {
                const numElements = Math.max(Number(dataType.slice(4)), 1);
                const bitsPerElement = bitLength / numElements;

                for (let k = 0; k < numElements; k++) {
                    let value = 0;
                    const bitStart = bitOffset + k * bitsPerElement;
                    const byteStart = Math.floor(bitStart / 8);
                    const bitPos = bitStart % 8;

                    for (let b = 0; b < bitsPerElement; b++) {
                        const byteIdx = byteStart + Math.floor((bitPos + b) / 8);
                        const bitInByte = (bitPos + b) % 8;
                        if (packedDataByte[i * componentLength + byteIdx] & (1 << bitInByte)) {
                            value |= (1 << b);
                        }
                    }
                    components[j][i * numElements + k] = value;
                }

            } else {
                // TODO: assert only quat is supported
            }
        });
    }

    const object = {};
    componentViews.forEach((view, j) => {
        let hasQuat = view.type.startsWith('quat');
        if (hasQuat) {
            let bl = buffers[view.quatBufferView].byteLength;
            if (bl < 4 || (bl & (bl - 1)) != 0)
                hasQuat = false;
        }
        if (hasQuat) {
            let numel = components[j].length;
            let component = new Float32Array(numel);
            let buffer = new Float32Array(buffers[view.quatBufferView]);
            for (var i = 0; i < numel; i++)
                component[i] = buffer[components[j][i]];
            object[view.key] = component;
        }
        else {
            object[view.key] = components[j];
        }
    });
    return object;
}

function unpackModel(header, buffer) {
    let model = {
        header: header
    };

    let buffers = [];
    header.bufferViews.forEach((view, j) => {
        buffers.push(buffer.slice(
            view.byteOffset, view.byteOffset+view.byteLength));
    });

    console.time("unpack model");
    for (var key in header.primitives) {
        let primitive = header.primitives[key];
        let bufferSlice = buffers[primitive.bufferView];
        // let components = unpackComponents(primitive, bufferSlice, buffers);
        let components = Module.unpackComponents(key, primitive, bufferSlice, buffers);
        model[key] = components;
    };
    console.timeEnd("unpack model");

    model.header = header;
    return model;
}


function createWorker(self) {
    let header;
    let base, harmonics;
    let vertexCount = 0;
    let viewProj;

    let lastProj = [];
    let depthIndex = new Uint32Array();
    let lastBaseVertexCount = 0,
        lastHarmonicsLength = 0;

    var _floatView = new Float32Array(1);
    var _int32View = new Int32Array(_floatView.buffer);

    function floatToHalf(float) {
        _floatView[0] = float;
        var f = _int32View[0];

        var sign = (f >> 31) & 0x0001;
        var exp = (f >> 23) & 0x00ff;
        var frac = f & 0x007fffff;

        var newExp;
        if (exp == 0) {
            newExp = 0;
        } else if (exp < 113) {
            newExp = 0;
            frac |= 0x00800000;
            frac = frac >> (113 - exp);
            if (frac & 0x01000000) {
                newExp = 1;
                frac = 0;
            }
        } else if (exp < 142) {
            newExp = exp - 112;
        } else {
            newExp = 31;
            frac = 0;
        }

        return (sign << 15) | (newExp << 10) | (frac >> 13);
    }

    function packHalf2x16(x, y) {
        return (floatToHalf(x) | (floatToHalf(y) << 16)) >>> 0;
    }

    function generateBaseTexture() {
        if (!base) return;
        let vertexCount = base.opacities.length;

        var basewidth = 1024 * 3;
        var baseheight = Math.ceil((3 * vertexCount) / basewidth);
        var basedata = new Uint32Array(basewidth * baseheight * 4);
        var basedataF = new Float32Array(basedata.buffer);
        for (let i = 0; i < vertexCount; i++) {
            let iTex = 12 * i;

            // position - 3 floats, info0 xyz
            basedataF[iTex + 0] = base.means[3*i+0];
            basedataF[iTex + 1] = base.means[3*i+1];
            basedataF[iTex + 2] = base.means[3*i+2];

            // scale - 2 floats, info0 w, info1 x
            basedataF[iTex + 3] = Math.exp(base.scales[2*i+0]);
            basedataF[iTex + 4] = Math.exp(base.scales[2*i+1]);

            // quat - 4 halfs, info1 yz
            let quat = [
                base.quats[4*i+0],
                base.quats[4*i+1],
                base.quats[4*i+2],
                base.quats[4*i+3],
            ];
            basedata[iTex + 5] = packHalf2x16(quat[0], quat[1]);
            basedata[iTex + 6] = packHalf2x16(quat[2], quat[3]);

            // anisotropy - 2 floats, info1 w, info2 x
            if (header.config.use_anisotropy !== false && base.anisotropies) {
                basedataF[iTex + 7] = Math.sinh(base.anisotropies[2*i+0]);
                basedataF[iTex + 8] = Math.sinh(base.anisotropies[2*i+1]);
            }
            else {
                basedataF[iTex + 7] = basedataF[iTex + 8] = 0.0;
            }

            // rgba - 4 halfs, info2 yz
            let m = header.config.sh_degree == 0 ? 0.28209479177387814 : 1.0;
            let b = header.config.sh_degree == 0 ? 0.5 : 0.0;
            let rgba = [
                b + m * base.features_dc[3*i+0],
                b + m * base.features_dc[3*i+1],
                b + m * base.features_dc[3*i+2],
                base.opacities[i]
            ];
            basedata[iTex + 9] = packHalf2x16(rgba[0], rgba[1]);
            basedata[iTex + 10] = packHalf2x16(rgba[2], rgba[3]);
        }

        self.postMessage({
            baseTexture: { basedata, basewidth, baseheight }
        }, [basedata.buffer]);
    }

    function generateCHTexture() {
        if (!harmonics) return;
        let degree_r = header.config.ch_degree_r;
        let degree_phi = header.config.ch_degree_phi;
        let dim_ch = degree_r * (2*degree_phi+1);
        if (dim_ch == 0)
            return;
        let pixelPerVert = Math.ceil(dim_ch/2);
        let vertexCount = Math.floor(harmonics.features_ch.length/(3*dim_ch));

        var chwidth = Math.floor(4096/pixelPerVert)*pixelPerVert;
        var chheight = Math.ceil(pixelPerVert*vertexCount/chwidth);
        var chdata = Module.packHarmonicTexture(
            "ch", dim_ch, chwidth, chheight, harmonics.features_ch);

        self.postMessage({
            chTexture: { chdata, chwidth, chheight }
        }, [chdata.buffer]);
    }

    function generateSHTexture() {
        if (!harmonics) return;
        let degree = header.config.sh_degree;
        let dim_sh = degree * (degree + 2);
        let pixelPerVert = Math.ceil(dim_sh/2);
        let vertexCount = Math.floor(harmonics.features_sh.length/(3*dim_sh));

        var shwidth = Math.floor(4096/pixelPerVert)*pixelPerVert;
        var shheight = Math.ceil(pixelPerVert*vertexCount/shwidth);
        var shdata = Module.packHarmonicTexture(
            "sh", dim_sh, shwidth, shheight, harmonics.features_sh);

        self.postMessage({
            shTexture: { shdata, shwidth, shheight }
        }, [shdata.buffer]);
    }

    function runSort(viewProj) {
        if (!base) return;
        if (harmonics && harmonics.features_sh.length > lastHarmonicsLength) {
            console.time("generate CH/SH texture");
            generateCHTexture();
            generateSHTexture();
            console.timeEnd("generate CH/SH texture");
            lastHarmonicsLength = harmonics.features_sh.length;
        }
        let vertexCount = base.opacities.length;
        if (lastBaseVertexCount == vertexCount) {
            let dot =
                lastProj[2] * viewProj[2] +
                lastProj[6] * viewProj[6] +
                lastProj[10] * viewProj[10];
            if (Math.abs(dot - 1) < 0.01) {
                return;
            }
        } else {
            console.time("generate base texture");
            generateBaseTexture();
            lastBaseVertexCount = vertexCount;
            console.timeEnd("generate base texture");
        }

        console.time("sort");
        let maxDepth = -Infinity;
        let minDepth = Infinity;
        let sizeList = new Int32Array(vertexCount);
        for (let i = 0; i < vertexCount; i++) {
            let depth =
                -(viewProj[2] * base.means[3*i+0] +
                    viewProj[6] * base.means[3*i+1] +
                    viewProj[10] * base.means[3*i+2]);
            depth = (depth * 4096) | 0;
            sizeList[i] = depth;
            if (depth > maxDepth) maxDepth = depth;
            if (depth < minDepth) minDepth = depth;
        }

        // This is a 16 bit single-pass counting sort
        let depthInv = (65536) / (maxDepth - minDepth);
        let counts0 = new Uint32Array(65536);
        for (let i = 0; i < vertexCount; i++) {
            sizeList[i] = ((sizeList[i] - minDepth) * depthInv) | 0;
            counts0[sizeList[i]]++;
        }
        let starts0 = new Uint32Array(65536);
        for (let i = 1; i < 65536; i++)
            starts0[i] = starts0[i - 1] + counts0[i - 1];
        depthIndex = new Uint32Array(vertexCount).fill(-1);
        for (let i = 0; i < vertexCount; i++)
            depthIndex[starts0[sizeList[i]]++] = i;

        console.timeEnd("sort");

        lastProj = viewProj;
        self.postMessage({ depthIndex, viewProj, vertexCount }, [
            depthIndex.buffer,
        ]);
    }

    const throttledSort = () => {
        if (!sortRunning) {
            sortRunning = true;
            let lastView = viewProj;
            runSort(lastView);
            setTimeout(() => {
                sortRunning = false;
                if (lastView !== viewProj) {
                    throttledSort();
                }
            }, 0);
        }
    };

    let sortRunning;
    window.addEventListener("message", (e) => {
        if (e.data.header) {
            base = e.data.base;
            harmonics = e.data.harmonics;
            header = e.data.header;
            postMessage({ base, harmonics });
        } else if (e.data.view) {
            viewProj = e.data.view;
            throttledSort();
        }
    });
}

let genericVertexShaderSource = `#version 300 es
precision highp float;
layout(location = 0) in vec2 vertexPosition;
void main() { gl_Position = vec4(vertexPosition, 0.0, 1.0); }
`
let vertexShaderSource = "";
let fragmentShaderSource = "";
let backgroundShaderSource = "";

let defaultViewMatrix = [
    1.0, 0.0, 0.0, 0.0,
    0.0, -0.199, 0.98, 0.0,
    0.0, -0.98, -0.199, 0.0,
    0.0, 0.0, 1.0, 1.0
];
// defaultViewMatrix = [-0.99396310067939, -0.0028527761475104297, 0.10967788056938267, 0, -0.10947236494874289, -0.04064304253875142, -0.9931590623621941, 0, 0.007290894500575474, -0.9991701580618988, 0.040085384110577396, 0, -0.008981957487634048, -0.019741792911632222, 0.80080398369808, 0.9999999999999599];
// defaultViewMatrix = [0.33168228758220725, 0.1352309888257432, 0.9336489201265321, 0, -0.9433629759026162, 0.039895202904019586, 0.3293551310485557, 0, 0.007290894500575899, -0.9900110870320623, 0.1408044397597277, 0, -0.13138587004540464, -0.07371462201299227, 0.49180482237051854, 0.9999999999998919];
let viewMatrix = defaultViewMatrix;
async function main() {
    let carousel = true;
    const params = new URLSearchParams(location.search);
    try {
        viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
        carousel = false;
    } catch (err) {}
    const url = new URL(
        params.get("url") || "model.ssplat",
        location.href
    );
    console.log(url);
    const req = await fetch(url, {
        mode: "cors", // no-cors, *cors, same-origin
        credentials: "omit", // include, *same-origin, omit
    });
    console.log(req);
    if (req.status != 200)
        throw new Error(req.status + " Unable to load " + req.url);

    const reader = req.body.getReader();
    let ssplatData = new Uint8Array(req.headers.get("content-length"));

    const downsample = 1;

    // const worker = new Worker(
    //     URL.createObjectURL(
    //         new Blob(["(", createWorker.toString(), ")(self)"], {
    //             type: "application/javascript",
    //         }),
    //     ),
    // );
    createWorker(window);
    let worker = window;

    const canvas = document.getElementById("canvas");
    const fps = document.getElementById("fps");
    let renderNeeded = true;

    let projectionMatrix;

    const gl = canvas.getContext("webgl2", {
        antialias: false,
    });

    // create shader programs
    function createShaderProgram(vsSource, fsSource) {
        const vertexShader = gl.createShader(gl.VERTEX_SHADER);
        gl.shaderSource(vertexShader, vsSource);
        gl.compileShader(vertexShader);
        if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS))
            console.error(gl.getShaderInfoLog(vertexShader));

        const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
        gl.shaderSource(fragmentShader, fsSource);
        gl.compileShader(fragmentShader);
        if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS))
            console.error(gl.getShaderInfoLog(fragmentShader));

        const program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        gl.useProgram(program);
        if (!gl.getProgramParameter(program, gl.LINK_STATUS))
            console.error(gl.getProgramInfoLog(program));

        return program;
    }
    const splatProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);
    const backgroundProgram = createShaderProgram(genericVertexShaderSource, backgroundShaderSource);

    // depth test and blend
    gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.blendEquation(gl.FUNC_ADD);

    // positions
    const triangleVertices = new Float32Array([-1, -1, 1, -1, 1, 1, -1, 1]);
    const vertexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);

    gl.useProgram(splatProgram);
    const a_position = gl.getAttribLocation(splatProgram, "vertexPosition");
    gl.enableVertexAttribArray(a_position);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(a_position, 2, gl.FLOAT, false, 0, 0);

    gl.useProgram(backgroundProgram);
    const bg_position = gl.getAttribLocation(backgroundProgram, "vertexPosition");
    gl.enableVertexAttribArray(bg_position);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(bg_position, 2, gl.FLOAT, false, 0, 0);

    // textures
    gl.useProgram(splatProgram);

    var baseTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, baseTexture);
    var uBaseTextureLocation = gl.getUniformLocation(splatProgram, "u_base_texture");
    gl.uniform1i(uBaseTextureLocation, 0);

    var chTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, chTexture);
    var uChTextureLocation = gl.getUniformLocation(splatProgram, "u_ch_texture");
    gl.uniform1i(uChTextureLocation, 1);

    var shTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, shTexture);
    var uShTextureLocation = gl.getUniformLocation(splatProgram, "u_sh_texture");
    gl.uniform1i(uShTextureLocation, 2);

    // indices
    const indexBuffer = gl.createBuffer();
    const a_index = gl.getAttribLocation(splatProgram, "index");
    gl.enableVertexAttribArray(a_index);
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
    gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
    gl.vertexAttribDivisor(a_index, 1);

    function setUniforms(program=null) {
        if (program === null) {
            setUniforms(splatProgram);
            setUniforms(backgroundProgram);
            return;
        }
        let camera = CameraPresets.camera;
        gl.useProgram(program);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "projection"),
            false, projectionMatrix);
        gl.uniform2fv(gl.getUniformLocation(program, "focal"),
            new Float32Array([camera.fx, camera.fy]));
        gl.uniform2fv(gl.getUniformLocation(program, "viewport"),
            new Float32Array([innerWidth, innerHeight]));
        gl.uniform1i(gl.getUniformLocation(program, "camera_model"), camera.model);
        gl.uniform4fv(gl.getUniformLocation(program, "distortion"),
            new Float32Array(camera.dist_coeffs));
        gl.uniform4fv(gl.getUniformLocation(program, "undistortion"),
            new Float32Array(camera.undist_coeffs));
    }

    const resize = () => {
        let camera = CameraPresets.camera;

        projectionMatrix = getProjectionMatrix(camera.fx, camera.fy, innerWidth, innerHeight);

        gl.canvas.width = Math.round(innerWidth / downsample);
        gl.canvas.height = Math.round(innerHeight / downsample);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

        setUniforms();
        renderNeeded = true;
    };

    window.addEventListener("resize", function() {
        setTimeout(resize, 1);
    });
    resize();

    window.addEventListener("message", (e) => {
        if (e.data.buffer) {
            ssplatData = new Uint8Array(e.data.buffer);
            const blob = new Blob([ssplatData.buffer], {
                type: "application/octet-stream",
            });
            const link = document.createElement("a");
            link.download = "model.ssplat";
            link.href = URL.createObjectURL(blob);
            document.body.appendChild(link);
            link.click();
        } else if (e.data.baseTexture) {
            const { basedata, basewidth, baseheight } = e.data.baseTexture;
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, baseTexture);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
            gl.texImage2D(
                gl.TEXTURE_2D, 0, gl.RGBA32UI,
                basewidth, baseheight, 0,
                gl.RGBA_INTEGER, gl.UNSIGNED_INT, basedata);
        } else if (e.data.chTexture) {
            const { chdata, chwidth, chheight } = e.data.chTexture;
            gl.activeTexture(gl.TEXTURE1);
            gl.bindTexture(gl.TEXTURE_2D, chTexture);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
            gl.texImage2D(
                gl.TEXTURE_2D, 0, gl.RGBA32UI,
                chwidth, chheight, 0,
                gl.RGBA_INTEGER, gl.UNSIGNED_INT, chdata);
        } else if (e.data.shTexture) {
            const { shdata, shwidth, shheight } = e.data.shTexture;
            gl.activeTexture(gl.TEXTURE2);
            gl.bindTexture(gl.TEXTURE_2D, shTexture);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
            gl.texImage2D(
                gl.TEXTURE_2D, 0, gl.RGBA32UI,
                shwidth, shheight, 0,
                gl.RGBA_INTEGER, gl.UNSIGNED_INT, shdata);
        } else if (e.data.depthIndex) {
            const { depthIndex, viewProj } = e.data;
            gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);
            vertexCount = e.data.vertexCount;
        }
    });

    let activeKeys = [];

    let lastFrame = 0;
    let avgFps = 0;
    let start = 0;
    function getVelocityScale() {
        return 120.0 / Math.min(Math.max(avgFps, 30.0), 240.0);
    }

    window.addEventListener("keydown", (e) => {
        // if (document.activeElement != document.body) return;
        if (/^Key[A-Z]$/.test(e.code) || /^Arrow/.test(e.code) || e.code == "Home" || e.code == "End")
            e.preventDefault();
        if (!activeKeys.includes(e.code)) activeKeys.push(e.code);
        renderNeeded = true;
    });
    window.addEventListener("keyup", (e) => {
        activeKeys = activeKeys.filter((k) => k !== e.code);
        renderNeeded = true;
    });
    window.addEventListener("blur", () => {
        activeKeys = [];
        renderNeeded = true;
    });

    window.addEventListener(
        "wheel",
        (e) => {
            carousel = false;
            e.preventDefault();
            const lineHeight = 10;
            const scale =
                e.deltaMode == 1
                    ? lineHeight
                    : e.deltaMode == 2
                    ? innerHeight
                    : 1;
            // let innerSize = Math.max(innerWidth, innerHeight);
            let innerSize = 840;
            let inv = invert4(viewMatrix);
            if (e.shiftKey) {
                inv = translate4(
                    inv,
                    (e.deltaX * scale) / innerSize,
                    (e.deltaY * scale) / innerSize,
                    0,
                );
            } else {
                let d = -(e.ctrlKey || e.metaKey ?  4.0 : 0.5)
                    * scale * e.deltaY / innerSize;
                inv = translate4(inv, 0, 0, d);
            }

            viewMatrix = invert4(inv);
            renderNeeded = true;
        },
        { passive: false },
    );

    let startX, startY, down;
    canvas.addEventListener("mousedown", (e) => {
        carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = e.ctrlKey || e.metaKey ? 2 : 1;
        renderNeeded = true;
    });
    canvas.addEventListener("contextmenu", (e) => {
        carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = 2;
        renderNeeded = true;
    });

    canvas.addEventListener("mousemove", (e) => {
        e.preventDefault();
        if (!down) return;
        if (e.ctrlKey || down == 2) {
            let inv = invert4(viewMatrix);
            inv = translate4(
                inv,
                (-5 * (e.clientX - startX)) / innerWidth,
                0,
                (5 * (e.clientY - startY)) / innerHeight,
            );
            viewMatrix = invert4(inv);
        }
        else if (e.shiftKey) {
            let inv = invert4(viewMatrix);
            inv = translate4(
                inv,
                (-5 * (e.clientX - startX)) / innerWidth,
                (-5 * (e.clientY - startY)) / innerHeight,
                0,
            );
            viewMatrix = invert4(inv);
        }
        else if (down == 1) {
            let dx = (5 * (e.clientX - startX)) / innerWidth;
            let dy = (5 * (e.clientY - startY)) / innerHeight;
            viewMatrix = rotate4(viewMatrix, dx, 0, 0, 1, true);
            viewMatrix = rotate4(viewMatrix, dy, viewMatrix[0], viewMatrix[4], viewMatrix[8], true);
        }
        startX = e.clientX;
        startY = e.clientY;
        renderNeeded = true;
    });
    canvas.addEventListener("mouseup", (e) => {
        e.preventDefault();
        down = false;
        startX = 0;
        startY = 0;
        renderNeeded = true;
    });

    let altX = 0,
        altY = 0;
    canvas.addEventListener(
        "touchstart",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1) {
                carousel = false;
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
                down = 1;
            } else if (e.touches.length === 2) {
                // console.log('beep')
                carousel = false;
                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
                down = 1;
            }
            renderNeeded = true;
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchmove",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1 && down) {
                let dx = (2 * (e.touches[0].clientX - startX)) / innerWidth;
                let dy = (2 * (e.touches[0].clientY - startY)) / innerHeight;
    
                viewMatrix = rotate4(viewMatrix, dx, 0, 0, 1, true);
                viewMatrix = rotate4(viewMatrix, dy, viewMatrix[0], viewMatrix[4], viewMatrix[8], true);
    
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
            } else if (e.touches.length === 2) {
                // alert('beep')
                const normalizeAngle = angle => (
                    (angle % (2 * Math.PI)) + 3 * Math.PI) % (2 * Math.PI) - Math.PI;
                const atan2_scale = (dx, dy) => {
                    const dtheta_r = 10.0;
                    return Math.tanh(1.5 *
                        (Math.sqrt(dx*dx + dy*dy + dtheta_r*dtheta_r) - dtheta_r) /
                        Math.sqrt(innerWidth * innerHeight)
                    );
                };
                let dy = startY - altY, dx = startX - altX;
                let dtx = e.touches[0].clientX - e.touches[1].clientX,
                    dty = e.touches[0].clientY - e.touches[1].clientY;
                const dtheta = normalizeAngle(Math.atan2(dy, dx) - Math.atan2(dty, dtx))
                    * Math.sqrt(atan2_scale(dx, dy) * atan2_scale(dtx, dty));
                const dscale_r = 5.0;
                const dscale = (Math.hypot(dx, dy) + dscale_r) /
                    (Math.hypot(dtx, dty) + dscale_r);
                dx =
                    2.0 * (e.touches[0].clientX +
                        e.touches[1].clientX -
                        (startX + altX)) /
                    2;
                dy =
                    2.0 * (e.touches[0].clientY +
                        e.touches[1].clientY -
                        (startY + altY)) /
                    2;
                let inv = invert4(viewMatrix);
                // inv = translate4(inv,  0, 0, d);
                inv = rotate4(inv, dtheta, 0, 0, 1, true);

                inv = translate4(inv, -dx / innerWidth, -dy / innerHeight, 0);

                // let preY = inv[13];
                inv = translate4(inv, 0, 0, 0.7 * (1 - dscale));
                // inv[13] = preY;

                viewMatrix = invert4(inv);

                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
            }
            renderNeeded = true;
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchend",
        (e) => {
            e.preventDefault();
            down = false;
            startX = 0;
            startY = 0;
            renderNeeded = true;
        },
        { passive: false },
    );

    let jumpDelta = 0;
    let vertexCount = 0;

    window.addEventListener("gamepadconnected", (e) => {
        const gp = navigator.getGamepads()[e.gamepad.index];
        console.log(
            `Gamepad connected at index ${gp.index}: ${gp.id}. It has ${gp.buttons.length} buttons and ${gp.axes.length} axes.`,
        );
        renderNeeded = true;
    });
    window.addEventListener("gamepaddisconnected", (e) => {
        console.log("Gamepad disconnected");
        renderNeeded = true;
    });

    let bytesRead = 0;
    let previousCamera = "";
    const frame = (now) => {
        let inv = invert4(viewMatrix);
        let shiftKey = activeKeys.includes("Shift") || activeKeys.includes("ShiftLeft") || activeKeys.includes("ShiftRight");
        let vsc = getVelocityScale();

        if (/^Key[ADQEWSJKLI]$/.test(activeKeys) ||
            /^Arrow/.test(activeKeys)) renderNeeded = true;

        if (activeKeys.includes("ArrowUp"))
            inv = translate4(inv, 0, -0.01*vsc, 0);
        if (activeKeys.includes("ArrowDown"))
            inv = translate4(inv, 0, 0.01*vsc, 0);
        if (activeKeys.includes("ArrowLeft"))
            inv = translate4(inv, -0.01*vsc, 0, 0);
        if (activeKeys.includes("ArrowRight"))
            inv = translate4(inv, 0.01*vsc, 0, 0);
        if (activeKeys.includes("KeyA")) inv = rotate4(inv, -0.01*vsc, 0, 1, 0);
        if (activeKeys.includes("KeyD")) inv = rotate4(inv, 0.01*vsc, 0, 1, 0);
        if (activeKeys.includes("KeyQ")) inv = rotate4(inv, 0.01*vsc, 0, 0, 1);
        if (activeKeys.includes("KeyE")) inv = rotate4(inv, -0.01*vsc, 0, 0, 1);
        if (activeKeys.includes("KeyW")) inv = rotate4(inv, 0.01*vsc, 1, 0, 0);
        if (activeKeys.includes("KeyS")) inv = rotate4(inv, -0.01*vsc, 1, 0, 0);

        const gamepads = navigator.getGamepads ? navigator.getGamepads() : [];
        let isJumping = activeKeys.includes("Space");
        for (let gamepad of gamepads) {
            if (!gamepad) continue;

            const axisThreshold = 0.1; // Threshold to detect when the axis is intentionally moved
            const moveSpeed = 0.02*vsc;
            const rotateSpeed = 0.02*vsc;

            let inv0 = JSON.stringify(inv);

            // Assuming the left stick controls translation (axes 0 and 1)
            if (Math.abs(gamepad.axes[0]) > axisThreshold) {
                inv = translate4(inv, moveSpeed * gamepad.axes[0], 0, 0);
                carousel = false;
            }
            if (Math.abs(gamepad.axes[1]) > axisThreshold) {
                inv = translate4(inv, 0, 0, -moveSpeed * gamepad.axes[1]);
                carousel = false;
            }
            if(gamepad.buttons[12].pressed || gamepad.buttons[13].pressed){
                inv = translate4(inv, 0, -moveSpeed*(gamepad.buttons[12].pressed - gamepad.buttons[13].pressed), 0);
                carousel = false;
            }
            if(gamepad.buttons[14].pressed || gamepad.buttons[15].pressed){
                inv = translate4(inv, -moveSpeed*(gamepad.buttons[14].pressed - gamepad.buttons[15].pressed), 0, 0);
                carousel = false;
            }

            // Assuming the right stick controls rotation (axes 2 and 3)
            if (Math.abs(gamepad.axes[2]) > axisThreshold) {
                inv = rotate4(inv, rotateSpeed * gamepad.axes[2], 0, 1, 0);
                carousel = false;
            }
            if (Math.abs(gamepad.axes[3]) > axisThreshold) {
                inv = rotate4(inv, -rotateSpeed * gamepad.axes[3], 1, 0, 0);
                carousel = false;
            }

            let tiltAxis = gamepad.buttons[6].value - gamepad.buttons[7].value;
            if (Math.abs(tiltAxis) > axisThreshold) {
                inv = rotate4(inv, rotateSpeed * tiltAxis, 0, 0, 1);
                carousel = false;
            }
            if (gamepad.buttons[0].pressed) {
                isJumping = true;
                carousel = false;
            }
            if(gamepad.buttons[3].pressed){
                carousel = true;
            }

            if (JSON.stringify(inv) != inv0)
                renderNeeded = true;
        }

        if (
            ["KeyJ", "KeyK", "KeyL", "KeyI"].some((k) => activeKeys.includes(k))
        ) {
            let d = Math.sqrt(0.5 * (
                inv[12]*inv[12] + inv[13]*inv[13] + inv[14]*inv[14] + inv[15]*inv[15]
            ));
            inv = translate4(inv, 0, 0, d, false);
            inv = rotate4(
                inv,
                activeKeys.includes("KeyJ")
                    ? -0.02*vsc
                    : activeKeys.includes("KeyL")
                    ? 0.02*vsc
                    : 0,
                0, 1, 0
            );
            inv = rotate4(
                inv,
                activeKeys.includes("KeyI")
                    ? 0.02*vsc
                    : activeKeys.includes("KeyK")
                    ? -0.02*vsc
                    : 0,
                1, 0, 0
            );
            inv = translate4(inv, 0, 0, -d, false);
        }

        viewMatrix = invert4(inv);

        if (carousel) {
            let inv = invert4(defaultViewMatrix);

            const t = Math.sin((Date.now() - start) / 5000);
            inv = translate4(inv, 0.2 * t, 0, 0.1 * (1 - Math.cos(t)));
            inv = rotate4(inv, -0.1 * t, 0, 1, 0);

            viewMatrix = invert4(inv);
            renderNeeded = true;
        }

        let jumpDeltaNew = jumpDelta;
        if (isJumping) {
            jumpDeltaNew = Math.min(1, jumpDelta + 0.05*vsc);
        } else {
            jumpDeltaNew = Math.max(0, jumpDelta - 0.05*vsc);
        }
        if (jumpDeltaNew != jumpDelta) {
            jumpDelta = jumpDeltaNew;
            renderNeeded = true;
        }

        let inv2 = invert4(viewMatrix);
        inv2 = translate4(inv2, 0, -jumpDelta, 0);
        inv2 = rotate4(inv2, -0.1 * jumpDelta, 1, 0, 0);
        let actualViewMatrix = invert4(inv2);

        const viewProj = multiply4(projectionMatrix, actualViewMatrix);
        worker.postMessage({ view: viewProj });

        const currentFps = 1000 / (now - lastFrame) || 0;
        avgFps = avgFps * 0.8 + currentFps * 0.2;

        let camera = CameraPresets.camera;
        if (camera !== previousCamera)
            renderNeeded = true;
        previousCamera = camera;

        // renderNeeded = true;
        if (vertexCount > 0 && renderNeeded) {
            document.getElementById("spinner").style.display = "none";

            let use_anisotropy = header.config.use_anisotropy !== false;
            document.getElementById("checkbox-aniso").disabled = !use_anisotropy;
            let sh_degree = header.config.sh_degree;
            document.getElementById("checkbox-sh").disabled = (sh_degree == 0);
            let ch_degree_r = header.config.ch_degree_r;
            document.getElementById("checkbox-ch").disabled = (ch_degree_r == 0);

            setUniforms();

            let background = header.config.background_color;
            document.getElementById("checkbox-bg").disabled = !(Math.max(background[0], background[1], background[2]) > 0);
            let bg = document.getElementById("checkbox-bg").checked;
            gl.clearColor(bg*background[0], bg*background[1], bg*background[2], 1.0);
            gl.clear(gl.COLOR_BUFFER_BIT);

            let background_sh_degree = header.config.background_sh_degree;
            if ((bg && background_sh_degree > 0) || camera.model == 1) {
                gl.useProgram(backgroundProgram);
                gl.uniform1f(gl.getUniformLocation(backgroundProgram, "sh_degree"),
                    background_sh_degree);
                let coeff = header.config.background_color;
                let location = gl.getUniformLocation(
                    backgroundProgram, "background_sh[0]");
                gl.uniform3f(location, coeff[0], coeff[1], coeff[2]);
                var dim_sh = (background_sh_degree+1)**2;
                for (var i = 1; i < dim_sh; i++) {
                    coeff = header.config.background_sh.slice(3*i-3, 3*i);
                    location = gl.getUniformLocation(
                        backgroundProgram, "background_sh["+i+"]");
                    gl.uniform3f(location, coeff[0], coeff[1], coeff[2]);
                }
                gl.uniformMatrix4fv(gl.getUniformLocation(backgroundProgram, "view"),
                    false, actualViewMatrix);
                gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
            }

            gl.useProgram(splatProgram);
            gl.uniformMatrix4fv(gl.getUniformLocation(splatProgram, "view"),
                false, actualViewMatrix);
            gl.uniform1i(gl.getUniformLocation(splatProgram, "u_use_aniso"),
                document.getElementById("checkbox-aniso").checked);
            gl.uniform2i(gl.getUniformLocation(splatProgram, "u_sh_config"),
                document.getElementById("checkbox-sh").checked,
                header.config.sh_degree);
            gl.uniform3i(gl.getUniformLocation(splatProgram, "u_ch_config"),
                document.getElementById("checkbox-ch").checked,
                header.config.ch_degree_r,
                header.config.ch_degree_phi);
            gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);

            renderNeeded = false;
        } else if (!(vertexCount > 0)) {
            gl.clear(gl.COLOR_BUFFER_BIT);
            document.getElementById("spinner").style.display = "";
            start = Date.now() + 2000;
        }

        const progress = (100 * bytesRead / ssplatData.length);
        if (progress < 100) {
            document.getElementById("progress").style.width = progress + "%";
        } else {
            document.getElementById("progress").style.display = "none";
        }
        fps.innerText = Math.round(avgFps) + " fps";
        lastFrame = now;
        requestAnimationFrame(frame);
    };

    frame();

    window.addEventListener("hashchange", (e) => {
        try {
            viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
            carousel = false;
        } catch (err) {}
    });

    let inputs = document.getElementsByTagName("input");
    for (var i = 0; i < inputs.length; i++) {
        inputs[i].addEventListener("input", () => { renderNeeded = true; });
    }

    const preventDefault = (e) => {
        e.preventDefault();
        e.stopPropagation();
    };

    let header = null;
    let headerLength = -1;
    let lastBytesRead = -1;
    let stopLoading = false;

    while (true) {
        const { done, value } = await reader.read();
        if (done || stopLoading) break;

        ssplatData.set(value, bytesRead);
        bytesRead += value.length;

        if (header === null) {
            if (bytesRead < 8)
                continue;
            let arrayint = new Uint32Array(ssplatData.buffer);
            let arraybyte = new Uint8Array(ssplatData.buffer);
            if (arrayint[0] != 1953263731)
                break;
            headerLength = arrayint[1]+8;
            if (bytesRead < headerLength)
                continue;
            header = new TextDecoder().decode(arraybyte.slice(8, headerLength));
            header = JSON.parse(header.trim());
        }

        if (bytesRead > lastBytesRead) {
            // worker.postMessage(unpackModel(
            //     header, ssplatData.buffer.slice(headerLength, bytesRead)));
            lastBytesRead = bytesRead;
            renderNeeded = true;
        }
    }
    if (!stopLoading) {
        worker.postMessage(unpackModel(
            header, ssplatData.buffer.slice(headerLength, bytesRead)));
        renderNeeded = true;
    }
}


// Function to load a shader file
async function loadShaderFile(url) {
    url = url + "?nocache=" + Math.floor(Date.now() / 1000);
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error(`Failed to load shader file: ${url}`);
    }
    return await response.text();
}

// Function to load both shaders and call main()
async function loadShadersAndInit() {
    try {
        const [
            utilsSource,
            vertexSource, fragmentSource,
            backgroundSource
        ] = await Promise.all([
            loadShaderFile('shader-utils.glsl'),
            loadShaderFile('shader-vert.glsl'),
            loadShaderFile('shader-frag.glsl'),
            loadShaderFile('shader-background.glsl')
        ]);

        function shaderInclude(s) {
            return s.replace('\n#include "shader-utils.glsl"\n', '\n'+utilsSource+'\n');
        }

        vertexShaderSource = shaderInclude(vertexSource);
        fragmentShaderSource = shaderInclude(fragmentSource);
        backgroundShaderSource = shaderInclude(backgroundSource);
  
        // Call the main function after both shaders are loaded
        main().catch((err) => {
            document.getElementById("spinner").style.display = "none";
            document.getElementById("message").innerText = err.toString();
            console.error(err);
        });

    } catch (err) {
        document.getElementById("spinner").style.display = "none";
        document.getElementById("message").innerText = err.toString();
        console.error(err);
    }
}


window.addEventListener("load", () => {
    CameraPresets.init();
    document.getElementById("camera-preset-container").appendChild(
        CameraPresets.createSelector());

    window.Console = {
        log: (msg) => {
            let container = document.getElementById("console");
            let span = document.createElement("pre");
            span.textContent = new String(msg);
            container.appendChild(span);
        }
    };
    window.addEventListener("error", (e) => {
        console.log(e.error);
        Console.log([e.lineno, e.colno, e.message].join(' ') + '\n' + e.error.stack);
    });
    createModule().then(module => {
        Module = module;
        loadShadersAndInit();
    });
});
