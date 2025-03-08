"use strict";

let Worker = {
    wasmModule: null
};

// config: JavaScript object
// packedData: ArrayBuffer
// buffers: Array of ArrayBuffers
Worker.unpackComponents = function(config, packedData, buffers) {
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

Worker.unpackModel = function(header, buffer) {
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
        let components = Worker.wasmModule.unpackComponents(key, primitive, bufferSlice, buffers);
        model[key] = components;
    };
    console.timeEnd("unpack model");

    model.header = header;
    return model;
}


Worker.createWorker = function(self) {
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

            // anisotropy - 2 floats, info1 w, info2 x, DEPRECATED
            basedataF[iTex + 7] = basedataF[iTex + 8] = 0.0;

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
        var chdata = Worker.wasmModule.packHarmonicTexture(
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
        var shdata = Worker.wasmModule.packHarmonicTexture(
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
