"use strict";


// misc help function - float/half packing

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


// misc help function - distance to ellipse from iq

function sdEllipse(x, y, a, b) {
    let px = Math.abs(x);
    let py = Math.abs(y);
    if (px > py) {
        [px, py] = [py, px];
        [a, b] = [b, a];
    }
    
    const l = b * b - a * a;
    const m = a * px / l;
    const m2 = m * m;
    const n = b * py / l;
    const n2 = n * n;
    const c = (m2 + n2 - 1.0) / 3.0;
    const c3 = c * c * c;
    const q = c3 + m2 * n2 * 2.0;
    const d = c3 + m2 * n2;
    const g = m + m * n2;
    
    let co;
    if (d < 0.0) {
        const p = Math.acos(q / c3) / 3.0;
        const s = Math.cos(p);
        const t = Math.sin(p) * Math.sqrt(3.0);
        const rx = Math.sqrt(-c * (s + t + 2.0) + m2);
        const ry = Math.sqrt(-c * (s - t + 2.0) + m2);
        co = (ry + Math.sign(l) * rx + Math.abs(g) / (rx * ry) - m) / 2.0;
    } else {
        const h = 2.0 * m * n * Math.sqrt(d);
        const s = Math.cbrt(q + h);
        const u = Math.cbrt(q - h);
        const rx = -s - u - c * 4.0 + 2.0 * m2;
        const ry = (s - u) * Math.sqrt(3.0);
        const rm = Math.hypot(rx, ry);
        const p = ry / Math.sqrt(rm - rx);
        co = (p + 2.0 * g / rm - m) / 2.0;
    }
    
    const si = Math.sqrt(1.0 - co * co);
    
    // Calculate closest point on the ellipse
    const closestX = a * co;
    const closestY = b * si;
    
    // Calculate distance and sign
    const distance = Math.hypot(closestX - px, closestY - py);
    const sign = py - closestY > 0 ? 1 : -1;
    
    return distance * sign;
}



// Worker


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
    console.log(header.primitives.base.length.toLocaleString() + " splats, "
        + header.buffer.byteLength.toLocaleString() + " bytes");
    console.log(header);

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


Worker.generateBaseTexture = function(header, base) {
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

Worker.generateCHTexture = function(header, harmonics) {
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

Worker.generateSHTexture = function(header, harmonics) {
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


Worker.createWorker = function(self) {
    let header;
    let base, harmonics;
    let vertexCount = 0;
    let viewProj;

    let lastProj = [];
    let depthIndex = new Uint32Array();
    let lastBaseVertexCount = 0,
        lastHarmonicsLength = 0;

    function runSort(viewProj) {
        if (!base) return;
        if (harmonics && harmonics.features_sh.length > lastHarmonicsLength) {
            console.time("generate CH/SH texture");
            Worker.generateCHTexture(header, harmonics);
            Worker.generateSHTexture(header, harmonics);
            console.timeEnd("generate CH/SH texture");
            lastHarmonicsLength = harmonics.features_sh.length;
        }
        vertexCount = base.opacities.length;
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
            Worker.generateBaseTexture(header, base);
            lastBaseVertexCount = vertexCount;
            console.timeEnd("generate base texture");
        }

        console.time("sort");
        let depthIndex = Worker.wasmModule.sortByDepth(
            vertexCount, base.means,
            viewProj[2], viewProj[6], viewProj[10]
        )
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


Worker.createPerPixelSortingWorker = function(self) {
    // modified from createWorker, not actually sorting anything

    let header;
    let base, harmonics;
    let vertexCount = 0;
    let viewProj;

    let lastProj = [];
    let depthIndex = new Uint32Array();
    let lastBaseVertexCount = 0,
        lastHarmonicsLength = 0;

    function runSort(viewProj) {
        if (!base) return;
        if (harmonics && harmonics.features_sh.length > lastHarmonicsLength) {
            console.time("generate CH/SH texture");
            Worker.generateCHTexture(header, harmonics);
            Worker.generateSHTexture(header, harmonics);
            console.timeEnd("generate CH/SH texture");
            lastHarmonicsLength = harmonics.features_sh.length;
        }

        vertexCount = base.opacities.length;
        if (vertexCount != lastBaseVertexCount) {
            console.time("generate base texture");
            Worker.generateBaseTexture(header, base);
            lastBaseVertexCount = vertexCount;
            console.timeEnd("generate base texture");

            depthIndex = new Uint32Array(vertexCount);
            for (var i = 0; i < vertexCount; i++)
                depthIndex[i] = i;
            lastBaseVertexCount = vertexCount;

            lastProj = viewProj;
            self.postMessage({ depthIndex, viewProj, vertexCount }, [
                depthIndex.buffer,
            ]);
        }
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


Worker.buildBVH = function(header, base) {
    let vertexCount = base.opacities.length;
    let numSplats = vertexCount;
    console.log(numSplats, "splats");

	let splatArray = new Float32Array(2048 * 2048 * 4);
    let triangleArrayUint32View = new Uint32Array(splatArray.buffer);

    const NT = 0;
    let tessellatedSplats = new Array(numSplats);
    let numAabb = 0;

	for (let i = 0; i < numSplats; i++) {

        let p = base.means.slice(3*i, 3*i+3);
        let q = base.quats.slice(4*i, 4*i+4);
        let R = quat2rotmat(q);
        let sx = Math.exp(base.scales[2*i+0]), sy = Math.exp(base.scales[2*i+1]);
        let pu = [R[0][0]*sx, R[0][1]*sx, R[0][2]*sx];
        let pv = [R[1][0]*sy, R[1][1]*sy, R[1][2]*sy];

        let m = header.config.sh_degree == 0 ? 0.28209479177387814 : 1.0;
        let b = header.config.sh_degree == 0 ? 0.5 : 0.0;
        m = 0.28209479177387814, b = 0.5;
        if (header.config.ch_degree_r > 0 || header.config.ch_degree_phi > 0)
            m *= 0.5, b *= 0.5;
        let rgba = [
            b + m * base.features_dc[3*i+0],
            b + m * base.features_dc[3*i+1],
            b + m * base.features_dc[3*i+2],
            base.opacities[i]
        ];

        let iTex = 8 * i;

        splatArray[iTex + 0] = p[0];
        splatArray[iTex + 1] = p[1];
        splatArray[iTex + 2] = p[2];
        triangleArrayUint32View[iTex + 3] = packHalf2x16(pu[0], pu[1]);
        triangleArrayUint32View[iTex + 4] = packHalf2x16(pu[2], pv[0]);
        triangleArrayUint32View[iTex + 5] = packHalf2x16(pv[1], pv[2]);
        triangleArrayUint32View[iTex + 6] = packHalf2x16(rgba[0], rgba[1]);
        triangleArrayUint32View[iTex + 7] = packHalf2x16(rgba[2], rgba[3]);

        let sr = [Math.hypot(pu[0], pv[0]), Math.hypot(pu[1], pv[1]), Math.hypot(pu[2], pv[2])];
        let r0 = 1 / (1+Math.sqrt(2)*NT);

        // tessellate splat
        // TODO: make this more compact without missing stuff
        var sc = Math.max(sx, sy);
        var boxes = [];
        for (var ti = -NT; ti <= NT; ti++) {
            for (var tj = -NT; tj <= NT; tj++) {
                var u = ti / (NT+0.5);
                var v = tj / (NT+0.5);
                // if (Math.hypot(u*sc/sx, v*sc/sy) >= 1)
                if (sdEllipse(u, v, sx/sc, sy/sc) >= r0)
                    continue;
                var c = [p[0]+u*pu[0]+v*pv[0], p[1]+u*pu[1]+v*pv[1], p[2]+u*pu[2]+v*pv[2]];
                var r = [sr[0]*r0, sr[1]*r0, sr[2]*r0];
                boxes.push([
                    c[0]-r[0], c[1]-r[1], c[2]-r[2],
                    c[0]+r[0], c[1]+r[1], c[2]+r[2],
                    c[0], c[1], c[2]
                ]);
            }
        }
        var boxesArray = new Float32Array(9*boxes.length);
        for (var k = 0; k < boxes.length; k++) {
            for (var _ = 0; _ < 9; _++)
                boxesArray[9*k+_] = boxes[k][_];
        }
        tessellatedSplats[i] = boxesArray;
        numAabb += boxes.length;
	}

	let aabbArray = new Float32Array(4096 * 4096 * 4);
	let totalWork = new Uint32Array(numAabb);
    let tmap = new Uint32Array(numAabb);
    for (var i = 0, ptr = 0; i < numSplats; i++) {
        var nt = tessellatedSplats[i].length / 9;
        // totalWork.fill(i, ptr, ptr+nt);
        aabbArray.set(tessellatedSplats[i], 9*ptr);
        for (var k = 0; k < nt; k++) {
            totalWork[ptr] = ptr;
            tmap[ptr] = i;
            ptr++;
        }
    }
    console.log(numAabb, "tesselated boxes");

	console.time("BvhGeneration");
	console.log("BvhGeneration...");
    BVH_Build_Iterative(totalWork, aabbArray);
	console.timeEnd("BvhGeneration");

    // compress AABB array
    let aabbArrayUint32View = new Uint32Array(aabbArray.buffer);
    for (var iTex = 0; iTex < aabbArray.length; iTex += 4) {
        var iTex0 = 2 * iTex;
        if (aabbArray.slice(iTex0, iTex0+8).every(x => x == 0))
            break;
        var idTriangle = aabbArray[iTex0];
        var aabbMin = aabbArray.slice(iTex0+1, iTex0+4);
        var idChild = aabbArray[iTex0+4];
        var aabbMax = aabbArray.slice(iTex0+5, iTex0+8);
        aabbArrayUint32View[iTex+0] = idTriangle < 0 ? -idChild : tmap[idTriangle];
        aabbArrayUint32View[iTex+1] = packHalf2x16(aabbMin[0], aabbMax[0]);
        aabbArrayUint32View[iTex+2] = packHalf2x16(aabbMin[1], aabbMax[1]);
        aabbArrayUint32View[iTex+3] = packHalf2x16(aabbMin[2], aabbMax[2]);
    }
    console.log(iTex/4, "bounding boxes");

    self.postMessage({
        bvhTexture: {
            splatArray: triangleArrayUint32View,
            aabbArray: aabbArrayUint32View,
        }
    }, []);
}


Worker.createRayTracingWorker = function(self) {
    let header;
    let base, harmonics;
    let vertexCount = 0;
    let viewProj;

    let lastProj = [];
    let depthIndex = new Uint32Array();
    let lastBaseVertexCount = 0,
        lastHarmonicsLength = 0;

    const throttledSort = () => {
        if (!sortRunning) {
            sortRunning = true;
            let lastView = viewProj;
            if (base && harmonics && harmonics.features_sh.length > lastHarmonicsLength) {
                lastHarmonicsLength = harmonics.features_sh.length;
            }
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
            vertexCount = base.opacities.length;
            postMessage({ base, harmonics, vertexCount });

            Worker.buildBVH(header, base);

        }
        
        else if (e.data.view) {
            viewProj = e.data.view;
            throttledSort();
        }

    });
}
