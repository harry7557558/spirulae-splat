"use strict";


let genericVertexShaderSource = `#version 300 es
precision highp float;
layout(location = 0) in vec2 vertexPosition;
void main() { gl_Position = vec4(vertexPosition, 0.0, 1.0); }
`;
let genericFragmentShaderSource = `#version 300 es
precision highp float;
flat in vec4 vColor;
void main() { gl_FragColor = vColor; }
`;
let backgroundShaderSource = "";
// default
let vertexShaderSource = "";
let fragmentShaderSource = "";
// ray tracing
let rtFragmentShaderSource = "";
// per pixel sorting
let ppsProjShaderSource = "";
let ppsProjShaderFragSource = "";
let ppsRasterShaderSource = "";


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
            utilsSource, utilsShSource, utilsChSource, backgroundSource,
            vertexSource, fragmentSource,
            rtFragmentSource,
            ppsProjSource, ppsProjFragSource, ppsRasterSource,
        ] = await Promise.all([
            loadShaderFile('src/shader-utils.glsl'),
            loadShaderFile('src/shader-utils-sh.glsl'),
            loadShaderFile('src/shader-utils-ch.glsl'),
            loadShaderFile('src/shader-background.glsl'),
            // default
            loadShaderFile('src/shader-vert.glsl'),
            loadShaderFile('src/shader-frag.glsl'),
            // ray tracing
            loadShaderFile('src/BVH_Ray_Tracing_Fragment.glsl'),
            // per pixel sorting
            loadShaderFile('src/shader-pps-projection.glsl'),
            loadShaderFile('src/shader-pps-projection-frag.glsl'),
            loadShaderFile('src/shader-pps-rasterize.glsl'),
        ]);

        function shaderInclude(s) {
            s = s.replace('\n#include "shader-utils.glsl"\n', '\n'+utilsSource+'\n');
            s = s.replace('\n#include "shader-utils-sh.glsl"\n', '\n'+utilsShSource+'\n');
            s = s.replace('\n#include "shader-utils-ch.glsl"\n', '\n'+utilsChSource+'\n');
            return s;
        }

        backgroundShaderSource = shaderInclude(backgroundSource);
        vertexShaderSource = shaderInclude(vertexSource);
        fragmentShaderSource = shaderInclude(fragmentSource);
        rtFragmentShaderSource = shaderInclude(rtFragmentSource);
        ppsProjShaderSource = shaderInclude(ppsProjSource);
        ppsProjShaderFragSource = shaderInclude(ppsProjFragSource);
        ppsRasterShaderSource = shaderInclude(ppsRasterSource);

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


function createShaderProgram(gl, vsSource, fsSource) {
    function onerror(msg) {
        console.error(msg);
        throw new Error(msg);
    }

    const vertexShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertexShader, vsSource);
    gl.compileShader(vertexShader);
    if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS)) {
        onerror(gl.getShaderInfoLog(vertexShader));
    }

    const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragmentShader, fsSource);
    gl.compileShader(fragmentShader);
    if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS))
        onerror(gl.getShaderInfoLog(fragmentShader));

    const program = gl.createProgram();
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    gl.useProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS))
        onerror(gl.getProgramInfoLog(program));

    return program;
}

function setDataTextureParameters(gl) {
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
}



function RasterRenderer(gl, viewportController) {
    const downsample = 1;

    // create shader programs
    const splatProgram = createShaderProgram(gl, vertexShaderSource, fragmentShaderSource);
    const backgroundProgram = createShaderProgram(gl, genericVertexShaderSource, backgroundShaderSource);

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

    this.setUniforms = function(program=null) {
        if (program === null) {
            this.setUniforms(splatProgram);
            this.setUniforms(backgroundProgram);
            return;
        }
        let camera = CameraPresets.camera;
        gl.useProgram(program);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "projection"),
            false, viewportController.projectionMatrix);
        gl.uniform2fv(gl.getUniformLocation(program, "focal"),
            new Float32Array([camera.fx/downsample, camera.fy/downsample]));
        gl.uniform2fv(gl.getUniformLocation(program, "viewport"),
            new Float32Array([innerWidth/downsample, innerHeight/downsample]));
        gl.uniform1i(gl.getUniformLocation(program, "camera_model"), camera.model);
        gl.uniform4fv(gl.getUniformLocation(program, "distortion"),
            new Float32Array(camera.dist_coeffs));
        gl.uniform4fv(gl.getUniformLocation(program, "undistortion"),
            new Float32Array(camera.undist_coeffs));
    }

    const resize = () => {
        let camera = CameraPresets.camera;

        viewportController.projectionMatrix = getProjectionMatrix(
            camera.fx, camera.fy, innerWidth, innerHeight);

        gl.canvas.width = Math.round(innerWidth / downsample);
        gl.canvas.height = Math.round(innerHeight / downsample);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

        this.setUniforms();
        viewportController.renderNeeded = true;
    };

    window.addEventListener("resize", function() {
        setTimeout(resize, 1);
    });
    resize();

    this.updateBaseTexture = function(basedata, basewidth, baseheight) {
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, baseTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            basewidth, baseheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, basedata);
    }

    this.updateChTexture = function(chdata, chwidth, chheight) {
        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, chTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            chwidth, chheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, chdata);
    }

    this.updateShTexture = function(shdata, shwidth, shheight) {
        gl.activeTexture(gl.TEXTURE2);
        gl.bindTexture(gl.TEXTURE_2D, shTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            shwidth, shheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, shdata);
    }

    this.updateDepthIndex = function(depthIndex, viewProj) {
        gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);
    }

    this.onFrame = function(header, vertexCount) {
        this.setUniforms();

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
                false, viewportController.actualViewMatrix);
            gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
        }

        gl.useProgram(splatProgram);
        gl.uniformMatrix4fv(gl.getUniformLocation(splatProgram, "view"),
            false, viewportController.actualViewMatrix);
        gl.uniform2i(gl.getUniformLocation(splatProgram, "u_sh_config"),
            document.getElementById("checkbox-sh").checked,
            header.config.sh_degree);
        gl.uniform3i(gl.getUniformLocation(splatProgram, "u_ch_config"),
            document.getElementById("checkbox-ch").checked,
            header.config.ch_degree_r,
            header.config.ch_degree_phi);
        gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);
    }
}



function RayTracingRenderer(gl, viewportController) {
    const downsample = 4.0;

    // create shader programs
    const splatProgram = createShaderProgram(gl, genericVertexShaderSource, rtFragmentShaderSource);
    const backgroundProgram = createShaderProgram(gl, genericVertexShaderSource, backgroundShaderSource);

    // depth test and blend
    gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    // gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
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

    var bvhTriTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, bvhTriTexture);
    var uBvhTriTextureLocation = gl.getUniformLocation(splatProgram, "tTriangleTexture");
    gl.uniform1i(uBvhTriTextureLocation, 0);

    var bvhAabbTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, bvhAabbTexture);
    var uBvhAabbTextureLocation = gl.getUniformLocation(splatProgram, "tAABBTexture");
    gl.uniform1i(uBvhAabbTextureLocation, 1);

    // indices
    const indexBuffer = gl.createBuffer();
    const a_index = gl.getAttribLocation(splatProgram, "index");
    gl.enableVertexAttribArray(a_index);
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
    gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
    gl.vertexAttribDivisor(a_index, 1);

    this.setUniforms = function(program=null) {
        if (program === null) {
            this.setUniforms(splatProgram);
            this.setUniforms(backgroundProgram);
            return;
        }
        let camera = CameraPresets.camera;
        gl.useProgram(program);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "projection"),
            false, viewportController.projectionMatrix);
        gl.uniform2fv(gl.getUniformLocation(program, "focal"),
            new Float32Array([camera.fx/downsample, camera.fy/downsample]));
        gl.uniform2fv(gl.getUniformLocation(program, "viewport"),
            new Float32Array([innerWidth/downsample, innerHeight/downsample]));
        gl.uniform1i(gl.getUniformLocation(program, "camera_model"), camera.model);
        gl.uniform4fv(gl.getUniformLocation(program, "distortion"),
            new Float32Array(camera.dist_coeffs));
        gl.uniform4fv(gl.getUniformLocation(program, "undistortion"),
            new Float32Array(camera.undist_coeffs));
    }

    const resize = () => {
        let camera = CameraPresets.camera;

        viewportController.projectionMatrix = getProjectionMatrix(
            camera.fx, camera.fy, innerWidth, innerHeight);

        gl.canvas.width = Math.round(innerWidth / downsample);
        gl.canvas.height = Math.round(innerHeight / downsample);
        console.log(gl.canvas);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

        this.setUniforms();
        viewportController.renderNeeded = true;
    };

    window.addEventListener("resize", function() {
        setTimeout(resize, 1);
    });
    resize();

    this.updateBvhTexture = function(splatArray, aabbArray) {
        var size = Math.sqrt(splatArray.length / 4);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, bvhTriTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            size, size, 0, gl.RGBA_INTEGER, gl.UNSIGNED_INT,
            splatArray
        );

        size = Math.sqrt(aabbArray.length / 4);
        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, bvhAabbTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            // gl.TEXTURE_2D, 0, gl.RGBA32F,
            // size, size, 0, gl.RGBA, gl.FLOAT,
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            size, size, 0, gl.RGBA_INTEGER, gl.UNSIGNED_INT,
            aabbArray
        );
    }

    this.onFrame = function(header, vertexCount) {
        this.setUniforms();

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
                false, viewportController.actualViewMatrix);
            gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
        }

        gl.useProgram(splatProgram);
        gl.uniformMatrix4fv(gl.getUniformLocation(splatProgram, "view"),
            false, viewportController.actualViewMatrix);
        gl.uniform2i(gl.getUniformLocation(splatProgram, "u_sh_config"),
            document.getElementById("checkbox-sh").checked,
            header.config.sh_degree);
        gl.uniform3i(gl.getUniformLocation(splatProgram, "u_ch_config"),
            document.getElementById("checkbox-ch").checked,
            header.config.ch_degree_r,
            header.config.ch_degree_phi);
        gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
    }
}



function PerPixelSortingRenderer(gl, viewportController) {
    const downsample = 2;

    // create shader programs
    const projProgram = createShaderProgram(gl, ppsProjShaderSource, ppsProjShaderFragSource);
    const rasterProgram = createShaderProgram(gl, genericVertexShaderSource, ppsRasterShaderSource);
    const backgroundProgram = createShaderProgram(gl, genericVertexShaderSource, backgroundShaderSource);

    // positions
    const triangleVertices = new Float32Array([-1, -1, 1, -1, 1, 1, -1, 1]);
    const vertexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);

    gl.useProgram(projProgram);
    const a_position = gl.getAttribLocation(projProgram, "vertexPosition");
    gl.enableVertexAttribArray(a_position);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(a_position, 2, gl.FLOAT, false, 0, 0);

    gl.useProgram(backgroundProgram);
    const bg_position = gl.getAttribLocation(backgroundProgram, "vertexPosition");
    gl.enableVertexAttribArray(bg_position);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(bg_position, 2, gl.FLOAT, false, 0, 0);

    // textures
    let baseTexture = gl.createTexture();
    let shTexture = gl.createTexture();
    let chTexture = gl.createTexture();
    let psaTexture = gl.createTexture();
    let intTexture = gl.createTexture();

    if (!gl.getExtension('EXT_color_buffer_float'))
        alert("Failed to get extension EXT_color_buffer_float");
    let projTexture = null, projFramebuffer = null;
    function createProjTarget(width, height) {
        if (projTexture === null)
            gl.deleteTexture(projTexture);
        if (projFramebuffer === null)
            gl.deleteFramebuffer(projFramebuffer);
        projTexture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, projTexture);
        // gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, width, height, 0, gl.RGBA, gl.FLOAT, null);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32UI, width, height, 0, gl.RGBA_INTEGER, gl.UNSIGNED_INT, null);
        setDataTextureParameters(gl);
        projFramebuffer = gl.createFramebuffer();
        gl.bindFramebuffer(gl.FRAMEBUFFER, projFramebuffer);
        gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, projTexture, 0);
    }

    // set uniforms
    gl.useProgram(projProgram);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, baseTexture);
    gl.uniform1i(gl.getUniformLocation(projProgram, "u_base_texture"), 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, shTexture);
    gl.uniform1i(gl.getUniformLocation(projProgram, "u_sh_texture"), 1);

    gl.useProgram(rasterProgram);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, psaTexture);
    gl.uniform1i(gl.getUniformLocation(rasterProgram, "u_psa_texture"), 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, intTexture);
    gl.uniform1i(gl.getUniformLocation(rasterProgram, "u_int_texture"), 1);
    gl.activeTexture(gl.TEXTURE3);
    gl.bindTexture(gl.TEXTURE_2D, chTexture);
    gl.uniform1i(gl.getUniformLocation(rasterProgram, "u_ch_texture"), 3);

    // indices
    const indexBuffer = gl.createBuffer();
    const a_index = gl.getAttribLocation(projProgram, "index");
    gl.enableVertexAttribArray(a_index);
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
    gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
    gl.vertexAttribDivisor(a_index, 1);

    this.setUniforms = function(program=null) {
        if (program === null) {
            this.setUniforms(projProgram);
            this.setUniforms(rasterProgram);
            this.setUniforms(backgroundProgram);
            return;
        }
        let camera = CameraPresets.camera;
        gl.useProgram(program);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "projection"),
            false, viewportController.projectionMatrix);
        gl.uniform2fv(gl.getUniformLocation(program, "focal"),
            new Float32Array([camera.fx/downsample, camera.fy/downsample]));
        gl.uniform2fv(gl.getUniformLocation(program, "viewport"),
            new Float32Array([innerWidth/downsample, innerHeight/downsample]));
        gl.uniform1i(gl.getUniformLocation(program, "camera_model"), camera.model);
        gl.uniform4fv(gl.getUniformLocation(program, "distortion"),
            new Float32Array(camera.dist_coeffs));
        gl.uniform4fv(gl.getUniformLocation(program, "undistortion"),
            new Float32Array(camera.undist_coeffs));
    }

    const resize = () => {
        let camera = CameraPresets.camera;

        viewportController.projectionMatrix = getProjectionMatrix(
            camera.fx, camera.fy, innerWidth, innerHeight);

        gl.canvas.width = Math.round(innerWidth / downsample);
        gl.canvas.height = Math.round(innerHeight / downsample);

        this.setUniforms();
        viewportController.renderNeeded = true;
    };

    window.addEventListener("resize", function() {
        setTimeout(resize, 1);
    });
    resize();

    this.updateBaseTexture = function(basedata, basewidth, baseheight) {
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, baseTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            basewidth, baseheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, basedata);
    }

    this.updateChTexture = function(chdata, chwidth, chheight) {
        gl.activeTexture(gl.TEXTURE3);
        gl.bindTexture(gl.TEXTURE_2D, chTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            chwidth, chheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, chdata);
    }

    this.updateShTexture = function(shdata, shwidth, shheight) {
        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, shTexture);
        setDataTextureParameters(gl);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            shwidth, shheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, shdata);
    }

    let depthIndex = [];
    let inverseDepthIndex = [];
    this.updateDepthIndex = function(depthIndex_, viewProj) {
        depthIndex = depthIndex_;
        gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);

        inverseDepthIndex = new Uint32Array(depthIndex.length);
        for (var i = 0; i < depthIndex.length; i++)
            inverseDepthIndex[depthIndex[i]] = i;
    }

    function glSynchronize(output) {
        if (output) {
            var data = new Uint8Array(4);
            gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, data);
        }
        else {
            var data = new Uint32Array(4);
            gl.readPixels(0, 0, 1, 1, gl.RGBA_INTEGER, gl.UNSIGNED_INT, data);
        }
    }

    let oldProjHeight = 0;
    this.onFrame = function(header, vertexCount) {
        this.setUniforms();
        gl.disable(gl.DEPTH_TEST);
        gl.disable(gl.BLEND);

        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

        /* projection */

        gl.useProgram(projProgram);
        let projWidth = 1536;
        let projHeight = Math.ceil(3 * vertexCount / projWidth);
        if (projHeight != oldProjHeight) {
            createProjTarget(projWidth, projHeight);
            oldProjHeight = projHeight;
            console.log("Update projTarget:", projWidth, projHeight);
        }
        gl.bindFramebuffer(gl.FRAMEBUFFER, projFramebuffer);
        gl.viewport(0, 0, projWidth, projHeight);
        gl.uniformMatrix4fv(gl.getUniformLocation(projProgram, "view"),
            false, viewportController.actualViewMatrix);
        gl.uniform2i(gl.getUniformLocation(projProgram, "u_sh_config"),
            document.getElementById("checkbox-sh").checked,
            header.config.sh_degree);
        gl.uniform1f(gl.getUniformLocation(projProgram, "wh_ratio"), projWidth / projHeight);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, baseTexture);

        // glSynchronize(false);
        // console.time("projection");
        gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);
        // glSynchronize(false);
        // console.timeEnd("projection");

        /* setup tiles */

        var projData = new Uint32Array(4 * projWidth * projHeight);
        // console.time("gl.readPixels");
        gl.readPixels(0, 0, projWidth, projHeight, gl.RGBA_INTEGER, gl.UNSIGNED_INT, projData);
        // console.timeEnd("gl.readPixels");

        // console.time("preparePPSTiles");
        let tiles = Worker.wasmModule.preparePPSTiles(
            vertexCount, innerWidth, innerHeight,
            projData, depthIndex
        );
        let numTilesX = tiles.numTilesX;
        let numTilesY = tiles.numTilesY;
        let intTexWidth = tiles.intTexWidth;
        let intTexHeight = tiles.intTexHeight;
        let numTilesPSA = tiles.numTilesPSA;
        let intersects = tiles.intersects;
        // console.timeEnd("preparePPSTiles");
        // kinda CPU bound, considering figuring out doing this on GPU

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, psaTexture);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI, numTilesX, numTilesY, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, numTilesPSA);
        setDataTextureParameters(gl);

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, intTexture);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI, intTexWidth, intTexHeight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, intersects);
        setDataTextureParameters(gl);

        /* rasterization */

        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);

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
                false, viewportController.actualViewMatrix);
            gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
        }

        gl.useProgram(rasterProgram);
        gl.uniformMatrix4fv(gl.getUniformLocation(rasterProgram, "view"),
            false, viewportController.actualViewMatrix);
        gl.uniform3i(gl.getUniformLocation(rasterProgram, "u_ch_config"),
            document.getElementById("checkbox-ch").checked,
            header.config.ch_degree_r,
            header.config.ch_degree_phi);
        gl.activeTexture(gl.TEXTURE2);
        gl.bindTexture(gl.TEXTURE_2D, projTexture);
        gl.uniform1i(gl.getUniformLocation(rasterProgram, "u_proj_texture"), 2);

        // glSynchronize(true);
        // console.time("rasterization");
        gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, 4);
        // glSynchronize(true);
        // console.timeEnd("rasterization");

    }
}

