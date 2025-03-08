"use strict";


let genericVertexShaderSource = `#version 300 es
precision highp float;
layout(location = 0) in vec2 vertexPosition;
void main() { gl_Position = vec4(vertexPosition, 0.0, 1.0); }
`
let vertexShaderSource = "";
let fragmentShaderSource = "";
let backgroundShaderSource = "";



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
            loadShaderFile('src/shader-utils.glsl'),
            loadShaderFile('src/shader-vert.glsl'),
            loadShaderFile('src/shader-frag.glsl'),
            loadShaderFile('src/shader-background.glsl')
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


function createShaderProgram(gl, vsSource, fsSource) {
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
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texImage2D(
            gl.TEXTURE_2D, 0, gl.RGBA32UI,
            basewidth, baseheight, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_INT, basedata);
    }

    this.updateChTexture = function(chdata, chwidth, chheight) {
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
    }

    this.updateShTexture = function(shdata, shwidth, shheight) {
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
