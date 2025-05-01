"use strict";

async function main() {
    let viewportController = new ViewportController(null);

    const params = new URLSearchParams(location.search);
    try {
        viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
        viewportController.carousel = false;
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

    let RenderModes = {
        default: {
            createWorker: Worker.createWorker,
            Renderer: RasterRenderer,
        },
        rt: {
            createWorker: Worker.createRayTracingWorker,
            Renderer: RayTracingRenderer
        },
        pps: {
            createWorker: Worker.createPerPixelSortingWorker,
            Renderer: PerPixelSortingRenderer
        },
    };
    let renderMode = RenderModes[params.get("renderer") || "default"];

    renderMode.createWorker(window);

    const canvas = document.getElementById("canvas");

    const gl = canvas.getContext("webgl2", {
        antialias: false,
    });
    let renderer = new renderMode.Renderer(gl, viewportController);

    let vertexCount = 0;

    window.addEventListener("message", (e) => {
        if (e.data.vertexCount) {
            vertexCount = e.data.vertexCount;
        }
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
            renderer.updateBaseTexture(basedata, basewidth, baseheight);
        } else if (e.data.chTexture) {
            const { chdata, chwidth, chheight } = e.data.chTexture;
            renderer.updateChTexture(chdata, chwidth, chheight);
        } else if (e.data.shTexture) {
            const { shdata, shwidth, shheight } = e.data.shTexture;
            renderer.updateShTexture(shdata, shwidth, shheight);
        } else if (e.data.depthIndex) {
            const { depthIndex, viewProj } = e.data;
            renderer.updateDepthIndex(depthIndex, viewProj);
        } else if (e.data.bvhTexture) {
            const { splatArray, aabbArray } = e.data.bvhTexture;
            renderer.updateBvhTexture(splatArray, aabbArray);
        }
    });

    let bytesRead = 0;
    const frame = (now) => {
        viewportController.onFrame(now);
        let viewProj = multiply4(viewportController.projectionMatrix, viewportController.actualViewMatrix);
        window.postMessage({ view: viewProj });

        // viewportController.renderNeeded = true;
        if (vertexCount > 0 && viewportController.renderNeeded) {
            document.getElementById("spinner").style.display = "none";

            let sh_degree = header.config.sh_degree;
            document.getElementById("checkbox-sh").disabled = (sh_degree == 0);
            let ch_degree_r = header.config.ch_degree_r;
            document.getElementById("checkbox-ch").disabled = (ch_degree_r == 0);

            renderer.onFrame(header, vertexCount);

            viewportController.renderNeeded = false;
        } else if (!(vertexCount > 0)) {
            gl.clear(gl.COLOR_BUFFER_BIT);
            document.getElementById("spinner").style.display = "";
        }

        const progress = (100 * bytesRead / ssplatData.length);
        if (progress < 100) {
            document.getElementById("progress").style.width = progress + "%";
        } else {
            document.getElementById("progress").style.display = "none";
        }
        requestAnimationFrame(frame);
    };

    frame();

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
            // window.postMessage(unpackModel(
            //     header, ssplatData.buffer.slice(headerLength, bytesRead)));
            lastBytesRead = bytesRead;
            viewportController.renderNeeded = true;
        }
    }
    if (!stopLoading) {
        window.postMessage(Worker.unpackModel(
            header, ssplatData.buffer.slice(headerLength, bytesRead)));
        viewportController.renderNeeded = true;

        // check deprecated features
        let comps = header.primitives.base.componentViews;
        for (var i in comps) {
            if (comps[i].key == "anisotropies") {
                document.getElementById("deprecation-warning").style.display = null;
            }
        }
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
    createWASMModule().then(module => {
        Worker.wasmModule = module;
        loadShadersAndInit();
    });
});
