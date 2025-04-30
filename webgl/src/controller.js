"use strict";

function quat2rotmat(quat) {
    const w = quat[0], x = quat[1], y = quat[2], z = quat[3];
    return [
        [1 - 2 * (y*y + z*z), 2 * (x*y + w*z), 2 * (x*z - w*y)],
        [2 * (x*y - w*z), 1 - 2 * (x*x + z*z), 2 * (y*z + w*x)],
        [2 * (x*z + w*y), 2 * (y*z - w*x), 1 - 2 * (x*x + y*y)]
    ];
}

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


let defaultViewMatrix = [
    1.0, 0.0, 0.0, 0.0,
    0.0, -0.199, 0.98, 0.0,
    0.0, -0.98, -0.199, 0.0,
    0.0, 0.0, 1.0, 1.0
];
// defaultViewMatrix = [-0.99396310067939, -0.0028527761475104297, 0.10967788056938267, 0, -0.10947236494874289, -0.04064304253875142, -0.9931590623621941, 0, 0.007290894500575474, -0.9991701580618988, 0.040085384110577396, 0, -0.008981957487634048, -0.019741792911632222, 0.80080398369808, 0.9999999999999599];
// defaultViewMatrix = [0.33168228758220725, 0.1352309888257432, 0.9336489201265321, 0, -0.9433629759026162, 0.039895202904019586, 0.3293551310485557, 0, 0.007290894500575899, -0.9900110870320623, 0.1408044397597277, 0, -0.13138587004540464, -0.07371462201299227, 0.49180482237051854, 0.9999999999998919];
let viewMatrix = defaultViewMatrix;


function ViewportController() {
    const canvas = document.getElementById("canvas");
    const fps = document.getElementById("fps");

    this.carousel = true;
    this.projectionMatrix = null;
    this.actualViewMatrix = null;
    this.renderNeeded = true;

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
        this.renderNeeded = true;
    });
    window.addEventListener("keyup", (e) => {
        activeKeys = activeKeys.filter((k) => k !== e.code);
        this.renderNeeded = true;
    });
    window.addEventListener("blur", () => {
        activeKeys = [];
        this.renderNeeded = true;
    });

    window.addEventListener(
        "wheel",
        (e) => {
            this.carousel = false;
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
            this.renderNeeded = true;
        },
        { passive: false },
    );

    let startX, startY, down;
    canvas.addEventListener("mousedown", (e) => {
        this.carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = e.ctrlKey || e.metaKey ? 2 : 1;
        this.renderNeeded = true;
    });
    canvas.addEventListener("contextmenu", (e) => {
        this.carousel = false;
        e.preventDefault();
        startX = e.clientX;
        startY = e.clientY;
        down = 2;
        this.renderNeeded = true;
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
        this.renderNeeded = true;
    });
    canvas.addEventListener("mouseup", (e) => {
        e.preventDefault();
        down = false;
        startX = 0;
        startY = 0;
        this.renderNeeded = true;
    });

    let altX = 0,
        altY = 0;
    canvas.addEventListener(
        "touchstart",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1) {
                this.carousel = false;
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
                down = 1;
            } else if (e.touches.length === 2) {
                // console.log('beep')
                this.carousel = false;
                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
                down = 1;
            }
            this.renderNeeded = true;
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
            this.renderNeeded = true;
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
            this.renderNeeded = true;
        },
        { passive: false },
    );

    let jumpDelta = 0;

    window.addEventListener("gamepadconnected", (e) => {
        const gp = navigator.getGamepads()[e.gamepad.index];
        console.log(
            `Gamepad connected at index ${gp.index}: ${gp.id}. It has ${gp.buttons.length} buttons and ${gp.axes.length} axes.`,
        );
        this.renderNeeded = true;
    });
    window.addEventListener("gamepaddisconnected", (e) => {
        console.log("Gamepad disconnected");
        this.renderNeeded = true;
    });


    let previousCamera = "";
    this.onFrame = function(now) {
        let oldRenderNeeded = this.renderNeeded;

        let inv = invert4(viewMatrix);
        let shiftKey = activeKeys.includes("Shift") || activeKeys.includes("ShiftLeft") || activeKeys.includes("ShiftRight");
        let vsc = getVelocityScale();

        if (/^Key[ADQEWSJKLI]$/.test(activeKeys) ||
            /^Arrow/.test(activeKeys)) this.renderNeeded = true;

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
                this.carousel = false;
            }
            if (Math.abs(gamepad.axes[1]) > axisThreshold) {
                inv = translate4(inv, 0, 0, -moveSpeed * gamepad.axes[1]);
                this.carousel = false;
            }
            if(gamepad.buttons[12].pressed || gamepad.buttons[13].pressed){
                inv = translate4(inv, 0, -moveSpeed*(gamepad.buttons[12].pressed - gamepad.buttons[13].pressed), 0);
                this.carousel = false;
            }
            if(gamepad.buttons[14].pressed || gamepad.buttons[15].pressed){
                inv = translate4(inv, -moveSpeed*(gamepad.buttons[14].pressed - gamepad.buttons[15].pressed), 0, 0);
                this.carousel = false;
            }

            // Assuming the right stick controls rotation (axes 2 and 3)
            if (Math.abs(gamepad.axes[2]) > axisThreshold) {
                inv = rotate4(inv, rotateSpeed * gamepad.axes[2], 0, 1, 0);
                this.carousel = false;
            }
            if (Math.abs(gamepad.axes[3]) > axisThreshold) {
                inv = rotate4(inv, -rotateSpeed * gamepad.axes[3], 1, 0, 0);
                this.carousel = false;
            }

            let tiltAxis = gamepad.buttons[6].value - gamepad.buttons[7].value;
            if (Math.abs(tiltAxis) > axisThreshold) {
                inv = rotate4(inv, rotateSpeed * tiltAxis, 0, 0, 1);
                this.carousel = false;
            }
            if (gamepad.buttons[0].pressed) {
                isJumping = true;
                this.carousel = false;
            }
            if(gamepad.buttons[3].pressed){
                this.carousel = true;
            }

            if (JSON.stringify(inv) != inv0)
                this.renderNeeded = true;
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

        if (this.carousel) {
            let inv = invert4(defaultViewMatrix);

            const t = Math.sin((Date.now() - start) / 5000);
            inv = translate4(inv, 0.2 * t, 0, 0.1 * (1 - Math.cos(t)));
            inv = rotate4(inv, -0.1 * t, 0, 1, 0);

            viewMatrix = invert4(inv);
            this.renderNeeded = true;
        }

        let jumpDeltaNew = jumpDelta;
        if (isJumping) {
            jumpDeltaNew = Math.min(1, jumpDelta + 0.05*vsc);
        } else {
            jumpDeltaNew = Math.max(0, jumpDelta - 0.05*vsc);
        }
        if (jumpDeltaNew != jumpDelta) {
            jumpDelta = jumpDeltaNew;
            this.renderNeeded = true;
        }

        let inv2 = invert4(viewMatrix);
        inv2 = translate4(inv2, 0, -jumpDelta, 0);
        inv2 = rotate4(inv2, -0.1 * jumpDelta, 1, 0, 0);
        this.actualViewMatrix = invert4(inv2);

        let camera = CameraPresets.camera;
        if (camera !== previousCamera)
            this.renderNeeded = true;
        previousCamera = camera;

        const currentFps = 1000 / (now - lastFrame) || 0;
        if (this.renderNeeded)
            avgFps = avgFps * 0.8 + currentFps * 0.2;
        fps.innerText = Math.round(avgFps) + " fps";
        lastFrame = now;
    }

    window.addEventListener("hashchange", (e) => {
        try {
            viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
            this.carousel = false;
        } catch (err) {}
    });

    let inputs = document.getElementsByTagName("input");
    for (var i = 0; i < inputs.length; i++) {
        inputs[i].addEventListener("input", () => { this.renderNeeded = true; });
    }

}
