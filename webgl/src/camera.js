"use strict";

let CameraPresets = {};

CameraPresets["s21"] = {
    name: "S21 Default",
    "w": 1440,
    "h": 1080,
    "fx": 1095.1245833946907,
    "fy": 1094.3580648863167,
    "cx": 720.0,
    "cy": 540.0,
    "k1": 0.05034023997067791,
    "k2": -0.043216212746082325,
    "p1": -0.0014314791258299502,
    "p2": 0.0008191675045178124,
    "model": "OPENCV",
};

CameraPresets["s21wide"] = {
    name: "S21 Wide",
    "w": 2344,
    "h": 1080,
    "fx": 1097.1793294408487,
    "fy": 1095.016791999343,
    "cx": 1153.0445385978421,
    "cy": 540.0554138732397,
    "k1": -0.0045764939109361295,
    "k2": -0.007780908612130821,
    "p1": -0.0003586612274147588,
    "p2": 0.002615958332449294,
    "model": "OPENCV"
};

CameraPresets["imx219"] = {
    name: "IMX219",
    "w": 1280,
    "h": 720,
    "fx": 575.681250922668,
    "fy": 576.1624750819685,
    "cx": 607.5612147329209,
    "cy": 361.0335080022911,
    "k1": 0.044675886593263776,
    "k2": -0.012454840182315909,
    "k3": 0.008726394633797787,
    "k4": -0.0032802897607147976,
    "model": "OPENCV_FISHEYE",
};

if (false)
CameraPresets["imx219ocv"] = {
    name: "IMX219 (OpenCV)",
    "w": 1280,
    "h": 720,
    "fx": 577.3855747940892,
    "fy": 578.4368064311739,
    "cx": 633.6586692524488,
    "cy": 360.55821768989074,
    "k1": -0.2161935881621547,
    "k2": 0.0350491696401335,
    "p1": 0.00027434574590604104,
    "p2": -0.00018872208213217676,
    "model": "OPENCV",
};

CameraPresets["t265"] = {
    name: "T265",
    "w": 848,
    "h": 800,
    "fx": 287.1207307480006,
    "fy": 287.1467554811063,
    "cx": 420.298272012518,
    "cy": 414.13518449424595,
    "k1": -0.017462005784106175,
    "k2": 0.05864943135824237,
    "k3": -0.05407212566738003,
    "k4": 0.011417049521949768,
    "model": "OPENCV_FISHEYE",
};

CameraPresets["debug"] = {
    name: "Debug",
    "w": 640,
    "h": 640,
    "fx": 200,
    "fy": 200,
    "cx": 320,
    "cy": 320,
    "k1": 0.2,
    "k2": -0.1,
    "p1": 0.05,
    "p2": -0.05,
    "model": "OPENCV",
};

CameraPresets["$"] = null;


CameraPresets.invFisheyeDistort = function(k1, k2, k3, k4) {
    `def inv_fisheye_distort(k1, k2, k3, k4):
        N = 128
        theta = np.pi/2 * (np.arange(N)+0.5)/N
        theta2 = theta*theta
        r = theta*(1+theta2*(k1+theta2*(k2+theta2*(k3+theta2*k4))))
        A = np.stack([r**3, r**5, r**7, r**9])
        l1, l2, l3, l4 = np.linalg.solve(A@A.T, A@(theta-r))
        return l1, l2, l3, l4
    `;

    const N = 128, n = 4;
    let theta = Array.from({length: N}, (_, i) => Math.PI/2 * (i + 0.5)/N);
    let r = theta.map(t => {
        let t2 = t*t, poly = 1 + t2*(k1 + t2*(k2 + t2*(k3 + t2*k4)));
        return t * poly;
    });
    
    // Create matrix A (4x128)
    let A = [[], [], [], []];
    for (let i = 0; i < N; i++) {
        let ri = r[i], r2 = ri*ri;
        A[0][i] = ri*r2;
        A[1][i] = A[0][i]*r2;
        A[2][i] = A[1][i]*r2;
        A[3][i] = A[2][i]*r2;
    }

    // Compute M = A*A^T and b = A*(theta - r)
    let M = Array(n).fill().map(() => Array(n).fill(0));
    let b = Array(n).fill(0);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) 
            for (let k = 0; k < N; k++) 
                M[i][j] += A[i][k] * A[j][k];
        for (let k = 0; k < N; k++)
            b[i] += A[i][k] * (theta[k] - r[k]);
    }

    // Cholesky decomposition
    let L = Array(n).fill().map(() => Array(n).fill(0));
    for (let i = 0; i < n; i++) 
        for (let j = 0; j <= i; j++) {
            let s = 0;
            for (let k = 0; k < j; k++) s += L[i][k] * L[j][k];
            L[i][j] = i === j ? Math.sqrt(M[i][i] - s) : (M[i][j] - s)/L[j][j];
        }
    
    // Forward substitution (Ly = b)
    let y = Array(n).fill(0);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < i; j++) y[i] -= L[i][j] * y[j];
        y[i] = (b[i] + y[i]) / L[i][i];
    }

    // Backward substitution (L^Tx = y)
    let x = Array(n).fill(0);
    for (let i = n-1; i >= 0; i--) {
        for (let j = i+1; j < n; j++) x[i] -= L[j][i] * x[j];
        x[i] = (y[i] + x[i]) / L[i][i];
    }

    return x;
}

CameraPresets.init = function() {
    for (var id in CameraPresets) {
        if (id == '$') break;
        let camera = CameraPresets[id];

        // add undistortion coefficients
        if (camera.model == "OPENCV") {
            camera.dist_coeffs = [camera.k1, camera.k2, camera.p1, camera.p2];
            camera.undist_coeffs = [0.0, 0.0, 0.0, 0.0];
        }
        else if (camera.model == "OPENCV_FISHEYE") {
            var l = CameraPresets.invFisheyeDistort(camera.k1, camera.k2, camera.k3, camera.k4);
            camera.dist_coeffs = [camera.k1, camera.k2, camera.k3, camera.k4];
            camera.undist_coeffs = l;
        }
        else {
            camera.dist_coeffs = [0.0, 0.0, 0.0, 0.0];
            camera.undist_coeffs = [0.0, 0.0, 0.0, 0.0];
        }
    }
};


CameraPresets.createSelector = function() {
    let select = document.createElement("select");
    for (var id in CameraPresets) {
        if (id == '$') break;
        let option = document.createElement("option");
        option.value = id;
        option.textContent = CameraPresets[id].name;
        select.appendChild(option);
    }
    function updateSelectedCamera() {
        let id = select.value;
        let camera = CameraPresets[id];
        let size0 = Math.sqrt(window.innerWidth*window.innerHeight);
        let size1 = Math.sqrt(camera.w*camera.h);
        let camera_s = {
            w: window.innerWidth,
            h: window.innerHeight,
            fx: camera.fx * size0/size1,
            fy: camera.fy * size0/size1,
            cx: camera.cx * window.innerWidth/camera.w,
            cy: camera.cy * window.innerHeight/camera.w,
            dist_coeffs: camera.dist_coeffs,
            undist_coeffs: camera.undist_coeffs,
            model: -1,
        };
        if (camera.model == "OPENCV") {
            camera_s.model = 0;
        }
        if (camera.model == "OPENCV_FISHEYE") {
            camera_s.model = 1;
        }
        CameraPresets.camera = camera_s;
    }
    select.addEventListener("input", updateSelectedCamera);
    window.addEventListener("resize", updateSelectedCamera);
    select.value = "s21";
    updateSelectedCamera();
    return select;
}
