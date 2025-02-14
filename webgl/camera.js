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

CameraPresets["$"] = null;


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
            dist_coeffs: [0.0, 0.0, 0.0, 0.0],
            model: -1,
        };
        if (camera.model == "OPENCV") {
            camera_s.model = 0;
            camera_s.dist_coeffs = [camera.k1, camera.k2, camera.p1, camera.p2];
        }
        if (camera.model == "OPENCV_FISHEYE") {
            camera_s.model = 1;
            camera_s.dist_coeffs = [camera.k1, camera.k2, camera.k3, camera.k4];
        }
        CameraPresets.camera = camera_s;
    }
    select.addEventListener("input", updateSelectedCamera);
    window.addEventListener("resize", updateSelectedCamera);
    select.value = "s21";
    updateSelectedCamera();
    return select;
}
