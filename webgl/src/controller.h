#include "../../spirulae_splat/splat/cuda/csrc/glm/glm/glm.hpp"
#include "../../spirulae_splat/splat/cuda/csrc/glm/glm/gtc/quaternion.hpp"

#include <set>
#include <string>
#include <array>

#include <functional>

#include <cstdio>

#include "collision.h"


template<typename T>
void copyToVector(const emscripten::val &typedArray, std::vector<T> &vec) {
    vec = emscripten::convertJSArrayToNumberVector<T>(typedArray);
    return;

    // https://github.com/emscripten-core/emscripten/issues/5519
    // doesn't work after memory growth
    unsigned int length = typedArray["length"].as<unsigned int>();
    emscripten::val heap = emscripten::val::module_property("HEAPU8");
    emscripten::val memory = heap["buffer"];

    emscripten::val memoryView = typedArray["constructor"].new_(memory, reinterpret_cast<uintptr_t>(vec.data()), length);
    
    memoryView.call<void>("set", typedArray);
}



struct ViewState {
    glm::vec3 p;
    glm::vec3 v;
    glm::quat q;
    glm::vec3 w;

    void fromMatrix(glm::mat4 m) {
        p = glm::vec3(m[3]);
        v = glm::vec3(0);
        q = glm::quat(m);
        w = glm::vec3(0);
    }

    glm::mat4 toMatrix() {
        glm::mat4 m = glm::mat4(this->q);
        m[3] = glm::vec4(this->p, 1.0f);
        return m;
    }

    ViewState operator+(const ViewState &d) {
        ViewState x;
        x.p = this->p + d.p;
        x.v = this->v + d.v;
        x.q = this->q + d.q;
        x.w = this->w + d.w;
        return x;
    }

    ViewState operator*(float dt) {
        ViewState x;
        x.p = this->p * dt;
        x.v = this->v * dt;
        x.q = this->q * dt;
        x.w = this->w * dt;
        return x;
    }

    void step(std::function<ViewState(ViewState)> f, float dt) {
        static constexpr float maxDt = 0.005f;
        if (dt > maxDt) {
            int n = int(ceilf(dt / maxDt));
            for (int i = 0; i < n; i++)
                step(f, dt / n);
            return;
        }
        // RK4 integration
        ViewState y = (*this);
        ViewState k1 = f(y);
        ViewState k2 = f(y + k1*0.5f*dt);
        ViewState k3 = f(y + k2*0.5f*dt);
        ViewState k4 = f(y + k3*dt);
        (*this) = y + (k1 + k2*2.0f + k3*2.0f + k4) * (dt/6.0f);
    }
};


class BaseController {

    static constexpr float kCollisionMinOpacity = 0.02f;
    static constexpr float kWheelTau = 0.1f;  // s

protected:
    float sceneScale = 1.0;

    std::set<std::string> activeKeys;
    glm::vec2 wheelDelta = glm::vec2(0.0f);
    glm::vec2 mousePos = glm::vec2(0.0f);
    glm::vec2 mouseDelta = glm::vec2(0.0f);
    std::array<float, 6> gamepadValues = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    std::vector<Primitive> prims;
    std::unique_ptr<BVHIntegral> collisionDetector = nullptr;

public:

    void setGaussians(int n, glm::vec3* means, glm::vec3* scales, glm::vec4* quats, float* opacities) {
        using namespace glm;

        // get scale
        vec3 sum_pos(0), sum_pos2(0);
        float sum_weight = 0.0f;
        for (int i = 0; i < n; i++) {
            float w = opacities[i];
            vec3 pos = means[i];
            sum_pos += w * pos;
            sum_pos2 += w * pos*pos;
            sum_weight += w;
        }
        sum_pos2 /= sum_weight, sum_pos /= sum_weight;
        vec3 std = sqrt(sum_pos2 - sum_pos*sum_pos);
        this->sceneScale = length(std);
        printf("sceneScale: %f\n", this->sceneScale);

        // convert to inverse cov and opacity
        prims.reserve(n);
        for (int i = 0; i < n; i++) {
            float opac = opacities[i];
            float extend = sqrtf(2.0f*logf(fmaxf(opac/kCollisionMinOpacity, 1.0f)));
            if (!(extend > 0.0f))
                continue;

            vec3 mean = means[i];
            quat q = quat(quats[i].x, quats[i].y, quats[i].z, quats[i].w);
            mat3 R = mat3(q);
            vec3 sc = scales[i];
            R[0] *= expf(sc.x); R[1] *= expf(sc.y); R[2] *= expf(sc.z);
            mat3 cov3d = R * transpose(R);
            mat3 inv_cov3d = inverse(cov3d);

            vec3 bound = extend * vec3(sqrtf(cov3d[0][0]), sqrtf(cov3d[1][1]), sqrtf(cov3d[2][2]));
            if (!isfinite(dot(bound, bound)))
                continue;

            Primitive prim = {
                mean, opac,
                {
                    -0.5f*inv_cov3d[0][0], -0.5f*inv_cov3d[1][1], -0.5f*inv_cov3d[2][2],
                    -inv_cov3d[0][1], -inv_cov3d[0][2], -inv_cov3d[1][2]
                },
                { mean-bound, mean+bound }
            };
            prims.push_back(prim);
        }
        // build collision detector
        collisionDetector = std::make_unique<BVHIntegral>(prims);
        printf("collisionDetector: %p\n", collisionDetector.get());
    }

    virtual void keydown(std::string keyCode) {
        activeKeys.insert(keyCode);
    }

    virtual void keyup(std::string keyCode) {
        activeKeys.erase(keyCode);
    }

    virtual void blur() {
        activeKeys.clear();
    }

    virtual void wheel(float deltaX, float deltaY) {
        if (deltaX != 0.0f)
            wheelDelta.x = deltaX;
        if (deltaY != 0.0f)
            wheelDelta.y = deltaY;
    }

    virtual ViewState getStateDerivative(ViewState x) = 0;

    void updateGamepadValues(float axes0, float axes1, float axes2, float axes3, float button6, float button7) {
        gamepadValues[0] = axes0;
        gamepadValues[1] = axes1;
        gamepadValues[2] = axes2;
        gamepadValues[3] = axes3;
        gamepadValues[4] = button6;
        gamepadValues[5] = button7;
    }
    
    void step(float dt) {
        wheelDelta *= exp(-1.0f / kWheelTau);

    }

};

class DefaultController : public BaseController {

    static constexpr float kSpeed = 0.5f;  // unit/s
    static constexpr float kOmegaR = 1.2f;  // rad/s
    static constexpr float kOmegaP = 1.0f;  // rad/s
    static constexpr float kOmegaY = 1.0f;  // rad/s
    static constexpr float kSpeedTau = 0.05f;  // s
    static constexpr float kOmegaTau = 0.05f;  // s
    static constexpr float kWheelSpeedX = 2.0f/20.0f;
    static constexpr float kWheelSpeedY = 2.5f/180.0f;

    static constexpr float kGamepadAxisThreshold = 0.1f;
    static constexpr float kGamepadSpeedX = kSpeed;
    static constexpr float kGamepadSpeedY = kSpeed;
    static constexpr float kGamepadAccelZ = 1.0f;
    static constexpr float kGamepadOmegaR = 2.0f*kOmegaR;
    static constexpr float kGamepadOmegaP = 2.0f*kOmegaP;
    static constexpr float kGamepadOmegaY = 2.0f*kOmegaY;

    static constexpr float kCollisionBoxRadius = 0.01f;  // unit
    static constexpr float kCollisionPenalty = 10.0f;  // unit^2 / (opac s^2)
    static constexpr float kCollisionVelocityDamping = 20.0f;  // unit / (opac s)
    static constexpr float kGravity = 0.0f;  // m/s^2

    float clipThreshold(float x) {
        return fabsf(x) > kGamepadAxisThreshold ? x : 0.0f;
    }

public:

    ViewState getStateDerivative(ViewState x) {
        using namespace glm;

        // x right, y forward, z up

        vec3 w = vec3(0);
        if (activeKeys.count("KeyA")) w.z += kOmegaY;
        if (activeKeys.count("KeyD")) w.z -= kOmegaY;
        if (activeKeys.count("KeyW")) w.x += kOmegaP;
        if (activeKeys.count("KeyS")) w.x -= kOmegaP;
        if (activeKeys.count("KeyQ")) w.y += kOmegaR;
        if (activeKeys.count("KeyE")) w.y -= kOmegaR;
        w.z -= kGamepadOmegaY * clipThreshold(gamepadValues[2]);
        w.x -= kGamepadOmegaP * clipThreshold(gamepadValues[3]);
        w.y += kGamepadOmegaR * (gamepadValues[4] - gamepadValues[5]);
        // w.z -= kGamepadOmegaY *  gamepadValues[0];
        // w.x -= kGamepadOmegaP *  gamepadValues[3];
        // w.y -= kGamepadOmegaR *  gamepadValues[2];

        vec3 v = vec3(0);
        if (activeKeys.count("ArrowUp")) v.y += kSpeed;
        if (activeKeys.count("ArrowDown")) v.y -= kSpeed;
        if (activeKeys.count("ArrowLeft")) v.x -= kSpeed;
        if (activeKeys.count("ArrowRight")) v.x += kSpeed;
        v.x += kSpeed * kWheelSpeedX * wheelDelta.x;
        v.y -= kSpeed * kWheelSpeedY * wheelDelta.y;
        v.x += kGamepadSpeedX * clipThreshold(gamepadValues[0]);
        v.y -= kGamepadSpeedY * clipThreshold(gamepadValues[1]);
        // v.y -= kGamepadSpeedY * gamepadValues[1];

        vec3 a = vec3(0);
        // a.z -= sceneScale * kGamepadAccelZ * gamepadValues[1];

        if (collisionDetector) {
            vec3 bound = vec3(kCollisionBoxRadius) * sceneScale;  // [m]
            vec3 grad = collisionDetector->query({ x.p-bound, x.p+bound });  // [opac m^2]
            grad /= 8.0f*bound.x*bound.y*bound.z;  // [opac / m]
            a -= (kCollisionPenalty*sceneScale*sceneScale) * grad;  // [m/s^2]
            a -= (kCollisionVelocityDamping*sceneScale) * length(grad) * x.v;
            a.z -= kGravity;
        }

        v *= sceneScale;  // unit/s -> m/s

        ViewState dx;
        if (kSpeedTau != 0.0) {
            dx.p = x.v;
            dx.v = a + (x.q * v - x.v) / kSpeedTau;  // also damps acceleration
        }
        else {
            dx.p = x.v + x.q * v;
            dx.v = a;
        }
        if (kOmegaTau != 0.0) {
            dx.q = 0.5f * x.q * quat(0.0f, x.w.x, x.w.y, x.w.z);
            dx.w = (w - x.w) / kOmegaTau;
        }
        else {
            x.w += w;
            dx.q = 0.5f * x.q * quat(0.0f, x.w.x, x.w.y, x.w.z);
            dx.w = vec3(0);
        }
        return dx;
    }

};


struct Controller {
    ViewState state;

    DefaultController controller;

    glm::mat4 b2w;

    static constexpr glm::mat4 b2c = glm::mat4(
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, -1, 0, 0,
        0, 0, 0, 1
    );

    static constexpr float kMaxStepDt = 1.0f;  // s
    static constexpr float kRenderNeededDelta = 1e-4f;

    static float matnorm(glm::mat4 a) {
        return sqrtf(glm::dot(a[0], a[0]) + glm::dot(a[1], a[1]) + glm::dot(a[2], a[2]) + glm::dot(a[3], a[3]));
    }

public:
    Controller() {
        using namespace glm;
        state.p = vec3(0);
        state.v = vec3(0);
        state.q = quat(1, 0, 0, 0);
        state.w = vec3(0);
        b2w = glm::mat4(0);
    }

    void setViewMatrix(emscripten::val m) {
        glm::mat4 w2c;
        for (int i = 0; i < 16; i++)
            w2c[i/4][i%4] = m[i].as<float>();
        glm::mat4 b2w = glm::inverse(w2c) * b2c;
        this->b2w = b2w;
        state.fromMatrix(b2w);
    }

    emscripten::val getViewMatrix() {
        glm::mat4 w2c = b2c * glm::inverse(b2w);
        return emscripten::val::global("Float32Array").new_(
            emscripten::typed_memory_view(16, (float*)&w2c));
    }

    void keydown(std::string keyCode) {
        controller.keydown(keyCode);
    }
    void keyup(std::string keyCode) {
        controller.keyup(keyCode);
    }
    void blur() {
        controller.blur();
    }

    void wheel(float deltaX, float deltaY) {
        controller.wheel(deltaX, deltaY);
    }
    void mousedown(float clientX, float clientY) {}
    void mouseup(float clientX, float clientY) {}
    void mousemove(float clientX, float clientY) {}
    void contextmenu(float clientX, float clientY) {}
    void touchstart() {}
    void touchmove() {}
    void touchend() {}
    void gamepadconnected() {
        controller.updateGamepadValues(0, 0, 0, 0, 0, 0);
    }
    void gamepaddisconnected() {
        controller.updateGamepadValues(0, 0, 0, 0, 0, 0);
    }

    void updateGamepadValues(emscripten::val gamepad) {
        controller.updateGamepadValues(
            gamepad["axes"][0].as<float>(),
            gamepad["axes"][1].as<float>(),
            gamepad["axes"][2].as<float>(),
            gamepad["axes"][3].as<float>(),
            gamepad["buttons"][6]["value"].as<float>(),
            gamepad["buttons"][7]["value"].as<float>()
        );
        
    }

    bool step(float dt) {
        if (!isfinite(dt)) return false;
        dt = glm::clamp(dt, 0.0f, kMaxStepDt);
        state.step(std::bind(&BaseController::getStateDerivative, &controller, std::placeholders::_1), dt);
        controller.step(dt);
        glm::mat4 b2w = state.toMatrix();
        float delta = matnorm(b2w - this->b2w);
        // printf("delta = %f\n", delta);
        if (delta > kRenderNeededDelta) {
            printf("delta = %f\n", delta);
            this->b2w = b2w;
            return true;
        }
        return false;
    }

    void setGaussians(emscripten::val base) {
        std::vector<float> means;
        copyToVector(base["means"], means);
        std::vector<float> scales;
        copyToVector(base["scales"], scales);
        std::vector<float> quats;
        copyToVector(base["quats"], quats);
        std::vector<float> opacs;
        copyToVector(base["opacities"], opacs);
        controller.setGaussians(
            (int)means.size()/3,
            (glm::vec3*)means.data(),
            (glm::vec3*)scales.data(),
            (glm::vec4*)quats.data(),
            opacs.data()
        );
    }
    
};
