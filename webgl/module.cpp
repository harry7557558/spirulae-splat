#include <emscripten/bind.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>

using namespace emscripten;

struct View {
    std::string type;
    int bitLength;
    int bitOffset;
    int quatBufferView;
    std::string key;
    int numElements;
};

struct Config {
    int componentLength;
    // std::vector<View> componentViews;
    val componentViews;
    int length;
};


template<typename T>
void copyToVector(const val &typedArray, std::vector<T> &vec) {
    vec = convertJSArrayToNumberVector<T>(typedArray);
    return;

    // https://github.com/emscripten-core/emscripten/issues/5519
    // doesn't work after memory growth
    unsigned int length = typedArray["length"].as<unsigned int>();
    val heap = val::module_property("HEAPU8");
    val memory = heap["buffer"];

    val memoryView = typedArray["constructor"].new_(memory, reinterpret_cast<uintptr_t>(vec.data()), length);
    
    memoryView.call<void>("set", typedArray);
}


struct Primitive {
    std::vector<float> output;
};


val unpackComponents(std::string key, const Config& config, const val& packedData, const val& buffers) {
    int componentLength = config.componentLength;
    // std::vector<View> componentViews = config.componentViews;
    std::vector<View> componentViews = vecFromJSArray<View>(config.componentViews);

    if (packedData.isNull() || packedData.isUndefined()) {
        return val::null();
    }

    int length = std::min(config.length, packedData["byteLength"].as<int>() / componentLength);
    std::vector<std::vector<uint32_t>> components;

    for (View& view : componentViews) {
        assert(view.type.find("quat") == 0 && "only quat type is supported");
        if (view.type.find("quat") == 0) {
            view.numElements = view.type.size() == 4 ? 1 :
                std::max(std::stoi(view.type.substr(4)), 1);
            components.push_back(std::vector<uint32_t>(length * view.numElements));
        }
    }

    val packedDataUint8 = val::global("Uint8Array").new_(packedData);
    std::vector<uint8_t> packedDataV;
    copyToVector(packedDataUint8, packedDataV);

    for (int i = 0; i < length; ++i) {

        for (size_t j = 0; j < componentViews.size(); ++j) {
            const View& view = componentViews[j];
            int bitLength = view.bitLength;
            int bitOffset = view.bitOffset;

            // assume quat type

            int numElements = view.numElements;
            int bitsPerElement = bitLength / numElements;

            for (int k = 0; k < numElements; ++k) {
                uint32_t value = 0;
                int bitStart = bitOffset + k * bitsPerElement;
                int byteStart = bitStart / 8;
                int bitPos = bitStart % 8;

                for (int b = 0; b < bitsPerElement; ++b) {
                    int byteIdx = byteStart + (bitPos + b) / 8;
                    int bitInByte = (bitPos + b) % 8;
                    if (packedDataV[i * componentLength + byteIdx] & (1 << bitInByte)) {
                        value |= (1 << b);
                    }
                }
                components[j][i * numElements + k] = value;
            }
        }
    }

    val object = val::object();
    for (size_t j = 0; j < componentViews.size(); ++j) {
        const auto& view = componentViews[j];
        bool hasQuat = view.type.find("quat") == 0;
        if (hasQuat) {
            int bl = buffers[view.quatBufferView]["byteLength"].as<int>();
            if (bl < 4 || (bl & (bl - 1)) != 0) {
                hasQuat = false;
            }
        }
        if (hasQuat) {
            int numel = components[j].size();
            std::vector<float> component(numel);
            val buffer = val::global("Float32Array").new_(buffers[view.quatBufferView]);
            std::vector<float> bufferV;
            copyToVector(buffer, bufferV);
            for (int i = 0; i < numel; ++i) {
                component[i] = bufferV[components[j][i]];
            }
            object.set(view.key, val::global("Float32Array").new_(typed_memory_view(component.size(), component.data())));
        } else {
            object.set(view.key, val::global("Uint32Array").new_(typed_memory_view(components[j].size(), components[j].data())));
        }
    }

    return object;
}



uint32_t floatToHalf(float floatValue) {
    union FloatInt {
        float f;
        uint32_t i;
    };

    FloatInt floatInt;
    floatInt.f = floatValue;
    uint32_t f = floatInt.i;

    uint32_t sign = (f >> 31) & 0x0001;
    uint32_t exp = (f >> 23) & 0x00ff;
    uint32_t frac = f & 0x007fffff;

    uint32_t newExp;
    if (exp == 0) {
        newExp = 0;
    } else if (exp < 113) {
        newExp = 0;
        frac |= 0x00800000; // Add implicit leading 1
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

uint32_t packHalf2x16(float x, float y) {
    return (floatToHalf(x) | (floatToHalf(y) << 16));
}


val packHarmonicTexture(std::string key, int dim, int width, int height, const val& features) {

    std::vector<uint32_t> data(width * height * 4);
    int pixelPerVert = (dim+1)/2;
    int vertexCount = int(features["length"].as<int>()/(3*dim));

    std::vector<float> featuresV;
    copyToVector(features, featuresV);

    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < dim; j++) {
            int iTex = 4 * pixelPerVert * i + 2 * j;
            int iBuffer = 3 * dim * i + 3 * j;
            data[iTex+0] = packHalf2x16(featuresV[iBuffer], featuresV[iBuffer+1]);
            data[iTex+1] = packHalf2x16(featuresV[iBuffer+2], 1.0);
        }
    }

    return val::global("Uint32Array").new_(typed_memory_view(data.size(), data.data()));
}



int main() {
    printf("WASM Module Initialized.\n");
    return 0;
}

EMSCRIPTEN_BINDINGS(module) {
    register_vector<View>("VectorView");
    register_vector<val>("VectorVal");

    value_object<View>("View")
        .field("type", &View::type)
        .field("bitLength", &View::bitLength)
        .field("bitOffset", &View::bitOffset)
        .field("quatBufferView", &View::quatBufferView)
        .field("key", &View::key);

    value_object<Config>("Config")
        .field("componentLength", &Config::componentLength)
        .field("componentViews", &Config::componentViews)
        .field("length", &Config::length);

    function("unpackComponents", &unpackComponents);
    function("packHarmonicTexture", &packHarmonicTexture);
}

