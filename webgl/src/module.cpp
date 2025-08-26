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


val sortByDepth(int vertexCount, const val& means_, bool is_fisheye, float n0, float n1, float n2) {
    static const int QN = 16384;

    std::vector<int> sizeList(vertexCount);

    std::vector<float> means;
    copyToVector(means_, means);

    for (int i = 0; i < vertexCount; i++) {
        float dx = means[3*i+0]-n0, dy = means[3*i+1]-n1, dz = means[3*i+2]-n2;
        // float z = sqrtf(dx*dx+dy*dy+dz*dz);
        float z = (dx*dx+dy*dy+dz*dz);
        sizeList[i] = int(((float)QN-0.1f)/(z+1.0f));
    }

    std::vector<int> counts0(QN, 0);
    for (int i = 0; i < vertexCount; i++) {
        counts0[sizeList[i]]++;
    }

    std::vector<int> starts0(QN);
    starts0[0] = 0;
    for (int i = 1; i < QN; i++)
        starts0[i] = starts0[i-1] + counts0[i-1];

    std::vector<uint32_t> depthIndex(vertexCount, -1);
    for (int i = 0; i < vertexCount; i++)
        depthIndex[starts0[sizeList[i]]++] = i;

    return val::global("Uint32Array").new_(typed_memory_view(depthIndex.size(), depthIndex.data()));
}


val preparePPSTiles(int vertexCount, int innerWidth, int innerHeight, const val& projData_, const val& depthIndex_) {
    const int TILE_SIZE = 12;
    int numTilesX = std::ceil(innerWidth / TILE_SIZE);
    int numTilesY = std::ceil(innerHeight / TILE_SIZE);
    int numTiles = numTilesX * numTilesY;

    std::vector<uint32_t> projData;
    copyToVector(projData_, projData);
    std::vector<uint32_t> depthIndex;
    copyToVector(depthIndex_, depthIndex);

    // Create 2D vector to store tile data
    std::vector<std::vector<int>> tiles(numTiles);
    
    // Process vertices and assign them to tiles
    for (int i0 = 0; i0 < vertexCount; i0++) {
        int i = (int)depthIndex[i0];
        int idx = (int)projData[12*i];
        
        if (idx == -1) {
            continue;
        }
        
        int x0 = (int)projData[12*i+8];
        int y0 = (int)projData[12*i+9];
        int x1 = (int)projData[12*i+10];
        int y1 = (int)projData[12*i+11];
        
        for (int x = x0; x < x1; x++)
            for (int y = y0; y < y1; y++) {
                tiles[y * numTilesX + x].push_back(idx);
            }
    }
    
    // Create PSA for tile indices
    std::vector<uint32_t> numTilesPSA(4 * numTiles);
    numTilesPSA[0] = 0;
    for (int i = 0; i < numTiles; i++) {
        numTilesPSA[4*i+1] = numTilesPSA[4*i] + tiles[i].size();
        numTilesPSA[4*(i+1)] = numTilesPSA[4*i+1];
    }
    
    int numIntersects = numTilesPSA[4*numTiles-3];
    int intTexWidth = 2048;
    int intTexHeight = std::ceil(numIntersects * 2.0 / intTexWidth);
    
    // Create the intersects array
    std::vector<uint32_t> intersects(4 * intTexWidth * intTexHeight);
    
    // Copy data from projData to intersects
    for (int i = 0; i < numTiles; i++) {
        for (size_t j = 0; j < tiles[i].size(); j++) {
            int sid = tiles[i][j];
            int iid = numTilesPSA[4*i] + j;
            for (int k = 0; k < 12; k++)
                intersects[8*iid+k] = projData[12*sid+k];
        }
    }
    
    // Create and return a result object with all the computed data
    val result = val::object();
    result.set("numTilesX", numTilesX);
    result.set("numTilesY", numTilesY);
    result.set("numTiles", numTiles);
    result.set("numIntersects", numIntersects);
    result.set("intTexWidth", intTexWidth);
    result.set("intTexHeight", intTexHeight);
    result.set("numTilesPSA", val::global("Uint32Array").new_(typed_memory_view(numTilesPSA.size(), numTilesPSA.data())));
    result.set("intersects", val::global("Uint32Array").new_(typed_memory_view(intersects.size(), intersects.data())));
    
    return result;
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

    function("sortByDepth", &sortByDepth);

    function("preparePPSTiles", &preparePPSTiles);
}

