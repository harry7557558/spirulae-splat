#pragma once

// A distorted camera into a monocular geometry model and back: the frame is
// resampled into pinhole faces on the way in and the prediction gathered into
// the source camera's own frame on the way out, planned once per camera.
//
// The faces OVERLAP and `gather` cross-fades them, because each is a separate
// call into a network that re-guesses its depth scale and its horizon every
// time. src/metric3d/README.md, "The camera side", is the rest of it.

#include "core/CameraModel.h"

#include <cstdint>
#include <string>
#include <vector>

namespace app {

// [sh, sw, channels] bytes -> [dh, dw, 3] floats in 0..1, what sampleFace
// takes. Area average, not bilinear: at a 2x downscale bilinear reads two of
// every three source pixels and aliases the structure the normals are made of.
std::vector<float> resize_area(const uint8_t* src, int sw, int sh, int channels,
                               int dw, int dh);

// One camera, as data/DatasetParser.h describes it.
struct GeometryCamera {
    int   model = 0;        // CameraModelType
    int   distortion = 0;   // CameraDistortionType
    int   width = 0, height = 0;
    float fx = 0, fy = 0, cx = 0, cy = 0;
    float dist[kCameraDistortionParams] = {};
    // A camera whose lens no tier represents exactly; the parser fitted it and
    // recorded the true projection. -1 when the tier is exact.
    int   source_model = -1;
    float source_params[16] = {};
};

class GeometryWarp {
public:
    // A split face is sized off cam.width x cam.height, NOT out_w x out_h, and
    // only then capped at `max_face`; sides land on a multiple of `patch`.
    // Throws std::runtime_error when the camera leaves no usable face.
    void plan(const GeometryCamera& cam, int out_w, int out_h, bool split, int patch,
              int max_face);

    bool split() const { return faces_ > 1; }
    int  faces() const { return faces_; }
    int  faceWidth() const { return fw_; }
    int  faceHeight() const { return fh_; }
    int  outWidth() const { return out_w_; }
    int  outHeight() const { return out_h_; }
    // A split keeps the wide capture's ray semantics; one face keeps z.
    bool defaultRayDepth() const { return split(); }

    // [faces(), 3, 3], or null for one face. Rows 0 and 1 span the face and
    // carry its half-extent in their LENGTH; row 2 is the unit optical axis.
    const double* faceAxes() const { return axes_; }

    // The pinhole a face IS, in its own pixels. Metric3D's depth is canonical
    // to a focal of 1000, so metres are the prediction times faceFocal() over
    // 1000; MoGe needs all four to resolve its point map's z shift.
    double faceFocal() const { return face_focal_; }
    double faceFocalY() const { return face_focal_y_; }
    double faceCx() const { return face_cx_; }
    double faceCy() const { return face_cy_; }

    // The resolution to hand `sampleFace`: full unless `max_face` capped the
    // faces, and NOT the output resolution, which is capped independently.
    int sampleWidth() const { return sw_; }
    int sampleHeight() const { return sh_; }

    // [sampleHeight(), sampleWidth(), 3] in 0..1 -> [faceHeight(),
    // faceWidth(), 3]. What the camera cannot see is filled mid-grey; a black
    // border reads as an edge to the network.
    void sampleFace(int k, const float* src, std::vector<float>& dst) const;

    // `depth[k]` is [fh*fw] canonical depth and `normal[k]` is [fh*fw*3] unit
    // normals in face k's frame; either list may be empty to skip it. Pixels
    // no face covers come back 0, the trainer's "no ground truth here".
    void gather(const std::vector<std::vector<float>>& depth,
                const std::vector<std::vector<float>>& normal, bool ray_depth,
                std::vector<float>* out_depth, std::vector<float>* out_normal) const;

private:
    // One face's claim on one output pixel.
    struct Contrib {
        int32_t face;
        float   px, py;    // where in that face, in its own pixel coordinates
        float   w;         // cross-fade weight, zero at the face's outer edge
        float   sz, sr;    // face z-depth -> this camera's z / distance
    };

    std::vector<double> alignFaces(const std::vector<std::vector<float>>& depth) const;

    const double* axes_ = nullptr;
    int faces_ = 1, fw_ = 0, fh_ = 0, out_w_ = 0, out_h_ = 0, sw_ = 0, sh_ = 0;
    // [fh*fw*2] coordinates in the sampleWidth() grid, NaN where the camera
    // has no ray at all, and [fh*fw] of how much of the sample is in frame.
    std::vector<std::vector<float>> to_src_;
    std::vector<std::vector<float>> face_valid_;
    // Which faces reach each output pixel. Row offsets rather than a fixed
    // slot per face: only the overlap names more than one.
    std::vector<int64_t> contrib_off_;
    std::vector<Contrib> contrib_;
    double face_focal_ = 0.0, face_focal_y_ = 0.0, face_cx_ = 0.0, face_cy_ = 0.0;
};

}  // namespace app
