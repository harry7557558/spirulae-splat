// Python bindings for the native dataset parsers (data/DatasetParser.h).
//
// These expose the SAME parsers the CLI trainer, the GUI and the WASM viewer
// already use, so Python stops carrying a second implementation of COLMAP /
// Nerfstudio / Metashape reading. See docs/restructure-proposal.md §4.1.
//
// Everything is returned as plain numpy arrays / lists; nothing here depends
// on torch, so this file also compiles in a torch-free build (it is simply not
// linked into one, since the extension module is CUDA+torch only).
//
// The Python-side adapter is spirulae_splat/modules/native_dataparser.py, and
// tests/python/test_dataparser_parity.py is the gate that must stay green
// until modules/dataparser.py is deleted.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "data/DatasetParser.h"

namespace py = pybind11;

namespace {

// Wrap a std::vector as a 2-D numpy array WITHOUT copying the vector twice:
// pybind copies once into the returned array, which is what we want (the
// ParsedDataset dies at the end of the call).
template <typename T>
py::array_t<T> as_array(const std::vector<T>& v, std::vector<py::ssize_t> shape) {
    py::ssize_t total = 1;
    for (auto s : shape) total *= s;
    if (total != (py::ssize_t)v.size()) {
        throw std::runtime_error("bind_data: shape does not match vector size");
    }
    py::array_t<T> out(shape);
    if (!v.empty()) {
        std::memcpy(out.mutable_data(), v.data(), v.size() * sizeof(T));
    }
    return out;
}

template <typename T>
py::array_t<T> as_array(const std::vector<T>& v) {
    return as_array(v, {(py::ssize_t)v.size()});
}

}  // namespace

void bind_data(py::module_& m) {
    // -----------------------------------------------------------------
    // DatasetParserConfig
    // -----------------------------------------------------------------
    // Field names match SpirulaeSplatDataParserConfig in
    // modules/dataparser.py wherever the two overlap, so the adapter is a
    // straight copy rather than a translation table.
    py::class_<DatasetParserConfig>(m, "DatasetParserConfig")
        .def(py::init<>())
        .def_readwrite("recon_dir", &DatasetParserConfig::recon_dir)
        .def_readwrite("image_dir", &DatasetParserConfig::image_dir)
        .def_readwrite("mask_dir", &DatasetParserConfig::mask_dir)
        .def_readwrite("depth_dir", &DatasetParserConfig::depth_dir)
        .def_readwrite("normal_dir", &DatasetParserConfig::normal_dir)
        .def_readwrite("require_image_files", &DatasetParserConfig::require_image_files)
        .def_readwrite("validation_fraction", &DatasetParserConfig::validation_fraction)
        .def_readwrite("eval_mode", &DatasetParserConfig::eval_mode)
        .def_readwrite("train_split_fraction", &DatasetParserConfig::train_split_fraction)
        .def_readwrite("eval_interval", &DatasetParserConfig::eval_interval)
        .def_readwrite("outlier_threshold", &DatasetParserConfig::outlier_threshold)
        .def_readwrite("rescale_camera_to_fit", &DatasetParserConfig::rescale_camera_to_fit)
        .def_readwrite("downscale_rounding_mode", &DatasetParserConfig::downscale_rounding_mode)
        .def_readwrite("metashape_xml", &DatasetParserConfig::metashape_xml)
        .def_readwrite("metashape_ply", &DatasetParserConfig::metashape_ply)
        .def_readwrite("metashape_psx", &DatasetParserConfig::metashape_psx)
        .def_readwrite("metashape_component", &DatasetParserConfig::metashape_component);

    // -----------------------------------------------------------------
    // ParsedDataset
    // -----------------------------------------------------------------
    py::class_<ParsedDataset>(m, "ParsedDataset")
        .def_readonly("num_cameras", &ParsedDataset::num_cameras)
        .def_readonly("image_filenames", &ParsedDataset::image_filenames)
        .def_readonly("mask_filenames", &ParsedDataset::mask_filenames)
        .def_readonly("depth_filenames", &ParsedDataset::depth_filenames)
        .def_readonly("normal_filenames", &ParsedDataset::normal_filenames)
        .def_readonly("train_frame_scale", &ParsedDataset::train_frame_scale)
        .def_property_readonly("camera_models", [](const ParsedDataset& d) {
            return as_array(d.camera_models);
        })
        .def_property_readonly("widths", [](const ParsedDataset& d) {
            return as_array(d.widths);
        })
        .def_property_readonly("heights", [](const ParsedDataset& d) {
            return as_array(d.heights);
        })
        // [N, 3, 4] camera-to-world, OpenGL/nerfstudio convention.
        .def_property_readonly("c2w", [](const ParsedDataset& d) {
            return as_array(d.c2w, {d.num_cameras, 3, 4});
        })
        // [N, 4] fx, fy, cx, cy
        .def_property_readonly("intrins", [](const ParsedDataset& d) {
            return as_array(d.intrins, {d.num_cameras, 4});
        })
        // [N, 10] k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2
        .def_property_readonly("dist_coeffs", [](const ParsedDataset& d) {
            return as_array(d.dist_coeffs, {d.num_cameras, 10});
        })
        .def_property_readonly("train_indices", [](const ParsedDataset& d) {
            return as_array(d.train_indices);
        })
        .def_property_readonly("val_indices", [](const ParsedDataset& d) {
            return as_array(d.val_indices);
        })
        .def_property_readonly("points_xyz", [](const ParsedDataset& d) {
            return as_array(d.points.xyz, {d.points.num(), 3});
        })
        .def_property_readonly("points_rgb", [](const ParsedDataset& d) {
            return as_array(d.points.rgb, {d.points.num(), 3});
        })
        // Row-major 4x4; identity when train_frame_scale == 1.
        .def_property_readonly("train_to_normalized", [](const ParsedDataset& d) {
            py::array_t<float> out({(py::ssize_t)4, (py::ssize_t)4});
            std::memcpy(out.mutable_data(), d.train_to_normalized.data(),
                        16 * sizeof(float));
            return out;
        })
        .def("__repr__", [](const ParsedDataset& d) {
            return "<ParsedDataset " + std::to_string(d.num_cameras) +
                   " cameras, " + std::to_string(d.points.num()) + " points>";
        });

    // -----------------------------------------------------------------
    // Parsers
    // -----------------------------------------------------------------
    m.def("parse_dataset",
          [](const std::string& dataset_dir, const DatasetParserConfig& cfg,
             const std::string& format) {
              py::gil_scoped_release release;  // pure C++, no Python touched
              return parse_dataset(dataset_dir, cfg, format);
          },
          py::arg("dataset_dir"), py::arg("config") = DatasetParserConfig(),
          py::arg("format") = "",
          "Parse a dataset directory. format: '' (auto-detect, tries "
          "nerfstudio then colmap then metashape), 'colmap', 'nerfstudio' or "
          "'metashape'.");

    m.def("parse_colmap_dataset",
          [](const std::string& d, const DatasetParserConfig& c) {
              py::gil_scoped_release release;
              return parse_colmap_dataset(d, c);
          },
          py::arg("dataset_dir"), py::arg("config") = DatasetParserConfig());

    m.def("parse_nerfstudio_dataset",
          [](const std::string& d, const DatasetParserConfig& c) {
              py::gil_scoped_release release;
              return parse_nerfstudio_dataset(d, c);
          },
          py::arg("dataset_dir"), py::arg("config") = DatasetParserConfig());

    m.def("parse_metashape_dataset",
          [](const std::string& d, const DatasetParserConfig& c) {
              py::gil_scoped_release release;
              return parse_metashape_dataset(d, c);
          },
          py::arg("dataset_dir"), py::arg("config") = DatasetParserConfig());
}
