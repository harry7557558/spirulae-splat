#include <vector>
#include <string>
#include <cassert>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

std::vector<uint8_t> pack_components(
    const py::dict& config,
    const std::vector<py::array_t<
        int32_t, py::array::c_style | py::array::forcecast
    >>& components
) {
    int length = config["length"].cast<int>();
    int component_length = config["componentLength"].cast<int>();
    auto component_views = config["componentViews"].cast<std::vector<py::dict>>();
    
    std::vector<uint8_t> packed_data(length * component_length, 0);
    
    for (int i = 0; i < length; ++i) {
        for (size_t j = 0; j < component_views.size(); ++j) {
            auto data = components[j].unchecked<2>();
            auto view = component_views[j];
            int bit_length = view["bitLength"].cast<int>();
            int bit_offset = view["bitOffset"].cast<int>();
            std::string data_type = view["type"].cast<std::string>();
            
            if (data_type.find("quat") == 0) {
                int num_elements = (data_type.size() > 4) ? std::stoi(data_type.substr(4)) : 1;
                if (num_elements == 0) continue;
                assert(bit_length % num_elements == 0);
                int bits_per_element = bit_length / num_elements;
                for (int k = 0; k < num_elements; ++k) {
                    uint32_t value = data(i, k);
                    assert(value >= 0 && value < (1U << bits_per_element));
                    int bit_start = bit_offset + k * bits_per_element;
                    int byte_start = bit_start / 8;
                    int bit_pos = bit_start % 8;
                    for (int b = 0; b < bits_per_element; ++b) {
                        int byte_idx = byte_start + (bit_pos + b) / 8;
                        int bit_in_byte = (bit_pos + b) % 8;
                        if (value & (1 << b)) {
                            packed_data[i*component_length + byte_idx] |= (1 << bit_in_byte);
                        }
                    }
                }
            }
            else {
                assert(false && "only quat type is supported");
            }
        }
    }
    
    return packed_data;
}

PYBIND11_MODULE(pack_components, m) {
    m.def("pack_components", &pack_components, "Pack components function");
}
