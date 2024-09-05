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


py::dict unpack_components(
    const py::dict& config,
    const py::buffer& packed_data_,
    const std::vector<py::buffer>& buffers
) {
    int component_length = config["componentLength"].cast<int>();
    auto component_views = config["componentViews"].cast<std::vector<py::dict>>();
    py::buffer_info packed_data_info = packed_data_.request();
    int length = std::min(config["length"].cast<int>(), static_cast<int>(packed_data_info.size / component_length));
    const uint8_t* packed_data = (uint8_t*)packed_data_info.ptr;

    std::vector<std::vector<uint32_t>> components;
    for (const auto& view : component_views) {
        std::string data_type = view["type"].cast<std::string>();
        if (data_type.find("quat") == 0) {
            int num_elements = (data_type == "quat") ? 1 : std::max(std::stoi(data_type.substr(4)), 1);
            components.emplace_back(length * num_elements);
        }
    }

    for (int i = 0; i < length; ++i) {
        for (size_t j = 0; j < component_views.size(); ++j) {
            const auto& view = component_views[j];
            int bit_length = view["bitLength"].cast<int>();
            int bit_offset = view["bitOffset"].cast<int>();
            std::string data_type = view["type"].cast<std::string>();

            if (data_type.find("quat") == 0) {
                int num_elements = (data_type == "quat") ? 1 : std::max(std::stoi(data_type.substr(4)), 1);
                int bits_per_element = bit_length / num_elements;

                for (int k = 0; k < num_elements; ++k) {
                    uint32_t value = 0;
                    int bit_start = bit_offset + k * bits_per_element;
                    int byte_start = bit_start / 8;
                    int bit_pos = bit_start % 8;

                    for (int b = 0; b < bits_per_element; ++b) {
                        int byte_idx = byte_start + (bit_pos + b) / 8;
                        int bit_in_byte = (bit_pos + b) % 8;
                        if (((uint8_t)packed_data[i * component_length + byte_idx] & (1 << bit_in_byte)) != 0) {
                            value |= (1U << b);
                        }
                    }
                    components[j][i * num_elements + k] = value;
                }
            }
        }
    }

    py::dict result;
    for (size_t j = 0; j < component_views.size(); ++j) {
        const auto& view = component_views[j];
        std::string data_type = view["type"].cast<std::string>();
        bool has_quat = data_type.find("quat") == 0;
        
        if (has_quat && view.contains("quatBufferView")) {
            int buffer_view = view["quatBufferView"].cast<int>();
            auto buffer = buffers[buffer_view].request();
            int bl = buffer.size;
            if (bl < 4 || (bl & (bl - 1)) != 0) {
                has_quat = false;
            }
        }

        if (has_quat) {
            int buffer_view = view["quatBufferView"].cast<int>();
            auto buffer = buffers[buffer_view].request();
            int numel = components[j].size();
            py::array_t<float> component(numel);
            auto component_ptr = component.mutable_data();
            auto buffer_ptr = (const float*)(buffer.ptr);
            for (int i = 0; i < numel; ++i) {
                component_ptr[i] = buffer_ptr[components[j][i]];
            }
            result[view["key"].cast<std::string>().c_str()] = component;
        } else {
            py::array_t<uint32_t> component(components[j].size());
            auto component_ptr = component.mutable_data();
            std::copy(components[j].begin(), components[j].end(), component_ptr);
            result[view["key"].cast<std::string>().c_str()] = component;
        }
    }

    return result;
}


PYBIND11_MODULE(pack_components, m) {
    m.def("pack_components", &pack_components, "Pack components function");
    m.def("unpack_components", &unpack_components, "Unpack components function");
}
