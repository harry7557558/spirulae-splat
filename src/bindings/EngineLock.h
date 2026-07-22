#pragma once

// A Python context manager over a `std::mutex` owned by C++.
//
// Both the web viewer (bind_viewer.cpp) and the trainer session
// (bind_trainer.cpp) hand Python the same engine mutex that their render
// worker threads take. Acquiring it must release the GIL: the worker holds
// the mutex *without* the GIL, so a Python thread that blocked on the mutex
// while holding the GIL would deadlock the pair.
//
// Bound once, as `_C._EngineLock` -- pybind registers a C++ type at most once
// per module, so `bind_engine_lock()` is called from exactly one place.

#include <pybind11/pybind11.h>

#include <mutex>

namespace ssplat_bindings {

class EngineLock {
public:
    explicit EngineLock(std::mutex* m) : _m(m) {}

    void enter() {
        pybind11::gil_scoped_release release;
        _m->lock();
        _held = true;
    }

    void exit(const pybind11::object&, const pybind11::object&,
              const pybind11::object&) {
        if (!_held) return;
        _held = false;
        _m->unlock();
    }

private:
    std::mutex* _m;
    bool _held = false;
};

inline void bind_engine_lock(pybind11::module_& m) {
    pybind11::class_<EngineLock>(m, "_EngineLock")
        .def("__enter__", &EngineLock::enter)
        .def("__exit__", &EngineLock::exit);
}

}  // namespace ssplat_bindings
