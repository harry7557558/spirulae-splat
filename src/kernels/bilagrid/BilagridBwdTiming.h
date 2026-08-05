#pragma once
// Async per-dispatch GPU timer for the bilateral-grid backward (Phase 2 of the
// bilagrid backward-selection plan; see BilagridBackwardSelection.md). It
// brackets each backward dispatch with a start/end timing event and reads the
// elapsed GPU time back ONE ITERATION BEHIND, so the hot path never stalls on
// a query wait. The harvested measurements feed the BilagridBwdSelector
// (Phase 5); in Phase 2 they are only logged.
//
// Backend-portable: uses only the backend:: event API. On CUDA an event record
// is a cheap cudaEventRecord; on Vulkan it flushes the stream so the two
// timestamps bracket exactly the bilagrid dispatch(es) in between.

#include "kernels/bilagrid/BilagridBwdSelector.h"          // ContextKey
#include "backend/api/BackendRuntime.h"   // backend::Event, Stream, event_*

#include <string>
#include <vector>

namespace spirula {
namespace bilagrid {

struct BwdMeasurement {
    ContextKey key;
    int arm;      // arm index that was dispatched (-1 = pre-selector baseline)
    double ms;    // measured GPU time
};

// Timed channels of the backward hook. RGB carries one of affine/ppisp/
// loglinear (distinguished by ContextKey.family); depth/normal are separate.
enum BwdChannel {
    BWD_RGB = 0,
    BWD_DEPTH = 1,
    BWD_NORMAL = 2,
    BWD_NUM_CHANNELS = 3,
};

inline int family_id(const std::string& type) {
    if (type == "affine") return 0;
    if (type == "ppisp") return 1;
    if (type == "loglinear") return 2;
    if (type == "depth") return 3;
    if (type == "normal") return 4;
    return -1;
}

inline const char* family_name(int id) {
    switch (id) {
        case 0: return "affine";
        case 1: return "ppisp";
        case 2: return "loglinear";
        case 3: return "depth";
        case 4: return "normal";
        default: return "?";
    }
}

// Fixed-slot async timer (one slot per BwdChannel). Timing events are created
// once and reused across iterations (on Vulkan each holds one query-pool slot).
class BwdTimingRing {
public:
    explicit BwdTimingRing(int num_slots = BWD_NUM_CHANNELS)
        : slots_(num_slots) {}

    ~BwdTimingRing() {
        for (Slot& s : slots_) {
            if (s.start) backend::event_destroy(s.start);
            if (s.end) backend::event_destroy(s.end);
        }
    }

    // Complete all pending measurements from PRIOR iterations and append them
    // to `out`. Call once at the top of the backward hook, BEFORE any begin()
    // this iteration -- every pending slot is then >= 1 iteration old, so the
    // synchronize is effectively free (the GPU finished it long ago).
    void harvest(std::vector<BwdMeasurement>& out) {
        for (Slot& s : slots_) {
            if (!s.pending) continue;
            backend::event_synchronize(s.end);
            double ms = (double)backend::event_elapsed_ms(s.start, s.end);
            out.push_back({s.key, s.arm, ms});
            s.pending = false;
        }
    }

    // Record the start event for `slot` (call immediately before the dispatch).
    void begin(int slot, const ContextKey& key, int arm,
               backend::Stream stream) {
        Slot& s = slots_[slot];
        if (!s.start) s.start = backend::event_create(/*enable_timing=*/true);
        if (!s.end) s.end = backend::event_create(/*enable_timing=*/true);
        s.key = key;
        s.arm = arm;
        backend::event_record(s.start, stream);
    }

    // Record the end event for `slot` (call immediately after the dispatch).
    void end(int slot, backend::Stream stream) {
        Slot& s = slots_[slot];
        backend::event_record(s.end, stream);
        s.pending = true;
    }

    int num_slots() const { return (int)slots_.size(); }

private:
    struct Slot {
        backend::Event* start = nullptr;
        backend::Event* end = nullptr;
        ContextKey key;
        int arm = -1;
        bool pending = false;
    };
    std::vector<Slot> slots_;
};

}  // namespace bilagrid
}  // namespace spirula
