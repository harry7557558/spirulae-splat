// Single-definition storage for the EngineState singleton.

#include "Engine.h"
#include "EngineState.h"

EngineState& engine() {
    static EngineState s;
    return s;
}

void engine_reset() {
    // All EngineState members are non-owning views into the global pool, so
    // overwriting with a default-constructed instance just drops the cached
    // pointers/flags; the freeAllDeviceMemory call below actually releases
    // the underlying device memory. Pool buffer keys will be re-acquired
    // (and resized) by the next scene's engine_* init calls.
    engine() = EngineState{};
    freeAllDeviceMemory();
}
