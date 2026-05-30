// Single-definition storage for the EngineState singleton.

#include "EngineState.h"

EngineState& engine() {
    static EngineState s;
    return s;
}
