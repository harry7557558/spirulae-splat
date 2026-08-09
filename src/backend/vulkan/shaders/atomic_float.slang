// Portable float atomic add (Vulkan backend). VK_EXT_shader_atomic_float is
// not baseline (and absent under MoltenVK), so the default atomic_add_f32 is
// a compare-exchange loop over the value's u32 bit image — same
// nondeterministic accumulation order as CUDA atomicAdd, so parity stays
// tolerance-based wherever this is used. The loop's exit is a successful
// CAS (well-defined; not a UB-terminated loop — see README loop rules).
//
// Under -DSS_NATIVE_F32_ATOMIC_ADD the emulation is replaced by the
// native OpAtomicFAddEXT (SPV_EXT_shader_atomic_float_add): spirv_tool.cpp
// compiles every entry whose include closure reaches this file a second
// time with the define, and the pipeline layer picks the native blob when
// the device enabled shaderBufferFloat32AtomicAdd — so the choice is made
// at pipeline creation, never as a runtime branch in the shader.

// Float atomic max mirroring Common.cuh's atomicMax(float*, float): signed
// int max for non-negative v, unsigned min for negative v (IEEE order in
// both halves). Single atomic — no CAS loop needed.
[ForceInline]
void atomic_max_f32(float* addr, float v) {
    if (asint(v) >= 0) {
        int* w = (int*)addr;
        InterlockedMax(w[0], asint(v));
    } else {
        uint* w = (uint*)addr;
        InterlockedMin(w[0], asuint(v));
    }
}

#ifdef SS_NATIVE_F32_ATOMIC_ADD

[ForceInline]
void atomic_add_f32(float* addr, float v) {
    if (v == 0.0f) return;
    InterlockedAdd(addr[0], v);
}

#else

[ForceInline]
void atomic_add_f32(float* addr, float v) {
    if (v == 0.0f) return;
    uint* w = (uint*)addr;
    uint expected = w[0];
    for (;;) {
        uint desired = asuint(asfloat(expected) + v);
        uint prev;
        InterlockedCompareExchange(w[0], expected, desired, prev);
        if (prev == expected) break;
        expected = prev;
    }
}

#endif
