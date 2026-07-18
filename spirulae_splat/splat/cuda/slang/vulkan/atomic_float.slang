// Portable float atomic add (Vulkan backend). VK_EXT_shader_atomic_float is
// not baseline (and absent under MoltenVK), so atomic float adds run as a
// compare-exchange loop over the value's u32 bit image — same
// nondeterministic accumulation order as CUDA atomicAdd, so parity stays
// tolerance-based wherever this is used. The loop's exit is a successful
// CAS (well-defined; not a UB-terminated loop — see README loop rules).

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
