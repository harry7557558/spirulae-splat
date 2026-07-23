// Adaptive byte-granular buffer access for the Vulkan compute entries.
//
// The baseline reads/writes uint8 buffer elements through u32 WORDS (shift +
// mask) because shaderInt8 is not in the device baseline; word reads past the
// last byte stay in bounds because device allocations are 16-byte rounded
// (VulkanRuntime). Byte writes must not clobber the neighboring bytes of the
// word, so the baseline store is an InterlockedAnd+InterlockedOr pair.
//
// Under -DSSPLAT_NATIVE_INT8 (the ".int8" blob variants, picked by
// VulkanPipelines.cpp on devices with shaderInt8 + storageBuffer8BitAccess)
// these lower to native byte loads/stores — notably turning the two-atomic
// store into a plain store.
//
// Buffers are passed as uint* (the word view) on BOTH paths so kernel
// signatures and push layouts never change between variants.

#ifndef SSPLAT_INT8_COMPAT
#define SSPLAT_INT8_COMPAT

#include "int64_compat.slang"

// Each helper has a uint overload alongside the gindex_t one so 32-bit call
// sites keep 32-bit index arithmetic — a lone gindex_t signature would
// promote them to int64 and drag OpCapability Int64 into base blobs that
// never had it.

#ifdef SSPLAT_NATIVE_INT8

// Unsigned byte load (0..255).
[ForceInline] uint u8_load(uint* buf, uint idx) {
    return uint(((uint8_t*)buf)[idx]);
}
[ForceInline] uint u8_load(uint* buf, gindex_t idx) {
    return uint(((uint8_t*)buf)[idx]);
}
// Sign-extended int8 load (-128..127).
[ForceInline] int s8_load(uint* buf, uint idx) {
    return int(((int8_t*)buf)[idx]);
}
[ForceInline] int s8_load(uint* buf, gindex_t idx) {
    return int(((int8_t*)buf)[idx]);
}
[ForceInline] void u8_store(uint* buf, uint idx, uint v) {
    ((uint8_t*)buf)[idx] = uint8_t(v);
}
[ForceInline] void u8_store(uint* buf, gindex_t idx, uint v) {
    ((uint8_t*)buf)[idx] = uint8_t(v);
}

#else  // !SSPLAT_NATIVE_INT8

[ForceInline] uint u8_load(uint* buf, uint idx) {
    uint w = buf[idx >> 2];
    return (w >> ((idx & 3u) * 8u)) & 0xffu;
}
[ForceInline] uint u8_load(uint* buf, gindex_t idx) {
    uint w = buf[idx >> 2];
    return (w >> ((uint(idx) & 3u) * 8u)) & 0xffu;
}
[ForceInline] int s8_load(uint* buf, uint idx) {
    return (int(u8_load(buf, idx)) << 24) >> 24;
}
[ForceInline] int s8_load(uint* buf, gindex_t idx) {
    return (int(u8_load(buf, idx)) << 24) >> 24;
}
[ForceInline] void u8_store(uint* buf, uint idx, uint v) {
    uint shift = 8u * (idx & 3u);
    InterlockedAnd(buf[idx >> 2], ~(0xffu << shift));
    InterlockedOr(buf[idx >> 2], (v & 0xffu) << shift);
}
[ForceInline] void u8_store(uint* buf, gindex_t idx, uint v) {
    uint shift = 8u * (uint(idx) & 3u);
    InterlockedAnd(buf[idx >> 2], ~(0xffu << shift));
    InterlockedOr(buf[idx >> 2], (v & 0xffu) << shift);
}

#endif  // SSPLAT_NATIVE_INT8

#endif  // SSPLAT_INT8_COMPAT
