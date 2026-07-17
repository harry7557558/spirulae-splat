// Device sort / scan / stream-compaction kernels behind
// backend/common/SortScan.h (Vulkan backend). See the design notes in
// csrc/backend/vulkan/README.md.
//
// Radix sort: classic 3-phase LSD — per-partition histogram, spine
// (per-digit exclusive scan across partitions + cross-digit base scan),
// rank-and-scatter — 8 bits per pass with begin_bit/end_bit masking, stable,
// 32-bit or 64-bit keys (keys are compared UNSIGNED within the masked bit
// range; every call site sorts non-negative keys). No device-scope atomics:
// per-digit totals fall out of the spine scans.
//
// Subgroup-size handling: ballot-based ranking uses uint4 masks (lanes up to
// 128) and never assumes a width; the host asserts subgroup size >= 8 so the
// per-subgroup scratch bound (MAX_NSG) holds.

// ---------------------------------------------------------------------------
// shared configuration (host mirrors: backend/vulkan/SortScanVulkan.cpp)
// ---------------------------------------------------------------------------

#define SORT_WG 128         // threads per sort workgroup
#define SORT_ELEMS 16       // keys per thread
#define SORT_PART 2048      // SORT_WG * SORT_ELEMS
#define RADIX 256
#define MAX_NSG 16          // SORT_WG / min supported subgroup size (8)

#define SCAN_WG 256         // threads per scan workgroup
#define SCAN_ELEMS 8        // elements per thread
#define SCAN_BLOCK 2048     // SCAN_WG * SCAN_ELEMS

// ---------------------------------------------------------------------------
// subgroup ballot helpers (uint4 masks, subgroup width <= 128)
// ---------------------------------------------------------------------------

uint popcount4(uint4 v) {
    return countbits(v.x) + countbits(v.y) + countbits(v.z) + countbits(v.w);
}

// Bits [0, lane) set.
uint4 lane_lt_mask(uint lane) {
    uint word = lane >> 5;
    uint bit = lane & 31;
    uint4 m = uint4(0);
    [ForceUnroll]
    for (uint w = 0; w < 4; w++) {
        if (w < word) m[w] = 0xffffffffu;
        else if (w == word) m[w] = (bit == 0) ? 0u : (0xffffffffu >> (32 - bit));
    }
    return m;
}

uint lsb4(uint4 m) {
    if (m.x != 0) return firstbitlow(m.x);
    if (m.y != 0) return 32 + firstbitlow(m.y);
    if (m.z != 0) return 64 + firstbitlow(m.z);
    if (m.w != 0) return 96 + firstbitlow(m.w);
    return 0xffffffffu;
}

// Mask of active lanes whose `digit` equals this lane's, restricted to
// `valid` lanes. All lanes must call this (uniform control flow).
uint4 match_digit(uint digit, bool valid) {
    uint4 mask = WaveActiveBallot(valid);
    [ForceUnroll]
    for (uint b = 0; b < 8; b++) {
        bool bit = ((digit >> b) & 1u) != 0;
        uint4 bal = WaveActiveBallot(bit);
        mask &= bit ? bal : ~bal;
    }
    return mask;
}

// ---------------------------------------------------------------------------
// radix sort: histogram
// ---------------------------------------------------------------------------

struct HistParams {
    uint64_t* keys;      // u32 or u64 elements depending on entry
    uint32_t* part_hist; // [RADIX][num_parts], digit-major
    uint32_t n;
    uint32_t shift;
    uint32_t digit_mask;
    uint32_t num_parts;
    uint32_t parts_per_row;  // grid fold: part = wg.y * parts_per_row + wg.x
};

groupshared uint _hist_local[RADIX];

// KEY_LOAD(params, idx) yields the uint digit source.
#define DEF_HIST(NAME, KTYPE)                                              \
[shader("compute")]                                                        \
[numthreads(SORT_WG, 1, 1)]                                                \
void NAME(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,                \
          uniform HistParams p) {                                          \
    uint part = wg.y * p.parts_per_row + wg.x;                             \
    if (part >= p.num_parts) return;                                       \
    for (uint i = tid; i < RADIX; i += SORT_WG) _hist_local[i] = 0;        \
    GroupMemoryBarrierWithGroupSync();                                     \
    KTYPE* keys = (KTYPE*)p.keys;                                          \
    uint base = part * SORT_PART;                                          \
    for (uint e = 0; e < SORT_ELEMS; e++) {                                \
        uint idx = base + e * SORT_WG + tid;                               \
        if (idx < p.n) {                                                   \
            uint d = uint(keys[idx] >> p.shift) & p.digit_mask;            \
            InterlockedAdd(_hist_local[d], 1u);                            \
        }                                                                  \
    }                                                                      \
    GroupMemoryBarrierWithGroupSync();                                     \
    for (uint d = tid; d < RADIX; d += SORT_WG)                            \
        p.part_hist[d * p.num_parts + part] = _hist_local[d];              \
}

DEF_HIST(sort_hist_u32, uint32_t)
DEF_HIST(sort_hist_u64, uint64_t)

// ---------------------------------------------------------------------------
// radix sort: spine
// ---------------------------------------------------------------------------

struct SpineLocalParams {
    uint32_t* part_hist;   // in: counts; out: per-digit exclusive prefixes
    uint32_t* digit_total; // [RADIX] out
    uint32_t num_parts;
};

groupshared uint _spinel_tmp[SORT_WG];

// One workgroup per digit: exclusive scan of part_hist[d][*] in place,
// total to digit_total[d].
[shader("compute")]
[numthreads(SORT_WG, 1, 1)]
void sort_spine_local(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                      uniform SpineLocalParams p) {
    uint d = wg.x;
    uint running = 0;
    for (uint chunk = 0; chunk < p.num_parts; chunk += SORT_WG) {
        uint idx = chunk + tid;
        uint v = idx < p.num_parts ? p.part_hist[d * p.num_parts + idx] : 0;
        // Hillis-Steele inclusive scan over the workgroup.
        _spinel_tmp[tid] = v;
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll]
        for (uint off = 1; off < SORT_WG; off <<= 1) {
            uint x = tid >= off ? _spinel_tmp[tid - off] : 0;
            GroupMemoryBarrierWithGroupSync();
            _spinel_tmp[tid] += x;
            GroupMemoryBarrierWithGroupSync();
        }
        uint incl = _spinel_tmp[tid];
        if (idx < p.num_parts)
            p.part_hist[d * p.num_parts + idx] = running + incl - v;
        uint chunk_total = _spinel_tmp[SORT_WG - 1];
        GroupMemoryBarrierWithGroupSync();
        running += chunk_total;
    }
    if (tid == 0) p.digit_total[d] = running;
}

struct SpineGlobalParams {
    uint32_t* digit_total;  // [RADIX] in
    uint32_t* global_base;  // [RADIX] out (exclusive scan)
};

groupshared uint _spineg_tmp[RADIX];

// Single workgroup: exclusive scan of the 256 digit totals.
[shader("compute")]
[numthreads(RADIX, 1, 1)]
void sort_spine_global(uint tid: SV_GroupThreadID,
                       uniform SpineGlobalParams p) {
    uint v = p.digit_total[tid];
    _spineg_tmp[tid] = v;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll]
    for (uint off = 1; off < RADIX; off <<= 1) {
        uint x = tid >= off ? _spineg_tmp[tid - off] : 0;
        GroupMemoryBarrierWithGroupSync();
        _spineg_tmp[tid] += x;
        GroupMemoryBarrierWithGroupSync();
    }
    p.global_base[tid] = _spineg_tmp[tid] - v;
}

// ---------------------------------------------------------------------------
// radix sort: rank and scatter (stable)
// ---------------------------------------------------------------------------

struct DownsweepParams {
    uint64_t* keys_in;
    uint64_t* keys_out;
    uint32_t* vals_in;
    uint32_t* vals_out;
    uint32_t* part_offsets;  // [RADIX][num_parts] per-digit exclusive prefix
    uint32_t* global_base;   // [RADIX]
    uint32_t n;
    uint32_t shift;
    uint32_t digit_mask;
    uint32_t num_parts;
    uint32_t parts_per_row;
};

// Ranking scratch. Iterations process SORT_WG contiguous keys; stability
// needs rank = (prior iterations) + (prior subgroups this iteration) +
// (rank within subgroup), all in memory order.
groupshared uint _ds_sg_hist[RADIX * MAX_NSG];  // [digit][subgroup]
groupshared uint _ds_running[RADIX];            // prior-iterations count

#define DEF_DOWNSWEEP(NAME, KTYPE)                                         \
[shader("compute")]                                                        \
[numthreads(SORT_WG, 1, 1)]                                                \
void NAME(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,                \
          uniform DownsweepParams p) {                                     \
    uint part = wg.y * p.parts_per_row + wg.x;                             \
    if (part >= p.num_parts) return;                                       \
    KTYPE* keys_in = (KTYPE*)p.keys_in;                                    \
    KTYPE* keys_out = (KTYPE*)p.keys_out;                                  \
    uint lane = WaveGetLaneIndex();                                        \
    uint sg = tid / WaveGetLaneCount();                                    \
    uint nsg = (SORT_WG + WaveGetLaneCount() - 1) / WaveGetLaneCount();    \
    uint4 lt = lane_lt_mask(lane);                                         \
    for (uint i = tid; i < RADIX; i += SORT_WG) _ds_running[i] = 0;        \
    for (uint i = tid; i < RADIX * MAX_NSG; i += SORT_WG)                  \
        _ds_sg_hist[i] = 0;                                                \
    GroupMemoryBarrierWithGroupSync();                                     \
    uint base = part * SORT_PART;                                          \
    for (uint e = 0; e < SORT_ELEMS; e++) {                                \
        uint idx = base + e * SORT_WG + tid;                               \
        bool valid = idx < p.n;                                            \
        KTYPE key = valid ? keys_in[idx] : KTYPE(0);                       \
        uint d = uint(key >> p.shift) & p.digit_mask;                      \
        uint4 match = match_digit(d, valid);                               \
        uint rank_sg = popcount4(match & lt);                              \
        bool leader = valid && (lsb4(match) == lane);                      \
        if (leader) _ds_sg_hist[d * MAX_NSG + sg] = popcount4(match);      \
        GroupMemoryBarrierWithGroupSync();                                 \
        if (valid) {                                                       \
            uint prior_sg = 0;                                             \
            for (uint s = 0; s < sg; s++)                                  \
                prior_sg += _ds_sg_hist[d * MAX_NSG + s];                  \
            uint pos = p.global_base[d] +                                  \
                       p.part_offsets[d * p.num_parts + part] +            \
                       _ds_running[d] + prior_sg + rank_sg;                \
            keys_out[pos] = key;                                           \
            p.vals_out[pos] = p.vals_in[idx];                              \
        }                                                                  \
        GroupMemoryBarrierWithGroupSync();                                 \
        /* fold this iteration's counts into the running totals and */     \
        /* clear the per-subgroup scratch for the next iteration */        \
        for (uint dd = tid; dd < RADIX; dd += SORT_WG) {                   \
            uint total = 0;                                                \
            for (uint s = 0; s < nsg; s++) {                               \
                total += _ds_sg_hist[dd * MAX_NSG + s];                    \
                _ds_sg_hist[dd * MAX_NSG + s] = 0;                         \
            }                                                              \
            _ds_running[dd] += total;                                      \
        }                                                                  \
        GroupMemoryBarrierWithGroupSync();                                 \
    }                                                                      \
}

DEF_DOWNSWEEP(sort_downsweep_u32, uint32_t)
DEF_DOWNSWEEP(sort_downsweep_u64, uint64_t)

// ---------------------------------------------------------------------------
// two-level scans (inclusive_sum / exclusive_sum)
// ---------------------------------------------------------------------------

// Spec constant 0: 1 = inclusive, 0 = exclusive.
[SpecializationConstant]
const int kInclusive = 0;

struct ScanParams {
    uint64_t* input;       // typed per entry
    uint64_t* output;
    uint64_t* block_sums;  // [num_blocks]
    uint32_t n;
    uint32_t blocks_per_row;
};

groupshared int64_t _scan_tmp_i64[SCAN_WG];
groupshared int _scan_tmp_i32[SCAN_WG];

// uint8 flag arrays are read as packed 32-bit words (byte extraction):
// 8-bit loads / 8-bit pointee types would drag in shaderInt8 +
// 8-bit-storage capabilities the baseline deliberately avoids. Reading the
// word containing byte n-1 may touch up to 3 bytes past n; device
// allocations are 16-byte rounded (VulkanRuntime), so this stays in bounds.
int load_flag_byte(uint32_t* words, uint idx) {
    uint w = words[idx >> 2];
    return int((w >> ((idx & 3u) * 8u)) & 0xffu);
}

// LOAD(idx) is an expression over `p` yielding T; blocked layout (each
// thread owns SCAN_ELEMS consecutive elements) keeps intra-thread order
// trivial.
#define DEF_SCAN_BLOCKS(NAME, T, TMP, LOAD)                                \
[shader("compute")]                                                        \
[numthreads(SCAN_WG, 1, 1)]                                                \
void NAME(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,                \
          uniform ScanParams p) {                                          \
    uint block = wg.y * p.blocks_per_row + wg.x;                           \
    uint base = block * SCAN_BLOCK + tid * SCAN_ELEMS;                     \
    T* output = (T*)p.output;                                              \
    T* block_sums = (T*)p.block_sums;                                      \
    T loaded[SCAN_ELEMS];                                                  \
    T v[SCAN_ELEMS];                                                       \
    [ForceUnroll]                                                          \
    for (uint e = 0; e < SCAN_ELEMS; e++) {                                \
        uint idx = base + e;                                               \
        loaded[e] = idx < p.n ? T(LOAD(idx)) : T(0);                       \
        v[e] = loaded[e];                                                  \
    }                                                                      \
    T thread_total = T(0);                                                 \
    [ForceUnroll]                                                          \
    for (uint e = 0; e < SCAN_ELEMS; e++) {                                \
        T x = v[e];                                                        \
        v[e] = thread_total; /* exclusive within thread */                 \
        thread_total += x;                                                 \
    }                                                                      \
    TMP[tid] = thread_total;                                               \
    GroupMemoryBarrierWithGroupSync();                                     \
    [ForceUnroll]                                                          \
    for (uint off = 1; off < SCAN_WG; off <<= 1) {                         \
        T x = tid >= off ? TMP[tid - off] : T(0);                          \
        GroupMemoryBarrierWithGroupSync();                                 \
        TMP[tid] += x;                                                     \
        GroupMemoryBarrierWithGroupSync();                                 \
    }                                                                      \
    T thread_base = TMP[tid] - thread_total; /* exclusive across threads */\
    [ForceUnroll]                                                          \
    for (uint e = 0; e < SCAN_ELEMS; e++) {                                \
        uint idx = base + e;                                               \
        if (idx < p.n) {                                                   \
            T excl = thread_base + v[e];                                   \
            output[idx] = kInclusive != 0 ? excl + loaded[e] : excl;       \
        }                                                                  \
    }                                                                      \
    if (tid == SCAN_WG - 1)                                                \
        block_sums[block] = thread_base + thread_total;                    \
}

#define LOAD_I32(idx) (((int*)p.input)[idx])
#define LOAD_I64(idx) (((int64_t*)p.input)[idx])
#define LOAD_FLAG(idx) (load_flag_byte((uint32_t*)p.input, idx) != 0 ? 1 : 0)

DEF_SCAN_BLOCKS(scan_blocks_i32, int, _scan_tmp_i32, LOAD_I32)
DEF_SCAN_BLOCKS(scan_blocks_i64, int64_t, _scan_tmp_i64, LOAD_I64)
// bool/uint8 flags scanned into an int32 output (select_flagged).
DEF_SCAN_BLOCKS(scan_blocks_flags, int, _scan_tmp_i32, LOAD_FLAG)

struct ScanSpineParams {
    uint64_t* block_sums;  // scanned exclusively in place
    uint32_t num_blocks;
};

#define DEF_SCAN_SPINE(NAME, T, TMP)                                       \
[shader("compute")]                                                        \
[numthreads(SCAN_WG, 1, 1)]                                                \
void NAME(uint tid: SV_GroupThreadID, uniform ScanSpineParams p) {         \
    T* sums = (T*)p.block_sums;                                            \
    T running = T(0);                                                      \
    for (uint chunk = 0; chunk < p.num_blocks; chunk += SCAN_WG) {         \
        uint idx = chunk + tid;                                            \
        T v = idx < p.num_blocks ? sums[idx] : T(0);                       \
        TMP[tid] = v;                                                      \
        GroupMemoryBarrierWithGroupSync();                                 \
        [ForceUnroll]                                                      \
        for (uint off = 1; off < SCAN_WG; off <<= 1) {                     \
            T x = tid >= off ? TMP[tid - off] : T(0);                      \
            GroupMemoryBarrierWithGroupSync();                             \
            TMP[tid] += x;                                                 \
            GroupMemoryBarrierWithGroupSync();                             \
        }                                                                  \
        if (idx < p.num_blocks) sums[idx] = running + TMP[tid] - v;        \
        T chunk_total = TMP[SCAN_WG - 1];                                  \
        GroupMemoryBarrierWithGroupSync();                                 \
        running += chunk_total;                                            \
    }                                                                      \
}

DEF_SCAN_SPINE(scan_spine_i32, int, _scan_tmp_i32)
DEF_SCAN_SPINE(scan_spine_i64, int64_t, _scan_tmp_i64)

struct ScanAddParams {
    uint64_t* output;
    uint64_t* block_sums;  // exclusive block bases
    uint32_t n;
    uint32_t blocks_per_row;
};

#define DEF_SCAN_ADD(NAME, T)                                              \
[shader("compute")]                                                        \
[numthreads(SCAN_WG, 1, 1)]                                                \
void NAME(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,                \
          uniform ScanAddParams p) {                                       \
    uint block = wg.y * p.blocks_per_row + wg.x;                           \
    T* output = (T*)p.output;                                              \
    T addend = ((T*)p.block_sums)[block];                                  \
    uint base = block * SCAN_BLOCK + tid * SCAN_ELEMS;                     \
    [ForceUnroll]                                                          \
    for (uint e = 0; e < SCAN_ELEMS; e++) {                                \
        uint idx = base + e;                                               \
        if (idx < p.n) output[idx] += addend;                              \
    }                                                                      \
}

DEF_SCAN_ADD(scan_add_i32, int)
DEF_SCAN_ADD(scan_add_i64, int64_t)

// ---------------------------------------------------------------------------
// stream compaction (select_flagged)
// ---------------------------------------------------------------------------

struct ScatterFlaggedParams {
    uint32_t* input;
    uint32_t* flags;  // uint8 array read as packed words (see load_flag_byte)
    int* scan;        // exclusive scan of flags
    uint32_t* output;
    uint32_t n;
    uint32_t blocks_per_row;  // SCAN_BLOCK-element blocks, matching the scan
};

[shader("compute")]
[numthreads(SCAN_WG, 1, 1)]
void scatter_flagged_u32(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                         uniform ScatterFlaggedParams p) {
    uint block = wg.y * p.blocks_per_row + wg.x;
    uint base = block * SCAN_BLOCK + tid * SCAN_ELEMS;
    [ForceUnroll]
    for (uint e = 0; e < SCAN_ELEMS; e++) {
        uint idx = base + e;
        if (idx < p.n && load_flag_byte(p.flags, idx) != 0)
            p.output[p.scan[idx]] = p.input[idx];
    }
}
