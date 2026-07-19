// Batched masked quantile radix-select (Vulkan backend), mirroring
// csrc/Quantile.cu. Selects the q-quantile of |x| over the strictly-positive
// finite elements of each row, excluding exact-max occurrences (pass 0
// subtracts them), via 4 radix passes over 8-bit digits of the float's bit
// image. All accumulation is integer (histograms + atomics), so results are
// bit-exact across backends and devices -> tight parity channel.
//
// Spec constants (IDs by declaration order):
//   0: kInvert -- 0: output the quantile value; 1: output its reciprocal
//      (the CUDA template's invert_quantile parameter).

[SpecializationConstant]
const int kInvert = 0;

static const int kTileElems = 4096;
static const uint kPosInfBits = 0x7f800000u;

struct QuantileParams {
    float* x;                // [B*N]
    uint* row_max_bits;      // [B]
    uint* row_max_count;     // [B]
    uint* row_prefix;        // [B]
    uint* row_rank;          // [B]
    uint* row_done;          // [B]
    uint* row_hist;          // [B*256]
    float* out;              // [B]
    int B, N;
    int tiles_per_row;
    int pass;                // 0..3
    float q;
    uint wgs_per_row;        // tile axis folds across (wg.x, wg.y)
};

[ForceInline]
bool _qt_pos_finite(uint u) {
    return (u > 0u) && (u < kPosInfBits);
}

groupshared uint _qt_red[256];

[ForceInline]
uint _qt_block_sum(uint v, uint tid) {
    _qt_red[tid] = v;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) _qt_red[tid] += _qt_red[tid + s];
        GroupMemoryBarrierWithGroupSync();
    }
    return _qt_red[0];
}

[ForceInline]
uint _qt_block_max(uint v, uint tid) {
    _qt_red[tid] = v;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) _qt_red[tid] = max(_qt_red[tid], _qt_red[tid + s]);
        GroupMemoryBarrierWithGroupSync();
    }
    return _qt_red[0];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void qt_find_row_max(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                     uniform QuantileParams p) {
    int tile_id = int(wg.y * p.wgs_per_row + wg.x);
    int row = tile_id / p.tiles_per_row;
    // No early return before the block reduce -- all threads must reach the
    // barriers. Out-of-range tiles contribute 0.
    bool valid = row < p.B;

    uint local_max = 0u;
    if (valid) {
        int tile_in_row = tile_id - row * p.tiles_per_row;
        int start = tile_in_row * kTileElems;
        int end = min(start + kTileElems, p.N);
        int64_t base = (int64_t)row * (int64_t)p.N;
        for (int i = start + int(tid); i < end; i += 256) {
            uint u = asuint(abs(p.x[base + i]));
            if (_qt_pos_finite(u) && u > local_max) local_max = u;
        }
    }
    uint block_max = _qt_block_max(local_max, tid);
    if (valid && tid == 0 && block_max != 0u)
        InterlockedMax(p.row_max_bits[row], block_max);
}

groupshared uint _qt_shist[256];

[shader("compute")]
[numthreads(256, 1, 1)]
void qt_hist_pass(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                  uniform QuantileParams p) {
    _qt_shist[tid] = 0u;
    GroupMemoryBarrierWithGroupSync();

    int tile_id = int(wg.y * p.wgs_per_row + wg.x);
    int row = tile_id / p.tiles_per_row;
    bool valid = row < p.B;
    if (valid && p.row_done[row] != 0u) valid = false;

    uint local_max_count = 0u;
    if (valid) {
        int tile_in_row = tile_id - row * p.tiles_per_row;
        int start = tile_in_row * kTileElems;
        int end = min(start + kTileElems, p.N);
        int64_t base = (int64_t)row * (int64_t)p.N;
        uint max_bits = p.row_max_bits[row];

        if (p.pass == 0) {
            for (int i = start + int(tid); i < end; i += 256) {
                uint u = asuint(abs(p.x[base + i]));
                if (!_qt_pos_finite(u)) continue;
                InterlockedAdd(_qt_shist[(u >> 24) & 255u], 1u);
                if (u == max_bits) local_max_count++;
            }
        } else {
            uint prefix = p.row_prefix[row];
            uint mask = 0xFFFFFFFFu << (32 - 8 * p.pass);
            uint shift = uint(24 - 8 * p.pass);
            for (int i = start + int(tid); i < end; i += 256) {
                uint u = asuint(abs(p.x[base + i]));
                if (!_qt_pos_finite(u)) continue;
                if ((u & mask) == prefix)
                    InterlockedAdd(_qt_shist[(u >> shift) & 255u], 1u);
            }
        }
    }
    // Reduce outside the valid-branch so every thread reaches the barriers.
    uint block_max_count = _qt_block_sum(local_max_count, tid);
    if (valid && p.pass == 0 && tid == 0 && block_max_count != 0u)
        InterlockedAdd(p.row_max_count[row], block_max_count);
    GroupMemoryBarrierWithGroupSync();
    if (valid) {
        uint c = _qt_shist[tid];
        if (c != 0u)
            InterlockedAdd(p.row_hist[uint(row) * 256u + tid], c);
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void qt_select_pass(uint3 wg: SV_GroupID, uint tid: SV_GroupThreadID,
                    uniform QuantileParams p) {
    int row = int(wg.x * 256 + tid);
    if (row >= p.B) return;
    if (p.row_done[row] != 0u) return;

    uint max_bits = p.row_max_bits[row];
    uint max_count = p.row_max_count[row];

    if (max_bits == 0u) {
        // No strictly positive finite values exist.
        p.out[row] = 1.0f;
        p.row_done[row] = 1u;
        return;
    }

    uint hbase = uint(row) * 256u;

    if (p.pass == 0) {
        uint64_t total = 0;
        for (int b = 0; b < 256; b++) total += p.row_hist[hbase + b];
        uint masked_count =
            (total > uint64_t(max_count)) ? uint(total - max_count) : 0u;
        if (masked_count == 0u) {
            p.out[row] = kInvert != 0 ? 2.0f / asfloat(max_bits)
                                      : asfloat(max_bits) * 0.5f;
            p.row_done[row] = 1u;
            return;
        }
        // Lower-rank quantile convention.
        uint rank = uint(min(1.0f, max(0.0f, p.q)) * float(masked_count - 1u));
        p.row_rank[row] = rank;
        p.row_prefix[row] = 0u;
    }

    uint prefix = p.row_prefix[row];
    uint rank = p.row_rank[row];

    uint shift = uint(24 - 8 * p.pass);
    uint mask_prefix =
        (p.pass == 0) ? 0u : (0xFFFFFFFFu << (32 - 8 * p.pass));
    bool exclude_max =
        (p.pass == 0) ? true : ((prefix & mask_prefix) == (max_bits & mask_prefix));
    uint max_bucket = (max_bits >> shift) & 255u;

    uint cum = 0u;
    uint chosen = 255u;
    for (uint b = 0; b < 256u; b++) {
        uint c = p.row_hist[hbase + b];
        if (exclude_max && b == max_bucket) c -= max_count;
        if (rank < cum + c) {
            chosen = b;
            break;
        }
        cum += c;
    }

    if (p.pass == 3) {
        prefix |= chosen;
        p.out[row] =
            kInvert != 0 ? 1.0f / asfloat(prefix) : asfloat(prefix);
        p.row_done[row] = 1u;
    } else {
        prefix |= (chosen << shift);
        p.row_prefix[row] = prefix;
        p.row_rank[row] = rank - cum;
    }
}
