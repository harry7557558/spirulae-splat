// DataManager — see DataManager.h for the public contract.

#include "DataManager.h"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <future>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>


// ===========================================================================
// Small utilities
// ===========================================================================
namespace {

// File suffix (lowercased extension). Returns "" if none.
std::string lower_suffix(const std::string& path) {
    auto dot = path.find_last_of('.');
    if (dot == std::string::npos) return {};
    std::string s = path.substr(dot);
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

// Bounded MPMC queue used by the disk-mode prefetch pipeline. One mutex
// guards both ends — this is for low-rate batch coordination (jobs/sec, not
// allocations/sec), so the contention cost is negligible.
template<typename T>
class BoundedQueue {
public:
    explicit BoundedQueue(size_t cap) : _cap(cap) {}

    // Push, blocking while full. Returns false if closed before space is
    // available.
    bool push(T v) {
        std::unique_lock<std::mutex> lk(_m);
        _cv_not_full.wait(lk, [&]{ return _closed || _q.size() < _cap; });
        if (_closed) return false;
        _q.push_back(std::move(v));
        _cv_not_empty.notify_one();
        return true;
    }

    // Pop, blocking while empty. Returns false if closed and drained.
    bool pop(T& out) {
        std::unique_lock<std::mutex> lk(_m);
        _cv_not_empty.wait(lk, [&]{ return _closed || !_q.empty(); });
        if (_q.empty()) return false;  // closed-and-drained
        out = std::move(_q.front());
        _q.pop_front();
        _cv_not_full.notify_one();
        return true;
    }

    // Wake everyone; subsequent pop()s drain remaining items then return false.
    void close() {
        {
            std::lock_guard<std::mutex> lk(_m);
            _closed = true;
        }
        _cv_not_empty.notify_all();
        _cv_not_full.notify_all();
    }

private:
    size_t                  _cap;
    std::deque<T>           _q;
    std::mutex              _m;
    std::condition_variable _cv_not_empty;
    std::condition_variable _cv_not_full;
    bool                    _closed = false;
};

// One element of a per-image grouping. Each group is homogeneous in (W,H),
// so a batch sampled from a single group has the contiguous shape the engine
// expects.
struct IndexGroup {
    int32_t              width  = 0;
    int32_t              height = 0;
    std::vector<int32_t> indices;       // dataset-global indices, shuffled
    size_t               cursor = 0;    // next index to emit

    // Refill from `indices` end-of-list. Reshuffles for non-eval datasets.
    void rewind(std::mt19937_64& rng, bool eval) {
        if (!eval) std::shuffle(indices.begin(), indices.end(), rng);
        cursor = 0;
    }
};

// Build (W,H)-keyed groups from a flat index list.
std::vector<IndexGroup> build_index_groups(
    const std::vector<int32_t>& flat_indices,
    const std::vector<int32_t>& widths,
    const std::vector<int32_t>& heights)
{
    std::unordered_map<uint64_t, IndexGroup> by_shape;
    for (int32_t i : flat_indices) {
        uint64_t key = ((uint64_t)(uint32_t)widths[i] << 32)
                     | (uint64_t)(uint32_t)heights[i];
        auto& g = by_shape[key];
        if (g.indices.empty()) {
            g.width  = widths[i];
            g.height = heights[i];
        }
        g.indices.push_back(i);
    }
    std::vector<IndexGroup> out;
    out.reserve(by_shape.size());
    for (auto& kv : by_shape) out.emplace_back(std::move(kv.second));
    return out;
}

// Group-size-weighted sampler — picks a group with probability proportional
// to its size, mirroring datamanager.py:IndexGroups::random_idx().
class GroupSampler {
public:
    GroupSampler() = default;
    explicit GroupSampler(const std::vector<IndexGroup>& groups) {
        _weights.reserve(groups.size());
        for (auto& g : groups) _weights.push_back((double)g.indices.size());
        _dist = std::discrete_distribution<size_t>(_weights.begin(), _weights.end());
    }
    size_t operator()(std::mt19937_64& rng) { return _dist(rng); }
    bool   empty() const { return _weights.empty(); }
private:
    std::vector<double>                _weights;
    std::discrete_distribution<size_t> _dist;
};

// ===========================================================================
// Image decoders
// ===========================================================================
//
// All decoders write directly into a pre-allocated batch slot `dst` (pointer
// to the start of row 0 of this image's slot in the [B,H,W,C] buffer).
//
// `expected_h`, `expected_w` are the IndexGroup's promised dimensions; if the
// decoded image disagrees we throw (the IndexGroup invariant was violated).

void decode_rgb_into(const std::string& path,
                     int expected_h, int expected_w,
                     PixelDType dtype,
                     uint8_t* dst)
{
    int w, h, ch;
    if (dtype == PixelDType::UINT16) {
        stbi_us* img = stbi_load_16(path.c_str(), &w, &h, &ch, 3);
        if (!img) throw std::runtime_error("DataManager: failed to load 16-bit RGB '" + path + "': " + stbi_failure_reason());
        if (w != expected_w || h != expected_h) {
            stbi_image_free(img);
            throw std::runtime_error("DataManager: rgb shape mismatch for '" + path + "'");
        }
        std::memcpy(dst, img, (size_t)w * h * 3 * sizeof(stbi_us));
        stbi_image_free(img);
    } else if (dtype == PixelDType::UINT8) {
        stbi_uc* img = stbi_load(path.c_str(), &w, &h, &ch, 3);
        if (!img) throw std::runtime_error("DataManager: failed to load 8-bit RGB '" + path + "': " + stbi_failure_reason());
        if (w != expected_w || h != expected_h) {
            stbi_image_free(img);
            throw std::runtime_error("DataManager: rgb shape mismatch for '" + path + "'");
        }
        std::memcpy(dst, img, (size_t)w * h * 3);
        stbi_image_free(img);
    } else {
        throw std::runtime_error("DataManager: float RGB inputs not supported in stb_image path");
    }
}

void decode_mask_into(const std::string& path,
                      int expected_h, int expected_w,
                      uint8_t* dst)
{
    int w, h, ch;
    stbi_uc* img = stbi_load(path.c_str(), &w, &h, &ch, 1);
    if (!img) throw std::runtime_error("DataManager: failed to load mask '" + path + "': " + stbi_failure_reason());
    if (w != expected_w || h != expected_h) {
        stbi_image_free(img);
        throw std::runtime_error("DataManager: mask shape mismatch for '" + path + "'");
    }
    for (size_t i = 0; i < (size_t)w * h; ++i)
        dst[i] = (uint8_t)(img[i] != 0);
    stbi_image_free(img);
}

void decode_depth_into(const std::string& path,
                       int expected_h, int expected_w,
                       PixelDType dtype,
                       uint8_t* dst)
{
    int w, h, ch;
    if (dtype == PixelDType::UINT16) {
        stbi_us* img = stbi_load_16(path.c_str(), &w, &h, &ch, 1);
        if (!img) throw std::runtime_error("DataManager: failed to load 16-bit depth '" + path + "': " + stbi_failure_reason());
        if (w != expected_w || h != expected_h) {
            stbi_image_free(img);
            throw std::runtime_error("DataManager: depth shape mismatch for '" + path + "'");
        }
        std::memcpy(dst, img, (size_t)w * h * sizeof(stbi_us));
        stbi_image_free(img);
    } else {
        throw std::runtime_error("DataManager: only 16-bit depth PNGs are supported in stb_image path");
    }
}

void decode_normal_into(const std::string& path,
                        int expected_h, int expected_w,
                        uint8_t* dst)
{
    int w, h, ch;
    stbi_uc* img = stbi_load(path.c_str(), &w, &h, &ch, 3);
    if (!img) throw std::runtime_error("DataManager: failed to load normal '" + path + "': " + stbi_failure_reason());
    if (w != expected_w || h != expected_h) {
        stbi_image_free(img);
        throw std::runtime_error("DataManager: normal shape mismatch for '" + path + "'");
    }
    std::memcpy(dst, img, (size_t)w * h * 3);
    stbi_image_free(img);
}

// Detect on-disk RGB dtype (8-bit vs 16-bit). Both are decided by a single
// stbi_is_16_bit() probe.
PixelDType probe_pixel_dtype(const std::string& path,
                             PixelDType fallback = PixelDType::UINT8)
{
    if (path.empty()) return fallback;
    if (stbi_is_16_bit(path.c_str())) return PixelDType::UINT16;
    return PixelDType::UINT8;
}

} // anonymous namespace


// ===========================================================================
// DecodedBatch::build_views
// ===========================================================================
void DecodedBatch::build_views() {
    auto mk = [](const void* data, uint32_t es,
                 std::vector<int64_t> shape) -> TorchTensorView {
        return TorchTensorView((uint64_t)data, es, std::move(shape));
    };
    int64_t B = num, H = height, W = width;
    viewmats_view    = mk(viewmats.data(),    4, {B, 4LL, 4LL});
    intrins_view     = mk(intrins.data(),     4, {B, 4LL});
    dist_coeffs_view = mk(dist_coeffs.data(), 4, {B, 10LL});

    if (B > 0 && H > 0 && W > 0 && !rgb_buffer.empty()) {
        rgb_view = mk(rgb_buffer.data(), pixel_dtype_size(rgb_dtype),
                      {B, H, W, 3LL});
    }
    if (!mask_buffer.empty()) {
        mask_view = mk(mask_buffer.data(), 1, {B, H, W, 1LL});
    }
    if (!depth_buffer.empty()) {
        depth_view = mk(depth_buffer.data(), pixel_dtype_size(depth_dtype),
                        {B, H, W, 1LL});
    }
    if (!normal_buffer.empty()) {
        normal_view = mk(normal_buffer.data(), pixel_dtype_size(normal_dtype),
                         {B, H, W, 3LL});
    }
}


// ===========================================================================
// DataManagerImpl — PIMPL
// ===========================================================================
class DataManagerImpl {
public:

    DataManagerImpl(
        DataManagerConfig config,
        CameraModelType   model,
        std::vector<std::string> image_filenames,
        std::vector<std::string> mask_filenames,
        std::vector<std::string> depth_filenames,
        std::vector<std::string> normal_filenames,
        std::vector<int32_t>     widths,
        std::vector<int32_t>     heights,
        std::vector<float>       viewmats,
        std::vector<float>       intrins,
        std::vector<float>       dist_coeffs,
        std::vector<int32_t>     train_indices,
        std::vector<int32_t>     val_indices);

    ~DataManagerImpl();

    const DecodedBatch& next_train_batch();
    const DecodedBatch* next_val_batch();

    int64_t num_train()    const { return (int64_t)_train_indices.size(); }
    int64_t num_val()      const { return (int64_t)_val_indices.size(); }
    bool    has_val()      const { return !_val_indices.empty(); }
    CacheMode cache_mode() const { return _cfg.cache_mode; }

    bool has_masks()   const { return !_mask_filenames.empty()   && _cfg.load_masks; }
    bool has_depths()  const { return !_depth_filenames.empty()  && _cfg.load_depths; }
    bool has_normals() const { return !_normal_filenames.empty() && _cfg.load_normals; }

private:

    // ---- Shared input data ------------------------------------------------
    DataManagerConfig         _cfg;
    CameraModelType           _model;
    std::vector<std::string>  _image_filenames;
    std::vector<std::string>  _mask_filenames;
    std::vector<std::string>  _depth_filenames;
    std::vector<std::string>  _normal_filenames;
    std::vector<int32_t>      _widths, _heights;
    std::vector<float>        _viewmats;     // [N,4,4] flat
    std::vector<float>        _intrins;      // [N,4]
    std::vector<float>        _dist_coeffs;  // [N,10]
    std::vector<int32_t>      _train_indices;
    std::vector<int32_t>      _val_indices;

    // Per-image detected RGB / depth / normal dtypes. Cached at construction
    // so the disk path doesn't re-probe headers on every fetch.
    std::vector<PixelDType>   _rgb_dtype;
    std::vector<PixelDType>   _depth_dtype;
    std::vector<PixelDType>   _normal_dtype;

    // ---- CPU mode preloaded buffers --------------------------------------
    //
    // Per-image, the decoded payload. Empty for images that don't have the
    // modality on disk. For DISK mode these stay empty.
    std::vector<std::vector<uint8_t>> _rgb_cache;
    std::vector<std::vector<uint8_t>> _mask_cache;
    std::vector<std::vector<uint8_t>> _depth_cache;
    std::vector<std::vector<uint8_t>> _normal_cache;

    // ---- Batch sampling state --------------------------------------------
    std::mutex                _sampling_mu;
    std::mt19937_64           _rng;
    std::vector<IndexGroup>   _train_groups;
    std::vector<IndexGroup>   _val_groups;
    GroupSampler              _train_sampler;
    GroupSampler              _val_sampler;

    // The currently-returned-to-caller batches. We hold one slot per kind so
    // the reference returned by next_*_batch() stays valid until the next
    // call to that same getter.
    DecodedBatch              _cpu_train_batch;
    DecodedBatch              _cpu_val_batch;

    // ---- DISK mode: prefetch pipeline ------------------------------------
    struct DecodeJob {
        int64_t              batch_id   = 0;
        int32_t              slot       = 0;   // row inside batch buffer
        int32_t              ds_index   = 0;   // dataset-global image index
        std::shared_ptr<DecodedBatch>     batch;
        std::shared_ptr<std::atomic<int>> remaining;  // decrement on done
        // The ready queue (train or val) this batch belongs to. The worker
        // that decrements `remaining` to zero takes ownership of publishing
        // the batch to this queue; no separate completer / scheduler-side
        // wait is involved, so all four worker pools stay saturated while
        // multiple batches are in flight.
        BoundedQueue<std::shared_ptr<DecodedBatch>>* ready_q = nullptr;
    };

    bool _disk_started = false;

    std::unique_ptr<BoundedQueue<DecodeJob>> _q_rgb;
    std::unique_ptr<BoundedQueue<DecodeJob>> _q_mask;
    std::unique_ptr<BoundedQueue<DecodeJob>> _q_depth;
    std::unique_ptr<BoundedQueue<DecodeJob>> _q_normal;
    std::unique_ptr<BoundedQueue<std::shared_ptr<DecodedBatch>>> _q_ready_train;
    std::unique_ptr<BoundedQueue<std::shared_ptr<DecodedBatch>>> _q_ready_val;

    std::vector<std::thread> _workers;
    std::thread              _scheduler;
    std::atomic<bool>        _stop{false};
    std::atomic<int64_t>     _next_batch_id{0};

    // The shared_ptr<DecodedBatch> currently held alive on behalf of the
    // most recent next_*_batch() return value.
    std::shared_ptr<DecodedBatch> _last_train_held;
    std::shared_ptr<DecodedBatch> _last_val_held;

    // ---- Setup helpers ---------------------------------------------------
    void probe_dtypes();
    void preload_cpu_cache();
    void allocate_batch(DecodedBatch& b,
                        const IndexGroup& g,
                        const std::vector<int32_t>& ds_indices);
    void fill_camera_params(DecodedBatch& b);

    // CPU-mode synchronous batch fetch helpers.
    DecodedBatch& next_batch_cpu(DecodedBatch& slot,
                                 std::vector<IndexGroup>& groups,
                                 GroupSampler& sampler,
                                 int batch_size);
    void fill_batch_from_cache(DecodedBatch& b);

    // Disk-mode helpers.
    void start_disk_pipeline();
    void stop_disk_pipeline();
    void scheduler_loop();
    void worker_loop_rgb();
    void worker_loop_mask();
    void worker_loop_depth();
    void worker_loop_normal();
    void enqueue_batch(IndexGroup& group, int batch_size,
                       BoundedQueue<std::shared_ptr<DecodedBatch>>& ready_q);
};


// ===========================================================================
// Construction
// ===========================================================================
DataManagerImpl::DataManagerImpl(
    DataManagerConfig          config,
    CameraModelType            model,
    std::vector<std::string>   image_filenames,
    std::vector<std::string>   mask_filenames,
    std::vector<std::string>   depth_filenames,
    std::vector<std::string>   normal_filenames,
    std::vector<int32_t>       widths,
    std::vector<int32_t>       heights,
    std::vector<float>         viewmats,
    std::vector<float>         intrins,
    std::vector<float>         dist_coeffs,
    std::vector<int32_t>       train_indices,
    std::vector<int32_t>       val_indices)
    : _cfg(config), _model(model),
      _image_filenames(std::move(image_filenames)),
      _mask_filenames(std::move(mask_filenames)),
      _depth_filenames(std::move(depth_filenames)),
      _normal_filenames(std::move(normal_filenames)),
      _widths(std::move(widths)), _heights(std::move(heights)),
      _viewmats(std::move(viewmats)),
      _intrins(std::move(intrins)),
      _dist_coeffs(std::move(dist_coeffs)),
      _train_indices(std::move(train_indices)),
      _val_indices(std::move(val_indices))
{
    int64_t N = (int64_t)_image_filenames.size();
    if ((int64_t)_widths.size()  != N ||
        (int64_t)_heights.size() != N ||
        (int64_t)_viewmats.size()    != N * 16 ||
        (int64_t)_intrins.size()     != N * 4 ||
        (int64_t)_dist_coeffs.size() != N * 10) {
        throw std::runtime_error("DataManager: per-camera array length mismatch");
    }
    if (!_mask_filenames.empty()   && (int64_t)_mask_filenames.size()   != N)
        throw std::runtime_error("DataManager: mask_filenames length mismatch");
    if (!_depth_filenames.empty()  && (int64_t)_depth_filenames.size()  != N)
        throw std::runtime_error("DataManager: depth_filenames length mismatch");
    if (!_normal_filenames.empty() && (int64_t)_normal_filenames.size() != N)
        throw std::runtime_error("DataManager: normal_filenames length mismatch");

    uint64_t seed = _cfg.seed != 0 ? _cfg.seed
        : (uint64_t)std::chrono::steady_clock::now().time_since_epoch().count();
    _rng.seed(seed);

    probe_dtypes();

    _train_groups = build_index_groups(_train_indices, _widths, _heights);
    _val_groups   = build_index_groups(_val_indices,   _widths, _heights);
    for (auto& g : _train_groups) g.rewind(_rng, /*eval=*/false);
    for (auto& g : _val_groups)   g.rewind(_rng, /*eval=*/true);
    _train_sampler = GroupSampler(_train_groups);
    _val_sampler   = GroupSampler(_val_groups);

    if (_cfg.cache_mode == CacheMode::CPU) {
        preload_cpu_cache();
    } else {
        start_disk_pipeline();
    }
}

DataManagerImpl::~DataManagerImpl() {
    if (_disk_started) stop_disk_pipeline();
}


// ---- Probe per-file dtypes (cheap; only reads PNG headers) -----------------
void DataManagerImpl::probe_dtypes() {
    int64_t N = (int64_t)_image_filenames.size();
    _rgb_dtype.assign((size_t)N, PixelDType::UINT8);
    for (int64_t i = 0; i < N; ++i) {
        _rgb_dtype[i] = probe_pixel_dtype(_image_filenames[i], PixelDType::UINT8);
    }
    if (has_depths()) {
        _depth_dtype.assign((size_t)N, PixelDType::UINT16);
        // depth assumed 16-bit grayscale PNG everywhere — we don't auto-detect
        // float because stb_image doesn't carry float PNG anyway.
    }
    if (has_normals()) {
        _normal_dtype.assign((size_t)N, PixelDType::UINT8);
    }
}


// ===========================================================================
// CPU mode — preload everything in parallel
// ===========================================================================
void DataManagerImpl::preload_cpu_cache() {
    int64_t N = (int64_t)_image_filenames.size();
    _rgb_cache.resize(N);
    if (has_masks())   _mask_cache.resize(N);
    if (has_depths())  _depth_cache.resize(N);
    if (has_normals()) _normal_cache.resize(N);

    // Pick a sensible worker count for the preload sweep — capped by the
    // number of images, so we don't spawn 32 threads to decode 4 images.
    int n_threads = (int)std::max(1u, std::thread::hardware_concurrency());
    n_threads = std::min(n_threads, (int)N);

    std::atomic<int64_t> next_idx{0};
    std::atomic<int64_t> done_count{0};
    std::mutex           print_mu;
    std::vector<std::thread> ts;
    ts.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        ts.emplace_back([&]() {
            for (;;) {
                int64_t i = next_idx.fetch_add(1);
                if (i >= N) return;
                int W = _widths[i], H = _heights[i];

                // RGB.
                {
                    PixelDType dt = _rgb_dtype[i];
                    size_t bytes = (size_t)W * H * 3 * pixel_dtype_size(dt);
                    _rgb_cache[i].assign(bytes, 0);
                    decode_rgb_into(_image_filenames[i], H, W, dt, _rgb_cache[i].data());
                }
                if (has_masks() && !_mask_filenames[i].empty()) {
                    _mask_cache[i].assign((size_t)W * H, 0);
                    decode_mask_into(_mask_filenames[i], H, W, _mask_cache[i].data());
                }
                if (has_depths() && !_depth_filenames[i].empty()) {
                    PixelDType dt = _depth_dtype[i];
                    size_t bytes = (size_t)W * H * pixel_dtype_size(dt);
                    _depth_cache[i].assign(bytes, 0);
                    decode_depth_into(_depth_filenames[i], H, W, dt, _depth_cache[i].data());
                }
                if (has_normals() && !_normal_filenames[i].empty()) {
                    _normal_cache[i].assign((size_t)W * H * 3, 0);
                    decode_normal_into(_normal_filenames[i], H, W, _normal_cache[i].data());
                }

                int64_t d = done_count.fetch_add(1) + 1;
                {
                    std::lock_guard<std::mutex> lk(print_mu);
                    std::fprintf(stderr, "\rLoading images %lld/%lld",
                                 (long long)d, (long long)N);
                    std::fflush(stderr);
                }
            }
        });
    }
    for (auto& th : ts) th.join();
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
}


// ===========================================================================
// Batch allocation + camera params + cache gather
// ===========================================================================
void DataManagerImpl::allocate_batch(
    DecodedBatch& b,
    const IndexGroup& g,
    const std::vector<int32_t>& ds_indices)
{
    int B = (int)ds_indices.size();
    int H = g.height, W = g.width;

    b.width  = W;
    b.height = H;
    b.num    = B;
    b.model  = _model;
    b.indices = ds_indices;

    b.viewmats.assign((size_t)B * 16, 0.0f);
    b.intrins.assign((size_t)B * 4, 0.0f);
    b.dist_coeffs.assign((size_t)B * 10, 0.0f);

    // RGB: pick the homogeneous dtype across the batch (the first image
    // wins; mixed dtypes within a (W,H) group are rare and treated as 8-bit).
    PixelDType rgb_dt = _rgb_dtype.empty() ? PixelDType::UINT8 : _rgb_dtype[ds_indices[0]];
    b.rgb_dtype = rgb_dt;
    b.rgb_buffer.assign((size_t)B * H * W * 3 * pixel_dtype_size(rgb_dt), 0);

    if (has_masks())
        b.mask_buffer.assign((size_t)B * H * W, 0);
    else
        b.mask_buffer.clear();

    if (has_depths()) {
        PixelDType d_dt = PixelDType::UINT16;
        b.depth_dtype = d_dt;
        b.depth_buffer.assign((size_t)B * H * W * pixel_dtype_size(d_dt), 0);
    } else {
        b.depth_buffer.clear();
    }

    if (has_normals()) {
        b.normal_dtype = PixelDType::UINT8;
        b.normal_buffer.assign((size_t)B * H * W * 3, 0);
    } else {
        b.normal_buffer.clear();
    }

    fill_camera_params(b);
}

void DataManagerImpl::fill_camera_params(DecodedBatch& b) {
    int B = b.num;
    for (int j = 0; j < B; ++j) {
        int32_t i = b.indices[j];
        std::memcpy(&b.viewmats[(size_t)j * 16],    &_viewmats[(size_t)i * 16],    16 * sizeof(float));
        std::memcpy(&b.intrins[(size_t)j * 4],      &_intrins[(size_t)i * 4],       4 * sizeof(float));
        std::memcpy(&b.dist_coeffs[(size_t)j * 10], &_dist_coeffs[(size_t)i * 10], 10 * sizeof(float));
    }
}

void DataManagerImpl::fill_batch_from_cache(DecodedBatch& b) {
    int B = b.num, H = b.height, W = b.width;
    size_t rgb_row = (size_t)H * W * 3 * pixel_dtype_size(b.rgb_dtype);
    size_t mask_row = (size_t)H * W;
    size_t depth_row = (size_t)H * W * pixel_dtype_size(b.depth_dtype);
    size_t normal_row = (size_t)H * W * 3;

    for (int j = 0; j < B; ++j) {
        int32_t i = b.indices[j];
        std::memcpy(b.rgb_buffer.data() + (size_t)j * rgb_row,
                    _rgb_cache[i].data(), rgb_row);
        if (!b.mask_buffer.empty() && !_mask_cache.empty() && !_mask_cache[i].empty())
            std::memcpy(b.mask_buffer.data() + (size_t)j * mask_row,
                        _mask_cache[i].data(), mask_row);
        if (!b.depth_buffer.empty() && !_depth_cache.empty() && !_depth_cache[i].empty())
            std::memcpy(b.depth_buffer.data() + (size_t)j * depth_row,
                        _depth_cache[i].data(), depth_row);
        if (!b.normal_buffer.empty() && !_normal_cache.empty() && !_normal_cache[i].empty())
            std::memcpy(b.normal_buffer.data() + (size_t)j * normal_row,
                        _normal_cache[i].data(), normal_row);
    }
}


// ===========================================================================
// CPU-mode batch fetch
// ===========================================================================
DecodedBatch& DataManagerImpl::next_batch_cpu(
    DecodedBatch& slot,
    std::vector<IndexGroup>& groups,
    GroupSampler& sampler,
    int batch_size)
{
    std::lock_guard<std::mutex> lk(_sampling_mu);
    if (sampler.empty()) {
        throw std::runtime_error("DataManager: no indices configured for this split");
    }
    size_t gi = sampler(_rng);
    auto& g = groups[gi];

    int eff_bs = std::min(batch_size, (int)g.indices.size());
    eff_bs = std::max(eff_bs, 1);
    std::vector<int32_t> picks;
    picks.reserve(eff_bs);
    for (int k = 0; k < eff_bs; ++k) {
        if (g.cursor >= g.indices.size()) g.rewind(_rng, /*eval=*/false);
        picks.push_back(g.indices[g.cursor++]);
    }

    allocate_batch(slot, g, picks);
    fill_batch_from_cache(slot);
    slot.build_views();
    return slot;
}


// ===========================================================================
// DISK mode pipeline
// ===========================================================================
void DataManagerImpl::start_disk_pipeline() {
    int prefetch = std::max(1, _cfg.prefetch_batches);
    // Job queues are sized generously — each batch enqueues B jobs per
    // modality so cap = prefetch * batch_size keeps the scheduler from
    // running ahead of the consumer.
    int job_cap = std::max(1, prefetch * std::max(_cfg.train_batch_size, _cfg.val_batch_size));
    _q_rgb    = std::make_unique<BoundedQueue<DecodeJob>>((size_t)job_cap);
    _q_mask   = std::make_unique<BoundedQueue<DecodeJob>>((size_t)job_cap);
    _q_depth  = std::make_unique<BoundedQueue<DecodeJob>>((size_t)job_cap);
    _q_normal = std::make_unique<BoundedQueue<DecodeJob>>((size_t)job_cap);
    _q_ready_train = std::make_unique<BoundedQueue<std::shared_ptr<DecodedBatch>>>((size_t)prefetch);
    _q_ready_val   = std::make_unique<BoundedQueue<std::shared_ptr<DecodedBatch>>>((size_t)prefetch);

    // Dedicated RGB workers + (optionally) dedicated mask workers — the
    // user asked for "rgb+mask" to share a thread group conceptually, but
    // a clean modality split (one pool per queue) is more robust than
    // tagging jobs through a shared queue. The mask pool shares the
    // workers_rgb knob: spawn ceil(N/2) mask workers up to a floor of 1.
    int n_rgb_workers = std::max(1, _cfg.workers_rgb);
    for (int t = 0; t < n_rgb_workers; ++t)
        _workers.emplace_back(&DataManagerImpl::worker_loop_rgb, this);

    if (has_masks()) {
        int n_mask = std::max(1, n_rgb_workers / 2);
        for (int t = 0; t < n_mask; ++t)
            _workers.emplace_back(&DataManagerImpl::worker_loop_mask, this);
    }

    if (has_depths()) {
        int n = std::max(1, _cfg.workers_depth);
        for (int t = 0; t < n; ++t)
            _workers.emplace_back(&DataManagerImpl::worker_loop_depth, this);
    }
    if (has_normals()) {
        int n = std::max(1, _cfg.workers_normal);
        for (int t = 0; t < n; ++t)
            _workers.emplace_back(&DataManagerImpl::worker_loop_normal, this);
    }

    _scheduler = std::thread(&DataManagerImpl::scheduler_loop, this);
    _disk_started = true;
}

void DataManagerImpl::stop_disk_pipeline() {
    _stop.store(true);
    _q_rgb->close();
    _q_mask->close();
    _q_depth->close();
    _q_normal->close();
    _q_ready_train->close();
    _q_ready_val->close();
    if (_scheduler.joinable()) _scheduler.join();
    for (auto& th : _workers) if (th.joinable()) th.join();
    _workers.clear();
}


void DataManagerImpl::worker_loop_rgb() {
    DecodeJob job;
    while (_q_rgb->pop(job)) {
        if (_stop.load()) return;
        DecodedBatch& b = *job.batch;
        int H = b.height, W = b.width;
        size_t row = (size_t)H * W * 3 * pixel_dtype_size(b.rgb_dtype);
        uint8_t* dst = b.rgb_buffer.data() + (size_t)job.slot * row;
        decode_rgb_into(_image_filenames[job.ds_index], H, W,
                        b.rgb_dtype, dst);
        if (job.remaining->fetch_sub(1) == 1) {
            job.batch->build_views();
            job.ready_q->push(job.batch);
        }
    }
}

void DataManagerImpl::worker_loop_mask() {
    DecodeJob job;
    while (_q_mask->pop(job)) {
        if (_stop.load()) return;
        DecodedBatch& b = *job.batch;
        int H = b.height, W = b.width;
        size_t row = (size_t)H * W;
        uint8_t* dst = b.mask_buffer.data() + (size_t)job.slot * row;
        decode_mask_into(_mask_filenames[job.ds_index], H, W, dst);
        if (job.remaining->fetch_sub(1) == 1) {
            job.batch->build_views();
            job.ready_q->push(job.batch);
        }
    }
}

void DataManagerImpl::worker_loop_depth() {
    DecodeJob job;
    while (_q_depth->pop(job)) {
        if (_stop.load()) return;
        DecodedBatch& b = *job.batch;
        int H = b.height, W = b.width;
        size_t row = (size_t)H * W * pixel_dtype_size(b.depth_dtype);
        uint8_t* dst = b.depth_buffer.data() + (size_t)job.slot * row;
        decode_depth_into(_depth_filenames[job.ds_index], H, W,
                          b.depth_dtype, dst);
        if (job.remaining->fetch_sub(1) == 1) {
            job.batch->build_views();
            job.ready_q->push(job.batch);
        }
    }
}

void DataManagerImpl::worker_loop_normal() {
    DecodeJob job;
    while (_q_normal->pop(job)) {
        if (_stop.load()) return;
        DecodedBatch& b = *job.batch;
        int H = b.height, W = b.width;
        size_t row = (size_t)H * W * 3;
        uint8_t* dst = b.normal_buffer.data() + (size_t)job.slot * row;
        decode_normal_into(_normal_filenames[job.ds_index], H, W, dst);
        if (job.remaining->fetch_sub(1) == 1) {
            job.batch->build_views();
            job.ready_q->push(job.batch);
        }
    }
}


// Scheduler: produces decoded-batch shells, hands jobs to per-modality
// queues, and routes the finished batch to the right ready queue.
void DataManagerImpl::scheduler_loop() {
    while (!_stop.load()) {
        // Alternate fill order: prefer training to validation, refill val
        // when ready_val is near-empty. (Simple heuristic.)
        if (!_train_groups.empty()) {
            std::unique_lock<std::mutex> lk(_sampling_mu);
            size_t gi = _train_sampler(_rng);
            auto& g = _train_groups[gi];
            lk.unlock();
            enqueue_batch(g, _cfg.train_batch_size, *_q_ready_train);
        }
        if (!_val_groups.empty()) {
            std::unique_lock<std::mutex> lk(_sampling_mu);
            size_t gi = _val_sampler(_rng);
            auto& g = _val_groups[gi];
            lk.unlock();
            enqueue_batch(g, _cfg.val_batch_size, *_q_ready_val);
        }
    }
}

void DataManagerImpl::enqueue_batch(
    IndexGroup& group, int batch_size,
    BoundedQueue<std::shared_ptr<DecodedBatch>>& ready_q)
{
    int eff_bs = std::min(batch_size, (int)group.indices.size());
    eff_bs = std::max(eff_bs, 1);

    std::vector<int32_t> picks;
    picks.reserve(eff_bs);
    {
        std::lock_guard<std::mutex> lk(_sampling_mu);
        for (int k = 0; k < eff_bs; ++k) {
            if (group.cursor >= group.indices.size()) group.rewind(_rng, /*eval=*/false);
            picks.push_back(group.indices[group.cursor++]);
        }
    }

    auto batch = std::make_shared<DecodedBatch>();
    allocate_batch(*batch, group, picks);

    // Number of decode jobs per image — RGB always, plus each enabled modality.
    int per_image = 1;
    if (has_masks())   per_image += 1;
    if (has_depths())  per_image += 1;
    if (has_normals()) per_image += 1;
    auto remaining = std::make_shared<std::atomic<int>>(eff_bs * per_image);

    // Helper: a single per-job completion decrement. When the counter flips
    // to 0 (returns 1 before the subtract), publish the batch ourselves —
    // no scheduler wait, no promise/future. Backpressure comes from the
    // bounded ready queue.
    auto try_publish = [&ready_q, batch, remaining]() {
        if (remaining->fetch_sub(1) == 1) {
            batch->build_views();
            ready_q.push(batch);
        }
    };

    // Stage all jobs into local lists first, so each modality queue is
    // pushed in one tight burst (more chance for workers to pick them up
    // in parallel) and we don't interleave producer-side contention with
    // worker wakeups.
    std::vector<DecodeJob> rgb_jobs, mask_jobs, depth_jobs, normal_jobs;
    rgb_jobs.reserve(eff_bs);
    if (has_masks())   mask_jobs.reserve(eff_bs);
    if (has_depths())  depth_jobs.reserve(eff_bs);
    if (has_normals()) normal_jobs.reserve(eff_bs);

    for (int j = 0; j < eff_bs; ++j) {
        int32_t i = picks[j];

        DecodeJob jb;
        jb.slot      = j;
        jb.ds_index  = i;
        jb.batch     = batch;
        jb.remaining = remaining;
        jb.ready_q   = &ready_q;
        jb.batch_id  = _next_batch_id.fetch_add(1);

        rgb_jobs.push_back(jb);
        if (has_masks()) {
            if (!_mask_filenames[i].empty()) mask_jobs.push_back(jb);
            else                              try_publish();  // skip slot
        }
        if (has_depths()) {
            if (!_depth_filenames[i].empty()) depth_jobs.push_back(jb);
            else                               try_publish();
        }
        if (has_normals()) {
            if (!_normal_filenames[i].empty()) normal_jobs.push_back(jb);
            else                                try_publish();
        }
    }

    // Hand off to the workers. The bounded ready queue (cap =
    // prefetch_batches) provides backpressure: once `prefetch_batches`
    // batches are decoded but unconsumed, the publishing worker blocks
    // here, while other workers keep draining their modality queues.
    // Modality job queues are also bounded (cap = prefetch_batches *
    // train_batch_size), so this scheduler call returns as soon as those
    // queues have room — and the next batch's jobs start enqueueing
    // immediately. This is what keeps all N workers busy when batch
    // size is 1.
    for (auto& jb : rgb_jobs)    if (!_q_rgb->push(std::move(jb))) return;
    for (auto& jb : mask_jobs)   if (!_q_mask->push(std::move(jb))) return;
    for (auto& jb : depth_jobs)  if (!_q_depth->push(std::move(jb))) return;
    for (auto& jb : normal_jobs) if (!_q_normal->push(std::move(jb))) return;
}


// ===========================================================================
// CPU vs DISK dispatch — public-facing fetch
// ===========================================================================
const DecodedBatch& DataManagerImpl::next_train_batch() {
    if (_cfg.cache_mode == CacheMode::CPU) {
        return next_batch_cpu(_cpu_train_batch, _train_groups,
                              _train_sampler, _cfg.train_batch_size);
    }
    std::shared_ptr<DecodedBatch> b;
    if (!_q_ready_train->pop(b))
        throw std::runtime_error("DataManager: train prefetch queue closed");
    _last_train_held = b;
    return *_last_train_held;
}

const DecodedBatch* DataManagerImpl::next_val_batch() {
    if (!has_val()) return nullptr;
    if (_cfg.cache_mode == CacheMode::CPU) {
        return &next_batch_cpu(_cpu_val_batch, _val_groups,
                               _val_sampler, _cfg.val_batch_size);
    }
    std::shared_ptr<DecodedBatch> b;
    if (!_q_ready_val->pop(b))
        throw std::runtime_error("DataManager: val prefetch queue closed");
    _last_val_held = b;
    return _last_val_held.get();
}


// ===========================================================================
// Public DataManager thin facade
// ===========================================================================
DataManager::DataManager(
    DataManagerConfig          config,
    CameraModelType            model,
    std::vector<std::string>   image_filenames,
    std::vector<std::string>   mask_filenames,
    std::vector<std::string>   depth_filenames,
    std::vector<std::string>   normal_filenames,
    std::vector<int32_t>       widths,
    std::vector<int32_t>       heights,
    std::vector<float>         viewmats,
    std::vector<float>         intrins,
    std::vector<float>         dist_coeffs,
    std::vector<int32_t>       train_indices,
    std::vector<int32_t>       val_indices)
{
    _impl = std::make_unique<DataManagerImpl>(
        std::move(config), model,
        std::move(image_filenames), std::move(mask_filenames),
        std::move(depth_filenames), std::move(normal_filenames),
        std::move(widths), std::move(heights),
        std::move(viewmats), std::move(intrins), std::move(dist_coeffs),
        std::move(train_indices), std::move(val_indices));
}

DataManager::~DataManager() = default;

const DecodedBatch& DataManager::next_train_batch()        { return _impl->next_train_batch(); }
const DecodedBatch* DataManager::next_val_batch()          { return _impl->next_val_batch(); }
int64_t   DataManager::num_train()      const              { return _impl->num_train(); }
int64_t   DataManager::num_val()        const              { return _impl->num_val(); }
bool      DataManager::has_val()        const              { return _impl->has_val(); }
CacheMode DataManager::cache_mode()     const              { return _impl->cache_mode(); }
bool      DataManager::has_masks()      const              { return _impl->has_masks(); }
bool      DataManager::has_depths()     const              { return _impl->has_depths(); }
bool      DataManager::has_normals()    const              { return _impl->has_normals(); }
