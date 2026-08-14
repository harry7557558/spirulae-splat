// PrepProgress.cpp -- see PrepProgress.h.

#include "app/gui/PrepProgress.h"

namespace gui {

void RunProgress::reset() {
    std::lock_guard<std::mutex> lk(_mu);
    for (StageProgress& s : _st) s = StageProgress{};
    _cur = Stage::Frames;
    _pending.clear();
}

void RunProgress::enter(Stage s, const std::string& detail) {
    std::lock_guard<std::mutex> lk(_mu);
    for (int i = 0; i < (int)s; i++)
        if (_st[i].status == StageStatus::Pending)
            _st[i].status = StageStatus::Skipped;
        else if (_st[i].status == StageStatus::Running)
            _st[i].status = StageStatus::Done;
    _st[(int)s].status = StageStatus::Running;
    _st[(int)s].detail = detail;
    _cur = s;
}

void RunProgress::count(Stage s, int64_t done, int64_t total) {
    std::lock_guard<std::mutex> lk(_mu);
    StageProgress& p = _st[(int)s];
    p.done = done;
    p.total = total;
    p.fraction = total > 0 ? (float)((double)done / (double)total) : -1.0f;
}

void RunProgress::fraction(Stage s, float f) {
    std::lock_guard<std::mutex> lk(_mu);
    _st[(int)s].fraction = f;
}

void RunProgress::detail(Stage s, const std::string& text) {
    std::lock_guard<std::mutex> lk(_mu);
    _st[(int)s].detail = text;
}

void RunProgress::mark(Stage s, StageStatus st) {
    std::lock_guard<std::mutex> lk(_mu);
    _st[(int)s].status = st;
    if (st == StageStatus::Done) _st[(int)s].fraction = 1.0f;
}

void RunProgress::finish(StageStatus st) {
    std::lock_guard<std::mutex> lk(_mu);
    for (StageProgress& p : _st)
        if (p.status == StageStatus::Running) p.status = st;
}

void RunProgress::note(const std::string& text, bool detail) {
    std::lock_guard<std::mutex> lk(_mu);
    _pending.push_back({_cur, detail, text});
}

void RunProgress::note(Stage s, const std::string& text, bool detail) {
    std::lock_guard<std::mutex> lk(_mu);
    _pending.push_back({s, detail, text});
}

std::vector<RunLine> RunProgress::drain() {
    std::lock_guard<std::mutex> lk(_mu);
    std::vector<RunLine> out;
    out.swap(_pending);
    return out;
}

StageProgress RunProgress::stage(Stage s) const {
    std::lock_guard<std::mutex> lk(_mu);
    return _st[(int)s];
}

Stage RunProgress::current() const {
    std::lock_guard<std::mutex> lk(_mu);
    return _cur;
}

bool RunProgress::ran(Stage s) const {
    std::lock_guard<std::mutex> lk(_mu);
    const StageStatus st = _st[(int)s].status;
    return st != StageStatus::Pending;
}

}  // namespace gui
