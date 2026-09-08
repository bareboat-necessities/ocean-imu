#pragma once

// Host-only observer. Templates read public state; callbacks in the overlay
// receive const views of the exact local measurement operands. No estimator
// pointer, writable view, shadow recursion or random-number call is retained.
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace ou3_brmm_trace {
struct Trace {
    std::ofstream samples, events;
    unsigned long index = 0, sequence = 0;
    bool pending = false, outer_pre = false, inner_pre = false, bias_pre = false;
    int predict = 0, acc_calls = 0, mag_calls = 0, mag_delivered = 0;
    int acc_applied = 0, mag_applied = 0, s_calls = 0, s_applied = 0;
    int resets = 0, tilt_resets = 0, reference_writes = 0, guard_calls = 0;
    bool guard_identical = true;
    float dt = 0, engagement = 0, excess = 0;
    Trace() {
        const char* prefix = std::getenv("OU3_BRMM_TRACE");
        if (!prefix || !*prefix) throw std::runtime_error("OU3_BRMM_TRACE is required");
        samples.open(std::string(prefix)+".samples.csv");
        events.open(std::string(prefix)+".events.csv");
        if (!samples || !events) throw std::runtime_error("cannot open BRMM trace");
        samples << std::setprecision(9);
        events << std::setprecision(9);
        samples << "index,dt,outer_pre,outer_post,inner_pre,inner_post,bias_pre,bias_post,"
                   "predict,acc_calls,acc_applied,mag_delivered,mag_calls,mag_applied,s_calls,s_applied,"
                   "resets,tilt_resets,reference_writes,guard_calls,guard_identical,engagement,excess,"
                   "tau,sigma,rs,pseudo_period,bias_norm,mag_lock,mag_refined\n";
        events << "index,sequence,kind,bx,by,bz,wx,wy,wz,R00,R01,R02,R10,R11,R12,R20,R21,R22\n";
    }
};
inline Trace& get() { static Trace t; return t; }

template<class Fusion>
void begin(const Fusion& f, float dt) {
    auto& t = get();
    if (t.pending) throw std::runtime_error("previous BRMM sample was not finished");
    t.pending = true;
    t.dt = dt;
    t.outer_pre = f.isLive();
    t.inner_pre = f.raw().isAdaptiveLive();
    t.bias_pre = f.raw().mekf().acc_bias_updates_enabled();
    t.predict = t.acc_calls = t.acc_applied = t.mag_delivered = t.mag_calls = t.mag_applied = 0;
    t.s_calls = t.s_applied = t.resets = t.tilt_resets = t.reference_writes = t.guard_calls = 0;
    t.guard_identical = true;
    t.engagement = t.excess = 0;
}

template<class V>
void guard(const V& raw, const V& conditioned, float engagement, float excess) {
    auto& t = get();
    ++t.guard_calls;
    t.guard_identical = (raw.array() == conditioned.array()).all();
    t.engagement = engagement;
    t.excess = excess;
}

template<class V, class Rotation, class Covariance>
void measurement(const char* kind, const V& body, const Rotation& rotation, const Covariance& r) {
    auto& t = get();
    if (!t.pending) throw std::runtime_error("measurement outside BRMM sample");
    if (std::string(kind) == "acc") ++t.acc_applied;
    else if (std::string(kind) == "mag") ++t.mag_applied;
    else ++t.s_applied;
    const auto world = (rotation.transpose()*body).eval();
    t.events << t.index << ',' << t.sequence++ << ',' << kind;
    for (int i = 0; i < 3; ++i) t.events << ',' << body(i);
    for (int i = 0; i < 3; ++i) t.events << ',' << world(i);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) t.events << ',' << r(i,j);
    t.events << '\n';
}

template<class Fusion>
void finish(const Fusion& f) {
    auto& t = get();
    if (!t.pending) return; // A repeated snapshot is not another IMU sample.
    const auto& a = f.raw();
    const auto& m = a.mekf();
    t.samples << t.index << ',' << t.dt << ',' << t.outer_pre << ',' << f.isLive()
              << ',' << t.inner_pre << ',' << a.isAdaptiveLive()
              << ',' << t.bias_pre << ',' << m.acc_bias_updates_enabled()
              << ',' << t.predict << ',' << t.acc_calls << ',' << t.acc_applied
              << ',' << t.mag_delivered << ',' << t.mag_calls << ',' << t.mag_applied
              << ',' << t.s_calls << ',' << t.s_applied << ',' << t.resets
              << ',' << t.tilt_resets << ',' << t.reference_writes << ',' << t.guard_calls
              << ',' << t.guard_identical << ',' << t.engagement << ',' << t.excess
              << ',' << a.getTauApplied() << ',' << a.getSigmaApplied() << ',' << a.getRSApplied()
              << ',' << m.get_pseudo_update_period_s() << ',' << m.get_acc_bias().norm()
              << ',' << f.hasMagNorthLock() << ',' << f.hasRefinedMagReference() << '\n';
    if (!t.samples || !t.events) throw std::runtime_error("BRMM trace write failed");
    t.pending = false;
    ++t.index;
}
} // namespace ou3_brmm_trace
