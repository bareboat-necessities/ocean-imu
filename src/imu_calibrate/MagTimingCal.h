#pragma once
// Copyright 2026, Mikhail Grushinskiy
// Relative delay identification against the delivered gyro clock. No absolute
// conversion epoch or independent-noise confidence interval is claimed.
#include "util/MagneticRotation.h"

namespace imu_cal {
enum class MagTimingReason : uint8_t {
    Ready, MorePairs, DelayRange, SegmentMismatch, WeakSignal, BroadProfile,
    FieldResidual, NoImprovement, TimingNoise
};
inline const char* magTimingReason(MagTimingReason reason) {
    switch (reason) {
    case MagTimingReason::Ready: return "Delay validated";
    case MagTimingReason::MorePairs: return "Need usable samples";
    case MagTimingReason::DelayRange: return "Delay outside range";
    case MagTimingReason::SegmentMismatch: return "Captures disagree";
    case MagTimingReason::WeakSignal: return "Reverse faster";
    case MagTimingReason::BroadProfile: return "Weak timing signal";
    case MagTimingReason::FieldResidual: return "Field/gyro mismatch";
    case MagTimingReason::NoImprovement: return "No timing improvement";
    case MagTimingReason::TimingNoise: return "Timing too noisy";
    }
    return "Timing invalid";
}
struct MagTimingFit {
    uint32_t delay_us = 0;
    float fit_width_s = 0, residual_time_s = 0;
    float train_rms_rad = 0, validation_rms_rad = 0, improvement = 0;
    uint32_t train_pairs = 0, validation_pairs = 0;
    MagTimingReason reason = MagTimingReason::MorePairs;
    bool ok = false;
};

class MagTimingCal {
    using V3 = ocean_imu::magnetic::V3;
    // 180 ms pair + 80 ms candidate lag + one 20 ms accepted gyro interval.
    using History = ocean_imu::magnetic::RotationHistory<128, 280000>;
    static constexpr unsigned N = 41; // 0..80 ms, 2 ms resolution
    struct Scores {
        std::array<double, N> cost{}; unsigned pairs = 0;
    };
public:
    void reset() { history_.reset(); scores_ = {}; have_ = false; seen_ = false; gaps_ = missing_ = 0; }
    void resetSegment(unsigned segment) {
        if (segment < 2) scores_[segment] = {};
        history_.reset(); have_ = seen_ = false;
    }
    unsigned pairs(unsigned segment) const { return segment < 2 ? scores_[segment].pairs : 0; }
    unsigned gaps() const { return gaps_; }
    unsigned missingHistory() const { return missing_; }
    // Capture progress only. Saving ALWAYS requires fit() with disjoint data.
    MagTimingFit segmentQuality(unsigned segment) const {
        return segment < 2 ? evaluate(scores_[segment], scores_[segment]) : MagTimingFit{};
    }
    // Two disjoint capture segments: segment 0 fits, segment 1 validates.
    // Use unaveraged spatially calibrated magnetic vectors and calibrated gyro.
    void add(uint32_t now, const V3& gyro, bool have, uint32_t sequence,
             uint32_t observed_us, const V3& field, unsigned segment) {
        using ocean_imu::magnetic::Status;
        const auto status = history_.push(now, gyro);
        if (status != Status::Ready && status != Status::First) { ++gaps_; have_ = false; return; }
        if (segment > 1 || !have || !field.allFinite() || field.norm() < 5 || field.norm() > 200) {
            have_ = false; return;
        }
        if (seen_ && sequence == sequence_) return;
        if (seen_ && (int32_t(sequence-sequence_) <= 0 || int32_t(observed_us-last_observed_) <= 0)) {
            have_ = false; return;
        }
        seen_ = true; sequence_ = sequence; last_observed_ = observed_us;
        if (int32_t(now-observed_us) < 0) { have_ = false; return; }
        if (!have_ || segment != segment_ || uint32_t(observed_us-anchor_us_) > 180000) {
            anchor(field, observed_us, segment); return;
        }
        if (uint32_t(observed_us-anchor_us_) < 120000) return;
        std::array<double, N> costs{};
        const float norm = field.norm();
        if (std::abs(norm-anchor_norm_) > .1f*anchor_norm_) { anchor(field, observed_us, segment); return; }
        bool complete = true;
        for (unsigned i=0; i<N; ++i) {
            History::Rotation r;
            if (history_.between(anchor_us_-i*2000, observed_us-i*2000, r) != Status::Ready) {
                ++missing_; complete = false; break;
            }
            const double e2 = (field/norm - r.q*anchor_unit_).squaredNorm();
            // Bounded loss limits isolated field disturbances; every candidate
            // sees exactly the same pairs and the same loss.
            costs[i] = std::min(e2, .04);
        }
        if (complete) {
            auto& s=scores_[segment];
            ++s.pairs;
            for (unsigned i=0; i<N; ++i) s.cost[i] += costs[i];
        }
        anchor(field, observed_us, segment);
    }
    MagTimingFit fit() const { return evaluate(scores_[0], scores_[1]); }
private:
    static MagTimingFit evaluate(const Scores& train, const Scores& validation) {
        MagTimingFit out;
        out.train_pairs=train.pairs; out.validation_pairs=validation.pairs;
        if (out.train_pairs < 50 || out.validation_pairs < 50) return out;
        // A scalar delay needs changing observable magnetic rotation, not a
        // rank-two gyro covariance. Profile separation below measures that
        // directly and rejects stillness, constant rate and field-axis spins.
        std::array<double,N> a{}, b{};
        for(unsigned i=0;i<N;++i) { a[i]=train.cost[i]/out.train_pairs; b[i]=validation.cost[i]/out.validation_pairs; }
        if (*std::max_element(a.begin(),a.end())-*std::min_element(a.begin(),a.end())<1e-6 ||
            *std::max_element(b.begin(),b.end())-*std::min_element(b.begin(),b.end())<1e-6) {
            out.reason=MagTimingReason::WeakSignal; return out;
        }
        const unsigned k=unsigned(std::min_element(a.begin(),a.end())-a.begin());
        const unsigned v=unsigned(std::min_element(b.begin(),b.end())-b.begin());
        // Boundary solutions, constant-rate ambiguity, and disagreement between
        // segments cannot authorize a saved correction.
        out.delay_us=k*2000;
        out.train_rms_rad=float(std::sqrt(a[k]));
        out.validation_rms_rad=float(std::sqrt(b[k]));
        out.improvement=float(1-b[k]/std::max(b[0],1e-12));
        if (k<5 || k>N-6 || v<5 || v>N-6) { out.reason=MagTimingReason::DelayRange; return out; }
        if (std::abs(int(k)-int(v))>4) { out.reason=MagTimingReason::SegmentMismatch; return out; }
        const double curvature=std::min(a[k-5]-a[k],a[k+5]-a[k]);
        const double validation_curve=std::min(b[v-5]-b[v],b[v+5]-b[v]);
        if (curvature<1e-6 || validation_curve<1e-6) { out.reason=MagTimingReason::WeakSignal; return out; }
        const double tolerance=std::max(1e-6, .05*a[k]);
        unsigned lo=k, hi=k;
        // Include disconnected competing minima, not just the local basin.
        for(unsigned i=0;i<N;++i) if(a[i]<=a[k]+tolerance) { lo=std::min(lo,i); hi=std::max(hi,i); }
        if(hi-lo>12) { out.reason=MagTimingReason::BroadProfile; return out; }
        out.fit_width_s=std::max(.002f, float(double(hi-lo)*.001));
        // Residual-to-time conversion is a model allowance, NOT a measured
        // jitter standard deviation. Keep a 25 ms floor for asynchronous phase.
        const double sensitivity=std::sqrt(validation_curve)/.01;
        out.residual_time_s=float(std::max(.025, std::sqrt(b[k])/sensitivity));
        if (out.train_rms_rad>=.05f || out.validation_rms_rad>=.05f) out.reason=MagTimingReason::FieldResidual;
        else if (out.improvement<.05f) out.reason=MagTimingReason::NoImprovement;
        else if (out.residual_time_s>.060f) out.reason=MagTimingReason::TimingNoise;
        else { out.reason=MagTimingReason::Ready; out.ok=true; }
        return out;
    }
    void anchor(const V3& m,uint32_t t,unsigned segment) {
        have_=true; anchor_norm_=m.norm(); anchor_unit_=m/anchor_norm_; anchor_us_=t; segment_=segment;
    }
    History history_{};
    std::array<Scores,2> scores_{};
    V3 anchor_unit_=V3::Zero();
    float anchor_norm_=0;
    uint32_t anchor_us_=0, sequence_=0, last_observed_=0;
    unsigned segment_=0, gaps_=0, missing_=0;
    bool have_=false, seen_=false;
};
} // namespace imu_cal
