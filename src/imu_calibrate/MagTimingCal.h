#pragma once
// Copyright 2026, Mikhail Grushinskiy
// Relative delay identification against the delivered gyro clock. No absolute
// conversion epoch or independent-noise confidence interval is claimed.
#include "util/MagneticRotation.h"

namespace imu_cal {
struct MagTimingFit {
    uint32_t delay_us = 0;
    float fit_width_s = 0, residual_time_s = 0;
    float train_rms_rad = 0, validation_rms_rad = 0, improvement = 0;
    uint32_t train_pairs = 0, validation_pairs = 0;
    bool ok = false;
};

class MagTimingCal {
    using V3 = ocean_imu::magnetic::V3;
    using History = ocean_imu::magnetic::RotationHistory<128, 250000>;
    static constexpr unsigned N = 41; // 0..80 ms, 2 ms resolution
    struct Scores {
        std::array<double, N> cost{}; unsigned pairs = 0, gyro_samples = 0;
        V3 gyro_sum=V3::Zero(); ocean_imu::magnetic::M3 gyro_outer=ocean_imu::magnetic::M3::Zero();
    };
public:
    void reset() { history_.reset(); scores_ = {}; have_ = false; seen_ = false; }
    // Two disjoint capture segments: segment 0 fits, segment 1 validates.
    // Use unaveraged spatially calibrated magnetic vectors and calibrated gyro.
    void add(uint32_t now, const V3& gyro, bool have, uint32_t sequence,
             uint32_t observed_us, const V3& field, unsigned segment) {
        using ocean_imu::magnetic::Status;
        const auto status = history_.push(now, gyro);
        if (status != Status::Ready && status != Status::First) { have_ = false; return; }
        if (segment > 1 || !have || !field.allFinite() || field.norm() < 5 || field.norm() > 200) {
            have_ = false; return;
        }
        auto& excitation=scores_[segment];
        ++excitation.gyro_samples; excitation.gyro_sum+=gyro;
        excitation.gyro_outer+=gyro*gyro.transpose();
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
                complete = false; break;
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
    MagTimingFit fit() const {
        MagTimingFit out;
        out.train_pairs=scores_[0].pairs; out.validation_pairs=scores_[1].pairs;
        if (out.train_pairs < 50 || out.validation_pairs < 50) return out;
        for (const auto& s:scores_) {
            if (s.gyro_samples<200) return out;
            const V3 mean=s.gyro_sum/float(s.gyro_samples);
            const ocean_imu::magnetic::M3 covariance=s.gyro_outer/float(s.gyro_samples)-mean*mean.transpose();
            Eigen::SelfAdjointEigenSolver<ocean_imu::magnetic::M3> es(covariance);
            if (es.info()!=Eigen::Success || es.eigenvalues()[1]<.02f) return out;
        }
        std::array<double,N> a{}, b{};
        for(unsigned i=0;i<N;++i) { a[i]=scores_[0].cost[i]/out.train_pairs; b[i]=scores_[1].cost[i]/out.validation_pairs; }
        const unsigned k=unsigned(std::min_element(a.begin(),a.end())-a.begin());
        const unsigned v=unsigned(std::min_element(b.begin(),b.end())-b.begin());
        // Boundary solutions, constant-rate ambiguity, and disagreement between
        // segments cannot authorize a saved correction.
        if (k<5 || k>N-6 || v<5 || v>N-6 || std::abs(int(k)-int(v))>4) return out;
        const double curvature=std::min(a[k-5]-a[k],a[k+5]-a[k]);
        const double validation_curve=std::min(b[v-5]-b[v],b[v+5]-b[v]);
        if (curvature<1e-6 || validation_curve<1e-6) return out;
        const double tolerance=std::max(1e-6, .05*a[k]);
        unsigned lo=k, hi=k;
        // Include disconnected competing minima, not just the local basin.
        for(unsigned i=0;i<N;++i) if(a[i]<=a[k]+tolerance) { lo=std::min(lo,i); hi=std::max(hi,i); }
        if(hi-lo>12) return out;
        out.delay_us=k*2000;
        out.fit_width_s=std::max(.002f, float(hi-lo)*.001f);
        out.train_rms_rad=float(std::sqrt(a[k]));
        out.validation_rms_rad=float(std::sqrt(b[k]));
        out.improvement=float(1-b[k]/std::max(b[0],1e-12));
        // Residual-to-time conversion is a model allowance, NOT a measured
        // jitter standard deviation. Keep a 25 ms floor for asynchronous phase.
        const double sensitivity=std::sqrt(validation_curve)/.01;
        out.residual_time_s=float(std::max(.025, std::sqrt(b[k])/sensitivity));
        out.ok=out.train_rms_rad<.05f && out.validation_rms_rad<.05f && out.improvement>=.05f && out.residual_time_s<=.060f;
        return out;
    }
private:
    void anchor(const V3& m,uint32_t t,unsigned segment) {
        have_=true; anchor_norm_=m.norm(); anchor_unit_=m/anchor_norm_; anchor_us_=t; segment_=segment;
    }
    History history_{};
    std::array<Scores,2> scores_{};
    V3 anchor_unit_=V3::Zero();
    float anchor_norm_=0;
    uint32_t anchor_us_=0, sequence_=0, last_observed_=0;
    unsigned segment_=0;
    bool have_=false, seen_=false;
};
} // namespace imu_cal
