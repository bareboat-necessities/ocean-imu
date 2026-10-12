#pragma once

// Copyright 2026, Mikhail Grushinskiy
// Shared magnetic preprocessing. Optional calibrated relative lag is an
// effective epoch with uncertainty, never an absolute conversion timestamp.
#include "MagneticRotation.h"

#ifndef SEA_STATE_MAG_HOST_ALIGNMENT
#define SEA_STATE_MAG_HOST_ALIGNMENT 1
#endif

namespace ocean_imu::magnetic {
class Input {
public:
    void reset(float gyro_density, const M3& measurement_covariance) {
        history_.reset(); consistency_.reset(); diagnostic_ = {};
        snapshot_count_ = snapshot_next_ = 0; snapshot_ready_ = false;
        delay_us_ = 0; timing_allowance_s_ = 0;
        observation_ = {};
        observation_.covariance = measurement_covariance;
        observation_.uncertainty.gyro_density = gyro_density;
        have_source_ = have_submitted_ = new_source_ = source_valid_ = bad_epoch_ = false;
        source_sequence_ = submitted_sequence_ = 0;
        duplicates_ = history_misses_ = 0;
        age_ms_ = latency_ms_ = min_latency_ms_ = max_latency_ms_ = NAN;
        applied_ = false;
    }

    // Timing has the acquisition source's have/sequence/frame_us/read_end_us.
    // A non-Kalman observer can supply its bias estimate without pretending
    // to possess a bias covariance: pass uncertainty_known=false in that case.
    template <class Timing>
    void advance(uint32_t now, const V3& corrected_gyro, float bias_variance,
                 const Timing& timing, const V3& calibrated, const V3& hard_iron,
                 bool uncertainty_known = true, float gyro_step_variance = 0) {
        now_ = now; field_ = calibrated;
        uncertainty_known_ = uncertainty_known;
        const auto status = history_.push(now, corrected_gyro, bias_variance, gyro_step_variance);
        bad_epoch_ = status == Status::Duplicate || status == Status::OutOfOrder || status == Status::Invalid;
        if (bad_epoch_) history_.reset();
        if (bad_epoch_ || status == Status::Gap) {
            snapshot_count_ = snapshot_next_ = 0; snapshot_ready_ = false;
        }
        new_source_ = timing.have && (!have_source_ ||
            (static_cast<int32_t>(timing.sequence - source_sequence_) > 0 &&
             static_cast<int32_t>(timing.frame_us - source_us_) > 0)) &&
            static_cast<int32_t>(now - timing.frame_us) >= 0 &&
            static_cast<int32_t>(timing.read_end_us - timing.frame_us) >= 0;
        if (new_source_) {
            source_sequence_ = timing.sequence;
            source_us_ = timing.frame_us; read_end_us_ = timing.read_end_us;
            have_source_ = true;
            delay_us_ = timing.delay_us;
            timing_allowance_s_ = timing.timing_allowance_s;
            timing_valid_ = delay_us_ <= 70000 && std::isfinite(timing_allowance_s_) &&
                timing_allowance_s_ >= 0 && timing_allowance_s_ <= .0721f &&
                (delay_us_ == 0 || timing_allowance_s_ >= .027f);
            measurement_us_ = source_us_ - delay_us_;
            snapshot_ready_ = false;
            Uncertainty u = observation_.uncertainty;
            u.bias_variance = bias_variance;
            // Unqualified conversion aperture/completion phase + accepted-frame gap.
            // Advisory diagnostics only; never injected as a guessed delay.
            u.from_time_sigma = u.to_time_sigma = .028295f + 1.0f/30.0f +
                float(RotationHistory<>::MaxGapUs) * 1e-6f;
            if (delay_us_ && timing_valid_) u.from_time_sigma = u.to_time_sigma = timing_allowance_s_;
            diagnostic_ = consistency_.observe(history_, source_sequence_, measurement_us_,
                calibrated, hard_iron, observation_.covariance, u);
            if (!uncertainty_known_) {
                diagnostic_.nis = diagnostic_.rotation_signal_to_noise = NAN;
                diagnostic_.unexplained = false;
            }
        }
        source_valid_ = timing.have && have_source_ && timing.sequence == source_sequence_;
        observation_.uncertainty.bias_variance = bias_variance;
    }

    bool pending() const {
        return source_valid_ && (!have_submitted_ ||
            static_cast<int32_t>(source_sequence_ - submitted_sequence_) > 0);
    }
    template <class Gate>
    bool due(Gate& gate, bool healthy, uint32_t now_ms) {
        if (!SEA_STATE_MAG_HOST_ALIGNMENT) return gate.update(healthy, now_ms);
        if (!healthy) { gate.update(false, now_ms); return false; }
        if (!pending()) { ++duplicates_; return false; }
        return gate.update(true, now_ms);
    }
    // Retain the independent proxy itself. Never reconstruct its historical
    // attitude with a MEKF bias estimate. Called on every accepted IMU frame.
    void capture(const Q& proxy_bw, const Q& attitude_bw, const V3& accel, const V3& gyro) {
        if (bad_epoch_ || !proxy_bw.coeffs().allFinite() || !attitude_bw.coeffs().allFinite() ||
            !accel.allFinite() || !gyro.allFinite()) return;
        const Snapshot current{now_, proxy_bw, attitude_bw, accel, gyro};
        snapshots_[snapshot_next_] = current;
        snapshot_next_ = (snapshot_next_ + 1) % snapshots_.size();
        snapshot_count_ = std::min(snapshot_count_ + 1, unsigned(snapshots_.size()));
        if (!pending() || snapshot_ready_) return;
        if (measurement_us_ == now_) { installSnapshot(current); return; }
        const unsigned oldest = (snapshot_next_ + snapshots_.size() - snapshot_count_) % snapshots_.size();
        for (unsigned j=1; j<snapshot_count_; ++j) {
            const auto& a = snapshots_[(oldest+j-1)%snapshots_.size()];
            const auto& b = snapshots_[(oldest+j)%snapshots_.size()];
            if (int32_t(measurement_us_-a.time)<0 || int32_t(b.time-measurement_us_)<0) continue;
            const uint32_t dt=b.time-a.time;
            if (!dt || dt>RotationHistory<>::MaxGapUs) return;
            const float f=float(uint32_t(measurement_us_-a.time))/float(dt);
            installSnapshot(Snapshot{measurement_us_, a.proxy.slerp(f,b.proxy),
                a.attitude.slerp(f,b.attitude), a.accel+f*(b.accel-a.accel), a.gyro+f*(b.gyro-a.gyro)});
            return;
        }
    }
    bool prepare() {
        if (!source_valid_ || bad_epoch_ || !timing_valid_ ||
            history_.between(measurement_us_, now_, observation_.rotation) != Status::Ready) {
            ++history_misses_; return false;
        }
        observation_.uncertainty.from_time_sigma = float(uint32_t(read_end_us_ - source_us_)) * 1e-6f + timing_allowance_s_;
        observation_.uncertainty.to_time_sigma = 0;
        return true;
    }
    const Observation& observation() const { return observation_; }
    const Consistency& diagnostic() const { return diagnostic_; }
    bool newSource() const { return new_source_; }
    bool lastApplied() const { return applied_; }

    // Continuous feedback intentionally reuses a held observation, rotated to
    // this gyro epoch each time. This is NOT repeated discrete information.
    bool currentField(const V3& hard_iron, V3& field, M3* covariance = nullptr) {
        if (!SEA_STATE_MAG_HOST_ALIGNMENT) {
            field = field_ - hard_iron;
            if (covariance) *covariance = observation_.covariance;
            return field.allFinite();
        }
        if (!prepare()) return false;
        if (covariance) return align(field_, hard_iron, observation_.covariance,
            observation_.rotation, observation_.uncertainty, field, *covariance);
        field = observation_.rotation.q * (field_ - hard_iron);
        return field.allFinite();
    }
    template <class Fusion>
    void correct(Fusion& fusion, uint32_t correction_start_us) {
        applied_ = false;
        if (SEA_STATE_MAG_HOST_ALIGNMENT && (!pending() || !snapshot_ready_ || !prepare())) return;
        if (!SEA_STATE_MAG_HOST_ALIGNMENT) fusion.updateMag(field_);
        else fusion.updateMagTransported(field_, observation_);
        // A callback may consume the observation for reference/bias learning
        // even before a core correction is applied. Offer it exactly once.
        complete(correction_start_us, fusion.lastMagUpdateApplied());
    }
    void complete(uint32_t correction_start_us, bool applied) {
        applied_ = applied; submitted_sequence_ = source_sequence_; have_submitted_ = true;
        age_ms_ = float(uint32_t(now_ - source_us_)) * .001f;
        latency_ms_ = float(uint32_t(correction_start_us - read_end_us_)) * .001f;
        if (!std::isfinite(min_latency_ms_)) min_latency_ms_ = max_latency_ms_ = latency_ms_;
        min_latency_ms_ = std::min(min_latency_ms_, latency_ms_);
        max_latency_ms_ = std::max(max_latency_ms_, latency_ms_);
    }
    // Used for a normalized-vector backend: the caller supplies its unchanged
    // nominal unit-vector covariance expressed at the current field magnitude.
    void setCovariance(const M3& covariance) { observation_.covariance = covariance; }
    template <class Output>
    void print(Output& out, const char* service_name = "applied") const {
        out.printf("[MAGTIME] host_age_ms=%.3f read_to_update_ms=%.3f jitter_range_ms=%.3f "
                   "delay_ms=%.1f timing_allowance_ms=%.1f rot_pred_rad=%.5f rot_seen_rad=%.5f residual_uT=%.4f d2=%.3f snr=%.3f "
                   "%s=%u uncertainty_known=%u status=%u duplicates=%lu history_miss=%lu\n",
            double(age_ms_), double(latency_ms_), double(max_latency_ms_-min_latency_ms_),
            double(delay_us_)*.001, double(timing_allowance_s_)*1000,
            double(diagnostic_.predicted_angle_rad), double(diagnostic_.observed_angle_rad),
            double(diagnostic_.residual_uT), double(diagnostic_.nis), double(diagnostic_.rotation_signal_to_noise),
            service_name, unsigned(applied_), unsigned(uncertainty_known_), unsigned(diagnostic_.status),
            static_cast<unsigned long>(duplicates_), static_cast<unsigned long>(history_misses_));
    }
private:
    struct Snapshot { uint32_t time; Q proxy, attitude; V3 accel, gyro; };
    void installSnapshot(const Snapshot& s) {
        observation_.proxy_bw=s.proxy; observation_.attitude_bw=s.attitude;
        observation_.accel=s.accel; observation_.gyro=s.gyro; snapshot_ready_=true;
    }
    std::array<Snapshot,32> snapshots_{};
    unsigned snapshot_count_=0, snapshot_next_=0;
    bool snapshot_ready_=false, timing_valid_=true;
    uint32_t measurement_us_=0, delay_us_=0;
    float timing_allowance_s_=0;
    RotationHistory<> history_{};
    RotationConsistency consistency_{};
    Consistency diagnostic_{};
    Observation observation_{};
    V3 field_ = V3::Zero();
    uint32_t now_ = 0, source_us_ = 0, read_end_us_ = 0;
    uint32_t source_sequence_ = 0, submitted_sequence_ = 0;
    uint32_t duplicates_ = 0, history_misses_ = 0;
    bool have_source_ = false, source_valid_ = false, have_submitted_ = false;
    bool new_source_ = false, bad_epoch_ = false, applied_ = false, uncertainty_known_ = true;
    float age_ms_ = NAN, latency_ms_ = NAN, min_latency_ms_ = NAN, max_latency_ms_ = NAN;
};
} // namespace ocean_imu::magnetic
