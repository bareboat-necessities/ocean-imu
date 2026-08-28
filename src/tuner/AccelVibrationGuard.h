#pragma once

/*
  Copyright 2025-2026, Mikhail Grushinskiy
*/

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

namespace seastate::tuner {

/*
  Front-end guard against out-of-band accelerometer vibration.

  Motivation.  A sea-state estimator reads the accelerometer for two purposes at
  once: it is the gravity reference that fixes tilt, and it is the only
  measurement of wave acceleration.  Both live below about 1 Hz.  Machinery
  vibration -- an inboard auxiliary diesel is the case this was built for --
  arrives one to two decades higher, and on a sensor sampled at a few hundred Hz
  a good part of it arrives folded, because the crank orders run past Nyquist.

  Left alone that energy does not simply add zero-mean noise.  The attitude
  correction wobbles at the vibration frequencies, and because the measurement
  model is nonlinear in attitude the wobble rectifies into a *static* tilt
  offset, which then leaks gravity into the horizontal acceleration and biases
  displacement.  Measured on this repository's engine-noise study, 0.45 m/s^2 of
  recorded vibration buys about 3 degrees of standing tilt error with almost no
  fluctuation on it, and the same sign on every record.

  What this does.  A cascade of one-pole low passes on the three accelerometer
  axes, placed ahead of every consumer, so the vibration never reaches the
  attitude loop.  The separation is what makes it cheap: the signal to keep is
  below ~1 Hz and the energy to remove is above ~5 Hz, so a corner in the gap
  costs the wave band little.

  What it costs, and why it is gated.  Group delay, tau = poles / (2 pi fc),
  roughly flat across the wave band.  Acceleration delayed by tau yields
  displacement delayed by tau, an error of A * 2 pi f * tau on a wave of
  amplitude A -- proportional to wave amplitude, so the price is paid hardest in
  exactly the big seas where the estimator is most valuable.  Running the filter
  unconditionally would therefore charge every quiet passage for a problem it
  does not have.

  So the guard measures the out-of-band content it would remove and engages only
  when there is something to remove, ramping in over a slow time constant so
  there is no transient at the throttle.  Below the lower threshold the output
  is the input, bit for bit; the low pass still runs, because that is what
  produces the measurement the decision is made on.

  Choosing a corner.  Prefer the highest corner that removes the machinery, not
  the lowest corner that clears the waves: the stopband requirement is soft
  (rejection is traded against delay smoothly) while the delay cost is not.  At
  matched group delay a longer cascade rejects less of what matters here, since
  the damaging content sits just above the corner rather than deep in the
  stopband, which is why the default is two poles.
*/
class AccelVibrationGuard {
public:
    static constexpr int kMaxPoles = 4;

    // Corner of the detector high pass.  This is deliberately well above the
    // conditioning corner and is not the same decision: conditioning asks what
    // to remove, detection asks whether there is machinery present at all.  A
    // detector placed at the conditioning corner would measure the top of the
    // sea spectrum, which grows with wave height, and would then engage the
    // guard hardest in big seas -- exactly backwards.
    static constexpr float kDetectHzDefault = 25.0f;

    // Engagement band, in detector-band RMS.  Defaults are set from this
    // repository's measurements.  Over the eight stationary records the clean
    // deployed configuration reads 0.00796 to 0.00805 m/s^2 here -- a one
    // percent spread across a 31:1 range of significant wave height, which is
    // what a detector above the sea should look like -- and it is essentially
    // all accelerometer white noise.  With the engine running the same reading
    // is 0.037 at the quietest level swept, 0.087 at idle, and 0.144 at the
    // nominal cruise condition.
    //
    // The lower rail is therefore placed at about four times the clean floor,
    // and the band is nearer that floor than the midpoint because the two
    // errors are not symmetric: engaging spuriously costs delay in exactly the
    // big seas where the estimator matters most, while failing to engage at the
    // bottom of the vibration range costs little, since that is also where the
    // damage is smallest.
    //
    // These are absolute levels referenced to the accelerometer white noise
    // this repository injects (1.51e-3 g per axis).  A noisier part raises the
    // clean floor proportionally and wants setEngagement() called to match.
    static constexpr float kEngageLoDefault = 0.03f;
    static constexpr float kEngageHiDefault = 0.08f;
    static constexpr float kSlewTauDefault  = 5.0f;

    // fc_hz <= 0 disables the guard entirely.
    void setCutoffHz(float fc_hz) noexcept {
        if (!std::isfinite(fc_hz)) return;
        cutoff_hz_ = (fc_hz > 0.0f) ? fc_hz : 0.0f;
    }

    void setPoles(int poles) noexcept {
        poles_ = std::max(1, std::min(kMaxPoles, poles));
    }

    // Corner of the detector high pass.  Raise it to keep a lively sea out of
    // the decision; lower it to catch machinery whose orders run unusually low.
    void setDetectHz(float hz) noexcept {
        if (std::isfinite(hz) && hz > 0.0f) detect_hz_ = hz;
    }

    [[nodiscard]] float detectHz() const noexcept { return detect_hz_; }

    // Out-of-band RMS at which the guard starts and finishes engaging, and the
    // time constant it ramps over.  Setting lo >= hi makes engagement a step at
    // lo; setting lo <= 0 engages the guard unconditionally.
    void setEngagement(float lo_mps2, float hi_mps2, float slew_tau_sec) noexcept {
        if (std::isfinite(lo_mps2) && lo_mps2 >= 0.0f) engage_lo_ = lo_mps2;
        if (std::isfinite(hi_mps2) && hi_mps2 >= 0.0f) engage_hi_ = hi_mps2;
        if (std::isfinite(slew_tau_sec) && slew_tau_sec > 0.0f) {
            slew_tau_sec_ = slew_tau_sec;
        }
    }

    [[nodiscard]] bool enabled() const noexcept { return cutoff_hz_ > 0.0f; }
    [[nodiscard]] float cutoffHz() const noexcept { return cutoff_hz_; }
    [[nodiscard]] int poles() const noexcept { return poles_; }

    // How much of the low-passed signal is currently being used, in [0, 1].
    [[nodiscard]] float engagement() const noexcept { return weight_; }

    // Group delay actually being applied, i.e. the pass-band delay of the low
    // pass scaled by how far the guard is engaged.
    [[nodiscard]] float groupDelaySec() const noexcept {
        if (!enabled()) return 0.0f;
        return weight_ * static_cast<float>(poles_) /
               (2.0f * static_cast<float>(M_PI) * cutoff_hz_);
    }

    void reset() noexcept {
        for (auto& stage : stages_) stage.setZero();
        for (auto& stage : detect_stages_) stage.setZero();
        removed_ms_.setZero();
        weight_ = 0.0f;
        initialized_ = false;
    }

    // Returns the accelerometer with its out-of-band content removed, to the
    // extent the guard is engaged.  A disabled or unengaged guard returns the
    // input unchanged, bit for bit.
    Eigen::Vector3f step(const Eigen::Vector3f& acc, float dt) noexcept {
        if (!enabled() || !(dt > 0.0f) || !std::isfinite(dt)) return acc;
        if (!acc.allFinite()) return acc;

        if (!initialized_) {
            // Seed every stage at the first sample so the guard does not begin
            // by slewing up from zero through the whole gravity vector.
            for (int i = 0; i < kMaxPoles; ++i) stages_[i] = acc;
            // Only the first detector stage sees the raw signal.  Every stage
            // after it sees the previous stage's high-passed output, which is
            // centred on zero, so seeding those at acc would make them decay
            // through a whole gravity vector and spike the detector at startup.
            detect_stages_[0] = acc;
            for (int i = 1; i < kDetectPoles; ++i) detect_stages_[i].setZero();
            initialized_ = true;
            removed_ms_.setZero();
            weight_ = 0.0f;
            return acc;
        }

        const float alpha =
            std::exp(-2.0f * static_cast<float>(M_PI) * cutoff_hz_ * dt);
        Eigen::Vector3f low_passed = acc;
        for (int i = 0; i < poles_; ++i) {
            stages_[i] = (1.0f - alpha) * low_passed + alpha * stages_[i];
            low_passed = stages_[i];
        }

        // Detector: cascaded one-pole high passes on the raw input, so the
        // decision does not feed back on itself.  Two of them, because a single
        // pole -- or the tempting x minus low_passed, which is also only
        // first order near DC -- rolls off too slowly to keep the sea out of
        // the reading.
        const float gamma =
            std::exp(-2.0f * static_cast<float>(M_PI) * detect_hz_ * dt);
        Eigen::Vector3f high_passed = acc;
        for (int i = 0; i < kDetectPoles; ++i) {
            detect_stages_[i] = (1.0f - gamma) * high_passed + gamma * detect_stages_[i];
            high_passed -= detect_stages_[i];
        }

        const float beta =
            1.0f - std::exp(-2.0f * static_cast<float>(M_PI) * kRemovedRmsHz * dt);
        removed_ms_ += beta * (high_passed.cwiseAbs2() - removed_ms_);

        // Engage on measured out-of-band level, slewed so a throttle change
        // does not step the measurement path.
        const float rms = removedRms();
        float target;
        if (!(engage_lo_ > 0.0f)) {
            target = 1.0f;
        } else if (engage_hi_ > engage_lo_) {
            target = (rms - engage_lo_) / (engage_hi_ - engage_lo_);
        } else {
            target = (rms >= engage_lo_) ? 1.0f : 0.0f;
        }
        target = std::max(0.0f, std::min(1.0f, target));

        const float slew = 1.0f - std::exp(-dt / slew_tau_sec_);
        weight_ += slew * (target - weight_);
        // Park at the rails so a dormant guard is exactly transparent and a
        // fully engaged one is exactly the low pass, instead of asymptotically
        // near either.
        if (weight_ < kWeightEpsilon) weight_ = 0.0f;
        else if (weight_ > 1.0f - kWeightEpsilon) weight_ = 1.0f;

        if (weight_ <= 0.0f) return acc;
        if (weight_ >= 1.0f) return low_passed;
        return acc + weight_ * (low_passed - acc);
    }

    // Per-axis RMS of the detector band, smoothed.
    [[nodiscard]] Eigen::Vector3f removedRmsAxes() const noexcept {
        return removed_ms_.cwiseSqrt();
    }

    // Vector RMS in the detector band, in m/s^2.  Measured whenever the guard
    // has a cutoff, engaged or not, so it is readable as a vibration
    // diagnostic in its own right.
    [[nodiscard]] float removedRms() const noexcept {
        return std::sqrt(std::max(0.0f, removed_ms_.sum()));
    }

private:
    // Averaging corner for the out-of-band statistic: slow enough to be a
    // steady reading, fast enough to follow a throttle change.
    static constexpr float kRemovedRmsHz = 0.05f;
    static constexpr float kWeightEpsilon = 1e-4f;
    static constexpr int kDetectPoles = 2;

    Eigen::Vector3f stages_[kMaxPoles] = {
        Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero(),
        Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero(),
    };
    Eigen::Vector3f detect_stages_[kDetectPoles] = {
        Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero(),
    };
    Eigen::Vector3f removed_ms_ = Eigen::Vector3f::Zero();
    float cutoff_hz_ = 0.0f;
    float detect_hz_ = kDetectHzDefault;
    float engage_lo_ = kEngageLoDefault;
    float engage_hi_ = kEngageHiDefault;
    float slew_tau_sec_ = kSlewTauDefault;
    float weight_ = 0.0f;
    int poles_ = 2;
    bool initialized_ = false;
};

} // namespace seastate::tuner
