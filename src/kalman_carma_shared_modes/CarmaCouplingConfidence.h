#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    CarmaCouplingConfidence -- the coupling-confidence state machine of the
    EXPERIMENTAL shared-modes QMEKF (KalmanCarmaSharedModes.h,
    docs/kalman-carma-shared-modes.md section 4).

    The shared directional prior is only as good as its cross-axis
    regression.  The research study found that a CONFIDENTLY wrong coupling
    (reversed direction, NEES 380) is far worse than no coupling, while the
    estimated coupling helps on every directional sea.  This class turns the
    available evidence into one number c in [0, 1] that scales the coupling
    the motion block applies (c = 1: the validated Option C prior; c = 0: the
    independent-axis prior, CARMA(4,1) in law on every axis):

      * consistency: the regression of the horizontal accelerations on the
        vertical one at the band, beta = S_hz / S_zz (complex, two-vector),
        of the coupling prior (the long 240 s estimate or an external one)
        against a short 60 s monitor estimate fed with the same world-frame
        accelerations,
            D^2 = sum_i |beta_f,i - beta_p,i|^2 / (v_f,i + v_p,i + (eps |beta_p,i|)^2),
            v_i = (1 - gamma_iz^2) S_ii / (N S_zz),
        (the large-sample variance of a transfer-function estimate from N
        independent envelope samples).  Two complex components: E[D^2] ~ 2
        when the prior is right.  Hysteresis: inconsistent above d2_high,
        consistent again below d2_low;
      * wave energy: the monitor's band power against the accelerometer
        noise floor sigma_n^2 dt pi b f / 2 of the same envelope filter;
      * manoeuvring: |world yaw rate| (low-passed) above a threshold, held
        for hold_s afterwards (the monitor must refill);
      * abrupt sea-state change: |d ln f_op / dt| (low-passed) above a
        threshold, held likewise;
      * innovation consistency: the 20 s mean accelerometer NIS over its own
        600 s mean (self-normalized: the configured R is inflated, so the
        absolute NIS level is not informative), with hysteresis;
      * warm-up: no coupling until the monitor holds one horizon of data.

    The evidence e = product of the gates is smoothed asymmetrically,
        c <- c + (1 - exp(-dt / tau)) (e - c),  tau = tau_down if e < c else tau_up,
    so the coupling is withdrawn within seconds and restored over a minute;
    with the hysteresis bands and hold times this prevents rapid switching.
    The caller applies c to the excitation at the next model commit and, when
    c falls, shrinks the unexplained horizontal-vertical motion
    cross-covariance by the same ratio (Kalman3D_Wave_CARMA3::
    shrink_motion_vertical_coupling), which keeps the covariance positive
    semidefinite and removes cross-covariance the reduced prior no longer
    sustains.  Nothing here touches the state: position, velocity and
    attitude are continuous through every change of c.

    Header-only, fixed size, no allocation; float on the target.
*/

#include <algorithm>
#include <cmath>
#include <complex>

namespace ocean_imu::carma_shared_modes {

enum class CouplingState : unsigned char {
    WarmUp,      // the monitor has not filled one horizon yet: c -> 0
    Coupled,     // all evidence present: c -> 1
    Reduced,     // some evidence missing (energy, manoeuvre, sea change): c -> fraction
    Fallback,    // the prior contradicts the data (consistency or innovations): c -> 0
    Disabled     // coupling off by configuration
};

inline const char* to_string(CouplingState s) {
    switch (s) {
        case CouplingState::WarmUp: return "warm-up";
        case CouplingState::Coupled: return "coupled";
        case CouplingState::Reduced: return "reduced";
        case CouplingState::Fallback: return "fallback";
        case CouplingState::Disabled: return "disabled";
    }
    return "?";
}

struct CouplingConfidenceConfig {
    bool  enabled = true;            // false: c = 1 always (the research SHM:est behaviour)
    float monitor_horizon_s = 60.0f; // short estimate the prior is checked against
    float d2_high = 16.0f;           // inconsistent above (E[D^2] ~ 2 when consistent)
    float d2_low = 6.0f;             // consistent again below
    float beta_rel_floor = 0.05f;    // eps: relative regression tolerance
    float snr_low = 10.0f;           // band power / noise floor: no coupling at or below
    float snr_high = 100.0f;         // full coupling at or above (log-linear in between)
    float accel_noise_std = 0.015f;  // m/s^2 per axis (white), the noise-floor model
    float yaw_rate_dps = 3.0f;       // manoeuvre threshold on the low-passed world yaw rate
    float yaw_tau_s = 2.0f;
    float freq_rate_per_s = 0.01f;   // abrupt sea change: |d ln f_op/dt| threshold
    float freq_tau_s = 10.0f;
    float hold_s = 30.0f;            // evidence held off after a manoeuvre / sea change
    float nis_ratio_high = 2.0f;     // short/long NIS ratio: fallback above
    float nis_ratio_low = 1.5f;      // ... recovered below
    float nis_short_tau_s = 20.0f;
    float nis_long_tau_s = 600.0f;
    float tau_down_s = 3.0f;         // withdraw quickly
    float tau_up_s = 60.0f;          // restore slowly
};

// Band cross-spectral evidence of one estimate (S_hz and the needed autos).
struct BandRegression {
    std::complex<float> beta[2]{};   // S_xz / S_zz, S_yz / S_zz
    float var[2] = {0.0f, 0.0f};     // large-sample variance of each component
    float szz = 0.0f;                // vertical band power
    float gamma2 = 0.0f;             // multiple coherence of z on (x, y)
    bool valid = false;
};

// re/im accessors of any Hermitian3-like type: S.re(i,j), S.im(i,j).
template<class H>
BandRegression band_regression(const H& S, float n_eff) {
    BandRegression r;
    const float szz = static_cast<float>(S.re(2,2));
    if (!(szz > 0.0f) || !std::isfinite(szz)) return r;
    r.szz = szz;
    float num = 0.0f;
    for (int i = 0; i < 2; ++i) {
        const std::complex<float> s(static_cast<float>(S.re(i,2)), static_cast<float>(S.im(i,2)));
        const float sii = std::max(static_cast<float>(S.re(i,i)), 0.0f);
        r.beta[i] = s / szz;
        const float g2 = sii > 0.0f ? std::min(1.0f, std::norm(s) / (sii * szz)) : 0.0f;
        r.var[i] = n_eff > 0.0f ? (1.0f - g2) * sii / (n_eff * szz) : 0.0f;
        num += std::norm(s);
    }
    const float sh = std::max(static_cast<float>(S.re(0,0)) + static_cast<float>(S.re(1,1)), 0.0f);
    r.gamma2 = sh > 0.0f ? std::min(1.0f, num / (sh * szz)) : 0.0f;
    r.valid = true;
    return r;
}

// Mismatch statistic D^2 of a prior regression against a monitor estimate.
inline float mismatch_d2(const BandRegression& prior, const BandRegression& monitor, float rel_floor) {
    if (!prior.valid || !monitor.valid) return 0.0f;
    float d2 = 0.0f;
    for (int i = 0; i < 2; ++i) {
        const float floor = rel_floor * rel_floor * std::max(std::norm(prior.beta[i]), std::norm(monitor.beta[i]));
        const float v = prior.var[i] + monitor.var[i] + floor + 1e-12f;
        d2 += std::norm(monitor.beta[i] - prior.beta[i]) / v;
    }
    return d2;
}

// Direction of travel of the dominant wave component from a regression:
// beta = e^{i phi} r u (u a real unit 2-vector, NED x north, y east), the
// sign of u fixed by the quadrature lag of linear waves (phase of the
// along-direction cross-spectrum in (-180, 0) deg).  Degrees clockwise from
// the filter's north; uncertainty from the coherence and sample count.
struct WaveDirectionEstimate {
    float to_deg = 0.0f;
    float uncertainty_deg = 180.0f;
    float phase_deg = 0.0f;          // x-z / y-z phase along the direction
    bool valid = false;
};

inline WaveDirectionEstimate wave_direction(const BandRegression& r, float n_eff) {
    WaveDirectionEstimate w;
    if (!r.valid) return w;
    const std::complex<float> q = r.beta[0] * r.beta[0] + r.beta[1] * r.beta[1];
    if (!(std::abs(q) > 0.0f)) return w;
    const float phi = 0.5f * std::arg(q);
    const std::complex<float> rot = std::polar(1.0f, -phi);
    float ux = (rot * r.beta[0]).real(), uy = (rot * r.beta[1]).real();
    float ph = phi;
    if (std::sin(ph) > 0.0f) { ux = -ux; uy = -uy; ph -= 3.14159265f; }
    if (ph < -3.14159265f) ph += 6.2831853f;
    w.to_deg = std::atan2(uy, ux) * 57.2957795f;
    if (w.to_deg < 0.0f) w.to_deg += 360.0f;
    w.phase_deg = ph * 57.2957795f;
    const float g2 = std::max(r.gamma2, 1e-6f);
    const float n = std::max(n_eff, 1.0f);
    w.uncertainty_deg = std::min(180.0f, std::sqrt((1.0f - std::min(g2, 0.999999f)) / (2.0f * n * g2)) * 57.2957795f);
    w.valid = true;
    return w;
}

struct CouplingConfidenceStatus {
    float confidence = 0.0f;         // c in [0, 1]
    float evidence = 0.0f;           // e (target of c)
    CouplingState state = CouplingState::WarmUp;
    float mismatch_d2 = 0.0f;
    float monitor_coherence = 0.0f;
    float snr = 0.0f;
    float nis_ratio = 1.0f;
    bool consistent = false;
    bool manoeuvre = false;
    bool sea_change = false;
};

class CarmaCouplingConfidence {
public:
    using Config = CouplingConfidenceConfig;

    void configure(const Config& c) { cfg_ = c; reset(); }
    const Config& config() const noexcept { return cfg_; }

    void reset() {
        st_ = CouplingConfidenceStatus{};
        st_.state = cfg_.enabled ? CouplingState::WarmUp : CouplingState::Disabled;
        st_.confidence = cfg_.enabled ? 0.0f : 1.0f;
        yaw_lp_ = 0.0f; dlnf_lp_ = 0.0f; last_lnf_ = NAN;
        hold_man_ = hold_sea_ = 0.0f;
        nis_s_ = nis_l_ = NAN;
        nis_bad_ = false;
        consistent_ = false;
        forced_ = -1.0f;
    }

    // Research / test hook: hold c at a value (negative: release).
    void force(float c) { forced_ = c; if (c >= 0.0f) st_.confidence = std::min(1.0f, c); }

    // One sample.  prior / monitor: regressions of the applied prior and of
    // the short monitor; monitor_full: the monitor has one horizon of data;
    // yaw_rate: world-frame yaw rate (rad/s); f_op: operating frequency (Hz);
    // nis: accelerometer NIS of this sample (NaN if none); dt in s.
    void update(float dt, const BandRegression& prior, const BandRegression& monitor, bool monitor_full,
                float yaw_rate, float f_op, float nis) {
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;
        if (!cfg_.enabled) { st_.confidence = 1.0f; st_.evidence = 1.0f; st_.state = CouplingState::Disabled; return; }
        // Manoeuvre and sea-change detectors.
        const float ay = -std::expm1(-dt / std::max(cfg_.yaw_tau_s, dt));
        if (std::isfinite(yaw_rate)) yaw_lp_ += ay * (yaw_rate - yaw_lp_);
        const bool man_now = std::fabs(yaw_lp_) * 57.2957795f > cfg_.yaw_rate_dps;
        hold_man_ = man_now ? cfg_.hold_s : std::max(0.0f, hold_man_ - dt);
        if (f_op > 0.0f && std::isfinite(f_op)) {
            const float lnf = std::log(f_op);
            if (std::isfinite(last_lnf_)) {
                const float af = -std::expm1(-dt / std::max(cfg_.freq_tau_s, dt));
                dlnf_lp_ += af * ((lnf - last_lnf_) / dt - dlnf_lp_);
            }
            last_lnf_ = lnf;
        }
        const bool sea_now = std::fabs(dlnf_lp_) > cfg_.freq_rate_per_s;
        hold_sea_ = sea_now ? cfg_.hold_s : std::max(0.0f, hold_sea_ - dt);
        st_.manoeuvre = hold_man_ > 0.0f;
        st_.sea_change = hold_sea_ > 0.0f;
        // Innovation consistency (self-normalized NIS).
        if (std::isfinite(nis) && nis >= 0.0f) {
            const float as = -std::expm1(-dt / cfg_.nis_short_tau_s), al = -std::expm1(-dt / cfg_.nis_long_tau_s);
            nis_s_ = std::isfinite(nis_s_) ? nis_s_ + as * (nis - nis_s_) : nis;
            nis_l_ = std::isfinite(nis_l_) ? nis_l_ + al * (nis - nis_l_) : nis;
        }
        st_.nis_ratio = (std::isfinite(nis_s_) && nis_l_ > 0.0f) ? nis_s_ / nis_l_ : 1.0f;
        if (st_.nis_ratio > cfg_.nis_ratio_high) nis_bad_ = true;
        else if (st_.nis_ratio < cfg_.nis_ratio_low) nis_bad_ = false;
        // Consistency of the prior with the monitor (hysteresis).
        st_.mismatch_d2 = mismatch_d2(prior, monitor, cfg_.beta_rel_floor);
        if (!prior.valid || !monitor.valid) consistent_ = false;
        else if (st_.mismatch_d2 > cfg_.d2_high) consistent_ = false;
        else if (st_.mismatch_d2 < cfg_.d2_low) consistent_ = true;
        st_.consistent = consistent_;
        st_.monitor_coherence = monitor.gamma2;
        // Wave energy against the accelerometer noise floor of the envelope
        // filter (sigma^2 dt pi b f / 2; b = band fraction).
        const float noise = cfg_.accel_noise_std * cfg_.accel_noise_std * dt * 3.14159265f * band_frac_ *
                            std::max(f_op, 1e-3f) * 0.5f;
        st_.snr = monitor.valid ? monitor.szz / std::max(noise, 1e-30f) : 0.0f;
        const float e_snr = st_.snr <= cfg_.snr_low ? 0.0f
                          : st_.snr >= cfg_.snr_high ? 1.0f
                          : std::log(st_.snr / cfg_.snr_low) / std::log(cfg_.snr_high / cfg_.snr_low);
        float e = e_snr;
        CouplingState s = CouplingState::Coupled;
        if (!monitor_full || !prior.valid) { e = 0.0f; s = CouplingState::WarmUp; }
        else if (!consistent_ || nis_bad_) { e = 0.0f; s = CouplingState::Fallback; }
        else if (st_.manoeuvre || st_.sea_change) { e = 0.0f; s = CouplingState::Reduced; }
        else if (e < 0.999f) s = CouplingState::Reduced;
        st_.evidence = e;
        st_.state = s;
        if (forced_ >= 0.0f) { st_.confidence = std::min(1.0f, forced_); return; }
        const float tau = e < st_.confidence ? cfg_.tau_down_s : cfg_.tau_up_s;
        st_.confidence += -std::expm1(-dt / std::max(tau, dt)) * (e - st_.confidence);
        st_.confidence = std::min(1.0f, std::max(0.0f, st_.confidence));
    }

    void set_band_fraction(float b) { if (b > 0.0f && std::isfinite(b)) band_frac_ = b; }
    float confidence() const noexcept { return st_.confidence; }
    const CouplingConfidenceStatus& status() const noexcept { return st_; }

private:
    Config cfg_{};
    CouplingConfidenceStatus st_{};
    float band_frac_ = 0.15f;
    float yaw_lp_ = 0.0f, dlnf_lp_ = 0.0f, last_lnf_ = NAN;
    float hold_man_ = 0.0f, hold_sea_ = 0.0f;
    float nis_s_ = NAN, nis_l_ = NAN;
    bool nis_bad_ = false, consistent_ = false;
    float forced_ = -1.0f;
};

}  // namespace ocean_imu::carma_shared_modes
