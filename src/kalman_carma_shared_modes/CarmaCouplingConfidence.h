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

    Header-only, fixed size, no allocation; templated on the scalar T
    (float on the target, double in the analysis tools).
*/

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

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

template<typename T = float>
struct CouplingConfidenceConfig {
    bool  enabled = true;            // false: c = 1 always (the estimated coupling, unchecked)
    T monitor_horizon_s = T(60.0); // short estimate the prior is checked against
    T d2_high = T(16.0);           // inconsistent above (E[D^2] ~ 2 when consistent)
    T d2_low = T(6.0);             // consistent again below
    T beta_rel_floor = T(0.05);    // eps: relative regression tolerance
    T snr_low = T(10.0);           // band power / noise floor: no coupling at or below
    T snr_high = T(100.0);         // full coupling at or above (log-linear in between)
    T accel_noise_std = T(0.015);  // m/s^2 per axis (white), the noise-floor model
    T yaw_rate_dps = T(3.0);       // manoeuvre threshold on the low-passed world yaw rate
    T yaw_tau_s = T(2.0);
    T freq_rate_per_s = T(0.01);   // abrupt sea change: |d ln f_op/dt| threshold
    T freq_tau_s = T(10.0);
    T hold_s = T(30.0);            // evidence held off after a manoeuvre / sea change
    T nis_ratio_high = T(2.0);     // short/long NIS ratio: fallback above
    T nis_ratio_low = T(1.5);      // ... recovered below
    T nis_short_tau_s = T(20.0);
    T nis_long_tau_s = T(600.0);
    T tau_down_s = T(3.0);         // withdraw quickly
    T tau_up_s = T(60.0);          // restore slowly
};

// Band cross-spectral evidence of one estimate (S_hz and the needed autos).
template<typename T = float>
struct BandRegression {
    std::complex<T> beta[2]{};   // S_xz / S_zz, S_yz / S_zz
    T var[2] = {T(0.0), T(0.0)};     // large-sample variance of each component
    T szz = T(0.0);                // vertical band power
    T gamma2 = T(0.0);             // multiple coherence of z on (x, y)
    bool valid = false;
};

// re/im accessors of any Hermitian3-like type: S.re(i,j), S.im(i,j).
template<typename T, class H>
BandRegression<T> band_regression(const H& S, T n_eff) {
    BandRegression<T> r;
    const T szz = static_cast<T>(S.re(2,2));
    if (!(szz > T(0.0)) || !std::isfinite(szz)) return r;
    r.szz = szz;
    T num = T(0.0);
    for (int i = 0; i < 2; ++i) {
        const std::complex<T> s(static_cast<T>(S.re(i,2)), static_cast<T>(S.im(i,2)));
        const T sii = std::max(static_cast<T>(S.re(i,i)), T(0.0));
        r.beta[i] = s / szz;
        const T g2 = sii > T(0.0) ? std::min(T(1.0), std::norm(s) / (sii * szz)) : T(0.0);
        r.var[i] = n_eff > T(0.0) ? (T(1.0) - g2) * sii / (n_eff * szz) : T(0.0);
        num += std::norm(s);
    }
    const T sh = std::max(static_cast<T>(S.re(0,0)) + static_cast<T>(S.re(1,1)), T(0.0));
    r.gamma2 = sh > T(0.0) ? std::min(T(1.0), num / (sh * szz)) : T(0.0);
    r.valid = true;
    return r;
}

// Mismatch statistic D^2 of a prior regression against a monitor estimate.
template<typename T>
T mismatch_d2(const BandRegression<T>& prior, const BandRegression<T>& monitor, T rel_floor) {
    if (!prior.valid || !monitor.valid) return T(0.0);
    T d2 = T(0.0);
    for (int i = 0; i < 2; ++i) {
        const T floor = rel_floor * rel_floor * std::max(std::norm(prior.beta[i]), std::norm(monitor.beta[i]));
        const T v = prior.var[i] + monitor.var[i] + floor + T(1e-12);
        d2 += std::norm(monitor.beta[i] - prior.beta[i]) / v;
    }
    return d2;
}

// Direction of travel of the dominant wave component from a regression:
// beta = e^{i phi} r u (u a real unit 2-vector, NED x north, y east), the
// sign of u fixed by the quadrature lag of linear waves (phase of the
// along-direction cross-spectrum in (-180, 0) deg).  Degrees clockwise from
// the filter's north; uncertainty from the coherence and sample count.
template<typename T = float>
struct WaveDirectionEstimate {
    T to_deg = T(0.0);
    T uncertainty_deg = T(180.0);
    T phase_deg = T(0.0);          // x-z / y-z phase along the direction
    bool valid = false;
};

template<typename T>
WaveDirectionEstimate<T> wave_direction(const BandRegression<T>& r, T n_eff) {
    WaveDirectionEstimate<T> w;
    if (!r.valid) return w;
    const std::complex<T> q = r.beta[0] * r.beta[0] + r.beta[1] * r.beta[1];
    if (!(std::abs(q) > T(0.0))) return w;
    const T phi = T(0.5) * std::arg(q);
    const std::complex<T> rot = std::polar(T(1.0), -phi);
    T ux = (rot * r.beta[0]).real(), uy = (rot * r.beta[1]).real();
    T ph = phi;
    if (std::sin(ph) > T(0.0)) { ux = -ux; uy = -uy; ph -= T(3.14159265); }
    if (ph < -T(3.14159265)) ph += T(6.2831853);
    w.to_deg = std::atan2(uy, ux) * T(57.2957795);
    if (w.to_deg < T(0.0)) w.to_deg += T(360.0);
    w.phase_deg = ph * T(57.2957795);
    const T g2 = std::max(r.gamma2, T(1e-6));
    const T n = std::max(n_eff, T(1.0));
    w.uncertainty_deg = std::min(T(180.0), std::sqrt((T(1.0) - std::min(g2, T(0.999999))) / (T(2.0) * n * g2)) * T(57.2957795));
    w.valid = true;
    return w;
}

template<typename T = float>
struct CouplingConfidenceStatus {
    T confidence = T(0.0);         // c in [0, 1]
    T evidence = T(0.0);           // e (target of c)
    CouplingState state = CouplingState::WarmUp;
    T mismatch_d2 = T(0.0);
    T monitor_coherence = T(0.0);
    T snr = T(0.0);
    T nis_ratio = T(1.0);
    bool consistent = false;
    bool manoeuvre = false;
    bool sea_change = false;
};

template<typename T = float>
class CarmaCouplingConfidence {
public:
    using Config = CouplingConfidenceConfig<T>;
    using Status = CouplingConfidenceStatus<T>;
    using Regression = BandRegression<T>;

    void configure(const Config& c) { cfg_ = c; reset(); }
    const Config& config() const noexcept { return cfg_; }

    void reset() {
        st_ = Status{};
        st_.state = cfg_.enabled ? CouplingState::WarmUp : CouplingState::Disabled;
        st_.confidence = cfg_.enabled ? T(0.0) : T(1.0);
        yaw_lp_ = T(0); dlnf_lp_ = T(0); last_lnf_ = std::numeric_limits<T>::quiet_NaN();
        hold_man_ = hold_sea_ = T(0.0);
        nis_s_ = nis_l_ = std::numeric_limits<T>::quiet_NaN();
        nis_bad_ = false;
        consistent_ = false;
        forced_ = -T(1.0);
    }

    // Research / test hook: hold c at a value (negative: release).
    void force(T c) { forced_ = c; if (c >= T(0.0)) st_.confidence = std::min(T(1.0), c); }

    // One sample.  prior / monitor: regressions of the applied prior and of
    // the short monitor; monitor_full: the monitor has one horizon of data;
    // yaw_rate: world-frame yaw rate (rad/s); f_op: operating frequency (Hz);
    // nis: accelerometer NIS of this sample (NaN if none); dt in s.
    void update(T dt, const Regression& prior, const Regression& monitor, bool monitor_full,
                T yaw_rate, T f_op, T nis) {
        if (!(dt > T(0.0)) || !std::isfinite(dt)) return;
        if (!cfg_.enabled) { st_.confidence = T(1.0); st_.evidence = T(1.0); st_.state = CouplingState::Disabled; return; }
        // Manoeuvre and sea-change detectors.
        const T ay = -std::expm1(-dt / std::max(cfg_.yaw_tau_s, dt));
        if (std::isfinite(yaw_rate)) yaw_lp_ += ay * (yaw_rate - yaw_lp_);
        const bool man_now = std::fabs(yaw_lp_) * T(57.2957795) > cfg_.yaw_rate_dps;
        hold_man_ = man_now ? cfg_.hold_s : std::max(T(0.0), hold_man_ - dt);
        if (f_op > T(0.0) && std::isfinite(f_op)) {
            const T lnf = std::log(f_op);
            if (std::isfinite(last_lnf_)) {
                const T af = -std::expm1(-dt / std::max(cfg_.freq_tau_s, dt));
                dlnf_lp_ += af * ((lnf - last_lnf_) / dt - dlnf_lp_);
            }
            last_lnf_ = lnf;
        }
        const bool sea_now = std::fabs(dlnf_lp_) > cfg_.freq_rate_per_s;
        hold_sea_ = sea_now ? cfg_.hold_s : std::max(T(0.0), hold_sea_ - dt);
        st_.manoeuvre = hold_man_ > T(0.0);
        st_.sea_change = hold_sea_ > T(0.0);
        // Innovation consistency (self-normalized NIS).
        if (std::isfinite(nis) && nis >= T(0.0)) {
            const T as = -std::expm1(-dt / cfg_.nis_short_tau_s), al = -std::expm1(-dt / cfg_.nis_long_tau_s);
            nis_s_ = std::isfinite(nis_s_) ? nis_s_ + as * (nis - nis_s_) : nis;
            nis_l_ = std::isfinite(nis_l_) ? nis_l_ + al * (nis - nis_l_) : nis;
        }
        st_.nis_ratio = (std::isfinite(nis_s_) && nis_l_ > T(0.0)) ? nis_s_ / nis_l_ : T(1.0);
        if (st_.nis_ratio > cfg_.nis_ratio_high) nis_bad_ = true;
        else if (st_.nis_ratio < cfg_.nis_ratio_low) nis_bad_ = false;
        // Consistency of the prior with the monitor (hysteresis).
        st_.mismatch_d2 = mismatch_d2<T>(prior, monitor, cfg_.beta_rel_floor);
        if (!prior.valid || !monitor.valid || st_.mismatch_d2 > cfg_.d2_high) consistent_ = false;
        else if (st_.mismatch_d2 < cfg_.d2_low) consistent_ = true;
        st_.consistent = consistent_;
        st_.monitor_coherence = monitor.gamma2;
        // Wave energy against the accelerometer noise floor of the envelope
        // filter (sigma^2 dt pi b f / 2; b = band fraction).
        const T noise = cfg_.accel_noise_std * cfg_.accel_noise_std * dt * T(3.14159265) * band_frac_ *
                            std::max(f_op, T(1e-3)) * T(0.5);
        st_.snr = monitor.valid ? monitor.szz / std::max(noise, T(1e-30)) : T(0.0);
        const T e_snr = st_.snr <= cfg_.snr_low ? T(0.0)
                          : st_.snr >= cfg_.snr_high ? T(1.0)
                          : std::log(st_.snr / cfg_.snr_low) / std::log(cfg_.snr_high / cfg_.snr_low);
        T e = e_snr;
        CouplingState s = CouplingState::Coupled;
        if (!monitor_full || !prior.valid) { e = T(0.0); s = CouplingState::WarmUp; }
        else if (!consistent_ || nis_bad_) { e = T(0.0); s = CouplingState::Fallback; }
        else if (st_.manoeuvre || st_.sea_change) { e = T(0.0); s = CouplingState::Reduced; }
        else if (e < T(0.999)) s = CouplingState::Reduced;
        st_.evidence = e;
        st_.state = s;
        if (forced_ >= T(0.0)) { st_.confidence = std::min(T(1.0), forced_); return; }
        const T tau = e < st_.confidence ? cfg_.tau_down_s : cfg_.tau_up_s;
        st_.confidence += -std::expm1(-dt / std::max(tau, dt)) * (e - st_.confidence);
        st_.confidence = std::min(T(1.0), std::max(T(0.0), st_.confidence));
    }

    void set_band_fraction(T b) { if (b > T(0.0) && std::isfinite(b)) band_frac_ = b; }
    T confidence() const noexcept { return st_.confidence; }
    const Status& status() const noexcept { return st_; }

private:
    Config cfg_{};
    Status st_{};
    T band_frac_ = T(0.15);
    T yaw_lp_ = T(0), dlnf_lp_ = T(0), last_lnf_ = std::numeric_limits<T>::quiet_NaN();
    T hold_man_ = T(0.0), hold_sea_ = T(0.0);
    T nis_s_ = std::numeric_limits<T>::quiet_NaN(), nis_l_ = std::numeric_limits<T>::quiet_NaN();
    bool nis_bad_ = false, consistent_ = false;
    T forced_ = -T(1.0);
};

}  // namespace ocean_imu::carma_shared_modes
