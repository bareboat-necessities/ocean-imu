#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
  Released under the MIT License

  CarmaWaveEvidence -- EXPERIMENTAL measurement evidence for wave excitation
  in the three-axis CARMA research wrappers (SeaStateFusionFilter_CARMA3 and,
  through it, SeaStateFusionFilter_TFG_CARMA3).  Nothing in the shipping OU-III
  filter includes it.  Derivation, calibration and measurements:
  docs/carma-wave-excitation.md.

  Inputs are the levelled vertical acceleration of the measurement-only front
  end and the configuration that front end is in force with: the corners,
  sample interval and exact white-noise variance gain of the period-scaled
  sigma band (AdaptiveWaveBandPass), the EW horizon of SeaStateAutoTuner, and
  the corner and moment horizon of WavePeriodEstimator.  The class runs its
  own copy of that band and of the period estimator's elevation proxy, with
  their EW variances, started from the first sample instead of from rest: the
  shared filters see the input's DC level (turn-on bias, levelling) as a step
  at every start or reset, and that deterministic transient, metres in the
  doubly integrated proxy, would otherwise read as wave evidence for minutes.
  Away from a start the copies equal the shared statistics.  No estimator
  state enters, so the result keeps the tuner's exogeneity and its one-sample
  staging (the wrapper commits it at the top of the next sample).

  Sensor model (vertical, after levelling): white accelerometer noise of
  per-sample variance s^2 at the sample interval dt (two-sided density
  c = s^2 dt) plus a bias random walk of intensity q [(m/s^2)^2/s].  Through
  the band H(s) = s b / ((s + a)(s + b)), a = 2 pi f_lo, b = 2 pi f_hi, the
  expected EW variance under "no waves" is

      E[V | 0] = s^2 g_band + q b / (2 a (a + b)) - q / (2 tau a^2),

  where g_band is the band's exact time-varying white-noise gain, the second
  term the random-walk variance and the third the part the tuner's EW mean
  removes (tau = EW horizon).  The EW estimate has effective sample size

      N = 2 tau (a + b) (c + q')^2 / (c^2 + 2 c q' + q'^2 (a^2 + 3ab + b^2)/(ab)),
      q' = q / (a b),

  (Isserlis: Var V = V^2 pi / (tau B), B = (int S)^2 / int S^2), N = 2 tau (a + b)
  for white noise alone.  The excess E = V - E[V | 0] and

      z = E / (E[V | 0] sqrt(2/N + u^2)),

  u being the relative 1-sigma uncertainty of the sensor-noise variance
  itself, are the wave-detection statistic.  p = smoothstep((z - z0)/(z1 - z0))
  is the evidence weight, and p max(E, 0) the moment-matched excitation
  variance of the two-hypothesis (no waves / waves) predictive distribution.

  The same sensor model through WavePeriodEstimator's elevation proxy
  s^2 / (s + l)^4 (l = 2 pi f_hp) has variance V_eta,0 = c / (32 l^3) +
  q / (32 l^5).  r = 1 - V_eta,0 / measured is the wave fraction of the proxy
  whose moment ratio is the zero-crossing frequency; a fraction r biases that
  frequency toward the noise's own (about f_hp) by about 1 - sqrt(r), and
  c_f = smoothstep((r - r0)/(r1 - r0)) is its confidence.  The proxy can also
  carry what neither the noise model nor the band contains (below-band
  motion, the slow levelling transient of the front end's observer).  A wave
  the acceleration band sees has frequency f >= f_lo, so its elevation is at
  most p E+ / (|H(f_lo)|^2 (2 pi f_lo)^4) = p E+ / (G_lo w_lo^4) (only
  significant band excess counts); a proxy variance beyond twice
  V_eta,0 + p E+ / (G_lo w_lo^4) is not explained by the band, and c_f is
  scaled by smoothstep(2 - rho), rho the ratio of the two.

  That bound is loose when the band's low corner is low, and the doubly
  integrated proxy is dominated by any slow non-wave content (the observer's
  levelling transient over the first minutes, slow bias drift), so a tiny sea
  can carry a period estimate pulled to a few hundredths of a hertz while its
  acceleration sits at the top of the band.  Acceleration itself is not
  weighted toward low frequency, so the band output decides: with one more
  pole at f_hi (y2 = b/(s + b) y, d = dy2/dt) its debiased zero-crossing
  frequency is

      f_c^2 = (V_d - V_d,0) / (V_y2 - V_y2,0) / (2 pi)^2,
      V_y2,0 = c b^3 / (4 (a + b)^2) + q b (a + 2b) / (4 a (a + b)^2) - q / (2 tau a^2),
      V_d,0  = c b^4 (2a + b) / (4 (a + b)^2) + q b^3 / (4 (a + b)^2).

  A band anchored on the period estimate centres the waves, f_c / sqrt(f_lo
  f_hi) of order one; a period estimate pulled far below the waves puts them
  at the top, and c_f is scaled by smoothstep((k1 - f_c / f_mid)/(k1 - k0)).

  All arithmetic is float: about sixty operations and four square roots per
  sample, fifteen floats of state.
*/

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace ocean_imu::carma {

// How the wrapper turns the measured sea state into CARMA driving noise.
enum class CarmaExcitationMode : uint8_t {
    // Published baseline: sigma = c_carma max(0.05 m/s^2, configured band
    // floor, 1e-4, sigma_aw / c_sigma) with the legacy sigma channel (its
    // configured-floor subtraction and stillness attenuation).
    Legacy = 0,
    // Ablation: the legacy sigma channel with the floors removed (numerical
    // floor only).
    FloorOff = 1,
    // Noise-debiased band variance max(0, V - E[V | 0]) from the sensor
    // model; no evidence weight; frequency follows the tuner.
    Debiased = 2,
    // Debiased variance times the evidence weight p; the frequency moves
    // only with p c_f and is held otherwise.
    Evidence = 3,
};

// Accelerometer noise as seen by the levelled vertical channel.  The defaults
// are the simulated sensor of the validated bench configuration (1.51 mg white
// at 200 Hz, 5e-4 m/s^2/sqrt(s) bias random walk); a device needs its own
// values from a stationary Allan-variance record.
struct CarmaSensorNoise {
    float white_sigma_mps2   = 0.0148f;
    float bias_rw_mps2_sqrts = 5.0e-4f;
    // Relative 1-sigma uncertainty of the modelled noise variance.
    float rel_var_uncertainty = 0.5f;
};

struct CarmaEvidenceConfig {
    // Detection statistic thresholds: p = 0 below z0, 1 above z1.
    float z0 = 2.0f;
    float z1 = 4.0f;
    // Wave fraction of the elevation proxy at which the period estimate
    // starts (r0) and stops (r1) being down-weighted.  r1 = 0.9 bounds the
    // noise bias of a fully trusted period by about 5 %; r0 = 0.75 needs a
    // proxy variance four times the noise model's, which neither its
    // fluctuations nor a 2x noise mis-calibration reach.
    float r0 = 0.75f;
    float r1 = 0.90f;
    // Band zero-crossing frequency over the band centre at which the period
    // estimate starts (k0) and ends (k1) being distrusted.  Correctly
    // anchored seas, broadband and bimodal included, measure 0.9-1.9 on the
    // tuning seeds; a period estimate pulled low by non-wave content, 2.4-6.
    float k0 = 2.0f;
    float k1 = 2.5f;
};

class CarmaWaveEvidence {
public:
    struct Inputs {
        float dt = 0.0f;              // sample interval [s]
        float accel = NAN;            // levelled vertical acceleration [m/s^2]
        float band_low_hz = NAN;      // sigma-band corners in force
        float band_high_hz = NAN;
        float band_white_gain = NAN;  // its exact white-noise variance gain
        float var_horizon_s = NAN;    // tuner EW horizon
        float elev_high_pass_hz = NAN;// period-estimator leak corner
        float elev_horizon_s = NAN;   // and its moment horizon
    };

    static float smoothstep(float x) noexcept {
        if (!(x > 0.0f)) return 0.0f;
        if (!(x < 1.0f)) return 1.0f;
        return x * x * (3.0f - 2.0f * x);
    }

    // Expected EW band variance and effective sample size under "no waves".
    static void bandNoise(const CarmaSensorNoise& n, float dt, float f_lo, float f_hi, float white_gain,
                          float tau, float& var0, float& n_eff) noexcept {
        constexpr float k2Pi = 6.28318530717958647692f;
        const float s2 = n.white_sigma_mps2 * n.white_sigma_mps2;
        const float q = n.bias_rw_mps2_sqrts * n.bias_rw_mps2_sqrts;
        const float a = k2Pi * f_lo, b = k2Pi * f_hi;
        var0 = s2 * white_gain;
        n_eff = NAN;
        if (!(a > 0.0f) || !(b > a) || !(tau > 0.0f) || !(dt > 0.0f)) return;
        if (q > 0.0f) var0 += q * b / (2.0f * a * (a + b)) - q / (2.0f * tau * a * a);
        const float c = s2 * dt;
        const float qp = q / (a * b);
        const float num = (c + qp) * (c + qp);
        const float den = c * c + 2.0f * c * qp + qp * qp * (a * a + 3.0f * a * b + b * b) / (a * b);
        if (den > 0.0f) n_eff = 2.0f * tau * (a + b) * num / den;
    }

    // Elevation-proxy variance of the same noise through WavePeriodEstimator.
    static float elevationNoise(const CarmaSensorNoise& n, float dt, float f_hp) noexcept {
        constexpr float k2Pi = 6.28318530717958647692f;
        const float l = k2Pi * f_hp;
        if (!(l > 0.0f) || !(dt > 0.0f)) return NAN;
        const float c = n.white_sigma_mps2 * n.white_sigma_mps2 * dt;
        const float q = n.bias_rw_mps2_sqrts * n.bias_rw_mps2_sqrts;
        const float l3 = l * l * l;
        return c / (32.0f * l3) + q / (32.0f * l3 * l * l);
    }

    void reset() noexcept {
        started_ = false;
        lp_ = band_ = 0.0f;
        hp1_ = hp2_ = vel_ = elev_ = x_prev_ = 0.0f;
        m1_ = m2_ = w_ = e1_ = e2_ = ew_ = 0.0f;
        y2_ = y2m_ = y2s_ = ds_ = 0.0f;
        band_freq_ratio_ = NAN;
        var_ = NAN; var0_ = NAN; n_eff_ = NAN; excess_ = 0.0f; z_ = NAN; p_ = 0.0f;
        elev_var_ = NAN; elev_frac_ = NAN; elev_explained_ = NAN; cf_ = 0.0f;
    }

    void update(const Inputs& in, const CarmaSensorNoise& n, const CarmaEvidenceConfig& cfg) noexcept {
        if (!(in.dt > 0.0f) || !std::isfinite(in.dt) || !std::isfinite(in.accel)) return;
        constexpr float k2Pi = 6.28318530717958647692f;
        const float dt = in.dt, x = in.accel;
        if (!started_) {
            // Start from the first sample: the copies carry no DC step.
            started_ = true;
            lp_ = x;
            x_prev_ = x;
        }
        // Band copy (AdaptiveWaveBandPass recursion at the corners in force).
        if (std::isfinite(in.band_low_hz) && std::isfinite(in.band_high_hz) && in.band_high_hz > in.band_low_hz &&
            std::isfinite(in.var_horizon_s) && in.var_horizon_s > 0.0f) {
            const float al = 1.0f - std::exp(-k2Pi * in.band_low_hz * dt);
            const float ah = 1.0f - std::exp(-k2Pi * in.band_high_hz * dt);
            const float lp_prev = lp_;
            lp_ = (1.0f - al) * lp_prev + al * x;
            band_ = (1.0f - ah) * band_ + ah * (1.0f - al) * (x - lp_prev);
            const float y2_prev = y2_;
            y2_ = (1.0f - ah) * y2_ + ah * band_;
            const float d = (y2_ - y2_prev) / dt;
            // EW moments at the tuner's horizon (SeaStateAutoTuner).
            const float a = 1.0f - std::exp(-dt / in.var_horizon_s);
            m1_ = (1.0f - a) * m1_ + a * band_;
            m2_ = (1.0f - a) * m2_ + a * band_ * band_;
            y2m_ = (1.0f - a) * y2m_ + a * y2_;
            y2s_ = (1.0f - a) * y2s_ + a * y2_ * y2_;
            ds_ = (1.0f - a) * ds_ + a * d * d;
            w_ = (1.0f - a) * w_ + a;
        }
        // Elevation-proxy copy (WavePeriodEstimator recursion).
        if (std::isfinite(in.elev_high_pass_hz) && in.elev_high_pass_hz > 0.0f) {
            const float lam = k2Pi * in.elev_high_pass_hz;
            const float decay = std::exp(-lam * dt);
            const float gain = (1.0f - decay) / lam;
            const float s1 = decay * (hp1_ + x - x_prev_);
            const float s2 = decay * (hp2_ + s1 - hp1_);
            x_prev_ = x;
            hp1_ = s1;
            hp2_ = s2;
            vel_ = decay * vel_ + gain * s2;
            elev_ = decay * elev_ + gain * vel_;
            if (std::isfinite(in.elev_horizon_s) && in.elev_horizon_s > 0.0f) {
                const float a = 1.0f - std::exp(-dt / in.elev_horizon_s);
                e1_ = (1.0f - a) * e1_ + a * elev_;
                e2_ = (1.0f - a) * e2_ + a * elev_ * elev_;
                ew_ = (1.0f - a) * ew_ + a;
            }
        }
        if (w_ > 1e-6f && std::isfinite(in.band_white_gain)) {
            const float mu = m1_ / w_;
            var_ = std::max(0.0f, m2_ / w_ - mu * mu);
            float v0, ne;
            bandNoise(n, dt, in.band_low_hz, in.band_high_hz, in.band_white_gain, in.var_horizon_s, v0, ne);
            var0_ = std::max(v0, 0.0f);
            n_eff_ = ne;
            excess_ = var_ - var0_;
            const float rel2 = ((std::isfinite(ne) && ne > 0.0f) ? 2.0f / ne : 0.0f) +
                               n.rel_var_uncertainty * n.rel_var_uncertainty;
            const float sd = var0_ * std::sqrt(rel2);
            if (sd > 0.0f) {
                z_ = excess_ / sd;
                p_ = smoothstep((z_ - cfg.z0) / std::max(cfg.z1 - cfg.z0, 1e-3f));
            } else {
                // A noiseless model: any excess is evidence.
                z_ = excess_ > 0.0f ? INFINITY : 0.0f;
                p_ = excess_ > 0.0f ? 1.0f : 0.0f;
            }
        }
        // The period confidence needs settled proxy moments.
        if (ew_ > 0.5f) {
            const float mu = e1_ / ew_;
            elev_var_ = std::max(0.0f, e2_ / ew_ - mu * mu);
            const float e0 = elevationNoise(n, dt, in.elev_high_pass_hz);
            if (std::isfinite(e0) && elev_var_ > 0.0f) {
                elev_frac_ = 1.0f - e0 / elev_var_;
                cf_ = smoothstep((elev_frac_ - cfg.r0) / std::max(cfg.r1 - cfg.r0, 1e-3f));
                // Elevation the band can explain (file header).
                const float a = k2Pi * in.band_low_hz, b = k2Pi * in.band_high_hz;
                float explained = e0;
                if (a > 0.0f && b > a) {
                    const float g_lo = 0.5f * b * b / (a * a + b * b);
                    explained += p_ * std::max(excess_, 0.0f) / (g_lo * a * a * a * a);
                }
                elev_explained_ = elev_var_ / explained;
                cf_ *= smoothstep(2.0f - elev_explained_);
            }
        }
        // Acceleration-domain check of the band's anchoring (file header).
        band_freq_ratio_ = NAN;
        if (w_ > 1e-6f) {
            const float a = k2Pi * in.band_low_hz, b = k2Pi * in.band_high_hz;
            if (a > 0.0f && b > a) {
                const float c = n.white_sigma_mps2 * n.white_sigma_mps2 * dt;
                const float q = n.bias_rw_mps2_sqrts * n.bias_rw_mps2_sqrts;
                const float ab2 = 4.0f * (a + b) * (a + b);
                const float vy0 = c * b * b * b / ab2 + q * b * (a + 2.0f * b) / (a * ab2) -
                                  q / (2.0f * in.var_horizon_s * a * a);
                const float vd0 = c * b * b * b * b * (2.0f * a + b) / ab2 + q * b * b * b / ab2;
                const float mu = y2m_ / w_;
                const float vy = std::max(0.0f, y2s_ / w_ - mu * mu) - std::max(vy0, 0.0f);
                const float vd = ds_ / w_ - vd0;
                if (vy > 0.0f && vd > 0.0f) {
                    const float fc = std::sqrt(vd / vy) / k2Pi;
                    band_freq_ratio_ = fc / std::sqrt(in.band_low_hz * in.band_high_hz);
                    cf_ *= smoothstep((cfg.k1 - band_freq_ratio_) / std::max(cfg.k1 - cfg.k0, 1e-3f));
                }
            }
        }
    }

    // Expected band variance with no waves, effective sample size, excess,
    // detection statistic and evidence weight.
    float bandVariance() const noexcept { return var_; }
    float noiseVariance() const noexcept { return var0_; }
    float effectiveDof() const noexcept { return n_eff_; }
    float excessVariance() const noexcept { return excess_; }
    float zScore() const noexcept { return z_; }
    float waveProbability() const noexcept { return p_; }
    // Elevation-proxy variance, its wave fraction and the period confidence.
    float elevationVariance() const noexcept { return elev_var_; }
    float elevationWaveFraction() const noexcept { return elev_frac_; }
    // Measured proxy variance over the noise model plus what the band
    // explains (above 2: not wave evidence).
    float elevationUnexplainedRatio() const noexcept { return elev_explained_; }
    // Debiased band zero-crossing frequency over the band centre.
    float bandFrequencyRatio() const noexcept { return band_freq_ratio_; }
    float frequencyConfidence() const noexcept { return cf_; }

    // Excitation variance: max(E, 0), times p when gated.
    float excitationVariance(bool gated) const noexcept {
        const float e = std::max(excess_, 0.0f);
        return gated ? p_ * e : e;
    }

private:
    bool started_ = false;
    // Band copy and its EW moments.
    float lp_ = 0.0f, band_ = 0.0f, m1_ = 0.0f, m2_ = 0.0f, w_ = 0.0f;
    // Extra pole at f_hi, its EW mean and square, and its derivative's.
    float y2_ = 0.0f, y2m_ = 0.0f, y2s_ = 0.0f, ds_ = 0.0f;
    // Elevation-proxy copy and its EW moments.
    float hp1_ = 0.0f, hp2_ = 0.0f, vel_ = 0.0f, elev_ = 0.0f, x_prev_ = 0.0f;
    float e1_ = 0.0f, e2_ = 0.0f, ew_ = 0.0f;
    // Outputs.
    float var_ = NAN;
    float var0_ = NAN;
    float n_eff_ = NAN;
    float excess_ = 0.0f;
    float z_ = NAN;
    float p_ = 0.0f;
    float elev_var_ = NAN;
    float elev_frac_ = NAN;
    float elev_explained_ = NAN;
    float band_freq_ratio_ = NAN;
    float cf_ = 0.0f;
};

}  // namespace ocean_imu::carma
