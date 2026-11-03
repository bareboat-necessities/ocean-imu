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
  sigma band (AdaptiveWaveBandPass) and the EW horizon of SeaStateAutoTuner.
  The class runs its own copy of that band and EW variance, started from the
  first sample instead of from rest: the shared band sees the input's DC level
  (turn-on bias, levelling) as a step at every start or reset, and that
  deterministic transient would otherwise read as wave evidence for a minute
  or two.  Away from a start the copy equals the shared statistic.  No
  estimator state enters, so the result keeps the tuner's exogeneity and its
  one-sample staging (the wrapper commits it at the top of the next sample).

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

  Period confidence.  The tuning period comes from WavePeriodEstimator's
  doubly integrated proxy, which weights the spectrum by 1/omega^4 and so is
  dominated, on small seas and in still water, by the sensor's slow errors
  and the front end's levelling transient; its period is then pulled far below
  the waves.  Acceleration is not weighted toward low frequency, so the band
  output decides: with one more pole at f_hi (y2 = b/(s + b) y, d = dy2/dt)
  its debiased zero-crossing frequency is

      f_c^2 = (V_d - V_d,0) / (V_y2 - V_y2,0) / (2 pi)^2,
      V_y2,0 = c b^3 / (4 (a + b)^2) + q b (a + 2b) / (4 a (a + b)^2) - q / (2 tau a^2),
      V_d,0  = c b^4 (2a + b) / (4 (a + b)^2) + q b^3 / (4 (a + b)^2).

  A band anchored on the period estimate centres the waves, f_c / sqrt(f_lo
  f_hi) of order one; a period pulled far below the waves puts them at the top
  of the band.  The ratio is meaningful only where the band is wave-dominated,
  so the period confidence is

      c_f = smoothstep((z - zf0)/(zf1 - zf0)) smoothstep((k1 - f_c/f_mid)/(k1 - k0)),

  zero without the ratio.

  All arithmetic is float: about forty operations and three square roots per
  sample, ten floats of state.
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
    // Detection statistic from which (zf0) and at which (zf1) the band is
    // wave-dominated enough to judge the period: zf1 = 16 is an excess of
    // eight times the noise variance at u = 0.5, a band 89 % waves.
    float zf0 = 8.0f;
    float zf1 = 16.0f;
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

    // Analysis helper (docs/carma-wave-excitation.md section 2.4): variance of
    // the same noise through WavePeriodEstimator's elevation proxy
    // s^2/(s + l)^4, l = 2 pi f_hp.  Not used by the filter.
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
        lp_ = band_ = m1_ = m2_ = w_ = 0.0f;
        y2_ = y2m_ = y2s_ = ds_ = 0.0f;
        var_ = NAN; var0_ = NAN; n_eff_ = NAN; excess_ = 0.0f; z_ = NAN; p_ = 0.0f;
        band_freq_ratio_ = NAN; cf_ = 0.0f;
    }

    void update(const Inputs& in, const CarmaSensorNoise& n, const CarmaEvidenceConfig& cfg) noexcept {
        if (!(in.dt > 0.0f) || !std::isfinite(in.dt) || !std::isfinite(in.accel)) return;
        if (!(std::isfinite(in.band_low_hz) && std::isfinite(in.band_high_hz) && in.band_high_hz > in.band_low_hz &&
              in.band_low_hz > 0.0f && std::isfinite(in.var_horizon_s) && in.var_horizon_s > 0.0f)) return;
        constexpr float k2Pi = 6.28318530717958647692f;
        const float dt = in.dt, x = in.accel;
        if (!started_) {
            started_ = true;   // from the first sample: no DC step
            lp_ = x;
        }
        // Band copy (AdaptiveWaveBandPass recursion at the corners in force),
        // one more pole at f_hi and its derivative.
        const float al = 1.0f - std::exp(-k2Pi * in.band_low_hz * dt);
        const float ah = 1.0f - std::exp(-k2Pi * in.band_high_hz * dt);
        const float lp_prev = lp_;
        lp_ = (1.0f - al) * lp_prev + al * x;
        band_ = (1.0f - ah) * band_ + ah * (1.0f - al) * (x - lp_prev);
        const float y2_prev = y2_;
        y2_ = (1.0f - ah) * y2_ + ah * band_;
        const float d = (y2_ - y2_prev) / dt;
        // EW moments at the tuner's horizon (SeaStateAutoTuner).
        const float alpha = 1.0f - std::exp(-dt / in.var_horizon_s);
        m1_ = (1.0f - alpha) * m1_ + alpha * band_;
        m2_ = (1.0f - alpha) * m2_ + alpha * band_ * band_;
        y2m_ = (1.0f - alpha) * y2m_ + alpha * y2_;
        y2s_ = (1.0f - alpha) * y2s_ + alpha * y2_ * y2_;
        ds_ = (1.0f - alpha) * ds_ + alpha * d * d;
        w_ = (1.0f - alpha) * w_ + alpha;
        if (!(w_ > 1e-6f) || !std::isfinite(in.band_white_gain)) return;

        // Detection.
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

        // Period confidence: band anchoring, judged where the band is waves.
        const float a = k2Pi * in.band_low_hz, b = k2Pi * in.band_high_hz;
        const float c = n.white_sigma_mps2 * n.white_sigma_mps2 * dt;
        const float q = n.bias_rw_mps2_sqrts * n.bias_rw_mps2_sqrts;
        const float ab2 = 4.0f * (a + b) * (a + b);
        const float vy0 = c * b * b * b / ab2 + q * b * (a + 2.0f * b) / (a * ab2) -
                          q / (2.0f * in.var_horizon_s * a * a);
        const float vd0 = c * b * b * b * b * (2.0f * a + b) / ab2 + q * b * b * b / ab2;
        const float my = y2m_ / w_;
        const float vy = std::max(0.0f, y2s_ / w_ - my * my) - std::max(vy0, 0.0f);
        const float vd = ds_ / w_ - vd0;
        band_freq_ratio_ = NAN;
        cf_ = 0.0f;
        if (vy > 0.0f && vd > 0.0f) {
            band_freq_ratio_ = std::sqrt(vd / vy) / (k2Pi * std::sqrt(in.band_low_hz * in.band_high_hz));
            cf_ = smoothstep((z_ - cfg.zf0) / std::max(cfg.zf1 - cfg.zf0, 1e-3f)) *
                  smoothstep((cfg.k1 - band_freq_ratio_) / std::max(cfg.k1 - cfg.k0, 1e-3f));
        }
    }

    // EW band variance, its expectation with no waves, effective sample size,
    // excess, detection statistic and evidence weight.
    float bandVariance() const noexcept { return var_; }
    float noiseVariance() const noexcept { return var0_; }
    float effectiveDof() const noexcept { return n_eff_; }
    float excessVariance() const noexcept { return excess_; }
    float zScore() const noexcept { return z_; }
    float waveProbability() const noexcept { return p_; }
    // Debiased band zero-crossing frequency over the band centre, and the
    // period confidence.
    float bandFrequencyRatio() const noexcept { return band_freq_ratio_; }
    float frequencyConfidence() const noexcept { return cf_; }

    // Excitation variance: max(E, 0), times p when gated.
    float excitationVariance(bool gated) const noexcept {
        const float e = std::max(excess_, 0.0f);
        return gated ? p_ * e : e;
    }

private:
    bool started_ = false;
    // Band copy and its EW mean and square; the extra pole's mean and
    // square, and its derivative's square.
    float lp_ = 0.0f, band_ = 0.0f, m1_ = 0.0f, m2_ = 0.0f, w_ = 0.0f;
    float y2_ = 0.0f, y2m_ = 0.0f, y2s_ = 0.0f, ds_ = 0.0f;
    // Outputs.
    float var_ = NAN;
    float var0_ = NAN;
    float n_eff_ = NAN;
    float excess_ = 0.0f;
    float z_ = NAN;
    float p_ = 0.0f;
    float band_freq_ratio_ = NAN;
    float cf_ = 0.0f;
};

}  // namespace ocean_imu::carma
