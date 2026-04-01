#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class SeaMetricsFromSpectrum {
public:
    enum class SpectrumType {
        Swell,
        Mixed,
        WindSea,
        Unknown
    };

    struct WaveCountEstimate {
        double expected   = 0.0;
        double ci_lower   = 0.0;
        double ci_upper   = 0.0;
        double confidence = 0.95;
    };

    struct Report {
        // Inputs / bookkeeping
        bool   valid = false;
        double lowfreq_guard_hz = 0.0;
        double observation_duration_s = 600.0;

        // Raw spectral moments in Hz-domain:
        //   m_n = ∫ f^n S_eta(f) df
        double m_neg1 = 0.0;
        double m0     = 0.0;
        double m1     = 0.0;
        double m2     = 0.0;
        double m3     = 0.0;
        double m4     = 0.0;

        // Central moments of the frequency-distribution induced by S(f)
        double mu2 = 0.0;
        double mu3 = 0.0;
        double mu4 = 0.0;

        // Peak / mean frequency
        double fp_hz  = 0.0;
        double fp_rad = 0.0;
        double tp_s   = 0.0;
        double mean_frequency_hz  = 0.0;
        double mean_frequency_rad = 0.0;

        // Standard periods
        double t1_s    = 0.0; // Tm01
        double tm02_s  = 0.0; // Tz
        double te_s    = 0.0; // Tm-10 / energy period
        double tm0m1_s = 0.0; // alias of Te
        double tm1m1_s = 0.0; // heuristic ratio-based period
        double mean_group_period_s = 0.0;

        // Heights
        double rms_displacement_m = 0.0;
        double hs_m               = 0.0; // Hm0 = 4*sqrt(m0)
        double h1_10_crest_m      = 0.0;
        double h1_100_crest_m     = 0.0;
        double most_probable_max_crest_m = 0.0;
        double expected_max_crest_m      = 0.0;

        // Rates / counts
        double upcrossing_rate_hz   = 0.0;
        double downcrossing_rate_hz = 0.0;
        WaveCountEstimate wave_count;

        // Bandwidth / regularity
        double rbw              = 0.0;
        double regularity_spec  = 0.0;
        double regularity_phase = std::numeric_limits<double>::quiet_NaN(); // not identifiable from spectrum alone
        double rbw_phase_increment = std::numeric_limits<double>::quiet_NaN(); // not identifiable from spectrum alone
        double narrowness_nu    = 0.0;

        double bandwidth_clh    = 0.0;
        double bandwidth_goda   = 0.0;
        double bandwidth_kuik   = 0.0;
        double bandwidth_epsilon = 0.0;
        double longuet_higgins_width = 0.0;
        double spectral_bandwidth_hz = 0.0;
        double spectral_narrowness_ratio = 0.0;

        // Shape / peakedness
        double spectral_skewness = 0.0;
        double spectral_kurtosis = 0.0;
        double spectral_excess_kurtosis = 0.0;
        double peakedness_ochi_q = 0.0;
        double benassai_parameter = 0.0;
        double peak_enhancement_gamma = 0.0;
        double peakedness_goda = 0.0;
        SpectrumType spectrum_type = SpectrumType::Unknown;

        // Steepness / nonlinearity / development
        double deep_water_wavelength_tp_m = 0.0;
        double wave_steepness = 0.0;
        double nonlinearity_parameter = 0.0;
        double ursell_number = std::numeric_limits<double>::quiet_NaN();
        double wave_age = std::numeric_limits<double>::quiet_NaN();

        // Extremes / probabilities
        double crest_exceed_prob_hs = 0.0;          // P(crest > Hs)
        double crest_exceed_prob_tayfun_hs = 0.0;   // heuristic
        double pot_mean_excess_m = 0.0;
        double weibull_return_height_m = 0.0;

        // Groupiness / instability
        double groupiness_factor = 0.0;
        double bfi = 0.0;
        double mean_group_length = 0.0;
        double group_height_factor = 0.0;
        double expected_run_length_above_hs = 0.0;

        // Energy / power
        double wave_energy_density_j_m2 = 0.0;
        double wave_power_w_m = 0.0;
        double energy_flux_period_s = 0.0;

        // Breaking / shallow water
        double breaking_probability = 0.0;
        double depth_limited_breaking_height_m = 0.0;
        double bottom_orbital_velocity_mps = std::numeric_limits<double>::quiet_NaN();

        // Radiation stress proxies
        double radiation_stress_xx_n_m = 0.0;
        double radiation_stress_yy_n_m = 0.0;
        double radiation_stress_xy_n_m = 0.0;

        // Comfort / sickness heuristics
        double msdv = 0.0;
        double seasickness_incidence_pct = 0.0;
        double motion_comfort_level_0_100 = 100.0;
        double vertical_motion_intensity = 0.0;
        double motion_character = 0.0;
        double time_to_onset_min = std::numeric_limits<double>::infinity();

        // Quality / stability
        double snr_db = std::numeric_limits<double>::quiet_NaN();
        double temporal_stability_hs = std::numeric_limits<double>::quiet_NaN();
        double temporal_stability_tp = std::numeric_limits<double>::quiet_NaN();
        double data_quality_0_1 = 0.0;

        // Application-specific
        double wec_capture_width_ratio = std::numeric_limits<double>::quiet_NaN();
        double sea_swell_partition = std::numeric_limits<double>::quiet_NaN();
        double wave_age_class = std::numeric_limits<double>::quiet_NaN();
    };

    SeaMetricsFromSpectrum(bool enable_extended = true,
                           bool enable_negative = true,
                           double observation_duration_s = 600.0,
                           double depth_m = std::numeric_limits<double>::quiet_NaN(),
                           double wind_speed_mps = std::numeric_limits<double>::quiet_NaN(),
                           double exposure_hours = 1.0,
                           double susceptibility = 0.5)
        : extended_metrics_(enable_extended),
          negative_moments_(enable_negative),
          observation_duration_s_(std::max(1.0, observation_duration_s)),
          depth_m_(depth_m),
          wind_speed_mps_(wind_speed_mps),
          exposure_hours_(std::max(0.0, exposure_hours)),
          susceptibility_(std::clamp(susceptibility, 0.0, 1.0)) {
        reset();
    }

    void reset() {
        report_ = Report{};
        report_.observation_duration_s = observation_duration_s_;

        has_history_ = false;
        hs_ema_ = tp_ema_ = 0.0;
        hs_var_ = tp_var_ = 0.0;
    }

    void setObservationDuration(double seconds) {
        observation_duration_s_ = std::max(1.0, seconds);
        report_.observation_duration_s = observation_duration_s_;
    }

    void setDepthMeters(double depth_m) {
        depth_m_ = depth_m;
    }

    void setWindSpeedMps(double wind_speed_mps) {
        wind_speed_mps_ = wind_speed_mps;
    }

    void setComfortExposureHours(double exposure_hours) {
        exposure_hours_ = std::max(0.0, exposure_hours);
    }

    void setComfortSusceptibility(double susceptibility) {
        susceptibility_ = std::clamp(susceptibility, 0.0, 1.0);
    }

    void setWECExtractedPowerPerMeter(double extracted_w_m) {
        wec_extracted_w_m_ = extracted_w_m;
    }

    template <std::size_t N>
    void updateFromSpectrum(const std::array<double, N>& freqs_hz,
                            const std::array<double, N>& S_eta_m2_per_hz,
                            double lowfreq_guard_hz = 0.0,
                            double dt_s = 0.0) {
        std::vector<double> f(freqs_hz.begin(), freqs_hz.end());
        std::vector<double> s(S_eta_m2_per_hz.begin(), S_eta_m2_per_hz.end());
        updateFromSpectrum(f, s, lowfreq_guard_hz, dt_s);
    }

    void updateFromSpectrum(const std::vector<double>& freqs_hz,
                            const std::vector<double>& S_eta_m2_per_hz,
                            double lowfreq_guard_hz = 0.0,
                            double dt_s = 0.0) {
        resetReportOnly_();

        const std::size_t N = freqs_hz.size();
        if (N == 0 || N != S_eta_m2_per_hz.size()) {
            return;
        }

        last_freqs_hz_ = freqs_hz;
        last_spectrum_ = S_eta_m2_per_hz;
        last_df_hz_ = computeDfFromCenters_(freqs_hz);

        report_.lowfreq_guard_hz = std::max(0.0, lowfreq_guard_hz);

        computeMoments_(freqs_hz, S_eta_m2_per_hz, last_df_hz_, report_.lowfreq_guard_hz);
        if (!(report_.m0 > EPS_)) {
            return;
        }

        computePeakMetrics_(freqs_hz, S_eta_m2_per_hz, report_.lowfreq_guard_hz);
        computeMomentDerivedMetrics_();
        computeHeightsAndExtremes_();
        computeBandwidthAndShape_();
        computeEnergetics_();
        computeNonlinearAndBreaking_();
        computeComfortMetrics_();
        computeQualityMetrics_(freqs_hz, S_eta_m2_per_hz);
        updateStability_(dt_s);

        report_.valid = true;
    }

    const Report& report() const {
        return report_;
    }

private:
    static constexpr double EPS_ = 1e-12;
    static constexpr double G_ = 9.80665;
    static constexpr double RHO_ = 1025.0;
    static constexpr double BETA_SPEC_ = 1.0;

    bool extended_metrics_ = true;
    bool negative_moments_ = true;

    double observation_duration_s_ = 600.0;
    double depth_m_ = std::numeric_limits<double>::quiet_NaN();
    double wind_speed_mps_ = std::numeric_limits<double>::quiet_NaN();
    double exposure_hours_ = 1.0;
    double susceptibility_ = 0.5;
    double wec_extracted_w_m_ = std::numeric_limits<double>::quiet_NaN();

    Report report_{};

    std::vector<double> last_freqs_hz_;
    std::vector<double> last_spectrum_;
    std::vector<double> last_df_hz_;

    bool has_history_ = false;
    double hs_ema_ = 0.0;
    double tp_ema_ = 0.0;
    double hs_var_ = 0.0;
    double tp_var_ = 0.0;

    static double twoPi_() { return 2.0 * M_PI; }
    static double safeSqrt_(double x) { return (x > 0.0) ? std::sqrt(x) : 0.0; }

    void resetReportOnly_() {
        const double obs = observation_duration_s_;
        report_ = Report{};
        report_.observation_duration_s = obs;
    }

    static std::vector<double> computeDfFromCenters_(const std::vector<double>& f) {
        const std::size_t N = f.size();
        std::vector<double> df(N, 0.0);
        if (N == 0) return df;
        if (N == 1) {
            df[0] = std::max(0.5 * f[0], 0.0);
            return df;
        }

        std::vector<double> edges(N + 1, 0.0);

        edges[0] = f[0] * std::sqrt(std::max(f[0], EPS_) / std::max(f[1], EPS_));
        for (std::size_t i = 1; i < N; ++i) {
            edges[i] = std::sqrt(std::max(f[i - 1], EPS_) * std::max(f[i], EPS_));
        }
        edges[N] = f[N - 1] * std::sqrt(std::max(f[N - 1], EPS_) / std::max(f[N - 2], EPS_));

        edges[0] = std::max(0.0, edges[0]);
        for (std::size_t i = 0; i < N; ++i) {
            if (!(edges[i + 1] > edges[i])) {
                edges[i + 1] = edges[i] + std::max(1e-6, 0.1 * std::max(f[i], EPS_));
            }
            df[i] = edges[i + 1] - edges[i];
        }
        return df;
    }

    void computeMoments_(const std::vector<double>& f_hz,
                         const std::vector<double>& S,
                         const std::vector<double>& df_hz,
                         double f_guard_hz) {
        double m_neg1 = 0.0, m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0;

        for (std::size_t i = 0; i < f_hz.size(); ++i) {
            const double f = f_hz[i];
            const double s = std::max(0.0, S[i]);
            const double df = std::max(0.0, df_hz[i]);

            if (!(f > 0.0) || !(df > 0.0) || f < f_guard_hz) continue;

            const double sf = s * df;
            m0 += sf;
            m1 += sf * f;
            m2 += sf * f * f;

            if (extended_metrics_) {
                m3 += sf * f * f * f;
                m4 += sf * f * f * f * f;
            }
            if (negative_moments_) {
                m_neg1 += sf / std::max(f, EPS_);
            }
        }

        report_.m0 = m0;
        report_.m1 = m1;
        report_.m2 = m2;
        report_.m3 = m3;
        report_.m4 = m4;
        report_.m_neg1 = m_neg1;
    }

    static double peakInterpolateHz_(const std::vector<double>& freqs_hz,
                                     const std::vector<double>& S,
                                     std::size_t k) {
        if (k == 0 || k + 1 >= freqs_hz.size()) return freqs_hz[k];

        const double f0 = freqs_hz[k - 1];
        const double f1 = freqs_hz[k];
        const double f2 = freqs_hz[k + 1];
        if (!(f0 > 0.0 && f1 > 0.0 && f2 > 0.0)) return f1;

        const double x0 = std::log(f0);
        const double x1 = std::log(f1);
        const double x2 = std::log(f2);

        const double y0 = std::log(std::max(S[k - 1], EPS_));
        const double y1 = std::log(std::max(S[k], EPS_));
        const double y2 = std::log(std::max(S[k + 1], EPS_));

        const double L0a = 1.0 / ((x0 - x1) * (x0 - x2));
        const double L1a = 1.0 / ((x1 - x0) * (x1 - x2));
        const double L2a = 1.0 / ((x2 - x0) * (x2 - x1));

        const double a = y0 * L0a + y1 * L1a + y2 * L2a;
        const double b =
            y0 * (-(x1 + x2) * L0a) +
            y1 * (-(x0 + x2) * L1a) +
            y2 * (-(x0 + x1) * L2a);

        if (a >= 0.0) return f1;

        const double x_peak = -b / (2.0 * a);
        const double f_peak = std::exp(x_peak);

        if (f_peak <= std::min(f0, f2) || f_peak >= std::max(f0, f2)) return f1;
        return f_peak;
    }

    void computePeakMetrics_(const std::vector<double>& freqs_hz,
                             const std::vector<double>& S,
                             double f_guard_hz) {
        std::size_t k0 = 0;
        while (k0 + 1 < freqs_hz.size() && freqs_hz[k0 + 1] <= std::max(f_guard_hz, 1.05 * freqs_hz.front())) {
            ++k0;
        }

        std::size_t k = k0;
        double vmax = -1.0;
        for (std::size_t i = k0; i < S.size(); ++i) {
            const double s = std::max(0.0, S[i]);
            if (s > vmax) {
                vmax = s;
                k = i;
            }
        }

        const double fp_hz = (vmax > 0.0) ? peakInterpolateHz_(freqs_hz, S, k) : 0.0;
        report_.fp_hz = fp_hz;
        report_.fp_rad = twoPi_() * fp_hz;
        report_.tp_s = (fp_hz > EPS_) ? (1.0 / fp_hz) : 0.0;
    }

    void computeMomentDerivedMetrics_() {
        const double m0 = report_.m0;
        const double m1 = report_.m1;
        const double m2 = report_.m2;
        const double m3 = report_.m3;
        const double m4 = report_.m4;
        const double mn1 = report_.m_neg1;

        report_.mean_frequency_hz  = (m0 > EPS_) ? (m1 / m0) : 0.0;
        report_.mean_frequency_rad = twoPi_() * report_.mean_frequency_hz;

        report_.t1_s   = (m1 > EPS_) ? (m0 / m1) : 0.0;             // Tm01
        report_.tm02_s = (m2 > EPS_) ? std::sqrt(m0 / m2) : 0.0;    // Tm02 / Tz
        report_.te_s   = (negative_moments_ && m0 > EPS_) ? (mn1 / m0) : 0.0;
        report_.tm0m1_s = report_.te_s;
        report_.tm1m1_s = (negative_moments_ && mn1 > EPS_) ? (m0 / mn1) : 0.0;
        report_.mean_group_period_s = report_.te_s;

        report_.upcrossing_rate_hz   = (m0 > EPS_) ? std::sqrt(m2 / m0) : 0.0;
        report_.downcrossing_rate_hz = report_.upcrossing_rate_hz;

        report_.wave_count = estimateWaveCountWithCI_(observation_duration_s_, 0.95);

        report_.rms_displacement_m = std::sqrt(std::max(m0, 0.0));
        report_.hs_m = 4.0 * report_.rms_displacement_m;

        report_.deep_water_wavelength_tp_m =
            (report_.tp_s > EPS_) ? (G_ * report_.tp_s * report_.tp_s / (2.0 * M_PI)) : 0.0;

        report_.wave_steepness =
            (report_.deep_water_wavelength_tp_m > EPS_)
                ? (report_.hs_m / report_.deep_water_wavelength_tp_m)
                : 0.0;

        report_.nonlinearity_parameter =
            (report_.deep_water_wavelength_tp_m > EPS_)
                ? 0.5 * report_.hs_m * (twoPi_() / report_.deep_water_wavelength_tp_m)
                : 0.0;

        if (extended_metrics_ && m0 > EPS_) {
            const double mu = m1 / m0;
            const double raw2 = m2 / m0;
            const double raw3 = m3 / m0;
            const double raw4 = m4 / m0;

            report_.mu2 = raw2 - mu * mu;
            report_.mu3 = raw3 - 3.0 * mu * raw2 + 2.0 * mu * mu * mu;
            report_.mu4 = raw4 - 4.0 * mu * raw3 + 6.0 * mu * mu * raw2 - 3.0 * mu * mu * mu * mu;
        }

        report_.rbw = (m1 > EPS_)
            ? safeSqrt_(std::max((m0 * m2) / (m1 * m1) - 1.0, 0.0))
            : 0.0;

        report_.regularity_spec = std::clamp(std::exp(-BETA_SPEC_ * report_.rbw), 0.0, 1.0);
        report_.narrowness_nu = report_.rbw;
    }

    void computeHeightsAndExtremes_() {
        const double sigma = report_.rms_displacement_m;
        if (!(sigma > EPS_)) return;

        report_.h1_10_crest_m  = sigma * std::sqrt(-2.0 * std::log(0.10));
        report_.h1_100_crest_m = sigma * std::sqrt(-2.0 * std::log(0.01));

        const double N = std::max(report_.wave_count.expected, 1.0);
        report_.most_probable_max_crest_m = sigma * std::sqrt(2.0 * std::log(N));

        const double aN = std::sqrt(2.0 * std::log(N));
        const double gamma_euler = 0.5772156649;
        report_.expected_max_crest_m = (aN > EPS_) ? sigma * (aN + gamma_euler / aN) : 0.0;

        report_.crest_exceed_prob_hs =
            std::exp(-(report_.hs_m * report_.hs_m) / (2.0 * sigma * sigma));

        const double Lambda = 0.5 * report_.wave_steepness;
        report_.crest_exceed_prob_tayfun_hs =
            report_.crest_exceed_prob_hs *
            std::exp(Lambda * std::pow(report_.hs_m / std::max(sigma, EPS_), 3.0));

        report_.pot_mean_excess_m = computePOTMeanExcess_(report_.hs_m);
        report_.weibull_return_height_m = computeWeibullReturnHeight_(100.0, observation_duration_s_);
    }

    void computeBandwidthAndShape_() {
        const double m0 = report_.m0;
        const double m1 = report_.m1;
        const double m2 = report_.m2;
        const double m4 = report_.m4;

        report_.bandwidth_clh =
            (m0 > EPS_ && m2 > EPS_)
                ? safeSqrt_(std::max(1.0 - (m1 * m1) / (m0 * m2), 0.0))
                : 0.0;

        report_.bandwidth_goda =
            (m1 > EPS_) ? safeSqrt_(std::max((m0 * m2) / (m1 * m1) - 1.0, 0.0)) : 0.0;

        report_.bandwidth_kuik =
            (m1 > EPS_)
                ? safeSqrt_(std::max((m0 * m2) - (m1 * m1), 0.0)) / m1
                : 0.0;

        report_.longuet_higgins_width = report_.bandwidth_goda;

        report_.spectral_bandwidth_hz =
            (m0 > EPS_)
                ? safeSqrt_(std::max((m2 / m0) - (m1 / m0) * (m1 / m0), 0.0))
                : 0.0;

        report_.spectral_narrowness_ratio =
            (m0 > EPS_ && m2 > EPS_) ? (m1 * m1) / (m0 * m2) : 0.0;

        if (extended_metrics_) {
            if (report_.mu2 > EPS_) {
                report_.spectral_skewness =
                    report_.mu3 / std::pow(report_.mu2, 1.5);
                report_.spectral_kurtosis =
                    report_.mu4 / (report_.mu2 * report_.mu2);
                report_.spectral_excess_kurtosis =
                    report_.spectral_kurtosis - 3.0;
            }

            if (m2 > EPS_) {
                report_.peakedness_ochi_q = (m0 * m4) / (m2 * m2);
                report_.benassai_parameter = report_.peakedness_ochi_q;
                report_.peakedness_goda = 2.0 * (report_.peakedness_ochi_q - 1.0);
                report_.peak_enhancement_gamma = std::clamp(
                    0.25 + 0.75 * std::sqrt(std::max(report_.peakedness_ochi_q, 0.0)),
                    1.0, 10.0);
            }

            report_.bandwidth_epsilon =
                (m0 > EPS_ && m4 > EPS_)
                    ? safeSqrt_(std::max(1.0 - (m2 * m2) / (m0 * m4), 0.0))
                    : 0.0;

            const double ratio = (report_.tm02_s > EPS_) ? (report_.tp_s / report_.tm02_s) : 1.0;
            if (report_.peak_enhancement_gamma < 1.5 && ratio > 1.2) {
                report_.spectrum_type = SpectrumType::Swell;
            } else if (report_.peak_enhancement_gamma > 3.0 && ratio < 0.9) {
                report_.spectrum_type = SpectrumType::WindSea;
            } else {
                report_.spectrum_type = SpectrumType::Mixed;
            }
        } else {
            report_.spectrum_type = SpectrumType::Unknown;
        }

        report_.groupiness_factor = 1.0 / std::max(1.0 + report_.rbw, 1.0);
        report_.mean_group_length = (report_.rbw > EPS_) ? (1.0 / report_.rbw) : 0.0;
        report_.group_height_factor = 1.0 / std::max(1.0 + 0.5 * report_.rbw, 1.0);

        report_.bfi =
            (report_.rbw > EPS_) ? (report_.wave_steepness / report_.rbw) : 0.0;

        const double p_exceed_hs = std::max(report_.crest_exceed_prob_hs, EPS_);
        report_.expected_run_length_above_hs = 1.0 / p_exceed_hs;
    }

    void computeEnergetics_() {
        report_.energy_flux_period_s = report_.te_s;
        report_.wave_energy_density_j_m2 = RHO_ * G_ * report_.m0;
        report_.wave_power_w_m =
            (RHO_ * G_ * G_ / (64.0 * M_PI)) *
            report_.hs_m * report_.hs_m * report_.te_s;

        // Deep-water normal-incidence scalar proxy
        report_.radiation_stress_xx_n_m = 0.5 * report_.wave_energy_density_j_m2;
        report_.radiation_stress_yy_n_m = 0.5 * report_.wave_energy_density_j_m2;
        report_.radiation_stress_xy_n_m = 0.0;

        if (std::isfinite(wec_extracted_w_m_) && report_.wave_power_w_m > EPS_) {
            report_.wec_capture_width_ratio = wec_extracted_w_m_ / report_.wave_power_w_m;
        }
    }

    void computeNonlinearAndBreaking_() {
        if (std::isfinite(depth_m_) && depth_m_ > EPS_) {
            report_.ursell_number =
                (report_.deep_water_wavelength_tp_m > EPS_)
                    ? (report_.hs_m *
                       report_.deep_water_wavelength_tp_m *
                       report_.deep_water_wavelength_tp_m) /
                      (depth_m_ * depth_m_ * depth_m_)
                    : 0.0;

            report_.depth_limited_breaking_height_m = 0.78 * depth_m_;

            if (report_.hs_m > EPS_) {
                const double ratio = 0.78 * depth_m_ / report_.hs_m;
                report_.breaking_probability =
                    std::exp(-0.5 * std::clamp(ratio, 0.0, 5.0) * std::clamp(ratio, 0.0, 5.0));
            }

            if (report_.tp_s > EPS_) {
                const double omega = twoPi_() / report_.tp_s;
                const double k = solveWavenumber_(omega, depth_m_);
                const double denom = std::sinh(std::max(k * depth_m_, EPS_));
                report_.bottom_orbital_velocity_mps =
                    (M_PI * report_.hs_m / std::max(report_.tp_s, EPS_)) / std::max(denom, EPS_);
            }
        } else {
            report_.breaking_probability = std::clamp(10.0 * report_.wave_steepness, 0.0, 1.0);
            report_.depth_limited_breaking_height_m = 0.0;
        }

        if (std::isfinite(wind_speed_mps_) && wind_speed_mps_ > EPS_) {
            const double cp =
                (report_.tp_s > EPS_) ? (G_ * report_.tp_s / (2.0 * M_PI)) : 0.0;
            report_.wave_age = cp / wind_speed_mps_;

            if (report_.wave_age < 0.8) report_.wave_age_class = 0.0;
            else if (report_.wave_age > 1.2) report_.wave_age_class = 2.0;
            else report_.wave_age_class = 1.0;

            report_.sea_swell_partition = (report_.wave_age > 1.2) ? 1.0 : 0.5;
        }
    }

    void computeComfortMetrics_() {
        const double accel_rms =
            (report_.m2 > EPS_) ? std::sqrt(report_.m2) : 0.0;
        report_.vertical_motion_intensity = accel_rms;

        if (!(accel_rms > EPS_) || !(report_.mean_frequency_hz > EPS_)) {
            report_.msdv = 0.0;
            report_.seasickness_incidence_pct = 0.0;
            report_.motion_comfort_level_0_100 = 100.0;
            report_.motion_character = 0.0;
            report_.time_to_onset_min = std::numeric_limits<double>::infinity();
            return;
        }

        report_.msdv = accel_rms * std::sqrt(std::max(exposure_hours_, 0.0));

        const double f_dom = report_.mean_frequency_hz;
        const double a = 0.5, b = 2.0, c = 0.4, k = 0.7;
        const double msdv2 = accel_rms * std::pow(f_dom, b / a) * std::pow(exposure_hours_, c);

        report_.seasickness_incidence_pct = std::clamp(
            100.0 * (1.0 - std::exp(-k * msdv2 * susceptibility_)),
            0.0, 100.0);

        if (report_.msdv <= 0.1) report_.motion_comfort_level_0_100 = 100.0;
        else if (report_.msdv <= 0.2) report_.motion_comfort_level_0_100 = 80.0;
        else if (report_.msdv <= 0.3) report_.motion_comfort_level_0_100 = 60.0;
        else if (report_.msdv <= 0.4) report_.motion_comfort_level_0_100 = 40.0;
        else if (report_.msdv <= 0.5) report_.motion_comfort_level_0_100 = 20.0;
        else report_.motion_comfort_level_0_100 = 0.0;

        const double optimal_freq = 0.2;
        report_.motion_character = std::clamp(
            std::exp(-std::pow((f_dom - optimal_freq) / 0.1, 2.0)),
            0.0, 1.0);

        report_.time_to_onset_min = std::max(
            1.0,
            30.0 / (accel_rms * std::pow(f_dom, 0.7) * std::max(susceptibility_, EPS_)));
    }

    void computeQualityMetrics_(const std::vector<double>& freqs_hz,
                                const std::vector<double>& S) {
        // Crude spectral SNR estimate: compare in-band median to tail median.
        const double fp = report_.fp_hz;
        if (!(fp > EPS_)) {
            report_.snr_db = std::numeric_limits<double>::quiet_NaN();
            report_.data_quality_0_1 = 0.0;
            return;
        }

        std::vector<double> band_vals;
        std::vector<double> tail_vals;

        for (std::size_t i = 0; i < freqs_hz.size(); ++i) {
            const double f = freqs_hz[i];
            const double s = std::max(0.0, S[i]);

            if (f >= 0.7 * fp && f <= 1.4 * fp) band_vals.push_back(s);
            if (f >= std::max(1.8 * fp, report_.lowfreq_guard_hz) || f <= 0.6 * report_.lowfreq_guard_hz) {
                tail_vals.push_back(s);
            }
        }

        const double sig = median_(band_vals);
        const double noi = std::max(median_(tail_vals), EPS_);
        report_.snr_db = (sig > 0.0) ? 10.0 * std::log10(sig / noi) : std::numeric_limits<double>::quiet_NaN();

        const double snr_score = std::isfinite(report_.snr_db)
            ? std::clamp((report_.snr_db - 0.0) / 20.0, 0.0, 1.0)
            : 0.0;

        const double reg_score = std::clamp(report_.regularity_spec, 0.0, 1.0);
        const double bw_score  = std::exp(-0.5 * std::max(0.0, report_.rbw));

        double stab_hs = std::isfinite(report_.temporal_stability_hs) ? report_.temporal_stability_hs : 0.75;
        double stab_tp = std::isfinite(report_.temporal_stability_tp) ? report_.temporal_stability_tp : 0.75;

        report_.data_quality_0_1 = std::clamp(
            0.30 * snr_score +
            0.25 * reg_score +
            0.20 * bw_score +
            0.125 * stab_hs +
            0.125 * stab_tp,
            0.0, 1.0);
    }

    void updateStability_(double dt_s) {
        if (!(dt_s > 0.0) || !(report_.valid || report_.m0 > EPS_)) return;

        const double alpha = 1.0 - std::exp(-dt_s / 180.0);

        if (!has_history_) {
            hs_ema_ = report_.hs_m;
            tp_ema_ = report_.tp_s;
            hs_var_ = 0.0;
            tp_var_ = 0.0;
            has_history_ = true;
        } else {
            const double dh = report_.hs_m - hs_ema_;
            hs_ema_ += alpha * dh;
            hs_var_ = (1.0 - alpha) * hs_var_ + alpha * dh * dh;

            const double dtp = report_.tp_s - tp_ema_;
            tp_ema_ += alpha * dtp;
            tp_var_ = (1.0 - alpha) * tp_var_ + alpha * dtp * dtp;
        }

        report_.temporal_stability_hs =
            std::clamp(1.0 - std::sqrt(hs_var_) / std::max(hs_ema_, EPS_), 0.0, 1.0);

        report_.temporal_stability_tp =
            std::clamp(1.0 - std::sqrt(tp_var_) / std::max(tp_ema_, EPS_), 0.0, 1.0);
    }

    WaveCountEstimate estimateWaveCountWithCI_(double duration_s, double confidence) const {
        WaveCountEstimate out;
        out.confidence = confidence;

        if (!(duration_s > 0.0)) return out;

        const double nu = report_.upcrossing_rate_hz;
        const double Nexp = std::max(0.0, nu * duration_s);
        const int N = static_cast<int>(std::round(Nexp));

        const double alpha = std::clamp(1.0 - confidence, 0.0, 1.0);
        double lower = 0.0, upper = 0.0;

        if (N == 0) {
            lower = 0.0;
            upper = 0.5 * chi2Quantile_(1.0 - alpha / 2.0, 2 * (N + 1));
        } else {
            lower = 0.5 * chi2Quantile_(alpha / 2.0, 2 * N);
            upper = 0.5 * chi2Quantile_(1.0 - alpha / 2.0, 2 * (N + 1));
        }

        out.expected = Nexp;
        out.ci_lower = lower;
        out.ci_upper = upper;
        return out;
    }

    double computeWeibullReturnHeight_(double T_return, double duration_s) const {
        if (!(report_.hs_m > EPS_) || !(duration_s > 0.0)) return 0.0;
        const double k = 2.0;
        const double lambda = report_.hs_m / std::sqrt(std::log(2.0));
        const double p = 1.0 - 1.0 / std::max(T_return, 1.0);
        return lambda * std::pow(-std::log(1.0 - p), 1.0 / k);
    }

    double computePOTMeanExcess_(double threshold) const {
        if (!(report_.hs_m > EPS_)) return 0.0;
        const double sigma = report_.hs_m / std::sqrt(2.0);
        const double thr = (threshold > EPS_) ? threshold : report_.hs_m;
        const double p = std::exp(-thr * thr / (2.0 * sigma * sigma));
        return (p > EPS_) ? sigma * p / (1.0 - p) : 0.0;
    }

    static double solveWavenumber_(double omega, double depth_m) {
        if (!(depth_m > EPS_)) return (omega * omega / G_);

        double k = omega * omega / G_;
        for (int i = 0; i < 5; ++i) {
            const double kh = k * depth_m;
            const double t = std::tanh(kh);
            const double sech2 = 1.0 / (std::cosh(kh) * std::cosh(kh));
            const double f = G_ * k * t - omega * omega;
            const double df = G_ * (t + k * depth_m * sech2);
            k -= f / std::max(df, EPS_);
        }
        return std::max(k, EPS_);
    }

    static double median_(std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        const std::size_t n = v.size();
        if (n & 1U) return v[n / 2];
        return 0.5 * (v[n / 2 - 1] + v[n / 2]);
    }

    static double erfinvApprox_(double x) {
        x = std::clamp(x, -0.999999, 0.999999);
        const double a = 0.147;
        const double ln = std::log(1.0 - x * x);
        const double tt1 = 2.0 / (M_PI * a) + 0.5 * ln;
        const double tt2 = (1.0 / a) * ln;
        const double sign = (x < 0.0) ? -1.0 : 1.0;
        double y = sign * std::sqrt(std::sqrt(tt1 * tt1 - tt2) - tt1);

        for (int i = 0; i < 2; ++i) {
            const double ey  = std::erf(y) - x;
            const double dy  = (2.0 / std::sqrt(M_PI)) * std::exp(-y * y);
            const double d2y = -2.0 * y * dy;
            const double denom = 2.0 * dy * dy - ey * d2y;
            if (std::fabs(denom) < 1e-12) break;
            y -= (2.0 * ey * dy) / denom;
        }
        return y;
    }

    static double normalQuantile_(double p) {
        p = std::clamp(p, 1e-7, 1.0 - 1e-7);
        return std::sqrt(2.0) * erfinvApprox_(2.0 * p - 1.0);
    }

    static double chi2Quantile_(double p, int k) {
        if (p <= 0.0) return 0.0;
        if (p >= 1.0) return std::numeric_limits<double>::infinity();
        const double z = normalQuantile_(p);
        const double v = static_cast<double>(k);
        const double h = 2.0 / (9.0 * v);
        return v * std::pow(1.0 - h + z * std::sqrt(h), 3.0);
    }
};

/*

// Usage:

SeaMetricsFromSpectrum sea(true, true);
sea.setDepthMeters(50.0);
sea.setWindSpeedMps(12.0);
sea.setObservationDuration(600.0);

sea.updateFromSpectrum(freqs_hz, S_eta, lowfreq_cut_hz, dt_s);
const auto& r = sea.report();

*/
