#pragma once

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <array>
#include <cmath>
#include <limits>
#include <algorithm>

#include "spectrum/SpectrumStats.h"
#include "spectrum/WaveSpectrumShared.h"

/*
    Classic Goertzel spectrum estimator.

    Hybrid rollback:
      - keeps decimation sanitization
      - keeps decimation-aware analysis_fmax_hz_
      - keeps guarded peak search from WaveSpectrumShared
      - restores older milder cutoff learning / smoothing / cleanup
      - avoids aggressive post-inversion low-frequency suppression
*/

template<int Nfreq = 32, int Nblock = 256>
class EIGEN_ALIGN_MAX WaveSpectrumEstimator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static constexpr double g = 9.80665;
    using Vec = Eigen::Matrix<double, Nfreq, 1>;
    using Biquad = WaveSpectrumShared::Biquad;

    struct PMFitResult {
        double alpha;
        double fp;
        double cost;
    };

    WaveSpectrumEstimator(double fs_raw_ = 200.0,
                          int decimFactor_ = 30,
                          bool hannEnabled_ = true)
        : fs_raw(fs_raw_),
          decimFactor(std::max(1, decimFactor_)),
          hannEnabled(hannEnabled_) {
        fs = fs_raw / double(decimFactor);

        const double fny_dec = fs_raw / (2.0 * double(decimFactor));
        const double lp_cutoffHz = 0.32 * fny_dec;

        // Less conservative than the recent regression-prone tuning.
        analysis_fmax_hz_ = WaveSpectrumShared::compute_safe_analysis_fmax(
            fs_raw, decimFactor, 1.2, 0.32, 0.90, 0.90, 0.04);

        WaveSpectrumShared::build_log_frequency_grid<Nfreq>(
            freqs_, f_edges_, df_, 0.04, analysis_fmax_hz_);

        for (int i = 0; i < Nfreq; ++i) {
            const double omega_rs = 2.0 * M_PI * freqs_[i] / fs;
            coeffs_[i] = 2.0 * std::cos(omega_rs);
        }

        double sumsq = 0.0;
        for (int n = 0; n < Nblock; ++n) {
            window_[n] = hannEnabled
                ? 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (Nblock - 1)))
                : 1.0;
            sumsq += window_[n] * window_[n];
        }
        window_sum_sq = sumsq;

        reset();

        WaveSpectrumShared::designLowpassBiquad(lp1_, lp_cutoffHz, fs_raw);
        WaveSpectrumShared::designLowpassBiquad(lp2_, lp_cutoffHz, fs_raw);
        WaveSpectrumShared::designHighpassBiquad(hp1_, hp_f0_hz, fs_raw);
        WaveSpectrumShared::designHighpassBiquad(hp2_, hp_f0_hz, fs_raw);
    }

    void reset() {
        buffer_.fill(0.0);
        writeIndex = 0;
        decimCounter = 0;
        filledSamples = 0;
        isWarm = false;
        lastSpectrum_.setZero();

        hp1_.reset();
        hp2_.reset();
        lp1_.reset();
        lp2_.reset();

        have_ema = false;
        psd_ema_.setZero();
        last_lowfreq_cut_hz_ = 0.0;
    }

    bool processSample(double x_raw) {
        const double x_hp = hp2_.process(hp1_.process(x_raw));
        const double y = lp2_.process(lp1_.process(x_hp));

        if (++decimCounter < decimFactor) return false;
        decimCounter = 0;

        buffer_[writeIndex] = y;
        writeIndex = (writeIndex + 1) % Nblock;
        filledSamples++;

        if (filledSamples >= Nblock) isWarm = true;

        if ((filledSamples % Nblock) == 0) {
            computeSpectrum();
            return true;
        }
        return false;
    }

    Vec getDisplacementSpectrum() const { return lastSpectrum_; }
    std::array<double, Nfreq> getFrequencies() const { return freqs_; }
    bool ready() const { return isWarm; }
    double getLowfreqCutHz() const { return last_lowfreq_cut_hz_; }
    double getAnalysisFmaxHz() const { return analysis_fmax_hz_; }

    double computeHs() const {
        return SpectrumStats::compute_hs<Nfreq>(lastSpectrum_, df_);
    }

    double estimateFp() const {
        return WaveSpectrumShared::estimate_fp_with_guard<Nfreq>(
            lastSpectrum_, freqs_, last_lowfreq_cut_hz_);
    }

    double estimateTp() const {
        const double fp = estimateFp();
        return (fp > 0.0) ? (1.0 / fp) : 0.0;
    }

    PMFitResult fitPiersonMoskowitz() const {
        Vec S_obs = lastSpectrum_;
        for (int i = 0; i < Nfreq; ++i) {
            if (S_obs[i] <= 0.0) S_obs[i] = 1e-12;
        }

        auto cost_fn = [&](double a, double fp) {
            const double omega_p = 2.0 * M_PI * fp;
            double cost = 0.0;
            constexpr double beta = 0.74;

            for (int i = 0; i < Nfreq - 1; ++i) {
                const double df = df_[i];
                const double f = freqs_[i];
                double model = a * g * g * std::pow(2.0 * M_PI * f, -5.0)
                             * std::exp(-beta * std::pow(omega_p / (2.0 * M_PI * f), 4.0));
                if (model <= 0.0) model = 1e-12;
                const double d = safe_log_(S_obs[i]) - safe_log_(model);
                cost += df * d * d;
            }
            return cost;
        };

        constexpr int N_fp_search = 32;
        constexpr double fp_min = 0.05;
        constexpr double fp_transition = 0.1;
        constexpr double fp_max = 1.0;

        std::array<double, N_fp_search> fp_grid{};
        const int n_log = static_cast<int>(N_fp_search * 0.4);
        const int n_lin = N_fp_search - n_log;

        for (int i = 0; i < n_log; ++i) {
            const double t = double(i) / double(std::max(1, n_log - 1));
            fp_grid[i] = fp_min * std::pow(fp_transition / fp_min, t);
        }
        for (int i = 0; i < n_lin; ++i) {
            const double t = double(i) / double(std::max(1, n_lin - 1));
            fp_grid[n_log + i] = fp_transition + t * (fp_max - fp_transition);
        }

        double bestA = 1e-5;
        double bestFp = fp_grid[0];
        double bestC = std::numeric_limits<double>::infinity();

        for (int ia = 0; ia < 8; ++ia) {
            const double a = 1e-5 + ia * (1.0 - 1e-5) / 7.0;
            for (int ifp = 0; ifp < N_fp_search; ++ifp) {
                const double fp = fp_grid[ifp];
                const double c = cost_fn(a, fp);
                if (c < bestC) {
                    bestC = c;
                    bestA = a;
                    bestFp = fp;
                }
            }
        }

        double alpha = bestA;
        double fp = bestFp;
        double stepA = 0.1;
        double stepFp = 0.1;

        for (int iter = 0; iter < 40; ++iter) {
            bool improved = false;
            double c = 0.0;

            c = cost_fn(alpha + stepA, fp);
            if (c < bestC) { bestC = c; alpha += stepA; improved = true; }

            c = cost_fn(alpha - stepA, fp);
            if (c < bestC) { bestC = c; alpha -= stepA; improved = true; }

            c = cost_fn(alpha, fp + stepFp);
            if (c < bestC) { bestC = c; fp += stepFp; improved = true; }

            c = cost_fn(alpha, fp - stepFp);
            if (c < bestC) { bestC = c; fp -= stepFp; improved = true; }

            if (!improved) {
                stepA *= 0.5;
                stepFp *= 0.5;
                if (stepA < 1e-12 && stepFp < 1e-12) break;
            }
        }

        return {alpha, fp, bestC};
    }

    void set_regularization_f0(double f0_hz) {
        reg_f0_hz = std::max(1e-6, f0_hz);
    }

    void set_highpass_f0(double f0_hz) {
        hp_f0_hz = std::max(0.0, f0_hz);
        WaveSpectrumShared::designHighpassBiquad(hp1_, hp_f0_hz, fs_raw);
        WaveSpectrumShared::designHighpassBiquad(hp2_, hp_f0_hz, fs_raw);
    }

private:
    static double safe_log_(double x) {
        return std::log(std::max(x, 1e-12));
    }

    inline double alpha_for_f(double f) const {
        return WaveSpectrumShared::ema_alpha_for_f(
            f, freqs_[0], freqs_[Nfreq - 1], ema_alpha_low, ema_alpha_high);
    }

    double estimate_lowfreq_cut_from_accel_(const std::array<double, Nfreq>& S_aa_true_arr,
                                            double Tblk) const {
        if (Nfreq < 4) return 0.0;

        std::array<double, Nfreq> E{};
        std::array<double, Nfreq> Es{};

        for (int i = 0; i < Nfreq; ++i) {
            E[i] = std::max(0.0, S_aa_true_arr[i]) * std::max(freqs_[i], 1e-12);
        }
        WaveSpectrumShared::smooth_3tap_array<Nfreq>(E, Es);

        const double f_floor_hz = std::max(
            std::max(1.04 * freqs_[0], (Tblk > 1e-12) ? (1.20 / Tblk) : 0.0),
            1.05 * hp_f0_hz);

        int i_floor = 0;
        while (i_floor + 1 < Nfreq && freqs_[i_floor + 1] <= f_floor_hz) ++i_floor;

        int i_peak = std::max(1, i_floor);
        double e_peak = Es[i_peak];
        for (int i = i_peak + 1; i < Nfreq - 1; ++i) {
            if (Es[i] > e_peak) {
                e_peak = Es[i];
                i_peak = i;
            }
        }

        if (!(e_peak > 0.0) || i_peak <= i_floor) {
            return f_floor_hz;
        }

        int i_valley = i_floor;
        double e_valley = Es[i_floor];
        for (int i = i_floor + 1; i < i_peak; ++i) {
            if (Es[i] < e_valley) {
                e_valley = Es[i];
                i_valley = i;
            }
        }

        const bool left_edge_raised =
            (E[0] > 1.15 * std::max(e_valley, 1e-18)) ||
            (Nfreq > 2 && E[1] > 1.06 * std::max(E[2], 1e-18));

        const bool valley_is_good =
            (i_valley > i_floor) &&
            (e_valley < 0.68 * e_peak) &&
            (freqs_[i_valley] < 0.72 * freqs_[i_peak]);

        if (left_edge_raised && valley_is_good) {
            return std::max(f_floor_hz, freqs_[i_valley]);
        }

        const double f_rel = 0.60 * freqs_[i_peak];
        const double f_cap = 0.90 * freqs_[i_peak];
        return std::max(f_floor_hz, std::min(f_rel, f_cap));
    }

    static double lowfreq_taper_(double f_hz, double f_cut_hz) {
        if (!(f_cut_hz > 0.0)) return 1.0;

        const double f_lo = 0.60 * f_cut_hz;
        const double f_hi = 1.35 * f_cut_hz;

        if (f_hz >= f_hi) return 1.0;
        if (f_hz <= f_lo) return 0.0;

        const double t = std::clamp((f_hz - f_lo) / (f_hi - f_lo), 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
    }

    void suppress_lowfreq_from_cut_inplace_mild_() {
        if (Nfreq < 4) return;
        if (!(last_lowfreq_cut_hz_ > 0.0)) return;

        const double f_eff = 1.00 * last_lowfreq_cut_hz_;

        std::array<double, Nfreq> E{};
        for (int i = 0; i < Nfreq; ++i) {
            E[i] = std::max(0.0, double(lastSpectrum_[i])) * std::max(freqs_[i], 1e-12);
        }

        int i_cut = 0;
        while (i_cut + 1 < Nfreq && freqs_[i_cut + 1] <= f_eff) ++i_cut;
        if (i_cut <= 0) return;

        const double f_ref = std::max(freqs_[i_cut], 1e-12);
        const double E_ref = std::max(E[i_cut], 1e-18);

        double prev = E_ref;
        const double shape_pow = 2.0;
        const double slope_cap = 1.08;

        for (int i = i_cut - 1; i >= 0; --i) {
            const double r = std::max(freqs_[i] / f_ref, 0.0);
            const double cap = E_ref * std::pow(r, shape_pow);

            double v = std::max(0.0, E[i]);
            v = std::min(v, cap);
            v = std::min(v, slope_cap * prev);

            E[i] = v;
            prev = v;
        }

        for (int i = 0; i < Nfreq; ++i) {
            lastSpectrum_[i] = E[i] / std::max(freqs_[i], 1e-12);
        }
    }

    void computeSpectrum() {
        const int blockSize = std::min(filledSamples, Nblock);
        const int startIdx = (writeIndex + Nblock - blockSize) % Nblock;
        const int N = blockSize;

        double sumx = 0.0, sumn = 0.0, sumn2 = 0.0, sumnx = 0.0;
        {
            int idx = startIdx;
            for (int n = 0; n < N; ++n) {
                const double xn = buffer_[idx];
                sumx += xn;
                sumn += n;
                sumn2 += double(n) * double(n);
                sumnx += double(n) * xn;
                idx = (idx + 1) % Nblock;
            }
        }

        const double denom = double(N) * sumn2 - sumn * sumn;
        const double a_lin = (denom != 0.0)
            ? (double(N) * sumnx - sumn * sumx) / denom
            : 0.0;
        const double b_lin = (sumx - a_lin * sumn) / double(N);

        const double base_scale =
            (window_sum_sq > 0.0 && fs > 0.0) ? (2.0 / (fs * window_sum_sq)) : 0.0;

        const double Tblk = (fs > 0.0) ? (double(N) / fs) : 0.0;
        const double f_blk = (Tblk > 0.0) ? (1.0 / (6.0 * Tblk)) : 0.0;

        std::array<double, Nfreq> S_aa_true_arr{};

        for (int i = 0; i < Nfreq; ++i) {
            const double f = freqs_[i];

            double s1 = 0.0, s2 = 0.0;
            int idx = startIdx;
            for (int n = 0; n < N; ++n) {
                const double xn = buffer_[idx];
                const double xdet = xn - (a_lin * n + b_lin);
                const double xw = xdet * window_[n];
                const double s_new = xw + coeffs_[i] * s1 - s2;
                s2 = s1;
                s1 = s_new;
                idx = (idx + 1) % Nblock;
            }

            const double power = s1 * s1 + s2 * s2 - s1 * s2 * coeffs_[i];
            const double S_aa_meas = std::max(0.0, power * base_scale);

            const double omega_raw = 2.0 * M_PI * f / fs_raw;
            const double H2_hp =
                WaveSpectrumShared::biquad_mag2_raw(hp1_, omega_raw) *
                WaveSpectrumShared::biquad_mag2_raw(hp2_, omega_raw);

            constexpr double hp_deconv_reg = 0.10;
            const double inv_hp =
                H2_hp / (H2_hp * H2_hp + hp_deconv_reg * hp_deconv_reg);

            S_aa_true_arr[i] = std::max(0.0, S_aa_meas * inv_hp);
        }

        last_lowfreq_cut_hz_ = estimate_lowfreq_cut_from_accel_(S_aa_true_arr, Tblk);

        const double f_knee = std::max({
            reg_f0_hz,
            0.60 * last_lowfreq_cut_hz_,
            0.95 * hp_f0_hz,
            0.90 * f_blk
        });

        const double lam = 2.0 * M_PI * f_knee;

        for (int i = 0; i < Nfreq; ++i) {
            const double f = freqs_[i];
            const double w = 2.0 * M_PI * f;
            const double den = (w * w + lam * lam);

            double S_eta = (den > 0.0) ? (S_aa_true_arr[i] / (den * den)) : 0.0;
            if (!std::isfinite(S_eta) || S_eta < 0.0) S_eta = 0.0;

            S_eta *= lowfreq_taper_(f, last_lowfreq_cut_hz_);

            if (use_psd_ema) {
                const double a = alpha_for_f(f);
                if (!have_ema) psd_ema_[i] = S_eta;
                else           psd_ema_[i] = (1.0 - a) * psd_ema_[i] + a * S_eta;
                lastSpectrum_[i] = psd_ema_[i];
            } else {
                lastSpectrum_[i] = S_eta;
            }
        }

        have_ema = true;

        WaveSpectrumShared::smooth_logfreq_3tap_inplace<Nfreq>(lastSpectrum_, freqs_);
        suppress_lowfreq_from_cut_inplace_mild_();
    }

    double fs_raw = 0.0;
    double fs = 0.0;
    int decimFactor = 1;
    bool hannEnabled = true;

    double reg_f0_hz = 0.008;
    double hp_f0_hz = 0.025;
    double analysis_fmax_hz_ = 1.2;

    std::array<double, Nfreq> freqs_{};
    std::array<double, Nfreq> coeffs_{};
    std::array<double, Nfreq + 1> f_edges_{};
    std::array<double, Nfreq> df_{};

    std::array<double, Nblock> buffer_{};
    std::array<double, Nblock> window_{};
    double window_sum_sq = 1.0;

    Biquad hp1_, hp2_, lp1_, lp2_;

    Eigen::Matrix<double, Nfreq, 1> lastSpectrum_;

    bool use_psd_ema = true;
    double ema_alpha_low = 0.06;
    double ema_alpha_high = 0.14;
    bool have_ema = false;
    Eigen::Matrix<double, Nfreq, 1> psd_ema_;

    int writeIndex = 0;
    int decimCounter = 0;
    int filledSamples = 0;
    bool isWarm = false;

    double last_lowfreq_cut_hz_ = 0.0;
};
