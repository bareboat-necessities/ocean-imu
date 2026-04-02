#pragma once

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <array>
#include <cmath>
#include <algorithm>
#include <limits>

#include "spectrum/SpectrumStats.h"
#include "spectrum/WaveSpectrumShared.h"

template<int Nfreq = 32, int Nblock = 256>
class EIGEN_ALIGN_MAX WaveSpectrumEstimator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Vec = Eigen::Matrix<double, Nfreq, 1>;
    using Biquad = WaveSpectrumShared::Biquad;

    WaveSpectrumEstimator(double fs_raw_ = 200.0,
                          int decimFactor_ = 30,
                          bool hannEnabled_ = true)
        : fs_raw(fs_raw_), decimFactor(decimFactor_), hannEnabled(hannEnabled_) {
        fs = fs_raw / decimFactor;

        WaveSpectrumShared::build_log_frequency_grid<Nfreq>(freqs_, f_edges_, df_);

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

        const double fny_dec = fs_raw_ / (2.0 * decimFactor_);
        const double lp_cutoffHz = 0.32 * fny_dec;

        WaveSpectrumShared::designLowpassBiquad(lp1_, lp_cutoffHz, fs_raw_);
        WaveSpectrumShared::designLowpassBiquad(lp2_, lp_cutoffHz, fs_raw_);
        WaveSpectrumShared::designHighpassBiquad(hp1_, hp_f0_hz, fs_raw_);
        WaveSpectrumShared::designHighpassBiquad(hp2_, hp_f0_hz, fs_raw_);
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

    void set_regularization_f0(double f0_hz) {
        reg_f0_hz = std::max(1e-6, f0_hz);
    }

    void set_highpass_f0(double f0_hz) {
        hp_f0_hz = std::max(0.0, f0_hz);
        WaveSpectrumShared::designHighpassBiquad(hp1_, hp_f0_hz, fs_raw);
        WaveSpectrumShared::designHighpassBiquad(hp2_, hp_f0_hz, fs_raw);
    }

private:
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

    bool lowfreq_edge_is_raised_() const {
        if (Nfreq < 4) return false;
        if (!(last_lowfreq_cut_hz_ > 0.0)) return false;

        std::array<double, Nfreq> E{};
        for (int i = 0; i < Nfreq; ++i) {
            E[i] = std::max(0.0, double(lastSpectrum_[i])) * std::max(freqs_[i], 1e-12);
        }

        int i_cut = 0;
        const double f_eff = 1.00 * last_lowfreq_cut_hz_;
        while (i_cut + 1 < Nfreq && freqs_[i_cut + 1] <= f_eff) ++i_cut;
        if (i_cut < 2) return false;

        const double E_ref = std::max(E[i_cut], 1e-18);

        return
            (E[0] > 1.08 * E_ref) ||
            (E[1] > 1.04 * std::max(E[2], 1e-18));
    }

    bool suppress_unsupported_lowfreq_spikes_(const std::array<double, Nfreq>& S_aa_true_arr,
                                              int k_peak_acc) {
        if (Nfreq < 4) return false;
        if (!(last_lowfreq_cut_hz_ > 0.0)) return false;
        if (k_peak_acc <= 3 || k_peak_acc >= Nfreq) return false;

        double f_soft = 1.35 * last_lowfreq_cut_hz_;
        f_soft = std::min(f_soft, 0.72 * freqs_[k_peak_acc]);

        int k_soft = 0;
        while (k_soft + 1 < Nfreq && freqs_[k_soft + 1] <= f_soft) ++k_soft;

        const int k_stop = std::min(k_soft, std::max(1, k_peak_acc - 3));
        if (k_stop < 2 || k_stop + 1 >= Nfreq) return false;

        std::array<double, Nfreq> Eeta{};
        std::array<double, Nfreq> Eaa{};
        std::array<double, Nfreq> Eaa_s{};

        for (int i = 0; i < Nfreq; ++i) {
            const double f = std::max(freqs_[i], 1e-12);
            Eeta[i] = std::max(0.0, double(lastSpectrum_[i])) * f;
            Eaa[i]  = std::max(0.0, S_aa_true_arr[i]) * f;
        }
        WaveSpectrumShared::smooth_3tap_array<Nfreq>(Eaa, Eaa_s);

        const int i_ref = std::min(k_stop + 1, Nfreq - 1);
        const double Eref = std::max(Eeta[i_ref], 1e-12);
        const double fref = std::max(freqs_[i_ref], 1e-12);
        const double Eaa_ref = std::max(Eaa_s[i_ref], 1e-12);

        bool changed = false;

        for (int i = 1; i <= k_stop; ++i) {
            const bool local_peak =
                (Eeta[i] > Eeta[i - 1]) &&
                (Eeta[i] > Eeta[i + 1]);

            if (!local_peak) continue;

            const double r = std::max(freqs_[i], 1e-12) / fref;
            const double env_cap = Eref * std::pow(r, 1.55);
            const double accel_support = Eaa_s[i] / Eaa_ref;

            if (accel_support < 0.18 && Eeta[i] > 1.9 * env_cap) {
                Eeta[i] = std::max(env_cap, 0.60 * Eeta[i]);
                changed = true;
            }
        }

        if (changed) {
            for (int i = 0; i < Nfreq; ++i) {
                lastSpectrum_[i] = Eeta[i] / std::max(freqs_[i], 1e-12);
            }
        }

        return changed;
    }

    void suppress_lowfreq_from_cut_inplace_() {
        if (Nfreq < 4) return;
        if (!(last_lowfreq_cut_hz_ > 0.0)) return;

        const double f_eff = 1.05 * last_lowfreq_cut_hz_;

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
        const double shape_pow = 2.60;
        const double slope_cap = 1.05;

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

    void finalize_displacement_spectrum_(const std::array<double, Nfreq>& S_aa_true_arr) {
        WaveSpectrumShared::smooth_logfreq_3tap_custom_inplace<Nfreq>(
            lastSpectrum_, freqs_, 0.20, 0.62, 1.08, 0.92);

        const int k_peak_acc =
            WaveSpectrumShared::find_accel_peak_index<Nfreq>(
                S_aa_true_arr, freqs_, last_lowfreq_cut_hz_, 1.03);

        const bool changed =
            suppress_unsupported_lowfreq_spikes_(S_aa_true_arr, k_peak_acc);

        const bool edge_raised = lowfreq_edge_is_raised_();

        if (changed || edge_raised) {
            suppress_lowfreq_from_cut_inplace_();
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

        const double f_knee = std::max(
            std::max(reg_f0_hz, 0.60 * last_lowfreq_cut_hz_),
            std::max(0.95 * hp_f0_hz, 0.90 * f_blk));

        const double lam = 2.0 * M_PI * f_knee;

        for (int i = 0; i < Nfreq; ++i) {
            const double f = freqs_[i];
            const double w = 2.0 * M_PI * f;
            const double den = w * w + lam * lam;

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
        finalize_displacement_spectrum_(S_aa_true_arr);
    }

    double fs_raw = 0.0;
    double fs = 0.0;
    int decimFactor = 1;
    bool hannEnabled = true;

    double reg_f0_hz = 0.008;
    double hp_f0_hz = 0.025;

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
    double ema_alpha_low  = 0.06;
    double ema_alpha_high = 0.20;
    bool have_ema = false;
    Eigen::Matrix<double, Nfreq, 1> psd_ema_;

    int writeIndex = 0;
    int decimCounter = 0;
    int filledSamples = 0;
    bool isWarm = false;

    double last_lowfreq_cut_hz_ = 0.0;
};
