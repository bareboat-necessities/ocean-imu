#pragma once

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <array>
#include <cmath>
#include <algorithm>

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

        last_lowfreq_cut_hz_ =
            WaveSpectrumShared::estimate_lowfreq_cut_from_accel<Nfreq>(
                S_aa_true_arr, freqs_, Tblk, hp_f0_hz);

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
            const double den = w * w + lam * lam;

            double S_eta = (den > 0.0) ? (S_aa_true_arr[i] / (den * den)) : 0.0;
            if (!std::isfinite(S_eta) || S_eta < 0.0) S_eta = 0.0;

            S_eta *= WaveSpectrumShared::lowfreq_taper(f, last_lowfreq_cut_hz_);

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

        WaveSpectrumShared::finalize_displacement_spectrum_goertzel_inplace<Nfreq>(
            lastSpectrum_, S_aa_true_arr, freqs_, last_lowfreq_cut_hz_); 
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
    double ema_alpha_low = 0.20;
    double ema_alpha_high = 0.06;
    bool have_ema = false;
    Eigen::Matrix<double, Nfreq, 1> psd_ema_;

    int writeIndex = 0;
    int decimCounter = 0;
    int filledSamples = 0;
    bool isWarm = false;

    double last_lowfreq_cut_hz_ = 0.0;
};
