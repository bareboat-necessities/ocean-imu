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
class EIGEN_ALIGN_MAX WaveSpectrumEstimator_Wavelets {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Vec = Eigen::Matrix<double, Nfreq, 1>;
    using Biquad = WaveSpectrumShared::Biquad;

    WaveSpectrumEstimator_Wavelets(double fs_raw_ = 200.0,
                                   int decimFactor_ = 30,
                                   bool hannEnabled_ = true)
        : fs_raw(fs_raw_), decimFactor(decimFactor_), hannEnabled(hannEnabled_) {
        fs = fs_raw / decimFactor;

        WaveSpectrumShared::build_log_frequency_grid<Nfreq>(freqs_, f_edges_, df_);

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

        buildWaveletBank_();
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

    inline int tapIndex_(int i, int n) const {
        return i * Nblock + n;
    }

    void orthogonalize_wavelet_to_dc_and_ramp_(int i, int L, int half) {
        double sum_re = 0.0, sum_im = 0.0;
        double sumk_re = 0.0, sumk_im = 0.0;
        double sumk2 = 0.0;

        for (int n = 0; n < L; ++n) {
            const double k = double(n - half);
            const double re = double(wave_re_[tapIndex_(i, n)]);
            const double im = double(wave_im_[tapIndex_(i, n)]);

            sum_re += re;
            sum_im += im;
            sumk_re += k * re;
            sumk_im += k * im;
            sumk2 += k * k;
        }

        const double mean_re = sum_re / double(L);
        const double mean_im = sum_im / double(L);
        const double slope_re = (sumk2 > 0.0) ? (sumk_re / sumk2) : 0.0;
        const double slope_im = (sumk2 > 0.0) ? (sumk_im / sumk2) : 0.0;

        for (int n = 0; n < L; ++n) {
            const double k = double(n - half);
            wave_re_[tapIndex_(i, n)] =
                float(double(wave_re_[tapIndex_(i, n)]) - (mean_re + slope_re * k));
            wave_im_[tapIndex_(i, n)] =
                float(double(wave_im_[tapIndex_(i, n)]) - (mean_im + slope_im * k));
        }
    }

    void buildWaveletBank_() {
        wave_re_.fill(0.0f);
        wave_im_.fill(0.0f);
        wave_half_.fill(0);
        wave_gain_onesided_hz_.fill(1.0);
        wave_valid_win_rms2_.fill(1.0);

        const double halfMax = double((Nblock - 1) / 2);
        const double sigma_t_max = (fs > 0.0)
            ? (halfMax / (WAVELET_SUPPORT_SIGMA * fs))
            : 0.0;

        for (int i = 0; i < Nfreq; ++i) {
            const double f0 = freqs_[i];
            if (!(f0 > 0.0) || !(fs > 0.0)) {
                wave_gain_onesided_hz_[i] = 1.0;
                wave_half_[i] = WAVELET_MIN_HALF;
                wave_valid_win_rms2_[i] = 1.0;
                continue;
            }

            double sigma_t = WAVELET_CYCLES_TARGET / (2.0 * M_PI * f0);
            sigma_t = std::min(sigma_t, std::max(1e-6, sigma_t_max));

            int half = int(std::ceil(WAVELET_SUPPORT_SIGMA * sigma_t * fs));
            half = std::clamp(half, WAVELET_MIN_HALF, int(halfMax));
            sigma_t = double(half) / (WAVELET_SUPPORT_SIGMA * fs);

            const int L = 2 * half + 1;
            wave_half_[i] = half;

            {
                const int n0 = half;
                const int n1 = Nblock - half;

                double win_sumsq_valid = 0.0;
                int M_valid = 0;
                for (int n = n0; n < n1; ++n) {
                    const double wv = window_[n];
                    win_sumsq_valid += wv * wv;
                    ++M_valid;
                }

                wave_valid_win_rms2_[i] =
                    (M_valid > 0) ? (win_sumsq_valid / double(M_valid)) : 1.0;
            }

            const double w0 = 2.0 * M_PI * f0;
            const double C0 = std::exp(-0.5 * (w0 * sigma_t) * (w0 * sigma_t));

            for (int n = 0; n < L; ++n) {
                const int k = n - half;
                const double t = double(k) / fs;
                const double u = t / std::max(1e-9, sigma_t);
                const double ga = std::exp(-0.5 * u * u);
                const double phi = 2.0 * M_PI * f0 * t;

                double re = ga * std::cos(phi);
                double im = ga * std::sin(phi);
                re -= (C0 * ga);

                wave_re_[tapIndex_(i, n)] = float(re);
                wave_im_[tapIndex_(i, n)] = float(im);
            }

            orthogonalize_wavelet_to_dc_and_ramp_(i, L, half);

            double H0r = 0.0;
            double H0i = 0.0;
            for (int n = 0; n < L; ++n) {
                const int k = n - half;
                const double t = double(k) / fs;
                const double phi = 2.0 * M_PI * f0 * t;
                const double c = std::cos(phi);
                const double s = std::sin(phi);
                const double re = double(wave_re_[tapIndex_(i, n)]);
                const double im = double(wave_im_[tapIndex_(i, n)]);

                H0r += re * c + im * s;
                H0i += im * c - re * s;
            }

            const double H0mag = std::sqrt(H0r * H0r + H0i * H0i);
            const double scale = (H0mag > 1e-12) ? (1.0 / H0mag) : 1.0;

            for (int n = 0; n < L; ++n) {
                wave_re_[tapIndex_(i, n)] = float(double(wave_re_[tapIndex_(i, n)]) * scale);
                wave_im_[tapIndex_(i, n)] = float(double(wave_im_[tapIndex_(i, n)]) * scale);
            }

            double tap_energy = 0.0;
            for (int n = 0; n < L; ++n) {
                const double re = double(wave_re_[tapIndex_(i, n)]);
                const double im = double(wave_im_[tapIndex_(i, n)]);
                tap_energy += re * re + im * im;
            }

            wave_gain_onesided_hz_[i] = std::max(0.5 * fs * tap_energy, 1e-12);
        }
    }

    void computeSpectrum() {
        constexpr int N = Nblock;

        std::array<double, Nblock> x{};
        int idx = writeIndex;
        for (int n = 0; n < N; ++n) {
            x[n] = buffer_[idx];
            idx = (idx + 1) % Nblock;
        }

        double sumx = 0.0, sumn = 0.0, sumn2 = 0.0, sumnx = 0.0;
        for (int n = 0; n < N; ++n) {
            sumx += x[n];
            sumn += n;
            sumn2 += double(n) * double(n);
            sumnx += double(n) * x[n];
        }

        const double denom = double(N) * sumn2 - sumn * sumn;
        const double a_lin = (denom != 0.0)
            ? (double(N) * sumnx - sumn * sumx) / denom
            : 0.0;
        const double b_lin = (sumx - a_lin * sumn) / double(N);

        std::array<double, Nblock> xw{};
        for (int n = 0; n < N; ++n) {
            const double xdet = x[n] - (a_lin * n + b_lin);
            xw[n] = xdet * window_[n];
        }

        const double Tblk = (fs > 0.0) ? (double(N) / fs) : 0.0;
        const double f_blk = (Tblk > 0.0) ? (1.0 / (6.0 * Tblk)) : 0.0;

        std::array<double, Nfreq> S_aa_true_arr{};

        for (int i = 0; i < Nfreq; ++i) {
            const double f = freqs_[i];
            const int half = wave_half_[i];
            const int L = 2 * half + 1;

            double pwr = 0.0;
            int M = 0;

            for (int n = half; n < N - half; ++n) {
                double yr = 0.0, yi = 0.0;

                for (int tn = 0; tn < L; ++tn) {
                    const int k = tn - half;
                    const double xv = xw[n - k];
                    const double re = double(wave_re_[tapIndex_(i, tn)]);
                    const double im = double(wave_im_[tapIndex_(i, tn)]);
                    yr += re * xv;
                    yi += im * xv;
                }

                pwr += (yr * yr + yi * yi);
                ++M;
            }

            const double var_out = (M > 0) ? (pwr / double(M)) : 0.0;
            const double inv_win_rms2 = 1.0 / std::max(wave_valid_win_rms2_[i], 1e-12);

            const double S_aa_meas =
                std::max(0.0,
                    (var_out * inv_win_rms2) /
                    std::max(wave_gain_onesided_hz_[i], 1e-12));

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

        WaveSpectrumShared::finalize_displacement_spectrum_wavelet_inplace<Nfreq>(
            lastSpectrum_, S_aa_true_arr, freqs_, last_lowfreq_cut_hz_); 
    }

    double fs_raw = 0.0;
    double fs = 0.0;
    int decimFactor = 1;
    bool hannEnabled = true;

    double reg_f0_hz = 0.008;
    double hp_f0_hz = 0.025;

    std::array<double, Nfreq> freqs_{};
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

    static constexpr double WAVELET_CYCLES_TARGET = 6.0;
    static constexpr double WAVELET_SUPPORT_SIGMA = 3.0;
    static constexpr int WAVELET_MIN_HALF = 8;

    std::array<float, Nfreq * Nblock> wave_re_{};
    std::array<float, Nfreq * Nblock> wave_im_{};
    std::array<int, Nfreq> wave_half_{};
    std::array<double, Nfreq> wave_gain_onesided_hz_{};
    std::array<double, Nfreq> wave_valid_win_rms2_{};
};
