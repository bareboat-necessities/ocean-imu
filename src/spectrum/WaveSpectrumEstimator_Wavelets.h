#pragma once

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "spectrum/SpectrumStats.h"

/*
    Copyright 2025-2026, Mikhail Grushinskiy

    Same interface as the original Goertzel estimator, but the block spectrum is computed
    using a complex Morlet wavelet filter bank (CWT-style) evaluated on the same fixed
    log-spaced frequency grid.

    Signal chain (per raw sample):
        raw a_z  ->  HP (2nd)  ->  HP (2nd)  ->  LP (2nd)  ->  downsample by D  ->  block buffer

    Wavelet spectrum estimation (once per full block at decimated fs):
        - Linear detrend the decimated block.
        - Apply Hann window (optional).
        - For each grid bin f_i:
            * Convolve the block with a complex Morlet wavelet centered at f_i.
            * Estimate coefficient power E{|c_i|^2} over valid samples.
            * Convert to one-sided acceleration PSD via:
                  S_aa(f_i) ~= Var_out / G_i
              where G_i is the one-sided gain integral:
                  G_i = integral_0^{fs/2} (|H(f)|^2 + |H(-f)|^2)/2 df
              computed numerically once in the constructor per bin.

        - Deconvolve ONLY the raw-rate HP chain magnitude at f_i.
          This version uses a SOFT inverse, not 1/max(H2, floor), to avoid a fake low-f shelf.
        - Convert to displacement PSD with Tikhonov-regularized 1/w^4:
              S_eta(f) = S_aa_true(f) / ( (w^2 + lambda^2)^2 )

    Improvements in this version:
      - wavelet taps are orthogonalized to both DC and ramp
      - HP deconvolution is softly regularized
      - displacement knee is tied to block length and HP corner
      - low-frequency/DC leak suppression is adaptive and uses E(f)=f*S(f)
*/

template<int Nfreq = 32, int Nblock = 256>
class EIGEN_ALIGN_MAX WaveSpectrumEstimator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static constexpr double g = 9.80665;
    using Vec = Eigen::Matrix<double, Nfreq, 1>;

    struct PMFitResult { double alpha, fp, cost; };

    struct Biquad {
        double b0 = 0.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;
        double z1 = 0.0, z2 = 0.0;

        inline double process(double x) {
            const double y = b0 * x + z1;
            z1 = b1 * x - a1 * y + z2;
            z2 = b2 * x - a2 * y;
            return y;
        }

        inline void reset() {
            z1 = 0.0;
            z2 = 0.0;
        }
    };

    WaveSpectrumEstimator(double fs_raw_ = 200.0,
                          int decimFactor_ = 30,
                          bool hannEnabled_ = true)
        : fs_raw(fs_raw_), decimFactor(decimFactor_), hannEnabled(hannEnabled_)
    {
        fs = fs_raw / decimFactor;

        buildFrequencyGrid();

        double sumsq = 0.0;
        for (int n = 0; n < Nblock; ++n) {
            if (hannEnabled) {
                window_[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (Nblock - 1)));
            } else {
                window_[n] = 1.0;
            }
            sumsq += window_[n] * window_[n];
        }
        window_sum_sq = sumsq;

        reset();

        const double fny_dec = fs_raw_ / (2.0 * decimFactor_);
        const double lp_cutoffHz = 0.32 * fny_dec;
        designLowpassBiquad(lp1_, lp_cutoffHz, fs_raw_);
        designLowpassBiquad(lp2_, lp_cutoffHz, fs_raw_);
        designHighpassBiquad(hp1_, hp_f0_hz, fs_raw_);
        designHighpassBiquad(hp2_, hp_f0_hz, fs_raw_);

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

    double computeHs() const {
        return SpectrumStats::compute_hs<Nfreq>(lastSpectrum_, df_);
    }

    double estimateFp() const {
        return SpectrumStats::estimate_fp<Nfreq>(lastSpectrum_, freqs_);
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
                const double d = safeLog(S_obs[i]) - safeLog(model);
                cost += df * d * d;
            }
            return cost;
        };

        constexpr int N_fp_search = 32;
        constexpr double fp_min = 0.05;
        constexpr double fp_transition = 0.10;
        constexpr double fp_max = 1.0;

        std::array<double, N_fp_search> fp_grid{};
        const int n_log = static_cast<int>(N_fp_search * 0.4);
        const int n_lin = N_fp_search - n_log;

        for (int i = 0; i < n_log; ++i) {
            const double t = double(i) / double(n_log - 1);
            fp_grid[i] = fp_min * std::pow(fp_transition / fp_min, t);
        }
        for (int i = 0; i < n_lin; ++i) {
            const double t = double(i) / double(n_lin - 1);
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
        designHighpassBiquad(hp1_, hp_f0_hz, fs_raw);
        designHighpassBiquad(hp2_, hp_f0_hz, fs_raw);
    }

private:
    inline double safeLog(double v) const {
        return std::log(std::max(v, 1e-18));
    }

    inline double alpha_for_f(double f) const {
        const double fmin = freqs_[0];
        const double fmax = freqs_[Nfreq - 1];
        double t = (f - fmin) / std::max(1e-12, (fmax - fmin));
        t = std::clamp(t, 0.0, 1.0);
        return ema_alpha_low + (ema_alpha_high - ema_alpha_low) * t;
    }

    void smooth_logfreq_3tap() {
        if (Nfreq < 3) return;

        Eigen::Matrix<double, Nfreq, 1> S = lastSpectrum_;
        Eigen::Matrix<double, Nfreq, 1> Sout = S;

        auto wpair = [&](int i) {
            const double eps = 1e-12;
            const double x_im1 = (i > 0)
                ? std::log(std::max(freqs_[i - 1], eps))
                : std::log(std::max(freqs_[i], eps));
            const double x_i = std::log(std::max(freqs_[i], eps));
            const double x_ip1 = (i < Nfreq - 1)
                ? std::log(std::max(freqs_[i + 1], eps))
                : std::log(std::max(freqs_[i], eps));

            const double dL = std::max(0.0, x_i - x_im1);
            const double dR = std::max(0.0, x_ip1 - x_i);

            const double k = 0.35;
            const double minC = 0.40;

            const double wL = k * dL;
            const double wR = k * dR;
            const double wC = std::max(minC, 1.0 - (wL + wR));
            const double W = wL + wC + wR;

            return std::array<double, 3>{ wL / W, wC / W, wR / W };
        };

        for (int i = 0; i < Nfreq; ++i) {
            const auto w = wpair(i);
            const double Sm1 = (i > 0) ? S[i - 1] : S[i];
            const double Sp1 = (i < Nfreq - 1) ? S[i + 1] : S[i];
            Sout[i] = w[0] * Sm1 + w[1] * S[i] + w[2] * Sp1;
        }

        lastSpectrum_ = Sout;
    }

    // Adaptive low-frequency/DC leak suppression using E(f)=f*S(f).
    // This is more robust than searching the raw peak of S(f), because S(f)
    // is exactly what gets exaggerated by 1/w^4 at the lowest bins.
    void suppress_lowfreq_dc_leak_() {
        if (Nfreq < 4 || fs <= 0.0) return;

        Eigen::Matrix<double, Nfreq, 1> E;
        Eigen::Matrix<double, Nfreq, 1> Es;

        for (int i = 0; i < Nfreq; ++i) {
            E[i] = std::max(0.0, (double)lastSpectrum_[i]) * std::max(freqs_[i], 1e-12);
        }

        Es[0] = 0.75 * E[0] + 0.25 * E[1];
        for (int i = 1; i < Nfreq - 1; ++i) {
            Es[i] = 0.25 * E[i - 1] + 0.50 * E[i] + 0.25 * E[i + 1];
        }
        Es[Nfreq - 1] = 0.25 * E[Nfreq - 2] + 0.75 * E[Nfreq - 1];

        const double Tblk = double(Nblock) / std::max(fs, 1e-12);
        const double f_floor_1 = 1.10 * freqs_[0];
        const double f_floor_2 = 2.0 / std::max(Tblk, 1e-12);
        const double f_floor_3 = 1.15 * hp_f0_hz;
        const double f_floor_hz = std::max(f_floor_1, std::max(f_floor_2, f_floor_3));

        int i_floor = 0;
        while (i_floor + 1 < Nfreq && freqs_[i_floor + 1] < f_floor_hz) {
            ++i_floor;
        }

        int i_peak = std::max(1, i_floor);
        if (i_peak >= Nfreq) i_peak = Nfreq - 1;

        double e_peak = Es[i_peak];
        for (int i = i_peak + 1; i < Nfreq - 1; ++i) {
            if (Es[i] > e_peak) {
                e_peak = Es[i];
                i_peak = i;
            }
        }

        if (!(e_peak > 0.0) || i_peak <= i_floor) return;

        int i_valley = i_floor;
        double e_valley = Es[i_floor];
        for (int i = i_floor + 1; i < i_peak; ++i) {
            if (Es[i] < e_valley) {
                e_valley = Es[i];
                i_valley = i;
            }
        }

        int i_cut = i_floor;
        const bool valley_is_good =
            (i_valley > i_floor) &&
            (e_valley < 0.72 * e_peak) &&
            (freqs_[i_valley] < 0.85 * freqs_[i_peak]);

        if (valley_is_good) {
            i_cut = i_valley;
        } else {
            const double f_rel = std::max(f_floor_hz, 0.42 * freqs_[i_peak]);
            for (int i = i_floor; i < i_peak; ++i) {
                if (freqs_[i] <= f_rel) i_cut = i;
            }
        }

        if (i_cut <= 0) return;

        const double f_ref = std::max(freqs_[i_cut], 1e-12);
        const double E_ref = std::max(E[i_cut], 1e-18);

        double left_width = std::log(std::max(freqs_[i_peak], 1e-12) / f_ref);
        left_width = std::clamp(left_width, 0.05, 1.5);

        const double shape_pow = std::clamp(3.6 - 1.1 * left_width, 2.2, 3.6);

        double prev = E_ref;
        for (int i = i_cut - 1; i >= 0; --i) {
            const double r = std::max(0.0, freqs_[i] / f_ref);
            const double shape_cap = E_ref * std::pow(r, shape_pow);

            double v = std::max(0.0, (double)E[i]);
            v = std::min(v, shape_cap);
            v = std::min(v, prev);

            E[i] = v;
            prev = v;
        }

        for (int i = 0; i < Nfreq; ++i) {
            lastSpectrum_[i] = E[i] / std::max(freqs_[i], 1e-12);
        }
    }

    void orthogonalize_wavelet_to_dc_and_ramp_(int i, int L, int half) {
        double sum_re = 0.0, sum_im = 0.0;
        double sumk_re = 0.0, sumk_im = 0.0;
        double sumk2 = 0.0;

        for (int n = 0; n < L; ++n) {
            const double k = double(n - half);
            const double re = (double)wave_re_[tapIndex_(i, n)];
            const double im = (double)wave_im_[tapIndex_(i, n)];

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
                (float)((double)wave_re_[tapIndex_(i, n)] - (mean_re + slope_re * k));
            wave_im_[tapIndex_(i, n)] =
                (float)((double)wave_im_[tapIndex_(i, n)] - (mean_im + slope_im * k));
        }
    }

    void buildFrequencyGrid() {
        constexpr double f_min = 0.04;
        constexpr double f_max = 1.2;

        for (int i = 0; i <= Nfreq; ++i) {
            const double t = double(i) / double(Nfreq);
            f_edges_[i] = f_min * std::pow(f_max / f_min, t);
        }

        for (int i = 0; i < Nfreq; ++i) {
            freqs_[i] = std::sqrt(f_edges_[i] * f_edges_[i + 1]);
            df_[i] = f_edges_[i + 1] - f_edges_[i];
        }
    }

    void designLowpassBiquad(Biquad& bq, double f_cut, double Fs) {
        const double Fc = f_cut / Fs;
        const double K = std::tan(M_PI * Fc);
        const double norm = 1.0 / (1.0 + K / Q + K * K);

        bq.b0 = K * K * norm;
        bq.b1 = 2.0 * bq.b0;
        bq.b2 = bq.b0;
        bq.a1 = 2.0 * (K * K - 1.0) * norm;
        bq.a2 = (1.0 - K / Q + K * K) * norm;
    }

    void designHighpassBiquad(Biquad& bq, double f_cut, double Fs) {
        if (f_cut <= 0.0) {
            bq.b0 = 1.0; bq.b1 = 0.0; bq.b2 = 0.0;
            bq.a1 = 0.0; bq.a2 = 0.0;
            bq.reset();
            return;
        }

        const double Fc = f_cut / Fs;
        const double K = std::tan(M_PI * Fc);
        const double norm = 1.0 / (1.0 + K / Q + K * K);

        bq.b0 = 1.0 * norm;
        bq.b1 = -2.0 * norm;
        bq.b2 = 1.0 * norm;
        bq.a1 = 2.0 * (K * K - 1.0) * norm;
        bq.a2 = (1.0 - K / Q + K * K) * norm;
        bq.reset();
    }

    inline double biquad_mag2_raw(const Biquad& bq, double Omega_raw) const {
        const double c1 = std::cos(Omega_raw);
        const double s1 = std::sin(Omega_raw);
        const double c2 = std::cos(2.0 * Omega_raw);
        const double s2 = std::sin(2.0 * Omega_raw);

        const double num_re = bq.b0 + bq.b1 * c1 + bq.b2 * c2;
        const double num_im = -(bq.b1 * s1 + bq.b2 * s2);

        const double den_re = 1.0 + bq.a1 * c1 + bq.a2 * c2;
        const double den_im = -(bq.a1 * s1 + bq.a2 * s2);

        const double num2 = num_re * num_re + num_im * num_im;
        const double den2 = den_re * den_re + den_im * den_im;

        return num2 / std::max(den2, 1e-16);
    }

    static constexpr int WAVELET_NFFT =
        (Nblock <= 256) ? 1024 :
        (Nblock <= 512) ? 2048 : 4096;

    static constexpr double WAVELET_CYCLES_TARGET = 6.0;
    static constexpr double WAVELET_SUPPORT_SIGMA = 3.0;
    static constexpr int    WAVELET_MIN_HALF = 8;

    std::array<float,  Nfreq * Nblock> wave_re_{};
    std::array<float,  Nfreq * Nblock> wave_im_{};
    std::array<int,    Nfreq>          wave_half_{};
    std::array<double, Nfreq>          wave_gain_onesided_hz_{};

    inline int tapIndex_(int i, int n) const {
        return i * Nblock + n;
    }

    void buildWaveletBank_() {
        wave_re_.fill(0.0f);
        wave_im_.fill(0.0f);
        wave_half_.fill(0);
        wave_gain_onesided_hz_.fill(1.0);

        const double halfMax = double((Nblock - 1) / 2);
        const double sigma_t_max = (fs > 0.0)
            ? (halfMax / (WAVELET_SUPPORT_SIGMA * fs))
            : 0.0;

        std::array<double, WAVELET_NFFT> Hre{};
        std::array<double, WAVELET_NFFT> Him{};

        for (int i = 0; i < Nfreq; ++i) {
            const double f0 = freqs_[i];

            if (!(f0 > 0.0) || !(fs > 0.0)) {
                wave_gain_onesided_hz_[i] = 1.0;
                wave_half_[i] = WAVELET_MIN_HALF;
                continue;
            }

            double sigma_t = WAVELET_CYCLES_TARGET / (2.0 * M_PI * f0);
            sigma_t = std::min(sigma_t, std::max(1e-6, sigma_t_max));

            int half = int(std::ceil(WAVELET_SUPPORT_SIGMA * sigma_t * fs));
            half = std::clamp(half, WAVELET_MIN_HALF, int(halfMax));
            sigma_t = double(half) / (WAVELET_SUPPORT_SIGMA * fs);

            const int L = 2 * half + 1;
            wave_half_[i] = half;

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

                wave_re_[tapIndex_(i, n)] = (float)re;
                wave_im_[tapIndex_(i, n)] = (float)im;
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
                const double re = (double)wave_re_[tapIndex_(i, n)];
                const double im = (double)wave_im_[tapIndex_(i, n)];

                H0r += re * c + im * s;
                H0i += im * c - re * s;
            }

            const double H0mag = std::sqrt(H0r * H0r + H0i * H0i);
            const double scale = (H0mag > 1e-12) ? (1.0 / H0mag) : 1.0;

            for (int n = 0; n < L; ++n) {
                wave_re_[tapIndex_(i, n)] = (float)((double)wave_re_[tapIndex_(i, n)] * scale);
                wave_im_[tapIndex_(i, n)] = (float)((double)wave_im_[tapIndex_(i, n)] * scale);
            }

            Hre.fill(0.0);
            Him.fill(0.0);

            const double df = fs / double(WAVELET_NFFT);

            for (int k = 0; k < WAVELET_NFFT; ++k) {
                const double w = 2.0 * M_PI * double(k) / double(WAVELET_NFFT);
                const double cw = std::cos(w);
                const double sw = std::sin(w);

                double cr = 1.0;
                double ci = 0.0;
                double sr = 0.0;
                double si = 0.0;

                for (int n = 0; n < L; ++n) {
                    const double re = (double)wave_re_[tapIndex_(i, n)];
                    const double im = (double)wave_im_[tapIndex_(i, n)];

                    sr += re * cr + im * ci;
                    si += im * cr - re * ci;

                    const double crn = cr * cw - ci * sw;
                    const double cin = cr * sw + ci * cw;
                    cr = crn;
                    ci = cin;
                }

                Hre[k] = sr;
                Him[k] = si;
            }

            double gain = 0.0;
            for (int k = 0; k <= WAVELET_NFFT / 2; ++k) {
                const int kneg = (k == 0 || k == WAVELET_NFFT / 2) ? k : (WAVELET_NFFT - k);
                const double mag2_pos = Hre[k] * Hre[k] + Him[k] * Him[k];
                const double mag2_neg = Hre[kneg] * Hre[kneg] + Him[kneg] * Him[kneg];
                gain += 0.5 * (mag2_pos + mag2_neg) * df;
            }

            wave_gain_onesided_hz_[i] = std::max(gain, 1e-12);
        }
    }

    void computeSpectrum() {
        const int N = Nblock;

        std::array<double, Nblock> x{};
        const int startIdx = writeIndex;
        int idx = startIdx;
        for (int n = 0; n < N; ++n) {
            x[n] = buffer_[idx];
            idx = (idx + 1) % Nblock;
        }

        double sumx = 0.0, sumn = 0.0, sumn2 = 0.0, sumnx = 0.0;
        for (int n = 0; n < N; ++n) {
            sumx  += x[n];
            sumn  += n;
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

        const double w_rms2 = (N > 0) ? (window_sum_sq / double(N)) : 1.0;
        const double inv_w_rms2 = 1.0 / std::max(w_rms2, 1e-12);

        const double Tblk = (fs > 0.0) ? (double(N) / fs) : 0.0;
        const double f_blk = (Tblk > 0.0) ? (1.0 / (6.0 * Tblk)) : 0.0;
        const double f_knee = std::max(reg_f0_hz, std::max(0.8 * hp_f0_hz, f_blk));
        const double lam = 2.0 * M_PI * f_knee;

        for (int i = 0; i < Nfreq; ++i) {
            const double f = freqs_[i];
            const int half = wave_half_[i];
            const int L = 2 * half + 1;

            double pwr = 0.0;
            int M = 0;

            for (int n = half; n < N - half; ++n) {
                double yr = 0.0;
                double yi = 0.0;

                for (int tn = 0; tn < L; ++tn) {
                    const int k = tn - half;
                    const double xv = xw[n - k];
                    const double re = (double)wave_re_[tapIndex_(i, tn)];
                    const double im = (double)wave_im_[tapIndex_(i, tn)];
                    yr += re * xv;
                    yi += im * xv;
                }

                pwr += (yr * yr + yi * yi);
                ++M;
            }

            const double var_out = (M > 0) ? (pwr / double(M)) : 0.0;
            const double S_aa_meas =
                std::max(0.0, (var_out * inv_w_rms2) / std::max(wave_gain_onesided_hz_[i], 1e-12));

            const double Omega_raw = 2.0 * M_PI * f / fs_raw;
            const double H2_hp =
                biquad_mag2_raw(hp1_, Omega_raw) *
                biquad_mag2_raw(hp2_, Omega_raw);

            // Soft inverse of HP power response.
            // Avoids the fake low-frequency shelf caused by 1/max(H2, floor).
            constexpr double hp_deconv_reg = 0.10;
            const double inv_hp =
                H2_hp / (H2_hp * H2_hp + hp_deconv_reg * hp_deconv_reg);

            const double S_aa_true = S_aa_meas * inv_hp;

            const double w = 2.0 * M_PI * f;
            const double den = (w * w + lam * lam);
            double S_eta = (den > 0.0) ? (S_aa_true / (den * den)) : 0.0;

            if (!std::isfinite(S_eta) || S_eta < 0.0) S_eta = 0.0;

            if (use_psd_ema) {
                const double a = alpha_for_f(f);
                if (!have_ema) {
                    psd_ema_[i] = S_eta;
                } else {
                    psd_ema_[i] = (1.0 - a) * psd_ema_[i] + a * S_eta;
                }
                lastSpectrum_[i] = psd_ema_[i];
            } else {
                lastSpectrum_[i] = S_eta;
            }
        }

        have_ema = true;
        smooth_logfreq_3tap();
        suppress_lowfreq_dc_leak_();
    }

    double fs_raw = 0.0;
    double fs = 0.0;
    int decimFactor = 1;
    bool hannEnabled = true;

    double reg_f0_hz = 0.008;
    double hp_f0_hz  = 0.025;

    std::array<double, Nfreq>   freqs_{};
    std::array<double, Nfreq + 1> f_edges_{};
    std::array<double, Nfreq>   df_{};

    std::array<double, Nblock> buffer_{};
    std::array<double, Nblock> window_{};
    double window_sum_sq = 1.0;

    Biquad hp1_, hp2_, lp1_, lp2_;

    Eigen::Matrix<double, Nfreq, 1> lastSpectrum_;

    bool use_psd_ema = true;
    double ema_alpha_low  = 0.20;
    double ema_alpha_high = 0.06;
    bool have_ema = false;
    Eigen::Matrix<double, Nfreq, 1> psd_ema_;

    int writeIndex = 0;
    int decimCounter = 0;
    int filledSamples = 0;
    bool isWarm = false;

    static constexpr double Q = 0.707;
};
