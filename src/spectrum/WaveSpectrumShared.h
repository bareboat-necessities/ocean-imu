#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace WaveSpectrumShared {

// -----------------------------
// Common small biquad
// -----------------------------
struct Biquad {
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;
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

static constexpr double BIQUAD_Q = 0.707;

// -----------------------------
// Basic helpers
// -----------------------------
inline double safe_pos(double x, double floor = 1e-12) {
    return std::max(x, floor);
}

inline double lerp(double a, double b, double t) {
    t = std::clamp(t, 0.0, 1.0);
    return a + (b - a) * t;
}

inline double peak_regime_from_fpk(double fpk_hz) {
    // 0 -> long / low-frequency sea
    // 1 -> shorter / higher-frequency sea
    return std::clamp((fpk_hz - 0.08) / (0.35 - 0.08), 0.0, 1.0);
}

inline double ema_alpha_for_f(double f,
                              double fmin,
                              double fmax,
                              double alpha_low,
                              double alpha_high) {
    double t = (f - fmin) / std::max(1e-12, (fmax - fmin));
    t = std::clamp(t, 0.0, 1.0);
    return alpha_low + (alpha_high - alpha_low) * t;
}

// -----------------------------
// Frequency grid
// -----------------------------
template<int Nfreq>
inline void build_log_frequency_grid(std::array<double, Nfreq>& freqs,
                                     std::array<double, Nfreq + 1>& f_edges,
                                     std::array<double, Nfreq>& df,
                                     double f_min = 0.04,
                                     double f_max = 1.2) {
    for (int i = 0; i <= Nfreq; ++i) {
        const double t = double(i) / double(Nfreq);
        f_edges[i] = f_min * std::pow(f_max / f_min, t);
    }

    for (int i = 0; i < Nfreq; ++i) {
        freqs[i] = std::sqrt(f_edges[i] * f_edges[i + 1]);
        df[i] = f_edges[i + 1] - f_edges[i];
    }
}

// -----------------------------
// Biquad design / response
// -----------------------------
inline void designLowpassBiquad(Biquad& bq, double f_cut, double Fs) {
    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / BIQUAD_Q + K * K);

    bq.b0 = K * K * norm;
    bq.b1 = 2.0 * bq.b0;
    bq.b2 = bq.b0;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / BIQUAD_Q + K * K) * norm;
    bq.reset();
}

inline void designHighpassBiquad(Biquad& bq, double f_cut, double Fs) {
    if (f_cut <= 0.0) {
        bq.b0 = 1.0; bq.b1 = 0.0; bq.b2 = 0.0;
        bq.a1 = 0.0; bq.a2 = 0.0;
        bq.reset();
        return;
    }

    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / BIQUAD_Q + K * K);

    bq.b0 = 1.0 * norm;
    bq.b1 = -2.0 * norm;
    bq.b2 = 1.0 * norm;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / BIQUAD_Q + K * K) * norm;
    bq.reset();
}

inline double biquad_mag2_raw(const Biquad& bq, double omega_raw) {
    const double c1 = std::cos(omega_raw);
    const double s1 = std::sin(omega_raw);
    const double c2 = std::cos(2.0 * omega_raw);
    const double s2 = std::sin(2.0 * omega_raw);

    const double num_re = bq.b0 + bq.b1 * c1 + bq.b2 * c2;
    const double num_im = -(bq.b1 * s1 + bq.b2 * s2);

    const double den_re = 1.0 + bq.a1 * c1 + bq.a2 * c2;
    const double den_im = -(bq.a1 * s1 + bq.a2 * s2);

    const double num2 = num_re * num_re + num_im * num_im;
    const double den2 = den_re * den_re + den_im * den_im;

    return num2 / std::max(den2, 1e-16);
}

// -----------------------------
// Accel-domain peak finder
// -----------------------------
template<int Nfreq>
inline int find_accel_peak_index(const std::array<double, Nfreq>& S_aa_true_arr,
                                 const std::array<double, Nfreq>& freqs,
                                 double lowfreq_guard_hz) {
    if (Nfreq < 3) return 0;

    std::array<double, Nfreq> Eaa{};
    std::array<double, Nfreq> Eaa_s{};

    for (int i = 0; i < Nfreq; ++i) {
        Eaa[i] = std::max(0.0, S_aa_true_arr[i]) * std::max(freqs[i], 1e-12);
    }

    Eaa_s[0] = 0.75 * Eaa[0] + 0.25 * Eaa[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Eaa_s[i] = 0.25 * Eaa[i - 1] + 0.50 * Eaa[i] + 0.25 * Eaa[i + 1];
    }
    Eaa_s[Nfreq - 1] = 0.25 * Eaa[Nfreq - 2] + 0.75 * Eaa[Nfreq - 1];

    int k0 = 0;
    const double f0 = std::max(lowfreq_guard_hz, 1.03 * freqs[0]);
    while (k0 + 1 < Nfreq && freqs[k0 + 1] <= f0) ++k0;

    int k_peak = k0;
    double v_peak = -1.0;
    for (int i = k0; i < Nfreq; ++i) {
        if (Eaa_s[i] > v_peak) {
            v_peak = Eaa_s[i];
            k_peak = i;
        }
    }
    return k_peak;
}

// -----------------------------
// Estimator-specific cutoff laws
// -----------------------------
template<int Nfreq>
inline double estimate_lowfreq_cut_from_accel_goertzel(
    const std::array<double, Nfreq>& S_aa_true,
    const std::array<double, Nfreq>& freqs,
    double Tblk,
    double hp_f0_hz)
{
    if (Nfreq < 4) return 0.0;

    std::array<double, Nfreq> E{};
    std::array<double, Nfreq> Es{};

    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, S_aa_true[i]) * std::max(freqs[i], 1e-12);
    }

    Es[0] = 0.75 * E[0] + 0.25 * E[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Es[i] = 0.25 * E[i - 1] + 0.50 * E[i] + 0.25 * E[i + 1];
    }
    Es[Nfreq - 1] = 0.25 * E[Nfreq - 2] + 0.75 * E[Nfreq - 1];

    const double f_floor_hz = std::max({
        1.05 * freqs[0],
        (Tblk > 1e-12) ? (1.25 / Tblk) : 0.0,
        1.05 * hp_f0_hz
    });

    int i_floor = 0;
    while (i_floor + 1 < Nfreq && freqs[i_floor + 1] < f_floor_hz) ++i_floor;

    int i_peak = std::max(1, i_floor);
    double e_peak = Es[i_peak];
    for (int i = i_peak + 1; i < Nfreq - 1; ++i) {
        if (Es[i] > e_peak) {
            e_peak = Es[i];
            i_peak = i;
        }
    }
    if (!(e_peak > 0.0) || i_peak <= i_floor) return f_floor_hz;

    const double fpk = freqs[i_peak];
    const double rpk = peak_regime_from_fpk(fpk);

    int i_valley = i_floor;
    double e_valley = Es[i_floor];
    for (int i = i_floor + 1; i < i_peak; ++i) {
        if (Es[i] < e_valley) {
            e_valley = Es[i];
            i_valley = i;
        }
    }

    const bool left_edge_raised =
        (E[0] > 1.18 * std::max(e_valley, 1e-18)) ||
        (Nfreq > 2 && E[1] > 1.08 * std::max(E[2], 1e-18));

    const bool valley_is_good =
        (i_valley > i_floor) &&
        (e_valley < 0.76 * e_peak) &&
        (freqs[i_valley] < 0.78 * fpk);

    double f_candidate;
    if (left_edge_raised && valley_is_good) {
        f_candidate = std::min(freqs[i_valley], lerp(0.52, 0.62, rpk) * fpk);
    } else {
        f_candidate = lerp(0.24, 0.34, rpk) * fpk;
    }

    return std::max(f_floor_hz, f_candidate);
}

template<int Nfreq>
inline double estimate_lowfreq_cut_from_accel_wavelet(
    const std::array<double, Nfreq>& S_aa_true,
    const std::array<double, Nfreq>& freqs,
    double Tblk,
    double hp_f0_hz)
{
    if (Nfreq < 4) return 0.0;

    std::array<double, Nfreq> E{};
    std::array<double, Nfreq> Es{};

    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, S_aa_true[i]) * std::max(freqs[i], 1e-12);
    }

    Es[0] = 0.75 * E[0] + 0.25 * E[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Es[i] = 0.25 * E[i - 1] + 0.50 * E[i] + 0.25 * E[i + 1];
    }
    Es[Nfreq - 1] = 0.25 * E[Nfreq - 2] + 0.75 * E[Nfreq - 1];

    const double f_floor_hz = std::max({
        1.04 * freqs[0],
        (Tblk > 1e-12) ? (1.20 / Tblk) : 0.0,
        1.03 * hp_f0_hz
    });

    int i_floor = 0;
    while (i_floor + 1 < Nfreq && freqs[i_floor + 1] < f_floor_hz) ++i_floor;

    int i_peak = std::max(1, i_floor);
    double e_peak = Es[i_peak];
    for (int i = i_peak + 1; i < Nfreq - 1; ++i) {
        if (Es[i] > e_peak) {
            e_peak = Es[i];
            i_peak = i;
        }
    }
    if (!(e_peak > 0.0) || i_peak <= i_floor) return f_floor_hz;

    const double fpk = freqs[i_peak];
    const double rpk = peak_regime_from_fpk(fpk);

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
        (e_valley < 0.80 * e_peak) &&
        (freqs[i_valley] < 0.70 * fpk);

    double f_candidate;
    if (left_edge_raised && valley_is_good) {
        f_candidate = std::min(freqs[i_valley], lerp(0.45, 0.55, rpk) * fpk);
    } else {
        f_candidate = lerp(0.18, 0.28, rpk) * fpk;
    }

    return std::max(f_floor_hz, f_candidate);
}

// -----------------------------
// Estimator-specific taper laws
// -----------------------------
inline double lowfreq_taper_goertzel(double f_hz,
                                     double f_cut_hz,
                                     double fpk_hz) {
    if (!(f_cut_hz > 0.0)) return 1.0;

    const double rpk = peak_regime_from_fpk(fpk_hz);
    const double start_mul = lerp(0.62, 0.54, rpk);
    const double full_mul  = lerp(1.08, 1.18, rpk);
    const double taper_pow = lerp(1.05, 1.25, rpk);

    if (f_hz >= full_mul * f_cut_hz) return 1.0;

    if (f_hz <= start_mul * f_cut_hz) {
        const double r = std::max(f_hz / std::max(start_mul * f_cut_hz, 1e-12), 0.0);
        return std::pow(r, taper_pow);
    }

    const double x = (f_hz - start_mul * f_cut_hz) /
                     std::max((full_mul - start_mul) * f_cut_hz, 1e-12);
    const double t = std::clamp(x, 0.0, 1.0);
    const double s = t * t * (3.0 - 2.0 * t);

    const double r0 = std::max(f_hz / std::max(start_mul * f_cut_hz, 1e-12), 0.0);
    const double g0 = std::pow(std::min(r0, 1.0), taper_pow);

    return g0 + (1.0 - g0) * s;
}

inline double lowfreq_taper_wavelet(double f_hz,
                                    double f_cut_hz,
                                    double fpk_hz) {
    if (!(f_cut_hz > 0.0)) return 1.0;

    const double rpk = peak_regime_from_fpk(fpk_hz);
    const double start_mul = lerp(0.68, 0.58, rpk);
    const double full_mul  = lerp(1.08, 1.16, rpk);
    const double taper_pow = lerp(1.00, 1.18, rpk);

    if (f_hz >= full_mul * f_cut_hz) return 1.0;

    if (f_hz <= start_mul * f_cut_hz) {
        const double r = std::max(f_hz / std::max(start_mul * f_cut_hz, 1e-12), 0.0);
        return std::pow(r, taper_pow);
    }

    const double x = (f_hz - start_mul * f_cut_hz) /
                     std::max((full_mul - start_mul) * f_cut_hz, 1e-12);
    const double t = std::clamp(x, 0.0, 1.0);
    const double s = t * t * (3.0 - 2.0 * t);

    const double r0 = std::max(f_hz / std::max(start_mul * f_cut_hz, 1e-12), 0.0);
    const double g0 = std::pow(std::min(r0, 1.0), taper_pow);

    return g0 + (1.0 - g0) * s;
}

// -----------------------------
// Peak-preserving smoother
// -----------------------------
template<int Nfreq, typename SpectrumLike>
inline void smooth_logfreq_3tap_inplace_custom(
    SpectrumLike& spectrum,
    const std::array<double, Nfreq>& freqs,
    double neighbor_k,
    double min_center,
    double peak_keep)
{
    if (Nfreq < 3) return;

    std::array<double, Nfreq> in{};
    std::array<double, Nfreq> out{};

    for (int i = 0; i < Nfreq; ++i) {
        in[i] = std::max(0.0, double(spectrum[i]));
        out[i] = in[i];
    }

    for (int i = 0; i < Nfreq; ++i) {
        const double xim1 = std::log(std::max((i > 0) ? freqs[i - 1] : freqs[i], 1e-12));
        const double xi   = std::log(std::max(freqs[i], 1e-12));
        const double xip1 = std::log(std::max((i < Nfreq - 1) ? freqs[i + 1] : freqs[i], 1e-12));

        const double dL = std::max(0.0, xi - xim1);
        const double dR = std::max(0.0, xip1 - xi);

        double wL = neighbor_k * dL;
        double wR = neighbor_k * dR;
        double wC = std::max(min_center, 1.0 - (wL + wR));
        const double W = wL + wC + wR;

        wL /= W; wC /= W; wR /= W;

        const double Sm1 = (i > 0) ? in[i - 1] : in[i];
        const double Si  = in[i];
        const double Sp1 = (i < Nfreq - 1) ? in[i + 1] : in[i];

        out[i] = wL * Sm1 + wC * Si + wR * Sp1;
    }

    for (int i = 1; i < Nfreq - 1; ++i) {
        const double li = in[i - 1];
        const double ci = in[i];
        const double ri = in[i + 1];

        const bool strong_local_peak =
            (ci > 1.10 * std::max(li, 1e-18)) &&
            (ci > 1.10 * std::max(ri, 1e-18));

        if (strong_local_peak) {
            out[i] = std::max(out[i], peak_keep * ci);
        }
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = std::max(0.0, out[i]);
    }
}

// -----------------------------
// Gentle low-frequency suppressor
// -----------------------------
template<int Nfreq, typename SpectrumLike>
inline void suppress_lowfreq_from_cut_inplace_custom(
    SpectrumLike& spectrum,
    const std::array<double, Nfreq>& freqs,
    double f_cut_hz,
    double f_eff_mul,
    double shape_pow,
    double allow_regrowth)
{
    if (Nfreq < 4) return;
    if (!(f_cut_hz > 0.0)) return;

    const double f_eff = f_eff_mul * f_cut_hz;

    std::array<double, Nfreq> E{};
    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, double(spectrum[i])) * std::max(freqs[i], 1e-12);
    }

    int i_cut = 0;
    while (i_cut + 1 < Nfreq && freqs[i_cut + 1] <= f_eff) ++i_cut;
    if (i_cut <= 0) return;

    const double f_ref = std::max(freqs[i_cut], 1e-12);
    const double E_ref = std::max(E[i_cut], 1e-18);

    double prev = E_ref;
    for (int i = i_cut - 1; i >= 0; --i) {
        const double r = std::max(freqs[i] / f_ref, 0.0);
        const double cap = E_ref * std::pow(r, shape_pow);

        double v = std::max(0.0, E[i]);
        v = std::min(v, cap);
        v = std::min(v, allow_regrowth * prev);

        E[i] = v;
        prev = v;
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = E[i] / std::max(freqs[i], 1e-12);
    }
}

// -----------------------------
// Goertzel-only spike veto
// -----------------------------
template<int Nfreq, typename SpectrumLike>
inline void suppress_unsupported_lowfreq_spikes_goertzel_inplace(
    SpectrumLike& spectrum,
    const std::array<double, Nfreq>& S_aa_true_arr,
    const std::array<double, Nfreq>& freqs,
    double lowfreq_cut_hz)
{
    if (Nfreq < 4) return;
    if (!(lowfreq_cut_hz > 0.0)) return;

    const int k_peak_acc = find_accel_peak_index<Nfreq>(S_aa_true_arr, freqs, lowfreq_cut_hz);
    if (k_peak_acc <= 2 || k_peak_acc >= Nfreq) return;

    const double fpk = freqs[k_peak_acc];
    const double rpk = peak_regime_from_fpk(fpk);

    double f_soft = lerp(1.12, 1.30, rpk) * lowfreq_cut_hz;
    f_soft = std::min(f_soft, lerp(0.62, 0.72, rpk) * fpk);

    int k_soft = 0;
    while (k_soft + 1 < Nfreq && freqs[k_soft + 1] <= f_soft) ++k_soft;

    const int k_stop = std::min(k_soft, std::max(1, k_peak_acc - 3));
    if (k_stop < 2 || k_stop + 1 >= Nfreq) return;

    std::array<double, Nfreq> Eeta{};
    std::array<double, Nfreq> Eaa{};
    std::array<double, Nfreq> Eaa_s{};

    for (int i = 0; i < Nfreq; ++i) {
        const double f = std::max(freqs[i], 1e-12);
        Eeta[i] = std::max(0.0, double(spectrum[i])) * f;
        Eaa[i]  = std::max(0.0, S_aa_true_arr[i]) * f;
    }

    Eaa_s[0] = 0.75 * Eaa[0] + 0.25 * Eaa[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Eaa_s[i] = 0.25 * Eaa[i - 1] + 0.50 * Eaa[i] + 0.25 * Eaa[i + 1];
    }
    Eaa_s[Nfreq - 1] = 0.25 * Eaa[Nfreq - 2] + 0.75 * Eaa[Nfreq - 1];

    const int i_ref = std::min(k_stop + 1, Nfreq - 1);
    const double Eref = std::max(Eeta[i_ref], 1e-12);
    const double fref = std::max(freqs[i_ref], 1e-12);
    const double Eaa_ref = std::max(Eaa_s[i_ref], 1e-12);

    const double accel_support_thr = lerp(0.10, 0.20, rpk);
    const double spike_ratio       = lerp(2.20, 1.80, rpk);
    const double env_pow           = lerp(1.60, 1.95, rpk);

    for (int i = 1; i <= k_stop; ++i) {
        const bool local_peak =
            (Eeta[i] > Eeta[i - 1]) &&
            (Eeta[i] > Eeta[i + 1]);

        const double r = std::max(freqs[i], 1e-12) / fref;
        const double env_cap = Eref * std::pow(r, env_pow);
        const double accel_support = Eaa_s[i] / Eaa_ref;

        if (local_peak &&
            accel_support < accel_support_thr &&
            Eeta[i] > spike_ratio * env_cap) {
            Eeta[i] = env_cap;
        }
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = Eeta[i] / std::max(freqs[i], 1e-12);
    }
}

// -----------------------------
// Guarded fp estimate
// -----------------------------
template<int Nfreq, typename SpectrumLike>
inline double estimate_fp_with_guard(const SpectrumLike& spectrum,
                                     const std::array<double, Nfreq>& freqs,
                                     double lowfreq_guard_hz) {
    if (Nfreq == 0) return 0.0;

    int k0 = 0;
    const double f_guard = std::max(lowfreq_guard_hz, freqs[0]);
    while (k0 + 1 < Nfreq && freqs[k0] < f_guard) ++k0;

    int k = k0;
    double vmax = -1.0;
    for (int i = k0; i < Nfreq; ++i) {
        const double s = std::max(0.0, double(spectrum[i]));
        if (s > vmax) {
            vmax = s;
            k = i;
        }
    }

    if (!(vmax > 0.0)) return freqs[k0];
    if (k <= 0 || k >= Nfreq - 1) return freqs[k];

    const double f0 = freqs[k - 1];
    const double f1 = freqs[k];
    const double f2 = freqs[k + 1];
    if (!(f0 > 0.0 && f1 > 0.0 && f2 > 0.0)) return f1;

    const double x0 = std::log(f0);
    const double x1 = std::log(f1);
    const double x2 = std::log(f2);

    const double y0 = std::log(std::max(double(spectrum[k - 1]), 1e-18));
    const double y1 = std::log(std::max(double(spectrum[k]),     1e-18));
    const double y2 = std::log(std::max(double(spectrum[k + 1]), 1e-18));

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

// -----------------------------
// Finalizers
// -----------------------------
template<int Nfreq, typename SpectrumLike>
inline void finalize_displacement_spectrum_goertzel_inplace(
    SpectrumLike& spectrum,
    const std::array<double, Nfreq>& S_aa_true_arr,
    const std::array<double, Nfreq>& freqs,
    double lowfreq_cut_hz)
{
    const int k_peak_acc = find_accel_peak_index<Nfreq>(S_aa_true_arr, freqs, lowfreq_cut_hz);
    const double fpk = freqs[std::clamp(k_peak_acc, 0, Nfreq - 1)];
    const double rpk = peak_regime_from_fpk(fpk);

    smooth_logfreq_3tap_inplace_custom<Nfreq>(
        spectrum, freqs,
        lerp(0.16, 0.24, rpk),
        lerp(0.62, 0.56, rpk),
        0.94);

    suppress_unsupported_lowfreq_spikes_goertzel_inplace<Nfreq>(
        spectrum, S_aa_true_arr, freqs, lowfreq_cut_hz);

    suppress_lowfreq_from_cut_inplace_custom<Nfreq>(
        spectrum, freqs, lowfreq_cut_hz,
        lerp(0.80, 0.88, rpk),
        lerp(1.25, 1.50, rpk),
        1.02);
}

template<int Nfreq, typename SpectrumLike>
inline void finalize_displacement_spectrum_wavelet_inplace(
    SpectrumLike& spectrum,
    const std::array<double, Nfreq>& S_aa_true_arr,
    const std::array<double, Nfreq>& freqs,
    double lowfreq_cut_hz)
{
    const int k_peak_acc = find_accel_peak_index<Nfreq>(S_aa_true_arr, freqs, lowfreq_cut_hz);
    const double fpk = freqs[std::clamp(k_peak_acc, 0, Nfreq - 1)];
    const double rpk = peak_regime_from_fpk(fpk);

    smooth_logfreq_3tap_inplace_custom<Nfreq>(
        spectrum, freqs,
        lerp(0.12, 0.18, rpk),
        lerp(0.72, 0.64, rpk),
        0.97);

    suppress_lowfreq_from_cut_inplace_custom<Nfreq>(
        spectrum, freqs, lowfreq_cut_hz,
        lerp(0.70, 0.80, rpk),
        lerp(1.15, 1.35, rpk),
        1.04);
}

} // namespace WaveSpectrumShared
