#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Literature baseline: the harmonic-mode heave EKF of

    S. Kuchler, C. Pregizer, J. K. Eberharter, K. Schneider, O. Sawodny,
    "Real-Time Estimation of a Ship's Attitude", IFAC World Congress 2011
    (heave part, Sections 3-4).

  Heave is a sum of N_m undamped modes z = sum_j z_j, each with state
  x_j = [z_j, z_j', w_j] and dynamics

      z_j'' = -w_j^2 z_j,   w_j' = 0 (random walk),

  plus one random-walk offset state that absorbs gravity and accelerometer
  bias.  The EKF starts at the first identification; the offset is seeded by
  the fitted constant of the identification buffer.  The measurement is the vertical acceleration

      y = sum_j (-w_j^2 z_j) + x_off + noise.

  The oscillator is discretised exactly for the current w_j (a rotation in the
  (z, z'/w) plane), which is the paper's exact discretisation.

  Identification: the vertical acceleration is block-averaged to fs_ident_hz
  into a ring buffer.  Every ident_period_s an FFT of the Hann-windowed buffer
  gives the acceleration amplitude spectrum; dividing by w^2 gives the heave
  spectrum, whose peaks set the number of modes and their frequencies.  Peaks
  that match no existing mode are appended to the EKF.  Only the new states
  are seeded; existing states keep their estimates.  The seeds (amplitude and
  phase at the current instant) come from a joint least-squares fit of the
  buffer at the peak frequencies.  The paper reads them off the FFT bins;
  the fit is the same quantity without the bin-phase ambiguity.  Observability
  maintenance follows the paper: a mode whose frequency leaves the admissible
  band is removed, and two modes whose frequencies converge are merged.

  Tuning follows the paper's pattern (higher process noise for higher
  frequency modes, slow frequency and offset states, R from the sensor
  datasheet).  Two choices are made here; the paper fixes neither:
  the velocity process-noise intensity of mode j is
  q_j = 4 zeta_q w_j^3 (A_j^2/2), the intensity that holds a damped
  oscillator of damping zeta_q at the identified amplitude A_j, so one
  setting serves every sea state; when all mode slots are occupied, a much
  stronger new peak replaces the weakest mode.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

namespace heave_baselines {

template <typename T = double, int MaxModes = 4>
class KuchlerHeaveEKF {
    static_assert(std::is_floating_point<T>::value,
                  "KuchlerHeaveEKF<T>: T must be a floating-point type.");

public:
    static constexpr int N = 3 * MaxModes + 1;
    static constexpr int OFF = 3 * MaxModes;
    using Vec = Eigen::Matrix<T, N, 1>;
    using Mat = Eigen::Matrix<T, N, N>;

    struct Config {
        // Identification.
        T fs_ident_hz = T(4);
        int fft_size = 512;              // 128 s at 4 Hz, 0.0078 Hz bins
        T ident_period_s = T(30);
        T f_min_hz = T(0.04);
        T f_max_hz = T(1.0);
        T peak_rel_threshold = T(0.25);  // of the largest heave-spectrum peak
        T match_tol = T(0.15);           // relative frequency
        T merge_tol = T(0.06);           // relative frequency
        T replace_ratio = T(2.0);

        // Noise model.
        T accel_noise_std = T(0.0148);   // datasheet white noise, m/s^2
        T zeta_q = T(0.01);
        T omega_rw = T(0.002);           // relative frequency random walk, 1/sqrt(s)
        T offset_rw = T(1e-3);           // m/s^2/sqrt(s)
        T offset_init_std = T(0.02);     // m/s^2, about the fitted buffer mean
        T seed_rel_std = T(0.5);         // initial amplitude uncertainty
        T seed_omega_rel_std = T(0.05);
    };

    KuchlerHeaveEKF() : KuchlerHeaveEKF(Config{}) {}
    explicit KuchlerHeaveEKF(const Config& cfg) : cfg_(cfg) { reset(); }

    void reset() {
        x_.setZero();
        P_.setZero();
        active_.fill(false);
        q_.fill(T(0));
        t_ = T(0);
        block_sum_ = T(0);
        block_n_ = 0;
        ring_.assign(static_cast<size_t>(cfg_.fft_size), T(0));
        ring_head_ = 0;
        ring_count_ = 0;
        since_ident_ = T(0);
        started_ = false;
        identifications_ = 0;
    }

    // a_up: levelled vertical specific force minus gravity, m/s^2, z up.
    T update(T a_up, T dt) {
        if (!(dt > T(0)) || !std::isfinite(a_up)) return heave();
        t_ += dt;

        if (started_) {
            predict_(dt);
            correct_(a_up);
            prune_();
        }

        // Identification buffer.
        block_sum_ += a_up;
        ++block_n_;
        const int block_len = std::max(1, static_cast<int>(std::lround(T(1) / (cfg_.fs_ident_hz * dt))));
        if (block_n_ >= block_len) {
            push_(block_sum_ / T(block_n_));
            block_sum_ = T(0);
            block_n_ = 0;
            block_dt_ = T(block_len) * dt;
            raw_dt_ = dt;
        }
        since_ident_ += dt;
        if (ring_count_ >= cfg_.fft_size &&
            since_ident_ >= cfg_.ident_period_s && block_n_ == 0)
        {
            identify_();
            since_ident_ = T(0);
        }
        return heave();
    }

    T heave() const {
        T z = T(0);
        for (int j = 0; j < MaxModes; ++j) if (active_[j]) z += x_(3 * j);
        return z;
    }
    T heaveRate() const {
        T v = T(0);
        for (int j = 0; j < MaxModes; ++j) if (active_[j]) v += x_(3 * j + 1);
        return v;
    }
    T heaveAccel() const {
        T a = T(0);
        for (int j = 0; j < MaxModes; ++j)
            if (active_[j]) a -= x_(3 * j + 2) * x_(3 * j + 2) * x_(3 * j);
        return a;
    }
    T offset() const { return x_(OFF); }
    int modeCount() const {
        return static_cast<int>(std::count(active_.begin(), active_.end(), true));
    }
    // Frequency of the mode with the largest heave amplitude, rad/s.
    T dominantFrequency() const {
        T best = T(-1), w = T(NAN);
        for (int j = 0; j < MaxModes; ++j) {
            if (!active_[j]) continue;
            const T a = modeAmplitude_(j);
            if (a > best) { best = a; w = x_(3 * j + 2); }
        }
        return w;
    }
    int identifications() const { return identifications_; }

private:
    struct Peak { T w; T amp; T zc; T zs; };

    T modeAmplitude_(int j) const {
        const T w = x_(3 * j + 2);
        const T z = x_(3 * j);
        const T v = x_(3 * j + 1);
        return std::sqrt(z * z + (v / w) * (v / w));
    }

    void predict_(T h) {
        Mat F = Mat::Identity();
        for (int j = 0; j < MaxModes; ++j) {
            if (!active_[j]) continue;
            const int i = 3 * j;
            const T z = x_(i), v = x_(i + 1), w = x_(i + 2);
            const T c = std::cos(w * h), s = std::sin(w * h);
            x_(i)     = z * c + v / w * s;
            x_(i + 1) = -z * w * s + v * c;
            F(i, i) = c;          F(i, i + 1) = s / w;
            F(i, i + 2) = -z * h * s + v * (h * c / w - s / (w * w));
            F(i + 1, i) = -w * s; F(i + 1, i + 1) = c;
            F(i + 1, i + 2) = -z * (s + w * h * c) - v * h * s;
        }
        P_ = F * P_ * F.transpose();
        for (int j = 0; j < MaxModes; ++j) {
            if (!active_[j]) continue;
            const int i = 3 * j;
            const T q = q_[j];
            P_(i, i)         += q * h * h * h / T(3);
            P_(i, i + 1)     += q * h * h / T(2);
            P_(i + 1, i)     += q * h * h / T(2);
            P_(i + 1, i + 1) += q * h;
            const T sw = cfg_.omega_rw * x_(i + 2);
            P_(i + 2, i + 2) += sw * sw * h;
        }
        P_(OFF, OFF) += cfg_.offset_rw * cfg_.offset_rw * h;
    }

    void correct_(T y) {
        Vec H = Vec::Zero();
        T yhat = x_(OFF);
        H(OFF) = T(1);
        for (int j = 0; j < MaxModes; ++j) {
            if (!active_[j]) continue;
            const int i = 3 * j;
            const T z = x_(i), w = x_(i + 2);
            yhat += -w * w * z;
            H(i) = -w * w;
            H(i + 2) = -T(2) * w * z;
        }
        const T R = cfg_.accel_noise_std * cfg_.accel_noise_std;
        const Vec PH = P_ * H;
        const T S = H.dot(PH) + R;
        const Vec K = PH / S;
        x_ += K * (y - yhat);
        P_ -= K * PH.transpose();
        P_ = T(0.5) * (P_ + P_.transpose());
    }

    void deactivate_(int j) {
        const int i = 3 * j;
        active_[j] = false;
        x_.template segment<3>(i).setZero();
        P_.template middleRows<3>(i).setZero();
        P_.template middleCols<3>(i).setZero();
        q_[j] = T(0);
    }

    void prune_() {
        const T wmin = T(2) * std::numbers::pi_v<T> * cfg_.f_min_hz;
        const T wmax = T(2) * std::numbers::pi_v<T> * cfg_.f_max_hz;
        for (int j = 0; j < MaxModes; ++j) {
            if (!active_[j]) continue;
            const T w = x_(3 * j + 2);
            if (!(w >= wmin && w <= wmax)) deactivate_(j);
        }
        for (int a = 0; a < MaxModes; ++a) {
            if (!active_[a]) continue;
            for (int b = a + 1; b < MaxModes; ++b) {
                if (!active_[b]) continue;
                const T wa = x_(3 * a + 2), wb = x_(3 * b + 2);
                if (std::abs(wa - wb) < cfg_.merge_tol * std::max(wa, wb)) merge_(a, b);
            }
        }
    }

    // Mode b folds into mode a: heave and heave rate add (two sinusoids of the
    // same frequency are one sinusoid); the frequency is amplitude weighted.
    void merge_(int a, int b) {
        const int ia = 3 * a, ib = 3 * b;
        const T Aa = modeAmplitude_(a), Ab = modeAmplitude_(b);
        const T wa = (Aa + Ab > T(0)) ? Aa / (Aa + Ab) : T(0.5);
        Mat Tm = Mat::Identity();
        Tm(ia, ib) = T(1);
        Tm(ia + 1, ib + 1) = T(1);
        Tm(ia + 2, ia + 2) = wa;
        Tm(ia + 2, ib + 2) = T(1) - wa;
        x_ = Tm * x_;
        P_ = Tm * P_ * Tm.transpose();
        q_[a] += q_[b];
        deactivate_(b);
    }

    void push_(T a) {
        ring_[static_cast<size_t>(ring_head_)] = a;
        ring_head_ = (ring_head_ + 1) % cfg_.fft_size;
        ring_count_ = std::min(ring_count_ + 1, cfg_.fft_size);
    }

    static void fft_(std::vector<std::complex<T>>& a) {
        const size_t n = a.size();
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(a[i], a[j]);
        }
        for (size_t len = 2; len <= n; len <<= 1) {
            const T ang = -T(2) * std::numbers::pi_v<T> / T(len);
            const std::complex<T> wl(std::cos(ang), std::sin(ang));
            for (size_t i = 0; i < n; i += len) {
                std::complex<T> w(1);
                for (size_t k = 0; k < len / 2; ++k) {
                    const auto u = a[i + k];
                    const auto v = a[i + k + len / 2] * w;
                    a[i + k] = u + v;
                    a[i + k + len / 2] = u - v;
                    w *= wl;
                }
            }
        }
    }

    void identify_() {
        const int n = cfg_.fft_size;
        std::vector<T> buf(static_cast<size_t>(n));
        T mean = T(0);
        for (int k = 0; k < n; ++k) {
            buf[static_cast<size_t>(k)] = ring_[static_cast<size_t>((ring_head_ + k) % n)];
            mean += buf[static_cast<size_t>(k)];
        }
        mean /= T(n);

        std::vector<std::complex<T>> X(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k) {
            const T w = T(0.5) - T(0.5) * std::cos(T(2) * std::numbers::pi_v<T> * T(k) / T(n));
            X[static_cast<size_t>(k)] = std::complex<T>(w * (buf[static_cast<size_t>(k)] - mean), T(0));
        }
        fft_(X);

        // Heave amplitude spectrum, lightly smoothed.
        const T df = cfg_.fs_ident_hz / T(n);
        const int k_lo = std::max(2, static_cast<int>(std::ceil(cfg_.f_min_hz / df)));
        const int k_hi = std::min(n / 2 - 2, static_cast<int>(std::floor(cfg_.f_max_hz / df)));
        if (k_hi <= k_lo + 2) return;
        std::vector<T> Z(static_cast<size_t>(n / 2), T(0));
        for (int k = 1; k < n / 2; ++k) {
            const T w = T(2) * std::numbers::pi_v<T> * df * T(k);
            Z[static_cast<size_t>(k)] = std::abs(X[static_cast<size_t>(k)]) / (w * w);
        }
        std::vector<T> Zs(Z);
        for (int k = 1; k < n / 2 - 1; ++k)
            Zs[static_cast<size_t>(k)] = (Z[static_cast<size_t>(k - 1)] + T(2) * Z[static_cast<size_t>(k)] + Z[static_cast<size_t>(k + 1)]) / T(4);

        T zmax = T(0);
        for (int k = k_lo; k <= k_hi; ++k) zmax = std::max(zmax, Zs[static_cast<size_t>(k)]);
        if (!(zmax > T(0))) return;

        std::vector<std::pair<T, T>> cand;  // (height, w)
        for (int k = k_lo; k <= k_hi; ++k) {
            const T c = Zs[static_cast<size_t>(k)];
            if (c > Zs[static_cast<size_t>(k - 1)] && c >= Zs[static_cast<size_t>(k + 1)] &&
                c >= cfg_.peak_rel_threshold * zmax)
            {
                // Parabolic refinement on log magnitude.
                const T l = std::log(Zs[static_cast<size_t>(k - 1)]);
                const T m = std::log(c);
                const T r = std::log(Zs[static_cast<size_t>(k + 1)]);
                const T den = l - T(2) * m + r;
                const T d = (den < T(0)) ? std::clamp(T(0.5) * (l - r) / den, T(-0.5), T(0.5)) : T(0);
                cand.emplace_back(c, T(2) * std::numbers::pi_v<T> * df * (T(k) + d));
            }
        }
        std::sort(cand.begin(), cand.end(), [](auto& a, auto& b) { return a.first > b.first; });
        if (static_cast<int>(cand.size()) > MaxModes) cand.resize(static_cast<size_t>(MaxModes));
        if (cand.empty()) return;

        T offset_fit = mean;
        std::vector<Peak> peaks = fit_(buf, mean, cand, offset_fit);
        ++identifications_;
        if (!started_) {
            // The EKF starts at the first identification, with the offset
            // seeded by the fitted constant of the buffer.
            x_(OFF) = offset_fit;
            P_(OFF, OFF) = cfg_.offset_init_std * cfg_.offset_init_std;
            started_ = true;
        }

        for (const auto& pk : peaks) {
            int match = -1;
            for (int j = 0; j < MaxModes; ++j) {
                if (!active_[j]) continue;
                if (std::abs(x_(3 * j + 2) - pk.w) < cfg_.match_tol * x_(3 * j + 2)) { match = j; break; }
            }
            if (match >= 0) {
                q_[match] = processNoise_(x_(3 * match + 2), pk.amp);
                continue;
            }
            int slot = -1;
            for (int j = 0; j < MaxModes; ++j) if (!active_[j]) { slot = j; break; }
            if (slot < 0) {
                int weakest = -1;
                T amin = std::numeric_limits<T>::infinity();
                for (int j = 0; j < MaxModes; ++j) {
                    const T a = modeAmplitude_(j);
                    if (a < amin) { amin = a; weakest = j; }
                }
                if (weakest >= 0 && pk.amp > cfg_.replace_ratio * amin) {
                    deactivate_(weakest);
                    slot = weakest;
                }
            }
            if (slot < 0) continue;
            seed_(slot, pk);
        }
    }

    T processNoise_(T w, T amp) const {
        return T(4) * cfg_.zeta_q * w * w * w * T(0.5) * amp * amp;
    }

    void seed_(int j, const Peak& pk) {
        const int i = 3 * j;
        active_[j] = true;
        x_(i) = pk.zc;
        x_(i + 1) = pk.zs;
        x_(i + 2) = pk.w;
        P_.template middleRows<3>(i).setZero();
        P_.template middleCols<3>(i).setZero();
        const T sa = cfg_.seed_rel_std * pk.amp;
        P_(i, i) = sa * sa;
        P_(i + 1, i + 1) = sa * sa * pk.w * pk.w;
        const T sw = cfg_.seed_omega_rel_std * pk.w;
        P_(i + 2, i + 2) = sw * sw;
        q_[j] = processNoise_(pk.w, pk.amp);
    }

    // Joint least squares of the buffer on cos/sin at each peak frequency plus
    // a constant, with time measured from the current instant.  Block averages
    // are stamped at their centres.  Returns heave and heave rate at the
    // current instant for each peak (z = -a/w^2).
    std::vector<Peak> fit_(const std::vector<T>& buf, T mean,
                           const std::vector<std::pair<T, T>>& cand,
                           T& offset) const
    {
        const int n = static_cast<int>(buf.size());
        const int m = static_cast<int>(cand.size());
        const int p = 2 * m + 1;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(n, p);
        Eigen::Matrix<T, Eigen::Dynamic, 1> b(n);
        for (int k = 0; k < n; ++k) {
            // Sample k is the block ending (n-1-k) blocks before now.
            const T tk = -(T(n - 1 - k) * block_dt_ + T(0.5) * (block_dt_ - raw_dt_));
            for (int c = 0; c < m; ++c) {
                const T w = cand[static_cast<size_t>(c)].second;
                A(k, 2 * c) = std::cos(w * tk);
                A(k, 2 * c + 1) = std::sin(w * tk);
            }
            A(k, p - 1) = T(1);
            b(k) = buf[static_cast<size_t>(k)] - mean;
        }
        const Eigen::Matrix<T, Eigen::Dynamic, 1> sol = A.colPivHouseholderQr().solve(b);
        offset = mean + sol(p - 1);
        std::vector<Peak> out;
        for (int c = 0; c < m; ++c) {
            const T w = cand[static_cast<size_t>(c)].second;
            const T ac = sol(2 * c), as = sol(2 * c + 1);
            // a(t) = ac cos(w t) + as sin(w t); z = -a / w^2.
            Peak pk;
            pk.w = w;
            pk.zc = -ac / (w * w);
            pk.zs = -as * w / (w * w);
            pk.amp = std::sqrt(ac * ac + as * as) / (w * w);
            out.push_back(pk);
        }
        return out;
    }

    Config cfg_;
    Vec x_ = Vec::Zero();
    Mat P_ = Mat::Zero();
    std::array<bool, MaxModes> active_{};
    std::array<T, MaxModes> q_{};
    T t_ = T(0);
    bool started_ = false;
    T block_sum_ = T(0);
    int block_n_ = 0;
    T block_dt_ = T(0.25);
    T raw_dt_ = T(0.005);
    std::vector<T> ring_;
    int ring_head_ = 0;
    int ring_count_ = 0;
    T since_ident_ = T(0);
    int identifications_ = 0;
};

} // namespace heave_baselines
