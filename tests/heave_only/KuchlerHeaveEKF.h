#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Heave-only filter: the harmonic-mode heave observer of

    S. Kuchler, J. K. Eberharter, K. Langer, K. Schneider, O. Sawodny,
    "Heave motion estimation of a vessel using acceleration measurements",
    18th IFAC World Congress, pp. 14742-14747, 2011.

  Equation numbers below are the paper's.

  Heave is a sum of N_m undamped modes (Eqs. 1, 3), each with state
  x_j = [z_j, z_j', w_j] and dynamics (Eq. 5)

      z_j'' = -w_j^2 z_j,   w_j random walk,   y_j = -w_j^2 z_j,

  plus a random-walk offset state for gravity and sensor bias (Eq. 8).  The
  measurement is y = sum_j y_j + y_off (Eqs. 7, 9).  The input here is the
  levelled vertical acceleration with gravity removed, so the paper's
  x_off,0 = -g becomes x_off,0 = 0.  Each mode is discretised exactly for the
  current w_j (Eq. 6): the frequency first, then the linear oscillator.

  Identification (Section 2.2, Fig. 2): an FFT of the buffered acceleration
  gives the heave amplitude spectrum |A_acc(w)| / w^2 (Eq. 2), whose peaks set
  the number of modes and their frequencies.  The EKF starts at the first
  identification.  The identification repeats at fixed intervals; a peak that
  matches no existing mode extends the model and only the new mode's states
  are initialised, all others keep their estimates (Section 2.3).  A mode
  whose frequency approaches zero is removed and two modes whose frequencies
  approach each other are merged (observability conditions, Eq. 11).
  Process noise penalises position and velocity of higher-frequency modes
  more; the frequency and offset states are slow; R is the datasheet sensor
  noise (Section 2.3).

  Choices the paper leaves open, made here:
    - new modes are seeded with amplitude and phase from a joint
      least-squares fit of the buffer at the peak frequencies (the paper
      reads them off the FFT bins, Eq. 2; the fit avoids the bin-phase error
      of an off-bin peak);
    - "approaches zero" is w_j below f_min_hz; modes above f_max_hz are also
      removed;
    - at most MaxModes modes; a new peak replace_ratio times stronger than
      the weakest mode replaces it;
    - velocity process noise q_j = 4 zeta_q w_j^3 (A_j^2/2), the intensity
      that holds a damped oscillator of damping zeta_q at the identified
      amplitude A_j, so one setting serves every sea state.
  The paper uses the raw body-z accelerometer and neglects roll and pitch
  (cos(phi) cos(theta) ~ 1); the replay feeds levelled acceleration, which
  only removes that approximation's error.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include "AccelSpectrum.h"

namespace heave_only {

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
        T offset_init_std = T(0.1);      // m/s^2
        // Seed the offset from the fitted buffer constant instead of the
        // paper's known-gravity value (0 after gravity removal).
        bool offset_from_fit = false;
        T seed_rel_std = T(0.5);         // initial amplitude uncertainty
        T seed_omega_rel_std = T(0.05);
    };

    KuchlerHeaveEKF() : KuchlerHeaveEKF(Config{}) {}
    explicit KuchlerHeaveEKF(const Config& cfg)
        : cfg_(cfg), spectrum_(cfg.fs_ident_hz, cfg.fft_size) { reset(); }

    void reset() {
        x_.setZero();
        P_.setZero();
        active_.fill(false);
        q_.fill(T(0));
        since_ident_ = T(0);
        started_ = false;
        identifications_ = 0;
        spectrum_ = AccelSpectrum<T>(cfg_.fs_ident_hz, cfg_.fft_size);
    }

    // a_up: vertical specific force minus gravity, m/s^2, z up.
    T update(T a_up, T dt) {
        if (!(dt > T(0)) || !std::isfinite(a_up)) return heave();

        if (started_) {
            predict_(dt);
            correct_(a_up);
            prune_();
        }

        const bool pushed = spectrum_.add(a_up, dt);
        since_ident_ += dt;
        if (pushed && spectrum_.full() && since_ident_ >= cfg_.ident_period_s) {
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
    struct Seed { T w; T amp; T z; T v; };

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

    void identify_() {
        const auto peaks = spectrum_.peaks(cfg_.f_min_hz, cfg_.f_max_hz,
                                           cfg_.peak_rel_threshold, MaxModes);
        if (peaks.empty()) return;
        T offset_fit = T(0);
        const std::vector<Seed> seeds = fit_(peaks, offset_fit);
        ++identifications_;
        if (!started_) {
            x_(OFF) = cfg_.offset_from_fit ? offset_fit : T(0);
            P_(OFF, OFF) = cfg_.offset_init_std * cfg_.offset_init_std;
            started_ = true;
        }

        for (const auto& sd : seeds) {
            int match = -1;
            for (int j = 0; j < MaxModes; ++j) {
                if (!active_[j]) continue;
                if (std::abs(x_(3 * j + 2) - sd.w) < cfg_.match_tol * x_(3 * j + 2)) { match = j; break; }
            }
            if (match >= 0) {
                q_[match] = processNoise_(x_(3 * match + 2), sd.amp);
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
                if (weakest >= 0 && sd.amp > cfg_.replace_ratio * amin) {
                    deactivate_(weakest);
                    slot = weakest;
                }
            }
            if (slot >= 0) seed_(slot, sd);
        }
    }

    T processNoise_(T w, T amp) const {
        return T(4) * cfg_.zeta_q * w * w * w * T(0.5) * amp * amp;
    }

    void seed_(int j, const Seed& sd) {
        const int i = 3 * j;
        active_[j] = true;
        x_(i) = sd.z;
        x_(i + 1) = sd.v;
        x_(i + 2) = sd.w;
        P_.template middleRows<3>(i).setZero();
        P_.template middleCols<3>(i).setZero();
        const T sa = cfg_.seed_rel_std * sd.amp;
        P_(i, i) = sa * sa;
        P_(i + 1, i + 1) = sa * sa * sd.w * sd.w;
        const T sw = cfg_.seed_omega_rel_std * sd.w;
        P_(i + 2, i + 2) = sw * sw;
        q_[j] = processNoise_(sd.w, sd.amp);
    }

    // Joint least squares of the buffer on cos/sin at each peak frequency plus
    // a constant, with time measured from the current instant.  Block averages
    // are stamped at their centres.  Returns heave and heave rate at the
    // current instant for each peak (z = -a/w^2).
    std::vector<Seed> fit_(const std::vector<typename AccelSpectrum<T>::Peak>& peaks,
                           T& offset) const
    {
        const std::vector<T> buf = spectrum_.samples();
        const int n = static_cast<int>(buf.size());
        const int m = static_cast<int>(peaks.size());
        const int p = 2 * m + 1;
        const T bdt = spectrum_.blockDt();
        const T rdt = spectrum_.rawDt();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(n, p);
        Eigen::Matrix<T, Eigen::Dynamic, 1> b(n);
        for (int k = 0; k < n; ++k) {
            const T tk = -(T(n - 1 - k) * bdt + T(0.5) * (bdt - rdt));
            for (int c = 0; c < m; ++c) {
                const T w = peaks[static_cast<size_t>(c)].w;
                A(k, 2 * c) = std::cos(w * tk);
                A(k, 2 * c + 1) = std::sin(w * tk);
            }
            A(k, p - 1) = T(1);
            b(k) = buf[static_cast<size_t>(k)];
        }
        const Eigen::Matrix<T, Eigen::Dynamic, 1> sol = A.colPivHouseholderQr().solve(b);
        offset = sol(p - 1);
        std::vector<Seed> out;
        for (int c = 0; c < m; ++c) {
            const T w = peaks[static_cast<size_t>(c)].w;
            const T ac = sol(2 * c), as = sol(2 * c + 1);
            // a(t) = ac cos(w t) + as sin(w t); z = -a / w^2.
            out.push_back({w, std::sqrt(ac * ac + as * as) / (w * w),
                           -ac / (w * w), -as / w});
        }
        return out;
    }

    Config cfg_;
    AccelSpectrum<T> spectrum_;
    Vec x_ = Vec::Zero();
    Mat P_ = Mat::Zero();
    std::array<bool, MaxModes> active_{};
    std::array<T, MaxModes> q_{};
    T since_ident_ = T(0);
    bool started_ = false;
    int identifications_ = 0;
};

} // namespace heave_only
