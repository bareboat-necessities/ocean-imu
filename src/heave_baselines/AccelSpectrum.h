#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Shared identification front end of the literature heave baselines: the
  vertical acceleration is block-averaged into a ring buffer, and the buffer's
  Hann-windowed FFT gives the heave amplitude spectrum |A_acc(w)| / w^2
  (Kuchler et al. 2011, Eq. 2) and its peaks.
*/

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <utility>
#include <vector>

namespace heave_baselines {

template <typename T>
class AccelSpectrum {
public:
    struct Peak {
        T w;       // rad/s, parabolic refinement on the log spectrum
        T height;  // smoothed heave amplitude spectrum at the peak
    };

    AccelSpectrum(T fs_hz, int size) : fs_hz_(fs_hz), n_(size) {
        ring_.assign(static_cast<size_t>(n_), T(0));
    }

    // Returns true when a block average was pushed.
    bool add(T a, T dt) {
        block_sum_ += a;
        ++block_n_;
        const int block_len = std::max(1, static_cast<int>(std::lround(T(1) / (fs_hz_ * dt))));
        if (block_n_ < block_len) return false;
        ring_[static_cast<size_t>(head_)] = block_sum_ / T(block_n_);
        head_ = (head_ + 1) % n_;
        count_ = std::min(count_ + 1, n_);
        block_sum_ = T(0);
        block_n_ = 0;
        block_dt_ = T(block_len) * dt;
        raw_dt_ = dt;
        return true;
    }

    bool full() const { return count_ >= n_; }
    int size() const { return n_; }
    T blockDt() const { return block_dt_; }
    T rawDt() const { return raw_dt_; }

    // Oldest-first copy of the buffer.
    std::vector<T> samples() const {
        std::vector<T> buf(static_cast<size_t>(n_));
        for (int k = 0; k < n_; ++k)
            buf[static_cast<size_t>(k)] = ring_[static_cast<size_t>((head_ + k) % n_)];
        return buf;
    }

    // Peaks of the heave amplitude spectrum in [f_min, f_max], strongest
    // first, at least rel_threshold of the largest, at most max_peaks.
    std::vector<Peak> peaks(T f_min_hz, T f_max_hz, T rel_threshold, int max_peaks) const {
        const std::vector<T> buf = samples();
        T mean = T(0);
        for (T v : buf) mean += v;
        mean /= T(n_);

        std::vector<std::complex<T>> X(static_cast<size_t>(n_));
        for (int k = 0; k < n_; ++k) {
            const T w = T(0.5) - T(0.5) * std::cos(T(2) * std::numbers::pi_v<T> * T(k) / T(n_));
            X[static_cast<size_t>(k)] = std::complex<T>(w * (buf[static_cast<size_t>(k)] - mean), T(0));
        }
        fft_(X);

        const T df = fs_hz_ / T(n_);
        const int half = n_ / 2;
        const int k_lo = std::max(2, static_cast<int>(std::ceil(f_min_hz / df)));
        const int k_hi = std::min(half - 2, static_cast<int>(std::floor(f_max_hz / df)));
        std::vector<Peak> out;
        if (k_hi <= k_lo + 2) return out;

        std::vector<T> Z(static_cast<size_t>(half), T(0));
        for (int k = 1; k < half; ++k) {
            const T w = T(2) * std::numbers::pi_v<T> * df * T(k);
            Z[static_cast<size_t>(k)] = std::abs(X[static_cast<size_t>(k)]) / (w * w);
        }
        std::vector<T> Zs(Z);
        for (int k = 1; k < half - 1; ++k) {
            Zs[static_cast<size_t>(k)] = (Z[static_cast<size_t>(k - 1)] + T(2) * Z[static_cast<size_t>(k)] +
                                          Z[static_cast<size_t>(k + 1)]) / T(4);
        }

        T zmax = T(0);
        for (int k = k_lo; k <= k_hi; ++k) zmax = std::max(zmax, Zs[static_cast<size_t>(k)]);
        if (!(zmax > T(0))) return out;

        for (int k = k_lo; k <= k_hi; ++k) {
            const T c = Zs[static_cast<size_t>(k)];
            if (c > Zs[static_cast<size_t>(k - 1)] && c >= Zs[static_cast<size_t>(k + 1)] &&
                c >= rel_threshold * zmax)
            {
                const T l = std::log(Zs[static_cast<size_t>(k - 1)]);
                const T m = std::log(c);
                const T r = std::log(Zs[static_cast<size_t>(k + 1)]);
                const T den = l - T(2) * m + r;
                const T d = (den < T(0)) ? std::clamp(T(0.5) * (l - r) / den, T(-0.5), T(0.5)) : T(0);
                out.push_back({T(2) * std::numbers::pi_v<T> * df * (T(k) + d), c});
            }
        }
        std::sort(out.begin(), out.end(), [](const Peak& a, const Peak& b) { return a.height > b.height; });
        if (static_cast<int>(out.size()) > max_peaks) out.resize(static_cast<size_t>(max_peaks));
        return out;
    }

private:
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

    T fs_hz_;
    int n_;
    std::vector<T> ring_;
    int head_ = 0;
    int count_ = 0;
    T block_sum_ = T(0);
    int block_n_ = 0;
    T block_dt_ = T(0.25);
    T raw_dt_ = T(0.005);
};

} // namespace heave_baselines
