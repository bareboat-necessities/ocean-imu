#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Geometry-independent shared-mode configuration, confidence and fallback.
// Core supplies the estimator interface, motion() and the coordinate-aware
// shrink_motion_vertical_coupling(). No stochastic or confidence code is copied.
#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include "CarmaSharedModesMotion.h"
#include "CarmaCouplingConfidence.h"

namespace ocean_imu::carma_shared_modes {

enum class CouplingMode : unsigned char {
    Estimated,   // deployable: the filter's own band cross-spectrum (default)
    External,    // research: a prior supplied with setExternalCoupling()
    Off          // independent axes (CARMA(4,1) law on every axis)
};

template<typename T = float>
struct SharedModesConfig {
    CouplingMode mode = CouplingMode::Estimated;
    CouplingConfidenceConfig<T> confidence{};
    T coupling_gain = T(1);       // upper bound of the applied coupling
    T horizon_s = T(240);         // coupling estimate averaging horizon
    T band_frac = T(0.15);        // envelope band / f_op
    T warmup_s = T(120);          // coupling estimate warm-up
};

template<typename T, class Core>
class SharedModesEstimator : public Core {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Base = Core;
    using Vector3 = typename Base::Vector3;
    using Motion = SharedModesMotion<T>;
    using Hermitian = kalman::carma::Hermitian3<T>;
    using CouplingConfig = kalman::carma::CarmaCouplingConfig;
    using Config = SharedModesConfig<T>;
    using Confidence = CarmaCouplingConfidence<T>;
    static constexpr int kStates = Base::NX;     // 27

    struct CouplingStatus {
        CouplingMode mode = CouplingMode::Estimated;
        CouplingConfidenceStatus<T> confidence{};
        T applied_coherence_xz = T(0), applied_coherence_yz = T(0);   // |Sigma_iz|^2/(Sigma_ii Sigma_zz)
        T applied_phase_xz_deg = T(0), applied_phase_yz_deg = T(0);
        unsigned long shrink_events = 0;
    };

    struct Health {
        T acc_nis = std::numeric_limits<T>::quiet_NaN();   // last accelerometer NIS (3 dof)
        T nis_ratio = T(1);                                 // short / long mean NIS
        bool covariance_finite = false;
        T covariance_asymmetry = std::numeric_limits<T>::quiet_NaN();   // max |P - P^T|
        T covariance_min_diag = std::numeric_limits<T>::quiet_NaN();
    };

    // The constructors of the CARMA3 core (noise densities, initial
    // covariances, ...); the shared-modes configuration defaults to the
    // deployable one until configureSharedModes() is called.
    template<class... Args>
    explicit SharedModesEstimator(Args&&... args) : Base(std::forward<Args>(args)...) {
        configureSharedModes(Config{});
    }

    void configureSharedModes(const Config& c) {
        cfg_ = c;
        auto& m = this->motion();
        CouplingConfig cc = m.coupling().config();
        cc.source = cfg_.mode == CouplingMode::Off ? CouplingConfig::Source::Off
                  : cfg_.mode == CouplingMode::Estimated ? CouplingConfig::Source::Estimated
                                                         : CouplingConfig::Source::External;
        cc.frame = CouplingConfig::Frame::World;
        cc.gain = static_cast<float>(cfg_.coupling_gain);
        cc.warmup_s = static_cast<float>(cfg_.warmup_s);
        cc.shrink = cfg_.mode != CouplingMode::External;
        m.coupling().set_config(cc);
        typename Motion::Estimator::Config ec;
        ec.horizon_s = cfg_.horizon_s;
        ec.band_frac = cfg_.band_frac;
        m.coupling().estimator().configure(ec);
        ec.horizon_s = cfg_.confidence.monitor_horizon_s;
        m.monitor().configure(ec);
        CouplingConfidenceConfig<T> k = cfg_.confidence;
        if (cfg_.mode == CouplingMode::Off) k.enabled = false;
        conf_.configure(k);
        conf_.set_band_fraction(cfg_.band_frac);
        shrunk_c_ = k.enabled ? T(0) : T(1);
        prior_ = BandRegression<T>{};
        prior_neff_ = T(0);
        shrink_events_ = 0;
        commit_clock_ = T(0);
        last_dt_ = T(0);
        apply_gain_(cfg_.mode == CouplingMode::Off ? T(0) : shrunk_c_);
    }
    const Config& sharedModesConfig() const noexcept { return cfg_; }

    // ---- per-sample interface (the CARMA3 core's, extended) ----
    void time_update(Vector3 const& gyr, T Ts) {
        Base::time_update(gyr, Ts);
        last_gyr_ = gyr;
        last_dt_ = Ts;
    }
    void measurement_update_acc_only(Vector3 const& acc, T tempC = Base::tempC_ref) {
        Base::measurement_update_acc_only(acc, tempC);
        step_confidence_();
    }

    // Research: a coupling prior from outside (CouplingMode::External).
    void setExternalCoupling(const Hermitian& S, T n_eff = T(0)) {
        this->motion().coupling().set_external(S, n_eff);
    }
    // Research / tests: hold the confidence (negative: release).
    void forceCouplingConfidence(T c) { conf_.force(c); }

    // ---- diagnostics ----
    T couplingConfidence() const noexcept { return conf_.confidence(); }

    // Direction of travel of the dominant waves from the coupling prior
    // (degrees clockwise from the filter's north) and its large-sample
    // uncertainty; invalid until the prior exists.
    WaveDirectionEstimate<T> waveDirection() const { return wave_direction<T>(prior_, prior_neff_); }

    CouplingStatus coupling() const {
        CouplingStatus c;
        c.mode = cfg_.mode;
        c.confidence = conf_.status();
        const Hermitian& E = this->motion().excitation();
        for (int i = 0; i < 2; ++i) {
            const std::complex<T> s(E.re(i,2), E.im(i,2));
            const T den = E.re(i,i) * E.re(2,2);
            const T coh = den > T(0) ? std::norm(s) / den : T(0);
            const T ph = std::abs(s) > T(0) ? std::arg(s) * T(57.29577951308232) : T(0);
            (i == 0 ? c.applied_coherence_xz : c.applied_coherence_yz) = coh;
            (i == 0 ? c.applied_phase_xz_deg : c.applied_phase_yz_deg) = ph;
        }
        c.shrink_events = shrink_events_;
        return c;
    }

    Health health() const {
        Health h;
        h.acc_nis = this->lastAccDiag().nis;
        h.nis_ratio = conf_.status().nis_ratio;
        const auto& P = this->covariance_full();
        h.covariance_finite = P.allFinite();
        h.covariance_asymmetry = (P - P.transpose()).cwiseAbs().maxCoeff();
        h.covariance_min_diag = P.diagonal().minCoeff();
        return h;
    }

private:
    void apply_gain_(T c) {
        auto& m = this->motion();
        CouplingConfig cc = m.coupling().config();   // only the gain changes
        cc.gain = static_cast<float>(cfg_.coupling_gain * std::min(T(1), std::max(T(0), c)));
        m.coupling().set_config(cc);
    }

    void step_confidence_() {
        if (cfg_.mode == CouplingMode::Off || !(last_dt_ > T(0))) return;
        auto& m = this->motion();
        // The candidate prior (what the coupling would apply at c = 1): the
        // coupling reports nothing while its gain is zero, so read it with
        // the gain lifted (configuration only; no state changes).
        Hermitian Sp;
        T np = T(0);
        BandRegression<T> prior;
        auto& cp = m.coupling();
        const CouplingConfig cur = cp.config();
        bool have = false;
        if (cur.gain == 0.0f) {
            CouplingConfig g1 = cur;
            g1.gain = 1.0f;
            cp.set_config(g1);
            have = cp.band_cross(Sp, np);
            cp.set_config(cur);
        } else {
            have = cp.band_cross(Sp, np);
        }
        const bool external = cfg_.mode == CouplingMode::External;
        if (have) prior = band_regression(Sp, external ? T(0) : np);
        prior_ = prior;
        prior_neff_ = external ? T(1e6) : np;
        const auto& mon = m.monitor();
        const BandRegression<T> monitor = band_regression(mon.cross(), mon.effective_samples());
        const bool full = mon.elapsed_s() >= cfg_.confidence.monitor_horizon_s;
        const Vector3 w_world = this->quaternion_boat() * (last_gyr_ - this->gyroscope_bias());
        conf_.update(last_dt_, prior, monitor, full, w_world.z(), m.operating_frequency_hz(), this->lastAccDiag().nis);
        // Commit at the model cadence (0.1 s): the excitation is rebuilt at
        // the schedule's next model commit; a falling confidence also shrinks
        // the unexplained horizontal-vertical cross-covariance by c_new / c_old.
        commit_clock_ += last_dt_;
        if (commit_clock_ < T(0.1)) return;
        commit_clock_ = T(0);
        const T c = conf_.confidence();
        if (c < shrunk_c_ && shrunk_c_ > T(0)) {
            const T ratio = c / shrunk_c_;
            if (ratio < T(0.999) && this->shrink_motion_vertical_coupling(ratio)) ++shrink_events_;
        }
        shrunk_c_ = c;
        apply_gain_(c);
    }

    Config cfg_{};
    Confidence conf_{};
    BandRegression<T> prior_{};
    T prior_neff_ = T(0);
    Vector3 last_gyr_ = Vector3::Zero();
    T last_dt_ = T(0);
    T shrunk_c_ = T(0);
    T commit_clock_ = T(0);
    unsigned long shrink_events_ = 0;
};

}  // namespace ocean_imu::carma_shared_modes
