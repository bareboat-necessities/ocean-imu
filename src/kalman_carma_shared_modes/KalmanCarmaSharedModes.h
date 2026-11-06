#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    KalmanCarmaSharedModes -- EXPERIMENTAL standalone QMEKF marine-motion
    filter with shared directional stochastic modes (Option C of the
    coupled-prior study; docs/kalman-carma-shared-modes.md,
    doc/kalman_carma_shared_modes/).  Not a deployed estimator and not a
    default: SeaStateFusionFilter_OU_III stays the shipping filter and the
    independent-axis CARMA3 filter (kalman_carma/) is unchanged.

    Composition (nothing re-implemented):
      * the QMEKF core Kalman3D_Wave_CARMA3 (attitude error, gyroscope bias,
        accelerometer bias, Joseph updates, multiplicative attitude reset)
        with the motion block SharedModesMotion: per axis one complexified
        CARMA(4,1) mode (8 real states), the axes sharing its dynamics and
        driven by circular complex noise with a Hermitian excitation, so the
        x-z / y-z quadrature of linear waves is representable (33 error
        states: 3 attitude, 3 gyro bias, 24 motion, 3 accelerometer bias);
      * the CARMA3 wrapper schedule SeaStateFusion_CARMA3 (proxy startup,
        magnetic acquisition, tuner, 0.1 s model commits, tilt watchdog);
      * the deployable coupling estimate (band cross-spectrum of the filter's
        own world-frame acceleration at the operating frequency, 240 s);
      * the coupling confidence (CarmaCouplingConfidence.h), which scales the
        applied coupling c in [0, 1] and, when c falls, shrinks the
        unexplained horizontal-vertical motion cross-covariance by the same
        ratio (Kalman3D_Wave_CARMA3::shrink_motion_vertical_coupling).  At
        c = 0 the prior is the independent-axis CARMA(4,1) law on every axis.

    Memory: begin() makes two heap allocations, the wrapper (schedule,
    front end, tuner) and, inside it, the 33-state QMEKF (covariance 33 x 33
    float, discretization and coupling state); update() allocates nothing.
    Outputs before begin() read the un-started wrapper (not Live).

    Ground truth, oracle coupling and RAO knowledge never enter this class.
    setExternalCoupling() exists for research (a coupling prior from another
    source, e.g. a wave buoy, or a deliberately wrong one in the robustness
    tests); the confidence checks it against the IMU like the estimated one.
*/

#include <algorithm>
#include <cmath>
#include <memory>

#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma_shared_modes/CarmaSharedModesMotion.h"
#include "kalman_carma_shared_modes/CarmaCouplingConfidence.h"

namespace ocean_imu::carma_shared_modes {

enum class CouplingMode : unsigned char {
    Estimated,   // deployable: the filter's own band cross-spectrum (default)
    External,    // research: a prior supplied with setExternalCoupling()
    Off          // independent axes (CARMA(4,1) law on every axis)
};

class KalmanCarmaSharedModes {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Motion = SharedModesMotion<float>;
    using Mekf = kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false, true, true, Motion>;
    using Fusion = carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, Mekf>;
    using Hermitian = kalman::carma::Hermitian3<float>;
    using CouplingConfig = kalman::carma::CarmaCouplingConfig;
    static constexpr int kStates = Mekf::NX;     // 33

    struct Config {
        carma::SeaStateFusionConfig_CARMA3 fusion{};   // sensors, magnetometer, startup
        CouplingMode mode = CouplingMode::Estimated;
        CouplingConfidenceConfig confidence{};
        float coupling_gain = 1.0f;     // upper bound of the applied coupling
        float horizon_s = 240.0f;       // coupling estimate (screened, seeds 1001-1020)
        float band_frac = 0.15f;        // envelope band / f_op (screened)
        float warmup_s = 120.0f;        // coupling estimate warm-up
        bool tuner = true;
        bool clamp = true;
    };

    struct SeaState {
        float frequency_hz = NAN;       // the tuner's wave frequency
        float operating_hz = NAN;       // CARMA operating frequency kappa f_z,base
        float kappa = NAN;              // vertical frequency scale
        Eigen::Vector3f accel_std = Eigen::Vector3f::Constant(NAN);   // committed prior (x, y, z)
    };

    struct CouplingStatus {
        CouplingMode mode = CouplingMode::Estimated;
        CouplingConfidenceStatus confidence{};
        float applied_coherence_xz = 0.0f, applied_coherence_yz = 0.0f;   // |Sigma_iz|^2/(Sigma_ii Sigma_zz)
        float applied_phase_xz_deg = 0.0f, applied_phase_yz_deg = 0.0f;
        unsigned long shrink_events = 0;
    };

    struct Health {
        bool live = false;
        float acc_nis = NAN;            // last accelerometer NIS (3 dof)
        float nis_ratio = 1.0f;         // short / long mean NIS
        bool covariance_finite = false;
        float covariance_asymmetry = NAN;   // max |P - P^T|
        float covariance_min_diag = NAN;
    };

    // An un-started wrapper exists from construction (as a plain wrapper
    // member would), so the outputs are safe to read before begin().
    KalmanCarmaSharedModes() : fusion_(std::make_unique<Fusion>()) {}
    explicit KalmanCarmaSharedModes(const Config& c) { begin(c); }

    void begin(const Config& c) {
        cfg_ = c;
        fusion_ = std::make_unique<Fusion>();
        fusion_->begin(cfg_.fusion);
        auto& f = fusion_->raw();
        f.enableTuner(cfg_.tuner);
        f.enableClamp(cfg_.clamp);
        auto& m = f.mekf().motion();
        CouplingConfig cc;
        cc.source = cfg_.mode == CouplingMode::Off ? CouplingConfig::Source::Off
                  : cfg_.mode == CouplingMode::Estimated ? CouplingConfig::Source::Estimated
                                                         : CouplingConfig::Source::External;
        cc.frame = CouplingConfig::Frame::World;
        cc.gain = cfg_.coupling_gain;
        cc.warmup_s = cfg_.warmup_s;
        if (cfg_.mode == CouplingMode::External) cc.shrink = false;
        m.coupling().set_config(cc);
        typename Motion::Estimator::Config ec;
        ec.horizon_s = cfg_.horizon_s;
        ec.band_frac = cfg_.band_frac;
        m.coupling().estimator().configure(ec);
        ec.horizon_s = cfg_.confidence.monitor_horizon_s;
        m.monitor().configure(ec);
        CouplingConfidenceConfig k = cfg_.confidence;
        if (cfg_.mode == CouplingMode::Off) k.enabled = false;
        conf_.configure(k);
        conf_.set_band_fraction(cfg_.band_frac);
        applied_c_ = shrunk_c_ = (k.enabled ? 0.0f : 1.0f);
        prior_ = BandRegression{};
        prior_neff_ = 0.0f;
        shrink_events_ = 0;
        commit_clock_ = 0.0f;
        apply_gain_(cfg_.mode == CouplingMode::Off ? 0.0f : applied_c_);
    }

    // Restart from scratch with the configuration of the last begin().
    void reset() { begin(cfg_); }

    // One IMU sample: rad/s and m/s^2 in the body NED frame (as the CARMA3
    // wrapper), temperature in deg C.
    void update(float dt, const Eigen::Vector3f& gyr, const Eigen::Vector3f& acc, float temperature_c = 35.0f) {
        fusion_->update(dt, gyr, acc, temperature_c);
        step_confidence_(dt, gyr);
    }
    // Magnetometer sample (uT, body NED), at its own rate.
    void updateMag(const Eigen::Vector3f& mag_body_ned) { fusion_->updateMag(mag_body_ned); }

    // Research: a coupling prior from outside (CouplingMode::External).
    void setExternalCoupling(const Hermitian& S, float n_eff = 0.0f) {
        fusion_->raw().mekf().motion().coupling().set_external(S, n_eff);
    }
    // Research / tests: hold the confidence (negative: release).
    void forceCouplingConfidence(float c) { conf_.force(c); }

    // ---- outputs (NED, filter's magnetic world frame) ----
    bool isLive() const { return fusion_->isLive(); }
    Eigen::Quaternionf attitude() const { return fusion_->raw().mekf().quaternion_boat(); }   // body -> world
    Eigen::Vector3f position() const { return fusion_->raw().mekf().get_position(); }
    Eigen::Vector3f velocity() const { return fusion_->raw().mekf().get_velocity(); }
    Eigen::Vector3f waveAcceleration() const { return fusion_->raw().mekf().get_world_accel(); }
    Eigen::Vector3f accelBias() const { return fusion_->raw().mekf().get_acc_bias(); }
    Eigen::Vector3f gyroBias() const { return fusion_->raw().mekf().gyroscope_bias(); }
    // Covariance of [p; v] (6 x 6) through the motion readout rows.
    Eigen::Matrix<float,6,6> positionVelocityCovariance() const {
        const auto& m = fusion_->raw().mekf();
        return m.pv_covariance6();
    }

    SeaState seaState() const {
        SeaState s;
        const auto& f = fusion_->raw();
        const auto& m = f.mekf().motion();
        s.frequency_hz = f.getFreqHz();
        s.operating_hz = m.operating_frequency_hz();
        s.kappa = f.getVerticalKappa();
        s.accel_std = m.accel_std();
        return s;
    }

    // Direction of travel of the dominant waves from the coupling prior
    // (degrees clockwise from the filter's magnetic north) and its
    // large-sample uncertainty; invalid until the prior exists.
    WaveDirectionEstimate waveDirection() const { return wave_direction(prior_, prior_neff_); }

    CouplingStatus coupling() const {
        CouplingStatus c;
        c.mode = cfg_.mode;
        c.confidence = conf_.status();
        const Hermitian& E = fusion_->raw().mekf().motion().excitation();
        for (int i = 0; i < 2; ++i) {
            const std::complex<float> s(E.re(i,2), E.im(i,2));
            const float den = E.re(i,i) * E.re(2,2);
            const float coh = den > 0.0f ? std::norm(s) / den : 0.0f;
            const float ph = std::abs(s) > 0.0f ? std::arg(s) * 57.2957795f : 0.0f;
            (i == 0 ? c.applied_coherence_xz : c.applied_coherence_yz) = coh;
            (i == 0 ? c.applied_phase_xz_deg : c.applied_phase_yz_deg) = ph;
        }
        c.shrink_events = shrink_events_;
        return c;
    }

    Health health() const {
        Health h;
        const auto& m = fusion_->raw().mekf();
        h.live = fusion_->isLive();
        h.acc_nis = m.lastAccDiag().nis;
        h.nis_ratio = conf_.status().nis_ratio;
        const auto& P = m.covariance_full();
        h.covariance_finite = P.allFinite();
        h.covariance_asymmetry = (P - P.transpose()).cwiseAbs().maxCoeff();
        h.covariance_min_diag = P.diagonal().minCoeff();
        return h;
    }

    float couplingConfidence() const noexcept { return conf_.confidence(); }
    const Config& config() const noexcept { return cfg_; }

    // Escape hatch for the simulators and tests (not a stable API).
    Fusion& fusion() noexcept { return *fusion_; }
    const Fusion& fusion() const noexcept { return *fusion_; }

private:
    void apply_gain_(float c) {
        auto& m = fusion_->raw().mekf().motion();
        CouplingConfig cc = m.coupling().config();   // only the gain changes
        cc.gain = cfg_.coupling_gain * std::min(1.0f, std::max(0.0f, c));
        m.coupling().set_config(cc);
    }

    void step_confidence_(float dt, const Eigen::Vector3f& gyr) {
        if (cfg_.mode == CouplingMode::Off) return;
        auto& mekf = fusion_->raw().mekf();
        auto& m = mekf.motion();
        // The candidate prior (what the coupling would apply at c = 1): the
        // coupling reports nothing while its gain is zero, so read it with
        // the gain lifted (configuration only; no state changes).
        Hermitian Sp;
        float np = 0.0f;
        BandRegression prior;
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
        if (have) prior = band_regression(Sp, cfg_.mode == CouplingMode::External ? 0.0f : np);
        prior_ = prior;
        prior_neff_ = cfg_.mode == CouplingMode::External ? 1.0e6f : np;
        const auto& mon = m.monitor();
        const BandRegression monitor = band_regression(mon.cross(), mon.effective_samples());
        const bool full = mon.elapsed_s() >= cfg_.confidence.monitor_horizon_s;
        const Eigen::Vector3f w_world = mekf.quaternion_boat() * (gyr - mekf.gyroscope_bias());
        conf_.update(dt, prior, monitor, full, w_world.z(), m.operating_frequency_hz(), mekf.lastAccDiag().nis);
        // Commit at the model cadence (0.1 s): the excitation is rebuilt at
        // the wrapper's next commit; a falling confidence also shrinks the
        // unexplained horizontal-vertical cross-covariance by c_new / c_old.
        commit_clock_ += dt;
        if (commit_clock_ < 0.1f) return;
        commit_clock_ = 0.0f;
        const float c = conf_.confidence();
        if (c < shrunk_c_ && shrunk_c_ > 0.0f && fusion_->isLive()) {
            const float ratio = c / shrunk_c_;
            if (ratio < 0.999f && mekf.shrink_motion_vertical_coupling(ratio)) ++shrink_events_;
        }
        shrunk_c_ = c;
        applied_c_ = c;
        apply_gain_(c);
    }

    Config cfg_{};
    std::unique_ptr<Fusion> fusion_;   // created by begin()
    CarmaCouplingConfidence conf_{};
    BandRegression prior_{};
    float prior_neff_ = 0.0f;
    float applied_c_ = 0.0f, shrunk_c_ = 0.0f;
    float commit_clock_ = 0.0f;
    unsigned long shrink_events_ = 0;
};

}  // namespace ocean_imu::carma_shared_modes
