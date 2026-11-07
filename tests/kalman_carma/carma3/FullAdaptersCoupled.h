#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Complete-INS adapters of the coupled three-axis CARMA priors
  (docs/carma-coupled-three-axis.md).  The wrapper, startup, magnetic
  acquisition, tuner and clamps are those of FullCARMA3 (FullAdapters.h); the
  inner filter is Kalman3D_Wave_CARMA3 with a motion block of
  kalman_carma/CarmaCoupledMotion.h.

  Names: <MODEL>:<source>[/opt...] with
    MODEL   MOD  the exact complex-OU modal form of the baseline (12 motion
                 states, real excitation)
            COR  Experiment A, real correlated excitation (12 motion states)
    source  off  no coupling (independent axes: the reduction check)
            est  DEPLOYABLE: the coupling estimator on the measured world-frame
                 acceleration (CarmaCrossSpectrumEstimator.h)
            orc  ORACLE: the same estimator fed with the true wave
                 acceleration (Sensors::a_wave_mag), no shrinkage
    opt     re (co-spectrum only), xy (horizontal anisotropy only, no
            horizontal-vertical term), an (baseline trace with the measured
            horizontal anisotropy), hm (ablation: the measured horizontal
            block scaled by c, the first form of the mapping),
            nz (MOD: directional modal coupling, per-axis numerator zeros),
            g<x> coupling gain, h<s> estimator horizon, b<x> envelope band
            fraction, w<s> warm-up, ns no shrinkage, fz<s> freeze at scenario
            time s (NED), keel (with fz: hold it in the keel frame, turning
            with the hull), flip (oracle: x-z and y-z cross terms negated,
            i.e. the wave direction reversed), rot<deg> (oracle: horizontal
            cross terms turned by deg), ph<deg> (oracle: x-z and y-z cross
            terms phase-shifted by deg), oh<s> oracle horizon (600 s).
  The shared directional modes (Experiment B) are the standalone filter of
  tests/kalman_carma_shared_modes (CsmBenchAdapter.h).
*/

#include "FullAdapters.h"

#include "kalman_carma/CarmaCoupledMotion.h"

#include <complex>
#include <cstdlib>
#include <string>
#include <vector>

namespace carma3_bench {

struct CoupledSpec {
    enum class Source { Off, Est, Oracle } source = Source::Off;
    ocean_imu::kalman::carma::CarmaCouplingConfig cfg{};
    float horizon_s = 240.0f, band_frac = 0.15f, oracle_horizon_s = 600.0f;   // the estimator defaults
    float freeze_at_s = -1.0f;
    bool keel = false, flip = false;
    float rot_deg = 0.0f, phase_deg = 0.0f;
};

inline bool parse_coupled_spec(const std::string& s, CoupledSpec& sp) {
    using Cfg = ocean_imu::kalman::carma::CarmaCouplingConfig;
    std::vector<std::string> tok;
    size_t a = 0;
    while (a <= s.size()) {
        const size_t b = s.find('/', a);
        tok.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) break;
        a = b + 1;
    }
    if (tok.empty()) return false;
    if (tok[0] == "off") sp.source = CoupledSpec::Source::Off;
    else if (tok[0] == "est") sp.source = CoupledSpec::Source::Est;
    else if (tok[0] == "orc") sp.source = CoupledSpec::Source::Oracle;
    else return false;
    for (size_t i = 1; i < tok.size(); ++i) {
        const std::string& o = tok[i];
        auto num = [&](size_t p) { return static_cast<float>(std::strtod(o.c_str() + p, nullptr)); };
        if (o == "re") sp.cfg.form = Cfg::Form::Real;
        else if (o == "an") sp.cfg.horizontal = Cfg::Horizontal::Anisotropic;
        else if (o == "hm") sp.cfg.horizontal = Cfg::Horizontal::Measured;
        else if (o == "xy") sp.cfg.form = Cfg::Form::RealXY;
        else if (o == "nz") sp.cfg.form = Cfg::Form::Numerator;
        else if (o == "ns") sp.cfg.shrink = false;
        else if (o == "keel") sp.keel = true;
        else if (o == "flip") sp.flip = true;
        else if (o.rfind("rot", 0) == 0) sp.rot_deg = num(3);
        else if (o.rfind("ph", 0) == 0) sp.phase_deg = num(2);
        else if (o.rfind("oh", 0) == 0) sp.oracle_horizon_s = num(2);
        else if (o.rfind("fz", 0) == 0) sp.freeze_at_s = num(2);
        else if (o[0] == 'g') sp.cfg.gain = num(1);
        else if (o[0] == 'h') sp.horizon_s = num(1);
        else if (o[0] == 'b') sp.band_frac = num(1);
        else if (o[0] == 'w') sp.cfg.warmup_s = num(1);
        else return false;
    }
    return true;
}

template<class MekfT>
class FullCoupled final : public IEstimator3 {
public:
    using Mekf = MekfT;
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, Mekf>;
    using H3 = ocean_imu::kalman::carma::Hermitian3<float>;
    using Cfg = ocean_imu::kalman::carma::CarmaCouplingConfig;

    FullCoupled(std::string label, const CoupledSpec& sp) : label_(std::move(label)), sp_(sp) {
        ocean_imu::carma::SeaStateFusionConfig_CARMA3 cfg;
        fill_config(cfg, MAG_DELAY_SEC);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.enableTuner(true);
        f.enableClamp(true);
        if constexpr (Mekf::kMotion) {
            auto& m = f.mekf().motion();
            Cfg c = sp_.cfg;
            c.source = sp_.source == CoupledSpec::Source::Off ? Cfg::Source::Off
                     : sp_.source == CoupledSpec::Source::Est ? Cfg::Source::Estimated : Cfg::Source::External;
            c.frame = sp_.keel ? Cfg::Frame::Keel : Cfg::Frame::World;
            if (sp_.source == CoupledSpec::Source::Oracle) c.shrink = false;
            m.coupling().set_config(c);
            typename ocean_imu::kalman::carma::CarmaCrossSpectrumEstimator<float>::Config ec;
            ec.horizon_s = sp_.horizon_s;
            ec.band_frac = sp_.band_frac;
            m.coupling().estimator().configure(ec);
            ec.horizon_s = sp_.oracle_horizon_s;
            oracle_.configure(ec);
        }
    }
    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        temp_ = s.tempC;
        t_ += s.dt;
        mag_.tick(s.dt);
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, s.tempC);
        if (s.mag_tick) {
            fusion_.updateMag(s.mag_ned);
            mag_.after_mag(fusion_, fusion_.raw().mekf().lastMagDiag().accepted);
        }
        if constexpr (Mekf::kMotion) {
            auto& m = fusion_.raw().mekf().motion();
            if (sp_.source == CoupledSpec::Source::Oracle) {
                oracle_.update(s.a_wave_mag.cast<float>(), s.dt, m.operating_frequency_hz());
                if (oracle_.elapsed_s() >= 30.0f) m.coupling().set_external(wrong_(oracle_.cross()), 0.0f);
            }
            if (sp_.freeze_at_s >= 0.0f && t_ >= sp_.freeze_at_s) m.coupling().freeze_now();
        }
    }
    // [dtheta; p; v] covariance in the filter's own coordinates through the
    // readout rows (unit rows: the plain index extraction).
    Eigen::Matrix<double,9,9> native9_() const {
        const auto& m = fusion_.raw().mekf();
        const Eigen::MatrixXd P = m.covariance_full().template cast<double>();
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(9, Mekf::NX);
        for (int i = 0; i < 3; ++i) S(i, i) = 1.0;
        S.bottomRows(6) = m.pv_readout().template cast<double>();
        return S * P * S.transpose();
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.p = m.get_position().template cast<double>();
        e.v = m.get_velocity().template cast<double>();
        e.a = m.get_world_accel().template cast<double>();
        const Eigen::Matrix<double,9,9> M9 = native9_();
        e.Ppv = M9.bottomRightCorner<6,6>();
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        e.has_ext = true;
        Eigen::Matrix<double,9,9> J = Eigen::Matrix<double,9,9>::Zero();
        J.topLeftCorner<3,3>() = e.q_bw.toRotationMatrix();
        J.block<3,3>(3,3) = -Eigen::Matrix3d::Identity();
        J.block<3,3>(6,6) = -Eigen::Matrix3d::Identity();
        e.P9 = J * M9 * J.transpose();
        e.b_g = m.gyroscope_bias().template cast<double>();
        e.b_a = m.get_acc_bias_body_at_temperature(temp_).template cast<double>();
        const auto& d = m.lastAccDiag();
        if (d.accepted) {
            Eigen::LLT<Eigen::Matrix3d> llt(d.S.template cast<double>());
            if (llt.info() == Eigen::Success)
                e.innov_white = llt.matrixL().solve(d.r.template cast<double>());
        }
        return e;
    }
    double native_nees9(const Truth9& tr) const override {
        const auto& m = fusion_.raw().mekf();
        const Eigen::Matrix3d R = m.quaternion_boat().template cast<double>().toRotationMatrix();
        Eigen::Matrix<double,9,1> err;
        const Eigen::AngleAxisd aa(Eigen::Matrix3d(tr.R_bw.transpose() * R));
        err.head<3>() = aa.angle() * aa.axis();
        err.segment<3>(3) = tr.p - m.get_position().template cast<double>();
        err.segment<3>(6) = tr.v - m.get_velocity().template cast<double>();
        Eigen::LDLT<Eigen::Matrix<double,9,9>> l(native9_());
        if (l.info() != Eigen::Success || !l.isPositive()) return NAN;
        return err.dot(l.solve(err));
    }
    double cond_corr() const override { return corr_condition(fusion_.raw().mekf().covariance_full()); }
    long mag_forwarded() const override { return mag_.forwarded; }
    long mag_applied() const override { return mag_.applied; }
    double op_fz_hz() const override {
        const auto& f = fusion_.raw();
        return static_cast<double>(f.getVerticalKappa()) * f.mekf().vertical_model().base_zero_crossing_hz();
    }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return Mekf::NX; }
    bool set_options(const std::string& o) override { return apply_carma_excitation_options(fusion_.raw(), o); }
    void set_label(const std::string& l) override { label_ = l; }
    bool excitation_diag(double (&d)[8]) const override { return carma_excitation_diag(fusion_.raw(), d); }
    bool coupling_diag(double (&d)[6]) const override {
        if constexpr (Mekf::kMotion) {
            const H3& S = fusion_.raw().mekf().motion().excitation();
            const int pr[3][2] = {{0, 2}, {1, 2}, {0, 1}};
            for (int k = 0; k < 3; ++k) {
                const int i = pr[k][0], j = pr[k][1];
                const double den = static_cast<double>(S.re(i,i)) * S.re(j,j);
                const double re = S.re(i,j), im = S.im(i,j);
                d[2 * k] = den > 0 ? (re * re + im * im) / den : NAN;
                d[2 * k + 1] = std::atan2(im, re) * 180.0 / 3.14159265358979323846;
            }
            return true;
        } else {
            (void)d;
            return false;
        }
    }
    Fusion& fusion() { return fusion_; }

private:
    // Deliberately wrong oracle coupling (flip, rot, ph).
    H3 wrong_(H3 S) const {
        if (sp_.flip)
            for (int i = 0; i < 2; ++i) {
                S.re(i,2) = -S.re(i,2); S.im(i,2) = -S.im(i,2);
                S.re(2,i) = -S.re(2,i); S.im(2,i) = -S.im(2,i);
            }
        if (sp_.rot_deg != 0.0f) {
            const float a = sp_.rot_deg * 3.14159265f / 180.0f, c = std::cos(a), sn = std::sin(a);
            Eigen::Matrix3f R;
            R << c, -sn, 0.0f, sn, c, 0.0f, 0.0f, 0.0f, 1.0f;
            S.re = R * S.re * R.transpose();
            S.im = R * S.im * R.transpose();
        }
        if (sp_.phase_deg != 0.0f) {
            const std::complex<float> e = std::polar(1.0f, sp_.phase_deg * 3.14159265f / 180.0f);
            for (int i = 0; i < 2; ++i) {
                const std::complex<float> v = e * std::complex<float>(S.re(i,2), S.im(i,2));
                S.re(i,2) = v.real(); S.im(i,2) = v.imag();
                S.re(2,i) = v.real(); S.im(2,i) = -v.imag();
            }
        }
        return S;
    }

    std::string label_;
    CoupledSpec sp_;
    Fusion fusion_;
    ocean_imu::kalman::carma::CarmaCrossSpectrumEstimator<float> oracle_{};
    MagCounter mag_;
    float temp_ = 35.0f;
    float t_ = 0.0f;
};

}  // namespace carma3_bench
