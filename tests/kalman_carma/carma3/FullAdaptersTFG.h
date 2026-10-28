#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Complete-INS adapter of the two-frame Lie-group CARMA filter
  (Kalman3D_Wave_TFG_CARMA3) for carma3_bench.  The wrapper is the CARMA3
  wrapper with the inner filter replaced (SeaStateFusionFilter_TFG_CARMA3.h),
  configured exactly as FullCARMA3, so the pair differs only in the error
  geometry.  The reset and covariance-write modes select the ablations.
*/

#include "FullAdapters.h"
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"

namespace carma3_bench {

template<int NC>
class FullTFGCARMA3 final : public IEstimator3 {
public:
    using Fusion = ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF, NC>;
    using M = typename Fusion::Impl::Mekf;

    FullTFGCARMA3(std::string label, typename M::ResetMode reset, typename M::CovarianceWriteMode write)
        : label_(std::move(label)) {
        ocean_imu::carma::SeaStateFusionConfig_CARMA3 cfg;
        fill_config(cfg, MAG_DELAY_SEC);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.enableTuner(true);
        f.enableClamp(true);
        // The filter is constructed at the identity state, where both write
        // modes coincide, so setting them after begin() is exact.
        f.mekf().set_reset_mode(reset);
        f.mekf().set_covariance_write_mode(write);
    }
    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        temp_ = s.tempC;
        mag_.tick(s.dt);
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, s.tempC);
        if (s.mag_tick) {
            fusion_.updateMag(s.mag_ned);
            mag_.after_mag(fusion_, fusion_.raw().mekf().lastMagDiag().accepted);
        }
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.p = m.get_position().template cast<double>();
        e.v = m.get_velocity().template cast<double>();
        e.a = m.get_world_accel().template cast<double>();
        // Additive world errors: dp = -rho_p + [phat]x phi (see the filter header).
        e.Ppv = m.physical_pv_covariance().template cast<double>();
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        e.has_ext = true;
        // e9 = [phi, phat - p, vhat - v] = J [phi, rho_p, rho_v], first order:
        // phat - p = rho_p - [phat]x phi.
        Eigen::Matrix<double,9,M::NX> J = Eigen::Matrix<double,9,M::NX>::Zero();
        J.template block<3,3>(0, M::OFF_PHI) = Eigen::Matrix3d::Identity();
        J.template block<3,3>(3, M::OFF_P) = Eigen::Matrix3d::Identity();
        J.template block<3,3>(3, M::OFF_PHI) = -skew_(e.p);
        J.template block<3,3>(6, M::OFF_V) = Eigen::Matrix3d::Identity();
        J.template block<3,3>(6, M::OFF_PHI) = -skew_(e.v);
        const Eigen::Matrix<double,M::NX,M::NX> P = m.covariance_full().template cast<double>();
        e.P9 = J * P * J.transpose();
        e.b_g = m.gyroscope_bias().template cast<double>();
        e.b_a = m.get_acc_bias_body_at_temperature(temp_).template cast<double>();
        return e;
    }
    // Native coordinates xi = Log_G(Xhat X^-1): phi, rho = J_l(phi)^-1 (xhat - Exp(phi) x).
    double native_nees9(const Truth9& tr) const override {
        const auto& m = fusion_.raw().mekf();
        const Eigen::Matrix3d R = m.state().R.template cast<double>();
        const Eigen::Vector3d phi = ocean_imu::lie::Log<double>(Eigen::Matrix3d(R * tr.R_bw.transpose()));
        const Eigen::Matrix3d dR = ocean_imu::lie::Exp<double>(phi);
        const Eigen::Matrix3d Jli = ocean_imu::lie::inv_left_jacobian<double>(phi);
        Eigen::Matrix<double,9,1> x;
        x.head<3>() = phi;
        x.segment<3>(3) = Jli * (m.get_position().template cast<double>() - dR * tr.p);
        x.segment<3>(6) = Jli * (m.get_velocity().template cast<double>() - dR * tr.v);
        const int idx[9] = {0, 1, 2, M::OFF_P, M::OFF_P + 1, M::OFF_P + 2, M::OFF_V, M::OFF_V + 1, M::OFF_V + 2};
        Eigen::Matrix<double,9,9> Pn;
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) Pn(i,j) = static_cast<double>(m.covariance_full()(idx[i], idx[j]));
        Eigen::LDLT<Eigen::Matrix<double,9,9>> l(Pn);
        if (l.info() != Eigen::Success || !l.isPositive()) return NAN;
        return x.dot(l.solve(x));
    }
    double cond_corr() const override { return corr_condition(fusion_.raw().mekf().covariance_full()); }
    long mag_forwarded() const override { return mag_.forwarded; }
    long mag_applied() const override { return mag_.applied; }
    double op_fz_hz() const override {
        const auto& f = fusion_.raw();
        return static_cast<double>(f.getVerticalKappa()) * f.mekf().vertical_model().base_zero_crossing_hz();
    }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return M::NX; }
    Fusion& fusion() { return fusion_; }

private:
    static Eigen::Matrix3d skew_(const Eigen::Vector3d& v) {
        Eigen::Matrix3d S;
        S << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
        return S;
    }
    std::string label_;
    Fusion fusion_;
    MagCounter mag_;
    float temp_ = 35.0f;
};

}  // namespace carma3_bench
