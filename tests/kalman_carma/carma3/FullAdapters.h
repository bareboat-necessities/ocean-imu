#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Complete-INS adapters of carma3_bench.  Every adapter runs the estimator
  behind the same proxy startup and magnetic acquisition, with the shipping
  simulator configuration (Bench3.h), the tuner and its clamps, and for
  OU-III its periodic a_w covariance sync.
*/

#include "Bench3.h"

#include "util/W3dSimCommon.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"

#include <Eigen/Eigenvalues>

namespace carma3_bench {

// Extended diagnostics shared by every adapter.

// Magnetometer samples that the wrapper forwarded to the estimator (the
// proxy's own gate: delay elapsed, north lock, Live) and the updates it
// applied.  Call after fusion.updateMag().
struct MagCounter {
    float t = 0.0f;
    long forwarded = 0, applied = 0;
    void tick(float dt) { t += dt; }
    template<class Fusion>
    void after_mag(const Fusion& fusion, bool accepted) {
        if (t < MAG_DELAY_SEC || !fusion.hasMagNorthLock() || !fusion.isLive()) return;
        ++forwarded;
        if (accepted) ++applied;
    }
};

template<class Mat>
double corr_condition(const Mat& Pin) {
    const Eigen::MatrixXd P = Pin.template cast<double>();
    const Eigen::VectorXd d = P.diagonal().cwiseMax(1e-300).cwiseSqrt();
    const Eigen::MatrixXd C = d.cwiseInverse().asDiagonal() * P * d.cwiseInverse().asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C, Eigen::EigenvaluesOnly);
    const double lo = es.eigenvalues().minCoeff(), hi = es.eigenvalues().maxCoeff();
    return lo > 0 ? hi / lo : INFINITY;
}

// Quaternion-MEKF error e = [dtheta_b (C_true = Exp(dtheta) Chat), p - phat,
// v - vhat]: the common e9 = [Log(Rhat R^T), phat - p, vhat - v] covariance is
// J P J^T with J = blkdiag(Rhat, -I, -I).
template<class Mat>
Eigen::Matrix<double,9,9> mekf_P9(const Mat& P, const Eigen::Matrix3d& R, const int (&ip)[3], const int (&iv)[3]) {
    int idx[9] = {0, 1, 2, ip[0], ip[1], ip[2], iv[0], iv[1], iv[2]};
    Eigen::Matrix<double,9,9> M;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) M(i,j) = static_cast<double>(P(idx[i], idx[j]));
    Eigen::Matrix<double,9,9> J = Eigen::Matrix<double,9,9>::Zero();
    J.topLeftCorner<3,3>() = R;
    J.block<3,3>(3,3) = -Eigen::Matrix3d::Identity();
    J.block<3,3>(6,6) = -Eigen::Matrix3d::Identity();
    return J * M * J.transpose();
}
template<class Mat>
double mekf_native_nees9(const Mat& P, const Eigen::Matrix3d& R, const Eigen::Vector3d& p, const Eigen::Vector3d& v,
                         const Truth9& tr, const int (&ip)[3], const int (&iv)[3]) {
    int idx[9] = {0, 1, 2, ip[0], ip[1], ip[2], iv[0], iv[1], iv[2]};
    Eigen::Matrix<double,9,9> M;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) M(i,j) = static_cast<double>(P(idx[i], idx[j]));
    Eigen::Matrix<double,9,1> e;
    const Eigen::AngleAxisd aa(Eigen::Matrix3d(tr.R_bw.transpose() * R));
    e.head<3>() = aa.angle() * aa.axis();
    e.segment<3>(3) = tr.p - p;
    e.segment<3>(6) = tr.v - v;
    Eigen::LDLT<Eigen::Matrix<double,9,9>> l(M);
    if (l.info() != Eigen::Success || !l.isPositive()) return NAN;
    return e.dot(l.solve(e));
}

// Shipping OU-III, unchanged.
class FullOU3 final : public IEstimator3 {
public:
    explicit FullOU3(std::string label = "OU3") : label_(std::move(label)) {
        SeaStateFusionConfig_OU_III cfg;
        fill_config(cfg, MAG_DELAY_SEC);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.setPeriodicAwCovarianceSync(true);
        f.setAwCovarianceSyncCongruent(false);
        f.enableTuner(true);
        f.enableClamp(true);
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
        e.p = m.get_position().cast<double>();
        e.v = m.get_velocity().cast<double>();
        e.a = m.get_world_accel().cast<double>();
        e.Ppv = extract_pv(m.covariance_full(), kIp, kIv);
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().cast<double>();
        e.has_ext = true;
        e.P9 = mekf_P9(m.covariance_full(), e.q_bw.toRotationMatrix(), kIp, kIv);
        e.b_g = m.gyroscope_bias().cast<double>();
        e.b_a = m.get_acc_bias_body_at_temperature(temp_).cast<double>();
        return e;
    }
    double native_nees9(const Truth9& tr) const override {
        const auto& m = fusion_.raw().mekf();
        return mekf_native_nees9(m.covariance_full(), m.quaternion_boat().cast<double>().toRotationMatrix(),
                                 m.get_position().cast<double>(), m.get_velocity().cast<double>(), tr, kIp, kIv);
    }
    double cond_corr() const override { return corr_condition(fusion_.raw().mekf().covariance_full()); }
    long mag_forwarded() const override { return mag_.forwarded; }
    long mag_applied() const override { return mag_.applied; }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return 21; }

private:
    // [dtheta, b_g, v, p, S, a_w, b_a]
    static constexpr int kIp[3] = {9, 10, 11}, kIv[3] = {6, 7, 8};
    std::string label_;
    SeaStateFusion_OU_III<TrackerType::KALMANF> fusion_;
    MagCounter mag_;
    float temp_ = 35.0f;
};

// Three-axis CARMA (Kalman3D_Wave_CARMA3).
template<int NHC, bool SH, int NZC, bool SZ>
class FullCARMA3 final : public IEstimator3 {
public:
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, NHC, SH, NZC, SZ>;
    using M = typename Fusion::Impl::Mekf;

    explicit FullCARMA3(std::string label, const ocean_imu::kalman::carma::CarmaShape* hshape = nullptr,
                        float h_corr = 0.0f)
        : label_(std::move(label)) {
        ocean_imu::carma::SeaStateFusionConfig_CARMA3 cfg;
        fill_config(cfg, MAG_DELAY_SEC);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.enableTuner(true);
        f.enableClamp(true);
        if (hshape) f.setHorizontalShape(*hshape);
        if (h_corr != 0.0f) f.setHorizontalCorrelation(h_corr);
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
        e.Ppv = extract_pv(m.covariance_full(), kIp, kIv);
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        e.has_ext = true;
        e.P9 = mekf_P9(m.covariance_full(), e.q_bw.toRotationMatrix(), kIp, kIv);
        e.b_g = m.gyroscope_bias().template cast<double>();
        e.b_a = m.get_acc_bias_body_at_temperature(temp_).template cast<double>();
        return e;
    }
    double native_nees9(const Truth9& tr) const override {
        const auto& m = fusion_.raw().mekf();
        return mekf_native_nees9(m.covariance_full(), m.quaternion_boat().template cast<double>().toRotationMatrix(),
                                 m.get_position().template cast<double>(), m.get_velocity().template cast<double>(),
                                 tr, kIp, kIv);
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
    static constexpr int kIp[3] = {M::pos_index(0), M::pos_index(1), M::pos_index(2)};
    static constexpr int kIv[3] = {M::vel_index(0), M::vel_index(1), M::vel_index(2)};
    std::string label_;
    Fusion fusion_;
    MagCounter mag_;
    float temp_ = 35.0f;
};

// Horizontal CARMA(4,1) shape fitted (tools/carma3_analysis.py, energetic
// band) to JONSWAP 0.29 Hz times the benchmark's own beam-sea sway RAO: it
// depends on that hull model.  (f1, z1, f2, z2, fl, fz, zz).
inline constexpr ocean_imu::kalman::carma::CarmaShape CARMA41_SWAY_SHAPE{
    0.283900, 0.072900, 0.453000, 0.255200, 0.0, 0.0, 0.0};

}  // namespace carma3_bench
