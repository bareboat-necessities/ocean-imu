#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Complete-INS adapters of carma3_bench.  Every adapter runs the estimator
  behind the same proxy startup and magnetic acquisition, with the shipping
  simulator configuration (Bench3.h), the periodic a_w covariance sync, the
  tuner and its clamps -- the configuration carma_bench uses.
*/

#include "Bench3.h"

#include "util/W3dSimCommon.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA.h"
#include "SeaStateFusionFilter_CARMA3.h"

namespace carma3_bench {

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
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, 35.0f);
        if (s.mag_tick) fusion_.updateMag(s.mag_ned);
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.p = m.get_position().cast<double>();
        e.v = m.get_velocity().cast<double>();
        e.a = m.get_world_accel().cast<double>();
        // [dtheta, b_g, v, p, S, a_w, b_a]
        constexpr int ip[3] = {9, 10, 11}, iv[3] = {6, 7, 8};
        e.Ppv = extract_pv(m.covariance_full(), ip, iv);
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().cast<double>();
        return e;
    }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return 21; }

private:
    std::string label_;
    SeaStateFusion_OU_III<TrackerType::KALMANF> fusion_;
};

// The PR #704 hybrid (OU-III horizontal chain + CARMA heave), unchanged.
template<int N, bool VI>
class FullHybrid final : public IEstimator3 {
public:
    explicit FullHybrid(std::string label) : label_(std::move(label)) {
        ocean_imu::carma::SeaStateFusionConfig_CARMA cfg;
        fill_config(cfg, MAG_DELAY_SEC);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.setPeriodicAwCovarianceSync(true);
        f.enableTuner(true);
        f.enableClamp(true);
    }
    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, 35.0f);
        if (s.mag_tick) fusion_.updateMag(s.mag_ned);
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = fusion_.raw().mekf();
        using M = typename Fusion::Impl::Mekf;
        e.live = fusion_.isLive();
        e.p = m.get_position().template cast<double>();
        e.v = m.get_velocity().template cast<double>();
        e.a = m.get_world_accel().template cast<double>();
        constexpr int ip[3] = {M::OFF_PH, M::OFF_PH + 1, M::OFF_Z};
        constexpr int iv[3] = {M::OFF_VH, M::OFF_VH + 1, M::OFF_Z + 1};
        e.Ppv = extract_pv(m.covariance_full(), ip, iv);
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        return e;
    }
    double op_fz_hz() const override {
        const auto& f = fusion_.raw();
        return static_cast<double>(f.getVerticalKappa()) * f.mekf().vertical_model().base_zero_crossing_hz();
    }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return Fusion::Impl::Mekf::NX; }

private:
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA<TrackerType::KALMANF, N, VI>;
    std::string label_;
    Fusion fusion_;
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
        f.setPeriodicAwCovarianceSync(true);
        f.enableTuner(true);
        f.enableClamp(true);
        if (hshape) f.setHorizontalShape(*hshape);
        if (h_corr != 0.0f) f.setHorizontalCorrelation(h_corr);
    }
    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, 35.0f);
        if (s.mag_tick) fusion_.updateMag(s.mag_ned);
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.p = m.get_position().template cast<double>();
        e.v = m.get_velocity().template cast<double>();
        e.a = m.get_world_accel().template cast<double>();
        constexpr int ip[3] = {M::pos_index(0), M::pos_index(1), M::pos_index(2)};
        constexpr int iv[3] = {M::vel_index(0), M::vel_index(1), M::vel_index(2)};
        e.Ppv = extract_pv(m.covariance_full(), ip, iv);
        e.nis = m.lastAccDiag().nis;
        e.has_att = true;
        e.q_bw = m.quaternion_boat().template cast<double>();
        return e;
    }
    double op_fz_hz() const override {
        const auto& f = fusion_.raw();
        return static_cast<double>(f.getVerticalKappa()) * f.mekf().vertical_model().base_zero_crossing_hz();
    }
    size_t memory_bytes() const override { return sizeof(fusion_) + sizeof(fusion_.raw().mekf()); }
    int states() const override { return M::NX; }
    Fusion& fusion() { return fusion_; }

private:
    std::string label_;
    Fusion fusion_;
};

// Horizontal CARMA(4,1) shape fitted (tools/carma3_analysis.py, energetic
// band) to JONSWAP 0.29 Hz times the benchmark's own beam-sea sway RAO: it
// depends on that hull model.  (f1, z1, f2, z2, fl, fz, zz).
inline constexpr ocean_imu::kalman::carma::CarmaShape CARMA41_SWAY_SHAPE{
    0.283900, 0.072900, 0.453000, 0.255200, 0.0, 0.0, 0.0};

}  // namespace carma3_bench
