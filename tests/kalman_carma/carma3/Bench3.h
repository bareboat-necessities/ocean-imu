#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Shared types of carma3_bench (three-axis CARMA benchmark): the sensor
  sample every estimator receives, the three-axis estimate it reports, the
  estimator interface and the factory of the complete-INS adapters (defined
  in carma3_filters_*.cpp so the template instantiations compile in
  parallel).
*/

#define EIGEN_NON_ARDUINO

#include <memory>
#include <string>

#include <Eigen/Dense>

namespace carma3_bench {

struct Sensors {
    float dt = 0.005f;
    Eigen::Vector3f gyr_ned = Eigen::Vector3f::Zero(), acc_ned = Eigen::Vector3f::Zero();
    bool mag_tick = false;
    Eigen::Vector3f mag_ned = Eigen::Vector3f::Zero();
    Eigen::Vector3d y_level = Eigen::Vector3d::Zero();   // level 1: R_nb f + g e_z (NED)
    float tempC = 35.0f;                                  // IMU temperature passed to every filter
};

struct Estimate3 {
    bool live = false;
    // NED, in the estimator's world frame (magnetic north for the complete
    // INS; the caller rotates into the truth's geographic frame).
    Eigen::Vector3d p = Eigen::Vector3d::Zero(), v = Eigen::Vector3d::Zero(), a = Eigen::Vector3d::Zero();
    // Covariance of [p_x, p_y, p_z, v_x, v_y, v_z].
    Eigen::Matrix<double,6,6> Ppv = Eigen::Matrix<double,6,6>::Constant(NAN);
    double nis = NAN;
    int nis_dof = 3;
    bool has_att = false;
    Eigen::Quaterniond q_bw = Eigen::Quaterniond::Identity();   // body -> world (NED)
    // Extended diagnostics (carma3_bench --extended), in the estimator's world
    // frame: covariance of e9 = [Log(Rhat R^T), phat - p, vhat - v] (common
    // coordinates for every filter) and the body-frame bias estimates (the
    // accelerometer bias at the sample temperature).
    bool has_ext = false;
    Eigen::Matrix<double,9,9> P9 = Eigen::Matrix<double,9,9>::Constant(NAN);
    Eigen::Vector3d b_g = Eigen::Vector3d::Constant(NAN), b_a = Eigen::Vector3d::Constant(NAN);
};

// Truth in an estimator's world frame, for the native-coordinate NEES.
struct Truth9 {
    Eigen::Matrix3d R_bw = Eigen::Matrix3d::Identity();
    Eigen::Vector3d p = Eigen::Vector3d::Zero(), v = Eigen::Vector3d::Zero();
};

class IEstimator3 {
public:
    virtual ~IEstimator3() = default;
    virtual std::string name() const = 0;
    virtual void step(const Sensors& s) = 0;
    virtual Estimate3 estimate() const = 0;
    virtual size_t memory_bytes() const = 0;
    virtual int states() const = 0;
    // CARMA operating frequency kappa * f_z,base (the tuner's f_tune); NaN
    // for OU-III, which has no CARMA model.  Diagnostic output only.
    virtual double op_fz_hz() const { return NAN; }
    // Joint NEES of [attitude, p, v] in the estimator's OWN error coordinates
    // (MEKF: dtheta_b, additive p, v; Lie group: phi, rho_p, rho_v).
    virtual double native_nees9(const Truth9&) const { return NAN; }
    // Correlation-scaled condition number of the full covariance.
    virtual double cond_corr() const { return NAN; }
    // Magnetometer samples that reached the estimator and the updates it
    // applied (rejected = forwarded - applied).
    virtual long mag_forwarded() const { return -1; }
    virtual long mag_applied() const { return -1; }
};

// Shipping simulator configuration (as kalman_ou_iii-sim).
constexpr float SIGMA_A_RESCALE = 0.71f;
constexpr float SIGMA_G_RESCALE = 0.01f;
constexpr float SIGMA_M_RESCALE = 8.0f;
constexpr float ACC_SIGMA = 1.51e-3f * 9.80665f;
constexpr float GYR_SIGMA = 0.00157f;
constexpr float MAG_SIGMA_UT = 0.80f;

template<class Cfg>
void fill_config(Cfg& cfg, float mag_delay_sec) {
    cfg.with_mag = true;
    cfg.sigma_a = Eigen::Vector3f::Constant(2.8f * ACC_SIGMA * SIGMA_A_RESCALE);
    cfg.sigma_g = Eigen::Vector3f::Constant(2.0f * GYR_SIGMA * SIGMA_G_RESCALE);
    cfg.sigma_m = Eigen::Vector3f::Constant(1.2f * MAG_SIGMA_UT * SIGMA_M_RESCALE);
    cfg.mag_delay_sec = mag_delay_sec;
}

// [p; v] covariance from a full covariance and index maps.
template<class Mat>
Eigen::Matrix<double,6,6> extract_pv(const Mat& P, const int (&ip)[3], const int (&iv)[3]) {
    int idx[6] = {ip[0], ip[1], ip[2], iv[0], iv[1], iv[2]};
    Eigen::Matrix<double,6,6> M;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) M(i,j) = static_cast<double>(P(idx[i], idx[j]));
    return M;
}

// Factories: nullptr for a name the unit does not define.
std::unique_ptr<IEstimator3> make_full_ou(const std::string& name);
std::unique_ptr<IEstimator3> make_full_c4(const std::string& name);
std::unique_ptr<IEstimator3> make_full_c5(const std::string& name);
std::unique_ptr<IEstimator3> make_full_tfg(const std::string& name);

inline std::unique_ptr<IEstimator3> make_full(const std::string& name) {
    if (auto e = make_full_ou(name)) return e;
    if (auto e = make_full_c4(name)) return e;
    if (auto e = make_full_tfg(name)) return e;
    return make_full_c5(name);
}

}  // namespace carma3_bench
