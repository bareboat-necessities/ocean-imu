// Read-only observer driver for the optional historical-readout diagnostic.
// The Python driver inserts taps into a temporary copy of the shipping header,
// compiles an untapped control, and requires identical terminal state/covariance.
// No setter/reseed is used after construction. This finite replay is not an
// all-time physical/magnetic-service or real-arithmetic certificate.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

static bool recording = false;
static std::vector<std::string> events;
static double physical_t = 0.0;
static Eigen::Vector3d physical_p = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_v = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_S = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_a = Eigen::Vector3d::Zero();
static Eigen::Vector4d physical_quat = Eigen::Vector4d(0,0,0,1);
static Eigen::Vector3d physical_bg = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_ba = Eigen::Vector3d::Zero();
static double physical_t_prev = 0.0;
static Eigen::Vector3d physical_p_prev = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_v_prev = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_S_prev = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_a_prev = Eigen::Vector3d::Zero();
static Eigen::Vector4d physical_quat_prev = Eigen::Vector4d(0,0,0,1);
static Eigen::Vector3d physical_bg_prev = Eigen::Vector3d::Zero();
static Eigen::Vector3d physical_ba_prev = Eigen::Vector3d::Zero();
static Eigen::Matrix<double,21,1> estimator_state = Eigen::Matrix<double,21,1>::Zero();
static Eigen::Vector4d estimator_quat = Eigen::Vector4d(0,0,0,1);
static Eigen::Matrix<double,21,1> estimator_state_after = Eigen::Matrix<double,21,1>::Zero();
static Eigen::Vector4d estimator_quat_after = Eigen::Vector4d(0,0,0,1);
template<class A> static std::string matrix_json(const A& a) {
    std::ostringstream out;
    out << std::setprecision(17) << '[';
    for (int i=0; i<a.rows(); ++i) {
        if (i) out << ',';
        out << '[';
        for (int j=0; j<a.cols(); ++j) {
            if (j) out << ',';
            out << static_cast<double>(a(i,j));
        }
        out << ']';
    }
    out << ']';
    return out.str();
}
static std::string estimator_json() {
    return std::string(",\"estimator_state\":")+matrix_json(estimator_state)
        +",\"estimator_quaternion\":"+matrix_json(estimator_quat)
        +",\"estimator_state_after\":"+matrix_json(estimator_state_after)
        +",\"estimator_quaternion_after\":"+matrix_json(estimator_quat_after);
}
static std::string physical_json() {
    return std::string(",\"physical_t\":")+std::to_string(physical_t)
        +",\"physical_p\":"+matrix_json(physical_p)
        +",\"physical_v\":"+matrix_json(physical_v)
        +",\"physical_S\":"+matrix_json(physical_S)
        +",\"physical_a\":"+matrix_json(physical_a)
        +",\"physical_quaternion\":"+matrix_json(physical_quat)
        +",\"physical_bg\":"+matrix_json(physical_bg)
        +",\"physical_ba\":"+matrix_json(physical_ba)+estimator_json();
}
static std::string prediction_physical_json() {
    return std::string(",\"physical_t\":")+std::to_string(physical_t_prev)
        +",\"physical_p\":"+matrix_json(physical_p_prev)
        +",\"physical_v\":"+matrix_json(physical_v_prev)
        +",\"physical_S\":"+matrix_json(physical_S_prev)
        +",\"physical_a\":"+matrix_json(physical_a_prev)
        +",\"physical_quaternion\":"+matrix_json(physical_quat_prev)
        +",\"physical_bg\":"+matrix_json(physical_bg_prev)
        +",\"physical_ba\":"+matrix_json(physical_ba_prev)
        +",\"physical_t_after\":"+std::to_string(physical_t)
        +",\"physical_p_after\":"+matrix_json(physical_p)
        +",\"physical_v_after\":"+matrix_json(physical_v)
        +",\"physical_S_after\":"+matrix_json(physical_S)
        +",\"physical_a_after\":"+matrix_json(physical_a)
        +",\"physical_quaternion_after\":"+matrix_json(physical_quat)
        +",\"physical_bg_after\":"+matrix_json(physical_bg)
        +",\"physical_ba_after\":"+matrix_json(physical_ba)
        +estimator_json();
}
template<class A, class B, class C, class D, class E>
static void readout_prediction(const A& f, const B& fl, const C& q,
                               const D& ql, float phi, const E& qb) {
    if (!recording) return;
    std::ostringstream out;
    out << std::setprecision(17) << "{\"kind\":\"prediction\",\"F_AG\":"
        << matrix_json(f) << ",\"F_LIN\":" << matrix_json(fl)
        << ",\"Q_AG\":" << matrix_json(q) << ",\"Q_LIN\":" << matrix_json(ql)
        << ",\"phi_BA\":" << phi << ",\"Q_BA\":" << matrix_json(qb)
        << prediction_physical_json() << '}';
    events.push_back(out.str());
}
template<class A> static void readout_sync(const A& q) {
    if (recording) events.push_back("{\"kind\":\"sync\",\"Q\":"+matrix_json(q)+physical_json()+'}');
}
template<class A, class B, class C>
static void readout_correction(const char* sensor, const A& h, const B& r, const C& k) {
    if (recording) events.push_back(std::string("{\"kind\":\"correction\",\"sensor\":\"")+sensor
        +"\",\"H\":"+matrix_json(h)+",\"R\":"+matrix_json(r)+",\"K\":"+matrix_json(k)
        +physical_json()+'}');
}
template<class A> static void readout_reset(const A& d) {
    if (recording) events.push_back("{\"kind\":\"reset\",\"d\":"+matrix_json(d)+physical_json()+'}');
}

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

int main(int argc, char** argv) {
    const bool wave = argc > 1 && std::string(argv[1])=="wave";
    const float heading = argc > 1 && !wave ? std::stof(argv[1]) : 0.0f;
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.12f);
    cfg.sigma_g.setConstant(.00135f);
    cfg.sigma_m.setConstant(.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    const Eigen::Vector3f field(75.0f*std::cos(heading),75.0f*std::sin(heading),0.0f);
    int live=-1, refined=-1, active=-1, applied=0;
    std::string root;
    for (int k=1; k<=45064; ++k) {
        if (k==45001) {
            if (active<0 || k-active<3400) return 2;
            root=matrix_json(filter.raw().mekf().covariance_full());
            recording=true;
        }
        const double t = static_cast<double>(k)*.005;
        physical_t_prev = physical_t;
        physical_p_prev = physical_p;
        physical_v_prev = physical_v;
        physical_S_prev = physical_S;
        physical_a_prev = physical_a;
        physical_quat_prev = physical_quat;
        physical_bg_prev = physical_bg;
        physical_ba_prev = physical_ba;
        physical_t = t;
        if (wave) {
            physical_p = Eigen::Vector3d(0.0,0.0,.4*std::sin(.6*t));
            physical_v = Eigen::Vector3d(0.0,0.0,.24*std::cos(.6*t));
            physical_a = Eigen::Vector3d(0.0,0.0,-.144*std::sin(.6*t));
            physical_S = Eigen::Vector3d(0.0,0.0,(.4/.6)*(1.0-std::cos(.6*t)));
        } else {
            physical_p.setZero(); physical_v.setZero(); physical_a.setZero(); physical_S.setZero();
        }
        const float roll = wave ? static_cast<float>(.02*std::sin(.5*t)) : 0.0f;
        const float rate = wave ? static_cast<float>(.01*std::cos(.5*t)) : 0.0f;
        const Eigen::Quaterniond qbw(Eigen::AngleAxisd(static_cast<double>(roll),Eigen::Vector3d::UnitX()));
        physical_quat = qbw.coeffs();
        const float az = wave ? static_cast<float>(-.144*std::sin(.6*t)) : 0.0f;
        const Eigen::Matrix3f rwb=Eigen::AngleAxisf(-roll,Eigen::Vector3f::UnitX()).toRotationMatrix();
        estimator_state = filter.raw().mekf().xext.template cast<double>();
        estimator_quat = filter.raw().mekf().qref.coeffs().template cast<double>();
        filter.update(.005f,Eigen::Vector3f(rate,0,0),rwb*Eigen::Vector3f(0,0,az-g_std));
        if (k%8==0) {
            const Eigen::Vector3f measured_field = wave ? (rwb*Eigen::Vector3f(60,0,30)).eval() : field;
            estimator_state = filter.raw().mekf().xext.template cast<double>();
            estimator_quat = filter.raw().mekf().qref.coeffs().template cast<double>();
            filter.updateMag(measured_field);
            if (recording && filter.raw().mekf().lastMagDiag().accepted) ++applied;
        }
        if (recording) {
            const auto& raw = filter.raw();
            estimator_state = raw.mekf().xext.template cast<double>();
            estimator_quat = raw.mekf().qref.coeffs().template cast<double>();
            std::ostringstream tune;
            tune << std::setprecision(17)
                 << "{\"kind\":\"adaptive_state\",\"physical_t\":" << physical_t
                 << ",\"tau\":" << raw.getTauApplied()
                 << ",\"sigma_aw\":" << raw.getSigmaApplied()
                 << ",\"R_S\":" << raw.getRSApplied()
                 << ",\"T_S\":" << raw.getPseudoUpdatePeriodSec()
                 << ",\"tau_target\":" << raw.getTauTarget()
                 << ",\"sigma_target\":" << raw.getSigmaTarget()
                 << ",\"variance\":" << raw.getAccelVariance()
                 << ",\"variance_horizon\":" << raw.getSigmaVarianceHorizonSec()
                 << ",\"freq\":" << raw.getFreqHz()
                 << ",\"period\":" << raw.getPeriodSec()
                 << ",\"proxy_q\":" << matrix_json(raw.startupProxyQuat().coeffs())
                 << physical_json() << "}";
            events.push_back(tune.str());
        }
        if (live<0 && filter.isLive()) live=k;
        if (refined<0 && filter.hasRefinedMagReference()) refined=k;
        if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;
    }
    const auto& m=filter.raw().mekf();
    if (!m.Pext.allFinite() || !m.xext.allFinite() || applied!=8) return 3;
    std::cout << "{\"live_step\":" << live << ",\"refined_step\":" << refined
              << ",\"active_step\":" << active << ",\"applied_magnetic_updates\":" << applied
              << ",\"root_covariance\":" << root
              << ",\"terminal_covariance\":" << matrix_json(m.Pext)
              << ",\"terminal_state\":" << matrix_json(m.xext)
              << ",\"terminal_quaternion\":" << matrix_json(m.qref.coeffs())
              << ",\"events\":[";
    for (std::size_t i=0; i<events.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << events[i];
    }
    std::cout << "]}\n";
}
