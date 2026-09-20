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
template<class A, class B, class C, class D, class E>
static void readout_prediction(const A& f, const B& fl, const C& q,
                               const D& ql, float phi, const E& qb) {
    if (!recording) return;
    std::ostringstream out;
    out << std::setprecision(17) << "{\"kind\":\"prediction\",\"F_AG\":"
        << matrix_json(f) << ",\"F_LIN\":" << matrix_json(fl)
        << ",\"Q_AG\":" << matrix_json(q) << ",\"Q_LIN\":" << matrix_json(ql)
        << ",\"phi_BA\":" << phi << ",\"Q_BA\":" << matrix_json(qb) << '}';
    events.push_back(out.str());
}
template<class A> static void readout_sync(const A& q) {
    if (recording) events.push_back("{\"kind\":\"sync\",\"Q\":"+matrix_json(q)+'}');
}
template<class A, class B, class C>
static void readout_correction(const char* sensor, const A& h, const B& r, const C& k) {
    if (recording) events.push_back(std::string("{\"kind\":\"correction\",\"sensor\":\"")+sensor
        +"\",\"H\":"+matrix_json(h)+",\"R\":"+matrix_json(r)+",\"K\":"+matrix_json(k)+'}');
}
template<class A> static void readout_reset(const A& d) {
    if (recording) events.push_back("{\"kind\":\"reset\",\"d\":"+matrix_json(d)+'}');
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
        const float roll = wave ? static_cast<float>(.02*std::sin(.5*t)) : 0.0f;
        const float rate = wave ? static_cast<float>(.01*std::cos(.5*t)) : 0.0f;
        const float az = wave ? static_cast<float>(-.144*std::sin(.6*t)) : 0.0f;
        const Eigen::Matrix3f rwb=Eigen::AngleAxisf(-roll,Eigen::Vector3f::UnitX()).toRotationMatrix();
        filter.update(.005f,Eigen::Vector3f(rate,0,0),rwb*Eigen::Vector3f(0,0,az-g_std));
        if (k%8==0) {
            const Eigen::Vector3f measured_field = wave ? (rwb*Eigen::Vector3f(60,0,30)).eval() : field;
            filter.updateMag(measured_field);
            if (recording && filter.raw().mekf().lastMagDiag().accepted) ++applied;
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
