/* Copyright (c) 2026 Mikhail Grushinskiy
 * Device-configured, deterministic ten-minute histories. No detrender is used:
 * the primary regression measures the actual estimator's raw world state.
 */
#define EIGEN_NON_ARDUINO
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
using V = Eigen::Vector3f;
using M = Eigen::Matrix3f;
constexpr float g = 9.80665f;
constexpr float pi = 3.14159265358979323846f;
struct Case {
    const char* name;
    float residual, dt;
    bool tilted=false, hold_bias=false, fixed=false, sync=true, moving=false, noise=false;
};

bool run(const Case& c) {
    Fusion f;
    Fusion::Config cfg;
    // Literal operating point from atomS3R_ins_tfg::resetFusion_.
    cfg.sigma_a = V::Constant(0.12f);
    cfg.sigma_m = V::Constant(0.80f);
    cfg.gyro_noise_density = 0.00135f;
    cfg.with_mag = true;
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    cfg.gravity_magnitude = g;
    cfg.online_tune_warmup_sec = 10.0f;
    if (c.hold_bias) cfg.acc_bias_unlock_sec = 1000.0f;
    f.begin(cfg);
    f.setPeriodicAwCovSync(c.sync);
    if (c.fixed && !f.setFixedTuning(2.5f, 0.25f, 1.5f))
        throw std::runtime_error("fixed configuration rejected");

    const M Rbase = (Eigen::AngleAxisf(c.tilted ? 2.9670597f : 0.0f, V::UnitZ()) *
                     Eigen::AngleAxisf(c.tilted ? 0.4363323f : 0.0f, V::UnitY()) *
                     Eigen::AngleAxisf(c.tilted ? 1.0471976f : 0.0f, V::UnitX())).toRotationMatrix();
    // Normal test runs must not put diagnostics among simulator input CSVs.
    // A focused replay may opt in to traces in an existing directory.
    std::ofstream csv;
    if (const char* directory = std::getenv("OCEAN_IMU_TRACE_DIR")) {
        csv.open(std::string(directory)+"/tfg-device-"+c.name+".csv");
        if (!csv) throw std::runtime_error("cannot open device trace");
    }
    if (csv.is_open()) csv << "t,live,bias_enabled,pz,vz,Sz,awz,baz,tau,sigma,rs,tilt_error_deg\n";
    float peak_p=0, peak_v=0, peak_S=0, peak_aw=0, peak_ba=0, peak_tilt=0;
    float first_live=-1, first_bias=-1, mag_clock=0, log_clock=0;
    bool finite=true;
    const int n = static_cast<int>(std::lround(600.0 / double(c.dt)));
    for (int k=0; k<n; ++k) {
        const float t=static_cast<float>(double(k)*double(c.dt));
        constexpr float w=2.0f*pi*0.17f;
        const float roll=c.moving ? 0.12f*std::sin(w*t) : 0.0f;
        const float roll_rate=c.moving ? 0.12f*w*std::cos(w*t) : 0.0f;
        const M R = Rbase * Eigen::AngleAxisf(roll,V::UnitX()).toRotationMatrix();
        const V aw(0.0f,0.0f,c.moving ? 0.20f*std::sin(w*t) : 0.0f);
        V acc=R.transpose()*(aw+V(0,0,-g+c.residual));
        V gyro(roll_rate,0,0);
        if (c.noise) {
            acc += V(0.01f*std::sin(2*pi*7.3f*t),0.013f*std::cos(2*pi*11.1f*t),
                     0.015f*std::sin(2*pi*5.7f*t));
            gyro += V(0.0005f,-0.0003f,0.0002f);
        }
        f.update(c.dt,gyro,acc,35.0f);
        mag_clock += c.dt;
        if (mag_clock >= 0.04f) {
            mag_clock=0.0f;
            f.updateMag(R.transpose()*V(20,0,43));
        }
        const auto& core=f.mekf();
        if (f.isLive() && first_live<0) {
            first_live=t;
            std::cout << "EVENT case=" << c.name << " live_s=" << t
                      << " timeout=" << f.handoffTimedOut() << '\n';
        }
        if (f.isLive() && core.acc_bias_updates_enabled() && first_bias<0) {
            first_bias=t;
            std::cout << "EVENT case=" << c.name << " bias_unlock_s=" << t << '\n';
        }
        const V p=core.get_position(), v=core.get_velocity();
        const V S=core.get_integral_displacement(), a=core.get_world_accel(), ba=core.get_acc_bias();
        const bool ok=p.allFinite() && v.allFinite() && S.allFinite() && a.allFinite() &&
                      ba.allFinite() && core.covariance_full().allFinite();
        if (!ok) { finite=false; std::cerr << "NONFINITE case=" << c.name << " t=" << t << '\n'; break; }
        const float cosine=std::clamp((core.R_bw()*R.transpose())(2,2),-1.0f,1.0f);
        const float tilt=std::acos(cosine)*180.0f/pi;
        if (f.isLive()) {
            peak_p=std::max(peak_p,std::abs(p.z())); peak_v=std::max(peak_v,std::abs(v.z()));
            peak_S=std::max(peak_S,std::abs(S.z())); peak_aw=std::max(peak_aw,std::abs(a.z()));
            peak_ba=std::max(peak_ba,ba.norm()); peak_tilt=std::max(peak_tilt,tilt);
        }
        log_clock+=c.dt;
        if (log_clock>=1.0f && csv.is_open()) {
            log_clock=0;
            csv << std::setprecision(9) << t << ',' << f.isLive() << ',' << core.acc_bias_updates_enabled()
                << ',' << p.z() << ',' << v.z() << ',' << S.z() << ',' << a.z() << ',' << ba.z()
                << ',' << f.getTauApplied() << ',' << f.getSigmaApplied() << ',' << f.getRSFilterInput()
                << ',' << tilt << '\n';
        }
    }
    std::cout << std::setprecision(9) << "RESULT case=" << c.name << " residual=" << c.residual
              << " dt=" << c.dt << " live=" << first_live << " unlock=" << first_bias
              << " maxp=" << peak_p << " maxv=" << peak_v << " maxS=" << peak_S
              << " maxaw=" << peak_aw << " maxba=" << peak_ba << " maxtilt=" << peak_tilt
              << " pz=" << f.get_position().z() << " tau=" << f.getTauApplied()
              << " sigma=" << f.getSigmaApplied() << " rs=" << f.getRSFilterInput()
              << " finite=" << finite << '\n';
    // These bounded, small-disturbance histories have peaks below 1.20 m.
    // This is a regression bound, not output saturation in the estimator.
    return finite && first_live>=0 && peak_p<=2.0f;
}
} // namespace
int main() {
    std::cout << std::unitbuf;
    const Case cases[] = {
        {"level-zero",0.0f,0.005f},
        {"level-001",0.01f,0.005f},
        {"level-003",0.03f,0.005f},
        {"level-005",0.05f,0.005f},
        {"tilted-005",0.05f,0.005f,true},
        {"tilted-held",0.05f,0.005f,true,true},
        {"tilted-fixed",0.05f,0.005f,true,false,true},
        {"tilted-no-sync",0.05f,0.005f,true,false,false,false},
        {"tilted-50hz",0.05f,0.020f,true},
        {"tilted-20hz",0.05f,0.050f,true},
        {"tilted-10hz",0.05f,0.100f,true},
        {"tilted-noisy",0.05f,0.005f,true,false,false,true,false,true},
        {"moving-noisy",0.05f,0.005f,true,false,false,true,true,true},
    };
    bool ok=true;
    try { for (const auto& c: cases) ok=run(c) && ok; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
    return ok ? 0 : 1;
}
