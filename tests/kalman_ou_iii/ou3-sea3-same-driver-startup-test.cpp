// Canonical non-promoting complete-SEA3 same-driver startup bridge.
//
// The mathematical source is tools/stability/ou3_sea3_fixed_history_source_core.py:
// one continuum Hilbert-ball coefficient field, one admissible SEA3 partition,
// one admissible continuum RAO member, no replay, no finite harmonic source,
// and no phase reseed. Simpson nodes below evaluate that continuum integral;
// they are not source modes. Keep every constant and transfer factor aligned
// with the Python source-core definition.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;

namespace {
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double H = 1.5;
constexpr double TP = 6.0;
constexpr double GAMMA = 3.3;
constexpr double FMIN = (1.0 / TP) / 64.0;
constexpr double FMAX = (1.0 / TP) * 256.0;
constexpr double RAO_G = 1.0;
constexpr double RAO_FC = 0.5;
constexpr double DRIVER_CENTER = 1.0 / TP;
constexpr double DRIVER_SIGMA = 0.005;
constexpr double DRIVER_BETA = 0.5;
constexpr double DT = 0.005;
constexpr double PREHISTORY_S = 60.0;
constexpr int WORD_SAMPLES = 601;
constexpr int COARSE_PANELS = 1024;
constexpr int FINE_PANELS = 2048;

struct SpectrumNorm { double scale{}; };

double unscaled_jonswap(double f) {
    const double fp = 1.0 / TP;
    const double x = f / fp;
    const double sigma = (x <= 1.0) ? 0.07 : 0.09;
    const double peak = std::exp(-((x - 1.0) * (x - 1.0)) /
                                 (2.0 * sigma * sigma));
    return std::pow(f, -5.0) * std::exp(-1.25 * std::pow(x, -4.0))
         * std::pow(GAMMA, peak);
}

double simpson(const std::vector<double>& y, double h) {
    if (y.size() < 3 || (y.size() % 2) == 0) std::abort();
    double s = y.front() + y.back();
    for (std::size_t i = 1; i + 1 < y.size(); ++i)
        s += (i % 2 ? 4.0 : 2.0) * y[i];
    return h * s / 3.0;
}

SpectrumNorm spectrum_norm(int panels) {
    if ((panels % 2) != 0) std::abort();
    const double u0 = std::log(FMIN), u1 = std::log(FMAX), du = (u1-u0)/panels;
    std::vector<double> vals(static_cast<std::size_t>(panels)+1);
    for (int i=0;i<=panels;++i) {
        const double f = std::exp(u0 + i*du);
        vals[static_cast<std::size_t>(i)] = unscaled_jonswap(f)*f;
    }
    const double integral = simpson(vals,du);
    return { (H*H/16.0)/integral };
}

double driver_norm_c(int panels) {
    const double u0 = std::log(FMIN), u1 = std::log(FMAX), du = (u1-u0)/panels;
    std::vector<double> vals(static_cast<std::size_t>(panels)+1);
    for (int i=0;i<=panels;++i) {
        const double f = std::exp(u0+i*du);
        vals[static_cast<std::size_t>(i)] =
            std::exp(-std::pow((f-DRIVER_CENTER)/DRIVER_SIGMA,2))*f;
    }
    const double norm2=simpson(vals,du);
    return 1.0/std::sqrt(norm2);
}

double source_accel(double t, int panels, double scale, double c) {
    const double u0 = std::log(FMIN), u1 = std::log(FMAX), du = (u1-u0)/panels;
    std::vector<double> vals(static_cast<std::size_t>(panels)+1);
    for (int i=0;i<=panels;++i) {
        const double f = std::exp(u0+i*du);
        const double S = scale*unscaled_jonswap(f);
        const double h = RAO_G*std::min(1.0,std::pow(RAO_FC/f,2.0));
        const double omega = 2.0*PI*f;
        const double acceleration_transfer = -(omega*omega)*h;
        const double a = DRIVER_BETA*c*
            std::exp(-0.5*std::pow((f-DRIVER_CENTER)/DRIVER_SIGMA,2));
        vals[static_cast<std::size_t>(i)] =
            acceleration_transfer*std::sqrt(S)*a*std::cos(omega*t)*f;
    }
    return simpson(vals,du);
}

Eigen::Vector3f specific_force_from_source(double a_z) {
    return Eigen::Vector3f(0.0f,0.0f,static_cast<float>(-g_std+a_z));
}

bool finite_positive(float x) { return std::isfinite(x) && x > 0.0f; }
} // namespace

int main() {
    const auto coarse_norm=spectrum_norm(COARSE_PANELS);
    const auto fine_norm=spectrum_norm(FINE_PANELS);
    const double coarse_c=driver_norm_c(COARSE_PANELS);
    const double fine_c=driver_norm_c(FINE_PANELS);

    double peak=0.0, max_delta=0.0;
    for(int k=0;k<WORD_SAMPLES;++k){
        const double t=PREHISTORY_S+k*DT;
        const double yc=source_accel(t,COARSE_PANELS,coarse_norm.scale,coarse_c);
        const double yf=source_accel(t,FINE_PANELS,fine_norm.scale,fine_c);
        peak=std::max(peak,std::abs(yf));
        max_delta=std::max(max_delta,std::abs(yf-yc));
    }
    const double quadrature_relative=max_delta/std::max(peak,1e-15);
    if (!(quadrature_relative < 2.0e-4)) {
        std::cerr << "continuum evaluation did not converge tightly enough\n";
        return 1;
    }

    using Filter=SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
    Filter f(false);
    f.setWithMag(false);
    f.setOnlineTuneWarmupSec(10.0f);
    f.initialize(Eigen::Vector3f::Constant(0.2f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.3f));
    const Eigen::Vector3f gyro=Eigen::Vector3f::Zero();
    f.initialize_from_acc(specific_force_from_source(
        source_accel(0.0,FINE_PANELS,fine_norm.scale,fine_c)));

    int first_ready=-1;
    double max_abs=0.0;
    float max_guard=0.0f;
    const int pre_samples=static_cast<int>(std::llround(PREHISTORY_S/DT));
    for(int k=0;k<pre_samples;++k){
        const double t=k*DT;
        const double y=source_accel(t,FINE_PANELS,fine_norm.scale,fine_c);
        max_abs=std::max(max_abs,std::abs(y));
        f.updateFrontEnd(static_cast<float>(DT),gyro,specific_force_from_source(y));
        max_guard=std::max(max_guard,f.accelVibrationGuardEngagement());
        if(first_ready<0 && f.isTunerReady()) first_ready=k;
    }
    if(first_ready<0 || !f.isTunerReady() || !f.wavePeriodUsable() || !f.startupProxyInitialized()) {
        std::cerr << "same continuum SEA3 history did not reach real shipping TunerReady before Live entry\n";
        return 1;
    }
    if(max_guard!=0.0f || max_abs>4.0) {
        std::cerr << "legal point left declared dormant-guard/acceleration branch\n";
        return 1;
    }

    const float tau_live=f.getTauApplied();
    const float sigma_live=f.getSigmaApplied();
    const float rs_live=f.getRSApplied();
    const float period_live=f.getWavePeriodSec();
    if(!finite_positive(tau_live)||!finite_positive(sigma_live)||
       !finite_positive(rs_live)||!finite_positive(period_live)){
        std::cerr << "shipping tuner schedule invalid at Live entry\n";
        return 1;
    }

    // Real shipping handoff at the declared t=60 s boundary. This is the H
    // (ungauged-yaw) startup mode; A21 release remains a separate hybrid event
    // in the complete word.
    f.goLive(f.startupProxyQuat(),0.035f,1.5708f);
    if(!f.isAdaptiveLive()) {
        std::cerr << "real goLive handoff failed\n";
        return 1;
    }
    const auto P0=f.mekf().covariance_full();
    const auto ldlt=P0.selfadjointView<Eigen::Lower>().ldlt();
    if(!P0.allFinite() || ldlt.info()!=Eigen::Success || ldlt.vectorD().minCoeff()<=0.0f){
        std::cerr << "shipping Live covariance seed not SPD\n";
        return 1;
    }

    float rs_min=std::numeric_limits<float>::infinity();
    float rs_max=0.0f;
    for(int k=0;k<WORD_SAMPLES;++k){
        const double t=PREHISTORY_S+k*DT;
        const double y=source_accel(t,FINE_PANELS,fine_norm.scale,fine_c);
        max_abs=std::max(max_abs,std::abs(y));
        f.updateTime(static_cast<float>(DT),gyro,specific_force_from_source(y),35.0f);
        const float rs=f.getRSApplied();
        rs_min=std::min(rs_min,rs);
        rs_max=std::max(rs_max,rs);
    }
    if(!finite_positive(rs_min)||!finite_positive(rs_max)||max_abs>4.0){
        std::cerr << "same-history Live continuation invalid\n";
        return 1;
    }

    std::cout<<std::setprecision(10)
             <<"COMPLETE_SEA3_SAME_DRIVER_STARTUP_PASS"
             <<" first_ready_s="<<(first_ready+1)*DT
             <<" wave_period_s="<<period_live
             <<" tau_s="<<tau_live
             <<" sigma_mps2="<<sigma_live
             <<" RS_entry="<<rs_live
             <<" RS_word_min="<<rs_min
             <<" RS_word_max="<<rs_max
             <<" max_abs_accel="<<max_abs
             <<" quadrature_rel="<<quadrature_relative
             <<"\n";
    return 0;
}