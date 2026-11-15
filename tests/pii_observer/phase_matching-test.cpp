// Independently check the scheduled observer's physical frequency response.
#include "pii_observer/VerticalPIIObserver.h"
#include <cmath>
#include <complex>
#include <iostream>

template<typename T> bool check_grid() {
    using Observer = marine_obs::VerticalPIIObserver<T, false>;
    for (const double f : {.03, .06, .12, .17, .25, .5, 1.0}) {
        for (const double sigma : {.001, .03, .1, .25, 1.1, 5.0}) {
            Observer old, corrected;
            typename Observer::AdaptConfig a;
            a.match_accel_filter_phase = false;
            old.configure_adaptation(a);
            a.match_accel_filter_phase = true;
            corrected.configure_adaptation(a);
            for (int i=0; i<1800; ++i) {
                old.update_adaptation(T(f), T(sigma), T(1), T(1));
                corrected.update_adaptation(T(f), T(sigma), T(1), T(1));
            }
            if (old.motion_pole_rate() != corrected.motion_pole_rate()) return false;
            const double tau0=old.accel_filter_tau(), tau1=corrected.accel_filter_tau();
            if (!std::isfinite(tau1) || tau1+1e-6<tau0 || tau1>a.tau_a_max+1e-6) return false;
            const std::complex<double> s(0, 6.283185307179586*f);
            const double r=corrected.motion_pole_rate();
            const auto core=std::pow(s/(s+r),3);
            const double e0=std::norm(core/(1.0+s*tau0)-1.0);
            const double e1=std::norm(core/(1.0+s*tau1)-1.0);
            if (e1>e0+2e-5) return false;
        }
    }
    return true;
}

int main() {
    if (!check_grid<float>() || !check_grid<double>()) {
        std::cerr << "PII phase/noise-smoothing guard failed\n";return 1;
    }
    using Observer=marine_obs::VerticalPIIObserver<double,false>;
    Observer old, corrected;
    Observer::Config c; c.r=.25;c.tau_a=.15;
    old.configure(c);corrected.configure(c);
    Observer::AdaptConfig a;
    a.r_freq_exp=a.r_sigma_exp=a.tau_a_freq_exp=a.tau_a_sigma_exp=0;
    a.match_accel_filter_phase=false;old.configure_adaptation(a);
    a.match_accel_filter_phase=true;corrected.configure_adaptation(a);
    constexpr double w=6.283185307179586*.17, dt=.005;
    for(int i=0;i<1800;++i) {
        old.update_adaptation(.17,1,1,1);corrected.update_adaptation(.17,1,1,1);
    }
    double e0=0,e1=0;
    for(int k=0;k<120000;++k) {
        const double truth=std::sin(w*k*dt), accel=-w*w*truth;
        const double p0=old.update(accel,dt),p1=corrected.update(accel,dt);
        if(k>=60000) {e0+=(p0-truth)*(p0-truth);e1+=(p1-truth)*(p1-truth);}
    }
    if (!(e1<.64*e0)) {std::cerr<<"PII long-wave reconstruction did not improve\n";return 1;}
    std::cout<<"PII phase response, bounds and long-wave reconstruction passed\n";
}
