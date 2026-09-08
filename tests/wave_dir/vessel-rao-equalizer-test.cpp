#define EIGEN_NON_ARDUINO
#include "wave_dir/VesselRaoEqualizer.h"
#include "wave_dir/VesselRaoNoiseWeighting.h"
#include <complex>
#include <iostream>

int main() {
    // Weighting must preserve strong signals, remain bounded near zero, and
    // follow physical channel assignments independently of direction on/off.
    wave_direction::VesselRaoNoiseWeighting noise;
    const Eigen::Vector3f nominal = Eigen::Vector3f::Constant(.02f);
    if (noise.scales(0, .3f, nominal) != Eigen::Vector3f::Ones()) return 4;
    noise.max_std_scale = 4;
    if (noise.scales(0, .3f, nominal) != Eigen::Vector3f(4,4,1)) return 5;
    if (noise.scales(100, .3f, nominal) != Eigen::Vector3f::Ones()) return 6;
    const auto original = noise.scales(.15f, .3f, nominal);
    std::swap(noise.response.horizontal_x_tau_s, noise.response.horizontal_y_tau_s);
    noise.response.enabled = true;
    const auto swapped = noise.scales(.15f, .3f, nominal);
    if (original.x() != swapped.y() || original.y() != swapped.x()) return 7;
    float previous = 4;
    for (int j=0; j<1000; ++j) {
        const auto scales = noise.scales(float(j)*.001f, .3f, nominal);
        if (!scales.allFinite() || scales.minCoeff()<1 || scales.maxCoeff()>4 ||
            scales.x()>previous || scales.z()!=1) return 8;
        previous=scales.x();
    }
    if (noise.scales(NAN, .3f, nominal) != Eigen::Vector3f::Ones()) return 9;
    using C=std::complex<double>;
    const double pi=std::acos(-1.0),dt=.005;
    wave_direction::VesselRaoEqualizer e;
    wave_direction::HeadingFrameAcceleration<float> in{.3f,-.2f,.1f,true};
    auto out=e.step(in,float(dt));
    if(out.forward_ms2!=in.forward_ms2 || out.starboard_ms2!=in.starboard_ms2 || out.up_ms2!=in.up_ms2)return 1;
    double worst=0;
    for(double hz : {.05,.1,.3,.5,.8}) for(double degrees : {-130.,-30.,15.,70.,160.}) {
        wave_direction::VesselRaoEqualizer::Config cfg;cfg.enabled=true;e.configure(cfg);
        const double w=2*pi*hz,angle=degrees*pi/180,r=hz*cfg.heave_period_s;
        const C i(0,1),qx=1.0/std::pow(1.0-i*w*.7,2),qy=1.0/std::pow(1.0-i*w,2);
        const C gz=1.0/C(1-r*r,-2*.45*r);
        C measured_x=0,measured_y=0,measured_z=0;
        const int n=int(200/hz/dt),start=n/2;
        for(int j=0;j<n;++j) {
            const C phase=std::exp(-i*w*(j*dt));
            in={float(std::real(i*std::cos(angle)*qx*phase)),
                float(std::real(i*std::sin(angle)*qy*phase)),float(std::real(gz*phase)),true};
            out=e.step(in,float(dt));
            if(j>=start){measured_x+=double(out.forward_ms2)*std::conj(phase);measured_y+=double(out.starboard_ms2)*std::conj(phase);measured_z+=double(out.up_ms2)*std::conj(phase);}
        }
        const double scale=2.0/(n-start);measured_x*=scale;measured_y*=scale;measured_z*=scale;
        worst=std::max({worst,std::abs(measured_x-i*std::cos(angle)*qy),std::abs(measured_y-i*std::sin(angle)*qy),std::abs(measured_z-qy)});
        // Common horizontal/vertical response restores orbital quadrature,
        // including above heave resonance, without changing sense thresholds.
        const C along=std::cos(angle)*measured_x+std::sin(angle)*measured_y;
        if(std::imag(along*std::conj(measured_z))<=0)return 2;
    }
    std::cout<<"Max complex response error: "<<worst<<'\n';
    return worst<8e-5 ? 0:3;
}
