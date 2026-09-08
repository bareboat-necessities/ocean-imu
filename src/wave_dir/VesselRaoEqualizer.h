#pragma once

#include "wave_dir/WaveDirectionFrame.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace wave_direction {

// Match all channels to a common low-pass response, rather than invert the
// vessel response to unity. Applies only to the direction branch. Translation
// footprint/depth attenuation remains; this does not reconstruct incident Hs.
class VesselRaoEqualizer {
public:
    struct Config {
        bool enabled = false;
        double horizontal_x_tau_s = 0.7;
        double horizontal_y_tau_s = 1.0;
        double heave_period_s = 2.4;
        double heave_damping = 0.45;
    };

    void configure(const Config& cfg) {
        for (double x : {cfg.horizontal_x_tau_s, cfg.horizontal_y_tau_s,
                         cfg.heave_period_s, cfg.heave_damping}) {
            if (!(std::isfinite(x) && x > 0))
                throw std::invalid_argument("Invalid direction RAO parameters");
        }
        config_ = cfg;
        reset();
    }
    void reset() { filters_ = {}; last_dt_ = 0; }
    const Config& config() const { return config_; }

    HeadingFrameAcceleration<float> step(
        const HeadingFrameAcceleration<float>& input, float dt) {
        if (!config_.enabled) return input;
        if (!(dt > 0 && std::isfinite(dt)) || !input.heading_valid ||
            !std::isfinite(input.forward_ms2) || !std::isfinite(input.starboard_ms2) ||
            !std::isfinite(input.up_ms2)) {
            reset();
            return {};
        }
        if (double(dt) != last_dt_) coefficients(double(dt));
        auto out = input;
        out.forward_ms2 = float(filters_[0].step(input.forward_ms2));
        out.starboard_ms2 = float(filters_[1].step(input.starboard_ms2));
        out.up_ms2 = float(filters_[2].step(input.up_ms2));
        return out;
    }

private:
    struct Biquad {
        double b0=1, b1=0, b2=0, a1=0, a2=0, z1=0, z2=0;
        bool initialized=false, identity=true;
        double step(double x) {
            if (identity) return x;
            if (!initialized) {
                // Unity DC response: initialize at a constant input, without
                // an artificial impulse when the direction branch starts.
                z1=(1-b0)*x; z2=(b2-a2)*x; initialized=true;
            }
            const double y=b0*x+z1;
            z1=b1*x-a1*y+z2; z2=b2*x-a2*y;
            return y;
        }
    };
    void coefficients(double dt) {
        const double wn=2*std::acos(-1.0)/config_.heave_period_s;
        // This target guarantees |E_j(iw)| <= 1 for each of the three
        // compensators at every frequency (coefficient-wise quadratic bound).
        const double th=std::max(1.0,std::sqrt(std::max(0.0,
                                   2*config_.heave_damping*config_.heave_damping-1)))/wn;
        const double t=std::max({config_.horizontal_x_tau_s,
                                 config_.horizontal_y_tau_s,th});
        const double k=2/dt, den=1+2*t*k+t*t*k*k;
        auto set=[&](int index,double n1,double n2,bool identity) {
            auto& f=filters_[std::size_t(index)];
            f.b0=(1+n1*k+n2*k*k)/den;
            f.b1=(2-2*n2*k*k)/den;
            f.b2=(1-n1*k+n2*k*k)/den;
            f.a1=(2-2*t*t*k*k)/den;
            f.a2=(1-2*t*k+t*t*k*k)/den;
            f.identity=identity;
        };
        const double tx=config_.horizontal_x_tau_s, ty=config_.horizontal_y_tau_s;
        set(0,2*tx,tx*tx,tx==t);
        set(1,2*ty,ty*ty,ty==t);
        set(2,2*config_.heave_damping/wn,1/(wn*wn),false);
        last_dt_=dt;
    }
    Config config_{};
    std::array<Biquad,3> filters_{};
    double last_dt_=0;
};

} // namespace wave_direction
