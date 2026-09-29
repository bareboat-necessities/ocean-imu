// Native driver for the non-promoting AW tracking feasibility diagnostic.
// The Python module inserts one read-only tap (the pre-correction AW mean)
// into a temporary copy of the shipping header and compiles an untapped
// control; both must print the identical terminal state. Only begin, update
// and updateMag are used after construction. The finite float replay is not
// an all-time physical, magnetic-service or real-arithmetic certificate.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#ifdef AW_TRACKING_TAP
static Eigen::Vector3d trace_aw_pre = Eigen::Vector3d::Zero();
// Event stream (prediction P, reset R, accelerometer row A, applied
// magnetic row M) for the literal world-frame six-column array.
static bool aw_event_rec = false;
#endif
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

// Even unit triangle wave of period one and its zero-mean integrals:
// tri=1-4x on [0,1/2]; u'=tri, w'=u, |u|<=1/8, |w-1/48|<=1/48.
static double frac(double x) { return x-std::floor(x); }
static double tri(double x) { x=frac(x); return x<.5 ? 1-4*x : 4*x-3; }
static double tri_u(double x) { x=frac(x); return x<.5 ? x-2*x*x : 2*x*x-3*x+1; }
static double tri_w(double x) {
    x=frac(x);
    const double w = x<.5 ? x*x/2-2*x*x*x/3
                          : 1.0/24+(2*x*x*x/3-1.5*x*x+x)-(2.0/24-.375+.5);
    return w-1.0/48;
}
// Quintic smoothstep envelope on [t0,t0+Tr]: C2, derivatives bounded by
// 15/8/Tr, 10/sqrt(3)/Tr^2 and 60/Tr^3.
static void envelope(double t, double t0, double tr, double e[4]) {
    const double x = (t-t0)/tr;
    if (x<=0) { e[0]=e[1]=e[2]=e[3]=0; return; }
    if (x>=1) { e[0]=1; e[1]=e[2]=e[3]=0; return; }
    e[0]=x*x*x*(10-15*x+6*x*x);
    e[1]=30*x*x*(1-x)*(1-x)/tr;
    e[2]=60*x*(1-x)*(1-2*x)/(tr*tr);
    e[3]=60*(1-6*x+6*x*x)/(tr*tr*tr);
}

int main(int argc, char** argv) {
    if (argc<11) return 64;
    // Horizontal (along the field's horizontal axis e_x) triangle amplitude and
    // frequency, horizontal cosine amplitude and frequency, vertical cosine
    // amplitude and frequency, onset t0 and ramp Tr, end time, record start.
    const double A1=atof(argv[1]), f1=atof(argv[2]), A2=atof(argv[3]), f2=atof(argv[4]);
    const double A3=atof(argv[5]), f3=atof(argv[6]);
    const double t0=atof(argv[7]), tr=atof(argv[8]), tend=atof(argv[9]), trec=atof(argv[10]);
    // Optional phase (cycles) of the horizontal triangle.
    const double ph = argc>11 ? atof(argv[11]) : 0.0;
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.12f);
    cfg.sigma_g.setConstant(.00135f);
    cfg.sigma_m.setConstant(.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    const double G=9.80665, pi=3.14159265358979323846;
    const int n = static_cast<int>(std::lround(tend/.005));
    int active=-1;
    for (int k=1; k<=n; ++k) {
        const double t=k*.005;
#ifdef AW_TRACKING_TAP
        aw_event_rec = t >= trec;
#endif
        // Horizontal displacement p=env*d; a=env d''+2 env' d'+env'' d.
        double e[4]; envelope(t,t0,tr,e);
        const double w2=2*pi*f2;
        const double disp = A1/(f1*f1)*tri_w(f1*t+ph) - A2/(w2*w2)*std::cos(w2*t);
        const double vel = A1/f1*tri_u(f1*t+ph) + A2/w2*std::sin(w2*t);
        const double acc = A1*tri(f1*t+ph) + A2*std::cos(w2*t);
        const double ax = e[0]*acc+2*e[1]*vel+e[2]*disp;
        const double az = A3*std::cos(2*pi*f3*t);
        const double roll=.01*std::sin(.5*t), rate=.005*std::cos(.5*t);
        const Eigen::Vector3d force(ax,0.0,az-G);
        const Eigen::Matrix3d rwb=Eigen::AngleAxisd(-roll,Eigen::Vector3d::UnitX()).toRotationMatrix();
        filter.update(.005f,Eigen::Vector3f(static_cast<float>(rate),0,0),(rwb*force).cast<float>());
        if (k%8==0) filter.updateMag((rwb*Eigen::Vector3d(21,0,72)).cast<float>());
        const auto& m=filter.raw().mekf();
        if (active<0 && m.acc_bias_updates_enabled()) active=k;
#ifdef AW_TRACKING_TAP
        if (t>=trec) {
            const Eigen::Matrix3d E=rwb*m.R_wb().template cast<double>().transpose();
            const double c=std::min(1.0,std::max(-1.0,(E.trace()-1)/2));
            std::printf("%d %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n", k, ax, az,
                trace_aw_pre.x(), trace_aw_pre.y(), trace_aw_pre.z(), std::acos(c),
                static_cast<double>(m.tau_aw));
        }
#endif
    }
    const auto& m=filter.raw().mekf();
#ifdef AW_TRACKING_TAP
    std::printf("B %.17g %.17g %.17g\n", static_cast<double>(m.v2ref.x()),
                static_cast<double>(m.v2ref.y()), static_cast<double>(m.v2ref.z()));
#endif
    std::printf("END %d", active);
    for (int i=0; i<m.xext.size(); ++i) std::printf(" %.9g", static_cast<double>(m.xext(i)));
    for (int i=0; i<4; ++i) std::printf(" %.9g", static_cast<double>(m.qref.coeffs()(i)));
    std::printf(" %.9g\n", static_cast<double>(m.Pext.trace()));
    return m.Pext.allFinite() && m.xext.allFinite() ? 0 : 3;
}
