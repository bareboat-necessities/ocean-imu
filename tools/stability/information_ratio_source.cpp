// Native driver for the non-promoting information-ratio feasibility diagnostic.
// The Python module inserts binary read-only taps into a temporary copy of the
// shipping header and compiles an untapped control; both must print the
// identical terminal state. Only begin, update and updateMag are used after
// construction. The finite float replay is not an all-time physical,
// magnetic-service or real-arithmetic certificate.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#ifdef INFO_RATIO_TAP
// Binary event stream: one tag byte, then a fixed number of doubles.
//   P prediction  F_AA(36) F_LL(144) Q_AA(36) Q_LL(144) phi_BA(1) Q_BA(9)
//   Y sync        Delta_AW(9)
//   A/M/S         correction H(63) R_eff(9)
//   G reset       d(3)
//   K marker      time(1) P(441) x(21) q(4)
static FILE* trace_file = nullptr;
static bool recording = false;
static void trace_put(double v) { std::fwrite(&v, sizeof v, 1, trace_file); }
template<class A> static void trace_mat(const A& a) {
    for (int i=0; i<a.rows(); ++i)
        for (int j=0; j<a.cols(); ++j) trace_put(static_cast<double>(a(i,j)));
}
static void trace_tag(char c) { std::fputc(c, trace_file); }
template<class A, class B, class C, class D, class E>
static void readout_prediction(const A& f, const B& fl, const C& q, const D& ql, double phi, const E& qb) {
    if (!recording) return;
    trace_tag('P'); trace_mat(f); trace_mat(fl); trace_mat(q); trace_mat(ql); trace_put(phi); trace_mat(qb);
}
template<class A> static void readout_sync(const A& q) {
    if (recording) { trace_tag('Y'); trace_mat(q); }
}
template<class A, class B>
static void readout_correction(char sensor, const A& h, const B& r) {
    if (recording) { trace_tag(sensor); trace_mat(h); trace_mat(r); }
}
template<class A> static void readout_reset(const A& d) {
    if (recording) { trace_tag('G'); trace_mat(d); }
}
#endif
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

// Even unit triangle wave of period one and its zero-mean integrals (as in
// aw_tracking_source.cpp): tri=1-4x on [0,1/2]; u'=tri, w'=u.
static double frac(double x) { return x-std::floor(x); }
static double tri(double x) { x=frac(x); return x<.5 ? 1-4*x : 4*x-3; }
static double tri_u(double x) { x=frac(x); return x<.5 ? x-2*x*x : 2*x*x-3*x+1; }
static double tri_w(double x) {
    x=frac(x);
    const double w = x<.5 ? x*x/2-2*x*x*x/3
                          : 1.0/24+(2*x*x*x/3-1.5*x*x+x)-(2.0/24-.375+.5);
    return w-1.0/48;
}
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
    // profile record_start_step record_stop_step marker_every trace_path
    if (argc<6) return 64;
    const std::string profile = argv[1];
    const int start = std::atoi(argv[2]), stop = std::atoi(argv[3]), every = std::atoi(argv[4]);
#ifdef INFO_RATIO_TAP
    trace_file = std::fopen(argv[5], "wb");
    if (!trace_file) return 65;
#endif
    const bool quiet = profile=="quiet", wave = profile=="wave";
    const bool collinear = profile=="collinear-service", locked = profile=="sync-locked";
    if (!(quiet || wave || collinear || locked)) return 66;
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.12f);
    cfg.sigma_g.setConstant(.00135f);
    cfg.sigma_m.setConstant(.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    const double G=9.80665, pi=3.14159265358979323846;
    // Collinear witness: a(t)=c cos(2 pi t), c=g(e_z-(e_z.b)b), b=(7,0,24)/25.
    const double cx = -168.0/625.0*G, cz = 49.0/625.0*G;
    int live=-1, refined=-1, active=-1, applied=0;
    for (int k=1; k<=stop; ++k) {
        const double t = static_cast<double>(k)*.005;
#ifdef INFO_RATIO_TAP
        recording = k>=start;
        if (recording && (k-start)%every==0) {
            // Marker at the root immediately before the prediction of step k.
            const auto& m=filter.raw().mekf();
            trace_tag('K'); trace_put(t-.005); trace_mat(m.Pext); trace_mat(m.xext); trace_mat(m.qref.coeffs());
        }
#endif
        if (k==start && (active<0 || k-active<3400)) return 2;
        float roll=0.0f, rate=0.0f;
        Eigen::Vector3d force(0.0,0.0,-G), field(75.0,0.0,0.0);
        if (wave) {
            roll = static_cast<float>(.02*std::sin(.5*t));
            rate = static_cast<float>(.01*std::cos(.5*t));
            force.z() = -.144*std::sin(.6*t)-G;
            field = Eigen::Vector3d(60.0,0.0,30.0);
        } else if (collinear) {
            roll = static_cast<float>(.01*std::sin(.5*t));
            rate = static_cast<float>(.005*std::cos(.5*t));
            const double c = std::cos(2.0*pi*t);
            force = Eigen::Vector3d(cx*c,0.0,cz*c-G);
            field = Eigen::Vector3d(21.0,0.0,72.0);
        } else if (locked) {
            // aw_tracking sync-locked-rectification profile (C2 onset 160--180 s).
            const double a1=2.6, f1=9.523809523809524, ph=.4;
            double e[4]; envelope(t,160.0,20.0,e);
            const double disp = a1/(f1*f1)*tri_w(f1*t+ph);
            const double vel = a1/f1*tri_u(f1*t+ph);
            const double acc = a1*tri(f1*t+ph);
            roll = static_cast<float>(.01*std::sin(.5*t));
            rate = static_cast<float>(.005*std::cos(.5*t));
            force = Eigen::Vector3d(e[0]*acc+2*e[1]*vel+e[2]*disp,0.0,-G);
            field = Eigen::Vector3d(21.0,0.0,72.0);
        }
        const Eigen::Matrix3d rwb=Eigen::AngleAxisd(-roll,Eigen::Vector3d::UnitX()).toRotationMatrix();
        filter.update(.005f,Eigen::Vector3f(rate,0,0),(rwb*force).cast<float>());
        if (k%8==0) {
            filter.updateMag((rwb*field).cast<float>());
            if (k>=start && filter.raw().mekf().lastMagDiag().accepted) ++applied;
        }
        if (live<0 && filter.isLive()) live=k;
        if (refined<0 && filter.hasRefinedMagReference()) refined=k;
        if (active<0 && filter.raw().mekf().acc_bias_updates_enabled()) active=k;
    }
    const auto& m=filter.raw().mekf();
#ifdef INFO_RATIO_TAP
    trace_tag('K'); trace_put(stop*.005); trace_mat(m.Pext); trace_mat(m.xext); trace_mat(m.qref.coeffs());
    std::fclose(trace_file);
#endif
    std::printf("END %d %d %d %d", live, refined, active, applied);
    for (int i=0; i<m.xext.size(); ++i) std::printf(" %.9g", static_cast<double>(m.xext(i)));
    for (int i=0; i<4; ++i) std::printf(" %.9g", static_cast<double>(m.qref.coeffs()(i)));
    std::printf(" %.9g\n", static_cast<double>(m.Pext.trace()));
    return m.Pext.allFinite() && m.xext.allFinite() ? 0 : 3;
}
