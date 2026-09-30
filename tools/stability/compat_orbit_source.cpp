// Native driver for the non-promoting compatibility-orbit feasibility
// diagnostic (compat_orbit_source_diagnostic.py). It runs the unchanged
// shipping SeaStateFusion_OU_III on one supplied physical history and prints
// the carried nominal AW, BA, S and attitude after every IMU sample. Only
// begin, update and updateMag are used after construction; nothing is written
// into the filter. The float replay is one execution, not an enclosure.
//
// argv: input file, roll amplitude [rad], roll frequency [Hz], pitch
// amplitude [rad], pitch frequency [Hz], field north/east/down [uT],
// magnetic decimation (samples), record start [s], constant body
// accelerometer bias x/y/z [m/s^2].
// The input file holds one line "ax ay az" (world NED physical acceleration,
// m/s^2) per 5-ms IMU sample.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
const float g_std = 9.80665f;
// Literal state offsets of Kalman3D_Wave_OU_III<float,true,true> (BASE_N=6).
constexpr int OFF_S = 12, OFF_AW = 15, OFF_BA = 18;

int main(int argc, char** argv) {
    if (argc < 14) return 64;
    std::ifstream in(argv[1]);
    std::vector<Eigen::Vector3d> acc;
    for (double x, y, z; in >> x >> y >> z;) acc.emplace_back(x, y, z);
    const double ra=atof(argv[2]), rf=atof(argv[3]), pa=atof(argv[4]), pf=atof(argv[5]);
    const Eigen::Vector3d B(atof(argv[6]), atof(argv[7]), atof(argv[8]));
    const int mag_every = atoi(argv[9]);
    const double trec = atof(argv[10]);
    const Eigen::Vector3d bias(atof(argv[11]), atof(argv[12]), atof(argv[13]));
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(.12f);
    cfg.sigma_g.setConstant(.00135f);
    cfg.sigma_m.setConstant(.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    const double G=9.80665, pi=3.14159265358979323846, h=.005;
    int active = -1;
    for (std::size_t k=1; k<=acc.size(); ++k) {
        const double t = k*h;
        // Body attitude: roll about x then pitch about y (world->body is
        // Rx(-roll) Ry(-pitch)); body rates from the exact derivative.
        const double wr=2*pi*rf, wp=2*pi*pf;
        const double roll=ra*std::sin(wr*t), droll=ra*wr*std::cos(wr*t);
        const double pitch=pa*std::sin(wp*t), dpitch=pa*wp*std::cos(wp*t);
        const Eigen::Matrix3d Rx=Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()).toRotationMatrix();
        const Eigen::Matrix3d Ry=Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()).toRotationMatrix();
        const Eigen::Matrix3d Rbw=Ry*Rx;                 // body->world
        const Eigen::Vector3d rate = Rx.transpose()*Eigen::Vector3d(0, dpitch, 0)
                                   + Eigen::Vector3d(droll, 0, 0);
        const Eigen::Vector3d force = acc[k-1]-Eigen::Vector3d(0, 0, G);
        filter.update(float(h), rate.cast<float>(), (Rbw.transpose()*force+bias).cast<float>());
        if (mag_every>0 && k%mag_every==0) filter.updateMag((Rbw.transpose()*B).cast<float>());
        const auto& m = filter.raw().mekf();
        if (active<0 && m.acc_bias_updates_enabled()) active = int(k);
        if (t >= trec) {
            // Attitude error E = R_true R_hat' (world frame), as rotation vector.
            const Eigen::Matrix3d Rhat = m.R_wb().template cast<double>().transpose();
            const Eigen::AngleAxisd e(Rbw*Rhat.transpose());
            const Eigen::Vector3d ev = e.angle()*e.axis();
            std::printf("%zu", k);
            for (int i=0; i<3; ++i) std::printf(" %.9g", double(m.xext(OFF_AW+i)));
            for (int i=0; i<3; ++i) std::printf(" %.9g", double(m.xext(OFF_BA+i)));
            for (int i=0; i<3; ++i) std::printf(" %.9g", double(m.xext(OFF_S+i)));
            std::printf(" %.9g %.9g %.9g %.9g %.9g\n", ev.x(), ev.y(), ev.z(),
                        double(m.tau_aw), double(m.get_pseudo_update_period_s()));
        }
    }
    const auto& m = filter.raw().mekf();
    std::printf("END %d", active);
    for (int i=0; i<m.xext.size(); ++i) std::printf(" %.9g", double(m.xext(i)));
    std::printf(" %.9g\n", double(m.Pext.trace()));
    return m.Pext.allFinite() && m.xext.allFinite() ? 0 : 3;
}
