#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

namespace {
using Filter = Kalman3D_Wave_OU_III<double, true, true>;
using Vec3 = Eigen::Vector3d;
constexpr int kNx = 21;
constexpr int kBg = 3;
constexpr int kV = 6;
constexpr int kP = 9;
constexpr int kS = 12;
constexpr int kAw = 15;
constexpr int kBa = 18;

struct Point {
    double dt;
    double tau;
    double sigma_aw;
    double sigma_s;
    double roll;
    double pitch;
    Vec3 omega;
};

Eigen::Quaterniond body_to_world(double roll, double pitch, double yaw = 0.0) {
    return Eigen::AngleAxisd(yaw, Vec3::UnitZ())
         * Eigen::AngleAxisd(pitch, Vec3::UnitY())
         * Eigen::AngleAxisd(roll, Vec3::UnitX());
}

Filter make_filter(const Point& p) {
    const Vec3 sigma_acc = Vec3::Constant(0.0148);
    const Vec3 sigma_gyro = Vec3::Constant(0.00157);
    const Vec3 sigma_mag = Vec3::Constant(0.25);
    Filter f(sigma_acc, sigma_gyro, sigma_mag);
    f.set_aw_time_constant(p.tau);
    f.set_aw_stationary_std(Vec3::Constant(p.sigma_aw));
    f.set_RS_noise(Vec3::Constant(p.sigma_s));
    f.set_Racc_std(sigma_acc);
    f.set_Rmag(sigma_mag);
    f.set_initial_linear_uncertainty(0.25, 0.25, 0.5);
    f.set_initial_acc_bias_std(0.02);
    f.set_Q_bacc_rw(Vec3::Constant(5.0e-4));
    f.set_mag_world_ref(Filter::ned_field_from_decl_incl(0.0, 1.05, 48.0));
    const Eigen::Quaterniond qbw = body_to_world(p.roll, p.pitch);
    f.qref = qbw.conjugate();
    f.qref.normalize();
    f.xext.setZero();
    f.Pext = 0.5 * (f.Pext + f.Pext.transpose());
    return f;
}

Vec3 predicted_acc(const Filter& f) {
    const Vec3 gravity_world(0.0, 0.0, 9.80665);
    return f.qref * (gravity_world + f.xext.segment<3>(kAw))
         + f.xext.segment<3>(kBa);
}

Vec3 predicted_mag(const Filter& f) {
    return f.qref * f.v2ref;
}

void one_cycle(Filter& f, const Point& p) {
    const Vec3 acc = predicted_acc(f);
    const Vec3 mag = predicted_mag(f);
    f.time_update(p.omega, p.dt);
    f.measurement_update_acc_only(acc, 35.0);
    f.measurement_update_mag_only(mag);
    f.applyIntegralZeroPseudoMeas();
}

void settle_covariance(Filter& f, const Point& p) {
    for (int i = 0; i < 240; ++i) one_cycle(f, p);
    f.xext.setZero();
}

void inject(Filter& f, int j, double h) {
    if (j < 3) {
        Vec3 d = Vec3::Zero();
        d(j) = h;
        Eigen::Quaterniond dq(1.0, 0.5*d.x(), 0.5*d.y(), 0.5*d.z());
        dq.normalize();
        f.qref = dq * f.qref;
        f.qref.normalize();
    } else {
        f.xext(j) += h;
    }
}

Eigen::Matrix<double, kNx, 1> error_between(const Filter& pert, const Filter& nominal) {
    Eigen::Matrix<double, kNx, 1> e;
    e.setZero();
    Eigen::Quaterniond dq = pert.qref * nominal.qref.conjugate();
    if (dq.w() < 0.0) dq.coeffs() *= -1.0;
    e.segment<3>(0) = 2.0 * dq.vec();
    e.segment<3>(kBg) = pert.xext.segment<3>(kBg) - nominal.xext.segment<3>(kBg);
    e.segment<3>(kV) = pert.xext.segment<3>(kV) - nominal.xext.segment<3>(kV);
    e.segment<3>(kP) = pert.xext.segment<3>(kP) - nominal.xext.segment<3>(kP);
    e.segment<3>(kS) = pert.xext.segment<3>(kS) - nominal.xext.segment<3>(kS);
    e.segment<3>(kAw) = pert.xext.segment<3>(kAw) - nominal.xext.segment<3>(kAw);
    e.segment<3>(kBa) = pert.xext.segment<3>(kBa) - nominal.xext.segment<3>(kBa);
    return e;
}

Eigen::Matrix<double, kNx, kNx> jacobian(const Point& p) {
    Filter base = make_filter(p);
    settle_covariance(base, p);
    Filter nominal = base;
    one_cycle(nominal, p);

    Eigen::Matrix<double, kNx, kNx> F;
    const double h_att = 2.0e-6;
    const double h_add = 1.0e-6;
    for (int j = 0; j < kNx; ++j) {
        const double h = j < 3 ? h_att : h_add;
        Filter plus = base;
        Filter minus = base;
        inject(plus, j, h);
        inject(minus, j, -h);
        one_cycle(plus, p);
        one_cycle(minus, p);
        F.col(j) = (error_between(plus, nominal) - error_between(minus, nominal)) / (2.0*h);
    }
    return F;
}

std::vector<Point> envelope() {
    std::vector<Point> out;
    const std::vector<double> dts{1.0/240.0, 1.0/200.0};
    const std::vector<double> taus{0.55, 1.8, 4.0};
    const std::vector<double> sigmas{0.12, 0.9, 2.8};
    const std::vector<double> sigma_s{0.35, 1.5, 6.0};
    const std::vector<double> rolls{-0.70, 0.0, 0.70};
    const std::vector<double> pitches{-0.45, 0.0, 0.45};
    const std::vector<Vec3> rates{
        Vec3(-0.55,-0.55,-0.30), Vec3::Zero(), Vec3(0.55,0.55,0.30)
    };
    for (double dt : dts)
      for (double tau : taus)
       for (double sa : sigmas)
        for (double ss : sigma_s)
         for (double r : rolls)
          for (double q : pitches)
           for (const Vec3& w : rates)
             out.push_back(Point{dt,tau,sa,ss,r,q,w});
    return out;
}
}

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "iss_jacobians.txt";
    std::ofstream os(path);
    if (!os) {
        std::cerr << "cannot open " << path << "\n";
        return 2;
    }
    os << std::setprecision(17);
    const auto points = envelope();
    os << "OCEAN_IMU_ISS_JACOBIANS_V1 " << kNx << " " << points.size() << "\n";
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        const auto F = jacobian(p);
        if (!F.allFinite()) {
            std::cerr << "non-finite Jacobian at point " << i << "\n";
            return 3;
        }
        os << "POINT " << i << " " << p.dt << " " << p.tau << " "
           << p.sigma_aw << " " << p.sigma_s << " " << p.roll << " "
           << p.pitch << " " << p.omega.x() << " " << p.omega.y() << " "
           << p.omega.z() << "\n";
        for (int r = 0; r < kNx; ++r) {
            for (int c = 0; c < kNx; ++c) {
                if (c) os << ' ';
                os << F(r,c);
            }
            os << '\n';
        }
    }
    std::cout << "wrote " << points.size() << " implementation Jacobians to " << path << "\n";
    return 0;
}
