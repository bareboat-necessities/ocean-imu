#define EIGEN_NON_ARDUINO
#include "ahrs/KalmanQMEKF.h"
#include <cstdio>
#include <cstdlib>

using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
using Quat = Eigen::Quaterniond;
static int failures = 0;
/** @brief Return the quaternion exponential of a rotation vector in radians. */
static Quat expq(const Vec& v) {
    const double a = v.norm();
    return a == 0 ? Quat::Identity() : Quat(Eigen::AngleAxisd(a, v/a));
}
/** @brief Return the shortest rotation vector in radians for a normalized quaternion. */
static Vec logq(Quat q) {
    q.normalize();
    if (q.w() < 0) q.coeffs() *= -1;
    const double n = q.vec().norm();
    return n < 1e-14 ? Vec(2*q.vec()) : Vec(2*std::atan2(n,q.w())/n*q.vec());
}
/** @brief Return the skew matrix representing cross product with v. */
static Mat cross(const Vec& v) {
    Mat m; m << 0,-v.z(),v.y(), v.z(),0,-v.x(), -v.y(),v.x(),0; return m;
}
/** @brief Record a test failure unless error is strictly below tolerance. */
static void check(double error, double tolerance, const char* name) {
    if (!(error < tolerance)) {
        std::fprintf(stderr, "%s: error %.12g >= %.12g\n", name, error, tolerance);
        ++failures;
    }
}

/** @brief Compare predicted covariance with finite differences of the physical right-error flow. */
template<bool bias> static void prediction(const Vec& rate, double dt) {
    constexpr int n = bias ? 6 : 3;
    using Matrix = Eigen::Matrix<double,n,n>;
    QuaternionMEKF<double,bias> f(Vec::Ones(),Vec::Zero(),Vec::Ones(),.7,.4,0);
    const Quat start = expq(Vec(-.4,.2,.3));
    const Vec gravity_body = start.conjugate()*Vec(0,0,-1);
    f.initialize_from_acc_mag(gravity_body,start.conjugate()*Vec(1,0,0));
    // Zero innovation makes the covariance anisotropic without changing truth.
    f.measurement_update_acc_only(gravity_body);
    const Matrix initial = f.covariance();
    Matrix jacobian = Matrix::Zero();
    const Quat nominal = start*expq(Vec(rate*dt));
    // Differentiate the actual right-error flow, not a second implementation
    // of Rodrigues' formula. Physical bias error changes the true body rate.
    const auto error = [&](int axis, double h) -> Eigen::Matrix<double,n,1> {
        Vec angle = Vec::Zero(), db = Vec::Zero();
        if (axis < 3) angle(axis) = h; else db(axis-3) = h;
        const Quat truth = start*expq(angle)*expq(Vec((rate-db)*dt));
        Eigen::Matrix<double,n,1> e;
        e.template head<3>() = logq(nominal.conjugate()*truth);
        if constexpr (bias) e.template tail<3>() = db;
        return e;
    };
    for (int j=0; j<n; ++j) jacobian.col(j) = (error(j,1e-6)-error(j,-1e-6))/2e-6;
    const Matrix expected = jacobian*initial*jacobian.transpose();
    f.time_update(rate,dt);
    check((f.covariance()-expected).norm(), 2e-9, "right-error prediction covariance");
    check(Quat(f.quaternion()).angularDistance(nominal), 1e-12, "nominal prediction");
}

/** @brief Compare covariance chart transport with finite differences for joint or accel-only updates. */
static void reset(bool joint) {
    using Matrix = Eigen::Matrix<double,6,6>;
    QuaternionMEKF<double,true> f(Vec::Constant(.2),Vec::Zero(),Vec::Constant(.2),.7,.4,0);
    f.initialize_from_acc_mag(Vec(0,0,-1),Vec(1,0,0));
    // Build attitude/bias cross-covariance before a deliberately finite update.
    f.time_update(Vec::Zero(),.3);
    const Matrix prior = f.covariance();
    Matrix H = Matrix::Zero();
    H.topLeftCorner<3,3>() = cross(Vec(0,0,-1));
    if (joint) H.block<3,3>(3,0) = cross(Vec(1,0,0));
    const Matrix R = .04*Matrix::Identity();
    const Quat target = expq(Vec(.6,-.4,.3));
    const Vec acc = target.conjugate()*Vec(0,0,-1);
    const Vec mag = target.conjugate()*Vec(1,0,0);
    Eigen::Matrix<double,6,1> residual = Eigen::Matrix<double,6,1>::Zero();
    residual.head<3>() = acc-Vec(0,0,-1);
    if (joint) residual.tail<3>() = mag-Vec(1,0,0);
    const Matrix gain = prior*H.transpose()*(H*prior*H.transpose()+R).inverse();
    const Vec correction = (gain*residual).head<3>();
    const Matrix a = Matrix::Identity()-gain*H;
    const Matrix posterior = a*prior*a.transpose()+gain*R*gain.transpose();
    Matrix reset_jacobian = Matrix::Identity();
    for (int j=0; j<3; ++j) {
        const Vec delta = 1e-6*Vec::Unit(j);
        // New error: Log(Exp(-correction) Exp(correction + epsilon)).
        reset_jacobian.block<3,1>(0,j) =
            (logq(expq(-correction)*expq(Vec(correction+delta)))-
             logq(expq(-correction)*expq(Vec(correction-delta))))/2e-6;
    }
    if (joint) f.measurement_update(acc,mag); else f.measurement_update_acc_only(acc);
    check(Quat(f.quaternion()).angularDistance(expq(correction)),1e-12,"measurement injection");
    check((f.covariance()-reset_jacobian*posterior*reset_jacobian.transpose()).norm(),
          2e-9,"right-error covariance reset");
}

/** @brief Run qMEKF propagation, spin-axis covariance, and finite-injection geometry checks. */
int main() {
    for (const Vec& rate : {Vec::Zero().eval(),Vec(1e-9,-2e-9,3e-9),Vec(0,0,4*M_PI),Vec(7,-9,4)})
        for (const double dt : {.005,.075}) {
            prediction<false>(rate,dt); prediction<true>(rate,dt);
        }
    QuaternionMEKF<double,false> f(Vec::Ones(),Vec::Zero(),Vec::Ones(),1,0,0);
    for (int i=0; i<1000; ++i) f.time_update(Vec(0,0,4*M_PI),.005);
    check((f.covariance()-Mat::Identity()).norm(),1e-10,"rotation cannot create information");
    reset(false); reset(true);
    QuaternionMEKF<double,true> bias_filter(Vec::Ones(),Vec::Zero(),Vec::Ones(),.7,.4,0);
    bias_filter.time_update(Vec(.4,.1,-.2),.1);
    const Quat q_before(bias_filter.quaternion());
    const Mat attitude_cov=bias_filter.covariance().topLeftCorner<3,3>();
    bias_filter.reacquire_gyro_bias(Vec(.01,-.02,.005));
    check(Quat(bias_filter.quaternion()).angularDistance(q_before),1e-12,"bias repair preserved attitude");
    check((bias_filter.covariance().topLeftCorner<3,3>()-attitude_cov).norm(),1e-12,"bias repair preserved attitude covariance");
    check(bias_filter.covariance().topRightCorner<3,3>().norm(),1e-12,"bias repair removed stale cross covariance");
    check((bias_filter.covariance().bottomRightCorner<3,3>()-.0009*Mat::Identity()).norm(),1e-12,"bias repair covariance floor");
    check((bias_filter.gyroscope_bias()-Vec(.01,-.02,.005)).norm(),1e-12,"bias repair mean");
    if (failures) return EXIT_FAILURE;
    std::puts("qMEKF prediction and reset agree with finite-differenced right-error geometry");
    return EXIT_SUCCESS;
}
