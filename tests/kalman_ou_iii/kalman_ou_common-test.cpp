#define EIGEN_NON_ARDUINO

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <array>
#include <cmath>
#include <iostream>

namespace detail = ocean_imu::kalman::ou_detail;

namespace {

using T = double;
using Matrix3 = Eigen::Matrix<T,3,3>;
using Matrix4 = Eigen::Matrix<T,4,4>;
using Vector3 = Eigen::Matrix<T,3,1>;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

template<int N>
Eigen::Matrix<T,N,N> van_loan_q(T tau, T dt, T sigma2) {
    Eigen::Matrix<T,N,N> F = Eigen::Matrix<T,N,N>::Zero();
    F(0,N-1) = T(1);
    F(1,0) = T(1);
    if constexpr (N == 4) F(2,1) = T(1);
    F(N-1,N-1) = -T(1) / tau;

    Eigen::Matrix<T,N,N> Qc = Eigen::Matrix<T,N,N>::Zero();
    Qc(N-1,N-1) = T(2) * sigma2 / tau;

    Eigen::Matrix<T,2*N,2*N> A = Eigen::Matrix<T,2*N,2*N>::Zero();
    A.template topLeftCorner<N,N>() = F;
    A.template topRightCorner<N,N>() = Qc;
    A.template bottomRightCorner<N,N>() = -F.transpose();
    const Eigen::Matrix<T,2*N,2*N> E = (A * dt).exp();
    const Eigen::Matrix<T,N,N> Phi = E.template topLeftCorner<N,N>();
    return E.template topRightCorner<N,N>() * Phi.transpose();
}

bool test_ou_covariance() {
    const std::array<T,4> taus = {T(0.05), T(0.5), T(2.1), T(10)};
    const std::array<T,6> ratios = {T(1e-5), T(1e-3), T(0.0099), T(0.0101), T(0.1), T(0.8)};
    const int map[3] = {0,1,3};

    for (T tau : taus) {
        for (T ratio : ratios) {
            const T dt = tau * ratio;
            const T sigma2 = T(0.83);
            Matrix3 q2;
            Matrix4 q3;
            detail::IntegratedOUChain<T,2>::process_covariance(tau, dt, sigma2, q2);
            detail::IntegratedOUChain<T,3>::process_covariance(tau, dt, sigma2, q3);
            const Matrix3 ref2 = van_loan_q<3>(tau, dt, sigma2);
            const Matrix4 ref3 = van_loan_q<4>(tau, dt, sigma2);
            const T scale2 = std::max(T(1), ref2.cwiseAbs().maxCoeff());
            const T scale3 = std::max(T(1), ref3.cwiseAbs().maxCoeff());
            if (!check((q2-ref2).cwiseAbs().maxCoeff() <= T(5e-10)*scale2, "OU-II Qd mismatch")) return false;
            if (!check((q3-ref3).cwiseAbs().maxCoeff() <= T(5e-10)*scale3, "OU-III Qd mismatch")) return false;
            for (int i=0; i<3; ++i) {
                for (int j=0; j<3; ++j) {
                    if (!check(std::abs(q2(i,j)-q3(map[i],map[j])) <= T(2e-13)*scale2,
                               "OU-II/OU-III marginal mismatch")) return false;
                }
            }
            Eigen::SelfAdjointEigenSolver<Matrix3> e2(q2);
            Eigen::SelfAdjointEigenSolver<Matrix4> e3(q3);
            if (!check(e2.eigenvalues().minCoeff() >= -T(1e-12)*scale2, "OU-II Qd is not PSD")) return false;
            if (!check(e3.eigenvalues().minCoeff() >= -T(1e-12)*scale3, "OU-III Qd is not PSD")) return false;
        }
    }
    return true;
}

template<typename Filter>
bool covariance_unchanged(Filter& filter, const Vector3& std_aw) {
    const auto before = filter.covariance_full();
    filter.set_aw_stationary_std(std_aw);
    return (filter.covariance_full() - before).cwiseAbs().maxCoeff() == T(0);
}

bool test_stationary_setters() {
    const Vector3 sigma_a = Vector3::Constant(T(0.02));
    const Vector3 gyro_density = Vector3::Constant(T(0.001));
    const Vector3 sigma_m = Vector3::Constant(T(0.5));
    const Vector3 zero = Vector3::Zero();
    const Vector3 std_aw(T(0.4), T(0.5), T(0.6));
    Matrix3 full;
    full << T(0.20), T(0.02), T(-0.01),
            T(0.02), T(0.30), T(0.015),
            T(-0.01), T(0.015), T(0.40);

    Kalman3D_Wave_OU_II<T> ou2(sigma_a, gyro_density, sigma_m);
    Kalman3D_Wave_OU_III<T> ou3(sigma_a, gyro_density, sigma_m);
    for (int i=0; i<8; ++i) {
        ou2.time_update(zero, T(0.005));
        ou3.time_update(zero, T(0.005));
    }

    if (!check(covariance_unchanged(ou2, std_aw), "OU-II std setter changed posterior covariance")) return false;
    if (!check(covariance_unchanged(ou3, std_aw), "OU-III std setter changed posterior covariance")) return false;

    auto cov_ou_ii = ou2.covariance_full();
    auto cov_ou_iii = ou3.covariance_full();
    ou2.set_aw_stationary_corr_std(std_aw, T(-0.2), T(0.15));
    ou3.set_aw_stationary_corr_std(std_aw, T(-0.2), T(0.15));
    if (!check((ou2.covariance_full()-cov_ou_ii).cwiseAbs().maxCoeff() == T(0), "OU-II corr setter changed posterior covariance")) return false;
    if (!check((ou3.covariance_full()-cov_ou_iii).cwiseAbs().maxCoeff() == T(0), "OU-III corr setter changed posterior covariance")) return false;

    cov_ou_ii = ou2.covariance_full();
    cov_ou_iii = ou3.covariance_full();
    ou2.set_aw_stationary_cov_full(full);
    ou3.set_aw_stationary_cov_full(full);
    if (!check((ou2.covariance_full()-cov_ou_ii).cwiseAbs().maxCoeff() == T(0), "OU-II full setter changed posterior covariance")) return false;
    if (!check((ou3.covariance_full()-cov_ou_iii).cwiseAbs().maxCoeff() == T(0), "OU-III full setter changed posterior covariance")) return false;

    ou2.reset_aw_covariance_to_stationary();
    ou3.reset_aw_covariance_to_stationary();
    const auto reset2 = ou2.covariance_full();
    const auto reset3 = ou3.covariance_full();
    constexpr int aw2 = 12;
    constexpr int aw3 = 15;
    if (!check((reset2.template block<3,3>(aw2,aw2)-full).cwiseAbs().maxCoeff() < T(1e-13), "OU-II reset covariance mismatch")) return false;
    if (!check((reset3.template block<3,3>(aw3,aw3)-full).cwiseAbs().maxCoeff() < T(1e-13), "OU-III reset covariance mismatch")) return false;
    for (int i=0; i<reset2.rows(); ++i) {
        if (i < aw2 || i >= aw2+3) {
            if (!check(reset2.template block<1,3>(i,aw2).cwiseAbs().maxCoeff() == T(0), "OU-II reset retained cross covariance")) return false;
        }
    }
    for (int i=0; i<reset3.rows(); ++i) {
        if (i < aw3 || i >= aw3+3) {
            if (!check(reset3.template block<1,3>(i,aw3).cwiseAbs().maxCoeff() == T(0), "OU-III reset retained cross covariance")) return false;
        }
    }
    return true;
}

bool test_noise_units_and_period() {
    const Vector3 sample_std = Vector3::Constant(T(0.02));
    const T dt = T(0.005);
    const Vector3 expected = sample_std * std::sqrt(dt);
    const Vector3 converted = Kalman3D_Wave_OU_II<T>::gyro_noise_density_from_sample_std(sample_std, dt);
    if (!check((converted-expected).cwiseAbs().maxCoeff() < T(1e-15), "gyro noise conversion mismatch")) return false;

    for (int rate : {100, 200, 400}) {
        T elapsed = T(0);
        int count = 0;
        const T step = T(1) / T(rate);
        for (int i=0; i<rate*3; ++i) {
            if (detail::periodic_update_due(step, T(0.015), elapsed)) ++count;
        }
        if (!check(count == 200, "pseudo update rate depends on sample rate")) return false;
        if (!check(elapsed >= T(0) && elapsed < T(0.015), "pseudo update remainder is out of range")) return false;
    }
    for (int rate : {100, 200, 400}) {
        float elapsed = 0.0f;
        int count = 0;
        const float step = 1.0f / static_cast<float>(rate);
        for (int sample=0; sample<rate*3; ++sample) {
            if (detail::periodic_update_due(step, 0.015f, elapsed)) ++count;
        }
        if (!check(count == 200, "float pseudo update rate depends on sample rate")) return false;
        if (!check(elapsed >= 0.0f && elapsed < 0.015f, "float pseudo update remainder is out of range")) return false;
    }
    return true;
}

bool test_period_retarget_progress() {
    // An unexpired deadline preserves accumulated service credit exactly.
    const float old_elapsed = 0.12f;
    const float longer_period = 0.14f;
    const float kept = detail::retarget_period_elapsed_progress_preserving(
        old_elapsed, longer_period);
    if (!check(kept == old_elapsed,
               "period retarget discarded unexpired elapsed credit")) return false;

    // A shortened deadline that is already overdue must not modulo the credit
    // away. Park immediately below the new deadline so the next sample fires.
    const float shorter_period = 0.13f;
    const float armed = detail::retarget_period_elapsed_progress_preserving(
        longer_period, shorter_period);
    if (!check(armed == std::nextafter(shorter_period, 0.0f),
               "overdue period retarget did not arm next-sample service")) return false;
    float due_elapsed = armed;
    if (!check(detail::periodic_update_due(0.005f, shorter_period, due_elapsed),
               "overdue retarget was not serviced on the next valid sample")) return false;

    // The largest deployed tau-coupled pseudo period fires on sample 33 from
    // zero elapsed under the exact float scheduler. This is the source-
    // independent worst case once retargets preserve progress.
    const float max_deployed_period = 0.16363636f;
    float elapsed = 0.0f;
    int first_fire = 0;
    for (int sample = 1; sample <= 64; ++sample) {
        elapsed = detail::retarget_period_elapsed_progress_preserving(
            elapsed, max_deployed_period);
        if (detail::periodic_update_due(0.005f, max_deployed_period, elapsed)) {
            first_fire = sample;
            break;
        }
    }
    if (!check(first_fire == 33,
               "largest deployed pseudo period no longer fires on sample 33")) return false;

    // Replay the former H,H,L,H,L,... 13-sample starvation timing. The new
    // retarget rule must produce updates rather than allowing a 635-sample
    // no-fire execution.
    constexpr float low = 0.1300000101327896f;
    constexpr float high = 0.13000193238258362f;
    elapsed = 0.0f;
    int fires = 0;
    int since_fire = 0;
    int worst_gap = 0;
    int samples = 0;
    int segment = 0;
    while (samples < 635) {
        const float period = (segment < 2 || (segment % 2) == 1) ? high : low;
        elapsed = detail::retarget_period_elapsed_progress_preserving(elapsed, period);
        const int used = std::min(13, 635 - samples);
        for (int k = 0; k < used; ++k) {
            ++samples;
            ++since_fire;
            if (detail::periodic_update_due(0.005f, period, elapsed)) {
                ++fires;
                worst_gap = std::max(worst_gap, since_fire);
                since_fire = 0;
            }
        }
        ++segment;
    }
    worst_gap = std::max(worst_gap, since_fire);
    if (!check(fires > 0, "former starvation cycle still has zero pseudo updates")) return false;
    if (!check(worst_gap <= 33, "former starvation cycle exceeds 33-sample gap")) return false;
    return true;
}

// The congruent a_w re-alignment must reach the stationary marginal without
// changing how a_w correlates with the rest of the state, and without pushing
// the joint covariance out of the positive semi-definite cone. Overwriting the
// marginal alone intentionally does not preserve this invariant.
bool test_aw_covariance_sync() {
    const Vector3 sigma_a = Vector3::Constant(T(0.35));
    const Vector3 gyro_density = Vector3::Constant(T(0.0015));
    const Vector3 sigma_m = Vector3::Constant(T(0.8));
    constexpr int aw = 15;

    Kalman3D_Wave_OU_III<T> ou3(sigma_a, gyro_density, sigma_m);
    ou3.set_aw_stationary_std(Vector3(T(2.36), T(2.36), T(1.26)));
    ou3.reset_aw_covariance_to_stationary();

    // Let the filter learn cross-covariances between a_w and [v, p, S].
    for (int i = 0; i < 4000; ++i) {
        const T t = T(i) * T(0.005);
        ou3.time_update(Vector3(T(0.02)*std::sin(T(0.7)*t), T(0.03)*std::cos(T(0.5)*t), T(0.01)), T(0.005));
        ou3.measurement_update_acc_only(Vector3(T(0.6)*std::sin(T(0.55)*t),
                                                T(0.4)*std::cos(T(0.61)*t),
                                                T(-9.80665) + T(1.4)*std::sin(T(0.57)*t)));
    }

    const auto before = ou3.covariance_full();
    const Matrix3 marginal_before = before.template block<3,3>(aw, aw);

    // Drive the stationary scale well away from the posterior in both
    // directions, which is what the adaptation cadence does in practice.
    for (T scale : {T(4.0), T(0.25)}) {
        const Vector3 target_std(T(2.36)*scale, T(2.36)*scale, T(1.26)*scale);
        Kalman3D_Wave_OU_III<T> probe = ou3;
        probe.set_aw_stationary_std(target_std);
        probe.synchronize_aw_covariance_to_stationary_congruent();
        const auto after = probe.covariance_full();

        const Matrix3 target = target_std.array().square().matrix().asDiagonal();
        const Matrix3 marginal_after = after.template block<3,3>(aw, aw);
        if (!check((marginal_after - target).cwiseAbs().maxCoeff() < T(1e-10),
                   "a_w sync did not reach the stationary marginal")) return false;

        // The quantity a congruence preserves is the cross-covariance measured
        // against the whitened a_w block, P_x,aw L^-T.
        const Eigen::LLT<Matrix3> chol_before(marginal_before);
        const Eigen::LLT<Matrix3> chol_after(marginal_after);
        if (!check(chol_before.info() == Eigen::Success &&
                   chol_after.info() == Eigen::Success,
                   "a_w marginal is not positive definite")) return false;
        const Matrix3 factor_before = chol_before.matrixL();
        const Matrix3 factor_after = chol_after.matrixL();

        for (int i = 0; i < after.rows(); ++i) {
            if (i >= aw && i < aw + 3) continue;
            if (!check(std::abs(after(i,i) - before(i,i)) < T(1e-10),
                       "a_w sync changed a non-a_w marginal")) return false;

            const Eigen::Matrix<T,1,3> whitened_before =
                factor_before.template triangularView<Eigen::Lower>()
                    .solve(before.template block<1,3>(i, aw).transpose()).transpose();
            const Eigen::Matrix<T,1,3> whitened_after =
                factor_after.template triangularView<Eigen::Lower>()
                    .solve(after.template block<1,3>(i, aw).transpose()).transpose();
            if (!check((whitened_after - whitened_before).cwiseAbs().maxCoeff() < T(1e-8),
                       "a_w sync changed the whitened cross-covariance")) return false;
        }

        const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> joint = after;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>
            solver(joint);
        const T smallest = solver.eigenvalues().minCoeff();
        if (!check(smallest > -T(1e-9),
                   "a_w sync produced an indefinite covariance")) return false;
    }
    return true;
}

bool test_attitude_helpers() {
    Vector3 rotation_vector;
    rotation_vector << T(0.01), T(-0.02), T(0.03);
    const auto q = detail::quat_from_delta_theta(rotation_vector);
    if (!check(std::abs(q.norm()-T(1)) <= T(1e-14), "quaternion increment is not unit length")) return false;
    Matrix3 R, B;
    detail::rot_and_B_from_wt(rotation_vector, T(0.005), R, B);
    return check((R*R.transpose()-Matrix3::Identity()).cwiseAbs().maxCoeff() <= T(1e-13),
                 "rotation transition is not orthogonal");
}

// Integrated-OU transition and process noise, entry by entry, against
// 60-digit references at the exact float inputs, in the precision the device
// runs (float) and in double.  Relative error per entry: the S rows are many
// orders of magnitude below the a-a entry, so a tolerance scaled by the matrix
// norm cannot see them.  tau from 0.0055 to 12 s spans x = h/tau on both sides
// of the series/closed-form switch at x = 1.
struct OUReference { float h, tau; double q[12]; };
// q: Qd (v,p,S,a) entries vv vp vS va pp pS pa SS Sa aa, then phi_pa phi_Sa.
const OUReference ou_reference[] = {
    {0.005f, 0.0055f, {6.7143097625450152268e-6, 1.3440736070324864007e-8, 1.851793478002004698e-11, 1.6276048791267115021e-3, 2.9460667358844910092e-11, 4.2816595965340214248e-14, 2.6400918980883641889e-6, 6.4394573140741892801e-17, 3.2551916368479006764e-9, 6.9527387459089979971e-1, 9.4374318459299975368e-6, 1.6844121682748924855e-8}},
    {0.005f, 0.05f, {1.2842569340985659308e-6, 2.4278133290391086705e-9, 3.2503809427211441229e-12, 3.7582052717081720153e-4, 4.9093159757997026055e-12, 6.8560310018223481572e-15, 6.2615881666545507109e-7, 9.855067438105146993e-18, 7.8256807223042950575e-10, 1.504534668570021214e-1, 1.2093544564044059247e-5, 2.0322744160952593561e-8}},
    {0.005f, 0.45f, {1.5242944550179065456e-7, 2.860694951971582438e-10, 3.816022911054479295e-13, 4.560206868874563154e-5, 5.7284494746930469338e-13, 7.9610859119663673363e-16, 7.6003133345759886981e-8, 1.1380870695886721409e-18, 9.5003719077734313185e-11, 1.8241015545666249321e-2, 1.2453831462355286766e-5, 2.0775589932663423404e-8}},
    {0.005f, 0.6f, {1.1456007074541199939e-7, 2.1494914952657447304e-10, 2.8669825471914423572e-13, 3.4296531277159524557e-5, 4.3029618517030234038e-13, 5.9791004984742291512e-16, 5.716075186777024745e-8, 8.5460178010069073563e-19, 7.1450855539906979503e-11, 1.371869220751742797e-2, 1.2465349440247035224e-5, 2.0790001401776150683e-8}},
    {0.005f, 2.0f, {3.4518562159589984702e-8, 6.4735782482186711225e-11, 8.6323362583617187913e-14, 1.0349099614870126822e-5, 1.2950752210884511922e-13, 1.7989653094851041489e-16, 1.7248495379146550897e-8, 2.5703518613337733061e-19, 2.1560616496118310178e-11, 4.1396420945383607236e-3, 1.2489589282400519197e-5, 2.0820317611884556814e-8}},
    {0.005f, 10.0f, {6.914072919174789676e-9, 1.2964426564099938167e-11, 1.7286261796318441418e-14, 2.073962668146755498e-6, 2.5930293006120054851e-14, 3.6015295388458892211e-17, 3.4566043408450942242e-9, 5.1452028596594156942e-20, 4.3207553114768118189e-12, 8.2958510308439560674e-4, 1.2497916368403430415e-5, 2.0830728030310604331e-8}},
    {0.0005f, 12.0f, {5.763709476245020937e-12, 1.0806993305233637068e-15, 1.4409350107426044636e-19, 1.7290947492905320269e-8, 2.1614087701640339703e-19, 3.0019637167611317971e-23, 2.8818247188634371932e-12, 4.2885309670514439114e-27, 3.6022810695742253804e-16, 6.9163786696524026667e-5, 1.2499827578108895814e-7, 2.0833119289802498547e-11}},
    {0.005f, 0.005f, {6.9757860379617905585e-6, 1.404103441194239584e-8, 1.9398786901652638151e-11, 1.6582419933011205025e-3, 3.1028311325545653129e-11, 4.5276084386145463825e-14, 2.6747958908690268058e-6, 6.840621779522758588e-17, 3.2885592845470394893e-9, 7.1767170048296214383e-1, 9.1969856181487486748e-6, 1.6515068746148003483e-8}},
    {0.02f, 0.01f, {1.2641110781395038499e-4, 1.0698584571985450978e-6, 6.0675756729800097043e-9, 6.2054538375692849525e-3, 1.0228382936822178467e-8, 6.2054532827584960664e-11, 3.6548485568932922554e-5, 3.9207071362796336123e-13, 1.7123225483423176402e-7, 8.1479800333872529519e-1, 1.1353352324831710654e-4, 8.6466465878310115001e-7}},
    {0.05f, 0.02f, {7.7077188191853620225e-4, 1.661987348812691323e-5, 2.3810597638859340512e-7, 1.3986627604980880495e-2, 4.1022967278704976823e-7, 6.3228386529579949936e-9, 1.9350189937042067505e-4, 1.0174461437950436689e-10, 2.1881237093741984924e-6, 8.2440748845559360937e-1, 6.3283400535467340246e-4, 1.2343320362069881454e-5}},
    {0.1f, 0.02f, {2.3329328916567723107e-3, 1.0659822175298168197e-4, 3.1688765560137249103e-6, 1.637705314449223485e-2, 5.8695648017422780009e-6, 1.9159189720561429766e-7, 3.096149264588804118e-4, 6.5884989993024642671e-9, 5.4323219855977297299e-6, 8.2996230138379190814e-1, 1.602695181157365566e-3, 6.7946097838371281134e-5}},
};

template<typename U>
bool test_ou_entries_precision(U tol_series, U tol_closed, const char* label) {
    const int r[10] = {0,0,0,0,1,1,1,2,2,3};
    const int c[10] = {0,1,2,3,1,2,3,2,3,3};
    for (const OUReference& ref : ou_reference) {
        const U h = U(ref.h), tau = U(ref.tau), sigma2 = U(0.83f);
        Eigen::Matrix<U,4,4> q;
        Eigen::Matrix<U,4,4> phi;
        detail::IntegratedOUChain<U,3>::process_covariance(tau, h, sigma2, q);
        detail::IntegratedOUChain<U,3>::transition(tau, h, phi);
        U got[12];
        for (int k = 0; k < 10; ++k) got[k] = q(r[k], c[k]);
        got[10] = phi(1,3);
        got[11] = phi(2,3);
        const U tol = (h / tau < U(1)) ? tol_series : tol_closed;
        for (int k = 0; k < 12; ++k) {
            const double rel = std::abs(double(got[k]) - ref.q[k]) / std::abs(ref.q[k]);
            if (!(rel <= double(tol))) {
                std::cerr << label << " OU entry " << k << " at h=" << ref.h << " tau=" << ref.tau
                          << " relative error " << rel << '\n';
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main() {
    if (!test_ou_covariance()) return 1;
    if (!test_ou_entries_precision<float>(2e-6f, 3e-5f, "float")) return 1;
    if (!test_ou_entries_precision<double>(1e-14, 1e-13, "double")) return 1;
    if (!test_stationary_setters()) return 1;
    if (!test_noise_units_and_period()) return 1;
    if (!test_period_retarget_progress()) return 1;
    if (!test_aw_covariance_sync()) return 1;
    if (!test_attitude_helpers()) return 1;
    return 0;
}
