#define main live_basin_euclidean_main
#include "live_basin_diagnostic.cpp"
#undef main

namespace {

double covariance_metric_norm(const Mat21& Psi,
                              const Mat21& P0,
                              const Mat21& P1)
{
    Eigen::LLT<Mat21> llt0(P0);
    Eigen::LLT<Mat21> llt1(P1);
    if (llt0.info() != Eigen::Success || llt1.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }

    const Mat21 L0 = llt0.matrixL();
    const Mat21 B = llt1.matrixL().solve(Psi * L0);
    return norm2(B);
}

std::pair<double,double> scaled_cov_eigen_bounds(const Mat21& P)
{
    const Mat21 S = state_scale();
    Mat21 D = Mat21::Zero();
    for (int i = 0; i < 21; ++i) D(i,i) = 1.0 / S(i,i);
    const Mat21 Pz = D * P * D;
    Eigen::SelfAdjointEigenSolver<Mat21> es(Pz, Eigen::EigenvaluesOnly);
    if (es.info() != Eigen::Success) {
        return {0.0, std::numeric_limits<double>::infinity()};
    }
    return {es.eigenvalues()(0), es.eigenvalues()(20)};
}

} // namespace

int main()
{
    std::cout << std::setprecision(12);
    std::cout << "name,chiV30,pz_min,pz_max,kappa_euclid\n";

    double worst_chi = 0.0;
    double global_pmin = std::numeric_limits<double>::infinity();
    double global_pmax = 0.0;
    bool all_contract = true;

    for (const auto& op : kReferencePoints) {
        const Vec3 sigma_a = Vec3::Constant(0.0294);
        const Vec3 sigma_g = Vec3::Constant(0.000157);
        const Vec3 sigma_m = Vec3::Constant(0.36);
        Core f(sigma_a, sigma_g, sigma_m);

        const Vec3 mag_world(20.0, 5.0, 44.0);
        f.set_mag_world_ref(mag_world);
        f.initialize_from_attitude(Eigen::Quaterniond::Identity(), 0.035, 0.087);
        f.set_linear_block_enabled(true);
        f.set_acc_bias_updates_enabled(true);
        f.set_aw_time_constant(op.tau);
        f.set_aw_stationary_std(Vec3::Constant(op.sigma_aw));
        f.set_RS_noise(Vec3::Constant(op.rS));
        const double period = std::clamp((0.015 / 1.1) * op.tau, 0.005, 0.25);
        f.set_pseudo_update_period_s(period);
        f.reset_aw_covariance_to_stationary();

        constexpr int warm_steps = 120 * 200;
        for (int k = 0; k < warm_steps; ++k) {
            (void)step(f, k, mag_world);
        }

        const Mat21 P0 = f.Pext;
        Mat21 Psi = Mat21::Identity();
        auto [pmin, pmax] = scaled_cov_eigen_bounds(P0);

        constexpr int H = 30 * 200;
        for (int k = 0; k < H; ++k) {
            const StepResult sr = step(f, warm_steps + k, mag_world);
            Psi = sr.A * Psi;
            const auto [lo, hi] = scaled_cov_eigen_bounds(f.Pext);
            pmin = std::min(pmin, lo);
            pmax = std::max(pmax, hi);
        }

        const Mat21 P1 = f.Pext;
        const double chi = covariance_metric_norm(Psi, P0, P1);
        const double kappa = (pmin > 0.0 && std::isfinite(pmax))
                           ? std::sqrt(pmax / pmin)
                           : std::numeric_limits<double>::infinity();

        std::cout << op.name << ',' << chi << ',' << pmin << ',' << pmax
                  << ',' << kappa << '\n';

        worst_chi = std::max(worst_chi, chi);
        global_pmin = std::min(global_pmin, pmin);
        global_pmax = std::max(global_pmax, pmax);
        all_contract = all_contract && std::isfinite(chi) && chi < 1.0;
    }

    const double global_kappa = (global_pmin > 0.0 && std::isfinite(global_pmax))
                              ? std::sqrt(global_pmax / global_pmin)
                              : std::numeric_limits<double>::infinity();
    std::cout << "PHASE2_METRIC_SUMMARY worst_chiV30=" << worst_chi
              << " global_pz_min=" << global_pmin
              << " global_pz_max=" << global_pmax
              << " global_kappa_euclid=" << global_kappa
              << " all_reference_contract=" << (all_contract ? 1 : 0)
              << '\n';

    return all_contract ? 0 : 1;
}
