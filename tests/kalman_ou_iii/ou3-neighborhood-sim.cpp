#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#include "util/W3dSimCommon.h"

// Host-only nonlinear certificate observer.  No production visibility or
// behavior is changed: this translation unit exposes the MEKF internals only so
// a controlled error can be injected into a second estimator that receives the
// exact same noisy measurements as the nominal estimator.
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

using Eigen::Matrix3f;
using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

constexpr float kSigmaARescale = 0.71f;
constexpr float kSigmaGRescale = 0.05f;
constexpr float kSigmaMRescale = 2.0f;
constexpr float kMagOdrHz = 25.0f;
constexpr float kDt = 1.0f / 200.0f;
constexpr int kNX = 21;
constexpr int kOffBg = 3;
constexpr int kOffV = 6;
constexpr int kOffP = 9;
constexpr int kOffS = 12;
constexpr int kOffAw = 15;
constexpr int kOffBa = 18;
using Vector21f = Eigen::Matrix<float, kNX, 1>;
using Matrix21f = Eigen::Matrix<float, kNX, kNX>;

float env_float(const char* name, float fallback)
{
    if (const char* raw = std::getenv(name)) {
        char* end = nullptr;
        const float x = std::strtof(raw, &end);
        if (end != raw && std::isfinite(x)) return x;
    }
    return fallback;
}

int env_positive_int(const char* name, int fallback)
{
    if (const char* raw = std::getenv(name)) {
        const int x = std::atoi(raw);
        if (x > 0) return x;
    }
    return fallback;
}

std::string env_string(const char* name, const std::string& fallback)
{
    if (const char* raw = std::getenv(name)) {
        if (*raw) return raw;
    }
    return fallback;
}

Vector21f parse_delta()
{
    Vector21f d = Vector21f::Zero();
    const char* raw = std::getenv("OU3_NEIGHBOR_DELTA");
    if (!raw || !*raw) return d;

    std::stringstream ss(raw);
    std::string token;
    int i = 0;
    while (std::getline(ss, token, ',')) {
        if (i >= kNX)
            throw std::runtime_error("OU3_NEIGHBOR_DELTA has more than 21 entries");
        char* end = nullptr;
        const float x = std::strtof(token.c_str(), &end);
        if (end == token.c_str() || !std::isfinite(x))
            throw std::runtime_error("OU3_NEIGHBOR_DELTA contains a non-finite/non-numeric entry");
        d(i++) = x;
    }
    if (i != kNX)
        throw std::runtime_error("OU3_NEIGHBOR_DELTA must contain exactly 21 comma-separated entries");
    return d;
}

Vector3f quaternion_log(const Quaternionf& q_in)
{
    Quaternionf q = q_in;
    if (!(q.norm() > 1.0e-12f) || !q.coeffs().allFinite())
        return Vector3f::Constant(std::numeric_limits<float>::quiet_NaN());
    q.normalize();
    if (q.w() < 0.0f) q.coeffs() *= -1.0f;
    const Vector3f v(q.x(), q.y(), q.z());
    const float nv = v.norm();
    if (nv < 1.0e-9f) return 2.0f * v;
    const float angle = 2.0f * std::atan2(nv, std::max(0.0f, q.w()));
    return (angle / nv) * v;
}

bool close_source_scalar(float a, float b)
{
    const float scale = std::max({1.0f, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 2.0e-6f * scale;
}

class NeighborhoodAdapter final : public IW3dFusionAdapter {
public:
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    using Mekf = std::remove_reference_t<decltype(std::declval<Fusion&>().raw().mekf())>;

    NeighborhoodAdapter(bool with_mag,
                        const Vector3f& sigma_a_init,
                        const Vector3f& sigma_g,
                        const Vector3f& sigma_m)
        : inject_requested_s_(env_float("OU3_NEIGHBOR_INJECT_TIME_S", 300.0f)),
          horizon_s_(env_float("OU3_NEIGHBOR_HORIZON_S", 4.0f)),
          requested_mode_(env_string("OU3_NEIGHBOR_MODE", "A")),
          trace_stride_(env_positive_int("OU3_NEIGHBOR_TRACE_STRIDE", 50)),
          delta_(parse_delta())
    {
        if (!(horizon_s_ > 0.0f))
            throw std::runtime_error("OU3_NEIGHBOR_HORIZON_S must be positive");
        if (requested_mode_ != "H" && requested_mode_ != "A")
            throw std::runtime_error("OU3_NEIGHBOR_MODE must be H or A");
        if (requested_mode_ == "H" && delta_.tail<3>().norm() > 0.0f)
            throw std::runtime_error("held-mode perturbation must have zero accelerometer-bias entries");

        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init * kSigmaARescale;
        cfg.sigma_g = sigma_g * kSigmaGRescale;
        cfg.sigma_m = sigma_m * kSigmaMRescale;
        cfg.mag_delay_sec = MAG_DELAY_SEC;

        nominal_.begin(cfg);
        perturbed_.begin(cfg);
        configure_filter(nominal_);
        configure_filter(perturbed_);

        const char* path = std::getenv("OU3_NEIGHBOR_TRACE");
        if (!path || !*path)
            throw std::runtime_error("OU3_NEIGHBOR_TRACE must name the pairwise trace CSV");
        trace_.open(path);
        if (!trace_) throw std::runtime_error("cannot open OU3_NEIGHBOR_TRACE");
        trace_ << std::setprecision(10);
        write_header();
    }

    void updateMag(const Vector3f& mag_body_ned) override
    {
        if (finished_) return;
        nominal_.updateMag(mag_body_ned);
        perturbed_.updateMag(mag_body_ned);
        const bool a = nominal_.raw().mekf().lastMagDiag().accepted;
        const bool b = perturbed_.raw().mekf().lastMagDiag().accepted;
        last_mag_accept_match_ = (a == b);
        if (injected_ && !last_mag_accept_match_) source_match_ = false;
        if (injected_) check_source_match();
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        if (finished_) {
            time_s_ += dt;
            return;
        }

        maybe_inject();

        nominal_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
        perturbed_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
        time_s_ += dt;
        ++sample_index_;

        const bool aa = nominal_.raw().mekf().lastAccDiag().accepted;
        const bool ab = perturbed_.raw().mekf().lastAccDiag().accepted;
        last_acc_accept_match_ = (aa == ab);
        if (injected_ && !last_acc_accept_match_) source_match_ = false;

        if (!injected_) return;
        check_source_match();
        update_prefix_stats();

        const bool stride_due =
            samples_since_injection_ % static_cast<unsigned>(trace_stride_) == 0u;
        const bool endpoint_due = time_s_ + 0.5f * dt >= inject_actual_s_ + horizon_s_;
        if (stride_due || endpoint_due) write_row(endpoint_due ? 1 : 0);
        ++samples_since_injection_;

        if (endpoint_due) {
            finished_ = true;
            std::cout << "OU3_NEIGHBOR_DONE"
                      << " inject_s=" << inject_actual_s_
                      << " horizon_s=" << horizon_s_
                      << " mode=" << injected_mode_
                      << " source_match=" << (source_match_ ? 1 : 0)
                      << " max_theta_deg=" << rad_to_deg(max_theta_rad_)
                      << " max_cov_rel=" << max_cov_rel_
                      << "\n";
        }
    }

    FilterSnapshot snapshot() const override
    {
        const auto& filter = nominal_.raw();
        const auto& mekf = filter.mekf();
        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(mekf.get_position());
        s.vel_est_zu = ned_to_zu(mekf.get_velocity());
        s.acc_est_zu = ned_to_zu(mekf.get_world_accel());
        const Quaternionf q_bw_ned = mekf.quaternion_boat().normalized();
        float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll, pitch, yaw);
        s.euler_nautical_deg = Vector3f(roll, pitch, wrapDeg(yaw));
        s.acc_bias_est_ned = mekf.get_acc_bias();
        s.gyro_bias_est_ned = mekf.gyroscope_bias();
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(mekf) + nominal_.magHardIronBodyUT();
        s.tau_target = filter.getTauTarget();
        s.sigma_target = filter.getSigmaTarget();
        s.tuning_target = filter.getRSTarget();
        s.tau_applied = filter.getTauApplied();
        s.sigma_applied = filter.getSigmaApplied();
        s.tuning_applied = filter.getRSApplied();
        s.freq_hz = filter.getFreqHz();
        s.wave_period_sec = filter.getWavePeriodSec();
        s.period_sec = filter.getPeriodSec();
        s.accel_variance = filter.getAccelVariance();
        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps = filter.getVerticalSpeedEnvelopeMps(true);
        return s;
    }

private:
    static void configure_filter(Fusion& f)
    {
        auto& filter = f.raw();
        filter.setPeriodicAwCovarianceSync(true);
        filter.setAwCovarianceSyncCongruent(false);
        filter.enableTuner(true);
        filter.enableClamp(true);
    }

    bool source_ready(const Fusion& f) const
    {
        const auto& m = f.raw().mekf();
        if (!f.isLive()) return false;
        const bool active = m.acc_bias_updates_enabled();
        return requested_mode_ == (active ? "A" : "H");
    }

    void maybe_inject()
    {
        if (injected_ || time_s_ + 1.0e-7f < inject_requested_s_) return;
        if (!source_ready(nominal_) || !source_ready(perturbed_)) return;
        check_source_match();
        if (!source_match_) return;

        inject_delta(perturbed_.raw().mekf(), delta_);
        injected_ = true;
        inject_actual_s_ = time_s_;
        injected_mode_ = nominal_.raw().mekf().acc_bias_updates_enabled() ? "A" : "H";
        samples_since_injection_ = 0;
        update_prefix_stats();
        write_row(0);
        std::cout << "OU3_NEIGHBOR_INJECT"
                  << " requested_s=" << inject_requested_s_
                  << " actual_s=" << inject_actual_s_
                  << " mode=" << injected_mode_
                  << " W0=" << information_energy(nominal_.raw().mekf())
                  << " theta_deg=" << rad_to_deg(pair_error().head<3>().norm())
                  << "\n";
    }

    static void inject_delta(Mekf& m, const Vector21f& d)
    {
        // Attitude uses the filter's exact left-error retraction and covariance
        // reset.  Remaining coordinates are ordinary additive states.
        m.xext.template head<3>() = d.template head<3>();
        m.applyQuaternionCorrectionFromErrorState();
        m.xext.template segment<3>(kOffBg) += d.template segment<3>(kOffBg);
        m.xext.template segment<3>(kOffV) += d.template segment<3>(kOffV);
        m.xext.template segment<3>(kOffP) += d.template segment<3>(kOffP);
        m.xext.template segment<3>(kOffS) += d.template segment<3>(kOffS);
        m.xext.template segment<3>(kOffAw) += d.template segment<3>(kOffAw);
        if (m.acc_bias_updates_enabled())
            m.xext.template segment<3>(kOffBa) += d.template segment<3>(kOffBa);
    }

    Vector21f pair_error() const
    {
        const auto& n = nominal_.raw().mekf();
        const auto& p = perturbed_.raw().mekf();
        Vector21f e = Vector21f::Zero();

        // qref is WORLD->BODY'.  The left-error quaternion that maps the
        // nominal local frame to the perturbed local frame is q_p q_n^{-1}; at
        // injection this is exactly Exp(delta_theta).
        const Quaternionf qrel = p.qref * n.qref.conjugate();
        e.template head<3>() = quaternion_log(qrel);
        e.template segment<3>(kOffBg) =
            p.xext.template segment<3>(kOffBg) - n.xext.template segment<3>(kOffBg);
        e.template segment<3>(kOffV) =
            p.xext.template segment<3>(kOffV) - n.xext.template segment<3>(kOffV);
        e.template segment<3>(kOffP) =
            p.xext.template segment<3>(kOffP) - n.xext.template segment<3>(kOffP);
        e.template segment<3>(kOffS) =
            p.xext.template segment<3>(kOffS) - n.xext.template segment<3>(kOffS);
        e.template segment<3>(kOffAw) =
            p.xext.template segment<3>(kOffAw) - n.xext.template segment<3>(kOffAw);
        e.template segment<3>(kOffBa) =
            p.xext.template segment<3>(kOffBa) - n.xext.template segment<3>(kOffBa);
        return e;
    }

    float information_energy(const Mekf& reference) const
    {
        const Vector21f e = pair_error();
        const int dim = injected_mode_ == "H" ? 18 : 21;
        const Matrix21f Pfull = reference.covariance_full();
        if (dim == 18) {
            const Eigen::Matrix<float,18,18> P =
                0.5f * (Pfull.topLeftCorner<18,18>() + Pfull.topLeftCorner<18,18>().transpose());
            Eigen::LDLT<Eigen::Matrix<float,18,18>> ldlt(P);
            if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
            const Eigen::Matrix<float,18,1> x = e.head<18>();
            return x.dot(ldlt.solve(x));
        }
        const Matrix21f P = 0.5f * (Pfull + Pfull.transpose());
        Eigen::LDLT<Matrix21f> ldlt(P);
        if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
        return e.dot(ldlt.solve(e));
    }

    float covariance_relative_difference() const
    {
        const Matrix21f Pn = nominal_.raw().mekf().covariance_full();
        const Matrix21f Pp = perturbed_.raw().mekf().covariance_full();
        return (Pp - Pn).norm() / std::max(1.0e-12f, Pn.norm());
    }

    void check_source_match()
    {
        if (!injected_) return;
        const auto& nf = nominal_.raw();
        const auto& pf = perturbed_.raw();
        const auto& n = nf.mekf();
        const auto& p = pf.mekf();
        const bool flags =
            nominal_.isLive() == perturbed_.isLive()
            && n.acc_bias_updates_enabled() == p.acc_bias_updates_enabled()
            && nominal_.hasMagNorthLock() == perturbed_.hasMagNorthLock()
            && nominal_.hasRefinedMagReference() == perturbed_.hasRefinedMagReference();
        const bool params =
            close_source_scalar(nf.getTauApplied(), pf.getTauApplied())
            && close_source_scalar(nf.getSigmaApplied(), pf.getSigmaApplied())
            && close_source_scalar(nf.getRSApplied(), pf.getRSApplied())
            && close_source_scalar(n.get_pseudo_update_period_s(), p.get_pseudo_update_period_s());
        if (!(flags && params && last_acc_accept_match_ && last_mag_accept_match_))
            source_match_ = false;
    }

    void update_prefix_stats()
    {
        const Vector21f e = pair_error();
        const float theta = e.head<3>().norm();
        max_theta_rad_ = std::max(max_theta_rad_, theta);
        const float c = covariance_relative_difference();
        if (std::isfinite(c)) max_cov_rel_ = std::max(max_cov_rel_, c);
    }

    void write_header()
    {
        trace_ << "time_s,time_from_injection_s,endpoint,mode,source_match,"
                  "acc_accept_match,mag_accept_match,W_nominal,W_perturbed,"
                  "theta_rad,error_norm,covariance_rel_fro,"
                  "nom_live,nom_active,nom_mag_lock,nom_mag_refined,"
                  "nom_tau,nom_sigma,nom_rs,pert_tau,pert_sigma,pert_rs";
        for (int i = 0; i < kNX; ++i) trace_ << ",e" << i;
        trace_ << '\n';
    }

    void write_row(int endpoint)
    {
        const auto& nf = nominal_.raw();
        const auto& pf = perturbed_.raw();
        const auto& n = nf.mekf();
        const auto& p = pf.mekf();
        const Vector21f e = pair_error();
        trace_ << time_s_ << ',' << (time_s_ - inject_actual_s_) << ',' << endpoint
               << ',' << injected_mode_ << ',' << (source_match_ ? 1 : 0)
               << ',' << (last_acc_accept_match_ ? 1 : 0)
               << ',' << (last_mag_accept_match_ ? 1 : 0)
               << ',' << information_energy(n)
               << ',' << information_energy(p)
               << ',' << e.head<3>().norm()
               << ',' << e.norm()
               << ',' << covariance_relative_difference()
               << ',' << (nominal_.isLive() ? 1 : 0)
               << ',' << (n.acc_bias_updates_enabled() ? 1 : 0)
               << ',' << (nominal_.hasMagNorthLock() ? 1 : 0)
               << ',' << (nominal_.hasRefinedMagReference() ? 1 : 0)
               << ',' << nf.getTauApplied() << ',' << nf.getSigmaApplied() << ',' << nf.getRSApplied()
               << ',' << pf.getTauApplied() << ',' << pf.getSigmaApplied() << ',' << pf.getRSApplied();
        for (int i = 0; i < kNX; ++i) trace_ << ',' << e(i);
        trace_ << '\n';
    }

    mutable Fusion nominal_;
    mutable Fusion perturbed_;
    std::ofstream trace_;
    float inject_requested_s_ = 300.0f;
    float inject_actual_s_ = std::numeric_limits<float>::quiet_NaN();
    float horizon_s_ = 4.0f;
    std::string requested_mode_ = "A";
    std::string injected_mode_ = "";
    int trace_stride_ = 50;
    Vector21f delta_ = Vector21f::Zero();
    float time_s_ = 0.0f;
    unsigned sample_index_ = 0;
    unsigned samples_since_injection_ = 0;
    bool injected_ = false;
    bool finished_ = false;
    bool source_match_ = true;
    bool last_acc_accept_match_ = true;
    bool last_mag_accept_match_ = true;
    float max_theta_rad_ = 0.0f;
    float max_cov_rel_ = 0.0f;
};

void process_one(const std::string& filename,
                 bool with_mag,
                 const W3dRandomSeeds& seeds)
{
    auto result = process_wave_file_for_tracker<NeighborhoodAdapter>(
        filename, kDt, with_mag, true, kMagOdrHz,
        "_fusion_ou3_neighbor", "_fusion_ou3_neighbor_nomag", seeds, false);
    if (!result) throw std::runtime_error("neighborhood simulation did not produce a result");
}

} // namespace

int main(int argc, char** argv)
{
    bool with_mag = true;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") with_mag = false;
        else if (arg == "--input" && i + 1 < argc) files.emplace_back(argv[++i]);
        else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [--nomag] --input PATH\n"
                      << "Environment: OU3_NEIGHBOR_TRACE, OU3_NEIGHBOR_DELTA (21 CSV values), "
                      << "OU3_NEIGHBOR_INJECT_TIME_S, OU3_NEIGHBOR_HORIZON_S, OU3_NEIGHBOR_MODE=H|A\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown or incomplete argument: " << arg << "\n";
            return 2;
        }
    }
    if (files.size() != 1u) {
        std::cerr << "ERROR: ou3-neighborhood-sim requires exactly one --input\n";
        return 2;
    }

    W3dRandomSeeds seeds;
    try { seeds = w3d_random_seeds_from_env(); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    try { process_one(files.front(), with_mag, seeds); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
