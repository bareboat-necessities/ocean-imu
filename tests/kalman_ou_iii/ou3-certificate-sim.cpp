#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

constexpr float kSigmaARescale = 0.71f;
constexpr float kSigmaGRescale = 0.05f;
constexpr float kSigmaMRescale = 2.0f;
constexpr float kMagOdrHz = 25.0f;
constexpr float kDt = 1.0f / 200.0f;

struct Limits {
    static constexpr W3dFailureLimits value{
        .err_limit_percent_z_jonswap = 4.489f,
        .err_limit_percent_z_pmstokes = 4.462f,
        .err_limit_yaw_deg = 0.9004f,
        .err_limit_roll_deg = 0.3513f,
        .err_limit_pitch_deg = 0.1975f,
        .err_limit_percent_3d_jonswap = 13.98f,
        .err_limit_percent_3d_pmstokes = 14.51f,
        .acc_z_bias_percent = 4.301f,
        .bias_3d_percent = 78.92f,
        .gyro_bias_3d_percent = 15.76f,
    };
};

int env_positive_int(const char* name, int fallback)
{
    if (const char* raw = std::getenv(name)) {
        const int value = std::atoi(raw);
        if (value > 0) return value;
    }
    return fallback;
}

class CertificateAdapter final : public IW3dFusionAdapter {
public:
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

    CertificateAdapter(bool with_mag,
                       const Vector3f& sigma_a_init,
                       const Vector3f& sigma_g,
                       const Vector3f& sigma_m)
        : trace_stride_(env_positive_int("OU3_CERT_TRACE_STRIDE", 10))
    {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init * kSigmaARescale;
        cfg.sigma_g = sigma_g * kSigmaGRescale;
        cfg.sigma_m = sigma_m * kSigmaMRescale;
        cfg.mag_delay_sec = MAG_DELAY_SEC;
        cfg.freeze_acc_bias_until_live = true;
        cfg.Racc_warmup_std = 0.5f;

        fusion_.begin(cfg);
        auto& filter = fusion_.raw();
        // Match the deployed/default test policy. The certificate runner does
        // not accept tuning ablations: its job is to test the implementation
        // that ships, not to find a more easily certifiable estimator.
        filter.setPeriodicAwCovarianceSync(true);
        filter.setAwCovarianceSyncCongruent(false);
        filter.enableLinearBlock(true);
        filter.enableTuner(true);
        filter.enableClamp(true);

        const char* trace_path = std::getenv("OU3_CERT_TRACE");
        if (!trace_path || !*trace_path) {
            throw std::runtime_error("OU3_CERT_TRACE must name a certificate trace CSV");
        }
        trace_.open(trace_path);
        if (!trace_) throw std::runtime_error("cannot open OU3_CERT_TRACE");
        trace_ << std::setprecision(9);
        write_header();
    }

    void updateMag(const Vector3f& mag_body_ned) override
    {
        fusion_.updateMag(mag_body_ned);
        if (fusion_.isLive() && fusion_.raw().mekf().lastMagDiag().accepted) {
            ++mag_accepted_since_trace_;
        }
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        auto& mekf = fusion_.raw().mekf();
        bool pseudo_due = false;
        if (mekf.linear_block_enabled()) {
            const float period = mekf.get_pseudo_update_period_s();
            pseudo_due = ocean_imu::kalman::ou_detail::periodic_update_due(
                dt, period, pseudo_elapsed_);
        } else {
            pseudo_elapsed_ = 0.0f;
        }
        if (pseudo_due) ++pseudo_since_trace_;

        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
        if (fusion_.isLive() && fusion_.raw().mekf().lastAccDiag().accepted) {
            ++acc_accepted_since_trace_;
        }
        time_s_ += dt;
        ++sample_index_;
    }

    FilterSnapshot snapshot() const override
    {
        const auto& filter = fusion_.raw();
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
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(mekf) + fusion_.magHardIronBodyUT();

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

        if (trace_ && sample_index_ % static_cast<unsigned>(trace_stride_) == 0u) {
            write_trace_row(q_bw_ned);
            acc_accepted_since_trace_ = 0;
            mag_accepted_since_trace_ = 0;
            pseudo_since_trace_ = 0;
        }
        return s;
    }

private:
    void write_header()
    {
        // The three *_accepted/pseudo columns are event COUNTS accumulated
        // since the previous trace row, not sticky last-diagnostic booleans.
        trace_ << "time_s,live,linear,bias_active,mag_lock,mag_refined,"
                  "acc_accepted,mag_accepted,pseudo_due_mirror,"
                  "tau_applied,sigma_applied,rs_applied,pseudo_period_s,"
                  "acc_nis,mag_nis,qw,qx,qy,qz,"
                  "bg_x,bg_y,bg_z,v_x,v_y,v_z,p_x,p_y,p_z,"
                  "S_x,S_y,S_z,aw_x,aw_y,aw_z,ba_x,ba_y,ba_z";
        for (int i = 0; i < 21; ++i) trace_ << ",Pdiag_" << i;
        // P(:,S) is the exact covariance information used by the S=0 gain.
        // The active implementation has 21 states and S occupies indices 12:14.
        for (int i = 0; i < 21; ++i)
            for (int j = 12; j < 15; ++j)
                trace_ << ",P_" << i << "_" << j;
        trace_ << "\n";
    }

    void write_vec(const Vector3f& v) const
    {
        trace_ << ',' << v.x() << ',' << v.y() << ',' << v.z();
    }

    void write_trace_row(const Quaternionf& q) const
    {
        const auto& filter = fusion_.raw();
        const auto& mekf = filter.mekf();
        const auto P = mekf.covariance_full();
        const auto& ad = mekf.lastAccDiag();
        const auto& md = mekf.lastMagDiag();

        trace_ << time_s_
               << ',' << (fusion_.isLive() ? 1 : 0)
               << ',' << (mekf.linear_block_enabled() ? 1 : 0)
               << ',' << (mekf.acc_bias_updates_enabled() ? 1 : 0)
               << ',' << (fusion_.hasMagNorthLock() ? 1 : 0)
               << ',' << (fusion_.hasRefinedMagReference() ? 1 : 0)
               << ',' << acc_accepted_since_trace_
               << ',' << mag_accepted_since_trace_
               << ',' << pseudo_since_trace_
               << ',' << filter.getTauApplied()
               << ',' << filter.getSigmaApplied()
               << ',' << filter.getRSApplied()
               << ',' << mekf.get_pseudo_update_period_s()
               << ',' << ad.nis << ',' << md.nis
               << ',' << q.w() << ',' << q.x() << ',' << q.y() << ',' << q.z();
        write_vec(mekf.gyroscope_bias());
        write_vec(mekf.get_velocity());
        write_vec(mekf.get_position());
        write_vec(mekf.get_integral_displacement());
        write_vec(mekf.get_world_accel());
        write_vec(mekf.get_acc_bias());
        for (int i = 0; i < 21; ++i) trace_ << ',' << P(i, i);
        for (int i = 0; i < 21; ++i)
            for (int j = 12; j < 15; ++j)
                trace_ << ',' << P(i, j);
        trace_ << '\n';
    }

    mutable Fusion fusion_;
    mutable std::ofstream trace_;
    mutable unsigned sample_index_ = 0;
    mutable float time_s_ = 0.0f;
    mutable float pseudo_elapsed_ = 0.0f;
    mutable int acc_accepted_since_trace_ = 0;
    mutable int mag_accepted_since_trace_ = 0;
    mutable int pseudo_since_trace_ = 0;
    int trace_stride_ = 10;
};

void process_one(const std::string& filename,
                 bool with_mag,
                 const W3dRandomSeeds& seeds,
                 bool write_timeseries,
                 float validation_window_sec)
{
    auto result = process_wave_file_for_tracker<CertificateAdapter>(
        filename, kDt, with_mag, true, kMagOdrHz,
        "_fusion_ou3_cert", "_fusion_ou3_cert_nomag", seeds, write_timeseries);
    if (!result) return;
    if (validation_window_sec > 0.0f)
        print_validation_metrics(*result, kDt, validation_window_sec, "OU_III_CERT");
    static constexpr W3dSummaryLabels labels{.target = "RS_target", .applied = "RS_applied"};
    print_summary_and_fail_if_needed(*result, kDt, Limits::value, labels);
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
            std::cout << "Usage: " << argv[0] << " [--nomag] [--input PATH]...\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown or incomplete argument: " << arg << "\n";
            return 2;
        }
    }
    if (files.empty()) files = collect_wave_data_files(".");

    W3dRandomSeeds seeds;
    try { seeds = w3d_random_seeds_from_env(); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    bool write_timeseries = true;
    if (const char* raw = std::getenv("W3D_WRITE_TIMESERIES"))
        write_timeseries = std::string(raw) != "0";
    float validation_window_sec = 900.0f;
    if (const char* raw = std::getenv("W3D_VALIDATION_WINDOW_SEC"))
        validation_window_sec = static_cast<float>(std::atof(raw));

    // A trace file corresponds to exactly one replay. The Python driver runs
    // the eight records one process at a time so provenance cannot be mixed.
    if (files.size() != 1u) {
        std::cerr << "ERROR: ou3-certificate-sim requires exactly one --input\n";
        return 2;
    }

    try { process_one(files.front(), with_mag, seeds, write_timeseries, validation_window_sec); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return w3d_any_quality_gate_failed() ? 1 : 0;
}
