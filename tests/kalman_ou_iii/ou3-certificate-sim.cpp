#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"

// Certificate-only visibility into the MEKF's already-computed scratch
// matrices.  This translation unit is a host test executable; production
// sources and estimator behavior are unchanged.  Include the MEKF once with
// private members visible, then include the wrapper normally (pragma once keeps
// the class definition from being repeated).
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

constexpr float kSigmaARescale = 0.71f;
constexpr float kSigmaGRescale = 0.05f;
constexpr float kSigmaMRescale = 2.0f;
constexpr float kMagOdrHz = 25.0f;
constexpr float kDt = 1.0f / 200.0f;
constexpr int kNX = 21;
constexpr int kOffV = 6;
constexpr int kOffS = 12;
constexpr int kOffAw = 15;
constexpr int kOffBa = 18;
using Matrix21f = Eigen::Matrix<float, kNX, kNX>;
using Matrix21x3f = Eigen::Matrix<float, kNX, 3>;
using Matrix3x21f = Eigen::Matrix<float, 3, kNX>;
using Vector12f = Eigen::Matrix<float, 12, 1>;

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

template<typename T>
void write_binary(std::ofstream& f, const T& value)
{
    f.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

Matrix21f reset_transport(const Vector3f& dtheta)
{
    Matrix21f G = Matrix21f::Identity();
    if (dtheta.allFinite() && dtheta.squaredNorm() > 0.0f) {
        G.block<3,3>(0,0) = Eigen::Matrix3f::Identity()
            + 0.5f * ocean_imu::kalman::ou_detail::skew(dtheta);
    }
    return G;
}

Matrix21f joseph_from_pct(const Matrix21f& P,
                          const Matrix21x3f& K,
                          const Eigen::Matrix3f& S,
                          const Matrix21x3f& PCt)
{
    Matrix21f out = P - K * PCt.transpose() - PCt * K.transpose()
                      + K * S * K.transpose();
    return 0.5f * (out + out.transpose());
}

class CertificateAdapter final : public IW3dFusionAdapter {
public:
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

    CertificateAdapter(bool with_mag,
                       const Vector3f& sigma_a_init,
                       const Vector3f& sigma_g,
                       const Vector3f& sigma_m)
        : trace_stride_(env_positive_int("OU3_CERT_TRACE_STRIDE", 10)),
          map_stride_(env_positive_int("OU3_CERT_MAP_STRIDE", 50))
    {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init * kSigmaARescale;
        cfg.sigma_g = sigma_g * kSigmaGRescale;
        cfg.sigma_m = sigma_m * kSigmaMRescale;
        cfg.mag_delay_sec = MAG_DELAY_SEC;

        fusion_.begin(cfg);
        auto& filter = fusion_.raw();
        // Match the deployed/default test policy.  The certificate executable
        // observes the shipping estimator; it does not retune it to pass.
        filter.setPeriodicAwCovarianceSync(true);
        filter.setAwCovarianceSyncCongruent(false);
        filter.enableTuner(true);
        filter.enableClamp(true);

        const char* trace_path = std::getenv("OU3_CERT_TRACE");
        if (!trace_path || !*trace_path)
            throw std::runtime_error("OU3_CERT_TRACE must name a certificate trace CSV");
        trace_.open(trace_path);
        if (!trace_) throw std::runtime_error("cannot open OU3_CERT_TRACE");
        trace_ << std::setprecision(9);
        write_header();

        const char* map_path = std::getenv("OU3_CERT_MAP_TRACE");
        if (!map_path || !*map_path)
            throw std::runtime_error("OU3_CERT_MAP_TRACE must name an exact-map binary trace");
        map_trace_.open(map_path, std::ios::binary);
        if (!map_trace_) throw std::runtime_error("cannot open OU3_CERT_MAP_TRACE");
        write_map_header();
    }

    ~CertificateAdapter() override
    {
        if (map_samples_ > 0) write_map_record();
    }

    void updateMag(const Vector3f& mag_body_ned) override
    {
        auto& mekf = fusion_.raw().mekf();
        ensure_map_block_started();
        const Matrix21f Ppre = mekf.covariance_full();
        const bool live_pre = fusion_.isLive();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refine_pre = fusion_.hasRefinedMagReference();

        fusion_.updateMag(mag_body_ned);

        const bool accepted = fusion_.raw().mekf().lastMagDiag().accepted;
        const bool lock_post = fusion_.hasMagNorthLock();
        const bool refine_post = fusion_.hasRefinedMagReference();
        if (lock_pre != lock_post || refine_pre != refine_post) {
            map_valid_ = false;
            map_hybrid_jump_ = true;
        }

        if (fusion_.isLive() && accepted) {
            ++mag_accepted_since_trace_;
            ++map_mag_count_;
            if (live_pre && active_pre == mekf.acc_bias_updates_enabled()
                    && lock_pre == lock_post && refine_pre == refine_post) {
                const Matrix21x3f PCt = mekf.PCt_scratch_;
                const Matrix21x3f K = mekf.K_scratch_;
                Eigen::LDLT<Matrix21f> ldlt(Ppre);
                if (ldlt.info() == Eigen::Success) {
                    const Matrix21x3f Ht = ldlt.solve(PCt);
                    const float denom = std::max(1.0f, PCt.norm());
                    map_linearization_residual_ = std::max(
                        map_linearization_residual_,
                        (Ppre * Ht - PCt).norm() / denom);
                    const Matrix3x21f H = Ht.transpose();
                    const Vector3f dtheta = K.topRows<3>() * mekf.lastMagDiag().r;
                    const Matrix21f C = reset_transport(dtheta)
                                      * (Matrix21f::Identity() - K * H);
                    map_accum_ = C * map_accum_;
                } else {
                    map_valid_ = false;
                }
            }
        }
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        auto& filter = fusion_.raw();
        auto& mekf = filter.mekf();
        ensure_map_block_started();

        const Matrix21f P0 = mekf.covariance_full();
        Vector12f xlin0;
        xlin0.segment<3>(0) = mekf.get_velocity();
        xlin0.segment<3>(3) = mekf.get_position();
        xlin0.segment<3>(6) = mekf.get_integral_displacement();
        xlin0.segment<3>(9) = mekf.get_world_accel();

        const bool live_pre = fusion_.isLive();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refine_pre = fusion_.hasRefinedMagReference();
        const bool pending_aw_floor = mekf.aw_covariance_floor_pending_;
        const Eigen::Matrix3f aw_floor_target = mekf.aw_covariance_floor_target_;
        const Eigen::Matrix3f R_S_pre = mekf.R_S;
        const float tau_b_pre = mekf.tau_bacc_;
        const Eigen::Matrix3f Q_bacc_pre = mekf.Q_bacc_;

        float pseudo_elapsed_copy = mekf.pseudo_update_elapsed_s_;
        const bool pseudo_due =
            ocean_imu::kalman::ou_detail::periodic_update_due(
                dt, mekf.pseudo_update_period_s_, pseudo_elapsed_copy);
        if (pseudo_due) {
            ++pseudo_since_trace_;
            ++map_pseudo_count_;
        }

        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);

        const bool live_post = fusion_.isLive();
        const bool active_post = mekf.acc_bias_updates_enabled();
        const bool lock_post = fusion_.hasMagNorthLock();
        const bool refine_post = fusion_.hasRefinedMagReference();

        if (live_pre != live_post
                || active_pre != active_post || lock_pre != lock_post
                || refine_pre != refine_post) {
            map_valid_ = false;
            map_hybrid_jump_ = true;
        }

        // The prediction matrix is not inferred from the trajectory.  These are
        // the exact matrices the MEKF just used in time_update().
        Matrix21f A = Matrix21f::Identity();
        A.block<6,6>(0,0) = mekf.F_AA_scratch_;
        A.block<12,12>(kOffV,kOffV) = mekf.F_LL_scratch_;
        const float tau_b = std::max(1e-3f, tau_b_pre);
        const float phi_b = active_pre ? std::exp(-dt / tau_b) : 1.0f;
        A.block<3,3>(kOffBa,kOffBa) = phi_b * Eigen::Matrix3f::Identity();

        Matrix21f Q = Matrix21f::Zero();
        Q.block<6,6>(0,0) = mekf.Q_AA_scratch_;
        Q.block<12,12>(kOffV,kOffV) = mekf.Q_LL_scratch_;
        if (active_pre) {
            const float qd_scale = -0.5f * tau_b * std::expm1(-2.0f * dt / tau_b);
            Q.block<3,3>(kOffBa,kOffBa) = Q_bacc_pre * qd_scale;
        }

        Matrix21f Pcur = A * P0 * A.transpose() + Q;
        Pcur = 0.5f * (Pcur + Pcur.transpose());

        // The proof-compatible covariance synchronization is a PSD increment
        // applied inside prediction. Reproduce it exactly because it changes
        // the subsequent S=0 gain, although it does not itself change the
        // deterministic error-state transition A.
        if (pending_aw_floor) {
            Eigen::Matrix3f Delta = aw_floor_target - Pcur.block<3,3>(kOffAw,kOffAw);
            Delta = 0.5f * (Delta + Delta.transpose());
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(Delta);
            if (es.info() == Eigen::Success) {
                Vector3f evals = es.eigenvalues().cwiseMax(0.0f);
                const Eigen::Matrix3f DeltaPlus =
                    es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
                Pcur.block<3,3>(kOffAw,kOffAw) += DeltaPlus;
            } else {
                map_valid_ = false;
            }
        }

        Matrix21f sample_map = A;
        const Vector12f xlin_pred = mekf.F_LL_scratch_ * xlin0;

        if (pseudo_due) {
            Matrix3x21f Hs = Matrix3x21f::Zero();
            Hs.block<3,3>(0,kOffS) = Eigen::Matrix3f::Identity();
            Matrix21x3f PCt = Pcur.block<kNX,3>(0,kOffS);
            if (!active_pre) PCt.block<3,3>(kOffBa,0).setZero();
            const Eigen::Matrix3f S = Pcur.block<3,3>(kOffS,kOffS) + R_S_pre;
            Eigen::LDLT<Eigen::Matrix3f> ldlt(S);
            if (ldlt.info() == Eigen::Success) {
                const Matrix21x3f K = ldlt.solve(PCt.transpose()).transpose();
                const Vector3f rS = -xlin_pred.segment<3>(6);
                const Vector3f dtheta = K.topRows<3>() * rS;
                const Matrix21f G = reset_transport(dtheta);
                const Matrix21f C = G * (Matrix21f::Identity() - K * Hs);
                sample_map = C * sample_map;
                Pcur = joseph_from_pct(Pcur, K, S, PCt);
                Pcur = G * Pcur * G.transpose();
                Pcur = 0.5f * (Pcur + Pcur.transpose());
            } else {
                map_valid_ = false;
            }
        }

        const auto& ad = mekf.lastAccDiag();
        if (live_post && ad.accepted) {
            ++acc_accepted_since_trace_;
            ++map_acc_count_;
        }

        // After fusion_.update(), K_scratch_/PCt_scratch_ belong to the
        // accelerometer correction (the last MEKF measurement in update()).
        // Recover H from the exact pre-accelerometer covariance and PCt=P H^T;
        // this avoids reimplementing the nonlinear measurement Jacobian.
        if (ad.accepted && live_pre && live_post) {
            const Matrix21x3f PCt = mekf.PCt_scratch_;
            const Matrix21x3f K = mekf.K_scratch_;
            Eigen::LDLT<Matrix21f> ldlt(Pcur);
            if (ldlt.info() == Eigen::Success) {
                const Matrix21x3f Ht = ldlt.solve(PCt);
                const float denom = std::max(1.0f, PCt.norm());
                map_linearization_residual_ = std::max(
                    map_linearization_residual_,
                    (Pcur * Ht - PCt).norm() / denom);
                const Matrix3x21f H = Ht.transpose();
                const Vector3f dtheta = K.topRows<3>() * ad.r;
                const Matrix21f C = reset_transport(dtheta)
                                  * (Matrix21f::Identity() - K * H);
                sample_map = C * sample_map;
            } else {
                map_valid_ = false;
            }
        }

        if (!(live_pre && live_post)) map_valid_ = false;
        map_accum_ = sample_map * map_accum_;

        time_s_ += dt;
        ++sample_index_;
        ++map_samples_;
        if (map_samples_ >= static_cast<unsigned>(map_stride_)) write_map_record();
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
        trace_ << "time_s,live,linear,bias_active,mag_lock,mag_refined,"
                  "acc_accepted,mag_accepted,pseudo_due_mirror,"
                  "tau_applied,sigma_applied,rs_applied,pseudo_period_s,"
                  "acc_nis,mag_nis,qw,qx,qy,qz,"
                  "bg_x,bg_y,bg_z,v_x,v_y,v_z,p_x,p_y,p_z,"
                  "S_x,S_y,S_z,aw_x,aw_y,aw_z,ba_x,ba_y,ba_z";
        for (int i = 0; i < 21; ++i) trace_ << ",Pdiag_" << i;
        for (int i = 0; i < 21; ++i)
            for (int j = 12; j < 15; ++j)
                trace_ << ",P_" << i << "_" << j;
        trace_ << "\n";
    }

    void write_map_header()
    {
        const char magic[8] = {'O','U','3','M','A','P','3','\0'};
        map_trace_.write(magic, sizeof(magic));
        const std::uint32_t version = 3;
        const std::uint32_t nx = kNX;
        const std::uint32_t stride = static_cast<std::uint32_t>(map_stride_);
        write_binary(map_trace_, version);
        write_binary(map_trace_, nx);
        write_binary(map_trace_, stride);
    }

    void ensure_map_block_started()
    {
        if (map_block_started_) return;
        const auto& filter = fusion_.raw();
        const auto& mekf = filter.mekf();
        map_block_started_ = true;
        map_start_time_ = time_s_;
        map_start_live_ = fusion_.isLive();
        map_start_active_ = mekf.acc_bias_updates_enabled();
        map_start_lock_ = fusion_.hasMagNorthLock();
        map_start_refined_ = fusion_.hasRefinedMagReference();
        map_start_tau_ = filter.getTauApplied();
        map_start_sigma_ = filter.getSigmaApplied();
        map_start_rs_ = filter.getRSApplied();
        map_valid_ = map_start_live_;
        map_hybrid_jump_ = false;
        map_linearization_residual_ = 0.0f;
        map_accum_.setIdentity();
        map_acc_count_ = map_mag_count_ = map_pseudo_count_ = 0;
    }

    void write_map_record()
    {
        if (!map_block_started_) return;
        const auto& filter = fusion_.raw();
        const auto& mekf = filter.mekf();
        const bool end_live = fusion_.isLive();
        const bool end_active = mekf.acc_bias_updates_enabled();
        const bool end_lock = fusion_.hasMagNorthLock();
        const bool end_refined = fusion_.hasRefinedMagReference();
        if (map_start_active_ != end_active || map_start_live_ != end_live
                || map_start_lock_ != end_lock || map_start_refined_ != end_refined)
            map_valid_ = false;

        std::uint32_t flags = 0;
        if (map_valid_) flags |= 1u << 0;
        if (map_start_live_) flags |= 1u << 1;
        if (end_live) flags |= 1u << 2;
        if (map_start_active_) flags |= 1u << 3;
        if (end_active) flags |= 1u << 4;
        if (map_start_lock_) flags |= 1u << 5;
        if (end_lock) flags |= 1u << 6;
        if (map_hybrid_jump_) flags |= 1u << 7;
        if (map_start_refined_) flags |= 1u << 8;
        if (end_refined) flags |= 1u << 9;

        const double t0 = static_cast<double>(map_start_time_);
        const double t1 = static_cast<double>(time_s_);
        write_binary(map_trace_, t0);
        write_binary(map_trace_, t1);
        write_binary(map_trace_, flags);
        write_binary(map_trace_, static_cast<std::int32_t>(map_acc_count_));
        write_binary(map_trace_, static_cast<std::int32_t>(map_mag_count_));
        write_binary(map_trace_, static_cast<std::int32_t>(map_pseudo_count_));
        write_binary(map_trace_, map_start_tau_);
        write_binary(map_trace_, map_start_sigma_);
        write_binary(map_trace_, map_start_rs_);
        const float tau1 = filter.getTauApplied();
        const float sigma1 = filter.getSigmaApplied();
        const float rs1 = filter.getRSApplied();
        write_binary(map_trace_, tau1);
        write_binary(map_trace_, sigma1);
        write_binary(map_trace_, rs1);
        write_binary(map_trace_, map_linearization_residual_);
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(map_trace_, map_accum_(i,j));

        map_block_started_ = false;
        map_samples_ = 0;
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
               << ",1"   // linear block: always on
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
    std::ofstream map_trace_;
    mutable unsigned sample_index_ = 0;
    mutable float time_s_ = 0.0f;
    mutable int acc_accepted_since_trace_ = 0;
    mutable int mag_accepted_since_trace_ = 0;
    mutable int pseudo_since_trace_ = 0;
    int trace_stride_ = 10;
    int map_stride_ = 50;

    bool map_block_started_ = false;
    bool map_valid_ = false;
    bool map_hybrid_jump_ = false;
    bool map_start_live_ = false;
    bool map_start_active_ = false;
    bool map_start_lock_ = false;
    bool map_start_refined_ = false;
    unsigned map_samples_ = 0;
    int map_acc_count_ = 0;
    int map_mag_count_ = 0;
    int map_pseudo_count_ = 0;
    float map_start_time_ = 0.0f;
    float map_start_tau_ = 0.0f;
    float map_start_sigma_ = 0.0f;
    float map_start_rs_ = 0.0f;
    float map_linearization_residual_ = 0.0f;
    Matrix21f map_accum_ = Matrix21f::Identity();
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
