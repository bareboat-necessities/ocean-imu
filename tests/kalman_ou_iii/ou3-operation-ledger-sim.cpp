#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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
#include "kalman_ou_common/KalmanOUCoreMath.h"

// Host-only complete-word observer. Production visibility/behavior is unchanged.
// One shipping estimator owns both point-map scanning and selected-direction
// event ledgers, so current-sample schedule data are read at exactly one
// shipping boundary. No second Riccati recursion or perturbed filter exists.
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

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
constexpr int kOffV = 6;
constexpr int kOffS = 12;
constexpr int kOffAw = 15;
constexpr int kOffBa = 18;
using Vector21f = Eigen::Matrix<float, kNX, 1>;
using Matrix21f = Eigen::Matrix<float, kNX, kNX>;
using Matrix21x3f = Eigen::Matrix<float, kNX, 3>;
using Matrix3x21f = Eigen::Matrix<float, 3, kNX>;
using Vector12f = Eigen::Matrix<float, 12, 1>;

template<typename T>
void write_binary(std::ofstream& f, const T& value)
{
    f.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

int env_positive_int(const char* name, int fallback)
{
    if (const char* raw = std::getenv(name)) {
        const int value = std::atoi(raw);
        if (value > 0) return value;
    }
    return fallback;
}

float env_float_required(const char* name)
{
    const char* raw = std::getenv(name);
    if (!raw || !*raw) throw std::runtime_error(std::string(name) + " is required");
    char* end = nullptr;
    const float x = std::strtof(raw, &end);
    if (end == raw || !std::isfinite(x))
        throw std::runtime_error(std::string(name) + " must be finite numeric");
    return x;
}

std::string env_string_required(const char* name)
{
    const char* raw = std::getenv(name);
    if (!raw || !*raw) throw std::runtime_error(std::string(name) + " is required");
    return std::string(raw);
}

Vector21f parse_direction()
{
    Vector21f d = Vector21f::Zero();
    const std::string raw = env_string_required("OU3_LEDGER_DIRECTION");
    std::stringstream ss(raw);
    std::string token;
    int i = 0;
    while (std::getline(ss, token, ',')) {
        if (i >= kNX) throw std::runtime_error("OU3_LEDGER_DIRECTION has more than 21 entries");
        char* end = nullptr;
        const float x = std::strtof(token.c_str(), &end);
        if (end == token.c_str() || !std::isfinite(x))
            throw std::runtime_error("OU3_LEDGER_DIRECTION contains invalid value");
        d(i++) = x;
    }
    if (i != kNX) throw std::runtime_error("OU3_LEDGER_DIRECTION must contain exactly 21 entries");
    return d;
}

Matrix21f reset_transport(const Vector3f& dtheta)
{
    Matrix21f G = Matrix21f::Identity();
    if (dtheta.allFinite() && dtheta.squaredNorm() > 0.0f) {
        G.block<3,3>(0,0) = Matrix3f::Identity()
            + 0.5f * ocean_imu::kalman::ou_detail::skew(dtheta);
    }
    return G;
}

Matrix21f joseph_from_pct(const Matrix21f& P,
                          const Matrix21x3f& K,
                          const Matrix3f& S,
                          const Matrix21x3f& PCt)
{
    Matrix21f out = P - K * PCt.transpose() - PCt * K.transpose()
                      + K * S * K.transpose();
    return 0.5f * (out + out.transpose());
}

float relative_matrix_difference(const Matrix21f& a, const Matrix21f& b)
{
    return (a - b).norm() / std::max(1.0e-12f, b.norm());
}

class OperationObserverAdapter final : public IW3dFusionAdapter {
public:
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

    OperationObserverAdapter(bool with_mag,
                             const Vector3f& sigma_a_init,
                             const Vector3f& sigma_g,
                             const Vector3f& sigma_m)
    {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init * kSigmaARescale;
        cfg.sigma_g = sigma_g * kSigmaGRescale;
        cfg.sigma_m = sigma_m * kSigmaMRescale;
        cfg.mag_delay_sec = MAG_DELAY_SEC;
        fusion_.begin(cfg);
        auto& filter = fusion_.raw();
        filter.setPeriodicAwCovarianceSync(true);
        filter.setAwCovarianceSyncCongruent(false);
        filter.enableTuner(true);
        filter.enableClamp(true);

        const char* map_path = std::getenv("OU3_LEDGER_MAP_TRACE");
        scan_mode_ = map_path && *map_path;
        if (scan_mode_) {
            scan_stride_ = env_positive_int("OU3_LEDGER_MAP_STRIDE", 600);
            const std::string cov_path = env_string_required("OU3_LEDGER_COV_TRACE");
            map_trace_.open(map_path, std::ios::binary);
            cov_trace_.open(cov_path, std::ios::binary);
            if (!map_trace_ || !cov_trace_)
                throw std::runtime_error("cannot open operation-observer map/covariance traces");
            write_scan_headers();
            return;
        }

        requested_t0_ = env_float_required("OU3_LEDGER_T0");
        requested_t1_ = env_float_required("OU3_LEDGER_T1");
        requested_mode_ = env_string_required("OU3_LEDGER_MODE");
        direction_ = parse_direction();
        if (!(requested_t1_ > requested_t0_))
            throw std::runtime_error("OU3_LEDGER_T1 must exceed OU3_LEDGER_T0");
        if (requested_mode_ != "H18" && requested_mode_ != "A21")
            throw std::runtime_error("OU3_LEDGER_MODE must be H18 or A21");
        if (requested_mode_ == "H18" && direction_.tail<3>().norm() > 1.0e-12f)
            throw std::runtime_error("H18 ledger direction must have zero b_a tail");

        const std::string path = env_string_required("OU3_LEDGER_TRACE");
        trace_.open(path);
        if (!trace_) throw std::runtime_error("cannot open OU3_LEDGER_TRACE");
        trace_ << std::setprecision(12);
        trace_ << "sequence,time_s,event,mode,V_before,V_after,delta_V,event_ratio,"
                  "rs_scalar,rs_period_s,rs_std_x,rs_std_y,rs_std_z,"
                  "rs_var_x,rs_var_y,rs_var_z,dtheta_norm,covariance_reconstruction_rel\n";
    }

    ~OperationObserverAdapter() override
    {
        if (scan_mode_ && scan_block_started_ && scan_samples_ > 0) write_scan_record();
    }

    void updateMag(const Vector3f& mag_body_ned) override
    {
        if (scan_mode_) {
            update_mag_scan(mag_body_ned);
            return;
        }
        update_mag_selected(mag_body_ned);
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        if (scan_mode_) {
            update_scan(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
            return;
        }
        update_selected(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
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
        return s;
    }

private:
    int dimension() const { return requested_mode_ == "H18" ? 18 : 21; }

    void write_scan_headers()
    {
        const char map_magic[8] = {'O','U','3','M','A','P','3','\0'};
        map_trace_.write(map_magic, sizeof(map_magic));
        const std::uint32_t map_version = 3;
        const std::uint32_t nx = kNX;
        const std::uint32_t stride = static_cast<std::uint32_t>(scan_stride_);
        write_binary(map_trace_, map_version);
        write_binary(map_trace_, nx);
        write_binary(map_trace_, stride);

        const char cov_magic[8] = {'O','U','3','C','O','V','1','\0'};
        cov_trace_.write(cov_magic, sizeof(cov_magic));
        const std::uint32_t cov_version = 1;
        write_binary(cov_trace_, cov_version);
        write_binary(cov_trace_, nx);
        write_binary(cov_trace_, stride);
    }

    void ensure_scan_block_started()
    {
        if (scan_block_started_) return;
        const auto& filter = fusion_.raw();
        const auto& mekf = filter.mekf();
        scan_block_started_ = true;
        scan_start_time_ = time_s_;
        scan_start_live_ = fusion_.isLive();
        scan_start_active_ = mekf.acc_bias_updates_enabled();
        scan_start_lock_ = fusion_.hasMagNorthLock();
        scan_start_refined_ = fusion_.hasRefinedMagReference();
        scan_start_tau_ = filter.getTauApplied();
        scan_start_sigma_ = filter.getSigmaApplied();
        scan_start_rs_ = filter.getRSApplied();
        scan_cov_start_ = mekf.covariance_full();
        scan_valid_ = scan_start_live_;
        scan_hybrid_jump_ = false;
        scan_linearization_residual_ = 0.0f;
        scan_accum_.setIdentity();
        scan_samples_ = 0;
        scan_acc_count_ = scan_mag_count_ = scan_pseudo_count_ = 0;
    }

    void invalidate_scan(bool hybrid = false)
    {
        scan_valid_ = false;
        scan_hybrid_jump_ = scan_hybrid_jump_ || hybrid;
    }

    void write_scan_record()
    {
        if (!scan_block_started_) return;
        const auto& filter = fusion_.raw();
        const auto& mekf = filter.mekf();
        const bool end_live = fusion_.isLive();
        const bool end_active = mekf.acc_bias_updates_enabled();
        const bool end_lock = fusion_.hasMagNorthLock();
        const bool end_refined = fusion_.hasRefinedMagReference();
        if (scan_start_active_ != end_active || scan_start_live_ != end_live
                || scan_start_lock_ != end_lock || scan_start_refined_ != end_refined)
            scan_valid_ = false;

        std::uint32_t flags = 0;
        if (scan_valid_) flags |= 1u << 0;
        if (scan_start_live_) flags |= 1u << 1;
        if (end_live) flags |= 1u << 2;
        if (scan_start_active_) flags |= 1u << 3;
        if (end_active) flags |= 1u << 4;
        if (scan_start_lock_) flags |= 1u << 5;
        if (end_lock) flags |= 1u << 6;
        if (scan_hybrid_jump_) flags |= 1u << 7;
        if (scan_start_refined_) flags |= 1u << 8;
        if (end_refined) flags |= 1u << 9;

        const double t0 = static_cast<double>(scan_start_time_);
        const double t1 = static_cast<double>(time_s_);
        write_binary(map_trace_, t0);
        write_binary(map_trace_, t1);
        write_binary(map_trace_, flags);
        write_binary(map_trace_, static_cast<std::int32_t>(scan_acc_count_));
        write_binary(map_trace_, static_cast<std::int32_t>(scan_mag_count_));
        write_binary(map_trace_, static_cast<std::int32_t>(scan_pseudo_count_));
        write_binary(map_trace_, scan_start_tau_);
        write_binary(map_trace_, scan_start_sigma_);
        write_binary(map_trace_, scan_start_rs_);
        const float tau1 = filter.getTauApplied();
        const float sigma1 = filter.getSigmaApplied();
        const float rs1 = filter.getRSApplied();
        write_binary(map_trace_, tau1);
        write_binary(map_trace_, sigma1);
        write_binary(map_trace_, rs1);
        write_binary(map_trace_, scan_linearization_residual_);
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(map_trace_, scan_accum_(i,j));

        write_binary(cov_trace_, t0);
        write_binary(cov_trace_, t1);
        const Matrix21f cov_end = mekf.covariance_full();
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(cov_trace_, scan_cov_start_(i,j));
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(cov_trace_, cov_end(i,j));

        scan_block_started_ = false;
        scan_samples_ = 0;
    }

    void scan_linear_event(const Matrix21f& C)
    {
        scan_accum_ = C * scan_accum_;
    }

    void update_mag_scan(const Vector3f& mag_body_ned)
    {
        ensure_scan_block_started();
        auto& mekf = fusion_.raw().mekf();
        const Matrix21f Ppre = mekf.covariance_full();
        const bool live_pre = fusion_.isLive();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refined_pre = fusion_.hasRefinedMagReference();

        fusion_.updateMag(mag_body_ned);

        const bool active_post = mekf.acc_bias_updates_enabled();
        const bool lock_post = fusion_.hasMagNorthLock();
        const bool refined_post = fusion_.hasRefinedMagReference();
        if (active_pre != active_post || lock_pre != lock_post || refined_pre != refined_post)
            invalidate_scan(true);

        const auto& md = mekf.lastMagDiag();
        if (!(fusion_.isLive() && md.accepted)) return;
        ++scan_mag_count_;
        if (!(live_pre && active_pre == active_post
                && lock_pre == lock_post && refined_pre == refined_post)) return;

        const Matrix21x3f PCt = mekf.PCt_scratch_;
        const Matrix21x3f K = mekf.K_scratch_;
        Eigen::LDLT<Matrix21f> ldlt(Ppre);
        if (ldlt.info() != Eigen::Success) {
            invalidate_scan();
            return;
        }
        const Matrix21x3f Ht = ldlt.solve(PCt);
        const float resid = (Ppre * Ht - PCt).norm() / std::max(1.0f, PCt.norm());
        scan_linearization_residual_ = std::max(scan_linearization_residual_, resid);
        const Matrix3x21f H = Ht.transpose();
        const Vector3f dtheta = K.topRows<3>() * md.r;
        scan_linear_event(reset_transport(dtheta) * (Matrix21f::Identity() - K * H));
    }

    void update_scan(float dt,
                     const Vector3f& gyr_meas_ned,
                     const Vector3f& acc_meas_ned,
                     float temperature_c)
    {
        ensure_scan_block_started();
        auto& filter = fusion_.raw();

        // This is exactly the first stateful shipping action in updateCore_.
        // Consume it once before observing the current Riccati schedule; the
        // subsequent wrapper call sees no pending commit and is trajectory-equivalent.
        filter.apply_pending_online_tune_();
        auto& mekf = filter.mekf();

        const Matrix21f P0 = mekf.covariance_full();
        Vector12f xlin0;
        xlin0.segment<3>(0) = mekf.get_velocity();
        xlin0.segment<3>(3) = mekf.get_position();
        xlin0.segment<3>(6) = mekf.get_integral_displacement();
        xlin0.segment<3>(9) = mekf.get_world_accel();

        const bool live_pre = fusion_.isLive();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refined_pre = fusion_.hasRefinedMagReference();
        const bool pending_aw_floor = mekf.aw_covariance_floor_pending_;
        const Matrix3f aw_floor_target = mekf.aw_covariance_floor_target_;
        const Matrix3f R_S_pre = mekf.R_S;
        const float pseudo_period_pre = mekf.pseudo_update_period_s_;
        const float tau_b_pre = mekf.tau_bacc_;
        const Matrix3f Q_bacc_pre = mekf.Q_bacc_;
        float pseudo_elapsed_copy = mekf.pseudo_update_elapsed_s_;
        const bool pseudo_due = ocean_imu::kalman::ou_detail::periodic_update_due(
            dt, pseudo_period_pre, pseudo_elapsed_copy);
        if (pseudo_due) ++scan_pseudo_count_;

        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);

        const bool live_post = fusion_.isLive();
        const bool active_post = mekf.acc_bias_updates_enabled();
        const bool lock_post = fusion_.hasMagNorthLock();
        const bool refined_post = fusion_.hasRefinedMagReference();
        if (live_pre != live_post || active_pre != active_post
                || lock_pre != lock_post || refined_pre != refined_post)
            invalidate_scan(true);

        Matrix21f A = Matrix21f::Identity();
        A.block<6,6>(0,0) = mekf.F_AA_scratch_;
        A.block<12,12>(kOffV,kOffV) = mekf.F_LL_scratch_;
        const float tau_b = std::max(1.0e-3f, tau_b_pre);
        const float phi_b = active_pre ? std::exp(-dt / tau_b) : 1.0f;
        A.block<3,3>(kOffBa,kOffBa) = phi_b * Matrix3f::Identity();

        Matrix21f Q = Matrix21f::Zero();
        Q.block<6,6>(0,0) = mekf.Q_AA_scratch_;
        Q.block<12,12>(kOffV,kOffV) = mekf.Q_LL_scratch_;
        if (active_pre) {
            const float qd_scale = -0.5f * tau_b * std::expm1(-2.0f * dt / tau_b);
            Q.block<3,3>(kOffBa,kOffBa) = Q_bacc_pre * qd_scale;
        }

        Matrix21f Pcur = A * P0 * A.transpose() + Q;
        Pcur = 0.5f * (Pcur + Pcur.transpose());
        scan_linear_event(A);

        if (pending_aw_floor) {
            Matrix3f Delta = aw_floor_target - Pcur.block<3,3>(kOffAw,kOffAw);
            Delta = 0.5f * (Delta + Delta.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix3f> es(Delta);
            if (es.info() != Eigen::Success) {
                invalidate_scan();
            } else {
                const Vector3f evals = es.eigenvalues().cwiseMax(0.0f);
                const Matrix3f DeltaPlus =
                    es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
                Pcur.block<3,3>(kOffAw,kOffAw) += DeltaPlus;
            }
        }

        const Vector12f xlin_pred = mekf.F_LL_scratch_ * xlin0;
        if (pseudo_due) {
            Matrix3x21f Hs = Matrix3x21f::Zero();
            Hs.block<3,3>(0,kOffS) = Matrix3f::Identity();
            Matrix21x3f PCt = Pcur.block<kNX,3>(0,kOffS);
            if (!active_pre) PCt.block<3,3>(kOffBa,0).setZero();
            const Matrix3f S = Pcur.block<3,3>(kOffS,kOffS) + R_S_pre;
            Eigen::LDLT<Matrix3f> ldlt(S);
            if (ldlt.info() != Eigen::Success) {
                invalidate_scan();
            } else {
                const Matrix21x3f K = ldlt.solve(PCt.transpose()).transpose();
                const Vector3f rS = -xlin_pred.segment<3>(6);
                const Vector3f dtheta = K.topRows<3>() * rS;
                const Matrix21f G = reset_transport(dtheta);
                scan_linear_event(G * (Matrix21f::Identity() - K * Hs));
                Pcur = joseph_from_pct(Pcur, K, S, PCt);
                Pcur = G * Pcur * G.transpose();
                Pcur = 0.5f * (Pcur + Pcur.transpose());
            }
        }

        const auto& ad = mekf.lastAccDiag();
        if (live_post && ad.accepted) ++scan_acc_count_;
        if (ad.accepted && live_pre && live_post) {
            const Matrix21x3f PCt = mekf.PCt_scratch_;
            const Matrix21x3f K = mekf.K_scratch_;
            Eigen::LDLT<Matrix21f> ldlt(Pcur);
            if (ldlt.info() != Eigen::Success) {
                invalidate_scan();
            } else {
                const Matrix21x3f Ht = ldlt.solve(PCt);
                const float lin_resid =
                    (Pcur * Ht - PCt).norm() / std::max(1.0f, PCt.norm());
                scan_linearization_residual_ = std::max(
                    scan_linearization_residual_, lin_resid);
                const Matrix3x21f H = Ht.transpose();
                const Vector3f dtheta = K.topRows<3>() * ad.r;
                const Matrix21f G = reset_transport(dtheta);
                scan_linear_event(G * (Matrix21f::Identity() - K * H));

                Matrix21f Precon = joseph_from_pct(Pcur, K, mekf.S_scratch_, PCt);
                Precon = G * Precon * G.transpose();
                Precon = 0.5f * (Precon + Precon.transpose());
                const float cov_resid = relative_matrix_difference(
                    Precon, mekf.covariance_full());
                if (!(std::isfinite(cov_resid) && cov_resid <= 2.0e-4f))
                    invalidate_scan();
            }
        }

        if (!(live_pre && live_post)) scan_valid_ = false;
        time_s_ += dt;
        ++scan_samples_;
        if (scan_samples_ >= static_cast<unsigned>(scan_stride_)) write_scan_record();
    }

    bool inside_requested_word() const
    {
        return started_ && !finished_ && time_s_ < requested_t1_ - 2.0e-6f;
    }

    void check_mode()
    {
        const auto& m = fusion_.raw().mekf();
        const std::string mode = m.acc_bias_updates_enabled() ? "A21" : "H18";
        if (!fusion_.isLive() || mode != requested_mode_)
            fail("requested complete-word mode not preserved");
    }

    void maybe_start()
    {
        if (started_ || finished_) return;
        if (time_s_ + 2.0e-6f < requested_t0_) return;
        if (std::abs(time_s_ - requested_t0_) > 2.0e-4f)
            fail("ledger could not align to requested word start");
        check_mode();
        started_ = true;
        actual_t0_ = time_s_;
        error_ = direction_;
        initial_energy_ = energy(fusion_.raw().mekf().covariance_full(), error_);
        if (!(std::isfinite(initial_energy_) && initial_energy_ > 0.0f))
            fail("invalid initial information energy");
        std::cout << "OU3_EVENT_LEDGER_START mode=" << requested_mode_
                  << " t0=" << actual_t0_ << " V0=" << initial_energy_ << "\n";
    }

    float energy(const Matrix21f& Pfull, const Vector21f& xfull) const
    {
        if (dimension() == 18) {
            Eigen::Matrix<float,18,18> P = Pfull.topLeftCorner<18,18>();
            P = 0.5f * (P + P.transpose());
            Eigen::LDLT<Eigen::Matrix<float,18,18>> ldlt(P);
            if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
            const Eigen::Matrix<float,18,1> x = xfull.head<18>();
            return x.dot(ldlt.solve(x));
        }
        Matrix21f P = 0.5f * (Pfull + Pfull.transpose());
        Eigen::LDLT<Matrix21f> ldlt(P);
        if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
        return xfull.dot(ldlt.solve(xfull));
    }

    void write_event(const char* event,
                     float event_time,
                     const Matrix21f& Pbefore,
                     const Matrix21f& Pafter,
                     const Matrix21f& C,
                     float rs_scalar,
                     const Matrix3f& R_S,
                     float dtheta_norm,
                     float covariance_reconstruction_rel)
    {
        if (failed_ || !started_ || finished_) return;
        const float vb = energy(Pbefore, error_);
        const Vector21f xafter = C * error_;
        const float va = energy(Pafter, xafter);
        if (!(std::isfinite(vb) && vb > 0.0f && std::isfinite(va) && va >= 0.0f))
            fail(std::string("invalid information energy at ") + event);
        const float dv = va - vb;
        const float ratio = va / vb;
        const float rsx = R_S(0,0) >= 0.0f ? std::sqrt(R_S(0,0)) : std::numeric_limits<float>::quiet_NaN();
        const float rsy = R_S(1,1) >= 0.0f ? std::sqrt(R_S(1,1)) : std::numeric_limits<float>::quiet_NaN();
        const float rsz = R_S(2,2) >= 0.0f ? std::sqrt(R_S(2,2)) : std::numeric_limits<float>::quiet_NaN();
        trace_ << sequence_++ << ',' << event_time << ',' << event << ',' << requested_mode_
               << ',' << vb << ',' << va << ',' << dv << ',' << ratio
               << ',' << rs_scalar << ',' << fusion_.raw().mekf().pseudo_update_period_s_
               << ',' << rsx << ',' << rsy << ',' << rsz
               << ',' << R_S(0,0) << ',' << R_S(1,1) << ',' << R_S(2,2)
               << ',' << dtheta_norm << ',' << covariance_reconstruction_rel << '\n';
        error_ = xafter;
        const std::string e(event);
        if (e == "prediction") prediction_delta_ += dv;
        else if (e == "aw_floor") floor_delta_ += dv;
        else if (e == "S_zero") s_delta_ += dv;
        else if (e == "accelerometer") acc_delta_ += dv;
        else if (e == "vector") vector_delta_ += dv;
    }

    void update_mag_selected(const Vector3f& mag_body_ned)
    {
        if (finished_) return;
        maybe_start();
        if (!started_) {
            fusion_.updateMag(mag_body_ned);
            return;
        }
        if (!inside_requested_word()) return;

        auto& mekf = fusion_.raw().mekf();
        const Matrix21f Ppre = mekf.covariance_full();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refined_pre = fusion_.hasRefinedMagReference();
        fusion_.updateMag(mag_body_ned);
        if (active_pre != mekf.acc_bias_updates_enabled()
                || lock_pre != fusion_.hasMagNorthLock()
                || refined_pre != fusion_.hasRefinedMagReference())
            fail("hybrid magnetic/mode event inside selected same-mode word");
        check_mode();

        const auto& md = mekf.lastMagDiag();
        if (!md.accepted) return;
        const Matrix21x3f PCt = mekf.PCt_scratch_;
        const Matrix21x3f K = mekf.K_scratch_;
        Eigen::LDLT<Matrix21f> ldlt(Ppre);
        if (ldlt.info() != Eigen::Success) fail("mag pre-covariance LDLT failed");
        const Matrix21x3f Ht = ldlt.solve(PCt);
        const Matrix3x21f H = Ht.transpose();
        const Vector3f dtheta = K.topRows<3>() * md.r;
        const Matrix21f C = reset_transport(dtheta)
                          * (Matrix21f::Identity() - K * H);
        write_event("vector", time_s_, Ppre, mekf.covariance_full(), C,
                    std::numeric_limits<float>::quiet_NaN(), Matrix3f::Zero(),
                    dtheta.norm(), 0.0f);
        ++vector_count_;
    }

    void update_selected(float dt,
                         const Vector3f& gyr_meas_ned,
                         const Vector3f& acc_meas_ned,
                         float temperature_c)
    {
        if (finished_) {
            time_s_ += dt;
            return;
        }
        maybe_start();
        if (!started_) {
            fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
            time_s_ += dt;
            return;
        }
        if (!inside_requested_word()) {
            time_s_ += dt;
            return;
        }

        auto& filter = fusion_.raw();
        filter.apply_pending_online_tune_();
        auto& mekf = filter.mekf();
        check_mode();
        const Matrix21f P0 = mekf.covariance_full();
        Vector12f xlin0;
        xlin0.segment<3>(0) = mekf.get_velocity();
        xlin0.segment<3>(3) = mekf.get_position();
        xlin0.segment<3>(6) = mekf.get_integral_displacement();
        xlin0.segment<3>(9) = mekf.get_world_accel();

        const bool live_pre = fusion_.isLive();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refined_pre = fusion_.hasRefinedMagReference();
        const bool pending_aw_floor = mekf.aw_covariance_floor_pending_;
        const Matrix3f aw_floor_target = mekf.aw_covariance_floor_target_;
        const Matrix3f R_S_pre = mekf.R_S;
        const float rs_scalar_pre = filter.getRSApplied();
        const float pseudo_period_pre = mekf.pseudo_update_period_s_;
        const float tau_b_pre = mekf.tau_bacc_;
        const Matrix3f Q_bacc_pre = mekf.Q_bacc_;
        float pseudo_elapsed_copy = mekf.pseudo_update_elapsed_s_;
        const bool pseudo_due = ocean_imu::kalman::ou_detail::periodic_update_due(
            dt, pseudo_period_pre, pseudo_elapsed_copy);

        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);

        if (!(live_pre && fusion_.isLive())
                || active_pre != mekf.acc_bias_updates_enabled()
                || lock_pre != fusion_.hasMagNorthLock()
                || refined_pre != fusion_.hasRefinedMagReference())
            fail("hybrid/live event inside selected same-mode word");
        check_mode();

        Matrix21f A = Matrix21f::Identity();
        A.block<6,6>(0,0) = mekf.F_AA_scratch_;
        A.block<12,12>(kOffV,kOffV) = mekf.F_LL_scratch_;
        const float tau_b = std::max(1.0e-3f, tau_b_pre);
        const float phi_b = active_pre ? std::exp(-dt / tau_b) : 1.0f;
        A.block<3,3>(kOffBa,kOffBa) = phi_b * Matrix3f::Identity();

        Matrix21f Q = Matrix21f::Zero();
        Q.block<6,6>(0,0) = mekf.Q_AA_scratch_;
        Q.block<12,12>(kOffV,kOffV) = mekf.Q_LL_scratch_;
        if (active_pre) {
            const float qd_scale = -0.5f * tau_b * std::expm1(-2.0f * dt / tau_b);
            Q.block<3,3>(kOffBa,kOffBa) = Q_bacc_pre * qd_scale;
        }

        Matrix21f Ppred = A * P0 * A.transpose() + Q;
        Ppred = 0.5f * (Ppred + Ppred.transpose());
        write_event("prediction", time_s_, P0, Ppred, A,
                    rs_scalar_pre, R_S_pre, 0.0f, 0.0f);
        ++prediction_count_;
        Matrix21f Pcur = Ppred;

        if (pending_aw_floor) {
            Matrix3f Delta = aw_floor_target - Pcur.block<3,3>(kOffAw,kOffAw);
            Delta = 0.5f * (Delta + Delta.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix3f> es(Delta);
            if (es.info() != Eigen::Success) fail("aw-floor eigensolve failed");
            const Vector3f evals = es.eigenvalues().cwiseMax(0.0f);
            const Matrix3f DeltaPlus =
                es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
            Matrix21f Pnext = Pcur;
            Pnext.block<3,3>(kOffAw,kOffAw) += DeltaPlus;
            write_event("aw_floor", time_s_, Pcur, Pnext, Matrix21f::Identity(),
                        rs_scalar_pre, R_S_pre, 0.0f, 0.0f);
            ++floor_count_;
            Pcur = Pnext;
        }

        const Vector12f xlin_pred = mekf.F_LL_scratch_ * xlin0;
        if (pseudo_due) {
            Matrix3x21f Hs = Matrix3x21f::Zero();
            Hs.block<3,3>(0,kOffS) = Matrix3f::Identity();
            Matrix21x3f PCt = Pcur.block<kNX,3>(0,kOffS);
            if (!active_pre) PCt.block<3,3>(kOffBa,0).setZero();
            const Matrix3f S = Pcur.block<3,3>(kOffS,kOffS) + R_S_pre;
            Eigen::LDLT<Matrix3f> ldlt(S);
            if (ldlt.info() != Eigen::Success) fail("S innovation LDLT failed");
            const Matrix21x3f K = ldlt.solve(PCt.transpose()).transpose();
            const Vector3f rS = -xlin_pred.segment<3>(6);
            const Vector3f dtheta = K.topRows<3>() * rS;
            const Matrix21f G = reset_transport(dtheta);
            const Matrix21f C = G * (Matrix21f::Identity() - K * Hs);
            Matrix21f Pnext = joseph_from_pct(Pcur, K, S, PCt);
            Pnext = G * Pnext * G.transpose();
            Pnext = 0.5f * (Pnext + Pnext.transpose());
            write_event("S_zero", time_s_, Pcur, Pnext, C,
                        rs_scalar_pre, R_S_pre, dtheta.norm(), 0.0f);
            ++s_count_;
            Pcur = Pnext;
        }

        const auto& ad = mekf.lastAccDiag();
        if (!ad.accepted) fail("missing accelerometer update inside complete Normal-Live word");
        const Matrix21x3f PCt = mekf.PCt_scratch_;
        const Matrix21x3f K = mekf.K_scratch_;
        Eigen::LDLT<Matrix21f> ldlt(Pcur);
        if (ldlt.info() != Eigen::Success) fail("accelerometer pre-covariance LDLT failed");
        const Matrix21x3f Ht = ldlt.solve(PCt);
        const float lin_resid = (Pcur * Ht - PCt).norm() / std::max(1.0f, PCt.norm());
        const Matrix3x21f H = Ht.transpose();
        const Vector3f dtheta = K.topRows<3>() * ad.r;
        const Matrix21f G = reset_transport(dtheta);
        const Matrix21f C = G * (Matrix21f::Identity() - K * H);
        Matrix21f Precon = joseph_from_pct(Pcur, K, mekf.S_scratch_, PCt);
        Precon = G * Precon * G.transpose();
        Precon = 0.5f * (Precon + Precon.transpose());
        const Matrix21f Ppost = mekf.covariance_full();
        const float cov_resid = std::max(lin_resid, relative_matrix_difference(Precon, Ppost));
        if (!(std::isfinite(cov_resid) && cov_resid <= 2.0e-4f))
            fail("accelerometer covariance reconstruction mismatch");
        write_event("accelerometer", time_s_, Pcur, Ppost, C,
                    rs_scalar_pre, R_S_pre, dtheta.norm(), cov_resid);
        ++acc_count_;

        time_s_ += dt;
        if (time_s_ + 2.0e-6f >= requested_t1_) finish();
    }

    void finish()
    {
        if (finished_ || failed_) return;
        if (std::abs(time_s_ - requested_t1_) > 2.0e-4f)
            fail("ledger endpoint does not match selected complete word");
        const float final_energy = energy(fusion_.raw().mekf().covariance_full(), error_);
        if (!(std::isfinite(final_energy) && final_energy >= 0.0f))
            fail("invalid final information energy");
        const float rho = final_energy / initial_energy_;
        const double telescoping = prediction_delta_ + floor_delta_ + s_delta_ + acc_delta_ + vector_delta_;
        const double total_delta = static_cast<double>(final_energy) - static_cast<double>(initial_energy_);
        const double telescope_error = telescoping - total_delta;
        std::cout << std::setprecision(12)
                  << "OU3_EVENT_LEDGER_DONE mode=" << requested_mode_
                  << " t0=" << actual_t0_ << " t1=" << time_s_
                  << " V0=" << initial_energy_ << " V1=" << final_energy
                  << " rho=" << rho
                  << " prediction_delta=" << prediction_delta_
                  << " floor_delta=" << floor_delta_
                  << " S_delta=" << s_delta_
                  << " accel_delta=" << acc_delta_
                  << " vector_delta=" << vector_delta_
                  << " telescope_error=" << telescope_error
                  << " prediction_count=" << prediction_count_
                  << " floor_count=" << floor_count_
                  << " S_count=" << s_count_
                  << " accel_count=" << acc_count_
                  << " vector_count=" << vector_count_ << "\n";
        finished_ = true;
    }

    [[noreturn]] void fail(const std::string& why)
    {
        if (!failed_) std::cerr << "OU3_EVENT_LEDGER_FAIL: " << why << "\n";
        failed_ = true;
        finished_ = true;
        throw std::runtime_error("operation-observer: " + why);
    }

    mutable Fusion fusion_;
    bool scan_mode_ = false;
    float time_s_ = 0.0f;

    // Scan mode: complete 21x21 map plus shipping covariance boundaries.
    std::ofstream map_trace_;
    std::ofstream cov_trace_;
    int scan_stride_ = 600;
    bool scan_block_started_ = false;
    bool scan_valid_ = false;
    bool scan_hybrid_jump_ = false;
    bool scan_start_live_ = false;
    bool scan_start_active_ = false;
    bool scan_start_lock_ = false;
    bool scan_start_refined_ = false;
    unsigned scan_samples_ = 0;
    int scan_acc_count_ = 0;
    int scan_mag_count_ = 0;
    int scan_pseudo_count_ = 0;
    float scan_start_time_ = 0.0f;
    float scan_start_tau_ = 0.0f;
    float scan_start_sigma_ = 0.0f;
    float scan_start_rs_ = 0.0f;
    float scan_linearization_residual_ = 0.0f;
    Matrix21f scan_accum_ = Matrix21f::Identity();
    Matrix21f scan_cov_start_ = Matrix21f::Zero();

    // Selected-direction event-ledger mode.
    std::ofstream trace_;
    float requested_t0_ = 0.0f;
    float requested_t1_ = 0.0f;
    std::string requested_mode_;
    Vector21f direction_ = Vector21f::Zero();
    Vector21f error_ = Vector21f::Zero();
    float actual_t0_ = 0.0f;
    float initial_energy_ = std::numeric_limits<float>::quiet_NaN();
    bool started_ = false;
    bool finished_ = false;
    bool failed_ = false;
    unsigned sequence_ = 0;
    int prediction_count_ = 0;
    int floor_count_ = 0;
    int s_count_ = 0;
    int acc_count_ = 0;
    int vector_count_ = 0;
    double prediction_delta_ = 0.0;
    double floor_delta_ = 0.0;
    double s_delta_ = 0.0;
    double acc_delta_ = 0.0;
    double vector_delta_ = 0.0;
};

void process_one(const std::string& filename,
                 bool with_mag,
                 const W3dRandomSeeds& seeds)
{
    auto result = process_wave_file_for_tracker<OperationObserverAdapter>(
        filename, kDt, with_mag, true, kMagOdrHz,
        "_fusion_ou3_operation_observer", "_fusion_ou3_operation_observer_nomag",
        seeds, false);
    if (!result) throw std::runtime_error("operation-observer simulation did not produce a result");
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
                      << "Scan env: OU3_LEDGER_MAP_TRACE, OU3_LEDGER_COV_TRACE, "
                         "OU3_LEDGER_MAP_STRIDE (default 600).\n"
                      << "Selected-ledger env: OU3_LEDGER_TRACE, OU3_LEDGER_T0, OU3_LEDGER_T1, "
                         "OU3_LEDGER_MODE=H18|A21, OU3_LEDGER_DIRECTION (21 CSV values).\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown or incomplete argument: " << arg << "\n";
            return 2;
        }
    }
    if (files.size() != 1u) {
        std::cerr << "ERROR: ou3-operation-ledger-sim requires exactly one --input\n";
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
