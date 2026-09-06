#include <algorithm>
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
#include "kalman_ou_common/KalmanOUCoreMath.h"

// Host-only nonlinear complete-word diagnostic. One shipping estimator owns
// the source, covariance, branch decisions, and gains. A covariance-free shadow
// state receives that frozen shipping word and recomputes only its nonlinear
// residuals. Production code and the Riccati recursion are unchanged.
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
constexpr float kTempRefC = 35.0f;
constexpr int kNX = 21;
constexpr int kOffBg = 3;
constexpr int kOffV = 6;
constexpr int kOffP = 9;
constexpr int kOffS = 12;
constexpr int kOffAw = 15;
constexpr int kOffBa = 18;
using Vector21f = Eigen::Matrix<float, kNX, 1>;
using Matrix21f = Eigen::Matrix<float, kNX, kNX>;
using Matrix21x3f = Eigen::Matrix<float, kNX, 3>;
using Matrix3x21f = Eigen::Matrix<float, 3, kNX>;
using Vector12f = Eigen::Matrix<float, 12, 1>;

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
    std::stringstream ss(env_string_required("OU3_SHADOW_DIRECTION"));
    std::string token;
    int i = 0;
    while (std::getline(ss, token, ',')) {
        if (i >= kNX) throw std::runtime_error("OU3_SHADOW_DIRECTION has more than 21 entries");
        char* end = nullptr;
        const float x = std::strtof(token.c_str(), &end);
        if (end == token.c_str() || !std::isfinite(x))
            throw std::runtime_error("OU3_SHADOW_DIRECTION contains invalid value");
        d(i++) = x;
    }
    if (i != kNX) throw std::runtime_error("OU3_SHADOW_DIRECTION must contain exactly 21 entries");
    return d;
}

Quaternionf quaternion_exp(const Vector3f& delta)
{
    const float angle = delta.norm();
    if (!std::isfinite(angle)) throw std::runtime_error("non-finite shadow attitude injection");
    const float half = 0.5f * angle;
    float scale = 0.5f;
    if (angle > 1.0e-7f) scale = std::sin(half) / angle;
    else scale = 0.5f - angle * angle / 48.0f;
    return Quaternionf(std::cos(half), scale * delta.x(), scale * delta.y(), scale * delta.z());
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

struct ShadowState {
    Quaternionf qref = Quaternionf::Identity();
    Vector21f x = Vector21f::Zero();
};

class FrozenGainShadowAdapter final : public IW3dFusionAdapter {
public:
    using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
    using Mekf = std::remove_reference_t<decltype(std::declval<Fusion&>().raw().mekf())>;

    FrozenGainShadowAdapter(bool with_mag,
                            const Vector3f& sigma_a_init,
                            const Vector3f& sigma_g,
                            const Vector3f& sigma_m)
        : requested_t0_(env_float_required("OU3_SHADOW_T0")),
          requested_t1_(env_float_required("OU3_SHADOW_T1")),
          requested_mode_(env_string_required("OU3_SHADOW_MODE")),
          scale_(env_float_required("OU3_SHADOW_SCALE")),
          direction_(parse_direction())
    {
        if (!(requested_t1_ > requested_t0_)) throw std::runtime_error("OU3_SHADOW_T1 must exceed T0");
        if (requested_mode_ != "H18" && requested_mode_ != "A21")
            throw std::runtime_error("OU3_SHADOW_MODE must be H18 or A21");
        if (!std::isfinite(scale_) || scale_ == 0.0f)
            throw std::runtime_error("OU3_SHADOW_SCALE must be finite and nonzero");
        if (requested_mode_ == "H18" && direction_.tail<3>().norm() > 1.0e-12f)
            throw std::runtime_error("H18 shadow direction must have zero b_a tail");

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

        const std::string path = env_string_required("OU3_SHADOW_TRACE");
        trace_.open(path);
        if (!trace_) throw std::runtime_error("cannot open OU3_SHADOW_TRACE");
        trace_ << std::setprecision(12)
               << "time_s,V,theta_norm,state_norm,reconstruction_error,prediction_count,S_count,acc_count,vector_count\n";
    }

    void updateMag(const Vector3f& mag_body_ned) override
    {
        if (finished_) return;
        maybe_start();
        if (!started_) {
            fusion_.updateMag(mag_body_ned);
            return;
        }
        if (!inside_word()) return;

        auto& mekf = fusion_.raw().mekf();
        ShadowState recon = copy_nominal(mekf);
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool lock_pre = fusion_.hasMagNorthLock();
        const bool refined_pre = fusion_.hasRefinedMagReference();

        fusion_.updateMag(mag_body_ned);

        if (active_pre != mekf.acc_bias_updates_enabled()
                || lock_pre != fusion_.hasMagNorthLock()
                || refined_pre != fusion_.hasRefinedMagReference())
            fail("hybrid magnetic/mode event inside selected word");
        check_mode();

        const auto& md = mekf.lastMagDiag();
        if (!md.accepted) return;
        const Matrix21x3f K = mekf.K_scratch_;
        const Vector3f v2ref_used = mekf.v2ref;
        const Vector3f nominal_pred = recon.qref.toRotationMatrix() * v2ref_used;
        const Vector3f measurement_used = md.r + nominal_pred;
        const Vector3f shadow_pred = shadow_.qref.toRotationMatrix() * v2ref_used;
        apply_correction(recon, K, md.r, mekf);
        apply_correction(shadow_, K, measurement_used - shadow_pred, mekf);
        verify_reconstruction(recon, mekf, "vector");
        ++vector_count_;
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
        maybe_start();
        if (!started_) {
            fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
            time_s_ += dt;
            return;
        }
        if (!inside_word()) {
            time_s_ += dt;
            return;
        }

        auto& filter = fusion_.raw();
        filter.apply_pending_online_tune_();
        auto& mekf = filter.mekf();
        check_mode();
        if (mekf.use_imu_lever_arm_)
            fail("frozen shadow is scoped to declared zero/disabled lever arm");

        ShadowState recon = copy_nominal(mekf);
        const Matrix21f P0 = mekf.covariance_full();
        const bool active_pre = mekf.acc_bias_updates_enabled();
        const bool live_pre = fusion_.isLive();
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

        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);

        if (!(live_pre && fusion_.isLive())
                || active_pre != mekf.acc_bias_updates_enabled()
                || lock_pre != fusion_.hasMagNorthLock()
                || refined_pre != fusion_.hasRefinedMagReference())
            fail("hybrid/live event inside selected word");
        check_mode();
        if (filter.accelVibrationGuardEngagement() > 1.0e-7f)
            fail("selected word left declared dormant-transparent vibration branch");

        predict_state(recon, mekf, gyr_meas_ned, dt, active_pre, tau_b_pre);
        predict_state(shadow_, mekf, gyr_meas_ned, dt, active_pre, tau_b_pre);
        ++prediction_count_;

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
        if (pending_aw_floor) {
            Matrix3f Delta = aw_floor_target - Pcur.block<3,3>(kOffAw,kOffAw);
            Delta = 0.5f * (Delta + Delta.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix3f> es(Delta);
            if (es.info() != Eigen::Success) fail("aw-floor eigensolve failed");
            const Vector3f evals = es.eigenvalues().cwiseMax(0.0f);
            Pcur.block<3,3>(kOffAw,kOffAw) +=
                es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
        }

        if (pseudo_due) {
            Matrix21x3f PCt = Pcur.block<kNX,3>(0,kOffS);
            if (!active_pre) PCt.block<3,3>(kOffBa,0).setZero();
            const Matrix3f S = Pcur.block<3,3>(kOffS,kOffS) + R_S_pre;
            Eigen::LDLT<Matrix3f> ldlt(S);
            if (ldlt.info() != Eigen::Success) fail("S innovation LDLT failed");
            const Matrix21x3f K = ldlt.solve(PCt.transpose()).transpose();
            const Vector3f r_nom = -recon.x.segment<3>(kOffS);
            const Vector3f r_shadow = -shadow_.x.segment<3>(kOffS);
            const Vector3f dtheta_nom = K.topRows<3>() * r_nom;
            apply_correction(recon, K, r_nom, mekf);
            apply_correction(shadow_, K, r_shadow, mekf);
            Matrix3x21f Hs = Matrix3x21f::Zero();
            Hs.block<3,3>(0,kOffS) = Matrix3f::Identity();
            const Matrix21f G = reset_transport(dtheta_nom);
            Pcur = joseph_from_pct(Pcur, K, S, PCt);
            Pcur = G * Pcur * G.transpose();
            Pcur = 0.5f * (Pcur + Pcur.transpose());
            ++s_count_;
        }

        const auto& ad = mekf.lastAccDiag();
        if (!ad.accepted) fail("missing accelerometer update inside complete Normal-Live word");
        const Matrix21x3f K = mekf.K_scratch_;
        const Vector3f nominal_pred = accel_prediction(recon, mekf, temperature_c);
        const Vector3f measurement_used = ad.r + nominal_pred;
        const Vector3f shadow_pred = accel_prediction(shadow_, mekf, temperature_c);
        apply_correction(recon, K, ad.r, mekf);
        apply_correction(shadow_, K, measurement_used - shadow_pred, mekf);
        verify_reconstruction(recon, mekf, "accelerometer");
        ++acc_count_;

        time_s_ += dt;
        if (time_s_ + 2.0e-6f >= requested_t1_) finish();
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

    ShadowState copy_nominal(const Mekf& m) const
    {
        ShadowState s;
        s.qref = m.qref;
        s.x = m.xext;
        s.x.head<3>().setZero();
        return s;
    }

    void check_mode()
    {
        const auto& m = fusion_.raw().mekf();
        const std::string mode = m.acc_bias_updates_enabled() ? "A21" : "H18";
        if (!fusion_.isLive() || mode != requested_mode_)
            fail("requested complete-word mode not preserved");
    }

    bool inside_word() const
    {
        return started_ && !finished_ && time_s_ < requested_t1_ - 2.0e-6f;
    }

    void maybe_start()
    {
        if (started_ || finished_) return;
        if (time_s_ + 2.0e-6f < requested_t0_) return;
        if (std::abs(time_s_ - requested_t0_) > 2.0e-4f)
            fail("shadow could not align to requested word start");
        check_mode();
        const auto& mekf = fusion_.raw().mekf();
        if (mekf.use_imu_lever_arm_)
            fail("frozen shadow is scoped to declared zero/disabled lever arm");
        shadow_ = copy_nominal(mekf);
        const Vector21f d = scale_ * direction_;
        shadow_.qref = (quaternion_exp(d.head<3>()) * shadow_.qref).normalized();
        shadow_.x.segment<3>(kOffBg) += d.segment<3>(kOffBg);
        shadow_.x.segment<3>(kOffV) += d.segment<3>(kOffV);
        shadow_.x.segment<3>(kOffP) += d.segment<3>(kOffP);
        shadow_.x.segment<3>(kOffS) += d.segment<3>(kOffS);
        shadow_.x.segment<3>(kOffAw) += d.segment<3>(kOffAw);
        if (requested_mode_ == "A21") shadow_.x.segment<3>(kOffBa) += d.segment<3>(kOffBa);
        started_ = true;
        actual_t0_ = time_s_;
        const Vector21f e0 = pair_error(shadow_, mekf);
        initial_energy_ = energy(mekf.covariance_full(), e0);
        if (!(std::isfinite(initial_energy_) && initial_energy_ > 0.0f))
            fail("invalid initial frozen-shadow energy");
        write_trace();
        std::cout << std::setprecision(12)
                  << "OU3_FROZEN_SHADOW_START mode=" << requested_mode_
                  << " scale=" << scale_ << " t0=" << actual_t0_
                  << " V0=" << initial_energy_ << "\n";
    }

    void predict_state(ShadowState& s,
                       const Mekf& nominal,
                       const Vector3f& gyr_meas_ned,
                       float dt,
                       bool active,
                       float tau_b) const
    {
        const Vector3f gyr = nominal.deheel_vector_(gyr_meas_ned);
        const Vector3f omega = gyr - s.x.segment<3>(kOffBg);
        s.qref = (ocean_imu::kalman::ou_detail::quat_from_delta_theta((-omega * dt).eval()) * s.qref).normalized();

        Vector12f lin;
        lin.segment<3>(0) = s.x.segment<3>(kOffV);
        lin.segment<3>(3) = s.x.segment<3>(kOffP);
        lin.segment<3>(6) = s.x.segment<3>(kOffS);
        lin.segment<3>(9) = s.x.segment<3>(kOffAw);
        lin = nominal.F_LL_scratch_ * lin;
        s.x.segment<3>(kOffV) = lin.segment<3>(0);
        s.x.segment<3>(kOffP) = lin.segment<3>(3);
        s.x.segment<3>(kOffS) = lin.segment<3>(6);
        s.x.segment<3>(kOffAw) = lin.segment<3>(9);
        if (active) s.x.segment<3>(kOffBa) *= std::exp(-dt / std::max(1.0e-3f, tau_b));
    }

    Vector3f accel_prediction(const ShadowState& s, const Mekf& nominal, float tempC) const
    {
        const Vector3f g_world(0.0f, 0.0f, nominal.gravity_magnitude_);
        Vector3f f = s.qref.toRotationMatrix() * (s.x.segment<3>(kOffAw) - g_world);
        f += s.x.segment<3>(kOffBa) + nominal.k_a_ * (tempC - kTempRefC);
        return f;
    }

    void apply_correction(ShadowState& s,
                          const Matrix21x3f& K,
                          const Vector3f& residual,
                          const Mekf& nominal) const
    {
        s.x.noalias() += K * residual;
        const Vector3f dtheta = s.x.head<3>();
        if (!dtheta.allFinite()) throw std::runtime_error("non-finite frozen-shadow correction");
        s.qref = (ocean_imu::kalman::ou_detail::quat_from_delta_theta(dtheta) * s.qref).normalized();
        s.x.head<3>().setZero();
        const float limit = nominal.accel_bias_limit();
        if (limit > 0.0f) {
            auto ba = s.x.segment<3>(kOffBa);
            if (!ba.allFinite()) throw std::runtime_error("non-finite frozen-shadow b_a");
            const float n = ba.norm();
            if (n > limit) ba *= limit / n;
        }
    }

    Vector21f pair_error(const ShadowState& s, const Mekf& nominal) const
    {
        Vector21f e = Vector21f::Zero();
        e.head<3>() = quaternion_log(s.qref * nominal.qref.conjugate());
        e.segment<3>(kOffBg) = s.x.segment<3>(kOffBg) - nominal.xext.segment<3>(kOffBg);
        e.segment<3>(kOffV) = s.x.segment<3>(kOffV) - nominal.xext.segment<3>(kOffV);
        e.segment<3>(kOffP) = s.x.segment<3>(kOffP) - nominal.xext.segment<3>(kOffP);
        e.segment<3>(kOffS) = s.x.segment<3>(kOffS) - nominal.xext.segment<3>(kOffS);
        e.segment<3>(kOffAw) = s.x.segment<3>(kOffAw) - nominal.xext.segment<3>(kOffAw);
        e.segment<3>(kOffBa) = s.x.segment<3>(kOffBa) - nominal.xext.segment<3>(kOffBa);
        if (dimension() == 18) e.tail<3>().setZero();
        return e;
    }

    float energy(const Matrix21f& Pfull, const Vector21f& e) const
    {
        if (dimension() == 18) {
            Eigen::Matrix<float,18,18> P = Pfull.topLeftCorner<18,18>();
            P = 0.5f * (P + P.transpose());
            Eigen::LDLT<Eigen::Matrix<float,18,18>> ldlt(P);
            if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
            return e.head<18>().dot(ldlt.solve(e.head<18>()));
        }
        Matrix21f P = 0.5f * (Pfull + Pfull.transpose());
        Eigen::LDLT<Matrix21f> ldlt(P);
        if (ldlt.info() != Eigen::Success) return std::numeric_limits<float>::quiet_NaN();
        return e.dot(ldlt.solve(e));
    }

    void verify_reconstruction(const ShadowState& recon, const Mekf& nominal, const char* event)
    {
        const Vector21f e = pair_error(recon, nominal);
        const float err = e.head(dimension()).norm();
        max_reconstruction_error_ = std::max(max_reconstruction_error_, err);
        if (!(std::isfinite(err) && err <= 5.0e-5f))
            fail(std::string("nominal state reconstruction mismatch after ") + event);
    }

    void write_trace()
    {
        if (!started_) return;
        const auto& mekf = fusion_.raw().mekf();
        const Vector21f e = pair_error(shadow_, mekf);
        const float V = energy(mekf.covariance_full(), e);
        trace_ << time_s_ << ',' << V << ',' << e.head<3>().norm() << ','
               << e.head(dimension()).norm() << ',' << max_reconstruction_error_ << ','
               << prediction_count_ << ',' << s_count_ << ',' << acc_count_ << ',' << vector_count_ << '\n';
    }

    void finish()
    {
        if (finished_ || failed_) return;
        if (std::abs(time_s_ - requested_t1_) > 2.0e-4f)
            fail("frozen shadow endpoint does not match selected word");
        const auto& mekf = fusion_.raw().mekf();
        const Vector21f ef = pair_error(shadow_, mekf);
        const float final_energy = energy(mekf.covariance_full(), ef);
        if (!(std::isfinite(final_energy) && final_energy >= 0.0f))
            fail("invalid final frozen-shadow energy");
        const float rho = final_energy / initial_energy_;
        write_trace();
        std::cout << std::setprecision(12)
                  << "OU3_FROZEN_SHADOW_DONE mode=" << requested_mode_
                  << " scale=" << scale_ << " t0=" << actual_t0_ << " t1=" << time_s_
                  << " V0=" << initial_energy_ << " V1=" << final_energy << " rho=" << rho
                  << " reconstruction_max=" << max_reconstruction_error_
                  << " prediction_count=" << prediction_count_
                  << " S_count=" << s_count_ << " accel_count=" << acc_count_
                  << " vector_count=" << vector_count_ << "\n";
        finished_ = true;
    }

    [[noreturn]] void fail(const std::string& why)
    {
        if (!failed_) std::cerr << "OU3_FROZEN_SHADOW_FAIL: " << why << "\n";
        failed_ = true;
        finished_ = true;
        throw std::runtime_error("frozen-gain-shadow: " + why);
    }

    mutable Fusion fusion_;
    std::ofstream trace_;
    float requested_t0_ = 0.0f;
    float requested_t1_ = 0.0f;
    std::string requested_mode_;
    float scale_ = 0.0f;
    Vector21f direction_ = Vector21f::Zero();
    ShadowState shadow_;
    float time_s_ = 0.0f;
    float actual_t0_ = 0.0f;
    float initial_energy_ = std::numeric_limits<float>::quiet_NaN();
    float max_reconstruction_error_ = 0.0f;
    bool started_ = false;
    bool finished_ = false;
    bool failed_ = false;
    int prediction_count_ = 0;
    int s_count_ = 0;
    int acc_count_ = 0;
    int vector_count_ = 0;
};

void process_one(const std::string& filename,
                 bool with_mag,
                 const W3dRandomSeeds& seeds)
{
    auto result = process_wave_file_for_tracker<FrozenGainShadowAdapter>(
        filename, kDt, with_mag, true, kMagOdrHz,
        "_fusion_ou3_frozen_shadow", "_fusion_ou3_frozen_shadow_nomag", seeds, false);
    if (!result) throw std::runtime_error("frozen-shadow simulation did not produce a result");
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
                      << "Env: OU3_SHADOW_TRACE, OU3_SHADOW_T0, OU3_SHADOW_T1, "
                         "OU3_SHADOW_MODE=H18|A21, OU3_SHADOW_DIRECTION, OU3_SHADOW_SCALE.\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown or incomplete argument: " << arg << "\n";
            return 2;
        }
    }
    if (files.size() != 1u) {
        std::cerr << "ERROR: ou3-frozen-gain-shadow-sim requires exactly one --input\n";
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
