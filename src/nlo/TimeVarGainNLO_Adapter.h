#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Reusable adapter for TimeVaryingGainNLO.

  Clients choose compile-time options:
      TimeVarGainNloAdapter<false, NloMagType::None, float>
      TimeVarGainNloAdapter<false, NloMagType::Magnetometer, float>
      TimeVarGainNloAdapter<true,  NloMagType::Compass, float>
      TimeVarGainNloAdapter<true,  NloMagType::Magnetometer, float>

  Frames:
    Inputs are BODY-frame NED/marine axes:
      x forward, y starboard, z down

    specific_force_b_mps2:
      level/still approximately [0, 0, -g]

    gyro_b_rad_s:
      body angular rate, rad/s

    State frame:
      NED, +Z down

    Convenience outputs:
      disp_zu / vel_zu / acc_zu use Z-up sign for charting/heave output.

  Startup behavior:
    - optional Mahony-like bootstrap is ON by default
    - it initializes q from accel/mag/compass if available
    - it then gyro-propagates during startup and uses gated accel correction
    - after bootstrap, it seeds TimeVaryingGainNLO and stops

  Important:
    Wave-time accelerometer tilt trim is OFF by default because it can chase
    horizontal wave acceleration and worsen roll/pitch.
*/

#if defined(ARDUINO)
  #include <Arduino.h>
  #if __has_include(<ArduinoEigenDense.h>)
    #include <ArduinoEigenDense.h>
  #else
    #include <Eigen/Dense>
    #include <Eigen/Geometry>
  #endif
#else
  #include <Eigen/Dense>
  #include <Eigen/Geometry>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "nlo/TimeVaryingGainNLO.h"
#include "freq/FrequencyTrackerPolicy.h"

template <bool WithGNSS = false,
          NloMagType Mag = NloMagType::None,
          typename R = float,
          TrackerType FreqTT = TrackerType::PLL>
class TimeVarGainNloAdapter {
public:
    using Filter = TimeVaryingGainNLO<WithGNSS, Mag, R>;
    using WaveFreqTracker = TrackerPolicy<FreqTT>;
    using WaveFreqTrackerConfig =
        marine_obs::detail::tracker_config_t<WaveFreqTracker, R>;

    using Vec2 = Eigen::Matrix<R, 2, 1>;
    using Vec3 = Eigen::Matrix<R, 3, 1>;
    using Mat3 = Eigen::Matrix<R, 3, 3>;
    using Quat = Eigen::Quaternion<R>;
    using Aux  = typename Filter::Aux;

    static constexpr bool kWithGNSS = WithGNSS;
    static constexpr NloMagType kMagType = Mag;

    struct Telemetry {
        R k1 = std::numeric_limits<R>::quiet_NaN();
        R k2 = std::numeric_limits<R>::quiet_NaN();
        R kI = std::numeric_limits<R>::quiet_NaN();
        R vartheta = std::numeric_limits<R>::quiet_NaN();
        R p0z_hat = std::numeric_limits<R>::quiet_NaN();

        Vec3 xi_n = Vec3::Constant(std::numeric_limits<R>::quiet_NaN());
        Vec3 fhat_n = Vec3::Constant(std::numeric_limits<R>::quiet_NaN());
        Vec3 sigma_b = Vec3::Constant(std::numeric_limits<R>::quiet_NaN());
        Vec3 gyro_bias_b = Vec3::Constant(std::numeric_limits<R>::quiet_NaN());

        R xi_norm = std::numeric_limits<R>::quiet_NaN();
        R fhat_norm = std::numeric_limits<R>::quiet_NaN();
        R sigma_norm = std::numeric_limits<R>::quiet_NaN();
        R gyro_bias_norm = std::numeric_limits<R>::quiet_NaN();

        R tilt_err_norm = std::numeric_limits<R>::quiet_NaN();
        R tilt_correction_norm = std::numeric_limits<R>::quiet_NaN();
        R tilt_down_lpf_valid = R(0);

        R mahony_boot_time_s = R(0);
        R mahony_boot_good_time_s = R(0);
        R mahony_boot_bias_norm = R(0);

        R theta = std::numeric_limits<R>::quiet_NaN();
        R wave_freq_hz = std::numeric_limits<R>::quiet_NaN();
        R wave_freq_confidence = std::numeric_limits<R>::quiet_NaN();
        bool wave_freq_locked = false;
    };

    struct Snapshot {
        bool initialized = false;
        bool warming = true;
        bool tilt_trim_active = false;
        bool mahony_bootstrap_active = false;

        R time_s = R(0);
        R time_since_init_s = R(0);
        R init_good_time_s = R(0);
        R init_total_time_s = R(0);

        Quat q_nb = Quat::Identity();

        Vec3 position_ned = Vec3::Zero();
        Vec3 velocity_ned = Vec3::Zero();

        Vec3 specific_force_ned = Vec3::Zero();
        Vec3 inertial_accel_ned = Vec3::Zero();

        Vec3 disp_zu = Vec3::Zero();
        Vec3 vel_zu = Vec3::Zero();
        Vec3 acc_zu = Vec3::Zero();

        Vec3 euler_rad = Vec3::Zero();

        Telemetry tvg;
    };

    struct Config {
        R gravity_mps2 = R(9.80665);

        typename Filter::Config filter{};

        /*
          Basic averaging fallback startup.
          Used only when mahony_bootstrap_enabled == false.
        */
        R init_required_good_time_s = R(2.0);
        R init_max_wait_s = R(4.0);
        R init_gyro_max_rad_s = R(0.08);
        R init_acc_norm_tol_frac = R(0.18);
        R yaw_seed_rad = R(0);

        /*
          Mahony-like bootstrap stage.
          Enabled by default. This is a startup stage only; it stops once the
          NLO is initialized.
        */
        bool mahony_bootstrap_enabled = true;
        R mahony_bootstrap_min_time_s = R(2.0);
        R mahony_bootstrap_max_time_s = R(6.0);

        R mahony_twoKp = R(0.35);
        R mahony_twoKi = R(0.002);

        R mahony_acc_norm_tol_frac = R(0.18);
        R mahony_gyro_max_rad_s = R(0.20);
        R mahony_down_lpf_tau_s = R(0.50);
        R mahony_bias_limit_rad_s = R(0.03);

        /*
          Optional roll/pitch-only tilt trim after NLO init.

          Default is OFF because in waves, accelerometer tilt trim can chase
          horizontal wave acceleration and make roll/pitch worse.

          Enable only for bench/calm/still startup or with a real stillness gate.

          tilt_trim_duration_s:
            > 0 : stop after this many seconds after initialization
            <=0 : run whenever enabled
        */
        bool tilt_trim_enabled = false;
        R tilt_trim_duration_s = R(0);
        R tilt_trim_tau_s = R(12.0);
        R tilt_trim_acc_lpf_tau_s = R(0.75);
        R tilt_trim_max_rate_rad_s = R(0.025);
        R tilt_trim_gyro_max_rad_s = R(0.08);
        R tilt_trim_acc_norm_tol_frac = R(0.10);

        /*
          Optional adapter-level roll/pitch gyro-bias learner for tilt trim.
          Default OFF for wave-time operation.
        */
        bool tilt_bias_enabled = false;
        R tilt_bias_ki = R(0);
        R tilt_bias_limit_rad_s = R(0);

        /*
          If false, update() returns before initialized and outputs remain zero.
          For real device startup, false is safer.
        */
        bool run_filter_before_initialized = false;

        /*
          Sea-state schedule for theta.

          theta is the paper's own scalar tuning parameter of the
          translational observer. It enters only through
          K(t) = vartheta(t)*theta*L_theta^-1*K0 (eq. (33)), which is an exact
          frequency scaling of the whole aiding loop by theta. Every published
          gain, and every ratio between them, is untouched by it.

          With the paper's K0 the aiding loop gives the heave path an
          effective high-pass corner of

              omega_c = theta * K_pz_p0z / K_p0z_p0z = 0.4125 * theta  rad/s

          A wave sensor needs that corner below the wave band: too high and
          the aiding loop attenuates and phase-shifts the wave itself, too low
          and double-integrated accelerometer drift is no longer suppressed.
          The paper's theta = 1 puts the corner at 0.41 rad/s, a 15 s period,
          which is inside the ocean wave band. That single fact accounts for
          essentially all of the observer's heave error at theta = 1.

          So theta is scheduled to hold the corner at a fixed fraction of the
          dominant wave frequency:

              theta = theta_from_omega_gain * omega_peak

          which places omega_c at 0.4125 * theta_from_omega_gain * omega_peak.
          The default 0.56 puts it a factor 4.3 below the peak, which is the
          flat optimum measured by a fixed-theta sweep across Hs = 0.27 m to
          8.5 m: the per-record optima came out at theta/omega_peak = 0.51,
          0.52, 0.53, 0.54, 0.58, 0.59, 0.61, 0.61.

          Scheduling a translational gain scalar on a measured condition
          signal is the paper's own construction; it schedules vartheta(t) on
          reported GNSS quality and says other laws may be used in its place.
          Here there is no receiver to report quality, and vartheta cannot
          move the corner anyway because it multiplies all four gains equally.
          theta is the parameter that can.

          Note this departs from the paper's theta >= 1. That bound comes from
          Theorem 1, where a large theta is what dominates the nonlinear
          coupling term; it is sufficient, not necessary, and the linear
          aiding loop is Hurwitz for any theta > 0 at these gains.

          Set auto_theta_from_wave_freq = false to run a fixed
          Config::filter.theta, e.g. the paper's 1.0.
        */
        bool auto_theta_from_wave_freq = true;

        R theta_from_omega_gain = R(0.56);
        R theta_min = R(0.15);
        R theta_max = R(1.60);

        // Time constant of the lag applied to theta itself, so the aiding
        // loop is never yanked when the tracker steps.
        R theta_smooth_tau_s = R(30.0);

        // Below this tracker confidence the schedule holds its last value.
        R theta_min_confidence = R(0.25);

        // theta is recomputed on this cadence, not every sample.
        R theta_update_period_s = R(0.5);

        // Corner of the fixed reference-heave channel the wave tracker runs
        // on. Below the longest swell of interest, above the drift band.
        R reference_heave_hp_hz = R(0.030);

        /*
          DC removal on the reported position and velocity.

          The absolute level of position is not observable here. The only
          thing aiding it is the virtual measurement p0n = 0, whose innovation
          is high-passed per eq. (27)-(28), and a high-pass has exactly zero
          gain at DC for any Th. So the observer can pin the wave-band motion
          and cannot pin the constant it sits on.

          Left alone that constant is not small. Measured against the
          simulator, reported heave carries a standing offset of -1.2 m on an
          Hs = 8.5 m sea. It comes from a rectified bias in (9d): the term
          -R*S(sigma)*f_b has a nonzero mean over the wave cycle, which biases
          xi and therefore the estimated vertical acceleration by about
          -0.012 m/s^2 at that sea state, and the aiding loop maps that into
          velocity and position through its own time constants. The bias
          scales with sea state, and the position offset it produces scales
          with 1/theta^2, so it grew when theta was moved below the wave band.

          Removing the acceleration bias at its source would be the better
          fix and is not done here. What is done is to report position and
          velocity through a first-order high-pass, so the output is the
          wave-band motion the observer can actually claim, with the
          unobservable constant taken out. This is what the zero-mean aiding
          assumption already asserts, so it is consistent with the observer
          rather than a cosmetic trim.

          The corner sits below the wave band: at tau = 50 s a 20 s swell sees
          0.2 % amplitude loss and 3.6 degrees of phase, and an 11 s sea sees
          0.3 % and 4.1 degrees.

          50 s is a measured optimum, not a round number. Sweeping it over the
          scored records gives mean raw Z RMS of 7.18 % at 20 s, 6.88 % at
          35 s, 6.83 % at 50 s, 6.97 % at 70 s, 7.51 % at 100 s and 9.97 % at
          400 s. Too short costs wave-band fidelity; too long leaves the
          observer's own low-frequency wander in the output.

          Reported acceleration is left alone. Its bias is the cause here, not
          a symptom of unobservable DC, and high-passing it would hide the one
          channel that shows the defect.

          Axis selection follows observability, not convenience: the vertical
          is always high-passed because the virtual measurement always aids it
          through the innovation high-pass, while the horizontal axes are
          high-passed only when no position reference is present. With GNSS
          the horizontal DC is observable and is reported as estimated.

          filter().positionNED() and filter().velocityNED() are unaffected and
          still return the raw observer states.
        */
        bool report_highpass_enabled = true;
        R report_highpass_tau_s = R(50.0);

        /*
          Same tracker the adaptive PII observer uses, with three settings
          changed because this one is driven by a displacement-shaped
          reference channel rather than by acceleration:

            f_max_hz     0.35 -> 0.50, so short wind seas near a 2.5 s period
                                 are still inside the search band
            pre_lp_hz    0.45 -> 0.60, to match the widened band
            lock_rms_min 0.012 -> 5e-4, because the policy normalizes its
                                 input by g and a displacement of a few
                                 centimetres would otherwise never lock
        */
        WaveFreqTrackerConfig wave_freq_tracker = [] {
            auto c = marine_obs::detail::
                make_default_tracker_config<WaveFreqTrackerConfig, R>();
            if constexpr (requires { c.f_max_hz; })     c.f_max_hz = R(0.50);
            if constexpr (requires { c.pre_lp_hz; })    c.pre_lp_hz = R(0.60);
            if constexpr (requires { c.lock_rms_min; }) c.lock_rms_min = R(5e-4);
            return c;
        }();
    };

    explicit TimeVarGainNloAdapter(const Config& cfg = makeDefaultConfig())
        : cfg_(cfg),
          filter_(cfg.filter)
    {
        reset();
    }

    static Config makeDefaultConfig() {
        Config cfg{};

        cfg.gravity_mps2 = R(9.80665);

        auto& f = cfg.filter;
        f.gravity_mps2 = cfg.gravity_mps2;

        f.gyro_bias_limit_rad_s = R(0.10);
        f.max_specific_force_mps2 = R(30.0);

        /*
          All gains below are the values published in Bryne/Fossen/Johansen,
          Sec. IV-C "Tuning and Gain Structure". Do not retune them here; use
          Config::filter.theta, which is the paper's own scalar tuning
          parameter, to move the translational loop in frequency.
        */

        // ga_dot = -(1/T) ga + (1/T) ka, T = 25 s,
        // ka = [20, 20, 1] for t <= 100 s, ka = [0.55, 1, 0.01] afterwards.
        f.use_time_varying_attitude_gains = true;
        f.attitude_gain_tau_s = R(25.0);
        f.attitude_gain_switch_s = R(100.0);

        f.k1_initial = R(20.0);
        f.k2_initial = R(20.0);

        f.k1_nominal = R(0.55);
        f.k2_nominal = R(1.0);
        f.kI_nominal = R(0.01);

        /*
          The paper ramps the integral gain from kI = 1 as well, which is safe
          there because Sigma2 is pinned by 5 Hz GNSS position from the start,
          so the injection term it integrates is meaningful immediately.

          Without GNSS, Sigma2's transient lasts far longer, and kI = 1 rad/s
          integrates that transient straight into the gyro bias estimate. The
          bias reaches gyro_bias_limit_rad_s, which rotates the attitude,
          which grows the injection term further, and the interconnection
          diverges: roll and pitch RMS run to tens of degrees. The paper's
          parameter projection does not stop this because the admissible ball
          used here is two orders of magnitude wider than the true biases.

          The proportional part of the ramp, k1 and k2 = 20, is harmless and
          is kept: it converges attitude quickly and is measurably fine.
        */
        f.kI_initial = WithGNSS ? R(1.0) : f.kI_nominal;

        // K0 = P*C', P from the ARE (12) with tau = 1/2 and
        // Q = blkdiag{50, 0.5*I3, 0.08*I3, 0.0025*I3}.
        // The four-state column is used by every axis aided by the virtual
        // zero-mean measurement; the three-state column by GNSS-aided axes.
        f.K_p0z_p0z = R(5.4295);
        f.K_pz_p0z  = R(2.2396);
        f.K_vz_p0z  = R(0.4454);
        f.K_xiz_p0z = R(0.0354);

        f.K_pp_scalar  = R(0.9513);
        f.K_vp_scalar  = R(0.3275);
        f.K_xip_scalar = R(0.0354);

        f.theta = R(1.0);

        // vartheta(t) = vartheta0 + vartheta1 + vartheta2, eq. (29)-(31),
        // with vartheta0 = 0.5, vartheta1 = b*exp(-a*ef), a = 2, b = 1.5,
        // ef a 125 s low-pass of the reported horizontal GNSS RMS error, and
        // vartheta2 a 25 s lag toward 1 for t <= 100 s and toward 0 after.
        f.use_time_varying_tmo_gain = true;
        f.vartheta0 = R(0.5);
        f.vartheta1_a = R(2.0);
        f.vartheta1_b = R(1.5);
        f.gnss_rms_lpf_tau_s = R(125.0);

        // No receiver is present to report a horizontal RMS error, so the
        // vartheta1 quality bonus is not claimed and vartheta settles on the
        // paper's floor vartheta0.
        f.vartheta1_without_gnss = R(0.0);

        f.vartheta2_tau_s = R(25.0);
        f.vartheta2_switch_s = R(100.0);

        f.p0z_highpass_tau_s = R(600.0);
        f.use_triad_style_force_injection = true;

        f.use_virtual_horizontal_position = !WithGNSS;

        return cfg;
    }

    void reset() {
        filter_ = Filter(cfg_.filter);
        filter_.reset();

        marine_obs::detail::tracker_configure(wave_freq_tracker_,
                                              cfg_.wave_freq_tracker);
        marine_obs::detail::tracker_reset(
            wave_freq_tracker_,
            marine_obs::detail::tracker_init_frequency_hz(cfg_.wave_freq_tracker,
                                                          R(0.12))
        );

        theta_ = cfg_.filter.theta;
        theta_sched_accum_s_ = R(0);
        theta_acquired_ = false;
        report_lpf_p_.setZero();
        report_lpf_v_.setZero();
        report_hp_primed_ = false;
        ref_heave_v_ = R(0);
        ref_heave_z_ = R(0);
        last_wave_freq_hz_ = R(0);
        last_wave_freq_conf_ = R(0);
        last_wave_freq_locked_ = false;

        aux_ = Aux{};
        resetAuxValidity_();

        initialized_ = false;
        time_s_ = R(0);
        time_since_init_s_ = R(0);

        init_total_time_s_ = R(0);
        init_good_time_s_ = R(0);

        init_good_count_ = 0;
        init_all_count_ = 0;
        init_good_mag_count_ = 0;
        init_all_mag_count_ = 0;

        init_good_acc_sum_.setZero();
        init_all_acc_sum_.setZero();
        init_good_mag_sum_.setZero();
        init_all_mag_sum_.setZero();

        down_lpf_b_.setZero();
        have_down_lpf_ = false;

        tilt_bias_b_.setZero();
        last_tilt_err_norm_ = R(0);
        last_tilt_correction_norm_ = R(0);

        boot_q_nb_.setIdentity();
        boot_bias_b_.setZero();
        boot_down_lpf_b_.setZero();

        boot_started_ = false;
        boot_have_down_lpf_ = false;
        boot_time_s_ = R(0);
        boot_good_time_s_ = R(0);
    }

    bool initialized() const { return initialized_; }

    bool warming() const {
        return !initialized_ ||
               (cfg_.tilt_trim_enabled &&
                cfg_.tilt_trim_duration_s > R(0) &&
                time_since_init_s_ < cfg_.tilt_trim_duration_s);
    }

    Filter& filter() { return filter_; }
    const Filter& filter() const { return filter_; }

    Config& config() { return cfg_; }
    const Config& config() const { return cfg_; }

    Aux& aux() { return aux_; }
    const Aux& aux() const { return aux_; }

    void setCompassHeadingRad(R heading_rad, bool valid = true) {
        if constexpr (Mag == NloMagType::Compass) {
            aux_.compass_heading_rad = heading_rad;
            aux_.compass_valid = valid;
        } else {
            (void)heading_rad;
            (void)valid;
        }
    }

    void clearCompass() {
        if constexpr (Mag == NloMagType::Compass) {
            aux_.compass_valid = false;
        }
    }

    void setMagBody(const Vec3& mag_b, bool valid = true) {
        if constexpr (Mag == NloMagType::Magnetometer) {
            aux_.mag_b = mag_b;
            aux_.mag_valid = valid;
        } else {
            (void)mag_b;
            (void)valid;
        }
    }

    void clearMag() {
        if constexpr (Mag == NloMagType::Magnetometer) {
            aux_.mag_valid = false;
        }
    }

    void setGnssXY(const Vec2& pos_n_m, R rms_xy_m, bool valid = true) {
        if constexpr (WithGNSS) {
            aux_.gnss.position_n_m = pos_n_m;
            aux_.gnss.rms_xy_m = rms_xy_m;
            aux_.gnss.valid = valid;
        } else {
            (void)pos_n_m;
            (void)rms_xy_m;
            (void)valid;
        }
    }

    void clearGnss() {
        if constexpr (WithGNSS) {
            aux_.gnss.valid = false;
        }
    }

    void setMagneticDeclinationRad(R declination_rad) {
        filter_.setMagneticDeclinationRad(declination_rad);
        cfg_.filter.magnetic_declination_rad = declination_rad;
    }

    void setAccelBiasBody(const Vec3& bias_mps2) {
        filter_.setAccelBiasBody(bias_mps2);
    }

    void setGyroBiasBody(const Vec3& bias_rad_s) {
        filter_.setGyroBiasBody(bias_rad_s);
        tilt_bias_b_ = bias_rad_s;
        boot_bias_b_ = bias_rad_s;

        if constexpr (Mag == NloMagType::None) {
            tilt_bias_b_.z() = R(0);
            boot_bias_b_.z() = R(0);
        }
    }

    void setQuaternionBodyToNED(const Quat& q_nb) {
        filter_.setQuaternionBodyToNED(q_nb);
        initialized_ = true;
        time_since_init_s_ = R(0);
        have_down_lpf_ = false;
        boot_started_ = false;
    }

    void setPositionNED(const Vec3& p_n_m) {
        filter_.setPositionNED(p_n_m);
    }

    void setVelocityNED(const Vec3& v_n_mps) {
        filter_.setVelocityNED(v_n_mps);
    }

    bool update(R dt,
                const Vec3& gyro_b_rad_s,
                const Vec3& specific_force_b_mps2)
    {
        return update(dt, gyro_b_rad_s, specific_force_b_mps2, aux_);
    }

    bool update(R dt,
                const Vec3& gyro_b_rad_s,
                const Vec3& specific_force_b_mps2,
                const Aux& aux)
    {
        if (!(dt > R(0)) || !isFinite_(dt)) {
            return false;
        }

        time_s_ += dt;

        if (!initialized_) {
            updateInitialization_(dt, gyro_b_rad_s, specific_force_b_mps2, aux);

            if (!initialized_ && !cfg_.run_filter_before_initialized) {
                return false;
            }
        } else {
            time_since_init_s_ += dt;
        }

        if (cfg_.tilt_bias_enabled) {
            Vec3 current = filter_.gyroBiasBody();
            current.x() = tilt_bias_b_.x();
            current.y() = tilt_bias_b_.y();

            if constexpr (Mag == NloMagType::None) {
                current.z() = R(0);
            } else {
                current.z() = tilt_bias_b_.z();
            }

            filter_.setGyroBiasBody(current);
        }

        updateThetaSchedule_(dt);

        filter_.update(dt, gyro_b_rad_s, specific_force_b_mps2, aux);

        updateReportHighpass_(dt, hasUsableGNSS_(aux));

        if (initialized_) {
            applyOptionalTiltTrim_(dt, gyro_b_rad_s, specific_force_b_mps2);
        }

        return initialized_;
    }

    Snapshot snapshot() const {
        Snapshot s;

        s.initialized = initialized_;
        s.warming = warming();
        s.tilt_trim_active = tiltTrimCurrentlyActive_();
        s.mahony_bootstrap_active =
            (!initialized_ && cfg_.mahony_bootstrap_enabled && boot_started_);

        s.time_s = time_s_;
        s.time_since_init_s = time_since_init_s_;
        s.init_good_time_s = init_good_time_s_;
        s.init_total_time_s = init_total_time_s_;

        s.q_nb = initialized_ ? filter_.quaternionBodyToNED() : boot_q_nb_;

        // Reported with the unobservable DC removed.
        // See Config::report_highpass_enabled.
        s.position_ned = filter_.positionNED() - report_lpf_p_;
        s.velocity_ned = filter_.velocityNED() - report_lpf_v_;
        s.specific_force_ned = filter_.specificForceNED();

        Vec3 gravity_n;
        gravity_n << R(0), R(0), cfg_.gravity_mps2;

        s.inertial_accel_ned = s.specific_force_ned + gravity_n;

        s.disp_zu << s.position_ned.x(), s.position_ned.y(), -s.position_ned.z();
        s.vel_zu  << s.velocity_ned.x(), s.velocity_ned.y(), -s.velocity_ned.z();
        s.acc_zu  << s.inertial_accel_ned.x(), s.inertial_accel_ned.y(), -s.inertial_accel_ned.z();

        s.euler_rad = initialized_ ? filter_.eulerRad() : eulerFromQuat_(boot_q_nb_);

        s.tvg.k1 = filter_.gainK1();
        s.tvg.k2 = filter_.gainK2();
        s.tvg.kI = filter_.gainKI();
        s.tvg.vartheta = filter_.gainVartheta();
        s.tvg.p0z_hat = filter_.integratedVerticalPositionState();

        s.tvg.xi_n = filter_.xiNED();
        s.tvg.fhat_n = filter_.specificForceNED();
        s.tvg.sigma_b = filter_.sigmaBody();
        s.tvg.gyro_bias_b = filter_.gyroBiasBody();

        s.tvg.xi_norm = s.tvg.xi_n.norm();
        s.tvg.fhat_norm = s.tvg.fhat_n.norm();
        s.tvg.sigma_norm = s.tvg.sigma_b.norm();
        s.tvg.gyro_bias_norm = s.tvg.gyro_bias_b.norm();

        s.tvg.tilt_err_norm = last_tilt_err_norm_;
        s.tvg.tilt_correction_norm = last_tilt_correction_norm_;
        s.tvg.tilt_down_lpf_valid = have_down_lpf_ ? R(1) : R(0);

        s.tvg.mahony_boot_time_s = boot_time_s_;
        s.tvg.mahony_boot_good_time_s = boot_good_time_s_;
        s.tvg.mahony_boot_bias_norm = boot_bias_b_.norm();

        s.tvg.theta = filter_.config().theta;
        s.tvg.wave_freq_hz = last_wave_freq_hz_;
        s.tvg.wave_freq_confidence = last_wave_freq_conf_;
        s.tvg.wave_freq_locked = last_wave_freq_locked_;

        return s;
    }

    R heaveUpM() const {
        return -filter_.positionNED().z();
    }

    R verticalVelocityUpMps() const {
        return -filter_.velocityNED().z();
    }

    R verticalAccelerationUpMps2() const {
        Vec3 gravity_n;
        gravity_n << R(0), R(0), cfg_.gravity_mps2;
        return -(filter_.specificForceNED() + gravity_n).z();
    }

private:
    Config cfg_;
    Filter filter_;
    Aux aux_{};

    WaveFreqTracker wave_freq_tracker_{};
    R theta_ = R(1);
    R theta_sched_accum_s_ = R(0);
    bool theta_acquired_ = false;
    Vec3 report_lpf_p_ = Vec3::Zero();
    Vec3 report_lpf_v_ = Vec3::Zero();
    bool report_hp_primed_ = false;
    R ref_heave_v_ = R(0);
    R ref_heave_z_ = R(0);
    R last_wave_freq_hz_ = R(0);
    R last_wave_freq_conf_ = R(0);
    bool last_wave_freq_locked_ = false;

    bool initialized_ = false;

    R time_s_ = R(0);
    R time_since_init_s_ = R(0);

    R init_total_time_s_ = R(0);
    R init_good_time_s_ = R(0);

    int init_good_count_ = 0;
    int init_all_count_ = 0;
    int init_good_mag_count_ = 0;
    int init_all_mag_count_ = 0;

    Vec3 init_good_acc_sum_ = Vec3::Zero();
    Vec3 init_all_acc_sum_ = Vec3::Zero();

    Vec3 init_good_mag_sum_ = Vec3::Zero();
    Vec3 init_all_mag_sum_ = Vec3::Zero();

    Vec3 down_lpf_b_ = Vec3::Zero();
    bool have_down_lpf_ = false;

    Vec3 tilt_bias_b_ = Vec3::Zero();
    R last_tilt_err_norm_ = R(0);
    R last_tilt_correction_norm_ = R(0);

    Quat boot_q_nb_ = Quat::Identity();
    Vec3 boot_bias_b_ = Vec3::Zero();
    Vec3 boot_down_lpf_b_ = Vec3::Zero();

    bool boot_started_ = false;
    bool boot_have_down_lpf_ = false;
    R boot_time_s_ = R(0);
    R boot_good_time_s_ = R(0);

    void resetAuxValidity_() {
        if constexpr (Mag == NloMagType::Compass) {
            aux_.compass_valid = false;
        }

        if constexpr (Mag == NloMagType::Magnetometer) {
            aux_.mag_valid = false;
        }

        if constexpr (WithGNSS) {
            aux_.gnss.valid = false;
        }
    }

    /*
      Hold the aiding loop a fixed factor below the dominant wave frequency by
      scheduling theta. See Config::auto_theta_from_wave_freq.
    */
    void updateThetaSchedule_(R dt) {
        if (!cfg_.auto_theta_from_wave_freq || !initialized_) {
            return;
        }

        /*
          Track on a fixed reference heave, never on the observer's own
          position or velocity estimate.

          What has to be located is the peak of the displacement spectrum.
          Three candidate inputs, two of which do not work:

            Acceleration. The tracker would see omega^4 * S_z, which for a
            JONSWAP sea falls only as omega^-1 above the peak, so it locks
            well above it. On the Hs = 4 m record it reports 0.32 Hz against
            a true peak of 0.12 Hz and theta comes out nearly 3x too high.

            The observer's own heave estimate. This closes a positive
            feedback loop: theta sets how much low-frequency drift the
            estimate admits, drift captures the tracker, the tracker lowers
            theta, which admits more drift. On the same record it walks down
            to 0.054 Hz and stays there.

          So the tracker gets its own fixed reference channel: vertical
          inertial acceleration double-integrated through two leaky
          integrators at reference_heave_hp_hz. That is a second-order
          high-pass on displacement, so it has the displacement spectrum's
          shape in the wave band and rejects drift at -40 dB/decade below the
          corner, and above all it does not depend on theta, so there is no
          loop to run away.

          Sign is irrelevant to a frequency tracker.
        */
        const R a_z_ned = filter_.specificForceNED().z() + cfg_.gravity_mps2;
        if (isFinite_(a_z_ned)) {
            const R tau = R(1) /
                (R(6.283185307179586) * cfg_.reference_heave_hp_hz);

            ref_heave_v_ += dt * (a_z_ned - ref_heave_v_ / tau);
            ref_heave_z_ += dt * (ref_heave_v_ - ref_heave_z_ / tau);

            if (isFinite_(ref_heave_z_)) {
                marine_obs::detail::tracker_step(wave_freq_tracker_,
                                                 ref_heave_z_, dt);
            } else {
                ref_heave_v_ = R(0);
                ref_heave_z_ = R(0);
            }
        }

        theta_sched_accum_s_ += dt;
        if (theta_sched_accum_s_ < cfg_.theta_update_period_s) {
            return;
        }

        const R step = theta_sched_accum_s_;
        theta_sched_accum_s_ = R(0);

        last_wave_freq_hz_ =
            marine_obs::detail::tracker_get_frequency_hz<WaveFreqTracker, R>(
                wave_freq_tracker_);
        last_wave_freq_conf_ =
            marine_obs::detail::tracker_get_confidence<WaveFreqTracker, R>(
                wave_freq_tracker_);
        last_wave_freq_locked_ =
            marine_obs::detail::tracker_is_locked(wave_freq_tracker_);

        if (!isFinite_(last_wave_freq_hz_) || !(last_wave_freq_hz_ > R(0))) {
            return;
        }
        if (!(last_wave_freq_conf_ >= cfg_.theta_min_confidence)) {
            return;
        }

        constexpr R kTwoPi = R(6.283185307179586);

        const R target = clamp_(
            cfg_.theta_from_omega_gain * kTwoPi * last_wave_freq_hz_,
            cfg_.theta_min,
            cfg_.theta_max
        );

        /*
          Snap on first acquisition, lag afterwards.

          theta starts at the configured fixed value, which is nowhere near
          the sea being measured. Lagging into the target from there costs
          several time constants during which the aiding loop is misplaced,
          and that transient is a real error in the output. Once the tracker
          is confident for the first time there is no reason to approach the
          answer slowly, so take it, and reserve the lag for tracking a sea
          state that actually changes.
        */
        if (!theta_acquired_) {
            theta_acquired_ = true;
            theta_ = target;
        } else {
            const R alpha = (cfg_.theta_smooth_tau_s > R(0))
                ? clamp_(step / cfg_.theta_smooth_tau_s, R(0), R(1))
                : R(1);
            theta_ += alpha * (target - theta_);
        }

        filter_.setTheta(theta_);
    }

    static bool hasUsableGNSS_(const Aux& aux) {
        if constexpr (!WithGNSS) {
            (void)aux;
            return false;
        } else {
            return aux.gnss.valid;
        }
    }

    /*
      Track the slow mean of position and velocity so the reported values can
      be given with the unobservable constant removed.
      See Config::report_highpass_enabled.
    */
    void updateReportHighpass_(R dt, bool have_gnss) {
        if (!cfg_.report_highpass_enabled || !initialized_) {
            return;
        }

        const R tau = cfg_.report_highpass_tau_s;
        if (!(tau > R(0))) {
            return;
        }

        const Vec3 p = filter_.positionNED();
        const Vec3 v = filter_.velocityNED();
        if (!p.allFinite() || !v.allFinite()) {
            return;
        }

        if (!report_hp_primed_) {
            report_hp_primed_ = true;
            report_lpf_p_ = p;
            report_lpf_v_ = v;
            return;
        }

        const R alpha = clamp_(dt / tau, R(0), R(1));
        report_lpf_p_ += alpha * (p - report_lpf_p_);
        report_lpf_v_ += alpha * (v - report_lpf_v_);

        // With GNSS the horizontal DC is observable, so leave it in.
        if (have_gnss) {
            report_lpf_p_.x() = R(0);
            report_lpf_p_.y() = R(0);
            report_lpf_v_.x() = R(0);
            report_lpf_v_.y() = R(0);
        }
    }

    static bool isFinite_(R x) {
        using std::isfinite;
        return isfinite(x);
    }

    static R clamp_(R x, R lo, R hi) {
        return (x < lo) ? lo : ((x > hi) ? hi : x);
    }

    static Vec3 clampNorm_(const Vec3& v, R max_norm) {
        const R n = v.norm();
        if (!(max_norm > R(0)) || n <= max_norm || n <= R(1e-12)) {
            return v;
        }
        return v * (max_norm / n);
    }

    static Vec3 nedDown_() {
        Vec3 d;
        d << R(0), R(0), R(1);
        return d;
    }

    static Vec3 eulerFromQuat_(const Quat& q_nb) {
        const Mat3 Rnb = q_nb.toRotationMatrix();

        Vec3 e;
        e.x() = std::atan2(Rnb(2, 1), Rnb(2, 2));
        e.y() = -std::asin(clamp_(Rnb(2, 0), R(-1), R(1)));
        e.z() = std::atan2(Rnb(1, 0), Rnb(0, 0));
        return e;
    }

    bool accelNormOk_(const Vec3& f_b, R tol_frac) const {
        const R n = f_b.norm();
        return isFinite_(n) &&
               std::abs(n - cfg_.gravity_mps2) <
                   tol_frac * cfg_.gravity_mps2;
    }

    bool gyroOk_(const Vec3& gyro_b, R max_norm) const {
        const R n = gyro_b.norm();
        return isFinite_(n) && n < max_norm;
    }

    bool tiltTrimCurrentlyActive_() const {
        if (!cfg_.tilt_trim_enabled || !initialized_) {
            return false;
        }

        if (cfg_.tilt_trim_duration_s <= R(0)) {
            return true;
        }

        return time_since_init_s_ <= cfg_.tilt_trim_duration_s;
    }

    void updateInitialization_(R dt,
                               const Vec3& gyro_b_rad_s,
                               const Vec3& specific_force_b_mps2,
                               const Aux& aux)
    {
        if (cfg_.mahony_bootstrap_enabled) {
            updateMahonyBootstrap_(dt, gyro_b_rad_s, specific_force_b_mps2, aux);
            return;
        }

        updateAveragingInitialization_(dt, gyro_b_rad_s, specific_force_b_mps2, aux);
    }

    void updateAveragingInitialization_(R dt,
                                        const Vec3& gyro_b_rad_s,
                                        const Vec3& specific_force_b_mps2,
                                        const Aux& aux)
    {
        init_total_time_s_ += dt;
        init_all_acc_sum_ += specific_force_b_mps2;
        ++init_all_count_;

        if constexpr (Mag == NloMagType::Magnetometer) {
            if (aux.mag_valid) {
                init_all_mag_sum_ += aux.mag_b;
                ++init_all_mag_count_;
            }
        }

        const bool good =
            accelNormOk_(specific_force_b_mps2, cfg_.init_acc_norm_tol_frac) &&
            gyroOk_(gyro_b_rad_s, cfg_.init_gyro_max_rad_s);

        if (good) {
            init_good_time_s_ += dt;
            init_good_acc_sum_ += specific_force_b_mps2;
            ++init_good_count_;

            if constexpr (Mag == NloMagType::Magnetometer) {
                if (aux.mag_valid) {
                    init_good_mag_sum_ += aux.mag_b;
                    ++init_good_mag_count_;
                }
            }
        }

        if (init_good_time_s_ >= cfg_.init_required_good_time_s &&
            init_good_count_ > 0) {
            const Vec3 acc0 =
                init_good_acc_sum_ / static_cast<R>(init_good_count_);

            if (tryInitializeFromAverage_(acc0, aux, true)) {
                return;
            }
        }

        if (init_total_time_s_ >= cfg_.init_max_wait_s &&
            init_all_count_ > 0) {
            const Vec3 acc0 =
                init_all_acc_sum_ / static_cast<R>(init_all_count_);

            (void)tryInitializeFromAverage_(acc0, aux, false);
        }
    }

    bool tryInitializeFromAverage_(const Vec3& acc0,
                                   const Aux& aux,
                                   bool prefer_good_window)
    {
        bool ok = false;

        if constexpr (Mag == NloMagType::Compass) {
            if (aux.compass_valid) {
                ok = filter_.initializeFromAccelCompass(
                    acc0,
                    aux.compass_heading_rad
                );
            } else {
                ok = filter_.initializeFromAccel(acc0, cfg_.yaw_seed_rad);
            }
        } else if constexpr (Mag == NloMagType::Magnetometer) {
            Vec3 mag0 = Vec3::Zero();
            bool have_mag = false;

            if (prefer_good_window && init_good_mag_count_ > 0) {
                mag0 = init_good_mag_sum_ / static_cast<R>(init_good_mag_count_);
                have_mag = true;
            } else if (!prefer_good_window && init_all_mag_count_ > 0) {
                mag0 = init_all_mag_sum_ / static_cast<R>(init_all_mag_count_);
                have_mag = true;
            } else if (aux.mag_valid) {
                mag0 = aux.mag_b;
                have_mag = true;
            }

            if (have_mag) {
                ok = filter_.initializeFromAccelMag(acc0, mag0);
            }

            if (!ok) {
                ok = filter_.initializeFromAccel(acc0, cfg_.yaw_seed_rad);
            }
        } else {
            ok = filter_.initializeFromAccel(acc0, cfg_.yaw_seed_rad);
        }

        if (ok) {
            initialized_ = true;
            time_since_init_s_ = R(0);
            have_down_lpf_ = false;

            tilt_bias_b_ = filter_.gyroBiasBody();
            if constexpr (Mag == NloMagType::None) {
                tilt_bias_b_.z() = R(0);
            }
        }

        return ok;
    }

    bool tryStartMahonyBootstrap_(const Vec3& specific_force_b_mps2,
                                  const Aux& aux)
    {
        Filter tmp(cfg_.filter);
        tmp.reset();

        bool ok = false;

        if constexpr (Mag == NloMagType::Compass) {
            if (aux.compass_valid) {
                ok = tmp.initializeFromAccelCompass(
                    specific_force_b_mps2,
                    aux.compass_heading_rad
                );
            } else {
                ok = tmp.initializeFromAccel(
                    specific_force_b_mps2,
                    cfg_.yaw_seed_rad
                );
            }
        } else if constexpr (Mag == NloMagType::Magnetometer) {
            if (aux.mag_valid) {
                ok = tmp.initializeFromAccelMag(
                    specific_force_b_mps2,
                    aux.mag_b
                );
            }

            if (!ok) {
                ok = tmp.initializeFromAccel(
                    specific_force_b_mps2,
                    cfg_.yaw_seed_rad
                );
            }
        } else {
            ok = tmp.initializeFromAccel(
                specific_force_b_mps2,
                cfg_.yaw_seed_rad
            );
        }

        if (!ok) {
            return false;
        }

        boot_q_nb_ = tmp.quaternionBodyToNED();
        boot_q_nb_.normalize();

        boot_bias_b_.setZero();
        if constexpr (Mag == NloMagType::None) {
            boot_bias_b_.z() = R(0);
        }

        boot_started_ = true;
        boot_have_down_lpf_ = false;
        boot_time_s_ = R(0);
        boot_good_time_s_ = R(0);

        return true;
    }

    void updateMahonyBootstrap_(R dt,
                                const Vec3& gyro_b_rad_s,
                                const Vec3& specific_force_b_mps2,
                                const Aux& aux)
    {
        if (!boot_started_) {
            if (!tryStartMahonyBootstrap_(specific_force_b_mps2, aux)) {
                return;
            }
        }

        boot_time_s_ += dt;
        init_total_time_s_ = boot_time_s_;

        Vec3 correction_b = Vec3::Zero();

        const bool accel_ok =
            accelNormOk_(specific_force_b_mps2, cfg_.mahony_acc_norm_tol_frac);

        const bool gyro_ok =
            gyroOk_(gyro_b_rad_s, cfg_.mahony_gyro_max_rad_s);

        if (accel_ok && gyro_ok) {
            const R fn = specific_force_b_mps2.norm();

            if (fn > R(1e-9) && isFinite_(fn)) {
                const Vec3 down_meas_b =
                    (-specific_force_b_mps2 / fn).normalized();

                const R alpha = clamp_(
                    dt / std::max(cfg_.mahony_down_lpf_tau_s, R(1e-3)),
                    R(0),
                    R(1)
                );

                if (!boot_have_down_lpf_) {
                    boot_down_lpf_b_ = down_meas_b;
                    boot_have_down_lpf_ = true;
                } else {
                    const Vec3 blended =
                        (R(1) - alpha) * boot_down_lpf_b_ +
                        alpha * down_meas_b;

                    const R bn = blended.norm();
                    if (bn > R(1e-9) && isFinite_(bn)) {
                        boot_down_lpf_b_ = blended / bn;
                    }
                }

                const Mat3 R_nb = boot_q_nb_.toRotationMatrix();
                const Vec3 down_hat_b =
                    (R_nb.transpose() * nedDown_()).normalized();

                /*
                  q_new = q * dq(delta_b)
                  err direction chosen to rotate predicted down toward measured down.
                */
                Vec3 err_b = boot_down_lpf_b_.cross(down_hat_b);

                if constexpr (Mag == NloMagType::None) {
                    err_b.z() = R(0);
                }

                correction_b = cfg_.mahony_twoKp * R(0.5) * err_b;

                if (cfg_.mahony_twoKi > R(0)) {
                    boot_bias_b_ +=
                        (-cfg_.mahony_twoKi * R(0.5) * dt) * err_b;

                    if constexpr (Mag == NloMagType::None) {
                        boot_bias_b_.z() = R(0);
                    }

                    boot_bias_b_ =
                        clampNorm_(boot_bias_b_, cfg_.mahony_bias_limit_rad_s);
                }

                boot_good_time_s_ += dt;
                init_good_time_s_ = boot_good_time_s_;
            }
        }

        Vec3 omega_b = gyro_b_rad_s - boot_bias_b_ + correction_b;

        if constexpr (Mag == NloMagType::None) {
            boot_bias_b_.z() = R(0);
        }

        integrateBootQuatRight_(dt, omega_b);

        const bool enough_good =
            boot_good_time_s_ >= cfg_.mahony_bootstrap_min_time_s;

        const bool forced =
            boot_time_s_ >= cfg_.mahony_bootstrap_max_time_s;

        if (enough_good || forced) {
            filter_ = Filter(cfg_.filter);
            filter_.reset();

            filter_.setQuaternionBodyToNED(boot_q_nb_);

            Vec3 b0 = boot_bias_b_;
            if constexpr (Mag == NloMagType::None) {
                b0.z() = R(0);
            }

            filter_.setGyroBiasBody(b0);
            tilt_bias_b_ = b0;

            initialized_ = true;
            time_since_init_s_ = R(0);
            have_down_lpf_ = false;
        }
    }

    void integrateBootQuatRight_(R dt, const Vec3& omega_b)
    {
        const R omega_norm = omega_b.norm();
        const R angle = omega_norm * dt;

        Quat dq;
        if (angle < R(1e-8) || omega_norm < R(1e-12)) {
            const Vec3 half = R(0.5) * dt * omega_b;
            dq = Quat(R(1), half.x(), half.y(), half.z());
        } else {
            const Vec3 axis = omega_b / omega_norm;
            const R half_angle = R(0.5) * angle;
            const R s = std::sin(half_angle);

            dq = Quat(std::cos(half_angle),
                      s * axis.x(),
                      s * axis.y(),
                      s * axis.z());
        }

        boot_q_nb_ = boot_q_nb_ * dq;
        boot_q_nb_.normalize();
    }

    void applyOptionalTiltTrim_(R dt,
                                const Vec3& gyro_b_rad_s,
                                const Vec3& specific_force_b_mps2)
    {
        last_tilt_err_norm_ = R(0);
        last_tilt_correction_norm_ = R(0);

        if (!tiltTrimCurrentlyActive_()) {
            return;
        }

        if (!accelNormOk_(specific_force_b_mps2, cfg_.tilt_trim_acc_norm_tol_frac)) {
            return;
        }

        if (!gyroOk_(gyro_b_rad_s, cfg_.tilt_trim_gyro_max_rad_s)) {
            return;
        }

        const R fn = specific_force_b_mps2.norm();
        if (!(fn > R(1e-9)) || !isFinite_(fn)) {
            return;
        }

        const Vec3 down_meas_b =
            (-specific_force_b_mps2 / fn).normalized();

        const R lpf_tau = std::max(cfg_.tilt_trim_acc_lpf_tau_s, R(1e-3));
        const R alpha = clamp_(dt / lpf_tau, R(0), R(1));

        if (!have_down_lpf_) {
            down_lpf_b_ = down_meas_b;
            have_down_lpf_ = true;
        } else {
            const Vec3 blended =
                (R(1) - alpha) * down_lpf_b_ + alpha * down_meas_b;
            const R bn = blended.norm();
            if (bn > R(1e-9) && isFinite_(bn)) {
                down_lpf_b_ = blended / bn;
            }
        }

        const Mat3 R_nb = filter_.rotationBodyToNED();
        const Vec3 down_hat_b =
            (R_nb.transpose() * nedDown_()).normalized();

        Vec3 err_b = down_lpf_b_.cross(down_hat_b);

        if constexpr (Mag == NloMagType::None) {
            err_b.z() = R(0);
        }

        const R en = err_b.norm();
        if (!(en > R(1e-12)) || !isFinite_(en)) {
            return;
        }

        last_tilt_err_norm_ = en;

        if (cfg_.tilt_bias_enabled) {
            tilt_bias_b_ += (-cfg_.tilt_bias_ki * dt) * err_b;

            if constexpr (Mag == NloMagType::None) {
                tilt_bias_b_.z() = R(0);
            }

            tilt_bias_b_ =
                clampNorm_(tilt_bias_b_, cfg_.tilt_bias_limit_rad_s);
        }

        const R gain = R(1) / std::max(cfg_.tilt_trim_tau_s, R(1e-3));
        const R max_step = cfg_.tilt_trim_max_rate_rad_s * dt;

        Vec3 delta_b = gain * dt * err_b;
        const R dn = delta_b.norm();
        if (dn > max_step && dn > R(1e-12)) {
            delta_b *= max_step / dn;
        }

        last_tilt_correction_norm_ = delta_b.norm();

        applyBodyDeltaQuaternion_(delta_b);
    }

    void applyBodyDeltaQuaternion_(const Vec3& delta_b) {
        const R angle = delta_b.norm();

        Quat dq;
        if (angle < R(1e-9)) {
            const Vec3 half = R(0.5) * delta_b;
            dq = Quat(R(1), half.x(), half.y(), half.z());
        } else {
            const Vec3 axis = delta_b / angle;
            const R half_angle = R(0.5) * angle;
            const R s = std::sin(half_angle);
            dq = Quat(std::cos(half_angle),
                      s * axis.x(),
                      s * axis.y(),
                      s * axis.z());
        }

        Quat q = filter_.quaternionBodyToNED() * dq;
        q.normalize();
        filter_.setQuaternionBodyToNED(q);
    }
};
