#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  LiveEntranceCertificate - constructive basin certificate for the
  proxy-to-Live handoff of the OU-III MEKF.

  What this is
  ------------
  doc/kalman_ou_iii/w3d-live-handoff-certificate.tex-part proves a block-local
  ISS bound for the Live error dynamics in the filter's own covariance
  (Riccati) metric, and gives a sufficient entrance inequality

      M_xi_xi * ||e_H||_{xi,H} + M_xi_ell * ||e_H||_{ell,H}
          + disturbance_margin  <  (1 - nu) * r,
      nu = c_eff * r,

  where the two norms are the sensitive and kinematic block norms in the
  filter's own metric at the handoff sample, r is the sensitive-tube radius,
  and c_eff is the small-gain slope assembled from the interval gain
  envelopes.  Only the sensitive block has to lie in the tube; the kinematic
  one is charged through its own transient gain, which is what Phase 1's block
  separation buys.  This header evaluates that inequality with the numbers a
  running system actually has at the handoff sample.

  What it is not
  --------------
  It is not a gate reverse-engineered from what startup happens to achieve,
  and it is not a boolean.  Every term is exposed, because the question a
  failed certification has to answer is *which* term spent the budget.  The
  four kinds of number are kept apart on purpose:

    * configured physical envelope - declared, not measured, not provable
      onboard (accelerometer calibration residual, magnetic reference error,
      gyro turn-on bias).  These are the assumptions the theorem is
      conditional on.
    * measured startup statistic - taken from the bootstrap's own observables
      (the gravity-residual hold and its branch, whether a magnetic gauge was
      accepted, and the seed marginals actually installed).
    * analytically constructed bound - what the paper's inequalities make of
      the two above.
    * certified margin - the final number, and the only one that decides.

  A handoff that fails this test is still allowed to go Live; it is simply not
  covered by the theorem, and the state machine records that difference.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <cstdint>
#include <limits>

namespace ocean_imu {
namespace kalman {
namespace ou3 {

// Where a handoff bound comes from.  The tag travels with the number because
// the theorem is exactly as strong as the weakest tag it rests on.
enum class BoundSource : uint8_t {
    Exact,     // zero, or otherwise fixed, by construction
    Envelope,  // deterministic given a configured physical operating envelope
    Measured,  // estimated from startup observables; not a truth bound
    External   // still an external hypothesis, not discharged onboard
};

inline const char* boundSourceName(BoundSource s) {
    switch (s) {
        case BoundSource::Exact:    return "exact";
        case BoundSource::Envelope: return "envelope";
        case BoundSource::Measured: return "measured";
        case BoundSource::External: return "external";
    }
    return "?";
}

struct HandoffBound {
    float       value  = std::numeric_limits<float>::quiet_NaN();  // physical units
    BoundSource source = BoundSource::External;

    bool valid() const { return std::isfinite(value) && value >= 0.0f; }
};

// Why a certificate was refused.  Ordered as the evaluator tests them, so the
// first failure is the one reported.
enum class LiveCertFailure : uint8_t {
    None = 0,
    NotAttempted,        // no handoff evaluated yet
    ConstantsMissing,    // interval constants or envelope not usable
    GateNotSatisfied,    // quality gate / aligned branch / tuner readiness
    AttitudeUnbounded,   // no finite yaw bound (no magnetic gauge)
    BoundNotConstructed, // some block bound is not finite
    SeedIncoherent,      // seed covariance is not usable as a metric
    BasinExceeded        // the inequality itself failed
};

inline const char* liveCertFailureName(LiveCertFailure f) {
    switch (f) {
        case LiveCertFailure::None:                return "certified";
        case LiveCertFailure::NotAttempted:        return "not-attempted";
        case LiveCertFailure::ConstantsMissing:    return "constants-missing";
        case LiveCertFailure::GateNotSatisfied:    return "gate-not-satisfied";
        case LiveCertFailure::AttitudeUnbounded:   return "attitude-unbounded";
        case LiveCertFailure::BoundNotConstructed: return "bound-not-constructed";
        case LiveCertFailure::SeedIncoherent:      return "seed-incoherent";
        case LiveCertFailure::BasinExceeded:       return "basin-exceeded";
    }
    return "?";
}

// Where a set of interval constants comes from, which decides what a passing
// certificate is evidence of.
enum class CertificateSource : uint8_t {
    // Verified arithmetic over the declared compact domain of
    // doc/kalman_ou_iii/w3d-computer-assisted-live-basin.tex-part, evaluated by
    // tools/ou_live_basin_interval_proof.py.  No sea spectrum, no replay, no
    // fitted margin.  A certificate issued against these is a theorem.
    Analytical,
    // Suprema measured by tests/kalman_ou_iii/live_basin_diagnostic.cpp over
    // the eight committed reference operating points.  Vastly tighter, and
    // conditional on the running trajectory actually satisfying the same
    // envelopes -- which nothing onboard checks.  A certificate issued against
    // these is evidence, not a proof.
    IntervalMeasured
};

inline const char* certificateSourceName(CertificateSource s) {
    switch (s) {
        case CertificateSource::Analytical:      return "analytical";
        case CertificateSource::IntervalMeasured: return "interval-measured";
    }
    return "?";
}

/*
  Interval constants of the metric block-local ISS theorem.

  Two sets exist and they differ by roughly 290 orders of magnitude, so which
  one is in force is part of the certificate rather than a detail.

  The default is the analytical set.  Its numbers come from the computer-
  assisted closure, which proves over a declared compact domain exactly the
  three things this header needs and Phase 3 previously had to assume: that the
  horizon contraction margin 1-chi_H is strictly positive, that the covariance
  is bounded above, and that the one-step nonlinear constants are finite.  They
  are quoted at the precision the proof program prints and are checked against
  it by tests/validation/test_ou_analytic_constants_match_proof.py.

  The measured set is available through intervalMeasured() for deployments that
  have verified their own interval envelopes.  It is the honest, useful, and
  unproved number; the analytical set is the honest, useless, and proved one.
  Neither is allowed to be reported as the other.
*/
struct LiveBasinConstants {
    CertificateSource source = CertificateSource::Analytical;

    // Metric transition constant, sensitive block against sensitive block.
    //
    // Exactly one under either source, and for the same reason: the Joseph
    // identity makes every intermediate Riccati step nonexpansive in the
    // covariance metric, so the metric norm over a horizon is largest at the
    // horizon's own start.  This is the constant that was ~566 in the Phase-2
    // Euclidean formulation.
    double M_xi_xi = 1.0;

    // Transient gain from the kinematic block into the sensitive one, bounded
    // by the same argument.  Phase 1's separation lives here: v and p are
    // charged through this gain rather than being required to lie in the tube.
    double M_xi_ell = 1.0;

    // Small-gain slope: nu = c_eff * r.
    //
    // The lifted analytical recursion is R_{n+1} <= chi_H R_n + C_H R_n^2 at
    // horizon boundaries, so the tube R <= r is invariant exactly when
    // C_H r <= 1 - chi_H.  Hence c_eff = C_H / (1 - chi_H), which is the
    // reciprocal of the proof program's riccati_nonlinear_radius_lower.
    //
    //   C_H       <= 2.413134956828e+202
    //   1 - chi_H >= 2.608762190650e-87
    //
    // The quotient is enormous, and that is the honest state of the broad-box
    // proof: the article's own conclusion is that the remaining obstacle is
    // conservatism of the global covariance and measurement bounds, not strict
    // stability.  See intervalMeasured() for what the deployed schedule
    // actually exhibits.
    double c_eff = 9.250114730570e+288;

    // Metric disturbance allowance already charged against the budget.
    double disturbance_margin = 0.0;

    // Largest tube radius the rotation-vector remainder expansion is valid on,
    // r_xi / sqrt(pbar) in metric units (the proof program's
    // riccati_chart_radius_lower).  Under the analytical source it does not
    // bind: the operating radius 1/(2 c_eff) is ~5.4e-290.
    double r_xi = 5.181025192391e-92;

    // Interval constants measured by live_basin_diagnostic.cpp over the eight
    // committed reference operating points, worst case.
    //
    // These are ~290 orders of magnitude tighter than the analytical set and
    // they are not a proof: they hold for the schedule the diagnostic replayed,
    // and nothing onboard verifies that a running trajectory satisfies the same
    // envelopes.  A deployment that has established its own envelopes may pin
    // them; the certificate then reports source = IntervalMeasured and must not
    // be described as theorem-certified.
    static LiveBasinConstants intervalMeasured() {
        LiveBasinConstants k;
        k.source = CertificateSource::IntervalMeasured;
        // Exactly one, as under the analytical source and for the same
        // reason.  The 1.02 that an earlier draft carried here was
        // rho_H^{-(H-1)}, the overshoot of insisting on an M rho^{k-j}
        // envelope; the tube argument needs only the supremum over the
        // interval, and metric monotonicity gives that as one.  Keeping both
        // sources at one also makes the inequality's left-hand side
        // independent of which constants are in force, so a single evaluated
        // handoff can be scored against both.
        k.M_xi_xi = 1.0;
        k.M_xi_ell = 1.0;
        // Worst of the eight points: 0.68 at the 0.27 m reference sea, 29.3 at
        // the 8.5 m one.  The worst is the default because a certificate has to
        // be sound over the declared envelope, not over the sea that happens to
        // be running.
        k.c_eff = 29.3;
        // Measured 88 to 327; does not bind against an operating radius of
        // at most 0.74.
        k.r_xi = 88.0;
        return k;
    }

    bool valid() const {
        return std::isfinite(M_xi_xi) && M_xi_xi >= 1.0 &&
               std::isfinite(M_xi_ell) && M_xi_ell >= 0.0 &&
               std::isfinite(c_eff) && c_eff > 0.0 &&
               std::isfinite(disturbance_margin) && disturbance_margin >= 0.0 &&
               std::isfinite(r_xi) && r_xi > 0.0;
    }
};

/*
  Declared physical operating envelope.

  Nothing here is measured or provable onboard.  Each entry is the hypothesis
  a deployment makes about its own hardware and sea, and the theorem's
  conclusion is conditional on it.  They are kept apart from the measured
  statistics so a reader of a certificate can see which is which.
*/
struct LiveEnvelope {
    // Angle between the low-passed measured rest-specific-force direction and
    // true world down during the accepted hold.  A small accelerometer
    // residual does NOT by itself prove a small true tilt in wave
    // acceleration: what the gate sees is the residual of an average, and the
    // part of the orbital acceleration that survives that average is exactly
    // what this number declares.  eta_g in the paper.
    float gravity_direction_rad = 0.05f;

    // Relative error of the learned horizontal magnetic reference.  A
    // self-consistent learned reference does not prove absolute heading
    // accuracy; this declares how far the reference may be from true north.
    float mag_reference_rel_error = 0.06f;

    // Residual gyro bias about the gravity direction.  The startup proxy
    // observes only the two components perpendicular to measured down, so the
    // third is a turn-on-bias declaration, not an estimate.
    float gyro_bias_axial_rad_s = 3.0e-4f;

    // Convergence residual of the proxy's perpendicular integral bias
    // estimate at handoff.
    float gyro_bias_perp_rad_s = 2.0e-4f;

    // Residual accelerometer bias after the deployed temperature correction,
    // i.e. a bound on the TRUE b_a error, not on the estimate.  Nothing
    // onboard measures this; it is a calibration-envelope declaration.
    float accel_bias_ms2 = 0.05f;

    // Accelerometer measurement error (noise plus scale-factor residual) used
    // when converting a levelled specific force into an a_w seed.
    float accel_meas_ms2 = 0.10f;

    // Wave-orbital velocity and displacement the platform may carry at the
    // handoff sample when no translational bootstrap is available.  Both scale
    // with the sea, so a single declaration has to cover the roughest sea in
    // the deployment envelope: the defaults cover the 8.5 m reference sea.
    //
    // Measured rather than assumed, these two turn out not to need a
    // bootstrap.  At handoff the MEKF's v and p are still at zero and the true
    // wave-orbital values are what the error is, and over the reference seas
    // that error stays inside the envelope below without any startup
    // integration at all -- while the kinematic block enters the basin
    // inequality through M_xi_ell rather than through the tube, so shrinking
    // it buys correspondingly little.  A translational bootstrap is therefore
    // not implemented; this is the measurement that says it is not needed, not
    // an omission.
    float translational_velocity_ms = 4.0f;
    float translational_position_m  = 8.0f;

    bool valid() const {
        return std::isfinite(gravity_direction_rad)   && gravity_direction_rad   >= 0.0f &&
               std::isfinite(mag_reference_rel_error) && mag_reference_rel_error >= 0.0f &&
               std::isfinite(gyro_bias_axial_rad_s)   && gyro_bias_axial_rad_s   >= 0.0f &&
               std::isfinite(gyro_bias_perp_rad_s)    && gyro_bias_perp_rad_s    >= 0.0f &&
               std::isfinite(accel_bias_ms2)          && accel_bias_ms2          >= 0.0f &&
               std::isfinite(accel_meas_ms2)          && accel_meas_ms2          >= 0.0f &&
               std::isfinite(translational_velocity_ms) &&
               translational_velocity_ms >= 0.0f &&
               std::isfinite(translational_position_m) &&
               translational_position_m >= 0.0f;
    }
};

// Startup observables the wrapper hands the evaluator.  These are measured,
// and they are labelled measured.
struct LiveHandoffObservables {
    bool  gate_satisfied    = false;  // quality hold, aligned branch, tuner ready
    bool  aligned_branch    = false;
    bool  yaw_gauged        = false;  // a magnetic reference was accepted
    float gravity_align_sin = std::numeric_limits<float>::quiet_NaN();

    // Marginal standard deviations of the seed covariance actually installed,
    // in physical units, block by block.  These define the metric the
    // certificate is evaluated in, so they are read back from the filter
    // rather than assumed.
    float sigma_theta = std::numeric_limits<float>::quiet_NaN();  // rad
    float sigma_bg    = std::numeric_limits<float>::quiet_NaN();  // rad/s
    float sigma_S     = std::numeric_limits<float>::quiet_NaN();  // m s
    float sigma_aw    = std::numeric_limits<float>::quiet_NaN();  // m/s^2
    float sigma_ba    = std::numeric_limits<float>::quiet_NaN();  // m/s^2
    float sigma_v     = std::numeric_limits<float>::quiet_NaN();  // m/s
    float sigma_p     = std::numeric_limits<float>::quiet_NaN();  // m

    // Bootstrap seeds actually installed, so the validation harness can
    // compare a predicted bound against the true error of the same quantity.
    bool  seeded_gyro_bias   = false;
    bool  seeded_world_accel = false;
    bool  seeded_translation = false;
};

/*
  The certificate itself.  Everything the inequality consumed is kept, so a
  rejected handoff can be read rather than guessed at.
*/
struct LiveEntranceCertificate {
    // ---- constructed physical bounds, each with its provenance ----
    HandoffBound attitude;     // rad
    HandoffBound gyro_bias;    // rad/s
    HandoffBound integral_S;   // m s
    HandoffBound world_accel;  // m/s^2
    HandoffBound accel_bias;   // m/s^2
    HandoffBound velocity;     // m/s
    HandoffBound position;     // m

    // ---- the two halves of the attitude bound, kept apart because they are
    //      two different assumptions ----
    float eta_gravity_rad = std::numeric_limits<float>::quiet_NaN();
    float eta_heading_rad = std::numeric_limits<float>::quiet_NaN();

    // ---- metric (dimensionless) block norms at the handoff sample ----
    float sensitive_metric_norm = std::numeric_limits<float>::quiet_NaN();
    float kinematic_metric_norm = std::numeric_limits<float>::quiet_NaN();

    // ---- the inequality ----
    //
    // Double, not float: under the analytical source c_eff is ~9.25e288 and
    // the tube radius ~5.4e-290, neither of which a float can hold.  The
    // physical bounds above stay float because they are physical.
    double nonlinear_radius        = std::numeric_limits<double>::quiet_NaN(); // r
    double small_gain_nu           = std::numeric_limits<double>::quiet_NaN(); // nu
    double linear_certificate_gain = std::numeric_limits<double>::quiet_NaN(); // M_xi_xi
    double kinematic_gain          = std::numeric_limits<double>::quiet_NaN(); // M_xi_ell
    double disturbance_margin      = std::numeric_limits<double>::quiet_NaN();
    double basin_lhs               = std::numeric_limits<double>::quiet_NaN();
    double basin_rhs               = std::numeric_limits<double>::quiet_NaN();
    double margin                  = std::numeric_limits<double>::quiet_NaN();

    // Which set of interval constants decided this certificate.  A pass under
    // IntervalMeasured is evidence about the deployed schedule; only a pass
    // under Analytical is a theorem.
    CertificateSource source = CertificateSource::Analytical;

    bool            certified = false;
    LiveCertFailure failure   = LiveCertFailure::NotAttempted;

    const char* failureName() const { return liveCertFailureName(failure); }
    const char* sourceName() const { return certificateSourceName(source); }

    // True only for a certificate that both passed and was decided by the
    // analytical constants.  This is the predicate anything reporting
    // "covered by the semiglobal theorem" has to use.
    bool theoremCertified() const {
        return certified && source == CertificateSource::Analytical;
    }

    // Weakest provenance among the bounds the inequality used.  A certificate
    // resting on an External bound is not a proof of anything the system
    // checked.
    BoundSource weakestSource() const {
        const HandoffBound* all[] = {&attitude, &gyro_bias, &integral_S,
                                     &world_accel, &accel_bias,
                                     &velocity, &position};
        BoundSource worst = BoundSource::Exact;
        for (const HandoffBound* b : all) {
            if (static_cast<uint8_t>(b->source) > static_cast<uint8_t>(worst)) {
                worst = b->source;
            }
        }
        return worst;
    }
};

// Right-hand side of the entrance inequality for a given constant set: the
// largest sensitive-block metric norm the tube can absorb.  Exposed because
// the left-hand side does not depend on the constants at all -- it is built
// from the handoff bounds and the installed seed -- so one evaluated handoff
// can be scored against both constant sets without being re-run.
inline double liveBasinBudget(const LiveBasinConstants& k) {
    if (!k.valid()) return std::numeric_limits<double>::quiet_NaN();
    const double r = std::fmin(0.5 / k.c_eff, k.r_xi);
    return (1.0 - k.c_eff * r) * r - k.disturbance_margin;
}

namespace cert_detail {

inline float safeRatio(float num, float den) {
    if (!std::isfinite(num) || !std::isfinite(den) || !(den > 0.0f)) {
        return std::numeric_limits<float>::infinity();
    }
    return num / den;
}

inline float asinClamped(float x) {
    if (!std::isfinite(x)) return std::numeric_limits<float>::quiet_NaN();
    if (x >= 1.0f) return 1.5707963f;
    if (x <= -1.0f) return -1.5707963f;
    return std::asin(x);
}

}  // namespace cert_detail

/*
  Build the certificate.

  gravity_align_max_sin is the deployed gate threshold, passed in rather than
  duplicated, so the certificate and the gate cannot drift apart.  The tilt
  half of the attitude bound is asin(threshold) + eta_g and is available only
  on the aligned branch, which is why the branch flag is a hard precondition
  rather than one more additive term.
*/
inline LiveEntranceCertificate evaluateLiveEntrance(
    const LiveHandoffObservables& obs,
    const LiveEnvelope&           env,
    const LiveBasinConstants&     k,
    float                         gravity_align_max_sin)
{
    LiveEntranceCertificate c;
    c.source                  = k.source;
    c.linear_certificate_gain = k.M_xi_xi;
    c.kinematic_gain          = k.M_xi_ell;
    c.disturbance_margin      = k.disturbance_margin;

    if (!k.valid() || !env.valid() ||
        !std::isfinite(gravity_align_max_sin) ||
        !(gravity_align_max_sin > 0.0f) || !(gravity_align_max_sin < 1.0f)) {
        c.failure = LiveCertFailure::ConstantsMissing;
        return c;
    }

    if (!obs.gate_satisfied || !obs.aligned_branch) {
        c.failure = LiveCertFailure::GateNotSatisfied;
        return c;
    }

    // ---- B1 attitude -------------------------------------------------
    // Tilt: the gate certifies the aligned branch and a sine residual, so the
    // angle between the predicted and the low-passed measured direction is at
    // most asin(threshold).  eta_g carries that measured direction the rest of
    // the way to true down and is a declared envelope, not an observation.
    c.eta_gravity_rad = env.gravity_direction_rad;
    const float tilt_bound =
        cert_detail::asinClamped(gravity_align_max_sin) + c.eta_gravity_rad;

    if (!obs.yaw_gauged) {
        // Without a magnetic gauge the heading error is bounded by nothing the
        // system observed.  The theorem needs a finite total attitude bound;
        // there is none, and saying so is the honest outcome.
        c.attitude.value  = std::numeric_limits<float>::infinity();
        c.attitude.source = BoundSource::External;
        c.failure         = LiveCertFailure::AttitudeUnbounded;
        return c;
    }

    // eta_m = asin(delta_B_h / B_h): the declared relative reference error.
    c.eta_heading_rad = cert_detail::asinClamped(env.mag_reference_rel_error);
    if (!std::isfinite(c.eta_heading_rad)) {
        c.failure = LiveCertFailure::BoundNotConstructed;
        return c;
    }

    c.attitude.value  = tilt_bound + c.eta_heading_rad;
    c.attitude.source = BoundSource::Envelope;

    // ---- B2 gyro bias -------------------------------------------------
    // The proxy's integral term supplies the two components perpendicular to
    // measured down; the axial one is a declaration.
    c.gyro_bias.value = std::sqrt(
        env.gyro_bias_perp_rad_s * env.gyro_bias_perp_rad_s +
        env.gyro_bias_axial_rad_s * env.gyro_bias_axial_rad_s);
    c.gyro_bias.source = BoundSource::Envelope;

    // ---- B4 integral state -------------------------------------------
    // Zero by the integration-epoch choice; see
    // Kalman3D_Wave_OU_III::reset_integral_epoch().
    c.integral_S.value  = 0.0f;
    c.integral_S.source = BoundSource::Exact;

    // ---- B6 accelerometer bias ---------------------------------------
    c.accel_bias.value  = env.accel_bias_ms2;
    c.accel_bias.source = BoundSource::Envelope;

    // ---- B3 world acceleration ---------------------------------------
    // The seed is the proxy-levelled specific force, so its error is the
    // levelling error times the specific-force magnitude, plus the
    // accelerometer's own error and its bias.  Without the seed the state
    // enters Live at zero and the bound is the acceleration envelope itself.
    constexpr float kG = 9.80665f;
    if (obs.seeded_world_accel) {
        c.world_accel.value =
            kG * std::sin(std::fmin(c.attitude.value, 1.5707963f)) +
            env.accel_meas_ms2 + c.accel_bias.value;
    } else {
        c.world_accel.value = std::numeric_limits<float>::infinity();
    }
    c.world_accel.source = BoundSource::Envelope;

    // ---- B5 velocity and position ------------------------------------
    c.velocity.value  = env.translational_velocity_ms;
    c.velocity.source = BoundSource::Envelope;
    c.position.value  = env.translational_position_m;
    c.position.source = BoundSource::Envelope;

    const HandoffBound* all[] = {&c.attitude, &c.gyro_bias, &c.integral_S,
                                 &c.world_accel, &c.accel_bias,
                                 &c.velocity, &c.position};
    for (const HandoffBound* b : all) {
        if (!b->valid()) {
            c.failure = LiveCertFailure::BoundNotConstructed;
            return c;
        }
    }

    // ---- metric norms at the handoff sample --------------------------
    // The seed covariance is block diagonal by construction (every seeder
    // drops its cross-covariances), so the metric norm of the handoff error
    // is the root sum of squares of per-block ratios.  This is a statement
    // about the tube coordinates, not a claim that the covariance bounds the
    // truth error: the numerators are the constructed bounds above.
    const float rt = cert_detail::safeRatio(c.attitude.value,    obs.sigma_theta);
    const float rg = cert_detail::safeRatio(c.gyro_bias.value,   obs.sigma_bg);
    const float rS = cert_detail::safeRatio(c.integral_S.value,  obs.sigma_S);
    const float ra = cert_detail::safeRatio(c.world_accel.value, obs.sigma_aw);
    const float rb = cert_detail::safeRatio(c.accel_bias.value,  obs.sigma_ba);
    const float rv = cert_detail::safeRatio(c.velocity.value,    obs.sigma_v);
    const float rp = cert_detail::safeRatio(c.position.value,    obs.sigma_p);

    const float sens[] = {rt, rg, rS, ra, rb};
    const float kin[]  = {rv, rp};
    float ss = 0.0f;
    for (float x : sens) {
        if (!std::isfinite(x)) { c.failure = LiveCertFailure::SeedIncoherent; return c; }
        ss += x * x;
    }
    float ks = 0.0f;
    for (float x : kin) {
        if (!std::isfinite(x)) { c.failure = LiveCertFailure::SeedIncoherent; return c; }
        ks += x * x;
    }
    c.sensitive_metric_norm = std::sqrt(ss);
    c.kinematic_metric_norm = std::sqrt(ks);

    // ---- the basin inequality -----------------------------------------
    // The budget (1 - c_eff r) r is maximised at r = 1/(2 c_eff); the tube is
    // capped at r_xi, where the remainder expansion stops being valid.
    const double r_opt = 0.5 / k.c_eff;
    c.nonlinear_radius = std::fmin(r_opt, k.r_xi);
    c.small_gain_nu    = k.c_eff * c.nonlinear_radius;

    // Only the sensitive block has to lie in the tube; the kinematic block is
    // charged through its own projected transient gain, which is Phase 1's
    // separation carried into the metric.
    c.basin_lhs = k.M_xi_xi * double(c.sensitive_metric_norm)
                + k.M_xi_ell * double(c.kinematic_metric_norm)
                + k.disturbance_margin;
    c.basin_rhs = (1.0 - c.small_gain_nu) * c.nonlinear_radius;
    c.margin    = c.basin_rhs - c.basin_lhs;

    if (!(c.small_gain_nu < 1.0) || !(c.margin > 0.0)) {
        c.failure = LiveCertFailure::BasinExceeded;
        return c;
    }

    c.certified = true;
    c.failure   = LiveCertFailure::None;
    return c;
}

}  // namespace ou3
}  // namespace kalman
}  // namespace ocean_imu
