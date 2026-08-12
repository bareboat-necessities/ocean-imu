/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    ContinuousMagHardIronEstimator.

    MagAutoTuner's own tests establish what a single startup window can and
    cannot do with a body-fixed magnetic offset.  These are about what changes
    when the accumulation is never closed, and about the one thing that does
    not change no matter how long it runs.

    A hull that only rolls and pitches never separates the offset from the
    world field cleanly, because the sensor's other errors -- soft iron and
    misalignment -- are modulated by the same tilt the offset is read from and
    are aliased into the fit.  That aliasing is not noise: it does not average
    away, and its size does not fall as the record lengthens.  So the tests
    below check the two properties that make a continuously running fit safe to
    leave switched on rather than merely accurate on a good day:

      - the regularisation returns a bounded fraction of the fit, and the
        fraction is set by the sensor rather than by the sea state, so a big
        sea does not license a confident wrong answer;
      - the reference belonging to a partly applied offset is available from
        the same statistics, so the caller can never pair one with the other.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <numbers>
#include <random>
#include <string>

#include <Eigen/Dense>

#include "tuner/ContinuousMagHardIronEstimator.h"

namespace {

using Eigen::Matrix3f;
using Eigen::Quaternionf;
using Eigen::Vector3f;

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::printf("%-62s %s\n", what.c_str(), ok ? "ok" : "FAILED");
    if (!ok) ++g_failures;
}

constexpr float kDt = 1.0f / 25.0f;

// NED: north, east, down.  Dip and magnitude of the reference field the OU
// simulator uses, so the numbers here are the ones the filter actually sees.
const Vector3f kFieldNed(52.0f * std::cos(66.5f * float(std::numbers::pi_v<float>) / 180.0f),
                         0.0f,
                         52.0f * std::sin(66.5f * float(std::numbers::pi_v<float>) / 180.0f));

Quaternionf tiltQuat(float roll_rad, float pitch_rad)
{
    return Quaternionf(Eigen::AngleAxisf(pitch_rad, Vector3f::UnitY()) *
                       Eigen::AngleAxisf(roll_rad, Vector3f::UnitX()));
}

struct Sea {
    float roll_amp_rad;
    float pitch_amp_rad;
    float roll_hz = 0.11f;
    float pitch_hz = 0.083f;
};

// One replay of a hull rolling and pitching at a fixed heading, with the
// sensor errors the simulator injects: a body-fixed offset, a soft-iron and
// misalignment matrix, and white noise.
struct Replay {
    Sea sea;
    Vector3f offset_body = Vector3f::Zero();
    Matrix3f distortion = Matrix3f::Identity();
    float noise_uT = 0.0f;
    float duration_sec = 900.0f;
    float heading_rate_dps = 0.0f;
    unsigned seed = 7u;
};

void run(ContinuousMagHardIronEstimator& est, const Replay& r)
{
    std::mt19937 rng(r.seed);
    std::normal_distribution<float> white(0.0f, r.noise_uT);

    const int n = int(r.duration_sec / kDt);
    for (int i = 0; i < n; ++i) {
        const float t = float(i) * kDt;
        const float roll = r.sea.roll_amp_rad *
                           std::sin(2.0f * float(std::numbers::pi_v<float>) * r.sea.roll_hz * t);
        const float pitch = r.sea.pitch_amp_rad *
                            std::sin(2.0f * float(std::numbers::pi_v<float>) * r.sea.pitch_hz * t + 0.7f);

        const Quaternionf q_tilt = tiltQuat(roll, pitch);

        // The hull's own heading rotates the world field within the level
        // frame; the estimator is not told about it, which is the point.
        const float heading =
            r.heading_rate_dps * t * float(std::numbers::pi_v<float>) / 180.0f;
        const Vector3f field_level =
            Eigen::AngleAxisf(-heading, Vector3f::UnitZ()) * kFieldNed;

        Vector3f mag = q_tilt.conjugate() * field_level;
        mag = r.distortion * mag + r.offset_body;
        if (r.noise_uT > 0.0f) {
            mag += Vector3f(white(rng), white(rng), white(rng));
        }

        est.update(kDt, q_tilt, mag);
    }
}

Matrix3f misalignment(float rx, float ry, float rz)
{
    return (Eigen::AngleAxisf(rz, Vector3f::UnitZ()) *
            Eigen::AngleAxisf(ry, Vector3f::UnitY()) *
            Eigen::AngleAxisf(rx, Vector3f::UnitX())).toRotationMatrix();
}

// Heading a level-frame field vector implies, in degrees.
float headingDeg(const Vector3f& v)
{
    return std::atan2(v.y(), v.x()) * 180.0f / float(std::numbers::pi_v<float>);
}

void testCleanOffsetIsRecovered()
{
    ContinuousMagHardIronEstimator::Config cfg;
    cfg.model_ridge = 0.0f;
    cfg.model_ridge_relative = 0.0f;
    ContinuousMagHardIronEstimator est(cfg);

    Replay r;
    r.sea = Sea{0.20f, 0.26f};          // a big sea, and no sensor error but the offset
    r.offset_body = Vector3f(1.4f, -0.9f, 2.1f);
    run(est, r);

    check(est.valid(), "clean model: the solve is offered");
    check((est.estimate().bias_body_uT - r.offset_body).norm() < 0.05f,
          "clean model: offset recovered to better than 0.05 uT");
}

void testFixedAttitudeDeclines()
{
    ContinuousMagHardIronEstimator est;

    Replay r;
    r.sea = Sea{0.0f, 0.0f};
    r.offset_body = Vector3f(1.4f, -0.9f, 2.1f);
    run(est, r);

    check(!est.valid(), "fixed attitude: the solve declines");
    check(est.information() < 1.0e-3f, "fixed attitude: reported information is ~0");
}

// The property the ridge exists for.  A hull working through a large sea does
// not fit a cleaner offset than one working through a small sea when the error
// is aliased distortion rather than noise -- it fits the same wrong one with
// more confidence -- so the fraction the solve returns must not grow with the
// sea.  A fixed prior cannot hold that: the excitation outruns it, and the big
// sea hands back an offset the sensor never had.
void testShrinkageDoesNotChaseTheSeaState()
{
    const Vector3f offset(1.4f, -0.9f, 2.1f);
    const Matrix3f distortion =
        misalignment(0.012f, -0.009f, 0.006f) *
        (Matrix3f::Identity() + 0.012f * Matrix3f::Identity());
    const Sea seas[2] = {Sea{0.035f, 0.045f}, Sea{0.16f, 0.21f}};

    float ratio[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};   // [relative ridge?][sea]

    for (int rel = 0; rel < 2; ++rel) {
        for (int k = 0; k < 2; ++k) {
            ContinuousMagHardIronEstimator::Config cfg;
            if (rel == 0) cfg.model_ridge_relative = 0.0f;
            ContinuousMagHardIronEstimator est(cfg);

            Replay r;
            r.sea = seas[k];
            r.offset_body = offset;
            r.distortion = distortion;
            r.noise_uT = 0.8f;
            run(est, r);

            if (rel == 1) {
                check(est.valid(),
                      std::string("aliased distortion: solve is offered, sea ") +
                          (k ? "large" : "small"));
            }
            ratio[rel][k] = est.estimate().bias_body_uT.norm() / offset.norm();
        }
    }

    std::printf("    returned fraction, fixed ridge only:  small %.3f  large %.3f\n",
                double(ratio[0][0]), double(ratio[0][1]));
    std::printf("    returned fraction, with relative:     small %.3f  large %.3f\n",
                double(ratio[1][0]), double(ratio[1][1]));

    const float spread_fixed = ratio[0][1] - ratio[0][0];
    const float spread_rel   = ratio[1][1] - ratio[1][0];

    check(spread_fixed > 0.4f,
          "a fixed ridge alone lets the fraction chase the sea state");
    check(spread_rel < 0.75f * spread_fixed,
          "the relative ridge materially narrows that spread");
    check(ratio[1][0] < 1.0f && ratio[1][1] < 1.0f,
          "neither sea returns more offset than the sensor carries");
}

void testHeadingChangeIsRejected()
{
    ContinuousMagHardIronEstimator est;

    Replay r;
    r.sea = Sea{0.10f, 0.13f};
    r.offset_body = Vector3f(1.4f, -0.9f, 2.1f);
    r.heading_rate_dps = 0.05f;   // 45 deg over the replay
    run(est, r);

    check(!est.valid(), "hull turning: the residual gate declines the window");
    check(est.residualRmsUT() > est.config().max_residual_rms_uT,
          "hull turning: the residual is what declines it");
}

void testReferenceTracksTheAppliedOffset()
{
    ContinuousMagHardIronEstimator::Config cfg;
    cfg.model_ridge = 0.0f;
    cfg.model_ridge_relative = 0.0f;
    ContinuousMagHardIronEstimator est(cfg);

    const Vector3f offset(1.4f, -0.9f, 2.1f);
    Replay r;
    r.sea = Sea{0.20f, 0.26f};
    r.offset_body = offset;
    run(est, r);

    Vector3f ref_none, ref_full, ref_half;
    check(est.levelReferenceForBias(Vector3f::Zero(), ref_none) &&
              est.levelReferenceForBias(offset, ref_full) &&
              est.levelReferenceForBias(0.5f * offset, ref_half),
          "reference: available at zero, half and full correction");

    check((ref_full - kFieldNed).norm() < 0.1f,
          "reference: correcting fully recovers the true field");

    // Half the offset applied puts the reference half way there, which is what
    // keeps a slewing correction consistent with what it is subtracting.
    check(((ref_half - 0.5f * (ref_none + ref_full)).norm() < 0.02f),
          "reference: a partial correction lands half way");

    const float err_none = std::fabs(headingDeg(ref_none) - headingDeg(kFieldNed));
    const float err_full = std::fabs(headingDeg(ref_full) - headingDeg(kFieldNed));
    std::printf("    heading implied by the reference: %.3f deg uncorrected, "
                "%.3f deg corrected\n", double(err_none), double(err_full));

    check(err_none > 1.0f && err_full < 0.05f,
          "reference: the standing heading error is what the correction removes");
}

} // namespace

int main()
{
    testCleanOffsetIsRecovered();
    testFixedAttitudeDeclines();
    testShrinkageDoesNotChaseTheSeaState();
    testHeadingChangeIsRejected();
    testReferenceTracksTheAppliedOffset();

    if (g_failures) {
        std::printf("\n%d continuous hard-iron check(s) failed\n", g_failures);
        return 1;
    }
    std::printf("\nall continuous hard-iron checks passed\n");
    return 0;
}
