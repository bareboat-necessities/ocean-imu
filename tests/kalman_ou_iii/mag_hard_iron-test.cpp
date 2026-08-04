/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    MagAutoTuner joint hard-iron estimation.

    The estimate is opt-in, and these tests are mostly about why.  Separating a
    body-fixed offset from the world field needs the accumulation frame to
    rotate, and how well it can be done is decided entirely by what rotates and
    by what else is wrong with the sensor at the same time:

      - a platform working through tens of degrees, heading and attitude both,
        sweeps the offset over a solid body of directions and the fit is exact;
      - a pure change of heading is not enough on its own, because rotating
        about the vertical never modulates the part of the offset that points
        along it;
      - tilt alone, of the few degrees a hull takes through a wave, constrains
        the offset so weakly that the model's blind spots -- soft iron, and
        error in the frame it is handed -- compete with the signal, and both of
        those are modulated by the same attitude the offset is read from;
      - one fixed attitude makes offset and field literally indistinguishable,
        and the solve has to say so rather than return a number.

    The arithmetic itself is exact whenever the model holds, which is the point:
    what fails in a wave record is the model, not the solver.
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "tuner/MagAutoTuner.h"

namespace {

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::printf("%-62s %s\n", what.c_str(), ok ? "ok" : "FAILED");
    if (!ok) ++g_failures;
}

constexpr float DEG = 3.14159265358979323846f / 180.0f;

// Field at the simulator's default site: 52 uT, 66.5 deg dip, in a NED frame
// whose x axis is magnetic north.
Eigen::Vector3f referenceField()
{
    const float incl = 66.5f * DEG;
    return 52.0f * Eigen::Vector3f(std::cos(incl), 0.0f, std::sin(incl));
}

Eigen::Quaternionf attitude(float roll_deg, float pitch_deg, float yaw_deg)
{
    return Eigen::Quaternionf(
        Eigen::AngleAxisf(yaw_deg * DEG, Eigen::Vector3f::UnitZ()) *
        Eigen::AngleAxisf(pitch_deg * DEG, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(roll_deg * DEG, Eigen::Vector3f::UnitX()));
}

struct Sample {
    Eigen::Quaternionf q_bw;   // accumulation frame, BODY->WORLD
    Eigen::Vector3f mag_body;  // what the magnetometer reads
};

// One synthetic startup window.
//
//   frame_error_deg  tilts the accumulation frame away from the attitude the
//                    body actually held, which is what an imperfect estimator
//                    hands the tuner.
//   soft_iron        a small symmetric distortion, the term the model omits.
std::vector<Sample> makeWindow(const Eigen::Vector3f& hard_iron_body,
                               float roll_amp_deg,
                               float pitch_amp_deg,
                               float yaw_amp_deg,
                               int count,
                               float noise_uT,
                               float frame_error_deg,
                               float soft_iron)
{
    const Eigen::Vector3f field = referenceField();

    Eigen::Matrix3f distortion = Eigen::Matrix3f::Identity();
    distortion(0, 0) += soft_iron;
    distortion(1, 1) -= soft_iron;
    distortion(0, 1) = distortion(1, 0) = 0.5f * soft_iron;

    std::mt19937 rng(20260804u);
    std::normal_distribution<float> white(0.0f, noise_uT);

    std::vector<Sample> samples;
    samples.reserve(static_cast<size_t>(count));

    for (int i = 0; i < count; ++i) {
        // Two incommensurate rates, so roll, pitch and heading do not trace a
        // single closed figure that would leave the fit degenerate anyway.
        const float phase = 6.2831853f * static_cast<float>(i) / 137.0f;

        const float roll = roll_amp_deg * std::sin(phase);
        const float pitch = pitch_amp_deg * std::sin(0.7f * phase + 1.1f);
        const float yaw = yaw_amp_deg * std::sin(0.31f * phase);

        const Eigen::Quaternionf q_true = attitude(roll, pitch, yaw);

        Sample sample;
        sample.mag_body =
            distortion * (q_true.conjugate() * field) + hard_iron_body;
        sample.mag_body +=
            Eigen::Vector3f(white(rng), white(rng), white(rng));

        // The frame the tuner is given, wrong by a fixed tilt.
        sample.q_bw =
            q_true * Eigen::Quaternionf(
                Eigen::AngleAxisf(frame_error_deg * DEG,
                                  Eigen::Vector3f::UnitY()));
        sample.q_bw.normalize();

        samples.push_back(sample);
    }
    return samples;
}

MagAutoTuner::Config baseConfig(int count)
{
    MagAutoTuner::Config cfg;
    cfg.estimate_hard_iron = true;
    cfg.min_samples = count;
    cfg.min_window_sec = 0.0f;
    cfg.sample_dt_sec = 1.0f / 25.0f;
    return cfg;
}

Eigen::Vector3f runWindow(MagAutoTuner& tuner,
                          const std::vector<Sample>& samples,
                          bool& valid)
{
    for (const auto& sample : samples) {
        tuner.addSampleWithTiltQuatDt(1.0f / 25.0f,
                                      sample.q_bw,
                                      Eigen::Vector3f(0.0f, 0.0f, -9.80665f),
                                      Eigen::Vector3f::Zero(),
                                      sample.mag_body);
    }
    Eigen::Vector3f estimate = Eigen::Vector3f::Zero();
    valid = tuner.getHardIronBodyUT(estimate);
    return estimate;
}

// A platform that really moves: heading swings and the hull works through tens
// of degrees, so the offset traces a solid body of directions in the world
// frame while the field does not, and the two come apart cleanly.
//
// Heading change alone is not enough, which is worth knowing: rotating about
// the vertical never modulates the part of the offset that points along it, so
// a pure turn leaves one direction as unconstrained as standing still does.
void test_real_motion_recovers_the_offset()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 20.0f, 20.0f, 70.0f, 2000, 0.8f,
                                    0.0f, 0.0f);

    MagAutoTuner tuner;
    tuner.setConfig(baseConfig(2000));

    bool valid = false;
    const Eigen::Vector3f estimate = runWindow(tuner, samples, valid);

    check(valid, "real motion: offset is estimated");
    check((estimate - truth).norm() < 0.25f,
          "real motion: offset recovered to better than 0.25 uT");
}

// The same window with the estimate not asked for.
void test_off_by_default()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 20.0f, 20.0f, 70.0f, 2000, 0.8f,
                                    0.0f, 0.0f);

    MagAutoTuner::Config cfg = baseConfig(2000);
    cfg.estimate_hard_iron = false;   // the shipped default

    MagAutoTuner tuner;
    tuner.setConfig(cfg);

    bool valid = false;
    const Eigen::Vector3f estimate = runWindow(tuner, samples, valid);

    check(!valid, "default: no offset is produced");
    check(estimate.isZero(), "default: nothing to subtract");
    check(tuner.isReady(), "default: the reference is still learned");
}

// One attitude for the whole window.  The offset and a rotated field are the
// same reading, and the solve must decline rather than split them arbitrarily.
void test_fixed_attitude_declines()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 0.0f, 0.0f, 0.0f, 2000, 0.8f,
                                    0.0f, 0.0f);

    MagAutoTuner tuner;
    tuner.setConfig(baseConfig(2000));

    bool valid = false;
    runWindow(tuner, samples, valid);

    check(!valid, "fixed attitude: solve declines");
    check(tuner.hardIronInformation() < 1.0f,
          "fixed attitude: reported information is ~0");
    check(tuner.isReady(),
          "fixed attitude: the reference is still learned");
}

// Wave tilt only, a few degrees, no heading change: this is the sea-state case,
// and it is the one the default is off for.  The information floor is what
// keeps a window this weak from being used.
void test_wave_tilt_alone_is_refused_by_the_information_floor()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 4.0f, 2.0f, 0.0f, 2000, 0.8f,
                                    0.0f, 0.0f);

    MagAutoTuner tuner;
    tuner.setConfig(baseConfig(2000));

    bool valid = false;
    runWindow(tuner, samples, valid);

    check(!valid, "wave tilt alone: solve declines at the default floor");
    check(tuner.hardIronInformation() <
              tuner.config().min_hard_iron_information,
          "wave tilt alone: information is below the floor");
}

// The floor is a bound on noise, not on bias.  Dropping it lets a weak window
// through, and what comes back is wrong for a reason no amount of further
// averaging fixes: the omitted soft-iron term is modulated by the same attitude
// the offset is being read from, so it lands in the answer.  Note the window
// here is noise-free -- more samples would not help.
void test_the_information_floor_does_not_bound_bias()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 4.0f, 2.0f, 0.0f, 20000, 0.0f,
                                    0.0f, 0.015f);

    MagAutoTuner::Config cfg = baseConfig(20000);
    cfg.min_hard_iron_information = 1.0e-6f;   // accept anything

    MagAutoTuner tuner;
    tuner.setConfig(cfg);

    bool valid = false;
    const Eigen::Vector3f estimate = runWindow(tuner, samples, valid);

    check(valid, "weak window: accepted once the floor is removed");
    check((estimate - truth).norm() > 0.5f,
          "weak window: and the answer is wrong despite no noise at all");
}

// An error in the frame handed to the tuner is likewise attitude-correlated,
// so it too is absorbed into the offset rather than averaged away.
void test_frame_error_biases_the_offset()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);

    MagAutoTuner::Config cfg = baseConfig(2000);

    MagAutoTuner clean;
    clean.setConfig(cfg);
    bool clean_valid = false;
    const Eigen::Vector3f clean_estimate = runWindow(
        clean,
        makeWindow(truth, 20.0f, 20.0f, 70.0f, 2000, 0.0f, 0.0f, 0.0f),
        clean_valid);

    MagAutoTuner tilted;
    tilted.setConfig(cfg);
    bool tilted_valid = false;
    const Eigen::Vector3f tilted_estimate = runWindow(
        tilted,
        makeWindow(truth, 20.0f, 20.0f, 70.0f, 2000, 0.0f, 0.5f, 0.0f),
        tilted_valid);

    check(clean_valid && (clean_estimate - truth).norm() < 0.05f,
          "exact frame: offset recovered essentially exactly");
    check(tilted_valid &&
              (tilted_estimate - truth).norm() >
                  10.0f * (clean_estimate - truth).norm(),
          "half-degree frame error: offset is materially worse");
}

// Whatever the offset comes out as, it must not silently corrupt the reference:
// the gauge-fixed field the MEKF consumes still has to look like a field.
void test_reference_stays_sane_when_the_offset_is_used()
{
    const Eigen::Vector3f truth(1.4f, -0.9f, 0.6f);
    const auto samples = makeWindow(truth, 20.0f, 20.0f, 70.0f, 2000, 0.8f,
                                    0.0f, 0.0f);

    MagAutoTuner tuner;
    tuner.setConfig(baseConfig(2000));

    bool valid = false;
    runWindow(tuner, samples, valid);

    Eigen::Vector3f reference;
    check(tuner.getMagWorldRef(reference), "reference is produced");
    check(std::fabs(reference.y()) < 1.0e-4f,
          "reference is gauge-fixed onto +X");
    check(std::fabs(reference.norm() - referenceField().norm()) < 1.0f,
          "reference norm matches the true field to within 1 uT");
}

} // namespace

int main()
{
    test_real_motion_recovers_the_offset();
    test_off_by_default();
    test_fixed_attitude_declines();
    test_wave_tilt_alone_is_refused_by_the_information_floor();
    test_the_information_floor_does_not_bound_bias();
    test_frame_error_biases_the_offset();
    test_reference_stays_sane_when_the_offset_is_used();

    if (g_failures != 0) {
        std::printf("\n%d check(s) failed\n", g_failures);
        return 1;
    }
    std::printf("\nall mag hard-iron checks passed\n");
    return 0;
}
