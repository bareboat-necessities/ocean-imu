// Bit-level characterization trace of SeaStateFusionFilter_TFG; see
// ou_wrapper_trace.cpp for the method.
#define EIGEN_NON_ARDUINO
#include <cstdio>
#include <string>
#include <Eigen/Dense>

#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
#include "SyntheticMarineHistory.h"

extern const float g_std = 9.80665f;

using namespace characterization;
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;

static void put_tfg(Trace& tr, const Fusion& f) {
    tr.put(static_cast<int>(f.stage())); tr.put(f.isLive());
    tr.put(f.magReferenceLearned()); tr.put(f.magReferenceRefined());
    tr.put(f.magRefineTimeSec()); tr.put(f.magNorthLockTimeSec());
    tr.put(f.magHardIronBodyUT()); tr.put(f.magContinuousHardIronAppliedUT());
    tr.put(f.magContinuousHardIron().estimate().valid);
    tr.put(f.isTunerReady()); tr.put(f.pseudoUpdatePeriodSec());
    tr.put(f.handoffTimedOut()); tr.put(f.getRSFilterInput());
    tr.put(f.accelVibrationGuardEngagement()); tr.put(f.accelVibrationRms());
    tr.put(f.accelVibrationRaccStd());
    tr.put(f.getTauApplied()); tr.put(f.getSigmaApplied()); tr.put(f.getRSApplied());
    tr.put(f.getTauTarget()); tr.put(f.getSigmaTarget()); tr.put(f.getRSTarget());
    tr.put(f.getWavePeriodSec()); tr.put(f.wavePeriodUsable()); tr.put(f.wavePeriodReady());
    tr.put(f.getAccelVariance()); tr.put(f.getSigmaVarianceHorizonSec());
    tr.put(f.quaternion());
    Eigen::Quaternionf qt = Eigen::Quaternionf::Identity();
    tr.put(f.startupTiltQuaternion(qt)); tr.put(qt);
    tr.put(f.get_velocity()); tr.put(f.get_position()); tr.put(f.get_world_accel());
    const auto& m = f.mekf();
    tr.put(m.get_integral_displacement());
    tr.put(m.gyroscope_bias()); tr.put(m.get_acc_bias());
    tr.put(m.acc_bias_updates_enabled());
    tr.put(m.lastMagDiag().accepted); tr.put(m.lastMagDiag().r);
    tr.put(m.covariance_full().diagonal());
}

static void run(const std::string& path, Fusion::Config cfg, const HistoryConfig& hc,
                int knobs, Trace& all) {
    Fusion f;
    f.begin(cfg);
    Trace tr;
    tr.fp = std::fopen(path.c_str(), "wb");
    SyntheticMarineHistory h(hc);
    Sample s;
    int k = 0;
    while (h.next(s)) {
        ++k;
        f.update(0.005f, s.gyro, s.acc, 30.0f);
        if (s.mag_due) f.updateMag(s.mag);
        if (knobs) {
            if (k == 120 * 200) f.setRSLaw(Fusion::RSLaw::LegacyCubic);
            if (k == 130 * 200) f.setTauScaledPseudoCadence(false);
            if (k == 140 * 200) f.setTauScaledPseudoCadence(true);
            if (k == 150 * 200) f.setRSLaw(Fusion::RSLaw::SpectralMSE);
            if (k == 155 * 200) f.setAccelVibrationRaccGain(1.2f);
            if (k == 160 * 200) f.setChannelFreeze(true, 3.0f, 0.4f, false, 0.0f);
            if (k == 170 * 200) f.setAdaptationTimeConstants(2.0f);
            if (k == 180 * 200) f.setFixedTuning(3.0f, 0.4f, 6.0f);
            if (k == 190 * 200) f.setPeriodicAwCovSync(false);
        }
        put_tfg(tr, f);
    }
    std::fclose(tr.fp);
    std::printf("%s values=%llu hash=%016llx\n", path.c_str(),
                (unsigned long long)tr.values, (unsigned long long)tr.hash);
    all.hash ^= tr.hash;
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    Trace all;
    HistoryConfig full;
    run(out + "/tfg_proxy_mag.bin", Fusion::Config{}, full, 0, all);

    Fusion::Config nomag;
    nomag.with_mag = false;
    HistoryConfig shortc = full;
    shortc.duration_sec = 260.0;
    run(out + "/tfg_proxy_nomag.bin", nomag, shortc, 0, all);

    Fusion::Config staged;
    staged.startup_init_policy = Fusion::StartupInitPolicy::StagedMekf;
    run(out + "/tfg_staged.bin", staged, shortc, 0, all);

    Fusion::Config hi;
    hi.mag_estimate_hard_iron = true;
    HistoryConfig capsize = full;
    capsize.capsize = true;
    capsize.seed = 0xc0ffeeu;
    run(out + "/tfg_proxy_hardiron_knobs.bin", hi, capsize, 1, all);

    std::printf("tfg combined=%016llx\n", (unsigned long long)all.hash);
    return 0;
}
