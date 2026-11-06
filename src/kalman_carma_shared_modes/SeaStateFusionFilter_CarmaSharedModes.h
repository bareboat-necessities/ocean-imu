#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    SeaStateFusion_CarmaSharedModes -- the EXPERIMENTAL shared-modes QMEKF
    (KalmanCarmaSharedModes<float>) behind the unchanged CARMA3 schedule
    (kalman_carma/SeaStateFusionFilter_CARMA3.h: proxy startup, magnetic
    acquisition and north lock, OU-III's tuner and envelope, 0.1 s model
    commits, tilt watchdog, wave-excitation modes, outputs).  It is
    SeaStateFusion_CARMA3 with the shared-modes core as its inner filter; the
    only addition is begin(cfg, shared_modes_cfg), which configures the core
    after the wrapper has created it.  Use it exactly as SeaStateFusion_CARMA3:

        SeaStateFusion_CarmaSharedModes f;
        f.begin(cfg);                       // estimated coupling + confidence, Evidence excitation
        f.update(dt, gyr, acc, temp_c);     // every IMU sample
        f.updateMag(mag);                   // at the magnetometer rate
        f.raw().mekf().coupling();          // confidence, state, applied coupling
        f.raw().mekf().waveDirection();     // direction of travel and uncertainty

    Not a default; see docs/kalman-carma-shared-modes.md.
*/

#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"

namespace ocean_imu::carma_shared_modes {

class SeaStateFusion_CarmaSharedModes
    : public carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, KalmanCarmaSharedModes<float>> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Base = carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, KalmanCarmaSharedModes<float>>;
    using Config = typename Base::Config;
    using Mekf = KalmanCarmaSharedModes<float>;
    using SharedModes = Mekf::Config;

    // The wave excitation defaults to Evidence mode (kalman_carma/CarmaWaveEvidence.h,
    // docs/carma-wave-excitation.md), the mode the AtomS3R sketches deploy:
    // the Legacy floor and stillness attenuation make the amplitude jump on
    // small seas, which the complexified realization tolerates worse than
    // the CARMA3 baseline (docs/kalman-carma-shared-modes.md section 2).
    void begin(const Config& cfg, const SharedModes& shared_modes = SharedModes{},
               carma::CarmaExcitationMode excitation = carma::CarmaExcitationMode::Evidence) {
        Base::begin(cfg);
        this->raw().setCarmaExcitationMode(excitation);
        this->raw().mekf().configureSharedModes(shared_modes);
    }

    Mekf& sharedModes() { return this->raw().mekf(); }
    const Mekf& sharedModes() const { return this->raw().mekf(); }
};

}  // namespace ocean_imu::carma_shared_modes
