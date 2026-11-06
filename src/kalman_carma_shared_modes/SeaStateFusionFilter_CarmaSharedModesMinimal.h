#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    SeaStateFusion_CarmaSharedModesMinimal<R> -- the EXPERIMENTAL 27-state
    shared-modes QMEKF (KalmanCarmaSharedModesMinimal<float, R>) behind the
    unchanged CARMA3 schedule, exactly as SeaStateFusion_CarmaSharedModes
    (SeaStateFusionFilter_CarmaSharedModes.h) puts the 33-state core there:

        SeaStateFusion_CarmaSharedModesMinimal<> f;   // Cascade realization
        f.begin(cfg);                                 // estimated coupling + confidence, Evidence excitation
        f.update(dt, gyr, acc, temp_c);
        f.updateMag(mag);
        f.raw().mekf().coupling();

    Not a default; see docs/kalman-carma-shared-modes-minimal.md.
*/

#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModesMinimal.h"

namespace ocean_imu::carma_shared_modes {

template<MinimalRealization R = MinimalRealization::Cascade>
class SeaStateFusion_CarmaSharedModesMinimal
    : public carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, KalmanCarmaSharedModesMinimal<float, R>> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Base = carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, KalmanCarmaSharedModesMinimal<float, R>>;
    using Config = typename Base::Config;
    using Mekf = KalmanCarmaSharedModesMinimal<float, R>;
    using SharedModes = typename Mekf::Config;

    // The 33-state wrapper's defaults: Evidence excitation, the deployable
    // shared-modes configuration.
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
