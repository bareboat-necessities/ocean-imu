#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Experimental only. Reuses the CARMA schedule, startup, Evidence excitation,
// magnetic acquisition, adaptation and coupling confidence of the QMEKF.
#include "KalmanTfgCarmaSharedModes.h"
#include "SeaStateFusionFilter_CarmaSharedModes.h"
namespace ocean_imu::carma_shared_modes {
using SeaStateFusion_TfgCarmaSharedModes = SharedModesFusion<KalmanTfgCarmaSharedModes<float>>;
}  // namespace ocean_imu::carma_shared_modes
