#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Experimental 27-state reduced tangent-bias geometry. All propagation,
// CARMA adaptation, physical covariance writes, Joseph updates and startup
// are the current shared TFG implementation. Only the finite correction
// retraction and its COMPLETE covariance reset change. See the article for
// the distinction from full TG-EqF and its virtual position-input bias.
#include "KalmanTfgCarmaSharedModes.h"
#include "SeaStateFusionFilter_CarmaSharedModes.h"
namespace ocean_imu::carma_shared_modes {
template<typename T=float>
class KalmanTangentCarmaSharedModes : public KalmanTfgCarmaSharedModes<T> {
public:
    using Base=KalmanTfgCarmaSharedModes<T>;
    template<class... Args>
        requires (!(sizeof...(Args)==1 && (std::is_base_of_v<KalmanTangentCarmaSharedModes,std::remove_cvref_t<Args>> && ...)))
    explicit KalmanTangentCarmaSharedModes(Args&&... args):Base(std::forward<Args>(args)...) {
        this->set_tangent_bias_geometry(true);
    }
};
using SeaStateFusion_TangentCarmaSharedModes=SharedModesFusion<KalmanTangentCarmaSharedModes<float>>;
}
