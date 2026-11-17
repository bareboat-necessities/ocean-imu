#define EIGEN_NON_ARDUINO
#include "../live_bias_consistency/StandardHeaders.h"
#include "ShippingBaseline.h"
// Read the acquired reference without adding a diagnostic API to production.
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
namespace startup_learning {
struct ShippingBaseline::Impl {
    ocean_imu::ou3::SeaStateFusion_OU_III<TrackerType::KALMANF> fusion;
    Handover result;
    double previous=0;
    bool have_previous=false;
    Impl() {
        ocean_imu::ou3::SeaStateFusionConfig_OU_III c;
        c.turn_separator.forward_body=Eigen::Vector3f::UnitY();
        fusion.begin(c);
    }
};
ShippingBaseline::ShippingBaseline():impl_(std::make_unique<Impl>()) {}
ShippingBaseline::~ShippingBaseline()=default;
bool ShippingBaseline::step(const Sample& s) {
    auto& p=*impl_;
    if (p.result.valid) return false;
    const double dt=p.have_previous?s.time-p.previous:0.005;
    p.previous=s.time; p.have_previous=true;
    p.fusion.update(float(dt),s.gyro.cast<float>(),s.accel.cast<float>());
    // Capture BEFORE a live magnetic correction; this is the seed itself.
    if (p.fusion.isLive()) {
        auto& h=p.result; const auto& core=p.fusion.raw().mekf();
        h.time=s.time; h.attitude=p.fusion.attitudeQuat().cast<double>();
        h.field=core.v2ref.cast<double>();
        h.down_body=h.attitude.conjugate()*V::UnitZ();
        h.gyro_bias=core.gyroscope_bias().cast<double>();
        h.gyro_bias_cov=core.covariance_full().template block<3,3>(3,3).cast<double>();
        h.accel_bias=core.get_acc_bias().cast<double>();
        h.angular_rate=s.gyro-h.gyro_bias; h.last_raw_gyro=s.gyro;
        h.valid=true; h.ready=true;
        return true;
    }
    if (s.magnetic) p.fusion.updateMag(s.mag.cast<float>());
    return false;
}
Handover ShippingBaseline::handover() const {return impl_->result;}
}
