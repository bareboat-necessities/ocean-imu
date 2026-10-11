#pragma once
#include <cstdint>
namespace ocean_imu {
class ImuFreshness {
public:
    void reset() { have_accel_ = pending_accel_ = false; }
    bool observe(bool accel, bool gyro, bool& accel_fresh) {
        if (accel) have_accel_ = pending_accel_ = true;
        accel_fresh = false;
        if (!gyro || !have_accel_) return false;
        accel_fresh = pending_accel_;
        pending_accel_ = false;
        return true;
    }
private:
    bool have_accel_ = false, pending_accel_ = false;
};
} // namespace ocean_imu
