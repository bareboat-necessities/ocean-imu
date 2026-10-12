#pragma once
#include <cstdint>
namespace ocean_imu {
class ImuFreshness {
public:
    /** @brief Forget acceleration availability and any sample pending gyro delivery. */
    void reset() { have_accel_ = pending_accel_ = false; }
    /** @brief Record sensor availability and return whether a gyro frame can be delivered.
     * Set accel_fresh only for acceleration not consumed by an earlier delivered
     * frame; acceleration-only observations remain pending until a gyro arrives.
     */
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
