#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Mahony update with a heading-only magnetometer correction.

  Mahony_AHRS::updateMag() adds the magnetometer error to the accelerometer
  error on all three axes with the same gain. The magnetometer carries little
  tilt information, so its calibration residuals (hard/soft iron,
  misalignment) and local field disturbances tilt the attitude estimate.

  Here the magnetometer error is computed exactly as in updateMag(), but only
  its component about the estimated vertical is kept, so it corrects heading
  and never roll or pitch. Tilt comes from the accelerometer alone, as in the
  IMU-only update.

  Mahony feeds every correction in as an extra body rate, so the heading term
  is added to the gyro rate (proportional part) and to the integral
  accumulator (integral part) before the unchanged IMU-only update runs. The
  result equals a single Mahony update whose error is the accelerometer error
  plus the projected magnetometer error.
*/

#include <cmath>

#include "ahrs/Mahony_AHRS.h"

template<typename T>
inline void mahony_AHRS_update_heading_mag(Mahony_AHRS<T>* m,
                                           T gx, T gy, T gz,
                                           T ax, T ay, T az,
                                           T mx, T my, T mz,
                                           T* pitch, T* roll, T* yaw,
                                           T delta_t_sec)
{
    const bool acc_ok = !((ax == T(0)) && (ay == T(0)) && (az == T(0)));
    const bool mag_ok = !((mx == T(0)) && (my == T(0)) && (mz == T(0)));
    if (acc_ok && mag_ok) {
        const T recipNorm = Mahony_AHRS<T>::invSqrt(mx * mx + my * my + mz * mz);
        mx *= recipNorm;
        my *= recipNorm;
        mz *= recipNorm;

        const T q0 = m->q0, q1 = m->q1, q2 = m->q2, q3 = m->q3;
        const T q0q0 = q0 * q0, q0q1 = q0 * q1, q0q2 = q0 * q2, q0q3 = q0 * q3;
        const T q1q1 = q1 * q1, q1q2 = q1 * q2, q1q3 = q1 * q3;
        const T q2q2 = q2 * q2, q2q3 = q2 * q3, q3q3 = q3 * q3;

        // Reference direction of Earth's magnetic field (as in updateMag).
        const T hx = T(2) * (mx * (T(0.5) - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
        const T hy = T(2) * (mx * (q1q2 + q0q3) + my * (T(0.5) - q1q1 - q3q3) + mz * (q2q3 - q0q1));
        const T hz = T(2) * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (T(0.5) - q1q1 - q2q2));
        const T bx = std::sqrt(hx * hx + hy * hy);
        const T bz = hz;

        // Estimated directions of gravity and magnetic field.
        const T halfvx = q1q3 - q0q2;
        const T halfvy = q0q1 + q2q3;
        const T halfvz = T(0.5) * (q0q0 - q1q1 - q2q2 + q3q3);
        const T halfwx = bx * (T(0.5) - q2q2 - q3q3) + bz * (q1q3 - q0q2);
        const T halfwy = bx * (q1q2 - q0q3)           + bz * (q0q1 + q2q3);
        const T halfwz = bx * (q0q2 + q1q3)           + bz * (T(0.5) - q1q1 - q2q2);

        // Magnetometer error, projected onto the estimated vertical.
        const T emx = my * halfwz - mz * halfwy;
        const T emy = mz * halfwx - mx * halfwz;
        const T emz = mx * halfwy - my * halfwx;
        const T vn2 = halfvx * halfvx + halfvy * halfvy + halfvz * halfvz;
        if (vn2 > T(0)) {
            const T k = (emx * halfvx + emy * halfvy + emz * halfvz) / vn2;
            const T ehx = k * halfvx, ehy = k * halfvy, ehz = k * halfvz;

            if (m->twoKi > T(0)) {
                m->integralFBx += m->twoKi * ehx * delta_t_sec;
                m->integralFBy += m->twoKi * ehy * delta_t_sec;
                m->integralFBz += m->twoKi * ehz * delta_t_sec;
            }
            gx += m->twoKp * ehx;
            gy += m->twoKp * ehy;
            gz += m->twoKp * ehz;
        }
    }
    m->update(gx, gy, gz, ax, ay, az, pitch, roll, yaw, delta_t_sec);
}
