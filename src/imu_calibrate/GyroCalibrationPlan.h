#pragma once
// Copyright 2026, Mikhail Grushinskiy
namespace imu_cal {
struct GyroCalibrationPlan {
  static constexpr int INITIAL=12,EXTRA=4,MAXIMUM=20;
  static constexpr int MIN_TRAIN=8,MIN_VALIDATION=4,MIN_ALIGNMENT_TRAIN=9;
  // Fixed before acquisition. Validation rotations never train the fit.
  static constexpr bool validation(int i) {return i<INITIAL?i%3==2:(i-INITIAL)%EXTRA==3;}
  static constexpr int validations(int n) {int count=0;for(int i=0;i<n;++i)count+=validation(i);return count;}
};
static_assert(GyroCalibrationPlan::validations(GyroCalibrationPlan::INITIAL)==GyroCalibrationPlan::MIN_VALIDATION);
static_assert(GyroCalibrationPlan::INITIAL-GyroCalibrationPlan::MIN_VALIDATION==GyroCalibrationPlan::MIN_TRAIN);
}
