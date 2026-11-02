.PHONY: all build test ensure-sim-data

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)

# Base flags
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

# Detect platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)   # macOS
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)    # Windows / Cygwin / others
endif

LDFLAGS  =

# Targets
TARGETS = calibration_frames-test sensor_noise-replay gyro_rotation-test sensor_noise-test calibration_accuracy-test mag_field_consistency-test mag_hand_motion-test calibration_workflow-test calibration_safety-test imu_calibrate-test accel_cal-test accel_cal-replay accel_cal_gravity_override-test mag_diagnostics-test

all: build test

TARGETS += mag_rotation-test mag_timing_cal-test
mag_timing_cal-test: mag_timing_cal-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
mag_timing_cal-test.o: ../../src/imu_calibrate/MagTimingCal.h ../../src/util/MagneticRotation.h ../../src/AtomS3R/AtomS3R_ImuCalBlob.h ../../src/AtomS3R/AtomS3R_MagAcquisition.h
mag_rotation-test: mag_rotation-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
mag_rotation-test.o: ../../src/util/MagneticRotation.h ../../src/AtomS3R/AtomS3R_MagAcquisition.h

MAG_FUSION_TARGETS = $(addprefix magnetic_fusion_,1 2 3 4 5 6)
TARGETS += magnetic_innovation-test $(MAG_FUSION_TARGETS)
build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

calibration_accuracy-test: calibration_accuracy-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mag_field_consistency-test: mag_field_consistency-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mag_hand_motion-test: mag_hand_motion-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

calibration_safety-test: calibration_safety-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

calibration_workflow-test: calibration_workflow-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

imu_calibrate-test: imu_calibrate-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

accel_cal-test: accel_cal-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

accel_cal-replay: accel_cal-replay.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

accel_cal_gravity_override-test: accel_cal_gravity_override-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mag_diagnostics-test: mag_diagnostics-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# The library headers under test (header-only), so edits there rebuild the tests.
SRC_DIR = ../../src
CAL_HEADERS = $(SRC_DIR)/imu_calibrate/GyroCalibrationPlan.h $(SRC_DIR)/imu_calibrate/CalibrateIMU.h $(SRC_DIR)/imu_calibrate/MagCalFit.h \
              $(SRC_DIR)/imu_calibrate/MagCalCapture.h $(SRC_DIR)/imu_calibrate/MagCalSampling.h $(SRC_DIR)/imu_calibrate/GyroCalCapture.h $(SRC_DIR)/imu_calibrate/AccelCalFit.h \
              $(SRC_DIR)/imu_calibrate/AccelCalCapture.h $(SRC_DIR)/AtomS3R/AtomS3R_ImuCalBlob.h \
              $(SRC_DIR)/AtomS3R/AtomS3R_ImuUnits.h $(SRC_DIR)/imu_calibrate/SensorNoise.h

calibration_accuracy-test.o: $(CAL_HEADERS)
mag_field_consistency-test.o: $(CAL_HEADERS)
mag_hand_motion-test.o: $(CAL_HEADERS)
calibration_safety-test.o: $(SRC_DIR)/imu_calibrate/MagCalCapture.h $(CAL_HEADERS)
calibration_workflow-test.o: $(CAL_HEADERS) $(SRC_DIR)/AtomS3R/AtomS3R_CalLog.h
imu_calibrate-test.o: $(CAL_HEADERS)
accel_cal-test.o: accel_cal_sim.h accel_cal_replay.h ../common/GravityChain.h $(CAL_HEADERS)
accel_cal-replay.o: accel_cal_replay.h $(CAL_HEADERS)
accel_cal_gravity_override-test.o: ../common/GravityChain.h $(CAL_HEADERS)
mag_diagnostics-test.o: $(SRC_DIR)/AtomS3R/AtomS3R_MagDiagnostics.h $(SRC_DIR)/AtomS3R/AtomS3R_MagAcquisition.h $(SRC_DIR)/AtomS3R/AtomS3R_Bmm150AuxPreset.h $(SRC_DIR)/AtomS3R/AtomS3R_Bmm150Compensation.h

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	find . -maxdepth 1 -type f -name '*.csv' ! -name 'calibrate*.csv' -delete
	rm -f calibrate_accel_replay_sample.log
	rm -f *.o $(TARGETS) *.exe

NEW_CAL_HEADERS = ../../src/imu_calibrate/GyroCalFit.h ../../src/imu_calibrate/SensorNoise.h
gyro_rotation-test: gyro_rotation-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
sensor_noise-test: sensor_noise-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
gyro_rotation-test.o sensor_noise-test.o: $(NEW_CAL_HEADERS) $(CAL_HEADERS)

magnetic_innovation-test: magnetic_innovation-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
magnetic_innovation-test.o: ../../src/kalman_common/MagneticInnovationGate.h

MAG_HEADERS = $(wildcard ../../src/kalman_common/*.h ../../src/kalman_carma/*.h ../../src/kalman_ou_iii/*.h ../../src/kalman_ou_ii/*.h ../../src/kalman_tfg/*.h ../../src/ahrs/*.h)
magnetic_fusion_%: magnetic_fusion-test.cpp $(MAG_HEADERS)
	$(CC) $(BASEFLAGS) -O1 -DEIGEN_UNROLLING_LIMIT=0 -DMAG_FILTER_KIND=$* $< -o $@

calibration_frames-test: calibration_frames-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
sensor_noise-replay: sensor_noise-replay.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
calibration_frames-test.o: $(CAL_HEADERS) ../../src/imu_calibrate/FrameAlignment.h ../../src/AtomS3R/AtomS3R_MagDiagnostics.h
sensor_noise-replay.o: $(NEW_CAL_HEADERS) $(CAL_HEADERS)
