.PHONY: all build test ensure-sim-data

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

LDFLAGS  =

# Tell make to look for source files in src/util
VPATH = $(REPO_ROOT)/src/util

TARGETS = kalman_ou_iii-sim ou3-certificate-sim ou3-information-certificate-sim kalman_ou_common-test aw_covariance_policy-test acc_bias_ou-test channel_freeze-test wave_period-test mag_hard_iron-test continuous_mag_hard_iron-test tuner_coupling-test tuner_schedule-test wave_band_sigma-test iss_contract-test rs_law-test startup_init-test

all: build test

build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

kalman_ou_iii-sim: kalman_ou_iii-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

ou3-certificate-sim: ou3-certificate-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

ou3-information-certificate-sim: ou3-information-certificate-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

kalman_ou_common-test: kalman_ou_common-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

aw_covariance_policy-test: aw_covariance_policy-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

acc_bias_ou-test: acc_bias_ou-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

channel_freeze-test: channel_freeze-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

wave_period-test: wave_period-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Needs W3dSimCommon for g_std, like the simulator target.
tuner_coupling-test: tuner_coupling-test.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mag_hard_iron-test: mag_hard_iron-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

continuous_mag_hard_iron-test: continuous_mag_hard_iron-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Needs W3dSimCommon for g_std, like the simulator target.
startup_init-test: startup_init-test.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

tuner_schedule-test: tuner_schedule-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

wave_band_sigma-test: wave_band_sigma-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

iss_contract-test: iss_contract-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

rs_law-test: rs_law-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# -MMD -MP records the headers each object depends on, so editing a filter
# header rebuilds the simulator instead of silently reusing a stale object.
%.o: %.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(wildcard *.d)

clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv
