.PHONY: all build test ensure-sim-data

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

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

TARGETS = pii_observer-adaptive local_gravity-test

all: build test

TARGETS += mag_continuous-test
mag_continuous-test: mag_continuous-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
mag_continuous-test.o: ../common/MagneticContinuousReplay.h ../common/MagneticMotion.h ../../src/util/MagneticInput.h ../../src/util/MagneticRotation.h

TARGETS += rapid_motion-test
rapid_motion-test: rapid_motion-test.o AttitudeProtection.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
rapid_motion-test.o: ../common/RapidMotionReplay.h ../common/RapidMotionBackends.h ../../src/util/AttitudeProtection.h ../../src/util/ImuFreshness.h

build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

pii_observer-adaptive: pii_observer-adaptive.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

local_gravity-test: local_gravity-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv

TARGETS += ins_heave-test
ins_heave-test: ins_heave-test.o AttitudeProtection.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
ins_heave-test.o: ins_heave-test.cpp ../common/InsHeaveReplay.h
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
AttitudeProtection.o: $(REPO_ROOT)/src/util/AttitudeProtection.h
build: ins_heave-test

-include $(wildcard *.d)
