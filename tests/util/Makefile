.PHONY: all build test clean

CXX ?= g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)

BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Werror -pedantic -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

HEADERS = \
	$(REPO_ROOT)/src/util/AngleUtils.h \
	$(REPO_ROOT)/src/util/ImuLoopHelpers.h \
	$(REPO_ROOT)/src/util/MagneticHeading.h \
	$(REPO_ROOT)/src/util/QuaternionUtils.h \
	$(REPO_ROOT)/src/wave_dir/WaveDirectionReport.h \
	$(REPO_ROOT)/src/detrend/AdaptiveWaveDetrender.h

TARGETS = ins-helpers-test

all: build test
build: $(TARGETS)
test: $(TARGETS)
	bash ./run_tests.sh

ins-helpers-test: ins-helpers-test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS) $(addsuffix .exe,$(TARGETS)) *.o
