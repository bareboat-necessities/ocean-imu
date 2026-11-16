.PHONY: all build test ensure-sim-data

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

UNAME_S := $(shell uname -s)
include $(REPO_ROOT)/tests/common/HostReplayFlags.mk

LDFLAGS  =

# Tell make to look for source files in src/util
VPATH = $(REPO_ROOT)/src/util

TARGETS = pii_observer-adaptive local_gravity-test phase_matching-test

all: build test

build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

pii_observer-adaptive: pii_observer-adaptive.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

phase_matching-test: phase_matching-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

local_gravity-test: local_gravity-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp Makefile $(REPO_ROOT)/tests/common/HostReplayFlags.mk
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(wildcard *.d)

clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv
