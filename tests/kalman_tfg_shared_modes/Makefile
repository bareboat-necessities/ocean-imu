.PHONY: all build test clean ensure-sim-data
CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)
UNAME_S := $(shell uname -s)
include $(REPO_ROOT)/tests/common/HostReplayFlags.mk
LDFLAGS =
VPATH = $(REPO_ROOT)/src/util
TARGETS = kalman_tfg_shared_modes-sim
all: build test
build: $(TARGETS)
test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh
ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data
kalman_tfg_shared_modes-sim: kalman_tfg_shared_modes-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.cpp Makefile $(REPO_ROOT)/tests/common/HostReplayFlags.mk
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(wildcard *.d)
clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv
