.PHONY: all build test clean ensure-sim-data
CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
BENCH_DIR := $(REPO_ROOT)/tests/kalman_carma
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src -I$(BENCH_DIR) $(EIGEN_CPPFLAGS) $(CPPFLAGS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif
LDFLAGS =
VPATH = $(REPO_ROOT)/src/util
TARGETS = carma_shared_modes-test kalman_carma_shared_modes-sim
# carma3_bench's engine and filter units (tests/kalman_carma), compiled here
# so the Monte Carlo mode pairs the new filter with the same truth and
# sensors; the engine's main() is renamed carma3_bench_main.
BENCH_OBJS = bench_engine.o bench_carma3_filters_ou.o bench_carma3_filters_c4.o bench_carma3_filters_c5.o \
             bench_carma3_filters_tfg.o bench_carma3_filters_coupled.o
all: build test
build: $(TARGETS)
test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh
ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data
carma_shared_modes-test: carma_shared_modes-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
kalman_carma_shared_modes-sim: kalman_carma_shared_modes-sim.o csm_filters.o $(BENCH_OBJS) W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
bench_engine.o: $(BENCH_DIR)/carma3_bench.cpp
	$(CC) $(CXXFLAGS) -Dmain=carma3_bench_main -MMD -MP -c $< -o $@
bench_%.o: $(BENCH_DIR)/%.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
%.o: %.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(wildcard *.d)
clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv
