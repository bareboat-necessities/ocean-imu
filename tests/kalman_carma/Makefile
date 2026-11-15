.PHONY: all build test clean
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
LDFLAGS =
# Published lever-arm ratios include amplified negative controls. Keep their
# arithmetic independent of the runner's SIMD/FMA instruction set; the other
# host targets retain their native performance profile.
LEVER_REPRO_FLAGS = $(BASEFLAGS) -ffp-contract=off -fno-tree-vectorize -DEIGEN_DONT_VECTORIZE
VPATH = $(REPO_ROOT)/src/util
TARGETS = carma_band_normalization-test carma_model-test carma_process-sim carma3-test tfg_carma3-test carma3_bench carma3-sim \
          carma_lever_arm-test carma_lever_arm-sim carma_excitation-test carma_coupled-test
all: build test
build: $(TARGETS)
test: $(TARGETS)
	bash ./run_tests.sh
carma_band_normalization-test: carma_band_normalization-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_model-test: carma_model-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_process-sim: carma_process-sim.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3-test: carma3-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
tfg_carma3-test: tfg_carma3-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3_bench: carma3_bench.o carma3_filters_ou.o carma3_filters_c4.o carma3_filters_c5.o carma3_filters_tfg.o carma3_filters_coupled.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3-sim: carma3-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_lever_arm-test: carma_lever_arm-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_lever_arm-sim: carma_lever_arm-repro.o W3dSimCommon-repro.o
	$(CC) $(LEVER_REPRO_FLAGS) -o $@ $^ $(LDFLAGS)
carma_lever_arm-repro.o: carma_lever_arm-sim.cpp Makefile
	$(CC) $(LEVER_REPRO_FLAGS) -MMD -MP -c $< -o $@
W3dSimCommon-repro.o: $(REPO_ROOT)/src/util/W3dSimCommon.cpp Makefile
	$(CC) $(LEVER_REPRO_FLAGS) -MMD -MP -c $< -o $@
carma_excitation-test: carma_excitation-test.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_coupled-test: carma_coupled-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(wildcard *.d)
clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv *.f64
