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
VPATH = $(REPO_ROOT)/src/util
TARGETS = carma_model-test carma_process-sim carma3-test tfg_carma3-test carma3_bench carma3-sim \
          carma_lever_arm-test carma_lever_arm-sim
all: build test
build: $(TARGETS)
test: $(TARGETS)
	bash ./run_tests.sh
carma_model-test: carma_model-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_process-sim: carma_process-sim.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3-test: carma3-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
tfg_carma3-test: tfg_carma3-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3_bench: carma3_bench.o carma3_filters_ou.o carma3_filters_c4.o carma3_filters_c5.o carma3_filters_tfg.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma3-sim: carma3-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_lever_arm-test: carma_lever_arm-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
carma_lever_arm-sim: carma_lever_arm-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(wildcard *.d)
clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv *.f64
