.PHONY: all build test ensure-sim-data clean

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

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

TARGETS = heave_baselines-sim heave_baselines-test

all: build test

build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

heave_baselines-sim: heave_baselines-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

heave_baselines-test: heave_baselines-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) *.exe *.csv
