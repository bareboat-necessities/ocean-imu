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
TARGETS = detrend-basic-test detrend-wave-test detrend-wave3d-test

all: build test

build: $(TARGETS)

test: ensure-sim-data $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

detrend-basic-test: detrend-basic-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	find . -maxdepth 1 -type f -name '*.csv' ! -name 'detrend*.csv' -delete
	rm -f *.o $(TARGETS) *.exe
