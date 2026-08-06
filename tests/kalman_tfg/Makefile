.PHONY: all build test ensure-sim-data clean

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
# Lie++ is the independent oracle for the world block of the group. Test-only:
# it never reaches the Arduino build. See third_party/Lie-plusplus/VENDORED.md.
LIEPP_CPPFLAGS := -isystem $(REPO_ROOT)/third_party/Lie-plusplus/include
BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(LIEPP_CPPFLAGS) $(CPPFLAGS)

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

TARGETS = lie_group-test convention-test

all: build test

build: $(TARGETS)

test: $(TARGETS)
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

lie_group-test: lie_group-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

convention-test: convention-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# -MMD -MP records the headers each object depends on, so editing a group or
# filter header rebuilds the test instead of silently reusing a stale object.
%.o: %.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(wildcard *.d)

clean:
	rm -f *.o *.d $(TARGETS) *.exe *.csv
