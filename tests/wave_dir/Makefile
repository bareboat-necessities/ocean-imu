.PHONY: all build test clean

CXX ?= g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)

BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Werror -pedantic -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

TARGET = wave-direction-test

all: build test
build: $(TARGET)
test: $(TARGET)
	./run_tests.sh

$(TARGET): wave-direction-test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET) $(TARGET).exe *.o
