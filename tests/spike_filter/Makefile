.PHONY: all build test clean

CXX ?= g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)

BASEFLAGS = -O3 -std=c++20 -Wall -Wextra -Werror -pedantic -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

TARGETS = spike-filter-test

all: build test
build: $(TARGETS)
test: $(TARGETS)
	bash ./run_tests.sh

spike-filter-test: spike-filter-test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS) $(addsuffix .exe,$(TARGETS)) *.o
