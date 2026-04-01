CC = g++
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I./../../src

ifneq (,$(wildcard ../../third_party/eigen/Eigen/Dense))
    EIGEN_INC = -I./../../third_party/eigen
else
    EIGEN_INC = -I/usr/include/eigen3
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) $(EIGEN_INC) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) $(EIGEN_INC) -march=native
else
    CXXFLAGS = $(BASEFLAGS) $(EIGEN_INC)
endif

TARGETS = spectrum-estimator-sim spectrum-wavelets-estimator-sim sea-metrics-from-spectrum-test

all: $(TARGETS)

spectrum-estimator-sim: spectrum-estimator-sim.o
	$(CC) $(CXXFLAGS) -o $@ $^

spectrum-wavelets-estimator-sim: spectrum-wavelets-estimator-sim.o
	$(CC) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) *.exe spectrum_*estimate*.csv

sea-metrics-from-spectrum-test: sea-metrics-from-spectrum-test.o
	$(CC) $(CXXFLAGS) -o $@ $^
