EIGEN_DIR ?= ../../third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
CXXFLAGS = -O3 -std=c++20 -fno-finite-math-only -DEIGEN_NON_ARDUINO -I../../src $(EIGEN_CPPFLAGS) -Wall -Wextra
all: live-bias-followup audit-sensors timing-test schmidt-test
build: all
live-bias-followup: main.o qmekf.o tfg.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@
common.o: ../../src/util/W3dSimCommon.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include *.d
.PHONY: all build clean
clean:
	rm -f *.o *.d live-bias-followup live-bias-baseline audit-sensors timing-test schmidt-test

audit-sensors: audit_sensors.o
	$(CXX) $(CXXFLAGS) $^ -o $@

timing-test: timing-test.cpp
	$(CXX) $(filter-out -O3,$(CXXFLAGS)) -O1 -DEIGEN_UNROLLING_LIMIT=0 -MMD -MP -MF $@.d -MT $@ $< -o $@

# Build exactly the same harness against an untouched baseline worktree.
# BASELINE_SRC is that worktree's absolute src directory; no production file
# is copied or modified. Both executables consume the same paired histories.
baseline-qmekf.o: qmekf.cpp Adapter.h Experiment.h
	$(CXX) -DBIAS_UNCHANGED_BASELINE -I$(BASELINE_SRC) $(CXXFLAGS) -MMD -MP -c $< -o $@
baseline-tfg.o: tfg.cpp Adapter.h Experiment.h
	$(CXX) -DBIAS_UNCHANGED_BASELINE -I$(BASELINE_SRC) $(CXXFLAGS) -MMD -MP -c $< -o $@
live-bias-baseline: main.o baseline-qmekf.o baseline-tfg.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@

schmidt-test: schmidt-test.cpp
	$(CXX) $(filter-out -O3,$(CXXFLAGS)) -O1 -DEIGEN_UNROLLING_LIMIT=0 -MMD -MP -MF $@.d -MT $@ $< -o $@
