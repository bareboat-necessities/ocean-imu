.PHONY: all build test clean
REPO_ROOT := $(abspath ../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
CXXFLAGS := -O3 -std=c++20 -fno-finite-math-only -DEIGEN_UNROLLING_LIMIT=0 -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) -Wall -Wextra
all: build test
build: geometry-test tangent-study tangent-replay
common.o: $(REPO_ROOT)/src/util/W3dSimCommon.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
tangent-study: study.o q.o independent.o tfg.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@
geometry-test: geometry-test.cpp
	$(CXX) $(filter-out -O3,$(CXXFLAGS)) -O1 -MMD -MP -MF $@.d -MT $@ $< -o $@
tangent-replay: tangent-replay.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@
# Same experiment source against an independently checked-out main baseline.
baseline-%.o: %.cpp
	$(CXX) -DTANGENT_BASELINE -I$(BASELINE_SRC) $(CXXFLAGS) -MMD -MP -c $< -o $@
tangent-baseline: baseline-study.o baseline-q.o baseline-independent.o baseline-tfg.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@
test: build
	bash ./run_tests.sh
clean:
	rm -f *.o *.d geometry-test tangent-study tangent-replay tangent-baseline
-include $(wildcard *.d)
