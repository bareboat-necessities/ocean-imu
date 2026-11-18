EIGEN_DIR ?= ../../third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
CXXFLAGS = -O2 -std=c++20 -fno-finite-math-only -DEIGEN_UNROLLING_LIMIT=0 -I../../src $(EIGEN_CPPFLAGS) -Wall -Wextra
.PHONY: all build clean
all: startup-learning startup-test
build: all
startup-learning: startup-learning.o StartupLearner.o ShippingBaseline.o common.o
	$(CXX) $(CXXFLAGS) $^ -o $@
startup-test: startup-test.o StartupLearner.o
	$(CXX) $(CXXFLAGS) $^ -o $@
# The shipping wrapper instantiates a large Eigen navigation filter even
# though this driver stops at handover. Keep its host build memory bounded.
ShippingBaseline.o: ShippingBaseline.cpp
	$(CXX) $(filter-out -O2,$(CXXFLAGS)) -O1 -MMD -MP -c $< -o $@
common.o: ../../src/util/W3dSimCommon.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include *.d
clean:
	rm -f *.o *.d startup-learning startup-test
