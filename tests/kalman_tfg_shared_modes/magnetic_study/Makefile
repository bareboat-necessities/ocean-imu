# Compile diagnostic TUs against an isolated instrumentation of baseline source.
# The ordinary simulators, evidence, firmware, and their build flags are unaffected.
ROOT := $(abspath ../../..)
OUT ?= build
OVERLAY := $(OUT)/overlay
EIGEN_DIR ?= $(ROOT)/third_party/eigen
EIGEN := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),$(EIGEN_DIR),/usr/include/eigen3)
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -funroll-loops -fno-finite-math-only
ARCH := $(if $(filter x86_64%,$(shell $(CXX) -dumpmachine)),-march=x86-64-v3,)
FLAGS = $(CXXFLAGS) $(ARCH) -I$(OVERLAY) -I$(OVERLAY)/src -I. -I$(ROOT)/src -isystem $(EIGEN)
INPUTS := $(shell sed -n 's|^diff --git .* b/|$(ROOT)/|p' instrumentation.patch)
.PHONY: all check
all: $(OUT)/magnetic-study $(OUT)/magnetic-reference-study $(OUT)/information-test
check: $(OUT)/information-test $(OUT)/packet-test
	$(OUT)/information-test
	$(OUT)/packet-test
$(OUT)/packet-test: $(OUT)/packet-test.o
	$(CXX) $(FLAGS) $^ -o $@
$(OVERLAY)/ready: prepare.py instrumentation.patch $(INPUTS)
	python3 prepare.py $(OVERLAY)
	touch $@
$(OUT)/%.o: %.cpp.in $(OVERLAY)/ready
	$(CXX) $(FLAGS) -MMD -MP -x c++ -c $< -o $@
$(OUT)/%.o: %.cpp $(OVERLAY)/ready
	$(CXX) $(FLAGS) -MMD -MP -c $< -o $@
$(OUT)/W3dSimCommon.o: $(OVERLAY)/ready
	$(CXX) $(FLAGS) -MMD -MP -c $(OVERLAY)/src/util/W3dSimCommon.cpp -o $@
$(OUT)/magnetic-study: $(OUT)/replay.o $(OUT)/SelectiveMag.o $(OUT)/W3dSimCommon.o
	$(CXX) $(FLAGS) $^ -o $@
$(OUT)/magnetic-reference-study: $(OUT)/reference.o $(OUT)/W3dSimCommon.o
	$(CXX) $(FLAGS) $^ -o $@
$(OUT)/information-test: $(OUT)/information-test.o $(OUT)/SelectiveMag.o
	$(CXX) $(FLAGS) $^ -o $@
-include $(wildcard $(OUT)/*.d)
