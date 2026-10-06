.PHONY: all build test clean

CC = g++
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)

# Same code-generating flags as the simulator Makefiles, so each heave_export
# build replays its estimator bit-for-bit like the shipped simulator does.
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif
WARN = -Wall -Wextra -Wshadow -Wconversion
HPDI_WARN = $(WARN) -Wdouble-promotion

EXPORTS = heave_export-ou3 heave_export-ou2 heave_export-tfg heave_export-pii heave_export-nlo
TARGETS = test_heave_hpdi hpdi_replay $(EXPORTS)

all: build test

build: $(TARGETS)

# Unit test, bilinear self-check and scorer self-check; no simulation data.
test: test_heave_hpdi hpdi_replay
	bash ./run_tests.sh

test_heave_hpdi: test_heave_hpdi.cpp $(REPO_ROOT)/src/hpdi/HeaveHPDI.h
	$(CC) $(CXXFLAGS) $(HPDI_WARN) -o $@ $<

hpdi_replay: $(REPO_ROOT)/tools/hpdi_replay.cpp $(REPO_ROOT)/src/hpdi/HeaveHPDI.h
	$(CC) $(CXXFLAGS) $(HPDI_WARN) -o $@ $<

W3dSimCommon.o: $(REPO_ROOT)/src/util/W3dSimCommon.cpp
	$(CC) $(CXXFLAGS) -MMD -MP -c $< -o $@

FAMILY_ou3 = OU_III
FAMILY_ou2 = OU_II
FAMILY_tfg = TFG
FAMILY_pii = PII
FAMILY_nlo = TVG_NLO

# Link-time hooks on the shared runner; see heave_export.cpp.
EXPORT_WRAPS = \
	-Wl,--wrap=_ZN19W3dSimulationRunnerC1E20W3dSimulationOptions21SimulationNoiseModelsR17IW3dFusionAdapter \
	-Wl,--wrap=_Z24print_validation_metricsRK22W3dSimulationRunResultffPKc

# The shipped unit compiles the simulator source exactly as its own Makefile
# does. Its anonymous-namespace helpers trigger GCC's subobject-linkage warning
# only because the source is included rather than compiled directly.
$(EXPORTS:%=%-shipped.o): heave_export-%-shipped.o: heave_export_shipped.cpp
	$(CC) $(CXXFLAGS) -DHEAVE_EXPORT_$(FAMILY_$*) -Wno-subobject-linkage -MMD -MP -c $< -o $@

$(EXPORTS:=.o): heave_export-%.o: heave_export.cpp
	$(CC) $(CXXFLAGS) $(WARN) -DHEAVE_EXPORT_$(FAMILY_$*) -MMD -MP -c $< -o $@

# The shipped unit links first, as the simulator object does in its own build.
$(EXPORTS): heave_export-%: heave_export-%-shipped.o heave_export-%.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(EXPORT_WRAPS)

-include $(wildcard *.d)

clean:
	rm -f *.o *.d $(TARGETS)
