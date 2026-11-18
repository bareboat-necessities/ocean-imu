.PHONY: all build test clean fetch-sim-data ensure-sim-data run-tests quality-gates

REPO_ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

TEST_DIRS := \
	$(REPO_ROOT)/tests/ahrs \
	$(REPO_ROOT)/tests/detrend \
	$(REPO_ROOT)/tests/freq \
	$(REPO_ROOT)/tests/heave_only \
	$(REPO_ROOT)/tests/imu_calibrate \
	$(REPO_ROOT)/tests/kalman_carma \
	$(REPO_ROOT)/tests/kalman_carma_shared_modes \
	$(REPO_ROOT)/tests/kalman_ou_ii \
	$(REPO_ROOT)/tests/kalman_ou_iii \
	$(REPO_ROOT)/tests/kalman_tfg \
	$(REPO_ROOT)/tests/kalman_tfg_carma \
	$(REPO_ROOT)/tests/kalman_tfg_shared_modes \
	$(REPO_ROOT)/tests/kalman_tangent_shared_modes \
	$(REPO_ROOT)/tests/live_bias_consistency \
	$(REPO_ROOT)/tests/nlo \
	$(REPO_ROOT)/tests/pii_observer \
	$(REPO_ROOT)/tests/spectrum \
	$(REPO_ROOT)/tests/spike_filter \
	$(REPO_ROOT)/tests/startup_learning \
	$(REPO_ROOT)/tests/util \
	$(REPO_ROOT)/tests/validation \
	$(REPO_ROOT)/tests/wave_dir \
	$(REPO_ROOT)/tests/wave_sim

# One active direction for ordinary unit/study runs. The full default target
# additionally runs the complete, isolated two-direction regression matrix.
export OCEAN_IMU_SIM_CASE ?= following
SIM_JOBS ?= 3

all: build
	@status=0; \
	$(MAKE) -C "$(REPO_ROOT)" test || status=$$?; \
	$(MAKE) -C "$(REPO_ROOT)" sim-regression || status=$$?; \
	exit $$status

build:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d build; \
	done

test: build ensure-sim-data
	@$(MAKE) -C "$(REPO_ROOT)" run-tests

# Sanitizer, static-analysis and coverage layer. Run one gate with
# `tools/quality_gates.sh <gate>`; see docs/quality-gates.md.
quality-gates:
	@"$(REPO_ROOT)/tools/quality_gates.sh" all

clean:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d clean >/dev/null 2>&1 || true; \
	done

# Both targets verify archive bytes and refresh every stale input directory.
# The shared immutable cache avoids duplicating the unpacked input payload.
fetch-sim-data ensure-sim-data:
	python3 "$(REPO_ROOT)/tools/sim_dataset.py" \
		--dest $(TEST_DIRS) $(patsubst $(REPO_ROOT)/tests/%,$(REPO_ROOT)/plots/%,$(TEST_DIRS))

run-tests:
	@set -e; \
	for d in $(TEST_DIRS); do \
		if [ -f "$$d/run_tests.sh" ]; then \
			echo "Running $$d/run_tests.sh"; \
			( cd "$$d" && bash -e ./run_tests.sh ); \
		fi; \
	done

.PHONY: sim-regression
sim-regression: build
	python3 "$(REPO_ROOT)/tools/sim_matrix.py" --jobs $(SIM_JOBS)
