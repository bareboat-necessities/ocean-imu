.PHONY: all build test clean fetch-sim-data ensure-sim-data run-tests quality-gates

REPO_ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

TEST_DIRS := \
	$(REPO_ROOT)/tests/ahrs \
	$(REPO_ROOT)/tests/detrend \
	$(REPO_ROOT)/tests/freq \
	$(REPO_ROOT)/tests/imu_calibrate \
	$(REPO_ROOT)/tests/kalman_ou_ii \
	$(REPO_ROOT)/tests/kalman_ou_iii \
	$(REPO_ROOT)/tests/kalman_tfg \
	$(REPO_ROOT)/tests/nlo \
	$(REPO_ROOT)/tests/pii_observer \
	$(REPO_ROOT)/tests/spectrum \
	$(REPO_ROOT)/tests/spike_filter \
	$(REPO_ROOT)/tests/validation \
	$(REPO_ROOT)/tests/wave_dir \
	$(REPO_ROOT)/tests/wave_sim

SIM_DATA_VERSION ?= v1.2.1
SIM_DATA_REPO ?= bareboat-necessities/oceanography-waves-lib
SIM_DATA_ZIP_NAME ?= sim-data-files-vessel-rao-28ft.zip
SIM_DATA_ZIP ?= $(REPO_ROOT)/$(SIM_DATA_ZIP_NAME)
SIM_DATA_URL ?= https://github.com/$(SIM_DATA_REPO)/releases/download/$(SIM_DATA_VERSION)/$(SIM_DATA_ZIP_NAME)

all: build test

build:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d build; \
	done

test: ensure-sim-data
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
# The shared immutable cache avoids duplicating the 823 MB archive payload.
fetch-sim-data ensure-sim-data:
	python3 "$(REPO_ROOT)/tools/sim_dataset.py" --archive "$(SIM_DATA_ZIP)" \
		--dest $(TEST_DIRS) $(patsubst $(REPO_ROOT)/tests/%,$(REPO_ROOT)/plots/%,$(TEST_DIRS))

run-tests:
	@set -e; \
	for d in $(TEST_DIRS); do \
		if [ -f "$$d/run_tests.sh" ]; then \
			echo "Running $$d/run_tests.sh"; \
			( cd "$$d" && bash -e ./run_tests.sh ); \
		fi; \
	done
