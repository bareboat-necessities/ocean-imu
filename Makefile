.PHONY: all build test clean fetch-sim-data ensure-sim-data run-tests

REPO_ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

TEST_DIRS := \
	$(REPO_ROOT)/tests/ahrs \
	$(REPO_ROOT)/tests/detrend \
	$(REPO_ROOT)/tests/freq \
	$(REPO_ROOT)/tests/imu_calibrate \
	$(REPO_ROOT)/tests/kalman_ou_ii \
	$(REPO_ROOT)/tests/kalman_ou_iii \
	$(REPO_ROOT)/tests/pii_observer \
	$(REPO_ROOT)/tests/spectrum \
	$(REPO_ROOT)/tests/wave_sim

SIM_DATA_VERSION ?= v1.1.3
SIM_DATA_REPO ?= bareboat-necessities/oceanography-waves-lib
SIM_DATA_ZIP_NAME ?= sim-data-files.zip
SIM_DATA_ZIP ?= $(REPO_ROOT)/$(SIM_DATA_ZIP_NAME)
SIM_DATA_URL ?= https://github.com/$(SIM_DATA_REPO)/releases/download/$(SIM_DATA_VERSION)/$(SIM_DATA_ZIP_NAME)
SIM_DATA_CHECK_FILE ?= $(REPO_ROOT)/tests/detrend/wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv

all: build test

build:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d build; \
	done

test: ensure-sim-data
	@$(MAKE) -C "$(REPO_ROOT)" run-tests

clean:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d clean >/dev/null 2>&1 || true; \
	done

fetch-sim-data:
	@set -e; \
	echo "Downloading $(SIM_DATA_URL)"; \
	curl -fL "$(SIM_DATA_URL)" -o "$(SIM_DATA_ZIP)"; \
	for d in $(TEST_DIRS); do \
		plot_dir="$(REPO_ROOT)/plots/$${d#$(REPO_ROOT)/tests/}"; \
		test_dir="$$d"; \
		mkdir -p "$$plot_dir"; \
		unzip -o "$(SIM_DATA_ZIP)" -d "$$plot_dir" >/dev/null; \
		echo "Unpacked $(SIM_DATA_ZIP) -> $$plot_dir"; \
		unzip -o "$(SIM_DATA_ZIP)" -d "$$test_dir" >/dev/null; \
		echo "Unpacked $(SIM_DATA_ZIP) -> $$test_dir"; \
	done

ensure-sim-data:
	@if [ ! -f "$(SIM_DATA_CHECK_FILE)" ]; then \
		echo "Wave simulation data not found ($(SIM_DATA_CHECK_FILE)); fetching..."; \
		$(MAKE) -C "$(REPO_ROOT)" fetch-sim-data; \
	else \
		echo "Wave simulation data found ($(SIM_DATA_CHECK_FILE))"; \
	fi

run-tests:
	@set -e; \
	for d in $(TEST_DIRS); do \
		if [ -f "$$d/run_tests.sh" ]; then \
			echo "Running $$d/run_tests.sh"; \
			( cd "$$d" && bash ./run_tests.sh ); \
		fi; \
	done
