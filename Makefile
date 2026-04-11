.PHONY: all clean fetch-sim-data

TEST_DIRS := \
	tests/ahrs \
	tests/detrend \
	tests/freq \
	tests/imu_calibrate \
	tests/kalman_ou_ii \
	tests/kalman_ou_iii \
	tests/pii_observer \
	tests/spectrum \
	tests/wave_sim

SIM_DATA_VERSION ?= v1.1.1
SIM_DATA_REPO ?= bareboat-necessities/oceanography-waves-lib
SIM_DATA_ZIP ?= sim-data-files.zip
SIM_DATA_URL ?= https://github.com/$(SIM_DATA_REPO)/releases/download/$(SIM_DATA_VERSION)/$(SIM_DATA_ZIP)

all:
	@set -e; \
	for d in $(TEST_DIRS); do \
		$(MAKE) -C $$d all; \
	done

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
		plot_dir="plots/$${d#tests/}"; \
		test_dir="$$d"; \
		mkdir -p "$$plot_dir"; \
		unzip -o "$(SIM_DATA_ZIP)" -d "$$plot_dir" >/dev/null; \
		echo "Unpacked $(SIM_DATA_ZIP) -> $$plot_dir"; \
		unzip -o "$(SIM_DATA_ZIP)" -d "$$test_dir" >/dev/null; \
		echo "Unpacked $(SIM_DATA_ZIP) -> $$test_dir"; \
	done
