.PHONY: all clean

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
