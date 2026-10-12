# Host run of the FIFO timing diagnostic sketch against a simulated BMI270
# (FifoTimingSketch.cpp), compiled unchanged against ins_sketch_stubs/.
# Invoke from a family's run_tests.sh, like InsSketchHeaveDrift.mk.
FIFO_SKETCH_FILE := $(REPO_ROOT)/sensors/imu_basic/atomS3R_imu_fifo_timing/atomS3R_imu_fifo_timing.ino
FIFO_SKETCH_FLAGS := -O2 -std=c++20 -fno-finite-math-only \
	-I$(REPO_ROOT)/tests/common/ins_sketch_stubs -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) \
	'-DFIFO_SKETCH_FILE="$(FIFO_SKETCH_FILE)"'

fifo_timing_sketch-test.o: $(REPO_ROOT)/tests/common/FifoTimingSketch.cpp $(FIFO_SKETCH_FILE) \
		$(wildcard $(REPO_ROOT)/tests/common/ins_sketch_stubs/*.h $(REPO_ROOT)/tests/common/ins_sketch_stubs/AtomS3R/*.h)
	$(CXX) $(FIFO_SKETCH_FLAGS) -MMD -MP -MF $@.d -MT $@ -c -o $@ $(REPO_ROOT)/tests/common/FifoTimingSketch.cpp

fifo_timing_sketch-test: fifo_timing_sketch-test.o
	$(CXX) $(FIFO_SKETCH_FLAGS) -o $@ $< $(LDFLAGS)

-include fifo_timing_sketch-test.o.d
