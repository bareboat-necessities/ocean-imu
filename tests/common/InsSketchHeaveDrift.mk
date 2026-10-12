# Heave/displacement drift regression for this family's shipped INS sketch.
# The unchanged .ino is compiled against host stand-ins for Arduino/M5 in
# ins_sketch_stubs/ and driven by InsSketchHeaveDrift.cpp.
# Invoke from a family's run_tests.sh, like the other shared regressions.
INS_SKETCH_ROOT := $(REPO_ROOT)/sensors/full_marine_ins
INS_FAMILY := $(notdir $(CURDIR))

INS_SKETCH_kalman_ou_ii  := atomS3R_ins_kalman_ou2
INS_SKETCH_kalman_ou_iii := atomS3R_ins_kalman_ou3
INS_SKETCH_kalman_tfg    := atomS3R_ins_tfg
INS_SKETCH_nlo           := atomS3R_ins_nlo
INS_SKETCH_pii_observer  := atomS3R_ins_pii_observer

INS_READY_kalman_ou_ii  := (app).fusion_.isLive()
INS_READY_kalman_ou_iii := (app).fusion_.isLive()
INS_READY_kalman_tfg    := (app).fusion_.isLive()
INS_READY_nlo           := (app).fusion_.initialized()
INS_READY_pii_observer  := (app).mahony_seeded_

INS_DISPLACEMENT_kalman_ou_ii  := (app).fusion_.displacementUpMeters()
INS_DISPLACEMENT_kalman_ou_iii := (app).fusion_.displacementUpMeters()
INS_DISPLACEMENT_kalman_tfg    := (app).displacement_up_m_
INS_DISPLACEMENT_nlo           := (app).fusion_.snapshot().disp_zu

INS_SKETCH_NAME := $(INS_SKETCH_$(INS_FAMILY))
INS_SKETCH_FILE := $(INS_SKETCH_ROOT)/$(INS_SKETCH_NAME)/$(INS_SKETCH_NAME).ino
INS_DISPLACEMENT := $(INS_DISPLACEMENT_$(INS_FAMILY))
INS_SKETCH_FLAGS := -O2 -std=c++20 -fno-finite-math-only -DSEA_STATE_ENABLE_WIZARD=0 \
	-I$(REPO_ROOT)/tests/common/ins_sketch_stubs -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) \
	'-DINS_SKETCH_FILE="$(INS_SKETCH_FILE)"' '-DINS_READY(app)=$(INS_READY_$(INS_FAMILY))' \
	$(if $(INS_DISPLACEMENT),'-DINS_DISPLACEMENT_UP(app)=$(INS_DISPLACEMENT)')

# The sketch includes the real estimator headers; the driver object records
# them so an edited header rebuilds the harness instead of a stale binary.
ins_sketch_heave-test.o: $(REPO_ROOT)/tests/common/InsSketchHeaveDrift.cpp $(INS_SKETCH_FILE) \
		$(wildcard $(REPO_ROOT)/tests/common/ins_sketch_stubs/*.h $(REPO_ROOT)/tests/common/ins_sketch_stubs/AtomS3R/*.h)
	$(if $(INS_SKETCH_NAME),,$(error no INS sketch is mapped to $(INS_FAMILY)))
	$(CXX) $(INS_SKETCH_FLAGS) -MMD -MP -MF $@.d -MT $@ -c -o $@ $(REPO_ROOT)/tests/common/InsSketchHeaveDrift.cpp

ins_sketch_heave-test: ins_sketch_heave-test.o
	$(CXX) $(INS_SKETCH_FLAGS) -o $@ $< $(REPO_ROOT)/src/util/W3dSimCommon.cpp $(LDFLAGS)

-include ins_sketch_heave-test.o.d
