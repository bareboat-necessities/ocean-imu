# Local-gravity regression target, added to a family's existing native rules
# the same way as StationaryDeviceRegression.mk (invoked from run_tests.sh).
local_gravity-test: local_gravity-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
