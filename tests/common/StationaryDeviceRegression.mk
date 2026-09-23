# Add only the device regression target to a family's existing native rules.
# Keeping this separate leaves all scientific replay build recipes unchanged.
# Invoke from that family's run_tests.sh; dependency files come from its %.o rule.
stationary_device-test: stationary_device-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
