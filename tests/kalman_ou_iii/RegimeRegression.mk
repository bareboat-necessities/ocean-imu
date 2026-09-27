# Additional native regression; the simulator build contract is unchanged.
regime_ambiguity-test: regime_ambiguity-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include regime_ambiguity-test.d
