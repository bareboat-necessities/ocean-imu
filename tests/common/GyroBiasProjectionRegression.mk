gyro_bias_projection-test: gyro_bias_projection-test.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
