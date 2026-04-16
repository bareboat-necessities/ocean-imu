CC = g++
EIGEN_DIR ?= ../../third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-I$(EIGEN_DIR),-I/usr/include/eigen3)
CFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I./../../src $(EIGEN_CPPFLAGS) $(CPPFLAGS)
TARGETS = wave-convention-check

.PHONY: all
all: $(TARGETS)

wave-convention-check: wave-convention-check.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(TARGETS) *.exe
