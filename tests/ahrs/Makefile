CC = g++
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I./../../src -I/usr/include/eigen3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)
endif

LDFLAGS  =

# Tell make to look for source files in ../../src
VPATH = ../../src/util

TARGETS = ahrs-qmekf-sim

all: $(TARGETS)

ahrs-qmekf-sim: ahrs-qmekf-sim.o W3dSimCommon.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) *.exe *.csv