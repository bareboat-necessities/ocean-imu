CC = g++

# Base flags
BASEFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I./../src -I/usr/include/eigen3

# Detect platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    CXXFLAGS = $(BASEFLAGS) -march=native
else ifeq ($(UNAME_S),Darwin)   # macOS
    CXXFLAGS = $(BASEFLAGS) -march=native
else
    CXXFLAGS = $(BASEFLAGS)    # Windows / Cygwin / others
endif

LDFLAGS  =

# Targets
TARGETS = freq_track

all: $(TARGETS)

freq_track: freq_track.o
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) *.exe *.csv