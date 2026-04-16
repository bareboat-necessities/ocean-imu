CC = g++
CFLAGS = -O3 -std=c++20 -funroll-loops -fno-finite-math-only -I./../../src $(CPPFLAGS)
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
