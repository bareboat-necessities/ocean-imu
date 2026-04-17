.PHONY: all build test ensure-sim-data
TEST_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TEST_DIR)/../..)

all: build test

build:
	@true

test: ensure-sim-data
	bash ./run_tests.sh

ensure-sim-data:
	@$(MAKE) -C "$(REPO_ROOT)" ensure-sim-data

clean:
	find . -maxdepth 1 -type f -name '*.csv' ! -name 'spectrum*.csv' -delete
	rm -f *.o $(TARGETS) *.exe
