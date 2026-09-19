.PHONY: all build h18-superword clean
CC = g++
TOOL_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
REPO_ROOT := $(abspath $(TOOL_DIR)/../..)
EIGEN_DIR ?= $(REPO_ROOT)/third_party/eigen
EIGEN_CPPFLAGS := $(if $(wildcard $(EIGEN_DIR)/Eigen/Dense),-isystem $(EIGEN_DIR),-isystem /usr/include/eigen3)
CXXFLAGS = -O2 -std=c++20 -Wall -Wextra -Wshadow -Wconversion -fno-finite-math-only \
           -I$(REPO_ROOT)/src $(EIGEN_CPPFLAGS) $(CPPFLAGS)
LDFLAGS =

# The export runs the shipping estimator, so it is built against the same
# sources the proof is pinned to. Optimization is left portable on purpose:
# the diagnostic reports the map it measured, and a machine-specific
# instruction selection is not part of the shipping contract.
EXPORT_SOURCE := $(TOOL_DIR)/ou3_theorem/h18_superword_export.cpp
EXPORT_BINARY := $(TOOL_DIR)/ou3_theorem/h18_superword_export

OUTPUT_DIR ?= /tmp
EXPORT_JSON ?= $(OUTPUT_DIR)/ou3-h18-superword.json
REPORT_JSON ?= $(OUTPUT_DIR)/ou3-h18-superword-feasibility.json

all: build

build: $(EXPORT_BINARY)

$(EXPORT_BINARY): $(EXPORT_SOURCE)
	$(CC) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# Measure one magnetically informed H18 service superword and evaluate the
# candidate storage on it at high precision. Neither step is a certificate.
h18-superword: $(EXPORT_BINARY)
	$(EXPORT_BINARY) --output $(EXPORT_JSON)
	python3 $(TOOL_DIR)/ou3_theorem/h18_superword.py \
	    --export $(EXPORT_JSON) --output $(REPORT_JSON)

clean:
	rm -f $(EXPORT_BINARY)
