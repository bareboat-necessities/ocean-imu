#!/bin/bash -e
set -e
./heave_only-test
./heave_only-sim
./heave_only-sim --frontend proxy
./heave_only-sim --frontend truth
