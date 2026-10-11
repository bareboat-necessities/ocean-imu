#!/bin/bash -e
set -e

./rapid_motion-test

./mag_continuous-test
./pii_observer-adaptive
./local_gravity-test
