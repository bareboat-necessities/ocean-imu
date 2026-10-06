#!/bin/bash -e
set -e

./test_heave_hpdi
python3 ../../tools/heave_baselines.py selfcheck --hpdi-bin ./hpdi_replay
