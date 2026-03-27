#!/bin/bash -e

python3 wave_sim-plots.py
python3 wave_sim-theories.py --formats png svg pgf --output-dir . --basename wave_sim-theories
