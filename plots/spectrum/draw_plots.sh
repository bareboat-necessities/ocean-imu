#!/bin/bash -e

python3 spectrum-plots.py
python3 spectrum-estimates-plots.py --mode classic
python3 spectrum-estimates-plots.py --mode wavelets
