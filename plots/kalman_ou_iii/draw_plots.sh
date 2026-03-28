#!/bin/bash -e

python3 ./kalman_ou_iii-plots.py
python3 ../spectrum/spectrum-plots.py  # .tex needs some spectrum .pgf later
