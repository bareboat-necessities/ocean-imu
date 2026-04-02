#!/bin/bash -e

python3 spectrum-plots.py
python3 spectrum-estimates-plots.py --mode classic
python3 spectrum-estimates-plots.py --mode wavelets
python3 spectrum-sea_metrics.py --csv sea_metrics_from_spectrum_report.csv --out-fragment ../../doc/spectrum/sea_metrics_from_spectrum_table_fragment.tex-part
