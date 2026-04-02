#!/bin/bash -e

python3 spectrum-plots.py
python3 spectrum-estimates-plots.py --mode classic
python3 spectrum-estimates-plots.py --mode wavelets
python3 spectrum-sea_metrics.py \
  --csv sea_metrics_from_spectrum_report.csv \
  --full-metrics-csv sea_metrics_from_spectrum_full_metrics.csv \
  --out-fragment ../../doc/spectrum/spectrum-sea_metrics.tex-part \
  --out-main ../../doc/spectrum/spectrum-sea_metrics.tex
