#!/bin/bash -e

./spectrum-estimator-sim
./spectrum-wavelets-estimator-sim

./sea-metrics-from-spectrum-test
python3 ../../plots/spectrum/sea-metrics-report.py --csv sea_metrics_from_spectrum_report.csv --out-fragment ../../doc/spectrum/sea_metrics_from_spectrum_table_fragment.tex-part --out-main ../../doc/spectrum/sea_metrics_from_spectrum_table.tex
