#!/bin/bash -e

python3 spectrum-plots.py
python3 spectrum-estimates-plots.py --mode combined
python3 spectrum-sea_metrics.py \
  --csv sea_metrics_from_spectrum_report_classic.csv \
  --full-metrics-name-value sea_metrics_from_spectrum_full_metrics_classic.txt \
  --out-fragment ../../doc/spectrum/spectrum-sea_metrics-classic.tex-part \
  --fragment-input spectrum-sea_metrics-classic.tex-part \
  --title "SeaMetricsFromSpectrum Report (Classic/Goertzel Estimator)" \
  --caption "SeaMetricsFromSpectrum validation from classic (Goertzel) estimator spectra." \
  --label "tab:sea_metrics_from_spectrum_classic"
python3 spectrum-sea_metrics.py \
  --csv sea_metrics_from_spectrum_report_wavelets.csv \
  --full-metrics-name-value sea_metrics_from_spectrum_full_metrics_wavelets.txt \
  --out-fragment ../../doc/spectrum/spectrum-sea_metrics-wavelets.tex-part \
  --fragment-input spectrum-sea_metrics-wavelets.tex-part \
  --title "SeaMetricsFromSpectrum Report (Wavelets Estimator)" \
  --caption "SeaMetricsFromSpectrum validation from wavelets estimator spectra." \
  --label "tab:sea_metrics_from_spectrum_wavelets"
