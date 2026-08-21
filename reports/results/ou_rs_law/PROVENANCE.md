# OU-III regularizer-law ablation figure provenance

`ou_rs_roundtrip_transition.svg` and its decimated series
`ou_rs_roundtrip_transition.csv` are one round-trip realization of the
ablation's non-stationary instrument, replayed through the deployed OU-III
filter and plotted by `tools/ou_validation.py:plot_transition_diagnostic`, the
same drawing code the general-validation transition figure uses.

Regenerate with:

    python3 tools/ou_rs_law_ablation.py --diagnostic-only \
        --output-dir reports/results/ou_rs_law

and mirror the SVG to `doc/kalman_ou_iii/ou_rs_roundtrip_transition.svg` for
the manuscript.  The plotted realization is seed triplet 11/101/1009, the first
triplet of the full validation ensemble, so it is the realization the scored
table also contains.  The two frozen operating points drawn as horizontal
reference lines are calibrated from the noise-free low-sea and scaled high-sea
reference records exactly as the fixed-tuning modes are.

The scored ablation tables are produced by the full driver
(`python3 tools/ou_rs_law_ablation.py --exponents 0,0.5,1,1.25`), which writes
this same figure at the end of its run.
