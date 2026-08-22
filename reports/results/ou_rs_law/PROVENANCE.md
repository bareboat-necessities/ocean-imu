# OU-III round-trip transition figure provenance

`ou_rs_roundtrip_transition.svg` and its decimated series
`ou_rs_roundtrip_transition.csv` are one bidirectional transition realization
replayed through the current deployed OU-III adaptive filter.  The figure uses
`tools/ou_validation.py:plot_transition_diagnostic`, the same drawing code as
the primary one-way transition diagnostic.

Regenerate with:

    python3 tools/ou_roundtrip_transition.py \
        --output-dir reports/results/ou_rs_law

and mirror the SVG to `doc/kalman_ou_iii/ou_rs_roundtrip_transition.svg` for the
manuscript.  The plotted realization uses seed triplet 11/101/1009.  The low and
high reference lines are independently calibrated frozen operating points for
the two endpoint records.

This artifact is retained only as a diagnostic of bidirectional response of the
current scheduler.  It is not evidence from, or a comparison among, historical
regularizer-law ablations.
