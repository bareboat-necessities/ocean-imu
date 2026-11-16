# Detailed heave comparison — Southbound, head / bow-quartering seas

All 136 full-record gates passed. [Metrics](heave_comparison_metrics.csv) use the
final 900 s of each 1200 s record. [Receipt](case-receipt.json) records inputs,
executable hashes, commands and completion; [summary](comparison.json) records
the trace and chart inventory.

The complete study produced 136 traces and 102 SVGs across shipped Mahony,
IMU-only proxy and oracle true-attitude frontends. The six production-default
overviews below are retained here; the full output is uploaded by the
`simulation-studies` workflow. Legends report RMS in the plotted final 60 s,
not the full scoring window. Traces preserve 200 Hz before 50 Hz plot sampling
and bounded vector export. All numeric gates use the original full-rate data.

| Spectrum | 0.27 m | 1.5 m | 8.5 m |
| --- | --- | --- | --- |
| jonswap | [Chart](charts/heave_only_jonswap_low.svg) | [Chart](charts/heave_only_jonswap_medium.svg) | [Chart](charts/heave_only_jonswap_high.svg) |
| pmstokes | [Chart](charts/heave_only_pmstokes_low.svg) | [Chart](charts/heave_only_pmstokes_medium.svg) | [Chart](charts/heave_only_pmstokes_high.svg) |

These verification runs started at the HEAD recorded in the receipt with the
implementation changes subsequently published in commit
`77ec3f40dd0a63d953920f40616e2a8367322697` present in the worktree. Their original
starting HEAD and actual executable hashes are retained; no replay provenance
has been relabelled. CI runs use their checked-out commit directly.
