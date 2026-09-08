# BRMM stability proof handover

Read `AGENTS.md` and `docs/ou3-proof-research-state.md` before continuing.
The current simulation source is the checksum-pinned v1.2.1 estimated 28 ft
sailboat RAO dataset, selected in `tools/sim_dataset.py`. All live replay
paths use vessel CG motion. Retained payload fixtures preserve their original
identities; they must not be restamped or described as new vessel replays.

BRMM is the only primary motion hypothesis. Keep same-history correlation,
actual-applied anisotropic R_S, full Q/covariance cross terms, finite resets,
physical -S_true forcing, zero proof lever arm, and the dormant vibration
branch. No filter or quality-gate changes are authorized by a dataset switch.
P3 remains delta=1e-18 for both H18/A21, conditional on the configured
Normal-Live matrices/caps and accepted-vector PE. Scalar BRMM recurrence does
not establish that PE. P4 and P5 are unclosed; favorable replay or forced
endpoint ratios do not promote them.

The full branch evidence workflow regenerates validation and robustness after
dataset changes. `simulation-studies.yml` runs noise-free model mismatch,
engine degradation, mitigation, cross-family guard checks and roundtrip motion.
The lever-arm workflow and build matrix regenerate installation studies,
ordinary simulations, comparison plots and LaTeX. Source/runtime audit reports
must record the actual archive, revision, phase partition and all violations.
Use their completed run identities; do not infer success from a queued job.

Physical-source qualification and a useful practical supply/storage bound
remain separate tasks. Continue from the ledger's falsifiable experiment,
without widening caps, changing proof premises, or promoting a theorem from
sampled extrema. The source-connected diagnostic uses the pinned RAO API and
an independent reconstruction of its shared kinematic primitive.
