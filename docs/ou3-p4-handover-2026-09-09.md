# OU-III P4 continuation status

Start with `docs/ou3-proof-research-state.md`,
`docs/ou3-live-entry-audit.md`, and the physical-wave-generator lemma in
`doc/kalman_ou_iii/w3d-brmm-stability-theorem.tex-part`.

## Corrected source, not an altered estimator

COMPLETE-BRMM now requires wave displacement relative to an equilibrium and a
qualified spectral or bounded shaping realization. Its potential derivative is
p, so `S_L(t)=phi(t)-phi(t_L)` follows by the fundamental theorem of calculus.
The exact backend is `ou3_brmm_physical_wave_source.py`; the correlated boundary
binding is `ou3_brmm_wave_primitive_binding.py`. The latter is a necessary outer
bridge, not a substitute for continuous generator/output qualification.

The source's numerical inverse-frequency amplitude or shaping-state budget is
not yet qualified for the entire family. Consequently D_S and the comparison
with the 300 m*s working radius remain unresolved. Do not choose D_S by searching
that working radius. The old `p=+/-d` obstruction stays in
`ou3_brmm_infinite_continuation.py`: B for the old source, E for the intended
physical omission, and excluded by the corrected physical theorem.

## Retained execution and proof graph

The fresh entry removes the common S origin once. Position is not re-anchored;
S is not reset at later words. Zero estimator means at fresh Live do not imply
zero physical displacement or velocity. The measured-period and prior-frequency
Live paths remain distinct and both need complete capture coverage.

The event attachment uses the code-owned same-signal JOINT frontend, including
prior-frequency evolution, after the current IMU events and before external
magnetometer handling. Every event cell retains its typed physical-potential
payload. The codec preserves absent prior WPE values as null; the legacy canonical
typed-window executor still requires the measured-period slice. Complete
prior-root windows and continuous physical output attachment remain open;
structural payload checks cannot close them.

Retain branch-complete production lineage, source-indexed H18/A21 prefixes,
separate BIAS0/1/2 recurrences, true-bias ancestry, same-history P/H/R/K,
Joseph/reset/projection splitting, the candidate/active tuner and scheduler,
and conditional finite-precision inputs. Use compatible joint24 storage; do not
resurrect the failed A21 marginal18 motion metric.

## Fail-closed continuation

P3 delta is `1e-18`. The corrected physical lemma removes the obsolete current
DC-source obstruction but does not close numerical source qualification, the
601-sample source-uniform cover, joint24 forcing attachment, endpoint and every
literal-prefix augmented LDLT, first-exit retention, a maximum P4 basin, finite
H18 capture, or H18/A21 basin transport. P4/P5 remain false.

Regenerate the source, centered-S and final-gate JSONs. Run the physical-source,
old obstruction, correlated boundary, prior-frequency attachment and exact
source-indexed rebase regressions. Keep any failed innovation enclosure separate
from a physical counterexample. The research ledger records the current limiter,
critic alternatives and next falsifiable complete-word experiment.
