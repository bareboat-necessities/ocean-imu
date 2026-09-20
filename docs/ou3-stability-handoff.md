# OU-III stability handoff

Continue only in PR #560, branch `ou3-explicit-regional-practical-stability`.
Read AGENTS.md and `ou3-proof-research-state.md`. Keep the one construction ->
capture -> H18 -> release -> A21 theorem and its same-history semantics.

The exact complete-word energy identity, full-factor representation and
singular nuisance range elimination are retained. The corrected LIN path
supplies a full-state singular lower comparison; convex combination with
fresh AG/BA injection closes the joint lower covariance at recurring regular
prediction roots after 16 s. See `root_covariance_certificate.py`.

`ou3-nuisance-upper-proof.md` now proves a recurring upper comparison for the
15 LIN/BA coordinates after 17 s of default regular A21 operation, including
all nuisance cross covariance. Three actual S observations eliminate the
unknown neutral root. The comparison keeps the same realized PSD-sync
increments while omitting acc/mag corrections. Do not recompute sync from
its auxiliary covariance. Constants are exact rational; no rho is inferred.

The remaining upper-covariance block is AG6. Full transported loss,
nonlinear radius, arithmetic supply, finite capture and retention remain
open. A finite trace cannot close these obligations.

`ou3-stationary-detectability.md` establishes stationary A21 detectability
with nonparallel gravity/magnetic field and recurring applied observations.
Do not exclude quiet water on an unproved assertion of instability. Any
necessary motion restriction must provide quantitative window excitation,
not merely nonzero amplitude. The active bias OU law must not be replaced
by a held-bias identity when assessing this question.

Reproduce with:

```
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py
```

Passing these checks is not end-to-end theorem completion. Do not revive
entrywise Riccati subdivision, endpoint-batch floors or restricted-information
lifting. The literature map and finite quiet-word diagnostic are removed.
