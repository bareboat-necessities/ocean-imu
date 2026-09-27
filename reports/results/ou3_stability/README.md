# OU-III stability evidence

The status remains fail-closed. Exact partial certificates do not establish
the end-to-end theorem.

- `lin-matrix-certificate.json`: full rational LIN endpoint precision and
  covariance comparison, plus a raw neutral S normalization bound.
- `root-covariance-certificate.json`: positive joint 21-state lower covariance
  from a convex combination at the same post-prediction root.
- `word-energy-audit.json`: exact counterexample to restricted-service lifting.
- `block-factor-status.json`: assembly status with no promoted full-state rho.
- `theorem-status.json`: controlling obligations and current next step.
- `provenance.json`: shipping and proof source hashes.

The evidence builder recomputes the matrix, audit, and assembly artifacts and
compares them exactly; a changed certificate cannot pass with stale evidence.

```
python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json
```

The article uses the exact complete-word covariance-energy identity. A source-
uniform full-state loss bound, nonlinear radius, finite capture, retention,
physical qualification and float32 supply remain open.

## Physical regimes and conditional moving pivots

`regime-certificate.json` exactly reproduces the smooth sin-cubed rest/motion
sensor-indistinguishability witness, its unchanged bias/rate margins and actual
magnetic-service reserve. It records the complete-window excitation quantifier,
conditional bridge algebra and failure of finite bridges alone under repetition.
`moving-pivots-certificate.json` checks the chronological gyro defect recurrence,
two-group six-column implication and all six greedy pivot floors using rational
arithmetic and an independent supplied matrix audit. These are conditional
lemmas, not source-uniform moving or stationary stability certificates.

The proof-side quiet-evidence monitor never certifies physical STILL. No
shipping mode switch, noise, cadence, magnetic semantics or quality gate is
changed. The native `regime_ambiguity-test` carries one 900-second construction/
rest/hidden-rocking/rest execution through Live, refinement and release. It is
a finite regression. Existing source replay artifacts remain byte-for-byte
unchanged because their shipping sources and diagnostic inputs are unchanged.
Only changed proof contracts, exact certificates and their bindings are renewed.
