# OU-III stability evidence

The status remains fail-closed. Exact partial certificates do not establish
the end-to-end theorem.

- `lin-matrix-certificate.json`: full rational LIN endpoint precision and
  covariance comparison, plus a raw neutral S normalization bound.
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
