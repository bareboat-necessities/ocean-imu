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

`stationary-covariance-certificate.json` supplies a uniform historical action
and every-operation full covariance upper comparison on the precisely stated
zero-residual quiet nominal subcase, plus qualitative homogeneous linear loss.
It does not certify nonlinear stability of the physical compatible class.
The updated moving certificate includes the literal reset-inverse identity,
an inverse-frame prefix budget, and a relaxed reset cancellation witness.
`regime-continuation-feasibility.json` is an 80-digit supplied-algebra audit,
not source-uniform coverage, reachability or a finite fitted contraction.

`moving-transport-source-feasibility.json` freshly replays the existing carried
quiet/moving source driver with read-only instrumentation and untapped control
parity. Exact rational exported operators factor both selected row groups with
E=0; 80-digit finite six-column budgets are about 1.20433 and 1.07449.
The evidence gate checks source/observer/driver fingerprints and fail-closed
scope; reproducing the numbers requires the optional native driver below.
Neither exported operands nor finite positive singular values enclose all
physical histories or certify magnetic service:

```
python3 -m tools.stability.ou3_theorem.moving_transport_source_diagnostic \
  --eigen /usr/include/eigen3 --output /tmp/ou3-moving-transport.json
```

`world-frame-certificate.json` exactly checks the world-frame row
factorization on a rational word, the reset Gram identity, attitude-invariant
same-cell geometry, the literal injection Loewner budget, the generalized quiet
floor, Corollary A's AW threshold, the jerk-limited collinear cadence and the
1-Hz collinear motion/bias witness (no MAGNETIC SERVICE admission claimed).
`aw-covariance-ceiling-certificate.json` exactly checks Lemma B: the
isotropic-sync AW ceiling `(1+eps)16`, its prediction invariance and excess
decay, the spectral-max sync witness, the S_factor=2 counterexample and the
Corollary A storage radii (.280954 at the clamp).
`world-frame-source-feasibility.json` replays quiet, wave and collinear source
words (1-Hz and 25-Hz magnetic cadence) through a derived observer with
untapped control parity. It records floor/actual ratios <=1 after the float
mean-injection charge, injection bound ratios, the collinear same-cell
collapse (.0301) beside aggregate sigma 28.50, the failed one-correction
service information (2.03e-5), the signed injection sum, the literal AW
block reconstruction and the failed uniform storage route (6.68 and 6.91).
CI reproduces the committed record exactly. These are non-promoting finite
audits:

```
python3 -m tools.stability.ou3_theorem.world_frame_source_diagnostic \
  --eigen /usr/include/eigen3 --output /tmp/ou3-world-frame.json \
  --expect reports/results/ou3_stability/world-frame-source-feasibility.json
```

`aw-tracking-certificate.json` exactly checks Corollary A* (attitude
columns from the nominal signed AW mean, threshold 1.96133 m/s^2), the
signed physical transfer, the weighted signed-correction representation of
the nominal mean, the world innovation factorization, the gain-weighted AW
loop identity and the AW increment cap. `aw-tracking-source-feasibility.json`
replays six admitted C2-onset histories through the unchanged estimator with
a one-tap observer and untapped control parity: the pointwise AW premise
fails by 6.80, the worst signed mean is .330 of 1.12383 (sync-locked
rectification) and the worst nominal transverse mean .177 of 1.96133. It
also rebuilds the literal world six-column array with every reset over a
100-s word: literal sigma_min 191.8--267.4, never below the injection-free
array, `|A~-I|<=.0045`, and G0's floor .072--.084 below both.
`signed-injection-certificate.json` checks Lemma I* (ordered injection
rotation bounded by endpoint attitude errors plus integrated rate error),
the third-order reset factor, the rotating-frame floor and the injection
feasibility table. `aggregate-floor-certificate.json` checks Lemma T and the
injection-free Theorem G0 (`s^2>=1.486786e-3` under supplied nominal window
premises), a synthetic falsification audit and the downstream feasibility.

```
python3 -m tools.stability.ou3_theorem.aw_tracking_source_diagnostic \
  --eigen /usr/include/eigen3 --output /tmp/ou3-aw-tracking.json \
  --expect reports/results/ou3_stability/aw-tracking-source-feasibility.json
```

The proof-side quiet-evidence monitor never certifies physical STILL. No
shipping mode switch, noise, cadence, magnetic semantics or quality gate is
changed. The native `regime_ambiguity-test` carries one 900-second construction/
rest/hidden-rocking/rest execution through Live, refinement and release. It is
a finite regression. Existing source replay artifacts remain byte-for-byte
unchanged because their shipping sources and diagnostic inputs are unchanged.
Only changed proof contracts, exact certificates and their bindings are renewed.
