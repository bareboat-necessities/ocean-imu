# OU-III proof research state

## Current hypothesis

Seek regional practical stability of the shipping estimator on one persistent
execution satisfying MARINE MOTION, IMU BIAS and MAGNETIC SERVICE. The same
path must carry history-dependent capture, magnetically informed H18, the
implemented release and magnetically informed A21. The theorem remains open.

## Evidence

The repository cleanup gate found 14 obsolete-architecture references across
seven files in CI run `35382314754`. Those remnants and the shared workflow
assertions have been corrected without weakening the cleanup gate. The replay
fingerprint is explicitly invalidated, not regenerated from un-replayed results.
Scientific simulation settings and shipping estimator behavior are unchanged.

The operation lemmas retain the joint bias state `[e_b; b_true]`, the single
shared physical-increment column, the separate Kalman correction and projection
defect, and the conditional Euclidean projection inequality. The integral
innovation retains `r_S=e_S-S_true` for one fixed physical-potential origin.
Rectangular magnetic sensitivity transport retains the original information
coordinates; its dimensional check does not certify the physical release map.

The focused Python suite has 63 passing tests, including the repository audit,
shared workflow and replay-fingerprint contracts, and finite-error
non-promotion regressions. Source and operation provenance validation passes.
The native shipping-contract, shipping-transition and common-math tests pass.
These results establish implementation regressions and subordinate algebra,
not a source-uniform finite-error stability certificate.

The local full-build command was
`make all EIGEN_DIR=/opt/pyvenv/lib/python3.13/site-packages/casadi/include/eigen3`.
Compilation completed, but `ensure-sim-data` failed when
`curl -fL --retry 3 https://github.com/bareboat-necessities/oceanography-waves-lib/releases/download/v1.2.1/sim-data-files-vessel-rao-28ft.zip -o /mnt/data/ocean-imu/sim-data-files-vessel-rao-28ft.zip.download`
returned exit 6: `Could not resolve host: github.com`.
Consequently the complete local simulation/data-dependent validation did not
run. This is an environment/data-acquisition failure, not evidence of a
mathematical or estimator regression.

## Current limiter

A complete source-uniform finite-error H18 service-superword inequality with
strict dissipation, inherited covariance and every-prefix retention is absent.
Capture, actual H18-to-A21 retention, A21 dissipation, recurring informative
service and target arithmetic remain open. Passing point-audit flags do not
close any of these obligations.

## Failed approaches / DEAD_ENDS

- A zero-radius proof projection erased the estimate, whereas the shipping
  operation disables projection. The proof now carries the actual branch.
  Invalid attitude injection can also bypass projection; successful finite
  injection remains an explicit domain premise.
- Squared-norm underflow classified tiny nonzero constant displacement as quiet
  water. Componentwise exact-zero admission preserves the mandatory DC rejection.
- Square-only sensitivity transport rejected dimension-changing error charts.
  Source and destination dimensions are now separate; no release certificate is
  inferred from this algebraic repair.
- A point-diagnostic constructor supplied infinite prefix bounds that its own
  validation rejected. Replacing those bounds by the observed values would be
  circular. Missing bounds are now explicit, empty prefixes are rejected, and
  finite endpoint/flag checks never set `certificate_complete`.
- Overflow of the scalar dissipation right-hand side cannot be used to pass an
  inequality. Nonfinite derived arithmetic makes the point audit fail.
- Strict full-state contraction on indefinitely ungauged intervals is blocked
  by the unit-spectral-radius heading/axial-gyro-bias unipotent block.
- Callback cadence, independent per-sample bias boxes, a finite replay, or
  imposing the estimator OU prior on physical bias cannot supply the missing
  magnetic, temporal or source-coverage premises.

## Retained facts

The bounded displacement potential and true IMU bias histories persist across
all estimator events. Physical reference acceleration remains in the specific
force or declared model disturbance. H18 and A21 share one bias prediction
relation; correction and projection do not reset truth. Actual applied magnetic
innovation covariance and complete preceding transport define information.
Euclidean bias projection compactness is not contraction in a full coupled
covariance metric. No estimator reset defines the certified-tail boundary.

## Alternatives

Before further enclosure work compare an inverse-free information storage, a
finite-error attitude-group/additive-state storage, and a bounded-real supply
construction retaining physical bias and model mismatch. These are candidate
methods for the same theorem, not parallel theorem paths.

## Next falsifiable experiment

Construct a literal magnetically informed H18 superword from an actually
reached state. Retain the full covariance, physical bias predecessor and every
shipping event. Export applied magnetic information and evaluate the complete
finite-error storage/supply residual and every-prefix retention margin at high
precision. Report the limiting direction and event contributions. Start a
source-uniform enclosure only after a useful feasible margin is demonstrated;
point consistency alone is not a certificate.
