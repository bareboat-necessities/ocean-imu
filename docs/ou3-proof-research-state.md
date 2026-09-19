# OU-III proof research state

## Current hypothesis

Seek regional practical stability of the shipping estimator on one persistent
execution satisfying MARINE MOTION, IMU BIAS and MAGNETIC SERVICE. The same
path must carry history-dependent capture, magnetically informed H18, the
implemented release and magnetically informed A21. The theorem remains open.

The H18 leg can no longer be posed as strict contraction of the full
21-coordinate finite error. The held accelerometer-bias coordinates are exactly
reproduced by the shipping H18 superword map, so the working hypothesis is now
that H18 dissipation holds on the complement of those coordinates, with the
held bias carried as a bounded input admitted by `|e_b| <= B_a + R_b`.

## Evidence

A literal magnetically informed H18 service superword was measured on an
actually reached shipping state, and the candidate storage `V(e)=e^T P^-1 e`
was evaluated on it at 60 decimal digits. Reproduce with

    make -C tools/stability h18-superword OUTPUT_DIR=/tmp

The export runs the deployed startup, hands over at the deployed `goLive()`,
keeps H18 through the shipping external accelerometer-bias hold, and takes
central differences of the complete closed-loop error map through the shipping
code itself, so the map is the deployed one rather than a re-derivation.

Root at 63.85 s, 30 s after Live, one second of superword, ten actually applied
magnetic corrections:

| quantity | value |
| --- | --- |
| worst admissible ratio, shipping covariance metric | 3.2595 |
| worst admissible ratio, root metric frozen at both ends | 3.3798 |
| worst admissible ratio at the coarse difference scale | 3.2496 |
| largest prefix retention ratio | 3.2595, reached at the endpoint |
| applied magnetic information, smallest eigenvalue | 1.4051 against `mu_M = 1` |
| held accelerometer-bias block of the map, deviation from the identity | 0 |
| held accelerometer-bias cross-covariance, either endpoint | 0 |
| held accelerometer-bias covariance block change | 0 |

The two difference scales agree to 0.6 percent entrywise and the two ratios to
0.3 percent, so the result is a property of the map and not of the differencing.
The frozen-metric ratio matches the time-varying one, so the growth is
amplification by the shipping map rather than drift of the covariance that
weights it. The limiting direction carries 43.7 percent of its energy in the
accelerometer bias and 37.0 percent in the integral state, with the attitude and
gyro-bias shares below 0.05 percent: the magnetic service condition is met with
margin and is not what blocks contraction. The same measurement at roots 120 s
and 300 s after Live gives 3.5375 and 3.6473, and a two-second window gives
9.0771, so the growth is secular rather than a startup transient or a window
artifact.

Classification: mathematical failure of the candidate storage formulation, not
an estimator regression, not a numerical failure and not an instability of the
deployed filter. The invalidated hypothesis was that a finite-error inequality
with `rho < 1` could be established for the full H18 error coordinate.

The failure is structural rather than quantitative. In H18 the shipping
estimator applies no accelerometer-bias mean dynamics, freezes the bias rows of
every gain, and leaves the bias cross-covariances at the zero the hold
installed. Writing the error as `(e_o, e_b)`, the measured superword map is
`(e_o, e_b) -> (A e_o + B e_b, e_b)` and the covariance is block diagonal
against the same split with `P_bb` unchanged, so for `e = (0, e_b)`

    V_end(Psi e) = (B e_b)^T P_oo,end^-1 (B e_b) + V_root(e) >= V_root(e).

No superword length, sea state or amount of magnetic information can make that
ratio smaller than one. Equivalently, the map is block triangular with an exact
identity block, so it carries eigenvalue one and no time-invariant quadratic
storage admits `rho < 1` either. The obstruction is recorded as
`held_bias_non_contraction` in `tools/stability/ou3_theorem/finite_error.py`,
its shipping premises are pinned by
`tests/kalman_ou_iii/shipping_transition-test.cpp`, and the diagnostic reports
it alongside every measurement.

Coercivity is a separate open concern on the same measurement: the reference
finite error sits at `V = 6512` at the root and `6622` at the endpoint, so the
reached error is far outside the ellipsoid the shipping covariance describes,
and the model mismatch carried as supply is large.

The architecture suite, the repository audit, the shared workflow and
replay-fingerprint contracts, the finite-error non-promotion regressions and
the native shipping-contract and shipping-transition tests pass. Source and
operation provenance validation passes. These establish implementation
regressions and subordinate algebra, not a source-uniform finite-error
certificate.

The local full build could not complete: `ensure-sim-data` failed because the
release archive host did not resolve in this environment. That is an
environment and data-acquisition failure, not evidence about the estimator.
The stability tooling itself needs no simulation archive.

## Current limiter

There is no restated H18 dissipation target yet. The full-state inequality is
refuted, and its replacement -- strict dissipation on the complement of the
held accelerometer-bias coordinates, with the held bias as a bounded input and
a supply constant that survives the `B e_b` coupling -- has neither a candidate
storage nor a feasibility measurement. Capture, actual H18-to-A21 retention,
A21 dissipation, recurring informative service and target arithmetic remain
open. Passing point-audit flags close none of these.

## Failed approaches / DEAD_ENDS

- Strict contraction of the full H18 error coordinate. The held
  accelerometer-bias block of the shipping superword map is exactly the
  identity and its covariance block is exactly frozen, so the worst admissible
  ratio is at least one by construction and was measured at 3.26 over one
  second. Sharpening enclosure bounds cannot cross that threshold; the
  formulation has to change.
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

Added by this measurement: the shipping H18 hold is exact, so the held
accelerometer-bias error is an invariant coordinate of the superword map, its
covariance block is constant, and it reaches the rest of the state only through
the coupling column `B`. Recurring magnetic service was available with margin
on every window measured, so heading information is not the H18 limiter. The
estimate projection bound `|e_b| <= B_a + R_b` remains the available handle on
the held bias, which is what makes carrying it as a bounded input plausible.

## Alternatives

1. Restate H18 dissipation on the complement coordinates with the held bias as
   a bounded input, using the projection bound as its admissible amplitude and
   the measured coupling column as its input matrix. This is the direct
   consequence of the obstruction and is the first candidate.
2. Enter the certified tail only at the H18-to-A21 release, where the shipping
   accelerometer-bias prediction multiplier is `phi_OU < 1` and the block is a
   genuine contraction, leaving informed H18 inside the finite pre-certified
   prefix. This changes where the tail begins and has to be reconciled with the
   declared proof path before it is adopted.
3. Inverse-free information storage, finite-error attitude-group/additive-state
   storage, and bounded-real supply constructions retaining physical bias and
   model mismatch remain candidate methods for the same theorem, not parallel
   theorem paths.

## Next falsifiable experiment

Measure the complement map. Take the same reached H18 superword, restrict the
finite-error coordinate to the complement of the held accelerometer bias, and
evaluate the worst admissible ratio of the restricted candidate storage with
the coupling column `B e_b` moved into the supply. Report the restricted ratio,
its limiting direction, the induced supply gain against the projection bound
`B_a + R_b`, and every-prefix retention. If the restricted ratio is not
usefully below one, alternative 2 is the next formulation to measure, not a
sharper bound on alternative 1.

Only after a useful feasible margin is demonstrated does a source-uniform
enclosure start. Point consistency, a single reached root and a finite window
are none of them a certificate.
