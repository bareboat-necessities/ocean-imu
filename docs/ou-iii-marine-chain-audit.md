# OU-III marine certificate audit

The source reference for this audit is main at
44d08d1c1f7c4ef1c7c7d4c7ad28e895f4db52dd, including the removal of
unreachable staged-MEKF startup code. The numerical changes are in proof tools;
filter gains, state updates and deployment bounds are unchanged.

**The complete P1–P5 stability certificate is not established.** A successful
test or a subcertificate's PASS must not be interpreted as finite inner capture.

| Stage | Established scope | Remaining issue for use by the current filter |
|---|---|---|
| P1 | Conditional startup/handoff norm bounds and heading branches | The narrower P5 half-Hs position entrance is an additional assumption, not a consequence of the 20 m P1 position ball. An ungauged timeout has only a tilt bound. |
| P2 | An overapproximating source language | The 800-node untimed graph admits all 640,000 transitions. It cannot justify restricting tuner movement in P4. |
| P3 | Linear information-word subcertificate | Regenerate with the actual horizontal S covariance; a positive information margin alone says nothing about a useful nonlinear radius or noise floor. |
| P4 | Finite-angle geometry and a microscopic local seed | The full 18/21-state nonlinear word dissipation, including cross blocks and source-consistent covariance evolution, remains open. |
| P5 | Entry into the outer attitude sector | No finite inner capture count is established. An outer entry count of zero is not an inner capture result. |

## Source-map repairs

The source stores WORLD-to-BODY attitude and injects the estimate correction on
the left: R_hat_plus = Q(dx) R_hat. For the physical error convention used by
the residual backend, E = R_true R_hat^T, this implies

    E_plus = E Q(dx)^T.

The corrected AD and full-matrix prefix maps therefore apply the inverse
correction on the **right**. The covariance reset continues to follow the
source's G = I + 0.5 [dx]_x. A noncommuting C++ comparison is essential;
collinear rotations cannot detect this side error.

For the pseudo-measurement, r = -S_hat. In true-minus-estimate coordinates,

    r = delta_S - S_true,    delta_S_plus = delta_S - K r.

The homogeneous map has r = delta_S. The former -delta_S made a scalar
P=R=1 update expand error 1 to 1.5 instead of contracting it to 0.5.
Physical S_true is an external forcing term, not zero merely because the
filter regularizes S towards zero. Prefix propagation now carries S_hat
separately from the physical error, including its zero initial value at first
handoff.

The finite gyro-bias prediction uses the source quaternion product
Q(-(omega-delta_b)h) E Q(omega h), rather than calling a first-order bias
approximation exact. The supplied omega is the estimate's bias-corrected rate.

The wrapper scales horizontal S standard deviations by its binary32
R_S_x_factor and R_S_y_factor (currently 0.72). Their variances use the squares
of those factors. Both the prefix gain calculation and the P3 lower noise
comparison now use these source values. The generic accelerometer H bound
also follows the estimated a_w mean, rather than treating a physical truth
force cap as an estimate cap. Effective-input norms use the entire a_w vector,
and attitude prediction consumes the propagated gyro-bias error rather than
reusing the startup bound indefinitely.

## Tuner timing and unnecessary widening

The source updates each EMA on every valid sample, tests a strict
time-last_stage > 0.1f condition, and applies a staged candidate before the
next sample. The clock helper keeps the binary64 absolute clock and last-stage
time separately and rounds dt to the source binary32 type. At an initial
zero clock and nominal 0.005f samples, the first staging is sample 21 and
its activation is sample 22. This must not be replaced by a 20-sample jump.

The EMA enclosure now retains the dependence between a and 1-a by evaluating
the multi-affine vertices. In particular a constant state equal to its target
cannot spuriously halve or double. A separately justified finite elapsed-time
label can restrict its image; the original untimed graph still covers arbitrary
late commits. Integrating the clock, candidate tuple and committed tuple into
the P4 source path remains necessary. The mathematical EMA enclosure does not
by itself qualify target-platform expf or floating-point roundoff.

Other conservative combinations must be tightened using an enclosure of the
same physical/source set, not by shrinking deployment assumptions:

* Keep force coupled to a_w and gravity, and use the nominal mean in H.
* Preserve declared vector norm balls. Their enclosing Cartesian cubes admit
  simultaneous axis extrema; dividing every component by sqrt(3) is not a
  valid replacement because it excludes permitted axis-aligned states.
* Keep a shared frequency/sea-time input for tau and sigma horizons, and the
  separate R_S target/horizon law. Neither full target clamps nor arbitrary
  tuner jumps establish a physically typical sea.
* Preserve covariance cross terms and accepted-packet recurrence. Requiring
  every individual correction to contract is stronger than word dissipation.

## Marine coverage and finite capture

The declared Live envelope permits non-gravitational CoG acceleration up to
4 m/s^2, giving a derived specific-force norm range 5.80665–13.80665 m/s^2,
and body rate up to 30 degrees/s. These are conditional motion assumptions,
not hardware saturation specifications or a universal rough-sea guarantee.
Hs <= 8.5 m alone does not bound acceleration, angular acceleration or slamming.

The domain covers the vibration guard only while it is transparent: zero
initial engagement and detector RMS <= 0.03 m/s^2. Active/transitioning guard
dynamics, nonzero lever arm, impacts/slams and saturated A-mode bias states
remain outside this certificate. A lever-arm extension needs the source
angular-acceleration and centripetal terms and their correlated uncertainties;
simply raising the CoG acceleration cap does not establish it.

A finite P5 capture proof needs a verified outer-word recurrence on the whole
entrance family and every prefix in its validity domain. For a pathwise bound

    V_next <= rho V + b,    0 <= rho < 1,

an inner target must lie strictly above b/(1-rho) for uniform finite capture
from above. An expected-value recurrence alone is not a pathwise capture
certificate. The stochastic domain currently supplies a failure budget of 0.05,
but no completed word-input budget, finite-horizon localization calculation
or validated inner capture count. The microscopic local seed cannot simply be
declared a practical stochastic target.

## Reproduction

The source-map audit workflow compiles the shipping C++ source and compares
its correction, covariance reset, forced S update, finite gyro-bias prediction
and tuner staging with the proof operations. It also runs the existing prefix
tests and regenerates P3. Numerical parity tolerances are diagnostic, not
theorem margins or a replacement for validated arithmetic.

On commit 1e3b49745a204a3ff690454a39c5522661b093ca, CI passed 18 source-map/AD
tests and 16 existing prefix tests. Regenerated P3 passed for H and A with
endpoint margin 2.2953997386276595e-20. The 30-degree full-word attempt remained
NOT_ESTABLISHED: its first H accelerometer enclosure at sample index 0 exceeded
the validated correction range [0,6] radians; A was not reached. This is an
enclosure obstruction, not evidence of an actual six-radian filter correction.
See runs 33567318416 and 33567318226. Later timing/convexity edits require their
own CI result.

Correct source covariance changed the reference-cell maximum gain tightening
from the old test's assumed >1.5 to 1.4088643079177852. The regression contract
now requires improvement (>1) and no worse bound on every cell (>=1), matching
the lemma; the analytic supremum and fixed-cell tests remain.
