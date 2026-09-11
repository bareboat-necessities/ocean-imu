# Conditional finite state/covariance composition

## Obligation and boundary

This lemma advances the finite-map stage of the ALT plan: consecutive accepted
core events must consume the SAME finite mean, full covariance and physical
predecessor. It does not close the source-uniform 600-step runtime word or
permit a storage search. The implementation is
`tools/stability/ou3_alt_contraction/finite_core.py`.

The supplying one-event identities are in
[the finite measurement proof](ou3-alt-finite-measurement-proof.md), including
physical prediction and the masked rank-three Joseph identity. They are reused,
not replaced by a second physical model or a derivative cocycle. No deployed
header, original P2/P3/P4/P5 premise or final gate is changed.

## Paired core state

Carry `X=(mode,z,P,q_hat,reference)`, where
`z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)` has 24 coordinates and P is the full
symmetric 21-state shipping covariance in BOTH H18 and A21. Require

```
c = Cayley(q_phys conjugate(q_hat)),
beta = reference.beta.
```

The reference contains the physical endpoint time, quaternion, v/p/centered-S/a,
gyro bias, accelerometer bias and one persistent Live origin. `Reference`
inherits `PhysicalKinematics`. `PhysicalSegment` still checks the SAME moment
and bias recurrence, and its before endpoint must equal this exact predecessor.
History/family labels are only consistency checks; neither labels nor finite
endpoint consistency establish COMPLETE-BRMM/BIAS admission. The all-time
primitive bound, admissible startup state and physical angular defect remain
source obligations.

A prediction derives `omega_hat=gyro_meas-(bg_true-e_bg)` from this state,
generates its nominal shipping polynomial step, and calls the existing physical
finite predictor with the actual physical endpoint increment. It also carries
`P_next=F P F'+Q`. The actual per-axis mean coefficients and F/Q are explicit
**conditional operands**, not independent admitted source boxes. Their literal
runtime/tuner/repair equations still have to be attached; floors and other
intervening events are not implicit in this prediction function.

## Accepted-event successor

For an accelerometer event reconstruct `a_hat=a_phys-e_aw`,
`b_hat=beta-e_ba` and the nominal rotated specific force from q_hat. For a
magnetometer event rotate the supplied applied reference by that same q_hat.
Form H and the exact finite secant Hbar from those nominal quantities and c.
For S retain `r=e_S-S_phys`, including nonzero physical S at zero motion error.
Sensor/model discrepancy is computed against that same physical endpoint; it
is not a newly independently bounded supply.

Using the current P and the supplied applied R, form the existing masked N and
full innovation Sigma, retaining H18's latent BA contribution. A conditional
safe-LDLT diagonal repair delta is included in the SAME innovation used in
both equations and covariance update:

```
Sigma_delta = Sigma + delta I,
K Sigma_delta = N,    Sigma_delta q = r,
d = K r = N q.
```

Here Sigma_delta is symmetric and nonsingular. An exact rational elimination
checks these relations; it does not certify the Eigen acceptance predicate,
positive definiteness, the value of delta or a floating-point solve residual.
An indefinite nonsingular matrix is not falsely treated as an SPD certificate.

The finite descriptor supplies the error successor and nominal quaternion
injection from the SAME d. The radial projection factor is checked against the
same pre-projection bias estimate `beta-(e_ba-d_ba)`, not an independent beta.
The covariance successor is

```
P_J = P - K N' - N K' + K Sigma_delta K' = P - K N',
G_d = diag(I_3 + [d_theta]x/2, I_18),
P_next = G_d P_J G_d'.
```

The first equality reuses the checked masked Joseph lemma; no `N=P H'`
assumption is introduced. The reset is the literal shipping COVARIANCE law,
not a linearization substituted for the nonlinear finite MEAN map. Its
attitude block and every attitude/non-attitude cross block are retained.
On the symmetric real branch upper-triangle copying and symmetrization do not
alter this congruence. Bias projection changes the mean but does not project P.
Numerical asymmetry and solve/accumulation defects are separate open obligations.

The rational evaluator implements the strict polynomial quaternion branch and
rejects other branches. The general finite w/k identity remains in the supplying
module. This evaluator limitation does not remove the axis-angle branch from
the shipping theorem or authorize a smaller physical domain. Rejected, not-due,
nonfinite, floor and mode-release events are not replaced by identity maps.

## Conditional composition theorem

For any finite event sequence, assume each event has its actual same-history
operands and branch premises attached and satisfies its supplying finite
identity. Starting from an actual paired state X_0, recursive substitution of
these event relations yields the actual paired state at every represented
prefix. No fresh covariance, nominal geometry or physical origin is inserted.

Proof: the empty prefix is X_0. If X_j is the actual pair, its reconstructed
nominal state and current P give the actual H/Hbar/N/Sigma/residual for the
next represented event. The inverse-free and finite reset/projection identities
above yield its mean successor; the Joseph and covariance-reset identities yield
its covariance successor. A prediction instead uses the existing physical
segment identity and its conditional covariance recurrence. Thus the next pair
is actual under the stated event premises. Induction proves the claim for any
finite number of such events. Substitution is associative, so splitting a
sequence after any represented prefix gives the same suffix result.

`compose_accepted` implements the accepted-measurement part and exposes the
empty prefix plus every represented successor. It imposes no arbitrary cap on
sequence length, but admits no schedule. Prediction composition uses the same
paired-state API. This theorem concerns equality, NOT chart retention,
coercivity, decay or a source-uniform event cover. Even an exact telescoping
storage identity supplies no sign or rho bound by itself.

## Validation and next binding

The exact-rational regressions compare the paired successor with direct nominal
subtraction and dense full-21 Joseph/reset arithmetic. They exercise held and
active masks, a common repair shift, active bias projection, nonhomogeneous S,
physical predecessor mutations and prefix/suffix equality. Their algebraic
roots and BIAS labels are explicitly unqualified: passing these tests is not
physical-source evidence. The existing passive shipping correspondence runner
remains the implementation check; no additional forced-runtime-root harness is
introduced.

The next attachment must eliminate the conditional F/Q/applied-R/repair inputs
by binding them to the actual same-history runtime memory. It must also supply
frontend/WPE/bandpass/sigma/tuner successors, all prediction floors and repairs,
axis-angle and rejection branches, actual S scheduling, asynchronous magnetic
history and every H18/A21 edge. The physical angular-defect and all-time
primitive/BIAS constraints remain joint with this graph. Local accepted-event
prefixes are not every literal runtime prefix.

No common metric, rho, source-uniform finite word, retained basin, finite
capture time or deployment-roundoff enclosure is certified. All ALT final gates
remain false, and the original route remains independently continuable.
