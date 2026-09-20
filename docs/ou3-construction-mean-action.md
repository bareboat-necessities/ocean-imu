# Construction-linked nominal mean action

This attempt enters the one finite-error tail inequality through the required
source-uniform historical AG reader. Before an action ceiling B_* can be
certified, the construction and subsequent nominal recursion must exclude
sustained force/field collinearity and gyro aliasing. The following is an
executed necessary separation test, not a replacement contraction argument.
It closes a finite recorded gyro barrier and rejects an energy-only route to
the simultaneous uniform exclusion. All source-uniform claims remain open.

## Exact joint recurrence

In the default profile let u=(b_hat_g,a_hat_w), a six-coordinate **mean**
vector, distinct from the six attitude/gyro error coordinates AG. Between
corrections its literal real-arithmetic predictor is

`u^- = A u^+`, `A=diag(I3,phi_aw I3)`.

At an accepted rank-three correction it is `u^+=u^-+K_u r`. Here K_u contains
the actual six rows of the full gain. In particular, nuisance cross covariance,
correlated process noise, varying parameters and previous resets remain in
K_u and the actual innovation covariance S. The attitude reset acts as the
identity on u. Bias projection changes BA only. A frame or other hard event
that changes u must be included explicitly; the recorded construction checks
every seam and fails on any unrecorded change. No estimator or physical state
is restarted at an audit boundary.

Start at the actual zero mean construction and set D=0, E=0, eta=0. Propagate

| Operation | D | E | eta |
|---|---|---|---|
| Prediction | A D A' | E | A eta + d_pred |
| Applied correction | D + K_u S K_u' | E + r' S^-1 r | eta + d_corr |

The d terms are the exact discrepancies between the literal recorded mean
update and the corresponding affine map with its exported operands. They
vanish in the ideal real-arithmetic recursion. Skipped observations add
nothing. An explicit factorization S=L diag(s) L' evaluates each addition
with three columns K_u L, without an inverse of the full covariance.

If S is symmetric SPD, then at every prefix

`[[D, u-eta], [(u-eta)', E]] >= 0`.

To prove this, stack the transported columns `A_(t<-i) K_ui S_i^(1/2)`
and the inputs `S_i^(-1/2) r_i`. Their Gram matrix is the displayed block.
The mean identity uses the same realized coefficients throughout. The full
6x6 D retains BG/AW cross terms; no independent coordinate disturbances or
diagonal covariance replacement is made. The action is accumulated directly
from the actual gains. It does not require a covariance-subtraction identity
across floating updates or hard covariance replacements, including events
that preserve the means but clear cross covariance.

This is a bound on nominal means, **not** the readout noise action B_W. An
upper bound on its scalar input energy E alone says nothing about uniform
observability. In particular, the innovations are not physical sensor noise.

## Executed construction and enclosure

The observer runs the existing diagonal-wave history from begin to 600 s,
through the actual frontend, handoff, H18, reference refinement and release.
All terminal means, quaternion, covariance, stage times and magnetic counts
are identical to the untapped control. It records 89998 predictions and
122251 actually applied rank-three corrections; it does not resample the
input or reseed means/covariance. The field in the test below is the actual
committed reference after refinement, not a substituted physical field.

First the complete action is evaluated with 80-decimal arithmetic. Then the
same binary operands are processed with outward 40-decimal intervals, using
only additions, products, divisions and positive-pivot LDL. The recorded S
is replaced by its exact symmetric part **only as the SPD action weight**;
the actual gain and innovation remain unchanged. This does not assert that
the literal covariance update is an exact Kalman update. All interval
endpoints are saved as exact dyadic rationals. A separate rational interval
LDL calculation rechecks the final six-coordinate action and exclusions.

BG has identity prediction, so its action block D_gg is Loewner-increasing,
and E is nondecreasing. Let eps bound every prefix of ||eta_g||. The full
matrix test

`(1-eps)^2 I3 - E_final D_gg,final > 0`

therefore proves ||b_hat_g||<1 rad/s at **every recorded prefix**. This ceiling
is derived from the recorded action, not assumed. Together with the existing
physical gyro and residual envelopes it gives

`h ||omega_measured-b_hat_g|| < .0099051914291880918 < .01`.

Thus this actual finite construction has a certified separation from a
complete-turn alias. It does not bound every admitted history or enclose the
floating implementation of the Rodrigues/integral helper.

An additional exact check uses the actual recorded pre-accelerometer means,
not a freely selected point in the ellipsoid. For all 40000 such observations
in the 400--600 s tail it verifies
`25 ||(aw-g) cross B||^2 - 4 ||aw-g||^2 ||B||^2 > 0`.
The smallest polynomial margin is
`317402241109776656710823995905/302231454903657293676544`.
Thus the actual recorded tail has sine >2/5 throughout. The energy enclosure
failure below must not be described as a nonphysical trajectory blocking the
theorem, or as observed sustained collinearity.

For force/field collinearity choose a full-row-rank 2x3 matrix C with CB=0
for the committed nonzero B. Collinearity is `C a_hat_w=C g`. If D is SPD,
the least energy in the **joint** action ellipsoid allowing this equality is

`E_col = [C(g-eta_aw)]' [C D_aw,aw C']^-1 [C(g-eta_aw)]`.

The other coordinates are minimized out by a full-covariance Schur identity;
their cross terms are retained, rather than set to zero. Separation would
require `E_col-E > 0`. The executed 80-digit values are

`E=817885.0623259853033`, `E_col=1.0138313892264883`,

so the margin is approximately **-817884.0484945961**. The rational interval
certificate confirms the strict negative sign. This energy ellipsoid admits
collinear means even on the recorded physical history whose actual tail
force/field sine stays above .443. It is the enclosure that fails, not a
constructed sustained-collinearity trajectory. All-time magnetic service
is still unproved for this run.

## Architecture decision

The roughly 8e5 deficit cannot be repaired by rounding precision, a smaller
pivot interval, or coefficient subdivision. Energy loses the signed and
time-linked restrictions in the literal innovations

`r_acc=f_measured-R(aw-g)-ba`, `r_mag=m_measured-R B`, `r_S=-S_hat`.

Do not refine this energy ellipsoid into an asserted common nominal-state
box. The next admissible technique must retain those signed relations with
the same physical velocity/primitive and nominal integral chain. Its
falsifiable target is a positive separation margin over a whole carried
window, while also controlling gyro transport, followed by the complete
historical reader action. Neither a physical Gramian nor the finite gyro
certificate supplies this missing implication. Uniform B_*, J, rho0, radius,
capture/release and every-prefix tail retention remain open.

Reproduce the diagnostic with
`python -m tools.stability.ou3_theorem.construction_mean_action --export --directory /tmp/ou3-mean --output /tmp/ou3-mean-diagnostic.json`
and then the enclosure with the same directory, `--intervals`, and a distinct
output. `--eigen` selects an alternate include directory. The exported trace
is scratch evidence bound by its hash; the compact factor certificate is
committed and checked independently by the normal evidence validator.
