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

## Signed time-linked exclusion construction

The next construction is **not** an unsigned innovation-energy bound.  Fix one
carried historical window and retain the chronological accepted acc/mag/S
updates.  At accelerometer event i write the literal signed relation

`r^a_i=f^m_i-R_i(\hat a_i-g)-\hat b_{a,i}`,

and at an applied magnetic event j write

`r^m_j=m^m_j-R_j B_j`.

The measured quantities on the left come from the *same* physical history:
`f^m_i=R^t_i(a_i-g)+b^t_{a,i}+n^a_i` and
`m^m_j=R^t_j B^t+n^m_j`.  Consequently no innovation may be replaced by a
free disturbance.  Substitution gives the signed defects

`R_i(\hat a_i-g)=R^t_i(a_i-g)+(b^t_{a,i}-\hat b_{a,i})+n^a_i-r^a_i`,

`R_i B_i=R^t_i B^t+n^m_i-r^m_i` at magnetic service events.

For a complete nonuniform physical cell use the already proved trapezoidal
weights alpha_i.  Since the displacement/velocity coordinates are carried
without restart,

`sum alpha_i a_i=(v_1-v_0)/T+q_a`,  `||q_a||<=J h_max/4`.

This is the required signed physical integral relation: the acceleration
samples cannot be chosen independently of their two endpoint velocities.
Likewise the bias samples are one rate-bounded history, so for any zero-sum
weights c_i, `||sum c_i b^t_{a,i}||` is bounded by the declared bias-rate
constant times the corresponding first absolute moment.  The same rule is
used for gyro bias.  No absolute values are taken before these cancellations.

The nominal side must be telescoped with the literal mean recursion.  For
`u=(\hat b_g,\hat a_w)`, predictions and corrections give, chronologically,

`u^-_{k+1}=A_k u^+_k+d^p_k`,
`u^+_k=u^-_k+K^u_k r_k+d^c_k`.

Thus for any terminal selector Z,

`Z u_N=Z Phi_{N,0}u_0 + sum_k Z Phi_{N,k+1}K^u_k r_k
       + sum_k Z Phi_{N,k+1}d_k`.

This identity keeps the **sign, time, sensor type, actual gain, OU transport
and hard-event order** of every innovation.  It is the relation that the
failed scalar energy ellipsoid discarded.  The same chronological products
must be used to transport the two vector equations above to a common root.

A sustained nominal force/field degeneracy over a service window means there
are unit directions z_k, transported by the actual nominal attitude, for
which both nominal cross products are small.  Contracting the two signed
sensor identities with those transported z_k and summing with the physical
cell weights gives one scalar/vector balance whose physical term contains
`(v_1-v_0)/T` and the nonvertical magnetic component, while its remaining
terms are exactly: physical sensor/bias residuals, the **signed** innovation
functional `sum W^a_k r^a_k+sum W^m_k r^m_k`, and literal mean/reset/rounding
defects.  The innovation functional is then eliminated with the telescoped
mean identity above; it is not bounded by `sum r' S^-1 r`.  The gyro-alias
case is treated in the same chronology: the Rodrigues bias-transport factor
`B(h,omega_measured-\hat b_g)` is paired with the signed magnetic/accelerometer
corrections that could have moved `\hat b_g` from its construction value.
A complete-turn value cannot simply be inserted as an independent nominal
coefficient.

The innovation elimination itself has a compatibility condition that must be
proved, not assumed.  For the chronological affine recursion

`u_{i+1}=A_i u_i+K_i r_i+d_i`,

and signed physical-balance weights `W_i`, introduce matrix multipliers
`Z_i`.  Exact summation by parts gives

`sum W_i r_i = Z_N u_N-Z_0 u_0
 +sum (Z_i-Z_(i+1)A_i)u_i
 +sum (W_i-Z_(i+1)K_i)r_i
 -sum Z_(i+1)d_i`.

Therefore endpoint telescoping eliminates the interior nominal-state and
innovation terms only if one constructs the same-history adjoint sequence

`Z_i=Z_(i+1)A_i`,  `W_i=Z_(i+1)K_i`.

Solving the second relation independently at each correction is insufficient:
it must be compatible with the first relation through every intervening
prediction and hard event.  If this exact compatibility cannot be obtained,
the residual two sums must remain in the separation inequality and be bounded
without unsigned innovation energy.  This is now an explicit prerequisite to
any claimed `Delta_col` enclosure.

Measurement corrections also jump the nominal integral chain.  For a smooth
normalized nonnegative proof weight psi with psi=psi'=psi''=0 at both window
endpoints, piecewise integration by parts gives

`int psi a_hat dt = -int psi''' S_hat dt - C_psi - D_psi`,

where

`C_psi=sum_j (psi_j K_v,j-psi'_j K_p,j+psi''_j K_S,j) r_j`

uses the three rows of the **same actual** correction `K_j r_j`, and
`D_psi` retains source/interpolation, arithmetic and applicable hard-event
defects.  Physical boundedness cannot be substituted for boundedness of
`S_hat` or this signed correction charge.  This identity supplies a concrete
whole-window target without inventing pseudo-observations between callbacks.

This yields the correct contradiction target.  For every admitted carried
window W define `Delta_col(W)` as the lower bound on the physical
force/field term minus the sensor/bias, quadrature, mean-defect and signed
innovation-transfer remainders, and `Delta_gyr(W)` analogously for the two
transverse singular values of the integrated gyro transport.  What is needed
is

`inf_W Delta_col(W)>0`,  `inf_W Delta_gyr(W)>0`.

Only after those two **same-history** margins are certified may compactness of
the remaining realized coefficient family be invoked.  On that separated
family the six-column raw observation matrix has full AG rank; exact
largest-residual factor pivoting is continuous on each pivot chart and a
finite chart cover supplies a reader.  The already bounded nuisance/process
factors and finite operation count then make its complete 6x6 action
continuous.  Therefore the common ceiling is

`B_* = sup_W lambda_max(B_W) < infinity`.

This is the route from signed temporal separation to the historical action
ceiling.  The supremum is **not yet numerically or rigorously enclosed** in
this PR: the missing executable step is an outward enclosure of the two
signed margins over the admitted same-history window family, followed by the
finite pivot-chart action maximum.  Until that succeeds,
`uniform_historical_AG_readout_action` remains false and no J/rho0 claim is
promoted.  The existing 400--600 s positive finite margins remain diagnostics
only.

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
