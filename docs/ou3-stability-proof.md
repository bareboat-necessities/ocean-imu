# OU-III stability theorem contract

## Scope

The shipping estimator is the object being analyzed. One persistent execution must simultaneously satisfy MARINE MOTION, IMU BIAS, and MAGNETIC SERVICE. Vessel motion, attitude, measurements, physical biases, frontend, tuner, covariance, acceptance state, scheduler, and estimator states belong to that same execution.

## MARINE MOTION

The wave coordinate `p` is displacement about a local equilibrium/reference, with `v=dp/dt` and `a=dv/dt`. The same history supplies attitude and angular rate. In addition to pointwise limits, admission requires

`|| integral_(t1)^(t2) p(tau) d tau || <= P_AC`

for every continuation. Equivalently, `qdot=p` has uniformly bounded potential differences. The same potential continuation persists through successive proof superwords; a one-time proof-coordinate origin at certified-tail entry does not authorize later physical resets or reanchoring. Therefore permanent `p(t)=p0 != 0` is outside the class, while quiet water `p=0` is inside. Bounded samples, finite means, spectra, or finite-window increments do not prove this all-time condition.

Physical translation is `p_CoG=p_eq+p`. Global origin, current, propulsion, leeway, and secular reference motion may be represented in `p_eq`, but nonzero `ddot(p_eq)` remains actual specific force/model disturbance unless the shipping implementation explicitly compensates it.

## IMU BIAS

One deterministic physical assumption covers total residual accelerometer and gyro bias after actual calibration:

`||b_a||<=B_a, ||dot b_a||<=D_a`

`||b_g||<=B_g, ||dot b_g||<=D_g`.

The sampled physical recurrence is predecessor constrained:

`b_a[k+1]=b_a[k]+w_a[k], ||w_a[k]||<=D_a dt[k]`

and analogously for gyro bias.

For accelerometer-bias prediction the proof uses one first-class relation in both shipping modes:

`e_b[k+1]- = phi_e[k] e_b[k]+ + (1-phi_e[k]) b[k] + w[k]`,

with `phi_e[k]=1` in H18/held prediction and `phi_e[k]=phi_OU[k]` in A21/active prediction. This is an estimator-mode coefficient, not a physical-bias law.

Correction is carried separately using the actual estimator increment:

`e_corr=e_minus-delta_bhat`.

Projection is then carried separately with the shipping Euclidean estimate projection:

`bhat_plus=Proj_R(bhat_corr)`, `e_plus=b_true-Proj_R(bhat_corr)`.

Estimator hold, release, correction, and projection never reset physical truth. The shipping gyro-bias mean predictor is identity, so its prediction error is `e_g_minus=e_g_plus+w_g`.

## MAGNETIC SERVICE

For each certified tail window rooted at `s`, let `Phi(k,s)` be the ordered differential of the complete preceding same-history shipping execution. The transport ends immediately before the magnetic correction, not after it. Let `H_m,k` be the literal pre-correction magnetic sensitivity and `S_m,k^act` the actual innovation covariance presented to the shipping factorization for a correction that was actually applied. In normalized heading/axial-gyro-bias root coordinates,

`G_k=W_k H_m,k Phi(k,s) E_hb`, with `W_k^T W_k=(S_m,k^act)^(-1)`.

Equivalently, if `S_m,k^act=L_k L_k^T`, use `W_k=L_k^(-1)`. An inverse-free LDLT form is equally valid. The essential point is that the same actual innovation covariance factored by the shipping update is used. Service requires

`sum G_k^T G_k >= mu_M I_2`

on every certified interval of length `T_M`.

Attempted callbacks, due events, packets, rejected/invalid measurements, saturation, and maximum gap alone do not establish service.

## One proof path

`construction -> startup/capture -> magnetically informed Live/H18 -> H18-to-A21 release -> magnetically informed A21 -> regional practical stability`.

Capture may be history dependent: `T_c=T_c(h,x0)<infinity`. No common startup deadline is assumed. If Live begins before sufficient magnetic information is established, that interval remains in the finite pre-certified prefix.

## Tail target

At informative-service superword boundaries, a completion certificate must establish coercivity, every-prefix retention, and a finite-error inequality

`V_(j+1) <= rho V_j + c_d ||d||^2_[j,j+1]`, with `rho<1`.

The complete relation must include prediction, accelerometer correction, integral pseudo-update, magnetic correction, bias correction and projection, quaternion reset, covariance evolution, tuner commits, scheduler/clocks, magnetic-reference changes, tilt relocks, and the H18-to-A21 release.

On the H18 leg this target is not available for the full error coordinate; see the next section. The H18 obligation is therefore the restated one: strict dissipation on the complement of the held accelerometer-bias coordinates, with the held bias carried as a bounded input.

## Loss of heading service

Indefinite absence of informative heading observations contains

`[theta+; e_bg+] = [[1,T],[0,1]] [theta; e_bg]`.

The block has spectral radius one and produces linear heading growth for nonzero axial gyro-bias error. It is a necessity result for MAGNETIC SERVICE, not a second stability theorem.

## Held accelerometer bias in H18

On a held segment with completed finite operations, an already feasible bias
estimate, and no frame/relock or projection-changing event, the estimator bias
is constant, its gain rows vanish, and its covariance block is frozen and
decoupled. Physical bias need not be constant: the absolute error obeys
`e_b[k+1]=e_b[k]+w[k]`. Only differences between executions with the **same
physical history** cancel this shared increment.

The same-history incremental differential has block form
`Psi=[[A,B],[0,I3]]`. It is not a global linear formula for the absolute
nonlinear finite error. With decoupled covariance and unchanged `P_bb`, a
direction `delta e=(0,u)` gives

`V_end(Psi delta e)=(B u)^T P_oo,end^-1 (B u)+V_root(delta e)>=V_root(delta e)`.

The identity block also precludes strict contraction in a fixed quadratic
metric. This is a held-mode necessity result, not instability, and does not
remove the separate projection, relock, frame-change or release obligations.
The absolute held bias is bounded by `|e_b|<=B_a+R_b` once the finite completed
projection/feasible-estimate premises hold. The single proof path therefore
uses the complement with held bias as input on H18, then the actual release
and A21 dynamics; it does not bypass the held leg.

## H18 complement and supply inequality

Write the exact finite-error relation on one retained superword as

`x_plus=A x+B e_b+r`, `V=x^T W x`, `W=P_oo^-1`.

This identity can define `r` for a chosen surrogate `A,B`; it does not bound
`r`. In particular, central differences about a nonzero reached error measure
local incremental sensitivity. They do not eliminate the reference forcing,
nonlinear terms, changing covariance/gains, frontend/tuner dependence, physical
bias increments, model mismatch or arithmetic. All remain in `r` unless
explicitly represented and bounded. Later superwords inherit the actual
auxiliary state rather than the fixed auxiliary root of this diagnostic.

**Conditional finite-error supply lemma.** If certified matrices satisfy
`A^T W_plus A <= alpha W` with `0<=alpha<rho<1`, then

`V_plus <= rho V + rho/(rho-alpha) |B e_b+r|^2_W_plus`.

Proof: use `|A x|_W_plus<=sqrt(alpha V)` and apply weighted Young's inequality
to the cross term with weight `(rho-alpha)/alpha`; the case `alpha=0` is
immediate. This is the location of the held-bias input in the controlling tail
inequality. It is conditional until the same-history matrix bound, finite-error
supply bound, coercivity and retained domain are certified.

For the isolated linear bias contribution, let `P_oo=L L^T` and
`P_oo,plus=L_plus L_plus^T`, `M=L_plus^-1 A L`, `N=L_plus^-1 B`.
If `Q=rho I-M^T M` is positive definite, completing the square gives the sharp
linear gain

`gamma_b=lambda_max(N^T N+(M^T N)^T Q^-1(M^T N))`.

It bounds `|M y+N u|^2<=rho |y|^2+gamma_b |u|^2` for the linear surrogate
only. It must not replace the unknown total finite-error supply. Every-prefix
retention additionally needs all-direction operator and supply bounds at each
prefix, not the trajectory of the endpoint-maximizing direction alone.

## Current closure

The machine-readable theorem status remains open. The complement diagnostic
finds a positive but small endpoint margin. Its fine/coarse operator
discrepancy is larger than the available norm margin; that discrepancy is not
a rigorous uncertainty bound. Neither the exact nonlinear matrix inequality
nor a bounded finite-error supply has been established. The next experiment
must resolve sensitivity conditioning and bound the same-history finite-error
remainder before any rigorous stability enclosure is claimed.

## Finite-error operation lemmas

In the complete finite-error storage inequality the bias prediction is a joint
map on `(e_b,b_true)`, with scalar blocks `[[phi_e,1-phi_e],[0,1]]` tensored
with `I3`. The physical increment has one shared column `[I3;I3]`. The held
and active modes differ only through the estimator coefficient; gyro-bias
self-prediction remains identity. Local absolute continuity gives the physical
sampled rate bound; an a.e. derivative bound without that regularity is not
sufficient to rule out jumps.

A correction uses its actual estimate increment. Projection contributes the
separate defect `d=bhat_corr-project(bhat_corr)`, so `e_plus=e_corr+d`.
For positive radius and finite completed projection, `|e_plus|<=B_a+R_b`.
When the true bias is inside the projection ball,
`|e_plus|^2+|d|^2<=|e_corr|^2`. This is a Euclidean component lemma, not a
full-state weighted-metric or covariance contraction result. Nonpositive
radius disables projection. The shipping invalid-attitude-injection return
bypasses projection and remains outside these finite-domain lemmas.

The persistent marine potential gives `S_true=q(t)-q(T_c)` with one fixed
origin, hence the actual integral innovation is `-S_hat=e_S-S_true`. The
physical forcing term cannot be dropped when composing the Kalman corrections.
Magnetic sensitivity transport accepts rectangular maps at coordinate changes;
this does not replace the still-open source-qualified release map.

Algebraic lemmas and native regression tests do not establish target arithmetic,
source-uniform finite-error dissipativity, capture, or retention.
