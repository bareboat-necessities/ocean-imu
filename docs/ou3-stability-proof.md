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

## Finite H18 bridge and A21 asymptotic tail

PR #557 rules out one full-state H18 contraction tactic. It does not create a
new proof obligation: shipping H18 is a finite bridge while the magnetic
reference is provisional. A finite recurrence

`V_{k+1} <= g_H V_k + d_H`

has a finite bound for every finite bridge length, even when `g_H >= 1`.
Therefore H18 needs finite retention until the implemented reference refinement
and bias release, not asymptotic contraction.

MAGNETIC SERVICE supplies at least one actually applied informative correction
per service window. Conditional on finite completion of the outer magnetic
reference refinement, the remaining internal accepted-update count and the
one-second guard therefore clear in finite time. Proving finite reference
refinement from the shipping state machine is the current release obligation.

After release, A21 is the recurring tail. Let a source-uniform linear A21
service word satisfy

`||F e||_W <= sqrt(rho_0) ||e||_W`,  `rho_0 < 1`.

If the complete finite nonlinear shipping word differs from its linear word by
a same-history remainder with Lipschitz storage gain `eta`, then

`||F_nl(e_1)-F_nl(e_2)||_W <= (sqrt(rho_0)+eta)||e_1-e_2||_W`.

Hence the finite-error storage ratio is

`rho=(sqrt(rho_0)+eta)^2`,

and strict contraction follows from the explicit small-gain condition
`sqrt(rho_0)+eta<1`. Bounded physical/model/arithmetic forcing then gives the
usual practical-stability radius after coercivity and every-prefix retention
are certified.

For the linear A21 word, use the covariance/information structure directly.
Prediction with positive-semidefinite process covariance is nonexpansive in the
covariance metric, and literal linear Kalman corrections are nonexpansive in
their updated metric. Define the **complete normalized A21 information floor**
to include the transported information of all actual accelerometer/gravity,
integral pseudo-, and magnetic corrections over the service word. If this full
floor satisfies `J_A21 >= mu I`, `mu>0`, the information comparison gives

`rho_0 <= 1/(1+mu)`.

MAGNETIC SERVICE establishes only the heading/axial-gyro-bias component; it
must not be substituted for the complete floor. On the inner bias domain
`||e_ba|| < R_b-B_a = 0.1748334 m/s^2`, the shipping estimate projection is
inactive and therefore does not contribute a local A21 nonlinearity. Release
retention must place the execution in this domain.

This is the controlling proof route. Finite-difference superwords and candidate
point ratios are not proof evidence and are not retained.

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


## Uniform time-varying A21 detectability

The A21 proof does not freeze the adaptive operating point. The literal
pseudo-update cadence and the compact bounds
`dt in [0.004,0.006]`, `tau in [0.02,12]` imply a positive uniform
time-varying OU translation observability minor. The proof uses the extended
Chebyshev property of the four S-response functions, not a sampled rank test.

Gravity has a uniform tilt sensitivity floor `g-A_max=1.00665 m/s^2`.
After quotienting the strictly stable active accelerometer-bias OU mode, two
gravity observations expose tilt/transverse gyro bias; MAGNETIC SERVICE covers
heading/axial gyro bias. The neutral quotient information Gramian is continuous
on the finite union of strict compact hybrid branch cells and pointwise
positive definite. Therefore

`mu_N = min lambda_min(J_N) > 0`.

Positive process-noise densities give uniform complete controllability. The
Riccati covariance is consequently bounded above and below and the literal
linear A21 error word has some source-uniform `rho_0<1`. This is an
existence theorem; a constructive numerical enclosure of `mu_N` is still
needed for an explicit capture radius.

Inside the strict inner domain (tilt <=6 deg, `||e_ba||<=0.15 m/s^2`) the
bias projection is inactive. Tuner and refined magnetic-reference schedules are
same-history exogenous. The remaining MEKF measurement/reset maps are smooth,
so their multiplicative remainder gain satisfies `eta(r)->0`. Hence some
positive `r_*` satisfies `sqrt(rho_0)+eta(r_*)<1`. Floating-point error is
an additive bounded supply and does not consume this derivative margin.


## Marine-forced attitude diversity

The attitude information floor follows from the existing physical contracts.
For `f=a-g` and true field `B`,

`integral_0^T f x B dt = Delta v x B - T g x B`.

Thus
`T^-1 ||integral f x B dt|| >= g B_h,min - 2 V_max B_max/T`.
With the declared limits this becomes positive after 5.60844 s and is
43.97475 (m/s^2) uT for T=8 s. Every 8 s marine history therefore contains
non-collinear gravity/specific-force and magnetic sensitivity. Recurring
MAGNETIC SERVICE transports an applied magnetic sensitivity into that interval.
Two 8 s subwindows form the 16 s A21 proof word and expose gyro bias through
attitude propagation. No additional attitude-excitation assumption is used.


## Quantitative enclosure status

A first fully analytic normalized translation enclosure is now constructive.
On the 16 s word, selecting S updates after times 0, 8 and 16 s with the
literal maximum scheduler delay 0.156 s, state scales
(V,p,S)=(5.5,8.1,1100), and worst declared S-noise standard deviation 100,
the determinant/Frobenius certificate gives

`mu_trans >= 2.04734e-3`.

This is a real lower bound, not a sampled singular value.

A deliberately sparse two-epoch attitude/gyro calculation gives a much weaker
candidate scale (~6.18e-7) and therefore is **not promoted as the final
shipping mu_N certificate**. The reason is important: MAGNETIC SERVICE is an
already-transported two-coordinate heading/axial-bias Gramian over a window,
not an instantaneous pure-heading row. A tight full certificate must compose
that actual 2-D service Gramian directly with the many accelerometer rows in
the same 16 s word; replacing it by a fictitious instantaneous attitude
measurement would be an invalid shortcut. The next quantitative proof step is
therefore the aggregate Schur/Gramian bound on the literal partition, using all
recurring accelerometer and S information rather than two sparse rows.

The float32 arithmetic path is likewise separated correctly. A straight-line
kernel with n rounded operations has the standard gamma_n bound
`gamma_n=n*u/(1-n*u)`, u=2^-24. This is implemented as a certificate
primitive, but no whole-word arithmetic supply is claimed until literal kernel
operation counts and magnitude envelopes are composed. Roundoff remains
additive supply and is not allowed to consume the nonlinear derivative margin.
