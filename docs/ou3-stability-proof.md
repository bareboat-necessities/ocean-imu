# OU-III stability theorem contract

## Scope

The shipping estimator is the object being analyzed. One persistent execution must simultaneously satisfy MARINE MOTION, IMU BIAS, and MAGNETIC SERVICE. Vessel motion, attitude, measurements, physical biases, frontend, tuner, covariance, acceptance state, scheduler, and estimator states belong to that same execution.

## MARINE MOTION

The wave coordinate `p` is displacement about a local equilibrium/reference, with `v=dp/dt` and `a=dv/dt`. The same history supplies attitude and angular rate. In addition to pointwise limits, admission requires

`|| integral_(t1)^(t2) p(tau) d tau || <= P_AC`

for every continuation. Equivalently, `qdot=p` has uniformly bounded potential differences. Therefore permanent `p(t)=p0 != 0` is outside the class, while quiet water `p=0` is inside. Bounded samples, finite means, spectra, or finite-window increments do not prove this all-time condition.

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

For each certified tail window rooted at `s`, let `Phi(k,s)` be the ordered differential of the complete preceding same-history shipping execution. Let `H_m,k` and `R_m,k` be the actual sensitivity and covariance of a correction that was actually applied. In normalized heading/axial-gyro-bias root coordinates,

`G_k=R_m,k^(-1/2) H_m,k Phi(k,s) E_hb`

or an algebraically equivalent inverse-free factorization. Service requires

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

## Loss of heading service

Indefinite absence of informative heading observations contains

`[theta+; e_bg+] = [[1,T],[0,1]] [theta; e_bg]`.

The block has spectral radius one and produces linear heading growth for nonzero axial gyro-bias error. It is a necessity result for MAGNETIC SERVICE, not a second stability theorem.

## Current closure

The machine-readable status in `reports/results/ou3_stability/theorem-status.json` intentionally remains open. The next controlling obligation is a finite-error magnetically informed H18 service-superword storage inequality with prefix retention.
