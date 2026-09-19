# OU-III proof research state

## Current hypothesis

The productive proof route does not require H18 asymptotic contraction. The
shipping wrapper holds accelerometer-bias learning while its magnetic reference
is provisional and releases that external hold after reference refinement. The
internal gate then clears after the accepted-update threshold and one-second
guard. H18 is therefore a finite bridge; the recurring asymptotic tail is A21.

The single route is
`construction -> finite capture -> finite H18 bridge -> finite reference
refinement/release -> recurring magnetically informed A21 -> regional practical
stability`.

## Negative result retained from PR #557

Strict full-state H18 incremental contraction in the covariance metric is
obstructed by the held accelerometer-bias identity block on a feasible held
segment. This invalidates that proof tactic only. The finite-difference export,
candidate-rho search, and associated measurement scripts were research
diagnostics, not proof machinery, and have been removed.

## Productive lemmas

`tail_stability.py` contains theorem algebra only. A finite H18 recurrence
`V+ <= g V + d` stays bounded for any finite number of bridge steps, including
`g >= 1`. Under MAGNETIC SERVICE, every service window contains at least one
actually applied informative correction; once magnetic-reference refinement
finishes at finite time, the remaining accepted-update count and one-second
guard therefore clear in finite time.

For the recurring A21 tail, if the source-uniform linear word satisfies
`||F e||_W <= sqrt(rho0)||e||_W` and the finite nonlinear remainder difference
has storage-norm Lipschitz gain `eta`, then the complete finite-error word has

`rho = (sqrt(rho0) + eta)^2`.

Thus the decisive small-gain condition is `sqrt(rho0)+eta < 1`. This is the
route to a real finite-error theorem; no sampled trajectory or candidate
numerical rho is used.

## Current limiter

Two implementation-linked results now control progress:

1. prove finite completion of the shipping magnetic-reference refinement. A
   sufficient lemma is now explicit: recurring usable samples close the literal
   min-sample/min-window, norm-ratio and horizontal-field gates in finite time.
   The deployed refinement has quality weighting disabled and hard-iron fitting
   disabled, so these are the controlling gates. The gate is now closed analytically on the captured domain rather than by a
   new assumption: rotation preserves the true field norm, so with true
   |B|>=20 uT and residual <=2 uT the running norm-ratio is at most
   4/(20-2)=0.222<0.35. With horizontal field >=15 uT, |B|<=75 uT and the same
   residual, the tuner's 5% horizontal gate is guaranteed whenever captured
   tilt-frame error is below about 7 degrees. MAGNETIC SERVICE then supplies at
   least one usable event per 1 s window, so the 128-sample/30-s window finishes
   within 128 s after refinement starts. The remaining release dependency is
   therefore capture into this <=7 degree tilt domain, not a separate magnetic
   norm-variation assumption. Then compose refinement with the internal
   accepted-update/guard logic to obtain finite H18 release;
2. establish a source-uniform **complete A21 information floor**. This is now
   the preferred linear route because it follows the Kalman covariance algebra
   instead of searching for a numerical rho. Prediction with Q>=0 is
   covariance-metric nonexpansive, and each literal Kalman correction is
   nonexpansive in its updated covariance metric. In normalized root
   coordinates a certified complete-word information floor mu>0 gives the
   comparison rho0<=1/(1+mu). The existing MAGNETIC SERVICE floor supplies the
   heading/axial-gyro-bias part; gravity/accelerometer and integral
   pseudo-measurements must supply the complementary directions. This **full
   information floor**, not magnetic service alone, is the next linear
   certificate. Then bound the nonlinear MEKF/reset/tuner remainder below
   1-sqrt(rho0).

   The A21 accelerometer-bias estimate projection need not enter that local
   nonlinear bound if the retained tail is kept inside its literal inactive
   region. With projection radius 0.4 m/s^2 and physical bias bound
   0.2251666 m/s^2, any bias-error norm below 0.1748334 m/s^2 keeps the estimate
   strictly inside the projection ball. Release/capture retention must prove
   entry into that inner domain.

Finite H18 bridge retention, the release operation itself, every-prefix A21
retention, physical qualification, capture and arithmetic closure remain open.

## Failed approaches / DEAD_ENDS

Do not reconstruct the PR #557 H18 full-state contraction search. H18 need not
be asymptotically contractive because it is a finite pre-tail bridge. Do not use
finite differences, a single simulated superword, or covariance-shaped point
ratios as theorem evidence. Do not narrow the physical motion/bias assumptions
to make a numerical contraction easier.

## Next proof step

Formalize the magnetic-reference refinement state machine and prove finite
release from recurring accepted informative corrections. In parallel, derive
the A21 linear error-transition energy inequality from the literal shipping
operations and the MAGNETIC SERVICE Gramian. The nonlinear small-gain remainder
is then bounded on the retained finite-error domain against that analytic
linear margin.


## A21 structural observability advance

For one scalar active-bias OU-III translation chain with state
`(v,p,S,a_w,b_a)`, exact constant-parameter OU transition, four successive
integral observations `H_S F^j`, `j=0..3`, and one attitude-normalized
accelerometer row `H_a=(0,0,0,1,1)`, a five-row observability minor has the
closed-form determinant

`det O = -Delta^3 tau^3 (1-exp(-Delta/tau))^3`.

It is nonzero for every finite `Delta,tau>0`. This removes a possible
structural rank obstruction: the integral chain plus active accelerometer-bias
state is observable per axis; the bias does not create an unobservable
translation mode after release.

This is not yet the shipping full-information certificate. The deployed tuner
can change `tau`, integral pseudo-updates occur on their actual scheduler, and
attitude is coupled rather than known. The next analytic step is to extend this
minor to bounded time-varying shipping transitions and combine its three-axis
floor with gravity/attitude and MAGNETIC SERVICE. That extension, rather than a
sampled matrix rank calculation, is the current route to the complete A21
information floor.


## Uniform time-varying A21 closure

The time-varying tuner/scheduler no longer requires freezing tau. The shipping
pseudo cadence is progress-preserving and satisfies
`T_S(tau)=clip((0.015/1.1)tau,0.005,0.15)`. On the admitted
`dt in [0.004,0.006]`, `tau in [0.02,12]` envelope, the accumulated OU
decay exponent between successive S updates is bounded by 0.55, including one
sample of scheduler overshoot. Across three gaps the latent-acceleration kernel
therefore retains at least `exp(-1.65)>0.19`.

For arbitrary positive time-varying `lambda(t)=1/tau(t)`, the S response to
initial acceleration has a kernel `K` with
`K'''(t)=exp(-integral lambda)>0`. Hence
`{1,t,t^2/2,K(t)}` is an extended complete Chebyshev system. Four successive
S observations have determinant equal to a Vandermonde factor times
`K'''(xi)/12`; the 4 ms sample-spacing floor makes this uniformly positive.
This proves time-varying translational observability for all three axes without
enumerating tuner histories.

The accelerometer attitude Jacobian has two nonzero singular values equal to
`|a_w-g|`, uniformly at least `9.80665-8.8=1.00665`. On the neutral
quotient (active accelerometer bias removed because its A21 OU predictor is
strictly stable), two separated gravity observations expose the two tilt/gyro
bias pairs. Known body rotation is orthogonal transport. The remaining
heading/axial-gyro-bias pair is exactly the subspace covered by MAGNETIC
SERVICE. Thus the non-decaying quotient is uniformly observable. With strict
inner guard margins, the finite hybrid branch cells form a compact union, so
the normalized quotient Gramian has an **existential uniform floor**

`mu_N := min_h lambda_min(J_N(h)) > 0`.

This is a proof of positivity, not yet a useful numerical enclosure of mu_N.

The A21 prediction is uniformly completely controllable on the same compact
envelope: the integrated-OU process covariance is SPD for positive dt, tau and
sigma; the exact attitude/gyro-bias covariance has positive gyro white-noise
and bias-RW densities; and active accelerometer bias has positive OU driving
density. The default periodic a_w covariance synchronization only queues a
bounded PSD increment, while bias release only raises finite P_ba diagonal
entries. These events preserve covariance compactness.

Uniform quotient observability/detectability plus uniform controllability gives
uniform upper/lower Riccati bounds and therefore a source-uniform linear A21
covariance-metric ratio `rho0<1`. No sampled rho is used.

On the strict inner domain (tilt error <=6 deg and accelerometer-bias error
<=0.15 m/s^2), projection is inactive and all quaternion measurement/reset
maps are smooth. The tuner and refined magnetic-reference paths are
measurement-only and therefore identical in a same-history comparison: they
belong to the bounded LTV schedule, not to the nonlinear remainder. Hence the
multiplicative nonlinear remainder satisfies `eta(r)->0` as `r->0`.
Since `rho0<1`, there exists `r_*>0` for which
`sqrt(rho0)+eta(r_*)<1`. Floating-point roundoff is retained as additive
supply rather than hidden in eta.

The remaining work is constructive rather than architectural: rigorously
enclose a usable numerical `mu_N`/`rho0`, derive an explicit `r_*`, bound
the additive arithmetic supply, and prove finite capture/H18 release enters and
retains that explicit inner domain.


## Marine motion supplies the missing attitude diversity

A separate attitude-excitation assumption is not required. Let
`f=a-g` be the world specific-force vector and `B` the true geomagnetic
field. On every interval of length `T`,

`integral f x B dt = (v(T)-v(0)) x B - T g x B`.

The marine contract gives `||v||<=5.5 m/s`; the magnetic contract gives
`B_h>=15 uT` and `||B||<=75 uT`. Therefore

`(1/T)||integral f x B dt|| >= g B_h - 2 V_max B_max/T`.

The right side is positive for `T>5.60844 s`. Choosing an 8 s subwindow gives
a uniform diversity floor `43.97475 (m/s^2) uT`. Hence every 8 s physical
marine history contains an instant where accelerometer and magnetic vector
sensitivities are non-collinear. MAGNETIC SERVICE supplies an applied magnetic
observation in every 1 s interval; transporting that sensitivity through the
known attitude transition to the diversity instant preserves its norm. The two
rank-two vector observations are therefore jointly full rank for attitude.
Two consecutive 8 s diversity windows expose gyro bias through its attitude
injection. The A21 proof word is consequently taken as 16 s.

This closes the earlier gap in the attitude/gyro-bias quotient argument using
existing MARINE MOTION and MAGNETIC SERVICE assumptions, rather than adding a
persistent-excitation assumption. The 16 s word also improves the active
accelerometer-bias homogeneous norm factor to `exp(-16/5000)=0.996805...`
(`rho_b=exp(-32/5000)<0.994`).
