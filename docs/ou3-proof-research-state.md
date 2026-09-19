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
