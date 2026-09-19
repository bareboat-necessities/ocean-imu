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

1. prove finite completion of the shipping magnetic-reference refinement under
   MAGNETIC SERVICE and the admitted sensor domain, then compose it with the
   internal accepted-update/guard logic to obtain finite H18 release;
2. establish a source-uniform linear A21 contraction margin from the actual
   prediction/correction sequence and MAGNETIC SERVICE, then bound the nonlinear
   MEKF/reset/projection/tuner remainder below that margin.

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
