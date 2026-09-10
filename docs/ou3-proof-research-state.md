# OU-III proof research state

## Current hypothesis and controlling inequality

The published finite-window COMPLETE-BRMM primitive qualification is insufficient
for a finite indefinite bound on all 18 motion errors because it bounds only
short-window `Delta S`, not the indefinite displacement primitive. The exact
quiet-source ambiguity certificate remains valid. Canonical P3 is unchanged at
`delta=1e-18`; P4/P5 remain false.

The positive continuation is an origin-invariant recurrence qualification, not an
independent fresh-entry S ball. For the real Live handoff time `t_L`, keep PR
#513's one-time coordinate `S_L(t)=S(t)-S(t_L)`. A physical sufficient source
condition is a uniform finite constant `D_S` satisfying

`||integral_u^t p(s) ds|| <= D_S` for every admitted history and all `t,u>=t_L`.

This preserves `e_S,L(t_L)=0`, does not re-anchor position, and does not re-zero S
at later words.

## Exact failure retained

`ou3_brmm_infinite_continuation.py` materializes the indistinguishable quiet pair
`p=+d` and `p=-d`, `v=a=omega=0`, fixed attitude and identical sensors. The same
shipping estimate gives

`e_S,L^+ - e_S,L^- = 2 h d`,

so `max(||e_S,L^+||,||e_S,L^-||) >= h||d||`. Uniform coercivity gives the
corresponding quadratic storage lower bound. This is classification **B** for the
finite-window-only indefinite target. If the intended source includes a bounded
indefinite S primitive, the missing numerical/source qualification is **E**.
No nominal filter instability, P3 failure or fresh-entry 300 m*s counterexample is
claimed.

## Materialized strengthened-source facts

`ou3_brmm_centered_S_recurrence.py` encodes the missing recurrence coordinate and
consumes the exact obstruction. It rejects an absolute S-origin bound, the legacy
300 m*s value as a fresh-entry premise, position re-anchoring, wordwise S
re-zeroing, and promotion from bounded position or bounded 3 s `Delta S` alone.

For the sign-symmetric strengthened source, the existing 300 m*s
**working/retention** radius gives a necessary search ceiling, not a source
assumption. Hence the positive search is `0 <= D_S <= 300 m*s`; `D_S_max` remains
unfrozen until the full same-history retention calculation closes.

`ou3_brmm_translation_observation_kernel.py` closes the exact integration-constant
kernel conditional on finite `D_S`. For two identical-acceleration histories,
per axis

`delta v=c_v`, `delta p=c_p+h c_v`,
`delta S_L=h c_p + h^2 c_v/2`.

Uniform bounded position forces `c_v=0`; finite centered-S recurrence then forces
`c_p=0`. Thus the two-dimensional translational zero-output kernel per axis drops
to one dimension under the position primitive and to zero under the centered-S
primitive. This uses no numerical value of `D_S`, no covariance membership and no
absolute S-origin bound. It removes the exact ambiguity mechanism but does not
prove contraction for nonzero source variations.

## Point feasibility after the repair

The non-promoting same-history joint24 diagnostic on the unmodified captured
shipping graph is strongly feasible in the S coordinate. At one correlated
covariance-root sigma it reports minimum extra Euclidean integral-displacement
headroom:

* H18: `299.2550971262734 m*s`, limiting at sample 6603 prediction;
* A21: `299.51436721585696 m*s`, limiting at sample 228612 prediction.

The point S retention ratios are only about `0.00248301` H18 and `0.00161878`
A21. The actual correlated-root limiter is instead H18 velocity, with critical
initial level `4.4551206340521015 sigma`; A21 is limited by latent acceleration at
`112.42280878895119 sigma`. The compatible-storage product rates remain
H18=`0.9996524355874086` and A21=`0.9959531012474646`.

These numbers justify continuing the joint24 route, but they do **not** qualify
`D_S`: the captured history has zero physical-bias driver and is not the universal
COMPLETE-BRMM family. An independent additive S port is still forbidden.

## Retained positive facts

Fresh v/p/S/aw/bg/ba estimator means are held zero on the fresh wrapper path; the
common S origin is removed exactly once. BIAS0/1/2 retain independent physical
driver recurrences and bounded-bias projection proofs; positive BIAS2 separation
is optional. The qualified runtime Live/H18 handoff remains at most 150 s and is
not P4 capture. Conditional binary32 mathematics remains separate from deployment.

The joint 24-state compatible-storage route remains preferred. The A21 18-state
marginal covariance storage and the full-product 3 s nonlinear route remain dead
ends; neither invalidates the joint error/physical-bias architecture.

## Current limiter

The limiting construction is now the **same-history nonzero-source continuation**,
not S point capacity. It must carry one physical `p/S_L` history and the same
acceleration/rotation sample into frontend/WPE, raw/effective sigma, tau, T_S,
anisotropic R_S, scheduler/commit state, reachable P/H/R/K, Joseph correction,
finite reset, BIAS0/1/2 driver/projection and successor compatible-storage cell.

The existing 601-sample typed execution kernel is downstream-ready, but the hard
source provider still lacks a validated correlated COMPLETE-BRMM window oracle.
For the primary BRMM theorem the indefinite centered-S recurrence also needs a
numerical/source qualification; for the retained directional-sea SEA0 extension
the hard spectral-driver/output oracle remains separately open. Neither may be
replaced by replay or independent sample boxes.

## Next falsifiable experiment

Construct the origin-invariant centered-S history as a retained coordinate of the
hard source transition, not as a per-word supply port. Parameterize the same-history
joint24 endpoint/prefix graph by `D_S` and search outward over `[0,300] m*s`.
For every candidate source cell require:

1. endpoint compatible-storage contraction including metric-transition penalty;
2. every literal prefix augmented positivity certificate;
3. correction/reset chart validity;
4. hard first-exit retention;
5. all BIAS0/1/2 driver and projection branches; and
6. conditional binary32 additive ISS enclosure.

Report the first limiting event/source lineage and the largest rigorously retained
`D_S`. Only that result may freeze the strengthened COMPLETE-BRMM indefinite-S
qualification. Afterward resume full source-uniform cover, P4 basin maximization,
finite H18 capture and H18->A21 transport. No universal endpoint margin, maximum
P4 basin or capture time is yet certified.
