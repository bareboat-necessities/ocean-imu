# OU-III proof research state

## Current hypothesis and controlling inequality

The published finite-window COMPLETE-BRMM primitive qualification is insufficient
for a finite indefinite bound on all 18 motion errors because it bounds only
short-window `Delta S`, not the indefinite displacement primitive. The exact
quiet-source ambiguity certificate remains valid. Canonical P3 is unchanged at
`delta=1e-18`; P4/P5 remain false.

The positive continuation is now an origin-invariant recurrence qualification,
not an independent fresh-entry S ball. For the real Live handoff time `t_L`, keep
PR #513's one-time coordinate `S_L(t)=S(t)-S(t_L)`. A physical sufficient source
condition is a uniform finite constant `D_S` satisfying

`||integral_u^t p(s) ds|| <= D_S` for every admitted history and all `t,u>=t_L`.

This is equivalent to a uniform bound on the diameter of each physical S
trajectory after removing only its constant origin. It preserves `e_S,L(t_L)=0`,
does not re-anchor position, and does not re-zero S at later words.

## Exact failure retained

`ou3_brmm_infinite_continuation.py` materializes the indistinguishable quiet pair
`p=+d` and `p=-d`, `v=a=omega=0`, fixed attitude and identical sensors. The same
shipping estimate then gives

`e_S,L^+ - e_S,L^- = 2 h d`,

hence `max(||e_S,L^+||,||e_S,L^-||) >= h||d||`. Uniform coercivity gives the
corresponding quadratic storage lower bound. This is classification **B** for the
finite-window-only indefinite target. If the intended source includes a bounded
indefinite S primitive, the missing numerical/source qualification is **E**.
No nominal filter instability, P3 failure or fresh-entry 300 m*s counterexample is
claimed.

## New materialized source condition

`ou3_brmm_centered_S_recurrence.py` now encodes the missing recurrence coordinate
and consumes the exact obstruction. It explicitly rejects:

* an absolute S-origin bound;
* the legacy 300 m*s value as a fresh-entry premise;
* position re-anchoring;
* wordwise S re-zeroing;
* promotion from bounded position or bounded 3 s `Delta S` alone.

The stronger source remains sign-symmetric. Therefore the existing 300 m*s
**working/retention** radius supplies a necessary search ceiling, not a source
assumption: a symmetric indistinguishable pair can separate by `2 D_S`, whereas
two errors inside a radius-300 S tube can separate by at most 600. Consequently
the positive retention search is restricted to

`0 <= D_S <= 300 m*s`.

The true admissible maximum can be smaller once Joseph/reset, compatible-storage,
every-prefix and first-exit margins are included. `D_S_max` therefore remains
unfrozen and `COMPLETE_BRMM_INDEFINITE_S_QUALIFIED=false`.

## Retained positive facts

Fresh v/p/S/aw/bg/ba estimator means are held zero on the fresh wrapper path; the
common S origin is removed exactly once. BIAS0/1/2 retain independent physical
driver recurrences and bounded-bias projection proofs; positive BIAS2 separation
is optional. The qualified runtime Live/H18 handoff remains at most 150 s and is
not P4 capture. Conditional binary32 mathematics remains separate from deployment.

The joint 24-state compatible-storage point candidate remains the useful route:
H18 product rate about 0.9996524356 and A21 about 0.9959531012. The A21 18-state
marginal covariance storage and the full-product 3 s nonlinear route remain dead
ends; neither invalidates the joint error/physical-bias architecture.

## Current limiter

The next controlling quantity is no longer whether an indefinite S premise is
necessary; that is settled. It is the largest `D_S` that the actual same-history
joint24 endpoint/prefix/first-exit construction can retain while preserving all
other declared source coordinates and literal shipping events.

The source transition must carry, on one history, physical `p/S_L`, frontend/WPE,
raw and effective sigma, tau, T_S, anisotropic R_S, scheduler/commit state,
reachable P/H/R/K, Joseph correction, finite reset, BIAS0/1/2 projection and the
successor storage cell. Independent per-word `Delta S` ports are forbidden.

## Next falsifiable experiment

Parameterize the existing joint24 reachable-prefix construction by one
source-history centered-S recurrence coordinate `D_S`, then search outward over
`[0,300] m*s`. For every candidate require the same source cell to close:

1. endpoint compatible-storage contraction including metric-transition penalty;
2. every literal prefix augmented positivity certificate;
3. correction/reset chart validity;
4. hard first-exit retention;
5. all BIAS0/1/2 driver and projection branches; and
6. conditional binary32 additive ISS enclosure.

Report the first limiting event/source lineage and the largest rigorously retained
`D_S`. Only that result may freeze the COMPLETE-BRMM indefinite-S qualification.
Afterward resume full source-uniform cover, P4 basin maximization, H18 capture and
H18->A21 transport. No endpoint margin, maximum P4 basin or capture time is yet
certified.
