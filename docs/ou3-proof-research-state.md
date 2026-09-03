# OU-III proof research state

This is the short current-state ledger required by the root `AGENTS.md`. Replace stale research history rather than accumulating it here.

## Current hypothesis

Canonical P3 is still open. P4 remains blocked by P3 and P5 remains blocked by P4. The canonical usefulness gate stays exactly `delta >= 1e-18`; the deployed filter and declared operating domain are unchanged.

The old P3 architecture is now rejected in two independent ways:

1. its covariance floor and ceiling lived on different horizons; and
2. its covariance ceiling discarded source order into four path maxima and assumed a finite S=0 pseudo-observation gap that the current P2->P3 source quotient does not prove.

The replacement P3 must use a whole-word covariance comparison and must preserve enough tuner/scheduler memory to justify whatever S-information recurrence it uses.

## Evidence

### Ordered covariance mechanism is large enough

The exact legal P2 path

`729 --16--> 568 --16--> 407 --18--> 246 --25--> 85 --13--> 74 --13--> 63`

reaches the full old four-max adverse label `[9,3,39,9]` in 101 samples. Therefore exact elapsed/Pareto/frontier enumeration cannot repair an upper that still feeds those four independent maxima to the old formula.

A non-promoting 635-sample point Riccati diagnostic that retains source order and the shipping pseudo-timer semantics starts from the same old covariance ceiling and obtains standard deviations

- old four-max ceiling: `[3.3454e4, 6.3004e4, 8.4732e4, 3.3302e1]`;
- ordered legal history: `[1.25985, 1.17849, 0.55534, 0.85421]`.

The variance reductions are about `7.05e8`, `2.86e9`, `2.33e10`, and `1.52e3`. This is feasibility evidence only, but it confirms that time ordering and repeated S updates are the missing mechanism rather than a small enclosure improvement.

### Existing finite S-gap premise is false for the current proof language

Shipping `set_pseudo_update_period_s()` reduces the retained pseudo elapsed time with `fmod(elapsed,new_period)`. The tuner commits the current smoothed tau on its finite source clock, and every commit reapplies the tau-scaled pseudo period.

The P2 clock certificate admits exact stage gaps `13..26` samples over the full binary64 lifetime. Node 720 has an exact gap-13 self edge. Inside that one P2 source cell,

- tau cell: `[8.3859254253, 12.0] s`;
- pseudo-period cell: `[0.1143535210, 0.1636363529] s`;
- `tau_low = 9.5333347321 s` -> `T_S = 0.1300000101 s`;
- `tau_high = 11.0 s` -> `T_S = 0.1499999911 s`.

The exact binary32 scheduler transcription has a repeatable `H,H,L,H,L,H,...` reset cycle with **zero S firings over all 635 samples**. CI run `33815085697` validates the exact gap-13 self edge, both tau/period values inside the same P2 cell, the repeatable no-fire cycle, and zero firings.

Therefore the old upper premise

`S observation gap <= max(pseudo period) + h`

is **not a theorem of `OU3_P2_CORRELATED_STAGE_TRANSFER_V1`**. Any P3 upper that uses that premise is invalid, irrespective of how tightly the covariance algebra is enclosed.

This result does **not** yet prove that the full deployed WavePeriodEstimator + tau EMA can realize the no-fire cycle. The current P2 quotient deliberately forgets within-cell tau evolution and upstream estimator memory. That distinction is the next research question.

## Failure classification / critic pass

**Failure type:** proof-method / source-abstraction failure.

**Invalidated strategy:** infer a finite S-packet gap from the maximum tau-coupled pseudo period while allowing source changes represented only by the current coarse P2 cell/pair state.

**Not invalidated:** the deployed filter, the declared physical operating envelope, the whole-word covariance feasibility signal, covariance monotonicity, or the possibility of proving P3 with a stronger source/scheduler representation.

**Current limiter:** guaranteed translation observability/information in the presence of source-dependent pseudo-period changes and retained timer phase.

Strongest reason to abandon the current architecture: the missing variable is not another covariance norm; it is an actual hybrid state (`pseudo_update_elapsed_s_`) whose reset map depends on the continuously evolving tau. Forgetting that state can manufacture a legal proof word with no S observations at all.

## Alternatives

1. **Shipping source + timer reachability.** Lift the P2 interface with the tau-EMA state (or a certified sufficient quotient) and pseudo timer phase, then propagate cumulative S information/covariance in time order. This keeps the theorem source complete and directly attacks the missing state.
2. **Finite-initialization / multiword covariance theorem.** Stop requiring every covariance word to wash out arbitrary initial covariance. Start from the actual certified finite filter initialization and prove bounded covariance using cumulative information over multiple words, allowing individual zero-S words. This changes the covariance theorem architecture, not the filter or domain.
3. **Derive estimator regularity from shipping code.** Prove a rate/dwell/total-variation property for the WavePeriodEstimator -> SeaStateAutoTuner -> tau EMA chain strong enough to exclude timer-reset starvation, then reuse a simpler scheduler quotient. This is acceptable only if derived from implementation and the existing physical source envelope; adding an artificial regularity assumption to make P3 pass is forbidden.
4. **Information-form recurrence instead of packet-gap selection.** Bound the complete ordered information contribution of whatever S packets actually occur rather than selecting three observations separated by a prescribed gap. This may combine naturally with either (1) or (2) and avoids making `max period + h` the controlling lemma.

## Next falsifiable experiment

**Test whether the shipping tau EMA itself can realize the scheduler reset cycle under the presently admitted tuning-frequency input family.**

Use the exact binary32 tau update, dynamic EMA horizon, 13-sample certified source clock, tau->pseudo-period map, `fmod` period setter, and `periodic_update_due` semantics. Search for legal tuning frequencies whose 13-sample EMA images alternate between two applied tau values around the critical reset period while producing no S firing from a post-fire timer state.

Interpretation:

- **cycle realizable:** coarse tau partitioning is not the root cause. The proof needs upstream WavePeriodEstimator/source regularity or a scheduler-aware theorem that tolerates such cycles; do not build a finer 800-state-style partition and call it closure.
- **cycle not realizable with tau EMA:** refine the P2/P3 interface to retain the tau-EMA relation and timer phase; the current cell quotient is the sole starvation source.

This experiment is diagnostic only. It cannot promote P3/P4/P5.

## Retained facts

- Canonical P3 usefulness threshold is exactly `1e-18`.
- `OU3_P2_CORRELATED_STAGE_TRANSFER_V1` remains the current source authority, but it is now known to be insufficient for any P3 lemma requiring a finite S firing gap.
- Physical P2 partition remains 800 states; arbitrary Cartesian tau/sigma/R_S switching remains forbidden.
- Whole-word lower construction is still required; `P=0` is a valid covariance lower start by Riccati monotonicity.
- H=18 and A=21 both remain required after translation P3 closes.
- No replay fitting, operating-domain shrink, gate tuning, or deployed-filter change is allowed for proof convenience.
- Lever arm remains disabled and the vibration guard remains on its dormant/transparent proof branch.
- P4 cannot promote before canonical P3; P5 cannot promote before strict canonical P4 contraction and must ultimately prove finite capture from the declared 45-degree entrance into the inner funnel.

## DEAD_ENDS / SHELVED

Do not resume these without a new mathematical fact:

- fixed-source 635-sample lower as a theorem proof;
- global whole-word lower against the old global four-max upper;
- 49-step gap-forgotten frontier;
- exact-elapsed/Pareto/minimum-cost enumeration while retaining the same four path maxima;
- independent Cartesian tau/sigma/R_S extrema;
- any covariance upper using `gap <= cadence_max + h` directly from the present P2 quotient;
- blind subdivision, scalar-norm tightening, or deeper search as a substitute for source/timer memory;
- sigma/R_S coefficient, filter, domain, or canonical-gate tuning;
- additional P4 micro-certificates before P3 has a real margin.
