# OU-III proof research state

This is the short current-state ledger required by the root `AGENTS.md`. Replace stale research history rather than accumulating it here.

## Current hypothesis

Canonical P3 is still open. P4 remains blocked by P3 and P5 remains blocked by P4. The canonical usefulness gate stays exactly `delta >= 1e-18`; the declared physical operating domain is unchanged.

PR #476 established that the current pseudo scheduler can erase already-earned elapsed time when tau changes its period:

`set_pseudo_update_period_s(): elapsed <- fmod(elapsed,new_period)`.

That is now the sharper limiter. Before lifting more upstream tuner state, test a progress-preserving period-retarget rule that never discards an owed S update. If the new period is already behind the retained elapsed time, park elapsed immediately below the new deadline so the next valid sample services the update; otherwise leave elapsed unchanged.

This is an implementation correction, not proof tuning. If adopted, all proof/evidence affected by the scheduler change must be rerun. No P3/P4/P5 promotion is permitted merely from the correction.

## Controlling inequality

For the deployed float scheduler, `h = 0.004999999888241291 s` and the largest tau-scaled pseudo period is about `T_S,max = 0.16363635 s`.

With progress preserved across period changes, the worst no-fire execution is simply the largest fixed period: period increases cannot remove elapsed credit, while a period decrease below the retained elapsed credit makes the next sample due. Exact binary32 `periodic_update_due<float>` reaches the largest period on sample 33 from zero elapsed.

Therefore the candidate source-independent liveness contract is

`N_S,max = 33 samples`,  i.e.  `T_gap,max <= 33 h ~= 0.165 s`.

This is the master quantity to validate before rebuilding any covariance certificate. It is stronger and simpler than carrying the complete `WavePeriodEstimator -> tuner -> tau` state merely to prevent timer-reset starvation.

## Evidence retained from #476

### Four path maxima are not a sufficient covariance-history quotient

The exact legal P2 path

`729 --16--> 568 --16--> 407 --18--> 246 --25--> 85 --13--> 74 --13--> 63`

reaches the full old adverse label `[9,3,39,9]` in 101 samples. Any covariance upper that keeps only those four path maxima loses the time order/residence duration needed by P3.

A non-promoting ordered 635-sample Riccati diagnostic reduced the old translation ceiling by roughly `1e9`, `1e9`, `1e10`, and `1e3` in variance across `[v,p,S,a_w]`. Ordered propagation therefore has ample numerical headroom once the scheduler contract is valid.

### The old fmod retarget rule admits exact starvation

The P2 clock language contains a legal gap-13 self edge. PR #476 certified an exact binary32 tau-EMA/timer execution inside one P2 source cell with:

- `tau_low = 9.533334732055664 s`;
- `tau_high = 9.533475875854492 s`;
- `T_S,low = 0.1300000101327896 s`;
- `T_S,high = 0.13000193238258362 s`;
- legal tuning frequencies `0.05247608572244644 / 0.05241831764578819 Hz`;
- zero S firings over 635 samples.

Current-head CI also repeated that result with the tau EMA exponential evaluated through float `expf`, so the counterexample is not an artifact of Python binary64 `math.exp`.

### Absolute-time source-clock aging is not the right repair

The dangerous source-clock gaps `13..16` arise only near binary64 uptime exponent 45 (about 1.11 million years), while normal-age spacing is 21 samples. That explains why deployed simulations never encounter the resonance, but silently imposing a deployment-lifetime bound would weaken the all-time theorem and is not accepted as the proof repair.

More importantly, nominal 21-sample commits do not by themselves create a source-independent liveness theorem while the setter may still erase elapsed progress. The timer retarget semantics, not another source partition, is the direct object to fix or certify.

## Failure classification / critic pass

**Failure type:** implementation-level hybrid scheduler liveness defect exposed by the proof.

**Invalidated approaches:**

- infer `S-gap <= cadence_max + h` while period changes use `fmod`;
- repair the defect only by finer tau/source partitions;
- replace the exact all-time source clock by nominal-age timing;
- lift `WavePeriodEstimator` state before testing whether the timer itself can be made source-independent.

**Current limiter:** establish and validate a progress-preserving period-retarget contract, then rebuild the ordered covariance upper against that actual shipping scheduler.

## Alternatives

1. **Progress-preserving scheduler retarget (selected next experiment).** Preserve elapsed credit when the period changes; if the new deadline has already passed, make the next sample due. Validate the 33-sample source-independent gap and rerun filter evidence.
2. **Full upstream source/timer reachability.** Keep the existing `fmod` behavior and lift `WavePeriodEstimator`, tuner, tau and timer state until starvation is excluded. This is much larger and remains the fallback if the scheduler correction materially harms filter behavior.
3. **Finite-initialization / multiword covariance theorem.** Allow zero-S words and prove boundedness from actual initialization plus later information. This is unattractive while the current timer admits an indefinite exact no-S cycle.

## Next falsifiable experiment

Implement only the period-retarget semantic correction and test three obligations:

1. elapsed below the new period is preserved exactly;
2. a shortened period that is already overdue causes an S update on the next valid sample rather than discarding elapsed time;
3. arbitrary legal period changes cannot extend the binary32 no-fire run beyond 33 samples, and the former 635-sample tau-EMA starvation witness is broken.

Then run the ordinary OU-III validation/evidence suite. If scored behavior changes materially, abandon or redesign the correction rather than tuning around it. If behavior remains acceptable, replace the obsolete starvation probes by a scheduler-progress certificate and resume the ordered whole-word P3 covariance construction.

## Retained facts

- Canonical P3 usefulness threshold is exactly `1e-18`.
- Physical P2 tuner partition remains 800 states; arbitrary Cartesian tau/sigma/R_S switching remains forbidden.
- Whole-word lower construction is still required; `P=0` remains a valid covariance lower start by Riccati monotonicity.
- H=18 and A=21 both remain required after translation P3 closes.
- No replay fitting, operating-domain shrink, gate tuning, or parameter tuning is allowed to make the proof pass.
- A deployed-filter change is acceptable only when it corrects an independently identified implementation defect and is followed by renewed evidence/proof validation; it cannot be credited as proof closure by itself.
- Lever arm remains disabled and the vibration guard remains on its dormant/transparent proof branch.
- P4 cannot promote before canonical P3; P5 cannot promote before strict canonical P4 contraction and must ultimately prove finite capture from the declared 45-degree entrance into the inner funnel.

## DEAD_ENDS / SHELVED

Do not resume these without a new mathematical fact:

- fixed-source 635-sample lower as a theorem proof;
- global whole-word lower against the old global four-max upper;
- 49-step gap-forgotten frontier;
- exact-elapsed/Pareto/minimum-cost enumeration while retaining the same four path maxima;
- independent Cartesian tau/sigma/R_S extrema;
- any covariance upper using `gap <= cadence_max + h` with the old `fmod` retarget rule;
- tau-cell subdivision or tau EMA alone as a proof of S recurrence;
- silently replacing all-time clock semantics by nominal deployment timing;
- blind subdivision, scalar-norm tightening, coefficient tuning, domain shrink, or gate tuning;
- additional P4 micro-certificates before P3 has a real margin.
