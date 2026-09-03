# OU-III proof research state

This is the short current-state ledger required by the root `AGENTS.md`. Replace stale research history rather than accumulating it here.

## Current hypothesis

Canonical P3 is still open. P4 remains blocked by P3 and P5 remains blocked by P4. The canonical usefulness gate stays exactly `delta >= 1e-18`; the deployed filter and declared operating domain are unchanged.

The old P3 architecture is rejected for three independent reasons:

1. its covariance floor and ceiling lived on different horizons;
2. its covariance ceiling collapsed source order into four path maxima; and
3. its finite S=0 observation-gap lemma is not implied by the shipping source/timer dynamics retained so far.

A replacement P3 must be whole-word and time ordered. Before constructing another covariance upper, it must establish a valid source-complete recurrence of S information, or explicitly prove a theorem architecture that remains useful without such recurrence.

## Evidence

### Four path maxima are not a sufficient history quotient

The exact legal P2 path

`729 --16--> 568 --16--> 407 --18--> 246 --25--> 85 --13--> 74 --13--> 63`

reaches the full old adverse label `[9,3,39,9]` in 101 samples. Therefore exact elapsed/Pareto/frontier enumeration cannot repair an upper that still feeds those four independent maxima to the old formula. Time order/residence duration is the missing information.

A non-promoting 635-sample point Riccati diagnostic retaining the legal source order and pseudo scheduler reduces the old translation ceiling by roughly `1e9`, `1e9`, `1e10`, and `1e3` in variance across `[v,p,S,a_w]`. This is feasibility evidence only; it shows that ordered propagation has ample numerical headroom if a valid scheduler/source theorem can be supplied.

### The current P2 quotient admits complete S starvation

Shipping `set_pseudo_update_period_s()` performs `elapsed = fmod(elapsed,new_period)` before the current sample's pseudo-update test. The P2 clock certificate admits exact stage gaps `13..26` samples over the full binary64 lifetime.

CI run `33815085697` certifies that P2 node 720 has an exact gap-13 self edge and that the current P2 quotient admits a repeatable 635-sample word with zero S firings. Hence

`S observation gap <= max(pseudo period) + h`

is not a theorem of `OU3_P2_CORRELATED_STAGE_TRANSFER_V1`.

### The starvation resonance survives the shipping tau EMA

CI run `33818178464`, artifact `ou3-p3-tau-ema-scheduler-cycle-diagnostic` (`9917480639`), validates an exact binary32 tau-EMA/timer cycle with:

- source node `720`;
- exact source-clock gap `13` samples;
- `tau_low = 9.533334732055664 s`;
- `tau_high = 9.533475875854492 s`;
- tau separation only `1.41143798828125e-4 s`;
- `T_S,low = 0.1300000101327896 s`;
- `T_S,high = 0.13000193238258362 s`;
- period separation `1.9222497940063477e-6 s`;
- scheduler tolerance `1.9073486328125e-6 s`;
- exact tau targets `9.528149604797363 / 9.538650512695312 s`;
- exact legal tuning frequencies `0.05247608572244644 / 0.05241831764578819 Hz`;
- zero S firings over all 635 samples.

The high and low tau/period values lie inside the same P2 source cell. The required frequency inputs lie inside the shipping `[0.05,1.5] Hz` effective tuning bounds. The 13-sample tau images and frequency-to-target round trips are exact at binary32 precision.

This does **not** yet prove that the full `WavePeriodEstimator` can synthesize the frequency sequence from an admissible physical acceleration history. It does prove that the tau EMA is not the missing regularizer and that tau-cell subdivision alone cannot establish S recurrence.

## Failure classification / critic pass

**Failure type:** proof-method / incomplete source-state abstraction.

**Newly invalidated strategy:** refine only the P2 tau partition or retain the continuous tau EMA while still treating its tuning-frequency input as an arbitrary legal bounded source. The real timer reset resonance survives that refinement.

**Not invalidated:** the deployed filter, the physical operating envelope, the ordered-covariance feasibility signal, or a proof that derives stronger source regularity from the actual upstream estimator dynamics.

**Current limiter:** prove enough regularity of the complete `WavePeriodEstimator -> SeaStateAutoTuner -> tau EMA -> pseudo timer` chain to exclude indefinite S starvation, or demonstrate that the theorem must explicitly tolerate such zero-S intervals.

Strongest reason to abandon another tau-only construction: the destructive reset requires only a `1.92 us` pseudo-period change and a `0.141 ms` tau change. Those are already exactly realizable through the shipping tau EMA. More tau subdivisions do not remove the mechanism.

## Alternatives

1. **Derive upstream estimator regularity from shipping code.** Lift the canonical log-period state and, if necessary, its slow moment states until a rate/dwell/total-variation theorem excludes the reset cycle. Any constraint must be derived from implementation plus the existing physical source envelope; adding a convenient external regularity assumption is forbidden.
2. **Full source/timer reachability.** Carry the necessary `WavePeriodEstimator`/tau/timer states and propagate cumulative S information directly instead of proving a scalar packet-gap lemma. This is more expensive but does not discard the hybrid state that caused the failure.
3. **Finite-initialization / multiword covariance theorem.** Allow individual zero-S words and prove bounded covariance from the actual finite initialization plus later cumulative information. This helps only if the full implementation still guarantees recurrent S information; it cannot repair a genuinely indefinite starvation execution.
4. **Separate mathematical estimator stability from machine-lifetime clock liveness only if justified.** The 13-sample source-clock regime occurs at an extreme binary64 absolute-time exponent, while normal-age spacing is 21 samples. A deployment-lifetime theorem would have to be explicit and independently justified; silently replacing the all-lifetime shipping clock by nominal timing would be an operating-domain reduction and is not allowed.

## Next falsifiable experiment

**Lift one layer upstream into the exact binary32 canonical log-period smoother.**

Construct the coupled state

`(log_period_sec_, tau_applied, pseudo_update_elapsed)`

under the certified 13-sample source clock. Treat the valid raw moment-ratio period as the input, transcribe the shipping log-period horizon/update exactly, and search for a periodic no-S cycle whose resulting tuning frequency drives the already-certified tau/timer resonance.

Interpretation:

- **cycle exists:** the canonical period smoother also fails to exclude starvation; the next and only remaining implementation regularizer is the slow positive second-moment/filter state inside `WavePeriodEstimator`.
- **cycle does not exist:** certify the log-period state as the missing source regularity and build the P3 scheduler/information quotient around it.

This remains diagnostic only and cannot promote P3/P4/P5.

## Retained facts

- Canonical P3 usefulness threshold is exactly `1e-18`.
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
- tau-cell subdivision or tau-EMA state alone as a proof of S recurrence;
- blind subdivision, scalar-norm tightening, or deeper search as a substitute for source/timer memory;
- sigma/R_S coefficient, filter, domain, or canonical-gate tuning;
- additional P4 micro-certificates before P3 has a real margin.
