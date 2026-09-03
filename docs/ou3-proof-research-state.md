# OU-III proof research state

This is the current research ledger required by the root `AGENTS.md`. Replace stale state rather than accumulating PR history.

## State

**REPLAN AFTER PR #475 — the four-max “same-history” covariance upper is now rejected.**

Canonical P3 still has no theorem PASS. P4 remains blocked by P3 and P5 remains blocked by P4. The canonical usefulness gate remains exactly `delta >= 1e-18`; the deployed filter, declared operating domain and P2 V1 source language are unchanged.

PR #475 correctly identified the first structural defect: P3 compared a zero-start 13–26-sample translation floor (65–130 ms) with a 3.02–3.17 s covariance-word upper. At the point cell `(tau=12 s, sigma=6, R_S=400)`, crediting the deployed S cadence in the upper alone leaves about 6.66 orders to the gate, whereas making the floor whole-word too gives the non-certified point diagnostic `delta ~= 3.79e-4`, binding on `p`. This remains the feasibility signal, not a certificate.

PR #476 then tested whether the existing same-history upper could preserve enough source correlation while the floor was lengthened. It cannot. The upper summarizes a source history by four independent path maxima:

1. maximum pseudo-update cadence;
2. maximum `sigma^2`;
3. maximum `q_c = 2 sigma^2/tau`;
4. maximum S-measurement variance.

A legal P2 V1 path attains the **full global adverse label `[9,3,39,9]` in only 101 samples = 0.505 s**, far inside the 635-sample covariance word:

`729 --16--> 568 --16--> 407 --18--> 246 --25--> 85 --13--> 74 --13--> 63`.

The supporting gaps are exact P2 V1 labelled transitions. Once the four maxima are attained, max-update is absorbing and the finite-clock graph has no dead continuation, so this history extends across the complete word while retaining the global label. Therefore the global four-max label is not merely a Cartesian artifact: it is an **actual legal history label**.

### Failure classification

**The four independent path maxima are not a sufficient source-correlation quotient for the covariance upper.**

Any exact-elapsed, Pareto, minimum-cost or endpoint refinement that feeds the same four maxima into `translation_upper_from_summary` must eventually admit the actual global label above. Improving enumeration cannot fix the loss because the loss occurs when time order and residence duration are discarded.

This is separate from the original horizon mismatch. The required P3 repair now has two structural parts:

- a whole-word covariance lower valid for every admissible PSD initial covariance and legal source history;
- a **time-ordered, duration-aware covariance upper** that does not replace a brief visit to each adverse source by holding every adverse statistic for the whole word.

### Measured falsifications on PR #476

- A rigorous-but-overconservative global whole-word lower projection gave arbitrary-phase `delta >= 1.3872e-27` against the old global upper. This is still about nine orders below the gate and is rejected as the canonical architecture because it pairs small-tau/frequent-S lower behavior with the large-tau/sparse-S global ceiling. The old probe also has an ulp-level metric reconstruction defect (`_metric_lower` divided by a nominal scale while the canonical gate multiplies by an outward interval square), so it is retained only as falsification evidence and must not be promoted.
- A 49-segment gap-forgotten source frontier collapsed to one label `[9,3,39,9]`, admitting all 800 physical source nodes and reproducing the global upper. Rejected.
- A finite-cost adverse frontier was attempted next, but its implementation did not produce a target-crossing terminal class. It is superseded by the exact 16-mask witness above; fixing the enumerator cannot repair a four-max summary whose global label is already legally reachable in 0.505 s.

## Next falsifiable experiment

**Build a time-ordered, cadence-aware covariance-upper diagnostic on legal P2 V1 words. Do not run another covariance lower against the four-max envelope.**

The diagnostic must preserve, in order:

- source-local `tau`, `sigma`, `q_c` and `R_S` coupling;
- duration of each applied source segment;
- the actual S pseudo-measurement scheduler state across source commits. Shipping `set_pseudo_update_period_s()` applies `elapsed = fmod(elapsed, new_period)`, so fixed-cell firing counts such as 19 at `tau=12 s` and 602 at `tau=0.333 s` cannot simply be pasted onto a changing-source word;
- the complete 635-sample covariance-word horizon and canonical 0–25-sample terminal phase.

The first diagnostic is non-promoting. It should compare the old four-max upper with an ordered upper on the explicit 101-sample global-label witness extended to the covariance word. This is deliberately the adversarial history that defeats the current summary. If the ordered upper drops by the expected orders of magnitude, implement a source-complete finite-state/monotone quotient of that ordered upper. If it does not, reconsider the covariance-upper theorem representation rather than tune parameters or enclosures.

A safe upper construction may start from a finite covariance/observability bound and propagate source-ordered Riccati prediction plus guaranteed S updates using monotonicity. Do **not** use a naive independent-observation `1/N` information gain: process-noise correlations and the source-varying scheduler must remain covered.

## Retained facts

- Canonical P3 gate remains exactly `1e-18`.
- `OU3_P2_CORRELATED_STAGE_TRANSFER_V1` remains the source-language authority. Legal stage-boundary transitions are `(c,s) --g--> (s,t)`, `g in 13..26`; the following segment uses source `s`.
- The physical source partition remains 800 states. Arbitrary Cartesian tuner switching is forbidden.
- Same-tau `(sigma_index=0, R_S_index=0)` source nodes remain valid covariance-lower dominators within their tau cell by Riccati monotonicity in process `Q` and measurement `R`. This does not justify collapsing tau or the ordered history.
- Starting a covariance lower from `P=0` is valid below every admissible PSD initial covariance by Riccati monotonicity; the difficulty is retaining a useful lower over the whole legal word, not validity of the zero start.
- Endpoint-only segment SPD certification remains feasible. Intermediate per-step strict SPD was surplus and is not a theorem obligation.
- H=18 and A=21 full-state branches are both still required after translation P3 closes.
- Existing exact Joseph, Cayley, co-rotated accelerometer and reset identities remain structural facts only; they do not authorize P4 promotion.
- Current proof scope remains zero/disabled lever arm and dormant/transparent vibration-guard branch, with no replay fitting or operating-domain shrink for PASS.
- P4 is blocked until canonical P3 passes. P5 is blocked until canonical P4 proves strict whole-word contraction and must ultimately prove finite capture from the declared 45-degree entrance into the inner funnel.

## DEAD_ENDS / SHELVED

Do not resume these as a route to the P3 gate unless a new mathematical fact invalidates the rejection:

- fixed-source 635-sample lower as a theorem proof;
- global whole-word lower against the global four-max upper;
- 49-step gap-forgotten source superset;
- exact-elapsed/Pareto/minimum-cost enumeration **with the same four path maxima**;
- independent Cartesian `tau/sigma/R_S` extrema;
- recursive natural interval covariance boxes that forget the repeated scalar source parameter;
- blind deeper subdivision, Taylor-model/enclosure tightening, or isotropic `rho I` repair as a substitute for the structural horizon/upper fix;
- sigma/R_S coefficient or canonical-gate tuning;
- additional P4 micro-certificates before P3 has a real margin.

Sound but shelved items such as univariate Taylor models, tighter Loewner enclosures and the anisotropic floor may be revisited only after the ordered covariance upper exists and a real remaining P3 deficit is measured.

## P4 / P5 sequencing

Before any new P4 theorem producer, obtain the canonical P3 numerical verdict with the ordered whole-word translation comparison and both H/A joins. If P3 passes, run the non-promoting complete-word P4 contraction diagnostic first. Only a clear source-complete contraction result justifies rebuilding the rigorous P4 certificate. P5 then proves finite startup-to-inner-funnel capture; outer-sector entry alone is not P5 closure.
