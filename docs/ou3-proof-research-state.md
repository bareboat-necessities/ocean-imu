# OU-III proof research state

This file is the current research ledger required by the root `AGENTS.md`. Keep it short and replace stale state rather than accumulating PR history.

## State

**REPLAN. Both the prior plan and its first replacement are superseded.**

Research question: **make the segment floor and the covariance ceiling refer to
the same horizon.** `delta` is certified through
`D_h L_segment D_h' - delta * diag(Sigma_upper) >= 0`, where `L_segment` is a
zero-start floor over one 13--26 sample segment (65--130 ms) and `Sigma_upper` is
a whole-word ceiling over 3.02--3.17 s. Those are quantities of different
horizons, and on the binding `S` channel the mismatch is worth about 21 orders.
Everything else queued -- univariate Taylor models, tighter Loewner enclosures,
deeper subdivision, the isotropic-floor repair, the sigma/R_S cost reduction --
stays **shelved**.

### Measured, at the cell that binds the ceiling (tau = 12 s, sigma = 6, R_S = 400)

| step | S-channel value | `delta` | vs gate |
| --- | --- | --- | --- |
| today: 130 ms zero-start floor vs 3-point whole-word ceiling | floor 1.37e-20, ceiling 7.17e9 | 1.90e-30 | 11.7 orders short |
| credit the real `S` cadence in the ceiling only | ceiling 6.18e4 | 2.21e-25 | **6.66 orders short** |
| make the floor a whole-word quantity too | floor 41.2 | 3.79e-4 | clears by 14.6 orders |

`delta` here is evaluated at the measured phase-0 identity floor
`rho = 8.742479e-07` (node 80 / gap 13), not the canonical minimum over phases,
x-subcells and all 800 nodes, which is 1.89e-31. Use the ratios, not the
absolute values.

The `S` channel binds `delta` in both of the first two rows: its per-channel
bracket is 1.90e-30, against 1.96e-20 for `v` and 1.60e-09 for `a_w`. `a_w` is
21 orders slack, so nothing that only improves `a_w` can move `delta` at all.

In the third row the binding channel **moves to `p`**, because a whole-word floor
lifts the three integrator channels to within a factor of a few thousand of the
ceiling and `p` is the tightest of them: floor/ceiling is 6.48e-4 on `v`,
**3.79e-4 on `p`**, 6.66e-4 on `S` and 3.88e-1 on `a_w`. `delta` is the minimum
over channels, so it is 3.79e-4. Expect the binding channel to change again once
these become rigorous uniform bounds rather than point recursions.

At this cell one 130 ms segment carries **zero** `S` firings. The floor is
therefore a pure free-integration quantity while the ceiling it is compared
against is a whole-word one.

### Why the previous replan was wrong

The prior version of this section claimed the ceiling was "12 to 13 orders
looser than the covariance the filter settles at" and predicted that crediting
the `S` cadence would raise `delta` by 12 orders and clear the gate. **Both
claims are withdrawn.** They came from comparing a Riccati fixed point computed
at one cell (tau = 3, sigma = 1, `rs.lo`) against `33459 m/s`, which is the
ceiling's maximum over the *whole* tuner range (attained at tau = 12, sigma = 6,
`rs.hi`). Mixing cells inflated the ratio by about eight orders.

Compared at matched cells the ceiling looseness is real but far smaller:

| cell (tau, sigma, R_S) | N firings in the word | 3-point v ceiling | N-firing v ceiling | variance gain |
| --- | --- | --- | --- | --- |
| 12, 6, 400 (**binds**) | 19 | 33434 m/s | 256.7 m/s | 1.70e4 |
| 12, 6, 0.15 | 19 | 875 m/s | 1.377 m/s | 4.04e5 |
| 1.1, 1, 400 | 201 | 5358 m/s | 83.3 m/s | 4.14e3 |
| 0.333, 6, 0.15 | 602 | 189.9 m/s | 1.609 m/s | 1.39e4 |

Over all matched cells the gain ranges 1.08e3 to 3.31e8 in variance, and at the
binding cell it is 1.70e4 on `v`, 2.59e4 on `p`, 1.16e5 on `S`. That is worth
5.06 orders on `delta` -- real, but it leaves 6.66.

The firing count was also wrong. The deployed cadence is tau-scaled,
`T_S = clamp(PSEUDO_UPDATE_TAU_RATIO_DEFAULT * tau, 0.005, 0.25)` with the ratio
`0.015 / 1.1`, and the word length is `2*max(1, 2*(T_S+h)) + (T_S+h) + 1`. The
two move together, so a word carries **19 firings at tau = 12 s and 602 at
tau = 0.333 s**; "about 211 per 3.17 s word" combined the longest word with a
different cell's cadence. `periodic_update_due` accumulates `dt` and does not
fire at `t = 0`, so these counts are `floor(T_word / T_S)`.

### Failure classification

**Proof-method failure: a horizon mismatch between the two sides of the master
inequality.** The previous entry called it a ceiling-only failure; that was
measured wrong and is corrected above. The ceiling is genuinely loose (5 orders)
but is not the limiter.

### What this invalidates

* `translation_upper`'s three-point inversion as the ceiling -- still true, now
  quantified at 5.06 orders rather than 12;
* the claim that fixing the ceiling alone reaches the gate;
* the earlier conclusion that "domain realism is not a route to the gate", which
  was measured against the mis-stated ceiling and is not evidence either way.

### What it does not invalidate

* the `1e-18` canonical gate;
* ~~same-history P2-V1 source correlation and the envelope machinery, which is
  not what broke~~ -- **withdrawn.** PR #476 showed the four-max summary admits
  the full global adverse label `[9,3,39,9]` after only 101 samples, far inside
  the 635-sample word, so at whole-word horizon the same-history *upper*
  degenerates to the global one. Verified independently here: all six
  transitions of `729 --16--> 568 --16--> 407 --18--> 246 --25--> 85 --13--> 74
  --13--> 63` are legal P2 V1 edges. The source *language* stands; the four-max
  quotient built on it does not;
* the H/A attitude and bias floors, which are separate from translation;
* the exact Joseph, Cayley, co-rotated accelerometer and reset identities;
* the endpoint-only SPD result, which remains correct and is not the limiter.

### Next falsifiable experiment

The horizon fix has two structural parts. PR #476 owns the second; this ledger
entry owns the first.

1. **whole-word covariance lower**, valid for every admissible PSD initial
   covariance and every legal source history;
2. **time-ordered, duration-aware covariance upper** -- #476's lane, because the
   four-max quotient it replaces is what broke.

The lower must be certified over *mixed* legal words, not one source held for a
word. The single-cell dwell measured below is a feasibility signal for it, not a
substitute: it says the geometry has room, and says nothing about whether a
finite-state or monotone argument can quantify over legal histories.

**Predicted:** a whole-word lower certified over mixed legal words retains a
`delta` within about two orders of the worst single-cell dwell, so above `1e-10`.
**Result: survives** -- greedy adversarial descent reaches 2.97e-8, a factor of
1.184 below the dwell. The remaining obligation is a *global* argument over legal
words rather than a local search: a finite-state or monotone quotient that bounds
every legal word, not the best one a greedy descent can find.

## Closed-form margin re-run under PR #478's clamp values

Applied #478's clamps to a scratch copy of the source, re-ran the full 800-node
tool, then restored `src/` (verified clean).

| | main clamps | #478 clamps | gain |
| --- | --- | --- | --- |
| per-node worst | 1.505e-05 | **2.802e-04** | 1.27 orders |
| mixed-word, all legal words | 3.809e-09 | **7.979e-08** | 1.32 orders |

800/800 nodes validated in both runs. Under #478's clamps the per-node margin
clears the `1e-18` gate by **14.45 orders** and the mixed-word bound by **10.90
orders**. The binding node moves from 780 to 790 and the max S gap falls from 33
to 31 samples.

### `MAX_SIGMA_A` does not reach the proof tooling at all

`ou3_source_reachable_matrix_p3.py:291` hardcodes
`"sigma_aw_applied_safety": [0.05, 6.0]` instead of parsing `MAX_SIGMA_A` from
the shipping header. So:

* #478's `MAX_SIGMA_A` 6 -> 4 change is **invisible** to the proof chain, which
  is part of why it measures at 0.00 orders;
* the "#478 clamps" run above still used `sigma` up to 6.0, so it is conservative
  -- the true figure under those clamps is slightly better than reported here;
* it is a source-of-truth gap of exactly the kind that produced the
  `PSEUDO_RATIO` defect recorded earlier in this ledger, where a hand-written
  constant silently disagreed with the shipping one.

Worth fixing on its own merits, independently of whether #478 lands.

## No clamp change can make canonical P3 pass

Measured directly against `translation_upper`, since the canonical floor does not
depend on the clamps. Canonical today is 11.72 orders short.

| clamp change | worst upper `S` | `delta` | still short by | gain |
| --- | --- | --- | --- | --- |
| none (today) | 7.171e+09 | 1.905e-30 | 11.72 | -- |
| `MAX_TUNE_FREQ_HZ` 1.5 -> 1.2 | 7.171e+09 | 1.905e-30 | 11.72 | **0.00** |
| `MAX_SIGMA_A` 6 -> 4 | 7.168e+09 | 1.906e-30 | 11.72 | **0.00** |
| `MAX_R_S` 400 -> 100 | 4.528e+08 | 3.017e-29 | 10.52 | 1.20 |
| `PSEUDO_UPDATE_PERIOD_MAX` 0.25 -> 0.15 | 3.442e+09 | 3.968e-30 | 11.40 | 0.32 |
| **all four (#478)** | 2.160e+08 | 6.324e-29 | **10.20** | 1.52 |

And driving the dominant lever to absurdity does not help either -- `MAX_R_S`
**saturates**:

| `MAX_R_S` | `delta` | short by |
| --- | --- | --- |
| 100 | 3.017e-29 | 10.52 |
| 10 | 1.456e-27 | 8.84 |
| 1 | 2.761e-27 | 8.56 |
| 0.15 (the minimum of its own range) | 2.786e-27 | **8.56** |

**Taking `R_S` to the bottom of its declared range still leaves 8.56 orders.**

### Why it saturates, and what that means

`translation_upper` builds `rstack = 3 (rmax + s_nuis + s_proc)`. As `R_S` falls,
`rmax` vanishes and the **process** corruption terms `s_nuis` and `s_proc`
dominate. They do not depend on `R_S` at all, so the ceiling floors out.

`translation_upper` is limited by process corruption, not by measurement noise.
That is the same fact that made endpoint referencing worth 2.9 orders while
crediting more observations was worth only 1.2: shortening the span over which
process noise corrupts each observation is what pays, not strengthening the
measurement.

**Consequence for #478.** Its clamp tightening is worth 1.52 orders on the
canonical margin and does not make P3 pass, so it has to be justified as
theorem-envelope hygiene on its own merits rather than as a route to P3. Two of
its four clamp moves -- `MAX_TUNE_FREQ_HZ` and `MAX_SIGMA_A` -- buy **nothing**
here, 0.00 orders each.

**Consequence generally: stop looking at clamps.** The deficit is structural.

## The closed-form route CANNOT be fed to the canonical gate as it stands

The chain exists on paper: `ou3_p3_p2_v1_full_state_join.build()` accepts a
`translation_candidate`, and `ou3_p3_canonical_gate.build()` accepts a
`p3_candidate`, so a parallel path could in principle let the **gate** render the
verdict rather than this ledger asserting one. It cannot, and the reason is worth
recording precisely.

`TRANS.validate` requires 16 flags true. Audited against the closed-form
construction, **6 can be claimed honestly, 1 is partial, and 9 cannot**:

| can claim | cannot claim |
| --- | --- |
| `source_only` | `full_word_history_sufficient_quotient_consumed` |
| `zero_lever_arm_branch` | `process_covariance_measurement_bounds_same_source_history` |
| `dormant_transparent_vibration_guard_branch` | `same_history_covariance_evaluated_before_uniform_endpoint_envelope` |
| `Riccati_order_monotonicity_used` | `full_4x4_translation_matrix_retained_inside_source_segments` |
| `strongest_accelerometer_and_S_measurements_applied_every_sample` | `one_complete_segment_common_boundary_floor_used` |
| `magnetometer_translation_jacobian_zero_on_declared_branch` | `older_process_covariance_discarded_at_each_boundary` |
| | `all_finite_clock_sample_phases_covered` |
| | `phase_26_is_next_stage_boundary_when_clock_advances` |
| | `frozen_clock_absorbing_hold_branch_included` |

plus the structural demands that `phase_samples == list(range(26))`,
`endpoint_phase_count == 800*26`, and `boundary_floor` carry
`all_800_sources_checked` with a strict `rho_z_identity_lower`.

**The nine are not deficiencies in the closed-form bound. They describe the
stepwise segment/phase architecture that it deliberately replaces.** There are no
segment boundaries to discard covariance at, no 0..25 phase sweep because
endpoint referencing removes that degree of freedom, and no same-history quotient
because the mixed-word bound is deliberately cross-worst -- which is *stronger*,
covering every legal word rather than one history at a time.

Two are real scope limits rather than architecture mismatches, and should not be
glossed with the rest: the construction is 3x3 on `(v, p, S)` and does not carry
`a_w`, and the frozen-clock absorbing hold is not modelled.

### Consequence

Feeding the gate would require either setting flags that do not describe the
construction, which is fabrication and is not an option, or revising the
translation-candidate contract so it admits a closed-form architecture. The
second is a proof-architecture decision for the project, not something this
branch should do unilaterally: those flags encode real conservatism guarantees,
and `older_process_covariance_discarded_at_each_boundary` in particular is a
safety property whose closed-form analogue has to be argued, not dropped.

**So this branch cannot produce a P3 verdict, and should not pretend the
remaining work is plumbing.** What it can offer is a measured alternative
construction with its own conservatism argument, for the project to accept or
reject on the merits.

## What replacing the canonical pair requires: both sides, not one

Swapping the closed-form ceiling into the canonical producer alone does **not**
reach the gate. The canonical floor is the limiter once the ceiling is repaired:

| floor | ceiling | `delta` | vs `1e-18` |
| --- | --- | --- | --- |
| canonical | canonical (today) | 1.905e-30 | short by 11.72 |
| canonical | closed form | 1.859e-25 | **short by 6.73** |
| closed form | canonical | 4.419e-14 | clears by 4.65 |
| closed form | closed form | 4.311e-09 | clears by 9.63 |

So this is not a matter of plugging a better `Sigma_upper` into
`phase_row`. The canonical floor is a stepwise zero-start segment quantity
(1.366e-20 on `S`) and the closed-form floor is an information-form endpoint
quantity (3.169e-04); they are different constructions and only the matched pair
clears. Replacing both is a change to the producer's core rather than wiring
into it, and it needs the canonical gate to render the verdict, not this module.

### Floor premise verified on both measurement channels

The floor rests on the other channels having no direct translation dependence.
Checked in the source, not argued:

* `measurement_update_acc_only` builds its innovation covariance from `OFF_TH`,
  `OFF_AW`, `OFF_BA` and `OFF_BG` only;
* `measurement_update_mag_only` (lines 2122-2198) references only
  `P_th_th = Pext.block<3,3>(0,0)`. The `OFF_AW`/`OFF_BA` uses just below it
  belong to `accelerometer_measurement_func`, which starts at 2199, not to the
  magnetometer.

Neither has `v`, `p` or `S` columns, so both reach translation solely through the
cross-covariance with `a_w`, which is the channel the floor already accounts for.

## Mixed-source words are covered without enumeration

The per-node rows pair one node's ceiling with its own floor, which is what
canonical `delta` does but says nothing about a word that changes source
mid-flight. Two monotonicity arguments close that without enumerating words:

* **Ceiling.** The certified S recurrence puts an observation in every max-gap
  window, so the worst a legal word can do is place them as late as possible in
  each -- exactly the max-cadence dwell. Any other legal word has at least as
  many observations, at least as close to the endpoint, hence a larger Gramian
  and a smaller ceiling.
* **Floor.** No word can carry more than one S observation per sample, and the
  min-cadence dwell already attains that, since the cadence clamps to `h` there.
  It is the most informative case and yields the smallest floor.

| | ceiling | floor |
| --- | --- | --- |
| `tau` | 12.0 | 0.3333 |
| `R_S` | 400 | 0.15 |
| S observations | 19 | 635 (one per sample, the hard limit) |

| channel | floor (var) | ceiling (var) | ratio |
| --- | --- | --- | --- |
| v | 2.510e-04 | 6.589e+04 | **3.809e-09** |
| p | 6.727e-04 | 1.664e+05 | 4.043e-09 |
| S | 3.169e-04 | 7.350e+04 | 4.311e-09 |

**Mixed-word margin 3.81e-09, valid for every legal word, clearing the `1e-18`
gate by 9.58 orders.** Binds on `v`. Strictly more pessimistic than the per-node
1.505e-05 because it pairs cross-worst, and strictly more general.

Four further tests pin it: the cross-worst pairing must not come out better than
the same-history one, the ceiling must take the fewest observations while the
floor takes the most, the floor must sit exactly at one observation per sample
or it is not conservative, and `validate()` must reject a dropped
source-changing claim.

## Per-node margin: 800/800 nodes, interval-certified

`tools/ou3_p3_closed_form_word_margin.py`. Both sides of the P3 comparison built
in closed form from one endpoint-referenced observability Gramian, so neither
steps a Riccati recursion -- which is what makes it survive interval arithmetic
where three previous attempts died between 343 and 503 of 635 samples.

| | value |
| --- | --- |
| source nodes evaluated / validated | **800 / 800** |
| worst margin | **1.505e-05** |
| clears the `1e-18` gate by | **13.18 orders** |
| binding node / channel | 780 (`tau` in [8.386, 12.0]) on `v` |
| S observations, ceiling / floor | 19 / 27 |

Run under **current main's** clamps, `MAX_R_S = 400`, not #478's tightened ones,
so the result does not depend on that PR landing. #478's clamps would improve it
by a further 1.29 orders.

Within each source cell the two sides take opposite endpoints -- the ceiling the
fewest observations and weakest S measurement, the floor the most and strongest
-- so both hold for every parameter value the cell admits. `R_S` comes from the
deployed clamped SpectralMSE law rather than the free P2 cell range, and `P0`
from the shipping constructor.

Nine tests pin the contracts: cadence and `R_S` agree with the repo's own source
parser and respect both clamps, the Gramian grows with observation count, the
uninflated floor Gramian dominates the inflated ceiling one, floor never exceeds
ceiling on any channel, a singular Gramian is reported rather than silently
inverted, and `validate()` catches tampering with the gate or any
non-promotion flag.

**This is a margin, not a P3 verdict.** The canonical producer and gate remain
the promotion authority; the artifact carries `certifies_theorem_stage: false`
and cannot set a theorem flag. What remains between this and a certificate is
the canonical producer wiring and mixed-source words, neither of which is
blocked by the interval wall now that both sides are closed form.

## How the floor is justified against the other measurement channels

The defect in the first attempt -- a floor that ignored accelerometer information
and so could be violated by the truth -- is closed by a structural argument
rather than by a tighter number.

**The accelerometer cannot drive v, p, S below their initial uncertainty.** Its
measurement Jacobian has zero columns for `v`, `p`, `S`: the innovation
covariance in `measurement_update_acc_only` references only `OFF_TH`, `OFF_AW`,
`OFF_BA` and `OFF_BG`. It touches the translation states purely through the
cross-covariance with `a_w`. And in the limiting case where `a_w` is known
*exactly* for the whole word, `v(T) = v(0) + int a_w` still carries `v(0)`'s
uncertainty, and likewise for `p` and `S`. The same holds for the magnetometer,
which informs attitude only.

So a floor safe against every other channel credits the `S` observations at their
**full** strength and keeps the raw shipping prior:

* `R_eff = rmax`, with **no** process inflation -- perfect `a_w` knowledge is the
  most informative case, hence the smallest and therefore most conservative
  floor;
* the raw shipping `P0` rather than a propagated one, since propagation only adds
  process noise and so only reduces information: `Y0_propagated <= Y0`, giving
  `(Y0 + G)^-1 <= (Y0' + G)^-1 <= truth`.

The ceiling keeps the pessimistic process-inflated `R_eff`. The two now differ by
a real mechanism instead of by a shared term, and the tell-tale `delta ~ 0.999`
is gone:

| tau | sigma | N | ceiling S | floor S | delta |
| --- | --- | --- | --- | --- | --- |
| 0.417 | 0.05 | 558 | 0.02041 | 0.01909 | 0.4744 |
| 0.417 | 4.0 | 558 | 0.03286 | 0.01909 | 0.003259 |
| 3.0 | 0.05 | 77 | 0.2369 | 0.2331 | 0.9169 |
| 6.0 | 4.0 | 38 | 50.37 | 23.54 | 0.0004839 |
| 12.0 | 4.0 | 21 | 61.70 | 26.97 | **0.0002803** |

**Worst `delta` = 2.80e-4, clearing the `1e-18` gate by 14.45 orders**, binding on
`v` at `tau = 12`, `sigma = 4`, `R_S = 100`. Both inversions validate in interval
arithmetic at all 15 cells.

### What is still missing before this is a certificate

* 15 cells, not the 800-node source partition;
* single-source dwell words, not mixed legal words;
* not run through the canonical P3 producer or its gate, which remain the
  promotion authority;
* the magnetometer argument is stated from the same structure but only the
  accelerometer Jacobian was checked in the source.

The construction is now closed-form on both sides, so none of these are blocked
by the interval wall. They are enumeration and plumbing rather than open
questions.

## Superseded first attempt: a floor that shared the ceiling's Gramian

Closed-form floor from the same Gramian as the ceiling: information adds, so
`Y_end <= Y0_propagated + G` and hence `P_end >= (Y0' + G)^-1`, with `P0` a
shipping constant so `Y0` is bounded. No recursion, so the interval wall does not
apply, and both inversions validate at all 15 cells.

| tau | sigma | N | ceiling S | floor S | delta | vs `1e-18` |
| --- | --- | --- | --- | --- | --- | --- |
| 0.417 | 0.05 | 558 | 0.02041 | 0.02041 | **0.99947** | +18.0 |
| 3.0 | 1.0 | 77 | 3.102 | 2.781 | 0.6639 | +17.8 |
| 6.0 | 1.0 | 38 | 34.27 | 18.63 | 0.01436 | +16.2 |
| 12.0 | 0.05 | 21 | 25.61 | 15.55 | **0.001683** | +15.2 |

Worst `delta` 1.68e-3, nominally clearing the gate by 15.23 orders.

**That number should not be trusted, and `delta ~ 0.999` at low `tau` is why.** A
floor and a ceiling that agree to 0.05 percent are not two independent bounds --
they are one construction differenced against itself. The floor here is the
ceiling's own Gramian minus a prior term, so the agreement is an artifact of
sharing `G`, not evidence that the covariance is pinned.

**The concrete defect:** `G` counts only `S` observations. The deployed filter
also runs accelerometer updates, which inform `v`, `p` and `S` through `a_w` and
the cross-covariances. Extra information means smaller covariance, so the true
`P` can fall **below** this floor, and a lower bound that the truth can violate
is not a lower bound. The canonical floor starts from `P = 0` precisely to be
safe against every information source, which is why it is 1.366e-20 and not
this.

So Step 2 is **not** done. What the computation does establish is narrower and
still useful: the closed-form route survives interval arithmetic on both sides,
and if a floor can be justified that credits the S observations at all, the
geometry has ample room. Fixing it means bounding the information the other
measurement channels can contribute, which is the real Step 2 and is not
attempted here.

## Stacked deficit accounting after PR #478

Re-ran the ceiling under #478's tightened safety clamps. All figures are the
binding `S` channel against the canonical segment floor 1.366e-20.

| step | `S` ceiling (variance) | `delta` | orders short of `1e-18` |
| --- | --- | --- | --- |
| canonical today: 3-point, start-referenced, old clamps | 7.171e+09 | 1.91e-30 | **11.72** |
| + endpoint referencing (#478's fix, and independently step 1 here) | 5.324e+07 | 2.57e-28 | 9.59 |
| + N-point Gramian instead of 3-point (step 1 here) | 7.350e+04 | 1.86e-25 | 6.73 |
| + #478's tightened safety clamps | 3.807e+03 | 3.59e-24 | **5.45** |

**6.27 orders recovered; 5.45 remain, and they must come from the floor.** The
earlier figure of 8.6 predates both #478's clamps and the ceiling work.

### What #478's clamps change

| constant | old | new | effect here |
| --- | --- | --- | --- |
| `MAX_TUNE_FREQ_HZ` | 1.5 | 1.2 | `tau` floor rises 0.3333 -> 0.4167 s |
| `MAX_SIGMA_A` | 6.0 | 4.0 | less process noise at the corner |
| `MAX_R_S` | 400 | 100 | 16x stronger worst-case `S` measurement |
| `PSEUDO_UPDATE_PERIOD_MAX_S` | 0.25 | 0.15 | max S gap 33 -> 31 samples |
| `kDynamicEmaHorizonMaxSec` | 30 | 35 | loosened, leaves the largest reference sea interior |

Worth **1.29 orders** on the worst endpoint `S` ceiling on its own, taking it
from 271.1 to 61.70 m.s. These are safety limits that no reference record
approaches, so this is a theorem-envelope tightening rather than a deployed
behaviour change -- but that claim is #478's to carry, and re-running the
deterministic sim against it is still outstanding here.

This ledger recorded the same corner as "unphysical but nearly free", measured at
4x on `a_w` and 1.002x on v/p/S. That measurement was taken against the
**broken** ceiling, where 1.29 orders was worth nothing against a 12.7-order
deficit. Against the repaired ceiling it is a fifth of what remains. **The
conclusion "not worth pursuing" was an artifact of the ceiling error**, and is
withdrawn.

### Endpoint referencing was found independently on both sides

#478's covariance upper now states that the reconstructed `[v,p,S]` covariance
"is already an endpoint covariance and MUST NOT be propagated forward through the
word again". That is the same fix as step 1 here, reached separately. #478 keeps
three selected observations; the N-point Gramian is worth a further 2.86 orders
on top of it.

## PR #478 supersedes the starvation work here and hands Step 2 its input

#478 repairs `set_pseudo_update_period_s` so a cadence retarget preserves
already-earned service credit instead of reducing elapsed time modulo the new
period:

```
if (elapsed < period) return elapsed;          // preserve
return std::nextafter(period, T(0));           // park just below the deadline
```

Verified here against this branch's own starvation witness, which #478's change
kills outright:

| retarget semantics | S firings over 635 samples |
| --- | --- |
| old `fmod` | **0** |
| new progress-preserving | **24**, worst gap 27 |

and the source-independent 33-sample bound reproduces exactly: the largest
deployed period 0.16363635659217834 s over `h` is 32.727, so 33 samples, attained
at `tau = 12`.

**Consequences for this branch.**

* The `tau` slew argument recorded below, and Step 0 of the plan built on it, are
  **obsolete**. #478 makes S recurrence a property of the implementation rather
  than something to be argued from the estimator's moment horizon. That is a
  better outcome than the argument it replaces.
* This branch concluded "no filter changes are warranted" on the evidence that
  the firing rate is nominal under realistic jitter and starvation is
  unreachable in practice. That was defensible but **it was the weaker call**:
  #478 removes the latent fragility instead of arguing it is unreachable, and
  gets a positive recurrence certificate for it.
* `tools/ou3_p3_whole_word_lower_feasibility.py` mirrors the **old** `fmod`
  retarget, and its four `StarvationWitness` tests assert that starvation
  happens. Both must be re-pointed when #478 lands, or they will pin behaviour
  the shipping code no longer has.

**What Step 2 gains.** The closed-form floor needs a rigorous bound on the time
since the last S firing. #478 certifies exactly that, source-independently, at
**33 samples = 0.165 s**. Step 2 no longer has to derive it.

## Step 1 SUCCEEDS: an interval-certified closed-form ceiling exists

The first rigorous object in this line of work. Closed form, so the stepwise
interval wall recorded below does not apply, and it **inverts cleanly in interval
arithmetic at every cell tested** -- the 3x3 determinant never straddles zero.

Construction: the S observations constrain the state at the **word endpoint**
directly, each through the backward map `[(t_k - T_w)^2/2, (t_k - T_w), 1]` with
its own process-corruption inflation, giving `G = sum_k h_k h_k' / R_eff,k` and
`P_end <= G^-1`. No recursion, no `P0`.

| construction | worst-case gain on the binding `S` channel |
| --- | --- |
| start-referenced, matching `translation_upper` | 1.24 orders |
| **endpoint-referenced** | **4.10 orders** (up to 8.28) |

**Endpoint-referencing is worth 2.9 orders on its own.** That is the more useful
finding: `translation_upper`'s looseness is not mainly its three-point count. It
builds the bound at word *start* and then free-runs the estimate through `F` for
a whole word, and that free-running amplification dominates everything else. Two
earlier comparisons in this ledger overstated the benefit of crediting the
cadence because they compared unequal quantities -- bare `R_S` against
`translation_upper`'s process-inflated `rstack`, and an at-start covariance
against an at-endpoint one. Both are corrected here.

**Where this leaves the deficit:** 12.7 orders needed, 4.10 delivered by the
ceiling. **8.6 orders must come from the floor**, which is Step 2 and is not yet
attempted. Step 1 retires the constructive risk on the ceiling only; it is not a
P3 verdict, has not been run through the canonical producer, and quantifies over
15 cells rather than 800.

## NO CERTIFICATE EXISTS YET. Everything else in this ledger is a diagnostic.

Stated plainly because the distinction was blurred: canonical P3 has never
produced a passing artifact, and nothing in this branch is interval-certified.
Every number recorded here is a double-precision or mpmath **point** evaluation,
at 15 of 800 cells, at one word length. The artifacts say so themselves --
`certifies_theorem_stage: false`, `interval_certified: false`,
`quantifies_over_legal_histories: false`. They are useful for locating the
limiter and for killing wrong hypotheses. They are not proof of anything.

### Attempted and FALSIFIED: validated invariant upper by interval Riccati

The natural rigorous construction was tried and does not work. The idea was
sound: the Riccati word map is Loewner monotone, so certifying
`D - W_c(D) >= 0` at the single box corner `D` certifies the whole box, meaning
no covariance *box* need be propagated and only rounding accumulates.

It still dies. Propagating the corner through one 635-sample word in interval
arithmetic loses positivity of the measurement denominator partway through, even
in **Joseph form** -- the covariance form dies within a few steps, since it
subtracts two nearly equal large quantities whenever `R` is small against
`P[S][S]`.

| `x = h/tau` cell width | samples survived of 635 |
| --- | --- |
| 1.796e-04 (full cell) | 343 |
| 2.806e-06 (64 subcells) | 434 |
| 4.384e-08 (4096 subcells) | 503 |

Narrowing `x` by 4096x buys 160 of the 292 missing steps, about 35 steps per 8x
narrowing. Reaching 635 needs roughly `1e7` subcells against a budget already
strained at 64. **This is the third independent time this project has hit the
same wall**, and it is the same one already frozen in DEAD_ENDS as "recursive
natural interval covariance boxes that forget the repeated scalar source
parameter at every Riccati operation". Evaluating at a single corner does not
escape it, because the transition itself is interval-valued through the `tau`
cell.

**Consequence: a rigorous whole-word upper cannot be built by stepwise interval
propagation.** It has to be closed-form. The tractable candidate is an
`N`-point observability-Gramian inversion -- structurally what
`translation_upper` already does with three points, evaluated instead over all
`N` firings the cadence actually delivers. That construction has no stepwise
recursion, so the dependency wall above does not apply to it. Measured worth of
crediting the real cadence, from the point diagnostics: 5.06 orders.

The probe was removed rather than retained, matching how #476 handled its own
falsified probes.

## Point diagnostics: the certificate needs no new theorem conditions

Earlier entries in this ledger asked whether a `P0`-bounded P3 would need a new
declared entrance bound, and treated that as a specification decision. **It does
not, and the framing was wrong.** `P0` is a shipping constant:

```
const T sigma_v0 = T(1.0);    // m/s
const T sigma_p0 = T(20.0);   // m
const T sigma_S0 = T(50.0);   // m.s
set_initial_linear_uncertainty(sigma_v0, sigma_p0, sigma_S0);
```

with `a_w` seeded at `Sigma_aw_stat`. The `S` entrance bound this ledger called
"undeclared" is declared: 50 m.s.

### P2 is a Cartesian product over parameters the source computes jointly

Two blockers, one cause. The P2 quotient ranges `(tau, sigma, R_S)` independently
inside their cells, but the deployed code **derives** `R_S` from `tau` and
`sigma` through the SpectralMSE law
`R_S = coeff * q_eff^(1/14) * (sigma*tau^4)^(6/7) / sqrt(T_S)`:

| tau | sigma | `R_S` the law gives | P2 cell max | reachable? |
| --- | --- | --- | --- | --- |
| 12 | 0.05 | **39.77** | 400 | no -- 10.1x below |
| 12 | 1.0 | 518 -> clamped 400 | 400 | yes |
| 3 | 0.05 | 0.686 | 400 | no |
| 0.333 | 0.05 | 0.0010 -> clamped 0.15 | 400 | no |

The adverse corner the covariance ceiling uses, `(tau=12, sigma=0.05, R_S=400)`,
is **not reachable**. Likewise `tau` is slew-limited by the estimator moment
horizon and the tuner EMA, which is what makes the starvation resonance
unreachable. Both are properties of the shipping code, so importing them is a
parsed source fact, not a domain shrink.

### Measured: the word map has an invariant set, from source constants alone

Iterating the joint word map `P <- max over law-consistent cells of
word_map(cell, P)` from the declared `P0` converges in 81 iterations:

| | v | p | S | a_w |
| --- | --- | --- | --- | --- |
| `P*` (variance) | 1578.4 | 10879 | 25699 | 35.80 |
| as std | 39.7 m/s | 104 m | 160 m.s | 5.98 m/s^2 |

`P*` is **invariant under every law-consistent cell** (verified, not argued), and

    worst delta over all law-consistent cells = 1.26e-8
    -> clears the canonical 1e-18 gate by 10.10 orders

So `Sigma_upper = P*`. There is no new hypothesis, no new declared bound, and no
weakening of what P3 asserts about initial covariance -- the theorem simply stops
demanding `P0`-freedom, which its own word map never needed, and uses the
invariant set instead. The `P0`-free Gramian inversion was a design choice in the
proof, not a requirement of it.

`P*` is large at the extreme corner because that corner is `sigma_aw = 6 m/s^2 =
0.61 g` with `tau = 12 s`, already recorded above as unphysical; at
`tau = 3, sigma = 1` the fixed point is v 2.3 m/s, p 3.2 m, S 2.6 m.s.

**Still required:** starved words must be excluded, because with zero `S` firings
the map has no fixed point at all -- `S` grows without bound. That is the `tau`
slew argument below. Everything else is now measured.

Not a certificate: double-precision point diagnostic over 15 cells at one word
length, no interval arithmetic.

## Supporting detail: P2's tau model is over-approximate

Starvation is not a physical regime. It is an artifact of the P2 quotient letting
`tau` jump anywhere inside a cell at every stage boundary, which the deployed
chain cannot do.

**Starvation is bit-exact.** Perturbing the resonant `tau` sequence by a single
binary32 ulp at 24 of the 49 stages restores `S` firing, and the magnitude of the
perturbation is irrelevant -- 1 ulp and 256 ulps break exactly the same 24
positions. Only *which* stage is perturbed matters.

**The deployed chain rate-limits `tau` by about 2.5 orders more than the
resonance can tolerate.** Starvation needs the tuning frequency to alternate by
1.1015e-3 relative every 13 samples = 65 ms, between #476's certified pair
`0.05241831764578819` and `0.05247608572244644` Hz. But the deployed
`WavePeriodEstimator` carries an exponentially weighted moment state with a
horizon of **at least 20 s** (`min_horizon_sec = 20.0`,
`moment_horizon_periods = 4.0`, up to `max_horizon_sec = 180.0`), plus a
canonical log-period smoothing state on top:

| estimator horizon | reachable fraction in one 65 ms stage | required moment-ratio swing |
| --- | --- | --- |
| 20 s (fastest allowed) | 3.24e-3 | 33.9 percent |
| 40 s | 1.62e-3 | 67.8 percent |
| 90 s | 7.22e-4 | 152.6 percent |
| 180 s | 3.61e-4 | 305.1 percent |

Even at the fastest allowed horizon the output moves only 3.2e-3 of the way to a
new input per stage, so driving the required alternation demands a **34 percent
relative swing in the moment ratio itself, every 65 ms, in exact antiphase, for
49 consecutive stages, landing on binary32-exact values at each one**. The moment
ratio is a ratio of spectral moments of wave acceleration and moves on wave
timescales. Downstream of it the tuner EMA adds its own limit,
`alpha = 1 - exp(-dt/adapt_sec) <= 0.095` per sample from the
`adapt_sec >= 0.05 s` clamp.

This is not a proof that the estimator cannot produce it, and it is not claimed
as one. It is a measured statement that the deployed chain low-passes the
required alternation by about 2.5 orders.

### The route that needs no weakening

1. Carry the deployed **`tau` slew limit** into the P2 source model as a parsed
   source fact, exactly as the cadence ratio and clamps already are -- not a
   domain shrink, because it is a property of the shipping code rather than a
   restriction on the sea state.
2. A slew-limited `tau` cannot hold the bit-exact 2-cycle, so zero-`S` words leave
   the legal set.
3. Every remaining legal word then has at least 19 `S` firings, so the `S`
   observability Gramian is nonsingular and **a `P0`-free `Sigma_upper` exists
   again**.
4. The whole-word lower feasibility result below then applies as stated:
   `delta ~ 3e-8`, clearing the `1e-18` gate by 10.5 orders.

This keeps P3 uniform in initial covariance and needs no new declared entrance
bound. It is preferred over the `P0`-bounded fallback recorded below, which
remains available and clears by 12.7 orders but weakens what P3 asserts.

**Next falsifiable experiment:** derive the per-stage `tau` slew bound from the
estimator horizon and tuner EMA, add it to the P2 transition relation, and re-run
the four-max witness and the starvation witness against it. **Predicted:** the
resonant 2-cycle becomes illegal while #476's `[9,3,39,9]` label witness stays
legal, since that one needs only ordinary cell-to-cell motion. **Falsified if**
the slew bound is loose enough to still admit the 2-cycle, in which case the
`P0`-bounded fallback is the remaining route.

## Why a P0-free certificate cannot be built on the current legal word set

`translation_upper` takes **no initial covariance**: it is `P0`-free by
construction, which is exactly why it inverts a three-point observability
Gramian and is four orders loose. The looseness buys `P0`-independence.

PR #476 certified, and this branch reproduced exactly, that **zero-`S` words are
legal**. The proof domain has no GNSS -- `position_error_norm_upper_m` and
`velocity_error_norm_upper_mps` are *entrance* bounds, not measurements -- so
over a proof word the only translation measurement is `S = 0`. On a zero-`S`
word, `v -> p -> S` is therefore a pure integrator chain **with no output at
all**: structurally unobservable, covariance unbounded above as `P0` grows.

**So no finite `P0`-free `Sigma_upper` exists over the current legal word set.**
No amount of tightening, subdivision, ordering or enclosure work produces one;
an unobservable state cannot be bounded without bounding where it started. This
is why the whole-word lower result below is conditional, and it is a structural
obstruction rather than a numerical one.

### Measured: a P0-bounded formulation clears the gate by 12.7 orders

Dropping `P0`-freedom and bounding the admissible initial covariance by the
domain's declared entrance set -- velocity 5.0 m/s, position component
`0.5*Hs <= 4.25` m, `||delta a_w|| <= 0.3 g` -- on the **starved** word itself:

| assumed `S` entrance | `delta` | vs `1e-18` |
| --- | --- | --- |
| 13.5 m.s (`p_component * T_word`) | 4.86e-6 | clears by **12.69 orders** |
| 63.5 m.s (`p_norm * T_word`, looser) | 1.16e-6 | clears by 12.06 orders |

The worst single-cell dwell gives 4.10e-6, so starvation costs essentially
nothing once `P0` is bounded. **The obstruction is the `P0`-freedom requirement,
not the starvation.**

The domain declares no entrance bound on `S`, the third integral, which is the
channel that binds `delta`. That bound is **not load-bearing** and was not
invented to make the proof close: the gate survives until an assumed `S`
entrance of about `1e8 m.s`, against a physically meaningful 10--60 m.s, so
there are six-plus orders of slack in the assumption itself.

| assumed `S` entrance std | `delta` | |
| --- | --- | --- |
| 1e3 m.s | 5.87e-9 | clears by 9.77 orders |
| 1e5 m.s | 5.87e-13 | clears by 5.77 orders |
| 1e7 m.s | 5.87e-17 | clears by 1.77 orders |
| 1e9 m.s | 5.87e-21 | **fails** by 2.23 orders |

### What this costs

A `P0`-bounded `Sigma_upper` is **weaker** than the present one: the theorem
becomes conditional on the declared entrance set rather than uniform in initial
covariance, and that must be stated wherever the result is claimed. It is not
ad hoc -- P1 already hands off through that entrance box and P5 must prove
capture into it -- but it is a real change to what P3 asserts, and it is a
specification decision rather than something a proof tool may adopt on its own.

Not a certificate: double-precision point diagnostic on two words, not interval
arithmetic, and it does not quantify over all legal words.

## The whole-word lower result is conditional on S firing, which is not a theorem

**PR #476 is ahead of this line of work and its result supersedes the headline
below.** #476 certified that alternating two applied `tau` values *inside one
source cell* drives the `S` firing count to **zero** across a whole 635-sample
word. Reproduced here exactly, using the shipping binary32 scheduler:

| word at node 720's cell, gap 13 | S firings |
| --- | --- |
| `tau` held at 9.533334732055664 | 24 |
| `tau` held at 9.533475875854492 | 24 |
| **alternating the two, same cell** | **0** |

The two periods are 0.1300000101327896 and 0.13000193238258362, differing by
1.92e-6 against a binary32 scheduler tolerance of 1.91e-6. Each stage boundary
rebases the timer by `fmod` against the other period and the remainder never
crosses.

A word with no `S` firing never forgets its initial covariance. Measured on that
word: `P0_independent = False` with a probe spread of 0.999, so the reported
ratio 3.11e-17 is a finite-probe artifact and the true value tends to zero as
`P0` grows. **No `P0`-independent whole-word lower exists for that word**, so the
dwell figures below hold only under the assumption that `S` fires at all -- an
assumption #476 showed is not implied by the current P2 quotient.

**Two defects this exposed in the diagnostic here, both fixed.** The scheduler
was transcribed in double precision, but the deployed MEKF is
`Kalman3D_Wave_OU_III<float>`; the `periodic_update_due` tolerance is 1.9e-6 in
binary32 against 3.6e-15 in double, and a double transcription cannot represent
the starvation at all. And `PSEUDO_RATIO` was hand-written as `0.015/1.1`, while
the C++ constant is a `constexpr float` quotient
`f32(f32(0.015)/f32(1.1)) = 0.013636362738907337`; the two differ in the eighth
digit, moving the period by 1.5e-8 s, which is decisive at this tolerance. The
schedule constants are now parsed from the shipping header instead of restated.

## Whole-word lower feasibility, conditional on S firing (measured)

`tools/ou3_p3_whole_word_lower_feasibility.py`, non-promoting point diagnostic.
Canonical `delta` is same-history matched -- `phase_row` pairs node `t`'s floor
image with `endpoint_phase_upper(t, ...)`, that same node's own upper, and then
minimises over nodes. The analogue measured is, per source configuration, the
horizon-matched ratio between the whole word propagated from `P = 0` and the same
word propagated from a large initial covariance.

Over all 800 physical source nodes at their adverse corner (slowest cadence,
largest `sigma`, weakest `S`), holding each for the full 635-sample word:

| | value |
| --- | --- |
| worst single-cell dwell | **3.51e-8** at node 729 (`tau = 12 s`, `sigma = 0.05`, `R_S = 400`, 19 firings, binds on `p`) |
| worst *mixed* legal word found | **2.97e-8** (greedy local search, 1.184x below the dwell) |
| **S-starved legal word (#476)** | **no P0-independent lower exists** -- the word never forgets `P0` |
| clears the `1e-18` gate by | **10.47 orders** |
| independent mpmath `dps=30` sweep | 3.50347e-08, agreeing to 0.24 percent |
| best node | about 1.0 (node 390) |

**Dwelling is not the adversary -- hypothesis falsified.** This ledger predicted
that mixing could only help a lower, because a mixed word necessarily visits
faster-cadence cells. False. Greedy single-swap descent from the pure dwell at
node 729 converged after four improving swaps:

| round | swap | ratio |
| --- | --- | --- |
| 0 | pure dwell at 729 | 3.51191e-08 |
| 1 | position 18 to node 649 | 3.14549e-08 |
| 2 | position 17 to node 569 | 3.12656e-08 |
| 3 | position 11 to node 649 | 2.96665e-08 |
| 4 | position 24 to node 798 | 2.96635e-08 |
| 5 | no legal single swap improves | converged |

The worst legal word found uses four sources -- 569, 649, 729, 798, all high-`tau`
high-`R_S` cells -- and is **1.184x** below the pure dwell, still clearing the
gate by **10.47 orders**. So the recorded prediction (mixed words stay within
about two orders, above `1e-10`) **survives**, while the stronger claim that the
dwell *is* the adversary does not. Treat the dwell as a proxy good to about 20
percent, not as the worst case.

Two limits on that search. It is greedy single-swap and returns a **local**
minimum, not the global worst legal word. And the legality constraint is far
tighter than the self-loop result suggests: every one of the 800 nodes self-loops
on all 14 gaps, so every dwell is a legal word, but only **75 of 3995**
single-segment substitutions around a 729 dwell are legal, because `729 -> cand
-> 729` usually is not.

On #476's own adversarial witness extended to the word, the ordered ratio is
2.00e-3 at the slowest-cadence corner and 3.88e-3 at the fastest, so the choice
of `tau` within a cell is worth about 2x, not orders. Ordering *helps*: the worst
fixed-source cell in that witness gives 3.5e-8 while the mixed word gives 2.0e-3,
because a mixed word necessarily visits faster-cadence cells that fire `S` more
often. The adversary for a lower is therefore dwelling, not mixing.

**Scheduler correction to this ledger's own earlier numbers.** The cadence is
tau-scaled *and* `set_pseudo_update_period_s` applies
`elapsed = fmod(elapsed, new_period)` on every source commit, so the fixed-cell
counts recorded above (19 at `tau = 12 s`, 602 at `tau = 0.333 s`) do not
transfer to a changing-source word. The witness word fires **106** times at the
fastest corner and **76** at the slowest -- neither equals any single cell's
count. Fixed-cell counts remain valid only for a word spent entirely in one cell.

**Numerical limits, surfaced rather than tuned away.** `P0 = infinity` is not
representable in double precision, and two walls bracket the usable range: below
the lower wall the word has not forgotten `P0`; above the upper wall the
covariance-form update loses the endpoint to cancellation, the tell being that
the endpoint *decreases* as `P0` grows, which the monotone Riccati map forbids.
The gap between the walls moves across cells by about eight orders -- the
accepted probe scale ranges 1e6 to 1e15 over the 800 nodes -- so no fixed probe
pair serves all of them. The probe escalates relative to each word's own
zero-start covariance until two consecutive probes agree, and reports the
accepted scale and achieved spread per row.

Not a certificate: point diagnostic at one corner per cell, double precision,
`quantifies_over_legal_histories: false`, `certifies_theorem_stage: false`.

## Failure classification (historical, both blockers now cleared)

The canonical run was blocked by `Loewner prediction lower lost strict SPD; split
x cell`, and separately by the `x=0.01` branch-partition defect. Both are
resolved: the branch clamping was an implementation defect, and the strict-SPD
demand was surplus -- the theorem certifies the segment endpoint, not every 5 ms
prediction. Neither was the mathematical limiter. See the State section for the
one that is.

## Evidence

For source node 137 / gap 13, the first uncertain prediction on the widest initial small-`x` cell has diagnostic center `lambda_min ~= 3.405e-10`, absolute row-radius `eps ~= 1.142e-5`, and `eps/lambda_min ~= 3.35e4`. Existing binary subdivision improves that ratio by about 2x per level, but at depth 12 it is still about `8.26`; roughly four additional binary levels would be needed merely to make this first 5-ms prediction plausible. That is not a viable source-complete architecture.

The theorem-relevant 13-sample point covariance is much healthier. On source 137 its diagnostic `lambda_min` is about `2.02e-5 ... 2.83e-5` across the tau endpoints. Across the ten physical tau cells, using the actual low-sigma/strongest-`R_S` source in each cell, the 13-sample point floor remains strict; the longest-tau cell still has `lambda_min ~= 1.20e-6`.

Point diagnostics also show the 13-sample covariance increasing in Loewner order as `x=h/tau` increases on all ten tau cells: a 101-point grid produced positive minimum eigenvalue increments, with the smallest observed increment about `4.3e-9`. This is feasibility evidence only, not a proof of monotonicity.

Two segment-level interval experiments have now failed:

1. **Post-hoc congruence normalization of the old natural interval Riccati recursion.** Dependency has already exploded before normalization. Even with 256 tau subcells, the normalized common lower remains strongly non-SPD (worst diagnostic normalized `lambda_min - eps` about `-80`).
2. **Natural-interval derivative/monotonicity propagation.** Propagating `P(x)` and `dP/dx` as ordinary interval boxes inherits the same recursive dependency loss; subdivision improves widths but does not produce a usable derivative-SPD certificate.

A third, independent experiment confirms the same conclusion from the opposite
direction. `one_step` imposed three strict-SPD obligations per sample (after the
prediction and after each scalar measurement update), 39 per 13-sample subcell,
where the theorem makes one: P3 consumes the segment endpoint floor, certified
downstream by `common_boundary_floor` via `certified_rho(posterior) > 0`.
Removing the 38 surplus obligations collapses the adaptive split tree and drops
the segment module from 1191 s to 434 s, and the failures stop being
`split x cell` exceptions. But the propagated endpoint lower is then **not SPD**
on the same cells, while the point diagnostics above put the true 13-sample
`lambda_min` at `2.02e-5 ... 2.83e-5`. So the per-step SPD demand was never the
limiter, and the collapse architecture loses the entire floor between the true
value and the propagated lower.

Therefore the rejected mechanism is now broader than the original `C-eps I` step: **recursive natural interval covariance boxes that forget the repeated scalar source parameter at every Riccati operation are frozen as a dead end.**

Canonical P3 still has not reached translation/H/A margin calculation, so none of these diagnostic values are P3 theorem margins.

## Endpoint-only certification succeeds; the rejection was against the wrong quantity

The dead-end call on recursive interval covariance boxes was measured against
the **intermediate per-step** SPD demand, where the first uncertain prediction
has `lambda_min ~= 3.4e-10` and about four more binary levels would be needed.
The theorem does not impose that demand; `common_boundary_floor` certifies the
**segment endpoint**, where `lambda_min ~= 2e-5`, about 5.9e4 times larger.
Against the endpoint quantity the existing machinery succeeds.

Measured at 64 uniform x subcells with the intermediate gates removed, 18 of 18
tested (node, gap) pairs certify a strictly positive endpoint floor:

| node (tau,sigma,R_S) | gap 13 | gap 20 | gap 26 |
| --- | --- | --- | --- |
| 0 (0,0,0) | 64 / 9.50e-7 | 64 / 1.18e-5 | 64 / 2.00e-5 |
| 137 (1,5,7) | 64 / 8.74e-7 | 64 / 9.04e-6 | 64 / 1.58e-5 |
| 399 (4,7,9) | 64 / 7.75e-5 | 128 / 4.72e-4 | 128 / 4.93e-4 |
| 555 (6,7,5) | 64 / 8.43e-5 | 128 / 3.70e-4 | 128 / 4.05e-4 |
| 729 (9,0,9) | 64 / 8.69e-8 | 64 / 5.47e-7 | 128 / - |
| 799 (9,7,9) | 64 / 6.08e-5 | 64 / 7.14e-5 | 128 / 2.68e-4 |

Subdivision now has the scaling argument it previously lacked: halving the cell
width flipped 33/33 failing to 65/65 passing, and the endpoint margin improves
roughly linearly in width. This is a different quantity from the rejected one,
not a third refinement of it.

**The blocker is now cost, not feasibility.** `common_boundary_floor` is 800
nodes x 14 gaps = 11200 segment scans; at the observed 2 s to 33 s each the full
run extrapolates to about 103 h against a 120 minute CI budget.

The floor is provably monotone in the sigma and R_S cell index, so node
`(tau,0,0)` dominates its whole tau cell and the scan collapses to 10 nodes:

- `Q` scales as `sigma_lo^2` and the Riccati map is Loewner monotone in `Q`;
- `d/dR [P - Pe(e'Pe+R)^-1 e'P] = Pe(e'Pe+R)^-2 e'P >= 0`.

This is domination among actual reachable nodes within one tau cell, not the
Cartesian tau/sigma/R_S extrema rejected above. Measured: for tau=1 the sigma
sweep gives 8.74e-7 for indices 0..5 (identical, because the committed filter
sigma is clamped at the 0.05 floor and the code uses `sigma.lo`), 2.47e-5 at
index 6 and 6.75e-5 at index 7. Reducing 800 nodes to 10 brings the estimate to
roughly 25-77 min, which fits the budget only just, and the dominating nodes are
the slow ones.

## The S=0 pseudo-measurement is inert over a *segment* (framing superseded)

Retained as a measurement. The conclusion drawn from it at the time was wrong:
inertness over 65 ms says nothing about `R_S` over a word, where the filter
fires it 19 to 602 times depending on the cell. See the State section. That
horizon difference is itself the current research question.

Across a 23-decade override sweep at node 80 the S update is live and strong
when `R_S_z` is small -- at `1e-4` it collapses `P[2][2]` from 4.277 to 8.95e-5
and moves the floor from 2.83e-5 to 4.62e-6 -- but it saturates by `R_S_z ~ 1e6`
and the deployed value is `7.46e11`, six orders past saturation. The control
sweep confirms the accelerometer update is in its active region and is doing the
work. So this is a real regime fact, not a defect.

The reason is the zero start: after 13 samples `P[2][2] = 4.277` in `z_S = S/h^3`
units, i.e. an `S` standard deviation of `2.6e-7 m.s`, against a pseudo
measurement standard deviation of `0.108 m.s`, a factor of 4.2e5. A triple
integral accumulates little uncertainty in 65 ms.

Consequence to test, not yet a claim: P3 translation observability is the
four-`S` spread argument over a 0.765 s window inside the 3.17 s word, but
`common_boundary_floor` builds its floor from a zero start over a single 13-26
sample segment and therefore never sees that mechanism. The bound stays valid --
`rs.lo` is the strongest measurement in the cell, which is the correct choice for
a lower bound -- but the segment floor may be a lossy proxy for the true relative
injection margin, and that is a candidate explanation for `delta` landing near
1e-18.

## Measured: the margin deficit is about 9 orders (cause since identified)

The measurements below stand. The open question they left -- what does explain
the deficit -- is answered by the covariance ceiling; see the State section.

`delta` was computed directly against the real `Sigma_upper` from
`HIST.endpoint_phase_upper`, at phase 0, for three source nodes. The target is
identical across them:

`Sigma_upper (v,p,S,a_w) = [1.119e9, 3.970e9, 7.179e9, 1.109e3]`

| node | delta with `rho*I` (global rho 8.69e-8) | delta with the true anisotropic segment floor | gain | binding group |
| --- | --- | --- | --- | --- |
| 0 | 1.890285e-31 | 1.610690e-27 | 8.5e3 | S |
| 137 | 1.890285e-31 | 1.518656e-27 | 8.0e3 | S |
| 729 | 1.890285e-31 | 1.607829e-28 | 8.5e2 | S |

The canonical gate is `1e-18`. The current isotropic floor is **13 orders**
short; repairing it to the full anisotropic floor still leaves **9 orders**.

Two consequences.

First, the isotropic `rho*I` collapse is a real defect worth about 8.5e3, but it
is not the limiter. An earlier estimate in this ledger put the gain at ~5e6 from
the anisotropy ratio `lambda_max/lambda_min ~ 4e6`; that was wrong. `delta`
comes from a matrix PSD binary search, not the diagonal bracket: for node 0 the
S per-group ratio `d_i^2 Pz_ii/upper_i` is 9.2e-24 while the certified `delta`
is 1.6e-27, about 5700 times lower.

Second, and more important: **no enclosure or scalarization work can close a 9
order gap.** Taylor models, tighter collapses, finer subdivision and the
sigma/R_S cost reduction are all worth factors of 10^3 or less against a deficit
of 10^9. They should not be pursued as a route to the gate. The cause was
subsequently identified as a horizon mismatch between the segment floor and the
whole-word ceiling; the ceiling itself is 5.06 orders loose. See the State
section.

### Where the deficit plausibly comes from

The S group binds everywhere, and its ratio is a tiny floor over a very large
ceiling. `Sigma_upper` for S is `7.18e9 (m.s)^2`, a standard deviation near
85000 m.s, against a declared entrance box of 300 m.s -- the ceiling is the
same-history covariance propagated over the word, some 8e4 times the entrance.
The floor is `h^6 * Pz[2][2] = 1.5625e-14 * 4.25 = 6.6e-14`.

`common_boundary_floor` discards all older process covariance and restarts every
segment from `P = 0`, which is what makes the segment bound easy to state. For a
triple-integrated OU state the accumulated floor grows steeply in horizon; over
the 635 sample word versus a 13 sample segment that is roughly
`(635/13)^7 ~ 1.4e12`, the right order to cover a 10^9 deficit. This is the same
conservatism that makes the S=0 pseudo-measurement inert over a segment.

Hypothesis to test before any further enclosure work: the zero-start
per-segment floor, not its enclosure, is what costs the margin. The falsifiable
version is to propagate a floor across a horizon long enough to contain the
four-`S` spread and recompute `delta` against the same `Sigma_upper`. If
`delta` moves by orders, the segment formulation is the limiter and the fix is
structural. If it does not, the theorem formulation itself is marginal and
should be reconsidered rather than re-enclosed.

## Physical-reality check on the delta comparison target

`Sigma_upper` is the ceiling `delta` is measured against. In physical units it
is far outside anything a marine IMU in ocean waves can reach:

| channel | std(Sigma_upper) | declared entrance box | ratio |
| --- | --- | --- | --- |
| v | 33459 m/s | 5 m/s | 6691x |
| p | 63000 m | 20 m | 3150x |
| S | 84730 m.s | 300 m.s | 282x |
| a_w | 33.3 m/s^2 (3.4 g) | 2.94 m/s^2 | 11x |

A covariance ceiling asserting 33 km/s of velocity uncertainty and 63 km of
position uncertainty describes a diverging filter, not this one. It is a
bounding looseness, not a physical claim.

Substituting the declared physical box as the ceiling moves the binding S ratio
from 9.25e-24 to 7.38e-19 and the matrix `delta` from 1.61e-27 to about
1.28e-22. So the 9 order deficit decomposes as roughly **5 orders from the loose
ceiling and 4 orders structural**, the latter being the segment/word horizon
mismatch: a 65 ms zero-start segment accumulates 2.6e-7 m.s of `S` uncertainty,
which against any realistic `S` ceiling is intrinsically near 1e-18 at best.

`Sigma_upper` may not simply be replaced by the entrance box -- `delta` is a
relative Riccati injection margin and the ceiling must stay a valid upper bound
on the filter covariance, not on the error. But tightening that bound so it
reflects covariance the filter can actually reach is worth about 5 orders,
which exceeds every enclosure avenue combined and is not currently being
pursued by any recorded plan.

Conditions that are already correctly excluded, so they are not the problem:
impacts and slams out of normal-Live scope, lever arm disabled, vibration guard
dormant and transparent, non-gravitational CoG acceleration 4.0 m/s^2 = 0.41 g.

Entrance realism: the 45 deg attitude entrance is **defensible, not inflated**.
The filter initialises with a flat-sea estimate while the vessel is already
under way in waves, so the initial attitude error is the true wave-induced
attitude, not a small number. A 30 deg roll in heavy seas against a zero
estimate is already a 30 deg tilt error before heading error is counted, and the
domain separately declares a 10 deg internal heading gauge error. An earlier
revision of this file claimed a vessel starts within about 1 deg of tilt; that
answered the wrong question -- it described the sea being flat rather than the
filter's initial estimate being flat -- and is withdrawn.

## Whole-chain physical-realism audit

Checked every declared bound against a marine IMU in ocean waves, and traced
what each costs. `Sigma_upper` is reproduced exactly by
`BASE.translation_upper` at the full declared tuner ranges, so that function is
the source.

**Correctly excluded already, not the problem:** impacts and slams out of
normal-Live scope, lever arm disabled, vibration guard dormant and transparent,
non-gravitational CoG acceleration 4.0 m/s^2 = 0.41 g.

**Attitude entrance is realistic.** The filter initialises with a flat-sea
estimate while the vessel is already under way in waves, so the 45 deg entrance
is the true wave-induced attitude error, not an inflated margin. The canonical
P4 gate's 0.8 rad requirement binds a candidate to the outer sector prerequisite
it consumes; the domain's `[30, 25, 20, 15]` list is the inner attitude cell and
is explicitly separate (`p4_outer_geometry_sector_remains_separate: true`).
There is no conflict between them.

**Rotation rate is realistic.** `body_rate_norm_upper_deg_s = 30`. A 20 degree
roll amplitude at a 5 s period peaks at `20 * 2pi / 5 = 25 deg/s`, so the bound
sits just above genuine heavy-sea motion rather than being inflated. (An earlier
version of this line used a 30 degree roll at 6 s, which peaks at 31.4 deg/s and
so exceeds the declared bound; that example was wrong, not the bound.)

**The tuner cells are unphysical but nearly free.** `sigma_aw` reaches
6.0 m/s^2 = 0.61 g and `tau` falls to 0.333 s, i.e. 0.61 g of process noise
decorrelating in a third of a second. That is vibration, which the domain
excludes elsewhere, and it exceeds the domain's own 0.41 g acceleration
envelope. Cost of restricting to a physics-consistent `tau >= 3 s`,
`sigma <= 4 m/s^2`: `a_w` improves 4.04x, and v, p, S improve by 1.002x. So the
unphysical corner is not what inflates the translation ceiling.

**`R_S` upper drives the ceiling, but tightening it is not enough.**
`Sigma_upper` scales as `R_S,max^2` through the four-`S` inversion:

| rs.hi | Sigma_upper v | std |
| --- | --- | --- |
| 0.15 | 2.086e6 | 1444 m/s |
| 100 | 7.190e7 | 8480 m/s |
| 400 (declared) | 1.119e9 | 33459 m/s |

The declared 400 costs 538x on v, p and S. But even the **lowest declared `R_S`
upper** (0.15; no physical lower bound is established here) leaves a 1444 m/s
velocity ceiling *in `translation_upper`* -- a statement about the three-point
proof ceiling, not about the filter, whose settled velocity error is millimetres
per second. **The ill-conditioning is intrinsic to inverting the triple
integrator v -> p -> S from three points over a roughly 0.8 s observation
span**, not a consequence of any tuner or domain choice.

### Deficit accounting on the binding S channel

| lever | worth | delta after |
| --- | --- | --- |
| current | - | 1.89e-31 |
| repair the isotropic `rho*I` collapse | 8.5e3 | 1.6e-27 |
| tightest `R_S` | 538 | ~8.7e-25 |
| both | 4.6e6 | ~4e-25 |
| **canonical gate** | | **1e-18** |

Every realism lever pulled together still leaves about **7 orders**. Physical
realism accounts for roughly 3 of the 9 order deficit; the remainder is
formulation. **This conclusion is withdrawn**: it was measured against the
three-point-inversion ceiling, which measurement now puts at 5.06 orders loose,
so it is not evidence about a correct ceiling. Re-evaluate domain realism only
after the floor and ceiling refer to the same horizon.

The two structural items that remain are the intrinsically ill-conditioned
four-`S` triple-integrator observability route, and the zero-start 65 ms segment
floor compared against a whole-word ceiling.

## The certified chain cannot detect a state-correction sign error

Checked in response to a direct challenge on the third-integral correction
sign. The sign is **correct** in both places:

* shipping `applyIntegralZeroPseudoMeas` sets `r = -xext.segment<3>(off_S)`,
  i.e. innovation toward the target `S = 0`, and `gain_from_ldlt3_` computes
  `K = P C' S^-1` with no sign flip, so `dS = K_S r ~ -P_SS S_mat^-1 S` is a
  negative multiple of `S` and pulls it to zero;
* the proof's `_transition` matches the shipping `Phi` exactly once mapped into
  `z = [v/h, p/h^2, S/h^3, a_w]`: rows `[1,0,0,k1]`, `[1,1,0,k2]`,
  `[0.5,1,1,k3]` with `k2 -> 1/2` and `k3 -> 1/6` as `x -> 0`, the triple
  integrator Taylor coefficients. All integration couplings are positive; the
  only negative entries are structural zeros outward-rounded to `-5e-324`.

**But the check that found this is one the proof chain cannot perform.** The
covariance update `P - P C'(C P C' + R)^-1 C P` does not reference the
innovation `r` at all. Had `r` carried the wrong sign, the filter would drive
`S` away from zero and diverge, while every P3 certificate passed identically:
the segment floor, `Sigma_upper` and `delta` would all be unchanged.

The source-marker contracts do not close the gap either. They pin the update
*shape* -- `xext.noalias() += K * r;` and `joseph_update3_(K, S_mat, PCt);` --
but nothing anywhere checks how `r` is defined.

P4 would catch it, since it bounds the actual nonlinear state map `F_w`. P4 has
never been reached because P3 blocks. So the certified chain does not verify
that the mean correction is stabilizing; P1-P3 certify covariance behaviour only,
and any claim that the deployed filter is certified should state that limitation.

**Closed by direct test**, independently of the P3/P4 work:
`tests/kalman_ou_iii/innovation_sign-test.cpp` pins that every shipping OU-III
measurement update corrects the state *toward* its measurement. The invariant is
sign-sensitive and gain-free: a correct linear Kalman update makes the posterior
a convex combination of prior and measurement, so the corrected component lies in
`[x-, z)` and can never move away from `z`. Covered: the `S = 0` integral
pseudo-measurement, the position, velocity and vertical-velocity pseudo-updates,
and the accelerometer attitude correction.

The test was validated by negative control -- flipping the innovation sign in
each of the four linear updates in turn, rebuilding, and confirming each flip is
caught. Note the test Makefiles track the `.cpp` only, so a header edit needs
`make -B`; without it the first negative-control pass silently ran a stale binary
and reported four false passes. This does not affect CI, which builds clean.

This closes the *sign* gap only. It does not make P1-P3 certify the mean, and it
is not a substitute for P4.

## The translation ceiling credits 3 of the 19 to 602 S firings in a word

Raised as a challenge: `R_S` is the filter's major stabilizing force by design,
so why does this ledger read as though it degrades the bound? The challenge is
right and the earlier framing was wrong.

| | S observations per word |
| --- | --- |
| proof covariance ceiling (`integrator_inverse`, 3x3, "three possibly correlated selected observations") | **3** |
| proof observability route (`FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO`, `aligned_firing_count`) | **4** |
| shipping filter at tau = 12 s, `T_S = 0.164 s` (the cell that **binds** the ceiling) | **19** |
| shipping filter at tau = 1.1 s, `T_S = 0.015 s` | 201 |
| shipping filter at tau = 0.333 s, `T_S = 0.005 s` (clamped floor) | 602 |

The cadence is tau-scaled and the word length depends on the same tau, so the
count is `floor(T_word / T_S)` evaluated **within one cell**; it is never valid
to take the word from one cell and the cadence from another.

`translation_upper` bounds the translation covariance by inverting a
three-point integrator observation. The filter applies the `S = 0`
pseudo-measurement `floor(T_word / T_S)` times across the word -- 19 at
tau = 12 s, 602 at tau = 0.333 s. A count ratio is not an information metric, so
this ledger states the consequence in covariance instead: at matched cells, a
ceiling that credits every firing is 1.08e3 to 3.31e8 times smaller in variance
than the three-point one, and 1.70e4 times smaller at the cell that binds.

Two earlier statements in this ledger were true of that model and false as
statements about the filter, and are withdrawn as framing:

* "the `S=0` pseudo-measurement is inert" -- true over a 65 ms segment carrying
  0 to 4 firings, and says nothing about `R_S` over a word;
* "even the tightest `R_S` leaves a 1444 m/s ceiling" -- that is the ceiling of
  the three-point inversion, not of the filter.

Widening the observation baseline does not rescue it. Sweeping the spread:

| Tpe | spacing | Tword | Sigma_upper v | v std |
| --- | --- | --- | --- | --- |
| 1.0 | 1.00 | 3.17 | 1.119e9 | 33459 m/s |
| 2.0 | 2.00 | 6.17 | 1.669e7 | **4085 m/s** |
| 4.0 | 4.00 | 12.17 | 2.836e7 | 5325 m/s |
| 8.0 | 8.00 | 24.17 | 5.492e8 | 23440 m/s |

The optimum near `Tpe = 2` is 67x better than the deployed configuration but
still three orders beyond any boat, because a longer word accumulates more
process noise. **The limiter is the firing count, not the baseline.**

This relocates the 9 order deficit again, and this time onto something that is
a modelling choice rather than a property of the filter. It is not the
enclosure, not the isotropic collapse (8.5e3), and not domain realism
(about 3 orders). Crediting the S firings the filter actually performs is the
first thing to try, and nothing else in this ledger should be attempted before
it.

## Master P3 quantity

The new backend is relevant only if it produces the complete-segment matrix lower used by canonical P3:

`P_segment(x) >= L_segment > 0`.

That lower enters the existing translation gate directly through

`D_h L_segment D_h^T - delta * Sigma_upper > 0`,

which is what `_certified_delta` tests before the H/A precision join. Improving unrelated one-step quantities is not useful unless it sharpens this complete-segment comparison.

## What the failures invalidate

- A useful absolute point Loewner lower after every uncertain 5-ms prediction.
- Blindly increasing `MAX_ADAPTIVE_X_DEPTH`.
- Post-hoc normalization of a covariance interval after natural Riccati dependency has already exploded.
- Natural-interval `P,dP/dx` recursion as a monotonicity proof engine.

## What the failures do not invalidate

- The `1e-18` canonical P3 usefulness threshold.
- The same-history P2-V1 source language.
- Existence of a rigorous complete 13--26-sample translation floor.
- H=18/A=21 full-state joining, which has not yet received a translation margin from this backend.
- P4 feasibility or infeasibility, because P4 remains blocked by canonical P3.

## Critic pass and alternatives

Assume the interval covariance recursion architecture is wrong. Its strongest defect is loss of the fact that the entire segment depends on **one repeated scalar** `x=h/tau`. Ordinary interval arithmetic replaces that one-dimensional curve by a new Cartesian matrix box after every operation.

Qualitatively different alternatives:

1. **Univariate centered Taylor-model / polynomial enclosure of the complete 13-sample Riccati map.** Preserve the same scalar symbol through all prediction and measurement operations and bound only a final remainder.
2. **Verified polynomial collocation/Bernstein or Chebyshev enclosure of each complete-segment matrix entry.** Again preserve one-dimensional dependence rather than recursively hulling matrices.
3. **Analytic complete-segment Gramian/Riccati lower in a different information/covariance representation.** Avoid intervalizing the covariance recursion itself. Note before attempting this: the natural candidate `P_k >= (G_c^-1 + G_o)^-1`, with `G_c` the reachability and `G_o` the observability Gramian, is **false** for segments started from `P_0 = 0` -- 223/352 violations on integrator-chain systems. The information-form derivation needs a finite `Y_0 = P_0^-1`, which a zero-start segment does not provide. A correct analytic lower has to come from the closed-loop reachability sum, whose transitions depend on the gains.
4. **Different Lyapunov representation** if the complete-segment covariance map remains unsuitable.

**Selection (historical, superseded -- do not execute):** pursue (1) once. It directly attacks the failed dependency mechanism, and the point map is smooth with generalized source-137 ratio `P(x)` versus the low-`x` endpoint close to `1 ... 1.08`, so a centered model should be proving a relative statement near unity rather than recovering five orders of lost SPD margin. If this Taylor-model attempt fails after one mathematically motivated refinement, freeze it and move to (3), not another interval subdivision variant.

## DEAD_ENDS

- **REJECTED: endpoint-only 800-node P2 ancestry.** It loses staged/committed path memory.
- **REJECTED: independent Cartesian `tau/sigma/R_S` extrema.** They destroy source-history correlation.
- **REJECTED: recursive absolute entrywise Loewner point lower plus deeper subdivision.** Tractable depth is insufficient.
- **REJECTED: one-step congruence normalization as the primary architecture.** It improves conditioning but still solves an unnecessary one-step property.
- **REJECTED: post-hoc complete-segment normalization of natural interval Riccati boxes.** Dependency has already exploded before normalization.
- **REJECTED: natural-interval derivative/monotonicity recursion.** It shares the same dependency loss.
- **REJECTED for #471: additional P4 micro-certificates before complete-word feasibility.**
- **REJECTED: `translation_upper`'s three-point integrator inversion as the
  covariance ceiling.** It credits three `S` observations where the filter
  applies 19 to 602 per word, and at matched cells is 1.08e3 to 3.31e8 looser in
  variance (1.70e4 at the binding cell), worth 5.06 orders on `delta`. Real, but
  not the limiter -- see the State section.
- **SHELVED, not rejected:** univariate Taylor model, tighter Loewner
  enclosures, deeper subdivision, the isotropic `rho*I` repair (worth 8.5e3) and
  the sigma/R_S cost reduction (worth 80x). All are sound but each is worth 10^4
  or less against a 10^12 ceiling error. Revisit only after the ceiling is
  fixed and a real `delta` exists.
- **PARKED: rigorous H=18/A=21 complete-word dissipation producer.** Written and
  unit-tested at commit `4d68493` (`tools/ou3_p4_complete_word_dissipation.py`),
  then removed: a rigorous certificate before the non-promoting `rho_w`
  diagnostic is the wrong order. Resurrect only after the diagnostic reports
  `rho_w` clearly below 1.

## Verified theorem steps

Randomised verification against `ou3_interval.symmetric_positive_definite_ldlt`.
All hold, so none of these are suspect as sources of the P3 blockage.

| Step | Result |
| --- | --- |
| `\|\|(R-I)v\|\|^2 = 4/(4+q^2) \|\|[c]x v\|\|^2` | exact, max err 4.7e-13 |
| `\|\|eta\|\| = sin(theta/2) \|\|h\|\|` | exact, max err 6.3e-13 |
| `\|\|R-I-[c]x\|\| <= (3/4) q^2` | 0 violations |
| `J <= n blockdiag(J_ii)` for PSD J | 0 violations |
| `K_theta S K_theta' <= P_theta_theta` | 0/1200 |
| `Phi_s' Sigma_s^-1 Phi_s <= Sigma_0^-1` | 0/1200 |
| `(1-3d/8)^2 <= 1-d/2` | holds for all `d <= 16/9` |

Prefix nonexpansiveness needs `Phi_s` invertible. It is, via the Joseph form and
invertible OU prediction, but the source-path document does not state it.

## Available relaxations, not yet warranted

These enlarge the certified P4 funnel `W_*` only; none changes
`rho = 1 - delta/2`, which is set entirely by the P3 margin. Do not spend
effort on them before a real `delta` exists.

- `lambda_max(Sigma) <= sum_g U_g` instead of `n_g max_g U_g`. With `U_S ~ 9e4`
  dominating this is 9.04e4 against 6.3e5, about 7x.
- Attitude corrections need only the attitude marginal,
  `|dtheta| <= sqrt(U_theta/lambda_min(R)) |y|`, so `sqrt(0.25)` replaces
  `sqrt(9e4)` -- about 600x on the correction gain.
- The sharp Cayley remainder `q^2(q+2)/(4+q^2)` is uniformly `0.805x` the
  `(3/4)q^2` bound.
- The `4` in `B_m = 4 N_op sqrt(m_+) C / m_-` comes from a `W_s <= 4 W_0`
  bootstrap; prefix nonexpansiveness gives `W_s <= (1+d/8)^2 W_0 <= 1.27 W_0`.

A rejected route may be resurrected only after recording the new mathematical fact that invalidates its rejection.

## Retained facts

- The scalar covariance-form measurement update `P - Pe(e'Pe+R)^-1e'P` requires
  only `e'Pe+R > 0`, never strict SPD, and is Loewner monotone on **all**
  symmetric arguments with that denominator positive. For `H >= 0`,
  `x'(dU)x = x'Hx - 2(x'He)(x'Pe)/d + (x'Pe)^2(e'He)/d^2` is a quadratic in
  `x'Pe` whose discriminant is `4((x'He)^2 - (x'Hx)(e'He))/d^2 <= 0` by
  Cauchy-Schwarz in the `H` semi-inner product, so `dU >= 0`; `d` is affine in
  `P` so it stays positive between ordered arguments. Any backend may therefore
  carry singular or indefinite intermediates without invalidating the order.
- Canonical P3 useful gate remains exactly `1e-18`.
- `OU3_P2_CORRELATED_STAGE_TRANSFER_V1` and same-history source correlation are retained.
- Both H=18 and A=21 remain required.
- Zero/disabled lever arm and dormant/transparent vibration-guard branch remain the current proof scope.
- No replay fitting, operating-domain shrink for PASS, or deployed-filter change is permitted.
- Existing exact Joseph, co-rotated accelerometer, reset, and finite-angle identities remain structural facts only; they do not authorize P4 promotion.
- P4 is blocked until canonical P3 passes; P5 is blocked until canonical P4 strictly contracts.

## Next falsifiable experiment

**The single current experiment is stated in the State section** (make the floor
and the ceiling whole-word quantities, and credit the `S` cadence). It is not
repeated here, so that this ledger carries exactly one next action.

The previously selected univariate centered Taylor model, and the
ceiling-only Riccati bound that briefly replaced it, are both **shelved**. The
Taylor model addresses enclosure dependency, which measurement shows is not the
limiter; the ceiling-only bound was measured at 5.06 orders, which leaves 6.66.

Before any new P4 proof producer, first obtain the canonical P3 numeric verdict.
If P3 passes, the only P4 work allowed is the non-promoting high-precision
complete-word ratio diagnostic `rho_w = V_after(F_w(x)) / V_before(x)`.
