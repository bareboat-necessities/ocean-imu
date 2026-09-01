# OU-III P1--P5 usable certificate envelope

This note defines what **usable** means for the deployed OU-III stability
certificate.  A certificate is not useful merely because every strict
inequality is positive.  The certified domains must overlap the actual source
handoffs and operating ranges in physical units, and the proof must not obtain a
large number by weakening the implementation model or the theorem claim.

The old P4/P5 perturbative route failed this criterion.  PR #441 proved a hard
attitude ceiling of about `1.24e-3 rad` at the shipping prefix factor, while the
source-faithful gauged P1 handoffs have Cayley norms about `0.272` (normal) and
`0.595` (timeout).  No further sharpening of the old scalar constants can bridge
that gap.  The route, not just its constants, must change.

PR #438 was therefore not wholly wasted, but it is not mergeable as a theorem
route.  Its direct/next-generation/third-generation/frontier scalar searches
still terminate in the retired whole-word scalar accounting.  The useful part
is the exact-rational **complete-word translation accumulation**.  This branch
reuses that core behind a repaired outward-rounded wrapper, then treats the
result only as one block of the still-missing complete full-state certificate.

## Stage semantics

| Stage | Usable meaning | Current route |
|---|---|---|
| **P1 startup/reset** | Source-faithful reset, normal handoff, and timeout handoff cover deployment-sized initial attitude and physical-coordinate errors. | Already physical rather than microscopic. Preserve these bounds; do not shrink P1 to make P5 easy. |
| **P2 source language/path** | Every admitted runtime branch is represented, but correlated adaptation states are propagated through shipping dynamics instead of independently re-choosing worst extrema on every sample. | `ou3_p4_source_path_reachability.py` models the deployed tuner EMA state. Raw tuner `sigma` is allowed below the separate `0.05 m/s^2` MEKF stationary-standard-deviation floor. EMA exponentials are validated and arbitrary late commits are over-approximated. The adaptive `R_S` slew horizon is covered. Until implementation `powf/sqrtf` error is separately enclosed, the `R_S` target uses the full deployed clamp instead of an unsafe tight target box. |
| **P3 linear word** | Global linear covariance/information statement over the declared source-parameter domain; no local state-amplitude ball is required. | Keep the validated source-reachable matrix certificate. Its small `delta` is a relative Riccati/noise comparison constant, **not a state radius** and must not be advertised as one. The useful PR #438 complete-word translation calculation is retained as a widened block result, not as a replacement for the full-state P3/P4 word. |
| **P4 nonlinear word** | The nonlinear geometry must cover a finite-angle domain that overlaps P1, and each nonlinear correction must be charged against its own information decrease with directional/block margins. | The operation-matched finite-angle sector is `0.80 rad` (`45.84 deg`) and contains both gauged P1 handoffs. Exact Cayley vector geometry, effective-vector reductions, Joseph identity, and exact reset congruence replace the old global `N x Lipschitz-defect` accounting. The retained PR #438 translation result is combined with the direct nontranslation margin only to emit a conservative Schur cross-block target; it is **not** promoted until the complete 18/21-state cross block is outward validated on the source graph. |
| **P5 startup capture** | P1 must enter a certified outer nonlinear domain without pretending ungauged yaw is observable. | The outer-sector certificate has `N_outer=0`: both gauged branches enter the full-SO(3) sector immediately; ungauged timeout enters the gravity-direction quotient sector with yaw treated as gauge. Its inclusion test compares the P1 gravity-cosine lower bound against the **upper** validated enclosure of the sector-boundary cosine. Finite capture from this outer sector to the final inner stochastic localization level remains downstream of complete-word P4 dissipation. |

## Non-regression rules

The P1--P5 stack must fail CI rather than regain a microscopic certificate by
any of the following shortcuts:

1. Do not shrink the P1 handoff or live operating domain to make a certificate
   pass.  A stronger theorem may widen the source domain; it may not silently
   delete previously admitted deployment states.
2. Do not treat the P3 relative Riccati margin as a nonlinear state radius.
3. Do not accumulate one global nonlinear defect once per packet and compare it
   with the weakest whole-word P3 direction.  PR #441 proves that route cannot
   reach the P1 handoff even with idealized constants.
4. Do not resurrect the PR #438 direct/next-generation/third-generation/frontier
   scalar optimizers as promotion inputs.  They may remain historical diagnostics
   but cannot establish a usable P4/P5 theorem.
5. Do not promote the PR #438 translation block, the minimum of translation and
   nontranslation block margins, or the Schur cross-block budget as a full-state
   certificate before the actual cross block has been outward enclosed.
6. Do not conflate the raw tuner `sigma_applied` state with the `0.05 m/s^2`
   filter-side stationary-standard-deviation floor.  The raw state can be below
   that floor after tuner variance readiness.
7. Do not gain source-path tightness by assuming a nominal-only `R_S` smoothing
   horizon, an upper bound on inter-commit delay that the implementation does
   not enforce, or unqualified `powf/sqrtf` implementation accuracy.
8. Do not charge the magnetometer radial finite-angle residual as a state
   disturbance when the shipping gain annihilates it exactly.
9. Do not charge the accelerometer finite-angle residual as independent
   measurement noise when `J_aw=R_wb` represents it exactly as an effective
   `a_w` input.
10. Do not multiply quaternion-reset transport by a covariance condition number;
    reset is an exact congruence and the proved inverse operator norm is one.
11. Do not assign a full-heading radius to the ungauged timeout branch.  Use the
    gravity/yaw quotient until a magnetic regauge event supplies a heading gauge.
12. Do not use the lower cosine enclosure to prove ungauged-sector inclusion.
    With a lower bound on the true cosine, the conservative comparison is against
    the **upper** enclosure of `cos(theta_sector)`.
13. Do not promote replay observations to a deployment theorem.  The widening
    remains source-generated and trajectory-independent.
14. Do not change filter gains, schedules, or estimator equations merely to make
    the proof pass.  A filter change is a separate design decision and requires
    new performance evidence.
15. Do not evaluate a gain row `m N / (m^2 p + lambda)` as an interval quotient.
    The specific-force magnitude `m` appears in the numerator and in the
    denominator, so the naive quotient overstates the accelerometer gain by up
    to the cell's magnitude ratio.  Use the exact unimodal supremum in
    `ou3_p4_shared_force_gain.py`.
16. Do not present the descending `30 -> 25 -> 20 -> 15` deg candidate ladder as
    the route that closes the first accelerometer operation.  The measured
    nuisance-over-budget ratio stays above one down to a 1 deg probe, so a
    narrower candidate cannot close that operation by itself.
17. Do not report `P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED` while the
    complete-word P4 dissipation or the finite P5 inner capture is open.  A
    stage with usable geometry and an open obligation is
    `USABLE_GEOMETRY_OPEN_OBLIGATION`, never `USABLE`.

## Quantitative outer target

The replacement outer nonlinear target is intentionally well above startup:

- design full-attitude sector: `0.80 rad` (`45.84 deg`);
- its Cayley boundary and two-sided cosine boundary are certified with the
  repository's validated transcendental enclosure and remain below `q=1`;
- exact vector strong-monotonicity factor is required to remain above `0.80`;
- exact nonlinear eta/residual information ratio is required to remain below
  `0.25`;
- normal and gauged-timeout P1 Cayley handoffs must both lie inside the sector;
- the ungauged timeout gravity-direction cosine must prove inclusion against the
  sector cosine **upper** enclosure without inventing a yaw bound.

These are minimum usability contracts, not tuned replay extrema.  Future proof
work may widen them further, but CI rejects regression below them.

## Useful result retained from PR #438

`ou3_p4_translation_full_word_rigorous.py` reuses the exact-rational complete-word
translation propagation from #438 while repairing the pieces that were not safe
to promote there:

- C++ `float` literals are interpreted as deployed binary32 values;
- conditioned scale products are intervalized from the first operation rather
  than formed in binary64 and wrapped afterward;
- the runtime-bounded dyadic Loewner compression is retained; and
- the endpoint is independently recertified.

`ou3_p4_post_translation_bottleneck.py` then emits, for each H/A mode, a
**lower-enclosed open budget**

`||C||_2 < sqrt(delta_translation * delta_nontranslation)`

for the missing normalized translation/nontranslation cross block.  The lower
enclosure is deliberate: rounding this threshold upward would make the future
full-state acceptance test fail-open.  The diagnostic recomputes all derived
fields during validation and remains `NOT_ESTABLISHED` until an actual cross
block is certified.

## Remaining theorem closure

The remaining hard theorem step is not another scalar-radius optimization.  It
is to propagate the complete 18/21-state word over every required source path,
outward-enclose the normalized full-state cross block below the emitted Schur
budget, and compose the operation-matched finite-angle correction dissipation on
those same cells.  That is the point at which the retained PR #438 linear
progress and the new finite-angle sector meet.

Once that outer complete-word decrease overlaps the existing inner stochastic
localization level, P5 can derive a finite word count from the P1/outer sector to
that inner level without the PR #441 uniform-transport obstruction.

## Measured first-accelerometer obstruction

`tools/ou3_p4_first_accel_sector_budget.py` measures why the complete-word P4
composition has not been emitted.  For each candidate angle it compares the
**budget** -- the largest correction norm whose worst-case signed Cayley
composition stays inside the `0.80` rad outer sector -- against the **nuisance
term**, the shared-force-magnitude accelerometer gain applied to the effective
`a_w` input.  The ratio falls from `5.05` at 30 deg to `2.12` at 15 deg, but the
non-candidate limit probes show it saturating at `1.34` for a 1 deg candidate.
The floor is source-faithful: the declared `0.3 g` latent-acceleration error
over the lowest admitted specific force is `2.941995 / 5.0 = 0.5884`, so the
accelerometer-implied gravity direction can be off by `0.6291` rad in the
low-force cells regardless of the candidate angle.

That producer reports a **distance and never a verdict**: the nuisance term is
an outward bound and the `a_w` covariance endpoints feeding the gain are
enclosure endpoints, not certified reachable points.  It does not retire the
ladder as a diagnostic; it retires the ladder as the closing route.

The complementary widening is `tools/ou3_p4_shared_force_gain.py`, which removes
the shared specific-force-magnitude dependency from the gain rows.  It is
uniformly between `1.17` and `1.78` tighter over the audited cells, and it moved
the signed 30 deg first-accelerometer stage from an abort after 13 of 40960
children -- with a composition denominator that could cross zero -- to a
complete evaluation of all 40960 with `max_d = 2.2251526515093487` rad and a
minimum denominator of `0.6356874937572935`.  The stage is still
`NOT_ESTABLISHED`, but for a measured reason rather than a lost enclosure.

`tools/ou3_p1_p5_certificate_numbers.py` collects every stage number and
re-applies the usability thresholds above to the recomputed values; see
`docs/ou-iii-p1-p5-certificate-numbers.md`.

Until that composition is emitted, the truthful status is:

- P1: source-faithful deployment-sized startup domain established;
- P2: source-complete conservative tuner path language established;
- P3: global linear source-word certificate established, with useful #438
  complete-word translation widening retained as partial block evidence;
- P4: physically useful finite-angle nonlinear sector established; conservative
  full-state cross-block target emitted; complete-word operation-matched
  full-state dissipation still to be composed;
- P5: source handoff into the useful outer sector established immediately;
  final finite capture to the inner stochastic level still downstream of P4.
