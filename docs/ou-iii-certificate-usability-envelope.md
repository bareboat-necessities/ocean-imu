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

## Stage semantics

| Stage | Usable meaning | Current route |
|---|---|---|
| **P1 startup/reset** | Source-faithful reset, normal handoff, and timeout handoff cover deployment-sized initial attitude and physical-coordinate errors. | Already physical rather than microscopic.  Preserve these bounds; do not shrink P1 to make P5 easy. |
| **P2 source language/path** | Every admitted runtime branch is represented, but correlated adaptation states are propagated through shipping dynamics instead of independently re-choosing worst extrema on every sample. | `ou3_p4_source_path_reachability.py` builds the deployed `tau/sigma_aw/R_S` EMA transition graph.  The historic weak corner remains explicit and its recurrence/residence is proved from source dynamics. |
| **P3 linear word** | Global linear covariance/information statement over the declared source-parameter domain; no local state-amplitude ball is required. | Keep the validated source-reachable matrix certificate.  Its small `delta` is a relative Riccati/noise comparison constant, **not a state radius** and must not be advertised as one. |
| **P4 nonlinear word** | The nonlinear geometry must cover a finite-angle domain that overlaps P1, and each nonlinear correction must be charged against its own information decrease with directional/block margins. | New operation-matched finite-angle sector is `0.80 rad` (`45.84 deg`) and contains both gauged P1 handoffs.  Exact Cayley vector geometry, effective-vector reductions, Joseph identity, and exact reset congruence replace the old global `N x Lipschitz-defect` accounting.  Complete-word path composition remains an explicit obligation. |
| **P5 startup capture** | P1 must enter a certified outer nonlinear domain without pretending ungauged yaw is observable. | New outer-sector certificate has `N_outer=0`: both gauged branches enter the full-SO(3) sector immediately; ungauged timeout enters the gravity-direction quotient sector with yaw treated as gauge and the gravity-parallel gyro-bias as a bounded neutral input.  Finite capture from this outer sector to the final inner stochastic localization level remains downstream of complete-word P4 dissipation. |

## Non-regression rules

The P1--P5 stack must fail CI rather than regain a microscopic certificate by
any of the following shortcuts:

1. Do not shrink the P1 handoff domain to fit the old P4 inner seed.
2. Do not treat the P3 relative Riccati margin as a nonlinear state radius.
3. Do not accumulate one global nonlinear defect once per packet and compare it
   with the weakest whole-word P3 direction.
4. Do not charge the magnetometer radial finite-angle residual as a state
   disturbance when the shipping gain annihilates it exactly.
5. Do not charge the accelerometer finite-angle residual as independent
   measurement noise when `J_aw=R_wb` represents it exactly as an effective
   `a_w` input.
6. Do not multiply quaternion-reset transport by a covariance condition number;
   reset is an exact congruence and the proved inverse operator norm is one.
7. Do not assign a full-heading radius to the ungauged timeout branch.  Use the
   gravity/yaw quotient until a magnetic regauge event supplies a heading gauge.
8. Do not promote replay observations to a deployment theorem.  The widening
   remains source-generated and trajectory-independent.
9. Do not change filter gains, schedules, or estimator equations merely to make
   the proof pass.  A filter change is a separate design decision and requires
   new performance evidence.

## Quantitative outer target

The replacement outer nonlinear target is intentionally well above startup:

- design full-attitude sector: `0.80 rad` (`45.84 deg`);
- its Cayley boundary is certified with the repository's validated
  transcendental enclosure and remains below `q=1`;
- exact vector strong-monotonicity factor is required to remain above `0.80`;
- exact nonlinear eta/residual information ratio is required to remain below
  `0.25`;
- normal and gauged-timeout P1 Cayley handoffs must both lie inside the sector;
- the ungauged timeout gravity-direction cosine must lie inside the same
  finite-angle gravity sector without inventing a yaw bound.

These are minimum usability contracts, not tuned replay extrema.  Future proof
work may widen them further, but CI should reject any regression below them.

## Remaining theorem closure

This change deliberately does **not** relabel the legacy microscopic P4 inner
level as a usable full certificate.  The remaining hard theorem step is to
compose the already-proved operation-level identities over every source-reachable
complete word, preserving directional information margins and the adaptation
path graph.  Once that outer-word decrease overlaps the existing inner
stochastic localization level, P5 can compute the finite word count from the
P1/outer sector to that inner level.

Until that final composition is emitted, the truthful status is:

- P1: source-faithful startup domain established;
- P2: source-dynamic path language established;
- P3: global linear source-word certificate established;
- P4: physically useful finite-angle nonlinear sector established, complete-word
  operation-matched dissipation still to be composed;
- P5: source handoff into the useful outer sector established immediately,
  final finite capture to the inner stochastic level still downstream of P4.
