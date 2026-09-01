# OU-III P1--P5 certificate numbers

This note answers two questions in one place: **what do the P1--P5 certificates
currently report**, and **are those numbers usable**?  Usable is defined in
`docs/ou-iii-certificate-usability-envelope.md`: the certified domains must
overlap the actual source handoffs and operating ranges in physical units, and
the proof must not obtain a large number by weakening the implementation model
or the theorem claim.

`tools/ou3_p1_p5_certificate_numbers.py` regenerates every number below and
re-applies each usability threshold to the recomputed value.  It does not treat
an upstream `PASS` field as a usability statement.  Run it with

```
python3 tools/ou3_p1_p5_certificate_numbers.py --output /tmp/ou3_p1_p5.json
```

and optionally `--first-accel` / `--translation` to fold in the two slow
producers without re-running them.

## Summary

| Stage | Status | Headline number | Usable? |
|---|---|---|---|
| **P1** startup/reset/handoff | PASS | gauged handoff Cayley norms `0.2721648148683776` (normal) and `0.5947333355555983` (timeout) | **yes** -- deployment-sized, not microscopic |
| **P2** source-word language | PASS | `800` source states, `640000` edges, one SCC, all states recurrent | **yes** -- source-complete, no replay |
| **P3** linear information word | PASS | H and A endpoint margins `2.2953997386276688e-20`, prefix information gain `<= 1.0` | **yes**, with the standing caveat that this margin is a *relative* Riccati/noise constant, not a state radius |
| **P4** nonlinear geometry | PASS | sector `0.80` rad (`45.836623610465864` deg), `q = 0.845586437476324 < 1`, monotonicity `0.8483533546735822`, eta ratio `0.17875410581097537` | geometry **yes**; complete-word dissipation **open** |
| **P5** startup capture | PASS | declared entrance `45` deg / `0.5 Hs`, `N_outer = 0`, every source branch enters immediately | outer entry **yes**; finite inner capture **open** |

So four of the five stages are physically sized rather than microscopic, and
none of them was obtained by shrinking a declared deployment domain.  What is
**not** established is the pair of theorem obligations that would turn the five
certificates into one stability proof:

1. complete 18/21-state source-correlated P4 word dissipation, and
2. the finite P5 word count from the 45 deg entrance to the inner stochastic
   localization level, which is downstream of (1).

The report therefore emits `P1_P5_COMPLETE_STABILITY_PROOF_ESTABLISHED = false`
and marks P4 and P5 as `USABLE_GEOMETRY_OPEN_OBLIGATION`.  It refuses to report
a complete proof while any obligation is open.

## Why the microscopic numbers are no longer the live ones

The retired uniform transported-defect route still reports an attitude capture
radius of `5.808479596010723e-32` rad, and its own ceiling -- independent of
every constant it contains -- is `1.23762376237624e-03` rad at the shipping
prefix factor.  Both are far below the `0.272` and `0.595` rad P1 handoffs, and
`tools/ou3_p4_p5_route_ceiling_certificate.py` proves no sharpening of that
accounting can close the gap.  Those numbers are retained as evidence that the
route is retired, not as the current certificate level.

The live P4/P5 geometry is the operation-matched finite-angle sector, which is
`0.80` rad and contains both gauged P1 handoffs directly.  That is the number to
quote for the current nonlinear domain.

## The remaining quantitative obstruction

The open P4 obligation now has a measured cause rather than a broken enclosure.
`tools/ou3_p4_first_accel_sector_budget.py` compares, for each P4 candidate
angle, two quantities on the same source/alignment/force children:

* the **budget** -- the largest correction norm whose worst-case signed Cayley
  composition with the post-prediction candidate norm still lands strictly
  inside the `0.80` rad outer sector;
* the **nuisance term** -- the shared-force-magnitude accelerometer gain applied
  to the effective `a_w` input (declared `0.3 g` startup latent-acceleration
  error, accelerometer bias, finite-angle force remainder).

| candidate | budget (rad) | nuisance (rad) | ratio |
|---|---|---|---|
| 15 deg | `0.5309502743982276` | `1.1270362643023522` | `2.12` |
| 20 deg | `0.4439819969818312` | `1.1909281418846758` | `2.68` |
| 25 deg | `0.3578622392054443` | `1.2735064196675971` | `3.56` |
| 30 deg | `0.27225152012902093` | `1.375114933176711` | `5.05` |
| 35 deg | `0.18677577048716126` | `1.4961948241502516` | `8.01` |
| 40 deg | `0.10101915914914625` | `1.6373388541150944` | `16.21` |
| 45 deg | `0.014514984270965614` | `1.7993361240597836` | `123.96` |

Shrinking the candidate angle does help -- the ratio falls from `5.05` at 30 deg
to `2.12` at 15 deg -- but it is bounded away from closing.  The producer's
non-candidate limit probes show the ratio saturating: at `1` deg the budget is
`0.7815378101554163` rad and the nuisance term is still `1.0456598581422538`
rad, a ratio of `1.34`.  The reason is visible in one source-faithful number:
the declared latent-acceleration error over the lowest admitted specific force
is `2.941995 / 5.0 = 0.5883990000000002`, so the accelerometer-implied gravity
direction can itself be off by `0.6291` rad in the low-force cells.  No choice
of candidate angle changes that.

**This is a distance, not a verdict.**  The nuisance term is an outward bound
and the `a_w` covariance endpoints feeding the gain are enclosure endpoints, not
certified reachable points, so a ratio above one does not prove that an
admissible state leaves the sector.  What it does establish is that the
`30 -> 25 -> 20 -> 15` deg search ladder cannot close the first deployed
accelerometer operation by angle alone, and that the two structural changes
already named by the route-ceiling certificate -- charging each defect against
its own operation's information decrease, and carrying a directional block
margin instead of one scalar whole-word `delta` -- are the ones that remain.

## Widening recorded with these numbers

`tools/ou3_p4_shared_force_gain.py` removes the shared specific-force-magnitude
dependency from the structured accelerometer gain rows.  Each row has the form
`m N / (m^2 p + lambda)` with `m` in both numerator and denominator; evaluating
it as an interval quotient overstates the gain by up to the cell's magnitude
ratio.  The rows are unimodal in `m`, so the exact supremum is the interior
value at `m* = sqrt(lambda/p)` -- `(1/2) sqrt(p/lambda)` for the self-`p` rows
and `C / (2 sqrt(lambda p))` for the independent-`C` axial row -- or the larger
endpoint.  The lemma reports a uniform tightening between `1.17` and `1.78` over
the audited cells, and is never looser.

The effect on the blocking stage is not marginal.  With the interval gain,
`ou3_p4_30deg_signed_first_accel_sector_v2` aborted after **13** of 40960
children with a correction norm of `2.6481953631664035` rad and a signed
composition denominator that reached `0.038994069387713985` and could cross
zero.  With the shared-force gain,
`ou3_p4_30deg_signed_first_accel_sector_v3` evaluates **all 40960** children:

| | v2 (interval gain) | v3 (shared-force gain) |
|---|---|---|
| children evaluated | `13` then abort | `40960`, complete |
| max correction norm | `2.6481953631664035` rad | `2.2251526515093487` rad |
| min composition denominator | `0.038994069387713985` | `0.6356874937572935` |
| max post-update `q` | not reachable | `6.583306237863824` |

The 30 deg family is therefore *measurable* for the first time.  It remains
`NOT_ESTABLISHED`, because `6.583` is still far outside the `0.8456` outer
sector -- which is exactly the gap the budget table above accounts for.

## Pricing the two remaining doors

`tools/ou3_p4_first_accel_aw_sigma_consistency.py` splits the nuisance term into
the part that is proof slack and the part that is a domain question, and puts a
number on each.

**Unconditional: joint force pairing.**  The nuisance is
`||K_theta|| * (a_w error + eta(q)|f| + b_a)`.  The gain falls roughly like
`1/|f|` while the finite-angle force remainder grows like `|f|`, so bounding the
two factors separately over one force cell charges the largest gain against the
largest force.  Maximising the product over subdivided magnitudes inside the
same cell tightens the nuisance term by `1.026` (15 deg) to `1.168` (45 deg)
with no domain, filter or candidate change.

**Conditional: the `a_w`/`sigma` pairing.**  The `a_w` covariance that sets the
gain is the tuner state, `P_aw = sigma_applied^2` with `sigma_applied` in the
deployed safety range `[0.05, 6.0] m/s^2`.  The `a_w` *error* is the separately
declared `0.3 g` startup envelope.  Nothing in the declared domain couples them,
so the worst cell pairs a tuner that believes the sea is flat with a
`2.9407694241234332 m/s^2` latent-acceleration error -- a ratio of `5.3691`
against that cell's `sigma` upper.  The shipping tuner does couple them
(`sigma_target = min(sigma_wave * sigma_coeff, max_sigma_a)`, with
`sigma_applied` an EMA toward it), but the EMA transient means no bound follows
from the update law alone.

The producer therefore assumes no coupling and instead measures the constant one
would have to supply.  With `||delta a_w|| <= c * sigma_applied` added to the
declared domain:

| candidate | budget (rad) | residual at zero `a_w` error (rad) | nuisance (rad) | ratio | required `c` |
|---|---|---|---|---|---|
| 15 deg | `0.5309502743982276` | `0.20587675702069888` | `1.0990056689852747` | `2.07` | **`1.9510667819413354`** |
| 20 deg | `0.4439819969818312` | `0.24791375638819854` | `1.1409879746372982` | `2.57` | **`1.176786821337373`** |
| 25 deg | `0.3578622392054443` | `0.30224527930294` | `1.195248807614594` | `3.34` | **`0.33380880686218`** |
| 30 deg | `0.27225152012902093` | `0.36909754878917767` | `1.2620140966060003` | `4.64` | none exists |
| 35 deg | `0.18677577048716126` | `0.44876080948633734` | `1.341573708617953` | `7.18` | none exists |
| 40 deg | `0.10101915914914625` | `0.5415985399383951` | `1.4342906492778942` | `14.20` | none exists |
| 45 deg | `0.014514984270965614` | `0.6480593087098117` | `1.5406129035243878` | `106.14` | none exists |

Two things follow, and they set the route.

First, **the 30 deg candidate cannot close the first accelerometer operation
under sector invariance at all.**  Its residual with a *perfect* `a_w` estimate,
`0.36909754878917767` rad -- accelerometer bias plus the finite-angle force
remainder alone -- already exceeds its `0.27225152012902093` rad budget.  No
consistency statement, no gain sharpening and no covariance enclosure can
recover that, for the same reason the route ceiling was independent of `kappa`.

Second, the candidates that can close need a consistency constant between `1.95`
(15 deg) and `0.33` (25 deg), against `5.3691` today.  Only the 15 deg rung asks
for something a filter plausibly satisfies -- an estimation error inside about
two applied sigma.  The 25 deg rung would require the latent-acceleration error
to stay below a third of one sigma, which is not a realistic filter property.

So there are exactly two doors, and the table prices both:

1. **Narrow the ladder and declare the coupling.**  Take the candidate to
   15--20 deg and add `||delta a_w|| <= c * sigma_applied` to the operating
   domain with `c` around `1.2`--`1.95`, justified from the tuner law and a
   bound on the EMA transient.  That is a new deployment theorem assumption and
   must be reviewed as one, not slipped in as a proof convenience.
2. **Stop testing sector invariance per operation.**  Charge each correction
   against its own Joseph information decrease and carry a directional block
   margin, so a transient excursion is allowed as long as the Lyapunov level
   decreases.  Under that criterion the 30 deg candidate is not excluded by the
   table above, because the table measures invariance, not dissipation.

Door 2 is the change the route-ceiling certificate already named, and it is the
one that does not add a domain assumption.  Door 1 is cheaper but narrows the
certified set and needs its own physical justification.

The producer reports a **distance and never a verdict**, and marks
`aw_sigma_consistency_declared_in_domain = false`: nothing here declares the
coupling, and nothing here is promoted.

## What would change these numbers

Only the two structural changes matter now; further sharpening of the gain, the
covariance enclosure or `delta` will not move the budget ratio, for the same
reason the route ceiling was independent of `kappa`:

1. pair each accepted correction with its own Joseph information decrease
   `W - W+ = ||H z||^2_{s S^-1}` on the step that injects it;
2. carry a directional/block margin so the attitude channel is not rate-limited
   by the slowest source channel of the word;
3. only then compose the finite-angle sector contraction over the complete
   18/21-state word and set a finite `N_inner` for P5.

Until those land, the honest status is the one the report emits: five usable
stage geometries, two open theorem obligations, and no complete P1--P5 proof.
