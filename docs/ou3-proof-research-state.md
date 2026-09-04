# OU-III proof research state

This file is the current proof ledger. It records the active theorem, completed
certificate facts, open obligations, and the next falsifiable step. Historical
routes belong in Git history.

## Active theorem

The target is **uniform regional practical stability with finite capture for
perturbations of physically admissible multimodal directional sea states**.

The fixed-dimensional sea class has `M_max = 3` and

`E(omega,theta) = sum_{r=1}^3 S_J(omega; H_r,T_p,r,gamma_r) D_r(theta; beta_r,s_r)`.

Inactive partitions use `H_r = 0`; PM is the `gamma_r = 1` boundary; modal
energies obey `H_s^2 = sum_r H_r^2`.  One-system seas, wind sea + swell,
crossing swell, and three-component seas therefore share one theorem class.

The physical peak periods are not identified with the deployed tuner period.
`WavePeriodEstimator` produces its zero-crossing/moment period after its own
leaky front end and finite-memory statistics.  In multimodal seas, moments
combine before the period ratio is formed.

The source chain remains

`L_actual_sea subset Lhat_SEA3 subset L_current_source`,

where the right-hand language is the existing implementation-correlated P2
language. SEA3 refines that language; it does not redefine the frozen P2/P3
interface.

Because active partitions need not be commensurate, the nonlinear stability
object is a finite **sea window**, not a common wave cycle.  P4 ultimately has
to certify both

`V_{k+N_W} <= rho V_k + gamma_s D_s,k + gamma_n D_n,k`, `rho < 1`,

and the intermediate-sample bound

`V_{k+l} <= kappa_V V_k + kappa_s D_s,k + kappa_n D_n,k`, `0 <= l < N_W`.

## Current proof architecture

1. **SEA0 — directional sea admissibility.** Three-partition sea/rate domain,
   robust vessel/IMU response family, finite-window realization enclosure, and
   exact estimator/tuner source dynamics.
2. **P1 — operating branch.** Startup/Normal-Live/hybrid assumptions, including
   no rejected accelerometer branch in declared Normal Live.
3. **P2 — SEA3 reachable language.** Preserve sea phase/frequency continuity,
   directional/cross-axis correlation, modal-energy coupling, prior/period
   source mode, tuner memory, stage/commit history, and pseudo-scheduler phase.
4. **P3 — finite-window linear certificate.** H=18 and A=21
   information/covariance bounds with the canonical usefulness gate exactly
   `1e-18`.
5. **P4 — nonlinear lifted dissipation.** Endpoint contraction, prefix gain,
   source/chart containment, and an invariant inner funnel.
6. **P5 — finite capture.** Entry from the declared 45-degree entrance into the
   P4 funnel and composition of consecutive SEA3 windows.

There is one general proof workflow: `.github/workflows/ou3-proof.yml`.
SEA3 response/inclusion is bound into the source-domain test already executed by
that workflow; the temporary `ou3-sea3-directional.yml` workflow is retired and
the cleanup test forbids its return.

## Shipping theorem envelope

Current source/parity values used by the proof are:

- `imu_dt = 0.005 s`;
- `0.03 <= f_tune <= 1.2 Hz`;
- `0.02 <= tau <= 12 s`;
- `sigma_aw <= 4 m/s^2`;
- `0.15 <= R_S <= 100`;
- `0.005 <= T_S <= 0.150 s`;
- dynamic EMA horizon `<= 35 s`;
- every valid Normal-Live IMU sample executes the accelerometer update;
- accelerometer rejection is outside the declared Normal-Live theorem branch;
- lever arm remains disabled;
- the accelerometer vibration guard is restricted to its dormant/transparent
  proof branch.

The pseudo-update scheduler uses the progress-preserving retarget helper.  The
retained scheduler certificate gives a maximum 30-sample recurrence at the
150 ms ceiling.

## Completed SEA0/source-bridge subcertificates

### Spectral-moment bridge

`tools/ou3_sea3_spectral_moment_bridge.py` is replay free and non-promoting.
Established analytical facts include:

- PM (`gamma=1`) surface-elevation `Tz/Tp = ((5/4)*pi)^(-1/4)`;
- for normalized uncorrelated partitions,
  `1/Tz_mix^2 = sum_r w_r/Tz_r^2`, `w_r=H_r^2/sum_j H_j^2`;
- ideal unbanded JONSWAP/PM *surface-acceleration* variance diverges because the
  `omega^4` weighting turns the `omega^-5` tail into `omega^-1`.

The padded numerical `1 <= gamma <= 7` screen remains a feasibility screen, not
an interval-integration theorem promotion.

### WavePeriodEstimator front end and startup causality

`tools/ou3_sea3_wave_period_frontend.py` remains replay free and non-promoting.
It source-certifies the fixed 0.2 Hz prior, the prior-to-first-finite-estimator
handoff, the one-sample tuner/period-estimator ordering edge, and the fact that
filter `TunerReady` is not WavePeriodEstimator readiness.  Its single-frequency
front-end certificate keeps discrete period warping below about 59 ppm on the
committed 5 ms / 0.03--1.2 Hz channel.

### Robust directional RAO envelope family

`tools/ou3_sea3_directional_response_domain.json` and
`tools/ou3_sea3_directional_p2_ha_feasibility.py` now certify a **continuous RAO
range**, not one hull and not one gain value.

For the complex three-axis CoG translational projection `h(f,theta)`, the family
is

`||h(f,theta)||_2 <= G min(1,(f_c/f)^p)` for `f>0`, with

- `0 <= G <= 4`;
- `0.03 <= f_c <= 1.2 Hz`;
- `p >= 2`;
- arbitrary complex phase;
- arbitrary frequency and heading dependence below the envelope;
- arbitrary cross-axis coupling consistent with the PSD outer product;
- an arbitrary six-DOF parent RAO, with rotational response handled by the
  separate P1 body-rate/attitude source bounds on the zero-lever-arm branch.

`G=4` is therefore **not one RAO**.  It is the upper face of a continuum of
three-axis response norms.  Likewise `f_c=1.2,p=2` is the monotone worst
*envelope corner*, not a representative vessel.

For every member of the family,

`tr M_disp <= G^2 H_s^2 / 16`,

`tr M_vel <= (2*pi*f_c)^2 G^2 H_s^2 / 16`,

`tr M_acc <= (2*pi*f_c)^4 G^2 H_s^2 / 16`.

Because these bounds increase with `G` and `f_c` and `p=2` is worst, verifying
`G=4, f_c=1.2, p=2` proves the whole parameter box.  No RAO or frequency grid is
used in that theorem.

The `p>=2` response roll-off also fixes the old high-frequency pathology:
above `f_c`, the squared `f^-2` response cancels the `f^4` acceleration
weighting.  The response-weighted acceleration moment is therefore finite even
without a hard 6 Hz truncation.  Compared with the old flat `G=4` outer corner
carried all the way to 6 Hz, the worst acceleration-moment coefficient improves
by approximately `(6/1.2)^4 = 625`.

This is a broad **declared deployment response envelope**, not a claim that every
conceivable vessel in the world satisfies it.  SEA0 still has to establish that
the physical vessel population admitted by the theorem lies inside it.  Jensen,
Mansour & Olsen, Ocean Engineering 31 (2004) 61--85,
doi:10.1016/S0029-8018(03)00108-2, is retained as engineering motivation for a
principal-dimension/operational-profile response family rather than a single
identified hull; it is not used as a universal envelope proof.

### SEA3-to-current-P2 right inclusion

The same producer mechanically certifies

`Lhat_SEA3 subset L_current_source`

for every Normal-Live SEA3 execution generated by any RAO inside the complete
parameter box.  The right inclusion itself is independent of which RAO member
is selected: shipping clamps `f_tune`, `tau`, `sigma_aw`, and `R_S` before the
same EMA/stage/one-sample-pending-apply semantics already over-approximated by
P2.

The result is deliberately **non-pruning** and still contains all 800 current P2
physical tuner cells.  That is sound inclusion, not yet a source-reachability
reduction.

The canonical source-domain test calls `sea3.build_inclusion()` directly, so the
right inclusion is now exercised by the existing `ou3-proof` workflow instead
of a parallel SEA3 workflow.

## Current certificate status

- **SEA0:** partial. Surface spectral bridge, estimator front-end/startup
  causality, robust continuum RAO-envelope moment theorem, and the right source
  inclusion now exist.  A hard finite-window physical sea realization/IQC
  enclosure and proof of physical-vessel coverage of the RAO envelope remain
  open.
- **P1:** existing Normal-Live assumptions retained.  No Gaussian pointwise
  acceleration bound is inferred from `H_s` alone.
- **P2:** `Lhat_SEA3 subset L_current_source` is mechanically closed for the
  whole RAO parameter box, but no SEA3-specific P2 cells are yet removed.
- **P3:** H=18/A=21 remains the unique canonical full-P2 route.  A canonical
  full-P2 PASS transfers to the entire SEA3 RAO family; a full-P2 FAIL is only
  inconclusive for SEA3 until response/finite-memory source pruning is built.
- **P4:** no SEA3 endpoint/prefix nonlinear certificate is promoted yet.
- **P5:** no SEA3 finite-capture certificate is promoted yet.

No SEA3 bridge changes the filter, shrinks the current P1 operating domain,
changes the `1e-18` gate, or creates alternate P3/P4 promotion authority.

## Current blocker

The immediate soundness blocker is the **left physical inclusion**:

`L_actual_sea subset Lhat_SEA3`.

SEA0 must supply a replay-free hard finite-window oscillator/IQC or equivalent
deterministic realization enclosure for the three-partition directional sea,
plus a physical-vessel qualification showing that the admitted ship/boat
response operators lie inside the declared RAO parameter box and that their
executions satisfy the existing P1 Normal-Live source bounds.  A random Gaussian
spectrum alone cannot establish an infinite-time hard pointwise bound.

The canonical full-P2 H/A calculation decides the next branch:

- if full-P2 P3 passes, its margin is inherited uniformly over the whole SEA3
  RAO box and source pruning is unnecessary for P3;
- if it fails, propagate the response-weighted WavePeriodEstimator/log-period,
  adaptive wave-band variance, tuner EMA, stage/commit state, and pseudo-clock to
  construct an actually narrower SEA3 history language, then rerun the same
  H/A theorem interface.

## Next proof increment

1. Finish the canonical full-P2 H/A run and record its actual H/A margins.
2. Add the replay-free finite-window sea realization/IQC contract and physical
   vessel qualification for the declared RAO envelope, closing the left
   inclusion.
3. If canonical P3 is not already positive, use the new response envelope in
   the finite-memory wave-period/variance/tuner reachability calculation to
   prune P2 source histories and recompute H/A without changing the canonical
   gate.
4. Once SEA3 P3 is unambiguous, continue to P4 endpoint + prefix dissipation and
   P5 finite capture.

## Retired / forbidden shortcuts

Do not reintroduce:

- one nominal hull RAO or a finite sampled RAO catalogue as a universal proof;
- flat high-frequency RAO gain carried to the 6 Hz band endpoint;
- independent Cartesian `tau/sigma/R_S` extrema;
- replay-selected sea/source words as a certificate domain;
- equating `T_p` with deployed `T_z` or averaging modal periods;
- unbanded JONSWAP surface-acceleration variance;
- treating `TunerReady` as WavePeriodEstimator readiness;
- common-period/Floquet reasoning for incommensurate multimodal seas;
- gate tuning, operating-domain shrink, or deployed-filter retuning to obtain a
  pass;
- parallel fallback OU-III proof workflows.
