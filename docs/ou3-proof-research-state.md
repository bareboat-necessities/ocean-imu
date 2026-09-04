# OU-III proof research state

This file is the current proof ledger. It records the active theorem, completed
certificate facts, open obligations, and the next falsifiable step. Historical
routes belong in Git history.

## Active theorem

The target is **uniform regional practical stability with finite capture for
perturbations of physically admissible multimodal directional JONSWAP sea
states**.

The fixed-dimensional sea class has `M_max = 3` and

`E(omega,theta) = sum_{r=1}^3 S_J(omega; H_r,T_p,r,gamma_r) D_r(theta; beta_r,s_r)`.

Inactive partitions use `H_r = 0`; `1 <= gamma_r <= 7`; PM is exactly the
`gamma_r = 1` boundary of this same JONSWAP family; modal energies obey
`H_s^2 = sum_r H_r^2`. One-system seas, wind sea + swell, crossing swell, and
three-component seas therefore share one theorem class.

The physical peak periods are not identified with the deployed tuner period.
`WavePeriodEstimator` produces a response-weighted zero-crossing/moment period
after its own leaky front end and finite-memory statistics. In multimodal seas,
moments combine before the period ratio is formed.

The source chain remains

`L_actual_sea subset Lhat_SEA3 subset L_current_source`,

where the right-hand language is the existing implementation-correlated P2
language. SEA3 refines that language; it does not redefine the frozen P2/P3
interface.

Because active partitions need not be commensurate, the nonlinear stability
object is a finite **sea window**, not a common wave cycle. P4 ultimately has to
certify both

`V_{k+N_W} <= rho V_k + gamma_s D_s,k + gamma_n D_n,k`, `rho < 1`,

and

`V_{k+l} <= kappa_V V_k + kappa_s D_s,k + kappa_n D_n,k`, `0 <= l < N_W`.

## Current proof architecture

1. **SEA0 — directional sea admissibility.** Three-partition JONSWAP sea/rate
   domain, robust vessel/IMU response family, finite-window realization
   enclosure, and exact estimator/tuner source dynamics.
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
PR #484 was closed unmerged as superseded by #482 because it used a different
provisional RAO family. Reusable response-independent mathematics from #484 is
being retained here; its conflicting RAO constants are not.

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

The pseudo-update scheduler uses the progress-preserving retarget helper. The
retained scheduler certificate gives a maximum 30-sample recurrence at the
150 ms ceiling.

## Completed SEA0/source-bridge subcertificates

### Spectral-moment bridge

`tools/ou3_sea3_spectral_moment_bridge.py` is replay free and non-promoting.
Established analytical facts include:

- PM is the `gamma=1` JONSWAP boundary and has
  `Tz/Tp = ((5/4)*pi)^(-1/4)`;
- for normalized uncorrelated partitions,
  `1/Tz_mix^2 = sum_r w_r/Tz_r^2`, `w_r=H_r^2/sum_j H_j^2`;
- ideal unbanded JONSWAP surface-acceleration variance diverges because
  `omega^4 * omega^-5 = omega^-1` in the tail.

The padded numerical `1 <= gamma <= 7` shape screen remains a feasibility
screen, not an interval-integration theorem promotion.

### Physical sea height/period admissibility

`tools/ou3_sea3_physical_admissibility.py` is replay free and non-promoting. It
removes the unphysical rectangular shortcut in which a partition can
simultaneously take an arbitrary maximum height and minimum peak period.

It uses the DNVGL-RP-C205 peak-steepness recommendation

`S_p = 2*pi*H/(g*T_p^2)`

with `S_p<=1/15` for `T_p<=8 s`, `S_p<=1/25` for `T_p>=15 s`, and linear
interpolation between 8 and 15 seconds. In SEA3 that single-sea rule is applied
separately to each active `(H_r,T_p,r)` as an explicit conservative theorem
design choice, while retaining

`H_s^2 = sum_r H_r^2`, `H_s<=8.5 m`.

The partitionwise multimodal extension is ours; it is not attributed to DNV as
a published multimodal theorem. This subcertificate does not close the finite
window realization or left language inclusion.

### WavePeriodEstimator front end and startup causality

`tools/ou3_sea3_wave_period_frontend.py` remains replay free and non-promoting.
It source-certifies the fixed 0.2 Hz prior, the same estimator's `3/lambda`
moment-start, `4/lambda` plus one-period startup-usable gate, strict
`6/lambda` readiness floor, and the one-sample tuner/period-estimator ordering
edge. OU-II, OU-III, and TFG now require that startup-usable measured period
before Live entry while keeping `WavePeriodEstimator::isReady()` as the stricter
diagnostic state. Its single-frequency discrete front-end certificate keeps
period warping below about 59 ppm on the committed 5 ms / 0.03--1.2 Hz channel.

### Exact arbitrary-spectrum leak identity

`tools/ou3_sea3_wave_period_spectral_identity.py` retains the useful
response-independent analytical result from superseded #484.

For the continuous-time steady-state counterpart of the shipping two-high-pass,
two-leaky-integration estimator,

`H_v(s)=s^2/(s+lambda)^3`,

`H_eta(s)=s^2/(s+lambda)^4`.

Let

`W(omega)=omega^4/(omega^2+lambda^2)^4`.

For **any** nonnegative input specific-force spectrum with finite weighted
moments,

`sigma_v^2/sigma_eta^2 - lambda^2 = int omega^2 W S_a / int W S_a`.

Thus the leak subtraction is exactly a weighted mean-square-frequency identity;
it is not a narrow-band or single-sinusoid approximation. This does not replace
the still-open exact discrete finite-EWMA/log-period transient enclosure and is
not itself a P2-pruning certificate.

### Robust directional RAO envelope family

`tools/ou3_sea3_directional_response_domain.json` and
`tools/ou3_sea3_directional_p2_ha_feasibility.py` certify a **continuous RAO
range**, not one hull and not one gain value.

For the complex three-axis CoG translational projection `h(f,theta)`,

`||h(f,theta)||_2 <= G min(1,(f_c/f)^p)`, `f>0`, with

- `0 <= G <= 4`;
- `0.03 <= f_c <= 1.2 Hz`;
- `p >= 2`;
- arbitrary complex phase;
- arbitrary frequency and heading dependence below the envelope;
- arbitrary cross-axis coupling consistent with the PSD outer product;
- an arbitrary six-DOF parent response, with rotational response handled by the
  separate P1 body-rate/attitude source bounds on the zero-lever-arm branch.

`G=4` is the upper face of a continuum, not one vessel RAO. Likewise
`f_c=1.2,p=2` is the monotone worst envelope corner, not a representative hull.

For every member,

`tr M_disp <= G^2 H_s^2 / 16`,

`tr M_vel <= (2*pi*f_c)^2 G^2 H_s^2 / 16`,

`tr M_acc <= (2*pi*f_c)^4 G^2 H_s^2 / 16`.

The `p>=2` roll-off makes the unbanded response-weighted acceleration moment
finite and removes the old flat-response dependence on the 6 Hz sigma-band
endpoint. Relative to that discarded flat outer corner, the worst
acceleration-moment coefficient improves by about `(6/1.2)^4 = 625`.

This is a broad declared deployment response envelope. SEA0 still has to show
which physical vessel population lies inside it; the envelope itself is not a
claim that every conceivable vessel has those RAOs.

### Coupled JONSWAP-sea / RAO / P1 compatibility

`tools/ou3_sea3_p1_compatibility.py` is replay free and fail closed. It proves
that the **independent Cartesian product**

`all physically admissible JONSWAP sea parameters x all response functions in the RAO envelope`

cannot itself imply the existing hard P1 condition

`||a_non-grav(t)|| <= 4 m/s^2`.

The analytical witness is one JONSWAP partition at `gamma=1`, i.e. the PM
boundary of the declared `gamma in [1,7]` family, with `T_p=8 s` on the DNV
peak-steepness boundary and the admitted envelope corner `G=4,f_c=1.2,p=2`.
PM is used only because its normalization is closed form; this is **not** a
PM-only theorem and no claim is made that `gamma=1` is the worst JONSWAP
member.

Using only positive spectral mass on `x=f/f_p in [1,3]`, the certificate obtains
an outward-rounded lower bound on acceleration mean square that is strictly
larger than `4^2`. Therefore the complete independent outer product cannot be
promoted into P1.

The correct physical theorem domain is consequently a **coupled SEA3 set**:
sea and RAO parameters are not independently maximized, and a physical pair is
admitted to P1 only after its finite-window response realization satisfies the
existing hard Normal-Live source bounds. An RMS or PSD bound alone is not a
pointwise deterministic P1 certificate.

This result closes the compatibility *logic* but deliberately leaves
`finite_window_realization_certificate_closed=false` and
`L_actual_sea_subset_Lhat_SEA3_closed=false`.

### SEA3-to-current-P2 right inclusion

The directional bridge mechanically certifies

`Lhat_SEA3 subset L_current_source`

for every Normal-Live SEA3 execution admitted after the P1/SEA3 physical-domain
conditions. Shipping clamps `f_tune`, `tau`, `sigma_aw`, and `R_S` before the
same EMA/stage/one-sample-pending-apply semantics already over-approximated by
P2.

The result is deliberately **non-pruning** and still contains all 800 current P2
physical tuner cells. That is sound right inclusion, not yet a source
reachability reduction.

## Current certificate status

- **SEA0:** partial. Surface spectral bridge, DNV-coupled sea-height/period
  admissibility, estimator front-end/startup causality, exact steady spectral
  leak identity, robust continuum RAO moment theorem, Cartesian-product
  incompatibility witness, and right source inclusion exist. The hard
  finite-window realization/IQC enclosure and physical-vessel qualification
  needed for the left inclusion remain open.
- **P1:** existing Normal-Live assumptions retained. The independent JONSWAP
  sea × RAO outer product is now explicitly rejected as a sufficient P1
  domain; admission must be checked on the coupled post-response realization.
- **P2:** `Lhat_SEA3 subset L_current_source` is mechanically closed, but no
  SEA3-specific P2 cells are yet removed.
- **P3:** H=18/A=21 remains the unique canonical full-P2 route. A canonical
  full-P2 PASS transfers to the admitted SEA3 subset; a full-P2 FAIL is only
  inconclusive for SEA3 until response/finite-memory pruning is built. A clean
  current-head CI verdict is still required before recording numerical H/A
  margins.
- **P4:** no SEA3 endpoint/prefix nonlinear certificate is promoted yet.
- **P5:** no SEA3 finite-capture certificate is promoted yet.

No SEA3 bridge changes the filter, shrinks the current P1 numerical bounds,
changes the `1e-18` P3 gate, or creates alternate P3/P4 promotion authority.

## Current blocker

The immediate physical soundness blocker is the **left inclusion on the coupled
domain**:

`L_actual_sea subset Lhat_SEA3`.

SEA0 must supply a replay-free hard finite-window oscillator/IQC or equivalent
realization enclosure for the three-partition directional JONSWAP sea and prove
that the post-RAO response satisfies the existing P1 Normal-Live source bounds.
A Gaussian spectrum alone cannot establish an infinite-time deterministic
pointwise bound.

Independently, the canonical full-P2 H/A calculation decides whether source
pruning is required for P3:

- if full-P2 P3 passes, its margin transfers to the admitted SEA3 subset and
  P2 pruning is unnecessary for P3;
- if it fails, propagate the response-weighted WavePeriodEstimator/log-period,
  adaptive wave-band variance, tuner EMA, stage/commit state, and pseudo-clock
  to construct an actually narrower SEA3 history language, then rerun the same
  H/A theorem interface.

## Next proof increment

1. Complete a clean current-head canonical full-P2 H/A run and record actual H/A
   margins and fail reasons.
2. Declare/certify the compact physical period/rate domain needed by the
   finite-window SEA3 realization.
3. Build the replay-free finite-window oscillator/IQC response certificate on
   the **coupled** JONSWAP-sea/RAO set and close the left inclusion.
4. If canonical P3 is not positive, propagate the exact discrete finite-memory
   period/variance/tuner state to prune P2 histories and recompute H/A without
   changing the canonical gate.
5. Once SEA3 P3 is unambiguous, continue to P4 endpoint+prefix dissipation and
   P5 finite capture.

## Retired / forbidden shortcuts

Do not reintroduce:

- one nominal hull RAO or a finite sampled RAO catalogue as a universal proof;
- a second OU-III/SEA3 proof PR or parallel proof workflow carrying a competing
  RAO theorem domain;
- flat high-frequency RAO gain carried to the 6 Hz band endpoint;
- treating PM as a separate theorem sea instead of the `gamma=1` JONSWAP
  boundary;
- independent maximum `H_r` and minimum `T_p,r` rectangular corners;
- three independently maximized partition heights while also fixing total Hs;
- an independent Cartesian product of sea extremes and RAO-envelope extremes
  as if it automatically satisfied P1;
- independent Cartesian `tau/sigma/R_S` extrema;
- replay-selected sea/source words as a certificate domain;
- equating physical `T_p` with deployed `T_z` or averaging modal periods;
- unbanded JONSWAP surface-acceleration variance;
- treating `TunerReady` as WavePeriodEstimator readiness;
- common-period/Floquet reasoning for incommensurate multimodal seas;
- gate tuning, operating-domain shrink, or deployed-filter retuning to obtain a
  pass;
- parallel fallback OU-III proof workflows.
