# OU-III ALT proof handover after PR #524

This file is the canonical continuation note for the ALT proof after PR #524.
After PR #524 is merged, sync the latest `main`, read `AGENTS.md`,
`docs/ou3-alt-proof-plan.md`, this file, `docs/ou3-alt-contraction-handover.md`,
the ALT section of `docs/ou3-proof-research-state.md`, and
`docs/ou3-brmm-main-handover.md`, then continue ALT in a new PR from `main`.

The original P2/P3/P4/P5 proof remains a separate independent proof track.
Do not delete it, weaken it, rewrite it to depend on ALT, or change frozen
`P3=1e-18`. Shipping implementation and quality gates remain authoritative.

## What PR #524 added

PR #524 substantially extends the finite-runtime ALT product rather than doing a
storage/rho experiment.

The retained product now carries one admitted corrected COMPLETE-BRMM physical
history, one BIAS0/BIAS1/BIAS2 history, one bounded sensor/disturbance history,
one persistent Live/S origin, full joint24 mean-error/true-bias coordinates,
full 21-state covariance, tuner/frontend/scheduler state, and asynchronous
magnetic/startup control state through the represented finite word.

The machine-side coefficient history is no longer detached from the actual
frontend path. The following relations are represented on the same admitted IMU
history:

- float API dt/gyro/accel projection into the shipping binary32 boundary;
- the persistent AccelVibrationGuard recurrence;
- the common private-Mahony vertical observer driven by that guarded sample;
- tracker LPF and stillness state derived from the same Mahony vertical output;
- machine band/statistics/sigma histories for separate/FMA compiler tracks;
- staged tuner state, pending commit, active tau/Sigma_aw/R_S/cadence state;
- prediction coefficient displacement, pending a_w covariance-floor snapshots,
  covariance hygiene and pseudo-S scheduler/service;
- the machine `Racc` vibration/RAO/restore branch derived from the same guarded
  and tuner history rather than held equal to the exact-shadow value;
- propagation of the machine-vs-exact coefficient displacement through the
  accelerometer measurement path;
- a joined wrapper requiring the frontend/stillness source and the Racc/TuneState
  event to be the same already executed admitted event, not two replayed events.

Guard configuration and tracker-LPF cutoff ancestry were also tightened. The
machine guard must be the binary32 image of the carried theorem runtime guard
configuration, and the tracker LPF uses the shipping default constructor cutoff
(`MAX_FREQ_HZ`, currently 6 Hz) rather than a freely supplied theorem value. The
private-Mahony configuration adapter uses the actual `vertical_cfg` fields
`two_kp`, `two_ki`, `gravity`, and `settle_sec` projected to binary32.

The Racc scalar projection bug encountered during continuation was corrected;
the joined machine frontend/Racc product now uses the proper runtime/config
ancestry. Do not reintroduce the earlier erroneous `accel_gate` read from
`vertical_cfg`.

## Important proof semantics retained

These machine histories are coefficient/supply relations against the one actual
shipping state/history. They are not second physical histories and not a second
persistent estimator. Separate/FMA branches may retain distinct compiler
arithmetic histories where the proof explicitly allows that, but they must share
the same source/event ancestry.

Do not replace same-history relations by independent parameter boxes, replay
successful events, restart S at word boundaries, assume covariance consistency,
shrink the physical domain merely to get a certificate, reuse exact-shadow
historical snapshots after machine/compiler histories diverge, or start a
storage/rho search before the finite-master guard closes.

ALT certified scope still excludes dynamic wind heel: `wind_heel_rad_==0` from
construction onward and no `update_wind_heel()` events.

## Current fail-closed status

The new finite maps are genuine progress, but there is still no end-to-end ALT
stability theorem.

Keep:

- `ALT_LIVE_PASS=false`
- `ALT_STARTUP_PASS=false`
- `ALT_END_TO_END_PASS=false`
- `storage_search_allowed=false`

The inherited ALT full-suite currently also exposes blockers from shared source
prerequisites. In particular, the continuous/private-Mahony source invariant is
still fail-closed (`continuous_all_live_PI_invariant_closed` / initial seed angle
not closed), and the inherited Riccati/source path can fail when the declared PE
certificate does not refine the vector certificate. Those are theorem/source
blockers; do not weaken them just to turn CI green.

There is also a pre-existing exact-test sensitivity in an endpoint-partition
assertion around 0.5 caused by outward interval rounding. Treat it separately
from theorem validity rather than silently narrowing the interval arithmetic.

## Next proof work

Continue from the joined finite event product, not from disconnected snapshots.
The next conversation should prioritize the following mathematical obligations:

1. Attach the now-bound guard -> private-Mahony -> tracker-LPF/stillness ->
   band/sigma -> tuner -> Racc chain through the actual startup/goLive history,
   so fresh Live entry carries these machine states instead of seeding them at
   Live.
2. Close the private-Mahony/startup capture and continuous all-Live invariant on
   the admitted corrected COMPLETE-BRMM source without fake numerical experiments
   or a reduced physical operating domain.
3. Finish raw WPE period/log ancestry and source-uniform deployment arithmetic
   supplies: exp/log/sqrt, band/statistics/tuner arithmetic, Q-axis operation
   roundoff, scheduler `nextafter`, a_w-sync clock arithmetic, trig/
   normalization, Eigen LDLT/eigensolver/floor decisions, comparisons and
   nonfinite branches.
4. Qualify every literal IMU/magnetic/hold/reset event of the exact 600-transition
   finite word with the same COMPLETE-BRMM + BIAS + disturbance history,
   including the corrected startup origin and H18/A21 semantics.
5. Only after the finite-master guard accepts that complete source-uniform word,
   search for common/compatible joint24 storage, prove every-prefix chart
   retention, and derive a disturbance-dependent ultimate bound.
6. For an indefinite theorem, additionally resolve finite-width wrapper-clock
   stall, signed magnetic-counter lifetime, and finite-word tiling without
   restarting the one-time Live/S origin or any frontend/covariance/bias state.

## Validation and handoff rule

Before claiming any newly added layer as closed, run the focused finite tests and
native shipping correspondence that exercise it, then inspect the full inherited
ALT suite. Distinguish new regressions from intentionally fail-closed inherited
proof prerequisites. Do not promote component regression success into universal
source admission or storage permission.

At the time this handoff was written, the current PR head already contained the
fix replacing the erroneous `vertical_cfg.accel_gate` reference with
`vertical_cfg.settle_sec`. The GitHub workflows for that head were still being
scheduled/executed, so the next conversation must inspect the merged-main CI
state rather than relying on an earlier failed run from a stale head.
