# OU-III startup initialization: staged MEKF vs Mahony proxy

`SeaStateFusion_OU_III::Config::startup_init_policy` selects which estimator
solves the attitude the filter starts from. `MahonyProxy` is the default;
`StagedMekf` restores the previous behaviour and is the matched ablation.
The simulator exposes both as `W3D_STARTUP_INIT=mahony_proxy|staged_mekf`.

## What the staged policy did

The MEKF was fed from the first sample and learned its own tilt while running
in a deliberately degraded configuration: linear block off, accelerometer bias
frozen, `Racc` inflated 2.5x. The wrapper then read that attitude back out for
two jobs:

- the gravity-agreement gate that decides when the magnetometer may start, and
- the tilt frame the magnetic world reference is averaged in.

Both reads happen while the MEKF is least able to level. Its accelerometer
update is fighting the orbital specific force with the linear block — the thing
that would otherwise absorb that acceleration — switched off.

That would be a bounded startup transient if the results were revisable, but
they are not. The world reference and the yaw gauge are locked exactly once, so
whatever tilt error survives the averaging window becomes a standing attitude
and heading bias for the rest of the run. It is scored: the executables report
the trailing 900 s of a 1200 s record, long after any transient has died.

## What the proxy policy does

The filter already runs a private Mahony observer (`VerticalAccelComplementary`)
that is a pure function of the raw gyro and accelerometer. Its correction corner
sits an order of magnitude below the wave band, so it rejects orbital
acceleration by construction instead of chasing it.

Under `MahonyProxy`:

1. The measurement-only front end — proxy, frequency tracker, wave-period
   estimator, sigma band, auto-tuner — runs from the first sample via
   `updateFrontEnd()`, with the MEKF untouched. None of those consumers reads a
   filter state, so this changes none of their outputs; `startup_init-test`
   pins that bit-for-bit.
2. The gravity gate and the magnetic accumulation frame read the proxy tilt --
   for both acquisition stages, not just the first.
3. When tilt has held agreement with gravity, north has been gauged, and the
   tuner is ready, `goLive()` seeds the MEKF with the finished attitude and it
   goes live in one step — linear block on, nominal `Racc`, operating point
   already converged. It never occupies the staged warmup at all.
4. The magnetic reference is re-learned once the observer has converged, and
   accelerometer-bias learning waits for that; see below.

The inner `SeaStateFusionFilter_OU_III` still defaults to `StagedMekf`. Only
something above it can perform the handoff, so a filter driven directly through
`updateTime()` would park at `TunerReady` forever if it defaulted otherwise.

## Two things the change forced

**The observer needs an integral term.** The vertical channel ran at
`two_ki = 0` for years and accepted the `~2b/two_kp` static tilt error a gyro
bias leaves, because everything downstream of it is high-passed. Nothing
high-passes an attitude *seed*: that error is a standing roll and pitch bias.
Measured, a 0.05 deg/s bias settles at 0.71 deg of tilt with no integral term,
and seeding from a zero-integral observer cost 0.17 deg of roll RMS on jonswap
H1.5.

This was briefly solved with a second observer instance so the vertical
channel's tuning could stay untouched. Turning the integral term on for both
and keeping one instance turns out to be free: over the eight records the
wave-period, sigma and r_S channels are unmoved and every scored metric is
identical to three decimals. Estimating the bias is also the better answer for
the vertical channel on its own terms -- the static tilt it tolerated was
leaking gravity into the levelled acceleration and being high-passed away
afterwards. So there is one observer, at `two_kp = 0.2`, `two_ki = 0.02`.

**Magnetic acquisition has to run twice.** The same low corner that makes the
observer reject waves makes it slow to converge from its accelerometer seed,
worst in the big seas where levelling matters most. Mean tilt error against
truth over a 23 s averaging window:

| window start | 7 s | 20 s | 40 s | 60 s | 120 s |
|---|---|---|---|---|---|
| best record  | 0.33 | 0.36 | 0.13 | 0.09 | 0.07 |
| worst record | 2.76 | 2.79 | 0.85 | 0.92 | 0.74 |

Waiting for the good frame before locking anything puts first heading around
105 s, which is not a usable device. So acquisition runs in two stages: a
provisional lock as soon as the gravity gate allows, giving heading and a live
filter in 22-52 s, and a refinement at `mag_refine_start_sec` that re-learns
the reference and re-gauges heading once the observer has converged.

Two things about the refinement are not optional:

- **It must be framed on the observer, not the MEKF.** By then the MEKF looks
  like the better frame, but it has been steering to the provisional reference
  the refinement exists to replace, so its tilt carries that reference's error
  and averaging the field in it re-derives the error. Framed that way the
  refinement is self-confirming: reference and yaw come back within 1e-3 deg of
  the provisional ones and the standing roll bias is untouched.
- **Accelerometer bias must be held until it lands** (`setAccBiasHold`). Bias
  and tilt error are barely separable in waves and the bias state has a 5000 s
  correlation time, so a value fitted to the provisional reference outlives the
  record. Without the hold the early-lock path fails the bias gate outright
  (187% against a 109% limit).

Refining *later* is worse, not better -- 90 s beats 120 and 150 -- because
damage accumulates for as long as the filter steers to the provisional
reference.

## Results

Eight reference records (jonswap and pmstokes, H 0.27 to 8.5 m), trailing 900 s,
mean across records:

| | first heading | Z %Hs | 3D % | roll deg | pitch deg | yaw deg | acc bias m/s^2 |
|---|---|---|---|---|---|---|---|
| `staged_mekf` | ~22 s | 4.451 | 20.318 | 0.373 | 0.274 | 1.795 | 0.0632 |
| `mahony_proxy`, two-stage (default) | 22-52 s | 4.435 | 20.286 | 0.340 | 0.266 | 1.834 | 0.0574 |
| `mahony_proxy`, single-stage | ~105 s | 4.441 | 20.242 | 0.310 | 0.259 | 1.839 | 0.0517 |

`staged_mekf` reproduces the pre-change numbers exactly, so the flag restores
the old path rather than approximating it.

The default improves roll 9%, accelerometer bias 9%, pitch 3%, 3D and Z
slightly, at the same time to first heading as before. All eight quality gates
pass in every configuration.

The single-stage row is the accuracy ceiling of this approach and is one
setting away (`mag_refine_enabled = false`, `proxy_mag_settle_sec = 90`): it
roughly doubles the roll and bias gain, and costs about 80 s of extra time
before the device reports a heading. Which of those matters is an application
decision, not a filter one, so both are supported and the fast one is the
default.

Per-record, the default is not uniformly better than the staged path: roll on
pmstokes H4 goes 0.240 -> 0.405 deg and on jonswap H0.27 0.230 -> 0.281, paid
for by 0.560 -> 0.256 on jonswap H1.5 and 0.408 -> 0.267 on pmstokes H0.27.
The single-stage configuration has its own smaller spread. Both improve the
mean; neither dominates record by record.

**Yaw is about 2% worse** in every proxy configuration. The residual is a
property of the learned reference *vector*, not of the one-time gauge: seeding
the handoff yaw variance anywhere from 5 to 90 degrees moves the scored yaw by
under 1e-4 deg, so the filter is not converging to the gauge, it is converging
to the reference. The refinement pass recovers part of it (2.008 -> 1.952 deg
on pmstokes H8.5) but not all, and re-framing the refinement on the MEKF to try
to do better is exactly the self-confirming loop described above. Closing the
gap needs a better exogenous tilt frame, not a later or longer average.

## Knobs

| setting | default | what it does |
|---|---|---|
| `startup_init_policy` | `MahonyProxy` | selects the policy |
| `proxy_mag_settle_sec` | 0 | observer settling before the *provisional* lock |
| `mag_refine_enabled` | true | run the second-stage acquisition |
| `mag_refine_start_sec` | 90 | when the refinement begins |
| `mag_refine_window_sec` | 30 | its averaging window |
| `proxy_startup_min_sec` | 8 | earliest handoff |
| `proxy_startup_timeout_sec` | 150 | latest handoff; raised internally so it can never cut the mag acquisition short |
| `proxy_handoff_tilt_sigma_rad` | 0.035 | seeded tilt variance |
| `proxy_handoff_yaw_sigma_rad` | 0.087 | seeded yaw variance once north is gauged |
| `acc_bias_unlock_mag_updates` | 250 | magnetometer updates after going live before the accelerometer-bias gate opens |

The last one is unchanged from the staged policy, and under this policy it is
inert: 0, 25, 50, 100 and 250 agree to the fourth decimal on every scored
metric, because by the time the filter goes live the attitude is already good
and there is nothing for the gate to protect against. It still matters at the
long end -- 2500 and above degrade 3D error badly -- so it stays at 250 rather
than being removed. The work it used to do is now done by `setAccBiasHold`,
which is keyed to the refinement rather than to a magnetometer-update count.
