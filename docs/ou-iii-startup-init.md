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
2. The gravity gate and the magnetic accumulation frame read the proxy tilt.
3. When tilt has held agreement with gravity, north has been gauged, and the
   tuner is ready, `goLive()` seeds the MEKF with the finished attitude and it
   goes live in one step — linear block on, nominal `Racc`, operating point
   already converged. It never occupies the staged warmup at all.

The inner `SeaStateFusionFilter_OU_III` still defaults to `StagedMekf`. Only
something above it can perform the handoff, so a filter driven directly through
`updateTime()` would park at `TunerReady` forever if it defaulted otherwise.

## Two things the change forced

**The bootstrap observer needs an integral term.** The vertical channel runs at
`two_ki = 0` and accepts the `~2b/two_kp` static tilt error a gyro bias leaves,
because everything downstream of it is high-passed. Nothing high-passes an
attitude *seed*: that error is a standing roll and pitch bias. Measured, a
0.05 deg/s bias settles at 0.71 deg of tilt with no integral term. Seeding from
a zero-integral observer cost 0.17 deg of roll RMS on jonswap H1.5. The
bootstrap therefore uses a second instance of the same class at
`two_ki = 0.02`, and the vertical channel's documented tuning is untouched.

**The magnetometer has to wait for the observer to settle.** The same low corner
that makes the proxy reject waves makes it slow to converge from its
accelerometer seed, and that seed is worst in the big seas where levelling
matters most. Mean tilt error against truth over a 23 s averaging window:

| window start | 7 s | 20 s | 40 s | 60 s | 120 s |
|---|---|---|---|---|---|
| best record  | 0.33 | 0.36 | 0.13 | 0.09 | 0.07 |
| worst record | 2.76 | 2.79 | 0.85 | 0.92 | 0.74 |

`proxy_mag_settle_sec` holds the whole magnetometer path off until then. Under
the old policy this wait was unaffordable because the MEKF was live and needed
the magnetometer; here it has not started and the tuner is converging anyway.

## Results

Eight reference records (jonswap and pmstokes, H 0.27 to 8.5 m), trailing 900 s,
mean across records:

| | Z %Hs | 3D % | roll deg | pitch deg | yaw deg | acc bias m/s^2 |
|---|---|---|---|---|---|---|
| `staged_mekf` | 4.451 | 20.318 | 0.373 | 0.274 | 1.795 | 0.0632 |
| `mahony_proxy` | 4.441 | 20.242 | 0.310 | 0.259 | 1.839 | 0.0517 |

`staged_mekf` reproduces the pre-change numbers exactly, so the flag restores
the old path rather than approximating it.

Roll improves 17%, accelerometer bias 18%, pitch 5%, 3D and Z slightly. All
eight quality gates pass under both.

**Yaw is 2.5% worse**, and it is worse at every settle time tried. The residual
is a property of the learned reference *vector*, not of the one-time gauge:
seeding the handoff yaw variance anywhere from 5 to 90 degrees moves the scored
yaw by under 1e-4 deg, so the filter is not converging to the gauge, it is
converging to the reference. Improving it means re-learning that reference in
the live MEKF's own tilt frame once it has converged — a runtime re-acquisition
rather than an initialization change, and deliberately not attempted here.

## Knobs

| setting | default | what it does |
|---|---|---|
| `startup_init_policy` | `MahonyProxy` | selects the policy |
| `proxy_mag_settle_sec` | 90 | observer settling before mag accumulation |
| `proxy_startup_min_sec` | 8 | earliest handoff |
| `proxy_startup_timeout_sec` | 150 | latest handoff; raised internally so it can never cut the mag acquisition short |
| `proxy_handoff_tilt_sigma_rad` | 0.035 | seeded tilt variance |
| `proxy_handoff_yaw_sigma_rad` | 0.087 | seeded yaw variance once north is gauged |
| `acc_bias_unlock_mag_updates` | 250 | magnetometer updates after going live before the accelerometer-bias gate opens |

The last one is unchanged from the staged policy. Accelerometer bias and tilt
error are only weakly separable in waves, so opening it early lets the bias
absorb the seed's tilt error; sweeping it from 250 to 25000 moves roll by under
0.02 deg while 3D error degrades badly at the long end, so 250 stands.
