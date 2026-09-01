# OU-III startup initialization: staged MEKF vs Mahony proxy

> The "22-52 s" and "roughly 20 s" figures below were what the gravity gate
> delivered when it happened to close. It usually did not: measured across the
> eight records, time to a live filter actually ran 22 s on the calm ones and
> 72 to 150 s on the big ones, the worst by timeout. That was a defect in the
> gate rather than in this policy, and `docs/ou-startup-gravity-gate.md` is the
> fix; with it the same policy reaches live in 22 to 33 s on all eight. The
> attitude-quality argument in this document is unchanged, but read its startup
> timings as the intent rather than as measurements.

This document records the comparison that made the Mahony proxy the startup
path for OU-II and OU-III. Both wrappers now have only that path: the
`startup_init_policy` selector, the `StagedMekf` arm and the simulator's
`W3D_STARTUP_INIT` variable have all been removed, so the numbers below are a
record rather than something re-runnable there. `SeaStateFusionFilter_TFG` still
carries both arms for anyone who wants to reproduce the shape of the result.

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
3. When tilt has held agreement with gravity *on the aligned branch*, north
   has been gauged, and the tuner is ready, `goLive()` seeds the MEKF with the
   finished attitude and it goes live in one step — linear block on, nominal `Racc`, operating point
   already converged. It never occupies the staged warmup at all.
4. The magnetic reference is re-learned once the observer has converged, and
   accelerometer-bias learning waits for that; see below.

The gate's residual is `||s_hat x s_meas||`, a sine, and a sine reads the same
at an angle and at its supplement: on its own it scores an attitude and its
180-degree flip alike. So the gate also tests `s_hat . s_meas > 0` and only
counts hold time on the aligned branch. The sign test is a condition, not a
threshold, so it has no knob; both timeouts -- the bootstrap gravity timeout and
`proxy_startup_timeout_sec` -- are held to it too, which delays a forced handoff
rather than letting it seed the filter upside down. The antipodal set is not
attracting for an accel-corrected observer, so that delay is bounded. This is
the aligned-branch condition the semiglobal startup theorem assumes of an
accepted handoff.

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
| `staged_mekf` | ~22 s | 4.301 | 19.872 | 0.372 | 0.275 | 1.796 | 0.0631 |
| `mahony_proxy`, two-stage (default) | 22-52 s | 4.287 | 19.851 | 0.339 | 0.266 | 1.835 | 0.0574 |
| `mahony_proxy`, single-stage | ~105 s | -- | -- | 0.310 | 0.259 | 1.839 | 0.0517 |

Measured after merging main's r_S adaptation-law work. That work moves
displacement and this one moves attitude, and they are independent: main alone
scores Z 4.302, 3D 19.876, roll 0.375, yaw 1.796, bias 0.0636, so the
`staged_mekf` row above reproduces it and the proxy row adds the attitude gain
on top without giving any of the displacement gain back. The single-stage row's
displacement figures are from the pre-merge tree and are not restated; its
attitude figures are the ones it was chosen for.

`staged_mekf` reproduces the pre-change numbers exactly, so the flag restores
the old path rather than approximating it.

The default improves roll 10%, accelerometer bias 10%, pitch 3%, 3D and Z
slightly, at the same time to first heading as before. All eight quality gates
pass in every configuration.

The single-stage row is the accuracy ceiling of this approach and is one
setting away (`mag_refine_enabled = false`, `proxy_mag_settle_sec = 90`): it
roughly doubles the roll and bias gain, and costs about 80 s of extra time
before the device reports a heading. Which of those matters is an application
decision, not a filter one, so both are supported and the fast one is the
default.

Per-record, the default is not uniformly better than the staged path: roll on
pmstokes H4 goes 0.242 -> 0.405 deg and on jonswap H0.27 0.234 -> 0.279, paid
for by 0.564 -> 0.256 on jonswap H1.5 and 0.396 -> 0.265 on pmstokes H0.27.
The single-stage configuration has its own smaller spread. Both improve the
mean; neither dominates record by record.

### Yaw

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
| `proxy_startup_timeout_sec` | 150 | latest handoff; raised internally so it can never cut the mag acquisition short, and still held to the aligned-branch test |
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

## Why the observer gains are not the yaw lever

The obvious next move on the yaw residual is to retune the Mahony observer that
frames the magnetic reference. It was swept, and it does not work. Recording
the result so the sweep is not repeated.

A full grid over `two_kp` in 0.05-0.4 and `two_ki` in 0.005-0.05 moves the
eight-record mean yaw only between 1.825 and 1.856 deg. On the calmest record
yaw is 2.07 deg in *every* cell of the grid, to three decimals.

The reason shows up immediately with the noise generators disabled: yaw RMS
falls from 2.071 to 0.153 deg on jonswap H0.27. About 93% of it is
magnetometer error, not attitude error. Heading is constant in these records,
so the hard-iron offset is not separable from the reference and is simply
absorbed into it; what is left is soft-iron distortion modulated by the wave
motion, which no attitude estimator can undo and no averaging window can
cancel.

Seed sensitivity says the same thing. The same configuration scores 1.83, 3.70
and 3.23 deg of yaw on three noise seeds, while the gains move it by less than
0.01 deg within each. Yaw here is a property of the magnetometer realization.

The sweep did turn up an apparent roll optimum at `two_kp = 0.1`,
`two_ki = 0.005` -- mean roll 0.340 -> 0.298 deg and bias 0.057 -> 0.053 on the
default seed. It was not adopted, because it does not replicate: on a second
seed it wins by a similar margin and on a third it loses by one, which is the
signature of fitting a single realization rather than a real improvement. The
optimum is also sharp -- `two_ki = 0.002` jumps roll to 0.560 deg -- which is
the other sign of the same thing. The gains stay at `two_kp = 0.2`, the
long-standing deployed corner, and `two_ki = 0.02`, the smallest change the
attitude seed actually requires.

Improving yaw further means calibrating the magnetometer, not retuning the
observer: either estimating the hard iron (`mag_estimate_hard_iron`, which
needs heading change during startup to be separable) or a soft-iron model,
neither of which is a startup-initialization question.

## Quality gates

The `kalman_ou_iii-sim` regression sentinels were re-derived once the numbers
above settled. They had accumulated 5-12% of slack against the filter that now
ships -- partly from this work, partly from the r_S adaptation-law work that
landed alongside it -- and slack that wide is somewhere a regression can hide.

Each limit follows the existing convention: the worst value across the eight
scored records, plus about half a percent, rounded up to the next tenth.

| gate | was | now | worst scored | margin |
|---|---|---|---|---|
| Z %Hs, jonswap | 5.4 | 4.8 | 4.69 (H0.27) | 2.3% |
| Z %Hs, pmstokes | 5.3 | 4.7 | 4.66 (H0.27) | 0.8% |
| 3D %, jonswap | 21.4 | 21.1 | 20.95 (H1.5) | 0.7% |
| 3D %, pmstokes | 22.0 | 20.6 | 20.42 (H4.0) | 0.9% |
| yaw deg | 2.2 | 2.2 | 2.16 (pmstokes H1.5) | 1.8% |
| acc Z bias % | 5.9 | 5.6 | 5.48 (pmstokes H8.5) | 2.3% |
| bias 3D % | 109.4 | 95.8 | 95.30 (jonswap H4.0) | 0.5% |

Yaw is unchanged because yaw did not improve; see the section above on why it
cannot be tuned here. `bias_3d_percent` moved the furthest because holding
accelerometer-bias learning until the magnetic reference is refined stops the
bias absorbing the provisional reference's tilt error.

Margins of 0.5-2.3% are small on purpose. The metrics are deterministic to
6e-6 relative across `-march=native`, `x86-64` and `x86-64-v2`, so these sit
about three orders of magnitude above the reproducibility floor.

These are fitted to the deployed default. The `W3D_STARTUP_INIT=staged_mekf`
ablation now exceeds two of them -- 3D on jonswap H1.5 (21.15 against 21.1) and
the bias aggregate (106.2 against 95.8) -- which is the point: those are the
numbers the default was changed to improve. Run the ablation under
`W3D_COLLECT_ALL_GATES=1` to get its metrics instead of an early exit.
