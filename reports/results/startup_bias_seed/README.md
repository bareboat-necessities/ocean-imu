# Startup vertical accelerometer-bias seed from a PII observer (OU-II, OU-III)

**Status: not adopted.** The prototype and the A/B tool that produced
`startup_bias_seed_runs.csv` are in commit `110fa93`; the next commit removes
them, and the shipping filters are unchanged.

## Motivation

In the 28 ft sailboat RAO replays, OU-II and OU-III go Live at about 32 s with a
zero accelerometer-bias estimate. Bias learning is then held until the magnetic
reference refinement completes, at about 120 s (`mag_refine_start_sec` 90 s plus
a 30 s window, then `setAccBiasHold(false)`). During that hold a few mg of
vertical bias drives the heave estimate to a one-sided offset of up to about 3 m
(see `../startup_transient/`). The heave state at go-live is already within
0.1–0.2 m of the truth, so seeding position, velocity, or acceleration would not
help. The quantity that needs a seed is the vertical bias.

## Prototype (commit `110fa93`, `SF_STARTUP_BIAS_SEED=1`)

- While the filter is not yet Live, a fixed-pole PII observer (triple pole at
  r = 0.3 rad/s, no adaptation) runs on the startup Mahony proxy's vertical
  acceleration.
- At go-live, the world-down component of the MEKF bias mean is set to the
  observer's implicit bias `−ks·S`, averaged over the last 8 s before go-live.
  The component is mapped into the MEKF's stored body frame and taken net of the
  temperature term.
- The covariance, the magnetic hold, and every other state are left exactly as
  the deployed hand-over leaves them.

A first version used the shipping adaptive PII. That schedule sets r ≈ 0.09 rad/s
for heave, so the observer's integral settles in about 75 s, which is longer
than the 32 s available before go-live. The seed was better than zero in only
67 % of runs, and the median startup heave improvement was 5–16 %. The fixed-pole
version is the one refinement allowed by the research protocol.

## Results (800 full 20 min replays)

The replays cover OU-II and OU-III, seed off and on, the engine at 2400 rpm and
off, the default seed plus four `W3D_INIT_SEED` initialisation seeds, and all 20
records.

| | OU-II off | OU-II on | OU-III off | OU-III on |
|---|---|---|---|---|
| Startup heave RMS 35–180 s, seed/deployed (median) | 0.77 | 0.78 | 0.74 | 0.76 |
| Peak heave error 0–300 s, seed/deployed (median) | 0.76 | 0.79 | 0.74 | 0.79 |
| Startup yaw RMS, seed/deployed (median) | 1.00 | 1.00 | 1.00 | 1.00 |

- **Share of runs with better startup heave:** 80 %. The median ratio is
  0.62 on Fenton, 0.70 on JONSWAP, and 0.94–1.00 on cnoidal.
- **Scored trailing 900 s window:** pooled x/y/z, roll, pitch, and yaw RMS change
  by at most 0.2 %. Individual runs move by up to ±30–100 % in pitch and yaw,
  in both directions, and never in heave.
- **Seed accuracy:** median error +0.034 m/s² against a median |true bias| of
  0.049 m/s². The error is positive in all 400 seeded runs. The seed beats the
  deployed zero seed in 81 % of runs but is within half the bias in only 21 %.

### Scored regression gates (`W3D_COLLECT_ALL_GATES=1`)

| | seed off (default) | seed on |
|---|---|---|
| OU-II  | 8/8 PASS | 8/8 PASS |
| OU-III | 8/8 PASS | **5/8**: JONSWAP 4 m yaw 0.967° > 0.819°; roll on PM-Stokes 8.5 m; 3D gyro-bias gate on JONSWAP 4 m, JONSWAP 8.5 m and PM-Stokes 8.5 m |

Seed off is logically identical to the deployed filters. Builds with
`-ffp-contract=off` give byte-identical output. With the default `-O3
-march=native` build, the differences are at most 0.1 mm in heave and 1e-5° in
yaw, from FMA contraction.

## Why the seed is inaccurate

The error is positive in every run. It does not come from the reference frame:
on the probed run (JONSWAP 4 m, seed 23), the true world-down bias averaged over
the observer window equals the body-z bias to 0.0004 m/s², and the true mean
vertical acceleration over that window is 0.004 m/s². The proxy's up-acceleration therefore carries a same-sign offset of
about −0.035 m/s², which is the gravity leak −g(1−cos ε) of its tilt error ε.
Its magnitude matches ε ≈ 5°. The proxy is a low-gain Mahony observer seeded
from a single accelerometer sample. That suits the wave-period estimator, which
ignores a DC offset, but not a bias measurement. With only about 25 s of data,
a wave-velocity residual Δv/T of ±0.01–0.04 m/s² is also comparable to the bias
itself.

## Answer to "does it help yaw?"

No. Startup yaw is unchanged to within 0.3 %, and the scored yaw is unchanged
when pooled. The one systematic yaw effect is a regression: the OU-III JONSWAP
4 m yaw gate fails.
