# Engine-noise degradation study

## Purpose

Every other noise source in this repository is electronic: accelerometer and
gyroscope white noise, initial bias, bias random walk, and magnetometer
calibration error.  A cruising sailboat carrying this estimator spends a
substantial share of its passages under power, and an inboard auxiliary diesel
puts mechanical vibration into the same two sensors.

This study measures how far the deployed OU-II, OU-III, and TFG estimators
degrade under that vibration, and identifies which part of the sensor path
causes it.  It is a **sensor-path** study: the vessel's rigid-body response to
the sea, the wave records, and the magnetometer are all unchanged.

## The vibration model

`src/util/W3dSimCommon.{h,cpp}` gains an engine vibration model for the
archetypal auxiliary of a 35-45 ft cruising boat: a naturally aspirated
three-cylinder four-stroke diesel on flexible mounts, a 2.6:1 reduction gear,
and a three-blade fixed propeller.  It injects

- crank orders: the half order, the first order, the firing order `n_cyl/2`,
  and the firing harmonics `2x`, `3x`, `4x`;
- driveline lines: shaft rate, propeller blade rate, and twice blade rate;
- a broadband structural floor, high-passed at 5 Hz so the engine contributes
  little directly in the wave band;
- governor hunting (an Ornstein-Uhlenbeck speed deviation of 0.4 % RMS with a
  2 s correlation time, plus a slow periodic hunt) and cycle-to-cycle
  combustion variability, which widen each order from a line into a band;
- accelerometer vibration rectification, specified in mg per g^2;
- gyroscope angular vibration via an effective lever arm, and the gyroscope's
  own linear-acceleration sensitivity.

Amplitudes are physical rather than fitted.  A line of order `k` is an
inertial excitation growing as `f^2`, transmitted through a mount of natural
frequency `f_n` and damping ratio `zeta` with

`T(f) = sqrt((1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2))`, `r = f/f_n`,

so one overall gain, calibrated once so the hull broadband RMS equals
`level_mps2` at `reference_rpm`, fixes the level at every other speed.  The
speed dependence that falls out is the familiar one rather than a plain `rpm^2`
law: a diesel auxiliary is rough at idle, where the low orders sit near the
mount resonance and `T > 1`.  At the default settings the hull RMS runs from
0.381 m/s^2 at 800 rpm to 0.808 m/s^2 at 3200 rpm, a factor of 2.1 over a
factor of 4 in speed.

The sensor is given a finite two-pole anti-alias bandwidth ahead of the 200 Hz
sample rate.  Each line's amplitude is attenuated at its **true** frequency and
its phase is then advanced at the sample rate, so orders above Nyquist fold on
their own; at 2400 rpm the 4.5 and 6.0 orders (180 and 240 Hz) are recorded at
20 and 40 Hz, and at 2000 rpm the 6.0 order lands on the sample rate exactly
and folds to DC.

The model is **off unless `W3D_ENGINE_RPM` is set**, so every existing noise
realization in the repository is bit-identical.  It is installed in the shared
`process_wave_file_for_tracker` runner, so all three families see the same
vibration.  Configuration:

| Variable | Default | Meaning |
| --- | ---: | --- |
| `W3D_ENGINE_RPM` | off | engine speed; 0 or unset disables the model |
| `W3D_ENGINE_CYLINDERS` | 3 | four-stroke cylinder count; firing order is `n/2` |
| `W3D_ENGINE_GEAR_RATIO` | 2.6 | engine rev per shaft rev |
| `W3D_ENGINE_BLADES` | 3 | propeller blade count |
| `W3D_ENGINE_LEVEL_MPS2` | 0.60 | hull broadband RMS at the reference speed |
| `W3D_ENGINE_REFERENCE_RPM` | 2400 | speed at which the level is stated |
| `W3D_ENGINE_MOUNT_HZ` | 10 | mount natural frequency |
| `W3D_ENGINE_MOUNT_ZETA` | 0.12 | mount damping ratio |
| `W3D_ENGINE_BANDWIDTH_HZ` | 80 | sensor anti-alias bandwidth |
| `W3D_ENGINE_GYRO_LEVER_M` | 1.5 | lever arm for hull angular vibration |
| `W3D_ENGINE_GYRO_G_SENS` | 1.78e-4 | gyro g-sensitivity, (rad/s)/(m/s^2) |
| `W3D_ENGINE_VRE_MG_PER_G2` | 1.0 | accelerometer vibration rectification |
| `W3D_ENGINE_SEED` | 20260828 | phase and modulation seed |

The simulator prints an `ENGINE_VIBRATION` banner with the resulting hull and
recorded RMS, followed by one `ENGINE_LINE` per harmonic giving its true
frequency, its aliased frequency, and its recorded amplitude.

### Is 0.60 m/s^2 realistic?

It is about 0.061 g.  At the 30-60 Hz carrying most of the energy that is
roughly 2 mm/s RMS velocity, which is *below* the 4 mm/s that ISO 6954 treats
as the threshold for adverse comment in accommodation spaces.  The nominal
point is therefore a routine cruise condition for a sensor on a cabin
bulkhead, not a worst case.  The level sweep runs from 0.15 m/s^2 (a quiet,
well-isolated installation) to 2.40 m/s^2 (a sensor near the engine bed).

## Protocol

`tools/engine_noise_degradation.py` replays the eight versioned stationary
JONSWAP and PM-Stokes cases from `oceanography-waves-lib` release `v1.1.3`, at
`H_s = {0.27, 1.5, 4.0, 8.5} m`, through each of OU-II, OU-III, and TFG.  The
ordinary sensor noise models stay on and the engine vibration is added on top;
the filters keep their deployed covariances, adaptation, pseudo-measurements,
startup, and regularization.  456 replays run across five arms and a common
engine-off baseline:

| Arm | What varies | Held fixed |
| --- | --- | --- |
| `baseline` | engine off | - |
| `speed` | 800 to 3200 rpm | level 0.60, bandwidth 80 Hz |
| `level` | 0.15 to 2.40 m/s^2 | 2400 rpm, bandwidth 80 Hz |
| `bandwidth` | 20 to 160 Hz | 2400 rpm, level 0.60 |
| `matched` | 20 to 160 Hz, level rescaled to a common recorded RMS | 2400 rpm |
| `path` | gyroscope coupling switched off | 2400 rpm, level 0.60, 80 Hz |

The `matched` arm is the control that separates *how much* vibration is
recorded from *where* the folded orders land: the recorded RMS is equalized at
0.454 m/s^2 while the fold frequencies still differ from cell to cell.  The
`path` arm attributes the effect between the two sensors.

Scoring uses the trailing 900 s of each 1200 s replay.  Pooling across the
eight equal-duration records is `RMS_pool = sqrt(mean(RMS_i^2))`, exactly the
RMS over their concatenation.

The window metrics gained a **mean error** alongside the RMS, emitted on a
separate `VALIDATION_METRICS_MEANS` line.  An RMS cannot distinguish a
zero-mean fluctuation from a constant offset, and here the offset is the whole
story.  The line is deliberately separate from `VALIDATION_METRICS`, which is
parsed into the committed validation and robustness evidence: adding columns
there would rewrite files the evidence fingerprint covers without any of their
numbers having changed.  Angles use a circular mean so a yaw error near the
+/-180 wrap does not average to zero.

## Results

### At the nominal cruise condition

| Family | Engine | 3-D [m] | 3-D offset [m] | Z [m] | Pitch RMS [deg] | Pitch offset [deg] | Yaw [deg] | 3-D / baseline |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.642 | 0.178 | 0.279 | 0.196 | 0.162 | 0.67 | 1.00 |
| OU-II | on | 3.363 | 2.862 | 0.740 | 2.436 | 2.303 | 58.59 | 5.24 |
| OU-III | off | 0.522 | 0.140 | 0.179 | 0.158 | 0.144 | 0.63 | 1.00 |
| OU-III | on | 4.232 | 3.875 | 0.828 | 3.084 | 3.075 | 45.38 | 8.10 |
| TFG | off | 0.745 | 0.206 | 0.190 | 0.386 | 0.353 | 0.80 | 1.00 |
| TFG | on | 14.026 | 12.903 | 1.610 | 16.205 | 16.127 | 74.35 | 18.84 |

"Offset" is the systematic part: the RMS across records of each record's mean
error.

All three families degrade heavily.  OU-III has the lowest engine-off error but
OU-II is lower at the cruise condition, so **the margin OU-III holds over OU-II
in the clean case does not survive vibration**.  Both OU families retain their
margin over TFG, whose cruise error is more than three times OU-III's and whose
degradation factor is 18.8 against 8.1.

### The error becomes an offset

With the engine off, OU-III's 3-D error is mostly fluctuation: 0.140 m of
0.522 m, 27 %, is systematic.  At cruise it is 3.875 m of 4.232 m, 92 %.  The
attitude channel is sharper still: pitch RMS 3.084 deg against a pitch offset
of 3.075 deg, i.e. a static tilt with essentially no fluctuation on it.
Vibration does not add noise to the estimate, it moves the estimate's
operating point.

The mechanism follows.  The accelerometer update is a joint measurement of
attitude and world acceleration, with predicted specific force
`R_wb (a_w - g)` and a covariance sized for a quiet sensor.  Wide-band
vibration cannot be attributed to `a_w`, whose OU prior only has wave-band
bandwidth, so the residual is taken up by attitude and by the accelerometer
bias state.  The resulting static tilt leaks gravity into the strapdown
horizontal acceleration (3 deg of tilt is 0.51 m/s^2), which the
integral-regularized displacement estimator turns into the bounded static
displacement offset that dominates the table.  OU-III's accelerometer
bias-state error correspondingly rises from 0.046 to 0.469 m/s^2.

Two candidates are ruled out by the same evidence:

- **Vibration rectification** contributes an offset of 2.3e-5 m/s^2, five
  orders of magnitude below the bias-state error.
- **The gyroscope path** contributes nothing measurable: with the model's
  angular-vibration and g-sensitivity terms switched off, every family is
  unchanged to four significant figures (OU-III 4.2318 vs 4.2316 m).  The
  degradation is entirely an accelerometer-path effect.

### There is no quiet regime

At 800 rpm OU-III already carries a 2.014 deg tilt offset against 0.144 deg
with the engine off.  Even at the quietest level swept - 0.15 m/s^2 hull,
0.114 m/s^2 recorded - the offset is 1.957 deg, more than thirteen times the
engine-off value.  The tilt offset appears almost fully formed as soon as any
vibration is present and then grows slowly; what grows quickly with level is
the displacement offset it drives.

The level sweep has a knee rather than a gentle trend.  Between 0.60 and
1.20 m/s^2 of hull vibration OU-III goes from 8.1x to 72x its baseline
(37.5 m), which is loss of lock rather than degradation.

### Power, not fold placement

| Family | Bandwidth as-is, 3-D [m] | Spread | Matched power, 3-D [m] | Spread |
| --- | --- | ---: | --- | ---: |
| OU-II | 0.87 / 1.89 / 3.36 / 8.21 | 9.4x | 2.47 / 4.21 / 3.36 / 4.79 | 1.9x |
| OU-III | 0.79 / 2.20 / 4.23 / 7.57 | 9.6x | 2.30 / 4.55 / 4.23 / 3.78 | 2.0x |
| TFG | 1.57 / 6.24 / 14.03 / 110.39 | 70.1x | 7.91 / 108.94 / 14.03 / 11.40 | 13.8x |

Cells are 20 / 40 / 80 / 160 Hz anti-alias bandwidth.

Equalizing the recorded vibration RMS collapses an 8:1 bandwidth range from a
factor of ~9.6 to a factor of ~2 for both OU families, while the folded line
frequencies still differ across cells.  Degradation is therefore set by **how
much out-of-band accelerometer power survives to the sample**, not by where the
aliased orders land - including the 2000 rpm cell where an order folds to
within a fraction of a Hz of DC, which does not stand out in the speed sweep.

TFG's matched cells do not collapse, but TFG has diverged in one of them
(>100 m), so its spread there is a statement about divergence, not about the
mechanism.

The practical consequence is that anti-alias filtering and mechanical isolation
act on the quantity that matters, while choosing a cruising rpm to move the
aliases does not.  Narrowing the sensor bandwidth from 80 to 20 Hz takes OU-III
from 4.232 m to 0.791 m.

## Interpretation boundary

This is a sensor-path result and is not the full cost of motoring.  The
magnetometer is unperturbed, the wave records are the same stationary seas used
elsewhere, and the vessel's rigid-body motion does not change; a real passage
under power also alters the encounter spectrum, adds propeller-induced surge,
and runs the engine at a speed that itself varies with the sea.  The estimator
is also not retuned for the condition, which is the point: these numbers are
what the shipped configuration does when the boat starts its engine.

They indicate that vibration-aware measurement covariance, or a motoring
detector gating the accelerometer update, is the natural next piece of work
rather than any refinement of the wave model.

## Reproducing

```
make -C tests/kalman_ou_ii kalman_ou_ii-sim
make -C tests/kalman_ou_iii kalman_ou_iii-sim
make -C tests/kalman_tfg kalman_tfg-sim
python3 tools/engine_noise_degradation.py --data-dir <wave-csv-dir> --jobs 4
```

Outputs land in `reports/results/engine_noise_degradation/`:
`engine_noise_runs.csv` (every cell), `engine_noise_summary.csv` (pooled),
`engine_noise_report.md`, `manifest.json`, and the two SVG figures, which are
mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.  The figures
use a fixed Matplotlib hash salt and carry no creation timestamp, so repeated
runs on the same evidence are byte-identical; `--no-plots` skips them where
Matplotlib is unavailable.
