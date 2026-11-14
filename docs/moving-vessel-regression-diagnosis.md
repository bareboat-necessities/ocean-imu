# Why the 5 kn replay inputs raise every filter's error

The oceanography-waves-lib v1.3.0 replay archive
(`sim-data-files-vessel-rao-28ft-5kn.zip`, SHA-256 `324e28cf…82baf3`) moves the
28 ft hull at 5 kn with heading `3° sin(2πt/120 s)`; v1.2.3
(`sim-data-files-vessel-rao-28ft.zip`, SHA-256 `87e3cee9…908d433c`) holds it on
station in the same seas. With unchanged estimators, heave, horizontal, attitude
and bias errors rise for every filter, most in the 0.27 m seas. This note
attributes that change. All numbers come from
`tools/moving_vessel_diagnosis.py` and are stored in
`reports/results/moving_vessel_diagnosis/`.

## Answer

1. **Every scored record is a following or quartering sea.** All eight records
   propagate within 30° of the +x track (elevation is in phase with +x orbital
   velocity in both archives). At 5 kn the encounter period lengthens: heave
   zero-crossing period ×1.26–1.60, vertical acceleration RMS ÷1.7–2.4, heave
   RMS −9 % to +3 %. The same metre of heave arrives with roughly half the
   acceleration.
2. **Steering is not a cause.** Adding the exact 5 kn / ±3° track to the
   stationary records, or removing it from the moving ones, changes heave error
   by at most ±0.25 % and bias error by −8 % to +4 %; horizontal error grows
   1–7 %.
3. **The lower encounter frequency is the cause.** Slowing each stationary record
   to the moving record's zero-crossing period (`stretch`: same spectral shape,
   same heave, accelerations ÷s²) reproduces the heave regression of every
   filter. It is at least as large as the moving one in 33 of 40
   filter/record cells; per filter, the median of `(stretch − S)/(M − S)` is
   1.4–1.9.
4. **It is an estimator property, not information loss.** A non-causal Wiener
   smoother on the simulated accelerometer reaches 0.11–0.26 % of H on the moving
   records, an order of magnitude below every filter. The regression remains
   with sensor noise and bias switched off.
5. **Two estimator mechanisms carry it** (quantified below). The fixed
   0.28 m/s² stillness detector reads the 0.27 m following seas as still water;
   it causes about 85 % of OU-II's worst-record regression and most of CARMA's
   from 1.5 m up. Every wave model's normalized heave error also grows with wave
   period at fixed heave; this carries most of the OU-III, TFG and shared-mode
   regression. TFG adds a tilt/bias excursion on long-period encounter spectra.
6. **Scoring at the CSV time stamp adds 0.43 cm RMS x-error to every record.**
   The time column has six significant digits, so samples after 1000 s are
   stamped to 0.01 s. An underway scorer that evaluates the declared track at
   that stamp misplaces it by up to 5 ms × 2.57 m/s. Evaluating it at
   `sample_index × dt` removes the artifact.

## Input physics, no estimator

| Record | heave T_z (s) | vertical acc RMS (m/s²) | heave RMS (m) | still (%) | Wiener bound (% H) |
| --- | ---: | ---: | ---: | ---: | ---: |
| JONSWAP 0.27 m | 2.97 → 4.76 | 0.246 → 0.106 | 0.053 → 0.049 | 41 → 96 | 0.11 → 0.25 |
| JONSWAP 1.5 m | 4.93 → 6.96 | 0.731 → 0.312 | 0.391 → 0.358 | 2 → 25 | 0.07 → 0.11 |
| JONSWAP 4 m | 6.83 → 8.78 | 1.101 → 0.583 | 1.043 → 0.976 | 0 → 4 | 0.06 → 0.12 |
| JONSWAP 8.5 m | 8.81 → 11.13 | 1.437 → 0.825 | 2.155 → 2.213 | 1 → 4 | 0.07 → 0.13 |
| PM-Stokes 0.27 m | 2.96 → 4.46 | 0.236 → 0.119 | 0.050 → 0.047 | 43 → 92 | 0.11 → 0.26 |
| PM-Stokes 1.5 m | 4.63 → 6.69 | 0.825 → 0.340 | 0.386 → 0.354 | 1 → 18 | 0.07 → 0.11 |
| PM-Stokes 4 m | 6.33 → 8.40 | 1.290 → 0.661 | 1.040 → 0.991 | 0 → 2 | 0.07 → 0.11 |
| PM-Stokes 8.5 m | 8.12 → 10.53 | 1.706 → 0.933 | 2.168 → 2.206 | 0 → 3 | 0.06 → 0.14 |

`still` is the fraction of the scored 900 s the shared `StillnessAdapter`
(`tuner/SeaStateFusionTunerCommon.h`) flags as still water: energy EMA of the
6 Hz-filtered vertical acceleration below `8e-4 g²` (0.28 m/s²) for 2 s. It is a
threshold on acceleration, so the same heave at half the acceleration reads as
still water. Both archives are internally consistent to their print
precision: body specific force matches `C_wb (a + g)` within 6e-5 m/s² and gyro
matches the attitude derivative within 2.3e-4 rad/s (`inputs.json`).

## Factor isolation

Five arms per record, replayed by the same executable with default seeds and
noise, scored on the trailing 900 s against double-precision input truth.
Counterfactual body specific force and rate are rebuilt from world kinematics,
so every arm is exactly self-consistent (tested in
`tests/validation/test_moving_vessel_diagnosis.py`).

| Arm | Content |
| --- | --- |
| `S` | v1.2.3 stationary record |
| `track` | `S` plus the declared track and heading: steering without Doppler |
| `stretch` | `S` slowed to the moving record's `T_z`: Doppler shift without steering or reshaping |
| `notrack` | `M` minus the declared track and heading: encounter spectrum without steering |
| `M` | v1.3.0 5 kn record |

Heave RMS error, % of nominal H:

| Filter / record | S | track | stretch | notrack | M |
| --- | ---: | ---: | ---: | ---: | ---: |
| **OU-II** worst of 8 | 6.63 | 6.63 | 23.24 | 20.85 | 20.85 |
| OU-II JONSWAP 0.27 m | 6.02 | 6.02 | 23.24 | 20.85 | 20.85 |
| OU-II PM-Stokes 0.27 m | 5.81 | 5.81 | 18.99 | 11.31 | 11.31 |
| OU-II JONSWAP 1.5 m | 6.63 | 6.63 | 8.35 | 7.59 | 7.60 |
| OU-II JONSWAP 8.5 m | 6.05 | 6.05 | 7.42 | 6.86 | 6.86 |
| **OU-III** worst of 8 | 4.25 | 4.25 | 7.37 | 8.22 | 8.24 |
| OU-III JONSWAP 0.27 m | 4.25 | 4.25 | 7.37 | 8.22 | 8.24 |
| OU-III PM-Stokes 0.27 m | 4.21 | 4.22 | 6.70 | 7.83 | 7.84 |
| OU-III JONSWAP 1.5 m | 4.11 | 4.11 | 6.09 | 5.50 | 5.50 |
| OU-III JONSWAP 8.5 m | 3.71 | 3.71 | 5.01 | 4.40 | 4.40 |
| **TFG** worst of 8 | 4.38 | 4.39 | 7.22 | 7.81 | 7.82 |
| TFG JONSWAP 0.27 m | 4.38 | 4.39 | 7.22 | 7.81 | 7.82 |
| TFG PM-Stokes 0.27 m | 4.33 | 4.34 | 6.77 | 7.47 | 7.48 |
| TFG JONSWAP 1.5 m | 4.19 | 4.20 | 6.04 | 5.52 | 5.52 |
| TFG JONSWAP 8.5 m | 3.73 | 3.73 | 4.99 | 4.40 | 4.40 |
| **CARMA** worst of 8 | 3.21 | 3.21 | 4.45 | 4.67 | 4.67 |
| CARMA JONSWAP 0.27 m | 3.06 | 3.05 | 4.30 | 4.19 | 4.19 |
| CARMA PM-Stokes 0.27 m | 3.21 | 3.21 | 4.08 | 4.67 | 4.67 |
| CARMA JONSWAP 1.5 m | 2.25 | 2.25 | 4.12 | 3.31 | 3.32 |
| CARMA JONSWAP 8.5 m | 2.09 | 2.09 | 2.71 | 3.00 | 3.00 |
| **Shared QMEKF** worst of 8 | 2.45 | 2.45 | 3.73 | 3.70 | 3.70 |
| Shared QMEKF JONSWAP 0.27 m | 2.45 | 2.45 | 3.73 | 3.55 | 3.55 |
| Shared QMEKF PM-Stokes 0.27 m | 2.40 | 2.40 | 3.46 | 3.70 | 3.70 |
| Shared QMEKF JONSWAP 1.5 m | 2.19 | 2.19 | 2.99 | 2.60 | 2.60 |
| Shared QMEKF JONSWAP 8.5 m | 2.02 | 2.02 | 2.53 | 2.26 | 2.25 |

`track ≈ S` and `notrack ≈ M` for every filter and record; `stretch` carries
the regression. Horizontal error (m RMS, exact track):

| Filter / record | S | track | stretch | notrack | M |
| --- | ---: | ---: | ---: | ---: | ---: |
| **OU-II** worst of 8 | 0.765 | 0.786 | 0.929 | 0.884 | 0.898 |
| **OU-III** worst of 8 | 0.665 | 0.689 | 0.916 | 0.802 | 0.819 |
| **TFG** worst of 8 | 0.674 | 0.683 | 0.945 | 0.817 | 0.832 |
| **CARMA** worst of 8 | 0.339 | 0.352 | 0.439 | 0.360 | 0.369 |
| **Shared QMEKF** worst of 8 | 0.310 | 0.326 | 0.386 | 0.293 | 0.310 |

Accelerometer-bias error (m/s² RMS) / worse of roll and pitch (° RMS) for the
two outlier groups:

| Filter / record (bias m/s² / tilt °) | S | track | stretch | notrack | M |
| --- | ---: | ---: | ---: | ---: | ---: |
| OU-II JONSWAP 0.27 m | 0.046 / 0.23 | 0.045 / 0.22 | 0.371 / 2.12 | 0.133 / 0.59 | 0.132 / 0.59 |
| OU-II PM-Stokes 0.27 m | 0.045 / 0.22 | 0.045 / 0.22 | 0.254 / 1.48 | 0.069 / 0.36 | 0.068 / 0.35 |
| TFG JONSWAP 8.5 m | 0.050 / 0.31 | 0.049 / 0.31 | 0.052 / 0.26 | 0.110 / 0.53 | 0.106 / 0.53 |
| TFG PM-Stokes 8.5 m | 0.057 / 0.35 | 0.056 / 0.35 | 0.047 / 0.30 | 0.162 / 0.69 | 0.157 / 0.67 |
| TFG JONSWAP 4 m | 0.047 / 0.28 | 0.046 / 0.27 | 0.106 / 0.55 | 0.049 / 0.29 | 0.053 / 0.30 |

## Mechanisms

Diagnostic runs on the four JONSWAP records; `S`, `M` and `stretch` arms.

| Filter | Run | 0.27 S | 0.27 M | 0.27 stretch | 1.5 S | 1.5 M | 1.5 stretch | 8.5 S | 8.5 M | 8.5 stretch |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | `shipping` | 6.02 | 20.85 | 23.24 | 6.63 | 7.60 | 8.35 | 6.05 | 6.86 | 7.42 |
| OU-II | `stillness_off` | 6.14 | 8.48 | 8.47 | 6.64 | 7.63 | 8.40 | 6.05 | 6.87 | 7.41 |
| OU-II | `no_noise` | 5.14 | 26.20 | 24.28 | 5.37 | 6.34 | 6.92 | 5.22 | 5.99 | 6.43 |
| OU-III | `shipping` | 4.25 | 8.24 | 7.37 | 4.11 | 5.50 | 6.09 | 3.71 | 4.40 | 5.01 |
| OU-III | `stillness_off` | 4.33 | 7.74 | 7.19 | 4.15 | 5.50 | 5.99 | 3.70 | 4.36 | 4.97 |
| OU-III | `no_noise` | 2.96 | 7.62 | 6.69 | 3.02 | 4.44 | 4.96 | 2.94 | 3.47 | 4.18 |
| OU-III | `ou3_floors_off` | 4.32 | 7.75 | 7.20 | 4.15 | 5.50 | 5.99 | 3.70 | 4.36 | 4.97 |
| OU-III | `ou3_floors_off_no_noise` | 3.05 | 5.50 | 4.96 | 2.96 | 3.74 | 4.28 | 2.91 | 3.35 | 4.05 |
| OU-III | `ou3_floors_off_R_scaled_no_noise` | — | — | 4.41 | — | — | 4.12 | — | — | 4.01 |
| TFG | `shipping` | 4.38 | 7.82 | 7.22 | 4.19 | 5.52 | 6.04 | 3.73 | 4.40 | 4.99 |
| TFG | `stillness_off` | 4.38 | 7.82 | 7.22 | 4.19 | 5.52 | 6.04 | 3.73 | 4.40 | 4.99 |
| TFG | `no_noise` | 3.10 | 5.92 | 5.35 | 3.07 | 3.85 | 4.40 | 2.95 | 3.41 | 4.07 |
| CARMA | `shipping` | 3.06 | 4.19 | 4.30 | 2.25 | 3.32 | 4.12 | 2.09 | 3.00 | 2.71 |
| CARMA | `stillness_off` | 2.46 | 3.57 | 3.75 | 2.21 | 2.65 | 3.01 | 2.08 | 2.33 | 2.62 |
| CARMA | `no_noise` | 2.95 | 4.14 | 3.94 | 1.87 | 3.35 | 3.90 | 1.71 | 2.73 | 2.38 |
| Shared QMEKF | `shipping` | 2.45 | 3.55 | 3.73 | 2.19 | 2.60 | 2.99 | 2.02 | 2.25 | 2.53 |
| Shared QMEKF | `stillness_off` | 2.45 | 3.55 | 3.73 | 2.19 | 2.60 | 2.99 | 2.02 | 2.25 | 2.53 |
| Shared QMEKF | `no_noise` | 2.06 | 3.19 | 3.20 | 1.78 | 2.12 | 2.52 | 1.60 | 1.73 | 2.11 |

**Stillness detector (OU-II, CARMA).** In the 0.27 m seas the detector reports
still water 92–96 % of the moving window (41–43 % stationary). The operating
point then multiplies the wave variance by `exp(-t_still/decay)`
(`measureOperatingPoint_`; decay 5 s in OU-II, 1 s in OU-III), pinning
`sigma_a` at its floor. Disabling the detector (`stillness_off`) cuts OU-II's
0.27 m increase from 14.8 to 2.3 points of H, removes all of its 0.27 m
tilt/bias excursion (bias 0.132 → 0.044 m/s², as stationary), and removes
59–72 % of CARMA's increase at 1.5–8.5 m. It is worth 0.6 point to OU-III.
TFG and the shared-mode filter do not use this path and are unchanged.

**Wave-model period dependence (all filters).** Without sensor noise every
filter is near scale-invariant on the stationary records (OU-III 2.9–3.0 % of
H in all seas) yet degrades on `M` and `stretch`, so the increase is
deterministic. OU-III's absolute floors (wave-acceleration prior
`max(0.05 m/s², band noise floor)`, `ACC_NOISE_FLOOR_SIGMA` = 0.12 m/s²,
`r_S ≥ 0.15`, stillness) explain about half of its noise-free 0.27 m and
1.5 m increase (`ou3_floors_off_no_noise`), but with the simulated noise
removing them is worth no more than disabling stillness alone. The rest
persists with floors off and with the assumed accelerometer noise rescaled by
`1/s²` (`ou3_floors_off_R_scaled_no_noise`: about 4 % of H on `stretch`
against 3 % on `S`), so it comes from constants fixed in seconds or from the
wave prior itself. In the band holding most heave energy, OU-III
over-estimates heave amplitude by 14 % on the 0.27 m stationary record
(0.2–0.35 Hz) and by 39 % on the moving one (0.1–0.2 Hz); OU-II
under-estimates it (gain 0.76) while `sigma_a` is floored. The shared-mode and
CARMA filters degrade least.

**TFG long-period tilt/bias excursion.** On the 8.5 m records TFG's bias error
grows ×2–3 and its pitch error ×6–7 (0.06 → 0.36° and 0.09 → 0.66°) on `M`
and `notrack`, but not on `stretch` or `track`; the stretched 4 m record
triggers the same excursion. It needs long-period encounter content, not
steering.

## Consequences for replay gates

- Gates rebaselined on these records measure the following-sea case only. A
  head sea raises the encounter frequency and would move each mechanism the
  other way; the archive has no such record.
- The largest single deterioration, OU-II's 0.27 m heave and bias, is the
  stillness detector classifying a 0.2 m following sea as still water.
- Horizontal and 3D scores of the low seas include the 0.43 cm time-stamp
  artifact until the scorer evaluates the track at sample time. A phase
  surrogate that re-derives displacement from velocity, and removes, restores
  and scores the track at the same rounded stamps, does not carry it; a
  wave-height-scaled record that scales the residual displacement carries it
  multiplied by the scale.

## Next falsifiable experiments

1. A head-sea counterpart of the eight records: predicts stillness fractions
   and OU-II/CARMA errors at or below the stationary values.
2. A stillness threshold in displacement-equivalent units (acceleration
   divided by the tracked `ω²`): predicts OU-II's 0.27 m moving error near the
   `stillness_off` value.
3. OU-III residual period dependence: bisect the constants fixed in seconds on
   the noise-free, floors-off `stretch` arm; the target is a `stretch/S` ratio
   near 1.
4. TFG: replay the stretched 4 m record to isolate the tilt/bias separation
   path that fails on long-period seas.

## Reproduce

```sh
# Unpack both archives into v123/ and v130/, then build the replays, e.g.
make -C tests/kalman_ou_iii kalman_ou_iii-sim   # likewise ou_ii, tfg, carma3, shared modes
python3 tools/moving_vessel_diagnosis.py inputs --stationary v123/ --moving v130/ \
    --json reports/results/moving_vessel_diagnosis/inputs.json
python3 tools/moving_vessel_diagnosis.py arms --stationary v123/ --moving v130/ --out arms/
python3 tools/moving_vessel_diagnosis.py run --tests tests --arms-dir arms/ \
    --label shipping --out reports/results/moving_vessel_diagnosis/shipping.csv
python3 tools/moving_vessel_diagnosis.py table reports/results/moving_vessel_diagnosis/*.csv
```

The scorer passes `W3D_REFERENCE_MOTION=stationary` and removes the track
itself, so it runs unchanged with executables with or without an underway
reference. `exit_code` 1 records a gate miss against the limits compiled into
the executable; the scores are complete. `no_noise` adds `--extra=--no-noise`.
The other diagnostic labels
need an experiment-only build with two environment hooks, never part of a
shipping build:

```diff
--- a/src/tuner/SeaStateFusionTunerCommon.h   (StillnessAdapter constructor)
+        if (const char* v = std::getenv("XP_STILL_ENERGY_THRESH")) energy_thresh = std::strtof(v, nullptr);
--- a/src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h   (apply_ou_tune_)
-        const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());
+        static const float xp_floor = std::getenv("XP_OU3_SIGMA_FLOOR")
+            ? std::strtof(std::getenv("XP_OU3_SIGMA_FLOOR"), nullptr) : 0.05f;
+        const float sigma_floor = std::max(xp_floor, band_noise_floor_sigma_());
```

`stillness_off` sets `XP_STILL_ENERGY_THRESH=0`. `ou3_floors_off*` also sets
`XP_OU3_SIGMA_FLOOR=1e-4`, `OU_III_ACC_NOISE_FLOOR_SIGMA=1e-6` and
`OU_III_R_S_MIN=1e-4`; the `R_scaled` rows add `SF_SIGMA_A_SCALE=1/s²` per
record. With the hooks unset the experiment build reproduces the shipping
scores exactly on all five filters (`manifest.json` records both executable sets).
