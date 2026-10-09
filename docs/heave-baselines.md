# Literature heave baselines

`src/heave_baselines/` implements two classic IMU heave estimators so the OU
filters can be compared against published methods on the same replay:

- `GodhavnHeaveFilter.h`: the adaptive heave filter of Godhavn (OCEANS'98),
  in the form analysed by Richter et al. (IFAC 2014).
- `KuchlerHeaveEKF.h`: the FFT-identified harmonic-mode EKF of Küchler et
  al. (IFAC 2011).

Both take the levelled vertical acceleration (z up, gravity removed) and
return heave. In `tests/heave_baselines` they share the Mahony attitude
front end of the PII observer (its tuned, sea-state-adaptive gains). The PII
observer runs alongside them as a third Mahony-levelled method.

## Godhavn 1998

Heave is the vertical acceleration through

    H(s) = s^2 / (s^2 + sqrt(2) w_c s + w_c^2)^2

This is a double integrator in series with a fourth-order high pass. It is
realised as two identical second-order sections, integrated with RK4. The
double zero at `s = 0` rejects a constant accelerometer bias. The cost is phase
lead in the wave band: `|1 - s^2 H| -> 2 sqrt(2) w_c / w_p`.

The cutoff adapts to the sea state. Every 5 s it minimises the error budget

    J(w_c) = (A^2/2) |1 - w_p^4 / D(j w_p)^2|^2 + q_w / (2^(7/2) w_c^3) + q_b / (2^(7/2) w_c^5)

The terms are:

- the error on a sinusoid of amplitude `A` at the dominant frequency `w_p`;
- white accelerometer noise of intensity `q_w`;
- an accelerometer bias random walk of intensity `q_b`.

Both noise integrals are exact for this `H`. `A` and `w_p` come from running
moments of the filter's own output. `q_w` comes from the input's first
difference. `q_b` is the sensor's datasheet figure. The applied cutoff relaxes
toward the target in `log w_c` with a 60 s time constant. A slow mean of the
input is subtracted ahead of the filter, so retuning does not excite the DC
state.

## Küchler et al. 2011

Heave is a sum of undamped modes `x_j = [z_j, z_j', w_j]` with
`z_j'' = -w_j^2 z_j`, where `w_j` is a random walk. A random-walk offset
state absorbs gravity and bias. The measurement is
`y = sum_j (-w_j^2 z_j) + x_off`. Propagation is the exact rotation for the
current `w_j`.

The steps are:

- **Identification.** The acceleration is block-averaged to 4 Hz into a
  128 s buffer. Every 30 s an FFT of the Hann-windowed buffer is divided by
  `w^2` to give the heave spectrum. Peaks above 25 % of the largest set the
  modes, up to 4.
- **Seeding.** New peaks are appended. Only their states are seeded, from a
  joint least-squares fit of the buffer at the peak frequencies. Existing
  modes keep their estimates.
- **Mode maintenance.** A mode that leaves 0.04–1 Hz is removed. Modes closer
  than 6 % in frequency are merged.
- **Start-up.** The EKF starts at the first identification. The offset is
  seeded from the fitted buffer constant.
- **Tuning.** `R` is the sensor datasheet white noise. Mode `j` gets
  velocity process noise `q_j = 4 zeta_q w_j^3 (A_j^2/2)`, with
  `zeta_q = 0.01`. This scales with the identified amplitude, so one setting
  serves every sea state. The frequency random walk is `0.002 w_j / sqrt(s)`.
  The scaling and the defaults are this implementation's choices, made by a
  sweep on the replay below. The paper reports neither.

`zeta_q` is the decisive tuning. With little process noise the filter locks
the frequency of a pure tone to 1 % (see `heave_baselines-test`). Larger
values let phase jumps explain the data instead, which biases `w_j` low.
Raising the mode count to 8 did not reduce the error on the replay.

## Comparison on the OU replay

`tests/heave_baselines/heave_baselines-sim` runs every method through
`process_wave_file_for_tracker`, the same noise template and seeds as the
OU-II, OU-III and TFG simulators. It uses:

- the 28 ft RAO dataset `v1.2.3`;
- an 8 mg accelerometer turn-on bias and a 5e-4 m/s²/√s bias random walk;
- a 25 Hz magnetometer.

It scores vertical displacement over the same trailing 900 s window. The
eight scored records give Z RMS error, as metres and as a percentage of Hs:

| Record | OU-III | OU-II | TFG | PII | Godhavn | Küchler |
|---|---|---|---|---|---|---|
| JONSWAP Hs 0.27 m | 1.1 cm (4.25%) | 1.6 cm (6.02%) | 1.2 cm (4.38%) | 1.0 cm (3.74%) | 1.8 cm (6.63%) | 1.3 cm (4.79%) |
| JONSWAP Hs 1.5 m | 6.1 cm (4.10%) | 9.9 cm (6.63%) | 6.3 cm (4.19%) | 8.4 cm (5.58%) | 15.6 cm (10.43%) | 16.2 cm (10.78%) |
| JONSWAP Hs 4 m | 16.1 cm (4.02%) | 25.4 cm (6.35%) | 16.1 cm (4.03%) | 25.2 cm (6.31%) | 34.9 cm (8.73%) | 56.0 cm (13.99%) |
| JONSWAP Hs 8.5 m | 31.5 cm (3.71%) | 51.4 cm (6.05%) | 31.7 cm (3.73%) | 65.5 cm (7.70%) | 125.1 cm (14.72%) | 137.5 cm (16.17%) |
| PM-Stokes Hs 0.27 m | 1.1 cm (4.21%) | 1.6 cm (5.81%) | 1.2 cm (4.33%) | 1.0 cm (3.83%) | 1.7 cm (6.45%) | 1.4 cm (5.27%) |
| PM-Stokes Hs 1.5 m | 6.0 cm (4.02%) | 9.6 cm (6.43%) | 6.1 cm (4.09%) | 8.3 cm (5.51%) | 14.8 cm (9.87%) | 17.9 cm (11.91%) |
| PM-Stokes Hs 4 m | 15.8 cm (3.96%) | 24.8 cm (6.20%) | 15.9 cm (3.98%) | 25.2 cm (6.31%) | 36.1 cm (9.03%) | 60.7 cm (15.17%) |
| PM-Stokes Hs 8.5 m | 31.5 cm (3.70%) | 50.8 cm (5.98%) | 31.7 cm (3.73%) | 76.0 cm (8.94%) | 152.5 cm (17.94%) | 156.6 cm (18.42%) |
| **Mean** | **4.00%** | **6.18%** | **4.06%** | **5.99%** | **10.48%** | **12.06%** |
| Worst | 4.25% | 6.63% | 4.38% | 8.94% | 17.94% | 18.42% |
| RMS error / OU-III (geometric mean) | 1.00× | 1.55× | 1.02× | 1.44× | 2.48× | 2.74× |

The adapted Godhavn cutoff runs from 0.20 rad/s (Hs 0.27 m) down to
0.047 rad/s (Hs 8.5 m). Its adaptive mean of 10.5 % beats every fixed cutoff
tried:

| Fixed `w_c` (rad/s) | 0.03 | 0.05 | 0.08 | 0.12 | 0.20 |
|---|---|---|---|---|---|---|
| Mean error (% Hs) | 128 | 35.4 | 15.4 | 12.0 | 15.3 |

Küchler's mean ranged from 12.1 % to 17.4 % over `zeta_q` 2e-4 to 3e-2 and
frequency random walk 5e-4 to 2e-3.

Both baselines are closest to OU-III on the smallest, shortest-period sea.
Their error grows with wave period. For Godhavn this is the fixed-structure
trade between phase lead and low-frequency noise. For Küchler it is a
broadband sea forced into a few undamped modes, with re-identification
transients between them.

Reproduce:

```bash
make ensure-sim-data
cd tests/heave_baselines && make build && ./heave_baselines-sim [--method godhavn|kuchler|pii]
```

`HB_GODHAVN_FIXED_WC`, `HB_KUCHLER_ZETA_Q` and `HB_KUCHLER_OMEGA_RW` override
the defaults for sweeps. The simulator gates each method on its own
regression sentinel.

## References

- J. M. Godhavn, "Adaptive tuning of heave filter in motion sensor", OCEANS'98, 1998.
- M. Richter, K. Schneider, D. Walser, O. Sawodny, "Real-time heave motion
  estimation using adaptive filtering techniques", IFAC World Congress, 2014.
- S. Küchler, C. Pregizer, J. K. Eberharter, K. Schneider, O. Sawodny,
  "Real-time estimation of a ship's attitude", IFAC World Congress, 2011.
