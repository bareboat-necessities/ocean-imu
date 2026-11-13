# Filter simulator regression gates

The estimator baseline is main `a053a1d2e571b4e4f39c566107c35154db1f3242`
(v2.7.6). The filter settings, noise, data, score windows and estimator equations
are unchanged by this gate update. The pinned input archive is
oceanography-waves-lib v1.2.3, vessel-RAO 28 ft at CG, verified by
`tools/sim_dataset.py`.

`tests/common/CurrentReplayLimits.h` contains the production-default W3d limits.
Each channel uses the worst of its eight scored records, with a 0.5% relative
margin and upward rounding to 0.00001. An existing, stricter production limit
is retained. These are deterministic regression sentinels for the default
seed and configuration, not population confidence intervals or accuracy goals.
Accelerometer and gyroscope 3D bias have separate gates. TFG also gates roll
and pitch. Historical interval-rate `--regression` limits remain separate;
OU-II now supports that same comparison protocol as the other native filters.

The W3d scores use the trailing 900 s of each 1200 s, 200 Hz record. Vertical
displacement is RMS/Hs; 3D displacement is RMS/max true displacement norm;
accelerometer Z bias is RMS/max absolute true Z bias; each 3D bias is vector
RMS/max true bias norm. The report also retains absolute vector bias RMS.
These denominators differ: a gate reduction is not itself an RMS improvement.

All 13 `*-sim.cpp` files are covered:

| Simulator | Gate/check |
| --- | --- |
| `ahrs-qmekf-sim.cpp` | Attitude, trailing 60 s; existing throughput floor |
| `heave_only-sim.cpp` | Five methods with Mahony/proxy; four with ideal attitude; raw heave RMS |
| `kalman_ou_ii-sim.cpp` | Production W3d limits; historical interval-rate comparison |
| `kalman_ou_iii-sim.cpp` | Production W3d limits; historical interval-rate comparison |
| `kalman_tfg-sim.cpp` | Production W3d limits; historical interval-rate comparison |
| `kalman_tfg-legacy-sim.cpp` | Current LegacyCubic law with its own measured W3d limits |
| `carma3-sim.cpp` | Default CARMA(4,1), production/interval profiles; other research variants diagnostic |
| `kalman_carma_shared_modes-sim.cpp` | Shared-mode QMEKF, production/interval profiles |
| `tfg_carma3-sim.cpp` | TFG CARMA(4,1), production/interval profiles |
| `kalman_tfg_shared_modes-sim.cpp` | Shared-mode TFG, production/interval profiles and protocol rejection checks |
| `nlo-sim.cpp` | Raw reported heave; separate low-wave PM-Stokes limit |
| `carma_process-sim.cpp` | Prior generator checked independently by `carma_process_validation.py --check`: PSD, variance, derivative and increment laws |
| `carma_lever_arm-sim.cpp` | Study instrument over several placements/models; no fitted accuracy envelope across ablations; `carma_lever_arm-test` checks its numerical model |

Kuchler's PM-Stokes heave gates tighten from 18.96% to 17.50% Hs with its
Mahony frontend, and from 20.23% to 20.13% with ideal attitude. The other
heave-only and stochastic-process tolerances retain their qualified margin
where no new improvement supports tightening. These
measurements do not qualify envelopes for other seeds, engine vibration,
lever-arm ablations or alternate research priors.
The heave test runner verifies and explicitly selects the eight canonical
records, so missing inputs fail instead of silently reducing test coverage.

Scored nonfinite physical histories fail the shared W3d gates, including
nonfinite bias truth and incomplete histories. NLO rejects nonfinite heave
instead of excluding it from RMS; AHRS rejects nonfinite attitude RMS. Optional
unavailable direction telemetry is unaffected. Zero true bias or zero
displacement in quiet diagnostic data does not manufacture a percentage score.

For reproducible measurements after building:

```sh
make all
python3 tools/filter_sim_gates.py --run --jobs 3 \
  --log-dir /tmp/filter-sim-gates --output /tmp/filter-sim-gates.json
```

The report tool verifies the pinned input hashes and requires every scored
record exactly once. Missing, duplicated, partial and nonfinite W3d scores
fail its checks. It records compiler, source commit, binary hashes and paired
legacy/production scores, and never changes a limit after a failed replay.
The versioned compact measurements are under
`reports/results/filter_sim_gates/`; dense output and runner logs are omitted.

To reuse the full native test run instead of running the same replays twice,
capture `W3D_WRITE_TIMESERIES=0 W3D_VALIDATION_WINDOW_SEC=900 make all` and pass
its output to the report tool with `--make-log LOG --log-dir /tmp/filter-gates`.
The importer keeps legacy and current sampling profiles separate and rejects
failed or incomplete numeric replays.

## Current 3D bias measurements

Absolute bias RMS uses the vector error over each trailing window. The table
shows the worst of the eight records. The legacy column is the interval-rate
sampling comparison on the same current estimator, inputs and seeds; it is not
an older main build. Positive reduction means lower error.

| Configuration | Legacy RMS (mm/s²) | Current RMS (mm/s²) | Reduction | Current accel 3D gate (% max true bias) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 52.128 | 52.041 | 0.17% | 61.99132 |
| OU-III | 58.632 | 58.014 | 1.05% | 69.10592 |
| TFG | 56.479 | 57.380 | -1.60% | 68.35106 |
| CARMA | 69.193 | 68.615 | 0.84% | 81.73414 |
| Shared-mode QMEKF | 70.501 | 69.911 | 0.84% | 83.27692 |
| TFG-CARMA | 69.281 | 68.692 | 0.85% | 81.82580 |
| Shared-mode TFG | 70.018 | 69.414 | 0.86% | 82.68527 |

The worst 3D bias improves modestly for six of the seven default native
families. TFG worsens in this matched sampling comparison, despite its
historical bias gate being much looser than either measured result. The
substantial gap between older gate limits and current scores must not be
attributed entirely to the recent sampling correction. Residual horizontal
accelerometer-bias error remains substantial.

AHRS maxima are roll 0.932953°, pitch 5.3324°, yaw 12.6579°; their new gates
are 0.9377°, 5.3591°, 12.7212°. NLO raw heave maxima are 6.93681% Hs
(JONSWAP), 7.32231% (low-wave PM-Stokes), and 6.98315% (other PM-Stokes);
the new gates are 6.98%, 7.36%, and 7.02%.
