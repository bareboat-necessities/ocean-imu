# OU-III horizontal wave displacement and velocity

The AtomS3R OU-III sketch and the replay simulator use the same causal horizontal wave estimate. In irregular waves the displacement primitive $S$ and velocity $v$ approximately satisfy $S + \beta v = 0$. For a sinusoid of angular frequency $\omega$, $\beta=1/\omega^2$. The shipping implementation estimates $\beta$ from the ratio of horizontal displacement and velocity variances over eight encounter periods. It uses the observed variance of $S+\beta v$ as the observation noise. The readout starts blending after three periods of live accelerometer-bias service and reaches full weight after four. The two-dimensional solve uses the current OU-III covariance, including any intervening magnetic correction. The recursive state and its covariance continue to supply the raw navigation estimate.

The eight JONSWAP and PM-Stokes records in each direction run at 200 Hz for 1,200 seconds. Each late RMS uses all 180,000 samples from 300 to 1,200 seconds. Pooled RMS is the square root of the mean squared per-record RMS. The alternate draws change the IMU noise and initial biases on the same physical wave records.

| Paired runs | Records | Parent 3-D RMS | OU-III 3-D RMS | Reduction |
| --- | ---: | ---: | ---: | ---: |
| Standard draw, both gyro profiles | 32 | 0.43360 m | 0.34928 m | 19.45% |
| Six draws, default gyro profile | 96 | 0.42728 m | 0.34472 m | 19.32% |
| Two reserved draws, default gyro profile | 32 | 0.43288 m | 0.34893 m | 19.39% |
| Rough seas and engine vibration | 8 | 0.87930 m | 0.71145 m | 19.09% |

The first two rows share 16 runs. All 152 distinct paired runs improved late 3-D displacement. On the two standard gyro profiles, horizontal X RMS fell from 0.29933 to 0.22766 m and Y RMS from 0.24946 to 0.18437 m. On the 16 default-profile records the 3-D velocity RMS fell from 0.07822 to 0.06608 m/s. Both profiles matched the parent on heave, roll, pitch, yaw and both bias RMS values. The 40–60 s startup displacement RMS also matched the parent.

All 32 standard gate decisions passed the previous limits. The additional seeds and engine tests retain their original gate outcomes in the CSV. Sixteen records with a zero-width, single-frequency observation gave a larger displacement reduction of 29.33%. The deployed readout allows finite bandwidth because an irregular sea contains several encounter frequencies.

Across the repository’s low–high–low wave transition, late 3-D RMS fell from 0.40659 to 0.35045 m in following seas and from 0.33329 to 0.26508 m in head seas. Both estimates improved in all seven scored phases. These runs use wave seed 7753, IMU seed 1063, initialization seed 3559 and the existing 1.5 → 4 → 1.5 m transition generator.

## Reproduction

Build the OU-III simulator at both revisions listed in the receipt using the repository replay compiler and the pinned v1.4.0 CSV inputs. Set `OCEAN_IMU_SIM_CASE` to the matching direction, `W3D_VALIDATION_WINDOW_SEC=900`, `W3D_VALIDATION_SEGMENTS=early:40:60,release:60:150`, `W3D_COLLECT_ALL_GATES=1` and `W3D_WRITE_TIMESERIES=0`. Run `kalman_ou_iii-sim --input <record>` for the default gyro profile and add `--regression` for the legacy profile. For velocity, enable time series and score the synchronized `vel_est_{x,y,z}` and `vel_ref_{x,y,z}` columns over the same trailing 900 seconds.

The six draws use the default seed and `W3D_SEED` values 7, 41, 101, 203 and 509. The two reserved seeds are 1063 and 3559. The engine cases use seed 203, both rough-sea spectra in both directions, `W3D_ENGINE_RPM=2400`, `W3D_ENGINE_BANDWIDTH_HZ=80` and `W3D_ENGINE_LEVEL_MPS2` of 0.6 or 1.2.

[Per-record measurements](paired-rms.csv) · [Source and executable hashes](receipt.json) · [Compiled replay limits](../../../../tests/common/CurrentReplayLimits.h)
