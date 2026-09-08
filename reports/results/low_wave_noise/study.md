# Low-wave filter study

Full 1200-second records, final 900 seconds; unchanged sensor injection and quality gates. Fresh validation uses paired IMU/initialization seeds 31013, 37003, 41011 and 47017 on all eight pinned records. These are new sensor draws, not independent stochastic seas. See the manifest for binary/source provenance and invalid arms.

| Family | Weighting | Violations / 32 | Low X acceleration RMS, m/s² | Low Y acceleration RMS, m/s² | Roll, deg | Pitch, deg | Yaw, deg | 3-D displacement, m |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| OU-II | off | 70 | 0.009421 | 0.009486 | 0.255251 | 0.328171 | 1.923037 | 0.382182 |
| OU-II | SNR-dependent | 69 | 0.006875 | 0.006750 | 0.255076 | 0.330511 | 1.928230 | 0.382119 |
| OU-III | off | 81 | 0.007941 | 0.007976 | 0.214166 | 0.327243 | 1.923412 | 0.294180 |
| OU-III | SNR-dependent | 80 | 0.006370 | 0.006081 | 0.212823 | 0.330782 | 1.933432 | 0.293226 |

Acceleration columns average only the two lowest-wave records across four draws; attitude/displacement columns average all 32 cases. OU-II still fails 30/32 cases and OU-III 32/32. Small pitch/yaw tradeoffs remain. Reduced acceleration noise is not a claim that low-wave attitude or direction is solved.

The separate heading/angle validation uses seed 53017: eight rigid world rotations and eight recomputed hull responses at incident angles −60° and 75°, retaining the prescribed harmonic phases. Weighting leaves violation counts unchanged (OU-II 73/16 cases; OU-III 65/16), while small attitude tradeoffs persist.

Direction bias correction is paired within a single executable. It does not alter upstream attitude/displacement/tuning. Zero axes are unavailable; conditional RMS cannot be interpreted without availability. At the default low-wave draw, corrected axis availability averages about 1% and travel direction is unresolved more than 99% of the time.

Rejected defaults include direct RAO scaling of the OU prior, longer tau alone, cadence-matched tau changes, globally enlarged accelerometer covariance, asymmetric R_S factors, increased bias priors, and magnetic-weight changes. OU-III R_S factors 0.5 improve training displacement but add a fresh-draw violation and worsen mean pitch. Global OU-II horizontal sigma ×4 adds a fresh-draw bias violation. The SNR policy preserves the original weighting in stronger seas.

Yaw screens and the magnetic-error diagnostic show sensitivity to both hard-iron offset and soft-iron distortion; single-heading low-wave motion does not identify a full magnetic calibration. No new yaw default is retained.

Applied OU-III tau/sigma/R_S telemetry now reports the active model. Earlier OU-III CSV fields with those names reported staged commands; they cannot locate actual covariance synchronization events. Tuning-runner checks reject variables absent from the executable and record the starting producer/source/binary hashes.

The shipping O3 checks retain 6/8 passing records in each OU family: six violations in OU-II and two in OU-III, all in the two 8.5 m records. O3 simulator copies are byte-identical to the Makefile outputs. Across the default low-wave spectral checks, horizontal peak gain is 0.916–0.980 after weighting; phase magnitude stays below 2.8 degrees. Broadband error falls, but X wave-band error increases and peak gain loses roughly 1–2.5 percentage points. The full per-axis comparison is in spectral-comparison.csv. The PNG shows the JONSWAP Hs=0.27 m error spectra over the final 900 seconds.

The additional Hs 0.27→1.5→0.27 surrogate uses kinematically closed 120-second transitions and two fresh sensor draws. Neither family adds a quality violation (OU-II 5, OU-III 6 across two cases each). Low-return 3-D acceleration RMS falls about 15% in OU-II and 10% in OU-III; mean pitch increases, and OU-III displacement during the rising transition increases about 3%. These are transition tradeoffs, not uniform improvement or a new gate certification.
