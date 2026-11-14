# Moving-vessel regressions and the stillness fix

The default replay dataset is the checksum-pinned oceanography-waves-lib
v1.3.1 28 ft vessel at 5 kn. Zero-speed controls also use v1.3.1, with exact
200 Hz timestamps. No v1.3.0 result is retained in this report.

## Low-wave failure mechanism

The shared stillness detector tests filtered vertical acceleration energy
against an absolute threshold, equivalent to about 0.28 m/s² RMS. Its
`isStill()` flag changes immediately; the two-second hold controls frequency
relaxation only. A small, long-period wave can fall below that acceleration
threshold while retaining substantial displacement. For sinusoidal motion,
acceleration amplitude is displacement amplitude times angular frequency
squared, so lower encounter frequency makes this ambiguity worse.

OU-II and the CARMA Legacy wrapper previously multiplied their estimated
wave variance by an exponential of continuous below-threshold time. On the
0.27 m moving seas this suppressed real wave energy and drove the amplitude
toward its floor. Both wrappers already subtract a sensor-noise floor in
the variance channel.

OU-II now disables this additional attenuation by default.
`setSigmaStillnessDecaySec(positive_seconds)` retains the historical study
setting; zero disables it again. Negative and nonfinite settings are ignored.
CARMA no longer passes stillness into variance attenuation. TFG-CARMA uses
that same wrapper. Frequency relaxation and the noise-floor subtraction
remain in place. The shared-mode Evidence excitation path, OU-III and TFG
retain their existing estimator behavior.

## Exact-clock measurements

The paired record-level values and input hashes are in
[`stillness_comparison.json`](../reports/results/filter_sim_gates/stillness_comparison.json).
The comparison uses the same v1.3.1 moving input records and trailing 900 s.
Values below are deterministic maxima over the eight default-seed records;
separate maxima need not come from the same record.

| Filter | Worst heave before → after (% Hs) | 0.27 m JONSWAP before → after (% Hs) |
| --- | ---: | ---: |
| OU-II | 20.84 → 8.48 | 20.84 → 8.48 |
| CARMA | 4.67 → 3.73 | 4.19 → 3.57 |
| TFG-CARMA | 4.67 → 3.73 | 4.19 → 3.57 |

The [zero-speed versus moving comparison](moving-vessel-comparison.md) uses
the corrected wrappers in both arms. OU-II worst zero-speed heave is about
6.65% Hs, slightly above the old 6.63%; CARMA and TFG-CARMA improve from about
3.21% to 2.46%. The complete position, attitude and bias measurements remain
visible; the fix does not improve every metric on every record. Production
and interval-rate limits are remeasured using the documented rule and the
previous production ceilings.

## Encounter-period phase correction

The PII library now retains the empirical noise-smoothing floor and adds an
analytical phase match between the acceleration filter and its PII restoring
core. Across all eight moving records and three sensor seeds, all 24 pairs
improve: mean heave 12.13 → 8.73% Hs; worst 15.39 → 10.49%. The zero-speed
worst case is unchanged; the largest small-wave cost is 0.316 mm RMS.
See the [scaling study](encounter-scaling-regressions.md) for the derivation,
per-record comparisons and rejected cutoff/clamp experiments.

## Remaining regressions

OU-III, TFG, Adaptive PII and NLO still show worse normalized heave on moving
inputs. TFG also has a substantial tilt/bias penalty in the largest sea.
The [current results](moving-vessel-results.md) and paired comparison report
those limitations explicitly. Removing OU-III attenuation and changing TFG
startup/bias-release timing are not adopted here: neither has established
acceptable behavior across both stationary and moving records.

The input switch changes encounter frequency and steering together. The
physics and the stillness mechanism support encounter-frequency sensitivity,
but the paired default runs alone cannot isolate all causes. The reusable
counterfactual driver below can remove/add the declared track and stretch a
self-consistent motion history. Its exploratory bounds are diagnostics,
not estimator-accuracy guarantees. No old-clock counterfactual score is
used as current evidence.

## Reproduction

Build the repository with `make all`. Replay the production and interval-rate
profiles with `python3 tools/filter_sim_gates.py --run --jobs 3 --log-dir
/tmp/filter-gates --output reports/results/filter_sim_gates/measurements.json`.
Use `tools/stationary_vessel_control.py` with the separately checksum-verified
zero-speed v1.3.1 archive to reproduce the paired control arm.

For new counterfactual studies, unpack the v1.3.1 zero-speed and 5 kn archives
into `stationary/` and `moving/`, then run:

```sh
python3 tools/moving_vessel_diagnosis.py inputs --stationary stationary --moving moving
python3 tools/moving_vessel_diagnosis.py arms --stationary stationary --moving moving --out /tmp/motion-arms
python3 tools/moving_vessel_diagnosis.py run --tests tests --arms-dir /tmp/motion-arms --out /tmp/motion-scores.csv
python3 tools/moving_vessel_diagnosis.py table /tmp/motion-scores.csv
```

The driver resolves executables before changing the child working directory,
uses immediate stillness for attenuation diagnostics, and groups tables by
full record identity. No experimental TFG handoff hook is compiled into the
shipping source.
