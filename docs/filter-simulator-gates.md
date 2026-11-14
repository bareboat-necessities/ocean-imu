# Filter simulator regression gates

External-data regression now uses oceanography-waves-lib v1.3.0,
`sim-data-files-vessel-rao-28ft-5kn.zip`: a 28 ft estimated hull response at
5 kn with heading oscillating ±3° about north every 120 s. The dataset change
has its own measured baseline. Estimator equations, tuning, default sensor
noise and seed schedules are unchanged.

See [current measurements](moving-vessel-results.md) and the
[track/scoring convention](moving-vessel-replays.md). All physical sensor
inputs and turning acceleration remain intact. Only reference displacement
and velocity are expressed about the declared cruise track. In particular,
3D error percentages are normalized by wave displacement, not total travel.

`tests/common/CurrentReplayLimits.h` contains separately qualified production
and interval-rate limits. Each channel uses the worst of eight full records,
plus 0.5%, rounded upward to 0.00001. Stationary-data limits are retained in
the baseline comparison as historical measurements; they do not cap limits
for a different physical experiment. The report exposes degradations rather
than attributing a dataset change to estimator improvement.

W3d scores use the trailing 900 s of each 1200 s, 200 Hz record. Heave is
RMS/incident Hs; 3D displacement is RMS/max reference wave-displacement norm.
Accelerometer Z bias uses max absolute true Z bias; each vector-bias score
uses max true vector-bias norm. Absolute displacement and bias RMS accompany
these different percentage normalizations. AHRS uses the trailing 60 s.
NLO gates raw heave, including its mean error. Adaptive PII has separately
measured vertical/yaw bars with the same 0.5% margin, rounded upward to
0.01 in its quoted units; it rejects nonfinite heave and attitude scores. Heave-only gates remain
separate for Mahony, proxy and ideal-attitude frontends and each method.

The inventory covers all native W3d filters, the legacy cubic TFG diagnostic,
AHRS, heave-only methods, NLO, the independently checked CARMA process prior,
and the lever-arm study instrument. Detrender inputs also use wave-relative
horizontal displacement from the moving records; cruise travel is not added
to their wave-reference target.

Nonfinite scored channels, missing/duplicate records, incomplete histories
and wrong input hashes remain failures. Full-rate samples determine RMS;
plot decimation never changes scoring. The conditional eight-record gates
are deterministic sentinels, not population confidence intervals, general
accuracy targets or guarantees for engine noise and alternate priors.

After building, reproduce the passing evidence with:

```sh
make all
python3 tools/filter_sim_gates.py --run --jobs 3 \
  --log-dir /tmp/filter-sim-gates \
  --output reports/results/filter_sim_gates/measurements.json
python3 tools/moving_vessel_publication.py
```

`--measure` permits collecting a changed experiment before gate qualification.
It retains failure exit statuses and still requires complete finite scores;
it neither changes limits nor produces a passing publication report. The
ordinary command requires every gate to pass. Both record input and binary
hashes and the source revision. Dense outputs and runner logs are omitted.
