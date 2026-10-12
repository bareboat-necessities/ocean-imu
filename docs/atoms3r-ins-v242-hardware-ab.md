# AtomS3R five-INS v2.4.2 A/B hardware qualification

This protocol answers one deployment question: does current magnetic
host-observation synchronization improve heading on real AtomS3R hardware
without starving valid magnetic service or degrading heave relative to v2.4.2?

## Builds

Use the same AtomS3R, calibration, USB/power arrangement and physical test
location. Test all five full-marine sketches: OU-II, OU-III, NLO, PII and TFG.

For current code use branch `validation/atoms3r-ins-v242-ab`. For the baseline
use `validation/v242-atoms3r-ins-ab`, which is v2.4.2 plus diagnostic-only
serial counters/rows. Neither validation branch changes estimator equations,
gates, gains or timing policy.

For capture set `SEA_STATE_SERIAL_NMEA=0` and `ARDUINO_PLOTTER=0`.
The existing 10 Hz debug stream then emits:

* `[INS_AB]`: time, filter, heading validity, heading, raw/clean heave and
  magnetic-use indicators.
* `[MAGTIME]`: host age/latency plus cumulative source/completion/service/
  rejection counts and maximum source/service gaps.

For NLO and PII, `service_total` counts continuous feedback offers, not
independent BMM150 observations. Use `source_total` and
`max_source_gap_ms` for source starvation. For OU-II, OU-III and TFG,
`service_total` is the count of actually applied discrete magnetic updates.

## Paired run order

Alternate firmware to reduce environment drift:

`v2.4.2 -> main -> v2.4.2 -> main`.

Do not recalibrate between a paired run unless the calibration is invalidated.
If recalibration is necessary, discard the pair and restart both sides from
the same newly saved calibration.

Run each sketch through these captures:

1. Still 60 s at four surveyed headings (0, 90, 180, 270 deg), one log per
   heading. Record the physical reference heading independently of the IMU.
2. Ten slow +/-90 deg yaw sweeps and ten rapid +/-180 deg yaw sweeps, returning
   to the same marked heading. Include combined roll/pitch during half.
3. Still 120 s after the dynamic sequence. This is the primary missed-service,
   north/south-flip and recovery capture.
4. Heave zero test: still 180 s on a fixed support.
5. Heave excursion test: repeat a mechanically measured vertical excursion
   (recommended 0.20-0.30 m peak) for at least 20 cycles, then return to the
   original support for 120 s.

Do not use hand-estimated displacement as an accuracy reference. A marked
vertical guide, linear scale, motion stage or independent displacement sensor
is sufficient.

## Quantitative acceptance

Analyze each fixed-heading or fixed-heave log with
`tools/hardware/atoms3r_ins_ab.py`. Example:

```sh
python3 tools/hardware/atoms3r_ins_ab.py \
  --run v242_ou3=ou3-v242-north.log \
  --run main_ou3=ou3-main-north.log \
  --heading-reference-deg 0 \
  --heave-reference-m 0 \
  --output ou3-north.csv
```

Current main passes only if all five sketches satisfy all of the following on
the repeated paired runs:

* no 180-degree heading flips or heading collapse to zero;
* heading validity is not lower than v2.4.2;
* fixed-heading RMS and peak error are improved or statistically unchanged,
  with no filter-specific regression hidden by averaging across filters;
* no unexplained multi-second magnetic source gap; compare
  `max_source_gap_ms` with the observed ~25 Hz host-source cadence;
* OU-II/OU-III/TFG do not lose valid discrete magnetic updates:
  `rejected_total`, `history_miss` and `max_service_gap_ms` must be
  inspected alongside source cadence;
* NLO/PII retain continuous magnetic feedback whenever a valid held source is
  available; source starvation is judged from `source_total` and source gap,
  not from the much higher continuous-offer count;
* zero-motion clean-heave drift and peak-to-peak noise are no worse than
  v2.4.2;
* measured heave-excursion amplitude/return-to-zero error are no worse than
  v2.4.2.

Do not call synchronization validated from replay alone. The BMM150 conversion
epoch remains physically unqualified; this protocol specifically measures the
real-device consequence of the deployed host-observation alignment.
