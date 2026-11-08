# Shared-mode CARMA method-article validation

This bundle supports `doc/kalman_carma_shared_modes/kalman_carma_shared_modes.tex`.
It evaluates the current complete estimators with the pinned
`oceanography-waves-lib` v1.2.3 28 ft vessel-response dataset. It is synthetic
development evidence for one hull model, not sea-trial accuracy.

## Experiment matrix

| Study | Methods | Motion / sensor draws | Replays |
|---|---|---|---:|
| Stationary comparison | Shared CARMA, uncoupled CARMA, OU-II, OU-III, TFG | 8 seas × 10 paired seed triplets | 400 |
| Common-output observer comparison | PII, TVG-NLO | 8 seas × first triplet | 16 |
| Low–high–low transition | Shared CARMA, uncoupled CARMA, OU-III | 10 paired triplets | 30 |
| Raw engine-vibration controls | Shared CARMA, uncoupled CARMA, OU-III | 14 conditions × 8 seas × first triplet | 336 |
| Guard-only and guard-plus-covariance arms | Shared CARMA, OU-III | 8 conditions × 8 seas × first triplet × 2 arms | 256 |
| Total | | | 1038 |

The stationary draws retain the common 0.02–1.60 Hz phase-surrogate procedure
of `tools/ou_validation.py`. The primary JONSWAP and separate PM-Stokes
summaries average the four height-normalized errors within each seed before
averaging or bootstrapping seeds. PII and TVG-NLO are single-draw comparisons;
their unavailable horizontal and bias outputs, and TVG-NLO yaw, remain unscored.
PII retains its 20 Hz magnetic path; the five Kalman methods use 25 Hz, and
TVG-NLO has no compass or horizontal GNSS aiding.

The uncoupled arm uses the same 27-coordinate cascade with coupling off. It
realizes the independent CARMA law and isolates sharing with a matched frontend;
it is not a cost measurement of a separate 21-state implementation.

All ordinary replay scores use the final 900 s (180,000 samples) of a 1200 s,
200 Hz run. Transition scores additionally retain seven full-rate intervals:
300–400, 400–520, 520–640, 640–800, 800–920, 920–1040 and 1040–1200 s. Both
crossfades are kinematically closed quintic blends, and no estimator is reset.

The engine conditions are off; 800/1600/2400/3200 rpm; 0.30/0.60/1.20 m/s²
hull vibration; 20/40/80/160 Hz sensor bandwidth; the same bandwidths at matched
recorded power; and the nominal condition with modeled gyro perturbation off.
The nominal 2400 rpm, 0.60 m/s², 80 Hz condition is reused across the controls.
The raw arm disables the guard and covariance inflation. Mitigation uses the
existing 14 Hz two-pole guard and covariance gain 0.75, without retuning.
Engine summaries pool eight per-record RMS errors quadratically, not arithmetically.

## Retained evidence

- `protocol.json`: input checksums, source revision/hashes, binary hashes,
  seed triplets, sample rate and transition intervals.
- `scores.csv`: one row per replay; scientific statistics read this file.
  `source_record_hs_m` retains the released source height. For the round trip,
  `hs_m=4` denotes the high endpoint after scaling the 8.5 m source, not the
  instantaneous height during the crossfade.
- `engine_summary.csv`: pooled engine cells derived from the full-rate scores.
- `publication_summary.json`: paired effects, publication-source hashes and
  rendered SVG sizes.
- `traces/*.csv.gz`: deterministic, extrema-preserving drawing traces capped
  at 2000 points each. They never enter a score, fit or bootstrap.

The original full-rate replay CSVs are temporary and are discarded. The shared
`plots/vector_export.py` exporter bounds SVG and PGF drawing geometry and shares
the point budget among overlaid curves. The raw signal and scoring resolution
are independent of the drawing budget.

## Reproduce

Build the simulators with the repository's ordinary `make all` environment,
then fetch the checksum-verified dataset with `tools/sim_dataset.py`. Run:

```sh
python3 tools/carma_shared_modes_validation.py run --data-dir tests/kalman_carma_shared_modes --jobs 4 --seeds 10 --results /tmp/carma-shared-modes-replay
python3 tools/carma_shared_modes_validation.py figures --results /tmp/carma-shared-modes-replay
```

The run stage resumes only against an identical protocol. Its append-only
`runs.jsonl` checkpoint is ignored by Git; `scores.csv` is the retained per-run
evidence. To start a changed protocol, use another `--results` directory.

The figures stage rebuilds all tables, SVG and PGF charts from the committed
scores and compact traces without replaying simulations. It is also called by
`plots/kalman_carma_shared_modes/draw_plots.sh`, so ordinary article CI rebuilds
the charts. Compile the main `.tex` with LuaLaTeX after figure generation.
To regenerate only the committed article figures, omit `--results` from the
figures command. Fresh simulations use a new output directory because the
source revision and binary checksums are intentionally guarded on resume.

Focused checks:

```sh
cd tests/validation
python3 -m unittest -v test_carma_shared_modes_validation test_vector_export test_figure_decimation test_ou_roundtrip_transition
```
