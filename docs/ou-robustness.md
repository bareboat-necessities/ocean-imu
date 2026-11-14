# OU-III robustness study

The full study uses the checksum-pinned v1.3.0 28 ft / 5 kn moving boat,
with ±3° heading over 120 s. See the [current measured baseline](moving-vessel-results.md)
and [reference-coordinate convention](moving-vessel-replays.md). The 310 new
simulator rows retain ten paired seed triplets, five sensitivity multipliers,
and all low-motion and transition conditions.

`tools/ou_robustness.py` implements the sensitivity and degradation-case study
that complements the paired OU-II/OU-III comparison in `tools/ou_validation.py`.
It does not retroactively tune the primary comparison.

The current committed robustness bundle was regenerated from source commit
`be5e5acd8cc599693a322fba9a7d2c17e3a85f7a`. It was produced outside GitHub
Actions, so `workflow_run` is null; the replay hashes and recorded build
environment carry the provenance.

## Design

The full profile uses the same ten predeclared wave-realization, IMU-noise, and
initialization seed triplets as the primary validation. Every comparison is
paired by all three seeds.

The local sensitivity experiment starts at the nominal noise-free OU-III
operating point for the $H_s=1.5$ m JONSWAP sea and perturbs the requested
parameters over the committed multiplier grid. It separates direct changes to
OU correlation time $\tau$, stationary acceleration scale $\sigma_{aw}$, and
integral pseudo-measurement scale $r_S$ from coupled perturbations that propagate
through the deployed scheduler. These are local sensitivity directions; they do
not estimate a global optimum.

The degradation cases include:

- a spectrally matched low-motion input scaled below the main study envelope;
- the controlled $H_s=1.5$ m, $T_p=5.7$ s to $H_s=4.0$ m, $T_p=11.4$ s
  transition under both the standard 120 s blend and a rapid 30 s blend. The
  standard case is the crossfade the paired validation study scores, so both
  studies mean the same record by "the controlled transition"; it was a 360 s
  blend until that crossfade was shortened;
- Adaptive and FixedNominal modes evaluated on identical paired realizations.

The study scores continuous displacement, attitude, bias, and tuning metrics.
As with the primary validation, deterministic simulator regression thresholds
are **not** statistical acceptance criteria and are not exported into the
committed CSV/JSON evidence as pass/fail fields.

## Current degradation results

- Reducing incident height from 0.27 m to 0.05 m increases mean normalized heave RMS from **8.33% to 25.86% Hs**. Absolute RMS is 0.0129 m in the smaller sea.
- Halving the coupled time scale and regularizer raises nominal-sea heave RMS from **6.94% to 13.92% Hs**. Halving only the time scale gives 7.02%; the coupling matters.
- Compressing the crossfade from 120 s to 30 s raises Adaptive RMS from **6.07% to 6.23% Hs**. The paired difference is +0.16 points, with bootstrap 95% interval **[+0.04, +0.29]**. Adaptive remains below FixedNominal in the rapid case by 3.84 points.

These are paired results within this development ensemble. The shorter crossfade also changes the scored mixture of sea conditions.

## Provenance

The robustness manifest uses the same two-layer contract as the primary study.
Immutable `replay_provenance` pins the full OU-III simulator/filter dependency
closure, the simulator build file, replay inputs, and the normalized raw-row
SHA-256. A later `restatement` block may identify different Python analysis code
and derived outputs, but it cannot replace the replay source commit or replay
hashes.

Before `--restat-from` writes anything, the current replay closure is compared
with the recorded one. A change in `SeaStateFusionFilter_OU_III.h`,
`Kalman3D_Wave_OU_III.h`, the simulator, shared replay code, or the recorded
Makefile is therefore a hard replay requirement rather than a warning. Analysis
or presentation changes may be restated only while the immutable replay and raw
rows still verify.

`tools/ou_evidence_contract.py --check` works both in a Git checkout and in a
source archive without `.git`; archive mode uses the committed manifest and
file hashes and does not weaken dependency verification. The protocol separately
records that all completed machine-readable replays are included and that
deterministic simulator regression gates are not Monte Carlo acceptance
criteria.

## Running

Fetch the versioned inputs and run the short integration profile:

```bash
make fetch-sim-data
python3 tools/ou_robustness.py --mode smoke \
  --sensitivity-scales 0.5,1.0,1.5
```

Run the completed publication profile:

```bash
python3 tools/ou_robustness.py --mode full
```

Useful controls include `--sensitivity-scales`, the three independent seed
lists, `--duration-sec`, `--window-sec`, `--bootstrap-resamples`, and
`--output-dir`. Smoke mode is an integration check and is not inferential
evidence.

The production CI path is `.github/workflows/ou-validation.yml`, which shards
validation and robustness separately, combines the full rows, verifies the
publication/evidence contract, and commits both bundles together.

### Restating a committed bundle

A pure change to statistical presentation can reuse the archived rows:

```bash
python3 tools/ou_robustness.py \
  --restat-from reports/results/ou_robustness/ou_robustness.json \
  --bootstrap-resamples 10000 \
  --output-dir reports/results/ou_robustness
```

Restatement cannot substitute for a replay after any file in the recorded
implementation closure changes. The evidence contract reports that mismatch and
requires a full regeneration. Source movement is therefore treated according to
what moved: formatting/statistical code can be restated when appropriate;
estimator/simulator implementation changes require new replays.

## Outputs

The versioned full result bundle under `reports/results/ou_robustness/` contains:

- raw completed-replay rows without legacy simulator gate fields;
- group summaries with bootstrap 95% confidence intervals;
- paired differences and paired standardized effects;
- a JSON record containing the protocol, frozen nominal tuning point, raw rows,
  summaries, and effects;
- a self-hashing manifest covering inputs, results, implementation dependencies,
  and analysis-pipeline sources;
- generated LaTeX and two SVG publication figures.

Differences are always `left - right`. Positive error differences therefore mean
that the left-hand sensitivity/degradation case has higher error.
