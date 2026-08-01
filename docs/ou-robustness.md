# OU-III robustness study

`tools/ou_robustness.py` implements the Phase 2 sensitivity and degradation-case
study. It complements, but does not replace or retroactively tune, the paired
OU-II/OU-III experiment in `tools/ou_validation.py`.

## Design

The full profile uses the same ten predeclared wave-phase, IMU-noise, and
initialization seed triplets as the primary validation. Every comparison is
paired by all three seeds.

The sensitivity experiment uses the nominal noise-free OU-III operating point
for the $H_s=1.5$ m JONSWAP sea. It varies exactly one quantity at a time over
the multipliers `0.5, 0.75, 1.0, 1.25, 1.5`:

- OU correlation time $\tau$;
- stationary acceleration scale $\sigma_{aw}$;
- integral pseudo-measurement standard deviation $r_S$.

The implementation accepts $r_S$, with
$R_S=\operatorname{diag}(r_S^2)$, so the corresponding covariance multiplier is
the square of the reported $r_S$ multiplier. This is a one-factor-at-a-time
study; it does not estimate parameter interactions or claim global robustness.

The degradation cases are:

- a spectrally matched low-motion input scaled from $H_s=0.27$ m to $0.05$ m,
  testing the fixed sensor-noise floor;
- the existing $H_s=1.5$ m, $T_p=5.7$ s to $H_s=4.0$ m, $T_p=11.4$ s
  transition run over both a controlled 360-second ramp and a rapid 30-second
  ramp, with Adaptive and FixedNominal evaluated on the same realizations.

## Running

Fetch the versioned input bundle and run the CI-sized integration profile:

```bash
make fetch-sim-data
python3 tools/ou_robustness.py --mode smoke
```

Run the completed publication profile:

```bash
python3 tools/ou_robustness.py --mode full
```

Useful controls include `--sensitivity-scales`, the three independent seed
lists, `--duration-sec`, `--window-sec`, `--bootstrap-resamples`, and
`--output-dir`. Custom sensitivity lists must retain the publication anchors
`0.5`, `1.0`, and `1.5`; intermediate scales may be added or removed.
Smoke mode is an integration check and is not inferential evidence.

## Outputs

The versioned full result bundle under `reports/results/ou_robustness/` contains:

- raw simulator rows;
- group summaries with bootstrap 95% confidence intervals;
- paired differences and paired standardized effects;
- a JSON record containing the protocol, frozen nominal tuning point, raw rows,
  summaries, and effects;
- a self-hashing manifest with source/input hashes;
- generated LaTeX and two SVG publication figures.

Differences are always `left - right`. Positive error differences therefore
mean that the left-hand sensitivity or degradation case has higher error.
