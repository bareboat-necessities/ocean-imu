# Bias-prior and reporting-filter feasibility

No candidate is promoted. The default numerical gates remain OU-II 6/8,
OU-III 6/8 and NLO 7/8 on the pinned v1.2.1 vessel records. The gate limits,
1200 s input histories, final 900 s scoring windows and injected sensor noise
are unchanged. Seed overrides are paired between baseline and candidate.

`studies.json` retains every aggregate row from 29 studies, including failed
records and candidates. NLO rows contain eight records each; OU rows contain
one. These are engineering parameter studies, not source-uniform stability
evidence. The source is commit `a2fbc8f86d02d1560d50725b16a41bd8e267d50c`
plus experimental hooks archived in `screening-source.patch`.

## Findings

Reducing horizontal accelerometer-bias driving noise improves OU-II's
high-sea roll/bias errors, but joint prior and magnetic-covariance changes
introduce low/medium-wave heave failures. Adjusting the pseudo-channel ratio
between .28 and .4 moves that failure between records. Small wave-time-
constant changes also fail the complete default set. OU-III can pass all
eight default records with a longer bias prior and reduced horizontal bias
noise, but that improvement does not generalize to separate sensor and
initialization draws.

| Paired validation | Default failing records, baseline → candidate | Fresh failing records | Fresh gate violations |
|---|---:|---:|---:|
| `ouII-final` | 2 → 1 | 26 → 28 | 81 → 78 |
| `ouIII-final` | 2 → 0 | 28 → 30 | 70 → 82 |
| `nlo-final` | 1 → 0 | 7 → 4 | 7 → 4 |

These final studies use fresh seeds 12713, 13901, 14923 and 15791, after
earlier separate checks using 4111/5237/6263/7307 and, for OU-II,
8191/9341/10427/11617. They are not repeatedly adjusted holdout draws.
Default records are screening data throughout.

The final OU-II candidate uses bias driving standard deviations
(.00015, .00015, .0004), bias time constant 20000 s, magnetic sigma scale 2
and pseudo-channel ratio .36. The OU-III candidate uses
(.0003, .0003, .0004), 50000 s, magnetic sigma scale 1.5 and sigma
coefficient .7. Both remain rejected because of the complete validation
tradeoffs; improving a fixed default gate is insufficient.

The NLO experiment adds an optional causal reporting high-pass horizon
`min(50 s, N/f_wave)`, using its independent reference-heave frequency
tracker. With theta gain .52 and N=4, all defaults pass, but the worst
fresh vertical error rises from 8.24759% to 12.5197% of incident Hs and
the fresh mean rises from 7.04781% to 7.07271%. A .54 gain also passes
defaults but has mixed earlier fresh-seed results. The optional reporting
change and both new settings are rejected rather than shipped.

## Reproduction

Use a separate checkout of the source commit above, apply
`screening-source.patch`, fetch the pinned data with `make ensure-sim-data`,
and build the OU-II, OU-III and NLO simulators using Eigen 3.4.0. The patch
contains the final experimental hooks, including scalar/per-axis bias
controls and the optional wave-relative reporting horizon; unset hooks
preserve the original default paths. It is not a proposed production patch.

For a named study, extract its unique `(config, env)` pairs from `runs`
into the runners' `{name, env}` configuration JSON. Run
`tools/rao_parameter_tuning.py` with the recorded family, record subset and
seed list, or `tools/observer_parameter_tuning.py --family nlo` with the
recorded seed list. Both accept `--configs`, `--jobs` and `--output-dir`.
Each archived manifest retains input provenance, binary/source hashes and
the execution arguments. Earlier screens used earlier subsets of these
hooks; their original hashes are preserved, not restamped as final-source
validation. The final paired studies used the complete archived patch.
