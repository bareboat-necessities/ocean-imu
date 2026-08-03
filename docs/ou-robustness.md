# OU-III robustness study

`tools/ou_robustness.py` implements the Phase 2 sensitivity and degradation-case
study. It complements, but does not replace or retroactively tune, the paired
OU-II/OU-III experiment in `tools/ou_validation.py`.

## Design

The full profile uses the same ten predeclared wave-realization, IMU-noise, and
initialization seed triplets as the primary validation. Every comparison is
paired by all three seeds.

The sensitivity experiment uses the nominal noise-free OU-III operating point
for the $H_s=1.5$ m JONSWAP sea and perturbs it over the multipliers
`0.5, 0.75, 1.0, 1.25, 1.5` along five directions.

One factor at a time, with the other two frozen:

- OU correlation time $\tau$;
- stationary acceleration scale $\sigma_{aw}$;
- integral pseudo-measurement standard deviation $r_S$.

Coupled, following the deployed regularization law
$r_S=\operatorname{clip}(1.2\,\sigma_{aw}\tau^3, 0.4, 35)$:

- `sigma_aw_rs`: $\sigma_{aw}\to c\,\sigma_{aw}$ with $r_S\to c\,r_S$;
- `tau_rs`: $\tau\to c\,\tau$ with $r_S\to c^3 r_S$.

The implementation accepts $r_S$, with
$R_S=\operatorname{diag}(r_S^2)$, so the corresponding covariance multiplier is
the square of the reported $r_S$ multiplier. All five directions are local
perturbations about one point; none estimates parameter interactions or claims
global robustness.

### Reading the frozen-companion columns

The $\tau$ and $\sigma_{aw}$ one-factor columns are nearly flat in displacement
error while $r_S$ and both coupled columns are not. That is a property of the
measurement geometry, not a parameter that fails to be applied:
$\sigma_{aw}^2$ does enter the discrete OU process covariance, and the raw rows
confirm the applied values match the requested ones. The accelerometer supplies
a direct 200 Hz observation of the latent $a_w$ with measurement noise orders of
magnitude below the stationary OU variance, so that gain is near saturation and
rescaling the prior barely changes the latent acceleration the kinematic chain
integrates; over one step $\Delta t/\tau \lesssim 10^{-2}$, so OU mean reversion
between updates is negligible too. Attitude *is* sensitive, because
$P_{a_w a_w}$ enters the accelerometer innovation covariance and therefore
governs how the residual is split between tilt error and wave acceleration --
halving $\sigma_{aw}$ roughly doubles pitch RMS.

So the supported statement is narrow: with $\tau$ and $r_S$ frozen at the
nominal point, displacement is locally insensitive to the *direct* OU
process-covariance effect of $\sigma_{aw}$, and attitude is not. The coupled
columns show the sensitivity the online tuner actually experiences.

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
python3 tools/ou_robustness.py --mode smoke \
  --sensitivity-scales 0.5,1.0,1.5
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

The generated LaTeX also defines `\OURobustnessTauSpan`,
`\OURobustnessSigmaSpan`, `\OURobustnessRSSpan`,
`\OURobustnessCoupledSigmaSpan`, and `\OURobustnessCoupledTauSpan`: the
peak-to-peak mean vertical error across each sweep, in points of $H_s$.

Differences are always `left - right`. Positive error differences therefore
mean that the left-hand sensitivity or degradation case has higher error.
