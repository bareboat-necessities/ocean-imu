# Conventional heave baselines: HPDI and FDDI

Two conventional comparators for OU–III's vertical channel. Both consume the
**same measurement-only inputs** OU–III's front end produces: the Mahony-proxy
vertical acceleration `a_up` (paper Eq. 73) and, for period-scaled cutoffs, the
canonical zero-crossing period `Tz_hat` (Eq. 80). Reference motion is used only
to score results and to select frozen parameters on a development split.

| Method | Causal | Where | Transfer from acceleration to heave |
|---|---|---|---|
| HPDI-classic | yes | `src/hpdi/HeaveHPDI.h` | `s^(n-2) / D(s)`, i.e. Butterworth HP of order n then 1/s² (m = n) |
| HPDI-compensated | yes | same, m < n | `[1 − R(s)/D(s)] / s²` |
| FDDI | no (offline) | `tools/heave_baselines.py` | `−W(f) / (2πf)²`, raised-cosine high-pass taper |

## Design summary

`D(s)` is the order-n Butterworth polynomial at cutoff ω_c, and `R(s)` is its
lowest m terms. The displacement gain relative to ideal double integration is
`G = 1 − R/D`.

- **Classic (m = n).** `G = s^n/D`. The error `|G − 1|` decays only as
  `c_{n−1} ω_c/ω`. This is the persistent phase lead of a causal high-pass, and
  is why textbook heave filters distort the wave band.
- **Compensated (m < n).** `|G − 1| ~ c_{m−1} (ω_c/ω)^(n−m+1)`. The cost is gain
  peaking near the corner (|G| up to ~2 for n = 4, m = 3; ~6 for n = 6, m = 3).
- **m ≥ 3.** Zero steady-state heave for a constant accelerometer bias, and
  bounded variance under a random-walk bias.
- **m = 2.** Leaves an offset of `c_2 b / ω_c²`. The member n = 3, m = 2 with
  ω_c = ω_R reproduces the reduced SpectralMSE response, paper Eq. (86), so the
  HPDI family contains OU–III's own scalar reduction.

**Discretization.** The trapezoidal rule (bilinear transform) is applied in delta
form on an observer-canonical realization with frequency-normalized states. This
stays float-accurate at ω_c·h ≈ 10⁻³. `selfcheck` verifies agreement with
`scipy.signal.bilinear` of H(s) to about 1e-8.

**Embedded footprint.** `HeaveHPDI<float, 6>`:
- 472 B of state, no heap.
- About 3 KB of code at `-Os` with `-fno-exceptions -fno-rtti`.
- Per sample: n² + 2n MACs plus one `expm1f`.
- Coefficient refresh about every 0.1 s (sample-and-hold, like OU–III).

## Build and test

```sh
mkdir -p build
g++ -std=c++14 -O2 -Wall -Wextra -I src tests/hpdi/test_heave_hpdi.cpp -o build/test_heave_hpdi
./build/test_heave_hpdi            # 67 checks
g++ -std=c++14 -O2 -I src tools/hpdi_replay.cpp -o build/hpdi_replay
python3 tools/heave_baselines.py selfcheck --hpdi-bin build/hpdi_replay
```

The header compiles cleanly from C++11 through C++20 under `-Wconversion
-Wdouble-promotion`.

## Harness export (one change in the replay harness)

For each scored record, write a CSV with these columns:

| Column | Content |
|---|---|
| `t` | time, s |
| `a_up` | Mahony-proxy up-positive vertical acceleration, gravity removed (Eq. 73), after the same vibration prefilter OU–III sees |
| `Tz_hat` | canonical period from the measurement-only estimator (Eq. 80); NaN before its startup gate |
| `z_ref` | reference heave; pass `--flip-z` if it is NED-down |

Use the same seed triplets as the OU–III runs so the comparison stays paired.

Then write a manifest listing every record:

```csv
path,family,hs,split,seed
rec/jonswap_hs0.27_s11.csv,JONSWAP,0.27,eval,11
rec/jonswap_hs0.27_d1.csv,JONSWAP,0.27,dev,d1
```

`split=dev` holds the additional draws used for comparator selection, as in
Table VIII. `split=eval` holds the ten paired seed triplets.

## Run

```sh
python3 tools/heave_baselines.py tune --manifest runs/manifest.csv --out runs/hpdi_fddi_frozen.json
python3 tools/heave_baselines.py evaluate --manifest runs/manifest.csv \
    --params runs/hpdi_fddi_frozen.json --out runs/hpdi_fddi_results.csv \
    --compare-csv runs/ou3_results.csv --compare-method OU-III
```

**What `tune` does.** It grid-searches both baselines on the dev split and
freezes three configurations:
- FDDI: fixed or period-scaled cutoff, crossed with the taper ratio.
- HPDI-classic: the best configuration with m = n.
- HPDI-compensated: the best configuration with m < n.

It warns when an optimum lands on a grid edge.

**What `evaluate` does.** It scores the frozen configurations on the eval split
and writes long-format per-record rows (method, family, hs, seed, pct_hs). With
`--compare-csv` it also reports the seed-level paired contrast:
- percentile bootstrap, 10,000 resamples, PCG64 seed 20260317;
- exact sign-flip p-value.

This matches the primary-endpoint protocol. The compare CSV needs the columns
`method,family,hs,seed,pct_hs`.

Add `--float32` to replay HPDI in single precision, as deployed on the ESP32-S3.

## Reporting notes

- FDDI is zero-phase and uses future samples. Report it as an offline reference
  bound, not a real-time competitor.
- Both baselines estimate heave only. Compare them on the vertical endpoint, and
  keep OU–III's 3-D and attitude outputs as a separate claim.
- `synth` writes idealized vertical-only records: no RAO, no tilt, no proxy, and
  the true T_z. Use them for plumbing tests only; they are not evidence.
