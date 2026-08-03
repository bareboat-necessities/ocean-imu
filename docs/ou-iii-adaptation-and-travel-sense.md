# OU-III adaptation loss and wave travel-sense reporting

Investigation record and remediation plan for the two reviewer findings on the
OU-III manuscript:

1. OU-III improves vertical RMS but degrades 3D RMS in moderate and large seas,
   and adapting `tau`/`sigma_aw` is worse than adapting `r_S` alone.
2. The reported travel-sense classes are defined against the estimator's own
   axis representative rather than against physical propagation.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_iii/kalman_ou_iii-sim`, scoring the trailing 900 s.

## 1. Finding 1: where the adaptation actually loses

### 1.1 The leading hypothesis was posterior/process covariance inconsistency

It is not the cause. Two independent checks:

- The committed study already prices the policy. `Adaptive` versus
  `AdaptiveHeldCovariance` (identical tuning, periodic re-alignment of the
  posterior `a_w` marginal switched off) differ by 0.002 m in 3D RMS at
  `H_s = 8.5` m — 0.07% of a 3.10 m error.
- A direct probe (40 000 steps, `sigma_aw` swept over its adaptation range,
  marginal re-aligned every 0.1 s) shows the minimum eigenvalue of `P` staying
  at float round-off both before and after the re-alignment
  (`-6.8e-05` before, `-2.3e-05` after, against diagonal entries of order 1).
  The re-alignment *inflates* the marginal, typically by 40-100x, and adding a
  PSD term to a diagonal block cannot destroy PSD.

The re-alignment does discard a large amount of posterior information: the
inflation factor is not "bounded" in any useful sense.

Making the operation *consistent* makes it worse, which settles the question
from the other side. `W3D_AW_COV_SYNC=congruent` re-aligns the marginal through
a congruence `x_aw -> A x_aw`, `A = L_new L_old^-1`, which reaches the same
marginal, leaves the whitened cross-covariance untouched and stays PSD by
construction. On `jonswap H8.500`:

| a_w sync policy | disp Z | 3D |
| --- | --- | --- |
| `periodic` (deployed overwrite) | 0.645 | 2.767 |
| `congruent` | 0.678 | 3.618 |
| `congruent`, heading rotated 45 deg | 3.655 | 11.849 |

The congruence propagates the marginal inflation into the cross-covariances,
which is what consistency demands and which the filter cannot absorb. The
deployed overwrite is not better because it is more correct; it is better
because it silently damps an intervention that is too large. The conclusion is
that the periodic re-alignment should be retired or bounded, not made
self-consistent - and the `HeldCovariance` ablation already shows retiring it
costs nothing.

### 1.2 The frequency that drives `tau` does not track the sea state

`tau_target = tau_coeff * 0.5 / f`, with `f` from the frequency tracker running
on the low-passed body-Z acceleration proxy. Acceleration spectra are
`(2*pi*f)^4` weighted, so their apparent frequency sits far above the elevation
peak, and the offset grows with the sea state:

| record | `f_p` (elevation) | `f_mean` (acceleration) | ratio |
| --- | --- | --- | --- |
| jonswap H0.270 L14.047 | 0.335 Hz | 0.548 Hz | 1.63 |
| jonswap H1.500 L50.710 | 0.177 Hz | 0.474 Hz | 2.68 |
| jonswap H4.000 L112.766 | 0.118 Hz | 0.440 Hz | 3.72 |
| jonswap H8.500 L202.839 | 0.085 Hz | 0.417 Hz | 4.90 |

The acceleration-band frequency is nearly constant across a 4x change in sea
state. The channel the filter uses to learn "how long the waves are" carries
almost no sea-state information — only estimation noise. Across the ten seed
triplets of the committed study at `H_s = 8.5` m, `tau_applied` is
`1.660 +/- 0.431` s against the `2.330` s the same tuner produces on the
noise-free record, and `r_S = clip(1.2 * sigma_aw * tau^3)` cubes that error:
`11.5 +/- 4.3` against `34.7`.

### 1.3 The integral regularizer is therefore far too strong in developed seas

The `S = 0` pseudo-measurement is a high-pass on displacement whose corner is
set by `r_S`. With `r_S` too small the corner sits inside the wave band and the
filter removes the signal it is supposed to estimate. Scoring
`jonswap H8.500` (900 s window, fixed `tau = 2.33`, `sigma_a = 2.283`):

| `r_S` (m*s) | `R_S_xy` | disp X | disp Y | disp Z | 3D |
| --- | --- | --- | --- | --- | --- |
| 3 | 0.36 | 3.478 | 2.170 | 1.489 | 4.362 |
| 8 | 0.36 | 2.857 | 1.834 | 0.841 | 3.498 |
| 15 | 0.36 | 2.221 | 1.476 | 0.590 | 2.732 |
| 35 | 0.36 | 1.475 | 1.080 | 0.405 | 1.873 |
| 3 | 1.00 | 2.835 | 1.798 | 1.343 | 3.616 |
| 8 | 1.00 | 1.848 | 1.266 | 0.782 | 2.372 |
| 15 | 1.00 | 1.353 | 1.023 | 0.568 | 1.789 |
| 35 | 1.00 | 0.964 | 0.902 | 0.402 | 1.380 |

The response is monotone up to `MAX_R_S = 35`, which is where `FixedOracle`
sits (34.66) and is exactly why `FixedOracle` wins: the clamp binds at the
operating point the filter needs.

### 1.4 The horizontal axes are pushed further into the same failure

`R_S_xy_factor = 0.36` scales the horizontal pseudo-measurement standard
deviation down, i.e. it makes the horizontal high-pass 2.8x *stronger* than the
vertical one, in the channel where the corner is already too high. Adaptive
mode, sweeping only that factor:

| record | 3D at `R_S_xy = 0.36` | 3D at `R_S_xy = 1.00` |
| --- | --- | --- |
| jonswap H0.270 | 0.068 | 0.107 |
| jonswap H1.500 | 0.289 | 0.291 |
| jonswap H4.000 | 1.188 | 0.891 |
| jonswap H8.500 | 2.767 | 1.829 |

0.36 is a small-sea optimum applied to every sea state. It costs 25% of the 3D
RMS at `H_s = 4` m and 34% at `H_s = 8.5` m.

By contrast `S_factor = 1.87` (the horizontal stationary acceleration scale) is
mis-specified relative to physics — the records give per-axis
`sigma_ax/sigma_az = 0.81` and `sigma_ay/sigma_az = 0.55`, with
`sqrt(sigma_ax^2 + sigma_ay^2)/sigma_az = 0.99` as deep-water theory requires —
but it is *not* what drives the loss: sweeping it 1.87 -> 0.80 -> 0.60 moves
3D RMS by under 2%.

### 1.5 Mechanism, stated once

The tuner learns the sea state through a frequency estimate taken in the wrong
spectral band. `tau` therefore stays roughly constant while the true wave period
grows, `r_S ~ sigma * tau^3` stays far below what the sea requires, the
`S = 0` regularizer's corner ends up inside the wave band, and the horizontal
axes get an extra 2.8x of the same medicine. The vertical channel still beats
the fixed reference because `r_S` at least moves with `sigma`; the horizontal
channel does not, so 3D RMS degrades exactly where the paper reports it.
Freezing `tau`/`sigma_aw` at the nominal point (`AdaptiveRSOnly`) helps because
it stops the OU model from following a frequency signal that carries no
information, while `r_S` keeps tracking `sigma`.

## 2. Finding 2: the travel-sense classes are a gauge label

### 2.1 Controlled experiment

`tools/rotate_record_heading.py` rewrites a record so the vessel heading is
rotated by a fixed angle about the world vertical: the wave field, the world
specific force and the seeds are untouched, the attitude trajectory becomes
`R_new(t) = Rz(psi) R_old(t)`, body specific force is re-resolved in the new
body frame, and body angular rate is left alone (it is invariant under a
constant world-frame rotation of the attitude). Running the same record at
`psi = 0` and `psi = 180` degrees:

| quantity | heading 0 deg | heading 180 deg |
| --- | --- | --- |
| disp Z RMS | 0.645 m | 0.654 m |
| disp 3D RMS | 2.767 m | 2.706 m |
| axis error vs generator | 4.70 deg | 5.19 deg |
| FORWARD share | 95.5% | 1.8% |
| REVERSE share | 1.8% | 94.6% |

Nothing physical changed and the axis estimate is unaffected, but the reported
class inverts completely. `KalmanWaveDirection` seeds its representative from
boat `+X`, so `FORWARD` means "along whichever end of the axis currently points
forward". The class is a gauge label.

### 2.2 The directed output, however, is correct

The estimator's directed propagation vector, with the vessel heading removed,
scored against the record's physical propagation direction (`A + 180`) by the
new `dir_travel_*` metrics:

| heading | travel error | correct | wrong | unresolved | FORWARD | REVERSE |
| --- | --- | --- | --- | --- | --- | --- |
| 0 deg | +4.56 deg | 95.5% | 1.1% | 3.4% | 95.5% | 1.1% |
| 45 deg | +5.43 deg | 95.5% | 1.2% | 3.3% | 95.5% | 1.2% |
| 180 deg | +4.91 deg | 94.6% | 1.8% | 3.6% | 1.8% | 94.6% |

The correctness rate is flat across headings while the class shares invert. The
propagation sense the estimator resolves is physical and heading-invariant; only
the exported enumeration and the metric built from it were not. The paper
therefore understates the result: 95% is a correctness rate, not just a
commitment rate.

Scoring the 45-degree record also exposed a latent bug in the *axis* metric: it
compared a boat-frame angle against a world-frame azimuth, so it reported a
39.5-degree axis error for an axis that was correct to 5.5 degrees. Every
shipped record has heading 0, which is why it never showed. The vessel heading
is now removed before the comparison.

### 2.3 Record convention

The generator azimuth `A` in the record filename is the direction the waves come
**from**: the propagation-to vector recovered from the truth channels (principal
axis of horizontal displacement, sense from
`sign(<xi_parallel * d(zeta)/dt>) < 0`) lies at `A + 180` in every record. This
was never asserted anywhere, which is why the 180-degree question could not be
settled from the existing metrics.

## 3. Remediation plan

Stages are ordered so each one is separately reviewable and separately
falsifiable. Stages A and B change no filter behaviour.

### Stage A - make the failure visible (no filter behaviour change)

- A1. `tools/rotate_record_heading.py`: the heading-rotation transform used in
  section 2.1, with the angular-rate invariance documented. (done)
- A2. Truth-referenced travel-sense metrics in the simulator harness:
  `dir_travel_error_deg`, `dir_travel_rmse_deg`, `dir_travel_correct_pct`,
  `dir_travel_wrong_pct`, `dir_travel_unresolved_pct`, scored against `A + 180`
  after removing the reference heading. `dir_sense_dominant_pct` is kept so the
  existing series stays readable. (done)
- A3. Assert the record convention (section 2.3) and the invariants of the
  heading transform in `tests/validation/test_record_conventions.py`. (done)
- A4. Fix the axis metric's own heading dependence: `dir_deg_hist` compared a
  boat-frame angle against a world-frame azimuth, which is only valid because
  every shipped record has heading 0. Scoring the 45-degree record gave a
  39.5-degree "error" for an axis that was in fact correct to 5.5 degrees. The
  vessel heading is now removed before the comparison. (done)
- A5. Still to do: adaptation diagnostics behind `W3D_WRITE_ADAPT_LOG` -
  `tau`, `sigma_aw`, `r_S`, accelerometer NIS, `trace(P)`, `trace(P_aw)` and the
  implied `S`-loop corner frequency - plus an offline splitter that separates
  displacement error into drift, wave-band and high-band parts.

### Stage B - price the covariance-synchronisation policy (done)

- B1. `synchronize_aw_covariance_to_stationary_congruent()` implements the
  consistent congruence re-alignment alongside the deployed overwrite, selected
  by `W3D_AW_COV_SYNC=congruent`. The deployed default is unchanged.
- B2. Unit test in `tests/kalman_ou_iii/kalman_ou_common-test.cpp`: the target
  marginal is reached exactly, the whitened cross-covariance is unchanged, the
  non-`a_w` marginals are untouched, and the joint covariance stays PSD, for
  marginal changes of both 16x and 1/16.
- B3. Measured (section 1.1): the consistent version is worse, so the remaining
  work is to retire or bound the periodic re-alignment rather than to fix its
  algebra. That decision belongs with the Stage E re-run, because
  `AdaptiveHeldCovariance` is already the mode that measures it.

### Stage C - estimate the wave period in the wave band

- C1. New estimator: band-limited double integration of the vertical
  acceleration proxy (fixed 0.02 Hz high-pass, well below every sea state of
  interest) into velocity and elevation proxies, long-horizon EWMA variances,
  and `omega_z = sigma_v / sigma_eta`, i.e. the standard zero-crossing period
  `T_z = 2*pi/omega_z`. No dependence on the filter's own states, so no new
  feedback loop.
- C2. `tau_target = tau_coeff * T_z / 2`, restoring the documented intent
  ("tau ~ half the dominant period"). Re-fit `tau_coeff` and raise
  `MAX_TAU_S` (currently 3.0 s, binding for `T_z ~ 8.5` s).
- C3. Re-fit `R_S_coeff` and the `r_S` clamps against the re-based `tau`.
  The `r_S ~ sigma * tau^3` form is physically right (`sigma_S ~ sigma_p/omega`
  and `sigma_p ~ sigma_a/omega^2`); only its frequency input and its clamps are
  wrong. Section 1.3 shows `MAX_R_S = 35` binding at `H_s = 8.5` m.
- C4. Validate `T_z` against the record truth (`sqrt(m0/m2)` of the reference
  elevation) across all four sea states and the transition record before any
  re-tuning is attempted.

### Stage D - per-axis regularisation instead of scalar anisotropy constants

- D1. Estimate world-horizontal acceleration variances from the world
  acceleration the direction stage already forms, with the same noise-floor
  subtraction as the vertical channel.
- D2. Derive `sigma_aw` and `r_S` per axis from those estimates, retiring
  `S_factor` and `R_S_xy_factor` as fallbacks used only until the per-axis
  estimator is ready. Anisotropy then follows the sea rather than a constant
  fitted on the smallest record.
- D3. Optional follow-up: feed the full 3x3 stationary covariance
  (`set_aw_stationary_cov_full` already exists) so the OU process is aligned
  with the estimated propagation axis.

### Stage E - re-validation and manuscript

Both committed evidence bundles pin a SHA-256 of the simulator sources they
were produced from, so *any* change under `src/util/` or `tests/*-sim.cpp` -
including the Stage A metric additions - invalidates them by design and
`tests/validation` fails until they are regenerated:

```
python3 tools/ou_robustness.py --mode full --bootstrap-resamples 10000 \
    --output-dir reports/results/ou_robustness
python3 tools/ou_validation.py --mode full
```

That is roughly 310 and 500 twenty-minute simulator runs respectively, so it is
a CI job, not a workstation job.

- E1. Re-fit coefficients with `tools/ou_tuning_sweep.py` on the stationary
  records, then re-run both bundles above and regenerate the generated `.tex`
  parts.
- E2. Update `doc/kalman_ou_iii/w3d-direction-methods.tex-part` and the results
  caption: report the travel-sense correctness rate and the heading-invariance
  experiment instead of the "not a correctness rate" caveat.
- E3. Update `docs/ou-validation.md` and the covariance-policy discussion to
  match the measured picture in section 1.1.

## 4. What this changes about the manuscript's claims

- "Adapting `tau` and `sigma_aw` makes things worse" is reproduced, but the
  cause is not adaptation. It is that the adaptation input is measured in the
  wrong spectral band, so the OU time constant is uninformative and `r_S`,
  which depends on its cube, lands far from what the sea requires.
- The horizontal-versus-vertical asymmetry is a tuning constant
  (`R_S_xy_factor`), not a property of the OU model or of the MEKF.
- The travel-sense result is stronger than reported, not weaker: the directed
  estimate is heading-invariant and correct to about 5 degrees. Only the
  exported class label was gauge-dependent.
