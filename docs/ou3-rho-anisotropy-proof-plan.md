# OU-III horizontal rho anisotropy: proof-track follow-up

The shipping filter now runs an anisotropic horizontal integral regularizer,
`rho_x = 0.72` and `rho_y = 0.50`. The proof track was written against the
isotropic pair and has not been moved. This file is the handover for that work.

Read `AGENTS.md`, `docs/ou3-brmm-main-handover.md` and
`docs/ou3-proof-research-state.md` first. Start from the latest `main` and open
a new PR. `P3 = 1e-18` stays frozen; nothing here licenses relaxing it.

## What shipped, and why

`SeaStateFusionFilter_OU_III.h` carries the measurement in full. In short: the
pinned v1.2.1 records are an estimated fin-keel sailboat response, and a keel
resists sway while leaving surge comparatively free.

A fixed `+/-30` degree heading projection would put the same horizontal RMS
ratio `cos30/sin30 = 1.732` in every record. The records instead run 2.33 at
`H1.5` down to 1.61 at `H8.5`, falling with wavelength and reaching the
geometric value only in the longest waves. A heading projection cannot vary
with period; a hull response can. Every record also holds `yaw` at exactly 0
(mean 0.000, sd 0.000, all eight), so the world horizontal axes are the
vessel's surge and sway here, and the knob pair is the surge/sway split.

Pooled over four fresh IMU draws and the eight records: pitch -17 percent,
lateral accelerometer bias -21 percent, 3D accelerometer bias -2 percent, yaw
unchanged, against roll +4 percent and x bias +5 percent. Both vertical
channels are flat at 1.000.

The body-frame form of this is **not** claimed. The keel's anisotropy is a
body-frame property and `R_S` is applied in world NED, so the two coincide only
while the vessel holds the heading these records hold. Rotating the split by
the estimated heading is the deployment-general statement and cannot be
measured on records that pin yaw.

## What this breaks

Three classes. The second is the one that will be missed, because it stays
green.

### A. Hard failure from one named premise

`tools/stability/ou3_brmm_complete_source.py:103` asserts the deployed source
text contains both factors at `0.72`:

```python
"horizontal_RS_factors_are_0p72": (
    "float R_S_x_factor_ = 0.72f;" in text
    and "float R_S_y_factor_ = 0.72f;" in text
),
```

That premise feeds `P3_source_contract_ready`, and its failure cascades. In the
ALT suite 51 of 55 errors reduce to it:

```
RuntimeError: correlated BRMM outer-enclosure prerequisites failed:
  {'complete': ['P3_source_contract_ready is not true',
                'R_S source parity failed: horizontal_RS_factors_are_0p72']}
```

Same file, `:363` emits and `:510` asserts `axis_std_factors == [0.72, 0.72,
1.0]`. `ou3_brmm_p3_full_preconditions.py:106` repeats the source-text match.
`ou3_brmm_full_normal_live_word_reset.py:71` runs its interval reset word with
`rs_std_xyz = [0.72, 0.72, 1.0]`.

### B. Silently stale certificates

These still build and validate, because each hard-codes `0.72` **and** checks
its own literal. They now certify bounds for a configuration that is not
shipped, and CI will not say so:

- `ou3_p4_innovation_binary32_bounds.py:52` computes
  `'S_zero': ((0.72*rslo)**2, rshi**2)`; `:82` emits
  `actual_RS_horizontal_factor: 0.72`; `:95` fails only if that emitted literal
  changes.
- `ou3_p4_marginal_correction_obstruction.py:86` computes
  `rs_std_x = down(0.72 * rs_lo)`; `:131` emits and `:168` guards the same way.

Both also carry `actual_RS_horizontal_factor` as a **single scalar**. An
anisotropic pair is not a value change to these modules, it is an interface
change: the `S_zero` innovation range and the marginal-correction obstruction
have to be re-derived with separate x and y horizontal standard deviations,
and the x branch is no longer the binding one by construction.

### C. Adapts on its own

`ou3_brmm_tuner_scheduler_step.py:227` regexes the deployed value out of the
header and needs no edit. This is the pattern the class A and class B modules
should adopt: read the deployed constants, do not restate them.

`ou3_brmm_riccati_tube_factored` is unaffected and still closes. It is what
`tests/validation/test_ou3_brmm_riccati_tube.py` actually exercises, and what
the P4 modules in class B import as `TUBE`, so do not confuse it with the
unfactored module in class D below.

### D. Pre-existing certificate failures, not caused by this change

These two regex the deployed factor correctly and the certificate then does not
close on it:

```
ou3_brmm_riccati_tube.py:254
  RuntimeError: cannot certify scaled OU process cell
                [0.009999999068167651, 0.010000000000000037] at depth 20
ou3_source_reachable_matrix_p3.py:247
  RuntimeError: cannot certify scaled OU process cell
                [0.00041666665735344007, 0.00041667021815050694]
```

Reproduce with `build()` on the default domain from `tools/stability`:

```
python3 -c "import ou3_brmm_riccati_tube as m; m.build()"
python3 -c "import ou3_source_reachable_matrix_p3 as m; m.build()"
```

**Attribution is settled: pre-existing.** Both builds were run at the parent
commit `238a70e`, where the header still carried the isotropic `rho_y = 0.72`,
and both fail there with byte-identical cell intervals. The obstruction is in
the scaled OU process cell itself and is independent of the horizontal factor,
so it is separable from the rho work and must not be folded into it. Fix or
write it up as its own change under the `AGENTS.md` taxonomy.

Neither is reached by `test_ou3_brmm_riccati_tube`, which exercises
`ou3_brmm_riccati_tube_factored` instead, so CI does not see either failure
today and no workflow turns red because of them.

## The work

1. Generalise the parity premise. Replace the `0p72` literal match with an
   extraction of both deployed factors, and carry them as a pair through
   `axis_std_factors`. Keep it fail-closed: a factor that is absent,
   non-finite, non-positive or outside the setter's `[0, 4]` clamp must still
   fail. Rename the premise; `horizontal_RS_factors_are_0p72` will be wrong
   again the next time the deployment moves.
2. Re-derive the two numeric bounds in class B for an anisotropic horizontal
   pair, and widen `actual_RS_horizontal_factor` to x and y. Do not simply
   substitute `min(rho_x, rho_y)`: check whether each bound is monotone in the
   factor before taking a worst case, and if it is not, carry both axes.
3. Leave class D alone. It is pre-existing and independent of `rho`; closing
   those two certificates is separate work with its own change. Do not let it
   block items 1 and 2, and do not let a green run of them be read as evidence
   that the rho work is incomplete.
4. Re-pin the six SHA gates on the filter header if it changes again. They are
   `finite_mahony_prefix_totality.py`, `finite_wpe_frequency_binary32.py`,
   `finite_mag_counter_saturation.py`,
   `finite_startup_disturbance_obstruction.py`,
   `ou3_brmm_infinite_continuation.py` and `ou3_brmm_magnetic_call_schedule.py`.
   They were re-pinned for this change and all 84 of their tests recompute and
   close, so they are current; they are listed here only because any further
   header edit trips them again.
5. Regenerate the OU validation and robustness evidence. The change moves two
   replay dependencies, `src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h` and
   `tests/kalman_ou_iii/kalman_ou_iii-sim.cpp`, so `ou_evidence_contract`
   reports both studies as stale until `ou-validation` regenerates and commits.

## Acceptance

- `ou3_brmm_complete_source.build()` validates with the deployed pair, and
  fails when either factor is edited away from what the header carries.
- The full `tests/ou3_alt_contraction` suite returns to its pre-change state:
  the only remaining errors are the four `mpmath` imports and the
  `test_carried_storage_rho_diagnostic` import, which are environment and
  predate this work.
- `ou3_p4_innovation_binary32_bounds` and
  `ou3_p4_marginal_correction_obstruction` report the deployed x and y factors
  and fail when either drifts from the header.
- Class D is unchanged: both builds still fail exactly as they do at
  `238a70e`, and no attempt to close them rides this change.
- `ou3-proof`, `ou3-complete-brmm` and the `ou3-p4-*` workflows are green.
- No theorem gate moves. `P4_PASS`, `P5_MAY_START` and every `ALT_LIVE_PASS`
  stay false.

## Do not

- Do not revert the filter to isotropic `rho` to make the proof close. The
  anisotropy is a measured hull property and the measurement is in the header.
- Do not weaken the parity premise to "any value is acceptable". The point of
  it is that the proof is derived for the configuration that ships.
- Do not claim the body-frame rotation. Records that pin yaw cannot measure it,
  and it is a separate obligation with its own dataset requirement: records at
  a second heading, `0` or `60` degrees, would separate hull from fixture.
