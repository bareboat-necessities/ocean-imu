# OU-III stability proof handoff — PR #639

## Controlling proof path

construction -> magnetically informed H18 -> compact A21 release -> outer A21
entry -> local V<=0.15^2 LaSalle strictness -> finite-error supply -> regime
composition -> float32.

Do not revive retracted scalar-reader/kernel-ceiling, coefficient-variation
Abel, AW-only passivity, source-uniform pre-entry pointwise AW tracking, or
"homogeneous zero action means base innovations are zero" arguments.

## Proved/retained

- Covariance-metric homogeneous operation nonexpansion and complete-word
  dissipation.
- Radius-local field-axis exclusion at r_FA=.15 and existence of local
  finite-window homogeneous strictness; no numerical eta_D.
- Direct H18/A21 release into the local ball is false as a theorem target.
- Correct BA storage elimination uses the outer marginal P_oo.
- Exact linked finite-error completed-square inequality
  V_N <= (1-gamma)V_0 + chi_gamma.
- EXCITED_MOVING candidate: 1 degree / 60 s, with a proposed assembled-device
  low-frequency residual qualification. Hardware qualification is OPEN.
- Continuous magnetic reference is canonical (h,0,z); exact reference-error
  decomposition is in the ledger. Reference cone remains OPEN.
- Default shipping sigma_aw floor is 9e-4 m/s2 after tuner readiness.
- Fresh 12-state LIN process covariance has a strictly positive qualitative
  lower floor by controllability + compactness; do not use its ~1e-29
  one-step scale as a contraction estimate.
- A21 release is first false->true acc_bias_updates_enabled transition, not
  Live handoff. Refinement has no forced timeout.

## Immediate blocker

Prove all-time held-BA LIN BIBO stability so arbitrarily delayed H18 release
still has a compact LIN mean.

A read-only diagnostic now exists:

    tools/stability/ou3_theorem/held_ba_lin_word_diagnostic.py

It propagates the exact 12x12 LIN root tangent over the final 17 s before BA
activation on carried quiet/wave histories, including prediction and literal
correction factors, with AG tangent columns treated as inputs. It is wired
into ou3-stability-proof as NON-PROMOTING evidence.

### Next calculation

1. Run/fix that diagnostic and record quiet/wave rho(M_LIN) and Euclidean
   singular norm.
2. Extend it to snapshot P_LL at the 17-s root and endpoint and compute the
   covariance-metric induced gain

       lambda_max(P0^(1/2) M' PN^(-1) M P0^(1/2)).

3. If carried gain is <1 with margin, prove source-uniform strictness using
   the fixed-factor nullspace/compactness argument over the bounded
   P_LL/tau/R_S/T_S/scheduler class. Do not promote the carried number.
4. If a LIN root has unit/noncontracting gain, inspect its eigenvector and
   identify the physical mode before doing interval refinement.
5. Once uniform LIN BIBO is proved, bound the affine AG/physical input over a
   complete word; attitude is compact on SO(3), BG estimate is projected, and
   physical/sensor inputs have declared bounds. This closes the principal
   noncompact H18->A21 release coordinate.

## Other open blockers after LIN release compactness

- Mahony proxy deterministic tilt bound and magnetic-reference cone E_B<15 uT.
- Real assembled AtomS3R/BMI270 LF residual qualification.
- Correct EXCITED_MOVING physical-to-nominal compatibility-line exclusion.
- Source-uniform outer A21 retention and finite entry below r=.15.
- Linked physical/model/arithmetic supply and prefix retention.
- Regime composition and float32 totality.

The stability-study LaTeX and docs/ou3-proof-research-state.md are controlling.
Finite carried diagnostics are evidence only and must not flip theorem-status
flags.
