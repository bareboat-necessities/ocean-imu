# Minimal (27-state) realization of the shared-modes QMEKF (experimental)

**Status: experimental, simulation only, not a default.**  The shipping
filter stays `SeaStateFusionFilter_OU_III`.  The 33-state shared-modes QMEKF
([kalman-carma-shared-modes.md](kalman-carma-shared-modes.md),
`KalmanCarmaSharedModes<T>`) is unchanged and bit-identical
(`carma_shared_modes-test` [2]); this is a research alternative to it.

* Motion block: `src/kalman_carma_shared_modes/CarmaSharedModesMinimalMotion.h`,
  `CarmaMinimalSharedModeMotion<T, R>`, six states per axis,
  `R in {Cascade (default), Quotient, QuotientHold}`.
* QMEKF core: `KalmanCarmaSharedModesMinimal.h`,
  `KalmanCarmaSharedModesMinimal<T = float, R = Cascade>`, the shared-modes
  core with this block (27 error states).
* Wrapper: `SeaStateFusionFilter_CarmaSharedModesMinimal.h`,
  `SeaStateFusion_CarmaSharedModesMinimal<R>`, the CARMA3 schedule.
* Tests and tools: `tests/kalman_carma_shared_modes/carma_shared_modes_minimal-test.cpp`
  (run by `make all`), the bench names `CSM27`, `CSM27Q`, `CSM27H`, the
  replay variants `W3D_CSM_VARIANT=CSM27...`, the host timing tool
  `kalman_carma_shared_modes_minimal-timing` (on demand), the ESP32-S3
  probe `tests/kalman_carma_shared_modes/arduino/csm_size_probe`.
* Study: `tools/carma_shared_modes_minimal_analysis.py` (mathematics),
  `tools/carma_shared_modes_minimal_study.py` (Monte Carlo, robustness,
  traces, dataset, timing, size, report) ->
  `reports/results/kalman_carma_shared_modes_minimal/`.
* Article: `doc/kalman_carma_shared_modes/kalman_carma_shared_modes-minimal.tex`.

RESULTS_SUMMARY

## 1. Why two directions per axis are unobservable

Per axis the 33-state block carries the complexified C41J mode
`zeta = B/(P+ D) xi` as eight real states `[Re(p, v, a, l), Im(p, v, a, l)]`,
and every accelerometer row reads `Re a`.  Each root `r` of the real
quadratic `D` (-7.93 and -1.25 rad/s at kappa = 1) is an eigenvalue of both
`A_c` and `conj(A_c)`, so the realified `A8` has it twice with a
two-dimensional eigenspace, which one scalar output cannot resolve.  Exactly:

    U = {0} x (span_C{w(r3), w(r4)} intersected with R^4),   w(r) = [B(r), r B(r), r^2 B(r), 1]

is the unobservable subspace: two real directions per axis, entirely in the
imaginary half, independent of the excitation `Sigma`, of the attitude and
of every QMEKF operation.  Its annihilator is the pair of rows whose
polynomial is `D(s)(alpha + beta s)`; they read `Im m`, `Im m'` of the
D-free mode `m = xi / P+`.  Two of the four imaginary directions stay
observable for every `Sigma` (a correction to the companion article's
observability paragraph, which calls `Im zeta` unobservable for real
`Sigma`).

## 2. Equivalence

* **Fixed parameters, any excitation sequence.**  `U` is invariant and
  annihilated by every measurement row (any rotation), the reset touches
  only attitude rows, and the fallback shrink
  `P_hz <- E + alpha (P_hz - E)` commutes with the projection.  The
  27-state QMEKF therefore gives the same attitude, bias, `p`, `v`, `a`,
  innovations and covariance as the 33-state filter, through amplitude
  commits, estimated coupling, confidence and fallback (test [3]:
  5-7e-12 relative in double precision over 400 s).
* **Three-dimensional law.**  Per axis `Re zeta = B/P * eta`,
  `eta = V [xi_r; xi_i]`, `V = [N_r, N_i] / D` with `P+ = N_r + i N_i`;
  `|D|^2 = |N_r|^2 + |N_i|^2` makes `V` lossless (co-inner).  Every
  realization of this `1 x 2` transfer function with the noise blocks
  `g2 [Re(Sigma_ab) I + Im(Sigma_ab) J2]` gives the same spectral density
  matrix: co-spectra, quadrature spectra, coherences and phases (test [2]:
  transfer functions 3e-14, stationary law and lag cross-covariance 9e-11).
* **Hermitian structure.**  `J` does not descend to the quotient, so the
  27-state block discretizes two intensities:
  `Q_ab = Re(Sigma_ab) Q_R + Im(Sigma_ab) Q_I`, `Q_R` symmetric from
  `g2 G G^T`, `Q_I` antisymmetric from `g2 G J2 G^T`.

## 3. Adaptation

The CARMA3 schedule commits a frequency scale `kappa` and amplitudes every
0.1 s, and the coupling changes `Sigma`.  A commit replaces `F` and `Q` and
keeps the state, so the adaptive model depends on the realization.

* **Excitation changes** (amplitudes, coupling, confidence, shrink) leave
  `U` and `F` alone: exact (above).
* **Frequency changes** are not.  Every realization scales exactly,
  `A(kappa) = kappa S A(1) S^-1` with a diagonal `S`, but
  `U(kappa) = S U(1)` moves: 0.45 deg per percent of `kappa`, 10 deg at
  +-25 %.  After a commit the 33-state filter's component in the old `U`,
  which it carries with a data-dependent posterior, leaks into the new
  observable coordinates.  No 27-state filter carries it, so none is
  exactly equivalent under frequency adaptation (motion-only experiment:
  `|p_33 - p_27Q| = 3e-5 m`, not rounding).
* **Quotient**: at a kappa commit the missing component is replaced by its
  stationary regression on the quotient (`L = P8 Pi^T (Pi P8 Pi^T)^-1`), and
  `x <- Pi(k_new) L(k_old) x` on every axis (`p, v, a` exactly continuous).
  With the scaling this is one double factorization at construction.
* **QuotientHold** holds the quotient coordinates (ablation).
* **Cascade** (default): the CARMA3 baseline's own `[p, v, a, l]` driven at
  `a` by `eta = xi_r + C_w w`, and the two colouring states `w` of `V` in
  input-normal coordinates (`A_w + A_w^T = -B_w B_w^T`, `C_w = -B_w e_1^T`,
  `B_w` upper triangular).  A commit holds `[p, v, a, l]` (as A41 does) and
  the normalized `w / sqrt(g2 sigma_a^2)`.  Because `V` is lossless with a
  constant Gramian, `eta` stays exactly white under any `kappa(t)`,
  `sigma(t)`: **with the coupling off the cascade is the CARMA3 baseline A41
  in law, adaptation included** (test [4]: 1.6e-13 relative; the 33-state
  filter and the quotient differ from A41 by up to 22 % under the same
  schedule).  This removes the "realization effect" documented for
  `CSM_OFF` in [kalman-carma-shared-modes.md](kalman-carma-shared-modes.md)
  section 2 by construction.  With the coupling on it is a third adaptive
  model, judged by the experiments below.

## 4. Is 27 the right dimension?

The minimal Markov order of the three-axis displacement law is
`4 rank(Sigma) + 2 rank(Im Sigma)` (residues of the spectral density: at the
roots of `D` the real part cancels and the residue is proportional to
`Im Sigma`, a real antisymmetric matrix of rank 0 or 2; Hankel singular
values confirm 16, 16, 12, 12, 8 for full Hermitian, quadrature, real,
diagonal and rank-one excitations).

* Real (in-phase) coupling or none: 12 motion states, i.e. the CARMA3
  baseline or Experiment A.
* Quadrature coupling (linear waves): 16 motion states, a 25-state filter,
  with the extra four in the vertical plane of propagation.
* 27 is the minimal realization that keeps the shared dynamics
  `F = I3 (x) F6` for every excitation.

A 25-state filter exists (C41J on every axis plus two colouring copies in
the plane of `range(Im Sigma)`), but its coordinates depend on `Sigma`: they
would move at every coupling, confidence and heading update, which loses
the exactness under excitation changes.  It is undefined at
`Im Sigma = 0`, which the confidence crosses at every warm-up and fallback,
and it loses the Kronecker structure.  All that for two states, about 14 %
of the 27-state covariance work.  Not implemented.

## 5. Implementation notes

* Paired discretization `discretize_lti_pair` / `StepDiscretizationPair`:
  `discretize_lti`'s scaled series, doublings and jitter series for
  `(Q_R, Q_I)`.
* Per-axis re-coordinatization hook `take_axis_recoordinatization` in the
  CARMA3 core (`Kalman3D_Wave_CARMA3::recoordinate_motion_`, applied as
  `S P S^T` right after a commit).
* Exact structural zeros.  Round-off residue of `Pi A8 L` (down to 1e-19)
  defeated the discretization's structural scaling and produced denormals,
  a 7x slower refresh on x86.  The cascade's `B_w` is exactly upper
  triangular and `C_w = -B_w e_1^T` exactly.
* Closed-form 2 x 2 algebra (Lyapunov, Cholesky, Givens) in the cascade's
  construction.  Eigen's generic decompositions added 44 kB of code on the
  ESP32-S3.
* Conditioning (stationary correlation per axis): 147 for the 33-state
  block, 128 for the quotient, 5.5 for the cascade (= A41: the colouring
  states are white and uncorrelated with `[p, v, a, l]`).

RESULTS_DETAIL
