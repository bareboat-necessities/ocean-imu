# Three-axis CARMA on the two-frame Lie group (TFG-CARMA3), experimental

**Status: experimental research.**  Nothing here changes the shipping OU-III
estimator or the existing CARMA3 and TFG filters.  `Kalman3D_Wave_OU_III` /
`SeaStateFusionFilter_OU_III` remain the deployed, authoritative
implementation; nothing in the deployed configuration includes the files
below.  All results are simulation results; no hardware was used.

Files:

* `src/kalman_carma/Kalman3D_Wave_TFG_CARMA3.h` -- the filter,
  `Kalman3D_Wave_TFG_CARMA3<T, NC>` (NC = 4: CARMA(4,1), 21 error states).
* `src/kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h` -- the wrapper: the
  CARMA3 wrapper with the inner filter replaced (a defaulted template
  parameter of `SeaStateFusionFilter_CARMA3`; CARMA3's own instantiations are
  unchanged, see section 4.4).
* `tests/kalman_carma/tfg_carma3-test.cpp` -- the mathematical verification
  (section 5).
* `tests/kalman_carma/carma3_bench.cpp --extended`, `carma3/FullAdaptersTFG.h`,
  `carma3_filters_tfg.cpp`, `carma3-sim` (`W3D_CARMA3_VARIANT=T41`) -- the
  benchmark.
* `tools/tfg_carma3_benchmark.py`, `tools/tfg_carma3_report.py` -- the Monte
  Carlo driver and the statistics/figures; results in
  `reports/results/tfg_carma3/`.

**Result in brief** (simulation; sections 7-13; generated tables in
`reports/results/tfg_carma3/summary.md`):

* TFG-CARMA3 is a genuinely different estimator -- 21 right-invariant
  error states on `TwoFrameGroup<T, 4, 2>`, estimate-independent
  measurement Jacobians, group injection, reset transport -- verified against
  independent references in 199 checks (section 5), with CARMA3's CARMA(4,1)
  prior reproduced exactly.
* It is first-order equivalent to CARMA3's quaternion MEKF (section 2.8),
  and translation shows it: in 34 of 36 conditions every displacement and
  velocity change is inside +-1 % with its whole 95 % interval, usually
  inside +-0.01 % (100 paired seeds on the principal seas, 30 or 20
  elsewhere, 80 independent dataset replays).  CARMA3's weaknesses on the
  still IMU and in the mixed sea are unchanged.
* Where the CARMA prior misfits (mixed sea, manoeuvring) the invariant
  error coordinates give better attitude and bias estimates -- tilt 2.6-6 %
  and accelerometer bias 3-7 % lower (significant after correction),
  heading up to 13 % lower (raw p = 7e-4; 24 % without reset transport,
  30 of 30 seeds) -- and less overconfidence (attitude NEES 7.6-12 %
  lower); with CARMA(5,2) the MEKF's mixed-sea
  attitude fails (tilt 1.6 deg) and the no-transport TFG keeps 0.34 deg.
  The cause is traced to the estimate dependence of the MEKF's coordinates
  (section 8); the reset transport gives back part of the gain.
* Through a magnetometer outage the TFG's heading covariance is far more
  honest (no-transport heading NEES 1.2 over 2400 s against CARMA3's 141 and
  OU-III's 2570) and it reacquires heading faster; its heading is not more
  accurate during long outages.
* Unfavourable: 3.4 times CARMA3's multiply-adds (1.13 times OU-III's;
  1.5 times CARMA3's host time); end-of-run gyro-bias error up to 0.8 %
  higher in regular seas and 17 % higher when manoeuvring (not
  significant after correction); no ESP32-S3 build evidence (section 11).
* Recommendation (section 13): do not replace CARMA3's Q-MEKF or deploy;
  continue only the no-transport variant as a targeted experiment against
  an observability-constrained CARMA3 MEKF.

## 1. Hypothesis and scope

CARMA3 (`docs/carma-three-axis.md`) runs a CARMA(4,1) prior on all three
world axes under OU-III's quaternion multiplicative EKF (Q-MEKF): additive
errors on the world-frame CARMA states and the body-frame biases, a
body-frame attitude error.  The two-frame Lie-group filter
(`docs/tfg-design.md`) carries the world states and the body biases in one
group with a right-invariant error.  The hypothesis under test:

> Representing attitude, the CARMA translational states and the IMU biases
> on the right-invariant two-frame group yields measurably better estimation
> (accuracy, robustness, bias tracking or covariance consistency) than the
> Q-MEKF with the same three-axis CARMA(4,1) prior.

The comparison holds everything else fixed: the same CARMA model library,
shape, frequency/amplitude mapping and exact discretization; the same
wrapper source code (front end, proxy startup, staged tune commit, tilt
watchdog, magnetometer-service bias gate, vibration Racc inflation); the
same sensor streams; the same physical initial uncertainty (section 4.3).

## 2. Mathematical model

### 2.1 State, group and error

Nominal state `X = (R, W, B)` in `TwoFrameGroup<T, 4, 2>`
(`src/lie/TwoFrameGroup.h`):

    R = R_bw in SO(3)                body -> world (NED, z down)
    W = [p v a l] in R^{3x4}         world frame; row k = CARMA state of axis k
                                     in the library coordinates [p, v, a, latent]
    B = [b_g b_a] in R^{3x2}         body frame

    (R1,W1,B1) o (R2,W2,B2) = (R1 R2, W1 + R1 W2, B2 + R2^T B1)

Convention contract: `R = R_bw`, right-invariant error
`eta = Xhat X^-1 = Exp_G(xi)`, correction on the left
(`Xhat+ = Exp_G(-xihat) Xhat`).  Local coordinates

    dR = Rhat R^T = Exp(phi),   rho_i = what_i - dR w_i,   beta = R (bhat - b)
    xi = [ phi | beta_g | rho_p rho_v rho_a rho_l | beta_a ]     21 states
           0     3        6     9     12    15      18

so the attitude error is world-referred and the body biases are carried as
world-referred vectors whose body-frame meaning is fixed by the true
attitude.  This is the same group as Kalman3D_Wave_TFG's (four world columns,
two body columns); only the world columns are relabelled.

### 2.2 Nominal dynamics

With measured rate `w_m`, gyro white noise `n_g`, bias walk `n_bg`, CARMA
driving noise `xi in R^3` (spatial covariance `Sigma_s`) and the residual
accelerometer-bias OU process of CARMA3:

    R'   = R [w_m - b_g - n_g]x
    W'   = W A^T + g e_a xi^T                x_k' = A x_k + g e_a xi_k per axis
    b_g' = n_bg
    b_a' = -b_a / tau_b + n_ba  (active),    b_a' = 0 (held: phi_e = 1)

`A` is CARMA3's CARMA(4,1) realization at the committed frequency scale
`kappa` with unit acceleration variance (`CarmaVerticalModel`), `g^2 =
kappa / Var_1(a)`, and `Sigma_s = [[s_x^2, r s_x s_y, 0], [r s_x s_y, s_y^2, 0],
[0, 0, s_z^2]]` the per-axis amplitudes (the wrapper's mapping gives
`s_x = s_y = s_z`, `r = 0`).  The discrete mean propagation over a step `h`
with zero-order-held rate is exact: `R <- R Exp(w h)`, `W <- W Phi_A^T`,
`b_a <- phi_b b_a` with `w = w_m - bhat_g`, `Phi_A = exp(A h)` from the
library's scaled-series discretization, `phi_b = exp(-h/tau_b)` (active)
or 1 (held).

### 2.3 Invariant error dynamics

Exactly (no linearization), with `u = Rhat (b_g - bhat_g + n_g)`:

    eta_R'  = [u]x eta_R
    eta_W'  = eta_W A^T - [u]x (What - eta_W) - eta_R g e_a xi^T
    eta_Bg' = [R (w_m - b_g)]x eta_Bg - R n_bg
    eta_Ba' = [R (w_m - b_g)]x eta_Ba - eta_Ba / tau_b - R n_ba

Linearized at `xi = 0` (`u = -beta_g + Rhat n_g`, `w_w = Rhat w`):

    phi'    = -beta_g + Rhat n_g
    beta_g' = [w_w]x beta_g - Rhat n_bg
    rho_i'  = sum_j A_ij rho_j - [what_i]x beta_g + [what_i]x Rhat n_g - g e_a,i xi
    beta_a' = [w_w]x beta_a - beta_a / tau_b - Rhat n_ba

Two structural facts follow.

* **No attitude -> translation term.**  `eta_W' ` contains `eta_W A^T` and
  `[u]x What`, but not `eta_R` acting on `W` through `A`: with one `A` on
  every axis the column mixing commutes with rotation,
  `(Q W) A^T = Q (W A^T)`, so the term `A (phi x w) - phi x (A w)` cancels
  identically.  `Phi_{rho,phi} = 0` (checked, < 1e-15).
* **Not group-affine.**  The `beta_g -> rho` coupling `-[what_i]x` depends on
  the estimate.  This filter is therefore a *right-invariant two-frame
  Lie-group error-state EKF*, not a log-linear invariant EKF; no
  estimate-independence of the error dynamics is claimed.

### 2.4 Transition over one step

With `R0`, `R1 = R0 Exp(w h)` the attitude at the start and end of the step
and `W1 = W0 Phi_A^T`:

    Phi_{phi,phi}   = I
    Phi_{phi,bg}    = F = -R0 h J_l(w h) R0^T          (int_0^h R(t) dt = R0 h J_l(w h))
    Phi_{bg,bg}     = E = R1 R0^T = Exp(w_w h)
    Phi_{rho,rho}   = Phi_A (x) I_3                     (exact CARMA transition, every axis)
    Phi_{rho_i,bg}  = [w_i(h)]x F
    Phi_{ba,ba}     = phi_b E

**The collapse identity.**  The structured impulse response
`D(s, h) = int_s^h (Phi_A(h-t) (x) I) [W(t)]x R(t) dt` that
Kalman3D_Wave_TFG integrates by five-point Gauss-Legendre quadrature
collapses for a common `A`: because the skew map is linear in each column
and `Phi_A(h-t) Phi_A(t) = Phi_A(h)`,

    (Phi_A(h-t) (x) I) [W(t)]x = [W(t) Phi_A(h-t)^T]x = [W(h)]x     for every t,

so `D(s,h) = [W(h)]x int_s^h R(t) dt` and every `rho` column of the
transition and of the noise injection is `[W(h)]x` times the attitude
column.  An independent derivation gives the same result: the additive
errors `dw = w - what` obey `dw' = A dw` with no attitude or bias term, and
`rho_i = -dw_i + [what_i]x phi` to first order, so
`rho_i(h) = (Phi_A rho(0))_i + [what_i(h)]x (phi(h) - phi(0))`.  The
quadrature of the existing TFG is therefore unnecessary for this prior; the
test checks the closed form against the quadrature (1e-9) and against a
fine-step integration of the continuous model (section 5).

### 2.5 Process noise

Noise inputs (body gyro density `Q_g`, bias walk `Q_bg`, CARMA `W_A (x)
Sigma_s`, bias OU `Q_ba`) enter as

    gyro:      phi += R(s) n_g,          rho_i += [w_i(s)]x R(s) n_g
    gyro bias: beta_g(h) += -R(h) n,     phi(h) += C(s) n,  C(s) = int_s^h R(t) dt
    CARMA:     rho -= (g e_a (x) I) xi
    accel bias: beta_a(h) += -phi_b(h-s) R(h) n_ba

and after the collapse the discrete covariance needs only 3 x 3 integrals:

    Q_pp = int R Q_g R^T + int C Q_bg C^T        (Q_g isotropic: q_g h I exactly)
    Q_pb = -(int C) Q_bg R1^T,   Q_bb = h R1 Q_bg R1^T
    Q_{rho_i, .}       = [w_i(h)]x Q_{phi, .}
    Q_{rho_i, rho_j}   = [w_i(h)]x Q_pp [w_j(h)]x^T + Q_A(i,j) Sigma_s
    Q_aa = R1 Q_ba R1^T (tau_b/2)(1 - exp(-2h/tau_b))   (active; 0 held)

`Q_A` is the library's exact CARMA process covariance at unit amplitude, so
the CARMA part is `Q_A (x) Sigma_s` exactly, including the horizontal
cross-covariance when `r != 0`; the gyro/bias integrals use five-point
Gauss-Legendre (exact to degree 9 in `s`).  There is no cross-covariance
between the CARMA noise and the gyro noise (independent sources); the
gyro-noise <-> translation and gyro-bias <-> translation cross-covariances
are the `[w_i(h)]x` blocks above and are all kept.

### 2.6 Measurements

Accelerometer (temperature-compensated bias as CARMA3) and magnetometer, in
the world frame:

    r_a = Rhat (f_m - bhat_a - k_a dT) - (ahat - g),    noise Rhat R_a Rhat^T
    r_m = Rhat m - B_w,                                  noise Rhat R_m Rhat^T

With the truth `X = Exp_G(-xi) Xhat` and a noise-free measurement,
`r_a(xi) = -rho_a + (I - Exp(phi)) g - eta_R eta_Ba`, so at `xi = 0`

    H_a = [ [g]x | 0 | 0 0 -I 0 | -I ],      H_m = [ -[B_w]x | 0 ... ]

exactly: these Jacobians do not depend on the estimate.  (The world residual
is `Rhat` times CARMA3's body residual, so the innovation information is
identical; only the error coordinates differ.)  CARMA3's semantics are kept:
no innovation gating; a held accelerometer bias stays in `S` but not in the
cross-covariance or the gain.

### 2.7 Correction and covariance transport

`a = K r` is the estimated error; the estimate is corrected by
`Xhat+ = Exp_G(-a) Xhat` (attitude correction rotates every world column and
re-expresses the biases), and the posterior covariance (CARMA3's Joseph
form) is transported into the new error coordinates by

    J_r(a) = d/dd Log_G( Exp_G(-a) Exp_G(a + d) ) |_{d=0}
           = sum_n (-ad_a)^n / (n+1)!
           = [[ J_l(-phi), 0 ], [ Q_r(phi, z), J_l(-phi) ]]   per vector block z

with `Q_r(phi, z) = Q_l(-phi, -z)` the SE(3) left-Jacobian block (Barfoot,
eq. 7.86), coefficients evaluated from their series below `|phi| = 1`.
Every vector block (world column or world-referred bias) has the same
`ad` structure because `B~ = R B` composes like a world column.  The
filter applies `P <- J P J^T` structurally (block diagonal plus the attitude
column).  Ablations: `AttitudeOnly` (J_l(-phi) on the attitude block, the
MEKF analogue) and `None`.

### 2.8 Relation to the quaternion MEKF

CARMA3's error coordinates are `e = [dtheta_b, db_g, dw, db_a]` with
`C_true = Exp(dtheta_b) Chat` (`C = R_wb`) and additive truth-minus-estimate
errors elsewhere.  To first order

    phi = Rhat dtheta_b,   beta = -Rhat db,   rho_i = -dw_i + [what_i]x Rhat dtheta_b,

a linear map `xi = T(Xhat) e` that depends on the estimate.  Both filters'
linear models are exact first-order linearizations of the same system, so

    Phi_TFG = T(Xhat1) Phi_MEKF T(Xhat0)^-1,  Q_TFG = T Q_MEKF T^T,  H_TFG = Rhat H_MEKF T^-1,

hence `K_TFG = T K_MEKF Rhat^T` and `a = T de` when the prior covariances
correspond.  **The two filters are first-order equivalent**: they can only
differ through (i) the injection, which agrees to first order in the
correction, (ii) the covariance transport, and (iii) the estimate at which
`T` is evaluated before and after a correction.  Working out (ii)-(iii) in
CARMA3's coordinates, the TFG's update implies, to first order in the
correction,

    dw_i+ = dw_i - (1/2)[phi_a]x dw_i + (1/2)[Delta what_i]x phi

against CARMA3's `dw_i+ = dw_i`.  The term is first order in the *state
correction* `Delta what` times the attitude error, so the filters separate
where the world states are corrected by large amounts -- the first seconds
after Live, when `ahat` jumps from 0 to the measured acceleration -- and
agree closely in steady state.  Section 5 measures it.

### 2.9 Unobservable heading

Without a magnetometer, rotating the world about the vertical is a symmetry
(gravity and the isotropic horizontal prior are invariant).  In the
invariant coordinates the unobservable direction is `e_phi_z` for every
estimate: `H_a e_phi_z = [g]x e_z = 0` and `Phi e_phi_z = e_phi_z`.  The
information along it, `e^T P^-1 e`, therefore cannot grow under prediction
(adding `Qd` only removes information) or under the accelerometer update
(`P+^-1 = P^-1 + H^T R^-1 H`, and also for CARMA3's held-bias update, by the
Woodbury identity).  In CARMA3's coordinates the direction is
`n(Xhat) = T^-1(Xhat) e_phi_z = [Rhat^T e_z, 0, [what]x e_z, 0]`, which moves
with the estimate; each update adds information along the new direction
(the classical source of false heading observability).  The full reset
`J_r(a)` re-introduces an estimate dependence in the TFG as well: the update
is linearized at `xi = 0`, the posterior sits at `a`, and `J_r(a)` maps the
symmetry direction at `a` (`J_r(a)^-1 e_phi_z`, which `H_a` does not
annihilate) onto `e_phi_z`.  Only the TFG without transport preserves the
direction exactly.  Section 5 quantifies all three.

### 2.10 Rotational equivariance and anisotropy

* **Deterministic model.**  For one `A` on every axis, `(Q W) A^T = Q (W A^T)`
  for every rotation `Q`: the CARMA transition commutes with any world
  rotation (checked, 1e-14).
* **Stochastic model.**  The driving noise transforms as
  `Sigma_s -> Q Sigma_s Q^T`, and `Qd(Q Sigma Q^T) = (I (x) Q) Qd(Sigma) (I (x) Q)^T`
  exactly.  The prior is invariant under `Q` only if `Q Sigma_s Q^T = Sigma_s`:
  under every rotation for isotropic `Sigma_s`, under yaw for
  `s_x = s_y, r = 0`, and under no rotation otherwise.
* **The filter.**  Gravity fixes the vertical, so the complete filter is a
  yaw-equivariant system: a run in a world frame rotated by `Rz(psi)` (state,
  covariance, magnetic reference and `Sigma_s` rotated) reproduces the
  original run rotated, to 1e-12.  With an anisotropic horizontal `Sigma_s`
  that is *not* rotated with the frame the runs differ (correlation units
  2.2): the prior then carries heading information.
* **Anisotropic dynamics.**  If the axes had different `A` (different
  horizontal and vertical shape or scale), the transition would not commute
  with rotation and the invariant error would acquire

      rho' = A_axis(rho) + C(W) phi,   C(W) phi = A_axis([phi]x W) - [phi]x A_axis(W),

  an attitude -> translation coupling with `Phi_{rho,phi} = int Phi(h-t)
  C(W(t)) dt` (verified against finite differences of the nonlinear error
  flow, 3e-10; 0.25 over 50 ms for a 20 % scale difference).  The latent
  coordinates would no longer form world vectors either.  This
  implementation requires one frequency scale and shape for all axes (the
  wrapper's mapping provides it) and counts, rather than approximates, any
  other configuration (`rejected_anisotropic_configurations()`, zero in
  every run).

## 3. Preservation of the CARMA prior

The filter uses `carma::CarmaVerticalModel<T, 4>` and
`carma::StepDiscretization<T, 4>` unchanged: the fitted CARMA(4,1) JONSWAP
shape (zero at DC, `f_z0 = 0`), the frequency scale `kappa = f_tune / f_z,model`,
the amplitude `g^2 = sigma_a^2 kappa / Var_1(a)`, the exact scaled-series
discretization with doublings and the reference-step Taylor extension for a
jittered clock.  The wrapper is the same code, so the mapping
`sigma_h = S_factor c_carma max(sigma_floor, sigma_aw / c_sigma)` and every
clamp are CARMA3's.  There is no `S` state and no `S = 0`
pseudo-measurement; the fourth CARMA state is the latent coordinate of the
cascade realization.  Section 5 checks stability, positive semidefiniteness,
stationary covariances and spectra of the three-axis process as the filter
discretizes it.

## 4. Implementation

### 4.1 Filter

`Kalman3D_Wave_TFG_CARMA3<T, NC>` stores the nominal state as a
`TwoFrameGroup<T, NC, 2>` and a dense symmetric `NX x NX` covariance
(NX = 21 for CARMA(4,1)).  Per sample:

* **Prediction** (`time_update`): exact mean propagation; `P <- Phi P Phi^T`
  as two structured passes of `Phi .` (the sparsity of section 2.4), then the
  structured `Qd` of section 2.5; exact symmetry by mirroring.
* **Update** (`measurement_update_acc_only`, `measurement_update_mag_only`):
  `P H^T` from the sparse `H`, `S`, `K = P H^T S^-1` through three solves,
  CARMA3's exactly symmetric Joseph form, the structured `J_r(a) P J_r(a)^T`,
  the `Exp_G(-a)` injection, CARMA3's gyro-bias ball and accelerometer-bias
  norm projections (means only).
* Every covariance-sized temporary is a member (no dynamic allocation after
  construction, nothing NX x NX on the stack), all scalar code is templated
  (float on the target, double in the tests) and `-Wconversion` clean.

### 4.2 Policy parity with CARMA3

The wrapper is the CARMA3 wrapper; the core exposes CARMA3's interface.
Every covariance *write* that CARMA3 performs in its own coordinates --
the attitude seeds of `initialize_from_attitude` / `initialize_from_acc`
(body-frame tilt/yaw block, attitude and gyro-bias decorrelation), the
stationary CARMA prior (`reset_aw_covariance_to_stationary`), freezing and
unfreezing the accelerometer bias, the absolute attitude writes of the
proxy's heading gauge (`set_quaternion_boat`) -- is applied to
`T^-1 P T^-T` and mapped back (`CovarianceWriteMode::Physical`), so both
filters represent the same physical uncertainty at every phase change.
`CovarianceWriteMode::Native` applies the same writes directly to the
invariant covariance (ablation `T41NW`).

### 4.3 What differs, structurally

Only the estimator: error coordinates, the transition/noise in them, the
group injection and the reset transport.  The measurement information, the
measurement-acceptance policy (no gating; held-bias semantics), startup,
adaptation, magnetic service and the CARMA prior are identical.  The
controlled ablations isolate the remaining design choices: reset transport
(`T41R0`, `T41RA`) and covariance-write coordinates (`T41NW`).

### 4.4 CARMA3 and OU-III unchanged

`SeaStateFusionFilter_CARMA3` gained a defaulted template parameter for the
inner filter; every existing instantiation names the same type as before.
Built with strict IEEE evaluation (`-O2 -ffp-contract=off`, no
`-march=native`), the benchmark rows of OU3 and A41 at this branch and at
the baseline `9676718b4591817692f8e67c8501158f6cf45339` are bit-identical in
all 63 columns.  With the repository's `-O3 -march=native` flags they
differ in the sixth significant digit because GCC contracts multiply-adds
differently once the adapter translation units change; that is a property
of the build, not of either filter.

## 5. Numerical verification

`tests/kalman_carma/tfg_carma3-test` (199 checks, all pass; part of
`make all`; full log `reports/results/tfg_carma3/math_checks.txt`; also
199/199 with no report when built with `-O1 -fsanitize=address,undefined`,
which the repository's sanitizer gate does not cover for this directory).  No
expected value is "what the implementation did": every reference is built
independently -- finite differences of the nonlinear group maps, fine-step
RK4 of the stated continuous model, the structured Gauss-Legendre quadrature
of Kalman3D_Wave_TFG, the reset Jacobian of the independent
Kalman3D_Wave_TFG implementation (same group), Monte Carlo of the nonlinear
stochastic model, and Kalman3D_Wave_CARMA3 run side by side.  Measured
discrepancies (worst over the cases of each family):

| Check | Cases | Worst measured | Limit |
|---|---|---|---|
| Group: `X X^-1 = I`, associativity, `Log_G(Exp_G(xi)) = xi` (`|phi|` 1e-6..2.5), right invariance of `eta`, retract/local | 50 random elements | 4.9e-15 | 1e-11 |
| Bias slot of `Xhat X^-1` equals `R (bhat - b)` exactly | 50 | 2.8e-17 | 1e-14 |
| `J_r(a)` vs FD of `Log_G(Exp_G(-a) Exp_G(a+d))`, `|a|` ~ 1e-3 / 0.1 / 1 | 30 | 9e-13 / 1.2e-10 / 8.0e-10 | 1e-8 |
| `J_r(a)` vs Kalman3D_Wave_TFG's series implementation | 30 | 1.8e-15 | 1e-12 |
| Transported covariance vs Monte Carlo of `Log_G(Exp_G(-a) Exp_G(a+d))`, correlated `P` with attitude-translation cross terms (20 000 draws) | 1 | 0.023 (corr. units; untransported 0.57) | 0.04 |
| Accelerometer and magnetometer `H` vs FD of the residual: identity, small (3/-2 deg) and large (40/-30 deg) roll/pitch, 170 deg yaw, 2 m CARMA states, nonzero biases | 10 | 1.7e-10 (relative) | 1e-8 |
| Complete update (gain, Joseph, `J_r`, `Exp_G` injection) vs a dense reference with FD `H` and FD `J_r`, active and held bias, accel then mag | 20 | state 2.4e-11, covariance 1.6e-9 (corr.) | 1e-7 / 1e-6 |
| `Phi` vs FD of the nonlinear error flow: the 5 states above x h = 5, 50 ms x rates 0, 0.4, 1.3 rad/s | 30 | 3.0e-9 | 2e-8 |
| `Phi` vs RK4 of the continuous error model | 30 | 3.6e-15 | 1e-10 |
| `Qd` vs RK4 Lyapunov integration of the continuous model | 30 | 3.7e-11 (corr.) | 1e-7 |
| `Phi_{rho,phi} = 0`; `Phi_{rho_i,bg} = [w_i(h)]x Phi_{phi,bg}` (collapse identity) | 1 | 0 / 0 | 1e-15 |
| Collapsed `Qd` vs Kalman3D_Wave_TFG-style structured GL-5 quadrature | 1 | 8.3e-16 | 1e-9 |
| Semigroup on non-uniform steps (4 ms + 7 ms): `Phi`, `Qd` | 1 | 4.4e-16 / 8.8e-16 | 1e-12 / 1e-8 |
| Structured `time_update` vs dense `Phi P Phi^T + Qd` | 2 steps | 2.7e-16 | 1e-12 |
| Monte Carlo of the nonlinear stochastic truth (4000 particles, 1 s, inflated gyro/bias noise) vs propagated `P` | 1 | 0.058 (corr.); worst variance ratio 5.8 %; mean 2.4 SE | 0.08 / 10 % / 4.5 |
| Stationary prior `Sigma_s (x) P_inf` invariant over 100 s of prediction, double / float | 2 | 3.8e-14 / 1.7e-5 | 1e-9 / 2e-4 |
| Three-axis prior through the filter's discretization, 5.8 h at 20 ms: sample covariance vs `Sigma_s (x) P_inf` (with `r = 0.3`) | 1 | 0.023 (corr.) | 0.05 |
| Welch PSD of `p_x, p_y, p_z` vs analytic (0.12-0.9 Hz): mean / max dB; `a_z` max dB | 1 | 0.19 / 0.51; 1.14 | 0.35 / 1.2; 1.2 |
| x-y co-spectrum / `(s_x s_y S_p)` - `r`, mean over frequency | 1 | 0.016 | 0.04 |
| DC zero: `S_p(1 mHz)/S_p(2 mHz) = 1/4`; largest pole real part | 1 | 2.4e-5; -0.168 rad/s | 1e-4; < 0 |
| Float vs double over 60 s with +-25 % dt jitter: attitude / position / covariance | 1 | 3.5e-5 rad / 0.085 mm / 1.5e-4 (corr.) | 1e-4 / 1 mm / 0.02 |
| CARMA transition commutes with any rotation; `Qd(Q S Q^T) = (I (x) Q) Qd(S) (I (x) Q)^T` | 10 rotations | 2.2e-16 / 6.6e-17 | 1e-14 / 1e-13 |
| Complete filter, world frame rotated by yaw (isotropic / anisotropic `Sigma_s` rotated with it) | 2 x 20 s | 1.3e-12 / 1.5e-12 | 1e-9 |
| Anisotropic `A`: derived `Phi_{rho,phi}` vs FD of the nonlinear error flow | 1 | 3.4e-10 (relative) | 1e-7 |
| One prediction from equal seeds: `T^-1 P_TFG T^-T` vs CARMA3's `P` | 1 | 6.5e-16 (corr.) | 1e-9 |
| Heading information `e^T P^-1 e` never rises, no reset transport, 120 s accelerometer only | 24 000 steps | 0 | 1e-9 |
| Held accelerometer bias: unchanged, zero cross-covariance after 2000 updates | 1 | 0 / 0 | 1e-15 |
| CARMA3's noiseless three-axis sinusoid (its own 10 % bar); rotating version, TFG/MEKF error ratio | 2 | 8.6 %; 4.3e-5 | 10 %; 1 % |
| Physical / native attitude write keeps its covariance | 2 | 1.2e-15 / 8.7e-17 | 1e-12 / 1e-15 |
| Wrapper on a still, noisy IMU (as carma3-test): live, finite, `|p|` after 240 s | 1 | 41 mm | 0.5 m |

Reported quantities (no pass/fail bar), the mechanism evidence of
section 8:

| Quantity | Value |
|---|---|
| Linearization error of one step, `|xi(h) - Phi xi(0)| / |xi|`, `|xi|` ~ 1e-4 / 1e-3 / 1e-2 / 1e-1 | 1.8e-7 / 1.5e-6 / 1.2e-5 / 1.5e-4 |
| First accelerometer update after a seed from identical states: state difference TFG vs MEKF | 5.4e-4 (rad, m) |
| First update: covariance difference (correlation units) | 0.14 |
| 60 s side by side, max attitude difference: full reset / attitude-only / none | 0.011 / 0.017 / 0.018 rad |
| 60 s side by side, max position difference: full / attitude-only / none | 6.8 / 10.2 / 10.4 mm |
| Heading information over 120 s accelerometer-only (final/initial): TFG full reset / attitude-only / none / CARMA3 MEKF | 273 / 15 / 0.99 / 751 |
| Heading std after 120 s accelerometer-only, from 28.6 deg: TFG full / none / CARMA3 MEKF | 2.6 / 29.8 / 1.5 deg |

### 5.1 Failures during the verification (A-E classification)

Following AGENTS.md's discipline for failed attempts (this study is not
OU-III proof work; the classes are used for clarity):

| Failure | Class | Cause and resolution |
|---|---|---|
| Anisotropic-dynamics check passed at `-O2` with the wrong commutator sign, failed at `-O3` with the right one | E | the test's RK4 held `auto` Eigen expressions referring to destroyed temporaries (undefined behaviour); fixed with explicit types; the derived sign agrees with FD to 3e-10 |
| Quadrature comparison reported 8e278 | E | division by a zeroed diagonal in the comparison; compared with the full diagonal instead |
| Rotating sinusoid: TFG 12.3 % vs MEKF 10.4 % | E | the two filters had different gyro-bias random walks in the test; identically configured they agree to 4e-5 |
| Benchmark CSV had 7 rows of a default run | E | `--header --extended` without `--scenario` ran the default scenario; header-only now, grid restarted from an empty CSV |

None is a filter or a mathematical failure.

## 6. Simulation methodology

`tests/kalman_carma/carma3_bench --extended` (one scenario and seed per
process): a random-phase directional sea through the vessel RAOs of
`docs/carma-three-axis.md` section 2, one IMU/magnetometer stream per seed
(200 Hz; accelerometer 1.5 mg white, turn-on bias up to 8 mg with a random
walk; gyro 0.0016 rad/s white, bias up to 0.05 deg/s; magnetometer 0.8 uT at
25 Hz with hard/soft-iron terms) fed identically to every filter.  Errors
are against the wave coordinate in the geographic frame (estimates rotated
by the declination); scored from 300 s; startup = the first 120 s after
Live.  No truth enters any estimator.  New scenario modifiers (every filter
sees the same stream): `turns` (90 deg heading changes in 6 s every 120 s,
peak rate 0.49 rad/s), `gbdrift` (gyro-bias ramps of 0.03-0.05 deg/s over
1200 s plus a 600 s oscillation), `temp` (25 -> 45 C warm-up with a 300 s
cycle, true accelerometer-bias temperature coefficients 0.0015-0.003
m/s^2/C against the filters' 0.002), `maggap<N>` (magnetometer absent for
N s from T/3).

**Allocation.**  100 paired seeds on the four principal seas (JONSWAP
0.29 Hz 0.8 m, PM 0.25 Hz 1 m, mixed swell + wind, changing f_p); 30 on
each of 29 further conditions (broadband, bimodal, abrupt change, two more
sea states, five headings, anchored, drifting, underway head and following,
manoeuvring, still IMU with and without calibrated bias, 5 cm and 2 cm
seas, large roll/pitch, heel/trim, rapid turns, 300 s and 900 s
magnetometer gaps, gyro-bias drift, accelerometer-bias drift, temperature
drift, +-25 % and +-1 % sample jitter); 20 seeds of 3600 s on JONSWAP, the
mixed sea and a 2400 s magnetometer loss.  Thirty seeds suffice for the
paired contrast: the paired TFG-MEKF differences have a seed-to-seed spread
one to three orders of magnitude below the RMS levels (section 7), so the
95 % interval of a paired change is far narrower than the +-1 % margin.
Filters per run: OU3, A41, T41 and the ablations T41R0, T41RA, T41NW (the
3600 s runs: OU3, A41, T41, T41R0).  The CARMA(5,2) follow-up (A52, T52,
T52R0 with A41/T41) is a separate grid.  The independent pinned
28 ft vessel-RAO dataset (8 records x 10 seeds) is replayed through
`carma3-sim` (A41, T41, T41R0) and `kalman_ou_iii-sim` (OU3).  Post hoc,
designed after the principal results and labelled as such (section 8,
summary section 11): 30 seeds each of the mixed sea with a 300 s
magnetometer gap and of JONSWAP with a calibrated accelerometer (OU3, A41,
T41, T41R0, T41RA), and the realized motion levels of every condition
(`--level truth`, 3 seeds, `truth_levels.csv`).

**Statistics.**  Paired change `1 - mean(T41)/mean(A41)` over the same
seeds with a 95 % bootstrap interval (10 000 resamples), win rate, two-sided
Wilcoxon signed-rank p with Holm correction within each table's family, and
an equivalence statement when the whole interval lies inside +-1 % (margin
fixed before the runs).  Significance and equivalence are reported
independently.

**Metrics.**  Per axis: RMS, mean, peak, 60 s running-mean RMS (low-frequency
drift), linear trend, velocity RMS, 2-dof `[p, v]` NEES, predicted std,
convergence time, startup RMS and peak, last-60 s RMS; horizontal and 3-D
RMS; position 3-dof and `[p, v]` 6-dof NEES; accelerometer NIS (mean and
fraction in the 95 % chi-square band); roll/pitch/heading RMS.  With
`--extended`: attitude NEES (3 dof) and heading NEES (1 dof) of
`phi_w = Log(Rhat R^T)` in the common world-referred coordinates; the joint
9-dof NEES of `[phi_w, phat - p, vhat - v]` in common coordinates and in each
filter's own coordinates; 95 % coverage per axis; tilt and heading-error
RMS/maximum; gyro- and accelerometer-bias RMS (scored window and last
60 s); maximum correlation-scaled condition number of `P`; magnetometer
samples forwarded to and applied by the estimator; gap metrics (heading RMS
and NEES in the gap, predicted heading std at its end, RMS after it,
recovery time); after each rapid turn, the peak attitude error and the
recovery time.

## 7. Results

All numbers are from `reports/results/tfg_carma3/summary.md` (generated;
section numbers below refer to it) and `paired_stats.json`.  Change =
`1 - mean(TFG)/mean(Q-MEKF)` over paired seeds, positive when TFG-CARMA3 is
lower; `[ ]` the 95 % bootstrap interval; `equiv.` the interval inside
+-1 %; `p` the Holm-corrected Wilcoxon p where below 0.05.

### 7.1 Principal seas (100 paired seeds each)

| metric | JONSWAP 0.29 Hz 0.8 m | PM 0.25 Hz 1 m | mixed swell + wind | changing f_p |
|---|---|---|---|---|
| heave RMS (m), Q-MEKF / TFG | 0.01461 / 0.01461 | 0.02104 / 0.02104 | 0.1491 / 0.1491 | 0.01401 / 0.01401 |
| heave change | +0.001 % equiv. | +0.001 % equiv. | +0.001 % equiv. | +0.002 % equiv. |
| horizontal RMS (m), Q-MEKF / TFG | 0.02127 / 0.02127 | 0.0291 / 0.0291 | 0.1939 / 0.1937 | 0.02017 / 0.02017 |
| horizontal change | -0.003 % equiv. | -0.000 % equiv. | **+0.091 % [+0.084, +0.098] equiv., p = 3e-16** (TFG lower in 100/100) | +0.002 % equiv. |
| sway change | -0.013 % equiv., p = 0.004 | -0.009 % equiv. | +0.44 % [+0.33, +0.54] equiv., p = 5e-9 | -0.006 % equiv. |
| roll / pitch RMS change | +0.15 % / +0.017 % equiv. | +0.10 % / +0.031 % equiv. | +0.004 % [-0.96, +0.94] equiv. / **+4.66 % [+2.59, +6.52], p = 0.003** | +0.11 % / +0.041 % equiv. |
| yaw RMS change | +0.17 % [+0.067, +0.28] equiv. | +0.15 % equiv. | +1.01 % [+0.38, +1.65] | +0.077 % equiv. |
| max heading excursion change | +0.36 % equiv., p = 0.02 | +0.34 % equiv. | **+2.77 % [+1.89, +3.65], p = 3e-6** | +0.34 % equiv. |
| gyro-bias RMS change | +0.15 % equiv. | +0.11 % equiv. | **+1.17 % [+0.82, +1.51], p = 2e-7** | +0.57 % equiv., p = 0.002 |
| accel-bias RMS change | +0.091 % equiv., p = 0.03 | +0.074 % equiv., p = 0.04 | **+3.17 % [+1.89, +4.34], p = 0.001** | +0.10 % equiv. |
| gyro-bias, last 60 s, change | **-0.79 % [-1.29, -0.29]** | **-0.78 % [-1.29, -0.28]** | -0.11 % equiv. | -0.29 % equiv. |
| accel-bias, last 60 s, change | +0.22 % equiv., p = 0.01 | +0.18 % equiv., p = 0.01 | **+7.47 % [+5.14, +9.57], p = 2e-5** | +0.18 % equiv., p = 0.007 |

* **Translation is equivalent** in every principal sea: every heave,
  surge, sway, horizontal, 3-D, velocity, low-frequency (60 s mean) and
  peak-error change lies inside +-1 % with its whole interval; most are
  below 0.02 %.  Several are statistically significant (the paired design
  resolves 0.01 %), none is practically relevant.  The largest is the
  mixed sea's sway, +0.44 %.
* **The mixed sea separates the pair** in attitude and bias, not in
  translation: pitch 4.7 %, maximum heading excursion 2.8 %, gyro bias
  1.2 %, accelerometer bias 3.2 % (7.5 % at the end of the run) lower with
  TFG-CARMA3, all significant after correction and broad across seeds (the
  median per-seed change is at least the mean; the five most favourable
  seeds carry about a fifth of the total).  Roll is unchanged on average
  with a wide interval.  Section 8 locates the cause.
* **One unfavourable signal:** the gyro-bias error over the last 60 s of
  the JONSWAP and PM runs is 0.8 % higher with TFG-CARMA3 (about
  0.17 deg/h on 21 deg/h; raw Wilcoxon p = 0.002 in both, Holm p = 0.12;
  TFG lower in 35 % of seeds).  Its interval crosses -1 %, so it is neither
  shown equivalent nor significant after correction.  The same filters'
  gyro-bias RMS over the whole scored window is 0.1-0.2 % *lower* with
  TFG-CARMA3, so the signal concerns the end-of-run estimate of the weakly
  observable bias components; it is not explained (section 12).
* **Accelerometer bias and tilt are one quantity here.**  In every condition
  the accelerometer-bias RMS error is 17.1-17.4 mg per degree of tilt RMS
  error, against `g pi/180 = 17.45 mg/deg`: the residual horizontal bias
  error is the tilt-equivalent bias that the accelerometer cannot separate
  from tilt
  without attitude excitation, so the tilt and bias rows move together
  (sections 2 and 2b of the summary).

### 7.2 Every condition

Summary sections 2 and 2b (one row per condition, Holm-corrected over the
whole matrix) and `paired_change_by_condition.svg`:

* **Translation:** equivalent (whole interval inside +-1 %) in 34 of the
  36 conditions, typically within +-0.01 %.  The two exceptions are the
  magnetometer outages of 900 s and 2400 s (7.4), where the heading error
  of both filters grows to several degrees and the horizontal intervals
  widen to -1.6..+0.4 % and -17.8..+0.8 % (point estimates unfavourable
  to TFG-CARMA3, neither significant).
* **Attitude and bias:** equivalent wherever the CARMA prior fits (tilt
  changes between -0.33 % and +0.44 %, mostly +0.1..+0.2 %); clearly better
  where it does not:
  mixed sea (tilt +3.0 %, 3600 s +6.0 %), manoeuvring (low-frequency
  translation outside the prior: yaw +13.3 % [+7.5, +18.9] (raw p = 7e-4,
  not significant after the matrix-wide correction), tilt +2.6 %
  (p = 2e-6), gyro bias +6.9 %), anchored (yaw +0.7 %), bimodal (yaw
  +0.7 %).  The
  gyro-bias drift case gains +0.9 % in yaw (p = 3e-7).
* **CARMA3's known weaknesses are untouched.**  On the still IMU the CARMA
  filters' horizontal error is 0.084 m against OU-III's 0.053 m and their
  heave 0.053 m against 0.027 m; in the mixed sea their heave is 0.149 m
  against OU-III's 0.132 m.  TFG-CARMA3 changes none of these by more than
  0.02 %: they are properties of the CARMA prior and the wrapper's mapping,
  not of the error geometry.
* Unfavourable points worth naming (none significant after correction):
  gyro-bias RMS in bow-quartering seas -5.8 % [-10.6, -1.2] (raw p = 0.04),
  tilt with head seas underway -0.33 % and in the abrupt change -0.30 %
  (both inside the margin).

### 7.3 Consistency

| | JONSWAP | PM | mixed | changing f_p | still IMU |
|---|---|---|---|---|---|
| position NEES (3): OU-III / Q-MEKF / TFG | 2.71 / 1.76 / 1.76 | 3.21 / 2.12 / 2.12 | 19.5 / 44.3 / 44.3 | 2.79 / 1.72 / 1.72 | |
| attitude NEES (3): OU-III / Q-MEKF / TFG / TFG no reset | 50.7 / 48.0 / 47.5 / 47.3 | 49.4 / 46.6 / 46.2 / 46.1 | 58.2 / 50.3 / 46.5 / 44.4 | 54.6 / 51.6 / 51.0 / 50.8 | 80.6 / 74.4 / 74.0 / 73.8 |
| yaw NEES (1): Q-MEKF / TFG / TFG no reset | 12.9 / 12.5 / 12.4 | 11.3 / 11.0 / 10.9 | 12.4 / 11.2 / 10.8 | 17.1 / 16.6 / 16.4 | 39.8 / 39.3 / 39.1 |
| attitude NEES in its 95 % band (ideal 0.95): Q-MEKF / TFG | 0.088 / 0.090 | 0.093 / 0.095 | 0.097 / 0.104 | 0.065 / 0.066 | |
| surge / heave 95 % coverage, Q-MEKF = TFG | 0.99 / 0.98 | 0.98 / 0.96 | 0.35 / 0.36 | 0.99 / 0.98 | |
| accelerometer NIS mean (3), Q-MEKF = TFG | 0.089 | 0.077 | 0.19 | 0.094 | |
| largest corr.-scaled cond(P): OU-III / Q-MEKF / TFG | 8600 / 4210 / 2330 | 9990 / 4910 / 2370 | 16 600 / 9040 / 4160 | 10 600 / 4540 / 3820 | |

* Translational consistency is identical (to 0.1 %) and inherited: slightly
  conservative position covariances in the regular seas, badly
  overconfident in the mixed sea (NEES 44, coverage 36 %), where the CARMA
  prior misfits.
* Attitude is overconfident by a factor of about 15 in every filter,
  including OU-III (section 8.4).  TFG-CARMA3 reduces attitude NEES by 1 %
  in regular seas and 8 % in the mixed sea, yaw NEES by 3 % and 10 %; the
  no-transport variant by a little more.  A real but small consistency
  improvement on top of a large shared inconsistency.
* The accelerometer NIS is far below its 3 degrees of freedom for every
  filter (the wrapper's accelerometer noise, with its vibration inflation,
  exceeds the simulated noise); it is the same for both error geometries
  to four digits, as first-order equivalence requires.

### 7.4 Magnetometer outage and reacquisition

| | 300 s gap | 900 s gap | lost at 1200 s (2400 s) |
|---|---|---|---|
| yaw RMS in the gap (deg): OU-III / Q-MEKF / TFG / TFG no reset | 3.69 / 2.13 / 2.13 / 2.10 | 22.5 / 5.17 / 5.45 / 4.60 | 58.7 / 16.5 / 18.4 / 17.3 |
| TFG vs Q-MEKF | -0.29 % [-3.1, +2.0] | -5.4 % [-14.4, +1.9] | -11.1 % [-28.6, +1.8] |
| yaw NEES in the gap (1): OU-III / Q-MEKF / TFG / TFG no reset | 13.7 / 5.3 / 4.3 / 3.9 | 395 / 16.4 / 7.3 / 2.2 | 2570 / 141 / 57 / 1.2 |
| predicted yaw std at the gap end (deg): same order | 1.16 / 1.65 / 2.31 / 2.89 | 1.16 / 1.71 / 2.83 / 10.8 | 1.16 / 1.72 / 2.84 / 41.5 |
| yaw RMS 60 s after the gap (deg) | 3.42 / 1.95 / 1.91 / 1.88 | 19.2 / 3.67 / 3.04 / 2.07 | -- |
| TFG vs Q-MEKF | +1.9 % [-2.4, +6.3] | **+17.2 % [+8.9, +23.6]** | -- |
| heading recovery (s) | 70 / 38 / 36 / 36 | 66 / 40 / 34 / 23 | -- |
| horizontal RMS (m) | 0.042 / 0.0217 / 0.0217 / 0.0217 | 0.060 / 0.0239 / 0.0241 / 0.0234 | 0.121 / 0.0382 / 0.0411 / 0.0402 |
| magnetometer samples forwarded = applied (every filter) | 21 700 | 21 700 | 29 200 |

* Without the magnetometer the heading is unobservable and its error is a
  random walk driven by the vertical gyro-bias error.  Only the TFG without
  reset transport reports that honestly: its heading NEES stays near 1-4
  and its predicted heading uncertainty grows to 41 deg over 2400 s.  The
  full-reset TFG is 1.2-2.5 times less overconfident than CARMA3, OU-III up
  to 2570 times overconfident (it predicts 1.2 deg at the end of a gap with
  a 59 deg error).
* The payoff of the honest covariance is reacquisition: after the 900 s
  gap TFG-CARMA3 returns to the pre-gap heading accuracy 6 s sooner and
  with 17 % less post-gap heading error, the no-transport variant 17 s
  sooner with 44 % less, because the first magnetometer updates are given
  the weight the actual error deserves.
* Accuracy *during* a long gap is not better: the point estimates of
  heading and horizontal error in the 900 s and 2400 s gaps favour CARMA3
  (-5 % and -11 % heading), with intervals that include zero at 30 and 20
  seeds.  An honest covariance does not make an unobservable state more
  accurate.
* Every magnetometer sample forwarded by the wrapper was applied by every
  filter (forwarded = applied in every run); there is no gating, so no
  rejected update.

### 7.5 Startup and rapid rotation

Equivalent: heave convergence time, startup heave/surge RMS and peak, and
startup tilt agree to better than 0.2 % in every principal sea; startup
heading is 0.1-0.8 % better with TFG-CARMA3 (mixed sea 2.99 against 3.02
deg).  Surge convergence differs by up to 2 % (JONSWAP 674 s against
687 s).  With 90-degree turns in 6 s (peak rate 0.49 rad/s) the peak
attitude error after a turn is 2.25 deg for both
and neither leaves its recovery band; the large roll/pitch scenario is
equivalent in translation and tilt and +0.45 % in heading.  Rapid
rotations expose neither an advantage nor a weakness: the invariant error
dynamics are exact in the rotation, and so is the MEKF's.

### 7.6 Independent vessel-RAO dataset

The pinned 28 ft vessel-RAO records (JONSWAP and Pierson-Moskowitz/Stokes,
H_s 0.27-8.5 m; 8 records x 10 seeds, `carma3-sim` for A41/T41/T41R0 and
`kalman_ou_iii-sim` for OU-III; summary section 8) are a second simulator
with its own sensor models.  Displacement errors in % of H_s:

| record | heave: OU-III / Q-MEKF / TFG | surge | sway | pitch RMS (deg): Q-MEKF / TFG / TFG no reset |
|---|---|---|---|---|
| JONSWAP 0.27 m | 4.21 / 3.12 / 3.12 | 3.20 / 1.21 / 1.21 | 3.31 / 0.878 / 0.878 | 0.259 / 0.257 / 0.255 |
| JONSWAP 1.5 m | 4.03 / 2.21 / 2.21 | 4.49 / 2.38 / 2.38 | 4.30 / 2.13 / 2.13 | 0.236 / 0.236 / 0.235 |
| JONSWAP 4.0 m | 3.93 / 2.17 / 2.17 | 5.60 / 2.62 / 2.63 | 4.47 / 2.41 / 2.42 | 0.239 / 0.239 / 0.239 |
| JONSWAP 8.5 m | 3.81 / 2.12 / 2.12 | 5.85 / 2.82 / 2.82 | 5.03 / 2.48 / 2.48 | 0.230 / 0.230 / 0.229 |
| PM/Stokes 0.27 m | 4.19 / 3.25 / 3.25 | 3.23 / 1.19 / 1.19 | 3.29 / 0.852 / 0.851 | 0.257 / 0.256 / 0.254 |
| PM/Stokes 1.5 m | 3.96 / 2.33 / 2.33 | 4.36 / 2.36 / 2.36 | 4.07 / 2.07 / 2.07 | 0.236 / 0.235 / 0.234 |
| PM/Stokes 4.0 m | 3.91 / 2.44 / 2.44 | 5.46 / 2.75 / 2.75 | 4.34 / 2.43 / 2.43 | 0.239 / 0.239 / 0.239 |
| PM/Stokes 8.5 m | 3.81 / 2.38 / 2.38 | 5.88 / 3.06 / 3.06 | 4.87 / 2.60 / 2.60 | 0.230 / 0.229 / 0.229 |

Over the 80 paired replays the translational changes are +0.006 %
(heave), -0.012 % (surge) and -0.018 % (sway), with no single replay
differing by more than 0.11 %; roll, pitch and heading change by +0.06 %,
+0.29 % and +0.25 % on average (single replays up to 1.3 %, 4.4 % and
5.5 %).  The dataset reproduces the synthetic result: no translational
difference, a small attitude difference in TFG-CARMA3's favour, and
CARMA3's advantage over OU-III unchanged.

## 8. Mechanism

### 8.1 Why translation does not change

Section 2.8 is the explanation, and the data follow it.  Both filters
linearize the same system exactly to first order, with the same
measurement information, the same CARMA prior and the same wrapper, so
their linear models are similar under `T(Xhat)` and their gains produce the
same first-order correction.  What remains is second order: the injection
(`Exp_G(-a) Xhat` against CARMA3's multiplicative/additive update) and the
coordinates in which the covariance is carried across a correction, whose
leading term in CARMA3's coordinates is `(1/2)[Delta what_i]x phi`, a
product of a state correction and the attitude error.  In the seas the
prior fits, both factors are small after convergence, so the filters agree
to 1e-4 in every translational RMS and to 2e-3 in attitude.  Against the
specific questions of the study:

* **Attitude-translation coupling.**  TFG-CARMA3 carries it in the
  dynamics (`Phi_{rho_i,bg} = [w_i(h)]x F` and the matching noise blocks)
  with a state-independent `H`; CARMA3 carries it in its accelerometer
  Jacobian `[Rhat^T (ahat - g)]x`.  To first order this is one coupling in
  two coordinate systems, and it produces no translational difference
  (below 0.02 % in the regular seas, below 0.15 % in the mixed sea).
* **Physically consistent CARMA corrections.**  The group injection rotates
  every world column by the attitude correction (`W+ = Exp(-a_phi) W +
  ...`); with per-update attitude corrections of order 1e-4 rad and world
  states below 1 m this is below 0.1 mm per update, and CARMA3's covariance
  carries the same term to first order (the `[what_i]x Rhat dtheta` part of
  `T`).  No accuracy effect is visible.
* **Gyro-bias error into translational uncertainty.**  Present in both;
  `[p, v]` NEES changes by 0.3-0.5 % (JONSWAP 5.44 -> 5.41), position
  NEES by less than 0.1 %.
* **Geometry or tuning.**  Tuning is identical by construction (one
  wrapper source, mapping and clamps; the committed operating point is
  driven by the shared front end).  The only horizontal differences above
  0.05 % (mixed sea: surge +0.07 %, sway +0.44 %) come with the attitude
  difference of 8.2 and vanish where it vanishes.
* **Low-frequency drift and startup.**  60 s running-mean errors change by
  at most 0.1 %, startup metrics by at most 0.8 % (heading) (7.5).

### 8.2 Where the pair separates

`mechanism_error_level.svg` plots, per condition, the paired tilt change
against CARMA3's own 3-D velocity error (left; the benchmark does not
record the acceleration-state error, and velocity is the closest recorded
state) and against the true horizontal wave acceleration (right).  All 32
conditions with a CARMA3 velocity error below 0.045 m/s show TFG-CARMA3 tilt changes
within +-0.5 %.  The four above 0.06 m/s are the mixed sea (1200 s and
3600 s: +3.0 % and +6.0 %), manoeuvring (+2.6 %) and the 2400 s
magnetometer loss, whose velocity error comes from the heading error
rotating the horizontal estimates and whose tilt is unchanged (+0.06 %).
Position error does not order the conditions (the still IMU and the 2.5 m
PM sea have CARMA3 position errors as large as manoeuvring and no change):
their errors are slow drift, not a misfit of the acceleration.  The motion
amplitude does not order them either: the mixed sea's horizontal
acceleration (0.28/0.47 m/s^2 RMS surge/sway) is no larger than JONSWAP's
(0.41/0.37 m/s^2).  What the three conditions share is that the
single-peak CARMA(4,1) prior misfits them -- a bimodal sea, and a
low-frequency translation the prior does not contain -- so the CARMA
acceleration state carries large errors and receives large corrections
every update.  That is the regime in which the second-order terms of 8.1
are not small.

### 8.3 Which ingredient produces it

The ablations separate the TFG's ingredients (change against CARMA3
Q-MEKF; NEES as means):

| | mixed: tilt RMS | mixed: accel-bias RMS | mixed: attitude NEES (Q-MEKF 50.3) | mixed: yaw NEES (12.4) | mixed: predicted yaw std (0.530 deg) | manoeuvring: yaw RMS | manoeuvring: tilt RMS |
|---|---|---|---|---|---|---|---|
| TFG-CARMA3 | +3.04 % | +3.17 % | 46.5 | 11.2 | 0.555 deg | +13.3 % | +2.56 % |
| no reset transport (T41R0) | +4.63 % | +4.77 % | 44.4 | 10.8 | 0.564 deg | +23.6 % | -0.20 % |
| attitude-only reset (T41RA) | +4.68 % | +4.83 % | 44.4 | 10.8 | 0.563 deg | +20.8 % | -0.83 % |
| invariant-coordinate writes (T41NW) | +3.06 % | +3.18 % | 46.5 | 11.1 | 0.555 deg | +13.3 % | +2.56 % |

* The covariance-write coordinates do not matter (T41NW = T41): the gain is
  not an artefact of how CARMA3's phase-change policies were mapped.
* The invariant error coordinates produce the gain; the full reset
  transport gives back about a third of it in the mixed sea and 44 % of
  the heading gain when manoeuvring (where it keeps a tilt gain the
  other variants do not have).  The transport is the correct first-order
  covariance map across the injection (section 5: 0.023 against 0.57
  untransported, correlation units, at a large correction), so removing it
  is not "more accurate" in the local sense.
* **The CARMA(5,2) follow-up amplifies the same ordering** (section 10):
  in the mixed sea that prior misfits even more (horizontal error 0.40 m),
  and the quaternion MEKF's tilt error grows to 1.64 deg (attitude NEES
  819), the full-reset TFG's to 0.99 deg (288), while the TFG without
  transport keeps 0.34 deg (NEES 50) -- the attitude it reaches with the
  better-fitting CARMA(4,1) prior, lower than the MEKF's in 30 of 30 seeds.

The interpretation consistent with all of this has two faces of one
property, the estimate dependence of the error coordinates:

* *Observable directions (tilt).*  CARMA3's accelerometer attitude Jacobian
  `[Rhat^T (ahat - g)]x` is evaluated at the CARMA acceleration estimate.
  When the prior misfits, `ahat` is wrong by several tenths of m/s^2, the
  MEKF linearizes the attitude about the wrong specific force, and the
  first-order equivalence of 2.8 holds only to `O(|da| |dtheta|)`.  In the
  invariant coordinates the accelerometer residual is exactly linear in the
  acceleration error (`H_a` does not depend on the estimate, 2.6), so the
  misfit cannot enter the attitude sensitivity.  The reset transport's
  attitude column `Q_r(phi, a_rho)` is built from the *world-state
  corrections* -- large exactly when the prior misfits -- and re-introduces
  a correction-dependent attitude-translation coupling.
* *The unobservable direction (heading).*  In CARMA3's coordinates the
  heading direction `n(Xhat) = [Rhat^T e_z, 0, [what]x e_z, 0]` contains
  the world states, so large CARMA corrections move it and each
  accelerometer update adds false heading information; the transport
  re-introduces the same dependence in the TFG (2.9).  The measured leakage
  over 120 s without magnetometer orders the variants as the results do:
  CARMA3 751x, full-reset TFG 273x, attitude-only 15x, none 0.99x (5).
  Every TFG variant predicts a larger heading uncertainty (0.555-0.564 deg
  against 0.530 deg in the mixed sea) while making a smaller heading error
  (1.535-1.538 against 1.551 deg): less overconfident and more accurate.

The post-hoc gap test (summary section 11; mixed sea, magnetometer absent
for 300 s, 30 seeds) is consistent with this: with heading uncorrectable,
the gains grow -- tilt +5.1 % (no transport +8.1 %) and accelerometer bias
+5.4 % (+8.4 %), both significant after correction; pitch +7.3 %
(+12.4 %) and heading after the gap +12.5 % (+19.8 %), with intervals
that exclude zero but do not survive the correction -- and the heading
NEES in the gap is 10.7 for CARMA3, 6.3 and 4.9 for the TFG variants.  It is not a proof: no experiment removed
the estimate dependence from CARMA3's MEKF itself (section 12).

### 8.4 Attitude consistency is mostly a shared accelerometer-bias effect

Every filter's attitude NEES is about 15 times its 3 degrees of freedom in
the principal seas (OU-III 49-58, CARMA3 47-52, TFG-CARMA3 46-51), with
7-10 % of samples inside the 95 % band.  TFG-CARMA3 lowers it by 1 %
(regular seas) to 8 % (mixed), far from removing it.  The cause is shared
input, not error geometry:

* the simulated accelerometer has a turn-on bias of up to 8 mg and a bias
  random walk of 5e-4 m/s^2/sqrt(s) (1.8 mg after 1200 s); the filters'
  accelerometer-bias model (OU-III's, inherited by both CARMA filters) has
  a 0.4 mg stationary standard deviation and an implied walk of
  8e-5 m/s^2/sqrt(s);
* the tilt error is the tilt-equivalent bias error (7.1), so the predicted
  tilt uncertainty is set by the 0.4 mg prior while the actual one is set
  by several mg;
* with a calibrated accelerometer (post-hoc test, JONSWAP, +-1 mg turn-on,
  the same walk) the tilt RMS falls from 0.309 to 0.124 deg and CARMA3's
  attitude NEES from 48 to 23, of which heading accounts for about 17; the
  tilt part falls from about 35 to about 6 (2 dof).  The still IMU shows
  the same (74 -> 48).

The remaining heading overconfidence (yaw NEES 11-17 with the magnetometer,
in every filter) is not attributed by this study; the simulated
magnetometer's hard/soft-iron and disturbance terms against the filters'
white-noise magnetometer model are the obvious candidate.

### 8.5 Numerical conditioning is not estimation accuracy

The largest correlation-scaled condition number of `P` is 1.2-2.2 times
smaller in TFG-CARMA3 (JONSWAP 2330 against 4210; mixed 4160 against 9040;
OU-III 8600-16 600).  This is a property of the coordinates: the invariant
`rho` absorbs the `[what]x Rhat dtheta` part of CARMA3's translational
error, so the attitude-translation correlations are weaker.  It costs
nothing and helps float arithmetic, but it is not an accuracy gain -- the
estimates agree to 1e-4 where the condition numbers differ by a factor of
two -- and it is reported separately from consistency and accuracy.

## 9. Figures

All figures are generated by `python3 tools/tfg_carma3_benchmark.py report`
from the CSV files in `reports/results/tfg_carma3/` (SVG, fixed hash salt,
no timestamps).  Time series are seed 1 of the condition, recorded at
20 Hz.

| figure | content | reading |
|---|---|---|
| `timeseries_xyz_jonswap.svg` | true vs estimated heave, surge, sway, 600-660 s | the two CARMA filters overlap to line width; OU-III overshoots the surge and sway peaks |
| `errors_xyz_jonswap.svg`, `errors_xyz_mixed.svg`, `errors_xyz_stationary.svg` | translational errors over time | identical traces for Q-MEKF and TFG; the mixed sea's large wave-band surge error, the still IMU's slow horizontal oscillation and both CARMA filters' heave transient before 150 s on the still IMU (-7.5 m, before the scored window; OU-III -2 m) are shared |
| `attitude_errors_jonswap.svg`, `_mixed.svg`, `_stationary.svg`, `_turns.svg` | roll, pitch, heading error | Q-MEKF and TFG nearly overlap in a single seed (the mixed-sea pitch difference of 7.1 is about 0.01 deg RMS and shows in the paired statistics); OU-III's roll and heading offsets differ |
| `error_psd_jonswap.svg`, `_mixed.svg`, `_stationary.svg` | error spectra | identical in the wave band; no low-frequency difference |
| `covariance_vs_error.svg` | predicted 2-sigma vs actual heave/surge error | position covariances identical in both geometries |
| `consistency_nees_nis.svg` | NEES/NIS distributions over seeds (regular seas as boxes, mixed-sea median as a diamond) | attitude NEES medians near 40-50 for every filter (8.4), the TFG's slightly lower; position NEES below 3 in the regular seas and about 45 in the mixed sea for both CARMA filters; NIS far below 3 |
| `startup_jonswap.svg` | first 120 s after Live | equivalent |
| `rms_distribution.svg` | per-seed RMS distributions, principal seas | overlapping |
| `paired_change_by_condition.svg` | paired change with interval, every condition | all translation inside the margin; tilt outside only for the mixed sea and manoeuvring |
| `mechanism_error_level.svg` | tilt change against CARMA3's error level and against the true motion | section 8.2 |
| `heading_maggap300.svg`, `_maggap900.svg`, `_maggap2400.svg`, `heading_mixed_maggap300.svg` | heading error and predicted 2-sigma through a gap | with the magnetometer every filter's heading error lies outside its own band much of the time (the shared yaw NEES of 11-17); in the gaps all bands widen, the no-transport TFG's most; over 2400 s only its band follows the error (to 80 deg) while the others stay within a few degrees |
| `cost.svg` | multiply-adds, host time, object size | section 11 |

Stationary and mixed-sea figures are separate files, as these were the
weak conditions of the earlier CARMA work.

## 10. Follow-up: CARMA(5,2) on all three axes

Run only after the primary comparison and kept separate from it (summary
section 7; 30 paired seeds on the four principal seas, the still IMU, the
300 s gap and the bimodal sea; `A52` CARMA3 Q-MEKF, `T52` TFG, `T52R0` TFG
without reset transport, with A41/T41 in the same rows).

**Mathematics.**  Nothing new is needed: `TwoFrameGroup<T, 5, 2>` carries
five world columns `[p v a l1 l2]`; one common `A` on all axes keeps
`Phi_{rho,phi} = 0` and the collapse identity of 2.4; `NX = 24`.  Mixed
orders (e.g. CARMA(4,1) horizontal, CARMA(5,2) vertical) are a different
matter: the latent coordinates of the two blocks are not columns of one
world matrix, so they would not form world vectors and the transition
would acquire the attitude coupling of 2.10.  They were not built.

**Computation.**  24 states; 19 083 multiply-adds per sample (TFG), 12 009
without reset transport, against 5697 for the CARMA(5,2) MEKF (3.3x and
2.1x).

**Results.**

| | JONSWAP | PM | changing f_p | still IMU | 300 s gap | bimodal | mixed |
|---|---|---|---|---|---|---|---|
| horizontal RMS (m): A41 / A52 / T52 / T52R0 | 0.0217 / 0.0168 / 0.0168 / 0.0168 | 0.0295 / 0.0268 / 0.0268 / 0.0268 | 0.0207 / 0.0177 / 0.0177 / 0.0177 | 0.084 / 0.087 / 0.087 / 0.087 | 0.0217 / 0.0168 / 0.0168 / 0.0168 | 0.079 / 0.116 / 0.116 / 0.116 | 0.194 / 0.402 / 0.397 / 0.394 |
| heave RMS (m): A41 / A52 / T52 / T52R0 | 0.0147 / 0.0164 / same / same | 0.0210 / 0.0265 / same / same | 0.0140 / 0.0165 / same / same | 0.053 / 0.046 / same / same | 0.0147 / 0.0164 / same / same | 0.053 / 0.100 / same / same | 0.149 / 0.346 / 0.346 / 0.346 |
| tilt RMS (deg): A52 / T52 / T52R0 | 0.306 / 0.306 / 0.306 | 0.316 / 0.316 / 0.315 | 0.302 / 0.301 / 0.301 | 0.351 / 0.350 / 0.350 | 0.297 / 0.297 / 0.296 | 0.313 / 0.309 / 0.307 | **1.644 / 0.991 / 0.338** |
| attitude NEES: A52 / T52 / T52R0 | 48.8 / 48.5 / 48.3 | 48.5 / 48.1 / 47.9 | 54.2 / 53.6 / 53.4 | 75.8 / 75.0 / 74.7 | 40.9 / 40.4 / 40.2 | 46.4 / 43.8 / 42.6 | **819 / 288 / 50** |

* Between the error geometries the CARMA(5,2) picture is the CARMA(4,1)
  one: translation equivalent (within 0.01 % in the regular seas, T52
  1.0 % better horizontally in the mixed sea), attitude and consistency
  slightly better.
* The prior matters far more than the geometry: CARMA(5,2) lowers the
  horizontal error by 9-23 % in the regular seas and raises heave by
  12-26 %, and in the bimodal and mixed seas it raises both by 48-132 % --
  the trade-off `docs/carma-three-axis.md` reports, unchanged here.
* **Mixed sea: the CARMA(5,2) MEKF's attitude fails** (tilt 1.64 deg
  against 0.32 deg with CARMA(4,1), attitude NEES 819), the full-reset TFG
  reduces it (0.99 deg, 288) and the TFG without reset transport
  removes it (0.34 deg, 50; lower than the MEKF in 30 of 30 seeds).  This
  is the clearest single observation of the mechanism of 8.3: the worse the
  prior fits, the more the estimate-dependent error coordinates feed the
  CARMA misfit into the attitude.

**Is TFG-CARMA(5,2) justified?**  Mathematically yes (no new structure);
computationally it is the most expensive variant studied.  It does not
change the CARMA(5,2) prior's translational trade-off, which is what
decides between (4,1) and (5,2).  Its one distinctive result -- an attitude
immune to the prior's misfit without reset transport -- is a property of
the invariant coordinates that T41R0 already shows at lower cost.

## 11. Computational feasibility

| | OU-III (shipping) | CARMA3 Q-MEKF | TFG-CARMA3 | TFG, no reset transport | TFG, attitude-only reset |
|---|---|---|---|---|---|
| error states | 21 | 21 | 21 | 21 | 21 |
| multiply-adds, prediction | 5913 | 1731 | 6465 | 6465 | 6465 |
| multiply-adds, accelerometer update | 6570 | 2439 | 7644 | 2217 | 4890 |
| multiply-adds per sample (prediction + accelerometer) | 12 483 | 4170 | 14 109 | 8682 | 11 355 |
| host update, fixed dt (us) | 22.31 | 7.42 | 11.29 | 8.08 | 9.64 |
| host update, +-25 % jitter (us) | 23.97 | 9.62 | 12.76 | 8.87 | 10.50 |
| filter object, x86-64 (B) | 5728 | 6128 | 5888 | as T41 | as T41 |
| code, object-level, x86-64 `-Os` (B) | 234 765 | 192 102 | 167 417 | as T41 | as T41 |
| largest stack frame in the filter / whole translation unit, x86-64 `-Os` (B) | 976 / 1456 | 1440 / 1456 | 1440 / 1456 | | |

* **Operation counts** are analytic loop counts of each implementation in
  one convention (`cost_model` in `tools/carma3_analysis.py` for OU-III and
  CARMA3, `tfg_cost_model` in `tools/tfg_carma3_report.py`).  TFG-CARMA3
  costs 3.4 times CARMA3 and 1.13 times OU-III.  Where it goes:
  * the prediction (6465 against 1731): CARMA3's error dynamics are block
    diagonal (attitude/biases and the three CARMA blocks never mix), while
    the invariant `Phi` couples every world column to the gyro bias
    (`Phi_{rho_i,bg} = [w_i]x F`, section 2.4) and the noise has the matching
    `[w_i]x` blocks; the implementation also applies `Phi` on both sides of
    the full `P` rather than forming one triangle in the second pass (about
    1200 of the 6465 could be saved that way; not done, the experiment does
    not need it);
  * the reset transport (5427 of the 7644 of the accelerometer update):
    `J_r(a) P J_r(a)^T` with the attitude column of `J_r` dense in every
    vector block; CARMA3 has no counterpart.  Without it (T41R0) the update
    is cheaper than CARMA3's (sparser `H`).
  The magnetometer update costs 7668 at 25 Hz (959 per 200 Hz sample).
* **Memory.**  One `NX x NX` covariance plus one covariance-sized scratch
  matrix and the `NX x 3` gain/cross-covariance scratch, all members; no
  dynamic allocation after construction; one CARMA model and one
  `StepDiscretization` (CARMA3 keeps two, horizontal and vertical), so the
  object is 240 B smaller than CARMA3's on the host.
* **Stack.**  The physical/invariant covariance mappings and the transport
  run in place through the member scratch with scalar loops (no `NX x NX`
  stack temporary, no Eigen `21 x 21` product).  The largest TFG-CARMA3 frame
  is `right_jacobian_cross` (1056 B, a leaf under the update), the next
  `time_update` (816 B); the library's `discretize_lti` (1440 B, shared with
  CARMA3) and the shared front end's Eigen triangular product (1456 B) are
  larger.
* **Float.**  The filter is templated; the target runs `float` (the ESP32-S3
  FPU is single precision).  Section 5's float-vs-double check over 60 s of
  jittered sampling: 3.5e-5 rad, 0.085 mm, covariance 1.5e-4 (correlation
  units).
* **Host time** (x86-64 `-O3`, one process on an otherwise idle host, mean
  update after Live over three JONSWAP seeds) ranks the filters and is not an
  ESP32-S3 time.  TFG-CARMA3 takes 11.3 us per sample against
  CARMA3's 7.4 us and OU-III's 22.3 us (1.5 and 0.51 times), 8.1 us
  without the reset transport.  Time grows less than the multiply-add
  count (1.5x against 3.4x) because both share the wrapper's front end,
  tuner and schedule, which the counts exclude.  This host is about three
  times slower than the one of `docs/carma-three-axis.md` (OU-III 7.07 us
  there): ratios carry over, absolute times do not.
* **ESP32-S3.**  `tests/kalman_carma/arduino/carma3_size_probe` gained the
  TFG-CARMA3 wrapper as variant 7 (`T41`), and `build.sh host` compiles the
  same eight translation units with the host compiler (the x86-64 rows above,
  `host_size_proxy.csv`).  **The Xtensa build was not run in this study**:
  executing the downloaded Espressif toolchain was not permitted in the
  session that produced it.  There is therefore no target compile, link,
  size, stack or timing evidence for TFG-CARMA3; CARMA3's target numbers
  (`reports/results/carma3_benchmark/esp32s3_size.csv`, 6040 B filter,
  135 694 B code, 1344 B largest frame) and the host proxies above (TFG
  smaller than CARMA3 in code and object size, no new frame above the
  existing ones) make a successful target build plausible but do not show
  it.  At 200-240 Hz the multiply-add count is the relevant comparison:
  14 109 per sample against the shipping OU-III's 12 483, i.e. about
  2.8-3.4 M multiply-adds per second (plus about 0.2 M for the magnetometer)
  -- if OU-III meets its deadline on the AtomS3R the TFG-CARMA3 arithmetic is
  of the same order, but that is an inference, not a measurement.  No sketch
  was added or changed.

## 12. Limitations and unresolved questions

* **No target evidence.**  No ESP32-S3 compile, size, stack or timing of
  TFG-CARMA3 exists (section 11); no hardware was used for anything.
* **Simulation only, one bench.**  The synthetic grid uses the CARMA3
  benchmark's vessel RAOs, sensor models and magnetometer disturbances;
  the pinned dataset is a second, independent simulator.  Neither contains
  a real vessel or a real IMU.  The IMU sits at the centre of gravity (no
  lever arm, which neither CARMA filter models).
* **The comparison inherits CARMA3's policies.**  Startup, staged commit,
  watchdog, magnetic-service gate, vibration inflation and every
  covariance write are CARMA3's, mapped into the invariant coordinates
  (Physical mode).  T41NW shows that mapping them natively changes nothing,
  but policies designed for the invariant filter from scratch were not
  studied.
* **Shared inconsistency.**  Attitude NEES of about 50 (3 dof) in every
  filter; section 8.4 attributes its tilt part to the accelerometer-bias
  model, which this study did not change (it is OU-III's), and leaves the
  heading part (yaw NEES 11-17) unattributed.
* **Unexplained: end-of-run gyro bias.**  The gyro-bias error over the
  last 60 s is 0.8 % higher with TFG-CARMA3 in JONSWAP, PM and large
  rotations (raw p about 0.002, not significant after correction) and
  16.7 % higher in the manoeuvring case (raw p = 0.01), while the
  whole-window gyro-bias RMS is lower in the same runs.  A per-axis
  (world-vertical against horizontal) decomposition of the bias error,
  which the benchmark does not record, is the next experiment.
* **Not group-affine.**  The `beta_g -> rho` coupling depends on the
  estimate, so the log-linear property and the estimate-independent error
  dynamics of invariant EKF theory do not apply; no stability, convergence
  or global invariance result is claimed.  Rotational equivariance holds
  for yaw with the isotropic horizontal prior the wrapper uses
  (section 2.10, checked to 1e-12) and fails for an anisotropic one.
* **Single shape and scale on all axes.**  Different horizontal and
  vertical CARMA models would add an attitude-to-translation coupling
  (section 2.10) and make the latent coordinates non-vectorial; this
  implementation refuses such a configuration rather than approximating
  it.
* **The reset trade-off is open.**  The reset transport is the correct
  first-order covariance map across a correction, yet the filter without
  it is more accurate in heading (manoeuvring, gaps) and more consistent
  (heading NEES 1.2 against 57 over 2400 s without magnetometer).  Whether
  a transport that preserves the unobservable heading direction (an
  observability-constrained reset) keeps both properties is untested.
* **The mechanism of section 8.3 is an interpretation.**  The evidence is
  the ablation ordering, the direction of the heading-covariance change,
  the dependence on the filters' own error level and the gap test; no
  experiment removed the false heading information from CARMA3's MEKF
  itself.  An observability-constrained CARMA3 MEKF is the falsifying
  experiment: if it recovers the mixed-sea and manoeuvring gains at
  CARMA3's cost, the gains belong to the unobservable-direction handling,
  not to the two-frame group.
* **Multiple comparisons.**  Holm correction is applied within each table;
  the several tables are not jointly corrected.  The post-hoc mechanism
  runs (summary section 11) were designed after the principal results and
  are labelled as such.

## 13. Recommendation

**Central question: does the coherent two-frame Lie-group representation
yield measurably better marine INS estimation than the quaternion MEKF with
the same three-axis CARMA(4,1) model?**  For the quantity the filter exists
to estimate -- wave displacement and velocity -- **no**: in 34 of 36
conditions every translational change is inside +-1 % with its whole
interval, usually inside +-0.01 %, and the two exceptions are long
magnetometer outages whose point estimates favour the MEKF without being
significant.  The derivation predicts this: the filters are first-order
equivalent, and the CARMA prior, not the error geometry, sets the
translational error.  For attitude, heading consistency and the
accelerometer bias, **yes, conditionally**: where the CARMA prior misfits
(mixed sea, manoeuvring) TFG-CARMA3 gives 2.6-6 % lower tilt error,
3-7 % lower accelerometer-bias error and 7.6-12 % lower attitude NEES
(significant after correction) and up to 13 % lower heading error (raw
p = 7e-4, not significant after the matrix-wide correction); it reacquires
heading faster after a magnetometer gap; elsewhere it is equivalent.  These gains do not touch CARMA3's known weaknesses (still
IMU, mixed-sea heave), and they come at 3.4 times CARMA3's multiply-adds.

* **Do not replace CARMA3's quaternion MEKF with TFG-CARMA3 as
  implemented**, and do not deploy it: there is no translational benefit,
  no target build evidence, and the cost is above the shipping OU-III's.
* **TFG-CARMA3 merits limited, targeted continuation, not general
  development.**  The most promising variant is the TFG *without* reset
  transport (T41R0): more consistent than the full filter everywhere
  (nearly consistent in heading through a 2400 s outage), more accurate in
  heading when manoeuvring (+23.6 % against CARMA3, 30 of 30 seeds) and
  after gaps, within 0.15 points of it elsewhere -- but without the full
  filter's tilt gain when manoeuvring (-0.2 % against +2.6 %) -- and at
  8682 multiply-adds per sample about 2.1 times CARMA3 and 0.7 times
  OU-III.  Two experiments decide whether to go further: (1) an observability-constrained CARMA3 MEKF (section 12) -- if
  it matches T41R0, the two-frame group is unnecessary; (2) a target build
  and timing of T41R0 on the ESP32-S3.
* **The verification framework stands on its own**: the independent
  references of section 5 (finite differences of the group maps, RK4
  Lyapunov `Qd`, Monte Carlo transport, the dense reference update, the
  yaw-information test) apply unchanged to any future invariant variant.

A LaTeX article was not prepared.  The primary result is an equivalence,
and the positive findings are conditional and partly post hoc; they do not
yet justify a publication-style claim.  This document and the generated
summary are the record.

## 14. Reproduction

    make -C tests/kalman_carma build                 # filters, tests, benches
    (cd tests/kalman_carma && ./tfg_carma3-test)     # section 5
    python3 tools/tfg_carma3_benchmark.py run --jobs 4        # main grid (resumable)
    python3 tools/tfg_carma3_benchmark.py followup --jobs 4   # CARMA(5,2), section 10
    python3 tools/tfg_carma3_benchmark.py mechanism --jobs 4  # post hoc, section 8
    python3 tools/tfg_carma3_benchmark.py truth               # motion levels per condition
    make ensure-sim-data && make -C tests/kalman_ou_iii kalman_ou_iii-sim
    python3 tools/tfg_carma3_benchmark.py dataset --jobs 4    # pinned dataset, 10 seeds
    python3 tools/tfg_carma3_benchmark.py traces --jobs 4
    python3 tools/tfg_carma3_benchmark.py timing              # one process, idle host
    python3 tools/tfg_carma3_benchmark.py report              # summary.md, paired_stats.json, figures
    bash tests/kalman_carma/arduino/carma3_size_probe/build.sh host   # host size proxy
    bash tests/kalman_carma/arduino/carma3_size_probe/build.sh <xtensa toolchain bin dir>   # target (not run here)

The main grid took about 3.2 h with four jobs on a 4-core x86-64 host.
Rows reproduce bit for bit (apart from `runtime_ns`) with this tree's
benches on the host that produced them (GCC `-O3 -march=native`); as for
CARMA3, a different compiler or CPU can move results in the fifth
significant digit through FMA contraction, which no conclusion depends on.
