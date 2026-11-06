# Coupled three-axis CARMA priors for the QMEKF CARMA3 filter (experimental)

**Status: experimental research, simulation only.**  Nothing here changes a
shipping estimator, a default, the wrapper, the OU-III filter or its
stability study, or anything TFG / TFG-CARMA3.  The reference is the
experimental QMEKF CARMA3 filter of [carma-three-axis.md](carma-three-axis.md)
(`Kalman3D_Wave_CARMA3<float, 4, false, 4, false>`, benchmark name `A41`),
which is preserved bit for bit (section 4).  The paper
`doc/kalman_carma/kalman_carma_coupled.tex` summarizes this study; the
generated tables are `reports/results/carma3_coupled/summary.md`.

**Result (simulation only; completed Monte Carlo, seeds 1-50, ten wave
seas plus calm water and a still IMU, complete INS).**  Shared directional
modes (`SHM:est`, 33 states, coupling estimated from the IMU alone) lower the
pooled wave-displacement RMS against the baseline A41 by **11.4 % (x), 4.8 %
(y) and 12.2 % (z)**, 95 % bootstrap intervals within +-0.6 points; the
oracle coupling adds only 0.5-1.4 points, and calm water and the still IMU
are unchanged.  Correlated excitation (`COR`) gains nothing (z +0.6 %,
horizontal -1.3 / -3.3 %): a real coupling can represent 28-41 % of the
truth's x-z cross-spectral energy, a constant phase 84-94 %.  The
constrained MCARMA (`MCA`, 35 states) is more accurate still (+14.9 / +9.5 /
+26.5 %), but its heave RAO is the benchmark hull's: with the RAO natural
frequency 27 % low it loses up to 25 % against the baseline, its worst host
update is 824 us (joint 26-state rediscretization) against 35 us, and the
unrestricted MCARMA is not identifiable (Fisher rank 51 of 54).  A wrong
coupling is harmful (reversed: -8 to -38 %).  **Recommendation: Option C**
(shared directional modes) as a separate experimental filter with a coupling
confidence and fallback, now `KalmanCarmaSharedModes`
([kalman-carma-shared-modes.md](kalman-carma-shared-modes.md)).  Nothing is
integrated into the production CARMA3 filter (section 8).

## 1. The baseline, as implemented

Verified against `src/kalman_carma/Kalman3D_Wave_CARMA3.h`,
`CarmaVerticalModel.h` and `SeaStateFusionFilter_CARMA3.h` (not against the
descriptions):

* **State (21).**  `[dtheta(3), b_g(3), x[p, v, a, l], y[p, v, a, l],
  z[h, v, a, l], b_a(3)]`, NED, z down; attitude error and gyro bias first,
  then one CARMA(4,1) block per axis, **axis-major**, then the accelerometer
  bias.  So the stacked motion block is `x_L = [x_x; x_y; x_z]` with
  `A_L = I_3 (x) A_C`, not `A_C (x) I_3` (the brief's `[p; v; a; r] (x) I_3`
  ordering is a permutation of this).
* **Per-axis model.**  `CARMA(4,1)`, `H(s) = g B(s)/P(s)`,
  `P = (s^2 + 2 zeta1 w1 s + w1^2)(s^2 + 2 zeta2 w2 s + w2^2)`,
  `B = s + w_z0`; the fitted JONSWAP shape C41J has `f1 0.306945 Hz,
  zeta1 0.104191, f2 0.815318 Hz, zeta2 0.478188, f_z0 = 0`, i.e. `B(s) = s`, a
  zero at DC.  Coordinates `[p, v, a, l]` with `p' = v`, `v' = a`, `l` the
  input of the second oscillator; the noise enters `a'` only
  (`W = g^2 e_a e_a'`), relative degree 3, so `p, v, a` are states and `a` is
  continuous with a white-noise derivative.
* **Discretization.**  `F = exp(A h)`, `Q = int_0^h exp(A t) W exp(A' t) dt`
  by the scaled Taylor series plus doublings of `discretize_lti`, and on a
  jittered clock the Taylor extension of `StepDiscretization` (exact in the
  sense of a controlled truncation below float rounding).
* **Noise injection.**  Horizontal: one unit-acceleration-variance model
  `A_h, W_1` and `Q_H = Sigma_h (x) Q_1`, `Sigma_h` the 2 x 2 horizontal
  covariance (x-y correlation `r`, zero by default; the earlier `A41C` run of
  carma-three-axis.md 9.7).  Vertical: its own model with `W = g_z^2`.  There
  is no x-z or y-z term.
* **Adaptation.**  The wrapper commits every 0.1 s (`ADAPT_EVERY_SECS`) an
  EMA of the tuner's operating point: `kappa = f_tune / f_z,base` (the same
  measured frequency for every axis) and `sigma_x = sigma_y = sigma_z =
  c_carma max(floor, sigma_aw / c_sigma)` (isotropic; no horizontal statistic
  is measured).  A commit changes `A`, `W` (and so `F`, `Q`); the posterior
  covariance is reset to the stationary prior only on entering Live.
* **Measurement.**  `f = R_wb (a - g e_z) + b_a (+ lever arm)`; the
  accelerometer rows of the CARMA states are the columns of `R_wb`; Joseph
  updates; the multiplicative attitude reset of `MekfAttitudeMath.h`.

## 2. Physics of the cross-axis relation (benchmark truth)

`tests/kalman_carma/CarmaBenchTruth.h` (physical horizontal model) is a
random-phase directional sea through explicit RAOs.  For one component of
direction `d` relative to north, `rel` relative to the bow:

    heave  z_up = H_z(w) eta,   H_z = wn^2 e^{-k 0.5} / (wn^2 - w^2 + 2 zeta wn i w)  (0.55 Hz, 0.35)
    north  x    = -i [surge cos psi - sway sin psi] eta,  surge = e^{-k(d + l_L|cos rel|)} cos rel
    east   y    = -i [surge sin psi + sway cos psi] eta,  sway  = e^{-k(d + l_B|sin rel|)} sin rel

The cross-spectral matrix `S(f) = sum_j S_eta(f) w_j t_j t_j^H` is Hermitian
PSD.  `tools/carma3_coupled_analysis.py` (section 3 of `analysis.json`):

* **x-z, y-z are in quadrature, shifted by the heave RAO.**  Phase -62 deg at
  the JONSWAP peak (not -90: the RAO lags 28 deg there), moving with
  frequency as the RAO phase does; energy-weighted coherence x-z 0.97.
* **Only 35-47 % of the x-z / y-z cross-spectral energy is real** (the
  co-spectrum) on every single-system sea (28 % x-z on the crossing sea);
  a constant complex coefficient captures 89-94 % of it (84 % x-z on the
  bimodal sea).  x-y is purely real (co-spectrum `cos th sin th`).
* Vessel response and water-particle motion differ: the heave RAO changes
  phase by 34 deg across 0.75-1.3 f_op; surge/sway are hull-averaged orbital
  motion.  No universal RAO is assumed anywhere below: Experiment C uses the
  benchmark hull's synthetic heave RAO (oracle) and mismatched ones.

## 3. Models

All motion blocks live in the filter's world (NED) frame: a wave field is
fixed in the world, so a world-frame coupling stays valid through heading
changes (tested in section 6); the wrapper's isotropic amplitudes are kept.

### 3.1 Experiment A -- correlated CARMA(4,1) (`COR`, 21 states)

Same blocks, joint driving covariance `Sigma_w` (3 x 3, real, PSD):

    F = I_3 (x) F_C,   Q = Sigma_w (x) Q_C,   P_inf = Sigma_w (x) P_C

(`Q_C`, `P_C` of the unit-acceleration-variance C41J).  Exact for identical
per-axis dynamics; checked against the joint Van Loan discretization to
1.6e-16 and against the double joint `discretize_lti` (float assembly
9e-7, scaled).  With `Sigma_w` diagonal it is the baseline: equal to it in
double to 7e-15 (positions) over 300 s of filtering (`carma_coupled-test` [1]).

**Structural limit.**  `S_p(w) = Sigma_w |H(iw)|^2` is real at every
frequency: zero (or 180 deg) cross-phase and frequency-independent coherence
`rho^2`.  Against the truth it can represent only the co-spectrum: 35-47 % of
the x-z energy (above).  The best it can do with an x-z quadrature pair is a
small zero-lag correlation `rho_xz ~ gamma cos(-62 deg)`.

### 3.2 Experiment B -- shared directional modes (`SHM`, 33 states)

**The complexified CARMA mode.**  A damped stochastic oscillator
`dq = [[-lambda, -w], [w, -lambda]] q dt + G dW` is the real form of a
complex OU process `zeta' = mu zeta + xi`, `mu = -lambda + i w`, whose
spectrum sits at positive frequency, so `Re zeta` and `Im zeta` are (nearly)
a Hilbert pair: the 90 deg relation the axes need.  One rotation has relative
degree 1 (`v` would be white); to keep `p, v, a` consistent the mode must have
relative degree 3.  The construction used here keeps the baseline's law
exactly:

> **Complexification.**  Let `P(s)` be the real CARMA(4,1) polynomial with
> poles `p_1, conj p_1, p_2, conj p_2`, `P+(s) = (s - p_1)(s - p_2)` its
> upper-half-plane factor and `D(s)` the real Hurwitz quadratic with
> `|D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2)/2` (closed form: `d_0 = sqrt(gamma)`,
> `d_1 = sqrt(beta + 2 d_0)` for `(|P+(iw)|^2 + |P+(-iw)|^2)/2 = w^4 + beta
> w^2 + gamma`).  Then the complex CARMA `zeta = B(s)/(P+(s) D(s)) xi`,
> `xi` circular complex white noise of intensity `2 g^2`, has
> `S_{Re zeta}(w) = g^2 |B(iw)|^2/|P(iw)|^2`: **its real part is exactly the
> C41J process** (Gaussian, same spectrum).  Its imaginary part is the
> quadrature companion, exact where the mirror term
> `|H(-iw)|^2/|H(iw)|^2` vanishes: 0.08 % at the peak, 0.6-2 % at
> 0.45-1 Hz, 6 % at 0.15 Hz and 17 % at 0.1 Hz (kappa = 1).

Numerically: the realification's real-part PSD equals C41J to 1.8e-14, its
`[p, v, a]` stationary covariance to 7.7e-15; poles strictly stable over
kappa 0.15-3.2.  Kinematic coordinates `[zeta, zeta', zeta'', l]` (complex),
real state `[Re(p, v, a, l), Im(p, v, a, l)]`: `p, v, a = Re(...)` are states
and the noise enters `a` only, so every accelerometer, position and velocity
row is a unit selector, as in the baseline.  (The 2021 representation of
MCARMA processes as sums of complex-valued OU processes, Fasen-Hartmann and
Scholz, is the closest published construction; it decomposes a real process,
whereas the complexification builds the single-sided companion of one.)

**Shared modes.**  Every axis carries a complexified mode; the axes share
its dynamics and are driven by circular complex noise with a Hermitian PSD
covariance `Sigma_c`.  Then

    F = I_3 (x) F_8,   Q_ab = Re(Sigma_ab) Q_8 + Im(Sigma_ab) J Q_8,   J = [[0, -I], [I, 0]]

(J commutes with `F_8`, `Q_8`; checked against the joint discretization to
1.7e-16, float 3e-6 scaled), and the displacement cross-spectrum is

    S_ij(w) = (g^2/2) [Sigma_ij |H(iw)|^2 + conj(Sigma_ij) |H(-iw)|^2]:

`Re Sigma` the co-spectrum, `Im Sigma` the quadrature spectrum, phase
`arg Sigma_ij` wherever the mirror term is small.  This is "one shared
directional wave component with axis-specific amplitude and phase" in the
form required by section 3.1 of the brief, plus the incoherent part
(spreading) carried by the full rank of `Sigma_c`.  Nesting, exact in law:
`Sigma_c` real gives Experiment A; diagonal real gives the baseline.  State
count: 8 per axis, 24 motion states.  A state-reduced form (`SHR`: the
complexified mode on z, x and y = `Re(beta zeta_z)` + C41J residuals, 16
motion states) reaches the same steady-state errors to 0.1-1 mm (analysis
section 5); it is not implemented because its coupling sits in the output
map, so adapting it would move the displacement estimate at every commit.

### 3.3 Experiment C -- constrained MCARMA (`MCA`, 35 states)

The general MCARMA(4,1) `P(s) = s^4 I + A_1 s^3 + ... + A_4`,
`B(s) = B_0 s + B_1`, `H = P^{-1} B`, `S_p = H Sigma H^*` (Marquardt and
Stelzer 2007; Schlemm and Stelzer 2012, corrected 2024) has 54 free
coefficients here.  At a generic point near the baseline its Whittle Fisher
information from 20 min of three-axis acceleration has rank 51 (the 3 null
directions are the `B -> B U` gauge), but 44 % of the coefficients have a
relative Cramer-Rao bound above 100 % (off-diagonal AR coefficients
60-115 %, `B_1` essentially unidentified): an unrestricted model is not
identifiable from one IMU in useful time.  The constrained form keeps the
physics: a complexified mode per axis, Hermitian excitation, and the vertical
output passed through a second-order heave RAO
`R_z(s) = wn^2/(s^2 + 2 zeta wn s + wn^2)`, i.e. the transfer matrix
`H(s) = diag(1, 1, R_z(s)) H_c(s)`.  The x-z phase then follows
`arg Sigma_xz - arg R_z(iw)` across frequency (test [6]: 34 deg between
0.75 and 1.3 f_op, to 0.06 deg).  12 parameters, all identifiable (RAO
`f_n` 0.8 %, `zeta` 2.9 % relative CRLB in 20 min).  The axes no longer share
dynamics, so Q is the exact discretization of the joint 26-state system
(`discretize_lti<26>`), redone when kappa or an excitation entry moves by more
than 1 % (the effect on the RMS error is below 0.1 %); the vertical mode
amplitude is the heave requirement divided by the RAO variance gain at the
current kappa.  The RAO is the benchmark hull's (labelled oracle where it
matters); `rao<fn>_<zeta>` runs mismatched ones.

### 3.4 Two shared modes, swell plus wind sea, crossing systems

Two complexified modes per axis at two frequencies with their own Hermitian
excitations (`K2`, 48 motion states) are evaluated analytically with oracle
frequencies (`analysis.json`, steady state): on the mixed sea the
horizontal/heave error falls from 89.7/93.6 mm (independent) to 27.6/29.7 mm
-- but 29.4/36.2 mm of that is the two-peak *marginal* shape (`K2_diag`, no
coupling); coupling adds 6 % / 18 %.  On two crossing systems of similar
frequency the coupling gain is ~0 (8.8 -> 8.9 mm with one mode).  K2 is not
implemented: no frequency estimator of a second spectral peak exists in the
wrapper, and its gain is dominated by the marginal shape, which is a
different question from coupling.

### 3.5 Coupling estimation (deployable)

`CarmaCrossSpectrumEstimator.h`: complex demodulation (Bingham, Godfrey and
Tukey 1967) of the three world-frame accelerations measured by the filter
itself (`R_wb^T (f - b_a) + g e_z`) at the operating frequency, a two-pole
envelope low-pass at `0.15 f_op` and a 240 s exponential average of
`A A^H` -> `S_hat`, Hermitian PSD by construction.  Mapping
(`CarmaCoupling::excitation`): the baseline diagonal (the wrapper's
isotropic `sigma_h^2 I` horizontally -- sensor noise, tilt and bias leakage in
the horizontal band are not wave energy -- and `sigma_z^2`), plus the cross
terms `b = s c S_hat_hz`, `c = |R_z|^2 sigma_z^2 / S_hat_zz`, so every
regression of a horizontal axis on z at `f_op` (co- and quadrature
spectrum) is the measured one.  `s = gain max(0, 1 - 1/(N_eff gamma^2))`
shrinks sampling noise, `gamma^2` being the debiased multiple coherence of z
on (x, y) and `N_eff` the number of independent envelope samples; `b` is
scaled down if `b^H H^-1 b` would exceed `0.999 sigma_z^2`; finally
`Sigma = R^-1 M R^-H` for a model whose axes respond with `R_i` at the band
(only C's vertical).  PSD by construction, the baseline diagonal without
coupling evidence, and equivariant under rotations about the vertical (every
factor is an invariant; test [7], 7e-16), so the prior does not depend on
how the waves sit relative to north.  With pure sensor noise the applied
coherence stays at 0.03 (test [12]).  Ablation forms: `an` adds the measured
horizontal anisotropy at the baseline trace (shrunk, capped at 90 %); `hm`
is the first form of this mapping (the measured horizontal block scaled by
`c`, floored at the prior), which read horizontal leakage as wave energy and
cost 25 % in calm water; it was replaced for that reason.  The band fraction
and horizon were chosen on screening seeds 1001-1020 only (section 5).

**Oracle versus deployable.**  `orc` runs feed the *same* estimator with the
true wave acceleration (labelled ORACLE everywhere); `est` runs use only the
IMU through the filter.  Ground-truth direction never enters a deployable
run.  The wrapper's own wave-direction estimator (`MarineWaveFrontEnd`) is not
used: the band cross-spectrum carries direction, its 180 deg sign (the sign
of the quadrature) and the coherence in one statistic.

### 3.6 Complex-OU modal form of CARMA(4,1) (`MOD`, 21 states)

**Theory.**  Fasen-Hartmann and Scholz (arXiv:2102.11681) write a d-variate
MCARMA(p, q) process as a sum of p multivariate complex-valued OU processes,
`dY_k = R_k Y_k dt + C_k dL`, with `R_k` right solvents of the AR polynomial
and `C_k` the residues of the transfer function there.  It requires an
irreducible fraction and a complete set of regular right solvents (for
instance, distinct eigenvalues of the block companion matrix); with repeated
solvents the process is not such a sum.  The components share ONE driving
process, so they are dependent, and the sampled process is VARMA(p, p-1).
The state-space characterization of Schlemm and Stelzer (2012) was corrected
in 2024 (arXiv:2411.06935).  The correction restates it through the observer
and controller canonical forms of the left and right matrix-fraction
descriptions with monic denominators, and this study uses it in that form.
All per-axis polynomials here are scalar, so left and right descriptions
coincide; the three-axis transfer matrices (`P^-1 B` with diagonal `P`,
`diag(1, 1, R_z) H_c` for C) are left fractions.

**Exact modal form of C41J.**  Scalar case (Brockwell, Davis and Yang 2011,
Prop. 2):

    p = 2 Re(z_1 + z_2),   dz_k = mu_k z_k dt + r_k g dW,   r_k = B(mu_k)/P'(mu_k)

with `mu_k` the upper poles and one common REAL Brownian motion `W`
(kappa = 1: `mu = -0.2009 + 1.9181i, -2.4497 + 4.4991i` rad/s,
`r = 0.02078 - 0.00587i, -0.02078 + 0.01289i`).  Per axis the real state is
`[Re z_1, Im z_1, Re z_2, Im z_2]` (12 motion states, as in the baseline):

* `A` = block-diag of damped rotations `[[Re mu, -Im mu], [Im mu, Re mu]]`;
  noise input `G = g [Re r_1, Im r_1, Re r_2, Im r_2]'` (residue-derived
  injection of the common noise);
* readouts `c_p = [2, 0, 2, 0]`, `c_v = c_p A`, `c_a = c_p A^2` (dense rows,
  4 terms each);
* `F = exp(A h)` closed form (rotation-decay blocks); `Q` closed form from
  `E[z_k conj z_l] = g^2 r_k conj r_l (e^{(mu_k + conj mu_l) h} - 1)/(mu_k + conj mu_l)`
  and `E[z_k z_l] = g^2 r_k r_l (e^{(mu_k + mu_l) h} - 1)/(mu_k + mu_l)`,
  including the cross-mode blocks; `h -> inf` gives `P_inf`;
* axes: `F = I_3 (x) F_4`, `Q = Sigma (x) Q_4` with a real `Sigma` (the
  baseline diagonal when coupling is off).

**Kinematics (the reason a naive modal model fails).**  The residue sums
are Markov parameters of `B/P`: `sum Re r_k = sum Re(r_k mu_k) = 0`
(`6.6e-17`, `6.0e-17`), so `c_p G = c_p A G = 0` and `p, v, a` are
noise-free readouts (relative degree 3, as in the baseline).

* Driving each mode by its own noise (same intensity) drops the cross-mode
  covariance and the cancellations: `Q_pp ~ h` (log-slope 1.00 against 4.98),
  white velocity, displacement variance x1.387, an acceleration spectrum
  rising at high frequency (+1.97 log2 per octave against -2.04): no finite
  acceleration variance.  Rejected.
* A readout rotated by a constant phase, `2 Re(e^{i phi}(z_1 + z_2))`, has
  `c G = -2 g sin(phi) sum Im r_k != 0` (`sum Im r_k = 0.00702`): relative
  degree 1.
* Per-mode complex gains `beta_k` on an axis driven by the same noise keep
  relative degree 3 iff `sum Re(beta_k r_k) = sum Re(beta_k r_k mu_k) = 0`,
  i.e. iff the axis reads `(b1 B(s) + b0)/P(s) xi` with real `b0, b1` (with
  distinct poles, the gains correspond one-to-one to real numerators of
  degree <= 3; relative degree 3 means degree <= 1).  The only admissible
  cross-axis phase is `-atan(b0/(b1 w))` (`w_z0 = 0`), never a constant, and
  `b0 = 0` leaves real coupling (Experiment A).

The complexification of 3.2 is consistent with this.  It is not a
decomposition of C41J but a complex CARMA with the upper poles and the two
real poles of `D` (-7.93, -1.25 rad/s), relative degree 3 by construction.
That is why arbitrary phases need its 8 states per axis.

**Equivalence** (`analysis.json` "modal", `carma_coupled-test` [16]): over
kappa 0.2-2.5 the transfer function matches C41J to 3e-14, the stationary
covariance to 9e-15, `F`/`Q` under the similarity `x_base = T x_modal` to
2e-15 / 7e-15, the displacement autocovariance (0.5-5 s lags) to 7e-15;
`cond(T) <= 196`.  Simulated `p, v, a` variances match `P_inf` within 2.1 %
(11 h at 50 Hz).

**Validity, adaptation, numerics.**  Distinct (regular) poles are required.
Shapes whose poles nearly coincide (`|mu_1 - mu_2| < 1e-3 max|mu|`) are
refused, and the baseline realization is the fallback.  The modal
coordinates depend on kappa (`v = c_p A(kappa) x`), so at every model commit
the filter applies `M = T_new^-1 T_old` to the state and, as `S P S'`, to the
covariance (`recoordinate_motion_`); the baseline coordinates and the
posterior are then continuous.  The modal QMEKF equals the baseline in
double to 2.7e-13 (fixed model) and 7.5e-13 (10 Hz kappa/amplitude
changes); without the re-coordinatization it does not.  Runtime arithmetic
is real; `std::complex<double>` appears only in the 10 Hz rebuild (poles,
residues, `T`, closed-form `F`, `Q`).

**QMEKF integration.**  The state ordering and the attitude, bias,
magnetometer, reset and Joseph code are the baseline's.  The accelerometer
Jacobian column block of axis `i` is `R_wb' e_i c_a` (4 terms instead of 1).
The attitude/bias-motion cross-covariances propagate through the same
congruence, and the multiplicative attitude reset does not touch motion
states.  Position, velocity and their covariance are read through the dense
rows (`pv_readout()`).

**Directional modal coupling (`MOD:est/nz`, Form `Numerator`).**  The only
kinematically admissible directional coupling of the 4-state modal form is
the per-axis numerator zero `(B + b0_i)/P` with a real excitation:
`b0_i = -w_op tan(arg S_iz) rho_iz` (`|arg| <= 70 deg`, `rho` the excitation
coherence) from the same Hermitian band estimate as B, excitation
`Re(D S D^H)` with `D = diag(1/u_x, 1/u_y, 1)`, `u_i = 1 - i b0_i / w_op`.  At
`f_op` the band gain and the phase `atan(rho tan arg S_iz)` match to 4e-16,
and the excitation stays PSD (test [17]).  The price is a nonzero
displacement density at DC on the coupled axes.  Steady state (analysis,
section 5, `MODN`): JONSWAP heave 10.2 -> 9.3 mm (SHM 9.5), horizontal
7.4 -> 9.7 mm; mixed sea horizontal 89.7 -> 68.6 mm, heave 93.6 -> 92.5 mm.
`MOD:est` (real coupling) is Experiment A in law and gives identical Monte
Carlo rows to `COR:est`.

## 4. Implementation

* `src/kalman_carma/CarmaCoupledMotion.h`: the motion blocks
  (`CarmaCorrelatedMotion`, `CarmaSharedModeMotion` =
  `CarmaKroneckerMotion<T, false/true>`, `CarmaModalMotion`,
  `CarmaRaoMcarmaMotion`).
* `src/kalman_carma/CarmaCrossSpectrumEstimator.h`: estimator, mapping,
  forms (Hermitian / Real / RealXY / Numerator), sources (Off / Estimated /
  External), frames (World / Keel), freeze.
* `Kalman3D_Wave_CARMA3<..., class Motion = void>`: one defaulted template
  parameter.  `void` is the filter as it was: every baseline line is
  unchanged inside the `!kMotion` branch, its members keep their types and
  order, and the motion block is an empty `[[no_unique_address]]` member.
  Checked bit for bit: every column (runtime excepted) of paired benchmark
  rows of A41, A52, A41S, M45, A41C, H41 on JONSWAP, of A41 and A52S on the
  mixed sea with turns and jitter, and of A41 and A52 at level 1, re-run after
  the last core change; and the ESP32-S3 object code of the A41 wrapper is
  instruction- and literal-identical to the original header (only symbol
  names gain the `void` argument).  A Motion instantiation shares the
  attitude, bias, magnetometer, reset and Joseph code line for line; only the
  linear block, its readouts (unit selectors or dense rows: `kDenseReadout`,
  `pv_readout()`), the accelerometer columns and the modal
  re-coordinatization dispatch.  The wrapper (`SeaStateFusionFilter_CARMA3`),
  the TFG headers and every default are untouched: a coupled filter is used
  as the wrapper's `MekfT`, and coupling is configured on the inner filter.
* `tests/kalman_carma/carma3/FullAdaptersCoupled.h`,
  `carma3_filters_coupled.cpp`: benchmark adapters `B41` (the baseline through
  the same adapter) and `MOD|COR|SHM|MCA:<off|est|orc>[/options]`
  (options in the header: `nz` directional modal coupling, `re`, `xy`, `an`,
  `hm`, `flip`, `rot<deg>`, `ph<deg>`, `fz<s>`, `keel`, `rao<fn>_<zeta>`, ...).
* `carma3_bench`: `--coupling` (acceleration error, whitened-innovation
  autocorrelation at 1 s and cross-correlation, applied coherence/phase,
  runtime p99/max), modifiers `sp<s>` (spreading), `magdist`, bases
  `regular_f<Hz>_h<m>`, `crossing`, `swell`, `dir_change`; the true wave
  acceleration reaches only `orc` adapters (`Sensors::a_wave_mag`).
* `tools/carma3_coupled_analysis.py` (exact analysis),
  `tools/carma3_coupled_benchmark.py` (Monte Carlo driver),
  `tools/carma3_coupled_report.py` (statistics, figures, LaTeX tables);
  data in `reports/results/carma3_coupled/`.

## 5. Tests

`tests/kalman_carma/carma_coupled-test.cpp` (102 checks, run by `make all`
through `run_tests.sh`; the same script checks `B41` = `A41` on every column
and a calm-water phantom bound for the estimated couplings), numbered as the
brief's section 9 plus [16]-[17] for the modal form:

| # | check | result |
|---|---|---|
| 1 | COR off = baseline (double / float); SHM off = baseline law; SHM real = COR | 1e-14 / 4e-6 / 4e-5 / 5e-5 relative |
| 2 | PSD excitation (2000 random Hermitian inputs, every form) and Q | min eig >= 0 |
| 3 | C2D: Kronecker / Hermitian-Kronecker / joint 26-state vs double joint; semigroup | 1e-6 / 3e-6 / 6e-6 scaled; 3e-14 |
| 4 | stable dynamics kappa 0.15-3.2; Re-part PSD = C41J | 2e-14 |
| 5 | stationary cross-covariance vs joint Lyapunov; float prediction invariance | 3e-14; 5e-4 |
| 6 | COR real; SHM phase = arg Sigma; MCA phase follows the RAO | 3e-17; 0.11 deg; 0.06 deg (of 34 deg) |
| 7 | rotation equivariance of the mapping; keel / world freeze under turns | 7e-16; exact |
| 8 | attitude/bias-motion cross-covariances through 12 deg tilt resets | COR 8e-5, SHM 2e-4 vs baseline |
| 9 | 10 Hz kappa / amplitude / coupling changes, 600 s | PSD, bounded |
| 10 | direction turned 90 deg: estimate re-converges, error recovers | yes (17.6 -> 17.2 mm) |
| 11 | excitation 1e-5 m/s^2 | finite, PSD, no displacement |
| 12 | sensor noise only: applied coherence; still-IMU phantom vs baseline | 0.03; equal |
| 13 | 3600 s SHM, 1800 s MCA, float | symmetric, PSD, bounded |
| 14 | layouts and state counts (21 / 21 / 33 / 35) | yes |
| 15 | model-matched: data from the SHM model itself | see below |
| 16 | modal form: \|H\|^2, F, Q, P_inf, autocovariance vs baseline; c G = c A G = 0; rotated readout; cond(T); independent driving; modal filter = baseline (double, fixed and 10 Hz changes); simulation; coalescing poles refused | 1e-14 / 4e-15 / 1e-13 / 3e-14 / 2e-13; 1e-16; lost; 196; Q_pp ~ h; 2.7e-13 / 7.5e-13; 2.1 %; yes |
| 17 | directional modal coupling: zeros on x, y only; relative degree 3; phase and band gain at f_op; zero phase = real coupling; PSD (300 random inputs) | yes; 3e-17; 2e-16 / 4e-16; exact; yes |

Test 15 is the decisive correctness check of Experiment B: data generated by
the shared-mode model itself (long-crested, isotropic horizontal marginals,
45 deg wave direction), with the magnetometer.  The correctly coupled filter
has RMS 16.8 / 18.1 mm (horizontal / z), against 19.1 / 22.5 mm for
independent axes, with NEES 6.77 (6 = consistent) against 7.24.  The gain
falls monotonically as the coupling gain goes 1 -> 0.5 -> 0 (18.0 / 20.9 mm
at 0.5).  A reversed coupling degrades monotonically in its gain
(21.3 / 25.3 mm at 0.5; 71.5 / 86.3 mm and NEES 380 at 1).  Every wrong
coupling is worse and more overconfident than the correct one (turned 45 deg:
46.9 mm, NEES 468; 90 deg: 85.4 mm, NEES 1785) and stays finite, and
confidently wrong coupling is flagged by the NEES.  Experiment A's form
(co-spectrum only) gains nothing on this pure-quadrature relation (identical
to independent axes).  Wrong coupling is an information source with the
wrong sign: through the accelerometer it biases heading and horizontal
position, so a deployable variant needs shrinkage and a validity monitor
(section 7).

## 6. Results (completed Monte Carlo)

Generated tables: `reports/results/carma3_coupled/summary.md` (per axis and
sea, consistency, ablations), `paired_stats.json`; figures in the same
directory.  Every row is one scenario and seed fed identically to every
filter (section 5); errors from 300 s to 1200 s; changes are
`1 - mean RMS(candidate) / mean RMS(A41)` with 95 % bootstrap intervals over
seeds (positive = better).

### 6.1 Summary table (brief, section 12)

Pooled over the ten wave seas (calm water and the still IMU excluded); RMS
in % of the incident H_s.  NEES of `[p, v]` (6 if consistent), accelerometer
NIS (3 if consistent).  Cross-spectral accuracy: energy-weighted relative
error of the estimated displacement cross-spectral matrix (10 seeds x 5 seas;
lower is better).  Runtime: host mean / 99th percentile per update (x86-64,
one job); ESP32-S3 object code / filter `sizeof` (Xtensa GCC 14.2.0, -Os).

| Model | States | X RMS | Y RMS | Z RMS | NEES / NIS | Cross-spectral error | Host us | ESP32-S3 code / filter (B) |
|---|---|---|---|---|---|---|---|---|
| Baseline CARMA3 (A41) | 21 | 4.18 | 2.47 | 3.64 | 15.42 / 0.183 | 0.226 | 3.6 / 7.9 | 141,209 / 6,136 |
| Complex-OU modal, exact | 21 | 4.18 (-0.0 %) | 2.47 (+0.0 %) | 3.64 (-0.0 %) | 15.42 / 0.183 | -- | 4.2 / 11.6 | -- |
| Modal, directional numerator, est. | 21 | 4.57 (-9.2 % [-10.7, -8.0]) | 3.03 (-22.6 % [-24.1, -21.3]) | 3.39 (+6.8 % [+6.6, +7.0]) | 12.35 / 0.220 | 0.160 | 4.3 / 12.0 | -- |
| Correlated (A), estimated | 21 | 4.24 (-1.3 % [-1.4, -1.2]) | 2.55 (-3.3 % [-3.6, -2.9]) | 3.62 (+0.6 % [+0.6, +0.6]) | 15.89 / 0.185 | 0.236 | 4.0 / 10.5 | 144,789 / 5,736 |
| Correlated (A), oracle | 21 | 4.24 (-1.4 %) | 2.57 (-3.8 %) | 3.62 (+0.6 %) | 15.96 / 0.185 | -- | -- | -- |
| **Shared modes (B), estimated** | 33 | **3.71 (+11.4 % [+11.2, +11.6])** | **2.36 (+4.8 % [+4.2, +5.4])** | **3.20 (+12.2 % [+12.0, +12.3])** | 13.91 / 0.181 | 0.182 | 7.0 / 34.9 | 153,210 / 15,288 |
| Shared modes (B), oracle | 33 | 3.66 (+12.6 %) | 2.34 (+5.3 %) | 3.14 (+13.6 %) | 13.84 / 0.182 | -- | -- | -- |
| Constrained MCARMA (C), estimated | 35 | 3.56 (+14.9 % [+14.7, +15.2]) | 2.24 (+9.5 % [+8.9, +10.0]) | 2.67 (+26.5 % [+26.3, +26.8]) | 13.57 / 0.366 | 0.143 | 26.4 / 824.0 | 160,062 / 19,352 |

Oracle rows feed the same estimator the true wave acceleration (labelled);
MCARMA's heave RAO is the benchmark hull's in every row of this table.

### 6.2 Per sea

* Shared modes, heave: long-period swell +23.7 %, two systems at different
  frequencies +16.1 %, mixed wind sea and swell +11.9 %, regular wave
  +10.1 %, Pierson-Moskowitz +8.0 %, JONSWAP +5.2 % (spread cos^8 +4.4 %),
  crossing seas +1.7 %.  The one loss is short-period wind waves, -1.2 %
  [-1.4, -1.0] (better on 2 % of seeds).  Horizontal: x +1.7 % (wind waves)
  to +27.4 % (swell); y +0.9 to +21.0 %, except the mixed sea (+0.5 %
  [-1.1, +2.0]) and the two systems at different frequencies (+0.1 %
  [-0.1, +0.3]), neither significant.
* Correlated excitation stays within 3 % of the baseline on every sea and
  axis except the mixed sea's y axis (-20.4 %): a real co-spectrum cannot
  carry the 90 deg x-z relation.
* The directional numerator form of the exact modal realization gains heave
  (+6.8 %) but loses 9-23 % horizontally (a constant real numerator puts
  displacement density at DC on the coupled axes).
* The exact complex-OU modal form reproduces A41 to the bootstrap resolution
  on every sea (the realization is exact; section 3.6).

### 6.3 Consistency

Pooled NEES is above 6 for every model because of the two-peak seas, which
a single CARMA mode cannot represent: the mixed wind sea and swell (NEES A41
62.9, SHM 60.1, MCA 56.1) and the two systems at different frequencies
(10.0, 8.4, 8.5).  On the single-system seas the coupled models are close to
consistent (NEES 5.5-6.1, 95 % position-interval coverage 0.97-0.99), the
shared modes slightly less overconfident than the baseline.  Pooled
accelerometer NIS: 0.18 for A41, COR and SHM, 0.22 for the numerator form,
0.37 for MCARMA.

### 6.4 Analysis (theoretical predictions, not Monte Carlo)

* Steady-state DARE with oracle parameters and white accelerometer noise
  only (`analysis.json`): JONSWAP heave RMS 10.2 mm independent, 9.5 mm
  shared modes, 5.5 mm MCARMA; the reversed shared-mode coupling 13.5 mm.
* Whittle Fisher information from 20 min of three-axis acceleration
  (JONSWAP): correlated 7 of 7 parameters, shared modes 10 of 10, constrained
  MCARMA 12 of 12 identifiable; unrestricted MCARMA(4,1) rank 51 of 54
  (section 3.3).

## 7. Robustness and deployability (ablations, seeds 1-30)

3-D RMS change against A41 (`summary.md`, Ablations):

* **Wrong coupling is harmful.**  Shared modes with the oracle coupling
  reversed: -8.3 to -38.5 %; turned 90 deg: -5.8 to -24.6 %; a 45 deg phase
  error still gains +1.8 to +3.3 % except on the mixed sea (-3.2 %).  A
  deployable variant therefore needs shrinkage of the estimate and a
  validity monitor with fallback to independent axes.
* **Coupling off.**  The complexified realization without coupling is
  -0.2 to -0.6 % against A41: equal in law at fixed parameters, slightly
  apart under adaptation.
* **Frame.**  A coupling frozen at 400 s keeps +1.5 to +12.9 % when held in
  NED, also through 90 deg heading zigzags (+3.2 %); frozen in the keel frame
  it loses the heading-change gain (-0.2 %).  The coupling belongs in the
  world frame (section 3).
* **Mapping.**  The first excitation mapping (measured horizontal block) cost
  22 % in calm water and was replaced (section 3.5); adding the measured
  horizontal anisotropy loses 4 % on the mixed sea.
* **MCARMA's RAO.**  With the hull RAO natural frequency 27 % low the
  constrained MCARMA loses 25.3 % (spread JONSWAP), 12.8 % (PM) and 24.4 %
  (heading changes), while still gaining 21.9 % on the mixed sea; +36 % and
  doubled damping keep gains of 8-23 %.  Its RAO prior alone, without
  coupling, already gains 2.4-7.4 %: part of the MCARMA gain is hull
  knowledge, not cross-axis coupling.

## 8. Recommendation

**Option C -- shared directional modes**, as a separate experimental filter.

* **A (keep independent CARMA3)** is not chosen: the shared-mode gains are
  consistent (every wave sea except short-period wind waves, -1.2 % heave),
  significant (intervals exclude zero) and leave calm water unchanged.
* **B (correlated CARMA3)** is rejected: no useful gain at any sea, and a
  real coupling cannot represent the quadrature x-z relation of waves.
* **C (shared directional modes)** is adopted: +11.4 / +4.8 / +12.2 %
  pooled, within 1.4 points of the oracle coupling, consistency no worse
  than the baseline, twice the host cost (7.0 against 3.6 us mean) and a
  15 KB filter object on the ESP32-S3.  It is fragile to a wrong coupling,
  so it needs the coupling confidence and fallback that the standalone
  filter adds.
* **D (full MCARMA)** is rejected: the unrestricted model is not
  identifiable; the constrained RAO form's extra gain rests on knowing the
  hull RAO (oracle here, and a 27 % frequency error turns it into a loss),
  it doubles the accelerometer NIS, and its worst host update (824 us, 16 %
  of the 5 ms step before any embedded slowdown) is unsuited to the AtomS3R.

**What goes into the production QMEKF CARMA3 filter: nothing.**  The
baseline `Kalman3D_Wave_CARMA3<float, 4, false, 4, false>` is unchanged.  The
shared-mode model is promoted to the standalone experimental filter
`KalmanCarmaSharedModes` ([kalman-carma-shared-modes.md](kalman-carma-shared-modes.md)),
which adds the coupling confidence and fallback and is not a default.
