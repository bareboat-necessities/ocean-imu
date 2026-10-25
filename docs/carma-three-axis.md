# Three-axis CARMA motion priors versus OU-III (experimental)

**Status: experimental research.**  Nothing here changes the shipping OU-III
estimator, its tuning, gates, magnetic handling or the OU-III stability study.
`Kalman3D_Wave_OU_III` / `SeaStateFusionFilter_OU_III` remain the deployed,
authoritative implementation.  The three-axis prototype lives under
`tests/kalman_carma/carma3/` so that the shipping source tree `src/` (hashed by
the OU-III stability provenance) is untouched; the PR #704 heave-only CARMA
code in `src/kalman_carma/` is used unchanged.  A CARMA prior is a statistical
model, not a law imposed on the vessel's physical history; no result below is
a stability, admissibility or hardware statement.

**Result (simulation only; no hardware).**  CARMA(4,1) on all three axes
(21 states, the same as OU-III, 25 % fewer multiply-adds per sample) is the
only configuration that beats OU-III horizontally on every simulated
condition.  In the paired complete-INS benchmark (100 seeds on the principal
seas, 30 on 27 further conditions) its horizontal RMS error is 43-73 % lower
on single-peaked, bimodal, broadband, underway, anchored, drifting and
manoeuvring cases and 17 % lower on the mixed swell + wind sea.  Its heave is
PR #704's CARMA(4,1): 27-46 % better, 13 % worse on the mixed sea.
CARMA(5,2) horizontally is better still on most single-peaked seas but 71 %
worse than OU-III on the mixed sea, because its spectral notch lands on the
swell; mixed orders inherit that failure.  A DC-zero CARMA needs no
horizontal `S` state, and OU-III's `S = 0` cannot be carried over to one (no
stabilizing Riccati solution); on CARMA(5,2), the derived `S` law changes
little.  Correlated x-y driving noise is free and gained 5 % with the
oracle correlation, untested with an estimated one.  Against this, every
CARMA variant is worse than OU-III on a still IMU (horizontal phantom 84
against 53 mm RMS, heave 53 against 27 mm) and more overconfident on the
mixed sea.  The models are re-discretized on every sample whose interval
changes.  **Recommendation: B, CARMA(4,1) on X/Y and Z, as the research
candidate; not shippable** (section 12).

## 1. Starting point

* **Commit.**  `main` at `1f9fa8c` (merge of PR #704, head `8b519ba`).  The
  PR #704 branch `claude/carma52-heave-filter` has one later commit,
  `ad740fb`, a CI bot regeneration of the OU validation evidence (no source
  change): the CARMA code is the one that produced the PR #704 benchmarks.
* **Reproduction.**  `make -C tests/kalman_carma build test` passes (801
  model checks, 30 filter checks, process validation, smoke replays).  Rows of
  `reports/results/carma_benchmark/raw_runs.csv` re-run on this host are bit
  identical (e.g. complete INS, JONSWAP 0.29 Hz and mixed sea, seed 4, OU-III,
  C41, C52: heave RMS and NEES identical to every printed digit).
* **Existing design.**  OU-III carries per axis the integrated OU chain
  `[v, p, S, a_w]`, `a_w` an OU process with time constant `tau` and
  stationary std `sigma_aw` (isotropic: `sigma_aw,H = S_factor sigma_aw,Z`,
  `S_factor = 1`), `S' = p`, and the pseudo-measurement `S = 0` every `T_S`
  with the SpectralMSE `r_S` (horizontal factors 0.72, 0.50).  PR #704 replaced
  only the vertical chain by CARMA(p, p-3) in `[h, v, a, latent]` coordinates
  (22 states for CARMA(5,2), no `S_z`), keeping OU-III's horizontal chain,
  attitude, bias, magnetometer, startup and adaptation.
* **Defect found: convergence time.**  `convergence_time()` of
  `carma_bench` summed squared errors over a rolling window that starts at the
  first sample, where the error is NaN before Live; once a NaN entered the sum
  every later window was NaN, and `std::max(0.0, NaN)` returned `0.0`.  Every
  tuner-driven run (all 6610 complete-INS rows and all level-1 tuner rows of
  `raw_runs.csv`) therefore reported `conv_s = 0`, and the PR #704 sentence
  "convergence time does not separate the filters" was an artifact.  Fixed in
  `carma_bench.cpp` (the window restarts after a NaN); with the fix, seed 4 of
  the complete INS gives OU-III 101 s, C41 95 s, C52 90 s on JONSWAP and
  OU-III 106 s, C41 1108 s, C52 never on the mixed sea.  The committed
  `conv_s` column of the PR #704 results is invalid; the three-axis benchmark
  reports the corrected metric.  No other committed PR #704 number changes
  (the fix and the refactors below are verified bit identical on every other
  field).
* **Frequency-range consistency.**  The CARMA amplitude mapping sets the
  model's *total* (0 to infinity) acceleration variance equal to an RMS
  measured in the adaptive band [0.5 f, 4 f] of the shipping sigma channel;
  the model's relative-degree-3 tail (`S_a ~ w^-2`) puts part of its variance
  outside that band, so the in-band prior is lower than measured.  The
  process-validation PSD checks of `tools/carma_process_validation.py`
  compare the acceleration PSD only on the *displacement* band (within 40 dB
  of the displacement peak, below 2 Hz), not on the band where acceleration
  carries its energy.  Both are reported, not changed: the mapping is the one
  every benchmarked CARMA filter uses (changing it would be re-tuning), and
  the steady-state analysis below shows that a *larger* horizontal amplitude
  than the true one is favourable under shape mismatch (section 6), so the
  inconsistency is not what limits CARMA.

## 2. Physics of three-axis wave-induced motion

Earth-fixed NED axes; `p_wave` is the oscillatory displacement about the local
equilibrium (the repository's wave coordinate).

**Heave.**  Hydrostatic restoring gives a resonant RAO (`f_n` = 0.55 Hz,
`zeta` 0.35 for the benchmark hull) times hull averaging `e^{-k d}`; the heave
spectrum is the incident JONSWAP/PM spectrum reshaped by it, single-peaked in
single seas, two-peaked in swell plus wind sea.

**Surge and sway.**  There is no hydrostatic restoring in the horizontal
plane.  A small free-floating hull in long waves follows the water-particle
orbit: in deep water the horizontal orbital displacement along the
propagation direction has the elevation amplitude and lags it by a quarter
period (`xi = H[eta]`, `xi ~ w int eta` for a narrow band).  Hull averaging
over draft, length and beam attenuates short waves; the benchmark uses
effective averaging lengths (the `e^{-k l}` form the truth already uses for
heave, roll and pitch, `l` = 0.09-0.14 L matches `sinc(kL/2)` at
`kL/2` = 1..pi/2):

    surge = e^{-k (d + l_L |cos b|)} cos b,   sway = e^{-k (d + l_B |sin b|)} sin b,
    d = 0.5 m, l_L = 0.75 m, l_B = 0.28 m,  b = propagation direction re bow.

So the horizontal spectra are the incident spectrum times `cos^2 b` /
`sin^2 b`, low-passed by the hull, with **no resonance** (the PR #704 truth
used 0.8 x the heave RAO, whose 0.55 Hz resonance is not physical for surge
and sway; it is kept as the `hlegacy` control).  Horizontal displacement is
about half the heave in these seas (JONSWAP 0.29 Hz, 0.8 m: std x/y/z
0.103/0.080/0.202 m) and its zero-crossing frequency is lower than heave's
(0.31 against 0.33 Hz), because heave's resonance lifts its shoulder.

**Directional coupling.**  For one direction component the three axes are
perfectly coherent: `x, y` in phase with each other (co-spectrum
`cos th sin th`), and in quadrature with the elevation.  The x-z phase at the
heave peak is -63 deg, not -90 deg: quadrature holds against the *incident*
elevation and the heave RAO adds its own lag (27 deg at 0.29 Hz).  Directional
spreading lowers coherence: JONSWAP with cos^8 spreading has energy-weighted
squared coherence x-z 0.97, y-z 0.76, x-y 0.64 (figure `coherence.svg`);
two systems 70 deg apart (bimodal) make x-y and y-z incoherent where the
other system dominates.  The cross-spectral matrix
`S_xyz(f) = sum_j S_eta(f) w_j t_j t_j^H` is Hermitian positive semidefinite
by construction; numerically the Hermitian defect is <= 4e-16 and the
smallest eigenvalue is >= 0 (relative to the trace) at every frequency of
every case.

**Roll and pitch coupling.**  Roll and pitch follow the beam and bow slopes;
the benchmark's IMU is at the centre of gravity, so roll/pitch reach the
translational estimate only through the attitude estimate: a tilt error
`dtheta` puts `g dtheta` into the horizontal specific force (0.1 deg is
0.017 m/s^2), which is why the horizontal channel is the attitude-coupled
one.  A lever arm (IMU above the roll centre) would add `r x alpha` to the
horizontal acceleration; it is studied elsewhere (`docs/ou-iii-imu-lever-arm.md`)
and not here.

**Low-frequency translation is not wave motion.**  Second-order
(difference-frequency) forces and wind drive horizontal motion far below the
wave band: an anchored hull sails on its rode with 30-300 s periods and metres
of amplitude (benchmark: 2.5 m at 1/90 Hz, 1.1 mg, with 16 deg yaw sheer); a
drifting hull moves with the mean drift and slowly varying velocity
(0.3 m/s mean, 0.08 m/s fluctuations, 0.7 mg); an underway hull makes way at
constant velocity, which no accelerometer sees, and changes speed and heading
(benchmark manoeuvre: up to 0.2 m/s^2 for 25-35 s).  These are motions of the
local equilibrium, not of the wave coordinate: the benchmark adds them to what
the IMU senses and scores the estimators against `p_wave` only.  Their
accelerations (1-30 mg) are of the order of accelerometer bias and of
`g dtheta`, so an IMU-only estimator cannot reconstruct them; what matters is
how much of them each prior turns into phantom wave displacement.

**Forward speed.**  Underway, every component is met at the encounter
frequency `w_e = w - k U cos b`: head seas move energy up (0.29 Hz to about
0.43 Hz at 2.5 m/s, with heave acceleration doubling), following seas move it
down towards and through zero (frequency folding), where the double
integration from acceleration has no information at all.

**Consequences for the prior.**  (i) The wave-induced horizontal spectrum of
a free-floating hull is, like heave, the incident spectrum reshaped by a
smooth RAO, with no extra resonance: forcing a JONSWAP-*fitted* shape onto
surge and sway is justified for the oscillatory part to the same extent as for
heave (a fixed shape scaled to the measured wave frequency), not more.
(ii) Both the elevation spectrum (`exp(-1.25 (f_p/f)^4)`) and every RAO vanish
at DC, and the displacement primitive of wave motion is bounded (a premise of
the repository's physical model), so a physically consistent prior has a
**zero at DC**.  (iii) Low-frequency translation is outside the model class;
a stationary prior should attribute it to bias/tilt or ignore it, not track
it.  (iv) Cross-axis coherence is high but its phase and magnitude depend on
direction, spreading and the hull's RAOs.

## 3. Candidate state-space models

A CARMA(p, q) displacement `p` with `H(s) = g B(s)/P(s)` (monic), driven by
unit-intensity white noise `xi`, is realized in the library coordinates of
`src/kalman_carma/CarmaVerticalModel.h`:

    x = [p, v, a, l_1 .. l_{N-3}],   dx = A x dt + G dW,   G = g e_a,
    p = e_0' x,  v = e_1' x = C A x,  a = e_2' x = C A^2 x,   W = G G' = g^2 e_a e_a'.

Relative degree `p - q = 3` gives `C G = C A G = 0`, `C A^2 G = g`:
displacement and velocity are differentiable, acceleration is continuous with
finite variance (`S_a = w^4 S_p ~ w^-2`), its derivative carries the white
noise.  Exact discretization `F = e^{A h}`, `Q = int_0^h e^{At} W e^{A't} dt`
uses the library's scaled Taylor series plus doublings (no stationary
difference `P - F P F'`, which cancels in float).  One shape at `kappa = 1`
is moved by the measured frequency (`kappa = f_tune / f_z,model`) and the
amplitude `g^2 = sigma_a^2 kappa / Var_1(a)` (closed form, no run-time
Lyapunov solve).

| Model | B(s) | P(s) | states | DC |
|---|---|---|---|---|
| CARMA(3,0) | 1 | (s + lambda) L1 | 3 | level |
| CARMA(4,1) | s + w_z0 | L1 L2 | 4 | zero (fitted `w_z0` = 0) |
| CARMA(5,2) | s^2 + 2 zeta_z w_z s + w_z^2 | (s + lambda) L1 L2 | 5 | level (-23 dB re peak) |
| CARMA(5,2), DC | s^2 | (s + lambda) L1 L2 | 5 | double zero |

`L_j = s^2 + 2 zeta_j w_j s + w_j^2`.  The library writes the latent
coordinates of the cascade realization (CARMA(4,1): the input `r` of `L2`;
CARMA(5,2): `(r, r')`); `w_z = 0` in the CARMA(5,2) branch gives `B = s^2`
exactly, so the double-DC-zero model needs no new code.

Values at the JONSWAP operating point (`f_z` = 0.37 Hz, `sigma_a`
= 0.5 m/s^2; `tools/carma3_analysis.py`, `carma3_analysis.json` holds A, G,
the output rows, W, F, Q at 200 Hz and P_inf of every model):

| | C30J | C41J | C52J | C52D | C41H |
|---|---|---|---|---|---|
| poles (rad/s) | -30.7, -0.63 +- 2.29i | -2.45 +- 4.49i, -0.20 +- 1.91i | -11.7, -0.17 +- 1.74i, -1.45 +- 2.37i | -30.9, -0.14 +- 1.76i, -1.33 +- 2.85i | -0.85 +- 3.23i, -0.15 +- 2.09i |
| controllable / observable from a / from p | 3/3/3 | 4/4/4 | 5/5/5 | 5/5/5 | 4/4/4 |
| `C A^2 G` (= g) | 3.75 | 1.27 | 2.27 | 3.65 | 0.48 |
| transfer-function error vs `g B/P` | 2e-12 | 5e-14 | 5e-13 | 4e-12 | 4e-14 |
| std p / v (m, m/s) | 0.032 / 0.075 | 0.058 / 0.134 | 0.047 / 0.109 | 0.036 / 0.084 | 0.081 / 0.189 |
| cond(P_inf), correlation-scaled / raw | 2.1 / 274 | 5.5 / 144 | 12 / 534 | 8.3 / 1426 | 17 / 181 |
| `Q_pp / Q_aa` at 200 Hz | 3.3e-11 | 3.2e-11 | 3.2e-11 | 3.4e-11 | 3.1e-11 |

C41J/C52J/C30J are the PR #704 JONSWAP fits; C52D is CARMA(5,2) with `B = s^2`
refitted to JONSWAP (1.12 dB, real pole on its 5 Hz bound); C41H is CARMA(4,1)
fitted to the benchmark's beam-sea sway spectrum (JONSWAP times its own
horizontal RAO, energetic band; **hull-dependent**).  Every model is strictly
stable, minimal, observable from acceleration alone and from displacement.
The correlation-scaled condition number is independent of `kappa` (scaling
is a similarity); the raw one grows with frequency spread (`Var p ~ kappa^-5`
against `Var a`), which the library's scaled discretization handles: float
`F`, `Q` agree with double to <= 1.5e-6 (relative to each entry's scale) over
f_z = 0.05-1 Hz for C41, C52 and both with S (`carma3-test`).

**Fits to horizontal targets.**  Log-spectral shape fits (objective and bounds
of `tools/carma_fit.py`) to the vessel-filtered targets need the energetic
band (within 30 dB of the peak): the hull filter makes the targets fall
30-75 dB across 0.17-1.30 Hz, faster than any relative-degree-3 rational
spectrum, and a full-band fit would be decided by the deep tail.  The
CARMA(4,1) fits to vessel heave, head-sea surge and beam-sea sway are close to
each other (f1 0.28-0.29 Hz, zeta1 0.073-0.075, f2 0.41-0.48 Hz,
zeta2 0.18-0.26) and all differ from the JONSWAP fit (f2 0.815 Hz,
zeta2 0.48): the hull filter, not the axis, changes the shape.

## 4. The auxiliary state S

OU-III's horizontal chain needs `S = int p dt` and `S = 0`: its `p` is a
double integral of an OU acceleration, so without `S` the pair `(v, p)` is not
detectable from acceleration.  A CARMA prior is Hurwitz, so it is detectable
(indeed observable from acceleration) without any pseudo-measurement.  The
question is whether the integral of CARMA displacement is controlled.

**Theory.**  `S(t) = int_0^t p` has transfer function `g B(s) / (s P(s))`.

* If `B(0) != 0` (CARMA(5,2) as fitted, CARMA(3,0)), the pole at 0 does not
  cancel: in the augmented system `[x; S]` the zero mode is **controllable
  from the noise and unobservable from acceleration** (PBH, both checked
  numerically).  `S` is a random walk with `Var S(t) ~ S_p(0) t`, `S_p(0)`
  the two-sided DC level; at f_z 0.37 Hz, 0.5 m/s^2: slope 5.497e-5
  m^2 s^2/s analytic, 5.497e-5 from the covariance recursion, 1.232e-4 both in
  the float filter at its test operating point (`carma3-test`).  The prior
  therefore *contradicts* the physical premise that the displacement
  primitive is bounded.
* If `B(0) = 0` (CARMA(4,1) as fitted, CARMA(5,2) with `B = s^2`), the zero
  mode is **uncontrollable**: `S(t) = S_inf(t) - S_inf(0)` with the stationary
  `S_inf = g B~(s)/P(s)`, `B = s B~`, which is a linear functional of the
  existing state (`c' x`, `c = A^{-T} e_p`).  `Var S(t) -> 2 Var S_inf`
  (C41J: 0.00182 = 2 x 0.000908 m^2 s^2; C52D: 0.000688 = 2 x 0.000344), with
  no growth over 3600 s.  The prior already encodes a bounded primitive.

**Option A -- CARMA without S.**  With a DC zero the stationary prior bounds
the displacement *and* its primitive; without one (C52J) it bounds the
displacement but lets its primitive wander as a random walk.  Drift control
of the displacement itself comes from the prior in both cases: the error at
DC is set by the prior's DC level and the bias/tilt confounding, not by S.

**Option B -- CARMA with S.**  `S' = p` appended adds an eigenvalue at 0.
* `B(0) != 0`: the mode is a controllable random walk; `S = 0` every `T_S`
  with variance `r^2` is an information-rate-`1/(r^2 T_S)` regularizer that
  closes a first-order loop at `w_R = sqrt(S_p(0)/(r^2 T_S))`; placing
  `w_R tau = kappa_R` (OU-III's pole constant 0.3627) gives
  `r_S = sqrt(S_p(0)) tau / (kappa_R sqrt(T_S))` per axis (the PR #704 `S_z`
  law, used for A52S).  It is a regularizer, not a measurement: under the
  prior `S` is unbounded, so `S = 0` asserts the physical premise the prior
  lacks.  The augmented pair is detectable (`[C_a; e_S]` observes the zero
  mode) and the steady-state filter exists.
* `B(0) = 0`: the mode is uncontrollable *and* on the unit circle, so the pair
  is **not stabilizable**: the discrete algebraic Riccati equation has no
  stabilizing solution (scipy refuses it).  The filter estimates the offset
  `S(0)` like a constant (gain -> 0), and the pseudo-measurement then keeps
  re-measuring `S_inf = c' x` -- prior information already in the model,
  counted again every `T_S` with temporally correlated "noise".  It is not a
  valid measurement model.  Carried over with OU-III's r_S (A41S, the negative
  control) it is still worse in steady state on PM, mixed and bimodal seas
  (table in section 6).
* Numerically, the zero mode integrates every rounding error of `Q`: a Van
  Loan step of 1 s with the 31 rad/s real pole of C52D leaks 3e-7 m^2 s^2 per
  step into the uncontrollable mode (the analysis therefore uses 0.1 s steps),
  and in float the `S` variance seeded at 50^2 m^2 s^2 has an ulp (2.4e-4)
  far above the per-sample random-walk increment (~1e-7), so the S-mode
  covariance of a float filter is only as good as its last pseudo-update.

**Option C -- the DC zero as the low-frequency constraint.**  A numerator
zero at DC makes the integral of displacement stationary with finite variance
(single zero: `S` stationary, `int S` a random walk; double zero: both
stationary).  It is the principled replacement of the S mechanism, and it is
what OU-III's S mechanism *is* in the frequency domain:

> **OU-III's S = 0 is a DC zero.**  Condition the integrated-OU prior
> (`S_p,OU = 2 sigma^2 tau / (w^4 (1 + w^2 tau^2))`) on `S = 0` observed at
> information rate `1/N`, `N = r_S^2 T_S`.  The conditional displacement
> spectrum is `S_eff = S_OU N w^2 / (S_OU + N w^2)
> = 2 sigma^2 tau N w^2 / (2 sigma^2 tau + N w^6 + N tau^2 w^8)`, a rational
> spectrum `(2 sigma^2/tau) |iw|^2 / |P(iw)|^2` with `P` the monic Hurwitz
> factor of `N tau^2 s^8 - N s^6 + 2 sigma^2 tau`: a **CARMA(4,1) with
> B(s) = s**.  (Spectral factorization error 4e-15..1e-14.)

At the JONSWAP 0.29 Hz operating point its poles are 0.293 Hz / zeta 0.39 and
0.299 Hz / 0.93 for the vertical r_S, 0.318 / 0.39 and 0.324 / 0.93 for the
horizontal (x factor 0.72): a broad, low-Q DC-zero prior, against the
C41J fit's sharp 0.307 Hz / 0.10 resonance.  The equivalence is exact for the
two-sided (smoothing) posterior.  The causal filters differ (the causal OU-III
filter conditions only on past S updates): their steady-state transfers differ
by up to 50 % in 0.1-1 Hz, and the causal OU-III error is 13.4 mm against
16.3 mm for its equivalent CARMA(4,1) on heave, 9.9 against 8.6 mm on x
(`ou3_s_equivalence.svg`).  So OU-III already *is* a broadband DC-zero
prior; CARMA(4,1) as fitted is a *narrow* DC-zero prior, and the CARMA(5,2) fit
gives up the DC zero for a notch at 0.15 Hz (x kappa).

**Recommendation on S.**  Do not carry OU-III's `S = 0` into a CARMA model.
Use a DC-zero CARMA (Option C; CARMA(4,1) has it as fitted).  If a model
without a DC zero is kept, its `S` needs the derived pole-placement law
(Option B, A52S), and it remains a regularizer of a mis-specified DC level.
The benchmark tests A (A41, A52), B (A52S; A41S as the negative control) and
C (A41's DC zero).

## 5. Multivariate (coupled) CARMA

**Independent axes.**  `A = diag(A_x, A_y, A_z)`, `W = diag(...)`.  The
prototype shares one `A_h` for x and y (the same waves, one measured
frequency) and gives z its own model and order.

**Correlated process noise.**  With `A_x = A_y = A_h`, a joint driving
intensity `Sigma_h (x) W_1` (`Sigma_h` a 2 x 2 directional covariance) is a
legitimate coupled vector CARMA: `F = I_2 (x) F_1`, `Q = Sigma_h (x) Q_1`
exactly (checked against the joint discretization to 6e-16 and a fine-step
reference to 1e-14 in `carma3-test`), at no extra state or cost.  Its
cross-spectrum is `S_xy(w) = Sigma_xy |H(iw)|^2`, **real at every frequency**:
correlated noise through identical dynamics can only represent the in-phase
x-y co-spectrum of directional seas.  The x-z and y-z relations of linear
waves are quadrature (phase-shifted), which correlated noise cannot produce;
they need different dynamics on the two axes.

**Shared-latent quadrature coupling.**  The minimal dynamic coupling uses the
narrow-band identity `xi ~ w int eta`: with a DC-zero heave prior,
`int z = S_inf = c' x_z` is a state functional, so

    p_N = alpha_N w_z c' x_z + r_N,   a_N = alpha_N w_z v_z + a(r_N),

with independent CARMA residuals `r_N, r_E`, is a linear Gaussian vector
CARMA (relative degree 4 for the coherent part).  Its cross-spectral matrix
is Hermitian PSD (smallest eigenvalue / trace >= 1e-5 on 0.02-2 Hz) and the
benchmark-truth cross-spectrum fixes `alpha` (oracle direction and coherence).

**How much can cross-axis information give?**  Noncausal Wiener bound with
the *true* `S_xyz` and white accelerometer noise (`tools/carma3_analysis.py`):

| sea | x gain, h-noise x1 / x3 / x10 | x-y only, x1 / x3 / x10 |
|---|---|---|
| JONSWAP 0.29 Hz (heading 0) | 6.3 / 13.5 / 20.8 % | 0.7 / 1.2 / 1.6 % |
| PM 0.25 Hz | 5.3 / 12.2 / 18.7 % | 0.5 / 1.1 / 1.4 % |
| mixed swell + wind | 16.5 / 39.8 / 64.9 % | 0.5 / 0.6 / 0.9 % |
| bimodal | 10.7 / 24.0 / 38.7 % | 0.0 / 0.0 / 0.0 % |

("h-noise x3" inflates the horizontal measurement noise threefold, a crude
stand-in for tilt-coupled error.)  The bound is relative to an *oracle prior*
whose own error is tiny (0.3 mm): it says what cross-axis information adds
once the spectrum is known exactly, not what a mismatched prior can recover.
The causal steady-state errors with real priors are 10-30 times larger and
**mismatch-dominated**.  The shared-latent coupled model, even with oracle
`alpha`, is worse than independent C41 axes on single-peaked seas
(x: 10.2 against 5.5 mm, JONSWAP; 14.1 against 8.8 mm, PM; inflating its
residual to the isotropic amplitude: 8.8 and 12.9 mm), because the narrow-band
quadrature `w_z / w` is off by 25 % at 0.8 `w_z` and the error budget is
dominated by exactly such mismatch.  It is better only on the mixed sea
(67 / 61 mm against 88 mm), where the independent C41 prior is itself badly
mismatched.  The x-y correlation alone is worth at most 0.5-1.6 % in the
bound; in the complete INS, with a mismatched prior, the oracle correlation
gained 5 % (section 9.7).

**Conclusion on coupling.**  A physically meaningful cross-spectral model
needs direction (with its 180 deg ambiguity), spreading and the hull's RAO
phase, all sea-state dependent; even with all three given exactly, the
achievable causal coupling structure loses to independent axes where CARMA
works, and the gain bound is small there.  A coupled model is not practical
for the real-time filter at this stage; correlated noise across identical x/y
dynamics is practical (zero cost) but can represent only the co-spectrum, has
a small bound, and needs an estimate of the correlation.

## 6. Steady-state error budgets (white noise only, oracle operating point)

Exact steady-state Kalman transfers (DARE at 200 Hz) applied to the true
displacement spectra; error = signal distortion / noise propagation, mm RMS.
Horizontal priors use the wrapper's mapping (frequency from the vertical
channel, amplitude = vertical band acceleration RMS, i.e. OU-III's isotropic
choice).

| x axis, mm RMS (signal / noise) | JONSWAP 0.29 Hz | PM 0.25 Hz | mixed | bimodal |
|---|---|---|---|---|
| OU-III | 10.0 (7.8 / 6.2) | 13.7 (11.7 / 7.1) | 70.2 (69.1 / 12.2) | 33.2 (29.9 / 14.4) |
| C30J | 18.4 (18.3 / 1.8) | 28.0 (27.9 / 2.0) | 113.7 (113.7 / 4.2) | 58.7 (58.4 / 5.6) |
| C41J | 5.6 (4.5 / 3.2) | 8.8 (8.0 / 3.7) | 87.6 (87.3 / 6.8) | 23.3 (21.8 / 8.4) |
| C52J | 6.8 (6.6 / 1.8) | 12.7 (12.5 / 2.0) | 283.7 (283.7 / 3.8) | 46.1 (45.8 / 4.8) |
| C52J+S | 6.5 (6.3 / 1.8) | 13.6 (13.5 / 2.0) | 283.1 (283.1 / 3.8) | 49.3 (49.0 / 4.8) |
| C41J+S (OU-III r_S) | 5.6 (4.7 / 3.1) | 10.6 (10.0 / 3.5) | 105.1 (104.9 / 6.4) | 30.1 (29.0 / 7.9) |
| C52D | 10.7 (10.6 / 1.3) | 19.0 (18.9 / 1.5) | 289.4 (289.4 / 2.8) | 56.2 (56.1 / 3.7) |
| C41H | 4.8 (2.9 / 3.8) | 7.5 (6.1 / 4.4) | 80.4 (80.0 / 7.6) | 23.2 (21.5 / 8.8) |
| C41J, true amplitude | 6.4 (5.9 / 2.5) | 10.4 (10.0 / 3.0) | 132.7 (132.6 / 4.6) | 33.3 (32.8 / 5.8) |
| *heave (z), same filters: OU-III / C41J / C52J* | 13.4 / 10.2 / 13.7 | 17.2 / 14.8 / 21.3 | 60.7 / 93.7 / 302.1 | 31.5 / 26.9 / 51.3 |

* Horizontal C41J halves OU-III's x error on single-peaked seas (signal
  distortion and noise propagation both lower), and is worse on the mixed
  sea.  C52J is worse than C41J on every sea, catastrophically on the mixed
  sea; C52D (double DC zero) and C30J are worse than OU-III.  The horizontal
  fit C41H is better than C41J on single-peaked seas.
* The *true* per-axis amplitude (`C41J_oracle_amp`) is worse than the larger
  isotropic one: under shape mismatch an inflated prior distorts less.
* This analysis has no bias, no tilt coupling and no adaptation; the Monte
  Carlo below has all three, and its ranking of C41 and C52 differs.

## 7. State dimension and cost

Every configuration shares nine states: attitude error (3), gyro bias (3) and
accelerometer bias (3).  Motion states per axis: OU-III `[v, p, S, a_w]` (4),
CARMA(4,1) `[p, v, a, r]` (4), CARMA(5,2) `[p, v, a, r, r']` (5), plus one
for an `S`.  `P` is a dense symmetric `NX x NX` float matrix.

| Configuration | id | states x / y / z | NX | P floats / bytes | multiply-adds per sample | ESP32-S3 `sizeof` filter (B) | ESP32-S3 code, object (B) | largest stack frame (B) | host update, fixed / jittered dt (us) |
|---|---|---|---|---|---|---|---|---|---|
| OU-III (shipping) | OU3 | 4 / 4 / 4 | 21 | 441 / 1764 | 12 483 | 5692 | 156 309 | 1024 | TIMING_OU3 |
| OU-III x/y + C41 z (PR #704) | Z41 | 4 / 4 / 4 | 21 | 441 / 1764 | 9315 | 4760 | 150 766 | 1344 | TIMING_Z41 |
| OU-III x/y + C52 z (PR #704) | Z52 | 4 / 4 / 5 | 22 | 484 / 1936 | 10 344 | 5224 | 153 358 | 2000 | TIMING_Z52 |
| C41 x/y/z | A41 | 4 / 4 / 4 | 21 | 441 / 1764 | 9315 | 5144 | 136 799 | 1344 | TIMING_A41 |
| C52 x/y/z | A52 | 5 / 5 / 5 | 24 | 576 / 2304 | 12 582 | 6448 | 141 055 | 2000 | TIMING_A52 |
| C41 x/y + C52 z | M45 | 4 / 4 / 5 | 22 | 484 / 1936 | 10 344 | 5608 | 147 121 | 2000 | TIMING_M45 |
| C52 x/y + C41 z | M54 | 5 / 5 / 4 | 23 | 529 / 2116 | 11 433 | 5952 | 148 058 | 2000 | TIMING_M54 |
| C52 x/y/z + S_x, S_y | A52S | 6 / 6 / 5 | 26 | 676 / 2704 | 15 232 | 7296 | 145 999 | 2032 | TIMING_A52S |
| C41 x/y/z + S_x, S_y | A41S | 5 / 5 / 4 | 23 | 529 / 2116 | 11 433 | 5880 | 143 966 | 1536 | -- |
| C41 x/y/z, correlated x-y noise | A41C | 4 / 4 / 4 | 21 | 441 / 1764 | 9315 | as A41 | as A41 | as A41 | -- |
| shared-latent coupled C41 (section 5) | -- | 4 / 4 / 4 | 21 | 441 / 1764 | about 9315 | not implemented in C++ | | | |

* **Operation counts** are the analytic multiply-adds of one prediction plus
  one accelerometer update of the structured implementation (block-diagonal
  `F_LL` applied per axis, symmetric Joseph update; `cost_model` in
  `tools/carma3_analysis.py`).  CARMA(4,1) everywhere is 25 % cheaper than
  OU-III, whose implementation propagates its 12-state linear block as one
  dense matrix; CARMA(5,2) everywhere costs the same as OU-III; the S
  states make CARMA the most expensive.  Correlated x-y noise and the
  shared-latent coupling add no state: `Q_H = Sigma_h (x) Q_1` and the
  coupling enters `C` and `A`, not `NX`.
* **ESP32-S3 numbers** are from the target compiler (Espressif Xtensa GCC
  14.2.0, `-Os`, arduino-esp32 flags) on one translation unit per variant
  holding the complete fusion wrapper's sample path
  (`tests/kalman_carma/arduino/carma3_size_probe`, `esp32s3_size.csv`); code
  bytes are object-level, not a linked image.  All variants are smaller in
  code than OU-III (one CARMA model class serves every axis); the filter
  object grows roughly with `NX^2` (+736 to +2152 B over C41 x/y/z for two to
  five added states); the largest stack frame is the
  library's `discretize_lti` (2000 B for a 5 x 5 model), against 1024 B for
  OU-III.  **No hardware was used: no on-target time, float behaviour or
  AtomS3R validation is claimed.**
* **Host time** (x86-64, `-O3`, one process on an idle host, mean update
  after Live over three JONSWAP seeds) ranks the variants; it is not an
  ESP32-S3 time.  An absolute target estimate would need a measured
  host/target ratio for this code, which does not exist; the defensible
  statement is relative: every C41 variant costs less per sample than the
  OU-III filter that the AtomS3R runs today.
* **Discretization.**  `F`, `Q` are recomputed only when the sample interval
  or the committed operating point changes (one horizontal model serves x and
  y through the Kronecker structure).  The library compares the interval
  exactly, so a jittered sample clock re-runs the scaled Taylor/doubling
  exponential of every CARMA block on every sample (the "jittered" column);
  OU-III's closed-form transition has no such cost.  A dt quantized to the
  nominal rate, or a cached table, would remove it at a controlled error;
  this is not implemented (PR #704's vertical model behaves the same way).
* **Conditioning** is section 3: correlation-scaled `cond(P_inf)` 5.5 (C41)
  to 12 (C52), raw up to 1426; float `F`, `Q` within 1.5e-6 of double over
  f_z 0.05-1 Hz; 100 s of float prediction keeps the stationary covariance
  to 1.5e-4.

## 8. Implementation

* `tests/kalman_carma/carma3/Kalman3D_Wave_CARMA3.h`:
  `Kalman3D_Wave_CARMA3<T, NHC, SH, NZC, SZ>`; `NHC = 0` is OU-III's
  horizontal chain exactly as in `Kalman3D_Wave_CARMA` (bit-for-bit
  equivalent on a tilted, rotating, noisy stream: every state and covariance
  difference 0; and in the complete-INS benchmark, Q41 = Z41 on every metric),
  `NHC = 3, 4, 5` a CARMA block per horizontal axis `[p, v, a, latent (, S)]`
  sharing `A_h`, with the Kronecker `Q`.  Attitude, gyro bias, accelerometer
  and magnetometer models, the residual accelerometer-bias OU and its freeze,
  Joseph updates and resets are those of `Kalman3D_Wave_CARMA` (OU-III's
  mathematics).  Not carried over: wind heel, lever arm, the commanded-parameter
  random walk (none enabled by the deployed wrapper).
* `tests/kalman_carma/carma3/SeaStateFusionFilter_CARMA3.h`: the PR #704
  wrapper schedule (shared front end, proxy startup, staged commit, watchdog,
  magnetometer-service bias gate, Racc inflation) plus the horizontal mapping
  `kappa_h = f_tune / f_z,h`, `sigma_h = S_factor sigma_a,B` (no horizontal
  statistic is measured), the horizontal S laws of section 4.  Mappings fixed
  before any Monte Carlo run.
* `tests/kalman_carma/carma3-test.cpp` (36 checks): layouts; bit equivalence
  with the PR #704 hybrid; Kronecker discretization; stationary covariance
  invariant over 100 s of float prediction (<= 1.5e-4); S growth = `S_p(0)`
  per second with a DC level, bounded with a DC zero; float vs double
  discretization over f_z 0.05-1 Hz; noiseless three-axis sinusoid recovered
  within 10 % on every axis with PSD covariance; correlated prior; wrapper on
  a still IMU.  Wired into `make all` through `run_tests.sh`, with paired
  smoke runs of `carma3_bench` at both levels.

## 9. Benchmark

**Design.**  `tests/kalman_carma/carma3_bench` runs one scenario and seed per
process: a random-phase directional sea (cos^2s spreading; waves 30 deg off
the bow unless stated) through the vessel RAOs of section 2, plus the
scenario's low-frequency translation.  One IMU/magnetometer stream (200 Hz;
accelerometer white noise 1.5 mg, turn-on bias up to 8 mg with a random walk;
gyro white noise 0.0016 rad/s, bias up to 0.05 deg/s; magnetometer 0.8 uT at
25 Hz with hard/soft-iron terms) is generated once per seed and fed
identically to every filter: every comparison is paired.  Errors are taken
against the wave coordinate in the geographic frame (the INS is magnetic;
estimates are rotated by the declination before scoring) over t >= 300 s of
a 1200 s run; startup metrics use the first 120 s after Live.  No ground
truth enters any estimator or tuner.  Two levels:

* **complete INS** (`full`): OU3 is `SeaStateFusion_OU_III` unchanged, Z41 and
  Z52 are PR #704's `SeaStateFusion_CARMA`, the three-axis variants
  `SeaStateFusion_CARMA3` (same schedule, startup, magnetic acquisition,
  tuner, bias handling);
* **level 1** (`l1`): prescribed attitude: per-axis reduced filters on
  `y = R_nb f + g e_z` with the true `R_nb` and the same tuner, which removes
  attitude and gyro coupling and isolates the motion prior.

100 paired seeds on the four principal seas at both levels; 30 on 27
secondary complete-INS conditions and 5 level-1 ones; 20 seeds of 3600 s on
JONSWAP and the mixed sea; RUNS_TOTAL filter runs in `raw_runs.csv`.  All
tables, per-axis metrics, confidence intervals and Wilcoxon tests are in the
generated `reports/results/carma3_benchmark/summary.md` and
`paired_stats.json`; the tables below are excerpts.

### 9.1 Principal seas (complete INS, 100 paired seeds)

HEADLINE_TABLE

![Paired improvement, horizontal](../reports/results/carma3_benchmark/paired_improvement_horizontal.svg)
![Paired improvement, heave](../reports/results/carma3_benchmark/paired_improvement_heave.svg)
![RMS distribution](../reports/results/carma3_benchmark/rms_distribution_full.svg)

### 9.2 All conditions

ALL_CONDITIONS_TABLE

* **CARMA(4,1) horizontally is better than OU-III on every condition**: 43-73 %
  on every single-peaked, bimodal, broadband, directional, underway, anchored,
  drifting, manoeuvring, low-amplitude, sensor-fault and long case, 17-18 % on
  the mixed sea (also with the legacy horizontal truth and over 3600 s).
  Heading changes it by a few per cent; the following sea at 2.5 m/s, whose
  encounter frequencies fold towards zero, is where OU-III is weakest
  (10.5 % H_s) and C41 gains most (4.0 %).
* **CARMA(5,2) horizontally** is better than C41 on single-peaked seas
  (JONSWAP 2.05 against 2.66 % H_s) and much better under low-frequency
  translation (manoeuvring 8.0 against 14.1, drifting 3.4 against 5.4), level
  with it on broadband and abrupt-change seas, worse on bimodal (9.2 against
  6.2) and PM 0.15 Hz, and **71 % worse than OU-III on the mixed sea**.
* **Heave** depends only on the vertical model: the horizontal model changes
  heave by at most 0.02 % H_s (Z41 against A41).  C41 z is better than C52 z
  on every wave condition (the PR #704 result, here with 30-100 seeds per
  condition); both lose to OU-III on the mixed sea (-13 %, -163 %) and C52 z
  also on the bimodal sea (-34 %).
* **Mixed orders.**  C52 x/y + C41 z (M54) combines the best single-peaked
  horizontal and heave results and inherits C52's horizontal failure on the
  mixed sea (M54_MIXED); C41 x/y + C52 z (M45) is dominated by C41 x/y/z.

### 9.3 Horizontal S and shape variants

S_TABLE

* **CARMA(5,2) + S (derived law, A52S)** is within a few per cent of
  CARMA(5,2) everywhere (JONSWAP 2.02 against 2.05, PM 2.77 against 2.68,
  mixed 31.8 against 31.7 % H_s; still IMU 80 against 87 mm): the drift
  control comes from the prior, and S does not repair the mixed-sea failure.
* **CARMA(4,1) + S with OU-III's r_S (A41S)**, the mathematically invalid
  carry-over of section 4, is better than C41 where the error is low-frequency
  (still IMU 42 against 84 mm, below OU-III's 53; manoeuvring 12.7 against
  14.1; drifting 4.9 against 5.4; following underway 3.6 against 4.0) and worse
  where the sea itself has broad or low-frequency energy (mixed 17.5 against
  15.1; PM 3.10 against 2.91; bimodal 6.9 against 6.2; broadband 8.3 against
  7.3; JONSWAP 0.10 Hz 5.2 against 4.8), and degrades heading in the
  magnetometer gap (yaw 3.36 against 2.09 deg).  It acts as extra low-frequency
  shrinkage that counts prior information twice.  The effect it captures is
  real; a valid way to obtain it would have to be part of the prior, and is
  not tested here.
* **Sway-fitted shape (H41; hull-dependent).**  The C41 shape fitted to this
  benchmark's own sway RAO is better than the JONSWAP shape on single-peaked
  seas (JONSWAP 2.42 against 2.66; following underway 2.78 against 4.03) and
  worse on the mixed and bimodal seas (16.4 against 15.1; 6.4 against 6.2).
  It shows that the horizontal shape is worth about 10 %, with a sign that
  depends on the hull and the sea; it is not adopted.

### 9.4 Still IMU

STILL_TABLE

Every CARMA variant is worse than OU-III on a still IMU, horizontally
(84-87 mm against 53 mm) and vertically (46-53 mm against 27 mm, the PR #704
result).  Calibrating the turn-on bias changes nothing (`still IMU,
calibrated bias` rows), and only 20 mm of C41's 84 mm is below 1/60 Hz: the
phantom is mostly an oscillation in the prior's own band (band table below),
not a drift.  The only variant below OU-III is the invalid A41S.

### 9.5 Level 1 (prescribed attitude) and where the error lives

L1_TABLE

With the attitude prescribed, every filter's horizontal error drops (OU-III
5.18 -> 3.70 % H_s on JONSWAP) and the ranking stays: the CARMA advantage is
not an artifact of the attitude coupling.  CARMA(3,0) is no better than
OU-III horizontally and 56-66 % worse in heave; it is not considered further.
The error by frequency band (seed 1, Welch over t >= 300 s) shows the
mechanism:

BAND_TABLES

Horizontally, OU-III's error is spread from 0 to 0.4 Hz; C41 cuts it in the
wave band (0.2-0.4 Hz) and below it, C52 cuts it most below the wave band
(0-0.05 Hz: bias, tilt and rotation leakage, which `g dtheta` makes larger
horizontally than vertically) at a slightly higher in-band error.  On the
mixed sea, C52's error sits at the swell (BAND_MIXED); on the still IMU the
CARMA error sits in the prior's band (BAND_STILL).

![Time series, JONSWAP](../reports/results/carma3_benchmark/timeseries_xyz_jonswap.svg)
![Error PSD, JONSWAP](../reports/results/carma3_benchmark/error_psd_xyz_jonswap.svg)
![Mixed sea](../reports/results/carma3_benchmark/mixed_xyz.svg)
![Still IMU](../reports/results/carma3_benchmark/stationary_xyz.svg)
![Manoeuvring](../reports/results/carma3_benchmark/maneuver_xyz.svg)

### 9.6 Consistency, startup, convergence, attitude

CONSISTENCY_TABLE

* **NEES** (6-dof position-velocity, 6 if consistent): on single-peaked seas
  CARMA is closer to consistent than OU-III, whose horizontal axes are
  overconfident (NEES 3.2-4.6 per horizontal axis); on the mixed sea every filter is
  overconfident and CARMA much more (C41 64, C52 324, OU-III 31).  The
  accelerometer NIS (3-dof, 3 if consistent) is 0.06-0.56 for every filter:
  the wrapper's accelerometer noise (with its inflation) far exceeds the
  realized innovations, more so with CARMA's better acceleration prediction,
  so NIS does not test the motion prior here.
* **Startup** (first 120 s after Live, JONSWAP): heave RMS OU-III 696 mm,
  C41 263 mm, C52 138 mm; horizontal x 51 / 20 / 18 mm (peaks 266 / 95 /
  100 mm), against 15-25 mm in steady state.  Every filter has a startup
  transient; CARMA's is smaller because it starts from its stationary
  covariance, OU-III from a large non-stationary one.
* **Convergence** (corrected metric: time after Live from which the 30 s
  rolling RMS stays within 1.5x the steady value; median): heave 126 / 96 /
  92 s; x 692 / 740 / 85 s.  The long x values are late low-frequency
  excursions crossing the threshold, which C52 suppresses; on the mixed sea
  the medians are 860-1030 s (or never) for every filter.
* **Attitude.**  Roll and pitch errors are similar for every filter
  (0.14-0.39 deg RMS over the conditions, CARMA slightly lower in roll),
  except C52 on all axes on the mixed sea (pitch 1.60 deg, roll 0.48 deg: the
  unexplained swell acceleration is absorbed as tilt) and, less, the mixed
  orders there.  The horizontal CARMA prior
  improves heading where the magnetometer is absent or the vessel turns: yaw
  RMS in the 300 s magnetometer gap 2.90 deg (OU-III) against 2.09 (C41, C52),
  manoeuvring 4.36 against 2.75 (C41) and 2.14 (C52).

### 9.7 Coupled x-y prior with the oracle correlation

COUPLED_PARAGRAPH  This exceeds the 0.7-1.6 % noncausal bound of
section 5, which assumed the true spectrum; with a mismatched prior the
correlation is an additional constraint.  The correlation came from the true
directional spectrum; no estimator of it was implemented or tested, so this
is an upper-bound diagnostic, not a candidate.

### 9.8 Pinned vessel-RAO dataset

DATASET_SECTION

## 10. The known CARMA problems

Causes are numbered as in the study brief: (1) spectral fitting, (2)
transfer-function zeros, (3) noise-model mismatch, (4) acceleration
observation sensitivity, (5) adaptation of the characteristic frequency,
(6) pseudo-measurement design, (7) covariance modelling, (8) observability,
(9) numerical implementation, (10) an overly restrictive model of the vessel
dynamics.

| Problem | Reproduced? | Evidence | Causes |
|---|---|---|---|
| C41 better than C52 in the complete INS | Heave: yes, on every wave condition.  Horizontal: **reversed** on single-peaked seas | Heave JONSWAP 1.83 / 2.04, PM 2.10 / 2.68 % H_s; identical at level 1 and in the white-noise steady state (10.2 / 13.7 mm): not attitude, bias or adaptation.  Horizontally C52 rejects sub-band error better (section 9.5) | 1, 2, 10 (heave); 2, 4 (horizontal reversal) |
| C52 collapse on the mixed sea | Yes, worse, and on all axes | Heave -163 %, horizontal -71 % against OU-III; 3600 s and legacy truth the same.  White-noise steady state at the oracle operating point (x): 284 mm signal distortion against 3.8 mm noise, so it is the prior's shape.  C52's numerator zeros sit at 0.151 Hz x kappa (zeta 0.2); the frequency mapping puts the mixed sea at f_z = 0.23 Hz (section 6), kappa = 0.62, so the notch falls at 0.094 Hz, on the swell peak (0.09 Hz) that carries about 60 % of the elevation variance.  The filter cannot represent the swell; the unexplained swell acceleration leaks into tilt (pitch 1.6 deg) | 1, 2, 5, 10; consequence 7 |
| Still-water phantom | Yes, vertical (C41 53, C52 46 mm against 27) and **also horizontal** (84-87 against 53 mm) | Unchanged by turn-on bias calibration; 20 of C41's 84 mm below 1/60 Hz, the rest in the prior's band.  In still water the wrapper keeps a nonzero operating point (amplitude floor 0.05 m/s^2, as for OU-III; the tuner's frequency bounds), and the resonant CARMA prior turns noise and slow tilt/bias leakage into motion at its resonance; OU-III's S = 0 is a broad, low-Q DC-zero prior (section 4) and a double-counted S (A41S) halves the CARMA phantom | 5, 6, 3, 10 |
| Overconfidence under mismatch | Yes on multi-peaked seas; **not** on single-peaked ones | NEES 6-dof mixed: OU-III 31, C41 64, C52 324; JONSWAP 7.7 / 5.4 / 4.4 (6 = consistent).  The Kalman covariance is conditional on the prior; with a prior outside the true family it is not an error bound, and the narrower the prior the larger the miss | 7 (as a consequence of 10); no model-error term is carried |
| Large startup transients | Present for every filter, **smaller** for CARMA than for OU-III | First 120 s after Live, JONSWAP: heave 696 (OU-III) / 263 (C41) / 138 (C52) mm RMS, against 15-25 mm in steady state; horizontal 51 / 20 / 18 mm.  CARMA starts at its stationary covariance; the remaining transient is the shared startup (proxy, magnetic acquisition, tuner commits) | 5 (operating point still settling) |
| `conv_s = 0` in PR #704 | Yes: metric defect | NaN before Live poisoned the rolling sum; fixed in `carma_bench.cpp` and in `carma3_bench` (section 1) | 9 (in the benchmark, not the filter) |
| Inconsistent frequency ranges | Yes, documented | Amplitude mapping matches total model variance to an in-band RMS; process-validation PSD compares acceleration on the displacement band (section 1) | 1 (mapping), not limiting |

Observability (8) is not a cause of any of them: every CARMA prior is
observable from acceleration alone, and the one non-stabilizable construction
(S on a DC-zero model) is not a candidate.  Numerical implementation (9) of
the filters is not a cause: float discretization, stationary-covariance
invariance and the bit-exact OU-III-horizontal equivalence are tested.

## 11. The fundamental tradeoff

1. **Is a stationary CARMA displacement model a better representation of
   horizontal marine motion than the integrated OU-III model?**  For the
   oscillatory, wave-induced part, yes in the sense that matters to the
   filter.  Wave-induced horizontal displacement is stationary with a
   spectrum that vanishes at DC (section 2).  A DC-zero CARMA(4,1)
   represents that directly; OU-III represents it only after conditioning on
   `S = 0`, which makes it an equivalent broad, low-Q DC-zero CARMA(4,1)
   (section 4).  With a fitted shape the horizontal error halves on every
   single-peaked sea.  The shape is the limit: with the white-noise steady
   state on the mixed sea, the fitted C41 is worse than OU-III's broad prior
   (88 against 70 mm, section 6).  In the complete INS, C41 is still better
   (+17 %), because CARMA rejects the low-frequency leakage that dominates
   the INS error.
2. **Can CARMA naturally suppress horizontal displacement drift?**  Yes.  A
   Hurwitz prior bounds the displacement, and with a DC zero (C41) also its
   primitive, without any pseudo-measurement.  Low-frequency error (60 s
   mean), x, JONSWAP: OU-III 7.6 mm, C41 3.9, C52 2.7; drifting 12.5 / 6.4 /
   4.0; manoeuvring 28.8 / 15.0 / 9.4.  The price is that the prior shrinks
   all low-frequency content, real swell included (the mixed sea).
3. **How does it distinguish oscillatory wave motion from low-frequency
   translation?**  It does not detect anything; it shrinks by frequency.
   Translation accelerations (1-30 mg, below 0.02 Hz) are indistinguishable
   in an IMU from accelerometer bias and `g dtheta`, so the prior attributes
   them to the bias/tilt states or to nothing.  The wave coordinate stays
   bounded, but its error under translation is 1.3-5 times that of the same
   sea without translation (section 9.2).  C52's real pole and zeros reject more of it than C41.
   No model here estimates the vessel's translation; that needs another
   sensor (GNSS).
4. **Would independent CARMA axes lose useful directional information?**  A
   bounded amount.  The noncausal bound with the true cross-spectrum is 5-17 %
   on x (single-peaked to mixed, white noise).  Most of it is x-z/y-z
   quadrature coherence, which needs wave direction (with its 180 deg
   ambiguity), spreading and the hull's RAO phase.  The x-y part is <= 1.6 %
   in the bound and 5 % in the complete INS with the oracle correlation.  A
   physically coupled x-z model with oracle coupling lost to independent
   axes on single-peaked seas (section 5).
5. **Could correlated process noise improve the three-axis model?**  For
   x-y, yes, a little, at zero state or compute cost (`Q_H = Sigma_h (x)
   Q_1`).  The x-y co-spectrum `cos th sin th` is invariant under the
   180 deg direction ambiguity, so its correlation could in principle be
   estimated in band from horizontal accelerations; that estimator is
   untested.  For x-z and y-z, no: correlated noise through identical
   dynamics has a real cross-spectrum, and the physical relation is
   quadrature.
6. **Is CARMA(4,1) preferable to CARMA(5,2) for horizontal axes?**  As a
   fixed choice, yes.  C52 is 8-26 % better than C41 on most single-peaked
   seas (worse on PM 0.15 Hz) and 27-43 % better under translation.  C41 is far better on the mixed sea
   (+17 % against -71 % relative to OU-III) and on the bimodal sea.  The sea
   state is not known in advance, C52's failure is catastrophic and C41's
   shortfall moderate.  Mathematically, C41 has the DC zero that the physics
   and the bounded-primitive premise ask for; C52's integral is a random
   walk.
7. **Is a mixed-order model mathematically preferable?**  No.  The
   wave-induced horizontal and vertical spectra of a free-floating hull are
   both the incident spectrum times a smooth RAO, so there is no physical
   reason for different orders.  The empirically best single-peaked
   combination (C52 x/y + C41 z) inherits C52's horizontal failure on the
   mixed sea; C41 x/y + C52 z is dominated.
8. **What additional observability or stability assumptions are required?**
   For the motion model, none: every CARMA prior is Hurwitz and observable
   from acceleration (section 3).  An `S` state requires `B(0) != 0` for the
   augmented pair to be stabilizable, so it must not be added to a DC-zero
   model.  For the INS, yaw observability (magnetic service) and
   bias/tilt/motion separability remain as for OU-III, but the CARMA prior
   changes how low-frequency acceleration is split among them.  The
   operating-point commits make the filter a switched linear system.  The
   premises and proofs of the OU-III stability study (MARINE MOTION, IMU
   BIAS, MAGNETIC SERVICE) are about the shipping OU-III filter; none of them
   transfers, and **nothing here is a stability statement**.
9. **What happens when the vessel's motion is not in the assumed CARMA
   family?**  Everything in this study is that case: the truth contains no
   CARMA model.  The error is then dominated by signal distortion rather than
   noise (section 6), the covariance stops being an error bound (NEES), and
   failures are sea-specific.  The more specific the prior's shape (C52's
   notch), the larger the failure on seas it does not fit.  OU-III's broad
   prior is the most robust on the mixed sea in heave.
10. **Does any measured improvement justify the increase in state dimension
    and computational cost?**  For C41 x/y/z there is no increase: 21 states
    like OU-III, 25 % fewer multiply-adds per sample, less code, a faster
    host update at a fixed sample interval (section 7).  The cost it does add
    is per-sample re-discretization when the interval jitters.  For C52 x/y/z
    (+3 states, OU-III's cost), the single-peaked gain over C41 does not
    justify the mixed-sea risk.  S states (+2) bring no consistent gain.
    Coupling adds no states but needs information the filter does not have.

## 12. Recommendation

**B: CARMA(4,1) on X/Y and CARMA(4,1) on Z**, as the research candidate.
**Not shippable.**

* Against **A** (OU-III horizontal, CARMA z): horizontal CARMA(4,1) is better
  than OU-III on every one of the 31 complete-INS sea conditions (43-73 %;
  mixed sea 17 %), at level 1 and under translation, with no extra states or
  cost; heave is the same as A.
* Against **C** (CARMA(5,2) everywhere): the mixed-sea collapse (horizontal
  -71 %, heave -163 %) and the bimodal heave loss.  Its mechanism (a notch
  that the frequency mapping can place on a swell) is structural, not a
  tuning issue.
* Against **D** (different orders or parameters per axis): the per-axis
  choice that wins on single-peaked seas (C52 x/y + C41 z) fails on the
  mixed sea; the sway-fitted C41 shape is hull-dependent and mixed in sign.
* Against **E** (coupled): the oracle gain is 5 % (x-y), a physically
  coupled x-z prior with oracle coupling lost to independent axes, and no
  estimator of the coupling exists.
* Against **F** (OU-III everywhere): the paired gains are large, consistent
  over 31 conditions and two levels, with tight intervals, at lower cost.

**Remaining obstacles before B could be considered for deployment**

1. Still water and very quiet seas: horizontal phantom 84 against 53 mm,
   heave 53 against 27 mm.  The cause is the resonant prior at a floored
   operating point; a valid low-frequency or "no waves" mechanism in the
   prior is needed (A41S shows the size of the effect, not a valid fix).
2. Mixed swell and wind sea: heave 13 % worse than OU-III (PR #704) and
   NEES twice OU-III's; a single fixed shape scaled by one frequency cannot
   represent two peaks.
3. Covariance under mismatch: no model-error term; the covariance is not an
   error bound on multi-peaked seas.
4. Sample-clock jitter re-runs the matrix exponential each sample; quantize
   dt or tabulate `F`, `Q` with a stated error before any target timing.
5. Hull and sea dependence: one synthetic hull and one replay hull; the
   horizontal shape is the JONSWAP fit; no lever arm, wind heel or
   nonlinear waves.
6. No hardware evidence: no AtomS3R timing, float behaviour or sea data.
7. No stability analysis: the OU-III stability study does not cover a CARMA
   motion model.

## Limitations

Synthetic truth (linear random-phase directional seas through explicit
RAOs, one hull) and one pinned replay dataset (one 28 ft hull, single-peaked
seas); the low-frequency translation models are simple; the IMU is at the
centre of gravity.  All results are simulation results; the analytic results
(sections 3-6) are exact for their stated models; the oracle diagnostics
(coupling, true-amplitude prior) use information the filter does not have and
are labelled.  H41 and every horizontal-RAO-dependent number hold for this
benchmark's hull only.

## Reproduction

    make -C tests/kalman_carma build test
    python3 tools/carma3_analysis.py                  # sections 2-7, figures
    python3 tools/carma3_benchmark.py run --jobs 4    # synthetic grid (resumable)
    make ensure-sim-data && make -C tests/kalman_ou_iii kalman_ou_iii-sim
    python3 tools/carma3_benchmark.py dataset          # pinned dataset, 10 seeds
    python3 tools/carma3_benchmark.py traces
    python3 tools/carma3_benchmark.py truth
    python3 tools/carma3_benchmark.py report           # summary.md, figures
