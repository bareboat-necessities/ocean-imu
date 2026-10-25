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

RESULT_PARAGRAPH

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
mismatched.  x-y correlation alone can give at most 0.5-3 % (and is tested
once in the complete INS with the oracle correlation, section 9).

**Conclusion on coupling.**  A physically meaningful cross-spectral model
needs direction (with its 180 deg ambiguity), spreading and the hull's RAO
phase, all sea-state dependent; even with all three given exactly, the
achievable causal coupling structure loses to independent axes where CARMA
works, and the gain bound is small there.  A coupled model is not practical
for the real-time filter at this stage; correlated noise across identical x/y
dynamics is practical (zero cost) but can represent only the co-spectrum and
has a negligible bound.

## 6. Steady-state error budgets (white noise only, oracle operating point)

Exact steady-state Kalman transfers (DARE at 200 Hz) applied to the true
displacement spectra; error = signal distortion / noise propagation, mm RMS.
Horizontal priors use the wrapper's mapping (frequency from the vertical
channel, amplitude = vertical band acceleration RMS, i.e. OU-III's isotropic
choice).

STEADY_STATE_TABLE

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

COST_TABLE

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

BENCHMARK_SECTION

## Reproduction

    make -C tests/kalman_carma build test
    python3 tools/carma3_analysis.py                  # sections 2-7, figures
    python3 tools/carma3_benchmark.py run --jobs 4    # synthetic grid (resumable)
    make ensure-sim-data && make -C tests/kalman_ou_iii kalman_ou_iii-sim
    python3 tools/carma3_benchmark.py dataset          # pinned dataset, 10 seeds
    python3 tools/carma3_benchmark.py traces
    python3 tools/carma3_benchmark.py truth
    python3 tools/carma3_benchmark.py report           # summary.md, figures
