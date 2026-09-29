# OU-III proof research state

## Current hypothesis

One persistent MARINE MOTION / IMU BIAS / MAGNETIC SERVICE execution follows
construction -> capture -> informed H18 -> refinement/release -> recurring
A21 -> regional practical stability. STILL / TRANSITION / MOVING qualify
portions of this same route. No shipping mode switch is enabled.

Excitation means Delta_g([t,t+T_E])>=theta_E on **every complete window
contained in one maximal physical moving episode**, between nondegenerate
rest intervals. Isolated zero-rate instants do not split an episode. Crossing
windows do not owe excitation. T_E/theta_E remain symbolic. The old global
rolling stillness-or-span rule is removed, not repaired by more assumptions.

The stationary target must allow physical attitude/BA ambiguity and cover
all measurement-compatible histories before quiet evidence changes estimator
behavior. The existing assumptions admit exact rest/motion indistinguishability;
a sound exact-rest detector with finite entry and universal finite exit is
impossible on that admitted parameter family. This does not establish that
no estimator change is needed. Full stationary practical boundedness is OPEN.
The moving signed margins Delta_col/Delta_gyr, six-pivot uniform premises,
B_*, J_AG, rho_0 and nonlinear every-prefix retention remain OPEN.
The exact zero-residual quiet nominal subfamily now has a root-independent
historical action/full covariance ceiling and qualitative homogeneous linear
loss.

In world coordinates every historical AG row is `-R_k[f_k]x[A~_k R_0',B~_k]`
(`ou3-world-frame-rows.md`): the attitude estimate and its error leave the
six-column geometry, except through world injections and the nominal
rotation integral of the gyro columns. Same-cell geometry also depends on the
applied magnetic cadence, so the MOVING premise is pursued as an **aggregate
world-frame** row statement. The attitude columns need only the nominal
signed mean of the AW state (Corollary A*: transverse mean below 1.96133
m/s^2); the pointwise physical tracking premise is false on admitted
histories. MAGNETIC SERVICE on every 1-s interval plus the implemented
nominal-rate bound make the field-axis coordinate of the gyro curve monotone
(Lemma T), which yields an explicit injection-free aggregate six-column floor
(Theorem G0) from two separated nominal accelerometer windows. OPEN: a
source bound on the nominal window statistics and G0 with injections.

**Contraction.** The complete A21 word contraction is one
information-ratio inequality. Its `k=0` form is the word's own Riccati
diameter `kappa_W=lambda_max(J^-1 A)=lambda_max(Pi^-1 P_diff)`, where `Pi` and
`P_diff` are the known-root and diffuse-root terminal covariances. It gives
`rho_W<=tanh(log(kappa_W)/4)` for every root covariance, and the bound is
sharp. On every carried word the only direction the data miss is the
one-dimensional physical tilt/BA kernel `nu` (tilt about the body field axis
with compensating BA). Along it Corollary K needs only the scalar
`nu'P_0 nu`, whose BA part is the proved `P_ba<=I/1600`. The remaining
source-uniform obligations (O1, O2 in `ou3-corrected-word-proof.md` §7)
therefore need no `B_*`, `I_eff` or excitation premise; `B_*`/`I_eff` remain
for coercivity.
Every first-prediction relative process comparison is capped at
`3.3741e-10` per prediction and is abandoned (DEAD_END 25).

**Historical reader.** It is the exact joint minimum-action reader
`B*=Pi+Tt I_eff^-1 Tt'`, the diffuse-AG-root Riccati limit with a full 21×21
domination. Carried windows of 16–64 s make it commensurate (≤12× the actual
AG covariance, ≤5× at 64 s). It now serves coercivity only (DEAD_END 27).

## Evidence

- Shared OU arithmetic uses cancellation-safe dimensionless SO(3) integral
  coefficients and alias-safe covariance symmetrization. The PSD checks use
  LDLT only for nonnegative-pivot acceptance; negative or failed pivots defer
  to the eigenvalues before applying the existing roundoff tolerance. The
  float/double regression retains zero-rate and branch-boundary checks,
  singular PSD invariance, indefinite repair and correlated negative blocks
  across three scales. These are implementation checks, not a finite-error
  contraction certificate. Real-arithmetic gyro-radius and reset lemmas are
  unchanged; native source diagnostics must be regenerated, not restamped.

- `ou3-regime-design.md` states observable stationary information, required
  practical theorem, detector obligations, carried transition state and the
  obstruction before any shipping change. Terminal gyro-bias averaging has
  exact supply N_g+Omega_tube+D_g sum(weight*sample_age); coherent deterministic
  noise does not average away. Gravity/BA retain a one-dimensional kernel
  about the known magnetic axis. Existing projections alone do not close it.
- The sin-cubed rest/motion/rest witness uses alpha=1/1000, nu=1/40,
  p=v=a=0, B=75e_x, b_g=-phi' e_x and b_a=g(R'e_z-e_z). Its bias/rate
  upper bounds are .00980665, .00073549875, .000075 and .000005625 in
  their respective SI units. All lie strictly inside unchanged limits;
  joins are C2, jerk/primitive vanish and complete 80*pi windows have span
  .002. The measured packets are exactly level/north for arbitrary duration.
  Existing actually-applied stationary service retains a floor >1 after
  the exact rational lower factor (1-alpha^2/2)^2. This is a witness parameter
  instance, not a deployment qualification of symbolic T_E/theta_E.
- `regimes.py` implements complete-window requirements, exact conditional
  every-prefix bridge composition and deterministic stationary gyro supply.
  Its constant-memory necessary-evidence monitor uses physical/noise bounds,
  a supplied continuous dwell and immediate evidence-failure exit. Invalid
  packets/gaps restart dwell. It never returns a physical-STILL certificate.
  Estimated wave states cannot distinguish identical input histories.
- `ou3-moving-six-pivots.md` continues the same-history reader: chronological
  prediction/reset recursion bounds |B-tI|, then the actual two-group row
  budget gives s=c/[1+(a+1)/b0]-epsilon. If s>0, every one of six greedy
  row pivots is >=s/sqrt(m). All asynchronous/reset/reference defects remain.
  Applied acc/mag groups within a prediction cell now give E=0 exactly by
  retaining every intervening reset in C_i. No simultaneous-event assumption
  is made. Source-uniform c,b0 and factor/count bounds remain unproved.
  The exact supplied matrix audit has s=19/100; its 80-digit singular value
  is .35070823886858529042926658460119746953531050023679. This is algebraic
  feasibility, not a contraction construction or a shipping certificate.
- The exact zero-residual quiet nominal record supplies a special two-group
  floor: actually applied acc/mag pairs eight qualified predictions apart
  have C'C>=81 I, A=I, B=(sum h_i)I with sum h_i>=4/125, and E=0.
  Thus s=18/127 and six pivots >=9/254. This is a
  nominal row theorem, compatible with physical tilt/BA ambiguity. It is
  not a physical identifiability or nonlinear theorem. The same historical
  reader theta_read=z1, bg_read=(z1-z0)/T now cancels arbitrary AG root and
  cross covariance. Inherited AW/BA, correlated process and actual observation
  bounds give B_q, an every-operation full upper C_q, and qualitative
  homogeneous linear loss on this exact subcase. The full compatible-class
  nonlinear theorem remains OPEN; covariance retention is not state retention.
- Literal reset G=I+[d]/2 and ideal Rodrigues prediction have inverse norm
  <=1; the shipping coefficient-series transfer remains qualified by
  DEAD_END 24. C=A^-1 B is unchanged by resets; later gyro injections retain A^-1.
  Exact prefix recurrences bound |C-tI| without amplifying its past defect
  by every reset. A terminal d=4e_x no longer erases an earlier .005-s floor.
  The 80-digit supplied noncommuting audit has floor .0396561974478783
  versus singular value .0399999504651045. The quiet correlated-root action
  audit ratio is .413929501952256 at root scales 1 and 10^12. These algebra
  audits precede enclosure and are non-promoting.
- A fresh carried-source audit reuses the existing read-only observer and
  unchanged driver, with exact untapped terminal parity on quiet and moving
  inputs. Each has eight complete same-cell groups. First/last group E=0
  is exact over exported rational operators; 80-digit six-column budgets
  are 1.204325434846729 and 1.074485913647320, below actual singular values
  1.922502992341225 and 1.717587283388896. These finite numbers are neither
  uniform certificates nor real-trajectory or magnetic-service enclosures.
- World-frame factorization (`world_frame.py`): predictions leave
  `A~=R'AR_0` invariant when the covariance rotation matches the mean
  increment, resets multiply it by `N=R+'GR` with
  `N'N=I+(|x|^2 I-xx')/4`, and `B~` integrates the nominal R_hat'. Exact on
  supplied rational words. Same-cell sigma(C) equals that of
  `[[f]x;[B_w]x N]`: world nominal force, reference and local injection only.
  The exact projector eigenvalue `(F+B)/2-sqrt((F-B)^2/4+FB kappa^2)` and the
  rational reset-angle bound `theta/2+theta^3/24+theta^2/(16-theta^2)` give c.
  The quiet floor extends to any constant attitude and admitted reference:
  `c^2>=1296/481`, s=16/635, six pivots >=4/635.
- Literal injections obey `dd'<=NIS K_theta S K_theta'<=NIS P_theta,theta`
  (Cauchy--Schwarz, then Joseph with `S-HPH'>=0`). The bound enters a, b0 and
  the `|d|^2|v|/6` reset supply.
- Exact collinear MARINE MOTION/IMU BIAS history: `a=c cos(2 pi t)`,
  `c=g(e_z-(e_z.b)b)`, `B=(21,0,72)` uT, makes force parallel to B at every
  integer second. With magnetic corrections only then, a carried shipping run
  gives same-cell sigma .0301--.0356 and sine >=.0032 but aggregate
  six-column sigma 28.50 over 4 s; AW tracking error reached .563 m/s^2. That
  cadence fails MAGNETIC SERVICE: one-correction 1-s windows have least
  service eigenvalue 2.03e-5 against mu_M=1. Jerk forbids force/field
  collinearity at every instant of a cadence with length-weighted mean gap
  below `4(g h-2V/L)/J` (.051 s at h=1/5, L=16 s); tent dips attain it.
- Carried world-frame audit (quiet, wave, collinear at 1 Hz and at 25 Hz;
  control parity): row factorization 1.7e-77, reset Gram 6.4e-81, prediction
  branch 2.1e-9. World floors/actual <=1 (quiet 1.0) once each reset charges
  its float mean-injection angle (<=1.8e-9 rad) on top of psi(theta); without
  that charge the 25-Hz word's near-collinear group exceeded its actual value
  by 1.7e-7 relative. Injection/NIS-prior ratios <=.134. Wave two-group
  budget 1.0656 versus actual 1.7176. `|A~-I|` matches half the signed world
  injection sum (.00112 versus norm sum .0135, collinear word). CI reproduces
  the committed record exactly and verifies every reported metric.
- AW covariance ceiling (`aw_covariance_ceiling.py`, Lemma B): with S_factor=1
  the pending sync is the spectral max of P_aw and sigma^2 I, predictions keep
  `(1+eps)16`, and corrections subtract `K_a S K_a'`, so
  `lambda_max(P_aw)<=(1+eps)16` (156^2 before) with excess decaying as
  `exp(-t/6)`. Every applied sync floors P_aw at sigma^2 I, so the ceiling is
  tight; an exact S_factor=2 witness shows isotropy is necessary. Corollary
  A's storage radius becomes .280954 at the clamp (.0072 before), >=1 for
  `sigma_max<=1.1238`. Carried: AW reconstruction 1.7e-6, sync isotropy
  8.2e-7, step ratio <=1+8.2e-7, ceiling ratio <=.999982.
- Corollary A: on 16-s windows the normalized nominal attitude-column Gram is
  `>=(1/5-e)^2/((1+e)^2+1)`, `e=(11/16+.15+epsilon_a)/g`, positive iff
  `epsilon_a<1.12383 m/s^2`; `gamma(0)` is about .00603. No attitude error or
  magnetic residual is charged.
- Corollary A* (`aw_tracking.py`, exact): the normalized attitude-column
  Gram is `>=(sigma_w-m_perp/g)^2/((1+m/g)^2+1)` for the nominal mean
  `mu=sum alpha_i a_hat_i`; positive iff `m_perp<1.96133 m/s^2`,
  `gamma*(0)=1/50`. The signed physical transfer recovers 1.12383 as a bound
  on `|sum alpha_i e_i|`, not on `sup|e|`. The nominal mean is a nonnegatively
  weighted signed sum of AW corrections; the world innovation factors as
  `r=R_hat(y-a_hat)`, and the AW loop obeys an exact gain-weighted identity.
- AW carried audit (`aw-tracking-source-feasibility.json`, six admitted
  C2-onset histories, exact envelopes, untapped-control parity): pointwise
  error up to 7.647 m/s^2 (6.80 x 1.12383) with attitude error .0062 rad;
  signed 16-s mean error <=.371 m/s^2 (ratio .330, sync-locked
  rectification); nominal transverse mean <=.348 (ratio .177); nominal L1
  force <=1.091. Rectification: Gamma is periodic in the 21-sample AW sync
  cycle, and a phase-locked jerk-limited triangle raises the signed mean
  from <=.05 to .371; a raised tuner sigma reduces it.
- Lemma I* (`signed_injection.py`): the ordered injection rotation J obeys
  `M_N M_0'=J K` exactly, so `angle(J)<=alpha_0+alpha_N+int|omega_tilde|`
  with no norm sum; one reset has `|x|<=alpha_-+alpha_+` (DEAD_END 17 needs
  attitude error >=1.1416 rad). `N=I-[x]/2+O(|x|^3)`. Corollary A** charges
  a signed mean of relative rotations.
- Lemma T and Theorem G0 (`aggregate_floor.py`, exact rationals): service
  gaps <=1 s and `Omega<=1.1508652` give `|b.T|>=cos(Omega)-Omega eps`
  (>=.4078-Omega eps); with m_perp<=2/5, u1<=6/5, L=16, G=64,
  `s^2>=1.486786e-3`, `q_I=16.81`. Synthetic falsification: actual
  sigma_min 37--60 versus floors .067--.070. Carried literal arrays (six
  100-s A21 words, every reset retained): sigma_min 191.8--267.4, at least the
  injection-free value, `|A~-I|<=.0045`, G0 floors .072--.084.
- The current excited moving forced-balance word already has roll
  .02 sin(t/2) and displacement .4 sin(.6t)e_z. It satisfies the moving span
  premise for T_E>=4*pi, 0<theta_E<=.04. Its 80-digit balance residual is
  about 3.38446651750e-82, with exact observer/control terminal parity.
  Its failed rotation and variation budgets remain failed after removing
  indefinite stillness; see DEAD_ENDS 13--14. No scalar refinement is retried.
- The exact carried-adjoint compatibility criterion is ker C subset ker W.
  Otherwise both residual sums remain. The forced-data identity retains
  root action, joint signed rotation/reference, literal gains, all Euclidean
  means and reset/projection/rounding defects. In the frozen-attitude map,
  F E_BG=E_BG and Z E_BG=0; that cancellation is not gyro feedback.
  The quiet word's exact seven-column null witness has |Wv|>
  754573637/5000000000000000, while its actual quiet innovations vanish.
- Joint minimum-action reader (`ag_readout.py`, `ou3-ag-readout-proof.md`):
  - **Construction.** Every action source is a unit column of one augmented
    chronological design `y=O_h h0+A s`: the nuisance-root factor, every
    process/sync factor and every applied noise factor.
  - **Exact results.** The Loewner-minimal feasible action is
    `B*=Pi+Tt I_eff^-1 Tt'`, and every feasible reader satisfies
    `B(L)=B*+(L-L*)Sigma(L-L*)'`. The diffuse-prior posterior equals the
    Riccati recursion from `diag(t I6,U_n)` at every scale, and
    `B* <= (1+1/g)TT'+(1+g)T_h I_eff^-1 T_h'`.
  - **Fixture checks.** All of these are exact on the supplied 21-state word,
    where the pivot reader's action is up to 4.3× larger.
  - **Carried replays** (real-arithmetic optimal-gain replay of the literal
    coefficients, corrected core): `B*` is within 6.3e8/1827/12.2/5.1× of
    the actual AG covariance on 0.32/4/16/64-s quiet windows, and 4.9e8/2560/
    12.4/5.4× on wave windows. The full 21×21 diffuse limit is within
    55–77× from 16 s on and dominates the literal terminal covariance.
- Contraction formulations (`corrected_word.py`, `word_energy.py`,
  `ou3-corrected-word-proof.md` §6):
  - **Exact identities.** `M=C_end'P_0^-1` and
    `M'P_end^-1 M=P_0^-1(Sigma_00|y-Sigma_00|y,x_end)P_0^-1`.
  - **Information-ratio lemma.** It is exact, with the proof in §6.
  - **First-prediction ceiling.** It is exact rational: `3.3741e-10` per
    prediction.
  - **Direct carried 0.32-s margins** (literal gains, corrected core):
    quiet `4.18722e-4` (unchanged from the pre-fix core), wave `6.929e-4`.
  - **Where the loss sits.** The slowest direction is translational (≈87%
    v/p, 10–13% BA). Its loss is 61–79% S pseudo-observation, 16–23% AW
    sync, 13% acc and 3% predictions.
  - **Other directions.** Bias directions have margins 1.7–2.2e-3.
    `lambda_min(D_NN)` is .26 quiet and 1.54 wave; the nuisance-eliminated AG
    Schur loss is ≥733 quiet and 1121 wave.
  - **Longer words.** The exact word margins are .024/.098/.39 at
    4/16/64 s.
  - **Separated bounds fail.** The information-only and forgetting-only
    bounds stay at 3e-26–3e-3 and 2e-10–6.6e-3 respectively.
- Word Riccati diameter (`word_diameter.py`, section 7 of
  `ou3-corrected-word-proof.md`):
  - **Exact.** Theorem D, the closed-form supremum and the sharp scalar word
    (`sup rho=3-2 sqrt2` at `p=1/sqrt2`) are checked exactly, as are
    composition, slow/fast factorization, compression/reader duals,
    Corollary K with invariance and the S-chain cancellation.
  - **Source-uniform cap.** Gyro-bias persistence gives
    `kappa_W>=sigma_g^2/(b0 T^2)` up to `1e-9` relative: 71.19, 17.80 and
    4.449 at 16, 32 and 64 s. No `k=0` certificate beats margin .212, .383
    or .643 there.
- Carried information-ratio feasibility (`information-ratio-source-feasibility.json`;
  float64 optimal-gain replays, 64-s history, 16/64-s words):
  - **Kill criterion.** The ideal `C=P_0` loses ≤1.3×. The joint-reader `C`
    loses at most 4.9× on MOVING words, so the criterion passes. On quiet
    words it loses 51.7× at the 289-s root, whose BA marginal is still an
    order of magnitude below its 5064-s level, but 5.7× at 5064 s.
  - **Diameter.** On MOVING words the joint-reader optimum is `k=0`. `kappa_W`
    is 19380/509 (wave), 570/37.5 (collinear) and 436/20.1 (sync-locked) at
    16/64 s; it is infinite on quiet words. The controlling direction is
    tilt/BA (67–93% BA share), except gyro bias on the 16-s sync-locked word.
    The fast-given-slow factor is ≤1.078.
  - **Kernel.** With the exact `nu'P_0 nu` the rank-one bound loses ≤1.9×. The
    ceiling `(sqrt(tau)+|nu_ba|/40)^2` with the actual tilt loses 1.2× on
    steady quiet words, and 1.5× with `tau=10^-3 rad^2`. The kernel set is
    invariant on every word, using the next root's kernel.
- Existing shipping projection evidence is retained unchanged: residual
  gyro-bias <=.5 rad/s, qualified prediction angle <.007, one-step transport
  floor .003999991833333333 s. A complete turn would require bias norm
  >1046.5466859583 rad/s on the qualified 4--6 ms family. None of these is
  a complete signed temporal gyro margin. No projection was active in the
  verified 840-case validation plus 18 calibration and 310-case robustness
  plus one calibration paired studies. Those prior source/input-bound
  replays do not certify the changed core; the existing full-evidence
  pipeline must regenerate the affected evidence, not restamp it.

## Current limiter

Exact physical rest is not identifiable from the current sensor/bias model.
A stationary practical theorem for the entire compatible class is required
before shipping detector/state/covariance changes can be justified. Conditional
finite bridge algebra does not supply detection liveness, a retained set or a
budget for arbitrarily repeated switches. These are explicit OPEN obligations.

On MOVING windows the six-column geometry is attitude-free. A same-cell
floor would need a coupling between applied magnetic cadence and the jerk
lemma; the aggregate premise avoids it. The injection-free aggregate floor
is now explicit (Theorem G0) under two nominal window statistics: the
transverse nominal AW mean (<1.96133 m/s^2 needed; carried worst .348) and
the L1 nominal force (carried 1.091). Controlling quantities, in order:
(i) a source bound on those nominal statistics from the literal AW loop,
whose average is gain-weighted and rectifies at the 21-sample sync cycle;
(ii) the injection frame `Q'b`, which the .02 rad/s gyro residual can rotate
by ~1 rad over a 100-s word, beyond the global fixed-b tube of G0;
(iii) the contraction factor, now one diameter of the word itself:
- **(O1)** a source-uniform ceiling on the kernel-bounded diameter
  `kappa_nu=lambda_max(Pi^-1 P_nu)`, i.e. observability of every root
  direction except the physical tilt/BA kernel;
- **(O2)** the scalar kernel ceiling `c_nu` and its invariance
  `nu_next'P_nu nu_next<=c_next`; the BA part is proved and a tilt ceiling
  about the body field axis of order `10^-3 rad^2` remains.

On carried 0.32-s words the slowest direction is translational, not AG.
Quiet water leaves the tilt/BA kernel about the magnetic axis outside `J`;
there contraction comes only through the scalar kernel variance and BA decay.
Norm-summed NIS/covariance injection bounds overcharge multi-second transport.
Quiet-subcase homogeneous decay does not control compatible physical mismatch.
Preserve actual chronological gains, resets, OU and bias histories. No
independent nominal boxes, unsigned energy, sampled Gramian, selected minor
or finite replay closes B_*; its exact remaining premise is `I_eff>=mu`, which
serves coercivity, not `rho_0`.
Every new lemma's role in V_next<=rho V+c_d|d|^2 is stated in the design:
stationary information sets a supply/ambiguity radius, bridge products and
supplies compose rho/c_d, and six pivots feed the historical covariance/loss
comparison. Source-uniform numerical contraction enclosure is not yet justified.

## Failed approaches / DEAD_ENDS

1. **Unsigned cumulative mean-action energy: DEAD_END.** E=817885.0623259853
   while collinearity cost is about 1.0138313892; exact E_col-E lies in
   [-817884.048495,-817884.048494]. Invalidated: the correlated energy
   ellipsoid separates all nominal forces from the field. Retain its finite
   old-domain BG<1 and force/field sine>2/5 certificates, actual factors and
   arithmetic bounds. Do not refine precision/subdivision of this relaxation.
2. **Independent nominal coefficients: exact nullspace.** f parallel B gives
   historical AG rank four; the complete-turn nominal h=.005,
   bhat_g=-400pi e_z example also has rank four and information-floor margin
   -mu. Invalidated: physical bounds or innovation bounds alone control a
   free nominal root. Construction reachability and all-time service for
   these relaxed examples were unproved. The implemented gyro projection now
   excludes the full-turn bias state; force/field collinearity and full temporal
   rank remain separate. Retain actual-history linkage.
3. **Endpoint-free forced adjoint: missing compatibility.** The exact failed
   equations are Z_i=Z_(i+1)A_i and W_i=Z_(i+1)K_i together. Zero-mean
   projection generally destroys them; the finite carried witness above
   also fails with unrestricted terminal multiplier. Invalidated: zero mean
   or spline jets alone leave only physical bias increments and defects.
   Retain both residual sums and signed physical Abel/velocity supplies.
4. **Pairwise unsigned tilt transfer: nonpositive margin.** The subtraction
   min(2A_max,(J_max+Omega_max A_max)T_E) already exceeds the maximum available
   gravity chord for every feasible theta_E<=pi/2. Invalidated: a pairwise
   triangle bound excludes nominal collinearity. Retain complete-window span
   and the exact signed pair relation; use whole-window physical integrals.
5. **Two margins to six pivots: unsupported implication.** No rank theorem
   proves that force/field separation and nonaliasing exhaust all varying
   six-column loss mechanisms. Coefficient compactness is also missing.
   Invalidated: simply taking p=min(Delta_col,Delta_gyr). Retain exact
   historical factor pivoting, root cancellation and the conditional p>0 bound.
6. **Universal V<=36 capture: refuted by admitted stillness.** For true
   R=Rx(2 atan(1/100)), p=v=a=omega=0, B=75e_x and
   b_a=g(R' e_z-e_z), ||b_a|| approximately .1961232<B_a, the measured record
   is level/north. The stationary service floor remains >1 after multiplying
   by (9999/10001)^2. At applicable regular boundaries P_ba,ba<=I/1600 gives
   V>=38468153689/625062500>61.54289>36. Invalidated: universal entry into
   the convenient projection-inactive ball or exact attitude/bias convergence
   in indefinite stillness. Retain the literal bias-projection sector;
   no larger/shaped invariant region has yet been proved.
7. **Old constant-attitude moving diagnostics: domain exclusion.** The smooth
   diagonal wave has finite 400--600-s mean tilt 8.211611 degrees, nominal
   acceleration 9.776391 and endpoint V>11492.6752 (guard margin <-11456.6752).
   These numbers are retained as diagnostics of the older domain, not current
   physical counterexamples. The 200-Hz stationary-looking moving witness is
   also excluded by jerk and tilt span. Its old-domain service proof remains.
8. **Local arithmetic factors: signed defect, not PSD noise.** The exported
   complete sync/symmetry defect has Rayleigh quotient -2^-44. Retain the
   full signed rank-one/two factor enclosure. Invalidated: treating that
   defect as a PSD process increment or restricting it to three AW coordinates.

9. **Future-only AG loss with a free prior: analytic obstruction.**
   P_root=diag(t I6,I15) gives D_AG,AG<=I6/t. At t=10^12 the margin against
   10^-6 I6 is <=-9.99999e-7. Retain historical root cancellation; future
   excitation alone cannot establish absolute J for that relaxed prior.
10. **Entrywise Riccati and cross ceilings: DEAD_END.** The former prediction
    enclosure [-18.7907040,2502.41442] loses acc/mag inverse verification;
    the cross-ceiling/Gershgorin margin is about -5.1223e8. Retain full matrix
    factors. LDL pivots are not eigenvalue floors and more entrywise
    subdivision is not the active approach.
11. **Restricted information lifting and endpoint process sums: invalid.**
    Heading/bias restriction I2 can coexist with full loss
    [[1,0,1],[0,1,0],[1,0,1]], cancelled by nuisance vector (1,0,-1).
    Interleaved corrections also defeat uncorrected accumulated-process
    lower bounds. Retain nuisance elimination and corrected endpoint paths.
12. **Automatic vanishing reset remainder: invalid.** At fixed nonzero
    injection d, the reset derivative J_l(d) differs from I+[d]/2. Smoothness
    alone does not prove eta(r)->0 for that comparison. Retain the explicit
    injection-dependent remainder and require a complete strict word margin.
13. **Forced-data rotation triangle: failed relaxation.** On the carried
    moving word, charge 39.22565272637 exceeds recorded projected gravity
    8.77133455729: margin -30.45431816909 before other supplies. Invalidated:
    separately norming the two rotations can close this budget. Retain the
    exact jointly signed rotation/reference action and actual magnetic gains.
    This does not falsify Delta_col or the physical assumptions. Current
    limiter: no uniform coupled rotation/reference bound. Next: include
    actual attitude dynamics and magnetic information before taking norms.
14. **Physical acceleration variation norm: stop after one refinement.**
    The cellwise velocity+jerk charge is 239.97053019538. Summing signed
    weights over actual S intervals first gives 23.98805965312+3.82238280532
    =27.81044245844, still margin -19.03910790115 against the same recorded
    threshold. Invalidated: either variation-norm bound closes the budget;
    no precision/subdivision retry is justified. Retain both exact integral
    identities, bias-rate Abel supply and the actual signed acceleration
    action (norm about .005 on this word). Architecture review: the single
    construction/H18/A21 path remains; the missing joint physical-integral,
    rotation/reference and gain bound cannot be replaced by separate norms.
    Finite small actions do not supply source-uniform ceilings.

15. **Exact-STILL entry/exit detection: analytic identifiability failure.**
    The new sin-cubed witness has exactly the same full IMU/magnetic record
    as true rest, including across C2 joins. Any causal detector entering on
    rest follows the identical state sequence through arbitrarily long hidden
    motion. Invalidated: finite dwell, gyro/acc variation or estimated v/aw
    alone guarantee both sound exact-rest entry and finite exit. Retain the
    stationary gyro information and physical ambiguity tube; no instability
    or failure of signed nominal margins follows. Limiter: no practical
    theorem for the whole compatible class. Next: derive that quotient/set
    storage before considering stationary mean or covariance changes.
16. **Finite bridges alone: switching-composition failure.** A moving norm
    gain 1/2 and finite bridge norm gain 3 give repeated gain (3/2)^n.
    Invalidated: individual finite retention establishes recurring stability.
    Retain exact every-prefix product/supply composition. Limiter: uniform
    product and transported-supply control through all switches. Next: derive
    that budget from actual operations or prove eventual regime retention;
    do not add a convenient physical switching/dwell assumption.
17. **Unrestricted reset inverse nonexpansion: exact algebraic obstruction.**
    With zero corrected rate, predict h0=3/625, apply d=4e_z,4e_z,(8/3)e_z,
    then eight h=1/200 predictions. Every literal reset is nonsingular with
    inverse norm <=1, yet B=diag(0,0,28/625). Invalidated: inverse nonexpansion
    and each one-prediction floor alone force full inter-anchor gyro rank.
    Retain exact inverse-frame recurrence and same-cell rows. This relaxed
    sequence has no shipping reachability or magnetic-service proof, and
    intermediate rows could restore full historical rank. Limiter: actual
    injection/geometry bounds. Next: bound those coupled to applied service,
    not by deleting resets or shrinking independent nominal boxes.
    Likewise the relaxed same-cell f=(1,0,1), b=e_x, d=2e_y gives Cf=0
    although raw f,b are nonparallel. The pulled-back field direction,
    not raw simultaneous separation, is the required geometry. Neither
    relaxed witness establishes source reachability or complete-row rank loss.
18. **Norm-summed NIS/covariance injection budget: overcharge.** Failed
    quantity: inverse-frame gyro floor over the 3-s inter-anchor interval of
    the collinear word. Summing `sqrt(NIS lambda_max(P_theta,theta))` gives
    4.98 rad against actual .0108 rad, so beta saturates and the budget is 0
    against actual 3.0. Classification: valid but quantitatively loose
    enclosure (per-correction ratio <=.134, typical far smaller). Invalidated:
    per-correction covariance/NIS norms, summed, control multi-second
    transport. Retained: the Loewner lemma and its use on short words (wave
    .2774 versus .28). One motivated refinement remains: `A~-I` is first
    order in the **signed** world injection sum (12 times below its norm sum).
19. **Covariance-normalized AW tracking: formulation failure.** Failed
    quantity: `sup_t(|e_aw|^2/lambda_max(P_aw)) sup_t lambda_max(P_aw)/1.12383^2`,
    6.68 on the carried 1-Hz collinear word and 6.91 at 25 Hz. A uniform
    retained radius must hold V>=106.6 (110.3) from the AW block alone, yet
    the storage route to Corollary A needs `r^2<1.12383^2/.0791`.
    Classification: high-precision feasibility ratio above one, structural:
    200-Hz acc corrections collapse lambda_max(P_aw) about 30 times between
    syncs while the jerk-driven lag error stays near .56 m/s^2, and every sync
    restores sigma^2. Invalidated: any AW covariance ceiling, measurement-aware
    or not, supplies `epsilon_a<1.12383` at a uniform retained radius.
    Retained: Lemma B and its tightness; the actual AW error .562 satisfies
    Corollary A on both words.
20. **Pointwise physical AW tracking: refuted on an admitted history.**
    Failed quantity: `sup|a_hat-a|` on 16-s A21 windows, 7.647 m/s^2 against
    1.12383 (ratio 6.80) for a C2-onset horizontal triangle 8.7 m/s^2 at 2.8 Hz
    (exact envelope: jerk <=99.9, |a|<=8.773, |v|<=.391). Classification:
    structural counterexample to a sufficient premise, not an estimator
    failure: the horizontal AW prior follows the vertical tuner at the .05
    sigma floor, and `Delta a_hat Delta a_hat'<=NIS P_aw` then needs NIS>=144
    per step to follow the jerk limit. Invalidated: any source-uniform bound
    on the pointwise AW error below 1.12383 (no refinement can cross it).
    Retained: Corollary A*, which needs only the signed nominal mean
    (carried worst .177 of its threshold). Do not strengthen MARINE MOTION.
21. **Perturbative injection charge over 16-s windows: infeasible.**
    Failed quantity: `delta^2(u_rms^2+1)/gamma*` with `delta` from Lemma I*,
    5.12 at the retained tilt, 1.79--7.17 at vanishing radius (only .96 for a
    quiet force RMS and an exact bias estimate), 700--1306 with the invariant.
    Classification: feasibility ratio above one; the deterministic .02 rad/s
    gyro residual is corrected by injections at that rate. Invalidated:
    charging `max|A~_k-A~_j|` against the attitude floor. Retained: Lemma I*
    (exact, physical), the third-order reset factor and the signed
    rotating-frame Corollary A** (feasible .22--.77 in the same table).
22. **Covariance-ceiling contraction from G0: existence only.** Failed
    quantity: `1-rho_0<=T b0/P_bg,max` with the least-squares reader ceiling
    from G0: 4.2e-9 per 100-s word (gyro block), 3.7e-13 (full action),
    3.4e-7 even at the actual floor 37; the scalar reader has
    `log10 B_*>=43367`. Classification: quantitatively vacuous but strict;
    b0=1e-11 makes the prediction comparison process-noise limited unless
    the ceiling is near the true P_bg. Invalidated: obtaining a useful rho_0
    from a least-singular-value floor through Q>=eps F C F'. Retained: the
    conditional implication chain. Next: blockwise reader action and an
    information-based (J_AG) contraction in the slow bias directions.

23. **LDLT pivot tolerance is not an eigenvalue tolerance: implementation
    failure.** For `N=6`, `tol=24*epsilon`, the symmetric matrix with
    `S_00=1` and lower `5x5` block `-(tol/2) ones` has an LDLT negative pivot
    `-tol/2` but eigenvalue `-5tol/2`. The accepted float result was
    -7.15256e-6 against tolerance 2.86102e-6; double was -1.33227e-14 against
    5.32907e-15. Invalidated: `min(D)>=-tol` guarantees
    `lambda_min(S)>=-tol`. Congruence preserves inertia, not eigenvalue
    magnitudes. Retained: the nonnegative-pivot fast path, the same numerical
    tolerance and unchanged harmless roundoff. Both the six-state projector
    and three/four-state OU regularizer now apply tolerance to eigenvalues.
    New scaled float/double regression checks pass. This does not supply a
    source-uniform arithmetic bound or close nonlinear retention.

24. **Ideal rotation identities are not a source-rounding certificate.** The
    shared SO(3) coefficients now switch on `x=|w|h` and use degree-18 Taylor
    series for |x|<1; the old `|w|<1e-7` branch description is invalid. Exact
    historical row factorization still uses actual Rs and Bs in W and Gamma.
    Simplifying `Rs^-1 Bs` to the exact rotation integral, or using a norm-one
    ideal rotation inverse, requires charging the finite Taylor remainder as
    well as mean-quaternion and floating-point defects. Classification:
    implementation-transfer obligation. Retained: the conditional ideal
    algebra, conservative root certificate and regenerated finite source
    audit; no source-uniform arithmetic transfer or theorem is promoted.
    The next falsifiable test is a directed remainder enclosure under the
    existing nominal-rate/step bounds, not changing the estimator for proof.

25. **First-prediction relative process comparison: structural ceiling.**
    - **Failed quantity.** `epsilon` in `Q>=epsilon F C_eta F'` with
      `C_eta>=P_root`. Testing it on one axis's S coordinate gives
      `epsilon<=Q_SS/(F L F')_SS<=3.3741e-10` for every admitted execution
      and every upper comparison (`B_*`, `U_n`, `eta`): the OU
      triple-integrator response is `<=t^3/6`, so one step has
      `Q_SS=O(h^7)`, and `L` is the certified S floor.
    - **What it can certify.** `N` predictions certify at most `N epsilon`,
      i.e. ≤1.73e-4 per 2048-s proof word.
    - **Carried confirmation.** The ideal ratio with the literal root
      covariance is 2e-23 (quiet) and 4e-22 (wave). The best
      `C_eta(B*,U_n)` chain is about 1e-36. The actual word margins are
      4e-4–7e-4 at 0.32 s and 0.39 at 64 s.
    - **Classification.** Structural and formulation-level; no reader,
      nuisance upper or Young split can lift it.
    - **Invalidated.** Obtaining a useful `rho_0` from any single-prediction
      (or per-prediction iterated) process comparison. This extends
      DEAD_END 22 from G0's scalar ceiling to every structured `C_eta`. The
      scalar root-precision cap (`J_root<=L_root^-1`, about 2e9) stays
      existence-only for the same reason.
    - **Retained.** The exact identity `D_pred=P^-1-(P+C)^-1`, the
      structured root upper `P<=diag((1+eta)B,(1+1/eta)U)` and the exact
      relative-Schur test of `D>=delta J_root` remain valid algebra.
26. **Separated information-only or forgetting-only word bounds: insufficient.**
    - **Failed quantities.** `1-lambda_max(P_0^-1/2 Sigma_00|y P_0^-1/2)`
      (information only) and `lambda_min(P_0^-1/2 Sigma_00|y,x_end P_0^-1/2)`
      (forgetting only).
    - **Carried results.** Quiet: 3e-26/2e-10 at 0.32 s up to 3e-12/5.2e-3
      at 64 s, against exact .389. Wave: 3e-10/2e-9 up to 3.0e-3/6.5e-3,
      against exact .393.
    - **Exact example.** `mixed_mechanism_example` has `rho=1/2` with both
      separated margins zero.
    - **Classification.** The feasibility ratio is 60–1e24 below exact.
      Slow directions contract through information and translation (or the
      quiet kernel) through forgetting, in different eigenvectors.
    - **Invalidated.** Any two scalar or eigenvalue-separated contraction
      bounds.
    - **Retained.** The exact smoother identity and the joint
      information-ratio lemma, which keeps both mechanisms in one matrix
      inequality.

27. **Root covariance matrix ceiling as the contraction input: redundant.**
    - **Failed quantity.** The joint-reader `C` loss (exact margin over
      certified margin) is 51.7× on the 16-s quiet word at 289 s, failing the
      10× kill criterion. On every MOVING word the optimum over `kappa` sits
      at `k=0`, where `C` does not enter.
    - **Young split.** The `U_n`-based `diag((1+eta)B*,(1+1/eta)U_n)`
      dominates the full joint `C` (λ_min 1.06–1.78). It cannot improve on
      `k=0`; its ~10^15 dynamic range also corrupts float64 values of `k`.
    - **Classification.** Formulation redundancy plus a transient. The quiet
      loss falls to 5.7× at 5064 s: the 289-s root's BA marginal (5.5e-5) is
      an order of magnitude below its 5064-s level, so its exact margin is
      transient.
    - **Invalidated.** That `I_eff>=mu`, `B_*` or a full 21×21 ceiling is on
      the critical path of `rho_0`. The lemma needs `P_0` only where `A-kappa J`
      is positive; on carried words that is the physical kernel `nu`
      (rank-one loss ≤1.9×).
    - **Retained.** The joint-reader identities (coercivity) and the lemma.
    - **Limiter.** The kernel-bounded diameter (O1, O2).

## Retained facts

World-frame row factorization, the reset Gram identity, attitude-invariant
same-cell geometry, the literal injection Loewner budget and the
isotropic-sync AW covariance ceiling are exact in real arithmetic, as are
Corollary A*, the AW-loop identities, Lemma I*, the third-order reset
factor, Corollary A**, Lemma T and the injection-free Theorem G0.
The sampled acceleration mean bound, constant-field joint physical vector
floor, full covariance-energy identity, LIN path action, full nuisance floor,
recurring nuisance upper covariance, historical reader factor/root algebra,
complete corrected-word matrix bootstrap, actual-gain finite-error composition,
reset remainder and projection-sector relations remain available.
Also exact in real arithmetic:

- the joint minimum-action reader, its diffuse-Riccati form, full 21×21
  domination and coercivity reduction;
- the smoother identity `M=C_end'P_0^-1`;
- the information-ratio word lemma;
- the relative-Schur characterization of `D>=delta J_root`;
- the first-prediction identity `D_pred=P^-1-(P+C)^-1`;
- the first-prediction ceiling;
- Theorem D (`kappa_W=lambda_max(J^-1 A)=lambda_max(Pi^-1 P_diff)`,
  `rho_W<=tanh(log(kappa_W)/4)`, sharp) and the closed-form supremum;
- composition (`kappa` shrinks under prefixing; `rho` multiplies), slow/fast
  factorization and the compression/reader duals;
- Corollary K and the invariance of `{P:nu'P nu<=c}` under `nu'P_nu nu<=c`;
- the S-chain cancellation of the `(v,p,S,a_w)` root, AW noise and syncs;
- the scalar kernel ceiling from `P_ba<=I/1600`;
- the source-uniform gyro-bias persistence cap on `kappa_W`.

The physical signed reserve is 27049050188592/625000000000000000, approximately
4.32784803017472e-5, under its stated physical field hypotheses. Its nominal
transfer ceiling remains unknown. The six-column-to-full contraction lemma
is conditional, not a certified numerical contraction for shipping histories.

The retained supplied-word 80-digit comparison has rho approximately
0.995037970315 and conditional decrement approximately 2.64851e-5; it is not
a source-uniform shipping result. Coarse existence quantities (nuisance ratio
approximately 1.02976e-18, scalar process floor approximately 2.49159e-26)
are not practical-radius certificates.

The theorem default profile uses g=9.80665, whereas the AtomS3R calibration-site
configuration uses default g=9.8025605 and its deployed S cadence. Other
configured profiles require their own source enclosure. No constants or
deployment behavior are changed to align them with the proof.

No H18 full-state contraction, independent nominal boxes, fitted/sampled rho,
scalar information lifting, scalar/Gershgorin process comparison, or finite
replay promotion is substituted for the retained proof path.

## Alternatives

For stationary operation, pursue a set/quotient practical storage retaining
physical tilt/BA ambiguity separately from estimator covariance. A gyro-only
observation must charge the entire compatible physical-rate tube. Do not hold
BA, force zero velocity, or shrink covariance solely from quiet-looking data.

For MOVING, work with the aggregate world-frame array: accelerometer rows
enter through nominal window statistics (Corollary A*), applied magnetic rows
at service gaps fix the components normal to B_w, and the gyro columns follow
the nominal rotation integral, whose field-axis coordinate Lemma T makes
monotone. No physical AW tracking, reference or gravity-mismatch transfer is
needed for the geometry; the body-frame rotation/reference action of
DEAD_END 13 does not arise. A same-cell route would instead have to derive
cadence from MAGNETIC SERVICE and couple it to the jerk lemma.

Downstream, certify the contraction of whole words (64 s or longer; the
persistence cap limits shorter words) through the kernel-bounded diameter,
not per prediction. The same formulation covers quiet and MOVING words, and
excitation only shrinks `kappa_nu`. G0, the nominal AW statistics, magnetic
service, the S-chain and the injection frame are the geometric inputs to
observability off the kernel. The joint reader remains the coercivity route.

## Analytical obstruction: the two proposed source-only scalar bounds

The literal recursions show that neither remaining scalar can be bounded by a
useful constant from MARINE MOTION / IMU BIAS / MAGNETIC SERVICE alone.

For the nominal AW mean, the exact loop identity is
`sum_acc Gamma(e-eta)=e_0-e_N+sum xi-sum_pred[(1-phi)a_hat+Delta a]`.
The physical increments `sum Delta a=a_N-a_0` telescope, but `eta`
contains attitude/BA/sensor residual and `xi` contains the AW increments of
other corrections.  The three physical assumptions do not bound these
estimator-error terms independently of the retained error storage.  Therefore
the desired source-only constant `m_perp_bar` does not follow from the
current premises.  This is distinct from DEAD_END 20: no pointwise AW
tracking is asserted.

For literal reset transport, Lemma I* gives the signed injection rotation in
terms of endpoint attitude errors plus the integrated gyro residual.  The
half-angle factor has the exact local expansion
`N=I-[x]/2+R_3`, `|R_3|<=|x|^3/6`, with product stretch controlled by
quadratic/cubic functions of the actual correction injections.  The physical
assumptions bound the gyro residual contribution but not the endpoint
attitude-error or correction-injection contribution independently of the
retained error storage.  Hence a source-only signed-reset constant is likewise
unavailable.

The correct formulation is radius-coupled.  On a candidate retained ball
`sqrt(V)<=r`, use covariance coercivity/projection guards to derive
`m_perp <= m_0+m_1 r+m_2 r^2` and
`delta_Q <= q_0+q_1 r+q_2 r^2+q_3 r^3`, with `m_0,q_0` the physical/source
terms.  Insert these directly into the local-tube quotient reader to obtain
`K(c,r)` and `d_bar(r)`.  The controlling closure is then the coupled
system
`d_bar(r) K(c,r)<=c` and the finite-error retained-radius inequality.
This ordering avoids falsely promoting an error-dependent geometry bound to
the linear source-uniform theorem.

## Coupled invariant assembly: current quantitative obstruction

The source-family action algebra now closes conditionally: S-chain removes
AW/root/sync sources, nuisance and measurement families have proved bounds,
and the AG family has the inverse-frame/leverage bound with literal whitening.
However the requested numerical two-dimensional solve is not yet well posed.

The reason is upstream of the action assembly.  The leverage estimate is

`M_C <= 2||F||^2/s(c,r)^2 [1+e_D^2 T h_max/q_T(c,r)]`.

Both `s(c,r)` and `q_T(c,r)` are local-tube G0 floors.  Their literal
values require the radius-dependent nominal signed-AW and relative-reset
bounds

`m_perp(c,r)<=A0+A1 sqrt(C(c,r))r`,
`delta_Q(c,r)<=Q0+Q1 sqrt(C(c,r))r+Q2 C(c,r)r^2+Q3 C(c,r)^(3/2)r^3`.

The exact identities defining these coefficients are proved, but finite
source-uniform numerical ceilings for `A0,A1,Q0..Q3` have not been derived.
The existing G0 numbers use supplied premises (for example m_perp=0.4 m/s^2)
and cannot be substituted into a theorem.  Therefore `q_T(c,r)`,
`s(c,r)`, `R_q^bar(c,r)`, `R_d^bar(c,r)`, and consequently
`K(c,r),D(c,r)` are still symbolic.

Likewise the finite-error composition is exact, but `E(c,r)` cannot yet be
evaluated because its actual-gain action coefficients for the complete word
and literal reset remainder are not all bounded numerically on the same
candidate rectangle.  The proved BA projection-inactive guard and AW
marginal are available but do not fill this gap.

Hence no rigorous positive or negative solution of

`D(c,r)<=c`,
`[1-sqrt(1-1/K(c,r))]r>E(c,r)`

can currently be asserted.  A grid search or substitution of carried-word
values would be fitted and is prohibited.

The next controlling analytical obligation is therefore not another
two-dimensional solve.  It is to derive finite radius-local coefficients for
the two signed-history maps:
(1) the AW-loop functional producing `A0,A1`, and
(2) the relative-reset partial-sum/product functional producing
`Q0..Q3`.
After those are inserted into G0, assemble the already derived factor-space
actions and only then solve the invariant inequalities.

## Attempt to derive A_i and Q_i: exact obstruction

The proposed summation-by-parts closure does not produce finite certified
`A0,A1` from the present assumptions.  The nominal mean representation is

`mu_hat=W0 a_hat_0+sum_c W_c Delta_c`, `0<=W_c<=1`.

Abel summation rewrites the correction term using partial sums of
`Delta_c`, but its coefficient is the total variation of the chronological
weights `W_c`.  Those weights depend on the literal OU prediction factors
and, through the correction increments, on the adaptive Kalman gains and AW
sync phase.  MARINE MOTION / IMU BIAS / MAGNETIC SERVICE do not bound that
variation.  The exact AW-loop identity replaces the raw increments by
gain-weighted `Gamma(e-eta)`, endpoint errors, other-correction increments
`xi` and prediction leakage; it does not remove the gain variation.  Thus a
radius-local `A1` requires an additional proved gain/weight-variation
inequality derived from the covariance recursion.  No such inequality is
currently in the proof.  Setting `A1` from carried sync-locked words would
be fitted.

The relative-reset coefficients have the analogous issue at second order.
Lemma I* bounds the **net rotation** of the ordered injection product by
endpoint attitude errors plus integrated gyro residual, so the first-order
signed term has a radius-local bound.  However the literal covariance factor
is `N=I-X/2+R3`, and its product remainder contains

`sum_l |x_l||S_(l-1)|/4 + exp(sum_l |x_l|^2/8)-1
 + sum_l |x_l|^3/6`.

Lemma I* bounds the final signed rotation, not `sum |x_l|^2` or the
partial-sum weighted quadratic term.  The Loewner injection lemma gives each
`x_l x_l'<=NIS_l P_theta,l`, but the current retained-storage argument has
no source-uniform bound on the sum of correction NIS/action over a complete
2048-s word before the contraction margin is known.  Consequently `Q1`
(first-order net rotation) can be expressed conditionally in terms of the
candidate radius and gyro residual, but finite `Q2,Q3` are not certified.

This is a genuine circularity in the present local-tube route:
the nonlinear storage/contraction would bound cumulative correction action,
while the literal G0 floor currently asks for that action to establish the
contraction.  The leverage bound removes accumulated gyro-history magnitude
but does not remove this reset-product remainder or AW gain-variation
dependence.

Therefore `q_T(c,r)`, `s(c,r)`, `K(c,r)`, `D(c,r)` and `E(c,r)`
cannot yet be assembled into a rigorous numerical two-dimensional solve.
The next proof must break the circle structurally, e.g. by formulating the
local-tube information directly with the exact reset factors `N_l` (using
`sigma_min(N_l)>=1` and their common signed action) so no Q2/Q3 product
remainder is needed, and by constructing the S-chain reader from literal
rows without first requiring a separate bound on the nominal AW mean.
No new physical assumption is implied by this diagnosis.

## Direct-information replacement of A_i/Q_i route

The A0/A1 and Q0..Q3 scalar route is retired.  The S-chain is now applied
directly to the raw auxiliary record, followed by whitening with its full
reduced covariance.  All literal nominal AW coefficients and all exact reset
factors N_l remain inside the reduced information matrix.  The new controlling
obligation is a same-history Loewner floor
`G_red(c;history)>=G_*(c,r)>0` on the quotient/kernel-augmented slow
coordinates.  Its quotient and gyro Schur eigenvalues are the `s(c,r)^2`
and `q_T(c,r)` needed by the already derived leverage/action bounds.

The exact identity `N_l'N_l>=I` guarantees inverse nonexpansion but does not
alone preserve force/field kernel angle, so it is not promoted to a G0 floor.
The proof must lower-bound the complete linked reduced matrix using MARINE
MOTION, MAGNETIC SERVICE and chronological gyro transport.  This avoids both
AW gain-variation and reset-product cumulative-action circularities.

## Quantitative variational-floor attempt: missing physical-to-nominal modulus

The optimal nuisance-annihilating Schur complement is exact, but the proposed
four-step contradiction is not yet justified by the current assumptions.

Magnetic service supplies a direct literal row modulus because the committed
world field is the coefficient of the magnetic attitude row and
`|B|>=B_min`; after whitening/nuisance projection this gives a transverse
distance proportional to the distance of attitude from the transported field
axis, modulo the exact shared-source metric.

The kernel row supplies the exact modulus `mu=1/c`.

Lemma T supplies a conditional chronological gyro modulus once the magnetic
residuals are small: its constant `q_I` is explicit in terms of service gap,
field floor, rate bound and the **accelerometer-window transverse geometry**.

The unresolved step is that accelerometer geometry.  MARINE MOTION constrains
the physical attitude/gravity direction, but the literal accelerometer
attitude row in `O_s` is built from the estimator nominal force
`a_hat-g e_z`.  Nuisance elimination allows `a_w,b_a,v,p,S` root
directions to mimic parts of this row.  The Schur complement removes those
directions optimally, but no existing theorem lower-bounds the remaining
distance of the literal nominal accelerometer row from the nuisance span
using physical attitude span alone.  Corollary A* did this through a nominal
signed-AW premise; that premise was exactly what the direct formulation was
intended to avoid.

Therefore compactness/contradiction cannot presently conclude
`G_red,mu>0` for every admissible MOVING history: a hypothetical sequence
may keep the physical attitude excitation while its estimator nominal force
approaches the magnetic-axis/nuisance-compatible geometry.  No current
assumption or proved estimator invariant excludes that sequence.

This is not repaired by `N_l'N_l>=I`: reset noncontraction preserves
invertibility, not the missing physical-to-nominal force separation.

Hence the direct variational route has reduced the gap to one precise
modulus:

`dist_(Sigma^-1)( O_acc,slow x_s,
                   range[O_f, magnetic-compatible nuisance] )
 >= a_phys(r)|x_tilt/BA|,quad a_phys(r)>0`                  (VF-A)

on every required moving excitation window, with all literal linked
coefficients retained.  A proof of (VF-A) may use the exact innovation
identity `f_measured-b_hat-r_acc = f_hat`, the physical sensor/bias bounds,
and the candidate storage radius to relate nominal force to physical force,
but it must not assume pointwise AW tracking.  Until (VF-A) is proved,
explicit magnetic/Lemma-T/kernel constants cannot be combined into a
positive `g_*(c,r)`.

## Physical-to-nominal projected modulus: innovation route fails

The exact accelerometer innovation identity does not supply (VF-A).

Shipping has
`r_acc=f_meas-[R_wb(a_w-g)+lever+b_a]` and
`J_att=-[R_wb(a_w-g)]x`, `J_aw=R_wb`, with the BA row identity when BA
updates are enabled.  Thus the physical measurement can be written as the
nominal prediction plus the realized innovation, but the attitude Jacobian
coefficient remains the estimator nuisance state `a_w`.

In the variational information
`min_xf ||Sigma^-1/2(O_s x_s+O_f x_f)||^2`, the nuisance minimizer is free
to vary exactly the root AW/BA coordinates whose historical transport enters
the accelerometer rows.  The realized innovation value is not a column of the
frozen root design; it is data.  Bounding it by sensor/bias/storage envelopes
therefore bounds a finite-error defect, not the distance of the slow Jacobian
column from `range(O_f)`.

Consequently the substitution
`f_hat=f_meas-b_hat-r_acc` cannot prove a positive nuisance-projected
attitude/BA information modulus without an additional relation restricting
the nominal AW/BA history.  Doing so would reintroduce, in another form, the
refuted pointwise AW-tracking premise or an unproved gain/history bound.

This identifies a structural obstruction to the proposed MOVING
observability theorem under the current assumptions: MARINE MOTION constrains
physical motion, but the linearized root-information matrix is built from
estimator nominal coefficients.  The current assumptions do not guarantee a
uniform separation between those coefficients and the magnetic-compatible
nuisance span.

Therefore (VF-A), DI-6 and a source-uniform positive `g_*(c,r)` are not
proved by the innovation route.  No numerical K,D,E or invariant solve may be
claimed from it.

A successful continuation needs a genuinely different mechanism already
present in shipping, for example an information argument using the actual
closed-loop correction/gain history rather than raw frozen Jacobian
observability, or a reachability theorem showing that the problematic nominal
coefficient histories cannot occur inside the candidate retained set.  Adding
a new physical assumption or changing estimator behavior is outside the
current task.

## Complete-word replacement after two-epoch failure

The controlling path now uses the complete corrected A21 word as one joint
Gaussian operator.  All root nuisance coordinates and every fresh
process/sync/noise factor are retained once; nuisance elimination is the full
Schur complement, and the kernel prior is appended afterward.  The target is

\`P_(nu,s)(W,c)<=K_MW(c,r) Pi_s(W)\`

uniformly over admissible MOVING words in the candidate retained region.
This lets magnetic service, all accelerometer epochs, S pseudo-observations,
chronological gyro transport, AW/BA process penalties, Joseph corrections,
resets and terminal forgetting cooperate.

The immediate analytical subproblem is the zero-action nullspace: prove that
a normalized complete-word sequence for which magnetic loss, accelerometer
loss after nuisance mimic, S loss, AG/BA/AW process action and terminal
forgetting all tend to zero converges only to the physical tilt/BA kernel.
Only after that qualitative coercivity statement is proved should explicit
moduli be extracted for K_MW.

## Shipping-invariant audit for the restricted detectability gap

All current shipping guards/service conditions were checked for a mechanism
that could lower-bound the local restricted gap needed by terminal retention.

- MAGNETIC SERVICE is genuinely coercive: on every certified T_M interval the
  sum of actually applied transported/whitened magnetic rows has a prescribed
  2-D information floor mu_M. It controls normalized heading/axial-gyro-bias
  root coordinates and is already the correct source-uniform service premise.
- The gyro-bias projection (0.5 rad/s) excludes complete-turn aliases and
  bounds transport rate, but does not impose an angle between the remaining
  tilt/BA compatibility line and accelerometer rows.
- The accelerometer-bias projection (0.4 m/s^2) and BA OU decay bound the BA
  component and keep the kernel family compact; they do not create
  transversality.
- S cadence plus the tau_aw clamp [0.02,12] gives the explicit four-S
  LIN/AW injectivity modulus, but S rows have no direct attitude/BA row.
- Racc/R_S/tuner clamps bound whitening and action constants. NIS/LDLT gates
  decide whether a correction is applied; they do not require a minimum
  attitude/BA information angle for an accepted accelerometer correction.
- Literal reset factors satisfy sigma_min(N)>=1 and inverse nonexpansion;
  this preserves invertibility but not force/field transversality.
- AW sync/floor operations bound covariance and nuisance action; they do not
  constrain nominal specific-force direction.

Therefore no unused shipping invariant supplies the missing restricted
tilt/BA gap. MAGNETIC SERVICE closes its intended 2-D heading/gyro sector,
and four-S closes LIN/AW, but the one-dimensional tilt/BA compatibility
sector can approach tangency continuously under the current MARINE MOTION
assumption.

Conclusion: with the present physical assumptions, a source-uniform
finite-horizon detectability constant C_det (and hence a source-uniform
strict recurring contraction margin obtained by this route) is not proved.
The obstruction is not an omitted implementation guard. To obtain such a
margin one needs either a different theorem that tolerates arbitrarily weak
tilt/BA detectability without a uniform ratio, or an additional quantitative
physical/service condition that excites that sector. Adding such a condition
would strengthen the assumptions and is outside the current task unless
explicitly authorized.
## Corrections to PR #625 analytical claims

Three earlier claims are corrected fail-closed:

1. The explicit four-S determinant/singular-value floor is RETRACTED. The
   inequality using tau<=12 had the exponential monotonicity reversed.
   Qualitative four-S injectivity remains valid; no quantitative four-S
   modulus is currently certified.
2. LE-4 must use the residualized gyro operator
   `(I-P_U)X_g`. The unprojected inequality was false. Consequently
   LE-5--LE-8 and the AG-process energy ceiling derived from them are
   RETRACTED pending a residualized rederivation.
3. Trivial pointwise augmented nullspace plus compactness does not imply a
   uniform positive g_MW across rank-changing/tangent kernel histories.
   The uniform-coercivity existence claim is RETRACTED. The later
   finite-horizon detectability/relative-action formulation is controlling.

Implication: K_MW(c,r), D(c,r), a strict recurring contraction margin and the
final invariant solve remain OPEN. No downstream certificate may cite the
retracted claims.
## Exact-compatible MOVING + strict MAGNETIC SERVICE construction test

To obtain an exact complete-word tilt/BA kernel, the simplest case sets the
kernel BA component to zero. In injection-free world coordinates the
magnetic-compatible attitude direction is the world field b. Zero
accelerometer loss then requires every applied nominal specific-force row to
satisfy `[f_hat_k]x b=0`, i.e. f_hat_k parallel b. With nonzero BA the
condition generalizes to a fixed/decaying transverse component
`[f_hat_k]x b = -R_ba,k phi_b,k b_a0`.

For the physical force, the existing jerk/velocity lemma rules out exact
field collinearity at every dense applied correction under the documented
marine bounds (at h=1/5 and L=16 s it gives a positive minimum weighted
collinear-sample gap, about .051 s, while regular accelerometer corrections
are much denser). Thus an exact-compatible construction cannot simply make
the physical force satisfy the kernel equation at every accelerometer row.

However the literal kernel equation uses estimator nominal force f_hat, not
physical force. The admitted assumptions constrain physical motion/bias and
MAGNETIC SERVICE, but there is currently no proved reachability invariant
forcing f_hat to inherit the physical jerk/velocity anti-collinearity. The
pointwise AW-tracking premise is explicitly refuted. Therefore the physical
jerk lemma cannot rule out an exact-compatible nominal history.

Conversely, constructing such a nominal history is not free: it must arise
from one actual shipping execution with the accelerometer correction, AW OU
prediction/sync, S corrections, tuner state and accepted measurements. No
existing theorem proves that an exactly collinear nominal-force sequence is
reachable while the physical force is not collinear. Carried sync-locked
examples are only approximate and cannot be promoted to an exact witness.

Strict MAGNETIC SERVICE itself is not the obstruction. It constrains the
transported magnetic heading/axial-bias rows and can remain strict under
small perturbations of translational acceleration/AW history. It does not
directly constrain the tilt/BA nominal-force compatibility equation.

Result: under the current proof state, existence of an exact-compatible
MOVING shipping execution with strict MAGNETIC SERVICE is neither constructed
nor ruled out. The question has reduced to a shipping reachability problem:
can the closed-loop AW/BA/S recursion realize the exact affine nominal-force
constraint at every applied accelerometer epoch while physical MARINE MOTION
remains admitted? A proof must use the literal mean recursion; physical
geometry alone cannot decide it.

Therefore the conditional kernel-disappearance counterexample cannot yet be
promoted to a disproof of source-uniform C_det, and finite C_det cannot be
proved by excluding the base word either. This reachability question is now
the controlling blocker.
## Next analytical step

Kernel-bounded observability certificate (O1, O2). Derive explicit symbolic
bounds before any new source run or enclosure.

1. **O1 reader.** Build a terminal reader that cancels every root direction
   except `nu`:
   - the `(v,p,S,a_w)` root through the S-chain identity;
   - the tilt/BA combination from S-chain accelerometer windows;
   - heading and tilt normal to the field from applied magnetic rows at
     service gaps;
   - the field-axis gyro bias from Lemma T with two separated windows (G0).

   Preserve the same-history signed coefficients. The kernel coordinate is
   supplied only by the fictitious precision `mu=1/c`.
2. **Known-root scalar action.** Bound directly
   `d_j=nu_(j+1)' Pi_j nu_(j+1)`; do not first seek a full 21-state ceiling.
   The S-chain should cancel the neutral/AW root and syncs before bounding
   process action.
3. **Close O2 as a fixed point.** From
   `P_end<=P_nu<=kappa_nu Pi`, prove linked bounds
   `kappa_nu(1/c)<=K(c)`, `d_j<=d_bar` and exhibit
   `d_bar K(c_bar)<=c_bar`. Prefer the sharper scalar Schur recursion
   `c_next<=d+a c/(1+b c)` if it exposes BA decay/kernel information without
   an independent attitude/BA split.
4. **Only then falsify constants.** If the symbolic reader produces explicit
   constants, evaluate the resulting `K(c)`, `d_bar` and fixed-point margin
   on the existing carried 64-s words as a non-promoting check. Reject the
   construction if its diameter is more than 10× the actual carried diameter
   or if no positive fixed-point margin exists.

The G0 extensions below remain the geometric inputs to these floors.

Extend Theorem G0 to the literal array: carry the injection frame as a
rotation Q (Lemma I*, half-angle remainder charged by signed partial sums),
apply Lemma T on local tubes of about pi/Omega with the kernel direction
`Q'b` frozen per tube, and chain the monotone field-axis coordinate across
tubes. A kernel turning at kappa=|omega_tilde|/2 drags the transverse residual
at kappa|eta|, so the chained speed is about `c0-kappa|eta|_max`: at .01 rad/s
this fails on 100-s words at the invariant rate (kappa T_w=1>c0) but leaves
room on ~40-s words with a physical rate (c0 about .6, gap G about 18 s),
provided the attitude-error jitter satisfies `alpha T_w<<eps`. Evaluate the resulting floor at 80 digits on the synthetic
falsification families with injection sequences at the deterministic
.02 rad/s residual and on the carried world-frame words; reject it if it is
nonpositive or exceeds an actual value. Separately, test whether the
nominal window statistics admit a loop bound: evaluate the gain-weighted AW
identity with the actual sync-cycle gains against the worst phase-locked
jerk-limited input (currently .371 signed, .348 transverse) and report the
worst admitted ratio against 1.96133. Passing finite tests remain
non-promoting. Do not repeat pairwise, unsigned-energy, variation-norm,
norm-summed or perturbative injection tactics.

## Validation and infrastructure

Based on main `ee0d46d`. This continuation adds proof tooling, evidence and
documentation only:

- `word_diameter.py` (exact) with `word-diameter-certificate.json`;
- `information_ratio_source_diagnostic.py` with its native driver
  `tools/stability/information_ratio_source.cpp` and
  `information-ratio-source-feasibility.json`.

No estimator source, tuning, gate or assumption is changed for the proof.

- **Native records.** They are generated from the shipping header through
  read-only taps with observer/control terminal parity. The
  information-ratio record keeps four significant digits of well-posed
  quantities only; unobservable quiet diameters are recorded as infinite,
  not as float64 values. CI reproduces it with `--expect` at relative
  tolerance `5e-3`, because float32 native replays differ at about `1e-6`
  between toolchains (below).
- **Checks.** `build_evidence.py` must reproduce every exact certificate and
  verify every diagnostic record. The focused `test_ou3_*.py` suite, ruff and
  `git diff --check` must pass, and `make all` is the primary validation.
- **Full evidence.** Full validation/robustness and TFG evidence are
  regenerated by the existing full-evidence pipeline, never certified by
  editing only their fingerprints.
- **Still to do.** The carried readout record
  (`ag-readout-source-feasibility.json`) predates the core fix and lacks the
  `contraction_feasibility` fields. A canonical record must come from the
  Ubuntu workflow artifact, which this environment cannot download.

### CI native replay binding (2026-09-29)

Run `36568721853`, job `109407034842`, failed the world-frame `--expect`
comparison on main `069c6b5301a79b2499476e6d07592d058aa8432c`. The quiet and
wave cases matched; the collinear cases differed in their trace hashes and
78 derived metrics. For example, the 1-Hz NIS maximum was
`17.4755170934018834186006` in CI versus `17.4755343260261939822796` in the
committed record. All source hashes, observer/control terminal parity and
`verify_diagnostic` invariants passed.

Classification: cross-environment native replay mismatch, not a failed
mathematical bound. The old record reproduces exactly with local GCC 14.2,
Eigen 3.4 and glibc 2.41. Independent Ubuntu 24.04 executions on main and
PR #622 (`36568698113`) produced byte-identical world-frame records. The
individual compiler/library contribution is not isolated. Invalidated:
identical source hashes imply bit-identical native traces across these
environments.

The canonical world-frame fixture is now the actual `ou3-world-frame.json`
from main's evidence artifact `11034245782`, independently matched against
artifact `11033851278`; its Git blob is
`2e6b651a638a7033466d2f144f491e902e4cd70c`. Only that fixture's provenance
binding is refreshed. The strict trace-hash comparison, decimal tolerance,
all invariant/quality gates, shipping source and false/open theorem flags
are unchanged. Use the workflow's Ubuntu environment for canonical replay;
this repair does not claim cross-toolchain bitwise portability.

Validation: `build_evidence.py` passes; 649 validation tests pass and one
simulation-record test is skipped because its inputs are unavailable.
The initial local compass setup failure was resolved by setting
`EIGEN_INCLUDE_DIR`, without a repository change. The noncanonical local AW
replay passes its source/invariant checks but differs in nine exact numeric
fields; its committed fixture is deliberately not replaced by local output.
Current limiter: native environment portability and the unchanged open
mathematical obligations above. Next falsifiable check: rerun the unchanged
world-frame and downstream AW `--expect` steps on Ubuntu CI. Do not weaken
comparison tolerances or promote a finite replay to fix an infrastructure
failure.
