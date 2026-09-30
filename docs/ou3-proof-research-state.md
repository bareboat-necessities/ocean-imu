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

**Causal terminal-AW correction (post-PR #630 audit).**  The claimed exact
AG-root cancellation was false: for one correction a terminal AW row q with
q E_h=0 pulls back to q(I-KH)E_h=-qKHE_h in general.  The useful sharp result
survives without that claim.  Starting at an actual post-sync boundary with
the full carried root covariance, Joseph telescoping makes the complete causal
reader action exactly u'P_Nu, and the literal AW sync/prediction/correction
chronology gives u'P_Nu<=16.  Thus sqrt(B_AW,*)<=4 remains source-uniform with
the AG root retained in the single action budget.  The physical primitive map
has zero root-uncertainty columns.  The unrestricted full action-source norm is also too coarse: process and
accelerometer source factors are orthogonal in covariance action, so taking
that norm before pairing destroys their exact same-history opposite-sign
cancellation.  The surviving constant is therefore the induced bilinear norm
restricted to the literal reachable causal-reader rows, not an arbitrary
reader cone.  The committed 17-s diagnostic previously
called C_port,W is not G_phys,W: it is only ||D2 h||_(2,1)/sqrt(action), so it
cannot be compared with the 8.45e-3 target.  The controlling obligation is to
construct the complete two-Abel physical-primitive -> normalized source
operator before taking norms, then enclose its induced norm source-uniformly.

**Exact sync-slab/17-s physical telescope.**  The physical LIN lift
chi=(v,p,S,a) with delta=chi_next-F_L chi has AW component
a_next-phi a exactly.  In physical-error coordinates S=0 contributes the real
forcing -K_S S_true with the single fixed capture origin, while deterministic
accelerometer elimination contributes no independent physical acceleration.
For one actual sync slab the prediction plus S forcing telescopes exactly to
lambda_L' E_L chi_L-lambda_R' E_L chi_R plus the missing accelerometer reader
jump sum q_a' R_wb a_phys.  Hence the S_true/P_AC term cancels the neutral
boundary created by the second Abel step and must not be charged eventwise.
The joint accelerometer+S completion is one chronological innovation square;
its positive q'Omega q term equals the covariance-storage decrement
lambda'(P^- - P^+)lambda and is retained with Joseph storage, not assigned a
second factor-four budget.  Summing all sync slabs over 17 s cancels every
internal physical boundary before norms.  OPEN: bound the resulting single
signed endpoint+accelerometer functional source-uniformly tightly enough for
the 0.0338 m/s^2 residual margin.

**Signed field-axis rotation x velocity telescope.**  On the retained local
field-axis chart, exact Stieltjes integration by parts gives one correction
cell as external M_b v endpoints minus the continuous gyro/gyro-bias rotation
pairing and the exact SO(3) correction jumps.  For each correction the exact
field-axis increment is retained; kappa-b' E_theta K r is charged as the
nonlinear twist/reset remainder.  The linear signed coefficient is
c_i=u'M_b[b]x v_i, so q_i=c_i E_theta' b remains inside the same-history
Joseph port.  Accelerometer AW and BA pieces are parts of the same
accelerometer innovation, S can rotate attitude through carried cross
covariance, and magnetic corrections remain sequential.  With
z_i=K_i' q_i, z_i'Omega_i z_i=q_i'(P_i^- - P_i^+)q_i exactly.  Summing every
17-s correction cell cancels all internal physical endpoints before norms and
produces one augmented full-root factor row.  OPEN: source-uniform enclosure
of that row including continuous gyro/gyro-bias transport and nonlinear
twist/reset remainder.  Raw TV, total NIS, V_max times unweighted covariance
loss, and the old 0.0338 residual margin are not promoted.

**Continuous field-axis enclosure and obstruction.**  The slow residual
gyro-bias term admits a second physical Abel bound
2 P_max B_g/T + P_max D_g = 0.01913992353 m/s^2 at T=17 s (apart from the
joint rotation/chart coefficient derivative).  The commissioned fast gyro
residual has only an amplitude contract, giving V_max N_g=0.11 m/s^2.
Moreover v=5.5 sin(t), p=-5.5 cos(t), a=5.5 cos(t), jerk=-5.5 sin(t), and
n_g=0.02 sin(t)b satisfy the declared primitive/fast-residual envelopes and
keep field-axis error amplitude at 0.02 rad while producing signed mean
0.055+O(0.0004) m/s^2.  Hence the old 0.0338084 margin cannot close the
continuous fast channel from independent envelopes.  This is not yet a
shipping counterexample because the literal correction chronology has not
been solved on the witness.  Exact SO(3) jump remainder is <=|v| kappa^2/2
plus the twist-coordinate remainder; no useful all-word sum of these terms is
currently proved.  NEXT: either prove same-history cancellation of fast gyro
against literal acc/S/mag jumps in the augmented full-root functional, or
extend the sinusoidal witness through those literal corrections and all
service/gate contracts to obtain a genuine shipping-reachable counterexample.

**Fast-gyro witness through literal measurements.**  Strengthened the
0.055-envelope witness by taking true attitude R=Exp(theta[b]x),
theta_dot=-0.02 sin(t), and fast gyro residual +0.02 sin(t)b.  The supplied
gyro is then exactly zero and R'B=B, so the level nominal gyro/magnetic
samples are exactly compatible; zero magnetic residual does not remove
magnetic covariance/service.  For any candidate nominal AW w, physical
a=g e_z+R(w-g e_z) makes the accelerometer sample exactly the level nominal
prediction.  A 5.5 cos(t) AW component needs only <=0.19613 m/s^2 additional
gravity compensation and remains within acceleration/jerk envelopes after
O(theta^2) DC centering.  The obstruction is therefore not the measurement
equations but self-consistency of the literal OU+S+acc covariance/gain orbit.
Over one period it is the finite-dimensional system
(I-A_per)x=B_per u plus sampled AW compatibility, exact innovation identity,
BA recurrence and the actual covariance-generated gains.  Net periodic
attitude balance does not imply the needed weighted cancellation because the
dangerous correction sum carries c_i=u'M_b[b]x v_i.  No structural identity
found forces that weighted sum to cancel.  NEXT: solve/exclude the coupled
periodic Riccati/mean compatibility system under the retained tuner chronology.

**Periodic Riccati/mean system analytically reduced.**  For a prescribed
periodic physical input, the stable front-end/tuner filters converge to a
periodic tau/sigma/R_S/T_S word.  Covariance evolution for that word is
independent of innovation values; its stabilizing periodic orbit fixes the
gains.  The mean period equation is (I-Phi)x0=G u+d.  If I-Phi is invertible
(as expected for a corrected contractive word), every periodic innovation
word has a unique periodic mean root.  For fixed gains, measurement-word to
innovation-word transport is block lower triangular with identity diagonal,
hence bijective.  Eliminating root and innovations leaves one physical
equation [I-R T_aw,a]a = g ez-R g ez+R t_aw+a_free.
Nonsingularity creates a unique periodic candidate rather than excluding it.
Therefore the periodic Riccati/mean architecture supplies no universal
same-history cancellation lemma.  Structural analytical exclusion has failed.
NEXT: interval-certify one actual coupled periodic tuner/covariance orbit and
its physical solution, then check MARINE/IMU/MAGNETIC/gates; or prove every
solution of that finite equation violates one existing condition.  Interval
arithmetic here would certify a constructed trajectory, not promote a fitted
diagnostic constant.

**Complete 17-s homogeneous linear word.**  Constructed the literal
21-state chronological mean map using prediction F, accepted acc/S/mag
I-KH, literal attitude resets, identity-mean PSD AW sync, active BA prediction,
and the carried covariance/source factors.  Exact expansion gives
P_W=M_W P_0 M_W'+S_W and D_W=P_0^-1-M_W'P_W^-1 M_W.  Therefore
rho_lin(W,P0)=lambda_max(P0^(1/2)M_W'P_W^-1 M_W P0^(1/2))
=1-lambda_min(P0^(1/2)D_WP0^(1/2)); innovations are eliminated exactly.
The existing smoother identity gives the equivalent information form and
Riccati diameter.  With the physical tilt/BA rank-one kernel, the surviving
source target is sup_W kappa_nu(W)<=K_17<infinity plus same-word scalar kernel
invariance, yielding rho_lin<=1-1/K_17.  Existing carried 16-s feasibility
words have exact rho about .884--.993 and difficult kappa_nu about 196--198,
so the plausible theorem scale is rho~.995, but these numbers remain
non-promoting.  OPEN: source-uniform K_17 over all retained histories.  Do not
claim a numerical rho until aggregate world-frame geometry, MAGNETIC SERVICE,
jerk/sample fidelity, S-chain and gyro persistence are enclosed jointly.

**K17 source-uniform target audited and NOT closed.**  The requested
combination of existing constants on one 17-s word is invalid: G0's explicit
s^2>=1.486786e-3 uses two 16-s windows separated by 64 s plus endpoint room
(~100 s total); the jerk fixed-attitude alias exclusion is 32 s; joint
measured-vector information is 64 s.  The 17-s nuisance comparison, S-chain
cancellation and P_ba<=I/1600 are compatible with 17 s but do not alone give
the missing slow quotient floor.  G0 additionally assumes nominal
m_perp<=0.4, u1<=1.2 and injection-free transport; these remain source-open.
Therefore carried kappa_nu~196--198 cannot be promoted and K17<=250 is not
proved.  Qualitative radius-local finiteness follows from the nuisance Schur
variational form G_red,mu=G_red+(1/c)nu nu' on a compact retained class if
its nullspace contradiction is made uniform.  Productive 17-s target:
derive an explicit same-word floor G_red,mu>=g17(c,r)I, then convert with the
exact S-chain/fast factor to K17(c,r).  Alternative: move contraction to a
>=64/100-s superword where existing geometry actually applies.

**Proof clock changed to 100-s contraction superword; estimator unchanged.**
17 s is retained only as the nuisance/root warm-up and intermediate-root
control horizon.  A complete moving 100-s word legitimately contains G0's
two 16-s windows with 64-s separation, the 32-s alias exclusion, 64-s joint
vector information and every 1-s magnetic-service interval on one history.
Define K100=sup lambda_max(Pi100^-1 Pnu,100); kernel invariance gives
rho100<=1-1/K100.  This removes the previous horizon mixing.  It does NOT
close K100 numerically: G0's nominal m_perp/u1 premises and literal
injection-frame extension remain source/radius-open, and adjacent superwords
still require exclusion of the unit-persistent exact-kernel equality in the
scalar return.  Once those close, use the existing exact every-prefix
composition to retain all intermediate 17-s roots and finite-error supplies.
No estimator schedule or assumption changed.

**Literal 100-s G0 radius-local attempt reduced to one modulus.**  Tried
the natural m_perp(r)=m0+A1 r+A2 r^2 and injection
delta_Q(r)=q0+Q1 r+Q2 r^2+Q3 r^3 route on the same 100-s word.  It remains
circular: the signed AW mean needs a source-uniform bound on variation of the
literal OU/Kalman/sync weights (time-varying gains can rectify); the reset
product needs complete-word bounds on sum|x_i|^2 and partial-sum quadratic
injection action, which are not available before contraction.  Therefore
these scalars are not promoted.  Adopt the exact nuisance-reduced information
G_red,mu=O_s' Sigma^-1/2(I-P_f)Sigma^-1/2 O_s+(1/c)nu nu'.
All literal AW coefficients, resets, S-chain cancellation and source
correlations remain inside it.  The 100-s target is G_red,mu>=g100(c,r)I.
The sole missing quantitative implication is physical-to-nominal
accelerometer separation: vanishing reduced nominal accelerometer information
must either force the true MARINE attitude/gravity span to vanish or incur
positive correction/process action already counted in Sigma.  This is now
the next controlling lemma.

**Physical-to-nominal separation lemma refuted in standalone form.**  Exact
two-epoch accelerometer/BA elimination gives
L2=C1 T_theta-phi_b B1 B0^-1 C0.  A nonzero zero-loss mode is exactly in
ker L2.  MARINE attitude span does not imply sigma_min^+(L2)>0 because the
literal C_i=-[fhat_i]x use nominal specific force; admitted translational
acceleration can compensate gravity-direction change while physical attitude
spans.  Existing collinear constructions realize this mechanism.  Therefore
do not seek a standalone accelerometer floor from attitude span.  On the
complete 100-s word, zero fresh action rigidifies nuisance trajectories,
four-S compatibility kills free homogeneous LIN/AW, magnetic service reduces
AG to at most one transported field-compatible line, and all accelerometer
rows leave at most one word-dependent compatibility kernel nu_W.  The
controlling quantitative target is now relative terminal/action coercivity
H_eff' Pi^-1 H_eff <= K_rel S_q on the quotient of nu_W, plus adjacent-word
kernel return.  A vanishing information eigenvalue is acceptable when the
terminal excess vanishes at the same rate.

**100-s relative quotient inequality derived exactly.**  After nuisance and
word-kernel shorting, write y=O_q v+A s, x_N=T_q v+T s, Sigma=AA'.  Then
J_q=O_q'Sigma^-1 O_q and Ttilde_q=T_q-T A'Sigma^-1 O_q.  Schur elimination
of the kernel line gives S_q and H_eff.  The sharp fixed-word constant is
K_rel-1=lambda_max(S_q^dagger/2 H_eff'Pi^-1 H_eff
S_q^dagger/2), with infinity exactly if Null(S_q) is not contained in
Null(H_eff); equivalently Null(O_q) subset Null(T_q).  Complete-word
zero-action classification supplies this for each fixed nondegenerate word,
but not uniformly through rank-changing sequences.  Remaining quantitative
lemma: a near-null quotient direction must satisfy terminal action <= C_det
times observation action uniformly.  A regularized literal backward reader
is an equivalent constructive route and gives K_rel<=1+C_det.  No Euclidean
information floor is needed.

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

**Nominal-AW complete reader (PR #630).** The raw beta-TV route has been
replaced by an exact complete covariance-weighted scalar reader.  For a
transverse multi-time AW readout q, augment the chronological state with its
readout accumulator and write the literal factor recursion once.  If M is the
root transition and Z the single signed coefficient on the stacked whitened
source vector, then

`q mu = q M e0 + Z s + d_impl`

and

`|q mu| <= C_root sqrt(V0) + C_src sqrt(A_W) + C_phys + C_impl`.

The stronger joint form is

`sqrt(C_root^2+C_src^2) sqrt(V0+A_W)+C_phys+C_impl < g sigma_w`.

Accelerometer and magnetic residual actions are square-summed in `A_W`
before taking norms; the nuisance-correlation term
`C H_n' K_L' lambda` is not a new source.  In the full covariance-weighted
adjoint every optimal correction obeys `P- p-=P+ p+`; the reduced
correlation term is exactly the bookkeeping residue from deleting nuisance
coordinates.  S corrections have zero weighted-adjoint jump.  Direct magnetic
mean forcing is its ordinary `r_mag' R_mag^-1 r_mag` action.  Prediction and
PSD AW sync remain chronological covariance factors; tuner lag qualifies the
coefficient trace and is not an independent mean forcing.

The scalar minimum source action is the dual of the existing joint
minimum-action reader
`B*=Pi+Ttilde I_eff^-1 Ttilde'`; no parallel observability architecture is
introduced.  However this does NOT numerically close the nominal-AW premise:
the currently useful source-uniform bounds for that reader still invoke G0
geometry whose premises include the nominal signed-AW statistic itself.
Substituting carried .348--.371 m/s^2 values would therefore be circular/fitted.
The new controlling obligation is a noncircular source-uniform evaluation of
the scalar reader coefficients/actions on the retained coupled
Riccati/tuner/physical trace class.

Keep the threshold symbolic as `g sigma_w`.  The value `g/5=1.96133`
requires the separate field-domain premise `sigma_w>=1/5`; an 80-degree
inclination premise alone gives `g cos(80 deg)~=1.7029069`.

**Noncircular nominal-AW reduction.** The complete reader shows that G0
cannot be used to prove its own nominal-AW premise.  Decomposing attitude by
the field direction removes most of the apparent box loss: a field-axis
rotation of angle theta changes gravity by at most
`2 g sigma_w sin(theta/2)`, only 0.2052962 m/s^2 at theta=6 deg and
sigma_w=1/5.  With the declared post-projection BA error and fast
accelerometer residual, the static kernel nuisance is 1.1304628 m/s^2,
below g/5.

The remaining translated physical term is exactly
`mean(M_b a)=[M_b v]_0^T/T-T^-1 int dot(M_b)v dt`, so its noncircular
control requires a signed temporal/TV bound on the SHIPPING field-axis
attitude-error loop.  The retained 6-degree angle bound, Lemma I* net
rotation, axial gyro-bias projection sector, and physical BA-rate bound do
not individually supply such a TV bound.  Using G0 to obtain it would be
circular because `f_hat parallel b` is precisely the boundary
`|a_hat x b|=g sigma_w`.

Equivalent target: prove a constrained minimum source-action inequality
`A_min(g sigma_w)>A_available` for the literal axial-bg + BA + AW/S
chronology.  OU leakage alone is far too weak: at tau=12 s replenishing a
g/5 DC AW component costs only 0.000654--0.000980 m/s2 per 4--6 ms step.
The weak-regularizer corner must also be retained: tau=12,sigma=4 drives the
SpectralMSE target into the r_S=100 m*s clamp with T_S=0.15 s.  Thus the
next analytical obligation is a reachability/passivity bound for the coupled
field-axis error/AW/S/bias loop, not another observability lemma and not a
stronger MARINE MOTION assumption.

**Axial coupled-loop reduction (PR #630).** A 32-s noncircular
field-axis calculation now gives a quantitative necessary condition for the
nominal-AW pathology.  At sigma_w=1/5 and the retained 6-degree tilt,
field-axis gravity mismatch is only 0.2052962 m/s2.  Adding the universal
post-projection BA error, fast accel residual, physical signed-mean/sampling
bound and the 0.02001 rad/s free gyro/bias rotation leaves a strict
0.2270622 m/s2 gap to g/5.  Thus any forbidden trajectory must obtain at
least that much rectification from correction-induced field-axis attitude
motion; the conservative equivalent average correction-induced angular
variation is 0.0412840 rad/s.

Total NIS/source action cannot close this gap: deterministic bounded sensor
residuals may be coherent, so square-summed measurement action grows with
word length.  The correct functional is the signed velocity/axial-correction
pairing.  For each correction,
`|delta theta_b|^2 <= NIS * q_b' K Omega K' q_b`, and the second factor is
exactly the Joseph decrement of axial attitude covariance.  The desired
certificate telescopes those decrements against prediction replenishment
before taking a norm.

A new structural obstruction is explicit: no source-uniform AG/axial
covariance ceiling is currently available without G0.  On the exactly
field-axis-degenerate nominal-force branch, magnetic and accelerometer rows
can both miss the axial attitude/bg pair, while bg process noise accumulates.
Therefore a generic Joseph-decrement bound cannot be promoted uniformly
without first excluding persistent degeneracy; using G0 for that exclusion
would be circular.  This is not a reachable counterexample.  It means the
next contradiction must come from the coupled AW/BA/S MEAN recursion itself,
showing that persistent `f_hat parallel b` cannot be maintained by the
declared physical/source histories.  Only after that escape-from-degeneracy
lemma may the AG covariance/information machinery be invoked.

**Persistent exact-degeneracy bridge CLOSED.** On the exact
zero-innovation branch `f_hat=a_hat_w-g e_z parallel b`, projection normal
to the committed field gives identically
`P_perp a_hat_w=P_perp(g e_z)`, magnitude `g sigma_w`.  Averaging the
literal projected accelerometer identity and using bounded physical velocity,
jerk sampling fidelity, the sharp field-axis gravity defect
`2 g sigma_w sin(theta/2)`, the universal post-projection BA-error bound and
the commissioned fast accelerometer residual yields

`g sigma_w <= 2 Vmax/T + J hmax/4
 +2 g sigma_w sin(theta_max/2)+B_ba,post+N_a`.

At sigma_w=1/5 and theta_max=6 deg the denominator margin is
0.6808672329 m/s2, so `Tcrit=2 Vmax/margin=16.1561 s`.  Every exact
zero-innovation degenerate interval of 17 s is therefore impossible, with
about 0.0338 m/s2 strict margin.  This uses no G0, AW covariance ceiling,
S-gain sign or carried AW statistic.

For a nonzero-innovation near-degenerate branch the same identity adds only
the SIGNED transverse accelerometer-innovation mean `Rbar_acc`.  On 17 s
the branch is excluded if `Rbar_acc<~0.0338 m/s2`.  The pointwise 0.3 m/s2
fast-residual box cannot supply that signed bound because deterministic
coherent residuals need not average away.  Thus the remaining robustification
is now: either prove the complete-reader/S-chain signed innovation mean below
this escape margin, or use a larger signed innovation directly as
measurement information/action forcing departure from degeneracy.  G0 may be
invoked only after this dichotomy; doing so before would be circular.

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
and qualitative four-S injectivity closes the zero-action LIN/AW nullspace.
The current physical assumptions do not themselves provide a proved
quantitative transversality bound for the remaining tilt/BA compatibility
sector. Algebraic tangent/near-tangent configurations exist, but subsequent
reachability analysis has **not** proved that literal shipping A21 executions
can approach them arbitrarily closely.

Conclusion: a source-uniform finite-horizon detectability constant C_det is
not proved by the current lemmas, but the present assumptions have **not been
disproved**. It remains open whether the literal shipping dynamics restrict
the reachable A21 histories enough to supply the missing quantitative gap.
A different theorem may avoid this gap; alternatively a new physical/service
condition could supply it, but such strengthening is not justified by the
reachability work completed so far. Adding such a condition
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

4. The later "Literal whitening and anchor identification" section formerly
   resurrected the false unprojected leverage inequality. It is corrected to
   use only `Xg_perp=(I-P_U)X_g`. The missing `P_U X_g` contribution remains
   OPEN; no AG-process bound follows from the bookkeeping identification.
5. Exact construction of `G_red(c)` is not a proof of the target uniform
   Loewner floor `G_red(c;history)>=G_*(c,r)>0`. That direct-information
   inequality remains an OPEN O1 target wherever it is invoked. Later
   rank-continuous/ordered-eigenvalue results must be cited explicitly if
   used instead; they do not retroactively prove the earlier DI-6 claim.

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
## Exact compatibility manifold: literal mean recursion

AW covariance synchronization does not change the AW mean. The mean
compatibility dynamics therefore consist of OU prediction plus actual
accelerometer/S corrections (and BA OU prediction/corrections when enabled).

In world coordinates, for the zero-BA compatibility subcase define
`fhat=a_hat-g=lambda b` immediately before an accelerometer correction.
Across a prediction of duration h with AW factor phi=exp(-h/tau),

`fhat^-_next=phi fhat +(phi-1)g`.

Hence its transverse component is

`P_b fhat^-_next=(phi-1)P_b g`.

Unless b is parallel to gravity or phi=1, the exact compatibility manifold is
not invariant under prediction. To return to it at the next accelerometer
epoch, the intervening mean corrections must supply exactly

`P_b Delta a_hat = (1-phi)P_b g`

plus the known contributions of any S correction and the nonzero-BA affine
compatibility term.

At an accepted accelerometer update the AW mean changes by
`Delta a_hat=K_aw r_acc`, where K_aw is the AW 3x3 block of the literal
Kalman gain. Therefore local exact reachability of the compatibility
manifold requires the transverse control-rank condition

`rank(P_b K_aw)=2`

at the relevant corrections, together with a residual r_acc whose implied
physical measurement remains inside MARINE MOTION / IMU BIAS bounds. S
updates add their literal AW gain times the S residual and must be included
in the same affine cycle equation; covariance sync adds no mean term.

This identifies a sharp reachability criterion but does not yet prove it.
The current proof has no source-uniform lower singular-value bound for
`P_b K_aw`; cross covariance can in principle make that block singular.
Conversely, no invariant forces it singular. Thus exact-compatible MOVING
reachability is reduced to the actual closed-loop gain-rank problem, not to
AW sync or autonomous OU dynamics.

A constructive existence proof can proceed from any strict-margin regular
state where `rank(P_b K_aw)=2`: the required transverse correction is O(h),
so by continuity sufficiently small h gives a small residual; physical
motion/service inequalities with strict margins persist under the resulting
small smooth input perturbation. To make this rigorous one still must exhibit
one reachable strict-margin state with that gain rank and verify the
longitudinal/S/BA cycle closure. No current analytical certificate supplies
that base state, so exact-compatible MOVING + strict MAGNETIC SERVICE remains
unresolved rather than ruled out.
## Full-rank AW correction base state

At the constructor/diagonal covariance state, all AW cross-covariances are
zero and `P_awaw=Sigma_aw_stat>0`. For an accepted accelerometer correction
the AW gain block is therefore exactly

`K_aw=Sigma_aw_stat R_wb' S^-1`.

`Sigma_aw_stat` and the innovation covariance S are positive definite and
`R_wb` is orthogonal, so K_aw is invertible. Hence for every nonzero magnetic
direction b,

`rank(P_b K_aw)=rank(P_b)=2`.

Thus the local two-component compatibility control rank exists analytically;
gain rank itself is not an obstruction.

This constructor covariance is not automatically a recurring A21 base state.
To use it for the exact-compatible MOVING counterexample one must connect it
to a strict-margin regular execution without invoking a reset/reinitialization
inside A21. A nearby diagonal-dominant covariance would suffice because rank
two is open. The remaining task is therefore to exhibit a reachable regular
covariance neighborhood with K_aw invertible (or prove that ordinary
prediction/corrections preserve invertibility long enough), while arranging
the mean cycle.

For the mean cycle, the transverse implicit equation at each accepted
accelerometer update has derivative P_b K_aw of rank two. The implicit
function theorem therefore solves the two transverse innovation components
locally for the O(h) OU drift. The longitudinal innovation remains free and
can be used together with S residual/control and, when enabled, BA innovation
to satisfy the scalar/integral closure conditions. A rigorous periodic or
finite-word construction still needs those longitudinal/S/BA equations and
the physical measurement realization checked against MARINE MOTION and IMU
BIAS.
## Regular-A21 AW-gain entry lemma: current result

For an accepted accelerometer correction the literal AW gain numerator is

`N_aw = P_aw,theta J_att' + P_aw,aw R_wb' + P_aw,ba`

(plus the AW/gyro-bias lever-arm cross term when that feature is active), and
`K_aw=N_aw S^-1`. Since S is positive definite, transverse gain rank is the
rank of `P_b N_aw`.

A positive AW marginal alone does not imply this rank. PSD covariance permits
cross blocks to cancel the direct `P_aw,aw R_wb'` term on one or more
transverse directions. Thus the constructor proof cannot be extended to all
A21 roots from the AW covariance floor alone.

Shipping prediction injects a favorable fresh term: the LL process covariance
adds `Q_aa=(1-phi^2) Sigma_aw_stat` to the AW marginal while adding no fresh
AW-attitude or AW-BA cross covariance. Hence immediately after prediction

`N_aw^- = N_aw,inherited^- + Q_aa R_wb'`.

A universal rank-two floor would follow if the fresh transverse singular
value exceeded the inherited cancellation norm. At regular 4--6 ms steps,
however, `1-phi^2=O(h/tau)` is small. Existing nuisance/cross-covariance
bounds do not prove

`sigma_min(P_b Q_aa R_wb') > ||P_b N_aw,inherited^-||`

or any signed variant preventing exact cancellation. Therefore no
source-uniform delta_K or finite entry time T_K is currently derivable from
the proved covariance bounds.

Conversely this algebra does not construct an actual rank-deficient A21
execution: PSD-compatible cancellation at one covariance matrix is not enough;
the matrix must be reachable under the literal Riccati recursion and strict
MAGNETIC SERVICE. No such execution has been analytically constructed.

Conclusion: the proposed regular-A21 AW-gain entry lemma is presently neither
proved nor disproved. Its controlling subproblem is covariance reachability:
can the actual Riccati recursion reach/approach the algebraic cancellation
manifold `det(P_b N_aw|_bperp)=0` under strict service? This is a lower
dimensional covariance invariant/reachability question. Until it is settled,
the exact-compatible IFT counterexample has no certified reachable base.
## Transverse AW-gain cancellation manifold under Riccati updates

Let C be the literal accepted accelerometer Jacobian at a fixed pre-correction
state and N=P C' its full gain numerator. The exact Kalman/Joseph covariance
update satisfies

`P+=P-P C' S^-1 C P`, `S=C P C'+Racc`,

and therefore

`P+ C' = P C' [I-S^-1 C P C']
        = P C' S^-1 Racc`

(equivalently with the invertible right factor written in the matching
order). Since S and Racc are positive definite, this right factor is
invertible. Hence for the same row C the rank of every row-block of `P C'`,
including the AW block, is preserved by its own accelerometer correction.
An accepted accelerometer update cannot create exact AW gain-numerator rank
loss from a full-rank pre-update AW numerator.

Literal attitude reset is an invertible covariance congruence and likewise
cannot create rank loss merely by coordinate change. AW covariance sync adds
a PSD AW-only increment on the default path; prediction adds the fresh
`Q_aa=(1-phi^2)Sigma_aw` AW term but also changes C through the nominal mean
and transports inherited cross covariance. S and magnetic corrections use
different rows and can change `P C_acc'` nontrivially.

Therefore the cancellation manifold is not invariant under the full A21
cycle, but neither is it reachable through an accelerometer correction alone.
Any approach to
`det(P_b N_aw|_bperp)=0` must be generated between accelerometer updates by
prediction, S correction, magnetic correction, BA mode changes, or the
change of the next accelerometer Jacobian C_acc itself.

Strict MAGNETIC SERVICE constrains the magnetic corrected rows but does not
bound their induced AW cross-covariance action relative to the next
accelerometer row. The current covariance bounds likewise do not supply a
positive distance from the cancellation manifold after those intervening
operations.

Thus exact rank loss is excluded across a single accepted accelerometer
correction when starting full rank, but the complete regular-A21 Riccati
reachability question remains open. To close it one needs a per-operation
distance-to-singularity inequality for prediction + S + magnetic operations,
or an invariant sign/determinant property of the 2x2 transverse numerator.
No such invariant is currently proved.
## Per-operation distance-to-singularity inequalities for the AW gain numerator

Fix the future accepted accelerometer row C_a and define its AW gain numerator
`N_a=(P C_a')_aw`. Let `d_a=sigma_min(P_b N_a|_bperp)`.

For any covariance correction with row H, innovation covariance
`S_H=H P H'+R_H`, the exact update is

`P+=P-P H' S_H^-1 H P`.

Against the future accelerometer row this gives

`N_a+=N_a- (P H')_aw S_H^-1 H P C_a'`.                    (OP-1)

Hence Weyl gives the exact safe inequality

`d_a+ >= d_a- - ||P_b(PH')_aw||
                    ||S_H^-1/2 H P C_a'||
                    ||S_H^-1/2||`.                          (OP-2)

Equivalently retain the linked correction matrix itself for a sharper bound:

`d_a+ >= d_a- - ||P_b(PH')_aw S_H^-1 H P C_a'||`.          (OP-3)

For H=C_a (the accelerometer's own row), the special identity proved above
replaces this subtraction: the AW row block is right-multiplied by an
invertible matrix, so exact rank is preserved.

For an S correction, H=E_S'. Therefore

`Delta N_a,S= -P_aw,S (P_SS+R_S)^-1 P_S,* C_a'`.            (OP-S)

For a magnetic correction,

`Delta N_a,M= -(P H_m')_aw S_m^-1 H_m P C_a'`.             (OP-M)

Strict MAGNETIC SERVICE lower-bounds cumulative magnetic information in its
heading/axial-bias root coordinates, but it gives no upper bound making
`||Delta N_a,M||<d_a`. Thus service alone does not prevent crossing.

Prediction has

`P-=F P+ F'+Q`.

For the future row C_a^- its AW numerator is

`N_a-= (F P+ F' C_a^-')_aw + (Q C_a^-')_aw`.               (OP-P)

The second term contains the favorable fresh AW block
`Q_aa R_wb'`, but Q also has correlated LIN/AW blocks and the first term
contains transported inherited cross covariance. Therefore

`d_a- >= sigma_min(P_b Q_aa R_wb'|_bperp)
       - ||P_b R_P||`,                                      (OP-second)

where R_P is the exact sum of all other transported/process contributions.
Current bounds do not make the right side positive.

AW covariance sync on the default path is an AW-only PSD increment Delta.
Against a fixed future accelerometer row it changes

`N_a -> N_a + Delta R_wb'`,                                 (OP-AW)

so

`d_new >= sigma_min(P_b Delta R_wb'|_bperp)-||P_b N_a||`

or, locally, `d_new>=d_old-||P_b Delta R_wb'||`; neither
inequality forbids a determinant crossing because adding a positive matrix
before an unrelated rotation/cross term is not sign preserving for the
2x2 transverse determinant.

A literal reset is an invertible covariance congruence, but the next
accelerometer row changes with the reset/mean attitude. For a fixed physical
row this is a coordinate transformation and preserves rank. The dangerous
piece is the change in nominal specific force, hence in J_att. If
`C_a,new=C_a,old+Delta C`, then

`N_new=N_old+(P Delta C')_aw`,
`d_new>=d_old-||(P Delta C')_aw||`.                         (OP-C)

These identities settle the structural question: none of prediction, S
correction, magnetic correction, AW sync, or Jacobian change has a
sign/determinant invariant that follows from PSD and strict MAGNETIC SERVICE
alone. Each can alter the transverse numerator by an additive matrix, and
the current assumptions provide no bound smaller than the incoming distance
`d_a`.

This does not yet exhibit a reachable crossing, but it rules out proving the
AW-gain entry lemma from per-operation rank preservation. A source-uniform
rank invariant would require a new quantitative dominance estimate on the
linked OP-S/OP-M/OP-P/OP-C terms. No such estimate is present in the current
proof assumptions.
## Can S/magnetic corrections cross the AW-gain cancellation manifold?

For a magnetic correction H_m=[J_m,0,...], the covariance update gives

`P_aw,*+ = P_aw,* - P_aw,theta J_m' S_m^-1 J_m P_theta,*`.

Thus its change of a future accelerometer AW numerator is

`Delta N_a,M = -P_aw,theta J_m' S_m^-1 J_m P_theta,* C_a'`.

If `P_aw,theta=0`, a magnetic correction cannot move N_a at all. Nonzero
AW-attitude correlation is therefore necessary. Accepted accelerometer
corrections generically create such correlation, so an acc->mag sequence is
the minimal candidate crossing mechanism.

However the covariance part of a Kalman correction is independent of the
measurement residual. For fixed pre-correction covariance and fixed magnetic
Jacobian/noise, there is no continuous residual parameter with which to tune
Delta N_a,M: varying the magnetic measurement changes the mean/reset and hence
future C_a, but not the pre-reset Joseph covariance map itself. Similarly the
S covariance correction has no residual parameter at all (`r_S=-S`) and its
covariance map is fixed by P and R_S.

Continuous crossing can therefore only be tuned through quantities that
change the pre-update covariance/Jacobians/noise along an actual execution:
elapsed prediction time, tuner parameters within clamps, attitude/field
geometry, or preceding correction chronology. Strict MAGNETIC SERVICE
constrains the accumulated transported magnetic rows but leaves these
parameters continuous.

At the algebraic PSD level, no sign invariant prevents crossing: choose a
positive-definite covariance with nonzero P_aw,theta and vary the strength of
a magnetic conditioning continuously from zero to its shipping value; the
rank-one/2-D downdate moves N_a continuously and can be arranged to cancel a
chosen transverse component while the Joseph covariance remains PSD. But
shipping does not expose magnetic-update strength as a free continuous
parameter: Rmag is fixed/adapted only through allowed configuration and the
update either occurs with its literal strength or not.

Therefore an algebraic crossing family is not yet a reachable shipping
family. To refute the universal entry lemma one must realize the required
conditioning-strength continuation through an allowed execution parameter
(most naturally the continuously varying pre-update attitude/cross covariance
generated by prediction and accelerometer history) while retaining strict
service. No existing reachability theorem provides that continuation.

Conclusion: PSD/Joseph structure does not forbid determinant crossing, and
magnetic/S updates can supply the necessary additive term, but actual
shipping reachability of a crossing remains unproved. There is no hidden
Riccati sign invariant found here; the blocker is again reachability of the
required cross covariance/Jacobian family.
## Attempted continuous sign-crossing family for transverse AW gain

Parameterizing one accepted accelerometer innovation alpha before a magnetic
correction gives a genuine continuous family of shipping executions as long
as acceptance/projection branches remain unchanged: the covariance Joseph
map at that accelerometer update is residual-independent, but the injected
attitude/reset, subsequent Jacobians, magnetic covariance update and next
accelerometer row depend continuously on alpha. Therefore

`D(alpha)=det(P_b N_aw(alpha)|_bperp)`

is continuous on such a branch.

Continuity alone is insufficient. A nonzero derivative D'(0) only proves
local variation, not opposite signs. The level/north/collinear geometry does
not provide an odd symmetry `D(-alpha)=-D(alpha)` because the Riccati
covariance path and reset Jacobians contain even and mixed terms. Existing
bounds also do not give a derivative lower bound large enough to force a
crossing before a gate/projection/service margin is reached.

The documented 25-Hz collinear carried word is a useful near-degenerate
diagnostic with dense magnetic corrections, but it reports a small positive
same-cell singular value, not the sign of this 2x2 determinant. It cannot be
promoted to an actual reachable sign-crossing theorem, and finite replay is
non-promoting under the research protocol.

Hence no continuous actual-execution family with rigorously opposite signs
has yet been constructed. Conversely no sign invariant was found: the exact
per-operation formulas permit additive changes capable of algebraic crossing.

The remaining exact condition for an intermediate-value disproof is now:
find one regular strict-service branch and two analytically certified
parameter values alpha_-<alpha_+ on that same branch such that
`D(alpha_-)D(alpha_+)<0`, with all MARINE MOTION/IMU BIAS/service inequalities
proved throughout the interval. This requires a signed determinant formula
or monotonicity estimate for the complete acc->reset->mag->prediction map;
norm bounds and singular values cannot establish it.

Until such a signed formula is derived, the universal AW-gain entry lemma is
neither proved nor refuted by reachability. The current proof should not
claim an intermediate-value crossing.
## New controlling route: soft-kernel ordered-eigenvalue compactness

The exact-kernel quotient/detectability route is no longer preferred because
its quotient changes discontinuously when a one-dimensional compatibility
kernel disappears. The corrected alternative keeps finite prior precision on
a continuously chosen least-information tilt/BA direction for every word.

The key spectral quantity is the second ordered eigenvalue lambda_2 of the
complete-word slow information (after the already justified nuisance
elimination), not lambda_min^+. Ordered eigenvalues are continuous. Therefore
if the complete-word nullspace theorem truly gives nullity <=1 on every
element of a compact closed retained word class, then lambda_2>0 pointwise
and compactness legitimately yields `inf lambda_2>0`, including through
rank-changing words.

This route avoids the epsilon^-2 kernel-disappearance pathology because the
finite rank-one prior remains on the weak direction before and after exact
rank loss. It also avoids requiring the unresolved AW-gain entry lemma.

Before promoting this to a theorem, compactness must be checked carefully:
word duration, dt/tuner/bias/state coefficients are bounded; however event
acceptance and scheduler patterns are discrete. Treat each regular event
pattern as a closed stratum and prove there are finitely many patterns on the
fixed word horizon, or include boundary patterns explicitly. J must be
continuous on each stratum. Any limit in which an applied correction becomes
rejected belongs to a neighboring stratum and must separately retain the
nullity<=1 conclusion under MAGNETIC SERVICE.

Thus the next exact task is topological rather than another covariance
reachability construction: prove the finite closed-stratum compactness and
nullity<=1 on every stratum closure. If it closes, the first legitimate
source-uniform positive spectral modulus follows nonconstructively as
`lambda_2,bar=inf lambda_2>0`; quantitative extraction can follow afterward.
## Controlling proof order after PR #625 consolidation

The controlling path is now:

`complete-word nullspace (Theorem A)`
` -> uniform qualitative O1 via fixed-factor lower semicontinuity`
` -> constructive g_under(c,r)`
` -> K(c,r)`
` -> same-history D(c,r)<=c`
` -> rho_0<1`
` -> nonlinear retained-radius inequality`
` -> every-prefix retention / H18-release entry / regime transitions / float32`.

The observation-only closed-range lemma is false across nuisance rank loss
(e.g. diag(1,epsilon)); it is not an obligation. The fixed-factor
nuisance+process variational operator is controlling because proved action
floors bound nuisance minimizers and provide lower semicontinuity.

Theorem A and the compactness contradiction now give a qualitative uniform
`g_mu>0` for every finite mu=1/c. This supersedes the need to use the
detectability, AW-gain reachability, determinant-sign, or soft-kernel routes
for O1. Those sections remain research history unless separately needed.

Next mathematical obligation: extract a constructive positive lower enclosure
`g_under(c,r)` from the literal zero-action implications while preserving
same-history correlations. Do not return to observation-only pseudoinverses
or the retracted unprojected leverage inequalities.
## Reachable-base dependency for the persistent-pair construction

The requested base cannot currently be certified from startup. The stability
proof explicitly leaves the H18/reference-refinement/BA-release map open:
actual release entry into the A21 retained region is a later obligation.
Therefore constructor covariance or a nominal A21 record cannot be promoted
to an **actual reachable recurring A21 base** without solving that release
obligation out of order.

The recurring A21 theorem is conditioned on coefficients from an actual
post-release execution. Consequently a counterexample/persistent-pair base
must likewise occur on an actual post-release history; an arbitrary frozen
state satisfying local A21 inequalities is insufficient.

Properties (1) strict physical/service/gate margins and (3) rank-two
transverse control are open properties once a suitable actual base exists;
property (2) exact compatibility is a closed codimension condition that the
IFT construction can preserve locally. Physical realization and strict
MAGNETIC SERVICE have already been shown locally open. But no current theorem
guarantees that the open post-release reachable set intersects the exact
compatibility manifold at a rank-two point.

Thus reachable-base existence is logically downstream of, or coupled to, the
still-open release reachability map. It cannot be solved from the current A21
tail lemmas alone. Conversely, failure to exhibit such a base does not prove
the compatibility manifold unreachable.

Proof discipline consequence: do not use the persistent-pair construction to
claim alpha_bar=1, and do not use its absence to claim alpha_bar<1. Complete
the A21 O1/O2 theorem conditionally on the actual retained A21 class first;
startup/H18/release entry remains a separate later composition obligation as
specified by the controlling proof order.
## Review corrections: compactness and controlling O2 lemma

The statement that finitely many possible matrix ranks yield finitely many
closed constant-rank strata is false and is withdrawn. Exact-rank sets are
not generally closed. O1 continuity now uses only the fixed-factor
variational action with nuisance/process coercivity; no nuisance-rank
pseudoinverse stratification is controlling.

The ordered second eigenvalue route is valid only with this continuous
quadratic-form representation: complete-word nullity<=1 on the closed compact
admissible class plus continuity of J gives `lambda2_bar>0` at the existence
level. No numerical lambda2 floor is claimed.

O2 must use a fixed positive physical root metric M. Compatibility generators
are normalized by `nu'Mnu=1`, and the exact same-history return is
`a_W=nu_+' M T_W nu_W`. Euclidean mixed-unit comparisons with one are not
theorem statements.

The controlling proof order is now:
`nullity<=1 -> lambda2_bar>0 -> soft/augmented O1 finite -> metric-normalized
compatibility return |a| (or finite product) <1 -> D<=c -> rho0<1 ->
nonlinear radius -> prefix/release/regime/float32`.

Exploratory detectability, AW-gain and sign-crossing sections are research
history and must not be cited as alternate controlling routes.
## Shipping-invariant audit of zero-dynamics continuation margins

Three candidate global margins were audited against shipping.

**Accelerometer gate margin.** In `measurement_update_acc_only` NIS is
computed for diagnostics but is not used as a rejection threshold. A finite
accelerometer sample is rejected only if the 3x3 safe LDLT fails. In real
arithmetic `S_acc=CPC'+Racc` is SPD on PSD covariance because `Racc>0`, so
the regular mathematical trajectory has no independent NIS/gate margin that
must decay. Thus accelerometer acceptance is not the global escape mechanism
(float32 LDLT totality is a later arithmetic obligation).

**Transverse authority D_perp.** This is generated by the actual correction
gain `K=P C' S^-1` composed with the compatibility-output derivative. Since
`S^-1` is nonsingular, loss of transverse authority comes from the relevant
projected rows of `P C'` and from compatibility geometry. AW process
covariance injection/sync supplies positive AW marginal variance but existing
cross covariances can algebraically cancel its contribution to `P C'`.
No shipping covariance floor or projection currently supplies a lower bound
on the projected gain numerator. Therefore no existing invariant proves
`sigma_min(D_perp)>=d0>0` along constrained trajectories.

Nor is there an invariant forcing D_perp to lose rank: accepted
accelerometer corrections preserve the rank of their own gain-numerator row
block under the same Jacobian, prediction injects fresh AW variance, and
literal resets are invertible coordinate changes. S/magnetic corrections and
the changing next accelerometer Jacobian can alter the numerator, but no
monotone determinant/sign law was found. Thus D_perp is genuinely undecided.

**MAGNETIC SERVICE margin.** The proof assumption supplies the non-strict
closed floor `lambda_min(G_M)>=mu_M` for actually applied magnetic service
windows. It does not supply a uniform surplus `>=mu_M+delta`. Consequently
strict service is locally open around a strict base, but an infinite
constrained trajectory is allowed by the assumptions to approach the service
boundary while still remaining admissible. There is no shipping invariant
forcing either a positive surplus or finite-time service failure.

Conclusion: after removing the spurious NIS-gate concern, the only substantive
global continuation quantities are D_perp rank and magnetic-service surplus.
Current shipping invariants control neither in the direction needed to decide
forward completeness versus escape. This is the exact residual O2 gap.
## Covariance/gain self-consistency of the persistent-mode construction

For accepted linearized Kalman corrections, the Joseph covariance recursion
depends on P, the measurement Jacobian and R, but not on the numerical
innovation residual. The S pseudo-update has the same property. Prediction
covariance likewise does not depend directly on physical translational
acceleration.

Hence physical acceleration chosen to realize the compatibility-maintaining
innovation changes the mean but does not directly change P or K. It couples
back only through mean-dependent Jacobians/resets, tuner variables and event
chronology.

There is no known covariance invariant forcing this coupled compatible
execution to escape: PSD is preserved, process covariance is injected, and
resets are invertible transports. Conversely no recurrent compatible
covariance orbit is yet proved. The remaining existence problem is a joint
Poincare/viability fixed point for nominal geometry+tuner+covariance with
D_perp nonsingular; physical residual realization is then supplied by the
already-derived PRDC equation.

Thus no residual-to-covariance contradiction closes O2. The unresolved
estimator-internal question is recurrence of the compatible nominal
geometry/covariance branch, not innovation magnitude.

## Periodic Riccati recurrence on a compatible coefficient cycle

Condition on a periodic compatible nominal geometry/event pattern and fixed
periodic tuner parameters. The covariance recursion is then a finite-period
discrete Riccati/Joseph map with positive measurement-noise floors and the
shipping process-noise injection. The complete-word nullity/detectability
results control all but the declared one-dimensional compatibility direction;
the added finite kernel prior used in O1 regularizes that direction for the
auxiliary Riccati diameter.

For the **actual** covariance recursion, standard periodic Riccati existence
cannot be invoked blindly because the physical compatibility direction may be
undetected and process noise may enter it. A bounded periodic covariance orbit
exists only if that neutral direction is dynamically stable or receives
sufficient recurring information. This is precisely O2. Therefore using a
periodic-Riccati theorem here would be circular.

The Poincare existence question cannot be reduced to O1. If an exact
compatibility mode is unit-persistent and receives process covariance, actual
P grows along it and no recurrent covariance orbit exists; if its deterministic
return is strictly below one, a bounded periodic covariance orbit is possible.
Thus covariance recurrence is mathematically equivalent to the scalar O2
return already under investigation.

Consequently no independent fixed-point theorem solves the blocker. The
compatible mean/tuner cycle plus D_perp nonsingularity is not enough; one must
also know the compatibility-mode covariance return. Conversely, assuming a
recurrent covariance orbit would assume the conclusion needed for O2.

This prevents a circular counterexample construction. A genuine persistent
physical/mean compatibility execution could coexist with **unbounded
covariance** in the invisible mode; that would refute the desired covariance
stability theorem even more directly, without requiring a recurrent P orbit.
Hence a recurrent covariance orbit is not necessary for an O2 counterexample.

The sharper counterexample target is therefore only a forward-complete
admissible mean/physical compatibility execution with D_perp nonzero. Along
it, either P stays bounded (recurrent/bounded counterexample) or P grows in the
unit-persistent mode (direct failure of the claimed uniform covariance
ceiling). Requiring covariance recurrence was unnecessarily strong.

## Forward continuation with evolving covariance

Augment the compatibility-manifold state by the actual covariance P. At each
accepted accelerometer epoch the two transverse compatibility equations solve
for the two transverse physical-acceleration components whenever the literal
control Jacobian D_perp(P,x) is nonsingular. The longitudinal physical
acceleration remains free and may be chosen to satisfy the kinematic moment
conditions.

The covariance update is defined for every finite PSD P because each shipping
innovation covariance has a positive measurement-noise floor. Joseph updates
preserve PSD and prediction adds finite PSD process covariance. Thus finite
covariance growth does not cause a finite-time algebraic singularity of the
Kalman update. Large P may change the gain and D_perp, but it does not by
itself terminate the recursion.

A standard stepwise continuation argument therefore gives: a compatible
execution can be extended through every finite horizon as long as
(1) D_perp remains nonsingular at the required correction epochs,
(2) the MI-2-selected physical acceleration/jerk remains inside MARINE MOTION
bounds, and (3) the closed MAGNETIC SERVICE condition continues to hold.
No separate bounded-P hypothesis is needed.

The current shipping invariants prove none of these three margins must fail in
finite time. Conversely they do not provide positive lower margins sufficient
for a global continuation theorem. Hence evolving/unbounded covariance does
not resolve the existence question; it only enters through D_perp and the
mean-dependent geometry.

In particular there is no estimator-internal finite-time blow-up obstruction:
for every finite horizon on which D_perp and physical/service conditions hold,
the coupled mean/covariance recursion is well defined and the compatibility-
maintaining physical input is obtained recursively. Infinite continuation is
equivalent to avoiding the three boundary events above for all time.

Thus the remaining global blocker is again D_perp/physical/service viability,
not covariance magnitude. Since service equality is admissible and physical
realizability has no structural amplitude/jerk contradiction, D_perp is the
only unresolved estimator-internal continuation boundary.

## Minimal covariance data for K_S, K_a and D_perp

At an S pseudo-update, `H_S` selects the three S coordinates. Therefore the
gain is determined by the full covariance column block
`P_:S=P(:,S)` and the 3x3 block `P_SS`:
`K_S=P_:S(P_SS+R_S)^-1`.

At an accelerometer update, the Jacobian has nonzero columns only in attitude,
AW, BA when enabled, and gyro bias when lever-arm terms are enabled. Hence
`K_a` is determined by the combined column image
`P C_a'=P_:theta J_att'+P_:aw R_wb'+P_:ba+P_:bg J_bg'`
and the innovation covariance. D_perp additionally needs the relevant rows of
K_a and the subsequent mean/reset transport.

These are the minimal **readout** blocks, but they do not form a closed
Riccati quotient. Prediction propagates the full LIN block: the S and AW
columns mix with v and p through the exact 4-state chain. Thus updating
P_:S and P_:aw requires P_:v and P_:p. Attitude prediction/reset couples
theta with gyro bias; BA prediction carries its cross blocks. Once these
columns are included, Joseph measurement updates
`P^+=P-PC'(CPC'+R)^-1CP`
modify every retained column through products involving the corresponding
rows, which by symmetry are the same retained column family.

For the LIN sector closure therefore requires all columns
`P_:{v,p,S,aw}`; for AG it requires all theta/bg columns when gyro bias is
present; BA requires its columns. Their union is every state block in the
shipping OU-III covariance. Cross blocks among these groups are needed by
prediction and measurement updates. Hence the closure of the gain-readout
column set under the literal Riccati recursion is the full Pext covariance.

There is no exact lower-dimensional covariance quotient that determines
K_S,K_a,D_perp and is autonomous under shipping. One can compress to
measurement-space Schur/information quantities for a **single** update, but
their next-step evolution depends on cross covariances discarded by that
compression.

Therefore an exactly periodic gain sequence generally requires recurrence of
the full covariance (or a special symmetry/invariant submanifold that reduces
it). The generic reduced Poincare-map shortcut does not close.

A special symmetric counterexample could still exploit an invariant covariance
submanifold (for example axis-decoupled diagonal/block-diagonal covariance
under a specially chosen attitude/field geometry). Establishing such an
invariant submanifold is now the only route to a genuinely lower-dimensional
periodic covariance construction; otherwise full covariance recurrence is
unavoidable.

## Symmetric covariance-submanifold test

Test the strongest natural symmetry: choose world axes so gravity and the
geomagnetic field lie in a coordinate plane, use diagonal/isotropic AW process
covariance, diagonal sensor-noise matrices, and a periodic principal-axis
rocking attitude. Start from block/axis-decoupled covariance.

LIN prediction and its process covariance preserve per-axis block structure.
BA OU prediction and the S=0 pseudo-update also preserve it when their gains
inherit the same diagonal axis symmetry.

The symmetry fails generically at the attitude measurements. The accelerometer
attitude Jacobian is `J_att=-[f_cog]x`; for a nonzero field-parallel
compatible specific force its skew matrix couples the two axes perpendicular
to f_cog. The magnetic Jacobian likewise contains a skew matrix of the body
magnetic vector. With gravity and magnetic field noncollinear, no fixed
coordinate basis diagonalizes both skew-induced measurement geometries over a
nontrivial rocking cycle. Joseph updates therefore create cross-axis
attitude/AW/BA covariance blocks even from a diagonal starting P.

Quaternion reset/attitude transport then propagates those cross blocks. Thus
the axis-decoupled/block-diagonal covariance family is not invariant under the
complete shipping word except in degenerate constant-attitude or collinear
gravity/magnetic geometries, both excluded by MOVING/nonvertical-field
premises.

No useful exact symmetric covariance submanifold was found from the physical
rotational symmetries. Consequently the lower-dimensional periodic
counterexample shortcut is unavailable for the nondegenerate theorem class.

The remaining exact periodic construction must use the full covariance
Poincare map. Any fixed-point existence argument must therefore be on the
full finite-dimensional PSD covariance together with compatible mean/tuner
variables; symmetry cannot be used to reduce it without leaving the admitted
MOVING geometry.

## Four regular-period seed obligations: resolved status

1. **Smooth tuner/event branch.** The OU-III class exposes tau_aw, AW
stationary covariance and S cadence as model parameters/setters; the Kalman
class does not autonomously retune them from residuals. A fixed admissible
interior parameter history and fixed scheduler phase is therefore a legitimate
smooth coefficient branch for the filter-level A21 analysis (subject to the
orchestrator theorem carrying those values). S events are periodic and acc/mag
events can be chosen strictly accepted. This obligation is mechanically
available; no tuner fixed-point equation exists inside the filter class.

2. **Endpoint controllability.** After transverse compatibility elimination,
the remaining longitudinal physical acceleration is one scalar input per acc
epoch. Kinematic velocity/displacement closure uses its zeroth/first sampled
moments. The proof already has the continuous-history identity
`T^-1 || integral (a-g) x B dt || >= g B_h,min-2 V_max B_max/T`,
which shows the physical forcing has nontrivial signed vector action on long
windows. However it does not prove full-rank endpoint controllability of the
reduced nominal mean recursion: the input-to-(v,p,S,AW,BA,attitude) endpoint
matrix still contains actual Kalman gains/cross covariances. No existing
shipping invariant gives its full row rank. Thus obligation 2 remains OPEN.

3. **Tuner closure.** At filter level choose fixed interior model parameters
and a period commensurate with the S scheduler; then parameter/scheduler
closure is exact. If the higher-level orchestrator is included as part of the
theorem state, its adaptation law must separately be shown periodic/fixed.
The current filter proof cannot assert that external tuner closure. Thus this
obligation is CLOSED conditionally at filter level, OPEN for the full
orchestrated execution unless a fixed admissible tuner mode is certified.

4. **Covariance seed with D_perp!=0.** Constructor/diagonal positive covariance
gives nonzero direct AW accelerometer gain, but it is not a certified recurring
A21 seed. For an arbitrary SPD covariance, D_perp singularity is a proper
algebraic condition unless the determinant is identically zero. The direct AW
term shows it is not identically zero, so nonsingular SPD covariances form an
open nonempty set algebraically. What remains unproved is intersection of that
set with the actual A21-reachable covariance set. Thus obligation 4 is
algebraically solved but reachability remains OPEN.

Conclusion: all four cannot honestly be declared solved. The hard obstruction
is still obligation 2 plus reachable intersection in 4; dimension counting
cannot replace those theorems. The next decisive lemma is full-rank endpoint
controllability of the reduced longitudinal-input mean recursion for one
actual regular A21 covariance/history.

## Quantitative H18/refinement/release rocking robustness radius

Existing branch provenance already supplies a quantitative captured-domain
release result: for tilt error within approximately 7 degrees, the physical
magnetic bounds |B| in [20,75] uT, horizontal field >=15 uT and residual <=2
uT imply the deployed MagAutoTuner 35% norm and 5% horizontal gates.
MAGNETIC SERVICE then supplies finite 128-sample/30-s refinement completion,
followed by the internal 250-accepted-update/1-s guard and finite A21 release.

Therefore the currently certified **reference-refinement/release** robustness
radius is
`theta_rel = about 7 deg`
in tilt-error space, conditional on already being in that captured domain and
retaining the stated magnetic margins. This is not a source-uniform capture
radius from construction.

For the perturbative reachable-covariance argument, a rocking history can be
kept inside the certified release branch whenever its induced estimator tilt
error plus the base captured error remains below theta_rel. To satisfy the
MOVING excitation requirement simultaneously by a small perturbation, one
needs a strict budget
`theta_E < theta_rel - theta_base - theta_defects`.        (RR-1)
If theta_E is left symbolic with no upper bound, RR-1 cannot be certified.

Thus the release mechanism itself has a concrete nonzero neighborhood; the
missing comparison is an explicit theorem-domain choice/bound on theta_E and
a capture margin theta_base. If the theorem declares theta_E below the
available residual budget, continuity of the finite release covariance map
and openness of D_perp!=0 allow a current-domain rocking release history near
a regular released history. If not, the old fixed-attitude replay cannot be
perturbed far enough by the present certificate.

This cleanly separates two issues: reference refinement/release is
quantitatively robust on the ~7-degree captured domain; source-uniform
construction/capture into that domain under the amended MOVING assumption
remains open.

## Nonperturbative moving-capture target with symbolic theta_E

Keep theta_E>0 symbolic. The ~7-degree quantity is an estimator tilt-error
entry tube for the already-certified magnetic-reference refinement mechanism;
it is not a bound on physical rocking amplitude.

Use the fixed LOCAL GRAVITY vector g0 and near-constant nonvertical magnetic
vector b0 as two noncollinear world references. For a candidate execution whose
attitude estimate remains outside the 7-degree tilt tube, compare physical and
nominal vector records over a complete MOVING window. Bounded physical
velocity gives the exact diversity identity already present in the proof,
`T^-1 || integral (a-g) x B dt ||
 >= g B_h,min - 2 V_max B_max/T`,
with explicit eps_g/eps_B degradation. Jerk/sampling fidelity transfers this
continuous action to the sampled applied-record level, while MAGNETIC SERVICE
guarantees recurring actually applied magnetic information.

The missing finite-error lemma is not local rank. It must show that for every
attitude error outside the 7-degree tube, the joint sampled acc+mag innovation
action over a sufficiently long source-qualified window has a positive lower
bound after minimizing over admissible AW/BA nuisance histories. If such a
bound c_cap(theta_E)>0 holds, repeated windows cannot leave the estimate
outside the tube indefinitely while the retained covariance/state remains
finite: each window supplies finite positive correction information/action.
This gives a history-dependent finite capture time T_c(h,x0), after which the
existing refinement/release certificate applies.

The principal unresolved issue is nuisance absorption: nominal AW and BA can
shift the accelerometer prediction. The physical vector-diversity identity
prevents physical translation from supplying a permanent rotated-gravity
surrogate, but one must prove the estimator's OU/S/BA mean dynamics cannot
cancel the joint finite-error innovation action on every window. This is the
same physical-to-nominal bridge encountered in O2, now at finite attitude
error.

Thus nonperturbative moving capture is not yet proved by the existing
identities. The exact next lemma is a finite-error variational inequality:
`inf joint innovation action >0`
over all admissible same-history nuisance trajectories and all attitude errors
with tilt >=7 degrees, under symbolic theta_E, LOCAL GRAVITY and MAGNETIC
SERVICE. Once proved, capture-to-refinement/release follows without any upper
restriction on theta_E.

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
