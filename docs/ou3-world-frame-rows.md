# World-frame historical rows, injections and the aggregate six-column floor

These results enter `V_next <= rho V + c_d |d|^2` only through the historical
reader: the six-column floor (c, a, b0, or directly s) feeds B_*, J_AG, the
full covariance upper bound and rho_0<1. The injection bounds also enter the
injection-squared term of the finite-angle reset remainder in the nonlinear
supply eta. None of them supplies a source-uniform premise by itself. The AW
covariance ceiling (section 7) enters through Corollary A's pointwise AW
premise and `|e_aw|^2<=lambda_max(P_aw) V`; section 8 replaces that premise by
a signed nominal mean. Exact checks are in `world_frame.py`,
`aw_covariance_ceiling.py`, `aw_tracking.py`, `signed_injection.py` and
`aggregate_floor.py`; carried audits are `world_frame_source_diagnostic.py`
and `aw_tracking_source_diagnostic.py`.

## 1. World-frame factorization

Work in regular real-arithmetic A21: fixed deployment frame, zero lever arm,
constant committed reference B_w, no frame/relock/reference event. Let R_k be
the mean world-to-body rotation (`qref`) after operation k and
`[[A_k,B_k],[0,I]]` the literal chronological AG transport from the root.

* Prediction: mean `R+=Rq(-w h)R`; covariance `[[Rs,Bs],[0,I]]` from
  `rot_and_B_from_wt`. Rq is the rotation of the normalized
  `quat_from_delta_theta`.
* Applied correction: identity on the latent AG state; its reset has mean
  `R+=Rq(d)R` and covariance rows `G=I+[d]/2`.
* Rows: acc `-[R f]x` with `f=a_hat-g_c e_z` (world nominal force before the
  correction, configured gravity), mag `-[R B_w]x`; S has no AG column.

Put `A~_k=R_k' A_k R_0` and `B~_k=R_k' B_k`.

**Lemma W.** Every historical row block equals
`-R_k [f_k]x [A~_k R_0', B~_k]`. Hence `O=diag(-R_k) O~ diag(R_0',I)`, with
O~ stacking `[f_k]x [A~_k,B~_k]`; `sigma(O)=sigma(O~)`, and `LO=T_h` is
equivalent to a reader for O~ with terminal `[[A~_N,B~_N],[0,I]]`. All reader
actions are congruent by the orthogonal `diag(R_N,I)`. Moreover

* prediction: `A~+=W A~`, `B~+=W(B~+Gamma)`, `W=R+' Rs R`,
  `Gamma=R' Rs^-1 Bs`. If Rs equals the mean increment, W=I; with exact
  integral coefficients Gamma is the world integral of the nominal
  body-to-world rotation over the step;
* reset: `A~+=N A~`, `B~+=N B~`, `N=R+' G R`. For any orthogonal mean
  injection `N'N=I+(|x|^2 I-xx')/4`, `x=R' d`, so sigma_min(N)=1,
  `|N|=sqrt(1+|x|^2/4)` and `|N^-1|=1`. For the exponential injection,
  `N=Exp(-[x])(I+[x]/2)` acts transversally as `(1+i theta/2)exp(-i theta)`.

*Proof.* `[Rv]x=R[v]x R'` for rotations; substitute the definitions. The closed
form with exact Rodrigues/integral coefficients gives
`Rs^-1 Bs=int_0^h Exp([w]u) du`; the finite coefficient-series defect must
be retained for the shipping implementation. Finally `N'N=R' G'G R` and
`G'G=I+(|d|^2 I-dd')/4`.

Consequently the attitude estimate and its error enter the six-column array
only through the world injections x (via N) and through Gamma, which
integrates the nominal rotation. With no injections and exact
rotation/integral coefficients the rows are
`[f_k]x [I, int R_hat']`: the world observability matrix of
`theta_w(t)=theta_w(0)+(int R_hat')b_g`. The source mean uses the polynomial
quaternion for |w h|<.01. The covariance/integral coefficients use alternating
Taylor series through x^18 for |x|<1, x=|w|h, and closed forms otherwise.
The exact factorization retains the actual Rs and Bs in W and Gamma.
|W-I| includes the deterministic coefficient/mean-quaternion approximation
and floating-point defects. Ideal orthogonal-rotation inverse bounds do not
silently apply to the truncated coefficients: their transfer is part of the
still-open implementation/arithmetic premise, not certified by finite replay.

## 2. Attitude-invariant same-cell geometry

For the same-cell group `C=[-[R f]x; -[R_2 B_w]x G_local]`,
`C R=diag(-R,-R_2)[[f]x; [B_w]x N]` with `N=R_2' G_local R`. Thus sigma(C)
depends only on the world nominal force, the committed reference and the
local world injection, not on attitude. With `k=N^-1 B_w/|N^-1 B_w|` and
`kappa=|u'k|`, `u=f/|f|`, the retained projector bound has exact least
eigenvalue

`c^2 >= (F+B)/2-sqrt((F-B)^2/4+F B kappa^2) >= F B (1-kappa^2)/(F+B)`,

`F=|f|^2`, `B=|B_w|^2`, sharper than `min(F,B)(1-kappa)`. The line angle
between `N^-1 B_w` and `B_w` is at most
`psi=theta-atan(theta/2)+atan((1-rho)/(2 sqrt rho))`, `rho=(1+theta^2/4)^-1/2`,
`theta=|d|`, and `psi<=theta/2+theta^3/24+theta^2/(16-theta^2)`: N^-1 rotates
the transverse part by `theta-atan(theta/2)` and shrinks it by rho, which turns
a line by at most `max_s(atan s-atan(rho s))`. With `s_0` the f/B_w sine,
`c^2 >= F B (s_0-psi)_+^2/(F+B)`. The polynomial injection adds its rotation
defect.

**Generalized quiet floor.** On the zero-residual quiet record with any
constant nominal attitude and any committed reference with |B|>=20 uT and
horizontal fraction >=1/5, and g>=9: f=-g e_z, `s_0` is the horizontal
fraction and N=I, so `c^2>=81*400/481/25=1296/481>(8/5)^2`. `A~=I` and
`B~=t R'` with t>=4/125 give `s=16/635` and all six pivots `>=4/635` for
twelve rows. The quiet covariance action of `stationary_covariance.py` is not
regenerated for general references.

## 3. Literal injection budget

**Lemma I.** For an applied correction with innovation r, factored innovation
covariance S and literal gain K, the attitude injection `d=K_theta r` obeys
`dd' <= (r'S^-1 r) K_theta S K_theta'`. For the Kalman gain `K=PH'S^-1` with
`S-HPH'>=0` (safety increments included), Joseph gives
`P-KSK'=(I-KH)P(I-KH)'+K(S-HPH')K'>=0`, so `dd'<=NIS P_theta,theta` (prior).
In the retained region `sqrt(NIS)<=sqrt(V)+||n||_(R_eff^-1)` since
`H'S^-1 H<=P^-1`.

The injections bound `a=|A~|<=prod sqrt(1+|d|^2/4)`, the inverse-frame gyro
floor b0 and the `|d|^2|v|/6` reset term. The world transport is
`A~=I-[sum x]/2+O((sum|x|)^2)`: to first order the **signed** world injection
sum, not the sum of norms, controls its drift.

## 4. Nominal attitude columns from physical transverse force

**Corollary A.** On a window of complete cells of duration T with every
accelerometer correction applied, one applied magnetic row, reference
horizontal fraction `sigma_w`, and `|a_hat_i-a(t_i)|<=epsilon_a` in the
estimator world frame, the normalized world attitude-column Gram (trapezoidal
acc weights, magnetic weight one, zero injections) satisfies `G>=gamma I`,

`gamma=(sigma_w-e)_+^2/((1+e)^2+1)`,
`e=(2 V_max/T+J_max h_max/4+epsilon_a)/g_c`.

*Proof.* Convexity gives `sum alpha_i [u_i]x'[u_i]x >= [u_bar]x'[u_bar]x`.
The sampling-fidelity mean bound applied to physical a, plus epsilon_a, gives
`|u_bar+e_z|<=e`. The least eigenvalue of `[u]x'[u]x+[b]x'[b]x` is at least
`|u x b|^2/(|u|^2+1)`, and `|u_bar x b|>=sigma_w-e`.

No attitude, attitude-error or magnetic-residual term appears: the nominal
magnetic rows use B_w exactly. Injections lower `sqrt(gamma)` by at most
`beta sqrt(u_max^2+1)`, `beta>=|A~-I|`. At T=16 s and `sigma_w=1/5`, positivity
needs `epsilon_a<1.12383 m/s^2` (1.811319 as T grows);
`gamma(0)=12629938689/2094683134450`, about .00603. Through the covariance,
`|e_aw|^2<=lambda_max(P_aw) V`; section 7 bounds `lambda_max(P_aw)` sharply
and shows why that route cannot be uniform. Section 8 drops the pointwise
premise (Corollary A*); section 12 adds the gyro columns through the time
structure of Gamma.

## 5. Same-cell geometry depends on magnetic cadence

Take `B=(21,0,72)` uT, `b=B/75`, `c=g(e_z-(e_z.b)b)=g(-168,0,49)/625`,
`a(t)=c cos(2 pi t)`, `p=-c cos(2 pi t)/(4 pi^2)`, roll `sin(t/2)/100`, zero
biases. Exact checks give |a|<=2.746, jerk <=17.26, |v|<=.438, bounded
primitive and zero displacement mean, with horizontal fraction 7/25;
complete `T_E>=4 pi` windows have gravity span 1/50. This history satisfies
MARINE MOTION and IMU BIAS. At integer seconds `a-g e_z=-(24/25)g b`, so an
ideal same-cell world group there has the exact kernel b; at half-integers
the transverse force is `14g/25`. For a tracking estimator,
`sigma_min(C~)<=|C~ b_w|` is at most the AW tracking, gravity-mismatch,
reference and `|B_w||N-I|` terms.

With applied magnetic corrections only at integer seconds, every same-cell
group degenerates. That cadence is **not** admitted by MAGNETIC SERVICE: a
single correction per 1-s window maps both normalized heading and axial-bias
columns nearly onto `H_m R e_z`, and the carried run below has least service
eigenvalue 2.03e-5 against `mu_M=1`. So no refutation of the same-cell route
under all three assumptions is claimed.

**Jerk-limited collinear cadence.** If the physical force is parallel to b at
instants with gaps D_k covering a window of length L, then `c.a(t_k)=|c|`,
the jerk bound gives `c.a>=|c|-J min(t-t_k,t_(k+1)-t)` between instants, and
`|v(t_n)-v(t_0)|<=2V` yields

`sum D_k^2/sum D_k >= 4(g h-2 V_max/L)/J_max`, `h=|e_z x b|`.

At h=1/5 and L=16 s the floor is `127383/2500000` s, about .051 s, so a
regular 25-Hz applied cadence cannot be collinear at every correction. Tent
dips attain the bound. A same-cell floor therefore needs a coupling between
the density of applied informative corrections and this jerk lemma; neither
motion nor bias bounds alone supply it. Corollary A needs only one applied
magnetic row and avoids that coupling.

## 6. Carried-source diagnostic

The observer and driver are derived at run time from the existing ones; an
untapped control reproduces the terminal state exactly. Floors use the
exported row operands and charge each reset's float mean-injection angle
delta (the rotation between the exported post-reset attitude and the exact
`Exp([d]x)` injection, at most 1.8e-9 rad) on top of psi(theta). Quiet and
wave reuse the existing .32-s windows. Both collinear words apply magnetic
corrections at 25 Hz through startup; during 225--229 s the first applies
them only at integer seconds, the second keeps 25 Hz.

| Quantity | Quiet | Wave | Collinear 1 Hz | Collinear 25 Hz |
|---|---:|---:|---:|---:|
| same-cell sigma_min (actual) | 9.80665 | 8.74939 | .030129 | .00037063 |
| world floor / actual, max | 1.0 | .99999999975 | .99999997 | .99999999999 |
| nominal force/field sine, min | 1 | .894 | .00322 | 4.0e-5 |
| injection / NIS P_theta bound, max | 0 | .108 | .134 | .061 |
| inter-anchor gyro sigma: actual / NIS budget | .28 / .28 | .28 / .2774 | 3.0 / 0 | 3.96 / 0 |
| injection sum: actual / NIS bound (rad) | 0 / 0 | 1.4e-5 / .0317 | .0108 / 4.98 | .0099 / 6.52 |
| two-group s: actual / world+NIS budget | 1.9225 / 1.2043 | 1.7176 / 1.0656 | .0289 / 0 | .0189 / 0 |
| aggregate six-column sigma_min | 7.152 | 6.278 | 28.50 | 40.62 |
| AW tracking error max (m/s^2) | 0 | .0136 | .563 | .562 |
| lambda_max(P_aw): at syncs / min at acc rows | .0025 / .0012 | .0025 / .0010 | .0791 / .0026 | .0791 / .0026 |
| AW storage route ratio | 0 | .00035 | 6.68 | 6.91 |
| one-correction service lambda_min | 1.6e-83 | 1.4e-10 | 2.03e-5 | 5.5e-11 |

Prediction world discrepancy is at most 2.1e-9 (float branch), the reset Gram
identity holds to 6.4e-81 and row factorization to 1.7e-77. `|A~-I|` matches
half the signed world injection sum: .00112 on the collinear word versus a
norm sum of .0135. The AW block is reconstructed through every correction
(`P_aw-K_a S K_a'`) and prediction to 1.7e-6 relative; applied syncs are the
isotropic spectral max to 8.2e-7, the Lemma B step ratio is at most
1+8.2e-7, the ceiling ratio at most .999982, and each post-sync least
eigenvalue is at least .99998 sigma^2 (38 syncs per 4 s). The one-correction
service value transports from just after the previous applied magnetic
correction: 40 ms at 25 Hz, the whole 1-s window for the 1-Hz word. These
finite values audit algebra on one run; they are neither enclosures nor
all-time service certificates. The committed record is reproduced exactly
in CI (`--expect`), and `verify_diagnostic` checks every reported metric.

## 7. Sharp AW covariance ceiling and the storage route

**Lemma B.** In the default profile of `ou3-nuisance-upper-proof.md` (real
arithmetic, S_factor=1 so `Sigma_aw=sigma^2 I`, `sigma<=4`, `tau<=12` s,
additive pending sync, no frame/relock reconfiguration), every operation
from construction on keeps

`lambda_max(P_aw)<=max(2.2^2,(1+epsilon)16)=(1+epsilon)16`,

epsilon the process-covariance defect of `small_x_source_defect()`. For any
level `s>=(1+epsilon)sigma_j^2` covering the predictions and applied sync
targets after `t_0`,

`(lambda_max(P_aw(t))-s)_+<=exp(-(t-t_0)/6)(lambda_max(P_aw(t_0))-s)_+`.

*Proof.* The AW rows of the LIN transition are `phi e_a`, so prediction maps
the block to `phi^2 P_aw+Q_aa` with `Q_aa<=(1+epsilon)(1-phi^2)sigma^2 I`;
hence `m'<=phi^2 m+(1+epsilon)(1-phi^2)sigma^2` and `m'-s<=phi^2(m-s)`. The
pending sync adds `Pi_+(sigma_t^2 I-P_aw)`; for an isotropic target the
result has the eigenvectors of P_aw and eigenvalues `max(beta_i,sigma_t^2)`,
so `m'=max(m,sigma_t^2)`. Every applied acc, magnetic or S correction uses
`K=PC'S^-1` on the AW rows with the same, possibly bumped, S that the
Joseph update receives, so its AW block is `P_aw-K_a S K_a'<=P_aw`, whether or
not BA rows are frozen. Attitude resets, heel reframing, gyro-bias
projection and AG/LIN cross-block zeroing leave the block unchanged, and
entering Live seats it at Sigma.

*Tightness.* After every applied isotropic sync each eigenvalue of P_aw is
at least `sigma_t^2`, and syncs run at the .1-s adaptation cadence. No
covariance ceiling below the stationary variance therefore exists, whatever
the acc/S rows do between syncs. Isotropy is necessary: with S_factor=2,
`Sigma=diag(16,16,4)` and `P=Sigma+q_1q_1'-2q_2q_2'`, `q_1=(3,0,4)/5`,
`q_2=(-4,0,3)/5`, satisfy `0<=P<=16 I`, yet the synchronized block has
(1,1) entry 16+9/25. `aw_covariance_ceiling.py` reproduces these exactly.

*Entry into the tail.* With `|e_aw|^2<=lambda_max(P_aw) V`, Corollary A holds
when `sqrt(V)<1.12383/sqrt((1+epsilon)sigma_max^2)`: at least .280954 at the
4 m/s^2 clamp, against .0072 from the inherited 156 m/s^2, and at least 1
when `sigma_max<=1.1238`. The ceiling is not propagated into the nuisance
upper comparison, where it does not control any threshold.

**The uniform storage route fails on carried words.** A uniform retained
radius r must contain `V(t)>=|e_aw(t)|^2/lambda_max(P_aw(t))` along the
history, while the route needs `r^2 sup_t lambda_max(P_aw)<1.12383^2`. On the
carried collinear motion, 200-Hz accelerometer corrections collapse
`lambda_max(P_aw)` to .0026 between syncs while the jerk-driven lag error
stays near .56 m/s^2, and syncs restore .0791. The ratio
`sup(|e_aw|^2/lambda_max) sup lambda_max/1.12383^2` is 6.68 with 1-Hz
magnetic corrections and 6.91 with the 25-Hz cadence of the quiet and
wave words, although the actual AW error .562 m/s^2 satisfies Corollary A.
Section 9 shows that the pointwise physical premise itself is false on
admitted histories, and section 8 removes it.

## 8. Corollary A*: the attitude columns need only the nominal signed mean

The rows are nominal: `u_i=(a_hat_i-g e_z)/g`. With convex weights alpha_i,
`mu=sum alpha_i a_hat_i`, `m=|mu|`, `m_perp=|mu x b|` and one applied magnetic
row, convexity and the eigenvalue bound of section 4 give

`G >= gamma* I`, `gamma*=(sigma_w-m_perp/g)_+^2/((1+m/g)^2+1)`.

No physical acceleration, attitude, attitude error or magnetic residual
appears, and no pointwise premise is needed: four rows with `a_hat=(+-8,0,0)`
and `(0,+-6,+-1)` have zero mean and exactly dominate `gamma*=49/1250` for
`b=(7,0,24)/25`. Positivity needs only `m_perp<g sigma_w=1.96133 m/s^2`, and
`gamma*(0)=1/50`, 3.3 times Corollary A's zero-error value. The physical
transfer `|mu|<=2V/T+J h/4+|sum alpha_i e_i|` shows that Corollary A's
threshold 1.12383 applies to the **signed mean** error, not to `sup|e|`.
The nominal mean is a nonnegatively weighted signed sum of AW corrections:
`sum alpha_i a_hat_i=W_0 a_hat_0+sum_c W_c Delta_c`, `0<=W_c<=1`
(exact, `aw_tracking.py`). A source-uniform bound on it is OPEN.

## 9. The literal AW correction loop

In the world frame the accelerometer innovation factors exactly as
`r=R_hat(y-a_hat)`, `y=R_hat'(f-b_hat)+g e_z=a+eta`, with attitude, bias and
residual only in `eta`. Hence `e+=(I-Gamma)e+Gamma eta` (Gamma=K_a R_hat),
`e+=e+xi` for S/mag corrections and `e+=phi e-(1-phi)a-Delta a` for the OU
prediction, and exactly

`sum_acc Gamma(e-eta)=e_0-e_N+sum xi-sum_pred[(1-phi)a_hat+Delta a]`,

with `sum Delta a=a_N-a_0`: the physical increments telescope, but the
average is **gain-weighted**, not trapezoidal. As for injections,
`Delta a_hat Delta a_hat'<=NIS K_a S K_a'<=NIS P_aw`; at the .05 sigma floor
(after the Lemma B excess decays) following a jerk-limit ramp needs
`NIS>=144` at every step.

**Pointwise premise refuted (DEAD_END 20).** The horizontal AW prior scale
is `S_factor sigma_Z`, and the tuner reads the vertical channel. On an admitted
history with a C2 onset (exact rational MARINE MOTION envelope) and a
horizontal triangle of amplitude 8.7 m/s^2 at 2.8 Hz (jerk bound 99.9),
the unchanged shipping estimator reaches A21 and keeps `sigma_aw=.05`:
`sup|a_hat-a|=7.647 m/s^2`, 6.80 times 1.12383, with attitude error
.0062 rad. Signed 16-s means stay small: error .0327 m/s^2, nominal
transverse mean .0225. The carried audit
(`aw-tracking-source-feasibility.json`) covers six profiles:

| Profile | sup error / 1.12383 | signed mean / 1.12383 | m_perp / 1.96133 |
|---|---:|---:|---:|
| horizontal triangle at jerk limit | 6.80 | .029 | .011 |
| triangle with vertical swell | 2.00 | .0074 | .081 |
| triangle near adaptation cadence | 2.28 | .050 | .103 |
| collinear witness motion | .216 | .00075 | .00035 |
| large swell with chop | 1.89 | .0090 | .036 |
| sync-locked rectification | 2.77 | .330 | .177 |

**Rectification.** The AW covariance sync fires every 21 samples (.105 s);
P_aw and hence Gamma are periodic in that cycle. A jerk-limited triangle
phase-locked to it (2.6 m/s^2, 200/21 Hz, phase .4) raises the signed mean
error from <=.05 to .371 m/s^2: the time-varying gain rectifies. A raised
tuner sigma reduces it (.10 and .025). It remains below both thresholds;
the worst nominal transverse mean is .348 m/s^2 and the L1 force mean 1.091.
These finite replays are neither enclosures nor service certificates.

## 10. Signed world-injection transport

**Lemma I\*.** Let `M=R_hat'R` be the world attitude error. An injection acts
by `M<-Exp(-x)M`, a prediction by `M<-P M` with `angle(P)<=h|omega_tilde|` up to
the branch defect. With J the ordered product of the injection rotations
alone, `M_N M_0'=J K` where K is the ordered product of the rate rotations,
each conjugated by the preceding injection product. By bi-invariance

`angle(J)<=angle(M_0)+angle(M_N)+int|omega_tilde|`.

This is exact on the rotation group and sums no injection norm. A single
reset gives `|x|<=alpha_-+alpha_+`, so the relaxed cancellation of DEAD_END
17 (|d|=4) needs a world attitude error of at least 1.1416 rad.

**Reset factor.** `N=Exp(-[x])(I+[x]/2)=I-[x]/2+R`,
`R=sum_{n>=3}(1-n/2)(-X)^n/n!`, `|R|<=|x|^3/6` for `|x|<=1`: the quadratic term
cancels, so A~ rotates by half the injection rotation to first order. The
half-angle product is not a group identity; its second-order remainder is
`Sum_l |x_l||S_(l-1)|/4` (signed partial sums S) plus the stretch
`exp(Sum|x|^2/8)-1`, and has no source bound (OPEN).

**Corollary A\*\*.** For rotation factors Q_k in the rows,
`|u x Q theta|=|Q'u x theta|` and convexity give the floor
`|Q_j v_bar x b|^2/(|v_bar|^2+1)` with
`Q_j v_bar=u_bar+sum alpha_k(Q_j Q_k'-I)u_k`: the loss is the **signed** mean of
relative rotations. Opposite rotations cancel exactly in the supplied check;
one-sided ones reduce the floor.

**Feasibility (DEAD_END 21).** Under the deterministic premises the
perturbative charge `|A~_k-A~_j|` against `sqrt(gamma*)/u_rms` fails on 16-s
windows: ratio^2 5.12 at the retained tilt .1047; at vanishing radius .96
only for a quiet force RMS with an exact bias estimate, otherwise 1.79--7.17;
700--1306 with only the .5 rad/s invariant. The gyro residual alone rotates
the injections by .02 rad/s. The signed rotating frame is feasible at
vanishing radius (.22--.64) and at the retained tilt for moderate
`mean|a_hat|` (.77; 1.006 at 3 m/s^2), and fails with only the invariant.
`signed-injection-certificate.json` records the table.

## 11. Tube lemma: what prevents gyro-column cancellation

Without injections, `c(t)=theta+Gamma(t)beta`, `Gamma(t)=int_0^t R_hat(s)'ds`,
is a curve of speed |beta| whose unit tangent `T=R_hat'beta/|beta|` turns at
most `Omega=|omega_hat|_max` per second. MAGNETIC SERVICE holds on **every**
interval of length T_M=1 s, so each contains an applied magnetic row and the
row gaps are Delta<=1 s. The implemented .5 rad/s invariant gives
`Omega<=.6108652381980153+.02+.02+.5`.

**Lemma T.** If `|P_b c(t_j)|<=eps|beta|` at every magnetic row, then on the
interior (distance >=pi/(2 Omega) from the word ends)

`|b.T(t)| >= c0 = cos(Omega Delta)-Omega eps`.

*Proof.* Let `s=|b.T(t*)|`, `c=arcsin s` and e the unit transverse direction
of `T(t*)`. For `|u-t*|<=R=(pi/2-c)/Omega`, `angle(T(u),e)<=Omega|u-t*|+c<=pi/2`.
If `c<=pi/2-Omega Delta`, samples `t_j in [t*-R,t*-R+Delta]`,
`t_l in [t*+R-Delta,t*+R]` exist and

`2 eps>=e.(c(t_l)-c(t_j))/|beta|>=int_{|u-t*|<=R-Delta}cos(Omega|u-t*|+c)du
 =(2/Omega)(cos(Omega Delta)-sin c)`.

Otherwise `s>cos(Omega Delta)` already. With less room l (possible only when
`Omega l<pi/2`), any `eps<r<=l-Delta` gives `2r cos(c+Omega r)<=2 eps`, i.e.
`s>=sqrt(1-x^2)cos(Omega r)-x Omega r`, `x=eps/r`; the lemma uses the smaller
of both bounds (`tube_constant`). Since T is continuous, `b.T` keeps one
sign: the field-axis coordinate `eta=b.c` is monotone with speed >=c0|beta|.

At the invariant rate `cos(Omega Delta)>=.4078` (rational alternating Taylor
bound). The condition matters: a circle of radius 1/Omega tangent to the
b-line returns to it every `2 pi/Omega=5.46 s` while `b.T` changes sign.
Monotonicity is what excludes cancellation: the relaxed DEAD_END 17 reset
sequence rotates the tangent by radians, which Lemma I* ties to attitude error.

## 12. Theorem G0: injection-free aggregate six-column floor

Take the injection-free world array (`A~=I`, `B~=Gamma`) on a word whose
interior contains accelerometer windows W1, W2 of length L separated by a gap
G, with at least n rows each, and magnetic row gaps <=Delta. Suppose each
window has nominal transverse mean `m_perp` and L1 force mean
`u1=mean|f_hat|/g`. Put `mu=sigma_w-m_perp/g`,
`rho=eps+(Delta/2)sqrt(1-c0^2)`.

* **Gyro part.** Either some magnetic residual exceeds `eps|beta|`, or Lemma T
  applies; then one window lies at distance >=G/2 from the zero of eta, so
  `|eta|>=c0 G|beta|/2` there, `|P_b c|<=rho|beta|` at its rows, and Jensen for
  `(.)_+^2` gives `Q>=q_I|beta|^2`,
  `q_I=min(B_min^2 eps^2, n g^2(c0 G mu/2-rho u1)_+^2)`.
* **Attitude part.** On W1 with centre t_c and one magnetic row within
  Delta/2, convexity, the reverse triangle inequality and `|Gamma(t)|<=t` give
  `Q>=K(sqrt(gamma*)|theta|-C|beta|)_+^2`, `K=min(n g^2,B_min^2)`,
  `gamma*=mu^2/(u1^2+1)`, `C=sqrt(gamma*) t_c+sqrt((u1 L/2)^2+(Delta/2)^2)`.
* **Combination.** For coordinates (theta, lambda beta),
  `s^2>=q_I k^2/(1+lambda^2 k^2)`, `k=sqrt(K gamma*)/(sqrt(q_I)+sqrt(K) C)`.

With the invariant Omega, Delta=1 s, L=16 s, G=64 s, B_min=20 uT,
n=16/.006, eps=41/200 and the premises `m_perp<=2/5`, `u1<=6/5` (carried worst
.348 and 1.091): `c0>=.171694`, `q_I=16.81` and `s^2>=1.486786e-3`
(`aggregate-floor-certificate.json`, exact rationals; 2 s of room before W1
and after W2, which exceeds pi/(2 Omega)). Disjoint 1-s service windows
(Delta=2) fail at the invariant rate, and at the physical rate they need more
than 2 s of room. An 80-digit tube audit and a synthetic falsification audit
(field-axis spin, transverse roll, frozen attitude, random rates) find actual
sigma_min 37--60 against recomputed floors .067--.070: valid and about 600
times conservative.

G0 is **not** a theorem about the literal array: the injection transport A~
is excluded, and on a 100-s word the deterministic .02 rad/s gyro residual
alone lets the injection frame, hence the kernel direction `Q'b`, rotate by
about 1 rad, which the global fixed-b tube does not tolerate. Its nominal
window premises are measured, not proved. On the six carried words the AW
audit rebuilds both arrays over 100 s from the literal event stream
(predictions, every reset, applied rows): the literal sigma_min is 191.8--267.4,
at least the injection-free value (ratio >=1.0000089), `|A~-I|<=.0045`, and
G0's floor from each word's own premises (service gap 1 s, B_min=20 uT) is
.072--.084. In practice the injection frame is benign; the proof gap is the
deterministic worst case.

## 13. Downstream feasibility

On a 100-s word (m<=150000 AG rows, >=33333 operations) G0 gives six pivots
>=9.96e-5. The existing scalar reader ceiling multiplies a coefficient bound
>=20 over every operation: `log10 B_*>=43367`. A least-squares reader has
action <=`R_max/s^2`=2690 (R_max=4, the magnetic residual bound) and gyro
block <=`R_max/q_I`=.238. Iterating the corrected-word prediction comparison
over every prediction bounds the contraction margin by `T b0/P_bg,max`:
3.7e-13 (full action), 4.2e-9 (gyro block) and 3.4e-7 even at the actual
floor 37 (b0=1e-11). Classification (DEAD_END 22): the chain is strict only
in the existence sense; the covariance-ceiling/process-noise route cannot
deliver a useful rho_0 unless the ceiling is near the true P_bg. The
contraction must come from the measurement loss J_AG in the slow directions,
with a blockwise (not least-singular-value) reader action.

The same product appears in information form: for a noise-free AG block
the homogeneous word factor is `1/(1+lambda_min(P_root^(1/2) I_W P_root^(1/2)))`
with `I_W=O'R^-1 O`, so the margin is set by the lower covariance times the
information in each slow direction. A lower bound built from process
accumulation alone (`P_bg>=T b0`) reproduces `T b0 q_I/R_max`; a Riccati
comparison with a maximal information rate would instead give
`P_bg>=sqrt(b0/i_max)` and a margin of order `T sqrt(b0 i_max)`. This is a
proposed formulation, not a certified bound.

## 14. What remains

1. A source-uniform bound on the nominal window statistics (signed mean,
   L1 force) of the literal AW loop; the pointwise physical premise is false.
2. Theorem G0 with injections: a local-tube version tolerating the kernel
   drift `Q'b` and the half-angle remainder of section 10.
3. A blockwise reader action and an information-based contraction in the
   gyro-bias and accelerometer-bias directions.
4. A same-cell floor is not needed; the aggregate route uses only the
   1-s service gap.
