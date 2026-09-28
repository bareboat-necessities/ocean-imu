# World-frame historical rows, injection budget and the same-cell limitation

These results enter `V_next <= rho V + c_d |d|^2` only through the historical
reader: the six-column floor (c, a, b0) feeds B_*, J_AG, the full covariance
upper bound and rho_0<1. The injection budget also bounds the injection-squared
term of the finite-angle reset remainder in the nonlinear supply eta. None of
them supplies a source-uniform premise by itself. The AW covariance ceiling
(section 7) enters through Corollary A's AW tracking premise and the
coercivity bound `|e_aw|^2<=lambda_max(P_aw) V`. Exact checks are in
`world_frame.py` and `aw_covariance_ceiling.py`; the carried audit is
`world_frame_source_diagnostic.py`.

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
  `Gamma=R' Rs^-1 Bs`. If Rs equals the mean increment, W=I and Gamma is the
  world integral of the nominal body-to-world rotation over the step;
* reset: `A~+=N A~`, `B~+=N B~`, `N=R+' G R`. For any orthogonal mean
  injection `N'N=I+(|x|^2 I-xx')/4`, `x=R' d`, so sigma_min(N)=1,
  `|N|=sqrt(1+|x|^2/4)` and `|N^-1|=1`. For the exponential injection,
  `N=Exp(-[x])(I+[x]/2)` acts transversally as `(1+i theta/2)exp(-i theta)`.

*Proof.* `[Rv]x=R[v]x R'` for rotations; substitute the definitions. The closed
form gives `Rs^-1 Bs=int_0^h Exp([w]u) du`. Finally `N'N=R' G'G R` and
`G'G=I+(|d|^2 I-dd')/4`.

Consequently the attitude estimate and its error enter the six-column array
only through the world injections x (via N) and through Gamma, which
integrates the nominal rotation. With no injections the rows are
`[f_k]x [I, int R_hat']`: the world observability matrix of
`theta_w(t)=theta_w(0)+(int R_hat')b_g`. The source mean uses the polynomial
quaternion for |w h|<.01 while the covariance uses Rodrigues for |w|>=1e-7;
|W-I| is that branch and arithmetic defect.

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
and shows why that route cannot be uniform. This supplies attitude columns
of the aggregate premise; the gyro columns need the time structure of Gamma
and are open.

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
eigenvalue 2.05e-5 against `mu_M=1`. So no refutation of the same-cell route
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
`Exp([d]x)` injection, at most 1.6e-9 rad) on top of psi(theta). Quiet and
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
| one-correction service lambda_min | 1.6e-83 | 1.8e-10 | 2.05e-5 | 6.9e-11 |

Prediction world discrepancy is at most 2.1e-9 (float branch), the reset Gram
identity holds to 4.3e-81 and row factorization to 9.9e-78. `|A~-I|` matches
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
The AW error must be bounded in physical units from the literal correction
loop, not through covariance-normalized storage.

## 8. What remains

1. A physical AW tracking bound `epsilon_a<1.12383 m/s^2` on 16-s windows
   from the literal acc/S correction loop driven by jerk, bias and attitude
   error. No covariance ceiling can supply it uniformly (section 7).
2. The aggregate gyro-column floor from transverse-force and magnetic rows at
   separated times, with nominal rotation up to 1.15 rad/s in Gamma.
3. The signed world injection sum in place of norm-summed per-correction
   covariance/NIS bounds, which overcharge the 3-s inter-anchor budget by about
   460 on the collinear word.
4. A same-cell floor, if pursued, from MAGNETIC SERVICE density coupled to the
   jerk-limited collinear cadence; the aggregate route does not need it.
